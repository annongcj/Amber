#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// This is included by kBWU.h with options for #define'd LOCAL_AFE and LOCAL_MBAR and an
// assumption that LOCAL_ENERGY is #define'd.  This code does not constitute an entire kernel,
// but rather a common refrain that is included in many pieces of a much larger kernel.  The
// vast majority of the code is wrapped inside a conditional dependent on there being at least
// one thread that has done a Thermodynamic Integration-related interaction, but the entire
// warp moves as one through this code (no divergence).
//
// There are two additional, optional #defines, NB14_CASE and CNST_CASE.  If NB14_CASE is
// present, then additional energy storage will be undertaken for a second accumulator,
// specifically named for non-bonded 1:4 electrostatics.  If CNST_CASE is NOT present, the
// usual storage of soft-core energy components will take place, but if it is #define'd, then
// soft-core energies are not applicable to positional constraints and anything related to
// this part is skipped.  In the special CNST_CASE, STORE_SCR1 and STORE_SCR2 need not be
// defined.
//
// Variables already defined before including this code:
//   nrgACC:         accumulator for energy terms
//   tgx:            thread index within the warp
//   warpIdx:        warp index within the block
//   eterm:          energy contribution accumulated from standard evaluation of the bonded
//                   interaction (bond, angle, dihedral, etc.)
//   edvdl:          energy contribution to dV/dLambda from each thread's bonded interaction
//                   (only needed if LOCAL_AFE is #define'd)
//   TIregion:       label for the Thermodynamic Integration region each thread's object
//                   belongs to (only needed if LOCAL_AFE is #define'd)
//   mbarTerm:       term central to Multi-state Bennet Acceptance Ratio summation
//   mbarRefLambda:  lambda at which each thread computed its interaction (this will vary with
//                   the mixing factor and TIregion)
//   sBarTot:        __shared__ memory BAR accumulators for each state
//
// Additional macros that must be #define'd immediately prior to including this file and
// should be immediately #undef'd afterwards:
//   EACC_OFFSET:    offset in the energy accumulator array nrgACC at which to begin merging
//   STORE_SCR1:     pointer to the first soft-core energy component's global accumulator
//   STORE_SCR2:     pointer to the second soft-core energy component's global accumulator
//---------------------------------------------------------------------------------------------
{
  for (unsigned int stride = warpSize >> 1; stride > 0; stride >>=1) {
    #ifdef NB14_CASE
    eel14 += __SHFL_DOWN(WARP_MASK, eel14, stride);
    #endif
    eterm += __SHFL_DOWN(WARP_MASK, eterm, stride);
  }
  if (tgx == 0) {
#ifdef NB14_CASE
    nrgACC[SCEE_EACC_OFFSET + warpIdx] += eel14;
#endif
    nrgACC[EACC_OFFSET + warpIdx] += eterm;
  }
#ifdef LOCAL_AFE
  PMEMask TIballot = __BALLOT(WARP_MASK, TIregion);
  if (TIballot) {
    for (unsigned int stride = warpSize >> 1; stride > 0; stride >>=1) {
      #ifdef NB14_CASE
      scEel14comp[0] += __SHFL_DOWN(WARP_MASK, scEel14comp[0], stride);
      scEel14comp[1] += __SHFL_DOWN(WARP_MASK, scEel14comp[1], stride);
      #endif
      #ifndef CNST_CASE
      scEcomp[0] += __SHFL_DOWN(WARP_MASK, scEcomp[0], stride);
      scEcomp[1] += __SHFL_DOWN(WARP_MASK, scEcomp[1], stride);
      #endif
      edvdl += __SHFL_DOWN(WARP_MASK, edvdl, stride);
    }
    if (tgx == 0) {

      // These atomics to global are not a bottleneck, as they happen inside
      // the warp vote function over TIballot (only warps with TI interactions
      // will reduce their sums and then contribute).
      atomicAdd(cSim.pDVDL, llitoulli(edvdl));
#  ifdef NB14_CASE
      atomicAdd(cSim.pSCEEL14R1, llitoulli(scEel14comp[0]));
      atomicAdd(cSim.pSCEEL14R2, llitoulli(scEel14comp[1]));
#  endif
#  ifndef CNST_CASE
      atomicAdd(STORE_SCR1, llitoulli(scEcomp[0]));
      atomicAdd(STORE_SCR2, llitoulli(scEcomp[1]));
#  endif
    }
#  ifdef LOCAL_MBAR
    // Do a second vote for CVterms.  A copy of this will get depleted
    // and refreshed each time a warp strides over the MBAR states.
    PMEMask CVballot = __BALLOT(WARP_MASK, CVterm);

    // cSim.bar_stride must be a multiple of GRID (the warp size) to avoid warp divergence
    int j = tgx;
    PMEMask mask1 = __BALLOT(WARP_MASK, j < cSim.bar_stride);
    while (j < cSim.bar_stride) {

      // Working on the jth MBAR state, sum up the contributions
      // from all relevant interactions in a register variable.
      // Dump the result into the __shared__ memory accumulator.
      PMEAccumulator bar_contrib = (PMEAccumulator)0;
      PMEMask CVreport = CVballot;
      for (int i = 0; i < GRID; i++) {

        // The warp will not diverge at this conditional over CVreport, nor
        // the one above over TIballot.  It is safe to use __shfl intrinsics.
        if (CVreport & 0x1) {
          int barLambdaOffset = (__SHFL(mask1, TIregion, i) - 1) * cSim.bar_stride;
          bar_contrib += llrint(ENERGYSCALE * __SHFL(mask1, mbarTerm, i) *
                                (cSim.pBarLambda[barLambdaOffset + j] -
                                 __SHFL(mask1, mbarRefLambda, i)));
        }
        CVreport >>= 1;
      }
      atomicAdd((unsigned long long int*)&sBarTot[j], llitoulli(bar_contrib));
      j += GRID;
      mask1 = __BALLOT(mask1, j < cSim.bar_stride);
    }
#  endif // LOCAL_MBAR
  }
#endif // LOCAL_AFE
}
