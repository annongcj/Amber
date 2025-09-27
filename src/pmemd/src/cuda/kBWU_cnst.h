#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// This is included by kBWU.h with options for #define'd LOCAL_AFE and LOCAL_MBAR.  A third
// define, LOCAL_ENERGY, applies at the level of kBWU.h but cascades down into this code as
// well.  This code does not constitute an entire kernel, but rather a section of one.
//
// Variables already defined before including this code:
//   atmcrd[x,y,z]:  coordinates for imported atoms
//   atmfrc[x,y,z]:  force accumulators for imported atoms
//   nrgACC:         accumulator for energy terms
//   startidx:       starting index for this warp's set of GRID (that is, 32) bonded terms
//   tgx:            thread index within the warp
//---------------------------------------------------------------------------------------------
{
  unsigned int atmI = cSim.pBwuCnstID[startidx + tgx];
#ifdef LOCAL_ENERGY
  if (atmI == 0xffffffff) {
    eterm = (PMEAccumulator)0;
  }
  else {
#else
  if (atmI != 0xffffffff) {
#endif
    PMEDouble2 constraint1 = cSim.pBwuCnst[2*startidx + tgx];
    PMEDouble2 constraint2 = cSim.pBwuCnst[2*startidx + GRID + tgx];
    PMEDouble ax    = atmcrdx[atmI] - constraint1.y;
    PMEDouble ay    = atmcrdy[atmI] - constraint2.x;
    PMEDouble az    = atmcrdz[atmI] - constraint2.y;
    PMEDouble wx    = constraint1.x * ax;
    PMEDouble wy    = constraint1.x * ay;
    PMEDouble wz    = constraint1.x * az;
#ifdef LOCAL_AFE
    unsigned int TIstatus = cSim.pBwuCnstStatus[startidx + tgx];
    TIregion = TIstatus >> 16;
    CVterm = TIstatus & 0xff;
    PMEDouble lambda = TISetLambda(CVterm, TIregion, cSim.AFElambdaSP[1]);
    wx *= lambda;
    wy *= lambda;
    wz *= lambda;
#endif
#ifdef use_DPFP
    PMEAccumulator ifx = llrint(wx * (PMEDouble)-2.0 * FORCESCALE);
    PMEAccumulator ify = llrint(wy * (PMEDouble)-2.0 * FORCESCALE);
    PMEAccumulator ifz = llrint(wz * (PMEDouble)-2.0 * FORCESCALE);
    atomicAdd((unsigned long long int*)&atmfrcx[atmI], llitoulli(ifx));
    atomicAdd((unsigned long long int*)&atmfrcy[atmI], llitoulli(ify));
    atomicAdd((unsigned long long int*)&atmfrcz[atmI], llitoulli(ifz));
#else
    int ifx = __float2int_rn(wx * (PMEDouble)-2.0 * BSCALE);
    int ify = __float2int_rn(wy * (PMEDouble)-2.0 * BSCALE);
    int ifz = __float2int_rn(wz * (PMEDouble)-2.0 * BSCALE);
    atomicAdd(&atmfrcx[atmI], ifx);
    atomicAdd(&atmfrcy[atmI], ify);
    atomicAdd(&atmfrcz[atmI], ifz);
#endif
#ifdef LOCAL_ENERGY
#  ifdef LOCAL_AFE
    SCterm = (TIstatus >> 8) & 0xff;
    if (SCterm == 0) {
      eterm = llrint(ENERGYSCALE * lambda * (wx*ax + wy*ay + wz*az));
    }
    else {

      // There is no soft-core energy decomposition for positional constraints
      eterm = (PMEAccumulator)0;
    }
    edvdl = llrint((PMEDouble)(CVterm * (2*TIregion - 3)) * ENERGYSCALE *
                   (wx*ax + wy*ay + wz*az));
#    ifdef LOCAL_MBAR
    // Set pieces of information that MBAR will need based on this thread's results
    mbarTerm = wx*ax + wy*ay + wz*az;
    mbarRefLambda = lambda;
#    endif
#  else
    eterm = llrint(ENERGYSCALE * (wx*ax + wy*ay + wz*az));
#  endif
  }
#  define CNST_CASE
#  define EACC_OFFSET   CNST_EACC_OFFSET
#include "kBWU_EnergyReduction.h"
#  undef EACC_OFFSET
#  undef CNST_CASE
#else
  }
#endif // LOCAL_ENERGY
}
