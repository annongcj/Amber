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
//   curr_wuidx:     CURRENT work unit index (as bonds write directly back to the global force
//                   arrays).  The counter wuidx will have previously been incremented, so rely
//                   on curr_wuidx instead.
//   tgx:            thread index within the warp
//   warpIdx:        warp index within the block
//   eterm:          PMEAccumulator into which the energy of the interaction will be cast
//
// Variables contingent on LOCAL_AFE and possibly LOCAL_MBAR #defines:
//   TIregion:       thermodynamic integration region code (0-2)
//   SCterm:         flag to indicate that the interaction is in the soft-core region
//   CVterm:         flag to indicate that the interaction is in the CV region
//   mbarTerm:       critical energy term for MBAR summation
//   mbarRefLambda:  reference mixing value for each interaction during MBAR summation
//---------------------------------------------------------------------------------------------
{
  unsigned int rawID = cSim.pBwuBondID[startidx + tgx];
#ifdef LOCAL_ENERGY
  if (rawID == 0xffffffff) {
    eterm = (PMEAccumulator)0;
  }
  else {
#else
  if (rawID != 0xffffffff) {
#endif
    unsigned int atmI  = rawID >> 8;
    unsigned int atmJ  = rawID & 0xff;
    PMEDouble xij      = atmcrdx[atmI] - atmcrdx[atmJ];
    PMEDouble yij      = atmcrdy[atmI] - atmcrdy[atmJ];
    PMEDouble zij      = atmcrdz[atmI] - atmcrdz[atmJ];
    PMEDouble rij      = sqrt(xij*xij + yij*yij + zij*zij);
    PMEDouble2 bond    = cSim.pBwuBond[startidx + tgx];
    PMEDouble da       = rij - bond.y;
    PMEDouble df       = bond.x * da;
    PMEDouble dfw      = (df + df) / rij;
#ifdef LOCAL_AFE
    unsigned int TIstatus = cSim.pBwuBondStatus[startidx + tgx];
    TIregion = TIstatus >> 16;
    CVterm = TIstatus & 0xff;
    PMEDouble lambda = TISetLambda(CVterm, TIregion, cSim.AFElambdaSP[1]);
    dfw *= lambda;
#endif
    PMEDouble fx = dfw * xij;
    PMEDouble fy = dfw * yij;
    PMEDouble fz = dfw * zij;
#ifdef LOCAL_NEIGHBORLIST
    int offset = (3*curr_wuidx + 2) * BOND_WORK_UNIT_THREADS_PER_BLOCK;
    int globalPosI = cSim.pBwuInstructions[offset + atmI];
    int globalPosJ = cSim.pBwuInstructions[offset + atmJ];
#else
    int offset = (3*curr_wuidx + 1) * BOND_WORK_UNIT_THREADS_PER_BLOCK;
    int globalPosI = cSim.pBwuInstructions[offset + atmI];
    int globalPosJ = cSim.pBwuInstructions[offset + atmJ];
#endif
    PMEAccumulator ifx = llrint(fx * FORCESCALE);
    PMEAccumulator ify = llrint(fy * FORCESCALE);
    PMEAccumulator ifz = llrint(fz * FORCESCALE);
    atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[globalPosJ],
              llitoulli(ifx));
    atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[globalPosJ],
              llitoulli(ify));
    atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[globalPosJ],
              llitoulli(ifz));
    atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[globalPosI],
              llitoulli(-ifx));
    atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[globalPosI],
              llitoulli(-ify));
    atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[globalPosI],
              llitoulli(-ifz));
#ifdef LOCAL_ENERGY
#  ifdef LOCAL_AFE
    SCterm = (TIstatus >> 8) & 0xff;
    if (SCterm == 0) {
      eterm = llrint(ENERGYSCALE * lambda * df * da);
    }
    else {
      scEcomp[TIregion - 1] = llrint(ENERGYSCALE * df * da);
      eterm = (PMEAccumulator)0;
    }
    edvdl = llrint((PMEDouble)(CVterm * (2*TIregion - 3)) * ENERGYSCALE * df * da);
#    ifdef LOCAL_MBAR
    // Set pieces of information that MBAR will need based on this thread's results
    mbarTerm = df * da;
    mbarRefLambda = lambda;
#    endif
#  else
    eterm = llrint(ENERGYSCALE * df * da);
#  endif
  }
#  define EACC_OFFSET   BOND_EACC_OFFSET
#  define STORE_SCR1    cSim.pSCBondR1
#  define STORE_SCR2    cSim.pSCBondR2
#include "kBWU_EnergyReduction.h"
#  undef EACC_OFFSET
#  undef STORE_SCR1
#  undef STORE_SCR2
#else
  }
#endif // LOCAL_ENERGY
}
