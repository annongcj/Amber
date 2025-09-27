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
  unsigned int rawID = cSim.pBwuAnglID[startidx + tgx];
#ifdef LOCAL_ENERGY
  if (rawID == 0xffffffff) {
    eterm = (PMEAccumulator)0;
  }
  else {
#else
  if (rawID != 0xffffffff) {
#endif
    unsigned int atmI = rawID >> 16;
    unsigned int atmJ = (rawID >> 8) & 0xff;
    unsigned int atmK = rawID & 0xff;
    PMEDouble xij     = atmcrdx[atmI] - atmcrdx[atmJ];
    PMEDouble yij     = atmcrdy[atmI] - atmcrdy[atmJ];
    PMEDouble zij     = atmcrdz[atmI] - atmcrdz[atmJ];
    PMEDouble xkj     = atmcrdx[atmK] - atmcrdx[atmJ];
    PMEDouble ykj     = atmcrdy[atmK] - atmcrdy[atmJ];
    PMEDouble zkj     = atmcrdz[atmK] - atmcrdz[atmJ];
    PMEDouble rij     = xij*xij + yij*yij + zij*zij;
    PMEDouble rkj     = xkj*xkj + ykj*ykj + zkj*zkj;
    PMEDouble rik     = sqrt(rij * rkj);
    PMEDouble2 angl   = cSim.pBwuAngl[startidx + tgx];
    PMEDouble cst     = min(pt999, max(-pt999, (xij*xkj + yij*ykj + zij*zkj) / rik));
    PMEDouble ant     = acos(cst);
    PMEDouble da      = ant - angl.y;
    PMEDouble df      = angl.x * da;
    PMEDouble dfw     = -(df + df) / sin(ant);
#ifdef LOCAL_AFE
    unsigned int TIstatus = cSim.pBwuAnglStatus[startidx + tgx];
    TIregion = TIstatus >> 16;
    CVterm = TIstatus & 0xff;
    PMEDouble lambda = TISetLambda(CVterm, TIregion, cSim.AFElambda[1]);
    dfw *= lambda;
#endif
    PMEDouble cik = dfw / rik;
    PMEDouble sth = dfw * cst;
    PMEDouble cii = sth / rij;
    PMEDouble ckk = sth / rkj;
    PMEDouble fx1 = cii * xij - cik * xkj;
    PMEDouble fy1 = cii * yij - cik * ykj;
    PMEDouble fz1 = cii * zij - cik * zkj;
    PMEDouble fx2 = ckk * xkj - cik * xij;
    PMEDouble fy2 = ckk * ykj - cik * yij;
    PMEDouble fz2 = ckk * zkj - cik * zij;
#ifdef use_DPFP
    PMEAccumulator ifx1 = llrint(fx1 * FORCESCALE);
    PMEAccumulator ify1 = llrint(fy1 * FORCESCALE);
    PMEAccumulator ifz1 = llrint(fz1 * FORCESCALE);
    PMEAccumulator ifx2 = llrint(fx2 * FORCESCALE);
    PMEAccumulator ify2 = llrint(fy2 * FORCESCALE);
    PMEAccumulator ifz2 = llrint(fz2 * FORCESCALE);
    atomicAdd((unsigned long long int*)&atmfrcx[atmI], llitoulli(ifx1));
    atomicAdd((unsigned long long int*)&atmfrcy[atmI], llitoulli(ify1));
    atomicAdd((unsigned long long int*)&atmfrcz[atmI], llitoulli(ifz1));
    atomicAdd((unsigned long long int*)&atmfrcx[atmK], llitoulli(ifx2));
    atomicAdd((unsigned long long int*)&atmfrcy[atmK], llitoulli(ify2));
    atomicAdd((unsigned long long int*)&atmfrcz[atmK], llitoulli(ifz2));
    atomicAdd((unsigned long long int*)&atmfrcx[atmJ], llitoulli(-(ifx1 + ifx2)));
    atomicAdd((unsigned long long int*)&atmfrcy[atmJ], llitoulli(-(ify1 + ify2)));
    atomicAdd((unsigned long long int*)&atmfrcz[atmJ], llitoulli(-(ifz1 + ifz2)));
#else
    int ifx1 = __float2int_rn(fx1 * BSCALE);
    int ify1 = __float2int_rn(fy1 * BSCALE);
    int ifz1 = __float2int_rn(fz1 * BSCALE);
    int ifx2 = __float2int_rn(fx2 * BSCALE);
    int ify2 = __float2int_rn(fy2 * BSCALE);
    int ifz2 = __float2int_rn(fz2 * BSCALE);
    atomicAdd(&atmfrcx[atmI], ifx1);
    atomicAdd(&atmfrcy[atmI], ify1);
    atomicAdd(&atmfrcz[atmI], ifz1);
    atomicAdd(&atmfrcx[atmK], ifx2);
    atomicAdd(&atmfrcy[atmK], ify2);
    atomicAdd(&atmfrcz[atmK], ifz2);
    atomicAdd(&atmfrcx[atmJ], -(ifx1 + ifx2));
    atomicAdd(&atmfrcy[atmJ], -(ify1 + ify2));
    atomicAdd(&atmfrcz[atmJ], -(ifz1 + ifz2));
#endif
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
#  define EACC_OFFSET   ANGL_EACC_OFFSET
#  define STORE_SCR1    cSim.pSCBondAngleR1
#  define STORE_SCR2    cSim.pSCBondAngleR2
#include "kBWU_EnergyReduction.h"
#  undef EACC_OFFSET
#  undef STORE_SCR1
#  undef STORE_SCR2
#else
  }
#endif // LOCAL_ENERGY
}
