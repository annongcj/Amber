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
  unsigned int rawID = cSim.pBwuQQxcID[startidx + tgx];
#ifdef LOCAL_VIRIAL
  PMEAccumulator v11 = (PMEAccumulator)0;
  PMEAccumulator v22 = (PMEAccumulator)0;
  PMEAccumulator v33 = (PMEAccumulator)0;
#endif
#ifdef LOCAL_ENERGY
  if (rawID == 0xffffffff) {
    eterm = (PMEAccumulator)0;
  }
  else {
#else
  if (rawID != 0xffffffff) {
#endif
    unsigned int atmI = rawID >> 8;
    unsigned int atmJ = rawID & 0xff;
    PMEFloat dx = atmcrdx[atmJ] - atmcrdx[atmI];
    PMEFloat dy = atmcrdy[atmJ] - atmcrdy[atmI];
    PMEFloat dz = atmcrdz[atmJ] - atmcrdz[atmI];
    PMEFloat qterm  = charges[atmJ] * charges[atmI];
    PMEFloat r2     = dx*dx + dy*dy + dz*dz;
    PMEFloat rinv   = rsqrt(r2);
    PMEFloat r2inv  = rinv * rinv;
#ifdef LOCAL_ENERGY
    PMEFloat g      = -rinv * qterm;
    PMEFloat df     = g * r2inv;
#else
    PMEFloat df     = -qterm * rinv * r2inv;
#endif
#ifdef use_DPFP
    PMEAccumulator ifx = llrint(df * dx * FORCESCALE);
    PMEAccumulator ify = llrint(df * dy * FORCESCALE);
    PMEAccumulator ifz = llrint(df * dz * FORCESCALE);
#  ifdef LOCAL_VIRIAL
    v11 -= llrint(df * dx * dx * FORCESCALE);
    v22 -= llrint(df * dy * dy * FORCESCALE);
    v33 -= llrint(df * dz * dz * FORCESCALE);
#  endif
    atomicAdd((unsigned long long int*)&atmfrcx[atmJ], llitoulli( ifx));
    atomicAdd((unsigned long long int*)&atmfrcy[atmJ], llitoulli( ify));
    atomicAdd((unsigned long long int*)&atmfrcz[atmJ], llitoulli( ifz));
    atomicAdd((unsigned long long int*)&atmfrcx[atmI], llitoulli(-ifx));
    atomicAdd((unsigned long long int*)&atmfrcy[atmI], llitoulli(-ify));
    atomicAdd((unsigned long long int*)&atmfrcz[atmI], llitoulli(-ifz));
#else  // use_DPFP
    int ifx = __float2int_rn(df * dx * BSCALEF);
    int ify = __float2int_rn(df * dy * BSCALEF);
    int ifz = __float2int_rn(df * dz * BSCALEF);
#  ifdef LOCAL_VIRIAL
    v11 -= fast_llrintf(df * dx * dx * FORCESCALEF);
    v22 -= fast_llrintf(df * dy * dy * FORCESCALEF);
    v33 -= fast_llrintf(df * dz * dz * FORCESCALEF);
#  endif
    atomicAdd(&atmfrcx[atmJ],  ifx);
    atomicAdd(&atmfrcy[atmJ],  ify);
    atomicAdd(&atmfrcz[atmJ],  ifz);
    atomicAdd(&atmfrcx[atmI], -ifx);
    atomicAdd(&atmfrcy[atmI], -ify);
    atomicAdd(&atmfrcz[atmI], -ifz);
#endif // use_DPFP
#ifdef LOCAL_ENERGY
    eterm = llrint(ENERGYSCALE * g);
  }
#  define EACC_OFFSET   QQXC_EACC_OFFSET

  // These two defines are placeholders
#  define STORE_SCR1    cSim.pSCVDW14R1
#  define STORE_SCR2    cSim.pSCVDW14R2
#include "kBWU_EnergyReduction.h"
#  undef STORE_SCR1
#  undef STORE_SCR2
  // The placeholder defines have now been undone

#  undef EACC_OFFSET
#else
  }
#endif // LOCAL_ENERGY
#ifdef LOCAL_VIRIAL
  for (unsigned int stride = warpSize >> 1; stride > 0; stride >>=1) {
    v11 += __SHFL_DOWN(WARP_MASK, v11, stride);
    v22 += __SHFL_DOWN(WARP_MASK, v22, stride);
    v33 += __SHFL_DOWN(WARP_MASK, v33, stride);
  }
  if (tgx == 0) {
    atomicAdd(cSim.pVirial_11, llitoulli(v11));
    atomicAdd(cSim.pVirial_22, llitoulli(v22));
    atomicAdd(cSim.pVirial_33, llitoulli(v33));
  }
#endif
}
