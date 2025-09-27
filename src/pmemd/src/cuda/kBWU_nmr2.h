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
  int offset         = 4*startidx + tgx;
  unsigned int rawID = cSim.pBwuNMR2ID[offset];
  int Inc       = cSim.pBwuNMR2ID[offset +   GRID];
  int stepInit  = cSim.pBwuNMR2ID[offset + 2*GRID];
  int stepFinal = cSim.pBwuNMR2ID[offset + 3*GRID];
  int active    = (NMRnstep >= stepInit && (NMRnstep <= stepFinal || stepFinal <= 0));
#ifdef LOCAL_ENERGY
  if (rawID == 0xffffffff || active == 0) {
    eterm = (PMEAccumulator)0;
  }
  else {
#else
  if (rawID != 0xffffffff && active == 1) {
#endif
    // Read timed data
    offset                 += 2 * startidx;
    PMEDouble2 R1R2Slope    = cSim.pBwuNMR2[offset         ];
    PMEDouble2 R3R4Slope    = cSim.pBwuNMR2[offset +   GRID];
    PMEDouble2 K2K3Slope    = cSim.pBwuNMR2[offset + 2*GRID];
    PMEDouble2 R1R2Intrcpt  = cSim.pBwuNMR2[offset + 3*GRID];
    PMEDouble2 R3R4Intrcpt  = cSim.pBwuNMR2[offset + 4*GRID];
    PMEDouble2 K2K3Intrcpt  = cSim.pBwuNMR2[offset + 5*GRID];

    // Calculate increment
    double dstep = NMRnstep - (double)((NMRnstep - stepInit) % abs(Inc));

    // Calculate restraint values
    PMEDouble2 R1R2, R3R4, K2K3;
    R1R2.x = dstep*R1R2Slope.x + R1R2Intrcpt.x;
    R1R2.y = dstep*R1R2Slope.y + R1R2Intrcpt.y;
    R3R4.x = dstep*R3R4Slope.x + R3R4Intrcpt.x;
    R3R4.y = dstep*R3R4Slope.y + R3R4Intrcpt.y;
    if (Inc > 0) {
      K2K3.x = K2K3Slope.x*dstep + K2K3Intrcpt.x;
      K2K3.y = K2K3Slope.y*dstep + K2K3Intrcpt.y;
    }
    else {
      int nstepu = (NMRnstep - stepInit) / abs(Inc);
      K2K3.x = K2K3Intrcpt.x * pow(K2K3Slope.x, nstepu);
      K2K3.y = K2K3Intrcpt.y * pow(K2K3Slope.y, nstepu);
    }

    // Compute displacement, force, and energy
    unsigned int atmI = rawID >> 8;
    unsigned int atmJ = rawID & 0xff;
    PMEDouble xij     = atmcrdx[atmI] - atmcrdx[atmJ];
    PMEDouble yij     = atmcrdy[atmI] - atmcrdy[atmJ];
    PMEDouble zij     = atmcrdz[atmI] - atmcrdz[atmJ];
    PMEDouble rij     = sqrt(xij*xij + yij*yij + zij*zij);
    PMEDouble df;
#ifdef LOCAL_ENERGY
    PMEDouble e;
#endif
    if (rij < R1R2.x) {
      PMEDouble dif = R1R2.x - R1R2.y;
      df            = (PMEDouble)2.0 * K2K3.x * dif;
#ifdef LOCAL_ENERGY
      e             = df*(rij - R1R2.x) + K2K3.x*dif*dif;
#endif
    }
    else if (rij < R1R2.y) {
      PMEDouble dif = rij - R1R2.y;
      df            = (PMEDouble)2.0 * K2K3.x * dif;
#ifdef LOCAL_ENERGY
      e             = K2K3.x * dif * dif;
#endif
    }
    else if (rij < R3R4.x) {
      df            = (PMEDouble)0.0;
#ifdef LOCAL_ENERGY
      e             = (PMEDouble)0.0;
#endif
    }
    else if (rij < R3R4.y) {
      PMEDouble dif = rij - R3R4.x;
      df            = (PMEDouble)2.0 * K2K3.y * dif;
#ifdef LOCAL_ENERGY
      e             = K2K3.y * dif * dif;
#endif
    }
    else {
      PMEDouble dif = R3R4.y - R3R4.x;
      df            = (PMEDouble)2.0 * K2K3.y * dif;
#ifdef LOCAL_ENERGY
      e             = df*(rij - R3R4.y) + K2K3.y*dif*dif;
#endif
    }
#ifdef LOCAL_AFE
    unsigned int TIstatus = cSim.pBwuNMR2Status[startidx + tgx];
    TIregion = TIstatus >> 16;
    CVterm = TIstatus & 0xff;
    PMEDouble lambda = TISetLambda(CVterm, TIregion, cSim.AFElambdaSP[1]);
#  ifdef LOCAL_ENERGY
    SCterm = (TIstatus >> 8) & 0xff;
    scEcomp[TIregion - SCterm] = llrint(ENERGYSCALE * e * (PMEDouble)SCterm);
    edvdl = llrint((PMEDouble)(CVterm * (2*TIregion - 3)) * ENERGYSCALE * e);
    //need to remove softcore terms from e but not scale by lambda
    //before MBAR accumulation
    e *= (PMEDouble)(1.0 - SCterm);
#    ifdef LOCAL_MBAR
    // Set pieces of information that MBAR will need based on
    // this thread's results
    mbarTerm = e;
    mbarRefLambda = lambda;
#    endif
    e *= lambda;
#  endif
    df *= lambda;
#endif
    if (cSim.bJar) {
      double fold  = cSim.pNMRJarData[2];
      double work  = cSim.pNMRJarData[3];
      double first = cSim.pNMRJarData[4];
      double fcurr = (PMEDouble)-2.0 * K2K3.x * (rij - R1R2.y);
      if (first == (PMEDouble)0.0) {
	fold                = -fcurr;
	cSim.pNMRJarData[4] = (PMEDouble)1.0;
      }
      work                += (PMEDouble)0.5 * (fcurr + fold) * cSim.drjar;
      cSim.pNMRJarData[0]  = R1R2.y;
      cSim.pNMRJarData[1]  = rij;
      cSim.pNMRJarData[2]  = fcurr;
      cSim.pNMRJarData[3]  = work;
    }
    df *= (PMEDouble)1.0 / rij;
    PMEDouble fx = df * xij;
    PMEDouble fy = df * yij;
    PMEDouble fz = df * zij;
#ifdef use_DPFP
    PMEAccumulator ifx = llrint(fx * FORCESCALE);
    PMEAccumulator ify = llrint(fy * FORCESCALE);
    PMEAccumulator ifz = llrint(fz * FORCESCALE);
    atomicAdd((unsigned long long int*)&atmfrcx[atmJ], llitoulli(ifx));
    atomicAdd((unsigned long long int*)&atmfrcy[atmJ], llitoulli(ify));
    atomicAdd((unsigned long long int*)&atmfrcz[atmJ], llitoulli(ifz));
    atomicAdd((unsigned long long int*)&atmfrcx[atmI], llitoulli(-ifx));
    atomicAdd((unsigned long long int*)&atmfrcy[atmI], llitoulli(-ify));
    atomicAdd((unsigned long long int*)&atmfrcz[atmI], llitoulli(-ifz));
#else
    int ifx = __float2int_rn(fx * BSCALEF);
    int ify = __float2int_rn(fy * BSCALEF);
    int ifz = __float2int_rn(fz * BSCALEF);
    atomicAdd(&atmfrcx[atmJ],  ifx);
    atomicAdd(&atmfrcy[atmJ],  ify);
    atomicAdd(&atmfrcz[atmJ],  ifz);
    atomicAdd(&atmfrcx[atmI], -ifx);
    atomicAdd(&atmfrcy[atmI], -ify);
    atomicAdd(&atmfrcz[atmI], -ifz);
#endif
#ifdef LOCAL_ENERGY
    eterm = llrint(ENERGYSCALE * e);
  }
#  define EACC_OFFSET   NMR2_EACC_OFFSET
#  define STORE_SCR1    cSim.pESCNMRDistanceR1
#  define STORE_SCR2    cSim.pESCNMRDistanceR2
#include "kBWU_EnergyReduction.h"
#  undef EACC_OFFSET
#  undef STORE_SCR1
#  undef STORE_SCR2
#else
  }
#endif // LOCAL_ENERGY
}
