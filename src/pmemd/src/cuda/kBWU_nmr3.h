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
  unsigned int rawID = cSim.pBwuNMR3ID[offset];
  int Inc       = cSim.pBwuNMR3ID[offset +   GRID];
  int stepInit  = cSim.pBwuNMR3ID[offset + 2*GRID];
  int stepFinal = cSim.pBwuNMR3ID[offset + 3*GRID];
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
    PMEDouble2 R1R2Slope    = cSim.pBwuNMR3[offset         ];
    PMEDouble2 R3R4Slope    = cSim.pBwuNMR3[offset +   GRID];
    PMEDouble2 K2K3Slope    = cSim.pBwuNMR3[offset + 2*GRID];
    PMEDouble2 R1R2Intrcpt  = cSim.pBwuNMR3[offset + 3*GRID];
    PMEDouble2 R3R4Intrcpt  = cSim.pBwuNMR3[offset + 4*GRID];
    PMEDouble2 K2K3Intrcpt  = cSim.pBwuNMR3[offset + 5*GRID];

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

    // Compute displacements, force, and energy
    unsigned int atmI = rawID >> 16;
    unsigned int atmJ = (rawID >> 8) & 0xff;
    unsigned int atmK = rawID & 0xff;
    PMEDouble xij     = atmcrdx[atmI] - atmcrdx[atmJ];
    PMEDouble yij     = atmcrdy[atmI] - atmcrdy[atmJ];
    PMEDouble zij     = atmcrdz[atmI] - atmcrdz[atmJ];
    PMEDouble xkj     = atmcrdx[atmK] - atmcrdx[atmJ];
    PMEDouble ykj     = atmcrdy[atmK] - atmcrdy[atmJ];
    PMEDouble zkj     = atmcrdz[atmK] - atmcrdz[atmJ];
    PMEDouble rij2    = xij*xij + yij*yij + zij*zij;
    PMEDouble rkj2    = xkj*xkj + ykj*ykj + zkj*zkj;
    PMEDouble rij     = sqrt(rij2);
    PMEDouble rkj     = sqrt(rkj2);
    PMEDouble rdenom  = rij * rkj;
    PMEDouble cst     = min(1.0, max(-1.0, (xij*xkj + yij*ykj + zij*zkj) / rdenom));
    PMEDouble theta   = acos(cst);
    PMEDouble df;
#ifdef LOCAL_ENERGY
    PMEDouble e;
#endif
    if (theta < R1R2.x) {
      PMEDouble dif = R1R2.x - R1R2.y;
      df            = (PMEDouble)2.0 * K2K3.x * dif;
#ifdef LOCAL_ENERGY
      e             = df*(theta - R1R2.x) + K2K3.x*dif*dif;
#endif
    }
    else if (theta < R1R2.y) {
      PMEDouble dif = theta - R1R2.y;
      df            = (PMEDouble)2.0 * K2K3.x * dif;
#ifdef LOCAL_ENERGY
      e             = K2K3.x * dif * dif;
#endif
    }
    else if (theta < R3R4.x) {
      df            = (PMEDouble)0.0;
#ifdef LOCAL_ENERGY
      e             = (PMEDouble)0.0;
#endif
    }
    else if (theta < R3R4.y) {
      PMEDouble dif = theta - R3R4.x;
      df            = (PMEDouble)2.0 * K2K3.y * dif;
#ifdef LOCAL_ENERGY
      e             = K2K3.y * dif * dif;
#endif
    }
    else {
      PMEDouble dif = R3R4.y - R3R4.x;
      df            = (PMEDouble)2.0 * K2K3.y * dif;
#ifdef LOCAL_ENERGY
      e             = df * (theta - R3R4.y) + K2K3.y * dif * dif;
#endif
    }
#ifdef LOCAL_AFE
    unsigned int TIstatus = cSim.pBwuNMR3Status[startidx + tgx];
    TIregion = TIstatus >> 16;
    CVterm = TIstatus & 0xff;
    PMEDouble lambda = TISetLambda(CVterm, TIregion, cSim.AFElambdaSP[1]);
#  ifdef LOCAL_ENERGY
    SCterm = (TIstatus >> 8) & 0xff;
    scEcomp[TIregion - SCterm] = llrint(ENERGYSCALE * e * (PMEDouble)SCterm);
    edvdl = llrint((PMEDouble)(CVterm * (2*TIregion - 3)) * ENERGYSCALE * e);
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
      double fold   = cSim.pNMRJarData[2];
      double work   = cSim.pNMRJarData[3];
      double first  = cSim.pNMRJarData[4];
      double fcurr  = (PMEDouble)-2.0 * K2K3.x * (theta - R1R2.y);
      if (first == (PMEDouble)0.0) {
	fold                = -fcurr;
	cSim.pNMRJarData[4] = (PMEDouble)1.0;
      }
      work                += (PMEDouble)0.5 * (fcurr + fold) * cSim.drjar;
      cSim.pNMRJarData[0]  = R1R2.y;
      cSim.pNMRJarData[1]  = theta;
      cSim.pNMRJarData[2]  = fcurr;
      cSim.pNMRJarData[3]  = work;
    }
    PMEDouble snt  = sin(theta);
    if (abs(snt) < (PMEDouble)1.0e-14) {
      snt = (PMEDouble)1.0e-14;
    }
    PMEDouble st  = -df / snt;
    PMEDouble sth = st * cst;
    PMEDouble cik = st / rdenom;
    PMEDouble cii = sth / rij2;
    PMEDouble ckk = sth / rkj2;
    PMEDouble fx1 = cii*xij - cik*xkj;
    PMEDouble fy1 = cii*yij - cik*ykj;
    PMEDouble fz1 = cii*zij - cik*zkj;
    PMEDouble fx2 = ckk*xkj - cik*xij;
    PMEDouble fy2 = ckk*ykj - cik*yij;
    PMEDouble fz2 = ckk*zkj - cik*zij;
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
    atomicAdd((unsigned long long int*)&atmfrcx[atmJ], llitoulli(-ifx1 - ifx2));
    atomicAdd((unsigned long long int*)&atmfrcy[atmJ], llitoulli(-ify1 - ify2));
    atomicAdd((unsigned long long int*)&atmfrcz[atmJ], llitoulli(-ifz1 - ifz2));
#else
    int ifx1 = fx1 * BSCALEF;
    int ify1 = fy1 * BSCALEF;
    int ifz1 = fz1 * BSCALEF;
    int ifx2 = fx2 * BSCALEF;
    int ify2 = fy2 * BSCALEF;
    int ifz2 = fz2 * BSCALEF;
    atomicAdd(&atmfrcx[atmI], llitoulli(ifx1));
    atomicAdd(&atmfrcy[atmI], llitoulli(ify1));
    atomicAdd(&atmfrcz[atmI], llitoulli(ifz1));
    atomicAdd(&atmfrcx[atmK], llitoulli(ifx2));
    atomicAdd(&atmfrcy[atmK], llitoulli(ify2));
    atomicAdd(&atmfrcz[atmK], llitoulli(ifz2));
    atomicAdd(&atmfrcx[atmJ], llitoulli(-ifx1 - ifx2));
    atomicAdd(&atmfrcy[atmJ], llitoulli(-ify1 - ify2));
    atomicAdd(&atmfrcz[atmJ], llitoulli(-ifz1 - ifz2));
#endif
#ifdef LOCAL_ENERGY
    eterm = llrint(ENERGYSCALE * e);
  }
#  define EACC_OFFSET   NMR3_EACC_OFFSET
#  define STORE_SCR1    cSim.pESCNMRAngleR1
#  define STORE_SCR2    cSim.pESCNMRAngleR2
#include "kBWU_EnergyReduction.h"
#  undef EACC_OFFSET
#  undef STORE_SCR1
#  undef STORE_SCR2
#else
  }
#endif // LOCAL_ENERGY
}
