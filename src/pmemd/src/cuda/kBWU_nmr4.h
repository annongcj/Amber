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
  unsigned int rawID = cSim.pBwuNMR4ID[offset];
  int Inc       = cSim.pBwuNMR4ID[offset +   GRID];
  int stepInit  = cSim.pBwuNMR4ID[offset + 2*GRID];
  int stepFinal = cSim.pBwuNMR4ID[offset + 3*GRID];
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
    PMEDouble2 R1R2Slope    = cSim.pBwuNMR4[offset         ];
    PMEDouble2 R3R4Slope    = cSim.pBwuNMR4[offset +   GRID];
    PMEDouble2 K2K3Slope    = cSim.pBwuNMR4[offset + 2*GRID];
    PMEDouble2 R1R2Intrcpt  = cSim.pBwuNMR4[offset + 3*GRID];
    PMEDouble2 R3R4Intrcpt  = cSim.pBwuNMR4[offset + 4*GRID];
    PMEDouble2 K2K3Intrcpt  = cSim.pBwuNMR4[offset + 5*GRID];

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

    // Compute displacements
    unsigned int atmI = rawID >> 24;
    unsigned int atmJ = (rawID >> 16) & 0xff;
    unsigned int atmK = (rawID >> 8) & 0xff;
    unsigned int atmL = rawID & 0xff;
    PMEDouble xij = atmcrdx[atmI] - atmcrdx[atmJ];
    PMEDouble yij = atmcrdy[atmI] - atmcrdy[atmJ];
    PMEDouble zij = atmcrdz[atmI] - atmcrdz[atmJ];
    PMEDouble xkj = atmcrdx[atmK] - atmcrdx[atmJ];
    PMEDouble ykj = atmcrdy[atmK] - atmcrdy[atmJ];
    PMEDouble zkj = atmcrdz[atmK] - atmcrdz[atmJ];
    PMEDouble xkl = atmcrdx[atmK] - atmcrdx[atmL];
    PMEDouble ykl = atmcrdy[atmK] - atmcrdy[atmL];
    PMEDouble zkl = atmcrdz[atmK] - atmcrdz[atmL];

    // Calculate ij X jk AND kl X jk:
    PMEDouble dx = yij*zkj - zij*ykj;
    PMEDouble dy = zij*xkj - xij*zkj;
    PMEDouble dz = xij*ykj - yij*xkj;
    PMEDouble gx = zkj*ykl - ykj*zkl;
    PMEDouble gy = xkj*zkl - zkj*xkl;
    PMEDouble gz = ykj*xkl - xkj*ykl;

    // Calculate the magnitudes of above vectors, and their dot product:
    PMEDouble fxi = dx*dx + dy*dy + dz*dz + tm24;
    PMEDouble fyi = gx*gx + gy*gy + gz*gz + tm24;
    PMEDouble ct  = dx*gx + dy*gy + dz*gz;

    // Branch if linear dihedral:
    PMEDouble z1  = rsqrt(fxi);
    PMEDouble z2  = rsqrt(fyi);
    PMEDouble z11 = z1 * z1;
    PMEDouble z22 = z2 * z2;
    PMEDouble z12 = z1 * z2;
    ct           *= z1 * z2;
    ct            = max((PMEDouble)-0.9999999999999, min(ct, (PMEDouble)0.9999999999999));
    PMEDouble ap  = acos(ct);
    PMEDouble s   = xkj*(dz*gy - dy*gz) + ykj*(dx*gz - dz*gx) + zkj*(dy*gx - dx*gy);
    if (s < (PMEDouble)0.0) {
      ap = -ap;
    }
    ap = PI - ap;
    PMEDouble sphi, cphi;
    faster_sincos(ap, &sphi, &cphi);

    // Translate the value of the torsion (by +- n*360) to bring it as close as
    // possible to one of the two central "cutoff" points (r2,r3). Use this as
    // the value of the torsion in the following comparison.
    PMEDouble apmean = (R1R2.y + R3R4.x) * (PMEDouble)0.5;
    if (ap - apmean > PI) {
      ap -= (PMEDouble)2.0 * (PMEDouble)PI;
    }
    if (apmean - ap > PI) {
      ap += (PMEDouble)2.0 * (PMEDouble)PI;
    }
    PMEDouble df;
#ifdef LOCAL_ENERGY
    PMEDouble e;
#endif
    if (ap < R1R2.x) {
      PMEDouble dif = R1R2.x - R1R2.y;
      df            = (PMEDouble)2.0 * K2K3.x * dif;
#ifdef LOCAL_ENERGY
      e             = df*(ap - R1R2.x) + K2K3.x*dif*dif;
#endif
    }
    else if (ap < R1R2.y) {
      PMEDouble dif = ap - R1R2.y;
      df            = (PMEDouble)2.0 * K2K3.x * dif;
#ifdef LOCAL_ENERGY
      e             = K2K3.x * dif * dif;
#endif
    }
    else if (ap < R3R4.x) {
      df            = (PMEDouble)0.0;
#ifdef LOCAL_ENERGY
      e             = (PMEDouble)0.0;
#endif
    }
    else if (ap < R3R4.y) {
      PMEDouble dif = ap - R3R4.x;
      df            = (PMEDouble)2.0 * K2K3.y * dif;
#ifdef LOCAL_ENERGY
      e             = K2K3.y * dif * dif;
#endif
    }
    else {
      PMEDouble dif = R3R4.y - R3R4.x;
      df            = (PMEDouble)2.0 * K2K3.y * dif;
#ifdef LOCAL_ENERGY
      e             = df*(ap - R3R4.y) + K2K3.y*dif*dif;
#endif
    }
#ifdef LOCAL_AFE
    unsigned int TIstatus = cSim.pBwuNMR4Status[startidx + tgx];
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
    e *= lambda * (PMEDouble)(1.0 - SCterm);
#  endif
    df *= lambda;
#endif
    df             *= -1.0 / sphi;
    if (cSim.bJar) {
      double fold   = cSim.pNMRJarData[2];
      double work   = cSim.pNMRJarData[3];
      double first  = cSim.pNMRJarData[4];
      double fcurr  = (PMEDouble)-2.0 * K2K3.x * (ap - R1R2.y);
      if (first == 0.0) {
	fold                = -fcurr;
	cSim.pNMRJarData[4] = (PMEDouble)1.0;
      }
      work                += (PMEDouble)0.5 * (fcurr + fold) * cSim.drjar;
      cSim.pNMRJarData[0]  = R1R2.y;
      cSim.pNMRJarData[1]  = ap;
      cSim.pNMRJarData[2]  = fcurr;
      cSim.pNMRJarData[3]  = work;
    }

    // dc = First derivative of cos(phi) w/respect to the cartesian differences t.
    PMEDouble dcdx = -gx*z12 - cphi*dx*z11;
    PMEDouble dcdy = -gy*z12 - cphi*dy*z11;
    PMEDouble dcdz = -gz*z12 - cphi*dz*z11;
    PMEDouble dcgx =  dx*z12 + cphi*gx*z22;
    PMEDouble dcgy =  dy*z12 + cphi*gy*z22;
    PMEDouble dcgz =  dz*z12 + cphi*gz*z22;
    PMEDouble dr1  =  df*( dcdz*ykj - dcdy*zkj);
    PMEDouble dr2  =  df*( dcdx*zkj - dcdz*xkj);
    PMEDouble dr3  =  df*( dcdy*xkj - dcdx*ykj);
    PMEDouble dr4  =  df*( dcgz*ykj - dcgy*zkj);
    PMEDouble dr5  =  df*( dcgx*zkj - dcgz*xkj);
    PMEDouble dr6  =  df*( dcgy*xkj - dcgx*ykj);
    PMEDouble drx  =  df*(-dcdy*zij + dcdz*yij + dcgy*zkl -  dcgz*ykl);
    PMEDouble dry  =  df*( dcdx*zij - dcdz*xij - dcgx*zkl +  dcgz*xkl);
    PMEDouble drz  =  df*(-dcdx*yij + dcdy*xij + dcgx*ykl -  dcgy*xkl);
#ifdef use_DPFP
    PMEAccumulator idr1 = llrint(dr1 * FORCESCALE);
    PMEAccumulator idr2 = llrint(dr2 * FORCESCALE);
    PMEAccumulator idr3 = llrint(dr3 * FORCESCALE);
    PMEAccumulator idr4 = llrint(dr4 * FORCESCALE);
    PMEAccumulator idr5 = llrint(dr5 * FORCESCALE);
    PMEAccumulator idr6 = llrint(dr6 * FORCESCALE);
    PMEAccumulator idrx = llrint(drx * FORCESCALE);
    PMEAccumulator idry = llrint(dry * FORCESCALE);
    PMEAccumulator idrz = llrint(drz * FORCESCALE);
    atomicAdd((unsigned long long int*)&atmfrcx[atmI], llitoulli(-idr1));
    atomicAdd((unsigned long long int*)&atmfrcy[atmI], llitoulli(-idr2));
    atomicAdd((unsigned long long int*)&atmfrcz[atmI], llitoulli(-idr3));
    atomicAdd((unsigned long long int*)&atmfrcx[atmJ], llitoulli(-idrx + idr1));
    atomicAdd((unsigned long long int*)&atmfrcy[atmJ], llitoulli(-idry + idr2));
    atomicAdd((unsigned long long int*)&atmfrcz[atmJ], llitoulli(-idrz + idr3));
    atomicAdd((unsigned long long int*)&atmfrcx[atmK], llitoulli(idrx + idr4));
    atomicAdd((unsigned long long int*)&atmfrcy[atmK], llitoulli(idry + idr5));
    atomicAdd((unsigned long long int*)&atmfrcz[atmK], llitoulli(idrz + idr6));
    atomicAdd((unsigned long long int*)&atmfrcx[atmL], llitoulli(-idr4));
    atomicAdd((unsigned long long int*)&atmfrcy[atmL], llitoulli(-idr5));
    atomicAdd((unsigned long long int*)&atmfrcz[atmL], llitoulli(-idr6));
#else
    int idr1 = dr1 * BSCALEF;
    int idr2 = dr2 * BSCALEF;
    int idr3 = dr3 * BSCALEF;
    int idr4 = dr4 * BSCALEF;
    int idr5 = dr5 * BSCALEF;
    int idr6 = dr6 * BSCALEF;
    int idrx = drx * BSCALEF;
    int idry = dry * BSCALEF;
    int idrz = drz * BSCALEF;
    atomicAdd(&atmfrcx[atmI], -idr1);
    atomicAdd(&atmfrcy[atmI], -idr2);
    atomicAdd(&atmfrcz[atmI], -idr3);
    atomicAdd(&atmfrcx[atmJ], -idrx + idr1);
    atomicAdd(&atmfrcy[atmJ], -idry + idr2);
    atomicAdd(&atmfrcz[atmJ], -idrz + idr3);
    atomicAdd(&atmfrcx[atmK], idrx + idr4);
    atomicAdd(&atmfrcy[atmK], idry + idr5);
    atomicAdd(&atmfrcz[atmK], idrz + idr6);
    atomicAdd(&atmfrcx[atmL], -idr4);
    atomicAdd(&atmfrcy[atmL], -idr5);
    atomicAdd(&atmfrcz[atmL], -idr6);
#endif
#ifdef LOCAL_ENERGY
    eterm = llrint(ENERGYSCALE * e);
  }
#  define EACC_OFFSET   NMR4_EACC_OFFSET
#  define STORE_SCR1    cSim.pESCNMRTorsionR1
#  define STORE_SCR2    cSim.pESCNMRTorsionR2
#include "kBWU_EnergyReduction.h"
#  undef EACC_OFFSET
#  undef STORE_SCR1
#  undef STORE_SCR2
#else
  }
#endif // LOCAL_ENERGY
}
