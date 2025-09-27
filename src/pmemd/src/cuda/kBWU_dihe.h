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
  unsigned int rawID = cSim.pBwuDiheID[startidx + tgx];
#ifdef LOCAL_ENERGY
  if (rawID == 0xffffffff) {
    eterm = (PMEAccumulator)0;
  }
  else {
#else
  if (rawID != 0xffffffff) {
#endif
    unsigned int atmI = rawID >> 24;
    unsigned int atmJ = (rawID >> 16) & 0xff;
    unsigned int atmK = (rawID >>  8) & 0xff;
    unsigned int atmL = rawID & 0xff;
    PMEFloat xij      = atmcrdx[atmI] - atmcrdx[atmJ];
    PMEFloat yij      = atmcrdy[atmI] - atmcrdy[atmJ];
    PMEFloat zij      = atmcrdz[atmI] - atmcrdz[atmJ];
    PMEFloat xkj      = atmcrdx[atmK] - atmcrdx[atmJ];
    PMEFloat ykj      = atmcrdy[atmK] - atmcrdy[atmJ];
    PMEFloat zkj      = atmcrdz[atmK] - atmcrdz[atmJ];
    PMEFloat xkl      = atmcrdx[atmK] - atmcrdx[atmL];
    PMEFloat ykl      = atmcrdy[atmK] - atmcrdy[atmL];
    PMEFloat zkl      = atmcrdz[atmK] - atmcrdz[atmL];
    PMEFloat2 dihe1   = cSim.pBwuDihe12[2*startidx + tgx];
    PMEFloat2 dihe2   = cSim.pBwuDihe12[2*startidx + GRID + tgx];
    PMEFloat  dihe3   = cSim.pBwuDihe3[startidx + tgx];

    // Get the normal vector
    PMEFloat dx      = yij*zkj - zij*ykj;
    PMEFloat dy      = zij*xkj - xij*zkj;
    PMEFloat dz      = xij*ykj - yij*xkj;
    PMEFloat gx      = zkj*ykl - ykj*zkl;
    PMEFloat gy      = xkj*zkl - zkj*xkl;
    PMEFloat gz      = ykj*xkl - xkj*ykl;
    PMEFloat fxi     = sqrt(dx*dx + dy*dy + dz*dz + tm24);
    PMEFloat fyi     = sqrt(gx*gx + gy*gy + gz*gz + tm24);
    PMEFloat ct      = dx*gx + dy*gy + dz*gz;

    // Branch if linear dihedral:
    PMEFloat z1     = (tenm3 <= fxi) ? (one / fxi) : zero;
    PMEFloat z2     = (tenm3 <= fyi) ? (one / fyi) : zero;
    PMEFloat z12    = z1 * z2;
    PMEFloat fzi    = (z12 != zero) ? one : zero;
    PMEFloat s      = xkj*(dz*gy - dy*gz) + ykj*(dx*gz - dz*gx) + zkj*(dy*gx - dx*gy);
    PMEFloat ap     = PI - (abs(acos(max(-one, min(one, ct * z12)))) *
			      (s >= (PMEFloat)0.0 ? (PMEFloat)1.0 : (PMEFloat)-1.0));
    PMEFloat sphi, cphi;
    #if defined(AMBER_PLATFORM_AMD) and defined(use_SPFP)
    __sincosf(ap, &sphi, &cphi);
    #else
    sincos(ap, &sphi, &cphi);
    #endif

    // Calculate the energy and the derivatives with respect to cosphi
    PMEFloat ct0 = dihe1.y * ap;
    PMEFloat sinnp, cosnp;
    #if defined(AMBER_PLATFORM_AMD) and defined(use_SPFP)
    __sincosf(ct0, &sinnp, &cosnp);
    #else
    sincos(ct0, &sinnp, &cosnp);
    #endif
    PMEFloat dums = sphi + tm24*(sphi >= (PMEFloat)0.0 ? (PMEFloat)1.0 : (PMEFloat)-1.0);
    PMEFloat df;
    if (tm06 > abs(dums)) {
      df = fzi * dihe2.y * (dihe1.y - dihe1.x + (dihe1.x * cphi));
    }
    else {
      df = fzi * dihe1.y * ((dihe2.y * sinnp) - dihe3*cosnp) / dums;
    }
#ifdef LOCAL_AFE
    unsigned int TIstatus = cSim.pBwuDiheStatus[startidx + tgx];
    TIregion = TIstatus >> 16;
    CVterm = TIstatus & 0xff;
    PMEFloat lambda = TISetLambda(CVterm, TIregion, cSim.AFElambdaSP[1]);
    df *= lambda;
#endif

    // Accelerated molecular dynamics: If iamd == 2 or 3, a dihedral boost will be added
    // to the dihedral forces, therefore df = df*fwgtd
    if ((cSim.iamd == 2) || (cSim.iamd == 3)) {
      PMEFloat fwgtd;
      fwgtd = cSim.pAMDfwgtd[0];
      df *= fwgtd;
    }

    // Gaussian accelerated MD: If igamd == 2, 3 or 5, a dihedral boost will be added to the
    // dihedral forces, therefore df = df*fwgtd
    // printf("Debug-GPU-kBWU_dihe) cSim.igamd, fwgtd = %12d, %12.5f\n", cSim.igamd, cSim.pGaMDfwgtd[0]);
    if ((cSim.igamd == 2) || (cSim.igamd == 3) || (cSim.igamd == 5)||(cSim.igamd == 19)) {
      PMEFloat fwgtd;
      fwgtd = cSim.pGaMDfwgtd[0];
      df *= fwgtd;
    }

    // Now do torsional first derivatives.  Set up array dc = 1st der. of cosphi
    // w/respect to cartesian differences:
    PMEFloat z11 = z1 * z1;
    z12          = z1 * z2;
    PMEFloat z22 = z2 * z2;
    PMEFloat dc1 = -gx*z12 - cphi*dx*z11;
    PMEFloat dc2 = -gy*z12 - cphi*dy*z11;
    PMEFloat dc3 = -gz*z12 - cphi*dz*z11;
    PMEFloat dc4 =  dx*z12 + cphi*gx*z22;
    PMEFloat dc5 =  dy*z12 + cphi*gy*z22;
    PMEFloat dc6 =  dz*z12 + cphi*gz*z22;

    // Update the first derivative array:
    PMEFloat dr1 = df * ( dc3*ykj - dc2*zkj);
    PMEFloat dr2 = df * ( dc1*zkj - dc3*xkj);
    PMEFloat dr3 = df * ( dc2*xkj - dc1*ykj);
    PMEFloat dr4 = df * ( dc6*ykj - dc5*zkj);
    PMEFloat dr5 = df * ( dc4*zkj - dc6*xkj);
    PMEFloat dr6 = df * ( dc5*xkj - dc4*ykj);
    PMEFloat drx = df * (-dc2*zij + dc3*yij + dc5*zkl - dc6*ykl);
    PMEFloat dry = df * ( dc1*zij - dc3*xij - dc4*zkl + dc6*xkl);
    PMEFloat drz = df * (-dc1*yij + dc2*xij + dc4*ykl - dc5*xkl);
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
    atomicAdd((unsigned long long int*)&atmfrcx[atmK], llitoulli( idrx + idr4));
    atomicAdd((unsigned long long int*)&atmfrcy[atmK], llitoulli( idry + idr5));
    atomicAdd((unsigned long long int*)&atmfrcz[atmK], llitoulli( idrz + idr6));
    atomicAdd((unsigned long long int*)&atmfrcx[atmL], llitoulli(-idr4));
    atomicAdd((unsigned long long int*)&atmfrcy[atmL], llitoulli(-idr5));
    atomicAdd((unsigned long long int*)&atmfrcz[atmL], llitoulli(-idr6));
#else
    int idr1 = __float2int_rn(dr1 * BSCALEF);
    int idr2 = __float2int_rn(dr2 * BSCALEF);
    int idr3 = __float2int_rn(dr3 * BSCALEF);
    int idr4 = __float2int_rn(dr4 * BSCALEF);
    int idr5 = __float2int_rn(dr5 * BSCALEF);
    int idr6 = __float2int_rn(dr6 * BSCALEF);
    int idrx = __float2int_rn(drx * BSCALEF);
    int idry = __float2int_rn(dry * BSCALEF);
    int idrz = __float2int_rn(drz * BSCALEF);
    atomicAdd(&atmfrcx[atmI], -idr1);
    atomicAdd(&atmfrcy[atmI], -idr2);
    atomicAdd(&atmfrcz[atmI], -idr3);
    atomicAdd(&atmfrcx[atmJ], -idrx + idr1);
    atomicAdd(&atmfrcy[atmJ], -idry + idr2);
    atomicAdd(&atmfrcz[atmJ], -idrz + idr3);
    atomicAdd(&atmfrcx[atmK],  idrx + idr4);
    atomicAdd(&atmfrcy[atmK],  idry + idr5);
    atomicAdd(&atmfrcz[atmK],  idrz + idr6);
    atomicAdd(&atmfrcx[atmL], -idr4);
    atomicAdd(&atmfrcy[atmL], -idr5);
    atomicAdd(&atmfrcz[atmL], -idr6);
#endif
#ifdef LOCAL_ENERGY
#  ifdef LOCAL_AFE
    SCterm = (TIstatus >> 8) & 0xff;
    if (SCterm == 0) {
      eterm = llrint(ENERGYSCALEF * lambda * (dihe2.x + cosnp*dihe2.y + sinnp*dihe3) * fzi);
    }
    else {
      scEcomp[TIregion - 1] = llrint(ENERGYSCALE * (dihe2.x + cosnp*dihe2.y +
                                                    sinnp*dihe3) * fzi);
      eterm = (PMEAccumulator)0;
    }
    edvdl = llrint((PMEDouble)(CVterm * (2*TIregion - 3)) * ENERGYSCALE *
                   (dihe2.x + cosnp*dihe2.y + sinnp*dihe3) * fzi);
#    ifdef LOCAL_MBAR
    // Set pieces of information that MBAR will need based on this thread's results
    mbarTerm = (dihe2.x + cosnp*dihe2.y + sinnp*dihe3) * fzi;
    mbarRefLambda = lambda;
#    endif
#  else
    eterm = llrint(ENERGYSCALEF * (dihe2.x + cosnp*dihe2.y + sinnp*dihe3) * fzi);
#  endif
  }
#  define EACC_OFFSET   DIHE_EACC_OFFSET
#  define STORE_SCR1    cSim.pSCDihedralR1
#  define STORE_SCR2    cSim.pSCDihedralR2
#include "kBWU_EnergyReduction.h"
#  undef EACC_OFFSET
#  undef STORE_SCR1
#  undef STORE_SCR2
#else
  }
#endif // LOCAL_ENERGY
}
