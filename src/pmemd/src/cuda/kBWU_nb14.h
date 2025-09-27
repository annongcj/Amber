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
  unsigned int rawID = cSim.pBwuNB14ID[startidx + tgx];
#ifdef LOCAL_VIRIAL
  PMEAccumulator v11 = (PMEAccumulator)0;
  PMEAccumulator v22 = (PMEAccumulator)0;
  PMEAccumulator v33 = (PMEAccumulator)0;
#endif
#ifdef LOCAL_ENERGY
  // In this special case, there are two accumulators to track: the accumulator inherited from
  // kBWU.h, scEcomp, will take the role of tracking Lennard-Jones soft-core interactions
  // (nb14), while a new accumulator, scEel14comp, will handle electrostatics 1-4 soft-core
  // interactions.  Similarly, eel14 will handle the overall accumulation of electrostatic
  // energy while the inherited eterm goes to Lennard-Jones interactions.
  PMEAccumulator eel14;
#  ifdef LOCAL_AFE
  PMEAccumulator scEel14comp[2] = {(PMEAccumulator)0, (PMEAccumulator)0};
#  endif
  if (rawID == 0xffffffff) {
    eterm = (PMEAccumulator)0;
    eel14 = (PMEAccumulator)0;
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
    PMEFloat2 scnb  = cSim.pBwuLJnb14[startidx + tgx];
#ifdef PHMD
    unsigned int first = indices[atmI];
    unsigned int second = indices[atmJ];
    PMEFloat qi = charges[atmI];
    PMEFloat qj = charges[atmJ];
    PMEFloat scee = cSim.pBwuEEnb14[startidx + tgx];
    PMEFloat qterm  = scee * qi * qj;
#else
    PMEFloat qterm  = cSim.pBwuEEnb14[startidx + tgx];
#endif
    PMEFloat r2     = dx*dx + dy*dy + dz*dz;
    PMEFloat rinv   = rsqrt(r2);
    PMEFloat r2inv  = rinv * rinv;
    PMEFloat r6     = r2inv * r2inv * r2inv;
    PMEFloat g      = rinv * qterm;
    PMEFloat f6     = scnb.y * r6;
    PMEFloat f12    = scnb.x * (r6 * r6);
    PMEFloat df     = (g + (PMEFloat)12.0*f12 - (PMEFloat)6.0*f6) * r2inv;
#ifdef PHMD
    PMEFloat x2h, x2g, qunprot, qprot, qxh, qh, qxg, qg, radh, lambdah, lambdah2, facth, factg, dxg, dxh, radg, lambdag, lambdag2, dudlh, dudlg, dudlhplus, dudlgplus;
    PMEFloat hpcharge_phmd, gpcharge_phmd, intdiel_inv;
    PMEFloat2 hpqstate1, hpqstate2, gpqstate1, gpqstate2, hpvstate1, gpvstate1;
    PMEFloat fact = 1.0;
#ifdef LOCAL_NEIGHBORLIST    
    int h = cSim.pImageGrplist[first] - 1;
    int k = cSim.pImageGrplist[second] - 1;
    hpqstate1 = cSim.pImageQstate1[first];
    hpqstate2 = cSim.pImageQstate2[first];
    gpqstate1 = cSim.pImageQstate1[second];
    gpqstate2 = cSim.pImageQstate2[second];
    hpcharge_phmd = cSim.pImageCharge_phmd[first];
    gpcharge_phmd = cSim.pImageCharge_phmd[second];
    hpvstate1 = cSim.pvstate1[cSim.pImageIndexPHMD[first]];
    gpvstate1 = cSim.pvstate1[cSim.pImageIndexPHMD[second]];
#else
    int h = cSim.pgrplist[first] - 1;
    int k = cSim.pgrplist[second] - 1;
    hpqstate1 = cSim.pqstate1[first];
    hpqstate2 = cSim.pqstate2[first];
    gpqstate1 = cSim.pqstate1[second];
    gpqstate2 = cSim.pqstate2[second];
    hpcharge_phmd = cSim.pcharge_phmd[first];
    gpcharge_phmd = cSim.pcharge_phmd[second];
    hpvstate1 = cSim.pvstate1[first];
    gpvstate1 = cSim.pvstate1[second];
#endif    
    int spgrp1 = cSim.psp_grp[h];
    int spgrp2 = cSim.psp_grp[k];
    intdiel_inv = cSim.intdiel_inv;
    dudlh = 0;
    dudlg = 0;
    dudlhplus = 0;
    dudlgplus = 0;
    if(h >= 0) {
      x2h = 1.0;
      if(spgrp1 > 0) {

        x2h = sin(cSim.pph_theta[h+1]);
        x2h *= x2h;
        lambdah = sin(cSim.pph_theta[h]);
        lambdah *= lambdah;
        qunprot = lambdah * (hpqstate2.x - hpqstate2.y);
        qprot = (1 - lambdah) * (hpqstate1.x - hpqstate1.y);
        qxh = gpcharge_phmd * (qunprot + qprot);
      }
      qunprot = x2h * hpqstate2.x + (1.0 - x2h) * hpqstate2.y;
      qprot = x2h * hpqstate1.x + (1.0 - x2h) * hpqstate1.y;
      qh = gpcharge_phmd * (qunprot - qprot);
      dudlh += qh * rinv * intdiel_inv * scee;
      if(spgrp1 > 0) {
        dudlhplus += qxh * rinv * intdiel_inv * scee;
      }
    }
    if(k >= 0) {
      x2g = 1.0;
      if(spgrp2 > 0) {
        x2g = sin(cSim.pph_theta[k+1]);
        x2g *= x2g;
        lambdag = sin(cSim.pph_theta[k]);
        lambdag *= lambdag;
        qunprot = lambdag * (gpqstate2.x - gpqstate2.y);
        qprot = (1.0 - lambdag) * (gpqstate1.x - gpqstate1.y);
        qxg = hpcharge_phmd * (qunprot + qprot);
      }
      qunprot = x2g * gpqstate2.x + (1.0 - x2g) * gpqstate2.y;
      qprot = x2g * gpqstate1.x + (1.0 - x2g) * gpqstate1.y;
      qg = hpcharge_phmd * (qunprot - qprot);
      dudlg += qg * rinv * intdiel_inv * scee;
      if(spgrp2 > 0) {
        dudlgplus += qxg * rinv * intdiel_inv * scee;
      }
    }
    bool lhtitr = false;
    bool lhtauto = false;
    PMEFloat xh = 1.0;
    int ihtaut = 0;
    if(h >= 0) {
      radh = hpvstate1.x;
      lambdah = sin(cSim.pph_theta[h]);
      lambdah *= lambdah;
      lambdah2 = 1.0 - lambdah;
      facth = lambdah2;
      if(spgrp1 == 1 || spgrp1 == 3) {
        lhtauto = true;
        radh += hpvstate1.y;
        if(hpvstate1.x > 0) {
          ihtaut = 1;
          xh = x2h;
        }
        else if(hpvstate1.y > 0) {
          ihtaut = -1;
          xh = 1.0 - x2h;
        }
        if(spgrp1 == 1) {
          facth = 1.0 - lambdah * xh;
          dxh = -ihtaut * lambdah;
        }
        else if(spgrp1 == 3) {
          facth = lambdah2 * xh;
          dxh = ihtaut * lambdah2;
        }
      }
      if(radh > 0) {
         lhtitr = true;
      }
    }
    bool lgtitr = false;
    bool lgtauto = false;
    PMEFloat xg = 1.0;
    int igtaut = 0;
    if(k >= 0) {
       radg = gpvstate1.x;
       lambdag = sin(cSim.pph_theta[k]);
       lambdag *= lambdag;
       lambdag2 = 1.0 - lambdag;
       factg = lambdag2;
       if(spgrp2 == 1 || spgrp2 == 3) {
         lgtauto = true;
         radg += gpvstate1.y;
         if(gpvstate1.x > 0) {
           igtaut = 1;
           xg = x2g;
         }
         else if(gpvstate1.y > 0) {
           igtaut = -1;
           xg = 1.0 - x2g;
         }
         if(spgrp2 == 1) {
           factg = 1.0 - lambdag * xg;
           dxg = -igtaut * lambdag;
         }
         else if(spgrp2 == 3) {
           factg = lambdag2 * xg;
           dxg = igtaut * lambdag2;
         }
       }
       if(radg > 0) {
         lgtitr = true;
       }
    }
    if(lhtitr && lgtitr) {
      if(h != k) {
        fact = facth * factg;
        dudlh += -xh * factg * (f12 - f6);
        dudlg += -xg * facth * (f12 - f6);
        if(lhtauto) {
          dudlhplus += dxh * factg * (f12 - f6);
        }
        if(lgtauto) {
          dudlgplus += dxg * facth * (f12 - f6);
        }
      }
      else if(spgrp1 == 1) {
        fact = 1.0 - lambdah;
        dudlh += -(f12 - f6);
      }
    }
    else if(lhtitr) {
      fact = facth;
      dudlh += -xh * (f12 - f6);
      if(lhtauto) {
        dudlhplus += dxh * (f12 - f6);
      }
    }
    else if(lgtitr) {
      fact = factg;
      dudlg += -xg * (f12 - f6);
      if(lgtauto) {
        dudlgplus += dxg * (f12 - f6);
      }
    }
    df     = (g + fact * (PMEFloat)12.0*f12 - fact * (PMEFloat)6.0*f6) * r2inv;
    #else
    df     = (g + (PMEFloat)12.0*f12 - (PMEFloat)6.0*f6) * r2inv;
    #endif
#ifdef LOCAL_AFE
    unsigned int TIstatus = cSim.pBwuNB14Status[startidx + tgx];
    TIregion = TIstatus >> 16;
    CVterm = TIstatus & 0xff;
    PMEFloat lambda = TISetLambda(CVterm, TIregion, cSim.AFElambdaSP[1]);
    df *= lambda;
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
  #ifdef PHMD
    if(h >= 0) {
      atomicAdd((unsigned long long int*)&atmdudl[atmI],
                llitoulli(llrint(dudlh * FORCESCALE)));
      if(cSim.psp_grp[h] > 0) {
        atomicAdd((unsigned long long int*)&atmdudlplus[atmI],
                  llitoulli(llrint(dudlhplus * FORCESCALE)));
      }
    }
    if(k >= 0) {
      atomicAdd((unsigned long long int*)&atmdudl[atmJ],
                llitoulli(llrint(dudlg * FORCESCALE)));
      if(cSim.psp_grp[k] > 0) {
        atomicAdd((unsigned long long int*)&atmdudlplus[atmJ],
                  llitoulli(llrint(dudlgplus * FORCESCALE)));
      }
    }
  #endif
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
  #ifdef PHMD
    if(h >= 0) {
      atomicAdd(&atmdudl[atmI], __float2int_rn(dudlh * BSCALEF));
      if(cSim.psp_grp[h] > 0) {
        atomicAdd(&atmdudlplus[atmI], __float2int_rn(dudlhplus * BSCALEF));
      }
    }
    if(k >= 0) {
      atomicAdd(&atmdudl[atmJ], __float2int_rn(dudlg * BSCALEF));
      if(cSim.psp_grp[k] > 0) {
        atomicAdd(&atmdudlplus[atmJ], __float2int_rn(dudlgplus * BSCALEF));
      }
    }
  #endif
#endif // use_DPFP
#ifdef LOCAL_ENERGY
#  ifdef LOCAL_AFE
    SCterm = (TIstatus >> 8) & 0xff;
    if (SCterm == 0) {
      eel14 = llrint(ENERGYSCALE * lambda * g);
      eterm = llrint(ENERGYSCALE * lambda * (f12 - f6));
    }
    else {
      scEel14comp[TIregion - 1] = llrint(ENERGYSCALE * g);
      scEcomp[TIregion - 1]     = llrint(ENERGYSCALE * (f12 - f6));
      eel14 = (PMEAccumulator)0;
      eterm = (PMEAccumulator)0;
    }
    edvdl = llrint((PMEDouble)(CVterm * (2*TIregion - 3)) * ENERGYSCALE * (g + (f12 - f6)));
#    ifdef LOCAL_MBAR
    // Set pieces of information that MBAR will need based on this thread's results
    mbarTerm = g + (f12 - f6);
    mbarRefLambda = lambda;
#    endif
#  else
    eel14 = llrint(ENERGYSCALE * g);
    eterm = llrint(ENERGYSCALE * (f12 - f6));
#  endif
  }
#  define NB14_CASE
#  define EACC_OFFSET   SCNB_EACC_OFFSET
#  define STORE_SCR1    cSim.pSCVDW14R1
#  define STORE_SCR2    cSim.pSCVDW14R2
#include "kBWU_EnergyReduction.h"
#  undef EACC_OFFSET
#  undef STORE_SCR1
#  undef STORE_SCR2
#  undef NB14_CASE
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
