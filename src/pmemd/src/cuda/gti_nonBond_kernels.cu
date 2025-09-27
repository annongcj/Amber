#ifdef GTI

#ifdef AMBER_PLATFORM_AMD
#include <hip/hip_runtime.h>
#endif

#include "gti_def.h"
#include "gti_const.cuh"
#include "gti_utils.cuh"
#include "gpuContext.h"
#include "ptxmacros.h"
#include "gti_schedule_functions.h"
#include "gti_nonBond_kernels.cuh"
#include "simulationConst.h"

CSIM_STO simulationConst cSim;

namespace GTI_NB_CALC {
  #if !(defined(__CUDACC_RDC__) || defined(__HIPCC_RDC__))
    #include "gti_localCUutils.inc"
  #endif
}

__forceinline__ __device__ PMEFloat test_alpha(PMEFloat q01) {
  return  max(0.01, (q01 / 332.0 + 0.8) / 3.0 + 0.05);
};

using namespace GTI_NB_CALC;
__device__ const PMEDouble delta = 1.0e-8;
//---------------------------------------------------------------------------------------------
// kgCalculateTINB_kernel:
//
// Arguments:
//   needEnergy:
//   needVirial:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgCalculateTINB_kernel(bool needEnergy, bool needVirial) {

  unsigned i = threadIdx.x;
  unsigned j = threadIdx.x;
  __shared__ PMEFloat recip[9], ucell[9];

  //PMEFloat4* cn[MaxLJIndex], sigEps[MaxLJIndex];
  PMEFloat4 *cn, *sigEps;
  PMEFloat *dvalue; //C4PairwiseCUDA2023
  int *dcoef; //C4PairwiseCUDA

  __shared__ PMEFloat4 cn0[MaxLJIndex];
  __shared__ PMEFloat4 sigEps0[MaxLJIndex]; 
  __shared__ int dcoef0[100];  //C4PairwiseCUDA2023
  __shared__ PMEFloat dvalue0[100]; //C4PairwiseCUDA2023
  //TODOC4PairwiseCUDA2023 the way of passing cn7/cn8 are different now, this line needs redone
  bool usingShare = MaxLJIndex > (cSim.TIVdwNTyes * (cSim.TIVdwNTyes + 1) / 2);

  if (usingShare) {
    cn = &cn0[0];
    sigEps = &sigEps0[0];
    dcoef = &dcoef0[0]; // C4PairwiseCUDA
    dvalue = &dvalue0[0]; //C4PairwiseCUDA2023
  } else {
    cn = &cSim.pTIcn[0];
    sigEps = &cSim.pTISigEps[0];
    dcoef = &cSim.pTIDcoef[0];
    dvalue = &cSim.pTIDvalue[0]; //C4PairwiseCUDA2023
  }

  __shared__ int my_iTIAtom;
  __shared__ PMEFloat3 myCoord;
  unsigned long long int myLocalForce[3] = { 0, 0, 0};
  unsigned long long int myLocalSCForce[3] = { 0, 0, 0 };

  __shared__ uint myNumberNB;

  if (cSim.pNTPData == NULL) {
    if (i < 3) {
      recip[i * 3] = cSim.recipf[i][0];
      recip[i * 3 + 1] = cSim.recipf[i][1];
      recip[i * 3 + 2] = cSim.recipf[i][2];
      ucell[i * 3] = cSim.ucellf[i][0];
      ucell[i * 3 + 1] = cSim.ucellf[i][1];
      ucell[i * 3 + 2] = cSim.ucellf[i][2];
    }
  } else {
    if (i < 9) {
      ucell[i] = cSim.pNTPData->ucellf[i];
      recip[i] = cSim.pNTPData->recipf[i];
    }
  }

  if (i == 0) {
    my_iTIAtom = cSim.pTINBList[blockIdx.x* gti_simulationConst::MaxNumberNBPerAtom].x;
    myNumberNB = cSim.pTINBList[(blockIdx.x+1) * gti_simulationConst::MaxNumberNBPerAtom -1].x;
    if (my_iTIAtom >= 0) {
      myCoord = { (PMEFloat)cSim.pImageX[my_iTIAtom] , (PMEFloat)cSim.pImageY[my_iTIAtom] , (PMEFloat)cSim.pImageZ[my_iTIAtom] };
    }
  }

  if (usingShare) {
    while (i < cSim.TIVdwNTyes * (cSim.TIVdwNTyes + 1) / 2) {
      cn[i] = cSim.pTIcn[i];
      sigEps[i] = cSim.pTISigEps[i];
      i += blockDim.x;
    }
    while (j < cSim.TIC4Pairwise) { //C4PairwiseCUDA2023
      dcoef[j*2+0] = cSim.pTIDcoef[j*2+0];
      dcoef[j*2+1] = cSim.pTIDcoef[j*2+1];
      dvalue[j] = cSim.pTIDvalue[j];
      j += blockDim.x;
    }
  }

  __syncthreads();
  
  if (myNumberNB == 0 || my_iTIAtom<0 ) return;

  static const uint shift = gti_simulationConst::MaxNumberNBPerAtom;
  uint bShift = blockIdx.x * shift;

  PMEFloat q1 = cSim.pOrigAtomCharge[cSim.pImageAtom[cSim.pTINBList[bShift].x]];

  int TIRegion = cSim.pTINBList[bShift].z & gti_simulationConst::Fg_TI_region;
  bool inReadMode = (cSim.reafMode == 1 - TIRegion);

  if (cSim.reafMode >= 0 && TIRegion == cSim.reafMode) return;

  int TISign = ((TIRegion == 0) ? -1 : 1);

  int& EleSCType = cSim.eleSC;
  PMEFloat lambda_e = (EleSCType == 1) ? smooth_step_func(cSim.TILambda[TIRegion], cSim.eleSmoothLambdaType) : cSim.TILambda[TIRegion];
  PMEFloat d_lambda_e = (EleSCType == 1) ? d_smooth_step_func(cSim.TILambda[TIRegion], cSim.eleSmoothLambdaType) : OneF;

  PMEFloat lambda_sc = (EleSCType == 1) ? smooth_step_func(cSim.TILambda[TIRegion], cSim.SCSmoothLambdaType) : cSim.TILambda[TIRegion];
  PMEFloat d_lambda_sc = (EleSCType == 1) ? d_smooth_step_func(cSim.TILambda[TIRegion], cSim.SCSmoothLambdaType) : OneF;

  int vdwSCType = cSim.vdwSC;
  PMEFloat lambda_v = (vdwSCType == 1)
    ? smooth_step_func(cSim.TILambda[TIRegion], cSim.vdwSmoothLambdaType)
    : cSim.TILambda[TIRegion];
  PMEFloat d_lambda_v = (vdwSCType == 1)
    ? d_smooth_step_func(cSim.TILambda[TIRegion], cSim.vdwSmoothLambdaType)
    : OneF;

  unsigned int pos = threadIdx.x;
  while (pos < myNumberNB) {

    int4 tt = cSim.pTINBList[pos+ bShift];
    if (tt.x >=0 ) {
      unsigned iatom0 = tt.x;
      unsigned iatom1 = tt.y;

      PMEFloat3 ddr = { (PMEFloat)cSim.pImageX[iatom1], (PMEFloat)cSim.pImageY[iatom1] ,(PMEFloat)cSim.pImageZ[iatom1] };
#ifdef AMBER_PLATFORM_AMD
      ddr = ddr - myCoord;
#else
      ddr -= myCoord;
#endif

      PMEFloat r2 = ddr.x * ddr.x + ddr.y * ddr.y + ddr.z * ddr.z;
      if (r2 > cSim.cut2)
        r2 = __image_dist2<PMEFloat>(ddr.x, ddr.y, ddr.z, recip, ucell);
      bool excluded = tt.z & gti_simulationConst::Fg_excluded;
      bool intTI = (tt.z & gti_simulationConst::Fg_int_TI) && cSim.tiCut > 0;

      if (r2 < cSim.cut2 || excluded || intTI) {

        bool addToDVDL_ele = tt.z & gti_simulationConst::ex_addToDVDL_ele;
        bool needSCEE = tt.z & gti_simulationConst::ex_SC_ELE;
        bool hasSC = tt.z & gti_simulationConst::Fg_has_SC;

        //aliases
        Schedule::InteractionType lambda_type = (hasSC) ? Schedule::TypeEleCC : Schedule::TypeEleSC;
        //PMEFloat& tiWeight_ele_rec = cSim.TIItemWeight[Schedule::TypeEleRec][TIRegion];
        PMEFloat& tiWeight_ele = cSim.TIItemWeight[lambda_type][TIRegion];
        PMEFloat& tiWeight_vdw = cSim.TIItemWeight[Schedule::TypeVDW][TIRegion];

        int& vdwIndex = tt.w;

        PMEFloat rinv = Rsqrt(r2);
        PMEFloat r = r2 * rinv;
        PMEFloat r2inv = rinv * rinv;

        unsigned long long int* pE, * pEDL;
        PMEFloat eEle = ZeroF;
        PMEFloat df = ZeroF;
        PMEFloat dfCorr = ZeroF;

        //Tail smoothing for SC
        PMEFloat r_ratio0 = cSim.cut_sc0 / cSim.cut_sc1;
        PMEFloat r_ratio = OneF - r_ratio0;
        PMEFloat rtp = r / cSim.cut_sc1;

        PMEFloat rtv, d_rtv; // vdw Tail smoothing
        if (cSim.cut_sc <= 0 || rtp < r_ratio0) {
          rtv = OneF;
          d_rtv = ZeroF;
        } else if (rtp > OneF) {
          rtv = ZeroF;
          d_rtv = ZeroF;
        } else {
          rtv = One - smooth_step_func((rtp - r_ratio0) / r_ratio, 2);
          d_rtv = -d_smooth_step_func((rtp - r_ratio0) / r_ratio, 2) / r_ratio;
        }

        PMEFloat rte = (cSim.cut_sc <= 1) ? OneF : rtv;
        PMEFloat d_rte = (cSim.cut_sc <= 1) ? ZeroF : d_rtv;

        //Tail smoothing for regular potentials
        if (cSim.fswitch > 0) {
          r_ratio0 = cSim.fswitch / cSim.cut;
          r_ratio = OneF - r_ratio0;
          rtp = r / cSim.cut;
        }

        PMEFloat rtr = (cSim.fswitch < 0 || rtp < r_ratio0) ? OneF :
          ((rtp > OneF) ? ZeroF : One - smooth_step_func((rtp - r_ratio0) / r_ratio, 2));
        PMEFloat d_rtr = (cSim.fswitch < 0 || rtp < r_ratio0) ? ZeroF :
          ((rtp > OneF) ? ZeroF : -d_smooth_step_func((rtp - r_ratio0) / r_ratio, 2) / r_ratio);

        uint EleScreenType = cSim.eleGauss;

        PMEFloat q1q2 = q1 * cSim.pOrigAtomCharge[cSim.pImageAtom[tt.y]];
        PMEFloat scalpha = (cSim.autoAlpha > 0) ? test_alpha(q1q2) : cSim.TISCAlpha;
        PMEFloat scbeta = (q1q2 < 0) ? cSim.TISCBeta : cSim.TISCGamma;

        PMEFloat erf_c = (excluded) ? -Erf(cSim.ew_coeffSP * r) : Erfc(cSim.ew_coeffSP * r);
        PMEFloat d_erf_c = -cSim.negTwoEw_coeffRsqrtPI * Exp(-cSim.ew_coeff2 * r2);
        if (inReadMode) {
          erf_c = (excluded) ? ZeroF : OneF;
          d_erf_c = ZeroF;
        }
        PMEFloat q1q2_erfc = q1q2 * erf_c;
        PMEFloat rsc_inv = ZeroF;
        PMEFloat rsc_3inv = ZeroF;
        PMEFloat d_rsc_inv_over_r = ZeroF;

        int& m = cSim.eleExp;

        if (cSim.scaleBeta > 0 && needSCEE) {
          PMEFloat sm = (vdwIndex < 0) ? (2 << (m - 1)) : cSim.pTISigMN[vdwIndex].x;
          scbeta *= max((PMEFloat)(2 << (m - 1)), sm);
        }

        if (!needSCEE) {
          eEle = q1q2_erfc * rinv;
          df = Fma(d_erf_c, q1q2, eEle) * r2inv;
          if (!addToDVDL_ele)
            dfCorr = q1q2 * rinv * r2inv * (OneF - tiWeight_ele);

        } else {
          switch (m) {
          case(1): rsc_inv = OneF / (r + scbeta * lambda_sc * rte); break;
          case(2): rsc_inv = Rsqrt(r2 + scbeta * lambda_sc * rte); break;
          case(4): rsc_inv = Rsqrt(Sqrt(r2 * r2 + scbeta * lambda_sc * rte)); break;
          case(6):default: rsc_inv = Rsqrt(Cbrt(r2 * r2 * r2 + scbeta * lambda_sc * rte)); break;
          }

          rsc_3inv = rsc_inv * rsc_inv * rsc_inv;

          d_rsc_inv_over_r = rsc_3inv;
          switch (m) {
          case(1): d_rsc_inv_over_r = rsc_inv * rsc_inv * (OneF + scbeta * lambda_sc * d_rte) * rinv; break;
          case(2): d_rsc_inv_over_r = rsc_3inv * (OneF + scbeta * lambda_sc * d_rte * rinv * HalfF); break;
          case(4): d_rsc_inv_over_r = QuarterF * rsc_3inv * rsc_inv * rsc_inv * (4 * r * r2 + scbeta * lambda_sc * d_rte) * rinv; break;
          case(6):default: d_rsc_inv_over_r = OneF / SixF * rsc_3inv * rsc_3inv * rsc_inv * (SixF * r2 * r2 * r + scbeta * lambda_sc * d_rte) * rinv; break;
          }

          switch (EleScreenType) {
            //case(0): df = (q1q2 * d_erf_c * rsc_inv * rinv) + (q1q2_erfc * d_rsc_inv); break;
          default:case(0): df = Fma(q1q2 * d_erf_c, rsc_inv * rinv, q1q2_erfc * d_rsc_inv_over_r); break;
            //case(1): df = (q1q2 * d_erf_c + q1q2_erfc * rinv) * r2inv + q1q2 *(d_rsc_inv - rinv * r2inv); break;
          case(1): df = (Fma(q1q2, d_erf_c, q1q2_erfc * rinv)) * r2inv + q1q2 * Fma(-rinv, r2inv, d_rsc_inv_over_r); break;
            break;
          }
        }

        if (needEnergy) {

          pE = (hasSC) ? cSim.pTIEESC[TIRegion] : cSim.pTIEECC[TIRegion];
          pEDL = (hasSC) ? cSim.pTISCEESC_DL[TIRegion] : cSim.pTISCEECC_DL[TIRegion];

          if (!needSCEE) {

            addEnergy(pE, eEle);
            if (!addToDVDL_ele) {
              addEnergy(cSim.pTISCEED[TIRegion], q1q2 * rinv);
              addEnergy(pE, -q1q2 * rinv);
            }

          } else {

            //lambda-"independent" contribution
            switch (EleScreenType) {
            default:case(0): eEle = q1q2_erfc * rsc_inv; break;
            case(1): eEle = q1q2_erfc * rinv + q1q2 * (rsc_inv - rinv); break;
            }
            addEnergy(pE, eEle);

            //lambda-dependent contribution
            PMEFloat qq = 0.0;
            switch (EleScreenType) {
            default:case(0): qq = q1q2_erfc; break;
            case(1): qq = q1q2; break;
            }
            PMEFloat dl = qq * TISign * scbeta * d_lambda_sc * rte;
            switch (m) {
            case(1): dl *= rsc_inv * rsc_inv; break;
            case(2):default: dl *= HalfF * rsc_3inv; break;
            case(4): dl *= QuarterF * rsc_3inv * rsc_inv * rsc_inv; break;
            case(6): dl *= OneF / SixF * rsc_3inv * rsc_3inv * rsc_inv; break;
            }
            addEnergy(pEDL, dl);

            // MBAR Stuff
            if (cSim.needMBAR) {
              unsigned i = 0;
              unsigned j = cSim.nMBARStates-1;
              if (cSim.currentMBARState >= 0 && cSim.rangeMBAR >= 0) {
                i = max(0, cSim.currentMBARState - cSim.rangeMBAR);
                j = min(cSim.nMBARStates-1, cSim.currentMBARState + cSim.rangeMBAR);
              }
              PMEFloat* lambda[2] = { cSim.pMBARLambda, &(cSim.pMBARLambda[cSim.nMBARStates]) };
              unsigned mbarShift = cSim.nMBARStates * (TIRegion + lambda_type * 3);
              unsigned long long int* MBAREnergy = &(cSim.pMBAREnergy[mbarShift]);
              for (unsigned l = i; l <= j; l++) {
                PMEFloat rsc_inv_l;
                switch (m) {
                case(1): rsc_inv_l = (EleSCType == 1)
                  ? OneF / (r + scbeta * rte * (smooth_step_func(lambda[TIRegion][l], cSim.SCSmoothLambdaType)))
                  : OneF / (r + scbeta * rte * lambda[TIRegion][l]); break;
                case(2):default: rsc_inv_l = (EleSCType == 1)
                  ? Rsqrt(r2 + scbeta * rte * (smooth_step_func(lambda[TIRegion][l], cSim.SCSmoothLambdaType)))
                  : Rsqrt(r2 + scbeta * rte * lambda[TIRegion][l]); break;
                case(4): rsc_inv_l = (EleSCType == 1)
                  ? Rsqrt(Sqrt(r2 * r2 + scbeta * rte * (smooth_step_func(lambda[TIRegion][l], cSim.SCSmoothLambdaType))))
                  : Rsqrt(Sqrt(r2 * r2 + scbeta * rte * lambda[TIRegion][l])); break;
                case(6): rsc_inv_l = (EleSCType == 1)
                  ? Rsqrt(Cbrt(r2 * r2 * r2 + scbeta * rte * (smooth_step_func(lambda[TIRegion][l], cSim.SCSmoothLambdaType))))
                  : Rsqrt(Cbrt(r2 * r2 * r2 + scbeta * rte * lambda[TIRegion][l])); break;
                }

                addEnergy(&(MBAREnergy[l]), qq * (rsc_inv_l - rsc_inv));
              }
            }  // if (cSim.needMBAR) {
          }
        }

        PMEFloat3 TIForce = { 0.0, 0.0, 0.0 }, TISCForce = { 0.0, 0.0, 0.0 };
        PMEFloat3 V = { 0.0, 0.0, 0.0 }, VSC = { 0.0, 0.0, 0.0 };

        PMEFloat3 fdd = df * ddr;

        //if (addToDVDL_ele) {
        TIForce = fdd * tiWeight_ele;
        addLocalForce(&myLocalForce[0], -TIForce.x);
        addLocalForce(&myLocalForce[1], -TIForce.y);
        addLocalForce(&myLocalForce[2], -TIForce.z);

        addForce(cSim.pTINBForceX[TIRegion] + iatom1, TIForce.x);
        addForce(cSim.pTINBForceY[TIRegion] + iatom1, TIForce.y);
        addForce(cSim.pTINBForceZ[TIRegion] + iatom1, TIForce.z);
        if (needVirial) V = -TIForce * ddr;

        if (!addToDVDL_ele) {
          TISCForce = dfCorr * ddr;
          if (needVirial) VSC = -TISCForce * ddr;
        }

        // EVDW
        int& n = cSim.vdwExp;
        bool addToDVDL_vdw = tt.z & gti_simulationConst::ex_addToDVDL_vdw;
        bool doLJ = (vdwIndex >= 0 && (tt.z & gti_simulationConst::ex_need_LJ) && !excluded);
        bool needSCLJ = tt.z & gti_simulationConst::ex_SC_LJ;

        df = ZeroF;
        dfCorr = Zero;
        TIForce = { 0.0, 0.0, 0.0 };
        fdd = { 0.0, 0.0, 0.0 };
        if (doLJ) {

          PMEFloat t0, t4, f4, t6, f6, f12, eLJ, u1_3, dvdu;
          PMEFloat& sn = cSim.pTISigMN[vdwIndex].y;

          PMEFloat rsc, rnsc;
          t0 = scalpha * lambda_v * rtv;
          switch (n) {
          case(1):
            rsc = (r + t0 * sn);
            rnsc = rsc;
            break;
          case(2): default:
            rnsc = r2 + t0 * sn;
            rsc = Sqrt(rnsc);
            break;
          case(4):
            rnsc = r2 * r2 + t0 * sn;
            rsc = Sqrt(Sqrt(rnsc));
            break;
          case(6):
            rnsc = r2 * r2 * r2 + t0 * sn;
            rsc = Sqrt(Cbrt(rnsc));
            break;
          }
          PMEFloat rsc_2 = rsc * rsc;
          if (needSCLJ && vdwSCType < 2) {

            f6 = OneF / (rsc_2 * rsc_2 * rsc_2 * sigEps[vdwIndex].x);
            f12 = f6 * f6;
            eLJ = sigEps[vdwIndex].y * (f12 - f6);
            // v is E_LJ and U is f6
            dvdu = sigEps[vdwIndex].y * (OneF - TwoF * f6);

            // C4 term for 12-6-4 potential 
            if (sigEps[vdwIndex].w > 0) {
              t4 = sigEps[vdwIndex].z;
              u1_3 = Cbrt(f6);
              f4 = t4 * u1_3 * u1_3;
              eLJ -= f4;
              dvdu += (TwoF * t4) / (ThreeF * u1_3);
            }
            if (cSim.TIC4Pairwise) {
              for (int i=0;i<cSim.TIC4Pairwise;i++) {
                if ((cSim.pImageAtom[tt.x] == (cSim.pTIDcoef[i*2]) && cSim.pImageAtom[tt.y] == (cSim.pTIDcoef[i*2+1])) || (cSim.pImageAtom[tt.x] == (cSim.pTIDcoef[i*2+1]) && cSim.pImageAtom[tt.y] == (cSim.pTIDcoef[i*2]))) {
                  f4 += cSim.pTIDvalue[i] * t4 / 3.0; // C4PairwiseCUDA2023
                  eLJ -= f4;
                  dvdu += (TwoF * t4) / (ThreeF * u1_3);
                }
              }
            }
            //if (cSim.TIC4Pairwise > 0) {  //C4PairwiseCUDA2023
            //  t4 = dvalue[0];
            //  u1_3 = Cbrt(f6);
            //  f4 = t4 * u1_3 * u1_3;
            //  eLJ -= f4;
            //  dvdu += (TwoF * t4) / (ThreeF * u1_3);
            //}
            PMEFloat dudr = OneF;
            switch (n) {
            case(1): break;
            case(2): dudr *= r; break;
            case(4): dudr *= r * r2; break;
            case(6):default: dudr *= r * r2 * r2; break;
            }

            dudr = -SixF * f6 / rnsc * rinv * (dudr + lambda_v * scalpha * sn * d_rtv / n);
            df = dvdu * dudr;
          } else {
            t4 = r2inv * r2inv;
            f4 = (sigEps[vdwIndex].w > 0) ? cn[vdwIndex].z * t4 : ZeroF;
            //f4 = (cSim.TIC4Pairwise > 0) ? dvalue[0] * t4 : ZeroF;  //C4PairwiseCUDA2023
            t6 = t4 * r2inv;
            f6 = cn[vdwIndex].y * t6;
            f12 = cn[vdwIndex].x * t6 * t6;
            if (vdwIndex == 55555) {
              for (int i=0;i<cSim.TIC4Pairwise;i++) {
                if ((cSim.pImageAtom[tt.x] == (cSim.pTIDcoef[i*2]) && cSim.pImageAtom[tt.y] == (cSim.pTIDcoef[i*2+1])) || (cSim.pImageAtom[tt.x] == (cSim.pTIDcoef[i*2+1]) && cSim.pImageAtom[tt.y] == (cSim.pTIDcoef[i*2]))) {
                  f4 += cSim.pTIDvalue[i] * t4 / 3.0; // C4PairwiseCUDA2023
                  eLJ -= f4;
                  dvdu += (TwoF * t4) / (ThreeF * u1_3);
                }
              }
            }
            df = (f12 * TwelveF - f6 * SixF - f4 * FourF) * r2inv;
            eLJ = (f12 - f6 - f4);
            if (needSCLJ) {  //Capped "SC"
              PMEFloat eCut = cSim.vdwCap;
              PMEFloat eMax = eCut + 10.0;;
              PMEFloat r2cut = Pow((cn[vdwIndex].x / eCut), OneF / SixF);
              PMEFloat qq = r2 / r2cut;
              if (qq < OneF) {
                eLJ = eMax - (eMax - eCut) * qq;
                df = TwoF * (eMax - eCut) * qq * r2inv;
              }
            }
          }
          df = df * rtr - eLJ * d_rtr / r;
          eLJ *= rtr;

          if (needEnergy) {
            pE = (addToDVDL_vdw) ? cSim.pTILJ[TIRegion] : cSim.pTISCLJ[TIRegion];
            addEnergy(pE, eLJ);

            if (needSCLJ) {
              if (vdwSCType < 2) {
                PMEFloat dudl = -scalpha * f6 * TISign * d_lambda_v * rtv;
                dudl *= SixF / (OneF * n) * sn / rnsc;
                addEnergy(cSim.pTISCLJ_DL[TIRegion], dvdu * dudl * rtr);

                // MBAR Stuff
                if (cSim.needMBAR && needEnergy) {
                  unsigned mbarShift = cSim.nMBARStates * (TIRegion + Schedule::TypeVDW * 3);
                  unsigned long long int* MBAREnergy = &(cSim.pMBAREnergy[mbarShift]);
                  PMEFloat* lambda[2] = { cSim.pMBARLambda , &(cSim.pMBARLambda[cSim.nMBARStates]) };

                  t0 = scalpha * rtv;
                  unsigned i = 0;
                  unsigned j = cSim.nMBARStates-1;
                  if (cSim.currentMBARState >= 0 && cSim.rangeMBAR >= 0) {
                    i = max(0, cSim.currentMBARState - cSim.rangeMBAR);
                    j = min(cSim.nMBARStates - 1, cSim.currentMBARState + cSim.rangeMBAR);
                  }
                  for (unsigned l = i; l <=j ; l++) {
                    t6 = (vdwSCType >= 1) ?
                      t0 * (smooth_step_func(lambda[TIRegion][l], cSim.vdwSmoothLambdaType))
                      : t0 * lambda[TIRegion][l];

                    f6 = OneF / sigEps[vdwIndex].x;
                    switch (n) {
                    case(1):
                      rsc = r + t6 * sn;
                      rsc_2 = rsc * rsc;
                      f6 /= (rsc_2 * rsc_2 * rsc_2);
                      break;
                    case(2):
                      rsc_2 = r2 + t6 * sn;
                      f6 /= (rsc_2 * rsc_2 * rsc_2);
                      break;
                    case(4):
                      rsc_2 = r2 * r2 + t6 * sn;
                      f6 /= (rsc_2 * Sqrt(rsc_2));
                      break;
                    case(6):default:
                      rsc_2 = r2 * r2 * r2 + t6 * sn;
                      f6 /= rsc_2;
                      break;
                    }

                    f12 = f6 * f6;
                    PMEFloat dl_bar = sigEps[vdwIndex].y * (f12 - f6);
                    if (sigEps[vdwIndex].w > 0) {
                      u1_3 = cbrt(f6);
                      dl_bar -= t4 * u1_3 * u1_3;
                    }
                    addEnergy(&(MBAREnergy[l]), dl_bar * rtr - eLJ);
                  }
                }
              }
            }

          }

          fdd = df * ddr;
        }

        if (addToDVDL_vdw) {
          TIForce = fdd * tiWeight_vdw;
          addLocalForce(&myLocalForce[0], -TIForce.x);
          addLocalForce(&myLocalForce[1], -TIForce.y);
          addLocalForce(&myLocalForce[2], -TIForce.z);

          addForce(cSim.pTINBForceX[TIRegion] + iatom1, TIForce.x);
          addForce(cSim.pTINBForceY[TIRegion] + iatom1, TIForce.y);
          addForce(cSim.pTINBForceZ[TIRegion] + iatom1, TIForce.z);
        } else {
#ifdef AMBER_PLATFORM_AMD
          TISCForce = TISCForce + fdd;
          if (needVirial) VSC = VSC - (ddr * TISCForce);
#else
          TISCForce += fdd;
          if (needVirial) VSC -= ddr * TISCForce;
#endif
        }

        if (needVirial) {
#ifdef AMBER_PLATFORM_AMD
          V = V - (ddr *  TIForce);
#else
          V -= ddr *  TIForce;
#endif
          addVirial(cSim.pV11[TIRegion], V.x);
          addVirial(cSim.pV22[TIRegion], V.y);
          addVirial(cSim.pV33[TIRegion], V.z);
        }

        if (!addToDVDL_vdw || !addToDVDL_ele) {
          addLocalForce(&myLocalSCForce[0], -TISCForce.x);
          addLocalForce(&myLocalSCForce[1], -TISCForce.y);
          addLocalForce(&myLocalSCForce[2], -TISCForce.z);
          addForce(cSim.pTINBSCForceX[TIRegion] + iatom1, TISCForce.x);
          addForce(cSim.pTINBSCForceY[TIRegion] + iatom1, TISCForce.y);
          addForce(cSim.pTINBSCForceZ[TIRegion] + iatom1, TISCForce.z);
          if (needVirial) {
            addVirial(cSim.pVirial_11, VSC.x);
            addVirial(cSim.pVirial_22, VSC.y);
            addVirial(cSim.pVirial_33, VSC.z);
          }
        }

      }

    }
    else break;

    pos += blockDim.x;
  }

  __syncthreads();

  for (int offset = warpSize / 2; offset > 0; offset /= 2) {
#ifndef AMBER_PLATFORM_AMD
    __syncwarp();
#endif
    myLocalForce[0] += __SHFL_DOWN(0xFFFFFFFF, myLocalForce[0], offset);
    myLocalForce[1] += __SHFL_DOWN(0xFFFFFFFF, myLocalForce[1], offset);
    myLocalForce[2] += __SHFL_DOWN(0xFFFFFFFF, myLocalForce[2], offset);

    myLocalSCForce[0] += __SHFL_DOWN(0xFFFFFFFF, myLocalSCForce[0], offset);
    myLocalSCForce[1] += __SHFL_DOWN(0xFFFFFFFF, myLocalSCForce[1], offset);
    myLocalSCForce[2] += __SHFL_DOWN(0xFFFFFFFF, myLocalSCForce[2], offset);

  }

  __syncthreads();

  if ((threadIdx.x & (warpSize - 1)) == 0) {
    atomicAdd((unsigned long long int*)cSim.pTINBForceX[TIRegion] + my_iTIAtom, myLocalForce[0]);
    atomicAdd((unsigned long long int*)cSim.pTINBForceY[TIRegion] + my_iTIAtom, myLocalForce[1]);
    atomicAdd((unsigned long long int*)cSim.pTINBForceZ[TIRegion] + my_iTIAtom, myLocalForce[2]);
    atomicAdd((unsigned long long int*)cSim.pTINBSCForceX[TIRegion] + my_iTIAtom, myLocalSCForce[0]);
    atomicAdd((unsigned long long int*)cSim.pTINBSCForceY[TIRegion] + my_iTIAtom, myLocalSCForce[1]);
    atomicAdd((unsigned long long int*)cSim.pTINBSCForceZ[TIRegion] + my_iTIAtom, myLocalSCForce[2]);
  }

}

//---------------------------------------------------------------------------------------------
// kgCalculateREAFNb_kernel:
//
// Arguments:
//   needEnergy:
//   needVirial:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgCalculateREAFNb_kernel(bool needEnergy, bool needVirial) {

  unsigned i = threadIdx.x;
  __shared__ PMEFloat recip[9], ucell[9];

  __shared__ PMEFloat4 cn[MaxLJIndex];
  __shared__ PMEFloat4 sigEps[MaxLJIndex];

  __shared__ int my_iREAFAtom;
  __shared__ PMEFloat3 myCoord;
  unsigned long long int myLocalForce[3] = { 0, 0, 0 };

  __shared__ uint myNumberNB;

  if (cSim.pNTPData == NULL) {
    if (i < 3) {
      recip[i * 3] = cSim.recipf[i][0];
      recip[i * 3 + 1] = cSim.recipf[i][1];
      recip[i * 3 + 2] = cSim.recipf[i][2];
      ucell[i * 3] = cSim.ucellf[i][0];
      ucell[i * 3 + 1] = cSim.ucellf[i][1];
      ucell[i * 3 + 2] = cSim.ucellf[i][2];
    }
  }
  else {
    if (i < 9) {
      ucell[i] = cSim.pNTPData->ucellf[i];
      recip[i] = cSim.pNTPData->recipf[i];
    }
  }

  if (i == 0) {
    my_iREAFAtom = cSim.pREAFNbList[blockIdx.x * gti_simulationConst::MaxNumberNBPerAtom].x;
    if (my_iREAFAtom >= 0) {
      myNumberNB = cSim.pREAFNbList[(blockIdx.x + 1) * gti_simulationConst::MaxNumberNBPerAtom - 1].x;
      myCoord = { (PMEFloat)cSim.pImageX[my_iREAFAtom] , (PMEFloat)cSim.pImageY[my_iREAFAtom] , (PMEFloat)cSim.pImageZ[my_iREAFAtom] };
    }
  }

  while (i < cSim.TIVdwNTyes * (cSim.TIVdwNTyes + 1) / 2) {
    cn[i] = cSim.pTIcn[i];
    sigEps[i] = cSim.pTISigEps[i];
    i += blockDim.x;
  }

  __syncthreads();

  if (myNumberNB == 0 || my_iREAFAtom < 0) return;

  static const uint shift = gti_simulationConst::MaxNumberNBPerAtom;
  uint bShift = blockIdx.x * shift;

  PMEFloat q1 = cSim.pOrigAtomCharge[cSim.pImageAtom[cSim.pREAFNbList[bShift].x]];

  unsigned int pos = threadIdx.x;
  while (pos < myNumberNB) {

    int4 tt = cSim.pREAFNbList[pos + bShift];
    if (tt.x >= 0) {

      unsigned iatom1 = tt.y;

      PMEFloat3 ddr = { (PMEFloat)cSim.pImageX[iatom1], (PMEFloat)cSim.pImageY[iatom1] , (PMEFloat)cSim.pImageZ[iatom1] };
#ifdef AMBER_PLATFORM_AMD
      ddr = ddr - myCoord;
#else
      ddr -= myCoord;
#endif

      PMEFloat r2 = ddr.x * ddr.x + ddr.y * ddr.y + ddr.z * ddr.z;
      if (r2 > cSim.cut2)
        r2 = __image_dist2<PMEFloat>(ddr.x, ddr.y, ddr.z, recip, ucell);

      bool excluded = tt.z & gti_simulationConst::Fg_excluded;
      bool intRE = (tt.z & gti_simulationConst::Fg_int_RE) && cSim.tiCut > 0;

      bool addRE_ele = tt.z & gti_simulationConst::ex_addRE_ele;
      bool addRE_vdw = tt.z & gti_simulationConst::ex_addRE_vdw;

      // Here the elec part is only a correction to the already scaled elec interaction
      Schedule::InteractionType tauType = (intRE) ? Schedule::TypeEleCC : Schedule::TypeEleSC;
      PMEFloat weight_ele = (intRE) ? -cSim.REAFItemWeight[Schedule::TypeEleRec][1]
        : -cSim.REAFItemWeight[Schedule::TypeEleRec][0];
      weight_ele += ((intRE) ? cSim.REAFItemWeight[tauType][1] : cSim.REAFItemWeight[tauType][0]);

      // Here the vdw part is a correction as well
      PMEFloat weight_vdw = (addRE_vdw) ? ( (intRE) ? cSim.REAFItemWeight[Schedule::TypeVDW][1] : cSim.REAFItemWeight[Schedule::TypeVDW][0]) - OneF
        : ZeroF;

      if (r2 < cSim.cut2 || excluded || intRE) {

        PMEFloat rinv = Rsqrt(r2);
        PMEFloat r = r2 * rinv;
        PMEFloat r2inv = rinv * rinv;

        PMEFloat df = ZeroF;
        PMEFloat3 READForce = { 0.0, 0.0, 0.0 };
        PMEFloat3 fdd = { 0.0, 0.0, 0.0 };
        PMEFloat3 V = { 0.0, 0.0, 0.0 };

        bool doEle = (addRE_ele && abs(weight_ele) > 1e-6 && !excluded);
        if (doEle){

          PMEFloat q1q2 = q1 * cSim.pOrigAtomCharge[cSim.pImageAtom[tt.y]];
          if (needEnergy) addEnergy(cSim.pEED, q1q2 * rinv * weight_ele);
          fdd = q1q2 * r2inv * rinv* ddr;

          READForce = fdd * weight_ele;
          addLocalForce(&myLocalForce[0], -READForce.x);
          addLocalForce(&myLocalForce[1], -READForce.y);
          addLocalForce(&myLocalForce[2], -READForce.z);

          addForce(cSim.pNBForceXAccumulator + iatom1, READForce.x);
          addForce(cSim.pNBForceYAccumulator + iatom1, READForce.y);
          addForce(cSim.pNBForceZAccumulator + iatom1, READForce.z);

          if (needVirial) V = -READForce * ddr;
        }

        // EVDW
        int& vdwIndex = tt.w;
        bool doLJ = (abs(weight_vdw) >1e06 && vdwIndex >= 0 && (tt.z & gti_simulationConst::ex_need_LJ) && !excluded);
        if (doLJ) {
          df = ZeroF;
          READForce = { 0.0, 0.0, 0.0 };
          fdd = { 0.0, 0.0, 0.0 };

          PMEFloat t4, f4, t6, f6, f12;

          t4 = r2inv * r2inv;
          f4 = (sigEps[vdwIndex].w > 0) ? cn[vdwIndex].z * t4 : ZeroF;
          t6 = t4 * r2inv;
          f6 = cn[vdwIndex].y * t6;
          f12 = cn[vdwIndex].x * t6 * t6;
          df = (f12 * TwelveF - f6 * SixF - f4 * FourF) * r2inv;
          if (vdwIndex == 55555) {
              for (int i=0;i<cSim.TIC4Pairwise;i++) {
                if ((cSim.pImageAtom[tt.x] == (cSim.pTIDcoef[i*2]) && cSim.pImageAtom[tt.y] == (cSim.pTIDcoef[i*2+1])) || (cSim.pImageAtom[tt.x] == (cSim.pTIDcoef[i*2+1]) && cSim.pImageAtom[tt.y] == (cSim.pTIDcoef[i*2]))) {
                  f4 += cSim.pTIDvalue[i] * t4 / 3.0; // C4PairwiseCUDA2023
                  //eLJ -= f4;
                  //dvdu += (TwoF * t4) / (ThreeF * u1_3);
                }
              }
            }
          if (needEnergy)  addEnergy(cSim.pEVDW, (f12 - f6 - f4) * weight_vdw);
          fdd = df * ddr;

          READForce = fdd * weight_vdw;
          addLocalForce(&myLocalForce[0], -READForce.x);
          addLocalForce(&myLocalForce[1], -READForce.y);
          addLocalForce(&myLocalForce[2], -READForce.z);

          addForce(cSim.pNBForceXAccumulator + iatom1, READForce.x);
          addForce(cSim.pNBForceYAccumulator + iatom1, READForce.y);
          addForce(cSim.pNBForceZAccumulator + iatom1, READForce.z);
        }

        if (needVirial) {
#ifdef AMBER_PLATFORM_AMD
          V = V - (ddr * READForce);
#else
          V -= ddr * READForce;
#endif
          addVirial(cSim.pVirial_11, V.x);
          addVirial(cSim.pVirial_22, V.y);
          addVirial(cSim.pVirial_33, V.z);
        }
      }

    }
    else break;

    pos += blockDim.x;
  }

  __syncthreads();

  for (int offset = warpSize / 2; offset > 0; offset /= 2) {
#ifndef AMBER_PLATFORM_AMD
    __syncwarp();
#endif
    myLocalForce[0] += __SHFL_DOWN(0xFFFFFFFF, myLocalForce[0], offset);
    myLocalForce[1] += __SHFL_DOWN(0xFFFFFFFF, myLocalForce[1], offset);
    myLocalForce[2] += __SHFL_DOWN(0xFFFFFFFF, myLocalForce[2], offset);

  }

  __syncthreads();

  if ((threadIdx.x & (warpSize - 1)) == 0) {
    atomicAdd((unsigned long long int*)cSim.pNBForceXAccumulator + my_iREAFAtom, myLocalForce[0]);
    atomicAdd((unsigned long long int*)cSim.pNBForceYAccumulator + my_iREAFAtom, myLocalForce[1]);
    atomicAdd((unsigned long long int*)cSim.pNBForceZAccumulator + my_iREAFAtom, myLocalForce[2]);
  }

}


//---------------------------------------------------------------------------------------------
// kgCalculateTI14NB_kernel:
//
// Arguments:
//   energy:
//   needVirial:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgCalculateTI14NB_kernel(bool needEnergy, bool needVirial) {
  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;

  unsigned long long int  *pE, *pEDL;
  PMEAccumulator* pFx, *pFy, *pFz;
  bool useSC14 = false;

  while (pos < cSim.numberTI14NBEntries) {
    unsigned atom0 = cSim.pTINb14ID[pos].x;
    unsigned atom1 = cSim.pTINb14ID[pos].y;
    unsigned TIRegion = cSim.pTINb14ID[pos].z;
    bool isIntSC = (cSim.pSCList[atom0] > 0 && cSim.pSCList[atom1] > 0);
    bool hasSC = (cSim.pSCList[atom0] > 0 || cSim.pSCList[atom1] > 0);

    int TISign = (TIRegion * 2) - 1;

    Schedule::InteractionType lambda_type = (hasSC) ? Schedule::TypeEleSC : Schedule::TypeEleCC;
    //PMEFloat& tiWeight_ele_rec = cSim.TIItemWeight[Schedule::TypeEleRec][TIRegion];
    PMEFloat& tiWeight_ele = cSim.TIItemWeight[lambda_type][TIRegion];
    PMEFloat& tiWeight_vdw = cSim.TIItemWeight[Schedule::TypeVDW][TIRegion];

    PMEFloat g1 = cSim.pTINb141[pos].x;


    PMEFloat f6 = ZeroF;
    PMEFloat f12 = ZeroF;
    unsigned iatom0 = cSim.pImageAtomLookup[atom0];
    unsigned iatom1 = cSim.pImageAtomLookup[atom1];
    PMEFloat dx = cSim.pImageX[iatom1] - cSim.pImageX[iatom0];
    PMEFloat dy = cSim.pImageY[iatom1] - cSim.pImageY[iatom0];
    PMEFloat dz = cSim.pImageZ[iatom1] - cSim.pImageZ[iatom0];
    PMEFloat r2 = dx*dx + dy*dy + dz*dz;
    PMEFloat rinv = Rsqrt(r2);
    PMEFloat r2inv = rinv*rinv;
    PMEFloat df = ZeroF;

    int& EleSCType = cSim.eleSC;
    PMEFloat lambda_e = (EleSCType == 1) ? smooth_step_func(cSim.TILambda[TIRegion], cSim.eleSmoothLambdaType) : cSim.TILambda[TIRegion];
    PMEFloat d_lambda_e = (EleSCType == 1) ? d_smooth_step_func(cSim.TILambda[TIRegion], cSim.eleSmoothLambdaType) : OneF;

    PMEFloat lambda_sc = (EleSCType == 1) ? smooth_step_func(cSim.TILambda[TIRegion], cSim.SCSmoothLambdaType) : cSim.TILambda[TIRegion];
    PMEFloat d_lambda_sc = (EleSCType == 1) ? d_smooth_step_func(cSim.TILambda[TIRegion], cSim.SCSmoothLambdaType) : OneF;

    int vdwSCType = cSim.vdwSC;
    PMEFloat lambda_v = (vdwSCType == 1)
      ? smooth_step_func(cSim.TILambda[TIRegion], cSim.vdwSmoothLambdaType)
      : cSim.TILambda[TIRegion];
    PMEFloat d_lambda_v = (vdwSCType == 1)
      ? d_smooth_step_func(cSim.TILambda[TIRegion], cSim.vdwSmoothLambdaType)
      : OneF;

    PMEFloat& scbeta = (g1 < 0) ? cSim.TISCBeta : cSim.TISCGamma;
    PMEFloat scalpha = (cSim.autoAlpha > 0) ? test_alpha(g1) : cSim.TISCAlpha;

    bool addToDVDL_ele = false, addToDVDL_vdw = false;

    switch (cSim.addSC) {
      case 0: addToDVDL_ele = (!hasSC); addToDVDL_vdw = (!hasSC); break;   // only non-SC part
      case 1: default:  addToDVDL_ele = (!isIntSC); addToDVDL_vdw = (!isIntSC); break; // add common-SC (ele+vdw) boundary
      case 2: addToDVDL_ele = true; addToDVDL_vdw = (!isIntSC); break;  // add SC internal terms (ele only)
      case 25: addToDVDL_ele = true; addToDVDL_vdw = true;  useSC14 = hasSC; break;
      case 3: case 4: case 5: case 6:addToDVDL_ele = true; addToDVDL_vdw = true; useSC14 = hasSC;  break;  //  add all NB SC terms
      case -1: addToDVDL_ele = (!isIntSC); addToDVDL_vdw = (!hasSC); break; // add common-SC (ele only) boundary
     }

    // EED: note that 1:4 electrostatics are double-counted
    PMEFloat de = 0.0;
    PMEFloat ddf = 0.0;
    PMEFloat den = 0.0, den_n=0.0;
    if (useSC14) {
      den = Rsqrt(r2 + (scbeta * lambda_sc));
      de = g1 * den;
      den_n = powf(den, 3.0);
      ddf = g1 * den_n;
    } else {
      de = g1 * rinv;
      ddf = de * r2inv;
    }

    if (needEnergy) {
      pE = (addToDVDL_ele) ? ((hasSC) ? cSim.pTIEE14SC[TIRegion] : cSim.pTIEE14CC[TIRegion]) : cSim.pTISCEE14[TIRegion];
      addEnergy(pE, de);

      if (useSC14 && addToDVDL_ele) {
        pEDL = (hasSC) ? cSim.pTIEE14SC_DL[TIRegion] : cSim.pTIEE14CC_DL[TIRegion];
        PMEFloat dl = 0.5 * g1 * den_n * TISign * scbeta * d_lambda_sc;
        addEnergy(pEDL, dl);
        // MBAR Stuff
        if (cSim.needMBAR) {
          PMEFloat* lambda[2] = { cSim.pMBARLambda, &(cSim.pMBARLambda[cSim.nMBARStates]) };
          unsigned mbarShift = cSim.nMBARStates * (TIRegion + lambda_type * 3);
          unsigned long long int* MBAREnergy = &(cSim.pMBAREnergy[mbarShift]);
          unsigned i = 0;
          unsigned j = cSim.nMBARStates - 1;
          if (cSim.currentMBARState >= 0 && cSim.rangeMBAR >= 0) {
            i = max(0, cSim.currentMBARState - cSim.rangeMBAR);
            j = min(cSim.nMBARStates - 1, cSim.currentMBARState + cSim.rangeMBAR);
          }
          for (unsigned l = i; l <= j; l++) {
            PMEFloat den_l = (EleSCType == 1)
              ? Rsqrt(r2 + scbeta * (smooth_step_func(lambda[TIRegion][l], cSim.SCSmoothLambdaType)))
              : Rsqrt(r2 + scbeta * lambda[TIRegion][l]);
            addEnergy(&(MBAREnergy[l]), g1 * (den_l - den));
          }
        }// if (cSim.needMBAR) {
      }

    }
    df += ddf;

    if (addToDVDL_ele) {
      pFx = cSim.pTINBForceX[TIRegion];
      pFy = cSim.pTINBForceY[TIRegion];
      pFz = cSim.pTINBForceZ[TIRegion];
    } else {
      pFx = cSim.pTINBSCForceX[TIRegion];
      pFy = cSim.pTINBSCForceY[TIRegion];
      pFz = cSim.pTINBSCForceZ[TIRegion];
    }

    PMEFloat ttx = df * dx;
    PMEFloat tty = df * dy;
    PMEFloat ttz = df * dz;
    if (addToDVDL_ele) {
      ttx *= tiWeight_ele;
      tty *= tiWeight_ele;
      ttz *= tiWeight_ele;
    }
    addForce(pFx + iatom0, -ttx);
    addForce(pFx + iatom1, ttx);
    addForce(pFy + iatom0, -tty);
    addForce(pFy + iatom1, tty);
    addForce(pFz + iatom0, -ttz);
    addForce(pFz + iatom1, ttz);
    if (needVirial && addToDVDL_ele) {
      addVirial(cSim.pV11[TIRegion], -ttx * dx);
      addVirial(cSim.pV22[TIRegion], -tty * dy);
      addVirial(cSim.pV33[TIRegion], -ttz * dz);
    }

    // VDW part
    PMEFloat g2 = cSim.pTINb141[pos].y;
    PMEFloat A = cSim.pTINb142[pos].x;
    PMEFloat B = cSim.pTINb142[pos].y;
    if (A < delta || B < delta) {
      pos += gridDim.x * blockDim.x;
      continue;
    }
    PMEFloat BoverA = cSim.pTINb142[pos].z;
    PMEFloat t6 = ZeroF;
    PMEFloat rt6 = ZeroF;

    df = OneF; // df must be 1.0 here
    // EVDW: note that 14 NB is not double-counted since it is zeroed-out
    t6 = r2inv * r2inv * r2inv;

    if (useSC14) {
      rt6 = r2 * r2 * r2 * BoverA;
      t6 = (BoverA) / ((scalpha * lambda_v) + rt6);
      df = t6 * r2 * r2 * r2;
    }
    f6 = g2 * B * t6;
    f12 = g2 * A * t6 * t6;

    de = f12 - f6;
    df *= (f12*SixF - f6*ThreeF) * TwoF * r2inv;

    if (needEnergy) {
      pE = (addToDVDL_vdw) ? cSim.pTILJ14[TIRegion] : cSim.pTISCLJ14[TIRegion];
      addEnergy(pE, de);

      if (addToDVDL_vdw && useSC14) {
        PMEFloat dt6_dl = scalpha * TISign * d_lambda_v * (BoverA) / (((scalpha * lambda_v) + rt6) * ((scalpha * lambda_v) + rt6));
        PMEFloat du_dt6 = g2 * (2 * A * t6 - B);
        addEnergy(cSim.pTISCLJ_DL[TIRegion], du_dt6 * dt6_dl);

        // MBAR Stuff
        if (cSim.needMBAR) {
          unsigned mbarShift = cSim.nMBARStates * (TIRegion + Schedule::TypeVDW * 3);
          unsigned long long int* MBAREnergy = &(cSim.pMBAREnergy[mbarShift]);
          PMEFloat* lambda[2] = { cSim.pMBARLambda , &(cSim.pMBARLambda[cSim.nMBARStates]) };

          unsigned i = 0;
          unsigned j = cSim.nMBARStates-1;
          if (cSim.currentMBARState >= 0 && cSim.rangeMBAR >= 0) {
            i = max(0, cSim.currentMBARState - cSim.rangeMBAR);
            j = min(cSim.nMBARStates - 1, cSim.currentMBARState + cSim.rangeMBAR);
          }
          for (unsigned l = i; l <= j; l++) {
            t6 = (vdwSCType == 1)
              ? BoverA / (scalpha * (smooth_step_func(lambda[TIRegion][l], cSim.vdwSmoothLambdaType)) + rt6)
              : BoverA / (scalpha * lambda[TIRegion][l] + rt6);
            PMEFloat dl_bar = g2 * t6 * (A * t6 - B);
            addEnergy(&(MBAREnergy[l]), dl_bar - de);
          }
        }
      }
    }

    // Add to the force accumulators
    if (addToDVDL_vdw) {
      pFx = cSim.pTINBForceX[TIRegion];
      pFy = cSim.pTINBForceY[TIRegion];
      pFz = cSim.pTINBForceZ[TIRegion];
    } else {
      pFx = cSim.pTINBSCForceX[TIRegion];
      pFy = cSim.pTINBSCForceY[TIRegion];
      pFz = cSim.pTINBSCForceZ[TIRegion];
    }

    ttx = df * dx;
    tty = df * dy;
    ttz = df * dz;
    if (addToDVDL_vdw) {
      ttx *= tiWeight_vdw;
      tty *= tiWeight_vdw;
      ttz *= tiWeight_vdw;
    }
    addForce(pFx + iatom0, -ttx);
    addForce(pFx + iatom1, ttx);
    addForce(pFy + iatom0, -tty);
    addForce(pFy + iatom1, tty);
    addForce(pFz + iatom0, -ttz);
    addForce(pFz + iatom1, ttz);
    if (needVirial && addToDVDL_vdw) {
      addVirial(cSim.pV11[TIRegion], -ttx * dx);
      addVirial(cSim.pV22[TIRegion], -tty * dy);
      addVirial(cSim.pV33[TIRegion], -ttz * dz);
    }

    pos += gridDim.x * blockDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kgCalculateREAF14NB_kernel:
//
// Arguments:
//   energy:
//   needVirial:
//---------------------------------------------------------------------------------------------F
_kPlainHead_ kgCalculateREAF14NB_kernel(bool needEnergy, bool needVirial){

  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;

  while (pos < cSim.numberREAF14NBEntries) {
    unsigned atom0 = cSim.pREAFNb14ID[pos].x;
    unsigned atom1 = cSim.pREAFNb14ID[pos].y;

    bool intRE = (cSim.pREAFList[atom0] > 0 && cSim.pREAFList[atom1] > 0);
    bool hasRE = (cSim.pREAFList[atom0] > 0 || cSim.pREAFList[atom1] > 0);

    Schedule::InteractionType tauType = (intRE) ? Schedule::TypeEleSC : Schedule::TypeEleCC;

    bool addRE_ele, addRE_vdw;
    switch (cSim.addRE) {
      case 0: addRE_ele = (!hasRE); addRE_vdw = (!hasRE); break;   // only non-RE part
      case 1: default:  addRE_ele = (!intRE); addRE_vdw = (!intRE); break; // add common-RE (ele+vdw) boundary
      case 2: addRE_ele = true; addRE_vdw = (!intRE); break;  // add SC internal terms (ele only)
      case 3: case 4: case 5: case 6: case 7: addRE_ele = true; addRE_vdw = true; break;  //  add all NB RE terms
      case -1: addRE_ele = (!intRE); addRE_vdw = (!hasRE); break; // add common-RE (ele only) boundary
    }

    // Note that the 1-4 terms aer already calculated by MD or TI; here are only corrections
    PMEFloat weight_ele = (addRE_ele) ? ((intRE) ? cSim.REAFItemWeight[tauType][1] : cSim.REAFItemWeight[tauType][0]) - OneF
      : ZeroF;
    PMEFloat weight_vdw = (addRE_vdw) ? ((intRE) ? cSim.REAFItemWeight[Schedule::TypeVDW][1] : cSim.REAFItemWeight[Schedule::TypeVDW][0]) -OneF
      : ZeroF;

    PMEFloat g1 = cSim.pREAFNb141[pos].x;

    PMEFloat f6 = ZeroF;
    PMEFloat f12 = ZeroF;
    unsigned iatom0 = cSim.pImageAtomLookup[atom0];
    unsigned iatom1 = cSim.pImageAtomLookup[atom1];
    PMEFloat dx = cSim.pImageX[iatom1] - cSim.pImageX[iatom0];
    PMEFloat dy = cSim.pImageY[iatom1] - cSim.pImageY[iatom0];
    PMEFloat dz = cSim.pImageZ[iatom1] - cSim.pImageZ[iatom0];
    PMEFloat r2 = dx * dx + dy * dy + dz * dz;
    PMEFloat rinv = Rsqrt(r2);
    PMEFloat r2inv = rinv * rinv;
    PMEFloat df = ZeroF;

    PMEFloat de = g1 * rinv;
    PMEFloat ddf = de * r2inv;

    if (needEnergy) addEnergy(cSim.pEEL14, de* weight_ele);
    df += ddf* weight_ele;

    // VDW part
    PMEFloat g2 = cSim.pREAFNb141[pos].y;
    PMEFloat A = cSim.pREAFNb142[pos].x;
    PMEFloat B = cSim.pREAFNb142[pos].y;

    if (A > delta && B > delta) {
      PMEFloat t6 = ZeroF;
      t6 = r2inv * r2inv * r2inv;
      f6 = g2 * B * t6;
      f12 = g2 * A * t6 * t6;
      de = f12 - f6;
      df += weight_vdw * (f12 * SixF - f6 * ThreeF) * TwoF * r2inv;

      if (needEnergy)  addEnergy(cSim.pENB14, de * weight_vdw);
    }

    PMEFloat ttx = df * dx;
    PMEFloat tty = df * dy;
    PMEFloat ttz = df * dz;

    addForce(cSim.pNBForceXAccumulator + iatom0, -ttx);
    addForce(cSim.pNBForceXAccumulator + iatom1, ttx);
    addForce(cSim.pNBForceYAccumulator + iatom0, -tty);
    addForce(cSim.pNBForceYAccumulator + iatom1, tty);
    addForce(cSim.pNBForceZAccumulator + iatom0, -ttz);
    addForce(cSim.pNBForceZAccumulator + iatom1, ttz);
    if (needVirial) {
      addVirial(cSim.pVirial_11, -ttx * dx);
      addVirial(cSim.pVirial_22, -tty * dy);
      addVirial(cSim.pVirial_33, -ttz * dz);
    }

    pos += gridDim.x * blockDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kgCalculate1264NB_kernel:
//
// Arguments:
//   needEnergy:
//   needVirial:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgCalculate1264NB_kernel(bool needEnergy, bool needVirial) {
  unsigned st = cSim.stride;
  PMEForceAccumulator* pFx = cSim.pNBForceAccumulator;
  PMEForceAccumulator* pFy = cSim.pNBForceAccumulator + st;
  PMEForceAccumulator* pFz = cSim.pNBForceAccumulator + st * 2;
  unsigned long long int* pV[3] = { cSim.pVirial_11, cSim.pVirial_22, cSim.pVirial_33 };

  __shared__ PMEFloat recip[9], ucell[9];
  __shared__ unsigned myNumberNB;
  unsigned i = threadIdx.x;

  if (cSim.pNTPData == NULL) {
    if (i < 3) {
      recip[i * 3] = cSim.recipf[i][0];
      recip[i * 3 + 1] = cSim.recipf[i][1];
      recip[i * 3 + 2] = cSim.recipf[i][2];
      ucell[i * 3] = cSim.ucellf[i][0];
      ucell[i * 3 + 1] = cSim.ucellf[i][1];
      ucell[i * 3 + 2] = cSim.ucellf[i][2];
    }
  } else {
    if (i < 9) {
      ucell[i] = cSim.pNTPData->ucellf[i];
      recip[i] = cSim.pNTPData->recipf[i];
    }
  }

  int bShift = blockIdx.x * gti_simulationConst::MaxNumberNBPerAtom;
  myNumberNB = cSim.pLJ1264NBList[(blockIdx.x + 1) * gti_simulationConst::MaxNumberNBPerAtom - 1].x;

  __syncthreads();

  unsigned int pos = threadIdx.x;

  PMEDouble*& pX = cSim.pImageX;
  PMEDouble*& pY = cSim.pImageY;
  PMEDouble*& pZ = cSim.pImageZ;

  while (pos < myNumberNB) {
    int4 tt = cSim.pLJ1264NBList[pos+ bShift];
    int& vdwIndex = tt.w;
    if (vdwIndex >= 0 && (tt.z & gti_simulationConst::ex_need_LJ)) {
      unsigned iatom0 = tt.x;
      unsigned iatom1 = tt.y;

      if (tt.x >= 0 && tt.y >= 0) {
        PMEFloat dx = pX[iatom1] - pX[iatom0];
        PMEFloat dy = pY[iatom1] - pY[iatom0];
        PMEFloat dz = pZ[iatom1] - pZ[iatom0];
        PMEFloat r2 = __image_dist2<PMEFloat>(dx, dy, dz, recip, ucell);

        if (r2 < cSim.cut2) {
          // Here we only deal with C4 terms
          PMEFloat rinv = Rsqrt(r2);
          PMEFloat r2inv = rinv * rinv;
          PMEFloat t4 = r2inv * r2inv;;
          PMEFloat f4 = cSim.pTIcn[vdwIndex].z * t4;

          PMEFloat df = -f4 * FourF * r2inv;
          if (needEnergy) {
            addEnergy((cSim.pEVDW), -f4);
          }
          PMEFloat ttx = df * dx;
          PMEFloat tty = df * dy;
          PMEFloat ttz = df * dz;

          // Add to the force accumulators
          addForce(pFx + iatom0, -ttx);
          addForce(pFx + iatom1, ttx);
          addForce(pFy + iatom0, -tty);
          addForce(pFy + iatom1, tty);
          addForce(pFz + iatom0, -ttz);
          addForce(pFz + iatom1, ttz);
          if (needVirial) {
            addVirial(pV[0], ttx * dx);
            addVirial(pV[1], tty * dy);
            addVirial(pV[2], ttz * dz);
          }
        }
      }
    }
    pos += blockDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kgCalculatep1264NB_kernel: C4PairwiseCUDA
//
// Arguments:
//   needEnergy:  
//   needVirial:  
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgCalculatep1264NB_kernel(bool needEnergy, bool needVirial) {
  unsigned st = cSim.stride;
  PMEForceAccumulator* pFx = cSim.pNBForceAccumulator;
  PMEForceAccumulator* pFy = cSim.pNBForceAccumulator + st;
  PMEForceAccumulator* pFz = cSim.pNBForceAccumulator + st * 2;
  unsigned long long int* pV[3] = { cSim.pVirial_11, cSim.pVirial_22, cSim.pVirial_33 };

  __shared__ PMEFloat recip[9], ucell[9];
  __shared__ unsigned myNumberNB;
  unsigned i = threadIdx.x;

  if (cSim.pNTPData == NULL) {
    if (i < 3) {
      recip[i * 3] = cSim.recipf[i][0];
      recip[i * 3 + 1] = cSim.recipf[i][1];
      recip[i * 3 + 2] = cSim.recipf[i][2];
      ucell[i * 3] = cSim.ucellf[i][0];
      ucell[i * 3 + 1] = cSim.ucellf[i][1];
      ucell[i * 3 + 2] = cSim.ucellf[i][2];
    }
  } else {
    if (i < 9) {
      ucell[i] = cSim.pNTPData->ucellf[i];
      recip[i] = cSim.pNTPData->recipf[i];
    }
  }

  int bShift = blockIdx.x * gti_simulationConst::MaxNumberNBPerAtom;
  myNumberNB = cSim.ppLJ1264NBList[(blockIdx.x + 1) * gti_simulationConst::MaxNumberNBPerAtom - 1].x;

  __syncthreads();

  unsigned int pos = threadIdx.x;

  PMEDouble*& pX = cSim.pImageX;
  PMEDouble*& pY = cSim.pImageY;
  PMEDouble*& pZ = cSim.pImageZ;

  while (pos < myNumberNB) {
    int4 tt = cSim.ppLJ1264NBList[pos+ bShift];
    int& vdwIndex = tt.w;
    if (vdwIndex >= 0 && (tt.z & gti_simulationConst::ex_need_LJ)) {
      unsigned iatom0 = tt.x;
      unsigned iatom1 = tt.y;

      if (tt.x >= 0 && tt.y >= 0) {
        PMEFloat dx = pX[iatom1] - pX[iatom0];
        PMEFloat dy = pY[iatom1] - pY[iatom0];
        PMEFloat dz = pZ[iatom1] - pZ[iatom0];
        PMEFloat r2 = __image_dist2<PMEFloat>(dx, dy, dz, recip, ucell);

        if (r2 < cSim.cut2) {
          // Here we only deal with C4 terms
          PMEFloat rinv = Rsqrt(r2);
          PMEFloat r2inv = rinv * rinv;
          PMEFloat t4 = r2inv * r2inv;;
          PMEFloat f4 = 0;
          if (vdwIndex == 11111) {
            f4 = 0;
            for (int i=0;i<=cSim.TIC4Pairwise;i++) {
              if ((cSim.pImageAtom[iatom0] == (cSim.pTIDcoef[i*2]) && cSim.pImageAtom[iatom1] == (cSim.pTIDcoef[i*2+1])) || (cSim.pImageAtom[iatom0] == (cSim.pTIDcoef[i*2+1]) && cSim.pImageAtom[iatom1] == (cSim.pTIDcoef[i*2]))) {
                f4 += cSim.pTIDvalue[i] * t4; // C4PairwiseCUDA2023
              }
            }
          }
          PMEFloat df = -f4 * FourF * r2inv;
          if (needEnergy) {
            addEnergy((cSim.pEVDW), -f4);
          }
          PMEFloat ttx = df * dx;
          PMEFloat tty = df * dy;
          PMEFloat ttz = df * dz;

          // Add to the force accumulators
          addForce(pFx + iatom0, -ttx);
          addForce(pFx + iatom1, ttx);
          addForce(pFy + iatom0, -tty);
          addForce(pFy + iatom1, tty);
          addForce(pFz + iatom0, -ttz);
          addForce(pFz + iatom1, ttz);
          if (needVirial) {
            addVirial(pV[0], ttx * dx);
            addVirial(pV[1], tty * dy);
            addVirial(pV[2], ttz * dz);
          }
        }
      }
    }
    pos += blockDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kgCalculate1264p1264NB_kernel: C4PairwiseCUDA
// TODO: C4PaireiseCUDA2023
// Arguments:
//   needEnergy:  
//   needVirial:  
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgCalculate1264p1264NB_kernel(bool needEnergy, bool needVirial) {
  unsigned st = cSim.stride;
  PMEForceAccumulator* pFx = cSim.pNBForceAccumulator;
  PMEForceAccumulator* pFy = cSim.pNBForceAccumulator + st;
  PMEForceAccumulator* pFz = cSim.pNBForceAccumulator + st * 2;
  unsigned long long int* pV[3] = { cSim.pVirial_11, cSim.pVirial_22, cSim.pVirial_33 };

  __shared__ PMEFloat recip[9], ucell[9];
  __shared__ unsigned myNumberNB;
  unsigned i = threadIdx.x;

  if (cSim.pNTPData == NULL) {
    if (i < 3) {
      recip[i * 3] = cSim.recipf[i][0];
      recip[i * 3 + 1] = cSim.recipf[i][1];
      recip[i * 3 + 2] = cSim.recipf[i][2];
      ucell[i * 3] = cSim.ucellf[i][0];
      ucell[i * 3 + 1] = cSim.ucellf[i][1];
      ucell[i * 3 + 2] = cSim.ucellf[i][2];
    }
  } else {
    if (i < 9) {
      ucell[i] = cSim.pNTPData->ucellf[i];
      recip[i] = cSim.pNTPData->recipf[i];
    }
  }

  int bShift = blockIdx.x * gti_simulationConst::MaxNumberNBPerAtom;
  myNumberNB = cSim.pLJ1264pLJ1264NBList[(blockIdx.x + 1) * gti_simulationConst::MaxNumberNBPerAtom - 1].x;

  __syncthreads();

  unsigned int pos = threadIdx.x;

  PMEDouble*& pX = cSim.pImageX;
  PMEDouble*& pY = cSim.pImageY;
  PMEDouble*& pZ = cSim.pImageZ;


  while (pos < myNumberNB) { //C4PairwiseCUDA2023
    int4 tt = cSim.pLJ1264pLJ1264NBList[pos+ bShift];
    int& vdwIndex = tt.w;
    if (vdwIndex >= 0 && (tt.z & gti_simulationConst::ex_need_LJ)) {
      unsigned iatom0 = tt.x;
      unsigned iatom1 = tt.y;

      if (tt.x >= 0 && tt.y >= 0) {
        PMEFloat dx = pX[iatom1] - pX[iatom0];
        PMEFloat dy = pY[iatom1] - pY[iatom0];
        PMEFloat dz = pZ[iatom1] - pZ[iatom0];
        PMEFloat r2 = __image_dist2<PMEFloat>(dx, dy, dz, recip, ucell);

        if (r2 < cSim.cut2) {
          // Here we only deal with C4 terms
          PMEFloat rinv = Rsqrt(r2);
          PMEFloat r2inv = rinv * rinv;
          PMEFloat t4 = r2inv * r2inv;
          PMEFloat f4 = cSim.pTIcn[vdwIndex].z * t4;
          if (vdwIndex == 55555) {
            //cSim.pLJ1264pLJ1264NBList[pos+ bShift].w = 55556; //C4PairwiseCUDA
            for (int i=0;i<cSim.TIC4Pairwise;i++) {
              if ((cSim.pImageAtom[iatom0] == (cSim.pTIDcoef[i*2]) && cSim.pImageAtom[iatom1] == (cSim.pTIDcoef[i*2+1])) || (cSim.pImageAtom[iatom0] == (cSim.pTIDcoef[i*2+1]) && cSim.pImageAtom[iatom1] == (cSim.pTIDcoef[i*2]))) {
                f4 += cSim.pTIDvalue[i] * t4 / 3.0; // C4PairwiseCUDA2023
              }
            }
          }
          PMEFloat df = -f4 * FourF * r2inv;
          if (needEnergy) {
            addEnergy((cSim.pEVDW), -f4);
          }
          PMEFloat ttx = df * dx;
          PMEFloat tty = df * dy;
          PMEFloat ttz = df * dz;

          // Add to the force accumulators
          addForce(pFx + iatom0, -ttx);
          addForce(pFx + iatom1, ttx);
          addForce(pFy + iatom0, -tty);
          addForce(pFy + iatom1, tty);
          addForce(pFz + iatom0, -ttz);
          addForce(pFz + iatom1, ttz);
          if (needVirial) {
            addVirial(pV[0], ttx * dx);
            addVirial(pV[1], tty * dy);
            addVirial(pV[2], ttz * dz);
          }
        }
      }
    }
    pos += blockDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kgCalculateTINB_gamd_kernel:
//
// Arguments:
//   needEnergy:
//   needVirial:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgCalculateTINB_gamd_kernel(bool needEnergy, bool needVirial){

  unsigned i = threadIdx.x;
  __shared__ PMEFloat recip[9], ucell[9];

  if (cSim.pNTPData == NULL) {
    if (i < 3) {
      recip[i * 3] = cSim.recipf[i][0];
      recip[i * 3 + 1] = cSim.recipf[i][1];
      recip[i * 3 + 2] = cSim.recipf[i][2];
      ucell[i * 3] = cSim.ucellf[i][0];
      ucell[i * 3 + 1] = cSim.ucellf[i][1];
      ucell[i * 3 + 2] = cSim.ucellf[i][2];
    }
  } else {
    if (i < 9) {
      ucell[i] = cSim.pNTPData->ucellf[i];
      recip[i] = cSim.pNTPData->recipf[i];
    }
  }

  __syncthreads();

  PMEFloat4* cn = cSim.pTIcn;
  PMEFloat4* sigEps = cSim.pTISigEps;

  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;

  while (pos < *cSim.pNumberTINBEntries) {

    int4 tt = cSim.pTINBList[pos];

    unsigned iatom0 = cSim.pImageAtomLookup[tt.x];
    unsigned iatom1 = cSim.pImageAtomLookup[tt.y];

    PMEFloat dx = cSim.pImageX[iatom1] - cSim.pImageX[iatom0];
    PMEFloat dy = cSim.pImageY[iatom1] - cSim.pImageY[iatom0];
    PMEFloat dz = cSim.pImageZ[iatom1] - cSim.pImageZ[iatom0];
    PMEFloat r2 = __image_dist2<PMEFloat>(dx, dy, dz, recip, ucell);

    bool excluded = tt.z & 4;

    if (r2 < cSim.cut2 || excluded) {

      int TIRegion = tt.z & gti_simulationConst::Fg_TI_region;
      int TISign = ((TIRegion == 0) ? -1 : 1);

      bool needLJ = (tt.w >= 0 && (tt.z & 2));
      bool needSCLJ = tt.z & 8;
      bool needSCEle = tt.z & 16;

      bool isIntSC = (cSim.pSCList[tt.x] > 0 && cSim.pSCList[tt.y] > 0);
      bool hasSC   = (cSim.pSCList[tt.x] > 0 || cSim.pSCList[tt.y] > 0);
      //bool addToDVDL_ele = (!hasSC);
      //bool addToDVDL_vdw = (!hasSC);

      unsigned long long int* pE, * pEDL;
      PMEFloat dE = ZeroF;
      PMEFloat df = ZeroF;
      PMEFloat dfCorr = ZeroF;
      PMEFloat df_gamd = ZeroF;

      int& vdwIndex = tt.w;
      bool forceSC = (isIntSC && !excluded);
      //bool simpleElec = (!needSCEle || isIntSC || excluded);

      // EED
      PMEFloat q01 = cSim.pOrigAtomCharge[tt.x] * cSim.pOrigAtomCharge[tt.y];
      PMEFloat rinv = Rsqrt(r2);
      PMEFloat r = r2 * rinv;
      PMEFloat swtch = excluded ? -Erf(cSim.ew_coeffSP * r) : Erfc(cSim.ew_coeffSP * r);
      PMEFloat d_swtch = -cSim.negTwoEw_coeffRsqrtPI * Exp(-cSim.ew_coeff2 * r2);
      PMEFloat g1 = q01 * swtch;
      PMEFloat r2inv = rinv*rinv;
      PMEFloat den = 0.0;
      PMEFloat den_d = 0.0;

      pE = (isIntSC) ? cSim.pTISCEED[TIRegion] : cSim.pTIEECC[TIRegion];
      pEDL = (isIntSC) ? cSim.pTISCEESC_DL[TIRegion] : cSim.pTISCEECC_DL[TIRegion];

      if (!needSCEle || excluded) {
        dE = g1 *rinv;
        df = Fma(d_swtch, q01, dE) * r2inv;
	addEnergy(cSim.pTIEECC[2], dE);  // EXCLUDED ELEC TERMS ARE NOT INCLUDED FOR ADDING GAMD BOOST POTENTIAL.
	// addEnergy(cSim.pTIEESC[TIRegion], dE);  // EXCLUDED ELEC TERMS ARE NOT INCLUDED FOR ADDING GAMD BOOST POTENTIAL.
      } else if (isIntSC) {
        dE = g1 *rinv;
        dfCorr = Fma(d_swtch, q01, dE) * r2inv;
        if (cSim.igamd == 6 || cSim.igamd == 7 || cSim.igamd == 10 || cSim.igamd == 11 || cSim.igamd == 100 || cSim.igamd == 103) {
          df_gamd = dfCorr;
        }
	addEnergy(pE, dE);
	// addEnergy(cSim.pTISCEED[TIRegion], dE);
      } else if (hasSC) {
	den = Rsqrt(r2 + cSim.TISCBeta*cSim.TILambda[TIRegion]);
	den_d = Pow(den, ThreeF);
	df = (q01 * d_swtch * den * rinv) + (g1 * den_d);
	df *= (1.0 - cSim.TILambda[TIRegion]);
	dE = (1.0 - cSim.TILambda[TIRegion]) * g1 * den;
        if (cSim.igamd == 6 || cSim.igamd == 7 || cSim.igamd == 10 || cSim.igamd == 11 || cSim.igamd == 100 || cSim.igamd == 103) {
          df_gamd = df;
	}
	addEnergy(pE, dE);
	// addEnergy(cSim.pTIEECC[TIRegion], dE);
      } else {
        dE = g1 *rinv;
        df = Fma(d_swtch, q01, dE) * r2inv;
        if (cSim.igamd == 6 || cSim.igamd == 7 || cSim.igamd == 10 || cSim.igamd == 11 || cSim.igamd == 100 || cSim.igamd == 103) {
          df_gamd = df;
	}
	addEnergy(pE, dE);
	// addEnergy(cSim.pTIEECC[TIRegion], dE);
      }

      if (!isIntSC) { // should be always true here...
        if (needEnergy) {
      		      //lambda-independent contribution
	              den = Rsqrt(r2 + cSim.TISCBeta*cSim.TILambda[TIRegion]);
	              den_d = Pow(den, ThreeF);
      		      PMEFloat ss = g1 * den;
      		      //lambda-dependent contribution
      		      PMEFloat dl = - cSim.TILambda[TIRegion] * ss ;
      		      dl += 0.5*(1.0 - cSim.TILambda[TIRegion]) * g1*den_d*cSim.TISCBeta * TISign;
      		      addEnergy(pEDL, dl*(1.0 - cSim.TILambda[TIRegion]));
		      // addEnergy(cSim.pTISCEESC_DL[TIRegion], dl*(1.0 - cSim.TILambda[TIRegion]));

          // MBAR Stuff
          if (cSim.needMBAR) {
		  PMEFloat* lambda[2] = { cSim.pMBARLambda, &(cSim.pMBARLambda[cSim.nMBARStates]) };
      			      for (unsigned l = 0; l < cSim.nMBARStates; l++) {
      				      PMEFloat den_n = Rsqrt(r2 + cSim.TISCBeta*lambda[TIRegion][l]);
      				      addEnergy(&(cSim.pMBAREnergy[l]),
      					      Pow((OneF - lambda[TIRegion][l]), cSim.TISCGamma)
      					      * g1 * (den_n - den));
      			      }
          }  // if (cSim.needMBAR) {
        }  // Need Energy
      } // if (!isIntSC)

      // Add to the force accumulators
      PMEFloat ttx = df * dx;
      PMEFloat tty = df * dy;
      PMEFloat ttz = df * dz;
      addForce(cSim.pTINBForceX[TIRegion] + iatom0, -ttx);
      addForce(cSim.pTINBForceX[TIRegion] + iatom1, ttx);
      addForce(cSim.pTINBForceY[TIRegion] + iatom0, -tty);
      addForce(cSim.pTINBForceY[TIRegion] + iatom1, tty);
      addForce(cSim.pTINBForceZ[TIRegion] + iatom0, -ttz);
      addForce(cSim.pTINBForceZ[TIRegion] + iatom1, ttz);
      if (needVirial) {
        addVirial(cSim.pV11[TIRegion], -ttx * dx);
        addVirial(cSim.pV22[TIRegion], -tty * dy);
        addVirial(cSim.pV33[TIRegion], -ttz * dz);
      }
      if (forceSC) {
        ttx = dfCorr * dx;
        tty = dfCorr * dy;
        ttz = dfCorr * dz;
        addForce(cSim.pTINBSCForceX[TIRegion] + iatom0, -ttx);
        addForce(cSim.pTINBSCForceX[TIRegion] + iatom1, ttx);
        addForce(cSim.pTINBSCForceY[TIRegion] + iatom0, -tty);
        addForce(cSim.pTINBSCForceY[TIRegion] + iatom1, tty);
        addForce(cSim.pTINBSCForceZ[TIRegion] + iatom0, -ttz);
        addForce(cSim.pTINBSCForceZ[TIRegion] + iatom1, ttz);
        if (needVirial) {
          addVirial(cSim.pVirial_11, -ttx * dx);
          addVirial(cSim.pVirial_22, -tty * dy);
          addVirial(cSim.pVirial_33, -ttz * dz);
        }
      }

      // record GaMDTINBForce
      ttx = df_gamd * dx;
      tty = df_gamd * dy;
      ttz = df_gamd * dz;
      addForce(cSim.pGaMDTINBForceX[TIRegion] + iatom0, -ttx);
      addForce(cSim.pGaMDTINBForceX[TIRegion] + iatom1, ttx);
      addForce(cSim.pGaMDTINBForceY[TIRegion] + iatom0, -tty);
      addForce(cSim.pGaMDTINBForceY[TIRegion] + iatom1, tty);
      addForce(cSim.pGaMDTINBForceZ[TIRegion] + iatom0, -ttz);
      addForce(cSim.pGaMDTINBForceZ[TIRegion] + iatom1, ttz);

      df=ZeroF;
      dfCorr=ZeroF;
      df_gamd=ZeroF;
      forceSC = (isIntSC && !excluded);
      // EVDW
      if (needLJ && !excluded) {
        PMEFloat t4, f4, t6, f6, f12, dd, dl, u1_3;
        pE = (!isIntSC) ? cSim.pTILJ[TIRegion] : cSim.pTISCLJ[TIRegion];
        if (!needSCLJ || isIntSC) {
          t4 = r2inv * r2inv;
          f4 = (sigEps[vdwIndex].w > 0) ? cn[vdwIndex].z * t4 : ZeroF;
          t6 = t4 * r2inv;
          f6 = cn[vdwIndex].y * t6;
          f12 = cn[vdwIndex].x * t6 * t6;
          dd = (f12*TwelveF - f6*SixF - f4*FourF) * r2inv;
          if (isIntSC) {
            dfCorr += dd;
          } else {
            df += dd;
          }
          if (needEnergy) {
            addEnergy(pE, (f12 - f6 - f4));
          }
          if (cSim.igamd == 6 || cSim.igamd == 7 || cSim.igamd == 10 || cSim.igamd == 11 || cSim.igamd == 100 || cSim.igamd == 104) {
              df_gamd = dd;
          }
        } else {
        if (vdwIndex == 55555) {
          for (int i=0;i<cSim.TIC4Pairwise;i++) {
            if ((cSim.pImageAtom[tt.x] == (cSim.pTIDcoef[i*2]) && cSim.pImageAtom[tt.y] == (cSim.pTIDcoef[i*2+1])) || (cSim.pImageAtom[tt.x] == (cSim.pTIDcoef[i*2+1]) && cSim.pImageAtom[tt.y] == (cSim.pTIDcoef[i*2]))) {
              f4 += cSim.pTIDvalue[i] * t4 / 3.0; // C4PairwiseCUDA2023
              //eLJ -= f4;
              //dvdu += (TwoF * t4) / (ThreeF * u1_3);
            }
          }
        }
          t4 = sigEps[vdwIndex].z;
          t6 = r2 * r2 * r2 * sigEps[vdwIndex].x;
          f6 = OneF / (cSim.TISCAlpha*cSim.TILambda[TIRegion] + t6);
          f12 = f6 * f6;
          dl = sigEps[vdwIndex].y * (f12 - f6);
          PMEFloat dvdu = sigEps[vdwIndex].y * (OneF - TwoF*f6);

          if (sigEps[vdwIndex].w > 0) {
            u1_3 = cbrt(f6);
            f4 = t4 * u1_3 * u1_3;
            dl -= f4;
            dvdu += (TwoF * t4) / (ThreeF * u1_3);
          } 
          if (vdwIndex == 55555) {
            f4 = 0;
            for (int i=0;i<cSim.TIC4Pairwise;i++) {
              if ((cSim.pImageAtom[tt.x] == (cSim.pTIDcoef[i*2]) && cSim.pImageAtom[tt.y] == (cSim.pTIDcoef[i*2+1])) || (cSim.pImageAtom[tt.x] == (cSim.pTIDcoef[i*2+1]) && cSim.pImageAtom[tt.y] == (cSim.pTIDcoef[i*2]))) {
                f4 += cSim.pTIDvalue[i] * t4 / 3.0; // C4PairwiseCUDA2023
                //eLJ -= f4;
                dvdu += (TwoF * t4) / (ThreeF * u1_3);
              }
            }
          }

          if (needEnergy) {
            addEnergy(pE, dl);
          }
          PMEFloat dudr = SixF * f12 * t6 * r2inv;
          dd = dvdu * dudr;
          df -= dd;
          if (cSim.igamd == 6 || cSim.igamd == 7 || cSim.igamd == 10 || cSim.igamd == 11 || cSim.igamd == 100 || cSim.igamd == 104) {
              df_gamd = dd;
          }
          if (!isIntSC) {
            PMEFloat dudl = -cSim.TISCAlpha * f12 * TISign;
            if (needEnergy) {
              addEnergy(cSim.pTISCLJ_DL[TIRegion], dvdu*dudl);
            }
          }

          // MBAR Stuff
          if (cSim.needMBAR && needEnergy) {
	    PMEFloat* lambda[2] = { cSim.pMBARLambda, &(cSim.pMBARLambda[cSim.nMBARStates]) };
            for (unsigned l = 0; l < cSim.nMBARStates; l++) {
              f6 = OneF / (cSim.TISCAlpha * lambda[TIRegion][l] + t6);
              f12 = f6 * f6;
              PMEFloat dl_bar = sigEps[vdwIndex].y * (f12 - f6);
              if (sigEps[vdwIndex].w > 0) {
                u1_3 = cbrt(f6);
                dl_bar -= t4 * u1_3 * u1_3;
              }
              addEnergy(&(cSim.pMBAREnergy[l]),
                (OneF - lambda[TIRegion][l]) * (dl_bar - dl));
            }
          }

        }
      }

      // record GaMDTINBForce
      ttx = df_gamd * dx;
      tty = df_gamd * dy;
      ttz = df_gamd * dz;
      addForce(cSim.pGaMDTINBForceX[TIRegion] + iatom0, -ttx);
      addForce(cSim.pGaMDTINBForceX[TIRegion] + iatom1, ttx);
      addForce(cSim.pGaMDTINBForceY[TIRegion] + iatom0, -tty);
      addForce(cSim.pGaMDTINBForceY[TIRegion] + iatom1, tty);
      addForce(cSim.pGaMDTINBForceZ[TIRegion] + iatom0, -ttz);
      addForce(cSim.pGaMDTINBForceZ[TIRegion] + iatom1, ttz);

      // Add to the force accumulators
      ttx = df * dx;
      tty = df * dy;
      ttz = df * dz;
      addForce(cSim.pTINBForceX[TIRegion] + iatom0, -ttx);
      addForce(cSim.pTINBForceX[TIRegion] + iatom1, ttx);
      addForce(cSim.pTINBForceY[TIRegion] + iatom0, -tty);
      addForce(cSim.pTINBForceY[TIRegion] + iatom1, tty);
      addForce(cSim.pTINBForceZ[TIRegion] + iatom0, -ttz);
      addForce(cSim.pTINBForceZ[TIRegion] + iatom1, ttz);
      if (needVirial) {
        addVirial(cSim.pV11[TIRegion], -ttx * dx);
        addVirial(cSim.pV22[TIRegion], -tty * dy);
        addVirial(cSim.pV33[TIRegion], -ttz * dz);
      }
      if (forceSC) {
        ttx = dfCorr * dx;
        tty = dfCorr * dy;
        ttz = dfCorr * dz;
        addForce(cSim.pTINBSCForceX[TIRegion] + iatom0, -ttx);
        addForce(cSim.pTINBSCForceX[TIRegion] + iatom1, ttx);
        addForce(cSim.pTINBSCForceY[TIRegion] + iatom0, -tty);
        addForce(cSim.pTINBSCForceY[TIRegion] + iatom1, tty);
        addForce(cSim.pTINBSCForceZ[TIRegion] + iatom0, -ttz);
        addForce(cSim.pTINBSCForceZ[TIRegion] + iatom1, ttz);
        if (needVirial) {
          addVirial(cSim.pVirial_11, -ttx * dx);
          addVirial(cSim.pVirial_22, -tty * dy);
          addVirial(cSim.pVirial_33, -ttz * dz);
        }
      }
    }
    pos += gridDim.x * blockDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kgCalculateTI14NB_gamd_kernel:
//
// Arguments:
//   energy:
//   needVirial:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgCalculateTI14NB_gamd_kernel(bool energy, bool needVirial) {
  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;

  unsigned long long int *pE0, *pE;
  PMEAccumulator* pFx, *pFy, *pFz;
  bool useSC14 = false;

  while (pos < cSim.numberTI14NBEntries) {
    unsigned atom0 = cSim.pTINb14ID[pos].x;
    unsigned atom1 = cSim.pTINb14ID[pos].y;
    unsigned TIRegion = cSim.pTINb14ID[pos].z;
    int TISign = (TIRegion * 2) - 1;
    PMEDouble g1 = cSim.pTINb141[pos].x;
    PMEDouble g2 = cSim.pTINb141[pos].y;
    PMEDouble f6 = cSim.pTINb142[pos].y;
    PMEDouble f12 = cSim.pTINb142[pos].x;
    unsigned iatom0 = cSim.pImageAtomLookup[atom0];
    unsigned iatom1 = cSim.pImageAtomLookup[atom1];
    PMEDouble dx = cSim.pImageX[iatom1] - cSim.pImageX[iatom0];
    PMEDouble dy = cSim.pImageY[iatom1] - cSim.pImageY[iatom0];
    PMEDouble dz = cSim.pImageZ[iatom1] - cSim.pImageZ[iatom0];
    PMEDouble r2 = dx*dx + dy*dy + dz*dz;
    PMEDouble rinv = Rsqrt(r2);
    PMEDouble r2inv = rinv*rinv;
    PMEDouble t6 = ZeroF;

    bool isIntSC = (cSim.pSCList[atom0] > 0 && cSim.pSCList[atom1] > 0);
    bool hasSC = (cSim.pSCList[atom0] > 0 || cSim.pSCList[atom1] > 0);
    //bool partSC = hasSC && (!isIntSC);
    bool addToDVDL_ele = false, addToDVDL_vdw = false;

    switch (cSim.addSC) {
      case 0: addToDVDL_ele = (!hasSC); addToDVDL_vdw = (!hasSC); break;   // only non-SC part
	  case 1: addToDVDL_ele = (!isIntSC); addToDVDL_vdw = (!isIntSC); break; // add common-SC (ele+vdw) boundary
	  case 2: addToDVDL_ele = true; addToDVDL_vdw = (!isIntSC); break;  // add SC internal terms (ele only)
	  case 3: addToDVDL_ele = true; addToDVDL_vdw = true; break;  //  add all NB SC terms
	  //case 3: addToDVDL_ele = true; addToDVDL_vdw = true; useSC14 = true;  break;  //  add all NB SC terms
	  case -1: addToDVDL_ele = (!isIntSC); addToDVDL_vdw = (!hasSC); break; // add common-SC (ele only) boundary
     }

    // EED: note that 1:4 electrostatics are double-counted
    PMEDouble de = 0.0;
    PMEDouble ddf = 0.0;
    PMEDouble df = ZeroF;
    PMEDouble df_gamd = ZeroF;
    if (useSC14) {
      PMEDouble den = Rsqrt(r2 + (cSim.TISCBeta * cSim.TILambda[TIRegion]));
      de = g1 * den;
      PMEDouble den_n = powf(den, 3.0);
      ddf = g1 * den_n;
      de = 0.5 * g1 * den_n * cSim.TISCBeta * TISign;
    } else {
      de = g1 * rinv;
      ddf = de * r2inv;
    }
    pE = (addToDVDL_ele) ? cSim.pTIEE14CC[TIRegion] : cSim.pTISCEE14[TIRegion];
    addEnergy(pE, de);
    pE0 = cSim.pTIEE14CC[2];
    // if (addToDVDL_ele) {
      addEnergy(pE0, -de);
    // }
    df += ddf;

    // Record GaMDTINBForce: 1-4 ele
    // if (cSim.igamd == 6 || cSim.igamd == 7 || cSim.igamd == 10 || cSim.igamd == 11 || cSim.igamd == 100 || cSim.igamd == 101) { 
    if (cSim.igamd == 101) { 
	df_gamd += ddf;
    }

    if (addToDVDL_ele) {
      pFx = cSim.pTINBForceX[TIRegion];
      pFy = cSim.pTINBForceY[TIRegion];
      pFz = cSim.pTINBForceZ[TIRegion];
    } else {
      pFx = cSim.pTINBSCForceX[TIRegion];
      pFy = cSim.pTINBSCForceY[TIRegion];
      pFz = cSim.pTINBSCForceZ[TIRegion];
    }

    PMEFloat ttx = df * dx;
    PMEFloat tty = df * dy;
    PMEFloat ttz = df * dz;
    addForce(pFx + iatom0, -ttx);
    addForce(pFx + iatom1, ttx);
    addForce(pFy + iatom0, -tty);
    addForce(pFy + iatom1, tty);
    addForce(pFz + iatom0, -ttz);
    addForce(pFz + iatom1, ttz);
    if (needVirial && addToDVDL_ele) {
      addVirial(cSim.pV11[TIRegion], -ttx * dx);
      addVirial(cSim.pV22[TIRegion], -tty * dy);
      addVirial(cSim.pV33[TIRegion], -ttz * dz);
    }

    // Record GaMDTINBForce
    ttx = df_gamd * dx;
    tty = df_gamd * dy;
    ttz = df_gamd * dz;
    addForce(cSim.pGaMDTINBForceX[TIRegion] + iatom0, -ttx);
    addForce(cSim.pGaMDTINBForceX[TIRegion] + iatom1, ttx);
    addForce(cSim.pGaMDTINBForceY[TIRegion] + iatom0, -tty);
    addForce(cSim.pGaMDTINBForceY[TIRegion] + iatom1, tty);
    addForce(cSim.pGaMDTINBForceZ[TIRegion] + iatom0, -ttz);
    addForce(cSim.pGaMDTINBForceZ[TIRegion] + iatom1, ttz);

    // Correction for double-counted 1-4 ele 
    ttx = df * dx;
    tty = df * dy;
    ttz = df * dz;
    addForce(cSim.pTINBForceX[2] + iatom0, ttx);
    addForce(cSim.pTINBForceX[2] + iatom1, -ttx);
    addForce(cSim.pTINBForceY[2] + iatom0, tty);
    addForce(cSim.pTINBForceY[2] + iatom1, -tty);
    addForce(cSim.pTINBForceZ[2] + iatom0, ttz);
    addForce(cSim.pTINBForceZ[2] + iatom1, -ttz);
    if (needVirial) {
      addVirial(cSim.pVirial_11, ttx * dx);
      addVirial(cSim.pVirial_22, tty * dy);
      addVirial(cSim.pVirial_33, ttz * dz);
    }

    // GaMDTINBForce: Correction for double-counted 1-4 ele ??  No CC atoms in LiGaMD.
    // check pTINBForceY[2], pTIEE14CC[2]

    df = OneF; // df must be 1.0 here
    // EVDW: note that 14 NB is not double-counted since it is zeroed-out
    t6 = r2inv * r2inv * r2inv;
    f6 *= g2 * t6;
    f12 *= g2 * t6 * t6;
    de = f12 - f6;

    if (useSC14) {
      PMEDouble sig_6 = f6 / f12;
      PMEDouble fourEps = f6 * sig_6;
      t6 = r2 * r2 * r2 * sig_6;
      f6 = OneF / ((cSim.TISCAlpha * cSim.TILambda[TIRegion]) + t6);
      f12 = f6 * f6;
      df = fourEps * f6 * t6;
    }
    df *= (f12*SixF - f6*ThreeF) * TwoF * r2inv;
    pE = (addToDVDL_vdw) ? cSim.pTILJ14[TIRegion] : cSim.pTISCLJ14[TIRegion];
    addEnergy(pE, de);

    // Add to the force accumulators
    if (addToDVDL_vdw) {
      pFx = cSim.pTINBForceX[TIRegion];
      pFy = cSim.pTINBForceY[TIRegion];
      pFz = cSim.pTINBForceZ[TIRegion];
    } else {
      pFx = cSim.pTINBSCForceX[TIRegion];
      pFy = cSim.pTINBSCForceY[TIRegion];
      pFz = cSim.pTINBSCForceZ[TIRegion];
    }

    ttx = df * dx;
    tty = df * dy;
    ttz = df * dz;
    addForce(pFx + iatom0, -ttx);
    addForce(pFx + iatom1, ttx);
    addForce(pFy + iatom0, -tty);
    addForce(pFy + iatom1, tty);
    addForce(pFz + iatom0, -ttz);
    addForce(pFz + iatom1, ttz);
    if (needVirial && addToDVDL_vdw) {
      addVirial(cSim.pV11[TIRegion], -ttx * dx);
      addVirial(cSim.pV22[TIRegion], -tty * dy);
      addVirial(cSim.pV33[TIRegion], -ttz * dz);
    }

    // Record GaMDTINBForce: 1-4 vdw
    // if (cSim.igamd == 6 || cSim.igamd == 7 || cSim.igamd == 10 || cSim.igamd == 11 || cSim.igamd == 100 || cSim.igamd == 102) { 
    if (cSim.igamd == 102) { 
	df_gamd = df;
    }
    ttx = df_gamd * dx;
    tty = df_gamd * dy;
    ttz = df_gamd * dz;
    addForce(cSim.pGaMDTINBForceX[TIRegion] + iatom0, -ttx);
    addForce(cSim.pGaMDTINBForceX[TIRegion] + iatom1, ttx);
    addForce(cSim.pGaMDTINBForceY[TIRegion] + iatom0, -tty);
    addForce(cSim.pGaMDTINBForceY[TIRegion] + iatom1, tty);
    addForce(cSim.pGaMDTINBForceZ[TIRegion] + iatom0, -ttz);
    addForce(cSim.pGaMDTINBForceZ[TIRegion] + iatom1, ttz);

    pos += gridDim.x * blockDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kgCalculate1264NB_gamd_kernel:
//
// Arguments:
//   needEnergy:  
//   needVirial:  
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgCalculate1264NB_gamd_kernel(bool needEnergy, bool needVirial) {
  unsigned st = cSim.stride;
  PMEForceAccumulator* pFx = cSim.pNBForceAccumulator;
  PMEForceAccumulator* pFy = cSim.pNBForceAccumulator + st;
  PMEForceAccumulator* pFz = cSim.pNBForceAccumulator + st * 2;
  unsigned long long int* pV[3] = { cSim.pVirial_11, cSim.pVirial_22, cSim.pVirial_33 };

  __shared__ PMEFloat recip[9], ucell[9];
  __shared__ unsigned myNumberNB;
  unsigned i = threadIdx.x;

  if (cSim.pNTPData == NULL) {
    if (i < 3) {
      recip[i * 3] = cSim.recipf[i][0];
      recip[i * 3 + 1] = cSim.recipf[i][1];
      recip[i * 3 + 2] = cSim.recipf[i][2];
      ucell[i * 3] = cSim.ucellf[i][0];
      ucell[i * 3 + 1] = cSim.ucellf[i][1];
      ucell[i * 3 + 2] = cSim.ucellf[i][2];
    }
  } else {
    if (i < 9) {
      ucell[i] = cSim.pNTPData->ucellf[i];
      recip[i] = cSim.pNTPData->recipf[i];
    }
  }

  int bShift = blockIdx.x * gti_simulationConst::MaxNumberNBPerAtom;
  myNumberNB = cSim.pLJ1264NBList[(blockIdx.x + 1) * gti_simulationConst::MaxNumberNBPerAtom - 1].x;

  __syncthreads();

  unsigned int pos = threadIdx.x;

  PMEDouble*& pX = cSim.pImageX;
  PMEDouble*& pY = cSim.pImageY;
  PMEDouble*& pZ = cSim.pImageZ;

  while (pos < myNumberNB) {
    int4 tt = cSim.pLJ1264NBList[pos+ bShift];
    int& vdwIndex = tt.w;
    if (vdwIndex >= 0 && (tt.z & gti_simulationConst::ex_need_LJ)) {
      unsigned iatom0 = tt.x;
      unsigned iatom1 = tt.y;

      if (tt.x >= 0 && tt.y >= 0) {
        PMEFloat dx = pX[iatom1] - pX[iatom0];
        PMEFloat dy = pY[iatom1] - pY[iatom0];
        PMEFloat dz = pZ[iatom1] - pZ[iatom0];
        PMEFloat r2 = __image_dist2<PMEFloat>(dx, dy, dz, recip, ucell);

        if (r2 < cSim.cut2) {
          // Here we only deal with C4 terms
          PMEFloat rinv = Rsqrt(r2);
          PMEFloat r2inv = rinv * rinv;
          PMEFloat t4 = r2inv * r2inv;;
          PMEFloat f4 = cSim.pTIcn[vdwIndex].z * t4;

          PMEFloat df = -f4 * FourF * r2inv;
          if (needEnergy) {
            addEnergy((cSim.pEVDW), -f4);
          }

	  /*
          // record GaMDTINBForce
          if (cSim.igamd == 6 || cSim.igamd == 7 || cSim.igamd == 10 || cSim.igamd == 11 || cSim.igamd == 100 || cSim.igamd == 104) {
              PMEFloat df_gamd = df;
          }

          ttx = df_gamd * dx;
          tty = df_gamd * dy;
          ttz = df_gamd * dz;
          addForce(cSim.pGaMDTINBForceX[TIRegion] + iatom0, -ttx);	// not yet ready for pGaMDTINBForceX[TIRegion] with 1264NB
          addForce(cSim.pGaMDTINBForceX[TIRegion] + iatom1, ttx);
          addForce(cSim.pGaMDTINBForceY[TIRegion] + iatom0, -tty);
          addForce(cSim.pGaMDTINBForceY[TIRegion] + iatom1, tty);
          addForce(cSim.pGaMDTINBForceZ[TIRegion] + iatom0, -ttz);
          addForce(cSim.pGaMDTINBForceZ[TIRegion] + iatom1, ttz);
	  */

          PMEFloat ttx = df * dx;
          PMEFloat tty = df * dy;
          PMEFloat ttz = df * dz;

          // Add to the force accumulators
          addForce(pFx + iatom0, -ttx);
          addForce(pFx + iatom1, ttx);
          addForce(pFy + iatom0, -tty);
          addForce(pFy + iatom1, tty);
          addForce(pFz + iatom0, -ttz);
          addForce(pFz + iatom1, ttz);
          if (needVirial) {
            addVirial(pV[0], ttx * dx);
            addVirial(pV[1], tty * dy);
            addVirial(pV[2], ttz * dz);
          }
        }
      }
    }
    pos += blockDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kgCalculateTI14NB_gamd_ppi_kernel:
//
// Arguments:
//   energy:
//   needVirial:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgCalculateTI14NB_ppi_gamd_kernel(bool energy, bool needVirial) {
  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;

  unsigned long long int *pE;
//  TIForce *pFx, *pFy, *pFz;
//  bool useSC14 = false;

  while (pos < cSim.numberTI14NBEntries) {
    unsigned atom0 = cSim.pTINb14ID[pos].x;
    unsigned atom1 = cSim.pTINb14ID[pos].y;
    unsigned TIRegion = cSim.pTINb14ID[pos].z;
    //int TISign = (TIRegion * 2) - 1;
    PMEDouble g1 = cSim.pTINb141[pos].x;
    PMEDouble g2 = cSim.pTINb141[pos].y;
    PMEDouble f6 = cSim.pTINb142[pos].y;
    PMEDouble f12 = cSim.pTINb142[pos].x;
    unsigned iatom0 = cSim.pImageAtomLookup[atom0];
    unsigned iatom1 = cSim.pImageAtomLookup[atom1];
    PMEDouble dx = cSim.pImageX[iatom1] - cSim.pImageX[iatom0];
    PMEDouble dy = cSim.pImageY[iatom1] - cSim.pImageY[iatom0];
    PMEDouble dz = cSim.pImageZ[iatom1] - cSim.pImageZ[iatom0];
    PMEDouble r2 = dx*dx + dy*dy + dz*dz;
    PMEDouble rinv = Rsqrt(r2);
    PMEDouble r2inv = rinv*rinv;
    PMEDouble t6 = ZeroF;


    bool isIntSC = (cSim.pSCList[atom0] > 0 && cSim.pSCList[atom1] > 0);
    //bool hasSC = (cSim.pSCList[atom0] > 0 || cSim.pSCList[atom1] > 0);

    // EED: note that 1:4 electrostatics are double-counted
    PMEDouble de = 0.0;
    PMEDouble ddf = 0.0;
    PMEDouble df = ZeroF;
    PMEDouble df_gamd = ZeroF;

    if (isIntSC) {
      de = g1 * rinv;
      ddf = de * r2inv;
    }
//    pE = (!isIntSC) ? cSim.pTISCEE14[TIRegion] : cSim.pTIEE14SC[TIRegion];
   pE = (isIntSC) ? cSim.pTISCEE14[TIRegion] : cSim.pTIEE14SC[TIRegion];
    addEnergy(pE, de);
    if ((cSim.igamd == 14 ||cSim.igamd == 15||cSim.igamd == 22||cSim.igamd == 24||cSim.igamd == 25||cSim.igamd == 27) && isIntSC){
        df_gamd += ddf;
    }

    // Record GaMDTINBForce
    PMEDouble ttx = df_gamd * dx;
    PMEDouble tty = df_gamd * dy;
    PMEDouble ttz = df_gamd * dz;


    if(cSim.igamd == 14 ||cSim.igamd == 15||cSim.igamd == 24||cSim.igamd == 25){
    addForce(cSim.pGaMDTINBForceX[TIRegion] + iatom0, -ttx);
    addForce(cSim.pGaMDTINBForceX[TIRegion] + iatom1, ttx);
    addForce(cSim.pGaMDTINBForceY[TIRegion] + iatom0, -tty);
    addForce(cSim.pGaMDTINBForceY[TIRegion] + iatom1, tty);
    addForce(cSim.pGaMDTINBForceZ[TIRegion] + iatom0, -ttz);
    addForce(cSim.pGaMDTINBForceZ[TIRegion] + iatom1, ttz);
    } else if (cSim.igamd == 22){
    addForce(cSim.pGaMDTIForceX[TIRegion] + iatom0, -ttx);
    addForce(cSim.pGaMDTIForceX[TIRegion] + iatom1, ttx);
    addForce(cSim.pGaMDTIForceY[TIRegion] + iatom0, -tty);
    addForce(cSim.pGaMDTIForceY[TIRegion] + iatom1, tty);
    addForce(cSim.pGaMDTIForceZ[TIRegion] + iatom0, -ttz);
    addForce(cSim.pGaMDTIForceZ[TIRegion] + iatom1, ttz);
    } //in igamd 22  included the 14 NB in bonded part


    df = OneF; // df must be 1.0 here
    df_gamd = ZeroF;
    // EVDW: note that 14 NB is not double-counted since it is zeroed-out
    t6 = r2inv * r2inv * r2inv;
    f6 *= g2 * t6;
    f12 *= g2 * t6 * t6;
    de = f12 - f6;
    df *= (f12*SixF - f6*ThreeF) * TwoF * r2inv;
//    pE = (!isIntSC) ? cSim.pTISCLJ14[TIRegion] : cSim.pTILJ14[TIRegion];
    pE = (isIntSC) ? cSim.pTISCLJ14[TIRegion] : cSim.pTILJ14[TIRegion];
    addEnergy(pE, de);

    // Record GaMDTINBForce
    if ((cSim.igamd == 14 ||cSim.igamd == 15||cSim.igamd == 22||cSim.igamd == 24||cSim.igamd == 25||cSim.igamd == 27) && isIntSC) {
        df_gamd = df;
    }

    ttx = df_gamd * dx;
    tty = df_gamd * dy;
    ttz = df_gamd * dz;


    if(cSim.igamd == 14||cSim.igamd == 15||cSim.igamd == 24||cSim.igamd == 25){
    addForce(cSim.pGaMDTINBForceX[TIRegion] + iatom0, -ttx);
    addForce(cSim.pGaMDTINBForceX[TIRegion] + iatom1, ttx);
    addForce(cSim.pGaMDTINBForceY[TIRegion] + iatom0, -tty);
    addForce(cSim.pGaMDTINBForceY[TIRegion] + iatom1, tty);
    addForce(cSim.pGaMDTINBForceZ[TIRegion] + iatom0, -ttz);
    addForce(cSim.pGaMDTINBForceZ[TIRegion] + iatom1, ttz);
    }  else if (cSim.igamd == 22){
    addForce(cSim.pGaMDTIForceX[TIRegion] + iatom0, -ttx);
    addForce(cSim.pGaMDTIForceX[TIRegion] + iatom1, ttx);
    addForce(cSim.pGaMDTIForceY[TIRegion] + iatom0, -tty);
    addForce(cSim.pGaMDTIForceY[TIRegion] + iatom1, tty);
    addForce(cSim.pGaMDTIForceZ[TIRegion] + iatom0, -ttz);
    addForce(cSim.pGaMDTIForceZ[TIRegion] + iatom1, ttz);
   } // igamd=22 not include the 14 NB
    pos += gridDim.x * blockDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kgCalculateTINB_ppi_gamd_kernel:
//
// Arguments:
//   needEnergy:
//   needVirial:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgCalculateTINB_ppi_gamd_kernel(bool needEnergy, bool needVirial) {

  unsigned i = threadIdx.x;
  __shared__ PMEFloat recip[9], ucell[9];

  if (cSim.pNTPData == NULL) {
    if (i < 3) {
      recip[i * 3] = cSim.recipf[i][0];
      recip[i * 3 + 1] = cSim.recipf[i][1];
      recip[i * 3 + 2] = cSim.recipf[i][2];
      ucell[i * 3] = cSim.ucellf[i][0];
      ucell[i * 3 + 1] = cSim.ucellf[i][1];
      ucell[i * 3 + 2] = cSim.ucellf[i][2];
    }
  } else {
    if (i < 9) {
      ucell[i] = cSim.pNTPData->ucellf[i];
      recip[i] = cSim.pNTPData->recipf[i];
    }
  }

  __syncthreads();

  PMEFloat4* cn = cSim.pTIcn;
  PMEFloat4* sigEps = cSim.pTISigEps;

  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;

  while (pos < *cSim.pNumberTINBEntries) {

    int4 tt = cSim.pTINBList[pos];

    unsigned iatom0 = cSim.pImageAtomLookup[tt.x];
    unsigned iatom1 = cSim.pImageAtomLookup[tt.y];
    PMEFloat dx = cSim.pImageX[iatom1] - cSim.pImageX[iatom0];
    PMEFloat dy = cSim.pImageY[iatom1] - cSim.pImageY[iatom0];
    PMEFloat dz = cSim.pImageZ[iatom1] - cSim.pImageZ[iatom0];
    PMEFloat r2 = __image_dist2<PMEFloat>(dx, dy, dz, recip, ucell);

    bool excluded = tt.z & 4;

    if (r2 < cSim.cut2 || excluded) {

      bool needLJ = (tt.w >= 0 && (tt.z & 2));
      //bool needSCLJ = tt.z & 8;
      bool needSCEle = tt.z & 16;
      bool isIntSC = tt.z & 32;

      bool hasSC = (cSim.pSCList[tt.x] > 0 || cSim.pSCList[tt.y] > 0);
      //bool addToDVDL_ele = (!hasSC);
      //bool addToDVDL_vdw = (!hasSC);

      //unsigned long long int* pE;
      PMEFloat dE = ZeroF;
      PMEFloat df = ZeroF;
      PMEFloat dfCorr = ZeroF;
      PMEFloat df_gamd = ZeroF;

      int& vdwIndex = tt.w;
      int TIRegion = tt.z & 1;
      //int TISign = ((TIRegion == 0) ? -1 : 1);
      //bool forceSC = (isIntSC && !excluded);
      //bool simpleElec = (!needSCEle || isIntSC || excluded);

      // EED
      PMEFloat q01 = cSim.pOrigAtomCharge[tt.x] * cSim.pOrigAtomCharge[tt.y];
      PMEFloat rinv = Rsqrt(r2);
      PMEFloat r = r2 * rinv;
      PMEFloat swtch = excluded ? -Erf(cSim.ew_coeffSP * r) : Erfc(cSim.ew_coeffSP * r);
      PMEFloat d_swtch = -cSim.negTwoEw_coeffRsqrtPI * Exp(-cSim.ew_coeff2 * r2);
      PMEFloat g1 = q01 * swtch;
      PMEFloat r2inv = rinv*rinv;
      //PMEFloat den = 0.0;
      //PMEFloat den_d = 0.0;

      if(!needSCEle || excluded){
//        dE = g1 *rinv;
//        df = Fma(d_swtch, q01, dE) * r2inv;
      } else if (isIntSC) {
        dE = g1 *rinv;
        dfCorr = Fma(d_swtch, q01, dE) * r2inv;
	if (cSim.igamd0 == 12 || cSim.igamd0 == 13 || cSim.igamd0 == 14 || cSim.igamd0 == 15 ||cSim.igamd0 == 18|| cSim.igamd0 == 19 ||cSim.igamd0 == 21 || cSim.igamd0 == 22||cSim.igamd0 == 23||cSim.igamd0 == 24||cSim.igamd0 == 25||cSim.igamd0 == 110 || cSim.igamd0 ==111|| cSim.igamd0 ==113|| cSim.igamd0 ==115) {
        df_gamd = dfCorr;
        addEnergy(cSim.pTISCEED[TIRegion], dE);
      }
      } else if (hasSC) {
        dE = g1 *rinv;
        df = Fma(d_swtch, q01, dE) * r2inv;
	if (cSim.igamd0 == 12 || cSim.igamd0 == 13 || cSim.igamd0 == 14 || cSim.igamd0 == 15 || cSim.igamd0 == 18 || cSim.igamd0 == 19 || cSim.igamd0 == 21||cSim.igamd0 == 22||cSim.igamd0 == 23||cSim.igamd0 == 24||cSim.igamd0 == 25||cSim.igamd0 == 27||cSim.igamd0 ==111|| cSim.igamd0 == 112|| cSim.igamd0 ==113|| cSim.igamd0 ==115) {
          df_gamd = df;
//          addEnergy(cSim.pTIEESC[TIRegion], dE);
	  addEnergy(cSim.pTISCEED[TIRegion], dE);
	}
      }
//recorded the force to GaMDTINBFoces
      PMEFloat    ttx = df_gamd * dx;
      PMEFloat    tty = df_gamd * dy;
      PMEFloat    ttz = df_gamd * dz;
      if(cSim.igamd0 ==25){
      addForce(cSim.pGaMDTINBForceX[TIRegion+1] + iatom0, -ttx);
      addForce(cSim.pGaMDTINBForceX[TIRegion+1] + iatom1, ttx);
      addForce(cSim.pGaMDTINBForceY[TIRegion+1] + iatom0, -tty);
      addForce(cSim.pGaMDTINBForceY[TIRegion+1] + iatom1, tty);
      addForce(cSim.pGaMDTINBForceZ[TIRegion+1] + iatom0, -ttz);
      addForce(cSim.pGaMDTINBForceZ[TIRegion+1] + iatom1, ttz);
      }else{
      addForce(cSim.pGaMDTINBForceX[TIRegion] + iatom0, -ttx);
      addForce(cSim.pGaMDTINBForceX[TIRegion] + iatom1, ttx);
      addForce(cSim.pGaMDTINBForceY[TIRegion] + iatom0, -tty);
      addForce(cSim.pGaMDTINBForceY[TIRegion] + iatom1, tty);
      addForce(cSim.pGaMDTINBForceZ[TIRegion] + iatom0, -ttz);
      addForce(cSim.pGaMDTINBForceZ[TIRegion] + iatom1, ttz);
      }
      df=ZeroF;
      dfCorr=ZeroF;
      df_gamd=ZeroF;
      //forceSC = (isIntSC && !excluded);
      // EVDW
      if (needLJ && !excluded) {
        PMEFloat t4, f4, t6, f6, f12, dd;
        //pE = (!isIntSC) ? cSim.pTILJ[TIRegion] : cSim.pTISCLJ[TIRegion];
        if (isIntSC) {
          t4 = r2inv * r2inv;
          f4 = (sigEps[vdwIndex].w > 0) ? cn[vdwIndex].z * t4 : ZeroF;
          t6 = t4 * r2inv;
          f6 = cn[vdwIndex].y * t6;
          f12 = cn[vdwIndex].x * t6 * t6;
          if (vdwIndex == 55555) {
            for (int i=0;i<cSim.TIC4Pairwise;i++) {
              if ((cSim.pImageAtom[tt.x] == (cSim.pTIDcoef[i*2]) && cSim.pImageAtom[tt.y] == (cSim.pTIDcoef[i*2+1])) || (cSim.pImageAtom[tt.x] == (cSim.pTIDcoef[i*2+1]) && cSim.pImageAtom[tt.y] == (cSim.pTIDcoef[i*2]))) {
                f4 += cSim.pTIDvalue[i] * t4 / 3.0; // C4PairwiseCUDA2023
                //eLJ -= f4;
                //dvdu += (TwoF * t4) / (ThreeF * u1_3);
              }
            }
          }
          dd = (f12*TwelveF - f6*SixF - f4*FourF) * r2inv;
	  if(cSim.igamd0 == 12 || cSim.igamd0 == 13 || cSim.igamd0 == 14 || cSim.igamd0 == 15 || cSim.igamd0 == 18 || cSim.igamd0 == 19 || cSim.igamd0 == 21 ||cSim.igamd0 == 22||cSim.igamd0 == 23 ||cSim.igamd0 == 24||cSim.igamd0 == 25||cSim.igamd0 == 27|| cSim.igamd0 ==111|| cSim.igamd0 ==113||cSim.igamd0 == 114 || cSim.igamd0 == 115){
          addEnergy(cSim.pTISCLJ[TIRegion], (f12 - f6 - f4));
          df_gamd += dd;
         }
         } else if (hasSC) {
             t4 = r2inv * r2inv;
             f4 = (sigEps[vdwIndex].w > 0) ? cn[vdwIndex].z * t4 : ZeroF;
             t6 = t4 * r2inv;
             f6 = cn[vdwIndex].y * t6;
             f12 = cn[vdwIndex].x * t6 * t6;
             if (vdwIndex == 55555) {
               for (int i=0;i<cSim.TIC4Pairwise;i++) {
                 if ((cSim.pImageAtom[tt.x] == (cSim.pTIDcoef[i*2]) && cSim.pImageAtom[tt.y] == (cSim.pTIDcoef[i*2+1])) || (cSim.pImageAtom[tt.x] == (cSim.pTIDcoef[i*2+1]) && cSim.pImageAtom[tt.y] == (cSim.pTIDcoef[i*2]))) {
                   f4 += cSim.pTIDvalue[i] * t4 / 3.0; // C4PairwiseCUDA2023
                   //eLJ -= f4;
                   //dvdu += (TwoF * t4) / (ThreeF * u1_3);
                 }
               }
             }
             dd = (f12*TwelveF - f6*SixF - f4*FourF) * r2inv;
	     if (cSim.igamd0 == 12 || cSim.igamd0 == 13 || cSim.igamd0 == 14 || cSim.igamd0 == 15 || cSim.igamd0 == 18||cSim.igamd0 == 19 || cSim.igamd0 == 21 || cSim.igamd0 == 22||cSim.igamd0 == 23 ||cSim.igamd0 == 24||cSim.igamd0 == 25||cSim.igamd0 == 27||cSim.igamd0 ==111|| cSim.igamd0 ==113|| cSim.igamd0 ==115||cSim.igamd0 == 116) {
                df_gamd += dd;
		addEnergy(cSim.pTISCLJ[TIRegion], (f12 - f6 - f4));
              }
             }
        }
   /*   else {
	  PMEFloat t4, f4, t6, f6, f12, dd, dl, u1_3;
         if(isIntSC){
             t4 = r2inv * r2inv;
             f4 = (sigEps[vdwIndex].w > 0) ? cn[vdwIndex].z * t4 : ZeroF;
             t6 = t4 * r2inv;
             f6 = cn[vdwIndex].y * t6;
             f12 = cn[vdwIndex].x * t6 * t6;
             dd = (f12*TwelveF - f6*SixF - f4*FourF) * r2inv;
             addEnergy(cSim.pTISCLJ[TIRegion], (f12 - f6 - f4));
	     if(cSim.igamd == 12 || cSim.igamd == 13 || cSim.igamd == 14 || cSim.igamd == 15 || cSim.igamd == 18 || cSim.igamd == 19 ||cSim.igamd == 114 || cSim.igamd == 115) {
                df_gamd += dd;
                }
            } else if (hasSC) {
             t4 = r2inv * r2inv;
             f4 = (sigEps[vdwIndex].w > 0) ? cn[vdwIndex].z * t4 : ZeroF;
             t6 = t4 * r2inv;
             f6 = cn[vdwIndex].y * t6;
             f12 = cn[vdwIndex].x * t6 * t6;
             dd = (f12*TwelveF - f6*SixF - f4*FourF) * r2inv;
             addEnergy(cSim.pTILJ[TIRegion], (f12 - f6 - f4));
	     if(cSim.igamd == 12 || cSim.igamd == 13 || cSim.igamd == 14 || cSim.igamd == 15 || cSim.igamd == 18 || cSim.igamd == 19 ||cSim.igamd == 114 || cSim.igamd == 115) {
             df_gamd += dd;
            }
	    }
	} */
      // record GaMDTINBForce
      ttx = df_gamd * dx;
      tty = df_gamd * dy;
      ttz = df_gamd * dz;
      if(cSim.igamd0 == 23 ||cSim.igamd0 == 24){
      addForce(cSim.pGaMDTINBForceX[TIRegion+1] + iatom0, -ttx);
      addForce(cSim.pGaMDTINBForceX[TIRegion+1] + iatom1, ttx);
      addForce(cSim.pGaMDTINBForceY[TIRegion+1] + iatom0, -tty);
      addForce(cSim.pGaMDTINBForceY[TIRegion+1] + iatom1, tty);
      addForce(cSim.pGaMDTINBForceZ[TIRegion+1] + iatom0, -ttz);
      addForce(cSim.pGaMDTINBForceZ[TIRegion+1] + iatom1, ttz);
      }else{
      addForce(cSim.pGaMDTINBForceX[TIRegion] + iatom0, -ttx);
      addForce(cSim.pGaMDTINBForceX[TIRegion] + iatom1, ttx);
      addForce(cSim.pGaMDTINBForceY[TIRegion] + iatom0, -tty);
      addForce(cSim.pGaMDTINBForceY[TIRegion] + iatom1, tty);
      addForce(cSim.pGaMDTINBForceZ[TIRegion] + iatom0, -ttz);
      addForce(cSim.pGaMDTINBForceZ[TIRegion] + iatom1, ttz);
      }
    }
    pos += gridDim.x * blockDim.x;
  }
}
//---------------------------------------------------------------------------------------------
// kgCalculateTINB_ppi2_gamd_kernel:
//
// Arguments: only for igamd.eq.16, 17,20,26 and 28
//   needEnergy:
//   needVirial:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgCalculateTINB_ppi2_gamd_kernel(bool needEnergy, bool needVirial, int bgpro2atm, int edpro2atm) {

  unsigned i = threadIdx.x;
  __shared__ PMEFloat recip[9], ucell[9];

  if (cSim.pNTPData == NULL) {
    if (i < 3) {
      recip[i * 3] = cSim.recipf[i][0];
      recip[i * 3 + 1] = cSim.recipf[i][1];
      recip[i * 3 + 2] = cSim.recipf[i][2];
      ucell[i * 3] = cSim.ucellf[i][0];
      ucell[i * 3 + 1] = cSim.ucellf[i][1];
      ucell[i * 3 + 2] = cSim.ucellf[i][2];
    }
  } else {
    if (i < 9) {
      ucell[i] = cSim.pNTPData->ucellf[i];
      recip[i] = cSim.pNTPData->recipf[i];
    }
  }

  __syncthreads();

  PMEFloat4* cn = cSim.pTIcn;
  PMEFloat4* sigEps = cSim.pTISigEps;

  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;

  while (pos < *cSim.pNumberTINBEntries) {

    int4 tt = cSim.pTINBList[pos];

    unsigned iatom0 = cSim.pImageAtomLookup[tt.x];
    unsigned iatom1 = cSim.pImageAtomLookup[tt.y];
    PMEFloat dx = cSim.pImageX[iatom1] - cSim.pImageX[iatom0];
    PMEFloat dy = cSim.pImageY[iatom1] - cSim.pImageY[iatom0];
    PMEFloat dz = cSim.pImageZ[iatom1] - cSim.pImageZ[iatom0];
    PMEFloat r2 = __image_dist2<PMEFloat>(dx, dy, dz, recip, ucell);

    bool excluded = tt.z & 4;

    if (r2 < cSim.cut2 || excluded) {

      bool needLJ = (tt.w >= 0 && (tt.z & 2));
      //bool needSCLJ = tt.z & 8;
      bool needSCEle = tt.z & 16;
      bool isIntSC = tt.z & 32;

      bool hasSC = (cSim.pSCList[tt.x] > 0 || cSim.pSCList[tt.y] > 0);
      bool is2region=((tt.x>=(bgpro2atm-1) && tt.x<edpro2atm) || (tt.y>=(bgpro2atm-1) && tt.y<edpro2atm));

      //unsigned long long int* pE;
      PMEFloat dE = ZeroF;
      PMEFloat df = ZeroF;
      //PMEFloat dfCorr = ZeroF;
      PMEFloat df_gamd = ZeroF;

      int& vdwIndex = tt.w;
      int TIRegion = tt.z & 1;
      //int TISign = ((TIRegion == 0) ? -1 : 1);
      //bool forceSC = (isIntSC && !excluded);
      //bool simpleElec = (!needSCEle || isIntSC || excluded);

      // EED
      PMEFloat q01 = cSim.pOrigAtomCharge[tt.x] * cSim.pOrigAtomCharge[tt.y];
      PMEFloat rinv = Rsqrt(r2);
      PMEFloat r = r2 * rinv;
      PMEFloat swtch = excluded ? -Erf(cSim.ew_coeffSP * r) : Erfc(cSim.ew_coeffSP * r);
      PMEFloat d_swtch = -cSim.negTwoEw_coeffRsqrtPI * Exp(-cSim.ew_coeff2 * r2);
      PMEFloat g1 = q01 * swtch;
      PMEFloat r2inv = rinv*rinv;
      //PMEFloat den = 0.0;
      //PMEFloat den_d = 0.0;

      if(!needSCEle || excluded){
//        dE = g1 *rinv;
//        df = Fma(d_swtch, q01, dE) * r2inv;
      } else if (isIntSC) {
         dE = g1 *rinv;
         df = Fma(d_swtch, q01, dE) * r2inv;
         if (cSim.igamd0 == 28){
	     df_gamd = df;
        addEnergy(cSim.pTISCEED[TIRegion], dE);
		}
      } else if (hasSC) {
//        bool is2region=((iatom0>=(bgpro2atm-1) && iatom0<edpro2atm) || (iatom1>=(bgpro2atm-1) && iatom1<edpro2atm));
        if(is2region){
         dE = g1 *rinv;
         df = Fma(d_swtch, q01, dE) * r2inv;
         if (cSim.igamd0 == 16 || cSim.igamd0 == 17||cSim.igamd0 == 20||cSim.igamd0 == 26||cSim.igamd0 == 28){
	     df_gamd = df;
        addEnergy(cSim.pTISCEED[TIRegion], dE);
	 }
	}
      }
//recorded the force to GaMDTINBFoces
      PMEFloat    ttx = df_gamd * dx;
      PMEFloat    tty = df_gamd * dy;
      PMEFloat    ttz = df_gamd * dz;
      addForce(cSim.pGaMDTINBForceX[TIRegion] + iatom0, -ttx);
      addForce(cSim.pGaMDTINBForceX[TIRegion] + iatom1, ttx);
      addForce(cSim.pGaMDTINBForceY[TIRegion] + iatom0, -tty);
      addForce(cSim.pGaMDTINBForceY[TIRegion] + iatom1, tty);
      addForce(cSim.pGaMDTINBForceZ[TIRegion] + iatom0, -ttz);
      addForce(cSim.pGaMDTINBForceZ[TIRegion] + iatom1, ttz);

      df=ZeroF;
      //dfCorr=ZeroF;
      df_gamd=ZeroF;
      // EVDW
      if (needLJ && !excluded) {
        PMEFloat t4, f4, t6, f6, f12, dd;
        //pE = (!isIntSC) ? cSim.pTILJ[TIRegion] : cSim.pTISCLJ[TIRegion];
		if(isIntSC){
		  t4 = r2inv * r2inv;
          f4 = (sigEps[vdwIndex].w > 0) ? cn[vdwIndex].z * t4 : ZeroF;
          t6 = t4 * r2inv;
          f6 = cn[vdwIndex].y * t6;
          f12 = cn[vdwIndex].x * t6 * t6;
          if (vdwIndex == 55555) {
            f4 = 0;
            for (int i=0;i<cSim.TIC4Pairwise;i++) {
              if ((cSim.pImageAtom[tt.x] == (cSim.pTIDcoef[i*2]) && cSim.pImageAtom[tt.y] == (cSim.pTIDcoef[i*2+1])) || (cSim.pImageAtom[tt.x] == (cSim.pTIDcoef[i*2+1]) && cSim.pImageAtom[tt.y] == (cSim.pTIDcoef[i*2]))) {
                f4 += cSim.pTIDvalue[i] * t4 / 3.0; // C4PairwiseCUDA2023
                //eLJ -= f4;
                //dvdu += (TwoF * t4) / (ThreeF * u1_3);
              }
            }
          }
          dd = (f12*TwelveF - f6*SixF - f4*FourF) * r2inv;
          if(cSim.igamd0 == 28) {
          df_gamd += dd;
	      addEnergy(cSim.pTISCLJ[TIRegion], (f12 - f6 - f4));
         }
		}
        if (hasSC && is2region) {
          t4 = r2inv * r2inv;
          f4 = (sigEps[vdwIndex].w > 0) ? cn[vdwIndex].z * t4 : ZeroF;
          t6 = t4 * r2inv;
          f6 = cn[vdwIndex].y * t6;
          f12 = cn[vdwIndex].x * t6 * t6;
          if (vdwIndex == 55555) {
            for (int i=0;i<cSim.TIC4Pairwise;i++) {
              if ((cSim.pImageAtom[tt.x] == (cSim.pTIDcoef[i*2]) && cSim.pImageAtom[tt.y] == (cSim.pTIDcoef[i*2+1])) || (cSim.pImageAtom[tt.x] == (cSim.pTIDcoef[i*2+1]) && cSim.pImageAtom[tt.y] == (cSim.pTIDcoef[i*2]))) {
                f4 += cSim.pTIDvalue[i] * t4 / 3.0; // C4PairwiseCUDA2023
                //eLJ -= f4;
                //dvdu += (TwoF * t4) / (ThreeF * u1_3);
              }
            }
          }
          dd = (f12*TwelveF - f6*SixF - f4*FourF) * r2inv;
          if(cSim.igamd0 == 16 || cSim.igamd0 == 17||cSim.igamd0 == 20||cSim.igamd0 == 26||cSim.igamd0 == 28) {
          df_gamd += dd;
	      addEnergy(cSim.pTISCLJ[TIRegion], (f12 - f6 - f4));
         }
	}
        }

      // record GaMDTINBForce
      ttx = df_gamd * dx;
      tty = df_gamd * dy;
      ttz = df_gamd * dz;
      if(cSim.igamd0 == 16 || cSim.igamd0 == 17||cSim.igamd0 == 26||cSim.igamd0 == 28){
      addForce(cSim.pGaMDTINBForceX[TIRegion] + iatom0, -ttx);
      addForce(cSim.pGaMDTINBForceX[TIRegion] + iatom1, ttx);
      addForce(cSim.pGaMDTINBForceY[TIRegion] + iatom0, -tty);
      addForce(cSim.pGaMDTINBForceY[TIRegion] + iatom1, tty);
      addForce(cSim.pGaMDTINBForceZ[TIRegion] + iatom0, -ttz);
      addForce(cSim.pGaMDTINBForceZ[TIRegion] + iatom1, ttz);
      }else if (cSim.igamd0 == 20) {
      addForce(cSim.pGaMDTINBForceX[TIRegion+1] + iatom0, -ttx);
      addForce(cSim.pGaMDTINBForceX[TIRegion+1] + iatom1, ttx);
      addForce(cSim.pGaMDTINBForceY[TIRegion+1] + iatom0, -tty);
      addForce(cSim.pGaMDTINBForceY[TIRegion+1] + iatom1, tty);
      addForce(cSim.pGaMDTINBForceZ[TIRegion+1] + iatom0, -ttz);
      addForce(cSim.pGaMDTINBForceZ[TIRegion+1] + iatom1, ttz);
       }
    }
    pos += gridDim.x * blockDim.x;
  }
}


#endif
