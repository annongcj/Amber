#ifdef GTI

#ifndef AMBER_PLATFORM_AMD
#include <cuda.h>
#endif
#include "gputypes.h"
#include "ptxmacros.h"
#include "gti_const.cuh"
#include "gti_general_kernels.cuh"
#include "gti_utils.cuh"

#include "simulationConst.h"
CSIM_STO simulationConst cSim;

namespace GTI_GENERAL_IMPL {
#  if !(defined(__CUDACC_RDC__) || defined(__HIPCC_RDC__))
#    include "gti_localCUutils.inc"
#  endif
}
using namespace GTI_GENERAL_IMPL;

template<class T>
__forceinline__ __device__  void vec_sync(T* pVector, unsigned& a0, unsigned& a1,
  int combinedMode) {

  if (combinedMode <= 2) {
    if (combinedMode == 0) {
      pVector[a0] += pVector[a1];
      pVector[a0 + cSim.stride] += pVector[a1 + cSim.stride];
      pVector[a0 + cSim.stride2] += pVector[a1 + cSim.stride2];
    } else if (combinedMode == 1) {
      pVector[a0] = pVector[a0] * cSim.TIWeight[0]
        + pVector[a1] * cSim.TIWeight[1];
      pVector[a0 + cSim.stride] = pVector[a0 + cSim.stride] * cSim.TIWeight[0]
        + pVector[a1 + cSim.stride] * cSim.TIWeight[1];
      pVector[a0 + cSim.stride2] = pVector[a0 + cSim.stride2] * cSim.TIWeight[0]
        + pVector[a1 + cSim.stride2] * cSim.TIWeight[1];
    }
    pVector[a1] = pVector[a0];
    pVector[a1 + cSim.stride] = pVector[a0 + cSim.stride];
    pVector[a1 + cSim.stride2] = pVector[a0 + cSim.stride2];
  } else {
    pVector[a0] = pVector[a1];
    pVector[a0+ cSim.stride] = pVector[a1 + cSim.stride];
    pVector[a0+ cSim.stride2] = pVector[a1 + cSim.stride2];
  }
}

//---------------------------------------------------------------------------------------------
// kgSyncVector_kernel:
//
// Arguments:
//   useImage:
//   mode: 0: force; 1: velocity; 2: coord; 3:coord+velocity
//   combinedMode:  0: direct combination; 1: weighted combination; 2: copy V0 to V1;  3: copy V1 to V0 ; 4: SHAKE mode
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgSyncVector_kernel(bool useImage, unsigned mode, int combinedMode)
{
  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  unsigned numberPairs = cSim.numberTICommonPairs;

  if (numberPairs == 0) {
    return;
  }
  if (combinedMode < 0) {
    combinedMode = mode;
  }
  while (pos < numberPairs) {
    unsigned atom0 = cSim.pTICommonPair[pos].x;
    unsigned atom1 = cSim.pTICommonPair[pos].y;
    unsigned a0 = (useImage) ? cSim.pImageAtomLookup[atom0] : atom0;
    unsigned a1 = (useImage) ? cSim.pImageAtomLookup[atom1] : atom1;

    int mycombinedMode = (combinedMode != 4) ? combinedMode :
      ( (cSim.pTICommonPair[pos].z == 1) ? 3 : 2);

    switch (mode) {
    case(0): vec_sync<PMEAccumulator>(cSim.pForceAccumulator, a0, a1, mycombinedMode);
      if(cSim.pNBForceAccumulator!= cSim.pForceAccumulator)
        vec_sync<PMEAccumulator>(cSim.pNBForceAccumulator, a0, a1, mycombinedMode);
      break;
    case(1): vec_sync<double>((useImage) ? cSim.pImageVelX :
      cSim.pVelX, a0, a1, mycombinedMode);
      break;
    case(2): vec_sync<double>((useImage) ? cSim.pImageX :
      cSim.pAtomX, a0, a1, mycombinedMode);
      break;
    case(3): vec_sync<double>((useImage) ? cSim.pImageX :
      cSim.pAtomX, a0, a1, mycombinedMode);
      vec_sync<double>((useImage) ? cSim.pImageVelX :
        cSim.pVelX, a0, a1, mycombinedMode);
      break;
    }
    pos += increment;
  }
}


//---------------------------------------------------------------------------------------------
// kgClearTIPotEnergy_kernel:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgClearTIPotEnergy_kernel()
{
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  unsigned long long int *pTIE = cSim.pTIPotEnergyBuffer;
  while (pos < cSim.GPUPotEnergyTerms * cSim.TIEnergyBufferMultiplier) {
    pTIE[pos] = Zero;
    pos += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kgClearTIKinEnergy_kernel:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgClearTIKinEnergy_kernel()
{
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  unsigned long long int *pTIE = cSim.pTIKinEnergyBuffer;
  while (pos < cSim.GPUKinEnergyTerms * 3) {
    pTIE[pos] = Zero;
    pos += increment;
  }

}

//---------------------------------------------------------------------------------------------
// kgClearMBAR_kernel:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgClearMBAR_kernel()
{
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  unsigned long long int *pE = cSim.pMBAREnergy;
  while (pos < cSim.nMBARStates* Schedule::TypeTotal *3) {
    pE[pos] = Zero;
    pos += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kgClearTIForce_kernel:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgClearTIForce_kernel()
{
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  pos = blockIdx.x * blockDim.x + threadIdx.x;

  PMEAccumulator* pTIF = cSim.pTIForce;
  unsigned multiplier = (cSim.needVirial) ? 10 : 5;
  while (pos < cSim.stride3 * multiplier) {
    pTIF[pos] = Zero;
    pos += increment;
  }

  pos = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned long long int *pTIE = cSim.pTIPotEnergyBuffer;
  while (pos < cSim.GPUPotEnergyTerms * cSim.TIEnergyBufferMultiplier) {
    pTIE[pos] = Zero;
    pos += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kgClearTIForce_gamd_kernel:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgClearTIForce_gamd_kernel()
{
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  pos = blockIdx.x * blockDim.x + threadIdx.x;

  PMEAccumulator* pTIF = cSim.pTIForce;
  unsigned multiplier = (cSim.needVirial) ? 14 : 7;
  while (pos < cSim.stride3 * multiplier) {
    pTIF[pos] = ZeroF;
    pos += increment;
  }

  pos = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned long long int *pTIE = cSim.pTIPotEnergyBuffer;
  while (pos < cSim.GPUPotEnergyTerms * 3) {
    pTIE[pos] = ZeroF;
    pos += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kgZeroTICharges_kernel:
//
// Arguments:
//
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgZeroTICharges_kernel(bool useImage, int mode)
{
  // mode=0: zero out TI region 1
  // mode=1: zero out TI region 0
  // mode=2: zero out both TI regions
  // mode=-1: restore to the original values; TBD

  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  PMEFloat charge, charge0=ZeroF;

  while (pos < cSim.numberTIAtoms) {
    unsigned int atom0 = cSim.pTIAtomList[pos].x;
    unsigned int myRegion = cSim.pTIAtomList[pos].y;
    unsigned int atomIndex = ( useImage ) ? cSim.pImageAtomLookup[atom0] : atom0;

    charge = cSim.pOrigAtomCharge[atom0];
    switch (mode) {
      case 0:
        if (myRegion==1) charge = charge0;  break;
      case 1:
        if (myRegion==0) charge = charge0;  break;
      case 2:
        charge = charge0;  break;
    }

    cSim.pAtomCharge[atomIndex] = charge;
    cSim.pAtomChargeSP[atomIndex] = charge;
    cSim.pAtomChargeSPLJID[atomIndex].x = charge;
    cSim.pImageCharge[atomIndex] = charge;

    pos += increment;
  }
}


//---------------------------------------------------------------------------------------------
// kgScaleRECharges_kernel:
//
// Arguments:
//
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgScaleRECharges_kernel(bool useImage, int mode) {
  // mode=0: scale RE charges of READ regtion
  // mode=2: restore charges

  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  PMEFloat charge, charge0 = ZeroF;
  PMEFloat scale = cSim.REAFItemWeight[Schedule::TypeEleRec][0];

  while (pos < cSim.numberREAFAtoms) {
    unsigned int atom0 = cSim.pREAFAtomList[pos];
    unsigned int atomIndex = (useImage) ? cSim.pImageAtomLookup[atom0] : atom0;

    charge = cSim.pOrigAtomCharge[atom0];
    switch (mode) {
    case 0: default: charge*= scale;  break;
    case 2:
      charge = charge0;
      break;
    }

    cSim.pAtomCharge[atomIndex] = charge;
    cSim.pAtomChargeSP[atomIndex] = charge;
    cSim.pAtomChargeSPLJID[atomIndex].x = charge;
    cSim.pImageCharge[atomIndex] = charge;

    pos += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kgZeroTIAtomForce_kernel:
//
// Arguments:
//   useImage:
//   mode:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgZeroTIAtomForce_kernel(bool useImage, int mode)
{
  // mode=0: zero out H0 region
  // mode=1: zero out H1 region
  // mode=2: zero out both H0 and H1 regions
  // mode=3: zero out both H0 and H1 regions AND main forces of TI atoms

  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  PMEAccumulator* pTarget = cSim.pTIForce;
  PMEAccumulator * pTargetMain = cSim.pForceAccumulator;

  while (pos < cSim.numberTIAtoms) {
    unsigned int realAtomindex = cSim.pTIAtomList[pos].x;
    unsigned int TIRegion = cSim.pTIAtomList[pos].y;
    unsigned int atomIndex = (useImage) ? cSim.pImageAtomLookup[realAtomindex] : realAtomindex;

    // normal TI region
    PMEAccumulator* pX = pTarget + TIRegion*cSim.stride3;
    PMEAccumulator* pY = pX + cSim.stride;
    PMEAccumulator* pZ = pX + cSim.stride * 2;
    pX[atomIndex] = 0; pY[atomIndex] = 0; pZ[atomIndex] = 0;

    // SC TI region
    PMEAccumulator* pSX = pTarget + (TIRegion+3)*cSim.stride3;
    PMEAccumulator* pSY = pSX + cSim.stride;
    PMEAccumulator* pSZ = pSX + cSim.stride * 2;
    pSX[atomIndex] = 0; pSY[atomIndex] = 0; pSZ[atomIndex] = 0;
    if (mode >= 2) {

      // normal TI region
      PMEAccumulator* pX1 = pTarget + ( 1 - TIRegion )*cSim.stride3;
      PMEAccumulator* pY1 = pX1 + cSim.stride;
      PMEAccumulator* pZ1 = pX1 + cSim.stride * 2;

      // SC TIRegion
      PMEAccumulator* pSX1 = pTarget + ( 4 - TIRegion )*cSim.stride3;
      PMEAccumulator* pSY1 = pSX1 + cSim.stride;
      PMEAccumulator* pSZ1 = pSX1 + cSim.stride * 2;

      pX1[atomIndex] = 0; pY1[atomIndex] = 0; pZ1[atomIndex] = 0;
      pSX1[atomIndex] = 0; pSY1[atomIndex] = 0; pSZ1[atomIndex] = 0;
    }
    if (mode == 3) {
      pTargetMain[atomIndex] = 0;
      (pTargetMain + cSim.stride)[atomIndex] = 0;
      (pTargetMain + (cSim.stride * 2))[atomIndex] = 0;
    }

    pos += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kgCopyToTIForce:
//
// Arguments:
//   TIRegion:
//   isNB:
//   keepSource:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgCopyToTIForce(int TIRegion, bool isNB, bool keepSource, PMEFloat weight)
{
  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  PMEAccumulator *pMain = (isNB) ? cSim.pNBForceAccumulator : cSim.pForceAccumulator;
  PMEAccumulator* pTI = (isNB) ? cSim.pTINBForce : cSim.pTIForce;
  pTI += TIRegion*cSim.stride3;

  if (TIRegion >= 0) {
    while (pos < cSim.stride3) {
      PMEFloat tt = ((PMEFloat)pMain[pos]) * weight * ONEOVERFORCESCALE;
      atomicAdd((unsigned long long*)&(pTI[pos]), (PMEAccumulator)(tt * FORCESCALE));
      if (!keepSource) {
        pMain[pos] = ZeroF;
      }
      pos += increment;
    }
  }
  else {
    while (pos < cSim.stride3) {
      pMain[pos] = (keepSource) ? (-pMain[pos]) : ZeroF;
      pos += increment;
    }
  }
}

//---------------------------------------------------------------------------------------------
// kgCombineTIForce_kernel:
//
// Arguments:
//   isNB:
//   useImage:
//   addOn:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgCombineTIForce_kernel(bool isNB, bool useImage, bool addOn, bool needWeight)
{
  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  PMEAccumulator *pMain = (isNB) ? cSim.pNBForceAccumulator : cSim.pForceAccumulator;
  PMEAccumulator* pTI0 = (isNB) ? cSim.pTINBForce : cSim.pTIForce;

  PMEAccumulator* pTI1 = pTI0 + cSim.stride3;
  PMEAccumulator* pTIT = pTI0 + cSim.stride3 * 2;

  double weight0 = cSim.TIWeight[0];
  double weight1 = cSim.TIWeight[1];

  PMEDouble4 tt;
  while (pos < cSim.stride3) {
    if (!addOn) pMain[pos] = 0;
    if (needWeight) {
      tt.x = ((PMEDouble)pTI0[pos]) * ONEOVERFORCESCALE;
      tt.y = ((PMEDouble)pTI1[pos]) * ONEOVERFORCESCALE;
      tt.z = ((PMEDouble)pTIT[pos]) * ONEOVERFORCESCALE;
      tt.w = (tt.x * weight0) + (tt.y * weight1) + tt.z;
      pMain[pos] += (PMEAccumulator)(tt.w * FORCESCALE);
    } else {
      pMain[pos] += pTI0[pos] + pTI1[pos] + pTIT[pos];
    }
    pos += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kgCombineSCForce_kernel:
//
// Arguments:
//   useImage:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgCombineSCForce_kernel(bool isNB, bool useImage)
{
  PMEAccumulator *pMain = (isNB) ? cSim.pNBForceAccumulator : cSim.pForceAccumulator;
  PMEAccumulator* pTI0 = (isNB) ? cSim.pTINBForce : cSim.pTIForce;

  PMEAccumulator* pSC0 = pTI0 + cSim.stride3 * 3;
  PMEAccumulator* pSC1 = pTI0 + cSim.stride3 * 4;

  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  while (pos < cSim.atoms) {
    //unsigned int realAtomindex = pos;
    //unsigned int atomindex = (useImage) ? cSim.pImageAtomLookup[realAtomindex] : realAtomindex;
    unsigned int atomindex = pos;
    unsigned int realAtomindex = (useImage) ? cSim.pImageAtom[atomindex] : atomindex;
    bool onlySC0 = (cSim.pTIList[realAtomindex] > 0 && cSim.pSCList[realAtomindex] > 0);
    bool onlySC1 = (cSim.pTIList[realAtomindex + cSim.stride] > 0 &&
                    cSim.pSCList[realAtomindex] > 0);
    if (!onlySC0) {
      pMain[atomindex] += pSC1[atomindex];
      pMain[atomindex + cSim.stride] += pSC1[atomindex + cSim.stride] ;
      pMain[atomindex + cSim.stride * 2] += pSC1[atomindex + cSim.stride * 2] ;
    }
    if (!onlySC1) {
      pMain[atomindex] += pSC0[atomindex] ;
      pMain[atomindex + cSim.stride] += pSC0[atomindex + cSim.stride] ;
      pMain[atomindex + cSim.stride * 2] += pSC0[atomindex + cSim.stride * 2] ;
    }

    pos += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kgCombineTIForce_gamd_kernel:
//
// Arguments:
//   isNB:
//   useImage:
//   addOn:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgCombineTIForce_gamd_kernel(bool isNB, bool useImage, bool addOn)
{
  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  PMEAccumulator *pMain = (isNB) ? cSim.pNBForceAccumulator : cSim.pForceAccumulator;
  PMEAccumulator*pTI0 = (isNB) ? cSim.pTINBForce : cSim.pTIForce;

  PMEAccumulator* pTI1 = pTI0 + cSim.stride3;
  PMEAccumulator* pTIT = pTI0 + cSim.stride3 * 2;

  double weight0 = cSim.TIWeight[0];
  double weight1 = cSim.TIWeight[1];

  PMEDouble4 tt;
  while (pos < cSim.stride3) {
    tt.x = ((PMEDouble)pTI0[pos]) * ONEOVERFORCESCALE;
    tt.y = ((PMEDouble)pTI1[pos]) * ONEOVERFORCESCALE;
    tt.z = ((PMEDouble)pTIT[pos]) * ONEOVERFORCESCALE;
    tt.w = (tt.x * weight0) + (tt.y * weight1) + tt.z;
    if (!addOn) {
      pMain[pos] = 0;
    }
    pMain[pos] += (PMEAccumulator)(tt.w * FORCESCALE);
    pos += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kgCombineSCForce_gamd_kernel:
//
// Arguments:
//   useImage:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgCombineSCForce_gamd_kernel(bool isNB, bool useImage)
{
  PMEAccumulator *pMain = (isNB) ? cSim.pNBForceAccumulator : cSim.pForceAccumulator;
  PMEAccumulator*pTI0 = (isNB) ? cSim.pTINBForce : cSim.pTIForce;

  PMEAccumulator* pSC0 = pTI0 + cSim.stride3 * 3;
  PMEAccumulator* pSC1 = pTI0 + cSim.stride3 * 4;

  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  while (pos < cSim.atoms) {
    unsigned int realAtomindex = pos;
    unsigned int atomindex = (useImage) ? cSim.pImageAtomLookup[realAtomindex] : realAtomindex;
    bool onlySC0 = (cSim.pTIList[realAtomindex] > 0 && cSim.pSCList[realAtomindex] > 0);
    bool onlySC1 = (cSim.pTIList[realAtomindex + cSim.stride] > 0 &&
                    cSim.pSCList[realAtomindex] > 0);
    if (!onlySC0) {
      pMain[atomindex] += pSC1[atomindex];
      pMain[atomindex + cSim.stride] += pSC1[atomindex + cSim.stride] ;
      pMain[atomindex + cSim.stride * 2] += pSC1[atomindex + cSim.stride * 2] ;
    }
    if (!onlySC1) {
      pMain[atomindex] += pSC0[atomindex] ;
      pMain[atomindex + cSim.stride] += pSC0[atomindex + cSim.stride] ;
      pMain[atomindex + cSim.stride * 2] += pSC0[atomindex + cSim.stride * 2] ;
    }

    pos += increment;
  }
}

//---------------------------------------------------------------------------------------------
// NetForce:
//---------------------------------------------------------------------------------------------
namespace NetForce {
  __device__ unsigned long long int netFrcX = 0;
  __device__ unsigned long long int netFrcY = 0;
  __device__ unsigned long long int netFrcZ = 0;
  __device__ unsigned ignored = 0;
  __device__ unsigned int blockCount = 0;
#ifdef use_DPFP
  __constant__ PMEFloat small = 1e-8;
  __constant__ PMEFloat forceCut = 500;
#else
  __constant__ PMEFloat small = 1e-6;
  __constant__ PMEFloat forceCut = 500;
#endif

}

//---------------------------------------------------------------------------------------------
// kgCalculateNetForce_kernel:
//
// Arguments:
//   TIRegion:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgCalculateNetForce_kernel(int TIRegion)
{
  using namespace NetForce;
  PMEFloat4 myNetFrc;
  unsigned int myWarpIgnored;
  bool p;
  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  myWarpIgnored =0;
  myNetFrc.x = 0.0;
  myNetFrc.y = 0.0;
  myNetFrc.z = 0.0;

  PMEAccumulator* pX, *pY, *pZ;

  if (TIRegion >= 0) {
    pX = cSim.pTINBForce + cSim.stride3*TIRegion;
  }
  else {
    pX = (PMEAccumulator*)cSim.pNBForceAccumulator;
  }
  pY = pX + cSim.stride;
  pZ = pX + cSim.stride2;

  PMEMask mask = __BALLOT(WARP_MASK, (pos < cSim.atoms) );
  while (pos < cSim.atoms) {
    PMEFloat fx = converter(pX[pos], ONEOVERFORCESCALE);
    PMEFloat fy = converter(pY[pos], ONEOVERFORCESCALE);
    PMEFloat fz = converter(pZ[pos], ONEOVERFORCESCALE);
    p = abs(fx) > small || abs(fy) > small || abs(fz) > small;
    if (p) {
      myNetFrc.x += (abs(fx) > forceCut) ? ((fx > 0) ? forceCut : -forceCut) : fx;
      myNetFrc.y += (abs(fy) > forceCut) ? ((fy > 0) ? forceCut : -forceCut) : fy;
      myNetFrc.z += (abs(fz) > forceCut) ? ((fz > 0) ? forceCut : -forceCut) : fz;
    }
    myWarpIgnored += maskPopc(__BALLOT(mask, !p));
    pos += increment;
    mask = __BALLOT(mask, (pos < cSim.atoms) );
  }
#pragma unroll
  for (int offset = warpSize / 2; offset > 0; offset /= 2) {
#ifndef AMBER_PLATFORM_AMD
    __syncwarp();
#endif
    myNetFrc.x += __SHFL_DOWN(0xFFFFFFFF, myNetFrc.x, offset);
    myNetFrc.y += __SHFL_DOWN(0xFFFFFFFF, myNetFrc.y, offset);
    myNetFrc.z += __SHFL_DOWN(0xFFFFFFFF, myNetFrc.z, offset);
  }

  if ((threadIdx.x & (warpSize - 1)) == 0){
    atomicAdd(&(netFrcX), ftoi(myNetFrc.x*eScale));
    atomicAdd(&(netFrcY), ftoi(myNetFrc.y*eScale));
    atomicAdd(&(netFrcZ), ftoi(myNetFrc.z*eScale));
    atomicAdd(&(ignored), myWarpIgnored);
  }
}

//---------------------------------------------------------------------------------------------
// kgRemoveNetForce_kernel:
//
// Arguments:
//   TIRegion:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgRemoveNetForce_kernel(int TIRegion)
{
  using namespace NetForce;
  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  PMEAccumulator* pX, *pY, *pZ;
  unsigned long long int  ifp = ftoi(forceCut *FORCESCALE);

  if (TIRegion >= 0) {
    pX = cSim.pTINBForce + (cSim.stride3 * TIRegion);
  }
  else {
    pX = cSim.pNBForceAccumulator;
  }
  pY = pX + cSim.stride;
  pZ = pX + cSim.stride * 2;

  PMEFloat normFactor = 1.0 / PMEFloat(cSim.atoms - ignored);
  PMEFloat ffx = converter(netFrcX, ONEOVERENERGYSCALE) * normFactor;
  PMEFloat ffy = converter(netFrcY, ONEOVERENERGYSCALE) * normFactor;
  PMEFloat ffz = converter(netFrcZ, ONEOVERENERGYSCALE) * normFactor;

  unsigned long long int nfX = ftoi(ffx * FORCESCALE);
  unsigned long long int nfY = ftoi(ffy * FORCESCALE);
  unsigned long long int nfZ = ftoi(ffz * FORCESCALE);

  while (pos < cSim.atoms) {
    PMEFloat fx = converter(pX[pos], ONEOVERFORCESCALE);
    PMEFloat fy = converter(pY[pos], ONEOVERFORCESCALE);
    PMEFloat fz = converter(pZ[pos], ONEOVERFORCESCALE);

    if (abs(fx) > small || abs(fy) > small || abs(fz) > small) {
      if (abs(fx) > forceCut) {
        pX[pos] = ((fx > 0) ? ifp : -ifp) - nfX;
      } else {
        pX[pos] -= nfX;
      }
      if (abs(fy) > forceCut) {
        pY[pos] = ((fy > 0) ? ifp : -ifp) - nfY;
      } else {
        pY[pos] -= nfY;
      }
      if (abs(fz) > forceCut) {
        pZ[pos] = ((fz > 0) ? ifp : -ifp) - nfZ;
      } else {
        pZ[pos] -= nfZ;
      }
    }
    pos += increment;
  }

  __syncthreads();

  if (threadIdx.x == 0) {
    unsigned int value = atomicInc(&blockCount, gridDim.x);
    if (value == (gridDim.x - 1)){
      netFrcX = 0;
      netFrcY = 0;
      netFrcZ = 0;
      ignored = 0;
      blockCount = 0;
    }
  }
}
//---------------------------------------------------------------------------------------------
// kgCopyToTIEnergy_kernel:
//
// Arguments:
//   TIRegion:
//   term1:
//   term2:
//   term3:
//   mode:
//   keepSource:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgCopyToTIEnergy_kernel(unsigned int TIRegion, int term1, int term2, int term3,
                                     bool isVirial, PMEFloat weight, bool addOn)
{
  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned long long int *pMain = cSim.pEnergyBuffer;
  unsigned long long int *pTI = cSim.pTIPotEnergyBuffer + (TIRegion * cSim.GPUPotEnergyTerms);


  if (pos == term1 || pos == term2 || pos == term3) {

    if (!addOn) pTI[pos] = Zero;
    if (TIRegion < 2) {
      PMEFloat tt = (isVirial) ? ONEOVERVIRIALSCALE : ONEOVERENERGYSCALE;
      PMEFloat ss = (isVirial) ? VIRIALSCALE : ENERGYSCALE;
      tt *= ((PMEFloat)pMain[pos]) * weight;
      pTI[pos] += (PMEAccumulator)(tt * ss);
    }
    else {
      pTI[pos] += pMain[pos];
    }
    pMain[pos] = ZeroF;
  }
}

//---------------------------------------------------------------------------------------------
// kgCorrectTIEnergy_kernel:
//
// Arguments:
//   beginTerm:
//   endTerm:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgCorrectTIEnergy_kernel(int beginTerm, int endTerm)
{
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned long long int *pMain = cSim.pEnergyBuffer;
  unsigned long long int *pTI0 = cSim.pTIPotEnergyBuffer;
  unsigned long long int *pTI1 = cSim.pTIPotEnergyBuffer + cSim.GPUPotEnergyTerms;
  unsigned long long int *pTI2 = cSim.pTIPotEnergyBuffer + cSim.GPUPotEnergyTerms * 2;

  if (pos >= beginTerm && pos <= endTerm) {
    pMain[pos] -= (pTI0[pos] + pTI1[pos]);
    pTI2[pos] += pMain[pos];
  }
}

//---------------------------------------------------------------------------------------------
// kgCorrectTIForce_kernel:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgCorrectTIForce_kernel()
{
  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  PMEAccumulator *pMain = cSim.pForceAccumulator;
  PMEAccumulator* pTI0 = cSim.pTIForce;
  PMEAccumulator* pTI1 = cSim.pTIForce + cSim.stride3;
  PMEAccumulator* pTIT = cSim.pTIForce + (cSim.stride3 * 2);

  while (pos < cSim.stride3) {
    pMain[pos] -= pTI0[pos] + pTI1[pos];
    pTIT[pos] += pMain[pos];
    pos += increment;
  }
}



#endif
