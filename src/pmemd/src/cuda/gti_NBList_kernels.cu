#ifdef GTI

#include <assert.h>
#include "gpuContext.h"
#include "gti_NBList_kernels.cuh"
#include "gti_utils.cuh"

#include "simulationConst.h"
CSIM_STO simulationConst cSim;

namespace GTI_NB_LIST {
  __device__ const int maxNumberTiles=1500;
  __device__ unsigned int possibleTITile[200][maxNumberTiles] = {0};
  __device__ unsigned int possible1264Tile[200][maxNumberTiles] = { 0 };
  __device__ unsigned int possiblep1264Tile[200][maxNumberTiles] = { 0 }; //C4PairwiseCUDA
  __device__ unsigned int possible1264p1264Tile[200][maxNumberTiles] = { 0 }; //C4PairwiseCUDA2023
  __device__ unsigned int numberPossibleTile[200] = {0};
  __constant__ bool useAllNB = true;
#  if !(defined(__CUDACC_RDC__) || defined(__HIPCC_RDC__))
#    include "gti_localCUutils.inc"
#  endif
}

using namespace GTI_NB;
using namespace GTI_NB_LIST;

//---------------------------------------------------------------------------------------------
// kgBuildSpecial2RestNBPreList_kernel:
//
// Arguments:
//   type:
//---------------------------------------------------------------------------------------------

_kPlainHead_
#if !defined(AMBER_PLATFORM_AMD)
#ifdef use_DPFP
__LAUNCH_BOUNDS__(128, 1)
#elif (__CUDA_ARCH__ >= 600)
__LAUNCH_BOUNDS__(256, 3)
#else
__LAUNCH_BOUNDS__(768, 1)
#endif
#endif
kgBuildSpecial2RestNBPreList_kernel(specialType type)
{
  if (useAllNB) return;

  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  PMEFloat cc2 = (cSim.cut + cSim.skinnb + 4.5) ;
  cc2 *= cc2;

  unsigned listSize = 0;
  unsigned *pList = NULL;
  unsigned (*pTiles)[maxNumberTiles] = NULL;
  switch(type) {
    case GTI_NB::TI : {
      listSize = cSim.numberUniqueTIAtoms;
      pList = cSim.pUniqueTIAtomList;
      pTiles = possibleTITile;
    }; break;
    case GTI_NB::LJ1264 : { 
      listSize = cSim.numberLJ1264Atoms;
      pList = cSim.pLJ1264AtomList;
      pTiles = possible1264Tile;
    }; break;
    case GTI_NB::pLJ1264 : { // C4PairwiseCUDA
      listSize = cSim.numberpLJ1264Atoms;
      pList = cSim.ppLJ1264AtomList;
      pTiles = possiblep1264Tile;
    }; break;
     case GTI_NB::LJ1264pLJ1264 : { // C4PairwiseCUDA2023
      listSize = cSim.numberLJ1264pLJ1264Atoms;
      pList = cSim.pLJ1264pLJ1264AtomList;
      pTiles = possible1264p1264Tile;
    }; break;
  }

  if (listSize == 0) {
    return;
  }

  unsigned th = threadIdx.x;
  while (th < listSize) {
    numberPossibleTile[th] = 0;
    th += increment;
  }
  __threadfence_block();

  __shared__ volatile float coordTI[gti_simulationConst::MaxNumberTIAtom][3];
  __shared__ PMEFloat recip[9], ucell[9];

  th = threadIdx.x;
  if (cSim.pNTPData == NULL) {
    if (th < 3) {
      recip[th*3    ] = cSim.recipf[th][0];
      recip[th*3 + 1] = cSim.recipf[th][1];
      recip[th*3 + 2] = cSim.recipf[th][2];
      ucell[th*3    ] = cSim.ucellf[th][0];
      ucell[th*3 + 1] = cSim.ucellf[th][1];
      ucell[th*3 + 2] = cSim.ucellf[th][2];
    }
  }
  else {
    if (th < 9) {
      ucell[th] = cSim.pNTPData->ucellf[th];
      recip[th] = cSim.pNTPData->recipf[th];
    }
  }

  th = threadIdx.x;
  while (th < listSize) {
    unsigned iatom = cSim.pImageAtomLookup[pList[th]];
    coordTI[th][0] = cSim.pImageX[iatom];
    coordTI[th][1] = cSim.pImageY[iatom];
    coordTI[th][2] = cSim.pImageZ[iatom];
    th += blockDim.x;
  }
  __syncthreads();

  while (pos < *(cSim.pNLEntries)) {

    NLEntry& p = cSim.pNLEntry[pos];
    if (p.NL.ymax & NLENTRY_HOME_CELL_MASK) {
      unsigned ypos = p.NL.ypos;
      unsigned ymax = p.NL.ymax >> NLENTRY_YMAX_SHIFT;

      PMEFloat x0 = cSim.pImageX[ypos];
      PMEFloat y0 = cSim.pImageY[ypos];
      PMEFloat z0 = cSim.pImageZ[ypos];
      PMEFloat x1 = cSim.pImageX[ymax-1];
      PMEFloat y1 = cSim.pImageY[ymax-1];
      PMEFloat z1 = cSim.pImageZ[ymax-1];

      for (unsigned i = 0; i < listSize; i++) {
        PMEFloat dx = coordTI[i][0] - x0;
        PMEFloat dy = coordTI[i][1] - y0;
        PMEFloat dz = coordTI[i][2] - z0;

        bool b1= __image_dist2_cut<PMEFloat>(dx, dy, dz, cc2, recip, ucell);
         dx = coordTI[i][0] - x1;
         dy = coordTI[i][1] - y1;
        dz = coordTI[i][2] - z1;
        bool b2 = __image_dist2_cut<PMEFloat>(dx, dy, dz, cc2, recip, ucell);
        if (b1 || b2) {
          unsigned tt = atomicAdd(&(numberPossibleTile[i]),1);
          pTiles[i][tt] = pos;
        }
      }

    }
    pos += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kgBuildSpecial2RestNBList_kernel:
//
// Arguments:
//   type:
//---------------------------------------------------------------------------------------------


_kPlainHead_
#if !defined(AMBER_PLATFORM_AMD)
#ifdef use_DPFP
__LAUNCH_BOUNDS__(128, 1)
#else
__LAUNCH_BOUNDS__(512, 2)
#endif
#endif
kgBuildSpecial2RestNBList_kernel(specialType type)
{
  using namespace GTI_NB;

  unsigned* pT = cSim.pTIList + cSim.stride * 2; // Index for non-TI atoms
  unsigned* pTI_REAF = cSim.pTIList + (1-cSim.reafMode)* cSim.stride; // Index for TI atoms
  PMEFloat cc =  (cSim.pNTPData == NULL) ? cSim.cutPlusSkin2 : cSim.pNTPData->cutPlusSkin2;

  unsigned listSize = 0;
  unsigned *pList = NULL;
  unsigned long long int *pNumberNBEntries = 0;
  int4* pNBList = NULL;
  switch (type) {
    case GTI_NB::TI: {
      listSize = cSim.numberUniqueTIAtoms;
      pList = cSim.pUniqueTIAtomList;
      pNumberNBEntries = cSim.pNumberTINBEntries;
      pNBList = cSim.pTINBList;
    }; break;
    case GTI_NB::LJ1264 : {
      listSize = cSim.numberLJ1264Atoms;
      pList = cSim.pLJ1264AtomList;
      pNumberNBEntries = cSim.pNumberLJ1264NBEntries;
      pNBList = cSim.pLJ1264NBList;
    }; break;
    case GTI_NB::pLJ1264 : { // C4PairwiseCUDA
      listSize = cSim.numberpLJ1264Atoms;
      pList = cSim.ppLJ1264AtomList;
      pNumberNBEntries = cSim.pNumberpLJ1264NBEntries;
      pNBList = cSim.ppLJ1264NBList;
    }; break;
    case GTI_NB::LJ1264pLJ1264 : { // C4PairwiseCUDA2023
      listSize = cSim.numberLJ1264pLJ1264Atoms;
      pList = cSim.pLJ1264pLJ1264AtomList;
      pNumberNBEntries = cSim.pNumberLJ1264pLJ1264NBEntries;
      pNBList = cSim.pLJ1264pLJ1264NBList;
    }; break;
    case GTI_NB::REAF: {
      listSize = cSim.numberREAFAtoms;
      pList = cSim.pREAFAtomList;
      pNumberNBEntries = cSim.pNumberREAFNbEntries;
      pNBList = cSim.pREAFNbList;
    }; break;
    assert(0); // Should not happen
  }

  __threadfence_block();

  if (listSize == 0) {
    return;
  }
  unsigned th = threadIdx.x;
  __shared__ PMEFloat recip[9], ucell[9];
  __shared__ float coordTI[gti_simulationConst::MaxNumberTIAtom][3];

  if (cSim.pNTPData == NULL) {
    if (th < 3) {
      recip[th*3    ] = cSim.recipf[th][0];
      recip[th*3 + 1] = cSim.recipf[th][1];
      recip[th*3 + 2] = cSim.recipf[th][2];
      ucell[th*3    ] = cSim.ucellf[th][0];
      ucell[th*3 + 1] = cSim.ucellf[th][1];
      ucell[th*3 + 2] = cSim.ucellf[th][2];
    }
  }
  else {
    if (th < 9) {
      ucell[th] = cSim.pNTPData->ucellf[th];
      recip[th] = cSim.pNTPData->recipf[th];
    }
  }

  __shared__ int tempList[gti_simulationConst::MaxNumberNBPerAtom];
  __shared__ unsigned counter;
  __shared__ unsigned int currentNLPos, currentNLPos2 ;
  __shared__ int iatom01, myregion;

  __syncthreads();

  unsigned bl = blockIdx.x;

  while (bl<listSize){

    if (threadIdx.x == 0) {
      counter = 0;
      tempList[0] = -1;
    }
    __syncthreads();

    bool isLinearTI = (type == GTI_NB::TI && bl < cSim.numberTICommonPairs);

    unsigned atom0;
    if (isLinearTI) {
      atom0 = cSim.pTICommonPair[bl].x;
    } else {
      atom0 = pList[bl];
    }

    unsigned iatom0 = cSim.pImageAtomLookup[atom0];
    float coord[3];
    coord[0] = cSim.pImageX[iatom0];
    coord[1] = cSim.pImageY[iatom0];
    coord[2] = cSim.pImageZ[iatom0];

    unsigned int limit = (useAllNB) ? cSim.atoms : numberPossibleTile[bl] * cSim.NLAtomsPerWarp;
    unsigned iatom1 = threadIdx.x;

    while (iatom1 < limit) {
      PMEFloat dx = cSim.pImageX[iatom1] - coord[0];
      PMEFloat dy = cSim.pImageY[iatom1] - coord[1];
      PMEFloat dz = cSim.pImageZ[iatom1] - coord[2];;
      if (__image_dist2_cut<PMEFloat>(dx, dy, dz, cc, recip, ucell)) {
        unsigned atom1 = cSim.pImageAtom[iatom1];
          switch (type) {
            case GTI_NB::TI: {
              if (cSim.numberTIAtoms > 0 && pT[atom1]>0) {
                unsigned mycount = atomicAdd(&counter, 1);
                tempList[mycount] = iatom1;
              }
            }; break;
            case GTI_NB::LJ1264 : {
              bool test=(atom1 != atom0);
              if (test && cSim.numberTIAtoms > 0) {
                if (pT[atom1] == 0 || pT[atom0] == 0) test = false;
              }
              if (test) {
                unsigned mycount = atomicAdd(&counter, 1);
                tempList[mycount] = iatom1;
              }
            }; break;
            case GTI_NB::pLJ1264 : {  // C4PairwiseCUDA
              bool test=(atom1 != atom0);
              if (test && cSim.numberTIAtoms > 0) {
                if (pT[atom1] == 0 || pT[atom0] == 0) test = false;
              }
              if (test) {
                unsigned mycount = atomicAdd(&counter, 1);
                tempList[mycount] = iatom1;
              }
            }; break;
            case GTI_NB::LJ1264pLJ1264 : {  // C4PairwiseCUDA
              bool test=(atom1 != atom0);
              if (test && cSim.numberTIAtoms > 0) {
                if (pT[atom1] == 0 || pT[atom0] == 0) test = false;
              }
              if (test) {
                unsigned mycount = atomicAdd(&counter, 1);
                tempList[mycount] = iatom1;
              }
            }; break;
            case GTI_NB::REAF : {
              bool test = (atom1 != atom0);
              if (test && cSim.numberTIAtoms > 0) { // exclude these already are in TI region
                if (pTI_REAF[atom1]>0) test = false;
              }
              if (test && cSim.pREAFList[atom1] > 0) test = atom0 > atom1;
              if (test) {
                unsigned mycount = atomicAdd(&counter, 1);
                tempList[mycount] = iatom1;
              }
            }; break;
          }  // switch
      }  // if (__image_dist2_c..
      iatom1 += blockDim.x;
    }

    __syncthreads();

    const int shift = gti_simulationConst::MaxNumberNBPerAtom;

    if (threadIdx.x == 0) {
      iatom01 = -1; myregion = -1;
      if (type == GTI_NB::TI) {
        if (isLinearTI) {
          currentNLPos = bl * shift;
          currentNLPos2 = (bl + cSim.numberTICommonPairs) * shift;
          iatom01 = cSim.pImageAtomLookup[cSim.pTICommonPair[bl].y];
          myregion = 0;
        } else {
          currentNLPos = (bl+ cSim.numberTICommonPairs)* shift;
          myregion = ((bl - cSim.numberTICommonPairs) < cSim.softcoreAtomShift) ? 0 : 1;
        }
      } else {
        currentNLPos = bl * shift;
      }

      pNBList[currentNLPos ].x = iatom0;
      pNBList[currentNLPos ].y = -1;
      pNBList[currentNLPos ].z = myregion;
      pNBList[currentNLPos + shift - 1].x = counter;
      pNBList[currentNLPos + shift - 1].y = -1;
      if (isLinearTI) {
        pNBList[currentNLPos2 ].x = iatom01;
        pNBList[currentNLPos2 ].y = -1;
        pNBList[currentNLPos2 ].z = 1;
        pNBList[currentNLPos2 + shift - 1].x = counter;
        pNBList[currentNLPos2 + shift - 1].y = -1;
      }

    }

    __syncthreads();

    if (counter > 0) {
      unsigned th = threadIdx.x;
      while (th < shift - 1) {
        if (th < counter) {
          pNBList[currentNLPos + th].x = iatom0;
          pNBList[currentNLPos + th].y = tempList[th];
          pNBList[currentNLPos + th].z = 0;
          if (isLinearTI) {
            pNBList[currentNLPos2 + th].x = iatom01;
            pNBList[currentNLPos2 + th].y = tempList[th];
            pNBList[currentNLPos2 + th].z = 1;
          }
          else {
            pNBList[currentNLPos + th].z = myregion;
          }
        }
        else {
          pNBList[currentNLPos + th] = { -1, -1, -1, -1 };
          if (isLinearTI) pNBList[currentNLPos2 + th] = { -1, -1, -1, -1 };
        }

        th += blockDim.x;
      }
    }

    __syncthreads();

    bl += gridDim.x;
  }
}


//---------------------------------------------------------------------------------------------
// kgBuildTI2TINBList_kernel:
//---------------------------------------------------------------------------------------------
_kPlainHead_
#ifdef use_DPFP
__launch_bounds__(128, 1)
#elif (__CDA_ARCH__ >= 600)
__launch_bounds__(64, 16)
#else
__launch_bounds__(1024, 1)
#endif
kgBuildTI2TINBList_kernel() {

  PMEDouble* pX = cSim.pImageX;
  PMEDouble* pY = cSim.pImageY;
  PMEDouble* pZ = cSim.pImageZ;

  const int shift = gti_simulationConst::MaxNumberNBPerAtom;

  __shared__ volatile float coordTI[gti_simulationConst::MaxNumberTIAtom][3];

  unsigned th = threadIdx.x;
  __shared__ PMEFloat recip[9], ucell[9];

  if (cSim.pNTPData == NULL) {
    if (th < 3) {
      recip[th * 3] = cSim.recipf[th][0];
      recip[th * 3 + 1] = cSim.recipf[th][1];
      recip[th * 3 + 2] = cSim.recipf[th][2];
      ucell[th * 3] = cSim.ucellf[th][0];
      ucell[th * 3 + 1] = cSim.ucellf[th][1];
      ucell[th * 3 + 2] = cSim.ucellf[th][2];
    }
  } else {
    if (th < 9) {
      ucell[th] = cSim.pNTPData->ucellf[th];
      recip[th] = cSim.pNTPData->recipf[th];
    }
  }
  th = threadIdx.x;
  while (th < cSim.numberTIAtoms) {
    int qq = th * shift;
    int iatom = cSim.pTINBList[qq].x;
    if (iatom >= 0) {
      coordTI[th][0] = pX[iatom];
      coordTI[th][1] = pY[iatom];
      coordTI[th][2] = pZ[iatom];
    }
    th += blockDim.x;
  }
  __syncthreads();

  unsigned totalPair = cSim.numberTIAtoms * (cSim.numberTIAtoms + 1) / 2;
  PMEFloat cc = (cSim.pNTPData == NULL) ? cSim.cutPlusSkin2 : cSim.pNTPData->cutPlusSkin2;
  PMEFloat r2 = ZeroF;
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  while (pos < totalPair) {
    unsigned index0 = (sqrt(pos * 8 + 1.0001) - 1) / 2;
    unsigned index1 = pos - index0 * (index0 + 1) / 2;

    int4 yy0 = cSim.pTINBList[index0 * shift];
    int4 yy1 = cSim.pTINBList[index1 * shift];
    unsigned region0 = yy0.z & gti_simulationConst::Fg_TI_region;
    unsigned region1 = yy1.z & gti_simulationConst::Fg_TI_region;
    PMEFloat r2 = ZeroF;
    if (region0 == region1 && index0 != index1) { // must be in the same TI region
      if (cSim.tiCut == 0) {
        PMEFloat dx = coordTI[index0][0] - coordTI[index1][0];
        PMEFloat dy = coordTI[index0][1] - coordTI[index1][1];
        PMEFloat dz = coordTI[index0][2] - coordTI[index1][2];
        r2 = __image_dist2<PMEFloat>(dx, dy, dz, recip, ucell);
      }
      if ((r2 < cc) || cSim.tiCut > 0) {
        unsigned iatom0 = yy0.x;
        unsigned iatom1 = yy1.x;
        int qq = (index0 + 1) * shift - 1;
        int currentNLPos = atomicAdd(&(cSim.pTINBList[qq].x), 1) + index0 * shift;
        cSim.pTINBList[currentNLPos].x = iatom0;
        cSim.pTINBList[currentNLPos].y = iatom1;
        cSim.pTINBList[currentNLPos].z = region0;
        cSim.pTINBList[currentNLPos].z |= gti_simulationConst::Fg_int_TI;
      }
    }
    pos += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kgTINBListFillAttribute_kernel:
//---------------------------------------------------------------------------------------------
_kPlainHead_
#if !defined(AMBER_PLATFORM_AMD)
#ifdef use_DPFP
__LAUNCH_BOUNDS__(128, 1)
#elif (__CUDA_ARCH__ >= 600)
__LAUNCH_BOUNDS__(64, 16)
#else
__LAUNCH_BOUNDS__(1024, 1)
#endif
#endif
kgTINBListFillAttribute_kernel()
{
  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  while (pos < *cSim.pNumberTINBEntries) {
    int4* nbEntry = cSim.pTINBList;
    int iatom0 = nbEntry[pos].x;
    int iatom1 = nbEntry[pos].y;
    if (iatom0 >= 0 && iatom1 >= 0) {
      int atom0 = cSim.pImageAtom[iatom0];
      int atom1 = cSim.pImageAtom[iatom1];
      nbEntry[pos].w = -1;

      bool found = false;  // Check if atom0 is on atom1's exclusion list
      for (unsigned kk = cSim.pNLExclusionStartCount[atom1].x;
        kk < cSim.pNLExclusionStartCount[atom1].x + cSim.pNLExclusionStartCount[atom1].y; kk++) {
        if (atom0 == cSim.pNLExclusionList[kk]) {
          found = true;
          break;
        }
      }
      if (found) {
        nbEntry[pos].z |= gti_simulationConst::Fg_excluded;
      } else {
        int tt = cSim.TIVdwNTyes * (cSim.pTIac[atom0] - 1) + (cSim.pTIac[atom1] - 1);
        int vdwindex = (tt >= 0) ? cSim.pTIico[tt] - 1 : -1;
        if (vdwindex >= 0) {
          nbEntry[pos].w = vdwindex;
          nbEntry[pos].z |= gti_simulationConst::ex_need_LJ;
        }
      }

      if (cSim.pSCList[atom0] > 0 || cSim.pSCList[atom1] > 0) {
        nbEntry[pos].z |= gti_simulationConst::Fg_has_SC;
        if (cSim.pSCList[atom0] > 0 && cSim.pSCList[atom1] > 0) {
          nbEntry[pos].z |= gti_simulationConst::Fg_int_SC;
        }
      }

      if (nbEntry[pos].z & gti_simulationConst::Fg_excluded) {
        nbEntry[pos].z |= gti_simulationConst::ex_addToDVDL_ele; // excluded pairs still have "reverse-PME" contribution to be corrected.
      } else {
        if (nbEntry[pos].z & gti_simulationConst::Fg_has_SC) {
          if (nbEntry[pos].z & gti_simulationConst::Fg_int_SC) {
            if (cSim.addSC == 2 || cSim.addSC == 3 || cSim.addSC == 25 || cSim.addSC == 35|| cSim.addSC == 5 || cSim.addSC == 6) {
              nbEntry[pos].z |= gti_simulationConst::ex_addToDVDL_ele;
              nbEntry[pos].z |= gti_simulationConst::ex_SC_ELE;
              if (cSim.addSC == 3 || cSim.addSC == 35 || cSim.addSC == 6) {
                nbEntry[pos].z |= gti_simulationConst::ex_addToDVDL_vdw;
                nbEntry[pos].z |= gti_simulationConst::ex_SC_LJ;
              }
            }
          } else {
            nbEntry[pos].z |= gti_simulationConst::ex_addToDVDL_ele;
            nbEntry[pos].z |= gti_simulationConst::ex_addToDVDL_vdw;
            nbEntry[pos].z |= gti_simulationConst::ex_SC_LJ;
            nbEntry[pos].z |= gti_simulationConst::ex_SC_ELE;
          }
        } else {
          nbEntry[pos].z |= gti_simulationConst::ex_addToDVDL_ele;
          nbEntry[pos].z |= gti_simulationConst::ex_addToDVDL_vdw;
        }
      }
    }
    pos += increment;
  }
}



//---------------------------------------------------------------------------------------------
// kgREAFNbListFillAttribute_kernel:
//---------------------------------------------------------------------------------------------

_kPlainHead_
#if !defined(AMBER_PLATFORM_AMD)
#ifdef use_DPFP
__LAUNCH_BOUNDS__(128, 1)
#elif (__CUDA_ARCH__ >= 600)
__LAUNCH_BOUNDS__(64, 16)
#else
__LAUNCH_BOUNDS__(1024, 1)
#endif
#endif
kgREAFNbListFillAttribute_kernel()
{
  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  while (pos < *cSim.pNumberREAFNbEntries) {
    int4* nbEntry = cSim.pREAFNbList;
    int iatom0 = nbEntry[pos].x;
    int iatom1 = nbEntry[pos].y;
    if (iatom0 >= 0 && iatom1 >= 0) {
      int atom0 = cSim.pImageAtom[iatom0];
      int atom1 = cSim.pImageAtom[iatom1];
      nbEntry[pos].w = -1;

      bool found = false;  // Check if atom0 is on atom1's exclusion list
      for (unsigned kk = cSim.pNLExclusionStartCount[atom1].x;
        kk < cSim.pNLExclusionStartCount[atom1].x + cSim.pNLExclusionStartCount[atom1].y; kk++) {
        if (atom0 == cSim.pNLExclusionList[kk]) {
          found = true;
          break;
        }
      }
      if (found) {
        nbEntry[pos].z |= gti_simulationConst::Fg_excluded;
      }
      else {
        int tt = cSim.TIVdwNTyes * (cSim.pTIac[atom0] - 1) + (cSim.pTIac[atom1] - 1);
        int vdwindex = (tt >= 0) ? cSim.pTIico[tt] - 1 : -1;
        if (vdwindex >= 0) {
          nbEntry[pos].w = vdwindex;
          nbEntry[pos].z |= gti_simulationConst::ex_need_LJ;
        }
      }

      if (cSim.pREAFList[atom0] > 0 || cSim.pREAFList[atom1] > 0) {
        nbEntry[pos].z |= gti_simulationConst::Fg_has_RE;
        if (cSim.pREAFList[atom0] > 0 && cSim.pREAFList[atom1] > 0) {
          nbEntry[pos].z |= gti_simulationConst::Fg_int_RE;
        }
      }

      if (nbEntry[pos].z & gti_simulationConst::Fg_excluded) {
        nbEntry[pos].z |= gti_simulationConst::ex_addRE_ele; // excluded pairs still have "reverse-PME" contribution to be corrected.
      }
      else {
        if (nbEntry[pos].z & gti_simulationConst::Fg_has_RE) {
          if (nbEntry[pos].z & gti_simulationConst::Fg_int_RE) {
            if (cSim.addRE == 2 || cSim.addRE == 3 || cSim.addRE == 5 || cSim.addRE == 6 || cSim.addRE == 7) {
              nbEntry[pos].z |= gti_simulationConst::ex_addRE_ele;
              if (cSim.addRE == 3 || cSim.addRE == 6 || cSim.addRE == 7 ) {
                nbEntry[pos].z |= gti_simulationConst::ex_addRE_vdw;
              }
            }
          }else {
            nbEntry[pos].z |= gti_simulationConst::ex_addRE_ele;
            nbEntry[pos].z |= gti_simulationConst::ex_addRE_vdw;
          }
        }
        else {
          nbEntry[pos].z |= gti_simulationConst::ex_addRE_ele;
          nbEntry[pos].z |= gti_simulationConst::ex_addRE_vdw;
        }
      }
    }
    pos += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kg1264NBListFillAttribute_kernel:
//---------------------------------------------------------------------------------------------
_kPlainHead_
#if !defined(AMBER_PLATFORM_AMD)
#ifdef use_DPFP
__launch_bounds__(128, 1)
#elif (__CDA_ARCH__ >= 600)
__launch_bounds__(64, 16)
#else
__LAUNCH_BOUNDS__(1024, 1)
#endif
#endif
kg1264NBListFillAttribute_kernel()
{
  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  int num = cSim.numberLJ1264Atoms * gti_simulationConst::MaxNumberNBPerAtom;
  int4* nbEntry = cSim.pLJ1264NBList;

  while (pos < num) {
    int iatom0 = nbEntry[pos].x;
    int iatom1 = nbEntry[pos].y;
    if (iatom0 >= 0 && iatom1 >= 0) {
      int atom0 = cSim.pImageAtom[iatom0];
      int atom1 = cSim.pImageAtom[iatom1];
      nbEntry[pos].w = -1;
      nbEntry[pos].z = 0;
      int tt = cSim.TIVdwNTyes * (cSim.pTIac[atom0] - 1) + (cSim.pTIac[atom1] - 1);
      int vdwindex = (tt >= 0) ? cSim.pTIico[tt] - 1 : -1;
      if (vdwindex >= 0) {
        nbEntry[pos].w = vdwindex;
        nbEntry[pos].z |= gti_simulationConst::ex_need_LJ;
      }
    }
    pos += increment;
  }
}

//========================================================

//---------------------------------------------------------------------------------------------
// kgp1264NBListFillAttribute_kernel: // C4PairwiseCUDA
//---------------------------------------------------------------------------------------------
// C4PairwiseCUDA
_kPlainHead_
#ifdef use_DPFP
__launch_bounds__(128, 1)
#elif (__CDA_ARCH__ >= 600)
__launch_bounds__(64, 16)
#else
__launch_bounds__(1024, 1)
#endif
kgp1264NBListFillAttribute_kernel() //C4PairwiseCUDA2023
{
  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  int num = cSim.TIC4Pairwise * 2 * gti_simulationConst::MaxNumberNBPerAtom;
  int4* nbEntry = cSim.ppLJ1264NBList;

  while (pos < num) {
    int iatom0 = nbEntry[pos].x;
    int iatom1 = nbEntry[pos].y;
    if (iatom0 >= 0 && iatom1 >= 0) {
      int atom0 = cSim.pImageAtom[iatom0];
      int atom1 = cSim.pImageAtom[iatom1];
      nbEntry[pos].w = 55555;
      nbEntry[pos].z |= gti_simulationConst::ex_need_LJ;
    }
    pos += increment;
  }
}

//C4PairwiseCUDA2023

_kPlainHead_
#ifdef use_DPFP
__launch_bounds__(128, 1)
#elif (__CDA_ARCH__ >= 600)
__launch_bounds__(64, 16)
#else
__launch_bounds__(1024, 1)
#endif
kg1264p1264NBListFillAttribute_kernel()
{
  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  int num1 = cSim.TIC4Pairwise * 2 * gti_simulationConst::MaxNumberNBPerAtom;
  int num2 = cSim.numberLJ1264pLJ1264Atoms * gti_simulationConst::MaxNumberNBPerAtom;
  int4* nbEntry = cSim.pLJ1264pLJ1264NBList;

  while (pos < num1) {
    int iatom0 = nbEntry[pos].x;
    int iatom1 = nbEntry[pos].y;
    if (iatom0 >= 0 && iatom1 >= 0) {
      int atom0 = cSim.pImageAtom[iatom0];
      int atom1 = cSim.pImageAtom[iatom1];
      nbEntry[pos].w = 55555;
      nbEntry[pos].z |= gti_simulationConst::ex_need_LJ;
    }
    pos += increment;
  }//C4PairwiseCUDA2023

  while (pos < num2) {
    int iatom0 = nbEntry[pos].x;
    int iatom1 = nbEntry[pos].y;
    if (iatom0 >= 0 && iatom1 >= 0) {
      int atom0 = cSim.pImageAtom[iatom0];
      int atom1 = cSim.pImageAtom[iatom1];
      nbEntry[pos].w = -1;
      nbEntry[pos].z = 0;
      int tt = cSim.TIVdwNTyes * (cSim.pTIac[atom0] - 1) + (cSim.pTIac[atom1] - 1);
      int vdwindex = (tt >= 0) ? cSim.pTIico[tt] - 1 : -1;
      if (vdwindex >= 0) {
        nbEntry[pos].w = vdwindex;
        nbEntry[pos].z |= gti_simulationConst::ex_need_LJ;
      }
    }
    pos += increment;
  }
}

//========================================================

//---------------------------------------------------------------------------------------------
// kgBuildSpecial2RestNBPreList_gamd_kernel:
//
// Arguments:
//   type:
//---------------------------------------------------------------------------------------------
_kPlainHead_
#if !defined(AMBER_PLATFORM_AMD)
#ifdef use_DPFP
__LAUNCH_BOUNDS__(128,1)
#elif (__CUDA_ARCH__ >= 600)
__LAUNCH_BOUNDS__(256, 3)
#else
__LAUNCH_BOUNDS__(768, 1)
#endif
#endif
kgBuildSpecial2RestNBPreList_gamd_kernel(specialType type)
{
  if (useAllNB) return;

  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  PMEFloat cc2 = (cSim.cut + cSim.skinnb + 4.5) ;
  cc2 *= cc2;

  unsigned listSize = 0;
  unsigned *pList = NULL;
  unsigned (*pTiles)[maxNumberTiles] = NULL;
  switch(type) {
    case GTI_NB::TI : {
      listSize = cSim.numberUniqueTIAtoms;
      pList = cSim.pUniqueTIAtomList;
      pTiles = possibleTITile;
    }; break;
    case GTI_NB::LJ1264 : {
      listSize = cSim.numberLJ1264Atoms;
      pList = cSim.pLJ1264AtomList;
      pTiles = possible1264Tile;
    }; break;
    case GTI_NB::pLJ1264 : { // C4PairwiseCUDA
      listSize = cSim.numberpLJ1264Atoms;
      pList = cSim.ppLJ1264AtomList;
      pTiles = possiblep1264Tile;
    }; break;
    case GTI_NB::LJ1264pLJ1264 : { // C4PairwiseCUDA2023
      listSize = cSim.numberLJ1264pLJ1264Atoms;
      pList = cSim.pLJ1264pLJ1264AtomList;
      pTiles = possible1264p1264Tile;
    }; break;
  }

  if (listSize == 0) {
    return;
  }

  unsigned i = threadIdx.x;
  while (i < listSize) {
    numberPossibleTile[i] = 0;
    i += increment;
  }
  __threadfence_block();

  float *coordTIX = NULL;
  float *coordTIY = NULL;
  float *coordTIZ = NULL;
  // __shared__ volatile float coordTI[cSim.atoms][3];
  __shared__ PMEFloat recip[9], ucell[9];

  i = threadIdx.x;
  if (cSim.pNTPData == NULL) {
    if (i < 3) {
      recip[i*3    ] = cSim.recipf[i][0];
      recip[i*3 + 1] = cSim.recipf[i][1];
      recip[i*3 + 2] = cSim.recipf[i][2];
      ucell[i*3    ] = cSim.ucellf[i][0];
      ucell[i*3 + 1] = cSim.ucellf[i][1];
      ucell[i*3 + 2] = cSim.ucellf[i][2];
    }
  }
  else {
    if (i < 9) {
      ucell[i] = cSim.pNTPData->ucellf[i];
      recip[i] = cSim.pNTPData->recipf[i];
    }
  }

  i = threadIdx.x;
  while (i < listSize) {
    unsigned iatom = cSim.pImageAtomLookup[pList[i]];
    coordTIX[i] = cSim.pImageX[iatom];
    coordTIY[i] = cSim.pImageY[iatom];
    coordTIZ[i] = cSim.pImageZ[iatom];
    i += blockDim.x;
  }
  __syncthreads();

  while (pos < *(cSim.pNLEntries)) {

    NLEntry& p = cSim.pNLEntry[pos];
    if (p.NL.ymax & NLENTRY_HOME_CELL_MASK) {
      unsigned ypos = p.NL.ypos;
      unsigned ymax = p.NL.ymax >> NLENTRY_YMAX_SHIFT;

      PMEFloat x0 = cSim.pImageX[ypos];
      PMEFloat y0 = cSim.pImageY[ypos];
      PMEFloat z0 = cSim.pImageZ[ypos];
      PMEFloat x1 = cSim.pImageX[ymax-1];
      PMEFloat y1 = cSim.pImageY[ymax-1];
      PMEFloat z1 = cSim.pImageZ[ymax-1];

      for (unsigned i = 0; i < listSize; i++) {
        PMEFloat dx = coordTIX[i] - x0;
        PMEFloat dy = coordTIY[i] - y0;
        PMEFloat dz = coordTIZ[i] - z0;

        bool b1= __image_dist2_cut<PMEFloat>(dx, dy, dz, cc2, recip, ucell);
         dx = coordTIX[i] - x1;
         dy = coordTIY[i] - y1;
         dz = coordTIZ[i] - z1;
        bool b2 = __image_dist2_cut<PMEFloat>(dx, dy, dz, cc2, recip, ucell);
        if (b1 || b2) {
          unsigned tt = atomicAdd(&(numberPossibleTile[i]),1);
          pTiles[i][tt] = pos;
        }
      }

    }
    pos += increment;
  }
}


//---------------------------------------------------------------------------------------------
// kgBuildSpecial2RestNBList_gamd_kernel:
//
// Arguments:
//   type:
//---------------------------------------------------------------------------------------------
_kPlainHead_
#if !defined(AMBER_PLATFORM_AMD)
#ifdef use_DPFP
__LAUNCH_BOUNDS__(128, 1)
#else
__LAUNCH_BOUNDS__(512, 2)
#endif
#endif
kgBuildSpecial2RestNBList_gamd_kernel(specialType type)
{
  using namespace GTI_NB;

  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;

  unsigned* pT = cSim.pTIList + cSim.stride * 2; // Index for non-TI atoms
  PMEFloat cc =  (cSim.pNTPData == NULL) ? cSim.cutPlusSkin2 : cSim.pNTPData->cutPlusSkin2;

  unsigned listSize = 0;
  unsigned *pList = NULL;
  unsigned long long int *pNumberNBEntries = 0;
  int4* pNBList = NULL;
  switch (type) {
    case GTI_NB::TI: {
      listSize = cSim.numberUniqueTIAtoms;
      pList = cSim.pUniqueTIAtomList;
      pNumberNBEntries = cSim.pNumberTINBEntries;
      pNBList = cSim.pTINBList;
    }; break;
    case GTI_NB::LJ1264: {
      listSize = cSim.numberLJ1264Atoms;
      pList = cSim.pLJ1264AtomList;
      pNumberNBEntries = cSim.pNumberLJ1264NBEntries;
      pNBList = cSim.pLJ1264NBList;
    }; break;
    case GTI_NB::pLJ1264 : { // C4PairwiseCUDA
      listSize = cSim.numberpLJ1264Atoms;
      pList = cSim.ppLJ1264AtomList;
      pNumberNBEntries = cSim.pNumberpLJ1264NBEntries;
      pNBList = cSim.ppLJ1264NBList;
    }; break;
    case GTI_NB::LJ1264pLJ1264 : { // C4PairwiseCUDA2023
      listSize = cSim.numberLJ1264pLJ1264Atoms;
      pList = cSim.pLJ1264pLJ1264AtomList;
      pNumberNBEntries = cSim.pNumberLJ1264pLJ1264NBEntries;
      pNBList = cSim.pLJ1264pLJ1264NBList;
    }; break;
  }

  if (pos == 0) {
    *(pNumberNBEntries) = 0;
  }
  __threadfence_block();

  if (listSize == 0) {
    return;
  }
  unsigned i = threadIdx.x;
  __shared__ PMEFloat recip[9], ucell[9];

  if (cSim.pNTPData == NULL) {
    if (i < 3) {
      recip[i*3    ] = cSim.recipf[i][0];
      recip[i*3 + 1] = cSim.recipf[i][1];
      recip[i*3 + 2] = cSim.recipf[i][2];
      ucell[i*3    ] = cSim.ucellf[i][0];
      ucell[i*3 + 1] = cSim.ucellf[i][1];
      ucell[i*3 + 2] = cSim.ucellf[i][2];
    }
  }
  else {
    if (i < 9) {
      ucell[i] = cSim.pNTPData->ucellf[i];
      recip[i] = cSim.pNTPData->recipf[i];
    }
  }
  __shared__ int tempList[2000];
  __shared__ unsigned counter;
  __shared__ volatile unsigned int currentNLPos, currentNLPos2 ;

  __syncthreads();

  i = blockIdx.x;

  while (i<listSize){

    if (threadIdx.x == 0) {
      counter = 0;
      currentNLPos = *pNumberNBEntries;
    }
    __syncthreads();

    unsigned atom1 = pList[i];
    unsigned iatom = cSim.pImageAtomLookup[atom1];
    float coord[3];
    coord[0] = cSim.pImageX[iatom];
    coord[1] = cSim.pImageY[iatom];
    coord[2] = cSim.pImageZ[iatom];

    unsigned int limit = (useAllNB) ? cSim.atoms : numberPossibleTile[i] * cSim.NLAtomsPerWarp;

    int atom11 = -1, myregion = -1;
    if (type == GTI_NB::TI) {
      if (i < cSim.numberTICommonPairs) atom11 = cSim.pTICommonPair[i].y;
      else myregion = ((i - cSim.numberTICommonPairs) < cSim.softcoreAtomShift) ? 0 : 1;
    }

    myregion |= 26;  // 2+8+16 ; 2 (need LJ) and 8 (need SC-LJ)  16 (need SC-ELE)

    unsigned iatom0 = threadIdx.x;
    while (iatom0 < limit) {
        PMEFloat dx = cSim.pImageX[iatom0] - coord[0];
        PMEFloat dy = cSim.pImageY[iatom0] - coord[1];
        PMEFloat dz = cSim.pImageZ[iatom0] - coord[2];;
        if (__image_dist2_cut<PMEFloat>(dx, dy, dz, cc, recip, ucell)) {
          unsigned atom0 = cSim.pImageAtom[iatom0];
          switch (type) {
            case GTI_NB::TI: {
              if (pT[atom0] > 0) {
                unsigned mycount = atomicAdd(&counter, 1);
                tempList[mycount] = atom0;
              }
            }; break;
            case GTI_NB::LJ1264: {
              if (atom1 != atom0) {
                unsigned mycount = atomicAdd(&counter, 1);
                tempList[mycount] = atom0;
              }
            }; break;
            case GTI_NB::pLJ1264 : { // C4PairwiseCUDA
              if (atom1 != atom0) {
                unsigned mycount = atomicAdd(&counter, 1);
                tempList[mycount] = atom0;
              }
            }; break;
            case GTI_NB::LJ1264pLJ1264 : { // C4PairwiseCUDA
              if (atom1 != atom0) {
                unsigned mycount = atomicAdd(&counter, 1);
                tempList[mycount] = atom0;
              }
            }; break;
          }  // switch
      }  // if (__image_dist2_c..
      iatom0 += blockDim.x;
    }

    __syncthreads();

    bool isLinearTI = (type == GTI_NB::TI && i < cSim.numberTICommonPairs);

    if (threadIdx.x == 0 && counter > 0) {
      currentNLPos = atomicAdd(pNumberNBEntries, counter);
      if (isLinearTI) {
        currentNLPos2 = atomicAdd(pNumberNBEntries, counter);
      }
    }
    __syncthreads();

    unsigned m = threadIdx.x;
    while (m < counter) {
      pNBList[currentNLPos + m].x = atom1;
      pNBList[currentNLPos + m].y = tempList[m];
      pNBList[currentNLPos + m].z = 0;
      if (isLinearTI) {
        pNBList[currentNLPos2 + m].x = atom11;
        pNBList[currentNLPos2 + m].y = tempList[m];
        pNBList[currentNLPos2 + m].z = 1;
      } else {
        pNBList[currentNLPos + m].z = myregion;
      }
      m += blockDim.x;
    }

    i += gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// TL:This kernel is out-of-date and is kept for gamd
// kgBuildTI2RestNBListFillAttribute_kernel:
//---------------------------------------------------------------------------------------------
_kPlainHead_
#if !defined(AMBER_PLATFORM_AMD)
#ifdef use_DPFP
__launch_bounds__(128, 1)
#elif (__CDA_ARCH__ >= 600)
__launch_bounds__(64, 16)
#else
__LAUNCH_BOUNDS__(1024, 1)
#endif
#endif
kgBuildTI2RestNBListFillAttribute_kernel()
{
  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  while (pos < *cSim.pNumberTINBEntries) {
    int4* nbEntry = cSim.pTINBList;
    int atom0 = nbEntry[pos].x;
    int atom1 = nbEntry[pos].y;
    nbEntry[pos].w = -1;
    bool found = false;  // Check if atom0 is on atom1's exclusion list
    for (unsigned kk = cSim.pNLExclusionStartCount[atom1].x;
      kk < cSim.pNLExclusionStartCount[atom1].x + cSim.pNLExclusionStartCount[atom1].y; kk++) {
      if (atom0 == cSim.pNLExclusionList[kk]) {
        found = true;
        break;
      }
    }
    if (found) {
      nbEntry[pos].z |= 4;
    } else {
      int tt = cSim.TIVdwNTyes * (cSim.pTIac[atom0] - 1) + (cSim.pTIac[atom1] - 1);
      int vdwindex = (tt >= 0) ? cSim.pTIico[tt] - 1 : -1;
      if (vdwindex >= 0) {
        nbEntry[pos].w = vdwindex;
        nbEntry[pos].z |= 2;
      }
    }
    pos += increment;

  }
}

//---------------------------------------------------------------------------------------------
// TL:This kernel is out-of-date and is kept for gamd
// kgBuild1264NBListFillAttribute_kernel:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgBuild1264NBListFillAttribute_kernel()
{
  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  while (pos < *cSim.pNumberLJ1264NBEntries) {
    int4* nbEntry = cSim.pLJ1264NBList;
    int iatom0 = nbEntry[pos].x;
    int iatom1 = nbEntry[pos].y;
    if (iatom0 >= 0 && iatom1 >= 0) {
      nbEntry[pos].w = -1;
      nbEntry[pos].z = 0;
      int atom0 = cSim.pImageAtom[iatom0];
      int atom1 = cSim.pImageAtom[iatom1];
      int tt = cSim.TIVdwNTyes * (cSim.pTIac[atom0] - 1) + (cSim.pTIac[atom1] - 1);
      int vdwindex = (tt >= 0) ? cSim.pTIico[tt] - 1 : -1;
      if (vdwindex >= 0) {
        nbEntry[pos].w = vdwindex;
        nbEntry[pos].z |= 2;
      }
    }
    pos += increment;
  }
}


//---------------------------------------------------------------------------------------------
// kgBuildTI2TINBList_gamd2_kernel:
//---------------------------------------------------------------------------------------------
_kPlainHead_
#ifdef use_DPFP
__launch_bounds__(128, 1)
#elif (__CUDA_ARCH__ >= 600)
__launch_bounds__(64, 16)
#else
__launch_bounds__(1024, 1)
#endif
kgBuildTI2TINBList_gamd2_kernel()
{
  PMEDouble* pX = cSim.pImageX;
  PMEDouble* pY = cSim.pImageY;
  PMEDouble* pZ = cSim.pImageZ;

  static const int MaxNumberTIAtom2 = 4000;

//  __shared__ volatile bool bondMask[(gti_simulationConst::MaxNumberTIAtom2 + 1) *
//                                    gti_simulationConst::MaxNumberTIAtom2];
//  __shared__ volatile float coordTI[gti_simulationConst::MaxNumberTIAtom2][3];

    __shared__ volatile bool bondMask[(MaxNumberTIAtom2 + 1) *
                                              MaxNumberTIAtom2];
    __shared__ volatile float coordTI[MaxNumberTIAtom2][3];


  unsigned i = threadIdx.x;
  __shared__ PMEFloat recip[9], ucell[9];

  if (cSim.pNTPData == NULL) {
    if (i < 3) {
      recip[i*3    ] = cSim.recipf[i][0];
      recip[i*3 + 1] = cSim.recipf[i][1];
      recip[i*3 + 2] = cSim.recipf[i][2];
      ucell[i*3    ] = cSim.ucellf[i][0];
      ucell[i*3 + 1] = cSim.ucellf[i][1];
      ucell[i*3 + 2] = cSim.ucellf[i][2];
    }
  }
  else {
    if (i < 9) {
      ucell[i] = cSim.pNTPData->ucellf[i];
      recip[i] = cSim.pNTPData->recipf[i];
    }
  }
  i = threadIdx.x;
  while (i<cSim.numberTIAtoms) {
    unsigned atom = cSim.pTIAtomList[i].x;
    unsigned iatom = cSim.pImageAtomLookup[atom];
    coordTI[i][0] = pX[iatom];
    coordTI[i][1] = pY[iatom];
    coordTI[i][2] = pZ[iatom];
    i += blockDim.x;
  }
  __syncthreads();

  unsigned totalPair = cSim.numberTIAtoms *( cSim.numberTIAtoms + 1 ) / 2;
  PMEFloat cc = (cSim.pNTPData == NULL) ? cSim.cutPlusSkin2 : cSim.pNTPData->cutPlusSkin2;
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  while (pos < totalPair) {
    unsigned index0 = (sqrt(pos*8+1.0001)-1)/2;
    unsigned index1 = pos-index0*(index0+1)/2;
    unsigned region0 = cSim.pTIAtomList[index0].y;
    unsigned region1 = cSim.pTIAtomList[index1].y;
    if (region0 == region1 && index0 != index1) { // must be in the same TI region
      PMEFloat dx = coordTI[index0][0] - coordTI[index1][0];
      PMEFloat dy = coordTI[index0][1] - coordTI[index1][1];
      PMEFloat dz = coordTI[index0][2] - coordTI[index1][2];
      PMEFloat r2 = __image_dist2<PMEFloat>(dx, dy, dz, recip, ucell);
      if (r2 < cc) {
        unsigned atom0 = cSim.pTIAtomList[index0].x;
        unsigned atom1 = cSim.pTIAtomList[index1].x;
        bool found = false;  //  Check if atom1 is on atom0's exclusion list
        unsigned start = cSim.pNLExclusionStartCount[atom0].x;
        unsigned count = cSim.pNLExclusionStartCount[atom0].y;
        for (unsigned k = start; k < start + count; k++) {
          if (atom1 == cSim.pNLExclusionList[k]) {
            found = true;
            break;
          }
        }
        unsigned int currentNLPos = atomicAdd(cSim.pNumberTINBEntries, 1);
        cSim.pTINBList[currentNLPos].x = atom0;
        cSim.pTINBList[currentNLPos].y = atom1;
        cSim.pTINBList[currentNLPos].z = region0;
//        printf("%d,%d,atom0,atom1 region\n",atom0,atom1,region0);
        if (found) {

          // This is a bonded or excluded pair
          cSim.pTINBList[currentNLPos].z |= 4;
        }
        else {
          int tt = cSim.TIVdwNTyes * ( cSim.pTIac[atom0] - 1 ) + ( cSim.pTIac[atom1] - 1 );
          int vdwindex = ( tt >= 0 ) ? cSim.pTIico[tt] - 1 : -1;
          if (vdwindex >= 0) {
            cSim.pTINBList[currentNLPos].w = vdwindex;
            cSim.pTINBList[currentNLPos].z |= 2;
          }
        }
        if (cSim.pSCList[atom0] > 0 || cSim.pSCList[atom1] > 0) {

          // Soft-core pair
          cSim.pTINBList[currentNLPos].z |= 24;
        }
        if (cSim.pSCList[atom0] > 0 && cSim.pSCList[atom1] > 0) {

          // Inter-region soft-core pair
          cSim.pTINBList[currentNLPos].z |= 32;
        }
      }
    }
    pos += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kgBuildTI2TINBList_gamd_kernel:
//---------------------------------------------------------------------------------------------
_kPlainHead_
#if !defined(AMBER_PLATFORM_AMD)
#ifdef use_DPFP
__LAUNCH_BOUNDS__(128, 1)
#elif (__CUDA_ARCH__ >= 600)
__LAUNCH_BOUNDS__(64, 16)
#else
__LAUNCH_BOUNDS__(1024, 1)
#endif
#endif
kgBuildTI2TINBList_gamd_kernel()
{
  PMEDouble* pX = cSim.pImageX;
  PMEDouble* pY = cSim.pImageY;
  PMEDouble* pZ = cSim.pImageZ;

  // __shared__ volatile bool bondMask[(cSim.atoms + 1) *
                                    // cSim.atoms];
  // __shared__ volatile float coordTI[cSim.atoms][3];
  // bool  *bondMask = NULL; // not used
  float *coordTIX = NULL;
  float *coordTIY = NULL;
  float *coordTIZ = NULL;

  unsigned i = threadIdx.x;
  __shared__ PMEFloat recip[9], ucell[9];

  if (cSim.pNTPData == NULL) {
    if (i < 3) {
      recip[i*3    ] = cSim.recipf[i][0];
      recip[i*3 + 1] = cSim.recipf[i][1];
      recip[i*3 + 2] = cSim.recipf[i][2];
      ucell[i*3    ] = cSim.ucellf[i][0];
      ucell[i*3 + 1] = cSim.ucellf[i][1];
      ucell[i*3 + 2] = cSim.ucellf[i][2];
    }
  }
  else {
    if (i < 9) {
      ucell[i] = cSim.pNTPData->ucellf[i];
      recip[i] = cSim.pNTPData->recipf[i];
    }
  }
  i = threadIdx.x;
  while (i<cSim.numberTIAtoms) {
    unsigned atom = cSim.pTIAtomList[i].x;
    unsigned iatom = cSim.pImageAtomLookup[atom];
    coordTIX[i] = pX[iatom];
    coordTIY[i] = pY[iatom];
    coordTIZ[i] = pZ[iatom];
    i += blockDim.x;
  }
  __syncthreads();

  unsigned totalPair = cSim.numberTIAtoms *( cSim.numberTIAtoms + 1 ) / 2;
  PMEFloat cc = (cSim.pNTPData == NULL) ? cSim.cutPlusSkin2 : cSim.pNTPData->cutPlusSkin2;
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  while (pos < totalPair) {
    unsigned index0 = (sqrt(pos*8+1.0001)-1)/2;
    unsigned index1 = pos-index0*(index0+1)/2;
    unsigned region0 = cSim.pTIAtomList[index0].y;
    unsigned region1 = cSim.pTIAtomList[index1].y;
    if (region0 == region1 && index0 != index1) { // must be in the same TI region
      PMEFloat dx = coordTIX[index0] - coordTIX[index1];
      PMEFloat dy = coordTIY[index0] - coordTIY[index1];
      PMEFloat dz = coordTIZ[index0] - coordTIZ[index1];
      PMEFloat r2 = __image_dist2<PMEFloat>(dx, dy, dz, recip, ucell);
      if (r2 < cc) {
        unsigned atom0 = cSim.pTIAtomList[index0].x;
        unsigned atom1 = cSim.pTIAtomList[index1].x;
        bool found = false;  //  Check if atom1 is on atom0's exclusion list
        unsigned start = cSim.pNLExclusionStartCount[atom0].x;
        unsigned count = cSim.pNLExclusionStartCount[atom0].y;
        for (unsigned k = start; k < start + count; k++) {
          if (atom1 == cSim.pNLExclusionList[k]) {
            found = true;
            break;
          }
        }

        unsigned int currentNLPos = atomicAdd(cSim.pNumberTINBEntries, 1);
        cSim.pTINBList[currentNLPos].x = atom0;
        cSim.pTINBList[currentNLPos].y = atom1;
        cSim.pTINBList[currentNLPos].z = region0;

        if (found) {

          // This is a bonded or excluded pair
          cSim.pTINBList[currentNLPos].z |= 4;
        }
        else {
          int tt = cSim.TIVdwNTyes * ( cSim.pTIac[atom0] - 1 ) + ( cSim.pTIac[atom1] - 1 );
          int vdwindex = ( tt >= 0 ) ? cSim.pTIico[tt] - 1 : -1;
          if (vdwindex >= 0) {
            cSim.pTINBList[currentNLPos].w = vdwindex;
            cSim.pTINBList[currentNLPos].z |= 2;
          }
        }
        if (cSim.pSCList[atom0] > 0 || cSim.pSCList[atom1] > 0) {

          // Soft-core pair
          cSim.pTINBList[currentNLPos].z |= 24;
        }
        if (cSim.pSCList[atom0] > 0 && cSim.pSCList[atom1] > 0) {

          // Inter-region soft-core pair
          cSim.pTINBList[currentNLPos].z |= 32;
        }
      }
    }
    pos += increment;
  }
}
#endif /* GTI */
