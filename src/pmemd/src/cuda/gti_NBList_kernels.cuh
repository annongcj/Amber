#ifndef GTI_NBList_KERNELS
#  define GTI_NBList_KERNELS
#  ifdef GTI

#    include "gpuContext.h"
#    include "gti_kernelAttribute.cuh"

// Real kernels
namespace GTI_NB {
  enum specialType {TI = 0, LJ1264 = 1, pLJ1264 = 3, LJ1264pLJ1264 = 4, REAF = 2 }; //C4PairwiseCUDA2023
}

#    include "gti_utils.cuh"
#    if !(defined(__CUDACC_RDC__) || defined(__HIPCC_RDC__))
  namespace GTI_NB_LIST {
#      include "gti_localCUutils.cuh"
  }
#    endif

_kPlainHead_ kgBuildSpecial2RestNBPreList_kernel(GTI_NB::specialType type);

_kPlainHead_ kgBuildSpecial2RestNBList_kernel(GTI_NB::specialType type);

_kPlainHead_ kgTINBListFillAttribute_kernel();

_kPlainHead_ kgREAFNbListFillAttribute_kernel();

_kPlainHead_ kg1264NBListFillAttribute_kernel();

_kPlainHead_ kgp1264NBListFillAttribute_kernel(); // C4PairwiseCUDA

_kPlainHead_ kg1264p1264NBListFillAttribute_kernel(); // C4PairwiseCUDA2023

/////////////////////////////////////////////////////////////////////////////
// TL:The following kernels are out-of-date and just kept for gamd

_kPlainHead_ kgBuildTI2RestNBListFillAttribute_kernel();

_kPlainHead_ kgBuild1264NBListFillAttribute_kernel();

_kPlainHead_ kgBuildp1264NBListFillAttribute_kernel(); //C4PairwiseCUDA2023

_kPlainHead_ kgBuild1264p1264NBListFillAttribute_kernel(); //C4PairwiseCUDA2023

_kPlainHead_ kgBuildTI2TINBList_kernel();

_kPlainHead_ kgBuildSpecial2RestNBPreList_gamd_kernel(GTI_NB::specialType type);

_kPlainHead_ kgBuildSpecial2RestNBList_gamd_kernel(GTI_NB::specialType type);

_kPlainHead_ kgBuildTI2TINBList_gamd_kernel();

_kPlainHead_ kgBuildTI2TINBList_gamd2_kernel();

#  endif /* GTI */
#endif /* GTI_NBList_KERNELS */
