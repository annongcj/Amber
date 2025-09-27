#ifndef GTI_KERNELS
#define GTI_KERNELS

#  ifdef GTI

#    include "gpuContext.h"
#    include "gti_kernelAttribute.cuh"

#    include "gti_utils.cuh"
#    if !(defined(__CUDACC_RDC__) || defined(__HIPCC_RDC__))
  namespace GTI_GENERAL_IMPL {
#      include "gti_localCUutils.cuh"
  }
#    endif

// Routines need to be on the device side code

// Real kernels
_kPlainHead_ kgClearTIForce_kernel();

_kPlainHead_ kgClearTIForce_gamd_kernel();

_kPlainHead_ kgClearTIKinEnergy_kernel();

_kPlainHead_ kgClearTIPotEnergy_kernel();

_kPlainHead_ kgClearMBAR_kernel();

_kPlainHead_ kgZeroTICharges_kernel(bool useImage, int mode);

_kPlainHead_ kgScaleRECharges_kernel(bool useImage, int mode);

_kPlainHead_ kgZeroTIAtomForce_kernel(bool useImage, int mode);

_kPlainHead_ kgCopyToTIForce(int TIRegion, bool isNB, bool keepSource, PMEFloat weight);

_kPlainHead_ kgCopyToTIEnergy_kernel(unsigned int TIRegion, int term1, int term2, int term3,
                                     bool isVirial, PMEFloat weight, bool addOn);

_kPlainHead_ kgCorrectTIEnergy_kernel(int beginTerm, int endTerm);

_kPlainHead_ kgCorrectTIForce_kernel();

_kPlainHead_ kgCombineTIForce_kernel(bool isNB, bool useImage, bool addOn, bool needWeight);

_kPlainHead_ kgCombineSCForce_kernel(bool isNB, bool useImage);

_kPlainHead_ kgCombineTIForce_gamd_kernel(bool isNB, bool useImage, bool addOn);

_kPlainHead_ kgCombineSCForce_gamd_kernel(bool isNB, bool useImage);

_kPlainHead_ kgCalculateNetForce_kernel(int TIRegion);

_kPlainHead_ kgRemoveNetForce_kernel(int TIRegion);

_kPlainHead_ kgSyncVector_kernel(bool useImage, unsigned mode, int combinedMode);
#  endif  /* GTI */

#endif  /* GTI_KERNELS */
