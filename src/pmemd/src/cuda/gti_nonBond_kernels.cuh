#ifndef GTI_NB_KERNELS
#  define GTI_NB_KERNELS
#  ifdef GTI

#    include "gpuContext.h"
#    include "gti_kernelAttribute.cuh"

#    include "gti_utils.cuh"
#    if !(defined(__CUDACC_RDC__) || defined(__HIPCC_RDC__))
  namespace GTI_NB_CALC {
#      include "gti_localCUutils.cuh"
  }
#    endif

// Real kernels
_kPlainHead_ kgCalculateTINB_kernel(bool needEnergy, bool needVirial);

_kPlainHead_ kgCalculateTI14NB_kernel(bool energy, bool virial);

_kPlainHead_ kgCalculateREAFNb_kernel(bool needEnergy, bool needVirial);

_kPlainHead_ kgCalculateREAF14NB_kernel(bool energy, bool virial);

_kPlainHead_ kgCalculate1264NB_kernel(bool energy, bool virial);

_kPlainHead_ kgCalculatep1264NB_kernel(bool energy, bool virial); // C4PairwiseCUDA

_kPlainHead_ kgCalculate1264p1264NB_kernel(bool energy, bool virial); // C4PairwiseCUDA

// gamd kernels
_kPlainHead_ kgCalculateTINB_gamd_kernel(bool needEnergy, bool needVirial);

_kPlainHead_ kgCalculateTI14NB_gamd_kernel(bool energy, bool virial);

_kPlainHead_ kgCalculate1264NB_gamd_kernel(bool energy, bool virial);

_kPlainHead_ kgCalculateTINB_ppi_gamd_kernel(bool needEnergy, bool needVirial);

_kPlainHead_ kgCalculateTINB_ppi2_gamd_kernel(bool energy, bool virial, int bgpro2atm, int edpro2atm);

_kPlainHead_ kgCalculateTI14NB_ppi_gamd_kernel(bool energy, bool virial);



#  endif /* GTI */
#endif /* GTI_NB_KERNELS */
