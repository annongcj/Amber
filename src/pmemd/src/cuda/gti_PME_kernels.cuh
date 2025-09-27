#ifndef GTI_PME_KERNELS
#  define GTI_PME_KERNELS
#  ifdef GTI

#    include "gti_def.h"
#    include "gpuContext.h"
#    include "gti_kernelAttribute.cuh"

#    include "gti_utils.cuh"
#    if !(defined(__CUDACC_RDC__) || defined(__HIPCC_RDC__))
  namespace GTI_PME_IMPL {
#      include "gti_localCUutils.cuh"
  }
#    endif

// Routines need to be on the device side code

// Real kernels
_kPlainHead_ kgTIPMEFillChargeGrid_kernel();

_kReduceFrcHead_ kgTIPMEReduceChargeGrid_kernel();

namespace GTI_PME_IMPL {
  struct FillChargeGridAtomData {
    int ix;
    int iy;
    int iz;
    PMEFloat tx[4];
    PMEFloat ty[4];
    PMEFloat tz[4];
  };

#    ifdef use_DPFP
    __constant__ const unsigned int LOADSIZE = 2;
    __constant__ const unsigned int FillY = 8;
#    else
#      if (__CUDA_ARCH__ >= 600)
      __constant__ const unsigned LOADSIZE = 8;
#      else
      __constant__ const unsigned LOADSIZE = 16;
#      endif /*  __CUDA_ARCH__  */
    __constant__ const unsigned int FillY = 16;
#    endif /* use_DPFP  */
  static __constant__ const PMEFloat ONEHALF = 0.5;
  static __constant__ const PMEFloat ONETHIRD = (1.0 / 3.0);
  static __constant__ const PMEFloat ONE = 1.0;
  static __constant__ const PMEFloat TWO = 2.0;
  static __constant__ const PMEFloat THREE = 3.0;
}

#  endif /* GTI */
#endif /* GTI_PME_KERNELS */
