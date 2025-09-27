#ifndef AMBER_XRAY_COMMON_H
#define AMBER_XRAY_COMMON_H
#ifdef AMBER_PLATFORM_AMD
     #define THRUST_DEVICE_SYSTEM 5
     #define __HIP_PLATFORM_AMD__
#else
     #define __HIP_PLATFORM_NVIDIA__
#endif
#include <thrust/complex.h>

using complex_double = thrust::complex<double>;

namespace xray {
  enum class KernelPrecision {
    Single,
    Double,
  };

  template<KernelPrecision>
  struct KernelConfig;

  template<>
  struct KernelConfig<KernelPrecision::Single>{
    using FloatType = float;
  };

  template<>
  struct KernelConfig<KernelPrecision::Double>{
    using FloatType = double;
  };
}

#endif //AMBER_XRAY_COMMON_H
