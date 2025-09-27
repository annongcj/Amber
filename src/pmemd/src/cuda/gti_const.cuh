#ifndef GTI_CONST
#define GTI_CONST
#ifdef GTI
#include "gti_def.h"

#ifdef GL_CONST
#undef GL_CONST
#endif

#define XSTR(x) STR(x)
#define STR(x) #x

#if defined(__CUDACC__) || defined(AMBER_PLATFORM_AMD)
#  if defined(AMBER_PLATFORM_AMD)
#    include <hip/hip_runtime.h>
#  endif
#  define GL_CONST  static const __constant__
#  define __TARGET __device__ inline
#else
#define GL_CONST  static const
#define __TARGET inline
#define Fma(a,b,c) (a*b+c)
#endif

#ifdef use_SPFP
const unsigned int MaxLJIndex = 1500;
const unsigned int MaxPLJIndex = 100; //C4PairwiseCUDA2023
#else 
const unsigned int MaxLJIndex = 700; 
const unsigned int MaxPLJIndex = 50; //C4PairwiseCUDA2023
#endif

GL_CONST unsigned Zero = 0;
GL_CONST unsigned One = 1;

GL_CONST Real ZeroF = 0.0;
GL_CONST Real QuarterF = 0.25;
GL_CONST Real HalfF = 0.5;
GL_CONST Real OneF = 1.0;
GL_CONST Real TwoF = 2.0;
GL_CONST Real ThreeF = 3.0;
GL_CONST Real FourF = 4.0;
GL_CONST Real SixF = 6.0;
GL_CONST Real TenF = 10.0;
GL_CONST Real TwelveF = 12.0;
GL_CONST Real FifteenF = 15.0;

#endif

#endif
