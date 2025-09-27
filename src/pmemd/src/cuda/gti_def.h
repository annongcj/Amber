
#ifndef __GTI_DEF_H__
#define __GTI_DEF_H__
#  ifdef GTI

#include "hip_definitions.h"

#define AMPERE (gpu->major == 8)
#define SM80 ( (gpu->major == 8) && (gpu->minor == 0) )

#define VOLTA (gpu->major == 7)
#define SM70 ( (gpu->major == 7) && (gpu->minor == 0) )

#define PASCAL (gpu->major == 6)
#define SM60 ( (gpu->major == 6) && (gpu->minor == 0) )
#define SM61 ( (gpu->major == 6) && (gpu->minor == 1) )
#define SM62 ( (gpu->major == 6) && (gpu->minor == 2) )

#define MAXWELL (gpu->major == 5)
#define SM50 ( (gpu->major == 5) && (gpu->minor == 0) )
#define SM52 ( (gpu->major == 5) && (gpu->minor == 2) )
#define SM53 ( (gpu->major == 5) && (gpu->minor == 3) )

#define KEPLER (gpu->major == 3)
#define SM30 ( (gpu->major == 3) && (gpu->minor == 0) )
#define SM35 ( (gpu->major == 3) && (gpu->minor == 5) )
#define SM37 ( (gpu->major == 3) && (gpu->minor == 7) )

//################################
#ifdef use_SPFP

#define Real float

#define Erf erff
#define Erfc erfcf

#define Exp expf
#define Sqrt sqrtf
#define Cbrt cbrtf
#define Rsqrt rsqrtf
#define Rcbrt rcbrtf

#define Pow powf
#define Fma fmaf
#define Tanh tanhf

#else

#define Real double

#define Erf erf
#define Erfc erfc
#define Exp exp
#define Sqrt sqrt
#define Cbrt cbrt
#define Rsqrt rsqrt
#define Rcbrt rcbrt
#define Pow pow
#define Fma fma
#define Tanh tanh

#endif

#define isDPFP ( sizeof(Real)==8)
#define isSPFP ( sizeof(Real)==4)

#undef CUDA_MEM_CLEAR

#undef TI_DP_ACC

#define SYNC_MASS

//################################

#define CLEAN(pointer) if (pointer) delete pointer

//################################
#endif  /* GTI */
#endif  /*  __GTI_DEF_H__ */

