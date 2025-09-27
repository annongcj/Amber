#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

#ifndef PTXMACROS_H
#define PTXMACROS_H

#if defined(CUDA_VERSION) && (CUDA_VERSION >= 9000)
//#define __shfl(a,b)  __shfl_sync(0xFFFFFFFF,a,b)
#define __shfl_up(a,b)          __shfl_up_sync(0xFFFFFFFF,a,b)
#define __shfl_down(a,b)        __shfl_down_sync(0xFFFFFFFF,a,b)
#define __shfl_xor(a,b)         __shfl_xor_sync(0xFFFFFFFF,a,b)
#define __ballot(a)             __ballot_sync(0xFFFFFFFF,a)
#define __any(a)                __any_sync(0xFFFFFFFF,a)
#define __all(a)                __all_sync(0xFFFFFFFF,a)
#endif /* CUDA_VERSION */

#if (__CUDACC_VER_MAJOR__ >= 9)
#  define __SHFL_XOR(m,v,l)       __shfl_xor_sync(m,v,l)
#  define __SHFL_UP(m,v,d)        __shfl_up_sync(m,v,d)
#  define __SHFL_DOWN(m,v,d)      __shfl_down_sync(m,v,d)
#  define __SHFL(m,v,l)           __shfl_sync(m,v,l)
#  define __ALL(m,p)              __all_sync(m,p)
#  define __ANY(m,p)              __any_sync(m,p)
#  define __BALLOT(m,p)           __ballot_sync(m,p)
#  if (__CUDACC_VER_MAJOR__ >= 11)
#    define __SYNCWARP(m)         __syncwarp(m)
#  else
#    define __SYNCWARP(m)
#  endif
#else
#  define __SHFL_XOR(m,v,l)       __shfl_xor(v,l)
#  define __SHFL_UP(m,v,d)        __shfl_up(v,d)
#  define __SHFL_DOWN(m,v,d)      __shfl_down(v,d)
#  define __SHFL(m,v,l)           __shfl(v,l)
#  define __ALL(m,p)              __all(p)
#  define __ANY(m,p)              __any(p)
#  define __BALLOT(m,p)           __ballot(p)
#  ifdef AMBER_PLATFORM_AMD
#    define __SYNCWARP(m)         __syncthreads()
#  else
#    define __SYNCWARP(m)
#  endif
#endif

#ifdef AMBER_PLATFORM_AMD
template<class T, int dpp_ctrl, int row_mask = 0xf, int bank_mask = 0xf, bool bound_ctrl = true>
__device__ inline
T WarpMoveDpp(const T& input) {
    constexpr int words_no = (sizeof(T) + sizeof(int) - 1) / sizeof(int);

    struct V { int words[words_no]; };
    V a = __builtin_bit_cast(V, input);

    #pragma unroll
    for (int i = 0; i < words_no; i++) {
        a.words[i] = __builtin_amdgcn_update_dpp(
          0, a.words[i],
          dpp_ctrl, row_mask, bank_mask, bound_ctrl
        );
    }

    return __builtin_bit_cast(T, a);
}

#if defined(__gfx1030__) || defined(__gfx1031__)
// RDNA
#define SUPPORT_DPP_WAVEFRONT_SHIFTS 0
#else
// CDNA
#define SUPPORT_DPP_WAVEFRONT_SHIFTS 1
#endif

template<int Subwarp, class T>
__device__ inline
typename std::enable_if<(Subwarp == 16), T>::type
WarpRotateLeft(const T& input) {
  // Row rotate right by 15 lanes == left rotate by 1 lane
  return WarpMoveDpp<T, 0x12f>(input);
}

template<int Subwarp, class T>
__device__ inline
typename std::enable_if<!(Subwarp == 16) && !(SUPPORT_DPP_WAVEFRONT_SHIFTS && Subwarp == warpSize), T>::type
WarpRotateLeft(const T& input) {
  int i = ((threadIdx.x + 1) & (Subwarp - 1)) | (threadIdx.x & ~(Subwarp - 1));
  return __SHFL(WARP_MASK, input, i);
}

template<int Subwarp, class T>
__device__ inline
typename std::enable_if<(SUPPORT_DPP_WAVEFRONT_SHIFTS && Subwarp == warpSize), T>::type
WarpRotateLeft(const T& input) {
    // Wavefront rotate left by 1 thread (supported on CDNA)
    return WarpMoveDpp<T, 0x134>(input);
}
#endif

// HIP-TODO: Remove these if the proper overloads are introduced in HIP later
#ifdef AMBER_PLATFORM_AMD
__host__ __device__ inline float rsqrt(float x)
{
#ifdef __HIP_DEVICE_COMPILE__
    return __frsqrt_rn(x);
#else
    return rsqrtf(x);
#endif
}
#endif

// HIP-TODO: Remove these once HIP -ffast-math is fixed to trigger
// the right hardware instructions
#if defined(use_SPFP) && defined(AMBER_PLATFORM_AMD)
__device__ inline float fast_exp(float x)
{
    return __expf(x);
}

__host__ __device__ inline float flog(float x)
{
#ifdef __HIP_DEVICE_COMPILE__
    return __logf(x);
#else
    return logf(x);
#endif

}

__host__ __device__ inline double flog(double x)
{
    return log(x);
}

__host__ __device__ inline float fsqrt(const float x)
{
#ifdef __HIP_DEVICE_COMPILE__
    return __fsqrt_rn(x);
#else
    return sqrtf(x);
#endif
}

#ifndef AMBER_PLATFORM_AMD
__host__ __device__ inline double fsqrt(const double x)
{
    return sqrt(x);
}
#endif

#define exp fast_exp
#define sqrt fsqrt
#define log flog
#endif

//---------------------------------------------------------------------------------------------
// llitoulli: converts long long int --to--> unsigned long long int
//---------------------------------------------------------------------------------------------

__device__ inline unsigned long long int llitoulli(long long int l)
{
#ifdef AMBER_PLATFORM_AMD
  return (unsigned long long int) l;
#else
  unsigned long long int u;
  asm("mov.b64    %0, %1;" : "=l"(u) : "l"(l));
  return u;
#endif
}

__forceinline__
__device__ unsigned int maskPopc(unsigned int mask)
{
  return __popc(mask);
};

__forceinline__
__device__ unsigned int maskPopc(unsigned long long int mask)
{
  return __popcll(mask);
};

__forceinline__
__device__ unsigned int maskFfs(unsigned int mask)
{
  return __ffs(mask);
};

__forceinline__
__device__ unsigned int maskFfs(unsigned long long int mask)
{
  return __ffsll(mask);
};


//---------------------------------------------------------------------------------------------
// ullitolli: converts unsigned long long int --to--> long long int
//---------------------------------------------------------------------------------------------
__device__ inline long long int ullitolli(unsigned long long int u)
{
#ifdef AMBER_PLATFORM_AMD
  return (long long int) u;
#else
  long long int l;
  asm("mov.b64    %0, %1;" : "=l"(l) : "l"(u));
  return l;
#endif
}

//---------------------------------------------------------------------------------------------
// float2todouble: converts float2 (two floats, back-to-back) --to--> double
//---------------------------------------------------------------------------------------------
__device__ inline double float2todouble(float2 f)
{
  double d;
#ifdef AMBER_PLATFORM_AMD
  __builtin_memcpy(&d, &f, sizeof(double));
#else
  asm("mov.b64         %0, {%1, %2};" : "=d"(d) : "f"(f.x), "f"(f.y));
#endif
  return d;
}

//---------------------------------------------------------------------------------------------
// doubletofloat2: converts double --to--> float2 (two floats, back-to-back)
//---------------------------------------------------------------------------------------------
__device__ inline float2 doubletofloat2(double d)
{
    float2 f;
#ifdef AMBER_PLATFORM_AMD
  __builtin_memcpy(&f, &d, sizeof(double));
#else
    asm("mov.b64    {%0, %1}, %2;" : "=f"(f.x), "=f"(f.y) : "d"(d));
#endif
    return f;
}

#ifndef AMBER_PLATFORM_AMD
//---------------------------------------------------------------------------------------------
// reduce: not sure what this does.
//---------------------------------------------------------------------------------------------
__device__ inline unsigned int reduce(unsigned int x)
{
  unsigned int s = x;
  asm("shfl.up %0, %1, 0x1, 0x0;\n\t"
   "   add.u32 %1, %0, %1;\n\t"
      "shfl.up %0, %1, 0x2, 0x0;\n\t"
   "   add.u32 %1, %0, %1;\n\t"
      "shfl.up %0, %1, 0x4, 0x0;\n\t"
   "   add.u32 %1, %0, %1;\n\t"
      "shfl.up %0, %1, 0x8, 0x0;\n\t"
   "   add.u32 %1, %0, %1;\n\t"
      "shfl.up %0, %1, 0x16, 0x0;\n\t"
   "   add.u32 %1, %0, %1;"
   : "=r"(s) : "r"(x));
  return s;
}
#if CUDA_VERSION < 6050

//---------------------------------------------------------------------------------------------
// __shfl: this is ALL over the place, but it is overloaded with different POD types for more
//         utility
//---------------------------------------------------------------------------------------------
static __device__ __inline__ unsigned int __shfl(unsigned int var, unsigned int srcLane,
                                                 int width=32)
{
  unsigned int ret, c;
  c = ((32-width) << 8) | 0x1f;
  asm volatile ("shfl.idx.b32 %0, %1, %2, %3;" : "=r"(ret) : "r"(var), "r"(srcLane), "r"(c));

  return ret;
}

//---------------------------------------------------------------------------------------------
static __device__ __inline__ unsigned int __shfl(unsigned int var, int srcLane, int width=32)
{
  unsigned int ret, c;

  c = ((32-width) << 8) | 0x1f;
  asm volatile ("shfl.idx.b32 %0, %1, %2, %3;" : "=r"(ret) : "r"(var), "r"(srcLane), "r"(c));

  return ret;
}

//---------------------------------------------------------------------------------------------
static __device__ __inline__ float __shfl(float var, unsigned int srcLane, int width=32)
{
  float ret;
  unsigned int c;

  c = ((32-width) << 8) | 0x1f;
  asm volatile ("shfl.idx.b32 %0, %1, %2, %3;" : "=f"(ret) : "f"(var), "r"(srcLane), "r"(c));

  return ret;
}

//---------------------------------------------------------------------------------------------
static __device__ __inline__ double __shfl(double var, unsigned int srcLane, int width=32)
{
  double ret;
  unsigned int c;
  c = ((32-width) << 8) | 0x1f;
  asm volatile ("{\n\t"
                ".reg .u32 dlo;\n\t"
                ".reg .u32 dhi;\n\t"
                " mov.b64 {dlo, dhi}, %1;\n\t"
                " shfl.idx.b32 dlo, dlo, %2, %3;\n\t"
                " shfl.idx.b32 dhi, dhi, %2, %3;\n\t"
                " mov.b64 %0, {dlo, dhi};\n\t"
                "}"
                : "=d"(ret) : "d"(var), "r"(srcLane), "r"(c));
  return ret;
}

//---------------------------------------------------------------------------------------------
static __device__ __inline__ double __shfl(double var, int srcLane, int width=32)
{
  double ret;
  unsigned int c;
  c = ((32-width) << 8) | 0x1f;
  asm volatile ("{\n\t"
                ".reg .u32 dlo;\n\t"
                ".reg .u32 dhi;\n\t"
                " mov.b64 {dlo, dhi}, %1;\n\t"
                " shfl.idx.b32 dlo, dlo, %2, %3;\n\t"
                " shfl.idx.b32 dhi, dhi, %2, %3;\n\t"
                " mov.b64 %0, {dlo, dhi};\n\t"
                "}"
                : "=d"(ret) : "d"(var), "r"(srcLane), "r"(c));
  return ret;
}

//---------------------------------------------------------------------------------------------
static __device__ __inline__ long long int __shfl(long long int var, int srcLane, int width=32)
{
  long long int ret;
  unsigned int c;

  c = ((32-width) << 8) | 0x1f;
  asm volatile ("{\n\t"
                ".reg .u32 dlo;\n\t"
                ".reg .u32 dhi;\n\t"
                " mov.b64 {dlo, dhi}, %1;\n\t"
                " shfl.idx.b32 dlo, dlo, %2, %3;\n\t"
                " shfl.idx.b32 dhi, dhi, %2, %3;\n\t"
                " mov.b64 %0, {dlo, dhi};\n\t"
                "}"
                : "=l"(ret) : "l"(var), "r"(srcLane), "r"(c));

  return ret;
}

//---------------------------------------------------------------------------------------------
static __device__ __inline__ long long int __shfl(long long int var, unsigned int srcLane,
                                                  int width=32)
{
  long long int ret;
  unsigned int c;

  c = ((32-width) << 8) | 0x1f;
  asm volatile ("{\n\t"
                ".reg .u32 dlo;\n\t"
                ".reg .u32 dhi;\n\t"
                " mov.b64 {dlo, dhi}, %1;\n\t"
                " shfl.idx.b32 dlo, dlo, %2, %3;\n\t"
                " shfl.idx.b32 dhi, dhi, %2, %3;\n\t"
                " mov.b64 %0, {dlo, dhi};\n\t"
                "}"
                : "=l"(ret) : "l"(var), "r"(srcLane), "r"(c));

  return ret;
}

//---------------------------------------------------------------------------------------------
static __device__ __inline__ unsigned long long int __shfl(unsigned long long int var,
                                                           int srcLane, int width=32)
{
  unsigned long long int ret;
  unsigned int c;
  c = ((32-width) << 8) | 0x1f;
  asm volatile ("{\n\t"
                ".reg .u32 dlo;\n\t"
                ".reg .u32 dhi;\n\t"
                " mov.b64 {dlo, dhi}, %1;\n\t"
                " shfl.idx.b32 dlo, dlo, %2, %3;\n\t"
                " shfl.idx.b32 dhi, dhi, %2, %3;\n\t"
                " mov.b64 %0, {dlo, dhi};\n\t"
                "}"
                : "=l"(ret) : "l"(var), "r"(srcLane), "r"(c));
  return ret;
}
#endif

//---------------------------------------------------------------------------------------------
// __add64: this inline function for 64-bit integer addition does not appear to be called
//          anywhere in the code.
//---------------------------------------------------------------------------------------------
static __device__ __inline__ long long int __add64(long long int a, long long int b)
{
  long long int ret;
  asm volatile ("{\n\t"
                ".reg .s32 alo;\n\t"
                ".reg .s32 ahi;\n\t"
                ".reg .s32 blo;\n\t"
                ".reg .s32 bhi;\n\t"
                " mov.b64 {alo, ahi}, %1;\n\t"
                " mov.b64 {blo, bhi}, %2;\n\t"
                " add.cc.s32 alo, alo, blo;\n\t"
                " addc.cc.s32 ahi, ahi, bhi;\n\t"
                " mov.b64 %0, {alo, ahi};\n\t"
                "}"
                : "=l"(ret) : "l"(a), "l"(b));

  return ret;
}
#endif

//---------------------------------------------------------------------------------------------
// This special version of llrintf is designed to provide better performance on Maxwell based
// cards.
//
// float -> int64 conversion
// 64-bit conversion operators run at reduced rate on Maxwell (4 threads per clock).
// This has a significant impact on AMBER performance.
//
// fast_llrintf provides a solution courtesy of Kate Clark:
//
// Inspired by Kahan summation, the approach breaks the conversion into two 32-bit conversion
// operations.  These operate at 32 threads per clock.  Conversion throughput is more than
// doubled.
//
// Publication: SPXP and other tricks - Walker, R.C., Le Grand, S., et al. 2016, in prep.
//
//   1) Requires CUDA >= 7.5
//
//   2) It looks like NVIDIA attempted to add something automatic along
//      these lines to CUDA 8.0 between RC and release version since having
//      Maxwell chips call regular llrintf with CUDA 8.0 has a lot less impact
//      on performance than it did with CUDA 7.5. But at the same time there was
//      a regression for Pascal chips on llrintf. So for now we will force
//      fastllrintf on Maxwell and 610 Pascal. Not optimal on GP100 but I
//      don't have access to the hardware to check such things so we force
//      llrintf for safety.
//---------------------------------------------------------------------------------------------
static __device__ __inline__ long long fast_llrintf(float x) {
#if defined(_WIN32) 
  return llrintf(x);
#else
  //  Maxwell hardware and GP102,104,107
  //  This solution also benefits AMD GPUs; this is faster than the default llrintf
  float z = x * (float)0x1.00000p-32;
  int hi = __float2int_rz( z );                         // First convert high bits
  float delta = x - ((float)0x1.00000p32*((float)hi));  // Check remainder sign
  int test = (__float_as_uint(delta) > 0xbf000000);
#  if defined(AMBER_PLATFORM_AMD)
  // HIP-TODO: implementation of __float2uint_rn is incorrect. Remove this when it's fixed.
  int lo = (unsigned int)rintf(fabsf(delta));           // Convert the (unsigned) remainder
#  else
  int lo = __float2uint_rn(fabsf(delta));               // Convert the (unsigned) remainder
#  endif
  lo = (test) ? -lo: lo;
  hi -= test;                                           // Two's complement correction
  long long res = __double_as_longlong(__hiloint2double(hi,lo)); // Return 64-bit result
  return res;
#endif
}

//---------------------------------------------------------------------------------------------
// TISetLambda: set the lambda value for TI functionality.  The lambda value changes based on
//              the TI region that each interaction resides in.
//
// Arguments:
//   CVterm:    flag to indicate that the interaction is within the CV region
//   TIregion:  identifier of the TI region to which the object belongs
//   lambda:    base value of the mixing factor
//---------------------------------------------------------------------------------------------
static __device__ __inline__ double TISetLambda(int CVterm, int TIregion, double lambda)
{
  if (CVterm) {
    if (TIregion == 1) {
      lambda = (double)1.0 - lambda;
    }
  }
  else {
    lambda = (double)1.0;
  }

  return lambda;
}

//---------------------------------------------------------------------------------------------
// TISetLambda: set the lambda value for TI functionality.  The lambda value changes based on
//              the TI region that each interaction resides in.
//
// Arguments:
//   CVterm:    flag to indicate that the interaction is within the CV region
//   TIregion:  identifier of the TI region to which the object belongs
//   lambda:    base value of the mixing factor
//---------------------------------------------------------------------------------------------
static __device__ __inline__ float TISetLambda(int CVterm, int TIregion, float lambda)
{
  if (CVterm) {
    if (TIregion == 1) {
      lambda = (float)1.0 - lambda;
    }
  }
  else {
    lambda = (float)1.0;
  }

  return lambda;
}

//TODO(dominic): Assumptions on warp size
//---------------------------------------------------------------------------------------------
// It can be useful to have prefix sums computed over warps, accumulated by __shfl instructions
// for maximum efficiency.  These macros handle inclusive and exclusive sums:
//
// Value:      1  2  4  3  1  5  7  4
// Inclusive:  1  3  7 10 11 16 23 27
// Exclusive:  0  1  3  7 10 11 16 23
//
// Arguments:
//   var:     name of the thread-local variable over which to take the prefix sum.  For
//            double prefix sums, two variables, varX, and varY, are formed into prefix sums
//            side-by-side.
//   tgx:     thread index within the warp
//   result:  Value to hold the total, inclusive sum of var over all threads
//            (exclusive sum with only)
//---------------------------------------------------------------------------------------------
#define inclusive_warp_prefixsum(var, tgx) \
{ \
  var += ((tgx &  1) ==  1) * __shfl_up(var, 1); \
  var += ((tgx &  3) ==  3) * __shfl_up(var, 2); \
  var += ((tgx &  7) ==  7) * __shfl_up(var, 4); \
  var += ((tgx & 15) == 15) * __shfl_up(var, 8); \
  var += (tgx == 31) * __shfl_up(var, 16); \
  var += ((tgx & 15) == 7 && tgx > 16) * __shfl_up(var, 8); \
  var += ((tgx &  7) == 3 && tgx > 8)  * __shfl_up(var, 4); \
  var += ((tgx &  3) == 1 && tgx > 4)  * __shfl_up(var, 2); \
  var += ((tgx &  1) == 0 && tgx >= 2) * __shfl_up(var, 1); \
}

#define inclusive_warp_double_prefixsum(varX, varY, tgx) \
{ \
  int mval = ((tgx & 1) == 1); \
  varX += mval * __shfl_up(varX, 1); \
  varY += mval * __shfl_up(varY, 1); \
  mval = ((tgx & 3) == 3); \
  varX += mval * __shfl_up(varX, 2); \
  varY += mval * __shfl_up(varY, 2); \
  mval = ((tgx & 7) == 7); \
  varX += mval * __shfl_up(varX, 4); \
  varY += mval * __shfl_up(varY, 4); \
  mval = ((tgx & 15) == 15); \
  varX += mval * __shfl_up(varX, 8); \
  varY += mval * __shfl_up(varY, 8); \
  mval = (tgx == 31); \
  varX += mval * __shfl_up(varX, 16); \
  varY += mval * __shfl_up(varY, 16); \
  mval = ((tgx & 15) == 7 && tgx > 16); \
  varX += mval * __shfl_up(varX, 8); \
  varY += mval * __shfl_up(varY, 8); \
  mval = ((tgx & 7) == 3 && tgx > 8); \
  varX += mval * __shfl_up(varX, 4); \
  varY += mval * __shfl_up(varY, 4); \
  mval = ((tgx & 3) == 1 && tgx > 4); \
  varX += mval * __shfl_up(varX, 2); \
  varY += mval * __shfl_up(varY, 2); \
  mval = ((tgx & 1) == 0 && tgx >= 2); \
  varX += mval * __shfl_up(varX, 1); \
  varY += mval * __shfl_up(varY, 1); \
}

#define exclusive_warp_prefixsum(var, tgx) \
{ \
  var += ((tgx &  1) ==  1) * __shfl_up(var, 1); \
  var += ((tgx &  3) ==  3) * __shfl_up(var, 2); \
  var += ((tgx &  7) ==  7) * __shfl_up(var, 4); \
  var += ((tgx & 15) == 15) * __shfl_up(var, 8); \
  var += (tgx == 31) * __shfl_up(var, 16); \
  var += ((tgx & 15) == 7 && tgx > 16) * __shfl_up(var, 8); \
  var += ((tgx &  7) == 3 && tgx > 8)  * __shfl_up(var, 4); \
  var += ((tgx &  3) == 1 && tgx > 4)  * __shfl_up(var, 2); \
  var += ((tgx &  1) == 0 && tgx >= 2) * __shfl_up(var, 1); \
  var = __shfl_up(var, 1); \
  if (tgx == 0) { \
    var = 0; \
  } \
}

#define exclusive_warp_prefixsum_savetotal(var, tgx, result) \
{ \
  var += ((tgx &  1) ==  1) * __shfl_up(var, 1); \
  var += ((tgx &  3) ==  3) * __shfl_up(var, 2); \
  var += ((tgx &  7) ==  7) * __shfl_up(var, 4); \
  var += ((tgx & 15) == 15) * __shfl_up(var, 8); \
  var += (tgx == 31) * __shfl_up(var, 16); \
  var += ((tgx & 15) == 7 && tgx > 16) * __shfl_up(var, 8); \
  var += ((tgx &  7) == 3 && tgx > 8)  * __shfl_up(var, 4); \
  var += ((tgx &  3) == 1 && tgx > 4)  * __shfl_up(var, 2); \
  var += ((tgx &  1) == 0 && tgx >= 2) * __shfl_up(var, 1); \
  if (tgx == 31) { \
    result = var; \
  } \
  var = __shfl_up(var, 1); \
  if (tgx == 0) { \
    var = 0; \
  } \
}

#define warp_reduction_down(var) \
{ \
  var += __SHFL_DOWN(0xffffffff, var, 16); \
  var += __SHFL_DOWN(0xffffffff, var,  8); \
  var += __SHFL_DOWN(0xffffffff, var,  4); \
  var += __SHFL_DOWN(0xffffffff, var,  2); \
  var += __SHFL_DOWN(0xffffffff, var,  1); \
}

#endif
