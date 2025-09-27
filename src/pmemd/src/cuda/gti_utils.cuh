#ifndef GTI_UTILS
#define GTI_UTILS
#ifdef GTI

#ifdef AMBER_PLATFORM_AMD
#include <hip/hip_runtime.h>
#else
#include <cuda.h>
#include <cuda_runtime.h>
#endif

template <class T>
inline __host__ __device__ bool __image_dist2_cut(T& dx, T& dy, T& dz, T& cut){
  return (__image_dist2(dx, dy, dz) < cut);
}

template <class T>
inline __host__ __device__ T __image_dist2(T& dx, T& dy, T& dz, T recip[], T ucell[])
{
  T x = recip[0]*dx + recip[3]*dy + recip[6]*dz;
  T y = recip[4]*dy + recip[7]*dz;
  T z = recip[8]*dz;

  x = x - rintf(x);
  y = y - rintf(y);
  z = z - rintf(z);

  dx = ucell[0]*x + ucell[1]*y + ucell[2]*z;
  dy = ucell[4]*y + ucell[5]*z;
  dz = ucell[8]*z;

  return dx*dx + dy*dy + dz*dz;
}

template <class T>
__forceinline__ __device__ bool __image_dist2_cut(T& dx, T& dy, T& dz, T& cut,
                                                  T recip[], T ucell[])
{
  return (__image_dist2(dx, dy, dz, recip, ucell) < cut );
}

inline __host__ __device__ float3 operator-(float3& a)
{
  return make_float3(-a.x, -a.y, -a.z);
}
inline __host__ __device__ void operator+=(float3& a, float3 b)
{
  a.x += b.x;
  a.y += b.y;
  a.z += b.z;
}

inline __host__ __device__ void operator-=(float3& a, float3 b)
{
  a.x -= b.x;
  a.y -= b.y;
  a.z -= b.z;
}

#ifndef AMBER_PLATFORM_AMD 
inline __host__ __device__ float3 operator*(float3 a, float3 b)
{
  return make_float3(a.x * b.x, a.y * b.y, a.z * b.z);
}

inline __host__ __device__ float3 operator*(float3 a, float b)
{
  return make_float3(a.x * b, a.y * b, a.z * b);
}

inline __host__ __device__ float3 operator*(float b, float3 a)
{
  return make_float3(b * a.x, b * a.y, b * a.z);
}

#endif

inline __host__ __device__ void operator*=(float3& a, float3 b)
{
  a.x *= b.x;
  a.y *= b.y;
  a.z *= b.z;
}
inline __host__ __device__ void operator*=(float3& a, float b)
{
  a.x *= b;
  a.y *= b;
  a.z *= b;
}
inline __host__ __device__ double3 operator-(double3& a)
{
  return make_double3(-a.x, -a.y, -a.z);
}

inline __host__ __device__ void operator+=(double3& a, double3 b)
{
  a.x += b.x;
  a.y += b.y;
  a.z += b.z;
}

inline __host__ __device__ void operator-=(double3& a, double3 b)
{
  a.x -= b.x;
  a.y -= b.y;
  a.z -= b.z;
}

inline __host__ __device__ void operator*=(double3& a, double3 b)
{
  a.x *= b.x;
  a.y *= b.y;
  a.z *= b.z;
}
inline __host__ __device__ void operator*=(double3& a, double b)
{
  a.x *= b;
  a.y *= b;
  a.z *= b;
}

#ifndef AMBER_PLATFORM_AMD
inline __host__ __device__ double3 operator*(double3 a, double3 b)
{
  return make_double3(a.x * b.x, a.y * b.y, a.z * b.z);
}

inline __host__ __device__ double3 operator*(double3 a, double b)
{
  return make_double3(a.x * b, a.y * b, a.z * b);
}

inline __host__ __device__ double3 operator*(double b, double3 a)
{
  return make_double3(b * a.x, b * a.y, b * a.z);
}
#endif

// GPU specific

#if defined(__CUDACC__) || defined(AMBER_PLATFORM_AMD)
#include "gputypes.h"
#include "ptxmacros.h"

#if defined(AMBER_PLATFORM_AMD)
// native implementation available
#elif (__CUDACC_VER_MAJOR__ >= 8) && (!defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600)
// native implementation available
#else
// Double precision atomicAdd, copied from CUDA_C_Programming_Guide.pdf (ver 5.0)
//
__device__ double atomicAdd(double* address, double val);
#endif

// Atomic inline accumulatings
static __forceinline__ __device__ void addEnergy(unsigned long long int* pE, PMEFloat energy)
{
  atomicAdd((unsigned long long int*) pE, ftoi(energy * eScale));
}

static __forceinline__ __device__ void addForce(PMEAccumulator* pF, PMEFloat force)
{
  atomicAdd((unsigned long long int*) pF, ftoi(force * fScale));
}

static __forceinline__ __device__ void addForce(unsigned long long int* pF, PMEFloat force)
{
  atomicAdd(pF, ftoi(force * fScale));
}

static __forceinline__ __device__ void addLocalForce(unsigned long long int* pF, PMEFloat force)
{
  (*pF) += ftoi(force * fScale);
}

static __forceinline__ __device__ void addLocalEnergy(unsigned long long int* pF, PMEFloat energy) {
  (*pF) += ftoi(energy * eScale);
}


static __forceinline__ __device__ void addVirial(unsigned long long int* pV, PMEFloat virial)
{
  atomicAdd((unsigned long long int*) pV, ftoi(virial * vScale));
}

static __forceinline__ __device__ void addForce(double* pF, PMEFloat force)
{
  atomicAdd(pF, force);
}

static __forceinline__ __device__ PMEFloat converter(unsigned long long int val, PMEFloat factor)
{
  PMEFloat result = (val >= 0x8000000000000000ull) ? -(PMEDouble)(val ^ 0xffffffffffffffffull) : (PMEDouble)val;

  return result * factor;
}

#if defined(__CUDACC_RDC__) || defined(__HIPCC_RDC__)
  /* declaration of global functions only visible when rdc=true*/
  void UpdateSimulationConst();

  /* declaration of already-defined global functions treated as extern functions when rdc=true*/
  extern __device__ void faster_sincos(PMEDouble a, PMEDouble *sptr, PMEDouble *cptr);
  extern __forceinline__ __device__ PMEFloat fasterfc(PMEFloat a);

#endif /* __CUDACC_RDC__ */

#endif

#endif  /* GTI */

#endif  /* GTI_UTILS */
