#ifndef _CUDA_SIMULATION
#define _CUDA_SIMULATION

#include "base_simulationConst.h"
#ifdef GTI
#  include "gti_simulationConst.h"
typedef gti_simulationConst simulationConst;
#else
typedef base_simulationConst simulationConst;
#endif

#if defined(__CUDACC_RDC__) || defined(__HIPCC_RDC__)
#define CSIM_STO extern __device__ __constant__
#else
#define CSIM_STO __device__ __constant__
#endif

#endif /* _CUDA_SIMULATION */
