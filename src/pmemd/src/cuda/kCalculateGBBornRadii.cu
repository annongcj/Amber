#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

#ifndef AMBER_PLATFORM_AMD
#include <cuda.h>
#endif
#include "gpu.h"
#include "ptxmacros.h"

// Use global instance instead of a local copy
#include "simulationConst.h"
CSIM_STO simulationConst cSim;

static __constant__ PMEFloat ta  = (PMEFloat)(1.0 / 3.0);
static __constant__ PMEFloat tb  = (PMEFloat)(2.0 / 5.0);
static __constant__ PMEFloat tc  = (PMEFloat)(3.0 / 7.0);
static __constant__ PMEFloat td  = (PMEFloat)(4.0 / 9.0);
static __constant__ PMEFloat tdd = (PMEFloat)(5.0 / 11.0);

#if !defined(__HIPCC_RDC__)

//---------------------------------------------------------------------------------------------
// SetkCalculateGBBornRadiiSim: establish parameters, notably cSim, on the GPU.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void SetkCalculateGBBornRadiiSim(gpuContext gpu)
{
  cudaError_t status;
  status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(simulationConst));
  RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

//---------------------------------------------------------------------------------------------
// GetkCalculateGBBornRadiiSim: download parameters from the GPU.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
void GetkCalculateGBBornRadiiSim(gpuContext gpu)
{
  cudaError_t status;
  status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(simulationConst));
  RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

#endif

//---------------------------------------------------------------------------------------------
// kCalculateGBBornRadii_kernel: general-purpose GB radii calculations
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GBBORNRADII_THREADS_PER_BLOCK, GBBORNRADII_BLOCKS_MULTIPLIER)
kCalculateGBBornRadii_kernel()
#include "kCalculateGBBornRadii.h"

//---------------------------------------------------------------------------------------------
// kCalculateGBBornRadiiIGB78_kernel: special cases for famous GB flavors
//---------------------------------------------------------------------------------------------
#define GB_IGB78
__global__ void
__LAUNCH_BOUNDS__(GBBORNRADII_THREADS_PER_BLOCK, GBBORNRADII_BLOCKS_MULTIPLIER)
kCalculateGBBornRadiiIGB78_kernel()
#include "kCalculateGBBornRadii.h"
#undef GB_IGB78

//---------------------------------------------------------------------------------------------
// kCalculateGBBornRadiiInitKernels: initialize kernels for GB radii calculations
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kCalculateGBBornRadiiInitKernels(gpuContext gpu)
{
  cudaFuncSetSharedMemConfig(kCalculateGBBornRadii_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalculateGBBornRadiiIGB78_kernel,
                             cudaSharedMemBankSizeEightByte);
    cudaFuncSetCacheConfig(kCalculateGBBornRadii_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kCalculateGBBornRadiiIGB78_kernel, cudaFuncCachePreferL1);
}

//---------------------------------------------------------------------------------------------
// kCalculateGBBornRadii: calculate the GB radii
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kCalculateGBBornRadii(gpuContext gpu)
{
  if ((gpu->sim.igb == 7) || (gpu->sim.igb == 8)) {
    kCalculateGBBornRadiiIGB78_kernel<<<gpu->GBBornRadiiBlocks,
                                        gpu->GBBornRadiiIGB78ThreadsPerBlock>>>();
  }
  else {
    kCalculateGBBornRadii_kernel<<<gpu->GBBornRadiiBlocks,
                                   gpu->GBBornRadiiThreadsPerBlock>>>();
  }
  LAUNCHERROR("kCalculateGBBornRadii");
}

//---------------------------------------------------------------------------------------------
// kReduceGBBornRadii_kernel: kernel for reducing GB radii computed across threads, called by
//                            the host function given directly below.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(REDUCEBUFFER_THREADS_PER_BLOCK, 1)
kReduceGBBornRadii_kernel()
{
  bool bIGB2578 = (cSim.igb == 2) || (cSim.igb == 5) || (cSim.igb == 7) || (cSim.igb == 8);
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  while (pos < cSim.atoms) {
    PMEDouble rborn_i = cSim.pAtomRBorn[pos];
    PMEDouble reff_i  = (PMEDouble)cSim.pReffAccumulator[pos] * ONEOVERFORCESCALE;

    // Process final Born Radii
    PMEDouble ri   = rborn_i - cSim.offset;
    PMEDouble ri1i = (PMEDouble)1.0 / ri;
    if (bIGB2578) {

      // Apply the new Onufriev "gbalpha, gbbeta, gbgamma" correction:
      PMEDouble psi_i = -ri * reff_i;
      if (cSim.igb == 8) {
        reff_i = ri1i - tanh((cSim.pgb_alpha[pos] + cSim.pgb_gamma[pos]*psi_i*psi_i -
                              cSim.pgb_beta[pos]*psi_i)*psi_i) / rborn_i;
      }
      else {
        reff_i = ri1i - tanh((cSim.gb_alpha + cSim.gb_gamma*psi_i*psi_i -
                              cSim.gb_beta*psi_i)*psi_i) / rborn_i;
      }
      reff_i = max(reff_i, (PMEDouble)1.0 / (PMEDouble)30.0);
      reff_i = (PMEDouble)1.0 / reff_i;
      cSim.pPsi[pos] = psi_i;
    }
    else {

      // "standard" GB, including the "diagonal" term here:
      reff_i = (PMEDouble)1.0 / (reff_i + ri1i);
    }

    cSim.pReffSP[pos] = reff_i;
    cSim.pReff[pos]   = reff_i;
    pos += blockDim.x * gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kReduceGBBornRadii: reduce the GB radii as pieces of them were calculate on many different
//                     threads.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kReduceGBBornRadii(gpuContext gpu)
{
  kReduceGBBornRadii_kernel<<<gpu->blocks, gpu->reduceBufferThreadsPerBlock>>>();
  LAUNCHERROR("kReduceGBBornRadii");
}

//---------------------------------------------------------------------------------------------
// kClearGBBuffers_kernel: clear buffers for GB computations
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(REDUCEBUFFER_THREADS_PER_BLOCK, 1)
kClearGBBuffers_kernel()
{
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  while (pos < cSim.atoms) {
    cSim.pReffAccumulator[pos]      = (PMEAccumulator)0;
    cSim.pSumdeijdaAccumulator[pos] = (PMEAccumulator)0;
    pos += blockDim.x * gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kClearGBBuffers: host function to launch the kernel above
//---------------------------------------------------------------------------------------------
void kClearGBBuffers(gpuContext gpu)
{
#ifdef AMBER_PLATFORM_AMD
  cudaMemsetAsync(gpu->sim.pReffAccumulator, 0, gpu->sim.atoms*sizeof(PMEAccumulator));
  cudaMemsetAsync(gpu->sim.pSumdeijdaAccumulator, 0, gpu->sim.atoms*sizeof(PMEAccumulator));
#else
  kClearGBBuffers_kernel<<<gpu->blocks, gpu->reduceBufferThreadsPerBlock>>>();
#endif
  LAUNCHERROR("kClearGBBuffers");
}
