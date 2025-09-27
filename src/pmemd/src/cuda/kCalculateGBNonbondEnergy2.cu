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

static __constant__ PMEFloat te  = (PMEFloat)(4.0 / 3.0);
static __constant__ PMEFloat tf  = (PMEFloat)(12.0 / 5.0);
static __constant__ PMEFloat tg  = (PMEFloat)(24.0 / 7.0);
static __constant__ PMEFloat th  = (PMEFloat)(40.0 / 9.0);
static __constant__ PMEFloat thh = (PMEFloat)(60.0 / 11.0);

#if !defined(__HIPCC_RDC__)

//---------------------------------------------------------------------------------------------
// SetkCalculateGBNonbondEnergy2Sim: upload information for the second stage of GB energy
//                                   calculation to the GPU.  Called by the most likely
//                                   function you'd suspect: gpuCopyConstants() in gpu.cpp.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void SetkCalculateGBNonbondEnergy2Sim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(simulationConst));
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

//---------------------------------------------------------------------------------------------
// GetkCalculateGBNonBondEnergy2Sim: download information for the second stage of GB energy
//                                   calculation from the GPU.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
void GetkCalculateGBNonBondEnergy2Sim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(simulationConst));
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

#endif

//---------------------------------------------------------------------------------------------
// Two kernels for calculating GB forces and energies in various flavors.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GBNONBONDENERGY2_THREADS_PER_BLOCK, GBNONBONDENERGY2_BLOCKS_LOAD)
kCalculateGBNonbondEnergy2_kernel()
#include "kCalculateGBNonbondEnergy2.h"

//---------------------------------------------------------------------------------------------
#define GB_IGB78
__global__ void
__LAUNCH_BOUNDS__(GBNONBONDENERGY2_THREADS_PER_BLOCK, GBNONBONDENERGY2_BLOCKS_LOAD)
kCalculateGBNonbondEnergy2IGB78_kernel()
#include "kCalculateGBNonbondEnergy2.h"
#undef IGB78

//---------------------------------------------------------------------------------------------
// kCalculateGBNonbondEnergy2InitKernels: initialize the kernels for GB force and energies.
//                                        Like its counterpart in the PME library, this
//                                        function serves to set the __shared__ memory stride
//                                        to eight bytes.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kCalculateGBNonbondEnergy2InitKernels(gpuContext gpu)
{
  cudaFuncSetSharedMemConfig(kCalculateGBNonbondEnergy2_kernel,
                             cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalculateGBNonbondEnergy2IGB78_kernel,
                             cudaSharedMemBankSizeEightByte);
  if (gpu->sm_version >= SM_3X) {
    cudaFuncSetCacheConfig(kCalculateGBNonbondEnergy2_kernel, cudaFuncCachePreferEqual);
  }
}

//---------------------------------------------------------------------------------------------
// kCalculateGBNonbondEnergy2: launch the appropriate kernel for computing GB forces and
//                             energies.
//---------------------------------------------------------------------------------------------
void kCalculateGBNonbondEnergy2(gpuContext gpu)
{
  if ((gpu->sim.igb == 7) || (gpu->sim.igb == 8)) {
    kCalculateGBNonbondEnergy2IGB78_kernel<<<gpu->GBNonbondEnergy2Blocks,
                                             gpu->GBNonbondEnergy2IGB78ThreadsPerBlock>>>();
  }
  else {
    kCalculateGBNonbondEnergy2_kernel<<<gpu->GBNonbondEnergy2Blocks,
                                        gpu->GBNonbondEnergy2ThreadsPerBlock>>>();
  }
  LAUNCHERROR("kCalculateGBNonbondEnergy2");
}
