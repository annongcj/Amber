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
//#include "cuda_profiler_api.h"

//#define PME_ENERGY

// Use global instance instead of a local copy
#include "simulationConst.h"
CSIM_STO simulationConst cSim;

#if !defined(__HIPCC_RDC__)

//---------------------------------------------------------------------------------------------
// SetkCalculateEFieldEnergySim: copy the contents of the host-side simulation command buffer
//                               to the device.
//
// Arguments:
//   gpu:      overarching type for storing all parameters, coordinates, and energy function
//             terms in a simulation
//---------------------------------------------------------------------------------------------
void SetkCalculateEFieldEnergySim(gpuContext gpu)
{
  cudaError_t status;
  status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(simulationConst));
  RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

//---------------------------------------------------------------------------------------------
// GetkCalculateEFieldEnergySim: download the contents of a CUDA simulation from the device to
//                               the host.
//
// Arguments:
//   gpu:      overarching type for storing all parameters, coordinates, and energy function
//             terms in a simulation
//---------------------------------------------------------------------------------------------
void GetkCalculateEFieldEnergySim(gpuContext gpu)
{
  cudaError_t status;
  status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(simulationConst));
  RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

#endif

// EField kernels

#define PME_ENERGY

//---------------------------------------------------------------------------------------------
// kCalculateEFieldEnergy_kernel: kernel for computing the electrostatic field (PME) energy of
//                                a system of particles
//
// Arguments:
//   nstep:     the step number (i.e. 250307 of 500000)
//   dt:        the time step
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kCalculateEFieldEnergy_kernel(PMEDouble nstep, PMEDouble dt)
#include "kEFE.h"

#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
// kCalculateEFieldForces_kernel: kernel for computing the electrostatic field (PME) forces on
//                                a system of particles
//
// Arguments:
//   nstep:     the step number (i.e. 250307 of 500000)
//   dt:        the time step
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kCalculateEFieldForces_kernel(PMEDouble nstep, PMEDouble dt)
#include "kEFE.h"

//---------------------------------------------------------------------------------------------
// kCalculateEFieldForces: wrapper function for the eponymous kernel
//
// Arguments:
//   gpu:      overarching type for storing all parameters, coordinates, and energy function
//             terms in a simulation
//   nstep:    the step number (i.e. 250307 of 500000)
//   dt:       the time step
//---------------------------------------------------------------------------------------------
extern "C" void kCalculateEFieldForces(gpuContext gpu, int nstep, PMEDouble dt)
{
  kCalculateEFieldForces_kernel<<<gpu->updateBlocks,
                                  gpu->updateThreadsPerBlock>>>((PMEDouble)nstep, dt);
  LAUNCHERROR("kCalculateEFieldForces");
}

//---------------------------------------------------------------------------------------------
// kCalculateEFieldEnergy: wrapper function for the eponymous kernel
//
// Arguments:
//   gpu:      overarching type for storing all parameters, coordinates, and energy function
//             terms in a simulation
//   nstep:    the step number (i.e. 250307 of 500000)
//   dt:       the time step
//---------------------------------------------------------------------------------------------
extern "C" void kCalculateEFieldEnergy(gpuContext gpu, int nstep, PMEDouble dt)
{
  kCalculateEFieldEnergy_kernel<<<gpu->updateBlocks,
                                  gpu->updateThreadsPerBlock>>>((PMEDouble)nstep, dt);
  LAUNCHERROR("kCalculatePMENonbondEnergy");
}

//---------------------------------------------------------------------------------------------
// kCalculateEFieldEnergyInitKernels: initialize PME kernels for a simulation
//
// Arguments:
//   gpu:      overarching type for storing all parameters, coordinates, and energy function
//             terms in a simulation
//---------------------------------------------------------------------------------------------
extern "C" void kCalculateEFieldEnergyInitKernels(gpuContext gpu)
{
  cudaFuncSetSharedMemConfig(kCalculateEFieldEnergy_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalculateEFieldForces_kernel, cudaSharedMemBankSizeEightByte);
}
