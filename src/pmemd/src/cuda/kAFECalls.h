#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// AFE: Alchemical Free Energy
// This could run as part of a set of streaming kernels since very little is actually being
// done here.  We would need to find something else that doesn't use the forces, velocities,
// or positions, though.
//---------------------------------------------------------------------------------------------
#define AFE_Crd
#define AFE_NEIGHBORLIST
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kPMEAFEExchangeCrd_kernel()
#include "kAFEKernels.h"
#undef AFE_NEIGHBORLIST

__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kAFEExchangeCrd_kernel()
#include "kAFEKernels.h"
#undef AFE_Crd

#define AFE_Vel
#define AFE_NEIGHBORLIST
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kPMEAFEExchangeVel_kernel()
#include "kAFEKernels.h"
#undef AFE_NEIGHBORLIST

__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kAFEExchangeVel_kernel()
#include "kAFEKernels.h"
#undef AFE_Vel

#define AFE_Frc
#define AFE_NEIGHBORLIST
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kPMEAFEExchangeFrc_kernel()
#include "kAFEKernels.h"
#undef AFE_NEIGHBORLIST

__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kAFEExchangeFrc_kernel()
#include "kAFEKernels.h"
#undef AFE_Frc

//---------------------------------------------------------------------------------------------
// kAFEExchange$VECTOR: exchange information between two variants of a system.
//
// Arguments:
//   gpu:     overarching type for storing parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kAFEExchangeVel(gpuContext gpu)
{
  if (gpu->bNeighborList) {
    kPMEAFEExchangeVel_kernel<<<gpu->AFEExchangeBlocks, gpu->AFEExchangeThreadsPerBlock>>>();
  }
  else {
    kAFEExchangeVel_kernel<<<gpu->AFEExchangeBlocks, gpu->AFEExchangeThreadsPerBlock>>>();
  }
  LAUNCHERROR("kAFEExchangeVel");
}

//---------------------------------------------------------------------------------------------
void kAFEExchangeFrc(gpuContext gpu)
{
  if (gpu->bNeighborList) {
    kPMEAFEExchangeFrc_kernel<<<gpu->AFEExchangeBlocks, gpu->AFEExchangeThreadsPerBlock>>>();
  }
  else {
    kAFEExchangeFrc_kernel<<<gpu->AFEExchangeBlocks, gpu->AFEExchangeThreadsPerBlock>>>();
  }
  LAUNCHERROR("kAFEExchangeFrc");
}

//---------------------------------------------------------------------------------------------
void kAFEExchangeCrd(gpuContext gpu)
{
  if (gpu->bNeighborList) {
    kPMEAFEExchangeCrd_kernel<<<gpu->AFEExchangeBlocks, gpu->AFEExchangeThreadsPerBlock>>>();
  }
  else {
    kAFEExchangeCrd_kernel<<<gpu->AFEExchangeBlocks, gpu->AFEExchangeThreadsPerBlock>>>();
  }
  LAUNCHERROR("kAFEExchangeCrd");
}
