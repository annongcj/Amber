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
#ifdef AMBER_PLATFORM_AMD
#include <hipcub/hipcub.hpp>
#include "hip_definitions.h"
#else
// CUDA toolkit includes CUB library since 11.0 release
#  if CUDA_VERSION >= 11000
#    define CUB_NS_QUALIFIER ::cub
#    include <cub/cub.cuh>
#  else
//   For CUDA < 11.0 use bundled cub library
#    include "../cub/cub.cuh"
#  endif
#endif

//---------------------------------------------------------------------------------------------
// kNLInitRadixSort: host function to initialize a neighbor list sort.  This is called by
//                   gpu_neighbor_list_setup_ in gpu.cpp, once during a simulation.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------

extern "C" void kNLInitRadixSort(gpuContext gpu)
{
  // Delete old temp data
  cudaDeviceSynchronize();
  delete gpu->pbSortTemp;
  gpu->pbSortTemp = NULL;

  // Create new sort for the Hilbert space neighbor list
  cudaError_t status =
  cub::DeviceRadixSort::SortPairs(NULL, gpu->sortTempBytes,
                                  gpu->sim.pImageHash2, gpu->sim.pImageHash,
                                  gpu->sim.pImageIndex2, gpu->sim.pImageIndex,
                                  gpu->sim.atoms, 0, gpu->neighborListBits);
  RTERROR(status, "kNLInitRadixSort: Sort initialization failed");

  gpu->pbSortTemp = new GpuBuffer<char>(gpu->sortTempBytes);

}

//---------------------------------------------------------------------------------------------
// kNLRadixSort: this is one of the first steps in refreshing the neighbor list, right after
//               creating the spatial hash grid.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kNLRadixSort(gpuContext gpu)
{
  cudaError_t status =
  cub::DeviceRadixSort::SortPairs(gpu->pbSortTemp->_pDevData, gpu->sortTempBytes,
                                  gpu->sim.pImageHash2, gpu->sim.pImageHash,
                                  gpu->sim.pImageIndex2, gpu->sim.pImageIndex,
                                  gpu->sim.atoms, 0, gpu->neighborListBits);
  RTERROR(status, "kNLRadixSort: Sort failed");
}
