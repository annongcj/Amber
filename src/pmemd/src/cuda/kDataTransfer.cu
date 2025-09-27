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

#if !defined(__HIPCC_RDC__)

//---------------------------------------------------------------------------------------------
// SetkDataTransferSim: called by gpuCopyConstants to orient the instance of cSim in this
//                      library
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void SetkDataTransferSim(gpuContext gpu)
{
  cudaError_t status;
  status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(simulationConst));
  RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

//---------------------------------------------------------------------------------------------
// GetkDataTransferSim: download information about the CUDA simulation data struct from the
//                      device.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//
// This is a debugging function.
//---------------------------------------------------------------------------------------------
void GetkDataTransferSim(gpuContext gpu)
{
  cudaError_t status;
  status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(simulationConst));
  RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy to cSim failed");
}

#endif

//---------------------------------------------------------------------------------------------
// Kernels for the "travel light" option of transferring coordinates and forces between the
// CPU and GPU
//---------------------------------------------------------------------------------------------
#define TAKE_IMAGE
#define PULL_COORD
__global__ void kRetrieveSimCoordNL_kernel()
#include "kShuttle.h"
#undef PULL_COORD

#define PULL_FORCE
__global__ void kRetrieveSimForceNL_kernel()
#include "kShuttle.h"
#undef PULL_FORCE

#define PULL_BOND_FORCE
__global__ void kRetrieveSimBondForceNL_kernel()
#include "kShuttle.h"
#undef PULL_BOND_FORCE

#define PULL_NONBOND_FORCE
__global__ void kRetrieveSimNonBondForceNL_kernel()
#include "kShuttle.h"
#undef PULL_NONBOND_FORCE

#define POST_FORCE
__global__ void kPostSimForceNL_kernel()
#include "kShuttle.h"
#undef POST_FORCE

#define POST_NONBOND_FORCE
__global__ void kPostSimNonBondForceNL_kernel()
#include "kShuttle.h"
#undef POST_NONBOND_FORCE

#define POST_FORCE_ADD
__global__ void kPostSimForceAddNL_kernel()
#include "kShuttle.h"
#undef POST_FORCE_ADD

#define POST_NONBOND_FORCE_ADD
__global__ void kPostSimNonBondForceAddNL_kernel()
#include "kShuttle.h"
#undef POST_NONBOND_FORCE_ADD

#define POST_COORD
__global__ void kPostSimCoordNL_kernel()
#include "kShuttle.h"
#undef POST_COORD
#undef TAKE_IMAGE

#define PULL_COORD
__global__ void kRetrieveSimCoord_kernel()
#include "kShuttle.h"
#undef PULL_COORD

#define PULL_FORCE
__global__ void kRetrieveSimForce_kernel()
#include "kShuttle.h"
#undef PULL_FORCE

#define PULL_BOND_FORCE
__global__ void kRetrieveSimBondForce_kernel()
#include "kShuttle.h"
#undef PULL_FORCE

#define PULL_NONBOND_FORCE
__global__ void kRetrieveSimNonBondForce_kernel()
#include "kShuttle.h"
#undef PULL_NONBOND_FORCE

#define POST_FORCE
__global__ void kPostSimForce_kernel()
#include "kShuttle.h"
#undef POST_FORCE

#define POST_NONBOND_FORCE
__global__ void kPostSimNonBondForce_kernel()
#include "kShuttle.h"
#undef POST_NONBOND_FORCE

#define POST_FORCE_ADD
__global__ void kPostSimForceAdd_kernel()
#include "kShuttle.h"
#undef POST_FORCE_ADD

#define POST_NONBOND_FORCE_ADD
__global__ void kPostSimNonBondForceAdd_kernel()
#include "kShuttle.h"
#undef POST_NONBOND_FORCE_ADD

#define POST_COORD
__global__ void kPostSimCoord_kernel()
#include "kShuttle.h"
#undef POST_COORD

//---------------------------------------------------------------------------------------------
// kRetrieveSimData: host function to download data (i.e. coordinates) on particular atoms
//                   from the GPU to the CPU RAM for processing by the main processor.
//
// Arguments:
//   gpu:       Overarching type for storing parameters, coordinates, and the energy function.
//              For this context, it stores a notable parameter: sim.Shuttletype...
//                0.) Information on atoms in a static list stored on the CPU and previously
//                    copied to the GPU
//                1.) Properties and identities of a number of atoms computed on the fly (on
//                    the GPU) that will participate in some sort of interaction to be computed
//                    on the CPU (this mode is not yet complete)
//   atm_data:  array holding atom data as handled by the Fortran code on the CPU (acceptable
//              inputs include atm_crd[][3], and atm_frc[][3], for example)
//   modifier:  specifies the type of information that will be shuttled in any of the modes
//              defined by gpu->sim.ShuttleType.  0: coordinates come down, 1: aggregate forces
//              come down, 2: forces due to bonded interactions come down, 3: forces due to
//              non-bonded interactions come down
//---------------------------------------------------------------------------------------------
extern "C" void kRetrieveSimData(gpuContext gpu, double atm_data[][3], int modifier)
{
  int i;

  if (gpu->sim.ShuttleType == 0) {

    // Load the data transfer buffer on the device
    if (modifier == 0) {
      if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
        kRetrieveSimCoordNL_kernel<<<gpu->blocks, GRID * 3>>>();
      }
      else {
        kRetrieveSimCoord_kernel<<<gpu->blocks, GRID * 3>>>();
      }
    }
    else if (modifier == 1) {
      if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
        kRetrieveSimForceNL_kernel<<<gpu->blocks, GRID * 3>>>();
      }
      else {
        kRetrieveSimForce_kernel<<<gpu->blocks, GRID * 3>>>();
      }
    }
    else if (modifier == 2) {
      if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
        kRetrieveSimBondForceNL_kernel<<<gpu->blocks, GRID * 3>>>();
      }
      else {
        kRetrieveSimBondForce_kernel<<<gpu->blocks, GRID * 3>>>();
      }
    }
    else if (modifier == 3) {
      if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
        kRetrieveSimNonBondForceNL_kernel<<<gpu->blocks, GRID * 3>>>();
      }
      else {
        kRetrieveSimNonBondForce_kernel<<<gpu->blocks, GRID * 3>>>();
      }
    }

    // Download from the device
    gpu->pbDataShuttle->Download();

    // Place atoms in the standard coordinate array
    int npt = gpu->sim.nShuttle;
    int *atmid = gpu->pbShuttleTickets->_pSysData;
    double *pShData = gpu->pbDataShuttle->_pSysData;
    for (i = 0; i < npt; i++) {
      atm_data[atmid[i]][0] = pShData[        i];
      atm_data[atmid[i]][1] = pShData[  npt + i];
      atm_data[atmid[i]][2] = pShData[2*npt + i];
    }
  }
  else if (gpu->sim.ShuttleType == 1) {

    // Print an error for now
    printf("| Error: this functionality (ShuttleType == 1) is not yet available.\n");
  }
}

//---------------------------------------------------------------------------------------------
// kPostSimData: host function to upload data (i.e. coordinates) on particular atoms to the
//               GPU from the CPU RAM for processing by the main processor.
//
// Arguments:
//   gpu:       overarching type for storing all parameters, coordinates, and energy functions.
//              For this context, it stores a notable parameter: sim.Shuttletype...
//                0.) Information on atoms in a static list stored on the CPU and previously
//                    copied to the GPU
//                1.) Properties and identities of a pre-defined number of atoms that will
//                    participate in some sort of interaction to be computed on the CPU
//   atm_data:  array holding atom data as handled by the Fortran code on the CPU (acceptable
//              inputs include atm_frc[][3], and atm_crd[][3], for example)
//   modifier:  specifies the type of information that will be shuttled in any of the modes
//              defined by gpu->sim.ShuttleType.
//
//              Mod   Data type      Destination on GPU                   Behavior
//              ---  -----------  ------------------------    ---------------------------------
//               0   Coordinates  pImage{X,Y,Z}               Set values on GPU with CPU data
//               1   Forces       pForceAccumulator{X,Y,Z}    Set values on GPU with CPU data
//               2   Bond forces  pForceAccumulator{X,Y,Z}    Set values on GPU with CPU data
//               3   NB forces    pNBForceAccumulator{X,Y,Z}  Set values on GPU with CPU data
//               4   Forces       pForceAccumulator{X,Y,Z}    Add CPU data to values on the GPU
//               5   Bond forces  pForceAccumulator{X,Y,Z}    Add CPU data to values on the GPU
//               4   NB forces    pNBForceAccumulator{X,Y,Z}  Add CPU data to values on the GPU
//---------------------------------------------------------------------------------------------
extern "C" void kPostSimData(gpuContext gpu, double atm_data[][3], int modifier)
{
  int i;

  if (gpu->sim.ShuttleType == 0) {

    // Load the data transfer buffer on the host
    int npt = gpu->sim.nShuttle;
    int *atmid = gpu->pbShuttleTickets->_pSysData;
    double *pShData = gpu->pbDataShuttle->_pSysData;
    for (i = 0; i < npt; i++) {
      pShData[        i] = atm_data[atmid[i]][0];
      pShData[  npt + i] = atm_data[atmid[i]][1];
      pShData[2*npt + i] = atm_data[atmid[i]][2];
    }

    // Upload data to the device and unpack
    gpu->pbDataShuttle->Upload();
    if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
      if (modifier == 0) {
        kPostSimCoordNL_kernel<<<gpu->blocks, GRID * 3>>>();
      }
      else if (modifier == 1 || modifier == 2) {
        kPostSimForceNL_kernel<<<gpu->blocks, GRID * 3>>>();
      }
      else if (modifier == 3) {
        kPostSimNonBondForceNL_kernel<<<gpu->blocks, GRID * 3>>>();
      }
      else if (modifier == 4 || modifier == 5) {
        kPostSimForceAddNL_kernel<<<gpu->blocks, GRID * 3>>>();
      }
      else if (modifier == 6) {
        kPostSimNonBondForceAddNL_kernel<<<gpu->blocks, GRID * 3>>>();
      }
    }
    else {
      if (modifier == 0) {
        kPostSimCoord_kernel<<<gpu->blocks, GRID * 3>>>();
      }
      else if (modifier == 1 || modifier == 2) {
        kPostSimForce_kernel<<<gpu->blocks, GRID * 3>>>();
      }
      else if (modifier == 3) {
        kPostSimNonBondForce_kernel<<<gpu->blocks, GRID * 3>>>();
      }
      else if (modifier == 4 || modifier == 5) {
        kPostSimForceAdd_kernel<<<gpu->blocks, GRID * 3>>>();
      }
      else if (modifier == 6) {
        kPostSimNonBondForceAdd_kernel<<<gpu->blocks, GRID * 3>>>();
      }
    }
  }
}
