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

// Use global instance instead of a local copy
#include "simulationConst.h"
CSIM_STO simulationConst cSim;

//---------------------------------------------------------------------------------------------
// Atom: structure to store SHAKE-critical atom position and mass information
//---------------------------------------------------------------------------------------------
struct Atom
{
  double invMassI;
  double xpl;
  double ypl;
  double zpl;
  double xil;
  double yil;
  double zil;
};

#if !defined(__HIPCC_RDC__)

//---------------------------------------------------------------------------------------------
// SetkShakeSim: upload critical SHAKE data to the GPU
//---------------------------------------------------------------------------------------------
void SetkShakeSim(gpuContext gpu)
{
  cudaError_t status;
  status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(simulationConst));
  RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

//---------------------------------------------------------------------------------------------
// GetkShakeSim: download critical SHAKE data from the GPU
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
void GetkShakeSim(gpuContext gpu)
{
  cudaError_t status;
  status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(simulationConst));
  RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

#endif

//---------------------------------------------------------------------------------------------
// Kernels for general, PME-specialized, and hardware-dependent flavors of SHAKE, with and
// without hydrogen mass repartitioning (HMR).
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
kShake_kernel()
#include "kShake.h"

#define TISHAKE2
__global__ void
__launch_bounds__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
kTIShake2_kernel()
#include "kShake.h"
#undef TISHAKE2

//---------------------------------------------------------------------------------------------
#define SHAKE_NEIGHBORLIST
#define TISHAKE2
__global__ void
__launch_bounds__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
kPMETIShake2_kernel()
#include "kShake.h"
#undef TISHAKE2
__global__ void
__LAUNCH_BOUNDS__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
kPMEShake_kernel()
#include "kShake.h"
#undef SHAKE_NEIGHBORLIST

//---------------------------------------------------------------------------------------------
#define NODPTEXTURE
#define TISHAKE2
__global__ void
__launch_bounds__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
kTIShake2NoDPTexture_kernel()
#include "kShake.h"
#undef TISHAKE2
__global__ void
__LAUNCH_BOUNDS__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
kShakeNoDPTexture_kernel()
#include "kShake.h"

//---------------------------------------------------------------------------------------------
#define SHAKE_NEIGHBORLIST
#define TISHAKE2
__global__ void
__launch_bounds__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
kPMETIShake2NoDPTexture_kernel()
#include "kShake.h"
#undef TISHAKE2
__global__ void
__LAUNCH_BOUNDS__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
kPMEShakeNoDPTexture_kernel()
#include "kShake.h"
#undef SHAKE_NEIGHBORLIST
#undef NODPTEXTURE

//---------------------------------------------------------------------------------------------
#define SHAKE_HMR
#define TISHAKE2
__global__ void
__launch_bounds__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
kTIShake2HMR_kernel()
#include "kShake.h"
#undef TISHAKE2
__global__ void
__LAUNCH_BOUNDS__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
kShakeHMR_kernel()
#include "kShake.h"

//---------------------------------------------------------------------------------------------
#define SHAKE_NEIGHBORLIST
#define TISHAKE2
__global__ void
__launch_bounds__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
kPMETIShake2HMR_kernel()
#include "kShake.h"
#undef TISHAKE2
__global__ void
__LAUNCH_BOUNDS__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
kPMEShakeHMR_kernel()
#include "kShake.h"
#undef SHAKE_NEIGHBORLIST

//---------------------------------------------------------------------------------------------
#define NODPTEXTURE
#define TISHAKE2
__global__ void
__launch_bounds__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
kTIShake2HMRNoDPTexture_kernel()
#include "kShake.h"
#undef TISHAKE2
__global__ void
__LAUNCH_BOUNDS__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
kShakeHMRNoDPTexture_kernel()
#include "kShake.h"

//---------------------------------------------------------------------------------------------
#define SHAKE_NEIGHBORLIST
#define TISHAKE2
__global__ void
__launch_bounds__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
kPMETIShake2HMRNoDPTexture_kernel()
#include "kShake.h"
#undef TISHAKE2
__global__ void
__LAUNCH_BOUNDS__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
kPMEShakeHMRNoDPTexture_kernel()
#include "kShake.h"
#undef SHAKE_NEIGHBORLIST
#undef NODPTEXTURE
#undef SHAKE_HMR

__global__ void
__launch_bounds__(UPDATE_THREADS_PER_BLOCK, 1)
kShakeOldPos_kernel()
#include "kOldShakePos.h"

//---------------------------------------------------------------------------------------------
#define SHAKE_NEIGHBORLIST
__global__ void
__launch_bounds__(UPDATE_THREADS_PER_BLOCK, 1)
kPMEShakeOldPos_kernel()
#include "kOldShakePos.h"
#undef SHAKE_NEIGHBORLIST

//---------------------------------------------------------------------------------------------
// Kernels for general, PME-specialized, and hardware-dependent flavors of RATTLE, with and
// without hydrogen mass repartitioning (HMR).
//---------------------------------------------------------------------------------------------
__global__ void
__launch_bounds__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
kRattle_kernel(PMEFloat dt)
#include "kRattle.h"

//---------------------------------------------------------------------------------------------
#define RATTLE_NEIGHBORLIST
__global__ void
__launch_bounds__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
kPMERattle_kernel(PMEFloat dt)
#include "kRattle.h"
#undef RATTLE_NEIGHBORLIST

//---------------------------------------------------------------------------------------------
//#define NODPTEXTURE
//__global__ void
//__launch_bounds__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
//kShakeNoDPTexture_kernel()
//#include "kShake.h"

//---------------------------------------------------------------------------------------------
//#define SHAKE_NEIGHBORLIST
//__global__ void
//__launch_bounds__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
//kPMEShakeNoDPTexture_kernel()
//#include "kShake.h"
//#undef SHAKE_NEIGHBORLIST
//#undef NODPTEXTURE

//---------------------------------------------------------------------------------------------
#define RATTLE_HMR
__global__ void
__launch_bounds__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
kRattleHMR_kernel(PMEFloat dt)
#include "kRattle.h"

//---------------------------------------------------------------------------------------------
#define RATTLE_NEIGHBORLIST
__global__ void
__launch_bounds__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
kPMERattleHMR_kernel(PMEFloat dt)
#include "kRattle.h"
#undef RATTLE_NEIGHBORLIST

//---------------------------------------------------------------------------------------------
//#define NODPTEXTURE
//__global__ void
//__launch_bounds__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
//kShakeHMRNoDPTexture_kernel()
//#include "kShake.h"

//---------------------------------------------------------------------------------------------
//#define SHAKE_NEIGHBORLIST
//__global__ void
//__launch_bounds__(SHAKE_THREADS_PER_BLOCK, SHAKE_BLOCKS)
//kPMEShakeHMRNoDPTexture_kernel()
//#include "kShake.h"
//#undef SHAKE_NEIGHBORLIST
//#undef NODPTEXTURE
#undef RATTLE_HMR

//---------------------------------------------------------------------------------------------
// kShakeInitKernels: initialize SHAKE-reltaed kernels
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kShakeInitKernels(gpuContext gpu)
{
  if (gpu->sm_version >= SM_3X) {
    cudaFuncSetCacheConfig(kShake_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kPMEShake_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kShakeNoDPTexture_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kPMEShakeNoDPTexture_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kShakeHMR_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kPMEShakeHMR_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kShakeHMRNoDPTexture_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kPMEShakeHMRNoDPTexture_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kTIShake2_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kPMETIShake2_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kTIShake2NoDPTexture_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kPMETIShake2NoDPTexture_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kTIShake2HMR_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kPMETIShake2HMR_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kTIShake2HMRNoDPTexture_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kPMETIShake2HMRNoDPTexture_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kRattle_kernel, cudaFuncCachePreferL1); //added by zhf, used for middle-scheme
    cudaFuncSetCacheConfig(kPMERattle_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kRattleHMR_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kPMERattleHMR_kernel, cudaFuncCachePreferL1);
    cudaFuncSetSharedMemConfig(kShake_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kPMEShake_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kShakeNoDPTexture_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kPMEShakeNoDPTexture_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kShakeHMR_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kPMEShakeHMR_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kShakeHMRNoDPTexture_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kPMEShakeHMRNoDPTexture_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kTIShake2_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kPMETIShake2_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kTIShake2NoDPTexture_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kPMETIShake2NoDPTexture_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kTIShake2HMR_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kPMETIShake2HMR_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kTIShake2HMRNoDPTexture_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kPMETIShake2HMRNoDPTexture_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kRattle_kernel, cudaSharedMemBankSizeEightByte); //added by zhf, used for middle-scheme
    cudaFuncSetSharedMemConfig(kPMERattle_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kRattleHMR_kernel, cudaSharedMemBankSizeEightByte);
    cudaFuncSetSharedMemConfig(kPMERattleHMR_kernel, cudaSharedMemBankSizeEightByte);
  }
}

//---------------------------------------------------------------------------------------------
// kShake: implement the various SHAKE kernels enumerated above.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kShake(gpuContext gpu)
{
  unsigned int totalBlocks = (gpu->sim.slowShakeOffset + gpu->shakeThreadsPerBlock - 1) /
                             gpu->shakeThreadsPerBlock;
  if(gpu->sim.tishake == 2) {
    totalBlocks = (gpu->sim.slowTIShakeOffset + gpu->shakeThreadsPerBlock - 1) /
                             gpu->shakeThreadsPerBlock;
  }
  int launchBlocks = 65535;

  if (gpu->bNeighborList) {
      kPMEShakeOldPos_kernel<<<gpu->updateBlocks, gpu->updateThreadsPerBlock>>>();
  }
  else {
      kShakeOldPos_kernel<<<gpu->updateBlocks, gpu->updateThreadsPerBlock>>>();
  }


  if (gpu->bNoDPTexture) {
    while (totalBlocks > 0) {
      int blocks = min(totalBlocks, launchBlocks);
      if (gpu->bUseHMR) {
        if (gpu->bNeighborList) {
          if (gpu->sim.tishake == 2) {
            kPMETIShake2HMRNoDPTexture_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>();
          }
          else {
            kPMEShakeHMRNoDPTexture_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>();
          }
        }
        else {
          if (gpu->sim.tishake == 2) {
            kTIShake2HMRNoDPTexture_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>();
          }
          else {
            kShakeHMRNoDPTexture_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>();
          }
        }
      }
      else {
        if (gpu->bNeighborList) {
          if (gpu->sim.tishake == 2) {
            kPMETIShake2NoDPTexture_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>();
          }
          else {
            kPMEShakeNoDPTexture_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>();
          }
        }
        else {
          if (gpu->sim.tishake == 2) {
            kTIShake2NoDPTexture_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>();
          }
          else {
            kShakeNoDPTexture_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>();
          }
        }
      }
      LAUNCHERROR("kShakeNoDPTexture");
      totalBlocks -= blocks;
    }
  }
  else {
    while (totalBlocks > 0) {
      int blocks = min(totalBlocks, launchBlocks);
      if (gpu->bUseHMR) {
        if (gpu->bNeighborList) {
          if (gpu->sim.tishake == 2) {
            kPMETIShake2HMR_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>();
          }
          else {
            kPMEShakeHMR_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>();
          }
        }
        else {
          if (gpu->sim.tishake == 2) {
            kTIShake2HMR_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>();
          }
          else {
            kShakeHMR_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>();
          }
        }
      }
      else {
        if (gpu->bNeighborList) {
          if (gpu->sim.tishake == 2) {
            kPMETIShake2_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>();
          }
          else {
            kPMEShake_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>();
          }
        }
        else {
          if (gpu->sim.tishake == 2) {
            kTIShake2_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>();
          }
          else {
            kShake_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>();
          }
        }
      }
      LAUNCHERROR("kShake");
      totalBlocks -= blocks;
    }
  }
}

//---------------------------------------------------------------------------------------------
// kRattle: implement the various RATTLE kernels enumerated above.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kRattle(gpuContext gpu, PMEFloat dt)
{
  unsigned int totalBlocks = (gpu->sim.slowShakeOffset + gpu->shakeThreadsPerBlock - 1) /
                             gpu->shakeThreadsPerBlock;
  int launchBlocks = 65535;
  //if (gpu->bNoDPTexture) {
  //  while (totalBlocks > 0) {
  //    int blocks = min(totalBlocks, launchBlocks);
  //    if (gpu->bUseHMR) {
  //      if (gpu->bNeighborList) {
  //        kPMEShakeHMRNoDPTexture_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>();
  //      }
  //      else {
  //        kShakeHMRNoDPTexture_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>();
  //      }
  //    }
  //    else {
  //      if (gpu->bNeighborList) {
  //        kPMEShakeNoDPTexture_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>();
  //      }
  //      else {
  //        kShakeNoDPTexture_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>();
  //      }
  //    }
  //    LAUNCHERROR("kShakeNoDPTexture");
  //    totalBlocks -= blocks;
  //  }
  //}
  //else {
    while (totalBlocks > 0) {
      int blocks = min(totalBlocks, launchBlocks);
      if (gpu->bUseHMR) {
        if (gpu->bNeighborList) {
          kPMERattleHMR_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>(dt);
        }
        else {
          kRattleHMR_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>(dt);
        }
      }
      else {
        if (gpu->bNeighborList) {
          kPMERattle_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>(dt);
        }
        else {
          kRattle_kernel<<<blocks, gpu->shakeThreadsPerBlock>>>(dt);
        }
      }
      LAUNCHERROR("kRattle");
      totalBlocks -= blocks;
    }
  //}
}
