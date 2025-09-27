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

// Use global instance instead of a local copy
#include "simulationConst.h"
CSIM_STO simulationConst cSim;

#ifndef use_DPFP
//Taken from kNLCPNE.h
//---------------------------------------------------------------------------------------------
// cudaERFC: struct to hold the ten constants required by the fasterfc function below
//---------------------------------------------------------------------------------------------
struct cudaERFC
{
  PMEFloat c0;
  PMEFloat c1;
  PMEFloat c2;
  PMEFloat c3;
  PMEFloat c4;
  PMEFloat c5;
  PMEFloat c6;
  PMEFloat c7;
  PMEFloat c8;
  PMEFloat c9;
};
__device__ __constant__ cudaERFC cCoarseGridERFC;

//---------------------------------------------------------------------------------------------
// SetkCalculateCoarseGridEnergyERFC: set constants for the fasterfc() function on the device.
//
// Arguments:
//   ewcoeff:    the Ewald coefficient, 1/(2*Gaussian sigma)
//---------------------------------------------------------------------------------------------
void SetkCalculateCoarseGridEnergyERFC(double ewcoeff)
{
  cudaERFC c;
  c.c0 = (PMEFloat)(-1.6488499458192755E-006 * pow(ewcoeff, 10.0));
  c.c1 = (PMEFloat)( 2.9524665006554534E-005 * pow(ewcoeff,  9.0));
  c.c2 = (PMEFloat)(-2.3341951153749626E-004 * pow(ewcoeff,  8.0));
  c.c3 = (PMEFloat)( 1.0424943374047289E-003 * pow(ewcoeff,  7.0));
  c.c4 = (PMEFloat)(-2.5501426008983853E-003 * pow(ewcoeff,  6.0));
  c.c5 = (PMEFloat)( 3.1979939710877236E-004 * pow(ewcoeff,  5.0));
  c.c6 = (PMEFloat)( 2.7605379075746249E-002 * pow(ewcoeff,  4.0));
  c.c7 = (PMEFloat)(-1.4827402067461906E-001 * pow(ewcoeff,  3.0));
  c.c8 = (PMEFloat)(-9.1844764013203406E-001 * ewcoeff * ewcoeff);
  c.c9 = (PMEFloat)(-1.6279070384382459E+000 * ewcoeff);
  cudaError_t status;
  status = cudaMemcpyToSymbol(cCoarseGridERFC, &c, sizeof(cudaERFC));
  RTERROR(status, "cudaMemcpyToSymbol: SetERFC copy to cCoarseGridERFC failed");
}

//---------------------------------------------------------------------------------------------
// __internal_fmad: encapsulates a call to __fmaf_rn (FMAD = Floating point Multiply-ADd), an
//                  intrinsic that returns a*b + c by the nomenclature given below.
//---------------------------------------------------------------------------------------------
static __forceinline__ __device__ float __internal_fmad(float a, float b, float c)
{
  return __fmaf_rn (a, b, c);
}

//---------------------------------------------------------------------------------------------
// Faster ERFC approximation courtesy of Norbert Juffa, NVIDIA Corporation
//
// Arguments:
//   a:     the argument to erfc(a)--take the complimentary error function of the number a
//---------------------------------------------------------------------------------------------
static __forceinline__ __device__ PMEFloat fasterfc(PMEFloat a)
{
  // Approximate log(erfc(a)) with rel. error < 7e-9
  PMEFloat t, x = a;
  t = cCoarseGridERFC.c0;
  t = __internal_fmad(t, x, cCoarseGridERFC.c1);
  t = __internal_fmad(t, x, cCoarseGridERFC.c2);
  t = __internal_fmad(t, x, cCoarseGridERFC.c3);
  t = __internal_fmad(t, x, cCoarseGridERFC.c4);
  t = __internal_fmad(t, x, cCoarseGridERFC.c5);
  t = __internal_fmad(t, x, cCoarseGridERFC.c6);
  t = __internal_fmad(t, x, cCoarseGridERFC.c7);
  t = __internal_fmad(t, x, cCoarseGridERFC.c8);
  t = __internal_fmad(t, x, cCoarseGridERFC.c9);
  t = t * x;
  return exp2f(t);
}

#endif

//---------------------------------------------------------------------------------------------
// SetkCalculateCoarseGridEnergySim: copy the contents of the host-side simulation command buffer
//                               to the device.
//
// Arguments:
//   gpu:      overarching type for storing all parameters, coordinates, and energy function
//             terms in a simulation
//---------------------------------------------------------------------------------------------
void SetkCalculateCoarseGridEnergySim(gpuContext gpu)
{
#if !defined(__HIPCC_RDC__)
  cudaError_t status;
  status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(simulationConst));
  RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
#endif
}

//---------------------------------------------------------------------------------------------
// GetkCalculateCoarseGridEnergySim: download the contents of a CUDA simulation from the device to
//                               the host.
//
// Arguments:
//   gpu:      overarching type for storing all parameters, coordinates, and energy function
//             terms in a simulation
//---------------------------------------------------------------------------------------------
void GetkCalculateCoarseGridEnergySim(gpuContext gpu)
{
#if !defined(__HIPCC_RDC__)
  cudaError_t status;
  status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(simulationConst));
  RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
#endif
}

//---------------------------------------------------------------------------------------------
// kCalculateCoarseGridEnergy_kernel: kernel for computing a pre-filter for mcwat
//
// Arguments:
//---------------------------------------------------------------------------------------------

#define TI
#define IS_ORTHOG
__global__ void
__launch_bounds__(UPDATE_THREADS_PER_BLOCK, 1)
kCalculateCoarseGridEnergyOrthogTI_kernel()
#include "kCCGE.h"
#undef IS_ORTHOG
__global__ void
__launch_bounds__(UPDATE_THREADS_PER_BLOCK, 1)
kCalculateCoarseGridEnergyNonOrthogTI_kernel()
#include "kCCGE.h"
#undef TI
#define IS_ORTHOG
__global__ void
__launch_bounds__(UPDATE_THREADS_PER_BLOCK, 1)
kCalculateCoarseGridEnergyOrthog_kernel()
#include "kCCGE.h"
#undef IS_ORTHOG
__global__ void
__launch_bounds__(UPDATE_THREADS_PER_BLOCK, 1)
kCalculateCoarseGridEnergyNonOrthog_kernel()
#include "kCCGE.h"

//---------------------------------------------------------------------------------------------
// kCalculateCoarseGridEnergy_orthog: wrapper function for the eponymous kernel
//
// Arguments:
//   gpu:      overarching type for storing all parameters, coordinates, and energy function
//             terms in a simulation
//---------------------------------------------------------------------------------------------
extern "C" void kCalculateCoarseGridEnergyOrthogTI(gpuContext gpu)
{
  kCalculateCoarseGridEnergyOrthogTI_kernel<<<gpu->updateBlocks,
                                  gpu->updateThreadsPerBlock>>>();
  LAUNCHERROR("kCalculateCoarseGridEnergy");
}

//---------------------------------------------------------------------------------------------
// kCalculateCoarseGridEnergy_nonorthog: wrapper function for the eponymous kernel
//
// Arguments:
//   gpu:      overarching type for storing all parameters, coordinates, and energy function
//             terms in a simulation
//---------------------------------------------------------------------------------------------
extern "C" void kCalculateCoarseGridEnergyNonOrthogTI(gpuContext gpu)
{
  kCalculateCoarseGridEnergyNonOrthogTI_kernel<<<gpu->updateBlocks,
                                  gpu->updateThreadsPerBlock>>>();
  LAUNCHERROR("kCalculateCoarseGridEnergy");
}

//---------------------------------------------------------------------------------------------
// kCalculateCoarseGridEnergy_orthog: wrapper function for the eponymous kernel
//
// Arguments:
//   gpu:      overarching type for storing all parameters, coordinates, and energy function
//             terms in a simulation
//---------------------------------------------------------------------------------------------
extern "C" void kCalculateCoarseGridEnergyOrthog(gpuContext gpu)
{
  kCalculateCoarseGridEnergyOrthog_kernel<<<gpu->updateBlocks,
                                  gpu->updateThreadsPerBlock>>>();
  LAUNCHERROR("kCalculateCoarseGridEnergy");
}

//---------------------------------------------------------------------------------------------
// kCalculateCoarseGridEnergy_nonorthog: wrapper function for the eponymous kernel
//
// Arguments:
//   gpu:      overarching type for storing all parameters, coordinates, and energy function
//             terms in a simulation
//---------------------------------------------------------------------------------------------
extern "C" void kCalculateCoarseGridEnergyNonOrthog(gpuContext gpu)
{
  kCalculateCoarseGridEnergyNonOrthog_kernel<<<gpu->updateBlocks,
                                  gpu->updateThreadsPerBlock>>>();
  LAUNCHERROR("kCalculateCoarseGridEnergy");
}

//---------------------------------------------------------------------------------------------
// kGetEmptyStericGridVoxels_kernel: kernel for computing move generator for coarse grid
//
// Arguments:
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kGetEmptyStericGridVoxels_kernel()
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  int xdim=cSim.stericMaxXVxl/cSim.stericBitCellWidth+1;
  int ydim=cSim.stericMaxYVxl/cSim.stericBitCellHeight+1;
  int zdim=cSim.stericMaxZVxl/cSim.stericBitCellDepth+1;
  unsigned int count = 0;
  if(pos < xdim*ydim*zdim) count=cSim.pStericThreadCount[pos];
  while(pos < xdim*ydim*zdim)
  {
    unsigned int gridval = cSim.pStericGrid[pos];
    int zbyte = pos/(xdim*ydim);
    int xypos = pos - (zbyte*xdim*ydim);
    int ybyte = xypos/xdim;
    int xbyte = xypos%xdim;
    for(int z=0; z<cSim.stericBitCellDepth; z++)
    {
      for(int y=0; y<cSim.stericBitCellHeight; y++)
      {
        for(int x=0; x<cSim.stericBitCellWidth; x++)
        {
          int xvxl = xbyte*cSim.stericBitCellWidth+x;
          int yvxl = ybyte*cSim.stericBitCellHeight+y;
          int zvxl = zbyte*cSim.stericBitCellDepth+z;
          bool xtrue=false;
          bool ytrue=false;
          bool ztrue=false;
          if(cSim.pVoxelOffset[1] < cSim.pVoxelOffset[0]) 
          {
            if (xvxl > cSim.pVoxelOffset[1] || xvxl < cSim.pVoxelOffset[0])
            {
              xtrue=true;
            }
          }
          else
          {
            if(xvxl >= cSim.pVoxelOffset[0] && xvxl < cSim.pVoxelOffset[1])
            {
              xtrue=true;
            }
          }
          if(cSim.pVoxelOffset[3] < cSim.pVoxelOffset[2]) 
          {
            if (yvxl > cSim.pVoxelOffset[3] || yvxl < cSim.pVoxelOffset[2])
            {
              ytrue=true;
            }
          }
          else
          {
            if(yvxl >= cSim.pVoxelOffset[2] && yvxl < cSim.pVoxelOffset[3])
            {
              ytrue=true;
            }
          }
          if(cSim.pVoxelOffset[5] < cSim.pVoxelOffset[4]) 
          {
            if (zvxl > cSim.pVoxelOffset[5] || zvxl < cSim.pVoxelOffset[4])
            {
              ztrue=true;
            }
          }
          else
          {
            if(zvxl >= cSim.pVoxelOffset[4] && zvxl < cSim.pVoxelOffset[5])
            {
              ztrue=true;
            }
          }
          if(xtrue && ytrue && ztrue)
          {
            if((gridval & (unsigned int)(1 << (x+y*cSim.stericBitCellWidth+z*cSim.stericBitCellWidth*cSim.stericBitCellHeight))) == 0)
            {
              cSim.pStericEmptyVoxels[count]=xvxl+yvxl*cSim.stericMaxXVxl+zvxl*cSim.stericMaxXVxl*cSim.stericMaxYVxl;
              count++;
            }
          }
        }
      }
    }
    pos+=increment;
  }
}

//---------------------------------------------------------------------------------------------
// kMakeStericGrid_kernel: kernel for computing move generator for coarse grid
//
// Arguments:
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCountStericGrid_kernel()
{
  unsigned int threadid= blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int pos = threadid;
  unsigned int increment = gridDim.x * blockDim.x;
  int count=0;
  int xdim=cSim.stericMaxXVxl/cSim.stericBitCellWidth+1;
  int ydim=cSim.stericMaxYVxl/cSim.stericBitCellHeight+1;
  int zdim=cSim.stericMaxZVxl/cSim.stericBitCellDepth+1;
  while(pos < xdim*ydim*zdim)
  {
    unsigned int gridval = cSim.pStericGrid[pos];
    int zbyte = pos/(xdim*ydim);
    int xypos = pos - (zbyte*xdim*ydim);
    int ybyte = xypos/xdim;
    int xbyte = xypos%xdim;
    for(int z=0; z<cSim.stericBitCellDepth; z++)
    {
      for(int y=0; y<cSim.stericBitCellHeight; y++)
      {
        for(int x=0; x<cSim.stericBitCellWidth; x++)
        {
          int xvxl = xbyte*cSim.stericBitCellWidth+x;
          int yvxl = ybyte*cSim.stericBitCellHeight+y;
          int zvxl = zbyte*cSim.stericBitCellDepth+z;
          bool xtrue=false;
          bool ytrue=false;
          bool ztrue=false;
          if(cSim.pVoxelOffset[1] < cSim.pVoxelOffset[0]) 
          {
            if (xvxl > cSim.pVoxelOffset[1] || xvxl < cSim.pVoxelOffset[0])
            {
              xtrue=true;
            }
          }
          else
          {
            if(xvxl >= cSim.pVoxelOffset[0] && xvxl < cSim.pVoxelOffset[1])
            {
              xtrue=true;
            }
          }
          if(cSim.pVoxelOffset[3] < cSim.pVoxelOffset[2]) 
          {
            if (yvxl > cSim.pVoxelOffset[3] || yvxl < cSim.pVoxelOffset[2])
            {
              ytrue=true;
            }
          }
          else
          {
            if(yvxl >= cSim.pVoxelOffset[2] && yvxl < cSim.pVoxelOffset[3])
            {
              ytrue=true;
            }
          }
          if(cSim.pVoxelOffset[5] < cSim.pVoxelOffset[4]) 
          {
            if (zvxl > cSim.pVoxelOffset[5] || zvxl < cSim.pVoxelOffset[4])
            {
              ztrue=true;
            }
          }
          else
          {
            if(zvxl >= cSim.pVoxelOffset[4] && zvxl < cSim.pVoxelOffset[5])
            {
              ztrue=true;
            }
          }
          if(xtrue && ytrue && ztrue)
          {
            if((gridval & (unsigned int)(1 << (x+y*cSim.stericBitCellWidth+z*cSim.stericBitCellWidth*cSim.stericBitCellHeight))) == 0)
            {
              count++;
            }
          }
        }
      }
    }
    pos+=increment;
  }
  cSim.pStericThreadCount[threadid]=count;
}

//---------------------------------------------------------------------------------------------
// kMakeStericGrid_kernel: kernel for computing move generator for coarse grid
//
// Arguments:
//---------------------------------------------------------------------------------------------
__global__ void
__launch_bounds__(UPDATE_THREADS_PER_BLOCK, 1)
kMakeStericGrid_kernel()
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  int cellwidth=cSim.stericBitCellWidth;
  int cellheight=cSim.stericBitCellHeight;
  int celldepth=cSim.stericBitCellDepth;
  int xdim=cSim.stericMaxXVxl/cellwidth+1;
  int ydim=cSim.stericMaxYVxl/cellheight+1;

  while(pos < cSim.atoms)
  {
    if(round(cSim.pAtomMass[pos]) > 3)
    {
      PMEFloat radius = cSim.pStericRadius[pos];
      int img_i=cSim.pImageAtomLookup[pos];
      PMEFloat crd_x=cSim.pImageX[img_i];
      PMEFloat crd_y=cSim.pImageY[img_i];
      PMEFloat crd_z=cSim.pImageZ[img_i];
      PMEFloat fx = cSim.recipf[0][0]*crd_x + cSim.recipf[1][0]*crd_y + cSim.recipf[2][0]*crd_z;
      PMEFloat fy =                           cSim.recipf[1][1]*crd_y + cSim.recipf[2][1]*crd_z;
      PMEFloat fz =                                                     cSim.recipf[2][2]*crd_z;
      fx = fx - round(fx);
      fy = fy - round(fy);
      fz = fz - round(fz);
      if(fx < 0.0) fx+=1.0;
      if(fy < 0.0) fy+=1.0;
      if(fz < 0.0) fz+=1.0;
      int vx = fx * cSim.stericMaxXVxl;
      int vy = fy * cSim.stericMaxYVxl;
      int vz = fz * cSim.stericMaxZVxl;
      PMEFloat zradius = radius/cSim.stericGridSpacing;
      int zmin=max(0,int(floor(vz-zradius))+1);
      int zmax=min(int(vz+zradius),cSim.stericMaxZVxl);
      for(int z=zmin; z <= zmax; z++)
      {
        PMEFloat yradius = sqrt(zradius*zradius-(vz-z)*(vz-z));
        int ymin=max(0,int(floor(vy-yradius))+1);
        int ymax=min(int(vy+yradius),cSim.stericMaxYVxl);
        for(int y=ymin; y <= ymax; y++)
        {
          PMEFloat xradius = sqrt(yradius*yradius-(vy-y)*(vy-y));
          if(xradius > 0.0)
          {
            int xmin=max(0,int(floor(vx-xradius))+1);
            int xmax=min(int(vx+xradius),cSim.stericMaxXVxl);

            for(int x=xmin; x <= xmax; x++)
            {
              int bytex=x/cellwidth;
              int bytey=y/cellheight;
              int bytez=z/celldepth;
              int bitx=x%cellwidth;
              int bity=y%cellheight;
              int bitz=z%celldepth;
              int bit=bitx+bity*cellwidth+bitz*cellwidth*cellheight;
              unsigned int mask=1<<bit;
              atomicOr(&cSim.pStericGrid[bytex+bytey*xdim+bytez*xdim*ydim], mask);
            } 
          }
        }
      }
    }
    pos += increment; 
  }
}

//---------------------------------------------------------------------------------------------
// kMakeStericGrid: wrapper function for the eponymous kernel
//
// Arguments:
//   gpu:      overarching type for storing all parameters, coordinates, and energy function
//             terms in a simulation
//---------------------------------------------------------------------------------------------
extern "C" void kMakeStericGrid(gpuContext gpu)
{
  kMakeStericGrid_kernel<<<gpu->updateBlocks,
                                  gpu->updateThreadsPerBlock>>>();
  LAUNCHERROR("kMakeStericGrid");
}


//---------------------------------------------------------------------------------------------
// kCountStericGrid: wrapper function for the eponymous kernel
//
// Arguments:
//   gpu:      overarching type for storing all parameters, coordinates, and energy function
//             terms in a simulation
//---------------------------------------------------------------------------------------------
extern "C" void kCountStericGrid(gpuContext gpu)
{
  kCountStericGrid_kernel<<<gpu->PMENonbondBlocks,
                                  gpu->PMENonbondForcesThreadsPerBlock>>>();
  LAUNCHERROR("kCountStericGrid");
}

//---------------------------------------------------------------------------------------------
// kGetEmptyStericVoxels: wrapper function for the eponymous kernel
//
// Arguments:
//   gpu:      overarching type for storing all parameters, coordinates, and energy function
//             terms in a simulation
//---------------------------------------------------------------------------------------------
extern "C" void kGetEmptyStericGridVoxels(gpuContext gpu)
{
  kGetEmptyStericGridVoxels_kernel<<<gpu->PMENonbondBlocks,
                                  gpu->PMENonbondForcesThreadsPerBlock>>>();
  LAUNCHERROR("kGetEmptyStericVoxels");
}

//---------------------------------------------------------------------------------------------
// kCalculateCoarseGridEnergyInitKernels: initialize PME kernels for a simulation
//
// Arguments:
//   gpu:      overarching type for storing all parameters, coordinates, and energy function
//             terms in a simulation
//---------------------------------------------------------------------------------------------
extern "C" void kCalculateCoarseGridEnergyInitKernels(gpuContext gpu)
{
  cudaFuncSetSharedMemConfig(kCalculateCoarseGridEnergyOrthogTI_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalculateCoarseGridEnergyNonOrthogTI_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalculateCoarseGridEnergyOrthog_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalculateCoarseGridEnergyNonOrthog_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kMakeStericGrid_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCountStericGrid_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kGetEmptyStericGridVoxels_kernel, cudaSharedMemBankSizeEightByte);
}
