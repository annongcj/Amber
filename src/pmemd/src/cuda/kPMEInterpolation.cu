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
// SetkPMEInterpolationSim: setup for PME interpolation, pushing cSim (including critical PME
//                          instructions) to the GPU
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void SetkPMEInterpolationSim(gpuContext gpu)
{
  cudaError_t status;
  status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(simulationConst));
  RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

//---------------------------------------------------------------------------------------------
// GetkPMEInterpolationSim: download critical information from the GPU.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
void GetkPMEInterpolationSim(gpuContext gpu)
{
  cudaError_t status;
  status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(simulationConst));
  RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

#endif

//---------------------------------------------------------------------------------------------
// kPMEGetGridWeights$OPTION_kernel: kernels for various flavors of the standard 4-4-4 PME
//                                   interpolation
//---------------------------------------------------------------------------------------------
#define COMP_FRACTIONAL_COORD
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kPMEGetGridWeights_kernel()
#include "kReImageCoord.h"

//---------------------------------------------------------------------------------------------
#define PME_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kPMEGetGridWeightsOrthogonal_kernel()
#include "kReImageCoord.h"
#undef PME_ORTHOGONAL

//---------------------------------------------------------------------------------------------
#define PME_NTP
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kPMEGetGridWeightsNTP_kernel()
#include "kReImageCoord.h"

//---------------------------------------------------------------------------------------------
#define PME_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kPMEGetGridWeightsOrthogonalNTP_kernel()
#include "kReImageCoord.h"
#undef PME_ORTHOGONAL
#undef PME_NTP

//---------------------------------------------------------------------------------------------
#define PME_SMALLBOX
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kPMEGetGridWeightsSmall_kernel()
#include "kReImageCoord.h"

//---------------------------------------------------------------------------------------------
#define PME_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kPMEGetGridWeightsSmallOrthogonal_kernel()
#include "kReImageCoord.h"
#undef PME_ORTHOGONAL

//---------------------------------------------------------------------------------------------
#define PME_NTP
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kPMEGetGridWeightsSmallNTP_kernel()
#include "kReImageCoord.h"

//---------------------------------------------------------------------------------------------
#define PME_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kPMEGetGridWeightsSmallOrthogonalNTP_kernel()
#include "kReImageCoord.h"
#undef PME_ORTHOGONAL
#undef PME_NTP
#undef PME_SMALLBOX
#undef COMP_FRACTIONAL_COORD

//---------------------------------------------------------------------------------------------
// kPMEGetGridWeights: host function for calling any of the kernels above
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kPMEGetGridWeights(gpuContext gpu)
{
  // Set tile size in local variables
  int nBlocks = gpu->updateBlocks;
  int genThreads = gpu->updateThreadsPerBlock;

  if (gpu->bSmallBox) {
    if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
      if (gpu->sim.is_orthog) {
        kPMEGetGridWeightsSmallOrthogonalNTP_kernel<<<nBlocks, genThreads>>>();
      }
      else {
        kPMEGetGridWeightsSmallNTP_kernel<<<nBlocks, genThreads>>>();
      }
    }
    else {
      if (gpu->sim.is_orthog) {
        kPMEGetGridWeightsSmallOrthogonal_kernel<<<nBlocks, genThreads>>>();
      }
      else {
        kPMEGetGridWeightsSmall_kernel<<<nBlocks, genThreads>>>();
      }
    }
  }
  else {
    if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
      if (gpu->sim.is_orthog) {
        kPMEGetGridWeightsOrthogonalNTP_kernel<<<nBlocks, genThreads>>>();
      }
      else {
        kPMEGetGridWeightsNTP_kernel<<<nBlocks, genThreads>>>();
      }
    }
    else {
      if (gpu->sim.is_orthog) {
        kPMEGetGridWeightsOrthogonal_kernel<<<nBlocks, genThreads>>>();
      }
      else {
        kPMEGetGridWeights_kernel<<<nBlocks, genThreads>>>();
      }
    }
  }
  LAUNCHERROR("kPMEGetGridWeights");
}

//---------------------------------------------------------------------------------------------
// kPMEClearChargeGridBuffer: set the entire charge grid to zero.  Called by the host, no need
//                            to do this manually through a kernel.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kPMEClearChargeGridBuffer(gpuContext gpu)
{
#ifdef use_DPFP
  cudaMemsetAsync(gpu->sim.plliXYZ_q, 0, gpu->sim.XYZStride * sizeof(long long int));
#else
  cudaMemsetAsync(gpu->sim.plliXYZ_q, 0, gpu->sim.XYZStride * sizeof(int));
#endif
  LAUNCHERROR("kPMEClearChargeGridBuffer");
}

//---------------------------------------------------------------------------------------------
// kPMEReduceChargeGridBuffer444_kernel: The charge grid buffer (integers, specially arranged
//                                       for ease of access during mapping) must now be
//                                       transposed into the normal arrangement for the FFT.
//                                       This routine is specific to 4th order interpolation.
//
// Arguments:
//   offset:    offset for successive kernel calls
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(REDUCEFORCES_THREADS_PER_BLOCK, 1)
kPMEReduceChargeGridBuffer444_kernel(unsigned int offset)
{
  unsigned int pos = (blockIdx.x + offset) * blockDim.x + threadIdx.x;

  if (pos < cSim.nfft1xnfft2xnfft3) {
    unsigned int iz   = pos / (cSim.nfft1xnfft2);
    unsigned int iy   = (pos - iz * cSim.nfft1xnfft2) / cSim.nfft1;
    unsigned int ix   = pos - iz * cSim.nfft1xnfft2 - iy * cSim.nfft1;
    unsigned int spos = ((ix & 0x3) + (iy & 0x3)*4) + (iz << 4) +
      (((ix >> 2) << 4) + ((iy >> 2) << 2)*cSim.nfft1)*cSim.nfft3;
#ifdef use_DPFP
    long long int value = 0;
    if (spos < cSim.nfft1xnfft2xnfft3) {
      value = cSim.plliXYZ_q[spos];
    }
    PMEFloat value1  = (PMEFloat)value * ONEOVERLATTICESCALEF;
    cSim.pXYZ_q[pos] = value1;
#else
    int value = 0;
    if (spos < cSim.nfft1xnfft2xnfft3) {
      value = cSim.plliXYZ_q[spos];
    }
    PMEFloat value1  = (PMEFloat)value * ONEOVERLATTICESCALEF;
    cSim.pXYZ_q[pos] = value1;
#endif
  }
}

//---------------------------------------------------------------------------------------------
// kPMEReduceChargeGridBuffer666_kernel: The charge grid buffer (integers, specially arranged
//                                       for ease of access during mapping) must now be
//                                       transposed into the normal arrangement for the FFT.
//                                       This routine is specific to 6th order interpolation.
//
// Arguments:
//   offset:    offset for successive kernel calls
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(REDUCEFORCES_THREADS_PER_BLOCK, 1)
kPMEReduceChargeGridBuffer666_kernel(unsigned int offset)
{
  unsigned int pos = (blockIdx.x + offset) * blockDim.x + threadIdx.x;

  if (pos < cSim.nfft1xnfft2xnfft3) {
    unsigned int iz   = pos / (cSim.nfft1xnfft2);
    unsigned int iy   = (pos - iz * cSim.nfft1xnfft2) / cSim.nfft1;
    unsigned int ix   = pos - iz * cSim.nfft1xnfft2 - iy * cSim.nfft1;
    unsigned int spos = (ix & 0x3) + (iy << 2) + ((ix >> 2) << 2)*cSim.nfft2 +
                        iz*cSim.nfft1xnfft2;
#ifdef use_DPFP
    long long int value = 0;
    if (spos < cSim.nfft1xnfft2xnfft3) {
      value = cSim.plliXYZ_q[spos];
    }
    PMEFloat value1  = (PMEFloat)value * ONEOVERLATTICESCALEF;
    cSim.pXYZ_q[pos] = value1;
#else
    int value = 0;
    if (spos < cSim.nfft1xnfft2xnfft3) {
      value = cSim.plliXYZ_q[spos];
    }
    PMEFloat value1  = (PMEFloat)value * ONEOVERLATTICESCALEF;
    cSim.pXYZ_q[pos] = value1;
#endif
  }
}

//---------------------------------------------------------------------------------------------
// kPMEReorderChargeGridBuffer444_kernel: reorder the charge grid buffer in an intelligent way,
//                                        making use of all the bandwidth possible.  Reads from
//                                        the charge grid will get buffered in a small
//                                        __shared__ memory grid.  Up to 32 points will get
//                                        read at a time.  At the end, up to 32 points will get
//                                        written at a time.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(TRANSPOSE_QMESH_THREADS_PER_BLOCK, 1)
kPMEReorderChargeGridBuffer444_kernel()
{
  __shared__ PMEFloat stage[288 * (TRANSPOSE_QMESH_THREADS_PER_BLOCK >> 5)];

  // Thread indices within each warp are required
  unsigned int tgx = threadIdx.x & 31;
  unsigned int warpIdx = threadIdx.x >> 5;

  // Compute the number of sub-sections the charge grid can break into
  unsigned int nxwrite = (cSim.nfft1 >> 5) + ((cSim.nfft1 & 31) > 0);
  unsigned int nysect = cSim.nfft2 >> 2;
  unsigned int nzsect = cSim.nfft3 >> 1;

  // Determine the starting sub-section.  {xyz}sidx track positions on the NORMAL grid.
  unsigned int pos = blockIdx.x * (TRANSPOSE_QMESH_THREADS_PER_BLOCK >> 5) + warpIdx;
  unsigned int zsidx = pos / (nxwrite * nysect);
  unsigned ysidx = (pos - zsidx*(nxwrite * nysect)) / nxwrite;
  unsigned xsidx = pos - (zsidx*nysect + ysidx)*nxwrite;

  // Useful constants: meshincr gives the increment on the BUFFER grid to move
  // between successive 4x4x2 sectors to read into the __shared__ stage memory.
  // stageBaseIdx gives the index of the stage memory at which each thread will
  // start.  Incrementing by four starting from stageBaseIdx will move along
  // the stage memory as meshincr moves along the BUFFER grid.
  unsigned int meshincr = (cSim.nfft3 << 4);
  unsigned int stageBaseIdx = (tgx & 0x03) + 36*(tgx >> 2) + 288*warpIdx;
  const int levelup = (cSim.nfft1xnfft2 - 4*cSim.nfft1);

  // Loop over all sub-sections, striding by the number of thread blocks
  while (xsidx < nxwrite && ysidx < nysect && zsidx < nzsect) {

    // Read up to eight sectors from the buffer grid
    const unsigned int imin = xsidx << 3;
    const unsigned int imax = min((xsidx+1) << 3, cSim.nfft1 >> 2);
    unsigned int bufferGpos = ((imin << 4) + (ysidx << 2)*cSim.nfft1)*cSim.nfft3 +
                     (zsidx << 5) + tgx;
    unsigned int stageIdx = stageBaseIdx;

    for (unsigned int i = imin; i < imax; i++) {
#ifdef use_DPFP
      long long int value = cSim.plliXYZ_q[bufferGpos];
      stage[stageIdx] = (PMEFloat)value * ONEOVERLATTICESCALEF;
#else
      int value = cSim.plliXYZ_q[bufferGpos];
      stage[stageIdx] = (PMEFloat)value * ONEOVERLATTICESCALEF;
#endif
      bufferGpos += meshincr;
      stageIdx += 4;
    }
    __SYNCWARP(0xFFFFFFFF);

    // Make eight writes to the normally arranged charge grid
    stageIdx = tgx + 288*warpIdx;
    unsigned int qgridIdx = (xsidx << 5) + tgx + (ysidx << 2)*cSim.nfft1 +
                   (zsidx << 1)*cSim.nfft1xnfft2;
    const int maxwrite = (imax - imin)*4;
    if (maxwrite == 32) {
      for (unsigned int i = 0; i < 8; i++) {
        cSim.pXYZ_q[qgridIdx] = stage[stageIdx];
        qgridIdx += cSim.nfft1 + (i == 3)*levelup;
        stageIdx += 36;
      }
    }
    if (maxwrite < 32) {
      for (unsigned int i = 0; i < 8; i++) {
        if (tgx < maxwrite) {
          cSim.pXYZ_q[qgridIdx] = stage[stageIdx];
        }
        qgridIdx += cSim.nfft1 + (i == 3)*levelup;
        stageIdx += 36;
      }
    }

    // Increment pos and move to the next section
    pos += gridDim.x * (TRANSPOSE_QMESH_THREADS_PER_BLOCK >> 5);
    zsidx = pos / (nxwrite * nysect);
    ysidx = (pos - zsidx*(nxwrite * nysect)) / nxwrite;
    xsidx = pos - (zsidx*nysect + ysidx)*nxwrite;
  }
}

//---------------------------------------------------------------------------------------------
// kPMEReduceChargeGridBuffer: host function to call the above kernel and reduce the charge
//                             grid buffer
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kPMEReduceChargeGridBuffer(gpuContext gpu)
{
  int blocks = (gpu->sim.nfft1 * gpu->sim.nfft2 * gpu->sim.nfft3 + 127) >> 7;
  int offset = 0;
  if (gpu->sim.pmeOrder == 4) {
    kPMEReorderChargeGridBuffer444_kernel<<<256, TRANSPOSE_QMESH_THREADS_PER_BLOCK>>>();
  }
  else {
    while (blocks > 0) {
      int lblocks = min(blocks, 65535);
      kPMEReduceChargeGridBuffer666_kernel<<<lblocks, 128>>>(offset);
      LAUNCHERROR("kPMEReduceChargeGridBuffer");
      offset += 65535;
      blocks -= 65535;
    }
  }
}


//---------------------------------------------------------------------------------------------
// kPMEFillChargeGridBuffer_kernel: kernel for performing particle --> mesh mapping at 4-4-4
//                                  interpolation
//---------------------------------------------------------------------------------------------
#define PME_ORDER 4
__global__ void
__LAUNCH_BOUNDS__(P2M_THREADS_PER_BLOCK, P2M_BLOCKS_MULTIPLIER)
kPMEFillQMeshBuffer444_kernel()
#include "kCQB.h"
#define AFE_REGION1
__global__ void
__LAUNCH_BOUNDS__(P2M_THREADS_PER_BLOCK, P2M_BLOCKS_MULTIPLIER)
kPMEFillQMeshAFERegion1Buffer444_kernel()
#include "kCQB.h"
#undef AFE_REGION1
#define AFE_REGION2
__global__ void
__LAUNCH_BOUNDS__(P2M_THREADS_PER_BLOCK, P2M_BLOCKS_MULTIPLIER)
kPMEFillQMeshAFERegion2Buffer444_kernel()
#include "kCQB.h"
#undef AFE_REGION2
#undef PME_ORDER
#define PME_ORDER 6
__global__ void
__LAUNCH_BOUNDS__(P2M_THREADS_PER_BLOCK, P2M_BLOCKS_MULTIPLIER)
kPMEFillQMeshBuffer666_kernel()
#include "kCQB.h"
#define AFE_REGION1
__global__ void
__LAUNCH_BOUNDS__(P2M_THREADS_PER_BLOCK, P2M_BLOCKS_MULTIPLIER)
kPMEFillQMeshAFERegion1Buffer666_kernel()
#include "kCQB.h"
#undef AFE_REGION1
#define AFE_REGION2
__global__ void
__LAUNCH_BOUNDS__(P2M_THREADS_PER_BLOCK, P2M_BLOCKS_MULTIPLIER)
kPMEFillQMeshAFERegion2Buffer666_kernel()
#include "kCQB.h"
#undef AFE_REGION2
#undef PME_ORDER

//---------------------------------------------------------------------------------------------
// kPMEInterpolationInitKernels: initialize particle mesh Ewald interpolation kernels.  This
//                               function IS called by gpu.cpp, even though it really does
//                               nothing at the moment.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kPMEInterpolationInitKernels(gpuContext gpu)
{

}

//---------------------------------------------------------------------------------------------
// kPMEFillChargeGridBuffer: fill the charge grid buffer, perhaps after clearing it with the
//                           kPMEClearChargeGridBuffer function.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kPMEFillChargeGridBuffer(gpuContext gpu)
{
  // ti_mode == 3 corresponds to vdw only change, no difference in charges between
  // the two regions.  Skip the extra reciprocal space calculation in that case.
#if !defined(AMBER_PLATFORM_AMD)
  int nblocks = P2M_BLOCKS_MULTIPLIER * gpu->blocks;
#else
  int nblocks = (gpu->sim.atoms + P2M_BATCHSIZE - 1) / P2M_BATCHSIZE;
#endif
  if ((gpu->sim.ti_mode == 0) || (gpu->sim.ti_mode == 3)) {
    if (gpu->sim.pmeOrder == 4) {
      kPMEFillQMeshBuffer444_kernel<<<nblocks, P2M_THREADS_PER_BLOCK>>>();
    }
    else if (gpu->sim.pmeOrder == 6) {
      kPMEFillQMeshBuffer666_kernel<<<nblocks, P2M_THREADS_PER_BLOCK>>>();
    }
  }
  else {
    if (gpu->sim.AFE_recip_region == 1) {
      if (gpu->sim.pmeOrder == 4) {
        kPMEFillQMeshAFERegion1Buffer444_kernel<<<nblocks, P2M_THREADS_PER_BLOCK>>>();
      }
      else if (gpu->sim.pmeOrder == 6) {
        kPMEFillQMeshAFERegion1Buffer666_kernel<<<nblocks, P2M_THREADS_PER_BLOCK>>>();
      }
    }
    else if (gpu->sim.AFE_recip_region == 2) {
      if (gpu->sim.pmeOrder == 4) {
        kPMEFillQMeshAFERegion2Buffer444_kernel<<<nblocks, P2M_THREADS_PER_BLOCK>>>();
      }
      else if (gpu->sim.pmeOrder == 6) {
        kPMEFillQMeshAFERegion2Buffer666_kernel<<<nblocks, P2M_THREADS_PER_BLOCK>>>();
      }
    }
  }
  LAUNCHERROR("kPMEFillChargeGridBuffer");
}

//---------------------------------------------------------------------------------------------
// Kernels for computing the scalar sum in several contexts: need for a virial and alchemical
// free energy computations.
//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kPMEScalarSumRCEnergy_kernel()
#include "kPSSE.h"

//---------------------------------------------------------------------------------------------
#define MBAR
#define SC_REGION_1
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kPMEScalarSumRCEnergyAFERegion1MBAR_kernel()
#include "kPSSE.h"
#undef SC_REGION_1

//---------------------------------------------------------------------------------------------
#define SC_REGION_2
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kPMEScalarSumRCEnergyAFERegion2MBAR_kernel()
#include "kPSSE.h"
#undef SC_REGION_2
#undef MBAR

//---------------------------------------------------------------------------------------------
#define SC_REGION_1
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kPMEScalarSumRCEnergyAFERegion1_kernel()
#include "kPSSE.h"
#undef SC_REGION_1

//---------------------------------------------------------------------------------------------
#define SC_REGION_2
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kPMEScalarSumRCEnergyAFERegion2_kernel()
#include "kPSSE.h"
#undef SC_REGION_2
#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kPMEScalarSumRC_kernel()
#include "kPSSE.h"

//---------------------------------------------------------------------------------------------
#define MBAR
#define SC_REGION_1
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kPMEScalarSumRCAFERegion1MBAR_kernel()
#include "kPSSE.h"
#undef SC_REGION_1

//---------------------------------------------------------------------------------------------
#define SC_REGION_2
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kPMEScalarSumRCAFERegion2MBAR_kernel()
#include "kPSSE.h"
#undef SC_REGION_2
#undef MBAR

//---------------------------------------------------------------------------------------------
#define SC_REGION_1
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kPMEScalarSumRCAFERegion1_kernel()
#include "kPSSE.h"
#undef SC_REGION_1

//---------------------------------------------------------------------------------------------
#define SC_REGION_2
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kPMEScalarSumRCAFERegion2_kernel()
#include "kPSSE.h"
#undef SC_REGION_2

//---------------------------------------------------------------------------------------------
#define PME_VIRIAL
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kPMEScalarSumRCEnergyVirial_kernel(PMEFloat pi_vol_inv)
#include "kPSSE.h"
#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kPMEScalarSumRCVirial_kernel(PMEFloat pi_vol_inv)
#include "kPSSE.h"
#undef PME_VIRIAL

//---------------------------------------------------------------------------------------------
// kPMEScalarSumRC: launch kernels to perform the "scalar sum," the element-wise multiplication
//                  of the charge grid with the B-spline adjusted erf(r)/r influence function
//                  when each is expressed in reciprocal space.  This will only compute
//                  quantities needed for force computations; the energy itself is an extra
//                  computation handled in a special case by the next function.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//   vol: periodic box volume
//---------------------------------------------------------------------------------------------
extern "C" void kPMEScalarSumRC(gpuContext gpu, PMEDouble vol)
{
  int nblocks = gpu->blocks;
  int nthreads = gpu->generalThreadsPerBlock;

  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
    kPMEScalarSumRCVirial_kernel<<<nblocks, nthreads>>>((PMEFloat)(1.0 / (PI*vol)));
  }
  else {

    // The following conditions denote "not alchemical free energy, or doing vdw only change"
    if ((gpu->sim.ti_mode == 0) || (gpu->sim.ti_mode == 3)) {
      kPMEScalarSumRC_kernel<<<nblocks, nthreads>>>();
    }
    else {
      if (gpu->sim.AFE_recip_region == 1) {
        if(gpu->sim.ifmbar > 0)
          kPMEScalarSumRCAFERegion1MBAR_kernel<<<nblocks, nthreads>>>();
        else
          kPMEScalarSumRCAFERegion1_kernel<<<nblocks, nthreads>>>();
      }
      else {
        if(gpu->sim.ifmbar > 0)
          kPMEScalarSumRCAFERegion2MBAR_kernel<<<nblocks, nthreads>>>();
        else
          kPMEScalarSumRCAFERegion2_kernel<<<nblocks, nthreads>>>();
      }
    }
  }
  LAUNCHERROR("kPMEScalarSumRC");
}

//---------------------------------------------------------------------------------------------
// kPMEScalarSumRCEnergy: launch kernels to perform the "scalar sum" and compute the associated
//                        energy.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//   vol: periodic box volume
//---------------------------------------------------------------------------------------------
extern "C" void kPMEScalarSumRCEnergy(gpuContext gpu, PMEDouble vol)
{
  int nblocks = gpu->blocks;
  int nthreads = gpu->generalThreadsPerBlock;

  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
    kPMEScalarSumRCEnergyVirial_kernel<<<nblocks, nthreads>>>((PMEFloat)(1.0 / (PI*vol)));
  }
  else {

    // The following conditions denote "not alchemical free energy, or doing vdw only change"
    if ((gpu->sim.ti_mode == 0) || (gpu->sim.ti_mode == 3)) {
      kPMEScalarSumRCEnergy_kernel<<<nblocks, nthreads>>>();
    }
    else {
      if (gpu->sim.AFE_recip_region == 1) {
        if(gpu->sim.ifmbar > 0)
          kPMEScalarSumRCEnergyAFERegion1MBAR_kernel<<<nblocks, nthreads>>>();
        else
          kPMEScalarSumRCEnergyAFERegion1_kernel<<<nblocks, nthreads>>>();
      }
      else {
        if(gpu->sim.ifmbar > 0)
          kPMEScalarSumRCEnergyAFERegion2MBAR_kernel<<<nblocks, nthreads>>>();
        else
          kPMEScalarSumRCEnergyAFERegion2_kernel<<<nblocks, nthreads>>>();
      }
    }
  }
  LAUNCHERROR("kPMEScalarSumRCEnergy");
}

//---------------------------------------------------------------------------------------------
// Kernels for computing the "gradient sum," that is mesh --> particle interpolation
// ("mapping") to complete the computation of particle forces based on mesh charges.  The
// various flavors mirror the kernels for scalar sum computation.
//---------------------------------------------------------------------------------------------

static const int GRADSUMTHREADS  = 64;
static const int GRADSUMLOADSIZE = 24;

//---------------------------------------------------------------------------------------------
__global__ void
#if !defined(AMBER_PLATFORM_AMD)
#ifdef use_DPFP
__LAUNCH_BOUNDS__(GRADSUMTHREADS, 8)
#else
__LAUNCH_BOUNDS__(GRADSUMTHREADS, 16)
#endif
#else
__LAUNCH_BOUNDS__(GRADSUMTHREADS, 1)
#endif
kPMEGradSum64_kernel(unsigned int offset)
#include "kPGS.h"

//---------------------------------------------------------------------------------------------
#define SC_REGION_1
__global__ void
#if !defined(AMBER_PLATFORM_AMD)
#ifdef use_DPFP
__LAUNCH_BOUNDS__(GRADSUMTHREADS, 8)
#else
__LAUNCH_BOUNDS__(GRADSUMTHREADS, 16)
#endif
#else
__LAUNCH_BOUNDS__(GRADSUMTHREADS, 1)
#endif
kPMEGradSum64AFERegion1_kernel(unsigned int offset)
#include "kPGS.h"
#undef SC_REGION_1

//---------------------------------------------------------------------------------------------
#define SC_REGION_2
__global__ void
#if !defined(AMBER_PLATFORM_AMD)
#ifdef use_DPFP
__LAUNCH_BOUNDS__(GRADSUMTHREADS, 8)
#else
__LAUNCH_BOUNDS__(GRADSUMTHREADS, 16)
#endif
#else
__LAUNCH_BOUNDS__(GRADSUMTHREADS, 1)
#endif
kPMEGradSum64AFERegion2_kernel(unsigned int offset)
#include "kPGS.h"
#undef SC_REGION_2

//---------------------------------------------------------------------------------------------
#define PME_VIRIAL
__global__ void
#if !defined(AMBER_PLATFORM_AMD)
#ifdef use_DPFP
__LAUNCH_BOUNDS__(GRADSUMTHREADS, 8)
#else
__LAUNCH_BOUNDS__(GRADSUMTHREADS, 8)
#endif
#else
__LAUNCH_BOUNDS__(GRADSUMTHREADS, 1)
#endif
kPMEGradSum64Virial_kernel(unsigned int offset)
#include "kPGS.h"
#undef PME_VIRIAL

//---------------------------------------------------------------------------------------------
#define PME_ORDER 4
__global__ void
__LAUNCH_BOUNDS__(128, 4)
kPMEGradientSum444_kernel()
#include "kPBGS.h"

#define SC_REGION_1
__global__ void
__LAUNCH_BOUNDS__(128, 4)
kPMEGradientSumAFERegion1_444_kernel()
#include "kPBGS.h"
#undef SC_REGION_1

#define SC_REGION_2
__global__ void
__LAUNCH_BOUNDS__(128, 4)
kPMEGradientSumAFERegion2_444_kernel()
#include "kPBGS.h"
#undef SC_REGION_2

#define PME_VIRIAL
__global__ void
__LAUNCH_BOUNDS__(128, 4)
kPMEGradientSumVirial444_kernel()
#include "kPBGS.h"
#undef PME_VIRIAL
#undef PME_ORDER

#define PME_ORDER 6
__global__ void
__LAUNCH_BOUNDS__(128, 4)
kPMEGradientSum666_kernel()
#include "kPBGS.h"

#define SC_REGION_1
__global__ void
__LAUNCH_BOUNDS__(128, 4)
kPMEGradientSumAFERegion1_666_kernel()
#include "kPBGS.h"
#undef SC_REGION_1

#define SC_REGION_2
__global__ void
__LAUNCH_BOUNDS__(128, 4)
kPMEGradientSumAFERegion2_666_kernel()
#include "kPBGS.h"
#undef SC_REGION_2

#define PME_VIRIAL
__global__ void
__LAUNCH_BOUNDS__(128, 4)
kPMEGradientSumVirial666_kernel()
#include "kPBGS.h"
#undef PME_VIRIAL
#undef PME_ORDER

//---------------------------------------------------------------------------------------------
// Kernel for computing grid titration forces
//---------------------------------------------------------------------------------------------
__global__
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
void kPMEGradSumPHMD_kernel()
#include "kPGSPHMD.h"

//---------------------------------------------------------------------------------------------
// kPMEGradSum: host function to launch the appropriate mesh --> particle interpolation kernel
//
// Arguments:
//---------------------------------------------------------------------------------------------
extern "C" void kPMEGradSum(gpuContext gpu)
{
  int blocks = (gpu->sim.atoms + GRADSUMLOADSIZE - 1) / GRADSUMLOADSIZE;
  int offset = 0;
  while (blocks > 0) {
    int lblocks = min(blocks, 65535);
    if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
      kPMEGradSum64Virial_kernel<<<lblocks, GRADSUMTHREADS>>>(offset);
    }
    else {

      // The following conditions denote "not alchemical free energy, or doing vdw only change"
      if ((gpu->sim.ti_mode == 0) || (gpu->sim.ti_mode == 3)) {
        kPMEGradSum64_kernel<<<lblocks, GRADSUMTHREADS>>>(offset);
      }
      else {

        // Need to call once for each region
        if (gpu->sim.AFE_recip_region == 1) {
          kPMEGradSum64AFERegion1_kernel<<<lblocks, GRADSUMTHREADS>>>(offset);
        }
        else {
          kPMEGradSum64AFERegion2_kernel<<<lblocks, GRADSUMTHREADS>>>(offset);
        }
      }
    }
    LAUNCHERROR("kPMEGradSum");
    blocks -= 65535;
    offset += 65535;
  }
  if (gpu->sim.iphmd == 3) {
    kPMEGradSumPHMD_kernel<<<gpu->PMENonbondBlocks, gpu->PMENonbondForcesThreadsPerBlock>>>();
  }
}

//---------------------------------------------------------------------------------------------
// kPMEForwardFFT: perform the forward 3D FFT, taking the charge grid into reciprocal space.
//---------------------------------------------------------------------------------------------
extern "C" void kPMEForwardFFT(gpuContext gpu)
{
  PRINTMETHOD("kPMEForwardFFT");
#ifdef VKFFT
  VkFFTAppend(gpu->appR2C, -1, NULL);
#else
#  ifdef use_DPFP
  cufftExecD2Z(gpu->forwardPlan, gpu->sim.pXYZ_q, gpu->sim.pXYZ_qt);
#  else
  cufftExecR2C(gpu->forwardPlan, gpu->sim.pXYZ_q, gpu->sim.pXYZ_qt);
#  endif
#endif
}

//---------------------------------------------------------------------------------------------
// kPMEBackwardFFT: perform the backward 3D FFT, obtaining the electrostatic potential in
//                  real space.
//---------------------------------------------------------------------------------------------
extern "C" void kPMEBackwardFFT(gpuContext gpu)
{
  PRINTMETHOD("kPMEBackwardFFT");
#ifdef VKFFT
  VkFFTAppend(gpu->appR2C, 1, NULL);
#else
#  ifdef use_DPFP
  cufftExecZ2D(gpu->backwardPlan, gpu->sim.pXYZ_qt, gpu->sim.pXYZ_q);
#  else
  cufftExecC2R(gpu->backwardPlan, gpu->sim.pXYZ_qt, gpu->sim.pXYZ_q);
#  endif
#endif
}
