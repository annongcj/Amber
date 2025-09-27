#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------
#include "gpu.h"
#include "matrix.h"
#include "ptxmacros.h"
#ifdef AMBER_PLATFORM_AMD
#  ifdef MPI
#    include <hip/hip_cooperative_groups.h>
#  endif
#else
#  include <cuda.h>
#endif


// Use global instance instead of a local copy
#include "simulationConst.h"
CSIM_STO simulationConst cSim;

#if !defined(__HIPCC_RDC__)

//---------------------------------------------------------------------------------------------
// SetkForcesUpdateSim: the first function called by gpuCopyConstants.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void SetkForcesUpdateSim(gpuContext gpu)
{
  cudaError_t status;
  status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(simulationConst));
  RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

//---------------------------------------------------------------------------------------------
// GetkForcesUpdateSim: download information about force computations from the device.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
void GetkForcesUpdateSim(gpuContext gpu)
{
  cudaError_t status;
  status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(simulationConst));
  RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

#endif

//---------------------------------------------------------------------------------------------
// kClearNBForces_kernel: set non-bonded interaction forces on all particle accumulators in
//                        the system to zero.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(CLEARFORCES_THREADS_PER_BLOCK, 1)
kClearNBForces_kernel()
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
#if !defined(AMBER_PLATFORM_AMD)
#pragma unroll 16
#endif
  while (pos < cSim.stride3) {
    cSim.pNBForceAccumulator[pos] = (PMEAccumulator)0;
    pos += blockDim.x * gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kClearNBForces: launch the kernel above to zero out forces.  This is only launched in MPI
//                 mode under constant pressure conditions, at the end of neighbor list
//                 building.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kClearNBForces(gpuContext gpu)
{
  kClearNBForces_kernel<<<gpu->blocks, gpu->NLClearForcesThreadsPerBlock>>>();
  LAUNCHERROR("kClearNBForces");
}
//---------------------------------------------------------------------------------------------
// Kernel for (re)initializing titration forces
//
//---------------------------------------------------------------------------------------------
__global__ void
kClearForces_PHMD_kernel()
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
#pragma unroll 16
  while (pos < cSim.stride) {
    cSim.pdph_accumulator[pos] = 0;
    cSim.pdph_plus_accumulator[pos] = 0;
    pos += blockDim.x * gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// Kernel for clearing titration forces
//
//---------------------------------------------------------------------------------------------
__global__ void
kClearDerivs_PHMD_kernel()
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  cSim.pdph_theta[pos] = 0;
}

//---------------------------------------------------------------------------------------------
// Kernel for model titration forces
//
//---------------------------------------------------------------------------------------------
__global__ void
kModelForces_PHMD_kernel()
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  int ntitr = cSim.ntitr;
  int psp_grp = cSim.psp_grp[pos];
  PMEFloat lambda,  x1, x2, r1, r2, r3, r4, r5, r6, dubarr, duph, dumod;
  PMEFloat dph_theta, dph_theta_minus, ppara, ppark, pparkminus, pparmod4;
  PMEFloat phmd_ene, umod, ubar, uph;
  PMEFloat2 psp_par;
  phmd_ene = 0.0;
  if( pos < ntitr ) {
    if( psp_grp == 0 ) {
      lambda = sin(cSim.pph_theta[pos]);
      lambda *= lambda;
      dph_theta = cSim.ppark[pos] + (0.5 - lambda) * 8.0 * cSim.pbarr[pos] -
                  2.0 * cSim.ppara[pos] * ( lambda - cSim.pparb[pos] );
      phmd_ene = cSim.ppark[pos] * lambda - 
                 cSim.ppara[pos] * ( lambda - cSim.pparb[pos] ) *
                 ( lambda - cSim.pparb[pos] ) -
                 4.0 * cSim.pbarr[pos] * ( 0.5 - lambda ) * ( 0.5 - lambda );
      if(cSim.use_forces) {
        cSim.pdph_theta[pos] += dph_theta;
      }
    }
    else if( psp_grp == 2 ) {
      lambda = sin(cSim.pph_theta[pos-1]);
      lambda *= lambda;
      x1 = sin(cSim.pph_theta[pos]);
      x1 *= x1;
      x2 = 1.0 - x1;
      ppara = cSim.ppara[pos];
      psp_par = cSim.psp_par[pos];
      ppark = cSim.ppark[pos];
      pparkminus = cSim.ppark[pos-1];
      r2 = -2.0 * ppara * cSim.pparb[pos];
      r1 = -2.0 * cSim.ppara[pos-1] * cSim.pparb[pos-1] - r2;
      r3 = -2.0 * psp_par.x * psp_par.y - r1;
      dubarr = 8.0 * cSim.pbarr[pos-1] * ( lambda - 0.5 );
      ubar = 4.0 * cSim.pbarr[pos-1] * ( lambda - 0.5 ) * ( lambda - 0.5 ) +
             4.0 * cSim.pbarr[pos] * ( x1 - 0.5 ) * ( x1 - 0.5 );
      uph = lambda * (pparkminus * x1 + ppark * x2);
      umod = psp_par.x * lambda * lambda * x1 * x1 + r1 * lambda * x1 +
             r2 * lambda + r3 * lambda * lambda * x1 + ppara * lambda * lambda;
      duph = pparkminus * x1 + ppark * x2;
      dumod = 2.0 * psp_par.x * lambda * x1 * x1 + r1 * x1 +
              r2 + 2.0 * r3 * lambda * x1 + 2.0 * ppara * lambda;
      dph_theta_minus = duph - dumod - dubarr;
      phmd_ene = uph - umod - ubar;
      dubarr = 8.0 * cSim.pbarr[pos] * ( x1 - 0.5 );
      duph = lambda * ( pparkminus - ppark );
      dumod = 2.0 * psp_par.x * lambda * lambda * x1 + r1 * lambda +
              r3 * lambda * lambda;
      dph_theta = duph - dumod - dubarr;

      if( cSim.use_forces ) {
        cSim.pdph_theta[pos] += dph_theta;
        cSim.pdph_theta[pos-1] += dph_theta_minus;
      }
    }
    else if( psp_grp == 4 ) {
      lambda = sin(cSim.pph_theta[pos-1]);
      lambda *= lambda;
      x1 = sin(cSim.pph_theta[pos]);
      x1 *= x1;
      x2 = 1.0 - x1;
      pparmod4 = cSim.pparmod[6*pos+4];
      r1 = cSim.pparmod[6*pos];
      r2 = cSim.pparmod[6*pos+1];
      r3 = cSim.pparmod[6*pos+2];
      r4 = cSim.pparmod[6*pos+3];
      r5 = pparmod4 - r1 * r4 * r4;
      r6 = -2.0 * pparmod4 * cSim.pparmod[6*pos+5] - r2 * r4 * r4;
      ubar = 4.0 * cSim.pbarr[pos] * (x1 - 0.5) * (x1 - 0.5) +
             4.0 * cSim.pbarr[pos-1] * (lambda - 0.5) * (lambda - 0.5);
      uph = lambda * cSim.ppark[pos-1];
      umod = (r1 * lambda * lambda + r2 * lambda + r3) * (x1 - r4) * (x1 - r4) +
             r5 * lambda * lambda + r6 * lambda;
      dubarr = 8.0 * cSim.pbarr[pos-1] * ( lambda - 0.5 );
      duph = cSim.ppark[pos-1];
      dumod = ( 2.0 * r1 * lambda + r2 ) * ( x1 - r4 ) * ( x1 - r4 ) +
              2.0 * r5 * lambda + r6;
      dph_theta_minus = duph - dumod - dubarr;
      phmd_ene = uph - umod - ubar;
      dubarr = 8.0 * cSim.pbarr[pos] * ( x1 - 0.5 );
      dumod = 2.0 * ( r1 * lambda * lambda + r2 * lambda + r3 ) * ( x1 - r4 );
      dph_theta = - dumod - dubarr;
      if( cSim.use_forces ) {
        cSim.pdph_theta[pos] += dph_theta;
        cSim.pdph_theta[pos-1] += dph_theta_minus;
      }
    }
    //atomicAdd(cSim.pphmd_ene, llitoulli(llrint(phmd_ene * ENERGYSCALE)));
  }
}
//---------------------------------------------------------------------------------------------
// Kernels for (re)initializing force accumulators.  The distinctions rest on whether there is
// a neighbor list in use and whether an alchemical free energy simulation is happening.
//
// Neighborlist-oriented kernels take the argument nonbondWarps, the number of warps devoted
// to non-bonded direct space interaction computations.  Each of these warps gets its own
// accumulators.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(CLEARFORCES_THREADS_PER_BLOCK, 1)
kClearForces_kernel(int bondWorkBlocks)
{
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

  // Clear GB NB kernel counters
  if (pos == 0) {
    cSim.pFrcBlkCounters[1] = bondWorkBlocks;
  }
  if (pos < 3) {
    cSim.pGBBRPosition[pos] = cSim.GBTotalWarps[pos];
  }
  if (pos < cSim.EnergyTerms) {
    cSim.pEnergyBuffer[pos] = (unsigned long long int)0;
  }
#if !defined(AMBER_PLATFORM_AMD)
#pragma unroll 16
#endif
  while (pos < cSim.forceBufferStride) {
    cSim.pForceAccumulator[pos] = (PMEAccumulator)0;
    pos += blockDim.x * gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(CLEARFORCES_THREADS_PER_BLOCK, 1)
kNLClearForces_kernel(int nonbondWarps, int bondWorkBlocks)
{
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  if (pos == 0) {
    cSim.pFrcBlkCounters[0] = nonbondWarps;
    cSim.pFrcBlkCounters[1] = bondWorkBlocks;
  }
  if (pos < cSim.EnergyTerms) {
    cSim.pEnergyBuffer[pos] = 0;
  }
#if !defined(AMBER_PLATFORM_AMD)
#pragma unroll 16
#endif
  while (pos < cSim.forceBufferStride) {
    cSim.pForceAccumulator[pos] = (PMEAccumulator)0;
    pos += blockDim.x * gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(CLEARFORCES_THREADS_PER_BLOCK, 1)
kAFEClearForces_kernel(int bondWorkBlocks)
{
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

  // Clear GB NB kernel counters
  if (pos == 0) {
    cSim.pFrcBlkCounters[1] = bondWorkBlocks;
  }
  if (pos < 3) {
    cSim.pGBBRPosition[pos] = cSim.GBTotalWarps[pos];
  }
  if (pos < cSim.EnergyTerms) {
    cSim.pEnergyBuffer[pos] = 0;
  }
  if (pos < cSim.AFETerms) {
    cSim.pAFEBuffer[pos] = 0;
  }
#if !defined(AMBER_PLATFORM_AMD)
#pragma unroll 16
#endif
  while (pos < cSim.forceBufferStride) {
    cSim.pForceAccumulator[pos] = (PMEAccumulator)0;
    pos += blockDim.x * gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(CLEARFORCES_THREADS_PER_BLOCK, 1)
kNLAFEClearForces_kernel(int nonbondWarps, int bondWorkBlocks)
{
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  if (pos == 0) {
    cSim.pFrcBlkCounters[0] = nonbondWarps;
    cSim.pFrcBlkCounters[1] = bondWorkBlocks;
  }
  if (pos < cSim.EnergyTerms) {
    cSim.pEnergyBuffer[pos] = 0;
  }
  if (pos < cSim.AFETerms) {
    cSim.pAFEBuffer[pos] = 0;
  }
#if !defined(AMBER_PLATFORM_AMD)
#pragma unroll 16
#endif
  // Will need to update this if we add the virial for berendsen barostat
  while (pos < cSim.forceBufferStride) {
    cSim.pForceAccumulator[pos] = (PMEAccumulator)0;
    pos += blockDim.x * gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(CLEARFORCES_THREADS_PER_BLOCK, 1)
kAFEMBARClearForces_kernel(int bondWorkBlocks)
{
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

  // Clear GB NB kernel counters
  if (pos == 0) {
    cSim.pFrcBlkCounters[1] = bondWorkBlocks;
  }
  if (pos < 3) {
    cSim.pGBBRPosition[pos] = cSim.GBTotalWarps[pos];
  }
  if (pos < cSim.EnergyTerms) {
    cSim.pEnergyBuffer[pos] = 0;
  }
  if (pos < cSim.AFETerms) {
    cSim.pAFEBuffer[pos] = 0;
  }
  if (pos < cSim.bar_states) {
    cSim.pBarTot[pos] = 0;
  }

#if !defined(AMBER_PLATFORM_AMD)
#pragma unroll 16
#endif
  while (pos < cSim.forceBufferStride) {
    cSim.pForceAccumulator[pos]  = (PMEAccumulator)0;
    pos                         += blockDim.x * gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(CLEARFORCES_THREADS_PER_BLOCK, 1)
kNLAFEMBARClearForces_kernel(int nonbondWarps, int bondWorkBlocks)
{
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  if (pos == 0) {
    cSim.pFrcBlkCounters[0] = nonbondWarps;
    cSim.pFrcBlkCounters[1] = bondWorkBlocks;
  }
  if (pos < cSim.EnergyTerms) {
    cSim.pEnergyBuffer[pos] = 0;
  }
  if (pos < cSim.AFETerms) {
    cSim.pAFEBuffer[pos] = 0;
  }
  if (pos < cSim.bar_states) {
    cSim.pBarTot[pos] = 0;
  }

#if !defined(AMBER_PLATFORM_AMD)
#pragma unroll 16
#endif
  // Will need to update this if we add the virial for berendsen barostat
  while (pos < cSim.forceBufferStride) {
    cSim.pForceAccumulator[pos]  = (PMEAccumulator)0;
    pos                         += blockDim.x * gridDim.x;
  }
}
//---------------------------------------------------------------------------------------------
// kClearDerivs_PHMD: host function to launch the appropriate kernel for clearing titration
//                    forces.
//
// Arguments:
//   gpu:          overarching type for storing all parameters, coordinates, and the energy
//                 function
//---------------------------------------------------------------------------------------------
void kClearDerivs_PHMD(gpuContext gpu)
{
  kClearDerivs_PHMD_kernel<<<1, gpu->sim.phmd_padded_ntitr>>>();
}
//---------------------------------------------------------------------------------------------
// kModelForces_PHMD: host function to launch the appropriate kernel for clearing titration
//                    forces.
//
// Arguments:
//   gpu:          overarching type for storing all parameters, coordinates, and the energy
//                 function
//---------------------------------------------------------------------------------------------
void kModelForces_PHMD(gpuContext gpu)
{
  kModelForces_PHMD_kernel<<<1, gpu->sim.phmd_padded_ntitr>>>();
}

//---------------------------------------------------------------------------------------------
// kClearForces: host function to launch the appropriate kernels for clearing forces.  The
//               choice dpeends on whether there is a neighborlist in use (and thus whether
//               the representation of atoms in the hash table is relevant).
//
// Arguments:
//   gpu:          overarching type for storing all parameters, coordinates, and the energy
//                 function
//   nonbondWarps: the number of warps performing non-bonded calculations (each warp has its
//                 own accumulators)
//---------------------------------------------------------------------------------------------
void kClearForces(gpuContext gpu, int nonbondWarps)
{
  if (gpu->bNeighborList) {
    if (gpu->sim.ti_mode == 0) {
      kNLClearForces_kernel<<<gpu->blocks,
                              gpu->clearForcesThreadsPerBlock>>>(nonbondWarps,
                                                                 gpu->bondWorkBlocks);
    }
    else {
      if (gpu->sim.ifmbar > 0) {
        kNLAFEMBARClearForces_kernel<<<gpu->blocks,
                                       gpu->clearForcesThreadsPerBlock>>>(nonbondWarps,
                                                                          gpu->bondWorkBlocks);
      }
      else {
        kNLAFEClearForces_kernel<<<gpu->blocks,
                                   gpu->clearForcesThreadsPerBlock>>>(nonbondWarps,
                                                                      gpu->bondWorkBlocks);
      }
    }
    if (gpu->iphmd > 0) {
      kClearForces_PHMD_kernel<<<gpu->blocks,
                                 gpu->clearForcesThreadsPerBlock>>>();
    }
  }
  else {
#ifdef AMBER_PLATFORM_AMD
    if (gpu->sim.ti_mode != 0) {
      if (gpu->sim.ifmbar > 0) {
        cudaMemsetAsync(gpu->sim.pBarTot, 0,
                        gpu->sim.bar_states * sizeof(unsigned long long int));
      }
      cudaMemsetAsync(gpu->sim.pAFEBuffer, 0,
                      gpu->sim.AFETerms * sizeof(unsigned long long int));
    }
    cudaMemsetAsync(gpu->sim.pEnergyBuffer, 0,
                    gpu->sim.EnergyTerms * sizeof(unsigned long long int));
    cudaMemsetAsync(gpu->sim.pForceAccumulator, 0,
                    gpu->sim.stride3 * gpu->sim.nonbondForceBuffers * sizeof(PMEAccumulator));
#else
    if (gpu->sim.ti_mode == 0) {
      kClearForces_kernel<<<gpu->blocks,
                            gpu->clearForcesThreadsPerBlock>>>(gpu->bondWorkBlocks);
    }
    else {
      if(gpu->sim.ifmbar > 0)
        kAFEMBARClearForces_kernel<<<gpu->blocks,
                                     gpu->clearForcesThreadsPerBlock>>>(gpu->bondWorkBlocks);
      else
        kAFEClearForces_kernel<<<gpu->blocks,
                                 gpu->clearForcesThreadsPerBlock>>>(gpu->bondWorkBlocks);
    }
    if (gpu->iphmd > 0) {
      kClearForces_PHMD_kernel<<<gpu->blocks,
                                 gpu->clearForcesThreadsPerBlock>>>();
    }
#endif
  }
  LAUNCHERROR("kClearForces");
}

//---------------------------------------------------------------------------------------------
// boltz2: one half the gas constant, in units of kcal/mol.
//---------------------------------------------------------------------------------------------
static __constant__ PMEDouble boltz2 = 0.00831441 * 0.5 / 4.184;

//---------------------------------------------------------------------------------------------
// COM: struct to hold bounds force center of mass computation
//---------------------------------------------------------------------------------------------
struct COM {
  PMEFloat xmin;
  PMEFloat ymin;
  PMEFloat zmin;
  PMEFloat xmax;
  PMEFloat ymax;
  PMEFloat zmax;
};

//---------------------------------------------------------------------------------------------
// kRecenter_Molecule1_kernel: compute the bounds on atoms in the system.  The comments below
//                             talk about SUMS, but there doesn't appear to be any addition
//                             going on.  Are the comments inside this function accurate?  This
//                             is only called by kRecenter_Molecule below.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kRecenter_Molecule1_kernel()
{
  __shared__ COM sA[GENERAL_THREADS_PER_BLOCK];
  PMEFloat xmin = (PMEFloat)999999999999.0;
  PMEFloat ymin = (PMEFloat)999999999999.0;
  PMEFloat zmin = (PMEFloat)999999999999.0;
  PMEFloat xmax = (PMEFloat)-999999999999.0;
  PMEFloat ymax = (PMEFloat)-999999999999.0;
  PMEFloat zmax = (PMEFloat)-999999999999.0;
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;

  // Perform individual sums
  while (pos < cSim.atoms) {
    PMEFloat2 xy = cSim.pAtomXYSP[pos];
    PMEFloat z   = cSim.pAtomZSP[pos];
    xmax         = max(xy.x, xmax);
    xmin         = min(xy.x, xmin);
    ymax         = max(xy.y, ymax);
    ymin         = min(xy.y, ymin);
    zmax         = max(z,    zmax);
    zmin         = min(z,    zmin);
    pos         += blockDim.x * gridDim.x;
  }

  // Perform local reduction to thread 0
  sA[threadIdx.x].xmin = xmin;
  sA[threadIdx.x].ymin = ymin;
  sA[threadIdx.x].zmin = zmin;
  sA[threadIdx.x].xmax = xmax;
  sA[threadIdx.x].ymax = ymax;
  sA[threadIdx.x].zmax = zmax;
  __syncthreads();
  for (uint stride = blockDim.x >> 1; stride > 0; stride >>= 1) {
    if (threadIdx.x < stride) {
      sA[threadIdx.x].xmin = min(sA[threadIdx.x].xmin ,sA[threadIdx.x + stride].xmin);
      sA[threadIdx.x].ymin = min(sA[threadIdx.x].ymin ,sA[threadIdx.x + stride].ymin);
      sA[threadIdx.x].zmin = min(sA[threadIdx.x].zmin ,sA[threadIdx.x + stride].zmin);
      sA[threadIdx.x].xmax = max(sA[threadIdx.x].xmax ,sA[threadIdx.x + stride].xmax);
      sA[threadIdx.x].ymax = max(sA[threadIdx.x].ymax ,sA[threadIdx.x + stride].ymax);
      sA[threadIdx.x].zmax = max(sA[threadIdx.x].zmax ,sA[threadIdx.x + stride].zmax);
    }
    __syncthreads();
  }

  // Output sum if thread 0
  if (threadIdx.x == 0) {
    cSim.pXMin[blockIdx.x] = sA[threadIdx.x].xmin;
    cSim.pYMin[blockIdx.x] = sA[threadIdx.x].ymin;
    cSim.pZMin[blockIdx.x] = sA[threadIdx.x].zmin;
    cSim.pXMax[blockIdx.x] = sA[threadIdx.x].xmax;
    cSim.pYMax[blockIdx.x] = sA[threadIdx.x].ymax;
    cSim.pZMax[blockIdx.x] = sA[threadIdx.x].zmax;
  }
}

//---------------------------------------------------------------------------------------------
// kRecenter_Molecule2_kernel: this appears to translate the system to put the midpoint between
//                             its extreme coordinates in x, y, and z at the origin.  Like
//                             kRecenter_Molecule1_kernel above, it is only called by
//                             kRecenter_Molecule.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kRecenter_Molecule2_kernel()
{
  __shared__ COM sA[GENERAL_THREADS_PER_BLOCK];
  // Read in local offsets
  unsigned int pos = threadIdx.x;
  while (pos < gridDim.x) {
    sA[pos].xmin = cSim.pXMin[pos];
    sA[pos].ymin = cSim.pYMin[pos];
    sA[pos].zmin = cSim.pZMin[pos];
    sA[pos].xmax = cSim.pXMax[pos];
    sA[pos].ymax = cSim.pYMax[pos];
    sA[pos].zmax = cSim.pZMax[pos];
    pos += blockDim.x;
  }
  __syncthreads();

 // Perform local reduction to thread 0
 for (uint stride = blockDim.x >> 1; stride > 0; stride >>= 1) {
   if (threadIdx.x < stride) {
     sA[threadIdx.x].xmin = min(sA[threadIdx.x].xmin ,sA[threadIdx.x + stride].xmin);
     sA[threadIdx.x].ymin = min(sA[threadIdx.x].ymin ,sA[threadIdx.x + stride].ymin);
     sA[threadIdx.x].zmin = min(sA[threadIdx.x].zmin ,sA[threadIdx.x + stride].zmin);
     sA[threadIdx.x].xmax = max(sA[threadIdx.x].xmax ,sA[threadIdx.x + stride].xmax);
     sA[threadIdx.x].ymax = max(sA[threadIdx.x].ymax ,sA[threadIdx.x + stride].ymax);
     sA[threadIdx.x].zmax = max(sA[threadIdx.x].zmax ,sA[threadIdx.x + stride].zmax);
   }
   __syncthreads();
 }
  PMEDouble xcenter = (PMEFloat)-0.5 * (sA[0].xmin + sA[0].xmax);
  PMEDouble ycenter = (PMEFloat)-0.5 * (sA[0].ymin + sA[0].ymax);
  PMEDouble zcenter = (PMEFloat)-0.5 * (sA[0].zmin + sA[0].zmax);

  // Perform individual sums
  pos = blockIdx.x * blockDim.x + threadIdx.x;
  while (pos < cSim.atoms) {
    PMEDouble x = cSim.pAtomX[pos];
    PMEDouble y = cSim.pAtomY[pos];
    PMEDouble z = cSim.pAtomZ[pos];
    x += xcenter;
    y += ycenter;
    z += zcenter;
    PMEFloat2 xy = {(PMEFloat)x, (PMEFloat)y};
    cSim.pAtomX[pos] = x;
    cSim.pAtomY[pos] = y;
    cSim.pAtomZ[pos] = z;
    cSim.pAtomXYSP[pos] = xy;
    cSim.pAtomZSP[pos]  = z;
    pos += blockDim.x * gridDim.x;
  }

  // Fix restraints: update the original arrays for atom position constraints
  pos = blockIdx.x*blockDim.x + threadIdx.x;
  while (pos < cSim.constraints) {
    PMEDouble2 constraint1 = cSim.pConstraint1[pos];
    PMEDouble2 constraint2 = cSim.pConstraint2[pos];
    constraint1.y += xcenter;
    constraint2.x += ycenter;
    constraint2.y += zcenter;
    cSim.pConstraint1[pos] = constraint1;
    cSim.pConstraint2[pos] = constraint2;
    pos += blockDim.x * gridDim.x;
  }

  // Fix restraints: update the bond work units
  pos = blockIdx.x*blockDim.x + threadIdx.x;
 unsigned int tgx = threadIdx.x & GRID_BITS_MASK;
  while (pos < cSim.BwuCnstCount) {
    int warpPos = 2 * (pos / GRID) * GRID;
    PMEDouble2 constraint1 = cSim.pBwuCnst[warpPos + tgx];
    PMEDouble2 constraint2 = cSim.pBwuCnst[warpPos + GRID + tgx];
    constraint1.y += xcenter;
    constraint2.x += ycenter;
    constraint2.y += zcenter;
    cSim.pBwuCnst[warpPos + tgx       ] = constraint1;
    cSim.pBwuCnst[warpPos + GRID + tgx] = constraint2;
    pos += blockDim.x * gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kRecenter_Molecule: this function removes the displacement of the system center, as measured
//                     by the midpoint of extrema in three cartesian directions.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kRecenter_Molecule(gpuContext gpu)
{
  kRecenter_Molecule1_kernel<<<gpu->blocks, gpu->generalThreadsPerBlock>>>();
  LAUNCHERROR("kRecenter_Molecule");
  kRecenter_Molecule2_kernel<<<gpu->blocks, gpu->generalThreadsPerBlock>>>();
  LAUNCHERROR("kRecenter_Molecule");
}

//---------------------------------------------------------------------------------------------
// kCalculateMaxGradient_kernel: find the maximum gradient out of all force accumulators.
//---------------------------------------------------------------------------------------------
__device__ unsigned long long int dMaxLength;
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kCalculateMaxGradient_kernel()
{
#if (__CUDA_ARCH__ >= 350) || defined(AMBER_PLATFORM_AMD)
  if (threadIdx.x == 0) {
    dMaxLength = 0;
  }
  __threadfence();
  __syncthreads();
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  double maxlength = 0.0;
  while (pos < cSim.atoms) {
    double fx = cSim.pForceAccumulator[pos];
    double fy = cSim.pForceAccumulator[pos + cSim.stride];
    double fz = cSim.pForceAccumulator[pos + cSim.stride2];
    double length = fx * fx + fy * fy + fz * fz;
    if (length > maxlength) {
      maxlength = length;
    }
    pos += increment;
  }
  maxlength = sqrt(maxlength);
  unsigned long long int ml = llitoulli(llround(maxlength));
  atomicMax(&dMaxLength, ml);
#endif
}

//---------------------------------------------------------------------------------------------
// Force update kernels.  The differentiators are the presence of a Langevin thermostat, a
// neighborlist, and constant pressure conditions. (Constant pressure conditions necessarily
// imply the presence of a neighbor list.) These kernels take code from kU.h.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kUpdate_kernel(PMEDouble dt)
#include "kU.h"

//---------------------------------------------------------------------------------------------
#define UPDATE_LANGEVIN
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kLangevinUpdate_kernel(PMEDouble dt, PMEDouble temp0, PMEDouble gamma_ln, int rpos)
#include "kU.h"
#undef UPDATE_LANGEVIN

//---------------------------------------------------------------------------------------------
#define UPDATE_NEIGHBORLIST
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kNLUpdate_kernel(PMEDouble dt)
#include "kU.h"

//---------------------------------------------------------------------------------------------
#define UPDATE_LANGEVIN
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kNLLangevinUpdate_kernel(PMEDouble dt, PMEDouble temp0, PMEDouble gamma_ln, int rpos)
#include "kU.h"
#undef UPDATE_LANGEVIN

//---------------------------------------------------------------------------------------------
#define UPDATE_NTP
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kNLNTPUpdate_kernel(PMEDouble dt)
#include "kU.h"

//---------------------------------------------------------------------------------------------
#define UPDATE_LANGEVIN
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kNLNTPLangevinUpdate_kernel(PMEDouble dt, PMEDouble temp0, PMEDouble gamma_ln, int rpos)
#include "kU.h"
#undef UPDATE_LANGEVIN
#undef UPDATE_NTP
#undef UPDATE_NEIGHBORLIST

__global__ void
kUpdateForces_PHMDNL_kernel()
{
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  int restype = cSim.psp_grp[pos];
  int atmX, lowerbound, upperbound, patmX;
  if(restype == 0 || restype == 1 || restype == 3) {
    lowerbound = cSim.presbounds[pos].x - 1;
    upperbound = cSim.presbounds[pos].y - 1;
  }
  else {
    lowerbound = cSim.presbounds[pos-1].x - 1;
    upperbound = cSim.presbounds[pos-1].y - 1;
  }
  for(atmX = lowerbound; atmX <= upperbound; atmX++) {
    patmX = cSim.pImageAtomLookup[atmX];
    if(restype == 0 || restype == 1 || restype == 3) {
      cSim.pdph_theta[pos] += (double)cSim.pdph_accumulator[patmX] * (double)ONEOVERFORCESCALE;
    }
    else {
      cSim.pdph_theta[pos] += (double)cSim.pdph_plus_accumulator[patmX] * (double)ONEOVERFORCESCALE;
    }
  }
}

__global__ void
kUpdateForces_PHMD_kernel()
{
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  int restype = cSim.psp_grp[pos];
  int atmX, lowerbound, upperbound, patmX;
  if(restype == 0 || restype == 1 || restype == 3) {
    lowerbound = cSim.presbounds[pos].x - 1;
    upperbound = cSim.presbounds[pos].y - 1;
  }
  else {
    lowerbound = cSim.presbounds[pos-1].x - 1;
    upperbound = cSim.presbounds[pos-1].y - 1;
  }
  for(atmX = lowerbound; atmX <= upperbound; atmX++) {
    if(restype == 0 || restype == 1 || restype == 3) {
      cSim.pdph_theta[pos] += (double)cSim.pdph_accumulator[atmX] * (double)ONEOVERFORCESCALE;
    }
    else {
      cSim.pdph_theta[pos] += (double)cSim.pdph_plus_accumulator[atmX] * (double)ONEOVERFORCESCALE;
    }
  }
}
__global__ void
kLangevinUpdatePHMD_kernel(PMEDouble dtx, int rpos)
{
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  PMEFloat kbt = KB_VAL * cSim.temp_phmd;
  PMEFloat gam = 0.0488882 * cSim.phbeta * dtx;
  PMEFloat mass = cSim.qmass_phmd;
  PMEFloat rfd = sqrt(2.0 * mass * gam * kbt) / dtx;
  PMEFloat dph_theta = cSim.pdph_theta[pos];
  cSim.pdudls[pos] = dph_theta;
  if(!cSim.use_forces) {
    dph_theta = 0;
  }
  dph_theta *= sin(2.0 * cSim.pph_theta[pos]);
  if( pos < cSim.ntitr ) {
    dph_theta += rfd * cSim.pRandomPHMD[rpos*cSim.phmd_padded_ntitr+pos];
  }
  PMEFloat dph_theta_old = dph_theta;
  cSim.pvph_theta[pos] = cSim.pvphold[pos] -
                         0.5 * dtx * dph_theta / mass -
                         0.5 * dtx * cSim.pdphold[pos] / mass;
  PMEFloat pvph_theta = (1.0 - gam ) * cSim.pvph_theta[pos];
  cSim.pph_theta[pos] += pvph_theta * dtx -
                        dph_theta * dtx * dtx / ( 2.0 * mass );
  cSim.pdphold[pos] = dph_theta_old;
  cSim.pvphold[pos] = pvph_theta;
  cSim.pdph_theta[pos] = dph_theta;
  cSim.pvph_theta[pos] = pvph_theta;
}

//---------------------------------------------------------------------------------------------
// kChargeUpdatePHMD_kernel: kernel to update phmd charges.
//
//---------------------------------------------------------------------------------------------
__global__ void
kChargeUpdatePHMD_kernel()
{
  unsigned int rpos = blockIdx.x;
  unsigned int pos = threadIdx.x;
  int i;
  double x1, qprot, qunprot;
  double x2 = 1.0;
  if( cSim.pbaseres[rpos] ) {
    x1 = sin( cSim.pph_theta[rpos] );
    x1 *= x1;
    if( cSim.psp_grp[rpos] > 0 ) {
      x2 = sin( cSim.pph_theta[rpos+1] );
      x2 *= x2;
    }
    i = cSim.presbounds[rpos].x + pos;
    if( i < cSim.presbounds[rpos].y ) {
      qprot = x2 * cSim.pqstate1[i].x + ( 1.0 - x2 ) * cSim.pqstate1[i].y;
      qunprot = x2 * cSim.pqstate2[i].x + ( 1.0 - x2 ) * cSim.pqstate2[i].y;
      cSim.pcharge_phmd[i] = ( 1.0 - x1 ) * qprot + x1 * qunprot; 
      qprot = x2 * cSim.pqstate1_md[i].x + ( 1.0 - x2 ) * cSim.pqstate1_md[i].y;
      qunprot = x2 * cSim.pqstate2_md[i].x + ( 1.0 - x2 ) * cSim.pqstate2_md[i].y;
      double charge = ( 1.0 - x1 ) * qprot + x1 * qunprot;
      cSim.pAtomChargeSPLJID[i].x = charge;
      cSim.pAtomCharge[i] = charge;
      cSim.pAtomChargeSP[i] = charge;
    }
  }
}

//---------------------------------------------------------------------------------------------
// kChargeUpdatePHMDNL_kernel: kernel to update phmd charges.
//
//---------------------------------------------------------------------------------------------
__global__ void
kChargeUpdatePHMDNL_kernel()
{
  unsigned int pos = (blockIdx.x * blockDim.x + threadIdx.x) * cSim.phmd_atom_stride;
  PMEFloat x1, x2, qprot, qunprot, charge;
  for (int i = 0; i < cSim.phmd_atom_stride; i++ ){
    int totpos = pos + i;
    if (totpos < cSim.atoms) {
      int h = cSim.pgrplist[totpos] - 1;
      int k = cSim.plinkgrplist[totpos] - 1;
      if (h >= 0 || k >= 0) {
        int atmx = cSim.pImageAtomLookup[pos + i];
        if (h >= 0) {
          x1 = sin(cSim.pph_theta[h]);
          x2 = 1.0;
          if (cSim.psp_grp[h] > 0) {
            x2 = sin(cSim.pph_theta[h+1]);
            x2 *= x2;
          }
          x1 *= x1;
          qprot = x2 * cSim.pqstate1[pos + i].x + (1.0 - x2) * cSim.pqstate1[pos + i].y;
          qunprot = x2 * cSim.pqstate2[pos + i].x + (1.0 - x2) * cSim.pqstate2[pos + i].y;
          charge = (1.0 - x1) * qprot + x1 * qunprot;
        }
        else {
          x1 = sin(cSim.pph_theta[k]);
          x1 *= x1;
          charge = x1 * cSim.pqstate1[pos + i].x + (1.0 - x1) * cSim.pqstate2[pos + i].x; 
        }
        cSim.pcharge_phmd[pos + i] = charge;
        cSim.pImageCharge_phmd[atmx] = charge;
        cSim.pAtomChargeSPLJID[atmx].x = charge;
        cSim.pAtomCharge[atmx] = charge;
        cSim.pAtomChargeSP[atmx] = charge;
        cSim.pImageCharge[atmx] = charge;
      }
    }
  }
}

//---------------------------------------------------------------------------------------------
// kUpdatePHMD: launch kernel to update titration forces.
//
// Arguments:
//   dtx:       time step
//   gpu:       overarching type for storing parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kUpdatePHMD(PMEDouble dtx, gpuContext gpu)
{
  if (gpu->randomCounterPHMD >= gpu->sim.randomSteps) {
    if (gpu->bCPURandoms) {
      cpu_kRandomPHMD(gpu);
    }
    else {
      kRandomPHMD(gpu);
    }
    gpu->randomCounterPHMD = 0;
  }
  kUpdateForces_PHMD(gpu);
  cudaDeviceSynchronize();
  kLangevinUpdatePHMD_kernel<<<1,gpu->sim.phmd_padded_ntitr>>>(dtx, gpu->randomCounterPHMD);
  cudaDeviceSynchronize();
  if (gpu->bNeighborList) {
    kChargeUpdatePHMDNL_kernel<<<gpu->PMENonbondBlocks, gpu->PMENonbondForcesThreadsPerBlock>>>();
    gpu->pbAtomChargeSP->Download();
    gpu->sumq = 0;
    gpu->sumq2 = 0;
    for (int i = 0; i < gpu->sim.atoms; i++) {
      double charge = gpu->pbAtomChargeSP->_pSysData[i];
      gpu->sumq += charge;
      gpu->sumq2 += charge * charge;
    }
  }
  else {
    kChargeUpdatePHMD_kernel<<<gpu->sim.ntitr,gpu->sim.padded_numch>>>();
  }
  cudaDeviceSynchronize();
  gpu->randomCounterPHMD++;
}
//---------------------------------------------------------------------------------------------
// kUpdateForces_PHMD: launch kernel to update titration forces.
//
// Arguments:
//   gpu:       overarching type for storing parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kUpdateForces_PHMD(gpuContext gpu)
{
  if (gpu->bNeighborList) {
    kUpdateForces_PHMDNL_kernel<<<1,gpu->sim.phmd_padded_ntitr>>>();
  }
  else {
    kUpdateForces_PHMD_kernel<<<1,gpu->sim.phmd_padded_ntitr>>>();
  }
  cudaDeviceSynchronize();
}
//---------------------------------------------------------------------------------------------
// kUpdate: launch the appropriate kernel to update forces.
//
// Arguments:
//   gpu:       overarching type for storing parameters, coordinates, and the energy function
//   dt:        the time step
//   temp0:     target temperature (for constant temperature runs)
//   gamma_ln:  Langevin collision frequency
//---------------------------------------------------------------------------------------------
void kUpdate(gpuContext gpu, PMEDouble dt, PMEDouble temp0, PMEDouble gamma_ln)
{
  // Local variables for the tile size
  int upBlocks = gpu->updateBlocks;
  int upThreads = gpu->updateThreadsPerBlock;

  // Choose Langevin update if necessary
  if (gpu->ntt == 3) {

    // Kernel names... too long... breaking my line format...
    int spNA = gpu->sim.paddedNumberOfAtoms;

    // Update random numbers if necessary
    if (gpu->randomCounter >= gpu->sim.randomSteps) {
      if (gpu->bCPURandoms) {
        cpu_kRandom(gpu);
      }
      else {
        kRandom(gpu);
      }
      gpu->randomCounter = 0;
    }
    if (gpu->bNeighborList) {
      if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
        kNLNTPLangevinUpdate_kernel<<<upBlocks, upThreads>>>(dt, temp0, gamma_ln,
                                                             gpu->randomCounter * spNA);
      }
      else {
        kNLLangevinUpdate_kernel<<<upBlocks, upThreads>>>(dt, temp0, gamma_ln,
                                                          gpu->randomCounter * spNA);
      }
    }
    else {
      kLangevinUpdate_kernel<<<upBlocks, upThreads>>>(dt, temp0, gamma_ln,
                                                      gpu->randomCounter * spNA);
    }
    gpu->randomCounter++;
  }
  else {
    if (gpu->bNeighborList) {
      if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
        kNLNTPUpdate_kernel<<<upBlocks, upThreads>>>(dt);
      }
      else {
        kNLUpdate_kernel<<<upBlocks, upThreads>>>(dt);
      }
    }
    else {
      kUpdate_kernel<<<upBlocks, upThreads>>>(dt);
    }
  }
  LAUNCHERROR("kUpdate");
}


//---------------------------------------------------------------------------------------------
// Kernels for updating forces with SGLD.  
//---------------------------------------------------------------------------------------------
#define UPDATE_SGLD

__global__ void
__launch_bounds__(UPDATE_THREADS_PER_BLOCK, 1)
kSGMDUpdate_kernel(PMEDouble dt, PMEDouble gamma_sg, PMEDouble dcomx1, PMEDouble dcomy1, PMEDouble dcomz1, PMEDouble dcomx2, PMEDouble dcomy2, PMEDouble dcomz2)
#include "kU.h"

//---------------------------------------------------------------------------------------------
#define UPDATE_LANGEVIN
__global__ void
__launch_bounds__(UPDATE_THREADS_PER_BLOCK, 1)
kSGLDUpdate_kernel(PMEDouble dt, PMEDouble temp0, PMEDouble gamma_ln, PMEDouble gamma_sg, PMEDouble dcomx1, PMEDouble dcomy1, PMEDouble dcomz1, PMEDouble dcomx2, PMEDouble dcomy2, PMEDouble dcomz2, int rpos)
#include "kU.h"
#undef UPDATE_LANGEVIN

//---------------------------------------------------------------------------------------------
#define UPDATE_NEIGHBORLIST
__global__ void
__launch_bounds__(UPDATE_THREADS_PER_BLOCK, 1)
kNLSGMDUpdate_kernel(PMEDouble dt, PMEDouble gamma_sg, PMEDouble dcomx1, PMEDouble dcomy1, PMEDouble dcomz1, PMEDouble dcomx2, PMEDouble dcomy2, PMEDouble dcomz2)
#include "kU.h"

//---------------------------------------------------------------------------------------------
#define UPDATE_LANGEVIN
__global__ void
__launch_bounds__(UPDATE_THREADS_PER_BLOCK, 1)
kNLSGLDUpdate_kernel(PMEDouble dt, PMEDouble temp0, PMEDouble gamma_ln, PMEDouble gamma_sg, PMEDouble dcomx1, PMEDouble dcomy1, PMEDouble dcomz1, PMEDouble dcomx2, PMEDouble dcomy2, PMEDouble dcomz2, int rpos)
#include "kU.h"
#undef UPDATE_LANGEVIN

//---------------------------------------------------------------------------------------------
#define UPDATE_NTP
__global__ void
__launch_bounds__(UPDATE_THREADS_PER_BLOCK, 1)
kNLNTPSGMDUpdate_kernel(PMEDouble dt, PMEDouble gamma_sg, PMEDouble dcomx1, PMEDouble dcomy1, PMEDouble dcomz1, PMEDouble dcomx2, PMEDouble dcomy2, PMEDouble dcomz2)
#include "kU.h"

//---------------------------------------------------------------------------------------------
#define UPDATE_LANGEVIN
__global__ void
__launch_bounds__(UPDATE_THREADS_PER_BLOCK, 1)
kNLNTPSGLDUpdate_kernel(PMEDouble dt, PMEDouble temp0, PMEDouble gamma_ln, PMEDouble gamma_sg, PMEDouble dcomx1, PMEDouble dcomy1, PMEDouble dcomz1, PMEDouble dcomx2, PMEDouble dcomy2, PMEDouble dcomz2, int rpos)
#include "kU.h"
#undef UPDATE_LANGEVIN
#undef UPDATE_NTP
#undef UPDATE_NEIGHBORLIST

#undef UPDATE_SGLD


//---------------------------------------------------------------------------------------------
// kUpdateSGLD: launch the appropriate SGLD kernel to update forces.
//
// Arguments:
//   gpu:       overarching type for storing parameters, coordinates, and the energy function
//   dt:        the time step
//   temp0:     target temperature (for constant temperature runs)
//   gamma_ln:  Langevin collision frequency
//---------------------------------------------------------------------------------------------
void kUpdateSGLD(gpuContext gpu, PMEDouble dt, PMEDouble temp0, PMEDouble gamma_ln, PMEDouble gamma_sg, PMEDouble com0sg[], PMEDouble com1sg[], PMEDouble com2sg[])
{
  // Local variables for the tile size
  int upBlocks = gpu->updateBlocks;
  int upThreads = gpu->updateThreadsPerBlock;
        double dcomx1                               = com0sg[0]-com1sg[0];
        double dcomy1                               = com0sg[1]-com1sg[1];
        double dcomz1                               = com0sg[2]-com1sg[2];
        double dcomx2                               = com0sg[0]-com2sg[0];
        double dcomy2                               = com0sg[1]-com2sg[1];
        double dcomz2                               = com0sg[2]-com2sg[2];

  // Choose Langevin update if necessary
  if (gpu->ntt == 3) {

    // Kernel names... too long... breaking my line format...
    int spNA = gpu->sim.paddedNumberOfAtoms;

    // Update random numbers if necessary
    if (gpu->randomCounter >= gpu->sim.randomSteps) {
      if (gpu->bCPURandoms) {
        cpu_kRandom(gpu);
      }
      else {
        kRandom(gpu);
      }
      gpu->randomCounter = 0;
    }
    if (gpu->bNeighborList) {
      if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
        kNLNTPSGLDUpdate_kernel<<<upBlocks, upThreads>>>(dt, temp0, gamma_ln, gamma_sg, dcomx1,dcomy1,dcomz1,dcomx2,dcomy2,dcomz2,
                             gpu->randomCounter * spNA);
      }
      else {
        kNLSGLDUpdate_kernel<<<upBlocks, upThreads>>>(dt, temp0, gamma_ln, gamma_sg, dcomx1,dcomy1,dcomz1,dcomx2,dcomy2,dcomz2,
                                     gpu->randomCounter * spNA);
      }
    }
    else {
      kSGLDUpdate_kernel<<<upBlocks, upThreads>>>(dt, temp0, gamma_ln, gamma_sg, dcomx1,dcomy1,dcomz1,dcomx2,dcomy2,dcomz2,
                                 gpu->randomCounter * spNA);
    }
    gpu->randomCounter++;
  }
  else {
    if (gpu->bNeighborList) {
      if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
        kNLNTPSGMDUpdate_kernel<<<upBlocks, upThreads>>>(dt, gamma_sg, dcomx1,dcomy1,dcomz1,dcomx2,dcomy2,dcomz2);
      }
      else {
        kNLSGMDUpdate_kernel<<<upBlocks, upThreads>>>(dt, gamma_sg, dcomx1,dcomy1,dcomz1,dcomx2,dcomy2,dcomz2);
      }
    }
    else {
      kSGMDUpdate_kernel<<<upBlocks, upThreads>>>(dt, gamma_sg, dcomx1,dcomy1,dcomz1,dcomx2,dcomy2,dcomz2);
    }
  }
  LAUNCHERROR("kUpdateSGLD");

}

//---------------------------------------------------------------------------------------------
// Kernels for updating forces with middle-schemed molecular dynamics.  The differentiators are,
// again, the use of a neighbor list, and constant pressure conditions.
// These kernels again take code from kMiddle.h.
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------

#define UPDATE_MIDDLE_SCHEME_1
__global__ void
__launch_bounds__(UPDATE_THREADS_PER_BLOCK, 1)
kUpdateMiddle1(PMEDouble dt)
#include "kMiddle.h"
#undef UPDATE_MIDDLE_SCHEME_1

#define UPDATE_MIDDLE_SCHEME_2
__global__ void
__launch_bounds__(UPDATE_THREADS_PER_BLOCK, 1)
kUpdateMiddle2(PMEDouble dt, PMEDouble temp0, PMEDouble gamma_ln, int rpos)
#include "kMiddle.h"
#undef UPDATE_MIDDLE_SCHEME_2

#define UPDATE_NEIGHBORLIST
//---------------------------------------------------------------------------------------------

#define UPDATE_MIDDLE_SCHEME_1
__global__ void
__launch_bounds__(UPDATE_THREADS_PER_BLOCK, 1)
kNLUpdateMiddle1(PMEDouble dt)
#include "kMiddle.h"
#undef UPDATE_MIDDLE_SCHEME_1

#define UPDATE_MIDDLE_SCHEME_2
__global__ void
__launch_bounds__(UPDATE_THREADS_PER_BLOCK, 1)
kNLUpdateMiddle2(PMEDouble dt, PMEDouble temp0, PMEDouble gamma_ln, int rpos)
#include "kMiddle.h"
#undef UPDATE_MIDDLE_SCHEME_2


//---------------------------------------------------------------------------------------------
#define UPDATE_NTP

#define UPDATE_MIDDLE_SCHEME_1
__global__ void
__launch_bounds__(UPDATE_THREADS_PER_BLOCK, 1)
kNLNTPUpdateMiddle1(PMEDouble dt)
#include "kMiddle.h"
#undef UPDATE_MIDDLE_SCHEME_1

//---------------------------------------------------------------------------------------------
#define UPDATE_MIDDLE_SCHEME_2
__global__ void
__launch_bounds__(UPDATE_THREADS_PER_BLOCK, 1)
kNLNTPUpdateMiddle2(PMEDouble dt, PMEDouble temp0, PMEDouble gamma_ln, int rpos)
#include "kMiddle.h"
#undef UPDATE_MIDDLE_SCHEME_2
#undef UPDATE_NTP
#undef UPDATE_NEIGHBORLIST

//---------------------------------------------------------------------------------------------
// kUpdateMiddle: launch the appropriate kernel to update forces.
//
// Arguments:
//   gpu:       overarching type for storing parameters, coordinates, and the energy function
//   dt:        the time step
//   temp0:     target temperature (for constant temperature runs)
//   gamma_ln:  Langevin collision frequency
//---------------------------------------------------------------------------------------------
void kUpdateMiddle(gpuContext gpu, PMEDouble dt, PMEDouble temp0, PMEDouble gamma_ln)
{
  // Local variables for the tile size
  int upBlocks = gpu->updateBlocks;
  int upThreads = gpu->updateThreadsPerBlock;

  if (gpu->ischeme == 1) { // Middle scheme procedure
    // Kernel names... too long... breaking my line format...
    int spNA = gpu->sim.paddedNumberOfAtoms;

    // Update random numbers if necessary
    if (gpu->randomCounter >= gpu->sim.randomSteps) {
      if (gpu->bCPURandoms) {
        cpu_kRandom(gpu);
      }
      else {
        kRandom(gpu);
      }
      gpu->randomCounter = 0;
    }
    if (gpu->bNeighborList) {
      if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          kNLNTPUpdateMiddle1<<<upBlocks, upThreads>>>(dt);
      }
      else {
          kNLUpdateMiddle1<<<upBlocks, upThreads>>>(dt);
      }
    }
    else {
        kUpdateMiddle1<<<upBlocks, upThreads>>>(dt);
    }
    if (gpu->ntc != 1) {
        kRattle(gpu, dt);
    }
    if (gpu->bNeighborList) {
      if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          kNLNTPUpdateMiddle2<<<upBlocks, upThreads>>>(dt, temp0, gamma_ln, gpu->randomCounter * spNA);
      }
      else {
          kNLUpdateMiddle2<<<upBlocks, upThreads>>>(dt, temp0, gamma_ln, gpu->randomCounter * spNA);
      }
    }
    else {
        kUpdateMiddle2<<<upBlocks, upThreads>>>(dt, temp0, gamma_ln, gpu->randomCounter * spNA);
    }
    gpu->randomCounter++;
  }
  LAUNCHERROR("kUpdateMiddle");

}


//---------------------------------------------------------------------------------------------
// Kernels for updating forces with relaxed molecular dynamics.  The differentiators are,
// again, the use of a Langevin thermostat, a neighbor list, and constant pressure conditions.
// These kernels again take code from kU.h.
//---------------------------------------------------------------------------------------------
#define UPDATE_RELAXMD
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kRelaxMDUpdate_kernel(PMEDouble dt)
#include "kU.h"

//---------------------------------------------------------------------------------------------
#define UPDATE_LANGEVIN
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kRelaxMDLangevinUpdate_kernel(PMEDouble dt, PMEDouble temp0, PMEDouble gamma_ln, int rpos)
#include "kU.h"
#undef UPDATE_LANGEVIN

//---------------------------------------------------------------------------------------------
#define UPDATE_NEIGHBORLIST
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kNLRelaxMDUpdate_kernel(PMEDouble dt)
#include "kU.h"

//---------------------------------------------------------------------------------------------
#define UPDATE_LANGEVIN
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kNLRelaxMDLangevinUpdate_kernel(PMEDouble dt, PMEDouble temp0, PMEDouble gamma_ln, int rpos)
#include "kU.h"
#undef UPDATE_LANGEVIN

//---------------------------------------------------------------------------------------------
#define UPDATE_NTP
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kNLNTPRelaxMDUpdate_kernel(PMEDouble dt)
#include "kU.h"

//---------------------------------------------------------------------------------------------
#define UPDATE_LANGEVIN
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kNLNTPRelaxMDLangevinUpdate_kernel(PMEDouble dt, PMEDouble temp0, PMEDouble gamma_ln, int rpos)
#include "kU.h"
#undef UPDATE_LANGEVIN
#undef UPDATE_NTP
#undef UPDATE_NEIGHBORLIST
#undef UPDATE_RELAXMD

//---------------------------------------------------------------------------------------------
// kRelaxMDUpdate: launch the appropriate kernel to update forces in the context of relaxed MD.
//
// Arguments:
//   gpu:       overarching type for storing parameters, coordinates, and the energy function
//   dt:        the time step
//   temp0:     target temperature (for constant temperature runs)
//   gamma_ln:  Langevin collision frequency
//---------------------------------------------------------------------------------------------
void kRelaxMDUpdate(gpuContext gpu, PMEDouble dt, PMEDouble temp0, PMEDouble gamma_ln)
{
  // Local variables for the tile size
  int upBlocks = gpu->updateBlocks;
  int upThreads = gpu->updateThreadsPerBlock;

  // Choose Langevin update if necessary
  if (gpu->ntt == 3) {

    // Kernel names... too long... breaking my line format...
    int spNA = gpu->sim.paddedNumberOfAtoms;

    // Update random numbers if necessary
    if (gpu->randomCounter >= gpu->sim.randomSteps) {
      if (gpu->bCPURandoms) {
        cpu_kRandom(gpu);
      }
      else {
        kRandom(gpu);
      }
      gpu->randomCounter = 0;
    }
    if (gpu->bNeighborList) {
      if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
        kNLNTPRelaxMDLangevinUpdate_kernel<<<upBlocks, upThreads>>>(dt, temp0, gamma_ln,
                                                                    gpu->randomCounter * spNA);
      }
      else {
        kNLRelaxMDLangevinUpdate_kernel<<<upBlocks, upThreads>>>(dt, temp0, gamma_ln,
                                                                 gpu->randomCounter * spNA);
      }
    }
    else {
      kRelaxMDLangevinUpdate_kernel<<<upBlocks, upThreads>>>(dt, temp0, gamma_ln,
                                                             gpu->randomCounter * spNA);
    }
    gpu->randomCounter++;
  }
  else {
    if (gpu->bNeighborList) {
      if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
        kNLNTPRelaxMDUpdate_kernel<<<upBlocks, upThreads>>>(dt);
      }
      else {
        kNLRelaxMDUpdate_kernel<<<upBlocks, upThreads>>>(dt);
      }
    }
    else {
      kRelaxMDUpdate_kernel<<<upBlocks, upThreads>>>(dt);
    }
  }
  LAUNCHERROR("kRelaxMDUpdate");
}

//---------------------------------------------------------------------------------------------
// kRefreshCharges_kernel: kernel that is ultimately called by a chain of functions set in
//                         motion for constant pH calculations.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kRefreshCharges_kernel()
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  while (pos < cSim.atoms) {
    PMEDouble charge              = cSim.pChargeRefreshBuffer[pos];
    cSim.pAtomCharge[pos]         = charge;
    cSim.pAtomChargeSP[pos]       = charge;
    cSim.pAtomChargeSPLJID[pos].x = charge;
    pos += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kNLRefreshCharges_kernel: kernel that is ultimately called by a chain of functions set in
//                           motion for constant pH calculations, in the context of a neighbor
//                           list.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kNLRefreshCharges_kernel()
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  while (pos < cSim.atoms) {
    PMEDouble charge                = cSim.pChargeRefreshBuffer[pos];
    unsigned int index              = cSim.pImageAtomLookup[pos];
    cSim.pImageCharge[index]        = charge;
    cSim.pAtomChargeSP[index]       = charge;
    cSim.pAtomChargeSPLJID[index].x = charge;
    pos += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kRefreshCharges: host function to launch one of the two kernels above, called by
//                  gpu_refresh_charges_ in gpu.cpp as part of constant pH calculations.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kRefreshCharges(gpuContext gpu)
{
  if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
    kNLRefreshCharges_kernel<<<gpu->blocks, gpu->generalThreadsPerBlock>>>();
  }
  else {
    kRefreshCharges_kernel<<<gpu->blocks, gpu->generalThreadsPerBlock>>>();
  }
  LAUNCHERROR("kRefreshCharges");
}

//---------------------------------------------------------------------------------------------
// kRefreshChargesGBCpH: refresh charges in the context of implicit solvent constant pH
//                       calculations.  It appears that this will never involve a neighbor
//                       list.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kRefreshChargesGBCpH(gpuContext gpu)
{
  kRefreshCharges_kernel<<<gpu->blocks, gpu->generalThreadsPerBlock>>>();
  LAUNCHERROR("kRefreshCharges");
}

//---------------------------------------------------------------------------------------------
// Kernels to reset velocities.  Critical details are the presence of a neighbor list or
// constant pressure conditions.  Draws on code in kRV.h.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kResetVelocities_kernel(PMEDouble temp, PMEDouble half_dtx, int rpos)
#include "kRV.h"

//---------------------------------------------------------------------------------------------
#define RV_NEIGHBORLIST
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kNLResetVelocities_kernel(PMEDouble temp, PMEDouble half_dtx, int rpos)
#include "kRV.h"

//---------------------------------------------------------------------------------------------
#define RV_NTP
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kNLNTPResetVelocities_kernel(PMEDouble temp, PMEDouble half_dtx, int rpos)
#include "kRV.h"

//---------------------------------------------------------------------------------------------
// Zero out the velocities.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(CLEARFORCES_THREADS_PER_BLOCK, 1)
kClearVelocities_kernel()
{
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  while(pos < cSim.atoms){
    cSim.pVelX[pos] = (double)0.0;
    cSim.pVelY[pos] = (double)0.0;
    cSim.pVelZ[pos] = (double)0.0;
    pos += blockDim.x * gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(CLEARFORCES_THREADS_PER_BLOCK, 1)
kNLClearVelocities_kernel()
{
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  if (pos < cSim.atoms) {
    cSim.pImageVelX[pos] = (double)0.0;
    cSim.pImageVelY[pos] = (double)0.0;
    cSim.pImageVelZ[pos] = (double)0.0;
    pos += blockDim.x * gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kClearVelocities : zero out the velocities
//---------------------------------------------------------------------------------------------
void kClearVelocities(gpuContext gpu)
{
  if (gpu->bNeighborList) {
    kNLClearVelocities_kernel<<<gpu->updateBlocks, gpu->updateThreadsPerBlock>>>();
  }
  else {
    kClearVelocities_kernel<<<gpu->updateBlocks, gpu->updateThreadsPerBlock>>>();
  }
  LAUNCHERROR("kClearVelocities");
}

//---------------------------------------------------------------------------------------------
// kResetVelocities: launch the appropriate kernel for resetting or updating velocities.
//
// Arguments:
//   gpu:       overarching type for storing parameters, coordinates, and the energy function
//   temp:      target temperature
//   half_dtx:  half the time step scaled by sqrt(418.4)
//---------------------------------------------------------------------------------------------
void kResetVelocities(gpuContext gpu, double temp, double half_dtx)
{
  // Update random numbers if necessary
  if (gpu->randomCounter >= gpu->sim.randomSteps) {
    if (gpu->bCPURandoms) {
      cpu_kRandom(gpu);
    }
    else {
      kRandom(gpu);
    }
    gpu->randomCounter = 0;
  }

  // Local variables for the tile size
  int upBlocks = gpu->updateBlocks;
  int upThreads = gpu->updateThreadsPerBlock;
  int spNA = gpu->sim.paddedNumberOfAtoms;

  // Choose the right kernel.
  if (gpu->bNeighborList) {
    if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
      kNLNTPResetVelocities_kernel<<<upBlocks, upThreads>>>(temp, half_dtx,
                                                            gpu->randomCounter * spNA);
    }
    else {
      kNLResetVelocities_kernel<<<upBlocks, upThreads>>>(temp, half_dtx,
                                                         gpu->randomCounter * spNA);
    }
  }
  else {
    kResetVelocities_kernel<<<upBlocks, upThreads>>>(temp, half_dtx,
                                                     gpu->randomCounter * spNA);
  }
  gpu->randomCounter++;
  LAUNCHERROR("kResetVelocities");
}

//---------------------------------------------------------------------------------------------
// Kernels for recalculating velocities (to satisfy the imposition of geometry constraints,
// which necessarily change the velocities as positions change).
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kRecalculateVelocities_kernel(PMEDouble dtx_inv)
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;

  if (pos < cSim.atoms) {
    PMEDouble oldAtomX = cSim.pOldAtomX[pos];
    PMEDouble atomX    = cSim.pAtomX[pos];
    PMEDouble oldAtomY = cSim.pOldAtomY[pos];
    PMEDouble atomY    = cSim.pAtomY[pos];
    PMEDouble oldAtomZ = cSim.pOldAtomZ[pos];
    PMEDouble atomZ    = cSim.pAtomZ[pos];
    PMEDouble velX     = (atomX - oldAtomX) * dtx_inv;
    PMEDouble velY     = (atomY - oldAtomY) * dtx_inv;
    PMEDouble velZ     = (atomZ - oldAtomZ) * dtx_inv;
    cSim.pVelX[pos]    = velX;
    cSim.pVelY[pos]    = velY;
    cSim.pVelZ[pos]    = velZ;
  }
}

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kPMERecalculateVelocities_kernel(PMEDouble dtx_inv)
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;

  if (pos < cSim.atoms) {
    PMEDouble oldAtomX   = cSim.pOldAtomX[pos];
    PMEDouble atomX      = cSim.pImageX[pos];
    PMEDouble oldAtomY   = cSim.pOldAtomY[pos];
    PMEDouble atomY      = cSim.pImageY[pos];
    PMEDouble oldAtomZ   = cSim.pOldAtomZ[pos];
    PMEDouble atomZ      = cSim.pImageZ[pos];
    PMEDouble velX       = (atomX - oldAtomX) * dtx_inv;
    PMEDouble velY       = (atomY - oldAtomY) * dtx_inv;
    PMEDouble velZ       = (atomZ - oldAtomZ) * dtx_inv;
    cSim.pImageVelX[pos] = velX;
    cSim.pImageVelY[pos] = velY;
    cSim.pImageVelZ[pos] = velZ;
  }
}

//---------------------------------------------------------------------------------------------
// kRecalculateVelocities: host function to launch the appropraite kernel for velocity
//                         rescaling due to geometry constraints.
//
// Arguments:
//   gpu:      overarching type for storing parameters, coordinates, and the energy function
//   dtx_inv:  inverse of [time step scaled by sqrt(418.4)]
//---------------------------------------------------------------------------------------------
void kRecalculateVelocities(gpuContext gpu, PMEDouble dtx_inv)
{
  if (gpu->bNeighborList) {
    kPMERecalculateVelocities_kernel<<<gpu->updateBlocks,
                                       gpu->updateThreadsPerBlock>>>(dtx_inv);
  }
  else {
    kRecalculateVelocities_kernel<<<gpu->updateBlocks, gpu->updateThreadsPerBlock>>>(dtx_inv);
  }
  LAUNCHERROR("kRecalculateVelocities");
}

//---------------------------------------------------------------------------------------------
// Kernels for recalculating velocities (to satisfy the imposition of geometry constraints,
// which necessarily change the velocities as positions change).
//---------------------------------------------------------------------------------------------
__global__ void
__launch_bounds__(UPDATE_THREADS_PER_BLOCK, 1)
kMiddleRecalculateVelocities_kernel(PMEDouble dtx_inv)
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;

  if (pos < cSim.atoms) {
    PMEDouble oldAtomX = cSim.pShakeOldAtomX[pos];
    PMEDouble atomX    = cSim.pAtomX[pos];
    PMEDouble oldAtomY = cSim.pShakeOldAtomY[pos];
    PMEDouble atomY    = cSim.pAtomY[pos];
    PMEDouble oldAtomZ = cSim.pShakeOldAtomZ[pos];
    PMEDouble atomZ    = cSim.pAtomZ[pos];
    PMEDouble velX     = (atomX - oldAtomX) * dtx_inv;
    PMEDouble velY     = (atomY - oldAtomY) * dtx_inv;
    PMEDouble velZ     = (atomZ - oldAtomZ) * dtx_inv;
    cSim.pVelX[pos]    += velX;
    cSim.pVelY[pos]    += velY;
    cSim.pVelZ[pos]    += velZ;
  }
}

//---------------------------------------------------------------------------------------------
__global__ void
__launch_bounds__(UPDATE_THREADS_PER_BLOCK, 1)
kMiddlePMERecalculateVelocities_kernel(PMEDouble dtx_inv)
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;

  if (pos < cSim.atoms) {
    PMEDouble oldAtomX   = cSim.pShakeOldAtomX[pos];
    PMEDouble atomX      = cSim.pImageX[pos];
    PMEDouble oldAtomY   = cSim.pShakeOldAtomY[pos];
    PMEDouble atomY      = cSim.pImageY[pos];
    PMEDouble oldAtomZ   = cSim.pShakeOldAtomZ[pos];
    PMEDouble atomZ      = cSim.pImageZ[pos];
    PMEDouble velX       = (atomX - oldAtomX) * dtx_inv;
    PMEDouble velY       = (atomY - oldAtomY) * dtx_inv;
    PMEDouble velZ       = (atomZ - oldAtomZ) * dtx_inv;
    cSim.pImageVelX[pos] += velX;
    cSim.pImageVelY[pos] += velY;
    cSim.pImageVelZ[pos] += velZ;
  }
}

//---------------------------------------------------------------------------------------------
// kRecalculateVelocities: host function to launch the appropraite kernel for velocity
//                         rescaling due to geometry constraints.
//
// Arguments:
//   gpu:      overarching type for storing parameters, coordinates, and the energy function
//   dtx_inv:  inverse of [time step scaled by sqrt(418.4)]
//---------------------------------------------------------------------------------------------
void kMiddleRecalculateVelocities(gpuContext gpu, PMEDouble dtx_inv)
{
  if (gpu->bNeighborList) {
    kMiddlePMERecalculateVelocities_kernel<<<gpu->updateBlocks,
                                       gpu->updateThreadsPerBlock>>>(dtx_inv);
  }
  else {
    kMiddleRecalculateVelocities_kernel<<<gpu->updateBlocks, gpu->updateThreadsPerBlock>>>(dtx_inv);
  }
  LAUNCHERROR("kRecalculateVelocities");
}

//---------------------------------------------------------------------------------------------
// Kernels for kinetic energy computation.  Once again there are distinct kernels depending on
// whether the context includes PME (which is equivalent to the presence of a neighbor list)
// and / or alchemical free energy calculations.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(THREADS_PER_BLOCK, 1)
kCalculateKineticEnergy_kernel(PMEFloat c_ave)
{
  extern __shared__ KineticEnergy sE[];

  PMEFloat eke   = (PMEFloat)0.0;
  PMEFloat ekph  = (PMEFloat)0.0;
  PMEFloat ekpbs = (PMEFloat)0.0;
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;

  // Sum up kinetic energies
  while (pos < cSim.atoms) {
    PMEFloat mass = cSim.pAtomMass[pos];
    PMEFloat vx   = cSim.pVelX[pos];
    PMEFloat lvx  = cSim.pLVelX[pos];
    PMEFloat vy   = cSim.pVelY[pos];
    PMEFloat lvy  = cSim.pLVelY[pos];
    PMEFloat vz   = cSim.pVelZ[pos];
    PMEFloat lvz  = cSim.pLVelZ[pos];
    PMEFloat svx  = vx + lvx;
    PMEFloat svy  = vy + lvy;
    PMEFloat svz  = vz + lvz;
    eke          += mass * (svx*svx + svy*svy + svz*svz);
    ekpbs        += mass * ( vx*lvx +  vy*lvy +  vz*lvz);
    ekph         += mass * ( vx*vx  +  vy*vy  +  vz*vz );
    pos          += blockDim.x * gridDim.x;
  }
  eke   *= (PMEFloat)0.125 * c_ave;
  ekph  *= (PMEFloat)0.5;
  ekpbs *= (PMEFloat)0.5;
  sE[threadIdx.x].KE.EKE   = eke;
  sE[threadIdx.x].KE.EKPH  = ekph;
  sE[threadIdx.x].KE.EKPBS = ekpbs;

  // Reduce per-thread kinetic energies
  __syncthreads();
  for (uint stride = blockDim.x >> 1; stride > 0; stride >>= 1) {
    if (threadIdx.x < stride) {
      sE[threadIdx.x].KE.EKE += sE[threadIdx.x + stride].KE.EKE;
      sE[threadIdx.x].KE.EKPH += sE[threadIdx.x + stride].KE.EKPH;
      sE[threadIdx.x].KE.EKPBS += sE[threadIdx.x + stride].KE.EKPBS;
    }
    __syncthreads();
  }

  // Save result
  if (threadIdx.x < 3) {
    cSim.pKineticEnergy[blockIdx.x].array[threadIdx.x] = sE[0].array[threadIdx.x];
  }
}

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(THREADS_PER_BLOCK, 1)
kPMECalculateKineticEnergy_kernel(PMEFloat c_ave)
{
  extern __shared__ KineticEnergy sE[];

  PMEFloat eke     = (PMEFloat)0.0;
  PMEFloat ekph    = (PMEFloat)0.0;
  PMEFloat ekpbs   = (PMEFloat)0.0;
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

  // Sum up kinetic energies
  while (pos < cSim.atoms) {
    PMEFloat mass = cSim.pImageMass[pos];
    PMEFloat vx   = cSim.pImageVelX[pos];
    PMEFloat lvx  = cSim.pImageLVelX[pos];
    PMEFloat vy   = cSim.pImageVelY[pos];
    PMEFloat lvy  = cSim.pImageLVelY[pos];
    PMEFloat vz   = cSim.pImageVelZ[pos];
    PMEFloat lvz  = cSim.pImageLVelZ[pos];
    PMEFloat svx  = vx + lvx;
    PMEFloat svy  = vy + lvy;
    PMEFloat svz  = vz + lvz;
    eke          += mass * (svx * svx + svy * svy + svz * svz);
    ekpbs        += mass * (vx * lvx + vy * lvy + vz * lvz);
    ekph         += mass * (vx * vx + vy * vy + vz * vz);
    pos          += blockDim.x * gridDim.x;
  }
  eke   *= (PMEFloat)0.125 * c_ave;
  ekph  *= (PMEFloat)0.5;
  ekpbs *= (PMEFloat)0.5;
  sE[threadIdx.x].KE.EKE   = eke;
  sE[threadIdx.x].KE.EKPH  = ekph;
  sE[threadIdx.x].KE.EKPBS = ekpbs;

  // Reduce per-thread kinetic energies
  __syncthreads();
  for (uint stride = blockDim.x >> 1; stride > 0; stride >>= 1) {
    if (threadIdx.x < stride) {
      sE[threadIdx.x].KE.EKE += sE[threadIdx.x + stride].KE.EKE;
      sE[threadIdx.x].KE.EKPH += sE[threadIdx.x + stride].KE.EKPH;
      sE[threadIdx.x].KE.EKPBS += sE[threadIdx.x + stride].KE.EKPBS;
    }
    __syncthreads();
  }

  // Save result
  if (threadIdx.x < 3) {
    cSim.pKineticEnergy[blockIdx.x].array[threadIdx.x] = sE[0].array[threadIdx.x];
  }
}

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(THREADS_PER_BLOCK, 1)
kCalculateKineticEnergyAFE_kernel(PMEFloat c_ave)
{
  extern __shared__ AFEKineticEnergy sAFEE[];

  PMEFloat eke          = (PMEFloat)0.0;
  PMEFloat ekpbs        = (PMEFloat)0.0;
  PMEFloat ekph         = (PMEFloat)0.0;
  PMEFloat ti_eke_R1    = (PMEFloat)0.0;
  PMEFloat ti_eke_R2    = (PMEFloat)0.0;
  PMEFloat ti_ekpbs_R1  = (PMEFloat)0.0;
  PMEFloat ti_ekpbs_R2  = (PMEFloat)0.0;
  PMEFloat ti_ekph_R1   = (PMEFloat)0.0;
  PMEFloat ti_ekph_R2   = (PMEFloat)0.0;
  PMEFloat ti_sc_eke_R1 = (PMEFloat)0.0;
  PMEFloat ti_sc_eke_R2 = (PMEFloat)0.0;
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

  // Sum up kinetic energies
  while (pos < cSim.atoms) {
    bool isSC     = (cSim.pImageTIRegion[pos] & 0x1);
    int TIRegion  = (cSim.pImageTIRegion[pos] >> 1);
    PMEFloat mass = cSim.pAtomMass[pos];
    PMEFloat vx   = cSim.pVelX[pos];
    PMEFloat lvx  = cSim.pLVelX[pos];
    PMEFloat vy   = cSim.pVelY[pos];
    PMEFloat lvy  = cSim.pLVelY[pos];
    PMEFloat vz   = cSim.pVelZ[pos];
    PMEFloat lvz  = cSim.pLVelZ[pos];
    PMEFloat svx  = vx + lvx;
    PMEFloat svy  = vy + lvy;
    PMEFloat svz  = vz + lvz;
    if (TIRegion == 1) {
      ti_eke_R1   += mass * (svx*svx + svy*svy + svz*svz);
      ti_ekpbs_R1 += mass * ( vx*lvx +  vy*lvy +  vz*lvz);
      ti_ekph_R1  += mass * ( vx*vx  +  vy*vy  +  vz*vz );
      if (isSC) {
        ti_sc_eke_R1 += mass * (svx*svx + svy*svy + svz*svz);
      }
    }
    if (TIRegion == 2) {
      ti_eke_R2   += mass * (svx*svx + svy*svy + svz*svz);
      ti_ekpbs_R2 += mass * ( vx*lvx +  vy*lvy +  vz*lvz);
      ti_ekph_R2  += mass * ( vx*vx  +  vy*vy  +  vz*vz );
      if (isSC) {
        ti_sc_eke_R2 += mass * (svx*svx + svy*svy + svz*svz);
      }
    }
    eke   += mass * (svx*svx + svy*svy + svz*svz);
    ekpbs += mass * ( vx*lvx +  vy*lvy +  vz*lvz);
    ekph  += mass * ( vx*vx  +  vy*vy  +  vz*vz );
    pos   += blockDim.x * gridDim.x;
  }
  eke          *= (PMEFloat)0.125 * c_ave;
  ti_eke_R1    *= (PMEFloat)0.125 * c_ave;
  ti_eke_R2    *= (PMEFloat)0.125 * c_ave;
  ekph         *= (PMEFloat)0.5;
  ti_ekph_R1   *= (PMEFloat)0.5;
  ti_ekph_R2   *= (PMEFloat)0.5;
  ekpbs        *= (PMEFloat)0.5;
  ti_ekpbs_R1  *= (PMEFloat)0.5;
  ti_ekpbs_R2  *= (PMEFloat)0.5;
  ti_sc_eke_R1 *= (PMEFloat)0.125 * c_ave;
  ti_sc_eke_R2 *= (PMEFloat)0.125 * c_ave;

  // We need to calculate the KE for essentially two different molecules:
  // The molecules share a common core, but in order to get the unique atom
  // contributions we subtract out the contributions from the common region.
  sAFEE[threadIdx.x].AFEKE.TI_EKER1    = eke - ti_eke_R2;
  sAFEE[threadIdx.x].AFEKE.TI_EKER2    = eke - ti_eke_R1;
  sAFEE[threadIdx.x].AFEKE.TI_EKPHR1   = ekph - ti_ekph_R2;
  sAFEE[threadIdx.x].AFEKE.TI_EKPHR2   = ekph - ti_ekph_R1;
  sAFEE[threadIdx.x].AFEKE.TI_EKPBSR1  = ekpbs - ti_ekpbs_R2;
  sAFEE[threadIdx.x].AFEKE.TI_EKPBSR2  = ekpbs - ti_ekpbs_R1;
  sAFEE[threadIdx.x].AFEKE.TI_SC_EKER1 = ti_sc_eke_R1;
  sAFEE[threadIdx.x].AFEKE.TI_SC_EKER2 = ti_sc_eke_R2;

  // Reduce per-thread kinetic energies
  __syncthreads();
  for (uint stride = blockDim.x >> 1; stride > 0; stride >>= 1) {
    if (threadIdx.x < stride) {
      sAFEE[threadIdx.x].AFEKE.TI_EKER1 += sAFEE[threadIdx.x + stride].AFEKE.TI_EKER1;
      sAFEE[threadIdx.x].AFEKE.TI_EKER2 += sAFEE[threadIdx.x + stride].AFEKE.TI_EKER2;
      sAFEE[threadIdx.x].AFEKE.TI_EKPHR1 += sAFEE[threadIdx.x + stride].AFEKE.TI_EKPHR1;
      sAFEE[threadIdx.x].AFEKE.TI_EKPHR2 += sAFEE[threadIdx.x + stride].AFEKE.TI_EKPHR2;
      sAFEE[threadIdx.x].AFEKE.TI_EKPBSR1 += sAFEE[threadIdx.x + stride].AFEKE.TI_EKPBSR1;
      sAFEE[threadIdx.x].AFEKE.TI_EKPBSR2 += sAFEE[threadIdx.x + stride].AFEKE.TI_EKPBSR2;
      sAFEE[threadIdx.x].AFEKE.TI_SC_EKER1 += sAFEE[threadIdx.x + stride].AFEKE.TI_SC_EKER1;
      sAFEE[threadIdx.x].AFEKE.TI_SC_EKER2 += sAFEE[threadIdx.x + stride].AFEKE.TI_SC_EKER2;
    }
    __syncthreads();
  }

  // Save result
  if (threadIdx.x < 8) {
    cSim.pAFEKineticEnergy[blockIdx.x].array[threadIdx.x] = sAFEE[0].array[threadIdx.x];
  }
}

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(THREADS_PER_BLOCK, 1)
kPMECalculateKineticEnergyAFE_kernel(PMEFloat c_ave)
{
  extern __shared__ AFEKineticEnergy sAFEE[];

  PMEFloat eke          = (PMEFloat)0.0;
  PMEFloat ekpbs        = (PMEFloat)0.0;
  PMEFloat ekph         = (PMEFloat)0.0;
  PMEFloat ti_eke_R1    = (PMEFloat)0.0;
  PMEFloat ti_eke_R2    = (PMEFloat)0.0;
  PMEFloat ti_ekpbs_R1  = (PMEFloat)0.0;
  PMEFloat ti_ekpbs_R2  = (PMEFloat)0.0;
  PMEFloat ti_ekph_R1   = (PMEFloat)0.0;
  PMEFloat ti_ekph_R2   = (PMEFloat)0.0;
  PMEFloat ti_sc_eke_R1 = (PMEFloat)0.0;
  PMEFloat ti_sc_eke_R2 = (PMEFloat)0.0;
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  while (pos < cSim.atoms) {
    bool isSC     = (cSim.pImageTIRegion[pos] & 0x1);
    int TIRegion  = (cSim.pImageTIRegion[pos] >> 1);
    PMEFloat mass = cSim.pImageMass[pos];
    PMEFloat vx   = cSim.pImageVelX[pos];
    PMEFloat lvx  = cSim.pImageLVelX[pos];
    PMEFloat vy   = cSim.pImageVelY[pos];
    PMEFloat lvy  = cSim.pImageLVelY[pos];
    PMEFloat vz   = cSim.pImageVelZ[pos];
    PMEFloat lvz  = cSim.pImageLVelZ[pos];
    PMEFloat svx  = vx + lvx;
    PMEFloat svy  = vy + lvy;
    PMEFloat svz  = vz + lvz;
    if (TIRegion == 1) {
      ti_eke_R1   += mass * (svx*svx + svy*svy + svz*svz);
      ti_ekpbs_R1 += mass * ( vx*lvx +  vy*lvy +  vz*lvz);
      ti_ekph_R1  += mass * ( vx*vx  +  vy*vy  +  vz*vz );
      if (isSC) {
        ti_sc_eke_R1 += mass * (svx*svx + svy*svy + svz*svz);
      }
    }
    if (TIRegion == 2) {
      ti_eke_R2   += mass * (svx*svx + svy*svy + svz*svz);
      ti_ekpbs_R2 += mass * ( vx*lvx +  vy*lvy +  vz*lvz);
      ti_ekph_R2  += mass * ( vx*vx  +  vy*vy  +  vz*vz );
      if (isSC) {
        ti_sc_eke_R2 += mass * (svx*svx + svy*svy + svz*svz);
      }
    }
    eke   += mass * (svx*svx + svy*svy + svz*svz);
    ekpbs += mass * ( vx*lvx +  vy*lvy +  vz*lvz);
    ekph  += mass * ( vx*vx  +  vy*vy  +  vz*vz );
    pos   += blockDim.x * gridDim.x;
  }
  eke          *= (PMEFloat)0.125 * c_ave;
  ti_eke_R1    *= (PMEFloat)0.125 * c_ave;
  ti_eke_R2    *= (PMEFloat)0.125 * c_ave;
  ekph         *= (PMEFloat)0.5;
  ti_ekph_R1   *= (PMEFloat)0.5;
  ti_ekph_R2   *= (PMEFloat)0.5;
  ekpbs        *= (PMEFloat)0.5;
  ti_ekpbs_R1  *= (PMEFloat)0.5;
  ti_ekpbs_R2  *= (PMEFloat)0.5;
  ti_sc_eke_R1 *= (PMEFloat)0.125 * c_ave;
  ti_sc_eke_R2 *= (PMEFloat)0.125 * c_ave;

  // We need to calculate the KE for essentially two different molecules.
  // The molecules share a common core, but in order to get the unique atom
  // contributions we subtract out the contributions from the common region.
  sAFEE[threadIdx.x].AFEKE.TI_EKER1    = eke - ti_eke_R2;
  sAFEE[threadIdx.x].AFEKE.TI_EKER2    = eke - ti_eke_R1;
  sAFEE[threadIdx.x].AFEKE.TI_EKPHR1   = ekph - ti_ekph_R2;
  sAFEE[threadIdx.x].AFEKE.TI_EKPHR2   = ekph - ti_ekph_R1;
  sAFEE[threadIdx.x].AFEKE.TI_EKPBSR1  = ekpbs - ti_ekpbs_R2;
  sAFEE[threadIdx.x].AFEKE.TI_EKPBSR2  = ekpbs - ti_ekpbs_R1;
  sAFEE[threadIdx.x].AFEKE.TI_SC_EKER1 = ti_sc_eke_R1;
  sAFEE[threadIdx.x].AFEKE.TI_SC_EKER2 = ti_sc_eke_R2;

  // Reduce per-thread kinetic energies
  __syncthreads();
  for (uint stride = blockDim.x >> 1; stride > 0; stride >>= 1) {
    if (threadIdx.x < stride) {
      sAFEE[threadIdx.x].AFEKE.TI_EKER1 += sAFEE[threadIdx.x + stride].AFEKE.TI_EKER1;
      sAFEE[threadIdx.x].AFEKE.TI_EKER2 += sAFEE[threadIdx.x + stride].AFEKE.TI_EKER2;
      sAFEE[threadIdx.x].AFEKE.TI_EKPHR1 += sAFEE[threadIdx.x + stride].AFEKE.TI_EKPHR1;
      sAFEE[threadIdx.x].AFEKE.TI_EKPHR2 += sAFEE[threadIdx.x + stride].AFEKE.TI_EKPHR2;
      sAFEE[threadIdx.x].AFEKE.TI_EKPBSR1 += sAFEE[threadIdx.x + stride].AFEKE.TI_EKPBSR1;
      sAFEE[threadIdx.x].AFEKE.TI_EKPBSR2 += sAFEE[threadIdx.x + stride].AFEKE.TI_EKPBSR2;
      sAFEE[threadIdx.x].AFEKE.TI_SC_EKER1 += sAFEE[threadIdx.x + stride].AFEKE.TI_SC_EKER1;
      sAFEE[threadIdx.x].AFEKE.TI_SC_EKER2 += sAFEE[threadIdx.x + stride].AFEKE.TI_SC_EKER2;
    }
    __syncthreads();
  }

  // Save result
  if (threadIdx.x < 8) {
    cSim.pAFEKineticEnergy[blockIdx.x].array[threadIdx.x] = sAFEE[0].array[threadIdx.x];
  }
}

//---------------------------------------------------------------------------------------------
// kCalculateKineticEnergy: compute the kinetic energy using the appropriate kernel for
//                          non-AFE applications.  This is separate from its AFE variant
//                          because there are separate host functions in gpu.cpp calling each
//                          one.
//
// Arguments:
//   gpu:      overarching type for storing parameters, coordinates, and the energy function
//   c_ave:    composite quantity referenced in Langevin dynamics
//---------------------------------------------------------------------------------------------
void kCalculateKineticEnergy(gpuContext gpu, PMEFloat c_ave)
{
  if (gpu->bNeighborList) {
    kPMECalculateKineticEnergy_kernel<<<gpu->blocks, gpu->threadsPerBlock,
                                        gpu->threadsPerBlock * sizeof(KineticEnergy)>>>(c_ave);
  }
  else {
    kCalculateKineticEnergy_kernel<<<gpu->blocks, gpu->threadsPerBlock,
                                     gpu->threadsPerBlock * sizeof(KineticEnergy)>>>(c_ave);
  }
  LAUNCHERROR("kCalculateKineticEnergy");
}

//---------------------------------------------------------------------------------------------
// kCalculateKineticEnergyAFE: compute kinetic energy for alchemical free energy simulations.
//
// Arguments:
//   gpu:      overarching type for storing parameters, coordinates, and the energy function
//   c_ave:    composite quantity referenced in Langevin dynamics
//---------------------------------------------------------------------------------------------
void kCalculateKineticEnergyAFE(gpuContext gpu, PMEFloat c_ave)
{
  if (gpu->bNeighborList) {
#ifdef use_DPFP
    // Use half as many threads: the dynamically allocated __shared__ array would
    // overload the hardware limits otherwise.
    kPMECalculateKineticEnergyAFE_kernel<<<gpu->blocks, gpu->threadsPerBlock/2,
                                           gpu->threadsPerBlock/2 *
                                           sizeof(AFEKineticEnergy)>>>(c_ave);
#else
    kPMECalculateKineticEnergyAFE_kernel<<<gpu->blocks, gpu->threadsPerBlock,
                                           gpu->threadsPerBlock *
                                           sizeof(AFEKineticEnergy)>>>(c_ave);
#endif
  }
  else {
#ifdef use_DPFP
    // Use half as many threads: the dynamically allocated __shared__ array would
    // overload the hardware limits otherwise.
    kCalculateKineticEnergyAFE_kernel<<<gpu->blocks * 2, gpu->threadsPerBlock/2,
                                        gpu->threadsPerBlock/2 *
                                        sizeof(AFEKineticEnergy)>>>(c_ave);
#else
    kCalculateKineticEnergyAFE_kernel<<<gpu->blocks, gpu->threadsPerBlock,
                                        gpu->threadsPerBlock *
                                        sizeof(AFEKineticEnergy)>>>(c_ave);
#endif
  }
  LAUNCHERROR("kCalculateKineticEnergyAFE");
}

//---------------------------------------------------------------------------------------------
// Kernels for velocity scaling in the context of constant temperature simulations.  The only
// distinction is whether PME (a neighbor list framework) is in effect.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kScaleVelocities_kernel(PMEDouble scale)
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;

  if (pos < cSim.atoms) {
    double vx = cSim.pVelX[pos];
    double vy = cSim.pVelY[pos];
    double vz = cSim.pVelZ[pos];
    vx *= scale;
    vy *= scale;
    vz *= scale;
    cSim.pVelX[pos] = vx;
    cSim.pVelY[pos] = vy;
    cSim.pVelZ[pos] = vz;
  }
}

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kPMEScaleVelocities_kernel(PMEDouble scale)
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;

  if (pos < cSim.atoms) {
    double vx = cSim.pImageVelX[pos];
    double vy = cSim.pImageVelY[pos];
    double vz = cSim.pImageVelZ[pos];
    vx *= scale;
    vy *= scale;
    vz *= scale;
    cSim.pImageVelX[pos] = vx;
    cSim.pImageVelY[pos] = vy;
    cSim.pImageVelZ[pos] = vz;
  }
}

//---------------------------------------------------------------------------------------------
// kScaleVelocities: host function to call the appropriate kernel for velocity rescaling.
//
// Arugments:
//   gpu:    overarching type for storing parameters, coordinates, and the energy function
//   scale:  scaling factor to apply
//---------------------------------------------------------------------------------------------
void kScaleVelocities(gpuContext gpu, PMEDouble scale)
{
  if (gpu->bNeighborList) {
    kPMEScaleVelocities_kernel<<<gpu->updateBlocks, gpu->updateThreadsPerBlock>>>(scale);
  }
  else {
    kScaleVelocities_kernel<<<gpu->updateBlocks, gpu->updateThreadsPerBlock>>>(scale);
  }
  LAUNCHERROR("kScaleVelocities");
}


//---------------------------------------------------------------------------------------------
// Kernels for velocity scaling in the context of constant temperature simulations.  The only
// distinction is whether PME (a neighbor list framework) is in effect.
//---------------------------------------------------------------------------------------------
__global__ void
__launch_bounds__(UPDATE_THREADS_PER_BLOCK, 1)
kScaleSGLD_kernel(PMEDouble scale,PMEDouble scalsg)
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;

  if (pos < cSim.atoms) {
    double vx = cSim.pVelX[pos];
    double vy = cSim.pVelY[pos];
    double vz = cSim.pVelZ[pos];
    vx *= scale;
    vy *= scale;
    vz *= scale;
    cSim.pVelX[pos] = vx;
    cSim.pVelY[pos] = vy;
    cSim.pVelZ[pos] = vz;
    double x0 = cSim.pX0sg[pos];
    double y0 = cSim.pY0sg[pos];
    double z0 = cSim.pZ0sg[pos];
    double x1 = cSim.pX1sg[pos];
    double y1 = cSim.pY1sg[pos];
    double z1 = cSim.pZ1sg[pos];
    x1 = x0 + scalsg*(x1-x0);
    y1 = y0 + scalsg*(y1-y0);
    z1 = z0 + scalsg*(z1-z0);
    cSim.pX1sg[pos] = x1;
    cSim.pY1sg[pos] = y1;
    cSim.pZ1sg[pos] = z1;
    double x2 = cSim.pX2sg[pos];
    double y2 = cSim.pY2sg[pos];
    double z2 = cSim.pZ2sg[pos];
    x2 = x0 + scalsg*(x2-x0);
    y2 = y0 + scalsg*(y2-y0);
    z2 = z0 + scalsg*(z2-z0);
    cSim.pX2sg[pos] = x2;
    cSim.pY2sg[pos] = y2;
    cSim.pZ2sg[pos] = z2;
    double rsgx = cSim.pRsgX[pos];
    double rsgy = cSim.pRsgY[pos];
    double rsgz = cSim.pRsgZ[pos];
    rsgx *= scalsg;
    rsgy *= scalsg;
    rsgz *= scalsg;
    cSim.pRsgX[pos] = rsgx;
    cSim.pRsgY[pos] = rsgy;
    cSim.pRsgZ[pos] = rsgz;
    double fpsg = cSim.pFPsg[pos];
    fpsg *= scalsg;
    cSim.pFPsg[pos] = fpsg;
    double ppsg = cSim.pPPsg[pos];
    ppsg *= scalsg*scalsg;
    cSim.pPPsg[pos] = ppsg;
  }
}

//---------------------------------------------------------------------------------------------
__global__ void
__launch_bounds__(UPDATE_THREADS_PER_BLOCK, 1)
kPMEScaleSGLD_kernel(PMEDouble scale,PMEDouble scalsg)
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;

  if (pos < cSim.atoms) {
    double vx = cSim.pImageVelX[pos];
    double vy = cSim.pImageVelY[pos];
    double vz = cSim.pImageVelZ[pos];
    vx *= scale;
    vy *= scale;
    vz *= scale;
    cSim.pImageVelX[pos] = vx;
    cSim.pImageVelY[pos] = vy;
    cSim.pImageVelZ[pos] = vz;
    double x0 = cSim.pImageX[pos];
    double y0 = cSim.pImageY[pos];
    double z0 = cSim.pImageZ[pos];
    double x1 = cSim.pImageX1sg[pos];
    double y1 = cSim.pImageY1sg[pos];
    double z1 = cSim.pImageZ1sg[pos];
    x1 = x0 + scalsg*(x1-x0);
    y1 = y0 + scalsg*(y1-y0);
    z1 = z0 + scalsg*(z1-z0);
    cSim.pImageX1sg[pos] = x1;
    cSim.pImageY1sg[pos] = y1;
    cSim.pImageZ1sg[pos] = z1;
    double x2 = cSim.pImageX2sg[pos];
    double y2 = cSim.pImageY2sg[pos];
    double z2 = cSim.pImageZ2sg[pos];
    x2 = x0 + scalsg*(x2-x0);
    y2 = y0 + scalsg*(y2-y0);
    z2 = z0 + scalsg*(z2-z0);
    cSim.pImageX2sg[pos] = x2;
    cSim.pImageY2sg[pos] = y2;
    cSim.pImageZ2sg[pos] = z2;
    double rsgx = cSim.pImageRsgX[pos];
    double rsgy = cSim.pImageRsgY[pos];
    double rsgz = cSim.pImageRsgZ[pos];
    rsgx *= scalsg;
    rsgy *= scalsg;
    rsgz *= scalsg;
    cSim.pImageRsgX[pos] = rsgx;
    cSim.pImageRsgY[pos] = rsgy;
    cSim.pImageRsgZ[pos] = rsgz;
    double fpsg = cSim.pImageFPsg[pos];
    fpsg *= scalsg;
    cSim.pImageFPsg[pos] = fpsg;
    double ppsg = cSim.pImagePPsg[pos];
    ppsg *= scalsg*scalsg;
    cSim.pImagePPsg[pos] = ppsg;
  }
}

//---------------------------------------------------------------------------------------------
// kScaleSGLD: host function to call the appropriate kernel for SGLD rescaling.
//
// Arugments:
//   gpu:    overarching type for storing parameters, coordinates, and the energy function
//   scale:  scaling factor to apply
//---------------------------------------------------------------------------------------------
void kScaleSGLD(gpuContext gpu, PMEDouble scale, PMEDouble scalsg)
{
  if (gpu->bNeighborList) {
    kPMEScaleSGLD_kernel<<<gpu->updateBlocks, gpu->updateThreadsPerBlock>>>(scale,scalsg);
  }
  else {
    kScaleSGLD_kernel<<<gpu->updateBlocks, gpu->updateThreadsPerBlock>>>(scale,scalsg);
  }
  LAUNCHERROR("kScaleSGLD");
}

//---------------------------------------------------------------------------------------------
// Two inclusions of kNTPKernels.h to perform constant pressure simulations with different
// numbers of molecules in the system.
//---------------------------------------------------------------------------------------------
#define KPMECALCULATECOM_KERNEL kPMECalculateCOM_kernel
#define KPMECALCULATESOLUTECOM_KERNEL kPMECalculateSoluteCOM_kernel
#define KPMECALCULATECOMKINETICENERGY_KERNEL kPMECalculateCOMKineticEnergy_kernel
#define KCALCULATEMOLECULARVIRIAL_KERNEL kCalculateMolecularVirial_kernel
#define KPRESSURESCALECOORDINATES_KERNEL kPressureScaleCoordinates_kernel
#define KPMECALCULATECOMAFE_KERNEL kPMECalculateCOMAFE_kernel
#include "kNTPKernels.h"
#undef KPMECALCULATECOM_KERNEL
#undef KPMECALCULATESOLUTECOM_KERNEL
#undef KPMECALCULATECOMKINETICENERGY_KERNEL
#undef KCALCULATEMOLECULARVIRIAL_KERNEL
#undef KPRESSURESCALECOORDINATES_KERNEL
#undef KPMECALCULATECOMAFE_KERNEL

//---------------------------------------------------------------------------------------------
#define NTP_LOTSOFMOLECULES
#define KPMECALCULATECOM_KERNEL kPMECalculateCOMLarge_kernel
#define KPMECALCULATESOLUTECOM_KERNEL kPMECalculateSoluteCOMLarge_kernel
#define KPMECALCULATECOMKINETICENERGY_KERNEL kPMECalculateCOMKineticEnergyLarge_kernel
#define KCALCULATEMOLECULARVIRIAL_KERNEL kCalculateMolecularVirialLarge_kernel
#define KPRESSURESCALECOORDINATES_KERNEL kPressureScaleCoordinatesLarge_kernel
#define KPMECALCULATECOMAFE_KERNEL kPMECalculateCOMLargeAFE_kernel
#include "kNTPKernels.h"
#undef KPMECALCULATECOM_KERNEL
#undef KPMECALCULATESOLUTECOM_KERNEL
#undef KPMECALCULATECOMKINETICENERGY_KERNEL
#undef KCALCULATEMOLECULARVIRIAL_KERNEL
#undef KPRESSURESCALECOORDINATES_KERNEL
#undef KPMECALCULATECOMAFE_KERNEL
#undef NTP_LOTSOFMOLECULES

#include "kNTPCalls.h"
#include "kRandom.h"
#include "kAFECalls.h"

//---------------------------------------------------------------------------------------------
// Kernels for handling extra points in the context of implicit solvent / isolated systems, a
// neighbor list framework, or constant pressure conditions.
//---------------------------------------------------------------------------------------------
#define EP_NEIGHBORLIST
__global__ void
__LAUNCH_BOUNDS__(VS_THREADS_PER_BLOCK, 1)
kNLOrientForces_kernel()
#include "kOrientForcesKernel.h"

//---------------------------------------------------------------------------------------------
#define EP_VIRIAL
__global__ void
__LAUNCH_BOUNDS__(VS_THREADS_PER_BLOCK, 1)
kNLOrientForcesVirial_kernel()
#include "kOrientForcesKernel.h"
#undef EP_VIRIAL
#undef EP_NEIGHBORLIST

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(VS_THREADS_PER_BLOCK, 1)
kOrientForces_kernel()
#include "kOrientForcesKernel.h"

//---------------------------------------------------------------------------------------------
// kOrientForces: host function to launch the appropriate kernel for handling extra point
//                force transmission (to put the forces on frame atoms which have mass).
//
// Arugments:
//   gpu:    overarching type for storing parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kOrientForces(gpuContext gpu)
{
  if (gpu->bNeighborList) {
    if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
      kNLOrientForcesVirial_kernel<<<gpu->blocks * VS_BLOCKS_MULTIPLIER,
                                     VS_THREADS_PER_BLOCK>>>();
    }
    else {
      kNLOrientForces_kernel<<<gpu->blocks * VS_BLOCKS_MULTIPLIER, VS_THREADS_PER_BLOCK>>>();
    }
  }
  else {
    kOrientForces_kernel<<<gpu->blocks * VS_BLOCKS_MULTIPLIER, VS_THREADS_PER_BLOCK>>>();
  }
  LAUNCHERROR("kOrientForces");
}

//---------------------------------------------------------------------------------------------
#define EP_NEIGHBORLIST
__global__ void
__LAUNCH_BOUNDS__(VS_THREADS_PER_BLOCK, VS_BLOCKS_MINIMUM)
kNLLocalToGlobal_kernel()
#include "kLocalToGlobalKernel.h"
#undef EP_NEIGHBORLIST

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(VS_THREADS_PER_BLOCK, VS_BLOCKS_MINIMUM)
kLocalToGlobal_kernel()
#include "kLocalToGlobalKernel.h"

//---------------------------------------------------------------------------------------------
// kLocalToGlobal: host function to call the appropriate kernel for transmitting EP forces to
//                 frame atoms.
//
// Arugments:
//   gpu:    overarching type for storing parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kLocalToGlobal(gpuContext gpu)
{
  if (gpu->bNeighborList) {
    kNLLocalToGlobal_kernel<<<gpu->blocks * VS_BLOCKS_MULTIPLIER, VS_THREADS_PER_BLOCK>>>();
  }
  else {
    kLocalToGlobal_kernel<<<gpu->blocks * VS_BLOCKS_MULTIPLIER, VS_THREADS_PER_BLOCK>>>();
  }
  LAUNCHERROR("kLocalToGlobal");
}

//---------------------------------------------------------------------------------------------
// kScaledMDScaleForces_kernel:
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kScaledMDScaleForces_kernel(PMEDouble pot_ene_tot, PMEDouble lambda)
{
  // Scale forces ScaledMD
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  while (pos < cSim.atoms) {
    PMEDouble forceX = cSim.pForceXAccumulator[pos];
    PMEDouble forceY = cSim.pForceYAccumulator[pos];
    PMEDouble forceZ = cSim.pForceZAccumulator[pos];
    forceX *= lambda;
    forceY *= lambda;
    forceZ *= lambda;
    cSim.pForceXAccumulator[pos] = forceX;
    cSim.pForceYAccumulator[pos] = forceY;
    cSim.pForceZAccumulator[pos] = forceZ;
    pos += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kScaledMDScaleForces: launch the one kernel to scale forces for Scaled MD applications.
//
// Arguments:
//   gpu:          overarching type for storing parameters, coordinates, and the energy
//                 function
//   pot_ene_tot:  the total (unscaled) potential energy of the system
//---------------------------------------------------------------------------------------------
void kScaledMDScaleForces(gpuContext gpu, PMEDouble pot_ene_tot, PMEDouble lambda)
{
  kScaledMDScaleForces_kernel<<<gpu->blocks,
                                gpu->generalThreadsPerBlock>>>(pot_ene_tot, lambda);
  LAUNCHERROR("kScaledMDScaleForces");
}

//---------------------------------------------------------------------------------------------
// kCheckGpuForces_kernel: check forces for each atom.
//---------------------------------------------------------------------------------------------
__global__ void __LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
  kCheckGpuForces_kernel(int *result, PMEFloat4 *frc)
{
  __shared__ volatile int nviolate, nNBviolate;
  __shared__ volatile PMEFloat sMaxFx[GRID];
  __shared__ volatile PMEFloat sMaxFy[GRID];
  __shared__ volatile PMEFloat sMaxFz[GRID];

  PMEFloat maxfx = 0.0;
  PMEFloat maxfy = 0.0;
  PMEFloat maxfz = 0.0;
  if (threadIdx.x == 0) {
    nviolate = 0;
    nNBviolate = 0;
  }
  __syncthreads();
  int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (cSim.ntp > 0 && cSim.barostat == 1) {
    while (pos < cSim.atoms) {
      PMEAccumulator llifx = cSim.pForceXAccumulator[pos];
      PMEAccumulator llify = cSim.pForceYAccumulator[pos];
      PMEAccumulator llifz = cSim.pForceZAccumulator[pos];
      PMEFloat fx = (PMEFloat)llifx / FORCESCALEF;
      PMEFloat fy = (PMEFloat)llify / FORCESCALEF;
      PMEFloat fz = (PMEFloat)llifz / FORCESCALEF;
      bool signal = false;
      if (fabs(fx) > (PMEFloat)500.0) {
        maxfx = max(maxfx, fabs(fx));
        atomicAdd((int*)&nviolate, 1);
	signal = true;
      }
      if (fabs(fy) > (PMEFloat)500.0) {
        maxfy = max(maxfy, fabs(fy));
        atomicAdd((int*)&nviolate, 1);
	signal = true;
      }
      if (fabs(fz) > (PMEFloat)500.0) {
        maxfz = max(maxfz, fabs(fz));
        atomicAdd((int*)&nviolate, 1);
	signal = true;
      }
      if (signal) {
	printf("Atom %8d    violation :: %14.4f %14.4f %14.4f\n", pos, fx, fy, fz);
      }
      llifx = cSim.pNBForceXAccumulator[pos];
      llify = cSim.pNBForceYAccumulator[pos];
      llifz = cSim.pNBForceZAccumulator[pos];
      fx = (PMEFloat)llifx / FORCESCALEF;
      fy = (PMEFloat)llify / FORCESCALEF;
      fz = (PMEFloat)llifz / FORCESCALEF;
      signal = false;
      if (fabs(fx) > (PMEFloat)500.0) {
        maxfx = max(maxfx, fabs(fx));
        atomicAdd((int*)&nNBviolate, 1);
	signal = true;
      }
      if (fabs(fy) > (PMEFloat)500.0) {
        maxfy = max(maxfy, fabs(fy));
        atomicAdd((int*)&nNBviolate, 1);
	signal = true;
      }
      if (fabs(fz) > (PMEFloat)500.0) {
        maxfz = max(maxfz, fabs(fz));
        atomicAdd((int*)&nNBviolate, 1);
	signal = true;
      }
      if (signal) {
	printf("Atom %8d NB violation :: %14.4f %14.4f %14.4f\n", pos, fx, fy, fz);
      }
      pos += gridDim.x * blockDim.x;
    }
  }
  else {
    while (pos < cSim.atoms) {
      PMEAccumulator llifx = cSim.pForceXAccumulator[pos];
      PMEAccumulator llify = cSim.pForceYAccumulator[pos];
      PMEAccumulator llifz = cSim.pForceZAccumulator[pos];
      PMEFloat fx = (PMEFloat)llifx / FORCESCALEF;
      PMEFloat fy = (PMEFloat)llify / FORCESCALEF;
      PMEFloat fz = (PMEFloat)llifz / FORCESCALEF;
      bool signal = false;
      if (fabs(fx) > (PMEFloat)500.0) {
        maxfx = max(maxfx, fabs(fx));
        atomicAdd((int*)&nviolate, 1);
	signal = true;
      }
      if (fabs(fy) > (PMEFloat)500.0) {
        maxfy = max(maxfy, fabs(fy));
        atomicAdd((int*)&nviolate, 1);
	signal = true;
      }
      if (fabs(fz) > (PMEFloat)500.0) {
        maxfz = max(maxfz, fabs(fz));
        atomicAdd((int*)&nviolate, 1);
	signal = true;
      }
      if (signal) {
	printf("Atom %8d    violation :: %14.4f %14.4f %14.4f\n", pos, fx, fy, fz);
      }
      pos += gridDim.x * blockDim.x;
    }
  }
  __syncthreads();
  for (unsigned int stride = warpSize >> 1; stride > 0; stride >>=1) {
    maxfx = max(maxfx, __SHFL_DOWN(WARP_MASK, maxfx, stride));
    maxfy = max(maxfy, __SHFL_DOWN(WARP_MASK, maxfy, stride));
    maxfz = max(maxfz, __SHFL_DOWN(WARP_MASK, maxfz, stride));
  }
  if ((threadIdx.x & GRID_BITS_MASK) == 0) {
    sMaxFx[threadIdx.x >> GRID_BITS] = maxfx;
    sMaxFy[threadIdx.x >> GRID_BITS] = maxfy;
    sMaxFz[threadIdx.x >> GRID_BITS] = maxfz;
  }
  __syncthreads();

  if (threadIdx.x < GRID && (nviolate > 0 || nNBviolate > 0)) {
    maxfx = sMaxFx[threadIdx.x];
    maxfy = sMaxFy[threadIdx.x];
    maxfz = sMaxFz[threadIdx.x];
    //There is no divergence within a warp, so using WARP_MASK mask
    for (unsigned int stride = warpSize >> 1; stride > 0; stride >>=1) {
      maxfx = max(maxfx, __SHFL_DOWN(WARP_MASK, maxfx, stride));
      maxfy = max(maxfy, __SHFL_DOWN(WARP_MASK, maxfy, stride));
      maxfz = max(maxfz, __SHFL_DOWN(WARP_MASK, maxfz, stride));
    }
  }
  if (threadIdx.x == 0 && (nviolate > 0 || nNBviolate > 0)) {
    atomicAdd(result, nviolate + nNBviolate);
    PMEFloat4 maxf;
    maxf.x = maxfx;
    maxf.y = maxfy;
    maxf.z = maxfz;
    *frc = maxf;
    printf("%4d violations :: %14.4f %14.4f %14.4f max\n", nviolate + nNBviolate, maxfx,
	   maxfy, maxfz);
  }
}

//---------------------------------------------------------------------------------------------
// kCheckGpuForces: launch the kernel to check that all forces on the GPU are within
//                  reasonable bounds.
//
// This is a debugging function.
//
// Arguments:
//   gpu:          overarching type for storing parameters, coordinates, and the energy
//                 function
//   frc_chk:      indication of the stage at which the force is being checked
//---------------------------------------------------------------------------------------------
void kCheckGpuForces(gpuContext gpu, int frc_chk)
{
  int hsignal = 0;
  int* dsignal;
  PMEFloat4 hfrc;
  PMEFloat4 *dfrc;
  cudaMalloc(&dsignal, sizeof(int));
  cudaMalloc(&dfrc, sizeof(PMEFloat4));
  cudaMemcpy(dsignal, &hsignal, sizeof(int), cudaMemcpyHostToDevice);
  kCheckGpuForces_kernel<<<gpu->blocks, gpu->generalThreadsPerBlock>>>(dsignal, dfrc);
  cudaMemcpy(&hsignal, dsignal, sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(&hfrc, dfrc, sizeof(PMEFloat4), cudaMemcpyDeviceToHost);
  cudaFree(dsignal);
  cudaFree(dfrc);
  cudaDeviceSynchronize();
  if (hsignal > 0) {
    printf("Signal detected after stage %2d, step %4d : %6d violations.\n", frc_chk, gpu->step,
           hsignal);
    exit(1);
  }
}

//---------------------------------------------------------------------------------------------
// kCheckGpuConsistency: launch a kernel to check the forces, virial, and energy.  This will
//                       repeat the same calculation many times, attempting to find anomalous
//                       results and reporting them.
//
// Arguments:
//   gpu:          overarching type for storing parameters, coordinates, and the energy
//                 function
//   iter:         the iteration counter
//   cchk:         the current checking stage to store and analyze
//   nchk:         the number of different stage checks to perform
//---------------------------------------------------------------------------------------------
void kCheckGpuConsistency(gpuContext gpu, int iter, int cchk, int nchk)
{
  int i;
  static dmat lcRefFrc, lcNewFrc, nbRefFrc, nbNewFrc;
  PMEAccumulator *pFrc, *pnbFrc;

  // Always download the imaging atom array
  if (gpu->bNeighborList && gpu->pbImageIndex != NULL) {
    gpu->pbImageIndex->Download();
  }
  unsigned int *pImageAtomLookup = &(gpu->pbImageIndex->_pSysData[gpu->sim.imageStride * 2]);

  // Download the forces
  gpu->pbForceAccumulator->Download();
  pFrc = gpu->pbForceAccumulator->_pSysData;
  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
    pnbFrc = gpu->pbForceAccumulator->_pSysData + gpu->sim.stride3;
  }

  // Setup if this is the first iteration
  if (iter == 0) {
    if (cchk == 0) {
      lcRefFrc = CreateDmat(nchk, 3 * gpu->sim.atoms);
      lcNewFrc = CreateDmat(nchk, 3 * gpu->sim.atoms);
      nbRefFrc = CreateDmat(nchk, 3 * gpu->sim.atoms);
      nbNewFrc = CreateDmat(nchk, 3 * gpu->sim.atoms);
    }
    double *dtmp = lcRefFrc.map[cchk];
    double *dtm2p = nbRefFrc.map[cchk];
    for (i = 0; i < gpu->sim.atoms; i++) {
      int i1 = pImageAtomLookup[i];
      dtmp[3*i    ] = (double)pFrc[i1                   ] * (double)ONEOVERFORCESCALE;
      dtmp[3*i + 1] = (double)pFrc[i1 + gpu->sim.stride ] * (double)ONEOVERFORCESCALE;
      dtmp[3*i + 2] = (double)pFrc[i1 + gpu->sim.stride2] * (double)ONEOVERFORCESCALE;
      if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
        dtm2p[3*i    ] = (double)pnbFrc[i1                   ] * (double)ONEOVERFORCESCALE;
        dtm2p[3*i + 1] = (double)pnbFrc[i1 + gpu->sim.stride ] * (double)ONEOVERFORCESCALE;
        dtm2p[3*i + 2] = (double)pnbFrc[i1 + gpu->sim.stride2] * (double)ONEOVERFORCESCALE;
      }
    }
  }
  else {
    double *dtmp = lcNewFrc.map[cchk];
    double *dtm2p = nbNewFrc.map[cchk];
    for (i = 0; i < gpu->sim.atoms; i++) {
      int i1 = pImageAtomLookup[i];
      dtmp[3*i    ] = (double)pFrc[i1                   ] * (double)ONEOVERFORCESCALE;
      dtmp[3*i + 1] = (double)pFrc[i1 + gpu->sim.stride ] * (double)ONEOVERFORCESCALE;
      dtmp[3*i + 2] = (double)pFrc[i1 + gpu->sim.stride2] * (double)ONEOVERFORCESCALE;
      if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
	dtm2p[3*i    ] = (double)pnbFrc[i1                   ] * (double)ONEOVERFORCESCALE;
	dtm2p[3*i + 1] = (double)pnbFrc[i1 + gpu->sim.stride ] * (double)ONEOVERFORCESCALE;
	dtm2p[3*i + 2] = (double)pnbFrc[i1 + gpu->sim.stride2] * (double)ONEOVERFORCESCALE;
      }
    }
  }

  // Compare the reference coordinates to the new coordinates
  if (iter > 0) {
    for (i = 0; i < gpu->sim.atoms; i++) {
      double rdx, rdy, rdz, ndx, ndy, ndz;
      rdx = lcRefFrc.map[cchk][3*i    ];
      rdy = lcRefFrc.map[cchk][3*i + 1];
      rdz = lcRefFrc.map[cchk][3*i + 2];
      ndx = lcNewFrc.map[cchk][3*i    ];
      ndy = lcNewFrc.map[cchk][3*i + 1];
      ndz = lcNewFrc.map[cchk][3*i + 2];
      if (cchk > 0) {
        rdx -= lcRefFrc.map[cchk - 1][3*i    ];
        rdy -= lcRefFrc.map[cchk - 1][3*i + 1];
        rdz -= lcRefFrc.map[cchk - 1][3*i + 2];
        ndx -= lcNewFrc.map[cchk - 1][3*i    ];
        ndy -= lcNewFrc.map[cchk - 1][3*i + 1];
        ndz -= lcNewFrc.map[cchk - 1][3*i + 2];
      }
      if (fabs(rdx - ndx) > 1.0e-6 || fabs(rdy - ndy) > 1.0e-6 || fabs(rdz - ndz) > 1.0e-6) {
	printf("Consistency check: Iter %4d check %2d atom %7d local forces do not match.\n",
	       iter, cchk, i);
	printf("   %9.4f %9.4f %9.4f != %9.4f %9.4f %9.4f\n", rdx, rdy, rdz, ndx, ndy, ndz);
      }
      rdx = nbRefFrc.map[cchk][3*i    ];
      rdy = nbRefFrc.map[cchk][3*i + 1];
      rdz = nbRefFrc.map[cchk][3*i + 2];
      ndx = nbNewFrc.map[cchk][3*i    ];
      ndy = nbNewFrc.map[cchk][3*i + 1];
      ndz = nbNewFrc.map[cchk][3*i + 2];
      if (cchk > 0) {
        rdx -= nbRefFrc.map[cchk - 1][3*i    ];
        rdy -= nbRefFrc.map[cchk - 1][3*i + 1];
        rdz -= nbRefFrc.map[cchk - 1][3*i + 2];
        ndx -= nbNewFrc.map[cchk - 1][3*i    ];
        ndy -= nbNewFrc.map[cchk - 1][3*i + 1];
        ndz -= nbNewFrc.map[cchk - 1][3*i + 2];
      }
      if (fabs(rdx - ndx) > 1.0e-6 || fabs(rdy - ndy) > 1.0e-6 || fabs(rdz - ndz) > 1.0e-6) {
	printf("Consistency check: Iter %4d check %2d atom %7d non-bonded forces do not "
               "match.\n", iter, cchk, i);
	printf("   %9.4f %9.4f %9.4f != %9.4f %9.4f %9.4f\n", rdx, rdy, rdz, ndx, ndy, ndz);
      }
    }
  }
}

__global__ void
__launch_bounds__(THREADS_PER_BLOCK, 1)
kCalculateSGLDaverages_kernel(double dt)
{
  extern __shared__ volatile SGLDAverage sG[];

  PMEFloat ekt   = (PMEFloat)0.0;
  PMEFloat eksg  = (PMEFloat)0.0;
  PMEFloat sumgam = (PMEFloat)0.0;
  PMEFloat x0sg = (PMEFloat)0.0;
  PMEFloat y0sg = (PMEFloat)0.0;
  PMEFloat z0sg = (PMEFloat)0.0;
  PMEFloat x1sg = (PMEFloat)0.0;
  PMEFloat y1sg = (PMEFloat)0.0;
  PMEFloat z1sg = (PMEFloat)0.0;
  PMEFloat x2sg = (PMEFloat)0.0;
  PMEFloat y2sg = (PMEFloat)0.0;
  PMEFloat z2sg = (PMEFloat)0.0;
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;

  // Sum up properties
  while (pos < cSim.atoms) {
    int wsg = cSim.pWsg[pos];
    if(wsg>0){
    PMEFloat mass =  cSim.pAtomMass[pos];
    PMEFloat vx   = cSim.pVelX[pos];
    PMEFloat lvx  = cSim.pLVelX[pos];
    PMEFloat vy   = cSim.pVelY[pos];
    PMEFloat lvy  = cSim.pLVelY[pos];
    PMEFloat vz   = cSim.pVelZ[pos];
    PMEFloat lvz  = cSim.pLVelZ[pos];
    PMEFloat svx  = vx + lvx;
    PMEFloat svy  = vy + lvy;
    PMEFloat svz  = vz + lvz;
    PMEFloat fpsg   = cSim.pFPsg[pos];
    PMEFloat ppsg   = cSim.pPPsg[pos];
    PMEFloat x0   = cSim.pX0sg[pos];
    PMEFloat y0   = cSim.pY0sg[pos];
    PMEFloat z0   = cSim.pZ0sg[pos];
    PMEFloat x1   = cSim.pX1sg[pos];
    PMEFloat y1   = cSim.pY1sg[pos];
    PMEFloat z1   = cSim.pZ1sg[pos];
    PMEFloat x2   = cSim.pX2sg[pos];
    PMEFloat y2   = cSim.pY2sg[pos];
    PMEFloat z2   = cSim.pZ2sg[pos];
    PMEFloat tsgfac   = cSim.sgavg/dt/20.455;
    PMEFloat vsgx   = tsgfac * (x0 -x1 );
    PMEFloat vsgy   = tsgfac * (y0 -y1 );
    PMEFloat vsgz   = tsgfac * (z0 -z1 );
    ekt          += mass * (svx*svx + svy*svy + svz*svz);
    eksg         += mass * ( vsgx*vsgx +  vsgy*vsgy +  vsgz*vsgz);
    sumgam       += - fpsg/ppsg;
    x0sg         += mass * x0;
    y0sg         += mass * y0;
    z0sg         += mass * z0;
    x1sg         += mass * x1;
    y1sg         += mass * y1;
    z1sg         += mass * z1;
    x2sg         += mass * x2;
    y2sg         += mass * y2;
    z2sg         += mass * z2;
    }
    pos          += blockDim.x * gridDim.x;      
    
  }
  ekt   *= (PMEFloat)0.125 ;
  eksg  *= (PMEFloat)0.5;
  
  sG[threadIdx.x].SGLDavg.EKT   = ekt;
  sG[threadIdx.x].SGLDavg.EKSG  = eksg;
  sG[threadIdx.x].SGLDavg.GAMSG = sumgam;
  sG[threadIdx.x].SGLDavg.COM0SG[0] = x0sg;
  sG[threadIdx.x].SGLDavg.COM0SG[1] = y0sg;
  sG[threadIdx.x].SGLDavg.COM0SG[2] = z0sg;
  sG[threadIdx.x].SGLDavg.COM1SG[0] = x1sg;
  sG[threadIdx.x].SGLDavg.COM1SG[1] = y1sg;
  sG[threadIdx.x].SGLDavg.COM1SG[2] = z1sg;
  sG[threadIdx.x].SGLDavg.COM2SG[0] = x2sg;
  sG[threadIdx.x].SGLDavg.COM2SG[1] = y2sg;
  sG[threadIdx.x].SGLDavg.COM2SG[2] = z2sg;

  // Reduce per-thread SGLD averages
  __syncthreads();
  unsigned int m = 1;
  while (m < blockDim.x) {
    int p = threadIdx.x + m;
    if(p < blockDim.x){
      ekt   = sG[p].SGLDavg.EKT;
      eksg  = sG[p].SGLDavg.EKSG ;
      sumgam =  sG[p].SGLDavg.GAMSG ;
      x0sg =  sG[p].SGLDavg.COM0SG[0] ;
      y0sg =  sG[p].SGLDavg.COM0SG[1] ;
      z0sg =  sG[p].SGLDavg.COM0SG[2] ;
      x1sg =  sG[p].SGLDavg.COM1SG[0] ;
      y1sg =  sG[p].SGLDavg.COM1SG[1] ;
      z1sg =  sG[p].SGLDavg.COM1SG[2] ;
      x2sg =  sG[p].SGLDavg.COM2SG[0] ;
      y2sg =  sG[p].SGLDavg.COM2SG[1] ;
      z2sg =  sG[p].SGLDavg.COM2SG[2] ;
    }
    else{
      ekt =(PMEFloat)0.0f ;
      eksg =(PMEFloat)0.0f ;
      sumgam =(PMEFloat)0.0f ;
      x0sg =(PMEFloat)0.0f ;
      y0sg =(PMEFloat)0.0f ;
      z0sg =(PMEFloat)0.0f ;
      x1sg =(PMEFloat)0.0f ;
      y1sg =(PMEFloat)0.0f ;
      z1sg =(PMEFloat)0.0f ;
      x2sg =(PMEFloat)0.0f ;
      y2sg =(PMEFloat)0.0f ;
      z2sg =(PMEFloat)0.0f ;
  }
    __syncthreads();
    sG[threadIdx.x].SGLDavg.EKT   += ekt;
    sG[threadIdx.x].SGLDavg.EKSG  += eksg;
    sG[threadIdx.x].SGLDavg.GAMSG += sumgam;
    sG[threadIdx.x].SGLDavg.COM0SG[0] += x0sg;
    sG[threadIdx.x].SGLDavg.COM0SG[1] += y0sg;
    sG[threadIdx.x].SGLDavg.COM0SG[2] += z0sg;
    sG[threadIdx.x].SGLDavg.COM1SG[0] += x1sg;
    sG[threadIdx.x].SGLDavg.COM1SG[1] += y1sg;
    sG[threadIdx.x].SGLDavg.COM1SG[2] += z1sg;
    sG[threadIdx.x].SGLDavg.COM2SG[0] += x2sg;
    sG[threadIdx.x].SGLDavg.COM2SG[1] += y2sg;
    sG[threadIdx.x].SGLDavg.COM2SG[2] += z2sg;
    __syncthreads();
    m *= 2;
  }

  // Save result
  if (threadIdx.x < 12) {
    cSim.pSGLDAverage[blockIdx.x].array[threadIdx.x] = sG[0].array[threadIdx.x];
  }
}

__global__ void
__launch_bounds__(THREADS_PER_BLOCK, 1)
kPMECalculateSGLDaverages_kernel(double dt)
{
  extern __shared__ volatile SGLDAverage sG[];

  PMEFloat ekt   = (PMEFloat)0.0;
  PMEFloat eksg  = (PMEFloat)0.0;
  PMEFloat sumgam = (PMEFloat)0.0;
  PMEFloat x0sg = (PMEFloat)0.0;
  PMEFloat y0sg = (PMEFloat)0.0;
  PMEFloat z0sg = (PMEFloat)0.0;
  PMEFloat x1sg = (PMEFloat)0.0;
  PMEFloat y1sg = (PMEFloat)0.0;
  PMEFloat z1sg = (PMEFloat)0.0;
  PMEFloat x2sg = (PMEFloat)0.0;
  PMEFloat y2sg = (PMEFloat)0.0;
  PMEFloat z2sg = (PMEFloat)0.0;
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;

  // Sum up properties
  while (pos < cSim.atoms) {
    int wsg = cSim.pImageWsg[pos];
    if(wsg>0){
    PMEFloat mass = cSim.pImageMass[pos];
    PMEFloat vx   = cSim.pImageVelX[pos];
    PMEFloat lvx  = cSim.pImageLVelX[pos];
    PMEFloat vy   = cSim.pImageVelY[pos];
    PMEFloat lvy  = cSim.pImageLVelY[pos];
    PMEFloat vz   = cSim.pImageVelZ[pos];
    PMEFloat lvz  = cSim.pImageLVelZ[pos];
    PMEFloat svx  = vx + lvx;
    PMEFloat svy  = vy + lvy;
    PMEFloat svz  = vz + lvz;
    PMEFloat fpsg   = cSim.pImageFPsg[pos];
    PMEFloat ppsg   = cSim.pImagePPsg[pos];
    PMEFloat x0   = cSim.pImageX0sg[pos];
    PMEFloat y0   = cSim.pImageY0sg[pos];
    PMEFloat z0   = cSim.pImageZ0sg[pos];
    PMEFloat x1   = cSim.pImageX1sg[pos];
    PMEFloat y1   = cSim.pImageY1sg[pos];
    PMEFloat z1   = cSim.pImageZ1sg[pos];
    PMEFloat x2   = cSim.pImageX2sg[pos];
    PMEFloat y2   = cSim.pImageY2sg[pos];
    PMEFloat z2   = cSim.pImageZ2sg[pos];
    PMEFloat tsgfac   = cSim.sgavg/dt/20.455;
    PMEFloat vsgx   = tsgfac * (x0-x1);
    PMEFloat vsgy   = tsgfac * (y0-y1);
    PMEFloat vsgz   = tsgfac * (z0-z1);
    ekt          += mass * (svx*svx + svy*svy + svz*svz);
    eksg         += mass * ( vsgx*vsgx +  vsgy*vsgy +  vsgz*vsgz);
    sumgam       += - fpsg/ppsg;
    x0sg         += mass * x0;
    y0sg         += mass * y0;
    z0sg         += mass * z0;
    x1sg         += mass * x1;
    y1sg         += mass * y1;
    z1sg         += mass * z1;
    x2sg         += mass * x2;
    y2sg         += mass * y2;
    z2sg         += mass * z2;
    }
    pos          += blockDim.x * gridDim.x;
        
}
  ekt   *= (PMEFloat)0.125 ;
  eksg  *= (PMEFloat)0.5;
  
  sG[threadIdx.x].SGLDavg.EKT   = ekt;
  sG[threadIdx.x].SGLDavg.EKSG  = eksg;
  sG[threadIdx.x].SGLDavg.GAMSG = sumgam;
  sG[threadIdx.x].SGLDavg.COM0SG[0] = x0sg;
  sG[threadIdx.x].SGLDavg.COM0SG[1] = y0sg;
  sG[threadIdx.x].SGLDavg.COM0SG[2] = z0sg;
  sG[threadIdx.x].SGLDavg.COM1SG[0] = x1sg;
  sG[threadIdx.x].SGLDavg.COM1SG[1] = y1sg;
  sG[threadIdx.x].SGLDavg.COM1SG[2] = z1sg;
  sG[threadIdx.x].SGLDavg.COM2SG[0] = x2sg;
  sG[threadIdx.x].SGLDavg.COM2SG[1] = y2sg;
  sG[threadIdx.x].SGLDavg.COM2SG[2] = z2sg;

  // Reduce per-thread kinetic energies
  __syncthreads();
  unsigned int m = 1;
  while (m < blockDim.x) {
    int p = threadIdx.x + m;
    if(p < blockDim.x){
      ekt   = sG[p].SGLDavg.EKT;
      eksg  = sG[p].SGLDavg.EKSG ;
      sumgam =  sG[p].SGLDavg.GAMSG ;
      x0sg =  sG[p].SGLDavg.COM0SG[0] ;
      y0sg =  sG[p].SGLDavg.COM0SG[1] ;
      z0sg =  sG[p].SGLDavg.COM0SG[2] ;
      x1sg =  sG[p].SGLDavg.COM1SG[0] ;
      y1sg =  sG[p].SGLDavg.COM1SG[1] ;
      z1sg =  sG[p].SGLDavg.COM1SG[2] ;
      x2sg =  sG[p].SGLDavg.COM2SG[0] ;
      y2sg =  sG[p].SGLDavg.COM2SG[1] ;
      z2sg =  sG[p].SGLDavg.COM2SG[2] ;
    }
    else{
      ekt =(PMEFloat)0.0f ;
      eksg =(PMEFloat)0.0f ;
      sumgam =(PMEFloat)0.0f ;
      x0sg =(PMEFloat)0.0f ;
      y0sg =(PMEFloat)0.0f ;
      z0sg =(PMEFloat)0.0f ;
      x1sg =(PMEFloat)0.0f ;
      y1sg =(PMEFloat)0.0f ;
      z1sg =(PMEFloat)0.0f ;
      x2sg =(PMEFloat)0.0f ;
      y2sg =(PMEFloat)0.0f ;
      z2sg =(PMEFloat)0.0f ;
  }
    __syncthreads();
    sG[threadIdx.x].SGLDavg.EKT   += ekt;
    sG[threadIdx.x].SGLDavg.EKSG  += eksg;
    sG[threadIdx.x].SGLDavg.GAMSG += sumgam;
    sG[threadIdx.x].SGLDavg.COM0SG[0] += x0sg;
    sG[threadIdx.x].SGLDavg.COM0SG[1] += y0sg;
    sG[threadIdx.x].SGLDavg.COM0SG[2] += z0sg;
    sG[threadIdx.x].SGLDavg.COM1SG[0] += x1sg;
    sG[threadIdx.x].SGLDavg.COM1SG[1] += y1sg;
    sG[threadIdx.x].SGLDavg.COM1SG[2] += z1sg;
    sG[threadIdx.x].SGLDavg.COM2SG[0] += x2sg;
    sG[threadIdx.x].SGLDavg.COM2SG[1] += y2sg;
    sG[threadIdx.x].SGLDavg.COM2SG[2] += z2sg;
    __syncthreads();
    m *= 2;
  }

  // Save result
  if (threadIdx.x < 12) {
    cSim.pSGLDAverage[blockIdx.x].array[threadIdx.x] = sG[0].array[threadIdx.x];
  }
}

void kCalculateSGLDaverages(gpuContext gpu, double dt)
{
        if (gpu->bNeighborList)
        {
#ifdef use_DPFP
    // Use half as many threads: the dynamically allocated __shared__ array would
    // overload the hardware limits otherwise.
          kPMECalculateSGLDaverages_kernel<<<gpu->blocks, gpu->threadsPerBlock/2, gpu->threadsPerBlock/2 * sizeof(SGLDAverage)>>>(dt);
#else
          kPMECalculateSGLDaverages_kernel<<<gpu->blocks, gpu->threadsPerBlock, gpu->threadsPerBlock * sizeof(SGLDAverage)>>>(dt);
#endif
        }
        else
        {
#ifdef use_DPFP
    // Use half as many threads: the dynamically allocated __shared__ array would
    // overload the hardware limits otherwise.
          kCalculateSGLDaverages_kernel<<<gpu->blocks, gpu->threadsPerBlock/2, gpu->threadsPerBlock/2 * sizeof(SGLDAverage)>>>(dt);
#else
          kCalculateSGLDaverages_kernel<<<gpu->blocks, gpu->threadsPerBlock, gpu->threadsPerBlock * sizeof(SGLDAverage)>>>(dt);
#endif
        }
    LAUNCHERROR("kCalculateSGLDaverages");
}
#ifdef MPI

//---------------------------------------------------------------------------------------------
// kCopyToAccumulator_kernel: kernel for accumulator copying.  This appears to be for debugging
//                            purposes only.
//---------------------------------------------------------------------------------------------
__global__ void
#if (__CUDA_ARCH__ != 750)
__LAUNCH_BOUNDS__(512, 3)
#else
__LAUNCH_BOUNDS__(512, 2)
#endif
kCopyToAccumulator_kernel(PMEAccumulator* p1, PMEAccumulator* p2, int size)
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  if (pos < size) {
    p1[pos] = p2[pos];
  }
}

//---------------------------------------------------------------------------------------------
// kCopyToAccumulator: copy the contents of one accumulator vector to another.  Think of this
//                     as strncpy() but on a GPU with strings of numbers.
//
// Arguments:
//   gpu:          overarching type for storing parameters, coordinates, and the energy
//                 function
//   p1:           the accumulator that will get the contents of p2
//   p2:           source of the information to copy
//   size:         length of the vectors p1 and p2
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
void kCopyToAccumulator(gpuContext gpu, PMEAccumulator* p1, PMEAccumulator* p2, int size)
{
  unsigned int blocks = (size + 511) / 512;
  kCopyToAccumulator_kernel<<<blocks, 512>>>(p1, p2, size);
  LAUNCHERROR("kCopyToAccumulator");
}

//---------------------------------------------------------------------------------------------
// kAddAccumulators_kernel: add two accumulators together.  This is used in reducing the force
//                          accumulators coming out of parallel GPU runs.
//
// Arguments:
//   p1:    accumulator that will ultimately hold the sum
//   p2:    contribution to add into p1
//   size:  size of p1 and p2
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kAddAccumulators_kernel(PMEAccumulator* p1, PMEAccumulator* p2, int size)
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  while (pos < size) {
    p1[pos] += p2[pos];
    pos += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kAddAccumulators: host function to launch accumulator addition.
//
// Arguments:
//   gpu:   overarching type for storing parameters, coordinates, and the energy function
//   p1:    accumulator that will ultimately hold the sum
//   p2:    contribution to add into p1
//   size:  size of p1 and p2
//---------------------------------------------------------------------------------------------
void kAddAccumulators(gpuContext gpu, PMEAccumulator* p1, PMEAccumulator* p2, int size)
{
  kAddAccumulators_kernel<<<gpu->blocks, gpu->generalThreadsPerBlock>>>(p1, p2, size);
  LAUNCHERROR("kAddAccumulators");
}

#  ifdef AMBER_PLATFORM_AMD

//---------------------------------------------------------------------------------------------
// kPersistentAddAccumulators_kernel: add nGpus accumulators together.  This is used in
//                                    reducing the force accumulators coming out of parallel
//                                    GPU runs.
//
// Arguments:
//   buffer:           buffer on the current GPU, it will contain the per-element sum of
//                     accumulators of all GPUs
//   accumulator:      temporary buffer visible on the current GPU visible for all other GPUs
//   peerAccumulators: temporary buffer (accumulator) on all GPUs
//   stages:           host buffer where device writes the current finished stage and host
//                     wait for it and then resets for further processing
//   nGpus:            number of GPUs
//   gpuID:            current GPU
//   start:            starting index of a segment for processing by the current GPU
//   end:              ending index of a segment for processing by the current GPU
//   size:             total size of buffer
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kPersistentAddAccumulators_kernel(PMEAccumulator* buffer,
                                  PMEAccumulator* accumulator,
                                  PMEAccumulator** peerAccumulators,
                                  unsigned int* stages,
                                  int nGpus, int gpuID, int start, int end, int size)
{
  cooperative_groups::grid_group grid = cooperative_groups::this_grid();

  const int stride = blockDim.x * gridDim.x;
  const int tid = threadIdx.x + blockIdx.x * blockDim.x;

  __shared__ PMEAccumulator * shAccumulators[16];

  if (threadIdx.x == 0) {
    for (int g = 0; g < nGpus; g++) {
      shAccumulators[g] = gpuID == g ? buffer : peerAccumulators[g];
    }
  }
  __syncthreads();

  for (int i = tid; i < start; i += stride) {
    accumulator[i] = buffer[i];
  }
  for (int i = end + tid; i < size; i += stride) {
    accumulator[i] = buffer[i];
  }
  __threadfence();

  grid.sync();
  if (tid == 0) {
    atomicExch(stages, 1);
    while (atomicMax(stages, 0) != 0) {}
  }
  grid.sync();

  for (int i = start + tid; i < end; i += stride){
    PMEAccumulator sum = 0;
    for (int g = 0; g < nGpus; g++) {
      sum += shAccumulators[g][i];
    }
    for (int g = 0; g < nGpus; g++) {
      shAccumulators[g][i] = sum;
    }
  }

  grid.sync();
  if (tid == 0) {
    atomicExch(stages, 2);
    while (atomicMax(stages, 0) != 0) {}
  }
  grid.sync();

  for (int i = tid; i < start; i += stride) {
    buffer[i] = accumulator[i];
  }
  for (int i = end + tid; i < size; i += stride) {
    buffer[i] = accumulator[i];
  }
}

static void kAddAccumulatorsWaitForStage(gpuContext gpu, unsigned int stage) 
{
  while (true) {
    unsigned int result = __atomic_load_n(&gpu->pbPersistentStages->_pSysData[0],
                                          __ATOMIC_RELAXED);
    if (result == stage) {
      MPI_Barrier(gpu->comm);
      __atomic_store_n(&gpu->pbPersistentStages->_pSysData[0], 0, __ATOMIC_RELAXED);
      break;
    }
  }
}

//---------------------------------------------------------------------------------------------
// kPersistentAddAccumulators: host function to launch accumulator addition among all GPUs.
//
// Arguments:
//   gpu:   overarching type for storing parameters, coordinates, and the energy function
//   pBuff: buffer on the current GPU, it will contain the per-element sum of accumulators
//          of all GPUs
//   size:  the number of element in pBuff
//---------------------------------------------------------------------------------------------
void kPersistentAddAccumulators(gpuContext gpu, PMEAccumulator* pBuff, int size)
{
  cudaError_t status;

  int blocks = 64;
  int threads = 64;
  int maxBlocksPerSM;
  status = cudaOccupancyMaxActiveBlocksPerMultiprocessor(&maxBlocksPerSM,
                                                         kPersistentAddAccumulators_kernel,
                                                         threads, 0);
  RTERROR(status, "cudaOccupancyMaxActiveBlocksPerMultiprocessor failed");
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, gpu->gpu_device_id);
  blocks = std::min(blocks, deviceProp.multiProcessorCount * maxBlocksPerSM);

  // Each GPU processes adds its own segment of accumulators
  int chunk = (size + gpu->nGpus - 1) / gpu->nGpus;
  int start = chunk * gpu->gpuID;
  int end = std::min(chunk * (gpu->gpuID + 1), size);
  void* params[] = {
    (void*)&pBuff,
    (void*)&gpu->pbPeerAccumulator->_pDevData,
    (void*)&gpu->pbPeerAccumulatorList->_pDevData,
    (void*)&gpu->pbPersistentStages->_pDevData,
    (void*)&gpu->nGpus,
    (void*)&gpu->gpuID,
    (void*)&start,
    (void*)&end,
    (void*)&size,
  };
  status = cudaLaunchCooperativeKernel(kPersistentAddAccumulators_kernel,
                                       dim3(blocks), dim3(threads), params, 0, 0);
  RTERROR(status, "cudaLaunchCooperativeKernel failed on kAddAccumulators2_kernel");

  // Barriers are implemented using atomic opertations on host memory and MPI_Barrier
  kAddAccumulatorsWaitForStage(gpu, 1);
  kAddAccumulatorsWaitForStage(gpu, 2);
}

#  endif // AMBER_PLATFORM_AMD

//---------------------------------------------------------------------------------------------
// kReduceForces_kernel: fold the non-bonded forces into the force accumulators containing
//                       bonded interactions.  This is called only in the context of MPI
//                       multi-GPU calculations.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(UPDATE_THREADS_PER_BLOCK, 1)
kReduceForces_kernel()
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  if (pos < cSim.stride3) {
    cSim.pForceAccumulator[pos] += cSim.pNBForceAccumulator[pos];

    // Zero out the non-bonded force values to avoid double-counting them at any time later.
    cSim.pNBForceAccumulator[pos] = 0.0;
  }
}

//---------------------------------------------------------------------------------------------
// kReduceForces: host function for launching force reduction in the context of parallel GPU
//                calculations.
//
// Arguments:
//   gpu:   overarching type for storing parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kReduceForces(gpuContext gpu)
{
  kReduceForces_kernel<<<3 * gpu->updateBlocks, gpu->updateThreadsPerBlock>>>();
  LAUNCHERROR("kReduceForces");
}

#endif // MPI
