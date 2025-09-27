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
// SetkCalculateGBNonbondEnergy1Sim: this is called by gpuCopyConstants (see gpu.cpp) and is
//                                   jused to port GB constants to the device.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void SetkCalculateGBNonbondEnergy1Sim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(simulationConst));
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

//---------------------------------------------------------------------------------------------
// SetkCalculateGBNonbondEnergy1Sim: this will download critical GB constants from the device.
//                                   It is not currently called anywhere in the code.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
void GetkCalculateGBNonBondEnergy1Sim(gpuContext gpu)
{
    cudaError_t status;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(simulationConst));
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

#endif

//---------------------------------------------------------------------------------------------
// Kernels for calculating forces and energies for various flavors of GB, even "gas phase"
//---------------------------------------------------------------------------------------------
// Forces for regular GB
__global__ void
__LAUNCH_BOUNDS__(GBNONBONDENERGY1_THREADS_PER_BLOCK, GBNONBONDENERGY1_BLOCKS_LOAD)
kCalculateGBNonbondForces1_kernel()
#include "kCalculateGBNonbondEnergy1.h"

//Forces for pwsasa + GB
#define GB_GBSA3
__global__ void
__LAUNCH_BOUNDS__(GBNONBONDENERGY1_THREADS_PER_BLOCK, GBNONBONDENERGY1_BLOCKS_LOAD)
kCalculateGBSA3Forces_kernel()
#include "kCalculateGBNonbondEnergy1.h"
#undef GB_GBSA3


//---------------------------------------------------------------------------------------------
// Gas phase GB
//---------------------------------------------------------------------------------------------
// Forces for gas phase GB
#define GB_IGB6
__global__ void
  __LAUNCH_BOUNDS__(GBNONBONDENERGY1_THREADS_PER_BLOCK, GBNONBONDENERGY1_BLOCKS_LOAD)
kCalculateGBNonbondForces1IGB6_kernel()
#include "kCalculateGBNonbondEnergy1.h"
#undef GB_IGB6

//For Energy
//---------------------------------------------------------------------------------------------
#define GB_ENERGY
__global__
__LAUNCH_BOUNDS__(GBNONBONDENERGY1_THREADS_PER_BLOCK, GBNONBONDENERGY1_BLOCKS_LOAD)
void kCalculateGBNonbondEnergy1_kernel()
#include "kCalculateGBNonbondEnergy1.h"
//---------------------------------------------------------------------------------------------
// pwsasa
#define GB_GBSA3
__global__
__LAUNCH_BOUNDS__(GBNONBONDENERGY1_THREADS_PER_BLOCK, GBNONBONDENERGY1_BLOCKS_LOAD)
void kCalculateGBSA3Energy_kernel()
#include "kCalculateGBNonbondEnergy1.h"
#undef GB_GBSA3

#undef GB_ENERGY

//---------------------------------------------------------------------------------------------
// Gas phase GB
//---------------------------------------------------------------------------------------------
#define GB_IGB6
#define GB_ENERGY
__global__
__LAUNCH_BOUNDS__(GBNONBONDENERGY1_THREADS_PER_BLOCK, GBNONBONDENERGY1_BLOCKS_LOAD)
void kCalculateGBNonbondEnergy1IGB6_kernel()
#include "kCalculateGBNonbondEnergy1.h"
#undef GB_ENERGY
#undef GB_IGB6

//---------------------------------------------------------------------------------------------
#define GB_MINIMIZATION
__global__ void
__LAUNCH_BOUNDS__(GBNONBONDENERGY1_THREADS_PER_BLOCK, GBNONBONDENERGY1_BLOCKS_LOAD)
kCalculateGBNonbondMinimizationForces1_kernel()
#include "kCalculateGBNonbondEnergy1.h"

//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// Gas phase GB
//---------------------------------------------------------------------------------------------
#define GB_IGB6
__global__ void
__LAUNCH_BOUNDS__(GBNONBONDENERGY1_THREADS_PER_BLOCK, GBNONBONDENERGY1_BLOCKS_LOAD)
kCalculateGBNonbondMinimizationForces1IGB6_kernel()
#include "kCalculateGBNonbondEnergy1.h"

//---------------------------------------------------------------------------------------------
#define GB_ENERGY
__global__
__LAUNCH_BOUNDS__(GBNONBONDENERGY1_THREADS_PER_BLOCK, GBNONBONDENERGY1_BLOCKS_LOAD)
void kCalculateGBNonbondMinimizationEnergy1IGB6_kernel()
#include "kCalculateGBNonbondEnergy1.h"
#undef GB_IGB6

//---------------------------------------------------------------------------------------------
__global__
__LAUNCH_BOUNDS__(GBNONBONDENERGY1_THREADS_PER_BLOCK, GBNONBONDENERGY1_BLOCKS_LOAD)
void kCalculateGBNonbondMinimizationEnergy1_kernel()
#include "kCalculateGBNonbondEnergy1.h"
#undef GB_ENERGY
#undef GB_MINIMIZATION


//---------------------------------------------------------------------------------------------
// kCalculateGBNonbondEnergy1InitKernels: what the name says.  Called by gpu_init_ in gpu.cpp.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kCalculateGBNonbondEnergy1InitKernels(gpuContext gpu)
{
  cudaFuncSetSharedMemConfig(kCalculateGBNonbondEnergy1_kernel,
                             cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalculateGBNonbondForces1_kernel,
                             cudaSharedMemBankSizeEightByte);
  //pwsasa set shared memory bank width to be eight bytes natively when launching the kernel
  cudaFuncSetSharedMemConfig(kCalculateGBSA3Energy_kernel,
                             cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalculateGBSA3Forces_kernel,
                             cudaSharedMemBankSizeEightByte);

  cudaFuncSetSharedMemConfig(kCalculateGBNonbondEnergy1IGB6_kernel,
                             cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalculateGBNonbondForces1IGB6_kernel,
                             cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalculateGBNonbondMinimizationEnergy1_kernel,
                             cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalculateGBNonbondMinimizationForces1_kernel,
                             cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalculateGBNonbondMinimizationEnergy1IGB6_kernel,
                             cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalculateGBNonbondMinimizationForces1IGB6_kernel,
                             cudaSharedMemBankSizeEightByte);

  if (gpu->sm_version >= SM_3X) {
    cudaFuncSetCacheConfig(kCalculateGBNonbondEnergy1_kernel, cudaFuncCachePreferEqual);
    cudaFuncSetCacheConfig(kCalculateGBNonbondForces1_kernel, cudaFuncCachePreferEqual);
    //pwsasa
    cudaFuncSetCacheConfig(kCalculateGBSA3Energy_kernel, cudaFuncCachePreferEqual);
    cudaFuncSetCacheConfig(kCalculateGBSA3Forces_kernel, cudaFuncCachePreferEqual);

    cudaFuncSetCacheConfig(kCalculateGBNonbondEnergy1IGB6_kernel, cudaFuncCachePreferEqual);
    cudaFuncSetCacheConfig(kCalculateGBNonbondForces1IGB6_kernel, cudaFuncCachePreferEqual);
    cudaFuncSetCacheConfig(kCalculateGBNonbondMinimizationEnergy1_kernel,
                           cudaFuncCachePreferEqual);
    cudaFuncSetCacheConfig(kCalculateGBNonbondMinimizationForces1_kernel,
                           cudaFuncCachePreferEqual);
    cudaFuncSetCacheConfig(kCalculateGBNonbondMinimizationEnergy1IGB6_kernel,
                           cudaFuncCachePreferEqual);
    cudaFuncSetCacheConfig(kCalculateGBNonbondMinimizationForces1IGB6_kernel,
                           cudaFuncCachePreferEqual);
  }
}

//---------------------------------------------------------------------------------------------
// kCalculateGBNonbondForces1: launch the appropriate kernel for computing GB forces.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kCalculateGBNonbondForces1(gpuContext gpu)
{
  if (gpu->imin == 0) {
    if (gpu->sim.igb != 6) {
       if(gpu->gbsa == 3){ //pwsasa
         kCalculateGBSA3Forces_kernel<<<gpu->GBNonbondEnergy1Blocks,
                                        gpu->GBNonbondEnergy1ThreadsPerBlock>>>();
       }
       else{
          kCalculateGBNonbondForces1_kernel<<<gpu->GBNonbondEnergy1Blocks,
                                          gpu->GBNonbondEnergy1ThreadsPerBlock>>>();
       }
    }
    else {
      kCalculateGBNonbondForces1IGB6_kernel<<<gpu->GBNonbondEnergy1Blocks,
                                              gpu->GBNonbondEnergy1ThreadsPerBlock>>>();
    }
  }
  else {
    if (gpu->sim.igb != 6) {
      kCalculateGBNonbondMinimizationForces1_kernel<<<gpu->GBNonbondEnergy1Blocks,
                                                      gpu->GBNonbondEnergy1ThreadsPerBlock>>>();
    }
    else {
      kCalculateGBNonbondMinimizationForces1IGB6_kernel<<<gpu->GBNonbondEnergy1Blocks,
                                                          gpu->GBNonbondEnergy1ThreadsPerBlock>>>();
    }
  }
  LAUNCHERROR("kCalculateGBNonbondForces1");
}

//---------------------------------------------------------------------------------------------
// kCalculateGBNonbondForces1: launch the appropriate kernel for computing GB energies.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kCalculateGBNonbondEnergy1(gpuContext gpu)
{
  if (gpu->imin == 0) {
    if (gpu->sim.igb != 6) {
       if (gpu->gbsa == 3) { //pwsasa
          kCalculateGBSA3Energy_kernel<<<gpu->GBNonbondEnergy1Blocks,
                                          gpu->GBNonbondEnergy1ThreadsPerBlock>>>();
       }
       else{
          kCalculateGBNonbondEnergy1_kernel<<<gpu->GBNonbondEnergy1Blocks,
                                          gpu->GBNonbondEnergy1ThreadsPerBlock>>>();
       }
    }
    else {
      kCalculateGBNonbondEnergy1IGB6_kernel<<<gpu->GBNonbondEnergy1Blocks,
                                              gpu->GBNonbondEnergy1ThreadsPerBlock>>>();
    }
  }
  else {
    if (gpu->sim.igb != 6) {
      kCalculateGBNonbondMinimizationEnergy1_kernel<<<gpu->GBNonbondEnergy1Blocks,
                                                      gpu->GBNonbondEnergy1ThreadsPerBlock>>>();
    }
    else {
      kCalculateGBNonbondMinimizationEnergy1IGB6_kernel<<<gpu->GBNonbondEnergy1Blocks,
                                                          gpu->GBNonbondEnergy1ThreadsPerBlock>>>();
    }
  }

  LAUNCHERROR("kCalculateGBNonbondEnergy1");
}

//--------------------------------------------------------------------------------------------
// kReduceMaxsasaEsurf_kernel: adding a kernel for reductions of atom loop to add
//                             maxsasa to Esurf from kCalculateGBNonbondEnergy1 kernel
// -------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(REDUCEBUFFER_THREADS_PER_BLOCK, 1)
kReduceMaxsasaEsurf_kernel()
{
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

  // Surface tension maxsasa correction--add the intecept only to the first thread
  PMEDouble tempes = (threadIdx.x == 0 && blockIdx.x == 0) * (PMEFloat)361.108307897 *
                     cSim.surften;
  while (pos < cSim.atoms) {
    tempes += (PMEFloat)0.681431329392 * cSim.pgbsa_maxsasa[pos] * cSim.surften;
    pos += blockDim.x * gridDim.x;
  }

  for (unsigned int stride = warpSize >> 1; stride > 0; stride >>=1) {
    tempes += __SHFL_DOWN(WARP_MASK, tempes, stride);
  }

  // Atomic add to global memory
  if ((threadIdx.x & GRID_BITS_MASK) == 0) {
    atomicAdd(cSim.pESurf, llitoulli(llrint(tempes * ENERGYSCALE)));
  }
}

//---------------------------------------------------------------------------------------------
// kReduceMaxsasaEsurf: launcher for the kernel above, to collect esurf data from all threads
//
// Arguments:
//   gpu:  overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kReduceMaxsasaEsurf(gpuContext gpu)
{
  kReduceMaxsasaEsurf_kernel<<<gpu->blocks, gpu->reduceBufferThreadsPerBlock>>>();
  LAUNCHERROR("kReduceMaxsasaEsurf");
}

//---------------------------------------------------------------------------------------------
// kReduceGBTemp7_kernel: GB forces reduction kernel
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(REDUCEBUFFER_THREADS_PER_BLOCK, 1)
kReduceGBTemp7_kernel()
{
  bool bIGB2578 = (cSim.igb == 2) || (cSim.igb == 5) || (cSim.igb == 7) || (cSim.igb == 8);
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  while (pos < cSim.atoms) {
    PMEFloat reff_i       = cSim.pReffSP[pos];
    PMEFloat psi_i        = cSim.pPsi[pos];
    PMEFloat rborn_i      = cSim.pAtomRBorn[pos];
    PMEFloat qi           = cSim.pAtomChargeSP[pos];
    PMEDouble sumdeijda_i = (PMEDouble)cSim.pSumdeijdaAccumulator[pos] * ONEOVERFORCESCALE;

    // Process Temp7 component
    PMEFloat expmkf = exp(-cSim.gb_kappa * reff_i) * cSim.extdiel_inv;
    PMEFloat dl     = cSim.intdiel_inv - expmkf;

    if (cSim.gbion == 2) {
      unsigned int ionmaski = cSim.pAtomIonMask[pos];
      PMEFloat intdiel_ion_2_inv = 1.0 / cSim.intdiel_ion_2;
      if (ionmaski>0){ // ion-ion pair (self)
        dl = cSim.intdiel_inv * intdiel_ion_2_inv - expmkf;
      }
      else { // solute-solute pair (self)
        // no change
      }
    }
    else if (cSim.gbion == 3) {
      unsigned int ionmaski = cSim.pAtomIonMask[pos];
      PMEFloat intdiel_ion_2_pp_inv = 1.0 / cSim.intdiel_ion_2_pp;
      PMEFloat intdiel_ion_2_nn_inv = 1.0 / cSim.intdiel_ion_2_nn;
      if (ionmaski == 3){ // anion-anion pair (self)
        dl = cSim.intdiel_inv * intdiel_ion_2_nn_inv - expmkf;
      }
      else if (ionmaski > 0){ // cation-cation pair (self)
        dl = cSim.intdiel_inv * intdiel_ion_2_pp_inv - expmkf;
      }
      else { // solute-solute pair (self)
        // no change
      }
    }

    PMEFloat qi2h   = (PMEFloat)0.50 * qi * qi;
    PMEFloat qid2h  = qi2h * dl;
    sumdeijda_i     = -sumdeijda_i + qid2h - cSim.gb_kappa * qi2h * expmkf * reff_i;
    if (cSim.alpb == 0) {

      // egb -= qid2h / reff_i;
    }
    else {

      // egb -= qid2h * (1.0 / reff_i + cSim.one_arad_beta);
      sumdeijda_i *= ((PMEFloat)1.0 + cSim.one_arad_beta * reff_i);
    }
    if (bIGB2578) {

      // New onufriev: scale values by alpha-, beta-, gamma- dependent factors later
      PMEFloat thi, thi2;
      if (cSim.igb == 8) {
        PMEFloat alpha = cSim.pgb_alpha[pos];
        PMEFloat gamma = cSim.pgb_gamma[pos];
        PMEFloat beta  = cSim.pgb_beta[pos];
        thi  = tanh((alpha + gamma * psi_i * psi_i - beta * psi_i) * psi_i);
        thi2 = (alpha + (PMEFloat)3.0 * gamma * psi_i * psi_i -
                (PMEFloat)2.0 * beta * psi_i) *
               ((PMEFloat)1.0 - thi * thi) * (rborn_i - cSim.offset) / rborn_i;
      }
      else {
        thi  = tanh((cSim.gb_alpha + cSim.gb_gamma * psi_i * psi_i -
                     cSim.gb_beta * psi_i) * psi_i);
        thi2 = (cSim.gb_alpha + (PMEFloat)3.0 * cSim.gb_gamma * psi_i * psi_i -
                (PMEFloat)2.0 * cSim.gb_beta * psi_i) *
               ((PMEFloat)1.0 - thi * thi) * (rborn_i - cSim.offset) / rborn_i;
      }
      sumdeijda_i *= thi2;
    }
    cSim.pTemp7[pos] = sumdeijda_i;
    pos += blockDim.x * gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kReduceGBTemp7: launch the eponymous kernel.  This is invoked in gpu_gb_ene_ (see gpu.cpp)
//                 on all GPUs but the master.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kReduceGBTemp7(gpuContext gpu)
{
  kReduceGBTemp7_kernel<<<gpu->blocks, gpu->reduceBufferThreadsPerBlock>>>();
  LAUNCHERROR("kReduceGBTemp7");
}

//---------------------------------------------------------------------------------------------
// kReduceGBTemp7Energy_kernel: GB energies reduction kernel
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(REDUCEBUFFER_THREADS_PER_BLOCK, 1)
kReduceGBTemp7Energy_kernel()
{
  //volatile __shared__ PMEDouble sE[1024];
  bool bIGB2578 = (cSim.igb == 2) || (cSim.igb == 5) || (cSim.igb == 7) || (cSim.igb == 8);
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  PMEDouble egb = (PMEDouble)0.0;
  while (pos < cSim.atoms) {
    PMEFloat reff_i  = cSim.pReffSP[pos];
    PMEFloat psi_i   = cSim.pPsi[pos];
    PMEFloat rborn_i = cSim.pAtomRBorn[pos];
    PMEFloat qi      = cSim.pAtomChargeSP[pos];
    PMEDouble sumdeijda_i = (PMEDouble)cSim.pSumdeijdaAccumulator[pos] * ONEOVERFORCESCALE;

    // Process Temp7 component
    PMEFloat expmkf = exp(-cSim.gb_kappa * reff_i) * cSim.extdiel_inv;
    PMEFloat dl     = cSim.intdiel_inv - expmkf;

    if (cSim.gbion == 2) {
      unsigned int ionmaski = cSim.pAtomIonMask[pos];
      PMEFloat intdiel_ion_2_inv = 1.0 / cSim.intdiel_ion_2;
      if (ionmaski>0){ // ion-ion pair (self)
        dl = cSim.intdiel_inv * intdiel_ion_2_inv - expmkf;
      }
      else { // solute-solute pair (self)
        // no change
      }
    }
    else if (cSim.gbion == 3) {
      unsigned int ionmaski = cSim.pAtomIonMask[pos];
      PMEFloat intdiel_ion_2_pp_inv = 1.0 / cSim.intdiel_ion_2_pp;
      PMEFloat intdiel_ion_2_nn_inv = 1.0 / cSim.intdiel_ion_2_nn;
      if (ionmaski == 3){ // anion-anion pair (self)
        dl = cSim.intdiel_inv * intdiel_ion_2_nn_inv - expmkf;
      }
      else if (ionmaski > 0){ // cation-cation pair (self)
        dl = cSim.intdiel_inv * intdiel_ion_2_pp_inv - expmkf;
      }
      else { // solute-solute pair (self)
        // no change
      }
    }

    PMEFloat qi2h   = (PMEFloat)0.50 * qi * qi;
    PMEFloat qid2h  = qi2h * dl;
    sumdeijda_i = -sumdeijda_i + qid2h - cSim.gb_kappa * qi2h * expmkf * reff_i;
    if (cSim.alpb == 0) {
      egb -= qid2h / reff_i;
    }
    else {
      egb -= qid2h * ((PMEFloat)1.0 / reff_i + cSim.one_arad_beta);
      sumdeijda_i *= ((PMEFloat)1.0 + cSim.one_arad_beta * reff_i);
    }
    if (bIGB2578) {

      // New onufriev: we later scale values by alpha-, beta-, and gamma-dependent factor:
      PMEFloat thi, thi2;
      if (cSim.igb == 8) {
        PMEFloat alpha = cSim.pgb_alpha[pos];
        PMEFloat gamma = cSim.pgb_gamma[pos];
        PMEFloat beta  = cSim.pgb_beta[pos];
        thi  = tanh((alpha + gamma * psi_i * psi_i - beta * psi_i) * psi_i);
        thi2 = (alpha + (PMEFloat)3.0 * gamma * psi_i * psi_i -
                (PMEFloat)2.0 * beta * psi_i) *
               ((PMEFloat)1.0 - thi * thi) * (rborn_i - cSim.offset) / rborn_i;
      }
      else {
        thi  = tanh((cSim.gb_alpha + cSim.gb_gamma * psi_i * psi_i - cSim.gb_beta * psi_i) *
                    psi_i);
        thi2 = (cSim.gb_alpha + (PMEFloat)3.0 * cSim.gb_gamma * psi_i * psi_i -
                (PMEFloat)2.0 * cSim.gb_beta * psi_i) *
               ((PMEFloat)1.0 - thi * thi) * (rborn_i - cSim.offset) / rborn_i;
      }
      sumdeijda_i *= thi2;
    }
    cSim.pTemp7[pos] = sumdeijda_i;
    pos += blockDim.x * gridDim.x;
  }

  // Reduce Generalized Born energy
  for (unsigned int stride = warpSize >> 1; stride > 0; stride >>=1) {
    egb += __SHFL_DOWN(WARP_MASK, egb, stride);
  }
  if ((threadIdx.x & GRID_BITS_MASK) == 0) {
    atomicAdd(cSim.pEGB, llitoulli(llrint(egb * ENERGYSCALE)));
  }
}

//---------------------------------------------------------------------------------------------
// kReduceGBTemp7Energy: launch the eponymous kernel
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kReduceGBTemp7Energy(gpuContext gpu)
{
  kReduceGBTemp7Energy_kernel<<<gpu->blocks, gpu->reduceBufferThreadsPerBlock>>>();
  LAUNCHERROR("kReduceGBTemp7Energy");
}
