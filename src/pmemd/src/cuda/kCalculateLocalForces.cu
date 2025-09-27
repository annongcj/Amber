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
#include "bondRemapDS.h"
#include "ptxmacros.h"
#include "hip_definitions.h"

// Use global instance instead of a local copy
#include "simulationConst.h"
CSIM_STO simulationConst cSim;

static __constant__ PMEDouble pt999           = (PMEDouble)(0.9990);
static __constant__ PMEDouble tm06            = (PMEDouble)(1.0e-06);
static __constant__ PMEDouble tenm3           = (PMEDouble)(1.0e-03);
static __constant__ PMEDouble tm24            = (PMEDouble)(1.0e-18);
static __constant__ PMEDouble one             = (PMEDouble)(1.0);
static __constant__ PMEDouble zero            = (PMEDouble)(0.0);
static __constant__ PMEFloat rad_to_deg_coeff = (PMEDouble)180.0 / ((PMEDouble)CMAP_STEP_SIZE *
                                                                    (PMEDouble)PI_VAL);

#if !defined(__HIPCC_RDC__)

//---------------------------------------------------------------------------------------------
// SetkCalculateLocalForcesSim: cast constants and pointers associated with local forces
//                              calculations to the device.  This is called by the oft-invoked
//                              gpuCopyConstants().
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void SetkCalculateLocalForcesSim(gpuContext gpu)
{
  cudaError_t status;
  status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(simulationConst));
  RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

//---------------------------------------------------------------------------------------------
// GetkCalculateLocalForcesSim: get constants and pointers associated with local forces
//                              calculations from the device.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
void GetkCalculateLocalForcesSim(gpuContext gpu)
{
  cudaError_t status;
  status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(simulationConst));
  RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

#endif

//---------------------------------------------------------------------------------------------
// This special version of sincos is designed for |a| < 6*PI. On a GTX 285 it is about 25%
// faster than sincos from the CUDA math library.  It also uses 8 fewer registers than the
// CUDA math library's sincos.  Maximum observed error is 2 ulps across the range stated above.
// Infinities and negative zero are not handled according to C99 specifications.  NaNs are
// handled fine.
//---------------------------------------------------------------------------------------------
__device__ void faster_sincos(PMEDouble a, PMEDouble *sptr, PMEDouble *cptr)
#include "kFastCosSin.h"

//---------------------------------------------------------------------------------------------
// Kernels for bond work units (these will perform all interactions involving pre-specified
// groups of atoms).  These kernels incorporate code in kBWU.h.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnits_kernel(int NMRnstep)
#include "kBWU.h"

//---------------------------------------------------------------------------------------------
#define LOCAL_ENERGY
__global__ void
__LAUNCH_BOUNDS__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsNrg_kernel(int NMRnstep)
#include "kBWU.h"
#define PHMD
__global__ void
__launch_bounds__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsNrgPHMD_kernel(int NMRnstep)
#include "kBWU.h"
#undef PHMD
#undef LOCAL_ENERGY

//---------------------------------------------------------------------------------------------
#define LOCAL_VIRIAL
__global__ void
__LAUNCH_BOUNDS__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsVir_kernel(int NMRnstep)
#include "kBWU.h"

//---------------------------------------------------------------------------------------------
#define LOCAL_ENERGY
__global__ void
__LAUNCH_BOUNDS__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsNrgVir_kernel(int NMRnstep)
#include "kBWU.h"
#undef LOCAL_ENERGY
#undef LOCAL_VIRIAL

//---------------------------------------------------------------------------------------------
#define LOCAL_NEIGHBORLIST
__global__ void
__LAUNCH_BOUNDS__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsNL_kernel(int NMRnstep)
#include "kBWU.h"

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsNLPHMD_kernel(int NMRnstep)
#include "kBWU.h"
#undef PHMD

//---------------------------------------------------------------------------------------------
#define LOCAL_ENERGY
__global__ void
__LAUNCH_BOUNDS__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsNrgNL_kernel(int NMRnstep)
#include "kBWU.h"
#undef LOCAL_ENERGY

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsNrgNLPHMD_kernel(int NMRnstep)
#include "kBWU.h"
#undef PHMD
#undef LOCAL_ENERGY

//---------------------------------------------------------------------------------------------
#define LOCAL_VIRIAL
__global__ void
__LAUNCH_BOUNDS__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsVirNL_kernel(int NMRnstep)
#include "kBWU.h"

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsVirNLPHMD_kernel(int NMRnstep)
#include "kBWU.h"
#undef PHMD

//---------------------------------------------------------------------------------------------
#define LOCAL_ENERGY
__global__ void
__LAUNCH_BOUNDS__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsNrgVirNL_kernel(int NMRnstep)
#include "kBWU.h"

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsNrgVirNLPHMD_kernel(int NMRnstep)
#include "kBWU.h"
#undef PHMD
#undef LOCAL_ENERGY
#undef LOCAL_VIRIAL
#undef LOCAL_NEIGHBORLIST

//---------------------------------------------------------------------------------------------
#define LOCAL_AFE
__global__ void
__LAUNCH_BOUNDS__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsTi_kernel(int NMRnstep)
#include "kBWU.h"

//---------------------------------------------------------------------------------------------
#define LOCAL_ENERGY
__global__ void
__LAUNCH_BOUNDS__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsNrgTi_kernel(int NMRnstep)
#include "kBWU.h"
#undef LOCAL_ENERGY

//---------------------------------------------------------------------------------------------
#define LOCAL_NEIGHBORLIST
__global__ void
__LAUNCH_BOUNDS__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsTiNL_kernel(int NMRnstep)
#include "kBWU.h"

//---------------------------------------------------------------------------------------------
#define LOCAL_ENERGY
__global__ void
__LAUNCH_BOUNDS__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsNrgTiNL_kernel(int NMRnstep)
#include "kBWU.h"
#undef LOCAL_ENERGY
#undef LOCAL_NEIGHBORLIST

#define LOCAL_MBAR
__global__ void
__LAUNCH_BOUNDS__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsTiMbar_kernel(int NMRnstep)
#include "kBWU.h"

//---------------------------------------------------------------------------------------------
#define LOCAL_ENERGY
__global__ void
__LAUNCH_BOUNDS__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsNrgTiMbar_kernel(int NMRnstep)
#include "kBWU.h"
#undef LOCAL_ENERGY

//---------------------------------------------------------------------------------------------
#define LOCAL_NEIGHBORLIST
__global__ void
__LAUNCH_BOUNDS__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsTiMbarNL_kernel(int NMRnstep)
#include "kBWU.h"

//---------------------------------------------------------------------------------------------
#define LOCAL_ENERGY
__global__ void
__LAUNCH_BOUNDS__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsNrgTiMbarNL_kernel(int NMRnstep)
#include "kBWU.h"
#undef LOCAL_ENERGY
#undef LOCAL_NEIGHBORLIST
#undef LOCAL_MBAR
#undef LOCAL_AFE

//---------------------------------------------------------------------------------------------
#define LOCAL_ENERGY
#define AMD_SUM
__global__ void
__LAUNCH_BOUNDS__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsAMD_kernel(int NMRnstep)
#include "kBWU.h"

//---------------------------------------------------------------------------------------------
#define LOCAL_NEIGHBORLIST
__global__ void
__LAUNCH_BOUNDS__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsAMDNL_kernel(int NMRnstep)
#include "kBWU.h"
#undef LOCAL_NEIGHBORLIST
#undef AMD_SUM

//---------------------------------------------------------------------------------------------
#define GAMD_SUM
__global__ void
__LAUNCH_BOUNDS__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsGaMD_kernel(int NMRnstep)
#include "kBWU.h"

//---------------------------------------------------------------------------------------------
#define LOCAL_NEIGHBORLIST
__global__ void
__LAUNCH_BOUNDS__(BOND_WORK_UNIT_THREADS_PER_BLOCK, BOND_WORK_UNIT_BLOCKS_MULTIPLIER)
kExecBondWorkUnitsGaMDNL_kernel(int NMRnstep)
#include "kBWU.h"
#undef LOCAL_NEIGHBORLIST
#undef GAMD_SUM
#undef LOCAL_ENERGY

//---------------------------------------------------------------------------------------------
// kResetBondWorkUnitsCounter_kernel: simply resets the bond work units' counter so that a
//                                    separate calculation, i.e. AMD, can be performed.
//
// Arguments:
//   bondWorkBlocks: the number of bond work blocks that will be launched
//---------------------------------------------------------------------------------------------
__global__ void __LAUNCH_BOUNDS__(32, 1) kResetBondWorkUnitsCounter_kernel(int bondWorkBlocks)
{
  if (threadIdx.x == 0 && blockIdx.x == 0) {
    cSim.pFrcBlkCounters[1] = bondWorkBlocks;
  }
}

//---------------------------------------------------------------------------------------------
// kExecuteBondWorkUnits: host function for calling the appropriate kernel to implement the
//                        pre-compiled bond work units.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, and the energy
//                function
//   calcEnergy:  flag to have the energy calculated.  0 by default, OFF--set to 1 for ON, set
//                to 2 for calculating energy only in select terms for accelerated MD, 3 for
//                calculating energy only in select terms for Gaussian accelerated MD.
//---------------------------------------------------------------------------------------------
extern "C" void kExecuteBondWorkUnits(gpuContext gpu, int calcEnergy)
{
#ifdef AMBER_PLATFORM_AMD
  int nblocks = gpu->sim.bondWorkUnits;
#else
  int nblocks = gpu->bondWorkBlocks;
#endif
  int nthreads = BOND_WORK_UNIT_THREADS_PER_BLOCK;

  // Bail out if there is no work to do--the bLocalInteractions flag is attuned to the presence
  // of a charge mesh, which is initialized to back-fill the bond work units' kernel.
  if (gpu->bLocalInteractions == false) {
    return;
  }
  if (calcEnergy == 1) {
#ifndef AMBER_PLATFORM_AMD
    kResetBondWorkUnitsCounter_kernel<<<1, 32>>>(nblocks);
#endif

    // Energy is calculated in addition to the force and perhaps virial trace
    if (gpu->bNeighborList) {
      if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
        if (gpu->iphmd > 0) {
          kExecBondWorkUnitsNrgVirNLPHMD_kernel<<<nblocks, nthreads>>>(gpu->NMRnstep);
        }
        else {
          kExecBondWorkUnitsNrgVirNL_kernel<<<nblocks, nthreads>>>(gpu->NMRnstep);
        }
      }
      else {
	if (gpu->sim.ti_mode == 0) {
          if (gpu->iphmd > 0) {
            kExecBondWorkUnitsNrgNLPHMD_kernel<<<nblocks, nthreads>>>(gpu->NMRnstep);
          }
          else {
            kExecBondWorkUnitsNrgNL_kernel<<<nblocks, nthreads>>>(gpu->NMRnstep);
          }
        }
        else {
          if (gpu->sim.ifmbar == 1) {
            kExecBondWorkUnitsNrgTiMbarNL_kernel<<<nblocks, nthreads>>>(gpu->NMRnstep);
          }
          else {
	    kExecBondWorkUnitsNrgTiNL_kernel<<<nblocks, nthreads>>>(gpu->NMRnstep);
          }
        }
      }
    }
    else {
      if (gpu->sim.ti_mode == 0) {
       if(gpu->iphmd > 0) {
          kExecBondWorkUnitsNrgPHMD_kernel<<<nblocks, nthreads>>>(gpu->NMRnstep);
       }
        else {
          kExecBondWorkUnitsNrg_kernel<<<nblocks, nthreads>>>(gpu->NMRnstep);
        }
      }
      else {
        if (gpu->sim.ifmbar == 1) {
          kExecBondWorkUnitsNrgTiMbar_kernel<<<nblocks, nthreads>>>(gpu->NMRnstep);
        }
        else {
          kExecBondWorkUnitsNrgTi_kernel<<<nblocks, nthreads>>>(gpu->NMRnstep);
        }
      }
    }
#ifdef AMBER_PLATFORM_AMD
#ifdef use_DPFP
    cudaMemsetAsync(gpu->sim.plliXYZ_q, 0,
                    gpu->sim.XYZStride * sizeof(long long int));
#else
    cudaMemsetAsync(gpu->sim.plliXYZ_q, 0,
                    gpu->sim.XYZStride * sizeof(int));
#endif
#endif
  }
  else if (calcEnergy == 2) {
#ifndef AMBER_PLATFORM_AMD
    kResetBondWorkUnitsCounter_kernel<<<1, 32>>>(nblocks);
#endif

    // Calculate only the energy associated with aMD-associated terms,
    // stashing the result in a special accumulator.
    if (gpu->bNeighborList) {
      kExecBondWorkUnitsAMDNL_kernel<<<nblocks, nthreads>>>(gpu->NMRnstep);
    }
    else {
      kExecBondWorkUnitsAMD_kernel<<<nblocks, nthreads>>>(gpu->NMRnstep);
    }
  }
  else if (calcEnergy == 3) {
#ifndef AMBER_PLATFORM_AMD
    kResetBondWorkUnitsCounter_kernel<<<1, 32>>>(nblocks);
#endif

    // Calculate only the energy associated with Gaussian aMD-associated terms,
    // stashing the result in a special accumulator.
    if (gpu->bNeighborList) {
      kExecBondWorkUnitsGaMDNL_kernel<<<nblocks, nthreads>>>(gpu->NMRnstep);
    }
    else {
      kExecBondWorkUnitsGaMD_kernel<<<nblocks, nthreads>>>(gpu->NMRnstep);
    }
  }
  else {

    // Only the force and perhaps virial trace will be computed
    if (gpu->bNeighborList) {
      if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
        if (gpu->iphmd > 0) {
          kExecBondWorkUnitsVirNLPHMD_kernel<<<nblocks, nthreads>>>(gpu->NMRnstep);
        }
        else {
          kExecBondWorkUnitsVirNL_kernel<<<nblocks, nthreads>>>(gpu->NMRnstep);
        }
      }
      else {
        if (gpu->sim.ti_mode == 0) {
          if (gpu->iphmd > 0) {
            kExecBondWorkUnitsNLPHMD_kernel<<<nblocks, nthreads>>>(gpu->NMRnstep);
          }
          else {
            kExecBondWorkUnitsNL_kernel<<<nblocks, nthreads>>>(gpu->NMRnstep);
          }
        }
        else {
          if (gpu->sim.ifmbar == 1) {
	    kExecBondWorkUnitsTiMbarNL_kernel<<<nblocks, nthreads>>>(gpu->NMRnstep);
          }
          else {
            kExecBondWorkUnitsTiNL_kernel<<<nblocks, nthreads>>>(gpu->NMRnstep);
          }
        }
      }
    }
    else {
      if (gpu->sim.ti_mode == 0) {
        kExecBondWorkUnits_kernel<<<nblocks, nthreads>>>(gpu->NMRnstep);
      }
      else {
        if (gpu->sim.ifmbar == 1) {
          kExecBondWorkUnitsTiMbar_kernel<<<nblocks, nthreads>>>(gpu->NMRnstep);
        }
        else {
          kExecBondWorkUnitsTi_kernel<<<nblocks, nthreads>>>(gpu->NMRnstep);
        }
      }
    }
#ifdef AMBER_PLATFORM_AMD
#ifdef use_DPFP
    cudaMemsetAsync(gpu->sim.plliXYZ_q, 0,
                    gpu->sim.XYZStride * sizeof(long long int));
#else
    cudaMemsetAsync(gpu->sim.plliXYZ_q, 0,
                    gpu->sim.XYZStride * sizeof(int));
#endif
#endif
  }
}

//---------------------------------------------------------------------------------------------
// NMR interactions kernels: there are separate instances for forces, forces and energies,
//                           with another twist about whether it happens in the context of a
//                           PME calculation (with NMR_NEIGHBORLIST defined).  A final detail
//                           distinguishing the kernels is whether double precision, accurate
//                           texture memory is available.  Code in kCNF.h is included here.
//---------------------------------------------------------------------------------------------
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcNMRFrc_kernel(int step)
#include "kCNF.h"

//---------------------------------------------------------------------------------------------
#define R6AV
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcNMRR6avFrc_kernel(int step)
#include "kCNF.h"
#undef R6AV

//---------------------------------------------------------------------------------------------
#define NMR_ENERGY
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcNMRNrg_kernel(int step)
#include "kCNF.h"

//---------------------------------------------------------------------------------------------
#define R6AV
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcNMRR6avNrg_kernel(int step)
#include "kCNF.h"
#undef R6AV
#undef NMR_ENERGY

//---------------------------------------------------------------------------------------------
#define NMR_NEIGHBORLIST
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcPMENMRFrc_kernel(int step)
#include "kCNF.h"

//---------------------------------------------------------------------------------------------
#define R6AV
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcPMENMRR6avFrc_kernel(int step)
#include "kCNF.h"
#undef R6AV

//---------------------------------------------------------------------------------------------
#define NMR_ENERGY
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcPMENMRNrg_kernel(int step)
#include "kCNF.h"

//---------------------------------------------------------------------------------------------
#define R6AV
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcPMENMRR6avNrg_kernel(int step)
#include "kCNF.h"
#undef R6AV
#undef NMR_ENERGY
#undef NMR_NEIGHBORLIST

//---------------------------------------------------------------------------------------------
#define NODPTEXTURE
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcNMRFrcNoDPTex_kernel(int step)
#include "kCNF.h"

//---------------------------------------------------------------------------------------------
#define R6AV
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcNMRR6avFrcNoDPTex_kernel(int step)
#include "kCNF.h"
#undef R6AV

//---------------------------------------------------------------------------------------------
#define NMR_ENERGY
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcNMRNrgNoDPTex_kernel(int step)
#include "kCNF.h"

//---------------------------------------------------------------------------------------------
#define R6AV
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcNMRR6avNrgNoDPTex_kernel(int step)
#include "kCNF.h"
#undef R6AV
#undef NMR_ENERGY

//---------------------------------------------------------------------------------------------
#define NMR_NEIGHBORLIST
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcPMENMRFrcNoDPTex_kernel(int step)
#include "kCNF.h"

//---------------------------------------------------------------------------------------------
#define R6AV
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcPMENMRR6avFrcNoDPTex_kernel(int step)
#include "kCNF.h"
#undef R6AV

//---------------------------------------------------------------------------------------------
#define NMR_ENERGY
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcPMENMRNrgNoDPTex_kernel(int step)
#include "kCNF.h"

//---------------------------------------------------------------------------------------------
#define R6AV
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcPMENMRR6avNrgNoDPTex_kernel(int step)
#include "kCNF.h"
#undef R6AV
#undef NMR_ENERGY
#undef NMR_NEIGHBORLIST
#undef NODPTEXTURE

//---------------------------------------------------------------------------------------------
// Alchemical free energy (AFE) NMR kernels: for now just TI, no MBAR
// we don't need specific force kernels because the restraint forces
// are independent of lambda, as is dvdl. we need these kernels to get
// region specific energies for AFE
//---------------------------------------------------------------------------------------------
#define NMR_AFE
#define NMR_ENERGY
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcNMRNrgAFE_kernel(int step)
#include "kCNF.h"

//---------------------------------------------------------------------------------------------
#define R6AV
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcNMRR6avNrgAFE_kernel(int step)
#include "kCNF.h"
#undef R6AV

//---------------------------------------------------------------------------------------------
#define NMR_NEIGHBORLIST
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcPMENMRNrgAFE_kernel(int step)
#include "kCNF.h"

//---------------------------------------------------------------------------------------------
#define R6AV
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcPMENMRR6avNrgAFE_kernel(int step)
#include "kCNF.h"
#undef R6AV
#undef NMR_NEIGHBORLIST

//---------------------------------------------------------------------------------------------
#define NODPTEXTURE
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcNMRNrgNoDPTexAFE_kernel(int step)
#include "kCNF.h"

//---------------------------------------------------------------------------------------------
#define R6AV
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcNMRR6avNrgNoDPTexAFE_kernel(int step)
#include "kCNF.h"
#undef R6AV

//---------------------------------------------------------------------------------------------
#define NMR_NEIGHBORLIST
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcPMENMRNrgNoDPTexAFE_kernel(int step)
#include "kCNF.h"

//---------------------------------------------------------------------------------------------
#define R6AV
__LAUNCH_BOUNDS__(NMRFORCES_THREADS_PER_BLOCK, NMRFORCES_BLOCKS)
__global__ void kCalcPMENMRR6avNrgNoDPTexAFE_kernel(int step)
#include "kCNF.h"
#undef R6AV

#undef NMR_ENERGY
#undef NMR_NEIGHBORLIST
#undef NODPTEXTURE
#undef NMR_AFE

//---------------------------------------------------------------------------------------------
// kCalculateNMRForces: compute the forces due to NMR constraints.  The vectorization of these
//                      things is a bit different than the vectorization of other aspects of
//                      the local interactions.
//
// Arguments:
//   gpu:    overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kCalculateNMRForces(gpuContext gpu)
{
  // Set the block dimension in a local variable.  The number of blocks for each
  // of a series of kernel launchs is determined inside the loop below.
  int nFrcThreads = gpu->NMRForcesThreadsPerBlock;
  int totalBlocks  = gpu->NMRForcesBlocks;
  int launchBlocks = 65535;
  if (gpu->bNMRInteractions) {
    if (gpu->bNoDPTexture) {
      while (totalBlocks > 0) {
        int blocks = min(totalBlocks, launchBlocks);
        if (gpu->bNeighborList)
          if(gpu->sim.NMRR6av) {
              kCalcPMENMRR6avFrcNoDPTex_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
          }
          else {
              kCalcPMENMRFrcNoDPTex_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
          }
        else
          if(gpu->sim.NMRR6av) {
              kCalcNMRR6avFrcNoDPTex_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
          }
          else {
              kCalcNMRFrcNoDPTex_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
          }
        LAUNCHERROR("kCalculateNMRForcesNoDPTexture");
        totalBlocks -= blocks;
      }
    }
    else {
      while (totalBlocks > 0) {
        int blocks = min(totalBlocks, launchBlocks);
        if (gpu->bNeighborList) {
          if(gpu->sim.NMRR6av) {
              kCalcPMENMRR6avFrc_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
          }
          else {
              kCalcPMENMRFrc_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
          }
        }
        else {
          if(gpu->sim.NMRR6av) {
              kCalcNMRR6avFrc_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
          }
          else {
              kCalcNMRFrc_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
          }
        }
        LAUNCHERROR("kCalculateNMRForces");
        totalBlocks -= blocks;
      }
    }
  }
}

//---------------------------------------------------------------------------------------------
// kCalculateNMREnergy: compute the energy due to NMR constraints.  The vectorization of these
//                      things is a bit different than the vectorization of other aspects of
//                      the local interactions.
//
// Arguments:
//   gpu:    overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kCalculateNMREnergy(gpuContext gpu)
{
  // Set the block dimension in a local variable.  The number of blocks for each
  // of a series of kernel launchs is determined inside the loop below.
  int nFrcThreads = gpu->NMRForcesThreadsPerBlock;
  int totalBlocks  = gpu->NMRForcesBlocks;
  int launchBlocks = 65535;
  if (gpu->bNMRInteractions) {
    if (gpu->bNoDPTexture) {
      while (totalBlocks > 0) {
        int blocks = min(totalBlocks, launchBlocks);
        if (gpu->bNeighborList) {
          if (gpu->sim.ti_mode == 0) {
            if(gpu->sim.NMRR6av) {
                kCalcPMENMRR6avNrgNoDPTex_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
            }
            else {
                kCalcPMENMRNrgNoDPTex_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
            }
          }
          else {
            if(gpu->sim.NMRR6av) {
                kCalcPMENMRR6avNrgNoDPTexAFE_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
            }
            else {
                kCalcPMENMRNrgNoDPTexAFE_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
            }
          }
        }
        else {
          if (gpu->sim.ti_mode == 0) {
            if(gpu->sim.NMRR6av) {
                kCalcNMRR6avNrgNoDPTex_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
            }
            else {
                kCalcNMRNrgNoDPTex_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
            }
          }
          else {
            if(gpu->sim.NMRR6av) {
                kCalcNMRR6avNrgNoDPTexAFE_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
            }
            else {
                kCalcNMRNrgNoDPTexAFE_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
            }
          }
        }
        LAUNCHERROR("kCalculateNMREnergyNoDPTexture");
        totalBlocks -= launchBlocks;
      }
    }
    else {
      while (totalBlocks > 0) {
        int blocks = min(totalBlocks, launchBlocks);
        if (gpu->bNeighborList) {
          if (gpu->sim.ti_mode == 0) {
            if(gpu->sim.NMRR6av) {
                kCalcPMENMRR6avNrg_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
            }
            else {
              kCalcPMENMRNrg_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
            }
          }
          else {
            if(gpu->sim.NMRR6av) {
                kCalcPMENMRR6avNrgAFE_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
            }
            else {
                kCalcPMENMRNrgAFE_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
            }
          }
        }
        else {
          if (gpu->sim.ti_mode == 0) {
            if(gpu->sim.NMRR6av) {
                kCalcNMRR6avNrg_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
            }
            else {
                kCalcNMRNrg_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
            }
          }
          else {
            if(gpu->sim.NMRR6av) {
                kCalcNMRR6avNrgAFE_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
            }
            else {
                kCalcNMRNrgAFE_kernel<<<blocks, nFrcThreads>>>(gpu->NMRnstep);
            }
          }
        }
        LAUNCHERROR("kCalculateNMREnergy");
        totalBlocks -= launchBlocks;
      }
    }
  }
}

//---------------------------------------------------------------------------------------------
// kCalcAMDWeightAndScaleFrc_kernel: see the host function below for a better description of
//                                   what this does.
//
// Arguments:
//   pot_ene_tot: total potential energy (passed in from gpu_calculate_and_apply_amd_weights_
//                in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the aMD weight
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kCalcAMDWeightAndScaleFrc_kernel(PMEDouble pot_ene_tot, PMEDouble dih_ene_tot, PMEDouble fwgt)
{
  // Calculate AMD weight, seting dihedral boost (tboost) to zero for now
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  if ((cSim.ntp > 0) && (cSim.barostat == 1)) {
    while (pos < cSim.atoms) {
      PMEDouble forceX = cSim.pForceXAccumulator[pos];
      PMEDouble forceY = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ = cSim.pForceZAccumulator[pos];
      PMEDouble NBForceX = cSim.pNBForceXAccumulator[pos];
      PMEDouble NBForceY = cSim.pNBForceYAccumulator[pos];
      PMEDouble NBForceZ = cSim.pNBForceZAccumulator[pos];
      forceX   *= fwgt;
      forceY   *= fwgt;
      forceZ   *= fwgt;
      NBForceX *= fwgt;
      NBForceY *= fwgt;
      NBForceZ *= fwgt;
      cSim.pForceXAccumulator[pos]   = forceX;
      cSim.pForceYAccumulator[pos]   = forceY;
      cSim.pForceZAccumulator[pos]   = forceZ;
      cSim.pNBForceXAccumulator[pos] = NBForceX;
      cSim.pNBForceYAccumulator[pos] = NBForceY;
      cSim.pNBForceZAccumulator[pos] = NBForceZ;
      pos += increment;
    }
  }
  else {
    while (pos < cSim.atoms) {
      PMEDouble forceX = cSim.pForceXAccumulator[pos];
      PMEDouble forceY = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ = cSim.pForceZAccumulator[pos];
      forceX *= fwgt;
      forceY *= fwgt;
      forceZ *= fwgt;
      cSim.pForceXAccumulator[pos] = forceX;
      cSim.pForceYAccumulator[pos] = forceY;
      cSim.pForceZAccumulator[pos] = forceZ;
      pos += increment;
    }
  }
}

//---------------------------------------------------------------------------------------------
// kCalculateAMDWeightAndScaleForces: launch the kernel to compute the weights for accelerated
//                                    MD and scale the forces appropriately.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, the energy
//                function, and the aMD settings
//   pot_ene_tot: total potential energy (passed in from the calling function
//                                        gpu_calculate_and_apply_amd_weights_ in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the aMD weight
//---------------------------------------------------------------------------------------------
void kCalculateAMDWeightAndScaleForces(gpuContext gpu, PMEDouble pot_ene_tot,
                                       PMEDouble dih_ene_tot, PMEDouble fwgt)
{
  // Get the tile size into shorter variable names: this is extremely picky but
  // I DO NOT like line wrapping, and these names may help to show that the tile
  // size for GAMD weighting and force scaling is indeed different than the
  // tile size for the local force computations themselves under GaMD.
  int nBlocks = gpu->blocks;
  int nGenThreads = gpu->generalThreadsPerBlock;

  kCalcAMDWeightAndScaleFrc_kernel<<<nBlocks, nGenThreads>>>(pot_ene_tot, dih_ene_tot, fwgt);
  LAUNCHERROR("kCalculateAMDWeightAndScaleForces");
}

//---------------------------------------------------------------------------------------------
// kCalcGAMDWeightAndScaleFrc_kernel: see the host function below for a better description of
//                                    what this does.
//
// Arguments:
//   pot_ene_tot: total potential energy (passed in from gpu_calculate_and_apply_gamd_weights_
//                in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
__global__ void
__launch_bounds__(GENERAL_THREADS_PER_BLOCK, 1)
kCalcGAMDWeightAndScaleFrc_kernel(PMEDouble pot_ene_tot, PMEDouble dih_ene_tot, PMEDouble fwgt)
{
  // Calculate GAMD weight, setting dihedral boost (tboost) to zero for now
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  if ((cSim.ntp > 0) && (cSim.barostat == 1)) {
    while (pos < cSim.atoms) {
      PMEDouble forceX    = cSim.pForceXAccumulator[pos];
      PMEDouble forceY    = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ    = cSim.pForceZAccumulator[pos];
      PMEDouble NBForceX  = cSim.pNBForceXAccumulator[pos];
      PMEDouble NBForceY  = cSim.pNBForceYAccumulator[pos];
      PMEDouble NBForceZ  = cSim.pNBForceZAccumulator[pos];
      forceX             *= fwgt;
      forceY             *= fwgt;
      forceZ             *= fwgt;
      NBForceX           *= fwgt;
      NBForceY           *= fwgt;
      NBForceZ           *= fwgt;
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
      cSim.pNBForceXAccumulator[pos]  = NBForceX;
      cSim.pNBForceYAccumulator[pos]  = NBForceY;
      cSim.pNBForceZAccumulator[pos]  = NBForceZ;
      pos += increment;
    }
  }
  else {
    while (pos < cSim.atoms) {
      PMEDouble forceX   = cSim.pForceXAccumulator[pos];
      PMEDouble forceY   = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ   = cSim.pForceZAccumulator[pos];
      forceX            *= fwgt;
      forceY            *= fwgt;
      forceZ            *= fwgt;
      cSim.pForceXAccumulator[pos] = forceX;
      cSim.pForceYAccumulator[pos] = forceY;
      cSim.pForceZAccumulator[pos] = forceZ;
      pos += increment;
    }
  }
}

//---------------------------------------------------------------------------------------------
// kCalcGAMDWeightAndScaleFrc_nb_kernel: see the host function below for a better description
//                                       of what this does.
//
// Arguments:
//   pot_ene_tot: total potential energy (passed in from gpu_calculate_and_apply_gamd_weights_
//                in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kCalcGAMDWeightAndScaleFrc_nb_kernel(PMEDouble pot_ene_nb, PMEDouble dih_ene_tot,
                                     PMEDouble fwgt)
{
  // Calculate GAMD weight, setting dihedral boost (tboost) to zero for now
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  if ((cSim.ntp > 0) && (cSim.barostat == 1)) {
    while (pos < cSim.atoms) {
      PMEDouble NBForceX  = cSim.pNBForceXAccumulator[pos];
      PMEDouble NBForceY  = cSim.pNBForceYAccumulator[pos];
      PMEDouble NBForceZ  = cSim.pNBForceZAccumulator[pos];
      NBForceX           *= fwgt;
      NBForceY           *= fwgt;
      NBForceZ           *= fwgt;
      cSim.pNBForceXAccumulator[pos]  = NBForceX;
      cSim.pNBForceYAccumulator[pos]  = NBForceY;
      cSim.pNBForceZAccumulator[pos]  = NBForceZ;
      pos += increment;
    }
  }
}

//---------------------------------------------------------------------------------------------
// kCalcGAMDWeightAndScaleFrc_gb_kernel: see the host function below for a better description of
//                                    what this does.
//
// Arguments:
//   pot_ene_tot: total potential energy (passed in from gpu_calculate_and_apply_gamd_weights_
//                in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
__global__ void
__launch_bounds__(GENERAL_THREADS_PER_BLOCK, 1)
kCalcGAMDWeightAndScaleFrc_gb_kernel(PMEDouble pot_ene_tot, PMEDouble dih_ene_tot, PMEDouble fwgt, bool useImage)
{
  // Calculate GAMD weight, setting dihedral boost (tboost) to zero for now
  unsigned int pos0 = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  while (pos0 < cSim.atoms) {
    unsigned int pos = ( useImage ) ? cSim.pImageAtomLookup[pos0] : pos0;
    PMEDouble forceX    = cSim.pForceXAccumulator[pos];
    PMEDouble forceY    = cSim.pForceYAccumulator[pos];
    PMEDouble forceZ    = cSim.pForceZAccumulator[pos];
    forceX             *= fwgt;
    forceY             *= fwgt;
    forceZ             *= fwgt;
    cSim.pForceXAccumulator[pos]    = forceX;
    cSim.pForceYAccumulator[pos]    = forceY;
    cSim.pForceZAccumulator[pos]    = forceZ;
    pos0 += increment;
  }
}

//---------------------------------------------------------------------------------------------
// NOTE: kCalcGAMDWeightAndScaleFrc_gb_nb_kernel should not be used for the moment
// as pNBForce is not saved separately with GB.
//
// kCalcGAMDWeightAndScaleFrc_gb_nb_kernel: see the host function below for a better description
//                                       of what this does.
//
// Arguments:
//   pot_ene_tot: total potential energy (passed in from gpu_calculate_and_apply_gamd_weights_
//                in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
__global__ void
__launch_bounds__(GENERAL_THREADS_PER_BLOCK, 1)
kCalcGAMDWeightAndScaleFrc_gb_nb_kernel(PMEDouble pot_ene_nb, PMEDouble dih_ene_tot,
                                     PMEDouble fwgt, bool useImage)
{
  // Calculate GAMD weight, setting dihedral boost (tboost) to zero for now
  unsigned int pos0 = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

    while (pos0 < cSim.atoms) {
      unsigned int pos = ( useImage ) ? cSim.pImageAtomLookup[pos0] : pos0;
      PMEDouble NBForceX  = cSim.pNBForceXAccumulator[pos];
      PMEDouble NBForceY  = cSim.pNBForceYAccumulator[pos];
      PMEDouble NBForceZ  = cSim.pNBForceZAccumulator[pos];
      NBForceX           *= fwgt;
      NBForceY           *= fwgt;
      NBForceZ           *= fwgt;
      cSim.pNBForceXAccumulator[pos]  = NBForceX;
      cSim.pNBForceYAccumulator[pos]  = NBForceY;
      cSim.pNBForceZAccumulator[pos]  = NBForceZ;
      pos0 += increment;
    }
}

//---------------------------------------------------------------------------------------------
// kCalcGAMDWeightAndScaleFrc_bond_kernel: see the host function below for a better description
//                                       of what this does.
//
// Arguments:
//   pot_ene_tot: total potential energy (passed in from gpu_calculate_and_apply_gamd_weights_
//                in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kCalcGAMDWeightAndScaleFrc_bond_kernel(PMEDouble pot_ene_bond, PMEDouble dih_ene_tot,
                                     PMEDouble fwgt,bool useImage)
{
  // Calculate GAMD weight, setting dihedral boost (tboost) to zero for now
  unsigned int pos0 = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  if ((cSim.ntp > 0) && (cSim.barostat == 1)) {
    while (pos0 < cSim.atoms) {
      unsigned int pos = ( useImage ) ? cSim.pImageAtomLookup[pos0] : pos0;
      PMEDouble BForceX  = cSim.pForceXAccumulator[pos];
      PMEDouble BForceY  = cSim.pForceYAccumulator[pos];
      PMEDouble BForceZ  = cSim.pForceZAccumulator[pos];
      BForceX           *= fwgt;
      BForceY           *= fwgt;
      BForceZ           *= fwgt;
      cSim.pForceXAccumulator[pos]  = BForceX;
      cSim.pForceYAccumulator[pos]  = BForceY;
      cSim.pForceZAccumulator[pos]  = BForceZ;
      pos0 += increment;
    }
  }
  else {
    while (pos0 < cSim.atoms) {
      unsigned int pos = ( useImage ) ? cSim.pImageAtomLookup[pos0] : pos0;
      PMEDouble forceX   = cSim.pForceXAccumulator[pos];
      PMEDouble forceY   = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ   = cSim.pForceZAccumulator[pos];
      forceX            *= fwgt;
      forceY            *= fwgt;
      forceZ            *= fwgt;
      cSim.pForceXAccumulator[pos] = forceX;
      cSim.pForceYAccumulator[pos] = forceY;
      cSim.pForceZAccumulator[pos] = forceZ;
      pos0 += increment;
    }
  }
}

//---------------------------------------------------------------------------------------------
// kCalcGAMDWeightAndScaleFrc_ti_kernel: see the host function below for a better description
//                                       of what this does.
//
// Arguments:
//   pot_ene_tot: total potential energy (passed in from gpu_calculate_and_apply_gamd_weights_
//                in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kCalcGAMDWeightAndScaleFrc_ti_kernel(PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgt, bool useImage)
{
  // Calculate GAMD weight, setting dihedral boost (tboost) to zero for now
  unsigned int pos0 = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  // scale atomic force in SC region
  while (pos0 < cSim.atoms) {
    unsigned int pos = ( useImage ) ? cSim.pImageAtomLookup[pos0] : pos0;
    if ((cSim.ntp > 0) && (cSim.barostat == 1)) {
      PMEDouble forceX    = cSim.pForceXAccumulator[pos];
      PMEDouble forceY    = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ    = cSim.pForceZAccumulator[pos];
      PMEDouble NBForceX  = cSim.pNBForceXAccumulator[pos];
      PMEDouble NBForceY  = cSim.pNBForceYAccumulator[pos];
      PMEDouble NBForceZ  = cSim.pNBForceZAccumulator[pos];
      PMEDouble TIForceX = cSim.pTIForceX[0][pos];
      PMEDouble TIForceY = cSim.pTIForceY[0][pos];
      PMEDouble TIForceZ = cSim.pTIForceZ[0][pos];
      PMEDouble TISCForceX = cSim.pTISCForceX[0][pos];
      PMEDouble TISCForceY = cSim.pTISCForceY[0][pos];
      PMEDouble TISCForceZ = cSim.pTISCForceZ[0][pos];
      PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
      PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
      PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
      forceX              = forceX - (TIForceX + TISCForceX) * (1.0-fwgt);
      forceY              = forceY - (TIForceY + TISCForceY) * (1.0-fwgt);
      forceZ              = forceZ - (TIForceZ + TISCForceZ) * (1.0-fwgt);
      NBForceX            = NBForceX - GaMDTINBForceX * (1.0-fwgt);
      NBForceY            = NBForceY - GaMDTINBForceY * (1.0-fwgt);
      NBForceZ            = NBForceZ - GaMDTINBForceZ * (1.0-fwgt);
      TIForceX           *= fwgt;
      TIForceY           *= fwgt;
      TIForceZ           *= fwgt;
      TISCForceX           *= fwgt;
      TISCForceY           *= fwgt;
      TISCForceZ           *= fwgt;
      GaMDTINBForceX           *= fwgt;
      GaMDTINBForceY           *= fwgt;
      GaMDTINBForceZ           *= fwgt;
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
      cSim.pNBForceXAccumulator[pos]  = NBForceX;
      cSim.pNBForceYAccumulator[pos]  = NBForceY;
      cSim.pNBForceZAccumulator[pos]  = NBForceZ;
      cSim.pTIForceX[0][pos]  = TIForceX;
      cSim.pTIForceY[0][pos]  = TIForceY;
      cSim.pTIForceZ[0][pos]  = TIForceZ;
      cSim.pTISCForceX[0][pos]  = TISCForceX;
      cSim.pTISCForceY[0][pos]  = TISCForceY;
      cSim.pTISCForceZ[0][pos]  = TISCForceZ;
      cSim.pGaMDTINBForceX[0][pos]  = GaMDTINBForceX;
      cSim.pGaMDTINBForceY[0][pos]  = GaMDTINBForceY;
      cSim.pGaMDTINBForceZ[0][pos]  = GaMDTINBForceZ;
    } else {
      PMEDouble forceX   = cSim.pForceXAccumulator[pos];
      PMEDouble forceY   = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ   = cSim.pForceZAccumulator[pos];
      PMEDouble TIForceX = cSim.pTIForceX[0][pos];
      PMEDouble TIForceY = cSim.pTIForceY[0][pos];
      PMEDouble TIForceZ = cSim.pTIForceZ[0][pos];
      PMEDouble TISCForceX = cSim.pTISCForceX[0][pos];
      PMEDouble TISCForceY = cSim.pTISCForceY[0][pos];
      PMEDouble TISCForceZ = cSim.pTISCForceZ[0][pos];
      PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
      PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
      PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
      forceX              = forceX - (TIForceX + TISCForceX + GaMDTINBForceX) * (1.0-fwgt);
      forceY              = forceY - (TIForceY + TISCForceY + GaMDTINBForceY) * (1.0-fwgt);
      forceZ              = forceZ - (TIForceZ + TISCForceZ + GaMDTINBForceZ) * (1.0-fwgt);
      TIForceX           *= fwgt;
      TIForceY           *= fwgt;
      TIForceZ           *= fwgt;
      TISCForceX           *= fwgt;
      TISCForceY           *= fwgt;
      TISCForceZ           *= fwgt;
      GaMDTINBForceX           *= fwgt;
      GaMDTINBForceY           *= fwgt;
      GaMDTINBForceZ           *= fwgt;
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
      cSim.pTIForceX[0][pos]  = TIForceX;
      cSim.pTIForceY[0][pos]  = TIForceY;
      cSim.pTIForceZ[0][pos]  = TIForceZ;
      cSim.pTISCForceX[0][pos]  = TISCForceX;
      cSim.pTISCForceY[0][pos]  = TISCForceY;
      cSim.pTISCForceZ[0][pos]  = TISCForceZ;
      cSim.pGaMDTINBForceX[0][pos]  = GaMDTINBForceX;
      cSim.pGaMDTINBForceY[0][pos]  = GaMDTINBForceY;
      cSim.pGaMDTINBForceZ[0][pos]  = GaMDTINBForceZ;
    }
    pos0 += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kCalcGAMDWeightAndScaleFrc_sc_others_ppi_kernel: see the host function below for a better description
//                                       of what this does.
//
// Arguments:
//   pot_ene_tot: total potential energy (passed in from gpu_calculate_and_apply_gamd_weights_
//                in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kCalcGAMDWeightAndScaleFrc_sc_others_ppi_kernel(PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgtd, bool useImage)
{
  // Calculate GAMD weight, setting dihedral boost (tboost) to zero for now
  unsigned int pos0 = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  // scale atomic force in SC region
  while (pos0 < cSim.atoms) {
    unsigned int pos = ( useImage ) ? cSim.pImageAtomLookup[pos0] : pos0;
    if ((cSim.ntp > 0) && (cSim.barostat == 1)) {
      PMEDouble forceX    = cSim.pForceXAccumulator[pos];
      PMEDouble forceY    = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ    = cSim.pForceZAccumulator[pos];
      PMEDouble NBForceX  = cSim.pNBForceXAccumulator[pos];
      PMEDouble NBForceY  = cSim.pNBForceYAccumulator[pos];
      PMEDouble NBForceZ  = cSim.pNBForceZAccumulator[pos];
      PMEDouble GaMDTIForceX = cSim.pGaMDTIForceX[0][pos];
      PMEDouble GaMDTIForceY = cSim.pGaMDTIForceY[0][pos];
      PMEDouble GaMDTIForceZ = cSim.pGaMDTIForceZ[0][pos];
      PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
      PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
      PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
      forceX              = (forceX - GaMDTIForceX) * fwgtd + GaMDTIForceX;
      forceY              = (forceY - GaMDTIForceY) * fwgtd + GaMDTIForceY;
      forceZ              = (forceZ - GaMDTIForceZ) * fwgtd + GaMDTIForceZ;
      NBForceX            = (NBForceX - GaMDTINBForceX) * fwgtd + GaMDTINBForceX;
      NBForceY            = (NBForceY - GaMDTINBForceY) * fwgtd + GaMDTINBForceY;
      NBForceZ            = (NBForceZ - GaMDTINBForceZ) * fwgtd + GaMDTINBForceZ;
      GaMDTIForceX       *= fwgtd;
      GaMDTIForceX       *= fwgtd;
      GaMDTIForceX       *= fwgtd;
      GaMDTINBForceX           *= fwgtd;
      GaMDTINBForceY           *= fwgtd;
      GaMDTINBForceZ           *= fwgtd;
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
      cSim.pNBForceXAccumulator[pos]  = NBForceX;
      cSim.pNBForceYAccumulator[pos]  = NBForceY;
      cSim.pNBForceZAccumulator[pos]  = NBForceZ;
//      cSim.pGaMDTIForceX[0][pos]  = GaMDTIForceX;
//      cSim.pGaMDTIForceY[0][pos]  = GaMDTIForceY;
//      cSim.pGaMDTIForceZ[0][pos]  = GaMDTIForceZ;
//      cSim.pGaMDTINBForceX[0][pos]  = GaMDTINBForceX;
//      cSim.pGaMDTINBForceY[0][pos]  = GaMDTINBForceY;
//      cSim.pGaMDTINBForceZ[0][pos]  = GaMDTINBForceZ;
    } else {
      PMEDouble forceX   = cSim.pForceXAccumulator[pos];
      PMEDouble forceY   = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ   = cSim.pForceZAccumulator[pos];
      PMEDouble GaMDTIForceX = cSim.pGaMDTIForceX[0][pos];
      PMEDouble GaMDTIForceY = cSim.pGaMDTIForceY[0][pos];
      PMEDouble GaMDTIForceZ = cSim.pGaMDTIForceZ[0][pos];
      PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
      PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
      PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
//      forceX              = forceX - (GaMDTIForceX + GaMDTINBForceX)* (1.0-fwgt);
//      forceY              = forceY - (GaMDTIForceY + GaMDTINBForceY)* (1.0-fwgt);
//      forceZ              = forceZ - (GaMDTIForceZ + GaMDTINBForceZ)* (1.0-fwgt);
      forceX              = (forceX - GaMDTIForceX) * fwgtd + GaMDTIForceX;
      forceY              = (forceY - GaMDTIForceY) * fwgtd + GaMDTIForceY;
      forceZ              = (forceZ - GaMDTIForceZ) * fwgtd + GaMDTIForceZ;
//      GaMDTIForceX       *= fwgtd;
//      GaMDTIForceX       *= fwgtd;
//      GaMDTIForceX       *= fwgtd;
//      GaMDTINBForceX           *= fwgtd;
//      GaMDTINBForceY           *= fwgtd;
//      GaMDTINBForceZ           *= fwgtd;
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
//      cSim.pGaMDTIForceX[0][pos]  = GaMDTIForceX;
//      cSim.pGaMDTIForceY[0][pos]  = GaMDTIForceY;
//      cSim.pGaMDTIForceZ[0][pos]  = GaMDTIForceZ;
//      cSim.pGaMDTINBForceX[0][pos]  = GaMDTINBForceX;
//     cSim.pGaMDTINBForceY[0][pos]  = GaMDTINBForceY;
//      cSim.pGaMDTINBForceZ[0][pos]  = GaMDTINBForceZ;
    }
    pos0 += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kCalcGAMDWeightAndScaleFrc_sc_ppi_kernel: see the host function below for a better description
//                                       of what this does.
//
// Arguments:
//   pot_ene_tot: total potential energy (passed in from gpu_calculate_and_apply_gamd_weights_
//                in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kCalcGAMDWeightAndScaleFrc_sc_ppi_kernel(PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgt, bool useImage)
{
  // Calculate GAMD weight, setting dihedral boost (tboost) to zero for now
  unsigned int pos0 = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  // scale atomic force in SC region
  while (pos0 < cSim.atoms) {
    unsigned int pos = ( useImage ) ? cSim.pImageAtomLookup[pos0] : pos0;
    if ((cSim.ntp > 0) && (cSim.barostat == 1)) {
      PMEDouble forceX    = cSim.pForceXAccumulator[pos];
      PMEDouble forceY    = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ    = cSim.pForceZAccumulator[pos];
      PMEDouble NBForceX  = cSim.pNBForceXAccumulator[pos];
      PMEDouble NBForceY  = cSim.pNBForceYAccumulator[pos];
      PMEDouble NBForceZ  = cSim.pNBForceZAccumulator[pos];
      PMEDouble GaMDTIForceX = cSim.pGaMDTIForceX[0][pos];
      PMEDouble GaMDTIForceY = cSim.pGaMDTIForceY[0][pos];
      PMEDouble GaMDTIForceZ = cSim.pGaMDTIForceZ[0][pos];
      PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
      PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
      PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
      forceX              = forceX - GaMDTIForceX * (1.0-fwgt);
      forceY              = forceY - GaMDTIForceY * (1.0-fwgt);
      forceZ              = forceZ - GaMDTIForceZ * (1.0-fwgt);
      NBForceX            = NBForceX - GaMDTINBForceX * (1.0-fwgt);
      NBForceY            = NBForceY - GaMDTINBForceY * (1.0-fwgt);
      NBForceZ            = NBForceZ - GaMDTINBForceZ * (1.0-fwgt);
      GaMDTIForceX       *= fwgt;
      GaMDTIForceX       *= fwgt;
      GaMDTIForceX       *= fwgt;
      GaMDTINBForceX           *= fwgt;
      GaMDTINBForceY           *= fwgt;
      GaMDTINBForceZ           *= fwgt;
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
      cSim.pNBForceXAccumulator[pos]  = NBForceX;
      cSim.pNBForceYAccumulator[pos]  = NBForceY;
      cSim.pNBForceZAccumulator[pos]  = NBForceZ;
      cSim.pGaMDTIForceX[0][pos]  = GaMDTIForceX;
      cSim.pGaMDTIForceY[0][pos]  = GaMDTIForceY;
      cSim.pGaMDTIForceZ[0][pos]  = GaMDTIForceZ;
      cSim.pGaMDTINBForceX[0][pos]  = GaMDTINBForceX;
      cSim.pGaMDTINBForceY[0][pos]  = GaMDTINBForceY;
      cSim.pGaMDTINBForceZ[0][pos]  = GaMDTINBForceZ;
    } else {
      PMEDouble forceX   = cSim.pForceXAccumulator[pos];
      PMEDouble forceY   = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ   = cSim.pForceZAccumulator[pos];
      PMEDouble GaMDTIForceX = cSim.pGaMDTIForceX[0][pos];
      PMEDouble GaMDTIForceY = cSim.pGaMDTIForceY[0][pos];
      PMEDouble GaMDTIForceZ = cSim.pGaMDTIForceZ[0][pos];
      PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
      PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
      PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
//      forceX              = forceX - (GaMDTIForceX + GaMDTINBForceX)* (1.0-fwgt);
//      forceY              = forceY - (GaMDTIForceY + GaMDTINBForceY)* (1.0-fwgt);
//      forceZ              = forceZ - (GaMDTIForceZ + GaMDTINBForceZ)* (1.0-fwgt);
      forceX              = forceX - GaMDTIForceX * (1.0-fwgt);
      forceY              = forceY - GaMDTIForceY * (1.0-fwgt);
      forceZ              = forceZ - GaMDTIForceZ * (1.0-fwgt);
      GaMDTIForceX       *= fwgt;
      GaMDTIForceX       *= fwgt;
      GaMDTIForceX       *= fwgt;
      GaMDTINBForceX           *= fwgt;
      GaMDTINBForceY           *= fwgt;
      GaMDTINBForceZ           *= fwgt;
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
      cSim.pGaMDTIForceX[0][pos]  = GaMDTIForceX;
      cSim.pGaMDTIForceY[0][pos]  = GaMDTIForceY;
      cSim.pGaMDTIForceZ[0][pos]  = GaMDTIForceZ;
      cSim.pGaMDTINBForceX[0][pos]  = GaMDTINBForceX;
      cSim.pGaMDTINBForceY[0][pos]  = GaMDTINBForceY;
      cSim.pGaMDTINBForceZ[0][pos]  = GaMDTINBForceZ;
    }
    pos0 += increment;
  }
}


//---------------------------------------------------------------------------------------------
// kCalcGAMDWeightAndScaleFrc_sc_bonded_ppi_kernel: see the host function below for a better description
//                                       of what this does.
//
// Arguments:
//   pot_ene_tot: total potential energy (passed in from gpu_calculate_and_apply_gamd_weights_
//                in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kCalcGAMDWeightAndScaleFrc_sc_bonded_ppi_kernel(PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgt, bool useImage)
{
  // Calculate GAMD weight, setting dihedral boost (tboost) to zero for now
  unsigned int pos0 = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  // scale SC bonded atomic forces
  while (pos0 < cSim.atoms) {
    unsigned int pos = ( useImage ) ? cSim.pImageAtomLookup[pos0] : pos0;
    PMEDouble forceX    = cSim.pForceXAccumulator[pos];
    PMEDouble forceY    = cSim.pForceYAccumulator[pos];
    PMEDouble forceZ    = cSim.pForceZAccumulator[pos];
    PMEDouble GaMDTIForceX = cSim.pGaMDTIForceX[0][pos];
    PMEDouble GaMDTIForceY = cSim.pGaMDTIForceY[0][pos];
    PMEDouble GaMDTIForceZ = cSim.pGaMDTIForceZ[0][pos];
    forceX              = forceX - GaMDTIForceX * (1.0-fwgt);
    forceY              = forceY - GaMDTIForceY * (1.0-fwgt);
    forceZ              = forceZ - GaMDTIForceZ * (1.0-fwgt);
    GaMDTIForceX           *= fwgt;
    GaMDTIForceY           *= fwgt;
    GaMDTIForceZ           *= fwgt;
    cSim.pForceXAccumulator[pos]    = forceX;
    cSim.pForceYAccumulator[pos]    = forceY;
    cSim.pForceZAccumulator[pos]    = forceZ;
    cSim.pGaMDTIForceX[0][pos]  = GaMDTIForceX;
    cSim.pGaMDTIForceY[0][pos]  = GaMDTIForceY;
    cSim.pGaMDTIForceZ[0][pos]  = GaMDTIForceZ;

    pos0 += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kCalcGAMDWeightAndScaleFrc_ti_bonded_kernel: see the host function below for a better description
//                                       of what this does.
//
// Arguments:
//   pot_ene_tot: total potential energy (passed in from gpu_calculate_and_apply_gamd_weights_
//                in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kCalcGAMDWeightAndScaleFrc_ti_bonded_kernel(PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgt, bool useImage)
{
  // Calculate GAMD weight, setting dihedral boost (tboost) to zero for now
  unsigned int pos0 = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  // scale SC bonded atomic forces
  while (pos0 < cSim.atoms) {
    unsigned int pos = ( useImage ) ? cSim.pImageAtomLookup[pos0] : pos0;
    PMEDouble forceX    = cSim.pForceXAccumulator[pos];
    PMEDouble forceY    = cSim.pForceYAccumulator[pos];
    PMEDouble forceZ    = cSim.pForceZAccumulator[pos];
    PMEDouble TIForceX = cSim.pTIForceX[0][pos];
    PMEDouble TIForceY = cSim.pTIForceY[0][pos];
    PMEDouble TIForceZ = cSim.pTIForceZ[0][pos];
    PMEDouble TISCForceX = cSim.pTISCForceX[0][pos];
    PMEDouble TISCForceY = cSim.pTISCForceY[0][pos];
    PMEDouble TISCForceZ = cSim.pTISCForceZ[0][pos];
    forceX              = forceX - (TIForceX + TISCForceX) * (1.0-fwgt);
    forceY              = forceY - (TIForceY + TISCForceY) * (1.0-fwgt);
    forceZ              = forceZ - (TIForceZ + TISCForceZ) * (1.0-fwgt);
    TIForceX           *= fwgt;
    TIForceY           *= fwgt;
    TIForceZ           *= fwgt;
    TISCForceX           *= fwgt;
    TISCForceY           *= fwgt;
    TISCForceZ           *= fwgt;
    cSim.pForceXAccumulator[pos]    = forceX;
    cSim.pForceYAccumulator[pos]    = forceY;
    cSim.pForceZAccumulator[pos]    = forceZ;
    cSim.pTIForceX[0][pos]  = TIForceX;
    cSim.pTIForceY[0][pos]  = TIForceY;
    cSim.pTIForceZ[0][pos]  = TIForceZ;
    cSim.pTISCForceX[0][pos]  = TISCForceX;
    cSim.pTISCForceY[0][pos]  = TISCForceY;
    cSim.pTISCForceZ[0][pos]  = TISCForceZ;
    pos0 += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kCalcGAMDWeightAndScaleFrc_sc_nonbonded_kernel: see the host function below for a better description
//                                       of what this does.
//
// Arguments:
//   pot_ene_tot: total potential energy (passed in from gpu_calculate_and_apply_gamd_weights_
//                in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kCalcGAMDWeightAndScaleFrc_sc_nonbonded_kernel(PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgt, bool useImage)
{
  // Calculate GAMD weight, setting dihedral boost (tboost) to zero for now
  unsigned int pos0 = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  // scale SC nonbonded atomic force
  while (pos0 < cSim.atoms) {
    unsigned int pos = ( useImage ) ? cSim.pImageAtomLookup[pos0] : pos0;
    if ((cSim.ntp > 0) && (cSim.barostat == 1)) {
      PMEDouble NBForceX  = cSim.pNBForceXAccumulator[pos];
      PMEDouble NBForceY  = cSim.pNBForceYAccumulator[pos];
      PMEDouble NBForceZ  = cSim.pNBForceZAccumulator[pos];
      PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
      PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
      PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
      NBForceX            = NBForceX - GaMDTINBForceX * (1.0-fwgt);
      NBForceY            = NBForceY - GaMDTINBForceY * (1.0-fwgt);
      NBForceZ            = NBForceZ - GaMDTINBForceZ * (1.0-fwgt);
      GaMDTINBForceX           *= fwgt;
      GaMDTINBForceY           *= fwgt;
      GaMDTINBForceZ           *= fwgt;
      cSim.pNBForceXAccumulator[pos]  = NBForceX;
      cSim.pNBForceYAccumulator[pos]  = NBForceY;
      cSim.pNBForceZAccumulator[pos]  = NBForceZ;
      cSim.pGaMDTINBForceX[0][pos]  = GaMDTINBForceX;
      cSim.pGaMDTINBForceY[0][pos]  = GaMDTINBForceY;
      cSim.pGaMDTINBForceZ[0][pos]  = GaMDTINBForceZ;
    } else {
      PMEDouble forceX   = cSim.pForceXAccumulator[pos];
      PMEDouble forceY   = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ   = cSim.pForceZAccumulator[pos];
      PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
      PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
      PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
      forceX              = forceX - GaMDTINBForceX * (1.0-fwgt);
      forceY              = forceY - GaMDTINBForceY * (1.0-fwgt);
      forceZ              = forceZ - GaMDTINBForceZ * (1.0-fwgt);
      GaMDTINBForceX           *= fwgt;
      GaMDTINBForceY           *= fwgt;
      GaMDTINBForceZ           *= fwgt;
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
      cSim.pGaMDTINBForceX[0][pos]  = GaMDTINBForceX;
      cSim.pGaMDTINBForceY[0][pos]  = GaMDTINBForceY;
      cSim.pGaMDTINBForceZ[0][pos]  = GaMDTINBForceZ;
    }
    pos0 += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kCalcGAMDWeightAndScaleFrc_ti_nonbonded_kernel: see the host function below for a better description
//                                       of what this does.
//
// Arguments:
//   pot_ene_tot: total potential energy (passed in from gpu_calculate_and_apply_gamd_weights_
//                in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kCalcGAMDWeightAndScaleFrc_ti_nonbonded_kernel(PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgt, bool useImage)
{
  // Calculate GAMD weight, setting dihedral boost (tboost) to zero for now
  unsigned int pos0 = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  // scale SC nonbonded atomic force
  while (pos0 < cSim.atoms) {
    unsigned int pos = ( useImage ) ? cSim.pImageAtomLookup[pos0] : pos0;
    if ((cSim.ntp > 0) && (cSim.barostat == 1)) {
      PMEDouble NBForceX  = cSim.pNBForceXAccumulator[pos];
      PMEDouble NBForceY  = cSim.pNBForceYAccumulator[pos];
      PMEDouble NBForceZ  = cSim.pNBForceZAccumulator[pos];
      PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
      PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
      PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
      NBForceX            = NBForceX - GaMDTINBForceX * (1.0-fwgt);
      NBForceY            = NBForceY - GaMDTINBForceY * (1.0-fwgt);
      NBForceZ            = NBForceZ - GaMDTINBForceZ * (1.0-fwgt);
      GaMDTINBForceX           *= fwgt;
      GaMDTINBForceY           *= fwgt;
      GaMDTINBForceZ           *= fwgt;
      cSim.pNBForceXAccumulator[pos]  = NBForceX;
      cSim.pNBForceYAccumulator[pos]  = NBForceY;
      cSim.pNBForceZAccumulator[pos]  = NBForceZ;
      cSim.pGaMDTINBForceX[0][pos]  = GaMDTINBForceX;
      cSim.pGaMDTINBForceY[0][pos]  = GaMDTINBForceY;
      cSim.pGaMDTINBForceZ[0][pos]  = GaMDTINBForceZ;
    } else {
      PMEDouble forceX   = cSim.pForceXAccumulator[pos];
      PMEDouble forceY   = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ   = cSim.pForceZAccumulator[pos];
      PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
      PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
      PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
      forceX              = forceX - GaMDTINBForceX * (1.0-fwgt);
      forceY              = forceY - GaMDTINBForceY * (1.0-fwgt);
      forceZ              = forceZ - GaMDTINBForceZ * (1.0-fwgt);
      GaMDTINBForceX           *= fwgt;
      GaMDTINBForceY           *= fwgt;
      GaMDTINBForceZ           *= fwgt;
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
      cSim.pGaMDTINBForceX[0][pos]  = GaMDTINBForceX;
      cSim.pGaMDTINBForceY[0][pos]  = GaMDTINBForceY;
      cSim.pGaMDTINBForceZ[0][pos]  = GaMDTINBForceZ;
    }
    pos0 += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kCalcGAMDWeightAndScaleFrc_ti_others_kernel: see the host function below for a better description
//                                       of what this does.
//
// Arguments:
//   pot_ene_tot: total potential energy (passed in from gpu_calculate_and_apply_gamd_weights_
//                in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kCalcGAMDWeightAndScaleFrc_ti_others_kernel(PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgtd, bool useImage)
{
  // Calculate GAMD weight, setting dihedral boost (tboost) to zero for now
  unsigned int pos0 = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  // scale atomic force in SC region
  while (pos0 < cSim.atoms) {
    unsigned int pos = ( useImage ) ? cSim.pImageAtomLookup[pos0] : pos0;
    if ((cSim.ntp > 0) && (cSim.barostat == 1)) {
      PMEDouble forceX    = cSim.pForceXAccumulator[pos];
      PMEDouble forceY    = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ    = cSim.pForceZAccumulator[pos];
      PMEDouble NBForceX  = cSim.pNBForceXAccumulator[pos];
      PMEDouble NBForceY  = cSim.pNBForceYAccumulator[pos];
      PMEDouble NBForceZ  = cSim.pNBForceZAccumulator[pos];
      PMEDouble TIForceX = cSim.pTIForceX[0][pos];
      PMEDouble TIForceY = cSim.pTIForceY[0][pos];
      PMEDouble TIForceZ = cSim.pTIForceZ[0][pos];
      PMEDouble TISCForceX = cSim.pTISCForceX[0][pos];
      PMEDouble TISCForceY = cSim.pTISCForceY[0][pos];
      PMEDouble TISCForceZ = cSim.pTISCForceZ[0][pos];
      PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
      PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
      PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
      forceX              = (forceX - TIForceX - TISCForceX) * fwgtd + TIForceX + TISCForceX;
      forceY              = (forceY - TIForceY - TISCForceY) * fwgtd + TIForceY + TISCForceY;
      forceZ              = (forceZ - TIForceZ - TISCForceZ) * fwgtd + TIForceZ + TISCForceZ;
      NBForceX            = (NBForceX - GaMDTINBForceX) * fwgtd + GaMDTINBForceX;
      NBForceY            = (NBForceY - GaMDTINBForceY) * fwgtd + GaMDTINBForceY;
      NBForceZ            = (NBForceZ - GaMDTINBForceZ) * fwgtd + GaMDTINBForceZ;
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
      cSim.pNBForceXAccumulator[pos]  = NBForceX;
      cSim.pNBForceYAccumulator[pos]  = NBForceY;
      cSim.pNBForceZAccumulator[pos]  = NBForceZ;
    } else {
      PMEDouble forceX   = cSim.pForceXAccumulator[pos];
      PMEDouble forceY   = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ   = cSim.pForceZAccumulator[pos];
      PMEDouble TIForceX = cSim.pTIForceX[0][pos];
      PMEDouble TIForceY = cSim.pTIForceY[0][pos];
      PMEDouble TIForceZ = cSim.pTIForceZ[0][pos];
      PMEDouble TISCForceX = cSim.pTISCForceX[0][pos];
      PMEDouble TISCForceY = cSim.pTISCForceY[0][pos];
      PMEDouble TISCForceZ = cSim.pTISCForceZ[0][pos];
      PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
      PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
      PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
      forceX            = (forceX - TIForceX - TISCForceX - GaMDTINBForceX) * fwgtd + TIForceX + TISCForceX + GaMDTINBForceX;
      forceY            = (forceY - TIForceY - TISCForceY - GaMDTINBForceY) * fwgtd + TIForceY + TISCForceY + GaMDTINBForceY;
      forceZ            = (forceZ - TIForceZ - TISCForceZ - GaMDTINBForceZ) * fwgtd + TIForceZ + TISCForceZ + GaMDTINBForceZ;
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
    }
    pos0 += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kCalcGAMDWeightAndScaleFrc_ti_bonded_others_kernel: see the host function below for a better description
//                                       of what this does.
//
// Arguments:
//   pot_ene_tot: total potential energy (passed in from gpu_calculate_and_apply_gamd_weights_
//                in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kCalcGAMDWeightAndScaleFrc_ti_bonded_others_kernel(PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgtd, bool useImage)
{
  // Calculate GAMD weight, setting dihedral boost (tboost) to zero for now
  unsigned int pos0 = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  // scale atomic force in SC region
  while (pos0 < cSim.atoms) {
    unsigned int pos = ( useImage ) ? cSim.pImageAtomLookup[pos0] : pos0;
    if ((cSim.ntp > 0) && (cSim.barostat == 1)) {
      PMEDouble forceX    = cSim.pForceXAccumulator[pos];
      PMEDouble forceY    = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ    = cSim.pForceZAccumulator[pos];
      PMEDouble NBForceX  = cSim.pNBForceXAccumulator[pos];
      PMEDouble NBForceY  = cSim.pNBForceYAccumulator[pos];
      PMEDouble NBForceZ  = cSim.pNBForceZAccumulator[pos];
      PMEDouble TIForceX = cSim.pTIForceX[0][pos];
      PMEDouble TIForceY = cSim.pTIForceY[0][pos];
      PMEDouble TIForceZ = cSim.pTIForceZ[0][pos];
      PMEDouble TISCForceX = cSim.pTISCForceX[0][pos];
      PMEDouble TISCForceY = cSim.pTISCForceY[0][pos];
      PMEDouble TISCForceZ = cSim.pTISCForceZ[0][pos];
      forceX              = (forceX - TIForceX - TISCForceX) * fwgtd + TIForceX + TISCForceX;
      forceY              = (forceY - TIForceY - TISCForceY) * fwgtd + TIForceY + TISCForceY;
      forceZ              = (forceZ - TIForceZ - TISCForceZ) * fwgtd + TIForceZ + TISCForceZ;
      NBForceX            = NBForceX * fwgtd;
      NBForceY            = NBForceY * fwgtd;
      NBForceZ            = NBForceZ * fwgtd;
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
      cSim.pNBForceXAccumulator[pos]  = NBForceX;
      cSim.pNBForceYAccumulator[pos]  = NBForceY;
      cSim.pNBForceZAccumulator[pos]  = NBForceZ;
    } else {
      PMEDouble forceX   = cSim.pForceXAccumulator[pos];
      PMEDouble forceY   = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ   = cSim.pForceZAccumulator[pos];
      PMEDouble TIForceX = cSim.pTIForceX[0][pos];
      PMEDouble TIForceY = cSim.pTIForceY[0][pos];
      PMEDouble TIForceZ = cSim.pTIForceZ[0][pos];
      PMEDouble TISCForceX = cSim.pTISCForceX[0][pos];
      PMEDouble TISCForceY = cSim.pTISCForceY[0][pos];
      PMEDouble TISCForceZ = cSim.pTISCForceZ[0][pos];
      forceX              = (forceX - TIForceX - TISCForceX) * fwgtd + TIForceX + TISCForceX;
      forceY              = (forceY - TIForceY - TISCForceY) * fwgtd + TIForceY + TISCForceY;
      forceZ              = (forceZ - TIForceZ - TISCForceZ) * fwgtd + TIForceZ + TISCForceZ;
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
    }
    pos0 += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kCalcGAMDWeightAndScaleFrc_ti_nonbonded_others_kernel: see the host function below for a better description
//                                       of what this does.
//
// Arguments:
//   pot_ene_tot: total potential energy (passed in from gpu_calculate_and_apply_gamd_weights_
//                in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kCalcGAMDWeightAndScaleFrc_ti_nonbonded_others_kernel(PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgtd, bool useImage)
{
  // Calculate GAMD weight, setting dihedral boost (tboost) to zero for now
  unsigned int pos0 = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  // scale atomic force in SC region
  while (pos0 < cSim.atoms) {
    unsigned int pos = ( useImage ) ? cSim.pImageAtomLookup[pos0] : pos0;
    if ((cSim.ntp > 0) && (cSim.barostat == 1)) {
      PMEDouble forceX    = cSim.pForceXAccumulator[pos];
      PMEDouble forceY    = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ    = cSim.pForceZAccumulator[pos];
      PMEDouble NBForceX  = cSim.pNBForceXAccumulator[pos];
      PMEDouble NBForceY  = cSim.pNBForceYAccumulator[pos];
      PMEDouble NBForceZ  = cSim.pNBForceZAccumulator[pos];
      PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
      PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
      PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
      forceX              = forceX * fwgtd;
      forceY              = forceY * fwgtd;
      forceZ              = forceZ * fwgtd;
      NBForceX            = (NBForceX - GaMDTINBForceX) * fwgtd + GaMDTINBForceX;
      NBForceY            = (NBForceY - GaMDTINBForceY) * fwgtd + GaMDTINBForceY;
      NBForceZ            = (NBForceZ - GaMDTINBForceZ) * fwgtd + GaMDTINBForceZ;
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
      cSim.pNBForceXAccumulator[pos]  = NBForceX;
      cSim.pNBForceYAccumulator[pos]  = NBForceY;
      cSim.pNBForceZAccumulator[pos]  = NBForceZ;
    } else {
      PMEDouble forceX   = cSim.pForceXAccumulator[pos];
      PMEDouble forceY   = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ   = cSim.pForceZAccumulator[pos];
      PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
      PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
      PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
      forceX            = (forceX - GaMDTINBForceX) * fwgtd + GaMDTINBForceX;
      forceY            = (forceY - GaMDTINBForceY) * fwgtd + GaMDTINBForceY;
      forceZ            = (forceZ - GaMDTINBForceZ) * fwgtd + GaMDTINBForceZ;
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
    }
    pos0 += increment;
  }
}
//---------------------------------------------------------------------------------------------
// kCalcGAMDWeightAndScaleFrc_sc_dual_nonbonded_kernel: see the host function below for a better description
//                                       of what this does.
//
// Arguments:
//   pot_ene_tot: total potential energy (passed in from gpu_calculate_and_apply_gamd_weights_
//                in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kCalcGAMDWeightAndScaleFrc_sc_dual_nonbonded_kernel(PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgt, PMEDouble fwgtd, bool useImage)
{
  // Calculate GAMD weight, setting dihedral boost (tboost) to zero for now
  unsigned int pos0 = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  // scale atomic force in SC region
  while (pos0 < cSim.atoms) {
    unsigned int pos = ( useImage ) ? cSim.pImageAtomLookup[pos0] : pos0;
    if ((cSim.ntp > 0) && (cSim.barostat == 1)) {
      PMEDouble forceX    = cSim.pForceXAccumulator[pos];
      PMEDouble forceY    = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ    = cSim.pForceZAccumulator[pos];
      PMEDouble NBForceX  = cSim.pNBForceXAccumulator[pos];
      PMEDouble NBForceY  = cSim.pNBForceYAccumulator[pos];
      PMEDouble NBForceZ  = cSim.pNBForceZAccumulator[pos];
      PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
      PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
      PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
      forceX              = forceX * fwgtd;
      forceY              = forceY * fwgtd;
      forceZ              = forceZ * fwgtd;
      NBForceX            = NBForceX - (NBForceX - GaMDTINBForceX) * (1-fwgtd) - GaMDTINBForceX * (1-fwgt);
      NBForceY            = NBForceY - (NBForceY - GaMDTINBForceY) * (1-fwgtd) - GaMDTINBForceY * (1-fwgt);
      NBForceZ            = NBForceZ - (NBForceZ - GaMDTINBForceZ) * (1-fwgtd) - GaMDTINBForceZ * (1-fwgt);
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
      cSim.pNBForceXAccumulator[pos]  = NBForceX;
      cSim.pNBForceYAccumulator[pos]  = NBForceY;
      cSim.pNBForceZAccumulator[pos]  = NBForceZ;
    } else {
      PMEDouble forceX   = cSim.pForceXAccumulator[pos];
      PMEDouble forceY   = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ   = cSim.pForceZAccumulator[pos];
      PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
      PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
      PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
      forceX              = forceX - (forceX - GaMDTINBForceX) * (1-fwgtd) - GaMDTINBForceX*(1-fwgt);
      forceY              = forceY - (forceY - GaMDTINBForceY) * (1-fwgtd) - GaMDTINBForceY*(1-fwgt);
      forceZ              = forceZ - (forceZ - GaMDTINBForceZ) * (1-fwgtd) - GaMDTINBForceZ*(1-fwgt);
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
    }
    pos0 += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kCalcGAMDWeightAndScaleFrc_sc_dual_bonded_kernel: see the host function below for a better description
//                                       of what this does.
//
// Arguments:
//   pot_ene_tot: total potential energy (passed in from gpu_calculate_and_apply_gamd_weights_
//                in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kCalcGAMDWeightAndScaleFrc_sc_dual_bonded_kernel(PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgt, PMEDouble fwgtd, bool useImage)
{
  // Calculate GAMD weight, setting dihedral boost (tboost) to zero for now
  unsigned int pos0 = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  // scale atomic force in SC region
  while (pos0 < cSim.atoms) {
    unsigned int pos = ( useImage ) ? cSim.pImageAtomLookup[pos0] : pos0;
    if ((cSim.ntp > 0) && (cSim.barostat == 1)) {
      PMEDouble forceX    = cSim.pForceXAccumulator[pos];
      PMEDouble forceY    = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ    = cSim.pForceZAccumulator[pos];
      PMEDouble NBForceX  = cSim.pNBForceXAccumulator[pos];
      PMEDouble NBForceY  = cSim.pNBForceYAccumulator[pos];
      PMEDouble NBForceZ  = cSim.pNBForceZAccumulator[pos];
      PMEDouble GaMDTIForceX = cSim.pGaMDTIForceX[0][pos];
      PMEDouble GaMDTIForceY = cSim.pGaMDTIForceY[0][pos];
      PMEDouble GaMDTIForceZ = cSim.pGaMDTIForceZ[0][pos];
      PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
      PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
      PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
      forceX              = forceX - (forceX-GaMDTIForceX) * (1-fwgtd) - GaMDTIForceX * (1-fwgt);
      forceY              = forceY - (forceY-GaMDTIForceY) * (1-fwgtd) - GaMDTIForceY * (1-fwgt);
      forceZ              = forceZ - (forceZ-GaMDTIForceZ) * (1-fwgtd) - GaMDTIForceZ * (1-fwgt);
      NBForceX            = NBForceX - (NBForceX - GaMDTINBForceX) *(1-fwgtd) - GaMDTINBForceX * (1-fwgt);
      NBForceY            = NBForceY - (NBForceY - GaMDTINBForceY) *(1-fwgtd) - GaMDTINBForceY * (1-fwgt);
      NBForceZ            = NBForceZ - (NBForceZ - GaMDTINBForceZ) *(1-fwgtd) - GaMDTINBForceZ * (1-fwgt);
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
      cSim.pNBForceXAccumulator[pos]  = NBForceX;
      cSim.pNBForceYAccumulator[pos]  = NBForceY;
      cSim.pNBForceZAccumulator[pos]  = NBForceZ;
    } else {
      PMEDouble forceX   = cSim.pForceXAccumulator[pos];
      PMEDouble forceY   = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ   = cSim.pForceZAccumulator[pos];
      PMEDouble GaMDTIForceX = cSim.pGaMDTIForceX[0][pos];
      PMEDouble GaMDTIForceY = cSim.pGaMDTIForceY[0][pos];
      PMEDouble GaMDTIForceZ = cSim.pGaMDTIForceZ[0][pos];
      PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
      PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
      PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
      forceX              = forceX - (forceX - GaMDTIForceX - GaMDTINBForceX) * (1-fwgtd) - (GaMDTIForceX+GaMDTINBForceX)*(1-fwgt);
      forceY              = forceY - (forceY - GaMDTIForceY - GaMDTINBForceY) * (1-fwgtd) - (GaMDTIForceY+GaMDTINBForceY)*(1-fwgt);
      forceZ              = forceZ - (forceZ - GaMDTIForceZ - GaMDTINBForceZ) * (1-fwgtd) - (GaMDTIForceZ+GaMDTINBForceZ)*(1-fwgt);
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
    }
    pos0 += increment;
  }
}
//---------------------------------------------------------------------------------------------
// kCalcGAMDWeightAndScaleFrc_sc_triple_kernel: see the host function below for a better description
//                                       of what this does.
//
// Arguments:
//   pot_ene_tot: total potential energy (passed in from gpu_calculate_and_apply_gamd_weights_
//                in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
// work with igamd 22
//---------------------------------------------------------------------------------------------
__global__ void
__launch_bounds__(GENERAL_THREADS_PER_BLOCK, 1)
kCalcGAMDWeightAndScaleFrc_sc_triple_kernel(PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgt, PMEDouble fwgtd, PMEDouble fwgtb,bool useImage)
{
  // Calculate GAMD weight, setting dihedral boost (tboost) to zero for now
  unsigned int pos0 = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  // scale atomic force in SC region
  while (pos0 < cSim.atoms) {
    unsigned int pos = ( useImage ) ? cSim.pImageAtomLookup[pos0] : pos0;
    if ((cSim.ntp > 0) && (cSim.barostat == 1)) {
      PMEDouble forceX    = cSim.pForceXAccumulator[pos];
      PMEDouble forceY    = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ    = cSim.pForceZAccumulator[pos];
      PMEDouble NBForceX  = cSim.pNBForceXAccumulator[pos];
      PMEDouble NBForceY  = cSim.pNBForceYAccumulator[pos];
      PMEDouble NBForceZ  = cSim.pNBForceZAccumulator[pos];
      PMEDouble GaMDTIForceX = cSim.pGaMDTIForceX[0][pos];
      PMEDouble GaMDTIForceY = cSim.pGaMDTIForceY[0][pos];
      PMEDouble GaMDTIForceZ = cSim.pGaMDTIForceZ[0][pos];
      PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
      PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
      PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
      forceX              = forceX - (forceX-GaMDTIForceX) * (1-fwgtd) - GaMDTIForceX * (1-fwgtb);
      forceY              = forceY - (forceY-GaMDTIForceY) * (1-fwgtd) - GaMDTIForceY * (1-fwgtb);
      forceZ              = forceZ - (forceZ-GaMDTIForceZ) * (1-fwgtd) - GaMDTIForceZ * (1-fwgtb);
      NBForceX            = NBForceX - (NBForceX - GaMDTINBForceX) *(1-fwgtd) - GaMDTINBForceX * (1-fwgt);
      NBForceY            = NBForceY - (NBForceY - GaMDTINBForceY) *(1-fwgtd) - GaMDTINBForceY * (1-fwgt);
      NBForceZ            = NBForceZ - (NBForceZ - GaMDTINBForceZ) *(1-fwgtd) - GaMDTINBForceZ * (1-fwgt);
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
      cSim.pNBForceXAccumulator[pos]  = NBForceX;
      cSim.pNBForceYAccumulator[pos]  = NBForceY;
      cSim.pNBForceZAccumulator[pos]  = NBForceZ;
    } else {
      PMEDouble forceX   = cSim.pForceXAccumulator[pos];
      PMEDouble forceY   = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ   = cSim.pForceZAccumulator[pos];
      PMEDouble GaMDTIForceX = cSim.pGaMDTIForceX[0][pos];
      PMEDouble GaMDTIForceY = cSim.pGaMDTIForceY[0][pos];
      PMEDouble GaMDTIForceZ = cSim.pGaMDTIForceZ[0][pos];
      PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
      PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
      PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
      forceX              = forceX - (forceX - GaMDTIForceX - GaMDTINBForceX) * (1-fwgtd) - GaMDTINBForceX*(1-fwgt)-GaMDTIForceX*(1-fwgtb);
      forceY              = forceY - (forceY - GaMDTIForceY - GaMDTINBForceY) * (1-fwgtd) - GaMDTINBForceY*(1-fwgt)-GaMDTIForceY*(1-fwgtb);
      forceZ              = forceZ - (forceZ - GaMDTIForceZ - GaMDTINBForceZ) * (1-fwgtd) - GaMDTINBForceZ*(1-fwgt)-GaMDTIForceZ*(1-fwgtb);
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
    }
    pos0 += increment;
  }
}
//---------------------------------------------------------------------------------------------
// kCalcGAMDWeightAndScaleFrc_sc_triple_kernel: see the host function below for a better description
//                                       of what this does.
// only works for igamd 21 and 26
// Arguments:
//   pot_ene_tot: total potential energy (passed in from gpu_calculate_and_apply_gamd_weights_
//                in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
__global__ void
__launch_bounds__(GENERAL_THREADS_PER_BLOCK, 1)
kCalcGAMDWeightAndScaleFrc_sc_triple2_kernel(PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgt, PMEDouble fwgtd, PMEDouble fwgtb,bool useImage)
{
  // Calculate GAMD weight, setting dihedral boost (tboost) to zero for now
  unsigned int pos0 = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  // scale atomic force in SC region
  while (pos0 < cSim.atoms) {
    unsigned int pos = ( useImage ) ? cSim.pImageAtomLookup[pos0] : pos0;
    if ((cSim.ntp > 0) && (cSim.barostat == 1)) {
      PMEDouble forceX    = cSim.pForceXAccumulator[pos];
      PMEDouble forceY    = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ    = cSim.pForceZAccumulator[pos];
      PMEDouble NBForceX  = cSim.pNBForceXAccumulator[pos];
      PMEDouble NBForceY  = cSim.pNBForceYAccumulator[pos];
      PMEDouble NBForceZ  = cSim.pNBForceZAccumulator[pos];
      PMEDouble GaMDTIForceX = cSim.pGaMDTIForceX[0][pos];
      PMEDouble GaMDTIForceY = cSim.pGaMDTIForceY[0][pos];
      PMEDouble GaMDTIForceZ = cSim.pGaMDTIForceZ[0][pos];
      PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
      PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
      PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
      forceX              = forceX*fwgtb;
      forceY              = forceY*fwgtb;
      forceZ              = forceZ*fwgtb;
      NBForceX            = NBForceX - (NBForceX - GaMDTINBForceX) *(1-fwgtd) - GaMDTINBForceX * (1-fwgt);
      NBForceY            = NBForceY - (NBForceY - GaMDTINBForceY) *(1-fwgtd) - GaMDTINBForceY * (1-fwgt);
      NBForceZ            = NBForceZ - (NBForceZ - GaMDTINBForceZ) *(1-fwgtd) - GaMDTINBForceZ * (1-fwgt);
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
      cSim.pNBForceXAccumulator[pos]  = NBForceX;
      cSim.pNBForceYAccumulator[pos]  = NBForceY;
      cSim.pNBForceZAccumulator[pos]  = NBForceZ;
    }
    pos0 += increment;
  }
}
//---------------------------------------------------------------------------------------------
// kCalcGAMDWeightAndScaleFrc_sc_triple3_kernel: see the host function below for a better description
//                                       of what this does.
//
// Arguments:
//   pot_ene_tot: total potential energy (passed in from gpu_calculate_and_apply_gamd_weights_
//                in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
// work with igamd eq 20 and 23
//---------------------------------------------------------------------------------------------
__global__ void
__launch_bounds__(GENERAL_THREADS_PER_BLOCK, 1)
kCalcGAMDWeightAndScaleFrc_sc_triple3_kernel(PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgt, PMEDouble fwgtd, PMEDouble fwgtb,bool useImage)
{
  // Calculate GAMD weight, setting dihedral boost (tboost) to zero for now
  unsigned int pos0 = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  // scale atomic force in SC region
  while (pos0 < cSim.atoms) {
    unsigned int pos = ( useImage ) ? cSim.pImageAtomLookup[pos0] : pos0;
    if ((cSim.ntp > 0) && (cSim.barostat == 1)) {
      PMEDouble forceX    = cSim.pForceXAccumulator[pos];
      PMEDouble forceY    = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ    = cSim.pForceZAccumulator[pos];
      PMEDouble NBForceX  = cSim.pNBForceXAccumulator[pos];
      PMEDouble NBForceY  = cSim.pNBForceYAccumulator[pos];
      PMEDouble NBForceZ  = cSim.pNBForceZAccumulator[pos];
      PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
      PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
      PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
      PMEDouble GaMDTINBForceX2 = cSim.pGaMDTINBForceX[1][pos];
      PMEDouble GaMDTINBForceY2 = cSim.pGaMDTINBForceY[1][pos];
      PMEDouble GaMDTINBForceZ2 = cSim.pGaMDTINBForceZ[1][pos];
      forceX              = forceX * fwgtd;
      forceY              = forceY * fwgtd;
      forceZ              = forceZ * fwgtd;
      NBForceX   = NBForceX-(NBForceX-GaMDTINBForceX-GaMDTINBForceX2) *(1-fwgtd)-GaMDTINBForceX*(1-fwgt)-GaMDTINBForceX2 * (1-fwgtb);
      NBForceY   = NBForceY-(NBForceY-GaMDTINBForceY-GaMDTINBForceY2) *(1-fwgtd)-GaMDTINBForceY*(1-fwgt)-GaMDTINBForceY2 * (1-fwgtb);
      NBForceZ   = NBForceZ-(NBForceZ-GaMDTINBForceZ-GaMDTINBForceZ2) *(1-fwgtd)-GaMDTINBForceZ*(1-fwgt)-GaMDTINBForceZ2 * (1-fwgtb);
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
      cSim.pNBForceXAccumulator[pos]  = NBForceX;
      cSim.pNBForceYAccumulator[pos]  = NBForceY;
      cSim.pNBForceZAccumulator[pos]  = NBForceZ;
    } else {
      PMEDouble forceX   = cSim.pForceXAccumulator[pos];
      PMEDouble forceY   = cSim.pForceYAccumulator[pos];
      PMEDouble forceZ   = cSim.pForceZAccumulator[pos];
      PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
      PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
      PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
      PMEDouble GaMDTINBForceX2 = cSim.pGaMDTINBForceX[1][pos];
      PMEDouble GaMDTINBForceY2 = cSim.pGaMDTINBForceY[1][pos];
      PMEDouble GaMDTINBForceZ2 = cSim.pGaMDTINBForceZ[1][pos];
      forceX   = forceX - (forceX-GaMDTINBForceX-GaMDTINBForceX2) * (1-fwgtd) - GaMDTINBForceX*(1-fwgt)-GaMDTINBForceX2*(1-fwgtb);
      forceY   = forceY - (forceY-GaMDTINBForceY-GaMDTINBForceY2) * (1-fwgtd) - GaMDTINBForceY*(1-fwgt)-GaMDTINBForceY2*(1-fwgtb);
      forceZ   = forceZ - (forceZ-GaMDTINBForceZ-GaMDTINBForceZ2) * (1-fwgtd) - GaMDTINBForceZ*(1-fwgt)-GaMDTINBForceZ2*(1-fwgtb);
      cSim.pForceXAccumulator[pos]    = forceX;
      cSim.pForceYAccumulator[pos]    = forceY;
      cSim.pForceZAccumulator[pos]    = forceZ;
    }
    pos0 += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kCalcGAMDWeightAndScaleFrc_sc_triple4_kernel: see the host function below for a better description
//                                       of what this does.
//
// Arguments:
// works only for igamd eq. 24 and 25
//   pot_ene_tot: total potential energy (passed in from gpu_calculate_and_apply_gamd_weights_
//                in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
__global__ void
__launch_bounds__(GENERAL_THREADS_PER_BLOCK, 1)
kCalcGAMDWeightAndScaleFrc_sc_triple4_kernel(PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgt, PMEDouble fwgtd, PMEDouble fwgtb,bool useImage)
{
 // Calculate GAMD weight, setting dihedral boost (tboost) to zero for now
  unsigned int pos0 = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  // scale atomic force in SC region
   while (pos0 < cSim.atoms) {
       unsigned int pos = ( useImage ) ? cSim.pImageAtomLookup[pos0] : pos0;
       if ((cSim.ntp > 0) && (cSim.barostat == 1)) {
        PMEDouble forceX    = cSim.pForceXAccumulator[pos];
        PMEDouble forceY    = cSim.pForceYAccumulator[pos];
        PMEDouble forceZ    = cSim.pForceZAccumulator[pos];
        PMEDouble NBForceX  = cSim.pNBForceXAccumulator[pos];
        PMEDouble NBForceY  = cSim.pNBForceYAccumulator[pos];
        PMEDouble NBForceZ  = cSim.pNBForceZAccumulator[pos];
        PMEDouble GaMDTIForceX = cSim.pGaMDTIForceX[0][pos];
        PMEDouble GaMDTIForceY = cSim.pGaMDTIForceY[0][pos];
        PMEDouble GaMDTIForceZ = cSim.pGaMDTIForceZ[0][pos];
        PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
        PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
        PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
        PMEDouble GaMDTINBForceX2 = cSim.pGaMDTINBForceX[1][pos];
        PMEDouble GaMDTINBForceY2 = cSim.pGaMDTINBForceY[1][pos];
        PMEDouble GaMDTINBForceZ2 = cSim.pGaMDTINBForceZ[1][pos];
        forceX = forceX-(forceX-GaMDTIForceX)*(1-fwgtd)-GaMDTIForceX*(1-fwgt);
        forceY = forceY-(forceY-GaMDTIForceY)*(1-fwgtd)-GaMDTIForceY*(1-fwgt);
        forceZ = forceZ-(forceZ-GaMDTIForceZ)*(1-fwgtd)-GaMDTIForceZ*(1-fwgt);
        NBForceX = NBForceX-(NBForceX-GaMDTINBForceX-GaMDTINBForceX2)*(1-fwgtd)-GaMDTINBForceX*(1-fwgt)-GaMDTINBForceX2*(1-fwgtb);
        NBForceY = NBForceY-(NBForceY-GaMDTINBForceY-GaMDTINBForceY2)*(1-fwgtd)-GaMDTINBForceY*(1-fwgt)-GaMDTINBForceY2*(1-fwgtb);
        NBForceZ = NBForceZ-(NBForceZ-GaMDTINBForceZ-GaMDTINBForceZ2)*(1-fwgtd)-GaMDTINBForceZ*(1-fwgt)-GaMDTINBForceZ2*(1-fwgtb);
        cSim.pForceXAccumulator[pos]    = forceX;
        cSim.pForceYAccumulator[pos]    = forceY;
        cSim.pForceZAccumulator[pos]    = forceZ;
        cSim.pNBForceXAccumulator[pos]  = NBForceX;
        cSim.pNBForceYAccumulator[pos]  = NBForceY;
        cSim.pNBForceZAccumulator[pos]  = NBForceZ;
    }
    pos0 += increment;
  }
}
//---------------------------------------------------------------------------------------------
// kCalcGAMDWeightAndScaleFrc_sc_triple4_kernel: see the host function below for a better description
//                                       of what this does.
//
// Arguments:
// works only for igamd eq. 27
//   pot_ene_tot: total potential energy (passed in from gpu_calculate_and_apply_gamd_weights_
//                in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
__global__ void
__launch_bounds__(GENERAL_THREADS_PER_BLOCK, 1)
kCalcGAMDWeightAndScaleFrc_sc_triple5_kernel(PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgt, PMEDouble fwgtd, PMEDouble fwgtb,bool useImage)
{
 // Calculate GAMD weight, setting dihedral boost (tboost) to zero for now
  unsigned int pos0 = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  // scale atomic force in SC region
   while (pos0 < cSim.atoms) {
       unsigned int pos = ( useImage ) ? cSim.pImageAtomLookup[pos0] : pos0;
       if ((cSim.ntp > 0) && (cSim.barostat == 1)) {
        PMEDouble forceX    = cSim.pForceXAccumulator[pos];
        PMEDouble forceY    = cSim.pForceYAccumulator[pos];
        PMEDouble forceZ    = cSim.pForceZAccumulator[pos];
        PMEDouble NBForceX  = cSim.pNBForceXAccumulator[pos];
        PMEDouble NBForceY  = cSim.pNBForceYAccumulator[pos];
        PMEDouble NBForceZ  = cSim.pNBForceZAccumulator[pos];
        PMEDouble GaMDTIForceX = cSim.pGaMDTIForceX[0][pos];
        PMEDouble GaMDTIForceY = cSim.pGaMDTIForceY[0][pos];
        PMEDouble GaMDTIForceZ = cSim.pGaMDTIForceZ[0][pos];
        PMEDouble GaMDTINBForceX = cSim.pGaMDTINBForceX[0][pos];
        PMEDouble GaMDTINBForceY = cSim.pGaMDTINBForceY[0][pos];
        PMEDouble GaMDTINBForceZ = cSim.pGaMDTINBForceZ[0][pos];
        forceX = forceX-(forceX-GaMDTIForceX)*(1-fwgtb)-GaMDTIForceX*(1-fwgt);
        forceY = forceY-(forceY-GaMDTIForceY)*(1-fwgtb)-GaMDTIForceY*(1-fwgt);
        forceZ = forceZ-(forceZ-GaMDTIForceZ)*(1-fwgtb)-GaMDTIForceZ*(1-fwgt);
        NBForceX = NBForceX-(NBForceX-GaMDTINBForceX)*(1-fwgtd)-GaMDTINBForceX*(1-fwgt);
        NBForceY = NBForceY-(NBForceY-GaMDTINBForceY)*(1-fwgtd)-GaMDTINBForceY*(1-fwgt);
        NBForceZ = NBForceZ-(NBForceZ-GaMDTINBForceZ)*(1-fwgtd)-GaMDTINBForceZ*(1-fwgt);
        cSim.pForceXAccumulator[pos]    = forceX;
        cSim.pForceYAccumulator[pos]    = forceY;
        cSim.pForceZAccumulator[pos]    = forceZ;
        cSim.pNBForceXAccumulator[pos]  = NBForceX;
        cSim.pNBForceYAccumulator[pos]  = NBForceY;
        cSim.pNBForceZAccumulator[pos]  = NBForceZ;
    }
    pos0 += increment;
  }
}


//---------------------------------------------------------------------------------------------
// kCalculateGAMDWeightAndScaleForces: launch the kernel to compute the weights for Gaussian
//                                     accelerated MD and scale the forces appropriately.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, the energy
//                function, and the GaMD settings
//   pot_ene_tot: total potential energy (passed in from the calling function
//                                        gpu_calculate_and_apply_gamd_weights_ in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
void kCalculateGAMDWeightAndScaleForces(gpuContext gpu, PMEDouble pot_ene_tot,
                                        PMEDouble dih_ene_tot, PMEDouble fwgt)
{
  // Get the tile size into shorter variable names: this is extremely picky but
  // I DO NOT like line wrapping, and these names may help to show that the tile
  // size for GAMD weighting and force scaling is indeed different than the
  // tile size for the local force computations themselves under GaMD.
  int nBlocks = gpu->blocks;
  int nGenThreads = gpu->generalThreadsPerBlock;

  kCalcGAMDWeightAndScaleFrc_kernel<<<nBlocks, nGenThreads>>>(pot_ene_tot, dih_ene_tot, fwgt);
  LAUNCHERROR("kCalculateGAMDWeightAndScaleForces");
}

//---------------------------------------------------------------------------------------------
// kCalculateGAMDWeightAndScaleForces_nb: launch the kernel to compute the weights for Gaussian
//                                     accelerated MD and scale the forces appropriately.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, the energy
//                function, and the GaMD settings
//   pot_ene_tot: total potential energy (passed in from the calling function
//                                        gpu_calculate_and_apply_gamd_weights_ in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
void kCalculateGAMDWeightAndScaleForces_nb(gpuContext gpu, PMEDouble pot_ene_nb,
                                           PMEDouble dih_ene_tot, PMEDouble fwgt)
{
  // Get the tile size into shorter variable names: this is extremely picky but
  // I DO NOT like line wrapping, and these names may help to show that the tile
  // size for GAMD weighting and force scaling is indeed different than the
  // tile size for the local force computations themselves under GaMD.
  int nBlocks = gpu->blocks;
  int nGenThreads = gpu->generalThreadsPerBlock;

  kCalcGAMDWeightAndScaleFrc_nb_kernel<<<nBlocks, nGenThreads>>>(pot_ene_nb, dih_ene_tot,
                                                                 fwgt);
  LAUNCHERROR("kCalculateGAMDWeightAndScaleForces_nb");
}

//---------------------------------------------------------------------------------------------
// kCalculateGAMDWeightAndScaleForces_gb: launch the kernel to compute the weights for Gaussian
//                                     accelerated MD and scale the forces appropriately.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, the energy
//                function, and the GaMD settings
//   pot_ene_tot: total potential energy (passed in from the calling function
//                                        gpu_calculate_and_apply_gamd_weights_ in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
void kCalculateGAMDWeightAndScaleForces_gb(gpuContext gpu, PMEDouble pot_ene_tot,
                                        PMEDouble dih_ene_tot, PMEDouble fwgt)
{
  // Get the tile size into shorter variable names: this is extremely picky but
  // I DO NOT like line wrapping, and these names may help to show that the tile
  // size for GAMD weighting and force scaling is indeed different than the
  // tile size for the local force computations themselves under GaMD.
  int nBlocks = gpu->blocks;
  int nGenThreads = gpu->generalThreadsPerBlock;
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));

  kCalcGAMDWeightAndScaleFrc_gb_kernel<<<nBlocks, nGenThreads>>>(pot_ene_tot, dih_ene_tot, fwgt, useImage);
  LAUNCHERROR("kCalculateGAMDWeightAndScaleForces_gb");
}

//---------------------------------------------------------------------------------------------
// NOTE: kCalculateGAMDWeightAndScaleForces_gb_nb should not be used for the moment
// as pNBForce is not saved separately with GB.
//
// kCalculateGAMDWeightAndScaleForces_gb_nb: launch the kernel to compute the weights for Gaussian
//                                     accelerated MD and scale the forces appropriately.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, the energy
//                function, and the GaMD settings
//   pot_ene_tot: total potential energy (passed in from the calling function
//                                        gpu_calculate_and_apply_gamd_weights_ in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
void kCalculateGAMDWeightAndScaleForces_gb_nb(gpuContext gpu, PMEDouble pot_ene_nb,
                                           PMEDouble dih_ene_tot, PMEDouble fwgt)
{
  // Get the tile size into shorter variable names: this is extremely picky but
  // I DO NOT like line wrapping, and these names may help to show that the tile
  // size for GAMD weighting and force scaling is indeed different than the
  // tile size for the local force computations themselves under GaMD.
  int nBlocks = gpu->blocks;
  int nGenThreads = gpu->generalThreadsPerBlock;
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));

  kCalcGAMDWeightAndScaleFrc_gb_nb_kernel<<<nBlocks, nGenThreads>>>(pot_ene_nb, dih_ene_tot,
                                                                 fwgt, useImage);
  LAUNCHERROR("kCalculateGAMDWeightAndScaleForces_gb_nb");
}

//---------------------------------------------------------------------------------------------
// kCalculateGAMDWeightAndScaleForces_ti_others: launch the kernel to compute the weights for Gaussian
//                                     accelerated MD and scale the forces appropriately.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, the energy
//                function, and the GaMD settings
//   pot_ene_tot: total potential energy (passed in from the calling function
//                                        gpu_calculate_and_apply_gamd_weights_ in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
void kCalculateGAMDWeightAndScaleForces_ti_others(gpuContext gpu, PMEDouble pot_ene_ti_region,
                                           PMEDouble dih_ene_tot, PMEDouble fwgtd)
{
  // Get the tile size into shorter variable names: this is extremely picky but
  // I DO NOT like line wrapping, and these names may help to show that the tile
  // size for GAMD weighting and force scaling is indeed different than the
  // tile size for the local force computations themselves under GaMD.
  int nBlocks = gpu->blocks;
  int nGenThreads = gpu->generalThreadsPerBlock;
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));

  kCalcGAMDWeightAndScaleFrc_ti_others_kernel<<<nBlocks, nGenThreads>>>(pot_ene_ti_region, dih_ene_tot, fwgtd, useImage);
  LAUNCHERROR("kCalculateGAMDWeightAndScaleForces_ti_others");
}

//---------------------------------------------------------------------------------------------
// kCalculateGAMDWeightAndScaleForces_bond: launch the kernel to compute the weights for Gaussian
//                                     accelerated MD and scale the forces appropriately.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, the energy
//                function, and the GaMD settings
//   pot_ene_tot: total potential energy (passed in from the calling function
//                                        gpu_calculate_and_apply_gamd_weights_ in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
void kCalculateGAMDWeightAndScaleForces_bond(gpuContext gpu, PMEDouble pot_ene_nb,
                                           PMEDouble dih_ene_tot, PMEDouble fwgt)
{
  // Get the tile size into shorter variable names: this is extremely picky but
  // I DO NOT like line wrapping, and these names may help to show that the tile
  // size for GAMD weighting and force scaling is indeed different than the
  // tile size for the local force computations themselves under GaMD.
  int nBlocks = gpu->blocks;
  int nGenThreads = gpu->generalThreadsPerBlock;
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));


  kCalcGAMDWeightAndScaleFrc_bond_kernel<<<nBlocks, nGenThreads>>>(pot_ene_nb, dih_ene_tot,
                                                                 fwgt,useImage);
  LAUNCHERROR("kCalculateGAMDWeightAndScaleForces_bond");
}

//---------------------------------------------------------------------------------------------
// kCalculateGAMDWeightAndScaleForces_ti_bonded_others: launch the kernel to compute the weights for Gaussian
//                                     accelerated MD and scale the forces appropriately.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, the energy
//                function, and the GaMD settings
//   pot_ene_tot: total potential energy (passed in from the calling function
//                                        gpu_calculate_and_apply_gamd_weights_ in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
void kCalculateGAMDWeightAndScaleForces_ti_bonded_others(gpuContext gpu, PMEDouble pot_ene_ti_region,
                                           PMEDouble dih_ene_tot, PMEDouble fwgtd)
{
  // Get the tile size into shorter variable names: this is extremely picky but
  // I DO NOT like line wrapping, and these names may help to show that the tile
  // size for GAMD weighting and force scaling is indeed different than the
  // tile size for the local force computations themselves under GaMD.
  int nBlocks = gpu->blocks;
  int nGenThreads = gpu->generalThreadsPerBlock;
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));

  kCalcGAMDWeightAndScaleFrc_ti_bonded_others_kernel<<<nBlocks, nGenThreads>>>(pot_ene_ti_region, dih_ene_tot, fwgtd, useImage);
  LAUNCHERROR("kCalculateGAMDWeightAndScaleForces_ti_bonded_others");
}

//---------------------------------------------------------------------------------------------
// kCalculateGAMDWeightAndScaleForces_ti_nonbonded_others: launch the kernel to compute the weights for Gaussian
//                                     accelerated MD and scale the forces appropriately.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, the energy
//                function, and the GaMD settings
//   pot_ene_tot: total potential energy (passed in from the calling function
//                                        gpu_calculate_and_apply_gamd_weights_ in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
void kCalculateGAMDWeightAndScaleForces_ti_nonbonded_others(gpuContext gpu, PMEDouble pot_ene_ti_region,
                                           PMEDouble dih_ene_tot, PMEDouble fwgtd)
{
  // Get the tile size into shorter variable names: this is extremely picky but
  // I DO NOT like line wrapping, and these names may help to show that the tile
  // size for GAMD weighting and force scaling is indeed different than the
  // tile size for the local force computations themselves under GaMD.
  int nBlocks = gpu->blocks;
  int nGenThreads = gpu->generalThreadsPerBlock;
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));

  kCalcGAMDWeightAndScaleFrc_ti_nonbonded_others_kernel<<<nBlocks, nGenThreads>>>(pot_ene_ti_region, dih_ene_tot, fwgtd, useImage);
  LAUNCHERROR("kCalculateGAMDWeightAndScaleForces_ti_nonbonded_others");
}

//---------------------------------------------------------------------------------------------
// kCalculateGAMDWeightAndScaleForces_sc_dual_nonbonded_ppi: launch the kernel to compute the weights for Gaussian
//                                     accelerated MD and scale the forces appropriately.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, the energy
//                function, and the GaMD settings
//   pot_ene_tot: total potential energy (passed in from the calling function
//                                        gpu_calculate_and_apply_gamd_weights_ in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
void kCalculateGAMDWeightAndScaleForces_sc_dual_nonbonded_ppi(gpuContext gpu, PMEDouble pot_ene_ti_region,
                                           PMEDouble dih_ene_tot, PMEDouble fwgt, PMEDouble fwgtd)
{
  // Get the tile size into shorter variable names: this is extremely picky but
  // I DO NOT like line wrapping, and these names may help to show that the tile
  // size for GAMD weighting and force scaling is indeed different than the
  // tile size for the local force computations themselves under GaMD.
  int nBlocks = gpu->blocks;
  int nGenThreads = gpu->generalThreadsPerBlock;
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));

  kCalcGAMDWeightAndScaleFrc_sc_dual_nonbonded_kernel<<<nBlocks, nGenThreads>>>(pot_ene_ti_region, dih_ene_tot,fwgt, fwgtd, useImage);
  LAUNCHERROR("kCalculateGAMDWeightAndScaleForces_sc_dual_nonbonded_ppi");
}

//---------------------------------------------------------------------------------------------
// kCalculateGAMDWeightAndScaleForces_sc_dual_bonded_ppi: launch the kernel to compute the weights for Gaussian
//                                     accelerated MD and scale the forces appropriately.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, the energy
//                function, and the GaMD settings
//   pot_ene_tot: total potential energy (passed in from the calling function
//                                        gpu_calculate_and_apply_gamd_weights_ in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
void kCalculateGAMDWeightAndScaleForces_sc_dual_bonded_ppi(gpuContext gpu, PMEDouble pot_ene_ti_region,
                                           PMEDouble dih_ene_tot, PMEDouble fwgt, PMEDouble fwgtd)
{
  // Get the tile size into shorter variable names: this is extremely picky but
  // I DO NOT like line wrapping, and these names may help to show that the tile
  // size for GAMD weighting and force scaling is indeed different than the
  // tile size for the local force computations themselves under GaMD.
  int nBlocks = gpu->blocks;
  int nGenThreads = gpu->generalThreadsPerBlock;
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));

  kCalcGAMDWeightAndScaleFrc_sc_dual_bonded_kernel<<<nBlocks, nGenThreads>>>(pot_ene_ti_region, dih_ene_tot,fwgt, fwgtd, useImage);
  LAUNCHERROR("kCalculateGAMDWeightAndScaleForces_sc_dual_nonbonded_ppi");
}

//---------------------------------------------------------------------------------------------
// kCalculateGAMDWeightAndScaleForces_sc_triple_ppi: launch the kernel to compute the weights for Gaussian
//                                     accelerated MD and scale the forces appropriately.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, the energy
//                function, and the GaMD settings
//   pot_ene_tot: total potential energy (passed in from the calling function
//                                        gpu_calculate_and_apply_gamd_weights_ in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
void kCalculateGAMDWeightAndScaleForces_sc_triple_ppi(gpuContext gpu, PMEDouble pot_ene_ti_region,
                                    PMEDouble dih_ene_tot, PMEDouble ppi_bond_tot,PMEDouble fwgt, PMEDouble fwgtd,PMEDouble fwgtb)
{
  // Get the tile size into shorter variable names: this is extremely picky but
  // I DO NOT like line wrapping, and these names may help to show that the tile
  // size for GAMD weighting and force scaling is indeed different than the
  // tile size for the local force computations themselves under GaMD.
  int nBlocks = gpu->blocks;
  int nGenThreads = gpu->generalThreadsPerBlock;
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));

  kCalcGAMDWeightAndScaleFrc_sc_triple_kernel<<<nBlocks, nGenThreads>>>(pot_ene_ti_region, dih_ene_tot,fwgt, fwgtd, fwgtb,useImage);
  LAUNCHERROR("kCalculateGAMDWeightAndScaleForces_sc_dual_nonbonded_ppi");
}

//---------------------------------------------------------------------------------------------
// kCalculateGAMDWeightAndScaleForces_sc_triple2_ppi: launch the kernel to compute the weights for Gaussian
//                                     accelerated MD and scale the forces appropriately.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, the energy
//                function, and the GaMD settings
//   pot_ene_tot: total potential energy (passed in from the calling function
//                                        gpu_calculate_and_apply_gamd_weights_ in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
void kCalculateGAMDWeightAndScaleForces_sc_triple2_ppi(gpuContext gpu, PMEDouble pot_ene_ti_region,
                                    PMEDouble dih_ene_tot, PMEDouble ppi_bond_tot,PMEDouble fwgt, PMEDouble fwgtd,PMEDouble fwgtb)
{
  // Get the tile size into shorter variable names: this is extremely picky but
  // I DO NOT like line wrapping, and these names may help to show that the tile
  // size for GAMD weighting and force scaling is indeed different than the
  // tile size for the local force computations themselves under GaMD.
  int nBlocks = gpu->blocks;
  int nGenThreads = gpu->generalThreadsPerBlock;
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));

  kCalcGAMDWeightAndScaleFrc_sc_triple2_kernel<<<nBlocks, nGenThreads>>>(pot_ene_ti_region, dih_ene_tot,fwgt, fwgtd, fwgtb,useImage);
  LAUNCHERROR("kCalculateGAMDWeightAndScaleForces_sc_triple2_nonbonded_ppi");
}

//---------------------------------------------------------------------------------------------
// kCalculateGAMDWeightAndScaleForces_sc_triple3_ppi: launch the kernel to compute the weights for Gaussian
//                                     accelerated MD and scale the forces appropriately.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, the energy
//                function, and the GaMD settings
//   pot_ene_tot: total potential energy (passed in from the calling function
//                                        gpu_calculate_and_apply_gamd_weights_ in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
void kCalculateGAMDWeightAndScaleForces_sc_triple3_ppi(gpuContext gpu, PMEDouble pot_ene_ti_region,
		                                    PMEDouble dih_ene_tot, PMEDouble ppi_bond_tot,PMEDouble fwgt, PMEDouble fwgtd,PMEDouble fwgtb)
{
  // Get the tile size into shorter variable names: this is extremely picky but
  // I DO NOT like line wrapping, and these names may help to show that the tile
  // size for GAMD weighting and force scaling is indeed different than the
  // tile size for the local force computations themselves under GaMD.
  int nBlocks = gpu->blocks;
  int nGenThreads = gpu->generalThreadsPerBlock;
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));

  kCalcGAMDWeightAndScaleFrc_sc_triple3_kernel<<<nBlocks, nGenThreads>>>(pot_ene_ti_region, dih_ene_tot,fwgt, fwgtd, fwgtb,useImage);
  LAUNCHERROR("kCalculateGAMDWeightAndScaleForces_sc_dual_nonbonded_ppi");
}
//---------------------------------------------------------------------------------------------
// kCalculateGAMDWeightAndScaleForces_sc_triple4_ppi: launch the kernel to compute the weights for Gaussian
//                                     accelerated MD and scale the forces appropriately.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, the energy
//                function, and the GaMD settings
//   pot_ene_tot: total potential energy (passed in from the calling function
//                                        gpu_calculate_and_apply_gamd_weights_ in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
void kCalculateGAMDWeightAndScaleForces_sc_triple4_ppi(gpuContext gpu, PMEDouble pot_ene_ti_region,
		                                    PMEDouble dih_ene_tot, PMEDouble ppi_bond_tot,PMEDouble fwgt, PMEDouble fwgtd,PMEDouble fwgtb)
{
  // Get the tile size into shorter variable names: this is extremely picky but
  // I DO NOT like line wrapping, and these names may help to show that the tile
  // size for GAMD weighting and force scaling is indeed different than the
  // tile size for the local force computations themselves under GaMD.
  int nBlocks = gpu->blocks;
  int nGenThreads = gpu->generalThreadsPerBlock;
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));

  kCalcGAMDWeightAndScaleFrc_sc_triple4_kernel<<<nBlocks, nGenThreads>>>(pot_ene_ti_region, dih_ene_tot,fwgt, fwgtd, fwgtb,useImage);
  LAUNCHERROR("kCalculateGAMDWeightAndScaleForces_sc_triple4_nonbonded_ppi");
}
//---------------------------------------------------------------------------------------------
// kCalculateGAMDWeightAndScaleForces_sc_triple4_ppi: launch the kernel to compute the weights for Gaussian
//                                     accelerated MD and scale the forces appropriately.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, the energy
//                function, and the GaMD settings
//   pot_ene_tot: total potential energy (passed in from the calling function
//                                        gpu_calculate_and_apply_gamd_weights_ in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
void kCalculateGAMDWeightAndScaleForces_sc_triple5_ppi(gpuContext gpu, PMEDouble pot_ene_ti_region,
		                                    PMEDouble dih_ene_tot, PMEDouble ppi_bond_tot,PMEDouble fwgt, PMEDouble fwgtd,PMEDouble fwgtb)
{
  // Get the tile size into shorter variable names: this is extremely picky but
  // I DO NOT like line wrapping, and these names may help to show that the tile
  // size for GAMD weighting and force scaling is indeed different than the
  // tile size for the local force computations themselves under GaMD.
  int nBlocks = gpu->blocks;
  int nGenThreads = gpu->generalThreadsPerBlock;
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));

  kCalcGAMDWeightAndScaleFrc_sc_triple5_kernel<<<nBlocks, nGenThreads>>>(pot_ene_ti_region, dih_ene_tot,fwgt, fwgtd, fwgtb,useImage);
  LAUNCHERROR("kCalculateGAMDWeightAndScaleForces_sc_triple5_nonbonded_ppi");
}

//---------------------------------------------------------------------------------------------
// kCalculateGAMDWeightAndScaleForces_ti: launch the kernel to compute the weights for Gaussian
//                                     accelerated MD and scale the forces appropriately.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, the energy
//                function, and the GaMD settings
//   pot_ene_tot: total potential energy (passed in from the calling function
//                                        gpu_calculate_and_apply_gamd_weights_ in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
void kCalculateGAMDWeightAndScaleForces_ti(gpuContext gpu, PMEDouble pot_ene_ti_region,
                                           PMEDouble dih_ene_tot, PMEDouble fwgt)
{
  // Get the tile size into shorter variable names: this is extremely picky but
  // I DO NOT like line wrapping, and these names may help to show that the tile
  // size for GAMD weighting and force scaling is indeed different than the
  // tile size for the local force computations themselves under GaMD.
  int nBlocks = gpu->blocks;
  int nGenThreads = gpu->generalThreadsPerBlock;
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));

  kCalcGAMDWeightAndScaleFrc_ti_kernel<<<nBlocks, nGenThreads>>>(pot_ene_ti_region, dih_ene_tot, fwgt, useImage);
  LAUNCHERROR("kCalculateGAMDWeightAndScaleForces_ti");
}

//---------------------------------------------------------------------------------------------
// kCalculateGAMDWeightAndScaleForces_ti_bonded: launch the kernel to compute the weights for Gaussian
//                                     accelerated MD and scale the forces appropriately.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, the energy
//                function, and the GaMD settings
//   pot_ene_tot: total potential energy (passed in from the calling function
//                                        gpu_calculate_and_apply_gamd_weights_ in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
void kCalculateGAMDWeightAndScaleForces_ti_bonded(gpuContext gpu, PMEDouble pot_ene_ti_region,
                                           PMEDouble dih_ene_tot, PMEDouble fwgt)
{
  // Get the tile size into shorter variable names: this is extremely picky but
  // I DO NOT like line wrapping, and these names may help to show that the tile
  // size for GAMD weighting and force scaling is indeed different than the
  // tile size for the local force computations themselves under GaMD.
  int nBlocks = gpu->blocks;
  int nGenThreads = gpu->generalThreadsPerBlock;
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));

  kCalcGAMDWeightAndScaleFrc_ti_bonded_kernel<<<nBlocks, nGenThreads>>>(pot_ene_ti_region, dih_ene_tot, fwgt, useImage);
  LAUNCHERROR("kCalculateGAMDWeightAndScaleForces_ti_bonded");
}

//---------------------------------------------------------------------------------------------
// kCalculateGAMDWeightAndScaleForces_sc_ppi: launch the kernel to compute the weights for Gaussian
//                                     accelerated MD and scale the forces appropriately.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, the energy
//                function, and the GaMD settings
//   pot_ene_tot: total potential energy (passed in from the calling function
//                                        gpu_calculate_and_apply_gamd_weights_ in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
void kCalculateGAMDWeightAndScaleForces_sc_ppi(gpuContext gpu, PMEDouble pot_ene_ti_region,
                                           PMEDouble dih_ene_tot, PMEDouble fwgt)
{
  // Get the tile size into shorter variable names: this is extremely picky but
  // I DO NOT like line wrapping, and these names may help to show that the tile
  // size for GAMD weighting and force scaling is indeed different than the
  // tile size for the local force computations themselves under GaMD.
  int nBlocks = gpu->blocks;
  int nGenThreads = gpu->generalThreadsPerBlock;
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));

  kCalcGAMDWeightAndScaleFrc_sc_ppi_kernel<<<nBlocks, nGenThreads>>>(pot_ene_ti_region, dih_ene_tot, fwgt, useImage);
  LAUNCHERROR("kCalculateGAMDWeightAndScaleForces_sc_ppi");
}
//---------------------------------------------------------------------------------------------
// kCalculateGAMDWeightAndScaleForces_sc_ppi: launch the kernel to compute the weights for Gaussian
//                                     accelerated MD and scale the forces appropriately.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, the energy
//                function, and the GaMD settings
//   pot_ene_tot: total potential energy (passed in from the calling function
//                                        gpu_calculate_and_apply_gamd_weights_ in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
void kCalculateGAMDWeightAndScaleForces_sc_others_ppi(gpuContext gpu, PMEDouble pot_ene_ti_region,
                                           PMEDouble dih_ene_tot, PMEDouble fwgtd)
{
  // Get the tile size into shorter variable names: this is extremely picky but
  // I DO NOT like line wrapping, and these names may help to show that the tile
  // size for GAMD weighting and force scaling is indeed different than the
  // tile size for the local force computations themselves under GaMD.
  int nBlocks = gpu->blocks;
  int nGenThreads = gpu->generalThreadsPerBlock;
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));

  kCalcGAMDWeightAndScaleFrc_sc_others_ppi_kernel<<<nBlocks, nGenThreads>>>(pot_ene_ti_region, dih_ene_tot, fwgtd, useImage);
  LAUNCHERROR("kCalculateGAMDWeightAndScaleForces_sc");
}

//---------------------------------------------------------------------------------------------
// kCalculateGAMDWeightAndScaleForces_sc_bonded_ppi: launch the kernel to compute the weights for Gaussian
//                                     accelerated MD and scale the forces appropriately.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, the energy
//                function, and the GaMD settings
//   pot_ene_tot: total potential energy (passed in from the calling function
//                                        gpu_calculate_and_apply_gamd_weights_ in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
void kCalculateGAMDWeightAndScaleForces_sc_bonded_ppi(gpuContext gpu, PMEDouble pot_ene_ti_region,
                                           PMEDouble dih_ene_tot, PMEDouble fwgt)
{
  // Get the tile size into shorter variable names: this is extremely picky but
  // I DO NOT like line wrapping, and these names may help to show that the tile
  // size for GAMD weighting and force scaling is indeed different than the
  // tile size for the local force computations themselves under GaMD.
  int nBlocks = gpu->blocks;
  int nGenThreads = gpu->generalThreadsPerBlock;
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));

  kCalcGAMDWeightAndScaleFrc_sc_bonded_ppi_kernel<<<nBlocks, nGenThreads>>>(pot_ene_ti_region, dih_ene_tot, fwgt, useImage);
  LAUNCHERROR("kCalculateGAMDWeightAndScaleForces_sc_bonded");
}

//---------------------------------------------------------------------------------------------
// kCalculateGAMDWeightAndScaleForces_sc_nonbonded: launch the kernel to compute the weights for Gaussian
//                                     accelerated MD and scale the forces appropriately.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, the energy
//                function, and the GaMD settings
//   pot_ene_tot: total potential energy (passed in from the calling function
//                                        gpu_calculate_and_apply_gamd_weights_ in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
void kCalculateGAMDWeightAndScaleForces_sc_nonbonded(gpuContext gpu, PMEDouble pot_ene_ti_region,
                                           PMEDouble dih_ene_tot, PMEDouble fwgt)
{
  // Get the tile size into shorter variable names: this is extremely picky but
  // I DO NOT like line wrapping, and these names may help to show that the tile
  // size for GAMD weighting and force scaling is indeed different than the
  // tile size for the local force computations themselves under GaMD.
  int nBlocks = gpu->blocks;
  int nGenThreads = gpu->generalThreadsPerBlock;
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));

  kCalcGAMDWeightAndScaleFrc_sc_nonbonded_kernel<<<nBlocks, nGenThreads>>>(pot_ene_ti_region, dih_ene_tot, fwgt, useImage);
  LAUNCHERROR("kCalculateGAMDWeightAndScaleForces_sc_nonbonded");
}

//---------------------------------------------------------------------------------------------
// kCalculateGAMDWeightAndScaleForces_ti_nonbonded: launch the kernel to compute the weights for Gaussian
//                                     accelerated MD and scale the forces appropriately.
//
// Arguments:
//   gpu:         overarching type for storing all parameters, coordinates, the energy
//                function, and the GaMD settings
//   pot_ene_tot: total potential energy (passed in from the calling function
//                                        gpu_calculate_and_apply_gamd_weights_ in gpu.cpp)
//   dih_ene_tot: total dihedral energy (passed in like pot_ene_tot)
//   fwgt:        set inside the calling function and passed in, relates to the GaMD weight
//---------------------------------------------------------------------------------------------
void kCalculateGAMDWeightAndScaleForces_ti_nonbonded(gpuContext gpu, PMEDouble pot_ene_ti_region,
                                           PMEDouble dih_ene_tot, PMEDouble fwgt)
{
  // Get the tile size into shorter variable names: this is extremely picky but
  // I DO NOT like line wrapping, and these names may help to show that the tile
  // size for GAMD weighting and force scaling is indeed different than the
  // tile size for the local force computations themselves under GaMD.
  int nBlocks = gpu->blocks;
  int nGenThreads = gpu->generalThreadsPerBlock;
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));

  kCalcGAMDWeightAndScaleFrc_ti_nonbonded_kernel<<<nBlocks, nGenThreads>>>(pot_ene_ti_region, dih_ene_tot, fwgt, useImage);
  LAUNCHERROR("kCalculateGAMDWeightAndScaleForces_ti_nonbonded");
}

//---------------------------------------------------------------------------------------------
// kCalculateLocalForcesInitKernels: initialize the kernels for local forces.  This is
//                                   analogous to kCalculatePMENonbondEnergyInitKernels (see
//                                   kCalculatePMENonbondEnergy.cu), but here there are two
//                                   things about the __shared__ memory cache to set, not one.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kCalculateLocalForcesInitKernels(gpuContext gpu)
{
  // Set the cache configuration for kernels handling special NMR interactions
  cudaFuncSetCacheConfig(kCalcNMRFrc_kernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(kCalcPMENMRFrc_kernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(kCalcNMRFrcNoDPTex_kernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(kCalcPMENMRFrcNoDPTex_kernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(kCalcNMRNrg_kernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(kCalcPMENMRNrg_kernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(kCalcNMRNrgNoDPTex_kernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(kCalcPMENMRNrgNoDPTex_kernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(kCalcNMRNrgAFE_kernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(kCalcPMENMRNrgAFE_kernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(kCalcNMRNrgNoDPTexAFE_kernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(kCalcPMENMRNrgNoDPTexAFE_kernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(kCalcNMRR6avFrc_kernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(kCalcPMENMRR6avFrc_kernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(kCalcNMRR6avFrcNoDPTex_kernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(kCalcPMENMRR6avFrcNoDPTex_kernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(kCalcNMRR6avNrg_kernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(kCalcPMENMRR6avNrg_kernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(kCalcNMRR6avNrgNoDPTex_kernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(kCalcPMENMRR6avNrgNoDPTex_kernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(kCalcNMRR6avNrgAFE_kernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(kCalcPMENMRR6avNrgAFE_kernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(kCalcNMRR6avNrgNoDPTexAFE_kernel, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(kCalcPMENMRR6avNrgNoDPTexAFE_kernel, cudaFuncCachePreferL1);

  // Set the __shared__ memory bank size to four bytes
  cudaFuncSetSharedMemConfig(kExecBondWorkUnits_kernel, cudaSharedMemBankSizeFourByte);
  cudaFuncSetSharedMemConfig(kExecBondWorkUnitsNrg_kernel, cudaSharedMemBankSizeFourByte);
  cudaFuncSetSharedMemConfig(kExecBondWorkUnitsVir_kernel, cudaSharedMemBankSizeFourByte);
  cudaFuncSetSharedMemConfig(kExecBondWorkUnitsNrgVir_kernel, cudaSharedMemBankSizeFourByte);
  cudaFuncSetSharedMemConfig(kExecBondWorkUnitsNL_kernel, cudaSharedMemBankSizeFourByte);
  cudaFuncSetSharedMemConfig(kExecBondWorkUnitsNrgNL_kernel, cudaSharedMemBankSizeFourByte);
  cudaFuncSetSharedMemConfig(kExecBondWorkUnitsVirNL_kernel, cudaSharedMemBankSizeFourByte);
  cudaFuncSetSharedMemConfig(kExecBondWorkUnitsNrgVirNL_kernel, cudaSharedMemBankSizeFourByte);

  // NMR special-purpose interactions: set the __shared__ memory bank size to eight bytes
  cudaFuncSetSharedMemConfig(kCalcNMRFrc_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalcPMENMRFrc_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalcNMRFrcNoDPTex_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalcPMENMRFrcNoDPTex_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalcNMRNrg_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalcPMENMRNrg_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalcNMRNrgNoDPTex_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalcPMENMRNrgNoDPTex_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalcNMRNrgAFE_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalcPMENMRNrgAFE_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalcNMRNrgNoDPTexAFE_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalcPMENMRNrgNoDPTexAFE_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalcNMRR6avFrc_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalcPMENMRR6avFrc_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalcNMRR6avFrcNoDPTex_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalcPMENMRR6avFrcNoDPTex_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalcNMRR6avNrg_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalcPMENMRR6avNrg_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalcNMRR6avNrgNoDPTex_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalcPMENMRR6avNrgNoDPTex_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalcNMRR6avNrgAFE_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalcPMENMRR6avNrgAFE_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalcNMRR6avNrgNoDPTexAFE_kernel, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kCalcPMENMRR6avNrgNoDPTexAFE_kernel, cudaSharedMemBankSizeEightByte);
}
