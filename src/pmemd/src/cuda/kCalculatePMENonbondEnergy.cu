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
#define GPU_CPP
#include "gputypes.h"
#include "gpu.h"
#undef GPU_CPP
#include "ptxmacros.h"
//#include "cuda_profiler_api.h"

//#define PME_VIRIAL
//#define PME_ENERGY

// Use global instance instead of a local copy
#include "simulationConst.h"
CSIM_STO simulationConst cSim;

#if !defined(__HIPCC_RDC__)

//---------------------------------------------------------------------------------------------
// SetkCalculatePMENonbondEnergySim: upload details of the PME non-bonded calculation to the
//                                   device.  This is called by the frequently invoked
//                                   gpuCopyConstants() function in gpu.cpp.
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void SetkCalculatePMENonbondEnergySim(gpuContext gpu)
{
  cudaError_t status;
  status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(simulationConst));
  RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

//---------------------------------------------------------------------------------------------
// GetkCalculatePMENonbondEnergySim: download details of the PME non-bonded calculation from
//                                   the device.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
void GetkCalculatePMENonBondEnergySim(gpuContext gpu)
{
  cudaError_t status;
  status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(simulationConst));
  RTERROR(status, "cudaMemcpyFromSymbol: SetSim copy to cSim failed");
}

#endif

#ifndef use_DPFP
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
__device__ __constant__ cudaERFC cERFC;

//---------------------------------------------------------------------------------------------
// SetkCalculatePMENonbondEnergyERFC: set constants for the fasterfc() function on the device.
//
// Arguments:
//   ewcoeff:    the Ewald coefficient, 1/(2*Gaussian sigma)
//---------------------------------------------------------------------------------------------
void SetkCalculatePMENonbondEnergyERFC(double ewcoeff)
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
  status = cudaMemcpyToSymbol(cERFC, &c, sizeof(cudaERFC));
  RTERROR(status, "cudaMemcpyToSymbol: SetERFC copy to cERFC failed");
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
  t = cERFC.c0;
  t = __internal_fmad(t, x, cERFC.c1);
  t = __internal_fmad(t, x, cERFC.c2);
  t = __internal_fmad(t, x, cERFC.c3);
  t = __internal_fmad(t, x, cERFC.c4);
  t = __internal_fmad(t, x, cERFC.c5);
  t = __internal_fmad(t, x, cERFC.c6);
  t = __internal_fmad(t, x, cERFC.c7);
  t = __internal_fmad(t, x, cERFC.c8);
  t = __internal_fmad(t, x, cERFC.c9);
  t = t * x;
  return exp2f(t);
}

//---------------------------------------------------------------------------------------------
// dOptimizeCoeffEntry: optimize an entry in the coefficients table with respect to 32 trial
//                      values over its applicable range.
//---------------------------------------------------------------------------------------------
__device__ PMEFloat4 dOptimizeCoeffEntry(double target, PMEFloat r2, PMEFloat4 basecoef)
{
  // Loop over perturbations of the coefficients
  PMEFloat4 best;
  double minr = 1.0e10;
  int i, j, k;
  for (i = -7; i <= 7; i++) {
    for (j = -7; j <= 7; j++) {
      for (k = -7; k <= 7; k++) {

        // Perturb the coefficients
        PMEFloat4 pcoef = basecoef;
        pcoef.x = __uint_as_float(__float_as_uint(pcoef.x) + i);
        pcoef.z = __uint_as_float(__float_as_uint(pcoef.z) + j);
        pcoef.w = __uint_as_float(__float_as_uint(pcoef.w) + k);

        // Compute the necessary change in the constant (y) coefficient
        double remain = target - (((double)pcoef.x * r2) + ((double)pcoef.z / r2) +
                                  ((double)pcoef.w / (r2 * r2)));
        double pcy = remain;
        // Average of 32 values of one segment (warp on NVIDIA or half-warp on AMD)
        pcy += __SHFL_XOR(WARP_MASK, pcy, 16);
        pcy += __SHFL_XOR(WARP_MASK, pcy,  8);
        pcy += __SHFL_XOR(WARP_MASK, pcy,  4);
        pcy += __SHFL_XOR(WARP_MASK, pcy,  2);
        pcy += __SHFL_XOR(WARP_MASK, pcy,  1);
        pcoef.y = (PMEFloat)(pcy * 0.03125);

        // Compute the mean error
        remain -= (double)pcoef.y;
        remain *= remain;
        // Average of 32 values of one segment (warp on NVIDIA or half-warp on AMD)
        remain += __SHFL_XOR(WARP_MASK, remain, 16);
        remain += __SHFL_XOR(WARP_MASK, remain,  8);
        remain += __SHFL_XOR(WARP_MASK, remain,  4);
        remain += __SHFL_XOR(WARP_MASK, remain,  2);
        remain += __SHFL_XOR(WARP_MASK, remain,  1);
        remain *= 0.03125;
        if (remain < minr) {
          minr = remain;
          best = pcoef;
        }
      }
    }
  }

  return best;
}

//---------------------------------------------------------------------------------------------
// kAdjustCoeffsTable_kernel: make an adjustment to the erfc spline coefficients table based
//                            on trial and error changes to the units of least place.
//---------------------------------------------------------------------------------------------
__global__ void __LAUNCH_BOUNDS__(1024, 1) kAdjustCoeffsTable_kernel()
{
  __shared__ PMEFloat4 basecoefs[64];

  int segmentIdx = threadIdx.x / 32;
  unsigned int tgx = threadIdx.x % 32;

  // Loop over all blocks
  int bpos = blockIdx.x;
  while (bpos < 144) {

    // Skip this for the first 121 blocks
    if (bpos < 121) {
      bpos += gridDim.x;
      continue;
    }

    // Each block will deal with one of the critical segments in the list
    if (threadIdx.x < 64) {
      basecoefs[threadIdx.x] = cSim.pErfcCoeffsTable[bpos*64 + threadIdx.x];
    }
    __syncthreads();

    // Each warp now takes care of its own segment on NVIDIA or two segments on AMD;
    // one block takes care of a full set of 32 equally spaced intervals for one power of two
    // in the r2 range
    double r2 = pow(2.0, bpos - 127);
    r2 += r2 * ((double)threadIdx.x + 0.5) / 1024.0;

    // Analytic result for the target in DPFP
    double target = (erfc(cSim.ew_coeff * sqrt(r2))/sqrt(r2) +
                     (2.0 * cSim.ew_coeff / sqrt(PI_VAL)) *
                     exp(-cSim.ew_coeff * cSim.ew_coeff * r2)) / r2;
    double targetEx = target - 1.0/(r2 * sqrt(r2));
    PMEFloat4 coefs;
    coefs = dOptimizeCoeffEntry(target, (PMEFloat)r2, basecoefs[2*segmentIdx]);
    __SYNCWARP(WARP_MASK);
    if (tgx == 0) {
      basecoefs[2*segmentIdx] = coefs;
    }
    coefs = dOptimizeCoeffEntry(targetEx, (PMEFloat)r2, basecoefs[2*segmentIdx + 1]);
    __SYNCWARP(WARP_MASK);
    if (tgx == 0) {
      basecoefs[2*segmentIdx + 1] = coefs;
    }
    __syncthreads();

    // Replace the coefficients with optimized values
    if (threadIdx.x < 64) {
      cSim.pErfcCoeffsTable[bpos*64 + threadIdx.x] = basecoefs[threadIdx.x];
    }
    __syncthreads();

    // Advance the block position
    bpos += gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kAdjustCoeffsTable: adjust the coefficients table to get as close as possible to the fp64
//                     results.  Launches the kernel above.
//
// Arguments:
//   gpu:       overarching struct storing information about the simulation, including atom
//              properties and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kAdjustCoeffsTable(gpuContext gpu)
{
  kAdjustCoeffsTable_kernel<<<gpu->blocks, 1024>>>();
}
#endif

//---------------------------------------------------------------------------------------------
// Particle Mesh Ewald nonbond kernels.  There are 120 in all.  A table of contents would be
// redundant, as the names are self-explanatory and so involved that finding the proper kernel
// in the pile would be no easier than searching for matching strings in the code itself.
// However, it is useful to provide some rules as to which features the various kernels have.
//
// - There are versions of every kernel that take 8, 16, and 32 atoms per warp.
// - Every kernel has special cases for computing forces alone and forces with energies
// - Every kernel has special cases for orthorhombic (orthogonal) and general unit cells
// - With those rules in mind, there are then ten types of kernels:
//   - Plain non-bonded interactions (good for standard constant volume simulations)
//   - Interactions that include a virial computation (for constant pressure simulations)
//   - Either of the above, with force switching (FSWITCH) (that's 4 so far)
//   - Minimization kernels for all three cases above (that's 8 so far)
//   - Soft core MBAR and soft core TI (just gives plain force and energies, no virial--that's
//     the final two cases)
//
// The main feature that differentiates minimization kernels is the presence of code to clamp
// the computed forces to ensure that they don't break the integer format for SPFP, or cause
// the system to distort wildly in its first moves.
//
// The organization below doesn't perhaps do the best job of separating #defines from
// particular kernels, but the annotation above should help to navigate the ~1100 lines of
// kernel declarations.
//---------------------------------------------------------------------------------------------

#if defined(AMBER_PLATFORM_AMD)
#  define INC_NLCPNE "kNLCPNE_AMD.h"
#else
#  define INC_NLCPNE "kNLCPNE.h"
#endif

// HIP-TODO: MI50 has better performance with 4
// HIP-TODO: DPFP on MI100 crashes with 5 because of a bug with v_accvgpr_read_b32 + ds_bpermute
#if defined(AMBER_PLATFORM_AMD) && defined(use_SPFP)
#define __LAUNCH_BOUNDS_PME_FORCES__ __launch_bounds__(PMENONBONDFORCES_THREADS_PER_BLOCK, 5)
#else
#define __LAUNCH_BOUNDS_PME_FORCES__ __LAUNCH_BOUNDS__(PMENONBONDFORCES_THREADS_PER_BLOCK, PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#endif
#define __LAUNCH_BOUNDS_PME_ENERGY__ __LAUNCH_BOUNDS__(PMENONBONDENERGY_THREADS_PER_BLOCK, PMENONBONDENERGY_BLOCKS_MULTIPLIER)
#define __LAUNCH_BOUNDS_IPS_ENERGY__ __LAUNCH_BOUNDS__(IPSNONBONDENERGY_THREADS_PER_BLOCK, IPSNONBONDENERGY_BLOCKS_MULTIPLIER)
#define __LAUNCH_BOUNDS_IPS_FORCES__ __LAUNCH_BOUNDS__(IPSNONBONDFORCES_THREADS_PER_BLOCK, IPSNONBONDENERGY_BLOCKS_MULTIPLIER)

#define PME_FSWITCH
#define PME_ATOMS_PER_WARP (32)
#define PME_VIRIAL
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBFrcVir32fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBFrcVir32fswitchPHMD_kernel()
#include INC_NLCPNE

#undef PHMD
//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBNrgVir32fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBNrgVir32fswitchPHMD_kernel()
#include INC_NLCPNE

#undef PHMD
#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBFrcVir32fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBFrcVir32fswitchPHMD_kernel()
#include INC_NLCPNE
#undef PHMD
#define PME_IS_ORTHOGONAL
#define PME_ENERGY

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBNrgVir32fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBNrgVir32fswitchPHMD_kernel()
#include INC_NLCPNE
#undef PHMD
#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_VIRIAL

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMENBFrc32fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDFORCES_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBFrc32fswitchPHMD_kernel()
#include INC_NLCPNE
#undef PHMD

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBNrg32fswitch_kernel()
#include INC_NLCPNE
//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBNrg32fswitchPHMD_kernel()
#include INC_NLCPNE
#undef PHMD

#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMEOrthoNBFrc32fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDFORCES_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBFrc32fswitchPHMD_kernel()
#include INC_NLCPNE
#undef PHMD
#define PME_ENERGY

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBNrg32fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBNrg32fswitchPHMD_kernel()
#include INC_NLCPNE
#undef PHMD
#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_ATOMS_PER_WARP

#define PME_ATOMS_PER_WARP (16)
#define PME_VIRIAL

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBFrcVir16fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBFrcVir16fswitchPHMD_kernel()
#include INC_NLCPNE
#undef PHMD
#define PME_ENERGY

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBNrgVir16fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBNrgVir16fswitchPHMD_kernel()
#include INC_NLCPNE
#undef PHMD
#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBFrcVir16fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBFrcVir16fswitchPHMD_kernel()
#include INC_NLCPNE
#undef PHMD
#define PME_ENERGY

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBNrgVir16fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBNrgVir16fswitchPHMD_kernel()
#include INC_NLCPNE
#undef PHMD
#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_VIRIAL

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMENBFrc16fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDFORCES_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBFrc16fswitchPHMD_kernel()
#include INC_NLCPNE
#undef PHMD
#define PME_ENERGY

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBNrg16fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBNrg16fswitchPHMD_kernel()
#include INC_NLCPNE
#undef PHMD
#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMEOrthoNBFrc16fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDFORCES_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBFrc16fswitchPHMD_kernel()
#include INC_NLCPNE
#undef PHMD
//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBNrg16fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBNrg16fswitchPHMD_kernel()
#include INC_NLCPNE
#undef PHMD
#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_ATOMS_PER_WARP

//---------------------------------------------------------------------------------------------
#define PME_ATOMS_PER_WARP (8)
#define PME_VIRIAL
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBFrcVir8fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBFrcVir8fswitchPHMD_kernel()
#include INC_NLCPNE
#undef PHMD

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBNrgVir8fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBNrgVir8fswitchPHMD_kernel()
#include INC_NLCPNE
#undef PHMD
#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBFrcVir8fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBFrcVir8fswitchPHMD_kernel()
#include INC_NLCPNE
#undef PHMD

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBNrgVir8fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBNrgVir8fswitchPHMD_kernel()
#include INC_NLCPNE
#undef PHMD

//---------------------------------------------------------------------------------------------

#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_VIRIAL

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMENBFrc8fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDFORCES_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBFrc8fswitchPHMD_kernel()
#include INC_NLCPNE
#undef PHMD

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBNrg8fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBNrg8fswitchPHMD_kernel()
#include INC_NLCPNE
#undef PHMD
#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMEOrthoNBFrc8fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDFORCES_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBFrc8fswitchPHMD_kernel()
#include INC_NLCPNE
#undef PHMD

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBNrg8fswitch_kernel()
#include INC_NLCPNE
//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBNrg8fswitchPHMD_kernel()
#include INC_NLCPNE

#undef PHMD
#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_ATOMS_PER_WARP

//---------------------------------------------------------------------------------------------
#define PME_MINIMIZATION
#define PME_ATOMS_PER_WARP (32)
#define PME_VIRIAL
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBMiniFrcVir32fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBMiniNrgVir32fswitch_kernel()
#include INC_NLCPNE

#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBMiniFrcVir32fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBMiniNrgVir32fswitch_kernel()
#include INC_NLCPNE

#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_VIRIAL

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMENBMiniFrc32fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBMiniNrg32fswitch_kernel()
#include INC_NLCPNE

#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMEOrthoNBMiniFrc32fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBMiniNrg32fswitch_kernel()
#include INC_NLCPNE

#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_ATOMS_PER_WARP

//---------------------------------------------------------------------------------------------
#define PME_ATOMS_PER_WARP (16)
#define PME_VIRIAL

__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBMiniFrcVir16fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBMiniNrgVir16fswitch_kernel()
#include INC_NLCPNE

#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBMiniFrcVir16fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBMiniNrgVir16fswitch_kernel()
#include INC_NLCPNE

#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_VIRIAL

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMENBMiniFrc16fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBMiniNrg16fswitch_kernel()
#include INC_NLCPNE

#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMEOrthoNBMiniFrc16fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBMiniNrg16fswitch_kernel()
#include INC_NLCPNE

#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_ATOMS_PER_WARP

//---------------------------------------------------------------------------------------------
#define PME_ATOMS_PER_WARP (8)
#define PME_VIRIAL
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBMiniFrcVir8fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBMiniNrgVir8fswitch_kernel()
#include INC_NLCPNE

#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBMiniFrcVir8fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBMiniNrgVir8fswitch_kernel()
#include INC_NLCPNE

#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_VIRIAL

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMENBMiniFrc8fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBMiniNrg8fswitch_kernel()
#include INC_NLCPNE

#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMEOrthoNBMiniFrc8fswitch_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBMiniNrg8fswitch_kernel()
#include INC_NLCPNE

#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_ATOMS_PER_WARP
#undef PME_MINIMIZATION
#undef PME_FSWITCH

//---------------------------------------------------------------------------------------------
#define PME_ATOMS_PER_WARP (32)
#define PME_VIRIAL
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBFrcVir32_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBFrcVir32PHMD_kernel()
#include INC_NLCPNE
#undef PHMD

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBNrgVir32_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBNrgVir32PHMD_kernel()
#include INC_NLCPNE
#undef PHMD
#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBFrcVir32_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBFrcVir32PHMD_kernel()
#include INC_NLCPNE
#undef PHMD

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBNrgVir32_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBNrgVir32PHMD_kernel()
#include INC_NLCPNE
#undef PHMD
#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_VIRIAL

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMENBFrc32_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDFORCES_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBFrc32PHMD_kernel()
#include INC_NLCPNE
#undef PHMD

//---------------------------------------------------------------------------------------------
#define PME_SCTI
#define PME_SCMBAR
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMENBFrc32SCMBAR_kernel()
#include "kNLCPNE_AFE.h"
#undef PME_SCMBAR

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMENBFrc32SCTI_kernel()
#include "kNLCPNE_AFE.h"

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBNrg32SCTI_kernel()
#include "kNLCPNE_AFE.h"

//---------------------------------------------------------------------------------------------
#define PME_SCMBAR
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBNrg32SCMBAR_kernel()
#include "kNLCPNE_AFE.h"
#undef PME_SCMBAR
#undef PME_SCTI


//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBNrg32_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBNrg32PHMD_kernel()
#include INC_NLCPNE
#undef PHMD
#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMEOrthoNBFrc32_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDFORCES_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBFrc32PHMD_kernel()
#include INC_NLCPNE
#undef PHMD
#define PME_SCTI

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMEOrthoNBFrc32SCTI_kernel()
#include "kNLCPNE_AFE.h"

//---------------------------------------------------------------------------------------------
#define PME_SCMBAR
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMEOrthoNBFrc32SCMBAR_kernel()
#include "kNLCPNE_AFE.h"

#undef PME_SCMBAR
#undef PME_SCTI

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBNrg32_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBNrg32PHMD_kernel()
#include INC_NLCPNE
#undef PHMD
//---------------------------------------------------------------------------------------------
#define PME_SCTI
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBNrg32SCTI_kernel()
#include "kNLCPNE_AFE.h"

//---------------------------------------------------------------------------------------------
#define PME_SCMBAR
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBNrg32SCMBAR_kernel()
#include "kNLCPNE_AFE.h"

#undef PME_SCMBAR
#undef PME_SCTI
#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_ATOMS_PER_WARP

//---------------------------------------------------------------------------------------------
#define PME_ATOMS_PER_WARP (16)
#define PME_VIRIAL

__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBFrcVir16_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBFrcVir16PHMD_kernel()
#include INC_NLCPNE
#undef PHMD
//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBNrgVir16_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBNrgVir16PHMD_kernel()
#include INC_NLCPNE

#undef PHMD
#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBFrcVir16_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBFrcVir16PHMD_kernel()
#include INC_NLCPNE
#undef PHMD

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBNrgVir16_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBNrgVir16PHMD_kernel()
#include INC_NLCPNE
#undef PHMD
#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_VIRIAL

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMENBFrc16_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDFORCES_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBFrc16PHMD_kernel()
#include INC_NLCPNE
#undef PHMD

//---------------------------------------------------------------------------------------------
#define PME_SCTI
#define PME_SCMBAR
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMENBFrc16SCMBAR_kernel()
#include "kNLCPNE_AFE.h"

#undef PME_SCMBAR

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMENBFrc16SCTI_kernel()
#include "kNLCPNE_AFE.h"

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBNrg16SCTI_kernel()
#include "kNLCPNE_AFE.h"

//---------------------------------------------------------------------------------------------
#define PME_SCMBAR
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBNrg16SCMBAR_kernel()
#include "kNLCPNE_AFE.h"

#undef PME_SCMBAR
#undef PME_SCTI

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBNrg16_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBNrg16PHMD_kernel()
#include INC_NLCPNE
#undef PHMD
#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMEOrthoNBFrc16_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDFORCES_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBFrc16PHMD_kernel()
#include INC_NLCPNE
#undef PHMD

//---------------------------------------------------------------------------------------------
#define PME_SCTI
#define PME_SCMBAR
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMEOrthoNBFrc16SCMBAR_kernel()
#include "kNLCPNE_AFE.h"

#undef PME_SCMBAR

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMEOrthoNBFrc16SCTI_kernel()
#include "kNLCPNE_AFE.h"

#undef PME_SCTI

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBNrg16_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBNrg16PHMD_kernel()
#include INC_NLCPNE
#undef PHMD

//---------------------------------------------------------------------------------------------
#define PME_SCTI
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBNrg16SCTI_kernel()
#include "kNLCPNE_AFE.h"


//---------------------------------------------------------------------------------------------
#define PME_SCMBAR
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBNrg16SCMBAR_kernel()
#include "kNLCPNE_AFE.h"

#undef PME_SCMBAR
#undef PME_SCTI
#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_ATOMS_PER_WARP

//---------------------------------------------------------------------------------------------
#define PME_ATOMS_PER_WARP (8)
#define PME_VIRIAL
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBFrcVir8_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBFrcVir8PHMD_kernel()
#include INC_NLCPNE
#undef PHMD

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBNrgVir8_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBNrgVir8PHMD_kernel()
#include INC_NLCPNE
#undef PHMD
#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBFrcVir8_kernel()
#include INC_NLCPNE
//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBFrcVir8PHMD_kernel()
#include INC_NLCPNE

#undef PHMD
//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBNrgVir8_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBNrgVir8PHMD_kernel()
#include INC_NLCPNE

#undef PHMD
#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_VIRIAL

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMENBFrc8_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDFORCES_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBFrc8PHMD_kernel()
#include INC_NLCPNE

#undef PHMD
//---------------------------------------------------------------------------------------------
#define PME_SCTI
#define PME_SCMBAR
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMENBFrc8SCMBAR_kernel()
#include "kNLCPNE_AFE.h"

#undef PME_SCMBAR

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMENBFrc8SCTI_kernel()
#include "kNLCPNE_AFE.h"

#undef PME_SCTI

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBNrg8_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMENBNrg8PHMD_kernel()
#include INC_NLCPNE

#undef PHMD

//---------------------------------------------------------------------------------------------
#define PME_SCTI
#define PME_SCMBAR
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBNrg8SCMBAR_kernel()
#include "kNLCPNE_AFE.h"

#undef PME_SCMBAR

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBNrg8SCTI_kernel()
#include "kNLCPNE_AFE.h"

#undef PME_SCTI
#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMEOrthoNBFrc8_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDFORCES_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBFrc8PHMD_kernel()
#include INC_NLCPNE

#undef PHMD
//---------------------------------------------------------------------------------------------
#define PME_SCTI
#define PME_SCMBAR
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMEOrthoNBFrc8SCMBAR_kernel()
#include "kNLCPNE_AFE.h"

#undef PME_SCMBAR

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMEOrthoNBFrc8SCTI_kernel()
#include "kNLCPNE_AFE.h"

#undef PME_SCTI

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBNrg8_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PHMD
__global__ void
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
kCalcPMEOrthoNBNrg8PHMD_kernel()
#include INC_NLCPNE

#undef PHMD

//---------------------------------------------------------------------------------------------
#define PME_SCTI
#define PME_SCMBAR
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBNrg8SCMBAR_kernel()
#include "kNLCPNE_AFE.h"

#undef PME_SCMBAR

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBNrg8SCTI_kernel()
#include "kNLCPNE_AFE.h"

#undef PME_SCTI
#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_ATOMS_PER_WARP

//---------------------------------------------------------------------------------------------
#define PME_MINIMIZATION
#define PME_ATOMS_PER_WARP (32)
#define PME_VIRIAL
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBMiniFrcVir32_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBMiniNrgVir32_kernel()
#include INC_NLCPNE

#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBMiniFrcVir32_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBMiniNrgVir32_kernel()
#include INC_NLCPNE

#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_VIRIAL

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__

kCalcPMENBMiniFrc32_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBMiniNrg32_kernel()
#include INC_NLCPNE

#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMEOrthoNBMiniFrc32_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBMiniNrg32_kernel()
#include INC_NLCPNE

#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_ATOMS_PER_WARP

//---------------------------------------------------------------------------------------------
#define PME_ATOMS_PER_WARP (16)
#define PME_VIRIAL
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBMiniFrcVir16_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBMiniNrgVir16_kernel()
#include INC_NLCPNE

#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBMiniFrcVir16_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBMiniNrgVir16_kernel()
#include INC_NLCPNE

#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_VIRIAL

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMENBMiniFrc16_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBMiniNrg16_kernel()
#include INC_NLCPNE

#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMEOrthoNBMiniFrc16_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBMiniNrg16_kernel()
#include INC_NLCPNE

#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_ATOMS_PER_WARP

//---------------------------------------------------------------------------------------------
#define PME_ATOMS_PER_WARP (8)
#define PME_VIRIAL
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBMiniFrcVir8_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBMiniNrgVir8_kernel()
#include INC_NLCPNE

#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBMiniFrcVir8_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBMiniNrgVir8_kernel()
#include INC_NLCPNE

#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_VIRIAL

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMENBMiniFrc8_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMENBMiniNrg8_kernel()
#include INC_NLCPNE

#undef PME_ENERGY

//---------------------------------------------------------------------------------------------
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_PME_FORCES__
kCalcPMEOrthoNBMiniFrc8_kernel()
#include INC_NLCPNE

//---------------------------------------------------------------------------------------------
#define PME_ENERGY
__global__ void
__LAUNCH_BOUNDS_PME_ENERGY__
kCalcPMEOrthoNBMiniNrg8_kernel()
#include INC_NLCPNE

#undef PME_ENERGY
#undef PME_IS_ORTHOGONAL
#undef PME_ATOMS_PER_WARP
#undef PME_MINIMIZATION

//---------------------------------------------------------------------------------------------
// kCalculatePMENonbondForces: host function to launch the non-bonded calculation and get
//                             forces.  The non-bonded direct space calculation is, by a good
//                             margin, the most expensive part of the simulation and the vast
//                             majority of steps will require only force computations.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kCalculatePMENonbondForces(gpuContext gpu)
{
  // Set local variables to the launch bounds that we need, to make the massive
  // case switch below easier to parse.  There are subtle differences between
  // the launch bounds that can get lost in the lengthy names--this will help
  // with that problem.
  int nbBlocks = gpu->PMENonbondBlocks;
  int nrgThreads = gpu->PMENonbondEnergyThreadsPerBlock;
  int frcThreads = gpu->PMENonbondForcesThreadsPerBlock;

  // Decide which non-bonded direct space force kernel to use.
  if (gpu->sim.fswitch < 0) {
    if (gpu->imin == 0) {
      if (gpu->sim.NLAtomsPerWarp == 32) {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->iphmd == 3) {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBFrcVir32PHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBFrcVir32PHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
          else {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBFrcVir32_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBFrcVir32_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
        }
        else {
          if (gpu->sim.is_orthog) {
            if (gpu->sim.ti_mode > 0) {
              if (gpu->sim.ifmbar > 0) {
                kCalcPMEOrthoNBFrc32SCMBAR_kernel<<<nbBlocks, frcThreads>>>();
              }
              else {
                kCalcPMEOrthoNBFrc32SCTI_kernel<<<nbBlocks, frcThreads>>>();
              }
            }
            else {
              if (gpu->iphmd == 3) {
                kCalcPMEOrthoNBFrc32PHMD_kernel<<<nbBlocks, frcThreads>>>();
              }
              else {
                kCalcPMEOrthoNBFrc32_kernel<<<nbBlocks, frcThreads>>>();
              }
            }
          }
          else {
            if (gpu->sim.ti_mode > 0) {
              if (gpu->sim.ifmbar > 0) {
                kCalcPMENBFrc32SCMBAR_kernel<<<nbBlocks, frcThreads>>>();
              }
              else {
                kCalcPMENBFrc32SCTI_kernel<<<nbBlocks, frcThreads>>>();
              }
            }
            else {
              if (gpu->iphmd == 3) {
                kCalcPMENBFrc32PHMD_kernel<<<nbBlocks, frcThreads>>>();
              }
              else {
                kCalcPMENBFrc32_kernel<<<nbBlocks, frcThreads>>>();
              }
            }
          }
        }
      }
      else if (gpu->sim.NLAtomsPerWarp == 16) {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->iphmd == 3) {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBFrcVir16PHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBFrcVir16PHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
          else {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBFrcVir16_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBFrcVir16_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
        }
        else {
          if (gpu->sim.is_orthog) {
            if (gpu->sim.ti_mode > 0) {
              if (gpu->sim.ifmbar > 0) {
                kCalcPMEOrthoNBFrc16SCMBAR_kernel<<<nbBlocks, frcThreads>>>();
              }
              else {
                kCalcPMEOrthoNBFrc16SCTI_kernel<<<nbBlocks, frcThreads>>>();
              }
            }
            else {
              if (gpu->iphmd == 3) {
                kCalcPMEOrthoNBFrc16PHMD_kernel<<<nbBlocks, frcThreads>>>();
              }
              else {
                kCalcPMEOrthoNBFrc16_kernel<<<nbBlocks, frcThreads>>>();
              }
            }
          }
          else {
            if (gpu->sim.ti_mode > 0) {
              if (gpu->sim.ifmbar > 0) {
                kCalcPMENBFrc16SCMBAR_kernel<<<nbBlocks, frcThreads>>>();
              }
              else {
                kCalcPMENBFrc16SCTI_kernel<<<nbBlocks, frcThreads>>>();
              }
            }
            else {
              if (gpu->iphmd == 3) {
                kCalcPMENBFrc16PHMD_kernel<<<nbBlocks, frcThreads>>>();
              }
              else {
                kCalcPMENBFrc16_kernel<<<nbBlocks, frcThreads>>>();
              }
            }
          }
        }
      }
      else {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->iphmd == 3) {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBFrcVir8PHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBFrcVir8PHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
          else {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBFrcVir8_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBFrcVir8_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
        }
        else {
          if (gpu->sim.is_orthog) {
            if (gpu->sim.ti_mode > 0) {
              if (gpu->sim.ifmbar > 0) {
                kCalcPMEOrthoNBFrc8SCMBAR_kernel<<<nbBlocks, frcThreads>>>();
              }
              else {
                kCalcPMEOrthoNBFrc8SCTI_kernel<<<nbBlocks, frcThreads>>>();
              }
            }
            else {
              if (gpu->iphmd == 3) {
                kCalcPMEOrthoNBFrc8PHMD_kernel<<<nbBlocks, frcThreads>>>();
              }
              else {
                kCalcPMEOrthoNBFrc8_kernel<<<nbBlocks, frcThreads>>>();
              }
            }
          }
          else {
            if (gpu->sim.ti_mode > 0) {
              if (gpu->sim.ifmbar > 0) {
                kCalcPMENBFrc8SCMBAR_kernel<<<nbBlocks, frcThreads>>>();
              }
              else {
                kCalcPMENBFrc8SCTI_kernel<<<nbBlocks, frcThreads>>>();
              }
            }
            else {
              if (gpu->iphmd == 3) {
                kCalcPMENBFrc8PHMD_kernel<<<nbBlocks, frcThreads>>>();
              }
              else {
                kCalcPMENBFrc8_kernel<<<nbBlocks, frcThreads>>>();
              }
            }
          }
        }
      }
    }
    else {
      if (gpu->sim.NLAtomsPerWarp == 32) {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniFrcVir32_kernel<<<nbBlocks, nrgThreads>>>();
          }
          else {
            kCalcPMENBMiniFrcVir32_kernel<<<nbBlocks, nrgThreads>>>();
          }
        }
        else {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniFrc32_kernel<<<nbBlocks, frcThreads>>>();
          }
          else {
            kCalcPMENBMiniFrc32_kernel<<<nbBlocks, frcThreads>>>();
          }
        }
      }
      else if (gpu->sim.NLAtomsPerWarp == 16) {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniFrcVir16_kernel<<<nbBlocks, nrgThreads>>>();
          }
          else {
            kCalcPMENBMiniFrcVir16_kernel<<<nbBlocks, nrgThreads>>>();
          }
        }
        else {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniFrc16_kernel<<<nbBlocks, frcThreads>>>();
          }
          else {
            kCalcPMENBMiniFrc16_kernel<<<nbBlocks, frcThreads>>>();
          }
        }
      }
      else {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniFrcVir8_kernel<<<nbBlocks, nrgThreads>>>();
          }
          else {
            kCalcPMENBMiniFrcVir8_kernel<<<nbBlocks, nrgThreads>>>();
          }
        }
        else {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniFrc8_kernel<<<nbBlocks, frcThreads>>>();
          }
          else {
            kCalcPMENBMiniFrc8_kernel<<<nbBlocks, frcThreads>>>();
          }
        }
      }
    }
  }
  else {
    if (gpu->imin == 0) {
      if (gpu->sim.NLAtomsPerWarp == 32) {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->iphmd == 3) {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBFrcVir32fswitchPHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBFrcVir32fswitchPHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
          else {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBFrcVir32fswitch_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBFrcVir32fswitch_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
        }
        else {
          if (gpu->iphmd == 3) {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBFrc32fswitchPHMD_kernel<<<nbBlocks, frcThreads>>>();
            }
            else {
              kCalcPMENBFrc32fswitchPHMD_kernel<<<nbBlocks, frcThreads>>>();
            }
          }
          else {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBFrc32fswitch_kernel<<<nbBlocks, frcThreads>>>();
            }
            else {
              kCalcPMENBFrc32fswitch_kernel<<<nbBlocks, frcThreads>>>();
            }
          }
        }
      }
      else if (gpu->sim.NLAtomsPerWarp == 16) {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->iphmd == 3) {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBFrcVir16fswitchPHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBFrcVir16fswitchPHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
          else {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBFrcVir16fswitch_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBFrcVir16fswitch_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
        }
        else {
          if (gpu->iphmd == 3) {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBFrc16fswitchPHMD_kernel<<<nbBlocks, frcThreads>>>();
            }
            else {
              kCalcPMENBFrc16fswitchPHMD_kernel<<<nbBlocks, frcThreads>>>();
            }
          }
          else {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBFrc16fswitch_kernel<<<nbBlocks, frcThreads>>>();
            }
            else {
              kCalcPMENBFrc16fswitch_kernel<<<nbBlocks, frcThreads>>>();
            }
          }
        }
      }
      else {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->iphmd == 3) {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBFrcVir8fswitchPHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBFrcVir8fswitchPHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
          else {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBFrcVir8fswitch_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBFrcVir8fswitch_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
        }
        else {
          if (gpu->iphmd == 3) {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBFrc8fswitchPHMD_kernel<<<nbBlocks, frcThreads>>>();
            }
            else {
              kCalcPMENBFrc8fswitchPHMD_kernel<<<nbBlocks, frcThreads>>>();
            }
          }
          else {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBFrc8fswitch_kernel<<<nbBlocks, frcThreads>>>();
            }
            else {
              kCalcPMENBFrc8fswitch_kernel<<<nbBlocks, frcThreads>>>();
            }
          }
        }
      }
    }
    else {
      if (gpu->sim.NLAtomsPerWarp == 32) {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniFrcVir32fswitch_kernel<<<nbBlocks, nrgThreads>>>();
          }
          else {
            kCalcPMENBMiniFrcVir32fswitch_kernel<<<nbBlocks, nrgThreads>>>();
          }
        }
        else {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniFrc32fswitch_kernel<<<nbBlocks, frcThreads>>>();
          }
          else {
            kCalcPMENBMiniFrc32fswitch_kernel<<<nbBlocks, frcThreads>>>();
          }
        }
      }
      else if (gpu->sim.NLAtomsPerWarp == 16) {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniFrcVir16fswitch_kernel<<<nbBlocks, nrgThreads>>>();
          }
          else {
            kCalcPMENBMiniFrcVir16fswitch_kernel<<<nbBlocks, nrgThreads>>>();
          }
        }
        else {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniFrc16fswitch_kernel<<<nbBlocks, frcThreads>>>();
          }
          else {
            kCalcPMENBMiniFrc16fswitch_kernel<<<nbBlocks, frcThreads>>>();
          }
        }
      }
      else {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniFrcVir8fswitch_kernel<<<nbBlocks, nrgThreads>>>();
          }
          else {
            kCalcPMENBMiniFrcVir8fswitch_kernel<<<nbBlocks, nrgThreads>>>();
          }
        }
        else {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniFrc8fswitch_kernel<<<nbBlocks, frcThreads>>>();
          }
          else {
            kCalcPMENBMiniFrc8fswitch_kernel<<<nbBlocks, frcThreads>>>();
          }
        }
      }
    }
  }
  LAUNCHERROR("kCalculatePMENonbondForces");
}

//---------------------------------------------------------------------------------------------
// kCalculatePMENonbondEnergy: similar to the function above for forces, but now invoking an
//                             energy calculation as the non-bonded interactions are examined.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kCalculatePMENonbondEnergy(gpuContext gpu)
{
  // Set local variables to the launch bounds that we need.  Note that
  // PMENonbondForcesThreadsPerBlock is used exclusively to launch kernels
  // for force calculations (see above), but if virials are needed as part
  // of the force calculation the energy thread count was used.
  int nbBlocks = gpu->PMENonbondBlocks;
  int nrgThreads = gpu->PMENonbondEnergyThreadsPerBlock;

  // Decide which non-bonded direct space energy kernel to use.
  if (gpu->sim.fswitch < 0) {
    if (gpu->imin == 0) {
      if (gpu->sim.NLAtomsPerWarp == 32) {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->iphmd == 3) {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBNrgVir32PHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBNrgVir32PHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
          else {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBNrgVir32_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBNrgVir32_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
        }
        else {
          if (gpu->sim.is_orthog) {
            if (gpu->sim.ti_mode > 0) {
              if (gpu->sim.ifmbar > 0) {
                kCalcPMEOrthoNBNrg32SCMBAR_kernel<<<nbBlocks, nrgThreads>>>();
              }
              else {
                kCalcPMEOrthoNBNrg32SCTI_kernel<<<nbBlocks, nrgThreads>>>();
              }
            }
            else {
              if (gpu->iphmd == 3) {
                kCalcPMEOrthoNBNrg32PHMD_kernel<<<nbBlocks, nrgThreads>>>();
              }
              else {
                kCalcPMEOrthoNBNrg32_kernel<<<nbBlocks, nrgThreads>>>();
              }
            }
          }
          else {
            if (gpu->sim.ti_mode > 0) {
              if (gpu->sim.ifmbar > 0) {
                kCalcPMENBNrg32SCMBAR_kernel<<<nbBlocks, nrgThreads>>>();
              }
              else {
                kCalcPMENBNrg32SCTI_kernel<<<nbBlocks, nrgThreads>>>();
              }
            }
            else {
              if (gpu->iphmd == 3) {
                kCalcPMENBNrg32PHMD_kernel<<<nbBlocks, nrgThreads>>>();
              }
              else {
                kCalcPMENBNrg32_kernel<<<nbBlocks, nrgThreads>>>();
              }
            }
          }
        }
      }
      else if (gpu->sim.NLAtomsPerWarp == 16) {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->iphmd == 3) {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBNrgVir16PHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBNrgVir16PHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
          else {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBNrgVir16_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBNrgVir16_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
        }
        else {
          if (gpu->sim.is_orthog) {
            if (gpu->sim.ti_mode > 0) {
              if (gpu->sim.ifmbar > 0) {
                kCalcPMEOrthoNBNrg16SCMBAR_kernel<<<nbBlocks, nrgThreads>>>();
              }
              else {
                kCalcPMEOrthoNBNrg16SCTI_kernel<<<nbBlocks, nrgThreads>>>();
              }
            }
            else {
              if (gpu->iphmd == 3) {
                kCalcPMEOrthoNBNrg16PHMD_kernel<<<nbBlocks, nrgThreads>>>();
              }
              else {
                kCalcPMEOrthoNBNrg16_kernel<<<nbBlocks, nrgThreads>>>();
              }
            }
          }
          else {
            if (gpu->sim.ti_mode > 0) {
              if (gpu->sim.ifmbar > 0) {
                kCalcPMENBNrg16SCMBAR_kernel<<<nbBlocks, nrgThreads>>>();
              }
              else {
                kCalcPMENBNrg16SCTI_kernel<<<nbBlocks, nrgThreads>>>();
              }
            }
            else {
              if (gpu->iphmd == 3) {
                kCalcPMENBNrg16PHMD_kernel<<<nbBlocks, nrgThreads>>>();
              }
              else {
                kCalcPMENBNrg16_kernel<<<nbBlocks, nrgThreads>>>();
              }
            }
          }
        }
      }
      else {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->iphmd == 3) {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBNrgVir8PHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBNrgVir8PHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
          else {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBNrgVir8_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBNrgVir8_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
        }
        else {
          if (gpu->sim.is_orthog) {
            if (gpu->sim.ti_mode > 0) {
              if (gpu->sim.ifmbar > 0) {
                kCalcPMEOrthoNBNrg8SCMBAR_kernel<<<nbBlocks, nrgThreads>>>();
              }
              else {
                kCalcPMEOrthoNBNrg8SCTI_kernel<<<nbBlocks, nrgThreads>>>();
              }
            }
            else {
              if (gpu->iphmd == 3) {
                kCalcPMEOrthoNBNrg8PHMD_kernel<<<nbBlocks, nrgThreads>>>();
              }
              else {
                kCalcPMEOrthoNBNrg8_kernel<<<nbBlocks, nrgThreads>>>();
              }
            }
          }
          else {
            if (gpu->sim.ti_mode > 0) {
              if (gpu->sim.ifmbar > 0) {
                kCalcPMENBNrg8SCMBAR_kernel<<<nbBlocks, nrgThreads>>>();
              }
              else {
                kCalcPMENBNrg8SCTI_kernel<<<nbBlocks, nrgThreads>>>();
              }
            }
            else {
              if (gpu->iphmd == 3) {
                kCalcPMENBNrg8PHMD_kernel<<<nbBlocks, nrgThreads>>>();
              }
              else {
                kCalcPMENBNrg8_kernel<<<nbBlocks, nrgThreads>>>();
              }
            }
          }
        }
      }
    }
    else {
      if (gpu->sim.NLAtomsPerWarp == 32) {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniNrgVir32_kernel<<<nbBlocks, nrgThreads>>>();
          }
          else {
            kCalcPMENBMiniNrgVir32_kernel<<<nbBlocks, nrgThreads>>>();
          }
        }
        else {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniNrg32_kernel<<<nbBlocks, nrgThreads>>>();
          }
          else {
            kCalcPMENBMiniNrg32_kernel<<<nbBlocks, nrgThreads>>>();
          }
        }
      }
      else if (gpu->sim.NLAtomsPerWarp == 16) {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniNrgVir16_kernel<<<nbBlocks, nrgThreads>>>();
          }
          else {
            kCalcPMENBMiniNrgVir16_kernel<<<nbBlocks, nrgThreads>>>();
          }
        }
        else {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniNrg16_kernel<<<nbBlocks, nrgThreads>>>();
          }
          else {
            kCalcPMENBMiniNrg16_kernel<<<nbBlocks, nrgThreads>>>();
          }
        }
      }
      else {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniNrgVir8_kernel<<<nbBlocks, nrgThreads>>>();
          }
          else {
            kCalcPMENBMiniNrgVir8_kernel<<<nbBlocks, nrgThreads>>>();
          }
        }
        else {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniNrg8_kernel<<<nbBlocks, nrgThreads>>>();
          }
          else {
            kCalcPMENBMiniNrg8_kernel<<<nbBlocks, nrgThreads>>>();
          }
        }
      }
    }
  }
  else {
    if (gpu->imin == 0) {
      if (gpu->sim.NLAtomsPerWarp == 32) {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->iphmd == 3) {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBNrgVir32fswitchPHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBNrgVir32fswitchPHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
          else {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBNrgVir32fswitch_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBNrgVir32fswitch_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
        }
        else {
          if (gpu->iphmd == 3) {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBNrg32fswitchPHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBNrg32fswitchPHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
          else {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBNrg32fswitch_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBNrg32fswitch_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
        }
      }
      else if (gpu->sim.NLAtomsPerWarp == 16) {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->iphmd == 3) {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBNrgVir16fswitchPHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBNrgVir16fswitchPHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
          else {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBNrgVir16fswitch_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBNrgVir16fswitch_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
        }
        else {
          if (gpu->iphmd == 3) {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBNrg16fswitchPHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBNrg16fswitchPHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
          else {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBNrg16fswitch_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBNrg16fswitch_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
        }
      }
      else {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->iphmd == 3) {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBNrgVir8fswitchPHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBNrgVir8fswitchPHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
          else {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBNrgVir8fswitch_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBNrgVir8fswitch_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
        }
        else {
          if (gpu->iphmd == 3) {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBNrg8fswitchPHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBNrg8fswitchPHMD_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
          else {
            if (gpu->sim.is_orthog) {
              kCalcPMEOrthoNBNrg8fswitch_kernel<<<nbBlocks, nrgThreads>>>();
            }
            else {
              kCalcPMENBNrg8fswitch_kernel<<<nbBlocks, nrgThreads>>>();
            }
          }
        }
      }
    }
    else {
      if (gpu->sim.NLAtomsPerWarp == 32) {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniNrgVir32_kernel<<<nbBlocks, nrgThreads>>>();
          }
          else {
            kCalcPMENBMiniNrgVir32_kernel<<<nbBlocks, nrgThreads>>>();
          }
        }
        else {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniNrg32_kernel<<<nbBlocks, nrgThreads>>>();
          }
          else {
            kCalcPMENBMiniNrg32_kernel<<<nbBlocks, nrgThreads>>>();
          }
        }
      }
      else if (gpu->sim.NLAtomsPerWarp == 16) {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniNrgVir16_kernel<<<nbBlocks, nrgThreads>>>();
          }
          else {
            kCalcPMENBMiniNrgVir16_kernel<<<nbBlocks, nrgThreads>>>();
          }
        }
        else {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniNrg16_kernel<<<nbBlocks, nrgThreads>>>();
          }
          else {
            kCalcPMENBMiniNrg16_kernel<<<nbBlocks, nrgThreads>>>();
          }
        }
      }
      else {
        if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniNrgVir8_kernel<<<nbBlocks, nrgThreads>>>();
          }
          else {
            kCalcPMENBMiniNrgVir8_kernel<<<nbBlocks, nrgThreads>>>();
          }
        }
        else {
          if (gpu->sim.is_orthog) {
            kCalcPMEOrthoNBMiniNrg8_kernel<<<nbBlocks, nrgThreads>>>();
          }
          else {
            kCalcPMENBMiniNrg8_kernel<<<nbBlocks, nrgThreads>>>();
          }
        }
      }
    }
  }
  LAUNCHERROR("kCalculatePMENonbondEnergy");
}

//---------------------------------------------------------------------------------------------
// Isotropic Periodic Sums nonbonded kernels.  There are not as many of these as PME kernels,
// but still a variety and lots of case specialization.  We still have the following:
//
// - Versions of every kernel taking 8, 16, and 32 atoms per warp
// - Special casing in every kernel for orthorhombic / orthogonal unit cells
// - Special casing in every kernel for energy computation
// - Special casing of virial computations in the context of forces or energies
//
// Aside from that, there is only the further option of specializing the kernel for energy
// minimizations--no alchemical free energy computations seem to be yet available with IPS.
//
// This gives a total of 16 kernels at each warp atom count, 48 IPS non-bonded kernels in all.
//---------------------------------------------------------------------------------------------
#define IPS_ATOMS_PER_WARP (32)
#define IPS_VIRIAL
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSNBFrcVir32_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSNBNrgVir32_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY

//---------------------------------------------------------------------------------------------
#define IPS_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSOrthoNBFrcVir32_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSOrthoNBNrgVir32_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY
#undef IPS_IS_ORTHOGONAL
#undef IPS_VIRIAL

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_IPS_FORCES__
kCalcIPSNBFrc32_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSNBNrg32_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY

//---------------------------------------------------------------------------------------------
#define IPS_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_IPS_FORCES__
kCalcIPSOrthoNBFrc32_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSOrthoNBNrg32_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY
#undef IPS_IS_ORTHOGONAL
#undef IPS_ATOMS_PER_WARP

//---------------------------------------------------------------------------------------------
#define IPS_ATOMS_PER_WARP (16)
#define IPS_VIRIAL
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSNBFrcVir16_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSNBNrgVir16_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY

//---------------------------------------------------------------------------------------------
#define IPS_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSOrthoNBFrcVir16_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSOrthoNBNrgVir16_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY
#undef IPS_IS_ORTHOGONAL
#undef IPS_VIRIAL

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_IPS_FORCES__
kCalcIPSNBFrc16_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSNBNrg16_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY

//---------------------------------------------------------------------------------------------
#define IPS_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_IPS_FORCES__
kCalcIPSOrthoNBFrc16_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSOrthoNBNrg16_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY
#undef IPS_IS_ORTHOGONAL
#undef IPS_ATOMS_PER_WARP

//---------------------------------------------------------------------------------------------
#define IPS_ATOMS_PER_WARP (8)
#define IPS_VIRIAL
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSNBFrcVir8_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSNBNrgVir8_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY

//---------------------------------------------------------------------------------------------
#define IPS_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSOrthoNBFrcVir8_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSOrthoNBNrgVir8_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY
#undef IPS_IS_ORTHOGONAL
#undef IPS_VIRIAL

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_IPS_FORCES__
kCalcIPSNBFrc8_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSNBNrg8_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY

//---------------------------------------------------------------------------------------------
#define IPS_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_IPS_FORCES__
kCalcIPSOrthoNBFrc8_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSOrthoNBNrg8_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY
#undef IPS_IS_ORTHOGONAL
#undef IPS_ATOMS_PER_WARP

//---------------------------------------------------------------------------------------------
#define IPS_MINIMIZATION
#define IPS_ATOMS_PER_WARP (32)
#define IPS_VIRIAL
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSNBMiniFrcVir32_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSNBMiniNrgVir32_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY

//---------------------------------------------------------------------------------------------
#define IPS_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSOrthoNBMiniFrcVir32_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSOrthoNBMiniNrgVir32_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY
#undef IPS_IS_ORTHOGONAL
#undef IPS_VIRIAL

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_IPS_FORCES__
kCalcIPSNBMiniFrc32_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSNBMiniNrg32_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY

//---------------------------------------------------------------------------------------------
#define IPS_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_IPS_FORCES__
kCalcIPSOrthoNBMiniFrc32_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSOrthoNBMiniNrg32_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY
#undef IPS_IS_ORTHOGONAL
#undef IPS_ATOMS_PER_WARP

//---------------------------------------------------------------------------------------------
#define IPS_ATOMS_PER_WARP (16)
#define IPS_VIRIAL
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSNBMiniFrcVir16_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSNBMiniNrgVir16_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY

//---------------------------------------------------------------------------------------------
#define IPS_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSOrthoNBMiniFrcVir16_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSOrthoNBMiniNrgVir16_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY
#undef IPS_IS_ORTHOGONAL
#undef IPS_VIRIAL

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_IPS_FORCES__
kCalcIPSNBMiniFrc16_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSNBMiniNrg16_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY

//---------------------------------------------------------------------------------------------
#define IPS_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_IPS_FORCES__
kCalcIPSOrthoNBMiniFrc16_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSOrthoNBMiniNrg16_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY
#undef IPS_IS_ORTHOGONAL
#undef IPS_ATOMS_PER_WARP

//---------------------------------------------------------------------------------------------
#define IPS_ATOMS_PER_WARP (8)
#define IPS_VIRIAL
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSNBMiniFrcVir8_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSNBMiniNrgVir8_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY

//---------------------------------------------------------------------------------------------
#define IPS_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSOrthoNBMiniFrcVir8_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSOrthoNBMiniNrgVir8_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY
#undef IPS_IS_ORTHOGONAL
#undef IPS_VIRIAL

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_IPS_FORCES__
kCalcIPSNBMiniFrc8_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY
__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSNBMiniNrg8_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY

//---------------------------------------------------------------------------------------------
#define IPS_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_IPS_FORCES__
kCalcIPSOrthoNBMiniFrc8_kernel()
#include "kNLCINE.h"

//---------------------------------------------------------------------------------------------
#define IPS_ENERGY

__global__ void
__LAUNCH_BOUNDS_IPS_ENERGY__
kCalcIPSOrthoNBMiniNrg8_kernel()
#include "kNLCINE.h"

#undef IPS_ENERGY
#undef IPS_IS_ORTHOGONAL
#undef IPS_ATOMS_PER_WARP
#undef IPS_MINIMIZATION

//---------------------------------------------------------------------------------------------
// kCalculatePMENonbondEnergyInitKernels: this is called by gpu_init_ and it does one thing--
//                                        set the __shared__ memory bank size to eight bytes,
//                                        not the default of four, to better align the stride
//                                        of our fast memory with what the non-bonded kernels
//                                        need most.
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kCalculatePMENonbondEnergyInitKernels(gpuContext gpu)
{
#ifdef use_SPFP
#  define PME_SHARED_BANK_SIZE cudaSharedMemBankSizeFourByte
#else
#  define PME_SHARED_BANK_SIZE cudaSharedMemBankSizeEightByte
#endif
#define IPS_SHARED_BANK_SIZE cudaSharedMemBankSizeEightByte

  // Kernels taking 32 atoms per warp
  cudaFuncSetSharedMemConfig(kCalcPMENBFrcVir32_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBNrgVir32_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBFrcVir32_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBNrgVir32_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBFrc32_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBFrc32SCTI_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBFrc32SCMBAR_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBNrg32_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBNrg32SCTI_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBNrg32SCMBAR_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBFrc32_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBFrc32SCTI_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBFrc32SCMBAR_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBNrg32_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBNrg32SCTI_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBNrg32SCMBAR_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBFrcVir32_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBNrgVir32_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBFrcVir32_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBNrgVir32_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBFrc32_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBNrg32_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBFrc32_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBNrg32_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBMiniFrcVir32_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBMiniNrgVir32_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBMiniFrcVir32_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBMiniNrgVir32_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBMiniFrc32_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBMiniNrg32_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBMiniFrc32_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBMiniNrg32_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBMiniFrcVir32_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBMiniNrgVir32_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBMiniFrcVir32_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBMiniNrgVir32_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBMiniFrc32_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBMiniNrg32_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBMiniFrc32_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBMiniNrg32_kernel, IPS_SHARED_BANK_SIZE);

  // Kernels taking 16 atoms per warp
  cudaFuncSetSharedMemConfig(kCalcPMENBFrcVir16_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBNrgVir16_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBFrcVir16_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBNrgVir16_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBFrc16_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBFrc16SCTI_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBFrc16SCMBAR_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBNrg16_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBNrg16SCTI_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBNrg16SCMBAR_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBFrc16_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBFrc16SCTI_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBFrc16SCMBAR_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBNrg16_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBNrg16SCTI_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBNrg16SCMBAR_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBFrcVir16_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBNrgVir16_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBFrcVir16_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBNrgVir16_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBFrc16_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBNrg16_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBFrc16_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBNrg16_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBMiniFrcVir16_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBMiniNrgVir16_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBMiniFrcVir16_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBMiniNrgVir16_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBMiniFrc16_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBMiniNrg16_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBMiniFrc16_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBMiniNrg16_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBMiniFrcVir16_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBMiniNrgVir16_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBMiniFrcVir16_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBMiniNrgVir16_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBMiniFrc16_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBMiniNrg16_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBMiniFrc16_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBMiniNrg16_kernel, IPS_SHARED_BANK_SIZE);

  // Kernels taking 8 atoms per warp
  cudaFuncSetSharedMemConfig(kCalcPMENBFrcVir8_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBNrgVir8_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBFrcVir8_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBNrgVir8_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBFrc8_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBFrc8SCTI_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBFrc8SCMBAR_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBNrg8_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBNrg8SCTI_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBNrg8SCMBAR_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBFrc8_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBFrc8SCTI_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBFrc8SCMBAR_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBNrg8_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBNrg8SCTI_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBNrg8SCMBAR_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBFrcVir8_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBNrgVir8_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBFrcVir8_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBNrgVir8_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBFrc8_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBNrg8_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBFrc8_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBNrg8_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBMiniFrcVir8_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBMiniNrgVir8_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBMiniFrcVir8_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBMiniNrgVir8_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBMiniFrc8_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMENBMiniNrg8_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBMiniFrc8_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcPMEOrthoNBMiniNrg8_kernel, PME_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBMiniFrcVir8_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBMiniNrgVir8_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBMiniFrcVir8_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBMiniNrgVir8_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBMiniFrc8_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSNBMiniNrg8_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBMiniFrc8_kernel, IPS_SHARED_BANK_SIZE);
  cudaFuncSetSharedMemConfig(kCalcIPSOrthoNBMiniNrg8_kernel, IPS_SHARED_BANK_SIZE);
#undef PME_SHARED_BANK_SIZE
#undef IPS_SHARED_BANK_SIZE
}

//---------------------------------------------------------------------------------------------
// kCalculateIPSNonbondForces: compute forces for Isotropic Periodic Sums using the appropriate
//                             kernel.  This follows the format of the corresponding PME
//                             function above.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kCalculateIPSNonbondForces(gpuContext gpu)
{
  // Copy the grid and block dimensions to local variables for brevity
  int nbBlocks = gpu->IPSNonbondBlocks;
  int nrgThreads = gpu->IPSNonbondEnergyThreadsPerBlock;
  int frcThreads = gpu->IPSNonbondForcesThreadsPerBlock;

  // Massive case switch to decide the appropriate non-bonded kernel
  if (gpu->imin == 0) {
    if (gpu->sim.NLAtomsPerWarp == 32) {
      if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
        if (gpu->sim.is_orthog) {
          kCalcIPSOrthoNBFrcVir32_kernel<<<nbBlocks, nrgThreads>>>();
        }
        else {
          kCalcIPSNBFrcVir32_kernel<<<nbBlocks, nrgThreads>>>();
        }
      }
      else {
        if (gpu->sim.is_orthog) {
          kCalcIPSOrthoNBFrc32_kernel<<<nbBlocks, frcThreads>>>();
        }
        else {
          kCalcIPSNBFrc32_kernel<<<nbBlocks, frcThreads>>>();
        }
      }
    }
    else if (gpu->sim.NLAtomsPerWarp == 16) {
      if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
        if (gpu->sim.is_orthog) {
          kCalcIPSOrthoNBFrcVir16_kernel<<<nbBlocks, nrgThreads>>>();
        }
        else {
          kCalcIPSNBFrcVir16_kernel<<<nbBlocks, nrgThreads>>>();
        }
      }
      else {
        if (gpu->sim.is_orthog) {
          kCalcIPSOrthoNBFrc16_kernel<<<nbBlocks, frcThreads>>>();
        }
        else {
          kCalcIPSNBFrc16_kernel<<<nbBlocks, frcThreads>>>();
        }
      }
    }
    else {
      if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
        if (gpu->sim.is_orthog) {
          kCalcIPSOrthoNBFrcVir8_kernel<<<nbBlocks, nrgThreads>>>();
        }
        else {
          kCalcIPSNBFrcVir8_kernel<<<nbBlocks, nrgThreads>>>();
        }
      }
      else {
        if (gpu->sim.is_orthog) {
          kCalcIPSOrthoNBFrc8_kernel<<<nbBlocks, frcThreads>>>();
        }
        else {
          kCalcIPSNBFrc8_kernel<<<nbBlocks, frcThreads>>>();
        }
      }
    }
  }
  else {
    if (gpu->sim.NLAtomsPerWarp == 32) {
      if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
        if (gpu->sim.is_orthog) {
          kCalcIPSOrthoNBMiniFrcVir32_kernel<<<nbBlocks, nrgThreads>>>();
        }
        else {
          kCalcIPSNBMiniFrcVir32_kernel<<<nbBlocks, nrgThreads>>>();
        }
      }
      else {
        if (gpu->sim.is_orthog) {
          kCalcIPSOrthoNBMiniFrc32_kernel<<<nbBlocks, frcThreads>>>();
        }
        else {
          kCalcIPSNBMiniFrc32_kernel<<<nbBlocks, frcThreads>>>();
        }
      }
    }
    else if (gpu->sim.NLAtomsPerWarp == 16) {
      if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
        if (gpu->sim.is_orthog) {
          kCalcIPSOrthoNBMiniFrcVir16_kernel<<<nbBlocks, nrgThreads>>>();
        }
        else {
          kCalcIPSNBMiniFrcVir16_kernel<<<nbBlocks, nrgThreads>>>();
        }
      }
      else {
        if (gpu->sim.is_orthog) {
          kCalcIPSOrthoNBMiniFrc16_kernel<<<nbBlocks, frcThreads>>>();
        }
        else {
          kCalcIPSNBMiniFrc16_kernel<<<nbBlocks, frcThreads>>>();
        }
      }
    }
    else {
      if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
        if (gpu->sim.is_orthog) {
          kCalcIPSOrthoNBMiniFrcVir8_kernel<<<nbBlocks, nrgThreads>>>();
        }
        else {
          kCalcIPSNBMiniFrcVir8_kernel<<<nbBlocks, nrgThreads>>>();
        }
      }
      else {
        if (gpu->sim.is_orthog) {
          kCalcIPSOrthoNBMiniFrc8_kernel<<<nbBlocks, frcThreads>>>();
        }
        else {
          kCalcIPSNBMiniFrc8_kernel<<<nbBlocks, frcThreads>>>();
        }
      }
    }
  }
  LAUNCHERROR("kCalculateIPSNonbondForces");
}

//---------------------------------------------------------------------------------------------
// kCalculateIPSNonbondEnergy: compute forces and energy for Isotropic Periodic Sums using the
//                             appropriate kernel.  This follows the format of the
//                             corresponding PME function above.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kCalculateIPSNonbondEnergy(gpuContext gpu)
{
  // Copy the grid and block dimensions to local variables for brevity
  int nbBlocks = gpu->IPSNonbondBlocks;
  int nrgThreads = gpu->IPSNonbondEnergyThreadsPerBlock;

  // Massive branch to decide the appropriate non-bonded kernel
  if (gpu->imin == 0) {
    if (gpu->sim.NLAtomsPerWarp == 32) {
      if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
        if (gpu->sim.is_orthog) {
          kCalcIPSOrthoNBNrgVir32_kernel<<<nbBlocks, nrgThreads>>>();
        }
        else {
          kCalcIPSNBNrgVir32_kernel<<<nbBlocks, nrgThreads>>>();
        }
      }
      else {
        if (gpu->sim.is_orthog) {
          kCalcIPSOrthoNBNrg32_kernel<<<nbBlocks, nrgThreads>>>();
        }
        else {
          kCalcIPSNBNrg32_kernel<<<nbBlocks, nrgThreads>>>();
        }
      }
    }
    else {
      if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
        if (gpu->sim.is_orthog) {
          kCalcIPSOrthoNBNrgVir16_kernel<<<nbBlocks, nrgThreads>>>();
        }
        else
          kCalcIPSNBNrgVir16_kernel<<<nbBlocks, nrgThreads>>>();
      }
      else {
        if (gpu->sim.is_orthog) {
          kCalcIPSOrthoNBNrg16_kernel<<<nbBlocks, nrgThreads>>>();
        }
        else {
          kCalcIPSNBNrg16_kernel<<<nbBlocks, nrgThreads>>>();
        }
      }
    }
  }
  else {
    if (gpu->sim.NLAtomsPerWarp == 32) {
      if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
        if (gpu->sim.is_orthog) {
          kCalcIPSOrthoNBMiniNrgVir32_kernel<<<nbBlocks, nrgThreads>>>();
        }
        else {
          kCalcIPSNBMiniNrgVir32_kernel<<<nbBlocks, nrgThreads>>>();
        }
      }
      else {
        if (gpu->sim.is_orthog) {
          kCalcIPSOrthoNBMiniNrg32_kernel<<<nbBlocks, nrgThreads>>>();
        }
        else {
          kCalcIPSNBMiniNrg32_kernel<<<nbBlocks, nrgThreads>>>();
        }
      }
    }
    else {
      if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
        if (gpu->sim.is_orthog) {
          kCalcIPSOrthoNBMiniNrgVir16_kernel<<<nbBlocks, nrgThreads>>>();
        }
        else {
          kCalcIPSNBMiniNrgVir16_kernel<<<nbBlocks, nrgThreads>>>();
        }
      }
      else {
        if (gpu->sim.is_orthog) {
          kCalcIPSOrthoNBMiniNrg16_kernel<<<nbBlocks, nrgThreads>>>();
        }
        else {
          kCalcIPSNBMiniNrg16_kernel<<<nbBlocks, nrgThreads>>>();
        }
      }
    }
  }
  LAUNCHERROR("kCalculateIPSNonbondEnergy");
}

//---------------------------------------------------------------------------------------------
// PHMD electrostatics
//---------------------------------------------------------------------------------------------
__global__
__launch_bounds__(PMENONBONDENERGY_THREADS_PER_BLOCK,
                  PMENONBONDENERGY_BLOCKS_MULTIPLIER)
void kCalculateSelfPHMD_kernel(PMEFloat vol, PMEFloat ewaldcof, PMEFloat sumq, PMEFloat sumq2)
{
  int titrindex = (blockIdx.x * blockDim.x + threadIdx.x) * cSim.phmd_stride;
  if (titrindex < cSim.ntitratoms) {
    for( int i = 0; i < cSim.phmd_stride; i++ ) {
      int temp_index = titrindex + i;
      if (temp_index < cSim.ntitratoms) {
        int atmi = cSim.ptitratoms[titrindex + i];
        int atmx = cSim.pImageAtomLookup[atmi];
        int h = cSim.pgrplist[atmi] - 1;
        PMEFloat dudl = 0;
        PMEFloat dudlplus = 0;
        PMEFloat qh = 0;
        PMEFloat qxh = 0;
        PMEFloat charge = cSim.pcharge_phmd[atmi];
        PMEFloat sqrt_pi = sqrt(PI);
        PMEFloat factor = -0.5 * PI / (ewaldcof * ewaldcof);
        PMEFloat qprot, qunprot, x2, lambda;
        int psp_grp = cSim.psp_grp[h];
        PMEFloat pph_theta = cSim.pph_theta[h];
        PMEFloat pluspph_theta = cSim.pph_theta[h+1];
        PMEFloat2 pqstate1 = cSim.pImageQstate1[atmx];
        PMEFloat2 pqstate2 = cSim.pImageQstate2[atmx];
        lambda = sin(pph_theta) * sin(pph_theta);
        x2 = 1.0;
        if (psp_grp > 0) {
          x2 = sin(pluspph_theta) * sin(pluspph_theta);
          qunprot = lambda * (pqstate2.x - pqstate2.y);
          qprot = (1 - lambda) * (pqstate1.x - pqstate1.y);
          qxh = qunprot + qprot;
          dudlplus += 2.0 * factor * sumq * qxh / vol;
          dudlplus -= 2.0 * charge * ewaldcof * qxh / sqrt_pi;
        }
        qunprot = x2 * pqstate2.x + (1 - x2) * pqstate2.y;
        qprot = x2 * pqstate1.x + (1 - x2) * pqstate1.y;
        qh = qunprot - qprot;
        dudl += 2.0 * factor * sumq * qh / vol;
        dudl -= 2.0 * charge * ewaldcof * qh / sqrt_pi;
#ifdef use_DPFP
        atomicAdd((unsigned long long int*)&cSim.pdph_accumulator[atmx],
                   llitoulli(llrint((PMEDouble)dudl * FORCESCALE)));
        if (psp_grp > 0) {
          atomicAdd((unsigned long long int*)&cSim.pdph_plus_accumulator[atmx],
                  llitoulli(llrint((PMEDouble)dudlplus * FORCESCALE)));
        }
#else
        atomicAdd((unsigned long long int*)&cSim.pdph_accumulator[atmx],
                llitoulli(fast_llrintf(FORCESCALEF * dudl)));
        if (psp_grp > 0) {
          atomicAdd((unsigned long long int*)&cSim.pdph_plus_accumulator[atmx],
                  llitoulli(fast_llrintf(FORCESCALEF * dudlplus)));
        }
#endif
      }
    }
  }
}

//---------------------------------------------------------------------------------------------
// kCalculateSelfPHMD: launch the kernel to compute the Self titrtion forces
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kCalculateSelfPHMD(double vol, double ewaldcof, gpuContext gpu)
{
  kCalculateSelfPHMD_kernel<<<gpu->PMENonbondBlocks, gpu->PMENonbondForcesThreadsPerBlock>>>(vol, ewaldcof,
                                                                                            gpu->sumq,
                                                                                            gpu->sumq2);
}
