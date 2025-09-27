#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------
#if defined(AMBER_PLATFORM_AMD)
#  include <hip/hip_runtime.h>
#  include <hip/hip_runtime_api.h>
#  include "hip_definitions.h"
#else
#  include <cuda.h>
#endif
#include "gpu.h"
#include "mdl_xray.h"
#include "ptxmacros.h"

// Use global instance instead of a local copy
#include "xrayDevConstants.h"
#include "simulationConst.h"
CSIM_STO simulationConst cSim;
__device__ __constant__ xrayDevConstants cXray;

// Constants for this module
GL_CONST PMEFloat TWOPI_F                 = (PMEFloat)6.283185482025146484375;
GL_CONST int XRAY_WORK_THREADS_PER_BLOCK  = 256;
GL_CONST int XRAY_WORK_BLOCKS_MULTIPLIER  = 6;
GL_CONST PMEFloat F_EPSILON               = (PMEFloat)1.0e-20;

#if !defined(__HIPCC_RDC__)

//---------------------------------------------------------------------------------------------
// SetkXrayUpdateSim: function to upload constants and pointers for the global instance of the
//                    simulationConst within this CUDA unit.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void SetkXrayUpdateSim(gpuContext gpu)
{
  cudaError_t status;
  status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(simulationConst));
  RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

#endif

//---------------------------------------------------------------------------------------------
// SetkXrayUpdateCalc: function to upload data to the GPU as constant memory through the
//                     xrayDevConstants symbol.
//
// Arguments:
//   xrd: set of constants and pointers for X-ray structure factor computations on the GPU
//---------------------------------------------------------------------------------------------
void SetkXrayUpdateCalc(xrayHostContext xrd)
{
  cudaError_t status;
  status = cudaMemcpyToSymbol(cXray, &xrd->devInfo, sizeof(xrayDevConstants));
  RTERROR(status, "cudaMemcpyToSymbol: Copying X-ray calculation bundle to devInfo failed");
}

//---------------------------------------------------------------------------------------------
// kXrayInitResults: kernel for initializing global accumulators needed by the X-ray structure
//                   factor calculations.
//---------------------------------------------------------------------------------------------
__global__ void kXrayInitResults_kernel()
{
  if (threadIdx.x < 32 && blockIdx.x == 0) {

    // This is pointer-array abuse.  pSumFoFc points to
    // the first index of an array of GRID indices.
    cXray.pSumFoFc[threadIdx.x] = (PMEAccumulator)0;
  }
  int pos = threadIdx.x + (blockIdx.x * blockDim.x);
  int pdatoms = 3 * ((cSim.atoms + GRID_BITS_MASK) / GRID) * GRID;
  while (pos < pdatoms) {
    cXray.pForceX[pos] = (PMEAccumulator)0.0;
    pos += blockDim.x * gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kXrayComputeMSS4: GPU kernel for computing the value of an intermediate quantity needed by
//                   structure factor calculations.  Most of this happens in double precision
//                   because this kernel is not called frequently.
//---------------------------------------------------------------------------------------------
__global__ void kXrayComputeMSS4_kernel()
{
  PMEDouble sina = sin(cSim.alpha);
  PMEDouble cosa = cos(cSim.alpha);
  PMEDouble sinb = sin(cSim.beta);
  PMEDouble cosb = cos(cSim.beta);
  PMEDouble sing = sin(cSim.gamma);
  PMEDouble cosg = cos(cSim.gamma);
  PMEDouble V    = (cSim.a * cSim.b * cSim.c) *
                   sqrt(1.0 - cosa*cosa - cosb*cosb - cosg*cosg + 2.0*cosa*cosb*cosg);
  PMEDouble astar = cSim.b * cSim.c * sina / V;
  PMEDouble bstar = cSim.a * cSim.c * sinb / V;
  PMEDouble cstar = cSim.a * cSim.b * sing / V;
  PMEDouble cosas = (cosb*cosg - cosa)/(sinb*sing);
  PMEDouble cosbs = (cosa*cosg - cosb)/(sina*sing);
  PMEDouble cosgs = (cosb*cosa - cosg)/(sinb*sina);
  
  int tidx = threadIdx.x + (blockIdx.x * blockDim.x);
  while (tidx < cXray.nHKL) {
    PMEDouble dh = cXray.pHKLaidx[tidx];
    PMEDouble dk = cXray.pHKLbidx[tidx];
    PMEDouble dl = cXray.pHKLcidx[tidx];
    PMEDouble S2 = (dh * dh * astar * astar) + (dk * dk * bstar * bstar) +
                   (dl * dl * cstar * cstar) +
                   2.0*((dk * dl * bstar * cstar * cosas) + (dh * dl * astar * cstar * cosbs) +
                        (dh * dk * astar * bstar * cosgs));
    cXray.pMSS4[tidx] = (PMEFloat)-0.25 * (PMEFloat)S2;
    tidx += gridDim.x * blockDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kXrayGetDerivative1_kernel: first kernel to compute the derivative of X-ray based structure
//                             factor restraints.  This mirrors xray_get_derivative() in
//                             xray_interface.F90.  This kernel MUST be launched with
//                             XRAY_WORK_THREADS_PER_BLOCK threads.
//---------------------------------------------------------------------------------------------
__global__ void kXrayGetDerivative1_kernel()
{
  int i;
  __shared__ volatile int nReflBatch, nAtomBatch;
  __shared__ volatile int scType[XRAY_WORK_THREADS_PER_BLOCK], sHKLMask[GRID];
  __shared__ volatile PMEFloat FoFcAcc[GRID], FcFcAcc[GRID], FobsAcc[GRID];
  __shared__ volatile PMEFloat sHKLa[GRID], sHKLb[GRID], sHKLc[GRID], sMSS4[GRID];
  __shared__ volatile PMEFloat sfcreal[GRID], sfcimag[GRID], sAbsFObs[GRID];
  __shared__ volatile PMEFloat crdx[XRAY_WORK_THREADS_PER_BLOCK];
  __shared__ volatile PMEFloat crdy[XRAY_WORK_THREADS_PER_BLOCK];
  __shared__ volatile PMEFloat crdz[XRAY_WORK_THREADS_PER_BLOCK];
  __shared__ volatile PMEFloat sBfactor[XRAY_WORK_THREADS_PER_BLOCK];
  __shared__ volatile PMEFloat soccp[XRAY_WORK_THREADS_PER_BLOCK];
  __shared__ volatile PMEFloat sSCfac[16*GRID];
  __shared__ volatile PMEFloat2 sFObs[GRID];
  __shared__ volatile int tileCounter;

  // Customary indexing
 unsigned int tgx     = (threadIdx.x & GRID_BITS_MASK);
  int warpIdx = threadIdx.x / GRID;
  
  // Initialize accumulators (only one warp will use these, but it's OK)
  PMEFloat sumFoFc = (PMEFloat)0.0;
  PMEFloat sumFcFc = (PMEFloat)0.0;
  //PMEFloat sumFobs = (PMEFloat)0.0;
  
  // Take in batches of 32 reflections at a time, store them in __shared__,
  // and deal with them one by one over batches of atoms
  int reflBase = blockIdx.x * GRID;
  while (reflBase < cXray.nHKL) {

    // Read HKL indices and reflection details for this tile
    int reflIdx = reflBase + tgx;
    if (reflIdx < cXray.nHKL) {
      if (warpIdx == 0) {
        sMSS4[tgx]     = cXray.pMSS4[reflIdx];
      }
      else if (warpIdx == 1) {
        sHKLa[tgx]     = TWOPI_F * (PMEFloat)(cXray.pHKLaidx[reflIdx]);
      }
      else if (warpIdx == 2) {
        sHKLb[tgx]     = TWOPI_F * (PMEFloat)(cXray.pHKLbidx[reflIdx]);
      }
      else if (warpIdx == 3) {
        sHKLc[tgx]     = TWOPI_F * (PMEFloat)(cXray.pHKLcidx[reflIdx]);
      }
      else if (warpIdx == 4) {
        sAbsFObs[tgx]  = cXray.pAbsFObs[reflIdx];
      }
      else if (warpIdx == 5) {

	// Initialize the sturcture factor calculations with
	// the solvent mask contribution, if there is any
	if (cXray.bsMaskModel == 0) {
          sfcreal[tgx]   = (PMEFloat)0.0;
          sfcimag[tgx]   = (PMEFloat)0.0;
	}
	else if (cXray.bsMaskModel == 1) {
	  PMEFloat  skm = cXray.pSolventKMask[reflIdx];
	  PMEFloat2 sfm = cXray.pSolventFMask[reflIdx];
          sfcreal[tgx]   = skm * sfm.x;
          sfcimag[tgx]   = skm * sfm.y;
	}
      }
      else if (warpIdx == 6) {
        sHKLMask[tgx]  = cXray.pHKLMask[reflIdx];
      }
      else if (warpIdx == 7) {
	PMEFloat2 xfer = cXray.pFObs[reflIdx];
	sFObs[tgx].x   = xfer.x;
	sFObs[tgx].y   = xfer.y;
      }
    }
    if (threadIdx.x == 0) {
      nReflBatch = (cXray.nHKL - reflBase >= GRID) ? GRID : cXray.nHKL - reflBase;
    }
    __syncthreads();

    // Pre-compute and tabulate the atomic scattering factors for these reflections
    int sccon = warpIdx;
    while (sccon < cXray.nScatterTypes) {
      int scoffset = 2 * sccon * cXray.nScatterCoeff;
      PMEFloat atmscfac = cXray.pScatterCoeffs[scoffset + cXray.nScatterCoeff - 1];
      for (i = 0; i < cXray.nScatterCoeff - 1; i++) {
        atmscfac += cXray.pScatterCoeffs[scoffset + i] *
                    exp(sMSS4[tgx] * cXray.pScatterCoeffs[scoffset + cXray.nScatterCoeff + i]);
      }
      sSCfac[sccon*GRID + tgx] = atmscfac;
      sccon += blockDim.x / GRID;
    }
    __syncthreads();
    
    // Loop over all atoms
    int atomBase = 0;
    while (atomBase < cXray.nSelectedAtom) {
      if (atomBase + threadIdx.x < cXray.nSelectedAtom) {
        int atmidx = cXray.pSelectedAtoms[atomBase + threadIdx.x];
        int imgidx = cSim.pImageAtomLookup[atmidx];
        scType[threadIdx.x] = cXray.pScatterType[atmidx];
        soccp[threadIdx.x]  = cXray.pOccupancy[atmidx];
        sBfactor[threadIdx.x] = cXray.pAtomBFactor[atmidx];
        PMEFloat fracx  = (cSim.pFractX[imgidx] / cSim.nfft1) + (PMEFloat)0.5;
        PMEFloat fracy  = (cSim.pFractY[imgidx] / cSim.nfft2) + (PMEFloat)0.5;
        PMEFloat fracz  = (cSim.pFractZ[imgidx] / cSim.nfft3) + (PMEFloat)0.5;
        fracx -= (fracx > (PMEFloat)1.0);
        fracy -= (fracy > (PMEFloat)1.0);
        fracz -= (fracz > (PMEFloat)1.0);
        crdx[threadIdx.x] = fracx;
        crdy[threadIdx.x] = fracy;
        crdz[threadIdx.x] = fracz;
      }
      if (threadIdx.x == 0) {
        int nremain = cXray.nSelectedAtom - atomBase;
        nAtomBatch = (nremain >= blockDim.x) ? blockDim.x : nremain;
      }
      __syncthreads();
      
      // Perform work
      int reflpos = warpIdx;
      while (reflpos < nReflBatch) {

        // Skip if this reflection has been masked out
#if 0
        if (sHKLMask[reflpos] == 0) {
          reflpos += blockDim.x / GRID;
          continue;
        }
#endif
        // Get the details for this reflection from __shared__ L1 cache
        PMEFloat tmss4 = sMSS4[reflpos];
        PMEFloat hkla  = sHKLa[reflpos];
        PMEFloat hklb  = sHKLb[reflpos];
        PMEFloat hklc  = sHKLc[reflpos];
        PMEFloat fcreal = (PMEFloat)0.0;
        PMEFloat fcimag = (PMEFloat)0.0;
        int atomBase2 = 0;
        while (atomBase2 < nAtomBatch) {
          int atmb2x = atomBase2 + tgx;
          PMEFloat toccp, angle, atmscfac;
          if (atmb2x < nAtomBatch) {
            atmscfac = sSCfac[(scType[atmb2x] * GRID) + reflpos];
            toccp  = soccp[atmb2x];
            angle  = hkla*crdx[atmb2x] + hklb*crdy[atmb2x] + hklc*crdz[atmb2x];
          }
          else {
            atmscfac = (PMEFloat)0.0;
            toccp    = (PMEFloat)0.0;
            angle    = (PMEFloat)0.0;
          }
          PMEFloat f = exp(tmss4 * sBfactor[atmb2x]) * toccp * atmscfac;
          fcreal    += f * cosf(angle);
          fcimag    += f * sinf(angle);

          // Increment the inner loop atom counter
          atomBase2 += GRID;
        }

        // Reduction of the structure factors for this reflection and batch of atoms
        warp_reduction_down(fcreal);
        warp_reduction_down(fcimag);
        if (tgx == 0) {
          sfcreal[reflpos]         += fcreal;
          sfcimag[reflpos]         += fcimag;
        }

        // Increment the reflection counter
        reflpos += blockDim.x / GRID;
      }
      __syncthreads();

      // Increment the atom window
      atomBase += blockDim.x;
    }
    __syncthreads();
    
    // Make some computations for the first stage of some massive
    // reductions that will be needed in the next kernel.
    if (threadIdx.x < nReflBatch) {
      PMEFloat2 fc;
      fc.x                      = sfcreal[threadIdx.x];
      fc.y                      = sfcimag[threadIdx.x];
      cXray.pFcalc[reflIdx]     = fc;
      PMEFloat absFcalc         = sqrt(fc.x*fc.x + fc.y*fc.y);
      cXray.pAbsFcalc[reflIdx]  = absFcalc;
      PMEFloat tabsfobs         = sAbsFObs[threadIdx.x];
      PMEFloat2 tfobs;
      tfobs.x                   = sFObs[threadIdx.x].x;
      tfobs.y                   = sFObs[threadIdx.x].y;

      // The sums of FoFc and FcFc denote quantities related to Rwork, not Rfree
      if (sHKLMask[threadIdx.x] != 0) {
        if (cXray.tarfunc == 0) {
          sumFoFc                += tabsfobs * absFcalc;
        }
        else {
          sumFoFc                += (tfobs.x * fc.x) + (tfobs.y * fc.y);
        }
        sumFcFc                  += absFcalc * absFcalc;
      }
    }
    __syncthreads();

    // Increment the reflection base
    reflBase += gridDim.x * GRID;
  }
  
  // Reduce the FObs and FCalc products across the warp, then the block
  if (threadIdx.x < GRID) {
    warp_reduction_down(sumFoFc);
    warp_reduction_down(sumFcFc);
    if (threadIdx.x == 0) {
      atomicAdd((unsigned long long int*)&cXray.pSumFoFc[0],
                llitoulli(llrint(sumFoFc * XRAYSCALEF)));
      atomicAdd((unsigned long long int*)&cXray.pSumFcFc[0],
                llitoulli(llrint(sumFcFc * XRAYSCALEF)));
    }
  }
}

//---------------------------------------------------------------------------------------------
// kXrayGetDerivative2_kernel: second kernel for obtaining forces from X-ray restraints.  This
//                             covers the operations in the dTargetML_dF and dTargetV_dF
//                             subroutine of xray_fourier.F90.
//---------------------------------------------------------------------------------------------
__global__ void kXrayGetDerivative2_kernel(int calcScalingValues)
{
  __shared__ volatile PMEFloat RsdAcc[GRID], RsdFreeAcc[GRID], NrgAcc[GRID];
  __shared__ volatile PMEFloat FcalcScale;

  int pos = threadIdx.x + (blockIdx.x * blockDim.x);
  int warpIdx = threadIdx.x / GRID;
 unsigned int tgx = (threadIdx.x & GRID_BITS_MASK);
  if (threadIdx.x < GRID) {
    RsdAcc[threadIdx.x]     = (PMEFloat)0.0;
    RsdFreeAcc[threadIdx.x] = (PMEFloat)0.0;
    NrgAcc[threadIdx.x]     = (PMEFloat)0.0;
  }
  
  // Get critical scaling factors
  if (threadIdx.x == GRID) {
    if (calcScalingValues == 1) {
      PMEAccumulator sumFoFc = *cXray.pSumFoFc;
      PMEAccumulator sumFcFc = *cXray.pSumFcFc;
      PMEFloat fsumFoFc = (PMEDouble)sumFoFc * ONEOVERXRAYSCALE;
      PMEFloat fsumFcFc = (PMEDouble)sumFcFc * ONEOVERXRAYSCALE;
      FcalcScale = fsumFoFc / fsumFcFc;
      if (blockIdx.x == 0) {
        cXray.pSFScaleFactors[0] = FcalcScale;
      }
    }
    else {
      FcalcScale = cXray.pSFScaleFactors[0];
    }
  }
  __syncthreads();

  PMEFloat rsdsum = (PMEFloat)0.0;
  PMEFloat rsdfreesum = (PMEFloat)0.0;
  PMEFloat nrgsum = (PMEFloat)0.0;
  while (pos < cXray.nHKL) {
    PMEFloat2 drv;
    if (cXray.pHKLMask[pos] == 0) {
      drv.x = (PMEFloat)0.0;
      drv.y = (PMEFloat)0.0;
      PMEFloat tfobs = cXray.pAbsFObs[pos];
      if (cXray.tarfunc == 0) {
        rsdfreesum += fabs(tfobs - (FcalcScale * cXray.pAbsFcalc[pos]));
      }
      else if (cXray.tarfunc == 1) {
	PMEFloat2 tvfobs = cXray.pFObs[pos];
	if (cXray.pAbsFcalc[pos] > F_EPSILON) {
          PMEFloat2 fc  = cXray.pFcalc[pos];
	  PMEFloat vdiffx = tvfobs.x - (FcalcScale * fc.x);
	  PMEFloat vdiffy = tvfobs.y - (FcalcScale * fc.y);
	  rsdfreesum += sqrt(vdiffx*vdiffx + vdiffy*vdiffy);
	}
      }
    }
    else {
      PMEFloat tfobs = cXray.pAbsFObs[pos];
      PMEFloat afc = cXray.pAbsFcalc[pos];
      if (cXray.tarfunc == 0) {
        if (afc > F_EPSILON) {
          PMEFloat2 fc  = cXray.pFcalc[pos];
          PMEFloat ascl = (PMEFloat)-2.0 * cXray.NormScale * (tfobs - FcalcScale*afc) *
                          (FcalcScale / afc);
          drv.x = ascl * fc.x * cXray.xrayGeneralWt;
          drv.y = ascl * fc.y * cXray.xrayGeneralWt;
        }
        else {
          drv.x = (PMEFloat)0.0;
          drv.y = (PMEFloat)0.0;
        }
        rsdsum += fabs(tfobs - (FcalcScale * afc));
        nrgsum += cXray.NormScale * (tfobs - FcalcScale*afc) * (tfobs - FcalcScale*afc);
      }
      else if (cXray.tarfunc == 1) {
	PMEFloat2 vdiff;
	PMEFloat2 tvfobs = cXray.pFObs[pos];
	if (afc > F_EPSILON) {
          PMEFloat2 fc  = cXray.pFcalc[pos];
	  vdiff.x = tvfobs.x - (FcalcScale * fc.x);
	  vdiff.y = tvfobs.y - (FcalcScale * fc.y);
	  drv.x = -(PMEFloat)2.0 * cXray.NormScale * FcalcScale * vdiff.x;
	  drv.y = -(PMEFloat)2.0 * cXray.NormScale * FcalcScale * vdiff.y;
	}
	else {
	  vdiff.x = (PMEFloat)0.0;
	  vdiff.y = (PMEFloat)0.0;
	  drv.x = (PMEFloat)0.0;
          drv.y = (PMEFloat)0.0;
        }
	rsdsum += sqrt(vdiff.x*vdiff.x + vdiff.y*vdiff.y);
	nrgsum += cXray.NormScale * (vdiff.x*vdiff.x + vdiff.y*vdiff.y);
      }
    }
    cXray.pReflDeriv[pos] = drv;
    pos += (blockDim.x * gridDim.x);
  }
  warp_reduction_down(rsdsum);
  warp_reduction_down(rsdfreesum);
  warp_reduction_down(nrgsum);
  if (tgx == 0) {
    RsdAcc[warpIdx] = rsdsum;
    RsdFreeAcc[warpIdx] = rsdfreesum;
    NrgAcc[warpIdx] = nrgsum * cXray.xrayGeneralWt;
  }
  __syncthreads();
  if (warpIdx == 0) {
    rsdsum = RsdAcc[tgx];
    warp_reduction_down(rsdsum);
    if (tgx == 0) {
      if (cXray.SumFObsWork > (PMEFloat)1.0e-8) {
        atomicAdd((unsigned long long int*)&cXray.pResidual[0],
                  llitoulli(llrint(rsdsum / cXray.SumFObsWork * ENERGYSCALEF)));
      }
    }
  }
  else if (warpIdx == 1) {
    nrgsum = NrgAcc[tgx];
    warp_reduction_down(nrgsum);
    if (tgx == 0) {
      atomicAdd((unsigned long long int*)&cXray.pXrayEnergy[0],
                llitoulli(llrint(nrgsum * ENERGYSCALEF)));
    }
  }
  else if (warpIdx == 2) {
    rsdfreesum = RsdFreeAcc[tgx];
    warp_reduction_down(rsdfreesum);
    if (tgx == 0) {
      if (cXray.SumFObsFree > (PMEFloat)1.0e-8) {
        atomicAdd((unsigned long long int*)&cXray.pFreeResidual[0],
                  llitoulli(llrint(rsdfreesum / cXray.SumFObsFree * ENERGYSCALEF)));
      }
    }
  }
}

//---------------------------------------------------------------------------------------------
// kXrayGetDerivative3_kernel: third kernel for obtaining forces from X-ray restraints.  This
//                             covers the operations of the subroutine fourier_dXYZBQ_dF from
//                             xray_fourier.F90.
//---------------------------------------------------------------------------------------------
__global__ void kXrayGetDerivative3_kernel()
{
  int i;
  __shared__ volatile int nReflBatch;
  __shared__ volatile int sHKLMask[GRID];
  __shared__ volatile PMEFloat FoFcAcc[GRID], FcFcAcc[GRID];
  __shared__ volatile PMEFloat sHKLa[GRID], sHKLb[GRID], sHKLc[GRID], sMSS4[GRID];
  __shared__ volatile PMEFloat sSCfac[16*GRID];
  __shared__ volatile PMEFloat2 sdrv[GRID];
  __shared__ volatile int tileCounter;

  // Customary indexing
 unsigned int tgx     = (threadIdx.x & GRID_BITS_MASK);
  int warpIdx = threadIdx.x / GRID;
  
  // Take in batches of 32 reflections at a time, store them in __shared__,
  // and deal with them one by one over batches of atoms
  int reflBase = blockIdx.x * GRID;
  while (reflBase < cXray.nHKL) {

    // Read HKL indices and reflection details for this tile
    int reflIdx = reflBase + tgx;
    if (reflIdx < cXray.nHKL) {
      if (warpIdx == 0) {
        sMSS4[tgx]  = cXray.pMSS4[reflIdx];
      }
      else if (warpIdx == 1) {
        sHKLa[tgx]  = TWOPI_F * (PMEFloat)(cXray.pHKLaidx[reflIdx]);
      }
      else if (warpIdx == 2) {
        sHKLb[tgx]  = TWOPI_F * (PMEFloat)(cXray.pHKLbidx[reflIdx]);
      }
      else if (warpIdx == 3) {
        sHKLc[tgx]  = TWOPI_F * (PMEFloat)(cXray.pHKLcidx[reflIdx]);
      }
      else if (warpIdx == 4) {
        sHKLMask[tgx] = cXray.pHKLMask[reflIdx];
      }
      else if (warpIdx == 5) {
        PMEFloat2 drv = cXray.pReflDeriv[reflIdx];
        sdrv[tgx].x = drv.x;
        sdrv[tgx].y = drv.y;
      }
    }
    if (threadIdx.x == 0) {
      nReflBatch = (cXray.nHKL - reflBase >= GRID) ? GRID : cXray.nHKL - reflBase;
    }
    __syncthreads();

    // Pre-compute and tabulate the atomic scattering factors for these reflections
    int sccon = warpIdx;
    while (sccon < cXray.nScatterTypes) {
      int scoffset = 2 * sccon * cXray.nScatterCoeff;
      PMEFloat atmscfac = cXray.pScatterCoeffs[scoffset + cXray.nScatterCoeff - 1];
      for (i = 0; i < cXray.nScatterCoeff - 1; i++) {
        atmscfac += cXray.pScatterCoeffs[scoffset + i] *
                    exp(sMSS4[tgx] * cXray.pScatterCoeffs[scoffset + cXray.nScatterCoeff + i]);
      }
      sSCfac[sccon*GRID + tgx] = atmscfac;
      sccon += blockDim.x / GRID;
    }
    __syncthreads();

    // Loop over all atoms in strides
    int atomBase = 0;
    while (atomBase < cXray.nSelectedAtom) {
      if (atomBase + threadIdx.x >= cXray.nSelectedAtom) {
        atomBase += blockDim.x;
        continue;
      }

      // Each thread takes one atom and loops over all reflections
      int atmidx        = cXray.pSelectedAtoms[atomBase + threadIdx.x];
      int imgidx        = cSim.pImageAtomLookup[atmidx];
      int scType        = cXray.pScatterType[atmidx];
      PMEFloat toccp    = cXray.pOccupancy[atmidx];
      PMEFloat tBfactor = cXray.pAtomBFactor[atmidx];
      PMEFloat fracx    = (cSim.pFractX[imgidx] / cSim.nfft1) + (PMEFloat)0.5;
      PMEFloat fracy    = (cSim.pFractY[imgidx] / cSim.nfft2) + (PMEFloat)0.5;
      PMEFloat fracz    = (cSim.pFractZ[imgidx] / cSim.nfft3) + (PMEFloat)0.5;
      fracx -= (fracx > (PMEFloat)1.0);
      fracy -= (fracy > (PMEFloat)1.0);
      fracz -= (fracz > (PMEFloat)1.0);
      PMEFloat crdx = fracx;
      PMEFloat crdy = fracy;
      PMEFloat crdz = fracz;
      PMEFloat frcx = (PMEFloat)0.0;
      PMEFloat frcy = (PMEFloat)0.0;
      PMEFloat frcz = (PMEFloat)0.0;
      int reflpos = 0;
      while (reflpos < nReflBatch) {

        // Skip if this reflection has been masked out
        if (sHKLMask[reflpos] == 0) {
          reflpos++;
          continue;
        }
        
        // Get the details for this reflection from __shared__ L1 cache
        PMEFloat tmss4 = sMSS4[reflpos];
        PMEFloat hkla  = sHKLa[reflpos];
        PMEFloat hklb  = sHKLb[reflpos];
        PMEFloat hklc  = sHKLc[reflpos];
        PMEFloat atmscfac = sSCfac[scType*GRID + reflpos];
        PMEFloat angle    = hkla*crdx + hklb*crdy + hklc*crdz;
        PMEFloat f        = exp(tmss4 * tBfactor) * toccp * atmscfac;
        PMEFloat fcreal   = f * cosf(angle);
        PMEFloat fcimag   = f * sinf(angle);
        PMEFloat cmpfac   = (fcimag * sdrv[reflpos].x) - (fcreal * sdrv[reflpos].y);
        frcx -= hkla * cmpfac;
        frcy -= hklb * cmpfac;
        frcz -= hklc * cmpfac;

        // Increment the reflection counter
        reflpos++;
      }

      // Write results back to global
#ifdef use_DPFP      
      atomicAdd((unsigned long long int*)&cXray.pForceX[atmidx],
                llitoulli(llrint(frcx * FORCESCALE)));
      atomicAdd((unsigned long long int*)&cXray.pForceY[atmidx],
                llitoulli(llrint(frcy * FORCESCALE)));
      atomicAdd((unsigned long long int*)&cXray.pForceZ[atmidx],
                llitoulli(llrint(frcz * FORCESCALE)));
#else
      atomicAdd((unsigned long long int*)&cXray.pForceX[atmidx],
                llitoulli(llrint(frcx * FORCESCALEF)));
      atomicAdd((unsigned long long int*)&cXray.pForceY[atmidx],
                llitoulli(llrint(frcy * FORCESCALEF)));
      atomicAdd((unsigned long long int*)&cXray.pForceZ[atmidx],
                llitoulli(llrint(frcz * FORCESCALEF)));
#endif
      // Increment the atom window
      atomBase += blockDim.x;
    }
    __syncthreads();

    // Increment the reflection base
    reflBase += gridDim.x * GRID;
  }
}

//---------------------------------------------------------------------------------------------
// kXrayGetDerivative4_kernel: final kernel to convert X-ray based forces back to cartesian
//                             coordinates and then add them to the non-bonded forces.
//---------------------------------------------------------------------------------------------
__global__ void kXrayGetDerivative4_kernel()
{
  __shared__ volatile PMEFloat sRecipf[9];
  
  int basepos = blockDim.x * blockIdx.x;

  if (cSim.ntp > 0 && cSim.barostat == 1) {
    if (threadIdx.x < 9) {
      sRecipf[threadIdx.x] = cSim.pNTPData->recipf[threadIdx.x];
    }
    __syncthreads();

    while (basepos < cSim.atoms) {

      // Check that this is a selected atom
      int ipos = basepos + threadIdx.x;
      if (ipos < cSim.atoms) {
        int addxray = cXray.pSelectionMask[ipos];

        // Read forces in direct order, then change coordinates
        // write back with scatter to the non-bonded forces
        if (addxray == 1) {
#ifdef use_DPFP
          PMEFloat sfx = (PMEDouble)cXray.pForceX[ipos] * ONEOVERFORCESCALE;
          PMEFloat sfy = (PMEDouble)cXray.pForceY[ipos] * ONEOVERFORCESCALE;
          PMEFloat sfz = (PMEDouble)cXray.pForceZ[ipos] * ONEOVERFORCESCALE;
#else
          PMEFloat sfx = (PMEDouble)cXray.pForceX[ipos] * ONEOVERFORCESCALEF;
          PMEFloat sfy = (PMEDouble)cXray.pForceY[ipos] * ONEOVERFORCESCALEF;
          PMEFloat sfz = (PMEDouble)cXray.pForceZ[ipos] * ONEOVERFORCESCALEF;
#endif
          PMEFloat nx  = sRecipf[0]*sfx;
          PMEFloat ny  = sRecipf[3]*sfx + sRecipf[4]*sfy;
          PMEFloat nz  = sRecipf[6]*sfx + sRecipf[7]*sfy + sRecipf[8]*sfz;        
          int imgpos = cSim.pImageAtomLookup[ipos];
#ifdef use_DPFP
          atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[imgpos],
                    llitoulli(llrint(-nx * FORCESCALE)));
          atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[imgpos],
                    llitoulli(llrint(-ny * FORCESCALE)));
          atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[imgpos],
                    llitoulli(llrint(-nz * FORCESCALE)));
#else
          atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[imgpos],
                    llitoulli(llrint(-nx * FORCESCALEF)));
          atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[imgpos],
                    llitoulli(llrint(-ny * FORCESCALEF)));
          atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[imgpos],
                    llitoulli(llrint(-nz * FORCESCALEF)));
#endif
          cXray.pForceX[ipos] = nx;
          cXray.pForceY[ipos] = ny;
          cXray.pForceZ[ipos] = nz;
        }
      }
      basepos += gridDim.x * blockDim.x;
    }
  }
  else {
    while (basepos < cSim.atoms) {

      // Check that this is a selected atom
      int ipos = basepos + threadIdx.x;
      if (ipos < cSim.atoms) {
        int addxray = cXray.pSelectionMask[ipos];

        // Read forces in direct order, then change coordinates,
        // write back with scatter to the non-bonded forces
        if (addxray == 1) {
#ifdef use_DPFP
          PMEFloat sfx = (PMEDouble)cXray.pForceX[ipos] * ONEOVERFORCESCALE;
          PMEFloat sfy = (PMEDouble)cXray.pForceY[ipos] * ONEOVERFORCESCALE;
          PMEFloat sfz = (PMEDouble)cXray.pForceZ[ipos] * ONEOVERFORCESCALE;
#else
          PMEFloat sfx = (PMEDouble)cXray.pForceX[ipos] * ONEOVERFORCESCALEF;
          PMEFloat sfy = (PMEDouble)cXray.pForceY[ipos] * ONEOVERFORCESCALEF;
          PMEFloat sfz = (PMEDouble)cXray.pForceZ[ipos] * ONEOVERFORCESCALEF;
#endif
          PMEFloat nx  = cSim.recipf[0][0]*sfx;
          PMEFloat ny  = cSim.recipf[1][0]*sfx + cSim.recipf[1][1]*sfy;
          PMEFloat nz  = cSim.recipf[2][0]*sfx + cSim.recipf[2][1]*sfy + cSim.recipf[2][2]*sfz;
          int imgpos = cSim.pImageAtomLookup[ipos];
#ifdef use_DPFP
          atomicAdd((unsigned long long int*)&cSim.pForceXAccumulator[imgpos],
                    llitoulli(llrint(-nx * FORCESCALE)));
          atomicAdd((unsigned long long int*)&cSim.pForceYAccumulator[imgpos],
                    llitoulli(llrint(-ny * FORCESCALE)));
          atomicAdd((unsigned long long int*)&cSim.pForceZAccumulator[imgpos],
                    llitoulli(llrint(-nz * FORCESCALE)));
#else
          atomicAdd((unsigned long long int*)&cSim.pForceXAccumulator[imgpos],
                    llitoulli(llrint(-nx * FORCESCALEF)));
          atomicAdd((unsigned long long int*)&cSim.pForceYAccumulator[imgpos],
                    llitoulli(llrint(-ny * FORCESCALEF)));
          atomicAdd((unsigned long long int*)&cSim.pForceZAccumulator[imgpos],
                    llitoulli(llrint(-nz * FORCESCALEF)));
#endif
          cXray.pForceX[ipos] = nx;
          cXray.pForceY[ipos] = ny;
          cXray.pForceZ[ipos] = nz;
        }
      }
      basepos += gridDim.x * blockDim.x;
    }
  }

}

//---------------------------------------------------------------------------------------------
// kXrayGetDerivative: function to launch the kernel for X-ray restraint force computations.
//
// Arguments:
//---------------------------------------------------------------------------------------------
extern "C" void kXrayGetDerivative(gpuContext gpu, xrayHostContext xrd, int calcScalingValues)
{
  // Set tile size in local variables
  int nBlocks = XRAY_WORK_BLOCKS_MULTIPLIER * gpu->blocks;

  // Decision on whether to recompute MSS4
  kXrayComputeMSS4_kernel<<<nBlocks, XRAY_WORK_THREADS_PER_BLOCK>>>();  
  kXrayInitResults_kernel<<<nBlocks, XRAY_WORK_THREADS_PER_BLOCK>>>();  
  kXrayGetDerivative1_kernel<<<nBlocks, XRAY_WORK_THREADS_PER_BLOCK>>>();
  kXrayGetDerivative2_kernel<<<nBlocks, XRAY_WORK_THREADS_PER_BLOCK>>>(calcScalingValues);
  kXrayGetDerivative3_kernel<<<nBlocks, XRAY_WORK_THREADS_PER_BLOCK>>>();
  kXrayGetDerivative4_kernel<<<nBlocks, XRAY_WORK_THREADS_PER_BLOCK>>>();
}
