#ifdef GTI
#ifndef AMBER_PLATFORM_AMD
#include <cuda.h>
#endif
#include "gputypes.h"
#include "ptxmacros.h"
#include "gti_const.cuh"
#include "gti_PME_kernels.cuh"

#include "gti_utils.cuh"

#include "simulationConst.h"
CSIM_STO simulationConst cSim;

namespace GTI_PME_IMPL {
#if !(defined(__CUDACC_RDC__) || defined(__HIPCC_RDC__))
#  include "gti_localCUutils.inc"
#endif
}
using namespace GTI_PME_IMPL;

//---------------------------------------------------------------------------------------------
// kgTIPMEFillChargeGrid_kernel:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgTIPMEFillChargeGrid_kernel()
{
  using namespace GTI_PME_IMPL;

  __shared__ volatile  FillChargeGridAtomData sAtom[FillY][3][LOADSIZE];
  __shared__ volatile unsigned int region[FillY][LOADSIZE];

  // Determine grid offsets
  const int iOffsetX = threadIdx.x & 0x03;
  const int iOffsetY = (threadIdx.x & 0x0f) >> 2;
  const int iOffsetZ = threadIdx.x >> 4;
  unsigned int cell = (blockIdx.x * blockDim.y) + threadIdx.y;

  while ((cell  * LOADSIZE) < cSim.atoms) {
    uint2 cellStartEnd = {cell * LOADSIZE, min((cell + 1) * LOADSIZE, cSim.atoms)};

    // Iterate through cell
    unsigned int pos = cellStartEnd.x;

    // Read Atom Data and procedurally generate spline weights
    unsigned int maxatom = min(pos + LOADSIZE, cellStartEnd.y);
    unsigned int pos1 = pos + threadIdx.x;

    if (pos1 < maxatom) {
      unsigned realAtomindex = cSim.pImageAtom[pos1];

      if (cSim.pTIList[realAtomindex] > 0) {
        region[threadIdx.y][threadIdx.x] = 0;
      }
      else if (cSim.pTIList[realAtomindex + cSim.stride] > 0) {
        region[threadIdx.y][threadIdx.x] = 1;
      }
      else {
        region[threadIdx.y][threadIdx.x] = 2;
      }
      unsigned myRegion = region[threadIdx.y][threadIdx.x];

      PMEFloat charge = cSim.pAtomChargeSP[pos1];
      PMEFloat fx = cSim.pFractX[pos1];
      PMEFloat fy = cSim.pFractY[pos1];
      PMEFloat fz = cSim.pFractZ[pos1];

      int ix = int(fx);
      int iy = int(fy);
      int iz = int(fz);
      fx -= ix;
      fy -= iy;
      fz -= iz;
      ix -= cSim.orderMinusOne;
      iy -= cSim.orderMinusOne;
      iz -= cSim.orderMinusOne;
      if (ix < 0) {
        ix += cSim.nfft1;
      }
      if (iy < 0) {
        iy += cSim.nfft2;
      }
      if (iz < 0) {
        iz += cSim.nfft3;
      }
      sAtom[threadIdx.y][myRegion][threadIdx.x].ix = ix;
      sAtom[threadIdx.y][myRegion][threadIdx.x].iy = iy;
      sAtom[threadIdx.y][myRegion][threadIdx.x].iz = iz;

      // Order 2
      PMEFloat4 tx;
      tx.x = ONE - fx;
      tx.y = fx;

      // Order 3
      tx.z = ONEHALF * fx * tx.y;
      tx.y = ONEHALF * ((fx + ONE)*tx.x + (TWO - fx)*tx.y);
      tx.x = ONEHALF * (ONE - fx) * tx.x;

      // Order 4
      PMEFloat tt = charge*ONETHIRD;
      tx.w = tt * fx * tx.z;
      tx.z = tt * ((fx + ONE)*tx.y + (THREE - fx)*tx.z);
      tx.y = tt * ((fx + TWO)*tx.x + (TWO - fx)*tx.y);
      tx.x = tt * (ONE - fx) * tx.x;

      sAtom[threadIdx.y][myRegion][threadIdx.x].tx[0] = tx.x;
      sAtom[threadIdx.y][myRegion][threadIdx.x].tx[1] = tx.y;
      sAtom[threadIdx.y][myRegion][threadIdx.x].tx[2] = tx.z;
      sAtom[threadIdx.y][myRegion][threadIdx.x].tx[3] = tx.w;

      // Order 2
      PMEFloat4 ty;
      ty.x = ONE - fy;
      ty.y = fy;

      // Order 3
      ty.z = ONEHALF * fy * ty.y;
      ty.y = ONEHALF * ((fy + ONE)*ty.x + (TWO - fy)*ty.y);
      ty.x = ONEHALF * (ONE - fy) * ty.x;

      // Order 4
      ty.w = ONETHIRD * fy * ty.z;
      ty.z = ONETHIRD * ((fy + ONE)*ty.y + (THREE - fy)*ty.z);
      ty.y = ONETHIRD * ((fy + TWO)*ty.x + (TWO - fy)*ty.y);
      ty.x = ONETHIRD * (ONE - fy) * ty.x;

      sAtom[threadIdx.y][myRegion][threadIdx.x].ty[0] = ty.x;
      sAtom[threadIdx.y][myRegion][threadIdx.x].ty[1] = ty.y;
      sAtom[threadIdx.y][myRegion][threadIdx.x].ty[2] = ty.z;
      sAtom[threadIdx.y][myRegion][threadIdx.x].ty[3] = ty.w;

      // Order 2
      PMEFloat4 tz;
      tz.x = ONE - fz;
      tz.y = fz;

      // Order 3
      tz.z = ONEHALF * fz * tz.y;
      tz.y = ONEHALF * ((fz + ONE)*tz.x + (TWO - fz)*tz.y);
      tz.x = ONEHALF * (ONE - fz) * tz.x;

      // Order 4
      tz.w = ONETHIRD * fz * tz.z;
      tz.z = ONETHIRD * ((fz + ONE)*tz.y + (THREE - fz)*tz.z);
      tz.y = ONETHIRD * ((fz + TWO)*tz.x + (TWO - fz)*tz.y);
      tz.x = ONETHIRD * (ONE - fz) * tz.x;

      sAtom[threadIdx.y][myRegion][threadIdx.x].tz[0] = tz.x;
      sAtom[threadIdx.y][myRegion][threadIdx.x].tz[1] = tz.y;
      sAtom[threadIdx.y][myRegion][threadIdx.x].tz[2] = tz.z;
      sAtom[threadIdx.y][myRegion][threadIdx.x].tz[3] = tz.w;
    }
    __syncthreads();

    // Interpolate onto grid
    pos1 = 0;
    unsigned int lastAtom = min(LOADSIZE, cellStartEnd.y - pos);
    while (pos1 < lastAtom) {

      // Calculate values
      unsigned myRegion = region[threadIdx.y][pos1];
      int ix = sAtom[threadIdx.y][myRegion][pos1].ix + iOffsetX;
      int iy = sAtom[threadIdx.y][myRegion][pos1].iy + iOffsetY;
      int iz = sAtom[threadIdx.y][myRegion][pos1].iz + iOffsetZ;

      // Insure coordinates stay in bounds
      if (ix >= cSim.nfft1) {
        ix -= cSim.nfft1;
      }
      if (iy >= cSim.nfft2) {
        iy -= cSim.nfft2;
      }
      if (iz >= cSim.nfft3) {
        iz -= cSim.nfft3;
      }

      // Calculate interpolation values and destinations
      unsigned int gpos = ((ix & 0x3) + (iy & 3)*4 + (iz & 1)*16) + (ix >> 2)*32 +
                          (((iy >> 2) << 3) + ((iz >> 1) << 1)*cSim.nfft2)*cSim.nfft1;
      gpos += myRegion * cSim.XYZStride;
#ifdef use_DPFP
      unsigned long long int value = llitoulli(llrint(LATTICESCALE *
        sAtom[threadIdx.y][myRegion][pos1].tx[iOffsetX] * sAtom[threadIdx.y][myRegion][pos1].ty[iOffsetY] * sAtom[threadIdx.y][myRegion][pos1].tz[iOffsetZ]));
#else
      unsigned long long int value = fast_llrintf(LATTICESCALEF *
        sAtom[threadIdx.y][myRegion][pos1].tx[iOffsetX] * sAtom[threadIdx.y][myRegion][pos1].ty[iOffsetY] * sAtom[threadIdx.y][myRegion][pos1].tz[iOffsetZ]);
#endif
      atomicAdd((cSim.pPMEChargeGrid + gpos), value);
      pos1++;
    }
    __syncthreads();

    cell += gridDim.x * blockDim.y;
  }
}

//---------------------------------------------------------------------------------------------
// kgTIPMEReduceChargeGrid_kernel:
//---------------------------------------------------------------------------------------------
_kReduceFrcHead_ kgTIPMEReduceChargeGrid_kernel()
{
  unsigned long long int *pP0 = cSim.pPMEChargeGrid;
  unsigned long long int *pP1 = cSim.pPMEChargeGrid + cSim.XYZStride;
  unsigned long long int *pP2 = cSim.pPMEChargeGrid + (cSim.XYZStride * 2);
  PMEFloat *pF0               = cSim.pFFTChargeGrid;
  PMEFloat *pF1               = cSim.pFFTChargeGrid + cSim.XYZStride;;

  unsigned int iz     = blockIdx.x;
  unsigned int zShift = iz*cSim.nfft1xnfft2;
  unsigned int z0     = (iz & 1) * 16;
  unsigned int z1     = ((iz >> 1) << 1) * cSim.nfft1xnfft2;
  unsigned int zz     = z0 + z1;
  unsigned int iy     = threadIdx.y;

  while (iy < cSim.nfft2 ) {
    unsigned shift = iy*cSim.nfft1 + zShift;
    unsigned y0    = (iy & 3) * 4;
    unsigned y1    = ((iy >> 2) << 3) * cSim.nfft1;
    unsigned yy    = y0 + y1;

    unsigned int ix = threadIdx.x;
    while (ix < cSim.nfft1) {

      unsigned int spos = (ix & 0x3) + (ix >> 2)*32 + yy + zz;
      if (spos < cSim.nfft1xnfft2xnfft3) {
        long long int value0 = pP0[spos] + pP2[spos];
        long long int value1 = pP1[spos] + pP2[spos];
        pF0[ix + shift] = (PMEFloat)(value0) * ONEOVERLATTICESCALEF;
        pF1[ix + shift] = (PMEFloat)(value1) * ONEOVERLATTICESCALEF;
        pP0[spos] = ZeroF;
        pP1[spos] = ZeroF;
        pP2[spos] = ZeroF;
      }
      ix += blockDim.x;
    }
    iy += blockDim.y;
  }
}

#endif /* GTI */
