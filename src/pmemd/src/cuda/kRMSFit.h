#include "copyright.i"
{
#if !defined(CALC_COM) && !defined(KABSCH)
  int pos = threadIdx.x + blockDim.x * blockIdx.x;
#endif

#ifdef MASS_CALC
  int tid = threadIdx.x;
  __shared__ PMEDouble totMass;
#endif

#if defined(CALC_COM) || defined(KABSCH)
  int warpIdx = threadIdx.x >> 6; //there are 64*3=7192 total threads per block and 3 warpIdxs: 0, 1, and 2
 unsigned int tgx = threadIdx.x & 63; // this goes from 0-63 in each warp
  int pos = tgx + blockIdx.x * 64;
  int tid = threadIdx.x;
  __shared__ PMEDouble mass[64];
  __shared__ double sXData[THREADS_PER_BLOCK*3/16];
  __shared__ double sYData[THREADS_PER_BLOCK*3/16];
  __shared__ double sZData[THREADS_PER_BLOCK*3/16];
#endif

  while (pos < cSim.nShuttle) {

#if defined(CALC_COM) || defined(KABSCH)
  sXData[tid] = 0.0;
  sYData[tid] = 0.0;
  sZData[tid] = 0.0;

  if (warpIdx == 0) {
  #ifdef NEIGHBOR_LIST
    int shidx = cSim.pShuttleTickets[pos];
    int atidx = cSim.pImageAtomLookup[shidx];
    mass[tgx] = cSim.pImageMass[atidx];
  #else
    int atidx = cSim.pShuttleTickets[pos];
    mass[tgx] = cSim.pAtomMass[atidx];
  #endif
  }
  __syncthreads();
#endif

#ifdef MASS_CALC
  #ifdef NEIGHBOR_LIST
    int shidx = cSim.pShuttleTickets[pos];
    int atidx = cSim.pImageAtomLookup[shidx];
    PMEDouble mass = cSim.pImageMass[atidx];
  #else 
    int atidx = cSim.pShuttleTickets[pos];
    PMEDouble mass = cSim.pAtomMass[atidx];
  #endif
#endif

#ifdef CALC_COM
    if ((cSim.pAtmIdx[pos] == 1) || (cSim.pAtmIdx[pos] == 2)) {
      if (warpIdx == 0) {
        sXData[tid] = cSim.pDataShuttle[pos] * mass[tgx];
        sYData[tid] = cSim.pDataShuttle[pos + cSim.nShuttle] * mass[tgx];
        sZData[tid] = cSim.pDataShuttle[pos + 2 * cSim.nShuttle] * mass[tgx];
      }
      else if (warpIdx == 1) {
        sXData[tid] = cSim.pPrevDataShuttle[pos] * mass[tgx];
        sYData[tid] = cSim.pPrevDataShuttle[pos + cSim.nShuttle] * mass[tgx];
        sZData[tid] = cSim.pPrevDataShuttle[pos + 2 * cSim.nShuttle] * mass[tgx];
      }
      else if (warpIdx == 2) {
        sXData[tid] = cSim.pNextDataShuttle[pos] * mass[tgx];
        sYData[tid] = cSim.pNextDataShuttle[pos + cSim.nShuttle] * mass[tgx];
        sZData[tid] = cSim.pNextDataShuttle[pos + 2 * cSim.nShuttle] * mass[tgx];
      }
    }
    __syncthreads();
#elif defined(KABSCH)
  #ifdef NEXT_NEIGHBOR
    if ((cSim.pAtmIdx[pos] == 1) || (cSim.pAtmIdx[pos] == 2)) {
      sXData[tid] = cSim.pNextDataShuttle[pos] * mass[tgx] *
                    cSim.pDataShuttle[pos + warpIdx * cSim.nShuttle];
      sYData[tid] = cSim.pNextDataShuttle[pos + cSim.nShuttle] * mass[tgx] *
                    cSim.pDataShuttle[pos + warpIdx * cSim.nShuttle];
      sZData[tid] = cSim.pNextDataShuttle[pos + 2 * cSim.nShuttle] * mass[tgx] *
                    cSim.pDataShuttle[pos + warpIdx * cSim.nShuttle];
    }
    __syncthreads();
  #elif defined(PREV_NEIGHBOR)
    if ((cSim.pAtmIdx[pos] == 1) || (cSim.pAtmIdx[pos] == 2)) {
      sXData[tid] = cSim.pPrevDataShuttle[pos] * mass[tgx] *
                    cSim.pDataShuttle[pos + warpIdx * cSim.nShuttle];
      sYData[tid] = cSim.pPrevDataShuttle[pos + cSim.nShuttle] * mass[tgx] *
                    cSim.pDataShuttle[pos + warpIdx * cSim.nShuttle];
      sZData[tid] = cSim.pPrevDataShuttle[pos + 2 * cSim.nShuttle] * mass[tgx] *
                    cSim.pDataShuttle[pos + warpIdx * cSim.nShuttle];
    }
    __syncthreads();
  #endif
#endif

#ifdef MASS_CALC
    if (tid == 0) {
      totMass = 0.0;
    }
    __syncthreads();
    if ((cSim.pAtmIdx[pos] == 1) || (cSim.pAtmIdx[pos] == 2)) {
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600
      voidatomicAdd(&totMass, mass);
#else
      atomicAdd(&totMass, mass);
#endif
    }
    __syncthreads();
#endif

#if defined(CALC_COM) || defined(KABSCH) 
    for (unsigned int i = blockDim.x/6; i > 0; i >>= 1) {
      if (tgx < i && (pos + i) < cSim.nShuttle) {
        sXData[tid] += sXData[tid + i];
        sYData[tid] += sYData[tid + i];
        sZData[tid] += sZData[tid + i];
      }
      __syncthreads();
    }
#endif

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600
  #ifdef CALC_COM
    if (tgx == 0) {
      if (warpIdx == 0) {
        voidatomicAdd(&cSim.pSelfCOM[0], sXData[0] / (*cSim.ptotFitMass));
        voidatomicAdd(&cSim.pSelfCOM[1], sYData[0] / (*cSim.ptotFitMass));
        voidatomicAdd(&cSim.pSelfCOM[2], sZData[0] / (*cSim.ptotFitMass));
      }
      else if (warpIdx == 1) {
        voidatomicAdd(&cSim.pPrevCOM[0], sXData[64] / (*cSim.ptotFitMass));
        voidatomicAdd(&cSim.pPrevCOM[1], sYData[64] / (*cSim.ptotFitMass));
        voidatomicAdd(&cSim.pPrevCOM[2], sZData[64] / (*cSim.ptotFitMass));
      }
      else if (warpIdx == 2) {
        voidatomicAdd(&cSim.pNextCOM[0], sXData[128] / (*cSim.ptotFitMass));
        voidatomicAdd(&cSim.pNextCOM[1], sYData[128] / (*cSim.ptotFitMass));
        voidatomicAdd(&cSim.pNextCOM[2], sZData[128] / (*cSim.ptotFitMass));
      }
    }
  #elif defined(KABSCH)
    if (tgx == 0) {
    #ifdef NEXT_NEIGHBOR
        voidatomicAdd(&cSim.pNextKabsch[0 + warpIdx * 3], sXData[warpIdx * 64]);
        voidatomicAdd(&cSim.pNextKabsch[1 + warpIdx * 3], sYData[warpIdx * 64]);
        voidatomicAdd(&cSim.pNextKabsch[2 + warpIdx * 3], sZData[warpIdx * 64]);
    #elif defined(PREV_NEIGHBOR)
        voidatomicAdd(&cSim.pPrevKabsch[0 + warpIdx * 3], sXData[warpIdx * 64]);
        voidatomicAdd(&cSim.pPrevKabsch[1 + warpIdx * 3], sYData[warpIdx * 64]);
        voidatomicAdd(&cSim.pPrevKabsch[2 + warpIdx * 3], sZData[warpIdx * 64]);
    #endif
    }
  #elif defined(MASS_CALC)
    if (tid == 0) {
      voidatomicAdd(cSim.ptotFitMass, totMass);
    }
  #endif
#else
  #ifdef CALC_COM
    if (tgx == 0) {
      if (warpIdx == 0) {
        atomicAdd(&cSim.pSelfCOM[0], sXData[0] / (*cSim.ptotFitMass));
        atomicAdd(&cSim.pSelfCOM[1], sYData[0] / (*cSim.ptotFitMass));
        atomicAdd(&cSim.pSelfCOM[2], sZData[0] / (*cSim.ptotFitMass));
      }
      else if (warpIdx == 1) {
        atomicAdd(&cSim.pPrevCOM[0], sXData[64] / (*cSim.ptotFitMass));
        atomicAdd(&cSim.pPrevCOM[1], sYData[64] / (*cSim.ptotFitMass));
        atomicAdd(&cSim.pPrevCOM[2], sZData[64] / (*cSim.ptotFitMass));
      }
      else if (warpIdx == 2) {
        atomicAdd(&cSim.pNextCOM[0], sXData[128] / (*cSim.ptotFitMass));
        atomicAdd(&cSim.pNextCOM[1], sYData[128] / (*cSim.ptotFitMass));
        atomicAdd(&cSim.pNextCOM[2], sZData[128] / (*cSim.ptotFitMass));
      }
    }
  #elif defined(KABSCH)
    if (tgx == 0) {
    #ifdef NEXT_NEIGHBOR
        atomicAdd(&cSim.pNextKabsch[0 + warpIdx * 3], sXData[warpIdx * 64]);
        atomicAdd(&cSim.pNextKabsch[1 + warpIdx * 3], sYData[warpIdx * 64]);
        atomicAdd(&cSim.pNextKabsch[2 + warpIdx * 3], sZData[warpIdx * 64]);
    #elif defined(PREV_NEIGHBOR)
        atomicAdd(&cSim.pPrevKabsch[0 + warpIdx * 3], sXData[warpIdx * 64]);
        atomicAdd(&cSim.pPrevKabsch[1 + warpIdx * 3], sYData[warpIdx * 64]);
        atomicAdd(&cSim.pPrevKabsch[2 + warpIdx * 3], sZData[warpIdx * 64]);
    #endif
    }
  #elif defined(MASS_CALC)
    if (tid == 0) {
      atomicAdd(cSim.ptotFitMass, totMass);
    }
  #endif
#endif

#ifdef CENTER_DATA
      cSim.pDataShuttle[                    pos] = cSim.pDataShuttle[                    pos] - cSim.pSelfCOM[0];
      cSim.pDataShuttle[    cSim.nShuttle + pos] = cSim.pDataShuttle[    cSim.nShuttle + pos] - cSim.pSelfCOM[1];
      cSim.pDataShuttle[2 * cSim.nShuttle + pos] = cSim.pDataShuttle[2 * cSim.nShuttle + pos] - cSim.pSelfCOM[2];
      cSim.pPrevDataShuttle[                    pos] = cSim.pPrevDataShuttle[                    pos] - cSim.pPrevCOM[0];
      cSim.pPrevDataShuttle[    cSim.nShuttle + pos] = cSim.pPrevDataShuttle[    cSim.nShuttle + pos] - cSim.pPrevCOM[1];
      cSim.pPrevDataShuttle[2 * cSim.nShuttle + pos] = cSim.pPrevDataShuttle[2 * cSim.nShuttle + pos] - cSim.pPrevCOM[2];
      cSim.pNextDataShuttle[                    pos] = cSim.pNextDataShuttle[                    pos] - cSim.pNextCOM[0];
      cSim.pNextDataShuttle[    cSim.nShuttle + pos] = cSim.pNextDataShuttle[    cSim.nShuttle + pos] - cSim.pNextCOM[1];
      cSim.pNextDataShuttle[2 * cSim.nShuttle + pos] = cSim.pNextDataShuttle[2 * cSim.nShuttle + pos] - cSim.pNextCOM[2];
#endif

#ifdef ROTATE_NEIGHBOR
  #ifdef NEXT_NEIGHBOR
    double xtmp = cSim.pNextDataShuttle[                    pos] * cSim.pRotAtm[0]
                + cSim.pNextDataShuttle[    cSim.nShuttle + pos] * cSim.pRotAtm[1]
                + cSim.pNextDataShuttle[2 * cSim.nShuttle + pos] * cSim.pRotAtm[2];
    double ytmp = cSim.pNextDataShuttle[                  + pos] * cSim.pRotAtm[3]
                + cSim.pNextDataShuttle[    cSim.nShuttle + pos] * cSim.pRotAtm[4]
                + cSim.pNextDataShuttle[2 * cSim.nShuttle + pos] * cSim.pRotAtm[5];
    double ztmp = cSim.pNextDataShuttle[                  + pos] * cSim.pRotAtm[6]
                + cSim.pNextDataShuttle[    cSim.nShuttle + pos] * cSim.pRotAtm[7]
                + cSim.pNextDataShuttle[2 * cSim.nShuttle + pos] * cSim.pRotAtm[8];
    cSim.pNextDataShuttle[                    pos] = xtmp;
    cSim.pNextDataShuttle[    cSim.nShuttle + pos] = ytmp;
    cSim.pNextDataShuttle[2 * cSim.nShuttle + pos] = ztmp;
  #elif defined(PREV_NEIGHBOR)
    double xtmp = cSim.pPrevDataShuttle[                    pos] * cSim.pRotAtm[0]
                + cSim.pPrevDataShuttle[    cSim.nShuttle + pos] * cSim.pRotAtm[1]
                + cSim.pPrevDataShuttle[2 * cSim.nShuttle + pos] * cSim.pRotAtm[2];
    double ytmp = cSim.pPrevDataShuttle[                  + pos] * cSim.pRotAtm[3]
                + cSim.pPrevDataShuttle[    cSim.nShuttle + pos] * cSim.pRotAtm[4]
                + cSim.pPrevDataShuttle[2 * cSim.nShuttle + pos] * cSim.pRotAtm[5];
    double ztmp = cSim.pPrevDataShuttle[                  + pos] * cSim.pRotAtm[6]
                + cSim.pPrevDataShuttle[    cSim.nShuttle + pos] * cSim.pRotAtm[7]
                + cSim.pPrevDataShuttle[2 * cSim.nShuttle + pos] * cSim.pRotAtm[8];
    cSim.pPrevDataShuttle[                    pos] = xtmp;
    cSim.pPrevDataShuttle[    cSim.nShuttle + pos] = ytmp;
    cSim.pPrevDataShuttle[2 * cSim.nShuttle + pos] = ztmp;
  #endif
#endif

#ifdef RECENTER_DATA
    cSim.pDataShuttle[                    pos] = cSim.pDataShuttle[                    pos] + cSim.pSelfCOM[0];
    cSim.pDataShuttle[    cSim.nShuttle + pos] = cSim.pDataShuttle[    cSim.nShuttle + pos] + cSim.pSelfCOM[1];
    cSim.pDataShuttle[2 * cSim.nShuttle + pos] = cSim.pDataShuttle[2 * cSim.nShuttle + pos] + cSim.pSelfCOM[2];
    cSim.pNextDataShuttle[                    pos] = cSim.pNextDataShuttle[                    pos] + cSim.pSelfCOM[0];
    cSim.pNextDataShuttle[    cSim.nShuttle + pos] = cSim.pNextDataShuttle[    cSim.nShuttle + pos] + cSim.pSelfCOM[1];
    cSim.pNextDataShuttle[2 * cSim.nShuttle + pos] = cSim.pNextDataShuttle[2 * cSim.nShuttle + pos] + cSim.pSelfCOM[2];
    cSim.pPrevDataShuttle[                    pos] = cSim.pPrevDataShuttle[                    pos] + cSim.pSelfCOM[0];
    cSim.pPrevDataShuttle[    cSim.nShuttle + pos] = cSim.pPrevDataShuttle[    cSim.nShuttle + pos] + cSim.pSelfCOM[1];
    cSim.pPrevDataShuttle[2 * cSim.nShuttle + pos] = cSim.pPrevDataShuttle[2 * cSim.nShuttle + pos] + cSim.pSelfCOM[2];
#endif


#if defined(CALC_COM) || defined(KABSCH)
    pos += gridDim.x*64;
#else
    pos += gridDim.x*256;
#endif
  }
}
