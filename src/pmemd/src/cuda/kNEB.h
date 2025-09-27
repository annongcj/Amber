#include "copyright.i"
{
  int warpIdx = threadIdx.x >> GRID_BITS;
 unsigned int tgx = threadIdx.x & GRID_BITS_MASK;
  int pos = tgx + (GRID * blockIdx.x);  // Position in the 1D array.
#ifdef REVISED_TAN
  #if defined(THIRD) || defined(FOURTH)
    double emax = 0.0;
    double emin = 0.0;
    int replica = cSim.beadid - 1;
  #endif
#endif

#ifdef NEB_FRC
  __shared__ int atidx[GRID];  
#endif

#ifdef NORMALIZE
  __shared__ double norm;
#endif

  while (pos < cSim.nShuttle) {
    int idx = pos + warpIdx * cSim.nShuttle;  // Position in the extended 3D array.

#ifdef NEB_FRC
    if (warpIdx == 0) {
#ifdef NEIGHBOR_LIST
      int shidx = cSim.pShuttleTickets[pos];
      atidx[tgx] = cSim.pImageAtomLookup[shidx];
#else
      atidx[tgx] = cSim.pShuttleTickets[pos];
#endif
    }
    __syncthreads();
#endif

    if ((cSim.pAtmIdx[pos] == 0) || (cSim.pAtmIdx[pos] == 2)) {
#ifdef REVISED_TAN
  #ifdef FIRST
      cSim.pTangents[idx] = cSim.pNextDataShuttle[idx] - cSim.pDataShuttle[idx];
  #elif defined(SECOND)
      cSim.pTangents[idx] = cSim.pDataShuttle[idx] - cSim.pPrevDataShuttle[idx];
  #elif defined(THIRD)
      emax = max(abs(cSim.pNEBEnergyAll[replica + 1] - cSim.pNEBEnergyAll[replica]),
                 abs(cSim.pNEBEnergyAll[replica - 1] - cSim.pNEBEnergyAll[replica]));
      emin = min(abs(cSim.pNEBEnergyAll[replica + 1] - cSim.pNEBEnergyAll[replica]),
                 abs(cSim.pNEBEnergyAll[replica - 1] - cSim.pNEBEnergyAll[replica]));
      cSim.pTangents[idx] = (cSim.pNextDataShuttle[idx] - cSim.pDataShuttle[idx]) * emax;
      cSim.pTangents[idx] += (cSim.pDataShuttle[idx] - cSim.pPrevDataShuttle[idx]) * emin;
  #elif defined(FOURTH)
      emax = max(abs(cSim.pNEBEnergyAll[replica + 1] - cSim.pNEBEnergyAll[replica]),
                 abs(cSim.pNEBEnergyAll[replica - 1] - cSim.pNEBEnergyAll[replica]));
      emin = min(abs(cSim.pNEBEnergyAll[replica + 1] - cSim.pNEBEnergyAll[replica]),
                 abs(cSim.pNEBEnergyAll[replica - 1] - cSim.pNEBEnergyAll[replica]));
      cSim.pTangents[idx] = (cSim.pNextDataShuttle[idx] - cSim.pDataShuttle[idx]) * emin;
      cSim.pTangents[idx] += (cSim.pDataShuttle[idx] - cSim.pPrevDataShuttle[idx]) * emax;
  #endif
#elif defined(BASIC_TAN)
      cSim.pTangents[idx] = cSim.pNextDataShuttle[idx] - cSim.pDataShuttle[idx];
#endif
    }
    else {
      cSim.pTangents[idx] = 0.0;
    }

#ifdef NORMALIZE
  if (threadIdx.x == 0) {
    if (*cSim.pNorm > 1.0e-06) {
      norm = 1.0/sqrt(*cSim.pNorm);
    }
    else {
      norm = 0.0;
    }
  }
  __syncthreads();

  if ((cSim.pAtmIdx[pos] == 0) || (cSim.pAtmIdx[pos] == 2)) {
    cSim.pTangents[idx] *= norm;
  }

  cSim.pSpringForce[idx] = ((cSim.pNextDataShuttle[idx] - cSim.pDataShuttle[idx]) * (*cSim.pSpring2)
  - (cSim.pDataShuttle[idx] - cSim.pPrevDataShuttle[idx]) * (*cSim.pSpring1)) * cSim.pTangents[idx];
#endif

#ifdef NEB_FRC
  //cSim.pNEBForce[idx] = ((*cSim.pDotProduct2) * cSim.pTangents[idx])
  //                    - ((*cSim.pDotProduct1) * cSim.pTangents[idx]);//comment here after test
  PMEAccumulator neb_frc = (((*cSim.pDotProduct2) * cSim.pTangents[idx])
                         - ((*cSim.pDotProduct1) * cSim.pTangents[idx])) * FORCESCALE;
  cSim.pForceAccumulator[atidx[tgx] + warpIdx * cSim.stride] += neb_frc;
#endif

    pos += gridDim.x * GRID;
  }
}
