#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// This code is included in kForcesUpdate.cu
//
// #defines: NODPTEXTURE, NTP_LOTSOFMOLECULES
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// KPMECALCULATECOM_KERNEL: the names of this and the other five kernels are simply stand-in
//                          macros, as the caps might suggest.  This kernel calculates the
//                          center of mass for each molecule in the system.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(THREADS_PER_BLOCK, 1)
KPMECALCULATECOM_KERNEL()
{
  struct COM1 {
    double COMX;
    double COMY;
    double COMZ;
  };

  struct UllCOM1 {
    PMEUllInt COMX;
    PMEUllInt COMY;
    PMEUllInt COMZ;
  };
  __shared__ volatile COM1 sC[THREADS_PER_BLOCK];
  __shared__ volatile UllCOM1 sCOM[MAXMOLECULES];
#ifdef NTP_LOTSOFMOLECULES
  const int maxMolecules = MAXMOLECULES;
#endif
  unsigned int pos;

  // Clear shared memory COM
  pos = threadIdx.x;
#ifdef NTP_LOTSOFMOLECULES
  while (pos < maxMolecules) {
#else
  while (pos < cSim.soluteMolecules) {
#endif
    sCOM[pos].COMX = (PMEUllInt)0;
    sCOM[pos].COMY = (PMEUllInt)0;
    sCOM[pos].COMZ = (PMEUllInt)0;
    pos += blockDim.x;
#ifdef NTP_LOTSOFMOLECULES
  }
#else
  }
#endif
  __syncthreads();

  // Calculate solute COM
  pos = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int tgx   = threadIdx.x & (GRID - 1);
  volatile COM1* psC = &sC[threadIdx.x - tgx];
  int oldMoleculeID  = -2;
  int atomID, moleculeID;
  PMEDouble CX = (PMEDouble)0.0;
  PMEDouble CY = (PMEDouble)0.0;
  PMEDouble CZ = (PMEDouble)0.0;
  PMEMask mask1 = __BALLOT(WARP_MASK, pos < cSim.soluteAtoms);
  while (pos < cSim.soluteAtoms) {
    atomID     = cSim.pImageSoluteAtomID[pos];
    moleculeID = cSim.pSoluteAtomMoleculeID[pos];
    PMEDouble mass = (PMEDouble)0.0;
    PMEDouble x, y, z;
    if (atomID != -1) {
      mass = cSim.pSoluteAtomMass[pos];
      x    = cSim.pImageX[atomID];
      y    = cSim.pImageY[atomID];
      z    = cSim.pImageZ[atomID];
    }

    // Output COM upon changed status
    PMEMask mask2 = __BALLOT(mask1, moleculeID != oldMoleculeID);
    if (moleculeID != oldMoleculeID) {
      psC[tgx].COMX = CX;
      psC[tgx].COMY = CY;
      psC[tgx].COMZ = CZ;
      __SYNCWARP(mask2);
      PMEMask mask3 = __BALLOT(mask2, oldMoleculeID >= 0);
      if (oldMoleculeID >= 0) {
#ifdef AMBER_PLATFORM_AMD_WARP64
        if (tgx < 32) {
          psC[tgx].COMX += psC[tgx + 32].COMX;
          psC[tgx].COMY += psC[tgx + 32].COMY;
          psC[tgx].COMZ += psC[tgx + 32].COMZ;
        }
        __SYNCWARP(mask3);
#endif
        if (tgx < 16) {
          psC[tgx].COMX += psC[tgx + 16].COMX;
          psC[tgx].COMY += psC[tgx + 16].COMY;
          psC[tgx].COMZ += psC[tgx + 16].COMZ;
        }
        __SYNCWARP(mask3);
        if (tgx < 8) {
          psC[tgx].COMX += psC[tgx + 8].COMX;
          psC[tgx].COMY += psC[tgx + 8].COMY;
          psC[tgx].COMZ += psC[tgx + 8].COMZ;
        }
        __SYNCWARP(mask3);
        if (tgx < 4) {
          psC[tgx].COMX += psC[tgx + 4].COMX;
          psC[tgx].COMY += psC[tgx + 4].COMY;
          psC[tgx].COMZ += psC[tgx + 4].COMZ;
        }
        __SYNCWARP(mask3);
        if (tgx < 2) {
          psC[tgx].COMX += psC[tgx + 2].COMX;
          psC[tgx].COMY += psC[tgx + 2].COMY;
          psC[tgx].COMZ += psC[tgx + 2].COMZ;
        }
        __SYNCWARP(mask3);
        if (tgx == 0) {
          psC->COMX += psC[1].COMX;
          psC->COMY += psC[1].COMY;
          psC->COMZ += psC[1].COMZ;
          unsigned long long int val1  = llitoulli(llrint(psC->COMX * ENERGYSCALE));
          unsigned long long int val2  = llitoulli(llrint(psC->COMY * ENERGYSCALE));
          unsigned long long int val3  = llitoulli(llrint(psC->COMZ * ENERGYSCALE));

#ifdef NTP_LOTSOFMOLECULES
          if (oldMoleculeID < maxMolecules) {
            atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
            atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
            atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
          }
          else {
            atomicAdd(&cSim.pSoluteUllCOMX[oldMoleculeID], val1);
            atomicAdd(&cSim.pSoluteUllCOMY[oldMoleculeID], val2);
            atomicAdd(&cSim.pSoluteUllCOMZ[oldMoleculeID], val3);
          }
#else  // NTP_LOTSOFMOLECULES
          atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
          atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
          atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
#endif // NTP_LOTSOFMOLECULES
        }
      }
      CX = (PMEDouble)0.0;
      CY = (PMEDouble)0.0;
      CZ = (PMEDouble)0.0;
    }
    oldMoleculeID = moleculeID;
    if (mass != 0.0) {
      CX += mass * x;
      CY += mass * y;
      CZ += mass * z;
    }
    pos += blockDim.x * gridDim.x;
    mask1 = __BALLOT(mask1, pos < cSim.soluteAtoms);
  }

  // Dump last batch of solute data to shared memory
  psC[tgx].COMX = CX;
  psC[tgx].COMY = CY;
  psC[tgx].COMZ = CZ;
  mask1 = __BALLOT(WARP_MASK, oldMoleculeID >= 0);
  if (oldMoleculeID >= 0) {
#ifdef AMBER_PLATFORM_AMD_WARP64
    if (tgx < 32) {
      psC[tgx].COMX += psC[tgx + 32].COMX;
      psC[tgx].COMY += psC[tgx + 32].COMY;
      psC[tgx].COMZ += psC[tgx + 32].COMZ;
    }
    __SYNCWARP(mask1);
#endif
    if (tgx < 16) {
      psC[tgx].COMX += psC[tgx + 16].COMX;
      psC[tgx].COMY += psC[tgx + 16].COMY;
      psC[tgx].COMZ += psC[tgx + 16].COMZ;
    }
    __SYNCWARP(mask1);
    if (tgx < 8) {
      psC[tgx].COMX += psC[tgx + 8].COMX;
      psC[tgx].COMY += psC[tgx + 8].COMY;
      psC[tgx].COMZ += psC[tgx + 8].COMZ;
    }
    __SYNCWARP(mask1);
    if (tgx < 4) {
      psC[tgx].COMX += psC[tgx + 4].COMX;
      psC[tgx].COMY += psC[tgx + 4].COMY;
      psC[tgx].COMZ += psC[tgx + 4].COMZ;
    }
    __SYNCWARP(mask1);
    if (tgx < 2) {
      psC[tgx].COMX += psC[tgx + 2].COMX;
      psC[tgx].COMY += psC[tgx + 2].COMY;
      psC[tgx].COMZ += psC[tgx + 2].COMZ;
    }
    __SYNCWARP(mask1);
    if (tgx == 0) {
      psC->COMX += psC[1].COMX;
      psC->COMY += psC[1].COMY;
      psC->COMZ += psC[1].COMZ;
      unsigned long long int val1 = llitoulli(llrint(psC->COMX * ENERGYSCALE));
      unsigned long long int val2 = llitoulli(llrint(psC->COMY * ENERGYSCALE));
      unsigned long long int val3 = llitoulli(llrint(psC->COMZ * ENERGYSCALE));
#ifdef NTP_LOTSOFMOLECULES
      if (oldMoleculeID < maxMolecules) {
        atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
        atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
        atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
      }
      else {
        atomicAdd(&cSim.pSoluteUllCOMX[oldMoleculeID], val1);
        atomicAdd(&cSim.pSoluteUllCOMY[oldMoleculeID], val2);
        atomicAdd(&cSim.pSoluteUllCOMZ[oldMoleculeID], val3);
      }
#else  // NTP_LOTSOFMOLECULES
      atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
      atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
      atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
#endif // NTP_LOTSOFMOLECULES
    }
  }
  // Dump solute data to main memory
  __syncthreads();
  pos = threadIdx.x;
#ifdef NTP_LOTSOFMOLECULES
  while (pos < maxMolecules) {
#else
  while (pos < cSim.soluteMolecules) {
#endif
    if (sCOM[pos].COMX != 0) {
      atomicAdd(&cSim.pSoluteUllCOMX[pos], sCOM[pos].COMX);
    }
    if (sCOM[pos].COMY != 0) {
      atomicAdd(&cSim.pSoluteUllCOMY[pos], sCOM[pos].COMY);
    }
    if (sCOM[pos].COMZ != 0) {
      atomicAdd(&cSim.pSoluteUllCOMZ[pos], sCOM[pos].COMZ);
    }
    pos += blockDim.x;
#ifdef NTP_LOTSOFMOLECULES
  }
#else
  }
#endif

  // Solvent atoms
  pos = blockIdx.x*blockDim.x + threadIdx.x;
  while (pos < cSim.solventMolecules) {
    int4 atomID       = cSim.pImageSolventAtomID[pos];
    PMEDouble invMass = cSim.pSolventInvMass[pos];
    PMEDouble mass = cSim.pSolventAtomMass1[pos];
    PMEDouble x    = cSim.pImageX[atomID.x];
    PMEDouble y    = cSim.pImageY[atomID.x];
    PMEDouble z    = cSim.pImageZ[atomID.x];
    PMEDouble CX   = mass * x;
    PMEDouble CY   = mass * y;
    PMEDouble CZ   = mass * z;
    if (atomID.y != -1) {
      PMEDouble mass = cSim.pSolventAtomMass2[pos];
      PMEDouble x    = cSim.pImageX[atomID.y];
      PMEDouble y    = cSim.pImageY[atomID.y];
      PMEDouble z    = cSim.pImageZ[atomID.y];
      CX            += mass * x;
      CY            += mass * y;
      CZ            += mass * z;
    }
    if (atomID.z != -1) {
      PMEDouble mass = cSim.pSolventAtomMass3[pos];
      PMEDouble x    = cSim.pImageX[atomID.z];
      PMEDouble y    = cSim.pImageY[atomID.z];
      PMEDouble z    = cSim.pImageZ[atomID.z];
      CX            += mass * x;
      CY            += mass * y;
      CZ            += mass * z;
    }
    if (atomID.w != -1) {
      PMEDouble mass = cSim.pSolventAtomMass4[pos];
      PMEDouble x    = cSim.pImageX[atomID.w];
      PMEDouble y    = cSim.pImageY[atomID.w];
      PMEDouble z    = cSim.pImageZ[atomID.w];
      CX            += mass * x;
      CY            += mass * y;
      CZ            += mass * z;
    }

    // Sum up results
    cSim.pSolventCOMX[pos]  = invMass * CX;
    cSim.pSolventCOMY[pos]  = invMass * CY;
    cSim.pSolventCOMZ[pos]  = invMass * CZ;
    pos                    += blockDim.x * gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// KPMECALCULATESOLUTECOM_KERNEL: calculate the solute center of mass.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(THREADS_PER_BLOCK, 1)
KPMECALCULATESOLUTECOM_KERNEL()
{
  struct COM1 {
    double COMX;
    double COMY;
    double COMZ;
  };

  struct UllCOM1 {
    PMEUllInt COMX;
    PMEUllInt COMY;
    PMEUllInt COMZ;
  };

  __shared__ volatile COM1 sC[THREADS_PER_BLOCK];
  __shared__ volatile UllCOM1 sCOM[MAXMOLECULES];
#ifdef NTP_LOTSOFMOLECULES
  const int maxMolecules = MAXMOLECULES;
#endif
  unsigned int pos;

  // Clear shared memory COM
  pos = threadIdx.x;
#ifdef NTP_LOTSOFMOLECULES
  while (pos < maxMolecules) {
#else
  while (pos < cSim.soluteMolecules) {
#endif
        sCOM[pos].COMX                              = (PMEUllInt)0;
        sCOM[pos].COMY                              = (PMEUllInt)0;
        sCOM[pos].COMZ                              = (PMEUllInt)0;
        pos                                        += blockDim.x;
#ifdef NTP_LOTSOFMOLECULES
  }
#else
  }
#endif
  __syncthreads();

  // Calculate solute COM
  pos = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int tgx   = threadIdx.x & (GRID - 1);
  volatile COM1* psC = &sC[threadIdx.x - tgx];
  int oldMoleculeID  = -2;
  int atomID, moleculeID;
  PMEDouble CX = (PMEDouble)0.0;
  PMEDouble CY = (PMEDouble)0.0;
  PMEDouble CZ = (PMEDouble)0.0;
  while (pos < cSim.soluteAtoms) {
    atomID     = cSim.pImageSoluteAtomID[pos];
    moleculeID = cSim.pSoluteAtomMoleculeID[pos];
    PMEDouble mass = (PMEDouble)0.0;
    PMEDouble x, y, z;
    if (atomID != -1) {
      mass = cSim.pSoluteAtomMass[pos];
      x    = cSim.pImageX[atomID];
      y    = cSim.pImageY[atomID];
      z    = cSim.pImageZ[atomID];
    }

    // Output COM upon changed status
    if (moleculeID != oldMoleculeID) {
      psC[tgx].COMX = CX;
      psC[tgx].COMY = CY;
      psC[tgx].COMZ = CZ;
      if (oldMoleculeID >= 0) {
#ifdef AMBER_PLATFORM_AMD_WARP64
        if (tgx < 32) {
          psC[tgx].COMX += psC[tgx + 32].COMX;
          psC[tgx].COMY += psC[tgx + 32].COMY;
          psC[tgx].COMZ += psC[tgx + 32].COMZ;
        }
#endif
        if (tgx < 16) {
          psC[tgx].COMX += psC[tgx + 16].COMX;
          psC[tgx].COMY += psC[tgx + 16].COMY;
          psC[tgx].COMZ += psC[tgx + 16].COMZ;
        }
        if (tgx < 8) {
          psC[tgx].COMX += psC[tgx + 8].COMX;
          psC[tgx].COMY += psC[tgx + 8].COMY;
          psC[tgx].COMZ += psC[tgx + 8].COMZ;
        }
        if (tgx < 4) {
          psC[tgx].COMX += psC[tgx + 4].COMX;
          psC[tgx].COMY += psC[tgx + 4].COMY;
          psC[tgx].COMZ += psC[tgx + 4].COMZ;
        }
        if (tgx < 2) {
          psC[tgx].COMX += psC[tgx + 2].COMX;
          psC[tgx].COMY += psC[tgx + 2].COMY;
          psC[tgx].COMZ += psC[tgx + 2].COMZ;
        }
        if (tgx == 0) {
          psC->COMX += psC[1].COMX;
          psC->COMY += psC[1].COMY;
          psC->COMZ += psC[1].COMZ;
          unsigned long long int val1 = llitoulli(llrint(psC->COMX * ENERGYSCALE));
          unsigned long long int val2 = llitoulli(llrint(psC->COMY * ENERGYSCALE));
          unsigned long long int val3 = llitoulli(llrint(psC->COMZ * ENERGYSCALE));
#ifdef NTP_LOTSOFMOLECULES
          if (oldMoleculeID < maxMolecules) {
            atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
            atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
            atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
          }
          else {
            atomicAdd(&cSim.pSoluteUllCOMX[oldMoleculeID], val1);
            atomicAdd(&cSim.pSoluteUllCOMY[oldMoleculeID], val2);
            atomicAdd(&cSim.pSoluteUllCOMZ[oldMoleculeID], val3);
          }
#else  // NTP_LOTSOFMOLECULES
          atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
          atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
          atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
#endif // NTP_LOTSOFMOLECULES
        }
      }
      CX = (PMEDouble)0.0;
      CY = (PMEDouble)0.0;
      CZ = (PMEDouble)0.0;
    }
    oldMoleculeID = moleculeID;
    if (atomID != -1) {
      CX += mass * x;
      CY += mass * y;
      CZ += mass * z;
    }
    pos += blockDim.x * gridDim.x;
  }

  // Dump last batch of solute data to shared memory
  psC[tgx].COMX = CX;
  psC[tgx].COMY = CY;
  psC[tgx].COMZ = CZ;
  if (oldMoleculeID >= 0) {
#ifdef AMBER_PLATFORM_AMD_WARP64
    if (tgx < 32) {
      psC[tgx].COMX += psC[tgx + 32].COMX;
      psC[tgx].COMY += psC[tgx + 32].COMY;
      psC[tgx].COMZ += psC[tgx + 32].COMZ;
    }
#endif
    if (tgx < 16) {
      psC[tgx].COMX += psC[tgx + 16].COMX;
      psC[tgx].COMY += psC[tgx + 16].COMY;
      psC[tgx].COMZ += psC[tgx + 16].COMZ;
    }
    if (tgx < 8) {
      psC[tgx].COMX += psC[tgx + 8].COMX;
      psC[tgx].COMY += psC[tgx + 8].COMY;
      psC[tgx].COMZ += psC[tgx + 8].COMZ;
    }
    if (tgx < 4) {
      psC[tgx].COMX += psC[tgx + 4].COMX;
      psC[tgx].COMY += psC[tgx + 4].COMY;
      psC[tgx].COMZ += psC[tgx + 4].COMZ;
    }
    if (tgx < 2) {
      psC[tgx].COMX += psC[tgx + 2].COMX;
      psC[tgx].COMY += psC[tgx + 2].COMY;
      psC[tgx].COMZ += psC[tgx + 2].COMZ;
    }
    if (tgx == 0) {
      psC->COMX += psC[1].COMX;
      psC->COMY += psC[1].COMY;
      psC->COMZ += psC[1].COMZ;
      unsigned long long int val1 = llitoulli(llrint(psC->COMX * ENERGYSCALE));
      unsigned long long int val2 = llitoulli(llrint(psC->COMY * ENERGYSCALE));
      unsigned long long int val3 = llitoulli(llrint(psC->COMZ * ENERGYSCALE));

#ifdef NTP_LOTSOFMOLECULES
      if (oldMoleculeID < maxMolecules) {
        atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
        atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
        atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
      }
      else {
        atomicAdd(&cSim.pSoluteUllCOMX[oldMoleculeID], val1);
        atomicAdd(&cSim.pSoluteUllCOMY[oldMoleculeID], val2);
        atomicAdd(&cSim.pSoluteUllCOMZ[oldMoleculeID], val3);
      }
#else  // NTP_LOTSOFMOLECULES
      atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
      atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
      atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
#endif // NTP_LOTSOFMOLECULES
    }
  }

  // Dump solute data to main memory
  __syncthreads();
  pos = threadIdx.x;
#ifdef NTP_LOTSOFMOLECULES
  while (pos < maxMolecules) {
#else
  while (pos < cSim.soluteMolecules) {
#endif
    if (sCOM[pos].COMX != 0) {
      atomicAdd(&cSim.pSoluteUllCOMX[pos], sCOM[pos].COMX);
    }
    if (sCOM[pos].COMY != 0) {
      atomicAdd(&cSim.pSoluteUllCOMY[pos], sCOM[pos].COMY);
    }
    if (sCOM[pos].COMZ != 0) {
      atomicAdd(&cSim.pSoluteUllCOMZ[pos], sCOM[pos].COMZ);
    }
    pos += blockDim.x;
#ifdef NTP_LOTSOFMOLECULES
  }
#else
  }
#endif
}

//---------------------------------------------------------------------------------------------
// KPMECALCULATECOMKINETICENERGY_KERNEL: calculate the kinetic energy of the center of mass
//                                       for each molecule in the system.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(THREADS_PER_BLOCK, 1)
KPMECALCULATECOMKINETICENERGY_KERNEL()
{
  struct COMKineticEnergy {
    double EKCOMX;
    double EKCOMY;
    double EKCOMZ;
  };

  struct UllCOMKineticEnergy {
    PMEUllInt EKCOMX;
    PMEUllInt EKCOMY;
    PMEUllInt EKCOMZ;
  };

  __shared__ volatile COMKineticEnergy sE[THREADS_PER_BLOCK];
  __shared__ volatile UllCOMKineticEnergy sEKCOM[MAXMOLECULES];
#ifdef NTP_LOTSOFMOLECULES
  const int maxMolecules = MAXMOLECULES;
#endif
  unsigned int pos;

  // Clear EKCOM
  pos = threadIdx.x;
#ifdef NTP_LOTSOFMOLECULES
  while (pos < maxMolecules) {
#else
  while (pos < cSim.soluteMolecules) {
#endif
    sEKCOM[pos].EKCOMX  = (PMEUllInt)0;
    sEKCOM[pos].EKCOMY  = (PMEUllInt)0;
    sEKCOM[pos].EKCOMZ  = (PMEUllInt)0;
    pos                += blockDim.x;
#ifdef NTP_LOTSOFMOLECULES
  }
#else
  }
#endif
  __syncthreads();

  // Calculate solute kinetic energy
  pos = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int tgx               = threadIdx.x & (GRID - 1);
  volatile COMKineticEnergy* psE = &sE[threadIdx.x - tgx];
  int oldMoleculeID = -2;
  PMEDouble EX = (PMEDouble)0.0;
  PMEDouble EY = (PMEDouble)0.0;
  PMEDouble EZ = (PMEDouble)0.0;
  while (pos < cSim.soluteAtoms) {
    int atomID     = cSim.pImageSoluteAtomID[pos];
    int moleculeID = cSim.pSoluteAtomMoleculeID[pos];
    PMEDouble mass = (PMEDouble)0.0;
    PMEDouble vx, vy, vz;
    if (atomID != -1) {
      mass = cSim.pSoluteAtomMass[pos];
      vx   = cSim.pImageVelX[atomID];
      vy   = cSim.pImageVelY[atomID];
      vz   = cSim.pImageVelZ[atomID];
    }

    // Output EKCOM upon changed status
    if ((moleculeID != oldMoleculeID) && (oldMoleculeID != -2)) {
      psE[tgx].EKCOMX = EX;
      psE[tgx].EKCOMY = EY;
      psE[tgx].EKCOMZ = EZ;
#ifdef AMBER_PLATFORM_AMD_WARP64
      if (tgx < 32) {
        psE[tgx].EKCOMX += psE[tgx + 32].EKCOMX;
        psE[tgx].EKCOMY += psE[tgx + 32].EKCOMY;
        psE[tgx].EKCOMZ += psE[tgx + 32].EKCOMZ;
      }
#endif
      if (tgx < 16) {
        psE[tgx].EKCOMX += psE[tgx + 16].EKCOMX;
        psE[tgx].EKCOMY += psE[tgx + 16].EKCOMY;
        psE[tgx].EKCOMZ += psE[tgx + 16].EKCOMZ;
      }
      if (tgx < 8) {
        psE[tgx].EKCOMX += psE[tgx + 8].EKCOMX;
        psE[tgx].EKCOMY += psE[tgx + 8].EKCOMY;
        psE[tgx].EKCOMZ += psE[tgx + 8].EKCOMZ;
      }
      if (tgx < 4) {
        psE[tgx].EKCOMX += psE[tgx + 4].EKCOMX;
        psE[tgx].EKCOMY += psE[tgx + 4].EKCOMY;
        psE[tgx].EKCOMZ += psE[tgx + 4].EKCOMZ;
      }
      if (tgx < 2) {
        psE[tgx].EKCOMX += psE[tgx + 2].EKCOMX;
        psE[tgx].EKCOMY += psE[tgx + 2].EKCOMY;
        psE[tgx].EKCOMZ += psE[tgx + 2].EKCOMZ;
      }
      if (tgx == 0) {
        psE->EKCOMX += psE[1].EKCOMX;
        psE->EKCOMY += psE[1].EKCOMY;
        psE->EKCOMZ += psE[1].EKCOMZ;
        unsigned long long int val1 = llitoulli(llrint(psE->EKCOMX * FORCESCALE));
        unsigned long long int val2 = llitoulli(llrint(psE->EKCOMY * FORCESCALE));
        unsigned long long int val3 = llitoulli(llrint(psE->EKCOMZ * FORCESCALE));
#ifdef NTP_LOTSOFMOLECULES
        if (oldMoleculeID < maxMolecules) {
          atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMX, val1);
          atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMY, val2);
          atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMZ, val3);
        }
        else {
          atomicAdd(&cSim.pSoluteUllEKCOMX[oldMoleculeID], val1);
          atomicAdd(&cSim.pSoluteUllEKCOMY[oldMoleculeID], val2);
          atomicAdd(&cSim.pSoluteUllEKCOMZ[oldMoleculeID], val3);
        }
#else  // NTP_LOTSOFMOLECULES
        atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMX, val1);
        atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMY, val2);
        atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMZ, val3);
#endif // NTP_LOTSOFMOLECULES
      }
      EX = (PMEDouble)0.0;
      EY = (PMEDouble)0.0;
      EZ = (PMEDouble)0.0;
    }
    oldMoleculeID = moleculeID;
    if (atomID != -1) {
      EX += mass * vx;
      EY += mass * vy;
      EZ += mass * vz;
    }
    pos += blockDim.x * gridDim.x;
  }

  // Dump last batch of solute data to shared memory
  psE[tgx].EKCOMX = EX;
  psE[tgx].EKCOMY = EY;
  psE[tgx].EKCOMZ = EZ;
#ifdef AMBER_PLATFORM_AMD_WARP64
  if (tgx < 32) {
    psE[tgx].EKCOMX += psE[tgx + 32].EKCOMX;
    psE[tgx].EKCOMY += psE[tgx + 32].EKCOMY;
    psE[tgx].EKCOMZ += psE[tgx + 32].EKCOMZ;
  }
#endif
  if (tgx < 16) {
    psE[tgx].EKCOMX += psE[tgx + 16].EKCOMX;
    psE[tgx].EKCOMY += psE[tgx + 16].EKCOMY;
    psE[tgx].EKCOMZ += psE[tgx + 16].EKCOMZ;
  }
  if (tgx < 8) {
    psE[tgx].EKCOMX += psE[tgx + 8].EKCOMX;
    psE[tgx].EKCOMY += psE[tgx + 8].EKCOMY;
    psE[tgx].EKCOMZ += psE[tgx + 8].EKCOMZ;
  }
  if (tgx < 4) {
    psE[tgx].EKCOMX += psE[tgx + 4].EKCOMX;
    psE[tgx].EKCOMY += psE[tgx + 4].EKCOMY;
    psE[tgx].EKCOMZ += psE[tgx + 4].EKCOMZ;
  }
  if (tgx < 2) {
    psE[tgx].EKCOMX += psE[tgx + 2].EKCOMX;
    psE[tgx].EKCOMY += psE[tgx + 2].EKCOMY;
    psE[tgx].EKCOMZ += psE[tgx + 2].EKCOMZ;
  }
  if ((tgx == 0) && (oldMoleculeID != -2)) {
    psE->EKCOMX += psE[1].EKCOMX;
    psE->EKCOMY += psE[1].EKCOMY;
    psE->EKCOMZ += psE[1].EKCOMZ;
    unsigned long long int val1 = llitoulli(llrint(psE->EKCOMX * FORCESCALE));
    unsigned long long int val2 = llitoulli(llrint(psE->EKCOMY * FORCESCALE));
    unsigned long long int val3 = llitoulli(llrint(psE->EKCOMZ * FORCESCALE));
#ifdef NTP_LOTSOFMOLECULES
    if (oldMoleculeID < maxMolecules) {
      atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMX, val1);
      atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMY, val2);
      atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMZ, val3);
    }
    else {
      atomicAdd(&cSim.pSoluteUllEKCOMX[oldMoleculeID], val1);
      atomicAdd(&cSim.pSoluteUllEKCOMY[oldMoleculeID], val2);
      atomicAdd(&cSim.pSoluteUllEKCOMZ[oldMoleculeID], val3);
    }
#else  // NTP_LOTSOFMOLECULES
    atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMX, val1);
    atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMY, val2);
    atomicAdd((PMEUllInt*)&sEKCOM[oldMoleculeID].EKCOMZ, val3);
#endif // NTP_LOTSOFMOLECULES
  }

  // Dump solute atoms to memory
  __syncthreads();
  pos = threadIdx.x;
#ifdef NTP_LOTSOFMOLECULES
  while (pos < maxMolecules) {
#else
  while (pos < cSim.soluteMolecules) {
#endif
    if (sEKCOM[pos].EKCOMX != 0) {
      atomicAdd(&cSim.pSoluteUllEKCOMX[pos], sEKCOM[pos].EKCOMX);
    }
    if (sEKCOM[pos].EKCOMY != 0) {
      atomicAdd(&cSim.pSoluteUllEKCOMY[pos], sEKCOM[pos].EKCOMY);
    }
    if (sEKCOM[pos].EKCOMZ != 0) {
      atomicAdd(&cSim.pSoluteUllEKCOMZ[pos], sEKCOM[pos].EKCOMZ);
    }
    pos += blockDim.x;
#ifdef NTP_LOTSOFMOLECULES
  }
#else
  }
#endif

  // Solvent atoms
  pos = blockIdx.x*blockDim.x + threadIdx.x;
  PMEDouble EKCOMX = (PMEDouble)0.0;
  PMEDouble EKCOMY = (PMEDouble)0.0;
  PMEDouble EKCOMZ = (PMEDouble)0.0;
  while (pos < cSim.solventMolecules) {
    int4 atomID       = cSim.pImageSolventAtomID[pos];
    PMEDouble invMass = cSim.pSolventInvMass[pos];
    PMEDouble mass    = cSim.pSolventAtomMass1[pos];
    PMEDouble vx      = cSim.pImageVelX[atomID.x];
    PMEDouble vy      = cSim.pImageVelY[atomID.x];
    PMEDouble vz      = cSim.pImageVelZ[atomID.x];
    PMEDouble EX      = mass * vx;
    PMEDouble EY      = mass * vy;
    PMEDouble EZ      = mass * vz;
    if (atomID.y != -1) {
      PMEDouble mass = cSim.pSolventAtomMass2[pos];
      PMEDouble vx   = cSim.pImageVelX[atomID.y];
      PMEDouble vy   = cSim.pImageVelY[atomID.y];
      PMEDouble vz   = cSim.pImageVelZ[atomID.y];
      EX            += mass * vx;
      EY            += mass * vy;
      EZ            += mass * vz;
    }
    if (atomID.z != -1) {
      PMEDouble mass = cSim.pSolventAtomMass3[pos];
      PMEDouble vx   = cSim.pImageVelX[atomID.z];
      PMEDouble vy   = cSim.pImageVelY[atomID.z];
      PMEDouble vz   = cSim.pImageVelZ[atomID.z];
      EX            += mass * vx;
      EY            += mass * vy;
      EZ            += mass * vz;
    }
    if (atomID.w != -1) {
      PMEDouble mass = cSim.pSolventAtomMass4[pos];
      PMEDouble vx   = cSim.pImageVelX[atomID.w];
      PMEDouble vy   = cSim.pImageVelY[atomID.w];
      PMEDouble vz   = cSim.pImageVelZ[atomID.w];
      EX            += mass * vx;
      EY            += mass * vy;
      EZ            += mass * vz;
    }

    // Sum up results
    EKCOMX += invMass * EX * EX;
    EKCOMY += invMass * EY * EY;
    EKCOMZ += invMass * EZ * EZ;
    pos    += blockDim.x * gridDim.x;
  }
  sE[threadIdx.x].EKCOMX = EKCOMX;
  sE[threadIdx.x].EKCOMY = EKCOMY;
  sE[threadIdx.x].EKCOMZ = EKCOMZ;
  __syncthreads();
  unsigned int m = 1;
  while (m < blockDim.x) {
    int p = threadIdx.x + m;

    // OPTIMIZE
    PMEDouble EX = ((p < blockDim.x) ? sE[p].EKCOMX : (PMEDouble)0.0);
    PMEDouble EY = ((p < blockDim.x) ? sE[p].EKCOMY : (PMEDouble)0.0);
    PMEDouble EZ = ((p < blockDim.x) ? sE[p].EKCOMZ : (PMEDouble)0.0);
    __syncthreads();
    sE[threadIdx.x].EKCOMX += EX;
    sE[threadIdx.x].EKCOMY += EY;
    sE[threadIdx.x].EKCOMZ += EZ;
    __syncthreads();
    m *= 2;
  }
  if (threadIdx.x == 0) {
    unsigned long long int val1 = llitoulli(llrint(sE[0].EKCOMX * FORCESCALE));
    unsigned long long int val2 = llitoulli(llrint(sE[0].EKCOMY * FORCESCALE));
    unsigned long long int val3 = llitoulli(llrint(sE[0].EKCOMZ * FORCESCALE));
    atomicAdd(cSim.pEKCOMX, val1);
    atomicAdd(cSim.pEKCOMY, val2);
    atomicAdd(cSim.pEKCOMZ, val3);
  }
}

//---------------------------------------------------------------------------------------------
// KCALCULATEMOLECULARVIRIAL_KERNEL: compute the molecular virial.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(THREADS_PER_BLOCK, 1)
KCALCULATEMOLECULARVIRIAL_KERNEL()
{
  struct COM1 {
    double COMX;
    double COMY;
    double COMZ;
  };

  struct Virial {
    PMEDouble vir_11;
    PMEDouble vir_22;
    PMEDouble vir_33;
  };

  __shared__ volatile COM1 sCOM[MAXMOLECULES];
  __shared__ volatile Virial sV[THREADS_PER_BLOCK];
#ifdef NTP_LOTSOFMOLECULES
  const int maxMolecules = MAXMOLECULES;
#endif

  // Read solute COMs into shared memory
  unsigned int pos = threadIdx.x;
  while (pos < cSim.soluteMolecules) {
#ifdef NTP_LOTSOFMOLECULES
    if (pos < maxMolecules) {
      sCOM[pos].COMX = cSim.pSoluteCOMX[pos];
      sCOM[pos].COMY = cSim.pSoluteCOMY[pos];
      sCOM[pos].COMZ = cSim.pSoluteCOMZ[pos];
    }
#else
    sCOM[pos].COMX = cSim.pSoluteCOMX[pos];
    sCOM[pos].COMY = cSim.pSoluteCOMY[pos];
    sCOM[pos].COMZ = cSim.pSoluteCOMZ[pos];
#endif
    if (blockIdx.x == 0) {
      cSim.pSoluteUllCOMX[pos] = 0;
      cSim.pSoluteUllCOMY[pos] = 0;
      cSim.pSoluteUllCOMZ[pos] = 0;
    }
    pos += blockDim.x;
  }
  __syncthreads();

  // Calculate per-thread virial
  PMEDouble vir_11 = 0.0;
  PMEDouble vir_22 = 0.0;
  PMEDouble vir_33 = 0.0;
  pos = blockIdx.x*blockDim.x + threadIdx.x;
  while (pos < cSim.soluteAtoms) {
    int atomID     = cSim.pImageSoluteAtomID[pos];
    int moleculeID = cSim.pSoluteAtomMoleculeID[pos];
    if (atomID != -1) {
      PMEDouble fx = cSim.pNBForceXAccumulator[atomID];
      PMEDouble x  = cSim.pImageX[atomID];
      PMEDouble fy = cSim.pNBForceYAccumulator[atomID];
      PMEDouble y  = cSim.pImageY[atomID];
      PMEDouble fz = cSim.pNBForceZAccumulator[atomID];
      PMEDouble z  = cSim.pImageZ[atomID];
#ifdef NTP_LOTSOFMOLECULES
      double cx, cy, cz;
      if (moleculeID < maxMolecules) {
        cx = sCOM[moleculeID].COMX;
        cy = sCOM[moleculeID].COMY;
        cz = sCOM[moleculeID].COMZ;
      }
      else {
        cx = cSim.pSoluteCOMX[moleculeID];
        cy = cSim.pSoluteCOMY[moleculeID];
        cz = cSim.pSoluteCOMZ[moleculeID];
      }
      vir_11 += fx * (x - cx);
      vir_22 += fy * (y - cy);
      vir_33 += fz * (z - cz);
#else
      vir_11 += fx * (x - sCOM[moleculeID].COMX);
      vir_22 += fy * (y - sCOM[moleculeID].COMY);
      vir_33 += fz * (z - sCOM[moleculeID].COMZ);
#endif
    }
    pos += blockDim.x * gridDim.x;
  }
  pos = blockIdx.x*blockDim.x + threadIdx.x;
  while (pos < cSim.solventMolecules) {
    int4 atomID    = cSim.pImageSolventAtomID[pos];
    PMEDouble COMX = cSim.pSolventCOMX[pos];
    PMEDouble COMY = cSim.pSolventCOMY[pos];
    PMEDouble COMZ = cSim.pSolventCOMZ[pos];
    PMEDouble fx   = cSim.pNBForceXAccumulator[atomID.x];
    PMEDouble x    = cSim.pImageX[atomID.x];
    PMEDouble fy   = cSim.pNBForceYAccumulator[atomID.x];
    PMEDouble y    = cSim.pImageY[atomID.x];
    PMEDouble fz   = cSim.pNBForceZAccumulator[atomID.x];
    PMEDouble z    = cSim.pImageZ[atomID.x];
    vir_11 += fx * (x - COMX);
    vir_22 += fy * (y - COMY);
    vir_33 += fz * (z - COMZ);
    if (atomID.y != -1) {
      fx      = cSim.pNBForceXAccumulator[atomID.y];
      x       = cSim.pImageX[atomID.y];
      fy      = cSim.pNBForceYAccumulator[atomID.y];
      y       = cSim.pImageY[atomID.y];
      fz      = cSim.pNBForceZAccumulator[atomID.y];
      z       = cSim.pImageZ[atomID.y];
      vir_11 += fx * (x - COMX);
      vir_22 += fy * (y - COMY);
      vir_33 += fz * (z - COMZ);
    }
    if (atomID.z != -1) {
      fx      = cSim.pNBForceXAccumulator[atomID.z];
      x       = cSim.pImageX[atomID.z];
      fy      = cSim.pNBForceYAccumulator[atomID.z];
      y       = cSim.pImageY[atomID.z];
      fz      = cSim.pNBForceZAccumulator[atomID.z];
      z       = cSim.pImageZ[atomID.z];
      vir_11 += fx * (x - COMX);
      vir_22 += fy * (y - COMY);
      vir_33 += fz * (z - COMZ);
    }
    if (atomID.w != -1) {
      fx      = cSim.pNBForceXAccumulator[atomID.w];
      x       = cSim.pImageX[atomID.w];
      fy      = cSim.pNBForceYAccumulator[atomID.w];
      y       = cSim.pImageY[atomID.w];
      fz      = cSim.pNBForceZAccumulator[atomID.w];
      z       = cSim.pImageZ[atomID.w];
      vir_11 += fx * (x - COMX);
      vir_22 += fy * (y - COMY);
      vir_33 += fz * (z - COMZ);
    }
    pos += blockDim.x * gridDim.x;
  }

  // Reduce virial
  sV[threadIdx.x].vir_11 = vir_11;
  sV[threadIdx.x].vir_22 = vir_22;
  sV[threadIdx.x].vir_33 = vir_33;
  __syncthreads();
  unsigned int m  = 1;
  while (m < blockDim.x) {
    int p = threadIdx.x + m;
    PMEDouble vir_11 = ((p < blockDim.x) ? sV[p].vir_11 : (PMEDouble)0.0);
    PMEDouble vir_22 = ((p < blockDim.x) ? sV[p].vir_22 : (PMEDouble)0.0);
    PMEDouble vir_33 = ((p < blockDim.x) ? sV[p].vir_33 : (PMEDouble)0.0);
    __syncthreads();
    sV[threadIdx.x].vir_11 += vir_11;
    sV[threadIdx.x].vir_22 += vir_22;
    sV[threadIdx.x].vir_33 += vir_33;
    __syncthreads();
    m *= 2;
  }
  if (threadIdx.x == 0) {
    unsigned long long int val1 = llitoulli(llrint(sV[0].vir_11));
    unsigned long long int val2 = llitoulli(llrint(sV[0].vir_22));
    unsigned long long int val3 = llitoulli(llrint(sV[0].vir_33));
    atomicAdd(cSim.pVirial_11, val1);
    atomicAdd(cSim.pVirial_22, val2);
    atomicAdd(cSim.pVirial_33, val3);
  }
}

//---------------------------------------------------------------------------------------------
// KPRESSURESCALECOORDINATES_KERNEL: scale coordinates in step with box resizing.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(THREADS_PER_BLOCK, 1)
KPRESSURESCALECOORDINATES_KERNEL()
{
  struct COM1 {
    double COMX;
    double COMY;
    double COMZ;
  };

  __shared__ volatile COM1 sCOM[MAXPSMOLECULES];
#ifdef NTP_LOTSOFMOLECULES
  const int maxMolecules = MAXPSMOLECULES;
#endif
  __shared__ volatile PMEDouble sLast_recip[9];
  __shared__ volatile PMEDouble sUcell[9];

  // Read transformation matrices
  if (threadIdx.x < 9) {
    sLast_recip[threadIdx.x] = cSim.pNTPData->last_recip[threadIdx.x];
    sUcell[threadIdx.x]      = cSim.pNTPData->ucell[threadIdx.x];
  }
  __syncthreads();

  // Read solute COMs into shared memory
  unsigned int pos = threadIdx.x;
  while (pos < cSim.soluteMolecules) {

    // Read solute COM
    PMEDouble invMass                           = cSim.pSoluteInvMass[pos];
    PMEUllInt ullCOMX                           = cSim.pSoluteUllCOMX[pos];
    PMEUllInt ullCOMY                           = cSim.pSoluteUllCOMY[pos];
    PMEUllInt ullCOMZ                           = cSim.pSoluteUllCOMZ[pos];

    invMass                                    *= ONEOVERENERGYSCALE;
    PMEDouble ox, oy, oz;
    if (ullCOMX >= 0x8000000000000000ull) {
      ox = -(PMEDouble)(ullCOMX ^ 0xffffffffffffffffull);
    }
    else {
      ox =  (PMEDouble)ullCOMX;
    }
    if (ullCOMY >= 0x8000000000000000ull) {
      oy = -(PMEDouble)(ullCOMY ^ 0xffffffffffffffffull);
    }
    else {
      oy =  (PMEDouble)ullCOMY;
    }
    if (ullCOMZ >= 0x8000000000000000ull) {
      oz = -(PMEDouble)(ullCOMZ ^ 0xffffffffffffffffull);
    }
    else {
      oz =  (PMEDouble)ullCOMZ;
    }
    ox *= invMass;
    oy *= invMass;
    oz *= invMass;

    // Calculate COM displacement
    PMEDouble cx = ox*sLast_recip[0] + oy*sLast_recip[3] + oz*sLast_recip[6];
    PMEDouble cy =                     oy*sLast_recip[4] + oz*sLast_recip[7];
    PMEDouble cz =                                         oz*sLast_recip[8];
    cx           = cx*sUcell[0] + cy*sUcell[1] + cz*sUcell[2];
    cy           =                cy*sUcell[4] + cz*sUcell[5];
    cz           =                               cz*sUcell[8];
    if (blockIdx.x == 0) {
      cSim.pSoluteCOMX[pos] = cx;
      cSim.pSoluteCOMY[pos] = cy;
      cSim.pSoluteCOMZ[pos] = cz;
    }

    cx -= ox;
    cy -= oy;
    cz -= oz;

#ifdef NTP_LOTSOFMOLECULES
    if (pos < maxMolecules) {
      sCOM[pos].COMX = cx;
      sCOM[pos].COMY = cy;
      sCOM[pos].COMZ = cz;
    }
    else {
      cSim.pSoluteDeltaCOMX[pos] = cx;
      cSim.pSoluteDeltaCOMY[pos] = cy;
      cSim.pSoluteDeltaCOMZ[pos] = cz;
    }
#else
    sCOM[pos].COMX = cx;
    sCOM[pos].COMY = cy;
    sCOM[pos].COMZ = cz;
#endif
    pos += blockDim.x;
  }
  __syncthreads();

  pos = blockIdx.x*blockDim.x + threadIdx.x;
  while (pos < cSim.soluteAtoms) {
    int atomID     = cSim.pImageSoluteAtomID[pos];
    int moleculeID = cSim.pSoluteAtomMoleculeID[pos];
    if (atomID != -1) {
      PMEDouble x = cSim.pImageX[atomID];
      PMEDouble y = cSim.pImageY[atomID];
      PMEDouble z = cSim.pImageZ[atomID];
#ifdef NTP_LOTSOFMOLECULES
      if (moleculeID < maxMolecules) {
        x += sCOM[moleculeID].COMX;
        y += sCOM[moleculeID].COMY;
        z += sCOM[moleculeID].COMZ;
      }
      else {
        x += cSim.pSoluteDeltaCOMX[moleculeID];
        y += cSim.pSoluteDeltaCOMY[moleculeID];
        z += cSim.pSoluteDeltaCOMZ[moleculeID];
      }
#else
      x += sCOM[moleculeID].COMX;
      y += sCOM[moleculeID].COMY;
      z += sCOM[moleculeID].COMZ;
#endif
      cSim.pImageX[atomID] = x;
      cSim.pImageY[atomID] = y;
      cSim.pImageZ[atomID] = z;
    }
    pos += blockDim.x * gridDim.x;
  }
  pos = blockIdx.x*blockDim.x + threadIdx.x;
  while (pos < cSim.solventMolecules) {
    int4 atomID  = cSim.pImageSolventAtomID[pos];
    PMEDouble m1 = cSim.pSolventAtomMass1[pos];
    PMEDouble invMass = cSim.pSolventInvMass[pos];
    PMEDouble x1 = cSim.pImageX[atomID.x];
    PMEDouble y1 = cSim.pImageY[atomID.x];
    PMEDouble z1 = cSim.pImageZ[atomID.x];
    PMEDouble ox = m1 * x1;
    PMEDouble oy = m1 * y1;
    PMEDouble oz = m1 * z1;
    PMEDouble x2, y2, z2;
    if (atomID.y != -1) {
      PMEDouble m2 = cSim.pSolventAtomMass2[pos];
      x2           = cSim.pImageX[atomID.y];
      y2           = cSim.pImageY[atomID.y];
      z2           = cSim.pImageZ[atomID.y];
      ox          += m2 * x2;
      oy          += m2 * y2;
      oz          += m2 * z2;
    }
    double x3, y3, z3;
    if (atomID.z != -1) {
      PMEDouble m3 = cSim.pSolventAtomMass3[pos];
      x3           = cSim.pImageX[atomID.z];
      y3           = cSim.pImageY[atomID.z];
      z3           = cSim.pImageZ[atomID.z];
      ox          += m3 * x3;
      oy          += m3 * y3;
      oz          += m3 * z3;
    }
    double x4, y4, z4;
    if (atomID.w != -1) {
      PMEDouble m4 = cSim.pSolventAtomMass4[pos];
      x4           = cSim.pImageX[atomID.w];
      y4           = cSim.pImageY[atomID.w];
      z4           = cSim.pImageZ[atomID.w];
      ox          += m4 * x4;
      oy          += m4 * y4;
      oz          += m4 * z4;
    }

    // Calculate change in COM
    ox *= invMass;
    oy *= invMass;
    oz *= invMass;
    PMEDouble cx = ox*sLast_recip[0] + oy*sLast_recip[3] + oz*sLast_recip[6];
    PMEDouble cy =                     oy*sLast_recip[4] + oz*sLast_recip[7];
    PMEDouble cz =                                         oz*sLast_recip[8];
    cx           = cx*sUcell[0] + cy*sUcell[1] + cz*sUcell[2];
    cy           =                cy*sUcell[4] + cz*sUcell[5];
    cz           =                               cz*sUcell[8];
    cSim.pSolventCOMX[pos] = cx;
    cSim.pSolventCOMY[pos] = cy;
    cSim.pSolventCOMZ[pos] = cz;
    cx -= ox;
    cy -= oy;
    cz -= oz;

    // Shift atomic coordinates
    cSim.pImageX[atomID.x] = x1 + cx;
    cSim.pImageY[atomID.x] = y1 + cy;
    cSim.pImageZ[atomID.x] = z1 + cz;
    if (atomID.y != -1) {
      cSim.pImageX[atomID.y] = x2 + cx;
      cSim.pImageY[atomID.y] = y2 + cy;
      cSim.pImageZ[atomID.y] = z2 + cz;
    }
    if (atomID.z != -1) {
      cSim.pImageX[atomID.z] = x3 + cx;
      cSim.pImageY[atomID.z] = y3 + cy;
      cSim.pImageZ[atomID.z] = z3 + cz;
    }
    if (atomID.w != -1) {
      cSim.pImageX[atomID.w] = x4 + cx;
      cSim.pImageY[atomID.w] = y4 + cy;
      cSim.pImageZ[atomID.w] = z4 + cz;
    }
    pos += blockDim.x * gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// KPMECALCULATECOMAFE_KERNEL: compute centers of mass when a soft core region is involved.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(THREADS_PER_BLOCK, 1)
KPMECALCULATECOMAFE_KERNEL()
{
  struct COM1 {
    double COMX;
    double COMY;
    double COMZ;
  };

  struct UllCOM1 {
    PMEUllInt COMX;
    PMEUllInt COMY;
    PMEUllInt COMZ;
  };
  __shared__ volatile COM1 sC[THREADS_PER_BLOCK];
  __shared__ volatile UllCOM1 sCOM[MAXMOLECULES];
#ifdef NTP_LOTSOFMOLECULES
  const int maxMolecules = MAXMOLECULES;
#endif
  unsigned int pos;

  // Clear shared memory COM
  pos = threadIdx.x;
#ifdef NTP_LOTSOFMOLECULES
  while (pos < maxMolecules) {
#else
  while (pos < cSim.soluteMolecules) {
#endif
    sCOM[pos].COMX  = (PMEUllInt)0;
    sCOM[pos].COMY  = (PMEUllInt)0;
    sCOM[pos].COMZ  = (PMEUllInt)0;
    pos            += blockDim.x;
#ifdef NTP_LOTSOFMOLECULES
  }
#else
  }
#endif
  __syncthreads();

  // Calculate solute center of mass
  pos = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int tgx   = threadIdx.x & (GRID - 1);
  volatile COM1* psC = &sC[threadIdx.x - tgx];
  int oldMoleculeID = -2;
  int atomID, moleculeID;
  PMEDouble CX = (PMEDouble)0.0;
  PMEDouble CY = (PMEDouble)0.0;
  PMEDouble CZ = (PMEDouble)0.0;
  PMEDouble SC_CY_region2 = (PMEDouble)0.0;
  PMEDouble SC_CZ_region2 = (PMEDouble)0.0;
  PMEDouble SC_CX_region2 = (PMEDouble)0.0;

  // Determine the ti_mol type for a given molecule
  //   0   : common region, do everything normally
  //   1,2 : fully sc (region 1,2), add the two soft core COMs together.
  //   3   : partially sc, need to average COM.
  //   4   : partially sc, need to find partner mol and average COM
  while (pos < cSim.soluteAtoms) {
    atomID     = cSim.pImageSoluteAtomID[pos];
    moleculeID = cSim.pSoluteAtomMoleculeID[pos];
    PMEDouble mass = 0;
    PMEDouble x, y, z;

    // Only have 1 entry per molecule for molecule type, but
    // can be grabbed by all threads with an atom in the molecule
    int AFE_mol_type = cSim.pAFEMolType[moleculeID];

    // Skip COM calculation for this molecule
    bool bFullySC    = ((AFE_mol_type + 1) >> 1) & 0x1;

    // Tells us to average the COM
    bool bPartSC     = (AFE_mol_type > 2);

    // If TRUE, need to find the partner
    bool bNeedExchg  = (!bFullySC) & (AFE_mol_type >> 2);

    // For bPartSC true, add to first com, second com, or both
    int region       = cSim.pImageTIRegion[pos];
    if (bPartSC) {
      if (bNeedExchg) {
        moleculeID = cSim.pAFEMolPartner[moleculeID];
      }
    }
    if (!bFullySC) {

      // Skip COM calculation for totally softcore molecules.  Because we
      // don't calculate COM for atoms with no mass, just leave the mass at 0.0.
      if (atomID != -1) {
        mass = cSim.pSoluteAtomMass[pos];
        x    = cSim.pImageX[atomID];
        y    = cSim.pImageY[atomID];
        z    = cSim.pImageZ[atomID];
      }
    }

    // Output COM upon changed status
    if (moleculeID != oldMoleculeID) {
      psC[tgx].COMX = CX;
      psC[tgx].COMY = CY;
      psC[tgx].COMZ = CZ;
      if (oldMoleculeID >= 0) {
#ifdef AMBER_PLATFORM_AMD_WARP64
        if (tgx < 32) {
          psC[tgx].COMX += psC[tgx + 32].COMX;
          psC[tgx].COMY += psC[tgx + 32].COMY;
          psC[tgx].COMZ += psC[tgx + 32].COMZ;
        }
#endif
        if (tgx < 16) {
          psC[tgx].COMX += psC[tgx + 16].COMX;
          psC[tgx].COMY += psC[tgx + 16].COMY;
          psC[tgx].COMZ += psC[tgx + 16].COMZ;
        }
        if (tgx < 8) {
          psC[tgx].COMX += psC[tgx + 8].COMX;
          psC[tgx].COMY += psC[tgx + 8].COMY;
          psC[tgx].COMZ += psC[tgx + 8].COMZ;
        }
        if (tgx < 4) {
          psC[tgx].COMX += psC[tgx + 4].COMX;
          psC[tgx].COMY += psC[tgx + 4].COMY;
          psC[tgx].COMZ += psC[tgx + 4].COMZ;
        }
        if (tgx < 2) {
          psC[tgx].COMX += psC[tgx + 2].COMX;
          psC[tgx].COMY += psC[tgx + 2].COMY;
          psC[tgx].COMZ += psC[tgx + 2].COMZ;
        }
        if (tgx == 0) {
          psC->COMX += psC[1].COMX;
          psC->COMY += psC[1].COMY;
          psC->COMZ += psC[1].COMZ;
          unsigned long long int val1 = llitoulli(llrint(psC->COMX * ENERGYSCALE));
          unsigned long long int val2 = llitoulli(llrint(psC->COMY * ENERGYSCALE));
          unsigned long long int val3 = llitoulli(llrint(psC->COMZ * ENERGYSCALE));

#ifdef NTP_LOTSOFMOLECULES
          if (oldMoleculeID < maxMolecules) {
            atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
            atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
            atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
          }
          else {
            atomicAdd(&cSim.pSoluteUllCOMX[oldMoleculeID], val1);
            atomicAdd(&cSim.pSoluteUllCOMY[oldMoleculeID], val2);
            atomicAdd(&cSim.pSoluteUllCOMZ[oldMoleculeID], val3);
          }
#else  // NTP_LOTSOFMOLECULES
          atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
          atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
          atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
#endif // NTP_LOTSOFMOLECULES
        }
      }
      CX                                      = (PMEDouble)0.0;
      CY                                      = (PMEDouble)0.0;
      CZ                                      = (PMEDouble)0.0;
    }
    // End contingency for oldMoleculeID != moleculeID

    oldMoleculeID = moleculeID;
    if (mass != 0.0) {
      if (bPartSC) {
        if (region != 2) {

          // Add to region 1
          CX += mass * x / 2.0;
          CY += mass * y / 2.0;
          CZ += mass * z / 2.0;
        }
        if (region != 1) {

          // Add to region 2.  This way, common adds to both region 1 and 2.
          SC_CX_region2 += mass * x / 2.0;
          SC_CY_region2 += mass * y / 2.0;
          SC_CZ_region2 += mass * z / 2.0;
        }
        CX += SC_CX_region2;
        CY += SC_CY_region2;
        CZ += SC_CZ_region2;
      }
      if (!bPartSC) {
        CX += mass * x;
        CY += mass * y;
        CZ += mass * z;
      }
    }
    pos += blockDim.x * gridDim.x;
  }

  // Dump last batch of solute data to shared memory
  psC[tgx].COMX = CX;
  psC[tgx].COMY = CY;
  psC[tgx].COMZ = CZ;
  if (oldMoleculeID >= 0) {
#ifdef AMBER_PLATFORM_AMD_WARP64
    if (tgx < 32) {
      psC[tgx].COMX += psC[tgx + 32].COMX;
      psC[tgx].COMY += psC[tgx + 32].COMY;
      psC[tgx].COMZ += psC[tgx + 32].COMZ;
    }
#endif
    if (tgx < 16) {
      psC[tgx].COMX += psC[tgx + 16].COMX;
      psC[tgx].COMY += psC[tgx + 16].COMY;
      psC[tgx].COMZ += psC[tgx + 16].COMZ;
    }
    if (tgx < 8) {
      psC[tgx].COMX += psC[tgx + 8].COMX;
      psC[tgx].COMY += psC[tgx + 8].COMY;
      psC[tgx].COMZ += psC[tgx + 8].COMZ;
    }
    if (tgx < 4) {
      psC[tgx].COMX += psC[tgx + 4].COMX;
      psC[tgx].COMY += psC[tgx + 4].COMY;
      psC[tgx].COMZ += psC[tgx + 4].COMZ;
    }
    if (tgx < 2) {
      psC[tgx].COMX += psC[tgx + 2].COMX;
      psC[tgx].COMY += psC[tgx + 2].COMY;
      psC[tgx].COMZ += psC[tgx + 2].COMZ;
    }
    if (tgx == 0) {
      psC->COMX += psC[1].COMX;
      psC->COMY += psC[1].COMY;
      psC->COMZ += psC[1].COMZ;
      unsigned long long int val1 = llitoulli(llrint(psC->COMX * ENERGYSCALE));
      unsigned long long int val2 = llitoulli(llrint(psC->COMY * ENERGYSCALE));
      unsigned long long int val3 = llitoulli(llrint(psC->COMZ * ENERGYSCALE));

#ifdef NTP_LOTSOFMOLECULES
      if (oldMoleculeID < maxMolecules) {
        atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
        atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
        atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
      }
      else {
        atomicAdd(&cSim.pSoluteUllCOMX[oldMoleculeID], val1);
        atomicAdd(&cSim.pSoluteUllCOMY[oldMoleculeID], val2);
        atomicAdd(&cSim.pSoluteUllCOMZ[oldMoleculeID], val3);
      }
#else  // NTP_LOTSOFMOLECULES
      atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMX, val1);
      atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMY, val2);
      atomicAdd((PMEUllInt*)&sCOM[oldMoleculeID].COMZ, val3);
#endif // NTP_LOTSOFMOLECULES
    }
  }

  // Dump solute data to main memory
  __syncthreads();
  pos = threadIdx.x;
#ifdef NTP_LOTSOFMOLECULES
  while (pos < maxMolecules) {
#else
  while (pos < cSim.soluteMolecules) {
#endif
    if (sCOM[pos].COMX != 0) {
      atomicAdd(&cSim.pSoluteUllCOMX[pos], sCOM[pos].COMX);
    }
    if (sCOM[pos].COMY != 0) {
      atomicAdd(&cSim.pSoluteUllCOMY[pos], sCOM[pos].COMY);
    }
    if (sCOM[pos].COMZ != 0) {
      atomicAdd(&cSim.pSoluteUllCOMZ[pos], sCOM[pos].COMZ);
    }
    pos += blockDim.x;
#ifdef NTP_LOTSOFMOLECULES
  }
#else
  }
#endif

  // Solvent atoms: no adjustment for AFE because, for now,
  // solvent molecules cannot be part of the soft core region
  pos = blockIdx.x*blockDim.x + threadIdx.x;
  while (pos < cSim.solventMolecules) {
    int4 atomID    = cSim.pImageSolventAtomID[pos];
    PMEDouble invMass = cSim.pSolventInvMass[pos];
    PMEDouble mass = cSim.pSolventAtomMass1[pos];
    PMEDouble x    = cSim.pImageX[atomID.x];
    PMEDouble y    = cSim.pImageY[atomID.x];
    PMEDouble z    = cSim.pImageZ[atomID.x];
    PMEDouble CX   = mass * x;
    PMEDouble CY   = mass * y;
    PMEDouble CZ   = mass * z;
    if (atomID.y != -1) {
      PMEDouble mass = cSim.pSolventAtomMass2[pos];
      PMEDouble x    = cSim.pImageX[atomID.y];
      PMEDouble y    = cSim.pImageY[atomID.y];
      PMEDouble z    = cSim.pImageZ[atomID.y];
      CX            += mass * x;
      CY            += mass * y;
      CZ            += mass * z;
    }
    if (atomID.z != -1) {
      PMEDouble mass = cSim.pSolventAtomMass3[pos];
      PMEDouble x    = cSim.pImageX[atomID.z];
      PMEDouble y    = cSim.pImageY[atomID.z];
      PMEDouble z    = cSim.pImageZ[atomID.z];
      CX            += mass * x;
      CY            += mass * y;
      CZ            += mass * z;
    }
    if (atomID.w != -1) {
      PMEDouble mass = cSim.pSolventAtomMass4[pos];
      PMEDouble x    = cSim.pImageX[atomID.w];
      PMEDouble y    = cSim.pImageY[atomID.w];
      PMEDouble z    = cSim.pImageZ[atomID.w];
      CX            += mass * x;
      CY            += mass * y;
      CZ            += mass * z;
    }

    // Sum up results
    cSim.pSolventCOMX[pos]  = invMass * CX;
    cSim.pSolventCOMY[pos]  = invMass * CY;
    cSim.pSolventCOMZ[pos]  = invMass * CZ;
    pos                    += blockDim.x * gridDim.x;
  }
}
