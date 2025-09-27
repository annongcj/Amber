#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// This is included by kneighborList.cu for placing ghost atoms.  Defines include MAP_ATOM_IDS.
//---------------------------------------------------------------------------------------------
{
#if defined(MAP_ATOM_IDS) && defined(PME_VIRIAL)
  __shared__ PMEFloat sRecipf[9];
  if (threadIdx.x < 9) {
    sRecip[threadIdx.x] = cSim.pNTPData->recipf[threadIdx.x];
  }
#endif
  // Set the counter for non-bonded Neutral-Territory domain regions.
  // The non-bonded Honeycomb kernel will increment this as needed to
  // get new assignments from the list.
  if (blockIdx.x == 0 && threadIdx.x == 0) {
    cSim.pFrcBlkCounters[7] = 2 * cSim.SMPCount;
  }
  int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  while (pos < cSim.atoms) {
    PMEFloat2 xy    = cSim.pHcmbXYSP[pos];
    PMEFloat z      = cSim.pHcmbZSP[pos];
    PMEFloat2 qljid = cSim.pHcmbChargeSPLJID[pos];
    int2 filter     = cSim.pSubImageFilter[pos];
    PMEFloat4 commit;
    commit.x = xy.x;
    commit.y = xy.y;
    commit.z = z;
#ifdef MAP_ATOM_IDS
#  ifdef PME_VIRIAL
    PMEFloat fx = (sRecipf[0] * xy.x) + (sRecipf[3] * xy.y) + (sRecipf[6] * z);
#  else
    PMEFloat fx = (cSim.recipf[0][0] * xy.x) + (cSim.recipf[1][0] * xy.y) +
                  (cSim.recipf[2][0] * z);
#  endif
#endif
    if (filter.x >= 0) {
      commit.w = qljid.x;
      cSim.pHcmbQQParticles[filter.x]  = commit;
#ifdef MAP_ATOM_IDS
      cSim.pHcmbQQIdentities[filter.x] = cSim.pHcmbIndex[pos];
      cSim.pHcmbQQFractX[filter.x]     = fx;
#endif
    }
    if (filter.y >= 0) {
      commit.w = qljid.y;
      cSim.pHcmbLJParticles[filter.y]  = commit;
#ifdef MAP_ATOM_IDS
      cSim.pHcmbLJIdentities[filter.y] = cSim.pHcmbIndex[pos];
      cSim.pHcmbLJFractX[filter.y]     = fx;
#endif
    }

    // Increment the position counter
    pos += blockDim.x*gridDim.x;
  }
}
