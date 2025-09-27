#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// This file is included by kNeighborList.cu multiple times with different pre-processor
// definitions (PME_VIRIAL, PME_IS_ORTHOGONAL, PME_ATOMS_PER_WARP) to generate multiple
// implementations of kNLBuildNeighborList_????_kernel().
//---------------------------------------------------------------------------------------------
{
#if !defined(AMBER_PLATFORM_AMD)
#  define VOLATILE volatile
#else
#  define VOLATILE
#endif
  struct BNLAtom {
    float x;     //
    float y;     // Coordinates of the atom
    float z;     //
  };

  union uintfloat {
    uint u;      // Makes it possible to read a floating point number
    float f;     //   as an unsigned integer
  };

  struct BNLWarp {
    uint offset;
    uint atomList[GRID];
    PMEMask exclusionMask[GRID];
    NLEntry nlEntry;
  };

  uint2 shCell;
  uint shCellID;
  VOLATILE uintfloat shNlRecord;
#if (PME_ATOMS_PER_WARP == 32)
  const int THREADS_PER_BLOCK = NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK;
#elif (PME_ATOMS_PER_WARP == 16)
  const int THREADS_PER_BLOCK = NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK;
#else
  const int THREADS_PER_BLOCK = NLBUILD_NEIGHBORLIST8_THREADS_PER_BLOCK;
#endif
  __shared__ VOLATILE BNLWarp sNLWarp[THREADS_PER_BLOCK / GRID];
#ifdef PME_VIRIAL
  __shared__ float sCutPlusSkin2;
  __shared__ float sUcellf[9];
#  ifdef PME_IS_ORTHOGONAL
  __shared__ float3 sCellOffset[NEIGHBOR_CELLS];
#  endif
#endif
#ifdef PME_VIRIAL
  if (threadIdx.x < 9) {
    sUcellf[threadIdx.x] = cSim.pNTPData->ucellf[threadIdx.x];
  }
  if (threadIdx.x == 32) {
    sCutPlusSkin2 = cSim.pNTPData->cutPlusSkin2;
  }
#ifdef PME_IS_ORTHOGONAL
  __syncthreads();
  if (threadIdx.x < NEIGHBOR_CELLS) {
    sCellOffset[threadIdx.x].x = sUcellf[0] * cSim.cellOffset[threadIdx.x][0];
    sCellOffset[threadIdx.x].y = sUcellf[4] * cSim.cellOffset[threadIdx.x][1];
    sCellOffset[threadIdx.x].z = sUcellf[8] * cSim.cellOffset[threadIdx.x][2];
  }
#endif
  __syncthreads();
  float cutPlusSkin2 = sCutPlusSkin2;
#else
  float cutPlusSkin2 = cSim.cutPlusSkin2;
#endif
  unsigned int warp = threadIdx.x >> GRID_BITS;
  VOLATILE BNLWarp* psWarp = &sNLWarp[warp];
  unsigned int globalWarp = warp + ((blockIdx.x * blockDim.x) >> GRID_BITS);
  unsigned int* psExclusion = &cSim.pBNLExclusionBuffer[globalWarp * cSim.NLExclusionBufferSize];
  if ((threadIdx.x & GRID_BITS_MASK) == (NEIGHBOR_CELLS + 2)) {
    shNlRecord.u = (blockIdx.x * blockDim.x + threadIdx.x) >> GRID_BITS;
  }
  PMEMask exclusionMask =
    cSim.NLAtomsPerWarpMask >> (threadIdx.x & cSim.NLAtomsPerWarpBitsMask);
#if (PME_ATOMS_PER_WARP == 16)
  exclusionMask = exclusionMask | (exclusionMask << 16);
#elif (PME_ATOMS_PER_WARP == 8)
  exclusionMask = exclusionMask | (exclusionMask << 8) |
                  (exclusionMask << 16) | (exclusionMask << 24);
#endif
#ifdef AMBER_PLATFORM_AMD_WARP64
  exclusionMask = exclusionMask | (exclusionMask << 32);
#endif

  // This is the beginning of a loop that extends to the end of the library.
  // The union shNlRecord is used for various pieces of information, both floats
  // and unsigned integers, and in only one critical case is it ever used to
  // interpret one as the other.
  while (__SHFL(WARP_MASK, shNlRecord.u, NEIGHBOR_CELLS + 2) < cSim.NLRecords) {
    unsigned int tgx = threadIdx.x & GRID_BITS_MASK;

    // Read NLRecord information
    unsigned int pos1 = __SHFL(WARP_MASK, shNlRecord.u, NEIGHBOR_CELLS + 2);
    if (tgx < NEIGHBOR_CELLS + 2) {
      shNlRecord.u = cSim.pNLRecord[pos1].array[tgx];
    }

    // Calculate Exclusion/neighbor list space required
    int atomOffset = __SHFL(WARP_MASK, shNlRecord.u, NEIGHBOR_CELLS + 1) >> NLRECORD_YOFFSET_SHIFT;
    uint2 homeCell = cSim.pNLNonbondCellStartEnd[__SHFL(WARP_MASK, shNlRecord.u, NEIGHBOR_CELLS)];

    if (tgx == (NEIGHBOR_CELLS + 4)) {
      shNlRecord.u = homeCell.x;
    }
    if (tgx == (NEIGHBOR_CELLS + 5)) {
      shNlRecord.u = homeCell.y;
    }
    int ysize = max(0, (int)(__SHFL(WARP_MASK, shNlRecord.u, NEIGHBOR_CELLS + 5) -
                             __SHFL(WARP_MASK, shNlRecord.u, NEIGHBOR_CELLS + 4) -
                             atomOffset * cSim.NLAtomsPerWarp));
    if (ysize > 0) {
      ysize = 1 + max(0, ysize - 1) / (cSim.NLAtomsPerWarp * cSim.NLYDivisor);
    }

    // Calculate maximum required space: shNlRecord's unsigned int has been assigned to
    // cSim.pNLRecord[pos1].array[tgx] for tgx < 16, homeCell.x for tgx == 18, and homeCell.y
    // for tgx == 19.  The value broadcast for setting cells on all threads is going to be
    // for tgx == 15, which is the final cell with which .
    unsigned int cells = __SHFL(WARP_MASK, shNlRecord.u, NEIGHBOR_CELLS + 1) & NLRECORD_CELL_COUNT_MASK;
    PMEMask mask1 = __BALLOT(WARP_MASK, (tgx < cells) );
    if (tgx < cells) {
      uint2 cell = cSim.pNLNonbondCellStartEnd[__SHFL(mask1, shNlRecord.u, tgx) >>
                                               NLRECORD_CELL_SHIFT];
      shCell.x = cell.x;
      shCell.y = cell.y;
      psWarp->atomList[tgx] = cell.y - cell.x;
    }
    else {
      psWarp->atomList[tgx] = 0;
    }
    __SYNCWARP(WARP_MASK);
    // Reduce xsize down to thread 0 of each warp
    uint temp;
    temp = psWarp->atomList[tgx];
#if (NEIGHBOR_CELLS > 16)
    temp += __SHFL_DOWN(WARP_MASK, temp, 16);
#endif
    for (int offset = 8; offset > 0; offset /= 2)
      temp += __SHFL_DOWN(WARP_MASK, temp, offset);
    psWarp->atomList[tgx]=temp;
    __SYNCWARP(WARP_MASK);

    if (tgx == (NEIGHBOR_CELLS + 2)) {
      uint totalXSize = ((psWarp->atomList[0] + GRID - 1) >> GRID_BITS);
      uint offset = atomicAdd(cSim.pNLTotalOffset,
                              totalXSize*ysize*cSim.NLOffsetPerWarp + cSim.NLAtomsPerWarp);
      psWarp->offset = offset;
    }
    __SYNCWARP(WARP_MASK);

    // Generate actual neighbor list entry
    uint ypos = __SHFL(WARP_MASK, shNlRecord.u, NEIGHBOR_CELLS + 4) +
                ((__SHFL(WARP_MASK, shNlRecord.u,
                         NEIGHBOR_CELLS + 1) >> NLRECORD_YOFFSET_SHIFT) * cSim.NLAtomsPerWarp);
    uint homecellY = __SHFL(WARP_MASK, shNlRecord.u, NEIGHBOR_CELLS + 5);
    while ( ypos < homecellY ) {

      // Calculate y bounds and set to calculate homecell interaction
      uint ymax = min(ypos + cSim.NLAtomsPerWarp, homecellY);
      if (tgx == 0) {
        psWarp->nlEntry.NL.ypos = ypos;
        psWarp->nlEntry.NL.ymax = (ymax << NLENTRY_YMAX_SHIFT) | NLENTRY_HOME_CELL_MASK;
        psWarp->nlEntry.NL.xatoms = 0;
        psWarp->nlEntry.NL.offset = psWarp->offset;
      }
      __SYNCWARP(WARP_MASK);

      // Read y atoms
      float xi;
      float yi;
      float zi;
      unsigned int index = ypos + (tgx & (cSim.NLAtomsPerWarpBitsMask));
      if (index < ymax) {
        PMEFloat2 xy = cSim.pAtomXYSP[index];
        zi = cSim.pAtomZSP[index];
        xi = xy.x;
        yi = xy.y;
      }
      else {
        xi = (float)10000.0 * index;
        yi = (float)10000.0 * index;
        zi = (float)10000.0 * index;
      }

#ifndef PME_IS_ORTHOGONAL
      // Transform into cartesian space
#ifdef PME_VIRIAL
      xi = sUcellf[0]*xi + sUcellf[1]*yi + sUcellf[2]*zi;
      yi = sUcellf[4]*yi + sUcellf[5]*zi;
      zi = sUcellf[8]*zi;
#else
      xi = cSim.ucellf[0][0]*xi + cSim.ucellf[0][1]*yi + cSim.ucellf[0][2]*zi;
      yi = cSim.ucellf[1][1]*yi + cSim.ucellf[1][2]*zi;
      zi = cSim.ucellf[2][2]*zi;
#endif
#endif

      // Calculate bounding box on SM 2.0 and up
      float bmin = (index < ymax) ? 0.5f * xi :  999999.0f;
      float bmax = (index < ymax) ? 0.5f * xi : -999999.0f;
      bmin = min(__SHFL(WARP_MASK, bmin, tgx ^ 1), bmin);
      bmin = min(__SHFL(WARP_MASK, bmin, tgx ^ 2), bmin);
      bmin = min(__SHFL(WARP_MASK, bmin, tgx ^ 4), bmin);
#if (PME_ATOMS_PER_WARP >= 16)
      bmin = min(__SHFL(WARP_MASK, bmin, tgx ^ 8), bmin);
#endif
#if (PME_ATOMS_PER_WARP == 32)
      bmin = min(__SHFL(WARP_MASK, bmin, tgx ^ 16), bmin);
#endif
      bmax = max(__SHFL(WARP_MASK, bmax, tgx ^ 1), bmax);
      bmax = max(__SHFL(WARP_MASK, bmax, tgx ^ 2), bmax);
      bmax = max(__SHFL(WARP_MASK, bmax, tgx ^ 4), bmax);

#if (PME_ATOMS_PER_WARP >= 16)
      bmax = max(__SHFL(WARP_MASK, bmax, tgx ^ 8), bmax);
#endif
#if (PME_ATOMS_PER_WARP == 32)
      bmax = max(__SHFL(WARP_MASK, bmax, tgx ^ 16), bmax);
#endif
      if (tgx == NEIGHBOR_CELLS + 6) {
        shNlRecord.f = bmax + bmin;
      }
      if (tgx == NEIGHBOR_CELLS + 7) {
       shNlRecord.f = bmax - bmin;
      }
      bmin = (index < ymax) ? 0.5f * yi :  999999.0f;
      bmax = (index < ymax) ? 0.5f * yi : -999999.0f;
      bmin = min(__SHFL(WARP_MASK, bmin, tgx ^ 1), bmin);
      bmin = min(__SHFL(WARP_MASK, bmin, tgx ^ 2), bmin);
      bmin = min(__SHFL(WARP_MASK, bmin, tgx ^ 4), bmin);
#if (PME_ATOMS_PER_WARP >= 16)
      bmin = min(__SHFL(WARP_MASK, bmin, tgx ^ 8), bmin);
#endif
#if (PME_ATOMS_PER_WARP == 32)
      bmin = min(__SHFL(WARP_MASK, bmin, tgx ^ 16), bmin);
#endif
      bmax = max(__SHFL(WARP_MASK, bmax, tgx ^ 1), bmax);
      bmax = max(__SHFL(WARP_MASK, bmax, tgx ^ 2), bmax);
      bmax = max(__SHFL(WARP_MASK, bmax, tgx ^ 4), bmax);
#if (PME_ATOMS_PER_WARP >= 16)
      bmax = max(__SHFL(WARP_MASK, bmax, tgx ^ 8), bmax);
#endif
#if (PME_ATOMS_PER_WARP == 32)
      bmax = max(__SHFL(WARP_MASK, bmax, tgx ^ 16), bmax);
#endif
      if (tgx == NEIGHBOR_CELLS + 8) {
        shNlRecord.f = bmax + bmin;
      }
      if (tgx == NEIGHBOR_CELLS + 9) {
        shNlRecord.f = bmax - bmin;
      }
      bmin = (index < ymax) ? 0.5f * zi :  999999.0f;
      bmax = (index < ymax) ? 0.5f * zi : -999999.0f;
      bmin = min(__SHFL(WARP_MASK, bmin, tgx ^ 1), bmin);
      bmin = min(__SHFL(WARP_MASK, bmin, tgx ^ 2), bmin);
      bmin = min(__SHFL(WARP_MASK, bmin, tgx ^ 4), bmin);
#if (PME_ATOMS_PER_WARP >= 16)
      bmin = min(__SHFL(WARP_MASK, bmin, tgx ^ 8), bmin);
#endif
#if (PME_ATOMS_PER_WARP == 32)
      bmin = min(__SHFL(WARP_MASK, bmin, tgx ^ 16), bmin);
#endif
      bmax = max(__SHFL(WARP_MASK, bmax, tgx ^ 1), bmax);
      bmax = max(__SHFL(WARP_MASK, bmax, tgx ^ 2), bmax);
      bmax = max(__SHFL(WARP_MASK, bmax, tgx ^ 4), bmax);
#if (PME_ATOMS_PER_WARP >= 16)
      bmax = max(__SHFL(WARP_MASK, bmax, tgx ^ 8), bmax);
#endif
#if (PME_ATOMS_PER_WARP == 32)
      bmax = max(__SHFL(WARP_MASK, bmax, tgx ^ 16), bmax);
#endif
      if (tgx == NEIGHBOR_CELLS + 10) {
        shNlRecord.f = bmax + bmin;
      }
      if (tgx == NEIGHBOR_CELLS + 11) {
        shNlRecord.f = bmax - bmin;
      }

      // Read exclusions into L1 or shared memory
      if (tgx < ymax - ypos) {
        uint atom = cSim.pImageAtom[index];
        uint2 exclusionStartCount = cSim.pNLExclusionStartCount[atom];
        psWarp->atomList[tgx] = exclusionStartCount.x;
        psWarp->exclusionMask[tgx] = exclusionStartCount.y;
      }
      else {
        psWarp->atomList[tgx] = 0;
        psWarp->exclusionMask[tgx] = 0;
      }
      __SYNCWARP(WARP_MASK);

      uint totalExclusions = 0;
      unsigned int minExclusion = cSim.atoms;
      unsigned int maxExclusion = 0;
      uint limit = ymax - ypos;
      for (int i = 0; i < limit; i++) {
        uint start = psWarp->atomList[i];
        uint count = psWarp->exclusionMask[i];
        for (int j = tgx; j < count; j += GRID) {
          int atom = cSim.pNLExclusionList[start + j];
          int imageAtom = cSim.pImageAtomLookup[atom];
          minExclusion = min(minExclusion, imageAtom);
          maxExclusion = max(maxExclusion, imageAtom);
          psExclusion[totalExclusions + j] = (imageAtom << NLEXCLUSION_SHIFT) | i;
        }
        totalExclusions += count;
      }
      minExclusion = min(__SHFL(WARP_MASK, minExclusion, tgx ^ 1), minExclusion);
      minExclusion = min(__SHFL(WARP_MASK, minExclusion, tgx ^ 2), minExclusion);
      minExclusion = min(__SHFL(WARP_MASK, minExclusion, tgx ^ 4), minExclusion);
      minExclusion = min(__SHFL(WARP_MASK, minExclusion, tgx ^ 8), minExclusion);
      minExclusion = min(__SHFL(WARP_MASK, minExclusion, tgx ^ 16), minExclusion);
      maxExclusion = max(__SHFL(WARP_MASK, maxExclusion, tgx ^ 1), maxExclusion);
      maxExclusion = max(__SHFL(WARP_MASK, maxExclusion, tgx ^ 2), maxExclusion);
      maxExclusion = max(__SHFL(WARP_MASK, maxExclusion, tgx ^ 4), maxExclusion);
      maxExclusion = max(__SHFL(WARP_MASK, maxExclusion, tgx ^ 8), maxExclusion);
      maxExclusion = max(__SHFL(WARP_MASK, maxExclusion, tgx ^ 16), maxExclusion);
#ifdef AMBER_PLATFORM_AMD_WARP64
      minExclusion = min(__SHFL(WARP_MASK, minExclusion, tgx ^ 32), minExclusion);
      maxExclusion = max(__SHFL(WARP_MASK, maxExclusion, tgx ^ 32), maxExclusion);
#endif
      if (tgx == NEIGHBOR_CELLS + 12) {
        shNlRecord.u = minExclusion;
      }
      if (tgx == NEIGHBOR_CELLS + 13) {
        shNlRecord.u = maxExclusion;
      }

      // Initialize Neighbor List variables for current line of entry
      __SYNCWARP(0xFFFFFFFF);
      psWarp->exclusionMask[tgx] = 0;
      __SYNCWARP(WARP_MASK);
      unsigned int cpos = 0;
      unsigned int atoms = 0;
      uint minAtom = cSim.atoms;
      uint maxAtom = 0;
      PMEMask threadmask = (PMEMask)1 << tgx;

      while (cpos < cells) {

        // Check for home cell
        shCellID = __SHFL(WARP_MASK, shNlRecord.u, cpos) & NLRECORD_CELL_TYPE_MASK;
        uint xpos;

        // Cell 0 always starts along force matrix diagonal
        if ((cpos == 0) && (psWarp->nlEntry.NL.ymax & NLENTRY_HOME_CELL_MASK)) {
          // Calculate exclusions assuming all atoms are in range of each other
          psWarp->exclusionMask[tgx] = 0;
          __SYNCWARP(WARP_MASK);
          for (int i = tgx; i < totalExclusions; i += GRID) {
            uint atom = psExclusion[i] >> NLEXCLUSION_SHIFT;
            if ((atom >= ypos) && (atom < ymax)) {
              unsigned int pos = atom - ypos;
              atomicOr((PMEMask*)&(psWarp->exclusionMask[pos]),
                       ((PMEMask)1 << (psExclusion[i] & NLEXCLUSION_ATOM_MASK)));
            }
          }
          __SYNCWARP(WARP_MASK);

          // Output exclusion masks
          if (tgx < cSim.NLAtomsPerWarp) {
            PMEMask mask = psWarp->exclusionMask[tgx];
            mask = ((mask >> (1 + tgx)) | (mask << (cSim.NLAtomsPerWarp - tgx - 1))) &
                   cSim.NLAtomsPerWarpMask;
            cSim.pNLAtomList[psWarp->offset + tgx] = mask;
          }
          __SYNCWARP(WARP_MASK);
          if (tgx == 0)
            psWarp->offset += cSim.NLAtomsPerWarp;
          __SYNCWARP(WARP_MASK);
          xpos = ypos + cSim.NLAtomsPerWarp;
        }
        else {
          xpos = __SHFL(WARP_MASK, shCell.x, cpos);
        }

        // Read x atoms
        unsigned int CELLY_cpos = __SHFL(WARP_MASK, shCell.y, cpos);
        while (xpos < CELLY_cpos) {

          // Calculate number of atoms in this iteration
          uint xmax = min(xpos + GRID, CELLY_cpos) - xpos;
          float sAtomx;
          float sAtomy;
          float sAtomz;

          // Read up to GRID atoms
          if (tgx < xmax) {
            PMEFloat2 xy = cSim.pAtomXYSP[xpos + tgx];
            sAtomz = cSim.pAtomZSP[xpos + tgx];
            sAtomx = xy.x;
            sAtomy = xy.y;
          }
          else {
            sAtomx = (float)-10000.0 * tgx;
            sAtomy = (float)-10000.0 * tgx;
            sAtomz = (float)-10000.0 * tgx;
          }

          // Translate all atoms into a local coordinate system within one unit
          // cell of the first atom read to avoid PBC handling within inner loops
          unsigned int cellID = shCellID;
#if defined(PME_VIRIAL) && defined(PME_IS_ORTHOGONAL)
          sAtomx += sCellOffset[cellID].x;
          sAtomy += sCellOffset[cellID].y;
          sAtomz += sCellOffset[cellID].z;
#else
          sAtomx += cSim.cellOffset[cellID][0];
          sAtomy += cSim.cellOffset[cellID][1];
          sAtomz += cSim.cellOffset[cellID][2];
#endif
#ifndef PME_IS_ORTHOGONAL
#  ifdef PME_VIRIAL
          sAtomx = sUcellf[0]*sAtomx + sUcellf[1]*sAtomy + sUcellf[2]*sAtomz;
          sAtomy = sUcellf[4]*sAtomy + sUcellf[5]*sAtomz;
          sAtomz = sUcellf[8]*sAtomz;
#  else
          sAtomx = cSim.ucellf[0][0]*sAtomx + cSim.ucellf[0][1]*sAtomy +
                   cSim.ucellf[0][2]*sAtomz;
          sAtomy = cSim.ucellf[1][1]*sAtomy + cSim.ucellf[1][2]*sAtomz;
          sAtomz = cSim.ucellf[2][2]*sAtomz;
#  endif
#endif
          // Bounding box test on SM 2.0+
          float trivialCut2 = 0.5625f * cutPlusSkin2;
          float bxc = __SHFL(WARP_MASK, shNlRecord.f, NEIGHBOR_CELLS + 6);
          float bxr = __SHFL(WARP_MASK, shNlRecord.f, NEIGHBOR_CELLS + 7);
          float byc = __SHFL(WARP_MASK, shNlRecord.f, NEIGHBOR_CELLS + 8);
          float byr = __SHFL(WARP_MASK, shNlRecord.f, NEIGHBOR_CELLS + 9);
          float bzc = __SHFL(WARP_MASK, shNlRecord.f, NEIGHBOR_CELLS + 10);
          float bzr = __SHFL(WARP_MASK, shNlRecord.f, NEIGHBOR_CELLS + 11);
          float tx = fabs(sAtomx - bxc);
          float ty = fabs(sAtomy - byc);
          float tz = fabs(sAtomz - bzc);
          tx = tx - min(tx, bxr);
          ty = ty - min(ty, byr);
          tz = tz - min(tz, bzr);
          float tr2 = tx*tx + ty*ty + tz*tz;
          PMEMask bpred = __BALLOT(WARP_MASK, (tgx < xmax) && (tr2 < cutPlusSkin2));
          PMEMask apred = __BALLOT(WARP_MASK, (tgx < xmax) && (tr2 < trivialCut2));

          // Bitwise exclusive OR plus assignment compound operator ^=
          // This is bpred = bpred ^ apred
          bpred ^= apred;

          // Perform tests on all non-trivial accepts in groups of (GRID / PME_ATOMS_PER_WARP)
          PMEMask mask = (PMEMask)cSim.NLAtomsPerWarpMask <<
                         (PME_ATOMS_PER_WARP * (tgx >> cSim.NLAtomsPerWarpBits));
          while (bpred) {
            int pos = maskFfs(bpred) - 1;
            bpred &= ~(PMEMask)1 << pos;
            // HIP-TODO: Support PME_ATOMS_PER_WARP = 32 and 8
#if (PME_ATOMS_PER_WARP == 16)
            int pos1 = maskFfs(bpred) - 1;
            if (pos1 != -1) {
              bpred &= ~(PMEMask)1 << pos1;
            }
            if (tgx >= 16) {
              pos = pos1;
            }
#ifdef AMBER_PLATFORM_AMD_WARP64
            pos1 = maskFfs(bpred) - 1;
            if (pos1 != -1) {
              bpred &= ~(PMEMask)1 << pos1;
            }
            if (tgx >= 32) {
              pos = pos1;
            }
            pos1 = maskFfs(bpred) - 1;
            if (pos1 != -1) {
              bpred &= ~(PMEMask)1 << pos1;
            }
            if (tgx >= 48) {
              pos = pos1;
            }
#endif
#elif (PME_ATOMS_PER_WARP == 8)
            int pos1 = maskFfs(bpred) - 1;
            if (pos1 != -1) {
              bpred &= ~(PMEMask)1 << pos1;
            }
            if (tgx >= 8) {
              pos = pos1;
            }
            pos1 = maskFfs(bpred) - 1;
            if (pos1 != -1) {
              bpred &= ~(PMEMask)1 << pos1;
            }
            if (tgx >= 16) {
              pos = pos1;
            }
            pos1 = maskFfs(bpred) - 1;
            if (pos1 != -1) {
              bpred &= ~(PMEMask)1 << pos1;
            }
            if (tgx >= 24) {
              pos = pos1;
            }
#endif

            float ax = __SHFL(WARP_MASK, sAtomx, pos);
            float ay = __SHFL(WARP_MASK, sAtomy, pos);
            float az = __SHFL(WARP_MASK, sAtomz, pos);
            int pred = 0;
            if (pos >= 0) {
              float dx = xi - ax;
              float dy = yi - ay;
              float dz = zi - az;
              float r2 = dx * dx + dy * dy + dz * dz;
              pred = (r2 < cutPlusSkin2);
            }

            // Signal acceptance or rejection of atoms
            if (__BALLOT(WARP_MASK, pred) & mask) {
              apred |= (PMEMask)1 << pos;
            }
          }
          __SYNCWARP(WARP_MASK);
          // HIP-TODO: Support PME_ATOMS_PER_WARP = 32 and 8
#if (PME_ATOMS_PER_WARP < 32)
          psWarp->exclusionMask[tgx] = apred;
          __SYNCWARP(WARP_MASK);
#if (PME_ATOMS_PER_WARP == 8)
          psWarp->exclusionMask[tgx] |= psWarp->exclusionMask[tgx ^ 8];
          __SYNCWARP(WARP_MASK);
          psWarp->exclusionMask[tgx] |= psWarp->exclusionMask[tgx ^ 16];
          __SYNCWARP(WARP_MASK);
          apred = psWarp->exclusionMask[tgx];
#else
#ifdef AMBER_PLATFORM_AMD_WARP64
          psWarp->exclusionMask[tgx] |= psWarp->exclusionMask[tgx ^ 16];
          __SYNCWARP(WARP_MASK);
          psWarp->exclusionMask[tgx] |= psWarp->exclusionMask[tgx ^ 32];
          __SYNCWARP(WARP_MASK);
          apred = psWarp->exclusionMask[tgx];
#else
          apred |= psWarp->exclusionMask[tgx ^ 16];
#endif
#endif
#endif
          bpred = apred;

          // Add all accepted atoms to atom list
          while (bpred) {
            int maxAccepts = min(GRID - atoms, maskPopc(bpred));
            int pos = -1;

            // Find number of predecessor bits and determine if thread can add atom
            if (bpred & threadmask) {
              PMEMask mask = threadmask - 1;
              pos = maskPopc(bpred & mask);
              if (pos >= maxAccepts) {
                pos = -1;
              }
            }

            // Accept each atom if there's room
            if (pos != -1) {
              unsigned int atom = xpos + tgx;
              minAtom = min(minAtom, atom);
              maxAtom = max(maxAtom, atom);
              psWarp->atomList[atoms + pos] = (atom << NLATOM_CELL_SHIFT) | shCellID;

              // Bitwise exclusive OR plus assignment compound operator
              bpred ^= threadmask;
            }
            __SYNCWARP(WARP_MASK);
            atoms += maxAccepts;

            // Output GRID atoms if ready
            if (atoms == GRID) {

              // Write swath of atoms to global memory
              cSim.pNLAtomList[psWarp->offset + tgx] = psWarp->atomList[tgx];
              __SYNCWARP(WARP_MASK);
              if (tgx == 0) {
                psWarp->offset += GRID;
                psWarp->nlEntry.NL.xatoms += GRID;
              }
              __SYNCWARP(WARP_MASK);

              // Clear used bits from bpred
              bpred &= __SHFL(WARP_MASK, bpred, tgx ^ 1);
              bpred &= __SHFL(WARP_MASK, bpred, tgx ^ 2);
              bpred &= __SHFL(WARP_MASK, bpred, tgx ^ 4);
              bpred &= __SHFL(WARP_MASK, bpred, tgx ^ 8);
              bpred &= __SHFL(WARP_MASK, bpred, tgx ^ 16);
#ifdef AMBER_PLATFORM_AMD_WARP64
              bpred &= __SHFL(WARP_MASK, bpred, tgx ^ 32);
#endif

              // Reduce minatom and maxatom
              minAtom = min(minAtom, __SHFL(WARP_MASK, minAtom, tgx ^ 1));
              minAtom = min(minAtom, __SHFL(WARP_MASK, minAtom, tgx ^ 2));
              minAtom = min(minAtom, __SHFL(WARP_MASK, minAtom, tgx ^ 4));
              minAtom = min(minAtom, __SHFL(WARP_MASK, minAtom, tgx ^ 8));
              minAtom = min(minAtom, __SHFL(WARP_MASK, minAtom, tgx ^ 16));
              maxAtom = max(maxAtom, __SHFL(WARP_MASK, maxAtom, tgx ^ 1));
              maxAtom = max(maxAtom, __SHFL(WARP_MASK, maxAtom, tgx ^ 2));
              maxAtom = max(maxAtom, __SHFL(WARP_MASK, maxAtom, tgx ^ 4));
              maxAtom = max(maxAtom, __SHFL(WARP_MASK, maxAtom, tgx ^ 8));
              maxAtom = max(maxAtom, __SHFL(WARP_MASK, maxAtom, tgx ^ 16));
#ifdef AMBER_PLATFORM_AMD_WARP64
              minAtom = min(minAtom, __SHFL(WARP_MASK, minAtom, tgx ^ 32));
              maxAtom = max(maxAtom, __SHFL(WARP_MASK, maxAtom, tgx ^ 32));
#endif

              // Search for y atom exclusions matching any x atom all at once (this should
              // reduce exclusion tests by a factor of approximately 100 overall).
              // But first, rule out skipping exclusion test.
              uint minExclusion = __SHFL(WARP_MASK, shNlRecord.u, NEIGHBOR_CELLS + 12);
              uint maxExclusion = __SHFL(WARP_MASK, shNlRecord.u, NEIGHBOR_CELLS + 13);

              psWarp->exclusionMask[tgx] = 0;
              __SYNCWARP(WARP_MASK);
              if ((minAtom <= maxExclusion) && (maxAtom >= minExclusion)) {
                uint atom = (psWarp->atomList[tgx] >> NLATOM_CELL_SHIFT);
                for (int j = 0; j < totalExclusions; j += GRID) {
                  int offset = j + tgx;
                  unsigned int eatom = 0xffffffff;
                  if (offset < totalExclusions) {
                    eatom = psExclusion[offset] >> NLEXCLUSION_SHIFT;
                  }
                  PMEMask vote = __BALLOT(WARP_MASK, (eatom >= minAtom) && (eatom <= maxAtom));

                  while (vote) {
                    unsigned int k = maskFfs(vote) - 1;
                    offset = k + j;
                    eatom = psExclusion[offset] >> NLEXCLUSION_SHIFT;
                    if (atom == eatom) {
                      psWarp->exclusionMask[psExclusion[offset] & NLEXCLUSION_ATOM_MASK] |=
                        threadmask;
                    }
                    __SYNCWARP(WARP_MASK);
                    vote ^= (PMEMask)1 << k;
                  }
                }
              }
              __SYNCWARP(WARP_MASK);

              // Output exclusion masks
              if (tgx < cSim.NLAtomsPerWarp) {
                PMEMask emask = psWarp->exclusionMask[tgx];
                ((PMEMask*)&cSim.pNLAtomList[psWarp->offset])[tgx] =
                  ((emask >> tgx) & exclusionMask) |
                  ((emask << (cSim.NLAtomsPerWarp - tgx)) & ~exclusionMask);
              }
              __SYNCWARP(WARP_MASK);
              if (tgx == 0)
                psWarp->offset += cSim.NLAtomsPerWarp * sizeof(PMEMask) / sizeof(unsigned int);
              atoms = 0;
              psWarp->atomList[tgx] = 0;
              minAtom = cSim.atoms;
              maxAtom = 0;
              __SYNCWARP(WARP_MASK);

              // Output neighbor list entry if width is sufficient
              if (psWarp->nlEntry.NL.xatoms >= cSim.NLXEntryWidth) {
                uint nlpos;
                if (tgx == 0) {
                  nlpos = atomicAdd(cSim.pNLEntries, 1);
                }
                nlpos = __SHFL(WARP_MASK, nlpos, 0);
                if (tgx < 4) {
                  cSim.pNLEntry[nlpos].array[tgx] = psWarp->nlEntry.array[tgx];
                }
                __SYNCWARP(0xFFFFFFFF);
                if (tgx == 0) {
                  psWarp->nlEntry.NL.xatoms  = 0;
                  psWarp->nlEntry.NL.ymax   &= ~NLENTRY_HOME_CELL_MASK;
                  psWarp->nlEntry.NL.offset  = psWarp->offset;
                }
              }
            }
            else {
              bpred = 0;
            }
            // End contingency for atoms == GRID
          }
          // End loop for adding atoms to the non-bonded list while bpred != 0

          // Move to next swath of atoms
          xpos += GRID;
        }

        // Move to next cell
        cpos++;
      }

      // Output last batch of atoms for this swath
      if (atoms > 0) {
        // Reduce minatom and maxatom
        minAtom = min(minAtom, __SHFL(WARP_MASK, minAtom, tgx ^ 1));
        minAtom = min(minAtom, __SHFL(WARP_MASK, minAtom, tgx ^ 2));
        minAtom = min(minAtom, __SHFL(WARP_MASK, minAtom, tgx ^ 4));
        minAtom = min(minAtom, __SHFL(WARP_MASK, minAtom, tgx ^ 8));
        minAtom = min(minAtom, __SHFL(WARP_MASK, minAtom, tgx ^ 16));
        maxAtom = max(maxAtom, __SHFL(WARP_MASK, maxAtom, tgx ^ 1));
        maxAtom = max(maxAtom, __SHFL(WARP_MASK, maxAtom, tgx ^ 2));
        maxAtom = max(maxAtom, __SHFL(WARP_MASK, maxAtom, tgx ^ 4));
        maxAtom = max(maxAtom, __SHFL(WARP_MASK, maxAtom, tgx ^ 8));
        maxAtom = max(maxAtom, __SHFL(WARP_MASK, maxAtom, tgx ^ 16));
#ifdef AMBER_PLATFORM_AMD_WARP64
        minAtom = min(minAtom, __SHFL(WARP_MASK, minAtom, tgx ^ 32));
        maxAtom = max(maxAtom, __SHFL(WARP_MASK, maxAtom, tgx ^ 32));
#endif
        if (tgx < atoms) {
          cSim.pNLAtomList[psWarp->offset + tgx] = psWarp->atomList[tgx];
        }
        else {
          cSim.pNLAtomList[psWarp->offset + tgx] = 0;
        }
        __SYNCWARP(WARP_MASK);
        if (tgx == 0)
          psWarp->offset += GRID;

        // Search for y atom exclusions matching any x atom all at once (this should
        // reduce exclusion tests by a factor of approximately 100 overall)
        uint minExclusion = __SHFL(WARP_MASK, shNlRecord.u, NEIGHBOR_CELLS + 12);
        uint maxExclusion = __SHFL(WARP_MASK, shNlRecord.u, NEIGHBOR_CELLS + 13);
        psWarp->exclusionMask[tgx] = 0;
        __SYNCWARP(WARP_MASK);
        if ((tgx < atoms) && ((minAtom <= maxExclusion) && (maxAtom >= minExclusion))) {
          unsigned int atom = (psWarp->atomList[tgx] >> NLATOM_CELL_SHIFT);
          for (int j = 0; j < totalExclusions; j++) {
            if ((psExclusion[j] >> NLEXCLUSION_SHIFT) == atom) {
              atomicOr((PMEMask*)(&psWarp->exclusionMask[psExclusion[j] &
                                  NLEXCLUSION_ATOM_MASK]), threadmask);
            }
          }
        }
        __SYNCWARP(WARP_MASK);

        // Output exclusion masks
        if (tgx < cSim.NLAtomsPerWarp) {
          PMEMask emask = psWarp->exclusionMask[tgx];
          ((PMEMask*)&cSim.pNLAtomList[psWarp->offset])[tgx] =
            ((emask >> tgx) & exclusionMask) |
            ((emask << (cSim.NLAtomsPerWarp - tgx)) & ~exclusionMask);
        }
        __SYNCWARP(WARP_MASK);
        if (tgx == 0) {
          psWarp->nlEntry.NL.xatoms += atoms;
          psWarp->offset += cSim.NLAtomsPerWarp * sizeof(PMEMask) / sizeof(unsigned int);
        }
        atoms  = 0;
        psWarp->atomList[tgx] = 0;
      }
      __SYNCWARP(0xFFFFFFFF);
      // End contingency for committing the final batch of atoms

      // Output final neighbor list entry if xatoms > 0
      if ((psWarp->nlEntry.NL.xatoms > 0) || (psWarp->nlEntry.NL.ymax & NLENTRY_HOME_CELL_MASK)) {
        uint nlpos;
        if (tgx == 0) {
          nlpos = atomicAdd(cSim.pNLEntries, 1);
        }
        nlpos = __SHFL(WARP_MASK, nlpos, 0);
        if (tgx < 4) {
          cSim.pNLEntry[nlpos].array[tgx] = psWarp->nlEntry.array[tgx];
        }
      }
      __SYNCWARP(WARP_MASK);

      // Advance to next swath of atoms
      ypos += cSim.NLYDivisor * cSim.NLAtomsPerWarp;
    }

    // Advance to next NLRecord entry
    if (tgx == (NEIGHBOR_CELLS + 2)) {
      shNlRecord.u = atomicAdd(&cSim.pFrcBlkCounters[0], 1);
    }
  }
#undef VOLATILE
}
