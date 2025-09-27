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
  struct BNLWarp {
    uint2 exclusionStartCount[PME_ATOMS_PER_WARP];
    uint atomList[GRID];
    uint homeCellExclusionMask[PME_ATOMS_PER_WARP];
    uint setBitPos[GRID];
    bool accepted[GRID];
    NLEntry nlEntry;

    // Bounding box
    float bxc;
    float bxr;
    float byc;
    float byr;
    float bzc;
    float bzr;
  };

#if (PME_ATOMS_PER_WARP == 32)
  const int THREADS_PER_BLOCK = NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK;
#elif (PME_ATOMS_PER_WARP == 16)
  const int THREADS_PER_BLOCK = NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK;
#else
  const int THREADS_PER_BLOCK = NLBUILD_NEIGHBORLIST8_THREADS_PER_BLOCK;
#endif
  __shared__ BNLWarp sNLWarp[THREADS_PER_BLOCK / GRID];
#ifdef PME_VIRIAL
  __shared__ float sUcellf[9];
#  ifdef PME_IS_ORTHOGONAL
  __shared__ float3 sCellOffset[NEIGHBOR_CELLS];
#  endif
#endif
#ifdef PME_VIRIAL
  if (threadIdx.x < 9) {
    sUcellf[threadIdx.x] = cSim.pNTPData->ucellf[threadIdx.x];
  }
#  ifdef PME_IS_ORTHOGONAL
  __syncthreads();
  if (threadIdx.x < NEIGHBOR_CELLS) {
    sCellOffset[threadIdx.x].x = sUcellf[0] * cSim.cellOffset[threadIdx.x][0];
    sCellOffset[threadIdx.x].y = sUcellf[4] * cSim.cellOffset[threadIdx.x][1];
    sCellOffset[threadIdx.x].z = sUcellf[8] * cSim.cellOffset[threadIdx.x][2];
  }
#  endif
  float cutPlusSkin2 = cSim.pNTPData->cutPlusSkin2;
#else
  float cutPlusSkin2 = cSim.cutPlusSkin2;
#endif

  unsigned int tgx = threadIdx.x & GRID_BITS_MASK;
  unsigned int warp = THREADS_PER_BLOCK == GRID ? 0 : (threadIdx.x >> GRID_BITS);
  unsigned int globalWarp = blockIdx.x * (THREADS_PER_BLOCK >> GRID_BITS) + warp;
  BNLWarp* psWarp = &sNLWarp[warp];

  uint* psExclusion = &cSim.pBNLExclusionBuffer[globalWarp * cSim.NLExclusionBufferSize];

  PMEMask exclusionMask =
    cSim.NLAtomsPerWarpMask >> (tgx & cSim.NLAtomsPerWarpBitsMask);
  for (int shift = PME_ATOMS_PER_WARP; shift < GRID; shift <<= 1) {
    exclusionMask = exclusionMask | (exclusionMask << shift);
  }
  PMEMask threadmask = (PMEMask)1 << tgx;

  // This is the beginning of a loop that extends to the end of the library.
  while (globalWarp < cSim.NLRecords) {
    uint2 shCell;
    uint shNeighborCell;

    // Calculate Exclusion/neighbor list space required
    uint neighborCells = cSim.pNLRecord[globalWarp].NL.neighborCells;
    int atomOffset = neighborCells >> NLRECORD_YOFFSET_SHIFT;
    uint2 homeCell = cSim.pNLNonbondCellStartEnd[cSim.pNLRecord[globalWarp].NL.homeCell];

    int ysize = max(0, (int)(homeCell.y - homeCell.x - atomOffset * cSim.NLAtomsPerWarp));
    if (ysize > 0) {
      ysize = 1 + max(0, ysize - 1) / (cSim.NLAtomsPerWarp * cSim.NLYDivisor);
    }

    // Load cell information and calculate maximum required space
    unsigned int cells = neighborCells & NLRECORD_CELL_COUNT_MASK;
    uint xsize = 0;
    if (tgx < cells) {
      shNeighborCell = cSim.pNLRecord[globalWarp].NL.neighborCell[tgx];
      shCell = cSim.pNLNonbondCellStartEnd[shNeighborCell >> NLRECORD_CELL_SHIFT];
      xsize = shCell.y - shCell.x;
    }

    // Reduce xsize
    for (int shift = 1; shift <= (NEIGHBOR_CELLS > 16 ? 16 : 8); shift <<= 1) {
      xsize += __SHFL(WARP_MASK, xsize, tgx ^ shift);
    }

    uint offset = 0;
    if (tgx == 0) {
      uint totalXSize = ((xsize + GRID - 1) >> GRID_BITS);
      offset = atomicAdd(cSim.pNLTotalOffset,
                         totalXSize*ysize*cSim.NLOffsetPerWarp + cSim.NLAtomsPerWarp);
    }
    offset = __SHFL(WARP_MASK, offset, 0);

    // Generate actual neighbor list entry
    uint ypos = homeCell.x + atomOffset * cSim.NLAtomsPerWarp;
    uint homecellY = homeCell.y;
    while (ypos < homecellY) {
      // Calculate y bounds and set to calculate homecell interaction
      uint ymax = min(ypos + cSim.NLAtomsPerWarp, homecellY);
      if (tgx == 0) {
        psWarp->nlEntry.NL.ypos = ypos;
        psWarp->nlEntry.NL.ymax = (ymax << NLENTRY_YMAX_SHIFT) | NLENTRY_HOME_CELL_MASK;
        psWarp->nlEntry.NL.xatoms = 0;
        psWarp->nlEntry.NL.offset = offset;
      }
      __SYNCWARP(WARP_MASK);

      // Read y atoms
      float xi;
      float yi;
      float zi;
      unsigned int index = ypos + (tgx & cSim.NLAtomsPerWarpBitsMask);
      if (index < ymax) {
#ifndef use_DPFP
        PMEFloat2 xy = tex1Dfetch<float2>(cSim.texAtomXYSP, index);
        zi = tex1Dfetch<float>(cSim.texAtomZSP, index);
#else
        PMEFloat2 xy = cSim.pAtomXYSP[index];
        zi = cSim.pAtomZSP[index];
#endif
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
#  ifdef PME_VIRIAL
      xi = sUcellf[0]*xi + sUcellf[1]*yi + sUcellf[2]*zi;
      yi = sUcellf[4]*yi + sUcellf[5]*zi;
      zi = sUcellf[8]*zi;
#  else
      xi = cSim.ucellf[0][0]*xi + cSim.ucellf[0][1]*yi + cSim.ucellf[0][2]*zi;
      yi = cSim.ucellf[1][1]*yi + cSim.ucellf[1][2]*zi;
      zi = cSim.ucellf[2][2]*zi;
#  endif
#endif

      // Calculate bounding box
      float bmin = (index < ymax) ? 0.5f * xi :  999999.0f;
      float bmax = (index < ymax) ? 0.5f * xi : -999999.0f;
      for (int shift = 1; shift < PME_ATOMS_PER_WARP; shift <<= 1) {
        bmin = min(__SHFL(WARP_MASK, bmin, tgx ^ shift), bmin);
        bmax = max(__SHFL(WARP_MASK, bmax, tgx ^ shift), bmax);
      }
      psWarp->bxc = bmax + bmin;
      psWarp->bxr = bmax - bmin;

      bmin = (index < ymax) ? 0.5f * yi :  999999.0f;
      bmax = (index < ymax) ? 0.5f * yi : -999999.0f;
      for (int shift = 1; shift < PME_ATOMS_PER_WARP; shift <<= 1) {
        bmin = min(__SHFL(WARP_MASK, bmin, tgx ^ shift), bmin);
        bmax = max(__SHFL(WARP_MASK, bmax, tgx ^ shift), bmax);
      }
      psWarp->byc = bmax + bmin;
      psWarp->byr = bmax - bmin;

      bmin = (index < ymax) ? 0.5f * zi :  999999.0f;
      bmax = (index < ymax) ? 0.5f * zi : -999999.0f;
      for (int shift = 1; shift < PME_ATOMS_PER_WARP; shift <<= 1) {
        bmin = min(__SHFL(WARP_MASK, bmin, tgx ^ shift), bmin);
        bmax = max(__SHFL(WARP_MASK, bmax, tgx ^ shift), bmax);
      }
      psWarp->bzc = bmax + bmin;
      psWarp->bzr = bmax - bmin;

      // Read exclusions
      uint limit = ymax - ypos;
      if (tgx < limit) {
        uint atom = cSim.pImageAtom[index];
        psWarp->exclusionStartCount[tgx] = cSim.pNLExclusionStartCount[atom];
      }
      __SYNCWARP(WARP_MASK);

      // Load all exclusions as-is (original atom indices)
      uint unfilteredExclusions = 0;
      for (int i = 0; i < limit; i++) {
        uint2 exclusionStartCount = psWarp->exclusionStartCount[i];
        uint start = exclusionStartCount.x;
        uint count = exclusionStartCount.y;
        for (int j = tgx; j < count; j += GRID) {
          uint atom = cSim.pNLExclusionList[start + j];
          psExclusion[unfilteredExclusions + j] = (atom << NLEXCLUSION_SHIFT) | i;
        }
        unfilteredExclusions += count;
      }

      if (tgx < cSim.NLAtomsPerWarp) {
        psWarp->homeCellExclusionMask[tgx] = 0;
      }
      __SYNCWARP(WARP_MASK);

      // Load atoms using their original indices and split exclusions in two parts:
      // 1) exclusions of current y atoms (for homecell interactions) -> build exclusion
      // mask immediately;
      // 2) other exclusions -> store them for future use into the same memory and count them.
      uint totalExclusions = 0;
      uint minExclusion = cSim.atoms;
      uint maxExclusion = 0;
      for (int j = 0; j < unfilteredExclusions; j += GRID) {
        uint atom = UINT_MAX;
        uint i = 0;
        bool inBound = j + tgx < unfilteredExclusions;
        if (inBound) {
          uint excl = psExclusion[j + tgx];
          i = (excl & NLEXCLUSION_ATOM_MASK);
          atom = cSim.pImageAtomLookup[excl >> NLEXCLUSION_SHIFT];
        }
        __SYNCWARP(WARP_MASK);
        bool selfExclusion = (atom >= ypos) && (atom < ymax);
        PMEMask m = __BALLOT(WARP_MASK, inBound && !selfExclusion);
        if (inBound) {
          if (selfExclusion) {
            // Calculate exclusions assuming all atoms are in range of each other
            // Only PME_ATOMS_PER_WARP bits are used, this means that we can use uint
            // instead of PMEMask for homeCellExclusionMask
            uint pos = atom - ypos;
            atomicOr(&psWarp->homeCellExclusionMask[pos], 1 << i);
          }
          else {
            minExclusion = min(minExclusion, atom);
            maxExclusion = max(maxExclusion, atom);
            // Write filtered exclusions
            uint pos = maskPopc(m & (threadmask - 1));
            psExclusion[totalExclusions + pos] = (atom << NLEXCLUSION_SHIFT) | i;
          }
        }
        totalExclusions += maskPopc(m);
      }
      for (int shift = 1; shift < GRID; shift <<= 1) {
        minExclusion = min(__SHFL(WARP_MASK, minExclusion, tgx ^ shift), minExclusion);
        maxExclusion = max(__SHFL(WARP_MASK, maxExclusion, tgx ^ shift), maxExclusion);
      }

      // Initialize Neighbor List variables for current line of entry
      unsigned int cpos = 0;
      unsigned int atoms = 0;

      while (cpos < cells) {

        // Check for home cell
        uint shCellID = __SHFL(WARP_MASK, shNeighborCell, cpos) & NLRECORD_CELL_TYPE_MASK;
        uint xpos;

        // Cell 0 always starts along force matrix diagonal
        if ((cpos == 0) && (psWarp->nlEntry.NL.ymax & NLENTRY_HOME_CELL_MASK)) {
          // Output exclusion masks
          if (tgx < cSim.NLAtomsPerWarp) {
            uint mask = psWarp->homeCellExclusionMask[tgx];
            mask = ((mask >> (1 + tgx)) | (mask << (cSim.NLAtomsPerWarp - tgx - 1))) &
                   cSim.NLAtomsPerWarpMask;
            cSim.pNLAtomList[offset + tgx] = mask;
          }
          offset += cSim.NLAtomsPerWarp;
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
#ifndef use_DPFP
            PMEFloat2 xy = tex1Dfetch<float2>(cSim.texAtomXYSP, xpos + tgx);
            sAtomz = tex1Dfetch<float>(cSim.texAtomZSP, xpos + tgx);
#else
            PMEFloat2 xy = cSim.pAtomXYSP[xpos + tgx];
            sAtomz = cSim.pAtomZSP[xpos + tgx];
#endif
            sAtomx = xy.x;
            sAtomy = xy.y;
          }
          else {
            sAtomx = INFINITY;
            sAtomy = INFINITY;
            sAtomz = INFINITY;
          }

          // Translate all atoms into a local coordinate system within one unit
          // cell of the first atom read to avoid PBC handling within inner loops
#if defined(PME_VIRIAL) && defined(PME_IS_ORTHOGONAL)
          sAtomx += sCellOffset[shCellID].x;
          sAtomy += sCellOffset[shCellID].y;
          sAtomz += sCellOffset[shCellID].z;
#else
          sAtomx += cSim.cellOffset[shCellID][0];
          sAtomy += cSim.cellOffset[shCellID][1];
          sAtomz += cSim.cellOffset[shCellID][2];
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
          // Bounding box test
          float trivialCut2 = 0.5625f * cutPlusSkin2;
          float tx = fabs(sAtomx - psWarp->bxc);
          float ty = fabs(sAtomy - psWarp->byc);
          float tz = fabs(sAtomz - psWarp->bzc);
          tx = tx - min(tx, psWarp->bxr);
          ty = ty - min(ty, psWarp->byr);
          tz = tz - min(tz, psWarp->bzr);
          float tr2 = tx*tx + ty*ty + tz*tz;
          bool between = (tr2 < cutPlusSkin2) && !(tr2 < trivialCut2);
          PMEMask bpred = __BALLOT(WARP_MASK, between);

          // Perform tests on all non-trivial accepts in groups of (GRID / PME_ATOMS_PER_WARP)
          int groupI = tgx >> cSim.NLAtomsPerWarpBits;
          // Find indices of all non-zero bits of bpred
          int setBitCount = maskPopc(bpred);
          if (between) {
            int pos = maskPopc(bpred & (threadmask - 1));
            psWarp->setBitPos[pos] = tgx;
          }
          psWarp->accepted[tgx] = (tr2 < trivialCut2);
          __SYNCWARP(WARP_MASK);
          for (int groupOffset = 0;
               groupOffset < setBitCount;
               groupOffset += GRID / PME_ATOMS_PER_WARP) {
            int pos = psWarp->setBitPos[groupOffset + groupI];

            float ax = __SHFL(WARP_MASK, sAtomx, pos);
            float ay = __SHFL(WARP_MASK, sAtomy, pos);
            float az = __SHFL(WARP_MASK, sAtomz, pos);
            float dx = xi - ax;
            float dy = yi - ay;
            float dz = zi - az;
            float r2 = dx * dx + dy * dy + dz * dz;

            // Signal acceptance or rejection of atoms
            if ((groupOffset + groupI < setBitCount) && (r2 < cutPlusSkin2)) {
              psWarp->accepted[pos] = true;
            }
          }
          __SYNCWARP(WARP_MASK);
          bpred = __BALLOT(WARP_MASK, psWarp->accepted[tgx]);

          // Add all accepted atoms to atom list
          while (bpred) {
            int maxAccepts = min(GRID - atoms, maskPopc(bpred));

            bool accepted = false;
            // Find number of predecessor bits and determine if thread can add atom
            if (bpred & threadmask) {
              int pos = maskPopc(bpred & (threadmask - 1));
              // Accept each atom if there's room
              if (pos < maxAccepts) {
                unsigned int atom = xpos + tgx;
                psWarp->atomList[atoms + pos] = (atom << NLATOM_CELL_SHIFT) | shCellID;
                accepted = true;
              }
            }
            __SYNCWARP(WARP_MASK);
            atoms += maxAccepts;

            // Clear used bits from bpred
            bpred = bpred ^ __BALLOT(WARP_MASK, accepted);

            // Output GRID atoms if ready
            if (atoms == GRID) {
              uint atomCell = psWarp->atomList[tgx];
              uint atom = (atomCell >> NLATOM_CELL_SHIFT);

              // Write swath of atoms to global memory
              cSim.pNLAtomList[offset + tgx] = atomCell;
              offset += GRID;
              if (tgx == 0) {
                psWarp->nlEntry.NL.xatoms += GRID;
              }

              // Search for y atom exclusions matching any x atom all at once (this should
              // reduce exclusion tests by a factor of approximately 100 overall).
              // But first, rule out skipping exclusion test.
              PMEMask emask = 0;
              if (__ANY(WARP_MASK, (atom <= maxExclusion) && (atom >= minExclusion))) {
                // Reduce minAtom and maxAtom
                uint minAtom = atom;
                uint maxAtom = atom;
                for (int shift = 1; shift < GRID; shift <<= 1) {
                  minAtom = min(minAtom, __SHFL(WARP_MASK, minAtom, tgx ^ shift));
                  maxAtom = max(maxAtom, __SHFL(WARP_MASK, maxAtom, tgx ^ shift));
                }
                for (int j = 0; j < totalExclusions; j += GRID) {
                  unsigned int excl = UINT_MAX;
                  unsigned int eatom = UINT_MAX;
                  if (j + tgx < totalExclusions) {
                    excl = psExclusion[j + tgx];
                    eatom = excl >> NLEXCLUSION_SHIFT;
                  }
                  PMEMask vote = __BALLOT(WARP_MASK, (eatom >= minAtom) && (eatom <= maxAtom));

                  while (vote) {
                    unsigned int k = maskFfs(vote) - 1;
                    unsigned int exclK = __SHFL(WARP_MASK, excl, k);
                    PMEMask m = __BALLOT(WARP_MASK, atom == (exclK >> NLEXCLUSION_SHIFT));
                    if ((exclK & NLEXCLUSION_ATOM_MASK) == tgx) {
                      emask |= m;
                    }
                    vote ^= (PMEMask)1 << k;
                  }
                }
              }

              // Output exclusion masks
              if (tgx < cSim.NLAtomsPerWarp) {
                ((PMEMask*)&cSim.pNLAtomList[offset])[tgx] =
                  ((emask >> tgx) & exclusionMask) |
                  ((emask << (cSim.NLAtomsPerWarp - tgx)) & ~exclusionMask);
              }
              offset += cSim.NLAtomsPerWarp * sizeof(PMEMask) / sizeof(unsigned int);
              atoms = 0;
              __SYNCWARP(WARP_MASK);

              // Output neighbor list entry if width is sufficient
              if (tgx == 0 && psWarp->nlEntry.NL.xatoms >= cSim.NLXEntryWidth) {
                uint nlpos = atomicAdd(cSim.pNLEntries, 1);
                cSim.pNLEntry[nlpos].NL = psWarp->nlEntry.NL;
                psWarp->nlEntry.NL.xatoms  = 0;
                psWarp->nlEntry.NL.ymax   &= ~NLENTRY_HOME_CELL_MASK;
                psWarp->nlEntry.NL.offset  = offset;
              }
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
        uint atomCell = psWarp->atomList[tgx];
        uint atom;
        uint minAtom;
        uint maxAtom;
        if (tgx < atoms) {
          cSim.pNLAtomList[offset + tgx] = atomCell;
          atom = (atomCell >> NLATOM_CELL_SHIFT);
          minAtom = atom;
          maxAtom = atom;
        }
        else {
          cSim.pNLAtomList[offset + tgx] = 0;
          atom = UINT_MAX;
          minAtom = cSim.atoms;
          maxAtom = 0;
        }
        offset += GRID;

        // Search for y atom exclusions matching any x atom all at once (this should
        // reduce exclusion tests by a factor of approximately 100 overall)
        PMEMask emask = 0;
        if (__ANY(WARP_MASK, (minAtom <= maxExclusion) && (maxAtom >= minExclusion))) {
          // Reduce minAtom and maxAtom
          for (int shift = 1; shift < GRID; shift <<= 1) {
            minAtom = min(minAtom, __SHFL(WARP_MASK, minAtom, tgx ^ shift));
            maxAtom = max(maxAtom, __SHFL(WARP_MASK, maxAtom, tgx ^ shift));
          }
          for (int j = 0; j < totalExclusions; j += GRID) {
            unsigned int excl = UINT_MAX;
            unsigned int eatom = UINT_MAX;
            if (j + tgx < totalExclusions) {
              excl = psExclusion[j + tgx];
              eatom = excl >> NLEXCLUSION_SHIFT;
            }
            PMEMask vote = __BALLOT(WARP_MASK, (eatom >= minAtom) && (eatom <= maxAtom));

            while (vote) {
              unsigned int k = maskFfs(vote) - 1;
              unsigned int exclK = __SHFL(WARP_MASK, excl, k);
              PMEMask m = __BALLOT(WARP_MASK, atom == (exclK >> NLEXCLUSION_SHIFT));
              if ((exclK & NLEXCLUSION_ATOM_MASK) == tgx) {
                emask |= m;
              }
              vote ^= (PMEMask)1 << k;
            }
          }
        }

        // Output exclusion masks
        if (tgx < cSim.NLAtomsPerWarp) {
          ((PMEMask*)&cSim.pNLAtomList[offset])[tgx] =
            ((emask >> tgx) & exclusionMask) |
            ((emask << (cSim.NLAtomsPerWarp - tgx)) & ~exclusionMask);
        }
        offset += cSim.NLAtomsPerWarp * sizeof(PMEMask) / sizeof(unsigned int);
        if (tgx == 0) {
          psWarp->nlEntry.NL.xatoms += atoms;
        }
        atoms  = 0;
      }
      // End contingency for committing the final batch of atoms

      // Output final neighbor list entry if xatoms > 0
      if (tgx == 0 &&
          ((psWarp->nlEntry.NL.xatoms > 0) ||
           (psWarp->nlEntry.NL.ymax & NLENTRY_HOME_CELL_MASK))) {
        uint nlpos = atomicAdd(cSim.pNLEntries, 1);
        cSim.pNLEntry[nlpos].NL = psWarp->nlEntry.NL;
      }
      __SYNCWARP(WARP_MASK);

      // Advance to next swath of atoms
      ypos += cSim.NLYDivisor * cSim.NLAtomsPerWarp;
    }

    // Advance to next NLRecord entry
    if (tgx == 0) {
      globalWarp = atomicAdd(&cSim.pFrcBlkCounters[0], 1);
    }
    globalWarp = __SHFL(WARP_MASK, globalWarp, 0);
  }
}
