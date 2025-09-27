#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// This code is #include'd by kNeighborList.cu as the kMapSubImageExpansion kernel.  The
// critical #define is YSHIFT_PENCILS, for situations in which the box angle made by the XY and
// YZ planes is less than 7/12 pi and the hexagonal prism cells are staggered depending on the
// value of Z.
//---------------------------------------------------------------------------------------------
{
  __shared__ int batchStart, bpos;
  __shared__ int qqHeapSums[1024], ljHeapSums[1024], focusSums[512], batchSums[32];
  __shared__ int lim1[1024], lim4[1024], ntMoves[22];

  // Stage I: take in the grouped total cell populations and compute a prefix sum to get
  //          the position at which to begin writing this batch's warp instructions
  qqHeapSums[threadIdx.x] = 0;
  ljHeapSums[threadIdx.x] = 0;
  int nheaps;
  int nbatch;
  nheaps = (cSim.nExpandedPencils >> 8);
  if (threadIdx.x < nheaps) {
    qqHeapSums[threadIdx.x] = cSim.pHcmbQQPopHeapSums[threadIdx.x];
    ljHeapSums[threadIdx.x] = cSim.pHcmbLJPopHeapSums[threadIdx.x];
  }
  nbatch = (cSim.nExpandedPencils + nCells - 1) / nCells;
  nheaps = 0x1 << (32 - __clz(nheaps));
  if (threadIdx.x == 0) {
    bpos = blockIdx.x;
  }
  __syncthreads();
  int maxstride = nheaps / 2;
  int stride = 1;
  while (stride <= maxstride) {
    int index = 2*(threadIdx.x + 1)*stride - 1;
    if (index < nheaps) {
      qqHeapSums[index] += qqHeapSums[index - stride];
      ljHeapSums[index] += ljHeapSums[index - stride];
    }
    stride *= 2;
    __syncthreads();
  }
  stride = maxstride / 2;
  while (stride > 0) {
    int index = (threadIdx.x + 1)*stride*2 - 1;
    if ((index + stride) < nheaps) {
      qqHeapSums[index + stride] += qqHeapSums[index];
      ljHeapSums[index + stride] += ljHeapSums[index];
    }
    stride = stride / 2;
    __syncthreads();
  }

  // Stage II: loop over batches: each block does a batch, possibly more than one
  //           if the system is huge.
  while (bpos < 2*nbatch) {

    // Get the focus sums surrounding the cells of interest--to bulletproof against
    // batches that span two heaps of 256 cells, two heaps will be taken.  Batches
    // can be no larger than 32 cells, but that's still huge!
    int bposCurr = bpos;
    int bpos1 = bposCurr - (bposCurr >= nbatch)*nbatch;
    int fstart = (bpos1 * nCells) >> 8;
    int popidx = 256*fstart + threadIdx.x;
    if (threadIdx.x < 512) {
      focusSums[threadIdx.x] = 0;
      if (popidx < cSim.nExpandedPencils) {
        if (bposCurr < nbatch) {
          focusSums[threadIdx.x] = cSim.pHcmbQQPopulations[popidx];
        }
        else {
          focusSums[threadIdx.x] = cSim.pHcmbLJPopulations[popidx];
        }
      }
    }
    __syncthreads();

    // Increment the value of bpos.  This happens between these two thread block
    // synchronizations to avoid having additional __synthreads() calls.  The
    // working value of bpos has been snapshotted as bposCurr in registers.
    if (threadIdx.x == 0) {
      bpos = atomicAdd(&cSim.pFrcBlkCounters[2], 1);
    }

    int ncells = 512;
    maxstride = 256;
    stride = 1;
    while (stride <= maxstride) {
      int index = 2*(threadIdx.x + 1)*stride - 1;
      if (index < ncells) {
        focusSums[index] += focusSums[index - stride];
      }
      stride *= 2;
      __syncthreads();
    }
    stride = maxstride / 2;
    while (stride > 0) {
      int index = (threadIdx.x + 1)*stride*2 - 1;
      if ((index + stride) < ncells) {
        focusSums[index + stride] += focusSums[index];
      }
      stride = stride / 2;
      __syncthreads();
    }

    // Advance the focus sums for cells in this batch.  The array focusSums contains
    // the population prefix sum for hash cells in a range that covers the cells in
    // this batch.  That range, and only that range, will be updated to reflect the
    // position of the group in the much larger expanded image.
    if (threadIdx.x >= nCells && threadIdx.x < 32) {
      batchSums[threadIdx.x] = 0;
    }
    if (popidx >= bpos1*nCells && popidx < (bpos1 + 1)*nCells) {
      int baseline;
      if (bposCurr < nbatch) {
        baseline = (fstart == 0) ? 0 : qqHeapSums[fstart - 1];
      }
      else {
        baseline = (fstart == 0) ? 0 : ljHeapSums[fstart - 1];
      }
      int bsumpos = threadIdx.x - (bpos1*nCells - fstart*256);
      batchSums[bsumpos] = focusSums[threadIdx.x] + baseline;
      if (bsumpos == 0) {
        batchStart = (threadIdx.x == 0) ? baseline : baseline + focusSums[threadIdx.x - 1];
      }
    }
    __syncthreads();

    // This block will now write warp instructions about cells in the batchSums array.
    // The batchsums array shows where the individual threads' instructions start--the
    // cell index to which the instructions pertain is derived from bpos1 * nCells,
    // and the cell index then indicates where in the expanded array to place the cells.
    // The challenge for each thread, then, is to determine which cell its instruction
    // is targeted at, which is done by figuring out which two elements of batchSums
    // the thread's instruction index lies between.
    int cellcount = (bpos1 == nbatch - 1) ? cSim.nExpandedPencils - bpos1*nCells : nCells;
    int endline = batchSums[cellcount - 1];
    if (threadIdx.x == 0 && bpos1 == nbatch - 1) {
      cSim.pExpansionInsrCount[(bposCurr >= nbatch)] = endline;
    }
    int insrpos = batchStart + threadIdx.x;
    while (insrpos < endline) {

      // Find the cell index by binary search over the array of limits
      int pivot  = nCells / 2;
      int dpivot = pivot / 2 + (pivot == 1);
      int search = 1;
      while (search != 0) {
        int move = (pivot < nCells && insrpos >= batchSums[pivot]) -
                   (pivot > 0 && insrpos < batchSums[pivot - 1]);
        pivot += move * dpivot;
        search *= move;
        dpivot = (dpivot / 2) + (dpivot == 1);
      }
      int cellidx = bpos1*nCells + pivot;

      // This is an important value that gives the particle's index in the range
      // [ lowest index particle in the expanded cell (cellpos = 0),
      //   highest index in the expanded cell (cellpos = 2*(lim4 - lim1) + lim2 - lim3) )
      int cellpos = (pivot == 0) ? insrpos - batchStart : insrpos - batchSums[pivot - 1];

      // The cell indices and positions within those cells to which each thread
      // pertains are now known.  Get the source cell and its limits one last time.
      int ywidth = cSim.nypencils + (2 * cSim.ypadding);
      int zpos = cellidx / ywidth;
      int ypos = cellidx - zpos*ywidth;
      int ysrc = ypos - cSim.ypadding;
      int zsrc = zpos - cSim.zpadding;
      ysrc += (((ysrc < 0) - (ysrc >= cSim.nypencils)) * cSim.nypencils);
      zsrc += ((zsrc < 0) * cSim.nzpencils);
      int srcidx = (zsrc * cSim.nypencils) + ysrc;
      if (bposCurr < nbatch) {
#if (__CUDA_ARCH__ >= 350) || defined(AMBER_PLATFORM_AMD)
        int lim1 = __ldg(&cSim.pHcmbQQCellLimits[srcidx]);
        int lim2 = __ldg(&cSim.pHcmbQQCellLimits[srcidx + cSim.npencils]);
        int lim3 = __ldg(&cSim.pHcmbQQCellLimits[srcidx + (2 * cSim.npencils)]);
        int lim4 = __ldg(&cSim.pHcmbQQCellLimits[srcidx + (3 * cSim.npencils)]);
#else
        int lim1 = cSim.pHcmbQQCellLimits[srcidx];
        int lim2 = cSim.pHcmbQQCellLimits[srcidx + cSim.npencils];
        int lim3 = cSim.pHcmbQQCellLimits[srcidx + (2 * cSim.npencils)];
        int lim4 = cSim.pHcmbQQCellLimits[srcidx + (3 * cSim.npencils)];
#endif
        lim2 += (lim2 < 0) * (lim1 + 1);
        lim3 += (lim3 < 0) * (lim4 + 1);

        // Decide whether to translate in X, Y, or Z
        int lowerpad = lim4 - lim3;
        int primpop = lim4 - lim1;
        int xtrans = 1 - (cellpos < lowerpad) + (cellpos >= primpop + lowerpad);
        int ytrans = 1 - (ypos < cSim.ypadding) + (ypos >= cSim.nypencils + cSim.ypadding);
        int ztrans = 1 - (zpos < cSim.zpadding);

        // Construct the translation / replication instruction
        int xsrc = cellpos - lowerpad + ((cellpos < lowerpad) -
                                         (cellpos >= primpop + lowerpad))*primpop;
        ysrc += cSim.ypadding;
        zsrc += cSim.zpadding;
        srcidx = (zsrc * ywidth) + ysrc;
        int2 insr;
        insr.x = (srcidx * cSim.QQCellSpace) + cSim.QQCapPadding + xsrc;
        insr.y = (cellidx * cSim.QQCellSpace) + cSim.QQCapPadding + cellpos - lowerpad;
        int tinsr = (ztrans << 16) | (ytrans << 8) | xtrans;
        cSim.pQQExpansions[insrpos] = insr;
        cSim.pQQTranslations[insrpos] = tinsr;
      }
      else {
#if (__CUDA_ARCH__ >= 350) || defined(AMBER_PLATFORM_AMD)
        int lim1 = __ldg(&cSim.pHcmbLJCellLimits[srcidx]);
        int lim2 = __ldg(&cSim.pHcmbLJCellLimits[srcidx + cSim.npencils]);
        int lim3 = __ldg(&cSim.pHcmbLJCellLimits[srcidx + (2 * cSim.npencils)]);
        int lim4 = __ldg(&cSim.pHcmbLJCellLimits[srcidx + (3 * cSim.npencils)]);
#else
        int lim1 = cSim.pHcmbLJCellLimits[srcidx];
        int lim2 = cSim.pHcmbLJCellLimits[srcidx + cSim.npencils];
        int lim3 = cSim.pHcmbLJCellLimits[srcidx + (2 * cSim.npencils)];
        int lim4 = cSim.pHcmbLJCellLimits[srcidx + (3 * cSim.npencils)];
#endif
        lim2 += (lim2 < 0) * (lim1 + 1);
        lim3 += (lim3 < 0) * (lim4 + 1);

        // Decide whether to translate in X, Y, or Z
        int lowerpad = lim4 - lim3;
        int primpop = lim4 - lim1;
        int xtrans = 1 - (cellpos < lowerpad) + (cellpos >= primpop + lowerpad);
        int ytrans = 1 - (ypos < cSim.ypadding) + (ypos >= cSim.nypencils + cSim.ypadding);
        int ztrans = 1 - (zpos < cSim.zpadding);

        // Construct the translation / replication instruction
        int xsrc = cellpos - lowerpad + ((cellpos < lowerpad) -
                                         (cellpos >= primpop + lowerpad))*primpop;
        ysrc += cSim.ypadding;
        zsrc += cSim.zpadding;
        srcidx = (zsrc * ywidth) + ysrc;
        int2 insr;
        insr.x = (srcidx * cSim.LJCellSpace) + cSim.LJCapPadding + xsrc;
        insr.y = (cellidx * cSim.LJCellSpace) + cSim.LJCapPadding + cellpos - lowerpad;
        int tinsr = (ztrans << 16) | (ytrans << 8) | xtrans;
        cSim.pLJExpansions[insrpos] = insr;
        cSim.pLJTranslations[insrpos] = tinsr;
      }

      // Advance the block until all atoms in this batch are done
      insrpos += blockDim.x;
    }
  }
  __syncthreads();

  // Stage III: the blocks of this kernel will be devoted to a different purpose once
  //            their work making particle expansion instructions is complete.  Continue
  //            from the earlier loop, packing the work of determining Neutral Territory
  //            import indices right behind the work of making expansion instructions.
  int yspan    = (cSim.ydimFull + 31 - 2*cSim.ypadding) / (32 - 2*cSim.ypadding);
  int zspan    = (cSim.zdimFull + 31 - cSim.zpadding) / (32 - cSim.zpadding);
  int nNTzones = yspan * zspan;
  while (bpos < 2*nbatch + 2*nNTzones) {

    // Again, take a secondary counter, this time to indicate
    // which zone of the hash cell grid to make NT indexing for
    int bpos1     = bpos - 2*nbatch;
    int znZidx, znYidx;
    if (bpos1 < nNTzones) {
      znZidx    = bpos1 / yspan;
      znYidx    = bpos1 - znZidx*yspan;
    }
    else {
      znZidx    = (bpos1 - nNTzones) / yspan;
      znYidx    = bpos1 - nNTzones - znZidx*yspan;
    }
    int lYabs     = znYidx * (32 - 2*cSim.ypadding);
    int hYabs     = min(lYabs + 32, cSim.ydimFull);
    int lZabs     = znZidx * (32 - cSim.zpadding);
    int hZabs     = min(lZabs + 32, cSim.zdimFull);
    int zoneYdim  = hYabs - lYabs;
    int zoneZdim  = hZabs - lZabs;
    int cellAbsZ  = threadIdx.x / zoneYdim;
    int cellAbsY  = threadIdx.x - zoneYdim*cellAbsZ;
    cellAbsY     += lYabs;
    cellAbsZ     += lZabs;

    // By construction, cellAbsY is within bounds, so only check cellAbsZ.  If it's
    // within the grid, assign qqHeapSums as if it were the lowest cap bound and
    // ljHeapSums as if it were the highest cap bound.
    if (cellAbsZ < cSim.zdimFull) {
      int cellSrcY  = cellAbsY - cSim.ypadding +
                      (((cellAbsY < cSim.ypadding) -
                        (cellAbsY >= cSim.ypadding + cSim.nypencils)) * cSim.nypencils);
      int cellSrcZ  = cellAbsZ - cSim.zpadding + ((cellAbsZ < cSim.zpadding) * cSim.nzpencils);
      int srcidx = (cellSrcZ * cSim.nypencils) + cellSrcY;
      if (bpos1 < nNTzones) {
	int cellidx  = (cellAbsZ * cSim.ydimFull) + cellAbsY;
        int tmplim1  = cSim.pHcmbQQCellLimits[srcidx];
        int tmplim2  = cSim.pHcmbQQCellLimits[srcidx + cSim.npencils];
        int tmplim3  = cSim.pHcmbQQCellLimits[srcidx + (2 * cSim.npencils)];
        int tmplim4  = cSim.pHcmbQQCellLimits[srcidx + (3 * cSim.npencils)];
        tmplim2     += (tmplim2 < 0) * (tmplim1 + 1);
        tmplim3     += (tmplim3 < 0) * (tmplim4 + 1);
	int plim1    = (cellidx * cSim.QQCellSpace) + cSim.QQCapPadding;
        lim1[threadIdx.x] =       plim1;
        lim4[threadIdx.x] =       plim1 + tmplim4 - tmplim1;
	qqHeapSums[threadIdx.x] = plim1 - tmplim4 + tmplim3;
	ljHeapSums[threadIdx.x] = plim1 + tmplim4 - tmplim1 + tmplim2 - tmplim1;
      }
      else {
	int cellidx  = (cellAbsZ * cSim.ydimFull) + cellAbsY;
        int tmplim1  = cSim.pHcmbLJCellLimits[srcidx];
        int tmplim2  = cSim.pHcmbLJCellLimits[srcidx + cSim.npencils];
        int tmplim3  = cSim.pHcmbLJCellLimits[srcidx + (2 * cSim.npencils)];
        int tmplim4  = cSim.pHcmbLJCellLimits[srcidx + (3 * cSim.npencils)];
        tmplim2     += (tmplim2 < 0) * (tmplim1 + 1);
        tmplim3     += (tmplim3 < 0) * (tmplim4 + 1);
	int plim1    = (cellidx * cSim.LJCellSpace) + cSim.LJCapPadding;
        lim1[threadIdx.x] =       plim1;
        lim4[threadIdx.x] =       plim1 + tmplim4 - tmplim1;
	qqHeapSums[threadIdx.x] = plim1 - tmplim4 + tmplim3;
	ljHeapSums[threadIdx.x] = plim1 + tmplim4 - tmplim1 + tmplim2 - tmplim1;
      }
    }

    // Last two warps read the table of Neutral Territory pencils (these warps
    // are most likely to have been idle during the previous reading cycle)
    if (threadIdx.x >= blockDim.x - 22) {
      int ntpos = threadIdx.x - blockDim.x + 22;
      ntMoves[ntpos] = cSim.pNTImports[ntpos];
    }
    __syncthreads();

    // Again, increment the value of bpos in between thread block synchronizations.
    // The critical information has already been snapshotted as bpos1 in registers.
    if (threadIdx.x == 0) {
      bpos = atomicAdd(&cSim.pFrcBlkCounters[2], 1);
    }

    // The table is made.  Now, each thread that controls a cell in the primary image
    // can write out a Neutral Territory decomposition.
    int cellidx = (threadIdx.x >> GRID_BITS);
    int tgx = (threadIdx.x & GRID_BITS_MASK);
    while (cellidx < zoneYdim * zoneZdim) {
      int znsrcz  = cellidx / zoneYdim;
      int znsrcy  = cellidx - zoneYdim*znsrcz;
      if (znsrcy >= cSim.ypadding && znsrcy < zoneYdim - cSim.ypadding &&
	  znsrcz >= cSim.zpadding) {
        cellAbsY     = znsrcy + lYabs;
        cellAbsZ     = znsrcz + lZabs;
        int cellSrcY  = cellAbsY - cSim.ypadding +
                        (((cellAbsY < cSim.ypadding) -
                          (cellAbsY >= cSim.ypadding + cSim.nypencils)) * cSim.nypencils);
        int cellSrcZ  = cellAbsZ - cSim.zpadding +
                        ((cellAbsZ < cSim.zpadding) * cSim.nzpencils);
        if (tgx < 11) {
          int mvY   = ntMoves[2*tgx];
	  int mvZ   = ntMoves[2*tgx + 1];
#ifdef YSHIFT_PENCILS
	  mvY      += ((cellSrcZ * -mvZ) & 0x1);
#endif
          int impY  = znsrcy + mvY;
          int impZ  = znsrcz + mvZ;
          int znidx = impZ*zoneYdim + impY;
	  int4 limits;
          limits.x = qqHeapSums[znidx];
	  limits.y = lim1[znidx];
	  limits.z = lim4[znidx];
	  limits.w = ljHeapSums[znidx];
	  int srcidx = (cellSrcZ * cSim.nypencils) + cellSrcY;
	  if (bpos1 < nNTzones) {
            cSim.pQQStencilLimits[11*srcidx + tgx] = limits;
          }
          else {
	    cSim.pLJStencilLimits[11*srcidx + tgx] = limits;
	  }
        }
      }
      cellidx += (blockDim.x >> GRID_BITS);
    }
    __syncthreads();
  }
}
