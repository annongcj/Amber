#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------
#ifdef use_DPFP
#  define BATCHSIZE        32
#  define BATCH_BITS        5
#  define BATCH_BITS_MASK  31
#else
#  define BATCHSIZE        32
#  define BATCH_BITS        5
#  define BATCH_BITS_MASK  31
#endif
{
  // The PME_VIRIAL define doesn't actually specify computation of virials--rather, it
  // specifies that forces be scaled by the real space box transformation matrix as found
  // in cSim.pNTPData.  If anything other than a Berendsen barostat and NTP is in place,
  // the constants in cSim.recipf will do.
#ifdef PME_VIRIAL
  __shared__ PMEFloat sRecipf[9];
  if (threadIdx.x < 9) {
    sRecipf[threadIdx.x] = cSim.pNTPData->recipf[threadIdx.x];
  }
  __syncthreads();
#endif

  // Partial charges of all atoms
  __shared__ volatile PMEFloat q[BATCHSIZE];

  // Forces on all atoms, computer by each warp and added atomically
  __shared__ volatile PMEAccumulator frc[3*BATCHSIZE];

  // The starting mesh indices (minimum x, y, and z stencil indices) are stored in the
  // order [ ix1, ix2, ix3, ..., ix(BATCHSIZE), iy1, iy2, ..., iy(BATCHSIZE), iz1, ... ],
  // mirroring crdq (also follows from kCQB.h).
  __shared__ volatile int sAtom_ixyz[3*BATCHSIZE];

  // The B-spline cofficients arrays are stored in padded arrays to make sure that
  // the 12 (or, 18, for 6th order interpolation) coefficients all fall on different
  // __shared__ memory banks.  For 4th order interpolation, the order is as follows:
  //
  // [ (atom 1, point 1x), (atom 2, point 1x), ... (atom BATCHSIZE, point 1x),
  //   (atom 1, point 2x), ... (atom BATCHSIZE, point 2x), (atom 1, point 3x), ...
  //   (atom BATCHSIZE, point 3x), (atom 1, point4x), ... (atom BATCHSIZE, point 4x), ...
  //   (atom1, point 1y), (atom 2, point 1y), ... (atom BATCHSIZE, point 4y), ...
  //   (atom1, point 1z), (atom 2, point 1z), ... (atom BATCHSIZE, point 4z) ]
  //
  // This ordering is vastly different from the arrays in kCQB.h, to facilitate vectorized
  // access to the array elements by each warp.
#if (PME_ORDER == 4)
  __shared__ volatile PMEFloat sAtom_txyz[12*BATCHSIZE];
  __shared__ volatile PMEFloat sAtom_dtxyz[12*BATCHSIZE];
#elif (PME_ORDER == 6)
  __shared__ volatile PMEFloat sAtom_txyz[18*BATCHSIZE];
  __shared__ volatile PMEFloat sAtom_dtxyz[18*BATCHSIZE];
#endif

  // There are multiple warps.  Prepare here for the final interpolation stage,
  // when each warp will read from the mesh and accumulate portions of forces on
  // one group of 32 atoms.  Four warps will operate on each group.
  const int warpIdx = threadIdx.x >> GRID_BITS;
#if (PME_ORDER == 4)
  const int zOffset = (warpIdx & 3);
#elif (PME_ORDER == 6)
  const int zOffset = (warpIdx & 3) + ((warpIdx & 3) >= 2);
  const int yOffset = 3*(warpIdx & 1);
#endif

  // Threads will initially cooperate to read atom coordinates and charges, then compute
  // B-spline coefficients in three dimensions.  Later, warps will cooperate to map each
  // atom to the mesh.  The variable ndim is only valid for threads 0-767.
  const int ndim = threadIdx.x / BATCHSIZE;
  const int bchidx = (threadIdx.x & BATCH_BITS_MASK);
  const int nfftval = (ndim == 0)*cSim.nfft1 + (ndim == 1)*cSim.nfft2 + (ndim == 2)*cSim.nfft3;

  // Constants for B-spline calculations
  const PMEFloat ONEHALF    = 0.5;
  const PMEFloat ONETHIRD   = (1.0 / 3.0);
  const PMEFloat ONE        = 1.0;
  const PMEFloat THREE      = 3.0;
#if (PME_ORDER == 6)
  const PMEFloat ONEQUARTER = 0.25;
  const PMEFloat ONEFIFTH   = 0.2;
  const PMEFloat TWO        = 2.0;
  const PMEFloat FOUR       = 4.0;
  const PMEFloat FIVE       = 5.0;
#endif

  // Iterate through the simulation cell according to the order of atoms in the topology
  int pos = blockIdx.x * BATCHSIZE;

  // Read atom data and generate spline weights
  while (pos < cSim.atoms) {
    int pos1 = pos + bchidx;

    // Read charges and intialize forces.  This utilizes the last 256 of 1024 threads.
    if (pos1 < cSim.atoms && threadIdx.x >= 3*BATCHSIZE) {
#ifdef AFE_REGION1
      int AFE_region = cSim.pImageTIRegion[pos1];
      bool AFE_Mask  = (AFE_region & 4);
      PMEFloat charge = (!AFE_Mask) * cSim.pAtomChargeSP[pos1];
#elif defined(AFE_REGION2)
      int AFE_region = cSim.pImageTIRegion[pos1];
      bool AFE_Mask  = (AFE_region & 2);
      PMEFloat charge = (!AFE_Mask) * cSim.pAtomChargeSP[pos1];
#else
      PMEFloat charge = cSim.pAtomChargeSP[pos1];
#endif
#ifdef use_DPFP
      q[threadIdx.x - 3*BATCHSIZE] = charge * FORCESCALE;
#else
      q[threadIdx.x - 3*BATCHSIZE] = charge * FORCESCALEF;
#endif
      frc[threadIdx.x - 3*BATCHSIZE] = cSim.pNBForceXAccumulator[pos1];
      frc[threadIdx.x - 2*BATCHSIZE] = cSim.pNBForceYAccumulator[pos1];
      frc[threadIdx.x -   BATCHSIZE] = cSim.pNBForceZAccumulator[pos1];
    }

    // Read coordinates, compute the B-spline coefficients procedurally, then store
    // results in __shared__ memory.  This utilizes the first 768 of 1024 threads.
    if (threadIdx.x < 3*BATCHSIZE) {
      PMEFloat fx;
      if (pos1 < cSim.atoms) {
        int segment = (threadIdx.x >> BATCH_BITS);
        if (segment == 0) {
          fx = cSim.pFractX[pos1];
        }
        else if (segment == 1) {
          fx = cSim.pFractY[pos1];
        }
        else if (segment == 2) {
          fx = cSim.pFractZ[pos1];
        }
      }
      int ix = int(fx);
      int ptstore = (ndim * BATCHSIZE * cSim.pmeOrder) + bchidx;
      fx -= ix;
      ix -= cSim.orderMinusOne;
      ix += (ix < 0) * nfftval;
      PMEFloat tx0, tx1, tx2, tx3;
      tx0 = ONE - fx;
      tx1 = fx;
      tx2 = ONEHALF * fx * tx1;
      tx0 = ONEHALF * (ONE - fx) * tx0;
      tx1 = ONE - tx0 - tx2;
#if (PME_ORDER == 4)
      sAtom_dtxyz[ptstore              ] = -tx0;
      sAtom_dtxyz[ptstore +   BATCHSIZE] = tx0 - tx1;
      sAtom_dtxyz[ptstore + 2*BATCHSIZE] = tx1 - tx2;
      sAtom_dtxyz[ptstore + 3*BATCHSIZE] = tx2;
#endif
      tx3 = ONETHIRD * fx * tx2;
      tx2 = ONETHIRD * ((fx + ONE) * tx1 + (THREE - fx) * tx2);
      tx0 = ONETHIRD *  (ONE - fx) * tx0;
      tx1 = ONE - tx0 - tx2 - tx3;
#if (PME_ORDER == 6)
      PMEFloat tx4, tx5;
      tx4 = ONEQUARTER * fx * tx3;
      tx3 = ONEQUARTER * ((fx + ONE)   * tx2 + (FOUR  - fx) * tx3);
      tx2 = ONEQUARTER * ((fx + TWO)   * tx1 + (THREE - fx) * tx2);
      tx1 = ONEQUARTER * ((fx + THREE) * tx0 + (TWO   - fx) * tx1);
      tx0 = ONEQUARTER * (ONE - fx) * tx0;
      sAtom_dtxyz[ptstore              ] = -tx0;
      sAtom_dtxyz[ptstore +   BATCHSIZE] = tx0 - tx1;
      sAtom_dtxyz[ptstore + 2*BATCHSIZE] = tx1 - tx2;
      sAtom_dtxyz[ptstore + 3*BATCHSIZE] = tx2 - tx3;
      sAtom_dtxyz[ptstore + 4*BATCHSIZE] = tx3 - tx4;
      sAtom_dtxyz[ptstore + 5*BATCHSIZE] = tx4;
      tx5 = ONEFIFTH * fx * tx4;
      tx4 = ONEFIFTH * ((fx + ONE)   * tx3 + (FIVE  - fx) * tx4);
      tx3 = ONEFIFTH * ((fx + TWO)   * tx2 + (FOUR  - fx) * tx3);
      tx2 = ONEFIFTH * ((fx + THREE) * tx1 + (THREE - fx) * tx2);
      tx1 = ONEFIFTH * ((fx + FOUR)  * tx0 + (TWO   - fx) * tx1);
      tx0 = ONEFIFTH * (ONE - fx) * tx0;
#endif
      sAtom_ixyz[threadIdx.x] = ix;
      sAtom_txyz[ptstore              ] = tx0;
      sAtom_txyz[ptstore +   BATCHSIZE] = tx1;
      sAtom_txyz[ptstore + 2*BATCHSIZE] = tx2;
      sAtom_txyz[ptstore + 3*BATCHSIZE] = tx3;
#if (PME_ORDER == 6)
      sAtom_txyz[ptstore + 4*BATCHSIZE] = tx4;
      sAtom_txyz[ptstore + 5*BATCHSIZE] = tx5;
#endif
    }
    __syncthreads();

    // Interpolate forces for a group of atoms.  This is the most
    // intensive part of the kernel and will make use of all threads.
    // The only check here is against the total number of atoms in
    // the simulation: BATCHSIZE is a multiple of 32 in this kernel.
    int workidx = ((threadIdx.x >> (GRID_BITS + 2)) << GRID_BITS) +
                 (threadIdx.x & GRID_BITS_MASK);
    if (pos + workidx < cSim.atoms) {

      // Calculate the starting mesh point
      int iMeshStart = sAtom_ixyz[              workidx];
#if (PME_ORDER == 4)
      int jMeshStart = sAtom_ixyz[  BATCHSIZE + workidx];
      int kMeshLevel = sAtom_ixyz[2*BATCHSIZE + workidx] + zOffset;
#elif (PME_ORDER == 6)
      int jMeshStart = sAtom_ixyz[  BATCHSIZE + workidx] + yOffset;
      int kMeshLevel = sAtom_ixyz[2*BATCHSIZE + workidx] + zOffset;
#endif
      kMeshLevel -= (kMeshLevel >= cSim.nfft3) * cSim.nfft3;
      PMEFloat fx = (PMEFloat)0.0;
      PMEFloat fy = (PMEFloat)0.0;
      PMEFloat fz = (PMEFloat)0.0;
      int j;
#if (PME_ORDER == 4)
      for (j = 0; j < 4; j++) {
#elif (PME_ORDER == 6)
      for (j = yOffset; j < 6; j++) {
#endif
        // Lay out the mesh indices iidx, jidx, and kidx, with iidx incrementing,
        // jidx fixed for this particular run along the X dimension, and kidx
        // holding the index into the linearized mesh data (less the changing iidx).
        int iidx = iMeshStart;
        int jidx = jMeshStart + j;
        jidx -= (jidx >= cSim.nfft2) * cSim.nfft2;
        int kidx = (kMeshLevel * cSim.nfft1xnfft2) + (jidx * cSim.nfft1);
        int k;
        PMEFloat xsum = (PMEFloat)0.0;
        PMEFloat dxsum = (PMEFloat)0.0;
#if (PME_ORDER == 4)
#pragma unroll 4
#elif (PME_ORDER == 6)
#pragma unroll 6
#endif
        for (k = 0; k < PME_ORDER; k++) {
#ifdef use_DPFP
          int2 i2term = tex1Dfetch<int2>(cSim.texXYZ_q, kidx + iidx);
          PMEFloat qterm = __hiloint2double(i2term.y, i2term.x);
#else
          PMEFloat qterm = tex1Dfetch<float>(cSim.texXYZ_q, kidx + iidx);
#endif
          xsum  +=  sAtom_txyz[k*BATCHSIZE + workidx] * qterm;
          dxsum += sAtom_dtxyz[k*BATCHSIZE + workidx] * qterm;
          iidx++;
          iidx -= (iidx >= cSim.nfft1) * cSim.nfft1;
        }

        // Now re-use jidx and kidx for folding in additional B-spline coefficients
        jidx = (cSim.pmeOrder + j)*BATCHSIZE + workidx;
        kidx = (2*cSim.pmeOrder + zOffset)*BATCHSIZE + workidx;
        fx -= dxsum * sAtom_txyz[jidx]  * sAtom_txyz[kidx];
        fy -= xsum  * sAtom_dtxyz[jidx] * sAtom_txyz[kidx];
        fz -= xsum  * sAtom_txyz[jidx]  * sAtom_dtxyz[kidx];
#if (PME_ORDER == 4)
      }
#elif (PME_ORDER == 6)
      }

      // Repeat the above loop for the second slab that this warp is charged to
      // handle within the 6-slab stack.  Each warp does 1.5 slabs (4*1.5 = 6).
      // If the first slab was a full slab, this one will be a half slab.  If
      // the first slab was a half slab, this will be a full slab.
      kMeshLevel++;
      kMeshLevel -= (kMeshLevel >= cSim.nfft3) * cSim.nfft3;
      for (j = 0; j < yOffset+3; j++) {
        int iidx = iMeshStart;
        int jidx = jMeshStart + j;
        jidx -= (jidx >= cSim.nfft2) * cSim.nfft2;
        int kidx = (kMeshLevel * cSim.nfft1xnfft2) + (jidx * cSim.nfft1);
        int k;
        PMEFloat xsum = (PMEFloat)0.0;
        PMEFloat dxsum = (PMEFloat)0.0;
#pragma unroll 6
        for (k = 0; k < 6; k++) {
#ifdef use_DPFP
          int2 i2term = tex1Dfetch<int2>(cSim.texXYZ_q, kidx + iidx);
          PMEFloat qterm = __hiloint2double(i2term.y, i2term.x);
#else
          PMEFloat qterm = tex1Dfetch<float>(cSim.texXYZ_q, kidx + iidx);
#endif
          xsum  +=  sAtom_txyz[k*BATCHSIZE + workidx] * qterm;
          dxsum += sAtom_dtxyz[k*BATCHSIZE + workidx] * qterm;
          iidx++;
          iidx -= (iidx >= cSim.nfft1) * cSim.nfft1;
        }

        // Now re-use jidx and kidx for folding in additional B-spline coefficients
        jidx = (cSim.pmeOrder + j)*BATCHSIZE + workidx;
        kidx = (2*cSim.pmeOrder + zOffset)*BATCHSIZE + workidx;
        fx -= dxsum * sAtom_txyz[jidx]  * sAtom_txyz[kidx];
        fy -= xsum  * sAtom_dtxyz[jidx] * sAtom_txyz[kidx];
        fz -= xsum  * sAtom_txyz[jidx]  * sAtom_dtxyz[kidx];
      }
#endif

      // Scale the force components with the unit cell dimensions and mesh size.  This
      // general-purpose code can transform the forces in-place thanks to the unit cell
      // transformation matrix being lower-triangular.
#ifdef PME_VIRIAL
      fz = sRecipf[6] * fx + sRecipf[7] * fy + sRecipf[8] * fz;
      fy = sRecipf[3] * fx + sRecipf[4] * fy;
      fx = sRecipf[0] * fx;
#else
      fz = cSim.recipf[2][0] * fx + cSim.recipf[2][1] * fy + cSim.recipf[2][2] * fz;
      fy = cSim.recipf[1][0] * fx + cSim.recipf[1][1] * fy;
      fx = cSim.recipf[0][0] * fx;
#endif
      fx *= cSim.nfft1 * q[workidx];
      fy *= cSim.nfft2 * q[workidx];
      fz *= cSim.nfft3 * q[workidx];
#ifdef use_DPFP
      atomicAdd((unsigned long long int*)&frc[              workidx], llrint(fx));
      atomicAdd((unsigned long long int*)&frc[  BATCHSIZE + workidx], llrint(fy));
      atomicAdd((unsigned long long int*)&frc[2*BATCHSIZE + workidx], llrint(fz));
#else
      atomicAdd((unsigned long long int*)&frc[              workidx], fast_llrintf(fx));
      atomicAdd((unsigned long long int*)&frc[  BATCHSIZE + workidx], fast_llrintf(fy));
      atomicAdd((unsigned long long int*)&frc[2*BATCHSIZE + workidx], fast_llrintf(fz));
#endif
    }

    // Dump results to global memory
    __syncthreads();
    if (pos1 < cSim.atoms) {
      if (threadIdx.x < BATCHSIZE) {
        cSim.pNBForceXAccumulator[pos1] = frc[threadIdx.x];
      }
      else if (threadIdx.x < 2*BATCHSIZE) {
        cSim.pNBForceYAccumulator[pos1] = frc[threadIdx.x];
      }
      else if (threadIdx.x < 3*BATCHSIZE) {
        cSim.pNBForceZAccumulator[pos1] = frc[threadIdx.x];
      }
    }

    // Advance to the next batch of atoms
    pos += gridDim.x * BATCHSIZE;
    __syncthreads();
  }
}
#undef BATCHSIZE
#undef BATCH_BITS_MASK
#undef BATCH_BITS
