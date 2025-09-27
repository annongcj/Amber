#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// The following code is included by kPMEInterpolation.cu to generate multiple forms of the
// kernel kPMEFillChargeGridBuffer(...)_kernel.
//---------------------------------------------------------------------------------------------
#define BATCHSIZE        30
#define BATCH_BITS        5
#define BATCH_BITS_MASK  31
#if !defined(AMBER_PLATFORM_AMD)
#  define VOLATILE volatile
#else
#  define VOLATILE
#endif
{
  // The ordering of crdq is as follows:
  // [ x1, x2, ..., x(BATCHSIZE), y1, y2, ..., z(BATCHSIZE), q1, q2, ..., q(BATCHSIZE) ]
  __shared__ VOLATILE PMEFloat crdq[4*BATCHSIZE];

  // The starting mesh indices (minimum x, y, and z stencil indices) are stored in the
  // order [ ix1, ix2, ix3, ..., ix(BATCHSIZE), iy1, iy2, ..., iy(BATCHSIZE), iz1, ... ],
  // mirroring crdq.
  __shared__ VOLATILE int sAtom_ixyz[3*BATCHSIZE];

  // The B-spline cofficients arrays are stored in padded arrays to make sure that
  // the 12 (or, 18, for 6th order interpolation) coefficients all fall on different
  // __shared__ memory banks.  For 4th order interpolation, the order is as follows:
  //
  // [ (atom 1, point 1x), (atom 1, point 2x), (atom1, point 3x), (atom1, point4x),
  //   (atom 2, point 1x), ... (atom BATCHSIZE, point 4x), (atom 1, point 1y),
  //   (atom 1, point 2y), (atom1, point 3y), (atom1, point 4y), (atom 2, point 1y), ...
  //   (atom BATCHSIZE, point 4y), ... (atom BATCHSIZE, point 4z) ]
#if (PME_ORDER == 4)
  __shared__ VOLATILE PMEFloat sAtom_txyz[12*BATCHSIZE];
#elif (PME_ORDER == 6)
  __shared__ VOLATILE PMEFloat sAtom_txyz[18*BATCHSIZE];
#endif

  // There are multiple warps
#ifdef AMBER_PLATFORM_AMD_WARP64
  // HIP-TODO: Support warp size = 64 (also see P2M_THREADS_PER_BLOCK, BATCHSIZE, etc)
  // As a workaround here we split work into logical "warps" of 32 threads (there are no
  // cross-lane operations so it is safe to use 3*32 = 96 = P2M_THREADS_PER_BLOCK on AMD GPUs)
  const int LOGICAL_GRID = 32;
  const int tgx = threadIdx.x % LOGICAL_GRID;
#else
  const int tgx = threadIdx.x & GRID_BITS_MASK;
#endif

  // Threads will initially cooperate to read atom coordinates and charges, then compute
  // B-spline coefficients in three dimensions.  Later, warps will cooperate to map each
  // atom to the mesh.
  const int ndim = threadIdx.x / BATCHSIZE;
  const int nfftval = (ndim == 0)*cSim.nfft1 + (ndim == 1)*cSim.nfft2 + (ndim == 2)*cSim.nfft3;
#ifdef AMBER_PLATFORM_AMD_WARP64
  // HIP-TODO: Remove after implementing support of warp size = 64
  const int nmap = BATCHSIZE / (P2M_THREADS_PER_BLOCK / LOGICAL_GRID);
#else
  const int nmap = BATCHSIZE / (P2M_THREADS_PER_BLOCK >> GRID_BITS);
#endif

  // Determine grid offsets and set constants for the B-spline calculations
  const PMEFloat ONEHALF    = 0.5;
  const PMEFloat ONETHIRD   = (1.0 / 3.0);
  const PMEFloat ONE        = 1.0;
  const PMEFloat THREE      = 3.0;
#if (PME_ORDER == 4)
  const unsigned int iOffsetX = tgx & 0x03;
  const unsigned int iOffsetY = (tgx & 0x0f) >> 2;
  const unsigned int iOffsetZ = tgx >> 4;
#elif (PME_ORDER == 6)
  const PMEFloat ONEQUARTER = 0.25;
  const PMEFloat ONEFIFTH   = 0.2;
  const PMEFloat TWO        = 2.0;
  const PMEFloat FOUR       = 4.0;
  const PMEFloat FIVE       = 5.0;
  const unsigned int iOffsetX = tgx % 6;
  const unsigned int iOffsetY = tgx / 6;
  const unsigned int iOffsetZ = 0;
  const unsigned int jOffsetX = (tgx & 0x03) + 2;
  const unsigned int jOffsetY = 5;
  const unsigned int jOffsetZ = tgx >> 2;
#endif

  // Iterate through the simulation cell according to the order of atoms in the topology
  int pos = blockIdx.x * BATCHSIZE;

  // Read atom data and generate spline weights
#if !defined(AMBER_PLATFORM_AMD)
  while (pos < cSim.atoms) {
#endif
    int atomid = (threadIdx.x & BATCH_BITS_MASK);
    int pos1 = pos + atomid;

    // Read coordinates and charges.  This attempts to make use of warp specialization.
    if ((pos1 < cSim.atoms) && (atomid < BATCHSIZE)) {
      int segment = (threadIdx.x >> BATCH_BITS);
      if (segment == 0) {
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
        PMEFloat fx = cSim.pFractX[pos1];
        crdq[              atomid] = fx;
        crdq[3*BATCHSIZE + atomid] = charge;
      }
      else if (segment == 1) {
        PMEFloat fy = cSim.pFractY[pos1];
        crdq[  BATCHSIZE + atomid] = fy;
      }
      else if (segment == 2) {
        PMEFloat fz = cSim.pFractZ[pos1];
        crdq[2*BATCHSIZE + atomid] = fz;
      }
    }
    __syncthreads();

    // Compute the B-spline coefficients procedurally, then store results
    // in __shared__ memory.  This utilizes 720 / 768 threads.
    if (threadIdx.x < 3*BATCHSIZE) {
      PMEFloat fx = crdq[threadIdx.x];
      int ix = int(fx);
      fx -= ix;
      ix -= cSim.orderMinusOne;
      ix += (ix < 0) * nfftval;
      PMEFloat tx0, tx1, tx2, tx3;
      tx0 = ONE - fx;
      tx1 = fx;
      tx2 = ONEHALF * fx * tx1;
      tx0 = ONEHALF * (ONE - fx) * tx0;
      tx1 = ONE - tx0 - tx2;
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
      tx5 = ONEFIFTH * fx * tx4;
      tx4 = ONEFIFTH * ((fx + ONE)   * tx3 + (FIVE  - fx) * tx4);
      tx3 = ONEFIFTH * ((fx + TWO)   * tx2 + (FOUR  - fx) * tx3);
      tx2 = ONEFIFTH * ((fx + THREE) * tx1 + (THREE - fx) * tx2);
      tx1 = ONEFIFTH * ((fx + FOUR)  * tx0 + (TWO   - fx) * tx1);
      tx0 = ONEFIFTH * (ONE - fx) * tx0;
#endif
      // Store results in shared memory
      sAtom_ixyz[threadIdx.x] = ix;
      int ptstore = cSim.pmeOrder * threadIdx.x;
      if (threadIdx.x < BATCHSIZE) {
        PMEFloat charge = crdq[threadIdx.x + 3*BATCHSIZE];
        tx0 *= charge;
        tx1 *= charge;
        tx2 *= charge;
        tx3 *= charge;
#if (PME_ORDER == 6)
        tx4 *= charge;
        tx5 *= charge;
#endif
      }
      sAtom_txyz[ptstore    ] = tx0;
      sAtom_txyz[ptstore + 1] = tx1;
      sAtom_txyz[ptstore + 2] = tx2;
      sAtom_txyz[ptstore + 3] = tx3;
#if (PME_ORDER == 6)
      sAtom_txyz[ptstore + 4] = tx4;
      sAtom_txyz[ptstore + 5] = tx5;
#endif
    }
    __syncthreads();

    // Interpolate up to nmap atoms onto the mesh.  This is the most intensive
    // part of the kernel and will make full use of all threads.
#ifdef AMBER_PLATFORM_AMD_WARP64
    // HIP-TODO: Remove after implementing support of warp size = 64
    pos1 = (threadIdx.x / LOGICAL_GRID) * nmap;
#else
    pos1 = (threadIdx.x >> GRID_BITS) * nmap;
#endif
    int lastAtom = min(pos1 + nmap, cSim.atoms - pos);
    int iOffsetXld = (cSim.pmeOrder *                 pos1) + iOffsetX;
    int iOffsetYld = (cSim.pmeOrder * (pos1 +   BATCHSIZE)) + iOffsetY;
    int iOffsetZld = (cSim.pmeOrder * (pos1 + 2*BATCHSIZE)) + iOffsetZ;
#if (PME_ORDER == 6)
    int jOffsetXld = (cSim.pmeOrder *                 pos1) + jOffsetX;
    int jOffsetYld = (cSim.pmeOrder * (pos1 +   BATCHSIZE)) + jOffsetY;
    int jOffsetZld = (cSim.pmeOrder * (pos1 + 2*BATCHSIZE)) + jOffsetZ;
#endif
    while (pos1 < lastAtom) {

      // Skip if no charge
      if (fabs(crdq[3*BATCHSIZE + pos1]) >= (PMEFloat)1.0e-8) {

        // Calculate values
        int ix = sAtom_ixyz[              pos1] + iOffsetX;
        int iy = sAtom_ixyz[  BATCHSIZE + pos1] + iOffsetY;
        int iz = sAtom_ixyz[2*BATCHSIZE + pos1] + iOffsetZ;

        // Insure coordinates stay in bounds
        ix -= (ix >= cSim.nfft1) * cSim.nfft1;
        iy -= (iy >= cSim.nfft2) * cSim.nfft2;

        // Calculate the XY component of the interpolation values and destinations
#if (PME_ORDER == 4)
        int gposxy = ((ix & 0x3) + (iy & 0x3)*4) +
                     (((ix >> 2) << 4) + ((iy >> 2) << 2)*cSim.nfft1)*cSim.nfft3;
#elif (PME_ORDER == 6)
        int gposxy = (ix & 0x3) + (iy << 2) + ((ix >> 2) << 2)*cSim.nfft2;
#endif
        PMEFloat sAtom_txy = sAtom_txyz[iOffsetXld] * sAtom_txyz[iOffsetYld] * LATTICESCALEF;
        int level;
#if (PME_ORDER == 4)
        for (level = 0; level < 4; level += 2) {
#elif (PME_ORDER == 6)
        for (level = 0; level < 6; level++) {
#endif
          iz += level;
          iz -= (iz >= cSim.nfft3) * cSim.nfft3;
#if (PME_ORDER == 4)
          int gposxyz = gposxy + (iz << 4);
#elif (PME_ORDER == 6)
          int gposxyz = gposxy + iz*cSim.nfft1xnfft2;
#endif
#ifdef use_DPFP
          unsigned long long int value = llitoulli(llrint(sAtom_txy *
                                                          sAtom_txyz[iOffsetZld + level]));
          atomicAdd((unsigned long long int*)&cSim.plliXYZ_q[gposxyz], value);
#else
          int value = __float2int_rn(sAtom_txy * sAtom_txyz[iOffsetZld + level]);
          atomicAdd(&cSim.plliXYZ_q[gposxyz], value);
#endif
#if (PME_ORDER == 4)
        }
#elif (PME_ORDER == 6)
        }
#endif

#if (PME_ORDER == 6)
        // Sixth-order interpolation has left out 24 pieces of the stencil.
        // Fill those in now.
        if (tgx < 24) {
          ix = sAtom_ixyz[              pos1] + jOffsetX;
          iy = sAtom_ixyz[  BATCHSIZE + pos1] + jOffsetY;
          iz = sAtom_ixyz[2*BATCHSIZE + pos1] + jOffsetZ;
          ix -= (ix >= cSim.nfft1) * cSim.nfft1;
          iy -= (iy >= cSim.nfft2) * cSim.nfft2;
          iz -= (iz >= cSim.nfft3) * cSim.nfft3;
          int gpos = (ix & 0x3) + (iy << 2) + ((ix >> 2) << 2)*cSim.nfft2 +
                     iz*cSim.nfft1xnfft2;
#ifdef use_DPFP
          unsigned long long int value = llitoulli(llrint(LATTICESCALE *
                                                          sAtom_txyz[jOffsetXld] *
                                                          sAtom_txyz[jOffsetYld] *
                                                          sAtom_txyz[jOffsetZld]));
          atomicAdd((unsigned long long int*)&cSim.plliXYZ_q[gpos], value);
#else
          int value = __float2int_rn(LATTICESCALEF * sAtom_txyz[jOffsetXld] *
                                     sAtom_txyz[jOffsetYld] * sAtom_txyz[jOffsetZld]);
          atomicAdd(&cSim.plliXYZ_q[gpos], value);
#endif
        }
        jOffsetXld += cSim.pmeOrder;
        jOffsetYld += cSim.pmeOrder;
        jOffsetZld += cSim.pmeOrder;
#endif
      }

      // Increment the mapping atom number within this block
      iOffsetXld += cSim.pmeOrder;
      iOffsetYld += cSim.pmeOrder;
      iOffsetZld += cSim.pmeOrder;
      pos1++;
    }

#if !defined(AMBER_PLATFORM_AMD)
    // Advance to the next batch of atoms
    pos += gridDim.x * BATCHSIZE;
    __syncthreads();
  }
#endif
}
#undef BATCHSIZE
#undef BATCH_BITS
#undef BATCH_BITS_MASK
#undef VOLATILE
