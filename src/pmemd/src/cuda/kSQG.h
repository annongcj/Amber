#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------
{
  // Each block does one of the pre-arranged sub-regions based on
  // data previously laid out in global memory.
  __shared__ double sRecip[9], sUcell[9], cellorig[3];
  __shared__ int istart, iend, xcells, ycells, zcells;
  __shared__ uint initgrid[32];
  __shared__ int stage[32];
#ifdef PME_NTP
  if (threadIdx.x < 9) {
    sRecip[threadIdx.x] = cSim.pNTPData->recip[threadIdx.x];
    sUcell[threadIdx.x] = cSim.pNTPData->ucell[threadIdx.x];
  }
#endif
  if (threadIdx.x < 3) {
#ifdef PME_IS_ORTHOGONAL
#  ifdef PME_NTP
    cellorig[threadIdx.x] = sRecip[3*threadIdx.x]*cSim.cellOffset[blockIdx.x][threadIdx.x];
#  endif
#else
    cellorig[threadIdx.x] = cSim.cellOffset[blockIdx.x][threadIdx.x];
#endif
  }
  if (threadIdx.x == 32) {
    uint2 clim = cSim.pNLNonbondCellStartEnd[blockIdx.x];
    istart = clim.x;
    iend = clim.y;
  }
  uint warpidx = threadIdx.x >> 5;
  uint tgx = threadIdx.x - (warpidx << 5);
#if PME_INTERPOLATION_ORDER == 4
  // The pattern for fourth order interpolation is fairly simple.  Four is a clean factor of
  // 32, so the only trick is to pad the regular grid in its two most rapidly varying
  // dimensions (x and y) to ensure that the memory accesses get staggered across all 32
  // banks during the application of each particle's charge onto the mesh.
  //
  //    0    1    2    3      400  401  402  403      800  801  802  803
  //   20   21   22   23 -->  420  421  422  423 -->  820  821  822  823 --> and once more
  //   40   41   42   43      440  441  442  443      840  841  842  843
  //   60   61   62   63      460  461  462  463      860  861  862  863
  //
  // This requires a staging mesh with dimension 20 in x and y.
  uint zdim = threadIdx.x >> 4;
  uint slabcount = (threadIdx.x - (zdim << 4));
  uint ydim = slabcount >> 2;
  uint allPassOffset = zdim*400 + ydim*20 + (slabcount - (ydim << 2));
#else  // PME_INTERPOLATION_ORDER == 6
  // The pattern for sixth order interpolation is more complex.  Starting from the initial
  // grid seed point (the minimum ix, iy, iz grid point affected by the particle, which will
  // be computed by the warps doing the B-spline coefficients), the assignment warp must
  // skip forward by the following amounts:
  //
  //    0    1    2    3    4    5      389  390  391  392  393  394
  //   19   20   21   22   23   24      408  409  410  411  412  413
  //   38   39   40   41   42   43 -->  427  428  429  430  431  432 --> four more times
  //   57   58   59   60   61   62      446  447  448  449  450  451
  //   76   77   78   79   80   81      465  466  467  468  469  470
  //   95                               484
  //
  // Following the initial six passes, the counts work differently to fill in gaps left behind
  // (only 31 of 36 spaces on each slab were filled).
  //
  //   96   97   98   99  100
  //  485  486  487  488  489
  //  874  875  876  877  878
  // 1263 1264 1265 1266 1267
  // 1652 1653 1654 1655 1656
  // 2041 2042 2043 2044 2045
  //
  // Each slab of the padded grid measures 19 in the fastest incrementing dimension (x) by 20
  // (y), and each slab is further padded by an additional 9 spaces so that the memory bank
  // owning position (i,j) on the next slab up in (z) will be five banks ahead.  If bank 0 owns
  // position (i,j,k) on slab k, bank 5 will own position (i,j,k+1).  The following is designed
  // to work with a 128 thread configuration.

  // Threads 32-63 will map out initial offsets and combos
  // for the first six passes and also the last pass.
  uint isblank = (tgx == 31)*388;
  uint notblank = (tgx != 31);
  uint sixPassOffset = isblank + notblank*(tgx + 13*(tgx / 6));
  isblank = (tgx >= 30)*388;
  notblank = (tgx < 30);
  uint lastPassOffset = isblank + notblank*(96 + tgx + 384*(tgx / 5));
  uint xSixPass = tgx % 6;
  uint ySixPass = tgx / 6;
  uint xLastPass = tgx % 5 + 1;
  uint yLastPass = tgx / 5;
#endif
  __syncthreads();
}
