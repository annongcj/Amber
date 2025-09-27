#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------
{
#ifdef PME_NTP
  __shared__ volatile PMEDouble sRecip[9];
#  ifdef PME_ORTHOGONAL
  __shared__ volatile PMEDouble sUcell[9];
  if (threadIdx.x < 9) {
    sUcell[threadIdx.x] = cSim.pNTPData->ucell[threadIdx.x];
  }
#  endif
  if (threadIdx.x < 9) {
    sRecip[threadIdx.x] = cSim.pNTPData->recip[threadIdx.x];
  }
  __syncthreads();
#endif
  unsigned int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  while (pos < cSim.atoms) {
    double x = cSim.pImageX[pos];
    double y = cSim.pImageY[pos];
    double z = cSim.pImageZ[pos];
    unsigned int cellID = cSim.pImageCellID[pos];
#if defined(PME_NTP)
#  if defined(PME_ORTHOGONAL)
    double dfx = x * sRecip[0];
    double dfy = y * sRecip[4];
    double dfz = z * sRecip[8];
#  else
    double dfx = x*sRecip[0] + y*sRecip[3] + z*sRecip[6];
    double dfy =               y*sRecip[4] + z*sRecip[7];
    double dfz =                             z*sRecip[8];
#  endif
#else
#  if defined(PME_ORTHOGONAL)
    double dfx = x * cSim.recip[0][0];
    double dfy = y * cSim.recip[1][1];
    double dfz = z * cSim.recip[2][2];
#  else
    double dfx = x*cSim.recip[0][0] + y*cSim.recip[1][0] + z*cSim.recip[2][0];
    double dfy =                      y*cSim.recip[1][1] + z*cSim.recip[2][1];
    double dfz =                                           z*cSim.recip[2][2];
#  endif
#endif

    // Account for minimum image convention
    dfx = (dfx - round(dfx) + 0.5);
    dfy = (dfy - round(dfy) + 0.5);
    dfz = (dfz - round(dfz) + 0.5);
    dfx = (dfx < 1.0 ? dfx : 0.0);
    dfy = (dfy < 1.0 ? dfy : 0.0);
    dfz = (dfz < 1.0 ? dfz : 0.0);

    // Can measure relative to cell edges for 3 or more cells on
    // each axis, otherwise measure relative to cell center
    double xCell = (double)(cellID & CELL_ID_MASK) * cSim.oneOverXcells;
    double yCell = (double)((cellID >> CELL_IDY_SHIFT) & CELL_ID_MASK) * cSim.oneOverYcells;
    double zCell = (double)((cellID >> CELL_IDZ_SHIFT) & CELL_ID_MASK) * cSim.oneOverZcells;
    x = dfx - xCell;
    y = dfy - yCell;
    z = dfz - zCell;
#ifdef PME_SMALLBOX
    if (x > cSim.maxCellX)  x -= 1.0;
    if (x < cSim.minCellX)  x += 1.0;
    if (y > cSim.maxCellY)  y -= 1.0;
    if (y < cSim.minCellY)  y += 1.0;
    if (z > cSim.maxCellZ)  z -= 1.0;
    if (z <= cSim.minCellZ) z += 1.0;
#else
    if (x > 0.5)   x -= 1.0;
    if (x <= -0.5) x += 1.0;
    if (y > 0.5)   y -= 1.0;
    if (y <= -0.5) y += 1.0;
    if (z > 0.5)   z -= 1.0;
    if (z <= -0.5) z += 1.0;
#endif
#ifdef PME_ORTHOGONAL
#  ifdef PME_NTP
    PMEFloat2 xy = {(PMEFloat)(sUcell[0] * x), (PMEFloat)(sUcell[4] * y)};
    cSim.pAtomXYSP[pos] = xy;
    cSim.pAtomZSP[pos]  = (PMEFloat)(sUcell[8] * z);
#  else
    PMEFloat2 xy = {(PMEFloat)(cSim.a * x), (PMEFloat)(cSim.b * y)};
    cSim.pAtomXYSP[pos] = xy;
    cSim.pAtomZSP[pos]  = (PMEFloat)(cSim.c * z);
#  endif
#else
    PMEFloat2 xy = {(PMEFloat)x, (PMEFloat)y};
    cSim.pAtomXYSP[pos] = xy;
    cSim.pAtomZSP[pos]  = z;
#endif
#ifdef COMP_FRACTIONAL_COORD
    PMEFloat fx = dfx;
    PMEFloat fy = dfy;
    PMEFloat fz = dfz;
    fx *= cSim.nfft1;
    fy *= cSim.nfft2;
    fz *= cSim.nfft3;
    cSim.pFractX[pos] = fx;
    cSim.pFractY[pos] = fy;
    cSim.pFractZ[pos] = fz;
#endif
    pos += increment;
  }
}
