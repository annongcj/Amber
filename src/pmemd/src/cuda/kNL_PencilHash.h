#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// This code is #include'd multiple times by kNeighborList.cu as the kNLGenerateSpatialHash
// kernel.  Critical #defines: PME_VIRIAL, for Berendsen barostat applications; YSHIFT_PENCILS,
// for situations in which the box angle made by the XY and XZ planes is less than 7/12 pi and
// the hexagonal prism cells are staggered depending on the value of Z. 
//---------------------------------------------------------------------------------------------
{
#ifdef PME_VIRIAL
  __shared__ double sRecip[9], sUcell[9];
  if (threadIdx.x < 9) {
    sRecip[threadIdx.x] = cSim.pNTPData->recip[threadIdx.x];
  }
  else if (threadIdx.x >= 32 && threadIdx.x < 41) {
    sUcell[threadIdx.x - 32] = cSim.pNTPData->ucell[threadIdx.x - 32];
  }
  __syncthreads();
#endif
  unsigned int pos       = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  while (pos < cSim.atoms) {

    // Double precision in global is reduced to fp32 in SP*P modes in
    // registers--this is fine because the spatial hash assignments can
    // tolerate this level of imprecision, even in very large numbers.
    double x = cSim.pHcmbX[pos];
    double y = cSim.pHcmbY[pos];
    double z = cSim.pHcmbZ[pos];

    // Orthogonal/nonorthogonal handled in the same code
#ifdef PME_VIRIAL
    double dfx = sRecip[0]*x + sRecip[3]*y + sRecip[6]*z;
    double dfy =               sRecip[4]*y + sRecip[7]*z;
    double dfz =                             sRecip[8]*z;
#else
    double dfx = cSim.recip[0][0]*x + cSim.recip[1][0]*y + cSim.recip[2][0]*z;
    double dfy =                      cSim.recip[1][1]*y + cSim.recip[2][1]*z;
    double dfz =                                           cSim.recip[2][2]*z;
#endif

    // Re-image all fractional coordinates to within [0.0, 1.0)
    dfx = dfx - round(dfx) + 0.5;
    dfy = dfy - round(dfy) + 0.5;
    dfz = dfz - round(dfz) + 0.5;
    dfx = (dfx < 1.0) ? dfx : 0.0;
    dfy = (dfy < 1.0) ? dfy : 0.0;
    dfz = (dfz < 1.0) ? dfz : 0.0;
#ifdef CALCULATE_HASH_CELL
    PMEFloat fx = dfx;
    PMEFloat fy = dfy;
    PMEFloat fz = dfz;

    // Decide on the Y and Z pencil cell locations
#  ifdef PME_VIRIAL
    double dry = sUcell[4]*dfy + sUcell[5]*dfz;
    double drz =                 sUcell[8]*dfz;
#  else
    double dry = cSim.ucell[1][1]*dfy + cSim.ucell[1][2]*dfz;
    double drz =                        cSim.ucell[2][2]*dfz;
#  endif
    PMEFloat ry = dry;
    PMEFloat rz = drz;
#  ifdef YSHIFT_PENCILS
    // Harder case: the pencils are arranged in staggered fashion, with each stack of
    // them offset by 0.5 lengths in Y depending on their position relative to the
    // stack starting at the origin.
    int iz           = fz * cSim.nzpencilsf;
    int iy           = (fy * cSim.nypencilsf) - 0.5*(iz & 1);
    PMEFloat izf     = iz;
    PMEFloat izp1f   = izf + (PMEFloat)1.0;
    PMEFloat iyfA    = (PMEFloat)iy + (PMEFloat)0.5 * (PMEFloat)(iz & 0x1);
    PMEFloat iyfB    = (PMEFloat)iy + (PMEFloat)0.5 * (PMEFloat)((iz + 1) & 0x1);
    PMEFloat iyp1fA  = iyfA + (PMEFloat)1.0;
    PMEFloat iyp1fB  = iyfB + (PMEFloat)1.0;

    // The four candidates are now (iyfA, izf), (iyp1fA, izf), (iyfB, izfp1), (iyp1fB, izfp1).
    // Up to this point, they have been calculated as fractional coordinates expressed in the
    // units of the hash cell width (i.e. fractional coordinate 0.3333 in a grid of 18 cells
    // is 6).  Express these coordinates strictly in terms of the box fractional coordinates.
    iyfA            *= cSim.invypencils;
    iyfB            *= cSim.invypencils;
    iyp1fA          *= cSim.invypencils;
    iyp1fB          *= cSim.invypencils;
    izf             *= cSim.invzpencils;
    izp1f           *= cSim.invzpencils;
#    ifdef PME_VIRIAL
    PMEFloat rdy  = sUcell[4]*iyfA + sUcell[5]*izf - ry;
    PMEFloat rdz  =                  sUcell[8]*izf - rz;
#    else
    PMEFloat rdy  = cSim.ucellf[1][1]*iyfA + cSim.ucellf[1][2]*izf - ry;
    PMEFloat rdz  =                          cSim.ucellf[2][2]*izf - rz;    
#    endif
    PMEFloat rmin = rdy*rdy + rdz*rdz;
    int besty = iy;
    int bestz = iz;
#    ifdef PME_VIRIAL
    rdy            = sUcell[4]*iyp1fA + sUcell[5]*izf - ry;
#    else
    rdy            = cSim.ucellf[1][1]*iyp1fA + cSim.ucellf[1][2]*izf - ry;
#    endif
    PMEFloat rtest = rdy*rdy + rdz*rdz;
    if (rtest < rmin) {
      rmin = rtest;
      besty = iy + 1;
    }
#    ifdef PME_VIRIAL
    rdy             = sUcell[4]*iyfB + sUcell[5]*izp1f - ry;
    rdz             =                  sUcell[8]*izp1f - rz;
#    else
    rdy             = cSim.ucellf[1][1]*iyfB + cSim.ucellf[1][2]*izp1f - ry;
    rdz             =                          cSim.ucellf[2][2]*izp1f - rz;
#    endif
    rtest = rdy*rdy + rdz*rdz;
    if (rtest < rmin) {
      rmin = rtest;
      besty = iy;
      bestz = iz + 1;
    }
#    ifdef PME_VIRIAL
    rdy             = sUcell[4]*iyp1fB + sUcell[5]*izp1f - ry;
#    else
    rdy             = cSim.ucellf[1][1]*iyp1fB + cSim.ucellf[1][2]*izp1f - ry;
#    endif
    rtest = rdy*rdy + rdz*rdz;
    if (rtest < rmin) {
      rmin = rtest;
      besty = iy + 1;
      bestz = iz + 1;
    }
    besty += cSim.nypencils * (besty < 0);
#  else
    // Easier case: the pencils are arranged on a regular grid (pencils are laid out
    // regularly in fractional space).  Find the nearest four, translate them into real
    // space, and see which of them is truly closest.
    int iy          = fy * cSim.nypencilsf;
    int iz          = fz * cSim.nzpencilsf;
    PMEFloat iyf    = iy;
    PMEFloat izf    = iz;
    PMEFloat iyp1f  = iyf + (PMEFloat)1.0;
    PMEFloat izp1f  = izf + (PMEFloat)1.0;
#    ifdef PME_VIRAL
    iyf             = (iyf   * cSim.invypencils * sUcell[4]) - ry;
    iyp1f           = (iyp1f * cSim.invypencils * sUcell[4]) - ry;
#    else
    iyf             = (iyf   * cSim.invypencils * cSim.ucellf[1][1]) - ry;
    iyp1f           = (iyp1f * cSim.invypencils * cSim.ucellf[1][1]) - ry;
#    endif
    izf            *= cSim.invzpencils;
    izp1f          *= cSim.invzpencils;
#    ifdef PME_VIRIAL
    PMEFloat rdy  = iyf + sUcell[5]*izf;
    PMEFloat rdz  =       sUcell[8]*izf - rz;
#    else
    PMEFloat rdy  = iyf + cSim.ucellf[1][2]*izf;
    PMEFloat rdz  =       cSim.ucellf[2][2]*izf - rz;
#    endif
    PMEFloat rmin = rdy*rdy + rdz*rdz;
    int besty = iy;
    int bestz = iz;
#    ifdef PME_VIRIAL
    rdy            = iyp1f + sUcell[5]*izf;
#    else
    rdy            = iyp1f + cSim.ucellf[1][2]*izf;
#    endif
    PMEFloat rtest = rdy*rdy + rdz*rdz;
    if (rtest < rmin) {
      rmin = rtest;
      besty = iy + 1;
    }
#    ifdef PME_VIRIAL
    rdy   = iyf + sUcell[5]*izp1f;
    rdz   =       sUcell[8]*izp1f - rz;
#    else
    rdy   = iyf + cSim.ucellf[1][2]*izp1f;
    rdz   =       cSim.ucellf[2][2]*izp1f - rz;
#    endif
    rtest = rdy*rdy + rdz*rdz;
    if (rtest < rmin) {
      rmin = rtest;
      besty = iy;
      bestz = iz + 1;
    }
#    ifdef PME_VIRIAL
    rdy   = iyp1f + sUcell[5]*izp1f;
#    else
    rdy   = iyp1f + cSim.ucellf[1][2]*izp1f;
#    endif
    rtest = rdy*rdy + rdz*rdz;
    if (rtest < rmin) {
      rmin = rtest;
      besty = iy + 1;
      bestz = iz + 1;
    }
#  endif
    besty -= cSim.nypencils * (besty == cSim.nypencils);
    bestz -= cSim.nzpencils * (bestz == cSim.nzpencils);

    // Record the hash cell and image index if fractional coords are not being
    // calculated--this concludes a long procedure to determine the proper cell
    // for each atom that happens during neighbor list construction.  It is
    // skipped when the neighbor list is actually used.
    cSim.pHcmbIndex[pos]  = pos;
    unsigned int ix       = fx * XDIR_HASH_CELLS;
    unsigned int cellid   = besty + (bestz * cSim.nypencils);
    cellid                = (cellid << 2) | ((fx >= cSim.UpperCapBoundary) << 1) |
                            (fx < cSim.LowerCapBoundary);
    cSim.pHcmbCellID[pos] = cellid;
    cSim.pHcmbHash[pos]   = ix + (besty << PENCIL_IDY_SHIFT) + (bestz << PENCIL_IDZ_SHIFT);
#else  // CALCULATE_HASH_CELL
    
    // Look up the hash cell (see code above for how it was computed),
    // then snap this atom to the nearest image of that hash cell.  The
    // cell cap bits, those indicating whether the atom is in the upper
    // or lower cap and must be replicated along the X direction, are
    // omitted for forming the primary image.
    int besty = cSim.pHcmbCellID[pos];
    besty >>= 2;
    int bestz = besty / cSim.nypencils;
    besty -= bestz * cSim.nypencils;

    // Commit the fractional coordinates--this is being done
    // in preparation for non-bonded calculations.
    cSim.pHcmbFractX[pos] = (PMEFloat)(dfx * cSim.dnfft1);
    cSim.pHcmbFractY[pos] = (PMEFloat)(dfy * cSim.dnfft2);
    cSim.pHcmbFractZ[pos] = (PMEFloat)(dfz * cSim.dnfft3);

    // Snap the atom into the image of its cell center, in preparation for
    // expanding the coordinates to cover any needed imaging calculations (this
    // must be done for neighbor list setup as well as non-bonded interactions).
#  ifdef YSHIFT_PENCILS
    PMEFloat fpy = ((PMEFloat)besty + (PMEFloat)0.5*(PMEFloat)(besty & 0x1)) *
                   cSim.invypencils;
    PMEFloat fpz = (PMEFloat)bestz * cSim.invzpencils;
#  else
    PMEFloat fpy = (PMEFloat)besty * cSim.invypencils;
    PMEFloat fpz = (PMEFloat)bestz * cSim.invzpencils;
#  endif
    // Demote dfy, dfz to work in single precision and minimize fp64 work
    // (the round function will only produce -1.0, 0.0, or 1.0)
    dfy -= round((PMEFloat)dfy - fpy);
    dfz -= round((PMEFloat)dfz - fpz);
    
    // Calculate the re-imaged, real-space location (dfy and dfz may have been updated)
#  ifdef PME_VIRIAL
    double drx = sUcell[0]*dfx + sUcell[1]*dfy + sUcell[2]*dfz;
    double dry =                 sUcell[4]*dfy + sUcell[5]*dfz;
    double drz =                                 sUcell[8]*dfz;
#  else
    double drx = cSim.ucell[0][0]*dfx + cSim.ucell[0][1]*dfy + cSim.ucell[0][2]*dfz;
    double dry =                        cSim.ucell[1][1]*dfy + cSim.ucell[1][2]*dfz;
    double drz =                                               cSim.ucell[2][2]*dfz;
#  endif
    // Real-space coordinates are committed regardless of the shape of the box
    PMEFloat2 xy = {(PMEFloat)drx, (PMEFloat)dry};
    cSim.pHcmbXYSP[pos] = xy;
    cSim.pHcmbZSP[pos] = (PMEFloat)drz;
#endif // End branch decided by CALCULATE_HASH_CELL

    // Increment the position counter
    pos += increment;
  }
}
