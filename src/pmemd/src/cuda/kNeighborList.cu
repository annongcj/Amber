#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

#ifndef AMBER_PLATFORM_AMD
#include <cuda.h>
#endif
#include "gpu.h"
#include "bondRemapDS.h"
#include "ptxmacros.h"

// Use global instance instead of a local copy
#include "simulationConst.h"
CSIM_STO simulationConst cSim;

struct Atom {
  PMEFloat xmin;
  PMEFloat xmax;
  PMEFloat ymin;
  PMEFloat ymax;
  PMEFloat zmin;
  PMEFloat zmax;
};

#if !defined(__HIPCC_RDC__)

//---------------------------------------------------------------------------------------------
// SetkNeighborListSim: function for copying the CUDA simulation image to the device
//
// Arguments:
//   gpu:  the memory data structure storing the remapped topology in its sim attribute, plus
//         lots of other memory allocated for simulating a collection of atoms with those
//         parameters
//---------------------------------------------------------------------------------------------
void SetkNeighborListSim(gpuContext gpu)
{
  cudaError_t status;
  status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(simulationConst));
  RTERROR(status, "kNeighborList-cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

//---------------------------------------------------------------------------------------------
// GetkNeighborListSim: function for establishing cSim in the context of this library
//
// Arguments:
//   See description before SetkNeighborListSim
//---------------------------------------------------------------------------------------------
void GetkNeighborListSim(gpuContext gpu)
{
  cudaError_t status;
  status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(simulationConst));
  RTERROR(status, "kNeighborList-cudaMemcpyFromSymbol: SetSim copy from cSim failed");
}

#endif

//---------------------------------------------------------------------------------------------
// stuff: the cell hash table grid dimensions
//---------------------------------------------------------------------------------------------
struct stuff {
  double fx;
  double fy;
  double fz;
  unsigned int ix;
  unsigned int iy;
  unsigned int iz;
  unsigned int ox;
  unsigned int oy;
  unsigned int oz;
  unsigned int cellHash;
  unsigned int hash;
};

//---------------------------------------------------------------------------------------------
// kNLResetCounter_kernel: this kernel merely resets the counter associated with asynchronous
//                         neighbor list building.
//
// This is a debugging kernel.
//---------------------------------------------------------------------------------------------
__global__ void __LAUNCH_BOUNDS__(GRID, 1) kNLResetCounter_kernel()
{
  if (threadIdx.x == 0) {
    *(cSim.pNLTotalOffset)  = 0;
    *(cSim.pNLEntries)      = 0;
    cSim.pFrcBlkCounters[0] = cSim.NLBuildWarps;
  }
}

//---------------------------------------------------------------------------------------------
// kNLResetCounter: function to launch the kernel above.
//
// This is a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void kNLResetCounter()
{
  kNLResetCounter_kernel<<<GRID, 1>>>();
}

//---------------------------------------------------------------------------------------------
// kNLGenerateSpatialHash_kernel: generates the spatial hash grid based on the fractional
//                                coordinates of each atom.  Fractional coordinates are taken
//                                from the Cartesian coordinates.
//
// The spatial decomposition occurs at two levels.  First, atoms are separated into cells of
// thickness at least the non-bonded cutoff plus the neighbor list buffer zone (nb skin).
// Second, the atoms' positions within each cell are localized to within 1/4 of the cell
// width (see CELL_HASH... enumerated constants in gputypes.h).  This defines 64 x the number
// of cells worth of bins into which atoms can be placed, a number that is constructed bitwise
// using the << CELL_HASH_BITS shifting and | (bitwise or) operations.  The sCellHash table is
// not a plain 3D stack like the cell grid itself: it is a Hilbert Space Filling curve given
// in gpu.cpp (search for cellHash).  By ordering atoms along this weaving, twisting curve,
// they can be accessed linearly to give a list of atoms that are near one another even within
// the confines of the cell interior.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(THREADS_PER_BLOCK, 1)
kNLGenerateSpatialHash_kernel()
{
  __shared__ unsigned int sCellHash[CELL_HASH_CELLS];
  __shared__ PMEFloat sRecipf[9];
  unsigned int pos       = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;

  // Clear atom list/exclusion mask space counter
  if (pos == 0) {
    *(cSim.pNLTotalOffset)  = 0;
    *(cSim.pNLEntries)      = 0;
    cSim.pFrcBlkCounters[0] = cSim.NLBuildWarps;
  }

  // Read cell hash
  if (threadIdx.x < CELL_HASH_CELLS) {
    sCellHash[threadIdx.x] = cSim.pNLCellHash[threadIdx.x];
  }
  if ((cSim.ntp > 0) && (cSim.barostat == 1)) {
    if (threadIdx.x < 9) {
      sRecipf[threadIdx.x] = cSim.pNTPData->recipf[threadIdx.x];
    }
    __syncthreads();
    while (pos < cSim.atoms) {
      PMEFloat x = cSim.pImageX[pos];
      PMEFloat y = cSim.pImageY[pos];
      PMEFloat z = cSim.pImageZ[pos];

      // Orthogonal/nonorthogonal handled in the same code (3 single precision
      // multiplies and adds? Who cares and why?)
      PMEFloat fx = sRecipf[0]*x + sRecipf[3]*y + sRecipf[6]*z;
      PMEFloat fy =                sRecipf[4]*y + sRecipf[7]*z;
      PMEFloat fz =                               sRecipf[8]*z;

      // Account for minimum image convention
      fx = fx - round(fx) + (PMEFloat)0.5;
      fy = fy - round(fy) + (PMEFloat)0.5;
      fz = fz - round(fz) + (PMEFloat)0.5;
      fx = (fx < (PMEFloat)1.0 ? fx : (PMEFloat)0.0);
      fy = (fy < (PMEFloat)1.0 ? fy : (PMEFloat)0.0);
      fz = (fz < (PMEFloat)1.0 ? fz : (PMEFloat)0.0);

      // Generate box coordinates
      cSim.pImageIndex2[pos]  = pos;
      unsigned int ix = fx * cSim.xcells;
      unsigned int iy = fy * cSim.ycells;
      unsigned int iz = fz * cSim.zcells;
      cSim.pImageCellID[pos] = ix + (iy << CELL_IDY_SHIFT) + (iz << CELL_IDZ_SHIFT);
      unsigned int ox = min(CELL_HASHX - 1,
                            (unsigned int)((PMEFloat)CELL_HASHX *
                                           ((fx - ix * cSim.oneOverXcellsf) * cSim.xcells)));
      unsigned int oy = min(CELL_HASHY - 1,
                            (unsigned int)((PMEFloat)CELL_HASHY *
                                           ((fy - iy * cSim.oneOverYcellsf) * cSim.ycells))) *
                        CELL_HASHX;
      unsigned int oz = min(CELL_HASHZ - 1,
                            (unsigned int)((PMEFloat)CELL_HASHZ *
                                           ((fz - iz * cSim.oneOverZcellsf) * cSim.zcells))) *
                        CELL_HASHXY;
      unsigned int cellHash = sCellHash[ox + oy + oz];
      unsigned int hash =
        (((iz * cSim.ycells + iy) * cSim.xcells + ix) << CELL_HASH_BITS) | cellHash;
      cSim.pImageHash2[pos] = hash;
      pos += increment;
    }
  }
  else {
    __syncthreads();
    while (pos < cSim.atoms) {
      PMEFloat x = cSim.pImageX[pos];
      PMEFloat y = cSim.pImageY[pos];
      PMEFloat z = cSim.pImageZ[pos];

      // Orthogonal/nonorthogonal handled in the same code (3 single precision
      // multiplies and adds? Who cares and why?)
      PMEFloat fx = cSim.recipf[0][0]*x + cSim.recipf[1][0]*y + cSim.recipf[2][0]*z;
      PMEFloat fy =                       cSim.recipf[1][1]*y + cSim.recipf[2][1]*z;
      PMEFloat fz =                                             cSim.recipf[2][2]*z;

      // Account for minimum image convention
      fx = fx - round(fx) + (PMEFloat)0.5;
      fy = fy - round(fy) + (PMEFloat)0.5;
      fz = fz - round(fz) + (PMEFloat)0.5;
      fx = (fx < (PMEFloat)1.0 ? fx : (PMEFloat)0.0);
      fy = (fy < (PMEFloat)1.0 ? fy : (PMEFloat)0.0);
      fz = (fz < (PMEFloat)1.0 ? fz : (PMEFloat)0.0);

      // Generate box coordinates
      cSim.pImageIndex2[pos]  = pos;
      unsigned int ix        = fx * cSim.xcells;
      unsigned int iy        = fy * cSim.ycells;
      unsigned int iz        = fz * cSim.zcells;
      cSim.pImageCellID[pos] = ix + (iy << CELL_IDY_SHIFT) + (iz << CELL_IDZ_SHIFT);
      unsigned int ox = min(CELL_HASHX - 1,
                            (unsigned int)((PMEFloat)CELL_HASHX *
                                           ((fx - ix * cSim.oneOverXcellsf) * cSim.xcells)));
      unsigned int oy = min(CELL_HASHY - 1,
                            (unsigned int)((PMEFloat)CELL_HASHY *
                                           ((fy - iy * cSim.oneOverYcellsf) * cSim.ycells))) *
                        CELL_HASHX;
      unsigned int oz = min(CELL_HASHZ - 1,
                            (unsigned int)((PMEFloat)CELL_HASHZ *
                                           ((fz - iz * cSim.oneOverZcellsf) * cSim.zcells))) *
                        CELL_HASHXY;
      unsigned int cellHash = sCellHash[ox + oy + oz];
      unsigned int hash     = (((iz*cSim.ycells + iy)*cSim.xcells + ix) << CELL_HASH_BITS) |
                              cellHash;
      cSim.pImageHash2[pos] = hash;
      pos += increment;
    }
  }
  // End branch for cases of constant pressure with a Berendsen barostat
  // (that's the first branch), or constant volume simulations (that's
  // one case the second branch serves), or constant pressure simulations
  // with a Monte-Carlo barostat (also the second branch).
}

//---------------------------------------------------------------------------------------------
// kNLGenerateSpatialHash: launches the kNLGenerateSpatialHash_kernel kernel above.
//
// Arguments:
//   gpu:    here, provides a tile size for the kernel
//---------------------------------------------------------------------------------------------
extern "C" void kNLGenerateSpatialHash(gpuContext gpu)
{
  kNLGenerateSpatialHash_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
  LAUNCHERROR("kNLGenerateSpatialHash");
}

//---------------------------------------------------------------------------------------------
// kNLRemapImage_kernel: changes all pointers to all atom properties (charge, position, etc) to
//                       follow the radix sort of the atoms that occurs when a new neighborlist
//                       is built
//
// Arguments:
//   pImageIndex: buffer containing pointers to all remapped atom indexes. Provides the
//                connection from the previous position in the atom list to the position in the
//                atom list following the radix sort.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(THREADS_PER_BLOCK, 1)
kNLRemapImage_kernel(unsigned int* pImageIndex)
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  unsigned int index = 0;
  unsigned int newindex;

  if (pos < cSim.atoms) {
    index = pImageIndex[pos];
  }

  // Amazing that this memory access pattern works
  while (pos < cSim.atoms) {
    unsigned int newpos = pos + increment;
    if (newpos < cSim.atoms) {
      newindex = pImageIndex[newpos];
    }

    // Read new data
    unsigned int atom   = cSim.pImageAtom[index];
    double x            = cSim.pImageX[index];
    double y            = cSim.pImageY[index];
    double z            = cSim.pImageZ[index];
    double vx           = cSim.pImageVelX[index];
    double vy           = cSim.pImageVelY[index];
    double vz           = cSim.pImageVelZ[index];
    double lvx          = cSim.pImageLVelX[index];
    double lvy          = cSim.pImageLVelY[index];
    double lvz          = cSim.pImageLVelZ[index];
    double q            = cSim.pImageCharge[index];
    double m            = cSim.pImageMass[index];
    PMEFloat2 sigeps    = cSim.pImageSigEps[index];
    int grplist, indexPHMD;
    PMEFloat qphmd;
    PMEFloat2 qstate1, qstate2, vstate1, vstate2;
    if (cSim.iphmd == 3) {
      grplist = cSim.pImageGrplist[index];
      qstate1 = cSim.pImageQstate1[index];
      qstate2 = cSim.pImageQstate2[index];
      vstate1 = cSim.pImageVstate1[index];
      vstate2 = cSim.pImageVstate2[index];
      qphmd = cSim.pImageCharge_phmd[index];
      indexPHMD = cSim.pImageIndexPHMD[index];
    }
    int wsg;
    double x0sg,y0sg,z0sg,x1sg,y1sg,z1sg,x2sg,y2sg,z2sg,rsgx,rsgy,rsgz,fpsg,ppsg;
    if(cSim.isgld > 0) {
      wsg=cSim.pImageWsg[index];
      x0sg=cSim.pImageX0sg[index];
      y0sg=cSim.pImageY0sg[index];
      z0sg=cSim.pImageZ0sg[index];
      x1sg=cSim.pImageX1sg[index];
      y1sg=cSim.pImageY1sg[index];
      z1sg=cSim.pImageZ1sg[index];
      x2sg=cSim.pImageX2sg[index];
      y2sg=cSim.pImageY2sg[index];
      z2sg=cSim.pImageZ2sg[index];
      rsgx=cSim.pImageRsgX[index];
      rsgy=cSim.pImageRsgY[index];
      rsgz=cSim.pImageRsgZ[index];
      fpsg=cSim.pImageFPsg[index];
      ppsg=cSim.pImagePPsg[index];
    }
#ifdef use_DPFP
    long long int LJID  = cSim.pImageLJID[index];
#else
    unsigned int LJID   = cSim.pImageLJID[index];
#endif
    unsigned int cellID = cSim.pImageCellID[index];
    cSim.pImageX2[pos] = x;
    cSim.pImageY2[pos] = y;
    cSim.pImageZ2[pos] = z;
    cSim.pImageAtom2[pos] = atom;
    cSim.pImageAtomLookup[atom] = pos;

    // Have to form the tuple separately before committing it to the array
    PMEFloat2 xy = {(PMEFloat)x, (PMEFloat)y};
    cSim.pAtomXYSaveSP[pos]  = xy;
    cSim.pAtomZSaveSP[pos]   = z;
    cSim.pImageVelX2[pos]    = vx;
    cSim.pImageVelY2[pos]    = vy;
    cSim.pImageVelZ2[pos]    = vz;
    cSim.pImageLVelX2[pos]   = lvx;
    cSim.pImageLVelY2[pos]   = lvy;
    cSim.pImageLVelZ2[pos]   = lvz;
    cSim.pImageCharge2[pos]  = q;
    cSim.pAtomChargeSP[pos]  = q;
    cSim.pImageMass2[pos]    = m;
    cSim.pImageInvMass2[pos] = (m != (PMEDouble)0.0 ? (PMEDouble)1.0 / m : (PMEDouble)0.0);
    cSim.pImageSigEps2[pos]  = sigeps;
    cSim.pImageLJID2[pos]    = LJID;
    cSim.pImageCellID2[pos]  = cellID;

    if(cSim.iphmd == 3) {
      cSim.pImageGrplist2[pos] = grplist;
      cSim.pImageQstate12[pos] = qstate1;
      cSim.pImageQstate22[pos] = qstate2;
      cSim.pImageVstate12[pos] = vstate1;
      cSim.pImageVstate22[pos] = vstate2;
      cSim.pImageCharge_phmd2[pos] = qphmd;
      cSim.pImageIndexPHMD2[pos] = indexPHMD;
    }

    if(cSim.isgld > 0) {
      cSim.pImageWsg2[pos] = wsg;
      cSim.pImageX0sg2[pos] = x0sg;
      cSim.pImageY0sg2[pos] = y0sg;
      cSim.pImageZ0sg2[pos] = z0sg;
      cSim.pImageX1sg2[pos] = x1sg;
      cSim.pImageY1sg2[pos] = y1sg;
      cSim.pImageZ1sg2[pos] = z1sg;
      cSim.pImageX2sg2[pos] = x2sg;
      cSim.pImageY2sg2[pos] = y2sg;
      cSim.pImageZ2sg2[pos] = z2sg;
      cSim.pImageRsgX2[pos] = rsgx;
      cSim.pImageRsgY2[pos] = rsgy;
      cSim.pImageRsgZ2[pos] = rsgz;
      cSim.pImageFPsg2[pos] = fpsg;
      cSim.pImagePPsg2[pos] = ppsg;
    }

    // Again, need to form the tuple separately before committing it to the array
#ifdef use_DPFP
    PMEFloat2 qljid = {q, __longlong_as_double(LJID)};
#else
    PMEFloat2 qljid = {(PMEFloat)q, __uint_as_float(LJID)};
#endif
    cSim.pAtomChargeSPLJID[pos] = qljid;

    // Advance to next atom
    index = newindex;
    pos = newpos;
  }
}

//---------------------------------------------------------------------------------------------
// kNLTIRemapImage_kernel: same as kNLRemapImage_kernel, but works in the context of alchemical
//                         free energy transformations or thermodynamic integration.
//
// Arguments:
//   pImageIndex: buffer containing pointers to all remapped atom indexes. Provides the
//                connection from the previous position in the atom list to the position in the
//                atom list following the radix sort.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(THREADS_PER_BLOCK, 1)
kNLTIRemapImage_kernel(unsigned int* pImageIndex)
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  unsigned int index = 0;
  unsigned int newindex;

  if (pos < cSim.atoms) {
    index = pImageIndex[pos];
  }
  while (pos < cSim.atoms) {
    unsigned int newpos = pos + increment;
    if (newpos < cSim.atoms) {
      newindex = pImageIndex[newpos];
    }

    // Read new data
    unsigned int atom           = cSim.pImageAtom[index];
    double x                    = cSim.pImageX[index];
    double y                    = cSim.pImageY[index];
    double z                    = cSim.pImageZ[index];
    double vx                   = cSim.pImageVelX[index];
    double vy                   = cSim.pImageVelY[index];
    double vz                   = cSim.pImageVelZ[index];
    double lvx                  = cSim.pImageLVelX[index];
    double lvy                  = cSim.pImageLVelY[index];
    double lvz                  = cSim.pImageLVelZ[index];
    double q                    = cSim.pImageCharge[index];
    double m                    = cSim.pImageMass[index];
    PMEFloat2 sigeps            = cSim.pImageSigEps[index];
    int grplist, indexPHMD;
    PMEFloat qphmd;
    PMEFloat2 qstate1, qstate2, vstate1, vstate2;
    if (cSim.iphmd == 3) {
      grplist = cSim.pImageGrplist[index];
      qstate1 = cSim.pImageQstate1[index];
      qstate2 = cSim.pImageQstate2[index];
      vstate1 = cSim.pImageVstate1[index];
      vstate2 = cSim.pImageVstate2[index];
      qphmd = cSim.pImageCharge_phmd2[index];
      indexPHMD = cSim.pImageIndexPHMD[index];
    }
    int wsg;
    double x0sg,y0sg,z0sg,x1sg,y1sg,z1sg,x2sg,y2sg,z2sg,rsgx,rsgy,rsgz,fpsg,ppsg;
    if(cSim.isgld > 0) {
      wsg=cSim.pImageWsg[index];
      x0sg=cSim.pImageX0sg[index];
      y0sg=cSim.pImageY0sg[index];
      z0sg=cSim.pImageZ0sg[index];
      x1sg=cSim.pImageX1sg[index];
      y1sg=cSim.pImageY1sg[index];
      z1sg=cSim.pImageZ1sg[index];
      x2sg=cSim.pImageX2sg[index];
      y2sg=cSim.pImageY2sg[index];
      z2sg=cSim.pImageZ2sg[index];
      rsgx=cSim.pImageRsgX[index];
      rsgy=cSim.pImageRsgY[index];
      rsgz=cSim.pImageRsgZ[index];
      fpsg=cSim.pImageFPsg[index];
      ppsg=cSim.pImagePPsg[index];
    }
#ifdef use_DPFP
    long long int LJID          = cSim.pImageLJID[index];
#else
    unsigned int LJID           = cSim.pImageLJID[index];
#endif
    unsigned int cellID         = cSim.pImageCellID[index];
    int TIRegion                = cSim.pImageTIRegion[index];
    cSim.pImageX2[pos]          = x;
    cSim.pImageY2[pos]          = y;
    cSim.pImageZ2[pos]          = z;
    cSim.pImageAtom2[pos]       = atom;
    cSim.pImageAtomLookup[atom] = pos;
    PMEFloat2 xy                = {(PMEFloat)x, (PMEFloat)y};
    cSim.pAtomXYSaveSP[pos]     = xy;
    cSim.pAtomZSaveSP[pos]      = z;
    cSim.pImageVelX2[pos]       = vx;
    cSim.pImageVelY2[pos]       = vy;
    cSim.pImageVelZ2[pos]       = vz;
    cSim.pImageLVelX2[pos]      = lvx;
    cSim.pImageLVelY2[pos]      = lvy;
    cSim.pImageLVelZ2[pos]      = lvz;
    cSim.pImageCharge2[pos]     = q;
    cSim.pAtomChargeSP[pos]     = q;
    cSim.pImageMass2[pos]       = m;
    cSim.pImageInvMass2[pos]    = (m != (PMEDouble)0.0 ? (PMEDouble)1.0 / m : (PMEDouble)0.0);
    cSim.pImageSigEps2[pos]     = sigeps;
    cSim.pImageLJID2[pos]       = LJID;
    cSim.pImageCellID2[pos]     = cellID;
    cSim.pImageTIRegion2[pos]   = TIRegion;

    if(cSim.iphmd == 3) {
      cSim.pImageGrplist2[pos] = grplist;
      cSim.pImageQstate12[pos] = qstate1;
      cSim.pImageQstate22[pos] = qstate2;
      cSim.pImageVstate12[pos] = vstate1;
      cSim.pImageVstate22[pos] = vstate2;
      cSim.pImageCharge_phmd2[pos] = qphmd;
      cSim.pImageIndexPHMD2[pos] = indexPHMD;
    }
    if(cSim.isgld > 0) {
      cSim.pImageWsg2[pos] = wsg;
      cSim.pImageX0sg2[pos] = x0sg;
      cSim.pImageY0sg2[pos] = y0sg;
      cSim.pImageZ0sg2[pos] = z0sg;
      cSim.pImageX1sg2[pos] = x1sg;
      cSim.pImageY1sg2[pos] = y1sg;
      cSim.pImageZ1sg2[pos] = z1sg;
      cSim.pImageX2sg2[pos] = x2sg;
      cSim.pImageY2sg2[pos] = y2sg;
      cSim.pImageZ2sg2[pos] = z2sg;
      cSim.pImageRsgX2[pos] = rsgx;
      cSim.pImageRsgY2[pos] = rsgy;
      cSim.pImageRsgZ2[pos] = rsgz;
      cSim.pImageFPsg2[pos] = fpsg;
      cSim.pImagePPsg2[pos] = ppsg;
    }
    // Need to store the value of pos so we can find it later for linear atom updating
    cSim.pUpdateIndex[index]    = pos;
#ifdef use_DPFP
    PMEFloat2 qljid             = {q, __longlong_as_double(LJID)};
#else
    PMEFloat2 qljid             = {(PMEFloat)q, __uint_as_float(LJID)};
#endif
    cSim.pAtomChargeSPLJID[pos] = qljid;

    // Advance to next atom
    index = newindex;
    pos   = newpos;
  }
}

//---------------------------------------------------------------------------------------------
// kNLTINoLinearAtmRemapImage_kernel: same as kNLTIRemapImage_kernel, but for cases with no
//                                    linear scaling atoms
//
// Arguments:
//   pImageIndex: buffer containing pointers to all remapped atom indexes. Provides the
//                connection from the previous position in the atom list to the position in the
//                atom list following the radix sort.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(THREADS_PER_BLOCK, 1)
kNLTINoLinearAtmRemapImage_kernel(unsigned int* pImageIndex)
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  unsigned int index = 0;
  unsigned int newindex;

  if (pos < cSim.atoms) {
    index = pImageIndex[pos];
  }
  while (pos < cSim.atoms) {
    unsigned int newpos = pos + increment;
    if (newpos < cSim.atoms) {
      newindex = pImageIndex[newpos];
    }

    // Read new data
    unsigned int atom           = cSim.pImageAtom[index];
    double x                    = cSim.pImageX[index];
    double y                    = cSim.pImageY[index];
    double z                    = cSim.pImageZ[index];
    double vx                   = cSim.pImageVelX[index];
    double vy                   = cSim.pImageVelY[index];
    double vz                   = cSim.pImageVelZ[index];
    double lvx                  = cSim.pImageLVelX[index];
    double lvy                  = cSim.pImageLVelY[index];
    double lvz                  = cSim.pImageLVelZ[index];
    double q                    = cSim.pImageCharge[index];
    double m                    = cSim.pImageMass[index];
    PMEFloat2 sigeps            = cSim.pImageSigEps[index];
    int grplist, indexPHMD;
    PMEFloat qphmd;
    PMEFloat2 qstate1, qstate2, vstate1, vstate2;
    if (cSim.iphmd == 3) {
      grplist = cSim.pImageGrplist[index];
      qstate1 = cSim.pImageQstate1[index];
      qstate2 = cSim.pImageQstate2[index];
      vstate1 = cSim.pImageVstate1[index];
      vstate2 = cSim.pImageVstate2[index];
      qphmd = cSim.pImageCharge_phmd2[index];
      indexPHMD = cSim.pImageIndexPHMD[index];
    }
    int wsg;
    double x0sg,y0sg,z0sg,x1sg,y1sg,z1sg,x2sg,y2sg,z2sg,rsgx,rsgy,rsgz,fpsg,ppsg;
    if(cSim.isgld > 0) {
      wsg=cSim.pImageWsg[index];
      x0sg=cSim.pImageX0sg[index];
      y0sg=cSim.pImageY0sg[index];
      z0sg=cSim.pImageZ0sg[index];
      x1sg=cSim.pImageX1sg[index];
      y1sg=cSim.pImageY1sg[index];
      z1sg=cSim.pImageZ1sg[index];
      x2sg=cSim.pImageX2sg[index];
      y2sg=cSim.pImageY2sg[index];
      z2sg=cSim.pImageZ2sg[index];
      rsgx=cSim.pImageRsgX[index];
      rsgy=cSim.pImageRsgY[index];
      rsgz=cSim.pImageRsgZ[index];
      fpsg=cSim.pImageFPsg[index];
      ppsg=cSim.pImagePPsg[index];
    }
#ifdef use_DPFP
    long long int LJID          = cSim.pImageLJID[index];
#else
    unsigned int LJID           = cSim.pImageLJID[index];
#endif
    unsigned int cellID         = cSim.pImageCellID[index];
    int TIRegion                = cSim.pImageTIRegion[index];
    cSim.pImageX2[pos]          = x;
    cSim.pImageY2[pos]          = y;
    cSim.pImageZ2[pos]          = z;
    cSim.pImageAtom2[pos]       = atom;
    cSim.pImageAtomLookup[atom] = pos;
    PMEFloat2 xy                = {(PMEFloat)x, (PMEFloat)y};
    cSim.pAtomXYSaveSP[pos]     = xy;
    cSim.pAtomZSaveSP[pos]      = z;
    cSim.pImageVelX2[pos]       = vx;
    cSim.pImageVelY2[pos]       = vy;
    cSim.pImageVelZ2[pos]       = vz;
    cSim.pImageLVelX2[pos]      = lvx;
    cSim.pImageLVelY2[pos]      = lvy;
    cSim.pImageLVelZ2[pos]      = lvz;
    cSim.pImageCharge2[pos]     = q;
    cSim.pAtomChargeSP[pos]     = q;
    cSim.pImageMass2[pos]       = m;
    cSim.pImageInvMass2[pos]    = (m != (PMEDouble)0.0 ? (PMEDouble)1.0 / m : (PMEDouble)0.0);
    cSim.pImageSigEps2[pos]     = sigeps;
    cSim.pImageLJID2[pos]       = LJID;
    cSim.pImageCellID2[pos]     = cellID;
    cSim.pImageTIRegion2[pos]   = TIRegion;
    if(cSim.iphmd == 3) {
      cSim.pImageGrplist2[pos] = grplist;
      cSim.pImageQstate12[pos] = qstate1;
      cSim.pImageQstate22[pos] = qstate2;
      cSim.pImageVstate12[pos] = vstate1;
      cSim.pImageVstate22[pos] = vstate2;
      cSim.pImageCharge_phmd2[pos] = qphmd;
      cSim.pImageIndexPHMD2[pos] = indexPHMD;
    }
    if(cSim.isgld > 0) {
      cSim.pImageWsg2[pos] = wsg;
      cSim.pImageX0sg2[pos] = x0sg;
      cSim.pImageY0sg2[pos] = y0sg;
      cSim.pImageZ0sg2[pos] = z0sg;
      cSim.pImageX1sg2[pos] = x1sg;
      cSim.pImageY1sg2[pos] = y1sg;
      cSim.pImageZ1sg2[pos] = z1sg;
      cSim.pImageX2sg2[pos] = x2sg;
      cSim.pImageY2sg2[pos] = y2sg;
      cSim.pImageZ2sg2[pos] = z2sg;
      cSim.pImageRsgX2[pos] = rsgx;
      cSim.pImageRsgY2[pos] = rsgy;
      cSim.pImageRsgZ2[pos] = rsgz;
      cSim.pImageFPsg2[pos] = fpsg;
      cSim.pImagePPsg2[pos] = ppsg;
    }
#ifdef use_DPFP
    PMEFloat2 qljid             = {q, __longlong_as_double(LJID)};
#else
    PMEFloat2 qljid             = {(PMEFloat)q, __uint_as_float(LJID)};
#endif
    cSim.pAtomChargeSPLJID[pos] = qljid;

    // Advance to next atom
    index = newindex;
    pos   = newpos;
  }
}

//---------------------------------------------------------------------------------------------
// kNLRemapLinearAtmID_kernel: Updates the indexes of linear atoms to allow for force mixing or
//                             position/velocity syncing of said atoms between region 1 and
//                             region 2 following force/position/velocity updates
//
// Arguments:
//   pImageIndex: buffer containing pointers to all remapped atom indexes. Provides the
//                connection from the previous position in the atom list to the position in the
//                atom list following the radix sort.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(THREADS_PER_BLOCK, 1)
kNLRemapLinearAtmID_kernel(unsigned int* pImageIndex)
{
  unsigned int pos                = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int increment          = gridDim.x * blockDim.x;
  unsigned int linearAtmCnt       = cSim.TIlinearAtmCnt;
  unsigned int paddedLinearAtmCnt = cSim.TIPaddedLinearAtmCnt;
  while (pos < linearAtmCnt) {

    // Region 1
    cSim.pImageTILinearAtmID[pos] = cSim.pUpdateIndex[cSim.pImageTILinearAtmID[pos]];

    // Region 2
    cSim.pImageTILinearAtmID[pos+paddedLinearAtmCnt] =
      cSim.pUpdateIndex[cSim.pImageTILinearAtmID[pos+paddedLinearAtmCnt]];
    pos += increment;
  }
}

//---------------------------------------------------------------------------------------------
// kNLRemapImage: this is called on the heels of the radix sort and will launch the appropriate
//                kernel (or kernels) for making the map to the neighbor list atom ordering
//                from the sequential atom ordering.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kNLRemapImage(gpuContext gpu)
{
  if (gpu->sim.ti_mode == 0) {
    kNLRemapImage_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>(gpu->sim.pImageIndex);
  }
  else {
    if (gpu->sim.TIlinearAtmCnt > 0) {
      kNLTIRemapImage_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>(gpu->sim.pImageIndex);
      kNLRemapLinearAtmID_kernel<<<gpu->AFEExchangeBlocks,
                                   gpu->AFEExchangeThreadsPerBlock>>>(gpu->sim.pImageIndex);
    }
    else {
      kNLTINoLinearAtmRemapImage_kernel<<<gpu->blocks,
                                          gpu->threadsPerBlock>>>(gpu->sim.pImageIndex);
    }
  }
  LAUNCHERROR("kNLRemapImage");
  unsigned int *pUint;
  int* pInt;
  double* pDouble;
  PMEFloat2* pPMEFloat2;
  PMEFloat* pPMEFloat;

  // Remap constant memory pointers
  pUint                           = gpu->sim.pImageAtom;
  gpu->sim.pImageAtom             = gpu->sim.pImageAtom2;
  gpu->sim.pImageAtom2            = pUint;
  pDouble                         = gpu->sim.pImageX;
  gpu->sim.pImageX                = gpu->sim.pImageX2;
  {
    cudaError_t status;
    status = cudaDestroyTextureObject(gpu->sim.texImageX);
    RTERROR(status, "cudaDestroyTextureObject gpu->sim.texImageX failed");
    cudaResourceDesc resDesc;
    cudaTextureDesc texDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceTypeLinear;
    resDesc.res.linear.devPtr = gpu->sim.pImageX;
    resDesc.res.linear.sizeInBytes = gpu->sim.stride3 * sizeof(int2);
    resDesc.res.linear.desc.f = cudaChannelFormatKindSigned;
    resDesc.res.linear.desc.x = 32;
    resDesc.res.linear.desc.y = 32;
    resDesc.res.linear.desc.z = 0;
    resDesc.res.linear.desc.w = 0;
    memset(&texDesc, 0, sizeof(texDesc));
    texDesc.normalizedCoords = 0;
    texDesc.filterMode = cudaFilterModePoint;
    texDesc.addressMode[0] = cudaAddressModeClamp;
    texDesc.readMode = cudaReadModeElementType;
    status = cudaCreateTextureObject(&(gpu->sim.texImageX), &resDesc, &texDesc, NULL);
    RTERROR(status, "cudaCreateTextureObject gpu->sim.texImageX failed");
  }
  gpu->sim.pImageX2               = pDouble;
  pDouble                         = gpu->sim.pImageY;
  gpu->sim.pImageY                = gpu->sim.pImageY2;
  gpu->sim.pImageY2               = pDouble;
  pDouble                         = gpu->sim.pImageZ;
  gpu->sim.pImageZ                = gpu->sim.pImageZ2;
  gpu->sim.pImageZ2               = pDouble;
  pDouble                         = gpu->sim.pImageVelX;
  gpu->sim.pImageVelX             = gpu->sim.pImageVelX2;
  gpu->sim.pImageVelX2            = pDouble;
  pDouble                         = gpu->sim.pImageVelY;
  gpu->sim.pImageVelY             = gpu->sim.pImageVelY2;
  gpu->sim.pImageVelY2            = pDouble;
  pDouble                         = gpu->sim.pImageVelZ;
  gpu->sim.pImageVelZ             = gpu->sim.pImageVelZ2;
  gpu->sim.pImageVelZ2            = pDouble;
  pDouble                         = gpu->sim.pImageLVelX;
  gpu->sim.pImageLVelX            = gpu->sim.pImageLVelX2;
  gpu->sim.pImageLVelX2           = pDouble;
  pDouble                         = gpu->sim.pImageLVelY;
  gpu->sim.pImageLVelY            = gpu->sim.pImageLVelY2;
  gpu->sim.pImageLVelY2           = pDouble;
  pDouble                         = gpu->sim.pImageLVelZ;
  gpu->sim.pImageLVelZ            = gpu->sim.pImageLVelZ2;
  gpu->sim.pImageLVelZ2           = pDouble;
  pDouble                         = gpu->sim.pImageCharge;
  gpu->sim.pImageCharge           = gpu->sim.pImageCharge2;
  gpu->sim.pImageCharge2          = pDouble;
  pDouble                         = gpu->sim.pImageMass;
  gpu->sim.pImageMass             = gpu->sim.pImageMass2;
  gpu->sim.pImageMass2            = pDouble;
  pDouble                         = gpu->sim.pImageInvMass;
  gpu->sim.pImageInvMass          = gpu->sim.pImageInvMass2;
  gpu->sim.pImageInvMass2         = pDouble;
  pPMEFloat2                      = gpu->sim.pImageSigEps;
  gpu->sim.pImageSigEps           = gpu->sim.pImageSigEps2;
  gpu->sim.pImageSigEps2          = pPMEFloat2;
  pUint                           = gpu->sim.pImageLJID;
  gpu->sim.pImageLJID             = gpu->sim.pImageLJID2;
  gpu->sim.pImageLJID2            = pUint;
  pUint                           = gpu->sim.pImageCellID;
  gpu->sim.pImageCellID           = gpu->sim.pImageCellID2;
  gpu->sim.pImageCellID2          = pUint;
  pInt                            = gpu->sim.pImageTIRegion;
  gpu->sim.pImageTIRegion         = gpu->sim.pImageTIRegion2;
  gpu->sim.pImageTIRegion2        = pInt;
  if (gpu->iphmd == 3) {
    pInt                          = gpu->sim.pImageGrplist;
    gpu->sim.pImageGrplist        = gpu->sim.pImageGrplist2;
    gpu->sim.pImageGrplist2       = pInt;
    pInt                          = gpu->sim.pImageIndexPHMD;
    gpu->sim.pImageIndexPHMD      = gpu->sim.pImageIndexPHMD2;
    gpu->sim.pImageIndexPHMD2     = pInt;
    pPMEFloat2                    = gpu->sim.pImageQstate1;
    gpu->sim.pImageQstate1        = gpu->sim.pImageQstate12;
    gpu->sim.pImageQstate12       = pPMEFloat2;
    pPMEFloat2                    = gpu->sim.pImageQstate2;
    gpu->sim.pImageQstate2        = gpu->sim.pImageQstate22;
    gpu->sim.pImageQstate22       = pPMEFloat2;
    pPMEFloat2                    = gpu->sim.pImageVstate1;
    gpu->sim.pImageVstate1        = gpu->sim.pImageVstate12;
    gpu->sim.pImageVstate12       = pPMEFloat2;
    pPMEFloat2                    = gpu->sim.pImageVstate2;
    gpu->sim.pImageVstate2        = gpu->sim.pImageVstate22;
    gpu->sim.pImageVstate22       = pPMEFloat2;
    pPMEFloat                     = gpu->sim.pImageCharge_phmd;
    gpu->sim.pImageCharge_phmd    = gpu->sim.pImageCharge_phmd2;
    gpu->sim.pImageCharge_phmd2   = pPMEFloat;
  }
  if(gpu->sim.isgld > 0) {
    pInt                       = gpu->sim.pImageWsg;
    gpu->sim.pImageWsg            = gpu->sim.pImageWsg2;
    gpu->sim.pImageWsg2           = pInt;
    pDouble                       = gpu->sim.pImageX0sg;
    gpu->sim.pImageX0sg            = gpu->sim.pImageX0sg2;
    gpu->sim.pImageX0sg2           = pDouble;
    pDouble                       = gpu->sim.pImageY0sg;
    gpu->sim.pImageY0sg            = gpu->sim.pImageY0sg2;
    gpu->sim.pImageY0sg2           = pDouble;
    pDouble                       = gpu->sim.pImageZ0sg;
    gpu->sim.pImageZ0sg            = gpu->sim.pImageZ0sg2;
    gpu->sim.pImageZ0sg2           = pDouble;
    pDouble                       = gpu->sim.pImageX1sg;
    gpu->sim.pImageX1sg            = gpu->sim.pImageX1sg2;
    gpu->sim.pImageX1sg2           = pDouble;
    pDouble                       = gpu->sim.pImageY1sg;
    gpu->sim.pImageY1sg            = gpu->sim.pImageY1sg2;
    gpu->sim.pImageY1sg2           = pDouble;
    pDouble                       = gpu->sim.pImageZ1sg;
    gpu->sim.pImageZ1sg            = gpu->sim.pImageZ1sg2;
    gpu->sim.pImageZ1sg2           = pDouble;
    pDouble                       = gpu->sim.pImageX2sg;
    gpu->sim.pImageX2sg            = gpu->sim.pImageX2sg2;
    gpu->sim.pImageX2sg2           = pDouble;
    pDouble                       = gpu->sim.pImageY2sg;
    gpu->sim.pImageY2sg            = gpu->sim.pImageY2sg2;
    gpu->sim.pImageY2sg2           = pDouble;
    pDouble                       = gpu->sim.pImageZ2sg;
    gpu->sim.pImageZ2sg            = gpu->sim.pImageZ2sg2;
    gpu->sim.pImageZ2sg2           = pDouble;
    pDouble                       = gpu->sim.pImageRsgX;
    gpu->sim.pImageRsgX            = gpu->sim.pImageRsgX2;
    gpu->sim.pImageRsgX2           = pDouble;
    pDouble                       = gpu->sim.pImageRsgY;
    gpu->sim.pImageRsgY            = gpu->sim.pImageRsgY2;
    gpu->sim.pImageRsgY2           = pDouble;
    pDouble                       = gpu->sim.pImageRsgZ;
    gpu->sim.pImageRsgZ            = gpu->sim.pImageRsgZ2;
    gpu->sim.pImageRsgZ2           = pDouble;
    pDouble                       = gpu->sim.pImageFPsg;
    gpu->sim.pImageFPsg            = gpu->sim.pImageFPsg2;
    gpu->sim.pImageFPsg2           = pDouble;
    pDouble                       = gpu->sim.pImagePPsg;
    gpu->sim.pImagePPsg            = gpu->sim.pImagePPsg2;
    gpu->sim.pImagePPsg2           = pDouble;
    }
}

//---------------------------------------------------------------------------------------------
// kNLRemapLocalInteractions_kernel: kernel for remapping local (read: bonded) interactions,
//                                   delegating them for processing within the groupings of the
//                                   hash table.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(THREADS_PER_BLOCK, 1)
kNLRemapLocalInteractions_kernel()
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;

  // Remapping for bond, angle, dihedral, CHARMM, and ordinary NMR
  // interactions is covered by kNLRemapBondWorkUnits below.

  // Remap Shake constraints
  while (pos < cSim.shakeOffset) {
    if (pos < cSim.shakeConstraints) {
      int4 atom = cSim.pShakeID[pos];
      if (cSim.tishake != 2) {
        if (atom.x != -1) {
          atom.x    = cSim.pImageAtomLookup[atom.x];
        }
        if (atom.y != -1) {
          atom.y = cSim.pImageAtomLookup[atom.y];
        }
        if (atom.z != -1) {
          atom.z  = cSim.pImageAtomLookup[atom.z];
        }
        if (atom.w != -1) {
          atom.w  = cSim.pImageAtomLookup[atom.w];
        }
      }
      else {
        if(atom.x < -1)
        {
          atom.x  = abs(atom.x)-2;
          atom.x  = -1*(cSim.pImageAtomLookup[atom.x]+1)-1;
        }
        else
        {
          atom.x  = cSim.pImageAtomLookup[atom.x];
        }
        if(atom.y < -1)
        {
          atom.y  = abs(atom.y)-2;
          atom.y  = -1*(cSim.pImageAtomLookup[atom.y]+1)-1;
        }
        else
        {
          atom.y  = cSim.pImageAtomLookup[atom.y];
        }
        if (atom.z != -1) {
          if(atom.z < -1)
          {
            atom.z  = abs(atom.z)-2;
            atom.z  = -1*(cSim.pImageAtomLookup[atom.z]+1)-1;
          }
          else
          {
            atom.z  = cSim.pImageAtomLookup[atom.z];
          }
        }
        if (atom.w != -1) {
          if(atom.w < -1)
          {
            atom.w  = abs(atom.w)-2;
            atom.w  = -1*(cSim.pImageAtomLookup[atom.w]+1)-1;
          }
          else
          {
            atom.w  = cSim.pImageAtomLookup[atom.w];
          }
        }
      }
      cSim.pImageShakeID[pos] = atom;
    }
    pos += blockDim.x * gridDim.x;
  }
  pos -= cSim.shakeOffset;
  while (pos < cSim.fastShakeOffset) {
    if (pos < cSim.fastShakeConstraints) {
      int4 atom = cSim.pFastShakeID[pos];
      atom.x    = cSim.pImageAtomLookup[atom.x];
      atom.y    = cSim.pImageAtomLookup[atom.y];
      atom.z    = cSim.pImageAtomLookup[atom.z];
      cSim.pImageFastShakeID[pos] = atom;
    }
    pos += blockDim.x * gridDim.x;
  }
  pos -= cSim.fastShakeOffset;

  if(cSim.tishake == 2) {
    while (pos < cSim.slowShakeOffset) {
      if (pos < cSim.slowShakeConstraints) {
        int atom1  = cSim.pSlowShakeID1[pos];
        int4 atom2 = cSim.pSlowShakeID2[pos];
        atom1      = cSim.pImageAtomLookup[atom1];
        if(atom2.x < -1)
        {
          atom2.x  = abs(atom2.x)-2;
          atom2.x  = -1*(cSim.pImageAtomLookup[atom2.x]+1)-1;
        }
        else
        {
          atom2.x  = cSim.pImageAtomLookup[atom2.x];
        }
        if(atom2.y < -1)
        {
          atom2.y  = abs(atom2.y)-2;
          atom2.y  = -1*(cSim.pImageAtomLookup[atom2.y]+1)-1;
        }
        else
        {
          atom2.y  = cSim.pImageAtomLookup[atom2.y];
        }
        if(atom2.z < -1)
        {
          atom2.z  = abs(atom2.z)-2;
          atom2.z  = -1*(cSim.pImageAtomLookup[atom2.z]+1)-1;
        }
        else
        {
          atom2.z  = cSim.pImageAtomLookup[atom2.z];
        }
        if(atom2.w < -1)
        {
          atom2.w  = abs(atom2.w)-2;
          atom2.w  = -1*(cSim.pImageAtomLookup[atom2.w]+1)-1;
        }
        else
        {
          atom2.w  = cSim.pImageAtomLookup[atom2.w];
        }
        cSim.pImageSlowShakeID1[pos] = atom1;
        cSim.pImageSlowShakeID2[pos] = atom2;
      }
      pos += blockDim.x * gridDim.x;
    }
  }
  else {
    while (pos < cSim.slowShakeOffset) {
      if (pos < cSim.slowShakeConstraints) {
        int atom1  = cSim.pSlowShakeID1[pos];
        int4 atom2 = cSim.pSlowShakeID2[pos];
        atom1      = cSim.pImageAtomLookup[atom1];
        atom2.x    = cSim.pImageAtomLookup[atom2.x];
        atom2.y    = cSim.pImageAtomLookup[atom2.y];
        atom2.z    = cSim.pImageAtomLookup[atom2.z];
        atom2.w    = cSim.pImageAtomLookup[atom2.w];
        cSim.pImageSlowShakeID1[pos] = atom1;
        cSim.pImageSlowShakeID2[pos] = atom2;
      }
      pos += blockDim.x * gridDim.x;
    }
  }
  pos -= cSim.slowShakeOffset;

  if(cSim.tishake == 2) {
    while ( pos < cSim.slowTIShakeOffset) {
      if(pos < cSim.slowTIShakeConstraints) {
        int atom1  = cSim.pSlowTIShakeID1[pos];
        int4 atom2 = cSim.pSlowTIShakeID2[pos];
        int4 atom3 = cSim.pSlowTIShakeID3[pos];
        atom1      = cSim.pImageAtomLookup[atom1];
        if(atom2.x < -1)
        {
          atom2.x  = abs(atom2.x)-2;
          atom2.x  = -1*(cSim.pImageAtomLookup[atom2.x]+1)-1;
        }
        else
        {
          atom2.x  = cSim.pImageAtomLookup[atom2.x];
        }
        if(atom2.y < -1)
        {
          atom2.y  = abs(atom2.y)-2;
          atom2.y  = -1*(cSim.pImageAtomLookup[atom2.y]+1)-1;
        }
        else
        {
          atom2.y  = cSim.pImageAtomLookup[atom2.y];
        }
        if(atom2.z < -1)
        {
          atom2.z  = abs(atom2.z)-2;
          atom2.z  = -1*(cSim.pImageAtomLookup[atom2.z]+1)-1;
        }
        else
        {
          atom2.z  = cSim.pImageAtomLookup[atom2.z];
        }
        if(atom2.w < -1)
        {
          atom2.w  = abs(atom2.w)-2;
          atom2.w  = -1*(cSim.pImageAtomLookup[atom2.w]+1)-1;
        }
        else
        {
          atom2.w  = cSim.pImageAtomLookup[atom2.w];
        }
        if(atom3.x < -1)
        {
          atom3.x  = abs(atom3.x)-2;
          atom3.x  = -1*(cSim.pImageAtomLookup[atom3.x]+1)-1;
        }
        else
        {
          atom3.x  = cSim.pImageAtomLookup[atom3.x];
        }
        if (atom3.y != -1) {
          if(atom3.y < -1)
          {
            atom3.y  = abs(atom3.y)-2;
            atom3.y  = -1*(cSim.pImageAtomLookup[atom3.y]+1)-1;
          }
          else
          {
            atom3.y  = cSim.pImageAtomLookup[atom3.y];
          }
        }
        if (atom3.z != -1) {
          if(atom3.z < -1)
          {
            atom3.z  = abs(atom3.z)-2;
            atom3.z  = -1*(cSim.pImageAtomLookup[atom3.z]+1)-1;
          }
          else
          {
            atom3.z  = cSim.pImageAtomLookup[atom3.z];
          }
        }
        if (atom3.w != -1) {
          if(atom3.w < -1)
          {
            atom3.w  = abs(atom3.w)-2;
            atom3.w  = -1*(cSim.pImageAtomLookup[atom3.w]+1)-1;
          }
          else
          {
            atom3.w  = cSim.pImageAtomLookup[atom3.w];
          }
        }
        cSim.pImageSlowTIShakeID1[pos] = atom1;
        cSim.pImageSlowTIShakeID2[pos] = atom2;
        cSim.pImageSlowTIShakeID3[pos] = atom3;
      }
      pos +=  blockDim.x * gridDim.x;
    }
    pos -= cSim.slowTIShakeOffset;
  }

  // Solute atoms already padded to warp width
  while (pos < cSim.soluteAtoms) {
    int atom = cSim.pSoluteAtomID[pos];
    if (atom != -1) {
      atom = cSim.pImageAtomLookup[atom];
    }
    cSim.pImageSoluteAtomID[pos] = atom;
    pos += blockDim.x * gridDim.x;
  }
  pos -= cSim.soluteAtoms;

  while (pos < cSim.solventMoleculeStride) {
    if (pos < cSim.solventMolecules) {
      int4 atom = cSim.pSolventAtomID[pos];
      atom.x    = cSim.pImageAtomLookup[atom.x];
      if (atom.y != -1) {
        atom.y  = cSim.pImageAtomLookup[atom.y];
      }
      if (atom.z != -1) {
        atom.z  = cSim.pImageAtomLookup[atom.z];
      }
      if (atom.w != -1) {
        atom.w  = cSim.pImageAtomLookup[atom.w];
      }
      cSim.pImageSolventAtomID[pos] = atom;
    }
    pos += blockDim.x * gridDim.x;
  }
  pos -= cSim.solventMoleculeStride;

  while (pos < cSim.EP11Offset) {
    if (pos < cSim.EP11s) {
      int4 frame = cSim.pExtraPoint11Frame[pos];
      int index  = cSim.pExtraPoint11Index[pos];
      frame.x    = cSim.pImageAtomLookup[frame.x];
      frame.y    = cSim.pImageAtomLookup[frame.y];
      frame.z    = cSim.pImageAtomLookup[frame.z];
      frame.w    = cSim.pImageAtomLookup[frame.w];
      index      = cSim.pImageAtomLookup[index];
      cSim.pImageExtraPoint11Frame[pos] = frame;
      cSim.pImageExtraPoint11Index[pos] = index;
    }
    pos += blockDim.x * gridDim.x;
  }
  while (pos < cSim.EP12Offset) {
    pos -= cSim.EP11Offset;
    if (pos < cSim.EP12s) {
      int4 frame = cSim.pExtraPoint12Frame[pos];
      int index  = cSim.pExtraPoint12Index[pos];
      frame.x    = cSim.pImageAtomLookup[frame.x];
      frame.y    = cSim.pImageAtomLookup[frame.y];
      frame.z    = cSim.pImageAtomLookup[frame.z];
      frame.w    = cSim.pImageAtomLookup[frame.w];
      index      = cSim.pImageAtomLookup[index];
      cSim.pImageExtraPoint12Frame[pos] = frame;
      cSim.pImageExtraPoint12Index[pos] = index;
    }
    pos += cSim.EP11Offset + blockDim.x * gridDim.x;
  }
  while (pos < cSim.EP21Offset) {
    pos -= cSim.EP12Offset;
    if (pos < cSim.EP21s) {
      int4 frame = cSim.pExtraPoint21Frame[pos];
      int2 index = cSim.pExtraPoint21Index[pos];
      frame.x    = cSim.pImageAtomLookup[frame.x];
      frame.y    = cSim.pImageAtomLookup[frame.y];
      frame.z    = cSim.pImageAtomLookup[frame.z];
      frame.w    = cSim.pImageAtomLookup[frame.w];
      index.x    = cSim.pImageAtomLookup[index.x];
      index.y    = cSim.pImageAtomLookup[index.y];
      cSim.pImageExtraPoint21Frame[pos] = frame;
      cSim.pImageExtraPoint21Index[pos] = index;
    }
    pos += cSim.EP12Offset + blockDim.x * gridDim.x;
  }
  while (pos < cSim.EP22Offset) {
    pos -= cSim.EP21Offset;
    if (pos < cSim.EP22s) {
      int4 frame = cSim.pExtraPoint22Frame[pos];
      int2 index = cSim.pExtraPoint22Index[pos];
      frame.x    = cSim.pImageAtomLookup[frame.x];
      frame.y    = cSim.pImageAtomLookup[frame.y];
      frame.z    = cSim.pImageAtomLookup[frame.z];
      frame.w    = cSim.pImageAtomLookup[frame.w];
      index.x    = cSim.pImageAtomLookup[index.x];
      index.y    = cSim.pImageAtomLookup[index.y];
      cSim.pImageExtraPoint22Frame[pos] = frame;
      cSim.pImageExtraPoint22Index[pos] = index;
    }
    pos += cSim.EP21Offset + blockDim.x * gridDim.x;
  }
  while (pos < cSim.EPCustomOffset) {
    pos -= cSim.EP22Offset;
    if (pos < cSim.EPCustomCount) {
      int4 frame = cSim.pExtraPointCustomFrame[pos];
      int  index = cSim.pExtraPointCustomIndex[pos];
      frame.x    = ((cSim.pImageAtomLookup[(frame.x & 0xffffff)]) | (frame.x & 0xff000000));
      frame.y    = cSim.pImageAtomLookup[frame.y];
      if (frame.z >= 0) {
        frame.z    = cSim.pImageAtomLookup[frame.z];
      }
      if (frame.w >= 0) {
        frame.w    = cSim.pImageAtomLookup[frame.w];
      }
      index      = cSim.pImageAtomLookup[index];
      cSim.pImageEPCustomFrame[pos] = frame;
      cSim.pImageEPCustomIndex[pos] = index;
    }
    pos += cSim.EP22Offset + blockDim.x * gridDim.x;
  }
  pos -= cSim.EPCustomOffset;

  // COM to COM distance restraints require fetching the Iatoms for each atom and re-sorting
  // into groups of atoms for each COM.  New COM groupings are placed in
  // pImageNMRCOMDistanceCOM and the indexing for each groups is placed into
  // pImageNMRCOMDistanceCOMGrp.  pImageNMRCOMDistanceCOM is of size distances * maxgrp to
  // avoid race conditions.  pImageNMRCOMDistanceCOMGrp is of size distances and indexes
  // pImageNMRCOMDistanceCOM for each COM distance restraint.
  while (pos < cSim.NMRCOMDistanceOffset) {
    pos -= cSim.NMRTorsionOffset;
    if (pos < cSim.NMRCOMDistances) {
      int tilespc = pos * cSim.NMRMaxgrp;
      int2 atom   = cSim.pNMRCOMDistanceID[pos];
      if (atom.x > 0)
        atom.x = cSim.pImageAtomLookup[atom.x];
      else {
        int2 nmrcom_range = cSim.pNMRCOMDistanceCOMGrp[pos * 2];
        int grpnum        = 0;
        for (int ip = nmrcom_range.x; ip < nmrcom_range.y; ip++) {
          for (int ip1 = cSim.pNMRCOMDistanceCOM[ip].x;
               ip1 < cSim.pNMRCOMDistanceCOM[ip].y + 1; ip1++) {
            atom.x = cSim.pImageAtomLookup[ip1];
            if (grpnum == 0) {
              cSim.pImageNMRCOMDistanceCOM[grpnum + tilespc].x = atom.x;
              cSim.pImageNMRCOMDistanceCOM[grpnum + tilespc].y = atom.x;
              cSim.pImageNMRCOMDistanceCOMGrp[pos * 2].x = grpnum + tilespc;
              grpnum = grpnum + 1;
              cSim.pImageNMRCOMDistanceCOMGrp[pos * 2].y = grpnum + tilespc;
            }
            else {
              if (cSim.pImageNMRCOMDistanceCOM[(grpnum - 1) + tilespc].y == atom.x - 1) {
                cSim.pImageNMRCOMDistanceCOM[grpnum - 1 + tilespc].y = atom.x;
              }
              else {
                cSim.pImageNMRCOMDistanceCOM[grpnum + tilespc].x = atom.x;
                cSim.pImageNMRCOMDistanceCOM[grpnum + tilespc].y = atom.x;
                grpnum = grpnum + 1;
                cSim.pImageNMRCOMDistanceCOMGrp[pos * 2].y = grpnum + tilespc;
              }
            }
          }
        }
        atom.x = cSim.pNMRCOMDistanceID[pos].x;
      }
      if (atom.y > 0) {
        atom.y = cSim.pImageAtomLookup[atom.y];
      }
      else {
        int2 nmrcom_range = cSim.pNMRCOMDistanceCOMGrp[pos * 2 + 1];
        int grpnum;
        if (atom.x > 0) {
          grpnum = 0;
        }
        else {
           grpnum = cSim.pImageNMRCOMDistanceCOMGrp[pos * 2].y - tilespc;
        }
        int grpnum2 = 0;
        for (int ip = nmrcom_range.x; ip < nmrcom_range.y; ip++) {
          for (int ip2 = cSim.pNMRCOMDistanceCOM[ip].x;
               ip2 < cSim.pNMRCOMDistanceCOM[ip].y + 1; ip2++) {
            atom.y = cSim.pImageAtomLookup[ip2];
            if (grpnum2 == 0) {
              grpnum2 = grpnum;
              cSim.pImageNMRCOMDistanceCOM[grpnum2 + tilespc].x = atom.y;
              cSim.pImageNMRCOMDistanceCOM[grpnum2 + tilespc].y = atom.y;
              cSim.pImageNMRCOMDistanceCOMGrp[pos*2 + 1].x = grpnum2 + tilespc;
              grpnum2 += 1;
              cSim.pImageNMRCOMDistanceCOMGrp[pos*2 + 1].y = grpnum2 + tilespc;
            }
            else {
              if (cSim.pImageNMRCOMDistanceCOM[(grpnum2 - 1) + tilespc].y == atom.y - 1) {
                cSim.pImageNMRCOMDistanceCOM[grpnum2 - 1 + tilespc].y = atom.y;
              }
              else {
                cSim.pImageNMRCOMDistanceCOM[grpnum2 + tilespc].x = atom.y;
                cSim.pImageNMRCOMDistanceCOM[grpnum2 + tilespc].y = atom.y;
                grpnum2 += 1;
                cSim.pImageNMRCOMDistanceCOMGrp[pos*2 + 1].y = grpnum2 + tilespc;
              }
            }
          }
        }
        atom.y = cSim.pNMRCOMDistanceID[pos].y;
      }
      cSim.pImageNMRCOMDistanceID[pos] = atom;
    }
    pos += cSim.NMRTorsionOffset + (blockDim.x * gridDim.x);
  }

  // r6av to r6av distance restraints require fetching the Iatoms for each atom and re-sorting
  // into groups of atoms for each r6av.  New r6av groupings are placed in
  // pImageNMRr6avDistancer6av and the indexing for the groups are placed into
  // pImageNMRr6avDistancer6avGrp.  pImageNMRr6avDistancer6av is of size distances * maxgrp to
  // avoid race conditions.  pImageNMRr6avDistancer6avGrp is of size distances and indexes
  // pImageNMRr6avDistancer6av for each r6av distance restraint.
  while (pos < cSim.NMRr6avDistanceOffset) {
    pos -= cSim.NMRCOMDistanceOffset;
    if (pos < cSim.NMRr6avDistances) {
      int tilespc = pos * cSim.NMRMaxgrp;
      int2 atom   = cSim.pNMRr6avDistanceID[pos];
      if (atom.x > 0) {
        atom.x = cSim.pImageAtomLookup[atom.x];
      }
      else {
        int2 nmrcom_range = cSim.pNMRr6avDistancer6avGrp[pos * 2];
        int grpnum        = 0;
        for (int ip = nmrcom_range.x; ip < nmrcom_range.y; ip++) {
          for (int ip1 = cSim.pNMRr6avDistancer6av[ip].x;
               ip1 < cSim.pNMRr6avDistancer6av[ip].y + 1; ip1++) {
            atom.x = cSim.pImageAtomLookup[ip1];
            if (grpnum == 0) {
              cSim.pImageNMRr6avDistancer6av[grpnum + tilespc].x = atom.x;
              cSim.pImageNMRr6avDistancer6av[grpnum + tilespc].y = atom.x;
              cSim.pImageNMRr6avDistancer6avGrp[pos * 2].x = grpnum + tilespc;
              grpnum += 1;
              cSim.pImageNMRr6avDistancer6avGrp[pos * 2].y = grpnum + tilespc;
            }
            else {
              if (cSim.pImageNMRr6avDistancer6av[(grpnum - 1) + tilespc].y ==
                  atom.x - 1) {
                cSim.pImageNMRr6avDistancer6av[grpnum - 1 + tilespc].y = atom.x;
              }
              else {
                cSim.pImageNMRr6avDistancer6av[grpnum + tilespc].x = atom.x;
                cSim.pImageNMRr6avDistancer6av[grpnum + tilespc].y = atom.x;
                grpnum += 1;
                cSim.pImageNMRr6avDistancer6avGrp[pos * 2].y = grpnum + tilespc;
              }
            }
          }
        }
        atom.x = cSim.pNMRr6avDistanceID[pos].x;
      }
      if (atom.y > 0) {
        atom.y = cSim.pImageAtomLookup[atom.y];
      }
      else {
        int2 nmrcom_range = cSim.pNMRr6avDistancer6avGrp[pos*2 + 1];
        int grpnum;
        if (atom.x > 0) {
          grpnum = 0;
        }
        else {
          grpnum = cSim.pImageNMRr6avDistancer6avGrp[pos * 2].y - tilespc;
        }
        int grpnum2 = 0;
        for (int ip = nmrcom_range.x; ip < nmrcom_range.y; ip++) {
          for (int ip2 = cSim.pNMRr6avDistancer6av[ip].x;
               ip2 < cSim.pNMRr6avDistancer6av[ip].y + 1; ip2++) {
            atom.y = cSim.pImageAtomLookup[ip2];
            if (grpnum2 == 0) {
              grpnum2 = grpnum;
              cSim.pImageNMRr6avDistancer6av[grpnum2 + tilespc].x = atom.y;
              cSim.pImageNMRr6avDistancer6av[grpnum2 + tilespc].y = atom.y;
              cSim.pImageNMRr6avDistancer6avGrp[pos*2 + 1].x = grpnum2 + tilespc;
              grpnum2 = grpnum2 + 1;
              cSim.pImageNMRr6avDistancer6avGrp[pos*2 + 1].y = grpnum2 + tilespc;
            }
            else {
              if (cSim.pImageNMRr6avDistancer6av[(grpnum2 - 1) + tilespc].y ==
                  atom.y - 1) {
                cSim.pImageNMRr6avDistancer6av[grpnum2 - 1 + tilespc].y = atom.y;
              }
              else {
                cSim.pImageNMRr6avDistancer6av[grpnum2 + tilespc].x = atom.y;
                cSim.pImageNMRr6avDistancer6av[grpnum2 + tilespc].y = atom.y;
                grpnum2 = grpnum2 + 1;
                cSim.pImageNMRr6avDistancer6avGrp[pos*2 + 1].y = grpnum2 + tilespc;
              }
            }
          }
        }
        atom.y = cSim.pNMRr6avDistanceID[pos].y;
      }
      cSim.pImageNMRr6avDistanceID[pos] = atom;
    }
    pos += cSim.NMRCOMDistanceOffset + blockDim.x * gridDim.x;
  }

  // pImageNMRCOMAngleCOM and the indexing for each groups is placed into
  // pImageNMRCOMAngleCOMGrp.  pImageNMRCOMAngleCOM is of size angles * maxgrp to
  // avoid race conditions.  pImageNMRCOMAngleCOMGrp is of size angles and indexes
  // pImageNMRCOMAngleCOM for each COM angle restraint.
  while (pos < cSim.NMRCOMAngleOffset) {
    pos -= cSim.NMRr6avDistanceOffset;
    if (pos < cSim.NMRCOMAngles) {
      int tilespc = pos * cSim.NMRMaxgrp;
      int2 atom1   = cSim.pNMRCOMAngleID1[pos];
      int  atom2   = cSim.pNMRCOMAngleID2[pos];
      if (atom1.x > 0)
        atom1.x = cSim.pImageAtomLookup[atom1.x];
      else {
        int2 nmrcom_range = cSim.pNMRCOMAngleCOMGrp[pos * 3];
        int grpnum        = 0;
        for (int ip = nmrcom_range.x; ip < nmrcom_range.y; ip++) {
          for (int ip1 = cSim.pNMRCOMAngleCOM[ip].x;
               ip1 < cSim.pNMRCOMAngleCOM[ip].y + 1; ip1++) {
            atom1.x = cSim.pImageAtomLookup[ip1];
            if (grpnum == 0) {
              cSim.pImageNMRCOMAngleCOM[grpnum + tilespc].x = atom1.x;
              cSim.pImageNMRCOMAngleCOM[grpnum + tilespc].y = atom1.x;
              cSim.pImageNMRCOMAngleCOMGrp[pos * 3].x = grpnum + tilespc;
              grpnum = grpnum + 1;
              cSim.pImageNMRCOMAngleCOMGrp[pos * 3].y = grpnum + tilespc;
            }
            else {
              if (cSim.pImageNMRCOMAngleCOM[(grpnum - 1) + tilespc].y == atom1.x - 1) {
                cSim.pImageNMRCOMAngleCOM[grpnum - 1 + tilespc].y = atom1.x;
              }
              else {
                cSim.pImageNMRCOMAngleCOM[grpnum + tilespc].x = atom1.x;
                cSim.pImageNMRCOMAngleCOM[grpnum + tilespc].y = atom1.x;
                grpnum = grpnum + 1;
                cSim.pImageNMRCOMAngleCOMGrp[pos * 3].y = grpnum + tilespc;
              }
            }
          }
        }
        atom1.x = cSim.pNMRCOMAngleID1[pos].x;
      }
      // start COM 2
      if (atom1.y > 0) {
        atom1.y = cSim.pImageAtomLookup[atom1.y];
      }
      else {
        int2 nmrcom_range = cSim.pNMRCOMAngleCOMGrp[pos * 3 + 1];
        int grpnum;
        if (atom1.x > 0) {
          grpnum = 0;
        }
        else {
           grpnum = cSim.pImageNMRCOMAngleCOMGrp[pos * 3].y - tilespc;
        }
        int grpnum2 = 0;
        for (int ip = nmrcom_range.x; ip < nmrcom_range.y; ip++) {
          for (int ip2 = cSim.pNMRCOMAngleCOM[ip].x;
               ip2 < cSim.pNMRCOMAngleCOM[ip].y + 1; ip2++) {
            atom1.y = cSim.pImageAtomLookup[ip2];
            if (grpnum2 == 0) {
              grpnum2 = grpnum;
              cSim.pImageNMRCOMAngleCOM[grpnum2 + tilespc].x = atom1.y;
              cSim.pImageNMRCOMAngleCOM[grpnum2 + tilespc].y = atom1.y;
              cSim.pImageNMRCOMAngleCOMGrp[pos*3 + 1].x = grpnum2 + tilespc;
              grpnum2 += 1;
              cSim.pImageNMRCOMAngleCOMGrp[pos*3 + 1].y = grpnum2 + tilespc;
            }
            else {
              if (cSim.pImageNMRCOMAngleCOM[(grpnum2 - 1) + tilespc].y == atom1.y - 1) {
                cSim.pImageNMRCOMAngleCOM[grpnum2 - 1 + tilespc].y = atom1.y;
              }
              else {
                cSim.pImageNMRCOMAngleCOM[grpnum2 + tilespc].x = atom1.y;
                cSim.pImageNMRCOMAngleCOM[grpnum2 + tilespc].y = atom1.y;
                grpnum2 += 1;
                cSim.pImageNMRCOMAngleCOMGrp[pos*3 + 1].y = grpnum2 + tilespc;
              }
            }
          }
        }
        atom1.y = cSim.pNMRCOMAngleID1[pos].y;
      }
      // start COM 3
      if (atom2 > 0) {
        atom2 = cSim.pImageAtomLookup[atom2];
      }
      else {
        int2 nmrcom_range = cSim.pNMRCOMAngleCOMGrp[pos * 3 + 2];
        int grpnum;
        int grpnum2;
        if (atom1.x > 0) {
          grpnum = 0;
        }
        else {
           grpnum = cSim.pImageNMRCOMAngleCOMGrp[pos * 3].y - tilespc;
        }
        if (atom1.y > 0) {
          grpnum2 = 0;
        }
        else {
           grpnum2 = cSim.pImageNMRCOMAngleCOMGrp[pos * 3 + 1].y - tilespc;
        }
        grpnum2 += grpnum;
        int grpnum3 = 0;
        for (int ip = nmrcom_range.x; ip < nmrcom_range.y; ip++) {
          for (int ip3 = cSim.pNMRCOMAngleCOM[ip].x;
               ip3 < cSim.pNMRCOMAngleCOM[ip].y + 1; ip3++) {
            atom2 = cSim.pImageAtomLookup[ip3];
            if (grpnum3 == 0) {
              grpnum3 = grpnum2;
              cSim.pImageNMRCOMAngleCOM[grpnum3 + tilespc].x = atom2;
              cSim.pImageNMRCOMAngleCOM[grpnum3 + tilespc].y = atom2;
              cSim.pImageNMRCOMAngleCOMGrp[pos*3 + 2].x = grpnum3 + tilespc;
              grpnum3 += 1;
              cSim.pImageNMRCOMAngleCOMGrp[pos*3 + 2].y = grpnum3 + tilespc;
            }
            else {
              if (cSim.pImageNMRCOMAngleCOM[(grpnum3 - 1) + tilespc].y == (atom2 - 1)) {
                cSim.pImageNMRCOMAngleCOM[grpnum3 - 1 + tilespc].y = atom2;
              }
              else {
                cSim.pImageNMRCOMAngleCOM[grpnum3 + tilespc].x = atom2;
                cSim.pImageNMRCOMAngleCOM[grpnum3 + tilespc].y = atom2;
                grpnum3 += 1;
                cSim.pImageNMRCOMAngleCOMGrp[pos*3 + 2].y = grpnum3 + tilespc;
              }
            }
          }
        }
        atom2 = cSim.pNMRCOMAngleID2[pos];
      }
      cSim.pImageNMRCOMAngleID1[pos] = atom1;
      cSim.pImageNMRCOMAngleID2[pos] = atom2;
    }
    pos += cSim.NMRr6avDistanceOffset + (blockDim.x * gridDim.x);
  }

  // pImageNMRCOMTorsionCOM and the indexing for each groups is placed into
  // pImageNMRCOMTorsionCOMGrp.  pImageNMRCOMTorsionCOM is of size torsions * maxgrp to
  // avoid race conditions.  pImageNMRCOMTorsionCOMGrp is of size torsions and indexes
  // pImageNMRCOMTorsionCOM for each COM torsion restraint.
  while (pos < cSim.NMRCOMTorsionOffset) {
    pos -= cSim.NMRCOMAngleOffset;
    if (pos < cSim.NMRCOMTorsions) {
      int tilespc = pos * cSim.NMRMaxgrp;
      int4 atom   = cSim.pNMRCOMTorsionID1[pos];
      // COM 1
      if (atom.x > 0)
        atom.x = cSim.pImageAtomLookup[atom.x];
      else {
        int2 nmrcom_range = cSim.pNMRCOMTorsionCOMGrp[pos * 4];
        int grpnum        = 0;
        for (int ip = nmrcom_range.x; ip < nmrcom_range.y; ip++) {
          for (int ip1 = cSim.pNMRCOMTorsionCOM[ip].x;
               ip1 < cSim.pNMRCOMTorsionCOM[ip].y + 1; ip1++) {
            atom.x = cSim.pImageAtomLookup[ip1];
            if (grpnum == 0) {
              cSim.pImageNMRCOMTorsionCOM[grpnum + tilespc].x = atom.x;
              cSim.pImageNMRCOMTorsionCOM[grpnum + tilespc].y = atom.x;
              cSim.pImageNMRCOMTorsionCOMGrp[pos * 4].x = grpnum + tilespc;
              grpnum += 1;
              cSim.pImageNMRCOMTorsionCOMGrp[pos * 4].y = grpnum + tilespc;
            }
            else {
              if (cSim.pImageNMRCOMTorsionCOM[(grpnum - 1) + tilespc].x == atom.x - 1) {
                cSim.pImageNMRCOMTorsionCOM[grpnum - 1 + tilespc].y = atom.x;
              }
              else {
                cSim.pImageNMRCOMTorsionCOM[grpnum + tilespc].x = atom.x;
                cSim.pImageNMRCOMTorsionCOM[grpnum + tilespc].y = atom.x;
                grpnum += 1;
                cSim.pImageNMRCOMTorsionCOMGrp[pos * 4].y = grpnum + tilespc;
              }
            }
          }
        }
        atom.x = cSim.pNMRCOMTorsionID1[pos].x;
      }
      // start COM 2
      if (atom.y > 0) {
        atom.y = cSim.pImageAtomLookup[atom.y];
      }
      else {
        int2 nmrcom_range = cSim.pNMRCOMTorsionCOMGrp[pos * 4 + 1];
        int grpnum;
        if (atom.x > 0) {
          grpnum = 0;
        }
        else {
           grpnum = cSim.pImageNMRCOMTorsionCOMGrp[pos * 4].y - tilespc;
        }
        int grpnum2 = 0;
        for (int ip = nmrcom_range.x; ip < nmrcom_range.y; ip++) {
          for (int ip2 = cSim.pNMRCOMTorsionCOM[ip].x;
               ip2 < cSim.pNMRCOMTorsionCOM[ip].y + 1; ip2++) {
            atom.y = cSim.pImageAtomLookup[ip2];
            if (grpnum2 == 0) {
              grpnum2 = grpnum;
              cSim.pImageNMRCOMTorsionCOM[grpnum2 + tilespc].x = atom.y;
              cSim.pImageNMRCOMTorsionCOM[grpnum2 + tilespc].y = atom.y;
              cSim.pImageNMRCOMTorsionCOMGrp[pos*4 + 1].x = grpnum2 + tilespc;
              grpnum2 += 1;
              cSim.pImageNMRCOMTorsionCOMGrp[pos*4 + 1].y = grpnum2 + tilespc;
            }
            else {
              if (cSim.pImageNMRCOMTorsionCOM[(grpnum2 - 1) + tilespc].y == atom.y - 1) {
                cSim.pImageNMRCOMTorsionCOM[grpnum2 - 1 + tilespc].y = atom.y;
              }
              else {
                cSim.pImageNMRCOMTorsionCOM[grpnum2 + tilespc].x = atom.y;
                cSim.pImageNMRCOMTorsionCOM[grpnum2 + tilespc].y = atom.y;
                grpnum2 += 1;
                cSim.pImageNMRCOMTorsionCOMGrp[pos*4 + 1].y = grpnum2 + tilespc;
              }
            }
          }
        }
        atom.y = cSim.pNMRCOMTorsionID1[pos].y;
      }
      // start COM 3
      if (atom.z > 0) {
        atom.z = cSim.pImageAtomLookup[atom.z];
      }
      else {
        int2 nmrcom_range = cSim.pNMRCOMTorsionCOMGrp[pos * 4 + 2];
        int grpnum;
        int grpnum2;
        if (atom.x > 0) {
          grpnum = 0;
        }
        else {
           grpnum = cSim.pImageNMRCOMTorsionCOMGrp[pos * 4].y - tilespc;
        }
        if (atom.y > 0) {
          grpnum2 = 0;
        }
        else {
           grpnum2 = cSim.pImageNMRCOMTorsionCOMGrp[pos * 4 + 1].y - tilespc;
        }
        grpnum2 += grpnum;
        int grpnum3 = 0;
        for (int ip = nmrcom_range.x; ip < nmrcom_range.y; ip++) {
          for (int ip3 = cSim.pNMRCOMTorsionCOM[ip].x;
               ip3 < cSim.pNMRCOMTorsionCOM[ip].y + 1; ip3++) {
            atom.z = cSim.pImageAtomLookup[ip3];
            if (grpnum3 == 0) {
              grpnum3 = grpnum2;
              cSim.pImageNMRCOMTorsionCOM[grpnum3 + tilespc].x = atom.z;
              cSim.pImageNMRCOMTorsionCOM[grpnum3 + tilespc].y = atom.z;
              cSim.pImageNMRCOMTorsionCOMGrp[pos*4 + 2].x = grpnum3 + tilespc;
              grpnum3 += 1;
              cSim.pImageNMRCOMTorsionCOMGrp[pos*4 + 2].y = grpnum3 + tilespc;
            }
            else {
              if (cSim.pImageNMRCOMTorsionCOM[(grpnum3 - 1) + tilespc].y == atom.z - 1) {
                cSim.pImageNMRCOMTorsionCOM[grpnum3 - 1 + tilespc].y = atom.z;
              }
              else {
                cSim.pImageNMRCOMTorsionCOM[grpnum3 + tilespc].x = atom.z;
                cSim.pImageNMRCOMTorsionCOM[grpnum3 + tilespc].y = atom.z;
                grpnum3 += 1;
                cSim.pImageNMRCOMTorsionCOMGrp[pos*4 + 2].y = grpnum3 + tilespc;
              }
            }
          }
        }
        atom.z = cSim.pNMRCOMTorsionID1[pos].z;
      }
      // start COM 4
      if (atom.w > 0) {
        atom.w = cSim.pImageAtomLookup[atom.w];
      }
      else {
        int2 nmrcom_range = cSim.pNMRCOMTorsionCOMGrp[pos * 4 + 3];
        int grpnum;
        int grpnum2;
        int grpnum3;
        if (atom.x > 0) {
          grpnum = 0;
        }
        else {
           grpnum = cSim.pImageNMRCOMTorsionCOMGrp[pos * 4].y - tilespc;
        }
        if (atom.y > 0) {
          grpnum2 = 0;
        }
        else {
           grpnum2 = cSim.pImageNMRCOMTorsionCOMGrp[pos * 4 + 1].y - tilespc;
        }
        if (atom.z > 0) {
          grpnum3 = 0;
        }
        else {
           grpnum3 = cSim.pImageNMRCOMTorsionCOMGrp[pos * 4 + 2].y - tilespc;
        }
        grpnum3 = grpnum3 + grpnum2 + grpnum;
        int grpnum4 = 0;
        for (int ip = nmrcom_range.x; ip < nmrcom_range.y; ip++) {
          for (int ip4 = cSim.pNMRCOMTorsionCOM[ip].x;
               ip4 < cSim.pNMRCOMTorsionCOM[ip].y + 1; ip4++) {
            atom.w = cSim.pImageAtomLookup[ip4];
            if (grpnum4 == 0) {
              grpnum4 = grpnum3;
              cSim.pImageNMRCOMTorsionCOM[grpnum4 + tilespc].x = atom.w;
              cSim.pImageNMRCOMTorsionCOM[grpnum4 + tilespc].y = atom.w;
              cSim.pImageNMRCOMTorsionCOMGrp[pos*4 + 3].x = grpnum4 + tilespc;
              grpnum4 += 1;
              cSim.pImageNMRCOMTorsionCOMGrp[pos*4 + 3].y = grpnum4 + tilespc;
            }
            else {
              if (cSim.pImageNMRCOMTorsionCOM[(grpnum4 - 1) + tilespc].y == atom.w - 1) {
                cSim.pImageNMRCOMTorsionCOM[grpnum4 - 1 + tilespc].y = atom.w;
              }
              else {
                cSim.pImageNMRCOMTorsionCOM[grpnum4 + tilespc].x = atom.w;
                cSim.pImageNMRCOMTorsionCOM[grpnum4 + tilespc].y = atom.w;
                grpnum4 += 1;
                cSim.pImageNMRCOMTorsionCOMGrp[pos*4 + 3].y = grpnum4 + tilespc;
              }
            }
          }
        }
        atom.w = cSim.pNMRCOMTorsionID1[pos].w;
      }
      cSim.pImageNMRCOMTorsionID1[pos] = atom;
    }
    pos += cSim.NMRCOMAngleOffset + (blockDim.x * gridDim.x);
  }
}

//---------------------------------------------------------------------------------------------
// kNLRemapLocalInteractions: this is called at the end of gpu_build_neighbor_list_() (see
//                            gpu.cpp) and assigns local (read: bonded) interactions to groups
//                            of atoms in the neighbor list for processing in the context of
//                            the hash table.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kNLRemapLocalInteractions(gpuContext gpu)
{
  kNLRemapLocalInteractions_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
  LAUNCHERROR("kNLRemapLocalInteractions");
}

//---------------------------------------------------------------------------------------------
// kNLRemapBondWorkUnits_kernel: kernel for remapping the atom imports needed by bond work
//                               units to the locally reordered neighbor list image.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(BOND_WORK_UNIT_THREADS_PER_BLOCK, 1)
kNLRemapBondWorkUnits_kernel()
{
  int wuidx = blockIdx.x;
  while (wuidx < cSim.bondWorkUnits) {
    int readpos = (3*wuidx + 1) * BOND_WORK_UNIT_THREADS_PER_BLOCK + threadIdx.x;
    int writepos = readpos + BOND_WORK_UNIT_THREADS_PER_BLOCK;
    int atomID = cSim.pBwuInstructions[readpos];
    if (atomID <= 0x7fffff) {
      cSim.pBwuInstructions[writepos] = cSim.pImageAtomLookup[atomID];
    }
    wuidx += gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kNLRemapBondWorkUnits:
//---------------------------------------------------------------------------------------------
extern "C" void kNLRemapBondWorkUnits(gpuContext gpu)
{
  if (gpu->bondWorkBlocks == 0) {
    return;
  }
  kNLRemapBondWorkUnits_kernel<<<gpu->bondWorkBlocks, BOND_WORK_UNIT_THREADS_PER_BLOCK>>>();
  LAUNCHERROR("kNLRemapBondWorkUnits");
}

//---------------------------------------------------------------------------------------------
// kNLClearCellBoundaries_kernel: clear all cell boundaries in case some are empty.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(THREADS_PER_BLOCK, 1)
kNLClearCellBoundaries_kernel()
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  uint2 nulldata   = {0, 0};
  while (pos < cSim.cells) {
    cSim.pNLNonbondCellStartEnd[pos]  = nulldata;
    pos += blockDim.x * gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kNLClearCellBoundaries: launch the kernel to do what the name says.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kNLClearCellBoundaries(gpuContext gpu)
{
  kNLClearCellBoundaries_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
  LAUNCHERROR("kNLClearCellBoundaries");
}

//---------------------------------------------------------------------------------------------
// kNLCalculateCellBoundaries_kernel:
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(THREADS_PER_BLOCK, 1)
kNLCalculateCellBoundaries_kernel(unsigned int* pImageHash)
{
  const int cSpan = 12000;
  __shared__ unsigned int sHash[cSpan];
  int pos = ((cSim.atoms + 1) * blockIdx.x) / gridDim.x - 1;
  int end = ((cSim.atoms + 1) * (blockIdx.x + 1)) / gridDim.x;

  while (pos < end) {

    // Read span to check for transitions
    int pos1  = pos + threadIdx.x;
    int spos  = threadIdx.x;
    int span  = min(end, pos + cSpan);
    int lSpan = min(end, pos + cSpan) - pos;
    while (pos1 < span) {

      // Read hash data or 0s on either end to force transitions
      if ((pos1 >= 0) && (pos1 < cSim.atoms))
        sHash[spos]= pImageHash[pos1];
      else
        sHash[spos]= 0;
      pos1 += blockDim.x;
      spos += blockDim.x;
    }
    __syncthreads();
    spos = threadIdx.x + 1;
    while (spos < lSpan) {
      int oldHash = sHash[spos - 1] >> CELL_HASH_BITS;
      int newHash = sHash[spos] >> CELL_HASH_BITS;
      if (oldHash != newHash) {
        if (pos + spos != cSim.atoms)
          cSim.pNLNonbondCellStartEnd[newHash].x  = pos + spos;
        if (pos + spos != 0)
          cSim.pNLNonbondCellStartEnd[oldHash].y  = pos + spos;
      }
      spos += blockDim.x;
    }
    __syncthreads();
    pos = min(end, pos + cSpan - 1);
  }
}

//---------------------------------------------------------------------------------------------
// kNLCalculateCellBoundaries:
//---------------------------------------------------------------------------------------------
extern "C" void kNLCalculateCellBoundaries(gpuContext gpu)
{
  kNLCalculateCellBoundaries_kernel<<<gpu->blocks,
                                      gpu->threadsPerBlock>>>(gpu->sim.pImageHash);
  LAUNCHERROR("kNLCalculateCellBoundaries");
}

//---------------------------------------------------------------------------------------------
// GEData: directives for the work units of the non-bonded interactions--neighbor list entries.
//---------------------------------------------------------------------------------------------
struct GEData {
  unsigned int workUnit;
  unsigned int xCellStart;
  unsigned int xCellEnd;
  unsigned int yCellStart;
  unsigned int yCellEnd;
  unsigned int exclusionMap;
  unsigned int imageAtom;
  unsigned int atom;
};

//---------------------------------------------------------------------------------------------
// Kernels to compute atomic coordinates as they will be needed within the cell grids.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kNLCalcCellCoordinates_kernel()
#include "kReImageCoord.h"

//---------------------------------------------------------------------------------------------
#define PME_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kNLCalcCellCoordinatesOrthogonal_kernel()
#include "kReImageCoord.h"
#undef PME_ORTHOGONAL

//---------------------------------------------------------------------------------------------
#define PME_NTP
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kNLCalcCellCoordinatesNTP_kernel()
#include "kReImageCoord.h"

//---------------------------------------------------------------------------------------------
#define PME_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kNLCalcCellCoordinatesOrthogonalNTP_kernel()
#include "kReImageCoord.h"
#undef PME_ORTHOGONAL
#undef PME_NTP

//---------------------------------------------------------------------------------------------
extern "C" void kNLCalculateCellCoordinates(gpuContext gpu)
{
  // Set tile size in local variables
  int nBlocks = gpu->blocks;
  int genThreads = gpu->generalThreadsPerBlock;

  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
    if (gpu->sim.is_orthog) {
      kNLCalcCellCoordinatesOrthogonalNTP_kernel<<<nBlocks, genThreads>>>();
    }
    else {
      kNLCalcCellCoordinatesNTP_kernel<<<nBlocks, genThreads>>>();
    }
  }
  else {
    if (gpu->sim.is_orthog) {
      kNLCalcCellCoordinatesOrthogonal_kernel<<<nBlocks, genThreads>>>();
    }
    else {
      kNLCalcCellCoordinates_kernel<<<nBlocks, genThreads>>>();
    }
  }
  LAUNCHERROR("kNLCalculateCellCoordinates");
}

//---------------------------------------------------------------------------------------------
// Kernels for neighbor list building.  Specialization depends on three factors: whether there
// are 8, 16, or 32 atoms per warp, whether the box is orthogonal, and whether constant
// pressure conditions are in effect.
//---------------------------------------------------------------------------------------------
#if defined(AMBER_PLATFORM_AMD)
#  define INC_BNL "kBNL_AMD.h"
#else
#  define INC_BNL "kBNL.h"
#endif

// HIP-TODO: Implement a better way (NLBUILD_NEIGHBORLIST_BLOCKS_MULTIPLIER is used for BNLBlocks
// But it can't be high in __launch_bounds__)
#ifdef AMBER_PLATFORM_AMD
#  define __LAUNCH_BOUNDS_BNL__(X, Y, Z) __launch_bounds__(X, Z)
#else
#  define __LAUNCH_BOUNDS_BNL__(X, Y, Z) __launch_bounds__(X, Y)
#endif

#define PME_ATOMS_PER_WARP (32)
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_BNL__(NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK,
                      NLBUILD_NEIGHBORLIST_BLOCKS_MULTIPLIER,
                      NLBUILD_NEIGHBORLIST_OCCUPANCY)
kNLMakeOrtho32_kernel()
#include INC_BNL

//---------------------------------------------------------------------------------------------
#define PME_VIRIAL
__global__ void
__LAUNCH_BOUNDS_BNL__(NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK,
                      NLBUILD_NEIGHBORLIST_BLOCKS_MULTIPLIER,
                      NLBUILD_NEIGHBORLIST_OCCUPANCY)
kNLMakeOrthoNTP32_kernel()
#include INC_BNL
#undef PME_VIRIAL
#undef PME_IS_ORTHOGONAL

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_BNL__(NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK,
                      NLBUILD_NEIGHBORLIST_BLOCKS_MULTIPLIER,
                      NLBUILD_NEIGHBORLIST_OCCUPANCY)
kNLMake32_kernel()
#include INC_BNL

//---------------------------------------------------------------------------------------------
#define PME_VIRIAL
__global__ void
__LAUNCH_BOUNDS_BNL__(NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK,
                      NLBUILD_NEIGHBORLIST_BLOCKS_MULTIPLIER,
                      NLBUILD_NEIGHBORLIST_OCCUPANCY)
kNLMakeNTP32_kernel()
#include INC_BNL
#undef PME_VIRIAL
#undef PME_ATOMS_PER_WARP

//---------------------------------------------------------------------------------------------
#define PME_ATOMS_PER_WARP (16)
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_BNL__(NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK,
                      NLBUILD_NEIGHBORLIST_BLOCKS_MULTIPLIER,
                      NLBUILD_NEIGHBORLIST_OCCUPANCY)
kNLMakeOrtho16_kernel()
#include INC_BNL

//---------------------------------------------------------------------------------------------
#define PME_VIRIAL
__global__ void
__LAUNCH_BOUNDS_BNL__(NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK,
                      NLBUILD_NEIGHBORLIST_BLOCKS_MULTIPLIER,
                      NLBUILD_NEIGHBORLIST_OCCUPANCY)
kNLMakeOrthoNTP16_kernel()
#include INC_BNL
#undef PME_VIRIAL
#undef PME_IS_ORTHOGONAL

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_BNL__(NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK,
                      NLBUILD_NEIGHBORLIST_BLOCKS_MULTIPLIER,
                      NLBUILD_NEIGHBORLIST_OCCUPANCY)
kNLMake16_kernel()
#include INC_BNL

//---------------------------------------------------------------------------------------------
#define PME_VIRIAL
__global__ void
__LAUNCH_BOUNDS_BNL__(NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK,
                      NLBUILD_NEIGHBORLIST_BLOCKS_MULTIPLIER,
                      NLBUILD_NEIGHBORLIST_OCCUPANCY)
kNLMakeNTP16_kernel()
#include INC_BNL
#undef PME_VIRIAL
#undef PME_ATOMS_PER_WARP

//---------------------------------------------------------------------------------------------
#define PME_ATOMS_PER_WARP (8)
#define PME_IS_ORTHOGONAL
__global__ void
__LAUNCH_BOUNDS_BNL__(NLBUILD_NEIGHBORLIST8_THREADS_PER_BLOCK,
                      NLBUILD_NEIGHBORLIST_BLOCKS_MULTIPLIER,
                      NLBUILD_NEIGHBORLIST_OCCUPANCY)
kNLMakeOrtho8_kernel()
#include INC_BNL

//---------------------------------------------------------------------------------------------
#define PME_VIRIAL
__global__ void
__LAUNCH_BOUNDS_BNL__(NLBUILD_NEIGHBORLIST8_THREADS_PER_BLOCK,
                      NLBUILD_NEIGHBORLIST_BLOCKS_MULTIPLIER,
                      NLBUILD_NEIGHBORLIST_OCCUPANCY)
kNLMakeOrthoNTP8_kernel()
#include INC_BNL
#undef PME_VIRIAL
#undef PME_IS_ORTHOGONAL

//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS_BNL__(NLBUILD_NEIGHBORLIST8_THREADS_PER_BLOCK,
                      NLBUILD_NEIGHBORLIST_BLOCKS_MULTIPLIER,
                      NLBUILD_NEIGHBORLIST_OCCUPANCY)
kNLMake8_kernel()
#include INC_BNL

//---------------------------------------------------------------------------------------------
#define PME_VIRIAL
__global__ void
__LAUNCH_BOUNDS_BNL__(NLBUILD_NEIGHBORLIST8_THREADS_PER_BLOCK,
                      NLBUILD_NEIGHBORLIST_BLOCKS_MULTIPLIER,
                      NLBUILD_NEIGHBORLIST_OCCUPANCY)
kNLMakeNTP8_kernel()
#include INC_BNL
#undef PME_VIRIAL
#undef PME_ATOMS_PER_WARP

#undef INC_BNL
#undef __LAUNCH_BOUNDS_BNL__

//---------------------------------------------------------------------------------------------
// kNeighborListInitKernels: initialize kernels for neighbor list building.  This function
//                           accomplishes one thing: sets the cache configuration to prefer L1,
//                           which prompts more register space and less shared memory space.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kNeighborListInitKernels(gpuContext gpu)
{
  if (gpu->sm_version >= SM_3X) {
    cudaFuncSetCacheConfig(kNLMakeOrthoNTP32_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kNLMakeNTP32_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kNLMakeOrtho32_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kNLMake32_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kNLMakeOrthoNTP16_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kNLMakeNTP16_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kNLMakeOrtho16_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kNLMake16_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kNLMakeOrthoNTP8_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kNLMakeNTP8_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kNLMakeOrtho8_kernel, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(kNLMake8_kernel, cudaFuncCachePreferL1);
  }
}

//---------------------------------------------------------------------------------------------
// kNLBuildNeighborList: launch the appropriate kernel to build the neighbor list.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kNLBuildNeighborList(gpuContext gpu)
{
  if (gpu->sim.NLAtomsPerWarp == 32) {
    if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
      if (gpu->sim.is_orthog)
        kNLMakeOrthoNTP32_kernel<<<gpu->BNLBlocks, gpu->NLBuildNeighborList32ThreadsPerBlock>>>();
      else
        kNLMakeNTP32_kernel<<<gpu->BNLBlocks, gpu->NLBuildNeighborList32ThreadsPerBlock>>>();
    }
    else {
      if (gpu->sim.is_orthog)
        kNLMakeOrtho32_kernel<<<gpu->BNLBlocks, gpu->NLBuildNeighborList32ThreadsPerBlock>>>();
      else
        kNLMake32_kernel<<<gpu->BNLBlocks, gpu->NLBuildNeighborList32ThreadsPerBlock>>>();
    }
  }
  else if (gpu->sim.NLAtomsPerWarp == 16) {
    if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
      if (gpu->sim.is_orthog)
        kNLMakeOrthoNTP16_kernel<<<gpu->BNLBlocks, gpu->NLBuildNeighborList16ThreadsPerBlock>>>();
      else
        kNLMakeNTP16_kernel<<<gpu->BNLBlocks, gpu->NLBuildNeighborList16ThreadsPerBlock>>>();
    }
    else {
      if (gpu->sim.is_orthog)
        kNLMakeOrtho16_kernel<<<gpu->BNLBlocks, gpu->NLBuildNeighborList16ThreadsPerBlock>>>();
      else
        kNLMake16_kernel<<<gpu->BNLBlocks, gpu->NLBuildNeighborList16ThreadsPerBlock>>>();
    }
  }
  else {
    if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
      if (gpu->sim.is_orthog)
        kNLMakeOrthoNTP8_kernel<<<gpu->BNLBlocks, gpu->NLBuildNeighborList8ThreadsPerBlock>>>();
      else
        kNLMakeNTP8_kernel<<<gpu->BNLBlocks, gpu->NLBuildNeighborList8ThreadsPerBlock>>>();
    }
    else {
      if (gpu->sim.is_orthog)
        kNLMakeOrtho8_kernel<<<gpu->BNLBlocks, gpu->NLBuildNeighborList8ThreadsPerBlock>>>();
      else
        kNLMake8_kernel<<<gpu->BNLBlocks, gpu->NLBuildNeighborList8ThreadsPerBlock>>>();
    }
  }
  LAUNCHERROR("kNLBuildNeighborList");
}

//---------------------------------------------------------------------------------------------
// kNLSkinTestCount_kernel: kernel to test all atom movements against the current pair list.
//                          If one has moved too far from its original position, the "skin"
//                          test fails.  This also tracks the number of atoms that have moved
//                          a certain portion of the pair list margin "skin" and prints the
//                          result at the end.
//
//                          Before any of that, this will loop over all atoms and check for
//                          violations.  Because the pair list violations do become more
//                          frequent in larger systems (albeit slowly), it could be useful to
//                          refresh the pair list less frequently, taking a slight risk of
//                          missing an interaction near the cutoff.  This will help test how
//                          much of a risk letting particles travel beyond half the pair list
//                          margin really is.
//
// This is a debugging function.
//---------------------------------------------------------------------------------------------
__global__ void __LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1) kNLSkinTestCount_kernel()
{
  __shared__ volatile bool sbFail;
  __shared__ volatile PMEFloat sOne_half_nonbond_skin_squared, sPermitted_skin_squared;
  __shared__ volatile int sMovement[512];
  __shared__ volatile int nviolation;
  __shared__ volatile PMEFloat sRecipf[9], sUcellf[9];
  __shared__ volatile int nPairObs[GENERAL_THREADS_PER_BLOCK / GRID];

  if (threadIdx.x < 9) {
    if (cSim.ntp > 0 && cSim.barostat == 1) {
      sRecipf[threadIdx.x] = cSim.pNTPData->recipf[threadIdx.x];
      sUcellf[threadIdx.x] = cSim.pNTPData->ucellf[threadIdx.x];
    }
    else {
      int tmod3 = threadIdx.x % 3;
      int tdiv3 = threadIdx.x / 3;
      sRecipf[threadIdx.x] = cSim.recipf[tmod3][tdiv3];
      sUcellf[threadIdx.x] = cSim.ucellf[tmod3][tdiv3];
    }
  }
  if (threadIdx.x == 0) {
    nviolation = 0;
  }
  __syncthreads();

  int i = blockIdx.x;
  int npair = 0;
  while (i < cSim.atoms) {
    int j = i + 1 + threadIdx.x;
    PMEFloat xi = cSim.pImageX[i];
    PMEFloat yi = cSim.pImageY[i];
    PMEFloat zi = cSim.pImageZ[i];
    PMEFloat2 oldxy = cSim.pAtomXYSaveSP[i];
    PMEFloat xiold  = oldxy.x;
    PMEFloat yiold  = oldxy.y;
    PMEFloat ziold  = cSim.pAtomZSaveSP[i];
    while (j < cSim.atoms) {
      PMEFloat dx = cSim.pImageX[j] - xi;
      PMEFloat dy = cSim.pImageY[j] - yi;
      PMEFloat dz = cSim.pImageZ[j] - zi;
      PMEFloat ndx = sRecipf[0]*dx + sRecipf[3]*dy + sRecipf[6]*dz;
      PMEFloat ndy =                 sRecipf[4]*dy + sRecipf[7]*dz;
      PMEFloat ndz =                                 sRecipf[8]*dz;
      ndx -= round(ndx);
      ndy -= round(ndy);
      ndz -= round(ndz);
      dx = sUcellf[0]*ndx + sUcellf[1]*ndy + sUcellf[2]*ndz;
      dy =                  sUcellf[4]*ndy + sUcellf[5]*ndz;
      dz =                                   sUcellf[8]*ndz;
      PMEFloat r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < cSim.cut2) {
        npair++;
        oldxy           = cSim.pAtomXYSaveSP[j];
        PMEFloat dxold  = oldxy.x - xiold;
        PMEFloat dyold  = oldxy.y - yiold;
        PMEFloat dzold  = cSim.pAtomZSaveSP[i] - ziold;
        PMEFloat ndxold = sRecipf[0]*dxold + sRecipf[3]*dyold + sRecipf[6]*dzold;
        PMEFloat ndyold =                    sRecipf[4]*dyold + sRecipf[7]*dzold;
        PMEFloat ndzold =                                       sRecipf[8]*dzold;
        ndxold -= round(ndxold);
        ndyold -= round(ndyold);
        ndzold -= round(ndzold);
        dxold   = sUcellf[0]*ndxold + sUcellf[1]*ndyold + sUcellf[2]*ndzold;
        dyold   =                     sUcellf[4]*ndyold + sUcellf[5]*ndzold;
        dzold   =                                         sUcellf[8]*ndzold;
        PMEFloat r2old = dxold*dxold + dyold*dyold + dzold*dzold;
        if (r2old > cSim.cutPlusSkin2) {
          atomicAdd((int*)&nviolation, 1);
        }
      }
      j += blockDim.x;
    }
    i += gridDim.x;
  }
  __syncthreads();

  // Reduce the pairs
#ifdef AMBER_PLATFORM_AMD_WARP64
  npair += __SHFL(WARP_MASK, npair, 32);
#endif
  npair += __SHFL(WARP_MASK, npair, 16);
  npair += __SHFL(WARP_MASK, npair, 8);
  npair += __SHFL(WARP_MASK, npair, 4);
  npair += __SHFL(WARP_MASK, npair, 2);
  npair += __SHFL(WARP_MASK, npair, 1);
  if ((threadIdx.x & GRID_BITS_MASK) == 0) {
    int warpIdx = threadIdx.x / GRID;
    nPairObs[warpIdx] = npair;
  }
  __syncthreads();
  npair = 0;
  if (threadIdx.x < GENERAL_THREADS_PER_BLOCK / GRID) {
    npair = nPairObs[threadIdx.x];
  }
#ifdef AMBER_PLATFORM_AMD_WARP64
  npair += __SHFL(WARP_MASK, npair, 32);
#endif
  npair += __SHFL(WARP_MASK, npair, 16);
  npair += __SHFL(WARP_MASK, npair, 8);
  npair += __SHFL(WARP_MASK, npair, 4);
  npair += __SHFL(WARP_MASK, npair, 2);
  npair += __SHFL(WARP_MASK, npair, 1);
  if (threadIdx.x == 0) {
    printf("%% kNLSkinTestCount :: %4d violations were detected out of %10d pairs.\n",
           nviolation, npair);
  }

  // Check for violations
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  sbFail = false;
  if (threadIdx.x == 0) {
    if (cSim.ntp > 0 && cSim.barostat == 1) {
       sOne_half_nonbond_skin_squared = cSim.pNTPData->one_half_nonbond_skin_squared;
    }
    else {
      sOne_half_nonbond_skin_squared = cSim.one_half_nonbond_skin_squared;
    }
  }
  for (int i = threadIdx.x; i < 501; i += blockDim.x) {
    sMovement[i] = 0;
  }
  __syncthreads();

  while (pos < cSim.atoms) {
    PMEFloat x      = cSim.pImageX[pos];
    PMEFloat y      = cSim.pImageY[pos];
    PMEFloat2 oldxy = cSim.pAtomXYSaveSP[pos];
    PMEFloat z      = cSim.pImageZ[pos];
    PMEFloat oldz   = cSim.pAtomZSaveSP[pos];
    PMEFloat dx     = x - oldxy.x;
    PMEFloat dy     = y - oldxy.y;
    PMEFloat dz     = z - oldz;
    PMEFloat r2     = dx*dx + dy*dy + dz*dz;
    int imvr        = (PMEFloat)100.0 * sqrt(r2 / sOne_half_nonbond_skin_squared);
    if (imvr < 500) {
      atomicAdd((int*)&sMovement[imvr], 1);
    }
    if (r2 >= sOne_half_nonbond_skin_squared) {
      sbFail = true;
    }
    pos += blockDim.x * gridDim.x;
  }
  __syncthreads();

  // Report results
  if ((threadIdx.x == 0) && sbFail) {
    *cSim.pNLbSkinTestFail = true;
    printf("%% Pair list was refreshed.  Permitted = %8.4f  Safe = %8.4f\n",
           sqrt(sOne_half_nonbond_skin_squared), sqrt(sOne_half_nonbond_skin_squared));
  }
  if ((threadIdx.x == 0) && nviolation > 0) {
    printf("movement = [\n");
    int i;
    for (i = 0; i < 200; i++) {
      printf("  %5.2f  %5d\n", (PMEFloat)0.01 * (PMEFloat)i, sMovement[i]);
    }
    printf("];\n");
  }
}

//---------------------------------------------------------------------------------------------
// kNLSkinTest_kernel: kernel to test all atom movements against the current pair list.  If
//                     one has moved too far from its original position, the "skin" test fails.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(GENERAL_THREADS_PER_BLOCK, 1)
kNLSkinTest_kernel()
{
  __shared__ bool sbFail;
  __shared__ PMEFloat sOne_half_nonbond_skin_squared;
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  sbFail = false;
  if (((cSim.ntp > 0) && (cSim.barostat == 1)) && (threadIdx.x == 0)) {
     sOne_half_nonbond_skin_squared = cSim.pNTPData->one_half_nonbond_skin_squared;
  }
  __syncthreads();

  if ((cSim.ntp > 0) && (cSim.barostat == 1)) {
    while (pos < cSim.atoms) {
      PMEFloat x      = cSim.pImageX[pos];
      PMEFloat y      = cSim.pImageY[pos];
      PMEFloat2 oldxy = cSim.pAtomXYSaveSP[pos];
      PMEFloat z      = cSim.pImageZ[pos];
      PMEFloat oldz   = cSim.pAtomZSaveSP[pos];
      PMEFloat dx     = x - oldxy.x;
      PMEFloat dy     = y - oldxy.y;
      PMEFloat dz     = z - oldz;
      PMEFloat r2     = dx*dx + dy*dy + dz*dz;
      if (r2 >= sOne_half_nonbond_skin_squared) {
        sbFail = true;
      }
      pos += blockDim.x * gridDim.x;
    }
  }
  else {
    while (pos < cSim.atoms) {
      PMEFloat x      = cSim.pImageX[pos];
      PMEFloat y      = cSim.pImageY[pos];
      PMEFloat2 oldxy = cSim.pAtomXYSaveSP[pos];
      PMEFloat z      = cSim.pImageZ[pos];
      PMEFloat oldz   = cSim.pAtomZSaveSP[pos];
      PMEFloat dx     = x - oldxy.x;
      PMEFloat dy     = y - oldxy.y;
      PMEFloat dz     = z - oldz;
      PMEFloat r2     = dx*dx + dy*dy + dz*dz;
      if (r2 >= cSim.one_half_nonbond_skin_squared) {
        sbFail = true;
      }
      pos += blockDim.x * gridDim.x;
    }
  }
  __syncthreads();
#if defined(AMBER_PLATFORM_AMD)
  if (threadIdx.x == 0) {
    // Pack total and failed block counts in one counter
    unsigned int delta = (sbFail ? (1 << 16) : 0) | 1;
    unsigned int count = atomicAdd(&cSim.pFinishedBlocksCounters[0], delta) + delta;
    unsigned int totalCount  = (count & 0xffff);
    unsigned int failedCount = (count >> 16);
    // The last block must signal completion to host
    if (totalCount == gridDim.x) {
      // Signal completion: 2 - failed, 1 - passed
      atomicExch(&cSim.pKernelCompletionFlags[0], (failedCount > 0) ? 2 : 1);
      atomicExch(&cSim.pFinishedBlocksCounters[0], 0);
    }
  }
#else
  if ((threadIdx.x == 0) && sbFail) {
    *cSim.pNLbSkinTestFail = true;
  }
#endif
}

//---------------------------------------------------------------------------------------------
// kNLSkinTest: host function to download tests on the neighbor list buffer skin.  If a
//              violation is detected, this will call for a pair list rebuild.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------

extern "C" void kNLSkinTest(gpuContext gpu)
{
  kNLSkinTest_kernel<<<gpu->blocks, gpu->generalThreadsPerBlock>>>();
  LAUNCHERROR("kNLSkinTest");

#if defined(AMBER_PLATFORM_AMD)
  while (true) {
    unsigned int result = __atomic_load_n(&gpu->pbKernelCompletionFlags->_pSysData[0],
                                          __ATOMIC_RELAXED);
    if (result != 0) {
      if (result == 2) {
        gpu->bNeedNewNeighborList = true;
      }
      // Clear the completion flag for future uses
      __atomic_store_n(&gpu->pbKernelCompletionFlags->_pSysData[0], 0, __ATOMIC_RELAXED);
      break;
    }
  }
#else
  cudaDeviceSynchronize();

  if (!(gpu->bCanMapHostMemory)) {
    gpu->pbNLbSkinTestFail->Download();
  }
  if (*(gpu->pbNLbSkinTestFail->_pSysData)) {
    gpu->bNeedNewNeighborList = true;
    *(gpu->pbNLbSkinTestFail->_pSysData) = false;
    if (!(gpu->bCanMapHostMemory)) {
      gpu->pbNLbSkinTestFail->Upload();
    }
  }
#endif
}
