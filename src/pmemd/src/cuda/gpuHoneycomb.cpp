#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------
#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include "matrix.h"
#include "gpuContext.h"
#include "gpuHoneycomb.h"

using namespace std;

//---------------------------------------------------------------------------------------------
// CountQQAtoms: count the number of atoms beraing partial charges that would allow them to
//               participate in electrostatic interactions.
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                object counts
//---------------------------------------------------------------------------------------------
int CountQQAtoms(gpuContext gpu)
{
  int i, nqq;

  nqq = 0;
  for (i = 0; i < gpu->sim.atoms; i++) {
    nqq += (fabs(gpu->pbAtomChargeSP->_pSysData[i]) > 1.0e-8);
  }

  return nqq;
}

//---------------------------------------------------------------------------------------------
// CountLJAtoms: count the number of atoms bearing Lennard-Jones properties that would allow
//               them to participate in non-trivial van-der Waals interactions.
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                object counts
//---------------------------------------------------------------------------------------------
int CountLJAtoms(gpuContext gpu)
{
  int i, nlj;

  nlj = 0;
  for (i = 0; i < gpu->sim.atoms; i++) {
    int ljt = gpu->pbAtomLJID->_pSysData[i] - 1;
    if (ljt >= 0) {
      nlj += ((gpu->sim.ljexist[ljt >> 5] & (0x1 << (ljt & 31))) > 0);
    }

  }

  return nlj;
}

//---------------------------------------------------------------------------------------------
// GetNTInstructionEst: calculate an estimate of the number of NT regions that will be needed
//                      in order to compute all of the non-bonded interactions.  This number
//                      will then be multiplied by the size of on NT instruction set to yield
//                      the total amount of space to allocate for the NT instructions.
//
// Arguments:
//   gpu:   overarching data structure containing simulation information, here used for atom
//          and cell counts
//---------------------------------------------------------------------------------------------
int GetNTInstructionEst(gpuContext gpu)
{
#ifdef use_DPFP
  double divfac    = gpu->sim.npencils * 350;
#else
  double divfac    = gpu->sim.npencils * 480;
#endif
  double cutfac    = (2.0 * gpu->sim.LowerCapBoundary) + 1.0;
  double meanQQreg = (double)(gpu->sim.nQQatoms * 11) * cutfac / divfac;
  double meanLJreg = (double)(gpu->sim.nLJatoms * 11) * cutfac / divfac;
  double qqlayers = gpu->sim.LowerCapBoundary / cutfac * meanQQreg;
  double ljlayers = gpu->sim.LowerCapBoundary / cutfac * meanLJreg;
  int ntre = 2.25 * ((meanQQreg + 1.0)*qqlayers + (meanLJreg + 1.0)*ljlayers + 4.0);
  ntre *= gpu->sim.npencils;

  return ntre;
}

//---------------------------------------------------------------------------------------------
// FindBisector: compute the bisector of a line segment determined by two points in a plane.
//               Provision is made for a bisector with infinite slope.
//
// Arguments:
//   y1, z1:      location of the first point
//   y2, z2:      location of the second point
//   slope:       slope of the bisecting line (returned)
//   intr:        intercept of the bisecting line (returned)
//---------------------------------------------------------------------------------------------
static void FindBisector(double y1, double z1, double y2, double z2, double *slope,
                         double *intr)
{
  // First, find the midpoint
  double mdpty = 0.5*(y1 + y2);
  double mdptz = 0.5*(z1 + z2);

  // Now find the slope of the line between the points
  double dy = y2 - y1;
  double dz = z2 - z1;

  // Slope and intercept of the bisector
  if (fabs(dz) < 1.0e-8) {
    *slope = 1.0e8;
    *intr  = mdpty;
  }
  else {

    // The slope of the bisector is the inverse
    double spty = mdpty + 1.0;
    double sptz = mdptz - dy/dz;
    *slope = (mdptz - sptz) / (mdpty - spty);
    *intr  = mdptz - (*slope)*mdpty;
  }
}

//---------------------------------------------------------------------------------------------
// FindBsctIntercept: finds the point of intersection for two lines.  It is assumed that the
//                    lines DO intersect.  There is also a provision for one of the lines
//                    having infinite slope (denoted by slope = 1.0e8).  In such a case, the
//                    intercept of that line refers to the horizontal axis intercept, and the
//                    intersection will be computed based on the non-infinite slope of the
//                    other line.
//
// Arguments:
//   slp1, intr1:   slope and intercept of the first line
//   slp2, intr2:   slope and intercept of the second line
//   y, z;          point of intersection (returned)
//---------------------------------------------------------------------------------------------
static inline void FindBsctIntercept(double slp1, double intr1, double slp2, double intr2,
                                     double *y, double *z)
{
  if (slp1 > 9.999e7) {
    *y = intr1;
    *z = intr1*slp2 + intr2;
  }
  else if (slp2 > 9.999e7) {
    *y = intr2;
    *z = intr2*slp1 + intr1;
  }
  else {

    // y*slp1 + intr1 = y*slp2 + intr2
    // y*(slp1 - slp2) = (intr2 - intr1)

    *y = (intr2 - intr1)/(slp1 - slp2);
    *z = (*y)*slp1 + intr1;  
  }
}

//---------------------------------------------------------------------------------------------
// GetPencilVertices: compute the vertices of a pencil's cross section.
//
// Arguments:
//   
//---------------------------------------------------------------------------------------------
static dmat GetPencilVertices(int ny, int nz, bool yshift, dmat *invUxs)
{
  int i;
  dmat bisectors, vertices;

  // Lay out the surrounding six points
  double invdy = 1.0 / (double)ny;
  double invdz = 1.0 / (double)nz;
  double *dtmp, *vtmp;
  dtmp = invUxs->data;
  vertices = CreateDmat(6, 2);

  // If the rows must shift, add 0.5 spaces in +/- Y to points above or below the Y axis
  double ytrans = (yshift) ? 0.5 / (double)ny : 0.0;
  
  // Point along Y-axis -1 spaces from the origin
  vertices.map[0][0] = dtmp[0] * (-invdy);

  // Point one layer above
  vertices.map[1][0] = dtmp[0]*(-ytrans) + dtmp[1]*invdz;
  vertices.map[1][1] = dtmp[2]*(-ytrans) + dtmp[3]*invdz;

  // Point one layer above, translated +1 spaces in Y
  vertices.map[2][0] = dtmp[0]*(invdy - ytrans) + dtmp[1]*invdz;
  vertices.map[2][1] = dtmp[2]*(invdy - ytrans) + dtmp[3]*invdz;

  // Point along Y-axis +1 spaces from the origin
  vertices.map[3][0] = dtmp[0] * invdy;

  // Point one layer below
  vertices.map[4][0] = dtmp[0]*ytrans + dtmp[1]*(-invdz);
  vertices.map[4][1] = dtmp[2]*ytrans + dtmp[3]*(-invdz);

  // Point one layer below, translated -1 spaces in Y
  vertices.map[5][0] = dtmp[0]*(ytrans - invdy) + dtmp[1]*(-invdz);
  vertices.map[5][1] = dtmp[2]*(ytrans - invdy) + dtmp[3]*(-invdz);

  // The surrounding vertices now determine bisectors
  bisectors = CreateDmat(6, 2);
  for (i = 0; i < 6; i++) {
    FindBisector(0.0, 0.0, vertices.map[i][0], vertices.map[i][1],
                 &bisectors.map[i][0], &bisectors.map[i][1]);
  }
  
  // The bisectors' intercepts now determine the new vertices
  for (i = 0; i < 6; i++) {
    int ip1 = i + 1 - 6*(i == 5);
    FindBsctIntercept(bisectors.map[i][0], bisectors.map[i][1], bisectors.map[ip1][0],
                      bisectors.map[ip1][1], &vertices.map[i][0], &vertices.map[i][1]);
  }

  // Free allocated memory
  DestroyDmat(&bisectors);

  return vertices;
}

//---------------------------------------------------------------------------------------------
// Point2LineSegmentDistance: find the distance between a point and a line segment.
//
// Arguments:
//   x, y:     coordinates of the point
//   x1, y1:   coordinates of the first end of the line segment
//   x2, y2:   coordinates of the second end of the line segment
//---------------------------------------------------------------------------------------------
static double Point2LineSegmentDistance(double x, double y, double x1, double y1, double x2,
                                        double y2)
{
  double A = x - x1;
  double B = y - y1;
  double C = x2 - x1;
  double D = y2 - y1;

  double dot = A*C + B*D;
  double len_sq = C*C + D*D;

  // Bullet-proof against zero length line
  double param = -1.0;
  if (len_sq != 0.0) {
    param = dot / len_sq;
  }
  double xx, yy;
  if (param < 0) {
    xx = x1;
    yy = y1;
  }
  else if (param > 1) {
    xx = x2;
    yy = y2;
  }
  else {
    xx = x1 + param*C;
    yy = y1 + param*D;
  }

  double dx = x - xx;
  double dy = y - yy;

  return sqrt(dx*dx + dy*dy);
}

//---------------------------------------------------------------------------------------------
// CheckAllPencilVertices: check the distances between pencils using all vertices and cross-
//                         section edges.
//---------------------------------------------------------------------------------------------
static double CheckAllPencilVertices(double ry, double rz, double range, dmat *vertices)
{
  int i, j;

  // Check all line segments in the home pencil against vertices of the other pencil
  double mindist = range;
  for (i = 0; i < 6; i++) {
    double vy = vertices->map[i][0];
    double vz = vertices->map[i][1];
    for (j = 0; j < 6; j++) {
      double l1y = vertices->map[j][0] + ry;
      double l1z = vertices->map[j][1] + rz;
      int jpidx = j + 1 - 6*(j == 5);
      double l2y = vertices->map[jpidx][0] + ry;
      double l2z = vertices->map[jpidx][1] + rz;
      double dist = Point2LineSegmentDistance(vy, vz, l1y, l1z, l2y, l2z);
      if (dist < mindist) {
        mindist = dist;
      }
    }
  }

  return mindist;
}

//---------------------------------------------------------------------------------------------
// ComputeCapBoundary: compute the needed padding for the expanded unit cell along the X axis.
//
// Arguments:
//   lj_cutoff:   the Lennard-Jones cutoff (only one cutoff is used, and this is the larger of
//                the two)
//   ucell:       matrix describing the transformation from fractional coordinates to real
//                space, in FORTRAN order (take this from gpu->sim.ucell or the like)
//---------------------------------------------------------------------------------------------
static double ComputeCapBoundary(double lj_cutoff, double nbskin, dmat *invU)
{
  double udata[9], cdepth[3];

  udata[0] = invU->map[0][0];
  udata[1] = invU->map[0][1];
  udata[2] = invU->map[0][2];
  udata[3] = invU->map[1][0];
  udata[4] = invU->map[1][1];
  udata[5] = invU->map[1][2];
  udata[6] = invU->map[2][0];
  udata[7] = invU->map[2][1];
  udata[8] = invU->map[2][2];
  HessianNorms(udata, cdepth);
  
  return (lj_cutoff + nbskin) / cdepth[0];
}

//---------------------------------------------------------------------------------------------
// ExtractTransformationMatrix: extract a transformation matrix from the gpuContext.
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                object counts
//   inverse:     flag to extract the inverse transformation matrix (ucell instead of recip)
//---------------------------------------------------------------------------------------------
static dmat ExtractTransformationMatrix(gpuContext gpu, bool inverse)
{
  int i, j;
  dmat M;

  M = CreateDmat(3, 3);
  if (gpu->sim.ntp > 0 && gpu->sim.barostat == 1) {
    gpu->pbNTPData->Download();
  }
  if (inverse) {
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        if (gpu->sim.ntp > 0 && gpu->sim.barostat == 1) {
#ifdef use_DPFP
          M.map[i][j] = gpu->pbNTPData->_pSysData->ucell[3*j + i];
#else
          M.map[i][j] = gpu->pbNTPData->_pSysData->ucellf[3*j + i];
#endif
        }
        else {
#ifdef use_DPFP
          M.map[i][j] = gpu->sim.ucell[i][j];
#else
          M.map[i][j] = gpu->sim.ucellf[i][j];
#endif
        }
      }
    }
  }
  else {
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        if (gpu->sim.ntp > 0 && gpu->sim.barostat == 1) {
#ifdef use_DPFP
          M.map[i][j] = gpu->pbNTPData->_pSysData->recip[3*j + i];
#else
          M.map[i][j] = gpu->pbNTPData->_pSysData->recipf[3*j + i];
#endif
        }
        else {
#ifdef use_DPFP
          M.map[i][j] = gpu->sim.recip[i][j];
#else
          M.map[i][j] = gpu->sim.recipf[i][j];
#endif
        }
      }
    }
  }

  return M;
}

//---------------------------------------------------------------------------------------------
// AssignAtomToPencil: assign an atom to a honeycomb pencil hash cell basedon its coordinates.
//
// Arguments:
//   x,y,z:       atomic coordinates in the YZ plane (X does not matter for cell assignment)
//   uxsdata:     linearized form of the matrix for taking coordinates into fractional space
//   invuxsdata:  linearized form of the matrix taking fractional coordinates into real space
//   yshift:      orientation of the pencil cells
//   ny,nz:       the number of pencil cells in Y or Z
//   reset:       flag to have the coordinates of the atom reset, and returned, along with
//                the cell Y and Z result
//---------------------------------------------------------------------------------------------
static cellcrd AssignAtomToPencil(double x, double y, double z, dmat *U, dmat *invU,
                                  bool yshift, int ny, int nz, bool reset)
{
  double dny = (double)ny;
  double dnz = (double)nz;
  double invny = 1.0 / dny;
  double invnz = 1.0 / dnz;
  double fy = U->data[4]*y + U->data[7]*z;
  double fz =                U->data[8]*z;
  fy = fy - round(fy) + 0.5;
  fz = fz - round(fz) + 0.5;
  fy -= (double)(fy >= 1.0);
  fz -= (double)(fz >= 1.0);    
  double ry = invU->data[4]*fy + invU->data[5]*fz;
  double rz =                    invU->data[8]*fz;
  
  // Branch based on the orientation of the hash cells
  int besty, bestz;
  if (yshift) {
    
    // Harder case: the pencils are arranged in staggered fashion, with each stack of
    // them offset by 0.5 lengths in Y depending on their position relative to the
    // stack starting at the origin.
    int iz          = fz * dnz;
    int iy          = (int)(fy*dny - 0.5*(double)(iz & 0x1));
    double izf      = iz;
    double izp1f    = izf + 1.0;
    double iyfA     = (double)iy + 0.5*(iz & 0x1);
    double iyfB     = (double)iy + 0.5*((iz + 1) & 0x1);
    double iyp1fA   = iyfA + 1.0;
    double iyp1fB   = iyfB + 1.0;

    // The four candidates are now (iyfA, izf), (iyp1fA, izf), (iyfB, izfp1), (iyp1fB, izfp1).
    // Up to this point, they have been calculated as fractional coordinates expressed in the
    // units of the hash cell width (i.e. fractional coordinate 0.3333 in a grid of 18 cells
    // is 6).  Express these coordinates strictly in terms of the box fractional coordinates.
    iyfA           *= invny;
    iyfB           *= invny;
    iyp1fA         *= invny;
    iyp1fB         *= invny;
    izf            *= invnz;
    izp1f          *= invnz;
    double dry      = invU->data[4]*iyfA + invU->data[5]*izf - ry;
    double drz      =                      invU->data[8]*izf - rz;
    double rmin     = dry*dry + drz*drz;
    besty = iy;
    bestz = iz;
    dry             = invU->data[4]*iyp1fA + invU->data[5]*izf - ry;
    double rtest    = dry*dry + drz*drz;
    if (rtest < rmin) {
      rmin = rtest;
      besty = iy + 1;
    }
    dry   = invU->data[4]*iyfB + invU->data[5]*izp1f - ry;
    drz   =                      invU->data[8]*izp1f - rz;
    rtest = dry*dry + drz*drz;
    if (rtest < rmin) {
      rmin = rtest;
      besty = iy;
      bestz = iz + 1;
    }
    dry   = invU->data[4]*iyp1fB + invU->data[5]*izp1f - ry;
    rtest = dry*dry + drz*drz;
    if (rtest < rmin) {
      rmin = rtest;
      besty = iy + 1;
      bestz = iz + 1;
    }
    besty += (besty < 0) * ny;
  }
  else {

    // Easier case: the pencils are arranged on a regular grid (pencil rows don't stagger
    // along the Y-axis when viewed in fractional coordinate space).  Find the nearest
    // four, translate them into real space, and see which of them is truly closest.
    int iz          = fz * dnz;
    int iy          = fy * dny;
    double iyf      = iy;
    double izf      = iz;
    double iyp1f    = iyf + 1.0;
    double izp1f    = izf + 1.0;
    iyf             = (iyf   * invny * invU->data[4]) - ry;
    iyp1f           = (iyp1f * invny * invU->data[4]) - ry;
    izf            *= invnz;
    izp1f          *= invnz;
    double dry      = iyf + invU->data[5]*izf;
    double drz      =       invU->data[8]*izf - rz;
    double rmin     = dry*dry + drz*drz;
    besty = iy;
    bestz = iz;
    dry             = iyp1f + invU->data[5]*izf;
    double rtest    = dry*dry + drz*drz;
    if (rtest < rmin) {
      rmin = rtest;
      besty = iy + 1;
    }
    dry   = iyf + invU->data[5]*izp1f;
    drz   =       invU->data[8]*izp1f - rz;
    rtest = dry*dry + drz*drz;
    if (rtest < rmin) {
      rmin = rtest;
      besty = iy;
      bestz = iz + 1;
    }
    dry   = iyp1f + invU->data[5]*izp1f;
    rtest = dry*dry + drz*drz;
    if (rtest < rmin) {
      rmin = rtest;
      besty = iy + 1;
      bestz = iz + 1;
    }
  }
  
  // Ensure that the pencil lies in the primary grid
  besty -= ny * (besty == ny);
  bestz -= nz * (bestz == nz);
  cellcrd result;
  result.ycell = besty;
  result.zcell = bestz;

  // Calculate the new coordinates of the atom when seated properly within its designated
  // cell (this will shift the atom to be as close as possible to the cell's central axis,
  // while keeping its X coordinate within the primary unit cell)
  if (reset) {
    double fx = U->data[0]*x + U->data[3]*y + U->data[6]*z;
    fx = fx - round(fx) + 0.5;
    double fpy, fpz;
    if (yshift) {
      fpy = ((double)besty + 0.5*(double)(bestz & 0x1)) * invny;
    }
    else {
      fpy = (double)besty * invny;
    }
    fpz = (double)bestz * invnz;
    fy -= round(fy - fpy);
    fz -= round(fz - fpz);
    double rx = invU->data[0]*fx + invU->data[1]*fy + invU->data[2]*fz;
    ry        =                    invU->data[4]*fy + invU->data[5]*fz;
    rz        =                                       invU->data[8]*fz;    
    result.newx = rx;
    result.newy = ry;
    result.newz = rz;
  }

  return result;
}

//---------------------------------------------------------------------------------------------
// CheckCellAssignmentMismatch: given two cells that could contain the named atom, detect
//                              whether the assignment to one cell or the other could hinge on
//                              roundoff error.  Return true if this is the case, false
//                              otherwise.
//
// Arguments:
//
//   gpu:         overarching data structure containing simulation information, here used for
//                box parameters
//   crd[x,y,z]:  atomic coordinates for the atom of interest
//   c[1,2]idx:   indices of each cell candidate to hold the atom
//---------------------------------------------------------------------------------------------
static bool CheckCellAssignmentMismatch(gpuContext gpu, double *crdx, double *crdy,
                                        double *crdz, int c1idx, int c2idx)
{
  // Unpack the GPU context (this would be inefficient, but for the fact that this routine
  // should hardly ever be called)
  dmat U, invU;
  U = ExtractTransformationMatrix(gpu, false);
  invU = ExtractTransformationMatrix(gpu, true);
  bool yshift  = gpu->bPencilYShift;
  int ny       = gpu->sim.nypencils;
  int nz       = gpu->sim.nzpencils;
  double invny = 1.0 / (double)ny;
  double invnz = 1.0 / (double)nz;
  
  // Transform into fractional coordinates
  double fx = U.data[0]*(*crdx) + U.data[3]*(*crdy) + U.data[6]*(*crdz);
  double fy =                     U.data[4]*(*crdy) + U.data[7]*(*crdz);
  double fz =                                         U.data[8]*(*crdz);

  // Get the Y and Z coordinates of each point
  int c1z     = c1idx / ny;
  int c1y     = c1idx - c1z*ny;
  int c2z     = c2idx / ny;
  int c2y     = c2idx - c2z*ny;
  double fp1y, fp1z, fp2y, fp2z;
  if (yshift) {
    fp1y = ((double)c1y + 0.5*(double)(c1z & 0x1)) * invny;
    fp2y = ((double)c2y + 0.5*(double)(c2z & 0x1)) * invny;
  }
  else {
    fp1y = (double)c1y * invny;
    fp2y = (double)c2y * invny;
  }
  fp1z = (double)c1z * invnz;
  fp2z = (double)c2z * invnz;

  // Image the fractional coordinates with respect to each candidate cell center:
  // they become displacements which can then be returned to real space.
  double df1y = fy - fp1y;
  double df1z = fz - fp1z;
  df1y -= round(df1y);
  df1z -= round(df1z);
  double df2y = fy - fp2y;
  double df2z = fz - fp2z;
  df2y -= round(df2y);
  df2z -= round(df2z);

  // Return to real space, then compute the actual displacement from each center
  double rf1y = invU.data[4]*df1y + invU.data[5]*df1z;
  double rf1z =                     invU.data[8]*df1z;
  double rf2y = invU.data[4]*df2y + invU.data[5]*df2z;
  double rf2z =                     invU.data[8]*df2z;
  bool isClose = (fabs(sqrt(rf1y*rf1y + rf1z*rf1z) - sqrt(rf2y*rf2y + rf2z*rf2z)) < 1.0e-4);

  // If this is a close call, the second cell candidate is going to be chosen,
  // but the coordinates are assumed to be imaged with respect to the first cell.
  // Re-image the coordinates with respect to the second cell instead.
  if (isClose) {
    fy = fp2y + df2y;
    fz = fp2z + df2z;
    *crdx = invU.data[0]*fx + invU.data[1]*fy + invU.data[2]*fz;
    *crdy =                   invU.data[4]*fy + invU.data[5]*fz;
    *crdz =                                     invU.data[8]*fz;
  }
  
  // Free allocated memory
  DestroyDmat(&U);
  DestroyDmat(&invU);

  return isClose;
}

//---------------------------------------------------------------------------------------------
// GetSubImagePadding: compute parameters for padding the non-bonded hash cell atom layout.
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                object counts
//   atm_crd:     atomic coordinates
//   Uxs:         transformation matrix for coordinates in the YZ plane, taking things into
//                fractional space
//   invUxs:      transformation matrix for coordinates in the YZ plane, taking things from
//                factional space into real space
//   ny, nz:      the Y or Z dimensions of the primary image hash cell table
//   yshift:      flag to indicate the orientation of the pencil hash cell grid (is the grid
//                staggered along the Y axis for different levels of Z?)
//   xlen:        the length of the unit cell along the X dimension
//   es_cutoff:   electrostatic cutoff non-bonded cutoff
//   lj_cutoff:   Lennard-Jones non-bonded cutoff (must be greater than or equal to es_cutoff)
//   nbskin:      the non-bonded buffer distance for pair list rebuilds
//   tNL:         collection of information about interactions between hash cells in the
//                neighbor list
//---------------------------------------------------------------------------------------------
static void GetSubImagePadding(gpuContext gpu, double atm_crd[][3], dmat *Uxs, dmat *invUxs,
                               int ny, int nz, bool yshift, double es_cutoff, double lj_cutoff,
                               double nbskin, nlkit *tNL)
{
  int i, j;
  double *uxsdata, *invuxsdata;
  dmat qqpop, ljpop, U, invU;

  // Allocate space to store hash cell populations in the initial snapshot
  qqpop = CreateDmat(ny, nz);
  ljpop = CreateDmat(ny, nz);

  // Loop over all atoms and accumulate cell populations
  U = CreateDmat(3, 3);
  invU = CreateDmat(3, 3);
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      U.map[i][j] = gpu->sim.recip[i][j];
      invU.map[i][j] = gpu->sim.ucell[i][j];
    }
  }
  for (i = 0; i < gpu->sim.atoms; i++) {

    // Compute fractional coordinates and re-imaged real coordinates
    cellcrd optcrd;
    optcrd = AssignAtomToPencil(atm_crd[i][0], atm_crd[i][1], atm_crd[i][2],
                                &U, &invU, yshift, ny, nz, false);
    unsigned int ljt = gpu->pbAtomLJID->_pSysData[i] - 1;
    int islj = ((gpu->sim.ljexist[ljt >> 5] & (0x1 << (ljt & 31))) > 0);
    int isqq = (fabs(gpu->pbAtomChargeSP->_pSysData[i]) > 1.0e-8);
    qqpop.map[optcrd.ycell][optcrd.zcell] += isqq;
    ljpop.map[optcrd.ycell][optcrd.zcell] += islj;
  }

  // Estimate the maximum number of atoms that any cell would have to be
  // padded with in order to hold atom images needed along the x direction.
  double xqq = gpu->sim.ucell[0][0] / (es_cutoff + nbskin);
  double xlj = gpu->sim.ucell[0][0] / (lj_cutoff + nbskin);

  // Find the maximum population.  Find the standard deviation.  Calculate the
  // largest number of atoms that any cell could reasonably be expected to hold
  // if these conditions persist.
  int nqqmax = DVecExtreme(qqpop.data, ny*nz);
  int nqqpad = ((int)(nqqmax / xqq) + 1) * 4;
  int nqqstd = DStDev(qqpop.data, ny*nz);
  nqqmax = (nqqmax*3/2 < nqqmax + nqqstd*10) ? nqqmax + nqqstd*10 : nqqmax * 3 / 2;
  int nljmax = DVecExtreme(ljpop.data, ny*nz);
  int nljpad = ((int)(nljmax / xlj) + 1) * 4;
  int nljstd = DStDev(ljpop.data, ny*nz);
  nljmax = (nljmax*3/2 < nljmax + nljstd*10) ? nljmax + nljstd*10 : nljmax * 3 / 2;

  // Set values in the gpuContext: these will be carried back up.  The Y-dimension padding
  // is decreased by 1 if the box is roughly orthorhombic in the YZ plane.
  int ypadding = 4;
  int zpadding = 4;
  if (yshift) {
    ypadding -= 1;
  }
  gpu->sim.ypadding = ypadding;
  gpu->sim.zpadding = zpadding;
  gpu->sim.QQCapPadding = nqqpad;
  gpu->sim.LJCapPadding = nljpad;
  gpu->sim.QQCellSpace = nqqmax + 2*nqqpad;
  gpu->sim.LJCellSpace = nljmax + 2*nljpad;
  gpu->sim.QQRowSpace = gpu->sim.QQCellSpace * (2*ypadding + ny);
  gpu->sim.LJRowSpace = gpu->sim.LJCellSpace * (2*ypadding + ny);
  gpu->sim.QQBaseOffset = (gpu->sim.QQRowSpace * zpadding) +
                          (gpu->sim.QQCellSpace * ypadding) + nqqpad;
  gpu->sim.LJBaseOffset = (gpu->sim.LJRowSpace * zpadding) +
                          (gpu->sim.LJCellSpace * ypadding) + nljpad;
  tNL->TotalQQSpace = (2*ypadding + ny) * (zpadding + nz) * gpu->sim.QQCellSpace;
  tNL->TotalLJSpace = (2*ypadding + ny) * (zpadding + nz) * gpu->sim.LJCellSpace;
  gpu->sim.ydimFull = 2*ypadding + ny;
  gpu->sim.zdimFull =   zpadding + nz;
  gpu->sim.nExpandedPencils = (2*ypadding + ny) * (zpadding + nz);

  // Determine the padding factors for the X dimension (needed for pre-imaging atoms)
  double capbound = ComputeCapBoundary(lj_cutoff, nbskin, &invU);
  gpu->sim.LowerCapBoundary = capbound;
  gpu->sim.UpperCapBoundary = 1.0 - capbound;

  // Free allocated memory
  DestroyDmat(&U);
  DestroyDmat(&invU);
  DestroyDmat(&qqpop);
  DestroyDmat(&ljpop);
}

//---------------------------------------------------------------------------------------------
// CullFlatPencilDesigns: function to mark pencils with abnormally flat or oblong
//                        cross-sections so that they will not be chosen for the final pair
//                        list building design.
//
// Arguments:
//   vertices:   the vertices of the pencil hash cell cross section
//   rstd:       the best (minimum) ratio of the longest to shortest distances between two
//               opposing vertices of the pencil cross section found thus far.  Only pencils
//               whose ratio falls within 15% of this cross section will be allowed.
//---------------------------------------------------------------------------------------------
static bool CullFlatPencilDesigns(dmat *vertices, double *rstd, double *rcurr)
{
  int i, j;

  // Find the aspect ratio
  double rval[3];
  for (i = 0; i < 3; i++) {
    double dy = vertices->map[i + 2][0] - vertices->map[i][0];
    double dz = vertices->map[i + 2][1] - vertices->map[i][1];
    rval[i] = sqrt(dy*dy + dz*dz);
  }
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      if (rval[j] > rval[j+1]) {
        double tmp  = rval[j];
        rval[j]     = rval[j + 1];
        rval[j + 1] = tmp;
      }
    }
  }
  double rnow = rval[2] / rval[0];
  if (*rstd < 0.0 || rnow < *rstd) {
    *rstd = rnow;
  }
  *rcurr = rnow;
  
  return (rnow >= 1.15 * (*rstd));
}

//---------------------------------------------------------------------------------------------
// CullStencilViolations: return true if the stated combination of pencils will violate the
//                        Neutral-Territory stencil
//---------------------------------------------------------------------------------------------
static bool CullStencilViolations(int k, int m, bool yshift, double separation, double maxsep)
{
  if (separation >= maxsep) {
    return false;
  }

  // Things that will always cause a problem
  if (k < -4 || m < -4) {
    return true;
  }
  
  // Trickier case: Y shifting for odd-numbered rows
  if (yshift) {
    if (k == -1 && m >= 4) {
      return true;
    }
    else if (k == -3 && (m >= 3 || m < -3)) {
      return true;
    }
    else if ((k == -2 || k == -4) && abs(m) > 4 + k/2) {
      return true;
    }
  }

  // Easier case: most problems have already been addressed if k >= -4,
  else {
    if (m > 4 + k) {
      return true;
    }
  }

  return false;
}

//---------------------------------------------------------------------------------------------
// GetHoneycombDesign: decide on the number of pencils (cells) to fill the simulation box and
//                     return the parameters for making the domain decomposition.
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                object counts
//   atm_crd:     atomic coordinates
//   es_cutoff:   electrostatic non-bonded cutoff
//   lj_cutoff:   Lennard-Jones non-bonded cutoff
//   nbskin:      the non-bonded pair list buffer distance
//---------------------------------------------------------------------------------------------
nlkit GetHoneycombDesign(gpuContext gpu, double atm_crd[][3], double es_cutoff,
                         double lj_cutoff, double nbskin)
{
  int i, j, k, m;

  // Determine the optimal count of pencils in each direction, then the degree
  // to which the box can shrink and still have this layout be relevant.  The
  // optimal pencil count is the one that minimizes the extra volume that must
  // be counted in order to cover everything.
  dmat U, invU, Uxs, invUxs;
  U = CreateDmat(3, 3);
  invU = CreateDmat(3, 3);
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      U.map[i][j] = gpu->sim.recip[i][j];
      invU.map[i][j] = gpu->sim.ucell[i][j];
    }
  }
  Uxs = CreateDmat(2, 2);
  invUxs = CreateDmat(2, 2);
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      Uxs.map[i][j] = U.map[i+1][j+1];
      invUxs.map[i][j] = invU.map[i+1][j+1];
    }
  }

  // Set yshift to true if the pencils will be staggered along the Y axis.
  double mvecy = sqrt((invUxs.map[0][0] * invUxs.map[0][0]) +
                      (invUxs.map[1][0] * invUxs.map[1][0]));
  double mvecz = sqrt((invUxs.map[0][1] * invUxs.map[0][1]) +
                      (invUxs.map[1][1] * invUxs.map[1][1]));
  double angle = acos(((invUxs.map[0][0] * invUxs.map[0][1]) +
                       (invUxs.map[1][0] * invUxs.map[1][1])) / (mvecy * mvecz));
  bool yshift = (angle < 7.0 * PI / 12.0);
  double boxv = invU.map[0][0] * invU.map[1][1] * invU.map[2][2];

  // Compute convenient quantities
  es_cutoff += nbskin;
  lj_cutoff += nbskin;
  double es_cutoff2 = es_cutoff * es_cutoff;
  double lj_cutoff2 = lj_cutoff * lj_cutoff;
  
  // Allocate space to record the pattern of pencil-to-pencil interactions
  imat qqPencils, ljPencils, bestqqPencils, bestljPencils;
  dmat qqRanges, ljRanges, bestqqRanges, bestljRanges;
  int nmaxp = 128;
  qqPencils = CreateImat(nmaxp, 2);
  ljPencils = CreateImat(nmaxp, 2);
  qqRanges = CreateDmat(1, nmaxp);
  ljRanges = CreateDmat(1, nmaxp);
  bestqqPencils = CreateImat(nmaxp, 2);
  bestljPencils = CreateImat(nmaxp, 2);
  bestqqRanges = CreateDmat(1, nmaxp);
  bestljRanges = CreateDmat(1, nmaxp);

  // Loop over different numbers of pencils to find the best fit.  The objective is to
  // fit everything in a stencil that is nine honeycomb cells wide--analogous, in fact,
  // to the 64-element Hilbert space filling curve that Scott used.  The 64-element
  // curve makes a grid 4 segments on a side, each of them as wide as one quarter of
  // the combined cutoff plus pair list margin.  Similarly, the honeycomb cells must be
  // at least as wide as one quarter (well, 7/2) the cutoff plus pair list margin.
  // This will make a list of up to 31 interactions between pencil cells that comprise
  // all non-bonded interactions.
  //
  //     O @ @ @ X             In the drawing at left, X is the home pencil, @
  //      @ @ @ @ @ @ @ @      represents a pencil with which it will certainly
  //       O @ @ @ @ @ O       interact, and O represents a pencil with which it
  //        @ @ @ @ @ @        may interact, depending on the exact relationship
  //         O @ O @ O         between the cutoff and box size.
  //
  int mintile = 0;
  int minPencilQQ = 0;
  int minPencilLJ = 0;
  int besti = 0;
  int bestj = 0;
  int imax = (invU.map[0][1] + invU.map[1][1]) / (lj_cutoff / 3.5);
  int jmax = invU.map[2][2] / (lj_cutoff / 3.5);
  int imin = imax - 3;
  int jmin = jmax - 3;
  imin = (imin < 7) ? 7 : imin;
  jmin = (jmin < 7) ? 7 : jmin;
  imax = (imax < imin) ? imin : imax;
  jmax = (jmax < jmin) ? jmin : jmax;
  double rstd = -1.0;
  double bestratio = -1.0;
  for (i = imin; i <= imax; i++) {
    double dbli = (double)i;
    for (j = jmin; j <= jmax; j++) {
      
      // The Z cell count must be a multiple of two if the pencils are staggered
      // along the Y axis.  The entire arrangement must tile and staggered
      // stacking implies a two-cell asymmetric unit to the packing.
      if (yshift && (j % 2) == 1) {
        continue;
      }
      double dblj = (double)j;

      // Find the length of the hexagonal prism inscribed in each pencil
      // That would be needed to hold four atoms, on average.
      double nQQatpen = (double)(gpu->sim.nQQatoms) / (dbli * dblj);
      double nLJatpen = (double)(gpu->sim.nLJatoms) / (dbli * dblj);
      double lsectQQ = invU.map[0][0] * (4.0 / nQQatpen);
      double lsectLJ = invU.map[0][0] * (4.0 / nLJatpen);
      int nQQsect = (int)((invU.map[0][0] / lsectQQ) + 0.999999);
      int nLJsect = (int)((invU.map[0][0] / lsectLJ) + 0.999999);

      // Compute the maximum distance from the center for each
      // hexagonal pencil, and get the locations of vertices
      // for the pencil at Y = Z = 0.
      dmat vertices;
      vertices = GetPencilVertices(i, j, yshift, &invUxs);

      // Prefer combinations with roughly round pencils
      double rcurr;
      if (CullFlatPencilDesigns(&vertices, &rstd, &rcurr)) {
        continue;
      }      
      double sqmaxdist = 0.0;
      for (k = 0; k < 6; k++) {
        double ry = vertices.map[k][0];
        double rz = vertices.map[k][1];
        if (ry*ry + rz*rz > sqmaxdist) {
          sqmaxdist = ry*ry + rz*rz;
        }
      }
      double sqmaxdistLJ = 2.0*sqrt(sqmaxdist) + lj_cutoff;
      sqmaxdistLJ *= sqmaxdistLJ;
      
      // Count the number of pencil-to-pencil interactions that will be
      // needed in electrostatic as well as Lennard-Jones contexts.
      int ntile = 0;
      int npencilQQ = 0;
      int npencilLJ = 0;
      for (k = -j; k <= 0; k++) {
        int mlim = (k == 0) ? 0 : i;
        for (m = -i; m <= mlim; m++) {
          double dy = (double)m / dbli;
          double dz = (double)k / dblj;
          if (yshift) {
            dy += 0.5 * (double)((abs(k) % 2) == 1) / dbli;
          }
          double ry = invUxs.map[0][0]*dy + invUxs.map[0][1]*dz;
          double rz =                       invUxs.map[1][1]*dz;
          double r2 = ry*ry + rz*rz;
          
          // Quick check that the two cells could possibly have atoms
          // within range of one another (lj_cutoff >= es_cutoff)
          if (r2 > sqmaxdistLJ) {
            continue;
          }

          // Check whether this cell:cell interaction will definitely
          // have electrostatic or Lennard-Jones tests to perform.
          double separation = CheckAllPencilVertices(ry, rz, lj_cutoff, &vertices);

          // Check that the two cells are indeed sufficiently spaced if
          // they fall outside the intended stencil
          if (CullStencilViolations(k, m, yshift, separation, lj_cutoff)) {
            continue;
          }

          // This tile configuration passes muster, so estimate how expensive it will be.
          int ntileQQ = 0;
          int ntileLJ = 0;
          if (separation < es_cutoff) {
            double esOverhang = sqrt(es_cutoff2 - (separation * separation));
            double esBreadth = 2.0*esOverhang;
            ntileQQ = (int)((esBreadth / lsectQQ) + 0.999999);
            qqPencils.map[npencilQQ][0] = m;
            qqPencils.map[npencilQQ][1] = k;
            qqRanges.map[0][npencilQQ] = separation;
          }
          if (separation < lj_cutoff) {
            double ljOverhang = sqrt(lj_cutoff2 - (separation * separation));
            double ljBreadth = 2.0*ljOverhang;
            ntileLJ = (int)((ljBreadth / lsectLJ) + 0.999999);            
            ljPencils.map[npencilLJ][0] = m;
            ljPencils.map[npencilLJ][1] = k;
            ljRanges.map[0][npencilLJ] = separation;
          }
          npencilQQ += (ntileQQ > 0);
          npencilLJ += (ntileLJ > 0);
          if (npencilQQ >= nmaxp || npencilLJ >= nmaxp) {
            nmaxp *= 2;
            qqPencils = ReallocImat(&qqPencils, nmaxp, 2);
            ljPencils = ReallocImat(&ljPencils, nmaxp, 2);
            qqRanges = ReallocDmat(&qqRanges, 1, nmaxp);
            ljRanges = ReallocDmat(&ljRanges, 1, nmaxp);
            bestqqPencils = ReallocImat(&bestqqPencils, nmaxp, 2);
            bestljPencils = ReallocImat(&bestljPencils, nmaxp, 2);
            bestqqRanges = ReallocDmat(&bestqqRanges, nmaxp, 2);
            bestljRanges = ReallocDmat(&bestljRanges, nmaxp, 2);
          }
          ntileQQ *= nQQsect;
          ntileLJ *= nLJsect;
          ntile += ntileQQ + ntileLJ;
        }
      }
      ntile *= i * j;
      if (besti == 0 || ntile < mintile) {
        besti = i;
        bestj = j;
        bestratio = rcurr;
        mintile = ntile;
        minPencilQQ = npencilQQ;
        minPencilLJ = npencilLJ;
        for (k = 0; k < npencilQQ; k++) {
          bestqqPencils.map[k][0] = qqPencils.map[k][0];
          bestqqPencils.map[k][1] = qqPencils.map[k][1];
          bestqqRanges.map[0][k]  = qqRanges.map[0][k];
        }
        for (k = 0; k < npencilLJ; k++) {
          bestljPencils.map[k][0] = ljPencils.map[k][0];
          bestljPencils.map[k][1] = ljPencils.map[k][1];
          bestljRanges.map[0][k]  = ljRanges.map[0][k];
        }
      }

      // Free allocated memory
      DestroyDmat(&vertices);
    }
  }
  
  // Construct the neighbor list kit.  Lay out the honeycomb cell (pencil) centers,
  // then record the map of pencil-to-pencil interactions for future use.
  nlkit tNL;
  tNL.yshift = yshift;
  tNL.ny = besti;
  tNL.nz = bestj;
  tNL.qqPencils = CreateImat(minPencilQQ, 2);
  tNL.ljPencils = CreateImat(minPencilLJ, 2);
  tNL.qqRanges = CreateDmat(1, minPencilQQ);
  tNL.ljRanges = CreateDmat(1, minPencilLJ);
  for (i = 0; i < minPencilQQ; i++) {
    tNL.qqPencils.map[i][0] = bestqqPencils.map[i][0];
    tNL.qqPencils.map[i][1] = bestqqPencils.map[i][1];
    tNL.qqRanges.map[0][i]  = bestqqRanges.map[0][i];
  }
  for (i = 0; i < minPencilLJ; i++) {
    tNL.ljPencils.map[i][0] = bestljPencils.map[i][0];
    tNL.ljPencils.map[i][1] = bestljPencils.map[i][1];
    tNL.ljRanges.map[0][i]  = bestljRanges.map[0][i];
  }
  double qqDens = (double)gpu->sim.nQQatoms / (invU.data[0] * invU.data[4] * invU.data[8]);
  double ljDens = (double)gpu->sim.nLJatoms / (invU.data[0] * invU.data[4] * invU.data[8]);
  double qqPairEst = (2.0 / 3.0) * PI * es_cutoff * es_cutoff * es_cutoff * qqDens;
  double ljPairEst = (2.0 / 3.0) * PI * lj_cutoff * lj_cutoff * lj_cutoff * ljDens;  
  qqPairEst *= (double)gpu->sim.nQQatoms;
  ljPairEst *= (double)gpu->sim.nLJatoms;
  tNL.TotalPairSpace  = 2.0 * (qqPairEst + ljPairEst);
  tNL.TotalPairSpace  = (tNL.TotalPairSpace + HCMB_PL_STORAGE_SIZE - 1);
  tNL.TotalPairSpace  = (tNL.TotalPairSpace / HCMB_PL_STORAGE_SIZE) + (2 * gpu->blocks);
  tNL.TotalPairSpace *= HCMB_PL_STORAGE_SIZE;
  GetSubImagePadding(gpu, atm_crd, &Uxs, &invUxs, besti, bestj, yshift, es_cutoff, lj_cutoff,
                     nbskin, &tNL);
  
  // Free allocated memory
  DestroyImat(&bestqqPencils);
  DestroyImat(&bestljPencils);
  DestroyImat(&qqPencils);
  DestroyImat(&ljPencils);
  DestroyDmat(&bestqqRanges);
  DestroyDmat(&bestljRanges);
  DestroyDmat(&qqRanges);
  DestroyDmat(&ljRanges);
  DestroyDmat(&U);
  DestroyDmat(&invU);
  DestroyDmat(&Uxs);
  DestroyDmat(&invUxs);

  return tNL;
}

//---------------------------------------------------------------------------------------------
// DestroyNLKit: destroy a neighbor list kit put together (allocated) by the above function.
//
// Arguments:
//   tNL:      the neighbor list kit that's no longer needed
//---------------------------------------------------------------------------------------------
void DestroyNLKit(nlkit *tNL)
{
  DestroyImat(&tNL->qqPencils);
  DestroyImat(&tNL->ljPencils);
  DestroyDmat(&tNL->qqRanges);
  DestroyDmat(&tNL->ljRanges);
}

//---------------------------------------------------------------------------------------------
// PackPencilRelationship: pack four integers describing the honeycomb pencil-to-pencil hash
//                         cell comparison that must be performed in order to cover all
//                         non-bonded interactions in a neutral territory relationship.  The
//                         hard-wired assumption is that the import region will be eleven
//                         pencils with up to 1020 atoms in all.
//
// Arguments:
//   gpu:       overarching data structure containing simulation information, including the
//              arrays to fill and non-bonded criteria of the simulation
//   P:         matrix containing all pencil relationships, (0, 0) to (y, z)
//   [qq,lj]R:  matrix containing all pencil to pencil minimum separations (ranges)
//   span:      the breadth of the Neutral Territory region to map
//   ntPencils: the list of pencils to import in order to fulfill the Neutral Territory
//              relationships (indexed relative ot the home cell, returned)
//---------------------------------------------------------------------------------------------
void PackPencilRelationship(gpuContext gpu, imat *P, dmat *R, int span, imat *ntPencils,
                            imat *ntInteractions)
{
  int i, j, k;

  imat occupancy, candidates, assigned;
  *ntPencils = CreateImat(3*span - 1, 2);
  *ntInteractions = CreateImat(P->row, 2);
  occupancy  = CreateImat(1, 3*span - 1);
  candidates = CreateImat(2, 3*span - 1);
  assigned   = CreateImat(1, P->row);
  
  if (gpu->bPencilYShift) {

    // Harder case: the imports are as in the case below, except every odd row is
    // shifted forward half a space (cell (0, 1) is one half space to the right of
    // cell (0, 0) and one space above it).  By construction, all the relationships
    // coming in via P start at (0, 0).  When moving these relationships to begin
    // on odd numbred rows, the offsets must again be considered--add one to the
    // necessary delta in Y.  If the central cell itself is located on an odd row,
    // then further adjustments need to happen.
    for (i = 0; i < span; i++) {
      ntPencils->map[i][0] = (i + 1) / 2;
      ntPencils->map[i][1] = -(i + 1);
      ntPencils->map[span + i][0] = -i;
      if (i > 0) {
        ntPencils->map[2*span + i - 1][0] = -(i + 2) / 2;
        ntPencils->map[2*span + i - 1][1] = -(i + 1);
      }
    }
  }
  else {

    // Easier case: the imported atoms will come from the original pencil,
    // then as many pencils as needed (N, the number needed) in the -Z
    // direction, then N-1 pencils in the -Y direction, and finally N-1
    // pencils along the -Y, -Z diagonal, skipping (-1, -1).
    for (i = 0; i < span; i++) {
      ntPencils->map[i][1] = -(i + 1);
      ntPencils->map[span + i][0] = -i;
      if (i > 0) {
        ntPencils->map[2*span + i - 1][0] = -(i + 1);
        ntPencils->map[2*span + i - 1][1] = -(i + 1);
      }
    }
  }
  
  // Search the list of neutral-territory pencils for a match to each comparison in P
  for (i = 0; i < P->row; i++) {

    // First find the comparison with the largest displacement
    // (this will probably have the least options for adding to
    // the neutral-territory decomposition)
    int designate = -1;
    int disp = -1;
    for (j = 0; j < P->row; j++) {
      if (assigned.data[j] == 0 && abs(P->map[j][0]) + abs(P->map[j][1]) > disp) {
        disp = abs(P->map[j][0]) + abs(P->map[j][1]);
        designate = j;
      }
    }
    assigned.data[designate] = 1;
    
    // Compute candidates
    SetIVec(candidates.data, 2 * ntPencils->row, 0);
    for (j = 0; j < ntPencils->row; j++) {
      for (k = 0; k < ntPencils->row; k++) {
        if (gpu->bPencilYShift) {
          if ((P->map[designate][1] & 1) == 0) {

            // If the relationship calls for two cells spaced by an even number of rows,
            // all of the alignments will take care of themselves.
            if (ntPencils->map[k][0] - ntPencils->map[j][0] == P->map[designate][0] &&
                ntPencils->map[k][1] - ntPencils->map[j][1] == P->map[designate][1]) {
              candidates.map[0][k] = 1;
              candidates.map[1][k] = j;
            }
          }
          else {

            // If the starting cell for this candidate is on an odd-numbered row
            // and the relationship spans an odd number of rows, then the origin
            // point needs to be promoted one space forward before checking the
            // relative positioning.  If, however, the putative neutral-territory
            // relationship starts on an even-numbered row, then there's nothing
            // to worry about.
            if (ntPencils->map[k][0] -
                (ntPencils->map[j][0] + (ntPencils->map[j][1] & 1)) == P->map[designate][0] &&
                ntPencils->map[k][1] - ntPencils->map[j][1] == P->map[designate][1]) {
              candidates.map[0][k] = 1;
              candidates.map[1][k] = j;
            }
          }
        }
        else {

          // No Y-shifting?  No problems.
          if (ntPencils->map[k][0] - ntPencils->map[j][0] == P->map[designate][0] &&
              ntPencils->map[k][1] - ntPencils->map[j][1] == P->map[designate][1]) {
            candidates.map[0][k] = 1;
            candidates.map[1][k] = j;
          }
        }
      }
    }

    // Select the best candidate
    int bestcand = -1;
    int bestocc = 0; 
    for (j = 0; j < ntPencils->row; j++) {
      if (candidates.map[0][j] == 1) {
        if (bestcand == -1 || occupancy.data[j] < bestocc) {
          bestcand = j;
          bestocc = occupancy.data[j];
        }
      }
    }
    if (bestcand == -1) {
      printf("PackPencilRelationship :: Error.  Failed to find a candidate for "
             "[ %2d %2d %2d %2d ]\n", 0, 0, P->map[designate][0], P->map[designate][1]);
      exit(1);
    }

    // Place the relationship, extending from candidates.map[1][bestcand]
    // to bestcand in the list of neutral-territory pencils
    occupancy.data[bestcand] += 1;
    int starty = ntPencils->map[candidates.map[1][bestcand]][0];
    int startz = ntPencils->map[candidates.map[1][bestcand]][1];
    int endy   = ntPencils->map[bestcand][0];
    int endz   = ntPencils->map[bestcand][1];
    gpu->pbPencils->_pSysData[i] = ((endy + 128) << 24) | ((endz + 128) << 16) |
                                   ((starty + 128) << 8) | (startz + 128);

    // Remember the interaction in terms of which imported pencils interact
    ntInteractions->map[i][0] = candidates.map[1][bestcand];
    ntInteractions->map[i][1] = bestcand;
    
    // Take the range into a 32-bit floating point representation expressible as an int
    double R0 = R->map[0][designate];
    double es2 = gpu->sim.esCutPlusSkin2 - R0*R0;
    gpu->pbPencilRanges->_pSysData[i     ] = (es2 > 1.0e-6) ? sqrt(es2) : -999.0;
    gpu->pbPencilRanges->_pSysData[i + 32] = sqrt(gpu->sim.ljCutPlusSkin2 - R0*R0);
  }

  // Free allocated memory
  DestroyImat(&occupancy);
  DestroyImat(&candidates);
  DestroyImat(&assigned);
}

//---------------------------------------------------------------------------------------------
// MakeSegmentRelays: make the array of segment relays based on the pencil relationships that
//                    were established in the PackPencilRelationship function above.
//
// Arguments:
//   gpu:             overarching data structure containing simulation information, including
//                    the Neutral Territory segment interactions that will be loaded
//   ntInteractions:  the CPU-computed Neutral Territory interaction list
//---------------------------------------------------------------------------------------------
void MakeSegmentRelays(gpuContext gpu, imat *ntInteractions)
{
  int i, j;

  // Initialize the array to reasonable but meaningless values (this will not protect against
  // bad memory accesses, but it will immediately indicate if arrays are overrun)
  for (i = 0; i < 128; i++) {
    gpu->pbSegmentRelays->_pSysData[i] = -1;
  }
  for (i = 0; i < ntInteractions->row; i++) {
    
    // Double-counting due to atom in the same pencil hash cell interacting
    // can be partly dealt with right here, by recognizing that atoms of the
    // second segment will necessarily come after atoms of the first. 
    if (ntInteractions->map[i][0] != ntInteractions->map[i][1]) {
      gpu->pbSegmentRelays->_pSysData[i     ] = ntInteractions->map[i][0];
      gpu->pbSegmentRelays->_pSysData[i + 32] = ntInteractions->map[i][1] + 11;
    }
    gpu->pbSegmentRelays->_pSysData[i + 64] = ntInteractions->map[i][0] + 11;
    gpu->pbSegmentRelays->_pSysData[i + 96] = ntInteractions->map[i][1];
  }
}

//---------------------------------------------------------------------------------------------
// SortInt2ByX: function for qsort to use when sorting int2 types by their X components into
//              ascending order.
//
// Arguments:
//   dual[A, B]:   the int2 vectors to compare
//---------------------------------------------------------------------------------------------
int SortInt2ByX(const void *dualA, const void *dualB)
{
  int intA = ((int2*)dualA)[0].x;
  int intB = ((int2*)dualB)[0].x;

  if (intA < intB) {
    return -1;
  }
  else if (intA > intB) {
    return 1;
  }
  else {
    return 0;
  }
}

//---------------------------------------------------------------------------------------------
// FindSubImgOrderingPartition: determine the best way to split kOrderSubImagePencils_kernel
//                              between blocks that will start a prefix sum over atoms in
//                              hash cells of electrostatic and Lennard-Jones sub-images and
//                              blocks that will update the primary image filtering.
//
// Arguments:
//   natom:      the number of atoms in the system
//   nsmp:       the number of streaming multi-processors on the GPU
//   maxcell:    the number of hash cells in the electrostatic sub-image
//---------------------------------------------------------------------------------------------
int FindSubImgOrderingPartition(int natom, int nsmp, int maxcell)
{
  // Default approach: give blocks to prefix sum computation such that each thread
  // will only have to run through the loop once.
  int nprfxblk = (maxcell + 255) >> 8;
  if (nprfxblk > 8) {
    nprfxblk = 8;
  }
  int bestblk = 0;
  int bestoption = 0;
  while (nprfxblk > 0) {

    // Valid for SM_3X and higher, but may need adjustment if
    // future architectures take more than 1024 threads per block
    int totalblk = nsmp / 256;
    int nfltrblk = totalblk - 2*nprfxblk;

    // Estimate that the prefix sum computation is about
    // three times as expensive as the filtering update.
    int prfxwork = 3 * ((maxcell + (256*nprfxblk - 1)) / (256 * nprfxblk));
    int fltrwork = (natom + 256*nfltrblk - 1) / (256 * nfltrblk);
    int worstcase = (prfxwork > fltrwork) ? prfxwork : fltrwork;
    if (bestblk == 0 || worstcase < bestoption) {
      bestoption = worstcase;
      bestblk = nprfxblk;
    }
    nprfxblk--;
  }

  return bestblk;
}

//---------------------------------------------------------------------------------------------
// FindExpansionBatchCount: find the optimal number of consecutive hash cells to assign to each
//                          thread block while laying out sub image expansion instructions.
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                object counts
//---------------------------------------------------------------------------------------------
int FindExpansionBatchCount(gpuContext gpu)
{
  int batchSize;

  // Calculate the average number of cells that each SM will need
  // to deal with (each SM is loaded with 1024 threads).
  batchSize = (gpu->sim.nExpandedPencils + (gpu->blocks - 1)) / gpu->blocks;
  if (batchSize > 32) {
    batchSize = 32;
  }

  return batchSize;
}

//---------------------------------------------------------------------------------------------
// CpuCalculateCellID: calculate cell indices for all atoms on the CPU, for comparison to GPU
//                     results.
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                object counts
//   cellID:      cell assignments computed on the CPU (should match HcmbCellID downloaded
//                from the GPU)
//   primaryCrd:  array to hold coordinates of the primary image, with pencil cell adjustments
//                that put each atom nearest the central axis of their appropriate pencil
//                cells in Y and Z (this may put them outside the primary unit cell by up to
//                half the thickness of a pencil) with X strictly inside the primary unit cell
//---------------------------------------------------------------------------------------------
static void CpuCalculateCellID(gpuContext gpu, int* cellID, double* primaryCrd)
{
  int i, j;

  // Extract the transformation matrices
  dmat U, invU;
  U = ExtractTransformationMatrix(gpu, false);
  invU = ExtractTransformationMatrix(gpu, true);

  // Loop over all atoms
  bool yshift = gpu->bPencilYShift;
  unsigned int *gpuCellID;
  double *xcrd, *ycrd, *zcrd;
  if (gpu->sim.pHcmbCellID == gpu->pbHcmbCellID->_pDevData) {
    gpuCellID = gpu->pbHcmbCellID->_pSysData;
    xcrd = gpu->pbHcmb->_pSysData;
    ycrd = gpu->pbHcmb->_pSysData + gpu->sim.stride;
    zcrd = gpu->pbHcmb->_pSysData + gpu->sim.stride2;
  }
  else {
    gpuCellID = gpu->pbHcmbCellID->_pSysData + gpu->sim.stride;
    xcrd = gpu->pbHcmb->_pSysData + gpu->sim.stride3;
    ycrd = gpu->pbHcmb->_pSysData + gpu->sim.stride4;
    zcrd = gpu->pbHcmb->_pSysData + (gpu->sim.stride * 5);
  }
  for (i = 0; i < gpu->sim.atoms; i++) {
    cellcrd optcrd;
    optcrd = AssignAtomToPencil(xcrd[i], ycrd[i], zcrd[i], &U, &invU, yshift,
                                gpu->sim.nypencils, gpu->sim.nzpencils, true);
    cellID[i] = (optcrd.zcell * gpu->sim.nypencils) + optcrd.ycell;
    primaryCrd[3*i    ] = optcrd.newx;
    primaryCrd[3*i + 1] = optcrd.newy;
    primaryCrd[3*i + 2] = optcrd.newz;

    // There may be a disagreement here if single precision arithmetic is in effect.
    // Detect mismatches due to roundoff and correct them.
    if (cellID[i] != (gpuCellID[i] >> 2) &&
        CheckCellAssignmentMismatch(gpu, &optcrd.newx, &optcrd.newy, &optcrd.newz, cellID[i],
                                    (gpuCellID[i] >> 2))) {
      cellID[i] = (gpuCellID[i] >> 2);
    }
    
    // Raise an error message if any cell assignments truly disagree
    if (cellID[i] != (gpuCellID[i] >> 2)) {
      int gpucellz = (gpuCellID[i] >> 2) / gpu->sim.nypencils;
      int gpucelly = (gpuCellID[i] >> 2) - (gpucellz * gpu->sim.nypencils);
      printf("CpuCalculateCellID :: Mismatch in atom %7d, CPU calculated cell %5d [ %3d %3d ] "
             "versus GPU %7d [ %3d %3d ]\n", i, cellID[i], optcrd.ycell, optcrd.zcell,
             gpuCellID[i] >> 2, gpucelly, gpucellz);
    }
  }

  // Free allocated memory
  DestroyDmat(&U);
  DestroyDmat(&invU);
}

//---------------------------------------------------------------------------------------------
// CpuCalculateCellLimits: compute the limits of cells in the primary image that will later be
//                         needed to expand the image
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                object counts
//   cellID:      cell assignments computed on the CPU (should match HcmbCellID downloaded
//                from the GPU)
//   primaryCrd:  array to hold coordinates of the primary image, with pencil cell adjustments
//                that put each atom nearest the central axis of their appropriate pencil
//                cells in Y and Z (this may put them outside the primary unit cell by up to
//                half the thickness of a pencil) with X strictly inside the primary unit cell
//---------------------------------------------------------------------------------------------
static int* CpuCalculateCellLimits(gpuContext gpu, int* cellID, double* primaryCrd)
{
  int i, j;
  int* limits;
  int3* populations;

  // Allocate the number of cells needed
  if (cellID[gpu->sim.atoms - 1] >= gpu->sim.npencils) {
    printf("CpuCalculateCellLimits :: Error.  There are %d pencils in the hash table but "
           "atoms in\nCpuCalculateCellLimits :: up to %d pencils.\n", gpu->sim.npencils,
           cellID[gpu->sim.atoms - 1]);
  }

  // Get the reciprocal space transformation matrix
  dmat U;
  U = ExtractTransformationMatrix(gpu, false);

  // Loop over all atoms and find the cell limits
  limits = (int*)malloc(8 * gpu->sim.npencils * sizeof(int));
  populations = (int3*)calloc(2 * gpu->sim.npencils, sizeof(int3));
  SetIVec(limits, 8 * gpu->sim.npencils, -1);
  int offset = 4 * gpu->sim.npencils;
  int qqcounter = 0;
  int ljcounter = 0;
  int lastqqID = -1;
  int lastljID = -1;
  for (i = 0; i < gpu->sim.atoms; i++) {
    PMEFloat2 qljid = gpu->pbHcmbChargeSPLJID->_pSysData[i];
    FloatShift ljt;
    ljt.f = qljid.y;
    ljt.ui -= 1;
    int islj = ((gpu->sim.ljexist[ljt.ui >> 5] & (0x1 << (ljt.ui & 31))) > 0);
    int isqq = (fabs(qljid.x) > 1.0e-8);

    // Compute fractional x coordinates for the atom
    double fx = U.data[0]*primaryCrd[3*i] + U.data[3]*primaryCrd[3*i + 1] +
                U.data[6]*primaryCrd[3*i + 2];

    // If the current cell ID does not match that of the last atom
    // bearing electrostatic properties, mark a new cell demarcation.
    if (isqq == 1) {

      // Mark the first and last atoms of each cell
      if (lastqqID != cellID[i]) {
        limits[cellID[i]] = qqcounter;
        if (cellID[i] > 0) {
          limits[lastqqID + (3 * gpu->sim.npencils)] = qqcounter;
        }
      }

      // Accumulate populations in each of the cell regions (low, middle, high)
      if (fx < gpu->sim.LowerCapBoundary) {
        populations[cellID[i]].x += 1;
      }
      else if (fx <= gpu->sim.UpperCapBoundary) {
        if (populations[cellID[i]].y == 0) {
          limits[cellID[i] + gpu->sim.npencils] = qqcounter;
        }
        populations[cellID[i]].y += 1;
      }
      else {
        if (populations[cellID[i]].z == 0) {
          limits[cellID[i] + (2 * gpu->sim.npencils)] = qqcounter;
        }
        populations[cellID[i]].z += 1;
      }
      qqcounter++;
      lastqqID = cellID[i];
    }
    if (islj == 1) {

      // Mark the first and last atoms of each cell
      if (lastljID != cellID[i]) {
        limits[cellID[i] + (4 * gpu->sim.npencils)] = ljcounter;
        if (cellID[i] > 0) {
          limits[lastljID + (7 * gpu->sim.npencils)] = ljcounter;
        }
      }

      // Accumulate populations in each of the cell regions (low, middle, high)
      if (fx < gpu->sim.LowerCapBoundary) {
        populations[cellID[i] + gpu->sim.npencils].x += 1;
      }
      else if (fx <= gpu->sim.UpperCapBoundary) {
        if (populations[cellID[i] + gpu->sim.npencils].y == 0) {
          limits[cellID[i] + (5 * gpu->sim.npencils)] = ljcounter;
        }
        populations[cellID[i] + gpu->sim.npencils].y += 1;
      }
      else {
        if (populations[cellID[i] + gpu->sim.npencils].z == 0) {
          limits[cellID[i] + (6 * gpu->sim.npencils)] = ljcounter;
        }
        populations[cellID[i] + gpu->sim.npencils].z += 1;
      }
      ljcounter++;
      lastljID = cellID[i];
    }

    // Special-case last atom
    if (i == gpu->sim.atoms - 1) {
      limits[lastqqID + (3 * gpu->sim.npencils)] = qqcounter;
      limits[lastljID + (7 * gpu->sim.npencils)] = ljcounter;
    }
  }

  // Check the cell populations
  int *gpulims = gpu->pbHcmbSubImgCellLimits->_pSysData;
  for (i = 0; i < gpu->sim.npencils; i++) {
    bool mismatch = false;
    for (j = 0; j < 4; j++) {
      int idx = i + (j * gpu->sim.npencils);
      if (limits[idx] != gpulims[idx]) {
        if (gpulims[idx] == -1) {
          if (j == 1 && limits[idx] == limits[i + ((j - 1) * gpu->sim.npencils)]) {
            continue;
          }
          if (j == 2 && limits[idx] == limits[i + ((j + 1) * gpu->sim.npencils)]) {
            continue;
          }
        }
        mismatch = true;
      }
    }
    if (mismatch) {
      int cellz = i / gpu->sim.nypencils;
      int celly = i - (cellz * gpu->sim.nypencils); 
      printf("CpuCalculateCellLimits :: Error.  CPU limits of cell [ %3d %3d ] do not "
             "match GPU limits\n                            [ %6d %6d %6d %6d ] (CPU) vs\n"
             "                            [ %6d %6d %6d %6d ] (GPU)\n", celly, cellz,
             limits[i], limits[i + gpu->sim.npencils], limits[i + (2 * gpu->sim.npencils)],
             limits[i + (3 * gpu->sim.npencils)], gpulims[i], gpulims[i + gpu->sim.npencils],
             gpulims[i + (2 * gpu->sim.npencils)], gpulims[i + (3 * gpu->sim.npencils)]);
      int cpurange = limits[i + (3 * gpu->sim.npencils)] - limits[i];
      int gpurange = gpulims[i + (3 * gpu->sim.npencils)] - gpulims[i];
      int jmax = (cpurange > gpurange) ? cpurange : gpurange;
      int cpubounds[4], gpubounds[4];
      for (j = 0; j < 4; j++) {
        cpubounds[j] = limits[i + (j * gpu->sim.npencils)];
        gpubounds[j] = gpulims[i + (j * gpu->sim.npencils)];
      }
      cpubounds[1] += (cpubounds[1] < 0) * (cpubounds[0] + 1);
      cpubounds[2] += (cpubounds[2] < 0) * (cpubounds[3] + 1);
      gpubounds[1] += (gpubounds[1] < 0) * (gpubounds[0] + 1);
      gpubounds[2] += (gpubounds[2] < 0) * (gpubounds[3] + 1);
      for (j = 0; j < jmax; j++) {
        int cpuidx = limits[i] + j;
        int gpuidx = gpulims[i] + j;
        if (j < cpurange) {
          printf("  CPU %9.4f %9.4f %9.4f ", primaryCrd[3*cpuidx], primaryCrd[3*cpuidx + 1],
                 primaryCrd[3*cpuidx + 2]);
          if (j < cpubounds[1] - cpubounds[0]) {
            printf("[ LCAP ] ");
          }
          else if (j < cpubounds[2] - cpubounds[0]) {
            printf("[ CORE ] ");
          }
          else {
            printf("[ HCAP ] ");
          }
        }
        if (j < gpurange) {
          printf("   GPU %9.4f %9.4f %9.4f ", primaryCrd[3*gpuidx], primaryCrd[3*gpuidx + 1],
                 primaryCrd[3*gpuidx + 2]);
          if (j < gpubounds[1] - gpubounds[0]) {
            printf("[ LCAP ] ");
          }
          else if (j < gpubounds[2] - gpubounds[0]) {
            printf("[ CORE ] ");
          }
          else {
            printf("[ HCAP ] ");
          }
        }
        printf("\n");
      }
    }
  }
  
  // Free allocated memory
  free(populations);

  return limits;
}

//---------------------------------------------------------------------------------------------
// AllocateExpandedCpuImage: allocate data to hold the expanded image for either electrostatic
//                           or Lennard-Jones interactions
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                numerous gpuBuffer attributes
//   subimg:      code for the sub-image to allocate (different cutoffs may imply different
//                levels of hash cell padding in Y and Z)
//---------------------------------------------------------------------------------------------
static eimage AllocateExpandedCpuImage(gpuContext gpu, const char* subimg)
{
  int i;

  // Codify the sub-image
  int isub;
  if (strcmp(subimg, "QQ") == 0) {
    isub = 0;
  }
  else if (strcmp(subimg, "LJ") == 0) {
    isub = 1;
  }

  int cellspace;
  eimage tIMG;
  if (isub == 0) {
    tIMG.xpadding  = gpu->sim.QQCapPadding;
    tIMG.ypadding  = gpu->sim.ypadding;
    tIMG.zpadding  = gpu->sim.zpadding;
    tIMG.cellspace = gpu->sim.QQCellSpace;
  }
  else if (isub == 1) {
    tIMG.xpadding  = gpu->sim.LJCapPadding;
    tIMG.ypadding  = gpu->sim.ypadding;
    tIMG.zpadding  = gpu->sim.zpadding;
    tIMG.cellspace = gpu->sim.LJCellSpace;
  }
  tIMG.ydimPrimary = gpu->sim.nypencils;
  tIMG.zdimPrimary = gpu->sim.nzpencils;
  tIMG.ydimFull = gpu->sim.nypencils + (2 * tIMG.ypadding);
  tIMG.zdimFull = gpu->sim.nzpencils + tIMG.zpadding;
  int ncell = tIMG.ydimFull * tIMG.zdimFull;
  int elemcount = ncell * tIMG.cellspace;
  tIMG.atoms = (int*)malloc(elemcount * sizeof(int));
  tIMG.absid = (int*)malloc(elemcount * sizeof(int));
  tIMG.crdq = (PMEFloat4*)malloc(elemcount * sizeof(PMEFloat4));
  tIMG.cells = (icontent*)malloc(ncell * sizeof(icontent));
  for (i = 0; i < ncell; i++) {
    int zcell = i / tIMG.ydimFull;
    int ycell = i - (tIMG.ydimFull * zcell);
    tIMG.cells[i].ypos = ycell;
    tIMG.cells[i].zpos = zcell;
    tIMG.cells[i].nntmap = 0;
    tIMG.cells[i].ysrc = ycell +
                         (((ycell < tIMG.ypadding) -
                           (ycell >= tIMG.ypadding + tIMG.ydimPrimary)) * tIMG.ydimPrimary);
    tIMG.cells[i].zsrc = zcell + ((zcell < tIMG.zpadding) * tIMG.zdimPrimary);
    tIMG.cells[i].atoms = &tIMG.atoms[i * tIMG.cellspace];
    tIMG.cells[i].absid = &tIMG.absid[i * tIMG.cellspace];
    tIMG.cells[i].crdq  = &tIMG.crdq[i * tIMG.cellspace];
  }

  return tIMG;
}

//---------------------------------------------------------------------------------------------
// DestroyEImage: free memory associated with an expanded image structure.
//
// Arguments:
//   tIMG:         the expanded sub-image of atoms to destroy
//---------------------------------------------------------------------------------------------
static void DestroyEImage(eimage *tIMG)
{
  int i;
  
  free(tIMG->atoms);
  free(tIMG->absid);
  free(tIMG->crdq);
  for (i = 0; i < tIMG->ydimPrimary * tIMG->zdimPrimary; i++) {
    if (tIMG->cells[i].nntmap > 0) {
      free(tIMG->cells[i].ntmaps);
      free(tIMG->cells[i].pairlimits);
      free(tIMG->cells[i].pairs);
    }
  }
  free(tIMG->cells);
}

//---------------------------------------------------------------------------------------------
// PencilNestedLoop: perform a nest loop over all atoms in each of two pencils, accounting for
//                   the need to count and then store atom neighbors in two separate cycles.
//
// Arguments:
//   tIMG:         the expanded sub-image of atoms to work from
//   hcmbPL:       the developing pair list
//   [i,j][1,2]:   indices of the cells (1 and 2) in Y (i) and Z (j)
//   nnb:          the number of non-bonded interactions (accumulated by this function)
//   xcut2:        cutoff distance (plus a pair list margin), squared
//---------------------------------------------------------------------------------------------
static int PencilNestedLoop(eimage *tIMG, pairlist *hcmbPL, int i1, int j1, int i2, int j2,
                            int nnb, double xcut2)
{
  int i, j;
  
  // Find the first pencil limits
  int Lcore  = tIMG->xpadding;
  int cidx1  = (tIMG->ydimFull * j1) + i1;
  int Hcore1 = tIMG->cells[cidx1].coreEnd;

  // Find the scond pencil limits
  int cidx2  = (tIMG->ydimFull * j2) + i2;  
  int Lcap2  = tIMG->cells[cidx2].lcapStart;
  int Hcap2  = tIMG->cells[cidx2].hcapEnd;

  // Nested loop over all atoms in the primary image of the first cell
  // and all atoms generally in the second cell
  int *atomid1, *atomid2;
  PMEFloat4 *crdq1, *crdq2;
  crdq1 = tIMG->cells[cidx1].crdq;
  crdq2 = tIMG->cells[cidx2].crdq;
  atomid1 = tIMG->cells[cidx1].atoms;
  atomid2 = tIMG->cells[cidx2].atoms;
  if (hcmbPL->ngbrRanges == NULL) {
    for (i = Lcore; i < Hcore1; i++) {
      PMEFloat4 packet;
      packet = crdq1[i];
      double atmx = packet.x;
      double atmy = packet.y;
      double atmz = packet.z;
      int atmid1 = atomid1[i];      
      for (j = Lcap2; j < Hcap2; j++) {
        if (cidx1 == cidx2 && j >= i) {
          continue;
        }
        packet = crdq2[j];
        double dx = packet.x - atmx;
        double dy = packet.y - atmy;
        double dz = packet.z - atmz;
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < xcut2) {
          int atmid2 = atomid2[j];
          if (atmid1 < atmid2) {
            hcmbPL->ngbrCounts[atmid1] += 1;
          }
          else {
            hcmbPL->ngbrCounts[atmid2] += 1;
          }
          nnb++;
        }
      }
    }
  }
  else {
    for (i = Lcore; i < Hcore1; i++) {
      PMEFloat4 packet;
      packet = crdq1[i];
      double atmx = packet.x;
      double atmy = packet.y;
      double atmz = packet.z;
      int atmid1 = atomid1[i];
      for (j = Lcap2; j < Hcap2; j++) {
        if (cidx1 == cidx2 && j >= i) {
          continue;
        }
        packet = crdq2[j];
        double dx = packet.x - atmx;
        double dy = packet.y - atmy;
        double dz = packet.z - atmz;
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < xcut2) {
          int atmid2 = atomid2[j];
          if (atmid1 < atmid2) {
            int plidx = hcmbPL->ngbrCounts[atmid1];
            hcmbPL->ngbrList[plidx] = atmid2;
            hcmbPL->ngbrRanges[plidx] = sqrt(r2);
            hcmbPL->ngbrCounts[atmid1] = plidx + 1;
          }
          else {
            int plidx = hcmbPL->ngbrCounts[atmid2];
            hcmbPL->ngbrList[plidx] = atmid1;
            hcmbPL->ngbrRanges[plidx] = sqrt(r2);
            hcmbPL->ngbrCounts[atmid2] = plidx + 1;
          }
        }
      }
    }
  }

  return nnb;
}

//---------------------------------------------------------------------------------------------
// FinalizePairlistAllocation: based on an array of all atoms' neighbor counts, allocate the
//                             necessary memory for storing the rest of the data.
//
// This is a debugging function.
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                object counts
//   thisPL:      the developing pair list   
//---------------------------------------------------------------------------------------------
static void FinalizePairlistAllocation(gpuContext gpu, pairlist *thisPL)
{
  int i, tnb;

  tnb = 0;
  for (i = 0; i < gpu->sim.atoms; i++) {
    int tmp = thisPL->ngbrCounts[i];
    thisPL->ngbrCounts[i] = tnb;
    tnb += tmp;
  }
  thisPL->natom = gpu->sim.atoms;
  thisPL->ngbrCounts[gpu->sim.atoms] = tnb;
  thisPL->ngbrList = (int*)malloc(tnb * sizeof(int));
  thisPL->ngbrRanges = (PMEFloat*)malloc(tnb * sizeof(PMEFloat));
}

//---------------------------------------------------------------------------------------------
// DestroyPairlist: free memory associated with a CPU-constructed, mock-up neighbor list.
//
// This is a debugging function.
//
// Arguments:
//   thisPL:     the pair list to free
//---------------------------------------------------------------------------------------------
static void DestroyPairlist(pairlist *thisPL)
{
  free(thisPL->ngbrCounts);
  free(thisPL->ngbrList);
  free(thisPL->ngbrRanges);
}

//---------------------------------------------------------------------------------------------
// SearchGridForAtom: search the expanded image for a particular atom and report all places
//                    that atom is found.
//
// This is a debugging function.
//
// Arguments:
//   tIMG:        the expanded sub-image representation
//   atmID:       the ID number of the atom to seek out
//   searchAll:   seek out the atom in every nook and cranny of the expanded image
//---------------------------------------------------------------------------------------------
static void SearchGridForAtom(eimage *tIMG, int srcID, int atmID, bool searchAll)
{
  int i, j, k;

  // Find the source atom in the primary image
  int yloc, zloc, kloc = -1;
  double atmx, atmy, atmz;
  for (i = tIMG->ypadding; i < tIMG->ydimPrimary + tIMG->ypadding; i++) {
    for (j = tIMG->zpadding; j < tIMG->zdimFull; j++) {
      int cidx = (j * tIMG->ydimFull) + i;
      for (k = tIMG->cells[cidx].coreStart; k < tIMG->cells[cidx].coreEnd; k++) {
        if (tIMG->cells[cidx].atoms[k] == srcID) {
          yloc = i;
          zloc = j;
          kloc = k;
          PMEFloat4 packet;
          packet = tIMG->cells[cidx].crdq[k];
          atmx = packet.x;
          atmy = packet.y;
          atmz = packet.z;
        }
      }
    }
  }
  if (kloc == -1) {
    printf("SearchGridForAtom :: failed to find atom %6d in the primary unit cell.\n", atmID);
    return;
  }
  if (searchAll) {

    // Find the atom in the primary
    bool found = false;
    for (i = 0; i < tIMG->ydimFull * tIMG->zdimFull; i++) {
      for (j = tIMG->cells[i].lcapStart; j < tIMG->cells[i].hcapEnd; j++) {
        if (tIMG->cells[i].atoms[j] == atmID) {
          if (found) {
            printf("      [                   ");
          }
          int zpos = i / tIMG->ydimFull;
          int ypos = i - (zpos * tIMG->ydimFull);
          printf("Found at %2d %2d - %3d ", ypos, zpos, j);
          if (j < tIMG->cells[i].coreStart) {        
            printf("LCAP ");
          }
          else if (j < tIMG->cells[i].coreEnd) {
            printf("CORE ");
          }
          else {
            printf("HCAP ");
          }
          PMEFloat4 packet = tIMG->cells[i].crdq[j];
          double dx = packet.x - atmx;
          double dy = packet.y - atmy;
          double dz = packet.z - atmz;
          double rdist = sqrt(dx*dx + dy*dy + dz*dz);
          printf("in %3d %3d %3d %3d : %8.4f ]\n", tIMG->cells[i].lcapStart,
                 tIMG->cells[i].coreStart, tIMG->cells[i].coreEnd, tIMG->cells[i].hcapEnd,
                 rdist);
          found = true;
        }
      }
    }
  }
  else {
    printf("      [ Atom %6d primary at     %2d %2d - %3d ", srcID, yloc, zloc, kloc);
    int cidx = (zloc * tIMG->ydimFull) + yloc;
    if (kloc < tIMG->cells[cidx].coreStart) {
      printf("LCAP ");
    }
    else if (kloc < tIMG->cells[cidx].coreEnd) {
      printf("CORE ");
    }
    else {
      printf("HCAP ");
    }
    printf("in %3d %3d %3d %3d ]\n", tIMG->cells[cidx].lcapStart, tIMG->cells[cidx].coreStart,
           tIMG->cells[cidx].coreEnd, tIMG->cells[cidx].hcapEnd);
  }
}

//---------------------------------------------------------------------------------------------
// MatchNeighbors: function for exhaustively checking the neighbor lists for all pairs and
//                 distances.  This is a rigorous test of the CPU neighbor list generator.
//
// This is a debugging function.
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                object counts
//   tIMG:        the expanded sub-image representation
//   thisPL:      the completed neighbor list derived by some complex means
//   simplePL:    the completed neighbor list derived by the simplest possible means (assumed
//                to be correct)
//   desc:        descriptor for the more complex neighbor list
//---------------------------------------------------------------------------------------------
static bool MatchNeighbors(gpuContext gpu, eimage *tIMG, pairlist *thisPL, pairlist *simplePL,
                           const char* desc)
{
  int i, j, k, m, n;
  int* missing;
  int* badrange;
  bool* duplicates;
  bool success;
  
  // Start out assuming success
  success = true;
  j = 0;
  for (i = 0; i < gpu->sim.atoms; i++) {
    if (simplePL->ngbrCounts[i+1] - simplePL->ngbrCounts[i] > j) {
      j = simplePL->ngbrCounts[i+1] - simplePL->ngbrCounts[i];
    }
    if (thisPL->ngbrCounts[i+1] - thisPL->ngbrCounts[i] > j) {
      j = thisPL->ngbrCounts[i+1] - thisPL->ngbrCounts[i];
    }
  }
  missing = (int*)malloc(j * sizeof(int));
  badrange = (int*)malloc(j * sizeof(int));
  duplicates = (bool*)calloc(j, sizeof(bool));
  for (i = 0; i < gpu->sim.atoms; i++) {

    // No need to cull based on the property of interest being present here--that will
    // already have been taken care of implicitly by the neighbor list generation
    int nmissing = 0;
    int nbadrange = 0;
    int nduplicate = 0;
    for (j = simplePL->ngbrCounts[i]; j < simplePL->ngbrCounts[i+1]; j++) {
      int atmid = simplePL->ngbrList[j];
      PMEFloat range = simplePL->ngbrRanges[j];
      bool found = false;
      bool wrongr2 = false;
      for (k = thisPL->ngbrCounts[i]; k < thisPL->ngbrCounts[i+1]; k++) {
        if (thisPL->ngbrList[k] == atmid) {
          if (found == false) {

            // The neighbor atom is now found, but is the distance right?
            found = true;
            if (fabs(thisPL->ngbrRanges[k] - range) > 1.0e-4) {
              wrongr2 = true;
            }
          }
          else {

            // The neighbor atom has already been found.  Mark the mistake.
            int2 dup;
            dup.x = j;
            dup.y = k;
            duplicates[j - simplePL->ngbrCounts[i]] = true;
            nduplicate++;
          }
        }
      }
      if (found == false) {
        missing[nmissing] = atmid;
        nmissing++;
      }
      if (found && wrongr2) {
        badrange[nbadrange] = atmid;
        nbadrange++;
      }
    }
    if (nmissing > 0) {
      printf("MatchNeighbors :: In neighbor list %s atom %6d, %d neighbors are missing\n",
             desc, i, nmissing);
      SearchGridForAtom(tIMG, i, 0, false);
      for (j = 0; j < nmissing; j++) {
        printf("      [ %6d ", missing[j]);
        for (k = simplePL->ngbrCounts[i]; k < simplePL->ngbrCounts[i+1]; k++) {
          if (simplePL->ngbrList[k] == missing[j]) {
            printf("(%8.4f) ", simplePL->ngbrRanges[k]);
            SearchGridForAtom(tIMG, i, missing[j], true);
          }
        }
      }
      printf("\n");
      success = false;
    }
    if (nbadrange > 0) {
      printf("MatchNeighbors :: In neighbor list %s atom %6d, %d neighbors are found with\n"
             "                  the wrong ranges\n                  [ ", desc, i, nbadrange);
      m = 0;
      for (j = 0; j < nbadrange; j++) {
        printf("%6d ", badrange[j]);
        for (k = simplePL->ngbrCounts[i]; k < simplePL->ngbrCounts[i+1]; k++) {
          if (simplePL->ngbrList[k] == badrange[j]) {
            printf("(%8.4f : ", simplePL->ngbrRanges[k]);
          }
        }
        for (k = thisPL->ngbrCounts[i]; k < thisPL->ngbrCounts[i+1]; k++) {
          if (thisPL->ngbrList[k] == badrange[j]) {
            printf("%8.4f) ", thisPL->ngbrRanges[k]);
          }
        }
        m++;
        if (m == 4) {
          printf("\n           ");
          m = 0;
        }
      }
      printf("]\n");
      success = false;
    }
    if (nduplicate > 0) {
      printf("MatchNeighbors :: In neighbor list %s atom %6d, %d neighbor atoms have "
             "duplicates:\n", desc, i, nduplicate);
      m = 0;
      for (j = simplePL->ngbrCounts[i]; j < simplePL->ngbrCounts[i+1]; j++) {
        if (duplicates[j - simplePL->ngbrCounts[i]]) {
          printf("                  [ Atom %6d (%8.4f) : ", simplePL->ngbrList[j],
                 simplePL->ngbrRanges[j]);
          for (k = thisPL->ngbrCounts[i]; k < thisPL->ngbrCounts[i+1]; k++) {
            if (thisPL->ngbrList[k] == simplePL->ngbrList[j]) {
              printf("%8.4f ", thisPL->ngbrRanges[k]);
            }
          }
          printf("]\n");
        }
      }
      success = false;
    }
  }

  // Free allocated memory
  free(missing);
  free(badrange);
  free(duplicates);

  return success;
}

//---------------------------------------------------------------------------------------------
// DefineActiveHashCellZones: define zones of each hash cell's memory in the expanded image,
//                            including the core and low or high end caps, which shall contain
//                            relevant atoms.  This can be done irrespective of whether atoms
//                            have been loaded into the hash cells or not, and can serve as a
//                            sanity check when the expanded image does get loaded.
//
// This is a debugging function.
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                object counts
//   tIMG:        the expanded sub-image representation
//   limits:      the limits of each cell in the primary image (lower index, lower cap index,
//                upper cap index, upper index)
//   subimg:      code for the sub-image to allocate (different cutoffs may imply different
//                levels of hash cell padding in Y and Z)
//---------------------------------------------------------------------------------------------
static void DefineActiveHashCellZones(gpuContext gpu, eimage *tIMG, int* limits,
                                      const char* subimg)
{
  int i, j;
  int *limptr;
  
  // Now, expand the primary image using a serial, CPU-based approach
  if (strcmp(subimg, "QQ") == 0) {
    limptr = limits;
  }
  else if (strcmp(subimg, "LJ") == 0) {
    limptr = &limits[4 * tIMG->ydimPrimary * tIMG->zdimPrimary];
  }
  int npencils = tIMG->ydimPrimary * tIMG->zdimPrimary;
  for (i = 0; i < tIMG->ydimFull; i++) {
    for (j = 0; j < tIMG->zdimFull; j++) {

      // Cell indexing (this cell expanded image cidx, primary image pcidx, source scidx)
      int cidx = (tIMG->ydimFull * j) + i;
      int pcidx = ((tIMG->cells[cidx].zsrc - tIMG->zpadding) * tIMG->ydimPrimary) +
                  tIMG->cells[cidx].ysrc - tIMG->ypadding;

      // Get the cell limits.  Record critical points in the atom list.
      int lim1 = limptr[pcidx];
      int lim2 = limptr[  npencils + pcidx];
      int lim3 = limptr[2*npencils + pcidx];
      int lim4 = limptr[3*npencils + pcidx];
      lim2 += (lim2 < 0)*(lim1 + 1);
      lim3 += (lim3 < 0)*(lim4 + 1);
      tIMG->cells[cidx].coreStart = tIMG->xpadding;
      tIMG->cells[cidx].coreEnd = tIMG->xpadding + lim4 - lim1;
      tIMG->cells[cidx].lcapStart = tIMG->xpadding - (lim4 - lim3);
      tIMG->cells[cidx].hcapEnd = tIMG->xpadding + lim4 + lim2 - 2*lim1;
    }
  }
}

//---------------------------------------------------------------------------------------------
// LoadPrimaryImage: load the primary image of the expanded cell grid using the first set of
//                   GPU-written instructions.
//
// This is a debugging function.
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                object counts
//   tIMG:        the expanded sub-image representation
//   limits:      the limits of each cell in the primary image (lower index, lower cap index,
//                upper cap index, upper index)
//   primaryCrd:  the coordinates of all atoms in the primary image
//   subimg:      code for the sub-image to allocate (different cutoffs may imply different
//                levels of hash cell padding in Y and Z)
//---------------------------------------------------------------------------------------------
static void LoadPrimaryImage(gpuContext gpu, eimage *tIMG, int* limits, double* primaryCrd,
                             const char* subimg)
{
  int i, isub;
  int *limptr;

  // Set pointers and cutoffs
  int npencils = tIMG->ydimPrimary * tIMG->zdimPrimary;
  unsigned int *HcmbAtomLookup = &(gpu->pbHcmbIndex->_pSysData[gpu->sim.imageStride * 2]);
  if (strcmp(subimg, "QQ") == 0) {
    limptr = limits;
    isub = 0;
  }
  else if (strcmp(subimg, "LJ") == 0) {
    limptr = &limits[4 * tIMG->ydimPrimary * tIMG->zdimPrimary];
    isub = 1;
  }
  for (i = 0; i < gpu->sim.atoms; i++) {
    int2 filter = gpu->pbSubImageFilter->_pSysData[i];
    PMEFloat4 packet;
    packet.x = primaryCrd[3*i    ];
    packet.y = primaryCrd[3*i + 1];
    packet.z = primaryCrd[3*i + 2];
    int cidx, cpos;
    if (isub == 0 && filter.x >= 0) {
      cidx = filter.x / tIMG->cellspace;
      cpos = filter.x - (cidx * tIMG->cellspace);
      packet.w = gpu->pbHcmbChargeSPLJID->_pSysData[i].x;
      tIMG->atoms[filter.x] = i;
      tIMG->absid[filter.x] = gpu->pbHcmbIndex->_pSysData[i];
      tIMG->crdq[filter.x] = packet;
    }
    else if (isub == 1 && filter.y >= 0) {
      cidx = filter.y / tIMG->cellspace;
      cpos = filter.y - (cidx * tIMG->cellspace);
      packet.w = gpu->pbHcmbChargeSPLJID->_pSysData[i].y;
      tIMG->atoms[filter.y] = i;
      tIMG->absid[filter.y] = gpu->pbHcmbIndex->_pSysData[i];
      tIMG->crdq[filter.y] = packet;
    }
    else {
      continue;
    }

    // Check the cell positioning of the filter indices
    int zpos = cidx / tIMG->ydimFull;
    int ypos = cidx - (zpos * tIMG->ydimFull);
    if (ypos < tIMG->ypadding || ypos >= tIMG->ypadding + tIMG->ydimPrimary ||
        zpos < tIMG->zpadding || zpos >= tIMG->zpadding + tIMG->zdimPrimary) {
      printf("LoadPrimaryImage :: (%s) Atom %6d of the re-ordered image is placed outside of "
             "the\n                    allowed cell boundaries [ placed in %2d %2d, needed "
             "%2d-%2d %2d-%2d ]\n", subimg, i, ypos, zpos, tIMG->ypadding,
             tIMG->ypadding + tIMG->ydimPrimary, tIMG->zpadding,
             tIMG->zpadding + tIMG->zdimPrimary);
    }
    if (cpos < tIMG->cells[cidx].coreStart || cpos >= tIMG->cells[cidx].coreEnd) {
      printf("LoadPrimaryImage :: (%s) Atom %6d of the re-ordered image is placed outside of "
             "the\n                    allowed cell boundaries [ placed at %3d, needed "
             "%3d-%3d ]\n", subimg, i, cpos, tIMG->cells[cidx].coreStart,
             tIMG->cells[cidx].coreEnd);
    }
  }
}

//---------------------------------------------------------------------------------------------
// TakeImageChecksum: take the sum of coordinates in a part or whole of the expanded image
//                    representation.  This routine is agnostic to how the image was made.
//
// This is a debugging function.
//
// Arguments:
//   tIMG:        the expanded sub-image representation
//   limits:      the limits of each cell in the primary image (lower index, lower cap index,
//                upper cap index, upper index)
//---------------------------------------------------------------------------------------------
static double4 TakeImageChecksum(gpuContext gpu, eimage *tIMG, const char* subimg,
                                 const char* region)
{
  int i, j, k, isub;
  int *limptr;
  
  // Set pointers and cutoffs
  int npencils = tIMG->ydimPrimary * tIMG->zdimPrimary;
  if (strcmp(subimg, "QQ") == 0) {
    isub = 0;
  }
  else if (strcmp(subimg, "LJ") == 0) {
    isub = 1;
  }
  double4 chksum;
  chksum.x = 0.0;
  chksum.y = 0.0;
  chksum.z = 0.0;
  chksum.w = 0.0;
  if (strcmp(region, "PRIMARY") == 0) {
    for (i = tIMG->ypadding; i < tIMG->ydimFull - tIMG->ypadding; i++) {
      for (j = tIMG->zpadding; j < tIMG->zdimFull; j++) {
        int cidx = (j * tIMG->ydimFull) + i;
        for (k = tIMG->cells[cidx].coreStart; k < tIMG->cells[cidx].coreEnd; k++) {
          PMEFloat4 packet = tIMG->cells[cidx].crdq[k];
          chksum.x += packet.x;
          chksum.y += packet.y;
          chksum.z += packet.z;
          if (isub == 0) {
            chksum.w += packet.w;
          }
          else {
            FloatShift xtyp;
            xtyp.f = packet.w;
            chksum.w += (PMEFloat)(xtyp.ui);
          }
        }
      }
    }
  }
  else if (strcmp(region, "ALL") == 0) {
    for (i = 0; i < tIMG->ydimFull; i++) {
      for (j = 0; j < tIMG->zdimFull; j++) {
        int cidx = (j * tIMG->ydimFull) + i;
        for (k = tIMG->cells[cidx].lcapStart; k < tIMG->cells[cidx].hcapEnd; k++) {
          PMEFloat4 packet = tIMG->cells[cidx].crdq[k];
          chksum.x += packet.x;
          chksum.y += packet.y;
          chksum.z += packet.z;
          chksum.w += packet.w;
          if (isub == 0) {
            chksum.w += packet.w;
          }
          else {
            FloatShift xtyp;
            xtyp.f = packet.w;
            chksum.w += (PMEFloat)(xtyp.ui);
          }
        }
      }
    }
  }

  return chksum;
}

//---------------------------------------------------------------------------------------------
// ExpandCpuImage: put all atoms into the primary image, in preparation to expand that to
//                   completely fill the pre-imaged space.
//
// This is a debugging function.
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                object counts
//   tIMG:        the expanded sub-image representation
//   limits:      the limits of each cell in the primary image (lower index, lower cap index,
//                upper cap index, upper index)
//   primaryCrd:  the coordinates of all atoms in the primary image
//   subimg:      code for the sub-image to allocate (different cutoffs may imply different
//                levels of hash cell padding in Y and Z)
//---------------------------------------------------------------------------------------------
static void ExpandCpuImage(gpuContext gpu, eimage *tIMG, int* limits, double* primaryCrd,
                             const char* subimg)
{
  int i, j, k, isub;
  int *limptr;
  unsigned int *HcmbIDs;
  
  // Make local copies of the transformation matrices and cutoff
  dmat U, invU;
  U = ExtractTransformationMatrix(gpu, false);
  invU = ExtractTransformationMatrix(gpu, true);
  int npencils = tIMG->ydimPrimary * tIMG->zdimPrimary;
  
  // Initialize the coordinate and property information
  j = tIMG->ydimFull * tIMG->zdimFull * tIMG->cellspace;
  for (i = 0; i < j; i++) {
    PMEFloat4 packet;
    packet.x = 1.0e8;
    packet.y = 1.0e8;
    packet.z = 1.0e8;
    FloatShift xtyp;
    xtyp.ui = 8192;
    packet.w = xtyp.f;
    tIMG->crdq[i] = packet;
  }
  
  // Set pointers and cutoffs
  double xcut2;
  if (strcmp(subimg, "QQ") == 0) {
    xcut2 = gpu->sim.es_cutoff + gpu->sim.skinnb;
    limptr = limits;
    isub = 0;
  }
  else if (strcmp(subimg, "LJ") == 0) {
    xcut2 = gpu->sim.lj_cutoff + gpu->sim.skinnb;
    limptr = &limits[4 * tIMG->ydimPrimary * tIMG->zdimPrimary];
    isub = 1;
  }
  xcut2 *= xcut2;
  if (gpu->sim.pHcmbIndex == gpu->pbHcmbIndex->_pDevData) {
    HcmbIDs = gpu->pbHcmbIndex->_pSysData;
  }
  else {
    HcmbIDs = gpu->pbHcmbIndex->_pSysData + (gpu->sim.imageStride * 4);
  }
  LoadPrimaryImage(gpu, tIMG, limits, primaryCrd, subimg);
  double4 chksum = TakeImageChecksum(gpu, tIMG, subimg, "PRIMARY");
  
  // Now, expand the primary image using a serial, CPU-based approach
  for (i = 0; i < tIMG->ydimFull; i++) {
    for (j = 0; j < tIMG->zdimFull; j++) {

      // Cell indexing (this cell expanded image cidx, primary image pcidx, source scidx)
      int cidx = (tIMG->ydimFull * j) + i;
      int pcidx = ((tIMG->cells[cidx].zsrc - tIMG->zpadding) * tIMG->ydimPrimary) +
                  tIMG->cells[cidx].ysrc - tIMG->ypadding;
      int scidx = (tIMG->cells[cidx].zsrc * tIMG->ydimFull) + tIMG->cells[cidx].ysrc;

      // Moves in Y and Z that each atom makes
      double ymove = (tIMG->cells[cidx].ypos > tIMG->cells[cidx].ysrc) ? 1.0 : 0.0;
      ymove = (tIMG->cells[cidx].ypos < tIMG->ypadding) ? -1.0 : ymove;
      double zmove = (tIMG->cells[cidx].zpos > tIMG->cells[cidx].zsrc) ? 1.0 : 0.0;
      zmove = (tIMG->cells[cidx].zpos < tIMG->zpadding) ? -1.0 : zmove;

      // Get the cell limits.  Record critical points in the atom list.
      int lim1 = limptr[pcidx];
      int lim2 = limptr[  npencils + pcidx];
      int lim3 = limptr[2*npencils + pcidx];
      int lim4 = limptr[3*npencils + pcidx];      
      lim2 += (lim2 < 0)*(lim1 + 1);
      lim3 += (lim3 < 0)*(lim4 + 1);

      // Expand the end caps
      for (k = lim1; k < lim2; k++) {
        int orig = tIMG->xpadding + k - lim1;
        int dest = orig + lim4 - lim1;
        tIMG->cells[cidx].atoms[dest] = tIMG->cells[scidx].atoms[orig];
        tIMG->cells[cidx].absid[dest] = tIMG->cells[scidx].absid[orig];
        PMEFloat4 packet;
        packet = tIMG->cells[scidx].crdq[orig];
        double xmove = 1.0;
        double newx = U.data[0]*packet.x + U.data[3]*packet.y + U.data[6]*packet.z + xmove;
        double newy = U.data[1]*packet.x + U.data[4]*packet.y + U.data[7]*packet.z + ymove;
        double newz = U.data[2]*packet.x + U.data[5]*packet.y + U.data[8]*packet.z + zmove;
        packet.x = invU.data[0]*newx + invU.data[1]*newy + invU.data[2]*newz;
        packet.y = invU.data[3]*newx + invU.data[4]*newy + invU.data[5]*newz;
        packet.z = invU.data[6]*newx + invU.data[7]*newy + invU.data[8]*newz;
        tIMG->cells[cidx].crdq[dest] = packet;
      }
      for (k = lim3; k < lim4; k++) {
        int orig = tIMG->xpadding + k - lim1;
        int dest = tIMG->xpadding + k - lim4;
        tIMG->cells[cidx].atoms[dest] = tIMG->cells[scidx].atoms[orig];
        tIMG->cells[cidx].absid[dest] = tIMG->cells[scidx].absid[orig];
        PMEFloat4 packet;
        packet = tIMG->cells[scidx].crdq[orig];
        double xmove = -1.0;
        double newx = U.data[0]*packet.x + U.data[3]*packet.y + U.data[6]*packet.z + xmove;
        double newy = U.data[1]*packet.x + U.data[4]*packet.y + U.data[7]*packet.z + ymove;
        double newz = U.data[2]*packet.x + U.data[5]*packet.y + U.data[8]*packet.z + zmove;
        packet.x = invU.data[0]*newx + invU.data[1]*newy + invU.data[2]*newz;
        packet.y = invU.data[3]*newx + invU.data[4]*newy + invU.data[5]*newz;
        packet.z = invU.data[6]*newx + invU.data[7]*newy + invU.data[8]*newz;
        tIMG->cells[cidx].crdq[dest] = packet;
      }

      // Expand cells outside the primary image
      if (scidx != cidx) {
        for (k = lim1; k < lim4; k++) {
          int pos = tIMG->xpadding + k - lim1;
          tIMG->cells[cidx].atoms[pos] = tIMG->cells[scidx].atoms[pos];
          tIMG->cells[cidx].absid[pos] = tIMG->cells[scidx].absid[pos];
          PMEFloat4 packet;
          packet = tIMG->cells[scidx].crdq[pos];
          double xmove = 0.0;
          double newx = U.data[0]*packet.x + U.data[3]*packet.y + U.data[6]*packet.z + xmove;
          double newy = U.data[1]*packet.x + U.data[4]*packet.y + U.data[7]*packet.z + ymove;
          double newz = U.data[2]*packet.x + U.data[5]*packet.y + U.data[8]*packet.z + zmove;
          packet.x = invU.data[0]*newx + invU.data[1]*newy + invU.data[2]*newz;
          packet.y = invU.data[3]*newx + invU.data[4]*newy + invU.data[5]*newz;
          packet.z = invU.data[6]*newx + invU.data[7]*newy + invU.data[8]*newz;
          tIMG->cells[cidx].crdq[pos] = packet;
        }
      }
    }
  }

  // Check over the primary image in the expanded representation: does it look all right?
  double4 primsum = TakeImageChecksum(gpu, tIMG, subimg, "PRIMARY");
  if (fabs(chksum.x - primsum.x) > 1.0e-4 || fabs(chksum.y - primsum.y) > 1.0e-4 ||
      fabs(chksum.z - primsum.z) > 1.0e-4 || fabs(chksum.w - primsum.w) > 1.0e-4) {
    printf("ExpandCpuImage :: (%s) Primary image mapping failure:\n", subimg);
    printf("ExpandCpuImage :: (%s) Check sum   X Y Z P : [ %12.4f %12.4f %12.4f %12.4f ]\n",
           subimg, chksum.x, chksum.y, chksum.z, chksum.w);
    printf("ExpandCpuImage :: (%s) Primary sum X Y Z P : [ %12.4f %12.4f %12.4f %12.4f ]\n",
           subimg, primsum.x, primsum.y, primsum.z, primsum.w);
  }
  else {
    printf("ExpandCpuImage :: (%s) Primary image coordinate checksum passed.\n", subimg);
  }
  for (i = 0; i < tIMG->ydimFull; i++) {
    for (j = 0; j < tIMG->zdimFull; j++) {
      int cidx = (tIMG->ydimFull * j) + i;
      for (k = tIMG->cells[cidx].lcapStart; k < tIMG->cells[cidx].hcapEnd; k++) {
        PMEFloat4 packet = tIMG->cells[cidx].crdq[k];
        FloatShift xtyp;
        xtyp.f = packet.w;
        if (packet.x > 0.9e8 || packet.y > 0.9e8 || packet.z > 0.9e8 || xtyp.ui == 8192) {
          printf("ExpandCpuImage :: (%s) Cell [ %2d %2d ] atom [ %3d ] is unwritten.\n",
                 subimg, i, j, k);
        }
      }
    }
  }

  // Free allocated memory
  DestroyDmat(&U);
  DestroyDmat(&invU);
}

//---------------------------------------------------------------------------------------------
// BuildSimplePairlist: build a pairlist based on the primary coordinates using an all-to-all
//                      comparison with explicit re-imaging.  This is terribly costly method
//                      but the simplest available.
//
// This is a debugging function.
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                object counts
//   primaryCrd:  the coordinates of all atoms in the primary image
//   subimg:      code for the sub-image to allocate (different cutoffs may imply different
//                levels of hash cell padding in Y and Z)
//---------------------------------------------------------------------------------------------
static pairlist BuildSimplePairlist(gpuContext gpu, double* primaryCrd, const char* subimg)
{
  int i, j, isub;
  
  // Determine the sub-image of interest
  double xcut2;
  if (strcmp(subimg, "QQ") == 0) {
    isub = 0;
    xcut2 = gpu->sim.es_cutoff + gpu->sim.skinnb;
  }
  else if (strcmp(subimg, "LJ") == 0) {
    isub = 1;
    xcut2 = gpu->sim.lj_cutoff + gpu->sim.skinnb;
  }
  xcut2 *= xcut2;
  
  // Make local copies of the transformation matrices and cutoff
  dmat U, invU;
  U = ExtractTransformationMatrix(gpu, false);
  invU = ExtractTransformationMatrix(gpu, true);

  // Loop over all atoms and make an explicit pair list with distance measurements.
  // This first construction performs a nested loop over all atoms in the system.
  pairlist simplePL;
  simplePL.ngbrCounts = (int*)calloc(gpu->sim.atoms + 1, sizeof(int));
  for (i = 0; i < gpu->sim.atoms - 1; i++) {
    int2 filter = gpu->pbSubImageFilter->_pSysData[i];
    if ((isub == 0 && filter.x == -1) || (isub == 1 && filter.y == -1)) {
      continue;
    }
    double atmx = primaryCrd[3*i    ];
    double atmy = primaryCrd[3*i + 1];
    double atmz = primaryCrd[3*i + 2];
    int nnb = 0;
    for (j = i + 1; j < gpu->sim.atoms; j++) {
      filter = gpu->pbSubImageFilter->_pSysData[j];
      if ((isub == 0 && filter.x == -1) || (isub == 1 && filter.y == -1)) {
        continue;
      }
      double dx = primaryCrd[3*j    ] - atmx;
      double dy = primaryCrd[3*j + 1] - atmy;
      double dz = primaryCrd[3*j + 2] - atmz;
      double ndx = U.data[0]*dx + U.data[3]*dy + U.data[6]*dz;
      double ndy = U.data[1]*dx + U.data[4]*dy + U.data[7]*dz;
      double ndz = U.data[2]*dx + U.data[5]*dy + U.data[8]*dz;
      ndx += (double)((ndx < -0.5) - (ndx >= 0.5));
      ndy += (double)((ndy < -0.5) - (ndy >= 0.5));
      ndz += (double)((ndz < -0.5) - (ndz >= 0.5));
      dx = invU.data[0]*ndx + invU.data[1]*ndy + invU.data[2]*ndz;
      dy = invU.data[3]*ndx + invU.data[4]*ndy + invU.data[5]*ndz;
      dz = invU.data[6]*ndx + invU.data[7]*ndy + invU.data[8]*ndz;
      double r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < xcut2) {
        nnb++;
      }
    }
    simplePL.ngbrCounts[i] = nnb;
  }  
  FinalizePairlistAllocation(gpu, &simplePL);
  int nnb = 0;
  for (i = 0; i < gpu->sim.atoms - 1; i++) {
    int2 filter = gpu->pbSubImageFilter->_pSysData[i];
    if ((isub == 0 && filter.x == -1) || (isub == 1 && filter.y == -1)) {
      continue;
    }
    double atmx = primaryCrd[3*i    ];
    double atmy = primaryCrd[3*i + 1];
    double atmz = primaryCrd[3*i + 2];
    for (j = i + 1; j < gpu->sim.atoms; j++) {
      filter = gpu->pbSubImageFilter->_pSysData[j];
      if ((isub == 0 && filter.x == -1) || (isub == 1 && filter.y == -1)) {
        continue;
      }
      double dx = primaryCrd[3*j    ] - atmx;
      double dy = primaryCrd[3*j + 1] - atmy;
      double dz = primaryCrd[3*j + 2] - atmz;
      double ndx = U.data[0]*dx + U.data[3]*dy + U.data[6]*dz;
      double ndy = U.data[1]*dx + U.data[4]*dy + U.data[7]*dz;
      double ndz = U.data[2]*dx + U.data[5]*dy + U.data[8]*dz;
      ndx += (double)((ndx < -0.5) - (ndx >= 0.5));
      ndy += (double)((ndy < -0.5) - (ndy >= 0.5));
      ndz += (double)((ndz < -0.5) - (ndz >= 0.5));
      dx = invU.data[0]*ndx + invU.data[1]*ndy + invU.data[2]*ndz;
      dy = invU.data[3]*ndx + invU.data[4]*ndy + invU.data[5]*ndz;
      dz = invU.data[6]*ndx + invU.data[7]*ndy + invU.data[8]*ndz;
      double r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < xcut2) {
        simplePL.ngbrList[nnb] = j;
        simplePL.ngbrRanges[nnb] = sqrt(r2);
        nnb++;
      }
    }
  }
  printf("BuildSimplePairlist :: There are %d interactions.\n", nnb);

  // Free allocated memory
  DestroyDmat(&U);
  DestroyDmat(&invU);

  return simplePL;
}

//---------------------------------------------------------------------------------------------
// BuildCorePairlist: build a neighbor list usign data in the core regions of primary image
//                    hash cells.  This is another costly method, as all primary image hash
//                    cells' core regions will be covered in the nested loop.  The number of
//                    atomic interactions tested will be the same as in BuildSimplePairlist()
//                    above, and again all tests will involve explicit re-imaging.
//
// This is a debugging function.
//
// Arguments:                    
//   gpu:         overarching data structure containing simulation information, here used for
//                object counts
//   tIMG:        the expanded sub-image representation
//   simplePL:    a reference pair list computed by the simplest possible method
//   subimg:      code for the sub-image to allocate (different cutoffs may imply different
//                levels of hash cell padding in Y and Z)
//   checkPairs:  flag to have all apir interactions checked against the reference
//---------------------------------------------------------------------------------------------
static void BuildCorePairlist(gpuContext gpu, eimage *tIMG, pairlist *simplePL,
                              const char* subimg, bool checkPairs)
{
  int i, j, k, m, ni, nj, isub;
  
  // Make local copies of the transformation matrices and cutoff
  double xcut2;
  if (strcmp(subimg, "QQ") == 0) {
    isub = 0;
    xcut2 = gpu->sim.es_cutoff + gpu->sim.skinnb;
  }
  else if (strcmp(subimg, "LJ") == 0) {
    isub = 1;
    xcut2 = gpu->sim.lj_cutoff + gpu->sim.skinnb;
  }
  xcut2 *= xcut2;
  dmat U, invU;
  U = ExtractTransformationMatrix(gpu, false);
  invU = ExtractTransformationMatrix(gpu, true);
  
  // Perform a nested loop over all atoms in the pencil cores, with re-imaging of
  // every interaction.
  pairlist pencilPL;
  int iLcore = tIMG->xpadding;
  pencilPL.ngbrCounts = (int*)calloc(gpu->sim.atoms + 1, sizeof(int));
  int nnb = 0;
  int* atmcover;
  atmcover = (int*)calloc(gpu->sim.atoms, sizeof(int));
  for (i = tIMG->ypadding; i < tIMG->ydimFull - tIMG->ypadding; i++) {
    for (j = tIMG->zpadding; j < tIMG->zdimFull; j++) {
      int cidx = (j * tIMG->ydimFull) + i;
      int ijlim = tIMG->cells[cidx].coreEnd;
      for (ni = tIMG->ypadding; ni < tIMG->ydimFull - tIMG->ypadding; ni++) {
        for (nj = tIMG->zpadding; nj < tIMG->zdimFull; nj++) {
          int cnidx = (nj * tIMG->ydimFull) + ni;
          if (cnidx > cidx) {
            continue;
          }
          int nijlim = tIMG->cells[cnidx].coreEnd;
          for (k = iLcore; k < ijlim; k++) {
            int katm = tIMG->cells[cidx].atoms[k];
            PMEFloat4 packet = tIMG->cells[cidx].crdq[k];
            double atmx = packet.x;
            double atmy = packet.y;
            double atmz = packet.z;
            atmcover[katm] = 1;
            for (m = iLcore; m < nijlim; m++) {
              if (cidx == cnidx && m >= k) {
                continue;
              }
              int matm = tIMG->cells[cnidx].atoms[m];
              packet = tIMG->cells[cnidx].crdq[m];
              double dx = atmx - packet.x;
              double dy = atmy - packet.y;
              double dz = atmz - packet.z;
              double ndx = U.data[0]*dx + U.data[3]*dy + U.data[6]*dz;
              double ndy = U.data[1]*dx + U.data[4]*dy + U.data[7]*dz;
              double ndz = U.data[2]*dx + U.data[5]*dy + U.data[8]*dz;
              ndx += (double)((ndx < -0.5) - (ndx >= 0.5));
              ndy += (double)((ndy < -0.5) - (ndy >= 0.5));
              ndz += (double)((ndz < -0.5) - (ndz >= 0.5));
              dx = invU.data[0]*ndx + invU.data[1]*ndy + invU.data[2]*ndz;
              dy = invU.data[3]*ndx + invU.data[4]*ndy + invU.data[5]*ndz;
              dz = invU.data[6]*ndx + invU.data[7]*ndy + invU.data[8]*ndz;
              double r2 = dx*dx + dy*dy + dz*dz;
              if (r2 < xcut2) {
                if (katm < matm) {
                  pencilPL.ngbrCounts[katm] += 1;
                }
                else {
                  pencilPL.ngbrCounts[matm] += 1;
                }
                nnb++;
              }
            }
          }
        }
      }
    }
  }
  FinalizePairlistAllocation(gpu, &pencilPL);  
  for (i = tIMG->ypadding; i < tIMG->ydimFull - tIMG->ypadding; i++) {
    for (j = tIMG->zpadding; j < tIMG->zdimFull; j++) {
      int cidx = (j * tIMG->ydimFull) + i;
      int ijlim = tIMG->cells[cidx].coreEnd;
      for (ni = tIMG->ypadding; ni < tIMG->ydimFull - tIMG->ypadding; ni++) {
        for (nj = tIMG->zpadding; nj < tIMG->zdimFull; nj++) {
          int cnidx = (nj * tIMG->ydimFull) + ni;
          if (cnidx > cidx) {
            continue;
          }
          int nijlim = tIMG->cells[cnidx].coreEnd;
          for (k = iLcore; k < ijlim; k++) {
            int katm = tIMG->cells[cidx].atoms[k];
            PMEFloat4 packet = tIMG->cells[cidx].crdq[k];
            double atmx = packet.x;
            double atmy = packet.y;
            double atmz = packet.z;
            atmcover[katm] = 1;
            for (m = iLcore; m < nijlim; m++) {
              if (cidx == cnidx && m >= k) {
                continue;
              }
              int matm = tIMG->cells[cnidx].atoms[m];
              packet = tIMG->cells[cnidx].crdq[m];
              double dx = atmx - packet.x;
              double dy = atmy - packet.y;
              double dz = atmz - packet.z;
              double ndx = U.data[0]*dx + U.data[3]*dy + U.data[6]*dz;
              double ndy = U.data[1]*dx + U.data[4]*dy + U.data[7]*dz;
              double ndz = U.data[2]*dx + U.data[5]*dy + U.data[8]*dz;
              ndx += (double)((ndx < -0.5) - (ndx >= 0.5));
              ndy += (double)((ndy < -0.5) - (ndy >= 0.5));
              ndz += (double)((ndz < -0.5) - (ndz >= 0.5));
              dx = invU.data[0]*ndx + invU.data[1]*ndy + invU.data[2]*ndz;
              dy = invU.data[3]*ndx + invU.data[4]*ndy + invU.data[5]*ndz;
              dz = invU.data[6]*ndx + invU.data[7]*ndy + invU.data[8]*ndz;
              double r2 = dx*dx + dy*dy + dz*dz;
              if (r2 < xcut2) {
                if (katm < matm) {
                  int plidx = pencilPL.ngbrCounts[katm];
                  pencilPL.ngbrList[plidx] = matm;
                  pencilPL.ngbrRanges[plidx] = sqrt(r2);
                  pencilPL.ngbrCounts[katm] = plidx + 1;
                }
                else {
                  int plidx = pencilPL.ngbrCounts[matm];
                  pencilPL.ngbrList[plidx] = katm;
                  pencilPL.ngbrRanges[plidx] = sqrt(r2);
                  pencilPL.ngbrCounts[matm] = plidx + 1;
                }
              }
            }
          }
        }
      }
    }
  }
  for (i = gpu->sim.atoms-1; i > 0; i--) {
    pencilPL.ngbrCounts[i] -= pencilPL.ngbrCounts[i] - pencilPL.ngbrCounts[i-1];
  }
  pencilPL.ngbrCounts[0] = 0;
  for (i = 0; i < gpu->sim.atoms; i++) {
    int2 filter = gpu->pbSubImageFilter->_pSysData[i];
    if (atmcover[i] == 0 && ((isub == 0 && filter.x >= 0) || (isub == 1 && filter.y >= 0))) {
      printf("BuildCorePairlist ::(%s) Atom %6d was not covered in the outer loop.\n",
             subimg, i);
    }    
  }

  // Check the individual atom neighbor counts
  if (checkPairs && MatchNeighbors(gpu, tIMG, &pencilPL, simplePL, "Pencil")) {
    printf("BuildCorePairlist :: (%s) Laboriously re-imaged honeycomb pair list passes.\n",
           subimg);
  }
  printf("BuildCorePairlist :: There are %d interactions.\n", nnb);

  // Free allocated memory
  free(atmcover);
  DestroyPairlist(&pencilPL);
  DestroyDmat(&U);
  DestroyDmat(&invU);
}  

//---------------------------------------------------------------------------------------------
// BuildHoneycombPairlist: build a neighbor list leveraging the honeycomb arrangement and its
//                         ghost atoms.  This is fast--there is no need to compute re-imaging,
//                         and there is a pre-defined list of interactions between hash cells
//                         that must be covered.  However, it needs to WORK, and that's what
//                         this function can test when the results are compared with those of
//                         BuildSimplePairlist().
//
// This is a debugging function.
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                object counts
//   tIMG:        the expanded sub-image representation
//   simplePL:    a reference pair list computed by the simplest possible method
//   subimg:      code for the sub-image to allocate (different cutoffs may imply different
//                levels of hash cell padding in Y and Z)
//   checkPairs:  flag to have all apir interactions checked against the reference
//---------------------------------------------------------------------------------------------
static pairlist BuildHoneycombPairlist(gpuContext gpu, eimage *tIMG, pairlist *simplePL,
                                       const char* subimg, bool checkPairs)
{
  int i, j, k;

  // Get the sub-image of interest
  int nNgbrPencils;
  int *ngbrPencils;
  double xcut2;
  nNgbrPencils = gpu->sim.nNgbrPencils;
  ngbrPencils = gpu->pbPencils->_pSysData;
  if (strcmp(subimg, "QQ") == 0) {
    xcut2 = gpu->sim.es_cutoff + gpu->sim.skinnb;
  }
  else if (strcmp(subimg, "LJ") == 0) {
    xcut2 = gpu->sim.lj_cutoff + gpu->sim.skinnb;
  }
  xcut2 *= xcut2;

  // Prepare an inventory of everything in the honeycomb expanded image.
  pairlist hcmbPL;
  hcmbPL.ngbrCounts = (int*)calloc(gpu->sim.atoms + 1, sizeof(int));
  hcmbPL.ngbrRanges = NULL;
  int nnb = 0;
  for (i = tIMG->ypadding; i < tIMG->ydimFull - tIMG->ypadding; i++) {
    for (j = tIMG->zpadding; j < tIMG->zdimFull; j++) {
      
      // Loop over all pencils that the present one could interact with
      for (k = 0; k < nNgbrPencils; k++) {
        int ni  = i + ((ngbrPencils[k] >> 8) & 0xff) - 128;
        int nj  = j + (ngbrPencils[k] & 0xff) - 128;
        int n2i = i + ((ngbrPencils[k] >> 24) & 0xff) - 128;
        int n2j = j + ((ngbrPencils[k] >> 16) & 0xff) - 128;
        if (gpu->bPencilYShift) {
          ni  += (((j - tIMG->zpadding)*(j - nj)) & 0x1);
          n2i += (((j - tIMG->zpadding)*(j - n2j)) & 0x1);
        }
        nnb = PencilNestedLoop(tIMG, &hcmbPL, ni, nj, n2i, n2j, nnb, xcut2);
      }
    }
  }
  printf("BuildHoneycombPairlist :: There are %d interactions.\n", nnb);
  FinalizePairlistAllocation(gpu, &hcmbPL);
  nnb = 0;
  for (i = tIMG->ypadding; i < tIMG->ydimFull - tIMG->ypadding; i++) {
    for (j = tIMG->zpadding; j < tIMG->zdimFull; j++) {

      // Loop over all pencils that the present one could interact with
      for (k = 0; k < nNgbrPencils; k++) {
        int ni  = i + ((ngbrPencils[k] >> 8) & 0xff) - 128;
        int nj  = j + (ngbrPencils[k] & 0xff) - 128;
        int n2i = i + ((ngbrPencils[k] >> 24) & 0xff) - 128;
        int n2j = j + ((ngbrPencils[k] >> 16) & 0xff) - 128;
        if (gpu->bPencilYShift) {
          ni  += (((j - tIMG->zpadding)*(j - nj)) & 0x1);
          n2i += (((j - tIMG->zpadding)*(j - n2j)) & 0x1);
        }
        nnb = PencilNestedLoop(tIMG, &hcmbPL, ni, nj, n2i, n2j, nnb, xcut2);
      }
    }
  }
  for (i = gpu->sim.atoms - 1; i > 0; i--) {
    hcmbPL.ngbrCounts[i] -= hcmbPL.ngbrCounts[i] - hcmbPL.ngbrCounts[i-1];
  }
  hcmbPL.ngbrCounts[0] = 0;

  // Check the individual atom neighbor counts
  if (checkPairs && MatchNeighbors(gpu, tIMG, &hcmbPL, simplePL, "Honeycomb")) {
    printf("BuildHoneycombPairlist :: (%s) Ghost atom honeycomb pair list passes.\n", subimg);
  }

  // Free allocated memory
  return hcmbPL;
}

//---------------------------------------------------------------------------------------------
// UnfoldExpansion: unfold the expansion of the coordinates using instructions written on the
//                  GPU.
//
// This is a debugging function.
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                object counts
//   tIMG:        the expanded sub-image representation
//   limits:      the limits of each cell in the primary image (lower index, lower cap index,
//                upper cap index, upper index)
//   subimg:      code for the sub-image to allocate (different cutoffs may imply different
//                levels of hash cell padding in Y and Z)
//---------------------------------------------------------------------------------------------
static void UnfoldExpansion(gpuContext gpu, eimage *tIMG, int* limits, const char* subimg)
{
  int i, isub;

  // Extract the real-space transformation matrix
  dmat invU;
  invU = ExtractTransformationMatrix(gpu, true);

  // Make shortcuts
  int cellSpace = tIMG->cellspace;
  
  int ninsr;
  int2 *insrList;
  int *transList;
  if (strcmp(subimg, "QQ") == 0) {
    isub = 0;
    ninsr = gpu->pbExpansionInsrCount->_pSysData[0];
    insrList = gpu->pbQQExpansions->_pSysData;
    transList = gpu->pbQQTranslations->_pSysData;
  }
  else if (strcmp(subimg, "LJ") == 0) {
    isub = 1;
    ninsr = gpu->pbExpansionInsrCount->_pSysData[1];
    insrList = gpu->pbLJExpansions->_pSysData;
    transList = gpu->pbLJTranslations->_pSysData;
  }
  int nNoTrans = 0;
  const int noTransMask = (1 << 16) + (1 << 8) + 1;
  for (i = 0; i < ninsr; i++) {

    // Record the number of particles that should not move
    if (insrList[i].x == insrList[i].y) {
      if (transList[i] != noTransMask) {
        printf("UnfoldExpansion :: Instruction %6d moves a particle but keeps it in the same "
               "place.\n", i);
      }
      nNoTrans++;
    }
    
    // Check the result against the known cell limits.  These "cellpos" like indices
    // are absolute with respect to the built-in padding that each cell has to
    // accommodate ghost atoms in the X direction, unlike the cellpos values in
    // kMapSubImageExpansion_kernel().
    int cOrigIdx = insrList[i].x / cellSpace;
    int cDestIdx = insrList[i].y / cellSpace;
    int cOrigPos = insrList[i].x - (cellSpace * cOrigIdx);
    int cDestPos = insrList[i].y - (cellSpace * cDestIdx);
    bool problem = false;
    if (cOrigPos < tIMG->cells[cOrigIdx].lcapStart ||
        cOrigPos >= tIMG->cells[cOrigIdx].hcapEnd) {
      printf("UnfoldExpansion :: Instruction %6d takes invalid atom %6d to position %6d.\n",
             i, insrList[i].x, insrList[i].y);
      problem = true;
    }
    if (cDestPos < tIMG->cells[cDestIdx].lcapStart ||
        cDestPos >= tIMG->cells[cDestIdx].hcapEnd) {
      printf("UnfoldExpansion :: Instruction %6d takes atom %6d to invalid position %6d.\n",
             i, insrList[i].x, insrList[i].y);
      problem = true;
    }
    if (problem) {
      int zcellOrig = cOrigIdx / tIMG->ydimFull;
      int ycellOrig = cOrigIdx - (tIMG->ydimFull * zcellOrig);
      int zcellDest = cDestIdx / tIMG->ydimFull;
      int ycellDest = cDestIdx - (tIMG->ydimFull * zcellDest);
      printf("                   This wants cell (%2d %2d, %3d) -> (%2d %2d, %3d)\n",
             ycellOrig, zcellOrig, cOrigPos, ycellDest, zcellDest, cDestPos);
      printf("                   Cell grid dimensions %2d x %2d (X padding = %2d)\n",
             tIMG->ydimFull, tIMG->zdimFull, tIMG->xpadding);
      printf("                   Valid indices %6d - %6d (orig) and %6d - %6d (dest)\n",
             tIMG->cells[cOrigIdx].lcapStart + (cellSpace * cOrigIdx),
             tIMG->cells[cOrigIdx].hcapEnd + (cellSpace * cOrigIdx),
             tIMG->cells[cDestIdx].lcapStart + (cellSpace * cDestIdx),
             tIMG->cells[cDestIdx].hcapEnd + (cellSpace * cDestIdx));
      tIMG->crdq[insrList[i].y].x = 0.0;
      tIMG->crdq[insrList[i].y].y = 0.0;
      tIMG->crdq[insrList[i].y].z = 0.0;
      tIMG->crdq[insrList[i].y].w = 0.0;
      continue;
    }
    
    // Copy and move the atoms
    PMEFloat4 packet = tIMG->crdq[insrList[i].x];
    int xtrans = transList[i] & 0xff;
    int ytrans = (transList[i] >> 8) & 0xff;
    int ztrans = (transList[i] >> 16);
    double xmove = (double)((xtrans == 2) - (xtrans == 0));
    double ymove = (double)((ytrans == 2) - (ytrans == 0));
    double zmove = (double)((ztrans == 2) - (ztrans == 0));
    packet.x += invU.data[0]*xmove + invU.data[1]*ymove + invU.data[2]*zmove;
    packet.y += invU.data[3]*xmove + invU.data[4]*ymove + invU.data[5]*zmove;
    packet.z += invU.data[6]*xmove + invU.data[7]*ymove + invU.data[8]*zmove;

    // Commit the result
    tIMG->crdq[insrList[i].y] = packet;
  }

  // Check on the number of particles that did not change array positions
  if (isub == 0 && gpu->sim.nQQatoms != nNoTrans) {
    printf("UnfoldExpansion :: With %d atoms in the primary image, %d did not change grid "
           "indices.\n", gpu->sim.nQQatoms, nNoTrans);
  }
  if (isub == 1 && gpu->sim.nLJatoms != nNoTrans) {
    printf("UnfoldExpansion :: With %d atoms in the primary image, %d did not change grid "
           "indices.\n", gpu->sim.nLJatoms, nNoTrans);
  }
  
  // Free allocated memory
  DestroyDmat(&invU);
}

//---------------------------------------------------------------------------------------------
// ExpandGpuImage: propagate the primary image using pre-computed instructions written on the
//                 GPU.  This follows on the work of ExpandCpuImage() above.
//
// This is a debugging function.
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                object counts
//   cpuIMG:      the expanded sub-image representation, computed on the CPU with alternative
//                code
//   limits:      the limits of each cell in the primary image (lower index, lower cap index,
//                upper cap index, upper index)
//   primaryCrd:  the coordinates of all atoms in the primary image
//   subimg:      code for the sub-image to allocate (different cutoffs may imply different
//                levels of hash cell padding in Y and Z)
//---------------------------------------------------------------------------------------------
static void ExpandGpuImage(gpuContext gpu, eimage *cpuIMG, eimage *gpuIMG, int* limits,
                           double* primaryCrd, const char* subimg)
{
  int i, j, k, isub;

  // Check on the number of atoms outside the primary image core in the CPU-built
  // expanded representation.
  int noncore = 0;
  int allcore = 0;
  for (i = 0; i < cpuIMG->ydimFull; i++) {
    for (j = 0; j < cpuIMG->zdimFull; j++) {
      int cidx = (j * cpuIMG->ydimFull) + i;
      if (i >= cpuIMG->ypadding && i < cpuIMG->ydimFull - cpuIMG->ypadding &&
          j >= cpuIMG->zpadding) {
        noncore += (cpuIMG->cells[cidx].coreStart - cpuIMG->cells[cidx].lcapStart) +
                   (cpuIMG->cells[cidx].hcapEnd - cpuIMG->cells[cidx].coreEnd);
        allcore += cpuIMG->cells[cidx].coreEnd - cpuIMG->cells[cidx].coreStart;
      }
      else {
        noncore += (cpuIMG->cells[cidx].hcapEnd - cpuIMG->cells[cidx].lcapStart);
      }
    }
  }
  printf("ExpandGpuImage :: (%s) There are %d atoms and %d ghost atoms.\n", subimg, allcore,
         noncore);
  
  // Step 1: unpack the primary image from coordinates, cell limits,
  //         and filter instructions onto the grid
  LoadPrimaryImage(gpu, gpuIMG, limits, primaryCrd, subimg);

  // Step 2: unpack the ghost atoms from the primary image already on
  //         the grid and the expansion instructions.  This will check
  //         cell limits to ensure that no expected bounds are violated.
  UnfoldExpansion(gpu, gpuIMG, limits, subimg);

  // Compare the expansion guided by GPU instructions to
  // that computed with independent code on the CPU
  bool problem = false;
  if (strcmp(subimg, "QQ") == 0) {
    isub = 0;
  }
  else if (strcmp(subimg, "LJ") == 0) {
    isub = 1;
  }
  for (i = 0; i < cpuIMG->ydimFull; i++) {
    for (j = 0; j < cpuIMG->zdimFull; j++) {
      int cidx = (j * cpuIMG->ydimFull) + i;
      PMEFloat4 *cpuCrdq = cpuIMG->cells[cidx].crdq;
      PMEFloat4 *gpuCrdq = gpuIMG->cells[cidx].crdq;
      for (k = cpuIMG->cells[cidx].lcapStart; k < cpuIMG->cells[cidx].hcapEnd; k++) {
        PMEFloat dx = cpuCrdq[k].x - gpuCrdq[k].x;
        PMEFloat dy = cpuCrdq[k].y - gpuCrdq[k].y;
        PMEFloat dz = cpuCrdq[k].z - gpuCrdq[k].z;
        PMEFloat dw;
        if (isub == 0) {
          dw = cpuCrdq[k].w - gpuCrdq[k].w;
        }
        else if (isub == 1) {
          FloatShift cxtyp, gxtyp;
          cxtyp.f = cpuCrdq[k].w;
          gxtyp.f = gpuCrdq[k].w;
          dw = cxtyp.ui - gxtyp.ui;
        }
        if (fabs(dx) > 1.0e-4 || fabs(dy) > 1.0e-4 || fabs(dz) > 1.0e-4 || fabs(dw) > 1.0e-4) {
          printf("ExpandGpuImage :: (%s) Mismatch in atom %3d of cell %3d %3d.\n"
                 "                       [ %9.4f %9.4f %9.4f %9.4f ] ->\n"
                 "                       [ %9.4f %9.4f %9.4f %9.4f ]\n", subimg, k, i, j,
                 cpuCrdq[k].x, cpuCrdq[k].y, cpuCrdq[k].z, cpuCrdq[k].w, gpuCrdq[k].x,
                 gpuCrdq[k].y, gpuCrdq[k].z, gpuCrdq[k].w);
          int ysrc = i +
                     (cpuIMG->ydimPrimary * ((i < cpuIMG->ypadding) -
                                             (i >= cpuIMG->ypadding + cpuIMG->ydimPrimary)));
          int zsrc = j + ((j < cpuIMG->zpadding) * cpuIMG->zdimPrimary);
          int xsrc = k + ((k < cpuIMG->cells[cidx].coreStart) -
                          (k >= cpuIMG->cells[cidx].coreEnd))*(cpuIMG->cells[cidx].coreEnd -
                                                               cpuIMG->cells[cidx].coreStart);
          int srcidx = (zsrc * cpuIMG->ydimFull) + ysrc;
          printf("                       Location: [ y z ] = %3d %3d  [ x ] = %3d (idx %6d) "
                 "[ %6d %6d %6d %6d ]\n"
                 "                       Origin:   [ y z ] = %3d %3d  [ x ] = %3d (idx %6d) "
                 "[ %6d %6d %6d %6d ]\n", i, j, k, (cidx * cpuIMG->cellspace) + k,
                 cpuIMG->cells[cidx].lcapStart + (cidx * cpuIMG->cellspace),
                 cpuIMG->cells[cidx].coreStart + (cidx * cpuIMG->cellspace),
                 cpuIMG->cells[cidx].coreEnd + (cidx * cpuIMG->cellspace),
                 cpuIMG->cells[cidx].hcapEnd + (cidx * cpuIMG->cellspace), ysrc, zsrc, xsrc,
                 (srcidx * cpuIMG->cellspace) + xsrc,
                 cpuIMG->cells[srcidx].lcapStart + (srcidx * cpuIMG->cellspace),
                 cpuIMG->cells[srcidx].coreStart + (srcidx * cpuIMG->cellspace),
                 cpuIMG->cells[srcidx].coreEnd + (srcidx * cpuIMG->cellspace),
                 cpuIMG->cells[srcidx].hcapEnd + (srcidx * cpuIMG->cellspace));
          problem = true;
        }
      }
    }
  }
  if (problem == false) {
    printf("ExpandGpuImage :: (%s) Comparison of CPU and GPU-derived expansions passes.\n",
           subimg);
  }
}

//---------------------------------------------------------------------------------------------
// CheckGpuExpandedImage: the particle expansion has now been computed on both the CPU and the
//                        GPU, then executed on the CPU using either set of instructions (any
//                        disagreement will have been reported).  The expansion has also been
//                        executed on the GPU, and can be downloaded for comparison in this
//                        function.
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                object counts
//   tIMG:        the expanded sub-image representation
//   limits:      the limits of each cell in the primary image (lower index, lower cap index,
//                upper cap index, upper index)
//   subimg:      code for the sub-image to allocate (different cutoffs may imply different
//                levels of hash cell padding in Y and Z)
//---------------------------------------------------------------------------------------------
static void CheckGpuExpandedImage(gpuContext gpu, eimage *tIMG, int* limits,
                                  const char* subimg)
{
  int i, j, k, nylim, nzlim;
  int *limptr, *gpulimptr, *idptr;
  bool isLJ;
  PMEFloat4 *atmptr;

  // Get the GPU cell limits, in case they are needed and not yet downloaded
  gpu->pbHcmbSubImgCellLimits->Download();

  // Set pointers
  if (strcmp(subimg, "QQ") == 0) {
    gpu->pbHcmbQQParticles->Download();
    gpu->pbHcmbQQIdentities->Download();
    atmptr = gpu->pbHcmbQQParticles->_pSysData;
    idptr = gpu->pbHcmbQQIdentities->_pSysData;
    limptr = limits;
    gpulimptr = gpu->pbHcmbSubImgCellLimits->_pSysData;
    isLJ = false;
  }
  else if (strcmp(subimg, "LJ") == 0) {
    gpu->pbHcmbLJParticles->Download();
    gpu->pbHcmbLJIdentities->Download();
    atmptr = gpu->pbHcmbLJParticles->_pSysData;
    idptr = gpu->pbHcmbLJIdentities->_pSysData;
    limptr = &limits[4 * gpu->sim.npencils];
    gpulimptr = &gpu->pbHcmbSubImgCellLimits->_pSysData[4 * gpu->sim.npencils];
    isLJ = true;
  }
  nylim = tIMG->ydimFull;
  nzlim = tIMG->zdimFull;
  for (i = 0; i < nylim; i++) {
    int ysrc = i + ((i < tIMG->ypadding) -
                    (i >= tIMG->ypadding + tIMG->ydimPrimary)) * tIMG->ydimPrimary;
    ysrc -= tIMG->ypadding;
    for (j = 0; j < nzlim; j++) {
      int zsrc = j + (j < tIMG->zpadding) * tIMG->zdimPrimary;
      zsrc -= tIMG->zpadding;
      int srcidx = (zsrc * tIMG->ydimPrimary) + ysrc;
      int cellidx = (j * nylim) + i;
      int lim1 = limptr[srcidx];
      int lim2 = limptr[srcidx +       gpu->sim.npencils];
      int lim3 = limptr[srcidx + (2 * gpu->sim.npencils)];
      int lim4 = limptr[srcidx + (3 * gpu->sim.npencils)];
      lim2 += (lim2 < 0) * (lim1 + 1);
      lim3 += (lim3 < 0) * (lim4 + 1);
      int lbound = (cellidx * tIMG->cellspace) + tIMG->xpadding - (lim4 - lim3);
      int hbound = (cellidx * tIMG->cellspace) + tIMG->xpadding + lim4 + lim2 - 2*lim1;

      // Verify the bounds
      if (hbound - lbound != tIMG->cells[cellidx].hcapEnd - tIMG->cells[cellidx].lcapStart) {
        printf("CheckGpuExpandedImage :: (%s) Bounds mismatch [ %6d %6d ] vs [ %6d %6d ] in "
               "cell [ %2d %2d ]\n", subimg, lbound, hbound,
               (cellidx * tIMG->cellspace) + tIMG->cells[cellidx].lcapStart,
               (cellidx * tIMG->cellspace) + tIMG->cells[cellidx].hcapEnd, i, j);
        printf("      CPU [ %6d %6d %6d %6d ]    GPU [ %6d %6d %6d %6d ]\n",
               lim1, lim2, lim3, lim4, gpulimptr[srcidx],
               gpulimptr[srcidx + gpu->sim.npencils],
               gpulimptr[srcidx + (2 * gpu->sim.npencils)],
               gpulimptr[srcidx + (3 * gpu->sim.npencils)]);
        printf("     Cell [ %6d %6d %6d %6d ]   tIMG [ %6d %6d %6d %6d ]\n",
               (cellidx * tIMG->cellspace) + tIMG->xpadding,
               (cellidx * tIMG->cellspace) + tIMG->xpadding + lim2 - lim1,
               (cellidx * tIMG->cellspace) + tIMG->xpadding + lim3 - lim1,
               (cellidx * tIMG->cellspace) + tIMG->xpadding + lim4 - lim1,
               (cellidx * tIMG->cellspace) + tIMG->cells[cellidx].lcapStart,
               (cellidx * tIMG->cellspace) + tIMG->cells[cellidx].coreStart,
               (cellidx * tIMG->cellspace) + tIMG->cells[cellidx].coreEnd,
               (cellidx * tIMG->cellspace) + tIMG->cells[cellidx].hcapEnd);
               
      }
      
      // Verify agreement between the bounds of the extended cell
      for (k = lbound; k < hbound; k++) {
        PMEFloat4 cpuAtom = tIMG->crdq[k];
        PMEFloat4 gpuAtom = atmptr[k];
        if (fabs(cpuAtom.x - gpuAtom.x) > 1.0e-4 || fabs(cpuAtom.y - gpuAtom.y) > 1.0e-4 ||
            fabs(cpuAtom.z - gpuAtom.z) > 1.0e-4 || fabs(cpuAtom.w - gpuAtom.w) > 1.0e-4) {
          printf("CheckGpuExpandedImage :: (%s) Mismatch in atom [ %3d %3d %5d -> %6d ] "
                 "coordinates and property:\n", subimg, i, j,
                 k - ((cellidx * tIMG->cellspace) + tIMG->xpadding), k);
          if (isLJ) {
            FloatShift fs;
            fs.f = cpuAtom.w;
#ifdef use_DPFP
            printf("      CPU : coord %12.7f %12.7f %12.7f, prop %4llu\n", cpuAtom.x,
                   cpuAtom.y, cpuAtom.z, fs.ui);
#else
            printf("      CPU : coord %12.7f %12.7f %12.7f, prop %4u\n", cpuAtom.x,
                   cpuAtom.y, cpuAtom.z, fs.ui);
#endif
            fs.f = gpuAtom.w;
#ifdef use_DPFP
            printf("      GPU : coord %12.7f %12.7f %12.7f, prop %4llu\n", gpuAtom.x,
                   gpuAtom.y, gpuAtom.z, fs.ui);
#else
            printf("      GPU : coord %12.7f %12.7f %12.7f, prop %4u\n", gpuAtom.x,
                   gpuAtom.y, gpuAtom.z, fs.ui);
#endif
          }
          else {
            printf("      CPU : coord %12.7f %12.7f %12.7f, prop %12.7f\n", cpuAtom.x,
                   cpuAtom.y, cpuAtom.z, cpuAtom.w);
            printf("      GPU : coord %12.7f %12.7f %12.7f, prop %12.7f\n", gpuAtom.x,
                   gpuAtom.y, gpuAtom.z, gpuAtom.w);
          }
        }
        if (tIMG->absid[k] != idptr[k]) {
          printf("CheckGpuExpandedImage :: (%s) Mismatch in atom [ %3d %3d %5d - %6d ] "
                 "identity:\n", subimg, i, j,
                 k - ((cellidx * tIMG->cellspace) + tIMG->xpadding), k);
          printf("      CPU : %6d    GPU : %6d\n", tIMG->absid[k], idptr[k]);
        }
      }
    }
  }
}

//---------------------------------------------------------------------------------------------
// CheckNTStencils: re-derive the Neutral Territory stencil import regions (called NT regions
//                  in kNTX.h and kNeighborList.cu) on the CPU.  Critical data for pencil
//                  relationships should already have been loaded into the CPU-side data of
//                  the gpuContext struct and will be accessed from there without downloading
//                  back from the GPU (it should not have changed on the GPU).
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                numerous gpuBuffer attributes
//   tIMG:        electrostatic or Lennard-Jones sub-image for which to re-derive NT limits
//   limits:      the limits of cells and capping regions for each Honeycomb hash cell
//   ntrAssigned: vector to track which NT regions found on the GPU have been assigned by the
//                CPU.  This assignment must take place over two iterations of this function,
//                as the NT regions list and the pair list are unified between electrostatic
//                and Lennard-Jones interactions.
//   subimg:      code for the sub-image to allocate (different cutoffs may imply different
//                levels of hash cell padding in Y and Z)
//---------------------------------------------------------------------------------------------
static void CheckNTStencils(gpuContext gpu, eimage *tIMG, int* limits, imat *ntrAssigned,
                            const char* subimg)
{
  int i, j, k, m, n;

  // Get the lists of pencil relationships.  Here, ntPencils will again be the list of 11
  // imported hash cells relative to some home cell, ntHandShakes will be the interactions
  // among each of two partitions of the 11 hash cells, and ntRanges the appropriate cutoff
  // less the minimum separation between any two points within each pair of hash cells.
  imat ntImports, ntHandShakes, ntCollection, climits, segpop, segloc, seglimits;
  dmat ntRanges, U, invU;
  ntImports = CreateImat(11, 2);
  ntCollection = CreateImat(11, 3);
  climits = CreateImat(11, 4);
  segpop = CreateImat(11, 64);
  segloc = CreateImat(11, 1024);
  seglimits = CreateImat(11, 64);
  for (i = 0; i < 11; i++) {
    ntImports.map[i][0] = gpu->pbNTImports->_pSysData[2*i    ];
    ntImports.map[i][1] = gpu->pbNTImports->_pSysData[2*i + 1];
  }
  ntHandShakes = CreateImat(2 * GRID, 2);
  ntRanges = CreateDmat(1, 2 * GRID);
  for (i = 0; i < gpu->sim.nNgbrPencils; i++) {
    ntHandShakes.map[i     ][0] = gpu->pbSegmentRelays->_pSysData[i     ];
    ntHandShakes.map[i     ][1] = gpu->pbSegmentRelays->_pSysData[i + 32];
    ntHandShakes.map[i + 32][0] = gpu->pbSegmentRelays->_pSysData[i + 64];
    ntHandShakes.map[i + 32][1] = gpu->pbSegmentRelays->_pSysData[i + 96];
    ntRanges.map[0][i     ] = gpu->pbPencilRanges->_pSysData[i     ];
    ntRanges.map[0][i + 32] = gpu->pbPencilRanges->_pSysData[i + 32];
  }
  U = ExtractTransformationMatrix(gpu, false);

  // Loop over all NT regions
  int gblntmap = 0;
  for (i = 0; i < gpu->sim.npencils; i++) {

    // Identify the cell of interest
    int srcz = i / tIMG->ydimPrimary;
    int srcy = i - (srcz * tIMG->ydimPrimary);
    int celly = srcy + tIMG->ypadding;
    int cellz = srcz + tIMG->zpadding;
    int cellidx = (tIMG->ydimFull * cellz) + celly;
    
    // Enumerate the imported cells
    for (j = 0; j < 11; j++) {
      ntCollection.map[j][0] = celly + ntImports.map[j][0];
      ntCollection.map[j][1] = cellz + ntImports.map[j][1];
      if (gpu->bPencilYShift) {
        ntCollection.map[j][0] += ((srcz * ntImports.map[j][1]) & 0x1);
      }
    }

    // Get the appropriate limits and atom count
    int ntAtomCount = 0;
    for (j = 0; j < 11; j++) {
      int pcy = ntCollection.map[j][0];
      int pcz = ntCollection.map[j][1];
      int psy = pcy + (((pcy < tIMG->ypadding) -
                        (pcy >= tIMG->ydimPrimary + tIMG->ypadding)) * tIMG->ydimPrimary);
      int psz = pcz + ((pcz < tIMG->zpadding) * tIMG->zdimPrimary);
      psy -= tIMG->ypadding;
      psz -= tIMG->zpadding;
      int psrc  = (psz * tIMG->ydimPrimary) + psy;
      int pcell = (pcz * tIMG->ydimFull) + pcy;
      ntCollection.map[j][2] = pcell;
      climits.map[j][0] = tIMG->cells[pcell].lcapStart;
      climits.map[j][1] = tIMG->cells[pcell].coreStart;
      climits.map[j][2] = tIMG->cells[pcell].coreEnd;
      climits.map[j][3] = tIMG->cells[pcell].hcapEnd;
      ntAtomCount += climits.map[j][3] - climits.map[j][0];
    }

    // Determine whether the atoms will all fit, or if multiple sections are needed.
    int nsegment = 1;
    if (ntAtomCount > NT_MAX_ATOMS) {
      nsegment = (ntAtomCount + NT_MAX_ATOMS - 1) / NT_MAX_ATOMS;
      bool success = false;
      while (!success) {
        success = true;
        PMEFloat invTWidth = (PMEFloat)nsegment / (1.0 + 2.0 * gpu->sim.LowerCapBoundary);
        SetIVec(segpop.data, segpop.row * segpop.col, 0);
        for (j = 0; j < 11; j++) {
          int pcell = ntCollection.map[j][2];
          for (k = tIMG->cells[pcell].lcapStart; k < tIMG->cells[pcell].hcapEnd; k++) {
            PMEFloat4 packet = tIMG->cells[pcell].crdq[k];
            PMEFloat atmx = packet.x;
            PMEFloat atmy = packet.y;
            PMEFloat atmz = packet.z;
            PMEFloat rbx = U.data[0]*atmx + U.data[3]*atmy + U.data[6]*atmz;
            int irbx = (rbx + gpu->sim.LowerCapBoundary) * invTWidth;
            segloc.map[j][k - tIMG->cells[pcell].lcapStart] = irbx;
          }
          for (k = 0; k < nsegment; k++) {
            int lbound = -1;
            int hbound = -1;
            for (m = tIMG->cells[pcell].lcapStart; m < tIMG->cells[pcell].hcapEnd; m++) {
              if (segloc.map[j][m - tIMG->cells[pcell].lcapStart] == k) {
                if (lbound < 0) {
                  lbound = m;
                }
                hbound = m + 1;
              }
            }
            segpop.map[j][k] = hbound - lbound;
          }
        }
        for (j = 0; j < nsegment; j++) {
          int tspop = 0;
          for (k = 0; k < 11; k++) {
            tspop += segpop.map[k][j];
          }
          if (tspop*2 > NT_MAX_ATOMS) {
            success = false;
          }
        }
        if (!success) {
          nsegment++;
        }
      }
    }

    // How many segment neighbors are needed to cover all of the interactions?
    int nspan;
    if (nsegment == 1) {
      nspan = 0;
    }
    else {
      nspan = ceil(gpu->sim.LowerCapBoundary * (double)nsegment /
                   (1.0 + (2.0 * gpu->sim.LowerCapBoundary)));
    }

    // Determine the compartmentalization of all hash cells.  Where does each
    // segment start and stop?  Store the results in terms of absolute
    // positions in the very large array of particles held in the eimage struct.
    SetIVec(seglimits.data, seglimits.row * seglimits.col, 0);
    for (j = 0; j < 11; j++) {
      int pcell = ntCollection.map[j][2];
      seglimits.map[j][0] = (tIMG->cellspace * pcell) + tIMG->cells[pcell].lcapStart;
      if (nsegment == 1) {
        seglimits.map[j][1] = (tIMG->cellspace * pcell) + tIMG->cells[pcell].hcapEnd;
      }
      else {
        for (k = 0; k < nsegment; k++) {
          for (m = tIMG->cells[pcell].lcapStart; m < tIMG->cells[pcell].hcapEnd; m++) {
            if (segloc.map[j][m - tIMG->cells[pcell].lcapStart] == k) {
              seglimits.map[j][k + 1] = (tIMG->cellspace * pcell) + m + 1;
            }
          }
        }
      }
    }

    // Allocate space for the neutral territory maps
    int maxntreg = nsegment * (nspan + (nspan == 0));
    tIMG->cells[cellidx].ntmaps = (int*)calloc(maxntreg * NTMAP_SPACE, sizeof(int));
    
    // Determine the limits for each pair of segments
    int ntmap[NTMAP_SPACE];
    SetIVec(ntmap, NTMAP_SPACE, 0);
    for (j = 0; j < nsegment - nspan; j++) {
      for (k = j + (nsegment > 1); k <= j + nspan; k++) {

        // Derive the limits for a prefix sum over the 22 sectors
        int atomPrefixSum = 0;
        if (nsegment == 1 || k == j + 1) {
          for (m = 0; m < 11; m++) {
            int rhalf = 0;
            if (nsegment == 1) {
              int rsize = seglimits.map[m][k + 1] - seglimits.map[m][j];
              rhalf = ((rsize + GRID_BITS_MASK) / GRID) / 2;
              rhalf *= GRID;
            }
            ntmap[m +  44] = atomPrefixSum;
            atomPrefixSum += seglimits.map[m][j + 1] - seglimits.map[m][j];
            ntmap[m +  55] = atomPrefixSum - (nsegment == 1)*rhalf;
            ntmap[m +  66] = atomPrefixSum;
            if (nsegment > 1) {
              atomPrefixSum += seglimits.map[m][k + 1] - seglimits.map[m][k];
            }
            ntmap[m +  77] = atomPrefixSum - (nsegment == 1)*rhalf;
          }
        }
        else {
          for (m = 0; m < 22; m++) {
            ntmap[m +  44] = atomPrefixSum;
            if (m < 11) {
              atomPrefixSum += seglimits.map[m][j + 1] - seglimits.map[m][j];
            }
            else {
              atomPrefixSum += seglimits.map[m - 11][k + 1] - seglimits.map[m - 11][k];
            }
            ntmap[m +  66] = atomPrefixSum;
          }
        }
        for (m = 0; m < 11; m++) {
          int readJs = seglimits.map[m][j];
          int readJf = seglimits.map[m][j + 1];
          int readKs = seglimits.map[m][k];
          int readKf = seglimits.map[m][k + 1];
          int pcell = ntCollection.map[m][2];
          
          // Determine the core atoms in segments J and K,
          // again expressing numbers as absolute indices
          int coreJs, coreJf, coreKs, coreKf;
          int coreS = tIMG->cells[pcell].coreStart + (tIMG->cellspace * pcell);
          int coreF = tIMG->cells[pcell].coreEnd + (tIMG->cellspace * pcell);
          if (readJf < coreS || readJs >= coreF) {
            coreJs = readJs;
            coreJf = readJs;
          }
          else if (readJs < coreS && readJf >= coreF) {
            coreJs = coreS;
            coreJf = coreF;
          }
          else if (readJs < coreS && readJf < coreF) {
            coreJs = coreS;
            coreJf = readJf;
          }
          else if (readJs >= coreS && readJs < coreF && readJf < coreF) {
            coreJs = readJs;
            coreJf = readJf;
          }
          else if (readJs >= coreS && readJs < coreF && readJf >= coreF) {
            coreJs = readJs;
            coreJf = coreF;
          }
          if (readKf < coreS || readKs >= coreF) {
            coreKs = readKs;
            coreKf = readKs;
          }
          else if (readKs < coreS && readKf >= coreF) {
            coreKs = coreS;
            coreKf = coreF;
          }
          else if (readKs < coreS && readKf < coreF) {
            coreKs = coreS;
            coreKf = readKf;
          }
          else if (readKs >= coreS && readKs < coreF && readKf < coreF) {
            coreKs = readKs;
            coreKf = readKf;
          }
          else if (readKs >= coreS && readKs < coreF && readKf >= coreF) {
            coreKs = readKs;
            coreKf = coreF;
          }
          
          // Load limits
          if (nsegment == 1) {
            ntmap[m      ] = ntmap[m +  44] + coreJs - readJs;
            ntmap[m +  11] = ntmap[m      ];
            ntmap[m +  22] = ntmap[m +  44] + coreJf - readJs;
            ntmap[m +  33] = ntmap[m +  22];
            int rsize = readKf - readJs;
            int rhalf = ((rsize + GRID_BITS_MASK) / GRID) / 2;
            rhalf *= GRID;
            readKs = readKf - rhalf;
            readJf = readJs + (rsize - rhalf);
            ntmap[m +  88] = readJs;
            ntmap[m +  99] = readJf;
            ntmap[m + 110] = readJf;
            ntmap[m + 121] = readKf;
          }
          else {

            // One of the hash cell segments may have zero atoms in it--find those cases
            // and deal with them appropriately.
            if (seglimits.map[m][j + 1] - seglimits.map[m][j] == 0) {
              seglimits.map[m][j    ] = seglimits.map[m][k];
              seglimits.map[m][j + 1] = seglimits.map[m][k];
            }
            if (seglimits.map[m][k + 1] - seglimits.map[m][k] == 0) {
              seglimits.map[m][k    ] = seglimits.map[m][j + 1];
              seglimits.map[m][k + 1] = seglimits.map[m][j + 1];              
            }
            
            // The read limits are not set in stone for the segmented case--it
            // remains to be seen whether this is one of two special cases.
            if (j == 0 && k == 1) {
              int rsize = readKf - readJs;
              int rhalf = ((rsize + GRID_BITS_MASK) / GRID) / 2;
              rhalf *= GRID;
              readKs = readKf - rhalf;
              readJf = readJs + (rsize - rhalf);
              ntmap[m      ]  = ntmap[m +  44] + coreJs - readJs;
              ntmap[m +  11]  = ntmap[m      ];
              ntmap[m +  22]  = ntmap[m +  44] + coreKf - readJs;
              ntmap[m +  33]  = ntmap[m +  22];
              ntmap[m +  66]  = ntmap[m +  77];
              ntmap[m +  77] -= rhalf;
              ntmap[m +  55]  = ntmap[m +  77];
              ntmap[m +  88] = readJs;
              ntmap[m +  99] = readJf;
              ntmap[m + 110] = readJf;
              ntmap[m + 121] = readKf;
            }
            else if (k == j + 1) {
              ntmap[m      ] = ntmap[m +  44] + coreJs - readJs;
              ntmap[m +  11] = ntmap[m +  55] + coreKs - readKs;
              ntmap[m +  22] = ntmap[m +  44] + coreJf - readJs;
              ntmap[m +  33] = ntmap[m +  55] + coreKf - readKs;
              ntmap[m +  66] = ntmap[m +  77];
              ntmap[m +  88] = readJs;
              ntmap[m +  99] = readJf;
              ntmap[m + 110] = readJf;
              ntmap[m + 121] = readKf;
            }
            else {
              ntmap[m      ] = ntmap[m +  44] + coreJs - readJs;
              ntmap[m +  11] = ntmap[m +  55] + coreKs - readKs;
              ntmap[m +  22] = ntmap[m +  44] + coreJf - readJs;
              ntmap[m +  33] = ntmap[m +  55] + coreKf - readKs;
              ntmap[m +  88] = readJs;
              ntmap[m +  99] = readKs;
              ntmap[m + 110] = readJf;
              ntmap[m + 121] = readKf;
            }
          }
        }

        // Check against the list of all GPU-generated Neutral-Territory maps
        bool matched = false;
        int nmatches = 0;
        for (m = 0; m < gpu->pbFrcBlkCounters->_pSysData[4]; m++) {
          bool matchCore = true;
          for (n = 0; n < 44; n++) {
            if (ntmap[n] != gpu->pbNTMaps->_pSysData[NTMAP_SPACE*m + n]) {
              matchCore = false;
            }
          }
          bool matchPos = true;
          for (n = 44; n < 88; n++) {
            if (ntmap[n] != gpu->pbNTMaps->_pSysData[NTMAP_SPACE*m + n]) {
              matchPos = false;
            }
          }
          bool matchReads = true;
          for (n = 88; n < 132; n++) {
            if (ntmap[n] != gpu->pbNTMaps->_pSysData[NTMAP_SPACE*m + n]) {
              matchReads = false;
            }
          }
          if (matchReads && matchPos && matchCore) {
            matched = true;
            ntrAssigned->data[m] += 1;
            nmatches++;
          }
          else if (matchReads || matchPos || matchCore) {
            printf("CheckNTStencils :: Error, incomplete match for %s NT stencil (%2d,%2d) of "
                   "%2d segments for cell [ %2d %2d / %2d %2d ]:\n", subimg, j, k, nsegment,
                   celly, cellz, srcz, srcy);
            for (n = 0; n < 11; n++) {
              printf("                   %6d %6d %6d %6d   %6d %6d %6d %6d   %6d %6d %6d %6d "
                     "-> %6d %6d %6d %6d   %6d %6d %6d %6d   %6d %6d %6d %6d\n", ntmap[n],
                     ntmap[11 + n], ntmap[22 + n], ntmap[33 + n], ntmap[44 + n], ntmap[55 + n],
                     ntmap[66 + n], ntmap[77 + n], ntmap[88 + n], ntmap[99 + n],
                     ntmap[110 + n], ntmap[121 + n],
                     gpu->pbNTMaps->_pSysData[NTMAP_SPACE*m + n      ],
                     gpu->pbNTMaps->_pSysData[NTMAP_SPACE*m + n +  11],
                     gpu->pbNTMaps->_pSysData[NTMAP_SPACE*m + n +  22],
                     gpu->pbNTMaps->_pSysData[NTMAP_SPACE*m + n +  33],
                     gpu->pbNTMaps->_pSysData[NTMAP_SPACE*m + n +  44],
                     gpu->pbNTMaps->_pSysData[NTMAP_SPACE*m + n +  55],
                     gpu->pbNTMaps->_pSysData[NTMAP_SPACE*m + n +  66],
                     gpu->pbNTMaps->_pSysData[NTMAP_SPACE*m + n +  77],
                     gpu->pbNTMaps->_pSysData[NTMAP_SPACE*m + n +  88],
                     gpu->pbNTMaps->_pSysData[NTMAP_SPACE*m + n +  99],
                     gpu->pbNTMaps->_pSysData[NTMAP_SPACE*m + n + 110],
                     gpu->pbNTMaps->_pSysData[NTMAP_SPACE*m + n + 121]);
            }
            printf("CheckNTStencils :: Relevant cell limits are:\n");
            for (n = 0; n < 11; n++) {
              int pcell = ntCollection.map[n][2];
              printf("                   %6d %6d %6d %6d\n", tIMG->cells[pcell].lcapStart,
                     tIMG->cells[pcell].coreStart, tIMG->cells[pcell].coreEnd,
                     tIMG->cells[pcell].hcapEnd);
            }
          }
        }
        if (!matched) {
          printf("CheckNTStencils :: Error, unable to find match for %s NT stencil (%2d,%2d) "
                 "of %2d segments for cell [ %2d %2d / %2d %2d -> %4d ]:\n", subimg, j, k,
                 nsegment, celly, cellz, srcz, srcy, i);
          for (m = 0; m < 11; m++) {
            printf("                   %6d %6d %6d %6d\n", ntmap[88 + m], ntmap[99 + m],
                   ntmap[110 + m], ntmap[121 + m]);
          }
          printf("CheckNTStencils :: Possible relatives, based on cell limits, include:\n");
          for (m = 0; m < gpu->pbFrcBlkCounters->_pSysData[4]; m++) {
            bool inbounds = true;
            for (n = 0; n < 11; n++) {
              int pcell = ntCollection.map[n][2];
              int lbound = (tIMG->cellspace * pcell) + tIMG->cells[pcell].lcapStart;
              int hbound = (tIMG->cellspace * pcell) + tIMG->cells[pcell].hcapEnd;
              if (gpu->pbNTMaps->_pSysData[NTMAP_SPACE*m + n +  88] <  lbound ||
                  gpu->pbNTMaps->_pSysData[NTMAP_SPACE*m + n + 121] >= hbound) {
                inbounds = false;
              }
            }
            if (inbounds) {
              for (n = 0; n < 11; n++) {
                int pcell = ntCollection.map[n][2];
                printf("                   %6d %6d %6d %6d |-> %6d %6d %6d %6d\n",
                       (tIMG->cellspace * pcell) + tIMG->cells[pcell].lcapStart,
                       (tIMG->cellspace * pcell) + tIMG->cells[pcell].coreStart,
                       (tIMG->cellspace * pcell) + tIMG->cells[pcell].coreEnd,
                       (tIMG->cellspace * pcell) + tIMG->cells[pcell].hcapEnd,
                       gpu->pbNTMaps->_pSysData[NTMAP_SPACE*m + n +  88],
                       gpu->pbNTMaps->_pSysData[NTMAP_SPACE*m + n +  99],
                       gpu->pbNTMaps->_pSysData[NTMAP_SPACE*m + n + 110],
                       gpu->pbNTMaps->_pSysData[NTMAP_SPACE*m + n + 121]);
              }
              printf("\n");
            }
          }
        }
        if (nmatches != 1) {
          printf("CheckNTStencils :: Error, %d matches have been found for NT stencil:\n",
                 nmatches);
          for (n = 0; n < 11; n++) {
            printf("                   %6d %6d %6d %6d   %6d %6d %6d %6d   %6d %6d %6d %6d\n",
                   ntmap[n], ntmap[11 + n], ntmap[22 + n], ntmap[33 + n], ntmap[44 + n],
                   ntmap[55 + n], ntmap[66 + n], ntmap[77 + n], ntmap[88 + n], ntmap[99 + n],
                   ntmap[110 + n], ntmap[121 + n]);
          }
        }

        // Log this Neutral Territory map
        ntmap[NTMAP_TYPE_IDX] = (strcmp(subimg, "QQ") == 0) ? 0 : 1;
        ntmap[NTMAP_ATOM_COUNT_IDX] = atomPrefixSum;
        for (m = 0; m < NTMAP_SPACE; m++) {
          tIMG->cells[cellidx].ntmaps[(NTMAP_SPACE *
                                       tIMG->cells[cellidx].nntmap) + m] = ntmap[m];
        }
        tIMG->cells[cellidx].nntmap += 1;
        if (tIMG->cells[cellidx].nntmap >= maxntreg) {
          tIMG->cells[cellidx].ntmaps = (int*)realloc(tIMG->cells[cellidx].ntmaps,
                                                      (tIMG->cells[cellidx].nntmap + 1) *
                                                      NTMAP_SPACE * sizeof(int));
          maxntreg++;
        }
        
        // Update the global counter (eponymous with each cell's counter)
        gblntmap++;
      }
    }
  }
}

//---------------------------------------------------------------------------------------------
// EnumerateNTInteractions: take one of the Neutral-Territory maps and enumerate interactions
//                          within it.
//
// This is a debugging function.
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                numerous gpuBuffer attributes
//   tIMG:        electrostatic or Lennard-Jones sub-image for which to re-derive NT limits
//   ntmap:       the Neutral Territory map, described in comments of kNTX.h.
//---------------------------------------------------------------------------------------------
static int EnumerateNTInteractions(gpuContext gpu, eimage *tIMG, int3 mapid, pairlist *ntPL,
                                   const char* subimg)
{
  int i, j, k, nnb;
  int* latmid;
  PMEFloat4* lcrd;

  // Allocate a local array to store coordinates
  latmid = (int*)malloc(1024 * sizeof(int));
  lcrd = (PMEFloat4*)malloc(1024 * sizeof(PMEFloat4));

  // Set constants based on the type of calculation
  double xcut2;
  if (strcmp(subimg, "QQ") == 0) {
    xcut2 = gpu->sim.es_cutoff + gpu->sim.skinnb;
  }
  else if (strcmp(subimg, "LJ") == 0) {
    xcut2 = gpu->sim.lj_cutoff + gpu->sim.skinnb;
  }
  xcut2 *= xcut2;

  // Get the Neutral-Territory map.  Set the pair list pointer if possible.
  int *ntmap, *plptr;
  i = (tIMG->ydimFull * mapid.y) + mapid.x;
  ntmap = &tIMG->cells[i].ntmaps[NTMAP_SPACE * mapid.z];
  if (ntPL->ngbrRanges != NULL) {
    plptr = &tIMG->cells[i].pairs[tIMG->cells[i].pairlimits[mapid.z]];
  }
  
  // Load the atoms from the expanded array
  for (i = 0; i < 22; i++) {
    for (j = ntmap[88 + i]; j < ntmap[110 + i]; j++) {
      int jbase = ntmap[44 + i] + j - ntmap[88 + i];
      lcrd[jbase] = tIMG->crdq[j];
      latmid[jbase] = tIMG->atoms[j];
    }
  }

  // Loop over all interactions between hash cells: core atoms of cell N to
  // all atoms of cell N + 11, then core atoms of cell N + 11 to all atoms of
  // cell N, where N is within [ 0, 10 ] inclusive.  Only work with the
  // prescribed list of cell-to-cell interactions.
  nnb = 0;
  for (i = 0; i < 31; i++) {

    // Get the limits according to the segment relays array.  This stores combinations
    // of imported hash cell segments that must interact with one another.  They will
    // be stored as hash cell A -> B, C -> D.
    int hcellA = gpu->pbSegmentRelays->_pSysData[i     ];
    int hcellB = gpu->pbSegmentRelays->_pSysData[i + 32];
    int hcellC = gpu->pbSegmentRelays->_pSysData[i + 64];
    int hcellD = gpu->pbSegmentRelays->_pSysData[i + 96];
    
    // Do the core of cell A to all of cell B
    if (hcellA >= 0 && hcellB >= 0) {
      for (j = ntmap[hcellA]; j < ntmap[hcellA + 22]; j++) {
        double atmx = lcrd[j].x;
        double atmy = lcrd[j].y;
        double atmz = lcrd[j].z;
        int atmidj = latmid[j];
        for (k = ntmap[hcellB + 44]; k < ntmap[hcellB + 66]; k++) {

          // Skip pairs that would double-count
          if (hcellB == hcellA + 11 && j <= k) {
            continue;
          }

          // Evaluate the pair distance
          double dx = (double)(lcrd[k].x) - atmx;
          double dy = (double)(lcrd[k].y) - atmy;
          double dz = (double)(lcrd[k].z) - atmz;
          double r2 = dx*dx + dy*dy + dz*dz;
          if (r2 < xcut2) {
            int atmidk = latmid[k];
            if (ntPL->ngbrRanges == NULL) {
              if (atmidj < atmidk) {
                ntPL->ngbrCounts[atmidj] += 1;
              }
              else {
                ntPL->ngbrCounts[atmidk] += 1;
              }
            }
            else {
              if (atmidj < atmidk) {
                int plidx = ntPL->ngbrCounts[atmidj];
                ntPL->ngbrList[plidx] = atmidk;
                ntPL->ngbrRanges[plidx] = sqrt(r2);
                plptr[nnb] = (j << 16) | k;
                ntPL->ngbrCounts[atmidj] = plidx + 1;
              }
              else {
                int plidx = ntPL->ngbrCounts[atmidk];
                ntPL->ngbrList[plidx] = atmidj;
                ntPL->ngbrRanges[plidx] = sqrt(r2);
                plptr[nnb] = (k << 16) | j;
                ntPL->ngbrCounts[atmidk] = plidx + 1;
              }
            }
            nnb++;
          }
        }
      }
    }

    // Do the core of cell C to all of cell D
    if (hcellC >= 0 && hcellD >= 0) {
      for (j = ntmap[hcellC]; j < ntmap[hcellC + 22]; j++) {
        double atmx = lcrd[j].x;
        double atmy = lcrd[j].y;
        double atmz = lcrd[j].z;
        int atmidj = latmid[j];
        for (k = ntmap[hcellD + 44]; k < ntmap[hcellD + 66]; k++) {

          // Skip pairs that would double-count
          if (hcellC == hcellD + 11 && j <= k) {
            continue;
          }

          // Evaluate the pair distance
          double dx = (double)(lcrd[k].x) - atmx;
          double dy = (double)(lcrd[k].y) - atmy;
          double dz = (double)(lcrd[k].z) - atmz;
          double r2 = dx*dx + dy*dy + dz*dz;
          if (r2 < xcut2) {
            int atmidk = latmid[k];
            if (ntPL->ngbrRanges == NULL) {
              if (atmidj < atmidk) {
                ntPL->ngbrCounts[atmidj] += 1;
              }
              else {
                ntPL->ngbrCounts[atmidk] += 1;
              }
            }
            else {
              if (atmidj < atmidk) {
                int plidx = ntPL->ngbrCounts[atmidj];
                ntPL->ngbrList[plidx] = atmidk;
                ntPL->ngbrRanges[plidx] = sqrt(r2);
                plptr[nnb] = (j << 16) | k;
                ntPL->ngbrCounts[atmidj] = plidx + 1;                
              }
              else {
                int plidx = ntPL->ngbrCounts[atmidk];
                ntPL->ngbrList[plidx] = atmidj;
                ntPL->ngbrRanges[plidx] = sqrt(r2);
                plptr[nnb] = (k << 16) | j;
                ntPL->ngbrCounts[atmidk] = plidx + 1;
              }
            }
            nnb++;
          }
        }
      }
    }    
  }
  
  // Free allocated memory
  free(latmid);
  free(lcrd);

  return nnb;
}

//---------------------------------------------------------------------------------------------
// EnumerateNTPairList: enumerate all interactions within a sub-image decomposition based on
//                      its Neutral Territory allocations.
//
// This is a debugging function.
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                numerous gpuBuffer attributes
//   tIMG:        electrostatic or Lennard-Jones sub-image for which to re-derive NT limits
//   subimg:      code for the sub-image to allocate (different cutoffs may imply different
//                levels of hash cell padding in Y and Z)
//---------------------------------------------------------------------------------------------
static void EnumerateNTPairList(gpuContext gpu, eimage *tIMG, pairlist *hcmbPL,
                                const char* subimg)
{
  int i, j;
  pairlist ntPL;

  // Set up the neutral territory pair list and get counts for each atom
  ntPL.ngbrCounts = (int*)calloc(gpu->sim.atoms + 1, sizeof(int));
  ntPL.ngbrRanges = NULL;
  for (i = 0; i < tIMG->ydimFull * tIMG->zdimFull; i++) {

    // Allocate space to hold pair list limits
    if (tIMG->cells[i].nntmap > 0) {
      tIMG->cells[i].pairlimits = (int*)calloc(tIMG->cells[i].nntmap + 1, sizeof(int));
      tIMG->cells[i].pairlimits[0] = 0;
    }
    int lnnb = 0;
    for (j = 0; j < tIMG->cells[i].nntmap; j++) {
      int3 mapid;
      mapid.y = i / tIMG->ydimFull;
      mapid.x = i - (mapid.y * tIMG->ydimFull);
      mapid.z = j;
      lnnb += EnumerateNTInteractions(gpu, tIMG, mapid, &ntPL, subimg);
      tIMG->cells[i].pairlimits[j + 1] = lnnb;
    }
    tIMG->cells[i].pairs = (int*)malloc(lnnb * sizeof(int));
  }
  int nnb = 0;
  for (i = 0; i < gpu->sim.atoms; i++) {
    nnb += ntPL.ngbrCounts[i];
  }
  printf("EnumerateNTPairlist :: There are %d interactions.\n", nnb);
  FinalizePairlistAllocation(gpu, &ntPL);
  for (i = 0; i < tIMG->ydimFull * tIMG->zdimFull; i++) {
    for (j = 0; j < tIMG->cells[i].nntmap; j++) {
      int3 mapid;
      mapid.y = i / tIMG->ydimFull;
      mapid.x = i - (mapid.y * tIMG->ydimFull);
      mapid.z = j;
      EnumerateNTInteractions(gpu, tIMG, mapid, &ntPL, subimg);
    }
  }
  for (i = gpu->sim.atoms - 1; i > 0; i--) {
    ntPL.ngbrCounts[i] -= ntPL.ngbrCounts[i] - ntPL.ngbrCounts[i-1];
  }
  ntPL.ngbrCounts[0] = 0;

  // Check the Neutral-Territory pair list
  printf("EnumerateNTPairlist :: Checking against simpler pair list method:\n");
  bool stt = MatchNeighbors(gpu, tIMG, &ntPL, hcmbPL, "Neutral Territory");
  printf("EnumerateNTPairlist :: Check complete.  ");
  if (stt) {
    printf("No errors were detected.\n");
  }
  else {
    printf("\n");
  }
}

//---------------------------------------------------------------------------------------------
// CheckNTAssignments: check that all NT assignments made on the GPU have been fulfilled on
//                     the CPU.
// 
// Arguments:
//   gpu:          overarching data structure containing simulation information, here used for
//                 the attribute pbNTMaps which is assumed to already have been downloaded
//   ntrAssigned:  
//---------------------------------------------------------------------------------------------
static void CheckNTAssignments(gpuContext gpu, imat *ntrAssigned)
{
  int i, j;
  
  // Check to see that all GPU-derived stencils have been accounted for.  If a CPU-derived
  // stencil matched more than one GPU-derived stencil, this will have been reported.
  for (i = 0; i < ntrAssigned->col; i++) {
    if (ntrAssigned->data[i] == 0) {
      printf("CheckNTStencils :: Failed to find a match for NT stencil %5d:\n", i);
      for (j = 0; j < 11; j++) {
        printf("                   %6d %6d %6d %6d   %6d %6d %6d %6d   %6d %6d %6d %6d\n",
               gpu->pbNTMaps->_pSysData[NTMAP_SPACE*i + j      ],
               gpu->pbNTMaps->_pSysData[NTMAP_SPACE*i + j +  11],
               gpu->pbNTMaps->_pSysData[NTMAP_SPACE*i + j +  22],
               gpu->pbNTMaps->_pSysData[NTMAP_SPACE*i + j +  33],
               gpu->pbNTMaps->_pSysData[NTMAP_SPACE*i + j +  44],
               gpu->pbNTMaps->_pSysData[NTMAP_SPACE*i + j +  55],
               gpu->pbNTMaps->_pSysData[NTMAP_SPACE*i + j +  66],
               gpu->pbNTMaps->_pSysData[NTMAP_SPACE*i + j +  77],
               gpu->pbNTMaps->_pSysData[NTMAP_SPACE*i + j +  88],
               gpu->pbNTMaps->_pSysData[NTMAP_SPACE*i + j +  99],
               gpu->pbNTMaps->_pSysData[NTMAP_SPACE*i + j + 110],
               gpu->pbNTMaps->_pSysData[NTMAP_SPACE*i + j + 121]);
      }
    }
  }
}

//---------------------------------------------------------------------------------------------
// LocateNTMapOnCPU: locate a Neutral-Territory map in the CPU list.
//
// Arguments:
//   gpuNTmap:   pointer to the list of integers constituting a Neutral-Territory map from the
//               GPU
//   tIMG:       the CPU-based image of the system, complete with limits, atom coordinates and
//               indexing, as well as Neutral Territory maps owned by each of its hash cells
//---------------------------------------------------------------------------------------------
static int2 LocateNTMapOnCPU(int *gpuNTmap, eimage *tIMG)
{
  int i, j, k;
  int2 result;

  result.x = -1;
  result.y = -1;
  for (i = 0; i < tIMG->ydimFull * tIMG->zdimFull; i++) {
    for (j = 0; j < tIMG->cells[i].nntmap; j++) {
      bool matched = true;
      k = 0;
      while (matched && k < 132) {
        matched = (gpuNTmap[k] == tIMG->cells[i].ntmaps[NTMAP_SPACE*j + k]);
        k++;
      }
      if (matched && k == 132) {
        result.x = i;
        result.y = j;
        return result;
      }
    }
  }

  return result;
}

//---------------------------------------------------------------------------------------------
// LocatePairInNTMap: locate the two sectors of a Neutral-Territory region that evaluate a
//                    given pair.
//
// Arguments:
//   gpu:        overarching data structure containing simulation information, here used for
//               numerous gpuBuffer attributes
//   atm[1,2]:   the atoms of the pair
//   ntmap:      definition of the NT map, 132 critical integers defining 22 segments
//---------------------------------------------------------------------------------------------
static int2 LocatePairInNTMap(gpuContext gpu, int *ntmap, int atm1, int atm2)
{
  int i;
  
  for (i = 0; i < 31; i++) {

    // Get the limits according to the segment relays array.  This stores combinations
    // of imported hash cell segments that must interact with one another.  They will
    // be stored as hash cell A -> B, C -> D.
    int hcellA = gpu->pbSegmentRelays->_pSysData[i     ];
    int hcellB = gpu->pbSegmentRelays->_pSysData[i + 32];
    int hcellC = gpu->pbSegmentRelays->_pSysData[i + 64];
    int hcellD = gpu->pbSegmentRelays->_pSysData[i + 96];

    if ((atm1 >= ntmap[hcellA]      && atm1 < ntmap[hcellA + 22] &&
         atm2 >= ntmap[hcellB + 44] && atm2 < ntmap[hcellB + 66]) ||
        (atm2 >= ntmap[hcellA]      && atm2 < ntmap[hcellA + 22] &&
         atm1 >= ntmap[hcellB + 44] && atm1 < ntmap[hcellB + 66])) {
      printf("[ %2d -> %2d ] ", hcellA, hcellB);
    }
    if ((atm1 >= ntmap[hcellC]      && atm1 < ntmap[hcellC + 22] &&
         atm2 >= ntmap[hcellD + 44] && atm2 < ntmap[hcellD + 66]) ||
        (atm2 >= ntmap[hcellC]      && atm2 < ntmap[hcellC + 22] &&
         atm1 >= ntmap[hcellD + 44] && atm1 < ntmap[hcellD + 66])) {
      printf("[ %2d -> %2d ] ", hcellC, hcellD);
    }
  }
}

//---------------------------------------------------------------------------------------------
// CheckPairForExclusion: check one of the CPU pairs to see whether it is really an exclusion
//                        and therefore properly excluded by the GPU.
//
// Arguments:
//   gpu:       overarching data structure containing simulation information, here used for
//              the Lennard-Jones exclusions arrays
//   tIMG:      the expanded particle image
//   latmid:    absolute atom ID numbers (as seen in the master topology) arranged as they
//              appear in the NT region
//   paircode:  the bit packed integer containing the indices of local data describing a pair
//              of interacting atoms
//---------------------------------------------------------------------------------------------
static bool CheckPairForExclusion(gpuContext gpu, eimage *tIMG, int* latmid, int paircode)
{
  int i, j;
  
  // See whether the pair code is an exclusion
  int atmA = (paircode >> 16);
  int atmB = (paircode & 0xffff);
  int aidx = latmid[atmA];
  int bidx = latmid[atmB];
  for (i = gpu->pbLJxcIndexing->_pSysData[aidx];
       i < gpu->pbLJxcIndexing->_pSysData[aidx + 1]; i++) {
    if (gpu->pbLJxcPartners->_pSysData[i] == bidx) {
      return true;
    }
  }
  
  return false;
}

//---------------------------------------------------------------------------------------------
// EvaluateGPUNTPairs: with the Neutral-Territory pair list created on the CPU and checked
//                     against a simpler method, the last thing is now to compare the CPU and
//                     GPU representations of the pair list.
//
// Arguments:
//   gpu:        overarching data structure containing simulation information, here used for
//               numerous gpuBuffer attributes
//   tIMG:       the CPU-based image of the system, complete with limits, atom coordinates and
//               indexing, as well as Neutral Territory maps owned by each of its hash cells
//---------------------------------------------------------------------------------------------
static void EvaluateGPUNTPairs(gpuContext gpu, eimage *tIMG, const char* subimg)
{
  int i, j, k, m, n, maptype;
  int* latmid;
  PMEFloat4* lcrd;

  // Allocate a local array to store coordinates
  latmid = (int*)malloc(1024 * sizeof(int));
  lcrd = (PMEFloat4*)malloc(1024 * sizeof(PMEFloat4));

  // Download the GPU-created pair list
  gpu->pbHcmbPairList->Download();
  gpu->pbHcmbTileList->Download();

  // Only a subset of the NT maps will be accessed on this pass through the function
  if (strcmp(subimg, "QQ") == 0) {
    maptype = 0;
  }
  else if (strcmp(subimg, "LJ") == 0) {
    maptype = 1;
  }
  
  // Step through each of the Neutral-Territory tiles as found on the GPU.
  bool problemFound = false;
  for (i = 0; i < gpu->pbFrcBlkCounters->_pSysData[4]; i++) {

    // Skip maps that won't be found in this image
    if (gpu->pbNTMaps->_pSysData[i*NTMAP_SPACE + NTMAP_TYPE_IDX] != maptype) {
      continue;
    }
    
    // Get the Neutral-Territory map.  Set the pair list pointer if possible.
    int *ntmap;
    ntmap = &gpu->pbNTMaps->_pSysData[i*NTMAP_SPACE];
  
    // Load the atoms from the expanded array (in case they are needed for debugging)
    for (j = 0; j < 22; j++) {
      for (k = ntmap[88 + j]; k < ntmap[110 + j]; k++) {
        int kbase     = ntmap[44 + j] + k - ntmap[88 + j];
        lcrd[kbase]   = tIMG->crdq[k];
        latmid[kbase] = tIMG->absid[k];
      }
    }
  
    // Lock on to the NT map on the CPU
    int2 mapid;
    mapid = LocateNTMapOnCPU(&gpu->pbNTMaps->_pSysData[i*NTMAP_SPACE], tIMG);

    // Lock on to the NT map's pairs and tiles as enumerated by the GPU
    int gPstart = ntmap[NTMAP_PAIR_START_IDX];
    int gnpair  = ntmap[NTMAP_PAIR_COUNT_IDX];
    int gTstart = ntmap[NTMAP_TILE_START_IDX];
    int gntile  = ntmap[NTMAP_TILE_COUNT_IDX];
    int *gpuPairList, *gpuTileList, *cpuPairList;
    gpuPairList = &gpu->pbHcmbPairList->_pSysData[gPstart];
    gpuTileList = &gpu->pbHcmbTileList->_pSysData[gTstart];

    // Lock on to the NT map's pairs as enumerated by the CPU.  There are no tiles on
    // the CPU, so any pairs that are not found explicitly must be located as parts of
    // some tile.
    int cPstart = tIMG->cells[mapid.x].pairlimits[mapid.y];
    int cnpair  = tIMG->cells[mapid.x].pairlimits[mapid.y + 1] - cPstart;
    cpuPairList = &tIMG->cells[mapid.x].pairs[cPstart];
    imat cpuPairFulfilled, gpuPairBySingleton, gpuPairByTile;
    cpuPairFulfilled   = CreateImat(1, cnpair);
    gpuPairBySingleton = CreateImat(1, cnpair);
    gpuPairByTile      = CreateImat(1, cnpair);

    // Check the CPU pairs--are there duplicates?
    for (j = 0; j < cnpair; j++) {
      int c1atm = cpuPairList[j];
      int c2atm = (c1atm & 0xffff);
      c1atm >>= 16;
      bool found = false;
      for (k = 0; k < gnpair; k++) {
        int g1atm = gpuPairList[k];
        if (g1atm == HCMB_PL_BLANK) {          
          continue;
        }
        int g2atm = (g1atm & 0xffff);
        g1atm >>= 16;
        if ((c1atm == g1atm && c2atm == g2atm) || (c1atm == g2atm && c2atm == g1atm)) {
          found = true;
          cpuPairFulfilled.data[j]   += 1;
          gpuPairBySingleton.data[j] += 1;
        }
      }
    }
    for (j = 0; j < gntile; j++) {
      int tcoreS = gpuTileList[j];
      int tothrS = (tcoreS & 0xffff);
      tcoreS >>= 16;
      for (k = tcoreS; k < tcoreS + 8; k++) {
        for (m = tothrS; m < tothrS + 4; m++) {
          for (n = 0; n < cnpair; n++) {
            int c1atm = cpuPairList[n];
            int c2atm = (c1atm & 0xffff);
            c1atm >>= 16;
            if ((k == c1atm && m == c2atm) || (k == c2atm && m == c1atm)) {
              cpuPairFulfilled.data[n] += 1;
              gpuPairByTile.data[n]    += 1;
            }
          }
        }
      }
    }

    // Check that no enumerated GPU pairs are Lennard-Jones exclusions
    if (strcmp(subimg, "QQ") != 0) {
      for (j = 0; j < gnpair; j++) {
        if (CheckPairForExclusion(gpu, tIMG, latmid, gpuPairList[j])) {
          int atm1 = (gpuPairList[j] & 0xffff);
          int atm2 = (gpuPairList[j] >> 16);
          int aidx = latmid[atm1];
          int bidx = latmid[atm2];
          printf("\nEvaluateGPUNTPairs :: %4d - %4d of NT region %4d (%s) evaluated by the "
                 "GPU is an exclusion between atoms %6d and %6d\n",  atm1, atm2, i, subimg,
                 aidx, bidx);
          problemFound = true;
        }
      }
    }
    
    // Check to see that all pairs are fulfilled by the GPU's singles and tiles, once
    bool pairmiss = false;
    bool doublecount = false;
    for (j = 0; j < cnpair; j++) {
      if (cpuPairFulfilled.data[j] == 0) {

        // Final check: this may not be a miss, just an exclusion
        if (strcmp(subimg, "QQ") == 0 ||
            CheckPairForExclusion(gpu, tIMG, latmid, cpuPairList[j]) == false) {
          pairmiss = true;
        }
      }
      if (cpuPairFulfilled.data[j] > 1) {
        doublecount = true;
      }
    }
    if (pairmiss || doublecount) {
      printf("EvaluateGPUNTPairs :: Error in NT tile %4d (%s):\n", i, subimg);
      problemFound = true;
      for (j = 0; j < 11; j++) {
        printf("                   %6d %6d %6d %6d   %6d %6d %6d %6d   %6d %6d %6d %6d\n",
               ntmap[j      ], ntmap[j +  11], ntmap[j +  22], ntmap[j +  33],
               ntmap[j +  44], ntmap[j +  55], ntmap[j +  66], ntmap[j +  77],
               ntmap[j +  88], ntmap[j +  99], ntmap[j + 110], ntmap[j + 121]);
      }
      for (j = 0; j < cnpair; j++) {
        if (cpuPairFulfilled.data[j] == 0) {

          // Identify the pair of atoms that is missing, then get the range between them,
          // and finally find the pair of segments containing the atoms that somehow was
          // not fully evaluated.
          int atm1 = (cpuPairList[j] & 0xffff);
          int atm2 = (cpuPairList[j] >> 16);
          double dx = lcrd[atm1].x - lcrd[atm2].x;
          double dy = lcrd[atm1].y - lcrd[atm2].y;
          double dz = lcrd[atm1].z - lcrd[atm2].z;
          double r = sqrt(dx*dx + dy*dy + dz*dz);
          printf("EvaluateGPUNTPairs :: %4d - %4d of NT region %4d (%s) is unmarked on the "
                 "GPU.  Range = %9.4f ", atm1, atm2, i, subimg, r);
          LocatePairInNTMap(gpu, ntmap, atm1, atm2);
          printf("\n");
        }
        if (cpuPairFulfilled.data[j] > 1) {
          int atm1 = (cpuPairList[j] & 0xffff);
          int atm2 = (cpuPairList[j] >> 16);
          printf("EvaluateGPUNTPairs :: %4d - %4d is duplicated on the GPU %2d times, %2d by "
                 "singletons and %2d by tiles. ", atm1, atm2, cpuPairFulfilled.data[j],
                 gpuPairBySingleton.data[j], gpuPairByTile.data[j]);
          LocatePairInNTMap(gpu, ntmap, atm1, atm2);
          printf("\n");
        }
      }
      printf("\n");
    }
    else {
      fprintf(stderr, "\rNT region %6d %s (of %6d) passes check.", i, subimg,
              gpu->pbFrcBlkCounters->_pSysData[4]);
      fflush(stderr);
    }

    
    // Free allocated memory
    DestroyImat(&cpuPairFulfilled);
  }
  if (problemFound) {
    printf("\n");
  }
  else {
    fprintf(stderr, "\r");
    fflush(stderr);
    printf("EvaluateGPUNTPairs :: All %s Neutral Territory regions passed checks.\n", subimg);
  }
  
  // Free allocated memory
  free(lcrd);
  free(latmid);
}

//---------------------------------------------------------------------------------------------
// CheckHoneycombExpansion: check the filtering and subsequent instructions by which the atom
//                          coordinates and properties will be taken from the pHcmbXYSP,
//                          pHcmbZSP, and pHcmbChargeSPLJID arrays into the expanded pencil
//                          representation on the device.
//
// This is a debugging function.
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                numerous gpuBuffer attributes
//   printAll:    flag to have raw data printed to the screen
//   checkPairs:  flag to have all pairs checked against independent, laboriously built lists
//---------------------------------------------------------------------------------------------
void CheckHoneycombExpansion(gpuContext gpu, bool printAll, bool checkPairs)
{
  int i, j;

  // Download the 'original' (that is, simply re-ordered) coordinates and atom
  // properties.  Recalculate the cell indices on the CPU.
  gpu->pbHcmb->Download();
  gpu->pbHcmbCellID->Download();
  gpu->pbHcmbChargeSPLJID->Download();
  int* cpu_cellID;
  double* cpu_primaryCrd;
  cpu_cellID = (int*)malloc(gpu->sim.stride * sizeof(int));
  cpu_primaryCrd = (double*)malloc(gpu->sim.stride3 * sizeof(double));
  CpuCalculateCellID(gpu, cpu_cellID, cpu_primaryCrd);

  // Report the grid size
  if (printAll) {
    double nnQQatpen = (double)(gpu->sim.nQQatoms) / ((double)gpu->sim.npencils);
    double nnLJatpen = (double)(gpu->sim.nLJatoms) / ((double)gpu->sim.npencils);
    double nlsectQQ = gpu->sim.ucell[0][0] * (4.0 / nnQQatpen);
    double nlsectLJ = gpu->sim.ucell[0][0] * (4.0 / nnLJatpen);
    printf("Mean lengths of QQ, LJ segments = %10.6f %10.6f (cutoff %10.6f %10.6f)\n",
           nlsectQQ, nlsectLJ, gpu->sim.es_cutoff, gpu->sim.lj_cutoff);
    printf("Count of neighbor pencils: [ %2d ]\n", gpu->sim.nNgbrPencils);
    printf("Hash cell grid dimensions: [ %2d %2d, padding %d %d (QQ) %d %d (LJ) ]\n",
           gpu->sim.nypencils, gpu->sim.nzpencils, gpu->sim.ypadding, gpu->sim.zpadding,
           gpu->sim.ypadding, gpu->sim.zpadding);
  }

  // Download the sub-image limits that dictate how the primary image will expand
  cudaDeviceSynchronize();
  gpu->pbHcmbSubImgCellLimits->Download();
  if (printAll) {
    printf("QQ Primary Image:\n");
    for (j = 0; j < gpu->sim.nzpencils; j++) {
      int k = 0;
      for (i = 0; i < gpu->sim.nypencils; i++) {
        int idx = (j * gpu->sim.nypencils) + i;
        int m1 = (gpu->pbHcmbSubImgCellLimits->_pSysData[idx + gpu->sim.npencils] >= 0) ?
                  gpu->pbHcmbSubImgCellLimits->_pSysData[idx + gpu->sim.npencils] -
                  gpu->pbHcmbSubImgCellLimits->_pSysData[idx] : -1;
        int m2 = (gpu->pbHcmbSubImgCellLimits->_pSysData[idx + (2 * gpu->sim.npencils)] >= 0) ?
                  gpu->pbHcmbSubImgCellLimits->_pSysData[idx + (3 * gpu->sim.npencils)] -
                  gpu->pbHcmbSubImgCellLimits->_pSysData[idx + (2 * gpu->sim.npencils)] : -1;
        printf(" %6d %3d %3d %6d   ", gpu->pbHcmbSubImgCellLimits->_pSysData[idx], m1, m2,
               gpu->pbHcmbSubImgCellLimits->_pSysData[idx + (3 * gpu->sim.npencils)]);
        k++;
        if (k == 10) {
          printf("\n");
          k = 0;
        }
      }
      if (k > 0) {
        printf("\n");
      }
    }
    printf("\n");
  }
  int paddedPencilCount = 4 * (((gpu->sim.npencils + 31) >> 5) << 5);
  if (printAll) {
    printf("LJ Primary Image:\n");
    for (j = 0; j < gpu->sim.nzpencils; j++) {
      int k = 0;
      for (i = 0; i < gpu->sim.nypencils; i++) {
        int idx = paddedPencilCount + (j * gpu->sim.nypencils) + i;
        int m1 = (gpu->pbHcmbSubImgCellLimits->_pSysData[idx + gpu->sim.npencils] >= 0) ?
                  gpu->pbHcmbSubImgCellLimits->_pSysData[idx + gpu->sim.npencils] -
                  gpu->pbHcmbSubImgCellLimits->_pSysData[idx] : -1;
        int m2 = (gpu->pbHcmbSubImgCellLimits->_pSysData[idx + (2 * gpu->sim.npencils)] >= 0) ?
                  gpu->pbHcmbSubImgCellLimits->_pSysData[idx + (3 * gpu->sim.npencils)] -
                  gpu->pbHcmbSubImgCellLimits->_pSysData[idx + (2 * gpu->sim.npencils)] : -1;
        printf(" %6d %3d %3d %6d   ", gpu->pbHcmbSubImgCellLimits->_pSysData[idx], m1, m2,
               gpu->pbHcmbSubImgCellLimits->_pSysData[idx + (3 * gpu->sim.npencils)]);
        k++;
        if (k == 10) {
          printf("\n");
          k = 0;
        }
      }
      if (k > 0) {
        printf("\n");
      }
    }
    printf("\n");
  }

  // Calculate the cell limits (for the primary image) on the CPU
  int* cpu_cellLimits;
  cpu_cellLimits = CpuCalculateCellLimits(gpu, cpu_cellID, cpu_primaryCrd);

  // Download cell populations
  gpu->pbHcmbSubImgPopulations->Download();
  int pny = gpu->sim.nypencils + (2 * gpu->sim.ypadding);
  int pnz = gpu->sim.nzpencils + gpu->sim.zpadding;
  if (printAll) {
    printf("QQ Expanded Image:\n");
    for (j = 0; j < pnz; j++) {
      for (i = 0; i < pny; i++) {
        int idx = j*pny + i;
        printf("%5d ", gpu->pbHcmbSubImgPopulations->_pSysData[idx]);
      }
      printf("\n");
    }
  }
  pny = gpu->sim.nypencils + 2 * gpu->sim.ypadding;
  pnz = gpu->sim.nzpencils + gpu->sim.zpadding;
  int paddedExpCount = ((gpu->sim.nExpandedPencils + 255) >> 8) << 8;
  if (printAll) {
    printf("LJ Expanded Image:\n");
    for (j = 0; j < pnz; j++) {
      for (i = 0; i < pny; i++) {
        int idx = j*pny + i + paddedExpCount;
        printf("%5d ", gpu->pbHcmbSubImgPopulations->_pSysData[idx]);
      }
      printf("\n");
    }
  }

  // Allocate arrays to hold various expanded images, based on GPU-calculated size parameters.
  eimage cpu_QQImage, cpu_LJImage, xgpu_QQImage, xgpu_LJImage;
  cpu_QQImage = AllocateExpandedCpuImage(gpu, "QQ");
  cpu_LJImage = AllocateExpandedCpuImage(gpu, "LJ");
  xgpu_QQImage = AllocateExpandedCpuImage(gpu, "QQ");
  xgpu_LJImage = AllocateExpandedCpuImage(gpu, "LJ");
  DefineActiveHashCellZones(gpu, &cpu_QQImage, cpu_cellLimits, "QQ");
  DefineActiveHashCellZones(gpu, &cpu_LJImage, cpu_cellLimits, "LJ");
  DefineActiveHashCellZones(gpu, &xgpu_QQImage, cpu_cellLimits, "QQ");
  DefineActiveHashCellZones(gpu, &xgpu_LJImage, cpu_cellLimits, "LJ");
  
  // Make a table of every cell's contents in the expanded image on the CPU.
  gpu->pbHcmbIndex->Download();
  gpu->pbSubImageFilter->Download();
  ExpandCpuImage(gpu, &cpu_QQImage, cpu_cellLimits, cpu_primaryCrd, "QQ");
  ExpandCpuImage(gpu, &cpu_LJImage, cpu_cellLimits, cpu_primaryCrd, "LJ");

  // Check the neighbor lists that result from the honeycomb decompositions
  pairlist simpleQQnl, simpleLJnl, hcmbQQnl, hcmbLJnl;
  if (checkPairs) {
    simpleQQnl = BuildSimplePairlist(gpu, cpu_primaryCrd, "QQ");
    simpleLJnl = BuildSimplePairlist(gpu, cpu_primaryCrd, "LJ");
    BuildCorePairlist(gpu, &cpu_QQImage, &simpleQQnl, "QQ", checkPairs);
    BuildCorePairlist(gpu, &cpu_LJImage, &simpleLJnl, "LJ", checkPairs);
  }
  hcmbQQnl = BuildHoneycombPairlist(gpu, &cpu_QQImage, &simpleQQnl, "QQ", checkPairs);
  hcmbLJnl = BuildHoneycombPairlist(gpu, &cpu_LJImage, &simpleLJnl, "LJ", checkPairs);

  // Download the instructions for making the sub-image expansions as calculated
  // on the device.  Replicate those instructions with CPU code and compare.
  gpu->pbHcmbChargeSPLJID->Download();
  gpu->pbQQExpansions->Download();
  gpu->pbLJExpansions->Download();
  gpu->pbQQTranslations->Download();
  gpu->pbLJTranslations->Download();
  gpu->pbExpansionInsrCount->Download();
  ExpandGpuImage(gpu, &cpu_QQImage, &xgpu_QQImage, cpu_cellLimits, cpu_primaryCrd, "QQ");
  ExpandGpuImage(gpu, &cpu_LJImage, &xgpu_LJImage, cpu_cellLimits, cpu_primaryCrd, "LJ");

  // Check the expanded image produced on the GPU in detail
  CheckGpuExpandedImage(gpu, &cpu_QQImage, cpu_cellLimits, "QQ");
  CheckGpuExpandedImage(gpu, &cpu_LJImage, cpu_cellLimits, "LJ");
  
  // Calculate the NT region decomposition
  gpu->pbNTMaps->Download();
  gpu->pbFrcBlkCounters->Download();
  imat ntrAssigned;
  ntrAssigned = CreateImat(1, gpu->pbFrcBlkCounters->_pSysData[4]);  
  CheckNTStencils(gpu, &cpu_QQImage, cpu_cellLimits, &ntrAssigned, "QQ");
  CheckNTStencils(gpu, &cpu_LJImage, cpu_cellLimits, &ntrAssigned, "LJ");
  CheckNTAssignments(gpu, &ntrAssigned);
  DestroyImat(&ntrAssigned);

  // Accumulate pair interactions in each NT map
  EnumerateNTPairList(gpu, &cpu_QQImage, &hcmbQQnl, "QQ");
  EnumerateNTPairList(gpu, &cpu_LJImage, &hcmbLJnl, "LJ");
  EvaluateGPUNTPairs(gpu, &cpu_QQImage, "QQ");
  EvaluateGPUNTPairs(gpu, &cpu_LJImage, "LJ");
  
  // Free allocated memory
  free(cpu_cellID);
  free(cpu_primaryCrd);
  if (checkPairs) {
    DestroyPairlist(&simpleQQnl);
    DestroyPairlist(&simpleLJnl);
  }
  DestroyPairlist(&hcmbQQnl);
  DestroyPairlist(&hcmbLJnl);
  DestroyEImage(&cpu_QQImage);
  DestroyEImage(&cpu_LJImage);
  DestroyEImage(&xgpu_QQImage);
  DestroyEImage(&xgpu_LJImage);
}
