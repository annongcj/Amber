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

//---------------------------------------------------------------------------------------------
// CreateDmat: create an M x N double-precision real matrix, initialized to zero.
//---------------------------------------------------------------------------------------------
dmat CreateDmat(int M, int N)
{
  int i;
  dmat A;

  A.row = M;
  A.col = N;
  A.map = (double**)malloc(M*sizeof(double*));
  A.data = (double*)calloc(M*N, sizeof(double));
  for (i = 0; i < M; i++) {
    A.map[i] = &A.data[N*i];
  }

  return A;
}

//---------------------------------------------------------------------------------------------
// ReallocDmat: re-allocate a double-precision real matrix.  Data in the original matrix will
//              remain "in place" in that indices of the new matrix, so long as they existed
//              in the original matrix, will contain the same data as before.  New indices will
//              contain zeros.
//
// Arguments:
//   A:     the original matrix
//   M:     the new number of rows
//   N:     the new number of columns
//---------------------------------------------------------------------------------------------
dmat ReallocDmat(dmat *A, int M, int N)
{
  int i, j, minM, minN;
  dmat Ap;

  Ap = CreateDmat(M, N);
  minM = (M < A->row) ? M : A->row;
  minN = (N < A->col) ? N : A->col;
  for (i = 0; i < minM; i++) {
    for (j = 0; j < minN; j++) {
      Ap.map[i][j] = A->map[i][j];
    }
  }
  DestroyDmat(A);

  return Ap;
}

//---------------------------------------------------------------------------------------------
// DestroyDmat: destroy a double-precision real matrix.
//---------------------------------------------------------------------------------------------
void DestroyDmat(dmat *A)
{
  free(A->data);
  free(A->map);
}

//---------------------------------------------------------------------------------------------
// CreateImat: create an M x N integer matrix, initialized to zero.
//---------------------------------------------------------------------------------------------
imat CreateImat(int M, int N)
{
  int i;
  imat A;

  A.row = M;
  A.col = N;
  A.map = (int**)malloc(M*sizeof(int*));
  A.data = (int*)calloc(M*N, sizeof(int));
  for (i = 0; i < M; i++) {
    A.map[i] = &A.data[N*i];
  }

  return A;
}

//---------------------------------------------------------------------------------------------
// ReallocImat: re-allocate an integer matrix.  Data in the original matrix will remain "in
//              place" in that indices of the new matrix, so long as they existed in the
//              original matrix, will contain the same data as before.  New indices will
//              contain zeros.
//
// Arguments:
//   A:     the original matrix
//   M:     the new number of rows
//   N:     the new number of columns
//---------------------------------------------------------------------------------------------
imat ReallocImat(imat *A, int M, int N)
{
  int i, j, minM, minN;
  imat Ap;

  Ap = CreateImat(M, N);
  minM = (M < A->row) ? M : A->row;
  minN = (N < A->col) ? N : A->col;
  for (i = 0; i < minM; i++) {
    for (j = 0; j < minN; j++) {
      Ap.map[i][j] = A->map[i][j];
    }
  }
  DestroyImat(A);

  return Ap;
}

//---------------------------------------------------------------------------------------------
// DestroyImat: destroy an integer matrix.
//---------------------------------------------------------------------------------------------
void DestroyImat(imat *A)
{
  free(A->data);
  free(A->map);
}

//---------------------------------------------------------------------------------------------
// DotP: Dot product of two double-precision real vectors V1 and V2 of length N.
//---------------------------------------------------------------------------------------------
double DotP(double* V1, double* V2, int N)
{
  int i;
  double dp;

  dp = 0.0;
  for (i = 0; i < N; i++) {
    dp += V1[i]*V2[i];
  }

  return dp;
}

//---------------------------------------------------------------------------------------------
// SetDVec: set all N elements of a double-precision real vector V to s.
//---------------------------------------------------------------------------------------------
void SetDVec(double* V, int N, double s)
{
  int i;

  for (i = 0; i < N; i++) {
    V[i] = s;
  }
}

//---------------------------------------------------------------------------------------------
// SetIVec: set all N elements of an integer vector V to s.
//---------------------------------------------------------------------------------------------
void SetIVec(int* V, int N, int s)
{
  int i;

  for (i = 0; i < N; i++) {
    V[i] = s;
  }
}

//---------------------------------------------------------------------------------------------
// IVecExtreme: take the extreme value of a vector of integers.  By default, the maximum is
//              taken.  If takemin is set to true, the minimum is found instead.
//
// Arguments:
//   V:    the vector to take the extreme value of
//   N:    length of the vector (or, length to consider)
//---------------------------------------------------------------------------------------------
int IVecExtreme(int* V, int N, bool takemin)
{
  int i;

  int extreme = V[0];
  if (takemin) {
    for (i = 1; i < N; i++) {
      if (V[i] < extreme) {
        extreme = V[i];
      }
    }
  }
  else {
    for (i = 1; i < N; i++) {
      if (V[i] > extreme) {
        extreme = V[i];
      }
    }
  }

  return extreme;
}

//---------------------------------------------------------------------------------------------
// iMaxAbsValue: find the maximum value of an integer matrix
//
// Arguments:
//   P:     the matrix to analyze
//---------------------------------------------------------------------------------------------
int iMaxAbsValue(imat *P)
{
  int i, maxd;

  maxd = 0;
  for (i = 0; i < P->col * P->row; i++) {
    if (abs(P->data[i]) > maxd) {
      maxd = abs(P->data[i]);
    }
  }

  return maxd;
}

//---------------------------------------------------------------------------------------------
// DVecExtreme: take the extreme value of a vector of integers.  By default, the maximum is
//              taken.  If takemin is set to true, the minimum is found instead.
//
// Arguments:
//   V:    the vector to take the extreme value of
//   N:    length of the vector (or, length to consider)
//---------------------------------------------------------------------------------------------
double DVecExtreme(double* V, int N, bool takemin)
{
  int i;

  double extreme = V[0];
  if (takemin) {
    for (i = 1; i < N; i++) {
      if (V[i] < extreme) {
        extreme = V[i];
      }
    }
  }
  else {
    for (i = 1; i < N; i++) {
      if (V[i] > extreme) {
        extreme = V[i];
      }
    }
  }

  return extreme;
}

//---------------------------------------------------------------------------------------------
// DSum: take the sum of a double-precision real vector.
//
// Arguments:
//   V:    the vector to sum
//   N:    length of the vector (or, length to consider)
//---------------------------------------------------------------------------------------------
double DSum(double* V, int N)
{
  int i;
  double sval = 0.0;

  for (i = 0; i < N; i++) {
    sval += V[i];
  }

  return sval;
}

//---------------------------------------------------------------------------------------------
// ISum: take the sum of an integer vector.
//
// Arguments:
//   V:    the vector to sum
//   N:    length of the vector (or, length to consider)
//---------------------------------------------------------------------------------------------
long long int ISum(int* V, int N)
{
  int i;
  long long int sval = 0;

  for (i = 0; i < N; i++) {
    sval += V[i];
  }

  return sval;
}

//---------------------------------------------------------------------------------------------
// DAverage: computes the mean value of a double-precision real vector.
//---------------------------------------------------------------------------------------------
double DAverage(double* V, int n)
{
  return DSum(V, n)/n;
}

//---------------------------------------------------------------------------------------------
// DStDev: computes the standard deviation of a double-precision real vector.
//---------------------------------------------------------------------------------------------
double DStDev(double* V, int n)
{
  int i;
  double m, m2;

  m = DAverage(V, n);
  m2 = 0.0;
  for (i = 0; i < n; i++) {
    m2 += (V[i]-m)*(V[i]-m);
  }
  if (m2 < 1.0e-8) {
    return 0.0;
  }

  return sqrt(m2/(n-1));
}

//---------------------------------------------------------------------------------------------
// AxbQRRxc: function for solving linear least-squares problems.  This performs the first stage
//           of the solution by taking a matrix problem Ax = b, where A { m by n, m >= n, and
//           decomposes A by the "modern classical" QR algorithm to recast the problem as
//           Rx = c, where R is the the upper-triangular component of A and the vector c is
//           implicitly computed as (Q*b).  This algorithm may be found in [REF]:
//
//           Trefethen, Lloyd N. and Bau, David III. "Numerical Linear
//           Algebra." pp.73.  Society for Industrial and Applied Mathematics,
//           Philadelphia.  1997.
//
// Arguments:
//   A:        the matrix to decompose and solve Ax = b
//   b:        solution vector in the linear system of equations
//---------------------------------------------------------------------------------------------
void AxbQRRxc(dmat A, double* b)
{
  int i, j, k;
  double tnm_v, tnm_v2, tempval, sign_v;
  double* v;
  double* vprime;
  double* tmp;

  v = (double*)malloc(A.row*sizeof(double));
  vprime = (double*)malloc(A.row*sizeof(double));
  for (k = 0; k < A.col; k++) {

    // Compute the kth column of Q*
    tnm_v2 = 0.0;
    for (i = 0; i < A.row-k; i++) {
      v[i] = A.map[i+k][k];
      tnm_v2 += v[i]*v[i];
    }
    sign_v = (v[0] >= 0.0) ? 1.0 : -1.0;
    tnm_v = sqrt(tnm_v2);
    tnm_v2 -= v[0]*v[0];
    v[0] += sign_v*tnm_v;
    tnm_v = 1.0/sqrt(tnm_v2 + v[0]*v[0]);
    for (i = 0; i < A.row-k; i++) {
      v[i] = v[i]*tnm_v;
    }

    // Update A as R evolves
    for (i = 0; i < A.col-k; i++) {
      vprime[i] = 0.0;
    }
    for (i = 0; i < A.row-k; i++) {
      tmp = &A.map[i+k][k];
      tempval = v[i];
      for (j = 0; j < A.col-k; j++) {
        vprime[j] += tempval*tmp[j];
      }
    }
    for (i = 0; i < A.row-k; i++) {
      tmp = &A.map[i+k][k];
      tempval = 2.0*v[i];
      for (j = 0; j < A.col-k; j++) {
        tmp[j] -= tempval*vprime[j];
      }
    }

    // Update b as Q* evolves
    tmp = &b[k];
    tempval = 2.0*DotP(v, tmp, A.row-k);
    for (i = 0; i < A.row-k; i++) {
      tmp[i] -= tempval*v[i];
    }
  }

  // Free Allocated Memory
  free(v);
  free(vprime);
}

//---------------------------------------------------------------------------------------------
// BackSub: solve the equation Rx = b to complete the work start by AxbQRRxc.
//
// Arguments:
//   R    : an upper triangular matrix of dimension n (this will contain a considerable
//          portion of garbage, the detritus of QR decomposition, below the N x N UT matrix)
//   b    : a vector of dimension n.  Results are returned in this vector.
//---------------------------------------------------------------------------------------------
void BackSub(dmat R, double* b)
{
  int i, j;
  double multval, temp_b, pivot;

  for (i = R.col-1; i > 0; i--) {
    pivot = 1.0/R.map[i][i];
    temp_b = b[i];
    for (j = i-1; j >= 0; j--) {
      multval = R.map[j][i]*pivot;
      b[j] -= multval*temp_b;
    }
    b[i] *= pivot;
  }
  b[0] /= R.map[0][0];
}

//---------------------------------------------------------------------------------------------
// CrossP: function for finding the cross-product cr of vectors p and q.  Note that vectors p
//         and q are assumed to be three-dimensional and only the first three components of
//         these vectors will be considered.
//---------------------------------------------------------------------------------------------
void CrossP(double* p, double* q, double* cr)
{
  cr[0] = p[1]*q[2] - p[2]*q[1];
  cr[1] = p[2]*q[0] - p[0]*q[2];
  cr[2] = p[0]*q[1] - p[1]*q[0];
}

//---------------------------------------------------------------------------------------------
// HessianNorms: this routine computes the distances between box faces given the inverse box
//               transformation matrix.
//
// Arguments:
//   invU:    the inverse box transformation matrix
//   cdepth:  upon return, stores the distances between box faces
//---------------------------------------------------------------------------------------------
void HessianNorms(double* invU, double* cdepth)
{
  int i;
  double mthx, mthy, mthz;
  double xyz0[3], thx[3], thy[3], thz[3], ic1[3], ic2[3], ic3[3];

  xyz0[0] = invU[0] + invU[1] + invU[2];
  xyz0[1] = invU[4] + invU[5];
  xyz0[2] = invU[8];
  for (i = 0; i < 3; i++) {
    ic1[i] = invU[3*i];
    ic2[i] = invU[3*i + 1];
    ic3[i] = invU[3*i + 2];
  }
  CrossP(ic2, ic3, thx);
  CrossP(ic1, ic3, thy);
  CrossP(ic1, ic2, thz);
  mthx = 1.0/sqrt(thx[0]*thx[0] + thx[1]*thx[1] + thx[2]*thx[2]);
  mthy = 1.0/sqrt(thy[0]*thy[0] + thy[1]*thy[1] + thy[2]*thy[2]);
  mthz = 1.0/sqrt(thz[0]*thz[0] + thz[1]*thz[1] + thz[2]*thz[2]);
  for (i = 0; i < 3; i++) {
    thx[i] *= mthx;
    thy[i] *= mthy;
    thz[i] *= mthz;
  }
  cdepth[0] = fabs(thx[0]*xyz0[0] + thx[1]*xyz0[1] + thx[2]*xyz0[2]);
  cdepth[1] = fabs(thy[0]*xyz0[0] + thy[1]*xyz0[1] + thy[2]*xyz0[2]);
  cdepth[2] = fabs(thz[0]*xyz0[0] + thz[1]*xyz0[1] + thz[2]*xyz0[2]);
}
