#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------
#ifndef MATRIX_FUNCS
#define MATRIX_FUNCS

#include "matrixDS.h"

dmat CreateDmat(int M, int N);

dmat ReallocDmat(dmat *A, int M, int N);

void DestroyDmat(dmat *A);

imat CreateImat(int M, int N);

imat ReallocImat(imat *A, int M, int N);

void DestroyImat(imat *A);

double DotP(double* V1, double* V2, int N);

void SetDVec(double* V, int N, double s);

void SetIVec(int* V, int N, int s);

int IVecExtreme(int* V, int N, bool takemin=false);

int iMaxAbsValue(imat *P);

double DVecExtreme(double* V, int N, bool takemin=false);

double DSum(double* V, int N);

long long int ISum(int* V, int N);

double DAverage(double* V, int n);

double DStDev(double* V, int n);

void AxbQRRxc(dmat A, double* b);

void BackSub(dmat R, double* b);

void CrossP(double* p, double* q, double* cr);

void HessianNorms(double* invU, double* cdepth);

#endif
