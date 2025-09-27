#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------
#ifndef MATRIX_DATA_STRUCTS
#define MATRIX_DATA_STRUCTS

//---------------------------------------------------------------------------------------------
// dmat: data type needed for the linear algebra, a double precision matrix
//---------------------------------------------------------------------------------------------
struct DMatrix {
  int row;
  int col;
  double* data;
  double** map;
};
typedef struct DMatrix dmat;

//---------------------------------------------------------------------------------------------
// imat: data type needed for storing two-dimensional tables of integers, or, umm, matrices
//---------------------------------------------------------------------------------------------
struct IMatrix {
  int row;
  int col;
  int* data;
  int** map;
};
typedef struct IMatrix imat;

#endif
