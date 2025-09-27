#include "copyright.i"
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
 unsigned int tgx = threadIdx.x & GRID_BITS_MASK;
#ifdef REDUCE_TAN
  double length = 0.0;
#elif defined(REDUCE_DOTPROD)
  double dotprod1 = 0.0;
  double dotprod2 = 0.0;
#endif

  for (int idx = index; idx < 3*cSim.nShuttle; idx += blockDim.x * gridDim.x) {
#ifdef REDUCE_TAN
    length += (cSim.pTangents[idx] * cSim.pTangents[idx]);
#elif defined(REDUCE_DOTPROD)
    dotprod1 += (cSim.pDataShuttle[idx] * cSim.pTangents[idx]); //pDataShuttle at this line should include forces, not the coordinates. 
    dotprod2 += (cSim.pSpringForce[idx]);
#endif
  }

#ifdef REDUCE_TAN
  length = warpReduceSum(length);
#elif defined(REDUCE_DOTPROD)
  dotprod1 = warpReduceSum(dotprod1);
  dotprod2 = warpReduceSum(dotprod2);
#endif

  if ( tgx == 0 ) {
#if (__CUDA_ARCH__ < 600) || (__CUDACC_VER_MAJOR__ < 8)
#ifdef REDUCE_TAN
    voidatomicAdd(cSim.pNorm,length);
#elif defined(REDUCE_DOTPROD)
    voidatomicAdd(cSim.pDotProduct1, dotprod1);
    voidatomicAdd(cSim.pDotProduct2, dotprod2);
#endif
#else
#ifdef REDUCE_TAN
    atomicAdd(cSim.pNorm,length);
#elif defined(REDUCE_DOTPROD)
    atomicAdd(cSim.pDotProduct1, dotprod1);
    atomicAdd(cSim.pDotProduct2, dotprod2);
#endif
#endif
  }
}
