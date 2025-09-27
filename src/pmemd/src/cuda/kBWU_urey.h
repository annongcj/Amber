#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// This is included by kBWU.h with options for #define'd LOCAL_AFE and LOCAL_MBAR.  A third
// define, LOCAL_ENERGY, applies at the level of kBWU.h but cascades down into this code as
// well.  This code does not constitute an entire kernel, but rather a section of one.
//
// Variables already defined before including this code:
//   atmcrd[x,y,z]:  coordinates for imported atoms
//   atmfrc[x,y,z]:  force accumulators for imported atoms
//   nrgACC:         accumulator for energy terms
//   startidx:       starting index for this warp's set of GRID (that is, 32) bonded terms
//   tgx:            thread index within the warp
//---------------------------------------------------------------------------------------------
{
  unsigned int rawID = cSim.pBwuUreyID[startidx + tgx];
#ifdef LOCAL_ENERGY
  PMEAccumulator eurey;
  if (rawID == 0xffffffff) {
    eurey = (PMEAccumulator)0;
  }
  else {
#else
  if (rawID != 0xffffffff) {
#endif
    unsigned int atmI = rawID >> 8;
    unsigned int atmJ = rawID & 0xff;
    PMEDouble xij     = atmcrdx[atmI] - atmcrdx[atmJ];
    PMEDouble yij     = atmcrdy[atmI] - atmcrdy[atmJ];
    PMEDouble zij     = atmcrdz[atmI] - atmcrdz[atmJ];
    PMEDouble rij = sqrt(xij*xij + yij*yij + zij*zij);
    PMEDouble2 UBAngle = cSim.pBwuUrey[startidx + tgx];
    double da = rij - UBAngle.y;
    double df = UBAngle.x * da;
    double dfw = (df + df) / rij;
    PMEDouble fx = dfw * xij;
    PMEDouble fy = dfw * yij;
    PMEDouble fz = dfw * zij;
#ifdef use_DPFP
    PMEAccumulator ifx = llrint(fx * FORCESCALE);
    PMEAccumulator ify = llrint(fy * FORCESCALE);
    PMEAccumulator ifz = llrint(fz * FORCESCALE);
    atomicAdd((unsigned long long int*)&atmfrcx[atmJ], llitoulli(ifx));
    atomicAdd((unsigned long long int*)&atmfrcy[atmJ], llitoulli(ify));
    atomicAdd((unsigned long long int*)&atmfrcz[atmJ], llitoulli(ifz));
    atomicAdd((unsigned long long int*)&atmfrcx[atmI], llitoulli(-ifx));
    atomicAdd((unsigned long long int*)&atmfrcy[atmI], llitoulli(-ify));
    atomicAdd((unsigned long long int*)&atmfrcz[atmI], llitoulli(-ifz));
#else
    int ifx = __float2int_rn(fx * BSCALE);
    int ify = __float2int_rn(fy * BSCALE);
    int ifz = __float2int_rn(fz * BSCALE);
    atomicAdd(&atmfrcx[atmJ], ifx);
    atomicAdd(&atmfrcy[atmJ], ify);
    atomicAdd(&atmfrcz[atmJ], ifz);
    atomicAdd(&atmfrcx[atmI], -ifx);
    atomicAdd(&atmfrcy[atmI], -ify);
    atomicAdd(&atmfrcz[atmI], -ifz);
#endif
#ifdef LOCAL_ENERGY
    eurey = llrint(ENERGYSCALE * df * da);
  }
  for (unsigned int stride = warpSize >> 1; stride > 0; stride >>=1) {
    eurey += __SHFL_DOWN(WARP_MASK, eurey, stride);
  }
  if (tgx == 0) {
    nrgACC[UREY_EACC_OFFSET + warpIdx] += eurey;
  }
#else
  }
#endif // LOCAL_ENERGY
}
