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
  unsigned int rawID = cSim.pBwuCImpID[startidx + tgx];
#ifdef LOCAL_ENERGY
  PMEAccumulator ecimp;
  if (rawID == 0xffffffff) {
    ecimp = (PMEAccumulator)0;
  }
  else {
#else
  if (rawID != 0xffffffff) {
#endif
    unsigned int atmI = rawID >> 24;
    unsigned int atmJ = (rawID >> 16) & 0xff;
    unsigned int atmK = (rawID >>  8) & 0xff;
    unsigned int atmL = rawID & 0xff;
    PMEDouble2 impDihedral = cSim.pBwuCImp[startidx + tgx];
    PMEDouble xij      = atmcrdx[atmI] - atmcrdx[atmJ];
    PMEDouble yij      = atmcrdy[atmI] - atmcrdy[atmJ];
    PMEDouble zij      = atmcrdz[atmI] - atmcrdz[atmJ];
    PMEDouble xkj      = atmcrdx[atmK] - atmcrdx[atmJ];
    PMEDouble ykj      = atmcrdy[atmK] - atmcrdy[atmJ];
    PMEDouble zkj      = atmcrdz[atmK] - atmcrdz[atmJ];
    PMEDouble xlk      = atmcrdx[atmL] - atmcrdx[atmK];
    PMEDouble ylk      = atmcrdy[atmL] - atmcrdy[atmK];
    PMEDouble zlk      = atmcrdz[atmL] - atmcrdz[atmK];

    // Calculate phi and quantities required for gradient
    PMEDouble oneOverRKJ = rsqrt(xkj*xkj + ykj*ykj + zkj*zkj);
    PMEDouble uxkj       = xkj * oneOverRKJ;
    PMEDouble uykj       = ykj * oneOverRKJ;
    PMEDouble uzkj       = zkj * oneOverRKJ;
    PMEDouble dotIJKJ = xij*uxkj + yij*uykj + zij*uzkj;
    PMEDouble upxij   = xij - dotIJKJ * uxkj;
    PMEDouble upyij   = yij - dotIJKJ * uykj;
    PMEDouble upzij   = zij - dotIJKJ * uzkj;
    dotIJKJ          *= oneOverRKJ;
    PMEDouble oneOverRUIJ  = rsqrt(upxij*upxij + upyij*upyij + upzij*upzij);
    upxij                 *= oneOverRUIJ;
    upyij                 *= oneOverRUIJ;
    upzij                 *= oneOverRUIJ;
    PMEDouble dotLKKJ = xlk * uxkj + ylk * uykj + zlk * uzkj;
    PMEDouble upxlk   = xlk - dotLKKJ * uxkj;
    PMEDouble upylk   = ylk - dotLKKJ * uykj;
    PMEDouble upzlk   = zlk - dotLKKJ * uzkj;
    dotLKKJ          *= oneOverRKJ;
    PMEDouble oneOverRULK  = rsqrt(upxlk * upxlk + upylk * upylk + upzlk * upzlk);
    upxlk                 *= oneOverRULK;
    upylk                 *= oneOverRULK;
    upzlk                 *= oneOverRULK;
    PMEDouble dot    = upxij*upxlk + upyij*upylk + upzij*upzlk;
    PMEDouble cosphi = min(max(dot, (PMEDouble)-1.0), (PMEDouble)1.0);
    PMEDouble cx     = upyij*upzlk - upzij*upylk;
    PMEDouble cy     = upzij*upxlk - upxij*upzlk;
    PMEDouble cz     = upxij*upylk - upyij*upxlk;
    dot              = cx*uxkj + cy*uykj + cz*uzkj;
    PMEDouble sinphi = min(max(dot, (PMEDouble)-1.0), (PMEDouble)1.0);
    PMEDouble phi    = acos(cosphi) *
      (sinphi >= (PMEDouble)0.0 ? (PMEDouble)1.0 : (PMEDouble)-1.0);
    PMEDouble df     = (PMEDouble)-2.0 * impDihedral.x * (phi - impDihedral.y);

    // Calculate gradient
    PMEDouble upxijk       = df * (uzkj*upyij - uykj*upzij) * oneOverRUIJ;
    PMEDouble upyijk       = df * (uxkj*upzij - uzkj*upxij) * oneOverRUIJ;
    PMEDouble upzijk       = df * (uykj*upxij - uxkj*upyij) * oneOverRUIJ;
    PMEDouble upxjkl       = df * (uykj*upzlk - uzkj*upylk) * oneOverRULK;
    PMEDouble upyjkl       = df * (uzkj*upxlk - uxkj*upzlk) * oneOverRULK;
    PMEDouble upzjkl       = df * (uxkj*upylk - uykj*upxlk) * oneOverRULK;
    PMEDouble vx           = dotIJKJ*upxijk + dotLKKJ*upxjkl;
    PMEDouble vy           = dotIJKJ*upyijk + dotLKKJ*upyjkl;
    PMEDouble vz           = dotIJKJ*upzijk + dotLKKJ*upzjkl;
#ifdef use_DPFP
    PMEAccumulator iupxijk = llrint(upxijk * FORCESCALE);
    PMEAccumulator iupyijk = llrint(upyijk * FORCESCALE);
    PMEAccumulator iupzijk = llrint(upzijk * FORCESCALE);
    PMEAccumulator iupxjkl = llrint(upxjkl * FORCESCALE);
    PMEAccumulator iupyjkl = llrint(upyjkl * FORCESCALE);
    PMEAccumulator iupzjkl = llrint(upzjkl * FORCESCALE);
    PMEAccumulator ivx     = llrint(vx * FORCESCALE);
    PMEAccumulator ivy     = llrint(vy * FORCESCALE);
    PMEAccumulator ivz     = llrint(vz * FORCESCALE);
    atomicAdd((unsigned long long int*)&atmfrcx[atmI], llitoulli(iupxijk));
    atomicAdd((unsigned long long int*)&atmfrcy[atmI], llitoulli(iupyijk));
    atomicAdd((unsigned long long int*)&atmfrcz[atmI], llitoulli(iupzijk));
    atomicAdd((unsigned long long int*)&atmfrcx[atmJ], llitoulli(ivx - iupxijk));
    atomicAdd((unsigned long long int*)&atmfrcy[atmJ], llitoulli(ivy - iupyijk));
    atomicAdd((unsigned long long int*)&atmfrcz[atmJ], llitoulli(ivz - iupzijk));
    atomicAdd((unsigned long long int*)&atmfrcx[atmK], llitoulli(-ivx - iupxjkl));
    atomicAdd((unsigned long long int*)&atmfrcy[atmK], llitoulli(-ivy - iupyjkl));
    atomicAdd((unsigned long long int*)&atmfrcz[atmK], llitoulli(-ivz - iupzjkl));
    atomicAdd((unsigned long long int*)&atmfrcx[atmL], llitoulli(iupxjkl));
    atomicAdd((unsigned long long int*)&atmfrcy[atmL], llitoulli(iupyjkl));
    atomicAdd((unsigned long long int*)&atmfrcz[atmL], llitoulli(iupzjkl));
#else
    int iupxijk = __float2int_rn(upxijk * BSCALEF);
    int iupyijk = __float2int_rn(upyijk * BSCALEF);
    int iupzijk = __float2int_rn(upzijk * BSCALEF);
    int iupxjkl = __float2int_rn(upxjkl * BSCALEF);
    int iupyjkl = __float2int_rn(upyjkl * BSCALEF);
    int iupzjkl = __float2int_rn(upzjkl * BSCALEF);
    int ivx     = __float2int_rn(vx * BSCALEF);
    int ivy     = __float2int_rn(vy * BSCALEF);
    int ivz     = __float2int_rn(vz * BSCALEF);
    atomicAdd(&atmfrcx[atmI],  iupxijk);
    atomicAdd(&atmfrcy[atmI],  iupyijk);
    atomicAdd(&atmfrcz[atmI],  iupzijk);
    atomicAdd(&atmfrcx[atmJ],  ivx - iupxijk);
    atomicAdd(&atmfrcy[atmJ],  ivy - iupyijk);
    atomicAdd(&atmfrcz[atmJ],  ivz - iupzijk);
    atomicAdd(&atmfrcx[atmK], -ivx - iupxjkl);
    atomicAdd(&atmfrcy[atmK], -ivy - iupyjkl);
    atomicAdd(&atmfrcz[atmK], -ivz - iupzjkl);
    atomicAdd(&atmfrcx[atmL],  iupxjkl);
    atomicAdd(&atmfrcy[atmL],  iupyjkl);
    atomicAdd(&atmfrcz[atmL],  iupzjkl);
#endif
#ifdef LOCAL_ENERGY
    ecimp = llrint(ENERGYSCALE * impDihedral.x * (phi - impDihedral.y) *
                   (phi - impDihedral.y));
  }

  for (unsigned int stride = warpSize >> 1; stride > 0; stride >>=1) {
    ecimp += __SHFL_DOWN(WARP_MASK, ecimp, stride);
  }
  if (tgx == 0) {
    nrgACC[CIMP_EACC_OFFSET + warpIdx] += ecimp;
  }
#else
  }
#endif // LOCAL_ENERGY
}
