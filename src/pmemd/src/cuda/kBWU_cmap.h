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
  unsigned int rawID = cSim.pBwuCmapID[2*startidx + tgx];
#ifdef LOCAL_ENERGY
  PMEAccumulator ecmap;
  if (rawID == 0xffffffff) {
    ecmap = (PMEAccumulator)0;
  }
  else {
#else
  if (rawID != 0xffffffff) {
#endif
    unsigned int atmI = rawID >> 24;
    unsigned int atmJ = (rawID >> 16) & 0xff;
    unsigned int atmK = (rawID >>  8) & 0xff;
    unsigned int atmL = rawID & 0xff;
    PMEFloat xij     = atmcrdx[atmI] - atmcrdx[atmJ];
    PMEFloat yij     = atmcrdy[atmI] - atmcrdy[atmJ];
    PMEFloat zij     = atmcrdz[atmI] - atmcrdz[atmJ];
    PMEFloat xkj     = atmcrdx[atmK] - atmcrdx[atmJ];
    PMEFloat ykj     = atmcrdy[atmK] - atmcrdy[atmJ];
    PMEFloat zkj     = atmcrdz[atmK] - atmcrdz[atmJ];
    PMEFloat xlk     = atmcrdx[atmL] - atmcrdx[atmK];
    PMEFloat ylk     = atmcrdy[atmL] - atmcrdy[atmK];
    PMEFloat zlk     = atmcrdz[atmL] - atmcrdz[atmK];

    // Calculate phi and quantities required for gradient
    PMEFloat oneOverRKJ = rsqrt(xkj*xkj + ykj*ykj + zkj*zkj);
    PMEFloat uxkj    = xkj * oneOverRKJ;
    PMEFloat uykj    = ykj * oneOverRKJ;
    PMEFloat uzkj    = zkj * oneOverRKJ;
    PMEFloat dotIJKJ = xij*uxkj + yij*uykj + zij*uzkj;
    PMEFloat upxij   = xij - dotIJKJ * uxkj;
    PMEFloat upyij   = yij - dotIJKJ * uykj;
    PMEFloat upzij   = zij - dotIJKJ * uzkj;
    dotIJKJ         *= oneOverRKJ;
    PMEFloat oneOverRUIJ = rsqrt(upxij*upxij + upyij*upyij + upzij*upzij);
    upxij   *= oneOverRUIJ;
    upyij   *= oneOverRUIJ;
    upzij   *= oneOverRUIJ;
    PMEFloat dotLKKJ = xlk*uxkj + ylk*uykj + zlk*uzkj;
    PMEFloat upxlk   = xlk - dotLKKJ * uxkj;
    PMEFloat upylk   = ylk - dotLKKJ * uykj;
    PMEFloat upzlk   = zlk - dotLKKJ * uzkj;
    dotLKKJ         *= oneOverRKJ;
    PMEFloat oneOverRULK = rsqrt(upxlk*upxlk + upylk*upylk + upzlk*upzlk);
    upxlk  *= oneOverRULK;
    upylk  *= oneOverRULK;
    upzlk  *= oneOverRULK;
    PMEFloat dot    = upxij*upxlk + upyij*upylk + upzij*upzlk;
    PMEFloat cosphi = min(max(dot, (PMEFloat)-1.0), (PMEFloat)1.0);
    PMEFloat cx     = upyij*upzlk - upzij*upylk;
    PMEFloat cy     = upzij*upxlk - upxij*upzlk;
    PMEFloat cz     = upxij*upylk - upyij*upxlk;
    dot             = cx*uxkj + cy*uykj + cz*uzkj;
    PMEFloat sinphi = min(max(dot, (PMEFloat)-1.0), (PMEFloat)1.0);
    PMEFloat phi    = acos(cosphi) * (sinphi >= (PMEFloat)0.0 ?
                                      (PMEFloat)1.0 : (PMEFloat)-1.0);
    PMEFloat upxijk = (uzkj*upyij - uykj*upzij) * oneOverRUIJ;
    PMEFloat upyijk = (uxkj*upzij - uzkj*upxij) * oneOverRUIJ;
    PMEFloat upzijk = (uykj*upxij - uxkj*upyij) * oneOverRUIJ;
    PMEFloat upxjkl = (uykj*upzlk - uzkj*upylk) * oneOverRULK;
    PMEFloat upyjkl = (uzkj*upxlk - uxkj*upzlk) * oneOverRULK;
    PMEFloat upzjkl = (uxkj*upylk - uykj*upxlk) * oneOverRULK;

    // Calculate psi
    unsigned int rawID2 = cSim.pBwuCmapID[2*startidx + GRID + tgx];
    unsigned int atmM = rawID2 & 0xff;
    unsigned int cmapType = rawID2 >> 8;
    PMEFloat oneOverRLK  = rsqrt(xlk*xlk + ylk*ylk + zlk*zlk);
    PMEFloat uxlk        = xlk * oneOverRLK;
    PMEFloat uylk        = ylk * oneOverRLK;
    PMEFloat uzlk        = zlk * oneOverRLK;
    PMEFloat dotJKLK     = xkj*uxlk + ykj*uylk + zkj*uzlk;
    PMEFloat upxjk       = -xkj + dotJKLK*uxlk;
    PMEFloat upyjk       = -ykj + dotJKLK*uylk;
    PMEFloat upzjk       = -zkj + dotJKLK*uzlk;
    dotJKLK             *= oneOverRLK;
    PMEFloat oneOverRUJK = rsqrt(upxjk*upxjk + upyjk*upyjk + upzjk*upzjk);
    upxjk               *= oneOverRUJK;
    upyjk               *= oneOverRUJK;
    upzjk               *= oneOverRUJK;
    PMEFloat xml         = atmcrdx[atmM] - atmcrdx[atmL];
    PMEFloat yml         = atmcrdy[atmM] - atmcrdy[atmL];
    PMEFloat zml         = atmcrdz[atmM] - atmcrdz[atmL];
    PMEFloat dotMLLK     = xml*uxlk + yml*uylk + zml*uzlk;
    PMEFloat upxml       = xml - dotMLLK*uxlk;
    PMEFloat upyml       = yml - dotMLLK*uylk;
    PMEFloat upzml       = zml - dotMLLK*uzlk;
    dotMLLK             *= oneOverRLK;
    PMEFloat oneOverRUML = rsqrt(upxml*upxml + upyml*upyml + upzml*upzml);
    upxml               *= oneOverRUML;
    upyml               *= oneOverRUML;
    upzml               *= oneOverRUML;
    dot                  = upxjk*upxml + upyjk*upyml + upzjk*upzml;
    PMEFloat cospsi      = min(max(dot, (PMEFloat)-1.0), (PMEFloat)1.0);
    cx                   = upyjk*upzml - upzjk*upyml;
    cy                   = upzjk*upxml - upxjk*upzml;
    cz                   = upxjk*upyml - upyjk*upxml;
    dot                  = cx*uxlk + cy*uylk + cz*uzlk;
    PMEFloat sinpsi      = min(max(dot, (PMEFloat)-1.0), (PMEFloat)1.0);
    PMEFloat psi         = acos(cospsi) * (sinpsi >= (PMEFloat)0.0 ?
                                           (PMEFloat)1.0 : (PMEFloat)-1.0);

    // Compute the cross terms
    PMEFloat upxjkl1  = (upyjk*uzlk - upzjk*uylk) * oneOverRUJK;
    PMEFloat upyjkl1  = (upzjk*uxlk - upxjk*uzlk) * oneOverRUJK;
    PMEFloat upzjkl1  = (upxjk*uylk - upyjk*uxlk) * oneOverRUJK;
    PMEFloat upxklm   = (uylk*upzml - uzlk*upyml) * oneOverRUML;
    PMEFloat upyklm   = (uzlk*upxml - uxlk*upzml) * oneOverRUML;
    PMEFloat upzklm   = (uxlk*upyml - uylk*upxml) * oneOverRUML;
    phi              += PI;
    psi              += PI;
    int x             = phi * (CMAP_RESOLUTION / ((PMEFloat)2.0 * PI));
    int y             = psi * (CMAP_RESOLUTION / ((PMEFloat)2.0 * PI));
    PMEFloat phifrac  = (phi - x*((PMEFloat)2.0 * PI / CMAP_RESOLUTION)) * CMAP_RESOLUTION /
                        ((PMEFloat)2.0*PI);
    PMEFloat psifrac  = (psi - y*((PMEFloat)2.0 * PI / CMAP_RESOLUTION)) * CMAP_RESOLUTION /
                        ((PMEFloat)2.0*PI);
    PMEFloat4* pSrc   = &(cSim.pCmapEnergy[cmapType + y*cSim.cmapRowStride +
                                           x*cSim.cmapTermStride]);
    PMEFloat4 E00     = *pSrc;
    PMEFloat4 E01     = pSrc[cSim.cmapRowStride];
    PMEFloat4 E10     = pSrc[cSim.cmapTermStride];
    PMEFloat4 E11     = pSrc[cSim.cmapRowStride + cSim.cmapTermStride];
#ifdef CHARMM_ENERGY
    PMEFloat a00      =                   E00.x;
#endif
    PMEFloat a10      =                   E00.y;
    PMEFloat a20      = (-(PMEFloat)3.0 * E00.x) + ( (PMEFloat)3.0 * E10.x) -
                        ( (PMEFloat)2.0 * E00.y) -                   E10.y;
    PMEFloat a30      = ( (PMEFloat)2.0 * E00.x) - ( (PMEFloat)2.0 * E10.x) +
                                          E00.y  +                   E10.y;
    PMEFloat a01      =                   E00.z;
    PMEFloat a11      =                   E00.w;
    PMEFloat a21      = (-(PMEFloat)3.0 * E00.z) + ( (PMEFloat)3.0 * E10.z) -
                        ( (PMEFloat)2.0 * E00.w) -                   E10.w;
    PMEFloat a31      = ( (PMEFloat)2.0 * E00.z) - ( (PMEFloat)2.0 * E10.z) +
                                          E00.w  +                   E10.w;
    PMEFloat a02      = (-(PMEFloat)3.0 * E00.x) + ( (PMEFloat)3.0 * E01.x) -
                        ( (PMEFloat)2.0 * E00.z) -                   E01.z;
    PMEFloat a12      = (-(PMEFloat)3.0 * E00.y) + ( (PMEFloat)3.0 * E01.y) -
                        ( (PMEFloat)2.0 * E00.w) -                   E01.w;
    PMEFloat a22      = ( (PMEFloat)9.0 * (E00.x - E10.x - E01.x + E11.x) ) +
                        ( (PMEFloat)6.0 * E00.y) + ( (PMEFloat)3.0 * E10.y) -
                        ( (PMEFloat)6.0 * E01.y) - ( (PMEFloat)3.0 * E11.y) +
                        ( (PMEFloat)6.0 * E00.z) - ( (PMEFloat)6.0 * E10.z) +
                        ( (PMEFloat)3.0 * E01.z) - ( (PMEFloat)3.0 * E11.z) +
                        ( (PMEFloat)4.0 * E00.w) + ( (PMEFloat)2.0 * E10.w) +
                        ( (PMEFloat)2.0 * E01.w) +                   E11.w;
    PMEFloat a32      = (-(PMEFloat)6.0 * (E00.x - E10.x - E01.x + E11.x) ) +
                        (-(PMEFloat)3.0 * (E00.y + E10.y - E01.y - E11.y) ) +
                        (-(PMEFloat)4.0 * E00.z) + ( (PMEFloat)4.0 * E10.z) -
                                          E01.w  -                   E11.w  +
                        (-(PMEFloat)2.0 * (E00.w + E10.w + E01.z - E11.z) );
    PMEFloat a03      = ( (PMEFloat)2.0 * E00.x) - ( (PMEFloat)2.0 * E01.x) +
                                          E00.z  +                   E01.z;
    PMEFloat a13      = ( (PMEFloat)2.0 * E00.y) - ( (PMEFloat)2.0 * E01.y) +
                                          E00.w  +                   E01.w;
    PMEFloat a23      = (-(PMEFloat)6.0 * (E00.x - E10.x - E01.x + E11.x) ) +
                        (-(PMEFloat)2.0 * (E00.w + E10.y + E01.w - E11.y) ) +
                        (-(PMEFloat)3.0 * (E00.z - E10.z + E01.z - E11.z) ) +
                        (-(PMEFloat)4.0 * E00.y) -                   E10.w +
                        ( (PMEFloat)4.0 * E01.y) -                   E11.w;
    PMEFloat a33     =  ( (PMEFloat)4.0 * (E00.x - E10.x - E01.x + E11.x) ) +
                        ( (PMEFloat)2.0 * (E00.y + E10.y - E01.y - E11.y) ) +
                        ( (PMEFloat)2.0 * (E00.z - E10.z + E01.z - E11.z) ) +
                                           E00.w  +                  E10.w  +
                                           E01.w  +                  E11.w;
    PMEFloat dPhi    = ((PMEFloat)3.0*a33*phifrac + (PMEFloat)2.0*a23)*phifrac + a13;
    dPhi             = ((PMEFloat)3.0*a32*phifrac + (PMEFloat)2.0*a22)*phifrac + a12 +
                       dPhi*psifrac;
    dPhi             = ((PMEFloat)3.0*a31*phifrac + (PMEFloat)2.0*a21)*phifrac + a11 +
                       dPhi*psifrac;
    dPhi             = ((PMEFloat)3.0*a30*phifrac + (PMEFloat)2.0*a20)*phifrac + a10 +
                       dPhi*psifrac;
    PMEFloat dPsi    = ((PMEFloat)3.0*a33*psifrac + (PMEFloat)2.0*a32)*psifrac + a31;
    dPsi             = ((PMEFloat)3.0*a23*psifrac + (PMEFloat)2.0*a22)*psifrac + a21 +
                       dPsi*phifrac;
    dPsi             = ((PMEFloat)3.0*a13*psifrac + (PMEFloat)2.0*a12)*psifrac + a11 +
                       dPsi*phifrac;
    dPsi             = ((PMEFloat)3.0*a03*psifrac + (PMEFloat)2.0*a02)*psifrac + a01 +
                       dPsi*phifrac;
    dPhi            *= rad_to_deg_coeff;
    dPsi            *= rad_to_deg_coeff;
    upxijk          *= dPhi;
    upyijk          *= dPhi;
    upzijk          *= dPhi;
    upxjkl          *= dPhi;
    upyjkl          *= dPhi;
    upzjkl          *= dPhi;
    upxjkl1         *= dPsi;
    upyjkl1         *= dPsi;
    upzjkl1         *= dPsi;
    upxklm          *= dPsi;
    upyklm          *= dPsi;
    upzklm          *= dPsi;

    // Calculate gradients
    PMEDouble vx =  dotIJKJ*upxijk  + dotLKKJ*upxjkl;
    PMEDouble vy =  dotIJKJ*upyijk  + dotLKKJ*upyjkl;
    PMEDouble vz =  dotIJKJ*upzijk  + dotLKKJ*upzjkl;
    PMEDouble wx = -dotJKLK*upxjkl1 + dotMLLK*upxklm;
    PMEDouble wy = -dotJKLK*upyjkl1 + dotMLLK*upyklm;
    PMEDouble wz = -dotJKLK*upzjkl1 + dotMLLK*upzklm;
#ifdef use_DPFP
    PMEAccumulator iupxijk  = llrint(upxijk * FORCESCALE);
    PMEAccumulator iupyijk  = llrint(upyijk * FORCESCALE);
    PMEAccumulator iupzijk  = llrint(upzijk * FORCESCALE);
    PMEAccumulator iupxjkl  = llrint(upxjkl * FORCESCALE);
    PMEAccumulator iupyjkl  = llrint(upyjkl * FORCESCALE);
    PMEAccumulator iupzjkl  = llrint(upzjkl * FORCESCALE);
    PMEAccumulator iupxjkl1 = llrint(upxjkl1 * FORCESCALE);
    PMEAccumulator iupyjkl1 = llrint(upyjkl1 * FORCESCALE);
    PMEAccumulator iupzjkl1 = llrint(upzjkl1 * FORCESCALE);
    PMEAccumulator iupxklm  = llrint(upxklm * FORCESCALE);
    PMEAccumulator iupyklm  = llrint(upyklm * FORCESCALE);
    PMEAccumulator iupzklm  = llrint(upzklm * FORCESCALE);
    PMEAccumulator ivx = llrint(vx * FORCESCALE);
    PMEAccumulator ivy = llrint(vy * FORCESCALE);
    PMEAccumulator ivz = llrint(vz * FORCESCALE);
    PMEAccumulator iwx = llrint(wx * FORCESCALE);
    PMEAccumulator iwy = llrint(wy * FORCESCALE);
    PMEAccumulator iwz = llrint(wz * FORCESCALE);
    atomicAdd((unsigned long long int*)&atmfrcx[atmI], llitoulli(-iupxijk));
    atomicAdd((unsigned long long int*)&atmfrcy[atmI], llitoulli(-iupyijk));
    atomicAdd((unsigned long long int*)&atmfrcz[atmI], llitoulli(-iupzijk));
    atomicAdd((unsigned long long int*)&atmfrcx[atmJ], llitoulli(iupxijk - ivx - iupxjkl1));
    atomicAdd((unsigned long long int*)&atmfrcy[atmJ], llitoulli(iupyijk - ivy - iupyjkl1));
    atomicAdd((unsigned long long int*)&atmfrcz[atmJ], llitoulli(iupzijk - ivz - iupzjkl1));
    atomicAdd((unsigned long long int*)&atmfrcx[atmK], llitoulli(ivx + iupxjkl - iwx +
                                                                 iupxjkl1));
    atomicAdd((unsigned long long int*)&atmfrcy[atmK], llitoulli(ivy + iupyjkl - iwy +
                                                                 iupyjkl1));
    atomicAdd((unsigned long long int*)&atmfrcz[atmK], llitoulli(ivz + iupzjkl - iwz +
                                                                 iupzjkl1));
    atomicAdd((unsigned long long int*)&atmfrcx[atmL], llitoulli(iupxklm + iwx - iupxjkl));
    atomicAdd((unsigned long long int*)&atmfrcy[atmL], llitoulli(iupyklm + iwy - iupyjkl));
    atomicAdd((unsigned long long int*)&atmfrcz[atmL], llitoulli(iupzklm + iwz - iupzjkl));
    atomicAdd((unsigned long long int*)&atmfrcx[atmM], llitoulli(-iupxklm));
    atomicAdd((unsigned long long int*)&atmfrcy[atmM], llitoulli(-iupyklm));
    atomicAdd((unsigned long long int*)&atmfrcz[atmM], llitoulli(-iupzklm));
#else
    int iupxijk  = llrint(upxijk * BSCALEF);
    int iupyijk  = llrint(upyijk * BSCALEF);
    int iupzijk  = llrint(upzijk * BSCALEF);
    int iupxjkl  = llrint(upxjkl * BSCALEF);
    int iupyjkl  = llrint(upyjkl * BSCALEF);
    int iupzjkl  = llrint(upzjkl * BSCALEF);
    int iupxjkl1 = llrint(upxjkl1 * BSCALEF);
    int iupyjkl1 = llrint(upyjkl1 * BSCALEF);
    int iupzjkl1 = llrint(upzjkl1 * BSCALEF);
    int iupxklm  = llrint(upxklm * BSCALEF);
    int iupyklm  = llrint(upyklm * BSCALEF);
    int iupzklm  = llrint(upzklm * BSCALEF);
    int ivx = llrint(vx * BSCALEF);
    int ivy = llrint(vy * BSCALEF);
    int ivz = llrint(vz * BSCALEF);
    int iwx = llrint(wx * BSCALEF);
    int iwy = llrint(wy * BSCALEF);
    int iwz = llrint(wz * BSCALEF);
    atomicAdd(&atmfrcx[atmI], -iupxijk);
    atomicAdd(&atmfrcy[atmI], -iupyijk);
    atomicAdd(&atmfrcz[atmI], -iupzijk);
    atomicAdd(&atmfrcx[atmJ], iupxijk - ivx - iupxjkl1);
    atomicAdd(&atmfrcy[atmJ], iupyijk - ivy - iupyjkl1);
    atomicAdd(&atmfrcz[atmJ], iupzijk - ivz - iupzjkl1);
    atomicAdd(&atmfrcx[atmK], ivx + iupxjkl - iwx + iupxjkl1);
    atomicAdd(&atmfrcy[atmK], ivy + iupyjkl - iwy + iupyjkl1);
    atomicAdd(&atmfrcz[atmK], ivz + iupzjkl - iwz + iupzjkl1);
    atomicAdd(&atmfrcx[atmL], iupxklm + iwx - iupxjkl);
    atomicAdd(&atmfrcy[atmL], iupyklm + iwy - iupyjkl);
    atomicAdd(&atmfrcz[atmL], iupzklm + iwz - iupzjkl);
    atomicAdd(&atmfrcx[atmM], -iupxklm);
    atomicAdd(&atmfrcy[atmM], -iupyklm);
    atomicAdd(&atmfrcz[atmM], -iupzklm);
#endif
#ifdef LOCAL_ENERGY
    PMEFloat a00 = E00.x;
    PMEFloat E = ((a33*psifrac + a32)*psifrac + a31)*psifrac + a30;
    E          = ((a23*psifrac + a22)*psifrac + a21)*psifrac + a20 + phifrac*E;
    E          = ((a13*psifrac + a12)*psifrac + a11)*psifrac + a10 + phifrac*E;
    E          = ((a03*psifrac + a02)*psifrac + a01)*psifrac + a00 + phifrac*E;
    ecmap      = llrint(ENERGYSCALE * E);
  }

  for (unsigned int stride = warpSize >> 1; stride > 0; stride >>=1) {
    ecmap += __SHFL_DOWN(WARP_MASK, ecmap, stride);
  }
  if (tgx == 0) {
    nrgACC[CMAP_EACC_OFFSET + warpIdx] += ecmap;
  }
#else
  }
#endif // LOCAL_ENERGY
}
