#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
#define EP_ONEPOINT
  while (pos < cSim.EP11Offset) {
#define EP_TYPE1
    if (pos < cSim.EP11s) {
#include "kLocalToGlobalLoop.h"
    }
#undef EP_TYPE1
    pos += blockDim.x * gridDim.x;
  }
  while (pos < cSim.EP12Offset) {
    pos -= cSim.EP11Offset;
#define EP_TYPE2
    if (pos < cSim.EP12s) {
#include "kLocalToGlobalLoop.h"
    }
#undef EP_TYPE2
    pos += cSim.EP11Offset + blockDim.x*gridDim.x;
  }
#undef EP_ONEPOINT

#define EP_TWOPOINTS
  while (pos < cSim.EP21Offset) {
    pos -= cSim.EP12Offset;
#define EP_TYPE1
    if (pos < cSim.EP21s) {
#include "kLocalToGlobalLoop.h"
    }
#undef EP_TYPE1
    pos += cSim.EP12Offset + blockDim.x*gridDim.x;
  }
  while (pos < cSim.EP22Offset) {
    pos -= cSim.EP21Offset;
#define EP_TYPE2
    if (pos < cSim.EP22s) {
#include "kLocalToGlobalLoop.h"
    }
#undef EP_TYPE2
    pos += cSim.EP21Offset + blockDim.x*gridDim.x;
  }
#undef EP_TWOPOINTS
  
  // Customized virtual sites
  while (pos < cSim.EPCustomOffset) {
    pos -= cSim.EP22Offset;
    if (pos < cSim.EPCustomCount) {
      double d1 = cSim.pExtraPointCustomD1[pos];
#ifdef EP_NEIGHBORLIST
      int  epidx = cSim.pImageEPCustomIndex[pos];
      int4 epfrm = cSim.pImageEPCustomFrame[pos];
#else
      int  epidx = cSim.pExtraPointCustomIndex[pos];
      int4 epfrm = cSim.pExtraPointCustomFrame[pos];
#endif
      int frmtype = (epfrm.x >> 24);
      epfrm.x = (epfrm.x & 0xffffff);
#ifdef EP_NEIGHBORLIST
      double parentx = cSim.pImageX[epfrm.x];
      double parenty = cSim.pImageY[epfrm.x];
      double parentz = cSim.pImageZ[epfrm.x];
      double fratm1x  = cSim.pImageX[epfrm.y];
      double fratm1y  = cSim.pImageY[epfrm.y];
      double fratm1z  = cSim.pImageZ[epfrm.y];
      double px0sg,py0sg,pz0sg,px1sg,py1sg,pz1sg,px2sg,py2sg,pz2sg;
      if(cSim.isgld>0){
        px0sg = cSim.pImageX0sg[epfrm.x];
        py0sg = cSim.pImageY0sg[epfrm.x];
        pz0sg = cSim.pImageZ0sg[epfrm.x];
        px1sg = cSim.pImageX1sg[epfrm.x];
        py1sg = cSim.pImageY1sg[epfrm.x];
        pz1sg = cSim.pImageZ1sg[epfrm.x];
        px2sg = cSim.pImageX2sg[epfrm.x];
        py2sg = cSim.pImageY2sg[epfrm.x];
        pz2sg = cSim.pImageZ2sg[epfrm.x];
      }
#else
      double parentx = cSim.pAtomX[epfrm.x];
      double parenty = cSim.pAtomY[epfrm.x];
      double parentz = cSim.pAtomZ[epfrm.x];
      double fratm1x  = cSim.pAtomX[epfrm.y];
      double fratm1y  = cSim.pAtomY[epfrm.y];
      double fratm1z  = cSim.pAtomZ[epfrm.y];
      double px0sg,py0sg,pz0sg,px1sg,py1sg,pz1sg,px2sg,py2sg,pz2sg;
      if(cSim.isgld>0){
        px0sg = cSim.pX0sg[epfrm.x];
        py0sg = cSim.pY0sg[epfrm.x];
        pz0sg = cSim.pZ0sg[epfrm.x];
        px1sg = cSim.pX1sg[epfrm.x];
        py1sg = cSim.pY1sg[epfrm.x];
        pz1sg = cSim.pZ1sg[epfrm.x];
        px2sg = cSim.pX2sg[epfrm.x];
        py2sg = cSim.pY2sg[epfrm.x];
        pz2sg = cSim.pZ2sg[epfrm.x];
      }
#endif
      if (frmtype == 4) {
#ifdef EP_NEIGHBORLIST
	cSim.pImageX[epidx] = parentx + d1*(fratm1x - parentx);
	cSim.pImageY[epidx] = parenty + d1*(fratm1y - parenty);
	cSim.pImageZ[epidx] = parentz + d1*(fratm1z - parentz);
        if(cSim.isgld>0){
          cSim.pImageX0sg[epidx] = px0sg + d1*(fratm1x - parentx);
          cSim.pImageY0sg[epidx] = py0sg + d1*(fratm1y - parenty);
          cSim.pImageZ0sg[epidx] = pz0sg + d1*(fratm1z - parentz);
          cSim.pImageX1sg[epidx] = px1sg + d1*(fratm1x - parentx);
          cSim.pImageY1sg[epidx] = py1sg + d1*(fratm1y - parenty);
          cSim.pImageZ1sg[epidx] = pz1sg + d1*(fratm1z - parentz);
          cSim.pImageX1sg[epidx] = px1sg + d1*(fratm1x - parentx);
          cSim.pImageY1sg[epidx] = py1sg + d1*(fratm1y - parenty);
          cSim.pImageZ1sg[epidx] = pz1sg + d1*(fratm1z - parentz);
          cSim.pImageX2sg[epidx] = px2sg + d1*(fratm1x - parentx);
          cSim.pImageY2sg[epidx] = py2sg + d1*(fratm1y - parenty);
          cSim.pImageZ2sg[epidx] = pz2sg + d1*(fratm1z - parentz);
        }
#else
	cSim.pAtomX[epidx] = parentx + d1*(fratm1x - parentx);
	cSim.pAtomY[epidx] = parenty + d1*(fratm1y - parenty);
	cSim.pAtomZ[epidx] = parentz + d1*(fratm1z - parentz);
        if(cSim.isgld>0){
          cSim.pX0sg[epidx] = px0sg + d1*(fratm1x - parentx);
          cSim.pY0sg[epidx] = py0sg + d1*(fratm1y - parenty);
          cSim.pZ0sg[epidx] = pz0sg + d1*(fratm1z - parentz);
          cSim.pX1sg[epidx] = px1sg + d1*(fratm1x - parentx);
          cSim.pY1sg[epidx] = py1sg + d1*(fratm1y - parenty);
          cSim.pZ1sg[epidx] = pz1sg + d1*(fratm1z - parentz);
          cSim.pX2sg[epidx] = px2sg + d1*(fratm1x - parentx);
          cSim.pY2sg[epidx] = py2sg + d1*(fratm1y - parenty);
          cSim.pZ2sg[epidx] = pz2sg + d1*(fratm1z - parentz);
        }
#endif
      }
      else if (frmtype == 5) {
#ifdef EP_NEIGHBORLIST
        double fratm2x  = cSim.pImageX[epfrm.z];
        double fratm2y  = cSim.pImageY[epfrm.z];
        double fratm2z  = cSim.pImageZ[epfrm.z];
#else
        double fratm2x  = cSim.pAtomX[epfrm.z];
        double fratm2y  = cSim.pAtomY[epfrm.z];
        double fratm2z  = cSim.pAtomZ[epfrm.z];
#endif
        double d2 = cSim.pExtraPointCustomD2[pos];
#ifdef EP_NEIGHBORLIST
	cSim.pImageX[epidx] = parentx + d1*(fratm1x - parentx) + d2*(fratm2x - parentx);
	cSim.pImageY[epidx] = parenty + d1*(fratm1y - parenty) + d2*(fratm2y - parenty);
	cSim.pImageZ[epidx] = parentz + d1*(fratm1z - parentz) + d2*(fratm2z - parentz);
        if(cSim.isgld>0){
          cSim.pImageX0sg[epidx] = px0sg + d1*(fratm1x - parentx) + d2*(fratm2x - parentx);
          cSim.pImageY0sg[epidx] = py0sg + d1*(fratm1y - parenty) + d2*(fratm2y - parenty);
          cSim.pImageZ0sg[epidx] = pz0sg + d1*(fratm1z - parentz) + d2*(fratm2z - parentz);
          cSim.pImageX1sg[epidx] = px1sg + d1*(fratm1x - parentx) + d2*(fratm2x - parentx);
          cSim.pImageY1sg[epidx] = py1sg + d1*(fratm1y - parenty) + d2*(fratm2y - parenty);
          cSim.pImageZ1sg[epidx] = pz1sg + d1*(fratm1z - parentz) + d2*(fratm2z - parentz);
          cSim.pImageX2sg[epidx] = px2sg + d1*(fratm1x - parentx) + d2*(fratm2x - parentx);
          cSim.pImageY2sg[epidx] = py2sg + d1*(fratm1y - parenty) + d2*(fratm2y - parenty);
          cSim.pImageZ2sg[epidx] = pz2sg + d1*(fratm1z - parentz) + d2*(fratm2z - parentz);
        }
#else
	cSim.pAtomX[epidx] = parentx + d1*(fratm1x - parentx) + d2*(fratm2x - parentx);
	cSim.pAtomY[epidx] = parenty + d1*(fratm1y - parenty) + d2*(fratm2y - parenty);
	cSim.pAtomZ[epidx] = parentz + d1*(fratm1z - parentz) + d2*(fratm2z - parentz);
        if(cSim.isgld>0){
          cSim.pX0sg[epidx] = px0sg + d1*(fratm1x - parentx) + d2*(fratm2x - parentx);
          cSim.pY0sg[epidx] = py0sg + d1*(fratm1y - parenty) + d2*(fratm2y - parenty);
          cSim.pZ0sg[epidx] = pz0sg + d1*(fratm1z - parentz) + d2*(fratm2z - parentz);
          cSim.pX1sg[epidx] = px1sg + d1*(fratm1x - parentx) + d2*(fratm2x - parentx);
          cSim.pY1sg[epidx] = py1sg + d1*(fratm1y - parenty) + d2*(fratm2y - parenty);
          cSim.pZ1sg[epidx] = pz1sg + d1*(fratm1z - parentz) + d2*(fratm2z - parentz);
          cSim.pX2sg[epidx] = px2sg + d1*(fratm1x - parentx) + d2*(fratm2x - parentx);
          cSim.pY2sg[epidx] = py2sg + d1*(fratm1y - parenty) + d2*(fratm2y - parenty);
          cSim.pZ2sg[epidx] = pz2sg + d1*(fratm1z - parentz) + d2*(fratm2z - parentz);
        }
#endif
      }
      else if (frmtype == 6) {
#ifdef EP_NEIGHBORLIST
        double fratm2x  = cSim.pImageX[epfrm.z];
        double fratm2y  = cSim.pImageY[epfrm.z];
        double fratm2z  = cSim.pImageZ[epfrm.z];
#else
        double fratm2x  = cSim.pAtomX[epfrm.z];
        double fratm2y  = cSim.pAtomY[epfrm.z];
        double fratm2z  = cSim.pAtomZ[epfrm.z];
#endif
        double d2 = cSim.pExtraPointCustomD2[pos];
	double dvx = (fratm1x - parentx) + d2*(fratm2x - fratm1x);
	double dvy = (fratm1y - parenty) + d2*(fratm2y - fratm1y);
	double dvz = (fratm1z - parentz) + d2*(fratm2z - fratm1z);
	double d1magdvec = d1 / sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
#ifdef EP_NEIGHBORLIST
        cSim.pImageX[epidx] = parentx + dvx*d1magdvec;
        cSim.pImageY[epidx] = parenty + dvy*d1magdvec;
        cSim.pImageZ[epidx] = parentz + dvz*d1magdvec;
        if(cSim.isgld>0){
          cSim.pImageX0sg[epidx] = px0sg + dvx*d1magdvec;
          cSim.pImageY0sg[epidx] = py0sg + dvy*d1magdvec;
          cSim.pImageZ0sg[epidx] = pz0sg + dvz*d1magdvec;
          cSim.pImageX1sg[epidx] = px1sg + dvx*d1magdvec;
          cSim.pImageY1sg[epidx] = py1sg + dvy*d1magdvec;
          cSim.pImageZ1sg[epidx] = pz1sg + dvz*d1magdvec;
          cSim.pImageX2sg[epidx] = px2sg + dvx*d1magdvec;
          cSim.pImageY2sg[epidx] = py2sg + dvy*d1magdvec;
          cSim.pImageZ2sg[epidx] = pz2sg + dvz*d1magdvec;
        }
#else
        cSim.pAtomX[epidx] = parentx + dvx*d1magdvec;
        cSim.pAtomY[epidx] = parenty + dvy*d1magdvec;
        cSim.pAtomZ[epidx] = parentz + dvz*d1magdvec;
        if(cSim.isgld>0){
          cSim.pX0sg[epidx] = px0sg + dvx*d1magdvec;
          cSim.pY0sg[epidx] = py0sg + dvy*d1magdvec;
          cSim.pZ0sg[epidx] = pz0sg + dvz*d1magdvec;
          cSim.pX1sg[epidx] = px1sg + dvx*d1magdvec;
          cSim.pY1sg[epidx] = py1sg + dvy*d1magdvec;
          cSim.pZ1sg[epidx] = pz1sg + dvz*d1magdvec;
          cSim.pX2sg[epidx] = px2sg + dvx*d1magdvec;
          cSim.pY2sg[epidx] = py2sg + dvy*d1magdvec;
          cSim.pZ2sg[epidx] = pz2sg + dvz*d1magdvec;
        }
#endif
      }
      else if (frmtype == 7) {
#ifdef EP_NEIGHBORLIST
        double fratm2x  = cSim.pImageX[epfrm.z];
        double fratm2y  = cSim.pImageY[epfrm.z];
        double fratm2z  = cSim.pImageZ[epfrm.z];
#else
        double fratm2x  = cSim.pAtomX[epfrm.z];
        double fratm2y  = cSim.pAtomY[epfrm.z];
        double fratm2z  = cSim.pAtomZ[epfrm.z];
#endif
        PMEFloat d2 = cSim.pExtraPointCustomD2[pos];

	// Need to switch to single precision with small differences in
	// local coordinates here, to conserve registers.
	PMEFloat rabx   = (PMEFloat)(fratm1x - parentx);
	PMEFloat raby   = (PMEFloat)(fratm1y - parenty);
	PMEFloat rabz   = (PMEFloat)(fratm1z - parentz);
	PMEFloat rbcx   = (PMEFloat)(fratm2x - fratm1x);
	PMEFloat rbcy   = (PMEFloat)(fratm2y - fratm1y);
	PMEFloat rbcz   = (PMEFloat)(fratm2z - fratm1z);
        PMEFloat invab2 = (PMEFloat)1.0 / (rabx*rabx + raby*raby + rabz*rabz);
	PMEFloat abbcOabab = (rabx*rbcx + raby*rbcy + rabz*rbcz) * invab2;
	PMEFloat rPerpx = rbcx - abbcOabab*rabx;
	PMEFloat rPerpy = rbcy - abbcOabab*raby;
	PMEFloat rPerpz = rbcz - abbcOabab*rabz;
	PMEFloat magrP2 = rPerpx*rPerpx + rPerpy*rPerpy + rPerpz*rPerpz;
#ifdef use_DPFP
	double sind2, cosd2;
	sincos(d2, &sind2, &cosd2);
#else
	float sind2, cosd2;
	sincosf(d2, &sind2, &cosd2);
#endif	
	PMEFloat rabfac = (PMEFloat)d1 * cosd2 * sqrt(invab2);
	PMEFloat rPfac  = (PMEFloat)d1 * sind2 / sqrt(magrP2);	
#ifdef EP_NEIGHBORLIST
        cSim.pImageX[epidx] = parentx + rabfac*rabx + rPfac*rPerpx;
        cSim.pImageY[epidx] = parenty + rabfac*raby + rPfac*rPerpy;
        cSim.pImageZ[epidx] = parentz + rabfac*rabz + rPfac*rPerpz;
        if(cSim.isgld>0){
          cSim.pImageX0sg[epidx] = px0sg + rabfac*rabx + rPfac*rPerpx;
          cSim.pImageY0sg[epidx] = py0sg + rabfac*raby + rPfac*rPerpy;
          cSim.pImageZ0sg[epidx] = pz0sg + rabfac*rabz + rPfac*rPerpz;
          cSim.pImageX1sg[epidx] = px1sg + rabfac*rabx + rPfac*rPerpx;
          cSim.pImageY1sg[epidx] = py1sg + rabfac*raby + rPfac*rPerpy;
          cSim.pImageZ1sg[epidx] = pz1sg + rabfac*rabz + rPfac*rPerpz;
          cSim.pImageX2sg[epidx] = px2sg + rabfac*rabx + rPfac*rPerpx;
          cSim.pImageY2sg[epidx] = py2sg + rabfac*raby + rPfac*rPerpy;
          cSim.pImageZ2sg[epidx] = pz2sg + rabfac*rabz + rPfac*rPerpz;
        }
#else
        cSim.pAtomX[epidx] = parentx + rabfac*rabx + rPfac*rPerpx;
        cSim.pAtomY[epidx] = parenty + rabfac*raby + rPfac*rPerpy;
        cSim.pAtomZ[epidx] = parentz + rabfac*rabz + rPfac*rPerpz;
        if(cSim.isgld>0){
          cSim.pX0sg[epidx] = px0sg + rabfac*rabx + rPfac*rPerpx;
          cSim.pY0sg[epidx] = py0sg + rabfac*raby + rPfac*rPerpy;
          cSim.pZ0sg[epidx] = pz0sg + rabfac*rabz + rPfac*rPerpz;
          cSim.pX1sg[epidx] = px1sg + rabfac*rabx + rPfac*rPerpx;
          cSim.pY1sg[epidx] = py1sg + rabfac*raby + rPfac*rPerpy;
          cSim.pZ1sg[epidx] = pz1sg + rabfac*rabz + rPfac*rPerpz;
          cSim.pX2sg[epidx] = px2sg + rabfac*rabx + rPfac*rPerpx;
          cSim.pY2sg[epidx] = py2sg + rabfac*raby + rPfac*rPerpy;
          cSim.pZ2sg[epidx] = pz2sg + rabfac*rabz + rPfac*rPerpz;
        }
#endif	
      }
      else if (frmtype == 8) {
#ifdef EP_NEIGHBORLIST
        double fratm2x  = cSim.pImageX[epfrm.z];
        double fratm2y  = cSim.pImageY[epfrm.z];
        double fratm2z  = cSim.pImageZ[epfrm.z];
#else
        double fratm2x  = cSim.pAtomX[epfrm.z];
        double fratm2y  = cSim.pAtomY[epfrm.z];
        double fratm2z  = cSim.pAtomZ[epfrm.z];
#endif
        PMEFloat d2 = cSim.pExtraPointCustomD2[pos];
        PMEFloat d3 = cSim.pExtraPointCustomD3[pos];

        // Again, it is best to switch to single precision with small
        // differences in local coordinates here, to conserve registers.
        PMEFloat rabx   = (PMEFloat)(fratm1x - parentx);
        PMEFloat raby   = (PMEFloat)(fratm1y - parenty);
        PMEFloat rabz   = (PMEFloat)(fratm1z - parentz);
        PMEFloat racx   = (PMEFloat)(fratm2x - parentx);
        PMEFloat racy   = (PMEFloat)(fratm2y - parenty);
        PMEFloat racz   = (PMEFloat)(fratm2z - parentz);
	PMEFloat abacx  = raby*racz - rabz*racy;
	PMEFloat abacy  = rabz*racx - rabx*racz;
	PMEFloat abacz  = rabx*racy - raby*racx;
#ifdef EP_NEIGHBORLIST
        cSim.pImageX[epidx] = parentx + (PMEFloat)d1*rabx + d2*racx + d3*abacx;
        cSim.pImageY[epidx] = parenty + (PMEFloat)d1*raby + d2*racy + d3*abacy;
        cSim.pImageZ[epidx] = parentz + (PMEFloat)d1*rabz + d2*racz + d3*abacz;
        if(cSim.isgld>0){
          cSim.pImageX0sg[epidx] = px0sg + (PMEFloat)d1*rabx + d2*racx + d3*abacx;
          cSim.pImageY0sg[epidx] = py0sg + (PMEFloat)d1*raby + d2*racy + d3*abacy;
          cSim.pImageZ0sg[epidx] = pz0sg + (PMEFloat)d1*rabz + d2*racz + d3*abacz;
          cSim.pImageX1sg[epidx] = px1sg + (PMEFloat)d1*rabx + d2*racx + d3*abacx;
          cSim.pImageY1sg[epidx] = py1sg + (PMEFloat)d1*raby + d2*racy + d3*abacy;
          cSim.pImageZ1sg[epidx] = pz1sg + (PMEFloat)d1*rabz + d2*racz + d3*abacz;
          cSim.pImageX2sg[epidx] = px2sg + (PMEFloat)d1*rabx + d2*racx + d3*abacx;
          cSim.pImageY2sg[epidx] = py2sg + (PMEFloat)d1*raby + d2*racy + d3*abacy;
          cSim.pImageZ2sg[epidx] = pz2sg + (PMEFloat)d1*rabz + d2*racz + d3*abacz;
        }
#else
        cSim.pAtomX[epidx] = parentx + (PMEFloat)d1*rabx + d2*racx + d3*abacx;
        cSim.pAtomY[epidx] = parenty + (PMEFloat)d1*raby + d2*racy + d3*abacy;
        cSim.pAtomZ[epidx] = parentz + (PMEFloat)d1*rabz + d2*racz + d3*abacz;
        if(cSim.isgld>0){
          cSim.pX0sg[epidx] = px0sg + (PMEFloat)d1*rabx + d2*racx + d3*abacx;
          cSim.pY0sg[epidx] = py0sg + (PMEFloat)d1*raby + d2*racy + d3*abacy;
          cSim.pZ0sg[epidx] = pz0sg + (PMEFloat)d1*rabz + d2*racz + d3*abacz;
          cSim.pX1sg[epidx] = px1sg + (PMEFloat)d1*rabx + d2*racx + d3*abacx;
          cSim.pY1sg[epidx] = py1sg + (PMEFloat)d1*raby + d2*racy + d3*abacy;
          cSim.pZ1sg[epidx] = pz1sg + (PMEFloat)d1*rabz + d2*racz + d3*abacz;
          cSim.pX2sg[epidx] = px2sg + (PMEFloat)d1*rabx + d2*racx + d3*abacx;
          cSim.pY2sg[epidx] = py2sg + (PMEFloat)d1*raby + d2*racy + d3*abacy;
          cSim.pZ2sg[epidx] = pz2sg + (PMEFloat)d1*rabz + d2*racz + d3*abacz;
        }
#endif	
      }
      else if (frmtype == 9) {
        PMEFloat rabx   = (PMEFloat)(fratm1x - parentx);
        PMEFloat raby   = (PMEFloat)(fratm1y - parenty);
        PMEFloat rabz   = (PMEFloat)(fratm1z - parentz);
        PMEFloat d1invab  = (PMEFloat)d1 * rsqrt(rabx*rabx + raby*raby + rabz*rabz);
#ifdef EP_NEIGHBORLIST
        cSim.pImageX[epidx] = parentx + d1invab*rabx;
        cSim.pImageY[epidx] = parenty + d1invab*raby;
        cSim.pImageZ[epidx] = parentz + d1invab*rabz;
        if(cSim.isgld>0){
          cSim.pImageX0sg[epidx] = px0sg + d1invab*rabx;
          cSim.pImageY0sg[epidx] = py0sg + d1invab*raby;
          cSim.pImageZ0sg[epidx] = pz0sg + d1invab*rabz;
          cSim.pImageX1sg[epidx] = px1sg + d1invab*rabx;
          cSim.pImageY1sg[epidx] = py1sg + d1invab*raby;
          cSim.pImageZ1sg[epidx] = pz1sg + d1invab*rabz;
          cSim.pImageX2sg[epidx] = px2sg + d1invab*rabx;
          cSim.pImageY2sg[epidx] = py2sg + d1invab*raby;
          cSim.pImageZ2sg[epidx] = pz2sg + d1invab*rabz;
        }
#else
        cSim.pAtomX[epidx] = parentx + d1invab*rabx;
        cSim.pAtomY[epidx] = parenty + d1invab*raby;
        cSim.pAtomZ[epidx] = parentz + d1invab*rabz;
        if(cSim.isgld>0){
          cSim.pX0sg[epidx] = px0sg + d1invab*rabx;
          cSim.pY0sg[epidx] = py0sg + d1invab*raby;
          cSim.pZ0sg[epidx] = pz0sg + d1invab*rabz;
          cSim.pX1sg[epidx] = px1sg + d1invab*rabx;
          cSim.pY1sg[epidx] = py1sg + d1invab*raby;
          cSim.pZ1sg[epidx] = pz1sg + d1invab*rabz;
          cSim.pX2sg[epidx] = px2sg + d1invab*rabx;
          cSim.pY2sg[epidx] = py2sg + d1invab*raby;
          cSim.pZ2sg[epidx] = pz2sg + d1invab*rabz;
        }
#endif	
      }
      else if (frmtype == 10) {
#ifdef EP_NEIGHBORLIST
        double fratm2x  = cSim.pImageX[epfrm.z];
        double fratm2y  = cSim.pImageY[epfrm.z];
        double fratm2z  = cSim.pImageZ[epfrm.z];
        double fratm3x  = cSim.pImageX[epfrm.w];
        double fratm3y  = cSim.pImageY[epfrm.w];
        double fratm3z  = cSim.pImageZ[epfrm.w];
#else
        double fratm2x  = cSim.pAtomX[epfrm.z];
        double fratm2y  = cSim.pAtomY[epfrm.z];
        double fratm2z  = cSim.pAtomZ[epfrm.z];
        double fratm3x  = cSim.pAtomX[epfrm.w];
        double fratm3y  = cSim.pAtomY[epfrm.w];
        double fratm3z  = cSim.pAtomZ[epfrm.w];
#endif
        PMEFloat rabx    = (PMEFloat)(fratm1x - parentx);
        PMEFloat raby    = (PMEFloat)(fratm1y - parenty);
        PMEFloat rabz    = (PMEFloat)(fratm1z - parentz);
        PMEFloat racx    = (PMEFloat)(fratm2x - parentx);
        PMEFloat racy    = (PMEFloat)(fratm2y - parenty);
        PMEFloat racz    = (PMEFloat)(fratm2z - parentz);
        PMEFloat radx    = (PMEFloat)(fratm3x - parentx);
        PMEFloat rady    = (PMEFloat)(fratm3y - parenty);
        PMEFloat radz    = (PMEFloat)(fratm3z - parentz);
        PMEFloat rjax    = (PMEFloat)d1*racx - rabx;
        PMEFloat rjay    = (PMEFloat)d1*racy - raby;
        PMEFloat rjaz    = (PMEFloat)d1*racz - rabz;
        PMEFloat d2      = cSim.pExtraPointCustomD2[pos];
        PMEFloat rjbx    = d2*radx - rabx;
        PMEFloat rjby    = d2*rady - raby;
        PMEFloat rjbz    = d2*radz - rabz;
        PMEFloat rmx     = rjay*rjbz - rjaz*rjby;
        PMEFloat rmy     = rjaz*rjbx - rjax*rjbz;
        PMEFloat rmz     = rjax*rjby - rjay*rjbx;
        PMEFloat d3      = cSim.pExtraPointCustomD3[pos];
	PMEFloat magdvec = d3 * rsqrt(rmx*rmx + rmy*rmy + rmz*rmz);
#ifdef EP_NEIGHBORLIST
        cSim.pImageX[epidx] = parentx + magdvec*rmx;
        cSim.pImageY[epidx] = parenty + magdvec*rmy;
        cSim.pImageZ[epidx] = parentz + magdvec*rmz;
#else
        cSim.pAtomX[epidx] = parentx + magdvec*rmx;
        cSim.pAtomY[epidx] = parenty + magdvec*rmy;
        cSim.pAtomZ[epidx] = parentz + magdvec*rmz;
        if(cSim.isgld>0){
          cSim.pX0sg[epidx] = px0sg + magdvec*rmx;
          cSim.pY0sg[epidx] = py0sg + magdvec*rmy;
          cSim.pZ0sg[epidx] = pz0sg + magdvec*rmz;
          cSim.pX1sg[epidx] = px1sg + magdvec*rmx;
          cSim.pY1sg[epidx] = py1sg + magdvec*rmy;
          cSim.pZ1sg[epidx] = pz1sg + magdvec*rmz;
          cSim.pX2sg[epidx] = px2sg + magdvec*rmx;
          cSim.pY2sg[epidx] = py2sg + magdvec*rmy;
          cSim.pZ2sg[epidx] = pz2sg + magdvec*rmz;
        }
#endif		
      }
    }
    pos += cSim.EP22Offset + blockDim.x*gridDim.x;
  }
}
