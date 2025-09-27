#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------
{
#ifdef EP_VIRIAL
  __shared__ PMEDouble sV[THREADS_PER_BLOCK];
  PMEDouble v11 = 0.0;
  PMEDouble v22 = 0.0;
  PMEDouble v33 = 0.0;
#endif

  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
#define EP_ONEPOINT
  while (pos < cSim.EP11Offset) {
#define EP_TYPE1
    if (pos < cSim.EP11s) {
#include "kOrientForcesLoop.h"
    }
#undef EP_TYPE1
    pos += blockDim.x * gridDim.x;
  }
  while (pos < cSim.EP12Offset) {
    pos -= cSim.EP11Offset;
#define EP_TYPE2
    if (pos < cSim.EP12s) {
#include "kOrientForcesLoop.h"
    }
#undef EP_TYPE2
    pos += cSim.EP11Offset + blockDim.x * gridDim.x;
  }
#undef EP_ONEPOINT

#define EP_TWOPOINTS
  while (pos < cSim.EP21Offset) {
    pos -= cSim.EP12Offset;
#define EP_TYPE1
    if (pos < cSim.EP21s) {
#include "kOrientForcesLoop.h"
    }
#undef EP_TYPE1
    pos += cSim.EP12Offset + blockDim.x * gridDim.x;
  }
  while (pos < cSim.EP22Offset) {
    pos -= cSim.EP21Offset;
#define EP_TYPE2
    if (pos < cSim.EP22s) {
#include "kOrientForcesLoop.h"
    }
#undef EP_TYPE2
    pos += cSim.EP21Offset + blockDim.x * gridDim.x;
  }
#undef EP_TWOPOINTS

  // Customized virtual sites
  while (pos < cSim.EPCustomOffset) {
    pos -= cSim.EP22Offset;
    if (pos < cSim.EPCustomCount) {
#ifdef EP_NEIGHBORLIST
      int  epidx = cSim.pImageEPCustomIndex[pos];
      int4 epfrm = cSim.pImageEPCustomFrame[pos];
#else
      int  epidx = cSim.pExtraPointCustomIndex[pos];
      int4 epfrm = cSim.pExtraPointCustomFrame[pos];
#endif

      // Get the accumulated force on the extra point.  The non-bonded force
      // accumulator will point to the standard force accumulator, unless
      // there is a virial calculation.
      PMEDouble forceX = cSim.pNBForceXAccumulator[epidx];
      PMEDouble forceY = cSim.pNBForceYAccumulator[epidx];
      PMEDouble forceZ = cSim.pNBForceZAccumulator[epidx];
      cSim.pForceXAccumulator[epidx] = (PMEAccumulator)0;
      cSim.pForceYAccumulator[epidx] = (PMEAccumulator)0;
      cSim.pForceZAccumulator[epidx] = (PMEAccumulator)0;
      cSim.pNBForceXAccumulator[epidx] = (PMEAccumulator)0;
      cSim.pNBForceYAccumulator[epidx] = (PMEAccumulator)0;
      cSim.pNBForceZAccumulator[epidx] = (PMEAccumulator)0;
#ifdef EP_VIRIAL
      cSim.pNBForceXAccumulator[epidx] = (PMEAccumulator)0;
      cSim.pNBForceYAccumulator[epidx] = (PMEAccumulator)0;
      cSim.pNBForceZAccumulator[epidx] = (PMEAccumulator)0;
#endif
      int frmtype = (epfrm.x >> 24);
      epfrm.x = (epfrm.x & 0xffffff);
      double d1 = cSim.pExtraPointCustomD1[pos];
      if (frmtype == 4) {
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[epfrm.x],
		  llitoulli(llrint((1.0 - d1) * forceX)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[epfrm.x],
		  llitoulli(llrint((1.0 - d1) * forceY)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[epfrm.x],
		  llitoulli(llrint((1.0 - d1) * forceZ)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[epfrm.y],
		  llitoulli(llrint(d1 * forceX)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[epfrm.y],
		  llitoulli(llrint(d1 * forceY)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[epfrm.y],
		  llitoulli(llrint(d1 * forceZ)));
      }
      else if (frmtype == 5) {
        double d2 = cSim.pExtraPointCustomD2[pos];
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[epfrm.x],
		  llitoulli(llrint((1.0 - d1 - d2) * forceX)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[epfrm.x],
		  llitoulli(llrint((1.0 - d1 - d2) * forceY)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[epfrm.x],
		  llitoulli(llrint((1.0 - d1 - d2) * forceZ)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[epfrm.y],
		  llitoulli(llrint(d1 * forceX)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[epfrm.y],
		  llitoulli(llrint(d1 * forceY)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[epfrm.y],
		  llitoulli(llrint(d1 * forceZ)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[epfrm.z],
		  llitoulli(llrint(d2 * forceX)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[epfrm.z],
		  llitoulli(llrint(d2 * forceY)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[epfrm.z],
		  llitoulli(llrint(d2 * forceZ)));
      }
      else if (frmtype == 6) {

	// Need to do this in single precision to keep the register count manageable
	int i;
	PMEFloat ap[3], ab[3], raEP[3], dv[3];
	PMEFloat fratm1[3], fratm2[3], parent[3], epfrc[3];
#ifdef EP_NEIGHBORLIST
        parent[0] = cSim.pImageX[epfrm.x];
        parent[1] = cSim.pImageY[epfrm.x];
        parent[2] = cSim.pImageZ[epfrm.x];
        fratm1[0] = cSim.pImageX[epfrm.y];
        fratm1[1] = cSim.pImageY[epfrm.y];
        fratm1[2] = cSim.pImageZ[epfrm.y];
        fratm2[0] = cSim.pImageX[epfrm.z];
        fratm2[1] = cSim.pImageY[epfrm.z];
        fratm2[2] = cSim.pImageZ[epfrm.z];
        raEP[0]   = cSim.pImageX[epidx] - parent[0];
        raEP[1]   = cSim.pImageY[epidx] - parent[1];
        raEP[2]   = cSim.pImageZ[epidx] - parent[2];
#else
        parent[0] = cSim.pAtomX[epfrm.x];
        parent[1] = cSim.pAtomY[epfrm.x];
        parent[2] = cSim.pAtomZ[epfrm.x];
        fratm1[0] = cSim.pAtomX[epfrm.y];
        fratm1[1] = cSim.pAtomY[epfrm.y];
        fratm1[2] = cSim.pAtomZ[epfrm.y];
        fratm2[0] = cSim.pAtomX[epfrm.z];
        fratm2[1] = cSim.pAtomY[epfrm.z];
        fratm2[2] = cSim.pAtomZ[epfrm.z];
        raEP[0]   = cSim.pAtomX[epidx] - parent[0];
        raEP[1]   = cSim.pAtomY[epidx] - parent[1];
        raEP[2]   = cSim.pAtomZ[epidx] - parent[2];
#endif
	PMEFloat d2 = cSim.pExtraPointCustomD2[pos];
	PMEFloat gamma = (PMEFloat)0.0;
	PMEFloat pbfac = (PMEFloat)0.0;
	epfrc[0] = (PMEFloat)forceX;
	epfrc[1] = (PMEFloat)forceY;
	epfrc[2] = (PMEFloat)forceZ;
        for (i = 0; i < 3; i++) {
	  ap[i] = fratm1[i] - parent[i];
	  ab[i] = fratm2[i] - fratm1[i];
	  dv[i] = ap[i] + d2 * ab[i];
	  gamma += dv[i] * dv[i];
	  pbfac += raEP[i] * epfrc[i];
	}
	gamma = (PMEFloat)d1 / sqrt(gamma);
	pbfac /= (raEP[0]*raEP[0] + raEP[1]*raEP[1] + raEP[2]*raEP[2]);
        PMEFloat gmpboldx = gamma * (epfrc[0] - (pbfac * raEP[0]));
        PMEFloat gmpboldy = gamma * (epfrc[1] - (pbfac * raEP[1]));
        PMEFloat gmpboldz = gamma * (epfrc[2] - (pbfac * raEP[2]));
#ifdef EP_VIRIAL
        v11 += raEP[0]*epfrc[0] + ap[0]*gmpboldx;
        v22 += raEP[1]*epfrc[1] + ap[1]*gmpboldy;
        v33 += raEP[2]*epfrc[2] + ap[2]*gmpboldz;
#endif
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[epfrm.x],
                  llitoulli(llrintf(epfrc[0] - gmpboldx)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[epfrm.x],
                  llitoulli(llrintf(epfrc[1] - gmpboldy)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[epfrm.x],
                  llitoulli(llrintf(epfrc[2] - gmpboldz)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[epfrm.y],
                  llitoulli(llrintf(((PMEFloat)1.0 - d2) * gmpboldx)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[epfrm.y],
                  llitoulli(llrintf(((PMEFloat)1.0 - d2) * gmpboldy)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[epfrm.y],
                  llitoulli(llrintf(((PMEFloat)1.0 - d2) * gmpboldz)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[epfrm.z],
                  llitoulli(llrintf(d2 * gmpboldx)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[epfrm.z],
                  llitoulli(llrintf(d2 * gmpboldy)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[epfrm.z],
                  llitoulli(llrintf(d2 * gmpboldz)));
      }
      else if (frmtype == 7) {

	// Need to do this in single precision to keep the register count manageable
#ifdef EP_NEIGHBORLIST
        PMEFloat parentx = cSim.pImageX[epfrm.x];
        PMEFloat parenty = cSim.pImageY[epfrm.x];
        PMEFloat parentz = cSim.pImageZ[epfrm.x];
        PMEFloat fratm1x = cSim.pImageX[epfrm.y];
        PMEFloat fratm1y = cSim.pImageY[epfrm.y];
        PMEFloat fratm1z = cSim.pImageZ[epfrm.y];
        PMEFloat fratm2x = cSim.pImageX[epfrm.z];
        PMEFloat fratm2y = cSim.pImageY[epfrm.z];
        PMEFloat fratm2z = cSim.pImageZ[epfrm.z];
#else
        PMEFloat parentx = cSim.pAtomX[epfrm.x];
        PMEFloat parenty = cSim.pAtomY[epfrm.x];
        PMEFloat parentz = cSim.pAtomZ[epfrm.x];
        PMEFloat fratm1x = cSim.pAtomX[epfrm.y];
        PMEFloat fratm1y = cSim.pAtomY[epfrm.y];
        PMEFloat fratm1z = cSim.pAtomZ[epfrm.y];
        PMEFloat fratm2x = cSim.pAtomX[epfrm.z];
        PMEFloat fratm2y = cSim.pAtomY[epfrm.z];
        PMEFloat fratm2z = cSim.pAtomZ[epfrm.z];
#endif
	PMEFloat d2 = cSim.pExtraPointCustomD2[pos];
	PMEFloat rabx = fratm1x - parentx;
	PMEFloat raby = fratm1y - parenty;
	PMEFloat rabz = fratm1z - parentz;
	PMEFloat rbcx = fratm2x - fratm1x;
	PMEFloat rbcy = fratm2y - fratm1y;
	PMEFloat rbcz = fratm2z - fratm1z;
	PMEFloat magab2  = rabx*rabx + raby*raby + rabz*rabz;
	PMEFloat magdvec = (rabx*rbcx + raby*rbcy + rabz*rbcz) / magab2;
	PMEFloat magab   = sqrt(magab2);
	PMEFloat invab   = (PMEFloat)1.0 / magab;
	PMEFloat invab2  = invab * invab;
	PMEFloat rPerpx = rbcx - magdvec*rabx;
	PMEFloat rPerpy = rbcy - magdvec*raby;
	PMEFloat rPerpz = rbcz - magdvec*rabz;
	PMEFloat magPr2 = rPerpx*rPerpx + rPerpy*rPerpy + rPerpz*rPerpz;
	PMEFloat magPr  = sqrt(magPr2);
        PMEFloat invPr  = (PMEFloat)1.0 / magPr;
	PMEFloat invPr2 = invPr * invPr;
        PMEFloat epfrcx = (PMEFloat)forceX;
        PMEFloat epfrcy = (PMEFloat)forceY;
        PMEFloat epfrcz = (PMEFloat)forceZ;
	PMEFloat f1fac  = (  rabx*epfrcx +   raby*epfrcy +   rabz*epfrcz) * invab2;
        PMEFloat f2fac  = (rPerpx*epfrcx + rPerpy*epfrcy + rPerpz*epfrcz) * invPr2;
#ifdef use_DPFP
	double cosd2, sind2;
	sincos(d2, &sind2, &cosd2);
#else
	float cosd2, sind2;
	sincosf(d2, &sind2, &cosd2);
#endif
        PMEFloat nf1fac = (PMEFloat)d1 * cosd2 * invab;
	PMEFloat nf2fac = (PMEFloat)d1 * sind2 * invPr;
        PMEFloat abbcOabab = (rabx*rbcx + raby*rbcy + rabz*rbcz) * invab2;
	PMEFloat F1 = epfrcx - f1fac*rabx;
	PMEFloat F2 = F1 - f2fac*rPerpx;
	PMEFloat F3 = f1fac * rPerpx;
#ifdef EP_VIRIAL
#  ifdef EP_NEIGHBORLIST
	PMEFloat raEPx = cSim.pImageX[epidx] - (PMEDouble)parentx;
#  else
	PMEFloat raEPx = cSim.pAtomX[epidx] - (PMEDouble)parentx;
#  endif
	PMEFloat racx  = fratm2x - parentx;
        v11 += raEPx*epfrcx - rabx*(nf1fac*F1 + nf2fac*(abbcOabab*F2 + F3)) - racx*nf2fac*F2;
#endif
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[epfrm.x],
                  llitoulli(llrintf(epfrcx - nf1fac*F1 + nf2fac*(abbcOabab*F2 + F3))));
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[epfrm.y],
                  llitoulli(llrintf(nf1fac*F1 - nf2fac*(F2 + abbcOabab*F2 + F3))));
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[epfrm.z],
                  llitoulli(llrintf(nf2fac * F2)));
	F1 = epfrcy - f1fac*raby;
	F2 = F1 - f2fac*rPerpy;
	F3 = f1fac * rPerpy;
#ifdef EP_VIRIAL
#  ifdef EP_NEIGHBORLIST
	PMEFloat raEPy = cSim.pImageY[epidx] - (PMEDouble)parenty;
#  else
	PMEFloat raEPy = cSim.pAtomY[epidx] - (PMEDouble)parenty;
#  endif
	PMEFloat racy  = fratm2y - parenty;
        v22 += raEPy*epfrcy - raby*(nf1fac*F1 + nf2fac*(abbcOabab*F2 + F3)) - racy*nf2fac*F2;
#endif
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[epfrm.x],
                  llitoulli(llrintf(epfrcy - nf1fac*F1 + nf2fac*(abbcOabab*F2 + F3))));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[epfrm.y],
                  llitoulli(llrintf(nf1fac*F1 - nf2fac*(F2 + abbcOabab*F2 + F3))));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[epfrm.z],
                  llitoulli(llrintf(nf2fac * F2)));
	F1 = epfrcz - f1fac*rabz;
	F2 = F1 - f2fac*rPerpz;
	F3 = f1fac * rPerpz;
#ifdef EP_VIRIAL
#  ifdef EP_NEIGHBORLIST
	PMEFloat raEPz = cSim.pImageZ[epidx] - (PMEDouble)parentz;
#  else
	PMEFloat raEPz = cSim.pAtomZ[epidx] - (PMEDouble)parentz;
#  endif
	PMEFloat racz  = fratm2z - parentz;
        v33 += raEPz*epfrcz - rabz*(nf1fac*F1 + nf2fac*(abbcOabab*F2 + F3)) - racz*nf2fac*F2;
#endif
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[epfrm.x],
                  llitoulli(llrintf(epfrcz - nf1fac*F1 + nf2fac*(abbcOabab*F2 + F3))));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[epfrm.y],
                  llitoulli(llrintf(nf1fac*F1 - nf2fac*(F2 + abbcOabab*F2 + F3))));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[epfrm.z],
                  llitoulli(llrintf(nf2fac * F2)));
      }
      else if (frmtype == 8) {
#ifdef EP_NEIGHBORLIST
        PMEFloat parentx = cSim.pImageX[epfrm.x];
        PMEFloat parenty = cSim.pImageY[epfrm.x];
        PMEFloat parentz = cSim.pImageZ[epfrm.x];
        PMEFloat fratm1x = cSim.pImageX[epfrm.y];
        PMEFloat fratm1y = cSim.pImageY[epfrm.y];
        PMEFloat fratm1z = cSim.pImageZ[epfrm.y];
        PMEFloat fratm2x = cSim.pImageX[epfrm.z];
        PMEFloat fratm2y = cSim.pImageY[epfrm.z];
        PMEFloat fratm2z = cSim.pImageZ[epfrm.z];
#  ifdef EP_VIRIAL
        PMEFloat raEPx   = cSim.pImageX[epidx] - parentx;
        PMEFloat raEPy   = cSim.pImageY[epidx] - parenty;
        PMEFloat raEPz   = cSim.pImageZ[epidx] - parentz;
#  endif
#else
        PMEFloat parentx = cSim.pAtomX[epfrm.x];
        PMEFloat parenty = cSim.pAtomY[epfrm.x];
        PMEFloat parentz = cSim.pAtomZ[epfrm.x];
        PMEFloat fratm1x = cSim.pAtomX[epfrm.y];
        PMEFloat fratm1y = cSim.pAtomY[epfrm.y];
        PMEFloat fratm1z = cSim.pAtomZ[epfrm.y];
        PMEFloat fratm2x = cSim.pAtomX[epfrm.z];
        PMEFloat fratm2y = cSim.pAtomY[epfrm.z];
        PMEFloat fratm2z = cSim.pAtomZ[epfrm.z];
#  ifdef EP_VIRIAL
        PMEFloat raEPx   = cSim.pAtomX[epidx] - parentx;
        PMEFloat raEPy   = cSim.pAtomY[epidx] - parenty;
        PMEFloat raEPz   = cSim.pAtomZ[epidx] - parentz;
#  endif
#endif
        PMEFloat epfrcx = (PMEFloat)forceX;
        PMEFloat epfrcy = (PMEFloat)forceY;
        PMEFloat epfrcz = (PMEFloat)forceZ;
        PMEFloat rabx   = (PMEFloat)(fratm1x - parentx);
        PMEFloat raby   = (PMEFloat)(fratm1y - parenty);
        PMEFloat rabz   = (PMEFloat)(fratm1z - parentz);
        PMEFloat racx   = (PMEFloat)(fratm2x - parentx);
        PMEFloat racy   = (PMEFloat)(fratm2y - parenty);
        PMEFloat racz   = (PMEFloat)(fratm2z - parentz);
	PMEFloat d2     = cSim.pExtraPointCustomD2[pos];
	PMEFloat d3     = cSim.pExtraPointCustomD3[pos];
        PMEFloat g01    = d3 * racz;
        PMEFloat g02    = d3 * racy;
        PMEFloat g12    = d3 * racx;
	PMEFloat fbx    =  (PMEFloat)d1*epfrcx - g01*epfrcy + g02*epfrcz;
	PMEFloat fby    =  g01*epfrcx + (PMEFloat)d1*epfrcy - g12*epfrcz;
	PMEFloat fbz    = -g02*epfrcx + g12*epfrcy + (PMEFloat)d1*epfrcz;
        g01             = d3 * rabz;
        g02             = d3 * raby;
        g12             = d3 * rabx;
	PMEFloat fcx    =   d2*epfrcx + g01*epfrcy - g02*epfrcz;
	PMEFloat fcy    = -g01*epfrcx +  d2*epfrcy + g12*epfrcz;
	PMEFloat fcz    =  g02*epfrcx - g12*epfrcy +  d2*epfrcz;
#ifdef EP_VIRIAL
        v11 += raEPx*epfrcx - rabx*fbx - racx*fcx;
        v22 += raEPy*epfrcy - raby*fby - racy*fcy;
        v33 += raEPz*epfrcz - rabz*fbz - racz*fcz;
#endif
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[epfrm.x],
		  llitoulli(llrintf(epfrcx - fbx - fcx)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[epfrm.x],
		  llitoulli(llrintf(epfrcy - fby - fcy)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[epfrm.x],
		  llitoulli(llrintf(epfrcz - fbz - fcz)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[epfrm.y],
		  llitoulli(llrintf(fbx)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[epfrm.y],
		  llitoulli(llrintf(fby)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[epfrm.y],
		  llitoulli(llrintf(fbz)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[epfrm.z],
		  llitoulli(llrintf(fcx)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[epfrm.z],
		  llitoulli(llrintf(fcy)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[epfrm.z],
		  llitoulli(llrintf(fcz)));
      }
      else if (frmtype == 9) {
#ifdef EP_NEIGHBORLIST
        PMEFloat parentx = cSim.pImageX[epfrm.x];
        PMEFloat parenty = cSim.pImageY[epfrm.x];
        PMEFloat parentz = cSim.pImageZ[epfrm.x];
        PMEFloat fratm1x = cSim.pImageX[epfrm.y];
        PMEFloat fratm1y = cSim.pImageY[epfrm.y];
        PMEFloat fratm1z = cSim.pImageZ[epfrm.y];
#else
        PMEFloat parentx = cSim.pAtomX[epfrm.x];
        PMEFloat parenty = cSim.pAtomY[epfrm.x];
        PMEFloat parentz = cSim.pAtomZ[epfrm.x];
        PMEFloat fratm1x = cSim.pAtomX[epfrm.y];
        PMEFloat fratm1y = cSim.pAtomY[epfrm.y];
        PMEFloat fratm1z = cSim.pAtomZ[epfrm.y];
#endif
        PMEFloat epfrcx  = (PMEFloat)forceX;
        PMEFloat epfrcy  = (PMEFloat)forceY;
        PMEFloat epfrcz  = (PMEFloat)forceZ;
        PMEFloat rabx    = (PMEFloat)(fratm1x - parentx);
        PMEFloat raby    = (PMEFloat)(fratm1y - parenty);
        PMEFloat rabz    = (PMEFloat)(fratm1z - parentz);
	PMEFloat invab   = rsqrt(rabx*rabx + raby*raby + rabz*rabz);
	PMEFloat fproj   = (rabx*epfrcx + raby*epfrcy + rabz*epfrcz) * invab * invab;
	invab           *= (PMEFloat)d1;
        PMEFloat fbx     = invab * (epfrcx - fproj*rabx);
        PMEFloat fby     = invab * (epfrcy - fproj*raby);
        PMEFloat fbz     = invab * (epfrcz - fproj*rabz);
#ifdef EP_VIRIAL
        v11 += (rabx * invab)*epfrcx - rabx*fbx;
        v22 += (raby * invab)*epfrcy - raby*fby;
        v33 += (rabz * invab)*epfrcz - rabz*fbz;
#endif
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[epfrm.x],
		  llitoulli(llrintf(epfrcx - fbx)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[epfrm.x],
		  llitoulli(llrintf(epfrcy - fby)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[epfrm.x],
		  llitoulli(llrintf(epfrcz - fbz)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[epfrm.y],
		  llitoulli(llrintf(fbx)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[epfrm.y],
		  llitoulli(llrintf(fby)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[epfrm.y],
		  llitoulli(llrintf(fbz)));
      }
      else if (frmtype == 10) {
#ifdef EP_NEIGHBORLIST
        PMEFloat parentx = cSim.pImageX[epfrm.x];
        PMEFloat parenty = cSim.pImageY[epfrm.x];
        PMEFloat parentz = cSim.pImageZ[epfrm.x];
        PMEFloat fratm1x = cSim.pImageX[epfrm.y];
        PMEFloat fratm1y = cSim.pImageY[epfrm.y];
        PMEFloat fratm1z = cSim.pImageZ[epfrm.y];
        PMEFloat fratm2x = cSim.pImageX[epfrm.z];
        PMEFloat fratm2y = cSim.pImageY[epfrm.z];
        PMEFloat fratm2z = cSim.pImageZ[epfrm.z];
        PMEFloat fratm3x = cSim.pImageX[epfrm.w];
        PMEFloat fratm3y = cSim.pImageY[epfrm.w];
        PMEFloat fratm3z = cSim.pImageZ[epfrm.w];
#else
        PMEFloat parentx = cSim.pAtomX[epfrm.x];
        PMEFloat parenty = cSim.pAtomY[epfrm.x];
        PMEFloat parentz = cSim.pAtomZ[epfrm.x];
        PMEFloat fratm1x = cSim.pAtomX[epfrm.y];
        PMEFloat fratm1y = cSim.pAtomY[epfrm.y];
        PMEFloat fratm1z = cSim.pAtomZ[epfrm.y];
        PMEFloat fratm2x = cSim.pAtomX[epfrm.z];
        PMEFloat fratm2y = cSim.pAtomY[epfrm.z];
        PMEFloat fratm2z = cSim.pAtomZ[epfrm.z];
        PMEFloat fratm3x = cSim.pAtomX[epfrm.w];
        PMEFloat fratm3y = cSim.pAtomY[epfrm.w];
        PMEFloat fratm3z = cSim.pAtomZ[epfrm.w];
#endif
        PMEFloat epfrcx   = (PMEFloat)forceX;
        PMEFloat epfrcy   = (PMEFloat)forceY;
        PMEFloat epfrcz   = (PMEFloat)forceZ;
        PMEFloat rabx     = (PMEFloat)(fratm1x - parentx);
        PMEFloat raby     = (PMEFloat)(fratm1y - parenty);
        PMEFloat rabz     = (PMEFloat)(fratm1z - parentz);
        PMEFloat racx     = (PMEFloat)(fratm2x - parentx);
        PMEFloat racy     = (PMEFloat)(fratm2y - parenty);
        PMEFloat racz     = (PMEFloat)(fratm2z - parentz);
        PMEFloat radx     = (PMEFloat)(fratm3x - parentx);
        PMEFloat rady     = (PMEFloat)(fratm3y - parenty);
        PMEFloat radz     = (PMEFloat)(fratm3z - parentz);
	PMEFloat fpd1     = (PMEFloat)d1;
        PMEFloat rjax     = fpd1*racx - rabx;
        PMEFloat rjay     = fpd1*racy - raby;
        PMEFloat rjaz     = fpd1*racz - rabz;
	PMEFloat d2       = cSim.pExtraPointCustomD2[pos];
        PMEFloat rjbx     = d2*radx - rabx;
        PMEFloat rjby     = d2*rady - raby;
        PMEFloat rjbz     = d2*radz - rabz;
        PMEFloat rjabx    = rjbx - rjax;
        PMEFloat rjaby    = rjby - rjay;
        PMEFloat rjabz    = rjbz - rjaz;
	PMEFloat rmx      = rjay*rjbz - rjaz*rjby;
	PMEFloat rmy      = rjaz*rjbx - rjax*rjbz;
	PMEFloat rmz      = rjax*rjby - rjay*rjbx;
	PMEFloat magdvec  = (PMEFloat)1.0 * rsqrt(rmx*rmx + rmy*rmy + rmz*rmz);
	PMEFloat magdvec2 = magdvec * magdvec;
	PMEFloat d3       = cSim.pExtraPointCustomD3[pos];
	PMEFloat cfx      = d3 * magdvec * epfrcx;
	PMEFloat cfy      = d3 * magdvec * epfrcy;
	PMEFloat cfz      = d3 * magdvec * epfrcz;
	PMEFloat rtx      = (rmy*rjabz - rmz*rjaby) * magdvec2;
        PMEFloat rty      = (rmz*rjabx - rmx*rjabz) * magdvec2;
        PMEFloat rtz      = (rmx*rjaby - rmy*rjabx) * magdvec2;
        PMEFloat fbx      = -rmx*rtx*cfx + (rjabz - rmy*rtx)*cfy - (rjaby + rmz*rtx)*cfz;
        PMEFloat fby      = -(rjabz + rmx*rty)*cfx - rmy*rty*cfy + (rjabx - rmz*rty)*cfz;
        PMEFloat fbz      = (rjaby - rmx*rtz)*cfx - (rjabx + rmy*rtz)*cfy - rmz*rtz*cfz;
        rtx               = (rjby*rmz - rjbz*rmy) * magdvec2 * fpd1;
        rty               = (rjbz*rmx - rjbx*rmz) * magdvec2 * fpd1;
        rtz               = (rjbx*rmy - rjby*rmx) * magdvec2 * fpd1;
        PMEFloat fcx      = -rmx*rtx*cfx - (fpd1*rjbz + rmy*rtx)*cfy +
                            (fpd1*rjby - rmz*rtx)*cfz;
        PMEFloat fcy      = (fpd1*rjbz - rmx*rty)*cfx - rmy*rty*cfy -
                            (fpd1*rjbx + rmz*rty)*cfz;
        PMEFloat fcz      = -(fpd1*rjby + rmx*rtz)*cfx + (fpd1*rjbx - rmy*rtz)*cfy -
                            rmz*rtz*cfz;
        rtx               = (rmy*rjaz - rmz*rjay) * magdvec2 * d2;
        rty               = (rmz*rjax - rmx*rjaz) * magdvec2 * d2;
        rtz               = (rmx*rjay - rmy*rjax) * magdvec2 * d2;
        PMEFloat fdx      = -rmx*rtx*cfx + (d2*rjaz - rmy*rtx)*cfy - (d2*rjay + rmz*rtx)*cfz;
        PMEFloat fdy      = -(d2*rjaz + rmx*rty)*cfx - rmy*rty*cfy + (d2*rjax - rmz*rty)*cfz;
        PMEFloat fdz      = (d2*rjay - rmx*rtz)*cfx - (d2*rjax + rmy*rtz)*cfy - rmz*rtz*cfz;
#ifdef EP_VIRIAL
#  ifdef EP_NEIGHBORLIST
        PMEFloat raEPx   = cSim.pImageX[epidx] - parentx;
        PMEFloat raEPy   = cSim.pImageY[epidx] - parenty;
        PMEFloat raEPz   = cSim.pImageZ[epidx] - parentz;
#  else
        PMEFloat raEPx  = cSim.pAtomX[epidx] - parentx;
        PMEFloat raEPy  = cSim.pAtomY[epidx] - parenty;
        PMEFloat raEPz  = cSim.pAtomZ[epidx] - parentz;
#  endif
        v11            += raEPx*epfrcx - rabx*fbx - racx*fcx - radx*fdx;
        v22            += raEPy*epfrcy - raby*fby - racy*fcy - rady*fdy;
        v33            += raEPz*epfrcz - rabz*fbz - racz*fcz - radz*fdz;
#endif
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[epfrm.x],
		  llitoulli(llrintf(epfrcx - fbx - fcx - fdx)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[epfrm.x],
		  llitoulli(llrintf(epfrcy - fby - fcy - fdy)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[epfrm.x],
		  llitoulli(llrintf(epfrcz - fbz - fcz - fdz)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[epfrm.y],
		  llitoulli(llrintf(fbx)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[epfrm.y],
		  llitoulli(llrintf(fby)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[epfrm.y],
		  llitoulli(llrintf(fbz)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[epfrm.z],
		  llitoulli(llrintf(fcx)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[epfrm.z],
		  llitoulli(llrintf(fcy)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[epfrm.z],
		  llitoulli(llrintf(fcz)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[epfrm.w],
		  llitoulli(llrintf(fdx)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[epfrm.w],
		  llitoulli(llrintf(fdy)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[epfrm.w],
		  llitoulli(llrintf(fdz)));
      }
    }
    pos += cSim.EP22Offset + blockDim.x*gridDim.x;
  }

#ifdef EP_VIRIAL
  sV[threadIdx.x]  = llrint(v11);
  sV[threadIdx.x] += sV[threadIdx.x ^ 1];
  sV[threadIdx.x] += sV[threadIdx.x ^ 2];
  sV[threadIdx.x] += sV[threadIdx.x ^ 4];
  sV[threadIdx.x] += sV[threadIdx.x ^ 8];
  sV[threadIdx.x] += sV[threadIdx.x ^ 16];
#ifdef AMBER_PLATFORM_AMD_WARP64
  sV[threadIdx.x] += sV[threadIdx.x ^ 32];
#endif
  if ((threadIdx.x & GRID_BITS_MASK) == 0) {
    atomicAdd(cSim.pVirial_11, llitoulli(sV[threadIdx.x]));
  }
  sV[threadIdx.x]  = llrint(v22);
  sV[threadIdx.x] += sV[threadIdx.x ^ 1];
  sV[threadIdx.x] += sV[threadIdx.x ^ 2];
  sV[threadIdx.x] += sV[threadIdx.x ^ 4];
  sV[threadIdx.x] += sV[threadIdx.x ^ 8];
  sV[threadIdx.x] += sV[threadIdx.x ^ 16];
#ifdef AMBER_PLATFORM_AMD_WARP64
  sV[threadIdx.x] += sV[threadIdx.x ^ 32];
#endif
  if ((threadIdx.x & GRID_BITS_MASK) == 0) {
    atomicAdd(cSim.pVirial_22, llitoulli(sV[threadIdx.x]));
  }
  sV[threadIdx.x]  = llrint(v33);
  sV[threadIdx.x] += sV[threadIdx.x ^ 1];
  sV[threadIdx.x] += sV[threadIdx.x ^ 2];
  sV[threadIdx.x] += sV[threadIdx.x ^ 4];
  sV[threadIdx.x] += sV[threadIdx.x ^ 8];
  sV[threadIdx.x] += sV[threadIdx.x ^ 16];
#ifdef AMBER_PLATFORM_AMD_WARP64
  sV[threadIdx.x] += sV[threadIdx.x ^ 32];
#endif
  if ((threadIdx.x & GRID_BITS_MASK) == 0) {
    atomicAdd(cSim.pVirial_33, llitoulli(sV[threadIdx.x]));
  }
#endif
}
