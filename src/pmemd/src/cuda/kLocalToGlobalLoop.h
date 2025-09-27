#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// This code is included by kLocalToGlobalKernel.h to handle extra points.
//---------------------------------------------------------------------------------------------
#ifdef EP_NEIGHBORLIST
#  define PATOMX(p) cSim.pImageX[p]
#  define PATOMY(p) cSim.pImageY[p]
#  define PATOMZ(p) cSim.pImageZ[p]
#  if defined(EP_ONEPOINT)
#    define EP1 index
#    if defined(EP_TYPE1)
  int4 frame = cSim.pImageExtraPoint11Frame[pos];
  int index  = cSim.pImageExtraPoint11Index[pos];
#      define PLOCALX1(p) cSim.pExtraPoint11X[p]
#      define PLOCALY1(p) cSim.pExtraPoint11Y[p]
#      define PLOCALZ1(p) cSim.pExtraPoint11Z[p]
#    elif defined(EP_TYPE2)
  int4 frame = cSim.pImageExtraPoint12Frame[pos];
  int index  = cSim.pImageExtraPoint12Index[pos];
#      define PLOCALX1(p) cSim.pExtraPoint12X[p]
#      define PLOCALY1(p) cSim.pExtraPoint12Y[p]
#      define PLOCALZ1(p) cSim.pExtraPoint12Z[p]
#    endif
#  elif defined(EP_TWOPOINTS)
#    define EP1 index.x
#    if defined(EP_TYPE1)
  int4 frame = cSim.pImageExtraPoint21Frame[pos];
  int2 index = cSim.pImageExtraPoint21Index[pos];
#      define PLOCALX1(p) cSim.pExtraPoint21X1[p]
#      define PLOCALY1(p) cSim.pExtraPoint21Y1[p]
#      define PLOCALZ1(p) cSim.pExtraPoint21Z1[p]
#      define PLOCALX2(p) cSim.pExtraPoint21X2[p]
#      define PLOCALY2(p) cSim.pExtraPoint21Y2[p]
#      define PLOCALZ2(p) cSim.pExtraPoint21Z2[p]
#    elif defined(EP_TYPE2)
  int4 frame = cSim.pImageExtraPoint22Frame[pos];
  int2 index = cSim.pImageExtraPoint22Index[pos];
#      define PLOCALX1(p) cSim.pExtraPoint22X1[p]
#      define PLOCALY1(p) cSim.pExtraPoint22Y1[p]
#      define PLOCALZ1(p) cSim.pExtraPoint22Z1[p]
#      define PLOCALX2(p) cSim.pExtraPoint22X2[p]
#      define PLOCALY2(p) cSim.pExtraPoint22Y2[p]
#      define PLOCALZ2(p) cSim.pExtraPoint22Z2[p]
#    endif // End branch over defined extra point types
#  endif // End branch over the number of defined extra points
#else  // EP_NEIGHBORLIST
#  define PATOMX(p) cSim.pAtomX[p]
#  define PATOMY(p) cSim.pAtomY[p]
#  define PATOMZ(p) cSim.pAtomZ[p]
#  if defined(EP_ONEPOINT)
#    define EP1 index
#    if defined(EP_TYPE1)
  int4 frame = cSim.pExtraPoint11Frame[pos];
  int index  = cSim.pExtraPoint11Index[pos];
#      define PLOCALX1(p) cSim.pExtraPoint11X[p]
#      define PLOCALY1(p) cSim.pExtraPoint11Y[p]
#      define PLOCALZ1(p) cSim.pExtraPoint11Z[p]
#    elif defined(EP_TYPE2)
  int4 frame = cSim.pExtraPoint12Frame[pos];
  int index  = cSim.pExtraPoint12Index[pos];
#      define PLOCALX1(p) cSim.pExtraPoint12X[p]
#      define PLOCALY1(p) cSim.pExtraPoint12Y[p]
#      define PLOCALZ1(p) cSim.pExtraPoint12Z[p]
#    endif
#  elif defined(EP_TWOPOINTS)
#    define EP1 index.x
#    if defined(EP_TYPE1)
  int4 frame = cSim.pExtraPoint21Frame[pos];
  int2 index = cSim.pExtraPoint21Index[pos];
#      define PLOCALX1(p) cSim.pExtraPoint21X1[p]
#      define PLOCALY1(p) cSim.pExtraPoint21Y1[p]
#      define PLOCALZ1(p) cSim.pExtraPoint21Z1[p]
#      define PLOCALX2(p) cSim.pExtraPoint21X2[p]
#      define PLOCALY2(p) cSim.pExtraPoint21Y2[p]
#      define PLOCALZ2(p) cSim.pExtraPoint21Z2[p]
#    elif defined(EP_TYPE2)
  int4 frame = cSim.pExtraPoint22Frame[pos];
  int2 index = cSim.pExtraPoint22Index[pos];
#      define PLOCALX1(p) cSim.pExtraPoint22X1[p]
#      define PLOCALY1(p) cSim.pExtraPoint22Y1[p]
#      define PLOCALZ1(p) cSim.pExtraPoint22Z1[p]
#      define PLOCALX2(p) cSim.pExtraPoint22X2[p]
#      define PLOCALY2(p) cSim.pExtraPoint22Y2[p]
#      define PLOCALZ2(p) cSim.pExtraPoint22Z2[p]
#    endif // End branch over defined extra point types
#  endif // End branch over the number of defined extra points
#endif // EP_NEIGHBORLIST

#ifdef EP_TYPE1
  PMEDouble ax = PATOMX(frame.y);
  PMEDouble ay = PATOMY(frame.y);
  PMEDouble az = PATOMZ(frame.y);
  PMEDouble bx = PATOMX(frame.z);
  PMEDouble by = PATOMY(frame.z);
  PMEDouble bz = PATOMZ(frame.z);
  PMEDouble cx = PATOMX(frame.w);
  PMEDouble cy = PATOMY(frame.w);
  PMEDouble cz = PATOMZ(frame.w);
#elif defined(EP_TYPE2)
  PMEDouble px = PATOMX(frame.z);
  PMEDouble py = PATOMY(frame.z);
  PMEDouble pz = PATOMZ(frame.z);
  PMEDouble ax = PATOMX(frame.y);
  PMEDouble ay = PATOMY(frame.y);
  PMEDouble az = PATOMZ(frame.y);
  PMEDouble cx = PATOMX(frame.w);
  PMEDouble cy = PATOMY(frame.w);
  PMEDouble cz = PATOMZ(frame.w);
  PMEDouble bx = PATOMX(frame.x);
  PMEDouble by = PATOMY(frame.x);
  PMEDouble bz = PATOMZ(frame.x);
  ax = (PMEDouble)0.5 * (ax + px);
  ay = (PMEDouble)0.5 * (ay + py);
  az = (PMEDouble)0.5 * (az + pz);
  cx = (PMEDouble)0.5 * (cx + px);
  cy = (PMEDouble)0.5 * (cy + py);
  cz = (PMEDouble)0.5 * (cz + pz);
#endif // End branch over defined extra point types
  PMEDouble uvecx = ax - bx;
  PMEDouble uvecy = ay - by;
  PMEDouble uvecz = az - bz;
  PMEDouble usiz  = rsqrt(uvecx*uvecx + uvecy*uvecy + uvecz*uvecz);
  PMEDouble vvecx = cx - bx;
  PMEDouble vvecy = cy - by;
  PMEDouble vvecz = cz - bz;
  PMEDouble vsiz  = rsqrt(vvecx*vvecx + vvecy*vvecy + vvecz*vvecz);
  uvecx *= usiz;
  uvecy *= usiz;
  uvecz *= usiz;
  vvecx *= vsiz;
  vvecy *= vsiz;
  vvecz *= vsiz;
  PMEDouble avex  = (PMEDouble)0.5 * (uvecx + vvecx);
  PMEDouble avey  = (PMEDouble)0.5 * (uvecy + vvecy);
  PMEDouble avez  = (PMEDouble)0.5 * (uvecz + vvecz);
  PMEDouble diffx = (PMEDouble)0.5 * (vvecx - uvecx);
  PMEDouble diffy = (PMEDouble)0.5 * (vvecy - uvecy);
  PMEDouble diffz = (PMEDouble)0.5 * (vvecz - uvecz);
  PMEDouble asiz  = rsqrt(avex*avex + avey*avey + avez*avez);
  PMEDouble dsiz  = rsqrt(diffx*diffx + diffy*diffy + diffz*diffz);
  PMEDouble f13 = -avex * asiz;
  PMEDouble f23 = -avey * asiz;
  PMEDouble f33 = -avez * asiz;
  PMEDouble f11 = diffx * dsiz;
  PMEDouble f21 = diffy * dsiz;
  PMEDouble f31 = diffz * dsiz;
  PMEDouble f12 = f23*f31 - f33*f21;
  PMEDouble f22 = f33*f11 - f13*f31;
  PMEDouble f32 = f13*f21 - f23*f11;

  // Handle 1st EP
  PMEDouble x1 = PLOCALX1(pos);
  PMEDouble y1 = PLOCALY1(pos);
  PMEDouble z1 = PLOCALZ1(pos);
#undef PLOCALX1
#undef PLOCALY1
#undef PLOCALZ1
#ifdef EP_TWOPOINTS

  // Handle second EP
  PMEDouble x2 = PLOCALX2(pos);
  PMEDouble y2 = PLOCALY2(pos);
  PMEDouble z2 = PLOCALZ2(pos);
#  undef PLOCALX2
#  undef PLOCALY2
#  undef PLOCALZ2
#endif // EP_TWOPOINTS
  PATOMX(EP1) = bx + f11*x1 + f12*y1 + f13*z1;
  PATOMY(EP1) = by + f21*x1 + f22*y1 + f23*z1;
  PATOMZ(EP1) = bz + f31*x1 + f32*y1 + f33*z1;
#undef EP1
#ifdef EP_TWOPOINTS
  PATOMX(index.y) = bx + f11*x2 + f12*y2 + f13*z2;
  PATOMY(index.y) = by + f21*x2 + f22*y2 + f23*z2;
  PATOMZ(index.y) = bz + f31*x2 + f32*y2 + f33*z2;
#endif
#undef PATOMX
#undef PATOMY
#undef PATOMZ
