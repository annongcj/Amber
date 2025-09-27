#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------
#ifdef EP_NEIGHBORLIST
#  define PATOMX(p) cSim.pImageX[p]
#  define PATOMY(p) cSim.pImageY[p]
#  define PATOMZ(p) cSim.pImageZ[p]
#  if defined(EP_ONEPOINT)
#    define EP1 index
#    if defined(EP_TYPE1)
  int index  = cSim.pImageExtraPoint11Index[pos];
  int4 frame = cSim.pImageExtraPoint11Frame[pos];
#    elif defined(EP_TYPE2)
  int index  = cSim.pImageExtraPoint12Index[pos];
  int4 frame = cSim.pImageExtraPoint12Frame[pos];
#    endif
#  elif defined(EP_TWOPOINTS)
#    define EP1 index.x
#    if defined(EP_TYPE1)
  int2 index = cSim.pImageExtraPoint21Index[pos];
  int4 frame = cSim.pImageExtraPoint21Frame[pos];
#    elif defined(EP_TYPE2)
  int2 index = cSim.pImageExtraPoint22Index[pos];
  int4 frame = cSim.pImageExtraPoint22Frame[pos];
#    endif
#  endif
#else
#  define PATOMX(p) cSim.pAtomX[p]
#  define PATOMY(p) cSim.pAtomY[p]
#  define PATOMZ(p) cSim.pAtomZ[p]
#  if defined(EP_ONEPOINT)
#    define EP1 index
#    if defined(EP_TYPE1)
  int index  = cSim.pExtraPoint11Index[pos];
  int4 frame = cSim.pExtraPoint11Frame[pos];
#    elif defined(EP_TYPE2)
  int index  = cSim.pExtraPoint12Index[pos];
  int4 frame = cSim.pExtraPoint12Frame[pos];
#    endif
#  elif defined(EP_TWOPOINTS)
#    define EP1 index.x
#    if defined(EP_TYPE1)
  int2 index = cSim.pExtraPoint21Index[pos];
  int4 frame = cSim.pExtraPoint21Frame[pos];
#    elif defined(EP_TYPE2)
  int2 index = cSim.pExtraPoint22Index[pos];
  int4 frame = cSim.pExtraPoint22Frame[pos];
#    endif
#  endif
#endif

  // Handle 1st or only EP
  PMEDouble x1 = PATOMX(EP1);
  PMEDouble y1 = PATOMY(EP1);
  PMEDouble z1 = PATOMZ(EP1);
  PMEDouble xp = PATOMX(frame.x);
  PMEDouble yp = PATOMY(frame.x);
  PMEDouble zp = PATOMZ(frame.x);
  PMEDouble forceX = cSim.pNBForceXAccumulator[EP1];
  PMEDouble forceY = cSim.pNBForceYAccumulator[EP1];
  PMEDouble forceZ = cSim.pNBForceZAccumulator[EP1];
  PMEDouble rx = x1 - xp;
  PMEDouble ry = y1 - yp;
  PMEDouble rz = z1 - zp;
  cSim.pForceXAccumulator[EP1] = (PMEAccumulator)0;
  cSim.pForceYAccumulator[EP1] = (PMEAccumulator)0;
  cSim.pForceZAccumulator[EP1] = (PMEAccumulator)0;
#ifdef EP_VIRIAL
  cSim.pNBForceXAccumulator[EP1] = (PMEAccumulator)0;
  cSim.pNBForceYAccumulator[EP1] = (PMEAccumulator)0;
  cSim.pNBForceZAccumulator[EP1] = (PMEAccumulator)0;
  v11 += forceX * rx;
  v22 += forceY * ry;
  v33 += forceZ * rz;
#endif
  PMEDouble torqueX = ry*forceZ - rz*forceY;
  PMEDouble torqueY = rz*forceX - rx*forceZ;
  PMEDouble torqueZ = rx*forceY - ry*forceX;
#ifdef EP_TWOPOINTS

  // Handle 2nd EP
  PMEDouble x2 = PATOMX(index.y);
  PMEDouble y2 = PATOMY(index.y);
  PMEDouble z2 = PATOMZ(index.y);
  PMEDouble forceX2 = cSim.pNBForceXAccumulator[index.y];
  PMEDouble forceY2 = cSim.pNBForceYAccumulator[index.y];
  PMEDouble forceZ2 = cSim.pNBForceZAccumulator[index.y];
  rx = x2 - xp;
  ry = y2 - yp;
  rz = z2 - zp;
  cSim.pForceXAccumulator[index.y] = (PMEAccumulator)0;
  cSim.pForceYAccumulator[index.y] = (PMEAccumulator)0;
  cSim.pForceZAccumulator[index.y] = (PMEAccumulator)0;
#ifdef EP_VIRIAL
  cSim.pNBForceXAccumulator[index.y] = (PMEAccumulator)0;
  cSim.pNBForceYAccumulator[index.y] = (PMEAccumulator)0;
  cSim.pNBForceZAccumulator[index.y] = (PMEAccumulator)0;
  v11 += forceX2 * rx;
  v22 += forceY2 * ry;
  v33 += forceZ2 * rz;
#endif
  forceX  += forceX2;
  forceY  += forceY2;
  forceZ  += forceZ2;
  torqueX += ry * forceZ2 - rz * forceY2;
  torqueY += rz * forceX2 - rx * forceZ2;
  torqueZ += rx * forceY2 - ry * forceX2;
#endif

#ifdef EP_TYPE1
  PMEDouble apx = PATOMX(frame.y);
  PMEDouble apy = PATOMY(frame.y);
  PMEDouble apz = PATOMZ(frame.y);
  PMEDouble bpx = PATOMX(frame.z);
  PMEDouble bpy = PATOMY(frame.z);
  PMEDouble bpz = PATOMZ(frame.z);
  PMEDouble cpx = PATOMX(frame.w);
  PMEDouble cpy = PATOMY(frame.w);
  PMEDouble cpz = PATOMZ(frame.w);
#elif defined(EP_TYPE2)
  PMEDouble apx = PATOMX(frame.y);
  PMEDouble apy = PATOMY(frame.y);
  PMEDouble apz = PATOMZ(frame.y);
  PMEDouble px  = PATOMX(frame.z);
  PMEDouble py  = PATOMY(frame.z);
  PMEDouble pz  = PATOMZ(frame.z);
  PMEDouble cpx = PATOMX(frame.w);
  PMEDouble cpy = PATOMY(frame.w);
  PMEDouble cpz = PATOMZ(frame.w);
  PMEDouble bpx = PATOMX(frame.x);
  PMEDouble bpy = PATOMY(frame.x);
  PMEDouble bpz = PATOMZ(frame.x);
  apx = (PMEDouble)0.5 * (apx + px);
  apy = (PMEDouble)0.5 * (apy + py);
  apz = (PMEDouble)0.5 * (apz + pz);
  cpx = (PMEDouble)0.5 * (cpx + px);
  cpy = (PMEDouble)0.5 * (cpy + py);
  cpz = (PMEDouble)0.5 * (cpz + pz);
#endif
  PMEDouble ux   = apx - bpx;
  PMEDouble uy   = apy - bpy;
  PMEDouble uz   = apz - bpz;
  PMEDouble usiz = rsqrt(ux*ux + uy*uy + uz*uz);
  PMEDouble vx   = cpx - bpx;
  PMEDouble vy   = cpy - bpy;
  PMEDouble vz   = cpz - bpz;
  PMEDouble vsiz = rsqrt(vx*vx + vy*vy + vz*vz);
  PMEDouble wx   = uy*vz - uz*vy;
  PMEDouble wy   = uz*vx - ux*vz;
  PMEDouble wz   = ux*vy - uy*vx;
  ux *= usiz;
  uy *= usiz;
  uz *= usiz;
  vx *= vsiz;
  vy *= vsiz;
  vz *= vsiz;
  PMEDouble wsiz = rsqrt(wx*wx + wy*wy + wz*wz);
  wx *= wsiz;
  wy *= wsiz;
  wz *= wsiz;
  PMEDouble dx  = vx - ux;
  PMEDouble dy  = vy - uy;
  PMEDouble dz  = vz - uz;
  PMEDouble dotdu  = ux*dx + uy*dy + uz*dz;
  PMEDouble dotdv  = vx*dx + vy*dy + vz*dz;
  PMEDouble upx = dx - dotdu*ux;
  PMEDouble upy = dy - dotdu*uy;
  PMEDouble upz = dz - dotdu*uz;
  PMEDouble vpx = dx - dotdv*vx;
  PMEDouble vpy = dy - dotdv*vy;
  PMEDouble vpz = dz - dotdv*vz;
  PMEDouble upsiz  = rsqrt(upx*upx + upy*upy + upz*upz);
  upx *= upsiz;
  upy *= upsiz;
  upz *= upsiz;
  PMEDouble vpsiz  = rsqrt(vpx*vpx + vpy*vpy + vpz*vpz);
  vpx *= vpsiz;
  vpy *= vpsiz;
  vpz *= vpsiz;
  PMEDouble c = ux*vx + uy*vy + uz*vz;
  PMEDouble s = rsqrt((PMEDouble)1.0 - c*c);
  PMEDouble uvdis  = usiz*s;
  PMEDouble vudis  = vsiz*s;
  PMEDouble dphidu = -(torqueX*ux + torqueY*uy + torqueZ*uz);
  PMEDouble dphidv = -(torqueX*vx + torqueY*vy + torqueZ*vz);
  PMEDouble dphidw = -(torqueX*wx + torqueY*wy + torqueZ*wz);
  dphidv *= uvdis;
  dphidu *= vudis;
  usiz *= (PMEDouble)0.5 * dphidw;
  vsiz *= (PMEDouble)0.5 * dphidw;
  PMEDouble dux = -wx*dphidv + upx*usiz;
  PMEDouble duy = -wy*dphidv + upy*usiz;
  PMEDouble duz = -wz*dphidv + upz*usiz;
  PMEDouble dvx =  wx*dphidu + vpx*vsiz;
  PMEDouble dvy =  wy*dphidu + vpy*vsiz;
  PMEDouble dvz =  wz*dphidu + vpz*vsiz;

#ifdef EP_VIRIAL
  // Get torque contribution to virial:
  v11 += dux*(apx - bpx) + dvx*(cpx - bpx);
  v22 += duy*(apy - bpy) + dvy*(cpy - bpy);
  v33 += duz*(apz - bpz) + dvz*(cpz - bpz);
#endif
#ifdef EP_TYPE1
  atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[frame.y],
            llitoulli(llrint(-dux)));
  atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[frame.y],
            llitoulli(llrint(-duy)));
  atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[frame.y],
            llitoulli(llrint(-duz)));
  atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[frame.w],
            llitoulli(llrint(-dvx)));
  atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[frame.w],
            llitoulli(llrint(-dvy)));
  atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[frame.w],
            llitoulli(llrint(-dvz)));
  atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[frame.z],
            llitoulli(llrint(dvx + dux + forceX)));
  atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[frame.z],
            llitoulli(llrint(dvy + duy + forceY)));
  atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[frame.z],
            llitoulli(llrint(dvz + duz + forceZ)));
#elif defined(EP_TYPE2)
  atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[frame.x],
            llitoulli(llrint(dvx + dux + forceX)));
  atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[frame.x],
            llitoulli(llrint(dvy + duy + forceY)));
  atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[frame.x],
            llitoulli(llrint(dvz + duz + forceZ)));
  atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[frame.y],
            llitoulli(llrint(-(PMEDouble)0.5 * dux)));
  atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[frame.y],
            llitoulli(llrint(-(PMEDouble)0.5 * duy)));
  atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[frame.y],
            llitoulli(llrint(-(PMEDouble)0.5 * duz)));
  atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[frame.w],
            llitoulli(llrint(-(PMEDouble)0.5 * dvx)));
  atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[frame.w],
            llitoulli(llrint(-(PMEDouble)0.5 * dvy)));
  atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[frame.w],
            llitoulli(llrint(-(PMEDouble)0.5 * dvz)));
  atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[frame.z],
            llitoulli(llrint(-(PMEDouble)0.5 * (dvx + dux))));
  atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[frame.z],
            llitoulli(llrint(-(PMEDouble)0.5 * (dvy + duy))));
  atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[frame.z],
            llitoulli(llrint(-(PMEDouble)0.5 * (dvz + duz))));
#endif
#undef EP1
#undef PATOMX
#undef PATOMY
#undef PATOMZ
