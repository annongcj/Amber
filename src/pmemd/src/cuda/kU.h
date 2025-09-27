#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// This is included by kForcesUpdate to integrate the equations of motion (move the atoms), and
// can be thought of as "k-Update.h."
//
// #defines: UPDATE_NEIGHBORLIST, UPDATE_LANGEVIN, UPDATE_RELAXMD
//---------------------------------------------------------------------------------------------
{
#ifdef UPDATE_NEIGHBORLIST
#define MASS(i) cSim.pImageMass[i]
#define INVMASS(i) cSim.pImageInvMass[i]
#define VELX(i) cSim.pImageVelX[i]
#define VELY(i) cSim.pImageVelY[i]
#define VELZ(i) cSim.pImageVelZ[i]
#define LVELX(i) cSim.pImageLVelX[i]
#define LVELY(i) cSim.pImageLVelY[i]
#define LVELZ(i) cSim.pImageLVelZ[i]
#define ATOMX(i) cSim.pImageX[i]
#define ATOMY(i) cSim.pImageY[i]
#define ATOMZ(i) cSim.pImageZ[i]
#ifdef UPDATE_SGLD
#define WTSG(i) cSim.pImageWsg[i]
#define X0SG(i) cSim.pImageX0sg[i]
#define Y0SG(i) cSim.pImageY0sg[i]
#define Z0SG(i) cSim.pImageZ0sg[i]
#define X1SG(i) cSim.pImageX1sg[i]
#define Y1SG(i) cSim.pImageY1sg[i]
#define Z1SG(i) cSim.pImageZ1sg[i]
#define X2SG(i) cSim.pImageX2sg[i]
#define Y2SG(i) cSim.pImageY2sg[i]
#define Z2SG(i) cSim.pImageZ2sg[i]
#ifdef UPDATE_LANGEVIN
#define RSGX(i) cSim.pImageRsgX[i]
#define RSGY(i) cSim.pImageRsgY[i]
#define RSGZ(i) cSim.pImageRsgZ[i]
#endif
#define PPSG(i) cSim.pImagePPsg[i]
#define FPSG(i) cSim.pImageFPsg[i]
#endif
#else
#define MASS(i) cSim.pAtomMass[i]
#define INVMASS(i) cSim.pAtomInvMass[i]
#define VELX(i) cSim.pVelX[i]
#define VELY(i) cSim.pVelY[i]
#define VELZ(i) cSim.pVelZ[i]
#define LVELX(i) cSim.pLVelX[i]
#define LVELY(i) cSim.pLVelY[i]
#define LVELZ(i) cSim.pLVelZ[i]
#define ATOMX(i) cSim.pAtomX[i]
#define ATOMY(i) cSim.pAtomY[i]
#define ATOMZ(i) cSim.pAtomZ[i]
#ifdef UPDATE_SGLD
#define WTSG(i) cSim.pWsg[i]
#define X0SG(i) cSim.pX0sg[i]
#define Y0SG(i) cSim.pY0sg[i]
#define Z0SG(i) cSim.pZ0sg[i]
#define X1SG(i) cSim.pX1sg[i]
#define Y1SG(i) cSim.pY1sg[i]
#define Z1SG(i) cSim.pZ1sg[i]
#define X2SG(i) cSim.pX2sg[i]
#define Y2SG(i) cSim.pY2sg[i]
#define Z2SG(i) cSim.pZ2sg[i]
#ifdef UPDATE_LANGEVIN
#define RSGX(i) cSim.pRsgX[i]
#define RSGY(i) cSim.pRsgY[i]
#define RSGZ(i) cSim.pRsgZ[i]
#endif
#define PPSG(i) cSim.pPPsg[i]
#define FPSG(i) cSim.pFPsg[i]
#endif
#endif

  double dtx       = dt * 20.455;
  double half_dtx = dtx * 0.5;
#ifdef UPDATE_LANGEVIN
  double gammai    = gamma_ln / 20.455;
  double c_implic  = 1.0 / (1.0 + gammai * half_dtx);
  double c_explic  = 1.0 - gammai * half_dtx;
  double sdfac     = 4.0 * gammai * boltz2 * temp0 / dtx;
#endif
#ifdef UPDATE_SGLD
        double tsgfac                               = cSim.sgavg/dtx;
        double sgavg1                               = cSim.sgavg;
        double sgavp1                               = cSim.sgavp;
        double sgavg0                               = 1.0-sgavg1;
        double sgavp0                               = 1.0-sgavp1;
        double wgam0                                = 0.0;
        if(cSim.isgld==1)wgam0           = 1.0;
        double wgam1                                = 1.0 - wgam0;
#endif
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  if (pos < cSim.atoms) {
    double aamass  = MASS(pos);

#ifdef UPDATE_LANGEVIN
    if (cSim.doTIHeating) {
      unsigned realAtomindex = cSim.pImageAtom[pos];
      if (cSim.pTIList[realAtomindex] > 0 || cSim.pTIList[realAtomindex + cSim.stride] > 0) {
        gammai *= 100.0;
        sdfac = 4.0 * gammai * boltz2 * cSim.TIHeatingTemp / dtx;  
      }
      else {
        gammai /= 100.0;
        sdfac = 4.0 * gammai * boltz2 * temp0 / dtx;  
      }
    }
    else {
      sdfac = 4.0 * gammai * boltz2 * temp0 / dtx;
    } 
    c_implic = 1.0 / (1.0 + gammai * half_dtx);
    c_explic = 1.0 - gammai * half_dtx;
#endif

    double atomX   = ATOMX(pos);
    double atomY   = ATOMY(pos);
    double atomZ   = ATOMZ(pos);
    double invMass = INVMASS(pos);
#ifdef UPDATE_RELAXMD
    unsigned int index = cSim.pImageAtom[pos];
#endif
    PMEAccumulator fx = cSim.pForceXAccumulator[pos];
    PMEAccumulator fy = cSim.pForceYAccumulator[pos];
    PMEAccumulator fz = cSim.pForceZAccumulator[pos];
#if defined(UPDATE_NTP) && !defined(MPI)
    PMEAccumulator nfx = cSim.pNBForceXAccumulator[pos];
    PMEAccumulator nfy = cSim.pNBForceYAccumulator[pos];
    PMEAccumulator nfz = cSim.pNBForceZAccumulator[pos];
    double forceX = (double)(fx + nfx) * (double)ONEOVERFORCESCALE;
    double forceY = (double)(fy + nfy) * (double)ONEOVERFORCESCALE;
    double forceZ = (double)(fz + nfz) * (double)ONEOVERFORCESCALE;
#else
    double forceX = (double)fx * (double)ONEOVERFORCESCALE;
    double forceY = (double)fy * (double)ONEOVERFORCESCALE;
    double forceZ = (double)fz * (double)ONEOVERFORCESCALE;
#endif
    double velX  = VELX(pos);
    double velY  = VELY(pos);
    double velZ  = VELZ(pos);
#ifdef UPDATE_LANGEVIN
    double gaussX = cSim.pRandomX[pos + rpos];
    double gaussY = cSim.pRandomY[pos + rpos];
    double gaussZ = cSim.pRandomZ[pos + rpos];
    double rsd    = sqrt(sdfac * aamass);
#endif
    double wfac = invMass * dtx;
    
    // Save previous velocities and positions
    LVELX(pos) = velX;
    LVELY(pos) = velY;
    LVELZ(pos) = velZ;
    cSim.pOldAtomX[pos] = atomX;
    cSim.pOldAtomY[pos] = atomY;
    cSim.pOldAtomZ[pos] = atomZ;
    // Update velocities
#ifdef UPDATE_RELAXMD
    if (index >= cSim.first_update_atom) {
#endif
#ifdef UPDATE_SGLD
      //if(WTSG(pos)){
        double sgwt  =WTSG(pos)>0?1.0:0.0;
        double gamsg                                = wgam0*gamma_sg - wgam1*FPSG(pos)/PPSG(pos);
        // Update velocities     
        // SGLD local averages
        double X0sg                                 = X0SG(pos);
        double Y0sg                                 = Y0SG(pos);
        double Z0sg                                 = Z0SG(pos);
        double X1sg                                 = X1SG(pos);
        double Y1sg                                 = Y1SG(pos);
        double Z1sg                                 = Z1SG(pos);
        double PsgX0                                 = tsgfac*aamass*(X0sg-X1sg);
        double PsgY0                                 = tsgfac*aamass*(Y0sg-Y1sg);
        double PsgZ0                                 = tsgfac*aamass*(Z0sg-Z1sg);
        double Xtran                                 = atomX-velX * dtx -X0sg;
        double Ytran                                 = atomY-velY * dtx -Y0sg;
        double Ztran                                 = atomZ-velZ * dtx -Z0sg;
        X1sg                                 = sgavg0*(X1sg+dcomx1+Xtran)+sgavg1*atomX;
        Y1sg                                 = sgavg0*(Y1sg+dcomy1+Ytran)+sgavg1*atomY;
        Z1sg                                 = sgavg0*(Z1sg+dcomz1+Ztran)+sgavg1*atomZ;
        double X2sg                                 = sgavg0*(X2SG(pos)+dcomx2+Xtran)+sgavg1*X1sg;
        double Y2sg                                 = sgavg0*(Y2SG(pos)+dcomy2+Ytran)+sgavg1*Y1sg;
        double Z2sg                                 = sgavg0*(Z2SG(pos)+dcomz2+Ztran)+sgavg1*Z1sg;
        //double PsgX                                 = sgavg0*PSGX(pos)+sgavg1*aamass*velX;
        //double PsgY                                 = sgavg0*PSGY(pos)+sgavg1*aamass*velY;
        //double PsgZ                                 = sgavg0*PSGZ(pos)+sgavg1*aamass*velZ;
        double PsgX                                 = tsgfac*aamass*(atomX-X1sg);
        double PsgY                                 = tsgfac*aamass*(atomY-Y1sg);
        double PsgZ                                 = tsgfac*aamass*(atomZ-Z1sg);
        double ptX                                  = (PsgX-sgavg0*PsgX0)/sgavg1;
        double ptY                                  = (PsgY-sgavg0*PsgY0)/sgavg1;
        double ptZ                                  = (PsgZ-sgavg0*PsgZ0)/sgavg1;

        //  guiding force
        //avgdfi3=tsgfac*(pi3t-2.0d0*avgpi3+tsgfac*amassi*(x1i3-x2i3))
	      //double fsgX                                 =tsgfac*(aamass*velX-2.0*PsgX+tsgfac*aamass*(X1sg-X2sg));
	      //double fsgY                                 =tsgfac*(aamass*velY-2.0*PsgY+tsgfac*aamass*(Y1sg-Y2sg));
	      //double fsgZ                                 =tsgfac*(aamass*velZ-2.0*PsgZ+tsgfac*aamass*(Z1sg-Z2sg));
	      double fsgX                                 =tsgfac*(ptX-2.0*PsgX+tsgfac*aamass*(X1sg-X2sg));
	      double fsgY                                 =tsgfac*(ptY-2.0*PsgY+tsgfac*aamass*(Y1sg-Y2sg));
	      double fsgZ                                 =tsgfac*(ptZ-2.0*PsgZ+tsgfac*aamass*(Z1sg-Z2sg));

        double wsgp                                 = (cSim.sgft*gamsg);
	      // calculate guiding forces
	      double guidX                                =  sgwt*( wsgp*PsgX+cSim.sgff*fsgX);
        double guidY                                =  sgwt*( wsgp*PsgY+cSim.sgff*fsgY);
        double guidZ                                =  sgwt*( wsgp*PsgZ+cSim.sgff*fsgZ);
        // update collision arrays
        double sumpp                           = (PsgX*PsgX+PsgY*PsgY+PsgZ*PsgZ);
        double sumfp                           =fsgX*PsgX+fsgY*PsgY+fsgZ*PsgZ;
        double sumpv                           =aamass*(velX*velX+velY*velY+velZ*velZ);
        double sumgv                           =guidX*velX+guidY*velY+guidZ*velZ;
        double PPsg                            =sgavp0*PPSG(pos)+sgavp1*sumpp;
        double FPsg                            =sgavp0*FPSG(pos)+sgavp1*sumfp;
        //double PPsg                            =sumpp;
        //double FPsg                            =sumfp;
	  // add guiding force 
	      forceX                                      += guidX;
	      forceY                                      += guidY;
	      forceZ                                      += guidZ;
#ifdef UPDATE_LANGEVIN
        double wsgg                                 = (cSim.sgfg*gammai);
	      // KSG is the local average random forces integrator
	      double RsgX                                 = sgavg0*RSGX(pos)+sgavg1*rsd * gaussX;
        double RsgY                                 = sgavg0*RSGY(pos)+sgavg1*rsd * gaussY;
        double RsgZ                                 = sgavg0*RSGZ(pos)+sgavg1*rsd * gaussZ;

	      // calculate GLE guiding forces
	      double gleX                                =  wsgg*PsgX+cSim.fsgldg*RsgX;
        double gleY                                =  wsgg*PsgY+cSim.fsgldg*RsgY;
        double gleZ                                =  wsgg*PsgZ+cSim.fsgldg*RsgZ;
	      forceX                                     += gleX;
	      forceY                                     += gleY;
	      forceZ                                     += gleZ;
        double sgbeta                          = sumgv/(sumpv/c_implic-sumgv*half_dtx);
        double fact                            =(gammai+sgbeta)*half_dtx;
        velX                                   =((1.0-fact)*velX+(forceX + rsd*gaussX)*wfac)/(1.0+fact);
        velY                                   =((1.0-fact)*velY+(forceY + rsd*gaussY)*wfac)/(1.0+fact);
        velZ                                   =((1.0-fact)*velZ+(forceZ + rsd*gaussZ)*wfac)/(1.0+fact);
        RSGX(pos)=RsgX;
        RSGY(pos)=RsgY;
        RSGZ(pos)=RsgZ;
#else
        double sgbeta                          =sumgv/(sumpv-sumgv*half_dtx);
        double fact                            = sumgv/(sumpv);
        forceX                                -=fact*aamass*(velX+0.5*forceX*wfac);
        forceY                                -=fact*aamass*(velY+0.5*forceY*wfac);
        forceZ                                -=fact*aamass*(velZ+0.5*forceZ*wfac);
        velX += forceX * wfac;
        velY += forceY * wfac;
        velZ += forceZ * wfac;
#endif
      X0SG(pos)=atomX;
      Y0SG(pos)=atomY;
      Z0SG(pos)=atomZ;
      X1SG(pos)=X1sg;
      Y1SG(pos)=Y1sg;
      Z1SG(pos)=Z1sg;
      X2SG(pos)=X2sg;
      Y2SG(pos)=Y2sg;
      Z2SG(pos)=Z2sg;
      PPSG(pos)=PPsg;
      FPSG(pos)=FPsg;
    //}
    //else{
#else  //UPDATE_SGLD
#ifdef UPDATE_LANGEVIN
      velX = (velX*c_explic + (forceX + rsd*gaussX)*wfac) * c_implic;
      velY = (velY*c_explic + (forceY + rsd*gaussY)*wfac) * c_implic;
      velZ = (velZ*c_explic + (forceZ + rsd*gaussZ)*wfac) * c_implic;
#else
      velX += forceX * wfac;
      velY += forceY * wfac;
      velZ += forceZ * wfac;
#endif
#endif  //UPDATE_SGLD
#ifdef UPDATE_SGLD
    //}
#endif
      if (cSim.vlimit > 0) {
        if (isnan(velX) || isinf(velX)) {
          velX = 0.0;
        }
        velX = max(min(velX, cSim.vlimit), -cSim.vlimit);
        if (isnan(velY) || isinf(velY)) {
          velY = 0.0;
        }
        velY = max(min(velY, cSim.vlimit), -cSim.vlimit);
        if (isnan(velZ) || isinf(velZ)) {
          velZ = 0.0;
        }
        velZ = max(min(velZ, cSim.vlimit), -cSim.vlimit);
      }

      // Save new velocity
      VELX(pos)= velX;
      VELY(pos)= velY;
      VELZ(pos)= velZ;

      // Update positions for SHAKE and kinetic energy kernel
      double newAtomX = atomX + velX*dtx;
      double newAtomY = atomY + velY*dtx;
      double newAtomZ = atomZ + velZ*dtx;
      ATOMX(pos) = newAtomX;
      ATOMY(pos) = newAtomY;
      ATOMZ(pos) = newAtomZ;
#ifndef UPDATE_NEIGHBORLIST
      PMEFloat2 xy;
      xy.x = newAtomX;
      xy.y = newAtomY;
      cSim.pAtomXYSP[pos] = xy;
      cSim.pAtomZSP[pos]  = newAtomZ;
#endif

#ifdef UPDATE_RELAXMD
    }
#endif
  }
#undef MASS
#undef INVMASS
#undef VELX
#undef VELY
#undef VELZ
#undef LVELX
#undef LVELY
#undef LVELZ
#undef ATOMX
#undef ATOMY
#undef ATOMZ
#ifdef UPDATE_SGLD
#undef WTSG
#undef X0SG
#undef Y0SG
#undef Z0SG
#undef X1SG
#undef Y1SG
#undef Z1SG
#undef X2SG
#undef Y2SG
#undef Z2SG
#undef RSGX
#undef RSGY
#undef RSGZ
#undef FPSG
#undef PPSG
#endif
}
