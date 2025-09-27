#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// Nov 2021, by zhf
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// This is included by kForcesUpdate to integrate the equations of motion (move the atoms), and
// can be thought of as "k-Update.h."
//
// #defines: UPDATE_NEIGHBORLIST, UPDATE_LANGEVIN, UPDATE_RELAXMD, UPDATE_MIDDLE_SCHEME
//---------------------------------------------------------------------------------------------
{
#ifdef UPDATE_NEIGHBORLIST
#  define MASS(i) cSim.pImageMass[i]
#  define INVMASS(i) cSim.pImageInvMass[i]
#  define VELX(i) cSim.pImageVelX[i]
#  define VELY(i) cSim.pImageVelY[i]
#  define VELZ(i) cSim.pImageVelZ[i]
#  define LVELX(i) cSim.pImageLVelX[i]
#  define LVELY(i) cSim.pImageLVelY[i]
#  define LVELZ(i) cSim.pImageLVelZ[i]
#  define ATOMX(i) cSim.pImageX[i]
#  define ATOMY(i) cSim.pImageY[i]
#  define ATOMZ(i) cSim.pImageZ[i]
#else
#  define MASS(i) cSim.pAtomMass[i]
#  define INVMASS(i) cSim.pAtomInvMass[i]
#  define VELX(i) cSim.pVelX[i]
#  define VELY(i) cSim.pVelY[i]
#  define VELZ(i) cSim.pVelZ[i]
#  define LVELX(i) cSim.pLVelX[i]
#  define LVELY(i) cSim.pLVelY[i]
#  define LVELZ(i) cSim.pLVelZ[i]
#  define ATOMX(i) cSim.pAtomX[i]
#  define ATOMY(i) cSim.pAtomY[i]
#  define ATOMZ(i) cSim.pAtomZ[i]
#endif

    double dtx       = dt * 20.455;
    double half_dtx  = dtx * 0.5;
#ifdef UPDATE_MIDDLE_SCHEME_2
    double lgv_c1    = exp(-gamma_ln * dt);
    double lgv_c2    = sqrt(1.0 - lgv_c1 * lgv_c1);
    double rtKT      = sqrt(2.0 * boltz2 * temp0);
#endif
    
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

    if (pos < cSim.atoms) {
#ifdef UPDATE_MIDDLE_SCHEME_2
        double atomX    = ATOMX(pos);
        double atomY    = ATOMY(pos);
        double atomZ    = ATOMZ(pos);
#endif

#ifdef UPDATE_MIDDLE_SCHEME_1
        double invMass  = INVMASS(pos);

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

        double wfac = invMass * dtx;

        LVELX(pos) = VELX(pos);
        LVELY(pos) = VELY(pos);
        LVELZ(pos) = VELZ(pos);
        cSim.pOldAtomX[pos] = ATOMX(pos);
        cSim.pOldAtomY[pos] = ATOMY(pos);
        cSim.pOldAtomZ[pos] = ATOMZ(pos);


        VELX(pos) += forceX * wfac;
        VELY(pos) += forceY * wfac;
        VELZ(pos) += forceZ * wfac;


#endif

#ifdef UPDATE_MIDDLE_SCHEME_2

        double velX = VELX(pos);
        double velY = VELY(pos);
        double velZ = VELZ(pos);
        double stdvel = rtKT * sqrt(INVMASS(pos));

        double newAtomX = atomX + velX * half_dtx;
        double newAtomY = atomY + velY * half_dtx;
        double newAtomZ = atomZ + velZ * half_dtx;

        double gaussX   = cSim.pRandomX[pos + rpos];
        double gaussY   = cSim.pRandomY[pos + rpos];
        double gaussZ   = cSim.pRandomZ[pos + rpos];


        velX = lgv_c1 * velX + lgv_c2 * gaussX * stdvel;
        velY = lgv_c1 * velY + lgv_c2 * gaussY * stdvel;
        velZ = lgv_c1 * velZ + lgv_c2 * gaussZ * stdvel;

        VELX(pos) = velX;
        VELY(pos) = velY;
        VELZ(pos) = velZ;

        ATOMX(pos) = newAtomX + velX * half_dtx;
        ATOMY(pos) = newAtomY + velY * half_dtx;
        ATOMZ(pos) = newAtomZ + velZ * half_dtx;

#ifndef UPDATE_NEIGHBORLIST
      PMEFloat2 xy;
      xy.x = newAtomX;
      xy.y = newAtomY;
      cSim.pAtomXYSP[pos] = xy;
      cSim.pAtomZSP[pos]  = newAtomZ;
#endif
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
}


