#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// This is included by kForcesUpdate to adjust velocities based on forces and Langevin
// collisions.  Say "k-Reset Velocities.h."
//
// #defines: RV_NEIGHBORLIST, RV_NTP
//---------------------------------------------------------------------------------------------
{
#ifdef RV_NEIGHBORLIST
#  define INVMASS(i) cSim.pImageInvMass[i]
#  define VELX(i) cSim.pImageVelX[i]
#  define VELY(i) cSim.pImageVelY[i]
#  define VELZ(i) cSim.pImageVelZ[i]
#else
#  define INVMASS(i) cSim.pAtomInvMass[i]
#  define VELX(i) cSim.pVelX[i]
#  define VELY(i) cSim.pVelY[i]
#  define VELZ(i) cSim.pVelZ[i]
#endif
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  double boltz     = 8.31441e-3 * temp / 4.184;
  if (pos < cSim.atoms) {
    double invMass    = INVMASS(pos);
    PMEAccumulator fx = cSim.pForceXAccumulator[pos];
    PMEAccumulator fy = cSim.pForceYAccumulator[pos];
    PMEAccumulator fz = cSim.pForceZAccumulator[pos];
#if defined(RV_NTP) && !defined(MPI)
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
    double velX, velY, velZ;

    // Zero velocities if it's really cold
    if (temp < 1.0e-6) {
      velX = (double)0.0;
      velY = (double)0.0;
      velZ = (double)0.0;
    }
    else {
      double gaussX = cSim.pRandomX[pos + rpos];
      double gaussY = cSim.pRandomY[pos + rpos];
      double gaussZ = cSim.pRandomZ[pos + rpos];
      double sd     = sqrt(boltz * invMass);
      velX = sd * gaussX;
      velY = sd * gaussY;
      velZ = sd * gaussZ;
    }

    // Back velocities up a half-step
    double wfac = invMass * half_dtx;
    velX -= forceX * wfac;
    velY -= forceY * wfac;
    velZ -= forceZ * wfac;

    // Write final velocities
    VELX(pos) = velX;
    VELY(pos) = velY;
    VELZ(pos) = velZ;
  }
#undef INVMASS
#undef VELX
#undef VELY
#undef VELZ
}
