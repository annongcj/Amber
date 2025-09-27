#ifdef GTI
#include "gputypes.h"
#include "gti_utils.cuh"
#include "gti_energies_kernels.cuh"

#include "simulationConst.h"
CSIM_STO simulationConst cSim;

namespace GTI_ENERGY_IMPL {
  __constant__ const PMEDouble tm06 = (PMEDouble)(1.0e-06);
  __constant__ const PMEDouble tenm3 = (PMEDouble)(1.0e-03);
  __constant__ const PMEDouble tm24 = (PMEDouble)(1.0e-18);
  __constant__ const PMEDouble one = (PMEDouble)(1.0);
  __constant__ const PMEDouble zero = (PMEDouble)(0.0);
  __constant__ const PMEDouble nlimit = -0.99999;
  __constant__ const PMEDouble plimit = 0.99999;

  #if !(defined(__CUDACC_RDC__) || defined(__HIPCC_RDC__))
    #include "gti_localCUutils.inc"
  #endif
}

using namespace GTI_ENERGY_IMPL;

//---------------------------------------------------------------------------------------------
// kgCalculateTIKineticEnergy_kernel:
//
// Arguments:
//   useImage:
//   c_ave:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgCalculateTIKineticEnergy_kernel(bool useImage, double c_ave)
{
  // Only use one block
  unsigned int pos = threadIdx.x;
  unsigned int increment = blockDim.x;

  __shared__ unsigned long long int seke[2], sekph[2], sekpbs[2], sekmh[2], sSCEke[2];
  if (pos == 0) {
    seke[0] = sekph[0] = sekpbs[0] = sekmh[0] = sSCEke[0] = Zero;
    seke[1] = sekph[1] = sekpbs[1] = sekmh[1] = sSCEke[1] = Zero;
  }

  unsigned long long int eke[2] = { Zero, Zero }, ekph[2] = { Zero, Zero }, ekpbs[2] = { Zero, Zero }, ekmh[2] = { Zero, Zero }, SCEke[2] = { Zero, Zero };

  __syncthreads();

  PMEDouble mass, vx, lvx, vy, lvy, vz, lvz;
  while (pos < cSim.numberTIAtoms) {
    unsigned int atom = cSim.pTIAtomList[pos].x;
    unsigned TIRegion = cSim.pTIAtomList[pos].y;
    bool isSC = (cSim.pTIAtomList[pos].z>0);
    unsigned long long int regionShift = TIRegion*cSim.GPUKinEnergyTerms;

    if (useImage) {
      unsigned iatom = cSim.pImageAtomLookup[atom];
      mass = cSim.pImageMass[iatom];
      vx   = cSim.pImageVelX[iatom];
      lvx  = cSim.pImageLVelX[iatom];
      vy   = cSim.pImageVelY[iatom];
      lvy  = cSim.pImageLVelY[iatom];
      vz   = cSim.pImageVelZ[iatom];
      lvz  = cSim.pImageLVelZ[iatom];
    }
    else {
      mass = cSim.pAtomMass[atom];
      vx   = cSim.pVelX[atom];
      lvx  = cSim.pLVelX[atom];
      vy   = cSim.pVelY[atom];
      lvy  = cSim.pLVelY[atom];
      vz   = cSim.pVelZ[atom];
      lvz  = cSim.pLVelZ[atom];
    }

    PMEDouble svx   = vx + lvx;
    PMEDouble svy   = vy + lvy;
    PMEDouble svz   = vz + lvz;
    PMEDouble ke   = c_ave * 0.125 * mass * (svx*svx + svy*svy + svz*svz);
    PMEDouble kpbs = 0.5 * mass * (vx*lvx + vy*lvy + vz*lvz);
    PMEDouble kph  = 0.5 * mass * (vx*vx + vy*vy + vz*vz);
    PMEDouble kmh  = 0.5 * mass * (lvx*lvx + lvy*lvy + lvz*lvz);

    addLocalEnergy(&eke[TIRegion], ke);
    addLocalEnergy(&ekpbs[TIRegion], kpbs);
    addLocalEnergy(&ekph[TIRegion], kph);
    addLocalEnergy(&ekmh[TIRegion], kmh);
    if (isSC) {
      addLocalEnergy(&SCEke[TIRegion], ke);
    }

    pos += increment;
  }

  __syncthreads();

  for (unsigned TIRegion = 0; TIRegion < 2; TIRegion++) {
    for (int offset = warpSize / 2; offset > 0; offset /= 2) {
#ifndef AMBER_PLATFORM_AMD 
      __syncwarp();
#endif
      eke[TIRegion] += __SHFL_DOWN(0xFFFFFFFF, eke[TIRegion], offset);
      ekpbs[TIRegion] += __SHFL_DOWN(0xFFFFFFFF, ekpbs[TIRegion], offset);
      ekph[TIRegion] += __SHFL_DOWN(0xFFFFFFFF, ekph[TIRegion], offset);
      ekmh[TIRegion] += __SHFL_DOWN(0xFFFFFFFF, ekmh[TIRegion], offset);
      SCEke[TIRegion] += __SHFL_DOWN(0xFFFFFFFF, SCEke[TIRegion], offset);
    }
  }

  __syncthreads();

  pos = threadIdx.x;
  if ((pos & (warpSize - 1)) == 0) {
    for (unsigned TIRegion = 0; TIRegion < 2; TIRegion++) {
      atomicAdd(&seke[TIRegion], eke[TIRegion]);
      atomicAdd(&sekph[TIRegion], ekph[TIRegion]);
      atomicAdd(&sekpbs[TIRegion], ekpbs[TIRegion]);
      atomicAdd(&sekmh[TIRegion], ekmh[TIRegion]);
      atomicAdd(&sSCEke[TIRegion], SCEke[TIRegion]);
    }
  }

  __syncthreads();

  if (pos == 0) {
    for (unsigned TIRegion = 0; TIRegion < 2; TIRegion++) {
      unsigned long long int regionShift = TIRegion * cSim.GPUKinEnergyTerms;
      cSim.pTIKinEnergyBuffer[regionShift] = seke[TIRegion];
      cSim.pTIKinEnergyBuffer[regionShift + 2] = sekph[TIRegion];
      cSim.pTIKinEnergyBuffer[regionShift + 3] = sekpbs[TIRegion];
      cSim.pTIKinEnergyBuffer[regionShift + 4] = sekmh[TIRegion];
      cSim.pTIKinEnergyBuffer[regionShift + 1] = sSCEke[TIRegion];
    }
  }
}

//---------------------------------------------------------------------------------------------
// kgBondedEnergy_kernel:
//
// Arguments:
//   needPotEnergy:
//   needvirial:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgBondedEnergy_kernel(bool needPotEnergy, bool needvirial)
{
  __shared__ PMEDouble  recip[9], ucell[9];

  int pos = blockIdx.x * blockDim.x + threadIdx.x;
  int i= threadIdx.x;
  if (cSim.pNTPData == NULL) {
    if (i < 3) {
      recip[i * 3] = cSim.recip[i][0];
      recip[i * 3 + 1] = cSim.recip[i][1];
      recip[i * 3 + 2] = cSim.recip[i][2];
      ucell[i * 3] = cSim.ucell[i][0];
      ucell[i * 3 + 1] = cSim.ucell[i][1];
      ucell[i * 3 + 2] = cSim.ucell[i][2];
    }
  }
  else {
    if (i < 9) {
      ucell[i] = cSim.pNTPData->ucell[i];
      recip[i] = cSim.pNTPData->recip[i];
    }
  }
  __syncthreads();

  //unsigned int totalTerms = cSim.numberTIBond + cSim.numberTIAngle + cSim.numberTIDihedral;

  int nwarp_bond = ((cSim.numberTIBond - 1) / GRID + 1) * GRID;
  int nwarp_angle = ((cSim.numberTIAngle - 1) / GRID + 1) * GRID;
  int nwarp_dihedral = ((cSim.numberTIDihedral - 1) / GRID + 1) * GRID;

  int nwarp_NMRDistance = ((cSim.numberTINMRDistance - 1) / GRID + 1) * GRID;
  int nwarp_NMRAngle = ((cSim.numberTINMRAngle - 1) / GRID + 1) * GRID;
  //int nwarp_NMRDihedral = ((cSim.numberTINMRDihedral - 1) / GRID + 1) * GRID;

  int bondBegin = 0;
  int bondEnd = cSim.numberTIBond - 1;

  int angleBegin = nwarp_bond;
  int angleEnd = angleBegin + cSim.numberTIAngle - 1;

  int dihedralBegin = angleBegin + nwarp_angle;
  int dihedralEnd = dihedralBegin + cSim.numberTIDihedral - 1;

  int NMRDistanceBegin = dihedralBegin + nwarp_dihedral;
  int NMRDistanceEnd = NMRDistanceBegin + cSim.numberTINMRDistance - 1;

  int NMRAngleBegin = NMRDistanceBegin + nwarp_NMRDistance;
  int NMRAngleEnd = NMRAngleBegin + cSim.numberTINMRAngle - 1;

  int NMRDihedralBegin = NMRAngleBegin + nwarp_NMRAngle;
  int NMRDihedralEnd = NMRDihedralBegin + cSim.numberTINMRDihedral - 1;

  PMEAccumulator* pTIFx[3] = { cSim.pTIForce,
                        cSim.pTIForce + cSim.stride3,
                        cSim.pTIForce + (cSim.stride3 * 2) };
  PMEAccumulator* pTIFy[3] = { cSim.pTIForce + cSim.stride,
                        cSim.pTIForce + cSim.stride + cSim.stride3,
                        cSim.pTIForce + cSim.stride + (cSim.stride3 * 2) };
  PMEAccumulator* pTIFz[3] = { cSim.pTIForce + (cSim.stride * 2),
                        cSim.pTIForce + (cSim.stride * 2) + cSim.stride3,
                        cSim.pTIForce + (cSim.stride * 2) +
                        (cSim.stride3 * 2) };

  PMEAccumulator* pSCFx[2] = { cSim.pTIForce + (cSim.stride3 * 3),
                        cSim.pTIForce + (cSim.stride3 * 4) };
  PMEAccumulator* pSCFy[2] = { cSim.pTIForce + (cSim.stride3 * 3) + cSim.stride,
                        cSim.pTIForce + (cSim.stride3 * 4) + cSim.stride };
  PMEAccumulator* pSCFz[2] = { cSim.pTIForce + (cSim.stride3 * 3) +
                       (cSim.stride * 2), cSim.pTIForce +
                       (cSim.stride3 * 4) + (cSim.stride * 2) };

  unsigned iatom[4], TIRegion;
  PMEDouble xij, xkj, xkl, yij, ykj, ykl, zij, zkj, zkl;
  PMEDouble ff, da, df, dfw, tt, tw;
  PMEDouble rij2, rij, rkj2, rik_i, cst;
  PMEDouble dx, dy, dz, gx, gy, gz;
  PMEDouble dcdx, dcdy, dcdz, dcgx, dcgy, dcgz;
  PMEDouble fxi, fyi, z1, z2, z11, z12, z22;
  PMEDouble sphi, cphi;
  PMEDouble dr1, dr2, dr3, dr4, dr5, dr6, drx, dry, drz;
  PMEDouble cik, sth, cii, ckk;
  PMEDouble fx1, fy1, fz1, fx2, fy2, fz2;

  unsigned long long int *pE0, *pE;
  unsigned long long int *pSE1[2] = { cSim.pTISpecialTerms[0][0], cSim.pTISpecialTerms[0][1] };
  unsigned long long int *pSE2[2] = { cSim.pTISpecialTerms[1][0], cSim.pTISpecialTerms[1][1] };

  PMEAccumulator* pFx, *pFy, *pFz;

  if (pos >= bondBegin && pos <= bondEnd) {
    unsigned mypos = pos - bondBegin;
    uint2& bondID = cSim.pTIBondID[mypos];
    PMEDouble2& bond = cSim.pTIBond[mypos];
    iatom[0] = cSim.pImageAtomLookup[bondID.x];
    iatom[1] = cSim.pImageAtomLookup[bondID.y];

    TIRegion = (cSim.pTIBondType[mypos] & gti_simulationConst::Fg_TI_region);
    bool hasSC = (cSim.pTIBondType[mypos] & gti_simulationConst::Fg_has_SC);
    bool inSC = (cSim.pTIBondType[mypos] & gti_simulationConst::Fg_int_SC);
    bool addToDVDL = (cSim.pTIBondType[mypos] & gti_simulationConst::ex_addToDVDL_bat);
    bool addToBATCorr = (cSim.pTIBondType[mypos] & gti_simulationConst::ex_addToBat_corr);
    PMEDouble myWeight = cSim.TIItemWeight[Schedule::TypeBAT][TIRegion];

    xij = cSim.pImageX[iatom[0]] - cSim.pImageX[iatom[1]];
    yij = cSim.pImageY[iatom[0]] - cSim.pImageY[iatom[1]];
    zij = cSim.pImageZ[iatom[0]] - cSim.pImageZ[iatom[1]];

    rij2 = __image_dist2<PMEDouble>(xij, yij, zij, recip, ucell);
    rij = sqrt(rij2);
    da = rij - bond.y;
    df = bond.x * da;
    dfw = (df + df) / rij;

    if (addToDVDL) {
      pE = cSim.pTIEBond[TIRegion];
      pFx = pTIFx[TIRegion];
      pFy = pTIFy[TIRegion];
      pFz = pTIFz[TIRegion];
    } else {
      pE = cSim.pTISCBond[TIRegion];
      pFx = pSCFx[TIRegion];
      pFy = pSCFy[TIRegion];
      pFz = pSCFz[TIRegion];
    }

    if (needPotEnergy) {
      tt = df * da;
      addEnergy(pE, tt);

      if (hasSC && !inSC) {
        addEnergy(pSE1[TIRegion], tt);
        if (!addToBATCorr) {
          addEnergy(pSE2[TIRegion], tt);
        }
      }
    }

    // Add to the force accumulators

    tw = addToDVDL ? myWeight : OneF;
    if (fabs(tw) > tm06) {
      tt = dfw * xij * tw;
      addForce(pFx + iatom[0], -tt);
      addForce(pFx + iatom[1], tt);

      tt = dfw * yij * tw;
      addForce(pFy + iatom[0], -tt);
      addForce(pFy + iatom[1], tt);

      tt = dfw * zij * tw;
      addForce(pFz + iatom[0], -tt);
      addForce(pFz + iatom[1], tt);
    }

  }
  if (pos >= angleBegin && pos <= angleEnd) {
    unsigned mypos = pos - angleBegin;
    uint4& angleID = cSim.pTIBondAngleID[mypos];
    PMEDouble2& angle = cSim.pTIBondAngle[mypos];
    iatom[0] = cSim.pImageAtomLookup[angleID.x];
    iatom[1] = cSim.pImageAtomLookup[angleID.y];
    iatom[2] = cSim.pImageAtomLookup[angleID.z];

    TIRegion = (cSim.pTIBondAngleType[mypos] & gti_simulationConst::Fg_TI_region);
    bool hasSC = (cSim.pTIBondAngleType[mypos] & gti_simulationConst::Fg_has_SC);
    bool inSC = (cSim.pTIBondAngleType[mypos] & gti_simulationConst::Fg_int_SC);
    bool addToDVDL = (cSim.pTIBondAngleType[mypos] & gti_simulationConst::ex_addToDVDL_bat);
    bool addToBATCorr = (cSim.pTIBondAngleType[mypos] & gti_simulationConst::ex_addToBat_corr);
    PMEDouble myWeight = cSim.TIItemWeight[Schedule::TypeBAT][TIRegion];

    xij = cSim.pImageX[iatom[0]] - cSim.pImageX[iatom[1]];
    yij = cSim.pImageY[iatom[0]] - cSim.pImageY[iatom[1]];
    zij = cSim.pImageZ[iatom[0]] - cSim.pImageZ[iatom[1]];
    xkj = cSim.pImageX[iatom[2]] - cSim.pImageX[iatom[1]];
    ykj = cSim.pImageY[iatom[2]] - cSim.pImageY[iatom[1]];
    zkj = cSim.pImageZ[iatom[2]] - cSim.pImageZ[iatom[1]];
    rij = xij * xij + yij * yij + zij * zij;
    PMEDouble rkj = xkj * xkj + ykj * ykj + zkj * zkj;
    PMEDouble rik = sqrt(rij * rkj);

    const PMEDouble nlimit = -0.99999;
    const PMEDouble plimit = 0.99999;
    PMEDouble cst = min(plimit, max(nlimit, (xij*xkj + yij * ykj + zij * zkj) / rik));
    PMEDouble ant = acos(cst);

    da = ant - angle.y;
    df = angle.x * da;
    dfw = -(df + df) / sin(ant);

    pE0 = cSim.pTIEAngle[2];
    if (addToDVDL) {
      pE = cSim.pTIEAngle[TIRegion];
      pFx = pTIFx[TIRegion];
      pFy = pTIFy[TIRegion];
      pFz = pTIFz[TIRegion];
    } else {
      pE = cSim.pTISCAngle[TIRegion];
      pFx = pSCFx[TIRegion];
      pFy = pSCFy[TIRegion];
      pFz = pSCFz[TIRegion];
    }

    if (needPotEnergy) {
      tt = df * da;
      addEnergy(pE, tt);

      if (hasSC && !inSC) {
        addEnergy(pSE1[TIRegion], tt);
        if (!addToBATCorr) {
          addEnergy(pSE2[TIRegion], tt);
        }
      }
    }

    // Calculation of the force
    cik = dfw / rik;
    sth = dfw * cst;
    cii = sth / rij;
    ckk = sth / rkj;
    fx1 = cii * xij - cik * xkj;
    fy1 = cii * yij - cik * ykj;
    fz1 = cii * zij - cik * zkj;
    fx2 = ckk * xkj - cik * xij;
    fy2 = ckk * ykj - cik * yij;
    fz2 = ckk * zkj - cik * zij;

    tw = addToDVDL ? myWeight : 1.0;
    if (fabs(tw) > tm06) {
      addForce(pFx + iatom[0], fx1 * tw);
      addForce(pFy + iatom[0], fy1 * tw);
      addForce(pFz + iatom[0], fz1 * tw);

      addForce(pFx + iatom[2], fx2 * tw);
      addForce(pFy + iatom[2], fy2 * tw);
      addForce(pFz + iatom[2], fz2 * tw);

      addForce(pFx + iatom[1], -tw * (fx1 + fx2));
      addForce(pFy + iatom[1], -tw * (fy1 + fy2));
      addForce(pFz + iatom[1], -tw * (fz1 + fz2));
    }
  }
  if (pos >= dihedralBegin && pos <= dihedralEnd) {
    unsigned mypos = pos - dihedralBegin;
    uint4& dihedralID = cSim.pTIDihedralID[mypos];
    iatom[0] = cSim.pImageAtomLookup[dihedralID.x];
    iatom[1] = cSim.pImageAtomLookup[dihedralID.y];
    iatom[2] = cSim.pImageAtomLookup[dihedralID.z];
    iatom[3] = cSim.pImageAtomLookup[dihedralID.w];

    PMEDouble2& dihedral1 = cSim.pTIDihedral1[mypos];
    PMEDouble2& dihedral2 = cSim. pTIDihedral2[mypos];
    PMEDouble& dihedral3 = cSim.pTIDihedral3[mypos];

    TIRegion = (cSim.pTIDihedralType[mypos] & gti_simulationConst::Fg_TI_region);
    bool addToDVDL = (cSim.pTIDihedralType[mypos] & gti_simulationConst::ex_addToDVDL_bat);
    PMEDouble myWeight = cSim.TIItemWeight[Schedule::TypeBAT][TIRegion];

    xij = cSim.pImageX[iatom[0]] - cSim.pImageX[iatom[1]];
    yij = cSim.pImageY[iatom[0]] - cSim.pImageY[iatom[1]];
    zij = cSim.pImageZ[iatom[0]] - cSim.pImageZ[iatom[1]];
    xkj = cSim.pImageX[iatom[2]] - cSim.pImageX[iatom[1]];
    ykj = cSim.pImageY[iatom[2]] - cSim.pImageY[iatom[1]];
    zkj = cSim.pImageZ[iatom[2]] - cSim.pImageZ[iatom[1]];
    xkl = cSim.pImageX[iatom[2]] - cSim.pImageX[iatom[3]];
    ykl = cSim.pImageY[iatom[2]] - cSim.pImageY[iatom[3]];
    zkl = cSim.pImageZ[iatom[2]] - cSim.pImageZ[iatom[3]];

    // Get the normal vector
    dx = yij * zkj - zij * ykj;
    dy = zij * xkj - xij * zkj;
    dz = xij * ykj - yij * xkj;
    gx = zkj * ykl - ykj * zkl;
    gy = xkj * zkl - zkj * xkl;
    gz = ykj * xkl - xkj * ykl;
    fxi = sqrt(dx*dx + dy * dy + dz * dz + tm24);
    fyi = sqrt(gx*gx + gy * gy + gz * gz + tm24);
    cst = dx * gx + dy * gy + dz * gz;

    // Branch if linear dihedral:
    PMEDouble z1 = (tenm3 <= fxi) ? (1.0 / fxi) : zero;
    PMEDouble z2 = (tenm3 <= fyi) ? (1.0 / fyi) : zero;
    PMEDouble z12 = z1 * z2;
    PMEDouble fzi = (z12 != zero) ? one : zero;
    PMEDouble s = xkj * (dz*gy - dy * gz) + ykj * (dx*gz - dz * gx) + zkj * (dy*gx - dx * gy);
    PMEDouble ff = PI - abs(acos(max(-1.0, min(1.0, cst * z12)))) * (s >= zero ? one : -one);
    faster_sincos(ff, &sphi, &cphi);

    // Calculate the energy and the derivatives with respect to cosphi
    PMEDouble ct0 = dihedral1.y * ff;
    PMEDouble sinnp, cosnp;
    faster_sincos(ct0, &sinnp, &cosnp);

    pE0 = cSim.pTIEDihedral[2];
    if (addToDVDL) {
      pE = cSim.pTIEDihedral[TIRegion];
      pFx = pTIFx[TIRegion];
      pFy = pTIFy[TIRegion];
      pFz = pTIFz[TIRegion];
    } else {
      pE = cSim.pTISCDihedral[TIRegion];
      pFx = pSCFx[TIRegion];
      pFy = pSCFy[TIRegion];
      pFz = pSCFz[TIRegion];
    }

    tt = (dihedral2.x + (cosnp * dihedral2.y) + (sinnp * dihedral3)) * fzi;
    if (needPotEnergy) addEnergy(pE, tt);

#ifdef use_DPFP
    PMEFloat delta = tm24;
#else
    PMEFloat delta = tm06;
#endif
    PMEDouble dums = sphi + (delta * (sphi >= zero ? one : -one));
    PMEDouble df;
    if (tm06 > abs(dums)) {
      df = fzi * dihedral2.y * (dihedral1.y - dihedral1.x + (dihedral1.x * cphi));
    } else {
      df = fzi * dihedral1.y * ((dihedral2.y * sinnp) - (dihedral3 * cosnp)) / dums;
    }

    // Now, set up array dc = 1st der. of cosphi w/respect to cartesian differences:
    z11 = z1 * z1;
    z12 = z1 * z2;
    z22 = z2 * z2;
    dcdx = -gx * z12 - cphi * dx*z11;
    dcdy = -gy * z12 - cphi * dy*z11;
    dcdz = -gz * z12 - cphi * dz*z11;
    dcgx = dx * z12 + cphi * gx*z22;
    dcgy = dy * z12 + cphi * gy*z22;
    dcgz = dz * z12 + cphi * gz*z22;

    // Update the first derivative array:
    dr1 = df * (dcdz*ykj - dcdy * zkj);
    dr2 = df * (dcdx*zkj - dcdz * xkj);
    dr3 = df * (dcdy*xkj - dcdx * ykj);
    dr4 = df * (dcgz*ykj - dcgy * zkj);
    dr5 = df * (dcgx*zkj - dcgz * xkj);
    dr6 = df * (dcgy*xkj - dcgx * ykj);
    drx = df * (-dcdy * zij + dcdz * yij + dcgy * zkl - dcgz * ykl);
    dry = df * (dcdx*zij - dcdz * xij - dcgx * zkl + dcgz * xkl);
    drz = df * (-dcdx * yij + dcdy * xij + dcgx * ykl - dcgy * xkl);

    tw = addToDVDL ? myWeight : 1.0;
    addForce(pFx + iatom[0], (-dr1*tw));
    addForce(pFy + iatom[0], (-dr2 * tw));
    addForce(pFz + iatom[0], (-dr3 * tw));

    addForce(pFx + iatom[1], (-drx + dr1)* tw);
    addForce(pFy + iatom[1], (-dry + dr2)* tw);
    addForce(pFz + iatom[1], (-drz + dr3)* tw);

    addForce(pFx + iatom[2], (drx + dr4)* tw);
    addForce(pFy + iatom[2], (dry + dr5)* tw);
    addForce(pFz + iatom[2], (drz + dr6)* tw);

    addForce(pFx + iatom[3], (-dr4 * tw));
    addForce(pFy + iatom[3], (-dr5 * tw));
    addForce(pFz + iatom[3], (-dr6 * tw));
  }

  if (pos >= NMRDistanceBegin && pos <= NMRDihedralEnd) {

    bool doDistance = (pos >= NMRDistanceBegin && pos <= NMRDistanceEnd);
    bool doAngle = (pos >= NMRAngleBegin && pos <= NMRAngleEnd);
    bool doDihedral = (pos >= NMRDihedralBegin && pos <= NMRDihedralEnd);

    if (!doDistance && !doAngle && !doDihedral) return;

      unsigned long long int* pTIRestEnergy = (
        (doDistance) ? cSim.pTIRestDist[2] : (
        (doAngle) ? cSim.pTIRestAng[2] : (
          (doDihedral) ? cSim.pTIRestTor[2] : 0))
        );

      unsigned long long int* pTIRestDer =  (
        (doDistance) ? cSim.pTIRestDist_DL[2] : (
        (doAngle) ? cSim.pTIRestAng_DL[2] : (
          (doDihedral) ? cSim.pTIRestTor_DL[2] : 0))
        );

      unsigned mypos0 = pos - (
        (doDistance) ? NMRDistanceBegin : (
        (doAngle) ? NMRAngleBegin : (
          (doDihedral) ? NMRDihedralBegin : 0)));

      if (doDistance && mypos0 > cSim.numberTINMRDistance) return;
      if (doAngle && mypos0 > cSim.numberTINMRAngle) return;
      if (doDihedral && mypos0 > cSim.numberTINMRDihedral) return;

      unsigned mypos1 = mypos0 + (
        (doDistance) ? cSim.numberTINMRDistance : (
        (doAngle) ? cSim.numberTINMRAngle : (
          (doDihedral) ? cSim.numberTINMRDihedral : 0)));

      PMEDouble4* pR = (
        (doDistance) ? cSim.pTINMRDistancesR : (
        (doAngle) ? cSim.pTINMRAngleR : (
          (doDihedral) ? cSim.pTINMRDihedralR : NULL)));

      PMEDouble4 R0 = pR[mypos0];
      PMEDouble4 R1 = pR[mypos1];
      PMEDouble4 R = {
        cSim.TIItemWeight[Schedule::TypeRestBA][0] * R0.x + cSim.TIItemWeight[Schedule::TypeRestBA][1] * R1.x,
        cSim.TIItemWeight[Schedule::TypeRestBA][0] * R0.y + cSim.TIItemWeight[Schedule::TypeRestBA][1] * R1.y,
        cSim.TIItemWeight[Schedule::TypeRestBA][0] * R0.z + cSim.TIItemWeight[Schedule::TypeRestBA][1] * R1.z,
        cSim.TIItemWeight[Schedule::TypeRestBA][0] * R0.w + cSim.TIItemWeight[Schedule::TypeRestBA][1] * R1.w
      };

      PMEDouble2* pK = (
        (doDistance) ? cSim.pTINMRDistancesK : (
        (doAngle) ? cSim.pTINMRAngleK : (
          (doDihedral) ? cSim.pTINMRDihedralK : NULL)));

      PMEDouble2 K0 = pK[mypos0];
      PMEDouble2 K1 = pK[mypos1];
      PMEDouble2 K{
      cSim.TIItemWeight[Schedule::TypeRestBA][0] * K0.x + cSim.TIItemWeight[Schedule::TypeRestBA][1] * K1.x,
      cSim.TIItemWeight[Schedule::TypeRestBA][0] * K0.y + cSim.TIItemWeight[Schedule::TypeRestBA][1] * K1.y,
      };

      uint4* pID =
        (doDistance) ? cSim.pTINMRDistancesID : (
        (doAngle) ? cSim.pTINMRAngleID : (
          (doDihedral) ? cSim.pTINMRDihedralID : NULL));

      uint4 ID = pID[mypos0];
      iatom[0] = cSim.pImageAtomLookup[ID.x];
      iatom[1] = cSim.pImageAtomLookup[ID.y];
      if (doAngle || doDihedral) iatom[2] = cSim.pImageAtomLookup[ID.z];
      if (doDihedral) iatom[3] = cSim.pImageAtomLookup[ID.w];

      xij = cSim.pImageX[iatom[0]] - cSim.pImageX[iatom[1]];
      yij = cSim.pImageY[iatom[0]] - cSim.pImageY[iatom[1]];
      zij = cSim.pImageZ[iatom[0]] - cSim.pImageZ[iatom[1]];
      if (doAngle || doDihedral) {
        xkj = cSim.pImageX[iatom[2]] - cSim.pImageX[iatom[1]];
        ykj = cSim.pImageY[iatom[2]] - cSim.pImageY[iatom[1]];
        zkj = cSim.pImageZ[iatom[2]] - cSim.pImageZ[iatom[1]];
        if (doDihedral) {
          xkl = cSim.pImageX[iatom[2]] - cSim.pImageX[iatom[3]];
          ykl = cSim.pImageY[iatom[2]] - cSim.pImageY[iatom[3]];
          zkl = cSim.pImageZ[iatom[2]] - cSim.pImageZ[iatom[3]];
        }
      }

      if (doDistance || doAngle) {
        rij2 = __image_dist2<PMEDouble>(xij, yij, zij, recip, ucell);
        rij = sqrt(rij2);
        if (doDistance) {
          ff = rij;
        } else if (doAngle) {
          rkj2 = xkj * xkj + ykj * ykj + zkj * zkj;
          rik_i = Rsqrt(rij2 * rkj2);
          cst = min(plimit, max(nlimit, (xij*xkj + yij * ykj + zij * zkj) * rik_i));
          ff = acos(cst);
        }
      } else if (doDihedral) {
        // Calculate ij X jk AND kl X jk:
        dx = yij * zkj - zij * ykj;
        dy = zij * xkj - xij * zkj;
        dz = xij * ykj - yij * xkj;
        gx = zkj * ykl - ykj * zkl;
        gy = xkj * zkl - zkj * xkl;
        gz = ykj * xkl - xkj * ykl;

        // Calculate the magnitudes of above vectors, and their dot product:
        fxi = dx * dx + dy * dy + dz * dz + tm24;
        fyi = gx * gx + gy * gy + gz * gz + tm24;
        z1 = Rsqrt(fxi);
        z2 = Rsqrt(fyi);
        z11 = z1 * z1;
        z22 = z2 * z2;
        z12 = z1 * z2;

        cst = min(plimit, max(nlimit, (dx * gx + dy * gy + dz * gz) * z1 * z2));
        ff = acos(cst);
        PMEDouble s = xkj * (dz*gy - dy * gz) + ykj * (dx*gz - dz * gx) + zkj * (dy*gx - dx * gy);
        if (s < (PMEDouble)0.0) ff *= -1.0;
        ff = PI - ff;
        faster_sincos(ff, &sphi, &cphi);

        // Translate the value of the torsion (by +- n*360) to bring it as close as
        // possible to (r2+r3)/2
        PMEDouble apmean = (R.y + R.z) * (PMEDouble)0.5;
        if (ff - apmean > PI) {
          ff -= (PMEDouble)2.0 * (PMEDouble)PI;
        }
        if (apmean - ff > PI) {
          ff += (PMEDouble)2.0 * (PMEDouble)PI;
        }
      }

      PMEDouble r0 = ZeroF, k0 = ZeroF, dr = ZeroF, dK = ZeroF, dr1= ZeroF, dudl= ZeroF;

      if (ff <= R.x) {  // < R1
        k0 = K.x;
        da = (R.y - R.x);
        df = -TwoF * k0 * da;

        if (needPotEnergy) {
          tt = k0 * da * (TwoF * (R.x - ff) + da);

          dK = cSim.TIItemdWeight[Schedule::TypeRestBA][0] * K0.x + cSim.TIItemdWeight[Schedule::TypeRestBA][1] * K1.x;
          dr = cSim.TIItemdWeight[Schedule::TypeRestBA][0] * (R0.y - R0.x) + cSim.TIItemdWeight[Schedule::TypeRestBA][1] * (R1.y - R1.x);
          dr1 = cSim.TIItemdWeight[Schedule::TypeRestBA][0] * R0.x + cSim.TIItemdWeight[Schedule::TypeRestBA][1] * R1.x;
          dudl = dK * da * (TwoF * (R.x - ff) + da) +
            k0 * (da * (dr - TwoF*dr1) + dr * (TwoF * (R.x - ff) + da));
        }

      } else if (ff <= R.y) {  // between R1 and R2
        k0 = K.x;
        da = ff - R.y;
        df = TwoF * k0 * da;

        if (needPotEnergy) {
          tt = k0 * da * da;
          dr = cSim.TIItemdWeight[Schedule::TypeRestBA][0] * R0.y + cSim.TIItemdWeight[Schedule::TypeRestBA][1] * R1.y;
          dK = cSim.TIItemdWeight[Schedule::TypeRestBA][0] * K0.x + cSim.TIItemdWeight[Schedule::TypeRestBA][1] * K1.x;
          dudl = dK * (da * da) + df * dr;
        }

      } else if (ff <= R.z) { // between R2 and R3 : flat region
         return;

      } else if (ff <= R.w) { // between R3 and R4
        k0 = K.y;
        da = ff - R.z;
        df = TwoF * k0 * da;

        if (needPotEnergy) {
          tt = k0 * da * da;
          dr = cSim.TIItemdWeight[Schedule::TypeRestBA][0] * R0.z + cSim.TIItemdWeight[Schedule::TypeRestBA][1] * R1.z;
          dK = cSim.TIItemdWeight[Schedule::TypeRestBA][0] * K0.y + cSim.TIItemdWeight[Schedule::TypeRestBA][1] * K1.y;
          dudl = dK * (da * da) + df * dr;
        }

      } else {  // >R4
        k0 = K.y;
        da = (R.w - R.z);
        df = TwoF * k0 * da;

        if (needPotEnergy) {
          tt = k0 * da * (TwoF * (ff - R.w) + da);

          dK = cSim.TIItemdWeight[Schedule::TypeRestBA][0] * K0.y + cSim.TIItemdWeight[Schedule::TypeRestBA][1] * K1.y;
          dr = cSim.TIItemdWeight[Schedule::TypeRestBA][0] * (R0.w - R0.z) + cSim.TIItemdWeight[Schedule::TypeRestBA][1] * (R1.w - R1.z);
          dr1 = cSim.TIItemdWeight[Schedule::TypeRestBA][0] * R0.w + cSim.TIItemdWeight[Schedule::TypeRestBA][1] * R1.w;
          dudl = dK * da * (TwoF * (ff - R.w) + da) +
            k0 * (da * (dr - TwoF*dr1) + dr * (TwoF * (ff - R.x) + da));
        }
      }

      //Restraint Energy
      if (needPotEnergy) {
        addEnergy(pTIRestEnergy, tt);

        // dU/dl
        addEnergy(pTIRestDer, dudl);

        // MBAR part--go through all listed lambda values
        if (cSim.needMBAR) {
          PMEFloat* weight[2] = { cSim.pMBARWeight , &(cSim.pMBARWeight[cSim.nMBARStates]) };
          unsigned mbarShift = cSim.nMBARStates * (2 + Schedule::TypeRestBA * 3); // Type BAT; put into the region-independent part
          unsigned long long int* MBAREnergy= &(cSim.pMBAREnergy[mbarShift]);
	  unsigned i = 0;
          unsigned j = cSim.nMBARStates-1;
          if (cSim.currentMBARState >= 0 && cSim.rangeMBAR >= 0) {
            i = max(0, cSim.currentMBARState - cSim.rangeMBAR);
            j = min(cSim.nMBARStates - 1, cSim.currentMBARState + cSim.rangeMBAR);
          }
          for (unsigned l = i; l <= j; l++) {
            PMEDouble2 tempK = {
              weight[0][l] * K0.x + weight[1][l] * K1.x ,
              weight[0][l] * K0.y + weight[1][l] * K1.y };
            PMEDouble4 tempR = {
              weight[0][l] * R0.x + weight[1][l] * R1.x ,
              weight[0][l] * R0.y + weight[1][l] * R1.y ,
              weight[0][l] * R0.z + weight[1][l] * R1.z , 
              weight[0][l] * R0.w + weight[1][l] * R1.w  };

            PMEDouble tt = ZeroF;
            if (ff <= tempR.x) {
              da = tempR.y- tempR.x;
              k0 = tempK.x;
              tt = (tempR.x - ff) * TwoF * k0* da + k0 * da * da;
            } else if (ff <= tempR.y) {
              da = ff - tempR.y;
              k0 = tempK.x;
              tt = k0 * da* da;
            } else if (ff <= tempR.z) {
              continue;
            } else if (ff <= tempR.w) {
              da = ff - tempR.z;
              k0 = tempK.y;
              tt = k0 * da * da;
            } else {
              da = tempR.w - tempR.z;
              k0 = tempK.y;
              tt = (ff - tempR.w) * TwoF * k0 * da + k0 * da * da;
            }
	    addEnergy(&(MBAREnergy[l]), tt);
          }
        }  // if (cSim.needMBAR) {
      }

      // Add to the force accumulators
      if (doDistance) {
        dfw = (df + df) / ff;
        tt = dfw * xij;
        addForce(pTIFx[2] + iatom[0], -tt);
        addForce(pTIFx[2] + iatom[1], tt);

        tt = dfw * yij;
        addForce(pTIFy[2] + iatom[0], -tt);
        addForce(pTIFy[2] + iatom[1], tt);

        tt = dfw * zij;
        addForce(pTIFz[2] + iatom[0], -tt);
        addForce(pTIFz[2] + iatom[1], tt);

      } else if (doAngle) {
        dfw = -(df + df) / sin(ff);

        cik = dfw * rik_i;
        sth = dfw * cst;
        cii = sth / rij2;
        ckk = sth / rkj2;
        fx1 = cii * xij - cik * xkj;
        fy1 = cii * yij - cik * ykj;
        fz1 = cii * zij - cik * zkj;
        fx2 = ckk * xkj - cik * xij;
        fy2 = ckk * ykj - cik * yij;
        fz2 = ckk * zkj - cik * zij;

        addForce(pTIFx[2] + iatom[0], fx1);
        addForce(pTIFy[2] + iatom[0], fy1);
        addForce(pTIFz[2] + iatom[0], fz1);

        addForce(pTIFx[2] + iatom[2], fx2);
        addForce(pTIFy[2] + iatom[2], fy2);
        addForce(pTIFz[2] + iatom[2], fz2);

        addForce(pTIFx[2] + iatom[1], -fx1 - fx2);
        addForce(pTIFy[2] + iatom[1], -fy1 - fy2);
        addForce(pTIFz[2] + iatom[1], -fz1 - fz2);

      } else if (doDihedral) {

        dfw = -(df + df) / sphi;

        dcdx = -gx * z12 - cphi * dx*z11;
        dcdy = -gy * z12 - cphi * dy*z11;
        dcdz = -gz * z12 - cphi * dz*z11;
        dcgx = dx * z12 + cphi * gx*z22;
        dcgy = dy * z12 + cphi * gy*z22;
        dcgz = dz * z12 + cphi * gz*z22;

        dr1 = dfw * (dcdz*ykj - dcdy * zkj);
        dr2 = dfw * (dcdx*zkj - dcdz * xkj);
        dr3 = dfw * (dcdy*xkj - dcdx * ykj);
        dr4 = dfw * (dcgz*ykj - dcgy * zkj);
        dr5 = dfw * (dcgx*zkj - dcgz * xkj);
        dr6 = dfw * (dcgy*xkj - dcgx * ykj);
        drx = dfw * (-dcdy * zij + dcdz * yij + dcgy * zkl - dcgz * ykl);
        dry = dfw * (dcdx*zij - dcdz * xij - dcgx * zkl + dcgz * xkl);
        drz = dfw * (-dcdx * yij + dcdy * xij + dcgx * ykl - dcgy * xkl);

        addForce(pTIFx[2] + iatom[0], (-dr1));
        addForce(pTIFy[2] + iatom[0], (-dr2));
        addForce(pTIFz[2] + iatom[0], (-dr3));
        addForce(pTIFx[2] + iatom[1], (-drx + dr1));
        addForce(pTIFy[2] + iatom[1], (-dry + dr2));
        addForce(pTIFz[2] + iatom[1], (-drz + dr3));
        addForce(pTIFx[2] + iatom[2], (drx + dr4));
        addForce(pTIFy[2] + iatom[2], (dry + dr5));
        addForce(pTIFz[2] + iatom[2], (drz + dr6));
        addForce(pTIFx[2] + iatom[3], (-dr4));
        addForce(pTIFy[2] + iatom[3], (-dr5));
        addForce(pTIFz[2] + iatom[3], (-dr6));
      }
   }
}


//---------------------------------------------------------------------------------------------
// kgREAFBondedEnergy_kernel:
//
// Arguments:
//   needPotEnergy:
//   needvirial:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgREAFBondedEnergy_kernel(bool needPotEnergy, bool needvirial) {
  int pos = blockIdx.x * blockDim.x + threadIdx.x;

  int nwarp_angle = ((cSim.numberREAFAngle - 1) / 32 + 1) * 32;
  int nwarp_dihedral = ((cSim.numberREAFDihedral - 1) / 32 + 1) * 32;

  int angleBegin = 0;
  int angleEnd = cSim.numberREAFAngle - 1;

  int dihedralBegin = angleBegin + nwarp_angle;
  int dihedralEnd = dihedralBegin + cSim.numberREAFDihedral - 1;


  unsigned iatom[4], READRegion;
  PMEDouble xij, xkj, xkl, yij, ykj, ykl, zij, zkj, zkl;
  PMEDouble ff, da, df, dfw, tt, tw;
  PMEDouble rij2, rij, rkj2, rik_i, cst;
  PMEDouble dx, dy, dz, gx, gy, gz;
  PMEDouble dcdx, dcdy, dcdz, dcgx, dcgy, dcgz;
  PMEDouble fxi, fyi, z1, z2, z11, z12, z22;
  PMEDouble sphi, cphi;
  PMEDouble dr1, dr2, dr3, dr4, dr5, dr6, drx, dry, drz;
  PMEDouble cik, sth, cii, ckk;
  PMEDouble fx1, fy1, fz1, fx2, fy2, fz2;

  PMEAccumulator* pFx= cSim.pNBForceXAccumulator, * pFy= cSim.pNBForceYAccumulator, * pFz= cSim.pNBForceZAccumulator;

  if (pos >= angleBegin && pos <= angleEnd) {
    unsigned mypos = pos - angleBegin;
    uint4& angleID = cSim.pREAFBondAngleID[mypos];
    PMEDouble2& angle = cSim.pREAFBondAngle[mypos];
    iatom[0] = cSim.pImageAtomLookup[angleID.x];
    iatom[1] = cSim.pImageAtomLookup[angleID.y];
    iatom[2] = cSim.pImageAtomLookup[angleID.z];

    bool hasRE = (cSim.pREAFBondAngleType[mypos] & gti_simulationConst::Fg_has_RE);
    bool inRE = (cSim.pREAFBondAngleType[mypos] & gti_simulationConst::Fg_int_RE);
    bool addRE = (cSim.pREAFBondAngleType[mypos] & gti_simulationConst::ex_addRE_bat);

    PMEDouble tw = addRE? (cSim.REAFItemWeight[Schedule::TypeBAT][inRE ? 1 : 0]) - OneF : ZeroF;
    if (abs(tw) > tm06) {

      xij = cSim.pImageX[iatom[0]] - cSim.pImageX[iatom[1]];
      yij = cSim.pImageY[iatom[0]] - cSim.pImageY[iatom[1]];
      zij = cSim.pImageZ[iatom[0]] - cSim.pImageZ[iatom[1]];
      xkj = cSim.pImageX[iatom[2]] - cSim.pImageX[iatom[1]];
      ykj = cSim.pImageY[iatom[2]] - cSim.pImageY[iatom[1]];
      zkj = cSim.pImageZ[iatom[2]] - cSim.pImageZ[iatom[1]];
      rij = xij * xij + yij * yij + zij * zij;
      PMEDouble rkj = xkj * xkj + ykj * ykj + zkj * zkj;
      PMEDouble rik = sqrt(rij * rkj);

      const PMEDouble nlimit = -0.99999;
      const PMEDouble plimit = 0.99999;
      PMEDouble cst = min(plimit, max(nlimit, (xij * xkj + yij * ykj + zij * zkj) / rik));
      PMEDouble ant = acos(cst);

      da = ant - angle.y;
      df = angle.x * da;
      dfw = -(df + df) / sin(ant);

      if (needPotEnergy) {
          tt = df * da;
          addEnergy(cSim.pEAngle, tt*tw);
      }

      // Calculation of the force
      cik = dfw / rik;
      sth = dfw * cst;
      cii = sth / rij;
      ckk = sth / rkj;
      fx1 = cii * xij - cik * xkj;
      fy1 = cii * yij - cik * ykj;
      fz1 = cii * zij - cik * zkj;
      fx2 = ckk * xkj - cik * xij;
      fy2 = ckk * ykj - cik * yij;
      fz2 = ckk * zkj - cik * zij;

      if (fabs(tw) > tm06) {
        addForce(pFx + iatom[0], fx1 * tw);
        addForce(pFy + iatom[0], fy1 * tw);
        addForce(pFz + iatom[0], fz1 * tw);

        addForce(pFx + iatom[2], fx2 * tw);
        addForce(pFy + iatom[2], fy2 * tw);
        addForce(pFz + iatom[2], fz2 * tw);

        addForce(pFx + iatom[1], -tw * (fx1 + fx2));
        addForce(pFy + iatom[1], -tw * (fy1 + fy2));
        addForce(pFz + iatom[1], -tw * (fz1 + fz2));
      }

    }
  }
  if (pos >= dihedralBegin && pos <= dihedralEnd) {
    unsigned mypos = pos - dihedralBegin;
    uint4& dihedralID = cSim.pREAFDihedralID[mypos];
    iatom[0] = cSim.pImageAtomLookup[dihedralID.x];
    iatom[1] = cSim.pImageAtomLookup[dihedralID.y];
    iatom[2] = cSim.pImageAtomLookup[dihedralID.z];
    iatom[3] = cSim.pImageAtomLookup[dihedralID.w];

    PMEDouble2& dihedral1 = cSim.pREAFDihedral1[mypos];
    PMEDouble2& dihedral2 = cSim.pREAFDihedral2[mypos];
    PMEDouble& dihedral3 = cSim.pREAFDihedral3[mypos];

    bool hasRE = (cSim.pREAFDihedralType[mypos] & gti_simulationConst::Fg_has_RE);
    bool inRE = (cSim.pREAFDihedralType[mypos] & gti_simulationConst::Fg_int_RE);
    bool addRE = (cSim.pREAFDihedralType[mypos] & gti_simulationConst::ex_addRE_bat);

    PMEDouble tw = addRE ? (cSim.REAFItemWeight[Schedule::TypeBAT][inRE ? 1 : 0]) - OneF: Zero;

    if (abs(tw) > tm06) {
      xij = cSim.pImageX[iatom[0]] - cSim.pImageX[iatom[1]];
      yij = cSim.pImageY[iatom[0]] - cSim.pImageY[iatom[1]];
      zij = cSim.pImageZ[iatom[0]] - cSim.pImageZ[iatom[1]];
      xkj = cSim.pImageX[iatom[2]] - cSim.pImageX[iatom[1]];
      ykj = cSim.pImageY[iatom[2]] - cSim.pImageY[iatom[1]];
      zkj = cSim.pImageZ[iatom[2]] - cSim.pImageZ[iatom[1]];
      xkl = cSim.pImageX[iatom[2]] - cSim.pImageX[iatom[3]];
      ykl = cSim.pImageY[iatom[2]] - cSim.pImageY[iatom[3]];
      zkl = cSim.pImageZ[iatom[2]] - cSim.pImageZ[iatom[3]];

      // Get the normal vector
      dx = yij * zkj - zij * ykj;
      dy = zij * xkj - xij * zkj;
      dz = xij * ykj - yij * xkj;
      gx = zkj * ykl - ykj * zkl;
      gy = xkj * zkl - zkj * xkl;
      gz = ykj * xkl - xkj * ykl;
      fxi = sqrt(dx * dx + dy * dy + dz * dz + tm24);
      fyi = sqrt(gx * gx + gy * gy + gz * gz + tm24);
      cst = dx * gx + dy * gy + dz * gz;

      // Branch if linear dihedral:
      PMEDouble z1 = (tenm3 <= fxi) ? (1.0 / fxi) : zero;
      PMEDouble z2 = (tenm3 <= fyi) ? (1.0 / fyi) : zero;
      PMEDouble z12 = z1 * z2;
      PMEDouble fzi = (z12 != zero) ? one : zero;
      PMEDouble s = xkj * (dz * gy - dy * gz) + ykj * (dx * gz - dz * gx) + zkj * (dy * gx - dx * gy);
      PMEDouble ff = PI - abs(acos(max(-1.0, min(1.0, cst * z12)))) * (s >= zero ? one : -one);
      faster_sincos(ff, &sphi, &cphi);

      // Calculate the energy and the derivatives with respect to cosphi
      PMEDouble ct0 = dihedral1.y * ff;
      PMEDouble sinnp, cosnp;
      faster_sincos(ct0, &sinnp, &cosnp);

      if (needPotEnergy) {
        tt = (dihedral2.x + (cosnp * dihedral2.y) + (sinnp * dihedral3)) * fzi;
        addEnergy(cSim.pEDihedral, tt*tw);
      }


#ifdef use_DPFP
      PMEFloat delta = tm24;
#else
      PMEFloat delta = tm06;
#endif
      PMEDouble dums = sphi + (delta * (sphi >= zero ? one : -one));
      PMEDouble df;
      if (tm06 > abs(dums)) {
        df = fzi * dihedral2.y * (dihedral1.y - dihedral1.x + (dihedral1.x * cphi));
      } else {
        df = fzi * dihedral1.y * ((dihedral2.y * sinnp) - (dihedral3 * cosnp)) / dums;
      }

      // Now, set up array dc = 1st der. of cosphi w/respect to cartesian differences:
      z11 = z1 * z1;
      z12 = z1 * z2;
      z22 = z2 * z2;
      dcdx = -gx * z12 - cphi * dx * z11;
      dcdy = -gy * z12 - cphi * dy * z11;
      dcdz = -gz * z12 - cphi * dz * z11;
      dcgx = dx * z12 + cphi * gx * z22;
      dcgy = dy * z12 + cphi * gy * z22;
      dcgz = dz * z12 + cphi * gz * z22;

      // Update the first derivative array:
      dr1 = df * (dcdz * ykj - dcdy * zkj);
      dr2 = df * (dcdx * zkj - dcdz * xkj);
      dr3 = df * (dcdy * xkj - dcdx * ykj);
      dr4 = df * (dcgz * ykj - dcgy * zkj);
      dr5 = df * (dcgx * zkj - dcgz * xkj);
      dr6 = df * (dcgy * xkj - dcgx * ykj);
      drx = df * (-dcdy * zij + dcdz * yij + dcgy * zkl - dcgz * ykl);
      dry = df * (dcdx * zij - dcdz * xij - dcgx * zkl + dcgz * xkl);
      drz = df * (-dcdx * yij + dcdy * xij + dcgx * ykl - dcgy * xkl);

      addForce(pFx + iatom[0], (-dr1 * tw));
      addForce(pFy + iatom[0], (-dr2 * tw));
      addForce(pFz + iatom[0], (-dr3 * tw));

      addForce(pFx + iatom[1], (-drx + dr1) * tw);
      addForce(pFy + iatom[1], (-dry + dr2) * tw);
      addForce(pFz + iatom[1], (-drz + dr3) * tw);

      addForce(pFx + iatom[2], (drx + dr4) * tw);
      addForce(pFy + iatom[2], (dry + dr5) * tw);
      addForce(pFz + iatom[2], (drz + dr6) * tw);

      addForce(pFx + iatom[3], (-dr4 * tw));
      addForce(pFy + iatom[3], (-dr5 * tw));
      addForce(pFz + iatom[3], (-dr6 * tw));
    }

  }

}

#ifndef AMBER_PLATFORM_AMD 
//---------------------------------------------------------------------------------------------
// kgRMSDPreparation_kernel:
//
// Arguments:
//   needPotEnergy:
//   needvirial:
//---------------------------------------------------------------------------------------------
__device__ PMEFloat com[5][3];
__device__ PMEFloat F[5][4][4];

_kPlainHead_ kgRMSDPreparation_kernel(bool needPotEnergy) {

  int mySet = blockIdx.x;
  int pos = threadIdx.x;
  if (cSim.rmsd_atom_count[mySet] == 0) return;

  if(pos<3)  com[mySet][pos] = 0;  
  __syncthreads();

  unsigned& max_count = cSim.rmsd_atom_max_count;

  // Calculate the center of geometry
  //
  PMEFloat myCom[3] = { 0.0, 0.0, 0.0 };
  while (pos < cSim.rmsd_atom_count[mySet]) {
      unsigned shift = mySet * max_count + pos;
      unsigned index = cSim.rmsd_atom_list[shift];
      unsigned iatom = cSim.pImageAtomLookup[index-1];
      myCom[0] += cSim.pImageX[iatom];
      myCom[1] += cSim.pImageY[iatom];
      myCom[2] += cSim.pImageZ[iatom];
      //TBT
      //printf("XXXXX %4d %4d %4d %10.6f %10.6f %10.6f \n", mySet
      //, index, iatom, cSim.pImageX[iatom], cSim.pImageY[iatom], cSim.pImageZ[iatom]);
      pos += blockDim.x;
  }
  __syncthreads();

#pragma unroll
  for (int offset = warpSize / 2; offset > 0; offset /= 2) {
      __syncwarp();
      myCom[0] += __SHFL_DOWN(0xFFFFFFFF, myCom[0], offset);
      myCom[1] += __SHFL_DOWN(0xFFFFFFFF, myCom[1], offset);
      myCom[2] += __SHFL_DOWN(0xFFFFFFFF, myCom[2], offset);
  }

  __syncthreads();

  pos = threadIdx.x;
  if (pos == 0) {
    com[mySet][0] = myCom[0] / cSim.rmsd_atom_count[mySet];
    com[mySet][1] = myCom[1] / cSim.rmsd_atom_count[mySet];
    com[mySet][2] = myCom[2] / cSim.rmsd_atom_count[mySet];
  }
  __syncthreads();

  // Calculate the R matrix
  //
  pos = threadIdx.x;
  PMEFloat R[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},{0.0, 0.0, 0.0} };
  while (pos < cSim.rmsd_atom_count[mySet]) {
      unsigned shift = mySet * max_count + pos;
      unsigned index = cSim.rmsd_atom_list[shift];
      unsigned iatom = cSim.pImageAtomLookup[index - 1];
      PMEFloat diffX[3] = {
          cSim.pImageX[iatom] - com[mySet][0] ,
          cSim.pImageY[iatom] - com[mySet][1] ,
          cSim.pImageZ[iatom] - com[mySet][2] };

      shift *= 3;
      for (unsigned i = 0; i < 3; i++)
        for (unsigned j = 0; j < 3; j++)
          R[i][j] += diffX[i] * cSim.rmsd_ref_crd[shift+j];
          
      pos += blockDim.x;
  }
  __syncthreads();

#pragma unroll
  for (int offset = warpSize / 2; offset > 0; offset /= 2) {
      __syncwarp();
      for (unsigned i = 0; i < 3; i++)
        for (unsigned j = 0; j < 3; j++)
          R[i][j] += __SHFL_DOWN(0xFFFFFFFF, R[i][j], offset);
  }

  __syncthreads();

  if (threadIdx.x == 0) {
      unsigned shift = mySet * 16;
      cSim.pMatrix[shift] = R[0][0] + R[1][1] + R[2][2];
      cSim.pMatrix[shift+1] = R[1][2] - R[2][1];
      cSim.pMatrix[shift+2] = R[2][0] - R[0][2];
      cSim.pMatrix[shift+3] = R[0][1] - R[1][0];

      cSim.pMatrix[shift+4] = R[1][2] - R[2][1];
      cSim.pMatrix[shift+5] = R[0][0] - R[1][1] - R[2][2];
      cSim.pMatrix[shift+6] = R[0][1] + R[1][0];
      cSim.pMatrix[shift+7] = R[0][2] + R[2][0];

      cSim.pMatrix[shift+8] = R[2][0] - R[0][2];
      cSim.pMatrix[shift+9] = R[0][1] + R[1][0];
      cSim.pMatrix[shift+10] = -R[0][0] + R[1][1] - R[2][2];
      cSim.pMatrix[shift+11] = R[1][2] + R[2][1];

      cSim.pMatrix[shift+12] = R[0][1] - R[1][0];
      cSim.pMatrix[shift+13] = R[0][2] + R[2][0];
      cSim.pMatrix[shift+14] = R[1][2] + R[2][1];
      cSim.pMatrix[shift+15] = -R[0][0] - R[1][1] + R[2][2];
  }
}

//---------------------------------------------------------------------------------------------
// kgRMSDEnergyForce_kernel:
//
// Arguments:
//   needPotEnergy:
//   needvirial:
//---------------------------------------------------------------------------------------------

_kPlainHead_ kgRMSDEnergyForce_kernel(bool needPotEnergy) {

  int mySet = blockIdx.x;
  int pos = threadIdx.x;
  if (cSim.rmsd_atom_count[mySet] == 0) return;

  PMEAccumulator* pTIFx[3] = { cSim.pTIForce,
                      cSim.pTIForce + cSim.stride3,
                      cSim.pTIForce + (cSim.stride3 * 2) };
  PMEAccumulator* pTIFy[3] = { cSim.pTIForce + cSim.stride,
                        cSim.pTIForce + cSim.stride + cSim.stride3,
                        cSim.pTIForce + cSim.stride + (cSim.stride3 * 2) };
  PMEAccumulator* pTIFz[3] = { cSim.pTIForce + (cSim.stride * 2),
                        cSim.pTIForce + (cSim.stride * 2) + cSim.stride3,
                        cSim.pTIForce + (cSim.stride * 2) +
                        (cSim.stride3 * 2) };

  unsigned& energyType = cSim.rmsd_type; //0: RMSD, 1:SD

  __shared__ PMEFloat U[3][3], rotatedRef[3], RMSD, RMSD2;
  PMEFloat mySD = 0.0, mySD2=0.0;

  unsigned shift = mySet * 4;
  PMEFloat eig[4] = { cSim.pResult[shift],  cSim.pResult[shift + 1], cSim.pResult[shift + 2], cSim.pResult[shift + 3] };

  int info = *(cSim.pInfo);
  unsigned& max_count = cSim.rmsd_atom_max_count;

  // Compute the rotation matrix.
  if (pos == 0) {
      unsigned shift = mySet * 16;
      PMEFloat q[] = {
         cSim.pMatrix[shift + 12],
         cSim.pMatrix[shift + 13],
         cSim.pMatrix[shift + 14],
         cSim.pMatrix[shift + 15] };
      PMEFloat q00 = q[0] * q[0], q01 = q[0] * q[1], q02 = q[0] * q[2], q03 = q[0] * q[3];
      PMEFloat q11 = q[1] * q[1], q12 = q[1] * q[2], q13 = q[1] * q[3];
      PMEFloat q22 = q[2] * q[2], q23 = q[2] * q[3];
      PMEFloat q33 = q[3] * q[3];

      U[0][0] = q00 + q11 - q22 - q33;
      U[0][1] = 2 * (q12 - q03);
      U[0][2] = 2 * (q13 + q02);
      U[1][0] = 2 * (q12 + q03);
      U[1][1] = q00 - q11 + q22 - q33 ;
      U[1][2] = 2 * (q23 - q01);
      U[2][0] = 2 * (q13 - q02); 
      U[2][1] = 2 * (q23 + q01);
      U[2][2] = q00 - q11 - q22 + q33;
  }

  __syncthreads();

  // Rotate the reference coordinates and compute forces
  //
  unsigned myTI = (cSim.rmsd_ti_region[mySet] > 0 && cSim.rmsd_ti_region[mySet] < 3) ? cSim.rmsd_ti_region[mySet] - 1 : 2;

  PMEFloat tiWeight= (myTI == 2) ? 1.0 : cSim.TIItemWeight[Schedule::TypeRMSD][cSim.rmsd_ti_region[mySet]-1];
  PMEFloat dlWeight = (myTI == 2) ? 1.0 : cSim.TIItemdWeight[Schedule::TypeRMSD][cSim.rmsd_ti_region[mySet] - 1];

  while (pos < cSim.rmsd_atom_count[mySet]) {
    unsigned shift = mySet * max_count + pos;  // atom index
    unsigned index = cSim.rmsd_atom_list[shift];
    unsigned iatom = cSim.pImageAtomLookup[index - 1];
    
    shift *= 3; // change to crd atom index
    PMEFloat p[3] = { 
      cSim.rmsd_ref_crd[shift], 
      cSim.rmsd_ref_crd[shift+1], 
      cSim.rmsd_ref_crd[shift+2] 
    };
    PMEFloat rotatedRef[3] = {
      U[0][0] * p[0] + U[1][0] * p[1] + U[2][0] * p[2],
      U[0][1] * p[0] + U[1][1] * p[1] + U[2][1] * p[2],
      U[0][2] * p[0] + U[1][2] * p[1] + U[2][2] * p[2]
    };

    PMEFloat diffX[3] = {
      cSim.pImageX[iatom] - com[mySet][0] ,
      cSim.pImageY[iatom] - com[mySet][1] ,
      cSim.pImageZ[iatom] - com[mySet][2] 
    };

    mySD +=
      (rotatedRef[0] - diffX[0]) * (rotatedRef[0] - diffX[0])
      + (rotatedRef[1] - diffX[1]) * (rotatedRef[1] - diffX[1])
      + (rotatedRef[2] - diffX[2]) * (rotatedRef[2] - diffX[2]);

    mySD2 += (p[0] * p[0]
      + p[1] * p[1]
      + p[2] * p[2]
      + diffX[0] * diffX[0]
      + diffX[1] * diffX[1]
      + diffX[2] * diffX[2] 
      );
     
    pos += blockDim.x;
  }

  __syncthreads();

#pragma unroll

  for (int offset = warpSize / 2; offset > 0; offset /= 2) {
    __syncwarp();
    mySD += __SHFL_DOWN(0xFFFFFFFF, mySD, offset);
    mySD2 += __SHFL_DOWN(0xFFFFFFFF, mySD2, offset);
  }

  __syncthreads();

  if (threadIdx.x == 0) {
    RMSD = Sqrt(mySD / cSim.rmsd_atom_count[mySet]);
    RMSD2 = Sqrt((mySD2 - 2 * cSim.pResult[mySet * 4 + 3]) / cSim.rmsd_atom_count[mySet]);
  }

  __syncthreads();

  pos = threadIdx.x;
  PMEFloat forcefactor = (energyType == 0)
    ? 1.0 /(RMSD *cSim.rmsd_atom_count[mySet])
    : 2.0;

  forcefactor *= (tiWeight* cSim.rmsd_weights[mySet]);

  //PMEFloat factor = 2*tiWeight * cSim.rmsd_weights[mySet];
  while (pos < cSim.rmsd_atom_count[mySet]) {
    unsigned shift = mySet * max_count + pos;
    unsigned index = cSim.rmsd_atom_list[shift];
    unsigned iatom = cSim.pImageAtomLookup[index - 1];

    shift *= 3;
    PMEFloat p[3] = {
      cSim.rmsd_ref_crd[shift],
      cSim.rmsd_ref_crd[shift + 1],
      cSim.rmsd_ref_crd[shift + 2]
    };
    PMEFloat rotatedRef[3] = {
      U[0][0] * p[0] + U[1][0] * p[1] + U[2][0] * p[2],
      U[0][1] * p[0] + U[1][1] * p[1] + U[2][1] * p[2],
      U[0][2] * p[0] + U[1][2] * p[1] + U[2][2] * p[2]
    };

    PMEFloat diffX[3] = {
      cSim.pImageX[iatom] - com[mySet][0] ,
      cSim.pImageY[iatom] - com[mySet][1] ,
      cSim.pImageZ[iatom] - com[mySet][2]
    };

    PMEFloat force[3] = {
      force[0] = forcefactor * (rotatedRef[0] - diffX[0]) ,
      force[1] = forcefactor * (rotatedRef[1] - diffX[1]) ,
      force[2] = forcefactor * (rotatedRef[2] - diffX[2])
    };

    addForce(pTIFx[2] + iatom, force[0]);
    addForce(pTIFy[2] + iatom, force[1]);
    addForce(pTIFz[2] + iatom, force[2]);

    pos += blockDim.x;
  }

  if (threadIdx.x == 0) {

    static int counter = 0;
    PMEFloat tt = (energyType == 0) 
      ? RMSD : mySD;

    tt *= cSim.rmsd_weights[mySet];

    addEnergy(cSim.pTIRestRMSD[2], tt* tiWeight);
    if (myTI != 2) {
        addEnergy(cSim.pTIRestRMSD_DL[2], tt * dlWeight);

        if (cSim.needMBAR) {
          PMEFloat* lambda[2] = { cSim.pMBARLambda , &(cSim.pMBARLambda[cSim.nMBARStates]) };
          unsigned mbarShift = cSim.nMBARStates * (2 + Schedule::TypeRMSD * 3); // Type BAT; put into the region-independent part
          unsigned long long int* MBAREnergy = &(cSim.pMBAREnergy[mbarShift]);
          unsigned i = 0;
          unsigned j = cSim.nMBARStates - 1;
          if (cSim.currentMBARState >= 0 && cSim.rangeMBAR >= 0) {
            i = max(Zero, cSim.currentMBARState - cSim.rangeMBAR);
            j = min(cSim.nMBARStates - 1, cSim.currentMBARState + cSim.rangeMBAR);
          }
          for (unsigned l = i; l <= j; l++) {
            addEnergy(&(MBAREnergy[l]), tt *(OneF- lambda[myTI][l]));
          }
        }  // if (cSim.needMBAR) 
    }
      
    //  for debug purpose--printout RMSD values
    /*
     if (counter++% 1
       0 == 0) {
        printf("RMSD Set # %1d  Iteration # %7d %12.8f %12.8f  %12.8f  %12.8f  %12.8f \n"
          , mySet+1
          , counter
          , RMSD
          , RMSD2
          , forcefactor
          , ((energyType == 0) ? RMSD : mySD )* cSim.rmsd_weights[mySet] * tiWeight
          , ((energyType == 0) ? RMSD : mySD )* cSim.rmsd_weights[mySet] * dlWeight);
    }
   */ 

  }
}
#endif

//---------------------------------------------------------------------------------------------
// kgBondedEnergy_ppi_kernel:
//
// Arguments:
//   needPotEnergy:
//   needvirial:
//---------------------------------------------------------------------------------------------
_kPlainHead_ kgBondedEnergy_ppi_kernel(bool needPotEnergy, bool needvirial)
{
  int pos = blockIdx.x * blockDim.x + threadIdx.x;
  //unsigned int totalTerms = cSim.numberTIBond + cSim.numberTIAngle + cSim.numberTIDihedral;

  int nwarp_bond = ((cSim.numberTIBond - 1) / GRID + 1) * GRID;
  int nwarp_angle = ((cSim.numberTIAngle - 1) / GRID + 1) * GRID;
  //int nwarp_dihedral = ((cSim.numberTIDihedral - 1) / GRID + 1) * GRID;

  //int nwarp_NMRDistance = ((cSim.numberTINMRDistance - 1) / GRID + 1) * GRID;
  //int nwarp_NMRAngle = ((cSim.numberTINMRAngle - 1) / GRID + 1) * GRID;
  //int nwarp_NMRDihedral = ((cSim.numberTINMRDihedral - 1) / GRID + 1) * GRID;

  int bondBegin = 0;
  int bondEnd = cSim.numberTIBond - 1;

  int angleBegin = nwarp_bond;
  int angleEnd = angleBegin + cSim.numberTIAngle - 1;

  int dihedralBegin = angleBegin + nwarp_angle;
  int dihedralEnd = dihedralBegin + cSim.numberTIDihedral - 1;

  //int NMRDistanceBegin = dihedralBegin + nwarp_dihedral;
  //int NMRDistanceEnd = NMRDistanceBegin + cSim.numberTINMRDistance - 1;

  //int NMRAngleBegin = NMRDistanceBegin + nwarp_NMRDistance;
  //int NMRAngleEnd = NMRAngleBegin + cSim.numberTINMRAngle - 1;

  //int NMRDihedralBegin = NMRAngleBegin + nwarp_NMRAngle;
  //int NMRDihedralEnd = NMRDihedralBegin + cSim.numberTINMRDihedral - 1;


  //unsigned long long int* energyPos[3] = { cSim.pTIPotEnergyBuffer,
  //                                        cSim.pTIPotEnergyBuffer + cSim.GPUPotEnergyTerms,
  //                                        cSim.pTIPotEnergyBuffer +
  //                                        (cSim.GPUPotEnergyTerms * 2) };
  //PMEAccumulator* pTIFx[3] = { cSim.pTIForce,
  //                      cSim.pTIForce + cSim.stride3,
  //                      cSim.pTIForce + (cSim.stride3 * 2) };
  //PMEAccumulator* pTIFy[3] = { cSim.pTIForce + cSim.stride,
  //                      cSim.pTIForce + cSim.stride + cSim.stride3,
  //                      cSim.pTIForce + cSim.stride + (cSim.stride3 * 2) };
  //PMEAccumulator* pTIFz[3] = { cSim.pTIForce + (cSim.stride * 2),
  //                      cSim.pTIForce + (cSim.stride * 2) + cSim.stride3,
  //                      cSim.pTIForce + (cSim.stride * 2) +
  //                      (cSim.stride3 * 2) };

  //PMEAccumulator* pSCFx[2] = { cSim.pTIForce + (cSim.stride3 * 3),
  //                      cSim.pTIForce + (cSim.stride3 * 4) };
  //PMEAccumulator* pSCFy[2] = { cSim.pTIForce + (cSim.stride3 * 3) + cSim.stride,
  //                      cSim.pTIForce + (cSim.stride3 * 4) + cSim.stride };
  //PMEAccumulator* pSCFz[2] = { cSim.pTIForce + (cSim.stride3 * 3) +
  //                     (cSim.stride * 2), cSim.pTIForce +
  //                     (cSim.stride3 * 4) + (cSim.stride * 2) };

  unsigned iatom[4], TIRegion;
  PMEDouble xij, xkj, xkl, yij, ykj, ykl, zij, zkj, zkl;
  PMEDouble da, df, dfw, tt, tw;
  PMEDouble rij, cst;
  PMEDouble dx, dy, dz, gx, gy, gz;
  PMEDouble dcdx, dcdy, dcdz, dcgx, dcgy, dcgz;
  PMEDouble fxi, fyi, z11, z22;
  PMEDouble sphi, cphi;
  PMEDouble dr1, dr2, dr3, dr4, dr5, dr6, drx, dry, drz;
  PMEDouble cik, sth, cii, ckk;
  PMEDouble fx1, fy1, fz1, fx2, fy2, fz2;

  unsigned long long int *pE;

  if (pos >= bondBegin && pos <= bondEnd) {
    unsigned mypos = pos - bondBegin;
    uint2& bondID = cSim.pTIBondID[mypos];
    PMEDouble2& bond = cSim.pTIBond[mypos];
    iatom[0] = cSim.pImageAtomLookup[bondID.x];
    iatom[1] = cSim.pImageAtomLookup[bondID.y];

    TIRegion = (cSim.pTIBondType[mypos] & gti_simulationConst::Fg_TI_region);
    //bool hasSC = (cSim.pTIBondType[mypos] & gti_simulationConst::Fg_has_SC);
    //bool inSC = (cSim.pTIBondType[mypos] & gti_simulationConst::Fg_int_SC);
//    bool addToDVDL = (cSim.pTIBondType[mypos] & gti_simulationConst::ex_addToDVDL_bat);
//    bool addToBATCorr = (cSim.pTIBondType[mypos] & gti_simulationConst::ex_addToBat_corr);
//    PMEDouble myWeight = cSim.TIItemWeight[Schedule::TypeBAT][TIRegion];

    xij = cSim.pImageX[iatom[0]] - cSim.pImageX[iatom[1]];
    yij = cSim.pImageY[iatom[0]] - cSim.pImageY[iatom[1]];
    zij = cSim.pImageZ[iatom[0]] - cSim.pImageZ[iatom[1]];
    rij = sqrt(xij*xij + yij * yij + zij * zij);
    da = rij - bond.y;
    df = bond.x * da;
    dfw = (df + df) / rij;

//    if (addToDVDL) {
//      pE = cSim.pTIEBond[TIRegion];
//      pFx = pTIFx[TIRegion];
//      pFy = pTIFy[TIRegion];
//      pFz = pTIFz[TIRegion];
//    } else {
//      pE = cSim.pTISCBond[TIRegion];
//      pFx = pSCFx[TIRegion];
//      pFy = pSCFy[TIRegion];
//      pFz = pSCFz[TIRegion];
//    }
    pE = cSim.pTISCBond[TIRegion];
    if (needPotEnergy) {
      tt = df * da;
      addEnergy(pE, tt);
//      if (hasSC && !inSC) {
//        addEnergy(pSE1[TIRegion], tt);
//        if (!addToBATCorr) {
//          addEnergy(pSE2[TIRegion], tt);
//        }
//      }
    }

    // Add to the force accumulators
//    tw = addToDVDL ? myWeight : OneF;
    tw=OneF;
    if (fabs(tw) > tm06 && (cSim.igamd == 111||cSim.igamd == 118||cSim.igamd == 14 ||cSim.igamd == 15||cSim.igamd == 21||cSim.igamd == 22||cSim.igamd == 24||cSim.igamd == 25||cSim.igamd == 27)) {
      tt = dfw * xij * tw;
//      addForce(pFx + iatom[0], -tt);
//      addForce(pFx + iatom[1], tt);
    addForce(cSim.pGaMDTIForceX[TIRegion] + iatom[0], -tt);
    addForce(cSim.pGaMDTIForceX[TIRegion] + iatom[1], tt);

      tt = dfw * yij * tw;
//      addForce(pFy + iatom[0], -tt);
//      addForce(pFy + iatom[1], tt);
    addForce(cSim.pGaMDTIForceY[TIRegion] + iatom[0], -tt);
    addForce(cSim.pGaMDTIForceY[TIRegion] + iatom[1], tt);

      tt = dfw * zij * tw;
//      addForce(pFz + iatom[0], -tt);
//      addForce(pFz + iatom[1], tt);
    addForce(cSim.pGaMDTIForceZ[TIRegion] + iatom[0], -tt);
    addForce(cSim.pGaMDTIForceZ[TIRegion] + iatom[1], tt);
    }

  }
  if (pos >= angleBegin && pos <= angleEnd) {
    unsigned mypos = pos - angleBegin;
    uint4& angleID = cSim.pTIBondAngleID[mypos];
    PMEDouble2& angle = cSim.pTIBondAngle[mypos];
    iatom[0] = cSim.pImageAtomLookup[angleID.x];
    iatom[1] = cSim.pImageAtomLookup[angleID.y];
    iatom[2] = cSim.pImageAtomLookup[angleID.z];

    TIRegion = (cSim.pTIBondAngleType[mypos] & gti_simulationConst::Fg_TI_region);
    //bool hasSC = (cSim.pTIBondAngleType[mypos] & gti_simulationConst::Fg_has_SC);
    //bool inSC = (cSim.pTIBondAngleType[mypos] & gti_simulationConst::Fg_int_SC);
//    bool addToDVDL = (cSim.pTIBondAngleType[mypos] & gti_simulationConst::ex_addToDVDL_bat);
//    bool addToBATCorr = (cSim.pTIBondAngleType[mypos] & gti_simulationConst::ex_addToBat_corr);
//    PMEDouble myWeight = cSim.TIItemWeight[Schedule::TypeBAT][TIRegion];

    xij = cSim.pImageX[iatom[0]] - cSim.pImageX[iatom[1]];
    yij = cSim.pImageY[iatom[0]] - cSim.pImageY[iatom[1]];
    zij = cSim.pImageZ[iatom[0]] - cSim.pImageZ[iatom[1]];
    xkj = cSim.pImageX[iatom[2]] - cSim.pImageX[iatom[1]];
    ykj = cSim.pImageY[iatom[2]] - cSim.pImageY[iatom[1]];
    zkj = cSim.pImageZ[iatom[2]] - cSim.pImageZ[iatom[1]];
    rij = xij * xij + yij * yij + zij * zij;

    PMEDouble rkj = xkj * xkj + ykj * ykj + zkj * zkj;
    PMEDouble rik = sqrt(rij * rkj);

    const PMEDouble nlimit = -0.99999;
    const PMEDouble plimit = 0.99999;
    PMEDouble cst = min(plimit, max(nlimit, (xij*xkj + yij * ykj + zij * zkj) / rik));
    PMEDouble ant = acos(cst);

    da = ant - angle.y;
    df = angle.x * da;
    dfw = -(df + df) / sin(ant);

    //pE0 = cSim.pTIEAngle[2];
//    if (addToDVDL) {
//      pE = cSim.pTIEAngle[TIRegion];
//      pFx = pTIFx[TIRegion];
//      pFy = pTIFy[TIRegion];
//      pFz = pTIFz[TIRegion];
//    } else {
//      pE = cSim.pTISCAngle[TIRegion];
//      pFx = pSCFx[TIRegion];
//      pFy = pSCFy[TIRegion];
//      pFz = pSCFz[TIRegion];
//    }
    pE = cSim.pTISCAngle[TIRegion];
    if (needPotEnergy) {
      tt = df * da;
      addEnergy(pE, tt);

//      if (hasSC && !inSC) {
//        addEnergy(pSE1[TIRegion], tt);
//        if (!addToBATCorr) {
//          addEnergy(pSE2[TIRegion], tt);
//        }
//      }
    }

    // Calculation of the force
    cik = dfw / rik;
    sth = dfw * cst;
    cii = sth / rij;
    ckk = sth / rkj;
    fx1 = cii * xij - cik * xkj;
    fy1 = cii * yij - cik * ykj;
    fz1 = cii * zij - cik * zkj;
    fx2 = ckk * xkj - cik * xij;
    fy2 = ckk * ykj - cik * yij;
    fz2 = ckk * zkj - cik * zij;
//    tw = addToDVDL ? myWeight : 1.0;
    tw=1.0;
    if (fabs(tw) > tm06 && (cSim.igamd == 14 ||cSim.igamd == 15 ||cSim.igamd == 21 ||cSim.igamd == 22 ||cSim.igamd == 111||cSim.igamd == 113||cSim.igamd == 119||cSim.igamd == 120||cSim.igamd == 24||cSim.igamd == 25||cSim.igamd == 27)) {
    addForce(cSim.pGaMDTIForceX[TIRegion] + iatom[0], fx1* tw);
    addForce(cSim.pGaMDTIForceY[TIRegion] + iatom[0], fy1* tw);
    addForce(cSim.pGaMDTIForceZ[TIRegion] + iatom[0], fz1* tw);

    addForce(cSim.pGaMDTIForceX[TIRegion] + iatom[2], fx2* tw);
    addForce(cSim.pGaMDTIForceY[TIRegion] + iatom[2], fy2 * tw);
    addForce(cSim.pGaMDTIForceZ[TIRegion] + iatom[2], fz2 * tw);

    addForce(cSim.pGaMDTIForceX[TIRegion] + iatom[1], -fx1 - fx2);
    addForce(cSim.pGaMDTIForceY[TIRegion] + iatom[1], -fy1 - fy2);
    addForce(cSim.pGaMDTIForceZ[TIRegion] + iatom[1], -fz1 - fz2);

//      addForce(pFx + iatom[0], fx1 * tw);
//      addForce(pFy + iatom[0], fy1 * tw);
//      addForce(pFz + iatom[0], fz1 * tw);

//      addForce(pFx + iatom[2], fx2 * tw);
//      addForce(pFy + iatom[2], fy2 * tw);
//      addForce(pFz + iatom[2], fz2 * tw);

//      addForce(pFx + iatom[1], -tw * (fx1 + fx2));
//      addForce(pFy + iatom[1], -tw * (fy1 + fy2));
//      addForce(pFz + iatom[1], -tw * (fz1 + fz2));
    }
  }
  if (pos >= dihedralBegin && pos <= dihedralEnd) {
    unsigned mypos = pos - dihedralBegin;
    uint4& dihedralID = cSim.pTIDihedralID[mypos];
    iatom[0] = cSim.pImageAtomLookup[dihedralID.x];
    iatom[1] = cSim.pImageAtomLookup[dihedralID.y];
    iatom[2] = cSim.pImageAtomLookup[dihedralID.z];
    iatom[3] = cSim.pImageAtomLookup[dihedralID.w];

    PMEDouble2& dihedral1 = cSim.pTIDihedral1[mypos];
    PMEDouble2& dihedral2 = cSim.pTIDihedral2[mypos];
    PMEDouble& dihedral3 = cSim.pTIDihedral3[mypos];

    TIRegion = (cSim.pTIDihedralType[mypos] & gti_simulationConst::Fg_TI_region);
    //bool hasSC = (cSim.pTIDihedralType[mypos] & gti_simulationConst::Fg_has_SC);
    //bool inSC = (cSim.pTIDihedralType[mypos] & gti_simulationConst::Fg_int_SC);
//    bool addToDVDL = (cSim.pTIDihedralType[mypos] & gti_simulationConst::ex_addToDVDL_bat);
//    bool addToBATCorr = (cSim.pTIDihedralType[mypos] & gti_simulationConst::ex_addToBat_corr);
    //PMEDouble myWeight = cSim.TIItemWeight[Schedule::TypeBAT][TIRegion];
    xij = cSim.pImageX[iatom[0]] - cSim.pImageX[iatom[1]];
    yij = cSim.pImageY[iatom[0]] - cSim.pImageY[iatom[1]];
    zij = cSim.pImageZ[iatom[0]] - cSim.pImageZ[iatom[1]];
    xkj = cSim.pImageX[iatom[2]] - cSim.pImageX[iatom[1]];
    ykj = cSim.pImageY[iatom[2]] - cSim.pImageY[iatom[1]];
    zkj = cSim.pImageZ[iatom[2]] - cSim.pImageZ[iatom[1]];
    xkl = cSim.pImageX[iatom[2]] - cSim.pImageX[iatom[3]];
    ykl = cSim.pImageY[iatom[2]] - cSim.pImageY[iatom[3]];
    zkl = cSim.pImageZ[iatom[2]] - cSim.pImageZ[iatom[3]];

    // Get the normal vector
    dx = yij * zkj - zij * ykj;
    dy = zij * xkj - xij * zkj;
    dz = xij * ykj - yij * xkj;
    gx = zkj * ykl - ykj * zkl;
    gy = xkj * zkl - zkj * xkl;
    gz = ykj * xkl - xkj * ykl;
    fxi = sqrt(dx*dx + dy * dy + dz * dz + tm24);
    fyi = sqrt(gx*gx + gy * gy + gz * gz + tm24);
    cst = dx * gx + dy * gy + dz * gz;

    // Branch if linear dihedral:
    PMEDouble z1 = (tenm3 <= fxi) ? (1.0 / fxi) : zero;
    PMEDouble z2 = (tenm3 <= fyi) ? (1.0 / fyi) : zero;
    PMEDouble z12 = z1 * z2;
    PMEDouble fzi = (z12 != zero) ? one : zero;
    PMEDouble s = xkj * (dz*gy - dy * gz) + ykj * (dx*gz - dz * gx) + zkj * (dy*gx - dx * gy);
    PMEDouble ff = PI - abs(acos(max(-1.0, min(1.0, cst * z12)))) * (s >= zero ? one : -one);
    faster_sincos(ff, &sphi, &cphi);

    // Calculate the energy and the derivatives with respect to cosphi
    PMEDouble ct0 = dihedral1.y * ff;
    PMEDouble sinnp, cosnp;
    faster_sincos(ct0, &sinnp, &cosnp);

    //pE0 = cSim.pTIEDihedral[2];
//    if (addToDVDL) {
//      pE = cSim.pTIEDihedral[TIRegion];
//      pFx = pTIFx[TIRegion];
//      pFy = pTIFy[TIRegion];
//      pFz = pTIFz[TIRegion];
//    } else {
//      pE = cSim.pTISCDihedral[TIRegion];
//      pFx = pSCFx[TIRegion];
//      pFy = pSCFy[TIRegion];
//      pFz = pSCFz[TIRegion];
//    }
    pE = cSim.pTISCDihedral[TIRegion];
    tt = (dihedral2.x + (cosnp * dihedral2.y) + (sinnp * dihedral3)) * fzi;
    if (needPotEnergy) {
      addEnergy(pE, tt);
//      if (addToDVDL) {
//        addEnergy(pE0, -tt);
//      }

      //if (hasSC && !inSC) {
      //  addEnergy(pSE1[TIRegion], tt);
      //  if (!addToBATCorr) {
      //    addEnergy(pSE2[TIRegion], tt);
      //  }
      //}
    }


#ifdef use_DPFP
    PMEFloat delta = tm24;
#else
    PMEFloat delta = tm06;
#endif
    PMEDouble dums = sphi + (delta * (sphi >= zero ? one : -one));
    PMEDouble df;
    if (tm06 > abs(dums)) {
      df = fzi * dihedral2.y * (dihedral1.y - dihedral1.x + (dihedral1.x * cphi));
    } else {
      df = fzi * dihedral1.y * ((dihedral2.y * sinnp) - (dihedral3 * cosnp)) / dums;
    }

    // Now, set up array dc = 1st der. of cosphi w/respect to cartesian differences:
    z11 = z1 * z1;
    z12 = z1 * z2;
    z22 = z2 * z2;
    dcdx = -gx * z12 - cphi * dx*z11;
    dcdy = -gy * z12 - cphi * dy*z11;
    dcdz = -gz * z12 - cphi * dz*z11;
    dcgx = dx * z12 + cphi * gx*z22;
    dcgy = dy * z12 + cphi * gy*z22;
    dcgz = dz * z12 + cphi * gz*z22;

    // Update the first derivative array:
    dr1 = df * (dcdz*ykj - dcdy * zkj);
    dr2 = df * (dcdx*zkj - dcdz * xkj);
    dr3 = df * (dcdy*xkj - dcdx * ykj);
    dr4 = df * (dcgz*ykj - dcgy * zkj);
    dr5 = df * (dcgx*zkj - dcgz * xkj);
    dr6 = df * (dcgy*xkj - dcgx * ykj);
    drx = df * (-dcdy * zij + dcdz * yij + dcgy * zkl - dcgz * ykl);
    dry = df * (dcdx*zij - dcdz * xij - dcgx * zkl + dcgz * xkl);
    drz = df * (-dcdx * yij + dcdy * xij + dcgx * ykl - dcgy * xkl);
if(cSim.igamd == 14 ||cSim.igamd == 15 || cSim.igamd == 18 ||cSim.igamd == 21 ||cSim.igamd == 22||cSim.igamd == 111||cSim.igamd == 113||cSim.igamd == 115||cSim.igamd == 117||cSim.igamd == 120||cSim.igamd == 24||cSim.igamd == 25||cSim.igamd == 27){
//    tw = addToDVDL ? myWeight : 1.0;
    tw = 1.0;
    addForce(cSim.pGaMDTIForceX[TIRegion] + iatom[0], (-dr1));
    addForce(cSim.pGaMDTIForceY[TIRegion] + iatom[0], (-dr2));
    addForce(cSim.pGaMDTIForceZ[TIRegion] + iatom[0], (-dr3));

    addForce(cSim.pGaMDTIForceX[TIRegion]+ iatom[1], (-drx + dr1));
    addForce(cSim.pGaMDTIForceY[TIRegion]+ iatom[1], (-dry + dr2));
    addForce(cSim.pGaMDTIForceZ[TIRegion] + iatom[1], (-drz + dr3));

    addForce(cSim.pGaMDTIForceX[TIRegion] + iatom[2], (drx + dr4));
    addForce(cSim.pGaMDTIForceY[TIRegion] + iatom[2], (dry + dr5));
    addForce(cSim.pGaMDTIForceZ[TIRegion] + iatom[2], (drz + dr6));

    addForce(cSim.pGaMDTIForceX[TIRegion] + iatom[3], (-dr4));
    addForce(cSim.pGaMDTIForceY[TIRegion] + iatom[3], (-dr5));
    addForce(cSim.pGaMDTIForceZ[TIRegion] + iatom[3], (-dr6));

//    addForce(pFx + iatom[0], (-dr1*tw));
//    addForce(pFy + iatom[0], (-dr2 * tw));
//    addForce(pFz + iatom[0], (-dr3 * tw));
//    addForce(pTIFx[2] + iatom[0], dr1);
//    addForce(pTIFy[2] + iatom[0], dr2);
//    addForce(pTIFz[2] + iatom[0], dr3);

//    addForce(pFx + iatom[1], (-drx + dr1)* tw);
//    addForce(pFy + iatom[1], (-dry + dr2)* tw);
//    addForce(pFz + iatom[1], (-drz + dr3)* tw);
//    addForce(pTIFx[2] + iatom[1], (drx - dr1));
//    addForce(pTIFy[2] + iatom[1], (dry - dr2));
//    addForce(pTIFz[2] + iatom[1], (drz - dr3));
//    addForce(pFx + iatom[2], (drx + dr4)* tw);
//    addForce(pFy + iatom[2], (dry + dr5)* tw);
//    addForce(pFz + iatom[2], (drz + dr6)* tw);
//    addForce(pTIFx[2] + iatom[2], -(drx + dr4));
//    addForce(pTIFy[2] + iatom[2], -(dry + dr5));
//    addForce(pTIFz[2] + iatom[2], -(drz + dr6));

//    addForce(pFx + iatom[3], (-dr4 * tw));
//    addForce(pFy + iatom[3], (-dr5 * tw));
//    addForce(pFz + iatom[3], (-dr6 * tw));
//    addForce(pTIFx[2] + iatom[3], dr4);
//    addForce(pTIFy[2] + iatom[3], dr5);
//    addForce(pTIFz[2] + iatom[3], dr6);
   }
  }
}

#endif
