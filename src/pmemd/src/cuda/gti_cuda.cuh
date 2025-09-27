// Interfaces to Kernels
#ifndef __GTI_CUDA_H__
#define __GTI_CUDA_H__

#  ifdef GTI

#    include "gpuContext.h"

// gpu-Kernel access

void ik_updateAllSimulationConst(gpuContext gpu);

void ik_ClearTIEnergyForce_kernel(gpuContext gpu);

void ik_ClearTIEnergyForce_gamd_kernel(gpuContext gpu);

void ik_ZeroTICharge(gpuContext gpu, unsigned int mode);

void ik_ScaleRECharge(gpuContext gpu, unsigned int mode);

void ik_ZeroTIAtomForce(gpuContext gpu, unsigned int mode);

void ik_PrepareChargeGrid(gpuContext gpu, int ti_mode, int TIRegion);

void ik_CopyToTIForce(gpuContext gpu, int TIRegion, bool isNB, bool keepSource, PMEFloat weight);

void ik_CopyToTIEnergy(gpuContext gpu, int TIRegion, int term1, int term2, int term3,
                       bool isVirial, PMEFloat weight, bool addon=false);

void ik_CorrectTIEnergy(gpuContext gpu, int beginTerm, int endTerm);

void ik_CorrectTIForce(gpuContext gpu);

void ik_CombineTIForce(gpuContext gpu, bool linear, bool needvirial);

void ik_CombineTIForce_gamd(gpuContext gpu, bool linear, bool needvirial);

void ik_RemoveTINetForce(gpuContext gpu, int TIRegion);

void ik_SyncVector_kernel(gpuContext gpu, unsigned int mode, int combinedMode);

void ik_CalculateTIKineticEnergy(gpuContext gpu, double c_ave);

void ik_CalculateTIBonded(gpuContext gpu, bool needEnergy, bool needVirial);

void ik_CalculateREAFBonded(gpuContext gpu, bool needEnergy, bool needVirial);

void ik_CalculateTIRestraint(gpuContext gpu, bool needEnergy, bool needVirial);

void ik_CalculateRMSD(gpuContext gpu, bool needPotEnergy);

void ik_CalculateTIBonded_ppi(gpuContext gpu, bool needEnergy, bool needVirial);

void ik_BuildTINBList(gpuContext gpu);

void ik_Build1264NBList(gpuContext gpu);

void ik_Buildp1264NBList(gpuContext gpu); // C4PairwiseCUDA

void ik_Build1264p1264NBList(gpuContext gpu); // C4PairwiseCUDA

void ik_BuildREAFNbList(gpuContext gpu);

void ik_BuildTINBList_gamd(gpuContext gpu);

void ik_CalculateTINB(gpuContext gpu, bool needEnergy, bool needVirial);

void ik_CalculateTI14NB(gpuContext gpu, bool needEnergy, bool needVirial);

void ik_CalculateREAFNb(gpuContext gpu, bool needEnergy, bool needVirial);

void ik_CalculateREAF14NB(gpuContext gpu, bool needEnergy, bool needVirial);

void ik_CalculateTINB_gamd(gpuContext gpu, bool needEnergy, bool needVirial);

void ik_CalculateTINB_ppi_gamd(gpuContext gpu, bool needEnergy, bool needVirial);

void ik_CalculateTINB_ppi2_gamd(gpuContext gpu, bool needEnergy, bool needVirial, int bgpro2atm, int edpro2atm);

void ik_CalculateTI14NB_gamd(gpuContext gpu, bool needEnergy, bool needVirial);

void ik_CalculateTI14NB_ppi_gamd(gpuContext gpu, bool needEnergy, bool needVirial);



void ik_Build1264NBList_gamd(gpuContext gpu);

void ik_Buildp1264NBList_gamd(gpuContext gpu); // C4PairwiseCUDA

void ik_Build1264p1264NBList_gamd(gpuContext gpu); // C4PairwiseCUDA

void ik_Calculate1264NB(gpuContext gpu, bool needEnergy, bool needVirial);

void ik_Calculatep1264NB(gpuContext gpu, bool needEnergy, bool needVirial); // C4PairwiseCUDA

void ik_Calculate1264p1264NB(gpuContext gpu, bool needEnergy, bool needVirial); // C4PairwiseCUDA

#  endif /* GTI */

#endif /* __GTI_CUDA_H__ */
