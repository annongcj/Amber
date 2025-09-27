#ifndef __GTI_GPU_H__
#define __GTI_GPU_H__

#  ifdef GTI
#    include "gpuContext.h"

// C++ routines / interfaces

void icc_updateAllSimulationConst(gpuContext gpu);

void icc_CalculateElecRecipForceEnergy(gpuContext gpu, double vol, bool needPotEnergy,
                                       int ti_mode, int TIRegion, int reaf_mode);

void icc_GetEnergyFromGPU(gpuContext gpu, double energy[EXTENDED_ENERGY_TERMS]);

void icc_GetTIPotEnergyFromGPU(gpuContext gpu,
                               double energy[3][gti_simulationConst::GPUPotEnergyTerms]);

void icc_GetTIKinEnergyFromGPU(gpuContext gpu,
                               double energy[3][gti_simulationConst::GPUKinEnergyTerms]);

void icc_GetMBAREnergyFromGPU(gpuContext gpu, double MBAREnergy[]);

#  endif /* GTI */
#endif /* __GTI_GPU_H__ */
