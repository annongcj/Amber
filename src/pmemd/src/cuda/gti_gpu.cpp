#ifdef GTI

#include <assert.h>
#include "gti_def.h"
#include "gpuContext.h"
#include "gpu.h"
#include "gti_gpu.h"
#include "gti_cuda.cuh"

//---------------------------------------------------------------------------------------------
// gti_gpuImpl: 
//---------------------------------------------------------------------------------------------
namespace gti_gpuImpl {
  inline double EnergyConverter(unsigned long long int val, bool useEnergyScale=true);
  inline double IntegerConverter(unsigned long long int val);  
}

using namespace gti_gpuImpl;

void icc_updateAllSimulationConst(gpuContext gpu) {
  ik_updateAllSimulationConst(gpu);
}

//---------------------------------------------------------------------------------------------
// icc_CalculateElecRecipForceEnergy:
//
// Arguments:
//
//---------------------------------------------------------------------------------------------
void icc_CalculateElecRecipForceEnergy(gpuContext gpu, double vol, bool needPotEnergy,
                                       int ti_mode, int TIRegion, int reaf_mode)
{
  if (gpu->ntf != 8) {

    // This is necessary for correct FFT force calc.
    if (ti_mode > 0) {
      (reaf_mode >= 0) ? ik_ZeroTICharge(gpu, reaf_mode)
        : ik_ZeroTICharge(gpu, TIRegion);
    } 
    
    if (reaf_mode >= 0) {
      ik_ScaleRECharge(gpu, reaf_mode);
    }
        
    // Zero the charge grid, if that has not already happened
    ik_PrepareChargeGrid(gpu, ti_mode, TIRegion);

    kPMEForwardFFT(gpu);
    if (needPotEnergy) {
      kPMEScalarSumRCEnergy(gpu, vol);
    }
    else {
      kPMEScalarSumRC(gpu, vol);
    }
    kPMEBackwardFFT(gpu);
    kPMEGradSum(gpu);
  }
}

//---------------------------------------------------------------------------------------------
// icc_GetEnergyFromGPU:
//
// Arguments:
//
//---------------------------------------------------------------------------------------------
void icc_GetEnergyFromGPU(gpuContext gpu, double energy[EXTENDED_ENERGY_TERMS])
{
  gpu->pbEnergyBuffer->Download();
  for (int i = 0; i < ENERGY_TERMS; i++) {
    bool useEnergyScale = (i<VIRIAL_OFFSET);
    energy[i] = EnergyConverter(gpu->pbEnergyBuffer->_pSysData[i], useEnergyScale);
  }

#ifdef use_DPFP
  if (gpu->imin != 0)
  {
    for (int i = ENERGY_TERMS; i < EXTENDED_ENERGY_TERMS; i++)
    {
      unsigned long long int val = gpu->pbEnergyBuffer->_pSysData[i];
      energy[i] = IntegerConverter(gpu->pbEnergyBuffer->_pSysData[i]);
    }
    energy[1] += energy[EVDWE_OFFSET];
    energy[10] += energy[EEDE_OFFSET];
    energy[VIRIAL_OFFSET + 0] += energy[VIRIALE_OFFSET + 0];
    energy[VIRIAL_OFFSET + 1] += energy[VIRIALE_OFFSET + 1];
    energy[VIRIAL_OFFSET + 2] += energy[VIRIALE_OFFSET + 2];
  }
#endif  
}

//---------------------------------------------------------------------------------------------
// icc_GetTIPotEnergyFromGPU:
//
// Arguments:
//
//---------------------------------------------------------------------------------------------
void icc_GetTIPotEnergyFromGPU(gpuContext gpu,
                               double energy[gti_simulationConst::TIEnergyBufferMultiplier][gti_simulationConst::GPUPotEnergyTerms])
{
  gpu->pbTIPotEnergyBuffer->Download();

  unsigned long long int *gputi;
    gputi = gpu->pbTIPotEnergyBuffer->_pSysData;
  for (int ti = 0; ti < gti_simulationConst::TIEnergyBufferMultiplier; ti++){
    for (int i = 0; i < gti_simulationConst::GPUPotEnergyTerms; i++) {
      bool useEnergyScale = (i<VIRIAL_OFFSET || i>(VIRIAL_OFFSET + 5));
      energy[ti][i] = EnergyConverter(gputi[i + ti*gti_simulationConst::GPUPotEnergyTerms],
                                      useEnergyScale);
    }
  }
}

//---------------------------------------------------------------------------------------------
// icc_GetTIKinEnergyFromGPU:
//
// Arguments:
//   gpu:            overarching structure holding all simulation information, critically in
//                   this case energy buffers
//---------------------------------------------------------------------------------------------
void icc_GetTIKinEnergyFromGPU(gpuContext gpu,
                               double energy[3][gti_simulationConst::GPUKinEnergyTerms])
{
  gpu->pbTIKinEnergyBuffer->Download();

  unsigned long long int *gputi;
  gputi = gpu->pbTIKinEnergyBuffer->_pSysData;
  for (int ti = 0; ti < 3; ti++) {
    for (int i = 0; i < gti_simulationConst::GPUKinEnergyTerms; i++) {
      energy[ti][i] = EnergyConverter(gputi[i + ti*gti_simulationConst::GPUKinEnergyTerms]);
    }
  }
}

//---------------------------------------------------------------------------------------------
// icc_GetMBAREnergyFromGPU:
//
// Arguments:
//   gpu:            overarching structure holding all simulation information, critically in
//                   this case energy buffers
//   energy:  
//---------------------------------------------------------------------------------------------
void icc_GetMBAREnergyFromGPU(gpuContext gpu, double energy[])
{

  if (gpu->sim.nMBARStates > 2001) return;

  const double MaxEnergyDiff=200.0;
  int n= gpu->sim.nMBARStates* Schedule::TypeTotal * 3;
  if (n <= 0) {
    return;
  }
  gpu->pbMBAREnergy->Download();
  for (int i = 0; i < n; i++) {

    double tt = EnergyConverter(gpu->pbMBAREnergy->_pSysData[i]);
    energy[i] = (std::isfinite(tt)) ? tt : MaxEnergyDiff;
  }
}

//---------------------------------------------------------------------------------------------
// EnergyConverter:
//
// Arguments:
//   val:             
//   useEnergyScale:  
//---------------------------------------------------------------------------------------------
double gti_gpuImpl::EnergyConverter(unsigned long long int val, bool useEnergyScale)
{
  double energy = 0.0;
  double factor = (useEnergyScale) ? 1.0/eScale : 1.0/vScale;
  if (val >= 0x8000000000000000ull) {
    energy = -(PMEDouble)(val ^ 0xffffffffffffffffull) * factor;
  }
  else {
    energy = (PMEDouble)val * factor;
  }

  return energy;
}

double gti_gpuImpl::IntegerConverter(unsigned long long int val)
{
  double energy = 0.0;
  energy = (double)(*((long long int*)(&val)));
  return energy;
}

#endif
