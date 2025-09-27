
#ifdef GTI

#include <assert.h>
#ifdef AMBER_PLATFORM_AMD
#include <hip/hip_runtime.h>
#else
#include <cuda_runtime.h>
#include <cusolverDn.h>
#endif
//#include <cusolver_utils.h>
#include "gputypes.h"
#include "gpu.h"
#include "gti_general_kernels.cuh"
#include "gti_energies_kernels.cuh"
#include "gti_NBList_kernels.cuh"
#include "gti_nonBond_kernels.cuh"
#include "gti_cuda.cuh"
#include <iostream> // C4PairwiseCUDA2023

#if defined(__CUDACC_RDC__) || defined(__HIPCC_RDC__)
#  include "simulationConst.h"
  __device__ __constant__ simulationConst cSim;
#  include "gti_utils.inc"
#endif

void ik_updateAllSimulationConst(gpuContext gpu) {

#if defined(__CUDACC_RDC__) || defined(__HIPCC_RDC__)
  UpdateSimulationConst();
#else
  // Original pmemd cu units
  SetkForcesUpdateSim(gpu);
  SetkDataTransferSim(gpu);
  SetkCalculateNEBForcesSim(gpu);
  SetkShakeSim(gpu);
  SetkCalculateLocalForcesSim(gpu);
  SetkXrayUpdateSim(gpu);
  if (gpu->ntb == 0) {
    SetkCalculateGBBornRadiiSim(gpu);
    SetkCalculateGBNonbondEnergy1Sim(gpu);
    SetkCalculateGBNonbondEnergy2Sim(gpu);
  } else {
    if (gpu->sim.icnstph == 2 || gpu->sim.icnste == 2) {
      SetkCalculateGBBornRadiiSim(gpu);
      SetkCalculateGBNonbondEnergy1Sim(gpu);
      SetkCalculateGBNonbondEnergy2Sim(gpu);
    }
    SetkPMEInterpolationSim(gpu);
    SetkNeighborListSim(gpu);
    SetkCalculateEFieldEnergySim(gpu);
    SetkCalculatePMENonbondEnergySim(gpu);
  }

  // GTI cu units
  GTI_GENERAL_IMPL::UpdateSimulationConst();
  GTI_ENERGY_IMPL::UpdateSimulationConst();
  GTI_NB_LIST::UpdateSimulationConst();
  GTI_NB_CALC::UpdateSimulationConst();

#endif
};
//---------------------------------------------------------------------------------------------
// ik_CopyToTIForce:
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                system size parameters
//   TIRegion:    numerical descriptor of the thermodynamic integration zone to handle
//   isNB:
//   keepSource:
//---------------------------------------------------------------------------------------------
void ik_CopyToTIForce(gpuContext gpu, int TIRegion, bool isNB, bool keepSource, PMEFloat weight)
{
  // Must use the default stream
  int nterms = gpu->sim.stride3;
  if (nterms <= 0) {
    return;
  }
  unsigned threadsPerBlock = 512;
  if (gpu->major == 6) {
    threadsPerBlock = 128;
  }
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);
  unsigned blocksToUse = (nterms / threadsPerBlock / 8) + 1;

  kgCopyToTIForce<<<blocksToUse, threadsPerBlock>>>(TIRegion, isNB, keepSource, weight);
  LAUNCHERROR("kCopyToTIForce");
}

//---------------------------------------------------------------------------------------------
// ik_CombineTIForce:
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                system size parameters
//   linear:
//   needvirial:
//---------------------------------------------------------------------------------------------
void ik_CombineTIForce(gpuContext gpu, bool linear, bool needvirial) {
  // Must use the default stream
  int nterms = gpu->sim.stride3;
  if (nterms <= 0) {
    return;
  }

  if (gpu->multiStream) cudaEventSynchronize(gpu->event_TIDone);

  unsigned threadsPerBlock = 512;
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);
  unsigned blocksToUse = (nterms / threadsPerBlock / 4) + 1;

  bool isNB = false;
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));
  bool addOn = true;
  bool needWeight = false;
  bool hasSC = ((gpu->sim.numberSoftcoreAtoms[0] + gpu->sim.numberSoftcoreAtoms[1]) > 0);

  kgCombineTIForce_kernel <<<blocksToUse, threadsPerBlock, 0,
    gpu->mainStream >>>(isNB, useImage, addOn, needWeight);
  LAUNCHERROR("kgCombineTIForce_kernel");

  if (hasSC) {
    kgCombineSCForce_kernel <<<blocksToUse, threadsPerBlock, 0, gpu->mainStream >>>(isNB, useImage);
    LAUNCHERROR("kgCombineSCForce_kernel");
  }

  if (needvirial) {
    bool isNB = true;
    kgCombineTIForce_kernel <<<blocksToUse, threadsPerBlock,
      0, gpu->mainStream >>>(isNB, useImage, addOn, needWeight);
    LAUNCHERROR("kgCombineTIForce_kernel");

    if (hasSC) {
    kgCombineSCForce_kernel <<<blocksToUse, threadsPerBlock, 0, gpu->mainStream >>>(isNB, useImage);
    LAUNCHERROR("kgCombineSCForce_kernel");
    }
  }

}

//---------------------------------------------------------------------------------------------
// ik_CombineTIForce_gamd:
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                system size parameters
//   linear:
//   needvirial:
//---------------------------------------------------------------------------------------------
void ik_CombineTIForce_gamd(gpuContext gpu, bool linear, bool needvirial) {
  // Must use the default stream
  int nterms = gpu->sim.stride3;
  if (nterms <= 0) {
    return;
  }

  if (gpu->multiStream) cudaEventSynchronize(gpu->event_TIDone);

  unsigned threadsPerBlock = 512;
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);
  unsigned blocksToUse = (nterms / threadsPerBlock / 4) + 1;

  bool isNB = false;
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));
  bool addOn = true;
  bool hasSC = ((gpu->sim.numberSoftcoreAtoms[0] + gpu->sim.numberSoftcoreAtoms[1]) > 0);

  kgCombineTIForce_gamd_kernel <<<blocksToUse, threadsPerBlock, 0,
    gpu->mainStream >>>(isNB, useImage, addOn);
  LAUNCHERROR("kgCombineTIForce_kernel");

  if (hasSC) {
    kgCombineSCForce_gamd_kernel <<<blocksToUse, threadsPerBlock, 0, gpu->mainStream >>>(isNB, useImage);
    LAUNCHERROR("kgCombineSCForce_kernel");
  }

  if (needvirial) {
    bool isNB = true;
    kgCombineTIForce_gamd_kernel <<<blocksToUse, threadsPerBlock,
      0, gpu->mainStream >>>(isNB, useImage, addOn);
    LAUNCHERROR("kgCombineTIForce_kernel");

    if (hasSC) {
    kgCombineSCForce_gamd_kernel <<<blocksToUse, threadsPerBlock, 0, gpu->mainStream >>>(isNB, useImage);
    LAUNCHERROR("kgCombineSCForce_kernel");
    }
  }

}

//---------------------------------------------------------------------------------------------
// ik_RemoveTINetForce:
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                kernel launch parameters
//   TIRegion:
//---------------------------------------------------------------------------------------------
void ik_RemoveTINetForce(gpuContext gpu, int TIRegion)
{
  cudaStream_t stream = (TIRegion<0) ? gpu->mainStream : gpu->TIStream;

  int threadsPerBlock = gpu->generalThreadsPerBlock;
#if ( defined(_WIN32) && defined(_DEBUG) )
  threadsPerBlock = 128;
#endif
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);

  kgCalculateNetForce_kernel<<<gpu->blocks, threadsPerBlock,
                               0, stream>>>(TIRegion);
  LAUNCHERROR("kgCalculateNetForce_kernel");


  kgRemoveNetForce_kernel<<<gpu->blocks, threadsPerBlock, 0, stream>>>(TIRegion);
  LAUNCHERROR("kgRemoveNetForce_kernel");
}

//---------------------------------------------------------------------------------------------
// ik_CopyToTIEnergy: overloaded form of the eponymous function above
//
// Arguments:
//   As above, but replacing being and end with...
//   term1:
//   term2:
//   term3:
//---------------------------------------------------------------------------------------------
void ik_CopyToTIEnergy(gpuContext gpu, int TIRegion, int term1, int term2, int term3,
                       bool isVirial, PMEFloat tiWeight, bool addon)
{
  // mode=0 : replace
  // mode=1 : add onto
  unsigned blocksToUse = 1;
  unsigned threadsPerBlock = GRID;
  kgCopyToTIEnergy_kernel<<<blocksToUse, threadsPerBlock>>>(TIRegion, term1, term2, term3,
                                                            isVirial, tiWeight, addon);
  LAUNCHERROR("kCopyToTIEnergy_kernel");
}

//---------------------------------------------------------------------------------------------
// ik_CorrectTIEnergy:
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                TI stream parameters
//   beginTerm:
//   endTerm:
//---------------------------------------------------------------------------------------------
void ik_CorrectTIEnergy(gpuContext gpu, int beginTerm, int endTerm)
{
  unsigned blocksToUse = 1;
  unsigned threadsPerBlock = GRID;
  kgCorrectTIEnergy_kernel<<<blocksToUse, threadsPerBlock, 0,
                             gpu->TIStream>>>(beginTerm, endTerm);
  LAUNCHERROR("kCorrectTIEnergy_kernel");
}

//---------------------------------------------------------------------------------------------
// ik_CorrectTIForce:
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                TI stream parameters
//---------------------------------------------------------------------------------------------
void ik_CorrectTIForce(gpuContext gpu)
{
  unsigned blocksToUse = gpu->blocks;
  unsigned threadsPerBlock = gpu->generalThreadsPerBlock;
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);
  kgCorrectTIForce_kernel<<<blocksToUse, threadsPerBlock, 0, gpu->TIStream>>>();
  LAUNCHERROR("kCorrectTIForce_kernel");
}

//---------------------------------------------------------------------------------------------
// ik_CalculateTIKineticEnergy:
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                TI stream parameters
//   c_ave:
//---------------------------------------------------------------------------------------------
void ik_CalculateTIKineticEnergy(gpuContext gpu, double c_ave)
{
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));
  unsigned blocksToUse = 1;
  unsigned threadsPerBlock = 128;

  kgCalculateTIKineticEnergy_kernel<<<blocksToUse, threadsPerBlock,
                                      0, gpu->TIStream>>>(useImage, c_ave);
  LAUNCHERROR("kCalculateTIKineticEnergy_kernel");
}

//---------------------------------------------------------------------------------------------
// ik_PrepareChargeGrid:
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                atom partial charge data and mesh parameters
//   ti_mode:     thermodynamic integration mode, taken from mdin
//---------------------------------------------------------------------------------------------
void ik_PrepareChargeGrid(gpuContext gpu, int ti_mode, int TIRegion)
{
  if (!gpu->bCalculateLocalForces || TIRegion > 0) {
    kPMEClearChargeGridBuffer(gpu);
  }
  kPMEFillChargeGridBuffer(gpu);
  kPMEReduceChargeGridBuffer(gpu);
}

//---------------------------------------------------------------------------------------------
// ik_ZeroTICharge:
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                CUDA stream and launch parameters
//   mode:
//---------------------------------------------------------------------------------------------
void ik_ZeroTICharge(gpuContext gpu, unsigned int mode)
{
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));

  int nterms = gpu->sim.numberTIAtoms;
  if (nterms <= 0) {
    return;
  }
  unsigned threadsPerBlock = (PASCAL || VOLTA || AMPERE) ? 512 : 128;
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);
  unsigned factor = 1;
  unsigned blocksToUse = min((nterms / threadsPerBlock) + 1, gpu->blocks*factor);

  kgZeroTICharges_kernel<<<blocksToUse, threadsPerBlock>>>(useImage, mode);
  LAUNCHERROR("kgZeroTICharges_kernel");
}

//---------------------------------------------------------------------------------------------
// ik_ScaleRECharge:
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                CUDA stream and launch parameters
//   mode:
//---------------------------------------------------------------------------------------------
void ik_ScaleRECharge(gpuContext gpu, unsigned int mode)
{
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));

  int nterms = gpu->sim.numberREAFAtoms;
  if (nterms <= 0) {
    return;
  }
  unsigned threadsPerBlock = (PASCAL || VOLTA || AMPERE) ? 512 : 128;
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);
  unsigned factor = 1;
  unsigned blocksToUse = min((nterms / threadsPerBlock) + 1, gpu->blocks * factor);

  kgScaleRECharges_kernel<<<blocksToUse, threadsPerBlock>>>(useImage, mode);
  LAUNCHERROR("kgScaleRECharges_kernel");
}

//---------------------------------------------------------------------------------------------
// ik_ZeroTIAtomForce:
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                CUDA stream and launch parameters
//   mode:
//---------------------------------------------------------------------------------------------
void ik_ZeroTIAtomForce(gpuContext gpu, unsigned int mode)
{
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));
  unsigned blocksToUse = 1;
  unsigned threadsPerBlock = 128;
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);
  kgZeroTIAtomForce_kernel<<<blocksToUse, threadsPerBlock,
                             0, gpu->TIStream>>>(useImage, mode);
  LAUNCHERROR("kZeroTIAtomForce_kernel");
}

//---------------------------------------------------------------------------------------------
// ik_ClearTIEnergyForce_kernel:
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                CUDA stream and launch parameters
//---------------------------------------------------------------------------------------------
void ik_ClearTIEnergyForce_kernel(gpuContext gpu)
{
  int nterms = gpu->sim.stride3 * ((gpu->sim.needVirial) ? 10 : 5);
  if (nterms <= 0) {
    return;
  }
  unsigned threadsPerBlock = (PASCAL || VOLTA || AMPERE) ? 64 : 1024;
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);
  unsigned factor = 1;
  unsigned blocksToUse = min((nterms / threadsPerBlock) + 1, gpu->blocks*factor);

  kgClearTIForce_kernel<<<blocksToUse, threadsPerBlock, 0, gpu->TIStream>>>();
  LAUNCHERROR("kClearTIForce_kernel");
  kgClearTIPotEnergy_kernel<<<1, GRID, 0, gpu->TIStream>>>();
  LAUNCHERROR("kClearTIPotEnergy_kernel");

  if (gpu->sim.nMBARStates > 0) {
    kgClearMBAR_kernel<<<1, GRID, 0, gpu->TIStream>>>();
    LAUNCHERROR("kgClearMBAR_kernel");
  }

  if (gpu->multiStream) {
    cudaEventRecord(gpu->event_TIClear, gpu->TIStream);
  }
}

//---------------------------------------------------------------------------------------------
// ik_ClearTIEnergyForce_gamd_kernel:
//
// Arguments:
//   gpu:         overarching data structure containing simulation information, here used for
//                CUDA stream and launch parameters
//---------------------------------------------------------------------------------------------
void ik_ClearTIEnergyForce_gamd_kernel(gpuContext gpu)
{
  int nterms = gpu->sim.stride3 * ((gpu->sim.needVirial) ? 14 : 7);
  if (nterms <= 0) {
    return;
  }
  unsigned threadsPerBlock = (PASCAL || VOLTA || AMPERE) ? 64 : 1024;
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);
  unsigned factor = 1;
  unsigned blocksToUse = min((nterms / threadsPerBlock) + 1, gpu->blocks*factor);

  kgClearTIForce_gamd_kernel<<<blocksToUse, threadsPerBlock, 0, gpu->TIStream>>>();
  LAUNCHERROR("kClearTIForce_kernel");
  kgClearTIPotEnergy_kernel<<<1, GRID, 0, gpu->TIStream>>>();
  LAUNCHERROR("kClearTIPotEnergy_kernel");

  if (gpu->sim.nMBARStates > 0) {
    kgClearMBAR_kernel<<<1, GRID, 0, gpu->TIStream>>>();
    LAUNCHERROR("kgClearMBAR_kernel");
  }


  if (gpu->multiStream) {
    cudaEventRecord(gpu->event_TIClear, gpu->TIStream);
  }
}

//---------------------------------------------------------------------------------------------
// ik_SyncVector_kernel:
//
// Arguments:
//   gpu:           overarching data structure containing simulation information, here used for
//                  neighbor list flags and kernel launch parameters
//   mode:
//   combinedMode:
//---------------------------------------------------------------------------------------------
void ik_SyncVector_kernel(gpuContext gpu, unsigned int mode, int combinedMode)
{
  int nterms = gpu->sim.numberTICommonPairs;
  if (nterms <= 0) {
    return;
  }
  int core = 128;
  bool useImage = (gpu->bNeighborList && (gpu->pbImageIndex != NULL));

  unsigned blocksToUse = min((nterms-1) / core + 1, gpu->blocks);
  unsigned threadsPerBlock = core;
  kgSyncVector_kernel<<<blocksToUse, threadsPerBlock>>>(useImage, mode, combinedMode);
  LAUNCHERROR("kgSyncVector_kernel");
}

//---------------------------------------------------------------------------------------------
// ik_BuildTINBList:
//
// Arguments:
//   gpu:           overarching data structure containing simulation information, here used for
//                  neighbor list flags and kernel launch parameters
//---------------------------------------------------------------------------------------------
void ik_BuildTINBList(gpuContext gpu)
{
  if (gpu->sim.numberTIAtoms == 0) {
    return;
  }
  int nterms = gpu->sim.atoms;

  unsigned threadsPerBlock =  (isDPFP) ? 128 : ((PASCAL || VOLTA || AMPERE) ? 256 : 768);
  unsigned factor = 1;
  unsigned blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock ) + 1,
                                                      gpu->blocks * factor);

  kgBuildSpecial2RestNBPreList_kernel<<<blocksToUse, threadsPerBlock,
                                        0, gpu->TIStream>>>(GTI_NB::TI);
  LAUNCHERROR("kgBuildSpecial2RestNBPreList");

  nterms = gpu->sim.atoms;

  threadsPerBlock = (isDPFP) ? 128 : 512;  // Tuned w/ M2000M
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);
  factor = (PASCAL || VOLTA || AMPERE) ? 1 : 2;
  blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock) + 1,
                                             gpu->blocks * factor);

  blocksToUse = gpu->sim.numberUniqueTIAtoms;
  kgBuildSpecial2RestNBList_kernel<<<blocksToUse, threadsPerBlock,
                                     0, gpu->TIStream>>>(GTI_NB::TI);
  LAUNCHERROR("kgBuildSpecial2RestNBList");

  nterms = gpu->sim.numberTIAtoms * 400;
  threadsPerBlock = (isDPFP) ? 128 : ((PASCAL || VOLTA || AMPERE) ? 64 : 1024);
  factor = 1;
  blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock) + 1,
                                             gpu->blocks*factor);

  kgBuildTI2TINBList_kernel <<<blocksToUse, threadsPerBlock, 0, gpu->TIStream >>> ();
  LAUNCHERROR("kgBuildTI2TINBList");

  kgTINBListFillAttribute_kernel<<<blocksToUse, threadsPerBlock,
                                             0, gpu->TIStream>>>();
  LAUNCHERROR("kgTINBListFillAttribute_kernel");

}

//---------------------------------------------------------------------------------------------
// ik_BuildTINBList_gamd:
//
// Arguments:
//   gpu:           overarching data structure containing simulation information, here used for
//                  neighbor list flags and kernel launch parameters
//---------------------------------------------------------------------------------------------
void ik_BuildTINBList_gamd(gpuContext gpu)
{
  if (gpu->sim.numberTIAtoms == 0) {
    return;
  }
  int nterms = gpu->sim.atoms;

  unsigned threadsPerBlock =  (isDPFP) ? 128 : ((PASCAL || VOLTA || AMPERE)? 256 : 768);
  unsigned factor = 1;
  unsigned blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock ) + 1,
                                                      gpu->blocks * factor);

  kgBuildSpecial2RestNBPreList_gamd_kernel<<<blocksToUse, threadsPerBlock,
                                        0, gpu->TIStream>>>(GTI_NB::TI);
  LAUNCHERROR("kgBuildSpecial2RestNBPreList_gamd");

  nterms = gpu->sim.atoms;

  threadsPerBlock = (isDPFP) ? 128: 512;  // Tuned w/ M2000M
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);
  factor = (PASCAL || VOLTA || AMPERE) ? 1 : 2;
  blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock) + 1,
                                             gpu->blocks * factor);

  kgBuildSpecial2RestNBList_gamd_kernel<<<blocksToUse, threadsPerBlock,
                                     0, gpu->TIStream>>>(GTI_NB::TI);
  LAUNCHERROR("kgBuildSpecial2RestNBList_gamd");

  nterms = gpu->sim.numberTIAtoms * 400;
  threadsPerBlock = (isDPFP) ? 128 : ((PASCAL || VOLTA || AMPERE) ? 64 : 1024);
  factor = 1;
  blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock) + 1,
                                             gpu->blocks*factor);

  kgBuildTI2RestNBListFillAttribute_kernel<<<blocksToUse, threadsPerBlock,
                                             0, gpu->TIStream>>>();
  LAUNCHERROR("kgBuildTI2RestNBListFillAttribute");

  // kgBuildTI2TINBList_gamd_kernel<<<blocksToUse, threadsPerBlock, 0, gpu->TIStream>>>();
//  kgBuildTI2TINBList_kernel<<<blocksToUse, threadsPerBlock, 0, gpu->TIStream>>>();
  kgBuildTI2TINBList_gamd2_kernel<<<blocksToUse, threadsPerBlock, 0, gpu->TIStream>>>();
  LAUNCHERROR("kgBuildTI2TINBList");
}

//---------------------------------------------------------------------------------------------
// ik_CalculateTINB:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//   needPotEnergy:
//   needvirial:
//---------------------------------------------------------------------------------------------
void ik_CalculateTINB(gpuContext gpu, bool needPotEnergy, bool needvirial)
{
  int nterms = gpu->sim.numberTIAtoms * gti_simulationConst::MaxNumberNBPerAtom;
  nterms /= 3;
  if (nterms <= 0) {
    return;
  }
  unsigned threadsPerBlock = (isSPFP) ? 768 : 256;

  if (PASCAL || VOLTA || AMPERE) {
    threadsPerBlock = 256;
  }
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);
  unsigned factor = 1;
  unsigned blocksToUse = min((nterms / threadsPerBlock ) + 1, gpu->blocks * factor);

  threadsPerBlock = 128;
  blocksToUse = gpu->sim.numberTIAtoms;
  kgCalculateTINB_kernel<<<blocksToUse, threadsPerBlock,
                           0, gpu->TIStream>>>(needPotEnergy, needvirial);
  LAUNCHERROR("kgCalculateTINB");
}

//---------------------------------------------------------------------------------------------
// ik_CalculateREAFNb:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//   needPotEnergy:
//   needvirial:
//---------------------------------------------------------------------------------------------
void ik_CalculateREAFNb(gpuContext gpu, bool needPotEnergy, bool needvirial)
{
  int nterms = gpu->sim.numberREAFAtoms * gti_simulationConst::MaxNumberNBPerAtom;
  nterms /= 3;
  if (nterms <= 0) {
    return;
  }
  unsigned threadsPerBlock = (isSPFP) ? 768 : 256;

  if (PASCAL || VOLTA || AMPERE) {
    threadsPerBlock = 256;
  }
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);
  unsigned factor = 1;
  unsigned blocksToUse = min((nterms / threadsPerBlock) + 1, gpu->blocks * factor);

  threadsPerBlock = 128;
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);
  blocksToUse = gpu->sim.numberREAFAtoms;

  kgCalculateREAFNb_kernel<<<blocksToUse, threadsPerBlock,
    0, gpu->TIStream>>>(needPotEnergy, needvirial);
  LAUNCHERROR("kgCalculateREAFNb");
}

//---------------------------------------------------------------------------------------------
// ik_CalculateTINB_gamd:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//   needPotEnergy:
//   needvirial:
//---------------------------------------------------------------------------------------------
void ik_CalculateTINB_gamd(gpuContext gpu, bool needPotEnergy, bool needvirial)
{
  int nterms = gpu->sim.numberTIAtoms * 400;
  if (nterms <= 0) {
    return;
  }
  unsigned threadsPerBlock = (isSPFP) ? 768 : 256;

  if (PASCAL || VOLTA || AMPERE) {
    threadsPerBlock = 256;
  }

  threadsPerBlock = 512;
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);
  unsigned factor = 1;
  unsigned blocksToUse = min((nterms / threadsPerBlock ) + 1, gpu->blocks * factor);

  kgCalculateTINB_gamd_kernel<<<blocksToUse, threadsPerBlock,
                           0, gpu->TIStream>>>(needPotEnergy, needvirial);
  LAUNCHERROR("kgCalculateTINB_gamd");
}

//---------------------------------------------------------------------------------------------
// ik_CalculateTINB_ppi_gamd:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//   needPotEnergy:
//   needvirial:
//---------------------------------------------------------------------------------------------
void ik_CalculateTINB_ppi_gamd(gpuContext gpu, bool needPotEnergy, bool needvirial)
{
  int nterms = gpu->sim.numberTIAtoms * 400;
  if (nterms <= 0) {
    return;
  }
  unsigned threadsPerBlock = (isSPFP) ? 768 : 256;

  if (PASCAL || VOLTA || AMPERE) {
    threadsPerBlock = 256;
  }
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);
  unsigned factor = 1;
  unsigned blocksToUse = min((nterms / threadsPerBlock ) + 1, gpu->blocks * factor);

  kgCalculateTINB_ppi_gamd_kernel<<<blocksToUse, threadsPerBlock,
                           0, gpu->TIStream>>>(needPotEnergy, needvirial);
  LAUNCHERROR("kgCalculateTINB_gamd");
}

//---------------------------------------------------------------------------------------------
// ik_CalculateTINB_ppi2_gamd:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//   needPotEnergy:
//   needvirial:
//---------------------------------------------------------------------------------------------
void ik_CalculateTINB_ppi2_gamd(gpuContext gpu, bool needPotEnergy, bool needvirial, int bgpro2atm, int edpro2atm)
{
  int nterms = gpu->sim.numberTIAtoms * 400;
  if (nterms <= 0) {
    return;
  }
  unsigned threadsPerBlock = (isSPFP) ? 768 : 256;

  if (PASCAL || VOLTA || AMPERE) {
    threadsPerBlock = 256;
  }
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);
  unsigned factor = 1;
  unsigned blocksToUse = min((nterms / threadsPerBlock ) + 1, gpu->blocks * factor);

  kgCalculateTINB_ppi2_gamd_kernel<<<blocksToUse, threadsPerBlock,
                           0, gpu->TIStream>>>(needPotEnergy, needvirial,bgpro2atm,edpro2atm);
  LAUNCHERROR("kgCalculateTINB_gamd");
}

//---------------------------------------------------------------------------------------------
// ik_CalculateTIBonded:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//   needPotEnergy:
//   needvirial:
//---------------------------------------------------------------------------------------------
void ik_CalculateTIBonded(gpuContext gpu, bool needPotEnergy, bool needvirial)
{
  unsigned int nwarp_bond = gpu->sim.numberTIBond>0 ? ((gpu->sim.numberTIBond - 1) / GRID + 1) : 0 ;
  unsigned int nwarp_angle = gpu->sim.numberTIAngle>0 ? ((gpu->sim.numberTIAngle - 1) / GRID + 1) : 0;
  unsigned int nwarp_dihedral = gpu->sim.numberTIDihedral>0 ? ((gpu->sim.numberTIDihedral - 1) / GRID + 1) :0 ;

  unsigned int nwarp_distranceRestraint = gpu->sim.numberTINMRDistance>0  ? ((gpu->sim.numberTINMRDistance - 1) / GRID + 1) : 0;
  unsigned int nwarp_angleRestraint = gpu->sim.numberTINMRAngle > 0 ? ((gpu->sim.numberTINMRAngle - 1) / GRID + 1) : 0;
  unsigned int nwarp_dihedralRestraint = gpu->sim.numberTINMRDihedral > 0 ? ((gpu->sim.numberTINMRDihedral - 1) / GRID + 1) : 0;

  int nterms = (nwarp_bond + nwarp_angle + nwarp_dihedral +
    nwarp_distranceRestraint + nwarp_angleRestraint + nwarp_dihedralRestraint) * GRID;
  if (nterms <= 0) {
    return;
  }
  unsigned threadsPerBlock = ((PASCAL || VOLTA || AMPERE) ? 128 : 256);
  unsigned factor = ((PASCAL || VOLTA || AMPERE) ? 2 : 1);
  unsigned blocksToUse = min((nterms / threadsPerBlock) + 1, gpu->blocks*factor);

  kgBondedEnergy_kernel<<<blocksToUse, threadsPerBlock,
                          0, gpu->TIStream>>>(needPotEnergy, needvirial);
  LAUNCHERROR("kgBondedEnergy_kernel");
}


//---------------------------------------------------------------------------------------------
// ik_CalculateREAFBonded:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//   needPotEnergy:
//   needvirial:
//---------------------------------------------------------------------------------------------
void ik_CalculateREAFBonded(gpuContext gpu, bool needPotEnergy, bool needvirial) {

  int nwarp_angle = gpu->sim.numberREAFAngle > 0 ? ((gpu->sim.numberREAFAngle - 1) / 32 + 1) : 0;
  int nwarp_dihedral = gpu->sim.numberREAFDihedral > 0 ? ((gpu->sim.numberREAFDihedral - 1) / 32 + 1) : 0;

  int nterms = ( nwarp_angle + nwarp_dihedral ) * 32;
  if (nterms <= 0) {
    return;
  }
  unsigned threadsPerBlock = ((PASCAL || VOLTA || AMPERE) ? 128 : 256);
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);
  unsigned factor = ((PASCAL || VOLTA || AMPERE) ? 2 : 1);
  unsigned blocksToUse = min((nterms / threadsPerBlock) + 1, gpu->blocks * factor);

  kgREAFBondedEnergy_kernel<<<blocksToUse, threadsPerBlock,
    0, gpu->TIStream>>>(needPotEnergy, needvirial);
  LAUNCHERROR("kgREAFBondedEnergy_kernel");
}

//---------------------------------------------------------------------------------------------
// ik_CalculateTIBonded_ppi:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//   needPotEnergy:
//   needvirial:
//---------------------------------------------------------------------------------------------
void ik_CalculateTIBonded_ppi(gpuContext gpu, bool needPotEnergy, bool needvirial)
{
  unsigned int nwarp_bond = gpu->sim.numberTIBond>0 ? ((gpu->sim.numberTIBond - 1) / GRID + 1) : 0 ;
  unsigned int nwarp_angle = gpu->sim.numberTIAngle>0 ? ((gpu->sim.numberTIAngle - 1) / GRID + 1) : 0;
  unsigned int nwarp_dihedral = gpu->sim.numberTIDihedral>0 ? ((gpu->sim.numberTIDihedral - 1) / GRID + 1) :0 ;

  unsigned int nwarp_distranceRestraint = gpu->sim.numberTINMRDistance>0  ? ((gpu->sim.numberTINMRDistance - 1) / GRID + 1) : 0;
  unsigned int nwarp_angleRestraint = gpu->sim.numberTINMRAngle > 0 ? ((gpu->sim.numberTINMRAngle - 1) / GRID + 1) : 0;
  unsigned int nwarp_dihedralRestraint = gpu->sim.numberTINMRDihedral > 0 ? ((gpu->sim.numberTINMRDihedral - 1) / GRID + 1) : 0;

  int nterms = (nwarp_bond + nwarp_angle + nwarp_dihedral +
    nwarp_distranceRestraint + nwarp_angleRestraint + nwarp_dihedralRestraint) * GRID;
  if (nterms <= 0) {
    return;
  }
  unsigned threadsPerBlock = ((PASCAL || VOLTA || AMPERE) ? 128 : 256);
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);
  unsigned factor = ((PASCAL || VOLTA || AMPERE) ? 2 : 1);
  unsigned blocksToUse = min((nterms / threadsPerBlock) + 1, gpu->blocks*factor);

  kgBondedEnergy_ppi_kernel<<<blocksToUse, threadsPerBlock,
                          0, gpu->TIStream>>>(needPotEnergy, needvirial);
  LAUNCHERROR("kgBondedEnergy_ppi_kernel");
}

#ifndef AMBER_PLATFORM_AMD
//---------------------------------------------------------------------------------------------
// ik_CalculateRMSD:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//   needPotEnergy:
//   needvirial:
//---------------------------------------------------------------------------------------------
void ik_CalculateRMSD(gpuContext gpu, bool needPotEnergy)
{
  
  if (gpu->sim.rmsd_atom_max_count == 0) return;

  //unsigned threadsPerBlock = min(1024, (gpu->sim.rmsd_atom_max_count / 32 + 1) * 32);
  unsigned threadsPerBlock = 32;

  kgRMSDPreparation_kernel << <5,  threadsPerBlock>> > (needPotEnergy);
  LAUNCHERROR("kgRMSDPreparation_kernel");

  cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR; // compute eigenvalues and eigenvectors.
  cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;

  static bool firstTime = true;
  static cusolverDnHandle_t cusolverH = NULL;
  static cudaStream_t stream = NULL;
  if (firstTime) {
    cusolverH = NULL;
    cusolverDnCreate(&cusolverH);
    cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    //cusolverDnSetStream(cusolverH, stream);
    firstTime = false;
  }
  
  const int m = 4;
  const int lda = m;
  unsigned shift = 16;
  int lwork = 65536;


  for (unsigned i = 0; i < gpu->sim.number_rmsd_set; i++) {
#ifdef use_DPFP
    cusolverDnDsyevd(
      cusolverH,
      jobz,
      uplo,
      m,
      (gpu->sim.pMatrix + i * shift),
      lda,
      (gpu->sim.pResult + i * lda),
      (gpu->sim.pWork),
      lwork,
      gpu->sim.pInfo
    );
#else
    cusolverDnSsyevd(
      cusolverH,
      jobz,
      uplo,
      m,
      (gpu->sim.pMatrix + i * shift),
      lda,
      (gpu->sim.pResult + i * lda),
      (gpu->sim.pWork),
      lwork,
      gpu->sim.pInfo
    );
#endif
  }
  kgRMSDEnergyForce_kernel << <5, threadsPerBlock >> > (needPotEnergy);
  LAUNCHERROR("kgRMSDEnergyForce_kernel");
}
#endif

//---------------------------------------------------------------------------------------------
// ik_CalculateTIRestraint:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//   needPotEnergy:
//   needvirial:
//---------------------------------------------------------------------------------------------
void ik_CalculateTIRestraint(gpuContext gpu, bool needPotEnergy, bool needvirial)
{

}

//---------------------------------------------------------------------------------------------
// ik_CalculateTI14NB:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//   needPotEnergy:
//   needvirial:
//---------------------------------------------------------------------------------------------
void ik_CalculateTI14NB(gpuContext gpu, bool needPotEnergy, bool needvirial)
{
  int nterms = gpu->sim.numberTI14NBEntries;
  if (nterms <= 0) {
    return;
  }
  int core = 128;
  unsigned blocksToUse = ((nterms - 1) / core) + 1;
  unsigned threadsPerBlock = core;
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);

  kgCalculateTI14NB_kernel<<<blocksToUse, threadsPerBlock,
                             0, gpu->TIStream>>>(needPotEnergy, needvirial);
  LAUNCHERROR("kgCalculateTI14NB");
}


//---------------------------------------------------------------------------------------------
// ik_CalculateREAF14NB:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//   needPotEnergy:
//   needvirial:
//---------------------------------------------------------------------------------------------
void ik_CalculateREAF14NB(gpuContext gpu, bool needPotEnergy, bool needvirial)
{
  int nterms = gpu->sim.numberREAF14NBEntries;
  if (nterms <= 0) {
    return;
  }
  int core = 128;
  unsigned blocksToUse = ((nterms - 1) / core) + 1;
  unsigned threadsPerBlock = core;

  kgCalculateREAF14NB_kernel<<<blocksToUse, threadsPerBlock,
    0, gpu->TIStream>>>(needPotEnergy, needvirial);
  LAUNCHERROR("kgCalculateREAF14NB");
}

//---------------------------------------------------------------------------------------------
// ik_CalculateTI14NB_gamd:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//   needPotEnergy:
//   needvirial:
//---------------------------------------------------------------------------------------------
void ik_CalculateTI14NB_gamd(gpuContext gpu, bool needPotEnergy, bool needvirial)
{
  int nterms = gpu->sim.numberTI14NBEntries;
  if (nterms <= 0) {
    return;
  }
  int core = 128;
  unsigned blocksToUse = ((nterms - 1) / core) + 1;
  unsigned threadsPerBlock = core;
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);

  kgCalculateTI14NB_gamd_kernel<<<blocksToUse, threadsPerBlock,
                             0, gpu->TIStream>>>(needPotEnergy, needvirial);
  LAUNCHERROR("kgCalculateTI14NB");
}

//---------------------------------------------------------------------------------------------
// ik_CalculateTI14NB_ppi_gamd:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//   needPotEnergy:
//   needvirial:
//---------------------------------------------------------------------------------------------
void ik_CalculateTI14NB_ppi_gamd(gpuContext gpu, bool needPotEnergy, bool needvirial)
{
  int nterms = gpu->sim.numberTI14NBEntries;
  if (nterms <= 0) {
    return;
  }
  int core = 128;
  unsigned blocksToUse = ((nterms - 1) / core) + 1;
  unsigned threadsPerBlock = core;
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);

  kgCalculateTI14NB_ppi_gamd_kernel<<<blocksToUse, threadsPerBlock,
                             0, gpu->TIStream>>>(needPotEnergy, needvirial);
  LAUNCHERROR("kgCalculateTI14NB");
}


//---------------------------------------------------------------------------------------------
// ik_Build1264NBList:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//---------------------------------------------------------------------------------------------
void ik_Build1264NBList(gpuContext gpu)
{
  if (gpu->sim.numberLJ1264Atoms == 0) {
    return;
  }
  int nterms = gpu->sim.atoms / GRID;

  //unsigned threadsPerBlock = (isDPFP) ? 128 : ((PASCAL) ? 256 : 1024);
  unsigned threadsPerBlock =  (isDPFP) ? 128 : ((PASCAL || VOLTA || AMPERE) ? 256 : 768);
  unsigned factor = 1;
  unsigned blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock ) + 1,
                                                      gpu->blocks * factor);

  kgBuildSpecial2RestNBPreList_kernel<<<blocksToUse, threadsPerBlock,
                                        0, gpu->mainStream>>>(GTI_NB::LJ1264);
  LAUNCHERROR("kgBuildSpecial2RestNBPreList");

  nterms = gpu->sim.atoms;


  threadsPerBlock = (isDPFP) ? 128 : 512;  // Tuned w/ M2000M
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);
  factor = (PASCAL || VOLTA || AMPERE) ? 1 : 2;
  blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock) + 1,
                                             gpu->blocks * factor);

  kgBuildSpecial2RestNBList_kernel<<<blocksToUse, threadsPerBlock,
                                     0, gpu->mainStream>>>(GTI_NB::LJ1264);
  LAUNCHERROR("kgBuildSpecial2RestNBList");

  nterms = gpu->sim.numberLJ1264Atoms * 400;
  threadsPerBlock = (isDPFP) ? 128 : ((PASCAL || VOLTA || AMPERE) ? 64 : 1024);
  factor = (PASCAL || VOLTA || AMPERE) ? 4 : 1;
  blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock) + 1,
                                             gpu->blocks * factor);

  kg1264NBListFillAttribute_kernel <<<blocksToUse, threadsPerBlock,
                                          0, gpu->mainStream>>>();
  LAUNCHERROR("kgBuild1264NBListFillAttribute");
}

//---------------------------------------------------------------------------------------------
// ik_Buildp1264NBList:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//---------------------------------------------------------------------------------------------
void ik_Buildp1264NBList(gpuContext gpu)
{
// C4PairwiseCUDA
  if (gpu->sim.numberpLJ1264Atoms == 0) {
    return;
  }
  int nterms = gpu->sim.atoms / 32;

  //unsigned threadsPerBlock = (isDPFP) ? 128 : ((PASCAL) ? 256 : 1024);
  unsigned threadsPerBlock =  (isDPFP) ? 128 : ((PASCAL || VOLTA || AMPERE) ? 256 : 768);
  unsigned factor = 1;
  unsigned blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock ) + 1,
                                                      gpu->blocks * factor);

  kgBuildSpecial2RestNBPreList_kernel<<<blocksToUse, threadsPerBlock,
                                        0, gpu->mainStream>>>(GTI_NB::pLJ1264);
  LAUNCHERROR("kgBuildSpecial2RestNBPreList");

  nterms = gpu->sim.atoms;


  threadsPerBlock = (isDPFP) ? 128 : 512;  // Tuned w/ M2000M
  factor = (PASCAL || VOLTA || AMPERE) ? 1 : 2;
  blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock) + 1,
                                             gpu->blocks * factor);

  kgBuildSpecial2RestNBList_kernel<<<blocksToUse, threadsPerBlock,
                                     0, gpu->mainStream>>>(GTI_NB::pLJ1264);
  LAUNCHERROR("kgBuildSpecial2RestNBList");

  nterms = gpu->sim.numberpLJ1264Atoms * 400;
  threadsPerBlock = (isDPFP) ? 128 : ((PASCAL || VOLTA || AMPERE) ? 64 : 1024);
  factor = (PASCAL || VOLTA || AMPERE) ? 4 : 1;
  blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock) + 1,
                                             gpu->blocks * factor);

  kgp1264NBListFillAttribute_kernel <<<blocksToUse, threadsPerBlock,
                                          0, gpu->mainStream>>>();
  LAUNCHERROR("kgBuild1264NBListFillAttribute");
}

//---------------------------------------------------------------------------------------------
// ik_Build1264p1264NBList:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//---------------------------------------------------------------------------------------------
void ik_Build1264p1264NBList(gpuContext gpu)
{
// C4PairwiseCUDA
  if (gpu->sim.numberLJ1264pLJ1264Atoms == 0) {
    return;
  }
  int nterms = gpu->sim.atoms; 

  //unsigned threadsPerBlock = (isDPFP) ? 128 : ((PASCAL) ? 256 : 1024);
  unsigned threadsPerBlock =  (isDPFP) ? 128 : ((PASCAL || VOLTA || AMPERE) ? 256 : 768);
  unsigned factor = 1;
  unsigned blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock ) + 1,
                                                      gpu->blocks * factor);

  kgBuildSpecial2RestNBPreList_kernel<<<blocksToUse, threadsPerBlock,
                                        0, gpu->mainStream>>>(GTI_NB::LJ1264pLJ1264); //C4PairwiseCUDA2023
  LAUNCHERROR("kgBuildSpecial2RestNBPreList");

  nterms = gpu->sim.atoms;


  threadsPerBlock = (isDPFP) ? 128 : 512;  // Tuned w/ M2000M
  factor = (PASCAL || VOLTA || AMPERE) ? 1 : 2;
  blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock) + 1,
                                             gpu->blocks * factor);

  kgBuildSpecial2RestNBList_kernel<<<blocksToUse, threadsPerBlock,
                                     0, gpu->mainStream>>>(GTI_NB::LJ1264pLJ1264);
  LAUNCHERROR("kgBuildSpecial2RestNBList");

  nterms = gpu->sim.numberLJ1264pLJ1264Atoms * 400; // C4PairwiseCUDA2023
  threadsPerBlock = (isDPFP) ? 128 : ((PASCAL || VOLTA || AMPERE) ? 64 : 1024);
  factor = (PASCAL || VOLTA || AMPERE) ? 4 : 1;
  blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock) + 1,
                                             gpu->blocks * factor);

  kg1264p1264NBListFillAttribute_kernel <<<blocksToUse, threadsPerBlock,
                                          0, gpu->mainStream>>>();
  LAUNCHERROR("kgBuild1264p1264NBListFillAttribute");
}
//---------------------------------------------------------------------------------------------
// ik_BuildREAFNbList:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//---------------------------------------------------------------------------------------------
void ik_BuildREAFNbList(gpuContext gpu)
{
  if (gpu->sim.numberREAFAtoms == 0) {
    return;
  }
  int nterms = gpu->sim.atoms / 32;

  //unsigned threadsPerBlock = (isDPFP) ? 128 : ((PASCAL) ? 256 : 1024);
  unsigned threadsPerBlock = (isDPFP) ? 128 : ((PASCAL || VOLTA || AMPERE) ? 256 : 768);
  unsigned factor = 1;
  unsigned blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock) + 1,
    gpu->blocks * factor);

  kgBuildSpecial2RestNBPreList_kernel<<<blocksToUse, threadsPerBlock,
    0, gpu->mainStream>>>(GTI_NB::REAF);
  LAUNCHERROR("kgBuildSpecial2RestNBPreList");

  nterms = gpu->sim.atoms;

  threadsPerBlock = (isDPFP) ? 128 : 512;  // Tuned w/ M2000M
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);
  factor = (PASCAL || VOLTA || AMPERE) ? 1 : 2;
  blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock) + 1,
    gpu->blocks * factor);

  kgBuildSpecial2RestNBList_kernel<<<blocksToUse, threadsPerBlock,
    0, gpu->mainStream>>>(GTI_NB::REAF);
  LAUNCHERROR("kgBuildSpecial2RestNBList");

  nterms = gpu->sim.numberREAFAtoms * 400;
  threadsPerBlock = (isDPFP) ? 128 : ((PASCAL || VOLTA || AMPERE) ? 64 : 1024);
  factor = (PASCAL || VOLTA || AMPERE) ? 4 : 1;
  blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock) + 1,
    gpu->blocks * factor);

  kgREAFNbListFillAttribute_kernel<<<blocksToUse, threadsPerBlock,
    0, gpu->mainStream>>>();
  LAUNCHERROR("kgBuildREAFNbListFillAttribute");
}

//---------------------------------------------------------------------------------------------
// ik_Build1264NBList_gamd:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//---------------------------------------------------------------------------------------------
void ik_Build1264NBList_gamd(gpuContext gpu)
{
  if (gpu->sim.numberLJ1264Atoms == 0) {
    return;
  }
  int nterms = gpu->sim.atoms / GRID;

  //unsigned threadsPerBlock = (isDPFP) ? 128 : ((PASCAL) ? 256 : 1024);
  unsigned threadsPerBlock = 768;
  unsigned factor = (PASCAL) ? 1 : 1;
  unsigned blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock) + 1,
                                                      gpu->blocks*factor);

  kgBuildSpecial2RestNBPreList_gamd_kernel<<<blocksToUse, threadsPerBlock,
                                        0, gpu->mainStream>>>(GTI_NB::LJ1264);
  LAUNCHERROR("kgBuildSpecial2RestNBPreList_gamd");

  nterms = gpu->sim.atoms;

  //threadsPerBlock = (isDPFP) ? 128 : ((PASCAL) ? 256 : 1024);
  threadsPerBlock = 512;  // Tuned w/ M2000M
  factor = (PASCAL) ? 4 : 1;
  blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock) + 1,
                                             gpu->blocks*factor);

  kgBuildSpecial2RestNBList_gamd_kernel<<<blocksToUse, threadsPerBlock,
                                     0, gpu->mainStream>>>(GTI_NB::LJ1264);
  LAUNCHERROR("kgBuildSpecial2RestNBList_gamd");

  nterms = gpu->sim.numberLJ1264Atoms * 400;
  threadsPerBlock = (isDPFP) ? 128 : ((PASCAL) ? 64 : 1024);
  threadsPerBlock = min(threadsPerBlock, MAX_THREADS_PER_BLOCK);
  factor = (PASCAL) ? 4 : 1;
  blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock) + 1,
                                             gpu->blocks * factor);

  LAUNCHERROR("kgBuild1264NBListFillAttribute");
}

//---------------------------------------------------------------------------------------------
// ik_Buildp1264NBList_gamd:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//---------------------------------------------------------------------------------------------
void ik_Buildp1264NBList_gamd(gpuContext gpu)
{
// C4PairwiseCUDA
  if (gpu->sim.numberpLJ1264Atoms == 0) {
    return;
  }
  int nterms = gpu->sim.atoms / 32;

  //unsigned threadsPerBlock = (isDPFP) ? 128 : ((PASCAL) ? 256 : 1024);
  unsigned threadsPerBlock = 768;
  unsigned factor = (PASCAL) ? 1 : 1;
  unsigned blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock) + 1,
                                                      gpu->blocks*factor);

  kgBuildSpecial2RestNBPreList_gamd_kernel<<<blocksToUse, threadsPerBlock,
                                        0, gpu->mainStream>>>(GTI_NB::pLJ1264);
  LAUNCHERROR("kgBuildSpecial2RestNBPreList_gamd");

  nterms = gpu->sim.atoms;

  //threadsPerBlock = (isDPFP) ? 128 : ((PASCAL) ? 256 : 1024);
  threadsPerBlock = 512;  // Tuned w/ M2000M
  factor = (PASCAL) ? 4 : 1;
  blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock) + 1,
                                             gpu->blocks*factor);

  kgBuildSpecial2RestNBList_gamd_kernel<<<blocksToUse, threadsPerBlock,
                                     0, gpu->mainStream>>>(GTI_NB::pLJ1264);
  LAUNCHERROR("kgBuildSpecial2RestNBList_gamd");

  nterms = gpu->sim.numberpLJ1264Atoms * 400;
  threadsPerBlock = (isDPFP) ? 128 : ((PASCAL) ? 64 : 1024);
  factor = (PASCAL) ? 4 : 1;
  blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock) + 1,
                                             gpu->blocks * factor);

  LAUNCHERROR("kgBuild1264NBListFillAttribute");
}

//---------------------------------------------------------------------------------------------
// ik_Build1264p1264NBList_gamd:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//---------------------------------------------------------------------------------------------
void ik_Build1264p1264NBList_gamd(gpuContext gpu)
{
// C4PairwiseCUDA
  if (gpu->sim.numberLJ1264pLJ1264Atoms == 0) {
    return;
  }
  int nterms = gpu->sim.atoms / 32;

  //unsigned threadsPerBlock = (isDPFP) ? 128 : ((PASCAL) ? 256 : 1024);
  unsigned threadsPerBlock = 768;
  unsigned factor = (PASCAL) ? 1 : 1;
  unsigned blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock) + 1,
                                                      gpu->blocks*factor);

  kgBuildSpecial2RestNBPreList_gamd_kernel<<<blocksToUse, threadsPerBlock,
                                        0, gpu->mainStream>>>(GTI_NB::pLJ1264);
  LAUNCHERROR("kgBuildSpecial2RestNBPreList_gamd");

  nterms = gpu->sim.atoms;

  //threadsPerBlock = (isDPFP) ? 128 : ((PASCAL) ? 256 : 1024);
  threadsPerBlock = 512;  // Tuned w/ M2000M
  factor = (PASCAL) ? 4 : 1;
  blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock) + 1,
                                             gpu->blocks*factor);

  kgBuildSpecial2RestNBList_gamd_kernel<<<blocksToUse, threadsPerBlock,
                                     0, gpu->mainStream>>>(GTI_NB::pLJ1264);
  LAUNCHERROR("kgBuildSpecial2RestNBList_gamd");

  nterms = gpu->sim.numberLJ1264pLJ1264Atoms * 400;
  threadsPerBlock = (isDPFP) ? 128 : ((PASCAL) ? 64 : 1024);
  factor = (PASCAL) ? 4 : 1;
  blocksToUse = (isDPFP) ? gpu->blocks : min((nterms / threadsPerBlock) + 1,
                                             gpu->blocks * factor);

  LAUNCHERROR("kgBuild1264p1264NBListFillAttribute");
}

//---------------------------------------------------------------------------------------------
// ik_Calculate1264NB:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//   needPotEnergy:
//   needvirial:
//---------------------------------------------------------------------------------------------
void ik_Calculate1264NB(gpuContext gpu, bool needPotEnergy, bool needvirial)
{
  int blocksToUse = gpu->sim.numberLJ1264Atoms;
  if (blocksToUse <= 0) {
    return;
  }
  unsigned threadsPerBlock = 256;
  
  kgCalculate1264NB_kernel<<<blocksToUse, threadsPerBlock,
                             0, gpu->mainStream>>>(needPotEnergy, needvirial);
  LAUNCHERROR("kgCalculateTINB");
}

//---------------------------------------------------------------------------------------------
// ik_Calculatep1264NB:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//   needPotEnergy:
//   needvirial:
//---------------------------------------------------------------------------------------------
void ik_Calculatep1264NB(gpuContext gpu, bool needPotEnergy, bool needvirial)
{
// C4PairwiseCUDA
  int blocksToUse = gpu->sim.numberpLJ1264Atoms;
  if (blocksToUse <= 0) {
    return;
  }
  unsigned threadsPerBlock = 256;

  kgCalculatep1264NB_kernel<<<blocksToUse, threadsPerBlock,
                             0, gpu->mainStream>>>(needPotEnergy, needvirial);
  LAUNCHERROR("kgCalculateTINB");
}

//---------------------------------------------------------------------------------------------
// ik_Calculate1264p1264NB:
//
// Arguments:
//   gpu:            overarching data structure containing simulation information, here used
//                   for stream directions and kernel launch parameters
//   needPotEnergy:
//   needvirial:
//---------------------------------------------------------------------------------------------
void ik_Calculate1264p1264NB(gpuContext gpu, bool needPotEnergy, bool needvirial)
{
// C4PairwiseCUDA
  int blocksToUse = gpu->sim.numberLJ1264pLJ1264Atoms; 
  if (blocksToUse <= 0) {
    return;
  }
  unsigned threadsPerBlock = 256;

  kgCalculate1264p1264NB_kernel<<<blocksToUse, threadsPerBlock,
                             0, gpu->mainStream>>>(needPotEnergy, needvirial);
  LAUNCHERROR("kgCalculateTINB");
}

#endif
