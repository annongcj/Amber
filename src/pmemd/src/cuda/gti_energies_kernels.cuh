#ifndef GTI_ENERGIES
#define GTI_ENERGIES
#  ifdef GTI

#    include "gpuContext.h"
#    include "gti_kernelAttribute.cuh"

#    include "gti_utils.cuh"
#    if !(defined(__CUDACC_RDC__) || defined(__HIPCC_RDC__))
  namespace GTI_ENERGY_IMPL {
#      include "gti_localCUutils.cuh"
  }
#    endif

// Routines need to be on the device side code

// Real kernels
_kPlainHead_ kgCalculateTIKineticEnergy_kernel(bool useImage, double c_ave);

_kPlainHead_ kgBondedEnergy_kernel(bool energy, bool virial);

_kPlainHead_ kgREAFBondedEnergy_kernel(bool energy, bool virial);

_kPlainHead_ kgRMSDPreparation_kernel(bool energy);

_kPlainHead_ kgRMSDEnergyForce_kernel(bool energy);


_kPlainHead_ kgBondedEnergy_ppi_kernel(bool energy, bool virial);
#  endif  /* GTI */
#endif  /* GTI_ENERGIES */
