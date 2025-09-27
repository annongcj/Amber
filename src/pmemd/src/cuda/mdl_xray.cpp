#include "copyright.i"

#include <vector>
#include "gpu.h"
#include "mdl_xray.h"
#include "xrayHostContext.h"
#include "gpuContext.h"

namespace xrayHostImpl {
  static xrayHostContext xrd = NULL;
}

using namespace std;
using namespace xrayHostImpl;

//---------------------------------------------------------------------------------------------
// gpu_xray_setup_i: called from xray_interface.F90 within xray_init() to initialize data
//                   structures for X-ray refinement calculations on the GPU.  This is part
//                   of the initialization.  The rest will come when the gpuContext has been
//                   established.
//
// Arguments:
// 
//---------------------------------------------------------------------------------------------
extern "C" void gpu_xray_setup_i_(int* atom_selection, int *natom, int *nres, double* pbc_box,
                                  double *pbc_alpha, double *pbc_beta, double *pbc_gamma,
                                  int hkl_index[][3], int *num_hkl, int *ivtarget,
				  double* realFObs, double* imagFObs, double* absFObs,
                                  double* sigFObs, double* mSS4, int* test_flag,
				  double *resolution_low, double *resolution_high,
				  double *xray_weight, double *solvent_mask_probe_radius,
				  double *solvent_mask_expand, double *solvent_scale,
				  double *solvent_bfactor, double *bfactor_min,
				  double *bfactor_max, int *bfactor_refinement_interval,
				  int *bs_mdlidx, double *k_mask, double* scatter_coefficients,
				  int *num_scatter_types, int *scatter_ncoeffs)
{
  PRINTMETHOD("gpu_xray_setup_i");
  
  int h, i, j, k, padded_natom, padded_nhkl, padded_ncoeff, nselected;
  double sumFObs, sumXFObs, sumFoFo;

  xrd = theXrayHostContext::GetPointer();

  padded_natom = ((*natom + GRID_BITS_MASK) / GRID) * GRID;
  xrd->pbSelectionMask = new GpuBuffer<int>(padded_natom);
  xrd->pbSelectedAtoms = new GpuBuffer<int>(padded_natom);

  // Perform the compaction that pack_index does for the CPU code.  Because this
  // mask does not change throughout the simulation, it will be made here, once.
  // This will also deprecate the indexing to convert to C++ from Fortran.
  nselected = 0;
  for (i = 0; i < *natom; i++) {
    xrd->pbSelectionMask->_pSysData[i] = atom_selection[i];      
    if (atom_selection[i] == 1) {
      xrd->pbSelectedAtoms->_pSysData[nselected] = i;
      nselected++;
    }
  }
  padded_nhkl = ((*num_hkl + GRID_BITS_MASK) / GRID) * GRID;
  xrd->pbHKLaidx          = new GpuBuffer<int>(padded_nhkl);
  xrd->pbHKLbidx          = new GpuBuffer<int>(padded_nhkl);
  xrd->pbHKLcidx          = new GpuBuffer<int>(padded_nhkl);
  xrd->pbFObs             = new GpuBuffer<PMEFloat2>(padded_nhkl);
  xrd->pbAbsFObs          = new GpuBuffer<PMEFloat>(padded_nhkl);
  xrd->pbSigFObs          = new GpuBuffer<PMEFloat>(padded_nhkl);
  xrd->pbHKLMask          = new GpuBuffer<int>(padded_nhkl);
  xrd->pbMSS4             = new GpuBuffer<PMEFloat>(padded_nhkl);
  xrd->pbSFScaleFactors   = new GpuBuffer<PMEFloat>(4);
  xrd->pbFcalc            = new GpuBuffer<PMEFloat2>(padded_nhkl);
  xrd->pbAbsFcalc         = new GpuBuffer<PMEFloat>(padded_nhkl);
  xrd->pbSolventKMask     = new GpuBuffer<PMEFloat>(padded_nhkl);
  xrd->pbSolventFMask     = new GpuBuffer<PMEFloat2>(padded_nhkl);
  sumFObs  = 0.0;
  sumXFObs = 0.0;
  sumFoFo  = 0.0;
  for (i = 0; i < *num_hkl; i++) {
    xrd->pbHKLaidx->_pSysData[i]         = hkl_index[i][0];
    xrd->pbHKLbidx->_pSysData[i]         = hkl_index[i][1];
    xrd->pbHKLcidx->_pSysData[i]         = hkl_index[i][2];
    xrd->pbFObs->_pSysData[i].x          = realFObs[i];
    xrd->pbFObs->_pSysData[i].y          = imagFObs[i];
    xrd->pbAbsFObs->_pSysData[i]         = absFObs[i];
    xrd->pbSigFObs->_pSysData[i]         = sigFObs[i];
    xrd->pbHKLMask->_pSysData[i]         = test_flag[i];
    xrd->pbMSS4->_pSysData[i]            = mSS4[i];
    if (*bs_mdlidx != 0) {
      xrd->pbSolventKMask->_pSysData[i]    = k_mask[i];
    }
    if (test_flag[i] != 0) {
      sumFObs += absFObs[i];
      sumFoFo += absFObs[i] * absFObs[i];
    }
    else {
      sumXFObs += absFObs[i];
    }
  }
  xrd->pbXrayResults = new GpuBuffer<PMEAccumulator>(GRID);
  xrd->pbReflDeriv   = new GpuBuffer<PMEFloat2>(padded_nhkl);

  // Store the scattering coefficients in a linearized format
  padded_ncoeff = 2 * (*num_scatter_types) * (*scatter_ncoeffs);
  padded_ncoeff = ((padded_ncoeff + GRID_BITS_MASK) / GRID) * GRID;
  xrd->pbScatterCoeffs = new GpuBuffer<PMEFloat>(padded_ncoeff);
  for (i = 0; i < 2 * (*num_scatter_types) * (*scatter_ncoeffs); i++) {
    xrd->pbScatterCoeffs->_pSysData[i] = scatter_coefficients[i];
  }

  // Forces are stored akin to the way they arrays in the rest of pmemd.cuda work
  xrd->pbXrayForces  = new GpuBuffer<PMEAccumulator>(3*padded_natom);

  // Store critical constants and pointers to data in xrd
  xrd->devInfo.nHKL                 = *num_hkl;
  xrd->devInfo.natom                = *natom;
  xrd->devInfo.nres                 = *nres;
  xrd->devInfo.tarfunc              = *ivtarget;
  xrd->devInfo.nScatterTypes        = *num_scatter_types;
  xrd->devInfo.nScatterCoeff        = *scatter_ncoeffs;
  xrd->devInfo.minBFactor           = *bfactor_min;
  xrd->devInfo.maxBFactor           = *bfactor_max;
  xrd->devInfo.BFacRefineInterval   = *bfactor_refinement_interval;
  xrd->devInfo.xrayGeneralWt        = *xray_weight;
  xrd->devInfo.nSelectedAtom        = nselected;
  xrd->devInfo.bsMaskModel          = *bs_mdlidx;
  xrd->devInfo.NormScale            = 1.0 / sumFoFo;
  xrd->devInfo.SumFObsWork          = sumFObs;
  xrd->devInfo.SumFObsFree          = sumXFObs;
  xrd->devInfo.pSelectionMask       = xrd->pbSelectionMask->_pDevData;
  xrd->devInfo.pSelectedAtoms       = xrd->pbSelectedAtoms->_pDevData;
  xrd->devInfo.pScatterCoeffs       = xrd->pbScatterCoeffs->_pDevData;
  xrd->devInfo.pHKLaidx             = xrd->pbHKLaidx->_pDevData;
  xrd->devInfo.pHKLbidx             = xrd->pbHKLbidx->_pDevData;
  xrd->devInfo.pHKLcidx             = xrd->pbHKLcidx->_pDevData;
  xrd->devInfo.pFObs                = xrd->pbFObs->_pDevData;
  xrd->devInfo.pAbsFObs             = xrd->pbAbsFObs->_pDevData;
  xrd->devInfo.pSigFObs             = xrd->pbSigFObs->_pDevData;
  xrd->devInfo.pHKLMask             = xrd->pbHKLMask->_pDevData;
  xrd->devInfo.pMSS4                = xrd->pbMSS4->_pDevData;
  xrd->devInfo.pSFScaleFactors      = xrd->pbSFScaleFactors->_pDevData;
  xrd->devInfo.pFcalc               = xrd->pbFcalc->_pDevData;
  xrd->devInfo.pAbsFcalc            = xrd->pbAbsFcalc->_pDevData;
  xrd->devInfo.pSolventKMask        = xrd->pbSolventKMask->_pDevData;
  xrd->devInfo.pSolventFMask        = xrd->pbSolventFMask->_pDevData;
  xrd->devInfo.pReflDeriv           = xrd->pbReflDeriv->_pDevData;
  xrd->devInfo.pSumFoFc             = &xrd->pbXrayResults->_pDevData[0];
  xrd->devInfo.pSumFoFo             = &xrd->pbXrayResults->_pDevData[1];
  xrd->devInfo.pSumFcFc             = &xrd->pbXrayResults->_pDevData[2];
  xrd->devInfo.pSumFobs             = &xrd->pbXrayResults->_pDevData[3];
  xrd->devInfo.pResidual            = &xrd->pbXrayResults->_pDevData[4];
  xrd->devInfo.pXrayEnergy          = &xrd->pbXrayResults->_pDevData[5];
  xrd->devInfo.pFreeResidual        = &xrd->pbXrayResults->_pDevData[6];
  xrd->devInfo.pForceX              = &xrd->pbXrayForces->_pDevData[               0];
  xrd->devInfo.pForceY              = &xrd->pbXrayForces->_pDevData[    padded_natom];
  xrd->devInfo.pForceZ              = &xrd->pbXrayForces->_pDevData[2 * padded_natom];

  // Upload constants, pointers, and array data
  SetkXrayUpdateCalc(xrd);
  xrd->pbSelectionMask->Upload();
  xrd->pbSelectedAtoms->Upload();
  xrd->pbHKLaidx->Upload();
  xrd->pbHKLbidx->Upload();
  xrd->pbHKLcidx->Upload();
  xrd->pbFObs->Upload();
  xrd->pbAbsFObs->Upload();
  xrd->pbSigFObs->Upload();
  xrd->pbHKLMask->Upload();
  xrd->pbMSS4->Upload();
  xrd->pbSolventKMask->Upload();
  xrd->pbScatterCoeffs->Upload();
}

//---------------------------------------------------------------------------------------------
// gpu_xray_setup_ii: called from pmemd.F90 after the gpuContext has been established.  This
//                    will draw on information in that class instance (more precisely a
//                    pointer it) to finish initializing data for X-ray structure factor
//                    computations.
//
// Arguments:
//---------------------------------------------------------------------------------------------
extern "C" void gpu_xray_setup_ii_(double* atom_bfactor, double* atom_occupancy,
				   int* atom_scatter_type)
{
  PRINTMETHOD("gpu_xray_setup_ii");

  int i, j, natom;
  gpuContext gpuPtr;

  // Get the pointer and then pull data out
  gpuPtr = GetGpuContextPtr();
  natom = gpuPtr->sim.atoms;
  
  // Allocate other arrays to hold information that is now known about the system
  xrd->pbScatterType                = new GpuBuffer<int>(natom);
  xrd->pbAtomBFactor                = new GpuBuffer<PMEFloat>(natom);
  xrd->pbOccupancy                  = new GpuBuffer<PMEFloat>(natom);
  xrd->devInfo.pAtomBFactor         = xrd->pbAtomBFactor->_pDevData;
  xrd->devInfo.pOccupancy           = xrd->pbOccupancy->_pDevData;
  xrd->devInfo.pScatterType         = xrd->pbScatterType->_pDevData;
  
  // Add information read from the prmtop, which should all
  // be in memory by the time this function is called.
  for (i = 0; i < natom; i++) {
    xrd->pbScatterType->_pSysData[i] = atom_scatter_type[i] - 1;
    xrd->pbOccupancy->_pSysData[i]   = atom_occupancy[i];
    xrd->pbAtomBFactor->_pSysData[i] = atom_bfactor[i];
  }

  // Update (upload for a second time) the X-ray device constants, then the arrays
  SetkXrayUpdateCalc(xrd);
  xrd->pbScatterType->Upload();
  xrd->pbAtomBFactor->Upload();
  xrd->pbOccupancy->Upload();
}

//---------------------------------------------------------------------------------------------
// gpu_xray_get_derivative: called from pme_force.F90 in place of the eponymous CPU subroutine.
//
// Arguments:
//   xray_e:  energy of X-ray based restraint penalties
//---------------------------------------------------------------------------------------------
extern "C" void gpu_xray_get_derivative_(double *xray_e, double *r_work, double *r_free,
                                         int *calcScalingValues, double *Fcalc_scale,
					 double *norm_scale)
{
  PRINTMETHOD("gpu_xray_get_derivative");

  gpuContext gpuPtr;  
  
  gpuPtr = GetGpuContextPtr();
  kXrayGetDerivative(gpuPtr, xrd, *calcScalingValues);
  xrd->pbXrayResults->Download();
  xrd->pbSFScaleFactors->Download();
  *xray_e = xrd->pbXrayResults->_pSysData[5];
  *xray_e /= ENERGYSCALE;
  *r_work = xrd->pbXrayResults->_pSysData[4];
  *r_work /= ENERGYSCALE;
  *r_free = xrd->pbXrayResults->_pSysData[6];
  *r_free /= ENERGYSCALE;
  *Fcalc_scale = xrd->pbSFScaleFactors->_pSysData[0];
  *norm_scale = xrd->devInfo.NormScale;
}

//---------------------------------------------------------------------------------------------
// gpu_download_sf: download the calculated structure factors.
//
// Arguments:
//   r_Fcalc:   array of calculated structure factor real components
//   i_Fcalc:   array of calculated structure factor imaginary components
//---------------------------------------------------------------------------------------------
extern "C" void gpu_download_sf_(double *r_Fcalc, double *i_Fcalc)
{
  int i;
  
  PRINTMETHOD("gpu_download_sf");
  xrd->pbFcalc->Download();
  for (i = 0; i < xrd->devInfo.nHKL; i++) {
    r_Fcalc[i] = xrd->pbFcalc->_pSysData[i].x;
    i_Fcalc[i] = xrd->pbFcalc->_pSysData[i].y;
  }
}

//---------------------------------------------------------------------------------------------
// gpu_upload_solvent_contribution: upload the solvent contributions calculated on the host.
//
// Arguments:
//   r_SolvFcalc:  array of calculated structure factor contributions' real components
//   i_SolvFcalc:  array of calculated structure factor contributions' imaginary components
//---------------------------------------------------------------------------------------------
extern "C" void gpu_upload_solvent_contribution_(double *r_SolvFcalc, double *i_SolvFcalc)
{
  int i;

  PRINTMETHOD("gpu_upload_solvent_contribution");
  for (i = 0; i < xrd->devInfo.nHKL; i++) {
    xrd->pbSolventFMask->_pSysData[i].x = r_SolvFcalc[i];
    xrd->pbSolventFMask->_pSysData[i].y = i_SolvFcalc[i];
  }
  xrd->pbSolventFMask->Upload();
}
