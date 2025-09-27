#ifndef MDLXRAY_FUNCS
#define MDLXRAY_FUNCS

#include "gpuContext.h"
#include "xrayHostContext.h"

// Fortran interface
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
				  int *bs_mdlidx, double *k_mask, double* scatter_coeffs,
				  int *num_scatter_types, int *scatter_ncoeffs);

extern "C" void gpu_xray_setup_ii_(double* atom_bfactor, double* atom_occupancy,
                                   int* atom_scatter_type);

extern "C" void gpu_xray_get_derivative_(double *xray_e, double *r_work, double *r_free,
                                         int *calcScalingValues, double *Fcalc_scale,
					 double *norm_scale);

extern "C" void gpu_download_sf_(double *r_Fcalc, double *i_Fcalc);

extern "C" void gpu_upload_solvent_contribution_(double *r_SolvFcalc, double *i_SolvFcalc);

// Kernel interface
extern "C" void SetkXrayUpdateCalc(xrayHostContext xrd);
extern "C" void kXrayGetDerivative(gpuContext gpu, xrayHostContext xrd, int calcScalingValues);

#endif
