#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

#ifndef __GPU_H__
#define __GPU_H__

#include "gpuContext.h"
#include "gputypes.h"
// F95 interface
extern "C" void gpu_init_(void);
extern "C" void gpu_shutdown_(void);
extern "C" void gpu_setup_system_(int *atoms, int *imin, double *tol, int *ntf, int *ntb,
                                  int *ips, int *ntp, int *barostat, int *ntt, int *gbsa,
                                  int *vrand, int *icnstph, int* icnste, int *hybridgb,
				                  int *ti_mode, double* surften, 
                                  int *ischeme, int *ithermostat, int *ntc);
extern "C" void gpu_setup_softcore_ti_(int ti_lst_repacked[], int ti_sc_lst[],
                                       int ti_latm_lst_repacked[], int *ti_latm_cnt);
extern "C" void gpu_set_ti_constants_(double* clambda, double* scalpha, double* scbeta,
                                      int* ifmbar, int* bar_intervall, int* bar_states,
                                      double bar_lambda[][2]);
extern "C" void gpu_ti_dynamic_lambda_(double* clambda);
extern "C" void gpu_upload_gamma_ln_(double* gamma_ln);
extern "C" void gpu_setup_shuttle_info_(int* nshatom, int* infocode, int atom_list[]);
extern "C" void gpu_shuttle_retrieve_data_(double atm_data[][3], int* flag);
extern "C" void gpu_shuttle_post_data_(double atm_data[][3], int* flag);
extern "C" void gpu_neb_exchange_crd_();
extern "C" void gpu_report_neb_ene_(int* master_size, double neb_nrg_all[]);
extern "C" void gpu_set_neb_springs_();
extern "C" void gpu_calculate_neb_frc_(double neb_force[][3]);
extern "C" void gpu_calculate_neb_frc_nstep_(double neb_force[][3]);
extern "C" void gpu_neb_rmsfit_();
extern "C" void gpu_upload_crd_(double atm_crd[][3]);
extern "C" void gpu_upload_crd_gb_cph_(double atm_crd[][3]);
extern "C" void gpu_download_crd_(double atm_crd[][3]);
extern "C" void gpu_upload_charges_(double charge[]);
extern "C" void gpu_upload_charges_gb_cph_(double charge[]);
extern "C" void gpu_upload_charges_pme_cph_(double charge[]);
extern "C" void gpu_download_charges_(double charge[]);
extern "C" void gpu_refresh_charges_(int cit_nb14[][3], double gbl_one_scee[], double qterm[]);
extern "C" void gpu_upload_fs_(double fs[]);
extern "C" void gpu_download_fs_(double fs[]);
extern "C" void gpu_upload_rborn_(double rborn[]);
extern "C" void gpu_download_rborn_(double rborn[]);
extern "C" void gpu_upload_reff_(double reff[]);
extern "C" void gpu_download_reff_(double reff[]);
extern "C" void gpu_upload_frc_(double atm_frc[][3]);
extern "C" void gpu_upload_frc_add_(double atm_frc[][3]);
extern "C" void gpu_download_frc_(double atm_frc[][3]);
extern "C" void gpu_upload_vel_(double atm_vel[][3]);
extern "C" void gpu_download_vel_(double atm_vel[][3]);
extern "C" void gpu_upload_last_vel_(double atm_last_vel[][3]);
extern "C" void gpu_download_last_vel_(double atm_last_vel[][3]);
extern "C" void gpu_clear_vel_();
extern "C" void gpu_bonds_setup_(int* cit_nbona, bond_rec cit_a_bond[], int* cit_nbonh,
                                 bond_rec cit_h_bond[], double gbl_req[], double gbl_rk[]);
extern "C" void gpu_angles_setup_(int* angle_cnt, int* ntheth, angle_rec cit_angle[],
                                  double gbl_teq[], double gbl_tk[]);
extern "C" void gpu_dihedrals_setup_(int* dihed_cnt, int* nhpih, dihed_rec cit_dihed[],
                                     int gbl_ipn[], double gbl_pn[], double gbl_pk[],
                                     double gbl_gamc[], double gbl_gams[]);
extern "C" void gpu_nb14_setup_(int *cit_nb14_cnt, int cit_nb14[][3], double gbl_one_scee[],
                                double gbl_one_scnb[], int *ntypes, int *lj1264, int *plj1264, int iac[],
                                int ico[], double cn1[], double cn2[], double cn6[],
                                int cn7[], double cn8[], int *C4Pairwise, double cn114[], double cn214[]); // C4PairwiseCUDA
extern "C" void gpu_nlinterpolation_setup_();
extern "C" void gpu_molecule_list_setup(int* molecules, listdata_rec listdata[]);
extern "C" void gpu_final_case_inhibitors_(char errmsg[], int *errlen, int *abortsig,
                                           double atm_crd[][3]);
extern "C" void gpu_create_outputbuffers_();
extern "C" void gpu_amd_setup_(int* iamd, int* w_amd, int* iamdlag, int* ntwx,double* EthreshP,
                               double* alphaP, double* EthreshD, double* alphaD,
                               double* EthreshP_w, double* alphaP_w, double* EthreshD_w,
                               double* alphaD_w, double* w_sign, double* temp0);
extern "C" void gpu_calculate_and_apply_amd_weights_(double* pot_ene_tot, double* dih_ene_tot,
                                                     double* amd_ene_tot, double* num_amd_lag);
extern "C" void gpu_download_amd_weights_(double amd_weights_and_energy[]);
extern "C" void gpu_gamd_setup_(int* igamd, int* igamdlag, int* ntwx, double* EthreshP,
                                double* kP, double* EthreshD, double* kD, double* temp0);
extern "C" void gpu_gamd_ppi_setup_(int* igamd, int* igamdlag, int* ntwx, double* EthreshP,
		                double* kP, double* EthreshD, double* kD,double* EthreshB, double* k_B, double* temp0);
extern "C" void gpu_igamd_update_(int* igamd);
extern "C" void gpu_gamd_update_(double* EthreshP, double* kP, double* EthreshD, double* kD);
extern "C" void gpu_gamd_triple_update_(double* EthreshP, double* kP, double* EthreshD, double* kD,double* EthreshB, double* k_B);
extern "C" void gpu_calculate_and_apply_gamd_weights_(double* pot_ene_tot, double* dih_ene_tot,
                                                      double* gamd_ene_tot,
                                                      double* num_gamd_lag);
extern "C" void gpu_calculate_and_apply_gamd_weights_ti_region_(double* pot_ene_nb,
                                                         double* dih_ene_tot,
                                                         double* gamd_ene_tot,
                                                         double* fwgt,
                                                         double* num_gamd_lag);
extern "C" void gpu_calculate_and_apply_gamd_weights_sc_dual_nonbonded_ppi_(double* pot_ene_nb,
                                                         double* dih_ene_tot,
                                                         double* gamd_ene_tot,
                                                         double* num_gamd_lag);
extern "C" void gpu_calculate_and_apply_gamd_weights_sc_triple_ppi_(double* pot_ene_nb,
		                                           double* dih_ene_tot,
							   double* ppi_bond_ene,
						           double* gamd_ene_tot,
							   double* num_gamd_lag);
extern "C" void gpu_calculate_and_apply_gamd_weights_sc_dual_bonded_ppi_(double* pot_ene_nb,
                                                         double* dih_ene_tot,
                                                         double* gamd_ene_tot,
                                                         double* num_gamd_lag);
extern "C" void gpu_calculate_and_apply_gamd_weights_ti_others_(double* pot_ene_nb,
                                                         double* dih_ene_tot,
                                                         double* gamd_ene_tot,
                                                         double* num_gamd_lag);
extern "C" void gpu_calculate_and_apply_gamd_weights_nb_(double* pot_ene_nb,
                                                         double* dih_ene_tot,
                                                         double* gamd_ene_tot,
                                                         double* num_gamd_lag);
extern "C" void gpu_calculate_and_apply_gamd_weights_bond_(double* pot_ene_tot,
                                                         double* dih_ene_tot,
                                                         double* gamd_ene_tot,
                                                         double* num_gamd_lag);
extern "C" void gpu_download_gamd_weights_(double gamd_weights_and_energy[]);
extern "C" void gpu_scaledmd_setup_(int* scaledMD, double* scaledMD_lambda);
extern "C" void gpu_scaledmd_scale_frc_(double* pot_ene_tot);
extern "C" void gpu_download_scaledmd_weights_(double* scaledMD_energy,
                                               double* scaledMD_weight,
                                               double* scaledMD_unscaled_energy);
extern "C" void gpu_download_cell_(double* ucell);
extern "C" void gpu_gbsa3_setup_(double sigma[], double epsilon[], double radius[], double maxsasa[]);
extern "C" void gpu_gbsa_setup_();
extern "C" void gpu_upload_iphmd_(int* iphmd);
extern "C" void gpu_reset_park_(double park[]);
extern "C" void gpu_upload_phmd_(bool* qsetlam, int* ntitr, double* qmass_phmd, double* temp_phmd,
                                 double* phbeta, double para[], double parb[], double park[],
				 double barr[], double ph_theta[], double dph_theta[],
				 double vph_theta[], double charge_phmd[], int grplist[],
				 int sp_grp[], double sp_par[][2], double qstate1[][2],
				 double qstate2[][2], double vstate1[][2], double vstate2[][2],
				 double qstate1_md[][2], double qstate2_md[][2], double parmod[][6],
				 int resbounds1[], int resbounds2[], int baseres[], int* max_numch,
				 int excludeatoms[], int* iphmd, double* solvph, int titratoms[],
				 int* ntitratoms, int numex[], int natex[], int linkgrplist[],
				 bool* use_linked, bool* use_forces);
extern "C" void gpu_download_phmd_theta_(double ph_theta[], double dph_theta[], double vph_theta[],
                                         double dudls[]);
extern "C" void gpu_update_phmd_(double* dtx);
void cpu_kRandomPHMD(gpuContext gpu);
// void kCalculateGBNonbondPHMD(gpuContext gpu);
extern "C" void kRandomPHMD(gpuContext gpu);
void kUpdatePHMD(PMEDouble dtx, gpuContext gpu);
extern "C" void kCalculateSelfPHMD(double vol, double ewaldcof, gpuContext gpu);
void kUpdateForces_PHMD(gpuContext gpu);
extern "C" void kModelForces_PHMD(gpuContext gpu);
extern "C" void kClearDerivs_PHMD(gpuContext gpu);
extern "C" void gpu_gbsa_frc_add_(double atm_frc[][3]);
extern "C" void gpu_force_new_neighborlist_();
extern "C" void gpu_map_honeycomb_(double ucell[][3], double recip[][3], double *es_cutoff,
                                   double *vdw_cutoff, double *skinnb);
extern "C" void gpu_allreduce(GpuBuffer<PMEAccumulator>* pbBuff, int size);
#  ifdef MPI
extern "C" void gpu_calculate_amd_dihedral_weight_(double* totdih);
extern "C" void gpu_calculate_amd_dihedral_energy_(double* totdih);
extern "C" void gpu_calculate_gb_amd_dihedral_energy_(double* totdih);
#  else
extern "C" void gpu_calculate_amd_dihedral_energy_weight_();
extern "C" void gpu_calculate_gb_amd_dihedral_energy_weight_();
#  endif
#  ifdef MPI
extern "C" void gpu_calculate_gamd_dihedral_weight_(double* totdih);
extern "C" void gpu_calculate_gamd_dihedral_energy_(double* totdih);
extern "C" void gpu_calculate_gb_gamd_dihedral_energy_(double* totdih);
#  else
extern "C" void gpu_calculate_gamd_dihedral_energy_weight_();
extern "C" void gpu_calculate_gamd_dihedral_energy_weight_ti_(double* totdih, double *dih_boost);
extern "C" void gpu_calculate_gb_gamd_dihedral_energy_weight_();
#  endif
extern "C" void gpu_set_first_update_atom_(int* first);
extern "C" void gpu_update_natoms(int* natoms, bool* NeighborList);
extern "C" void gpu_ti_exchange_frc_();
extern "C" void gpu_ti_exchange_vel_();
extern "C" void gpu_ti_exchange_crd_();
//MCRES
extern "C" void gpu_steric_grid_(int *atm_cnt, int *maxXVxl, int *maxYVxl, int *maxZVxl,
                                 double *stericGSP, double steric_radius[], int* nmEmpVxls,
                                 int voxel_offset[6]);
extern "C" void gpu_steric_grid_empty_voxels_(int *numvoxels, int listEmptyVoxels[]);
extern "C" void gpu_coarsegrid_setup_(int *watatomcnt, double *ew_coeff, int *mcressize);
extern "C" void gpu_calc_coarsegrid_(int coarseGrid[], double *coarseGSP, int *maxCoarseAtoms,
                                     int *maxCoarseXVxl, int *maxCoarseYVxl,
                                     int *maxcoarseZVxl, int *mcressize, int rndmWatArray[],
                                     double translateArray[][3], double acceptance_array[], int *is_orthog,
                                     int *ti_mode);


extern "C" void gpu_sgld_update_(double* dt, double* temp0, double* gamma_ln, double* gamma_sg, double com0sg[], double com1sg[], double com2sg[]);
extern "C" void gpu_sgld_setup_(int* isgld, int* nsgatom, int sgatoms[], double* dt, double* tsgset, double* sgmsum, double* tsgavg, double* tsgavp, double* sgft, double* sgff, double* sgfg, double* fsgldg, double* tempsg, double sg_x0[][3], double sg_x1[][3], double sg_x2[][3], double sg_r[][3], double sg_fps[], double sg_pps[]);
extern "C" void gpu_sgld_rndf_(double* sgrndf);

void cpu_amrset(int seed);
void cpu_kRandom(gpuContext gpu);

extern "C" void gpu_build_neighbor_list_();

// Local interface
extern "C" void gpuCopyConstants();

// Spline computations
extern "C" void kApplyErfcSplines(PMEFloat* r2, PMEFloat2* results, int *npts,
                                  PMEFloat4* tCdata, int texsize);
extern "C" void kAdjustCoeffsTable(gpuContext gpu);

// Kernel interfaces
extern "C" void kCalculateGBBornRadiiInitKernels(gpuContext gpu);
extern "C" void kCalculateGBNonbondEnergy1InitKernels(gpuContext gpu);
extern "C" void kCalculateGBNonbondEnergy2InitKernels(gpuContext gpu);
extern "C" void kNeighborListInitKernels(gpuContext gpu);
extern "C" void kCalculatePMENonbondEnergyInitKernels(gpuContext gpu);
extern "C" void kPMEInterpolationInitKernels(gpuContext gpu);
extern "C" void kPMEBuildChargeGrid(gpuContext gpu);
extern "C" void kCalculateLocalForcesInitKernels(gpuContext gpu);
extern "C" void kShakeInitKernels(gpuContext gpu);
extern "C" void SetkForcesUpdateSim(gpuContext gpu);
extern "C" void GetkForcesUpdateSim(gpuContext gpu);
extern "C" void SetkXrayUpdateSim(gpuContext gpu);
extern "C" void SetkDataTransferSim(gpuContext gpu);
extern "C" void GetkDataTransferSim(gpuContext gpu);
extern "C" void SetkCalculateNEBForcesSim(gpuContext gpu);
extern "C" void GetkCalculateNEBForcesSim(gpuContext gpu);
extern "C" void SetkCalculateEFieldEnergySim(gpuContext gpu);
extern "C" void GetkCalculateEFieldEnergySim(gpuContext gpu);
extern "C" void kCalculateCoarseGridEnergyInitKernels(gpuContext gpu);
extern "C" void SetkCalculateCoarseGridEnergySim(gpuContext gpu);
extern "C" void GetkCalculateCoarseGridEnergySim(gpuContext gpu);
extern "C" void kMakeStericGrid(gpuContext gpu);
extern "C" void kCountStericGrid(gpuContext gpu);
extern "C" void kGetEmptyStericGridVoxels(gpuContext gpu);
extern "C" void SetkCalculateLocalForcesSim(gpuContext gpu);
extern "C" void GetkCalculateLocalForcesSim(gpuContext gpu);
extern "C" void SetkCalculateGBBornRadiiSim(gpuContext gpu);
extern "C" void GetkCalculateGBBornRadiiSim(gpuContext gpu);
extern "C" void SetkCalculateGBNonbondEnergy1Sim(gpuContext gpu);
extern "C" void GetkCalculateGBNonbondEnergy1Sim(gpuContext gpu);
extern "C" void SetkCalculateGBNonbondEnergy2Sim(gpuContext gpu);
extern "C" void GetkCalculateGBNonbondEnergy2Sim(gpuContext gpu);
extern "C" void SetkShakeSim(gpuContext gpu);
extern "C" void GetkShakeSim(gpuContext gpu);
extern "C" void SetkNeighborListSim(gpuContext gpu);
extern "C" void GetkNeighborListSim(gpuContext gpu);
extern "C" void SetkPMEInterpolationSim(gpuContext gpu);
extern "C" void GetkPMEInterpolationSim(gpuContext gpu);
extern "C" void SetkCalculatePMENonbondEnergySim(gpuContext gpu);
extern "C" void GetkCalculatePMENonbondEnergySim(gpuContext gpu);
extern "C" void SetkCalculatePMENonbondEnergyERFC(double ewcoeff);
extern "C" void SetkCalculateCoarseGridEnergyERFC(double ewcoeff);

// Kernel drivers
extern "C" void SetNLClearForcesKernel(gpuContext gpu);
extern "C" void kClearForces(gpuContext gpu, int totalWarps = 0);
extern "C" void kClearNBForces(gpuContext gpu);
extern "C" void SetNLReduceForcesKernel(gpuContext gpu);
extern "C" void kReduceForces(gpuContext gpu);
extern "C" void kReduceNBForces(gpuContext gpu);
extern "C" void kOrientForces(gpuContext gpu);
extern "C" void kLocalToGlobal(gpuContext gpu);
extern "C" void kExecuteBondWorkUnits(gpuContext gpu, int calcEnergy = 0);
extern "C" void kCalculateLocalForces(gpuContext gpu);
extern "C" void kCalculateLocalEnergy(gpuContext gpu);
extern "C" void kCalculateCHARMMForces(gpuContext gpu);
extern "C" void kCalculateCHARMMEnergy(gpuContext gpu);
extern "C" void kCalculateNMRForces(gpuContext gpu);
extern "C" void kCalculateNMREnergy(gpuContext gpu);
extern "C" void kCalculateGBBornRadii(gpuContext gpu);
extern "C" void kReduceGBBornRadii(gpuContext gpu);
extern "C" void kProcessGBBornRadii(gpuContext gpu);
extern "C" void kClearGBBuffers(gpuContext gpu);
extern "C" void kCalculateGBNonbondEnergy1(gpuContext gpu);
extern "C" void kCalculateGBNonbondForces1(gpuContext gpu);
extern "C" void kReduceGBTemp7(gpuContext gpu);
extern "C" void kReduceMaxsasaEsurf(gpuContext gpu);
extern "C" void kProcessGBTemp7(gpuContext gpu);
extern "C" void kReduceGBTemp7Energy(gpuContext gpu);
extern "C" void kProcessGBTemp7Energy(gpuContext gpu);
extern "C" void kCalculateGBNonbondEnergy2(gpuContext gpu);
extern "C" void kUpdate(gpuContext gpu, PMEDouble dt, PMEDouble temp0, PMEDouble gamma_ln);
extern "C" void kUpdateMiddle(gpuContext gpu, PMEDouble dt, PMEDouble temp0, PMEDouble gamma_ln); //added by zhf
extern "C" void kRefreshCharges(gpuContext gpu);
extern "C" void kRefreshChargesGBCpH(gpuContext gpu);
extern "C" void kRelaxMDUpdate(gpuContext gpu, PMEDouble dt, PMEDouble temp0,
                               PMEDouble gamma_ln);
extern "C" void kShake(gpuContext gpu);
extern "C" void kRattle(gpuContext gpu, PMEFloat dt); // added by zhf
extern "C" void kFastShake(gpuContext gpu);
extern "C" void kCalculateKineticEnergy(gpuContext gpu, PMEFloat c_ave);
extern "C" void kCalculateKineticEnergyAFE(gpuContext gpu, PMEFloat c_ave);
extern "C" void kCalculateEFieldForces(gpuContext gpu, int nstep, double dt);
extern "C" void kCalculateEFieldEnergy(gpuContext gpu, int nstep, double dt);
extern "C" void kCalculateCoarseGridEnergyOrthogTI(gpuContext gpu);
extern "C" void kCalculateCoarseGridEnergyNonOrthogTI(gpuContext gpu);
extern "C" void kCalculateCoarseGridEnergyOrthog(gpuContext gpu);
extern "C" void kCalculateCoarseGridEnergyNonOrthog(gpuContext gpu);
extern "C" void kCalculateCOM(gpuContext gpu);
extern "C" void kCalculateSoluteCOM(gpuContext gpu);
extern "C" void kReduceSoluteCOM(gpuContext gpu);
extern "C" void kClearSoluteCOM(gpuContext gpu);
extern "C" void kCalculateCOMKineticEnergy(gpuContext gpu);
extern "C" void kReduceCOMKineticEnergy(gpuContext gpu);
extern "C" void kCalculateMolecularVirial(gpuContext gpu);
extern "C" void kCalculateSoluteCOM(gpuContext gpu);
extern "C" void kPressureScaleCoordinates(gpuContext gpu);
extern "C" void kPressurScaleConstraints(gpuContext gpu);
extern "C" void kPressureScaleConstraintCoordinates(gpuContext gpu);
extern "C" void kResetVelocities(gpuContext gpu, double temp, double half_dtx);
extern "C" void kClearVelocities(gpuContext gpu);
extern "C" void kRecalculateVelocities(gpuContext gpu, PMEDouble dtx_inv);
extern "C" void kMiddleRecalculateVelocities(gpuContext gpu, PMEDouble dtx_inv); //added by zhf
extern "C" void kScaleVelocities(gpuContext gpu, PMEDouble scale);
extern "C" void kRecenter_Molecule(gpuContext gpu);
extern "C" void kRandom(gpuContext gpu);
extern "C" void kPMEGetGridWeights(gpuContext gpu);
extern "C" void kPMEClearChargeGrid(gpuContext gpu);
extern "C" void kPMEFillChargeGrid(gpuContext gpu);
extern "C" void kPMEConvertChargeGrid(gpuContext gpu);
extern "C" void kPMEClearChargeGridBuffer(gpuContext gpu);
extern "C" void kPMEFillChargeGridBuffer(gpuContext gpu);
extern "C" void kPMEReduceChargeGridBuffer(gpuContext gpu);
extern "C" void kPMEScalarSumRC(gpuContext gpu, PMEDouble vol);
extern "C" void kPMEScalarSumRCEnergy(gpuContext gpu, PMEDouble vol);
extern "C" void kPMEC2CScalarSumRC(gpuContext gpu, PMEDouble vol);
extern "C" void kPMEGradSum(gpuContext gpu);
extern "C" void kPMEForwardFFT(gpuContext gpu);
extern "C" void kPMEBackwardFFT(gpuContext gpu);
extern "C" void kNLResetCounter();
extern "C" void kNLGenerateSpatialHash(gpuContext gpu);
extern "C" void kFilterImage(gpuContext gpu, bool map_ids);
extern "C" void kExpandImage(gpuContext gpu, bool map_ids);
extern "C" void kNLInitRadixSort(gpuContext gpu);
extern "C" void kNLDeleteRadixSort(gpuContext gpu);
extern "C" void kNLRadixSort(gpuContext gpu);
extern "C" void kNLRemapLocalInteractions(gpuContext gpu);
extern "C" void kNLRemapBondWorkUnits(gpuContext gpu);
extern "C" void kNLRemapImage(gpuContext gpu);
extern "C" void kNLCalculateOffsets(gpuContext gpu);
extern "C" void kNLCalculateCellCoordinates(gpuContext gpu);
extern "C" void kNLBuildNeighborList(gpuContext gpu);
extern "C" void kNLClearCellBoundaries(gpuContext gpu);
extern "C" void kNLCalculateCellBoundaries(gpuContext gpu);
extern "C" void kCalculatePMENonbondEnergy(gpuContext gpu);
extern "C" void kCalculatePMENonbondForces(gpuContext gpu);
extern "C" void kCalculatePMELocalForces(gpuContext gpu);
extern "C" void kCalculatePMELocalEnergy(gpuContext gpu);
extern "C" void kCalculateIPSNonbondEnergy(gpuContext gpu);
extern "C" void kCalculateIPSNonbondForces(gpuContext gpu);
extern "C" void kLongRangePairR4(gpuContext gpu);
extern "C" void kLongRangeNrgPairR4(gpuContext gpu);
extern "C" void kNLSkinTest(gpuContext gpu);
extern "C" void kCalculateAMDWeights(gpuContext gpu);
extern "C" void kCalculateAMDWeightAndScaleForces(gpuContext gpu, PMEDouble pot_ene_tot,
                                                  PMEDouble dih_ene_tot, PMEDouble fwgt);
extern "C" void kCalculateAMDDihedralEnergy(gpuContext gpu);
extern "C" void kCalculateGAMDWeights(gpuContext gpu);
extern "C" void kCalculateGAMDWeightAndScaleForces(gpuContext gpu, PMEDouble pot_ene_tot,
                                                   PMEDouble dih_ene_tot, PMEDouble fwgt);
extern "C" void kCalculateGAMDWeightAndScaleForces_nb(gpuContext gpu, PMEDouble pot_ene_nb,
                                                   PMEDouble dih_ene_tot, PMEDouble fwgt);
extern "C" void kCalculateGAMDWeightAndScaleForces_gb(gpuContext gpu, PMEDouble pot_ene_tot,
                                                   PMEDouble dih_ene_tot, PMEDouble fwgt);
extern "C" void kCalculateGAMDWeightAndScaleForces_gb_nb(gpuContext gpu, PMEDouble pot_ene_nb,
                                                   PMEDouble dih_ene_tot, PMEDouble fwgt);

extern "C" void kCalculateGAMDWeightAndScaleForces_bond(gpuContext gpu, PMEDouble pot_ene_nb,
                                                   PMEDouble dih_ene_tot, PMEDouble fwgt);
extern "C" void kCalculateGAMDWeightAndScaleForces_nb_ppi(gpuContext gpu, PMEDouble pot_ene_nb,
                                                   PMEDouble dih_ene_tot, PMEDouble fwgt);
extern "C" void kCalculateGAMDWeightAndScaleForces_ti(gpuContext gpu,
                                                   PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgt);
extern "C" void kCalculateGAMDWeightAndScaleForces_ti_bonded(gpuContext gpu,
                                                   PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgt);
extern "C" void kCalculateGAMDWeightAndScaleForces_ti_nonbonded(gpuContext gpu,
                                                   PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgt);
extern "C" void kCalculateGAMDWeightAndScaleForces_ti_others(gpuContext gpu,
                                                   PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgtd);
extern "C" void kCalculateGAMDWeightAndScaleForces_ti_bonded_others(gpuContext gpu,
                                                   PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgtd);
extern "C" void kCalculateGAMDWeightAndScaleForces_ti_nonbonded_others(gpuContext gpu,
                                                   PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgtd);
extern "C" void kCalculateGAMDWeightAndScaleForces_sc(gpuContext gpu,
                                                   PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgt);
extern "C" void kCalculateGAMDWeightAndScaleForces_sc_ppi(gpuContext gpu,
                                                   PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgt);
extern "C" void kCalculateGAMDWeightAndScaleForces_sc_others_ppi(gpuContext gpu,
                                                   PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgtd);
extern "C" void kCalculateGAMDWeightAndScaleForces_sc_bonded_ppi(gpuContext gpu,
                                                   PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgt);
extern "C" void kCalculateGAMDWeightAndScaleForces_sc_nonbonded(gpuContext gpu,
                                                   PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgt);
extern "C" void kCalculateGAMDWeightAndScaleForces_sc_dual_nonbonded_ppi(gpuContext gpu,
	                                           PMEDouble pot_ene_ti_region,	PMEDouble dih_ene_tot, PMEDouble fwgt, PMEDouble fwgtd);
extern "C" void kCalculateGAMDWeightAndScaleForces_sc_dual_bonded_ppi(gpuContext gpu,
                                                   PMEDouble pot_ene_ti_region, PMEDouble dih_ene_tot, PMEDouble fwgt, PMEDouble fwgtd);
extern "C" void kCalculateGAMDWeightAndScaleForces_sc_triple_ppi(gpuContext gpu, PMEDouble pot_ene_ti_region,
		                        PMEDouble dih_ene_tot,PMEDouble ppi_bond_tot, PMEDouble fwgt, PMEDouble fwgtd,PMEDouble fwgtb);
extern "C" void kCalculateGAMDWeightAndScaleForces_sc_triple2_ppi(gpuContext gpu, PMEDouble pot_ene_ti_region,
		                       PMEDouble dih_ene_tot,PMEDouble ppi_bond_tot, PMEDouble fwgt, PMEDouble fwgtd,PMEDouble fwgtb);
extern "C" void kCalculateGAMDWeightAndScaleForces_sc_triple3_ppi(gpuContext gpu, PMEDouble pot_ene_ti_region,
		                       PMEDouble dih_ene_tot,PMEDouble ppi_bond_tot, PMEDouble fwgt, PMEDouble fwgtd,PMEDouble fwgtb);
extern "C" void kCalculateGAMDWeightAndScaleForces_sc_triple4_ppi(gpuContext gpu, PMEDouble pot_ene_ti_region,
	                      PMEDouble dih_ene_tot,PMEDouble ppi_bond_tot, PMEDouble fwgt, PMEDouble fwgtd,PMEDouble fwgtb);
extern "C" void kCalculateGAMDWeightAndScaleForces_sc_triple5_ppi(gpuContext gpu, PMEDouble pot_ene_ti_region,
		                      PMEDouble dih_ene_tot,PMEDouble ppi_bond_tot, PMEDouble fwgt, PMEDouble fwgtd,PMEDouble fwgtb);
extern "C" void kCalculateGAMDDihedralEnergy(gpuContext gpu);
extern "C" void kScaledMDScaleForces(gpuContext gpu, PMEDouble pot_ene_tot, PMEDouble lambda);
extern "C" void kAFEExchangeFrc(gpuContext gpu);
extern "C" void kAFEExchangeVel(gpuContext gpu);
extern "C" void kAFEExchangeCrd(gpuContext gpu);

// Data transfer
extern "C" void kRetrieveSimData(gpuContext gpu, double atm_data[][3], int modifier);
extern "C" void kPostSimData(gpuContext gpu, double atm_data[][3], int modifier);
extern "C" void DownloadAllCoord(gpuContext gpu, double atm_crd[][3]);
//NEB
extern "C" void kNEBSendRecv(gpuContext gpu, int buff_size);
extern "C" void NEB_report_energy(gpuContext gpu, int* master_size, double neb_nrg_all[]);
extern "C" void kNEBspr(gpuContext gpu);
extern "C" void kNEBfrc(gpuContext gpu, double neb_force[][3]);
extern "C" void kNEBfrc_nstep(gpuContext gpu, double neb_force[][3]);
extern "C" void kFitCOM(gpuContext gpu);

extern "C" void kCheckGpuForces(gpuContext gpu, int frc_chk);
extern "C" void kCheckGpuConsistency(gpuContext gpu, int iter, int cchk, int nchk);
extern "C" void gpu_sumq_(double* sumq, double* sumq2);

#  ifdef MPI
extern "C" void kCopyToAccumulator(gpuContext gpu, PMEAccumulator* p1, PMEAccumulator* p2,
                                   int size);
extern "C" void kAddAccumulators(gpuContext gpu, PMEAccumulator* p1, PMEAccumulator* p2,
                                 int size);
#    ifdef AMBER_PLATFORM_AMD
extern "C" void kPersistentAddAccumulators(gpuContext gpu, PMEAccumulator* pBuff, int size);
#    endif
extern "C" void kReduceForces(gpuContext gpu);
#  endif

// Window into the main library's global gpuContext class instance
gpuContext GetGpuContextPtr();
#endif

void kCalculateSGLDaverages(gpuContext gpu, PMEDouble dt);
extern "C" void kUpdateSGLD(gpuContext gpu, PMEDouble dt, PMEDouble temp0, PMEDouble gamma_ln, PMEDouble gamma_sg, PMEDouble com0sg[], PMEDouble com1sg[], PMEDouble com2sg[]);
void kScaleSGLD(gpuContext gpu, PMEDouble scale, PMEDouble scalsg);
