#ifndef __GTI_F95_H__
#define __GTI_F95_H__

#include "gputypes.h"
#include "gti_simulationConst.h"

// "Global variables"

#  ifdef MPI
static bool useMPI = true;
#  else
static bool useMPI = false;
#  endif
enum RestraintType { RDist, RAngle, RDihedral };

// F95 interface

extern "C" void gti_init_md_parameters_(double* vlimit);

extern "C" void gti_init_ti_parameters_(int ti_latm_lst[][3], int ti_lst[][3], int ti_sc_lst[], int ti_sc_bat_lst[], gti_simulationConst::TypeCtrlVar* ctrlVar);

extern "C" void gti_init_reaf_parameters_(int* reaf_mode, double* reaf_tau, int* addRE, int reaf_atom_list[][2]);

extern "C" void gti_build_nl_list_(int* ti_moode, int* lj1264, int* plj1264, int* read = NULL); // C4PairwiseCUDA

extern "C" void gti_clear_(bool* need_pot_enes, int* ti_mode, double atm_crd[][3]);

extern "C" void gti_update_md_ene_(pme_pot_ene_rec * pEnergy, double enmr[3], double virial[3],
  double ekcmt[3], int* ineb, int* nebfreq, int* nstep);

extern "C" void gti_ele_recip_(bool* needPotEnergy, bool* need_virials, double* vol,
  int* ti_mode, int* reaf_mode);

extern "C" void gti_bonded_setup_(int* bat_type, int* reaf_mode, int* error);

extern "C" void gti_reaf_bonded_setup_(int* error);

extern "C" void gti_restraint_setup_(int* bat_type);

extern "C" void gti_bonded_(bool* needPotEnergy, bool* need_virials, int* ti_moode, int* reaf_mode);

extern "C" void gti_others_(bool* needEnergy, int* nstep, double* dt, double* vol, double* ewaldcof);

extern "C" void gti_kinetic_(double* c_ave);

extern "C" void gti_read_lambda_schedule_(char* fileName);

extern "C" void gti_get_lambda_weights_(int* pType, int* pN, double lambdas[], double weights[], double dWeights[]);

extern "C" void gti_update_lambda_(double* pLambda, int* pIndex, int* k, double glb_weight[2], double glb_dweight[2], double item_weight[2][Schedule::TypeTotal], double item_dweights[2][Schedule::TypeTotal], bool* pUse);

extern "C" void gti_update_lambda_weights_(double* pLambda, int* k, double glb_weight[2], double glb_dweight[2], double item_weight[2][Schedule::TypeTotal], double item_dweights[2][Schedule::TypeTotal], bool* pUse);

extern "C" void gti_reaf_tau_schedule_(char* schFileName);

extern "C" void gti_update_tau_(double* pTau, double reaf_weights[2], double read_item_weights[2][Schedule::TypeTotal], bool* pUse);

extern "C" void gti_update_tau_weights_(double* pTau, double reaf_weights[2], double read_item_weights[2][Schedule::TypeTotal], bool* pUseSP);

extern "C" void gti_setup_update_(double charge[], int* sync_mass, double mass[]);

extern "C" void gti_nb_setup_(int* pNtypes, int iac[], int ico[], double cn1[], double cn2[],
                              double cn6[], int cn7[], double cn8[], int* pC4Pairwise); // C4PairwiseCUDA

extern "C" void gti_ti_nb_setup_(int* pNtypes, int iac[], int ico[], int* pNb14_cnt,
  int cit_nb14[][3], double gbl_one_scee[],
  double gbl_one_scnb[], double cn114[], double cn214[]);

extern "C" void gti_reaf_nb_setup_(int* pNtypes, int iac[], int ico[], int* pNb14_cnt,
  int cit_nb14[][3], double gbl_one_scee[],
  double gbl_one_scnb[], double cn114[], double cn214[]);

extern "C" void gti_ti_nb_gamd_setup_(int* pNtypes, int iac[], int ico[], int* pNb14_cnt,
		                      int cit_nb14[][3], double gbl_one_scee[],
				      double gbl_one_scnb[], double cn114[], double cn214[]);

extern "C" void gti_lj1264_nb_setup_(char atm_igraph[][4], int* ntypes, int atm_iac[],
                                     int ico[], double gbl_cn6[]); 

extern "C" void gti_plj1264_nb_setup_(char atm_igraph[][4], int* ntypes, int atm_iac[],
                                     int ico[], int gbl_cn7[], double gbl_cn8[], int* pC4Pairwise); //C4PairwiseCUDA

extern "C" void gti_lj1264plj1264_nb_setup_(char atm_igraph[][4], int* ntypes, int atm_iac[],
                                     int ico[], double gbl_cn6[], int gbl_cn7[], double gbl_cn8[], int* pC4Pairwise); //C4PairwiseCUDA

extern "C" void gti_pme_setup_();

extern "C" void gti_nb_(bool* needPotEnergy, bool* needvirial, int* ti_mode, int* lj1264, int* plj1264, int* reaf_mode); // C4PairwiseCUDA

extern "C" void gti_finalize_force_virial_(int* numextra, bool* needvirial, int* ti_mode,
                                           int*net_force, double frc[][3]);

extern "C" void gti_get_pot_energy_(double
                                    TIPotEnergy[3][gti_simulationConst::GPUPotEnergyTerms],
                                    int* ti_mode);
extern "C" void gti_get_virial_(double TIWeights[2], double virial[3], double ekcmt[3],
                                int* ti_mode);

extern "C" void gti_get_kin_energy_(double
                                    TIKinEnergy[3][gti_simulationConst::GPUKinEnergyTerms]);

extern "C" void gti_get_mbar_energy_(double MBAREnergy[]);

extern "C" void gti_sync_vector_(int* mode, int* combined);

extern "C" void gti_download_frc_(int* pTIRegion, double atm_frc[][3], bool* needVirialbool,
                                  bool* fresh);

extern "C" void gti_setup_localheating_(double* pTI_tempi);

extern "C" void gti_init_rmsd_(int* pNumber_rmsd_set, int* pRmsd_type, int rmsd_atom_count[],
  int rmsd_ti_region[], int rmsd_atom_list[], double rmsd_ref_crd[], double rmsd_ref_com[], double rmsd_weights[]);

extern "C" void gti_turnoff_localheating_();

extern "C" void gti_update_simulation_const_();

/////   GAMD///
extern "C" void gti_get_virial_gamd_(double TIWeights[2], double virial[3], double ekcmt[3],
  int* ti_mode);

extern "C" void gti_nb_gamd_(bool* needPotEnergy, bool* needvirial, int* ti_mode, int* lj1264, int* plj1264); // C4PairwiseCUDA

extern "C" void gti_finalize_force_virial_gamd_(int* numextra, bool* needvirial, int* ti_mode,
  int* net_force, double frc[][3]);

extern "C" void gti_finalize_force_gamd_(int* numextra, bool* needvirial, int* ti_mode,
  int* net_force, double frc[][3]);

extern "C" void gti_virial_gamd_(int* numextra, bool* needvirial, int* ti_mode,
  int* net_force, double frc[][3]);

extern "C" void gti_init_gamd_ti_parameters_(int ti_latm_lst[][3], int ti_lst[][3],
  int ti_sc_lst[], int ti_sc_bat_lst[],
  double* lambda, unsigned* klambda, double* scalpha,
  double* scbeta, double* gamma, int* pCut,
  int* pAddSC, int* pEleGauss, int* pEleSC, int* pVdwSC, double* pVdwCap);

extern "C" void gti_build_nl_list_gamd_(int* ti_moode, int* lj1264, int* plj1264); // C4PairwiseCUDA

extern "C" void gti_clear_gamd_(bool* need_pot_enes, int* ti_mode, double atm_crd[][3]);

extern "C" void gti_bonded_gamd_setup_(int* pBat_type, int* pError);


#endif
