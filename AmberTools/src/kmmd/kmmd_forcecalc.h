#include "kmmd_context.h"

void measure_dh(  KMMDCrd_t* dh_values,  int* dh_atoms, KMMDCrd_t* crd, int n_dh );
void measure_dof( KMMDCrd_t* dof_values, int* dh_atoms, KMMDCrd_t* crd, int n_dh );
void measure_dof_host( KMMDCrd_t* dof_values, int* dh_atoms, KMMDCrd_t* crd, int n_dh );
void map_dh_forces( int* dh_atoms, KMMDCrd_t* dcrd, KMMDFrcAcc_t* dfrc, KMMDCrd_t* crd, KMMDFrcAcc_t* frc, int n_dh );

 KMMDEneAcc_t KMMD_point_ene_host( kmmdHostContext *KMMD,  KMMDCrd_t *crd,  KMMDCrd_t *frcq);

 KMMDEneAcc_t point_ene_host( KMMDCrd_t     *dof_values,
                                      KMMDCrd_t     *trainpt_preWeights,
                                      KMMDCrd_t     *trainpt_dWeights_dLambda,
                                      KMMDEneAcc_t  *dVdL,
                                      KMMDEneTrn_t  *enes,
                                      int           *dh_atoms,
                                      KMMDCrd_t     *crd,
                                      KMMDFrcAcc_t  *frc,
                                      int            n_dh,
                                      int            n_trainpts,
                                      KMMDCrd_t      sigma2,
                                      char          *calc_id );




 void saveVec_host(KMMDCrd_t **data, int n_pts, const char *filename); //debug function.

//exported for debug purposes, shouldn't be considered part of API
 KMMDEneAcc_t point_dx_vecs_host(KMMDCrd_t *dof_values, 
                                    KMMDCrd_t      *preWeights, 
                                    KMMDCrd_t      *test_dof,
                                    KMMDCrd_t      *dof_deltas,
                                    KMMDEneAcc_t   *trainVec_weights,
                                    KMMDCrd_t       sigma2,
                                    int             n_dh,
                                    int             n_trainpts);

//exported for debug purposes, shouldn't be considered part of API

 KMMDEneAcc_t reduce_dof_forces_host(KMMDCrd_t      *dof_deltas,
                                             KMMDCrd_t      *trainVec_weights, 
                                             KMMDCrd_t      *trainVec_dWeights_dLambda, 
                                             KMMDEneAcc_t   *dEne_dLambda,
                                             KMMDCrd_t      *crds_dof,
                                             KMMDEneAcc_t   *frcs_dof,
                                             KMMDCrd_t       sigma2,
                                             KMMDEneTrn_t   *enes,
                                             KMMDEneAcc_t    sum_weight,
                                             int             n_dh,
                                             int             n_trainpts );

//exported for debug purposes, shouldn't be considered part of API
 void map_dh_forces_host( int*          dh_atoms, 
                                  KMMDCrd_t*    dcrd, 
                                  KMMDFrcAcc_t* dfrc, 
                                  KMMDCrd_t*    crd, 
                                  KMMDFrcAcc_t* frc, 
                                  int           n_dh );
