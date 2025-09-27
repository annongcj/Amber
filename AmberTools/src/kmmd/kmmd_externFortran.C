#ifndef _KMMD_INIT
#define _KMMD_INIT

#include "kmmd_context.h"

#include "hash.h"
#include "scanDB.h"

#include "pdb_io.h"
#include "kmmd_forcecalc.h"

#include "stdio.h"

//what are we interfacing with? Maybe a #ifdef CUDA in here.
typedef double EXTERN_CRD_T;
typedef double EXTERN_FRC_T;
typedef double EXTERN_ENE_T;



//static context to be allocated at initialisation,
//calling program never sees it.
kmmdHostContext *KMMD;


/*trailing underscore for import to fortran.. */
extern "C" void kmmd_init_( char* kmmd_inpFileName ){
    /* this function to be called to setup an KMMD calculation for embedding in pmemd. */

    //start sigma, default is 0.1
    KMMD = new kmmdHostContext( kmmd_inpFileName );
    printf("Initialised KMMD, was called from fortran code: start sigma2=%.4e\n", KMMD->sigma2);
    
}


/*trailing underscore for import to fortran.. */
extern "C" void kmmd_frccalc_( EXTERN_CRD_T *crd, EXTERN_FRC_T *frc, EXTERN_ENE_T *ene ){
    
    KMMDEneAcc_t point_ene, dVdL;
    
   // fprintf(stderr, "CRD size extern: %i   KMMD: %i\n", (int)sizeof(EXTERN_CRD_T), (int)sizeof(KMMDCrd_t));
   // fprintf(stderr, "%i training points of %i dihedrals\n", KMMD->n_trainpts, KMMD->n_dh);
    
    dVdL = 0.0;
    
#ifdef KMMD_CUDA
    point_ene   = point_ene_host( KMMD->dh_trainpts->_pSysData,
                                  KMMD->preWeights->_pSysData, KMMD->dweights_dlambda_scale->_pSysData, &dVdL,
                                  KMMD->trainenes->_pSysData,
                                  KMMD->dh_atoms->_pSysData,    
                                  crd,
                                  frc,
                                  KMMD->n_dh,
                                  KMMD->n_trainpts,
                                  KMMD->sigma2,
                                  (char *)"test_Cuda" );
#else
    point_ene   = point_ene_host( KMMD->dh_trainpts,
                                  KMMD->preWeights, KMMD->dweights_dlambda_scale, &dVdL,
                                  KMMD->trainenes,
                                  KMMD->dh_atoms, 
                                  crd,
                                  frc,
                                  KMMD->n_dh,
                                  KMMD->n_trainpts,
                                  KMMD->sigma2,
                                  (char *)"test_Ext" );
                                  
#endif
    //fprintf(stderr, "KMMD returned point_ene %e  , adding to %e\n", point_ene, *ene);
                                  
   *ene += (EXTERN_ENE_T)point_ene;

    //feeding TI energy gradients back into amber is complicated, just write and flush from here:
    if( KMMD->dweights_dlambda_scale ){
        KMMD->acc_dvdl           += dVdL;
        KMMD->acc_dvdl_stepcount += 1;
        if ( KMMD->acc_dvdl_stepcount % KMMD->acc_dvdl_steps == 0 ){
            fprintf(KMMD->dvdl_logfile, "%.12f\n", KMMD->acc_dvdl / KMMD->acc_dvdl_stepcount);
            fflush(KMMD->dvdl_logfile);
            KMMD->acc_dvdl           = 0.0;
            KMMD->acc_dvdl_stepcount = 0;
        }
    }
}



#endif
