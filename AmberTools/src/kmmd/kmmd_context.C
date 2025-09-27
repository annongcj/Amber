#ifndef _KMMD_HOST_CONTEXT_FUNCS
#define _KMMD_HOST_CONTEXT_FUNCS

#include "kmmd_types.h"
#include "kmmd_context.h"

#include "hash.h"
#include "scanDB.h"
#include "DB_vs_forcefield.h"
#include "kmmd_forcecalc.h"

#include "pdb_io.h"

void kmmdHostContext::Init(char *kmmd_inpFileName) {

    //training data
    KMMDCrd_t     *dh_coords,  *coords, *new_crd, *test_dof, *preWeights, *weightGradientsTI, sigma2;
    KMMDFrcTrn_t  *frc;
    KMMDEneTrn_t  *enes, *classEnes, *deltaEnes, point_ene;
    
    //basic system info
    int         N_atoms;
    char      **atnames;
    int        *resids, *dh_atoms;
    int         i, ii, n_db_strucs, n_tracked_dh, d, acc_dvdl_steps;
    //i / o.
    FILE       *output_pdb, *output_denes, *output_at2, *kmmd_inpFile;
     
    //load a pdb structure, identify dihedrals
    //these values in  a json file, eg: "GGCC_setup.json":
    char *DB_fileList_fname;       // = "../../DYADS_RESULTS_DB/DYADS_DB_05_2021";
    char *ref_pdb;//        = "../../DYADS_RESULTS_DB/DYADS_DB_05_2021/GG/GG_nW20.pdb";
    char *dvdl_logfile_name;
     
    /* data structures to record names of dihedrals and ids of defining atoms */
    HashTable  *dh_atnames, *resdel, *use_dh;
    
    //setup all variables needed and check that they were read correctly.
    kmmd_inpFile = fopen(kmmd_inpFileName, "r");
    if( !kmmd_inpFile ){
        fprintf(stderr, "Couldn't open KMMD control file, path: %s\n",
                kmmd_inpFileName);
        exit( 1 );
    }else{
        fprintf(stderr, "Opened KMMD control file, path: %s\n",
                kmmd_inpFileName);
    }
    DB_fileList_fname = NULL;
    ref_pdb    = NULL;
    dh_atnames = NULL;
    resdel     = NULL;
    use_dh     = NULL;
    weightGradientsTI = NULL;
    sigma2     = 0.1; //sensible-ish (conservative) default.
    dvdl_logfile_name = NULL;    
    acc_dvdl_steps    = 1;

    fprintf(stderr, "scanning KMMD input JSON file %s\n", kmmd_inpFileName);
    scan_json_tags( kmmd_inpFile, &DB_fileList_fname, &ref_pdb, &dvdl_logfile_name, &acc_dvdl_steps,
                      &dh_atnames,   &resdel,  &use_dh, &sigma2);
    fclose(kmmd_inpFile);
    
    
    // load reference PDB file for topology.
    N_atoms = get_pdb(ref_pdb, &coords, &atnames, &resids);
    for( i = 0; i < N_atoms; i++ ){
        printf("%i atname: %s x,y,z: %.2f, %.2f, %.2f resids: %i\n",
                    i, atnames[i], coords[i*3], coords[i*3+1], coords[i*3+2], resids[i]);
        if( i >= 40 ){
              if( i == 40 ){printf("    ...etc\n");}
              break;
        }
    }
    
    printf("building some atom indices...\n");

    //list of atoms taking part in each dihedral.
    dh_atoms = assign_dh_atoms(use_dh, resdel, dh_atnames, atnames, resids, N_atoms);
 
    //load training data points from the DB
    n_db_strucs = scan_db((char  *)DB_fileList_fname, &dh_coords, &enes, &classEnes, &preWeights, &weightGradientsTI, NULL, dh_atoms, use_dh->count);
 
    //get the classical forcefield energies, shift by those, and re-centre on zero.
    //d = DB_try_load_forcefield_enes((char  *)DB_subdir_path, &deltaEnes, n_db_strucs);
    deltaEnes = (KMMDEneTrn_t  *)malloc(sizeof(KMMDEneTrn_t)*n_db_strucs);
    for( i = 0; i < n_db_strucs; i++ ){
         deltaEnes[i] = enes[i] - classEnes[i];
    }
    free(enes); //only need the deltas, these are precalculated for the given forcefield.
    free(classEnes); 


    { 
       //repack data, discarding outlier structures.
       int n_new;
       double minEne = 9e99, meanEne = 0.0;
     
 
       n_new = DB_clean_highEneSnaps( &dh_coords,
                                      &preWeights,
                                      &weightGradientsTI, 
                                      &deltaEnes, 
                                       NULL,
                                       n_db_strucs,
                                       use_dh->count );
       fprintf(stderr, "discarding outliers, reduced number of training points from %i to %i\n", n_db_strucs, n_new);

       n_db_strucs = n_new;
       for( i = 0; i < n_db_strucs; i++ ){
          meanEne += deltaEnes[i];
          if ( deltaEnes[i] < minEne ) minEne = deltaEnes[i]; 
       }
       meanEne /= n_db_strucs;
       //for( i = 0; i < n_db_strucs; i++ ){
       //   deltaEnes[i] -= meanEne;
       //}
       printf("mean energy correction: %e   minimum: %e\n \n", meanEne, minEne);
    }

    output_denes = fopen("trimmed_ene_deltas.dat", "w");
    for(i = 0; i < n_db_strucs; i++){
       fprintf(output_denes, "%i %.12f\n", i, deltaEnes[i]);
    }
    fclose(output_denes);
    

    this->n_trainpts = n_db_strucs;
    this->n_dh       = use_dh->count;
    
    //scale the point'point distance for stability.
    //higher dimensionality means higher expected distance between points
    this->sigma2     = sigma2 * sqrt(use_dh->count * 2);

    if( weightGradientsTI ){
        this->acc_dvdl = 0.0;
        this->acc_dvdl_stepcount = 0;
        this->acc_dvdl_steps     = acc_dvdl_steps;
        if( dvdl_logfile_name )
           this->dvdl_logfile = fopen(dvdl_logfile_name, "w");
        else
           this->dvdl_logfile = fopen("KMMD_ti_gradients.txt", "w");
    }



#ifdef KMMD_CUDA
    //marshal data into GPU-friendly structured types.
    this->dh_atoms    = new GpuBuffer<int>(use_dh->count * 4);
    this->dh_trainpts = new GpuBuffer<KMMDCrd_t>(use_dh->count * this->n_trainpts * 2);
    this->preWeights  = new GpuBuffer<KMMDCrd_t>(this->n_trainpts);
    this->trainenes   = new GpuBuffer<KMMDEneTrn_t>(this->n_trainpts);
    this->kmmd_frc    = new GpuBuffer<KMMDFrcAcc_t>(N_atoms * 3);
    
    if( weightGradientsTI ){
       this->dweights_dlambda_scale = new GpuBuffer<KMMDCrd_t>(this->n_trainpts);
       for( i = 0; i < this->n_trainpts; i++ )
           this->dweights_dlambda_scale->_pSysData[i] = weightGradientsTI[i];
    }else{
       this->dweights_dlambda_scale = NULL;
    }

    for( i = 0; i < use_dh->count*4; i++ ){
       this->dh_atoms->_pSysData[i] = dh_atoms[i];
    }
    for( i = 0; i < use_dh->count * this->n_trainpts; i++ ){
       this->dh_trainpts->_pSysData[2*i]   = dh_coords[2*i];
       this->dh_trainpts->_pSysData[2*i+1] = dh_coords[2*i+1];
    }
    for( i = 0; i < this->n_trainpts; i++ ){
       this->trainenes->_pSysData[i]       = deltaEnes[i]; 
       this->preWeights->_pSysData[i]      = preWeights[i]; 
    }

#else
    this->dh_atoms    = new int[ use_dh->count * 4 ];
    this->dh_trainpts = new KMMDCrd_t[use_dh->count * this->n_trainpts * 2];
    this->preWeights  = new KMMDCrd_t[this->n_trainpts];
    this->trainenes   = new KMMDEneTrn_t[this->n_trainpts];
    this->kmmd_frc    = new KMMDFrcAcc_t[N_atoms * 3];

    if( weightGradientsTI ){
       this->dweights_dlambda_scale = new KMMDCrd_t[this->n_trainpts];
       for( i = 0; i < this->n_trainpts; i++ )
           this->dweights_dlambda_scale[i] = weightGradientsTI[i];
    }else{
       this->dweights_dlambda_scale = NULL;
    }
    for( i = 0; i < use_dh->count*4; i++ ){
       this->dh_atoms[i] = dh_atoms[i];
    }
    for( i = 0; i < use_dh->count * this->n_trainpts; i++ ){
       this->dh_trainpts[2*i]   = dh_coords[2*i];
       this->dh_trainpts[2*i+1] = dh_coords[2*i+1];
    }
    for( i = 0; i < this->n_trainpts; i++ ){
       this->trainenes[i]       = deltaEnes[i]; 
       this->preWeights[i]      = preWeights[i]; 
    }
#endif
    free(dh_atoms);
    free(dh_coords);
    free(deltaEnes);
    free(preWeights);
    if( weightGradientsTI ) free( weightGradientsTI );


    //analyse some properties for debug purposes:
    KMMDCrd_t    *test_frc, *dof_deltas, max_frc, this_frc;
    KMMDEneAcc_t  sum_weight;
    KMMDEneAcc_t *trainVec_weights, *frcs_dof, dVdL;
    FILE         *smoothEneFile;
    
    test_frc         = (KMMDCrd_t *)malloc(N_atoms * 3 * sizeof(KMMDCrd_t));
    frcs_dof         = (KMMDCrd_t *)malloc(this->n_dh  * 2 * sizeof(KMMDCrd_t));
    dof_deltas       = (KMMDCrd_t *)malloc(this->n_trainpts * 2 * this->n_dh * sizeof(KMMDCrd_t));
    trainVec_weights = (KMMDCrd_t *)malloc(this->n_trainpts * sizeof(KMMDEneAcc_t));
    smoothEneFile    = fopen("trainpt_smooth_enes.dat", "w");
    
    max_frc  = 0.0;
    for( i = 0; i < this->n_trainpts; i++ ){
     
         //estimate energy and force at each of the training points.
         //force component due to the point itself should be zero, so 
         //this is a good sampling of the range of forces?
         test_dof  = &(this->dh_trainpts[i*this->n_dh*2]);
        
         
         ////////////////////////////////////////////////////////////////////////////
         //cannot test the whole of point_ene_host() because it is too much effort to recover Cartesian coords
         //but can check the important three next-level-down function calls:
         sum_weight = point_dx_vecs_host(this->dh_trainpts, this->preWeights,
                                         test_dof, 
                                         dof_deltas, 
                                         trainVec_weights, 
                                         this->sigma2, 
                                         this->n_dh, 
                                         this->n_trainpts);
                                                  
         point_ene  = reduce_dof_forces_host(dof_deltas, \
                                         trainVec_weights, this->dweights_dlambda_scale, &dVdL,\
                                         test_dof, \
                                         frcs_dof, \
                                         this->sigma2, \
                                         this->trainenes, sum_weight, this->n_dh, this->n_trainpts );
                                         
//         map_dh_forces_host( this->dh_atoms, test_dof, frcs_dof, crd, frc, n_dh );
         ////////////////////////////////////////////////////////////////////////////
         
         
         this_frc = 0.0;
         for(ii = 0; ii < this->n_dh*2; ii++){
            this_frc += frcs_dof[ii]*frcs_dof[ii];
         }
         if( this_frc > max_frc ){
            max_frc = this_frc;
         }
         fprintf(smoothEneFile, "%i true_deltaene: %e  smooth_deltaene: %e |frc|: %.12e\n", 
                  i,  this->trainenes[i], point_ene, sqrt(this_frc));
         
    }
    free(test_frc);    
    free(frcs_dof);    
    free(dof_deltas);  
    free(trainVec_weights);    
    fclose(smoothEneFile);

    printf("sigma2: %.4f scaled to: %.4e gives max force component: %.12e\n",
                    sigma2, this->sigma2, max_frc );
    
    
    
}




kmmdHostContext::kmmdHostContext(char *kmmd_inpFileName) {
    this->Init(kmmd_inpFileName);
}

kmmdHostContext::~kmmdHostContext(){
   
}

#endif
