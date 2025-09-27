//test wrapper for port of cythonised functions from the forcecalc to CUDA.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#include <dirent.h> //Linux/posix: directory listing

#ifdef CUDA
   #include <cuda_runtime.h>
   #include "gpuContext.h"
   #include "gpu.h"
   #ifdef GTI
      #include "gti_cuda.cuh"
   #endif
   #include "kmmd_types.h"

   namespace gppImpl {
     static gpuContext gpu = NULL;
   }
   using namespace std;
   using namespace gppImpl;
#endif

#include "DB_vs_forcefield.h"

#define DB_EXCLUDE_OUTLIERS_NSIGMAS 2

int cmpfunc_qsrtKMMDEneTrn(const void * a, const void * b) {
   if      ( *(KMMDEneTrn_t*)a  < *(KMMDEneTrn_t*)b ) return -1;
   else if ( *(KMMDEneTrn_t*)a == *(KMMDEneTrn_t*)b ) return  0;
   return  1;
}

int  DB_clean_highEneSnaps(KMMDCrd_t    **dh_coords, 
                           KMMDCrd_t    **preWeights, 
                           KMMDCrd_t    **weightGradientsTI,
                           KMMDEneTrn_t **denes, 
                           char        ***DB_fnames,
                           int            n_DB_snaps,
                           int            n_dh_per_snap){

//Go over the energy corrections, if anything is wildly deviated from the mean
//then it is unphysical and should be stripped from the dataset.

     int           i, i_new, n_snaps_new;
     double        mean, stdd, dev;
     KMMDEneTrn_t *denes_new;
     KMMDCrd_t    *coords_new, *preWeights_new, *gradients_new;
     char        **fnames_new;
     double        Q1,Q3;
    


     //find mean and stddev of energies based only on the cetral two quartiles, 
     //so fitting a Gaussian subject to exclusion of outliers
     Q1  = -9e99;
     Q3  = -9e99;
     
     //buffer for sorted energies.  Need sorted version but do not want to change original (yet).
     denes_new  = (KMMDEneTrn_t *)malloc(n_DB_snaps * sizeof(KMMDEneTrn_t));
     memcpy( denes_new, *denes, n_DB_snaps * sizeof(KMMDEneTrn_t) );
     qsort( denes_new, n_DB_snaps, sizeof(KMMDEneTrn_t), cmpfunc_qsrtKMMDEneTrn);
     Q1   = denes_new[(int)(n_DB_snaps/4)];
     Q3   = denes_new[(int)((3*n_DB_snaps)/4)];
     stdd = (Q3 - Q1)/1.333333;                //estimate stddev via interquartile range.
     mean = denes_new[(int)(n_DB_snaps/2)];    //estimate mean as median.
     
     printf("estimated mean, stddev of delta_e (via interquartile range): %e %e\n", mean, stdd);
     printf("Q0,Q1,Q2,Q3: %e  %e  %e  %e\n", (double)denes_new[0], (double)Q1, (double)mean, (double)Q3);    
     printf("first three enes, unsorted: %.12e %.12e %.12e\n", (double)(*denes)[0], (double)(*denes)[1], (double)(*denes)[2]);
  

     //find new size of dataset: how many outliers?
     n_snaps_new = 0;
     for( i = 0; i < n_DB_snaps; i++ ){
         dev = (*denes)[i] - mean;
         if ( dev >  DB_EXCLUDE_OUTLIERS_NSIGMAS * stdd ) continue; 
         if ( dev < -DB_EXCLUDE_OUTLIERS_NSIGMAS * stdd ) continue; 
         n_snaps_new++;
     }
     
     //allocate new smaller arrays and pack down.
     printf("discarding outliers, to %i points from %i\n", n_snaps_new, n_DB_snaps);
     free(denes_new); //will reallocate with smaller size to contain cleaned data.
     denes_new      = (KMMDEneTrn_t *)malloc(n_snaps_new * sizeof(KMMDEneTrn_t));
     coords_new     = (KMMDCrd_t    *)malloc(n_snaps_new * n_dh_per_snap * sizeof(KMMDCrd_t) * 2);
     preWeights_new = (KMMDCrd_t    *)malloc(n_snaps_new * sizeof(KMMDCrd_t));
     if( *weightGradientsTI ){
         gradients_new = (KMMDCrd_t    *)malloc(n_snaps_new * sizeof(KMMDCrd_t));
     }else{
         gradients_new = NULL;
     }
     if( DB_fnames ){
         fnames_new = (char     **)malloc(n_snaps_new * sizeof(char *));
     }
     i_new = 0;
     for( i = 0; i < n_DB_snaps; i++ ){
         dev = (*denes)[i] - mean;
         if ( dev >  DB_EXCLUDE_OUTLIERS_NSIGMAS * stdd ||
              dev < -DB_EXCLUDE_OUTLIERS_NSIGMAS * stdd ) {
              if( DB_fnames ) free( (*DB_fnames)[i] );    
              continue; 
         }
         denes_new[i_new]         = (*denes)[i];
         preWeights_new[i_new]    = (*preWeights)[i];
         if( *weightGradientsTI ) 
             gradients_new[i_new] = (*weightGradientsTI)[i];


         memcpy( &(coords_new[i_new*n_dh_per_snap*2]),
                 &((*dh_coords)[i*n_dh_per_snap*2]),
                  n_dh_per_snap*2*sizeof(KMMDCrd_t) );
                  
         if( DB_fnames ){
            fnames_new[i_new] = (char *)malloc((strlen((*DB_fnames)[i])+1)*sizeof(char));
            memcpy( fnames_new[i_new], (*DB_fnames)[i], (strlen((*DB_fnames)[i])+1)*sizeof(char) );
            free( (*DB_fnames)[i] );
         }
         i_new++;
     }
     free( *dh_coords );
     free( *denes );
     free( *preWeights );
         
     //location in mem has changed now.
    *dh_coords  = coords_new;
    *denes      = denes_new;
    *preWeights = preWeights_new;
     if( *weightGradientsTI ){
         free(*weightGradientsTI);
        *weightGradientsTI = gradients_new;
     }

     FILE *f = fopen("tidied_DB_check.txt", "w");
     if( *weightGradientsTI ){
         for( i = 0; i < n_snaps_new; i++ )
             fprintf(f, "%e    %f %f\n", (*denes)[i], (*preWeights)[i], (*weightGradientsTI)[i] );
     }else{
         for( i = 0; i < n_snaps_new; i++ )
             fprintf(f, "%e    %f\n", (*denes)[i], (*preWeights)[i] );
     }
     fclose( f );
 
     return( n_snaps_new );
}
#undef DB_EXCLUDE_OUTLIERS_NSIGMAS











