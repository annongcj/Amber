#include <stdio.h>
#include <math.h>

#ifdef CUDA
   #include <cuda_runtime.h>
   #include "gputypes.h"
   #include "gpu.h"
   #ifdef GTI
      #include "gti_cuda.cuh"
   #endif
   #include "gpuContext.h"
#endif

#include "kmmd_types.h"

#ifdef CUDA
   namespace gppImpl {
      static gpuContext gpu = NULL;
   }
   using namespace std;
   using namespace gppImpl;
#endif

#include "kmmd_forcecalc.h"


#define INLINE_CROSS(a,b,c) a.x = b.y*c.z-c.y*b.z; a.y = b.z*c.x-c.z*b.x; a.z = b.x*c.y-c.x*b.y;

#define INLINE_DOT(a, b) (a.x*b.x+a.y*b.y+a.z*b.z)

#define INLINE_VECSCALE(a, s) (a.x) *= (s); (a.y) *= (s); (a.z) *= (s);

void saveVec_host(float *data, int n_pts, const char *filename){
    /* debug function to save a list of datapoints */
    FILE *f;
    int   i;

    f = fopen(filename, "w");
    for( i = 0; i < n_pts; i++ ){
       fprintf(f, "%i %.6f\n", i, data[i]);
    }    
    fclose(f);
}



void measure_dh( KMMDCrd_t* dh_values, int* dh_atoms, KMMDCrd_t* crd, int n_dh ){

     /*Measure values of the dihedrals defined by atoms dh_atoms, takes global atom coordinates in *crd.
       To stripe this function as a kernel, need to make sure that coordinates for the dh_atoms 
       are in the corresponding crds.  Annoying pointer-chase.
       
       After integration, this should be taken care of by existing AMBER code.
      */
        
     int         i, at0, at1, at2,  at3;
     KMMDFloat4  u1, u2, u3, u1Xu2, u2Xu3;  
     KMMDCrd_t   modU2, cPhi, sPhi;         //I guess it makes sense, powers of two and all that... 

     //loop over dihedrals, 4 dh_atoms per dihedral
     for(i = 0; i < n_dh*4; i+=4){
     
         at0 = dh_atoms[i];
         at1 = dh_atoms[i+1];
         at2 = dh_atoms[i+2];
         at3 = dh_atoms[i+3];
         
         u1.x = crd[at1]-crd[at0]; u1.y = crd[at1+1]-crd[at0+1]; u1.z = crd[at1+2]-crd[at0+2];
         u2.x = crd[at2]-crd[at1]; u2.y = crd[at2+1]-crd[at1+1]; u2.z = crd[at2+2]-crd[at1+2];
         u3.x = crd[at3]-crd[at2]; u3.y = crd[at3+1]-crd[at2+1]; u3.z = crd[at3+2]-crd[at2+2];
        
         INLINE_CROSS(u1Xu2, u1, u2);
         INLINE_CROSS(u2Xu3, u2, u3);

         modU2  = sqrt( u2.x*u2.x + u2.y*u2.y + u2.z*u2.z );
         cPhi   = modU2 * (u1.x*u2Xu3.x + u1.y*u2Xu3.y + u1.z*u2Xu3.z);
         sPhi   = u1Xu2.x*u2Xu3.x + u1Xu2.y*u2Xu3.y +u1Xu2.z*u2Xu3.z;

         dh_values[i/4] = (KMMDCrd_t)atan2(sPhi, cPhi);
     }
}

void measure_dof_host( KMMDCrd_t* dof_values, int* dh_atoms, KMMDCrd_t* crd, int n_dh ){

     /*Measure values of the dihedrals defined by atoms dh_atoms, takes global atom coordinates in *crd.
       To stripe this function as a kernel, need to make sure that coordinates for the dh_atoms 
       are in the corresponding crds.  Annoying pointer-chase.

       After measuring dihedrals, no atan2(), just save the sin and cos values adjacent.
       
       After integration, (some of?) this should be taken care of by existing AMBER code.
      */
        
     int        i, at0, at1, at2, at3;
    //changed these from PMEFloat3 to KMMDFloat4because Float3 is not defined in the gputypes.h that I have.
     KMMDFloat4 u1, u2, u3, u1Xu2, u2Xu3;  
     KMMDCrd_t  modU2, cPhi, sPhi;         //I guess it makes sense, powers of two and all that... 

     #define DBG_DOF_HOST

     //loop over dihedrals, 4 dh_atoms per dihedral
     for(i = 0; i < n_dh*4; i+=4){
     
         at0 = 3*dh_atoms[i];
         at1 = 3*dh_atoms[i+1];
         at2 = 3*dh_atoms[i+2];
         at3 = 3*dh_atoms[i+3];

         u1.x = crd[at1]-crd[at0]; u1.y = crd[at1+1]-crd[at0+1]; u1.z = crd[at1+2]-crd[at0+2];
         u2.x = crd[at2]-crd[at1]; u2.y = crd[at2+1]-crd[at1+1]; u2.z = crd[at2+2]-crd[at1+2];
         u3.x = crd[at3]-crd[at2]; u3.y = crd[at3+1]-crd[at2+1]; u3.z = crd[at3+2]-crd[at2+2];
        
         INLINE_CROSS(u1Xu2, u1, u2);
         INLINE_CROSS(u2Xu3, u2, u3);

         modU2  = sqrt( u2.x*u2.x + u2.y*u2.y + u2.z*u2.z );
         cPhi   = modU2 * (u1.x*u2Xu3.x + u1.y*u2Xu3.y + u1.z*u2Xu3.z);
         sPhi   = u1Xu2.x*u2Xu3.x + u1Xu2.y*u2Xu3.y +u1Xu2.z*u2Xu3.z;

         modU2  = 1./sqrt(cPhi*cPhi + sPhi*sPhi); //reuse the variable name for a quick normalisation

         dof_values[i/2]     = cPhi * modU2;
         dof_values[1+(i/2)] = sPhi * modU2;
         #ifdef DBG_DOF_HOST
         if( dof_values[i/2]   != dof_values[i/2] ||
             dof_values[1+i/2] != dof_values[1+i/2] ){
              fprintf(stderr, "ERROR measured a Nan angle for dihedral %i atoms %i %i %i %i\n",
                      i/4, at0, at1, at2, at3);
              fprintf(stderr, "atom 0 coords: %f %f %f\n", crd[at0],crd[at0+1],crd[at0+2]);
              fprintf(stderr, "atom 1 coords: %f %f %f\n", crd[at1],crd[at1+1],crd[at1+2]);
              fprintf(stderr, "atom 2 coords: %f %f %f\n", crd[at2],crd[at2+1],crd[at2+2]);
              fprintf(stderr, "atom 3 coords: %f %f %f\n", crd[at3],crd[at3+1],crd[at3+2]);
              exit( 1 );
         }
         #endif
         
         
     }
     #undef DBG_DOF_HOST
}

KMMDEneAcc_t point_dx_vecs_host(KMMDCrd_t  *dof_values, 
                                    KMMDCrd_t       *trainVec_preWeights,
                                    KMMDCrd_t       *test_dof,
                                    KMMDCrd_t       *dof_deltas,
                                    KMMDEneAcc_t    *trainVec_weights,
                                    KMMDCrd_t        sigma2,
                                    int              n_dh,
                                    int              n_trainpts){
                              
     /*Measure dihedrals defined by dh_atoms[0..4*n_dh] in the coordinates crd[0..3*n_atoms].
      
       After measuring dihedrals, no atan2(), just save the sin and cos values adjacent.
       
       Compare each training point dof_values[0..2*n_dh*n_trainpts] and save the deltas to dof_deltas[0..] 
       
       While we are in the loop, calculate the (un-normalised) weight of each training point.
       
       Return the sum of weights for later normalisation.
       
      */
     int          i, j, ii;
     KMMDCrd_t    inv_sigma2;     
     KMMDEneAcc_t sum_weight;
          
     //loop over training points
     sum_weight =  0.;
     inv_sigma2 = -1./sigma2;
     
#define  _DBG_FWEIGHTS
#ifdef   DBG_FWEIGHTS
     FILE *f;
     f = fopen("values_deltas.dat", "w");
#endif
     
     
     for( j = 0; j < n_trainpts; j++ ){
     
         //loop over degrees of freedom, 2 dof per dihedral.
         trainVec_weights[j] = 0.;
         for(i = 0; i < n_dh*2; i++){
             ii = j*n_dh*2 + i;
             dof_deltas[ii]       = dof_values[ii] - test_dof[i];
             trainVec_weights[j] += dof_deltas[ii] * dof_deltas[ii]; //accumulate the r2 values for each point
             
#ifdef DBG_FWEIGHTS
             fprintf(f, "%i %i %.3f %.3f %.6e\n", j, i, dof_values[ii], test_dof[i], trainVec_weights[j]);
         }
         fprintf(f, "\n");
#else
     }
#endif
             
         trainVec_weights[j] = exp(inv_sigma2 * trainVec_weights[j]) * trainVec_preWeights[j];
         sum_weight         += trainVec_weights[j];
     }                   
#ifdef DBG_FWEIGHTS
     fprintf(f, "\n");  
     fclose(f);
     //exit( 1 ) ;  
#endif
     
     return( sum_weight );
}

KMMDEneAcc_t KMMD_point_dx_vecs_host(KMMDCrd_t     *dof_values, 
                                              KMMDCrd_t     *trainVec_preWeights, 
                                              KMMDCrd_t     *test_dof,
                                              KMMDCrd_t     *dof_deltas,
                                              KMMDEneAcc_t  *trainVec_weights,
                                              KMMDCrd_t      sigma2,
                                              int            n_dh,
                                              int            n_trainpts){
                              
     /*Measure dihedrals defined by dh_atoms[0..4*n_dh] in the coordinates crd[0..3*n_atoms].
      
       After measuring dihedrals, no atan2(), just save the sin and cos values adjacent.
       
       Compare each training point dof_values[0..2*n_dh*n_trainpts] and save the deltas to dof_deltas[0..] 
       
       While we are in the loop, calculate the (un-normalised) weight of each training point.
       
       Return the sum of weights for later normalisation.
       
      */
     int          i, j, i_trainpt, i_delta;
     KMMDCrd_t    inv_sigma2;     
     KMMDEneAcc_t sum_weight;
          
     //Halfway through porting this function to cuda: annoying mix of Float2 and Float types.

     //loop over training points
     sum_weight =  0.;
     inv_sigma2 = -1./sigma2;
     for( j = 0; j < n_trainpts; j++ ){
     
         //loop over degrees of freedom, 2 dof per dihedral.
         trainVec_weights[j] = 0.;

         for(i = 0; i < n_dh; i++){
             i_delta    = (j*n_trainpts*n_dh + i)*2;
             i_trainpt  = (j*n_dh            + i)*2;
             dof_deltas[i_delta]    = dof_values[i_trainpt]   - test_dof[i*2];
             dof_deltas[i_delta+1]  = dof_values[i_trainpt+1] - test_dof[i*2+1];

             //accumulate the r2 values for each point
             trainVec_weights[j]   += dof_deltas[i_delta]   * dof_deltas[i_delta]; 
             trainVec_weights[j]   += dof_deltas[i_delta+1] * dof_deltas[i_delta+1];
         }
         trainVec_weights[j] = exp(inv_sigma2 * trainVec_weights[j]) * trainVec_preWeights[j];
         sum_weight         += trainVec_weights[j];
     }                       
     
     return( sum_weight );
}

KMMDEneAcc_t reduce_dof_forces_host(KMMDCrd_t *dof_deltas,
                                     KMMDCrd_t         *trainVec_weights,
                                     KMMDCrd_t         *trainVec_dWeights_dLambda, 
                                     KMMDEneAcc_t      *dEne_dLambda,
                                     KMMDCrd_t         *crds_dof,
                                     KMMDCrd_t         *frcs_dof,
                                     KMMDCrd_t          sigma2,
                                     KMMDCrd_t         *enes, 
                                     KMMDEneAcc_t       sum_weight,
                                     int                n_dh,
                                     int                n_trainpts ){

     /* collect forces in dof-space, also total energy */
     int            i, j, ii;
     KMMDEneAcc_t   total_ene, iZ2, inv_sum_weight;
     KMMDCrd_t      inv_sigma2, scale;
     
     inv_sigma2     = 1./sigma2;
     inv_sum_weight = 1./sum_weight;
     iZ2            = inv_sum_weight * inv_sum_weight;
     total_ene      = 0.;
     
#define  _DBG_EWEIGHTS
#ifdef   DBG_EWEIGHTS
     FILE *f;
     f = fopen("values_enew.dat", "w");
     fprintf(f, "#sumW: %e  %ix2x%i\n", sum_weight, n_trainpts, n_dh);
#endif
     
     //collect total ene but do not normalise yet.
     if( trainVec_dWeights_dLambda == NULL ){
          for( j = 0; j < n_trainpts; j++ ){
              total_ene += enes[j] * trainVec_weights[j];
          }
     } else {
          KMMDEneAcc_t term;
         *dEne_dLambda = 0.;
          for( j = 0; j < n_trainpts; j++ ){
              term          = enes[j] * trainVec_weights[j];
              total_ene    += term;
             *dEne_dLambda += term * trainVec_dWeights_dLambda[j];
          }
         *dEne_dLambda *= inv_sum_weight;
     }
     
     
     //find forces via chain rule
     for( j = 0; j < n_trainpts; j++ ){
          for( i = 0; i < n_dh*2; i++ ){
              ii           =  j*n_dh*2 + i;
              frcs_dof[i] += 2*dof_deltas[ii]*trainVec_weights[j]*inv_sigma2 *\
                                        ( total_ene - enes[j] * sum_weight ) * iZ2;
              
          }
     }
     
     

     //project forces onto hypertorus tangent plane
     for( i = 0; i < n_dh*2; i+=2 ){
         scale          = frcs_dof[i]*crds_dof[i] + frcs_dof[i+1]*crds_dof[i+1];
//         printf("scale %.4f \n", scale);

#ifdef   DBG_EWEIGHTS
          fprintf(f, "%i  %.12e  %.12e  %.12e\n", i, frcs_dof[i], crds_dof[i], frcs_dof[i]-scale*crds_dof[i]);
          fprintf(f, "%i  %.12e  %.12e  %.12e\n", i, frcs_dof[i+1], crds_dof[i+1], frcs_dof[i+1]-scale*crds_dof[i+1]);
#endif
         
         frcs_dof[i]   -= scale * crds_dof[i];
         frcs_dof[i+1] -= scale * crds_dof[i+1];
         
     }


#ifdef   DBG_EWEIGHTS
     fprintf(f, "#total ene: %e\n", total_ene * inv_sum_weight);
     fclose( f );
#endif
     
     //return normalised energy 
     return( total_ene * inv_sum_weight );                          

}
                      
     

void map_dh_forces_host( int*          dh_atoms,
                                  KMMDCrd_t*    dcrd, 
                                  KMMDFrcTrn_t* dfrc, 
                                  KMMDCrd_t*    crd, 
                                  KMMDFrcAcc_t* frc, 
                                  int           n_dh ){
     
     /* update Cartesian force in place given dihedral force. 
     
        dh_atoms[n_dh*4] : indices of atoms defining each dihedral
        dcrd[n_dh*2]     : dihedrals expressed as cos(Phi), sin(Phi) pairs
        dfrc[n_dh*2]     : forces in the dcrd space
        crd[n_atoms*3]   : Cartesian coordinates
        frc[n_atoms*3]   : forces in Cartesian space
        n_dh             : number of dihedrals.
     
     */
     
     int           i, at0, at1, at2, at3;
     KMMDFloat4    v12, v23, v34;
     KMMDFloat4    v12Xv23, v23Xv34;
     KMMDCrd_t     sinPhi, cosPhi;
     KMMDCrd_t     dUdCos, dUdSin, ccPhi, dUdPhi;
     KMMDFrcAcc_t  mod12, inv12, mod23, inv23, mod34, inv34, cosb, cb2;
     KMMDFrcAcc_t  isinb2, cosc, cc2, isinc2, inv123, inv234; 
     KMMDFrcAcc_t  fa, fb1, fb2, fc1, fc2, fd;


#ifdef DBG_UNITTESTS
     FILE *save_inpcrd, *save_dhats, *save_dudphi, *save_frcdels;

     printf("Sizeof KMMDFloat, int: %lu, %lu\n", sizeof(KMMDFloat), sizeof(int));
     printf("Number of inp dihedrals: %i\n", n_dh);

     save_inpcrd  = fopen("save_inpcrd_50ats.dat", "w");
     save_dudphi  = fopen("save_dudphi.dat",  "w");
     save_frcdels = fopen("save_frcdels.dat",  "w");
#endif



     for( i=0; i < n_dh; i++){
         
         at0 = 3*dh_atoms[4*i];
         at1 = 3*dh_atoms[4*i+1];
         at2 = 3*dh_atoms[4*i+2];
         at3 = 3*dh_atoms[4*i+3];
         
         //vectors defining the dihedrals
         v12.x = crd[at1]-crd[at0]; v12.y = crd[at1+1]-crd[at0+1]; v12.z = crd[at1+2]-crd[at0+2];
         v23.x = crd[at2]-crd[at1]; v23.y = crd[at2+1]-crd[at1+1]; v23.z = crd[at2+2]-crd[at1+2];
         v34.x = crd[at3]-crd[at2]; v34.y = crd[at3+1]-crd[at2+1]; v34.z = crd[at3+2]-crd[at2+2];

         INLINE_CROSS(v12Xv23, v12, v23);
         //l_v12Xv23 = sqrt(v12Xv23.x*v12Xv23.x+v12Xv23.y*v12Xv23.y+v12Xv23.z*v12Xv23.z);

         INLINE_CROSS(v23Xv34, v23, v34);
         //l_v23Xv34 = sqrt(v23Xv34.x*v23Xv34.x+v23Xv34.y*v23Xv34.y+v23Xv34.z*v23Xv34.z);
 
         ccPhi  = INLINE_DOT(v12Xv23, v23Xv34);
         ccPhi /= sqrt(INLINE_DOT(v12Xv23, v12Xv23) * INLINE_DOT(v23Xv34, v23Xv34) );
         if      ( ccPhi < -1.0 )  ccPhi = -1.0;
         else if ( ccPhi >  1.0 )  ccPhi =  1.0;

         //INLINE_CROSS(scr, v12Xv23, v23Xv34);
         //if ( INLINE_DOT(scr, v23) > 0.0 ) theta =  acos( ccPhi );
         //else                              theta = -acos( ccPhi );
         
         //########################################
         //  actual torsion, based on chain rule
         //
         //  dUdPhi = dU_dh0 * dh0_dPhi
         //  
         //
         //
         //chain  rule dUdPhi =  -sin phi * dUdx
         //                    +  cos phi * dUdy
         cosPhi = dcrd[2*i];
         sinPhi = dcrd[2*i+1];
         dUdCos = dfrc[2*i];
         dUdSin = dfrc[2*i+1];

         //Check that everything is parallel to the manifold
         //PMEFloat scale;
         //scale = dfrc[2*i]*dcrd[2*i] + dfrc[2*i+1]*dcrd[2*i+1];
         //printf("scale:  %f\n", scale);
         //dUdPhi = dUdSin*cos(phi) - dUdCos*sin(phi);
         dUdPhi = dUdSin*cosPhi - dUdCos*sinPhi;

    //phi    = atan2(sinPhi, cosPhi);
   // cosPhi = cos(phi - 0.5 * M_PI);
   // sinPhi = sin(phi - 0.5 * M_PI);
   // dUdPhi = dUdCos*cosPhi + dUdSin*sinPhi;  //getting fucking desperate here
    
         mod12 = sqrt(INLINE_DOT(v12, v12));
         inv12 = 1./mod12;
         mod23 = sqrt(INLINE_DOT(v23, v23));
         inv23 = 1./mod23;
         mod34 = sqrt(INLINE_DOT(v34, v34));
         inv34 = 1./mod34;

         cosb   = -INLINE_DOT(v12, v23)*inv12*inv23;
         cb2    =  cosb*cosb;
         isinb2 =  0.0;
         if( cb2 < 0.9999 ) isinb2 = dUdPhi/(cb2 - 1.0);

         cosc   = -INLINE_DOT(v23, v34)*inv23*inv34;
         cc2    =  cosc*cosc;
         isinc2 =  0.0;
         if( cc2 < 0.9999 ) isinc2 = dUdPhi/(cc2 - 1.0);

         inv123   = inv12*inv23;
         inv234   = inv23*inv34;
         INLINE_VECSCALE( v12Xv23, inv123 );
         INLINE_VECSCALE( v23Xv34, inv234 );

         fa  = -inv12 * isinb2;
         fb1 =  (mod23 - mod12*cosb) * inv123 * isinb2;
         fb2 =  cosc * inv23 * isinc2;
         fc1 =  (mod23 - mod34*cosc) * inv234 * isinc2;
         fc2 =  cosb * inv23 * isinb2;
         fd  =  inv34 * isinc2;


         //#define SC (1.0/4184.0)
         #define SC 1.0
         fa*=SC; fb1*=SC; fb2*=SC; fc1*=SC; fc2*=SC; fd*=SC;
         #undef  SC

         frc[at0]   +=  v12Xv23.x * fa;
         frc[at0+1] +=  v12Xv23.y * fa;
         frc[at0+2] +=  v12Xv23.z * fa;

         frc[at1]   +=  fb1 * v12Xv23.x - fb2 * v23Xv34.x;
         frc[at1+1] +=  fb1 * v12Xv23.y - fb2 * v23Xv34.y;
         frc[at1+2] +=  fb1 * v12Xv23.z - fb2 * v23Xv34.z;

         frc[at2]   += -fc1 * v23Xv34.x + fc2 * v12Xv23.x;
         frc[at2+1] += -fc1 * v23Xv34.y + fc2 * v12Xv23.y;
         frc[at2+2] += -fc1 * v23Xv34.z + fc2 * v12Xv23.z;

         frc[at3]   +=  v23Xv34.x * fd;
         frc[at3+1] +=  v23Xv34.y * fd;
         frc[at3+2] +=  v23Xv34.z * fd;

#ifdef DBG_UNITTESTS
         fprintf(save_dudphi, "%.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",\
                       dUdPhi, fa, fb1, fb2, fc1, fc2, fd);

         fprintf(save_frcdels, "%.6f %.6f %.6f  %.6f %.6f %.6f  %.6f %.6f %.6f  %.6f %.6f %.6f\n",
                       v12Xv23.x * fa, v12Xv23.y * fa, v12Xv23.z * fa,
                       fb1 * v12Xv23.x - fb2 * v23Xv34.x,
                       fb1 * v12Xv23.y - fb2 * v23Xv34.y,
                       fb1 * v12Xv23.z - fb2 * v23Xv34.z,
                       -fc1 * v23Xv34.x + fc2 * v12Xv23.x,
                       -fc1 * v23Xv34.y + fc2 * v12Xv23.y,
                       -fc1 * v23Xv34.z + fc2 * v12Xv23.z,
                       v23Xv34.x * fd, v23Xv34.y * fd, v23Xv34.z * fd);
#endif

    }
#ifdef DBG_UNITTESTS
    fclose(save_dudphi);
    fclose(save_frcdels);
#endif
}


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
                                      char          *calc_id ){

     /*Measure dihedrals defined by dh_atoms[0..4*n_dh] in the coordinates crd[0..3*n_atoms].
      
       After measuring dihedrals, no atan2(), just save the sin and cos values adjacent.
       
       Compare each training point dof_values[0..2*n_dh*n_trainpts] using the Gaussian kernel and take a weighted mean 
       of the training point energies enes[0..n_trainpts].
       
      */
        
     KMMDFrcAcc_t *frcs_dof;
     KMMDEneAcc_t *trainVec_weights;
     KMMDCrd_t    *test_dof, *dof_deltas;     
     KMMDEneAcc_t  sum_weight, point_ene;
          
     //workspaces.
     //printf("allocating workspaces for %i dihedrals, %i training points\n", n_dh, n_trainpts);
     //calloc: allocate and set to zero.
     test_dof         = (KMMDCrd_t    *)calloc(n_dh*2,     sizeof(KMMDCrd_t));
     dof_deltas       = (KMMDCrd_t    *)calloc(n_trainpts*n_dh*2, sizeof(KMMDCrd_t)); 
     trainVec_weights = (KMMDEneAcc_t *)calloc(n_trainpts, sizeof(KMMDEneAcc_t)); 
     frcs_dof         = (KMMDFrcAcc_t *)calloc(n_dh*2,     sizeof(KMMDFrcAcc_t)); 
    
     if( !test_dof || !dof_deltas || !trainVec_weights || !frcs_dof ){
        fprintf(stderr,
           "some problems allocating memory. Your system is too big or you have too many training pts.\n");
        exit( 1 );
     }

     //measure dihedral values for the test point only
     measure_dof_host( test_dof, dh_atoms, crd, n_dh );
     if ( calc_id != NULL ){  
            char dbg_fname[128];
            sprintf(dbg_fname, "%s_testpoint_dof.dat", calc_id);
        //    saveVec_host(test_dof, n_dh*2, dbg_fname); 
            
            //sprintf(dbg_fname, "%s_testpoint_crd.dat", calc_id);
            //saveVec_host(crd, 124*3, dbg_fname);  //hardode number of atoms, not needed for non-debug code.
     }

     //get un-normalised deltas to each training point, and also sum of weights.
     sum_weight = point_dx_vecs_host(dof_values,\
                                     trainpt_preWeights,\
                                     test_dof,\
                                     dof_deltas,\
                                     trainVec_weights,\
                                     sigma2,\
                                     n_dh,\
                                     n_trainpts);
     
     
     //normalise and reduce to get force in dof-space, plus total energy
     point_ene  = reduce_dof_forces_host(dof_deltas, \
                                         trainVec_weights, trainpt_dWeights_dLambda, dVdL, \
                                         test_dof, \
                                         frcs_dof, \
                                         sigma2, \
                                         enes, sum_weight, n_dh, n_trainpts );
                                         
     //printf("KMMD sum_w, point ene: %e %e\n", sum_weight, point_ene);
     if( sum_weight != sum_weight || point_ene != point_ene){
         fprintf(stderr, "Error! nans.  Stopping in func point_ene_host().\n");
         fprintf(stderr, "point ene:  %e\n", point_ene);
         fprintf(stderr, "sum weight: %e\n", sum_weight);
         exit( 1 );
     }
     
     
     if ( calc_id != NULL ){  
            char dbg_fname[128];
            sprintf(dbg_fname, "%s_testpoint_dof_forces.dat", calc_id);
         //   saveVec_host(frcs_dof, n_dh*2, dbg_fname); 
            
            //sprintf(dbg_fname, "%s_testpoint_crd.dat", calc_id);
            //saveVec_host(crd, 124*3, dbg_fname);  //hardode number of atoms, not needed for non-debug code.
     }
     
     //convert dof-space forces to cartesian forces on atoms
     map_dh_forces_host( dh_atoms, test_dof, frcs_dof, crd, frc, n_dh );
     
     free(test_dof);
     free(dof_deltas);
     free(trainVec_weights);
     free(frcs_dof);
     
     return( point_ene );
}

KMMDEneAcc_t KMMD_point_ene_host( kmmdHostContext *KMMD,  KMMDCrd_t *crd, KMMDFrcAcc_t *frc){

     /*Measure dihedrals defined by dh_atoms[0..4*n_dh] in the coordinates crd[0..3*n_atoms].
      
       After measuring dihedrals, no atan2(), just save the sin and cos values adjacent.
       
       Compare each training point dof_values[0..2*n_dh*n_trainpts] using the Gaussian kernel and take a weighted mean 
       of the training point energies enes[0..n_trainpts].
       
      */
        
     KMMDFrcAcc_t  *frcs_dof; 
     KMMDCrd_t     *test_dof, *dof_deltas;     
     KMMDEneAcc_t   sum_weight, point_ene, *trainVec_weights, dV_dL;
          
     //workspaces.
     //printf("allocating workspaces for %i dihedrals, %i training points\n", n_dh, n_trainpts);
     //calloc: allocate and set to zero.
     test_dof         = (KMMDCrd_t    *)calloc(KMMD->n_dh*2,                  sizeof(KMMDCrd_t));
     dof_deltas       = (KMMDCrd_t    *)calloc(KMMD->n_trainpts*KMMD->n_dh*2, sizeof(KMMDCrd_t)); 
     trainVec_weights = (KMMDEneAcc_t *)calloc(KMMD->n_trainpts, sizeof(KMMDEneAcc_t)); 
     frcs_dof         = (KMMDFrcAcc_t *)calloc(KMMD->n_dh*2,     sizeof(KMMDFrcAcc_t)); 
     if( !test_dof || !dof_deltas || !trainVec_weights || !frcs_dof ){
        fprintf(stderr,
           "some problems allocating memory. Your system is too big or you have too many training pts.\n");
        exit( 1 );
     }

#ifdef KMMD_CUDA
     //measure dihedral values for the test point only
     measure_dof_host( test_dof, KMMD->dh_atoms->_pSysData, crd, KMMD->n_dh );

     //get un-normalised deltas to each training point, and also sum of weights.
     sum_weight = KMMD_point_dx_vecs_host(KMMD->dh_trainpts->_pSysData, KMMD->preWeights->_pSysData,\
                                          test_dof, dof_deltas, trainVec_weights,\
                                          KMMD->sigma2, KMMD->n_dh, KMMD->n_trainpts);
     
     //normalise and reduce to get force in dof-space, plus total energy
     point_ene  = reduce_dof_forces_host(dof_deltas, \
                                         trainVec_weights, KMMD->dweights_dlambda_scale->_pSysData, &dV_dL,\
                                         test_dof, \
                                         frcs_dof, \
                                         KMMD->sigma2, \
                                         KMMD->trainenes->_pSysData, \
                                         sum_weight,\
                                         KMMD->n_dh,\
                                         KMMD->n_trainpts );
     
     //convert dof-space forces to cartesian forces on atoms
     map_dh_forces_host( KMMD->dh_atoms->_pSysData,\
                         test_dof, frcs_dof, crd, frc, KMMD->n_dh );
#else
     //measure dihedral values for the test point only
     measure_dof_host( test_dof, KMMD->dh_atoms, crd, KMMD->n_dh );

     //get un-normalised deltas to each training point, and also sum of weights.
     sum_weight = KMMD_point_dx_vecs_host(KMMD->dh_trainpts, KMMD->preWeights,\
                                          test_dof, dof_deltas, trainVec_weights,\
                                          KMMD->sigma2, KMMD->n_dh, KMMD->n_trainpts);
     
     //normalise and reduce to get force in dof-space, plus total energy
     point_ene  = reduce_dof_forces_host(dof_deltas, \
                                         trainVec_weights,  KMMD->dweights_dlambda_scale, &dV_dL,\
                                         test_dof, \
                                         frcs_dof, \
                                         KMMD->sigma2, \
                                         KMMD->trainenes, \
                                         sum_weight,\
                                         KMMD->n_dh,\
                                         KMMD->n_trainpts );
     
     //convert dof-space forces to cartesian forces on atoms
     map_dh_forces_host( KMMD->dh_atoms,\
                         test_dof, frcs_dof, crd, frc, KMMD->n_dh );
#endif
     
     free(test_dof);
     free(dof_deltas);
     free(trainVec_weights);
     free(frcs_dof);
     
     return( point_ene );
}




