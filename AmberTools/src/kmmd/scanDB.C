//test wrapper for port of cythonised functions from the forcecalc to CUDA.
#include <stdio.h>
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

   namespace gppImpl {
      static gpuContext gpu = NULL;
   }
   using namespace std;
   using namespace gppImpl;
#endif


extern "C" {
#include "hash.h" //rough and ready hash table, use for building DH list in wrapper code.
}

#include "kmmd_forcecalc.h"
#include "scanDB.h"

int scan_db(char     *DB_fileList_fname, 
                KMMDCrd_t    **dh_coords, 
                KMMDEneTrn_t **dh_enes,  
                KMMDEneTrn_t **class_enes,  
                KMMDCrd_t    **preWeights,
                KMMDCrd_t    **weightGradientsTI, 
                char        ***DB_fnames, 
                int           *dh_atoms,
                int            n_dh){

    int            n_ats, n_read = 0, n_files = 0; 
    int            i_at, fields_read;   
    double        *snap_coords, mean_ene, test_weight, test_gradient, test_qene, test_ene;
    int            chars_read, line_count;
    char          *line = NULL, element_name[12], *test_fname;
    size_t         llen = 0;
    FILE          *f, *g;


    //loading up the training data.
    f = fopen(DB_fileList_fname, "r");
    if( !f ){
        fprintf(stderr, "Couldn't open KMMD database of reference structures, qenes, classenes, weights\n");
        fprintf(stderr, "input path was: _ %s _\n", DB_fileList_fname);
        exit( 1 );
    }
    
    line       = (char *)malloc(1024*sizeof(char));
    test_fname = (char *)malloc(1024*sizeof(char));
    
    //first pass through the file, to check number of structures
    n_files   = 0;
    while( getline(&line, &llen, f) != -1 ) {
         //require xyz format input structures.  Doesn't matter if we slightly over-allocate space.
         if( strstr(line, "xyz") == NULL ) continue;

         fields_read = sscanf(line, "%s %lf %lf %lf %lf",
                       test_fname, &test_qene, &test_ene, &test_weight, &test_gradient);
         n_files   += 1;
    }
    rewind( f );
    

    //allocate all training pt coords in a single block
    *dh_coords  = (KMMDCrd_t    *)malloc( 2 * n_files * n_dh * sizeof(KMMDCrd_t) ); 
    *dh_enes    = (KMMDEneTrn_t *)malloc( n_files * sizeof(KMMDEneTrn_t) ); 
    *preWeights = (KMMDCrd_t    *)malloc( n_files * sizeof(KMMDCrd_t) ); 

    //don't really need to save classical energies, but will return un-processed from this function if found
    *class_enes = (KMMDEneTrn_t *)malloc( n_files * sizeof(KMMDEneTrn_t) ); 
    
    //allow fragmentation of filenames, could be long.
    if(  DB_fnames != NULL ){
        *DB_fnames  = (char **)malloc( n_files * sizeof(char *) ); 
    }



    //load each file, convert into internal coordinates.
    while( getline(&line, &llen, f) != -1 ) {


        //line format:
        // filename.xyz   weight  classical_ene
        // OR:
        // filename.xyz

        if( strstr(line, "xyz") == NULL ) continue;
        fields_read = sscanf(line, "%s %lf %lf %lf %lf", 
                       test_fname, &test_qene, &test_ene, &test_weight, &test_gradient);
        if( fields_read < 4 ) continue;
      (*dh_enes)[n_read]    = (KMMDEneTrn_t) test_qene; 
      (*class_enes)[n_read] = (KMMDEneTrn_t) test_ene; 
      (*preWeights)[n_read] = (KMMDCrd_t)    test_weight;
        if( fields_read == 5 ){
            if( *weightGradientsTI == NULL ) *weightGradientsTI = (KMMDCrd_t *)malloc(n_files * sizeof(KMMDCrd_t));
          (*weightGradientsTI)[n_read] = test_gradient; 
        }

        //record the filename (if asked for)
        if(  DB_fnames != NULL ){
            (*DB_fnames)[n_read] = (char *)malloc( sizeof(char) * (strlen(test_fname)+1) );
            strcpy(((*DB_fnames)[n_read]), test_fname);
        }

        //read the structure
        g = fopen(test_fname, "r");
        if( !g ){
           fprintf(stderr, "Couldn't open KMMD reference file!\n");
           fprintf(stderr, "input path was: _ %s _\n", test_fname);
           exit( 1 );
        }
        
        
        line_count  = 0; 
        snap_coords = NULL;
        while( 1 ){
            chars_read = getline(&line, &llen, g);
            if( chars_read < 0 ){ break; }
            if( line_count == 0 ){
                //workspace for snap coords, to calc dihedrals from.
                n_ats       = (int)atoi(line);
                snap_coords = (double *)malloc( n_ats * 3 * sizeof(double) );
            }
            //else if( line_count == 1 ){
                //read reference energy for this snap
            //    (*dh_enes)[n_read]  = (KMMDEneTrn_t)strtod( &(line[6]), NULL );
            //    if( n_read < 6 ){
            //         printf("file %i %s has ene: %.12f eV\n", n_read, test_fname, (double)(*dh_enes)[n_read]);
            //    }
            //}
            else if(line_count > 1 ){ 
                i_at = line_count - 2;
                if( sscanf(line, "%s %lf %lf %lf", element_name, 
                             &snap_coords[i_at*3], &snap_coords[i_at*3+1], &snap_coords[i_at*3+2] ) < 4){
                    
                    fprintf(stderr, "stopping read at line: %i : %s\n", line_count, line);
                    break;
                }
            }
            line_count += 1;
            if( line_count - 2 >= n_ats ){break;}
        }
        if( g ) fclose(g);
        if( snap_coords == NULL ){
            fprintf(stderr, "couldn't read snap %s\n", test_fname);
            exit(1);
        }

        //measure the degrees of freedom used to represent the point in kernel space.
        measure_dof_host( &((*dh_coords)[2*n_read*n_dh]), dh_atoms, snap_coords, n_dh );
        mean_ene += (*dh_enes)[n_read];
        n_read   += 1;

        if( snap_coords != NULL ) free(snap_coords);

    }
    fclose( f );
    free( test_fname );
    free( line );

    if( n_read < 1 ){
        fprintf(stderr, "ERROR could not parse DB file %s\n", DB_fileList_fname);
        fprintf(stderr, "line format should be:\n");
        fprintf(stderr, "filename.xyz  Reference_Ene_kcal/mol  FF_Ene_kcal/mol weight <optional: weight gradient>\n");
        fprintf(stderr, "weight is typically equal to 1 unless doing TI or unless optimised as a training parameter.\n");
        fprintf(stderr, "weight gradient need not be specified, unless doing TI, in which case -1 or 1.\n");
    } 


    /*
     Assume units of DB are in standard AMBER kcal/mol
    { 
     int i_snap;
     for( i_snap = 0; i_snap < n_read; i_snap++ ){
        (*dh_enes)[i_snap] *= EV_TO_KCAL_MOL;
    }}
    */
    fprintf(stderr, "first three quantum enes, assumed in kcal/mol: %.6e %.6e %.6f\n",\
                 (double)(*dh_enes)[0], (double)(*dh_enes)[1], (double)(*dh_enes)[2]);

    if( n_read > 10 ){
        if( *weightGradientsTI )
             fprintf(stderr, "example line data 10: xx_fname qene:%e classEne:%e weight:%f sign:%f \n", 
                      (*dh_enes)[10], (*class_enes)[10], (*preWeights)[10],  (*weightGradientsTI)[10] );
        else
             fprintf(stderr, "example line data 10: xx_fname qene:%e classEne:%e weight:%f \n", 
                      (*dh_enes)[10], (*class_enes)[10], (*preWeights)[10] );
             
    }

    return( n_files );
}

