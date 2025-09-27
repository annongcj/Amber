#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef CUDA
   #include "gpuContext.h"
   #include "gpu.h"
   #ifdef GTI
     #include "gti_cuda.cuh"
   #endif
#endif

extern "C" {
#include "hash.h" //rough and ready hash table, use for building DH list in wrapper code.
}

#include "pdb_io.h"

/***************
Start the file with some local-scope utility functions
***************/

char **split4_chars( char *space_list ){
   //utility function to split a space-separated list into an array of 4 strings
   char **x;
   x    = (char **)malloc(4*sizeof(char *));
   x[0] = (char  *)malloc(16*sizeof(char));
   x[1] = (char  *)malloc(16*sizeof(char));
   x[2] = (char  *)malloc(16*sizeof(char));
   x[3] = (char  *)malloc(16*sizeof(char));
   sscanf(space_list, "%s %s %s %s", x[0],x[1],x[2],x[3]);
   return( x );
}

int *split4_ints( char *space_list ){
   //utility function to split a space-separated list into an array of 4 integers
   int *x, got;
   
   x   = (int *)malloc(4*sizeof(int));
   got = sscanf(space_list, "%d %d %d %d", &(x[0]),&(x[1]),&(x[2]),&(x[3]));
   if( got != 4 ){
       fprintf(stderr, "Error! could not scan four ints from json entry: _%s_\n", space_list);
       exit( 8 );
   }
   return( x );
}

void **split2_chars_int( char *space_list ){
   //utility function to split a space-separated list into a string and an int.
   void  **x;
   int     got;
   
   x    = (void **)malloc(sizeof(void*)*2);
   x[0] = (void  *)malloc(sizeof(char)*32);
   x[1] = (void  *)malloc(sizeof(int));
   got  = sscanf(space_list, "%s %d", (char *)(x[0]), (int *)(x[1]));
   if( got != 2 ){
       fprintf(stderr, "Error! could not scan str+int from json entry: _%s_\n", space_list);
       exit( 8 );
   }
   return( x );
}


int  get_pdb(const char *path, KMMDCrd_t **coords, char ***atnames, int **resids){
    /* read of a pdb file, just grabbing info needed for KMMD calculation */
    FILE   *f;
    size_t  len = 0;
    int     chars_read;
    char   *line = NULL;
    char   *pdb_field_charbuf;
    int     n_ats = 0, i_at = 0, from_p, to_p, d;
    
    f = fopen(path, "r");
    
    if(!f){fprintf(stderr, "no file %s\n", path); exit(-1);}
    else{fprintf(stderr, "reading pdb: %s\n", path);}
    while( 1 ){
        chars_read = getline(&line, &len, f);
        if( chars_read < 0 ){ break; }
        if( line[0] == 'A' && line[1] == 'T' && line[2] == 'O' && line[3] == 'M' ){
           n_ats += 1;
        }
    }
    fclose(f);
    fprintf(stderr, "   : %i atoms\n", n_ats);
    
    *coords  = (KMMDCrd_t *)    malloc( n_ats * sizeof(KMMDCrd_t) * 3);
    *atnames = (char **)        malloc( n_ats * sizeof(char*));
    *resids  = (int *)          malloc( n_ats * sizeof(int));
    pdb_field_charbuf = (char *)malloc( 12 * sizeof(char) );
    
    f = fopen(path, "r");
    if(!f){fprintf(stderr, "no file %s\n", path); exit(-1);}
    while ( 1 ){
        chars_read = getline(&line, &len, f);
        if( chars_read < 0 ){ break; }
        if( line[0] != 'A' || line[1] != 'T' || line[2] != 'O' || line[3] != 'M' )
           continue;
           
        /*get atnames*/
        (*atnames)[i_at] = (char *)calloc( 8, sizeof(char) ); //calloc inits to zero == NULL = string terminator
        from_p = 12;
        to_p   = 0;
        while(from_p < chars_read){
           if( from_p > 15 ){ break; }
           if( line[from_p] == ' ' ){ 
               from_p++;
           }else{
               (*atnames)[i_at][to_p++] = line[from_p++];
           }
        }
        
        /* get atom coords */
        for( d = 0; d < 3; d++ ){
            memset(pdb_field_charbuf, 0, 12*sizeof(char) );
            memcpy(pdb_field_charbuf, &line[30+8*d], 8*sizeof(char));
            (*coords)[3*i_at+d] = atof(pdb_field_charbuf);
        }
        
        /* residue IDs */
        memset(pdb_field_charbuf, 0, 12*sizeof(char) );
        memcpy(pdb_field_charbuf, &line[22], 4*sizeof(char));
        (*resids)[i_at] = (int)atoi(pdb_field_charbuf);
        
        i_at += 1;
    }
    fclose(f);
    free(pdb_field_charbuf);

    return( n_ats );
}

void   write_pdb(char *path, FILE *g, KMMDCrd_t *coords, int n_atoms){

    FILE   *f;
    size_t  len = 0;
    int     chars_read, atoms_read;
    char   *line = NULL;
    char    pdb_field_charbuf[128];
    
    f = fopen(path, "r");
    if( !f ){
        fprintf(stderr, "problem opening %s\n", path); exit(-1);
    }

    atoms_read = 0;
    while ( atoms_read < n_atoms ){
        chars_read = getline(&line, &len, f);
        if( chars_read < 0 ){ break; }
        if( line[0] != 'A' || line[1] != 'T' || line[2] != 'O' || line[3] != 'M' )
           continue;

        //cut the template line off at the start of the coordinates field.
        line[30] = '\0';
        sprintf(pdb_field_charbuf, "%s%8.3f%8.3f%8.3f\n", line, 
                          coords[atoms_read*3],coords[atoms_read*3+1],coords[atoms_read*3+2]);
        fprintf(g, "%s", pdb_field_charbuf);
        atoms_read += 1;
   }
   fclose(f);
}



int*  assign_dh_atoms(HashTable  *dh_names,\
                              HashTable  *dh_resdels,\
                              HashTable  *dh_atnames,\
                                   char **atnames,\
                                   int   *resids,\
                                   int    n_ats){
 /*
     Wrapper code to build list of dihedrals to treat, for standalone test calculation on GG.CC duplex.
 */                                  
    int         i_dh, ii_dh, n_dh, i_at, *dh_atoms;
    LinkedList *dhItemList;
    void      **chars_int;
    char       *dhname, **dhSet_atnames, *dummy;
    int         dhresid, *dhSet_deltas, found_OK;
    
    fprintf(stderr, "getting n_dh...\n");
    
    n_dh = dh_names->count; //how many dihedrals define our system
    fprintf(stderr, "    ... %i\n", n_dh);
    
    dh_atoms   = (int *)malloc(n_dh * sizeof(int) * 4);
    dhItemList = dh_names->keyList;
    i_dh       = n_dh - 1; //fill the array backwards, we are reading from a LIFO.
    while( dhItemList ){
    
       fprintf(stderr, "processing dihedral with key: _%s_\n", dhItemList->item->key);
    
       chars_int  = split2_chars_int( dhItemList->item->key );
       dhname     =  (char *)chars_int[0];
       dhresid    = *((int*)chars_int[1]);

       fprintf(stderr, "name: _%s_  resid: _%i_\n", dhname, dhresid);


       //get chain offsets for member residues of this dihedral type
       char *got_value;
       got_value     = ht_search(dh_resdels,  dhname);
       if( got_value == NULL ){
           fprintf(stderr, "Error, could not load matching resdel key to _%s_\n", dhname);
           fprintf(stderr, "Check the KMMD input json file.\n");
           exit( 8 );
       }
       dhSet_deltas  = split4_ints( got_value );
       
       got_value     = ht_search(dh_atnames, dhname);
       if( got_value == NULL ){
           fprintf(stderr, "Error, could not load matching dh_atnames key to _%s_\n", dhname);
           fprintf(stderr, "Check the KMMD input json file.\n");
           exit( 8 );
       }
       dhSet_atnames = split4_chars( got_value );
       
       for(ii_dh = 0; ii_dh < 4; ii_dh+=1 ){

          //get atom id based on residue atnames and resindex
          found_OK = 0;
          for( i_at = 0; i_at < n_ats; i_at++ ){
              if( resids[i_at] != dhresid + dhSet_deltas[ii_dh] ) continue; 
              if(  strcmp( atnames[i_at], dhSet_atnames[ii_dh]) != 0 ) {
                  //fprintf(stderr, "  res: %i at: %i %s\n", resids[i_at], i_at, atnames[i_at]);
                  continue;
              }
              dh_atoms[i_dh*4 + ii_dh] = i_at;  //found the atom: resid and atom name match.
              found_OK = 1;
              break; 
          }
          if( found_OK != 1 ){
              fprintf(stderr, "Error, could not find atom %s in resid %i dhname %s\n",
                                dhSet_atnames[ii_dh], dhresid + dhSet_deltas[ii_dh], dhname   );

              for( i_at = 0; i_at < n_ats; i_at++ ){
                 fprintf(stderr, "  compared: %i:%s to %i:%s and got %i:%i \n", 
                        resids[i_at], atnames[i_at], dhresid + dhSet_deltas[ii_dh], dhSet_atnames[ii_dh], 
                        strcmp( atnames[i_at], dhSet_atnames[ii_dh]), strcmp( atnames[i_at], dhSet_atnames[ii_dh]));
              }

              exit( 1 );
          }
       } 
       free(dhSet_deltas);       
       for( ii_dh = 0; ii_dh < 4; ii_dh++ ) free(dhSet_atnames[ii_dh]);
       free(dhSet_atnames);

       printf("  %i %i:%s ats: %i %i %i %i\n",
                  i_dh, dhresid, dhname, dh_atoms[i_dh*4], dh_atoms[i_dh*4+1], dh_atoms[i_dh*4+2], dh_atoms[i_dh*4+3] );

       dhItemList = dhItemList->next;
       i_dh      -= 1;
    }
    return( dh_atoms ); 
}




