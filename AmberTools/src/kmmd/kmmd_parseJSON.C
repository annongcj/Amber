#ifndef _KMMD_PARSE_JSON
#define _KMMD_PARSE_JSON

#include <string.h>

extern "C" {
#include "hash.h"
}

#include "kmmd_types.h"
#include "kmmd_context.h"

#include "scanDB.h"
#include "DB_vs_forcefield.h"
#include "kmmd_forcecalc.h"

#include "pdb_io.h"

int getStrTag(char* line_buf, const char* tag, size_t read_len, char **got_str){
    char  *pos;
    int    i;
    int    read_start, read_stop;
    int    seen_quotes, write_posn;
    int    actual_line_len;
    
    //check line for occurence of tag in the form: "tag_name": "string"
    pos = strstr(line_buf, tag);
    
    if( pos == NULL ) return -1; //-1: tag not found
    
    actual_line_len = (int)strlen(line_buf);
   // if( pos )printf("  checked line %s length %i strlen %i for tag %s got: %p\n", 
   //          line_buf, read_len, actual_line_len, tag, pos);
    read_start = -1;
    read_stop  = -1;
    for( i = int((pos - line_buf)/sizeof(char)); i < (int)read_len && i < actual_line_len;i++ ){
        if( line_buf[i] == ':' ){
           read_start = i+1;
           continue;
        }
        if( read_start > -1 &&\
                ( line_buf[i]=='}' || line_buf[i]==',' ) ){
           read_stop = i-1;
           break;
        }
    }
    if( read_start != -1 && read_stop == -1 ){
        printf("found tag _%s_ but no value. This can be OK.\n", tag);
        return 0; //tag found but no data 
    }
    
    //copy data if found.
   *got_str = (char *)malloc( sizeof(char)*(2 + read_stop - read_start));
    memset(*got_str, '\0', sizeof(char)*(2 + read_stop - read_start));
   
    seen_quotes = 0; 
    write_posn   = 0;
    for( i = read_start; i<= read_stop; i++ ){
       if( line_buf[i]  == '"' ){
          seen_quotes += 1; 
       }else if( seen_quotes == 1 ){
          (*got_str)[write_posn++] = line_buf[i];
       }
    }
  (*got_str)[write_posn] = '\0';
  
    return write_posn; //tag found with data.
}
int getIntTag(char* line_buf, const char* tag, size_t read_len, int *got_int){

    char   *pos;
    int     i;
    int     read_start, read_stop;
    int     got_value_count;
    int     actual_line_len;

    //check line for occurence of tag in the form: "tag_name": int
    pos = strstr(line_buf, tag);
    if( pos == NULL ) return -1; //-1: tag not found
    read_start = -1;
    read_stop  = -1;
    actual_line_len = (int)strlen(line_buf);
    for( i = int((pos - line_buf)/sizeof(char)); i < (int)read_len && i < actual_line_len;i++ ){
        if( line_buf[i] == ':' ){
            read_start = i+1;
            continue;
        }
        if( read_start > -1 &&\
            ( line_buf[i]=='}' || line_buf[i]==',' ) ){
            read_stop = i-1;
            break;
        }
    }
    if( read_start != -1 && read_stop == -1 ){
        printf("found tag _%s_ but no value. This can be OK.\n", tag);
        return 0; //tag found but no data 
    }

    got_value_count = sscanf( (char *)(&line_buf[read_start]), "%i", got_int );

    return got_value_count;
}


int getFloatTag(char* line_buf, const char* tag, size_t read_len, double *got_float){
    char   *pos;
    int     i;
    int     read_start, read_stop;
    int     got_value_count;
    int     actual_line_len;
    
    //check line for occurence of tag in the form: "tag_name": "string"
    pos = strstr(line_buf, tag);
    if( pos == NULL ) return -1; //-1: tag not found
    read_start = -1;
    read_stop  = -1;
    actual_line_len = (int)strlen(line_buf);
    for( i = int((pos - line_buf)/sizeof(char)); (int)read_len && i < actual_line_len;i++ ){
        if( line_buf[i] == ':' ){
           read_start = i+1;
           continue;
        }
        if( read_start > -1 &&\
                ( line_buf[i]=='}' || line_buf[i]==',' ) ){
           read_stop = i-1;
           break;
        }
    }
    if( read_start != -1 && read_stop == -1 ){
        printf("found tag _%s_ but no value. This can be OK.\n", tag);
        return 0; //tag found but no data 
    }
    
    got_value_count = sscanf( (char *)(&line_buf[read_start]), "%lf", got_float ); 
    
    return got_value_count;
}

int getHashEntry( char* line_buf, char** got_key, char** got_value){
    
    char *pos;
    int   i;
    int   read_start, read_stop;
    int   status;
    
    //check line for occurence of tag in the form: {"something": "something"}
    pos = strstr(line_buf, "\"");
    if( pos == NULL ) return 0; //no entry on this line.
    
    i = (int)((pos - line_buf)/sizeof(char));
    read_start = -1;
    read_stop  = -1;
    while( line_buf[i] ){
       if( read_start == -1 && line_buf[i] == '"' ){
          read_start = i+1;
          i += 1;
          continue;
       }
       if( read_start != -1 && line_buf[i] == '"'){
          read_stop = i-1;
          break;
       }
       i += 1;
    }
    if (read_stop < read_start ) return 0; //empty string.
    if (read_start == -1 || read_stop == -1) return 0;//no string or unterminated
    
    //copy data if found.
   *got_key = (char *)malloc( sizeof(char)*(2 + read_stop - read_start) );
    for( i = read_start; i<= read_stop; i++ ){
       (*got_key)[i-read_start] = line_buf[i];
    }
  (*got_key)[read_stop-read_start+1] = '\0';
  
  
    //now check for value, can return OK with only a key, and set empty value.
    pos = strstr(line_buf, ":");
    if( pos == NULL ){
       *got_value = (char *)malloc( sizeof(char) );
      **got_value = '\0';
       return 1; //have key but not value
    }
    i = (int)(pos - line_buf);
    read_start = -1;
    read_stop  = -1;
    while( line_buf[i] != '\0' ){
       if( read_start == -1 && line_buf[i] == '"' ){
          read_start = i+1;
          i += 1;
          continue;
       }
       if( read_start != -1 && line_buf[i] == '"'){
          read_stop = i-1;
          break;
       }
       i += 1;
    }
    if( (read_stop < read_start )
     || (read_start == -1 || read_stop == -1) ){
      *got_value = (char *)malloc( sizeof(char) );
     **got_value = '\0';
       return 1; 
    }
    
    //copy data if found.
   *got_value = (char *)malloc( sizeof(char)*(2 + read_stop - read_start) );
    for( i = read_start; i<= read_stop; i++ ){
       (*got_value)[i-read_start] = line_buf[i];
    }
  (*got_value)[read_stop-read_start+1] = '\0';
    
    return 2;
}

int check_for_closeTags(char *line){

   int i, count_close;
   i           = 0;
   count_close = 0;
   while( line[i] != '\0' ){
      if( line[i] == ']' ) count_close+=1; 
      i++;
   }
   return count_close;
}


int scan_json_tags( FILE *f, char **DB_fileList, char **ref_pdb,
                    char **dvdl_logfile_name,    int   *acc_dvdl_steps,
                    HashTable **dh_atnames, HashTable **dh_resdel, HashTable **dh_byres, 
                    KMMDCrd_t  *sigma2){

    char *line_buf;
    size_t  read_len = 1024;
    char *pos;
    char *got_key   = NULL;
    char *got_value = NULL;
    
    HashTable *active_hash = NULL;
    
    line_buf = (char *)malloc(read_len*sizeof(char));
    while (getline(&line_buf, &read_len, f) != -1) {
        
        if( *DB_fileList == NULL ){
	        if( getStrTag(line_buf, "DB_fileList", read_len, DB_fileList) > 0 ){
	           printf("got list of struct files, optionally weights: '%s\n'", *DB_fileList);
	           continue;
	        }
        }
        if( *ref_pdb == NULL){
	        if( getStrTag(line_buf, "ref_pdb", read_len, ref_pdb) > 0 ){
	           printf("got ref pdb: '%s'\n", *ref_pdb);
	           continue;
	        }
        }
        if( *dvdl_logfile_name == NULL ){
                if( getStrTag(line_buf, "dvdl_logfile", read_len, dvdl_logfile_name) > 0 ){
                   printf("defined dvdl_logfile: '%s'\n", *dvdl_logfile_name);
                   continue;
                }
        } 

        if( getIntTag( line_buf, "acc_dvdl_steps", read_len, acc_dvdl_steps) > 0 ){
               printf("got acc_dvdl_steps : %i\n", *acc_dvdl_steps);
        }        

	double hold_double;
	if( getFloatTag(line_buf, "sigma2", read_len, &hold_double) > 0 ){
	      *sigma2 = (KMMDCrd_t)hold_double;
	       printf("got sigma2: '%f'\n", (float)(*sigma2));
	       continue;
	}

	    
	//if we find a tag indicating start of a hash table, then create 
	//one and add subsequent hash entries to it.
	if( getStrTag(line_buf, "dh_atnames", read_len, &pos) == 0){
	       printf("creating hash: dh_atnames\n");
	      *dh_atnames = create_table(CAPACITY);
	       active_hash = *dh_atnames;
	       continue;
	}
	if( getStrTag(line_buf, "dh_resdel", read_len, &pos) == 0){
	       printf("creating hash: dh_resdel\n");
	      *dh_resdel = create_table(CAPACITY);
	       active_hash = *dh_resdel;
	       continue;
	}
	if( getStrTag(line_buf, "dh_byres", read_len, &pos) == 0){
	       printf("creating hash: dh_byres\n");
	      *dh_byres = create_table(CAPACITY);
	       active_hash = *dh_byres;
	       continue;
	}
	    
	    //whatever key-value pair we find, add it to the active hash table.
	if( active_hash != NULL ){
	    getHashEntry( line_buf, &got_key, &got_value);
    	    if( got_key != NULL ){
    	       ht_insert(active_hash, got_key, got_value);
    	       free(got_key);
    	       free(got_value);
    	       got_key   = NULL;
    	       got_value = NULL;
    	    }
    	    if( check_for_closeTags(line_buf) > 0 ){
    	       active_hash = NULL;
    	    }
    	}
    }
    if( line_buf != NULL ) free(line_buf);
    
    if( *DB_fileList == NULL ){
        fprintf(stderr, "MISSING KEY FROM KMMD JSON FILE: 'DB_fileList'\n");
        exit(8);
    }
    if( *ref_pdb == NULL ){
        fprintf(stderr, "MISSING KEY FROM KMMD JSON FILE: 'ref_pdb'\n");
        exit(8);
    }
    if( *dh_atnames == NULL || *dh_resdel == NULL || *dh_byres == NULL ){
        fprintf(stderr, "MISSING ONE OR MORE JSON ARRAYS FROM KMMD CONTROL FILE\n");
        
        if( *dh_atnames == NULL ) { fprintf(stderr, "  ..'dh_atnames' is NULL\n");}
        if( *dh_resdel  == NULL ) { fprintf(stderr, "  ..'dh_resdel'  is NULL\n");}
        if( *dh_byres   == NULL ) { fprintf(stderr, "  ..'dh_byres'   is NULL\n");}
        
        exit(8);
    }
    
    return 0;
}

#ifdef DBG_JSON
int main(){
    FILE *f;
         
    //load a pdb structure, identify dihedrals
    char *DB_fileList;  // = "../../DYADS_RESULTS_DB/DYADS_DB_05_2021/GG";
    char *ref_pdb;//        = "../../DYADS_RESULTS_DB/DYADS_DB_05_2021/GG/GG_nW20.pdb";
     
    /* data structures to record names of dihedrals and ids of defining atoms */
    HashTable  *dh_atnames, *resdel, *use_dh;
    
    KMMDCrd_t   sig2;
    
    
    //f = fopen("/home/josh/stretch_2/AIMD_clean_checkout/sluk/TEANMD/TEAN_AIMD_CUDA/test_amber/kmmd_test.json", "r");
    
    f = fopen("/home/josh/stretch_2/AIMD_clean_checkout/sluk/AIMD/STRETCH_GG/GGCC_setup.json", "r");
    
    if( !f ){
        fprintf(stderr, "Couldn't open KMMD control file\n");
        exit( 1 );
    }
    DB_fileList = NULL;
    ref_pdb     = NULL;
    dh_atnames  = NULL;
    resdel      = NULL;
    use_dh      = NULL;
    scan_json_tags( f, &DB_fileList, &ref_pdb,
                  &dh_atnames,   &resdel,  &use_dh, &sig2);
    fclose(f);

    return 0;
}
#endif

#endif
