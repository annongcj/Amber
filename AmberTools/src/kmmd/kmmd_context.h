#ifndef _KMMD_HOST_CONTEXT
#define _KMMD_HOST_CONTEXT

#include <stdio.h> 
#include <cstdlib>
#include <math.h>

#ifdef CUDA
   #include "gputypes.h"
   #include "gpuBuffer.h"
#endif

#include "kmmd_types.h"
#include "hash.h"

//utility function to read MLMLD control data
//..defined in KMMDParseJSON.cpp.
int scan_json_tags( FILE *f, 
                    char **DB_fileList, 
                    char **ref_pdb,
                    char **dvdl_logfile_name,
                    int   *acc_dvdl_steps,
                    HashTable **dh_atnames,
                    HashTable **dh_resdel,
                    HashTable **dh_byres, 
                    KMMDCrd_t  *sigma2 );


class kmmdHostContext {

public:
  virtual void Init(char *input_filename);
  
  kmmdHostContext(char *input_filename);
  virtual ~kmmdHostContext();
  
public:
  int           n_kmmd_atoms;
  int           n_dh;
  int           n_trainpts;
  KMMDCrd_t     sigma2;

  int           acc_dvdl_steps, acc_dvdl_stepcount;
  FILE         *dvdl_logfile;
  KMMDEneAcc_t  acc_dvdl; 


#ifdef KMMD_CUDA
  GpuBuffer<int>*            dh_atoms;      // 4*n_dh atoms to take part in KMMD dihedrals
  GpuBuffer<KMMDCrd_t>*      dh_trainpts;   // n_dh * n_trainpts * 2 dihedral cos/sin pairs
  GpuBuffer<KMMDCrd_t>*      preWeights;    // n_trainpts
  GpuBuffer<KMMDCrd_t>*      dweights_dlambda_scale; // n_trainpts or NULL if not used. 
  GpuBuffer<KMMDEneTrn_t>*   trainenes;     // n_trainpts 
  GpuBuffer<KMMDFrcAcc_t>*   kmmd_frc;      // force delta
#else
  int*            dh_atoms;      // 4*n_dh atoms to take part in KMMD dihedrals
  KMMDCrd_t*      dh_trainpts;   // n_dh * n_trainpts * 2 dihedral cos/sin pairs
  KMMDCrd_t*      preWeights;    // n_trainpts 
  KMMDCrd_t*      dweights_dlambda_scale; // n_trainpts or NULL if not used. 
  KMMDEneTrn_t*   trainenes;     // n_trainpts 
  KMMDFrcAcc_t*   kmmd_frc;      // force delta
#endif
};

#endif
