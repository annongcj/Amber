#ifndef HAVE_KMMD_TYPES
#define HAVE_KMMD_TYPES


#ifdef CUDA
#include "gpuTypes.h"

typedef PMEFloat       KMMDCrd_t;     //type for coordinates
typedef PMEDouble      KMMDFrcAcc_t;  //type for total force
typedef PMEDouble      KMMDEneAcc_t;  //type for total ene
typedef PMEFloat       KMMDFrcTrn_t;  //type for training point individual force
typedef PMEDouble      KMMDEneTrn_t;  //type for training point individual energy
typedef PMEFloat4      KMMDFloat4;    //type for four floats x,y,z,w.  Use for 3-vectors so have nice alignment to cache boundaries.

#else

typedef double  KMMDCrd_t;          
typedef double  KMMDFrcAcc_t;
typedef double  KMMDEneAcc_t;
typedef double  KMMDFrcTrn_t;
typedef double  KMMDEneTrn_t;
typedef struct  KMMDFloat4 {
    float x;
    float y;
    float z;
    float w;
} KMMDFloat4;

#endif //CUDA
#endif //HAVE_KMMD_TYPES

