#include "kmmd_types.h"

typedef struct NbData {
  int    id;
  float  distance;
} NbData;

void buildNeighbourGraph(KMMDCrd_t    *dh_coords, 
                           KMMDEneTrn_t *denes, 
                           int           n_snaps, 
                           int           n_dims,
                           KMMDCrd_t     test_sig2);
