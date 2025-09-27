
#include "kmmd_types.h"

int    get_pdb(const char *path, KMMDCrd_t **coords, char ***atnames, int **resids);
void   write_pdb(char *path, FILE *g, KMMDCrd_t *coords, int n_atoms);
int*   assign_dh_atoms(HashTable  *dh_names,\
                       HashTable  *dh_resdels,\
                       HashTable  *dh_atnames,\
                                   char **atnames,\
                                   int   *resids,\
                                   int    n_ats);
