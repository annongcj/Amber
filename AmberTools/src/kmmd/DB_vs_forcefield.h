#include "kmmd_types.h"

int  DB_snap_forcefield_enes(char  *DB_subdir_path, 
                     KMMDEneTrn_t **snap_ff_enes, 
                     char         **DB_fnames,
                     char          *topfile_path,
                     int            n_DB_snaps);




int  DB_clean_highEneSnaps(KMMDCrd_t    **dh_coords,
                           KMMDCrd_t    **preWeights,
                           KMMDCrd_t    **weight_gradients_TI,
                           KMMDEneTrn_t **denes,
                           char        ***DB_fnames,
                           int            n_DB_snaps,
                           int            n_dh_per_snap);
          


