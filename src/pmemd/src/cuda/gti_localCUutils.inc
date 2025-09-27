#if defined(GTI) && !(defined(__CUDACC_RDC__) || defined(__HIPCC_RDC__))

#include "gti_utils.inc"

__device__ void faster_sincos(PMEDouble a, PMEDouble *sptr, PMEDouble *cptr)
#include "kFastCosSin.h"

#endif
