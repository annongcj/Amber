/* The program body will have implementation only when rdc=true  */
#if defined(GTI)  && defined(__CUDACC_RDC__)
#  include "gti_utils.cuh"
#  include "gpuContext.h"

#  include "simulationConst.h"
CSIM_STO simulationConst cSim;

#  include "gti_utils.inc"

#endif

