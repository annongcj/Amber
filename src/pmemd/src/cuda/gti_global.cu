//Place to hold globally accessible variables
// Note: only applicable to the rdc=true case

#if defined(GTI) && defined(__CUDACC_RDC__)
#  include "simulationConst.h"
  __device__ __constant__ simulationConst cSim;
#endif
