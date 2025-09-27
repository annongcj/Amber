//---------------------------------------------------------------------------------------------
// pmemd.cuda TI extension (pmemdGTI)
//
// Proprietary and confidential
// Copyright (C) Tai-Sung Lee & Darrin York
// Laboratory for Biomolecular Simulation Research, Rutgers University
//
// All usage should be under appropriate license terms.  Any copying via any medium,
// modification, extraction, or distribution of this file, is strictly prohibited without the
// express written permission of the authors.
//---------------------------------------------------------------------------------------------

#ifdef GTI

///////////////////////////////////////////////////////////////
//  Functions to be included for both rdc=true and rdc=false cases
///////////////////////////////////////////////////////////////

#include "gpuContext.h"

void UpdateSimulationConst() {

  gpuContext gpu = theGPUContext::GetPointer();

  cudaError_t status;

  status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(simulationConst));
  RTERROR(status, "cudaMemcpyToSymbol: UpdateSimulationConst failed");
}

#endif  /* GTI */
