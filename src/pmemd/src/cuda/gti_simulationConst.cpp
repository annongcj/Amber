#ifdef GTI
#include "gti_simulationConst.h"

//---------------------------------------------------------------------------------------------
// InitSimulationConst: initialize the simulation constants structure within the gpuContext
//---------------------------------------------------------------------------------------------
void gti_simulationConst::InitSimulationConst()
{
  base_simulationConst::InitSimulationConst();

  needVirial = false;
  needMBAR = false;

  pNTPData = NULL;
  pTIList = NULL;
  pSCList = NULL;
  numberTIAtoms = 0;
  numberTICommonPairs = 0;
  doTIHeating = false;

  numberTI14NBEntries = 0;
  numberTIBond = 0;
  numberTIAngle = 0;
  numberTIDihedral = 0;

  numberTINMRDistance = 0;
  numberTINMRAngle = 0;
  numberTINMRDihedral = 0;

  nMBARStates = 0;
  currentMBARState = -1;
  rangeMBAR = -1;
  needMBAR = false;

  numberLJ1264Atoms = 0;
  numberpLJ1264Atoms = 0;  //C4PairwiseCUDA2023
  numberLJ1264pLJ1264Atoms = 0; 
  numberREAFAtoms = 0;
  numberREAF14NBEntries = 0;
  numberREAFAngle = 0;
  numberREAFDihedral = 0;

  eleExp = 2;
  vdwExp = 6;
  
  reafMode=-1;
  addRE=-999;
  eleTauType=-999;
  vdwTauType=-999;
  current_tau=1.0e99;

  number_rmsd_set = 0;
  rmsd_atom_max_count = 0;
}

#endif

