#ifdef GTI

#include "gti_def.h"
#include <stdexcept>   
#include <assert.h>
#include <math.h>
#include <string>
#include "gti_gpu.h"
#include "gti_gpuContext.h"
#include "hip_definitions.h"

//---------------------------------------------------------------------------------------------
// gti_gpuContext: constructor function for the central gpuContext class objects
//---------------------------------------------------------------------------------------------
gti_gpuContext::gti_gpuContext()
{
  pbTLNeedNewNeighborList = NULL;
  pbTIList = NULL;               // Orig TI list
  pbSCList = NULL;               // Orig TI list
  pbSCBATList = NULL;            // Orig TI list
  pbTIAtomList = NULL;           // TI atom list
  pbTISoftcoreList = NULL;       // TI softcore list
  pbUniqueTIAtomList = NULL;     // Unique TI atom list 
  pbTICommonPair = NULL;         // TI atom pair (region #1, region #2)

  pbTIPotEnergyBuffer = NULL;  // Energy accumulation buffer
  pbTIKinEnergyBuffer = NULL;  // Energy accumulation buffer
  pbTIForce = NULL; // Force accumulators = NULL; 5 sets

  pbPMEChargeGrid = NULL;      // Atom charge grid
  pbOrigAtomCharge = NULL;       // Atom charge storage
  pbFFTChargeGrid = NULL;        // Atom charge on FFT grid

  pbNumberTINBEntries = NULL;
  pbTINBList = NULL;

  pbTINb141 = NULL;
  pbTINb142 = NULL;
  pbTINb14ID = NULL;

  pbTILJTerm = NULL;

  pbTIac = NULL;
  pbTIico = NULL;
  pbTIcn = NULL;
  pbTIDcoef = NULL; // C4PairwiseCUDA
  pbTIDvalue = NULL; // C4PairwiseCUDA
  pbTISigEps = NULL;
  pbTISigMN = NULL;

  pbTIBond = NULL;         // Bond Kr, Req
  pbTIBondID = NULL;       // Bond i, j
  pbTIBondType = NULL;     // Bond TI Type  
  pbTIBondAngle = NULL;    // Bond Angle Ka, Aeq
  pbTIBondAngleID = NULL;  // Bond Angle i, j, k = NULL; regin
  pbTIBondAngleType = NULL;// Dihedral TI Type  
  pbTIDihedral1 = NULL;    // Dihedral gmul, pn
  pbTIDihedral2 = NULL;    // Dihedral pk, gamc
  pbTIDihedral3 = NULL;    // Dihedral gams
  pbTIDihedralID = NULL;   // Dihedral i, j, k, l
  pbTIDihedralType = NULL; // Dihedral TI Type  

  pbTINMRDistancesR = NULL;
  pbTINMRDistancesK = NULL;
  pbTINMRDistancesID = NULL;

  pbTINMRAngleR = NULL;
  pbTINMRAngleK = NULL;
  pbTINMRAngleID = NULL;

  pbTINMRDihedralR = NULL;
  pbTINMRDihedralK = NULL;
  pbTINMRDihedralID = NULL;

  // READ Data section
  pbREAFList = NULL;
  pbREAFAtomList = NULL;

  pbNumberREAFNbEntries = NULL;

  pbREAFNbList = NULL;
  pbREAFNb141 = NULL;
  pbREAFNb142 = NULL;
  pbREAFNb14ID = NULL;

  pbREAFLJTerm = NULL;

  pbREAFac = NULL;
  pbREAFico = NULL;
  pbREAFcn = NULL;
  pbREAFSigEps = NULL;
  pbREAFSigMN = NULL;

  pbREAFBondAngle = NULL;    // Bond Angle Ka, Aeq
  pbREAFBondAngleID = NULL;  // Bond Angle i, j, k = NULL; regin
  pbREAFBondAngleType = NULL;// Dihedral TI Type  
  pbREAFDihedral1 = NULL;    // Dihedral gmul, pn
  pbREAFDihedral2 = NULL;    // Dihedral pk, gamc
  pbREAFDihedral3 = NULL;    // Dihedral gams
  pbREAFDihedralID = NULL;   // Dihedral i, j, k, l
  pbREAFDihedralType = NULL; // Dihedral TI Type  

  // MBAR
  pbMBAREnergy = NULL;  // energy values for states
  pbMBARLambda = NULL;  // lambda values for states
  pbMBARWeight = NULL;  // weight values for states 

  // LJ1264 Data section
  pbNumberLJ1264NBEntries = NULL;
  pbLJ1264AtomList = NULL;
  pbLJ1264NBList = NULL;

  //RMSD fit
  pbMatrix = NULL;
  pbResult = NULL;
  pbWork = NULL;
  pbInfo = NULL;

  // pLJ1264 Data section C4PairwiseCUDA
  pbNumberpLJ1264NBEntries = NULL;
  pbpLJ1264AtomList = NULL;
  pbpLJ1264NBList = NULL;

  // LJ1264pLJ1264 Data section C4PairwiseCUDA2023
  pbNumberLJ1264pLJ1264NBEntries = NULL;
  pbLJ1264pLJ1264AtomList = NULL;
  pbLJ1264pLJ1264NBList = NULL;

}

//---------------------------------------------------------------------------------------------
// ~gti_gpuContext: destructor function for the central gpuContext class objects
//---------------------------------------------------------------------------------------------
gti_gpuContext::~gti_gpuContext()
{
  CleanUp();
}

//---------------------------------------------------------------------------------------------
// UpdateSimulationConst: update the gpuContext information on the device with new host input
//---------------------------------------------------------------------------------------------
void gti_gpuContext::UpdateSimulationConst()
{
  icc_updateAllSimulationConst(this);
}

//---------------------------------------------------------------------------------------------
// CleanUp: garbage dump function
//---------------------------------------------------------------------------------------------
void gti_gpuContext::CleanUp(){
}

unsigned gti_gpuContext::GetMatchAtomID(const uint myID) const {

  int result = myID;

  if (sim.numberTICommonPairs > 0 && pbTICommonPair && pbTICommonPair->_pSysData) {

    for (unsigned i = 0; i < sim.numberTICommonPairs; i++) {
      if (pbTICommonPair->_pSysData[i].x == myID) {
        result = pbTICommonPair->_pSysData[i].y; break;
      }
      if (pbTICommonPair->_pSysData[i].y == myID) {
        result = pbTICommonPair->_pSysData[i].x; break;
      }
    }
  }

  return result;
}

//---------------------------------------------------------------------------------------------
// InitMDParameters: 
//
// Arguments:
//   vlimit:  
//---------------------------------------------------------------------------------------------
void gti_gpuContext::InitMDParameters(double vlimit)
{

  multiStream = false;

  mainStream = NULL; MDStream = NULL; memStream = NULL; TIStream = NULL;

  if (multiStream) {
    cudaStreamCreateWithFlags(&(TIStream), cudaStreamNonBlocking);

    cudaEventCreate(&event_IterationBegin, CU_EVENT_DISABLE_TIMING);
    cudaEventCreate(&event_IterationDone, CU_EVENT_DISABLE_TIMING);
    cudaEventCreate(&event_PMEDone, CU_EVENT_DISABLE_TIMING);
    cudaEventCreate(&event_TIClear, CU_EVENT_DISABLE_TIMING);
    cudaEventCreate(&event_TIDone, CU_EVENT_DISABLE_TIMING);
  }

  CleanUp();
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, gpu_device_id);
  major = deviceProp.major;
  minor = deviceProp.minor;
  sim.needVirial = ((sim.ntp > 0) && (sim.barostat == 1));
  pbTLNeedNewNeighborList = std::unique_ptr<GpuBuffer<bool>>(new GpuBuffer<bool>(1));
  pbTLNeedNewNeighborList->_pSysData[0] = false;
  pbTLNeedNewNeighborList->Upload();
  pbTLNeedNewNeighborList->Download();
  sim.pTLNeedNewNeighborList = pbTLNeedNewNeighborList->_pDevData;
  if (vlimit > 0) {
    sim.bUseVlimit = true;
    sim.vlimit = vlimit;
  }
}

//---------------------------------------------------------------------------------------------
// InitAccumulators:
//---------------------------------------------------------------------------------------------
void gti_gpuContext::InitAccumulators()
{
  unsigned NBMultiplier = sim.needVirial ? 2 : 1;

  if (sim.igamd > 0 ) {
    // the second 7-set is for NB
    pbTIForce = std::unique_ptr< GpuBuffer<PMEAccumulator> > (new GpuBuffer<PMEAccumulator>(sim.stride3 * 7 * NBMultiplier));
    for (unsigned i = 0; i < sim.stride3 * 7 * NBMultiplier; i++) pbTIForce->_pSysData[i] = 0;
  } else {
  // the second 5-set is for NB
    pbTIForce = std::unique_ptr< GpuBuffer<PMEAccumulator> > (new GpuBuffer<PMEAccumulator>(sim.stride3 * 5 * NBMultiplier));
  for (unsigned i = 0; i < sim.stride3 * 5 * NBMultiplier; i++) pbTIForce->_pSysData[i] = 0;
  }

  pbTIForce->Upload();

  sim.pTIForce = pbTIForce->_pDevData;
  if (sim.igamd > 0 ) {
    sim.pTINBForce = (sim.needVirial) ? sim.pTIForce + sim.stride3 * 7 :
                                              sim.pTIForce;
  } else {
  sim.pTINBForce = (sim.needVirial) ? sim.pTIForce + sim.stride3 * 5 :
                                              sim.pTIForce;
  }

  for (unsigned i = 0; i < 3; i++) {
    sim.pTIForceX[i] = sim.pTIForce + sim.stride3*i;
    sim.pTIForceY[i] = sim.pTIForce + sim.stride3*i + sim.stride;
    sim.pTIForceZ[i] = sim.pTIForce + sim.stride3*i + sim.stride2;

    if (i<2){
      sim.pTISCForceX[i] = sim.pTIForceX[i] + sim.stride3 * 3;
      sim.pTISCForceY[i] = sim.pTIForceY[i] + sim.stride3 * 3;
      sim.pTISCForceZ[i] = sim.pTIForceZ[i] + sim.stride3 * 3;

     if (sim.igamd > 0 ) {
      sim.pGaMDTIForceX[i] = sim.pTIForceX[i] + sim.stride3 * 5;
      sim.pGaMDTIForceY[i] = sim.pTIForceY[i] + sim.stride3 * 5;
      sim.pGaMDTIForceZ[i] = sim.pTIForceZ[i] + sim.stride3 * 5;
     }
    }
  }

  for (unsigned i = 0; i < 3; i++) {
    sim.pTINBForceX[i] = sim.pTINBForce + sim.stride3*i;
    sim.pTINBForceY[i] = sim.pTINBForce + sim.stride3*i + sim.stride;
    sim.pTINBForceZ[i] = sim.pTINBForce + sim.stride3*i + sim.stride2;

    if (i<2) {
      sim.pTINBSCForceX[i] = sim.pTINBForceX[i] + sim.stride3 * 3;
      sim.pTINBSCForceY[i] = sim.pTINBForceY[i] + sim.stride3 * 3;
      sim.pTINBSCForceZ[i] = sim.pTINBForceZ[i] + sim.stride3 * 3;

     if (sim.igamd > 0 ) {
      sim.pGaMDTINBForceX[i] = sim.pTINBForceX[i] + sim.stride3 * 5;
      sim.pGaMDTINBForceY[i] = sim.pTINBForceY[i] + sim.stride3 * 5;
      sim.pGaMDTINBForceZ[i] = sim.pTINBForceZ[i] + sim.stride3 * 5;
     }
    }
  }
  // TIEnergyBuffer
  // 0, 1 : linear energy of region 0 and 1
  // 2: the common energy
  // 3, 4 : lambda-dependent of region 0 and 1
  // 5, 6 : the SC-only energy of region 0 and 1
  pbTIPotEnergyBuffer = std::unique_ptr< GpuBuffer<unsigned long long int> > (new GpuBuffer<unsigned long long int>(sim.GPUPotEnergyTerms * sim.TIEnergyBufferMultiplier));
  
  pbTIPotEnergyBuffer->Upload();
  sim.pTIPotEnergyBuffer = pbTIPotEnergyBuffer->_pDevData;

  unsigned long long int* energyPos[sim.TIEnergyBufferMultiplier];
  for (unsigned i = 0; i < sim.TIEnergyBufferMultiplier; i++) {
    energyPos[i] = sim.pTIPotEnergyBuffer + (sim.GPUPotEnergyTerms * i);
  }

  static const int DLShift = gti_simulationConst::TIEnergyDLShift;
  static const int SCShift = gti_simulationConst::TIEnergySCShift;
  static const int ExShift = gti_simulationConst::TIExtraShift;

  for (unsigned i = 0; i < 3; i++) {
    sim.pTILJ[i] = energyPos[i] + 1;
    sim.pTIEBond[i] = energyPos[i] + 3;
    sim.pTIEAngle[i] = energyPos[i] + 4;
    sim.pTIEDihedral[i] = energyPos[i] + 5;
    sim.pTIRestDist[i] = energyPos[i] + 14;
    sim.pTIRestAng[i] = energyPos[i] + 15;
    sim.pTIRestTor[i] = energyPos[i] + 16;
    sim.pTIRestRMSD[i] = energyPos[i] + 17;

    sim.pTIEE14CC[i] = energyPos[i] + 6;
    sim.pTILJ14[i] = energyPos[i] + 7;
    sim.pTIEER[i] = energyPos[i] + 9;
    sim.pTIEECC[i] = energyPos[i] + 10;
    sim.pTIEESC[i] = energyPos[i] + ExShift;
    sim.pTIEE14SC[i] = energyPos[i] + ExShift +1;

    sim.pTISCLJ_DL[i] = energyPos[i + DLShift] + 1;
    sim.pTISCEBond_DL[i] = energyPos[i + DLShift] + 3;
    sim.pTISCEAngle_DL[i] = energyPos[i + DLShift] + 4;
    sim.pTISCEDihedral_DL[i] = energyPos[i + DLShift] + 5;
    sim.pTISCEECC_DL[i] = energyPos[i+ DLShift] + 10;
    sim.pTISCEESC_DL[i] = energyPos[i + DLShift] + ExShift;
    sim.pTIEE14CC_DL[i] = energyPos[i + DLShift] + 7;
    sim.pTIEE14SC_DL[i] = energyPos[i + DLShift] + ExShift +1;
    sim.pTIRestDist_DL[i] = energyPos[i + DLShift] + 14;
    sim.pTIRestAng_DL[i] = energyPos[i+ DLShift] + 15;
    sim.pTIRestTor_DL[i] = energyPos[i+ DLShift] + 16;
    sim.pTIRestRMSD_DL[i] = energyPos[i + DLShift] + 17;

    sim.pTISCLJ[i] = energyPos[i+ SCShift] + 1;
    sim.pTISCBond[i] = energyPos[i+ SCShift] + 3;
    sim.pTISCAngle[i] = energyPos[i+ SCShift] + 4;
    sim.pTISCDihedral[i] = energyPos[i+ SCShift] + 5;
    sim.pTISCEE14[i] = energyPos[i + SCShift] + 6;
    sim.pTISCLJ14[i] = energyPos[i + SCShift] + 7;
    sim.pTISCEED[i] = energyPos[i+ SCShift] + 10;

    for (unsigned j = 0; j < gti_simulationConst::TINumberSpecialTerms; j++) {
      sim.pTISpecialTerms[j][i] = energyPos[i] + gti_simulationConst::TISpecialShift+j;
    }

    sim.pV11[i] = energyPos[i] + VIRIAL_OFFSET;
    sim.pV22[i] = energyPos[i] + VIRIAL_OFFSET+1;
    sim.pV33[i] = energyPos[i] + VIRIAL_OFFSET+2;
  }

  pbTIKinEnergyBuffer = std::unique_ptr< GpuBuffer<unsigned long long int> > (new GpuBuffer<unsigned long long int>(sim.GPUKinEnergyTerms * 3));
  pbTIKinEnergyBuffer->Upload();
  sim.pTIKinEnergyBuffer = pbTIKinEnergyBuffer->_pDevData;
}

//---------------------------------------------------------------------------------------------
// InitTIParameters:
//
// Arguments:
//   ti_latm_lst:  
//   ti_lst:       
//   ti_sc_lst:    
//   lambda:       
//   klambda:      
//   scalpha:      
//   scbeta:  
//   addSC14:     
//---------------------------------------------------------------------------------------------
void gti_gpuContext::InitTIParameters(int ti_latm_lst[][3], int ti_lst[][3],
  int ti_sc_lst[], int ti_sc_bat_lst[], gti_simulationConst::TypeCtrlVar ctrlVar)
{
  // Setup SC parameters
  sim.TISCAlpha = ctrlVar.scalpha;
  sim.TISCBeta = ctrlVar.scbeta;
  sim.TISCGamma = ctrlVar.scgamma;
  sim.TIk = ctrlVar.klambda;
  sim.addSC = ctrlVar.addSC;
  sim.eleGauss = ctrlVar.eleGauss;
  sim.tiCut = ctrlVar.tiCut;
  sim.eleSC = ctrlVar.eleSC;
  sim.vdwSC = ctrlVar.vdwSC;
  sim.autoAlpha = ctrlVar.autoAlpha;
  sim.scaleBeta = ctrlVar.scaleBeta;
  sim.eleExp = ctrlVar.eleExp;
  sim.vdwExp = ctrlVar.vdwExp;
  sim.vdwCap = ctrlVar.vdwCap;
  sim.cut_sc = ctrlVar.cut_sc;
  sim.cut_sc0 = ctrlVar.cut_sc0;
  sim.cut_sc1 = ctrlVar.cut_sc1;

  // TI part
  pbTIList = std::unique_ptr< GpuBuffer<unsigned> > (new GpuBuffer<unsigned>(sim.stride3));
  pbSCList = std::unique_ptr< GpuBuffer<unsigned> > (new GpuBuffer<unsigned>(sim.stride));
  pbSCBATList = std::unique_ptr< GpuBuffer<unsigned> > (new GpuBuffer<unsigned>(sim.stride));

  uint4 *temp = new uint4[gti_simulationConst::MaxNumberTIAtom];
  unsigned* tempSC = new unsigned[gti_simulationConst::MaxNumberTIAtom*2];

  unsigned counterCommon[2] = {0, 0};
  unsigned counterSC[2] = {0, 0};
  unsigned tempShift = gti_simulationConst::MaxNumberTIAtom;

  for (unsigned i = 0; i < sim.atoms; i++) {
    pbSCList->_pSysData[i] = ti_sc_lst[i];
    pbSCBATList->_pSysData[i] = ti_sc_bat_lst[i];
  }
  unsigned TISize = 0;
  for (unsigned j = 0; j < 3; j++){
    for (unsigned i = 0; i < sim.atoms; i++){
      pbTIList->_pSysData[i + j*sim.stride] = ti_lst[i][j];
      if (j < 2 && ti_lst[i][j]>0) {
        temp[TISize].x = i;
        temp[TISize].y = j;
        if (ti_sc_lst[i] > 0) {
          temp[TISize].z = 1;
          tempSC[counterSC[j] + tempShift*j] = i;
          counterSC[j]++;
        } else {
          temp[TISize].z = 0;
          counterCommon[j]++;
        }
        TISize++;
        if (TISize > gti_simulationConst::MaxNumberTIAtom) {
          std::string e("Too many TI atoms--current limit: 500 atoms");
          throw std::runtime_error(e);
        }
      }
    }
  }

  // pmemd orig list
  pbTIList->Upload();
  pbSCList->Upload();
  pbSCBATList->Upload();
  sim.pTIList = pbTIList->_pDevData;
  sim.pSCList = pbSCList->_pDevData;
  sim.pSCBATList = pbSCList->_pDevData;

  // TI atom list
  pbTIAtomList = std::unique_ptr< GpuBuffer<uint4> > (new GpuBuffer<uint4>(TISize));
  for (unsigned i = 0; i < TISize; i++) pbTIAtomList->_pSysData[i] = temp[i];
  delete[] temp;
  pbTIAtomList->Upload();
  sim.numberTIAtoms = TISize;
  sim.pTIAtomList = pbTIAtomList->_pDevData;

  // TI common pairs
  assert(counterCommon[0] == counterCommon[1]);
  pbTICommonPair = std::unique_ptr< GpuBuffer<int4> > (new GpuBuffer<int4>(counterCommon[0]));
  for (unsigned i = 0; i < counterCommon[0]; i++) {
    pbTICommonPair->_pSysData[i].x = ti_latm_lst[i][0]-1;
    pbTICommonPair->_pSysData[i].y = ti_latm_lst[i][1]-1;
    pbTICommonPair->_pSysData[i].z = -1;
  }
  pbTICommonPair->Upload();
  sim.numberTICommonPairs = counterCommon[0];
  sim.pTICommonPair = pbTICommonPair->_pDevData;

  // TI SC list
  pbTISoftcoreList = std::unique_ptr< GpuBuffer<unsigned> > (new GpuBuffer<unsigned>(counterSC[0] + counterSC[1]));
  for (unsigned j = 0; j < 2; j++) {
    for (unsigned i = 0; i < counterSC[j]; i++)
      pbTISoftcoreList->_pSysData[i + j* counterSC[0] ] = tempSC[i + tempShift*j];
  }
  delete[] tempSC;
  pbTISoftcoreList->Upload();
  sim.softcoreAtomShift = counterSC[0];
  sim.numberSoftcoreAtoms[0] = counterSC[0];
  sim.numberSoftcoreAtoms[1] = counterSC[1];
  sim.pTISoftcoreList = pbTISoftcoreList->_pDevData;

  // Unique list
  sim.numberUniqueTIAtoms = (sim.numberTICommonPairs + sim.numberSoftcoreAtoms[0] +
                             sim.numberSoftcoreAtoms[1]);
  pbUniqueTIAtomList = std::unique_ptr< GpuBuffer<unsigned> > (new GpuBuffer<unsigned>(sim.numberUniqueTIAtoms));

  for (unsigned i = 0; i < sim.numberTICommonPairs; i++) {
    pbUniqueTIAtomList->_pSysData[i] = pbTICommonPair->_pSysData[i].x;
  }
  for (unsigned i = sim.numberTICommonPairs;
       i < (sim.numberTICommonPairs + sim.numberSoftcoreAtoms[0] +
            sim.numberSoftcoreAtoms[1]) ; i++) {
    pbUniqueTIAtomList->_pSysData[i] = pbTISoftcoreList->_pSysData[i- sim.numberTICommonPairs];
  }
  pbUniqueTIAtomList->Upload();
  sim.pUniqueTIAtomList = pbUniqueTIAtomList->_pDevData;

  
  if (!pbOrigAtomCharge) {
    pbOrigAtomCharge = std::unique_ptr< GpuBuffer<PMEFloat> >(new GpuBuffer<PMEFloat>(sim.stride));
    sim.pOrigAtomCharge = pbOrigAtomCharge->_pDevData;
  }
}

void gti_gpuContext::InitREADParameters(int reaf_mode, double reaf_tau, int addRE, int read_Atom_list[][2]) {

  assert(reaf_mode == 0 || reaf_mode == 1);

  sim.reafMode = reaf_mode;
  sim.current_tau = reaf_tau;
  sim.addRE = addRE;
  unsigned readSize=0;
  int tempList[gti_simulationConst::MaxNumberTIAtom];

  pbREAFList = std::unique_ptr< GpuBuffer<unsigned> >(new GpuBuffer<unsigned>(sim.stride));

  for (unsigned i = 0; i < sim.atoms; i++) {
    pbREAFList->_pSysData[i] = read_Atom_list[i][reaf_mode];
    if (read_Atom_list[i][reaf_mode] > 0) {
      tempList[readSize++] = i;
    }
  }

  if (readSize == 0) {
    sim.reafMode = -1;
    return;
  } 

  sim.numberREAFAtoms = readSize;

  pbREAFAtomList = std::unique_ptr< GpuBuffer<unsigned> >(new GpuBuffer<unsigned>(readSize));
  for (unsigned i = 0; i < readSize; i++) {
    pbREAFAtomList->_pSysData[i] = tempList[i];
  }

  pbREAFList->Upload();
  sim.pREAFList = pbREAFList->_pDevData;
  
  pbREAFAtomList->Upload();
  sim.pREAFAtomList = pbREAFAtomList->_pDevData;

  if (!pbOrigAtomCharge) {
    pbOrigAtomCharge = std::unique_ptr< GpuBuffer<PMEFloat> >(new GpuBuffer<PMEFloat>(sim.stride));
    sim.pOrigAtomCharge = pbOrigAtomCharge->_pDevData;
  }

}


//---------------------------------------------------------------------------------------------
// InitGaMDTIParameters:
//
// Arguments:
//   ti_latm_lst:  
//   ti_lst:       
//   ti_sc_lst:    
//   lambda:       
//   klambda:      
//   scalpha:      
//   scbeta:  
//   addSC14:     
//---------------------------------------------------------------------------------------------
void gti_gpuContext::InitGaMDTIParameters(int ti_latm_lst[][3], int ti_lst[][3], int ti_sc_lst[], int ti_sc_bat_lst[],
                                      double lambda, unsigned klambda, double scalpha,
                                      double scbeta, double scgamma, int tiCut,
                                      int addSC, int eleGauss, int eleSC, int vdwSC, double vdwCap)
{
  // TI part
  pbTIList = std::unique_ptr< GpuBuffer<unsigned> > (new GpuBuffer<unsigned>(sim.stride3));
  pbSCList = std::unique_ptr< GpuBuffer<unsigned> > (new GpuBuffer<unsigned>(sim.stride));
  pbSCBATList = std::unique_ptr< GpuBuffer<unsigned> > (new GpuBuffer<unsigned>(sim.stride));

  uint4 *temp = new uint4[sim.atoms];
  unsigned* tempSC = new unsigned[sim.atoms*2];

  unsigned counterCommon[2] = {0, 0};
  unsigned counterSC[2] = {0, 0};
  unsigned tempShift = sim.atoms;

  for (unsigned i = 0; i < sim.atoms; i++) {
    pbSCList->_pSysData[i] = ti_sc_lst[i];
    pbSCBATList->_pSysData[i] = ti_sc_bat_lst[i];
  }
  unsigned TISize = 0;
  for (unsigned j = 0; j < 3; j++){
    for (unsigned i = 0; i < sim.atoms; i++){
      pbTIList->_pSysData[i + j*sim.stride] = ti_lst[i][j];
      if (j < 2 && ti_lst[i][j]>0) {
        temp[TISize].x = i;
        temp[TISize].y = j;
        if (ti_sc_lst[i] > 0) {
          temp[TISize].z = 1;
          tempSC[counterSC[j] + tempShift*j] = i;
          counterSC[j]++;
        } else {
          temp[TISize].z = 0;
          counterCommon[j]++;
        }
        TISize++;
      }
    }
  }

  // pmemd orig list
  pbTIList->Upload();
  pbSCList->Upload();
  pbSCBATList->Upload();
  sim.pTIList = pbTIList->_pDevData;
  sim.pSCList = pbSCList->_pDevData;
  sim.pSCBATList = pbSCList->_pDevData;

  // TI atom list
  pbTIAtomList = std::unique_ptr< GpuBuffer<uint4> > (new GpuBuffer<uint4>(TISize));
  for (unsigned i = 0; i < TISize; i++) pbTIAtomList->_pSysData[i] = temp[i];
  delete[] temp;
  pbTIAtomList->Upload();
  sim.numberTIAtoms = TISize;
  sim.pTIAtomList = pbTIAtomList->_pDevData;

  // TI common pairs
  assert(counterCommon[0] == counterCommon[1]);
  pbTICommonPair = std::unique_ptr< GpuBuffer<int4> > (new GpuBuffer<int4>(counterCommon[0]));
  for (unsigned i = 0; i < counterCommon[0]; i++) {
    pbTICommonPair->_pSysData[i].x = ti_latm_lst[i][0]-1;
    pbTICommonPair->_pSysData[i].y = ti_latm_lst[i][1]-1;
  }
  pbTICommonPair->Upload();
  sim.numberTICommonPairs = counterCommon[0];
  sim.pTICommonPair = pbTICommonPair->_pDevData;

  // TI SC list
  pbTISoftcoreList = std::unique_ptr< GpuBuffer<unsigned> > (new GpuBuffer<unsigned>(counterSC[0] + counterSC[1]));
  for (unsigned j = 0; j < 2; j++) {
    for (unsigned i = 0; i < counterSC[j]; i++)
      pbTISoftcoreList->_pSysData[i + j* counterSC[0] ] = tempSC[i + tempShift*j];
  }
  delete[] tempSC;
  pbTISoftcoreList->Upload();
  sim.softcoreAtomShift = counterSC[0];
  sim.numberSoftcoreAtoms[0] = counterSC[0];
  sim.numberSoftcoreAtoms[1] = counterSC[1];
  sim.pTISoftcoreList = pbTISoftcoreList->_pDevData;

  // Unique list
  sim.numberUniqueTIAtoms = (sim.numberTICommonPairs + sim.numberSoftcoreAtoms[0] +
                             sim.numberSoftcoreAtoms[1]);
  pbUniqueTIAtomList = std::unique_ptr< GpuBuffer<unsigned> > (new GpuBuffer<unsigned>(sim.numberUniqueTIAtoms));

  for (unsigned i = 0; i < sim.numberTICommonPairs; i++) {
    pbUniqueTIAtomList->_pSysData[i] = pbTICommonPair->_pSysData[i].x;
  }
  for (unsigned i = sim.numberTICommonPairs;
       i < (sim.numberTICommonPairs + sim.numberSoftcoreAtoms[0] +
            sim.numberSoftcoreAtoms[1]) ; i++) {
    pbUniqueTIAtomList->_pSysData[i] = pbTISoftcoreList->_pSysData[i- sim.numberTICommonPairs];
  }
  pbUniqueTIAtomList->Upload();
  sim.pUniqueTIAtomList = pbUniqueTIAtomList->_pDevData;

  // Setup SC parameters
  sim.TISCAlpha = scalpha;
  sim.TISCBeta = scbeta;
  sim.TISCGamma = abs(scgamma);
  sim.addSC = addSC;
  sim.eleGauss = eleGauss;
  sim.eleSC = (eleSC<0) ? 0 : eleSC;
  sim.tiCut = tiCut;
  sim.vdwSC = (vdwSC<0) ? 0 : vdwSC;
  sim.vdwCap = vdwCap;
  
  pbOrigAtomCharge = std::unique_ptr< GpuBuffer<PMEFloat> > (new GpuBuffer<PMEFloat>(sim.stride));
  sim.pOrigAtomCharge = pbOrigAtomCharge->_pDevData;
}

//---------------------------------------------------------------------------------------------
// InitMBARParameters:
//
// Arguments:
//   nMBAR_state:  
//   MBAR_lambda:  
//---------------------------------------------------------------------------------------------
void gti_gpuContext::InitMBARParameters(unsigned nMBAR_state, double MBAR_lambda[][2], int rangeMBAR)
{  
  sim.needMBAR = true;
  sim.nMBARStates = nMBAR_state;  // Number of MBAR states 
  sim.rangeMBAR = rangeMBAR;

  double* Lam[2] = { new double[sim.nMBARStates], new double[sim.nMBARStates] };
  double* weight = new double[sim.nMBARStates];
  double* dWeight = new double[sim.nMBARStates];
  for (unsigned l = 0; l < sim.nMBARStates; l++) {
    Lam[0][l] = 1.0 - MBAR_lambda[l][0];
    Lam[1][l] = 1.0 - MBAR_lambda[l][1];
  }

  pbMBAREnergy = std::unique_ptr< GpuBuffer<unsigned long long int> >(new GpuBuffer<unsigned long long int>(sim.nMBARStates*3* LambdaSchedule::TypeTotal)); // 3=2+1: two TI regions and one region-indep;
  pbMBARLambda = std::unique_ptr< GpuBuffer<PMEFloat> > (new GpuBuffer<PMEFloat>(2*sim.nMBARStates));
  pbMBARWeight = std::unique_ptr< GpuBuffer<PMEFloat> > (new GpuBuffer<PMEFloat>(2 * sim.nMBARStates));
  for (unsigned l = 0; l < sim.nMBARStates; l++) {
    pbMBARLambda->_pSysData[l] = Lam[0][l];
    pbMBARLambda->_pSysData[l+ sim.nMBARStates] = Lam[1][l];
  }

  LambdaSchedule& schd = LambdaSchedule::GetReference();
  schd.GetWeight(Schedule::TypeRestBA, sim.nMBARStates, Lam[0], weight, dWeight);
  for (unsigned l = 0; l < sim.nMBARStates; l++) {
    pbMBARWeight->_pSysData[l] = weight[l];
  }
  schd.GetWeight(Schedule::TypeRestBA, sim.nMBARStates, Lam[1], weight, dWeight);
  for (unsigned l = 0; l < sim.nMBARStates; l++) {
    pbMBARWeight->_pSysData[l + sim.nMBARStates] = weight[l];
  }

  for (unsigned l = 0; l < sim.nMBARStates * 2 * Schedule::TypeTotal ; l++) {
    pbMBAREnergy->_pSysData[l] = 0;
  }

  pbMBARLambda->Upload();
  sim.pMBARLambda = pbMBARLambda->_pDevData;

  pbMBARWeight->Upload();
  sim.pMBARWeight = pbMBARWeight->_pDevData;

  pbMBAREnergy->Upload();
  sim.pMBAREnergy = pbMBAREnergy->_pDevData;

  delete[] Lam[0];
  delete[] Lam[1];
  delete[] weight;
  delete[] dWeight;
}

//---------------------------------------------------------------------------------------------
// UpdateTICharges:
// 
// Arguments:
//   charge:  
//---------------------------------------------------------------------------------------------
void gti_gpuContext::UpdateCharges(double charge[])
{
  if (!pbAtomCharge || !pbOrigAtomCharge) return;

  for (unsigned i = 0; i < sim.atoms; i++) {
    pbAtomCharge->_pSysData[i] = charge[i];
    pbOrigAtomCharge->_pSysData[i] = charge[i];
  }
  pbAtomCharge->Upload();
  pbOrigAtomCharge->Upload();

  // Add updated charges to the ee-14 term holder
  if (sim.numberTI14NBEntries > 0) {
    for (unsigned i = 0; i < sim.numberTI14NBEntries; i++) {
      double q0 = pbAtomCharge->_pSysData[pbTINb14ID->_pSysData[i].x];
      double q1 = pbAtomCharge->_pSysData[pbTINb14ID->_pSysData[i].y];
      pbTINb141->_pSysData[i].x *= q0 * q1;
    }
    if (pbTINb141) pbTINb141->Upload();
  }

  if (sim.numberREAF14NBEntries > 0) {
    for (unsigned i = 0; i < sim.numberREAF14NBEntries; i++) {
      double q0 = pbAtomCharge->_pSysData[pbREAFNb14ID->_pSysData[i].x];
      double q1 = pbAtomCharge->_pSysData[pbREAFNb14ID->_pSysData[i].y];
      pbREAFNb141->_pSysData[i].x *= q0 * q1;
    }
    if (pbREAFNb141) pbREAFNb141->Upload();
  }
}

//---------------------------------------------------------------------------------------------
// UpdateTIMasses:
// 
// Arguments:
//   mass:  
//---------------------------------------------------------------------------------------------
void gti_gpuContext::UpdateTIMasses(double mass[])
{
  for (int i = 0; i < sim.atoms; i++) {
    pbAtomMass->_pSysData[i] = mass[i];
    if (mass[i] > 0.0) {
      pbAtomMass->_pSysData[i +sim.stride] = (PMEDouble)(1.0 / mass[i]);
    }
    else {
      pbAtomMass->_pSysData[i +sim.stride] = (PMEDouble)0.0;
    }
    pbImageMass->_pSysData[i] = pbAtomMass->_pSysData[i];
    pbImageMass->_pSysData[i + sim.stride] = pbAtomMass->_pSysData[i + sim.stride];
  }
  pbAtomMass->Upload(); 
  pbImageMass->Upload();
}


//void gti_gpuContext::InitRMSD(int RMSDMask[5], int ti_lst[][3],
#endif
