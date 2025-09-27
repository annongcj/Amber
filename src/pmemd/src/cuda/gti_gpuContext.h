#ifndef _GTI_GPU_Context
#define _GTI_GPU_Context
#include <memory>
#include "gputypes.h"
#include "gpuBuffer.h"
#include "base_gpuContext.h"
#include "gti_simulationConst.h"

class gti_gpuContext : public base_gpuContext {

public:

  gti_gpuContext();
  ~gti_gpuContext();
  
  void UpdateSimulationConst();

  void InitMDParameters(double vlimit);
  void InitAccumulators();
  void InitTIParameters(int ti_latm_lst[][3], int ti_lst[][3],
    int ti_sc_lst[], int ti_sc_bat_lst[], gti_simulationConst::TypeCtrlVar ctrlVar);
  void InitREADParameters(int reaf_mode, double reaf_tau, int addRE, int reaf_atom_list[][2]);

  void InitGaMDTIParameters(int ti_latm_lst[][3], int ti_lst[][3], 
    int ti_sc_lst[], int ti_sc_bat_lst[], double lambda,
    unsigned klambda, double scalpha, double scbeta, double scgamma,
    int tiCut, int addSC, int eleGauss, int eleSC, int vdwSC, double dwCap);
  void InitMBARParameters(unsigned nMAR_state, double MBAR_lambda[][2], int rangeMBAR);

  void UpdateCharges(double charge[]);
  void UpdateTIMasses(double mass[]);

  unsigned GetMatchAtomID(const uint myID) const;

  cudaStream_t mainStream, MDStream, memStream;
  cudaStream_t TIStream;
  cudaEvent_t event_PMEDone, event_IterationBegin, event_IterationDone, event_TIClear, event_TIDone;

  // GPU Info
  unsigned int major;
  unsigned int minor;

  // Flags
  bool needVirial;
  bool multiStream;

  // TI Data section
  std::unique_ptr< GpuBuffer<bool>>       pbTLNeedNewNeighborList;
  std::unique_ptr< GpuBuffer<unsigned>>   pbTIList;               // Orig TI list
  std::unique_ptr< GpuBuffer<unsigned>>   pbSCList;               // Orig TI list
  std::unique_ptr< GpuBuffer<unsigned>>   pbSCBATList;            // Orig TI list
  std::unique_ptr< GpuBuffer<uint4>>      pbTIAtomList;           // TI atom list
  std::unique_ptr< GpuBuffer<unsigned>>   pbTISoftcoreList;       // TI softcore list
  std::unique_ptr< GpuBuffer<unsigned>>   pbUniqueTIAtomList;     // Unique TI atom list 
  std::unique_ptr< GpuBuffer<int4>>      pbTICommonPair;         // TI atom pair (region #1, region #2)

  std::unique_ptr< GpuBuffer<unsigned long long int>> pbTIPotEnergyBuffer;  // Energy accumulation buffer
  std::unique_ptr< GpuBuffer<unsigned long long int>> pbTIKinEnergyBuffer;  // Energy accumulation buffer
  std::unique_ptr< GpuBuffer<PMEAccumulator>>     pbTIForce; // Force accumulators; 5 sets
                                                           //   separated by sim.stride*3;
  // 0: TI region 0; 1: TI region 1; 2: TI common; 3: SC region 0; 4: SC reion 1
  // 5: GaMD SC region 0; 6: GaMD SC region 1

  std::unique_ptr< GpuBuffer<unsigned long long int>> pbPMEChargeGrid;      // Atom charge grid
  std::unique_ptr< GpuBuffer<PMEFloat>>   pbOrigAtomCharge;       // Atom charge storage
  std::unique_ptr< GpuBuffer<PMEFloat>>   pbFFTChargeGrid;        // Atom charge on FFT grid

  std::unique_ptr< GpuBuffer<unsigned long long int>> pbNumberTINBEntries;
  std::unique_ptr< GpuBuffer<int4>>       pbTINBList;

  std::unique_ptr< GpuBuffer<PMEDouble2>> pbTINb141;
  std::unique_ptr< GpuBuffer<PMEDouble3>> pbTINb142;
  std::unique_ptr< GpuBuffer<uint4>>      pbTINb14ID;

  std::unique_ptr< GpuBuffer<PMEFloat2>>  pbTILJTerm;

  std::unique_ptr< GpuBuffer<int>>        pbTIac;
  std::unique_ptr< GpuBuffer<int>>        pbTIico;
  std::unique_ptr< GpuBuffer<PMEFloat4>>  pbTIcn;
  std::unique_ptr< GpuBuffer<int>>        pbTIDcoef; // C4PairwiseCUDA
  std::unique_ptr< GpuBuffer<PMEFloat>>  pbTIDvalue; // C4PairwiseCUDA
  std::unique_ptr< GpuBuffer<PMEFloat4>>  pbTISigEps;
  std::unique_ptr< GpuBuffer<PMEFloat2>>  pbTISigMN;

  std::unique_ptr< GpuBuffer<PMEDouble2>> pbTIBond;         // Bond Kr, Req
  std::unique_ptr< GpuBuffer<uint2>>      pbTIBondID;       // Bond i, j
  std::unique_ptr< GpuBuffer<uint>>       pbTIBondType;     // Bond TI Type  
  std::unique_ptr< GpuBuffer<PMEDouble2>> pbTIBondAngle;    // Bond Angle Ka, Aeq
  std::unique_ptr< GpuBuffer<uint4>>      pbTIBondAngleID;  // Bond Angle i, j, k; regin
  std::unique_ptr< GpuBuffer<uint>>       pbTIBondAngleType;// Dihedral TI Type  
  std::unique_ptr< GpuBuffer<PMEDouble2>> pbTIDihedral1;    // Dihedral gmul, pn
  std::unique_ptr< GpuBuffer<PMEDouble2>> pbTIDihedral2;    // Dihedral pk, gamc
  std::unique_ptr< GpuBuffer<PMEDouble>>  pbTIDihedral3;    // Dihedral gams
  std::unique_ptr< GpuBuffer<uint4>>      pbTIDihedralID;   // Dihedral i, j, k, l
  std::unique_ptr< GpuBuffer<uint>>       pbTIDihedralType; // Dihedral TI Type  

  std::unique_ptr< GpuBuffer<PMEDouble4>> 		pbTINMRDistancesR;
  std::unique_ptr< GpuBuffer<PMEDouble2>> 		pbTINMRDistancesK;
  std::unique_ptr< GpuBuffer<uint4>>        pbTINMRDistancesID;

  std::unique_ptr< GpuBuffer<PMEDouble4>> 		pbTINMRAngleR;
  std::unique_ptr< GpuBuffer<PMEDouble2>> 		pbTINMRAngleK;
  std::unique_ptr< GpuBuffer<uint4>>        pbTINMRAngleID;

  std::unique_ptr< GpuBuffer<PMEDouble4>> 		pbTINMRDihedralR;
  std::unique_ptr< GpuBuffer<PMEDouble2>> 		pbTINMRDihedralK;
  std::unique_ptr< GpuBuffer<uint4>>        pbTINMRDihedralID;

  // READ Data section
  std::unique_ptr< GpuBuffer<unsigned>>   pbREAFList;
  std::unique_ptr< GpuBuffer<unsigned>>   pbREAFAtomList;

  std::unique_ptr< GpuBuffer<unsigned long long int>> pbNumberREAFNbEntries;

  std::unique_ptr< GpuBuffer<int4>>       pbREAFNbList;
  std::unique_ptr< GpuBuffer<PMEDouble2>> pbREAFNb141;
  std::unique_ptr< GpuBuffer<PMEDouble3>> pbREAFNb142;
  std::unique_ptr< GpuBuffer<uint4>>      pbREAFNb14ID;

  std::unique_ptr< GpuBuffer<PMEFloat2>>  pbREAFLJTerm;

  std::unique_ptr< GpuBuffer<int>>        pbREAFac;
  std::unique_ptr< GpuBuffer<int>>        pbREAFico;
  std::unique_ptr< GpuBuffer<PMEFloat4>>  pbREAFcn;
  std::unique_ptr< GpuBuffer<PMEFloat4>>  pbREAFSigEps;
  std::unique_ptr< GpuBuffer<PMEFloat2>>  pbREAFSigMN;

  std::unique_ptr< GpuBuffer<PMEDouble2>> pbREAFBondAngle;    // Bond Angle Ka, Aeq
  std::unique_ptr< GpuBuffer<uint4>>      pbREAFBondAngleID;  // Bond Angle i, j, k; regin
  std::unique_ptr< GpuBuffer<uint>>       pbREAFBondAngleType;// Dihedral TI Type  
  std::unique_ptr< GpuBuffer<PMEDouble2>> pbREAFDihedral1;    // Dihedral gmul, pn
  std::unique_ptr< GpuBuffer<PMEDouble2>> pbREAFDihedral2;    // Dihedral pk, gamc
  std::unique_ptr< GpuBuffer<PMEDouble>>  pbREAFDihedral3;    // Dihedral gams
  std::unique_ptr< GpuBuffer<uint4>>      pbREAFDihedralID;   // Dihedral i, j, k, l
  std::unique_ptr< GpuBuffer<uint>>       pbREAFDihedralType; // Dihedral TI Type  

  // MBAR
  std::unique_ptr< GpuBuffer<unsigned long long int>> pbMBAREnergy;  // energy values for states
  std::unique_ptr< GpuBuffer<PMEFloat>> pbMBARLambda;  // lambda values for states
  std::unique_ptr< GpuBuffer<PMEFloat>> pbMBARWeight;  // weight values for states 

  // LJ1264 Data section
  std::unique_ptr< GpuBuffer<unsigned long long int>> pbNumberLJ1264NBEntries;
  std::unique_ptr< GpuBuffer<unsigned>>   pbLJ1264AtomList;    
  std::unique_ptr< GpuBuffer<int4>>       pbLJ1264NBList;

  //RMSD fit
  std::unique_ptr< GpuBuffer<PMEFloat>> pbMatrix, pbResult, pbWork;
  std::unique_ptr< GpuBuffer<int>> pbInfo;

  // pLJ1264 Data section //C4PairwiseCUDA
  std::unique_ptr< GpuBuffer<unsigned long long int>> pbNumberpLJ1264NBEntries;
  std::unique_ptr< GpuBuffer<unsigned>>   pbpLJ1264AtomList;
  std::unique_ptr< GpuBuffer<int4>>       pbpLJ1264NBList;

  // LJ1264pLJ1264 Data section //C4PairwiseCUDA2023
  std::unique_ptr< GpuBuffer<unsigned long long int>> pbNumberLJ1264pLJ1264NBEntries;
  std::unique_ptr< GpuBuffer<unsigned>>   pbLJ1264pLJ1264AtomList;
  std::unique_ptr< GpuBuffer<int4>>       pbLJ1264pLJ1264NBList;

protected:

  void CleanUp();
};

#endif /*  _GTI_GPU_Context */
