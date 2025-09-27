#ifndef _GTI_simulationConst
#define _GTI_simulationConst
#ifdef GTI

#include "gti_def.h"
#include "gti_schedule_functions.h"
#include "base_simulationConst.h"

class gti_simulationConst : public base_simulationConst {

public:
 
  void InitSimulationConst();

// TI control variables for both C and FORTRAN are put in gti_controlVariable.i
// to enforce consistency
#define C_COMPILER
#include "../gti_controlVariable.i"
#undef  C_COMPILER

//data section

  //Flags/Types of various interactions
  static const int Fg_TI_region = 1;
  static const int Fg_excluded = 1 << 2;
  static const int Fg_has_SC = 1 << 3;
  static const int Fg_int_SC = 1 << 4;
  static const int Fg_int_TI = 1 << 6;

  static const int ex_need_LJ = 1 << 11;
  static const int ex_SC_LJ = 1 << 12;
  static const int ex_SC_ELE = 1 << 13;
  static const int ex_addToDVDL_ele = 1 << 14;
  static const int ex_addToDVDL_vdw = 1 << 15;
  static const int ex_addToDVDL_bat = 1 << 16;
  static const int ex_addToBat_corr = 1 << 17;
  static const int ex_userDefined_bat = 1 << 18;

  static const int Fg_has_RE = 1 << 21;
  static const int Fg_int_RE = 1 << 22;

  static const int ex_addRE_ele = 1 << 25;
  static const int ex_addRE_vdw = 1 << 26;
  static const int ex_addRE_bat = 1 << 27;

// Added and modified stuff
	bool needVirial;
	bool* pTLNeedNewNeighborList;

// TI stuff
	unsigned TIk;
  double TILambda[2], TIWeight[2], TIdWeight[2];
  PMEFloat TISCAlpha, TISCBeta, TIPrefactor;
  PMEFloat TIItemWeight[Schedule::TypeTotal][2];
  PMEFloat TIItemdWeight[Schedule::TypeTotal][2];
  PMEFloat TISCGamma, TISCGamma_factor[2], TISCGamma_dFactor[2];
  

  int addSC; //ctrlVar.addSC;
  int batSC; //ctrlVar.batSC;
  int eleGauss; //ctrlVar.eleGauss;
  int tiCut; //ctrlVar.tiCut;
  int eleSC; //ctrlVar.eleSC;
  int vdwSC; //ctrlVar.vdwSC;
  int autoAlpha; //ctrlVar.autoAlpha;
  int scaleBeta; //ctrlVar.scaleBeta;
  int eleExp; //ctrlVar.eleExp;
  int vdwExp; //ctrlVar.vdwExp;

  PMEFloat vdwCap; //ctrlVar.vdwCap;

  int cut_sc;// = ctrlVar.cut_sc;
  PMEFloat cut_sc0; //ctrlVar.cut_sc0;
  PMEFloat cut_sc1; //ctrlVar.cut_sc1;


  bool doTIHeating;
  PMEDouble TIHeatingTemp;

  int eleSmoothLambdaType;
  int vdwSmoothLambdaType;
  int SCSmoothLambdaType;


	unsigned *pTIList;	// Original TI list from pmemd
	unsigned *pSCList;	// Original TI list from pmemd
  unsigned* pSCBATList;	// Original TI list from pmemd

	unsigned numberTIAtoms;
	uint4 *pTIAtomList;  // .x: atom index
						// .y: TI region
						// .z: softcore or not

	unsigned numberUniqueTIAtoms; // numberUniqueTIAtoms =  number of TI common atoms + number of SC atoms
	unsigned *pUniqueTIAtomList;

	unsigned numberTICommonPairs;
	int4 *pTICommonPair;	// .x index of the common atom in the region #0		
							// .y index of the common atom in the region #1	
              // .z index of the SHAKE atom: 0: the region #0 is a SHAKE atom, 1: the region #1 is a SHAKE atom,; -1 none

	unsigned numberSoftcoreAtoms[2];	// number of softcore atoms for each region
	unsigned softcoreAtomShift;	// the shift between softcore atom list 
	unsigned *pTISoftcoreList;	// index of softcore atoms for each region (shifted by softcoreAtomStride) 

	unsigned long long int* pNumberTINBEntries;
	int4* pTINBList;		
    // .x: the first atom: always in TI region
		// .y: the second atom
    // .w: the LJ index in pTIcn; 
		// .z: the type of this non-bond interaction and the execution flag

	int4* pTINBCellList;

	PMEFloat* pOrigAtomCharge;

	unsigned long long int *pPMEChargeGrid;
	PMEFloat *pFFTChargeGrid;

  // Accumulators
  PMEAccumulator* pTIForce;
  PMEAccumulator* pTIForceX[3], *pTIForceY[3], *pTIForceZ[3];
  PMEAccumulator* pTISCForceX[2], *pTISCForceY[2], *pTISCForceZ[2];
  PMEAccumulator* pGaMDTIForceX[2], *pGaMDTIForceY[2], *pGaMDTIForceZ[2];

  PMEAccumulator* pTINBForce;
  PMEAccumulator* pTINBForceX[3], *pTINBForceY[3], *pTINBForceZ[3];
  PMEAccumulator* pTINBSCForceX[2], *pTINBSCForceY[2], *pTINBSCForceZ[2];
  PMEAccumulator* pGaMDTINBForceX[2], *pGaMDTINBForceY[2], *pGaMDTINBForceZ[2];

  unsigned long long int *pTIPotEnergyBuffer;
  unsigned long long int *pTILJ[3], *pTIEBond[3], *pTIEAngle[3], *pTIEDihedral[3], *pTIEER[3], *pTIEECC[3], *pTIEESC[3], *pTILJ14[3];
  unsigned long long int* pTIRestDist[3], * pTIRestAng[3], * pTIRestTor[3];
  unsigned long long int* pTIRestRMSD[3];
  unsigned long long int* pTIEE14CC[3], * pTIEE14SC[3], * pTIEE14CC_DL[3], * pTIEE14SC_DL[3];
  unsigned long long int *pTISCLJ[3], *pTISCBond[3], *pTISCAngle[3], *pTISCDihedral[3], *pTISCEED[3], *pTISCEE14[3], *pTISCLJ14[3];
  unsigned long long int *pTISCEBond_DL[3], *pTISCEAngle_DL[3], *pTISCEDihedral_DL[3];
  unsigned long long int *pTISCLJ_DL[3], *pTISCEECC_DL[3], *pTISCEESC_DL[3], *pTISCLJ14_DL[3];
  unsigned long long int *pTIRestDist_DL[3], *pTIRestAng_DL[3], *pTIRestTor_DL[3], * pTIRestRMSD_DL[3];
  unsigned long long int* pTISpecialTerms[TINumberSpecialTerms][3];
  unsigned long long int *pV11[3], *pV22[3], *pV33[3];
  
	unsigned long long int *pTIKinEnergyBuffer;


	int            numberTI14NBEntries;
  PMEDouble2* pTINb141;
  PMEDouble3 *pTINb142;
	uint4*          pTINb14ID;

	unsigned int  TIVdwNTyes;
        unsigned int  TIC4Pairwise; //C4PairwiseCUDA2023
	int*  pTIac;
	int*  pTIico;
	PMEFloat4*     pTIcn;
        int*  pTIDcoef;             //C4PairwiseCUDA
        PMEFloat*     pTIDvalue;   //C4PairwiseCUDA
	PMEFloat4*     pTISigEps;
        PMEFloat2*     pTISigMN;

	unsigned int	numberTIBond;
	PMEDouble2*		pTIBond;                              // Bond rk, req
	uint2*        pTIBondID;                            // Bond i, j
  uint*         pTIBondType;

	unsigned int	numberTIAngle;
	PMEDouble2*   pTIBondAngle;                         // Bond Angle Kt, teq
	uint4*        pTIBondAngleID;                      // Bond Angle i, j
  uint*         pTIBondAngleType;

	unsigned int	numberTIDihedral;
	PMEDouble2*   pTIDihedral1;                         // Dihedral Ipn, pn
	PMEDouble2*   pTIDihedral2;                         // Dihedral pk, gamc
	PMEDouble*    pTIDihedral3;                         // Dihedral gams
	uint4*        pTIDihedralID;                       // Dihedral i, j, k, l
  uint*         pTIDihedralType;

  // NMR restraint terms
  unsigned int numberTINMRDistance;
  PMEDouble4*		pTINMRDistancesR;    
  PMEDouble2*		pTINMRDistancesK;     
  uint4*        pTINMRDistancesID;                          

  unsigned int numberTINMRAngle;
  PMEDouble4*		pTINMRAngleR;
  PMEDouble2*		pTINMRAngleK;
  uint4*        pTINMRAngleID;

  unsigned int numberTINMRDihedral;
  PMEDouble4*		pTINMRDihedralR;
  PMEDouble2*		pTINMRDihedralK;
  uint4*        pTINMRDihedralID;
  
//  MBAR
	bool needMBAR;   
	unsigned nMBARStates;  // Number of MBAR states 	
  int currentMBARState;
  int rangeMBAR;

	unsigned long long int* pMBAREnergy;	// Energies for states
  PMEFloat* pMBARLambda; // lambdas for MBAR
  PMEFloat* pMBARWeight; // lambdas for MBAR

// REAF section
  int reafMode;
  int addRE;
  int eleTauType;
  int vdwTauType;
  PMEFloat current_tau;

  PMEFloat REAFWeight[2], REAFItemWeight[Schedule::TypeTotal][2];

  unsigned numberREAFAtoms;
  unsigned* pREAFList;
  unsigned* pREAFAtomList;

  unsigned long long int* pNumberREAFNbEntries;
  int4* pREAFNbList;

  int         numberREAF14NBEntries;
  PMEDouble2* pREAFNb141;
  PMEDouble3* pREAFNb142;
  uint4* pREAFNb14ID;

  unsigned int	numberREAFAngle;
  PMEDouble2* pREAFBondAngle;                         // Bond Angle Kt, teq
  uint4* pREAFBondAngleID;                      // Bond Angle i, j
  uint* pREAFBondAngleType;

  unsigned int	numberREAFDihedral;
  PMEDouble2* pREAFDihedral1;                         // Dihedral Ipn, pn
  PMEDouble2* pREAFDihedral2;                         // Dihedral pk, gamc
  PMEDouble* pREAFDihedral3;                         // Dihedral gams
  uint4* pREAFDihedralID;                       // Dihedral i, j, k, l
  uint* pREAFDihedralType;

  // LJ-12-6-4 section
	unsigned numberLJ1264Atoms; 
	unsigned *pLJ1264AtomList;
	unsigned long long int *pNumberLJ1264NBEntries;
	int4* pLJ1264NBList;		// .x: the metal ion
							// .y: the LJ pair
              
					// .y: the LJ pair
  // pLJ-12-6-4 section // C4PairwiseCUDA
        unsigned numberpLJ1264Atoms;
        unsigned *ppLJ1264AtomList;
        unsigned long long int *pNumberpLJ1264NBEntries;
        int4* ppLJ1264NBList;            // .x: the metal ion
                                         // .y: the LJ pair

  // LJ-12-6-4 and pLJ-12-6-4 coexist section // C4PairwiseCUDA2023
        unsigned numberLJ1264pLJ1264Atoms;
        unsigned *pLJ1264pLJ1264AtomList;
        unsigned long long int *pNumberLJ1264pLJ1264NBEntries;
        int4* pLJ1264pLJ1264NBList;  

  unsigned rmsd_type;
  unsigned number_rmsd_set;
  unsigned rmsd_atom_count[5];
  unsigned rmsd_atom_max_count;
  unsigned rmsd_ti_region[MaxNumberRMSDRegion];
  unsigned rmsd_atom_list[MaxNumberRMSDRegion * MaxNumberRMSDAtom];
  PMEFloat rmsd_weights[MaxNumberRMSDRegion];
  PMEFloat rmsd_ref_crd[MaxNumberRMSDRegion * MaxNumberRMSDAtom * 3];
  PMEFloat rmsd_ref_com[MaxNumberRMSDRegion * 3];

  PMEFloat *pMatrix, *pResult, *pWork;
  int* pInfo;
 
};


#endif
#endif
