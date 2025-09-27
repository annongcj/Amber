#ifdef GTI

#include <vector>
#include <algorithm>
#include <stdexcept> 
#include "gpu.h"
#include "gti_gpu.h"
#include "gti_cuda.cuh"
#include "gti_schedule_functions.h"
#include "gti_f95.h"

//---------------------------------------------------------------------------------------------
// gti_sync_vector_:
//
// Arguments:
//   mode:          
//   combinedMode:  
//---------------------------------------------------------------------------------------------
extern "C" void gti_sync_vector_(int* mode, int* combinedMode)
{
  PRINTMETHOD("gti_sync_vector");

  // mode: 0 for force; 1 for velocity; 2 for coord
  gpuContext gpu = theGPUContext::GetPointer();
  ik_SyncVector_kernel(gpu, *mode, *combinedMode);
}

//---------------------------------------------------------------------------------------------
// gti_init_md_parameters_:
//
// Arguments:
//   vlimit:
//---------------------------------------------------------------------------------------------
extern "C" void gti_init_md_parameters_(double* vlimit)
{
  PRINTMETHOD("gti_init_md_parameters");

  gpuContext gpu = theGPUContext::GetPointer();
  gpu->InitMDParameters(*vlimit);

}

//---------------------------------------------------------------------------------------------
// gti_init_ti_parameters_:
//
// Arguments:
//   ti_latm_lst:  
//   ti_lst:       
//   ti_sc_lst:    
//   plambda:      
//   pklambda:     
//   pscalpha:     
//   pscbeta:      
//---------------------------------------------------------------------------------------------
extern "C" void gti_init_ti_parameters_(int ti_latm_lst[][3], int ti_lst[][3], int ti_sc_lst[], int ti_sc_bat_lst[], gti_simulationConst::TypeCtrlVar* pCtrlVar)
{
  PRINTMETHOD(__func__);

  gpuContext gpu = theGPUContext::GetPointer();
  gpu->InitAccumulators();
  gpu->InitTIParameters(ti_latm_lst, ti_lst, ti_sc_lst, ti_sc_bat_lst, *pCtrlVar);
}

extern "C" void gti_init_reaf_parameters_(int* reaf_mode, double* reaf_tau, int* addRE, int reaf_atom_list[][2]) {
  PRINTMETHOD(__func__);

  gpuContext gpu = theGPUContext::GetPointer();
  gpu->InitREADParameters(*reaf_mode, *reaf_tau, *addRE, reaf_atom_list);
}

//---------------------------------------------------------------------------------------------
// gti_init_gamd_ti_parameters_:
//
// Arguments:
//   ti_latm_lst:  
//   ti_lst:       
//   ti_sc_lst:    
//   plambda:      
//   pklambda:     
//   pscalpha:     
//   pscbeta:      
//---------------------------------------------------------------------------------------------
extern "C" void gti_init_gamd_ti_parameters_(int ti_latm_lst[][3], int ti_lst[][3], int ti_sc_lst[], int ti_sc_bat_lst[],
                                        double* plambda, unsigned* pklambda, double* pscalpha,
                                        double* pscbeta, double* pscgamma, int* pCut,
                                        int* pAddSC, int* pEleGauss, int*  pEleSC, int* pVdwSC, double* pVdwCap)
{
  PRINTMETHOD("gti_init_gamd_ti_parameters");

  gpuContext gpu = theGPUContext::GetPointer();
  gpu->InitAccumulators(); 
  gpu->InitGaMDTIParameters(ti_latm_lst, ti_lst, ti_sc_lst, ti_sc_bat_lst , *plambda, *pklambda, *pscalpha,
    *pscbeta, *pscgamma, *pCut, *pAddSC, *pEleGauss, *pEleSC, *pVdwSC, *pVdwCap);
}

//---------------------------------------------------------------------------------------------
// gti_init_mbar_parameters_:
//
// Arguments:
//   pMBAR_states:  
//   MBAR_lambda:   
//---------------------------------------------------------------------------------------------
extern "C" void gti_init_mbar_parameters_(int* pMBAR_states, double MBAR_lambda[][2], int* pRangeMBAR)
{
  PRINTMETHOD("gti_init_mbar_parameters");

  gpuContext gpu = theGPUContext::GetPointer();
  if (*pMBAR_states <= 0) {
    return;
  }
  gpu->InitMBARParameters((unsigned)*pMBAR_states, MBAR_lambda, *pRangeMBAR);
}

//---------------------------------------------------------------------------------------------
// gti_setup_update_:
//
// Arguments:
//   charge:     
//   sync_mass:  
//   mass:       
//---------------------------------------------------------------------------------------------
extern "C" void gti_setup_update_(double charge[], int* sync_mass, double mass[])
{
  PRINTMETHOD("gti_setup_update");

  gpuContext gpu = theGPUContext::GetPointer();
  gpu->UpdateCharges(charge);

  if (gpu->sim.numberTIAtoms>0 && *sync_mass>0) {
    gpu->UpdateTIMasses(mass);
  }
}

//---------------------------------------------------------------------------------------------
// gti_nb_setup_:
//
// Arguments:
//   pNtypes:  
//   iac:      
//   ico:      
//   cn1:      
//   cn2:      
//   cn6:      
//   cn7:
//   cn8: C4PairwiseCUDA
//---------------------------------------------------------------------------------------------
extern "C" void gti_nb_setup_(int* pNtypes, int iac[], int ico[], double cn1[], double cn2[],
                              double cn6[], int cn7[], double cn8[], int* pC4Pairwise)
{
  PRINTMETHOD("gti_nb_setup");

  gpuContext gpu = theGPUContext::GetPointer();

  unsigned int ntypes = *pNtypes;
  unsigned int C4Pairwise = *pC4Pairwise;
  gpu->sim.TIVdwNTyes = ntypes;
  gpu->sim.TIC4Pairwise = C4Pairwise; //C4PairwiseCUDA2023
  gpu->pbTIac = std::unique_ptr< GpuBuffer<int> > (new GpuBuffer<int>(gpu->sim.atoms));
  gpu->pbTIico = std::unique_ptr< GpuBuffer<int> > (new GpuBuffer<int>(ntypes*ntypes));
  gpu->pbTIcn = std::unique_ptr< GpuBuffer<PMEFloat4> > (new GpuBuffer<PMEFloat4>(ntypes*(ntypes + 1) / 2));
  gpu->pbTIDcoef = std::unique_ptr< GpuBuffer<int> > (new GpuBuffer<int>(C4Pairwise*2)); // C4PairwiseCUDA
  gpu->pbTIDvalue = std::unique_ptr< GpuBuffer<PMEFloat> > (new GpuBuffer<PMEFloat>(C4Pairwise)); // C4PairwiseCUDA
  gpu->pbTISigEps = std::unique_ptr< GpuBuffer<PMEFloat4> > (new GpuBuffer<PMEFloat4>(ntypes*(ntypes + 1) / 2));
  gpu->pbTISigMN = std::unique_ptr< GpuBuffer<PMEFloat2> >(new GpuBuffer<PMEFloat2>(ntypes * (ntypes + 1) / 2));

  for (unsigned i = 0; i < C4Pairwise; i++) {
    gpu->pbTIDcoef->_pSysData[i*2] = (cn7[i*3])/3;
    gpu->pbTIDcoef->_pSysData[i*2+1] = (cn7[i*3+1])/3;
  } // C4PairwiseCUDA
  
  for (unsigned i = 0; i < C4Pairwise; i++) {
    gpu->pbTIDvalue->_pSysData[i] = (PMEFloat)cn8[i];
  } // C4PairwiseCUDA

  for (unsigned i = 0; i < gpu->sim.atoms; i++) {
    gpu->pbTIac->_pSysData[i] = iac[i];
  }
  for (unsigned i = 0; i < ntypes*ntypes; i++) {
    gpu->pbTIico->_pSysData[i] = ico[i];
  }

  for (unsigned i = 0; i < ntypes*(ntypes + 1) / 2; i++) {
    gpu->pbTIcn->_pSysData[i].x = (PMEFloat)cn1[i];
    gpu->pbTIcn->_pSysData[i].y = (PMEFloat)cn2[i];
    gpu->pbTIcn->_pSysData[i].z = (PMEFloat)cn6[i];

    gpu->pbTISigEps->_pSysData[i].x = 0.0;
    gpu->pbTISigEps->_pSysData[i].y = 0.0;
    gpu->pbTISigEps->_pSysData[i].z = 0;
    gpu->pbTISigEps->_pSysData[i].w = -1.0;

    if (fabs(cn1[i]) > 1e-5 && fabs(cn2[i]) > 1e-5) {
      double sigma6 = cn1[i] / cn2[i];
      gpu->pbTISigEps->_pSysData[i].x = (PMEFloat)(1.0/sigma6); // signma^-6
      gpu->pbTISigEps->_pSysData[i].y = (PMEFloat)(cn2[i] * cn2[i] / cn1[i]); // 4*epslon

      gpu->pbTISigMN->_pSysData[i].x = (PMEFloat)(pow(sigma6, gpu->sim.eleExp/6.0));
      gpu->pbTISigMN->_pSysData[i].y = (PMEFloat)(pow(sigma6, gpu->sim.vdwExp/6.0));

      if (fabs(cn6[i]) > 1e-5) {
        double ss = pow(cn2[i] / cn1[i], 1.0 / 3.0);
        gpu->pbTISigEps->_pSysData[i].z = (PMEFloat)(ss*ss*cn6[i]);
        gpu->pbTISigEps->_pSysData[i].w = 1.0;
      }
    }
  }

  gpu->pbTIac->Upload();
  gpu->pbTIico->Upload();
  gpu->pbTIcn->Upload();
  gpu->pbTIDcoef->Upload(); // C4PairwiseCUDA
  gpu->pbTIDvalue->Upload(); // C4PairwiseCUDA
  gpu->pbTISigEps->Upload();
  gpu->pbTISigMN->Upload();
  gpu->sim.pTIac = gpu->pbTIac->_pDevData;
  gpu->sim.pTIico = gpu->pbTIico->_pDevData;
  gpu->sim.pTIcn = gpu->pbTIcn->_pDevData;
  gpu->sim.pTIDcoef = gpu->pbTIDcoef->_pDevData; // C4PairwiseCUDA
  gpu->sim.pTIDvalue = gpu->pbTIDvalue->_pDevData; // C4PairwiseCUDA
  gpu->sim.pTISigEps = gpu->pbTISigEps->_pDevData;
  gpu->sim.pTISigMN = gpu->pbTISigMN->_pDevData;

  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gti_ti_nb_setup_:
//
// Arguments:       
//   pNtypes:       
//   iac:           
//   ico:           
//   pNb14_cnt:     
//   cit_nb14:      
//   gbl_one_scee:  
//   gbl_one_scnb:  
//   cn114:         
//   cn214:         
//---------------------------------------------------------------------------------------------
extern "C" void gti_ti_nb_setup_(int* pNtypes, int iac[], int ico[], int* pNb14_cnt,
                                 int cit_nb14[][3], double gbl_one_scee[],
                                 double gbl_one_scnb[], double cn114[], double cn214[])
{
  PRINTMETHOD(__func__);

  gpuContext gpu = theGPUContext::GetPointer();
  unsigned int ntypes = *pNtypes;
  unsigned maxNumberTINBEntries = gpu->sim.numberTIAtoms * gti_simulationConst::MaxNumberNBPerAtom;
  gpu->pbTINBList = std::unique_ptr< GpuBuffer<int4> > (new GpuBuffer<int4>(maxNumberTINBEntries));

  for (unsigned i = 0; i < maxNumberTINBEntries; i++) {
    gpu->pbTINBList->_pSysData[i] = { -1, -1, -1, -1 };
  }
  gpu->pbTINBList->Upload();
  gpu->sim.pTINBList = gpu->pbTINBList->_pDevData;
  gpu->pbNumberTINBEntries = std::unique_ptr< GpuBuffer<unsigned long long int> > (new GpuBuffer<unsigned long long int>(1));
  gpu->pbNumberTINBEntries->_pSysData[0] = maxNumberTINBEntries;
  gpu->pbNumberTINBEntries->Upload();
  gpu->sim.pNumberTINBEntries = gpu->pbNumberTINBEntries->_pDevData;

  // 1:4 NB part
  unsigned int nb14s = (gpu->ntf < 8) ? *pNb14_cnt : 0;
  int(*temp)[4] = new int[nb14s][4];
  bool inRegion[2][2];
  unsigned counter = 0;
  for (unsigned int i = 0; i < nb14s; i++) {
    for (unsigned int j = 0; j < 2; j++) {
      unsigned atom = cit_nb14[i][j] - 1;
      inRegion[0][j] = (gpu->pbTIList->_pSysData[atom])>0;
      inRegion[1][j] = (gpu->pbTIList->_pSysData[atom + gpu->sim.stride])>0;
    }
    temp[counter][3] = -1;
    if (inRegion[0][0] || inRegion[0][1]) {
      temp[counter][3] = 0;
    }
    if (inRegion[1][0] || inRegion[1][1]) {
      temp[counter][3] = 1;
    }
    if (temp[counter][3] >= 0) {
      for (unsigned int j = 0; j < 3; j++) {
        temp[counter][j] = cit_nb14[i][j];
      }
      cit_nb14[i][2] = -1;
      counter++;
    }
  }
  nb14s = counter;

  // Copy 1-4 interactions
  gpu->pbTINb141 = std::unique_ptr< GpuBuffer<PMEDouble2> > (new GpuBuffer<PMEDouble2>(nb14s));
  gpu->pbTINb142 = std::unique_ptr< GpuBuffer<PMEDouble3> > (new GpuBuffer<PMEDouble3>(nb14s));
  gpu->pbTINb14ID = std::unique_ptr< GpuBuffer<uint4> > (new GpuBuffer<uint4>(nb14s));
  if (gpu->ntf < 8) {
    for (unsigned int i = 0; i < nb14s; i++) {
      int parm_idx = temp[i][2] - 1;
      unsigned atom0 = abs(temp[i][0]) - 1;
      unsigned atom1 = abs(temp[i][1]) - 1;

      gpu->pbTINb141->_pSysData[i].x = gbl_one_scee[parm_idx];
      gpu->pbTINb141->_pSysData[i].y = gbl_one_scnb[parm_idx];
      int tt = ntypes * (iac[abs(temp[i][0]) - 1] - 1) + (iac[abs(temp[i][1]) - 1] - 1);
      int nbtype = tt >= 0 ? ico[tt] - 1 : -1;
      if (nbtype >= 0) {
        gpu->pbTINb142->_pSysData[i].x = cn114[nbtype];
        gpu->pbTINb142->_pSysData[i].y = cn214[nbtype];
        gpu->pbTINb142->_pSysData[i].z = cn214[nbtype]/cn114[nbtype]; // B/A
      }
      else {
        gpu->pbTINb142->_pSysData[i].x = ZeroF;
        gpu->pbTINb142->_pSysData[i].y = ZeroF;
        gpu->pbTINb142->_pSysData[i].z = ZeroF;
      }
      gpu->pbTINb14ID->_pSysData[i].x = atom0;
      gpu->pbTINb14ID->_pSysData[i].y = atom1;
      gpu->pbTINb14ID->_pSysData[i].z = temp[i][3];
    }
  }
  delete[] temp;
  gpu->pbTINb141->Upload();
  gpu->pbTINb142->Upload();
  gpu->pbTINb14ID->Upload();

  // Set constants
  gpu->sim.numberTI14NBEntries = nb14s;
  gpu->sim.pTINb141 = gpu->pbTINb141->_pDevData;
  gpu->sim.pTINb142 = gpu->pbTINb142->_pDevData;
  gpu->sim.pTINb14ID = gpu->pbTINb14ID->_pDevData;

  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gti_reaf_nb_setup_:
//
// Arguments:
//   atm_isymbl:  
//   ntypes:      
//   atm_iac:     
//   ico:         
//   gbl_cn6:     
//---------------------------------------------------------------------------------------------
extern "C" void gti_reaf_nb_setup_(int* pNtypes, int iac[], int ico[], int* pNb14_cnt,
  int cit_nb14[][3], double gbl_one_scee[],
  double gbl_one_scnb[], double cn114[], double cn214[])
{
  PRINTMETHOD(__func__);

  if (*pNtypes <= 0 ) return;

  gpuContext gpu = theGPUContext::GetPointer();

  if (gpu->sim.reafMode < 0) return;
  unsigned nlistatoms = gpu->sim.numberREAFAtoms;
  if (nlistatoms == 0) return;

  unsigned ntypes = *pNtypes;

  // regular NB
    unsigned maxNumberREAFNbEntries = gpu->sim.numberREAFAtoms * gti_simulationConst::MaxNumberNBPerAtom;
    gpu->pbREAFNbList = std::unique_ptr< GpuBuffer<int4> >(new GpuBuffer<int4>(maxNumberREAFNbEntries));
    for (unsigned i = 0; i < maxNumberREAFNbEntries; i++) {
      gpu->pbREAFNbList->_pSysData[i] = { -1, -1, -1, -1 };
    }
    gpu->pbREAFNbList->Upload();
    gpu->sim.pREAFNbList = gpu->pbREAFNbList->_pDevData;
    gpu->pbNumberREAFNbEntries = std::unique_ptr< GpuBuffer<unsigned long long int> >(new GpuBuffer<unsigned long long int>(1));
    gpu->pbNumberREAFNbEntries->_pSysData[0] = maxNumberREAFNbEntries;
    gpu->pbNumberREAFNbEntries->Upload();
    gpu->sim.pNumberREAFNbEntries = gpu->pbNumberREAFNbEntries->_pDevData;

  if (*pNb14_cnt <= 0) return;
  // 1:4 NB part
  unsigned int nb14s = (gpu->ntf < 8) ? *pNb14_cnt : 0;
  int(*temp)[4] = new int[nb14s][4];

  unsigned counter = 0;
  for (unsigned int i = 0; i < nb14s; i++) {
    bool inRegion = false;
    temp[counter][3] = -1;

    for (unsigned int j = 0; j < 2; j++) {
      unsigned atom = cit_nb14[i][j] - 1;
      if (gpu->pbREAFList->_pSysData[atom] > 0) {
        inRegion = true;
        break;
      }
    }

    if (inRegion) {
      temp[counter][3] = 0;
      for (unsigned int j = 0; j < 3; j++) {
        temp[counter][j] = cit_nb14[i][j];
      }
      counter++;
    }
  }
  nb14s = counter;

  // Copy 1-4 interactions
  gpu->pbREAFNb141 = std::unique_ptr< GpuBuffer<PMEDouble2> >(new GpuBuffer<PMEDouble2>(nb14s));
  gpu->pbREAFNb142 = std::unique_ptr< GpuBuffer<PMEDouble3> >(new GpuBuffer<PMEDouble3>(nb14s));
  gpu->pbREAFNb14ID = std::unique_ptr< GpuBuffer<uint4> >(new GpuBuffer<uint4>(nb14s));
  if (gpu->ntf < 8) {
    for (unsigned int i = 0; i < nb14s; i++) {
      int parm_idx = temp[i][2] - 1;
      unsigned atom0 = abs(temp[i][0]) - 1;
      unsigned atom1 = abs(temp[i][1]) - 1;

      gpu->pbREAFNb141->_pSysData[i].x = gbl_one_scee[parm_idx];
      gpu->pbREAFNb141->_pSysData[i].y = gbl_one_scnb[parm_idx];
      int tt = ntypes * (iac[abs(temp[i][0]) - 1] - 1) + (iac[abs(temp[i][1]) - 1] - 1);
      int nbtype = tt >= 0 ? ico[tt] - 1 : -1;
      if (nbtype >= 0) {
        gpu->pbREAFNb142->_pSysData[i].x = cn114[nbtype];
        gpu->pbREAFNb142->_pSysData[i].y = cn214[nbtype];
        gpu->pbREAFNb142->_pSysData[i].z = cn214[nbtype] / cn114[nbtype]; // B/A
      }
      else {
        gpu->pbREAFNb142->_pSysData[i].x = ZeroF;
        gpu->pbREAFNb142->_pSysData[i].y = ZeroF;
        gpu->pbREAFNb142->_pSysData[i].z = ZeroF;
      }
      gpu->pbREAFNb14ID->_pSysData[i].x = atom0;
      gpu->pbREAFNb14ID->_pSysData[i].y = atom1;
      gpu->pbREAFNb14ID->_pSysData[i].z = temp[i][3];
    }
  }
  delete[] temp;
  gpu->pbREAFNb141->Upload();
  gpu->pbREAFNb142->Upload();
  gpu->pbREAFNb14ID->Upload();

  // Set constants
  gpu->sim.numberREAF14NBEntries = nb14s;
  gpu->sim.pREAFNb141 = gpu->pbREAFNb141->_pDevData;
  gpu->sim.pREAFNb142 = gpu->pbREAFNb142->_pDevData;
  gpu->sim.pREAFNb14ID = gpu->pbREAFNb14ID->_pDevData;

  gpuCopyConstants();
}


//---------------------------------------------------------------------------------------------
// gti_ti_nb_gamd_setup_:
//
// Arguments:       
//   pNtypes:       
//   iac:           
//   ico:           
//   pNb14_cnt:     
//   cit_nb14:      
//   gbl_one_scee:  
//   gbl_one_scnb:  
//   cn114:         
//   cn214:         
//---------------------------------------------------------------------------------------------
extern "C" void gti_ti_nb_gamd_setup_(int* pNtypes, int iac[], int ico[], int* pNb14_cnt,
                                 int cit_nb14[][3], double gbl_one_scee[],
                                 double gbl_one_scnb[], double cn114[], double cn214[])
{
  PRINTMETHOD("gti_ti_nb_gamd_setup");

  gpuContext gpu = theGPUContext::GetPointer();
  unsigned int ntypes = *pNtypes;
  unsigned maxNumberTINBEntries = gpu->sim.numberTIAtoms * 3000;
  gpu->pbTINBList = std::unique_ptr< GpuBuffer<int4> > (new GpuBuffer<int4>(maxNumberTINBEntries));

  for (unsigned i = 0; i < maxNumberTINBEntries; i++) {
    gpu->pbTINBList->_pSysData[i].x = -1;
    gpu->pbTINBList->_pSysData[i].y = -1;
    gpu->pbTINBList->_pSysData[i].z = -1;
    gpu->pbTINBList->_pSysData[i].w = -1;
  }
  gpu->pbTINBList->Upload();
  gpu->sim.pTINBList = gpu->pbTINBList->_pDevData;
  gpu->pbNumberTINBEntries = std::unique_ptr< GpuBuffer<unsigned long long int> > (new GpuBuffer<unsigned long long int>(1));
  gpu->pbNumberTINBEntries->_pSysData[0] = 0;
  gpu->pbNumberTINBEntries->Upload();
  gpu->sim.pNumberTINBEntries = gpu->pbNumberTINBEntries->_pDevData;

  // 1:4 NB part
  unsigned int nb14s = (gpu->ntf < 8) ? *pNb14_cnt : 0;
  int(*temp)[4] = new int[nb14s][4];
  bool inRegion[2][2];
  unsigned counter = 0;
  for (unsigned int i = 0; i < nb14s; i++) {
    for (unsigned int j = 0; j < 2; j++) {
      unsigned atom = cit_nb14[i][j] - 1;
      inRegion[0][j] = (gpu->pbTIList->_pSysData[atom])>0;
      inRegion[1][j] = (gpu->pbTIList->_pSysData[atom + gpu->sim.stride])>0;
    }
    temp[counter][3] = -1;
    if (inRegion[0][0] || inRegion[0][1]) {
      temp[counter][3] = 0;
    }
    if (inRegion[1][0] || inRegion[1][1]) {
      temp[counter][3] = 1;
    }
    if (temp[counter][3] >= 0) {
      for (unsigned int j = 0; j < 3; j++) {
        temp[counter][j] = cit_nb14[i][j];
      }
//      cit_nb14[i][2]=-1;
      counter++;
    }
  }
  nb14s = counter;

  // Copy 1-4 interactions
  gpu->pbTINb141 = std::unique_ptr< GpuBuffer<PMEDouble2> > (new GpuBuffer<PMEDouble2>(nb14s));
  gpu->pbTINb142 = std::unique_ptr< GpuBuffer<PMEDouble3> > (new GpuBuffer<PMEDouble3>(nb14s));
  gpu->pbTINb14ID = std::unique_ptr< GpuBuffer<uint4> > (new GpuBuffer<uint4>(nb14s));
  if (gpu->ntf < 8) {
    for (unsigned int i = 0; i < nb14s; i++) {
      int parm_idx = temp[i][2] - 1;
      unsigned atom0 = abs(temp[i][0]) - 1;
      unsigned atom1 = abs(temp[i][1]) - 1;

      gpu->pbTINb141->_pSysData[i].x = gbl_one_scee[parm_idx];
      gpu->pbTINb141->_pSysData[i].y = gbl_one_scnb[parm_idx];
      int tt = ntypes * (iac[abs(temp[i][0]) - 1] - 1) + (iac[abs(temp[i][1]) - 1] - 1);
      int nbtype = tt >= 0 ? ico[tt] - 1 : -1;
      if (nbtype >= 0) {
        gpu->pbTINb142->_pSysData[i].x = cn114[nbtype];
        gpu->pbTINb142->_pSysData[i].y = cn214[nbtype];
        gpu->pbTINb142->_pSysData[i].z = cn214[nbtype]/cn114[nbtype]; // B/A
      }
      else {
        gpu->pbTINb142->_pSysData[i].x = ZeroF;
        gpu->pbTINb142->_pSysData[i].y = ZeroF;
        gpu->pbTINb142->_pSysData[i].z = ZeroF;
      }
      gpu->pbTINb14ID->_pSysData[i].x = atom0;
      gpu->pbTINb14ID->_pSysData[i].y = atom1;
      gpu->pbTINb14ID->_pSysData[i].z = temp[i][3];
    }
  }
  delete[] temp;
  gpu->pbTINb141->Upload();
  gpu->pbTINb142->Upload();
  gpu->pbTINb14ID->Upload();

  // Set constants
  gpu->sim.numberTI14NBEntries = nb14s;
  gpu->sim.pTINb141 = gpu->pbTINb141->_pDevData;
  gpu->sim.pTINb142 = gpu->pbTINb142->_pDevData;
  gpu->sim.pTINb14ID = gpu->pbTINb14ID->_pDevData;
}

//---------------------------------------------------------------------------------------------
// gti_lj1264_nb_setup_:
//
// Arguments:
//   atm_isymbl:  
//   ntypes:      
//   atm_iac:     
//   ico:         
//   gbl_cn6:     
//---------------------------------------------------------------------------------------------
extern "C" void gti_lj1264_nb_setup_(char atm_isymbl[][4], int* ntypes, int atm_iac[],
  int ico[], double gbl_cn6[])
{
  PRINTMETHOD(__func__);

  gpuContext gpu = theGPUContext::GetPointer();
  unsigned natoms = gpu->sim.atoms;
  unsigned nlistatoms = 0;
  unsigned* list = new unsigned[natoms];

  for (unsigned i = 0; i < natoms; i++) {
    if (atm_iac[i] > 0) {
      for (unsigned j = 1; j < 4; j++) {
        char tt = atm_isymbl[i][j];
        if (tt == '+' || tt == '-') {
          int index = ico[(*ntypes) * (atm_iac[i] - 1) + atm_iac[i] - 1] - 1;
          if (index >= 0) {
            if (gbl_cn6[index] > 1e-10) {
              list[nlistatoms] = i;
              nlistatoms++;
              break;
            }
          }
        }
      }
    }
  }
  gpu->sim.numberLJ1264Atoms = nlistatoms;
  if (nlistatoms > 0) {
    gpu->pbLJ1264AtomList = std::unique_ptr< GpuBuffer<unsigned> >(new GpuBuffer<unsigned>(nlistatoms));
    for (unsigned i = 0; i < nlistatoms; i++) {
      gpu->pbLJ1264AtomList->_pSysData[i] = list[i];
    }
    gpu->pbLJ1264AtomList->Upload();
    gpu->sim.pLJ1264AtomList = gpu->pbLJ1264AtomList->_pDevData;
    unsigned maxNumberLJ1264NBEntries = gpu->sim.numberLJ1264Atoms * gti_simulationConst::MaxNumberNBPerAtom;
    gpu->pbLJ1264NBList = std::unique_ptr< GpuBuffer<int4> >(new GpuBuffer<int4>(maxNumberLJ1264NBEntries));
    for (unsigned i = 0; i < maxNumberLJ1264NBEntries; i++) {
      gpu->pbLJ1264NBList->_pSysData[i] = { -1, -1, -1, -1 };
    }
    gpu->pbLJ1264NBList->Upload();
    gpu->sim.pLJ1264NBList = gpu->pbLJ1264NBList->_pDevData;
    gpu->pbNumberLJ1264NBEntries = std::unique_ptr< GpuBuffer<unsigned long long int> >(new GpuBuffer<unsigned long long int>(1));
    gpu->pbNumberLJ1264NBEntries->_pSysData[0] = maxNumberLJ1264NBEntries;
    gpu->pbNumberLJ1264NBEntries->Upload();
    gpu->sim.pNumberLJ1264NBEntries = gpu->pbNumberLJ1264NBEntries->_pDevData;
  }

  delete[] list;
  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gti_plj1264_nb_setup_:
//
// Arguments:
//   atm_isymbl:
//   ntypes:
//   atm_iac:
//   ico:
//   gbl_cn7: // C4PairwiseCUDA
//   gbl_cn8:
//   C4Pairwise: 
//---------------------------------------------------------------------------------------------
extern "C" void gti_plj1264_nb_setup_(char atm_isymbl[][4], int* ntypes, int atm_iac[],
  int ico[], int gbl_cn7[], double gbl_cn8[], int* pC4Pairwise)
{
  PRINTMETHOD(__func__);

  gpuContext gpu = theGPUContext::GetPointer();
  unsigned natoms = gpu->sim.atoms;
  unsigned nlistatoms = *pC4Pairwise*2; //C4PairwiseCUDA
  unsigned* list = new unsigned[*pC4Pairwise*2];
  for (unsigned k=0;k<*pC4Pairwise;k++) { //C4PairwiseCUDA2023
    //if (i==(gbl_cn7[k*3]+3)/3 || i==(gbl_cn7[k*3+1]+3)/3) { //if any pairwise matches, check whether already counted
      list[k*2] = (gbl_cn7[k*3])/3;
      list[k*2+1]=(gbl_cn7[k*3+1])/3;
    //}
  }

  gpu->sim.numberpLJ1264Atoms = nlistatoms;
  if (nlistatoms > 0) {
    gpu->pbpLJ1264AtomList = std::unique_ptr< GpuBuffer<unsigned> >(new GpuBuffer<unsigned>(nlistatoms));
    for (unsigned i = 0; i < nlistatoms; i++) {
      gpu->pbpLJ1264AtomList->_pSysData[i] = list[i];
    }
    gpu->pbpLJ1264AtomList->Upload();
    gpu->sim.ppLJ1264AtomList = gpu->pbpLJ1264AtomList->_pDevData;
    unsigned maxNumberpLJ1264NBEntries = gpu->sim.numberpLJ1264Atoms * gti_simulationConst::MaxNumberNBPerAtom;
    gpu->pbpLJ1264NBList = std::unique_ptr< GpuBuffer<int4> >(new GpuBuffer<int4>(maxNumberpLJ1264NBEntries));
    for (unsigned i = 0; i < maxNumberpLJ1264NBEntries; i++) {
      gpu->pbpLJ1264NBList->_pSysData[i] = { -1, -1, -1, -1 };
    }
    gpu->pbpLJ1264NBList->Upload();
    gpu->sim.ppLJ1264NBList = gpu->pbpLJ1264NBList->_pDevData;
    gpu->pbNumberpLJ1264NBEntries = std::unique_ptr< GpuBuffer<unsigned long long int> >(new GpuBuffer<unsigned long long int>(1));
    gpu->pbNumberpLJ1264NBEntries->_pSysData[0] = maxNumberpLJ1264NBEntries;
    gpu->pbNumberpLJ1264NBEntries->Upload();
    gpu->sim.pNumberpLJ1264NBEntries = gpu->pbNumberpLJ1264NBEntries->_pDevData;
  }

  delete[] list;
  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gti_lj1264plj1264_nb_setup_:
//
// Arguments:
//   atm_isymbl:
//   ntypes:
//   atm_iac:
//   ico:
//   gbl_cn6:
//   gbl_cn7: // C4PairwiseCUDA
//   gbl_cn8:
//   C4Pairwise:
//---------------------------------------------------------------------------------------------
extern "C" void gti_lj1264plj1264_nb_setup_(char atm_isymbl[][4], int* ntypes, int atm_iac[],
  int ico[], double gbl_cn6[], int gbl_cn7[], double gbl_cn8[], int* pC4Pairwise)
{
  PRINTMETHOD(__func__);

  gpuContext gpu = theGPUContext::GetPointer();
  unsigned natoms = gpu->sim.atoms;
  unsigned nlistatoms = *pC4Pairwise*2; //C4PairwiseCUDA
  unsigned* list = new unsigned[*pC4Pairwise*2+natoms];
  for (unsigned k=0;k<*pC4Pairwise;k++) { //C4PairwiseCUDA2023
    //if (i==(gbl_cn7[k*3]+3)/3 || i==(gbl_cn7[k*3+1]+3)/3) { //if any pairwise matches, check whether already counted
      list[k*2] = (gbl_cn7[k*3])/3;
      list[k*2+1]=(gbl_cn7[k*3+1])/3;
    //}
  }
  for (unsigned i = 0; i < natoms; i++) {
    if (atm_iac[i] > 0) {
      for (unsigned j = 1; j < 4; j++) {
        char tt = atm_isymbl[i][j]; 
        if (tt == '+' || tt == '-') { 
          int index = ico[(*ntypes) * (atm_iac[i] - 1) + atm_iac[i] - 1] - 1; 
          if (index >= 0)
            if (gbl_cn6[index] > 1e-10) {
              list[nlistatoms] = i;
              nlistatoms++;
              break;
            }
        }
      }
    }
  }
  gpu->sim.numberLJ1264pLJ1264Atoms = nlistatoms;
  if (nlistatoms > 0) {
    gpu->pbLJ1264pLJ1264AtomList = std::unique_ptr< GpuBuffer<unsigned> >(new GpuBuffer<unsigned>(nlistatoms));
    for (unsigned i = 0; i < nlistatoms; i++) {
      gpu->pbLJ1264pLJ1264AtomList->_pSysData[i] = list[i];
    }
    gpu->pbLJ1264pLJ1264AtomList->Upload();
    gpu->sim.pLJ1264pLJ1264AtomList = gpu->pbLJ1264pLJ1264AtomList->_pDevData;
    unsigned maxNumberLJ1264pLJ1264NBEntries = gpu->sim.numberLJ1264pLJ1264Atoms * gti_simulationConst::MaxNumberNBPerAtom;
    gpu->pbLJ1264pLJ1264NBList = std::unique_ptr< GpuBuffer<int4> >(new GpuBuffer<int4>(maxNumberLJ1264pLJ1264NBEntries));
    for (unsigned i = 0; i < maxNumberLJ1264pLJ1264NBEntries; i++) {
      gpu->pbLJ1264pLJ1264NBList->_pSysData[i] = { -1, -1, -1, -1 };
    }
    gpu->pbLJ1264pLJ1264NBList->Upload();
    gpu->sim.pLJ1264pLJ1264NBList = gpu->pbLJ1264pLJ1264NBList->_pDevData;
    gpu->pbNumberLJ1264pLJ1264NBEntries = std::unique_ptr< GpuBuffer<unsigned long long int> >(new GpuBuffer<unsigned long long int>(1));
    gpu->pbNumberLJ1264pLJ1264NBEntries->_pSysData[0] = maxNumberLJ1264pLJ1264NBEntries;
    gpu->pbNumberLJ1264pLJ1264NBEntries->Upload();
    gpu->sim.pNumberLJ1264pLJ1264NBEntries = gpu->pbNumberLJ1264pLJ1264NBEntries->_pDevData;
  }
  delete[] list;
  gpuCopyConstants();
}
//---------------------------------------------------------------------------------------------
// gti_pme_setup_:
//---------------------------------------------------------------------------------------------
extern "C" void gti_pme_setup_()
{
  PRINTMETHOD("gti_pme_setup");

  // PME part
  gpuContext gpu = theGPUContext::GetPointer();
  gpu->pbPMEChargeGrid = std::unique_ptr< GpuBuffer<unsigned long long int> > (new GpuBuffer<unsigned long long int>(gpu->sim.XYZStride * 3));
  for (unsigned i = 0; i < gpu->sim.XYZStride * 3; i++) {
    gpu->pbPMEChargeGrid->_pSysData[i] = 0;
  }
  gpu->pbPMEChargeGrid->Upload();
  gpu->sim.pPMEChargeGrid = gpu->pbPMEChargeGrid->_pDevData;


  gpu->pbFFTChargeGrid = std::unique_ptr< GpuBuffer<PMEFloat> > (new GpuBuffer<PMEFloat>(gpu->sim.XYZStride * 2));
  for (unsigned i = 0; i < gpu->sim.XYZStride * 2; i++) {
    gpu->pbFFTChargeGrid->_pSysData[i] = 0;
  }
  gpu->pbFFTChargeGrid->Upload();
  gpu->sim.pFFTChargeGrid = gpu->pbFFTChargeGrid->_pDevData;
}

//---------------------------------------------------------------------------------------------
// gti_bonded_setup_:
//---------------------------------------------------------------------------------------------
extern "C" void gti_bonded_setup_(int* pBat_type, int* reaf_mode, int* pError)
{
  gpuContext gpu = theGPUContext::GetPointer();
  unsigned aa[4];
  unsigned currentTIRegion;
  bool inTI = false;

  std::vector<std::pair<unsigned, unsigned> > rotatableBond; rotatableBond.clear();

  // Bonds
  if (true) {
    unsigned numberBond = gpu->sim.bonds;
    std::vector<PMEDouble2> tempBond; tempBond.clear();
    std::vector<uint2> tempID; tempID.clear();
    std::vector<uint> tempType; tempType.clear();

    for (unsigned int i = 0; i < numberBond; i++) {
      aa[0] = gpu->pbBondID->_pSysData[i].x;
      aa[1] = gpu->pbBondID->_pSysData[i].y;
      inTI = false;
      for (unsigned int j = 0; j < 2; j++) {
        bool b0 = (gpu->pbTIList->_pSysData[aa[j]])>0;
        bool b1 = (gpu->pbTIList->_pSysData[aa[j] + gpu->sim.stride]) > 0;
        if (b0 || b1) {
          inTI = true;
          currentTIRegion = (b0) ? 0 : 1; 
          break;
        }
      }
      if (inTI) {
        PMEDouble2 tb = gpu->pbBond->_pSysData[i];
        tempBond.push_back(tb);
        uint2 tbid;
        tbid.x = aa[0];
        tbid.y = aa[1];
        tempID.push_back(tbid);
        tempType.push_back(currentTIRegion);

        // Zero out the force constant
        gpu->pbBond->_pSysData[i].x = 0.0;
      }
    }
    gpu->sim.numberTIBond = tempBond.size(); 
    uint regionCounter[2] = { 0, 0 };
    if (*pBat_type == 2) {
      for (unsigned int i = 0; i < gpu->sim.numberTIBond; i++) {
        aa[0] = tempID[i].x;
        aa[1] = tempID[i].y;
        uint currentTIRegion = (tempType[i] & 1);
        uint& type = tempType[i];
        if ( ((gpu->pbSCBATList->_pSysData[aa[0]] & 2) && (gpu->pbSCBATList->_pSysData[aa[1]] & 2) && currentTIRegion == 0) ||
          ((gpu->pbSCBATList->_pSysData[aa[0]] & 16) && (gpu->pbSCBATList->_pSysData[aa[1]] & 16) && currentTIRegion == 1)  ){
          type |= gti_simulationConst::ex_userDefined_bat;
          regionCounter[currentTIRegion]++;
        }
      }
    }
    for (unsigned int i = 0; i < gpu->sim.numberTIBond; i++) {
      aa[0] = tempID[i].x;
      aa[1] = tempID[i].y;
      uint currentTIRegion = (tempType[i] & 1);
      uint& type = tempType[i];

      bool hasSC = (gpu->pbSCList->_pSysData[aa[0]] > 0 || gpu->pbSCList->_pSysData[aa[1]] > 0);
      bool inSC = (gpu->pbSCList->_pSysData[aa[0]] > 0 && gpu->pbSCList->_pSysData[aa[1]] > 0);
      if (hasSC) type |= gti_simulationConst::Fg_has_SC;
      if (inSC) type |= gti_simulationConst::Fg_int_SC;
      if (!hasSC) type |= gti_simulationConst::ex_addToDVDL_bat;

      if (hasSC && !inSC) {
        if ((*pBat_type) >= 1) {
          if (hasSC && !inSC) {
            if (regionCounter[currentTIRegion] == 0) {
              regionCounter[currentTIRegion]++;  // pick up the first CC-SC bond
            }
            else {
              if (!(type & gti_simulationConst::ex_userDefined_bat)) {
                type |= gti_simulationConst::ex_addToBat_corr;
                if (*pBat_type > 0) type |= gti_simulationConst::ex_addToDVDL_bat;
              }
            }
          }
        }
      }

      // search for rotatable bonds
      // ignore any bonds involving hydrogens
      if (gpu->pbAtomMass->_pSysData[aa[0]] > 3.1 && gpu->pbAtomMass->_pSysData[aa[1]] > 3.1) {

        std::vector<unsigned> c0, c1;
        c0.clear(); c1.clear();
        c0.push_back(aa[0]);
        c1.push_back(aa[1]);

        int numberDihedral = gpu->sim.dihedrals;
        bool inTorsion = false;
        for (unsigned int i = 0; i < numberDihedral; i++) {
          //a[0] is the termial atom
          if (aa[0] == gpu->pbDihedralID1->_pSysData[i].x
            && aa[1] != gpu->pbDihedralID1->_pSysData[i].y
            && aa[1] != gpu->pbDihedralID1->_pSysData[i].z
            && aa[1] != gpu->pbDihedralID1->_pSysData[i].w
            ) {
            c0.push_back(gpu->pbDihedralID1->_pSysData[i].y);
            c0.push_back(gpu->pbDihedralID1->_pSysData[i].z);
            c0.push_back(gpu->pbDihedralID1->_pSysData[i].w);
          }
          if (aa[0] == gpu->pbDihedralID1->_pSysData[i].w
            && aa[1] != gpu->pbDihedralID1->_pSysData[i].x
            && aa[1] != gpu->pbDihedralID1->_pSysData[i].y
            && aa[1] != gpu->pbDihedralID1->_pSysData[i].z
            ) {
            c0.push_back(gpu->pbDihedralID1->_pSysData[i].x);
            c0.push_back(gpu->pbDihedralID1->_pSysData[i].y);
            c0.push_back(gpu->pbDihedralID1->_pSysData[i].z);
          }

          //a[1] is the terminal atom 
          if (aa[1] == gpu->pbDihedralID1->_pSysData[i].x
            && aa[0] != gpu->pbDihedralID1->_pSysData[i].y
            && aa[0] != gpu->pbDihedralID1->_pSysData[i].z
            && aa[0] != gpu->pbDihedralID1->_pSysData[i].w
            ) {
            c1.push_back(gpu->pbDihedralID1->_pSysData[i].y);
            c1.push_back(gpu->pbDihedralID1->_pSysData[i].z);
            c1.push_back(gpu->pbDihedralID1->_pSysData[i].w);
          }
          if (aa[1] == gpu->pbDihedralID1->_pSysData[i].w
            && aa[0] != gpu->pbDihedralID1->_pSysData[i].x
            && aa[0] != gpu->pbDihedralID1->_pSysData[i].y
            && aa[0] != gpu->pbDihedralID1->_pSysData[i].z
            ) {
            c1.push_back(gpu->pbDihedralID1->_pSysData[i].x);
            c1.push_back(gpu->pbDihedralID1->_pSysData[i].y);
            c1.push_back(gpu->pbDihedralID1->_pSysData[i].z);
          }

          if (aa[0] == gpu->pbDihedralID1->_pSysData[i].y
            && aa[1] == gpu->pbDihedralID1->_pSysData[i].z)  inTorsion = true;
          if (aa[0] == gpu->pbDihedralID1->_pSysData[i].z
            && aa[1] == gpu->pbDihedralID1->_pSysData[i].y)  inTorsion = true;

        }

        bool rotatable = true;
        for (unsigned i = 0; i < c0.size(); i++) {
          rotatable = (std::find(begin(c1), end(c1), c0[i]) == end(c1));
          if (!rotatable) break;
        }

        if (rotatable && inTorsion) rotatableBond.push_back({ aa[0], aa[1] });
      }
    }

    gpu->pbTIBond = std::unique_ptr< GpuBuffer<PMEDouble2> > (new GpuBuffer<PMEDouble2>(gpu->sim.numberTIBond));
    gpu->pbTIBondID = std::unique_ptr< GpuBuffer<uint2> > (new GpuBuffer<uint2>(gpu->sim.numberTIBond));
    gpu->pbTIBondType = std::unique_ptr< GpuBuffer<uint> > (new GpuBuffer<uint>(gpu->sim.numberTIBond));
    for (unsigned int i = 0; i < gpu->sim.numberTIBond; i++) {
      gpu->pbTIBond->_pSysData[i] = tempBond[i];
      gpu->pbTIBondID->_pSysData[i] = tempID[i];
      gpu->pbTIBondType->_pSysData[i] = tempType[i];
    }
    gpu->pbTIBond->Upload();
    gpu->pbTIBondID->Upload();
    gpu->pbTIBondType->Upload();
    gpu->sim.pTIBond = gpu->pbTIBond->_pDevData;
    gpu->sim.pTIBondID = gpu->pbTIBondID->_pDevData;
    gpu->sim.pTIBondType = gpu->pbTIBondType->_pDevData;
  }

  // Angles
  if (true) {
    int numberAngle = gpu->sim.bondAngles;
    std::vector<PMEDouble2> tempAngle; tempAngle.clear();
    std::vector<uint4> tempID; tempID.clear();
    std::vector<uint> tempType; tempType.clear();

    for (unsigned int i = 0; i < numberAngle; i++) {
      aa[0] = gpu->pbBondAngleID1->_pSysData[i].x;
      aa[1] = gpu->pbBondAngleID1->_pSysData[i].y;
      aa[2] = gpu->pbBondAngleID2->_pSysData[i];
      inTI = false;
      for (unsigned int j = 0; j < 3; j++) {
        bool b0 = (gpu->pbTIList->_pSysData[aa[j]])>0;
        bool b1 = (gpu->pbTIList->_pSysData[aa[j] + gpu->sim.stride]) > 0;
        if (b0 || b1) {
          inTI = true;
          currentTIRegion = (b0) ? 0 : 1;
          break;
        }
      }
      if (*reaf_mode < 0 || (*reaf_mode != currentTIRegion)) {
        if (inTI) {
          PMEDouble2 ta = gpu->pbBondAngle->_pSysData[i];
          tempAngle.push_back(ta);
          uint4 tbid;
          tbid.x = aa[0];
          tbid.y = aa[1];
          tbid.z = aa[2];
          tempType.push_back(currentTIRegion);
          tempID.push_back(tbid);

          // Zero out the force constant
          gpu->pbBondAngle->_pSysData[i].x = 0.0;
        }
      }
    }

    gpu->sim.numberTIAngle = tempAngle.size();
    uint regionCounter[2] = { 0, 0 };
    if (*pBat_type == 2) {
      for (unsigned int i = 0; i < gpu->sim.numberTIAngle; i++) {
        aa[0] = tempID[i].x;
        aa[1] = tempID[i].y;
        aa[2] = tempID[i].z;
        uint currentTIRegion = (tempType[i] & 1);
        uint& type = tempType[i];
        if (((gpu->pbSCBATList->_pSysData[aa[0]] & 4) && (gpu->pbSCBATList->_pSysData[aa[1]] & 4) && (gpu->pbSCBATList->_pSysData[aa[2]] & 4) && currentTIRegion == 0) ||
          ((gpu->pbSCBATList->_pSysData[aa[0]] & 32) && (gpu->pbSCBATList->_pSysData[aa[1]] & 32) && (gpu->pbSCBATList->_pSysData[aa[2]] & 32) && currentTIRegion == 1)) {
          type |= gti_simulationConst::ex_userDefined_bat;
          regionCounter[currentTIRegion]++;
        }
      }
    }

    for (unsigned int i = 0; i < gpu->sim.numberTIAngle; i++) {
      aa[0] = tempID[i].x;
      aa[1] = tempID[i].y;
      aa[2] = tempID[i].z;
      uint currentTIRegion = (tempType[i] & 1);
      uint& type = tempType[i];

      bool b0 = gpu->pbSCList->_pSysData[aa[0]] > 0;
      bool b1 = gpu->pbSCList->_pSysData[aa[1]] > 0;
      bool b2 = gpu->pbSCList->_pSysData[aa[2]] > 0;
      bool hasSC = (b0 || b1 || b2);
      bool inSC = (b0 && b1 && b2);
      if (hasSC) type |= gti_simulationConst::Fg_has_SC;
      if (inSC) type |= gti_simulationConst::Fg_int_SC;
      if (!hasSC) type |= gti_simulationConst::ex_addToDVDL_bat;
      if (hasSC && !inSC) {
        if ( (*pBat_type) >= 1) {
            // RRD/DRR case
            if ((!b0 && !b1 && b2) || (b0 && !b1 && !b2)) {
              if (regionCounter[currentTIRegion] == 0) {
                regionCounter[currentTIRegion]++;
              } else 
              {
                if (!(type & gti_simulationConst::ex_userDefined_bat)) {
                  type |= gti_simulationConst::ex_addToBat_corr;
                  if (*pBat_type > 0) type |= gti_simulationConst::ex_addToDVDL_bat;
                }
              }
            } else {
              bool RDR = ( !b0 && b1 && !b2 );
              bool DRD = (b0 && !b1 && b2);
              if (RDR || DRD) {
                type |= gti_simulationConst::ex_addToBat_corr;
                if (*pBat_type > 0) type |= gti_simulationConst::ex_addToDVDL_bat;
              }
            }
        }
      }
    }

    gpu->pbTIBondAngle = std::unique_ptr< GpuBuffer<PMEDouble2> > (new GpuBuffer<PMEDouble2>(gpu->sim.numberTIAngle));
    gpu->pbTIBondAngleID = std::unique_ptr< GpuBuffer<uint4> > (new GpuBuffer<uint4>(gpu->sim.numberTIAngle));
    gpu->pbTIBondAngleType = std::unique_ptr< GpuBuffer<uint> > (new GpuBuffer<uint>(gpu->sim.numberTIAngle));
    for (unsigned int i = 0; i < gpu->sim.numberTIAngle; i++) {
      gpu->pbTIBondAngle->_pSysData[i] = tempAngle[i];
      gpu->pbTIBondAngleID->_pSysData[i] = tempID[i];
      gpu->pbTIBondAngleType->_pSysData[i] = tempType[i];
    }
    gpu->pbTIBondAngle->Upload();
    gpu->pbTIBondAngleID->Upload();
    gpu->pbTIBondAngleType->Upload();
    gpu->sim.pTIBondAngle = gpu->pbTIBondAngle->_pDevData;
    gpu->sim.pTIBondAngleID = gpu->pbTIBondAngleID->_pDevData;
    gpu->sim.pTIBondAngleType = gpu->pbTIBondAngleType->_pDevData;
  }

  // Dihedrals
  if (true) {
    int numberDihedral = gpu->sim.dihedrals;
    std::vector<PMEDouble2> tempDihedral1; tempDihedral1.clear();
    std::vector<PMEDouble2> tempDihedral2; tempDihedral2.clear();
    std::vector<PMEDouble> tempDihedral3; tempDihedral3.clear();
    std::vector<uint4> tempID; tempID.clear();
    std::vector<uint> tempDihedralRegion; tempDihedralRegion.clear();
    std::vector<uint> tempType; tempType.clear();

    for (unsigned int i = 0; i < numberDihedral; i++) {
      aa[0] = gpu->pbDihedralID1->_pSysData[i].x;
      aa[1] = gpu->pbDihedralID1->_pSysData[i].y;
      aa[2] = gpu->pbDihedralID1->_pSysData[i].z;
      aa[3] = gpu->pbDihedralID1->_pSysData[i].w;
      inTI = false;
      for (unsigned int j = 0; j < 4; j++) {
        bool b0 = (gpu->pbTIList->_pSysData[aa[j]])>0;
        bool b1 = (gpu->pbTIList->_pSysData[aa[j] + gpu->sim.stride]) > 0;
        if (b0 || b1) {
          inTI = true;
          currentTIRegion = (b0) ? 0 : 1;
          break;
        }
      }
      if (*reaf_mode < 0 || (*reaf_mode != currentTIRegion)) {
        if (inTI) {
          PMEDouble2 td1 = { gpu->pbDihedral1->_pSysData[i].x, gpu->pbDihedral1->_pSysData[i].y };
          tempDihedral1.push_back(td1);
          PMEDouble2 td2 = { gpu->pbDihedral2->_pSysData[i].x, gpu->pbDihedral2->_pSysData[i].y };
          tempDihedral2.push_back(td2);
          PMEDouble td3 = gpu->pbDihedral3->_pSysData[i];
          tempDihedral3.push_back(td3);
          uint4 tt = { aa[0], aa[1], aa[2], aa[3] };
          tempID.push_back(tt);
          tempType.push_back(currentTIRegion);

          // Zero out the force constant
          gpu->pbDihedral2->_pSysData[i].x = ZeroF;
          gpu->pbDihedral2->_pSysData[i].y = ZeroF;
          gpu->pbDihedral3->_pSysData[i] = ZeroF;
        }
      }
    }

    gpu->sim.numberTIDihedral = tempDihedral1.size();
    uint regionCounter[2] = { 0, 0 };
    if (*pBat_type == 2) {
      for (unsigned int i = 0; i < gpu->sim.numberTIDihedral; i++) {
        aa[0] = tempID[i].x;
        aa[1] = tempID[i].y;
        aa[2] = tempID[i].z;
        aa[3] = tempID[i].w;
        uint currentTIRegion = (tempType[i] & 1);
        uint& type = tempType[i];
        if (((gpu->pbSCBATList->_pSysData[aa[0]] & 8) && (gpu->pbSCBATList->_pSysData[aa[1]] & 8)
          && (gpu->pbSCBATList->_pSysData[aa[2]] & 8) && (gpu->pbSCBATList->_pSysData[aa[3]] & 8)  && currentTIRegion == 0) ||
            ((gpu->pbSCBATList->_pSysData[aa[0]] & 64) && (gpu->pbSCBATList->_pSysData[aa[1]] & 64) 
          && (gpu->pbSCBATList->_pSysData[aa[2]] & 64) && (gpu->pbSCBATList->_pSysData[aa[3]] & 64) && currentTIRegion == 1)) {
          type |= gti_simulationConst::ex_userDefined_bat;
          regionCounter[currentTIRegion]++;
        }
      }
    }

    for (unsigned int i = 0; i < gpu->sim.numberTIDihedral; i++) {
      aa[0] = tempID[i].x;
      aa[1] = tempID[i].y;
      aa[2] = tempID[i].z;
      aa[3] = tempID[i].w;
      uint currentTIRegion = (tempType[i] & 1);
      uint& type = tempType[i];

      bool b0 = gpu->pbSCList->_pSysData[aa[0]] > 0;
      bool b1 = gpu->pbSCList->_pSysData[aa[1]] > 0;
      bool b2 = gpu->pbSCList->_pSysData[aa[2]] > 0;
      bool b3 = gpu->pbSCList->_pSysData[aa[3]] > 0;
      bool hasSC = (b0 || b1 || b2 || b3);
      bool inSC = (b0 && b1 && b2 && b3);
      if (hasSC) type |= gti_simulationConst::Fg_has_SC;
      if (inSC) type |= gti_simulationConst::Fg_int_SC;

      // All torsions should be scaled when gti_add_sc>=4, =25,35
      if (gpu->sim.addSC >= 4 && (gpu->sim.addSC!=25) && (gpu->sim.addSC != 35)) {
        type |= gti_simulationConst::ex_addToDVDL_bat;
      } else {

        if (!hasSC) { // All CC torsions should be scaled
          type |= gti_simulationConst::ex_addToDVDL_bat;

        } else {// Treatment of SC 

          // look for Rotatable bonds
          if ( (gpu->sim.addSC== 25 || gpu->sim.addSC == 35) && (b1||b2) ) {
            std::pair <unsigned, unsigned> qq = { unsigned(aa[1]), unsigned(aa[2]) };
            if (std::find(begin(rotatableBond), end(rotatableBond),qq) != end(rotatableBond) ) {
              type |= gti_simulationConst::ex_addToDVDL_bat;
            }
            qq = { unsigned(aa[2]), unsigned(aa[1]) };
            if (std::find(begin(rotatableBond), end(rotatableBond), qq) != end(rotatableBond)) {
              type |= gti_simulationConst::ex_addToDVDL_bat;
            }
          }

          // CC-SC boundar
          if (!inSC) {  
            if ( (*pBat_type) >= 1) {
            if (hasSC && !inSC) {
              // RRDD/DDRR case
              if ((!b0 && !b1 && b2 && b3) || (b0 && b1 && !b2 && !b3)) {
                if (regionCounter[currentTIRegion] == 0) {
                  regionCounter[currentTIRegion]++;
                } else {
                  if (!(type & gti_simulationConst::ex_userDefined_bat)) {
                    type |= gti_simulationConst::ex_addToBat_corr;
                    if (*pBat_type > 0) type |= gti_simulationConst::ex_addToDVDL_bat;
                  }
                }
              } else {
                // RDDD/DDDR case
                bool RDDD = ((!b0 && b1 && b2 && b3) || (b0 && b1 && b2 && b3));
                if (!RDDD) {
                  type |= gti_simulationConst::ex_addToBat_corr;
                  if (*pBat_type > 0) type |= gti_simulationConst::ex_addToDVDL_bat;
                }
              }
            }
          }
        }

        }

      }
    }


    gpu->pbTIDihedral1 = std::unique_ptr< GpuBuffer<PMEDouble2> > (new GpuBuffer<PMEDouble2>(gpu->sim.numberTIDihedral));
    gpu->pbTIDihedral2 = std::unique_ptr< GpuBuffer<PMEDouble2> > (new GpuBuffer<PMEDouble2>(gpu->sim.numberTIDihedral));
    gpu->pbTIDihedral3 = std::unique_ptr< GpuBuffer<PMEDouble> > (new GpuBuffer<PMEDouble>(gpu->sim.numberTIDihedral));
    gpu->pbTIDihedralID = std::unique_ptr< GpuBuffer<uint4> > (new GpuBuffer<uint4>(gpu->sim.numberTIDihedral));
    gpu->pbTIDihedralType = std::unique_ptr< GpuBuffer<uint> > (new GpuBuffer<uint>(gpu->sim.numberTIDihedral));
    for (unsigned int i = 0; i < gpu->sim.numberTIDihedral; i++) {
      gpu->pbTIDihedral1->_pSysData[i] = tempDihedral1[i];
      gpu->pbTIDihedral2->_pSysData[i] = tempDihedral2[i];
      gpu->pbTIDihedral3->_pSysData[i] = tempDihedral3[i];
      gpu->pbTIDihedralID->_pSysData[i] = tempID[i];
      gpu->pbTIDihedralType->_pSysData[i] = tempType[i];
    }
    gpu->pbTIDihedral1->Upload();
    gpu->pbTIDihedral2->Upload();
    gpu->pbTIDihedral3->Upload();
    gpu->pbTIDihedralID->Upload();
    gpu->pbTIDihedralType->Upload();
    gpu->sim.pTIDihedral1 = gpu->pbTIDihedral1->_pDevData;
    gpu->sim.pTIDihedral2 = gpu->pbTIDihedral2->_pDevData;
    gpu->sim.pTIDihedral3 = gpu->pbTIDihedral3->_pDevData;
    gpu->sim.pTIDihedralID = gpu->pbTIDihedralID->_pDevData;
    gpu->sim.pTIDihedralType = gpu->pbTIDihedralType->_pDevData;
  }
}

//---------------------------------------------------------------------------------------------
// gti_reaf_bonded_setup_:
//---------------------------------------------------------------------------------------------
extern "C" void gti_reaf_bonded_setup_(int* pError) {
  gpuContext gpu = theGPUContext::GetPointer();
  unsigned aa[4];
  bool hasRE = false;
  bool intRE = false;

  // Angles
  if (true) {
    int numberAngle = gpu->sim.bondAngles;
    std::vector<PMEDouble2> tempAngle; tempAngle.clear();
    std::vector<uint4> tempID; tempID.clear();
    std::vector<uint> tempType; tempType.clear();

    for (unsigned int i = 0; i < numberAngle; i++) {
      aa[0] = gpu->pbBondAngleID1->_pSysData[i].x;
      aa[1] = gpu->pbBondAngleID1->_pSysData[i].y;
      aa[2] = gpu->pbBondAngleID2->_pSysData[i];
      hasRE = false;
      for (unsigned int j = 0; j < 3; j++) {
        hasRE = ((gpu->pbREAFList->_pSysData[aa[j]]) > 0);
        if (hasRE) break;
      }
      if (hasRE) {
        PMEDouble2 ta = gpu->pbBondAngle->_pSysData[i];
        tempAngle.push_back(ta);
        uint4 tbid;
        tbid.x = aa[0];
        tbid.y = aa[1];
        tbid.z = aa[2];
        tempType.push_back(0);
        tempID.push_back(tbid);

        // Zero out theforce constant
        // gpu->pbBondAngle->_pSysData[i].x = 0.0;
      }
    }

    gpu->sim.numberREAFAngle = tempAngle.size();
    uint regionCounter[2] = { 0, 0 };

    for (unsigned int i = 0; i < gpu->sim.numberREAFAngle; i++) {
      aa[0] = tempID[i].x;
      aa[1] = tempID[i].y;
      aa[2] = tempID[i].z;
      uint currentTIRegion = (tempType[i] & 1);
      uint& type = tempType[i];

      bool b0 = gpu->pbREAFList->_pSysData[aa[0]] > 0;
      bool b1 = gpu->pbREAFList->_pSysData[aa[1]] > 0;
      bool b2 = gpu->pbREAFList->_pSysData[aa[2]] > 0;
      bool hasRE = (b0 || b1 || b2);
      bool intRE = (b0 && b1 && b2);
      if (hasRE) type |= gti_simulationConst::Fg_has_RE;
      if (intRE) type |= gti_simulationConst::Fg_int_RE;
      if (gpu->sim.addRE >= 7) {
        type |= gti_simulationConst::ex_addRE_bat;
      }
    }

    gpu->pbREAFBondAngle = std::unique_ptr< GpuBuffer<PMEDouble2> >(new GpuBuffer<PMEDouble2>(gpu->sim.numberREAFAngle));
    gpu->pbREAFBondAngleID = std::unique_ptr< GpuBuffer<uint4> >(new GpuBuffer<uint4>(gpu->sim.numberREAFAngle));
    gpu->pbREAFBondAngleType = std::unique_ptr< GpuBuffer<uint> >(new GpuBuffer<uint>(gpu->sim.numberREAFAngle));
    for (unsigned int i = 0; i < gpu->sim.numberREAFAngle; i++) {
      gpu->pbREAFBondAngle->_pSysData[i] = tempAngle[i];
      gpu->pbREAFBondAngleID->_pSysData[i] = tempID[i];
      gpu->pbREAFBondAngleType->_pSysData[i] = tempType[i];
    }
    gpu->pbREAFBondAngle->Upload();
    gpu->pbREAFBondAngleID->Upload();
    gpu->pbREAFBondAngleType->Upload();
    gpu->sim.pREAFBondAngle = gpu->pbREAFBondAngle->_pDevData;
    gpu->sim.pREAFBondAngleID = gpu->pbREAFBondAngleID->_pDevData;
    gpu->sim.pREAFBondAngleType = gpu->pbREAFBondAngleType->_pDevData;
  }

  // Dihedrals
  if (true) {
    int numberDihedral = gpu->sim.dihedrals;
    std::vector<PMEDouble2> tempDihedral1; tempDihedral1.clear();
    std::vector<PMEDouble2> tempDihedral2; tempDihedral2.clear();
    std::vector<PMEDouble> tempDihedral3; tempDihedral3.clear();
    std::vector<uint4> tempID; tempID.clear();
    std::vector<uint> tempDihedralRegion; tempDihedralRegion.clear();
    std::vector<uint> tempType; tempType.clear();

    for (unsigned int i = 0; i < numberDihedral; i++) {
      aa[0] = gpu->pbDihedralID1->_pSysData[i].x;
      aa[1] = gpu->pbDihedralID1->_pSysData[i].y;
      aa[2] = gpu->pbDihedralID1->_pSysData[i].z;
      aa[3] = gpu->pbDihedralID1->_pSysData[i].w;
      for (unsigned int j = 0; j < 4; j++) {
        hasRE = (gpu->pbREAFList->_pSysData[aa[j]]) > 0;
        if (hasRE) break;
      }
      if (hasRE) {
        PMEDouble2 td1 = { gpu->pbDihedral1->_pSysData[i].x, gpu->pbDihedral1->_pSysData[i].y };
        tempDihedral1.push_back(td1);
        PMEDouble2 td2 = { gpu->pbDihedral2->_pSysData[i].x, gpu->pbDihedral2->_pSysData[i].y };
        tempDihedral2.push_back(td2);
        PMEDouble td3 = gpu->pbDihedral3->_pSysData[i];
        tempDihedral3.push_back(td3);
        uint4 tt = { aa[0], aa[1], aa[2], aa[3] };
        tempID.push_back(tt);
        tempType.push_back(0);
      }
    }

    gpu->sim.numberREAFDihedral = tempDihedral1.size();

    for (unsigned int i = 0; i < gpu->sim.numberREAFDihedral; i++) {
      aa[0] = tempID[i].x;
      aa[1] = tempID[i].y;
      aa[2] = tempID[i].z;
      aa[3] = tempID[i].w;
      uint currentTIRegion = (tempType[i] & 1);
      uint& type = tempType[i];

      bool b0 = gpu->pbREAFList->_pSysData[aa[0]] > 0;
      bool b1 = gpu->pbREAFList->_pSysData[aa[1]] > 0;
      bool b2 = gpu->pbREAFList->_pSysData[aa[2]] > 0;
      bool b3 = gpu->pbREAFList->_pSysData[aa[3]] > 0;
      bool hasRE = (b0 || b1 || b2 || b3);
      bool intRE = (b0 && b1 && b2 && b3);
      if (hasRE) type |= gti_simulationConst::Fg_has_RE;
      if (intRE) type |= gti_simulationConst::Fg_int_RE;

      // All torsions should be scaled when gti_add_re >=4
      if (gpu->sim.addRE >= 4) {
        type |= gti_simulationConst::ex_addRE_bat;
      } 
    }


    gpu->pbREAFDihedral1 = std::unique_ptr< GpuBuffer<PMEDouble2> >(new GpuBuffer<PMEDouble2>(gpu->sim.numberREAFDihedral));
    gpu->pbREAFDihedral2 = std::unique_ptr< GpuBuffer<PMEDouble2> >(new GpuBuffer<PMEDouble2>(gpu->sim.numberREAFDihedral));
    gpu->pbREAFDihedral3 = std::unique_ptr< GpuBuffer<PMEDouble> >(new GpuBuffer<PMEDouble>(gpu->sim.numberREAFDihedral));
    gpu->pbREAFDihedralID = std::unique_ptr< GpuBuffer<uint4> >(new GpuBuffer<uint4>(gpu->sim.numberREAFDihedral));
    gpu->pbREAFDihedralType = std::unique_ptr< GpuBuffer<uint> >(new GpuBuffer<uint>(gpu->sim.numberREAFDihedral));
    for (unsigned int i = 0; i < gpu->sim.numberREAFDihedral; i++) {
      gpu->pbREAFDihedral1->_pSysData[i] = tempDihedral1[i];
      gpu->pbREAFDihedral2->_pSysData[i] = tempDihedral2[i];
      gpu->pbREAFDihedral3->_pSysData[i] = tempDihedral3[i];
      gpu->pbREAFDihedralID->_pSysData[i] = tempID[i];
      gpu->pbREAFDihedralType->_pSysData[i] = tempType[i];
    }
    gpu->pbREAFDihedral1->Upload();
    gpu->pbREAFDihedral2->Upload();
    gpu->pbREAFDihedral3->Upload();
    gpu->pbREAFDihedralID->Upload();
    gpu->pbREAFDihedralType->Upload();
    gpu->sim.pREAFDihedral1 = gpu->pbREAFDihedral1->_pDevData;
    gpu->sim.pREAFDihedral2 = gpu->pbREAFDihedral2->_pDevData;
    gpu->sim.pREAFDihedral3 = gpu->pbREAFDihedral3->_pDevData;
    gpu->sim.pREAFDihedralID = gpu->pbREAFDihedralID->_pDevData;
    gpu->sim.pREAFDihedralType = gpu->pbREAFDihedralType->_pDevData;
  }
}


extern "C" void gti_restraint_setup_(int* pBat_type)
{

  gpuContext gpu = theGPUContext::GetPointer();

  int* pNRestraint=NULL;
  unsigned currentTIRegion = 0;
  bool inTI = false;

  std::vector<RestraintType> tempType;
  std::vector<PMEDouble4> tempR, R[2];
  std::vector<PMEDouble2> tempK, K[2];
  std::vector<uint4> tempID, ID[2];
  std::vector<int> tempRegion;  // TI region
  std::vector<int> tempInfo;  // Aux info for marking entry being selected

  GpuBuffer <PMEDouble2>* pR1R2;
  GpuBuffer <PMEDouble2>* pR3R4;
  GpuBuffer <PMEDouble2>* pK2K3;
  GpuBuffer <int2>* pID12;
  GpuBuffer <int>* pID3;
  GpuBuffer <int4>* pID1234;
  unsigned len;
  
  for (RestraintType type = RDist; type <= RDihedral; type=RestraintType(type+1)) {
    tempR.clear();
    tempK.clear();
    tempID.clear();
    tempRegion.clear();
    tempInfo.clear();
    len = 0;

    switch(type) {
      case(RDist): 
        pNRestraint = &(gpu->sim.NMRDistances);
        pR1R2 = gpu->pbNMRDistanceR1R2Int;
        pR3R4 = gpu->pbNMRDistanceR3R4Int;
        pK2K3 = gpu->pbNMRDistanceK2K3Int;
        pID12 = gpu->pbNMRDistanceID;
        len = 2;
        break;
      case(RAngle): 
        pNRestraint= &(gpu->sim.NMRAngles);
        pR1R2 = gpu->pbNMRAngleR1R2Int;
        pR3R4 = gpu->pbNMRAngleR3R4Int;
        pK2K3 = gpu->pbNMRAngleK2K3Int;
        pID12 = gpu->pbNMRAngleID1;
        pID3 = gpu->pbNMRAngleID2;
        len = 3;
        break;
      case(RDihedral): 
        pNRestraint = &(gpu->sim.NMRTorsions);
        pR1R2 = gpu->pbNMRTorsionR1R2Int;
        pR3R4 = gpu->pbNMRTorsionR3R4Int;
        pK2K3 = gpu->pbNMRTorsionK2K3Int;
        pID1234 = gpu->pbNMRTorsionID1;
        len = 4;
        break;
    }

    tempR.clear(); tempK.clear(); tempID.clear(); tempInfo.clear(); tempRegion.clear();
    if (*pNRestraint > 0) {
      for (unsigned int i = 0; i < *pNRestraint; i++) {

        int tID[4] = { -1, -1, -1, -1 };
        switch (type) {
        case(RDist):
          tID[0] = pID12->_pSysData[i].x; tID[1] = pID12->_pSysData[i].y;
          break;
        case(RAngle):
          tID[0] = pID12->_pSysData[i].x; tID[1] = pID12->_pSysData[i].y;
          tID[2] = pID3->_pSysData[i];
          break;
        case(RDihedral):
          tID[0] = pID1234->_pSysData[i].x; tID[1] = pID1234->_pSysData[i].y;
          tID[2] = pID1234->_pSysData[i].z; tID[3] = pID1234->_pSysData[i].w;
          break;
        }

        bool scFlag = false;
        for (unsigned i = 0; i < len; i++)
          scFlag = (scFlag || (gpu->pbSCList->_pSysData[tID[i]] > 0));

        //Any restraint involved SC atoms will not be considered if gti_bat_sc=0
        if (*pBat_type == 0 && scFlag) continue;

        inTI = false;
        for (unsigned int j = 0; j < len; j++) {
          bool b0 = (gpu->pbTIList->_pSysData[tID[j]]) > 0;
          bool b1 = (gpu->pbTIList->_pSysData[tID[j] + gpu->sim.stride]) > 0;
          if (b0 || b1) {
            inTI = true;
            currentTIRegion = (b0) ? 0 : 1;
            break;
          }
        }
        
        if (inTI) {
          tempType.push_back(type);

          PMEDouble4 tr = {
            pR1R2->_pSysData[i].x,
            pR1R2->_pSysData[i].y,
            pR3R4->_pSysData[i].x,
            pR3R4->_pSysData[i].y };
          tempR.push_back(tr);

          PMEDouble2 tk = {
            pK2K3->_pSysData[i].x,
            pK2K3->_pSysData[i].y };
          tempK.push_back(tk);
          pK2K3->_pSysData[i] = { 0.0, 0.0 };

          uint4 tbid = {
            (uint)tID[0], (uint)tID[1], (uint)tID[2], (uint)tID[3] };
          tempID.push_back(tbid);

          tempRegion.push_back(currentTIRegion);
          if (scFlag) {
            tempInfo.push_back(1);
          } else {
            tempInfo.push_back(0);
          }
        }

      }
    }

    int tempNumberRestraint = tempR.size();
    R[0].clear(); K[0].clear(); ID[0].clear();
    R[1].clear(); K[1].clear(); ID[1].clear();

    if (tempNumberRestraint > 0) {
      int nTIPair = gpu->sim.numberTICommonPairs;
      int i0 = -1, i1 = -1;
      for (unsigned int i = 0; i < tempNumberRestraint; i++) {
        for (unsigned int j = 0; (j < tempNumberRestraint && j!=i && tempType[i]==tempType[j] ); j++) {

          bool found = false;
          if (tempInfo[i] == 0 && tempInfo[j] == 0 && ((tempRegion[i] + tempRegion[j]) == 1)) {
            if (gpu->GetMatchAtomID(tempID[i].x) == tempID[j].x) {
              if (gpu->GetMatchAtomID(tempID[i].y) == tempID[j].y) {
                if (tempType[i] == RDist) {
                  found = true;
                  break;
                } else if (gpu->GetMatchAtomID(tempID[i].z) == tempID[j].z) {
                  if (tempType[i] == RAngle) {
                    found = true;
                    break;
                  } else if (gpu->GetMatchAtomID(tempID[i].w) == tempID[j].w) {
                    if (tempType[i] == RDihedral) {
                      found = true;
                      break;
                    }
                  }
                }
              }
                 
            }
                     
          }
            // found the match pair: create a matched paired R, K, and ID records.
          if (found) {
            R[tempRegion[i]].push_back(tempR[i]); R[tempRegion[j]].push_back(tempR[j]);
            K[tempRegion[i]].push_back(tempK[i]); K[tempRegion[j]].push_back(tempK[j]);
            ID[tempRegion[i]].push_back(tempID[i]); ID[tempRegion[j]].push_back(tempID[j]);
            tempInfo[i] = -1;
            tempInfo[j] = -1;
          }
        }
      }

      for (unsigned int i = 0; i < tempNumberRestraint; i++) {
        // create an artifical matched partner
        if (tempInfo[i] >= 0) {
          unsigned myTI = tempRegion[i];
          R[myTI].push_back(tempR[i]); R[1 - myTI].push_back(tempR[i]);
          K[myTI].push_back(tempK[i]); K[1 - myTI].push_back(PMEDouble2({ 0.0, 0.0 }));

          uint4 tID = tempID[i];
          if (tempInfo[i] == 0) {
            tID = {
               gpu->GetMatchAtomID(tempID[i].x),
               gpu->GetMatchAtomID(tempID[i].y),
               gpu->GetMatchAtomID(tempID[i].z),
               gpu->GetMatchAtomID(tempID[i].w)
            };
          }
          ID[myTI].push_back(tempID[i]); ID[1 - myTI].push_back(tID);
        }
      }
    }

    unsigned n = R[0].size();
    if (n>0) {

      std::unique_ptr<GpuBuffer<PMEDouble4>> *pTIR = NULL;
      std::unique_ptr<GpuBuffer<PMEDouble2>> * pTIK = NULL;
      std::unique_ptr<GpuBuffer<uint4>> *pTIID = NULL;

      PMEDouble4** pSTIR = NULL;
      PMEDouble2** pSTIK = NULL;
      uint4** pSTIID = NULL;

      switch (type) {
      case(RDist):
        gpu->sim.numberTINMRDistance=n;
        pTIR = &(gpu->pbTINMRDistancesR);
        pTIK = &(gpu->pbTINMRDistancesK);
        pTIID = &(gpu->pbTINMRDistancesID);
        pSTIR = &(gpu->sim.pTINMRDistancesR);
        pSTIK = &(gpu->sim.pTINMRDistancesK);
        pSTIID = &(gpu->sim.pTINMRDistancesID);
        break;
      case(RAngle):
        gpu->sim.numberTINMRAngle = n;
        pTIR = &(gpu->pbTINMRAngleR);
        pTIK = &(gpu->pbTINMRAngleK);
        pTIID = &(gpu->pbTINMRAngleID);
        pSTIR = &(gpu->sim.pTINMRAngleR);
        pSTIK = &(gpu->sim.pTINMRAngleK);
        pSTIID = &(gpu->sim.pTINMRAngleID);
        break;
      case(RDihedral):
        gpu->sim.numberTINMRDihedral = n;
        pTIR = &(gpu->pbTINMRDihedralR);
        pTIK = &(gpu->pbTINMRDihedralK);
        pTIID = &(gpu->pbTINMRDihedralID);
        pSTIR = &(gpu->sim.pTINMRDihedralR);
        pSTIK = &(gpu->sim.pTINMRDihedralK);
        pSTIID = &(gpu->sim.pTINMRDihedralID);
        break;
      }

      (*pTIR) = std::unique_ptr< GpuBuffer<PMEDouble4> > (new GpuBuffer<PMEDouble4>(2 * n));
      (*pTIK) = std::unique_ptr< GpuBuffer<PMEDouble2> > (new GpuBuffer<PMEDouble2>(2 * n));
      (*pTIID) = std::unique_ptr< GpuBuffer<uint4> > (new GpuBuffer<uint4>(2 * n));

      for (unsigned int i = 0; i < n; i++) {
        (*pTIR)->_pSysData[i] = R[0][i];
        (*pTIR)->_pSysData[i + n] = R[1][i];

        (*pTIK)->_pSysData[i] = K[0][i];
        (*pTIK)->_pSysData[i + n] = K[1][i];

        (*pTIID)->_pSysData[i] = ID[0][i];
        (*pTIID)->_pSysData[i + n] = ID[1][i];
      }
      (*pTIR)->Upload();
      (*pTIK)->Upload();
      (*pTIID)->Upload();

      (*pSTIR) = (*pTIR)->_pDevData;
      (*pSTIK) = (*pTIK)->_pDevData;
      (*pSTIID) = (*pTIID)->_pDevData;

    }

  }

  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gti_build_nl_list_:
//
// Arguments:
//   ti_mode:  
//   lj1264:   
//   plj1264 C4PairwiseCUDA
//---------------------------------------------------------------------------------------------
extern "C" void gti_build_nl_list_(int* ti_mode, int* lj1264, int* plj1264,  int* reaf) 
{
  PRINTMETHOD(__func__);
  
  // This function should be called only after the new NB has been updated
  static bool firstTime = true;
  gpuContext gpu = theGPUContext::GetPointer();

  bool buildSpecial = gpu->bNeedNewNeighborList;
  gpu_build_neighbor_list_();
  if (gpu->multiStream) {
    cudaEventRecord(gpu->event_IterationBegin, gpu->mainStream);
  }
  if (buildSpecial || firstTime)  {
    if (*ti_mode > 0) {
      ik_BuildTINBList(gpu);
    }
    if (*lj1264 > 0 && *plj1264 == 0) {
      ik_Build1264NBList(gpu);
    }
    if (*plj1264 > 0 && *lj1264 == 0) { // C4PairwiseCUDA
      ik_Buildp1264NBList(gpu);
    }
    if (*plj1264 > 0 && *lj1264 > 0) { // C4PairwiseCUDA
      ik_Build1264p1264NBList(gpu);
    }
    if (*reaf >= 0) {
      ik_BuildREAFNbList(gpu);
    }

    firstTime = false;
  }  
}

//---------------------------------------------------------------------------------------------
// gti_build_nl_list_gamd_:
//
// Arguments:
//   ti_mode:  
//   lj1264:   
//   plj1264 C4PairwiseCUDA
//---------------------------------------------------------------------------------------------
extern "C" void gti_build_nl_list_gamd_(int* ti_mode, int* lj1264, int* plj1264)
{
  PRINTMETHOD("gti_build_nl_list_gamd");
  
  // This function should be called only after the new NB has been updated
  static bool firstTime = true;
  gpuContext gpu = theGPUContext::GetPointer();

  bool buildSpecial = gpu->bNeedNewNeighborList;
  gpu_build_neighbor_list_();
  if (gpu->multiStream) {
    cudaEventRecord(gpu->event_IterationBegin, gpu->mainStream);
  }
  if (buildSpecial || firstTime)  {
    if (gpu->multiStream) {
      cudaEventSynchronize(gpu->event_IterationBegin);
    }
    if (*ti_mode > 0) {
      ik_BuildTINBList_gamd(gpu);
    }
    if (*lj1264 > 0 && *plj1264 == 0 ) {
      ik_Build1264NBList_gamd(gpu);
    }
    if (*plj1264 > 0 && *lj1264 == 0) { // C4PairwiseCUDA
      ik_Buildp1264NBList_gamd(gpu);
    }
    if (*plj1264 > 0 && *lj1264 > 0) { // C4PairwiseCUDA
      ik_Build1264p1264NBList_gamd(gpu);
    }
    firstTime = false;
  }  
}

//---------------------------------------------------------------------------------------------
// gti_clear_gamd_:
//
// Arguments:
//   needEnergy:  
//   ti_mode:     
//   crd:         
//---------------------------------------------------------------------------------------------
extern "C" void gti_clear_gamd_(bool* needEnergy, int* ti_mode, double crd[][3])
{
  PRINTMETHOD("gti_clear");

  gpuContext gpu = theGPUContext::GetPointer();
  if (gpu->multiStream) cudaEventSynchronize(gpu->event_IterationBegin);

  // Clear forces        
  if ((*needEnergy) || ( (gpu->sim.ntp > 0) && (gpu->sim.barostat == 1) ) ) {
    kClearForces(gpu, gpu->sim.NLNonbondEnergyWarps);
  }
  else {
    kClearForces(gpu, gpu->sim.NLNonbondForcesWarps);
  }
  if (*ti_mode>0) {
    ik_ClearTIEnergyForce_gamd_kernel(gpu);
  }

  // Download critical coordinates for force calculations to occur on the host side
  if (gpu->sim.nShuttle > 0) {
    kRetrieveSimData(gpu, crd, 0);
  }
  kPMEGetGridWeights(gpu);
}

//---------------------------------------------------------------------------------------------
// gti_clear_:
//
// Arguments:
//   needEnergy:  
//   ti_mode:     
//   crd:         
//---------------------------------------------------------------------------------------------
extern "C" void gti_clear_(bool* needEnergy, int* ti_mode, double crd[][3])
{
  PRINTMETHOD("gti_clear");

  gpuContext gpu = theGPUContext::GetPointer();

  // Make sure the TIStream is catching up
  if (gpu->multiStream) cudaStreamSynchronize(gpu->TIStream);

  // Clear forces        
  if ((*needEnergy) || ( (gpu->sim.ntp > 0) && (gpu->sim.barostat == 1) ) ) {
    kClearForces(gpu, gpu->sim.NLNonbondEnergyWarps);
  }
  else {
    kClearForces(gpu, gpu->sim.NLNonbondForcesWarps);
  }
  if (*ti_mode>0) {
    ik_ClearTIEnergyForce_kernel(gpu);
  }

  // Download critical coordinates for force calculations to occur on the host side
  if (gpu->sim.nShuttle > 0) {
    kRetrieveSimData(gpu, crd, 0);
  }
  if (gpu->sim.iphmd > 0) {
    kClearDerivs_PHMD(gpu);
  }
  kPMEGetGridWeights(gpu);
}
//---------------------------------------------------------------------------------------------
// gti_bonded_:
//
// Arguments:
//   needPotEnergy:  
//   need_virials:   
//   ti_mode:        
//---------------------------------------------------------------------------------------------
extern "C" void gti_bonded_(bool* needPotEnergy, bool* need_virials, int* ti_mode, int* reaf_mode)
{
  PRINTMETHOD(__func__);

  gpuContext gpu = theGPUContext::GetPointer();

  if (gpu->iphmd > 0) {
    kExecuteBondWorkUnits(gpu, 1);
  }
  else if (!useMPI || gpu->bCalculateLocalForces) {
     kExecuteBondWorkUnits(gpu, (*needPotEnergy ? 1 : 0));
  }
  if (*ti_mode > 0) {
    ik_CalculateTIBonded(gpu, *needPotEnergy, *need_virials);
    ik_CalculateTIRestraint(gpu, *needPotEnergy, *need_virials);
  }

  if (*reaf_mode >= 0) {
    ik_CalculateREAFBonded(gpu, *needPotEnergy, *need_virials);
  }
}

//---------------------------------------------------------------------------------------------
// gti_bonded_ppi_:
//
// Arguments:
//   needPotEnergy:
//   need_virials:
//   ti_mode:
//---------------------------------------------------------------------------------------------
extern "C" void gti_bonded_ppi_(bool* needPotEnergy, bool* need_virials, int* ti_mode)
{
  PRINTMETHOD("gti_bonded");

  gpuContext gpu = theGPUContext::GetPointer();

  if (!useMPI || gpu->bCalculateLocalForces) {
    kExecuteBondWorkUnits(gpu, (*needPotEnergy ? 1 : 0));
  }
  if (*ti_mode > 0) {
    ik_CalculateTIBonded_ppi(gpu, *needPotEnergy, *need_virials);
  }
}


//---------------------------------------------------------------------------------------------
// gti_others_:
//
// Arguments:
//   needEnergy:  
//   nstep:       
//   dt:          
//---------------------------------------------------------------------------------------------
extern "C" void gti_others_(bool* needEnergy, int* nstep, double* dt, double* vol, double* ewaldcof)
{
  PRINTMETHOD(__func__);
  
  gpuContext gpu = theGPUContext::GetPointer();

  if (useMPI && !gpu->bCalculateLocalForces) {
    return;
  }
  if (gpu->ntf != 8) {
    *needEnergy ? kCalculateNMREnergy(gpu) : kCalculateNMRForces(gpu);
    if (gpu->sim.efx != 0 || gpu->sim.efy != 0 || gpu->sim.efz != 0) {
      *needEnergy ?
      kCalculateEFieldEnergy(gpu, *nstep, *dt) :
      kCalculateEFieldForces(gpu, *nstep, *dt);
    }
  }
  if(gpu->iphmd > 0) {
    if(!gpu->sim.qsetlam) {
        kModelForces_PHMD(gpu);
    }
  }
  if (gpu->iphmd == 3) {
    kCalculateSelfPHMD(*vol, *ewaldcof, gpu);
  }

#ifndef AMBER_PLATFORM_AMD  
  if (gpu->sim.number_rmsd_set > 0) {
    ik_CalculateRMSD(gpu, *needEnergy);
  }
#endif
}

//---------------------------------------------------------------------------------------------
// gti_update_md_ene_:
//
// Arguments:
//   pEnergy:  
//   enmr:     
//   virial:   
//   ekcmt:    
//---------------------------------------------------------------------------------------------
extern "C" void gti_update_md_ene_(pme_pot_ene_rec* pEnergy, double enmr[3], double virial[3],
                                   double ekcmt[3], int* ineb, int* nebfreq, int* nstep)
{
  PRINTMETHOD("gti_update_md_ene");

  gpuContext gpu = theGPUContext::GetPointer();
  double energy[EXTENDED_ENERGY_TERMS];

  bool masterNode= (!useMPI || gpu->gpuID == 0);

  icc_GetEnergyFromGPU(gpu, energy);
  pEnergy->vdw_dir = energy[1];
  pEnergy->vdw_recip = gpu->vdw_recip;
  pEnergy->vdw_tot = pEnergy->vdw_dir;
  if (masterNode) pEnergy->vdw_tot += pEnergy->vdw_recip;
    
  pEnergy->elec_dir = energy[10];
  pEnergy->elec_recip = energy[9];
  pEnergy->elec_nb_adjust = 0.0;
  pEnergy->elec_tot = pEnergy->elec_dir + pEnergy->elec_recip;
  if (masterNode) pEnergy->elec_tot += pEnergy->elec_self;

  pEnergy->hbond = 0.0;

  pEnergy->bond = energy[3];
  pEnergy->angle = energy[4];
  pEnergy->dihedral = energy[5];

  pEnergy->vdw_14 = energy[7];
  pEnergy->elec_14 = energy[6];

  pEnergy->restraint = energy[8] + energy[14] + energy[15] + energy[16];
  pEnergy->angle_ub = energy[11];
  pEnergy->imp = energy[12];
  pEnergy->cmap = energy[13];

  enmr[0] = energy[14];
  enmr[1] = energy[15];
  enmr[2] = energy[16];
  pEnergy->efield = energy[17];
//  printf("efield in md_energ %8.3f\n",energy[17]);

  pEnergy->total = pEnergy->vdw_dir + pEnergy->elec_dir + pEnergy->elec_recip +
    pEnergy->hbond + pEnergy->bond +
    pEnergy->angle + pEnergy->dihedral + pEnergy->vdw_14 + pEnergy->elec_14 +
    pEnergy->restraint + pEnergy->imp + pEnergy->angle_ub + pEnergy->efield +
    pEnergy->cmap;

  if (masterNode) pEnergy->total += (pEnergy->vdw_recip + pEnergy->elec_self);
//NEB
#ifdef MPI
  if ((*ineb > 0) && masterNode && (*nstep%*nebfreq ==0)) {
    MPI_Allgather(&pEnergy->total, 1, MPI_PMEDOUBLE, gpu->pbNEBEnergyAll->_pSysData,
                                   1, MPI_PMEDOUBLE, MPI_COMM_WORLD);
    gpu->pbNEBEnergyAll->Upload();
  }
#endif

  // Grab virial if needed
  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
    virial[0] = 0.5 * energy[VIRIAL_OFFSET + 0];
    virial[1] = 0.5 * energy[VIRIAL_OFFSET + 1];
    virial[2] = 0.5 * energy[VIRIAL_OFFSET + 2];
    if (masterNode) {
      virial[0] -= 0.5 * (gpu->ee_plasma + 2.0 * gpu->vdw_recip);
      virial[1] -= 0.5 * (gpu->ee_plasma + 2.0 * gpu->vdw_recip);
      virial[2] -= 0.5 * (gpu->ee_plasma + 2.0 * gpu->vdw_recip);
    }
    ekcmt[0] = energy[VIRIAL_OFFSET + 3];
    ekcmt[1] = energy[VIRIAL_OFFSET + 4];
    ekcmt[2] = energy[VIRIAL_OFFSET + 5];
  }
}

//---------------------------------------------------------------------------------------------
// gti_ele_recip_:
//
// Arguments:
//   needPotEnergy:
//   need_virials:
//   vol:
//   ti_weights:
//   ti_mode:
//---------------------------------------------------------------------------------------------
extern "C" void gti_ele_recip_(bool* needPotEnergy, bool* need_virials, double* vol, 
                               int* ti_mode, int* reaf_mode)
{
  PRINTMETHOD("gti_ele_recip");

  gpuContext gpu = theGPUContext::GetPointer(); 

  if (useMPI && !gpu->bCalculateReciprocalSum) {
    return;
  }
  if (*ti_mode == 0) {
    icc_CalculateElecRecipForceEnergy(gpu, *vol, *needPotEnergy, 0, 0, *reaf_mode);
  }
  else {
    for (int TIRegion = 0; TIRegion < 2; TIRegion++) {

      if (*reaf_mode < 0 || TIRegion == (*reaf_mode)) {
        icc_CalculateElecRecipForceEnergy(gpu, *vol, *needPotEnergy, *ti_mode, TIRegion, *reaf_mode);

        ik_CopyToTIForce(gpu, TIRegion, true, false, gpu->sim.TIItemWeight[Schedule::TypeEleRec][TIRegion]);
        // Virials are pre-weighted
        if (*need_virials) {
          ik_CopyToTIEnergy(gpu, TIRegion, VIRIAL_OFFSET, VIRIAL_OFFSET + 1, VIRIAL_OFFSET + 2, true, gpu->sim.TIItemWeight[Schedule::TypeEleRec][TIRegion], true);
        }
        // PME energy is not pre-weighted
        if (*needPotEnergy) {
          ik_CopyToTIEnergy(gpu, TIRegion, 9, -1, -1, false, 1.0);
        }
      }
    }
  }
}

//---------------------------------------------------------------------------------------------
// gti_nb_:
//
// Arguments:
//   needPotEnergy:  
//   needvirial:     
//   ti_mode:        
//   lj1264:  
//   plj1264 C4PairwiseCUDA       
//---------------------------------------------------------------------------------------------
extern "C" void gti_nb_(bool* needPotEnergy, bool* needvirial, int* ti_mode, int* lj1264, int* plj1264, int* reaf_mode)
{
  PRINTMETHOD("gti_nb");

  gpuContext gpu = theGPUContext::GetPointer();
  
  if (useMPI && !gpu->bCalculateDirectSum) {
    return;
  }

  if (*reaf_mode >= 0) {
    ik_ScaleRECharge(gpu, *reaf_mode);
  }
  if (*ti_mode > 0 && *reaf_mode<0 ) {
    ik_ZeroTICharge(gpu, 2); // Zero out all TI region charges
  }

  (*needPotEnergy) ? kCalculatePMENonbondEnergy(gpu) : kCalculatePMENonbondForces(gpu);
  
  if (*lj1264 > 0 && *plj1264 == 0) {
    ik_Calculate1264NB(gpu, *needPotEnergy, *needvirial);
  }
  if (*plj1264 > 0 && *lj1264 == 0) { // C4PairwiseCUDA
    ik_Calculatep1264NB(gpu, *needPotEnergy, *needvirial);
  }
  if (*plj1264 > 0 && *lj1264 > 0) { // C4PairwiseCUDA
    ik_Calculate1264p1264NB(gpu, *needPotEnergy, *needvirial);
  }
  if (*reaf_mode >= 0) {
    ik_CalculateREAFNb(gpu, *needPotEnergy, *needvirial);
    ik_CalculateREAF14NB(gpu, *needPotEnergy, *needvirial);
  }

  if (*ti_mode) {
    ik_CalculateTINB(gpu, *needPotEnergy, *needvirial);
    ik_CalculateTI14NB(gpu, *needPotEnergy, *needvirial);
  }

}

//---------------------------------------------------------------------------------------------
// gti_nb_gamd_:
//
// Arguments:
//   needPotEnergy:  
//   needvirial:     
//   ti_mode:        
//   lj1264:      
//   plj1264: C4PairwiseCUDA   
//---------------------------------------------------------------------------------------------
extern "C" void gti_nb_gamd_(bool* needPotEnergy, bool* needvirial, int* ti_mode, int* lj1264, int* plj1264)
{
  PRINTMETHOD("gti_nb_gamd");

  gpuContext gpu = theGPUContext::GetPointer();
  
  if (useMPI && !gpu->bCalculateDirectSum) {
    return;
  }
  if (*ti_mode > 0) {
    ik_ZeroTICharge(gpu, 2); // Zero out all TI region charges
  }
  (*needPotEnergy) ? kCalculatePMENonbondEnergy(gpu) : kCalculatePMENonbondForces(gpu);
  if (*lj1264 > 0 && *plj1264 == 0) {
    ik_Calculate1264NB(gpu, *needPotEnergy, *needvirial);
  }
  if (*plj1264 > 0 && *lj1264 == 0) { // C4PairwiseCUDA
    ik_Calculatep1264NB(gpu, *needPotEnergy, *needvirial);
  }
  if (*plj1264 > 0 && *lj1264 > 0) { // C4PairwiseCUDA
    ik_Calculate1264p1264NB(gpu, *needPotEnergy, *needvirial);
  }
  if (*ti_mode > 0) {
    ik_CalculateTINB_gamd(gpu, *needPotEnergy, *needvirial);
    // ik_CalculateTI14NB(gpu, *needPotEnergy, *needvirial);
    ik_CalculateTI14NB_gamd(gpu, *needPotEnergy, *needvirial);
  }
}
//---------------------------------------------------------------------------------------------
// gti_nb_ppi_gamd_:
//
// Arguments:
//   needPotEnergy:
//   needvirial:
//   ti_mode:
//   lj1264:
//   plj1264: C4PairwiseCUDA
//---------------------------------------------------------------------------------------------
extern "C" void gti_nb_ppi_gamd_(bool* needPotEnergy, bool* needvirial, int* ti_mode, int* lj1264, int* plj1264)
{
  PRINTMETHOD("gti_nb_gamd");

  gpuContext gpu = theGPUContext::GetPointer();

  if (useMPI && !gpu->bCalculateDirectSum) {
    return;
  }
  if (*ti_mode > 0) {
    ik_ZeroTICharge(gpu, 2); // Zero out all TI region charges
  }
  (*needPotEnergy) ? kCalculatePMENonbondEnergy(gpu) : kCalculatePMENonbondForces(gpu);
  if (*lj1264 > 0 && *plj1264 == 0) {
    ik_Calculate1264NB(gpu, *needPotEnergy, *needvirial);
  }
  if (*plj1264 > 0 && *lj1264 == 0) { // C4PairwiseCUDA
    ik_Calculatep1264NB(gpu, *needPotEnergy, *needvirial);
  }
  if (*plj1264 > 0 && *lj1264 > 0) { // C4PairwiseCUDA
    ik_Calculate1264p1264NB(gpu, *needPotEnergy, *needvirial);
  }
  if (*ti_mode > 0) {
    ik_CalculateTINB_ppi_gamd(gpu, *needPotEnergy, *needvirial);
//    ik_CalculateTI14NB_gamd(gpu, *needPotEnergy, *needvirial); !! first not use the 14NB
  }
  ik_ZeroTICharge(gpu, -1);//restore origin charges
  //update the energy from gti calculations
  for (int TIRegion = 0; TIRegion < 2; TIRegion++) {
      if (*needPotEnergy) {
        ik_CopyToTIEnergy(gpu, TIRegion, 9, -1, -1, 0, false);
      }
    }
}

//---------------------------------------------------------------------------------------------
// gti_nb_ppi2_gamd_:
//
// Arguments:
//   needPotEnergy:
//   needvirial:
//   ti_mode:
//   lj1264:
//   plj1264:
//-------------------------------------------------------------------------------
extern "C" void gti_nb_ppi2_gamd_(bool* needPotEnergy, bool* needvirial, int* ti_mode, int* lj1264, int* plj1264, int* bgpro2atm,int* edpro2atm)
{
  PRINTMETHOD("gti_nb_gamd");

  gpuContext gpu = theGPUContext::GetPointer();

  if (useMPI && !gpu->bCalculateDirectSum) {
    return;
  }
  if (*ti_mode > 0) {
    ik_ZeroTICharge(gpu, 2); // Zero out all TI region charges
  }
  (*needPotEnergy) ? kCalculatePMENonbondEnergy(gpu) : kCalculatePMENonbondForces(gpu);
  if (*lj1264 > 0 && *plj1264 == 0) {
    ik_Calculate1264NB(gpu, *needPotEnergy, *needvirial);
  }
  if (*plj1264 > 0 && *lj1264 == 0) { // C4PairwiseCUDA
    ik_Calculatep1264NB(gpu, *needPotEnergy, *needvirial);
  }
  if (*plj1264 > 0 && *lj1264 > 0) { // C4PairwiseCUDA
    ik_Calculate1264p1264NB(gpu, *needPotEnergy, *needvirial);
  }
  if (*ti_mode > 0) {
    ik_CalculateTINB_ppi2_gamd(gpu, *needPotEnergy, *needvirial, *bgpro2atm,*edpro2atm);
  }
  ik_ZeroTICharge(gpu, -1);//restore origin charges

  for (int TIRegion = 0; TIRegion < 2; TIRegion++) {
      if (*needPotEnergy) {
        ik_CopyToTIEnergy(gpu, TIRegion, 9, -1, -1, 0, false);
      }
    }
}
//---------------------------------------------------------------------------------------------
// gti_nb_ppi3_gamd_:
//
// Arguments:
//   needPotEnergy:
//   needvirial:
//   ti_mode:
//   lj1264:
//   plj1264: // C4PairwiseCUDA
//---------------------------------------------------------------------------------------------
extern "C" void gti_nb_ppi3_gamd_(bool* needPotEnergy, bool* needvirial, int* ti_mode, int* lj1264, int* plj1264)
{
  PRINTMETHOD("gti_nb_gamd");

  gpuContext gpu = theGPUContext::GetPointer();

  if (useMPI && !gpu->bCalculateDirectSum) {
    return;
  }
  if (*ti_mode > 0) {
    ik_ZeroTICharge(gpu, 2); // Zero out all TI region charges
  }
  (*needPotEnergy) ? kCalculatePMENonbondEnergy(gpu) : kCalculatePMENonbondForces(gpu);
  if (*lj1264 > 0 && *plj1264 == 0) {
    ik_Calculate1264NB(gpu, *needPotEnergy, *needvirial);
  }
  if (*plj1264 > 0 && *lj1264 == 0) { // C4PairwiseCUDA
    ik_Calculatep1264NB(gpu, *needPotEnergy, *needvirial);
  }
  if (*plj1264 > 0 && *lj1264 > 0) { // C4PairwiseCUDA
    ik_Calculate1264p1264NB(gpu, *needPotEnergy, *needvirial);
  }
  if (*ti_mode > 0) {
    ik_CalculateTINB_ppi_gamd(gpu, *needPotEnergy, *needvirial);
    ik_CalculateTI14NB_ppi_gamd(gpu, *needPotEnergy, *needvirial); // first not use the 14NB
  }
  ik_ZeroTICharge(gpu, -1);//restore origin charges
  //update the energy from gti calculations
  for (int TIRegion = 0; TIRegion < 2; TIRegion++) {
      if (*needPotEnergy) {
        ik_CopyToTIEnergy(gpu, TIRegion, 9, -1, -1, 0, false);
      }
    }
}



//---------------------------------------------------------------------------------------------
// gti_finalize_force_virial_:
//
// Arguments:
//   numextra:    
//   needvirial:  
//   ti_mode:     
//   net_force:   
//   frc:         
//---------------------------------------------------------------------------------------------
extern "C" void gti_finalize_force_virial_(int* numextra, bool* needvirial, int* ti_mode,
                                           int* net_force, double frc[][3])
{
  PRINTMETHOD(__func__);
  gpuContext gpu = theGPUContext::GetPointer();

  // make everyone is finished
  if (gpu->multiStream) cudaDeviceSynchronize();

  // Download forces on critical atoms for additional CPU calculations
  if (gpu->sim.nShuttle > 0) {
    kRetrieveSimData(gpu, frc, 2);
  }

  // Upload forces computed on the CPU to the device
  if (gpu->sim.nShuttle > 0) {
    kPostSimData(gpu, frc, 2);
  }
  if (*ti_mode > 0) {
    if (gpu->multiStream) cudaEventRecord(gpu->event_TIDone, gpu->TIStream);
    ik_CombineTIForce(gpu, true, *needvirial);
  }
  if (gpu->sim.EPs > 0) {
    kOrientForces(gpu);
  }
  if ((*net_force > 0) && (!useMPI || gpu->bCalculateReciprocalSum)) {
    ik_RemoveTINetForce(gpu, -1);
  }
  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
    if (!useMPI || gpu->gpuID == 0) {
      kCalculateCOMKineticEnergy(gpu);
      kReduceCOMKineticEnergy(gpu);
    }
  }

  if (*needvirial) {
    kCalculateMolecularVirial(gpu);
    if (*ti_mode > 0) {
      ik_CopyToTIEnergy(gpu, 2, VIRIAL_OFFSET, VIRIAL_OFFSET+1, VIRIAL_OFFSET + 2, true, 1.0, true);
    }
  }
  if (useMPI) gpu_allreduce(gpu->pbForceAccumulator, gpu->sim.stride3);
}

//---------------------------------------------------------------------------------------------
// gti_finalize_force_virial_gamd_:
//
// Arguments:
//   numextra:    
//   needvirial:  
//   ti_mode:     
//   net_force:   
//   frc:         
//---------------------------------------------------------------------------------------------
extern "C" void gti_finalize_force_virial_gamd_(int* numextra, bool* needvirial, int* ti_mode,
                                           int* net_force, double frc[][3])
{
  PRINTMETHOD("gti_finalize_force_virial_gamd");

  gpuContext gpu = theGPUContext::GetPointer();

  // Download forces on critical atoms for additional CPU calculations
  if (gpu->sim.nShuttle > 0) {
    kRetrieveSimData(gpu, frc, 2);
  }

  // Upload forces computed on the CPU to the device
  if (gpu->sim.nShuttle > 0) {
    kPostSimData(gpu, frc, 2);
  }
  if (*ti_mode > 0) {
    if (gpu->multiStream) cudaEventRecord(gpu->event_TIDone, gpu->TIStream);
    ik_CombineTIForce_gamd(gpu, true, *needvirial);
  }
  if (gpu->sim.EPs > 0) {
    kOrientForces(gpu);
  }
  if ((*net_force > 0) && (!useMPI || gpu->bCalculateReciprocalSum)) {
    ik_RemoveTINetForce(gpu, -1);
  }
  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
    if (!useMPI || gpu->gpuID == 0) {
      kCalculateCOMKineticEnergy(gpu);
      kReduceCOMKineticEnergy(gpu);
    }
  }

  if (*needvirial) {
    kCalculateMolecularVirial(gpu);
    if (*ti_mode > 0) {
      //TBF ik_CopyToTIEnergy(gpu, 2, VIRIAL_OFFSET, VIRIAL_OFFSET + 2, 1, false);
    }
  }
  if (useMPI) gpu_allreduce(gpu->pbForceAccumulator, gpu->sim.stride3);
}

//---------------------------------------------------------------------------------------------
// gti_finalize_force_gamd_:
//
// Arguments:
//   numextra:    
//   needvirial:  
//   ti_mode:     
//   net_force:   
//   frc:         
//---------------------------------------------------------------------------------------------
extern "C" void gti_finalize_force_gamd_(int* numextra, bool* needvirial, int* ti_mode,
                                           int* net_force, double frc[][3])
{
  PRINTMETHOD("gti_finalize_force_gamd");

  gpuContext gpu = theGPUContext::GetPointer();

  // Download forces on critical atoms for additional CPU calculations
  if (gpu->sim.nShuttle > 0) {
    kRetrieveSimData(gpu, frc, 2);
  }

  // Upload forces computed on the CPU to the device
  if (gpu->sim.nShuttle > 0) {
    kPostSimData(gpu, frc, 2);
  }
  if (*ti_mode > 0) {
    if (gpu->multiStream) cudaEventRecord(gpu->event_TIDone, gpu->TIStream);
    ik_CombineTIForce_gamd(gpu, true, *needvirial);
  }
  if (gpu->sim.EPs > 0) {
    kOrientForces(gpu);
  }
  if ((*net_force > 0) && (!useMPI || gpu->bCalculateReciprocalSum)) {
    ik_RemoveTINetForce(gpu, -1);
  }
  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
    if (!useMPI || gpu->gpuID == 0) {
      kCalculateCOMKineticEnergy(gpu);
      kReduceCOMKineticEnergy(gpu);
    }
  }

  if (useMPI) gpu_allreduce(gpu->pbForceAccumulator, gpu->sim.stride3);
}

//---------------------------------------------------------------------------------------------
// gti_virial_gamd_:
//
// Arguments:
//   numextra:    
//   needvirial:  
//   ti_mode:     
//   net_force:   
//   frc:         
//---------------------------------------------------------------------------------------------
extern "C" void gti_virial_gamd_(int* numextra, bool* needvirial, int* ti_mode,
                                           int* net_force, double frc[][3])
{
  PRINTMETHOD("gti_finalize_force_virial_gamd");

  gpuContext gpu = theGPUContext::GetPointer();

  if (*needvirial) {
    kCalculateMolecularVirial(gpu);
    if (*ti_mode > 0) {
      //TBF ik_CopyToTIEnergy(gpu, 2, VIRIAL_OFFSET, VIRIAL_OFFSET + 2, 1, false);
    }
  }
}

//---------------------------------------------------------------------------------------------
// gti_get_virial_:
//
// Arguments:
//   TIWeights:  
//   virial:     
//   ekcmt:      
//   ti_mode:    
//---------------------------------------------------------------------------------------------
extern "C" void gti_get_virial_(double TIWeights[2], double virial[3], double ekcmt[3],
                                int* ti_mode)
{
  PRINTMETHOD("gti_get_virial");

  gpuContext gpu = theGPUContext::GetPointer();
  
  const int EKCMT_OFFSET = VIRIAL_OFFSET + 3;
  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
    double TIPotEnergy[gti_simulationConst::TIEnergyBufferMultiplier][gti_simulationConst::GPUPotEnergyTerms];
    icc_GetTIPotEnergyFromGPU(gpu, TIPotEnergy);
    for (unsigned i = 0; i < 3; i++) {
      virial[i] = TIPotEnergy[2][VIRIAL_OFFSET + i] - gpu->ee_plasma - 2*gpu->vdw_recip;
      ekcmt[i] += TIPotEnergy[2][EKCMT_OFFSET + i];
      if (*ti_mode > 0) {
        virial[i] += TIPotEnergy[0][VIRIAL_OFFSET + i] +
          TIPotEnergy[1][VIRIAL_OFFSET + i] ;
        ekcmt[i] += (TIPotEnergy[0][EKCMT_OFFSET + i] * TIWeights[0]) +
                    (TIPotEnergy[1][EKCMT_OFFSET + i] * TIWeights[1]);
      }
      virial[i] *= 0.5;
    }
  }
}

//---------------------------------------------------------------------------------------------
// gti_get_virial_gamd_:
//
// Arguments:
//   TIWeights:  
//   virial:     
//   ekcmt:      
//   ti_mode:    
//---------------------------------------------------------------------------------------------
extern "C" void gti_get_virial_gamd_(double TIWeights[2], double virial[3], double ekcmt[3],
                                int* ti_mode)
{
  PRINTMETHOD("gti_get_virial_gamd");

  gpuContext gpu = theGPUContext::GetPointer();
  
  const int EKCMT_OFFSET = VIRIAL_OFFSET + 3;
  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
    double TIPotEnergy[gti_simulationConst::TIEnergyBufferMultiplier][gti_simulationConst::GPUPotEnergyTerms];
    icc_GetTIPotEnergyFromGPU(gpu, TIPotEnergy);
    for (unsigned i = 0; i < 3; i++) {
      virial[i] += 0.5 * TIPotEnergy[2][VIRIAL_OFFSET + i];
      ekcmt[i] += TIPotEnergy[2][EKCMT_OFFSET + i];
      if (*ti_mode > 0) {
        virial[i] += 0.5 * (TIPotEnergy[0][VIRIAL_OFFSET + i] * TIWeights[0]) +
                     0.5 * (TIPotEnergy[1][VIRIAL_OFFSET + i] * TIWeights[1]);
        ekcmt[i] += (TIPotEnergy[0][EKCMT_OFFSET + i] * TIWeights[0]) +
                    (TIPotEnergy[1][EKCMT_OFFSET + i] * TIWeights[1]);
      }
    }
  }
}

//---------------------------------------------------------------------------------------------
// gti_kinetic_:
//
// Arguments:
//   c_ave:  
//---------------------------------------------------------------------------------------------
extern "C" void gti_kinetic_(double* c_ave)
{
  PRINTMETHOD("gti_kinetic");

  gpuContext gpu = theGPUContext::GetPointer();

  ik_CalculateTIKineticEnergy(gpu, *c_ave);
}

//---------------------------------------------------------------------------------------------
// gti_get_pot_energy_:
//
// Arguments:
//   TIPotEnergy:
//   ti_mode:      
//---------------------------------------------------------------------------------------------
extern "C" void gti_get_pot_energy_(double
                                    TIPotEnergy[3][gti_simulationConst::GPUPotEnergyTerms],
                                    int* ti_mode)
{  
  PRINTMETHOD("gti_get_pot_energy");

  gpuContext gpu = theGPUContext::GetPointer();  
  //cudaStreamSynchronize(gpu->TIStream);
  icc_GetTIPotEnergyFromGPU(gpu, TIPotEnergy);
}

//---------------------------------------------------------------------------------------------
// gti_get_mbar_energy_:
//
// Arguments:
//   MBAREnergy:  
//---------------------------------------------------------------------------------------------
extern "C" void gti_get_mbar_energy_(double MBAREnergy[])
{
  PRINTMETHOD("gti_get_mbar_energy");

  gpuContext gpu = theGPUContext::GetPointer();
  //cudaStreamSynchronize(gpu->TIStream);
  icc_GetMBAREnergyFromGPU(gpu, MBAREnergy);
}

//---------------------------------------------------------------------------------------------
// gti_get_kin_energy_:
//
// Arguments:
//   TIKinEnergy:  
//---------------------------------------------------------------------------------------------
extern "C" void gti_get_kin_energy_(double
                                    TIKinEnergy[3][gti_simulationConst::GPUKinEnergyTerms])
{
  PRINTMETHOD("gti_get_kin_energy");

  gpuContext gpu = theGPUContext::GetPointer();
  cudaStreamSynchronize(gpu->TIStream);
  //cudaDeviceSynchronize();

  icc_GetTIKinEnergyFromGPU(gpu, TIKinEnergy);
}

//---------------------------------------------------------------------------------------------
// gti_download_frc_:
//
// Arguments:
//   pTIRegion:      
//   atm_frc:       
//   need_virials:  
//   fresh:         
//---------------------------------------------------------------------------------------------
extern "C" void gti_download_frc_(int* pTIRegion, double atm_frc[][3], bool* need_virials,
                                  bool* fresh)
{
  PRINTMETHOD(__func__);

  gpuContext gpu = theGPUContext::GetPointer();
  
  if (*pTIRegion < 0) {
    if (!gpu->pbForceAccumulator) return;
  } else {
    if (!gpu->pbTIForce) return;
  }

  if (*fresh){
    if (*pTIRegion < 0) {
      gpu->pbForceAccumulator->Download();
    } else {
      gpu->pbTIForce->Download();
    }
  }
  double scale = 1.0;
  scale = ONEOVERFORCESCALE;
  PMEAccumulator *pForce= (*pTIRegion >= 0) ? gpu->pbTIForce->_pSysData +
                                              gpu->sim.stride3*(*pTIRegion) :
                                              gpu->pbForceAccumulator->_pSysData;
  PMEAccumulator *pNBForce;
  if (*need_virials) {
    pNBForce = (*pTIRegion >= 0) ? pForce + (gpu->sim.stride3 * 5) : pForce + gpu->sim.stride3;
  }
  if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
    gpu->pbImageIndex->Download();
    unsigned int* pImageAtomLookup = &(gpu->pbImageIndex->_pSysData[gpu->sim.imageStride * 2]);
    for (int i = 0; i < gpu->sim.atoms; i++){
      int i1 = pImageAtomLookup[i];
      atm_frc[i][0] = (double)(pForce[i1])* scale;
      atm_frc[i][1] = (double)(pForce[i1 + gpu->sim.stride])* scale;
      atm_frc[i][2] = (double)(pForce[i1 + gpu->sim.stride2])* scale;
      if (*need_virials) {
        atm_frc[i][0] += (double)(pNBForce[i1])* scale;
        atm_frc[i][1] += (double)(pNBForce[i1 + gpu->sim.stride])* scale;
        atm_frc[i][2] += (double)(pNBForce[i1 + gpu->sim.stride2])* scale;
      }
    }
  }
  else {
    for (int i = 0; i < gpu->sim.atoms; i++) {
      atm_frc[i][0] = (double)(pForce[i])* scale;
      atm_frc[i][1] = (double)(pForce[i + gpu->sim.stride])* scale;
      atm_frc[i][2] = (double)(pForce[i + gpu->sim.stride2])* scale;
      if (*need_virials) {
        atm_frc[i][0] += (double)(pNBForce[i])* scale;
        atm_frc[i][1] += (double)(pNBForce[i + gpu->sim.stride])* scale;
        atm_frc[i][2] += (double)(pNBForce[i + gpu->sim.stride2])* scale;
      }
    }
  }
}

extern "C" void gti_reaf_lambda_schedule_(char* schFileName) {

  LambdaSchedule& schd = LambdaSchedule::GetReference();
  schFileName[256] = 0;
  schd.Init(schFileName);
}

//---------------------------------------------------------------------------------------------
// gti_update_lambda_:
//
// Arguments:
//   pLambda:  
//---------------------------------------------------------------------------------------------
extern "C" void gti_update_lambda_(double* pLambda, int* pIndex, int* pK, double ti_weights[2], double ti_dweights[2], double ti_item_weights[2][Schedule::TypeTotal], double ti_item_dweights[2][Schedule::TypeTotal], bool* pUse)
{
  PRINTMETHOD(__func__);
  
  gti_update_lambda_weights_(pLambda, pK, ti_weights, ti_dweights,
    ti_item_weights, ti_item_dweights, pUse);

  gpuContext gpu = theGPUContext::GetPointer();

  gpu->sim.TIk = *pK;
  if (*pIndex >=0) gpu->sim.currentMBARState = *pIndex;
  gpu->sim.TILambda[0] = *pLambda;
  gpu->sim.TILambda[1] = 1.0 - *pLambda;
  gpu->sim.TIWeight[0] = ti_weights[0];
  gpu->sim.TIWeight[1] = ti_weights[1];
  gpu->sim.TIdWeight[0] = ti_dweights[0];
  gpu->sim.TIdWeight[1] = ti_dweights[1];
  gpu->sim.TIPrefactor = (*pK) * pow((1.0 - *pLambda), *pK - 1);

  for (unsigned i = 0; i < Schedule::TypeTotal ; i++) {
    gpu->sim.TIItemWeight[i][0] = ti_item_weights[0][i];
    gpu->sim.TIItemWeight[i][1] = ti_item_weights[1][i];
    gpu->sim.TIItemdWeight[i][0] = ti_item_dweights[0][i];
    gpu->sim.TIItemdWeight[i][1] = ti_item_dweights[1][i];
  }

  LambdaSchedule& schd = LambdaSchedule::GetReference();
  gpu->sim.eleSmoothLambdaType = schd.GetFunctionType(Schedule::TypeEleSC);
  gpu->sim.vdwSmoothLambdaType = schd.GetFunctionType(Schedule::TypeVDW);
  gpu->sim.SCSmoothLambdaType = schd.GetFunctionType(Schedule::TypeEleSSC);

  gpu->UpdateSimulationConst();
}

extern "C" void gti_update_lambda_weights_(double* pLambda, int* pK, double ti_weights[2], double ti_dweights[2], double ti_item_weights[2][Schedule::TypeTotal], double ti_item_dweights[2][Schedule::TypeTotal], bool* pUseSP) {
  
  double w[2], dW[2];

  LambdaSchedule& schd = LambdaSchedule::GetReference();
  schd.Init(*pUseSP);

  // TL: This part is not right for k!=1, need to re-visit--need to add as a special case in LambdaSchedule
  if (*pUseSP) {

    for (LambdaSchedule::InteractionType type = Schedule::TypeGen; type < Schedule::TypeTotal; type = LambdaSchedule::InteractionType(type + 1)) {
      schd.GetWeight(type, *pLambda, w[0], dW[0], LambdaSchedule::forward);
      ti_item_weights[0][type] = w[0];
      ti_item_dweights[0][type] = dW[0];
      schd.GetWeight(type, *pLambda, w[1], dW[1], LambdaSchedule::backward);
      ti_item_weights[1][type] = w[1];
      ti_item_dweights[1][type] = dW[1];
    }
    
    ti_weights[0] = ti_item_weights[0][Schedule::TypeGen];
    ti_weights[1] = ti_item_weights[1][Schedule::TypeGen];
    ti_dweights[0] = ti_item_dweights[0][Schedule::TypeGen];
    ti_dweights[1] = ti_item_dweights[1][Schedule::TypeGen];

  } else {

    w[0] = pow((1.0 - *pLambda), *pK);
    w[1] = 1.0 - w[0];

    ti_weights[0] = w[0];
    ti_weights[1] = w[1];
    ti_dweights[0] = -1.0;
    ti_dweights[1] = 1.0;

    for (unsigned i = 0; i < Schedule::TypeTotal ; i++) {
      ti_item_weights[0][i] = ti_weights[0];
      ti_item_weights[1][i] = ti_weights[1];
      ti_item_dweights[0][i] = ti_dweights[0];
      ti_item_dweights[1][i] = ti_dweights[1];
    }
  }

}


//---------------------------------------------------------------------------------------------
// gti_get_lambda_weights_: Get the weights and the derivatives of lambda values (for the current 
// set lambda scheduling)
//
// Arguments:
//   *pType: type of the interaction
//   *pN size of the input/output arrays
//   lambdas: an arrary contains multiple global lambda values
//   weights[]: the weights
//   dWeights[]: the weight derivatives wrt lambda
//---------------------------------------------------------------------------------------------
extern "C" void gti_get_lambda_weights_(int* pType, int* pN, double lambdas[], double weights[], double dWeights[]) {

  LambdaSchedule& schd = LambdaSchedule::GetReference();
  Schedule::InteractionType type = (Schedule::InteractionType)*pType;
  schd.GetWeight(type, *pN, lambdas, weights, dWeights);
}


extern "C" void gti_reaf_tau_schedule_(char* schFileName) {

  TauSchedule& schd = TauSchedule::GetReference();
  schFileName[256] = 0;
  schd.Init(schFileName);
}

//---------------------------------------------------------------------------------------------
// gti_update_tau_:
//
// Arguments:
//   pLambda:  
//---------------------------------------------------------------------------------------------
extern "C" void gti_update_tau_(double* pTau, double reaf_weights[2], double read_item_weights[2][Schedule::TypeTotal], bool* pUse)
{
  PRINTMETHOD(__func__);

  gti_update_tau_weights_(pTau, reaf_weights, read_item_weights, pUse);

  gpuContext gpu = theGPUContext::GetPointer();

  gpu->sim.current_tau = *pTau;
  gpu->sim.REAFWeight[0] = reaf_weights[0];
  gpu->sim.REAFWeight[1] = reaf_weights[1];

  for (unsigned i = 0; i < Schedule::TypeTotal; i++) {
    gpu->sim.REAFItemWeight[i][0] = read_item_weights[0][i];
    gpu->sim.REAFItemWeight[i][1] = read_item_weights[1][i];
  }

  gpu->sim.REAFItemWeight[Schedule::TypeEleRec][1] = read_item_weights[0][Schedule::TypeEleRec]* read_item_weights[0][Schedule::TypeEleRec];

  TauSchedule& schd = TauSchedule::GetReference();
  gpu->sim.eleTauType = schd.GetFunctionType(Schedule::TypeEleSC);
  gpu->sim.vdwTauType = schd.GetFunctionType(Schedule::TypeVDW);

  gpu->UpdateSimulationConst();
}

extern "C" void gti_update_tau_weights_(double* pTau, double reaf_weights[2], double read_item_weights[2][Schedule::TypeTotal], bool* pUseSP) {

  PRINTMETHOD(__func__);

  TauSchedule& schd = TauSchedule::GetReference();
  schd.Init(*pUseSP);

  if (*pUseSP) {

    double w, dW;
    for (Schedule::InteractionType type = Schedule::TypeGen; type < Schedule::TypeTotal; type = Schedule::InteractionType(type + 1)) {
      schd.GetWeight(type, *pTau, w, dW, Schedule::external);
      read_item_weights[0][type] = w;
      schd.GetWeight(type, *pTau, w, dW, Schedule::internal);
      read_item_weights[1][type] = w;
    }
    reaf_weights[0] = read_item_weights[0][Schedule::TypeGen];
    reaf_weights[1] = read_item_weights[1][Schedule::TypeGen];

  } else {

    reaf_weights[0] = OneF - *pTau;
    reaf_weights[1] = reaf_weights[0] * reaf_weights[0];

    for (Schedule::InteractionType type = Schedule::TypeGen; type < Schedule::TypeTotal; type = Schedule::InteractionType(type + 1)) {
      read_item_weights[0][type] = reaf_weights[0];
      read_item_weights[1][type] = reaf_weights[1];
    }
  }
}


//---------------------------------------------------------------------------------------------
// gti_setup_RMSD:
//---------------------------------------------------------------------------------------------
extern "C" void gti_init_rmsd_(int* pNumber_rmsd_set, int* pRmsd_type, int rmsd_atom_count[],
  int rmsd_ti_region[], int rmsd_atom_list[], double rmsd_ref_crd[], double rmsd_ref_com[],
  double rmsd_weights[])
{
  PRINTMETHOD(__func__);

  gpuContext gpu = theGPUContext::GetPointer();
  gpu->sim.rmsd_type = *pRmsd_type;

  for (unsigned i = 0; i < 5; i++) {
    gpu->sim.rmsd_atom_count[i] = 0;
    gpu->sim.rmsd_ti_region[i] = -1;
    gpu->sim.rmsd_weights[i] = 0.0;
  }

  unsigned n= *pNumber_rmsd_set;
  if (n > gti_simulationConst::MaxNumberRMSDRegion) {
    std::string e("Too many RMSD Regions--current limit: 5 regions");
    throw std::runtime_error(e);
  }
  unsigned max = 0;
  for (unsigned i = 0; i < n; i++) 
    if (rmsd_atom_count[i] > max) max = rmsd_atom_count[i];
  gpu->sim.rmsd_atom_max_count = max;

  if (max*n > gti_simulationConst::MaxNumberRMSDAtom * gti_simulationConst::MaxNumberRMSDRegion) {
    std::string e("Too many RMSD atoms--current limit: total 1000 atoms");
    throw std::runtime_error(e);
  }

  gpu->sim.number_rmsd_set = n;
  for (unsigned i = 0; i < n; i++) {
    gpu->sim.rmsd_atom_count[i] = rmsd_atom_count[i];
    gpu->sim.rmsd_ti_region[i] = rmsd_ti_region[i];
    gpu->sim.rmsd_weights[i] = rmsd_weights[i];
    for (unsigned j = 0; j < rmsd_atom_count[i]; j++) {
      unsigned k = i * max + j;
      gpu->sim.rmsd_atom_list[k] = rmsd_atom_list[k];
      for (unsigned l=0; l<3; l++)
        gpu->sim.rmsd_ref_crd[k*3+l] = rmsd_ref_crd[k*3+l];
    }
    for (unsigned l = 0; l < 3; l++)
      gpu->sim.rmsd_ref_com[i* 3 + l] = rmsd_ref_com[i * 3 + l];
  }
 
  gpu->pbMatrix = std::unique_ptr < GpuBuffer<PMEFloat>> (new GpuBuffer<PMEFloat>(5*16) );
  gpu->pbResult = std::unique_ptr < GpuBuffer<PMEFloat>>(new GpuBuffer<PMEFloat>(5*4));

  unsigned workSize = 65536;
  gpu->pbWork = std::unique_ptr < GpuBuffer<PMEFloat>> (new GpuBuffer<PMEFloat>(workSize) );
  gpu->pbInfo = std::unique_ptr < GpuBuffer<int>>(new GpuBuffer<int>(1));
  gpu->sim.pMatrix = gpu->pbMatrix->_pDevData;
  gpu->sim.pResult = gpu->pbResult->_pDevData;
  gpu->sim.pWork = gpu->pbWork->_pDevData;
  gpu->sim.pInfo = gpu->pbInfo->_pDevData;
  gpu->UpdateSimulationConst();

}

//---------------------------------------------------------------------------------------------
// gti_setup_localheating_:
//---------------------------------------------------------------------------------------------
extern "C" void gti_setup_localheating_(double* pTI_tempi) {

  gpuContext gpu = theGPUContext::GetPointer();
  gpu->sim.TIHeatingTemp = *pTI_tempi;
  gpu->sim.doTIHeating = true;

  gpu->UpdateSimulationConst();
}

//---------------------------------------------------------------------------------------------
// gti_update_localheating_:
//---------------------------------------------------------------------------------------------
extern "C" void gti_turnoff_localheating_() {

  gpuContext gpu = theGPUContext::GetPointer();
  gpu->sim.doTIHeating = false;

  gpu->UpdateSimulationConst();
}

//---------------------------------------------------------------------------------------------
// gti_update_simulation_const_:
//---------------------------------------------------------------------------------------------
extern "C" void gti_update_simulation_const_()
{
  PRINTMETHOD("gti_update_simulation_const");

  gpuContext gpu = theGPUContext::GetPointer();

  gpu->UpdateSimulationConst();
}

//---------------------------------------------------------------------------------------------
// gti_bonded_gamd_setup_:
//---------------------------------------------------------------------------------------------
extern "C" void gti_bonded_gamd_setup_(int* pBat_type, int* pError)
{
  gpuContext gpu = theGPUContext::GetPointer();
  unsigned aa[4];
  unsigned currentTIRegion;
  bool inTI = false;

  // Bonds
  if (true) {
    unsigned numberBond = gpu->sim.bonds;
    std::vector<PMEDouble2> tempBond; tempBond.clear();
    std::vector<uint2> tempID; tempID.clear();
    std::vector<uint> tempType; tempType.clear();

    for (unsigned int i = 0; i < numberBond; i++) {
      aa[0] = gpu->pbBondID->_pSysData[i].x;
      aa[1] = gpu->pbBondID->_pSysData[i].y;
      inTI = false;
      for (unsigned int j = 0; j < 2; j++) {
        bool b0 = (gpu->pbTIList->_pSysData[aa[j]])>0;
        bool b1 = (gpu->pbTIList->_pSysData[aa[j] + gpu->sim.stride]) > 0;
        if (b0 || b1) {
          inTI = true;
          currentTIRegion = (b0) ? 0 : 1;
          break;
        }
      }
      if (inTI) {
        PMEDouble2 tb = gpu->pbBond->_pSysData[i];
        tempBond.push_back(tb);
        uint2 tbid;
        tbid.x = aa[0];
        tbid.y = aa[1];
        tempID.push_back(tbid);
        tempType.push_back(currentTIRegion);

        // Zero out theforce constant
//        gpu->pbBond->_pSysData[i].x = 0.0;
      }
    }
    gpu->sim.numberTIBond = tempBond.size();
    uint regionCounter[2] = { 0, 0 };
    if (*pBat_type == 2) {
      for (unsigned int i = 0; i < gpu->sim.numberTIBond; i++) {
        aa[0] = tempID[i].x;
        aa[1] = tempID[i].y;
        uint currentTIRegion = (tempType[i] & 1);
        uint& type = tempType[i];
        if ( ((gpu->pbSCBATList->_pSysData[aa[0]] & 2) && (gpu->pbSCBATList->_pSysData[aa[1]] & 2) && currentTIRegion == 0) ||
          ((gpu->pbSCBATList->_pSysData[aa[0]] & 16) && (gpu->pbSCBATList->_pSysData[aa[1]] & 16) && currentTIRegion == 1)  ){
          type |= gti_simulationConst::ex_userDefined_bat;
          regionCounter[currentTIRegion]++;
        }
      }
    }
    for (unsigned int i = 0; i < gpu->sim.numberTIBond; i++) {
      aa[0] = tempID[i].x;
      aa[1] = tempID[i].y;
      uint currentTIRegion = (tempType[i] & 1);
      uint& type = tempType[i];

      bool hasSC = (gpu->pbSCList->_pSysData[aa[0]] > 0 || gpu->pbSCList->_pSysData[aa[1]] > 0);
      bool inSC = (gpu->pbSCList->_pSysData[aa[0]] > 0 && gpu->pbSCList->_pSysData[aa[1]] > 0);
      if (hasSC) type |= gti_simulationConst::Fg_has_SC;
      if (inSC) type |= gti_simulationConst::Fg_int_SC;
      if (!hasSC) type |= gti_simulationConst::ex_addToDVDL_bat;

      if (hasSC && !inSC) {
        if (abs(*pBat_type) == 1) {
          if (hasSC && !inSC) {
            if (regionCounter[currentTIRegion] == 0) {
              regionCounter[currentTIRegion]++;  // pick up the first CC-SC bond
            } else {
              if (!(type & gti_simulationConst::ex_userDefined_bat)) {
                type |= gti_simulationConst::ex_addToBat_corr;
                if (*pBat_type > 0) type |= gti_simulationConst::ex_addToDVDL_bat;
              }
            }
          }
        } else if (abs(*pBat_type == 2) && (! (type & gti_simulationConst::ex_userDefined_bat)) ) {
          type |= gti_simulationConst::ex_addToBat_corr;
          if (*pBat_type > 0) type |= gti_simulationConst::ex_addToDVDL_bat;
        }
      }
    }

    gpu->pbTIBond = std::unique_ptr< GpuBuffer<PMEDouble2> > (new GpuBuffer<PMEDouble2>(gpu->sim.numberTIBond));
    gpu->pbTIBondID = std::unique_ptr< GpuBuffer<uint2> > (new GpuBuffer<uint2>(gpu->sim.numberTIBond));
    gpu->pbTIBondType = std::unique_ptr< GpuBuffer<uint> > (new GpuBuffer<uint>(gpu->sim.numberTIBond));
    for (unsigned int i = 0; i < gpu->sim.numberTIBond; i++) {
      gpu->pbTIBond->_pSysData[i] = tempBond[i];
      gpu->pbTIBondID->_pSysData[i] = tempID[i];
      gpu->pbTIBondType->_pSysData[i] = tempType[i];
    }
    gpu->pbTIBond->Upload();
    gpu->pbTIBondID->Upload();
    gpu->pbTIBondType->Upload();
    gpu->sim.pTIBond = gpu->pbTIBond->_pDevData;
    gpu->sim.pTIBondID = gpu->pbTIBondID->_pDevData;
    gpu->sim.pTIBondType = gpu->pbTIBondType->_pDevData;
  }

  // Angles
  if (true) {
    int numberAngle = gpu->sim.bondAngles;
    std::vector<PMEDouble2> tempAngle; tempAngle.clear();
    std::vector<uint4> tempID; tempID.clear();
    std::vector<uint> tempType; tempType.clear();

    for (unsigned int i = 0; i < numberAngle; i++) {
      aa[0] = gpu->pbBondAngleID1->_pSysData[i].x;
      aa[1] = gpu->pbBondAngleID1->_pSysData[i].y;
      aa[2] = gpu->pbBondAngleID2->_pSysData[i];
      inTI = false;
      for (unsigned int j = 0; j < 3; j++) {
        bool b0 = (gpu->pbTIList->_pSysData[aa[j]])>0;
        bool b1 = (gpu->pbTIList->_pSysData[aa[j] + gpu->sim.stride]) > 0;
        if (b0 || b1) {
          inTI = true;
          currentTIRegion = (b0) ? 0 : 1;
          break;
        }
      }
      if (inTI) {
        PMEDouble2 ta = gpu->pbBondAngle->_pSysData[i];
        tempAngle.push_back(ta);
        uint4 tbid;
        tbid.x = aa[0];
        tbid.y = aa[1];
        tbid.z = aa[2];
        tempType.push_back(currentTIRegion);
        tempID.push_back(tbid);

        // Zero out theforce constant
//        gpu->pbBondAngle->_pSysData[i].x =0.0;
      }
    }

    gpu->sim.numberTIAngle = tempAngle.size();
    uint regionCounter[2] = { 0, 0 };
    if (*pBat_type == 2) {
      for (unsigned int i = 0; i < gpu->sim.numberTIAngle; i++) {
        aa[0] = tempID[i].x;
        aa[1] = tempID[i].y;
        aa[2] = tempID[i].z;
        uint currentTIRegion = (tempType[i] & 1);
        uint& type = tempType[i];
        if (((gpu->pbSCBATList->_pSysData[aa[0]] & 4) && (gpu->pbSCBATList->_pSysData[aa[1]] & 4) && (gpu->pbSCBATList->_pSysData[aa[2]] & 4) && currentTIRegion == 0) ||
          ((gpu->pbSCBATList->_pSysData[aa[0]] & 32) && (gpu->pbSCBATList->_pSysData[aa[1]] & 32) && (gpu->pbSCBATList->_pSysData[aa[2]] & 32) && currentTIRegion == 1)) {
          type |= gti_simulationConst::ex_userDefined_bat;
          regionCounter[currentTIRegion]++;
        }
      }
    }
    for (unsigned int i = 0; i < gpu->sim.numberTIAngle; i++) {
      aa[0] = tempID[i].x;
      aa[1] = tempID[i].y;
      aa[2] = tempID[i].z;
      uint currentTIRegion = (tempType[i] & 1);
      uint& type = tempType[i];

      bool b0 = gpu->pbSCList->_pSysData[aa[0]] > 0;
      bool b1 = gpu->pbSCList->_pSysData[aa[1]] > 0;
      bool b2 = gpu->pbSCList->_pSysData[aa[2]] > 0;
      bool hasSC = (b0 || b1 || b2);
      bool inSC = (b0 && b1 && b2);
      if (hasSC) type |= gti_simulationConst::Fg_has_SC;
      if (inSC) type |= gti_simulationConst::Fg_int_SC;
      if (!hasSC) type |= gti_simulationConst::ex_addToDVDL_bat;
      if (hasSC && !inSC) {
        if (abs(*pBat_type) == 1) {
          if (hasSC && !inSC) {
            // RRD/DRR case
            if ((!b0 && !b1 && b2) || (b0 && !b1 && !b2)) {
              if (regionCounter[currentTIRegion] == 0) {
                regionCounter[currentTIRegion]++;
              } else {
                type |= gti_simulationConst::ex_addToBat_corr;
                if (*pBat_type > 0) type |= gti_simulationConst::ex_addToDVDL_bat;
              }
            } else {
              bool RDD = ((!b0 && !b1 && b2) || (b0 && !b1 && !b2));
              if (!RDD) {
                type |= gti_simulationConst::ex_addToBat_corr;
                if (*pBat_type > 0) type |= gti_simulationConst::ex_addToDVDL_bat;
              }
            }
          }
        } else if (abs(*pBat_type) == 2 && (!(type & gti_simulationConst::ex_userDefined_bat))) {
          type |= gti_simulationConst::ex_addToBat_corr;
          if (*pBat_type > 0) type |= gti_simulationConst::ex_addToDVDL_bat;
        }
      }
    }
    gpu->pbTIBondAngle = std::unique_ptr< GpuBuffer<PMEDouble2> > (new GpuBuffer<PMEDouble2>(gpu->sim.numberTIAngle));
    gpu->pbTIBondAngleID = std::unique_ptr< GpuBuffer<uint4> > (new GpuBuffer<uint4>(gpu->sim.numberTIAngle));
    gpu->pbTIBondAngleType = std::unique_ptr< GpuBuffer<uint> > (new GpuBuffer<uint>(gpu->sim.numberTIAngle));
    for (unsigned int i = 0; i < gpu->sim.numberTIAngle; i++) {
      gpu->pbTIBondAngle->_pSysData[i] = tempAngle[i];
      gpu->pbTIBondAngleID->_pSysData[i] = tempID[i];
      gpu->pbTIBondAngleType->_pSysData[i] = tempType[i];
    }
    gpu->pbTIBondAngle->Upload();
    gpu->pbTIBondAngleID->Upload();
    gpu->pbTIBondAngleType->Upload();
    gpu->sim.pTIBondAngle = gpu->pbTIBondAngle->_pDevData;
    gpu->sim.pTIBondAngleID = gpu->pbTIBondAngleID->_pDevData;
    gpu->sim.pTIBondAngleType = gpu->pbTIBondAngleType->_pDevData;
  }

  // Dihedrals
  if (true) {
    int numberDihedral = gpu->sim.dihedrals;
    std::vector<PMEDouble2> tempDihedral1; tempDihedral1.clear();
    std::vector<PMEDouble2> tempDihedral2; tempDihedral2.clear();
    std::vector<PMEDouble> tempDihedral3; tempDihedral3.clear();
    std::vector<uint4> tempID; tempID.clear();
    std::vector<uint> tempDihedralRegion; tempDihedralRegion.clear();
    std::vector<uint> tempType; tempType.clear();

    for (unsigned int i = 0; i < numberDihedral; i++) {
      aa[0] = gpu->pbDihedralID1->_pSysData[i].x;
      aa[1] = gpu->pbDihedralID1->_pSysData[i].y;
      aa[2] = gpu->pbDihedralID1->_pSysData[i].z;
      aa[3] = gpu->pbDihedralID1->_pSysData[i].w;
      inTI = false;
      for (unsigned int j = 0; j < 4; j++) {
        bool b0 = (gpu->pbTIList->_pSysData[aa[j]])>0;
        bool b1 = (gpu->pbTIList->_pSysData[aa[j] + gpu->sim.stride]) > 0;
        if (b0 || b1) {
          inTI = true;
          currentTIRegion = (b0) ? 0 : 1;
          break;
        }
      }
      if (inTI) {
        PMEDouble2 td1 = {gpu->pbDihedral1->_pSysData[i].x, gpu->pbDihedral1->_pSysData[i].y};
        tempDihedral1.push_back(td1);
        PMEDouble2 td2 = {gpu->pbDihedral2->_pSysData[i].x, gpu->pbDihedral2->_pSysData[i].y};
        tempDihedral2.push_back(td2);
        PMEDouble td3 = gpu->pbDihedral3->_pSysData[i];
        tempDihedral3.push_back(td3);
        uint4 tt = {aa[0], aa[1], aa[2], aa[3]};
        tempID.push_back(tt);
        tempType.push_back(currentTIRegion);
      }
    }

    gpu->sim.numberTIDihedral = tempDihedral1.size();
    uint regionCounter[2] = { 0, 0 };
    if (*pBat_type == 2) {
      for (unsigned int i = 0; i < gpu->sim.numberTIDihedral; i++) {
        aa[0] = tempID[i].x;
        aa[1] = tempID[i].y;
        aa[2] = tempID[i].z;
        aa[3] = tempID[i].w;
        uint currentTIRegion = (tempType[i] & 1);
        uint& type = tempType[i];
        if (((gpu->pbSCBATList->_pSysData[aa[0]] & 8) && (gpu->pbSCBATList->_pSysData[aa[1]] & 8)
          && (gpu->pbSCBATList->_pSysData[aa[2]] & 8) && (gpu->pbSCBATList->_pSysData[aa[3]] & 8)  && currentTIRegion == 0) ||
            ((gpu->pbSCBATList->_pSysData[aa[0]] & 64) && (gpu->pbSCBATList->_pSysData[aa[1]] & 64)
          && (gpu->pbSCBATList->_pSysData[aa[2]] & 64) && (gpu->pbSCBATList->_pSysData[aa[3]] & 64) && currentTIRegion == 1)) {
          type |= gti_simulationConst::ex_userDefined_bat;
          regionCounter[currentTIRegion]++;
        }
      }
    }
    for (unsigned int i = 0; i < gpu->sim.numberTIDihedral; i++) {
      aa[0] = tempID[i].x;
      aa[1] = tempID[i].y;
      aa[2] = tempID[i].z;
      aa[3] = tempID[i].w;
      uint currentTIRegion = (tempType[i] & 1);
      uint& type = tempType[i];

      bool b0 = gpu->pbSCList->_pSysData[aa[0]] > 0;
      bool b1 = gpu->pbSCList->_pSysData[aa[1]] > 0;
      bool b2 = gpu->pbSCList->_pSysData[aa[2]] > 0;
      bool b3 = gpu->pbSCList->_pSysData[aa[3]] > 0;
      bool hasSC = (b0 || b1 || b2 || b3);
      bool inSC = (b0 && b1 && b2 && b3);
      if (hasSC) type |= gti_simulationConst::Fg_has_SC;
      if (inSC) type |= gti_simulationConst::Fg_int_SC;
      if (!hasSC) type |= gti_simulationConst::ex_addToDVDL_bat;
      if (hasSC && !inSC) {
        if (abs(*pBat_type) == 1) {
          if (hasSC && !inSC) {
            // RRDD/DDRR case
            if ((!b0 && !b1 && b2 && b3) || (b0 && b1 && !b2 && !b3)) {
              if (regionCounter[currentTIRegion] == 0) {
                regionCounter[currentTIRegion]++;
              } else {
                type |= gti_simulationConst::ex_addToBat_corr;
                if (*pBat_type > 0) type |= gti_simulationConst::ex_addToDVDL_bat;
              }
            } else {
              // RDDD/DDDR case
              bool RDDD = ((!b0 && b1 && b2 && b3) || (b0 && b1 && b2 && b3));
              if (!RDDD) {
                type |= gti_simulationConst::ex_addToBat_corr;
                if (*pBat_type > 0) type |= gti_simulationConst::ex_addToDVDL_bat;
              }
            }
          }
        } else if (abs(*pBat_type) == 2 && (!(type & gti_simulationConst::ex_userDefined_bat))) {
          type |= gti_simulationConst::ex_addToBat_corr;
          if (*pBat_type > 0) type |= gti_simulationConst::ex_addToDVDL_bat;
        }
      }
    }

    gpu->pbTIDihedral1 = std::unique_ptr< GpuBuffer<PMEDouble2> > (new GpuBuffer<PMEDouble2>(gpu->sim.numberTIDihedral));
    gpu->pbTIDihedral2 = std::unique_ptr< GpuBuffer<PMEDouble2> > (new GpuBuffer<PMEDouble2>(gpu->sim.numberTIDihedral));
    gpu->pbTIDihedral3 = std::unique_ptr< GpuBuffer<PMEDouble> > (new GpuBuffer<PMEDouble>(gpu->sim.numberTIDihedral));
    gpu->pbTIDihedralID = std::unique_ptr< GpuBuffer<uint4> > (new GpuBuffer<uint4>(gpu->sim.numberTIDihedral));
    gpu->pbTIDihedralType = std::unique_ptr< GpuBuffer<uint> > (new GpuBuffer<uint>(gpu->sim.numberTIDihedral));
    for (unsigned int i = 0; i < gpu->sim.numberTIDihedral; i++) {
      gpu->pbTIDihedral1->_pSysData[i] = tempDihedral1[i];
      gpu->pbTIDihedral2->_pSysData[i] = tempDihedral2[i];
      gpu->pbTIDihedral3->_pSysData[i] = tempDihedral3[i];
      gpu->pbTIDihedralID->_pSysData[i] = tempID[i];
      gpu->pbTIDihedralType->_pSysData[i] = tempType[i];
    }
    gpu->pbTIDihedral1->Upload();
    gpu->pbTIDihedral2->Upload();
    gpu->pbTIDihedral3->Upload();
    gpu->pbTIDihedralID->Upload();
    gpu->pbTIDihedralType->Upload();
    gpu->sim.pTIDihedral1 = gpu->pbTIDihedral1->_pDevData;
    gpu->sim.pTIDihedral2 = gpu->pbTIDihedral2->_pDevData;
    gpu->sim.pTIDihedral3 = gpu->pbTIDihedral3->_pDevData;
    gpu->sim.pTIDihedralID = gpu->pbTIDihedralID->_pDevData;
    gpu->sim.pTIDihedralType = gpu->pbTIDihedralType->_pDevData;
  }
}
#endif /* GTI */
