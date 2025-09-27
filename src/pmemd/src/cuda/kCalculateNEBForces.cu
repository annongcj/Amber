#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------
#ifndef AMBER_PLATFORM_AMD
#include <cuda.h>
#endif

#include "gpu.h"
#include "ptxmacros.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

// Use global instance instead of a local copy
#include "simulationConst.h"
CSIM_STO simulationConst cSim;

bool first = true;

#if !defined(__HIPCC_RDC__)

//---------------------------------------------------------------------------------------------
// SetkCalculateNEBForcesSim: called by gpuCopyConstants to orient the instance of cSim in this
//                     library
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void SetkCalculateNEBForcesSim(gpuContext gpu)
{
  cudaError_t status;
  status = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(simulationConst));
  RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

//---------------------------------------------------------------------------------------------
// GetkCalculateNEBForcesSim: download information about the CUDA simulation data struct from the
//                     device.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//
// This is a debugging function.
//---------------------------------------------------------------------------------------------
void GetkCalculateNEBForcesSim(gpuContext gpu)
{
  cudaError_t status;
  status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(simulationConst));
  RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cSim failed");
}

#endif

//----------------------------------------------------------------------------------------------
// Shuffle Warp Reduce: https://devblogs.nvidia.com/faster-parallel-reductions-kepler/
//----------------------------------------------------------------------------------------------
__inline__ __device__ double warpReduceSum(double val) {
  for (int offset = GRID/2; offset > 0; offset /= 2)
    val += __shfl_down(val, offset);
  return val;
}

//----------------------------------------------------------------------------------------------
// atomicAdd function for double data types
//----------------------------------------------------------------------------------------------
//#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600
#if (__CUDA_ARCH__ < 600) || (__CUDACC_VER_MAJOR__ < 8)
__device__ void voidatomicAdd(double* ptr, double val) {
  unsigned long long int* ulli_ptr = (unsigned long long int*)ptr;
  unsigned long long int old = *ulli_ptr;
  unsigned long long int assumed;

  do {
    assumed = old;
    old = atomicCAS(ulli_ptr, assumed,  __double_as_longlong(val + __longlong_as_double(assumed)));
  } while (old != assumed);
}
#endif

//---------------------------------------------------------------------------------------------
// dgeev from lapack for computing eigenvalues and left/right enigenvectors of a square matrix
//---------------------------------------------------------------------------------------------
extern "C" {
extern int dgeev_(char*,char*,int*,double*,int*,double*, double*, double*, int*, double*, int*, double*, int*, int*);
}

//---------------------------------------------------------------------------------------------
// Kernel for clearing the NEB buffers at each step
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// Kernel for shuttling the coordinates for sending and receiving between the neighboring
// replicas in NEB.
//---------------------------------------------------------------------------------------------

#define TAKE_IMAGE
#define PULL_COORD
__global__ void kRetrieveCoordNL_kernel()
#include "kShuttle.h"
#undef PULL_COORD

#define PULL_FORCE
__global__ void kRetrieveForceNL_kernel()
#include "kShuttle.h"
#undef PULL_FORCE
#undef TAKE_IMAGE

#define PULL_COORD
__global__ void kRetrieveCoord_kernel()
#include "kShuttle.h"
#undef PULL_COORD

#define PULL_FORCE
__global__ void kRetrieveForce_kernel()
#include "kShuttle.h"
#undef PULL_FORCE

//----------------------------------------------------------------------------------------------
// kernels for RMSfit calculations.
//----------------------------------------------------------------------------------------------

#define NEIGHBOR_LIST
#define MASS_CALC
__global__ void kTotMassFitNL_kernel()
#include "kRMSFit.h"
#undef MASS_CALC

#define CALC_COM
__global__ void kCALCCOMNL_kernel()
#include "kRMSFit.h"
#undef CALC_COM

#define KABSCH
#define NEXT_NEIGHBOR
__global__ void kKabschNextNL_kernel()
#include "kRMSFit.h"
#undef NEXT_NEIGHBOR

#define PREV_NEIGHBOR
__global__ void kKabschPrevNL_kernel()
#include "kRMSFit.h"
#undef PREV_NEIGHBOR
#undef KABSCH
#undef NEIGHBOR_LIST

#define MASS_CALC
__global__ void kTotMassFit_kernel()
#include "kRMSFit.h"
#undef MASS_CALC

#define CALC_COM
__global__ void kCALCCOM_kernel()
#include "kRMSFit.h"
#undef CALC_COM

#define CENTER_DATA
__global__ void kCenterData_kernel()
#include "kRMSFit.h"
#undef CENTER_DATA

#define RECENTER_DATA
__global__ void kReCenterData_kernel()
#include "kRMSFit.h"
#undef RECENTER_DATA

#define ROTATE_NEIGHBOR
#define NEXT_NEIGHBOR
__global__ void kRotateNext_kernel()
#include "kRMSFit.h"
#undef NEXT_NEIGHBOR

#define PREV_NEIGHBOR
__global__ void kRotatePrev_kernel()
#include "kRMSFit.h"
#undef PREV_NEIGHBOR
#undef ROTATE_NEIGHBOR

#define KABSCH
#define NEXT_NEIGHBOR
__global__ void kKabschNext_kernel()
#include "kRMSFit.h"
#undef NEXT_NEIGHBOR

#define PREV_NEIGHBOR
__global__ void kKabschPrev_kernel()
#include "kRMSFit.h"
#undef PREV_NEIGHBOR
#undef KABSCH

//----------------------------------------------------------------------------------------------
// Kernel for calculating NEB forces.
//----------------------------------------------------------------------------------------------

#define REVISED_TAN
#define FIRST
__global__ void kNEBRevisedTanFIRST_kernel()
#include "kNEB.h"
#undef FIRST

#define SECOND
__global__ void kNEBRevisedTanSECOND_kernel()
#include "kNEB.h"
#undef SECOND

#define THIRD
__global__ void kNEBRevisedTanTHIRD_kernel()
#include "kNEB.h"
#undef THIRD

#define FOURTH
__global__ void kNEBRevisedTanFOURTH_kernel()
#include "kNEB.h"
#undef FOURTH
#undef REVISED_TAN

#define BASIC_TAN
__global__ void kNEBBasicTan_kernel()
#include "kNEB.h"
#undef BASIC_TAN

#define NORMALIZE
__global__ void kNormalizeTan_kernel()
#include "kNEB.h"
#undef NORMALIZE

#define NEIGHBOR_LIST
#define NEB_FRC
__global__ void kNEBFRCNL_kernel()
#include "kNEB.h"
#undef NEB_FRC
#undef NEIGHBOR_LIST

#define NEB_FRC
__global__ void kNEBFRC_kernel()
#include "kNEB.h"
#undef NEB_FRC

#define REDUCE_TAN
__global__ void kNormalize_kernel()
#include "kNorm.h"
#undef REDUCE_TAN

#define REDUCE_DOTPROD
__global__ void kDotProduct_kernel()
#include "kNorm.h"
#undef REDUCE_DOTPROD

//-----------------------------------------------------------------------------------------------
extern "C" void kNEBSendRecv(gpuContext gpu, int buff_size)
{
#ifdef MPI
  if (gpu->sim.ShuttleType == 0) {
    int next_node = gpu->sim.beadid;  //replica id, goes from 0 to neb_nbead-1
    if (next_node >= gpu->sim.neb_nbead) {
      next_node = 0;
    }
    int prev_node = gpu->sim.beadid - 2;
    if (prev_node < 0) {
      prev_node = gpu->sim.neb_nbead -1;
    }

    if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
      kRetrieveCoordNL_kernel<<<gpu->blocks, GRID * 3>>>();
    }
    else {
      kRetrieveCoord_kernel<<<gpu->blocks, GRID * 3>>>();
    }
    cudaDeviceSynchronize();
    // Download from the device
    if (!(gpu->bCanMapHostMemory)) {
      gpu->pbDataShuttle->Download();
    }
    int tag1 = 101010;
    int tag2 = 202020;
    MPI_Status status;
    if ( (gpu->sim.beadid & 1) == 0) {
      MPI_Sendrecv(gpu->pbDataShuttle->_pSysData, buff_size, MPI_DOUBLE, next_node, tag1,
                   gpu->pbNextDataShuttle->_pSysData, buff_size, MPI_DOUBLE, next_node, tag1,
                   MPI_COMM_WORLD, &status);
    }
    else {
      MPI_Sendrecv(gpu->pbDataShuttle->_pSysData, buff_size, MPI_DOUBLE, prev_node, tag1,
                   gpu->pbPrevDataShuttle->_pSysData, buff_size, MPI_DOUBLE, prev_node, tag1,
                   MPI_COMM_WORLD, &status);
    }
    if ( (gpu->sim.beadid & 1) == 0) {
      MPI_Sendrecv(gpu->pbDataShuttle->_pSysData, buff_size, MPI_DOUBLE, prev_node, tag2,
                   gpu->pbPrevDataShuttle->_pSysData, buff_size, MPI_DOUBLE, prev_node, tag2,
                   MPI_COMM_WORLD, &status);
    }
    else {
      MPI_Sendrecv(gpu->pbDataShuttle->_pSysData, buff_size, MPI_DOUBLE, next_node, tag2,
                   gpu->pbNextDataShuttle->_pSysData, buff_size, MPI_DOUBLE, next_node, tag2,
                   MPI_COMM_WORLD, &status);
    }

    //Upload to the device
    if (!(gpu->bCanMapHostMemory)) {
      gpu->pbNextDataShuttle->Upload();
      gpu->pbPrevDataShuttle->Upload();
    }
  }
  else if (gpu->sim.ShuttleType == 1) {
    // Print an error for now
    printf("| Error: this functionality (ShuttleType == 1) is not yet available.\n");
  }
#endif
}

extern "C" void NEB_report_energy(gpuContext gpu, int* master_size, double neb_nrg_all[])
{
    //this is needed for energy output
    // at this point the pbNEBEnergyAll is filled with the replica's energies
    for (int i = 0; i < *master_size; i++) {
      neb_nrg_all[i] = gpu->pbNEBEnergyAll->_pSysData[i];
    }
}

// Set the springs:
// at this point the pbNEBEnergyAll is filled with the replica's energies
extern "C" void kNEBspr(gpuContext gpu)
{
    PMEDouble ezero = MAX(gpu->pbNEBEnergyAll->_pSysData[0],
                       gpu->pbNEBEnergyAll->_pSysData[gpu->sim.neb_nbead-1]);
    PMEDouble emax = gpu->pbNEBEnergyAll->_pSysData[0];
    for (int rep = 1; rep < gpu->sim.neb_nbead; rep++) {
      emax = MAX(emax, gpu->pbNEBEnergyAll->_pSysData[rep]);
    }

    for (int i = 0; i < 5; i++) {
      gpu->pbDataSPR->_pSysData[i] = 0.0;
    }

    if (gpu->sim.skmin == gpu->sim.skmax) {
      gpu->pbDataSPR->_pSysData[4] = gpu->sim.skmax;
    }
    else if ((gpu->pbNEBEnergyAll->_pSysData[1] > ezero) && (emax != ezero)) {
      gpu->pbDataSPR->_pSysData[4] = gpu->sim.skmax - gpu->sim.skmin * (emax -
         MAX(gpu->pbNEBEnergyAll->_pSysData[0], gpu->pbNEBEnergyAll->_pSysData[1])/(emax - ezero));
    }
    else {
      gpu->pbDataSPR->_pSysData[4] = gpu->sim.skmax - gpu->sim.skmin;
    }
    int rep = gpu->sim.beadid - 1;  //replica id, goes from 0 to neb_nbead-1
    if ((rep == 0) || (rep == (gpu->sim.neb_nbead - 1))){
      goto exit3;
    }
    gpu->pbDataSPR->_pSysData[3] = gpu->pbDataSPR->_pSysData[4];
    if (gpu->sim.skmin == gpu->sim.skmax) {
      gpu->pbDataSPR->_pSysData[4] = gpu->sim.skmax;
    }
    else if ((gpu->pbNEBEnergyAll->_pSysData[rep + 1] > ezero) && (emax != ezero)) {
      gpu->pbDataSPR->_pSysData[4] = gpu->sim.skmax - gpu->sim.skmin * (emax -
         MAX(gpu->pbNEBEnergyAll->_pSysData[rep - 1], gpu->pbNEBEnergyAll->_pSysData[rep])/(emax - ezero));
    }
    else {
      gpu->pbDataSPR->_pSysData[4] = gpu->sim.skmax - gpu->sim.skmin;
    }

    if (!(gpu->bCanMapHostMemory)) {
      gpu->pbDataSPR->Upload();
    }

    exit3: {}
}

extern "C" void kFitCOM(gpuContext gpu)
{

  double det;
  double small;
  double b[9];
  double norm;
  int ismall;

  //vars needed for lapack diagonalization
  double Kabsch2[9], eigReal[3], eigImag[3];
  double vl[1];
  double vr[9];
  int info;
  int ldvl = 1, ldvr = 3, n = 3;
  const int worksize = 48;
  int lwork = worksize;
  double work[worksize];
  char Nchar='N', Vchar='V';

  if (first) {
    // Calculate total mass of the fit region
    if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
     kTotMassFitNL_kernel<<<gpu->blocks, gpu->threadsPerBlock/4>>>();
      LAUNCHERROR("kFitCOM");
    }
    else {
      kTotMassFit_kernel<<<gpu->blocks, gpu->threadsPerBlock/4>>>();
      LAUNCHERROR("kFitCOM");
    }
    first = false;
  }

  cudaDeviceSynchronize();

    for (int i = 0; i < 27; i++) {
      gpu->pbKabschCOM->_pSysData[i] = 0.0;
    }
    gpu->pbKabschCOM->Upload();

    // Calculate center of mass of the self and neighboring replicas
    if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
      kCALCCOMNL_kernel<<<gpu->blocks, gpu->threadsPerBlock*3/16>>>();
      LAUNCHERROR("kFitCOM");
    }
    else {
      kCALCCOM_kernel<<<gpu->blocks, gpu->threadsPerBlock*3/16>>>();
      LAUNCHERROR("kFitCOM");
    }

    // Center atoms wrt COM of fit regions.
    kCenterData_kernel<<<gpu->blocks, gpu->threadsPerBlock/4>>>();
    LAUNCHERROR("kFitCOM");


    // Calculate Kabsch matrix
    if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
      kKabschPrevNL_kernel<<<gpu->blocks, gpu->threadsPerBlock*3/16>>>();
      LAUNCHERROR("kFitCOM");
      kKabschNextNL_kernel<<<gpu->blocks, gpu->threadsPerBlock*3/16>>>();
      LAUNCHERROR("kFitCOM");
    }
    else {
      kKabschPrev_kernel<<<gpu->blocks, gpu->threadsPerBlock*3/16>>>();
      LAUNCHERROR("kFitCOM");
      kKabschNext_kernel<<<gpu->blocks, gpu->threadsPerBlock*3/16>>>();
      LAUNCHERROR("kFitCOM");
    }
    gpu->pbKabschCOM->Download();

    det = gpu->pbKabschCOM->_pSysData[0] * gpu->pbKabschCOM->_pSysData[4] * gpu->pbKabschCOM->_pSysData[8] -
          gpu->pbKabschCOM->_pSysData[0] * gpu->pbKabschCOM->_pSysData[5] * gpu->pbKabschCOM->_pSysData[7] -
          gpu->pbKabschCOM->_pSysData[1] * gpu->pbKabschCOM->_pSysData[3] * gpu->pbKabschCOM->_pSysData[8] +
          gpu->pbKabschCOM->_pSysData[1] * gpu->pbKabschCOM->_pSysData[5] * gpu->pbKabschCOM->_pSysData[6] +
          gpu->pbKabschCOM->_pSysData[2] * gpu->pbKabschCOM->_pSysData[3] * gpu->pbKabschCOM->_pSysData[7] -
          gpu->pbKabschCOM->_pSysData[2] * gpu->pbKabschCOM->_pSysData[4] * gpu->pbKabschCOM->_pSysData[6];

    if (abs(det) < 1.0e-07) {
      printf("small determinant in rmsfit(): %f\n", abs(det));
      goto exit1;
    }

    for(int i = 0; i < 9; i++) {
      Kabsch2[i] = 0.0;
      b[i] = 0.0;
    }

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          Kabsch2[3*i+j] += gpu->pbKabschCOM->_pSysData[3*k+i]
                          * gpu->pbKabschCOM->_pSysData[3*k+j];
        }
      }
    }

    // Calculate eigenvalues using the DGEEV subroutine
    dgeev_(&Nchar,&Vchar,&n,Kabsch2,&n,eigReal,eigImag,
              vl,&ldvl,vr,&ldvr,work,&lwork,&info);

    // Check for errors
    if (info!=0){
      printf("Error in diagonalization routine dgeev %d\n", info);
      goto exit1;
    }

    // Find the smallest eigenvalue
    small = 1.0e20;
    for (int i = 0; i < 3; i++) {
      if (eigReal[i] < small) {
        ismall = i;
        small = eigReal[i];
      }
    }

    // Calculate the b vectors
    for (int j = 0; j < 3; j++) {
      norm = 1.0/sqrt(abs(eigReal[j]));
      for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
          b[3*i+j] += gpu->pbKabschCOM->_pSysData[3*i+k] * vr[3*j+k] * norm;
        }
      }
    }

    CONTINUE1:

    for(int i = 0; i < 9; i++) {
      gpu->pbRotAtm->_pSysData[i] = 0.0;
    }

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          gpu->pbRotAtm->_pSysData[3*i+j] += b[3*i+k] * vr[3*k+j];
        }
      }
    }

    det = gpu->pbRotAtm->_pSysData[0] * gpu->pbRotAtm->_pSysData[4] * gpu->pbRotAtm->_pSysData[8] -
          gpu->pbRotAtm->_pSysData[0] * gpu->pbRotAtm->_pSysData[5] * gpu->pbRotAtm->_pSysData[7] -
          gpu->pbRotAtm->_pSysData[1] * gpu->pbRotAtm->_pSysData[3] * gpu->pbRotAtm->_pSysData[8] +
          gpu->pbRotAtm->_pSysData[1] * gpu->pbRotAtm->_pSysData[5] * gpu->pbRotAtm->_pSysData[6] +
          gpu->pbRotAtm->_pSysData[2] * gpu->pbRotAtm->_pSysData[3] * gpu->pbRotAtm->_pSysData[7] -
          gpu->pbRotAtm->_pSysData[2] * gpu->pbRotAtm->_pSysData[4] * gpu->pbRotAtm->_pSysData[6];

    if (abs(det) < 1.0e-10) {
      printf("small determinant in rmsfit(): %f\n", abs(det));
      goto exit1;
    }

    if (det < 0) {
      for (int i = 0; i < 3; i++) {
        b[3*i+ismall] = -b[3*i+ismall];
      }
      goto CONTINUE1;
    }

    if (!(gpu->bCanMapHostMemory)) {
      gpu->pbRotAtm->Upload();
    }

    kRotatePrev_kernel<<<gpu->blocks, gpu->threadsPerBlock/4>>>();
    LAUNCHERROR("kFitCOM");

    cudaDeviceSynchronize();

    exit1: {}

    det = gpu->pbKabschCOM->_pSysData[9] * gpu->pbKabschCOM->_pSysData[13] * gpu->pbKabschCOM->_pSysData[17] -
          gpu->pbKabschCOM->_pSysData[9] * gpu->pbKabschCOM->_pSysData[14] * gpu->pbKabschCOM->_pSysData[16] -
          gpu->pbKabschCOM->_pSysData[10] * gpu->pbKabschCOM->_pSysData[12] * gpu->pbKabschCOM->_pSysData[17] +
          gpu->pbKabschCOM->_pSysData[10] * gpu->pbKabschCOM->_pSysData[14] * gpu->pbKabschCOM->_pSysData[15] +
          gpu->pbKabschCOM->_pSysData[11] * gpu->pbKabschCOM->_pSysData[12] * gpu->pbKabschCOM->_pSysData[16] -
          gpu->pbKabschCOM->_pSysData[11] * gpu->pbKabschCOM->_pSysData[13] * gpu->pbKabschCOM->_pSysData[15];

    if (abs(det) < 1.0e-07) {
      printf("small determinant in rmsfit(): %f\n", abs(det));
      goto exit2;
    }

    for(int i = 0; i < 9; i++) {
      Kabsch2[i] = 0.0;
      b[i] = 0.0;
    }

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          Kabsch2[3*i+j] += gpu->pbKabschCOM->_pSysData[3*k+i+9]
                          * gpu->pbKabschCOM->_pSysData[3*k+j+9];
        }
      }
    }

    // Calculate eigenvalues using the DGEEV subroutine
    dgeev_(&Nchar,&Vchar,&n,Kabsch2,&n,eigReal,eigImag,
              vl,&ldvl,vr,&ldvr,work,&lwork,&info);

    // Check for errors
    if (info!=0){
      printf("Error in diagonalization routine dgeev");
      goto exit2;
    }

    // Find the smallest eigenvalue
    small = 1.0e20;
    for (int i = 0; i < 3; i++) {
      if (eigReal[i] < small) {
        ismall = i;
        small = eigReal[i];
      }
    }

    // Calculate the b vectors
    for (int j = 0; j < 3; j++) {
      norm = 1.0/sqrt(abs(eigReal[j]));
      for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
          b[3*i+j] += gpu->pbKabschCOM->_pSysData[3*i+k+9] * vr[3*j+k] * norm;
        }
      }
    }

    CONTINUE2:

    for(int i = 0; i < 9; i++) {
      gpu->pbRotAtm->_pSysData[i] = 0.0;
    }

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          gpu->pbRotAtm->_pSysData[3*i+j] += b[3*i+k] * vr[3*k+j];
        }
      }
    }

    det = gpu->pbRotAtm->_pSysData[0] * gpu->pbRotAtm->_pSysData[4] * gpu->pbRotAtm->_pSysData[8] -
          gpu->pbRotAtm->_pSysData[0] * gpu->pbRotAtm->_pSysData[5] * gpu->pbRotAtm->_pSysData[7] -
          gpu->pbRotAtm->_pSysData[1] * gpu->pbRotAtm->_pSysData[3] * gpu->pbRotAtm->_pSysData[8] +
          gpu->pbRotAtm->_pSysData[1] * gpu->pbRotAtm->_pSysData[5] * gpu->pbRotAtm->_pSysData[6] +
          gpu->pbRotAtm->_pSysData[2] * gpu->pbRotAtm->_pSysData[3] * gpu->pbRotAtm->_pSysData[7] -
          gpu->pbRotAtm->_pSysData[2] * gpu->pbRotAtm->_pSysData[4] * gpu->pbRotAtm->_pSysData[6];

    if (abs(det) < 1.0e-10) {
      printf("small determinant in rmsfit(): %f\n", abs(det));
      goto exit2;
    }

    if (det < 0) {
      for (int i = 0; i < 3; i++) {
        b[3*i+ismall] = -b[3*i+ismall];
      }
      goto CONTINUE2;
    }

    if (!(gpu->bCanMapHostMemory)) {
      gpu->pbRotAtm->Upload();
    }

    kRotateNext_kernel<<<gpu->blocks, gpu->threadsPerBlock/4>>>();
    LAUNCHERROR("kFitCOM");

    kReCenterData_kernel<<<gpu->blocks, gpu->threadsPerBlock/4>>>();
    LAUNCHERROR("kFitCOM");
    cudaDeviceSynchronize();

    exit2: {}

}

extern "C" void kNEBfrc(gpuContext gpu, double neb_force[][3])
{
    if (gpu->sim.tmode == 1) {

      int rep = gpu->sim.beadid - 1;
      if ((gpu->pbNEBEnergyAll->_pSysData[rep + 1] > gpu->pbNEBEnergyAll->_pSysData[rep]) &&
          (gpu->pbNEBEnergyAll->_pSysData[rep] > gpu->pbNEBEnergyAll->_pSysData[rep - 1])) {
        kNEBRevisedTanFIRST_kernel<<<gpu->blocks, GRID * 3>>>();
      }
      else if ((gpu->pbNEBEnergyAll->_pSysData[rep + 1] < gpu->pbNEBEnergyAll->_pSysData[rep]) &&
               (gpu->pbNEBEnergyAll->_pSysData[rep] < gpu->pbNEBEnergyAll->_pSysData[rep - 1])) {
        kNEBRevisedTanSECOND_kernel<<<gpu->blocks, GRID * 3>>>();
      }
      else if (gpu->pbNEBEnergyAll->_pSysData[rep + 1] > gpu->pbNEBEnergyAll->_pSysData[rep - 1]) {
        kNEBRevisedTanTHIRD_kernel<<<gpu->blocks, GRID * 3>>>();
      }
      else {
        kNEBRevisedTanFOURTH_kernel<<<gpu->blocks, GRID * 3>>>();
      }
    }
    else {
      kNEBBasicTan_kernel<<<gpu->blocks, GRID * 3>>>();
    }

    kNormalize_kernel<<<gpu->blocks, GRID * 3>>>();
    cudaDeviceSynchronize();

    kNormalizeTan_kernel<<<gpu->blocks, GRID * 3>>>();
    cudaDeviceSynchronize();

    if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
      kRetrieveForceNL_kernel<<<gpu->blocks, GRID * 3>>>();
    }
    else {
      kRetrieveForce_kernel<<<gpu->blocks, GRID * 3>>>();
    }
    cudaDeviceSynchronize();

      kDotProduct_kernel<<<gpu->blocks, GRID * 3>>>();

    if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
      kNEBFRCNL_kernel<<<gpu->blocks, GRID * 3>>>();
    }
    else {
      kNEBFRC_kernel<<<gpu->blocks, GRID * 3>>>();
    }
}

extern "C" void kNEBfrc_nstep(gpuContext gpu, double neb_force[][3])
{
    if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
      kNEBFRCNL_kernel<<<gpu->blocks, GRID * 3>>>();
    }
    else {
      kNEBFRC_kernel<<<gpu->blocks, GRID * 3>>>();
    }
}
