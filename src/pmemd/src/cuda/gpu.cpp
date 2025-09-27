#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------
#include <stdio.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>
#include <math.h>
#include <string.h>
#include <algorithm>
#include <vector>

#ifdef _WIN32
#include "extra.h"
#endif

#include "gpuContext.h"
#include "gpu.h"
#ifdef GTI
#include "gti_cuda.cuh"
#endif
namespace gppImpl {
  static gpuContext gpu = NULL;
}

using namespace std;
using namespace gppImpl;

#include "matrix.h"
#include "bondRemap.h"

//---------------------------------------------------------------------------------------------
// gpu_startup_: establish the identity of the GPU to use and set up the MPI communicator if
//               running in parallel.
//
// Arguments (MPI only):
//   gpuID:       ID number of the device when multiple GPUs are involved (0 == master)
//   nGpus:       the number of GPUs to use in the calculation
//   comm_number: identifier for the local communicator among GPUs
//---------------------------------------------------------------------------------------------
#ifdef MPI
extern "C" void gpu_startup_(int* gpuID, int* nGpus, int* comm_number)
#else
extern "C" void gpu_startup_(void)
#endif
{
#ifdef GTI
  gpu = theGPUContext::Initialize();
#else
  gpu = new _gpuContext;
#endif
  PRINTMETHOD("gpu_startup");
#ifdef MPI
  gpu->gpuID = *gpuID;
  gpu->nGpus = *nGpus;
  gpu->comm = MPI_COMM_WORLD;
#ifdef GVERBOSE
  printf("GPU %d out of %d processes\n", gpu->gpuID, gpu->nGpus);
#endif

  // Create local communicator from comm number
  MPI_Comm_split(MPI_COMM_WORLD, *comm_number, *gpuID, &gpu->comm);
#endif
}

//---------------------------------------------------------------------------------------------
// gpu_set_device_: set the device to use in the calculation
//
// Arguments:
//   device:   ID number of the GPU (CUDA_VISIBLE_DEVICES environment variable)
//---------------------------------------------------------------------------------------------
extern "C" void gpu_set_device_(int* device)
{
  PRINTMETHOD("gpu_set_device");
  gpu->gpu_device_id = *device;
}

//---------------------------------------------------------------------------------------------
// gpu_init_: master function for initializing the device.  This calls setup for all types of
//            simulations and detects essential information about the GPU, like Streaming
//            Multiprocessor (SM) version.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_init_(void)
{
  PRINTMETHOD("gpu_init");
  int LRFSize = 0;
  int SMCount = 0;
  int SMMajor = 0;
  int SMMinor = 0;

#ifdef MPI
  if (getenv("CUDA_PROFILE") != 0) {
    char profile_log[80];
    if(getenv("CUDA_PROFILE_LOG")) {
      sprintf(profile_log, "../%s%d", getenv("CUDA_PROFILE_LOG"), gpu->gpuID);
    }
    else {
      sprintf(profile_log, "../cu%d.csv", gpu->gpuID);
    }
#ifdef _WIN32
    _putenv_s("CUDA_PROFILE_LOG", profile_log);
    _putenv_s("CUDA_PROFILE_CSV", "1");
#else
    setenv("CUDA_PROFILE_LOG", profile_log, 1);
    setenv("CUDA_PROFILE_CSV", "1", 1);
#endif
  }
#endif

  int device = -1;
  int gpuCount = 0;
  cudaError_t status;
  cudaDeviceProp deviceProp;
  status = cudaGetDeviceCount(&gpuCount);
  RTERROR(status, "cudaGetDeviceCount failed");
  if (gpuCount == 0) {
    printf("No CUDA-capable devices found, exiting.\n");
    cudaDeviceReset();
    exit(-1);
  }

#ifdef MPI
  // Grab node names from all other processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  int length;
  char myName[MPI_MAX_PROCESSOR_NAME + 1];
  char* pName = new char[world_size * (MPI_MAX_PROCESSOR_NAME + 1)];
  int* pNameCount = new int[world_size];
  int* pNameDisp = new int[world_size];
  MPI_Get_processor_name(myName, &length);
  strcpy(&pName[world_rank * (MPI_MAX_PROCESSOR_NAME + 1)], myName);
  for (int i = 0; i < world_size; i++) {
    pNameCount[i] = MPI_MAX_PROCESSOR_NAME + 1;
    pNameDisp[i] = i*(MPI_MAX_PROCESSOR_NAME + 1);
  }
  MPI_Allgatherv(myName, MPI_MAX_PROCESSOR_NAME + 1, MPI_CHAR, pName, pNameCount, pNameDisp,
                 MPI_CHAR, MPI_COMM_WORLD);

  // Test for single node run
  bool bSingleNode = true;
  bool bP2P = false;
  bool bNCCL = false;
  for (int i = 0; i < gpu->nGpus; i++) {
    if (strcmp(&pName[i * (MPI_MAX_PROCESSOR_NAME + 1)], myName)) {
      bSingleNode = false;
    }
  }
#endif

  // Activate zero-copy
  cudaSetDeviceFlags(cudaDeviceMapHost);

  // If the device id is -1 this means it was left at the default so we
  // or the user specified -1 on the command line meaning autoselect GPU.
  // choose CUDA device with the most memory (flesh out later for multi-gpu)
  // otherwise we use the device id specified on the command line.
  if (gpu->gpu_device_id == -1) {
#ifdef MPI
    // Check for duplicate processes on current node
    int localCount = 0;
    int offset = 1;
    for (int i = 0; i < world_size; i++) {
      if (!strcmp(&pName[i * (MPI_MAX_PROCESSOR_NAME + 1)], myName)) {
        localCount++;
        if (i < world_rank) {
          offset++;
        }
      }
    }
    if (localCount > 1) {

      // Choose nth gpu that can run AMBER
      int pos = 0;
      while (offset > 0) {
        cudaGetDeviceProperties(&deviceProp, pos);
        if (deviceProp.canMapHostMemory &&
            ((deviceProp.major >= 3) ||
             ((deviceProp.major == 1) && (deviceProp.minor == 3)))) {
          device = pos;
          offset--;
        }
        pos++;
        if (pos == gpuCount) {
          pos = 0;
        }
      }
#  ifdef GVERBOSE
      char hostname[128];
      gethostname(hostname, 127);
      printf("Node %d running on device %d out of %d GPUs on %s\n", gpu->gpuID, device,
             gpuCount, hostname);
#  endif
    }
    else {
#endif

    // Generate list of compatible GPUs scored by GPU revision first and total memory second
    int* pGPUList = new int[gpuCount];
    unsigned int* pGPUScore = new unsigned int[gpuCount];
    int gpus = 0;
    for (int i = 0; i < gpuCount; i++) {
      cudaGetDeviceProperties(&deviceProp, i);
#ifdef MPI
      if (deviceProp.canMapHostMemory) {
#endif
      if (((deviceProp.major >= 3) ||
           ((deviceProp.major == 1) && (deviceProp.minor == 3)))) {
        pGPUList[gpus] = i;
        pGPUScore[gpus] = (deviceProp.major << 24) + (deviceProp.totalGlobalMem >> 20);
        gpus += (deviceProp.major >= 3);
      }
#ifdef MPI
      }
#endif
    }
    if (gpus == 0) {
      printf("Error searching for compatible GPU");
      exit(-1);
    }

    // Select best GPU according to score
    if (gpus > 0) {

      // Bubble sort (go bubblesort!) device list by score
      bool done = true;
      do {
        done = true;
        for (int i = 0; i < gpus - 1; i++) {
          if (pGPUScore[i] < pGPUScore[i + 1]) {
            done = false;
            int gpu = pGPUList[i];
            unsigned int score = pGPUScore[i];
            pGPUList[i] = pGPUList[i + 1];
            pGPUScore[i] = pGPUScore[i + 1];
            pGPUList[i + 1] = gpu;
            pGPUScore[i + 1] = score;
          }
        }
      } while (!done);
    }

    // Let CUDA select any device from this list
    status = cudaSetValidDevices(pGPUList, gpus);
    delete[] pGPUList;
    delete[] pGPUScore;
    RTERROR(status, "Error searching for compatible GPU");
    // Trick driver into creating a context on an available and valid GPU
    status = cudaFree(0);
    RTERROR(status, "Error selecting compatible GPU");

    // Get device
    status = cudaGetDevice(&device);
    RTERROR(status, "Error fetching current GPU");
#ifdef MPI
  }
#endif
  }
  else {
    cudaGetDeviceProperties(&deviceProp, gpu->gpu_device_id);
#ifdef MPI
    if (deviceProp.canMapHostMemory && (deviceProp.major >= 3)) {
#else
    if (deviceProp.major >= 3) {
#endif
      device = gpu->gpu_device_id;
    }
    else {
#ifdef MPI
      printf("Selected GPU does not support both zero-copy and SM 3.0, exiting.\n");
#else
      printf("Selected GPU lacks SM 3.0 or better support, exiting.\n");
#endif
      cudaDeviceReset();
      exit(-1);
    }
  }

#ifdef MPI
  // Release list of gpu names
  delete[] pName;
  delete[] pNameCount;
  delete[] pNameDisp;
#endif

  if (device == -1) {
#ifdef MPI
    printf("No zero-copy and double-precision capable gpu located, exiting.\n");
#else
    printf("No double-precision capable gpu located, exiting.\n");
#endif
    cudaDeviceReset();
    exit(-1);
  }


  // Finally set CUDA device
  status = cudaSetDevice(device);
  RTERROR(status, "Error setting CUDA device");
  cudaDeviceSynchronize();

  // Test for universal P2P access
#if defined(MPI)
#ifdef NCCL
  if (bSingleNode && gpu->nGpus > 2) {
    bNCCL = true;
  } else
#endif
#ifdef AMBER_PLATFORM_AMD
  if (bSingleNode) {
#else
  if (bSingleNode && ((gpu->nGpus & (gpu->nGpus - 1)) == 0)) {
#endif
    bP2P = true;
    int* pDevice = new int[gpu->nGpus];
    pDevice[gpu->gpuID] = device;
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, pDevice, sizeof(int), MPI_BYTE,
                  gpu->comm);
    int* pUnifiedAddressing = new int[gpu->nGpus];
    cudaGetDeviceProperties(&deviceProp, device);
#ifdef AMBER_PLATFORM_AMD
    pUnifiedAddressing[gpu->gpuID] = true;
#else
    pUnifiedAddressing[gpu->gpuID] = deviceProp.unifiedAddressing;
#endif
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, pUnifiedAddressing, sizeof(int),
                  MPI_BYTE, gpu->comm);
    for (int i = 0; i < gpu->nGpus; i++) {
      if (pDevice[i] != device) {
        int canAccessPeer;
        cudaError_t status = cudaDeviceCanAccessPeer(&canAccessPeer, device, pDevice[i]);
        RTERROR(status, "cudaDeviceCanAccessPeer");
        if (canAccessPeer == 0) {
          bP2P = false;
        }
        else {
          status = cudaDeviceEnablePeerAccess(pDevice[i], 0);

          // Ignore error that really isn't an error, just bad API design
          if (status != cudaErrorPeerAccessAlreadyEnabled) {
            RTERROR(status, "cudaDeviceEnablePeerAccess");
          }
          else {
            cudaGetLastError();
          }
        }
      }
      if (!pUnifiedAddressing[i]) {
        bSingleNode = false;
      }
    }
    delete[] pDevice;
  }
  gpu->bSingleNode = bSingleNode;
  gpu->bP2P = bP2P;
#ifdef AMBER_PLATFORM_AMD
  gpu->bPersistentP2P = gpu->bP2P;
#else
  gpu->bPersistentP2P = false;
#endif
  gpu->bNCCL = bNCCL;

#ifdef NCCL
  // Allocate NCCL communicator
  if (gpu->bNCCL) {
    int nrank, rank;
    MPI_Comm_size(gpu->comm, &nrank);
    MPI_Comm_rank(gpu->comm, &rank);
    ncclUniqueId ncclId;

    if (rank == 0) {
      ncclGetUniqueId(&ncclId);
    }
    MPI_Bcast(&ncclId, NCCL_UNIQUE_ID_BYTES, MPI_CHAR, 0, gpu->comm);
    ncclResult_t ret = ncclCommInitRank(&(gpu->ncclComm), nrank, ncclId, rank);
    if (ret != ncclSuccess) {
      printf("NCCL Init failed : %d\n", ret);
      exit(-1);
    }
  }
#endif
#endif

  // Determine kernel call configuration and grab desired additional GPU properties
  cudaGetDeviceProperties(&deviceProp, device);
  //gpu->bECCSupport = deviceProp.ECCEnabled || deviceProp.tccDriver ||
  //                   (strcasestr(deviceProp.name, "tesla") != NULL);
  // TL: comment out "telsa"--ECC support can be on non-tesla GPUs
  gpu->bECCSupport = deviceProp.ECCEnabled || deviceProp.tccDriver ;
  gpu->bCanMapHostMemory = deviceProp.canMapHostMemory;
  gpu->totalMemory = deviceProp.totalGlobalMem;

#ifdef GVERBOSE
  double memsize = (double)deviceProp.totalGlobalMem / (1024.0 * 1024.0);
  printf("Using GPU %d, %s, SM %d.%d, %.1f MBytes of memory\n", device, deviceProp.name,
         deviceProp.major, deviceProp.minor, memsize);
#endif

  // Store GPU Device ID for later use
  gpu->gpu_device_id = device;
  gpu->blocks = deviceProp.multiProcessorCount;

  // Determine SM version
  unsigned int blocksPerSM;
  gpu->sm_version = SM_3X;
  gpu->threadsPerBlock = THREADS_PER_BLOCK;
  gpu->NLBuildNeighborList32ThreadsPerBlock = NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK;
  gpu->NLBuildNeighborList16ThreadsPerBlock = NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK;
  gpu->NLBuildNeighborList8ThreadsPerBlock = NLBUILD_NEIGHBORLIST8_THREADS_PER_BLOCK;
  gpu->BNLBlocks = NLBUILD_NEIGHBORLIST_BLOCKS_MULTIPLIER * gpu->blocks;
  gpu->localForcesThreadsPerBlock = LOCALFORCES_THREADS_PER_BLOCK;
  gpu->CHARMMForcesThreadsPerBlock = CHARMMFORCES_THREADS_PER_BLOCK;
  gpu->AFEExchangeThreadsPerBlock = AFE_EXCHANGE_THREADS_PER_BLOCK;
  gpu->NMRForcesThreadsPerBlock = NMRFORCES_THREADS_PER_BLOCK;
  gpu->clearForcesThreadsPerBlock = CLEARFORCES_THREADS_PER_BLOCK;
  gpu->NLClearForcesThreadsPerBlock = NLCLEARFORCES_THREADS_PER_BLOCK;
  gpu->reduceForcesThreadsPerBlock = REDUCEFORCES_THREADS_PER_BLOCK;
  gpu->reduceBufferThreadsPerBlock = REDUCEBUFFER_THREADS_PER_BLOCK;
  gpu->GBBornRadiiThreadsPerBlock = GBBORNRADII_THREADS_PER_BLOCK;
  gpu->GBBornRadiiIGB78ThreadsPerBlock = GBBORNRADII_THREADS_PER_BLOCK;
  gpu->GBBornRadiiBlocks = GBBORNRADII_BLOCKS_MULTIPLIER * gpu->blocks;
  gpu->GBNonbondEnergy1ThreadsPerBlock = GBNONBONDENERGY1_THREADS_PER_BLOCK;
  gpu->GBNonbondEnergy1Blocks = GBNONBONDENERGY1_BLOCKS_MULTIPLIER * gpu->blocks;
  gpu->GBNonbondEnergy2ThreadsPerBlock = GBNONBONDENERGY2_THREADS_PER_BLOCK;
  gpu->GBNonbondEnergy2IGB78ThreadsPerBlock = GBNONBONDENERGY2_THREADS_PER_BLOCK;
  gpu->GBNonbondEnergy2Blocks = GBNONBONDENERGY2_BLOCKS_MULTIPLIER * gpu->blocks;
  gpu->PMENonbondEnergyThreadsPerBlock = PMENONBONDENERGY_THREADS_PER_BLOCK;
  gpu->PMENonbondForcesThreadsPerBlock = PMENONBONDFORCES_THREADS_PER_BLOCK;
  gpu->PMENonbondBlocks = PMENONBONDENERGY_BLOCKS_MULTIPLIER * gpu->blocks;
  gpu->IPSNonbondEnergyThreadsPerBlock = IPSNONBONDENERGY_THREADS_PER_BLOCK;
  gpu->IPSNonbondForcesThreadsPerBlock = IPSNONBONDFORCES_THREADS_PER_BLOCK;
  gpu->IPSNonbondBlocks = IPSNONBONDENERGY_BLOCKS_MULTIPLIER * gpu->blocks;
  gpu->updateThreadsPerBlock = UPDATE_THREADS_PER_BLOCK;
  gpu->shakeThreadsPerBlock = SHAKE_THREADS_PER_BLOCK;
  gpu->generalThreadsPerBlock = GENERAL_THREADS_PER_BLOCK;
  gpu->readSize = READ_SIZE;
  gpu->maxSoluteMolecules = MAXMOLECULES;
  gpu->maxPSSoluteMolecules = MAXPSMOLECULES;
  gpu->bNeighborList = false;
#ifdef GVERBOSE
  printf("LRF size is %d\n", deviceProp.regsPerBlock);
#endif

  // Support for architectures < 300 was dropped: all GPUs are assumed to support DP textures.
  gpu->bNoDPTexture       = false;

  // Set up kernel cache and shared memory preferences
  kCalculateGBBornRadiiInitKernels(gpu);
  kCalculateGBNonbondEnergy1InitKernels(gpu);
  kCalculateGBNonbondEnergy2InitKernels(gpu);
  kNeighborListInitKernels(gpu);
  kCalculatePMENonbondEnergyInitKernels(gpu);
  kPMEInterpolationInitKernels(gpu);
  kCalculateLocalForcesInitKernels(gpu);
  kShakeInitKernels(gpu);
  return;
}

//---------------------------------------------------------------------------------------------
// gpu_shutdown_: function for shutting down GPUs and deleting memory allocated on the host
//                associated with them.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_shutdown_(void)
{
  PRINTMETHOD("gpu_shutdown");
#ifdef MPI
  if (gpu->bP2P) {
    for (int i = 0; i < gpu->nGpus; i++) {
      if (i != gpu->gpuID) {
        cudaError_t status;
        status = cudaIpcCloseMemHandle((void*)gpu->pbPeerAccumulatorList->_pSysData[i]);
        RTERROR(status, "cudaIpcCloseMemHandle failed on gpu->pbPeerAccumulator");
      }
    }

    // Make sure all processes reach here before deleting GPU memory
    MPI_Barrier(gpu->comm);
  }
#ifdef NCCL
  if (gpu->bNCCL) {
    ncclCommDestroy(gpu->ncclComm);
  }
#endif
#endif

#if defined(GTI)
  theGPUContext::Shutdown();
#else
  delete gpu;
  cudaDeviceReset();
#endif
  return;
}

//---------------------------------------------------------------------------------------------
// gpu_get_device_info_: obtains information on the GPUs installed and what they are called.
//
// Arguments:
//   gpu_dev_count:  the number of GPU devices
//   gpu_dev_id:     identifications of the GPU device for which to obtain information
//   gpu_dev_mem:    total device memory of this GPU
//   gpu_num_proc:   the number of CUDA cores in this GPU
//   gpu_core_freq:  cycles/second of each core (i.e. 1.4GHz)
//   gpu_name:       the name of this GPU device
//   name_len:       length of the device name
//---------------------------------------------------------------------------------------------
extern "C" void gpu_get_device_info_(int* gpu_dev_count, int* gpu_dev_id, int* gpu_dev_mem,
                                     int* gpu_num_proc, double* gpu_core_freq, char* gpu_name,
                                     int* name_len)
{
  PRINTMETHOD("gpu_get_device_info");
  cudaError_t status;
  cudaDeviceProp deviceProp;
  size_t device_mem;

  status = cudaGetDeviceCount(gpu_dev_count);
  RTERROR(status, "cudaGetDeviceCount failed");

  *gpu_dev_id = gpu->gpu_device_id;
  cudaGetDeviceProperties(&deviceProp, *gpu_dev_id);

  device_mem = deviceProp.totalGlobalMem/(1024*1024);
  *gpu_dev_mem = (int )device_mem;

  *gpu_num_proc = (int )deviceProp.multiProcessorCount;
  *gpu_core_freq = (double )(deviceProp.clockRate * 1e-6f);

  strcpy(gpu_name,deviceProp.name);
  *name_len=strlen(deviceProp.name);
}

//---------------------------------------------------------------------------------------------
// gpu_check_titan_: function to check Amber compatibility with NVIDIA GTX TITAN or 780
//---------------------------------------------------------------------------------------------
extern "C" void gpu_check_titan_()
{
  PRINTMETHOD("gpu_check_titan");

  // Variables to get GPU name
  int gpu_dev_id = 0;
  cudaDeviceProp deviceProp;
  char* gpu_name;
  gpu_name = (char *) malloc(128 * sizeof(char));

  // Variables to get driver version
  float driver_version = 0.0f;
  FILE* nvidia_version_file;
  char* version_line;
  version_line = (char *) malloc(128 * sizeof(char));

  // Get GPU name
  gpu_dev_id = gpu->gpu_device_id;
  cudaGetDeviceProperties(&deviceProp, gpu_dev_id);
  strcpy(gpu_name,deviceProp.name);

  // Check for GeForce GTX TITAN or 780
  if (strcmp(gpu_name, "GeForce GTX TITAN") == 0 || strcmp(gpu_name, "GeForce GTX 780") == 0) {

    // Get driver version from /proc/driver/nvidia/version
    // Note that this may be linux-specific
    nvidia_version_file = fopen("/proc/driver/nvidia/version", "r");
    if (nvidia_version_file != NULL) {
      char* throwaway = fgets(version_line, 128, nvidia_version_file);
      if (version_line != NULL) {
        sscanf(version_line, "%*s %*s %*s %*s %*s %*s %*s %f", &driver_version);
      }
    }
    fclose(nvidia_version_file);

    // Check for NVIDIA driver version 325.15.
    // Note on NVIDIA drivers: driver_branch.build_in_branch (not necessarily ordered in time).
    // This works for now, but because of the way NVIDIA drivers are organized,
    // this may not always be correct.
    if (driver_version < 325.15f) {
      printf("Error: NVIDIA graphics card and NVIDIA driver version "
             "incompatible with Amber\n");
      printf("GeForce GTX TITAN or GeForce GTX 780 need at least driver "
             "version 325.15 for Amber\n");
      exit(-1);
    }
  }
  free(gpu_name);
  free(version_line);
}

#ifdef MPI
//---------------------------------------------------------------------------------------------
// gpu_get_p2p_info_: get peer-to-peer information for a GPU device
//
// Arguments:
//   b_single_node:
//   b_p2p_enabled:
//   b_nccl_enabled:
//---------------------------------------------------------------------------------------------
extern "C" void gpu_get_p2p_info_(bool* b_single_node, bool* b_p2p_enabled, bool* b_nccl_enabled)
{
  PRINTMETHOD("gpu_get_p2p_info");
  *b_single_node = gpu->bSingleNode;
  *b_p2p_enabled = gpu->bP2P;
  *b_nccl_enabled = gpu->bNCCL;
}

//---------------------------------------------------------------------------------------------
// gpu_get_slave_device_info_: get information on a GPU slave device
//
// Arguments:
//   taskID:         the task thread ID (each GPU is a massively parallel architecture, but
//                   from the point of view of the program each GPU device will essentially
//                   get one thread)
//   gpu_dev_count:  the number of GPU devices in use
//   gpu_dev_id:     the ID number of this particular slave device
//   gpu_dev_mem:    the amount of memory in this device
//   gpu_num_proc:   the number of CUDA cores in this slave device
//   gpu_core_freq:  cycles/second of each core (i.e. 1.4GHz)
//   gpu_name:       the name of this GPU device
//   name_len:       length of the device name
//---------------------------------------------------------------------------------------------
extern "C" void gpu_get_slave_device_info_(int* taskID, int* gpu_dev_count, int* gpu_dev_id,
                                           int* gpu_dev_mem, int* gpu_num_proc,
                                           double* gpu_core_freq, char* gpu_name,
                                           int* name_len)
{
  PRINTMETHOD("gpu_get_slave_device_info");
  MPI_Status status;
  MPI_Recv(gpu_dev_count, 1, MPI_INT, *taskID, 0, gpu->comm, &status);
  MPI_Recv(gpu_dev_id, 1, MPI_INT, *taskID, 0, gpu->comm, &status);
  MPI_Recv(gpu_dev_mem, 1, MPI_INT, *taskID, 0, gpu->comm, &status);
  MPI_Recv(gpu_num_proc, 1, MPI_INT, *taskID, 0, gpu->comm, &status);
  MPI_Recv(gpu_core_freq, 1, MPI_DOUBLE, *taskID, 0, gpu->comm, &status);
  MPI_Recv(gpu_name, 80, MPI_CHAR, *taskID, 0, gpu->comm, &status);
  MPI_Recv(name_len, 1, MPI_INT, *taskID, 0, gpu->comm, &status);
}

//---------------------------------------------------------------------------------------------
// gpu_send_slave_device_info_: send information on the slave device to the master process
//---------------------------------------------------------------------------------------------
extern "C" void gpu_send_slave_device_info_()
{
  PRINTMETHOD("gpu_send_slave_device_info");
  int gpu_dev_count;
  int gpu_dev_id;
  int gpu_dev_mem;
  int gpu_num_proc;
  double gpu_core_freq;
  char gpu_name[81];
  int name_len;
  gpu_get_device_info_(&gpu_dev_count, &gpu_dev_id, &gpu_dev_mem, &gpu_num_proc,
                       &gpu_core_freq, gpu_name, &name_len);
  MPI_Send(&gpu_dev_count, 1, MPI_INT, 0, 0, gpu->comm);
  MPI_Send(&gpu_dev_id, 1, MPI_INT, 0, 0, gpu->comm);
  MPI_Send(&gpu_dev_mem, 1, MPI_INT, 0, 0, gpu->comm);
  MPI_Send(&gpu_num_proc, 1, MPI_INT, 0, 0, gpu->comm);
  MPI_Send(&gpu_core_freq, 1, MPI_DOUBLE, 0, 0, gpu->comm);
  MPI_Send(gpu_name, 80, MPI_CHAR, 0, 0, gpu->comm);
  MPI_Send(&name_len, 1, MPI_INT, 0, 0, gpu->comm);
}
#endif

//---------------------------------------------------------------------------------------------
// gpu_get_memory_info_: returns KB of memory in use on CPU and GPU
//
// Arguments:
//   gpumemory: KB of memory on the GPU device
//   cpumemory: KB of memory in the host RAM
//---------------------------------------------------------------------------------------------
extern "C" void gpu_get_memory_info_(int* gpumemory, int* cpumemory)
{
  // The "ll" suffix denotes that the integer be formulated as a long long int
  //*gpumemory  = (int)(gpu->totalGPUMemory / 1024ll);
  //*cpumemory  = (int)(gpu->totalCPUMemory / 1024ll);

  *gpumemory = gpuMemoryInfo::Instance().totalGPUMemory / 1024ll;
  *cpumemory = gpuMemoryInfo::Instance().totalCPUMemory / 1024ll;
  return;
}

//---------------------------------------------------------------------------------------------
// get_shuttle_size_: function for obtaining the size of the information shuttle based on the
//                   type of shuttling and the number of atoms involved.
//
// Arguments:
//   mult:        multiplier on the number of array elements per atom
//   elemsize:    size of each data element
//   natom:       the number of atoms involved
//---------------------------------------------------------------------------------------------
static int get_shuttle_size(int mult, size_t elemsize, int natom)
{
  int dsize;

  // The array must be sized as a multiple of GRID and large enough to hold the required
  // number of atoms.  That means rounding up to the nearest multiple of GRID. 512 bytes
  // will be downloaded at a time, so the array should also be rounded up to the nearest
  // multiple of 512 / (size of data).
  dsize = (mult*natom / GRID) + ((mult*natom % GRID) > 0 ? 1 : 0);
  dsize *= GRID;
  dsize *= elemsize;
  dsize = (dsize / 512) + ((dsize % 512) > 0 ? 1 : 0);
  dsize *= 512;

  // Implicit assumption that elemsize is a factor of 512
  dsize /= elemsize;

  return dsize;
}

//---------------------------------------------------------------------------------------------
// gpu_setup_system_: this is the gateway to the GPU, called from the main program in
//                    pmemd.F90.
//
// Arguments:
//   atoms:     the number of atoms / particles in the system
//   imin:      from &cntrl namelist
//   tol:       goal for SHAKE convergence
//   ntf:       force calculation setting (from &cntrl namelist)
//   ntb:       box type (0 = isolated system, 1 = periodic, 2 = periodic with variable size)
//   ntc:       SHAKE setting (1 = no constraint, 2 = bonds involving hydrogen constraint)
//   ips:       flag to activate Isotropic Periodic Sums
//   ntp:       constant pressure setting (0 = constant volume, 1 = isotropic rescaling)
//   barostat:  barostat setting (1 = Berendsen, 2 = Monte Carlo)
//   ntt:       thermostat setting (0 = no thermostat, implies NVE, 1 = Berendsen,
//              2 = Andersen, 3 = Langevin)
//   gbsa:      Generalized Born Surface Area active flag
//   vrand:     velocity randomization frequency for Andersen thermostating (number of steps)
//   icnstph:   flag to activate constant pH dynamics
//   icnste:    flag to activate constant Redox Potential dynamics
//   hybridgb:  flag to activate Hybrid Solvent REMD
//   ti_mode:   Thermodynamic Integration mode
//   surften:   Surface tension for GBSA
//   ischeme:   Middle-scheme parameter
//   ithermostat: Middle-scheme parameter
//---------------------------------------------------------------------------------------------
extern "C" void gpu_setup_system_(int *atoms, int *imin, double* tol, int* ntf, int* ntb,
                                  int *ips, int *ntp, int *barostat, int *ntt, int *gbsa,
                                  int *vrand, int *icnstph, int* icnste, int *hybridgb,
				                  int *ti_mode, double* surften,
                                  int *ischeme, int *ithermostat, int *ntc)
{
  PRINTMETHOD("gpu_setup_system");

  // Grab simulation parameters;
  gpu->imin   = *imin;
  gpu->ntf    = *ntf;
  gpu->ntb    = *ntb;
  gpu->ips    = *ips;
  gpu->ntt    = *ntt;
  gpu->gbsa   = *gbsa;
  gpu->sim.surften = *surften; //pwsasa need surften
  gpu->sim.ntp = *ntp;
  gpu->sim.barostat = *barostat;
  gpu->sim.icnstph = *icnstph;
  gpu->sim.icnste  = *icnste;
  gpu->sim.hybridgb = *hybridgb;
  gpu->sim.ti_mode = *ti_mode;

  // middle-schem parameters
  gpu->ischeme = *ischeme;
  gpu->ithermostat = *ithermostat;
  gpu->ntc = *ntc;

  // Allocate system based on atom count
  gpu->sim.atoms = *atoms;
  gpu->sim.paddedNumberOfAtoms = ((*atoms + gpu->sim.grid - 1) >> gpu->sim.gridBits) <<
                                 gpu->sim.gridBits;
  gpu->updateBlocks = (*atoms + gpu->updateThreadsPerBlock - 1) / gpu->updateThreadsPerBlock;
  gpu->sim.SMPCount = gpu->blocks;

  // Calculate stride to insure buffers begin on safe texture boundaries
  gpu->sim.stride = ((*atoms + 63) >> 6) << 6;
  gpu->sim.stride2 = 2 * gpu->sim.stride;
  gpu->sim.stride3 = 3 * gpu->sim.stride;
  gpu->sim.stride4 = 4 * gpu->sim.stride;
  gpu->sim.stride5 = 5 * gpu->sim.stride;
  gpu->sim.tol = *tol;

  // Determine number of randoms to generate
  if ((gpu->ntt == 3) || (gpu->ntt == 2) || (gpu->iphmd > 0) || (gpu->ischeme == 1)) {
    gpu->sim.randomSteps = MAX_RANDOM_STEPS;

    // - 1 here to avoid overflowing a signed int. = 2GB.
    if (gpu->totalMemory < aligned_uli(2u * 1024u * 1024u * 1024u )) {
      gpu->sim.randomSteps /= 2;
    }
    if (*atoms >= 131072) {
      gpu->sim.randomSteps /= 2;
    }
    if (*atoms >= 262144) {
      gpu->sim.randomSteps /= 2;
    }
    if (*atoms >= 524288) {
      gpu->sim.randomSteps /= 2;
    }
    gpu->sim.randomNumbers = gpu->sim.paddedNumberOfAtoms * 3 * gpu->sim.randomSteps;
  }
  else {
    gpu->sim.randomNumbers = 1;
  }

  // Clear any previous stuff
  delete gpu->pbAtom;
  delete gpu->pbAtomXYSP;
  delete gpu->pbAtomZSP;
  delete gpu->pbAtomSigEps;
  delete gpu->pbAtomLJID;
  delete gpu->pbAtomRBorn;
  delete gpu->pbAtomS;
  delete gpu->pbAtomCharge;
  delete gpu->pbAtomChargeSP;
  delete gpu->pbAtomChargeSPLJID;
  delete gpu->pbAtomMass;
  delete gpu->pbAtomIonMask;
  delete gpu->pbReff;
  delete gpu->pbFrcBlkCounters;
  gpu->pbReff = NULL;
  delete gpu->pbReffSP;
  delete gpu->pbPsi;
  delete gpu->pbTemp7;
  delete gpu->pbVel;
  delete gpu->pbLVel;
  delete gpu->pbCenter;
  delete gpu->pbRandom;
  delete gpu->pbTIRegion;
  delete gpu->pbTILinearAtmID;
  delete gpu->pbBarLambda;
  delete gpu->pbBarTot;
  gpu->pbAtom               = new GpuBuffer<double>(gpu->sim.stride3 * 3); //fixed to 3 to store SHAKEOLDATOM for middle-scheme
  gpu->pbAtomXYSP           = new GpuBuffer<PMEFloat2>(gpu->sim.stride);
  gpu->pbAtomZSP            = new GpuBuffer<PMEFloat>(gpu->sim.stride);
  gpu->pbAtomCharge         = new GpuBuffer<double>(gpu->sim.stride);
  gpu->pbAtomChargeSP       = new GpuBuffer<PMEFloat>(gpu->sim.stride);
  gpu->pbAtomSigEps         = new GpuBuffer<PMEFloat2>(gpu->sim.stride);
  gpu->pbAtomChargeSPLJID   = new GpuBuffer<PMEFloat2>(gpu->sim.stride);
  gpu->pbAtomLJID           = new GpuBuffer<unsigned int>(gpu->sim.stride);
  gpu->pbAtomRBorn          = new GpuBuffer<PMEFloat>(gpu->sim.stride);
  gpu->pbAtomS              = new GpuBuffer<PMEFloat>(gpu->sim.stride);
  gpu->pbAtomMass           = new GpuBuffer<double>(gpu->sim.stride2);
  gpu->pbAtomIonMask        = new GpuBuffer<unsigned int>(gpu->sim.stride);
  gpu->pbReff               = new GpuBuffer<PMEDouble>(gpu->sim.stride);
  gpu->pbFrcBlkCounters     = new GpuBuffer<unsigned int>(8);
  gpu->sim.pFrcBlkCounters  = gpu->pbFrcBlkCounters->_pDevData;
  if (gpu->ntb == 0 || gpu->sim.icnstph == 2 || gpu->sim.icnste == 2 || gpu->sim.hybridgb > 0) {
    gpu->pbReffSP = new GpuBuffer<PMEFloat>(gpu->sim.stride);
    gpu->pbPsi = new GpuBuffer<PMEFloat>(gpu->sim.stride);
    gpu->pbTemp7 = new GpuBuffer<PMEFloat>(gpu->sim.stride);
  }
  gpu->pbVel = new GpuBuffer<double>(gpu->sim.stride3);
  gpu->pbLVel = new GpuBuffer<double>(gpu->sim.stride3);
  gpu->pbCenter = new GpuBuffer<PMEFloat>(gpu->blocks * 6);
  gpu->pbRandom = new GpuBuffer<double>(gpu->sim.randomNumbers, bShadowedOutputBuffers);

  // Allocate force accumulators
  unsigned int nonbondForceBuffers;
  unsigned int maxForceBuffers;

  // The thermodynamic integration mode no longer requires separate force buffers
  if (gpu->bNeighborList) {
    if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
      nonbondForceBuffers = 2;
    }
    else {
      nonbondForceBuffers = 1;
    }
    maxForceBuffers = nonbondForceBuffers;
  }
  else {
    nonbondForceBuffers = 1;
    maxForceBuffers = 1;
  }

#ifdef use_DPFP
  // Extra buffers for minimization?
  unsigned int intForceBuffers = (gpu->imin != 0) ? 1 : 0;
  gpu->sim.maxForceBuffers = maxForceBuffers + intForceBuffers;
#else
  gpu->sim.maxForceBuffers = maxForceBuffers;
#endif
  gpu->sim.nonbondForceBuffers = nonbondForceBuffers;
  gpu->sim.forceBufferStride = gpu->sim.maxForceBuffers * gpu->sim.stride3;

  // Delete any existing accumulators
  delete gpu->pbForceAccumulator;
  delete gpu->pbReffAccumulator;
#ifdef MPI
  delete gpu->pbPeerAccumulator;
  delete gpu->pbPeerAccumulatorList;
  delete gpu->pPeerAccumulatorMemHandle;
  delete gpu->pbPersistentStages;
#endif
  delete gpu->pbSumdeijdaAccumulator;
  delete gpu->pblliXYZ_q;
  delete gpu->pbXYZ_q;
  delete gpu->pbXYZ_qt;
  gpu->pbXYZ_q = NULL;
  gpu->pbXYZ_qt = NULL;
  delete gpu->pbEnergyBuffer;
  delete gpu->pbAFEBuffer;
  delete gpu->pbKineticEnergyBuffer;
  delete gpu->pbAFEKineticEnergyBuffer;
  
  // Only allocate GB accumulators if they're needed
  if (gpu->ntb == 0 || gpu->sim.icnstph == 2 || gpu->sim.icnste == 2 || gpu->sim.hybridgb > 0) {

    // Allocate Born Radius accumulator
    gpu->pbReffAccumulator = new GpuBuffer<PMEAccumulator>(gpu->sim.stride);
    gpu->sim.pReffAccumulator = gpu->pbReffAccumulator->_pDevData;

    // Allocate Born Force accumulator
    gpu->pbSumdeijdaAccumulator = new GpuBuffer<PMEAccumulator>(gpu->sim.stride);
    gpu->sim.pSumdeijdaAccumulator = gpu->pbSumdeijdaAccumulator->_pDevData;
  }

  // Allocate new force accumulators
  gpu->pbForceAccumulator = new GpuBuffer<PMEAccumulator>(gpu->sim.stride3 * gpu->sim.maxForceBuffers);
  gpu->sim.pForceAccumulator = gpu->pbForceAccumulator->_pDevData;
  gpu->sim.pForceXAccumulator = gpu->pbForceAccumulator->_pDevData;
  gpu->sim.pForceYAccumulator = gpu->sim.pForceAccumulator + gpu->sim.stride;
  gpu->sim.pForceZAccumulator = gpu->sim.pForceAccumulator + gpu->sim.stride2;
  gpu->sim.pBondedForceAccumulator = gpu->sim.pForceAccumulator;
  gpu->sim.pBondedForceXAccumulator = gpu->sim.pForceXAccumulator;
  gpu->sim.pBondedForceYAccumulator = gpu->sim.pForceYAccumulator;
  gpu->sim.pBondedForceZAccumulator = gpu->sim.pForceZAccumulator;
  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
    gpu->sim.pNBForceAccumulator = gpu->pbForceAccumulator->_pDevData + gpu->sim.stride3;
    gpu->sim.pNBForceXAccumulator = gpu->sim.pNBForceAccumulator;
    gpu->sim.pNBForceYAccumulator = gpu->sim.pNBForceAccumulator + gpu->sim.stride;
    gpu->sim.pNBForceZAccumulator = gpu->sim.pNBForceAccumulator + gpu->sim.stride2;
  }
  else {
    gpu->sim.pNBForceAccumulator = gpu->pbForceAccumulator->_pDevData;
    gpu->sim.pNBForceXAccumulator = gpu->sim.pNBForceAccumulator;
    gpu->sim.pNBForceYAccumulator = gpu->sim.pNBForceAccumulator + gpu->sim.stride;
    gpu->sim.pNBForceZAccumulator = gpu->sim.pNBForceAccumulator + gpu->sim.stride2;
  }
#ifdef use_DPFP
  if (gpu->imin != 0)
  {
    gpu->sim.pIntForceAccumulator = gpu->sim.pNBForceAccumulator + gpu->sim.stride3;
    gpu->sim.pIntForceXAccumulator = gpu->sim.pIntForceAccumulator;
    gpu->sim.pIntForceYAccumulator = gpu->sim.pIntForceAccumulator + gpu->sim.stride;
    gpu->sim.pIntForceZAccumulator = gpu->sim.pIntForceAccumulator + gpu->sim.stride2;
  }
  else
#endif
  {
    gpu->sim.pIntForceAccumulator = NULL;
    gpu->sim.pIntForceXAccumulator = NULL;
    gpu->sim.pIntForceYAccumulator = NULL;
    gpu->sim.pIntForceZAccumulator = NULL;
  }

#ifdef MPI
  // Allocate peer accumulator buffer if P2P is active
  if (gpu->bP2P) {
    int levels = 1 << (int)log2(gpu->nGpus * 2 - 1);
    size_t length = gpu->sim.stride3 * levels;
#  ifdef AMBER_PLATFORM_AMD
    // HIP-TODO: We assume that buffers for IPC must be aligned by 2 MiB.
    // Remove when https://github.com/ROCm-Developer-Tools/HIP/issues/2257 is fixed.
    size_t alignment = 2 * 1024 * 1024 / sizeof(PMEAccumulator);
    length = (length + alignment - 1) / alignment * alignment;
#  endif
    gpu->pbPeerAccumulator = new GpuBuffer<PMEAccumulator>(length);
    cudaError_t status;
#  ifdef AMBER_PLATFORM_AMD
    if (gpu->bPersistentP2P) {
      // Recreate with different flags
      status = cudaFree(gpu->pbPeerAccumulator->_pDevData);
      RTERROR(status, "cudaFree failed on gpu->pbPeerAccumulator->_pDevData");
      status = hipExtMallocWithFlags((void**)&gpu->pbPeerAccumulator->_pDevData,
                                     length * sizeof(PMEAccumulator),
                                     hipDeviceMallocFinegrained);
      RTERROR(status, "hipExtMallocWithFlags failed on gpu->pbPeerAccumulator->_pDevData");
      gpu->pbPersistentStages = new GpuBuffer<unsigned int>(1, false, true);
    }
#  endif
    gpu->pPeerAccumulatorMemHandle = new cudaIpcMemHandle_t[gpu->nGpus];
    status = cudaIpcGetMemHandle(&(gpu->pPeerAccumulatorMemHandle[gpu->gpuID]),
                                 gpu->pbPeerAccumulator->_pDevData);
    RTERROR(status, "cudaIpcGetMemHandle failed on gpu->pbPeerAccumulator");
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, gpu->pPeerAccumulatorMemHandle,
                  sizeof(cudaIpcMemHandle_t), MPI_BYTE, gpu->comm);
    gpu->pbPeerAccumulatorList = new GpuBuffer<PMEAccumulator*>(gpu->nGpus);
    for (int i = 0; i < gpu->nGpus; i++) {
      if (i != gpu->gpuID) {
        status = cudaIpcOpenMemHandle((void**)&(gpu->pbPeerAccumulatorList->_pSysData[i]),
                                      gpu->pPeerAccumulatorMemHandle[i],
                                      cudaIpcMemLazyEnablePeerAccess);
        RTERROR(status, "cudaIpcOpenMemHandle failed on gpu->pbPeerAccumulator Handle");
      } else {
        gpu->pbPeerAccumulatorList->_pSysData[i] = NULL;
      }
    }
    gpu->pbPeerAccumulatorList->Upload();
  }
#endif

  // Link up the remaining buffer arrays in the GPU context
  gpu->sim.pAtomX = gpu->pbAtom->_pDevData;
  {
    cudaError_t status;
    if (gpu->sim.texAtomX != 0) {
      status = cudaDestroyTextureObject(gpu->sim.texAtomX);
      RTERROR(status, "cudaDestroyTextureObject gpu->sim.texAtomX failed");
    }
    cudaResourceDesc resDesc;
    cudaTextureDesc texDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceTypeLinear;
    resDesc.res.linear.devPtr = (int2*)gpu->sim.pAtomX;
    resDesc.res.linear.sizeInBytes = gpu->sim.stride3 * sizeof(int2);
    resDesc.res.linear.desc.f = cudaChannelFormatKindSigned;
    resDesc.res.linear.desc.x = 32;
    resDesc.res.linear.desc.y = 32;
    resDesc.res.linear.desc.z = 0;
    resDesc.res.linear.desc.w = 0;
    memset(&texDesc, 0, sizeof(texDesc));
    texDesc.normalizedCoords = 0;
    texDesc.filterMode = cudaFilterModePoint;
    texDesc.addressMode[0] = cudaAddressModeClamp;
    texDesc.readMode = cudaReadModeElementType;
    status = cudaCreateTextureObject(&(gpu->sim.texAtomX), &resDesc, &texDesc, NULL);
    RTERROR(status, "cudaCreateTextureObject gpu->sim.texAtomX failed");
  }
  gpu->sim.pAtomY = gpu->pbAtom->_pDevData + gpu->sim.stride;
  gpu->sim.pAtomXYSP = gpu->pbAtomXYSP->_pDevData;
#if !defined(use_DPFP)
  {
    cudaError_t status;
    if (gpu->sim.texAtomXYSP != 0) {
      status = cudaDestroyTextureObject(gpu->sim.texAtomXYSP);
      RTERROR(status, "cudaDestroyTextureObject gpu->sim.texAtomXYSP failed");
    }
    cudaResourceDesc resDesc;
    cudaTextureDesc texDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceTypeLinear;
    resDesc.res.linear.devPtr = gpu->sim.pAtomXYSP;
    resDesc.res.linear.sizeInBytes = gpu->sim.stride * sizeof(float2);
    resDesc.res.linear.desc.f = cudaChannelFormatKindFloat;
    resDesc.res.linear.desc.x = 32;
    resDesc.res.linear.desc.y = 32;
    resDesc.res.linear.desc.z = 0;
    resDesc.res.linear.desc.w = 0;
    memset(&texDesc, 0, sizeof(texDesc));
    texDesc.normalizedCoords = 0;
    texDesc.filterMode = cudaFilterModePoint;
    texDesc.addressMode[0] = cudaAddressModeClamp;
    texDesc.readMode = cudaReadModeElementType;
    status = cudaCreateTextureObject(&(gpu->sim.texAtomXYSP), &resDesc, &texDesc, NULL);
    RTERROR(status, "cudaCreateTextureObject gpu->sim.texAtomXYSP failed");
  }
#endif
  gpu->sim.pAtomZ = gpu->pbAtom->_pDevData + gpu->sim.stride2;
  gpu->sim.pAtomZSP = gpu->pbAtomZSP->_pDevData;
#if !defined(use_DPFP)
  {
    cudaError_t status;
    if (gpu->sim.texAtomZSP != 0) {
      status = cudaDestroyTextureObject(gpu->sim.texAtomZSP);
      RTERROR(status, "cudaDestroyTextureObject gpu->sim.texAtomZSP failed");
    }
    cudaResourceDesc resDesc;
    cudaTextureDesc texDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceTypeLinear;
    resDesc.res.linear.devPtr = gpu->sim.pAtomZSP;
    resDesc.res.linear.sizeInBytes = gpu->sim.stride * sizeof(float);
    resDesc.res.linear.desc.f = cudaChannelFormatKindFloat;
    resDesc.res.linear.desc.x = 32;
    resDesc.res.linear.desc.y = 0;
    resDesc.res.linear.desc.z = 0;
    resDesc.res.linear.desc.w = 0;
    memset(&texDesc, 0, sizeof(texDesc));
    texDesc.normalizedCoords = 0;
    texDesc.filterMode = cudaFilterModePoint;
    texDesc.addressMode[0] = cudaAddressModeClamp;
    texDesc.readMode = cudaReadModeElementType;
    status = cudaCreateTextureObject(&(gpu->sim.texAtomZSP), &resDesc, &texDesc, NULL);
    RTERROR(status, "cudaCreateTextureObject gpu->sim.texAtomZSP failed");
  }
#endif
  gpu->sim.pOldAtomX = gpu->pbAtom->_pDevData + gpu->sim.stride3;
  {
    cudaError_t status;
    if (gpu->sim.texOldAtomX != 0) {
      status = cudaDestroyTextureObject(gpu->sim.texOldAtomX);
      RTERROR(status, "cudaDestroyTextureObject gpu->sim.texOldAtomX failed");
    }
    cudaResourceDesc resDesc;
    cudaTextureDesc texDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceTypeLinear;
    resDesc.res.linear.devPtr = (int2*)gpu->sim.pOldAtomX;
    resDesc.res.linear.sizeInBytes = gpu->sim.stride3 * sizeof(int2);
    resDesc.res.linear.desc.f = cudaChannelFormatKindSigned;
    resDesc.res.linear.desc.x = 32;
    resDesc.res.linear.desc.y = 32;
    resDesc.res.linear.desc.z = 0;
    resDesc.res.linear.desc.w = 0;
    memset(&texDesc, 0, sizeof(texDesc));
    texDesc.normalizedCoords = 0;
    texDesc.filterMode = cudaFilterModePoint;
    texDesc.addressMode[0] = cudaAddressModeClamp;
    texDesc.readMode = cudaReadModeElementType;
    status = cudaCreateTextureObject(&(gpu->sim.texOldAtomX), &resDesc, &texDesc, NULL);
    RTERROR(status, "cudaCreateTextureObject gpu->sim.texOldAtomX failed");
  }
  gpu->sim.pOldAtomY = gpu->pbAtom->_pDevData + gpu->sim.stride3 + gpu->sim.stride;
  gpu->sim.pOldAtomZ = gpu->pbAtom->_pDevData + gpu->sim.stride3 + gpu->sim.stride2;
  gpu->sim.pShakeOldAtomX = gpu->pbAtom->_pDevData + gpu->sim.stride3 * 2; // added by zhf
  gpu->sim.pShakeOldAtomY = gpu->pbAtom->_pDevData + gpu->sim.stride3 * 2 + gpu->sim.stride;
  gpu->sim.pShakeOldAtomZ = gpu->pbAtom->_pDevData + gpu->sim.stride3 * 2 + gpu->sim.stride2;
  gpu->sim.pAtomSigEps = gpu->pbAtomSigEps->_pDevData;
  gpu->sim.pAtomLJID = gpu->pbAtomLJID->_pDevData;
  gpu->sim.pAtomRBorn = gpu->pbAtomRBorn->_pDevData;
  gpu->sim.pAtomS = gpu->pbAtomS->_pDevData;
  gpu->sim.pAtomCharge = gpu->pbAtomCharge->_pDevData;
  gpu->sim.pAtomChargeSP = gpu->pbAtomChargeSP->_pDevData;
  gpu->sim.pAtomChargeSPLJID = gpu->pbAtomChargeSPLJID->_pDevData;
  gpu->sim.pAtomIonMask = gpu->pbAtomIonMask->_pDevData;
#if !defined(use_DPFP)
  {
    cudaError_t status;
    if (gpu->sim.texAtomChargeSPLJID != 0) {
      status = cudaDestroyTextureObject(gpu->sim.texAtomChargeSPLJID);
      RTERROR(status, "cudaDestroyTextureObject gpu->sim.texAtomChargeSPLJID failed");
    }
    cudaResourceDesc resDesc;
    cudaTextureDesc texDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceTypeLinear;
    resDesc.res.linear.devPtr = gpu->sim.pAtomChargeSPLJID;
    resDesc.res.linear.sizeInBytes = gpu->sim.stride * sizeof(float2);
    resDesc.res.linear.desc.f = cudaChannelFormatKindFloat;
    resDesc.res.linear.desc.x = 32;
    resDesc.res.linear.desc.y = 32;
    resDesc.res.linear.desc.z = 0;
    resDesc.res.linear.desc.w = 0;
    memset(&texDesc, 0, sizeof(texDesc));
    texDesc.normalizedCoords = 0;
    texDesc.filterMode = cudaFilterModePoint;
    texDesc.addressMode[0] = cudaAddressModeClamp;
    texDesc.readMode = cudaReadModeElementType;
    status = cudaCreateTextureObject(&(gpu->sim.texAtomChargeSPLJID), &resDesc, &texDesc, NULL);
    RTERROR(status, "cudaCreateTextureObject gpu->sim.texAtomChargeSPLJID failed");
  }
#endif
  gpu->sim.pAtomMass = gpu->pbAtomMass->_pDevData;
  gpu->sim.pAtomInvMass = gpu->pbAtomMass->_pDevData + gpu->sim.stride;
  if (gpu->ntb == 0 || gpu->sim.icnstph == 2 || gpu->sim.icnste == 2 || gpu->sim.hybridgb > 0) {
    gpu->sim.pReff = gpu->pbReff->_pDevData;
    gpu->sim.pReffSP = gpu->pbReffSP->_pDevData;
    gpu->sim.pPsi = gpu->pbPsi->_pDevData;
    gpu->sim.pTemp7 = gpu->pbTemp7->_pDevData;
  }
  gpu->sim.pVelX = gpu->pbVel->_pDevData;
  gpu->sim.pVelY = gpu->pbVel->_pDevData + gpu->sim.stride;
  gpu->sim.pVelZ = gpu->pbVel->_pDevData + gpu->sim.stride2;
  gpu->sim.pLVelX = gpu->pbLVel->_pDevData;
  gpu->sim.pLVelY = gpu->pbLVel->_pDevData + gpu->sim.stride;
  gpu->sim.pLVelZ = gpu->pbLVel->_pDevData + gpu->sim.stride2;
  gpu->sim.pXMin = gpu->pbCenter->_pDevData + 0 * gpu->blocks;
  gpu->sim.pYMin = gpu->pbCenter->_pDevData + 1 * gpu->blocks;
  gpu->sim.pZMin = gpu->pbCenter->_pDevData + 2 * gpu->blocks;
  gpu->sim.pXMax = gpu->pbCenter->_pDevData + 3 * gpu->blocks;
  gpu->sim.pYMax = gpu->pbCenter->_pDevData + 4 * gpu->blocks;
  gpu->sim.pZMax = gpu->pbCenter->_pDevData + 5 * gpu->blocks;
  gpu->sim.pRandom = gpu->pbRandom->_pDevData;
  gpu->sim.pRandomX = gpu->sim.pRandom;
  gpu->sim.pRandomY = gpu->sim.pRandomX + gpu->sim.randomSteps * gpu->sim.paddedNumberOfAtoms;
  gpu->sim.pRandomZ = gpu->sim.pRandomY + gpu->sim.randomSteps * gpu->sim.paddedNumberOfAtoms;
  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gpu_setup_softcore_ti_: setup for softcore Thermodynamic integration routines.
//
//  Arguments:
//    ti_lst_repacked:      1-D array containing corresponding to ti_lst in ti.F90
//    ti_sc_lst:            see ti.F90
//    ti_latm_lst_repacked: 1-D array corresponding to ti_latm_lst in ti.F90
//    ti_latm_cnt:          number of linear scaling atoms in each region
//---------------------------------------------------------------------------------------------
extern "C" void gpu_setup_softcore_ti_(int ti_lst_repacked[], int ti_sc_lst[],
                                       int ti_latm_lst_repacked[], int *ti_latm_cnt)
{
  PRINTMETHOD("gpu_setup_softcore_ti");

  // Delete any existing data
  delete gpu->pbTIRegion;
  delete gpu->pbTILinearAtmID;

  // Allocate TI Region data
  gpu->pbTIRegion = new GpuBuffer<int>(gpu->sim.stride);
  for (int i = 0; i < gpu->sim.atoms; i++) {

    // Needs to be the first two values in ti_lst_repacked, the third is common
    gpu->pbTIRegion->_pSysData[i] = (ti_sc_lst[i] != 0) | (ti_lst_repacked[i] << 1) |
                                    (ti_lst_repacked[i + gpu->sim.atoms] << 2);
  }
  for (int i = gpu->sim.atoms; i < gpu->sim.stride; i++) {

    // This will work up to a 10 million atom system. then the ints will blow up
    gpu->pbTIRegion->_pSysData[i]   = 99999999 + i * 200;
  }
  gpu->pbTIRegion->Upload();

  int padded_latm_cnt;
  padded_latm_cnt = ((*ti_latm_cnt + 31) >> 5) << 5;
  gpu->pbTILinearAtmID = new GpuBuffer<int>(2 * padded_latm_cnt);
  for (int i = 0; i < *ti_latm_cnt; i++) {
    for (int j = 0; j < 2; j++) {

      // Use the number of ti_latms as the offset allowing us to have one array
      // with both sets of atom IDs
      gpu->pbTILinearAtmID->_pSysData[i + j * padded_latm_cnt] =
        (ti_latm_lst_repacked[i + j*(*ti_latm_cnt)] - 1);
    }
  }
  for (int i = *ti_latm_cnt; i < padded_latm_cnt; i++) {
    for (int j = 0; j < 2; j++) {

      // Used 999999 because it is unlikely anyone will ever
      // want to run something with >999999 linear scaling atoms
      gpu->pbTILinearAtmID->_pSysData[i + j * padded_latm_cnt] = (999999 * i + 100 + j);
    }
  }
  gpu->pbTILinearAtmID->Upload();

  // We don't need to know the number of threads on the gpu, just on the cpu,
  // so don't bother setting gpu->sim
  gpu->AFEExchangeBlocks = (*ti_latm_cnt + gpu->AFEExchangeThreadsPerBlock - 1) /
                           gpu->AFEExchangeThreadsPerBlock;

  // This will work even if *ti_latm_cnt = threads per block
  // source: Number Conversion, Roland Backhouse 3/5/2001
  gpu->sim.TIlinearAtmCnt       = *ti_latm_cnt;
  gpu->sim.TIPaddedLinearAtmCnt = padded_latm_cnt;
  gpu->sim.pTIRegion            = gpu->pbTIRegion->_pDevData;
  gpu->sim.pTILinearAtmID       = gpu->pbTILinearAtmID->_pDevData;
  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gpu_set_ti_constants_: set lambda and 1-lambda for thermodynamic integration simulations, in
//                        both single and double precision.  This could be combined with
//                        gpu_setup_softcore_ti_ above.
//
// Arguments:
//   clambda:       the mixing factor to set
//   scalpha:       softcore vdw "fudge factor"
//   scbeta:        softcore eel "fudge factor"
//   ifmbar:        flag for performing MBAR calculation
//   bar_intervall: spacing between lambda values for MBAR. will soon be deprecated
//---------------------------------------------------------------------------------------------
extern "C" void gpu_set_ti_constants_(double* clambda, double* scalpha, double* scbeta,
                                      int* ifmbar, int* bar_intervall, int* bar_states,
                                      double bar_lambda[][2])
{
  PRINTMETHOD("gpu_set_ti_constants");
  gpu->sim.AFElambda[0]   = 1.0 - *clambda;
  gpu->sim.AFElambda[1]   = *clambda;
  gpu->sim.AFElambdaSP[0] = 1.0 - *clambda;
  gpu->sim.AFElambdaSP[1] = *clambda;
  gpu->sim.scalpha        = *scalpha;
  gpu->sim.scbeta         = *scbeta;
  gpu->sim.scalphaSP      = *scalpha;
  gpu->sim.scbetaSP       = *scbeta;
  gpu->sim.TIsigns[0]     = -1.0;
  gpu->sim.TIsigns[1]     = 1.0;
  gpu->sim.ifmbar         = *ifmbar;
  gpu->sim.bar_intervall  = *bar_intervall;
  gpu->sim.bar_states  = *bar_states;

  int bar_stride          = ((*bar_states + 63) >> 6) << 6;
  gpu->sim.bar_stride     = bar_stride;
  delete gpu->pbBarLambda;
  delete gpu->pbBarTot;
  gpu->pbBarLambda        = new GpuBuffer<double>(2 * bar_stride);
  gpu->pbBarTot           = new GpuBuffer<unsigned long long int>(bar_stride);
  for (int i = 0; i < *bar_states; i++) {
    gpu->pbBarLambda->_pSysData[i]              = bar_lambda[i][0];
    gpu->pbBarLambda->_pSysData[i + bar_stride] = bar_lambda[i][1];
  }
  int barsize             = gpu->sim.bar_states * gpu->sim.atoms;
  if (barsize % 4 != 0) {
    barsize              += barsize % 4;
  }
  gpu->pbBarLambda->Upload();
  gpu->pbBarTot->Upload();
  gpu->sim.pBarLambda     = gpu->pbBarLambda->_pDevData;
  gpu->sim.pBarTot        = gpu->pbBarTot->_pDevData;
  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gpu_ti_dynamic_lambda_(clambda): set lambda and 1-lambda following ntave steps.  Could
//                                  combine this with the above function.
//
// Arguments:
//   clambda:       the mixing factor to set
//---------------------------------------------------------------------------------------------
extern "C" void gpu_ti_dynamic_lambda_(double* clambda)
{
  PRINTMETHOD("gpu_ti_dynamic_lambda");
  gpu->sim.AFElambda[0]   = 1.0 - *clambda;
  gpu->sim.AFElambda[1]   = *clambda;
  gpu->sim.AFElambdaSP[0] = 1.0 - *clambda;
  gpu->sim.AFElambdaSP[1] = *clambda;
  gpuCopyConstants();

}

//---------------------------------------------------------------------------------------------
// gpu_steric_grid_: generate steric grid returns size of empty voxels
//
// Arguments:
//---------------------------------------------------------------------------------------------
extern "C" void gpu_steric_grid_(int *atm_cnt, int *maxXVxl, int *maxYVxl, int *maxZVxl,
                                 double *stericGSP, double steric_radius[], int *nmEmpVxls,
                                 int voxel_offset[6])
{
  gpu->sim.stericMaxXVxl = *maxXVxl;
  gpu->sim.stericMaxYVxl = *maxYVxl;
  gpu->sim.stericMaxZVxl = *maxZVxl;
  gpu->sim.stericGridSpacing = *stericGSP;
  delete gpu->pbStericRadius;
  delete gpu->pbVoxelOffset;
  delete gpu->pbStericGrid;
  delete gpu->pbStericThreadCount;
  delete gpu->pbStericEmptyVoxels;
  // 32 bits in int32
  int xdim=*maxXVxl/gpu->sim.stericBitCellWidth+1;
  int ydim=*maxYVxl/gpu->sim.stericBitCellHeight+1;
  int zdim=*maxZVxl/gpu->sim.stericBitCellDepth+1;
  gpu->pbStericGrid         = new GpuBuffer<unsigned int>(xdim*ydim*zdim);
  gpu->sim.pStericGrid     = gpu->pbStericGrid->_pDevData;
  for(int i=0; i < xdim*ydim*zdim; i++)
  {
    gpu->pbStericGrid->_pSysData[i]=0;
  }
  gpu->pbStericRadius       = new GpuBuffer<PMEFloat>(*atm_cnt);
  gpu->sim.pStericRadius   = gpu->pbStericRadius->_pDevData;
  for(int i=0; i < *atm_cnt; i++)
  {
     gpu->pbStericRadius->_pSysData[i]=steric_radius[i];
  }
  gpu->pbVoxelOffset        = new GpuBuffer<int>(6);
  gpu->sim.pVoxelOffset     = gpu->pbVoxelOffset->_pDevData;
  for(int i=0; i<6; i++)
  {
    gpu->pbVoxelOffset->_pSysData[i]=voxel_offset[i]-1; //Fortran indexes to C
  }
  gpu->pbStericThreadCount = new GpuBuffer<int>(gpu->PMENonbondForcesThreadsPerBlock*gpu->PMENonbondBlocks);
  gpu->sim.pStericThreadCount   = gpu->pbStericThreadCount->_pDevData;
  gpu->pbStericRadius->Upload();
  gpu->pbVoxelOffset->Upload();
  gpu->pbStericThreadCount->Upload();
  gpu->pbStericGrid->Upload();
  gpuCopyConstants();
  SetkCalculateCoarseGridEnergySim(gpu);
  kMakeStericGrid(gpu);
  cudaDeviceSynchronize();
  kCountStericGrid(gpu);
  gpu->pbStericThreadCount->Download();
  int emptyvoxels=0;
  for(int i=0; i< gpu->PMENonbondForcesThreadsPerBlock*gpu->PMENonbondBlocks; i++)
  {
    int stericcount = gpu->pbStericThreadCount->_pSysData[i];
    gpu->pbStericThreadCount->_pSysData[i] = emptyvoxels;
    emptyvoxels += stericcount;
  }
  gpu->pbStericThreadCount->Upload();
  *nmEmpVxls=emptyvoxels;
  gpu->pbStericEmptyVoxels         = new GpuBuffer<int>(*nmEmpVxls);
  gpu->sim.pStericEmptyVoxels      = gpu->pbStericEmptyVoxels->_pDevData;
  gpu->pbStericEmptyVoxels->Upload();
  gpuCopyConstants();
  SetkCalculateCoarseGridEnergySim(gpu);
}

//---------------------------------------------------------------------------------------------
// gpu_steric_grid_empty_voxels_: generate steric grid returns size of empty voxels
//
// Arguments:
//---------------------------------------------------------------------------------------------
extern "C" void  gpu_steric_grid_empty_voxels_(int *numvoxels, int listEmptyVoxels[])
{
  kGetEmptyStericGridVoxels(gpu);
  cudaDeviceSynchronize();
  cudaError_t status = cudaDeviceSynchronize();
  RTERROR(status, "kGetEmptyStericGridVoxels cudaDeviceSynchronize failed");
  gpu->pbStericEmptyVoxels->Download();
  for(int i=0; i < *numvoxels; i++)
  {
    listEmptyVoxels[i]=gpu->pbStericEmptyVoxels->_pSysData[i];
  }
}

//---------------------------------------------------------------------------------------------
// gpu_coarsegrid_setup_: setup coarse grid calculations
//
// Arguments:
//---------------------------------------------------------------------------------------------
extern "C" void gpu_coarsegrid_setup_(int *watatomcnt, double *ew_coeff, int *mcressize)
{
  gpu->sim.waterAtomCount = *watatomcnt;
  gpu->sim.mcressize = *mcressize;
#ifndef use_DPFP
  SetkCalculateCoarseGridEnergyERFC(*ew_coeff);
#endif
  delete gpu->pbCoarseGridEne;
  gpu->pbCoarseGridEne     = new GpuBuffer<PMEFloat>(*mcressize);
  gpu->sim.pCoarseGridEne  = gpu->pbCoarseGridEne->_pDevData;
  gpu->pbCoarseGridEne->Upload();
  gpuCopyConstants();
  kCalculateCoarseGridEnergyInitKernels(gpu);
}

//---------------------------------------------------------------------------------------------
// gpu_calc_coarsegrid_: calculates coarse grid energy
//
// Arguments:
//   coarseGrid: 4-D array with x,y,z and indicie for each atom
//   coarseGSP: coarsegrid spacing double
//   maxCoarseXVxl, maxCoarseYVxl, maxCoarseZVxl, length of x,y,z of coarse grid
//   mcressize: number of atoms going through coarse grid calculation
//   rndmVxlArray: array of which voxel a water atom is going to
//   rndmWatArray: selection of random waters
//   acceptance_array: output array for data
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calc_coarsegrid_(int coarseGrid[], double *coarseGSP, int *maxCoarseAtoms,
                                     int *maxCoarseXVxl, int *maxCoarseYVxl, int *maxCoarseZVxl,
                                     int *mcressize, int rndmWatArray[], double translateArray[][3],
                                     double acceptance_array[], int *is_orthog, int *ti_mode)
{
   gpu->sim.coarseMaxXvxl = *maxCoarseXVxl;
   gpu->sim.coarseMaxYvxl = *maxCoarseYVxl;
   gpu->sim.coarseMaxZvxl = *maxCoarseZVxl;
   gpu->sim.coarseGridSpacing = *coarseGSP;
   gpu->sim.coarseMaxAtoms = *maxCoarseAtoms;
   delete gpu->pbRandomWater;
   delete gpu->pbCoarseGrid;
   delete gpu->pbCoarseTrans;
   int GridSize = *maxCoarseAtoms * (*maxCoarseXVxl) * (*maxCoarseYVxl) * (*maxCoarseZVxl);
   gpu->pbRandomWater       = new GpuBuffer<int>(*mcressize);
   gpu->pbCoarseGrid        = new GpuBuffer<int>(GridSize);
   gpu->pbCoarseTrans       = new GpuBuffer<PMEFloat>(3 * *mcressize);
   gpu->sim.pRandomWater    = gpu->pbRandomWater->_pDevData;
   gpu->sim.pCoarseGrid     = gpu->pbCoarseGrid->_pDevData;
   gpu->sim.pCoarseTrans    = gpu->pbCoarseTrans->_pDevData;
   for(int i=0; i < *mcressize; i++)
   {
       gpu->pbCoarseTrans->_pSysData[i]=translateArray[i][0];
       gpu->pbCoarseTrans->_pSysData[*mcressize + i]=translateArray[i][1];
       gpu->pbCoarseTrans->_pSysData[*mcressize * 2 + i]=translateArray[i][2];
       gpu->pbRandomWater->_pSysData[i]=rndmWatArray[i];
   }
   int atomsize=int(*maxCoarseAtoms);
   for(int i=0; i < *maxCoarseXVxl; i++)
   {
       for(int j=0; j < *maxCoarseYVxl; j++)
       {
           for(int k=0; k < *maxCoarseZVxl; k++)
           {
                for(int l=0; l < atomsize; l++)
                {
                    // Not sure why these are padded by 3. atom size is padded by 2x
                    gpu->pbCoarseGrid->_pSysData[l + k * atomsize +
                                                 j * (*maxCoarseZVxl) * atomsize +
                                                 i * (*maxCoarseZVxl) * (*maxCoarseYVxl) * atomsize] =
                                                               coarseGrid[i+j*(*maxCoarseXVxl+3)+
                                                               k*(*maxCoarseYVxl+3)*(*maxCoarseXVxl+3)+
                                                               l*(*maxCoarseXVxl+3)*(*maxCoarseYVxl+3)*(*maxCoarseZVxl+3)];
                }
           }
       }
   }
   gpu->pbRandomWater->Upload();
   gpu->pbCoarseGrid->Upload();
   gpu->pbCoarseTrans->Upload();
   gpuCopyConstants();
   SetkCalculateCoarseGridEnergySim(gpu);
   if(*ti_mode>0)
   {
       if(*is_orthog>0)
       {
           kCalculateCoarseGridEnergyOrthogTI(gpu);
       }
       else
       {
           kCalculateCoarseGridEnergyNonOrthogTI(gpu);
       }
   }
   else
   {
       if(*is_orthog>0)
       {
           kCalculateCoarseGridEnergyOrthog(gpu);
       }
       else
       {
           kCalculateCoarseGridEnergyNonOrthog(gpu);
       }
   }
   cudaDeviceSynchronize();
   gpu->pbCoarseGridEne->Download();
   for(int i=0; i<*mcressize; i++)
   {
     acceptance_array[i]= gpu->pbCoarseGridEne->_pSysData[i];
   }
}

//---------------------------------------------------------------------------------------------
// gpu_upload_gamma_ln_: upload gamma_ln to the device. This is needed for neb call to fix the
//                       endpoint structures.
// Arguments:
//   gamma_ln: Langevin collision frequency
//---------------------------------------------------------------------------------------------
extern "C" void gpu_upload_gamma_ln_(double* gamma_ln)
{
  gpu->sim.gamma_ln                   = *gamma_ln;
  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gpu_upload_neb_info_: Upload nattgtrms, nattgtfit, rmsmask, fitmask for neb calculation.
//
// Arguments:
//   neb_nbead:    number of NEB replicas.
//   beadid:       NEB replica ID.
//   nattgtrms:    number of atoms in rmsmask.
//   nattgtfit:    number of atoms in fitmask.
//   last_neb_atm: last atom included in either rmsmask or fitmask.
//   skmin:        NEB min spring constant.
//   skmax:        NEB max spring constant.
//   tmode:        NEB tangent mode.
//   rmsmask:      atom selection mask used for NEB force application.
//   fitmask:      atom selection mask used for NEB structure fitting.
//   idx_mask:     combined ordered rmsmask and fitmask atoms.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_upload_neb_info_(int* neb_nbead, int* beadid, int* nattgtrms, int* nattgtfit,
                                     int* last_neb_atm, int* skmin, int* skmax, int* tmode,
                                     int rmsmask[], int fitmask[], int idx_mask[])
{
    PRINTMETHOD("gpu_upload_neb_info");
    gpu->sim.neb_nbead                  = *neb_nbead;
    gpu->sim.beadid                     = *beadid;
    gpu->sim.nattgtrms                  = *nattgtrms;
    gpu->sim.nattgtfit                  = *nattgtfit;
    gpu->sim.last_neb_atm               = *last_neb_atm;
    gpu->sim.skmin                      = *skmin;
    gpu->sim.skmax                      = *skmax;
    gpu->sim.tmode                      = *tmode;

    delete gpu->pbRMSMask;
    delete gpu->pbFitMask;
    delete gpu->pbAtmIdx;
    int RMSStride = ((gpu->sim.nattgtrms + gpu->sim.grid - 1) >>
                        gpu->sim.gridBits) << gpu->sim.gridBits;
    int FitStride = ((gpu->sim.nattgtfit + gpu->sim.grid - 1) >>
                        gpu->sim.gridBits) << gpu->sim.gridBits;
    gpu->pbRMSMask = new GpuBuffer<int>(RMSStride);
    gpu->pbFitMask = new GpuBuffer<int>(FitStride);
    gpu->pbAtmIdx = new GpuBuffer<int>(get_shuttle_size(1, sizeof(int),
                                                   gpu->sim.nShuttle));

    gpu->sim.pRMSMask = gpu->pbRMSMask->_pDevData;
    gpu->sim.pFitMask = gpu->pbFitMask->_pDevData;
    gpu->sim.pAtmIdx = gpu->pbAtmIdx->_pDevData;

    for (int i = 0; i < gpu->sim.nattgtrms; i++){
      gpu->pbRMSMask->_pSysData[i] = rmsmask[i] - 1;
    }

    for (int i = 0; i < gpu->sim.nattgtfit; i++){
      gpu->pbFitMask->_pSysData[i] = fitmask[i] - 1;
    }

    for (int i = 0; i < gpu->sim.nShuttle; i++){
      gpu->pbAtmIdx->_pSysData[i] = idx_mask[i];
    }

    gpu->pbRMSMask->Upload();
    gpu->pbFitMask->Upload();
    gpu->pbAtmIdx->Upload();
    gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gpu_setup_neb_: Setup for NEB, called from pmemd.F90.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_setup_neb_()
{
    PRINTMETHOD("gpu_setup_neb");
    delete gpu->pbNEBEnergyAll;
    delete gpu->pbTangents;
    delete gpu->pbSpringForce;
    delete gpu->pbNEBForce;
    delete gpu->pbNextDataShuttle;
    delete gpu->pbPrevDataShuttle;
    delete gpu->pbKabschCOM;
    delete gpu->pbDataSPR;
    delete gpu->pbRotAtm;
    delete gpu->pbtotFitMass;

    int NEBStride = ((gpu->sim.neb_nbead + gpu->sim.grid - 1) >>
                        gpu->sim.gridBits) << gpu->sim.gridBits;
    gpu->pbNEBEnergyAll = new GpuBuffer<PMEDouble>(NEBStride);
    gpu->pbTangents = new GpuBuffer<double>(get_shuttle_size(3, sizeof(double),
                                                   gpu->sim.nShuttle));
    gpu->pbSpringForce = new GpuBuffer<double>(get_shuttle_size(3, sizeof(double),
                                                   gpu->sim.nShuttle));
    gpu->pbNEBForce = new GpuBuffer<double>(get_shuttle_size(3, sizeof(double),
                                                   gpu->sim.nShuttle));
    if (gpu->bCanMapHostMemory) {
      gpu->pbNextDataShuttle = new GpuBuffer<double>(get_shuttle_size(3, sizeof(double),
                                                   gpu->sim.nShuttle), false, true);
      gpu->pbPrevDataShuttle = new GpuBuffer<double>(get_shuttle_size(3, sizeof(double),
                                                   gpu->sim.nShuttle), false, true);
    }
    else {
      gpu->pbNextDataShuttle = new GpuBuffer<double>(get_shuttle_size(3, sizeof(double),
                                                   gpu->sim.nShuttle));
      gpu->pbPrevDataShuttle = new GpuBuffer<double>(get_shuttle_size(3, sizeof(double),
                                                   gpu->sim.nShuttle));
    }
    gpu->pbKabschCOM = new GpuBuffer<double>(27);
    if (gpu->bCanMapHostMemory) {
      gpu->pbDataSPR = new GpuBuffer<double>(5, false, true);
    }
    else {
      gpu->pbDataSPR = new GpuBuffer<double>(5);
    }
    if (gpu->bCanMapHostMemory) {
      gpu->pbRotAtm = new GpuBuffer<double>(9, false, true);
    }
    else {
      gpu->pbRotAtm = new GpuBuffer<double>(9);
    }
    gpu->pbtotFitMass = new GpuBuffer<double>(1);

    gpu->sim.pNEBEnergyAll = gpu->pbNEBEnergyAll->_pDevData;
    gpu->sim.pTangents = gpu->pbTangents->_pDevData;
    gpu->sim.pSpringForce = gpu->pbSpringForce->_pDevData;
    gpu->sim.pNEBForce = gpu->pbNEBForce->_pDevData;
    gpu->sim.pNextDataShuttle = gpu->pbNextDataShuttle->_pDevData;
    gpu->sim.pPrevDataShuttle = gpu->pbPrevDataShuttle->_pDevData;
    gpu->sim.pKabschCOM = gpu->pbKabschCOM->_pDevData;
    gpu->sim.pPrevKabsch = gpu->pbKabschCOM->_pDevData + 0;
    gpu->sim.pNextKabsch = gpu->pbKabschCOM->_pDevData + 9;
    gpu->sim.pSelfCOM = gpu->pbKabschCOM->_pDevData + 18;
    gpu->sim.pPrevCOM = gpu->pbKabschCOM->_pDevData + 21;
    gpu->sim.pNextCOM = gpu->pbKabschCOM->_pDevData + 24;
    gpu->sim.pDataSPR = gpu->pbDataSPR->_pDevData;
    gpu->sim.pNorm = gpu->pbDataSPR->_pDevData + 0;
    gpu->sim.pDotProduct1 = gpu->pbDataSPR->_pDevData + 1;
    gpu->sim.pDotProduct2 = gpu->pbDataSPR->_pDevData + 2;
    gpu->sim.pSpring1 = gpu->pbDataSPR->_pDevData + 3;
    gpu->sim.pSpring2 = gpu->pbDataSPR->_pDevData + 4;
    gpu->sim.ptotFitMass = gpu->pbtotFitMass->_pDevData;
    gpu->sim.pRotAtm = gpu->pbRotAtm->_pDevData;

    gpuCopyConstants();
}

//--------------------------------------------------------------------------------------------
// gpu_neb_exchange_crd_: NEB routine for exchanging the neighboring replicas coordinates
//                        through the shuttle transfer.
//--------------------------------------------------------------------------------------------
extern "C" void gpu_neb_exchange_crd_()
{
  PRINTMETHOD("gpu_neb_exchange_crd");
  int buff_size = get_shuttle_size(3, sizeof(double),gpu->sim.nShuttle);
  kNEBSendRecv(gpu, buff_size);
}

//--------------------------------------------------------------------------------------------
// gpu_report_neb_ene_: NEB routine for the replicas' energy print-out.
//
// Arguments:
//   master_size:
//   neb_nrg_all:
//--------------------------------------------------------------------------------------------
extern "C" void gpu_report_neb_ene_(int* master_size, double neb_nrg_all[])
{
  PRINTMETHOD("gpu_report_neb_ene");
  NEB_report_energy(gpu, master_size, neb_nrg_all);
}

//--------------------------------------------------------------------------------------------
// gpu_set_neb_springs_: NEB routine for setting the spring constants.
//--------------------------------------------------------------------------------------------
extern "C" void gpu_set_neb_springs_()
{
  PRINTMETHOD("gpu_set_neb_springs");
  kNEBspr(gpu);
}

//--------------------------------------------------------------------------------------------
// gpu_calculate_neb_frc_: Routine for calculating the NEB forces every nebfreq steps.
//
// Arguments:
//   neb_force: NEB forces.
//--------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_neb_frc_(double neb_force[][3])
{
  PRINTMETHOD("gpu_calculate_neb_frc");
  kNEBfrc(gpu, neb_force);
}

//--------------------------------------------------------------------------------------------
// gpu_calculate_neb_frc_nstep_: Routine for calculating the NEB forces at the off steps.
//
// Arguments:
//   neb_force: NEB forces.
//--------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_neb_frc_nstep_(double neb_force[][3])
{
  PRINTMETHOD("gpu_calculate_neb_frc_nstep");
  kNEBfrc_nstep(gpu, neb_force);
}

//--------------------------------------------------------------------------------------------
// gpu_neb_rmsfit_: Routine for calculating the NEB structure fitting.
//--------------------------------------------------------------------------------------------
extern "C" void gpu_neb_rmsfit_()
{
  PRINTMETHOD("gpu_neb_rmsfit");
  kFitCOM(gpu);
}

//---------------------------------------------------------------------------------------------
// gpu_setup_shuttle_info_: Upload atom_list and nshatom for partial coordinates download.
//
// Arguments:
//   nshatom:   Number of partial coordinates to be downloaded.
//   infocode:  Determines shuttle type (PMEFloat or integer).
//   atom_list: List of atom numbers for partial coordinates download.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_setup_shuttle_info_(int* nshatom, int* infocode, int atom_list[])
{
  PRINTMETHOD("gpu_setup_shuttle_info");
  gpu->sim.nShuttle                   = *nshatom;
  gpu->sim.ShuttleType                = *infocode;
  if (gpu->sim.ShuttleType == 0) {
    delete gpu->pbDataShuttle;
    delete gpu->pbShuttleTickets;
    gpu->pbDataShuttle = new GpuBuffer<double>(get_shuttle_size(3, sizeof(double),
                                                                gpu->sim.nShuttle),
					       true, true);
    gpu->pbShuttleTickets = new GpuBuffer<int>(get_shuttle_size(1, sizeof(int),
								gpu->sim.nShuttle),
                                               true, true);
  }
  else {

    // Print an error for now
    printf("| Error: this functionality (ShuttleType == 1) is not yet available.\n");
  }
  gpu->sim.pDataShuttle = gpu->pbDataShuttle->_pDevData;
  gpu->sim.pShuttleTickets = gpu->pbShuttleTickets->_pDevData;
  for (int i = 0; i < gpu->sim.nShuttle; i++){
    gpu->pbShuttleTickets->_pSysData[i] = atom_list[i] - 1;
  }
  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gpu_shuttle_retrieve_data_: shuttle download of the data.
//
// Arguments:
//   atm_data: the atomic array for the data
//   int* flag: type of shuttle download
//---------------------------------------------------------------------------------------------
extern "C" void gpu_shuttle_retrieve_data_(double atm_data[][3], int* flag)
{
  PRINTMETHOD("gpu_shuttle_retrieve_data");
  if (gpu->sim.nShuttle > 0) {
    kRetrieveSimData(gpu, atm_data, *flag);
    cudaError_t status = cudaDeviceSynchronize();
    RTERROR(status, "gpu_shuttle_retrieve_data cudaDeviceSynchronize failed");
  }
}

//---------------------------------------------------------------------------------------------
// gpu_shuttle_post_data_: shuttle upload of the data.
//
// Arguments:
//   atm_data: the atomic data
//   int* flag: type of shuttle upload
//---------------------------------------------------------------------------------------------
extern "C" void gpu_shuttle_post_data_(double atm_data[][3], int* flag)
{
  PRINTMETHOD("gpu_shuttle_post_data");
  if (gpu->sim.nShuttle > 0) {
    kPostSimData(gpu, atm_data, *flag);
    cudaError_t status = cudaDeviceSynchronize();
    RTERROR(status, "gpu_shuttle_post_data cudaDeviceSynchronize failed");
  }
}

//---------------------------------------------------------------------------------------------
// gpu_upload_crd_: upload coordinates to the device.  This routine is called by all Fortran
//                  routines seeking to interface with the CUDA code.
//
// Arguments:
//    atm_crd:  the atomic coordinates
//---------------------------------------------------------------------------------------------
extern "C" void gpu_upload_crd_(double atm_crd[][3])
{
  PRINTMETHOD("gpu_upload_crd");
  if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
    if (gpu->bNewNeighborList) {
      gpu->pbImageIndex->Download();
      gpu->bNewNeighborList = false;
    }
    gpu->pbImage->Download();
    unsigned int* pImageAtomLookup = &(gpu->pbImageIndex->_pSysData[gpu->sim.imageStride * 2]);
    double *pCrd = gpu->pbImage->_pSysData;
    if (gpu->sim.pImageX != gpu->pbImage->_pDevData) {
      pCrd = gpu->pbImage->_pSysData + gpu->sim.stride3;
    }
    for (int i = 0; i < gpu->sim.atoms; i++) {
      int i1 = pImageAtomLookup[i];
      pCrd[i1] = atm_crd[i][0];
      pCrd[i1 + gpu->sim.stride] = atm_crd[i][1];
      pCrd[i1 + gpu->sim.stride2] = atm_crd[i][2];
    }
    gpu->pbImage->Upload();
  }
  else {
    for (int i = 0; i < gpu->sim.atoms; i++) {
      gpu->pbAtom->_pSysData[i] = atm_crd[i][0];
      gpu->pbAtom->_pSysData[i + gpu->sim.stride] = atm_crd[i][1];
      gpu->pbAtom->_pSysData[i + gpu->sim.stride2] = atm_crd[i][2];
      gpu->pbAtomXYSP->_pSysData[i].x = atm_crd[i][0];
      gpu->pbAtomXYSP->_pSysData[i].y = atm_crd[i][1];
      gpu->pbAtomZSP->_pSysData[i] = atm_crd[i][2];
    }
    for (int i = gpu->sim.atoms; i < gpu->sim.stride; i++) {
      gpu->pbAtom->_pSysData[i] = 9999990000.0 + i*2000.0;
      gpu->pbAtom->_pSysData[i + gpu->sim.stride] = 9999990000.0 + i*2000.0;
      gpu->pbAtom->_pSysData[i + gpu->sim.stride2] = 9999990000.0 + i*2000.0;
      gpu->pbAtomXYSP->_pSysData[i].x = 9999990000.0f + i*2000.0;
      gpu->pbAtomXYSP->_pSysData[i].y = 9999990000.0f + i*2000.0;
      gpu->pbAtomZSP->_pSysData[i] = 9999990000.0f + i*2000.0;
    }
    gpu->pbAtom->Upload();
    gpu->pbAtomXYSP->Upload();
    gpu->pbAtomZSP->Upload();
  }
}

//---------------------------------------------------------------------------------------------
// gpu_upload_crd_gb_cph_: upload coordinates to the device in the context of constant PH
//                         simulations or GB implicit solvent calculations.  This routine is
//                         called in addition to the standard gpu_upload_crd_ in constantph.F90
//                         of the pmemd Fortran code.
//
// Arguments:
//   atm_crd:  the atomic coordinates
//---------------------------------------------------------------------------------------------
extern "C" void gpu_upload_crd_gb_cph_(double atm_crd[][3])
{
  PRINTMETHOD("gpu_upload_crd_gb_cph");
  for (int i = 0; i < gpu->sim.atoms; i++) {
    gpu->pbAtom->_pSysData[i] = atm_crd[i][0];
    gpu->pbAtom->_pSysData[i + gpu->sim.stride] = atm_crd[i][1];
    gpu->pbAtom->_pSysData[i + gpu->sim.stride2] = atm_crd[i][2];
    gpu->pbAtomXYSP->_pSysData[i].x = atm_crd[i][0];
    gpu->pbAtomXYSP->_pSysData[i].y = atm_crd[i][1];
    gpu->pbAtomZSP->_pSysData[i] = atm_crd[i][2];
  }
  for (int i = gpu->sim.atoms; i < gpu->sim.stride; i++) {
    gpu->pbAtom->_pSysData[i] = 9999990000.0 + i * 2000.0;
    gpu->pbAtom->_pSysData[i + gpu->sim.stride] = 9999990000.0 + i * 2000.0;
    gpu->pbAtom->_pSysData[i + gpu->sim.stride2] = 9999990000.0 + i * 2000.0;
    gpu->pbAtomXYSP->_pSysData[i].x = 9999990000.0f + i * 2000.0;
    gpu->pbAtomXYSP->_pSysData[i].y = 9999990000.0f + i * 2000.0;
    gpu->pbAtomZSP->_pSysData[i] = 9999990000.0f + i * 2000.0;
  }
  gpu->pbAtom->Upload();
  gpu->pbAtomXYSP->Upload();
  gpu->pbAtomZSP->Upload();
}

//---------------------------------------------------------------------------------------------
// gpu_download_crd_: download coordinates from the device.  This routine is called by all
//                    routines seeking to interface with the GPU code.
//
// Arguments:
//   atm_crd:  the atomic coordinates
//---------------------------------------------------------------------------------------------
extern "C" void gpu_download_crd_(double atm_crd[][3])
{
  PRINTMETHOD("gpu_download_crd");
  if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
    if (gpu->bNewNeighborList) {
      gpu->pbImageIndex->Download();
      gpu->bNewNeighborList = false;
    }
    gpu->pbImage->Download();
    unsigned int* pImageAtomLookup = &(gpu->pbImageIndex->_pSysData[gpu->sim.imageStride * 2]);
    double *pCrd = gpu->pbImage->_pSysData;
    if (gpu->sim.pImageX != gpu->pbImage->_pDevData) {
      pCrd = gpu->pbImage->_pSysData + gpu->sim.stride3;
    }
    int i;
    for (i = 0; i < gpu->sim.atoms; i++) {
      int i1 = pImageAtomLookup[i];
      atm_crd[i][0] = pCrd[i1];
      atm_crd[i][1] = pCrd[i1 + gpu->sim.stride];
      atm_crd[i][2] = pCrd[i1 + gpu->sim.stride2];
    }
  }
  else {
    gpu->pbAtom->Download();
    for (int i = 0; i < gpu->sim.atoms; i++) {
      atm_crd[i][0] = gpu->pbAtom->_pSysData[i];
      atm_crd[i][1] = gpu->pbAtom->_pSysData[i + gpu->sim.stride];
      atm_crd[i][2] = gpu->pbAtom->_pSysData[i + gpu->sim.stride2];
    }
  }
}

//-----------------------------
// gpu_upload_iphmd_: upload phmd flag.
//
// Arguments:
//   iphmd:   phmd flag
//---------------------------------------------------------------------------------------------
extern "C" void gpu_upload_iphmd_(int* iphmd)
{
  gpu->iphmd = *iphmd;
}

//---------------------------------------------------------------------------------------------
// gpu_reset_park_: upload phmd flag.
//
// Arguments:
//   park:         The pH differences for the residues.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_reset_park_(double park[])
{
  for(int i = 0; i < gpu->sim.ntitr; i++) {
    gpu->park->_pSysData[i] = park[i];
  }
  gpu->park->Upload();
}

//---------------------------------------------------------------------------------------------
// gpu_upload_phmd_: upload phmd information.
//
// Arguments:
//   qsetlam:      Is lambda fixed?
//   ntitr:        The number of titration variables.
//   qmass_phmd:   The mass of the lambda particles.
//   temp_phmd:    The temperature of the lambda particles.
//   phbeta:       The friction constant for the lambda particle integrator.
//   para:         The A parameters for the residues.
//   parb:         The B parameters for the residues.
//   park:         The pH differences for the residues.
//   barr:         The barrier potentials.
//   ph_theta:     The titration coordinates.
//   dph_theta:    The forces on the titration coordinates.
//   vph_theta:    The velocities of the titration coordinates.
//   charge_phmd:  The charges of the atoms.
//   grplist:      The titratable residue each atom belongs to.
//   sp_grp:       What type of titratable residue each residue is.
//   sp_par:       The special paramters for each titratable residue.
//   qstate1:      The charges on the atoms when lambda=1.
//   qstate2:      The charges on the atoms when lambda=0.
//   vstate1:      The vdw states of the atoms when lambda=1.
//   vstate2:      The vdw states of the atoms when lambda=0.
//   qstate1_md:   The charges on the atoms in MD when lambda=1.
//   qstate2_md:   The charges on the atoms in MD when lambda=0.
//   parmod:       The six extra paramters for ASP and GLU residues.
//   resbounds1:   The atoms beginning the titratable residues.
//   resbounds2:   The atoms ending the titratable residues.
//   baseres:      Is a residue the root titratable residue?
//   max_numch:    The maximum number of charges in a titratable residue.
//   excludeatoms: Identifies atoms excluded from the Born radii calculations.
//   iphmd:        What version of phmd to run.
//   solvph:       Solution pH.
//   titratoms:    An array of the indices of those atoms that are titrating
//   ntitratoms:   The number of titrating atoms
//   numex:        Exclusion bounds
//   natex:        Exclusion list
//   linkgrplist:  List of residues linked to each titratable residue (for PME)
//   use_linked:   Use linked water?
//   use_forces:   Use the forces in the motion of the titration variables?
//---------------------------------------------------------------------------------------------
extern "C" void gpu_upload_phmd_(bool* qsetlam, int* ntitr, double* qmass_phmd,
double* temp_phmd, double* phbeta, double para[], double parb[], double park[],
double barr[], double ph_theta[], double dph_theta[], double vph_theta[],
double charge_phmd[], int grplist[], int sp_grp[], double sp_par[][2],
double qstate1[][2], double qstate2[][2], double vstate1[][2],
double vstate2[][2], double qstate1_md[][2], double qstate2_md[][2],
double parmod[][6], int resbounds1[], int resbounds2[], int baseres[],
int* max_numch, int excludeatoms[], int* iphmd,
double* solvph, int titratoms[],
int* ntitratoms, int numex[], int natex[], int linkgrplist[],
bool* use_linked, bool* use_forces)
{
  gpu->use_linked = *use_linked;
  gpu->iphmd = *iphmd;
  gpu->sim.iphmd = *iphmd;
  gpu->sim.qsetlam = *qsetlam;
  gpu->sim.use_forces = *use_forces;
  gpu->sim.ntitr = *ntitr;
  gpu->sim.qmass_phmd = *qmass_phmd;
  gpu->sim.temp_phmd = *temp_phmd;
  gpu->sim.phbeta = *phbeta;
  gpu->sim.solvph = *solvph;
  gpu->sim.phmd_padded_ntitr = ((*ntitr + 31) >> 5) << 5;
  gpu->sim.padded_numch = ((*max_numch + 31) >> 5) << 5;
  gpu->sim.randomNumbersPHMD = gpu->sim.phmd_padded_ntitr * gpu->sim.randomSteps;
  gpu->sim.ntitratoms = *ntitratoms;
  gpu->sim.phmd_stride = (*ntitratoms + gpu->PMENonbondBlocks * gpu->PMENonbondForcesThreadsPerBlock - 1) /
                         (gpu->PMENonbondBlocks * gpu->PMENonbondForcesThreadsPerBlock);
  if (gpu->sim.atoms > gpu->PMENonbondBlocks * gpu->PMENonbondForcesThreadsPerBlock) {
    gpu->sim.phmd_atom_stride = gpu->sim.atoms /
                           (gpu->PMENonbondBlocks * gpu->PMENonbondForcesThreadsPerBlock);
  }
  else {
    gpu->sim.phmd_atom_stride = 1;
  }
  delete gpu->para;
  delete gpu->parb;
  delete gpu->park;
  delete gpu->barr;
  delete gpu->ph_theta;
  delete gpu->dph_theta;
  delete gpu->vph_theta;
  delete gpu->dphold;
  delete gpu->vphold;
  delete gpu->charge_phmd;
  delete gpu->grplist;
  delete gpu->linkgrplist;
  delete gpu->ImageGrplist;
  delete gpu->ImageIndexPHMD;
  delete gpu->sp_grp;
  delete gpu->sp_par;
  delete gpu->qstate1;
  delete gpu->qstate2;
  delete gpu->ImageQstate1;
  delete gpu->ImageQstate2;
  delete gpu->vstate1;
  delete gpu->vstate2;
  delete gpu->ImageVstate1;
  delete gpu->ImageVstate2;
  delete gpu->qstate1_md;
  delete gpu->qstate2_md;
  delete gpu->parmod;
  delete gpu->pbRandomPHMD;
  delete gpu->baseres;
  delete gpu->excludeatoms;
  delete gpu->dph_accumulator;
  delete gpu->dph_plus_accumulator;
  delete gpu->phmd_ene;
  delete gpu->titratoms;
  delete gpu->dudls;
  //Code borrowed for exclusions
  delete gpu->pbNLExclusionStartCount;
  gpu->parb = new GpuBuffer<PMEFloat>(gpu->sim.phmd_padded_ntitr);
  gpu->para = new GpuBuffer<PMEFloat>(gpu->sim.phmd_padded_ntitr);
  gpu->park = new GpuBuffer<PMEFloat>(gpu->sim.phmd_padded_ntitr);
  gpu->barr = new GpuBuffer<PMEFloat>(gpu->sim.phmd_padded_ntitr);
  gpu->ph_theta = new GpuBuffer<PMEFloat>(gpu->sim.phmd_padded_ntitr);
  gpu->dph_theta = new GpuBuffer<PMEFloat>(gpu->sim.phmd_padded_ntitr);
  gpu->vph_theta = new GpuBuffer<PMEFloat>(gpu->sim.phmd_padded_ntitr);
  gpu->dphold = new GpuBuffer<PMEFloat>(gpu->sim.phmd_padded_ntitr);
  gpu->vphold = new GpuBuffer<PMEFloat>(gpu->sim.phmd_padded_ntitr);
  gpu->charge_phmd = new GpuBuffer<PMEFloat>(gpu->sim.paddedNumberOfAtoms);
  gpu->ImageCharge_phmd = new GpuBuffer<PMEFloat>(2 * gpu->sim.paddedNumberOfAtoms);
  gpu->grplist = new GpuBuffer <int>(gpu->sim.paddedNumberOfAtoms);
  gpu->linkgrplist = new GpuBuffer <int>(gpu->sim.paddedNumberOfAtoms);
  gpu->ImageGrplist = new GpuBuffer <int>(2 * gpu->sim.paddedNumberOfAtoms);
  gpu->ImageIndexPHMD = new GpuBuffer<int>(2 * gpu->sim.paddedNumberOfAtoms);
  gpu->sp_grp = new GpuBuffer<int>(gpu->sim.phmd_padded_ntitr);
  gpu->sp_par = new GpuBuffer<PMEFloat2>(gpu->sim.phmd_padded_ntitr);
  gpu->qstate1 = new GpuBuffer<PMEFloat2>(gpu->sim.paddedNumberOfAtoms);
  gpu->qstate2 = new GpuBuffer<PMEFloat2>(gpu->sim.paddedNumberOfAtoms);
  gpu->ImageQstate1 = new GpuBuffer<PMEFloat2>(2 * gpu->sim.paddedNumberOfAtoms);
  gpu->ImageQstate2 = new GpuBuffer<PMEFloat2>(2 * gpu->sim.paddedNumberOfAtoms);
  gpu->vstate1 = new GpuBuffer<PMEFloat2>(gpu->sim.paddedNumberOfAtoms);
  gpu->vstate2 = new GpuBuffer<PMEFloat2>(gpu->sim.paddedNumberOfAtoms);
  gpu->ImageVstate1 = new GpuBuffer<PMEFloat2>(2 * gpu->sim.paddedNumberOfAtoms);
  gpu->ImageVstate2 = new GpuBuffer<PMEFloat2>(2 * gpu->sim.paddedNumberOfAtoms);
  gpu->qstate1_md = new GpuBuffer<PMEFloat2>(gpu->sim.paddedNumberOfAtoms);
  gpu->qstate2_md = new GpuBuffer<PMEFloat2>(gpu->sim.paddedNumberOfAtoms);
  gpu->parmod = new GpuBuffer<PMEFloat>(6*gpu->sim.phmd_padded_ntitr);
  gpu->pbRandomPHMD = new GpuBuffer<double>(gpu->sim.randomNumbersPHMD);
  gpu->resbounds = new GpuBuffer<int2>(gpu->sim.phmd_padded_ntitr);
  gpu->baseres = new GpuBuffer<bool>(gpu->sim.phmd_padded_ntitr);
  gpu->excludeatoms = new GpuBuffer<int>(gpu->sim.paddedNumberOfAtoms);
  gpu->dph_accumulator = new GpuBuffer<PMEAccumulator>(gpu->sim.paddedNumberOfAtoms);
  gpu->dph_plus_accumulator = new GpuBuffer<PMEAccumulator>(gpu->sim.paddedNumberOfAtoms);
  gpu->titratoms = new GpuBuffer<int>(*ntitratoms);
  gpu->dudls = new GpuBuffer<PMEFloat>(gpu->sim.phmd_padded_ntitr);
  for(int i = 0; i < *ntitratoms; i++) {
    gpu->titratoms->_pSysData[i] = titratoms[i];
  }
  gpu->titratoms->Upload();
  gpu->sim.ptitratoms = gpu->titratoms->_pDevData;
  gpu->phmd_ene = new GpuBuffer<unsigned long long int>(1);
  gpu->phmd_ene->_pSysData[0] = 0;
  for(int i = 0; i < *ntitr; i++) {
    gpu->para->_pSysData[i] = para[i];
    gpu->parb->_pSysData[i] = parb[i];
    gpu->park->_pSysData[i] = park[i];
    gpu->barr->_pSysData[i] = barr[i];
    gpu->ph_theta->_pSysData[i] = ph_theta[i];
    gpu->dph_theta->_pSysData[i] = dph_theta[i];
    gpu->vph_theta->_pSysData[i] = vph_theta[i];
    gpu->dphold->_pSysData[i] = dph_theta[i];
    gpu->vphold->_pSysData[i] = vph_theta[i];
    gpu->sp_grp->_pSysData[i] = sp_grp[i];
    gpu->sp_par->_pSysData[i].x = sp_par[i][0];
    gpu->sp_par->_pSysData[i].y = sp_par[i][1];
    gpu->resbounds->_pSysData[i].x = resbounds1[i];
    gpu->resbounds->_pSysData[i].y = resbounds2[i];
    gpu->baseres->_pSysData[i] = baseres[i];
    gpu->dudls->_pSysData[i] = 0;
    for(int j = 0; j < 6; j++ ) {
      gpu->parmod->_pSysData[i*6+j] = parmod[i][j];
    }
  }
  for(int i = *ntitr; i < gpu->sim.phmd_padded_ntitr; i++) {
    gpu->para->_pSysData[i] = (PMEFloat)0.0;
    gpu->parb->_pSysData[i] = (PMEFloat)0.0;
    gpu->park->_pSysData[i] = (PMEFloat)0.0;
    gpu->barr->_pSysData[i] = (PMEFloat)0.0;
    gpu->ph_theta->_pSysData[i] = (PMEFloat)0.0;
    gpu->dph_theta->_pSysData[i] = (PMEFloat)0.0;
    gpu->vph_theta->_pSysData[i] = (PMEFloat)0.0;
    gpu->dphold->_pSysData[i] = (PMEFloat)0.0;
    gpu->vphold->_pSysData[i] = (PMEFloat)0.0;
    gpu->sp_grp->_pSysData[i] = 0;
    gpu->sp_par->_pSysData[i].x = (PMEFloat)0.0;
    gpu->sp_par->_pSysData[i].y = (PMEFloat)0.0;
    gpu->resbounds->_pSysData[i].x = 0;
    gpu->resbounds->_pSysData[i].y = 0;
    gpu->baseres->_pSysData[i] = false;
    gpu->dudls->_pSysData[i] = 0;
    for(int j = 0; j < 6; j++ ) {
      gpu->parmod->_pSysData[i*6+j] = (PMEFloat)0.0;
    }
  }
  gpu->sumq = 0;
  gpu->sumq2 = 0;
  for(int i = 0; i < gpu->sim.atoms; i++) {
    gpu->sumq += charge_phmd[i];
    gpu->sumq2 += charge_phmd[i] * charge_phmd[i];
    gpu->charge_phmd->_pSysData[i] = charge_phmd[i];
    gpu->ImageCharge_phmd->_pSysData[i] = charge_phmd[i];
    gpu->grplist->_pSysData[i] = grplist[i];
    gpu->ImageGrplist->_pSysData[i] = grplist[i];
    gpu->ImageIndexPHMD->_pSysData[i] = i;
    gpu->qstate1->_pSysData[i].x = qstate1[i][0];
    gpu->qstate1->_pSysData[i].y = qstate1[i][1];
    gpu->qstate2->_pSysData[i].x = qstate2[i][0];
    gpu->qstate2->_pSysData[i].y = qstate2[i][1];
    gpu->ImageQstate1->_pSysData[i].x = qstate1[i][0];
    gpu->ImageQstate1->_pSysData[i].y = qstate1[i][1];
    gpu->ImageQstate2->_pSysData[i].x = qstate2[i][0];
    gpu->ImageQstate2->_pSysData[i].y = qstate2[i][1];
    gpu->vstate1->_pSysData[i].x = vstate1[i][0];
    gpu->vstate1->_pSysData[i].y = vstate1[i][1];
    gpu->vstate2->_pSysData[i].x = vstate2[i][0];
    gpu->vstate2->_pSysData[i].y = vstate2[i][1];
    gpu->ImageVstate1->_pSysData[i].x = vstate1[i][0];
    gpu->ImageVstate1->_pSysData[i].y = vstate1[i][1];
    gpu->ImageVstate2->_pSysData[i].x = vstate2[i][0];
    gpu->ImageVstate2->_pSysData[i].y = vstate2[i][1];
    gpu->qstate1_md->_pSysData[i].x = qstate1_md[i][0];
    gpu->qstate1_md->_pSysData[i].y = qstate1_md[i][1];
    gpu->qstate2_md->_pSysData[i].x = qstate2_md[i][0];
    gpu->qstate2_md->_pSysData[i].y = qstate2_md[i][1];
    gpu->excludeatoms->_pSysData[i] = excludeatoms[i];
    gpu->dph_accumulator->_pSysData[i] = 0;
    gpu->dph_plus_accumulator->_pSysData[i] = 0;
    gpu->linkgrplist->_pSysData[i] = linkgrplist[i];
  }
  for(int i = gpu->sim.atoms; i < gpu->sim.paddedNumberOfAtoms; i++) {
    gpu->charge_phmd->_pSysData[i] = (PMEFloat)0.0;
    gpu->ImageCharge_phmd->_pSysData[i] = (PMEFloat)0.0;
    gpu->grplist->_pSysData[i] = 0;
    gpu->ImageGrplist->_pSysData[i] = 0;
    gpu->ImageIndexPHMD->_pSysData[i] = 0;
    gpu->qstate1->_pSysData[i].x = (PMEFloat)0.0;
    gpu->qstate1->_pSysData[i].y = (PMEFloat)0.0;
    gpu->qstate2->_pSysData[i].x = (PMEFloat)0.0;
    gpu->qstate2->_pSysData[i].y = (PMEFloat)0.0;
    gpu->ImageQstate1->_pSysData[i].x = (PMEFloat)0.0;
    gpu->ImageQstate1->_pSysData[i].y = (PMEFloat)0.0;
    gpu->ImageQstate2->_pSysData[i].x = (PMEFloat)0.0;
    gpu->ImageQstate2->_pSysData[i].y = (PMEFloat)0.0;
    gpu->vstate1->_pSysData[i].x = (PMEFloat)0.0;
    gpu->vstate1->_pSysData[i].y = (PMEFloat)0.0;
    gpu->vstate2->_pSysData[i].x = (PMEFloat)0.0;
    gpu->vstate2->_pSysData[i].y = (PMEFloat)0.0;
    gpu->ImageVstate1->_pSysData[i].x = (PMEFloat)0.0;
    gpu->ImageVstate1->_pSysData[i].y = (PMEFloat)0.0;
    gpu->ImageVstate2->_pSysData[i].x = (PMEFloat)0.0;
    gpu->ImageVstate2->_pSysData[i].y = (PMEFloat)0.0;
    gpu->qstate1_md->_pSysData[i].x = (PMEFloat)0.0;
    gpu->qstate1_md->_pSysData[i].y = (PMEFloat)0.0;
    gpu->qstate2_md->_pSysData[i].x = (PMEFloat)0.0;
    gpu->qstate2_md->_pSysData[i].y = (PMEFloat)0.0;
    gpu->excludeatoms->_pSysData[i] = 0;
    gpu->dph_accumulator->_pSysData[i] = 0;
    gpu->dph_plus_accumulator->_pSysData[i] = 0;
    gpu->linkgrplist->_pSysData[i] = 0;
  }
  gpu->para->Upload();
  gpu->parb->Upload();
  gpu->park->Upload();
  gpu->barr->Upload();
  gpu->ph_theta->Upload();
  gpu->dph_theta->Upload();
  gpu->vph_theta->Upload();
  gpu->dphold->Upload();
  gpu->vphold->Upload();
  gpu->charge_phmd->Upload();
  gpu->ImageCharge_phmd->Upload();
  gpu->grplist->Upload();
  gpu->linkgrplist->Upload();
  gpu->ImageGrplist->Upload();
  gpu->ImageIndexPHMD->Upload();
  gpu->sp_grp->Upload();
  gpu->sp_par->Upload();
  gpu->qstate1->Upload();
  gpu->qstate2->Upload();
  gpu->ImageQstate1->Upload();
  gpu->ImageQstate2->Upload();
  gpu->vstate1->Upload();
  gpu->vstate2->Upload();
  gpu->ImageVstate1->Upload();
  gpu->ImageVstate2->Upload();
  gpu->qstate1_md->Upload();
  gpu->qstate2_md->Upload();
  gpu->parmod->Upload();
  gpu->resbounds->Upload();
  gpu->baseres->Upload();
  gpu->excludeatoms->Upload();
  gpu->dph_accumulator->Upload();
  gpu->dph_plus_accumulator->Upload();
  gpu->phmd_ene->Upload();
  gpu->dudls->Upload();
  if (gpu->bCPURandoms) {
    cpu_kRandomPHMD(gpu);
  }
  else {
    kRandomPHMD(gpu);
  }
  gpu->sim.pRandomPHMD = gpu->pbRandomPHMD->_pDevData;
  gpu->sim.pdph_theta = gpu->dph_theta->_pDevData;
  gpu->sim.pvph_theta = gpu->vph_theta->_pDevData;
  gpu->sim.pvphold = gpu->vphold->_pDevData;
  gpu->sim.pph_theta = gpu->ph_theta->_pDevData;
  gpu->sim.pdphold = gpu->dphold->_pDevData;
  gpu->sim.presbounds = gpu->resbounds->_pDevData;
  gpu->sim.psp_grp = gpu->sp_grp->_pDevData;
  gpu->sim.pqstate1 = gpu->qstate1->_pDevData;
  gpu->sim.pqstate2 = gpu->qstate2->_pDevData;
  gpu->sim.pImageQstate1 = gpu->ImageQstate1->_pDevData;
  gpu->sim.pImageQstate2 = gpu->ImageQstate2->_pDevData;
  gpu->sim.pImageQstate12 = gpu->ImageQstate1->_pDevData + gpu->sim.paddedNumberOfAtoms;
  gpu->sim.pImageQstate22 = gpu->ImageQstate2->_pDevData + gpu->sim.paddedNumberOfAtoms;
  gpu->sim.pqstate1_md = gpu->qstate1_md->_pDevData;
  gpu->sim.pqstate2_md = gpu->qstate2_md->_pDevData;
  gpu->sim.pcharge_phmd = gpu->charge_phmd->_pDevData;
  gpu->sim.pImageCharge_phmd = gpu->ImageCharge_phmd->_pDevData;
  gpu->sim.pImageCharge_phmd2 = gpu->ImageCharge_phmd->_pDevData + gpu->sim.paddedNumberOfAtoms;
  gpu->sim.pbaseres = gpu->baseres->_pDevData;
  gpu->sim.pparmod = gpu->parmod->_pDevData;
  gpu->sim.pbarr = gpu->barr->_pDevData;
  gpu->sim.ppark = gpu->park->_pDevData;
  gpu->sim.ppara = gpu->para->_pDevData;
  gpu->sim.pparb = gpu->parb->_pDevData;
  gpu->sim.psp_par = gpu->sp_par->_pDevData;
  gpu->sim.pexcludeatoms = gpu->excludeatoms->_pDevData;
  gpu->sim.pgrplist = gpu->grplist->_pDevData;
  gpu->sim.plinkgrplist = gpu->linkgrplist->_pDevData;
  gpu->sim.pImageGrplist = gpu->ImageGrplist->_pDevData;
  gpu->sim.pImageIndexPHMD = gpu->ImageIndexPHMD->_pDevData;
  gpu->sim.pImageIndexPHMD2 = gpu->sim.pImageIndexPHMD + gpu->sim.paddedNumberOfAtoms;
  gpu->sim.pImageGrplist2 = gpu->sim.pImageGrplist + gpu->sim.paddedNumberOfAtoms;
  gpu->sim.pvstate1 = gpu->vstate1->_pDevData;
  gpu->sim.pvstate2 = gpu->vstate2->_pDevData;
  gpu->sim.pImageVstate1 = gpu->ImageVstate1->_pDevData;
  gpu->sim.pImageVstate2 = gpu->ImageVstate2->_pDevData;
  gpu->sim.pImageVstate12 = gpu->ImageVstate1->_pDevData + gpu->sim.paddedNumberOfAtoms;
  gpu->sim.pImageVstate22 = gpu->ImageVstate2->_pDevData + gpu->sim.paddedNumberOfAtoms;
  gpu->sim.pdph_accumulator = gpu->dph_accumulator->_pDevData;
  gpu->sim.pdph_plus_accumulator = gpu->dph_plus_accumulator->_pDevData;
  gpu->sim.pphmd_ene = gpu->phmd_ene->_pDevData;
  gpu->sim.pdudls = gpu->dudls->_pDevData;
  //Borrowed code for exclusions
  gpu->pbNLExclusionStartCount = new GpuBuffer<uint2>(gpu->sim.atoms);
  // Double count exclusions
  unsigned int* pExclusionCount = new unsigned int[gpu->sim.atoms];
  memset(pExclusionCount, 0, gpu->sim.atoms * sizeof(unsigned int));
  int offset = 0;
  for (int i = 0; i < gpu->sim.atoms; i++) {
    for (int j = 0; j < numex[i]; j++) {
      if (natex[offset + j] > 0) {
        pExclusionCount[i] += 1;
        pExclusionCount[natex[offset + j] - 1] += 1;
      }
    }
    offset += numex[i];
  }

  // Count total exclusions
  int exclusions = 0;
  for (int i = 0; i < gpu->sim.atoms; i++) {
    gpu->pbNLExclusionStartCount->_pSysData[i].x = exclusions;
    gpu->pbNLExclusionStartCount->_pSysData[i].y = 0;
    exclusions += ((pExclusionCount[i] + (GRID - 1)) >> GRID_BITS) << GRID_BITS;
  }
  gpu->pbNLExclusionList = new GpuBuffer<unsigned int>(exclusions);
  memset(gpu->pbNLExclusionList->_pSysData, 0xff, exclusions * sizeof(unsigned int));
  offset = 0;
  for (int i = 0; i < gpu->sim.atoms; i++) {
    int pos_i = gpu->pbNLExclusionStartCount->_pSysData[i].x +
                gpu->pbNLExclusionStartCount->_pSysData[i].y;
    for (int j = 0; j < numex[i]; j++) {
      if (natex[offset + j] > 0) {
        gpu->pbNLExclusionList->_pSysData[pos_i++] = natex[offset + j] - 1;
        gpu->pbNLExclusionStartCount->_pSysData[i].y++;
        int pos_j = gpu->pbNLExclusionStartCount->_pSysData[natex[offset + j] - 1].x +
                    gpu->pbNLExclusionStartCount->_pSysData[natex[offset + j] - 1].y++;
        gpu->pbNLExclusionList->_pSysData[pos_j] = i;
      }
    }
    offset += numex[i];
  }
  gpu->pbNLExclusionList->Upload();
  gpu->pbNLExclusionStartCount->Upload();
  gpu->sim.pNLExclusionList       = gpu->pbNLExclusionList->_pDevData;
  gpu->sim.pNLExclusionStartCount = gpu->pbNLExclusionStartCount->_pDevData;
}

//---------------------------------------------------------------------------------------------
// gpu_download_phmd_theta_: Download titration variables.
//
// Arguments:
//   phmd_theta:   the titration variables.
//   dph_theta:    the forces on the titration variables.
//   vph_theta:    the velocities of the titration variables
//---------------------------------------------------------------------------------------------
extern "C" void gpu_download_phmd_theta_(double ph_theta[], double dph_theta[], double vph_theta[], double dudls[])
{
  if(gpu->sim.qsetlam) {
      kUpdateForces_PHMD(gpu);
  }
  gpu->ph_theta->Download();
  gpu->dph_theta->Download();
  gpu->vph_theta->Download();
  gpu->dudls->Download();
  for (int i = 0; i < gpu->sim.ntitr; i++) {
    ph_theta[i] = gpu->ph_theta->_pSysData[i];
    dph_theta[i] = gpu->dph_theta->_pSysData[i];
    vph_theta[i] = gpu->vph_theta->_pSysData[i];
    dudls[i] = gpu->dudls->_pSysData[i];
  }
}

//---------------------------------------------------------------------------------------------
// gpu_upload_charges_: upload the charges for atoms in the simulation.
//
// Arguments:
//   charge:   the atomic partial charges
//---------------------------------------------------------------------------------------------
extern "C" void gpu_upload_charges_(double charge[])
{
  PRINTMETHOD("gpu_upload_charges");
  if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
    if (gpu->bNewNeighborList) {
      gpu->pbImageIndex->Download();
      gpu->bNewNeighborList   = false;
    }
    gpu->pbImage->Download();
    gpu->pbAtomChargeSPLJID->Download();
    unsigned int* pImageAtomLookup = &(gpu->pbImageIndex->_pSysData[gpu->sim.imageStride * 2]);
    double *pCharge = gpu->pbImageCharge->_pSysData;
    if (gpu->sim.pImageCharge != gpu->pbImageCharge->_pDevData) {
      pCharge = gpu->pbImageCharge->_pSysData + gpu->sim.stride;
    }
    for (int i = 0; i < gpu->sim.atoms; i++) {
      int i1 = pImageAtomLookup[i];
      pCharge[i1] = charge[i];
      gpu->pbAtomChargeSP->_pSysData[i1] = charge[i];
      gpu->pbAtomChargeSPLJID->_pSysData[i1].x = charge[i];
    }
    gpu->pbImageCharge->Upload();
  }
  else {
    for (int i = 0; i < gpu->sim.atoms; i++) {
      gpu->pbAtomCharge->_pSysData[i] = charge[i];
      gpu->pbAtomChargeSP->_pSysData[i] = charge[i];
      gpu->pbAtomChargeSPLJID->_pSysData[i].x = charge[i];
    }
    for (int i = gpu->sim.atoms; i < gpu->sim.paddedNumberOfAtoms; i++) {
      gpu->pbAtomCharge->_pSysData[i] = (PMEDouble)0.0;
      gpu->pbAtomChargeSP->_pSysData[i] = (PMEFloat)0.0;
      gpu->pbAtomChargeSPLJID->_pSysData[i].x = (PMEFloat)0.0;
    }
    gpu->pbAtomCharge->Upload();
  }
  gpu->pbAtomChargeSP->Upload();
  gpu->pbAtomChargeSPLJID->Upload();
}

//---------------------------------------------------------------------------------------------
// gpu_upload_charges_pme_cph_: upload atomic partial charges for PME and constant pH
//                              simulations.
//
// Arguments:
//   charge:   the atomic partial charges
//---------------------------------------------------------------------------------------------
extern "C" void gpu_upload_charges_pme_cph_(double charge[])
{
  PRINTMETHOD("gpu_upload_charges_pme_cph");

  FloatShift ljid;
  if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
    if (gpu->bNewNeighborList) {
      gpu->pbImageIndex->Download();
      gpu->bNewNeighborList   = false;
    }
    gpu->pbImage->Download();
    gpu->pbAtomChargeSPLJID->Download();
    gpu->pbImageLJID->Download();
    unsigned int* pImageAtomLookup = &(gpu->pbImageIndex->_pSysData[gpu->sim.imageStride * 2]);
    double *pCharge = gpu->pbImageCharge->_pSysData;
    if (gpu->sim.pImageCharge != gpu->pbImageCharge->_pDevData) {
      pCharge = gpu->pbImageCharge->_pSysData + gpu->sim.stride;
    }
    for (int i = 0; i < gpu->sim.atoms; i++) {
      int i1 = pImageAtomLookup[i];
      pCharge[i1] = charge[i];
      gpu->pbAtomChargeSP->_pSysData[i1] = charge[i];
      gpu->pbAtomChargeSPLJID->_pSysData[i1].x = charge[i];
      ljid.ui = gpu->pbAtomLJID->_pSysData[i];
      gpu->pbAtomChargeSPLJID->_pSysData[i1].y = ljid.f;
    }
    gpu->pbImageCharge->Upload();
  }
  else {
    for (int i = 0; i < gpu->sim.atoms; i++) {
      gpu->pbAtomCharge->_pSysData[i] = charge[i];
      gpu->pbAtomChargeSP->_pSysData[i] = charge[i];
      gpu->pbAtomChargeSPLJID->_pSysData[i].x = charge[i];
      ljid.ui = gpu->pbAtomLJID->_pSysData[i];
      gpu->pbAtomChargeSPLJID->_pSysData[i].y = ljid.f;
    }
    for (int i = gpu->sim.atoms; i < gpu->sim.paddedNumberOfAtoms; i++) {
      gpu->pbAtomCharge->_pSysData[i] = (PMEDouble)0.0;
      gpu->pbAtomChargeSP->_pSysData[i] = (PMEFloat)0.0;
      gpu->pbAtomChargeSPLJID->_pSysData[i].x = (PMEFloat)0.0;
      ljid.ui = 0;
      gpu->pbAtomChargeSPLJID->_pSysData[i].y = ljid.f;
    }
    gpu->pbAtomCharge->Upload();
  }
  gpu->pbAtomChargeSP->Upload();
  gpu->pbAtomChargeSPLJID->Upload();
}

//---------------------------------------------------------------------------------------------
// gpu_upload_charges_gb_cph_: this routine is a method for uploading charges in constant PH
//                             (implicit or explicit solvent) and GB simulations.
//
// Arguments:
//   charge:   the atomic partial charges
//---------------------------------------------------------------------------------------------
extern "C" void gpu_upload_charges_gb_cph_(double charge[])
{
  PRINTMETHOD("gpu_upload_charges_gb_cph");

  FloatShift ljid;
  for (int i = 0; i < gpu->sim.atoms; i++) {
    gpu->pbAtomCharge->_pSysData[i] = charge[i];
    gpu->pbAtomChargeSP->_pSysData[i] = charge[i];
    gpu->pbAtomChargeSPLJID->_pSysData[i].x = charge[i];
    ljid.ui = gpu->pbAtomLJID->_pSysData[i];
    gpu->pbAtomChargeSPLJID->_pSysData[i].y = ljid.f;
  }
  for (int i = gpu->sim.atoms; i < gpu->sim.paddedNumberOfAtoms; i++) {
    gpu->pbAtomCharge->_pSysData[i] = (PMEDouble)0.0;
    gpu->pbAtomChargeSP->_pSysData[i] = (PMEFloat)0.0;
    gpu->pbAtomChargeSPLJID->_pSysData[i].x = (PMEFloat)0.0;
    ljid.ui = 0;
    gpu->pbAtomChargeSPLJID->_pSysData[i].y = ljid.f;
  }
  gpu->pbAtomCharge->Upload();
  gpu->pbAtomChargeSP->Upload();
  gpu->pbAtomChargeSPLJID->Upload();
}

//---------------------------------------------------------------------------------------------
// gpu_download_charges_: this routine is the sole method for download charge information from
//                        the device for all simulations.
//
// Arguments:
//   charge:   the atomic partial charges
//---------------------------------------------------------------------------------------------
extern "C" void gpu_download_charges_(double charge[])
{
  PRINTMETHOD("gpu_download_charges");
  if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
    if (gpu->bNewNeighborList) {
      gpu->pbImageIndex->Download();
      gpu->bNewNeighborList = false;
    }
    gpu->pbImageCharge->Download();
    unsigned int* pImageAtomLookup = &(gpu->pbImageIndex->_pSysData[gpu->sim.imageStride * 2]);
    double *pCharge = gpu->pbImageCharge->_pSysData;
    if (gpu->sim.pImageCharge != gpu->pbImageCharge->_pDevData) {
      pCharge = gpu->pbImageCharge->_pSysData + gpu->sim.stride;
    }
    for (int i = 0; i < gpu->sim.atoms; i++) {
      int i1 = pImageAtomLookup[i];
      charge[i] = pCharge[i1];
    }
  }
  else {
    gpu->pbAtomCharge->Download();
    for (int i = 0; i < gpu->sim.atoms; i++) {
      charge[i] = gpu->pbAtomCharge->_pSysData[i];
    }
  }
}

//---------------------------------------------------------------------------------------------
// gpu_upload_sigeps_: upload sigma and epsilon parameters for simulations.
//
// Arguments:
//   sig:   Lennard-Jones sigma parameters
//   eps:   Lennard-Jones epsilon parameters
//---------------------------------------------------------------------------------------------
extern "C" void gpu_upload_sigeps_(double sig[], double eps[])
{
  PRINTMETHOD("gpu_upload_sigeps");
  for (int i = 0; i < gpu->sim.atoms; i++) {
    gpu->pbAtomSigEps->_pSysData[i].x = sig[i];
    gpu->pbAtomSigEps->_pSysData[i].y = eps[i];
  }
  for (int i = gpu->sim.atoms; i < gpu->sim.paddedNumberOfAtoms; i++) {
    gpu->pbAtomSigEps->_pSysData[i].x = (PMEFloat)0.0;
    gpu->pbAtomSigEps->_pSysData[i].y = (PMEFloat)0.0;
  }
  gpu->pbAtomSigEps->Upload();
}

//---------------------------------------------------------------------------------------------
// gpu_download_sigeps_: the sole function for downloading sigma
//
// Arguments:
//   sig:   Lennard-Jones sigma parameters
//   eps:   Lennard-Jones epsilon parameters
//---------------------------------------------------------------------------------------------
extern "C" void gpu_download_sigeps_(double sig[], double eps[])
{
  PRINTMETHOD("gpu_download_sigeps");
  gpu->pbAtomSigEps->Download();
  for (int i = 0; i < gpu->sim.atoms; i++) {
    sig[i] = gpu->pbAtomSigEps->_pSysData[i].x;
    eps[i] = gpu->pbAtomSigEps->_pSysData[i].y;
  }
}

//---------------------------------------------------------------------------------------------
// gpu_upload_fs_: upload screening terms for GB calculations to the device
//
// Arguments:
//   fs:    atomic screening parameters
//---------------------------------------------------------------------------------------------
extern "C" void gpu_upload_fs_(double fs[])
{
  PRINTMETHOD("gpu_upload_fs");
  for (int i = 0; i < gpu->sim.atoms; i++) {
    gpu->pbAtomS->_pSysData[i] = fs[i];
  }
  for (int i = gpu->sim.atoms; i < gpu->sim.paddedNumberOfAtoms; i++) {
    gpu->pbAtomS->_pSysData[i] = (PMEFloat)0.0;
  }
  gpu->pbAtomS->Upload();
}

//---------------------------------------------------------------------------------------------
// gpu_download_fs_: download screening terms for GB calculations from the device
//
// Arguments:
//   fs:    atomic screening parameters
//---------------------------------------------------------------------------------------------
extern "C" void gpu_download_fs_(double fs[])
{
  PRINTMETHOD("gpu_download_fs");
  gpu->pbAtomS->Download();
  for (int i = 0; i < gpu->sim.atoms; i++) {
    fs[i]   = gpu->pbAtomS->_pSysData[i];
  }
}

//---------------------------------------------------------------------------------------------
// gpu_upload_rborn_: upload atomic Born radii to the device.
//
// Arguments:
//   rborn:   atomic Born radii
//---------------------------------------------------------------------------------------------
extern "C" void gpu_upload_rborn_(double rborn[])
{
  PRINTMETHOD("gpu_upload_rborn");
  for (int i = 0; i < gpu->sim.atoms; i++) {
    gpu->pbAtomRBorn->_pSysData[i] = rborn[i];
  }
  for (int i = gpu->sim.atoms; i < gpu->sim.paddedNumberOfAtoms; i++) {
    gpu->pbAtomRBorn->_pSysData[i] = (PMEFloat)0.0001;
    gpu->pbReff->_pSysData[i] = (PMEFloat)0.0001;
    gpu->pbReffSP->_pSysData[i] = (PMEFloat)0.0001;
    gpu->pbTemp7->_pSysData[i] = (PMEFloat)0.0;
  }
  gpu->pbAtomRBorn->Upload();
  gpu->pbReff->Upload();
  gpu->pbReffSP->Upload();
  gpu->pbTemp7->Upload();
}

//---------------------------------------------------------------------------------------------
// gpu_download_rborn_: download atomic Born radii from the device in GB simulations
//
// Arguments:
//   rborn:   atomic Born radii
//---------------------------------------------------------------------------------------------
extern "C" void gpu_download_rborn_(double rborn[])
{
  PRINTMETHOD("gpu_download_rborn");
  gpu->pbAtomRBorn->Download();
  for (int i = 0; i < gpu->sim.atoms; i++) {
    rborn[i] = gpu->pbAtomRBorn->_pSysData[i];
  }
}

//---------------------------------------------------------------------------------------------
// gpu_upload_masses_: upload the atomic masses to the device
//
// Arguments:
//   mass:    atomic masses in Daltons
//---------------------------------------------------------------------------------------------
extern "C" void gpu_upload_masses_(double mass[])
{
  PRINTMETHOD("gpu_upload_masses");
  for (int i = 0; i < gpu->sim.atoms; i++) {
    gpu->pbAtomMass->_pSysData[i] = mass[i];
    if (mass[i] > 0.0) {
      gpu->pbAtomMass->_pSysData[i + gpu->sim.stride] = (PMEDouble)(1.0 / mass[i]);
    }
    else {
      gpu->pbAtomMass->_pSysData[i + gpu->sim.stride] = (PMEDouble)0.0;
    }
  }
  for (int i = gpu->sim.atoms; i < gpu->sim.paddedNumberOfAtoms; i++) {
    gpu->pbAtomMass->_pSysData[i] = (PMEDouble)0.0;
    gpu->pbAtomMass->_pSysData[i + gpu->sim.stride] = (PMEDouble)10000000.0;
  }
  gpu->pbAtomMass->Upload();
}

//---------------------------------------------------------------------------------------------
// gpu_download_masses_: download the atomic masses to the device
//
// Arguments:
//   mass:    atomic masses in Daltons
//---------------------------------------------------------------------------------------------
extern "C" void gpu_download_masses_(double mass[])
{
  PRINTMETHOD("gpu_download_masses");
  gpu->pbAtomMass->Download();
  for (int i = 0; i < gpu->sim.atoms; i++) {
    mass[i] = gpu->pbAtomMass->_pSysData[i];
  }
}

//---------------------------------------------------------------------------------------------
// gpu_upload_reff_: upload the effective Born radii for GB simulations
//
// Arguments:
//   reff:    the effective Born radii
//---------------------------------------------------------------------------------------------
extern "C" void gpu_upload_reff_(double reff[])
{
  PRINTMETHOD("gpu_upload_reff");
  for (int i = 0; i < gpu->sim.atoms; i++) {
    gpu->pbReff->_pSysData[i] = reff[i];
  }
  for (int i = gpu->sim.atoms; i < gpu->sim.paddedNumberOfAtoms; i++) {
    gpu->pbReff->_pSysData[i] = (PMEFloat)0.0001;
  }
  gpu->pbReff->Upload();
}

//---------------------------------------------------------------------------------------------
// gpu_download_reff_: download the effective Born radii for GB simulations
//
// Arguments:
//   reff:    the effective Born radii
//---------------------------------------------------------------------------------------------
extern "C" void gpu_download_reff_(double reff[])
{
  PRINTMETHOD("gpu_download_reff");
  gpu->pbReff->Download();
  for (int i = 0; i < gpu->sim.atoms; i++) {
    reff[i] = gpu->pbReff->_pSysData[i];
  }
}

//---------------------------------------------------------------------------------------------
// gpu_upload_frc_: upload the forces on all atoms to the device.  While it is anathema to
//                  transfer so much data between the host and the device, this is needed in
//                  Monte Carlo barostating (even when a change in box volume succeeds, the
//                  old forces must still be reloaded) and in REMD (when a move fails, the old
//                  forces again must be reloaded).
//
// Arguments:
//   atm_frc:    the forces on all atoms
//---------------------------------------------------------------------------------------------
extern "C" void gpu_upload_frc_(double atm_frc[][3])
{
  PRINTMETHOD("gpu_upload_frc");
  if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
    if (gpu->bNewNeighborList) {
      gpu->pbImageIndex->Download();
      gpu->bNewNeighborList = false;
    }
    unsigned int* pImageAtomLookup = &(gpu->pbImageIndex->_pSysData[gpu->sim.imageStride * 2]);
    PMEAccumulator *pForce = gpu->pbForceAccumulator->_pSysData;
    for (int i = 0; i < gpu->sim.atoms; i++) {
      int i1 = pImageAtomLookup[i];
      pForce[i1] = (PMEAccumulator)(FORCESCALE * atm_frc[i][0]);
      pForce[i1 + gpu->sim.stride] = (PMEAccumulator)(FORCESCALE * atm_frc[i][1]);
      pForce[i1 + gpu->sim.stride2] = (PMEAccumulator)(FORCESCALE * atm_frc[i][2]);
    }
    gpu->pbForceAccumulator->Upload();
  }
  else {
    for (int i = 0; i < gpu->sim.atoms; i++) {
      gpu->pbForceAccumulator->_pSysData[i] = (PMEAccumulator)(FORCESCALE * atm_frc[i][0]);
      gpu->pbForceAccumulator->_pSysData[i + gpu->sim.stride] =
        (PMEAccumulator)(FORCESCALE * atm_frc[i][1]);
      gpu->pbForceAccumulator->_pSysData[i + gpu->sim.stride2] =
        (PMEAccumulator)(FORCESCALE * atm_frc[i][2]);
    }
    gpu->pbForceAccumulator->Upload();
  }
}

//---------------------------------------------------------------------------------------------
// gpu_upload_frc_add_: this routine is needed by methods that add biasing forces to the atoms
//                      which are not (yet) calculated on the GPU.  Specifically, the method
//                      is adaptively biased MD, in both GB and PME situations.
//
// Arguments:
//   atm_frc:    forces on all atoms
//---------------------------------------------------------------------------------------------
extern "C" void gpu_upload_frc_add_(double atm_frc[][3])
{
  PRINTMETHOD("gpu_upload_frc_add");
  gpu->pbForceAccumulator->Download();
  if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
    if (gpu->bNewNeighborList) {
      gpu->pbImageIndex->Download();
      gpu->bNewNeighborList = false;
    }
    unsigned int* pImageAtomLookup = &(gpu->pbImageIndex->_pSysData[gpu->sim.imageStride * 2]);
    PMEAccumulator *pForce  = gpu->pbForceAccumulator->_pSysData;
    for (int i = 0; i < gpu->sim.atoms; i++) {
      int i1 = pImageAtomLookup[i];
      pForce[i1] += (PMEAccumulator)(FORCESCALE * atm_frc[i][0]);
      pForce[i1 + gpu->sim.stride] += (PMEAccumulator)(FORCESCALE * atm_frc[i][1]);
      pForce[i1 + gpu->sim.stride2] += (PMEAccumulator)(FORCESCALE * atm_frc[i][2]);
    }
    gpu->pbForceAccumulator->Upload();
  }
  else {
    for (int i = 0; i < gpu->sim.atoms; i++) {
      gpu->pbForceAccumulator->_pSysData[i] += (PMEAccumulator)(FORCESCALE * atm_frc[i][0]);
      gpu->pbForceAccumulator->_pSysData[i + gpu->sim.stride] +=
        (PMEAccumulator)(FORCESCALE * atm_frc[i][1]);
      gpu->pbForceAccumulator->_pSysData[i + gpu->sim.stride2] +=
        (PMEAccumulator)(FORCESCALE * atm_frc[i][2]);
    }
    gpu->pbForceAccumulator->Upload();
  }
}

//---------------------------------------------------------------------------------------------
// gpu_gbsa_frc_add_: update forces on the GPU device based on GBSA contributions calculated
//                    on the host.
//
// Arguments:
//   atm_frc:    forces on all atoms
//---------------------------------------------------------------------------------------------
extern "C" void gpu_gbsa_frc_add_(double atm_frc[][3])
{
  PRINTMETHOD("gpu_gbsa_frc_add");
  gpu->pbForceAccumulator->Download();
  if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
    if (gpu->bNewNeighborList) {
      gpu->pbImageIndex->Download();
      gpu->bNewNeighborList = false;
    }
    unsigned int* pImageAtomLookup = &(gpu->pbImageIndex->_pSysData[gpu->sim.stride2]);
    PMEAccumulator *pForce  = gpu->pbForceAccumulator->_pSysData;
    for (int i = 0; i < gpu->sim.atoms; i++) {
      int i1 = pImageAtomLookup[i];
      pForce[i1] += (PMEAccumulator)(FORCESCALE * atm_frc[i][0]);
      pForce[i1 + gpu->sim.stride] += (PMEAccumulator)(FORCESCALE * atm_frc[i][1]);
      pForce[i1 + gpu->sim.stride2] += (PMEAccumulator)(FORCESCALE * atm_frc[i][2]);
    }
  }
  else {
    for (int i = 0; i < gpu->sim.atoms; i++) {
      gpu->pbForceAccumulator->_pSysData[i] += (PMEAccumulator)(FORCESCALE * atm_frc[i][0]);
      gpu->pbForceAccumulator->_pSysData[i + gpu->sim.stride] +=
        (PMEAccumulator)(FORCESCALE * atm_frc[i][1]);
      gpu->pbForceAccumulator->_pSysData[i + gpu->sim.stride2] +=
        (PMEAccumulator)(FORCESCALE * atm_frc[i][2]);
    }
  }
  gpu->pbForceAccumulator->Upload();
}

//---------------------------------------------------------------------------------------------
// gpu_download_frc_: download forces from the GPU device.  Generally, this is no longer used
//                    when updating forces on the GPU with contributions calculated on the
//                    host--those contributions are simply uploaded.  However, it is called
//                    in REMD and MC barostating, when failed moves mean that system must
//                    revert to its old state and use the old forces.
//
// Arguments:
//   atm_frc:   the forces on all atoms
//---------------------------------------------------------------------------------------------
extern "C" void gpu_download_frc_(double atm_frc[][3])
{
  PRINTMETHOD("gpu_download_frc");
  gpu->pbForceAccumulator->Download();
  PMEAccumulator *pForce  = gpu->pbForceAccumulator->_pSysData;
  if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
    if (gpu->bNewNeighborList) {
      gpu->pbImageIndex->Download();
      gpu->bNewNeighborList = false;
    }
    unsigned int* pImageAtomLookup = &(gpu->pbImageIndex->_pSysData[gpu->sim.imageStride * 2]);
#ifndef MPI
    if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
      PMEAccumulator *pNBForce  = gpu->pbForceAccumulator->_pSysData + gpu->sim.stride3;
      for (int i = 0; i < gpu->sim.atoms; i++) {
        int i1 = pImageAtomLookup[i];
        atm_frc[i][0] = (double)(pForce[i1] + pNBForce[i1]) * (double)ONEOVERFORCESCALE;
        atm_frc[i][1] = (double)(pForce[i1 + gpu->sim.stride] +
                                 pNBForce[i1 + gpu->sim.stride]) * (double)ONEOVERFORCESCALE;
        atm_frc[i][2] = (double)(pForce[i1 + gpu->sim.stride2] +
                                 pNBForce[i1 + gpu->sim.stride2]) * (double)ONEOVERFORCESCALE;
#ifdef use_DPFP
        if (gpu->imin != 0)
        {
          atm_frc[i][0]    += (double)pNBForce[i1 + gpu->sim.stride3];
          atm_frc[i][1]    += (double)pNBForce[i1 + gpu->sim.stride4];
          atm_frc[i][2]    += (double)pNBForce[i1 + gpu->sim.stride5];
        }
#endif
      }
    }
    else {
#endif
    for (int i = 0; i < gpu->sim.atoms; i++) {
      int i1              = pImageAtomLookup[i];
      atm_frc[i][0]       = (double)pForce[i1] * (double)ONEOVERFORCESCALE;
      atm_frc[i][1]       = (double)pForce[i1 + gpu->sim.stride] * (double)ONEOVERFORCESCALE;
      atm_frc[i][2]       = (double)pForce[i1 + gpu->sim.stride2] * (double)ONEOVERFORCESCALE;
#ifdef use_DPFP
      if (gpu->imin != 0)
      {
        atm_frc[i][0]    += (double)pForce[i1 + gpu->sim.stride3];
        atm_frc[i][1]    += (double)pForce[i1 + gpu->sim.stride4];
        atm_frc[i][2]    += (double)pForce[i1 + gpu->sim.stride5];
      }
#endif
    }
#ifndef MPI
    }
#endif
  }
  else {
    for (int i = 0; i < gpu->sim.atoms; i++) {
      atm_frc[i][0] = (double)pForce[i] * (double)ONEOVERFORCESCALE;
      atm_frc[i][1] = (double)pForce[i + gpu->sim.stride] * (double)ONEOVERFORCESCALE;
      atm_frc[i][2] = (double)pForce[i + gpu->sim.stride2] * (double)ONEOVERFORCESCALE;
#ifdef use_DPFP
      if (gpu->imin != 0)
      {
        atm_frc[i][0]    += (double)pForce[i + gpu->sim.stride3];
        atm_frc[i][1]    += (double)pForce[i + gpu->sim.stride4];
        atm_frc[i][2]    += (double)pForce[i + gpu->sim.stride5];
      }
#endif
    }
  }
}

//---------------------------------------------------------------------------------------------
// gpu_upload_vel_: upload velocities from the host to the device.  Invoked at the times one
//                  might expect, to set up the simulation after a checkpoint file, to write a
//                  checkpoint file, to perform a REMD swap.  Also invoked to zero the center
//                  of mass velocity--that's done on the CPU.
//
// Arguments:
//   atm_vel:  the velocities of all atoms
//---------------------------------------------------------------------------------------------
extern "C" void gpu_upload_vel_(double atm_vel[][3])
{
  PRINTMETHOD("gpu_upload_vel");
  if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
    if (gpu->bNewNeighborList) {
      gpu->pbImageIndex->Download();
      gpu->bNewNeighborList = false;
    }
    gpu->pbImageVel->Download();
    unsigned int* pImageAtomLookup = &(gpu->pbImageIndex->_pSysData[gpu->sim.imageStride * 2]);
    double *pVel = gpu->pbImageVel->_pSysData;
    if (gpu->sim.pImageVelX != gpu->pbImageVel->_pDevData) {
      pVel = gpu->pbImageVel->_pSysData + gpu->sim.stride3;
    }
    for (int i = 0; i < gpu->sim.atoms; i++) {
      int i1 = pImageAtomLookup[i];
      pVel[i1] = atm_vel[i][0];
      pVel[i1 + gpu->sim.stride] = atm_vel[i][1];
      pVel[i1 + gpu->sim.stride2] = atm_vel[i][2];
    }
    gpu->pbImageVel->Upload();
  }
  else {
    for (int i = 0; i < gpu->sim.atoms; i++) {
      gpu->pbVel->_pSysData[i] = atm_vel[i][0];
      gpu->pbVel->_pSysData[i + gpu->sim.stride] = atm_vel[i][1];
      gpu->pbVel->_pSysData[i + gpu->sim.stride2] = atm_vel[i][2];
    }
    gpu->pbVel->Upload();
  }
}

//---------------------------------------------------------------------------------------------
// gpu_download_vel_: download the velocities of all atoms from the device to the host.
//
// Arguments:
//   atm_vel:  the velocities of all atoms
//---------------------------------------------------------------------------------------------
extern "C" void gpu_download_vel_(double atm_vel[][3])
{
  PRINTMETHOD("gpu_download_vel");
  if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
    if (gpu->bNewNeighborList) {
      gpu->pbImageIndex->Download();
      gpu->bNewNeighborList = false;
    }
    gpu->pbImageVel->Download();
    unsigned int* pImageAtomLookup = &(gpu->pbImageIndex->_pSysData[gpu->sim.imageStride * 2]);
    double *pVel = gpu->pbImageVel->_pSysData;
    if (gpu->sim.pImageVelX != gpu->pbImageVel->_pDevData) {
      pVel = gpu->pbImageVel->_pSysData + gpu->sim.stride3;
    }
    for (int i = 0; i < gpu->sim.atoms; i++) {
      int i1 = pImageAtomLookup[i];
      atm_vel[i][0] = pVel[i1];
      atm_vel[i][1] = pVel[i1 + gpu->sim.stride];
      atm_vel[i][2] = pVel[i1 + gpu->sim.stride2];
    }
  }
  else {
    gpu->pbVel->Download();
    for (int i = 0; i < gpu->sim.atoms; i++) {
      atm_vel[i][0] = gpu->pbVel->_pSysData[i];
      atm_vel[i][1] = gpu->pbVel->_pSysData[i + gpu->sim.stride];
      atm_vel[i][2] = gpu->pbVel->_pSysData[i + gpu->sim.stride2];
    }
  }
}

//---------------------------------------------------------------------------------------------
// gpu_upload_last_vel_: upload the previous step's velocities.  This is called during system
//                       initialization.
//
// Arguments:
//   atm_last_vel:  the atomic velocities from a previous step
//---------------------------------------------------------------------------------------------
extern "C" void gpu_upload_last_vel_(double atm_last_vel[][3])
{
  PRINTMETHOD("gpu_upload_last_vel");
  if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
    if (gpu->bNewNeighborList) {
      gpu->pbImageIndex->Download();
      gpu->bNewNeighborList = false;
    }
    gpu->pbImageLVel->Download();
    unsigned int* pImageAtomLookup = &(gpu->pbImageIndex->_pSysData[gpu->sim.imageStride * 2]);
    double *pLVel = gpu->pbImageLVel->_pSysData;
    if (gpu->sim.pImageLVelX != gpu->pbImageLVel->_pDevData) {
      pLVel = gpu->pbImageLVel->_pSysData + gpu->sim.stride3;
    }
    for (int i = 0; i < gpu->sim.atoms; i++) {
      int i1 = pImageAtomLookup[i];
      pLVel[i1] = atm_last_vel[i][0];
      pLVel[i1 + gpu->sim.stride] = atm_last_vel[i][1];
      pLVel[i1 + gpu->sim.stride2] = atm_last_vel[i][2];
    }
    gpu->pbImageLVel->Upload();
  }
  else {
    for (int i = 0; i < gpu->sim.atoms; i++) {
      gpu->pbLVel->_pSysData[i] = atm_last_vel[i][0];
      gpu->pbLVel->_pSysData[i + gpu->sim.stride] = atm_last_vel[i][1];
      gpu->pbLVel->_pSysData[i + gpu->sim.stride2]= atm_last_vel[i][2];
    }
    gpu->pbLVel->Upload();
  }
}

//---------------------------------------------------------------------------------------------
// gpu_download_last_vel_: download the velocities from a previous step.  Not sure where this
//                         is called.
//
// Arguments:
//   atm_last_vel:  the atomic velocities from a previous step
//---------------------------------------------------------------------------------------------
extern "C" void gpu_download_last_vel_(double atm_last_vel[][3])
{
  PRINTMETHOD("gpu_download_last_vel");
  if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
    if (gpu->bNewNeighborList) {
      gpu->pbImageIndex->Download();
      gpu->bNewNeighborList = false;
    }
    gpu->pbImageLVel->Download();
    unsigned int* pImageAtomLookup = &(gpu->pbImageIndex->_pSysData[gpu->sim.imageStride * 2]);
    double *pLVel = gpu->pbImageLVel->_pSysData;
    if (gpu->sim.pImageLVelX != gpu->pbImageLVel->_pDevData) {
      pLVel = gpu->pbImageLVel->_pSysData + gpu->sim.stride3;
    }
    for (int i = 0; i < gpu->sim.atoms; i++) {
      int i1 = pImageAtomLookup[i];
      atm_last_vel[i][0] = pLVel[i1];
      atm_last_vel[i][1] = pLVel[i1 + gpu->sim.stride];
      atm_last_vel[i][2] = pLVel[i1 + gpu->sim.stride2];
    }
  }
  else {
    gpu->pbLVel->Download();
    for (int i = 0; i < gpu->sim.atoms; i++) {
      atm_last_vel[i][0] = gpu->pbLVel->_pSysData[i];
      atm_last_vel[i][1] = gpu->pbLVel->_pSysData[i + gpu->sim.stride];
      atm_last_vel[i][2] = gpu->pbLVel->_pSysData[i + gpu->sim.stride2];
    }
  }
}

//---------------------------------------------------------------------------------------------
// gpu_clear_vel_: zero out the velocities.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_clear_vel_()
{
  PRINTMETHOD("gpu_clear_vel");
  kClearVelocities(gpu);
}

//---------------------------------------------------------------------------------------------
// gpu_bonds_setup_: setup for bonded interactions on the GPU.
//
// Arguments:
//   cit_nbona:   the number of bonds between heavy atoms
//   cit_nbonh:   the number of bonds to hydrogen
//   cit_a_bond:  array of bond parameters for bonds between heavy atoms
//   cit_h_bond:  array of bond parameters for bonds to hydrogen
//   gbl_req:     equilbrium lengths for all bond types
//   gbl_rk:      stiffness constants for all bond types
//---------------------------------------------------------------------------------------------
extern "C" void gpu_bonds_setup_(int* cit_nbona, bond_rec cit_a_bond[], int* cit_nbonh,
                                 bond_rec cit_h_bond[], double gbl_req[], double gbl_rk[])
{
  PRINTMETHOD("gpu_bonds_setup");

  // Count non-zero bonds
  int bonds = 0;
  if (gpu->ntf < 3) {
    for (int i = 0; i < *cit_nbona; i++) {
      if (gbl_rk[cit_a_bond[i].parm_idx - 1]) {
        bonds++;
      }
    }
  }
  if (gpu->ntf < 2) {
    for (int i = 0; i < *cit_nbonh; i++) {
      if (gbl_rk[cit_h_bond[i].parm_idx - 1]) {
        bonds++;
      }
    }
  }

  // Allocate/reallocate GPU bond buffers
  delete gpu->pbBond;
  delete gpu->pbBondID;

#ifdef GVERBOSE
  printf("%d bonds, %d active\n", *cit_nbona + *cit_nbonh, bonds);
#endif
  gpu->pbBond = new GpuBuffer<PMEDouble2>(bonds);
  gpu->pbBondID = new GpuBuffer<int2>(bonds);
  bonds = 0;
  if (gpu->ntf < 3) {
    for (int i = 0; i < *cit_nbona; i++) {
      if (gbl_rk[cit_a_bond[i].parm_idx - 1]) {
        gpu->pbBond->_pSysData[bonds].x = gbl_rk[cit_a_bond[i].parm_idx - 1];
        gpu->pbBond->_pSysData[bonds].y = gbl_req[cit_a_bond[i].parm_idx - 1];
        gpu->pbBondID->_pSysData[bonds].x = abs(cit_a_bond[i].atm_i) - 1;
        gpu->pbBondID->_pSysData[bonds].y = abs(cit_a_bond[i].atm_j) - 1;
        bonds++;
      }
    }
  }
  if (gpu->ntf < 2) {
    for (int i = 0; i < *cit_nbonh; i++) {
      if (gbl_rk[cit_h_bond[i].parm_idx - 1]) {
        gpu->pbBond->_pSysData[bonds].x = gbl_rk[cit_h_bond[i].parm_idx - 1];
        gpu->pbBond->_pSysData[bonds].y = gbl_req[cit_h_bond[i].parm_idx - 1];
        gpu->pbBondID->_pSysData[bonds].x = abs(cit_h_bond[i].atm_i) - 1;
        gpu->pbBondID->_pSysData[bonds].y = abs(cit_h_bond[i].atm_j) - 1;
        bonds++;
      }
    }
  }
  gpu->pbBond->Upload();
  gpu->pbBondID->Upload();

  // Set constants
  gpu->sim.bonds = bonds;
  gpu->sim.pBond = gpu->pbBond->_pDevData;
  gpu->sim.pBondID = gpu->pbBondID->_pDevData;
  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gpu_angles_setup_: setup for angle interactions on the GPU.
//
// Arguments:
//   angle_cnt:   the number of angles in the system
//   ntheth:      the number of angle parameter types involving hydrogens
//   cit_angle:   parameter arrays for all angles
//   gbl_teq:     angle equilibrium values
//   gbl_tk:      angle stiffness constants
//---------------------------------------------------------------------------------------------
extern "C" void gpu_angles_setup_(int* angle_cnt, int* ntheth, angle_rec cit_angle[],
                                  double gbl_teq[], double gbl_tk[])
{
  PRINTMETHOD("gpu_angles_setup");

  // Allocate/reallocate GPU bond angle buffers
  delete gpu->pbBondAngle;
  delete gpu->pbBondAngleID1;
  delete gpu->pbBondAngleID2;
  int angles = 0;
  if (gpu->ntf < 4) {
    angles = *angle_cnt;
  }
  else if (gpu->ntf < 5) {
    angles = *angle_cnt - *ntheth;
  }

#ifdef GVERBOSE
  printf("%d bond angles, %d active\n", *angle_cnt, angles);
#endif
  gpu->pbBondAngle = new GpuBuffer<PMEDouble2>(angles);
  gpu->pbBondAngleID1 = new GpuBuffer<int2>(angles);
  gpu->pbBondAngleID2 = new GpuBuffer<int>(angles);

  // Copy bond angles
  angles = 0;
  if (gpu->ntf < 5) {
    for (int i = 0; i < *angle_cnt; i++) {
      if ((gpu->ntf < 4) || (i >= *ntheth)) {
        gpu->pbBondAngle->_pSysData[angles].x = gbl_tk[cit_angle[i].parm_idx -1];
        gpu->pbBondAngle->_pSysData[angles].y = gbl_teq[cit_angle[i].parm_idx - 1];
        gpu->pbBondAngleID1->_pSysData[angles].x = abs(cit_angle[i].atm_i) - 1;
        gpu->pbBondAngleID1->_pSysData[angles].y = abs(cit_angle[i].atm_j) - 1;
        gpu->pbBondAngleID2->_pSysData[angles] = abs(cit_angle[i].atm_k) - 1;
        angles++;
      }
    }
  }
  gpu->pbBondAngle->Upload();
  gpu->pbBondAngleID1->Upload();
  gpu->pbBondAngleID2->Upload();

  // Set constants
  gpu->sim.bondAngles = angles;
  gpu->sim.pBondAngle = gpu->pbBondAngle->_pDevData;
  gpu->sim.pBondAngleID1 = gpu->pbBondAngleID1->_pDevData;
  gpu->sim.pBondAngleID2 = gpu->pbBondAngleID2->_pDevData;
  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gmul: array of even numbers interspersed with zeros for dihedral computations
//---------------------------------------------------------------------------------------------
static PMEDouble gmul[] = {0.0, 2.0, 0.0, 4.0, 0.0, 6.0, 0.0, 8.0, 0.0, 10.0};

//---------------------------------------------------------------------------------------------
// gpu_dihedrals_setup_: setup for dihedral computations on the GPU device
//
// Arguments:
//   dihed_cnt:    the number of dihedrals in the system
//   nphih:        the number of dihedrals involving hydrogen (with NTF = 6 these are omitted)
//   cit_dihed:    arrays of dihedral parameters and atom indices
//   gbl_ipn:      integer representations of the periodicity for each dihedral type
//   gbl_pn:       floating-point representation of the periodicity for each dihedral type
//   gbl_pk:       amplitudes for each dihedral (all dihedrals are
//   gbl_gamc:
//   gbl_gams:
//---------------------------------------------------------------------------------------------
extern "C" void gpu_dihedrals_setup_(int* dihed_cnt, int* nphih, dihed_rec cit_dihed[],
                                     int gbl_ipn[], double gbl_pn[], double gbl_pk[],
                                     double gbl_gamc[], double gbl_gams[])
{
  PRINTMETHOD("gpu_dihedrals_setup");

  // Allocate/reallocate GPU dihedral buffers
  delete gpu->pbDihedral1;
  delete gpu->pbDihedral2;
  delete gpu->pbDihedral3;
  delete gpu->pbDihedralID1;
  int dihedrals = 0;
  if (gpu->ntf < 6) {
    dihedrals = *dihed_cnt;
  }
  else if (gpu->ntf < 7) {
    dihedrals = *dihed_cnt - *nphih;
  }

#ifdef GVERBOSE
  printf("%d dihedrals, %d active\n", *dihed_cnt, dihedrals);
#endif
  gpu->pbDihedral1 = new GpuBuffer<PMEFloat2>(dihedrals);
  gpu->pbDihedral2 = new GpuBuffer<PMEFloat2>(dihedrals);
  gpu->pbDihedral3 = new GpuBuffer<PMEFloat>(dihedrals);
  gpu->pbDihedralID1 = new GpuBuffer<int4>(dihedrals);

  // Copy dihedrals
  dihedrals = 0;
  if (gpu->ntf < 7) {
    for (int i = 0; i < *dihed_cnt; i++) {
      if ((gpu->ntf < 6) || (i >= *nphih)) {
        gpu->pbDihedral1->_pSysData[dihedrals].x = gmul[gbl_ipn[cit_dihed[i].parm_idx -1]];
        gpu->pbDihedral1->_pSysData[dihedrals].y = gbl_pn[cit_dihed[i].parm_idx - 1];
        gpu->pbDihedral2->_pSysData[dihedrals].x = gbl_pk[cit_dihed[i].parm_idx - 1];
        gpu->pbDihedral2->_pSysData[dihedrals].y = gbl_gamc[cit_dihed[i].parm_idx - 1];
        gpu->pbDihedral3->_pSysData[dihedrals] = gbl_gams[cit_dihed[i].parm_idx - 1];
        gpu->pbDihedralID1->_pSysData[dihedrals].x = abs(cit_dihed[i].atm_i) - 1;
        gpu->pbDihedralID1->_pSysData[dihedrals].y = abs(cit_dihed[i].atm_j) - 1;
        gpu->pbDihedralID1->_pSysData[dihedrals].z = abs(cit_dihed[i].atm_k) - 1;
        gpu->pbDihedralID1->_pSysData[dihedrals].w = abs(cit_dihed[i].atm_l) - 1;
        dihedrals++;
      }
    }
  }

  gpu->pbDihedral1->Upload();
  gpu->pbDihedral2->Upload();
  gpu->pbDihedral3->Upload();
  gpu->pbDihedralID1->Upload();

  // Set constants
  gpu->sim.dihedrals = dihedrals;
  gpu->sim.pDihedral1 = gpu->pbDihedral1->_pDevData;
  gpu->sim.pDihedral2 = gpu->pbDihedral2->_pDevData;
  gpu->sim.pDihedral3 = gpu->pbDihedral3->_pDevData;
  gpu->sim.pDihedralID1 = gpu->pbDihedralID1->_pDevData;
  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gpu_nmr_get_jar_data_: for calculations with Jarzynski's equality.  This takes data FROM the
//                        GPU and downloads it to the host.
//
// Arguments:
//   r2:
//   rint:
//   fcurr:
//   work:
//---------------------------------------------------------------------------------------------
extern "C" void gpu_nmr_get_jar_data_(double* r2, double* rint, double* fcurr, double* work)
{
  PRINTMETHOD("gpu_nmr_get_jar_data");
  gpu->pbNMRJarData->Download();
  *r2 = gpu->pbNMRJarData->_pSysData[0];
  *rint = gpu->pbNMRJarData->_pSysData[1];
  *fcurr = gpu->pbNMRJarData->_pSysData[2];
  *work = gpu->pbNMRJarData->_pSysData[3];
}

//---------------------------------------------------------------------------------------------
// gpu_nmr_set_nstep_: this function is called by nmr_calls.F90 at various times to set the
//                     time step as the NMR calculation is to see it.  The NMR restraints can
//                     be programmed to evolve over time, and this function is critical to
//                     orchestrating that evolution.
//
// Arguments:
//   nstep:    the step count to set
//---------------------------------------------------------------------------------------------
extern "C" void gpu_nmr_set_nstep_(int* nstep)
{
  PRINTMETHOD("gpu_nmr_set_nstep");
  gpu->NMRnstep = *nstep;
}

//---------------------------------------------------------------------------------------------
// gpu_nmr_setup_: setup for GPU-based NMR calculations of all kinds
// Added COM angle and torsions
// Arguments:
//   nmrnum:  number of restaints
//   resttype: Determines which type of restraints e.g distance, angles, torsionas
//---------------------------------------------------------------------------------------------
extern "C" void gpu_nmr_setup_(int* nmrnum, int nmrat[][16], int nmrst[][3], double r1nmr[][2],
                               double r2nmr[][2], double r3nmr[][2], double r4nmr[][2],
                               double rk2nmr[][2], double rk3nmr[][2], int* jar,
                               double* drjar, int nmrcom[][2], int* maxgrp, int igravt[],
                               int ifxyz[], double idxyz[], int resttype[])
{
  PRINTMETHOD("gpu_nmr_setup");


  // Count restraints of each type
  int distances = 0;
  int COMdistances = 0;
  int r6avdistances = 0;
  int angles = 0;
  int torsions = 0;
  int COMangles = 0;
  int COMtorsions = 0;

  // let us figure how much of a type of restraints
  for (int i = 0; i < *nmrnum; i++) {
    // distances, resttype = 1
    if (resttype[i] == 1) {
      if (nmrat[i][0] >= 0 && nmrat[i][1] >= 0) {
        distances++;
      }
      else if (igravt[i] == 1) {
        r6avdistances++;
      }
      else {
        COMdistances++;
      }
    }
    // angles, resttype = 2
    else if (resttype[i] == 2) {
      if (nmrat[i][0] >= 0 && nmrat[i][1] >= 0 && nmrat[i][2] >= 0) {
        angles++;
      }
      else {
        COMangles++;
      }
    }
    // torsions, resttype = 3
    else if (resttype[i] == 3) {
      if (nmrat[i][0] >= 0 && nmrat[i][1] >= 0 && nmrat[i][2] >= 0 && nmrat[i][3] >= 0) {
        torsions++;
      }
      else {
        COMtorsions++;
      }
    }
  } // end of nmr


  bool r6restraints = false;
  if(r6avdistances > 0) {
       r6restraints = true;
  }
  gpu->sim.NMRR6av = r6restraints;
  gpu->sim.NMRDistances = distances;
  gpu->sim.NMRAngles = angles;
  gpu->sim.NMRTorsions = torsions;
  gpu->sim.NMRCOMDistances = COMdistances;
  gpu->sim.NMRr6avDistances = r6avdistances;
  gpu->sim.NMRCOMAngles = COMangles;
  gpu->sim.NMRCOMTorsions = COMtorsions;

  gpu->sim.NMRDistanceOffset =
    (((distances + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits);
  gpu->sim.NMRAngleOffset = gpu->sim.NMRDistanceOffset +
    (((angles + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits);
  gpu->sim.NMRTorsionOffset = gpu->sim.NMRAngleOffset +
    (((torsions + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits);
  gpu->sim.NMRCOMDistanceOffset = gpu->sim.NMRTorsionOffset +
    (((COMdistances  + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits);
  gpu->sim.NMRr6avDistanceOffset = gpu->sim.NMRCOMDistanceOffset +
    (((r6avdistances + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits);
  gpu->sim.NMRCOMAngleOffset = gpu->sim.NMRr6avDistanceOffset +
    (((COMangles + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits);
  gpu->sim.NMRCOMTorsionOffset = gpu->sim.NMRCOMAngleOffset +
    (((COMtorsions + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits);
  gpu->sim.NMRMaxgrp = *maxgrp;

  // Delete existing restraints
  delete gpu->pbNMRJarData;
  delete gpu->pbNMRDistanceID;
  delete gpu->pbNMRDistanceR1R2;
  delete gpu->pbNMRDistanceR3R4;
  delete gpu->pbNMRDistanceK2K3;
  delete gpu->pbNMRDistanceK4;
  delete gpu->pbNMRDistanceAve;
  delete gpu->pbNMRDistanceTgtVal;
  delete gpu->pbNMRDistanceStep;
  delete gpu->pbNMRDistanceInc;
  delete gpu->pbNMRDistanceR1R2Slp;
  delete gpu->pbNMRDistanceR3R4Slp;
  delete gpu->pbNMRDistanceK2K3Slp;
  delete gpu->pbNMRDistanceK4Slp;
  delete gpu->pbNMRDistanceR1R2Int;
  delete gpu->pbNMRDistanceR3R4Int;
  delete gpu->pbNMRDistanceK2K3Int;
  delete gpu->pbNMRDistanceK4Int;

  delete gpu->pbNMRCOMDistanceID;
  delete gpu->pbNMRCOMDistanceCOM;
  delete gpu->pbNMRCOMDistanceCOMGrp;
  delete gpu->pbNMRCOMDistanceR1R2;
  delete gpu->pbNMRCOMDistanceR3R4;
  delete gpu->pbNMRCOMDistanceK2K3;
  delete gpu->pbNMRCOMDistanceK4;
  delete gpu->pbNMRCOMDistanceAve;
  delete gpu->pbNMRCOMDistanceTgtVal;
  delete gpu->pbNMRCOMDistanceStep;
  delete gpu->pbNMRCOMDistanceInc;
  delete gpu->pbNMRCOMDistanceR1R2Slp;
  delete gpu->pbNMRCOMDistanceR3R4Slp;
  delete gpu->pbNMRCOMDistanceK2K3Slp;
  delete gpu->pbNMRCOMDistanceK4Slp;
  delete gpu->pbNMRCOMDistanceR1R2Int;
  delete gpu->pbNMRCOMDistanceR3R4Int;
  delete gpu->pbNMRCOMDistanceK2K3Int;
  delete gpu->pbNMRCOMDistanceK4Int;
  delete gpu->pbNMRCOMDistanceWeights;
  delete gpu->pbNMRCOMDistanceXYZ;

  delete gpu->pbNMRr6avDistanceID;
  delete gpu->pbNMRr6avDistancer6av;
  delete gpu->pbNMRr6avDistancer6avGrp;
  delete gpu->pbNMRr6avDistanceR1R2;
  delete gpu->pbNMRr6avDistanceR3R4;
  delete gpu->pbNMRr6avDistanceK2K3;
  delete gpu->pbNMRr6avDistanceK4;
  delete gpu->pbNMRr6avDistanceAve;
  delete gpu->pbNMRr6avDistanceTgtVal;
  delete gpu->pbNMRr6avDistanceStep;
  delete gpu->pbNMRr6avDistanceInc;
  delete gpu->pbNMRr6avDistanceR1R2Slp;
  delete gpu->pbNMRr6avDistanceR3R4Slp;
  delete gpu->pbNMRr6avDistanceK2K3Slp;
  delete gpu->pbNMRr6avDistanceK4Slp;
  delete gpu->pbNMRr6avDistanceR1R2Int;
  delete gpu->pbNMRr6avDistanceR3R4Int;
  delete gpu->pbNMRr6avDistanceK2K3Int;
  delete gpu->pbNMRr6avDistanceK4Int;
  delete gpu->pbNMRr6avDistanceWeights;
  delete gpu->pbNMRr6avDistanceXYZ;

  delete gpu->pbNMRAngleID1;
  delete gpu->pbNMRAngleID2;
  delete gpu->pbNMRAngleR1R2;
  delete gpu->pbNMRAngleR3R4;
  delete gpu->pbNMRAngleK2K3;
  delete gpu->pbNMRAngleK4;
  delete gpu->pbNMRAngleAve;
  delete gpu->pbNMRAngleTgtVal;
  delete gpu->pbNMRAngleStep;
  delete gpu->pbNMRAngleInc;
  delete gpu->pbNMRAngleR1R2Slp;
  delete gpu->pbNMRAngleR3R4Slp;
  delete gpu->pbNMRAngleK2K3Slp;
  delete gpu->pbNMRAngleK4Slp;
  delete gpu->pbNMRAngleR1R2Int;
  delete gpu->pbNMRAngleR3R4Int;
  delete gpu->pbNMRAngleK2K3Int;
  delete gpu->pbNMRAngleK4Int;

  delete gpu->pbNMRCOMAngleID1;
  delete gpu->pbNMRCOMAngleID2;
  delete gpu->pbNMRCOMAngleCOM;
  delete gpu->pbNMRCOMAngleCOMGrp;
  delete gpu->pbNMRCOMAngleR1R2;
  delete gpu->pbNMRCOMAngleR3R4;
  delete gpu->pbNMRCOMAngleK2K3;
  delete gpu->pbNMRCOMAngleK4;
  delete gpu->pbNMRCOMAngleAve;
  delete gpu->pbNMRCOMAngleTgtVal;
  delete gpu->pbNMRCOMAngleStep;
  delete gpu->pbNMRCOMAngleInc;
  delete gpu->pbNMRCOMAngleR1R2Slp;
  delete gpu->pbNMRCOMAngleR3R4Slp;
  delete gpu->pbNMRCOMAngleK2K3Slp;
  delete gpu->pbNMRCOMAngleK4Slp;
  delete gpu->pbNMRCOMAngleR1R2Int;
  delete gpu->pbNMRCOMAngleR3R4Int;
  delete gpu->pbNMRCOMAngleK2K3Int;
  delete gpu->pbNMRCOMAngleK4Int;

  delete gpu->pbNMRTorsionID1;
  delete gpu->pbNMRTorsionR1R2;
  delete gpu->pbNMRTorsionR3R4;
  delete gpu->pbNMRTorsionK2K3;
  delete gpu->pbNMRTorsionK4;
  delete gpu->pbNMRTorsionAve1;
  delete gpu->pbNMRTorsionAve2;
  delete gpu->pbNMRTorsionTgtVal;
  delete gpu->pbNMRTorsionStep;
  delete gpu->pbNMRTorsionInc;
  delete gpu->pbNMRTorsionR1R2Slp;
  delete gpu->pbNMRTorsionR3R4Slp;
  delete gpu->pbNMRTorsionK2K3Slp;
  delete gpu->pbNMRTorsionK4Slp;
  delete gpu->pbNMRTorsionR1R2Int;
  delete gpu->pbNMRTorsionR3R4Int;
  delete gpu->pbNMRTorsionK2K3Int;
  delete gpu->pbNMRTorsionK4Int;

  delete gpu->pbNMRCOMTorsionID1;
  delete gpu->pbNMRCOMTorsionCOM;
  delete gpu->pbNMRCOMTorsionCOMGrp;
  delete gpu->pbNMRCOMTorsionR1R2;
  delete gpu->pbNMRCOMTorsionR3R4;
  delete gpu->pbNMRCOMTorsionK2K3;
  delete gpu->pbNMRCOMTorsionK4;
  delete gpu->pbNMRCOMTorsionAve;
  delete gpu->pbNMRCOMTorsionTgtVal;
  delete gpu->pbNMRCOMTorsionStep;
  delete gpu->pbNMRCOMTorsionInc;
  delete gpu->pbNMRCOMTorsionR1R2Slp;
  delete gpu->pbNMRCOMTorsionR3R4Slp;
  delete gpu->pbNMRCOMTorsionK2K3Slp;
  delete gpu->pbNMRCOMTorsionK4Slp;
  delete gpu->pbNMRCOMTorsionR1R2Int;
  delete gpu->pbNMRCOMTorsionR3R4Int;
  delete gpu->pbNMRCOMTorsionK2K3Int;
  delete gpu->pbNMRCOMTorsionK4Int;

  // Allocate restraint info
  gpu->pbNMRJarData             = new GpuBuffer<double>(5);
  gpu->pbNMRDistanceID          = new GpuBuffer<int2>(distances);
  gpu->pbNMRDistanceR1R2        = new GpuBuffer<PMEDouble2>(distances);
  gpu->pbNMRDistanceR3R4        = new GpuBuffer<PMEDouble2>(distances);
  gpu->pbNMRDistanceK2K3        = new GpuBuffer<PMEDouble2>(distances);
  gpu->pbNMRDistanceK4          = new GpuBuffer<PMEDouble>(distances);
  gpu->pbNMRDistanceAve         = new GpuBuffer<PMEDouble>(distances);
  gpu->pbNMRDistanceTgtVal      = new GpuBuffer<PMEDouble2>(distances);
  gpu->pbNMRDistanceStep        = new GpuBuffer<int2>(distances);
  gpu->pbNMRDistanceInc         = new GpuBuffer<int>(distances);
  gpu->pbNMRDistanceR1R2Slp     = new GpuBuffer<PMEDouble2>(distances);
  gpu->pbNMRDistanceR3R4Slp     = new GpuBuffer<PMEDouble2>(distances);
  gpu->pbNMRDistanceK2K3Slp     = new GpuBuffer<PMEDouble2>(distances);
  gpu->pbNMRDistanceK4Slp       = new GpuBuffer<PMEDouble>(distances);
  gpu->pbNMRDistanceR1R2Int     = new GpuBuffer<PMEDouble2>(distances);
  gpu->pbNMRDistanceR3R4Int     = new GpuBuffer<PMEDouble2>(distances);
  gpu->pbNMRDistanceK2K3Int     = new GpuBuffer<PMEDouble2>(distances);
  gpu->pbNMRDistanceK4Int       = new GpuBuffer<PMEDouble>(distances);

  gpu->pbNMRCOMDistanceID       = new GpuBuffer<int2>(COMdistances);
  gpu->pbNMRCOMDistanceCOM      = new GpuBuffer<int2>((*maxgrp) * (COMdistances));
  gpu->pbNMRCOMDistanceCOMGrp   = new GpuBuffer<int2>(COMdistances * 2);
  gpu->pbNMRCOMDistanceR1R2     = new GpuBuffer<PMEDouble2>(COMdistances);
  gpu->pbNMRCOMDistanceR3R4     = new GpuBuffer<PMEDouble2>(COMdistances);
  gpu->pbNMRCOMDistanceK2K3     = new GpuBuffer<PMEDouble2>(COMdistances);
  gpu->pbNMRCOMDistanceK4       = new GpuBuffer<PMEDouble>(COMdistances);
  gpu->pbNMRCOMDistanceAve      = new GpuBuffer<PMEDouble>(COMdistances);
  gpu->pbNMRCOMDistanceTgtVal   = new GpuBuffer<PMEDouble2>(COMdistances);
  gpu->pbNMRCOMDistanceStep     = new GpuBuffer<int2>(COMdistances);
  gpu->pbNMRCOMDistanceInc      = new GpuBuffer<int>(COMdistances);
  gpu->pbNMRCOMDistanceWeights  = new GpuBuffer<int>(COMdistances*5);
  gpu->pbNMRCOMDistanceXYZ      = new GpuBuffer<PMEDouble>(COMdistances*3);
  gpu->pbNMRCOMDistanceR1R2Slp  = new GpuBuffer<PMEDouble2>(COMdistances);
  gpu->pbNMRCOMDistanceR3R4Slp  = new GpuBuffer<PMEDouble2>(COMdistances);
  gpu->pbNMRCOMDistanceK2K3Slp  = new GpuBuffer<PMEDouble2>(COMdistances);
  gpu->pbNMRCOMDistanceK4Slp    = new GpuBuffer<PMEDouble>(COMdistances);
  gpu->pbNMRCOMDistanceR1R2Int  = new GpuBuffer<PMEDouble2>(COMdistances);
  gpu->pbNMRCOMDistanceR3R4Int  = new GpuBuffer<PMEDouble2>(COMdistances);
  gpu->pbNMRCOMDistanceK2K3Int  = new GpuBuffer<PMEDouble2>(COMdistances);
  gpu->pbNMRCOMDistanceK4Int    = new GpuBuffer<PMEDouble>(COMdistances);

  gpu->pbNMRr6avDistanceID      = new GpuBuffer<int2>(r6avdistances);
  gpu->pbNMRr6avDistancer6av    = new GpuBuffer<int2>(*maxgrp * r6avdistances);
  gpu->pbNMRr6avDistancer6avGrp = new GpuBuffer<int2>(r6avdistances * 2);
  gpu->pbNMRr6avDistanceR1R2    = new GpuBuffer<PMEDouble2>(r6avdistances);
  gpu->pbNMRr6avDistanceR3R4    = new GpuBuffer<PMEDouble2>(r6avdistances);
  gpu->pbNMRr6avDistanceK2K3    = new GpuBuffer<PMEDouble2>(r6avdistances);
  gpu->pbNMRr6avDistanceK4      = new GpuBuffer<PMEDouble>(r6avdistances);
  gpu->pbNMRr6avDistanceAve     = new GpuBuffer<PMEDouble>(r6avdistances);
  gpu->pbNMRr6avDistanceTgtVal  = new GpuBuffer<PMEDouble2>(r6avdistances);
  gpu->pbNMRr6avDistanceStep    = new GpuBuffer<int2>(r6avdistances);
  gpu->pbNMRr6avDistanceInc     = new GpuBuffer<int>(r6avdistances);
  gpu->pbNMRr6avDistanceWeights  = new GpuBuffer<int>(r6avdistances*5);
  gpu->pbNMRr6avDistanceXYZ      = new GpuBuffer<PMEDouble>(r6avdistances*3);
  gpu->pbNMRr6avDistanceR1R2Slp = new GpuBuffer<PMEDouble2>(r6avdistances);
  gpu->pbNMRr6avDistanceR3R4Slp = new GpuBuffer<PMEDouble2>(r6avdistances);
  gpu->pbNMRr6avDistanceK2K3Slp = new GpuBuffer<PMEDouble2>(r6avdistances);
  gpu->pbNMRr6avDistanceK4Slp   = new GpuBuffer<PMEDouble>(r6avdistances);
  gpu->pbNMRr6avDistanceR1R2Int = new GpuBuffer<PMEDouble2>(r6avdistances);
  gpu->pbNMRr6avDistanceR3R4Int = new GpuBuffer<PMEDouble2>(r6avdistances);
  gpu->pbNMRr6avDistanceK2K3Int = new GpuBuffer<PMEDouble2>(r6avdistances);
  gpu->pbNMRr6avDistanceK4Int   = new GpuBuffer<PMEDouble>(r6avdistances);

  gpu->pbNMRAngleID1       = new GpuBuffer<int2>(angles);
  gpu->pbNMRAngleID2       = new GpuBuffer<int>(angles);
  gpu->pbNMRAngleR1R2      = new GpuBuffer<PMEDouble2>(angles);
  gpu->pbNMRAngleR3R4      = new GpuBuffer<PMEDouble2>(angles);
  gpu->pbNMRAngleK2K3      = new GpuBuffer<PMEDouble2>(angles);
  gpu->pbNMRAngleK4        = new GpuBuffer<PMEDouble>(angles);
  gpu->pbNMRAngleAve       = new GpuBuffer<PMEDouble>(angles);
  gpu->pbNMRAngleTgtVal    = new GpuBuffer<PMEDouble2>(angles);
  gpu->pbNMRAngleStep      = new GpuBuffer<int2>(angles);
  gpu->pbNMRAngleInc       = new GpuBuffer<int>(angles);
  gpu->pbNMRAngleR1R2Slp   = new GpuBuffer<PMEDouble2>(angles);
  gpu->pbNMRAngleR3R4Slp   = new GpuBuffer<PMEDouble2>(angles);
  gpu->pbNMRAngleK2K3Slp   = new GpuBuffer<PMEDouble2>(angles);
  gpu->pbNMRAngleK4Slp     = new GpuBuffer<PMEDouble>(angles);
  gpu->pbNMRAngleR1R2Int   = new GpuBuffer<PMEDouble2>(angles);
  gpu->pbNMRAngleR3R4Int   = new GpuBuffer<PMEDouble2>(angles);
  gpu->pbNMRAngleK2K3Int   = new GpuBuffer<PMEDouble2>(angles);
  gpu->pbNMRAngleK4Int     = new GpuBuffer<PMEDouble>(angles);

  gpu->pbNMRCOMAngleID1      = new GpuBuffer<int2>(COMangles);
  gpu->pbNMRCOMAngleID2      = new GpuBuffer<int>(COMangles);
  gpu->pbNMRCOMAngleCOM      = new GpuBuffer<int2>(*maxgrp * COMangles);
  gpu->pbNMRCOMAngleCOMGrp   = new GpuBuffer<int2>(COMangles * 3);
  gpu->pbNMRCOMAngleR1R2     = new GpuBuffer<PMEDouble2>(COMangles);
  gpu->pbNMRCOMAngleR3R4     = new GpuBuffer<PMEDouble2>(COMangles);
  gpu->pbNMRCOMAngleK2K3     = new GpuBuffer<PMEDouble2>(COMangles);
  gpu->pbNMRCOMAngleK4       = new GpuBuffer<PMEDouble>(COMangles);
  gpu->pbNMRCOMAngleAve      = new GpuBuffer<PMEDouble>(COMangles);
  gpu->pbNMRCOMAngleTgtVal   = new GpuBuffer<PMEDouble2>(COMangles);
  gpu->pbNMRCOMAngleStep     = new GpuBuffer<int2>(COMangles);
  gpu->pbNMRCOMAngleInc      = new GpuBuffer<int>(COMangles);
  gpu->pbNMRCOMAngleR1R2Slp  = new GpuBuffer<PMEDouble2>(COMangles);
  gpu->pbNMRCOMAngleR3R4Slp  = new GpuBuffer<PMEDouble2>(COMangles);
  gpu->pbNMRCOMAngleK2K3Slp  = new GpuBuffer<PMEDouble2>(COMangles);
  gpu->pbNMRCOMAngleK4Slp    = new GpuBuffer<PMEDouble>(COMangles);
  gpu->pbNMRCOMAngleR1R2Int  = new GpuBuffer<PMEDouble2>(COMangles);
  gpu->pbNMRCOMAngleR3R4Int  = new GpuBuffer<PMEDouble2>(COMangles);
  gpu->pbNMRCOMAngleK2K3Int  = new GpuBuffer<PMEDouble2>(COMangles);
  gpu->pbNMRCOMAngleK4Int    = new GpuBuffer<PMEDouble>(COMangles);

  gpu->pbNMRTorsionID1     = new GpuBuffer<int4>(torsions);
  gpu->pbNMRTorsionR1R2    = new GpuBuffer<PMEDouble2>(torsions);
  gpu->pbNMRTorsionR3R4    = new GpuBuffer<PMEDouble2>(torsions);
  gpu->pbNMRTorsionK2K3    = new GpuBuffer<PMEDouble2>(torsions);
  gpu->pbNMRTorsionK4      = new GpuBuffer<PMEDouble>(torsions);
  gpu->pbNMRTorsionAve1    = new GpuBuffer<PMEDouble>(torsions);
  gpu->pbNMRTorsionAve2    = new GpuBuffer<PMEDouble>(torsions);
  gpu->pbNMRTorsionTgtVal  = new GpuBuffer<PMEDouble2>(torsions);
  gpu->pbNMRTorsionStep    = new GpuBuffer<int2>(torsions);
  gpu->pbNMRTorsionInc     = new GpuBuffer<int>(torsions);
  gpu->pbNMRTorsionR1R2Slp = new GpuBuffer<PMEDouble2>(torsions);
  gpu->pbNMRTorsionR3R4Slp = new GpuBuffer<PMEDouble2>(torsions);
  gpu->pbNMRTorsionK2K3Slp = new GpuBuffer<PMEDouble2>(torsions);
  gpu->pbNMRTorsionK4Slp   = new GpuBuffer<PMEDouble>(torsions);
  gpu->pbNMRTorsionR1R2Int = new GpuBuffer<PMEDouble2>(torsions);
  gpu->pbNMRTorsionR3R4Int = new GpuBuffer<PMEDouble2>(torsions);
  gpu->pbNMRTorsionK2K3Int = new GpuBuffer<PMEDouble2>(torsions);
  gpu->pbNMRTorsionK4Int   = new GpuBuffer<PMEDouble>(torsions);

  gpu->pbNMRCOMTorsionID1      = new GpuBuffer<int4>(COMtorsions);
  gpu->pbNMRCOMTorsionCOM      = new GpuBuffer<int2>(*maxgrp * COMtorsions);
  gpu->pbNMRCOMTorsionCOMGrp   = new GpuBuffer<int2>(COMtorsions * 4);
  gpu->pbNMRCOMTorsionR1R2     = new GpuBuffer<PMEDouble2>(COMtorsions);
  gpu->pbNMRCOMTorsionR3R4     = new GpuBuffer<PMEDouble2>(COMtorsions);
  gpu->pbNMRCOMTorsionK2K3     = new GpuBuffer<PMEDouble2>(COMtorsions);
  gpu->pbNMRCOMTorsionK4       = new GpuBuffer<PMEDouble>(COMtorsions);
  gpu->pbNMRCOMTorsionAve      = new GpuBuffer<PMEDouble>(COMtorsions);
  gpu->pbNMRCOMTorsionTgtVal   = new GpuBuffer<PMEDouble2>(COMtorsions);
  gpu->pbNMRCOMTorsionStep     = new GpuBuffer<int2>(COMtorsions);
  gpu->pbNMRCOMTorsionInc      = new GpuBuffer<int>(COMtorsions);
  gpu->pbNMRCOMTorsionR1R2Slp  = new GpuBuffer<PMEDouble2>(COMtorsions);
  gpu->pbNMRCOMTorsionR3R4Slp  = new GpuBuffer<PMEDouble2>(COMtorsions);
  gpu->pbNMRCOMTorsionK2K3Slp  = new GpuBuffer<PMEDouble2>(COMtorsions);
  gpu->pbNMRCOMTorsionK4Slp    = new GpuBuffer<PMEDouble>(COMtorsions);
  gpu->pbNMRCOMTorsionR1R2Int  = new GpuBuffer<PMEDouble2>(COMtorsions);
  gpu->pbNMRCOMTorsionR3R4Int  = new GpuBuffer<PMEDouble2>(COMtorsions);
  gpu->pbNMRCOMTorsionK2K3Int  = new GpuBuffer<PMEDouble2>(COMtorsions);
  gpu->pbNMRCOMTorsionK4Int    = new GpuBuffer<PMEDouble>(COMtorsions);

  // Clear Jarzynski data
  gpu->sim.bJar = (*jar == 1);
  gpu->sim.drjar = *drjar;
  for (int i = 0; i < 5; i++) {
    gpu->pbNMRJarData->_pSysData[i] = 0.0;
  }

  // Add restraints
  distances = 0;
  COMdistances = 0;
  r6avdistances = 0;
  angles = 0;
  torsions = 0;
  COMangles = 0;
  COMtorsions = 0;
  for (int i = 0; i < *nmrnum; i++) {
    if (resttype[i] == 1) {
      if (nmrat[i][0] >= 0 && nmrat[i][1] >= 0) {
        gpu->pbNMRDistanceID->_pSysData[distances].x = nmrat[i][0] / 3;
        gpu->pbNMRDistanceID->_pSysData[distances].y = nmrat[i][1] / 3;
        gpu->pbNMRDistanceR1R2->_pSysData[distances].x = r1nmr[i][1];
        gpu->pbNMRDistanceR1R2->_pSysData[distances].y = r2nmr[i][1];
        gpu->pbNMRDistanceR3R4->_pSysData[distances].x = r3nmr[i][1];
        gpu->pbNMRDistanceR3R4->_pSysData[distances].y = r4nmr[i][1];
        gpu->pbNMRDistanceK2K3->_pSysData[distances].x = rk2nmr[i][1];
        gpu->pbNMRDistanceK2K3->_pSysData[distances].y = rk3nmr[i][1];
        gpu->pbNMRDistanceAve->_pSysData[distances] = (PMEDouble)0.0;
        gpu->pbNMRDistanceStep->_pSysData[distances].x = nmrst[i][0];
        gpu->pbNMRDistanceStep->_pSysData[distances].y = nmrst[i][1];
        gpu->pbNMRDistanceInc->_pSysData[distances] = nmrst[i][2];
        gpu->pbNMRDistanceR1R2Slp->_pSysData[distances].x = r1nmr[i][0];
        gpu->pbNMRDistanceR1R2Slp->_pSysData[distances].y = r2nmr[i][0];
        gpu->pbNMRDistanceR3R4Slp->_pSysData[distances].x = r3nmr[i][0];
        gpu->pbNMRDistanceR3R4Slp->_pSysData[distances].y = r4nmr[i][0];
        gpu->pbNMRDistanceK2K3Slp->_pSysData[distances].x = rk2nmr[i][0];
        gpu->pbNMRDistanceK2K3Slp->_pSysData[distances].y = rk3nmr[i][0];
        gpu->pbNMRDistanceR1R2Int->_pSysData[distances].x = r1nmr[i][1];
        gpu->pbNMRDistanceR1R2Int->_pSysData[distances].y = r2nmr[i][1];
        gpu->pbNMRDistanceR3R4Int->_pSysData[distances].x = r3nmr[i][1];
        gpu->pbNMRDistanceR3R4Int->_pSysData[distances].y = r4nmr[i][1];
        gpu->pbNMRDistanceK2K3Int->_pSysData[distances].x = rk2nmr[i][1];
        gpu->pbNMRDistanceK2K3Int->_pSysData[distances].y = rk3nmr[i][1];
        distances++;
      }
      else if (igravt[i] == 1) {
        if (nmrat[i][0] >= 0) {
          gpu->pbNMRr6avDistanceID->_pSysData[r6avdistances].x = nmrat[i][0] / 3;
        }
        else {

          // This code is for r6av restraints. By keeping it a multiple of 3 the number of
          // atoms in the group of atoms can be indexed in the same manner as the CPU code.
          // nmrcom stores the atom ranges ivolved in r6av restraints.
          gpu->pbNMRr6avDistanceID->_pSysData[r6avdistances].x = nmrat[i][0];
          gpu->pbNMRr6avDistancer6avGrp->_pSysData[r6avdistances * 2].x = -nmrat[i][0]-1;
          gpu->pbNMRr6avDistancer6avGrp->_pSysData[r6avdistances * 2].y = -nmrat[i][8];
          for (int ip1 = -nmrat[i][0]-1 ; ip1 < -nmrat[i][8] ; ip1++) {
            gpu->pbNMRr6avDistancer6av->_pSysData[ip1].x = nmrcom[ip1][0] - 1;
            gpu->pbNMRr6avDistancer6av->_pSysData[ip1].y = nmrcom[ip1][1] - 1;
          }
        }
        if (nmrat[i][1] >= 0) {
          gpu->pbNMRr6avDistanceID->_pSysData[r6avdistances].y = nmrat[i][1] / 3;
        }
        else {

          // This code is for r6av restraints. By keeping it a multiple of 3 the number of
          // atoms in the group of atoms can be indexed in the same manner as the CPU code.
          // nmrcom stores the atom ranges ivolved in r6av restraints.
          gpu->pbNMRr6avDistanceID->_pSysData[r6avdistances].y = nmrat[i][1];
          gpu->pbNMRr6avDistancer6avGrp->_pSysData[r6avdistances*2 + 1].x = -nmrat[i][1]-1;
          gpu->pbNMRr6avDistancer6avGrp->_pSysData[r6avdistances*2 + 1].y = -nmrat[i][9];
          for (int ip2 = -nmrat[i][1]-1; ip2 < -nmrat[i][9]; ip2++) {
            gpu->pbNMRr6avDistancer6av->_pSysData[ip2].x = nmrcom[ip2][0] - 1;
            gpu->pbNMRr6avDistancer6av->_pSysData[ip2].y = nmrcom[ip2][1] - 1;
          }
        }
        gpu->pbNMRr6avDistanceR1R2->_pSysData[r6avdistances].x = r1nmr[i][1];
        gpu->pbNMRr6avDistanceR1R2->_pSysData[r6avdistances].y = r2nmr[i][1];
        gpu->pbNMRr6avDistanceR3R4->_pSysData[r6avdistances].x = r3nmr[i][1];
        gpu->pbNMRr6avDistanceR3R4->_pSysData[r6avdistances].y = r4nmr[i][1];
        gpu->pbNMRr6avDistanceK2K3->_pSysData[r6avdistances].x = rk2nmr[i][1];
        gpu->pbNMRr6avDistanceK2K3->_pSysData[r6avdistances].y = rk3nmr[i][1];
        gpu->pbNMRr6avDistanceAve->_pSysData[r6avdistances] = (PMEDouble)0.0;
        gpu->pbNMRr6avDistanceStep->_pSysData[r6avdistances].x = nmrst[i][0];
        gpu->pbNMRr6avDistanceStep->_pSysData[r6avdistances].y = nmrst[i][1];
        gpu->pbNMRr6avDistanceInc->_pSysData[r6avdistances] = nmrst[i][2];
        gpu->pbNMRr6avDistanceR1R2Slp->_pSysData[r6avdistances].x = r1nmr[i][0];
        gpu->pbNMRr6avDistanceR1R2Slp->_pSysData[r6avdistances].y = r2nmr[i][0];
        gpu->pbNMRr6avDistanceR3R4Slp->_pSysData[r6avdistances].x = r3nmr[i][0];
        gpu->pbNMRr6avDistanceR3R4Slp->_pSysData[r6avdistances].y = r4nmr[i][0];
        gpu->pbNMRr6avDistanceK2K3Slp->_pSysData[r6avdistances].x = rk2nmr[i][0];
        gpu->pbNMRr6avDistanceK2K3Slp->_pSysData[r6avdistances].y = rk3nmr[i][0];
        gpu->pbNMRr6avDistanceR1R2Int->_pSysData[r6avdistances].x = r1nmr[i][1];
        gpu->pbNMRr6avDistanceR1R2Int->_pSysData[r6avdistances].y = r2nmr[i][1];
        gpu->pbNMRr6avDistanceR3R4Int->_pSysData[r6avdistances].x = r3nmr[i][1];
        gpu->pbNMRr6avDistanceR3R4Int->_pSysData[r6avdistances].y = r4nmr[i][1];
        gpu->pbNMRr6avDistanceK2K3Int->_pSysData[r6avdistances].x = rk2nmr[i][1];
        gpu->pbNMRr6avDistanceK2K3Int->_pSysData[r6avdistances].y = rk3nmr[i][1];
        gpu->pbNMRr6avDistanceWeights->_pSysData[r6avdistances*5] = ifxyz[i*5];
        gpu->pbNMRr6avDistanceWeights->_pSysData[5*r6avdistances+1] = ifxyz[5*i+1];
        gpu->pbNMRr6avDistanceWeights->_pSysData[5*r6avdistances+2] = ifxyz[5*i+2];
        gpu->pbNMRr6avDistanceWeights->_pSysData[5*r6avdistances+3] = ifxyz[5*i+3];
        gpu->pbNMRr6avDistanceWeights->_pSysData[5*r6avdistances+4] = ifxyz[5*i+4];
        gpu->pbNMRr6avDistanceXYZ->_pSysData[r6avdistances*3] = idxyz[i*3];
        gpu->pbNMRr6avDistanceXYZ->_pSysData[3*r6avdistances+1] = idxyz[3*i+1];
        gpu->pbNMRr6avDistanceXYZ->_pSysData[3*r6avdistances+2] = idxyz[3*i+2];
        r6avdistances++;
      }
      else {  //com distances
        if (nmrat[i][0] >= 0) { //check if only one atom
          gpu->pbNMRCOMDistanceID->_pSysData[COMdistances].x = nmrat[i][0] / 3;
        }
        else {

          // This code is for COM restraints. By keeping it a multiple of 3 the number of
          // atoms in the group of atoms can be indexed in the same manner as the CPU code.
          // nmrcom stores the atom ranges ivolved in COM restraints.
          // last index is stored in nrmat[9-1, 10-1]
          gpu->pbNMRCOMDistanceID->_pSysData[COMdistances].x = nmrat[i][0];
          gpu->pbNMRCOMDistanceCOMGrp->_pSysData[COMdistances * 2].x = -nmrat[i][0]-1;
          gpu->pbNMRCOMDistanceCOMGrp->_pSysData[COMdistances * 2].y = -nmrat[i][8];
          for (int ip1 = -nmrat[i][0]-1 ; ip1 < -nmrat[i][8] ; ip1++) {
            gpu->pbNMRCOMDistanceCOM->_pSysData[ip1].x = nmrcom[ip1][0] - 1;
            gpu->pbNMRCOMDistanceCOM->_pSysData[ip1].y = nmrcom[ip1][1] - 1;
          }
        }
        if (nmrat[i][1] >= 0) {
          gpu->pbNMRCOMDistanceID->_pSysData[COMdistances].y = nmrat[i][1] / 3;
        }
        else {

          // This code is for COM restraints. By keeping it a multiple of 3 the number of
          // atoms in the group of atoms can be indexed in the same manner as the CPU code.
          // nmrcom stores the atom ranges ivolved in COM restraints.
          gpu->pbNMRCOMDistanceID->_pSysData[COMdistances].y = nmrat[i][1];
          gpu->pbNMRCOMDistanceCOMGrp->_pSysData[COMdistances * 2 + 1].x = -nmrat[i][1]-1;
          gpu->pbNMRCOMDistanceCOMGrp->_pSysData[COMdistances * 2 + 1].y = -nmrat[i][9];
          for (int ip2 = -nmrat[i][1]-1; ip2 < -nmrat[i][9]; ip2++) {
            gpu->pbNMRCOMDistanceCOM->_pSysData[ip2].x = nmrcom[ip2][0] - 1;
            gpu->pbNMRCOMDistanceCOM->_pSysData[ip2].y = nmrcom[ip2][1] - 1;
          }
        }
        gpu->pbNMRCOMDistanceR1R2->_pSysData[COMdistances].x = r1nmr[i][1];
        gpu->pbNMRCOMDistanceR1R2->_pSysData[COMdistances].y = r2nmr[i][1];
        gpu->pbNMRCOMDistanceR3R4->_pSysData[COMdistances].x = r3nmr[i][1];
        gpu->pbNMRCOMDistanceR3R4->_pSysData[COMdistances].y = r4nmr[i][1];
        gpu->pbNMRCOMDistanceK2K3->_pSysData[COMdistances].x = rk2nmr[i][1];
        gpu->pbNMRCOMDistanceK2K3->_pSysData[COMdistances].y = rk3nmr[i][1];
        gpu->pbNMRCOMDistanceAve->_pSysData[COMdistances] = (PMEDouble)0.0;
        gpu->pbNMRCOMDistanceStep->_pSysData[COMdistances].x = nmrst[i][0];
        gpu->pbNMRCOMDistanceStep->_pSysData[COMdistances].y = nmrst[i][1];
        gpu->pbNMRCOMDistanceInc->_pSysData[COMdistances] = nmrst[i][2];
        gpu->pbNMRCOMDistanceR1R2Slp->_pSysData[COMdistances].x = r1nmr[i][0];
        gpu->pbNMRCOMDistanceR1R2Slp->_pSysData[COMdistances].y = r2nmr[i][0];
        gpu->pbNMRCOMDistanceR3R4Slp->_pSysData[COMdistances].x = r3nmr[i][0];
        gpu->pbNMRCOMDistanceR3R4Slp->_pSysData[COMdistances].y = r4nmr[i][0];
        gpu->pbNMRCOMDistanceK2K3Slp->_pSysData[COMdistances].x = rk2nmr[i][0];
        gpu->pbNMRCOMDistanceK2K3Slp->_pSysData[COMdistances].y = rk3nmr[i][0];
        gpu->pbNMRCOMDistanceR1R2Int->_pSysData[COMdistances].x = r1nmr[i][1];
        gpu->pbNMRCOMDistanceR1R2Int->_pSysData[COMdistances].y = r2nmr[i][1];
        gpu->pbNMRCOMDistanceR3R4Int->_pSysData[COMdistances].x = r3nmr[i][1];
        gpu->pbNMRCOMDistanceR3R4Int->_pSysData[COMdistances].y = r4nmr[i][1];
        gpu->pbNMRCOMDistanceK2K3Int->_pSysData[COMdistances].x = rk2nmr[i][1];
        gpu->pbNMRCOMDistanceK2K3Int->_pSysData[COMdistances].y = rk3nmr[i][1];
        gpu->pbNMRCOMDistanceWeights->_pSysData[COMdistances*5] = ifxyz[i*5];
        gpu->pbNMRCOMDistanceWeights->_pSysData[5*COMdistances+1] = ifxyz[5*i+1];
        gpu->pbNMRCOMDistanceWeights->_pSysData[5*COMdistances+2] = ifxyz[5*i+2];
        gpu->pbNMRCOMDistanceWeights->_pSysData[5*COMdistances+3] = ifxyz[5*i+3];
        gpu->pbNMRCOMDistanceWeights->_pSysData[5*COMdistances+4] = ifxyz[5*i+4];
        gpu->pbNMRCOMDistanceXYZ->_pSysData[COMdistances*3] = idxyz[i*3];
        gpu->pbNMRCOMDistanceXYZ->_pSysData[3*COMdistances+1] = idxyz[3*i+1];
        gpu->pbNMRCOMDistanceXYZ->_pSysData[3*COMdistances+2] = idxyz[3*i+2];
        COMdistances++;
      }
    }
    else if (resttype[i] == 2) {
      if (nmrat[i][0] >= 0 && nmrat[i][1] >= 0 && nmrat[i][2] >= 0) {
        gpu->pbNMRAngleID1->_pSysData[angles].x = nmrat[i][0] / 3;
        gpu->pbNMRAngleID1->_pSysData[angles].y = nmrat[i][1] / 3;
        gpu->pbNMRAngleID2->_pSysData[angles] = nmrat[i][2] / 3;
        gpu->pbNMRAngleR1R2->_pSysData[angles].x = r1nmr[i][1];
        gpu->pbNMRAngleR1R2->_pSysData[angles].y = r2nmr[i][1];
        gpu->pbNMRAngleR3R4->_pSysData[angles].x = r3nmr[i][1];
        gpu->pbNMRAngleR3R4->_pSysData[angles].y = r4nmr[i][1];
        gpu->pbNMRAngleK2K3->_pSysData[angles].x = rk2nmr[i][1];
        gpu->pbNMRAngleK2K3->_pSysData[angles].y = rk3nmr[i][1];
        gpu->pbNMRAngleAve->_pSysData[angles] = (PMEDouble)0.0;
        gpu->pbNMRAngleStep->_pSysData[angles].x = nmrst[i][0];
        gpu->pbNMRAngleStep->_pSysData[angles].y = nmrst[i][1];
        gpu->pbNMRAngleInc->_pSysData[angles] = nmrst[i][2];
        gpu->pbNMRAngleR1R2Slp->_pSysData[angles].x = r1nmr[i][0];
        gpu->pbNMRAngleR1R2Slp->_pSysData[angles].y = r2nmr[i][0];
        gpu->pbNMRAngleR3R4Slp->_pSysData[angles].x = r3nmr[i][0];
        gpu->pbNMRAngleR3R4Slp->_pSysData[angles].y = r4nmr[i][0];
        gpu->pbNMRAngleK2K3Slp->_pSysData[angles].x = rk2nmr[i][0];
        gpu->pbNMRAngleK2K3Slp->_pSysData[angles].y = rk3nmr[i][0];
        gpu->pbNMRAngleR1R2Int->_pSysData[angles].x = r1nmr[i][1];
        gpu->pbNMRAngleR1R2Int->_pSysData[angles].y = r2nmr[i][1];
        gpu->pbNMRAngleR3R4Int->_pSysData[angles].x = r3nmr[i][1];
        gpu->pbNMRAngleR3R4Int->_pSysData[angles].y = r4nmr[i][1];
        gpu->pbNMRAngleK2K3Int->_pSysData[angles].x = rk2nmr[i][1];
        gpu->pbNMRAngleK2K3Int->_pSysData[angles].y = rk3nmr[i][1];
        angles++;
      }
      else {
        if (nmrat[i][0] >= 0) { //check if only one atom for first atom index
          gpu->pbNMRCOMAngleID1->_pSysData[COMangles].x = nmrat[i][0] / 3;
        }
        else {
          // This code is for COM restraints. By keeping it a multiple of 3 the number of
          // atoms in the group of atoms can be indexed in the same manner as the CPU code.
          // nmrcom stores the atom ranges ivolved in COM restraints.
          gpu->pbNMRCOMAngleID1->_pSysData[COMangles].x = nmrat[i][0];
          gpu->pbNMRCOMAngleCOMGrp->_pSysData[COMangles * 3].x = -nmrat[i][0]-1;
          gpu->pbNMRCOMAngleCOMGrp->_pSysData[COMangles * 3].y = -nmrat[i][8];
          for (int ip1 = -nmrat[i][0]-1 ; ip1 < -nmrat[i][8] ; ip1++) {
            gpu->pbNMRCOMAngleCOM->_pSysData[ip1].x = nmrcom[ip1][0] - 1;
            gpu->pbNMRCOMAngleCOM->_pSysData[ip1].y = nmrcom[ip1][1] - 1;
          }
        }
        if (nmrat[i][1] >= 0) { //check if only one atom for second atom index
          gpu->pbNMRCOMAngleID1->_pSysData[COMangles].y = nmrat[i][1] / 3;
        }
        else {
          // This code is for COM restraints. By keeping it a multiple of 3 the number of
          // atoms in the group of atoms can be indexed in the same manner as the CPU code.
          // nmrcom stores the atom ranges ivolved in COM restraints.
          gpu->pbNMRCOMAngleID1->_pSysData[COMangles].y = nmrat[i][1];
          gpu->pbNMRCOMAngleCOMGrp->_pSysData[COMangles * 3 + 1].x = -nmrat[i][1]-1;
          gpu->pbNMRCOMAngleCOMGrp->_pSysData[COMangles * 3 + 1].y = -nmrat[i][9];
          for (int ip2 = -nmrat[i][1]-1; ip2 < -nmrat[i][9]; ip2++) {
            gpu->pbNMRCOMAngleCOM->_pSysData[ip2].x = nmrcom[ip2][0] - 1;
            gpu->pbNMRCOMAngleCOM->_pSysData[ip2].y = nmrcom[ip2][1] - 1;
          }
        }
        if (nmrat[i][2] >= 0) { // check if only one atom in 3rd atom index
          gpu->pbNMRCOMAngleID2->_pSysData[COMangles] = nmrat[i][2] / 3;
        }
        else {
          // This code is for COM restraints. By keeping it a multiple of 3 the number of
          // atoms in the group of atoms can be indexed in the same manner as the CPU code.
          // nmrcom stores the atom ranges ivolved in COM group
          gpu->pbNMRCOMAngleID2->_pSysData[COMangles] = nmrat[i][2];
          gpu->pbNMRCOMAngleCOMGrp->_pSysData[COMangles * 3 + 2].x = -nmrat[i][2]-1;
          gpu->pbNMRCOMAngleCOMGrp->_pSysData[COMangles * 3 + 2].y = -nmrat[i][10];
          for (int ip3 = -nmrat[i][2]-1; ip3 < -nmrat[i][10]; ip3++) {
            gpu->pbNMRCOMAngleCOM->_pSysData[ip3].x = nmrcom[ip3][0] - 1;
            gpu->pbNMRCOMAngleCOM->_pSysData[ip3].y = nmrcom[ip3][1] - 1;
          }
        }
        gpu->pbNMRCOMAngleR1R2->_pSysData[COMangles].x = r1nmr[i][1];
        gpu->pbNMRCOMAngleR1R2->_pSysData[COMangles].y = r2nmr[i][1];
        gpu->pbNMRCOMAngleR3R4->_pSysData[COMangles].x = r3nmr[i][1];
        gpu->pbNMRCOMAngleR3R4->_pSysData[COMangles].y = r4nmr[i][1];
        gpu->pbNMRCOMAngleK2K3->_pSysData[COMangles].x = rk2nmr[i][1];
        gpu->pbNMRCOMAngleK2K3->_pSysData[COMangles].y = rk3nmr[i][1];
        gpu->pbNMRCOMAngleAve->_pSysData[COMangles] = (PMEDouble)0.0;
        gpu->pbNMRCOMAngleStep->_pSysData[COMangles].x = nmrst[i][0];
        gpu->pbNMRCOMAngleStep->_pSysData[COMangles].y = nmrst[i][1];
        gpu->pbNMRCOMAngleInc->_pSysData[COMangles] = nmrst[i][2];
        gpu->pbNMRCOMAngleR1R2Slp->_pSysData[COMangles].x = r1nmr[i][0];
        gpu->pbNMRCOMAngleR1R2Slp->_pSysData[COMangles].y = r2nmr[i][0];
        gpu->pbNMRCOMAngleR3R4Slp->_pSysData[COMangles].x = r3nmr[i][0];
        gpu->pbNMRCOMAngleR3R4Slp->_pSysData[COMangles].y = r4nmr[i][0];
        gpu->pbNMRCOMAngleK2K3Slp->_pSysData[COMangles].x = rk2nmr[i][0];
        gpu->pbNMRCOMAngleK2K3Slp->_pSysData[COMangles].y = rk3nmr[i][0];
        gpu->pbNMRCOMAngleR1R2Int->_pSysData[COMangles].x = r1nmr[i][1];
        gpu->pbNMRCOMAngleR1R2Int->_pSysData[COMangles].y = r2nmr[i][1];
        gpu->pbNMRCOMAngleR3R4Int->_pSysData[COMangles].x = r3nmr[i][1];
        gpu->pbNMRCOMAngleR3R4Int->_pSysData[COMangles].y = r4nmr[i][1];
        gpu->pbNMRCOMAngleK2K3Int->_pSysData[COMangles].x = rk2nmr[i][1];
        gpu->pbNMRCOMAngleK2K3Int->_pSysData[COMangles].y = rk3nmr[i][1];
        COMangles++;
      }
    }
    else if (resttype[i] == 3) {
      if (nmrat[i][0] >= 0 && nmrat[i][1] >= 0 && nmrat[i][2] >= 0 && nmrat[i][3] >= 0) { //reg Torsions
        gpu->pbNMRTorsionID1->_pSysData[torsions].x = nmrat[i][0] / 3;
        gpu->pbNMRTorsionID1->_pSysData[torsions].y = nmrat[i][1] / 3;
        gpu->pbNMRTorsionID1->_pSysData[torsions].z = nmrat[i][2] / 3;
        gpu->pbNMRTorsionID1->_pSysData[torsions].w = nmrat[i][3] / 3;
        gpu->pbNMRTorsionR1R2->_pSysData[torsions].x = r1nmr[i][1];
        gpu->pbNMRTorsionR1R2->_pSysData[torsions].y = r2nmr[i][1];
        gpu->pbNMRTorsionR3R4->_pSysData[torsions].x = r3nmr[i][1];
        gpu->pbNMRTorsionR3R4->_pSysData[torsions].y = r4nmr[i][1];
        gpu->pbNMRTorsionK2K3->_pSysData[torsions].x = rk2nmr[i][1];
        gpu->pbNMRTorsionK2K3->_pSysData[torsions].y = rk3nmr[i][1];
        gpu->pbNMRTorsionAve1->_pSysData[torsions] = (PMEDouble)0.0;
        gpu->pbNMRTorsionAve2->_pSysData[torsions] = (PMEDouble)0.0;
        gpu->pbNMRTorsionStep->_pSysData[torsions].x = nmrst[i][0];
        gpu->pbNMRTorsionStep->_pSysData[torsions].y = nmrst[i][1];
        gpu->pbNMRTorsionInc->_pSysData[torsions] = nmrst[i][2];
        gpu->pbNMRTorsionR1R2Slp->_pSysData[torsions].x = r1nmr[i][0];
        gpu->pbNMRTorsionR1R2Slp->_pSysData[torsions].y = r2nmr[i][0];
        gpu->pbNMRTorsionR3R4Slp->_pSysData[torsions].x = r3nmr[i][0];
        gpu->pbNMRTorsionR3R4Slp->_pSysData[torsions].y = r4nmr[i][0];
        gpu->pbNMRTorsionK2K3Slp->_pSysData[torsions].x = rk2nmr[i][0];
        gpu->pbNMRTorsionK2K3Slp->_pSysData[torsions].y = rk3nmr[i][0];
        gpu->pbNMRTorsionR1R2Int->_pSysData[torsions].x = r1nmr[i][1];
        gpu->pbNMRTorsionR1R2Int->_pSysData[torsions].y = r2nmr[i][1];
        gpu->pbNMRTorsionR3R4Int->_pSysData[torsions].x = r3nmr[i][1];
        gpu->pbNMRTorsionR3R4Int->_pSysData[torsions].y = r4nmr[i][1];
        gpu->pbNMRTorsionK2K3Int->_pSysData[torsions].x = rk2nmr[i][1];
        gpu->pbNMRTorsionK2K3Int->_pSysData[torsions].y = rk3nmr[i][1];
        torsions++;
      }
      else { //COM torsions
        if (nmrat[i][0] >= 0) { //check if only one atom for first atom index
          gpu->pbNMRCOMTorsionID1->_pSysData[COMtorsions].x = nmrat[i][0] / 3;
        }
        else {
          // This code is for COM restraints. By keeping it a multiple of 3 the number of
          // atoms in the group of atoms can be indexed in the same manner as the CPU code.
          // nmrcom stores the atom ranges ivolved in COM restraints.
          gpu->pbNMRCOMTorsionID1->_pSysData[COMtorsions].x = nmrat[i][0];
          gpu->pbNMRCOMTorsionCOMGrp->_pSysData[COMtorsions * 4].x = -nmrat[i][0]-1;
          gpu->pbNMRCOMTorsionCOMGrp->_pSysData[COMtorsions * 4].y = -nmrat[i][8];
          for (int ip1 = -nmrat[i][0]-1 ; ip1 < -nmrat[i][8] ; ip1++) {
            gpu->pbNMRCOMTorsionCOM->_pSysData[ip1].x = nmrcom[ip1][0] - 1;
            gpu->pbNMRCOMTorsionCOM->_pSysData[ip1].y = nmrcom[ip1][1] - 1;
          }
        }
        if (nmrat[i][1] >= 0) { //check if only one atom for second atom index
          gpu->pbNMRCOMTorsionID1->_pSysData[COMtorsions].y = nmrat[i][1] / 3;
        }
        else {
          // This code is for COM restraints. By keeping it a multiple of 3 the number of
          // atoms in the group of atoms can be indexed in the same manner as the CPU code.
          // nmrcom stores the atom ranges ivolved in COM restraints.
          gpu->pbNMRCOMTorsionID1->_pSysData[COMtorsions].y = nmrat[i][1];
          gpu->pbNMRCOMTorsionCOMGrp->_pSysData[COMtorsions * 4 + 1].x = -nmrat[i][1]-1;
          gpu->pbNMRCOMTorsionCOMGrp->_pSysData[COMtorsions * 4 + 1].y = -nmrat[i][9];
          for (int ip2 = -nmrat[i][1]-1; ip2 < -nmrat[i][9]; ip2++) {
            gpu->pbNMRCOMTorsionCOM->_pSysData[ip2].x = nmrcom[ip2][0] - 1;
            gpu->pbNMRCOMTorsionCOM->_pSysData[ip2].y = nmrcom[ip2][1] - 1;
          }
        }
        if (nmrat[i][2] >= 0) {
          gpu->pbNMRCOMTorsionID1->_pSysData[COMtorsions].z = nmrat[i][2] / 3;
        }
        else {
          // This code is for COM restraints. By keeping it a multiple of 3 the number of
          // atoms in the group of atoms can be indexed in the same manner as the CPU code.
          // nmrcom stores the atom ranges ivolved in COM restraints.
          gpu->pbNMRCOMTorsionID1->_pSysData[COMtorsions].z = nmrat[i][2];
          gpu->pbNMRCOMTorsionCOMGrp->_pSysData[COMtorsions * 4 + 2].x = -nmrat[i][2]-1;
          gpu->pbNMRCOMTorsionCOMGrp->_pSysData[COMtorsions * 4 + 2].y = -nmrat[i][10];
          for (int ip3 = -nmrat[i][2]-1; ip3 < -nmrat[i][10]; ip3++) {
            gpu->pbNMRCOMTorsionCOM->_pSysData[ip3].x = nmrcom[ip3][0] - 1;
            gpu->pbNMRCOMTorsionCOM->_pSysData[ip3].y = nmrcom[ip3][1] - 1;
          }
        }
        if (nmrat[i][3] >= 0) { //check if only one atom for first atom index
          gpu->pbNMRCOMTorsionID1->_pSysData[COMtorsions].w = nmrat[i][3] / 3;
        }
        else {
          // This code is for COM restraints. By keeping it a multiple of 3 the number of
          // atoms in the group of atoms can be indexed in the same manner as the CPU code.
          // nmrcom stores the atom ranges ivolved in COM restraints.
          gpu->pbNMRCOMTorsionID1->_pSysData[COMtorsions].w = nmrat[i][3];
          gpu->pbNMRCOMTorsionCOMGrp->_pSysData[COMtorsions * 4 + 3].x = -nmrat[i][3]-1;
          gpu->pbNMRCOMTorsionCOMGrp->_pSysData[COMtorsions * 4 + 3].y = -nmrat[i][11];
          for (int ip4 = -nmrat[i][3]-1; ip4 < -nmrat[i][11]; ip4++) {
            gpu->pbNMRCOMTorsionCOM->_pSysData[ip4].x = nmrcom[ip4][0] - 1;
            gpu->pbNMRCOMTorsionCOM->_pSysData[ip4].y = nmrcom[ip4][1] - 1;
          }
        }
        gpu->pbNMRCOMTorsionR1R2->_pSysData[COMtorsions].x = r1nmr[i][1];
        gpu->pbNMRCOMTorsionR1R2->_pSysData[COMtorsions].y = r2nmr[i][1];
        gpu->pbNMRCOMTorsionR3R4->_pSysData[COMtorsions].x = r3nmr[i][1];
        gpu->pbNMRCOMTorsionR3R4->_pSysData[COMtorsions].y = r4nmr[i][1];
        gpu->pbNMRCOMTorsionK2K3->_pSysData[COMtorsions].x = rk2nmr[i][1];
        gpu->pbNMRCOMTorsionK2K3->_pSysData[COMtorsions].y = rk3nmr[i][1];
        gpu->pbNMRCOMTorsionAve->_pSysData[COMtorsions] = (PMEDouble)0.0;
        gpu->pbNMRCOMTorsionStep->_pSysData[COMtorsions].x = nmrst[i][0];
        gpu->pbNMRCOMTorsionStep->_pSysData[COMtorsions].y = nmrst[i][1];
        gpu->pbNMRCOMTorsionInc->_pSysData[COMtorsions] = nmrst[i][2];
        gpu->pbNMRCOMTorsionR1R2Slp->_pSysData[COMtorsions].x = r1nmr[i][0];
        gpu->pbNMRCOMTorsionR1R2Slp->_pSysData[COMtorsions].y = r2nmr[i][0];
        gpu->pbNMRCOMTorsionR3R4Slp->_pSysData[COMtorsions].x = r3nmr[i][0];
        gpu->pbNMRCOMTorsionR3R4Slp->_pSysData[COMtorsions].y = r4nmr[i][0];
        gpu->pbNMRCOMTorsionK2K3Slp->_pSysData[COMtorsions].x = rk2nmr[i][0];
        gpu->pbNMRCOMTorsionK2K3Slp->_pSysData[COMtorsions].y = rk3nmr[i][0];
        gpu->pbNMRCOMTorsionR1R2Int->_pSysData[COMtorsions].x = r1nmr[i][1];
        gpu->pbNMRCOMTorsionR1R2Int->_pSysData[COMtorsions].y = r2nmr[i][1];
        gpu->pbNMRCOMTorsionR3R4Int->_pSysData[COMtorsions].x = r3nmr[i][1];
        gpu->pbNMRCOMTorsionR3R4Int->_pSysData[COMtorsions].y = r4nmr[i][1];
        gpu->pbNMRCOMTorsionK2K3Int->_pSysData[COMtorsions].x = rk2nmr[i][1];
        gpu->pbNMRCOMTorsionK2K3Int->_pSysData[COMtorsions].y = rk3nmr[i][1];
        COMtorsions++;
      }
    }
  } // end of nmrnum loop

  // Upload data
  gpu->pbNMRJarData->Upload();
  gpu->pbNMRDistanceID->Upload();
  gpu->pbNMRDistanceR1R2->Upload();
  gpu->pbNMRDistanceR3R4->Upload();
  gpu->pbNMRDistanceK2K3->Upload();
  gpu->pbNMRDistanceK4->Upload();
  gpu->pbNMRDistanceAve->Upload();
  gpu->pbNMRDistanceStep->Upload();
  gpu->pbNMRDistanceInc->Upload();
  gpu->pbNMRDistanceR1R2Slp->Upload();
  gpu->pbNMRDistanceR3R4Slp->Upload();
  gpu->pbNMRDistanceK2K3Slp->Upload();
  gpu->pbNMRDistanceK4Slp->Upload();
  gpu->pbNMRDistanceR1R2Int->Upload();
  gpu->pbNMRDistanceR3R4Int->Upload();
  gpu->pbNMRDistanceK2K3Int->Upload();
  gpu->pbNMRDistanceK4Int->Upload();

  gpu->pbNMRCOMDistanceID->Upload();
  gpu->pbNMRCOMDistanceCOM->Upload();
  gpu->pbNMRCOMDistanceCOMGrp->Upload();
  gpu->pbNMRCOMDistanceR1R2->Upload();
  gpu->pbNMRCOMDistanceR3R4->Upload();
  gpu->pbNMRCOMDistanceK2K3->Upload();
  gpu->pbNMRCOMDistanceK4->Upload();
  gpu->pbNMRCOMDistanceAve->Upload();
  gpu->pbNMRCOMDistanceStep->Upload();
  gpu->pbNMRCOMDistanceInc->Upload();
  gpu->pbNMRCOMDistanceR1R2Slp->Upload();
  gpu->pbNMRCOMDistanceR3R4Slp->Upload();
  gpu->pbNMRCOMDistanceK2K3Slp->Upload();
  gpu->pbNMRCOMDistanceK4Slp->Upload();
  gpu->pbNMRCOMDistanceR1R2Int->Upload();
  gpu->pbNMRCOMDistanceR3R4Int->Upload();
  gpu->pbNMRCOMDistanceK2K3Int->Upload();
  gpu->pbNMRCOMDistanceK4Int->Upload();
  gpu->pbNMRCOMDistanceWeights->Upload();
  gpu->pbNMRCOMDistanceXYZ->Upload();

  gpu->pbNMRr6avDistanceID->Upload();
  gpu->pbNMRr6avDistancer6av->Upload();
  gpu->pbNMRr6avDistancer6avGrp->Upload();
  gpu->pbNMRr6avDistanceR1R2->Upload();
  gpu->pbNMRr6avDistanceR3R4->Upload();
  gpu->pbNMRr6avDistanceK2K3->Upload();
  gpu->pbNMRr6avDistanceK4->Upload();
  gpu->pbNMRr6avDistanceAve->Upload();
  gpu->pbNMRr6avDistanceStep->Upload();
  gpu->pbNMRr6avDistanceInc->Upload();
  gpu->pbNMRr6avDistanceR1R2Slp->Upload();
  gpu->pbNMRr6avDistanceR3R4Slp->Upload();
  gpu->pbNMRr6avDistanceK2K3Slp->Upload();
  gpu->pbNMRr6avDistanceK4Slp->Upload();
  gpu->pbNMRr6avDistanceR1R2Int->Upload();
  gpu->pbNMRr6avDistanceR3R4Int->Upload();
  gpu->pbNMRr6avDistanceK2K3Int->Upload();
  gpu->pbNMRr6avDistanceK4Int->Upload();
  gpu->pbNMRr6avDistanceWeights->Upload();
  gpu->pbNMRr6avDistanceXYZ->Upload();

  gpu->pbNMRAngleID1->Upload();
  gpu->pbNMRAngleID2->Upload();
  gpu->pbNMRAngleR1R2->Upload();
  gpu->pbNMRAngleR3R4->Upload();
  gpu->pbNMRAngleK2K3->Upload();
  gpu->pbNMRAngleK4->Upload();
  gpu->pbNMRAngleAve->Upload();
  gpu->pbNMRAngleStep->Upload();
  gpu->pbNMRAngleInc->Upload();
  gpu->pbNMRAngleR1R2Slp->Upload();
  gpu->pbNMRAngleR3R4Slp->Upload();
  gpu->pbNMRAngleK2K3Slp->Upload();
  gpu->pbNMRAngleK4Slp->Upload();
  gpu->pbNMRAngleR1R2Int->Upload();
  gpu->pbNMRAngleR3R4Int->Upload();
  gpu->pbNMRAngleK2K3Int->Upload();
  gpu->pbNMRAngleK4Int->Upload();

  gpu->pbNMRCOMAngleID1->Upload();
  gpu->pbNMRCOMAngleID2->Upload();
  gpu->pbNMRCOMAngleCOM->Upload();
  gpu->pbNMRCOMAngleCOMGrp->Upload();
  gpu->pbNMRCOMAngleR1R2->Upload();
  gpu->pbNMRCOMAngleR3R4->Upload();
  gpu->pbNMRCOMAngleK2K3->Upload();
  gpu->pbNMRCOMAngleK4->Upload();
  gpu->pbNMRCOMAngleAve->Upload();
  gpu->pbNMRCOMAngleStep->Upload();
  gpu->pbNMRCOMAngleInc->Upload();
  gpu->pbNMRCOMAngleR1R2Slp->Upload();
  gpu->pbNMRCOMAngleR3R4Slp->Upload();
  gpu->pbNMRCOMAngleK2K3Slp->Upload();
  gpu->pbNMRCOMAngleK4Slp->Upload();
  gpu->pbNMRCOMAngleR1R2Int->Upload();
  gpu->pbNMRCOMAngleR3R4Int->Upload();
  gpu->pbNMRCOMAngleK2K3Int->Upload();
  gpu->pbNMRCOMAngleK4Int->Upload();

  gpu->pbNMRTorsionID1->Upload();
  gpu->pbNMRTorsionR1R2->Upload();
  gpu->pbNMRTorsionR3R4->Upload();
  gpu->pbNMRTorsionK2K3->Upload();
  gpu->pbNMRTorsionK4->Upload();
  gpu->pbNMRTorsionAve1->Upload();
  gpu->pbNMRTorsionAve2->Upload();
  gpu->pbNMRTorsionStep->Upload();
  gpu->pbNMRTorsionInc->Upload();
  gpu->pbNMRTorsionR1R2Slp->Upload();
  gpu->pbNMRTorsionR3R4Slp->Upload();
  gpu->pbNMRTorsionK2K3Slp->Upload();
  gpu->pbNMRTorsionK4Slp->Upload();
  gpu->pbNMRTorsionR1R2Int->Upload();
  gpu->pbNMRTorsionR3R4Int->Upload();
  gpu->pbNMRTorsionK2K3Int->Upload();
  gpu->pbNMRTorsionK4Int->Upload();

  gpu->pbNMRCOMTorsionID1->Upload();
  gpu->pbNMRCOMTorsionCOM->Upload();
  gpu->pbNMRCOMTorsionCOMGrp->Upload();
  gpu->pbNMRCOMTorsionR1R2->Upload();
  gpu->pbNMRCOMTorsionR3R4->Upload();
  gpu->pbNMRCOMTorsionK2K3->Upload();
  gpu->pbNMRCOMTorsionK4->Upload();
  gpu->pbNMRCOMTorsionAve->Upload();
  gpu->pbNMRCOMTorsionStep->Upload();
  gpu->pbNMRCOMTorsionInc->Upload();
  gpu->pbNMRCOMTorsionR1R2Slp->Upload();
  gpu->pbNMRCOMTorsionR3R4Slp->Upload();
  gpu->pbNMRCOMTorsionK2K3Slp->Upload();
  gpu->pbNMRCOMTorsionK4Slp->Upload();
  gpu->pbNMRCOMTorsionR1R2Int->Upload();
  gpu->pbNMRCOMTorsionR3R4Int->Upload();
  gpu->pbNMRCOMTorsionK2K3Int->Upload();
  gpu->pbNMRCOMTorsionK4Int->Upload();

  // Set constant bank data
  gpu->sim.pNMRJarData             = gpu->pbNMRJarData->_pDevData;
  gpu->sim.pNMRDistanceID          = gpu->pbNMRDistanceID->_pDevData;
  gpu->sim.pNMRDistanceR1R2        = gpu->pbNMRDistanceR1R2->_pDevData;
  gpu->sim.pNMRDistanceR3R4        = gpu->pbNMRDistanceR3R4->_pDevData;
  gpu->sim.pNMRDistanceK2K3        = gpu->pbNMRDistanceK2K3->_pDevData;
  gpu->sim.pNMRDistanceAve         = gpu->pbNMRDistanceAve->_pDevData;
  gpu->sim.pNMRDistanceTgtVal      = gpu->pbNMRDistanceTgtVal->_pDevData;
  gpu->sim.pNMRDistanceStep        = gpu->pbNMRDistanceStep->_pDevData;
  gpu->sim.pNMRDistanceInc         = gpu->pbNMRDistanceInc->_pDevData;
  gpu->sim.pNMRDistanceR1R2Slp     = gpu->pbNMRDistanceR1R2Slp->_pDevData;
  gpu->sim.pNMRDistanceR3R4Slp     = gpu->pbNMRDistanceR3R4Slp->_pDevData;
  gpu->sim.pNMRDistanceK2K3Slp     = gpu->pbNMRDistanceK2K3Slp->_pDevData;
  gpu->sim.pNMRDistanceR1R2Int     = gpu->pbNMRDistanceR1R2Int->_pDevData;
  gpu->sim.pNMRDistanceR3R4Int     = gpu->pbNMRDistanceR3R4Int->_pDevData;
  gpu->sim.pNMRDistanceK2K3Int     = gpu->pbNMRDistanceK2K3Int->_pDevData;

  gpu->sim.pNMRCOMDistanceID       = gpu->pbNMRCOMDistanceID->_pDevData;
  gpu->sim.pNMRCOMDistanceCOM      = gpu->pbNMRCOMDistanceCOM->_pDevData;
  gpu->sim.pNMRCOMDistanceCOMGrp   = gpu->pbNMRCOMDistanceCOMGrp->_pDevData;
  gpu->sim.pNMRCOMDistanceR1R2     = gpu->pbNMRCOMDistanceR1R2->_pDevData;
  gpu->sim.pNMRCOMDistanceR3R4     = gpu->pbNMRCOMDistanceR3R4->_pDevData;
  gpu->sim.pNMRCOMDistanceK2K3     = gpu->pbNMRCOMDistanceK2K3->_pDevData;
  gpu->sim.pNMRCOMDistanceAve      = gpu->pbNMRCOMDistanceAve->_pDevData;
  gpu->sim.pNMRCOMDistanceTgtVal   = gpu->pbNMRCOMDistanceTgtVal->_pDevData;
  gpu->sim.pNMRCOMDistanceStep     = gpu->pbNMRCOMDistanceStep->_pDevData;
  gpu->sim.pNMRCOMDistanceInc      = gpu->pbNMRCOMDistanceInc->_pDevData;
  gpu->sim.pNMRCOMDistanceR1R2Slp  = gpu->pbNMRCOMDistanceR1R2Slp->_pDevData;
  gpu->sim.pNMRCOMDistanceR3R4Slp  = gpu->pbNMRCOMDistanceR3R4Slp->_pDevData;
  gpu->sim.pNMRCOMDistanceK2K3Slp  = gpu->pbNMRCOMDistanceK2K3Slp->_pDevData;
  gpu->sim.pNMRCOMDistanceR1R2Int  = gpu->pbNMRCOMDistanceR1R2Int->_pDevData;
  gpu->sim.pNMRCOMDistanceR3R4Int  = gpu->pbNMRCOMDistanceR3R4Int->_pDevData;
  gpu->sim.pNMRCOMDistanceK2K3Int  = gpu->pbNMRCOMDistanceK2K3Int->_pDevData;
  gpu->sim.pNMRCOMDistanceWeights  = gpu->pbNMRCOMDistanceWeights->_pDevData;
  gpu->sim.pNMRCOMDistanceXYZ      = gpu->pbNMRCOMDistanceXYZ->_pDevData;

  gpu->sim.pNMRr6avDistanceID      = gpu->pbNMRr6avDistanceID->_pDevData;
  gpu->sim.pNMRr6avDistancer6av    = gpu->pbNMRr6avDistancer6av->_pDevData;
  gpu->sim.pNMRr6avDistancer6avGrp = gpu->pbNMRr6avDistancer6avGrp->_pDevData;
  gpu->sim.pNMRr6avDistanceR1R2    = gpu->pbNMRr6avDistanceR1R2->_pDevData;
  gpu->sim.pNMRr6avDistanceR3R4    = gpu->pbNMRr6avDistanceR3R4->_pDevData;
  gpu->sim.pNMRr6avDistanceK2K3    = gpu->pbNMRr6avDistanceK2K3->_pDevData;
  gpu->sim.pNMRr6avDistanceAve     = gpu->pbNMRr6avDistanceAve->_pDevData;
  gpu->sim.pNMRr6avDistanceTgtVal  = gpu->pbNMRr6avDistanceTgtVal->_pDevData;
  gpu->sim.pNMRr6avDistanceStep    = gpu->pbNMRr6avDistanceStep->_pDevData;
  gpu->sim.pNMRr6avDistanceInc     = gpu->pbNMRr6avDistanceInc->_pDevData;
  gpu->sim.pNMRr6avDistanceR1R2Slp = gpu->pbNMRr6avDistanceR1R2Slp->_pDevData;
  gpu->sim.pNMRr6avDistanceR3R4Slp = gpu->pbNMRr6avDistanceR3R4Slp->_pDevData;
  gpu->sim.pNMRr6avDistanceK2K3Slp = gpu->pbNMRr6avDistanceK2K3Slp->_pDevData;
  gpu->sim.pNMRr6avDistanceR1R2Int = gpu->pbNMRr6avDistanceR1R2Int->_pDevData;
  gpu->sim.pNMRr6avDistanceR3R4Int = gpu->pbNMRr6avDistanceR3R4Int->_pDevData;
  gpu->sim.pNMRr6avDistanceK2K3Int = gpu->pbNMRr6avDistanceK2K3Int->_pDevData;
  gpu->sim.pNMRr6avDistanceWeights = gpu->pbNMRr6avDistanceWeights->_pDevData;
  gpu->sim.pNMRr6avDistanceXYZ     = gpu->pbNMRr6avDistanceXYZ->_pDevData;

  gpu->sim.pNMRAngleID1            = gpu->pbNMRAngleID1->_pDevData;
  gpu->sim.pNMRAngleID2            = gpu->pbNMRAngleID2->_pDevData;
  gpu->sim.pNMRAngleR1R2           = gpu->pbNMRAngleR1R2->_pDevData;
  gpu->sim.pNMRAngleR3R4           = gpu->pbNMRAngleR3R4->_pDevData;
  gpu->sim.pNMRAngleK2K3           = gpu->pbNMRAngleK2K3->_pDevData;
  gpu->sim.pNMRAngleAve            = gpu->pbNMRAngleAve->_pDevData;
  gpu->sim.pNMRAngleTgtVal         = gpu->pbNMRAngleTgtVal->_pDevData;
  gpu->sim.pNMRAngleStep           = gpu->pbNMRAngleStep->_pDevData;
  gpu->sim.pNMRAngleInc            = gpu->pbNMRAngleInc->_pDevData;
  gpu->sim.pNMRAngleR1R2Slp        = gpu->pbNMRAngleR1R2Slp->_pDevData;
  gpu->sim.pNMRAngleR3R4Slp        = gpu->pbNMRAngleR3R4Slp->_pDevData;
  gpu->sim.pNMRAngleK2K3Slp        = gpu->pbNMRAngleK2K3Slp->_pDevData;
  gpu->sim.pNMRAngleR1R2Int        = gpu->pbNMRAngleR1R2Int->_pDevData;
  gpu->sim.pNMRAngleR3R4Int        = gpu->pbNMRAngleR3R4Int->_pDevData;
  gpu->sim.pNMRAngleK2K3Int        = gpu->pbNMRAngleK2K3Int->_pDevData;

  gpu->sim.pNMRCOMAngleID1       = gpu->pbNMRCOMAngleID1->_pDevData;
  gpu->sim.pNMRCOMAngleID2       = gpu->pbNMRCOMAngleID2->_pDevData;
  gpu->sim.pNMRCOMAngleCOM      = gpu->pbNMRCOMAngleCOM->_pDevData;
  gpu->sim.pNMRCOMAngleCOMGrp   = gpu->pbNMRCOMAngleCOMGrp->_pDevData;
  gpu->sim.pNMRCOMAngleR1R2     = gpu->pbNMRCOMAngleR1R2->_pDevData;
  gpu->sim.pNMRCOMAngleR3R4     = gpu->pbNMRCOMAngleR3R4->_pDevData;
  gpu->sim.pNMRCOMAngleK2K3     = gpu->pbNMRCOMAngleK2K3->_pDevData;
  gpu->sim.pNMRCOMAngleAve      = gpu->pbNMRCOMAngleAve->_pDevData;
  gpu->sim.pNMRCOMAngleTgtVal   = gpu->pbNMRCOMAngleTgtVal->_pDevData;
  gpu->sim.pNMRCOMAngleStep     = gpu->pbNMRCOMAngleStep->_pDevData;
  gpu->sim.pNMRCOMAngleInc      = gpu->pbNMRCOMAngleInc->_pDevData;
  gpu->sim.pNMRCOMAngleR1R2Slp  = gpu->pbNMRCOMAngleR1R2Slp->_pDevData;
  gpu->sim.pNMRCOMAngleR3R4Slp  = gpu->pbNMRCOMAngleR3R4Slp->_pDevData;
  gpu->sim.pNMRCOMAngleK2K3Slp  = gpu->pbNMRCOMAngleK2K3Slp->_pDevData;
  gpu->sim.pNMRCOMAngleR1R2Int  = gpu->pbNMRCOMAngleR1R2Int->_pDevData;
  gpu->sim.pNMRCOMAngleR3R4Int  = gpu->pbNMRCOMAngleR3R4Int->_pDevData;
  gpu->sim.pNMRCOMAngleK2K3Int  = gpu->pbNMRCOMAngleK2K3Int->_pDevData;
  gpu->sim.pNMRTorsionID1          = gpu->pbNMRTorsionID1->_pDevData;
  gpu->sim.pNMRTorsionR1R2         = gpu->pbNMRTorsionR1R2->_pDevData;
  gpu->sim.pNMRTorsionR3R4         = gpu->pbNMRTorsionR3R4->_pDevData;
  gpu->sim.pNMRTorsionK2K3         = gpu->pbNMRTorsionK2K3->_pDevData;
  gpu->sim.pNMRTorsionAve1         = gpu->pbNMRTorsionAve1->_pDevData;
  gpu->sim.pNMRTorsionAve2         = gpu->pbNMRTorsionAve2->_pDevData;
  gpu->sim.pNMRTorsionTgtVal       = gpu->pbNMRTorsionTgtVal->_pDevData;
  gpu->sim.pNMRTorsionStep         = gpu->pbNMRTorsionStep->_pDevData;
  gpu->sim.pNMRTorsionInc          = gpu->pbNMRTorsionInc->_pDevData;
  gpu->sim.pNMRTorsionR1R2Slp      = gpu->pbNMRTorsionR1R2Slp->_pDevData;
  gpu->sim.pNMRTorsionR3R4Slp      = gpu->pbNMRTorsionR3R4Slp->_pDevData;
  gpu->sim.pNMRTorsionK2K3Slp      = gpu->pbNMRTorsionK2K3Slp->_pDevData;
  gpu->sim.pNMRTorsionR1R2Int      = gpu->pbNMRTorsionR1R2Int->_pDevData;
  gpu->sim.pNMRTorsionR3R4Int      = gpu->pbNMRTorsionR3R4Int->_pDevData;
  gpu->sim.pNMRTorsionK2K3Int      = gpu->pbNMRTorsionK2K3Int->_pDevData;

  gpu->sim.pNMRCOMTorsionID1       = gpu->pbNMRCOMTorsionID1->_pDevData;
  gpu->sim.pNMRCOMTorsionCOM      = gpu->pbNMRCOMTorsionCOM->_pDevData;
  gpu->sim.pNMRCOMTorsionCOMGrp   = gpu->pbNMRCOMTorsionCOMGrp->_pDevData;
  gpu->sim.pNMRCOMTorsionR1R2     = gpu->pbNMRCOMTorsionR1R2->_pDevData;
  gpu->sim.pNMRCOMTorsionR3R4     = gpu->pbNMRCOMTorsionR3R4->_pDevData;
  gpu->sim.pNMRCOMTorsionK2K3     = gpu->pbNMRCOMTorsionK2K3->_pDevData;
  gpu->sim.pNMRCOMTorsionAve      = gpu->pbNMRCOMTorsionAve->_pDevData;
  gpu->sim.pNMRCOMTorsionTgtVal   = gpu->pbNMRCOMTorsionTgtVal->_pDevData;
  gpu->sim.pNMRCOMTorsionStep     = gpu->pbNMRCOMTorsionStep->_pDevData;
  gpu->sim.pNMRCOMTorsionInc      = gpu->pbNMRCOMTorsionInc->_pDevData;
  gpu->sim.pNMRCOMTorsionR1R2Slp  = gpu->pbNMRCOMTorsionR1R2Slp->_pDevData;
  gpu->sim.pNMRCOMTorsionR3R4Slp  = gpu->pbNMRCOMTorsionR3R4Slp->_pDevData;
  gpu->sim.pNMRCOMTorsionK2K3Slp  = gpu->pbNMRCOMTorsionK2K3Slp->_pDevData;
  gpu->sim.pNMRCOMTorsionR1R2Int  = gpu->pbNMRCOMTorsionR1R2Int->_pDevData;
  gpu->sim.pNMRCOMTorsionR3R4Int  = gpu->pbNMRCOMTorsionR3R4Int->_pDevData;
  gpu->sim.pNMRCOMTorsionK2K3Int  = gpu->pbNMRCOMTorsionK2K3Int->_pDevData;
  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gpu_refresh_charges_: refresh the charges based on things that happened on the CPU.  This
//                       is called mainly by constant pH calculations.
//
// Arguments:
//   cit_nb14:     list of 1-4 non-bonded pairs and their parameter indices (different 1-4
//                 scaling factors can be indexed)
//   gbl_one_scee: the different 1-4 electrostatic scaling factors in use
//   qterm:        atomic partial charges for Amber PME
//---------------------------------------------------------------------------------------------
extern "C" void gpu_refresh_charges_(int cit_nb14[][3], double gbl_one_scee[], double qterm[])
{
  PRINTMETHOD("gpu_refresh_charges");

  // Copy 1-4 interactions
  if (gpu->ntf < 8) {
    int i, j;
    for (i = 0; i < gpu->sim.nb14s; i++) {
      int parm_idx = cit_nb14[i][2] - 1;
      gpu->pbNb141->_pSysData[i].x = (parm_idx>=0)? gbl_one_scee[parm_idx] * qterm[abs(cit_nb14[i][0]) - 1] *
                                     qterm[abs(cit_nb14[i][1]) - 1] : 0.0;
    }

    // Updating the Nb141 array on the device keeps things consistent between the device and
    // host memory, but is now legacy behavior--the real update for 1:4 charge interactions
    // happens in the bond work units below.
    gpu->pbNb141->Upload();
    gpu->pbChargeRefreshBuffer->Upload(qterm);
    kRefreshCharges(gpu);

    // Update bond work units
    PMEFloat *ftmp = gpu->pbBondWorkUnitPFLOAT->_pSysData;
    PMEDouble2 *d2tmp = gpu->pbNb141->_pSysData;
    for (i = 0; i < gpu->sim.bondWorkUnits; i++) {
      int offset = gpu->bondWorkUnitRecord[i].nb14PFloatIdx;
      int *itmp;
      itmp = gpu->bondWorkUnitRecord[i].nb14List;
      for (j = 0; j < gpu->bondWorkUnitRecord[i].nnb14; j++) {
        if (itmp[j] >= 0) {
          ftmp[offset + j] = d2tmp[itmp[j]].x;
        }
      }
      itmp = gpu->bondWorkUnitRecord[i].atomList;
      offset = gpu->bondWorkUnitRecord[i].qPFloatIdx;
      for (j = 0; j < gpu->bondWorkUnitRecord[i].natom; j++) {
        if (itmp[j] >= 0) {
          ftmp[offset + j] = qterm[itmp[j]];
        }
      }
    }
    gpu->pbBondWorkUnitPFLOAT->Upload();
  }
}

//---------------------------------------------------------------------------------------------
// DetectLennardJonesPotentials: scan through the Lennard-Jones interaction matrix to determine
//                               whether a particular atom type will interact with nonzero
//                               potentials in combination with any others.  A bitmask is
//                               constructed to record this information for rapid access.
//                               C4PairwiseCUDA modified to add cn7, cn8 and C4Pairwise
//---------------------------------------------------------------------------------------------
static void DetectLennardJonesPotentials(int ico[], double cn1[], double cn2[], double cn6[],
                                         int cn7[], double cn8[], int *C4Pairwise, int *lj1264, int *plj1264, int *ntypes)
{
  PRINTMETHOD("DetectLennardJonesPotentials");
  int i, j;
  dmat LJmat, r4mat, r4arr; // C4PairwiseCUDA2023

  // Compute a matrix of the energy that two particles of either
  // type would have if separated by a distance of 3A.
  LJmat = CreateDmat(*ntypes, *ntypes);
  r4mat = CreateDmat(*ntypes, *ntypes);
  r4arr = CreateDmat(1,*C4Pairwise * 4); //C4PairwiseCUDA2023 array to keep both cn7 and cn8
  for (i = 1; i <= *ntypes; i++) {
    for (j = 1; j <= i; j++) {
      int pos = ico[*ntypes * (i - 1) + j - 1] - 1;
      if (pos >= 0) {
        LJmat.map[i-1][j-1] = (cn1[pos]/729.0 + cn2[pos])/729.0;
        LJmat.map[j-1][i-1] = LJmat.map[i-1][j-1];
      }
    }
  }
  if (*lj1264 == 1 && *plj1264 == 0) {
    for (i = 1; i <= *ntypes; i++) {
      for (j = 1; j <= i; j++) {
        int pos = ico[*ntypes * (i - 1) + j - 1] - 1;
        if (pos >= 0) {
          r4mat.map[i-1][j-1] = cn6[pos] / 81.0;
        }
      }
    }
  }
  if (*lj1264 == 0 && *plj1264 == 1) {
    for (i = 1; i <= *C4Pairwise * 4; i+=4) {
      r4arr.map[0][i-1] = cn7[i-1];
      r4arr.map[0][i] = cn7[i];
      r4arr.map[0][i+1] = cn7[i+1];
      r4arr.map[0][i+2] = cn8[(i-1) / 4] / 81.0;
    }
  }
  if (*lj1264 == 1 && *plj1264 == 1) {
    for (i = 1; i <= *ntypes; i++) {
      for (j = 1; j <= i; j++) {
        int pos = ico[*ntypes * (i - 1) + j - 1] - 1;
        if (pos >= 0) {
          r4mat.map[i-1][j-1] = cn6[pos] / 81.0;
        }
      }
    }
    for (i = 1; i <= *C4Pairwise; i++) {
      std::cout << "Assigning values to r4arr" << std::endl; //C4PairwiseCUDA2023debug
      r4arr.map[0][i*4-4] = cn7[i*3-3];
      r4arr.map[0][i*4-3] = cn7[i*3-2];
      r4arr.map[0][i*4-2] = cn7[i*3-1];
      r4arr.map[0][i*4-1] = cn8[i] / 81.0;
    }
  }
  
  // Loop over both matrices and detect whether any atom types experience non-zero
  // potentials in the presence of any others.
  unsigned int T1 = 0;
  unsigned int T2 = 0;
  unsigned int T3 = 0;
  unsigned int T4 = 0;
  for (i = 0; i < *ntypes; i++) {
    int ljpo = 0;
    for (j = 0; j < *ntypes; j++) {
      ljpo += (fabs(LJmat.map[i][j]) > 1.0e-8) + (fabs(r4mat.map[i][j]) > 1.0e-8);
    }
    for (j = 0; j < *C4Pairwise; j++) { // C4PairwiseCUDA2023
      //std::cout << "before ljpo" << std::endl; //C4PairwiseCUDA2023debug
      ljpo += (fabs(r4mat.map[0][j]) > 1.0e-8);  //Maybe adding by atom, not atom type is OK
      //std::cout << "after ljpo" << std::endl; //C4PairwiseCUDA2023debug
    }
    if (ljpo > 0) {
      if (i < 32) {
        T1 += 0x1 << i;
      }
      else if (i < 64) {
        T2 += 0x1 << (i - 32);
      }
      else if (i < 96) {
        T3 += 0x1 << (i - 64);
      }
      else if (i < 128) {
        T4 += 0x1 << (i - 96);
      }
    }
  }
  gpu->sim.ljexist[0] = T1;
  gpu->sim.ljexist[1] = T2;
  gpu->sim.ljexist[2] = T3;
  gpu->sim.ljexist[3] = T4;

  // Free allocated memory
  DestroyDmat(&LJmat);
  DestroyDmat(&r4mat);
  DestroyDmat(&r4arr); //C4PairwiseCUDA2023
}

//---------------------------------------------------------------------------------------------
// gpu_nb14_setup_: setup for 1-4 interactions on the GPU
//
// Arguments:
//   nb14_cnt:     the number of 1-4 interactions
//   cit_nb14:     list of 1-4 non-bonded pairs and their parameter indices (different 1-4
//                 scaling factors can be indexed)
//   gbl_one_scee: the different 1-4 electrostatic scaling factors in use
//   gbl_one_scnb: the different 1-4 Lennard-Jones scaling factors in use
//   ntypes:       the number of Lennard-Jones parameter types
//   iac:          oh, no--something from sander nomenclature... it's the Lennard-Jones
//                 parameter type index of each atom
//   ico:          the non-bonded parameter type index of each atom
//   cn1:          Lennard-Jones A parameters
//   cn2:          Lennard-Jones B parameters
//   cn6:          Lennard-Jones C parameters
//   cn7:          Lennard-Jones D indices C4PairwiseCUDA
//   cn8:          Lennard-Jones D parameters C4PairwiseCUDA
//   C4Pairwise:   Count of pairwise LJ interactions
//   cn114:        Lennard-Jones A 1-4 parameters
//   cn114:        Lennard-Jones B 1-4 parameters
//---------------------------------------------------------------------------------------------
extern "C" void gpu_nb14_setup_(int *nb14_cnt, int cit_nb14[][3], double gbl_one_scee[],
                                double gbl_one_scnb[], int *ntypes, int *lj1264, int *plj1264, int iac[],
                                int ico[], double cn1[], double cn2[], double cn6[],
                                int cn7[], double cn8[], int *C4Pairwise, double cn114[], double cn214[])
{
  PRINTMETHOD("gpu_nb14_setup");
  int i, j, k;

  // Allocate/reallocate GPU nb14 buffers
  int nb14s = 0;
  if (gpu->ntf < 8) {
    nb14s = *nb14_cnt;
  }
  else {
    nb14s = 0;
  }
  delete gpu->pbNb141;
  delete gpu->pbNb142;
  delete gpu->pbNb14ID;

#ifdef GVERBOSE
  printf("%d 1-4 nonbonds, %d active\n", *nb14_cnt, nb14s);
#endif
  if (nb14s > 0) {
    gpu->pbNb141 = new GpuBuffer<PMEDouble2>(nb14s);
    gpu->pbNb142 = new GpuBuffer<PMEDouble2>(nb14s);
    gpu->pbNb14ID = new GpuBuffer<int2>(nb14s);
  }

#ifdef GVERBOSE
  printf("%d atom types\n", *ntypes);
#endif

  // Generate sigma and epsilon from atom types
  double* sigma   = new double[*ntypes];
  double* epsilon = new double[*ntypes];
  for (i = 0; i < *ntypes; i++) {
    int nbtype = ico[*ntypes * i + i] - 1;
    if (nbtype >= 0) {
      double c1  = cn1[nbtype];
      double c2  = cn2[nbtype];
      double sig = pow(c1 / c2, 1.0 / 6.0);
      double eps = (c2 * c2) / (4.0 * c1);
      sigma[i]   = sig / 2.0;
      epsilon[i] = sqrt(eps);
    }
    else {
      sigma[i]   = 0.5;
      epsilon[i] = 0.0;
    }
  }
  // Detect Lennard-Jones atom types that have zero potential
  // with themselves and all other atom types.
  DetectLennardJonesPotentials(ico, cn1, cn2, cn6, cn7, cn8, C4Pairwise, lj1264, plj1264, ntypes); // C4PairwiseCUDA
  // Set up and upload Lennard-Jones and NBFIX data
  int LJMultiplier  = 1 + (gpu->sim.ti_mode >= 2) + (*lj1264 == 1);
  gpu->sim.LJTypes  = *ntypes + 1;
  gpu->sim.LJOffset = (*ntypes + 1) * (*ntypes + 1);
  gpu->sim.LJTerms  = gpu->sim.LJOffset * LJMultiplier + *C4Pairwise*4; //C4PairwiseCUDA2023
  gpu->pbLJTerm     = new GpuBuffer<PMEFloat2>(gpu->sim.LJTerms);
  memset(gpu->pbLJTerm->_pSysData, 0, gpu->sim.LJTerms * sizeof(PMEFloat2)); 
  for (i = 1; i <= *ntypes; i++) { 
    for (j = 1; j <= i; j++) {
      int pos = ico[*ntypes * (i - 1) + j - 1] - 1;
      int o1  = (i * gpu->sim.LJTypes) + j;
      int o2  = (j * gpu->sim.LJTypes) + i;

      // pbLJTerm is initialized to zero.  Don't set it to anything if there is no
      // Lennard-Jones interaction for two atom types.
      if (pos >= 0) {
        if (gpu->ips == 0) {
          gpu->pbLJTerm->_pSysData[o1].x  = cn1[pos] * 12.0;
          gpu->pbLJTerm->_pSysData[o2].x  = cn1[pos] * 12.0;
          gpu->pbLJTerm->_pSysData[o1].y  = cn2[pos] * 6.0;
          gpu->pbLJTerm->_pSysData[o2].y  = cn2[pos] * 6.0;
        }
        else {
          gpu->pbLJTerm->_pSysData[o1].x  = cn1[pos];
          gpu->pbLJTerm->_pSysData[o2].x  = cn1[pos];
          gpu->pbLJTerm->_pSysData[o1].y  = cn2[pos];
          gpu->pbLJTerm->_pSysData[o2].y  = cn2[pos];
        }
        if (*lj1264 == 1) {
          int ooffset = gpu->sim.LJOffset*(LJMultiplier - 1);
          gpu->pbLJTerm->_pSysData[o1 + ooffset].x = cn6[pos] * 4.0;
          gpu->pbLJTerm->_pSysData[o2 + ooffset].x = cn6[pos] * 4.0;
        }
      }

      // If doing TI, calculate 1/sigma^6 and 4*epsilon
      if (gpu->sim.ti_mode >= 2) {

        // Contingency for pos >= 0: pbLJTerm is zeroed beforehand,
        // which would lead to division by 0 for softcore hydrogens
        if (pos >= 0) {
          double oneOverSigma6            = cn2[pos] / cn1[pos];
          double foureps                  = cn2[pos] * oneOverSigma6;
          o1                             += gpu->sim.LJOffset;
          o2                             += gpu->sim.LJOffset;
          gpu->pbLJTerm->_pSysData[o1].x  = oneOverSigma6;
          gpu->pbLJTerm->_pSysData[o2].x  = oneOverSigma6;
          gpu->pbLJTerm->_pSysData[o1].y  = foureps;
          gpu->pbLJTerm->_pSysData[o2].y  = foureps;
        }
        else {
          o1                             += gpu->sim.LJOffset;
          o2                             += gpu->sim.LJOffset;
          gpu->pbLJTerm->_pSysData[o1].x  = 1.0;
          gpu->pbLJTerm->_pSysData[o2].x  = 1.0;
          gpu->pbLJTerm->_pSysData[o1].y  = 0.0;
          gpu->pbLJTerm->_pSysData[o2].y  = 0.0;
        }
      }
    }
  }
  if (*plj1264 == 1) {// C4PairwiseCUDA2023 add cn7 and cn8 to the end of pbLJTerm
    int oooffset = gpu->sim.LJOffset*(LJMultiplier);
    for (k=1; k<=*C4Pairwise; k++) {
      gpu->pbLJTerm->_pSysData[oooffset+k*4-4].x = cn7[k*3-3];
      gpu->pbLJTerm->_pSysData[oooffset+k*4-3].x = cn7[k*3-2];
      gpu->pbLJTerm->_pSysData[oooffset+k*4-2].x = cn7[k*3-1];
      gpu->pbLJTerm->_pSysData[oooffset+k*4-1].x = cn8[k];
    }
  }
  

  gpu->pbLJTerm->Upload();
  gpu->sim.pLJTerm = gpu->pbLJTerm->_pDevData;
#if !defined(use_DPFP)
  {
    cudaError_t status;
    if (gpu->sim.texLJTerm != 0) {
      status = cudaDestroyTextureObject(gpu->sim.texLJTerm);
      RTERROR(status, "cudaDestroyTextureObject gpu->sim.texLJTerm failed");
    }
    cudaResourceDesc resDesc;
    cudaTextureDesc texDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceTypeLinear;
    resDesc.res.linear.devPtr = gpu->sim.pLJTerm;
    resDesc.res.linear.sizeInBytes = gpu->sim.LJTerms * sizeof(float2);
    resDesc.res.linear.desc.f = cudaChannelFormatKindFloat;
    resDesc.res.linear.desc.x = 32;
    resDesc.res.linear.desc.y = 32;
    resDesc.res.linear.desc.z = 0;
    resDesc.res.linear.desc.w = 0;
    memset(&texDesc, 0, sizeof(texDesc));
    texDesc.normalizedCoords = 0;
    texDesc.filterMode = cudaFilterModePoint;
    texDesc.addressMode[0] = cudaAddressModeClamp;
    texDesc.readMode = cudaReadModeElementType;
    status = cudaCreateTextureObject(&(gpu->sim.texLJTerm), &resDesc, &texDesc, NULL);
    RTERROR(status, "cudaCreateTextureObject gpu->sim.texLJTerm failed");
  }
#endif
  // Copy atomic nonbond parameters
  for (i = 0; i < gpu->sim.atoms; i++) {
    if (iac[i] >= 1) {
      // Safe-guard the original Lennard-Jones typing
      gpu->pbAtomLJID->_pSysData[i]           = iac[i];
      FloatShift ljid;
      ljid.ui = iac[i];
      gpu->pbAtomChargeSPLJID->_pSysData[i].y = ljid.f;
      gpu->pbAtomSigEps->_pSysData[i].x       = sigma[iac[i] - 1];
      gpu->pbAtomSigEps->_pSysData[i].y       = (PMEFloat)(2.0 * epsilon[iac[i] - 1]);
    }
    else {
      gpu->pbAtomLJID->_pSysData[i] = 0;
      gpu->pbAtomChargeSPLJID->_pSysData[i].y = 0;
      gpu->pbAtomSigEps->_pSysData[i].x = (PMEFloat)0.0;
      gpu->pbAtomSigEps->_pSysData[i].y = (PMEFloat)0.0;
    }
    gpu->pbAtomLJID->_pSysData[i]           = iac[i];
    FloatShift ljid;
    ljid.ui = iac[i];
    gpu->pbAtomChargeSPLJID->_pSysData[i].y = ljid.f;
    if (iac[i] >= 1) {
      gpu->pbAtomSigEps->_pSysData[i].x       = sigma[iac[i] - 1];
      gpu->pbAtomSigEps->_pSysData[i].y       = (PMEFloat)(2.0 * epsilon[iac[i] - 1]);
    }
  }
  for (i = gpu->sim.atoms; i < gpu->sim.paddedNumberOfAtoms; i++) {
    gpu->pbAtomLJID->_pSysData[i]           = (PMEFloat)0;
    gpu->pbAtomChargeSPLJID->_pSysData[i].y = 0;
    gpu->pbAtomSigEps->_pSysData[i].x       = (PMEFloat)0.0;
    gpu->pbAtomSigEps->_pSysData[i].y       = (PMEFloat)0.0;
  }
  gpu->pbAtomSigEps->Upload();
  gpu->pbAtomLJID->Upload();
  gpu->pbAtomChargeSPLJID->Upload();
  // Search for 1/r4 sources
#if !defined(GTI)
  if (*lj1264 == 1 && *plj1264 == 0) {
    // Step 1: make tables of the number of atoms with each Lennard-Jones type
    //         and the total number of interactions that this would be if the
    //         cutoff were infinite.
    int* typcounts = (int*)calloc(*ntypes, sizeof(int));
    for (i = 0; i < gpu->sim.atoms; i++) {
      typcounts[iac[i] - 1] += 1;
    }
    imat r4intr = CreateImat(*ntypes, *ntypes);
    for (i = 1; i <= *ntypes; i++) {
      for (j = 1; j <= i; j++) {
        int pos = ico[*ntypes * (i - 1) + j - 1] - 1;
        if (pos < 0) {
          continue;
        }
        r4intr.map[i-1][j-1] = (cn6[pos] > 1.0e-8) * typcounts[i-1] * typcounts[j-1];
        r4intr.map[j-1][i-1] = (cn6[pos] > 1.0e-8) * typcounts[i-1] * typcounts[j-1];
      }
    }
    int totcov = 0;
    for (i = 0; i < *ntypes; i++) {
      for (j = 0; j < *ntypes; j++) {
        totcov += r4intr.map[i][j];
      }
    }

    // Step 2: take the most "bang for the buck" by iteratively including the atom types
    //         with the most coverage on r4intr and the fewest members on typcounts, until
    //         all interactions involving 1/r4 terms can be accounted for.
    int nacc = 0;
    int* r4list = (int*)malloc((*ntypes) * sizeof(int));
    int nlist = 0;
    while (nacc < totcov) {
      int tbest = 0;
      int bestcov = 0;
      double bestratio;
      for (i = 0; i < *ntypes; i++) {
        double cratio = 0.0;
        int ncov = 0;
        for (j = 0; j < *ntypes; j++) {
          if (j != i) {
            cratio += r4intr.map[i][j] + r4intr.map[j][i];
            ncov += r4intr.map[i][j] + r4intr.map[j][i];
          }
          else {
            cratio += r4intr.map[i][j];
            ncov += r4intr.map[i][j];
          }
        }
        if (typcounts[i] > 0) {
          cratio /= (double)typcounts[i];
        }
        else {
          cratio = 0.0;
        }
        if (i == 0 || cratio > bestratio) {
          bestratio = cratio;
          bestcov = ncov;
          tbest = i;
        }
      }
      nacc += bestcov;
      r4list[nlist] = tbest;
      nlist++;
      for (i = 0; i < *ntypes; i++) {
        r4intr.map[i][tbest] = 0;
        r4intr.map[tbest][i] = 0;
      }
    }

    // Step 3: enumerate all atoms of the critical types
    int nR4 = 0;
    for (i = 0; i < nlist; i++) {
      for (j = 0; j < gpu->sim.atoms; j++) {
        nR4 += (iac[j] == r4list[i]);
      }
    }
    nR4 = (nR4/GRID + 1) * GRID;
    gpu->pbR4sources = new GpuBuffer<int2>(nR4);
    nR4 = 0;
    for (i = 0; i < gpu->sim.atoms; i++) {
      if (fabs(cn6[iac[i]]) > 1.0e-8) {
        gpu->pbR4sources->_pSysData[nR4].x = i;
        gpu->pbR4sources->_pSysData[nR4].y = iac[i];
        nR4++;
      }
    }
    gpu->sim.nR4sources = nR4;
    gpu->pbR4sources->Upload();
    gpu->sim.pR4sources = gpu->pbR4sources->_pDevData;
  }
  // C4PairwiseCUDA2023
  if (*lj1264 == 0 && *plj1264 == 1) {
    
    
    int* typcounts = (int*)calloc(*ntypes, sizeof(int));
    for (i = 0; i < gpu->sim.atoms; i++) {
      typcounts[iac[i] - 1] += 1;
    }
    imat r4intr = CreateImat(*ntypes, *ntypes);
    for (i = 1; i <= *ntypes; i++) {
      for (j = 1; j <= i; j++) {
        int pos = ico[*ntypes * (i - 1) + j - 1] - 1;
        if (pos < 0) {
          continue;
        }
        r4intr.map[i-1][j-1] = (cn6[pos] > 1.0e-8) * typcounts[i-1] * typcounts[j-1];
        r4intr.map[j-1][i-1] = (cn6[pos] > 1.0e-8) * typcounts[i-1] * typcounts[j-1];
      }
    }
    int totcov = 0;
    for (i = 0; i < *ntypes; i++) {
      for (j = 0; j < *ntypes; j++) {
        totcov += r4intr.map[i][j];
      }
    }
    
    
    int nacc = 0;
    int* r4list = (int*)malloc((*ntypes) * sizeof(int));
    int nlist = 0;
    while (nacc < totcov) {
      int tbest = 0;
      int bestcov = 0;
      double bestratio;
      for (i = 0; i < *ntypes; i++) {
        double cratio = 0.0;
        int ncov = 0;
        for (j = 0; j < *ntypes; j++) {
          if (j != i) {
            cratio += r4intr.map[i][j] + r4intr.map[j][i];
            ncov += r4intr.map[i][j] + r4intr.map[j][i];
          }
          else {
            cratio += r4intr.map[i][j];
            ncov += r4intr.map[i][j];
          }
        }
        if (typcounts[i] > 0) {
          cratio /= (double)typcounts[i];
        }
        else {
          cratio = 0.0;
        }
        if (i == 0 || cratio > bestratio) {
          bestratio = cratio;
          bestcov = ncov;
          tbest = i;
        }
      }
      nacc += bestcov;
      r4list[nlist] = tbest;
      nlist++;
      for (i = 0; i < *ntypes; i++) {
        r4intr.map[i][tbest] = 0;
        r4intr.map[tbest][i] = 0;
      }
    }
    
    
    int nR4 = 0;
    for (i = 0; i < nlist; i++) {
      for (j = 0; j < gpu->sim.atoms; j++) {
        nR4 += (iac[j] == r4list[i]);
      }
    }
    nR4 = (nR4/32 + 1) * 32;
    gpu->pbR4sources = new GpuBuffer<int2>(nR4);
    nR4 = 0;
    for (i = 0; i < gpu->sim.atoms; i++) {
      if (fabs(cn6[iac[i]]) > 1.0e-8) {
        gpu->pbR4sources->_pSysData[nR4].x = i;
        gpu->pbR4sources->_pSysData[nR4].y = iac[i];
        nR4++;
      }
    }
    gpu->sim.nR4sources = nR4;
    gpu->pbR4sources->Upload();
    gpu->sim.pR4sources = gpu->pbR4sources->_pDevData;
  }
  // C4PairwiseCUDA2023
  if (*lj1264 == 1 && *plj1264 == 1) {


    int* typcounts = (int*)calloc(*ntypes, sizeof(int));
    for (i = 0; i < gpu->sim.atoms; i++) {
      typcounts[iac[i] - 1] += 1;
    }
    imat r4intr = CreateImat(*ntypes, *ntypes);
    for (i = 1; i <= *ntypes; i++) {
      for (j = 1; j <= i; j++) {
        int pos = ico[*ntypes * (i - 1) + j - 1] - 1;
        if (pos < 0) {
          continue;
        }
        r4intr.map[i-1][j-1] = (cn6[pos] > 1.0e-8) * typcounts[i-1] * typcounts[j-1];
        r4intr.map[j-1][i-1] = (cn6[pos] > 1.0e-8) * typcounts[i-1] * typcounts[j-1];
      }
    }
    int totcov = 0;
    for (i = 0; i < *ntypes; i++) {
      for (j = 0; j < *ntypes; j++) {
        totcov += r4intr.map[i][j];
      }
    }


    int nacc = 0;
    int* r4list = (int*)malloc((*ntypes) * sizeof(int));
    int nlist = 0;
    while (nacc < totcov) {
      int tbest = 0;
      int bestcov = 0;
      double bestratio;
      for (i = 0; i < *ntypes; i++) {
        double cratio = 0.0;
        int ncov = 0;
        for (j = 0; j < *ntypes; j++) {
          if (j != i) {
            cratio += r4intr.map[i][j] + r4intr.map[j][i];
            ncov += r4intr.map[i][j] + r4intr.map[j][i];
          }
          else {
            cratio += r4intr.map[i][j];
            ncov += r4intr.map[i][j];
          }
        }
        if (typcounts[i] > 0) {
          cratio /= (double)typcounts[i];
        }
        else {
          cratio = 0.0;
        }
        if (i == 0 || cratio > bestratio) {
          bestratio = cratio;
          bestcov = ncov;
          tbest = i;
        }
      }
      nacc += bestcov;
      r4list[nlist] = tbest;
      nlist++;
      for (i = 0; i < *ntypes; i++) {
        r4intr.map[i][tbest] = 0;
        r4intr.map[tbest][i] = 0;
      }
    }


    int nR4 = 0;
    for (i = 0; i < nlist; i++) {
      for (j = 0; j < gpu->sim.atoms; j++) {
        nR4 += (iac[j] == r4list[i]);
      }
    }
    nR4 = (nR4/32 + 1) * 32;
    gpu->pbR4sources = new GpuBuffer<int2>(nR4);
    nR4 = 0;
    for (i = 0; i < gpu->sim.atoms; i++) {
      if (fabs(cn6[iac[i]]) > 1.0e-8) {
        gpu->pbR4sources->_pSysData[nR4].x = i;
        gpu->pbR4sources->_pSysData[nR4].y = iac[i];
        nR4++;
      }
    }
    gpu->sim.nR4sources = nR4;
    gpu->pbR4sources->Upload();
    gpu->sim.pR4sources = gpu->pbR4sources->_pDevData;
  }
#endif /* GTI */

  // Copy 1-4 interactions
  gpu->sim.nb14s = nb14s;
  if (nb14s > 0) {
    if (gpu->ntf < 8) {
      for (i = 0; i < nb14s; i++) {
        int parm_idx = cit_nb14[i][2] - 1;
        if(gpu->iphmd > 0) {
          gpu->pbNb141->_pSysData[i].x = gbl_one_scee[parm_idx];
        }
        else {
          gpu->pbNb141->_pSysData[i].x = (parm_idx>=0) ? gbl_one_scee[parm_idx] *
                                       gpu->pbAtomCharge->_pSysData[abs(cit_nb14[i][0]) - 1] *
                                       gpu->pbAtomCharge->_pSysData[abs(cit_nb14[i][1]) - 1] :0.0;
        }
        gpu->pbNb141->_pSysData[i].y = (parm_idx>=0) ? gbl_one_scnb[parm_idx] :0.0 ;
        int tt = (*ntypes)*(iac[abs(cit_nb14[i][0]) - 1] - 1) +
                 (iac[abs(cit_nb14[i][1]) - 1] - 1);

        // Make sure that nbtype is well-behaved
        int nbtype = tt >= 0 ? ico[tt] - 1 : -1;
        if (nbtype >= 0) {
          gpu->pbNb142->_pSysData[i].x = cn114[nbtype];
          gpu->pbNb142->_pSysData[i].y = cn214[nbtype];
        }
        else {
          gpu->pbNb142->_pSysData[i].x = (PMEDouble)0.0;
          gpu->pbNb142->_pSysData[i].y = (PMEDouble)0.0;
        }
        gpu->pbNb14ID->_pSysData[i].x = abs(cit_nb14[i][0]) - 1;
        gpu->pbNb14ID->_pSysData[i].y = abs(cit_nb14[i][1]) - 1;
      }
    }
    gpu->pbNb141->Upload();
    gpu->pbNb142->Upload();
    gpu->pbNb14ID->Upload();

    // Set constants
    gpu->sim.pNb141 = gpu->pbNb141->_pDevData;
    gpu->sim.pNb142 = gpu->pbNb142->_pDevData;
    gpu->sim.pNb14ID = gpu->pbNb14ID->_pDevData;
  }

  // Set up charge refresh buffer for constant pH or constant Redox potential
  if (gpu->sim.icnstph != 0 || gpu->sim.icnste != 0) {
    delete gpu->pbChargeRefreshBuffer;
    gpu->pbChargeRefreshBuffer = new GpuBuffer<double>(gpu->sim.atoms, false);
    gpu->sim.pChargeRefreshBuffer = gpu->pbChargeRefreshBuffer->_pDevData;
  }

  // Upload constant data
  gpuCopyConstants();

  // Delete temporary arrays
  delete[] sigma;
  delete[] epsilon;
}

//---------------------------------------------------------------------------------------------
// gpu_init_extra_pnts_nb14_: initialize 1-4 interactions for extra points on the GPU.  Extra
//                            points inherit the 1-4 interactions of their parent atoms.
//
// Arguments:
//   frame_cnt:    the number of extra points frames
//   ep_frames:    extra point frame parameters
//   ep_lcl_crd:   the local coordinates for placing extra points
//---------------------------------------------------------------------------------------------
extern "C" void gpu_init_extra_pnts_nb14_(int* frame_cnt, ep_frame_rec* ep_frames,
                                          double ep_lcl_crd[][2][3])
{
  PRINTMETHOD("gpu_init_extra_pnts_nb14");
  int EP11s = 0;
  int EP12s = 0;
  int EP21s = 0;
  int EP22s = 0;
  int EPCustomCount = 0;

  // Count each type of extra point frame
  for (int i = 0; i < *frame_cnt; i++) {
    if (ep_frames[i].ep_cnt == 1 && ep_frames[i].type <= 3) {
      if (ep_frames[i].type == 1) {
        EP11s++;
      }
      else if (ep_frames[i].type == 2) {
        EP12s++;
      }
    }
    else if (ep_frames[i].ep_cnt == 2) {
      if (ep_frames[i].type == 1) {
        EP21s++;
      }
      else if (ep_frames[i].type == 2) {
        EP22s++;
      }
    }
    else {
      EPCustomCount++;
    }
  }
  int EP11Stride = ((EP11s + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits;
  int EP12Stride = ((EP12s + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits;
  int EP21Stride = ((EP21s + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits;
  int EP22Stride = ((EP22s + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits;
  int EPCustomStride = ((EPCustomCount + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) <<
                       gpu->sim.gridBits;

  // Delete existing extra point frames
  delete gpu->pbExtraPoint11Frame;
  delete gpu->pbExtraPoint11Index;
  delete gpu->pbExtraPoint11;
  delete gpu->pbExtraPoint12Frame;
  delete gpu->pbExtraPoint12Index;
  delete gpu->pbExtraPoint12;
  delete gpu->pbExtraPoint21Frame;
  delete gpu->pbExtraPoint21Index;
  delete gpu->pbExtraPoint21;
  delete gpu->pbExtraPoint22Frame;
  delete gpu->pbExtraPoint22Index;
  delete gpu->pbExtraPoint22;

  // Allocate extra point frames
  gpu->pbExtraPoint11Frame     = new GpuBuffer<int4>(EP11s);
  gpu->pbExtraPoint11Index     = new GpuBuffer<int>(EP11s);
  gpu->pbExtraPoint11          = new GpuBuffer<double>(3 * EP11Stride);
  gpu->pbExtraPoint12Frame     = new GpuBuffer<int4>(EP12s);
  gpu->pbExtraPoint12Index     = new GpuBuffer<int>(EP12s);
  gpu->pbExtraPoint12          = new GpuBuffer<double>(3 * EP12Stride);
  gpu->pbExtraPoint21Frame     = new GpuBuffer<int4>(EP21s);
  gpu->pbExtraPoint21Index     = new GpuBuffer<int2>(EP21s);
  gpu->pbExtraPoint21          = new GpuBuffer<double>(6 * EP21Stride);
  gpu->pbExtraPoint22Frame     = new GpuBuffer<int4>(EP22s);
  gpu->pbExtraPoint22Index     = new GpuBuffer<int2>(EP22s);
  gpu->pbExtraPoint22          = new GpuBuffer<double>(6 * EP22Stride);
  gpu->pbExtraPointCustomFrame = new GpuBuffer<int4>(EPCustomCount);
  gpu->pbExtraPointCustomIndex = new GpuBuffer<int>(EPCustomCount);
  gpu->pbExtraPointCustom      = new GpuBuffer<double>(6 * EPCustomStride);

  // Copy Extra point data
  EP11s = 0;
  EP12s = 0;
  EP21s = 0;
  EP22s = 0;
  int epcc = 0;
  for (int i = 0; i < *frame_cnt; i++) {
    if (ep_frames[i].ep_cnt == 1 && ep_frames[i].type <= 3) {
      if (ep_frames[i].type == 1) {
        gpu->pbExtraPoint11Index->_pSysData[EP11s] = ep_frames[i].extra_pnt[0] - 1;
        gpu->pbExtraPoint11Frame->_pSysData[EP11s].x = ep_frames[i].parent_atm - 1;
        gpu->pbExtraPoint11Frame->_pSysData[EP11s].y = ep_frames[i].frame_atm1 - 1;
        gpu->pbExtraPoint11Frame->_pSysData[EP11s].z = ep_frames[i].frame_atm2 - 1;
        gpu->pbExtraPoint11Frame->_pSysData[EP11s].w = ep_frames[i].frame_atm3 - 1;
        gpu->pbExtraPoint11->_pSysData[EP11s] = ep_lcl_crd[i][0][0];
        gpu->pbExtraPoint11->_pSysData[EP11s + EP11Stride] = ep_lcl_crd[i][0][1];
        gpu->pbExtraPoint11->_pSysData[EP11s + EP11Stride * 2] = ep_lcl_crd[i][0][2];
        EP11s++;
      }
      else if (ep_frames[i].type == 2) {
        gpu->pbExtraPoint12Index->_pSysData[EP12s] = ep_frames[i].extra_pnt[0] - 1;
        gpu->pbExtraPoint12Frame->_pSysData[EP12s].x = ep_frames[i].parent_atm - 1;
        gpu->pbExtraPoint12Frame->_pSysData[EP12s].y = ep_frames[i].frame_atm1 - 1;
        gpu->pbExtraPoint12Frame->_pSysData[EP12s].z = ep_frames[i].frame_atm2 - 1;
        gpu->pbExtraPoint12Frame->_pSysData[EP12s].w = ep_frames[i].frame_atm3 - 1;
        gpu->pbExtraPoint12->_pSysData[EP12s] = ep_lcl_crd[i][0][0];
        gpu->pbExtraPoint12->_pSysData[EP12s + EP12Stride] = ep_lcl_crd[i][0][1];
        gpu->pbExtraPoint12->_pSysData[EP12s + EP12Stride * 2] = ep_lcl_crd[i][0][2];
        EP12s++;
      }
    }
    else if (ep_frames[i].ep_cnt == 2) {
      if (ep_frames[i].type == 1) {
        gpu->pbExtraPoint21Index->_pSysData[EP21s].x = ep_frames[i].extra_pnt[0] - 1;
        gpu->pbExtraPoint21Index->_pSysData[EP21s].y = ep_frames[i].extra_pnt[1] - 1;
        gpu->pbExtraPoint21Frame->_pSysData[EP21s].x = ep_frames[i].parent_atm - 1;
        gpu->pbExtraPoint21Frame->_pSysData[EP21s].y = ep_frames[i].frame_atm1 - 1;
        gpu->pbExtraPoint21Frame->_pSysData[EP21s].z = ep_frames[i].frame_atm2 - 1;
        gpu->pbExtraPoint21Frame->_pSysData[EP21s].w = ep_frames[i].frame_atm3 - 1;
        gpu->pbExtraPoint21->_pSysData[EP21s] = ep_lcl_crd[i][0][0];
        gpu->pbExtraPoint21->_pSysData[EP21s + EP21Stride] = ep_lcl_crd[i][0][1];
        gpu->pbExtraPoint21->_pSysData[EP21s + EP21Stride * 2] = ep_lcl_crd[i][0][2];
        gpu->pbExtraPoint21->_pSysData[EP21s + EP21Stride * 3] = ep_lcl_crd[i][1][0];
        gpu->pbExtraPoint21->_pSysData[EP21s + EP21Stride * 4] = ep_lcl_crd[i][1][1];
        gpu->pbExtraPoint21->_pSysData[EP21s + EP21Stride * 5] = ep_lcl_crd[i][1][2];
        EP21s++;
      }
      else if (ep_frames[i].type == 2) {
        gpu->pbExtraPoint22Index->_pSysData[EP22s].x = ep_frames[i].extra_pnt[0] - 1;
        gpu->pbExtraPoint22Index->_pSysData[EP22s].y = ep_frames[i].extra_pnt[1] - 1;
        gpu->pbExtraPoint22Frame->_pSysData[EP22s].x = ep_frames[i].parent_atm - 1;
        gpu->pbExtraPoint22Frame->_pSysData[EP22s].y = ep_frames[i].frame_atm1 - 1;
        gpu->pbExtraPoint22Frame->_pSysData[EP22s].z = ep_frames[i].frame_atm2 - 1;
        gpu->pbExtraPoint22Frame->_pSysData[EP22s].w = ep_frames[i].frame_atm3 - 1;
        gpu->pbExtraPoint22->_pSysData[EP22s] = ep_lcl_crd[i][0][0];
        gpu->pbExtraPoint22->_pSysData[EP22s + EP22Stride] = ep_lcl_crd[i][0][1];
        gpu->pbExtraPoint22->_pSysData[EP22s + EP22Stride * 2] = ep_lcl_crd[i][0][2];
        gpu->pbExtraPoint22->_pSysData[EP22s + EP22Stride * 3] = ep_lcl_crd[i][1][0];
        gpu->pbExtraPoint22->_pSysData[EP22s + EP22Stride * 4] = ep_lcl_crd[i][1][1];
        gpu->pbExtraPoint22->_pSysData[EP22s + EP22Stride * 5] = ep_lcl_crd[i][1][2];
        EP22s++;
      }
    }
    else {

      // Configure custom virtual site (EP) data for all type numbers >= 4.
      // Encode the frame type in the x field of the int4 for compact memory access.
      gpu->pbExtraPointCustomIndex->_pSysData[epcc]   = ep_frames[i].extra_pnt[0] - 1;
      gpu->pbExtraPointCustomFrame->_pSysData[epcc].x = ((ep_frames[i].parent_atm - 1) |
                                                         (ep_frames[i].type << 24));
      gpu->pbExtraPointCustomFrame->_pSysData[epcc].y = ep_frames[i].frame_atm1 - 1;
      gpu->pbExtraPointCustomFrame->_pSysData[epcc].z = -1;
      gpu->pbExtraPointCustomFrame->_pSysData[epcc].w = -1;
      gpu->pbExtraPointCustom->_pSysData[epcc] = ep_lcl_crd[i][0][0];
      if (ep_frames[i].type == 4 || ep_frames[i].type == 9) {
	epcc++;
      }
      else if (ep_frames[i].type == 5 || ep_frames[i].type == 6 || ep_frames[i].type == 7) {
        gpu->pbExtraPointCustomFrame->_pSysData[epcc].z = ep_frames[i].frame_atm2 - 1;
        gpu->pbExtraPointCustom->_pSysData[epcc +   EPCustomStride] = ep_lcl_crd[i][0][1];
	epcc++;
      }
      else if (ep_frames[i].type == 8) {
	gpu->pbExtraPointCustomFrame->_pSysData[epcc].z = ep_frames[i].frame_atm2 - 1;
	gpu->pbExtraPointCustom->_pSysData[epcc +   EPCustomStride] = ep_lcl_crd[i][0][1];
	gpu->pbExtraPointCustom->_pSysData[epcc + 2*EPCustomStride] = ep_lcl_crd[i][0][2];
	epcc++;
      }
      else if (ep_frames[i].type == 10) {
	gpu->pbExtraPointCustomFrame->_pSysData[epcc].z = ep_frames[i].frame_atm2 - 1;
	gpu->pbExtraPointCustomFrame->_pSysData[epcc].w = ep_frames[i].frame_atm3 - 1;
	gpu->pbExtraPointCustom->_pSysData[epcc +   EPCustomStride] = ep_lcl_crd[i][0][1];
	gpu->pbExtraPointCustom->_pSysData[epcc + 2*EPCustomStride] = ep_lcl_crd[i][0][2];
	epcc++;
      }
    }
  }

  // Upload constants and data
  gpu->pbExtraPoint11Frame->Upload();
  gpu->pbExtraPoint11Index->Upload();
  gpu->pbExtraPoint11->Upload();
  gpu->pbExtraPoint12Frame->Upload();
  gpu->pbExtraPoint12Index->Upload();
  gpu->pbExtraPoint12->Upload();
  gpu->pbExtraPoint21Frame->Upload();
  gpu->pbExtraPoint21Index->Upload();
  gpu->pbExtraPoint21->Upload();
  gpu->pbExtraPoint22Frame->Upload();
  gpu->pbExtraPoint22Index->Upload();
  gpu->pbExtraPoint22->Upload();
  gpu->pbExtraPointCustomFrame->Upload();
  gpu->pbExtraPointCustomIndex->Upload();
  gpu->pbExtraPointCustom->Upload();
  gpu->sim.EPs    = *frame_cnt;
  gpu->sim.EP11s  = EP11s;
  gpu->sim.EP12s  = EP12s;
  gpu->sim.EP21s  = EP21s;
  gpu->sim.EP22s  = EP22s;
  gpu->sim.EPCustomCount = EPCustomCount;
  gpu->sim.EP11Offset     =                       EP11Stride;
  gpu->sim.EP12Offset     = gpu->sim.EP11Offset + EP12Stride;
  gpu->sim.EP21Offset     = gpu->sim.EP12Offset + EP21Stride;
  gpu->sim.EP22Offset     = gpu->sim.EP21Offset + EP22Stride;
  gpu->sim.EPCustomOffset = gpu->sim.EP22Offset + EPCustomStride;
  gpu->sim.pExtraPoint11Frame = gpu->pbExtraPoint11Frame->_pDevData;
  gpu->sim.pExtraPoint11Index = gpu->pbExtraPoint11Index->_pDevData;
  gpu->sim.pExtraPoint11X = gpu->pbExtraPoint11->_pDevData;
  gpu->sim.pExtraPoint11Y = gpu->pbExtraPoint11->_pDevData + EP11Stride;
  gpu->sim.pExtraPoint11Z = gpu->pbExtraPoint11->_pDevData + EP11Stride * 2;
  gpu->sim.pExtraPoint12Frame = gpu->pbExtraPoint12Frame->_pDevData;
  gpu->sim.pExtraPoint12Index = gpu->pbExtraPoint12Index->_pDevData;
  gpu->sim.pExtraPoint12X = gpu->pbExtraPoint12->_pDevData;
  gpu->sim.pExtraPoint12Y = gpu->pbExtraPoint12->_pDevData + EP12Stride;
  gpu->sim.pExtraPoint12Z = gpu->pbExtraPoint12->_pDevData + EP12Stride * 2;
  gpu->sim.pExtraPoint21Frame = gpu->pbExtraPoint21Frame->_pDevData;
  gpu->sim.pExtraPoint21Index = gpu->pbExtraPoint21Index->_pDevData;
  gpu->sim.pExtraPoint21X1 = gpu->pbExtraPoint21->_pDevData;
  gpu->sim.pExtraPoint21Y1 = gpu->pbExtraPoint21->_pDevData + EP21Stride;
  gpu->sim.pExtraPoint21Z1 = gpu->pbExtraPoint21->_pDevData + EP21Stride * 2;
  gpu->sim.pExtraPoint21X2 = gpu->pbExtraPoint21->_pDevData + EP21Stride * 3;
  gpu->sim.pExtraPoint21Y2 = gpu->pbExtraPoint21->_pDevData + EP21Stride * 4;
  gpu->sim.pExtraPoint21Z2 = gpu->pbExtraPoint21->_pDevData + EP21Stride * 5;
  gpu->sim.pExtraPoint22Frame = gpu->pbExtraPoint22Frame->_pDevData;
  gpu->sim.pExtraPoint22Index = gpu->pbExtraPoint22Index->_pDevData;
  gpu->sim.pExtraPoint22X1 = gpu->pbExtraPoint22->_pDevData;
  gpu->sim.pExtraPoint22Y1 = gpu->pbExtraPoint22->_pDevData + EP22Stride;
  gpu->sim.pExtraPoint22Z1 = gpu->pbExtraPoint22->_pDevData + EP22Stride * 2;
  gpu->sim.pExtraPoint22X2 = gpu->pbExtraPoint22->_pDevData + EP22Stride * 3;
  gpu->sim.pExtraPoint22Y2 = gpu->pbExtraPoint22->_pDevData + EP22Stride * 4;
  gpu->sim.pExtraPoint22Z2 = gpu->pbExtraPoint22->_pDevData + EP22Stride * 5;
  gpu->sim.pExtraPointCustomFrame = gpu->pbExtraPointCustomFrame->_pDevData;
  gpu->sim.pExtraPointCustomIndex = gpu->pbExtraPointCustomIndex->_pDevData;
  gpu->sim.pExtraPointCustomD1 = gpu->pbExtraPointCustom->_pDevData;
  gpu->sim.pExtraPointCustomD2 = gpu->pbExtraPointCustom->_pDevData + EPCustomStride;
  gpu->sim.pExtraPointCustomD3 = gpu->pbExtraPointCustom->_pDevData + EPCustomStride * 2;
  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gpu_constraints_setup_: setup for position constraints on the GPU.  This is called at
//                         initialization in pmemd.F90.
//
//  Arguments:
//    natc:        the number of atom constraints
//    atm_jrc:     restraint group atom selection array
//    atm_weight:  weights (stiffnesses) by which to constrain
//    atm_xc:      atom position coordinates for constraints array
//---------------------------------------------------------------------------------------------
extern "C" void gpu_constraints_setup_(int* natc, int* atm_jrc, double* atm_weight,
                                       double* atm_xc)
{
  PRINTMETHOD("gpu_constraints_setup");

  // Delete old constraints
  delete gpu->pbConstraint1;
  delete gpu->pbConstraint2;
  delete gpu->pbConstraintID;

  // Punt if any constraint pointer is NULL
  if ((atm_jrc == NULL) || (atm_weight == NULL) || (atm_xc == NULL)) {
    return;
  }

  int constraints = *natc;
#ifdef GVERBOSE
  printf("%d constraints\n", constraints);
#endif

  // Allocate/reallocate constraint buffers
  gpu->pbConstraint1  = new GpuBuffer<PMEDouble2>(constraints);
  gpu->pbConstraint2  = new GpuBuffer<PMEDouble2>(constraints);
  gpu->pbConstraintID = new GpuBuffer<int>(constraints);

  // Copy constraints
  for (int i = 0; i < constraints; i++) {
    int j = atm_jrc[i] - 1;
    gpu->pbConstraint1->_pSysData[i].x = atm_weight[i];
    gpu->pbConstraint1->_pSysData[i].y = atm_xc[3*j    ];
    gpu->pbConstraint2->_pSysData[i].x = atm_xc[3*j + 1];
    gpu->pbConstraint2->_pSysData[i].y = atm_xc[3*j + 2];
    gpu->pbConstraintID->_pSysData[i] = j;
  }
  gpu->pbConstraint1->Upload();
  gpu->pbConstraint2->Upload();
  gpu->pbConstraintID->Upload();

  // Set constants
  gpu->sim.constraints    = constraints;
  gpu->sim.pConstraint1   = gpu->pbConstraint1->_pDevData;
  gpu->sim.pConstraint2   = gpu->pbConstraint2->_pDevData;
  gpu->sim.pConstraintID  = gpu->pbConstraintID->_pDevData;
  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gpu_angles_ub_setup_: setup for Urey-Bradley angles on the GPU
//
// Arguments:
//   angle_ub_cnt:  the number of Urey-Bradley angles
//   angle_ub:      parameters for all Urey-Bradley angles (affected atom indices + index into
//                  parameter arrays)
//   ub_r0:         equilibrium length for the Urey-Bradley spring
//   ub_rk:         spring constants for the Urey-Bradley potentials
//---------------------------------------------------------------------------------------------
extern "C" void gpu_angles_ub_setup_(int* angle_ub_cnt, angle_ub_rec angle_ub[],
                                     double ub_r0[], double ub_rk[])
{
  PRINTMETHOD("gpu_angles_ub_setup");

  delete gpu->pbUBAngle;
  delete gpu->pbUBAngleID;

#ifdef GVERBOSE
  printf("%d Urey Bradley angles\n", *angle_ub_cnt);
#endif
  int UBAngles = *angle_ub_cnt;

  // Allocate/reallocate Urey Bradley angles
  gpu->pbUBAngle = new GpuBuffer<PMEDouble2>(UBAngles);
  gpu->pbUBAngleID = new GpuBuffer<int2>(UBAngles);

  // Add Urey Bradley angles
  UBAngles = 0;
  for (int i = 0; i < *angle_ub_cnt; i++) {
    gpu->pbUBAngle->_pSysData[UBAngles].x = ub_rk[angle_ub[i].parm_idx - 1];
    gpu->pbUBAngle->_pSysData[UBAngles].y = ub_r0[angle_ub[i].parm_idx - 1];
    gpu->pbUBAngleID->_pSysData[UBAngles].x = abs(angle_ub[i].atm_i) - 1;
    gpu->pbUBAngleID->_pSysData[UBAngles].y = abs(angle_ub[i].atm_j) - 1;
    UBAngles++;
  }
  gpu->pbUBAngle->Upload();
  gpu->pbUBAngleID->Upload();

  // Set constants
  gpu->sim.UBAngles = UBAngles;
  gpu->sim.pUBAngle = gpu->pbUBAngle->_pDevData;
  gpu->sim.pUBAngleID = gpu->pbUBAngleID->_pDevData;
  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gpu_dihedrals_imp_setup_: setup for imporper dihedrals on the GPU
//
// Arguments:
//   dihed_imp_cnt:   the number of improper dihedrals in the system
//   dihed_imp:       indices of affected atoms and an index into the parameter arrays for
//                    each improper
//   imp_pk:          amplitudes for improper cosine terms
//   imp_phase:       phase angles for improper cosine terms
//---------------------------------------------------------------------------------------------
extern "C" void  gpu_dihedrals_imp_setup_(int* dihed_imp_cnt, dihed_imp_rec dihed_imp[],
                                          double imp_pk[], double imp_phase[])
{
  PRINTMETHOD("gpu_dihedrals_imp_setup");

  // Delete old improper dihedrals
  delete gpu->pbImpDihedral;
  delete gpu->pbImpDihedralID1;

#ifdef GVERBOSE
  printf("%d improper dihedrals\n", *dihed_imp_cnt);
#endif
  int impDihedrals = *dihed_imp_cnt;

  // Allocate/reallocate improper dihedrals
  gpu->pbImpDihedral = new GpuBuffer<PMEDouble2>(impDihedrals);
  gpu->pbImpDihedralID1 = new GpuBuffer<int4>(impDihedrals);

  // Add improper dihedrals
  impDihedrals = 0;
  for (int i = 0; i < *dihed_imp_cnt; i++) {
    gpu->pbImpDihedral->_pSysData[impDihedrals].x = imp_pk[dihed_imp[i].parm_idx - 1];
    gpu->pbImpDihedral->_pSysData[impDihedrals].y = imp_phase[dihed_imp[i].parm_idx - 1];
    gpu->pbImpDihedralID1->_pSysData[impDihedrals].x = abs(dihed_imp[i].atm_i) - 1;
    gpu->pbImpDihedralID1->_pSysData[impDihedrals].y = abs(dihed_imp[i].atm_j) - 1;
    gpu->pbImpDihedralID1->_pSysData[impDihedrals].z = abs(dihed_imp[i].atm_k) - 1;
    gpu->pbImpDihedralID1->_pSysData[impDihedrals].w = abs(dihed_imp[i].atm_l) - 1;
    impDihedrals++;
  }
  gpu->pbImpDihedral->Upload();
  gpu->pbImpDihedralID1->Upload();

  // Set constants
  gpu->sim.impDihedrals = impDihedrals;
  gpu->sim.pImpDihedral = gpu->pbImpDihedral->_pDevData;
  gpu->sim.pImpDihedralID1 = gpu->pbImpDihedralID1->_pDevData;
  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gpu_cmap_setup_: setup for CMAP terms on the GPU
//
// Arguments:
//   cmap_cnt:         the number of CMAP terms in the system
//   cmap:             parameters for all types of CMAP interactions (length cmap_cnt)
//   cmap_type_count:  the number of CMAP types
//   cmap_grid:        the contents of every CMAP grid (grids are interlaced in this packed
//                     array)
//   cmap_dPhi:        partial derivative along the first dimension ("Phi" if the CMAP applies
//                     to the usual protein backbone atoms)
//   cmap_dPsi:        partial derivative along the second dimension ("Psi" if the CMAP applies
//                     to the usual protein backbone atoms)
//   cmap_dPhi_dPsi:   mixed partial derivative of the CMAP term
//---------------------------------------------------------------------------------------------
extern "C" void gpu_cmap_setup_(int* cmap_cnt, cmap_rec cmap[], int* cmap_type_count,
                                double cmap_grid[], double cmap_dPhi[], double cmap_dPsi[],
                                double cmap_dPhi_dPsi[])
{
  PRINTMETHOD("gpu_cmap_setup");

  // Delete old cmap data
  delete gpu->pbCmapID1;
  delete gpu->pbCmapID2;
  delete gpu->pbCmapType;
  delete gpu->pbCmapEnergy;

#ifdef GVERBOSE
  printf("%d cmaps\n", *cmap_cnt);
#endif
  int cmaps = *cmap_cnt;

  // Allocate/reallocate cmaps and cmap data
  gpu->pbCmapID1 = new GpuBuffer<int4>(cmaps);
  gpu->pbCmapID2 = new GpuBuffer<int>(cmaps);
  gpu->pbCmapType = new GpuBuffer<int>(cmaps);

  // Add cmaps
  cmaps = 0;
  for (int i = 0; i < *cmap_cnt; i++) {
    gpu->pbCmapID1->_pSysData[cmaps].x = abs(cmap[i].atm_i) - 1;
    gpu->pbCmapID1->_pSysData[cmaps].y = abs(cmap[i].atm_j) - 1;
    gpu->pbCmapID1->_pSysData[cmaps].z = abs(cmap[i].atm_k) - 1;
    gpu->pbCmapID1->_pSysData[cmaps].w = abs(cmap[i].atm_l) - 1;
    gpu->pbCmapID2->_pSysData[cmaps] = abs(cmap[i].atm_m) - 1;
    gpu->pbCmapType->_pSysData[cmaps] = cmap[i].parm_idx - 1;
    cmaps++;
  }
  gpu->pbCmapID1->Upload();
  gpu->pbCmapID2->Upload();
  gpu->pbCmapType->Upload();

  // Allocate energy map
  int termReads = (*cmap_type_count * sizeof(PMEDouble) + gpu->readSize - 1) / gpu->readSize;
  int cmapTermStride = (termReads * gpu->readSize) / sizeof(PMEDouble);
  int cmapRowStride = cmapTermStride * CMAP_DIMENSION;
  gpu->pbCmapEnergy = new GpuBuffer<PMEFloat4>(cmapTermStride * CMAP_DIMENSION *
                                               CMAP_DIMENSION);

  // Copy cmap terms
  for (int i = 0; i < *cmap_type_count; i++) {
    PMEFloat4* pCmap = gpu->pbCmapEnergy->_pSysData + i + cmapRowStride + cmapTermStride;
    for (int x = -1; x < (int)CMAP_RESOLUTION + 3; x++) {
      int ix = x;
      if (ix < 0) {
        ix += CMAP_RESOLUTION;
      }
      else if (ix >= (int)CMAP_RESOLUTION) {
        ix -= CMAP_RESOLUTION;
      }
      for (int y = -1; y < (int)CMAP_RESOLUTION + 3; y++) {
        int iy = y;
        if (iy < 0) {
          iy += CMAP_RESOLUTION;
        }
        else if (iy >= (int)CMAP_RESOLUTION) {
          iy -= CMAP_RESOLUTION;
        }
        pCmap[y*cmapRowStride + x*cmapTermStride].x =
          cmap_grid[i + iy*(*cmap_type_count) + ix*CMAP_RESOLUTION*(*cmap_type_count)];
        pCmap[y*cmapRowStride + x*cmapTermStride].y =
          cmap_dPhi[i + iy*(*cmap_type_count) + ix*CMAP_RESOLUTION*(*cmap_type_count)] *
          CMAP_STEP_SIZE;
        pCmap[y*cmapRowStride + x*cmapTermStride].z =
          cmap_dPsi[i + iy*(*cmap_type_count) + ix*CMAP_RESOLUTION*(*cmap_type_count)] *
          CMAP_STEP_SIZE;
        pCmap[y*cmapRowStride + x*cmapTermStride].w =
          cmap_dPhi_dPsi[i + iy*(*cmap_type_count) + ix*CMAP_RESOLUTION*(*cmap_type_count)] *
          CMAP_STEP_SIZE * CMAP_STEP_SIZE;
      }
    }
  }
  gpu->pbCmapEnergy->Upload();

  // Set constants
  gpu->sim.cmaps          = cmaps;
  gpu->sim.pCmapID1       = gpu->pbCmapID1->_pDevData;
  gpu->sim.pCmapID2       = gpu->pbCmapID2->_pDevData;
  gpu->sim.pCmapType      = gpu->pbCmapType->_pDevData;
  gpu->sim.cmapTermStride = cmapTermStride;
  gpu->sim.cmapRowStride  = cmapRowStride;
  gpu->sim.pCmapEnergy    = gpu->pbCmapEnergy->_pDevData + cmapTermStride + cmapRowStride;
  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gpu_shake_setup_: setup for SHAKE constraints on the GPU
//
// Arguments:
//   atm_mass:               masses of all atoms
//   my_nonfastwat_bond_cnt: count of bonds not in fast waters
//   my_nonfastwat_bond_dat: parameters for bonds not in fast waters (indices of affected atoms
//                           plus an index into the bond parameter arrays)
//   my_fastwat_res_cnt:     the number of fast waters to SETTLE / RATTLE
//   iorwat:                 indicator of how the waters are arranged in the topology
//   my_fastwat_res_lst:
//---------------------------------------------------------------------------------------------
extern "C" void gpu_shake_setup_(double atm_mass[], int* my_nonfastwat_bond_cnt,
                                 shake_bond_rec my_nonfastwat_bond_dat[],
                                 int* my_fastwat_res_cnt, int* iorwat,
                                 int my_fastwat_res_lst[], int* ti_shake)
{
  PRINTMETHOD("gpu_shake_setup");

  // Delete any existing Shake constraints
  delete gpu->pbShakeID;
  delete gpu->pbShakeParm;
  delete gpu->pbShakeInvMassH;
  delete gpu->pbFastShakeID;
  delete gpu->pbSlowShakeID1;
  delete gpu->pbSlowShakeID2;
  delete gpu->pbSlowShakeParm;
  delete gpu->pbSlowShakeInvMassH;
  delete gpu->pbSlowTIShakeID1;
  delete gpu->pbSlowTIShakeID2;
  delete gpu->pbSlowTIShakeID3;
  delete gpu->pbSlowTIShakeParm;
  delete gpu->pbSlowTIShakeInvMassH;
  delete gpu->pbShakeParm2;
  delete gpu->pbSlowShakeParm2;
  delete gpu->pbSlowTIShakeParm2;

#ifdef GVERBOSE
  printf("%d nonfast SHAKE constraints, %d fast SHAKE constraints\n", *my_nonfastwat_bond_cnt,
         *my_fastwat_res_cnt);
#endif

  // Determine number of H-bond networks in traditional Shake constraint list
  bool bHMR                  = false;
  int shakeConstraints       = 0;
  int slowShakeConstraints   = 0;
  bool* shakeMap             = new bool[*my_nonfastwat_bond_cnt];
  int4* pShakeID             = new int4[*my_nonfastwat_bond_cnt];
  int* pSlowShakeID1         = new int[*my_nonfastwat_bond_cnt];
  int4* pSlowShakeID2        = new int4[*my_nonfastwat_bond_cnt];
  double2* pShakeParm        = new double2[*my_nonfastwat_bond_cnt];
  double* pShakeInvMassH     = new double[*my_nonfastwat_bond_cnt];
  double2* pSlowShakeParm    = new double2[*my_nonfastwat_bond_cnt];
  double* pSlowShakeInvMassH = new double[*my_nonfastwat_bond_cnt];
  int slowTIShakeConstraints;
  int* pSlowTIShakeID1;
  int4* pSlowTIShakeID2;
  int4* pSlowTIShakeID3;
  double2* pShakeParm2;
  double2* pSlowShakeParm2;
  double2* pSlowTIShakeParm;
  double2* pSlowTIShakeParm2;
  double* pSlowTIShakeInvMassH;
  if(*ti_shake == 2) {
    slowTIShakeConstraints    = 0;
    pSlowTIShakeID1          = new int[*my_nonfastwat_bond_cnt];
    pSlowTIShakeID2         = new int4[*my_nonfastwat_bond_cnt];
    pSlowTIShakeID3         = new int4[*my_nonfastwat_bond_cnt];
    pShakeParm2          = new double2[*my_nonfastwat_bond_cnt];
    pSlowShakeParm2      = new double2[*my_nonfastwat_bond_cnt];
    pSlowTIShakeParm     = new double2[*my_nonfastwat_bond_cnt];
    pSlowTIShakeParm2    = new double2[*my_nonfastwat_bond_cnt];
    pSlowTIShakeInvMassH  = new double[*my_nonfastwat_bond_cnt];
  }

  for (int i = 0; i < *my_nonfastwat_bond_cnt; i++) {
    shakeMap[i] = false;
  }
  bool fail = false;
  for (int i = 0; i < *my_nonfastwat_bond_cnt; i++) {
    if (!shakeMap[i]) {
      shakeMap[i]         = true;
      unsigned int atm_i;
      unsigned int atm_j;
      double mass_i;
      double mass_j;
      double r = my_nonfastwat_bond_dat[i].parm;
      if(*ti_shake == 2)
      {
        double tmp_mass_i;
        double tmp_mass_j;
        if(my_nonfastwat_bond_dat[i].atm_i < -1) {
            tmp_mass_i=atm_mass[abs(my_nonfastwat_bond_dat[i].atm_i + 1)];
        }
        else {
            tmp_mass_i=atm_mass[my_nonfastwat_bond_dat[i].atm_i - 1];
        }
        if(my_nonfastwat_bond_dat[i].atm_j < -1) {
            tmp_mass_j=atm_mass[abs(my_nonfastwat_bond_dat[i].atm_j + 1)];
        }
        else {
            tmp_mass_j=atm_mass[my_nonfastwat_bond_dat[i].atm_j - 1];
        }
        if (tmp_mass_i > tmp_mass_j) {
           atm_i = my_nonfastwat_bond_dat[i].atm_i;
           atm_j = my_nonfastwat_bond_dat[i].atm_j;
           mass_i = tmp_mass_i;
           mass_j = tmp_mass_j;
        }
        else {
           atm_i = my_nonfastwat_bond_dat[i].atm_j;
           atm_j = my_nonfastwat_bond_dat[i].atm_i;
           mass_i = tmp_mass_j;
           mass_j = tmp_mass_i;
        }
      }
      else {
        // Determine central atom in network
        if (atm_mass[my_nonfastwat_bond_dat[i].atm_i - 1] >
            atm_mass[my_nonfastwat_bond_dat[i].atm_j - 1]) {
          atm_i = my_nonfastwat_bond_dat[i].atm_i;
          atm_j = my_nonfastwat_bond_dat[i].atm_j;
          mass_i = atm_mass[my_nonfastwat_bond_dat[i].atm_i - 1];
          mass_j = atm_mass[my_nonfastwat_bond_dat[i].atm_j - 1];
        }
        else {
          atm_i = my_nonfastwat_bond_dat[i].atm_j;
          atm_j = my_nonfastwat_bond_dat[i].atm_i;
          mass_i = atm_mass[my_nonfastwat_bond_dat[i].atm_j - 1];
          mass_j = atm_mass[my_nonfastwat_bond_dat[i].atm_i - 1];
        }
      }

      // Check for presence hydrogen mass repartitioning (HMR)
      if ((mass_j / gpu->sim.massH > 1.0001) || (gpu->sim.massH / mass_j > 1.0001)) {
        bHMR = true;
      }

      int max_h_atoms=4;
      if(*ti_shake == 2) {
        max_h_atoms = 8;
      }

      const int MAX_h_atoms = 8;
      unsigned int atm_i_cnt = 0;
      int atom_jj[MAX_h_atoms];
      double mass_ii[MAX_h_atoms];
      double mass_jj[MAX_h_atoms];
      double r_jj[MAX_h_atoms];
      for (int j = 0; j < max_h_atoms; j++) {
        atom_jj[j] = 0;
        mass_ii[j] = 0.0;
        mass_jj[j] = 0.0;
        r_jj[j] = -10.0;
      }
      for (int j = i + 1; j <  *my_nonfastwat_bond_cnt; j++) {
        if (!shakeMap[j]) {

          // Check for additional members of network
          if (my_nonfastwat_bond_dat[j].atm_i == atm_i) {
            atom_jj[atm_i_cnt] = my_nonfastwat_bond_dat[j].atm_j;
            mass_jj[atm_i_cnt] = atm_mass[my_nonfastwat_bond_dat[j].atm_j - 1];
            if(*ti_shake == 2) mass_ii[atm_i_cnt] = my_nonfastwat_bond_dat[j].mass_i;
            r_jj[atm_i_cnt] = my_nonfastwat_bond_dat[j].parm;
            shakeMap[j] = true;
          }
          else if (my_nonfastwat_bond_dat[j].atm_j == atm_i) {
            atom_jj[atm_i_cnt] = my_nonfastwat_bond_dat[j].atm_i;
            mass_jj[atm_i_cnt] = atm_mass[my_nonfastwat_bond_dat[j].atm_i - 1];
            if(*ti_shake == 2) mass_ii[atm_i_cnt] = my_nonfastwat_bond_dat[j].mass_i;
            r_jj[atm_i_cnt] = my_nonfastwat_bond_dat[j].parm;
            shakeMap[j] = true;
          }

          // Check for inconsistent network
          if ((my_nonfastwat_bond_dat[j].atm_i == atm_j) ||
              (my_nonfastwat_bond_dat[j].atm_j == atm_j)) {
            printf("%d: %d - %d and %d: %d - %d\n", i, my_nonfastwat_bond_dat[i].atm_i,
                   my_nonfastwat_bond_dat[i].atm_j, j, my_nonfastwat_bond_dat[j].atm_i,
                   my_nonfastwat_bond_dat[j].atm_j);
            printf("Hydrogen atom %d appears to have multiple bonds to atoms %d and %d which "
                   "is illegal for SHAKEH.\n", atm_j, atm_i, my_nonfastwat_bond_dat[j].atm_j);
            fail = true;
          }
          if (shakeMap[j]) {
            atm_i_cnt++;
            if ((*ti_shake != 2 && atm_i_cnt == 4) || (*ti_shake == 2 && atm_i_cnt == 8)) {
              printf("Too many hydrogens for a hydrogen network, exiting.\n");
              gpu_shutdown_();
              exit(-1);
            }
          }
        }
      }

      double r_ti2=r_jj[0];
      double mass_ti2=mass_i;
      if(*ti_shake == 2) {
        for(int k=0; k<atm_i_cnt; k++)
        {
          if(r != r_jj[k])
          {
            r_ti2=r_jj[k];
          }
          if(mass_i != mass_ii[k])
          {
            mass_ti2=mass_ii[k];
          }
        }

      }
      else {
        for(int k=0; k<atm_i_cnt; k++)
        {
          if(r != r_jj[k])
          {
            r_ti2=r_jj[k];
          }
        }
      }

      // Construct constraint
      if (atm_i_cnt < 3) {
        pShakeID[shakeConstraints].x = atm_i - 1;
        pShakeID[shakeConstraints].y = atm_j - 1;
        pShakeID[shakeConstraints].z = atom_jj[0] - 1;
        pShakeID[shakeConstraints].w = atom_jj[1] - 1;
        pShakeParm[shakeConstraints].x = 1.0 / mass_i;
        pShakeParm[shakeConstraints].y = r;
        if(*ti_shake == 2) {
          pShakeParm2[shakeConstraints].x = 1.0 / mass_ti2;
          if(r_ti2 < 0.0) {
              r_ti2 = r;
          }
          pShakeParm2[shakeConstraints].y = r_ti2;
        }
        pShakeInvMassH[shakeConstraints] = 1.0 / mass_j;
        shakeConstraints++;
      }
      else if (atm_i_cnt < 4) {
        pSlowShakeID1[slowShakeConstraints] = atm_i - 1;
        pSlowShakeID2[slowShakeConstraints].x = atm_j - 1;
        pSlowShakeID2[slowShakeConstraints].y = atom_jj[0] - 1;
        pSlowShakeID2[slowShakeConstraints].z = atom_jj[1] - 1;
        pSlowShakeID2[slowShakeConstraints].w = atom_jj[2] - 1;
        pSlowShakeParm[slowShakeConstraints].x = 1.0 / mass_i;
        pSlowShakeParm[slowShakeConstraints].y = r;
        if(*ti_shake == 2) {
          pSlowShakeParm2[slowShakeConstraints].x = 1.0 / mass_ti2;
          if(r_ti2 < 0.0) {
              r_ti2 = r;
          }
          pSlowShakeParm2[slowShakeConstraints].y = r_ti2;
        }
        pSlowShakeInvMassH[slowShakeConstraints]= 1.0 / mass_j;
        slowShakeConstraints++;
      }
      else {
        if(*ti_shake == 2) {
          pSlowTIShakeID1[slowTIShakeConstraints] = atm_i - 1;
          pSlowTIShakeID2[slowTIShakeConstraints].x = atm_j - 1;
          pSlowTIShakeID2[slowTIShakeConstraints].y = atom_jj[0] - 1;
          pSlowTIShakeID2[slowTIShakeConstraints].z = atom_jj[1] - 1;
          pSlowTIShakeID2[slowTIShakeConstraints].w = atom_jj[2] - 1;
          pSlowTIShakeID3[slowTIShakeConstraints].x = atom_jj[3] - 1;
          pSlowTIShakeID3[slowTIShakeConstraints].y = atom_jj[4] - 1;
          pSlowTIShakeID3[slowTIShakeConstraints].z = atom_jj[5] - 1;
          pSlowTIShakeID3[slowTIShakeConstraints].w = atom_jj[6] - 1;
          pSlowTIShakeParm[slowTIShakeConstraints].x = 1.0 / mass_i;
          pSlowTIShakeParm[slowTIShakeConstraints].y = r;
          pSlowTIShakeParm2[slowTIShakeConstraints].x = 1.0 / mass_ti2;
          if(r_ti2 < 0.0) {
              r_ti2 = r;
          }
          pSlowTIShakeParm2[slowTIShakeConstraints].y = r_ti2;
          pSlowTIShakeInvMassH[slowTIShakeConstraints]= 1.0 / mass_j;
          slowTIShakeConstraints++;
        }
      }
    }

    // Exit if one or more illegal SHAKEH networks detected
    if (fail) {
      printf("Exiting due to the presence of inconsistent SHAKEH hydrogen clusters.\n");
      gpu_shutdown_();
      exit(-1);
    }
  }

  // Allocate and copy standard SHAKE constraints to gpu RAM
  gpu->sim.shakeConstraints = shakeConstraints;
  gpu->pbShakeID          = new GpuBuffer<int4>(shakeConstraints);
  gpu->pbShakeParm        = new GpuBuffer<double2>(shakeConstraints);
  gpu->sim.pShakeID       = gpu->pbShakeID->_pDevData;
  gpu->sim.pShakeParm     = gpu->pbShakeParm->_pDevData;
  gpu->sim.tishake        = *ti_shake;
  if(*ti_shake == 2)
  {
    gpu->pbShakeParm2       = new GpuBuffer<double2>(shakeConstraints);
    gpu->sim.pShakeParm2    = gpu->pbShakeParm2->_pDevData;
  }
  gpu->bUseHMR            = bHMR;
  if (bHMR) {
    gpu->pbShakeInvMassH = new GpuBuffer<double>(shakeConstraints);
    gpu->sim.pShakeInvMassH = gpu->pbShakeInvMassH->_pDevData;
  }
  else {
    gpu->pbShakeInvMassH = NULL;
    gpu->sim.pShakeInvMassH = NULL;
  }
  for (int i = 0; i < shakeConstraints; i++) {
    gpu->pbShakeID->_pSysData[i] = pShakeID[i];
    gpu->pbShakeParm->_pSysData[i] = pShakeParm[i];
    if(*ti_shake == 2) {
      gpu->pbShakeParm2->_pSysData[i] = pShakeParm2[i];
    }
    if (bHMR) {
      gpu->pbShakeInvMassH->_pSysData[i] = pShakeInvMassH[i];
    }
  }
  gpu->pbShakeID->Upload();
  gpu->pbShakeParm->Upload();
  if(*ti_shake == 2) {
    gpu->pbShakeParm2->Upload();
  }
  if (bHMR) {
    gpu->pbShakeInvMassH->Upload();
  }

  // Allocate and copy rarely occurring sp3-hybridized SHAKE constraints to GPU RAM
  gpu->sim.slowShakeConstraints = slowShakeConstraints;
  gpu->pbSlowShakeID1     = new GpuBuffer<int>(slowShakeConstraints + 10); // extra 10 elements reserved for possible TI addition
  gpu->pbSlowShakeID2     = new GpuBuffer<int4>(slowShakeConstraints + 10);
  gpu->pbSlowShakeParm    = new GpuBuffer<double2>(slowShakeConstraints + 10 );
  gpu->sim.pSlowShakeID1  = gpu->pbSlowShakeID1->_pDevData;
  gpu->sim.pSlowShakeID2  = gpu->pbSlowShakeID2->_pDevData;
  gpu->sim.pSlowShakeParm = gpu->pbSlowShakeParm->_pDevData;
  if(*ti_shake == 2) {
    gpu->pbSlowShakeParm2    = new GpuBuffer<double2>(slowShakeConstraints + 10 );
    gpu->sim.pSlowShakeParm2 = gpu->pbSlowShakeParm2->_pDevData;
  }
  if (bHMR) {
    gpu->pbSlowShakeInvMassH = new GpuBuffer<double>(slowShakeConstraints+10);
    gpu->sim.pSlowShakeInvMassH = gpu->pbSlowShakeInvMassH->_pDevData;
  }
  else {
    gpu->pbSlowShakeInvMassH = NULL;
    gpu->sim.pSlowShakeInvMassH = NULL;
  }

  for (int i = 0; i < slowShakeConstraints; i++) {
    gpu->pbSlowShakeID1->_pSysData[i] = pSlowShakeID1[i];
    gpu->pbSlowShakeID2->_pSysData[i] = pSlowShakeID2[i];
    gpu->pbSlowShakeParm->_pSysData[i] = pSlowShakeParm[i];
    if(*ti_shake == 2) {
      gpu->pbSlowShakeParm2->_pSysData[i] = pSlowShakeParm2[i];
    }
    if (bHMR) {
      gpu->pbSlowShakeInvMassH->_pSysData[i] = pSlowShakeInvMassH[i];
    }
  }
  gpu->pbSlowShakeID1->Upload();
  gpu->pbSlowShakeID2->Upload();
  gpu->pbSlowShakeParm->Upload();
  if(*ti_shake == 2) {
    gpu->pbSlowShakeParm2->Upload();
  }
  if (bHMR) {
    gpu->pbSlowShakeInvMassH->Upload();
  }

  if(*ti_shake == 2) {
    // Allocate and copy even more rare TI SHAKE constraints to GPU RAM
    gpu->sim.slowTIShakeConstraints = slowTIShakeConstraints;
    gpu->pbSlowTIShakeID1      = new GpuBuffer<int>(slowTIShakeConstraints);
    gpu->pbSlowTIShakeID2      = new GpuBuffer<int4>(slowTIShakeConstraints);
    gpu->pbSlowTIShakeID3      = new GpuBuffer<int4>(slowTIShakeConstraints);
    gpu->pbSlowTIShakeParm     = new GpuBuffer<double2>(slowTIShakeConstraints);
    gpu->pbSlowTIShakeParm2    = new GpuBuffer<double2>(slowTIShakeConstraints);
    gpu->sim.pSlowTIShakeID1   = gpu->pbSlowTIShakeID1->_pDevData;
    gpu->sim.pSlowTIShakeID2   = gpu->pbSlowTIShakeID2->_pDevData;
    gpu->sim.pSlowTIShakeID3   = gpu->pbSlowTIShakeID3->_pDevData;
    gpu->sim.pSlowTIShakeParm  = gpu->pbSlowTIShakeParm->_pDevData;
    gpu->sim.pSlowTIShakeParm2 = gpu->pbSlowTIShakeParm2->_pDevData;
    if (bHMR) {
      gpu->pbSlowTIShakeInvMassH = new GpuBuffer<double>(slowTIShakeConstraints);
      gpu->sim.pSlowTIShakeInvMassH = gpu->pbSlowTIShakeInvMassH->_pDevData;
    }
    else {
      gpu->pbSlowTIShakeInvMassH = NULL;
      gpu->sim.pSlowTIShakeInvMassH = NULL;
    }

    for (int i = 0; i < slowTIShakeConstraints; i++) {
      gpu->pbSlowTIShakeID1->_pSysData[i] = pSlowTIShakeID1[i];
      gpu->pbSlowTIShakeID2->_pSysData[i] = pSlowTIShakeID2[i];
      gpu->pbSlowTIShakeID3->_pSysData[i] = pSlowTIShakeID3[i];
      gpu->pbSlowTIShakeParm->_pSysData[i] = pSlowTIShakeParm[i];
      gpu->pbSlowTIShakeParm2->_pSysData[i] = pSlowTIShakeParm2[i];
      if (bHMR) {
        gpu->pbSlowTIShakeInvMassH->_pSysData[i] = pSlowTIShakeInvMassH[i];
      }
    }
    gpu->pbSlowTIShakeID1->Upload();
    gpu->pbSlowTIShakeID2->Upload();
    gpu->pbSlowTIShakeID3->Upload();
    gpu->pbSlowTIShakeParm->Upload();
    gpu->pbSlowTIShakeParm2->Upload();
    if (bHMR) {
      gpu->pbSlowTIShakeInvMassH->Upload();
    }
  }

  int ind1, ind2, ind3;
  switch (*iorwat) {
    case 1:
      ind1 = 0;
      ind2 = 1;
      ind3 = 2;
      break;
    case 2:
      ind1 = 1;
      ind2 = 2;
      ind3 = 0;
      break;
    default:
      ind1 = 2;
      ind2 = 0;
      ind3 = 1;
      break;
  }

  // Set up fast SHAKE constraints
  gpu->pbFastShakeID = new GpuBuffer<int4>(*my_fastwat_res_cnt);
  for (int i = 0; i < *my_fastwat_res_cnt; i++) {
    int atm_1 = my_fastwat_res_lst[i] + ind1 - 1;
    int atm_2 = my_fastwat_res_lst[i] + ind2 - 1;
    int atm_3 = my_fastwat_res_lst[i] + ind3 - 1;
    gpu->pbFastShakeID->_pSysData[i].x = atm_1;
    gpu->pbFastShakeID->_pSysData[i].y = atm_2;
    gpu->pbFastShakeID->_pSysData[i].z = atm_3;
    gpu->pbFastShakeID->_pSysData[i].w = -1;
  }
  gpu->sim.fastShakeConstraints = *my_fastwat_res_cnt;
  gpu->pbFastShakeID->Upload();
  gpu->sim.pFastShakeID   = gpu->pbFastShakeID->_pDevData;

  // Calculate SHAKE offsets
  gpu->sim.shakeOffset = (((gpu->sim.shakeConstraints + 63) >> 6) << 6);
  gpu->sim.fastShakeOffset = (((gpu->sim.fastShakeConstraints + 63) >> 6) << 6) +
                             gpu->sim.shakeOffset;
  gpu->sim.slowShakeOffset = (((gpu->sim.slowShakeConstraints + 63) >> 6) << 6) +
                             gpu->sim.fastShakeOffset;
  if(*ti_shake == 2) {
    gpu->sim.slowTIShakeOffset = (((gpu->sim.slowTIShakeConstraints + 63) >> 6) << 6) +
                               gpu->sim.slowShakeOffset;
  }
  gpuCopyConstants();

  // Release temporary space
  delete[] shakeMap;
  delete[] pShakeID;
  delete[] pShakeParm;
  delete[] pShakeInvMassH;
  delete[] pSlowShakeID1;
  delete[] pSlowShakeID2;
  delete[] pSlowShakeParm;
  delete[] pSlowShakeInvMassH;
  if(*ti_shake == 2)
  {
    delete[] pShakeParm2;
    delete[] pSlowShakeParm2;
    delete[] pSlowTIShakeID1;
    delete[] pSlowTIShakeID2;
    delete[] pSlowTIShakeID3;
    delete[] pSlowTIShakeParm;
    delete[] pSlowTIShakeParm2;
    delete[] pSlowTIShakeInvMassH;
  }
}


//---------------------------------------------------------------------------------------------
// gpu_shake_ti_setup_: setup for modifying TI SHAKE constraints on the GPU
//
// Arguments:

//---------------------------------------------------------------------------------------------
extern "C" void gpu_shake_ti_setup_(double atm_mass[], int ti_lst[][3], int ti_sc_lst[], double ti_weights[]) {

  double t0 = ti_weights[0];
  double t1 = ti_weights[1];
  std::vector<bool> picked(gpu->sim.shakeConstraints,false);

  for (int i = 0; i < gpu->sim.shakeConstraints; i++) {
    if (!picked[i]) {
      int iatom = gpu->pbShakeID->_pSysData[i].x;
      if (ti_lst[iatom][2] > 0) continue;

      int iatom0 = gpu->pbShakeID->_pSysData[i].y;
      int iatom1 = gpu->pbShakeID->_pSysData[i].z;
      int iatom2 = gpu->pbShakeID->_pSysData[i].w;

      bool b0 = (iatom0 >= 0 && ti_sc_lst[iatom0] > 0);
      bool b1 = (iatom1 >= 0 && ti_sc_lst[iatom1] > 0);
      bool b2 = (iatom2 >= 0 && ti_sc_lst[iatom2] > 0);

      if (ti_sc_lst[iatom] == 0 && (b0 || b1 || b2)) {
        unsigned int myRegion = (ti_lst[iatom][0] > 0) ? 0 : ((ti_lst[iatom][1] > 0) ? 1 : -1);
        int pos = -1;
        for (int k = 0; k < gpu->sim.numberTICommonPairs; k++) {
          bool f0 = ((myRegion == 0) && (gpu->pbTICommonPair->_pSysData[k].x == iatom));
          bool f1 = ((myRegion == 1) && (gpu->pbTICommonPair->_pSysData[k].y == iatom));
          if ( f0 || f1){
            pos = k;
            break;
          }
        }
        if (pos >= 0) {
          gpu->pbTICommonPair->_pSysData[pos].z = myRegion;
          int j_index = -1;
          for (int j = i + 1; j < gpu->sim.shakeConstraints; j++) {
            if (!picked[j]) {
              bool f0 = ((myRegion == 0) && (gpu->pbShakeID->_pSysData[j].x == gpu->pbTICommonPair->_pSysData[pos].y));
              bool f1 = ((myRegion == 1) && (gpu->pbShakeID->_pSysData[j].x == gpu->pbTICommonPair->_pSysData[pos].x));
              if (f0 || f1) {
                j_index = j;
                break;
              }
            }
          }
          if (j_index > 0) {  // Merge i and j_index
            picked[j_index] = true;
            unsigned int counter = 0;
            int tt[6] = { -1,-1,-1,-1,-1,-1 };

            double r = ti_weights[myRegion] * gpu->pbShakeParm->_pSysData[i].y
              + ti_weights[1 - myRegion] * gpu->pbShakeParm->_pSysData[j_index].y;
            gpu->pbShakeParm->_pSysData[i].y = r;

            //Re-assign the eq bond lenght to the new SHAKE length
            for (unsigned iBond = 0; iBond < gpu->sim.numberTIBond; iBond++) {

              bool b0 = gpu->pbTIBondID->_pSysData[iBond].x == iatom;
              bool b1 = gpu->pbTIBondID->_pSysData[iBond].y == iatom;

              bool b01 = (gpu->pbTIBondID->_pSysData[iBond].y == iatom0) || (gpu->pbTIBondID->_pSysData[iBond].y == iatom1) || (gpu->pbTIBondID->_pSysData[iBond].y == iatom2);
              bool b11 = (gpu->pbTIBondID->_pSysData[iBond].x == iatom0) || (gpu->pbTIBondID->_pSysData[iBond].x == iatom1) || (gpu->pbTIBondID->_pSysData[iBond].x == iatom2);

              if ( ( b0 &&  b01) || (b01 && b11) )   gpu->pbTIBond->_pSysData[iBond].y = sqrt(r);
            }

            int jatom = gpu->pbShakeID->_pSysData[j_index].x;
            int jatom0 = gpu->pbShakeID->_pSysData[j_index].y;
            int jatom1 = gpu->pbShakeID->_pSysData[j_index].z;
            int jatom2 = gpu->pbShakeID->_pSysData[j_index].w;

            for (unsigned iBond = 0; iBond < gpu->sim.numberTIBond; iBond++) {

              bool b0 = gpu->pbTIBondID->_pSysData[iBond].x == jatom;
              bool b1 = gpu->pbTIBondID->_pSysData[iBond].y == jatom;

              bool b01 = (gpu->pbTIBondID->_pSysData[iBond].y == jatom0) || (gpu->pbTIBondID->_pSysData[iBond].y == jatom1) || (gpu->pbTIBondID->_pSysData[iBond].y == jatom2);
              bool b11 = (gpu->pbTIBondID->_pSysData[iBond].x == jatom0) || (gpu->pbTIBondID->_pSysData[iBond].x == jatom1) || (gpu->pbTIBondID->_pSysData[iBond].x == jatom2);

              if ((b0 && b01) || (b01 && b11))   gpu->pbTIBond->_pSysData[iBond].y = sqrt(r);
            }
            gpu->pbTIBond->Upload();

            if (iatom0 >= 0) tt[counter++] = iatom0;
            if (iatom1 >= 0) tt[counter++] = iatom1;
            if (iatom2 >= 0) tt[counter++] = iatom2;

            if (jatom0 >= 0 && ti_sc_lst[jatom0] > 0) {
              tt[counter++] = jatom0;
              gpu->pbShakeID->_pSysData[j_index].y = -1;
            }
            if (jatom1 >= 0 && ti_sc_lst[jatom1] > 0) {
              tt[counter++] = jatom1;
              gpu->pbShakeID->_pSysData[j_index].z = -1;
            }
            if (jatom2 >= 0 && ti_sc_lst[jatom2] > 0) {
              tt[counter++] = jatom2;
              gpu->pbShakeID->_pSysData[j_index].w = -1;
            }

            if (gpu->pbShakeID->_pSysData[j_index].y < 0 && gpu->pbShakeID->_pSysData[j_index].z >= 0) {
              gpu->pbShakeID->_pSysData[j_index].y = gpu->pbShakeID->_pSysData[j_index].z;
              gpu->pbShakeID->_pSysData[j_index].z = -1;
            }
            if (gpu->pbShakeID->_pSysData[j_index].y < 0 && gpu->pbShakeID->_pSysData[j_index].w >= 0) {
              gpu->pbShakeID->_pSysData[j_index].y = gpu->pbShakeID->_pSysData[j_index].w;
              gpu->pbShakeID->_pSysData[j_index].w = -1;
            }

            if (counter <= 3) {
              gpu->pbShakeID->_pSysData[i].y = tt[0];
              gpu->pbShakeID->_pSysData[i].z = tt[1];
              gpu->pbShakeID->_pSysData[i].w = tt[2];
            } else if (counter==4) {
              unsigned k = gpu->sim.slowShakeConstraints++;
              gpu->sim.slowShakeOffset = (((gpu->sim.slowShakeConstraints + 63) >> 6) << 6) +
                gpu->sim.fastShakeOffset;
              gpu->pbSlowShakeID1->_pSysData[k] = iatom;
              gpu->pbSlowShakeID2->_pSysData[k].x = tt[0];
              gpu->pbSlowShakeID2->_pSysData[k].y = tt[1];
              gpu->pbSlowShakeID2->_pSysData[k].z = tt[2];
              gpu->pbSlowShakeID2->_pSysData[k].w = tt[3];
              gpu->pbSlowShakeParm->_pSysData[k].x = gpu->pbShakeParm->_pSysData[i].x;
              gpu->pbSlowShakeParm->_pSysData[k].y = gpu->pbShakeParm->_pSysData[i].y;
              gpu->pbSlowShakeInvMassH->_pSysData[k] = gpu->pbShakeInvMassH->_pSysData[i];
              gpu->pbShakeID->_pSysData[i].y = -1;
              gpu->pbShakeID->_pSysData[i].z = -1;
              gpu->pbShakeID->_pSysData[i].w = -1;
            }

            // make sure all TC Hydrogen positions will be updated according to the right region
            for (int t = 0; t < 4; t++) {
              for (int k = 0; k < gpu->sim.numberTICommonPairs; k++) {
                if ((gpu->pbTICommonPair->_pSysData[k].x == tt[t]) || (gpu->pbTICommonPair->_pSysData[k].y == tt[t]))
                  gpu->pbTICommonPair->_pSysData[k].z = myRegion;
              }
            }

          }
        } //  if (pos >= 0) {
        picked[i] = true;
      } //  if (ti_sc_lst[iatom] == 0 && (b0 || b1 || b2)) {
    }
  }

  if (gpu->pbTICommonPair) gpu->pbTICommonPair->Upload();

  if (gpu->pbShakeID) gpu->pbShakeID->Upload();
  if (gpu->pbShakeParm) gpu->pbShakeParm->Upload();
  if (gpu->pbSlowShakeID1) gpu->pbSlowShakeID1->Upload();
  if (gpu->pbSlowShakeID2) gpu->pbSlowShakeID2->Upload();
  if (gpu->pbSlowShakeParm) gpu->pbSlowShakeParm->Upload();
  if (gpu->pbSlowShakeInvMassH) gpu->pbSlowShakeInvMassH->Upload();

  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gpu_final_case_inhibitors_: function for checking data previously loaded into the
//                             simulationConst struct.  This is called directly by pmemd.F90
//                             and must be called after all other setup is finished.
//
// Arguments:
//   errmsg:     error message allocated by pmemd.F90 to store errors and post them back to
//               the mdout
//   errlen:     length of the error message (returned)
//   abortsig:   flag to have pmemd abort the run after finding fatal problems (returned)
//---------------------------------------------------------------------------------------------
extern "C" void gpu_final_case_inhibitors_(char errmsg[], int *errlen, int *abortsig,
                                           double atm_crd[][3])
{
  PRINTMETHOD("gpu_final_case_inhibitors");
  errmsg[0] = '\0';

#ifndef use_DPFP
  int i;
  bool errspc = true, fatality = false;

  // Check charges
  const double maxq = 18.0 * 18.2223;
  for (i = 0; i < gpu->sim.atoms; i++) {
    if (fabs(gpu->pbAtomCharge->_pSysData[i]) >= maxq) {
      sprintf(&errmsg[strlen(errmsg)],
              "| ERROR :: Charge of atom %d is %8.4fe, which exceeds the maximum\n"
              "|          tolerance of 18.0e.  This topology must be run using the DPFP\n"
              "|          version of the GPU code, or with the CPU code.\n", i,
              gpu->pbAtomCharge->_pSysData[i] / 18.2223);
      fatality = true;
      if (strlen(errmsg) > 32500) {
        sprintf(&errmsg[strlen(errmsg)], "| Maximum number of warnings reached.\n");
        errspc = false;
      }
    }
  }

  // Check on bond angles
  for (i = 0; i < gpu->sim.bondAngles; i++) {
    double maxang = gpu->pbBondAngle->_pSysData[i].y;
    if (maxang < 0.5*PI) {
      maxang = PI - maxang;
    }
    if (gpu->pbBondAngle->_pSysData[i].x >= 7200.0 / (maxang * maxang) && errspc) {
      sprintf(&errmsg[strlen(errmsg)],
              "| WARNING :: An angle between atoms %d, %d, and %d has an\n"
              "|            extremely high stiffness constant (%10.4f), which could cause\n"
              "|            an overflow in the single-precision force accumulation if this\n"
              "|            angle is stretched very far from equilibrium.\n",
              gpu->pbBondAngleID1->_pSysData[i].x + 1, gpu->pbBondAngleID1->_pSysData[i].y + 1,
              gpu->pbBondAngleID2->_pSysData[i] + 1, gpu->pbBondAngle->_pSysData[i].x);
      if (strlen(errmsg) > 32500) {
        sprintf(&errmsg[strlen(errmsg)], "| Maximum number of warnings reached.\n");
        errspc = false;
      }
    }
  }

  // Check on dihedrals
  for (i = 0; i < gpu->sim.dihedrals; i++) {
    PMEFloat2 dihe1 = gpu->pbDihedral1->_pSysData[i];
    PMEFloat2 dihe2 = gpu->pbDihedral2->_pSysData[i];
    if (dihe2.x * dihe1.y > 7200.0 && errspc) {
      sprintf(&errmsg[strlen(errmsg)],
              "| WARNING :: A dihedral between atoms %d, %d, %d, and %d\n"
              "|            has an extremely high stiffness constant (%10.4f) and\n"
              "|            periodicity (%5.2f).  If this dihedral rotates completely, it\n"
              "|            could cause an overflow in the single-precision force\n"
              "|            accumulation.\n",
              gpu->pbDihedralID1->_pSysData[i].x + 1, gpu->pbDihedralID1->_pSysData[i].y + 1,
              gpu->pbDihedralID1->_pSysData[i].z + 1, gpu->pbDihedralID1->_pSysData[i].w + 1,
              dihe2.x, dihe1.y);
      if (strlen(errmsg) > 32500) {
        sprintf(&errmsg[strlen(errmsg)], "| Maximum number of warnings reached.\n");
        errspc = false;
      }
    }
  }

  // Check on NMR restraints
  for (i = 0; i < gpu->sim.NMRDistances; i++) {
    double2 r1r2 = gpu->pbNMRDistanceR1R2Int->_pSysData[i];
    double2 r3r4 = gpu->pbNMRDistanceR3R4Int->_pSysData[i];
    double2 k2k3 = gpu->pbNMRDistanceK2K3Int->_pSysData[i];
    int2 atmid   = gpu->pbNMRDistanceID->_pSysData[i];
    double dx    = atm_crd[atmid.x][0] - atm_crd[atmid.y][0];
    double dy    = atm_crd[atmid.x][1] - atm_crd[atmid.y][1];
    double dz    = atm_crd[atmid.x][2] - atm_crd[atmid.y][2];
    double r     = sqrt(dx*dx + dy*dy + dz*dz);
    double dl2   = fabs(r - r1r2.y);
    double dl3   = fabs(r - r3r4.x);
    if ((2.0 * k2k3.x * dl2 >= 7200.0 || 2.0 * k2k3.y * dl3 >= 7200.0) && errspc) {
      sprintf(&errmsg[strlen(errmsg)],
              "| WARNING :: The NMR distance restraint between atoms %d and %d could\n"
              "|            be stretched to an unsafe degree.  The run will proceed, but\n"
              "|            overflow in the force accumulation is possible.\n", atmid.x + 1,
              atmid.y + 1);
      if (strlen(errmsg) > 32500) {
        sprintf(&errmsg[strlen(errmsg)], "| Maximum number of warnings reached.\n");
        errspc = false;
      }
    }
  }
  for (i = 0; i < gpu->sim.NMRAngles; i++) {
    double2 k2k3 = gpu->pbNMRAngleK2K3Int->_pSysData[i];
    int2 atmid1  = gpu->pbNMRAngleID1->_pSysData[i];
    int atmid2   = gpu->pbNMRAngleID2->_pSysData[i];
    if ((k2k3.x >= 800.0 || k2k3.y >= 800.0) && errspc) {
      sprintf(&errmsg[strlen(errmsg)],
              "| WARNING :: The NMR angle restraint between atoms %d, %d, and %d\n"
              "|            may be unsafe if stretched very far.  The run will proceed, but\n"
              "|            overflow in the force accumulation is possible.\n", atmid1.x + 1,
              atmid1.y + 1, atmid2 + 1);
      if (strlen(errmsg) > 32500) {
        sprintf(&errmsg[strlen(errmsg)], "| Maximum number of warnings reached.\n");
        errspc = false;
      }
    }
  }
  for (i = 0; i < gpu->sim.NMRTorsions; i++) {
    double2 k2k3 = gpu->pbNMRTorsionK2K3Int->_pSysData[i];
    int4 atmid  = gpu->pbNMRTorsionID1->_pSysData[i];
    if ((k2k3.x >= 800.0 || k2k3.y >= 800.0) && errspc) {
      sprintf(&errmsg[strlen(errmsg)],
              "| WARNING :: The NMR torsion restraint between atoms %d, %d, %d,\n"
              "|            and %d may be unsafe if stretched very far.  The run will\n"
              "|            proceed, but overflow in the force accumulation is possible.\n",
	      atmid.x + 1, atmid.y + 1, atmid.z + 1, atmid.w + 1);
      if (strlen(errmsg) > 32500) {
        sprintf(&errmsg[strlen(errmsg)], "| Maximum number of warnings reached.\n");
        errspc = false;
      }
    }
  }

  if (strlen(errmsg) != 0) {
    printf("%s", errmsg);
  }
  if (fatality) {
    *abortsig = 1;
  }
  else {
    *abortsig = 0;
  }
  *errlen = strlen(errmsg);
#else
  *errlen = 0;
  *abortsig = 0;
#endif
}

//---------------------------------------------------------------------------------------------
// gpu_get_water_distances_: set up fast SHAKE parameters
//
// Arguments:
//   rbtarg:   array of constants needed by SETTLE / RATTLE
//---------------------------------------------------------------------------------------------
extern "C" void gpu_get_water_distances_(double* rbtarg)
{
  PRINTMETHOD("gpu_get_water_distances");
  gpu->sim.ra          = rbtarg[0];
  gpu->sim.ra_inv      = 1.0 / gpu->sim.ra;
  gpu->sim.rb          = rbtarg[1];
  gpu->sim.rc          = rbtarg[2];
  gpu->sim.rc2         = rbtarg[3];
  gpu->sim.hhhh        = rbtarg[4];
  gpu->sim.wo_div_wohh = rbtarg[5] / rbtarg[7];
  gpu->sim.wh_div_wohh = rbtarg[6] / rbtarg[7];
  gpu->sim.wo          = rbtarg[5];
  gpu->sim.wh          = rbtarg[6];
  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gpu_create_outputbuffers_: create buffers for GPU data to land in before it gets processed
//                            by the CPU for printing to file output.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_create_outputbuffers_()
{
  PRINTMETHOD("gpu_create_outputbuffers");

  // Allocate energy buffer
  gpu->sim.EnergyTerms        = EXTENDED_ENERGY_TERMS;
  gpu->pbEnergyBuffer         = new GpuBuffer<unsigned long long int>(gpu->sim.EnergyTerms);
  gpu->sim.pEnergyBuffer      = gpu->pbEnergyBuffer->_pDevData;
  gpu->sim.pEELT              = gpu->sim.pEnergyBuffer + 0;
  gpu->sim.pEVDW              = gpu->sim.pEnergyBuffer + 1;
  gpu->sim.pEGB               = gpu->sim.pEnergyBuffer + 2;
  gpu->sim.pEBond             = gpu->sim.pEnergyBuffer + 3;
  gpu->sim.pEAngle            = gpu->sim.pEnergyBuffer + 4;
  gpu->sim.pEDihedral         = gpu->sim.pEnergyBuffer + 5;
  gpu->sim.pEEL14             = gpu->sim.pEnergyBuffer + 6;
  gpu->sim.pENB14             = gpu->sim.pEnergyBuffer + 7;
  gpu->sim.pEConstraint       = gpu->sim.pEnergyBuffer + 8;
  gpu->sim.pEER               = gpu->sim.pEnergyBuffer + 9;
  gpu->sim.pEED               = gpu->sim.pEnergyBuffer + 10;
  gpu->sim.pEAngle_UB         = gpu->sim.pEnergyBuffer + 11;
  gpu->sim.pEImp              = gpu->sim.pEnergyBuffer + 12;
  gpu->sim.pECmap             = gpu->sim.pEnergyBuffer + 13;
  gpu->sim.pENMRDistance      = gpu->sim.pEnergyBuffer + 14;
  gpu->sim.pENMRCOMDistance   = gpu->sim.pEnergyBuffer + 14;
  gpu->sim.pENMRr6avDistance  = gpu->sim.pEnergyBuffer + 14;
  gpu->sim.pENMRAngle         = gpu->sim.pEnergyBuffer + 15;
  gpu->sim.pENMRCOMAngle      = gpu->sim.pEnergyBuffer + 15;
  gpu->sim.pENMRTorsion       = gpu->sim.pEnergyBuffer + 16;
  gpu->sim.pENMRCOMTorsion    = gpu->sim.pEnergyBuffer + 16;
  gpu->sim.pEEField           = gpu->sim.pEnergyBuffer + 17;
  gpu->sim.pESurf             = gpu->sim.pEnergyBuffer + 18; //pwsasa ESurf
  gpu->sim.pAMDEDihedral      = gpu->sim.pEnergyBuffer + AMD_E_DIHEDRAL_OFFSET;
  gpu->sim.pGaMDEDihedral     = gpu->sim.pEnergyBuffer + GAMD_E_DIHEDRAL_OFFSET;
  gpu->sim.pVirial            = gpu->sim.pEnergyBuffer + VIRIAL_OFFSET;
  gpu->sim.pVirial_11         = gpu->sim.pEnergyBuffer + VIRIAL_OFFSET;
  gpu->sim.pVirial_22         = gpu->sim.pEnergyBuffer + VIRIAL_OFFSET + 1;
  gpu->sim.pVirial_33         = gpu->sim.pEnergyBuffer + VIRIAL_OFFSET + 2;
  gpu->sim.pEKCOMX            = gpu->sim.pEnergyBuffer + VIRIAL_OFFSET + 3;
  gpu->sim.pEKCOMY            = gpu->sim.pEnergyBuffer + VIRIAL_OFFSET + 4;
  gpu->sim.pEKCOMZ            = gpu->sim.pEnergyBuffer + VIRIAL_OFFSET + 5;
  gpu->sim.pVirial_11E        = gpu->sim.pEnergyBuffer + VIRIALE_OFFSET;
  gpu->sim.pVirial_22E        = gpu->sim.pEnergyBuffer + VIRIALE_OFFSET + 1;
  gpu->sim.pVirial_33E        = gpu->sim.pEnergyBuffer + VIRIALE_OFFSET + 2;
  gpu->sim.pEVDWE             = gpu->sim.pEnergyBuffer + EVDWE_OFFSET;
  gpu->sim.pEEDE              = gpu->sim.pEnergyBuffer + EEDE_OFFSET;
  gpu->sim.pEELTE             = gpu->sim.pEnergyBuffer + EELTE_OFFSET;

  if (gpu->sim.ti_mode > 0) {
    gpu->sim.AFETerms       = AFE_TERMS;
    gpu->pbAFEBuffer        = new GpuBuffer<unsigned long long int>(AFE_TERMS);
    gpu->sim.pAFEBuffer     = gpu->pbAFEBuffer->_pDevData;
    gpu->sim.pDVDL          = gpu->sim.pAFEBuffer + 0;
    gpu->sim.pSCBondR1      = gpu->sim.pAFEBuffer + 1;
    gpu->sim.pSCBondR2      = gpu->sim.pAFEBuffer + 2;
    gpu->sim.pSCBondAngleR1 = gpu->sim.pAFEBuffer + 3;
    gpu->sim.pSCBondAngleR2 = gpu->sim.pAFEBuffer + 4;
    gpu->sim.pSCDihedralR1  = gpu->sim.pAFEBuffer + 5;
    gpu->sim.pSCDihedralR2  = gpu->sim.pAFEBuffer + 6;
    gpu->sim.pESCNMRDistanceR1 = gpu->sim.pAFEBuffer + 7;
    gpu->sim.pESCNMRDistanceR2 = gpu->sim.pAFEBuffer + 8;
    gpu->sim.pESCNMRCOMDistanceR1   = gpu->sim.pAFEBuffer + 7;
    gpu->sim.pESCNMRCOMDistanceR2   = gpu->sim.pAFEBuffer + 8;
    gpu->sim.pESCNMRr6avDistanceR1  = gpu->sim.pAFEBuffer + 7;
    gpu->sim.pESCNMRr6avDistanceR2  = gpu->sim.pAFEBuffer + 8;
    gpu->sim.pESCNMRAngleR1    = gpu->sim.pAFEBuffer + 9;
    gpu->sim.pESCNMRAngleR2    = gpu->sim.pAFEBuffer + 10;
    gpu->sim.pESCNMRTorsionR1  = gpu->sim.pAFEBuffer + 11;
    gpu->sim.pESCNMRTorsionR2  = gpu->sim.pAFEBuffer + 12;
    gpu->sim.pSCVDWDirR1    = gpu->sim.pAFEBuffer + 13;
    gpu->sim.pSCVDWDirR2    = gpu->sim.pAFEBuffer + 14;
    gpu->sim.pSCEELDirR1    = gpu->sim.pAFEBuffer + 15;
    gpu->sim.pSCEELDirR2    = gpu->sim.pAFEBuffer + 16;
    gpu->sim.pSCVDW14R1     = gpu->sim.pAFEBuffer + 17;
    gpu->sim.pSCVDW14R2     = gpu->sim.pAFEBuffer + 18;
    gpu->sim.pSCEEL14R1     = gpu->sim.pAFEBuffer + 19;
    gpu->sim.pSCEEL14R2     = gpu->sim.pAFEBuffer + 20;
    gpu->sim.pSCVDWDerR1    = gpu->sim.pAFEBuffer + 21;
    gpu->sim.pSCVDWDerR2    = gpu->sim.pAFEBuffer + 22;
    gpu->sim.pSCEELDerR1    = gpu->sim.pAFEBuffer + 23;
    gpu->sim.pSCEELDerR2    = gpu->sim.pAFEBuffer + 24;
  }

  // Allocate kinetic energy buffer
  if (gpu->bCanMapHostMemory) {
    gpu->pbKineticEnergyBuffer = new GpuBuffer<KineticEnergy>(gpu->blocks, false, true);
    if (gpu->sim.ti_mode > 0) {
      gpu->pbAFEKineticEnergyBuffer = new GpuBuffer<AFEKineticEnergy>(gpu->blocks,
                                                                      false, true);
    }
  }
  else {
    gpu->pbKineticEnergyBuffer = new GpuBuffer<KineticEnergy>(gpu->blocks);
    if (gpu->sim.ti_mode > 0) {
      gpu->pbAFEKineticEnergyBuffer = new GpuBuffer<AFEKineticEnergy>(gpu->blocks);
    }
  }
  gpu->sim.pKineticEnergy = gpu->pbKineticEnergyBuffer->_pDevData;
  if (gpu->sim.ti_mode > 0) {
    gpu->sim.pAFEKineticEnergy = gpu->pbAFEKineticEnergyBuffer->_pDevData;
  }

    
 // Calculate offsets
  gpu->sim.bondOffset = (((gpu->sim.bonds +
                           (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits);
  gpu->sim.bondAngleOffset = gpu->sim.bondOffset +
                             (((gpu->sim.bondAngles + (gpu->sim.grid - 1)) >>
                               gpu->sim.gridBits) << gpu->sim.gridBits);
  gpu->sim.dihedralOffset = gpu->sim.bondAngleOffset +
                            (((gpu->sim.dihedrals + (gpu->sim.grid - 1)) >>
                              gpu->sim.gridBits) << gpu->sim.gridBits);
  gpu->sim.nb14Offset = gpu->sim.dihedralOffset +
                        (((gpu->sim.nb14s + (gpu->sim.grid - 1)) >>
                          gpu->sim.gridBits) << gpu->sim.gridBits);
  gpu->sim.constraintOffset = gpu->sim.nb14Offset +
                              (((gpu->sim.constraints + (gpu->sim.grid - 1)) >>
                                gpu->sim.gridBits) << gpu->sim.gridBits);
  gpu->localForcesBlocks = (gpu->sim.constraintOffset + gpu->localForcesThreadsPerBlock - 1) /
                           gpu->localForcesThreadsPerBlock;

  // Check for zero local interations and skip local force kernel in that case.
  if (gpu->sim.bonds + gpu->sim.bondAngles + gpu->sim.dihedrals + gpu->sim.nb14s +
      gpu->sim.constraints + gpu->sim.UBAngles + gpu->sim.impDihedrals + gpu->sim.cmaps +
      gpu->sim.NMRDistances + gpu->sim.NMRAngles + gpu->sim.NMRTorsions > 0) {
    gpu->bLocalInteractions = true;
  }
  else {
    gpu->bLocalInteractions = false;
  }

  // One more possibility: the force kernel is tasked with initializing the non-bonded
  // reciprocal space charge grid for PME.  Make sure to flag "local interactions" so
  // that this step is not skipped!
  if (gpu->sim.XYZStride > 0) {
    gpu->bLocalInteractions = true;
  }

  // Calculate CHARMM offsets
  gpu->sim.UBAngleOffset = (((gpu->sim.UBAngles +
                              (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits);
  gpu->sim.impDihedralOffset = gpu->sim.UBAngleOffset +
                               (((gpu->sim.impDihedrals + (gpu->sim.grid - 1)) >>
                                 gpu->sim.gridBits) << gpu->sim.gridBits);
  gpu->sim.cmapOffset = gpu->sim.impDihedralOffset +
                        (((gpu->sim.cmaps + (gpu->sim.grid - 1)) >>
                          gpu->sim.gridBits) << gpu->sim.gridBits);
  gpu->CHARMMForcesBlocks = (gpu->sim.cmapOffset + gpu->CHARMMForcesThreadsPerBlock - 1) /
                            gpu->CHARMMForcesThreadsPerBlock;
  if (gpu->sim.UBAngles + gpu->sim.impDihedrals + gpu->sim.cmaps > 0) {
    gpu->bCharmmInteractions = true;
  }
  else {
    gpu->bCharmmInteractions = false;
  }

  // Calculate NMR offsets added COM angle and torsion offset
  gpu->sim.NMRDistanceOffset = (((gpu->sim.NMRDistances + (gpu->sim.grid - 1)) >>
                                 gpu->sim.gridBits) << gpu->sim.gridBits);
  gpu->sim.NMRAngleOffset = gpu->sim.NMRDistanceOffset +
                            (((gpu->sim.NMRAngles + (gpu->sim.grid - 1)) >>
                              gpu->sim.gridBits) << gpu->sim.gridBits);
  gpu->sim.NMRTorsionOffset = gpu->sim.NMRAngleOffset +
                              (((gpu->sim.NMRTorsions + (gpu->sim.grid - 1)) >>
                                gpu->sim.gridBits) << gpu->sim.gridBits);
  gpu->sim.NMRCOMDistanceOffset = gpu->sim.NMRTorsionOffset +
                                  (((gpu->sim.NMRCOMDistances + (gpu->sim.grid - 1)) >>
                                    gpu->sim.gridBits) << gpu->sim.gridBits);
  gpu->sim.NMRr6avDistanceOffset = gpu->sim.NMRCOMDistanceOffset +
                                   (((gpu->sim.NMRr6avDistances + (gpu->sim.grid - 1)) >>
                                     gpu->sim.gridBits) << gpu->sim.gridBits);
  gpu->sim.NMRCOMAngleOffset = gpu->sim.NMRr6avDistanceOffset +
                               (((gpu->sim.NMRCOMAngles + (gpu->sim.grid - 1)) >>
                                     gpu->sim.gridBits) << gpu->sim.gridBits);
  gpu->sim.NMRCOMTorsionOffset = gpu->sim.NMRCOMAngleOffset +
                                 (((gpu->sim.NMRCOMTorsions + (gpu->sim.grid - 1)) >>
                                     gpu->sim.gridBits) << gpu->sim.gridBits);


  gpu->NMRForcesBlocks = (gpu->sim.NMRCOMTorsionOffset + gpu->NMRForcesThreadsPerBlock - 1) /
                         gpu->NMRForcesThreadsPerBlock;

  // Added COMangles and COMtorsions
  if (gpu->sim.NMRCOMDistances + gpu->sim.NMRr6avDistances + gpu->sim.NMRCOMAngles + gpu->sim.NMRCOMTorsions > 0) {
    gpu->bNMRInteractions = true;
  }
  else {
    gpu->bNMRInteractions = false;
  }
  gpuCopyConstants();
  if (!gpu->bNeighborList) {
    kClearGBBuffers(gpu);
  }
}

//---------------------------------------------------------------------------------------------
// gpuCopyConstants: calls a series of functions, depending on the type of simulation, to
//                   establish critical constants on the GPU.  This is invoked in many places!
//---------------------------------------------------------------------------------------------
extern "C" void gpuCopyConstants()
{
  PRINTMETHOD("gpuCopyConstants");
#ifdef GTI
  gpu->UpdateSimulationConst();
#else
  SetkForcesUpdateSim(gpu);
  SetkDataTransferSim(gpu);
  SetkCalculateNEBForcesSim(gpu);
  SetkShakeSim(gpu);
  SetkCalculateLocalForcesSim(gpu);
  SetkXrayUpdateSim(gpu);
  if (gpu->ntb == 0) {
    SetkCalculateGBBornRadiiSim(gpu);
    SetkCalculateGBNonbondEnergy1Sim(gpu);
    SetkCalculateGBNonbondEnergy2Sim(gpu);
  }
  else {
    if (gpu->sim.icnstph == 2 || gpu->sim.icnste == 2 || gpu->sim.hybridgb > 0) {
      SetkCalculateGBBornRadiiSim(gpu);
      SetkCalculateGBNonbondEnergy1Sim(gpu);
      SetkCalculateGBNonbondEnergy2Sim(gpu);
    }
    SetkPMEInterpolationSim(gpu);
    SetkNeighborListSim(gpu);
    SetkCalculateEFieldEnergySim(gpu);
    SetkCalculatePMENonbondEnergySim(gpu);
  }
#endif
}

//---------------------------------------------------------------------------------------------
// gpu_init_pbc_: initialize periodic boundary conditions on the GPU
//
// Arguments:
//   a,b,c:      the three lengths of the simulation box
//   alpha:
//   beta:       box angles for the simulation box
//   gamma:
//   uc_volume:  volume of the simulation unit cell
//   uc_sphere:  largest sphere that will fit within the unit cell
//   max_cutoff: the maximum cutoff that could be imposed based on uc_sphere and the
//               non-bonded skin (buffer region for diffusion before pairlist rebuilds)
//   pbc_box:
//   reclng:
//   cut_factor:
//   ucell:      transformation matrix to take fractional coordinates into real space
//   recip:      transformation matrix to take Cartesian coordinates into fractional
//---------------------------------------------------------------------------------------------
extern "C" void gpu_init_pbc_(double* a, double* b, double* c, double* alpha, double* beta,
                              double* gamma, double* uc_volume, double* uc_sphere,
                              double* max_cutoff, double pbc_box[], double reclng[],
                              double cut_factor[], double ucell[][3], double recip[][3])
{
  PRINTMETHOD("gpu_init_pbc");
  gpu->bNeighborList = true;
  gpu->bNeedNewNeighborList = true;
  gpu->bNewNeighborList = false;

  // Check for orthogonal system
  if ((*alpha == 90.0) && (*beta == 90.0) && (*gamma == 90.0)) {
    gpu->sim.is_orthog = true;
  }
  else {
    gpu->sim.is_orthog = false;
  }
  gpu->sim.af    = (PMEFloat)*a;
  gpu->sim.bf    = (PMEFloat)*b;
  gpu->sim.cf    = (PMEFloat)*c;
  gpu->sim.a     = (PMEDouble)gpu->sim.af;
  gpu->sim.b     = (PMEDouble)gpu->sim.bf;
  gpu->sim.c     = (PMEDouble)gpu->sim.cf;
  gpu->sim.alpha = *alpha * PI / 180.0;
  gpu->sim.beta  = *beta * PI / 180.0;
  gpu->sim.gamma = *gamma * PI / 180.0;

  for (int i = 0; i < 3; i++) {
    gpu->sim.pbc_box[i] = pbc_box[i];
    gpu->sim.reclng[i] = reclng[i];
    gpu->sim.cut_factor[i] = cut_factor[i];
    for (int j = 0; j < 3; j++) {
      gpu->sim.ucellf[i][j] = ucell[j][i];
      gpu->sim.ucell[i][j] = gpu->sim.ucellf[i][j];
      gpu->sim.recipf[i][j] = recip[j][i];
      gpu->sim.recip[i][j] = gpu->sim.recipf[i][j];
    }
  }
  gpu->sim.uc_volume = *uc_volume;
  gpu->sim.uc_sphere = *uc_sphere;
  gpu->sim.pi_vol_inv = 1.0 / (PI * gpu->sim.uc_volume);

  gpuCopyConstants();
  return;
}

//---------------------------------------------------------------------------------------------
// cellOffset:
//---------------------------------------------------------------------------------------------
static const int cellOffset[][3] = {
  { 0,  0,  0},
  { 1,  0,  0},
  {-1,  1,  0},
  { 0,  1,  0},
  { 1,  1,  0},
  {-1, -1,  1},
  { 0, -1,  1},
  { 1, -1,  1},
  {-1,  0,  1},
  { 0,  0,  1},
  { 1,  0,  1},
  {-1,  1,  1},
  { 0,  1,  1},
  { 1,  1,  1}
};


//---------------------------------------------------------------------------------------------
// cellHash: so, this is the hash table--a static 4 x 4 x 4 cube of cells following a Hilbert
//           space-filling curve.
//---------------------------------------------------------------------------------------------
static const unsigned int cellHash[CELL_HASH_CELLS] = {
  34, 35, 60, 61,
  33, 32, 63, 62,
  30, 31,  0,  1,
  29, 28,  3,  2,

  37, 36, 59, 58,
  38, 39, 56, 57,
  25, 24,  7,  6,
  26, 27,  4,  5,

  42, 43, 54, 53,
  41, 40, 55, 52,
  20, 23,  8,  9,
  21, 22, 11, 10,

  45, 44, 49, 50,
  46, 47, 48, 51,
  19, 16, 15, 14,
  18, 17, 12, 13,
};

//---------------------------------------------------------------------------------------------
// gpu_neighbor_list_setup_: this function manages neighbor list setup, called by
//                           final_pme_setup once during a periodic simulation.
//
// Arguments:
//   numex:        the number of non-bonded exclusions applicable to each atom in the system
//   natex:        indices of excluded atoms
//   vdw_cutoff:   the Lennard-Jones cutoff (this is the larger of two cutoffs in the CPU code,
//                 and the GPU code will just take the electrostatic cutoff to be equal to the
//                 vdW cutoff)
//   skinnb:       the non-bonded skin for neighbor list production
//   skin_permit:  the permitted factor of skinnb which particles may travel before triggering
//                 a pair list rebuild
//---------------------------------------------------------------------------------------------
extern "C" void gpu_neighbor_list_setup_(int numex[], int natex[], double *vdw_cutoff,
                                         double *skinnb, double *skin_permit, int *smbx_permit,
                                         int *ti_shake)
{
  PRINTMETHOD("gpu_neighbor_list_setup");

  int i, j;

  // Delete any existing neighbor list data
  delete gpu->pbImageIndex;
  delete gpu->pbSubImageLookup;
  delete gpu->pbAtomXYSaveSP;
  delete gpu->pbAtomZSaveSP;
  delete gpu->pbImage;
  delete gpu->pbImageVel;
  delete gpu->pbImageLVel;
  delete gpu->pbImageMass;
  delete gpu->pbImageCharge;
  delete gpu->pbImageSigEps;
  delete gpu->pbImageLJID;
  delete gpu->pbImageCellID;
  delete gpu->pbImageTIRegion;
  delete gpu->pbImageTILinearAtmID;
  delete gpu->pbNLCellHash;
  delete gpu->pbNLNonbondCellStartEnd;
  delete gpu->pbBNLExclusionBuffer;
  delete gpu->pbNLExclusionList;
  delete gpu->pbNLExclusionStartCount;
  delete gpu->pbNLAtomList;
  delete gpu->pbNLTotalOffset;
  delete gpu->pbNLEntries;
  delete gpu->pbFinishedBlocksCounters;
  delete gpu->pbKernelCompletionFlags;
  delete gpu->pbNLbSkinTestFail;
  delete gpu->pbImageShakeID;
  delete gpu->pbImageFastShakeID;
  delete gpu->pbImageSlowShakeID1;
  delete gpu->pbImageSlowShakeID2;
  delete gpu->pbImageSlowTIShakeID1;
  delete gpu->pbImageSlowTIShakeID2;
  delete gpu->pbImageSlowTIShakeID3;
  delete gpu->pbImageSoluteAtomID;
  delete gpu->pbImageSolventAtomID;
  delete gpu->pbImageExtraPoint11Frame;
  delete gpu->pbImageExtraPoint11Index;
  delete gpu->pbImageExtraPoint12Frame;
  delete gpu->pbImageExtraPoint12Index;
  delete gpu->pbImageExtraPoint21Frame;
  delete gpu->pbImageExtraPoint21Index;
  delete gpu->pbImageExtraPoint22Frame;
  delete gpu->pbImageExtraPoint22Index;
#ifdef MPI
  delete[] gpu->pMinLocalCell;
  delete[] gpu->pMaxLocalCell;
  delete[] gpu->pMinLocalAtom;
  delete[] gpu->pMaxLocalAtom;
#endif

  // Copy constants
  gpu->sim.cut                    = *vdw_cutoff;
  gpu->sim.cut2                   = (*vdw_cutoff) * (*vdw_cutoff);
  gpu->sim.cut3                   = gpu->sim.cut * gpu->sim.cut2;
  gpu->sim.cut6                   = gpu->sim.cut3 * gpu->sim.cut3;
  gpu->sim.cutinv                 = 1/(*vdw_cutoff);
  gpu->sim.cut2inv                = gpu->sim.cutinv * gpu->sim.cutinv;
  gpu->sim.cut3inv                = gpu->sim.cut2inv * gpu->sim.cutinv;
  gpu->sim.cut6inv                = gpu->sim.cut3inv * gpu->sim.cut3inv;
  gpu->sim.invfswitch6cut6        = 1 / (gpu->sim.fswitch6 * gpu->sim.cut6);
  gpu->sim.invfswitch3cut3        = 1 / (gpu->sim.fswitch3 * gpu->sim.cut3);
  gpu->sim.cut6invcut6minfswitch6 = gpu->sim.cut6 / (gpu->sim.cut6 - gpu->sim.fswitch6);
  gpu->sim.cut3invcut3minfswitch3 = gpu->sim.cut3 / (gpu->sim.cut3 - gpu->sim.fswitch3);
  gpu->sim.skinnb                 = *skinnb;
  gpu->sim.skinPermit             = *skin_permit;
  gpu->sim.cutPlusSkin            = *vdw_cutoff + *skinnb;
  gpu->sim.cutPlusSkin2           = gpu->sim.cutPlusSkin * gpu->sim.cutPlusSkin;

  gpu->pbFinishedBlocksCounters = new GpuBuffer<unsigned int>(1);
  if (gpu->bCanMapHostMemory) {
    gpu->pbKernelCompletionFlags = new GpuBuffer<unsigned int>(1, false, true);
  }
  else {
    gpu->pbKernelCompletionFlags = new GpuBuffer<unsigned int>(1);
  }

  // Allocate skin test buffer
  if (gpu->bCanMapHostMemory) {
    gpu->pbNLbSkinTestFail = new GpuBuffer<bool>(1, false, true);
  }
  else {
    gpu->pbNLbSkinTestFail = new GpuBuffer<bool>(1);
  }

  // Allocate cell hash
  gpu->pbNLCellHash = new GpuBuffer<unsigned int>(CELL_HASH_CELLS);

  // Fill cell hash array
  for (i = 0; i < CELL_HASH_CELLS; i++) {
    gpu->pbNLCellHash->_pSysData[i] = cellHash[i];
  }
  gpu->pbNLCellHash->Upload();

  // Allocate new exclusion lists
  gpu->pbNLExclusionStartCount = new GpuBuffer<uint2>(gpu->sim.atoms);

  // Double count exclusions
  unsigned int* pExclusionCount = new unsigned int[gpu->sim.atoms];
  memset(pExclusionCount, 0, gpu->sim.atoms * sizeof(unsigned int));
  int offset = 0;
  for (i = 0; i < gpu->sim.atoms; i++) {
    for (j = 0; j < numex[i]; j++) {
      if (natex[offset + j] > 0) {
        pExclusionCount[i] += 1;
        pExclusionCount[natex[offset + j] - 1] += 1;
      }
    }
    offset += numex[i];
  }

  // Count total exclusions
  int exclusions = 0;
  for (i = 0; i < gpu->sim.atoms; i++) {
    gpu->pbNLExclusionStartCount->_pSysData[i].x = exclusions;
    gpu->pbNLExclusionStartCount->_pSysData[i].y = 0;
    exclusions += ((pExclusionCount[i] + (GRID - 1)) >> GRID_BITS) << GRID_BITS;
  }
  gpu->pbNLExclusionList = new GpuBuffer<unsigned int>(exclusions);
  memset(gpu->pbNLExclusionList->_pSysData, 0xff, exclusions * sizeof(unsigned int));
  offset = 0;
  for (i = 0; i < gpu->sim.atoms; i++) {
    int pos_i = gpu->pbNLExclusionStartCount->_pSysData[i].x +
                gpu->pbNLExclusionStartCount->_pSysData[i].y;
    for (j = 0; j < numex[i]; j++) {
      if (natex[offset + j] > 0) {
        int pos_j = gpu->pbNLExclusionStartCount->_pSysData[natex[offset + j] - 1].x +
                    gpu->pbNLExclusionStartCount->_pSysData[natex[offset + j] - 1].y;
        gpu->pbNLExclusionList->_pSysData[pos_i] = natex[offset + j] - 1;
        gpu->pbNLExclusionList->_pSysData[pos_j] = i;
        pos_i++;
        gpu->pbNLExclusionStartCount->_pSysData[i].y++;
        gpu->pbNLExclusionStartCount->_pSysData[natex[offset + j] - 1].y++;
      }
    }
    offset += numex[i];
  }

  // Sort exclusion count for determining maximum observable exclusions per warp
  std::sort(pExclusionCount, pExclusionCount + gpu->sim.atoms);
  gpu->pbNLExclusionStartCount->Upload();
  gpu->pbNLExclusionList->Upload();
  gpu->sim.pNLExclusionList       = gpu->pbNLExclusionList->_pDevData;
  gpu->sim.pNLExclusionStartCount = gpu->pbNLExclusionStartCount->_pDevData;

  // Use Hessian normal form to compute the distance between box faces
  double udata[9], cdepth[3];
  udata[0] = gpu->sim.ucell[0][0];
  udata[1] = gpu->sim.ucell[0][1];
  udata[2] = gpu->sim.ucell[0][2];
  udata[3] = gpu->sim.ucell[1][0];
  udata[4] = gpu->sim.ucell[1][1];
  udata[5] = gpu->sim.ucell[1][2];
  udata[6] = gpu->sim.ucell[2][0];
  udata[7] = gpu->sim.ucell[2][1];
  udata[8] = gpu->sim.ucell[2][2];
  HessianNorms(udata, cdepth);

  // Calculate nonbond cell size
  PMEFloat cell = gpu->sim.cut + gpu->sim.skinnb;
  gpu->sim.xcells = int(cdepth[0] / cell);
  if (gpu->sim.xcells < 1) {
    gpu->sim.xcells = 1;
  }
  gpu->sim.ycells = int(cdepth[1] / cell);
  if (gpu->sim.ycells < 1) {
    gpu->sim.ycells = 1;
  }
  gpu->sim.zcells = int(cdepth[2] / cell);
  if (gpu->sim.zcells < 1) {
    gpu->sim.zcells = 1;
  }
  gpu->sim.xcell = (PMEFloat)(gpu->sim.a / (PMEDouble)gpu->sim.xcells);
  gpu->sim.ycell = (PMEFloat)(gpu->sim.b / (PMEDouble)gpu->sim.ycells);
  gpu->sim.zcell = (PMEFloat)(gpu->sim.c / (PMEDouble)gpu->sim.zcells);
  gpu->sim.minCellX = -0.5 / gpu->sim.xcells;
  gpu->sim.minCellY = -0.5 / gpu->sim.ycells;
  gpu->sim.minCellZ = -0.5 / gpu->sim.zcells;
  gpu->sim.maxCellX =  1.5 / gpu->sim.xcells;
  gpu->sim.maxCellY =  1.5 / gpu->sim.ycells;
  gpu->sim.maxCellZ =  1.5 / gpu->sim.zcells;
  gpu->sim.oneOverXcellsf = (PMEFloat)1.0 / (PMEFloat)(gpu->sim.xcells);
  gpu->sim.oneOverYcellsf = (PMEFloat)1.0 / (PMEFloat)(gpu->sim.ycells);
  gpu->sim.oneOverZcellsf = (PMEFloat)1.0 / (PMEFloat)(gpu->sim.zcells);
  gpu->sim.oneOverXcells = gpu->sim.oneOverXcellsf;
  gpu->sim.oneOverYcells = gpu->sim.oneOverYcellsf;
  gpu->sim.oneOverZcells = gpu->sim.oneOverZcellsf;

  // Test for small nonbond cell count on any dimension
  if ((gpu->sim.xcells <= 2) || (gpu->sim.ycells <= 2) || (gpu->sim.zcells <= 2)) {
    gpu->bSmallBox = true;
#ifdef MPI
    if (gpu->gpuID == 0 && *smbx_permit == 0) {
#else
    if (*smbx_permit == 0) {
#endif
      printf("gpu_neighbor_list_setup :: Small box detected, with <= 2 cells in one or more\n"
             "                           dimensions.  The current GPU code has been deemed\n"
             "                           unsafe for these situations.  Please alter the\n"
             "                           cutoff to increase the number of hash cells, make\n"
             "                           use of the CPU code, or (if absolutely necessary)\n"
             "                           run pmemd.cuda with the -AllowSmallBox flag.  This\n"
             "                           behavior will be corrected in a forthcoming "
             "patch.\n");
      exit(-1);
#ifdef MPI
    }
#else
    }
#endif
  }
  else {
    gpu->bSmallBox = false;
  }
  gpu->sim.cell = gpu->sim.xcell;
  if (gpu->sim.ycell > gpu->sim.cell) {
    gpu->sim.cell = gpu->sim.ycell;
  }
  if (gpu->sim.zcell > gpu->sim.cell) {
    gpu->sim.cell = gpu->sim.zcell;
  }
  double xskin = cdepth[0]/gpu->sim.xcells - gpu->sim.cut;
  double yskin = cdepth[1]/gpu->sim.ycells - gpu->sim.cut;
  double zskin = cdepth[2]/gpu->sim.zcells - gpu->sim.cut;
  gpu->sim.nonbond_skin = xskin;
  if (yskin < gpu->sim.nonbond_skin) {
    gpu->sim.nonbond_skin = yskin;
  }
  if (zskin < gpu->sim.nonbond_skin) {
    gpu->sim.nonbond_skin = zskin;
  }
  gpu->sim.one_half_nonbond_skin_squared  = (gpu->sim.nonbond_skin * (*skin_permit));
  gpu->sim.one_half_nonbond_skin_squared *= gpu->sim.one_half_nonbond_skin_squared;
  gpu->sim.cutPlusSkin = gpu->sim.cut + gpu->sim.nonbond_skin;
  gpu->sim.cutPlusSkin2 = gpu->sim.cutPlusSkin * gpu->sim.cutPlusSkin;
  gpu->sim.cells = gpu->sim.xcells * gpu->sim.ycells * gpu->sim.zcells;
  gpu->pbNLNonbondCellStartEnd = new GpuBuffer<uint2>(gpu->sim.cells);
  gpu->neighborListBits = int(log((double)gpu->sim.cells)/log((double)2.0) + 1 +
                              CELL_HASH_BITS);

  // Calculate local cell offsets
  for (i = 0; i < NEIGHBOR_CELLS; i++) {
    if ((gpu->sim.is_orthog) &&
        ((gpu->sim.ntp == 0) || ((gpu->sim.ntp > 0) && (gpu->sim.barostat != 1)))) {
      gpu->sim.cellOffset[i][0] = (PMEDouble)cellOffset[i][0] * gpu->sim.xcell;
      gpu->sim.cellOffset[i][1] = (PMEDouble)cellOffset[i][1] * gpu->sim.ycell;
      gpu->sim.cellOffset[i][2] = (PMEDouble)cellOffset[i][2] * gpu->sim.zcell;
    }
    else {
      gpu->sim.cellOffset[i][0] = (PMEDouble)cellOffset[i][0] / (PMEDouble)gpu->sim.xcells;
      gpu->sim.cellOffset[i][1] = (PMEDouble)cellOffset[i][1] / (PMEDouble)gpu->sim.ycells;
      gpu->sim.cellOffset[i][2] = (PMEDouble)cellOffset[i][2] / (PMEDouble)gpu->sim.zcells;
    }
  }
#ifdef MPI
  bool** pCellMap = new bool*[gpu->nGpus];
  for (i = 0; i < gpu->nGpus; i++) {
    pCellMap[i] = new bool[gpu->sim.cells];
    for (j = 0; j < gpu->sim.cells; j++) {
      pCellMap[i][j] = false;
    }
  }

  // Determine work allocation
  gpu->bCalculateLocalForces      = false;
  gpu->bCalculateDirectSum        = false;
  gpu->bCalculateReciprocalSum    = false;
  if (gpu->gpuID == 0) {
    gpu->bCalculateLocalForces  = true;
  }
  if (gpu->ips == 0) {
    if ((gpu->gpuID == 0) && (gpu->nGpus < 5)) {
        gpu->bCalculateReciprocalSum = true;
    }
    else if ((gpu->gpuID == 1) && (gpu->nGpus >= 5)) {
      gpu->bCalculateReciprocalSum = true;
    }
    if (gpu->nGpus < 3) {
      gpu->bCalculateDirectSum = true;
#if defined(AMBER_PLATFORM_AMD) && defined(use_SPFP)
      if (gpu->nGpus == 2 && gpu->gpuID == 0) {
        gpu->bCalculateDirectSum = false;
      }
#endif
    }
    else if ((gpu->nGpus < 5) && (gpu->gpuID >= 1)) {
      gpu->bCalculateDirectSum = true;
    }
    else if ((gpu->nGpus >= 5) && (gpu->gpuID >= 2)) {
      gpu->bCalculateDirectSum = true;
    }
  }
  else {
    gpu->bCalculateDirectSum    = true;
  }

  // Determine global cell allocation
  gpu->pMinLocalCell = new int[gpu->nGpus];
  gpu->pMaxLocalCell = new int[gpu->nGpus];
  gpu->pMinLocalAtom = new int[gpu->nGpus];
  gpu->pMaxLocalAtom = new int[gpu->nGpus];

  // Now allocate cells
  for (i = 0; i < gpu->nGpus; i++) {
    int minCell = 0;
    int maxCell = 0;
    if (gpu->sim.bIPSActive) {
      minCell = (i * gpu->sim.cells) / gpu->nGpus;
      maxCell = ((i + 1) * gpu->sim.cells) / gpu->nGpus;
    }
    else {
      if (gpu->nGpus == 1) {
        maxCell = gpu->sim.cells;
      }
#if defined(AMBER_PLATFORM_AMD) && defined(use_SPFP)
      else if (false) {
#else
      else if (gpu->nGpus == 2) {
#endif

        // Load-balance direct sum for this case only
        unsigned int pmeweight = 2200;
        unsigned int totalweight = 10000;

        // Adjust for larger systems
        if (gpu->sim.atoms > 150000) {
          pmeweight = 3150;
        }
        else if (gpu->sim.atoms > 24000) {
          pmeweight = 2200 + (gpu->sim.atoms - 24000) / 133;
        }

        // Adjust for FFT volume normalization
        float volume    = gpu->sim.a * gpu->sim.b * gpu->sim.c;
        float fftvolume = gpu->sim.nfft1 * gpu->sim.nfft2 * gpu->sim.nfft3;
        float reference = 0.92f;
        float actual    = volume / fftvolume;
        float relative  = actual / reference;
        pmeweight *= pow((actual / reference), 0.4f);

        // Adjust for nonbond cutoff
        if (gpu->gpuID == 0) {
          maxCell = (pmeweight * gpu->sim.cells) / totalweight;
        }
        else {
          minCell = (pmeweight * gpu->sim.cells) / totalweight;
          maxCell = gpu->sim.cells;
        }
      }
      else if ((gpu->nGpus < 5) && (i > 0)) {
        minCell = (i - 1) * gpu->sim.cells / (gpu->nGpus - 1);
        maxCell = i * gpu->sim.cells / (gpu->nGpus - 1);
      }
      else if ((gpu->nGpus >= 5) && (i > 1)) {
        minCell = (i - 2) * gpu->sim.cells / (gpu->nGpus - 2);
        maxCell = (i - 1) * gpu->sim.cells / (gpu->nGpus - 2);
      }
    }
    for (j = minCell; j < maxCell; j++) {
      pCellMap[i][j]          = true;
    }
    gpu->pMinLocalCell[i]       = minCell;
    gpu->pMaxLocalCell[i]       = maxCell;
    gpu->pMinLocalAtom[i]       = 0;
    gpu->pMaxLocalAtom[i]       = 0;
  }

  // Determine local cell allocation
  gpu->minLocalCell = gpu->pMinLocalCell[gpu->gpuID];
  gpu->maxLocalCell = gpu->pMaxLocalCell[gpu->gpuID];
#  ifdef GVERBOSE
  printf("Node %d %d %d\n", gpu->gpuID, gpu->minLocalCell, gpu->maxLocalCell);
#  endif
#else
  // Determine work allocation
  gpu->bCalculateLocalForces      = true;
#endif // MPI

  // Build work unit list
  int workUnits           = 0;
  NLRecord* pNLCellRecord = new NLRecord[gpu->sim.cells];
  unsigned int homeCells  = 0;
  int offUnits = 0;
  for (i = 0; i < gpu->sim.cells; i++) {
    int x = i % gpu->sim.xcells;
    int y = ((i - x) / gpu->sim.xcells) % gpu->sim.ycells;
    int z = (i - x - y * gpu->sim.xcells) / (gpu->sim.xcells * gpu->sim.ycells);

    // Clear neighbor list data
    pNLCellRecord[i].NL.homeCell = i;
    pNLCellRecord[i].NL.neighborCells = 0;

    for (j = 0; j < NEIGHBOR_CELLS; j++) {
      int x1 = x + cellOffset[j][0];
      int y1 = y + cellOffset[j][1];
      int z1 = z + cellOffset[j][2];
      if (x1 < 0) {
        x1 += gpu->sim.xcells;
      }
      if (y1 < 0) {
        y1 += gpu->sim.ycells;
      }
      if (z1 < 0) {
        z1 += gpu->sim.zcells;
      }
      if (x1 >= gpu->sim.xcells) {
        x1 -= gpu->sim.xcells;
      }
      if (y1 >= gpu->sim.ycells) {
        y1 -= gpu->sim.ycells;
      }
      if (z1 >= gpu->sim.zcells) {
        z1 -= gpu->sim.zcells;
      }
      int cell1 = (z1 * gpu->sim.ycells + y1) * gpu->sim.xcells + x1;

#ifdef MPI
      // Add to local work list if central cell is controlled by node
      if ((i >= gpu->minLocalCell) && (i < gpu->maxLocalCell)) {
#endif
      workUnits++;

      // Add to neighbor list records
      if (pNLCellRecord[i].NL.neighborCells == 0) {
        homeCells++;
      }
      pNLCellRecord[i].NL.neighborCell[pNLCellRecord[i].NL.neighborCells++] =
        (cell1 << NLRECORD_CELL_SHIFT) | j;
#ifdef MPI
      }
#endif
    }
  }

  // Determine cell divisors for neighbor list generation
  unsigned int totalEnergyWarps;
  unsigned int totalForcesWarps;
  if (gpu->ips == 0) {
    totalEnergyWarps = (gpu->PMENonbondEnergyThreadsPerBlock * gpu->PMENonbondBlocks) / GRID;
    totalForcesWarps = (gpu->PMENonbondForcesThreadsPerBlock * gpu->PMENonbondBlocks) / GRID;
  }
  else {
    totalEnergyWarps = (gpu->IPSNonbondEnergyThreadsPerBlock * gpu->IPSNonbondBlocks) / GRID;
    totalForcesWarps = (gpu->IPSNonbondForcesThreadsPerBlock * gpu->IPSNonbondBlocks) / GRID;
  }
  unsigned int yDivisor     = 1;
  unsigned int xDivisor     = 1;
  unsigned int xEntryWidth  = 128;
  unsigned int atomsPerCell = (gpu->sim.atoms + GRID) / gpu->sim.cells;
  unsigned int atomsPerWarp = 32;

  // Adjust for crazy vacuum-filled solvent boxes, can waste memory like crazy if we hit this,
  // but we *shouldn't* be hitting this, period.
  unsigned int minAtomsPerCell = 0.09 * gpu->sim.xcell * gpu->sim.ycell * gpu->sim.zcell;
  if (atomsPerCell < minAtomsPerCell) {
    atomsPerCell = minAtomsPerCell;
  }
#ifdef GVERBOSE
  printf("%d home cells, %d and %d total warps\n", homeCells, totalEnergyWarps,
         totalForcesWarps);
#endif

  // Adjust Y Divisor
  if (homeCells < 8) {
    yDivisor = 8;
  }
  else if (homeCells < totalForcesWarps * 2) {
    yDivisor = 4;
  }
  else {
    yDivisor = 3;
  }

  // Adjust X divisor further for really really small cell counts and high multi-GPU node
  // counts
  if (homeCells < 2 * totalForcesWarps / 3) {
    xDivisor = 4;
  }
  else if (homeCells < 5 * totalForcesWarps / 2) {
    xDivisor = 3;
  }
  else if (homeCells < 3 * totalForcesWarps) {
    xDivisor = 2;
  }
  else {
    xDivisor = 1;
  }

  // Determine maximum possible observable exclusions
  int maxExclusionsPerWarp8 = 0;
  int maxExclusionsPerWarp16 = 0;
  int maxExclusionsPerWarp32 = 0;
  int maxExclusionsPerWarp = 0;

  for (int i = gpu->sim.atoms - 1; i >= max(0, gpu->sim.atoms - 8); i--) {
    maxExclusionsPerWarp8 += pExclusionCount[i];
  }
  for (int i = gpu->sim.atoms - 1; i >= max(0, gpu->sim.atoms - 16); i--) {
    maxExclusionsPerWarp16 += pExclusionCount[i];
  }
  for (int i = gpu->sim.atoms - 1; i >= max(0, gpu->sim.atoms - 32); i--) {
    maxExclusionsPerWarp32 += pExclusionCount[i];
  }
#ifdef GVERBOSE
  printf("XM %d %d %d\n", maxExclusionsPerWarp8, maxExclusionsPerWarp16,
         maxExclusionsPerWarp32);
#endif

  // Adjust warp atoms based on cell volume
  if (atomsPerCell < 256) {
    atomsPerWarp                = 16;
    maxExclusionsPerWarp        = maxExclusionsPerWarp16;
  }
  else {
    maxExclusionsPerWarp        = maxExclusionsPerWarp32;
  }

  // Determine y divisor for building neighbor lists
  yDivisor = (atomsPerCell + 2 * atomsPerWarp - 1) / atomsPerWarp;

  // Manual overrides
#ifdef AMBER_PLATFORM_AMD
  xEntryWidth          = 256;
#else
  xEntryWidth          = 128;
#endif
  yDivisor             = 8;
  xDivisor             = 1;
  atomsPerWarp         = 16;
  maxExclusionsPerWarp = maxExclusionsPerWarp16;
  unsigned int NLMaxEntries = 4 * gpu->sim.cells *
                              ((atomsPerCell + atomsPerWarp) / atomsPerWarp) * NEIGHBOR_CELLS *
                              ((atomsPerCell + GRID) / GRID);

  // Determine exclusion buffer size
  if (maxExclusionsPerWarp <= 256) {
    gpu->sim.NLExclusionBufferSize = 256;
  }
  else {
    gpu->sim.NLExclusionBufferSize = 256 + ((maxExclusionsPerWarp - 256 + 511) / 512) * 512;
  }

#ifdef GVERBOSE
  printf("APW MEX EBS YD XD ME %d %d %d %d %d %d\n", atomsPerWarp, maxExclusionsPerWarp,
         gpu->sim.NLExclusionBufferSize, yDivisor, xDivisor, NLMaxEntries);
#endif

  // Count neighbor list records
  unsigned int NLRecords = gpu->sim.cells * yDivisor;

#ifdef GVERBOSE
  printf("%d neighbor list records %d %d %d %d %d\n", NLRecords, gpu->sim.xcells,
         gpu->sim.ycells, gpu->sim.zcells, yDivisor, xDivisor);
#endif
  if (atomsPerWarp == 32) {
    gpu->sim.NLBuildWarps = (gpu->BNLBlocks *
                             gpu->NLBuildNeighborList32ThreadsPerBlock) / GRID;
  }
  else if (atomsPerWarp == 16) {
    gpu->sim.NLBuildWarps = (gpu->BNLBlocks *
                             gpu->NLBuildNeighborList16ThreadsPerBlock) / GRID;
  }
  else {
    gpu->sim.NLBuildWarps = (gpu->BNLBlocks *
                             gpu->NLBuildNeighborList8ThreadsPerBlock) / GRID;
  }

  // Allocate buffers for building neighbor list
  gpu->pbBNLExclusionBuffer       = new GpuBuffer<uint>(gpu->sim.NLExclusionBufferSize *
                                                        gpu->sim.NLBuildWarps);
  gpu->sim.pBNLExclusionBuffer    = gpu->pbBNLExclusionBuffer->_pDevData;
  gpu->sim.NLNonbondEnergyWarps   = totalEnergyWarps;
  gpu->sim.NLNonbondForcesWarps   = totalForcesWarps;
  gpu->sim.NLMaxEntries           = NLMaxEntries;
  gpu->pbNLRecord                 = new GpuBuffer<NLRecord>(NLRecords);
  gpu->pbNLEntry                  = new GpuBuffer<NLEntry>(NLMaxEntries);
  gpu->sim.NLYDivisor             = yDivisor;
  gpu->sim.NLXEntryWidth          = xEntryWidth;

  // Placeholder, not redundant upon the return of xDivisor
  gpu->sim.NLEntryTypes           = yDivisor;
  gpu->sim.NLYStride              = yDivisor * atomsPerWarp;
  gpu->sim.NLAtomsPerWarp         = atomsPerWarp;
  gpu->sim.NLAtomsPerWarpBits     = ffs(atomsPerWarp) - 1;
  gpu->sim.NLAtomsPerWarpBitsMask = (1 << gpu->sim.NLAtomsPerWarpBits) - 1;
  gpu->sim.NLAtomsPerWarpMask     =
    (unsigned long)((1ull << (1ull << (gpu->sim.NLAtomsPerWarpBits))) - 1ull);
  NLRecord* pNLRecord             = gpu->pbNLRecord->_pSysData;
  gpu->sim.pNLRecord              = gpu->pbNLRecord->_pDevData;
  gpu->sim.pNLEntry               = gpu->pbNLEntry->_pDevData;
  NLRecords                       = 0;

  // Build Neighbor list records from raw neighbor cell data
  for (i = 0; i < gpu->sim.cells; i++) {
    if (pNLCellRecord[i].NL.neighborCells != 0) {
      int divisor = min(xDivisor, pNLCellRecord[i].NL.neighborCells);
      for (int y = 0; y < yDivisor; y++) {
        for (int x = 0; x < divisor; x++) {
          int start = (x * pNLCellRecord[i].NL.neighborCells) / divisor;
          int end = ((x + 1) * pNLCellRecord[i].NL.neighborCells) / divisor;
          if (end > start) {
            pNLRecord[NLRecords].NL.homeCell = i;
            pNLRecord[NLRecords].NL.neighborCells = (end - start) |
                                                    (y << NLRECORD_YOFFSET_SHIFT);
            for (j = start; j < end; j++) {
              pNLRecord[NLRecords].NL.neighborCell[j] = pNLCellRecord[i].NL.neighborCell[j];
            }
            NLRecords++;
          }
        }
      }
    }
  }
  gpu->sim.NLRecords = NLRecords;
  gpu->pbNLRecord->Upload();

  // Map excessive exclusion mask space
  unsigned int warpsPerCell   = (atomsPerCell + GRID - 1) / GRID;
  unsigned int multiplier     = GRID / gpu->sim.NLAtomsPerWarp;
  unsigned int offsetPerWarp  = GRID +                     // atom ID
                                (gpu->sim.NLAtomsPerWarp * sizeof(PMEMask) /
                                 sizeof(unsigned int));    // masks stored as uints
  unsigned int maxTotalOffset = 2 * workUnits * multiplier * warpsPerCell * warpsPerCell *
                                offsetPerWarp;
  gpu->pbNLAtomList         = new GpuBuffer<unsigned int>(maxTotalOffset,
                                                          bShadowedOutputBuffers);
  gpu->sim.pNLAtomList      = gpu->pbNLAtomList->_pDevData;
  gpu->sim.NLMaxTotalOffset = maxTotalOffset;
  gpu->sim.NLOffsetPerWarp  = offsetPerWarp;
  gpu->pbNLTotalOffset      = new GpuBuffer<unsigned int>(1);
  gpu->sim.pNLTotalOffset   = gpu->pbNLTotalOffset->_pDevData;
  gpu->pbNLEntries          = new GpuBuffer<unsigned int>(1);
  gpu->sim.pNLEntries       = gpu->pbNLEntries->_pDevData;
#ifdef MPI
#  ifdef GVERBOSE
  printf("Node %d, %d work units\n", gpu->gpuID, workUnits);
#  endif
#endif

  // Allocate remapped local interactions
  gpu->pbImageNMRCOMDistanceID       = new GpuBuffer<int2>(gpu->sim.NMRCOMDistances);
  gpu->pbImageNMRCOMDistanceCOM      = new GpuBuffer<int2>((gpu->sim.NMRMaxgrp) *
                                                           (gpu->sim.NMRCOMDistances));
  gpu->pbImageNMRCOMDistanceCOMGrp   = new GpuBuffer<int2>(gpu->sim.NMRCOMDistances * 2);
  gpu->pbImageNMRr6avDistanceID      = new GpuBuffer<int2>(gpu->sim.NMRr6avDistances);
  gpu->pbImageNMRr6avDistancer6av    = new GpuBuffer<int2>((gpu->sim.NMRMaxgrp) *
                                                           (gpu->sim.NMRr6avDistances));
  gpu->pbImageNMRr6avDistancer6avGrp = new GpuBuffer<int2>(gpu->sim.NMRr6avDistances * 2);

  gpu->pbImageNMRCOMAngleID1         = new GpuBuffer<int2>(gpu->sim.NMRCOMAngles);
  gpu->pbImageNMRCOMAngleID2         = new GpuBuffer<int>(gpu->sim.NMRCOMAngles);
  gpu->pbImageNMRCOMAngleCOM      = new GpuBuffer<int2>((gpu->sim.NMRMaxgrp) *
                                                           (gpu->sim.NMRCOMAngles));
  gpu->pbImageNMRCOMAngleCOMGrp   = new GpuBuffer<int2>(gpu->sim.NMRCOMAngles * 3);

  gpu->pbImageNMRCOMTorsionID1       = new GpuBuffer<int4>(gpu->sim.NMRCOMTorsions);
  gpu->pbImageNMRCOMTorsionCOM      = new GpuBuffer<int2>((gpu->sim.NMRMaxgrp) *
                                                           (gpu->sim.NMRCOMTorsions));
  gpu->pbImageNMRCOMTorsionCOMGrp   = new GpuBuffer<int2>(gpu->sim.NMRCOMTorsions * 4);


  gpu->pbImageShakeID                = new GpuBuffer<int4>(gpu->sim.shakeConstraints);
  gpu->pbImageFastShakeID            = new GpuBuffer<int4>(gpu->sim.fastShakeConstraints);
  gpu->pbImageSlowShakeID1           = new GpuBuffer<int>(gpu->sim.slowShakeConstraints);
  gpu->pbImageSlowShakeID2           = new GpuBuffer<int4>(gpu->sim.slowShakeConstraints);
  if(*ti_shake == 2) {
    gpu->pbImageSlowTIShakeID1         = new GpuBuffer<int>(gpu->sim.slowTIShakeConstraints);
    gpu->pbImageSlowTIShakeID2         = new GpuBuffer<int4>(gpu->sim.slowTIShakeConstraints);
    gpu->pbImageSlowTIShakeID3         = new GpuBuffer<int4>(gpu->sim.slowTIShakeConstraints);
  }
  if ((gpu->sim.solventMolecules > 0) || (gpu->sim.soluteAtoms > 0)) {
    gpu->pbImageSolventAtomID = new GpuBuffer<int4>(gpu->sim.solventMolecules);
    gpu->pbImageSoluteAtomID = new GpuBuffer<int>(gpu->sim.soluteAtoms);
  }

  // Allocate new lists
  gpu->sim.imageStride  = ((gpu->sim.stride + 255) >> 8) << 8;
  gpu->pbImageIndex     = new GpuBuffer<unsigned int>(7 * gpu->sim.imageStride);
  gpu->pbSubImageLookup = new GpuBuffer<unsigned int>(gpu->sim.stride);
  gpu->pbAtomXYSaveSP   = new GpuBuffer<PMEFloat2>(gpu->sim.stride);
  gpu->pbAtomZSaveSP    = new GpuBuffer<PMEFloat>(gpu->sim.stride);
  gpu->pbImage          = new GpuBuffer<double>(6 * gpu->sim.stride);
  gpu->pbImageVel       = new GpuBuffer<double>(6 * gpu->sim.stride);
  gpu->pbImageLVel      = new GpuBuffer<double>(6 * gpu->sim.stride);
  gpu->pbImageMass      = new GpuBuffer<double>(4 * gpu->sim.stride);
  gpu->pbImageCharge    = new GpuBuffer<double>(2 * gpu->sim.stride);
  gpu->pbImageSigEps    = new GpuBuffer<PMEFloat2>(2 * gpu->sim.stride);
  gpu->pbImageLJID      = new GpuBuffer<unsigned int>(2 * gpu->sim.stride);
  gpu->pbImageCellID    = new GpuBuffer<unsigned int>(2 * gpu->sim.stride);
  if (gpu->sim.ti_mode != 0) {
    gpu->pbImageTIRegion = new GpuBuffer<int>(2 * gpu->sim.stride);
    if (gpu->sim.TIPaddedLinearAtmCnt > 0) {
      gpu->pbImageTILinearAtmID = new GpuBuffer<int>(4 * gpu->sim.TIPaddedLinearAtmCnt);
      gpu->pbUpdateIndex        = new GpuBuffer<unsigned int>(gpu->sim.stride);
    }
  }

  // Allocate new extra point arrays
  gpu->pbImageExtraPoint11Frame = new GpuBuffer<int4>(gpu->sim.EP11s);
  gpu->pbImageExtraPoint11Index = new GpuBuffer<int>(gpu->sim.EP11s);
  gpu->pbImageExtraPoint12Frame = new GpuBuffer<int4>(gpu->sim.EP12s);
  gpu->pbImageExtraPoint12Index = new GpuBuffer<int>(gpu->sim.EP12s);
  gpu->pbImageExtraPoint21Frame = new GpuBuffer<int4>(gpu->sim.EP21s);
  gpu->pbImageExtraPoint21Index = new GpuBuffer<int2>(gpu->sim.EP21s);
  gpu->pbImageExtraPoint22Frame = new GpuBuffer<int4>(gpu->sim.EP22s);
  gpu->pbImageExtraPoint22Index = new GpuBuffer<int2>(gpu->sim.EP22s);
  gpu->pbImageEPCustomFrame = new GpuBuffer<int4>(gpu->sim.EPCustomCount);
  gpu->pbImageEPCustomIndex = new GpuBuffer<int>(gpu->sim.EPCustomCount);

  // Copy data
  for (i = 0; i < gpu->sim.stride; i++) {
    int istr = i + gpu->sim.stride;
    int istr2 = i + gpu->sim.stride2;
    gpu->pbImageIndex->_pSysData[i] = i;
    gpu->pbImageIndex->_pSysData[i + gpu->sim.imageStride] = i;
    gpu->pbImageIndex->_pSysData[i + gpu->sim.imageStride * 2] = i;
    gpu->pbImageIndex->_pSysData[i + gpu->sim.imageStride * 3] = 0;
    gpu->pbImageIndex->_pSysData[i + gpu->sim.imageStride * 4] = i;
    gpu->pbImageIndex->_pSysData[i + gpu->sim.imageStride * 5] = i;
    gpu->pbImageIndex->_pSysData[i + gpu->sim.imageStride * 6] = i;
    gpu->pbImage->_pSysData[i]          = gpu->pbAtom->_pSysData[i];
    gpu->pbImage->_pSysData[istr]       = gpu->pbAtom->_pSysData[istr];
    gpu->pbImage->_pSysData[istr2]      = gpu->pbAtom->_pSysData[istr2];
    gpu->pbAtomXYSaveSP->_pSysData[i].x = gpu->pbAtom->_pSysData[i];
    gpu->pbAtomXYSaveSP->_pSysData[i].y = gpu->pbAtom->_pSysData[istr];
    gpu->pbAtomZSaveSP->_pSysData[i]    = gpu->pbAtom->_pSysData[istr2];
    gpu->pbImageVel->_pSysData[i]       = gpu->pbVel->_pSysData[i];
    gpu->pbImageVel->_pSysData[istr]    = gpu->pbVel->_pSysData[istr];
    gpu->pbImageVel->_pSysData[istr2]   = gpu->pbVel->_pSysData[istr2];
    gpu->pbImageLVel->_pSysData[i]      = gpu->pbLVel->_pSysData[i];
    gpu->pbImageLVel->_pSysData[istr]   = gpu->pbLVel->_pSysData[istr];
    gpu->pbImageLVel->_pSysData[istr2]  = gpu->pbLVel->_pSysData[istr2];
    gpu->pbImageMass->_pSysData[i]      = gpu->pbAtomMass->_pSysData[i];
    gpu->pbImageMass->_pSysData[istr]   = gpu->pbAtomMass->_pSysData[istr];
    gpu->pbImageCharge->_pSysData[i]    = gpu->pbAtomCharge->_pSysData[i];
    gpu->pbImageSigEps->_pSysData[i]    = gpu->pbAtomSigEps->_pSysData[i];
    gpu->pbImageLJID->_pSysData[i]      = gpu->pbAtomLJID->_pSysData[i];
    gpu->pbImageCellID->_pSysData[i]    = 0;

    // In it's own loop to avoid having to evaluate the conditional sim.stride times
    if (gpu->sim.ti_mode != 0) {
      for (i = 0; i < gpu->sim.stride; i++) {
        gpu->pbImageTIRegion->_pSysData[i] = gpu->pbTIRegion->_pSysData[i];
        if (gpu->sim.TIPaddedLinearAtmCnt > 0) {
          gpu->pbUpdateIndex->_pSysData[i] = 0;
        }
      }
    }
    if (gpu->sim.ti_mode != 0) {
      for (i = 0; i < gpu->sim.TIPaddedLinearAtmCnt; i++) {
        gpu->pbImageTILinearAtmID->_pSysData[i] = gpu->pbTILinearAtmID->_pSysData[i];
        gpu->pbImageTILinearAtmID->_pSysData[i+gpu->sim.TIPaddedLinearAtmCnt] =
          gpu->pbTILinearAtmID->_pSysData[i+gpu->sim.TIPaddedLinearAtmCnt];
      }
    }
  }
  gpu->pbImageIndex->Upload();
  gpu->pbImage->Upload();
  gpu->pbImageVel->Upload();
  gpu->pbImageLVel->Upload();
  gpu->pbImageMass->Upload();
  gpu->pbImageCharge->Upload();
  gpu->pbImageSigEps->Upload();
  gpu->pbImageLJID->Upload();
  gpu->pbImageCellID->Upload();
  if (gpu->sim.ti_mode != 0) {
    gpu->pbImageTIRegion->Upload();
    if (gpu->sim.TIPaddedLinearAtmCnt > 0) {
      gpu->pbImageTILinearAtmID->Upload();
      gpu->pbUpdateIndex->Upload();
    }
  }

  // Copy bonded interactions
  for (i = 0; i < gpu->sim.NMRCOMDistances; i++) {
    gpu->pbImageNMRCOMDistanceID->_pSysData[i] = gpu->pbNMRCOMDistanceID->_pSysData[i];
  }
  for (j = 0; j < gpu->sim.NMRCOMDistances * 2; j++) {
    int lbound = gpu->pbNMRCOMDistanceCOMGrp->_pSysData[j].x;
    int ubound = gpu->pbNMRCOMDistanceCOMGrp->_pSysData[j].y;
    int counter=0;
    if (j % 2 == 0) {
      counter = 0;
    }
    for (i = lbound; i < ubound; i++) {
      gpu->pbImageNMRCOMDistanceCOM->_pSysData[counter + (gpu->sim.NMRMaxgrp * (j/2))] =
        gpu->pbNMRCOMDistanceCOM->_pSysData[i];
      counter = counter + 1;
    }
  }
  for (i = 0; i < gpu->sim.NMRCOMDistances * 2; i++) {
    gpu->pbImageNMRCOMDistanceCOMGrp->_pSysData[i] = gpu->pbNMRCOMDistanceCOMGrp->_pSysData[i];
  }

  //r6 average
  for (i = 0; i < gpu->sim.NMRr6avDistances; i++) {
    gpu->pbImageNMRr6avDistanceID->_pSysData[i] = gpu->pbNMRr6avDistanceID->_pSysData[i];
  }
  for (j = 0; j < gpu->sim.NMRr6avDistances * 2; j++) {
    int lbound = gpu->pbNMRr6avDistancer6avGrp->_pSysData[j].x;
    int ubound = gpu->pbNMRr6avDistancer6avGrp->_pSysData[j].y;
    int counter=0;
    if (j % 2 == 0) {
      counter = 0;
    }
    for (i = lbound; i < ubound; i++) {
      gpu->pbImageNMRr6avDistancer6av->_pSysData[counter + (gpu->sim.NMRMaxgrp * (j/2))] =
        gpu->pbNMRr6avDistancer6av->_pSysData[i];
      counter = counter + 1;
    }
  }
  for (i = 0; i < gpu->sim.NMRr6avDistances * 2; i++) {
    gpu->pbImageNMRr6avDistancer6avGrp->_pSysData[i] =
      gpu->pbNMRr6avDistancer6avGrp->_pSysData[i];
  }

  // COM Angles
  for (i = 0; i < gpu->sim.NMRCOMAngles; i++) {

    // Get the first three elements of nmrat array, which is the first index into each group
    gpu->pbImageNMRCOMAngleID1->_pSysData[i] = gpu->pbNMRCOMAngleID1->_pSysData[i];
    gpu->pbImageNMRCOMAngleID2->_pSysData[i] = gpu->pbNMRCOMAngleID2->_pSysData[i];
  }

  // Get the upper and lower bound for each igr (nmrat[iat], nmrat[8*iat])
  for (j = 0; j < gpu->sim.NMRCOMAngles * 3; j++)  {
    int lbound = gpu->pbNMRCOMAngleCOMGrp->_pSysData[j].x;
    int ubound = gpu->pbNMRCOMAngleCOMGrp->_pSysData[j].y;
    int counter=0;
    if ( j % 3 == 0){
      counter = 0;
    }

    // Get the atom index within the upper lower bound for an igr
    for (i = lbound; i < ubound; i++){
      gpu->pbImageNMRCOMAngleCOM->_pSysData[counter + (gpu->sim.NMRMaxgrp * (j/3))] =
        gpu->pbNMRCOMAngleCOM->_pSysData[i];
      counter = counter + 1;
    }
  }
  for (i = 0; i < gpu->sim.NMRCOMAngles * 3; i++){
    gpu->pbImageNMRCOMAngleCOMGrp->_pSysData[i]      = gpu->pbNMRCOMAngleCOMGrp->_pSysData[i];
  }

  // COM Torions
  for (i = 0; i < gpu->sim.NMRCOMTorsions; i++){
    gpu->pbImageNMRCOMTorsionID1->_pSysData[i]       = gpu->pbNMRCOMTorsionID1->_pSysData[i];
  }
  for (j = 0; j < gpu->sim.NMRCOMTorsions * 4; j++) {
    int lbound = gpu->pbNMRCOMTorsionCOMGrp->_pSysData[j].x;
    int ubound = gpu->pbNMRCOMTorsionCOMGrp->_pSysData[j].y;
    int counter=0;
    if ( j % 4 == 0){
      counter = 0;
    }
    for (i = lbound; i < ubound; i++){
      gpu->pbImageNMRCOMTorsionCOM->_pSysData[counter + (gpu->sim.NMRMaxgrp * (j/4))]  =
        gpu->pbNMRCOMTorsionCOM->_pSysData[i];
      counter = counter + 1;
    }
  }
  for (i = 0; i < gpu->sim.NMRCOMTorsions * 4; i++){
    gpu->pbImageNMRCOMTorsionCOMGrp->_pSysData[i] = gpu->pbNMRCOMTorsionCOMGrp->_pSysData[i];
  }
  for (i = 0; i < gpu->sim.shakeConstraints; i++) {
    gpu->pbImageShakeID->_pSysData[i] = gpu->pbShakeID->_pSysData[i];
  }
  for (i = 0; i < gpu->sim.fastShakeConstraints; i++) {
    gpu->pbImageFastShakeID->_pSysData[i] = gpu->pbFastShakeID->_pSysData[i];
  }
  for (i = 0; i < gpu->sim.slowShakeConstraints; i++) {
    gpu->pbImageSlowShakeID1->_pSysData[i] = gpu->pbSlowShakeID1->_pSysData[i];
    gpu->pbImageSlowShakeID2->_pSysData[i] = gpu->pbSlowShakeID2->_pSysData[i];
  }
  if(*ti_shake == 2) {
    for (i = 0; i < gpu->sim.slowTIShakeConstraints; i++) {
      gpu->pbImageSlowTIShakeID1->_pSysData[i] = gpu->pbSlowTIShakeID1->_pSysData[i];
      gpu->pbImageSlowTIShakeID2->_pSysData[i] = gpu->pbSlowTIShakeID2->_pSysData[i];
      gpu->pbImageSlowTIShakeID3->_pSysData[i] = gpu->pbSlowTIShakeID3->_pSysData[i];
    }
  }
  for (i = 0; i < gpu->sim.solventMolecules; i++) {
    gpu->pbImageSolventAtomID->_pSysData[i] = gpu->pbSolventAtomID->_pSysData[i];
  }
  for (i = 0; i < gpu->sim.soluteAtoms; i++) {
    gpu->pbImageSoluteAtomID->_pSysData[i] = gpu->pbSoluteAtomID->_pSysData[i];
  }
  for (i = 0; i < gpu->sim.EP11s; i++) {
    gpu->pbImageExtraPoint11Frame->_pSysData[i] = gpu->pbExtraPoint11Frame->_pSysData[i];
    gpu->pbImageExtraPoint11Index->_pSysData[i] = gpu->pbExtraPoint11Index->_pSysData[i];
  }
  for (i = 0; i < gpu->sim.EP12s; i++) {
    gpu->pbImageExtraPoint12Frame->_pSysData[i] = gpu->pbExtraPoint12Frame->_pSysData[i];
    gpu->pbImageExtraPoint12Index->_pSysData[i] = gpu->pbExtraPoint12Index->_pSysData[i];
  }
  for (i = 0; i < gpu->sim.EP21s; i++) {
    gpu->pbImageExtraPoint21Frame->_pSysData[i] = gpu->pbExtraPoint21Frame->_pSysData[i];
    gpu->pbImageExtraPoint21Index->_pSysData[i] = gpu->pbExtraPoint21Index->_pSysData[i];
  }
  for (i = 0; i < gpu->sim.EP22s; i++) {
    gpu->pbImageExtraPoint22Frame->_pSysData[i] = gpu->pbExtraPoint22Frame->_pSysData[i];
    gpu->pbImageExtraPoint22Index->_pSysData[i] = gpu->pbExtraPoint22Index->_pSysData[i];
  }
  if (gpu->pbNMRCOMDistanceID) {
    gpu->pbImageNMRCOMDistanceID->Upload();
  }
  if (gpu->pbNMRCOMDistanceCOM) {
    gpu->pbImageNMRCOMDistanceCOM->Upload();
  }
  if (gpu->pbNMRCOMDistanceCOMGrp) {
    gpu->pbImageNMRCOMDistanceCOMGrp->Upload();
  }
  if (gpu->pbNMRr6avDistanceID) {
    gpu->pbImageNMRr6avDistanceID->Upload();
  }
  if (gpu->pbNMRr6avDistancer6av) {
    gpu->pbImageNMRr6avDistancer6av->Upload();
  }
  if (gpu->pbNMRr6avDistancer6avGrp) {
    gpu->pbImageNMRr6avDistancer6avGrp->Upload();
  }
  if (gpu->pbNMRCOMAngleID1) {
    gpu->pbImageNMRCOMAngleID1->Upload();
    gpu->pbImageNMRCOMAngleID2->Upload();
  }
  if (gpu->pbNMRCOMAngleCOM) {
    gpu->pbImageNMRCOMAngleCOM->Upload();
  }
  if (gpu->pbNMRCOMAngleCOMGrp) {
    gpu->pbImageNMRCOMAngleCOMGrp->Upload();
  }
  if (gpu->pbNMRCOMTorsionID1) {
    gpu->pbImageNMRCOMTorsionID1->Upload();
  }
  if (gpu->pbNMRCOMTorsionCOM) {
    gpu->pbImageNMRCOMTorsionCOM->Upload();
  }
  if (gpu->pbNMRCOMTorsionCOMGrp) {
    gpu->pbImageNMRCOMTorsionCOMGrp->Upload();
  }
  gpu->pbImageShakeID->Upload();
  gpu->pbImageFastShakeID->Upload();
  gpu->pbImageSlowShakeID1->Upload();
  gpu->pbImageSlowShakeID2->Upload();
  if(*ti_shake == 2) {
    gpu->pbImageSlowTIShakeID1->Upload();
    gpu->pbImageSlowTIShakeID2->Upload();
    gpu->pbImageSlowTIShakeID3->Upload();
  }
  if ((gpu->sim.solventMolecules > 0) || (gpu->sim.soluteAtoms > 0)) {
    gpu->pbImageSolventAtomID->Upload();
    gpu->pbImageSoluteAtomID->Upload();
  }
  if (gpu->sim.EP11s > 0) {
    gpu->pbImageExtraPoint11Frame->Upload();
    gpu->pbImageExtraPoint11Index->Upload();
  }
  if (gpu->sim.EP12s > 0) {
    gpu->pbImageExtraPoint12Frame->Upload();
    gpu->pbImageExtraPoint12Index->Upload();
  }
  if (gpu->sim.EP21s > 0) {
    gpu->pbImageExtraPoint21Frame->Upload();
    gpu->pbImageExtraPoint21Index->Upload();
  }
  if (gpu->sim.EP22s > 0) {
    gpu->pbImageExtraPoint22Frame->Upload();
    gpu->pbImageExtraPoint22Index->Upload();
  }

  // Set up pointers
  gpu->sim.pImageIndex            = gpu->pbImageIndex->_pDevData;
  gpu->sim.pImageAtom             = gpu->pbImageIndex->_pDevData + gpu->sim.imageStride;
  gpu->sim.pImageAtomLookup       = gpu->pbImageIndex->_pDevData + gpu->sim.imageStride * 2;
  gpu->sim.pImageHash             = gpu->pbImageIndex->_pDevData + gpu->sim.imageStride * 3;
  gpu->sim.pImageIndex2           = gpu->pbImageIndex->_pDevData + gpu->sim.imageStride * 4;
  gpu->sim.pImageAtom2            = gpu->pbImageIndex->_pDevData + gpu->sim.imageStride * 5;
  gpu->sim.pImageHash2            = gpu->pbImageIndex->_pDevData + gpu->sim.imageStride * 6;
  gpu->sim.pAtomXYSaveSP          = gpu->pbAtomXYSaveSP->_pDevData;
  gpu->sim.pAtomZSaveSP           = gpu->pbAtomZSaveSP->_pDevData;
  gpu->sim.pImageX                = gpu->pbImage->_pDevData;
  gpu->sim.pImageY                = gpu->pbImage->_pDevData + gpu->sim.stride;
  gpu->sim.pImageZ                = gpu->pbImage->_pDevData + gpu->sim.stride2;
  gpu->sim.pImageX2               = gpu->pbImage->_pDevData + gpu->sim.stride3;
  gpu->sim.pImageY2               = gpu->pbImage->_pDevData + gpu->sim.stride4;
  gpu->sim.pImageZ2               = gpu->pbImage->_pDevData + gpu->sim.stride * 5;
  {
    cudaError_t status;
    if (gpu->sim.texImageX != 0) {
      status = cudaDestroyTextureObject(gpu->sim.texImageX);
      RTERROR(status, "cudaDestroyTextureObject gpu->sim.texImageX failed");
    }
    cudaResourceDesc resDesc;
    cudaTextureDesc texDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceTypeLinear;
    resDesc.res.linear.devPtr = gpu->sim.pImageX;
    resDesc.res.linear.sizeInBytes = gpu->sim.stride3 * sizeof(int2);
    resDesc.res.linear.desc.f = cudaChannelFormatKindSigned;
    resDesc.res.linear.desc.x = 32;
    resDesc.res.linear.desc.y = 32;
    resDesc.res.linear.desc.z = 0;
    resDesc.res.linear.desc.w = 0;
    memset(&texDesc, 0, sizeof(texDesc));
    texDesc.normalizedCoords = 0;
    texDesc.filterMode = cudaFilterModePoint;
    texDesc.addressMode[0] = cudaAddressModeClamp;
    texDesc.readMode = cudaReadModeElementType;
    status = cudaCreateTextureObject(&(gpu->sim.texImageX), &resDesc, &texDesc, NULL);
    RTERROR(status, "cudaCreateTextureObject gpu->sim.texImageX failed");
  }
  {
    cudaError_t status;
    if (gpu->sim.texImageX2 != 0) {
      status = cudaDestroyTextureObject(gpu->sim.texImageX2);
      RTERROR(status, "cudaDestroyTextureObject gpu->sim.texImageX2 failed");
    }
    cudaResourceDesc resDesc;
    cudaTextureDesc texDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceTypeLinear;
    resDesc.res.linear.devPtr = gpu->sim.pImageX2;
    resDesc.res.linear.sizeInBytes = gpu->sim.stride3 * sizeof(int2);
    resDesc.res.linear.desc.f = cudaChannelFormatKindSigned;
    resDesc.res.linear.desc.x = 32;
    resDesc.res.linear.desc.y = 32;
    resDesc.res.linear.desc.z = 0;
    resDesc.res.linear.desc.w = 0;
    memset(&texDesc, 0, sizeof(texDesc));
    texDesc.normalizedCoords = 0;
    texDesc.filterMode = cudaFilterModePoint;
    texDesc.addressMode[0] = cudaAddressModeClamp;
    texDesc.readMode = cudaReadModeElementType;
    status = cudaCreateTextureObject(&(gpu->sim.texImageX2), &resDesc, &texDesc, NULL);
    RTERROR(status, "cudaCreateTextureObject gpu->sim.texImageX2 failed");
  }
  gpu->sim.pImageVelX             = gpu->pbImageVel->_pDevData;
  gpu->sim.pImageVelY             = gpu->pbImageVel->_pDevData + gpu->sim.stride;
  gpu->sim.pImageVelZ             = gpu->pbImageVel->_pDevData + gpu->sim.stride2;
  gpu->sim.pImageVelX2            = gpu->pbImageVel->_pDevData + gpu->sim.stride3;
  gpu->sim.pImageVelY2            = gpu->pbImageVel->_pDevData + gpu->sim.stride4;
  gpu->sim.pImageVelZ2            = gpu->pbImageVel->_pDevData + gpu->sim.stride * 5;
  gpu->sim.pImageLVelX            = gpu->pbImageLVel->_pDevData;
  gpu->sim.pImageLVelY            = gpu->pbImageLVel->_pDevData + gpu->sim.stride;
  gpu->sim.pImageLVelZ            = gpu->pbImageLVel->_pDevData + gpu->sim.stride2;
  gpu->sim.pImageLVelX2           = gpu->pbImageLVel->_pDevData + gpu->sim.stride3;
  gpu->sim.pImageLVelY2           = gpu->pbImageLVel->_pDevData + gpu->sim.stride4;
  gpu->sim.pImageLVelZ2           = gpu->pbImageLVel->_pDevData + gpu->sim.stride * 5;
  gpu->sim.pImageMass             = gpu->pbImageMass->_pDevData;
  gpu->sim.pImageInvMass          = gpu->pbImageMass->_pDevData + gpu->sim.stride;
  gpu->sim.pImageMass2            = gpu->pbImageMass->_pDevData + gpu->sim.stride2;
  gpu->sim.pImageInvMass2         = gpu->pbImageMass->_pDevData + gpu->sim.stride3;
  gpu->sim.pImageCharge           = gpu->pbImageCharge->_pDevData;
  gpu->sim.pImageCharge2          = gpu->pbImageCharge->_pDevData + gpu->sim.stride;
  gpu->sim.pImageSigEps           = gpu->pbImageSigEps->_pDevData;
  gpu->sim.pImageSigEps2          = gpu->pbImageSigEps->_pDevData + gpu->sim.stride;
  gpu->sim.pImageLJID             = gpu->pbImageLJID->_pDevData;
  gpu->sim.pImageLJID2            = gpu->pbImageLJID->_pDevData + gpu->sim.stride;
  gpu->sim.pImageCellID           = gpu->pbImageCellID->_pDevData;
  gpu->sim.pImageCellID2          = gpu->pbImageCellID->_pDevData + gpu->sim.stride;
  if (gpu->sim.ti_mode != 0) {
    gpu->sim.pImageTIRegion  = gpu->pbImageTIRegion->_pDevData;
    gpu->sim.pImageTIRegion2 = gpu->pbImageTIRegion->_pDevData + gpu->sim.stride;
    if (gpu->sim.TIPaddedLinearAtmCnt > 0) {
      gpu->sim.pImageTILinearAtmID  = gpu->pbImageTILinearAtmID->_pDevData;
      gpu->sim.pImageTILinearAtmID2 = gpu->pbImageTILinearAtmID->_pDevData +
                                      2*gpu->sim.TIPaddedLinearAtmCnt;
      gpu->sim.pUpdateIndex         = gpu->pbUpdateIndex->_pDevData;
    }
  }
  if ((gpu->sim.solventMolecules > 0) || (gpu->sim.soluteAtoms > 0)) {
    gpu->sim.pImageSolventAtomID = gpu->pbImageSolventAtomID->_pDevData;
    gpu->sim.pImageSoluteAtomID = gpu->pbImageSoluteAtomID->_pDevData;
  }
  gpu->sim.pImageNMRCOMDistanceID       = gpu->pbImageNMRCOMDistanceID->_pDevData;
  gpu->sim.pImageNMRCOMDistanceCOM      = gpu->pbImageNMRCOMDistanceCOM->_pDevData;
  gpu->sim.pImageNMRCOMDistanceCOMGrp   = gpu->pbImageNMRCOMDistanceCOMGrp->_pDevData;
  gpu->sim.pImageNMRr6avDistanceID      = gpu->pbImageNMRr6avDistanceID->_pDevData;
  gpu->sim.pImageNMRr6avDistancer6av    = gpu->pbImageNMRr6avDistancer6av->_pDevData;
  gpu->sim.pImageNMRr6avDistancer6avGrp = gpu->pbImageNMRr6avDistancer6avGrp->_pDevData;
  gpu->sim.pImageNMRCOMAngleID1       = gpu->pbImageNMRCOMAngleID1->_pDevData;
  gpu->sim.pImageNMRCOMAngleID2       = gpu->pbImageNMRCOMAngleID2->_pDevData;
  gpu->sim.pImageNMRCOMAngleCOM      = gpu->pbImageNMRCOMAngleCOM->_pDevData;
  gpu->sim.pImageNMRCOMAngleCOMGrp   = gpu->pbImageNMRCOMAngleCOMGrp->_pDevData;
  gpu->sim.pImageNMRCOMTorsionID1       = gpu->pbImageNMRCOMTorsionID1->_pDevData;
  gpu->sim.pImageNMRCOMTorsionCOM      = gpu->pbImageNMRCOMTorsionCOM->_pDevData;
  gpu->sim.pImageNMRCOMTorsionCOMGrp   = gpu->pbImageNMRCOMTorsionCOMGrp->_pDevData;
  gpu->sim.pImageShakeID                = gpu->pbImageShakeID->_pDevData;
  gpu->sim.pImageFastShakeID            = gpu->pbImageFastShakeID->_pDevData;
  gpu->sim.pImageSlowShakeID1           = gpu->pbImageSlowShakeID1->_pDevData;
  gpu->sim.pImageSlowShakeID2           = gpu->pbImageSlowShakeID2->_pDevData;
  if(*ti_shake == 2) {
    gpu->sim.pImageSlowTIShakeID1         = gpu->pbImageSlowTIShakeID1->_pDevData;
    gpu->sim.pImageSlowTIShakeID2         = gpu->pbImageSlowTIShakeID2->_pDevData;
    gpu->sim.pImageSlowTIShakeID3         = gpu->pbImageSlowTIShakeID3->_pDevData;
  }
  gpu->sim.pImageExtraPoint11Frame      = gpu->pbImageExtraPoint11Frame->_pDevData;
  gpu->sim.pImageExtraPoint11Index      = gpu->pbImageExtraPoint11Index->_pDevData;
  gpu->sim.pImageExtraPoint12Frame      = gpu->pbImageExtraPoint12Frame->_pDevData;
  gpu->sim.pImageExtraPoint12Index      = gpu->pbImageExtraPoint12Index->_pDevData;
  gpu->sim.pImageExtraPoint21Frame      = gpu->pbImageExtraPoint21Frame->_pDevData;
  gpu->sim.pImageExtraPoint21Index      = gpu->pbImageExtraPoint21Index->_pDevData;
  gpu->sim.pImageExtraPoint22Frame      = gpu->pbImageExtraPoint22Frame->_pDevData;
  gpu->sim.pImageExtraPoint22Index      = gpu->pbImageExtraPoint22Index->_pDevData;
  gpu->sim.pImageEPCustomFrame          = gpu->pbImageEPCustomFrame->_pDevData;
  gpu->sim.pImageEPCustomIndex          = gpu->pbImageEPCustomIndex->_pDevData;
  gpu->sim.pNLNonbondCellStartEnd       = gpu->pbNLNonbondCellStartEnd->_pDevData;
  gpu->sim.pNLAtomList                  = gpu->pbNLAtomList->_pDevData;
  gpu->sim.pFinishedBlocksCounters      = gpu->pbFinishedBlocksCounters->_pDevData;
  gpu->sim.pKernelCompletionFlags       = gpu->pbKernelCompletionFlags->_pDevData;
  gpu->sim.pNLbSkinTestFail             = gpu->pbNLbSkinTestFail->_pDevData;
  gpu->sim.pNLCellHash                  = gpu->pbNLCellHash->_pDevData;

  // Set up radix sort
  kNLInitRadixSort(gpu);
  gpuCopyConstants();
  gpu_create_outputbuffers_();

  // Free allocated memory
  delete[] pNLCellRecord;
  delete[] pExclusionCount;
}

//---------------------------------------------------------------------------------------------
// gpu_skin_test_: run the GPU neighbor list skin test (this will kill the simulation if a
//                 violation is detected)
//---------------------------------------------------------------------------------------------
extern "C" void gpu_skin_test_()
{
  PRINTMETHOD("gpu_skin_test");
  kNLSkinTest(gpu);
}

//---------------------------------------------------------------------------------------------
// gpu_build_neighbor_list_: create the neighbor list!
//---------------------------------------------------------------------------------------------
extern "C" void gpu_build_neighbor_list_()
{
  PRINTMETHOD("gpu_build_neighbor_list");
  if (gpu->bNeedNewNeighborList) {
    gpu->bNeedNewNeighborList = false;
    gpu->bNewNeighborList = true;
    kNLGenerateSpatialHash(gpu);
    kNLRadixSort(gpu);
    kNLRemapImage(gpu);
    kNLClearCellBoundaries(gpu);
    kNLCalculateCellBoundaries(gpu);
#ifdef MPI
    // Regenerate force send/receives
    gpu->pbNLNonbondCellStartEnd->Download();

    // Calculate overall limits and PME scatter info
    for (int i = 0; i < gpu->nGpus; i++) {
      if (gpu->pMinLocalCell[i] != gpu->pMaxLocalCell[i]) {

        // We have to detect crazy users who have large areas of vacuum in their simulation.
        // Technically speaking, I suspect we ought to crash the simulation with a nasty
        // warning message but in the interest of world peace, here's a complicated loop to
        // handle this (1 or more empty nonbond cells).
        int j = gpu->pMinLocalCell[i];
        while ((j < gpu->pMaxLocalCell[i]) && (gpu->pbNLNonbondCellStartEnd->_pSysData[j].x ==
                                               gpu->pbNLNonbondCellStartEnd->_pSysData[j].y)) {
          j++;
        }
        if (j < gpu->pMaxLocalCell[i]) {
          gpu->pMinLocalAtom[i] = gpu->pbNLNonbondCellStartEnd->_pSysData[j].x;
          int k = gpu->pMaxLocalCell[i] - 1;

          // We'll never get here unless there's at least one nonempty
          // cell so no need to check cell bounds
          while (gpu->pbNLNonbondCellStartEnd->_pSysData[k].x ==
                 gpu->pbNLNonbondCellStartEnd->_pSysData[k].y) {
            k--;
          }
          gpu->pMaxLocalAtom[i] = gpu->pbNLNonbondCellStartEnd->_pSysData[k].y;
        }
        else {
          gpu->pMinLocalAtom[i] = 0;
          gpu->pMaxLocalAtom[i] = 0;
        }
      }
      else {
        gpu->pMinLocalAtom[i] = 0;
        gpu->pMaxLocalAtom[i] = 0;
      }
    }

    // Calculate local limits
    if (gpu->minLocalCell != gpu->maxLocalCell) {
      gpu->sim.minLocalAtom = gpu->pMinLocalAtom[gpu->gpuID];
      gpu->sim.maxLocalAtom = gpu->pMaxLocalAtom[gpu->gpuID];
    }
    else {
      gpu->sim.minLocalAtom = 0;
      gpu->sim.maxLocalAtom = 0;
    }
#endif
    gpuCopyConstants();
    kNLCalculateCellCoordinates(gpu);
    kNLBuildNeighborList(gpu);
    kNLRemapLocalInteractions(gpu);
    kNLRemapBondWorkUnits(gpu);
#ifdef MPI
    if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
      kClearNBForces(gpu);
    }
#endif
  }
}

//---------------------------------------------------------------------------------------------
// gpu_molecule_list_setup_: set up the list of molecules.  Solute and solvent molecules are
//                           delineated here.
//
// Arguments:
//   molecules:   the total number of molecules (known from the topology)
//   listdata:    a list of all the different types of molecules (i.e. biotin, streptavidin,
//                water, sodium, chloride)
//---------------------------------------------------------------------------------------------
extern "C" void gpu_molecule_list_setup_(int* molecules, listdata_rec listdata[])
{
  PRINTMETHOD("gpu_molecule_list_setup");

  // Delete previous molecule list
  delete gpu->pbSoluteAtomID;
  delete gpu->pbSoluteAtomMass;
  delete gpu->pbSolute;
  delete gpu->pbUllSolute;
  delete gpu->pbSolventAtomID;
  delete gpu->pbSolvent;

  // Count molecule parameters
  int soluteAtoms      = 0;
  int soluteMolecules  = 0;
  int solventMolecules = 0;
  for (int i = 0; i < *molecules; i++) {

    // Distinguish between solute and solvent
    if (listdata[i].cnt < 5) {
      solventMolecules++;
    }
    else {
      int offset =
        ((listdata[i].cnt + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits;
      soluteAtoms += offset;
      soluteMolecules++;
    }
  }
  int soluteMoleculeStride =
    ((soluteMolecules + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits;
  int solventMoleculeStride =
    ((solventMolecules + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits;

  // Allocate solvent/solute data
  gpu->pbSoluteAtomID = new GpuBuffer<int>(2 * soluteAtoms);
  gpu->pbSoluteAtomMass = new GpuBuffer<PMEDouble>(soluteAtoms);
  if (soluteMolecules <= gpu->maxPSSoluteMolecules) {
    gpu->pbSolute = new GpuBuffer<PMEDouble>(4 * soluteMoleculeStride);
  }
  else {
    gpu->pbSolute = new GpuBuffer<PMEDouble>(7 * soluteMoleculeStride);
  }
  gpu->pbUllSolute = new GpuBuffer<PMEUllInt>(6 * soluteMoleculeStride);
  gpu->pbSolventAtomID = new GpuBuffer<int4>(solventMoleculeStride);
  gpu->pbSolvent = new GpuBuffer<PMEDouble>(8 * solventMoleculeStride);

  // Fill in data
  int* pSoluteAtomID         = gpu->pbSoluteAtomID->_pSysData;
  int* pSoluteAtomMoleculeID = gpu->pbSoluteAtomID->_pSysData + soluteAtoms;
  PMEDouble* pSoluteAtomMass = gpu->pbSoluteAtomMass->_pSysData;
  PMEDouble* pSoluteInvMass  = gpu->pbSolute->_pSysData + soluteMoleculeStride * 3;
  int4* pSolventAtomID       = gpu->pbSolventAtomID->_pSysData;
  PMEDouble* pSolventMass1   = gpu->pbSolvent->_pSysData;
  PMEDouble* pSolventMass2   = gpu->pbSolvent->_pSysData + solventMoleculeStride;
  PMEDouble* pSolventMass3   = gpu->pbSolvent->_pSysData + solventMoleculeStride * 2;
  PMEDouble* pSolventMass4   = gpu->pbSolvent->_pSysData + solventMoleculeStride * 3;
  PMEDouble* pSolventInvMass = gpu->pbSolvent->_pSysData + solventMoleculeStride * 7;
  soluteAtoms      = 0;
  soluteMolecules  = 0;
  solventMolecules = 0;
  for (int i = 0; i < *molecules; i++) {

    // Distinguish between solute and solvent
    double totalMass = 0.0;
    if (listdata[i].cnt < 5) {
      pSolventAtomID[solventMolecules].y = -1;
      pSolventAtomID[solventMolecules].z = -1;
      pSolventAtomID[solventMolecules].w = -1;
      pSolventMass2[solventMolecules] = 0.0;
      pSolventMass3[solventMolecules] = 0.0;
      pSolventMass4[solventMolecules] = 0.0;
      double totalMass = 0.0;
      switch (listdata[i].cnt) {
        case 4:
          pSolventAtomID[solventMolecules].w = listdata[i].offset + 3;
          pSolventMass4[solventMolecules] = gpu->pbAtomMass->_pSysData[listdata[i].offset + 3];
          totalMass += gpu->pbAtomMass->_pSysData[listdata[i].offset + 3];
        case 3:
          pSolventAtomID[solventMolecules].z = listdata[i].offset + 2;
          pSolventMass3[solventMolecules] = gpu->pbAtomMass->_pSysData[listdata[i].offset + 2];
          totalMass += gpu->pbAtomMass->_pSysData[listdata[i].offset + 2];
        case 2:
          pSolventAtomID[solventMolecules].y = listdata[i].offset + 1;
          pSolventMass2[solventMolecules] = gpu->pbAtomMass->_pSysData[listdata[i].offset + 1];
          totalMass += gpu->pbAtomMass->_pSysData[listdata[i].offset + 1];
        case 1:
          pSolventAtomID[solventMolecules].x = listdata[i].offset;
          pSolventMass1[solventMolecules] = gpu->pbAtomMass->_pSysData[listdata[i].offset];
          totalMass += gpu->pbAtomMass->_pSysData[listdata[i].offset];
      }
      pSolventInvMass[solventMolecules] = (totalMass > 1.0e-8) ? 1.0 / totalMass : 0.0;
      solventMolecules++;
    }
    else {
      double totalMass = 0.0;
      for (int j = 0; j < listdata[i].cnt; j++) {
        pSoluteAtomID[soluteAtoms + j] = listdata[i].offset + j;
        pSoluteAtomMoleculeID[soluteAtoms + j] = soluteMolecules;
        pSoluteAtomMass[soluteAtoms + j] = gpu->pbAtomMass->_pSysData[listdata[i].offset + j];
        totalMass += gpu->pbAtomMass->_pSysData[listdata[i].offset + j];
      }
      int offset = ((listdata[i].cnt + (GRID - 1)) >> GRID_BITS) << GRID_BITS;
      for (int j = listdata[i].cnt; j < offset; j++) {
        pSoluteAtomID[soluteAtoms + j] = -1;
        pSoluteAtomMoleculeID[soluteAtoms + j] = soluteMolecules;
        pSoluteAtomMass[soluteAtoms + j] = 0.0;
      }
      soluteAtoms += offset;
      pSoluteInvMass[soluteMolecules] = 1.0 / totalMass;
      soluteMolecules++;
    }
  }

  // Upload data
  gpu->pbSoluteAtomID->Upload();
  gpu->pbSoluteAtomMass->Upload();
  gpu->pbSolute->Upload();
  gpu->pbUllSolute->Upload();
  gpu->pbSolventAtomID->Upload();
  gpu->pbSolvent->Upload();

  // Set up constant pointers
  gpu->sim.soluteMolecules = soluteMolecules;
  gpu->sim.soluteMoleculeStride = soluteMoleculeStride;
  gpu->sim.soluteAtoms = soluteAtoms;
  gpu->sim.solventMolecules = solventMolecules;
  gpu->sim.solventMoleculeStride = solventMoleculeStride;

  // Solute
  gpu->sim.pSoluteCOMX           = gpu->pbSolute->_pDevData;
  gpu->sim.pSoluteCOMY           = gpu->pbSolute->_pDevData + soluteMoleculeStride;
  gpu->sim.pSoluteCOMZ           = gpu->pbSolute->_pDevData + soluteMoleculeStride * 2;
  gpu->sim.pSoluteInvMass        = gpu->pbSolute->_pDevData + soluteMoleculeStride * 3;
  gpu->sim.pSoluteDeltaCOMX      = gpu->pbSolute->_pDevData + soluteMoleculeStride * 4;
  gpu->sim.pSoluteDeltaCOMY      = gpu->pbSolute->_pDevData + soluteMoleculeStride * 5;
  gpu->sim.pSoluteDeltaCOMZ      = gpu->pbSolute->_pDevData + soluteMoleculeStride * 6;
  gpu->sim.pSoluteUllCOMX        = gpu->pbUllSolute->_pDevData;
  gpu->sim.pSoluteUllCOMY        = gpu->pbUllSolute->_pDevData + soluteMoleculeStride;
  gpu->sim.pSoluteUllCOMZ        = gpu->pbUllSolute->_pDevData + soluteMoleculeStride * 2;
  gpu->sim.pSoluteUllEKCOMX      = gpu->pbUllSolute->_pDevData + soluteMoleculeStride * 3;
  gpu->sim.pSoluteUllEKCOMY      = gpu->pbUllSolute->_pDevData + soluteMoleculeStride * 4;
  gpu->sim.pSoluteUllEKCOMZ      = gpu->pbUllSolute->_pDevData + soluteMoleculeStride * 5;
  gpu->sim.pSoluteAtomID         = gpu->pbSoluteAtomID->_pDevData;
  gpu->sim.pSoluteAtomMoleculeID = gpu->pbSoluteAtomID->_pDevData + soluteAtoms;
  gpu->sim.pSoluteAtomMass       = gpu->pbSoluteAtomMass->_pDevData;

  // Solvent
  gpu->sim.pSolventAtomMass1 = gpu->pbSolvent->_pDevData;
  gpu->sim.pSolventAtomMass2 = gpu->pbSolvent->_pDevData + solventMoleculeStride;
  gpu->sim.pSolventAtomMass3 = gpu->pbSolvent->_pDevData + solventMoleculeStride * 2;
  gpu->sim.pSolventAtomMass4 = gpu->pbSolvent->_pDevData + solventMoleculeStride * 3;
  gpu->sim.pSolventCOMX      = gpu->pbSolvent->_pDevData + solventMoleculeStride * 4;
  gpu->sim.pSolventCOMY      = gpu->pbSolvent->_pDevData + solventMoleculeStride * 5;
  gpu->sim.pSolventCOMZ      = gpu->pbSolvent->_pDevData + solventMoleculeStride * 6;
  gpu->sim.pSolventInvMass   = gpu->pbSolvent->_pDevData + solventMoleculeStride * 7;
  gpu->sim.pSolventAtomID    = gpu->pbSolventAtomID->_pDevData;

  // Another of many invocations of gpuCopyConstants()
  gpuCopyConstants();
#ifdef GVERBOSE
    printf("Found %d molecules, %d solvent, %d solute, and %d effective solute atoms\n",
           *molecules, solventMolecules, soluteMolecules, soluteAtoms);
#endif
}

//---------------------------------------------------------------------------------------------
// gpu_constraint_molecule_list_setup_: setup positional constraints for molecules on the GPU
//
// Arguments:
//   molecules:   the total number of molecules (known from the topology)
//   listdata:    a list of all the different types of molecules (i.e. biotin, streptavidin,
//                water, sodium, chloride)
//   atm_xc:      locations at which constrained atoms are pegged to the board
//---------------------------------------------------------------------------------------------
extern "C" void gpu_constraint_molecule_list_setup_(int* molecules, listdata_rec listdata[],
                                                    double atm_xc[])
{
  PRINTMETHOD("gpu_constraint_molecule_list_setup");

  // Delete previous constained molecule list
  delete gpu->pbConstraintAtomX;
  delete gpu->pbConstraintAtomY;
  delete gpu->pbConstraintAtomZ;
  delete gpu->pbConstraintCOMX;
  delete gpu->pbConstraintCOMY;
  delete gpu->pbConstraintCOMZ;

  // Allocate and clear molecule IDs for NTP constaints
  // (assumes regular molecule list is already generated!)
  int* pMoleculeID       = new int[gpu->sim.atoms];
  int* pAtomConstraintID = new int[gpu->sim.atoms];
  for (int i = 0; i < gpu->sim.atoms; i++) {
    pMoleculeID[i] = -1;
    pAtomConstraintID[i] = -1;
  }

  // Calculate centers of mass and assign molecule IDs to all atoms
  int soluteAtoms      = 0;
  int soluteMolecules  = 0;
  int solventMolecules = 0;
  double *pCOMX  = new double[*molecules];
  double *pCOMY  = new double[*molecules];
  double *pCOMZ  = new double[*molecules];
  double *pFCOMX = new double[*molecules];
  double *pFCOMY = new double[*molecules];
  double *pFCOMZ = new double[*molecules];
  for (int i = 0; i < *molecules; i++) {
    double comX = 0.0;
    double comY = 0.0;
    double comZ = 0.0;
    double totalMass = 0.0;
    int offset = ((listdata[i].cnt +
                   (gpu->sim.grid - 1)) >> gpu->sim.gridBits) << gpu->sim.gridBits;
    for (int j = listdata[i].offset; j < listdata[i].offset + listdata[i].cnt; j++) {
      pMoleculeID[j] = i;
      double mass = gpu->pbAtomMass->_pSysData[j];
      comX += mass * atm_xc[j * 3];
      comY += mass * atm_xc[j * 3 + 1];
      comZ += mass * atm_xc[j * 3 + 2];
      totalMass += mass;
    }
    comX /= totalMass;
    comY /= totalMass;
    comZ /= totalMass;
    pCOMX[i] = comX;
    pCOMY[i] = comY;
    pCOMZ[i] = comZ;

    // Now convert to fractional coordinates
    pFCOMX[i] = gpu->sim.recip[0][0] * comX + gpu->sim.recip[1][0] * comY +
                gpu->sim.recip[2][0] * comZ;
    pFCOMY[i] = gpu->sim.recip[1][1] * comY + gpu->sim.recip[2][1] * comZ;
    pFCOMZ[i] = gpu->sim.recip[2][2] * comZ;
  }

  // Create NTP data for each constraint
  gpu->pbConstraintAtomX    = new GpuBuffer<PMEDouble>(gpu->sim.constraints);
  gpu->pbConstraintAtomY    = new GpuBuffer<PMEDouble>(gpu->sim.constraints);
  gpu->pbConstraintAtomZ    = new GpuBuffer<PMEDouble>(gpu->sim.constraints);
  gpu->pbConstraintCOMX     = new GpuBuffer<PMEDouble>(gpu->sim.constraints);
  gpu->pbConstraintCOMY     = new GpuBuffer<PMEDouble>(gpu->sim.constraints);
  gpu->pbConstraintCOMZ     = new GpuBuffer<PMEDouble>(gpu->sim.constraints);
  gpu->sim.pConstraintAtomX = gpu->pbConstraintAtomX->_pDevData;
  gpu->sim.pConstraintAtomY = gpu->pbConstraintAtomY->_pDevData;
  gpu->sim.pConstraintAtomZ = gpu->pbConstraintAtomZ->_pDevData;
  gpu->sim.pConstraintCOMX  = gpu->pbConstraintCOMX->_pDevData;
  gpu->sim.pConstraintCOMY  = gpu->pbConstraintCOMY->_pDevData;
  gpu->sim.pConstraintCOMZ  = gpu->pbConstraintCOMZ->_pDevData;
  for (int i = 0; i < gpu->sim.constraints; i++) {
    int atom = gpu->pbConstraintID->_pSysData[i];
    int moleculeID = pMoleculeID[atom];
    gpu->pbConstraintAtomX->_pSysData[i] = atm_xc[3 * atom]     - pCOMX[moleculeID];
    gpu->pbConstraintAtomY->_pSysData[i] = atm_xc[3 * atom + 1] - pCOMY[moleculeID];
    gpu->pbConstraintAtomZ->_pSysData[i] = atm_xc[3 * atom + 2] - pCOMZ[moleculeID];
    gpu->pbConstraintCOMX->_pSysData[i] = pFCOMX[moleculeID];
    gpu->pbConstraintCOMY->_pSysData[i] = pFCOMY[moleculeID];
    gpu->pbConstraintCOMZ->_pSysData[i] = pFCOMZ[moleculeID];
  }
  gpu->pbConstraintAtomX->Upload();
  gpu->pbConstraintAtomY->Upload();
  gpu->pbConstraintAtomZ->Upload();
  gpu->pbConstraintCOMX->Upload();
  gpu->pbConstraintCOMY->Upload();
  gpu->pbConstraintCOMZ->Upload();
  gpuCopyConstants();

  // Free temporary arrays
  delete[] pMoleculeID;
  delete[] pAtomConstraintID;
  delete[] pCOMX;
  delete[] pCOMY;
  delete[] pCOMZ;
  delete[] pFCOMX;
  delete[] pFCOMY;
  delete[] pFCOMZ;
}

//---------------------------------------------------------------------------------------------
// gpu_ti_molecule_list_setup_: set up the list of molecules for Thermodynamic Integration
//
// Arguments:
//   molecules:       the number of molecules in the system
//   listdata:        critical information on the identity of each molecule
//   ti_mol_type:
//   ti_sc_partner:   soft core partners for each molecule involved in TI
//---------------------------------------------------------------------------------------------
extern "C" void gpu_ti_molecule_list_setup_(int* molecules, listdata_rec listdata[],
                                            int ti_mol_type[], int ti_sc_partner[])
{
  PRINTMETHOD("gpu_ti_molecule_list_setup");

  // Delete previous molecule list
  delete gpu->pbSoluteAtomID;
  delete gpu->pbSoluteAtomMass;
  delete gpu->pbSolute;
  delete gpu->pbUllSolute;
  delete gpu->pbSolventAtomID;
  delete gpu->pbSolvent;
  delete gpu->pbAFEMolType;
  delete gpu->pbAFEMolPartner;

  // Count molecule parameters
  int soluteAtoms      = 0;
  int soluteMolecules  = 0;
  int solventMolecules = 0;
  for (int i = 0; i < *molecules; i++) {

    // Distinguish between solute and solvent
    if (listdata[i].cnt < 5) {
      solventMolecules++;
    }
    else {
      int offset = ((listdata[i].cnt + (gpu->sim.grid - 1)) >> gpu->sim.gridBits) <<
                   gpu->sim.gridBits;
      soluteAtoms += offset;
      soluteMolecules++;
    }
  }
  int soluteMoleculeStride  = ((soluteMolecules + (gpu->sim.grid - 1)) >>
                               gpu->sim.gridBits) << gpu->sim.gridBits;
  int solventMoleculeStride = ((solventMolecules + (gpu->sim.grid - 1)) >>
                               gpu->sim.gridBits) << gpu->sim.gridBits;

  // Allocate solvent/solute data
  gpu->pbSoluteAtomID   = new GpuBuffer<int>(2 * soluteAtoms);
  gpu->pbAFEMolType     = new GpuBuffer<int>(*molecules);
  gpu->pbAFEMolPartner  = new GpuBuffer<int>(*molecules);
  gpu->pbSoluteAtomMass = new GpuBuffer<PMEDouble>(soluteAtoms);
  if (soluteMolecules <= gpu->maxPSSoluteMolecules) {
    gpu->pbSolute = new GpuBuffer<PMEDouble>(4 * soluteMoleculeStride);
  }
  else {
    gpu->pbSolute = new GpuBuffer<PMEDouble>(7 * soluteMoleculeStride);
  }
  gpu->pbUllSolute     = new GpuBuffer<PMEUllInt>(6 * soluteMoleculeStride);
  gpu->pbSolventAtomID = new GpuBuffer<int4>(solventMoleculeStride);
  gpu->pbSolvent       = new GpuBuffer<PMEDouble>(8 * solventMoleculeStride);

  // Fill in data
  int* pSoluteAtomID         = gpu->pbSoluteAtomID->_pSysData;
  int* pSoluteAtomMoleculeID = gpu->pbSoluteAtomID->_pSysData + soluteAtoms;
  PMEDouble* pSoluteAtomMass = gpu->pbSoluteAtomMass->_pSysData;
  PMEDouble* pSoluteInvMass  = gpu->pbSolute->_pSysData + soluteMoleculeStride * 3;
  int4* pSolventAtomID       = gpu->pbSolventAtomID->_pSysData;
  PMEDouble* pSolventMass1   = gpu->pbSolvent->_pSysData;
  PMEDouble* pSolventMass2   = gpu->pbSolvent->_pSysData + solventMoleculeStride;
  PMEDouble* pSolventMass3   = gpu->pbSolvent->_pSysData + solventMoleculeStride * 2;
  PMEDouble* pSolventMass4   = gpu->pbSolvent->_pSysData + solventMoleculeStride * 3;
  PMEDouble* pSolventInvMass = gpu->pbSolvent->_pSysData + solventMoleculeStride * 7;
  soluteAtoms         = 0;
  soluteMolecules     = 0;
  solventMolecules    = 0;
  int* pAFEMolType    = gpu->pbAFEMolType->_pSysData;
  int* pAFEMolPartner = gpu->pbAFEMolPartner->_pSysData;
  for (int i = 0; i < *molecules; i++) {

    // Distinguish between solute and solvent
    double totalMass = 0.0;
    if (listdata[i].cnt < 5) {
      pSolventAtomID[solventMolecules].y = -1;
      pSolventAtomID[solventMolecules].z = -1;
      pSolventAtomID[solventMolecules].w = -1;
      pSolventMass2[solventMolecules]    = 0.0;
      pSolventMass3[solventMolecules]    = 0.0;
      pSolventMass4[solventMolecules]    = 0.0;
      double totalMass                   = 0.0;
      switch (listdata[i].cnt) {
        case 4:
          pSolventAtomID[solventMolecules].w = listdata[i].offset + 3;
          pSolventMass4[solventMolecules]    = gpu->pbAtomMass->_pSysData[listdata[i].offset +
                                                                          3];
          totalMass                         += gpu->pbAtomMass->_pSysData[listdata[i].offset +
                                                                          3];
        case 3:
          pSolventAtomID[solventMolecules].z = listdata[i].offset + 2;
          pSolventMass3[solventMolecules]    = gpu->pbAtomMass->_pSysData[listdata[i].offset +
                                                                          2];
          totalMass                         += gpu->pbAtomMass->_pSysData[listdata[i].offset +
                                                                          2];
        case 2:
          pSolventAtomID[solventMolecules].y = listdata[i].offset + 1;
          pSolventMass2[solventMolecules]    = gpu->pbAtomMass->_pSysData[listdata[i].offset +
                                                                          1];
          totalMass                         += gpu->pbAtomMass->_pSysData[listdata[i].offset +
                                                                          1];
        case 1:
          pSolventAtomID[solventMolecules].x = listdata[i].offset;
          pSolventMass1[solventMolecules]    = gpu->pbAtomMass->_pSysData[listdata[i].offset];
          totalMass                         += gpu->pbAtomMass->_pSysData[listdata[i].offset];
      }
      pSolventInvMass[solventMolecules]    = 1.0 / totalMass;
      solventMolecules++;
    }
    else {
      double totalMass = 0.0;
      for (int j = 0; j < listdata[i].cnt; j++) {
        pSoluteAtomID[soluteAtoms + j]         = listdata[i].offset + j;
        pSoluteAtomMoleculeID[soluteAtoms + j] = soluteMolecules;
        pSoluteAtomMass[soluteAtoms + j]       =
          gpu->pbAtomMass->_pSysData[listdata[i].offset + j];
        totalMass                             +=
          gpu->pbAtomMass->_pSysData[listdata[i].offset + j];
      }
      int offset = ((listdata[i].cnt + (GRID - 1)) >> GRID_BITS) << GRID_BITS;
      for (int j = listdata[i].cnt; j < offset; j++) {
        pSoluteAtomID[soluteAtoms + j] = -1;
        pSoluteAtomMoleculeID[soluteAtoms + j] = soluteMolecules;
        pSoluteAtomMass[soluteAtoms + j]       = 0.0;
      }
      soluteAtoms                             += offset;
      pSoluteInvMass[soluteMolecules]          = 1.0 / totalMass;
      pAFEMolType[soluteMolecules]             = ti_mol_type[i];
      pAFEMolPartner[soluteMolecules]          = ti_sc_partner[i];
      soluteMolecules++;
    }
  }

  // Upload data
  gpu->pbSoluteAtomID->Upload();
  gpu->pbSoluteAtomMass->Upload();
  gpu->pbSolute->Upload();
  gpu->pbUllSolute->Upload();
  gpu->pbSolventAtomID->Upload();
  gpu->pbSolvent->Upload();

  // Set up constant pointers
  gpu->sim.soluteMolecules       = soluteMolecules;
  gpu->sim.soluteMoleculeStride  = soluteMoleculeStride;
  gpu->sim.soluteAtoms           = soluteAtoms;
  gpu->sim.solventMolecules      = solventMolecules;
  gpu->sim.solventMoleculeStride = solventMoleculeStride;
  gpu->sim.pSoluteCOMX           = gpu->pbSolute->_pDevData;
  gpu->sim.pSoluteCOMY           = gpu->pbSolute->_pDevData + soluteMoleculeStride;
  gpu->sim.pSoluteCOMZ           = gpu->pbSolute->_pDevData + soluteMoleculeStride * 2;
  gpu->sim.pSoluteInvMass        = gpu->pbSolute->_pDevData + soluteMoleculeStride * 3;
  gpu->sim.pSoluteDeltaCOMX      = gpu->pbSolute->_pDevData + soluteMoleculeStride * 4;
  gpu->sim.pSoluteDeltaCOMY      = gpu->pbSolute->_pDevData + soluteMoleculeStride * 5;
  gpu->sim.pSoluteDeltaCOMZ      = gpu->pbSolute->_pDevData + soluteMoleculeStride * 6;
  gpu->sim.pSoluteUllCOMX        = gpu->pbUllSolute->_pDevData;
  gpu->sim.pSoluteUllCOMY        = gpu->pbUllSolute->_pDevData + soluteMoleculeStride;
  gpu->sim.pSoluteUllCOMZ        = gpu->pbUllSolute->_pDevData + soluteMoleculeStride * 2;
  gpu->sim.pSoluteUllEKCOMX      = gpu->pbUllSolute->_pDevData + soluteMoleculeStride * 3;
  gpu->sim.pSoluteUllEKCOMY      = gpu->pbUllSolute->_pDevData + soluteMoleculeStride * 4;
  gpu->sim.pSoluteUllEKCOMZ      = gpu->pbUllSolute->_pDevData + soluteMoleculeStride * 5;
  gpu->sim.pSoluteAtomID         = gpu->pbSoluteAtomID->_pDevData;
  gpu->sim.pSoluteAtomMoleculeID = gpu->pbSoluteAtomID->_pDevData + soluteAtoms;
  gpu->sim.pSoluteAtomMass       = gpu->pbSoluteAtomMass->_pDevData;
  gpu->sim.pSolventAtomMass1     = gpu->pbSolvent->_pDevData;
  gpu->sim.pSolventAtomMass2     = gpu->pbSolvent->_pDevData + solventMoleculeStride;
  gpu->sim.pSolventAtomMass3     = gpu->pbSolvent->_pDevData + solventMoleculeStride * 2;
  gpu->sim.pSolventAtomMass4     = gpu->pbSolvent->_pDevData + solventMoleculeStride * 3;
  gpu->sim.pSolventCOMX          = gpu->pbSolvent->_pDevData + solventMoleculeStride * 4;
  gpu->sim.pSolventCOMY          = gpu->pbSolvent->_pDevData + solventMoleculeStride * 5;
  gpu->sim.pSolventCOMZ          = gpu->pbSolvent->_pDevData + solventMoleculeStride * 6;
  gpu->sim.pSolventInvMass       = gpu->pbSolvent->_pDevData + solventMoleculeStride * 7;
  gpu->sim.pSolventAtomID        = gpu->pbSolventAtomID->_pDevData;
  gpu->sim.pAFEMolType           = gpu->pbAFEMolType->_pDevData;
  gpu->sim.pAFEMolPartner        = gpu->pbAFEMolPartner->_pDevData;
  gpuCopyConstants();

#ifdef GVERBOSE
  printf("Found %d molecules, %d solvent, %d solute, and %d effective solute atoms\n",
         *molecules, solventMolecules, soluteMolecules, soluteAtoms);
#endif
}

//---------------------------------------------------------------------------------------------
// gpu_ntp_setup_: setup for constant pressure simulations
//---------------------------------------------------------------------------------------------
extern "C" void gpu_ntp_setup_()
{
  PRINTMETHOD("gpu_pme_ntp_setup");

  gpu->pbNTPData      = new GpuBuffer<NTPData>(1);
  NTPData* pNTPData   = gpu->pbNTPData->_pSysData;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pNTPData->ucell[i*3 + j] = gpu->sim.ucell[i][j];
      pNTPData->ucellf[i*3 + j] = gpu->sim.ucellf[i][j];
      pNTPData->recip[i*3 + j] = gpu->sim.recip[i][j];
      pNTPData->recipf[i*3 + j] = gpu->sim.recipf[i][j];
      pNTPData->last_recip[i*3 + j] = gpu->sim.recip[i][j];
    }
  }
  pNTPData->one_half_nonbond_skin_squared = gpu->sim.one_half_nonbond_skin_squared;
  pNTPData->cutPlusSkin2 = (gpu->sim.cut + gpu->sim.nonbond_skin) *
                           (gpu->sim.cut + gpu->sim.nonbond_skin);
  gpu->sim.pNTPData = gpu->pbNTPData->_pDevData;
  gpu->pbNTPData->Upload();
  gpuCopyConstants();

  // Calculate initial COM
  kCalculateCOM(gpu);
  kReduceSoluteCOM(gpu);
}

//---------------------------------------------------------------------------------------------
// Here begin functions for erfc() spline setup
//---------------------------------------------------------------------------------------------

#ifndef use_DPFP
//---------------------------------------------------------------------------------------------
// GenSplines: generate a batch of splines for evenly spaced intervals over a given range.  The
//             splines take the form
//
//             f(r2) = A*r2 + B + C/r2 + D/(r2*r2)
//
// Arguments:
//   lbound:  lower limit of the range
//   hbound:  upper limit of the range
//   nspl:    the number of spline intervals
//   alpha:   the Ewald coefficient, 0.5/(Gaussian sigma for charge spreading)
//   coef:    pointer for dumping results into the coefficients array
//   offset:  offset for placing spline coefficients in the array coef.  This will determine
//            the behavior of the data fitting target: 0 for standard d/dr [erfc(r)/r], 1 for
//            d/dr [erfc(r)/r - 1/r]
//---------------------------------------------------------------------------------------------
static void GenSplines(double lbound, double hbound, int nspl, double alpha, PMEFloat4* coef,
                       int offset)
{
  int i, j, imin;
  double lint, gspc, vspc, r2;
  double* b;
  double* bsave;
  dmat A;

  // Allocate for the matrix equation
  A = CreateDmat(21, 4);
  b = (double*)malloc(21 * sizeof(double));
  bsave = (double*)malloc(21 * sizeof(double));
  gspc = (hbound - lbound) / nspl;
  imin = (lbound < 1.0e-8) ? 1 : 0;
  for (i = imin; i < nspl; i++) {
    lint = lbound + i*gspc;

    // Initialize the matrix equation
    SetDVec(A.data, A.row * A.col, 0.0);
    SetDVec(b, A.row, 0.0);

    // Take 21 samples over the interval, inclusive of the limits
    vspc = 0.05*gspc;
    for (j = 0; j <= 20; j++) {
      r2 = lint + j*vspc;
      A.map[j][0] = r2;
      A.map[j][1] = 1.0;
      A.map[j][2] = 1.0 / r2;
      A.map[j][3] = 1.0 / (r2 * r2);
      b[j] = (erfc(alpha * sqrt(r2))/sqrt(r2) +
              (2.0 * alpha / sqrt(PI_VAL)) * exp(-alpha * alpha * r2)) / r2;
      if (offset == 1) {
        b[j] -= 1.0 / (r2 * sqrt(r2));
      }
      bsave[j] = b[j];
    }

    // Solve the matrix equation
    AxbQRRxc(A, b);
    BackSub(A, b);

    // Commit the result (and convert double --> float)
    coef[2*i + offset].x = b[0];
    coef[2*i + offset].y = b[1];
    coef[2*i + offset].z = b[2];
    coef[2*i + offset].w = b[3];
  }

  // Fill in the spline containing the origin with an extension of the
  // next nearest spline.  Two non-bonded point charges should never get this
  // close, but if they do the spline will apply a finite force.
  if (imin == 1) {
    coef[offset].x = coef[2+offset].x;
    coef[offset].y = coef[2+offset].y + coef[2+offset].z/gspc + coef[2+offset].w/(gspc*gspc);
    coef[offset].z = 0.0;
    coef[offset].z = 0.0;
  }

  // Free allocated memory
  DestroyDmat(&A);
  free(b);
  free(bsave);
}

//---------------------------------------------------------------------------------------------
// MakeCoeffsTable: build the coefficients table for interpolating d/dr [erfc(r)/r].  The
//                  table is indexed by another ErfcIndexTable, which itself is indexed by r2.
//                  The splines are given by f(r2) = A*r2 + B + C/r2 + D/(r2*r2).  Because
//                  1/r2 needs to be calculated anyway, this spine is very easy to handle and
//                  fits the target function well.
//
// Arguments:
//   gpu:       overarching struct storing information about the simulation, including atom
//              properties and the energy function
//---------------------------------------------------------------------------------------------
static void MakeCoeffsTable()
{
  int i, j, splmin;
  double range, rint, yfac, yfacEX;
  PMEFloat4 coef;

  // Splines for each interval, with decreasing density
  // as the function becomes easier to interpolate.
  splmin = 117;
  for (i = splmin; i < 144; i++) {
    range = 1.0 / pow(2.0, 127-i);
    rint = range / 32.0;
    GenSplines(range, range + 32.0*rint, 32, gpu->sim.ew_coeff,
               &gpu->pbErfcCoeffsTable->_pSysData[64*i], 0);
    GenSplines(range, range + 32.0*rint, 32, gpu->sim.ew_coeff,
               &gpu->pbErfcCoeffsTable->_pSysData[64*i], 1);
    if (i == splmin) {
      coef = gpu->pbErfcCoeffsTable->_pSysData[64*i];
      yfac = (coef.x * range) + coef.y + (coef.z / range) + (coef.w / (range * range));
      coef = gpu->pbErfcCoeffsTable->_pSysData[64*i + 1];
      yfacEX = (coef.x * range) + coef.y + (coef.z / range) + (coef.w / (range * range));
    }
  }
  coef.x = 0.0;
  coef.z = 0.0;
  coef.w = 0.0;
  for (i = 0; i < splmin; i++) {
    for (j = 0; j < 32; j++) {
      coef.y = yfac;
      gpu->pbErfcCoeffsTable->_pSysData[64*i + 2*j] = coef;
      coef.y = yfacEX;
      gpu->pbErfcCoeffsTable->_pSysData[64*i + 2*j + 1] = coef;
    }
  }

  // Special case: the zeroth array element is a redirect for anything that overflows the
  // cutoff.  Set its coefficients to all zeros.
  coef.y = 0.0;
  gpu->pbErfcCoeffsTable->_pSysData[0] = coef;
  gpu->pbErfcCoeffsTable->_pSysData[1] = coef;

  // Assign device pointer and upload
  gpu->sim.pErfcCoeffsTable = gpu->pbErfcCoeffsTable->_pDevData;
}

//---------------------------------------------------------------------------------------------
// CalculateErfcSpline: driver function for creating a wicked d/dr erfc(r)/r spline--for
//                      32-bit floating point numbers the accuracy rivals analytic
//                      computation.  The texture memory requirement is only 4kb for a 10A
//                      cutoff and increases as log2(cutoff squared).  This is useful in force
//                      computations, but for energies erfc() itself is needed and should be
//                      computed directly using fasterfc() in this library.
//
// Arguments:
//   gpu:       overarching struct storing information about the simulation, including atom
//              properties and the energy function
//---------------------------------------------------------------------------------------------
static void CalculateErfcSpline()
{
  // Allocate the coefficients table with room for cutoffs up to 360A
  int nset = 144 * 32;
  gpu->pbErfcCoeffsTable = new GpuBuffer<PMEFloat4>(2 * nset, true, false);

  // Fill and upload the coefficients table
  MakeCoeffsTable();
  gpu->pbErfcCoeffsTable->Upload();
  gpuCopyConstants();

  // Adjust the coefficients if possible
  kAdjustCoeffsTable(gpu);

#if !defined(use_DPFP)
  cudaError_t status;
  if (gpu->sim.texErfcCoeffsTable != 0) {
    status = cudaDestroyTextureObject(gpu->sim.texErfcCoeffsTable);
    RTERROR(status, "cudaDestroyTextureObject gpu->sim.texErfcCoeffsTable failed");
  }
  cudaResourceDesc resDesc;
  cudaTextureDesc texDesc;
  memset(&resDesc, 0, sizeof(resDesc));
  resDesc.resType = cudaResourceTypeLinear;
  resDesc.res.linear.devPtr = gpu->pbErfcCoeffsTable->_pDevData;
  resDesc.res.linear.sizeInBytes = gpu->pbErfcCoeffsTable->_length*sizeof(PMEFloat4);
  resDesc.res.linear.desc.f = cudaChannelFormatKindFloat;
  resDesc.res.linear.desc.x = 32;
  resDesc.res.linear.desc.y = 32;
  resDesc.res.linear.desc.z = 32;
  resDesc.res.linear.desc.w = 32;
  memset(&texDesc, 0, sizeof(texDesc));
  texDesc.normalizedCoords = 0;
  texDesc.filterMode = cudaFilterModePoint;
  texDesc.addressMode[0] = cudaAddressModeClamp;
  texDesc.readMode = cudaReadModeElementType;
  status = cudaCreateTextureObject(&(gpu->sim.texErfcCoeffsTable), &resDesc, &texDesc, NULL);
  RTERROR(status, "cudaCreateTextureObject gpu->sim.texErfcCoeffsTable failed");
#endif
}
#endif // use_DPFP not defined
// Here end declarations of functions for erfc() spline setup

//---------------------------------------------------------------------------------------------
// gpu_pme_alltasks_setup_: setup for all tasks involving (Smooth) Particle Mesh Ewald.  Note
//                          where this is called by the Fortran code, in pme_alltasks_setup,
//                          at the end of pme_alltasks_setup, in which bonded interactions
//                          and PME parameters (i.e. grid dimensions) have all been laid out.
//                          This is a good place to stash extra preparatory work for the GPU,
//                          as most of the CPU-side data has been established.
//
// Arguments:
//   nfft(1,2,3):    the FFT grid dimensions
//   bspl_order:     particle-mesh interpolation order
//   prefac(1,2,3):  lists of prefactors for scalar sum computation that can be pre-tabulated
//                   (in the case of orthorhombic unit cells)
//   ew_coeff:       the Ewald coefficient, 1 / (2 x Gaussian sigma)
//   ips:            flag to indicate that Isotropic Periodic Sums are in effect, not PME
//   fswitch:        the non-bonded cutoff
//   ef(x,y,z):      time-based electric field charges (these are NOT atom partial charges--
//                   rather, things that will exert a spooky force on all particles in the
//                   simulation as it progresses)
//   efn:            flag to have the electric field normalized (only works in an orthorhombic
//                   unit cell--the CPU code still needs a trap to ensure this)
//   efphase:        phase of the electric field (it waxes and wanes in a cosine waveform as a
//                   function of time)
//   effreq:         frequency of electric field oscillation
//---------------------------------------------------------------------------------------------
extern "C" void gpu_pme_alltasks_setup_(int* nfft1, int* nfft2, int* nfft3, int* bspl_order,
                                        double* prefac1, double* prefac2, double* prefac3,
                                        double* ew_coeff, int* ips, double* fswitch,
                                        double* efx, double* efy, double* efz, int* efn,
                                        double* efphase, double* effreq, double *vdw_cutoff)
{
  PRINTMETHOD("gpu_pme_alltasks_setup");

  // Delete any existing PME state
  delete gpu->pbPrefac;
  delete gpu->pbFract;

  // Determine whether running IPS or PME
  if (*ips != 0) {
    gpu->sim.bIPSActive = true;
  }
  else {
    gpu->sim.bIPSActive = false;
  }

  // Allocate GPU data
  int n1                 = ((*nfft1 + 1) + PADDING) & PADDING_MASK;
  int n2                 = ((*nfft2 + 1) + PADDING) & PADDING_MASK;
  int n3                 = ((*nfft3 + 1) + PADDING) & PADDING_MASK;
  gpu->sim.efx           = *efx;
  gpu->sim.efy           = *efy;
  gpu->sim.efz           = *efz;
  gpu->sim.efn           = *efn;
  gpu->sim.efphase       = *efphase;
  gpu->sim.effreq        = *effreq;
  gpu->sim.fswitch       = *fswitch;
  gpu->sim.fswitch2      = gpu->sim.fswitch * gpu->sim.fswitch;
  gpu->sim.fswitch3      = gpu->sim.fswitch * gpu->sim.fswitch2;
  gpu->sim.fswitch6      = gpu->sim.fswitch3 * gpu->sim.fswitch3;
  gpu->sim.n2Offset      = n1;
  gpu->sim.n3Offset      = n1 + n2;
  gpu->sim.nSum          = n1 + n2 + n3;
  gpu->sim.nfft1         = *nfft1;
  gpu->sim.nfft2         = *nfft2;
  gpu->sim.nfft3         = *nfft3;
  gpu->sim.dnfft1        = (double)(*nfft1);
  gpu->sim.dnfft2        = (double)(*nfft2);
  gpu->sim.dnfft3        = (double)(*nfft3);
  gpu->sim.pmeOrder      = *bspl_order;
  gpu->sim.orderMinusOne = *bspl_order - 1;
  gpu->sim.fft_x_dim     = *nfft1 / 2 + 1;
  gpu->sim.fft_y_dim     = *nfft2;
  gpu->sim.fft_z_dim     = *nfft3;
  gpu->sim.nf1           = *nfft1 / 2;
  gpu->sim.nfft1xnfft2   = gpu->sim.nfft1 * gpu->sim.nfft2;
  gpu->sim.nfft1xnfft2xnfft3 = gpu->sim.nfft1 * gpu->sim.nfft2 * gpu->sim.nfft3;
  gpu->sim.fft_y_dim_times_x_dim = gpu->sim.fft_y_dim * gpu->sim.fft_x_dim;
  gpu->sim.fft_x_y_z_dim         = gpu->sim.fft_y_dim_times_x_dim * gpu->sim.fft_z_dim;
  gpu->sim.fft_x_y_z_quarter_dim = gpu->sim.fft_x_y_z_dim / 4;
  gpu->sim.fft_quarter_z_dim_m1  = (*nfft3 / 4) - 1;
  if (2 * gpu->sim.nf1 < *nfft1) {
    gpu->sim.nf1++;
  }
  gpu->sim.nf2 = *nfft2 / 2;
  if (2 * gpu->sim.nf2 < *nfft2) {
    gpu->sim.nf2++;
  }
  gpu->sim.nf3 = *nfft3 / 2;
  if (2 * gpu->sim.nf3 < *nfft3) {
    gpu->sim.nf3++;
  }

  // Set up PME charge buffer
  gpu->sim.ew_coeff   = *ew_coeff;
  gpu->sim.ew_coeffSP = *ew_coeff;
  gpu->sim.ew_coeff2  = (*ew_coeff) * (*ew_coeff);
#ifndef use_DPFP
  SetkCalculatePMENonbondEnergyERFC(*ew_coeff);
  CalculateErfcSpline();
#endif
  gpu->sim.negTwoEw_coeffRsqrtPI = -2.0 * gpu->sim.ew_coeffSP / sqrt(PI);
  gpu->sim.fac       = (PI * PI) / ((*ew_coeff) * (*ew_coeff));
  gpu->sim.fac2      = 2.0 * gpu->sim.fac;
  gpu->sim.XYZStride = ((2 * gpu->sim.fft_x_dim * gpu->sim.fft_y_dim * gpu->sim.fft_z_dim) +
                        31) & 0xffffffe0;
  gpu->pbPrefac      = new GpuBuffer<PMEFloat>(n1 + n2 + n3);
  gpu->pbFract       = new GpuBuffer<PMEFloat>(gpu->sim.stride3);

  // Copy Prefac data
  for (int i = 0; i < *nfft1; i++) {
    gpu->pbPrefac->_pSysData[i + 1] = prefac1[i];
  }
  for (int i = 0; i < *nfft2; i++) {
    gpu->pbPrefac->_pSysData[i + 1 + n1] = prefac2[i];
  }
  for (int i = 0; i < *nfft3; i++) {
    gpu->pbPrefac->_pSysData[i + 1 + n1 + n2] = prefac3[i];
  }
  gpu->pbPrefac->Upload();

  // Allocate PME buffers if needed
  if (gpu->ntb != 0) {
#ifdef use_DPFP
    gpu->pblliXYZ_q    = new GpuBuffer<long long int>(gpu->sim.XYZStride);
#else
    gpu->pblliXYZ_q    = new GpuBuffer<int>(gpu->sim.XYZStride);
#endif
    gpu->sim.plliXYZ_q = gpu->pblliXYZ_q->_pDevData;
    gpu->pbXYZ_q       = new GpuBuffer<PMEFloat>(gpu->sim.XYZStride);
    gpu->sim.pXYZ_q    = gpu->pbXYZ_q->_pDevData;
    {
      cudaError_t status;
      if (gpu->sim.texXYZ_q != 0) {
        status = cudaDestroyTextureObject(gpu->sim.texXYZ_q);
        RTERROR(status, "cudaDestroyTextureObject gpu->sim.texXYZ_q failed");
      }
      cudaResourceDesc resDesc;
      cudaTextureDesc texDesc;
      memset(&resDesc, 0, sizeof(resDesc));
      resDesc.resType = cudaResourceTypeLinear;
      resDesc.res.linear.devPtr = gpu->sim.pXYZ_q;
#ifdef use_DPFP
      resDesc.res.linear.sizeInBytes = gpu->sim.nfft1 * gpu->sim.nfft2 * gpu->sim.nfft3 * sizeof(int2);
      resDesc.res.linear.desc.f = cudaChannelFormatKindSigned;
      resDesc.res.linear.desc.y = 32;
#else
      resDesc.res.linear.sizeInBytes = gpu->sim.nfft1 * gpu->sim.nfft2 * gpu->sim.nfft3 * sizeof(float);
      resDesc.res.linear.desc.f = cudaChannelFormatKindFloat;
      resDesc.res.linear.desc.y = 0;
#endif
      resDesc.res.linear.desc.x = 32;
      resDesc.res.linear.desc.z = 0;
      resDesc.res.linear.desc.w = 0;
      memset(&texDesc, 0, sizeof(texDesc));
      texDesc.normalizedCoords = 0;
      texDesc.filterMode = cudaFilterModePoint;
      texDesc.addressMode[0] = cudaAddressModeClamp;
      texDesc.readMode = cudaReadModeElementType;
      status = cudaCreateTextureObject(&(gpu->sim.texXYZ_q), &resDesc, &texDesc, NULL);
      RTERROR(status, "cudaCreateTextureObject gpu->sim.texXYZ_q failed");
    }
    gpu->pbXYZ_qt      = new GpuBuffer<PMEComplex>(gpu->sim.fft_x_dim * gpu->sim.fft_y_dim *
                                                   gpu->sim.fft_z_dim);
    gpu->sim.pXYZ_qt   = gpu->pbXYZ_qt->_pDevData;
  }

  // Set up pointers
  gpu->sim.pPrefac1 = gpu->pbPrefac->_pDevData;
  gpu->sim.pPrefac2 = gpu->pbPrefac->_pDevData + n1;
  gpu->sim.pPrefac3 = gpu->pbPrefac->_pDevData + n1 + n2;
  gpu->sim.pFractX  = gpu->pbFract->_pDevData;
  gpu->sim.pFractY  = gpu->pbFract->_pDevData + gpu->sim.stride;
  gpu->sim.pFractZ  = gpu->pbFract->_pDevData + gpu->sim.stride2;

  // Set up FFT plans
#ifdef use_DPFP
  cufftPlan3d(&(gpu->forwardPlan), gpu->sim.nfft3, gpu->sim.nfft2, gpu->sim.nfft1, CUFFT_D2Z);
  cufftPlan3d(&(gpu->backwardPlan),  gpu->sim.nfft3, gpu->sim.nfft2, gpu->sim.nfft1,
              CUFFT_Z2D);
#else
  cufftPlan3d(&(gpu->forwardPlan), gpu->sim.nfft3, gpu->sim.nfft2, gpu->sim.nfft1, CUFFT_R2C);
  cufftPlan3d(&(gpu->backwardPlan),  gpu->sim.nfft3, gpu->sim.nfft2, gpu->sim.nfft1,
              CUFFT_C2R);
#endif

#ifdef VKFFT
  gpu->inputBufferSize = gpu->sim.nfft1 * gpu->sim.nfft2 * gpu->sim.nfft3 * sizeof(PMEFloat);
  gpu->outputBufferSize = gpu->sim.fft_x_dim * gpu->sim.fft_y_dim * gpu->sim.fft_z_dim *
                          sizeof(PMEComplex);

  VkFFTConfiguration configuration = {};
  configuration.performR2C = true;
  configuration.device = &gpu->gpu_device_id;
#  ifdef use_DPFP
  configuration.doublePrecision = true;
#  endif

  configuration.FFTdim = 3;
  configuration.size[0] = gpu->sim.nfft1;
  configuration.size[1] = gpu->sim.nfft2;
  configuration.size[2] = gpu->sim.nfft3;

  configuration.inverseReturnToInputBuffer = true;
  configuration.isInputFormatted = true;
  configuration.inputBufferSize = &gpu->inputBufferSize;
  configuration.inputBuffer = (void**)&gpu->sim.pXYZ_q;
  configuration.inputBufferStride[0] = gpu->sim.nfft1;
  configuration.inputBufferStride[1] = configuration.inputBufferStride[0] * gpu->sim.nfft2;
  configuration.inputBufferStride[2] = configuration.inputBufferStride[1] * gpu->sim.nfft3;

  configuration.bufferSize = &gpu->outputBufferSize;
  configuration.buffer = (void**)&gpu->sim.pXYZ_qt;
  configuration.bufferStride[0] = gpu->sim.fft_x_dim;
  configuration.bufferStride[1] = configuration.bufferStride[0] * gpu->sim.fft_y_dim;
  configuration.bufferStride[2] = configuration.bufferStride[1] * gpu->sim.fft_z_dim;

  gpu->appR2C = new VkFFTApplication();
  VkFFTResult fftResult = initializeVkFFT(gpu->appR2C, configuration);
  if (fftResult != VKFFT_SUCCESS) {
    printf("VkFFT error: %d\n", fftResult);
    cudaDeviceReset();
    exit(-1);
  }
#endif

  // Upload information to the device
  gpuCopyConstants();

  return;
}

//---------------------------------------------------------------------------------------------
// gpu_remap_bonded_setup_: function to encapsulate a lot of effort at the end of
//                          pme_alltasks_setup in the Fortran code.  This is intended to take
//                          place after most other work has been done to lay out data for the
//                          GPU, and builds on a great deal of that data.
//
// Arguments:
//   numex:     array holding the number of exclusions for each atom in the system
//              (taken directly from the Fortran code)
//   natex:     the atomic indices of all exclusions, packed together one after another (see
//              the Amber file formats desciption of the topology for a detailed description)
//   ico:       the non-bonded parameter type index of each atom
//   cn1:       Lennard-Jones A parameters
//   cn2:       Lennard-Jones B parameters
//---------------------------------------------------------------------------------------------
extern "C" void gpu_remap_bonded_setup_(int ico[], double cn1[], double cn2[], int numex[],
                                        int natex[], int *ntypes)
{
  int nbwunits;
  bondwork* bwunits;

  // Lay out bonded interaction blocks
  bwunits = AssembleBondWorkUnits(gpu, &nbwunits);
  gpu->sim.bondWorkUnits = nbwunits;
#ifdef AMBER_PLATFORM_AMD
  gpu->sim.clearQBWorkUnits = ((gpu->sim.XYZStride + CHARGE_BUFFER_STRIDE - 1) /
                               CHARGE_BUFFER_STRIDE);
#else
  gpu->sim.clearQBWorkUnits = nbwunits + ((gpu->sim.XYZStride + CHARGE_BUFFER_STRIDE - 1) /
                                          CHARGE_BUFFER_STRIDE);
#endif
  int ncluster = (nbwunits + 65535) / 65536;
  ncluster += (ncluster == 0);
  gpu->bondWorkBlocks = gpu->blocks * BOND_WORK_UNIT_BLOCKS_MULTIPLIER;
  gpu->bondWorkUnitRecord = bwunits;

  // Allocate space for the bond work units.  Set pointers to operate on the GPU.
  bwalloc bwdims = CalcBondBlockSize(gpu, bwunits, nbwunits);
  gpu->pbBondWorkUnitUINT    = new GpuBuffer<unsigned int>(bwdims.nUint, true, false);
  gpu->pbBondWorkUnitDBL2    = new GpuBuffer<PMEDouble2>(bwdims.nDbl2, true, false);
  gpu->pbBondWorkUnitPFLOAT  = new GpuBuffer<PMEFloat>(bwdims.nPFloat, true, false);
  gpu->pbBondWorkUnitPFLOAT2 = new GpuBuffer<PMEFloat2>(bwdims.nPFloat2, true, false);
  int insrUintPos    = 0;
  int bondUintPos    = insrUintPos +   (3 * BOND_WORK_UNIT_THREADS_PER_BLOCK) * nbwunits;
  int anglUintPos    = bondUintPos +   (GRID * bwdims.nbondwarps);
  int diheUintPos    = anglUintPos +   (GRID * bwdims.nanglwarps);
  int cmapUintPos    = diheUintPos +   (GRID * bwdims.ndihewarps);
  int qqxcUintPos    = cmapUintPos + 2*(GRID * bwdims.ncmapwarps);
  int nb14UintPos    = qqxcUintPos +   (GRID * bwdims.nqqxcwarps);
  int nmr2UintPos    = nb14UintPos +   (GRID * bwdims.nnb14warps);
  int nmr3UintPos    = nmr2UintPos + 4*(GRID * bwdims.nnmr2warps);
  int nmr4UintPos    = nmr3UintPos + 4*(GRID * bwdims.nnmr3warps);
  int ureyUintPos    = nmr4UintPos + 4*(GRID * bwdims.nnmr4warps);
  int cimpUintPos    = ureyUintPos +   (GRID * bwdims.nureywarps);
  int cnstUintPos    = cimpUintPos +   (GRID * bwdims.ncimpwarps);
  int cnstUpdatePos  = cnstUintPos +   (GRID * bwdims.ncnstwarps);
  int bondDbl2Pos    = 0;
  int anglDbl2Pos    = bondDbl2Pos +   (GRID * bwdims.nbondwarps);
  int nmr2Dbl2Pos    = anglDbl2Pos +   (GRID * bwdims.nanglwarps);
  int nmr3Dbl2Pos    = nmr2Dbl2Pos + 6*(GRID * bwdims.nnmr2warps);
  int nmr4Dbl2Pos    = nmr3Dbl2Pos + 6*(GRID * bwdims.nnmr3warps);
  int ureyDbl2Pos    = nmr4Dbl2Pos + 6*(GRID * bwdims.nnmr4warps);
  int cimpDbl2Pos    = ureyDbl2Pos +   (GRID * bwdims.nureywarps);
  int cnstDbl2Pos    = cimpDbl2Pos +   (GRID * bwdims.ncimpwarps);
  int qPmefPos       = 0;
  int dihePmefPos    = nbwunits * BOND_WORK_UNIT_THREADS_PER_BLOCK;
  int nb14PmefPos    = dihePmefPos +    (GRID * bwdims.ndihewarps);
  int dihePmef2Pos   = 0;
  int nb14Pmef2Pos   = dihePmef2Pos + 2*(GRID * bwdims.ndihewarps);
  int bondStatusPos, anglStatusPos, diheStatusPos, cmapStatusPos, qqxcStatusPos, nb14StatusPos,
      nmr2StatusPos, nmr3StatusPos, nmr4StatusPos, ureyStatusPos, cimpStatusPos, cnstStatusPos;
  if (gpu->sim.ti_mode > 0) {
    bondStatusPos = cnstUpdatePos + (GRID * bwdims.ncnstwarps);
    anglStatusPos = bondStatusPos + (GRID * bwdims.nbondwarps);
    diheStatusPos = anglStatusPos + (GRID * bwdims.nanglwarps);
    cmapStatusPos = diheStatusPos + (GRID * bwdims.ndihewarps);
    qqxcStatusPos = cmapStatusPos + (GRID * bwdims.ncmapwarps);
    nb14StatusPos = qqxcStatusPos + (GRID * bwdims.nqqxcwarps);
    nmr2StatusPos = nb14StatusPos + (GRID * bwdims.nnb14warps);
    nmr3StatusPos = nmr2StatusPos + (GRID * bwdims.nnmr2warps);
    nmr4StatusPos = nmr3StatusPos + (GRID * bwdims.nnmr3warps);
    ureyStatusPos = nmr4StatusPos + (GRID * bwdims.nnmr4warps);
    cimpStatusPos = ureyStatusPos + (GRID * bwdims.nureywarps);
    cnstStatusPos = cimpStatusPos + (GRID * bwdims.ncimpwarps);
  }

  // Set pointers to data on the device--these will be used by kernels
  gpu->sim.pBwuInstructions  = &gpu->pbBondWorkUnitUINT->_pDevData[insrUintPos];
  gpu->sim.pBwuBondID        = &gpu->pbBondWorkUnitUINT->_pDevData[bondUintPos];
  gpu->sim.pBwuAnglID        = &gpu->pbBondWorkUnitUINT->_pDevData[anglUintPos];
  gpu->sim.pBwuDiheID        = &gpu->pbBondWorkUnitUINT->_pDevData[diheUintPos];
  gpu->sim.pBwuCmapID        = &gpu->pbBondWorkUnitUINT->_pDevData[cmapUintPos];
  gpu->sim.pBwuQQxcID        = &gpu->pbBondWorkUnitUINT->_pDevData[qqxcUintPos];
  gpu->sim.pBwuNB14ID        = &gpu->pbBondWorkUnitUINT->_pDevData[nb14UintPos];
  gpu->sim.pBwuNMR2ID        = &gpu->pbBondWorkUnitUINT->_pDevData[nmr2UintPos];
  gpu->sim.pBwuNMR3ID        = &gpu->pbBondWorkUnitUINT->_pDevData[nmr3UintPos];
  gpu->sim.pBwuNMR4ID        = &gpu->pbBondWorkUnitUINT->_pDevData[nmr4UintPos];
  gpu->sim.pBwuUreyID        = &gpu->pbBondWorkUnitUINT->_pDevData[ureyUintPos];
  gpu->sim.pBwuCImpID        = &gpu->pbBondWorkUnitUINT->_pDevData[cimpUintPos];
  gpu->sim.pBwuCnstID        = &gpu->pbBondWorkUnitUINT->_pDevData[cnstUintPos];
  gpu->sim.pBwuCnstUpdateIdx = &gpu->pbBondWorkUnitUINT->_pDevData[cnstUpdatePos];
  gpu->sim.pBwuBond          = &gpu->pbBondWorkUnitDBL2->_pDevData[bondDbl2Pos];
  gpu->sim.pBwuAngl          = &gpu->pbBondWorkUnitDBL2->_pDevData[anglDbl2Pos];
  gpu->sim.pBwuNMR2          = &gpu->pbBondWorkUnitDBL2->_pDevData[nmr2Dbl2Pos];
  gpu->sim.pBwuNMR3          = &gpu->pbBondWorkUnitDBL2->_pDevData[nmr3Dbl2Pos];
  gpu->sim.pBwuNMR4          = &gpu->pbBondWorkUnitDBL2->_pDevData[nmr4Dbl2Pos];
  gpu->sim.pBwuUrey          = &gpu->pbBondWorkUnitDBL2->_pDevData[ureyDbl2Pos];
  gpu->sim.pBwuCImp          = &gpu->pbBondWorkUnitDBL2->_pDevData[cimpDbl2Pos];
  gpu->sim.pBwuCnst          = &gpu->pbBondWorkUnitDBL2->_pDevData[cnstDbl2Pos];
  gpu->sim.pBwuCharges       = &gpu->pbBondWorkUnitPFLOAT->_pDevData[qPmefPos];
  gpu->sim.pBwuDihe3         = &gpu->pbBondWorkUnitPFLOAT->_pDevData[dihePmefPos];
  gpu->sim.pBwuEEnb14        = &gpu->pbBondWorkUnitPFLOAT->_pDevData[nb14PmefPos];
  gpu->sim.pBwuDihe12        = &gpu->pbBondWorkUnitPFLOAT2->_pDevData[dihePmef2Pos];
  gpu->sim.pBwuLJnb14        = &gpu->pbBondWorkUnitPFLOAT2->_pDevData[nb14Pmef2Pos];
  if (gpu->sim.ti_mode > 0) {
    gpu->sim.pBwuBondStatus  = &gpu->pbBondWorkUnitUINT->_pDevData[bondStatusPos];
    gpu->sim.pBwuAnglStatus  = &gpu->pbBondWorkUnitUINT->_pDevData[anglStatusPos];
    gpu->sim.pBwuDiheStatus  = &gpu->pbBondWorkUnitUINT->_pDevData[diheStatusPos];
    gpu->sim.pBwuCmapStatus  = &gpu->pbBondWorkUnitUINT->_pDevData[cmapStatusPos];
    gpu->sim.pBwuQQxcStatus  = &gpu->pbBondWorkUnitUINT->_pDevData[qqxcStatusPos];
    gpu->sim.pBwuNB14Status  = &gpu->pbBondWorkUnitUINT->_pDevData[nb14StatusPos];
    gpu->sim.pBwuNMR2Status  = &gpu->pbBondWorkUnitUINT->_pDevData[nmr2StatusPos];
    gpu->sim.pBwuNMR3Status  = &gpu->pbBondWorkUnitUINT->_pDevData[nmr3StatusPos];
    gpu->sim.pBwuNMR4Status  = &gpu->pbBondWorkUnitUINT->_pDevData[nmr4StatusPos];
    gpu->sim.pBwuUreyStatus  = &gpu->pbBondWorkUnitUINT->_pDevData[ureyStatusPos];
    gpu->sim.pBwuCimpStatus  = &gpu->pbBondWorkUnitUINT->_pDevData[cimpStatusPos];
    gpu->sim.pBwuCnstStatus  = &gpu->pbBondWorkUnitUINT->_pDevData[cnstStatusPos];
  }

  // Record the lengths of array segments which will require updating.
  gpu->sim.BwuCnstCount = bwdims.ncnstwarps * GRID;

  // Set pointers to data on the host--these will be used in the event that host data
  // needs to change (i.e. in constant pH simulations).  Any data changed by referencing
  // these pointers will need to be uploaded to the device by taking the entire array
  // these point into.
  gpu->ptrBwuInstructions  = &gpu->pbBondWorkUnitUINT->_pSysData[insrUintPos];
  gpu->ptrBwuBondID        = &gpu->pbBondWorkUnitUINT->_pSysData[bondUintPos];
  gpu->ptrBwuAnglID        = &gpu->pbBondWorkUnitUINT->_pSysData[anglUintPos];
  gpu->ptrBwuDiheID        = &gpu->pbBondWorkUnitUINT->_pSysData[diheUintPos];
  gpu->ptrBwuCmapID        = &gpu->pbBondWorkUnitUINT->_pSysData[cmapUintPos];
  gpu->ptrBwuQQxcID        = &gpu->pbBondWorkUnitUINT->_pSysData[qqxcUintPos];
  gpu->ptrBwuNB14ID        = &gpu->pbBondWorkUnitUINT->_pSysData[nb14UintPos];
  gpu->ptrBwuNMR2ID        = &gpu->pbBondWorkUnitUINT->_pSysData[nmr2UintPos];
  gpu->ptrBwuNMR3ID        = &gpu->pbBondWorkUnitUINT->_pSysData[nmr3UintPos];
  gpu->ptrBwuNMR4ID        = &gpu->pbBondWorkUnitUINT->_pSysData[nmr4UintPos];
  gpu->ptrBwuUreyID        = &gpu->pbBondWorkUnitUINT->_pSysData[ureyUintPos];
  gpu->ptrBwuCimpID        = &gpu->pbBondWorkUnitUINT->_pSysData[cimpUintPos];
  gpu->ptrBwuCnstID        = &gpu->pbBondWorkUnitUINT->_pSysData[cnstUintPos];
  gpu->ptrBwuCnstUpdateIdx = &gpu->pbBondWorkUnitUINT->_pSysData[cnstUintPos];
  gpu->ptrBwuBond          = &gpu->pbBondWorkUnitDBL2->_pSysData[bondDbl2Pos];
  gpu->ptrBwuAngl          = &gpu->pbBondWorkUnitDBL2->_pSysData[anglDbl2Pos];
  gpu->ptrBwuNMR2          = &gpu->pbBondWorkUnitDBL2->_pSysData[nmr2Dbl2Pos];
  gpu->ptrBwuNMR3          = &gpu->pbBondWorkUnitDBL2->_pSysData[nmr3Dbl2Pos];
  gpu->ptrBwuNMR4          = &gpu->pbBondWorkUnitDBL2->_pSysData[nmr4Dbl2Pos];
  gpu->ptrBwuUrey          = &gpu->pbBondWorkUnitDBL2->_pSysData[ureyDbl2Pos];
  gpu->ptrBwuCimp          = &gpu->pbBondWorkUnitDBL2->_pSysData[cimpDbl2Pos];
  gpu->ptrBwuCnst          = &gpu->pbBondWorkUnitDBL2->_pSysData[cnstDbl2Pos];
  gpu->ptrBwuCharges       = &gpu->pbBondWorkUnitPFLOAT->_pSysData[qPmefPos];
  gpu->ptrBwuDihe3         = &gpu->pbBondWorkUnitPFLOAT->_pSysData[dihePmefPos];
  gpu->ptrBwuEEnb14        = &gpu->pbBondWorkUnitPFLOAT->_pSysData[nb14PmefPos];
  gpu->ptrBwuDihe12        = &gpu->pbBondWorkUnitPFLOAT2->_pSysData[dihePmef2Pos];
  gpu->ptrBwuLJnb14        = &gpu->pbBondWorkUnitPFLOAT2->_pSysData[nb14Pmef2Pos];
  if (gpu->sim.ti_mode > 0) {
    gpu->ptrBwuBondStatus  = &gpu->pbBondWorkUnitUINT->_pSysData[bondStatusPos];
    gpu->ptrBwuAnglStatus  = &gpu->pbBondWorkUnitUINT->_pSysData[anglStatusPos];
    gpu->ptrBwuDiheStatus  = &gpu->pbBondWorkUnitUINT->_pSysData[diheStatusPos];
    gpu->ptrBwuCmapStatus  = &gpu->pbBondWorkUnitUINT->_pSysData[cmapStatusPos];
    gpu->ptrBwuQQxcStatus  = &gpu->pbBondWorkUnitUINT->_pSysData[qqxcStatusPos];
    gpu->ptrBwuNB14Status  = &gpu->pbBondWorkUnitUINT->_pSysData[nb14StatusPos];
    gpu->ptrBwuNMR2Status  = &gpu->pbBondWorkUnitUINT->_pSysData[nmr2StatusPos];
    gpu->ptrBwuNMR3Status  = &gpu->pbBondWorkUnitUINT->_pSysData[nmr3StatusPos];
    gpu->ptrBwuNMR4Status  = &gpu->pbBondWorkUnitUINT->_pSysData[nmr4StatusPos];
    gpu->ptrBwuUreyStatus  = &gpu->pbBondWorkUnitUINT->_pSysData[ureyStatusPos];
    gpu->ptrBwuCimpStatus  = &gpu->pbBondWorkUnitUINT->_pSysData[cimpStatusPos];
    gpu->ptrBwuCnstStatus  = &gpu->pbBondWorkUnitUINT->_pSysData[cnstStatusPos];
  }

  // Set the indexing pointers to 0xff (all bits 1), as on the GPU a blaring bit string
  // (unsigned int as a whole == 0xffffffff) will indicate there is no term to compute.
  // Can't do this with zero as there are terms which involve only one unique atom and
  // atom number 0 is a thing.
  memset(gpu->pbBondWorkUnitUINT->_pSysData, 0xff, bwdims.nUint * sizeof(unsigned int));

  // Write the GPU instructions for each bond work unit
  MakeBondedUnitDirections(gpu, bwunits, nbwunits, &bwdims);

  // Upload work units
  gpu->pbBondWorkUnitUINT->Upload();
  gpu->pbBondWorkUnitDBL2->Upload();
  gpu->pbBondWorkUnitPFLOAT->Upload();
  gpu->pbBondWorkUnitPFLOAT2->Upload();

  // Upload constants and pointers to the device
  gpuCopyConstants();

  return;
}

//---------------------------------------------------------------------------------------------
// gpu_get_self: test whether the electrostatic energy of the plasma is consistent on the host
//               and device.  If not, make it consistent and reupload other constants.
//
// Arguments:
//   ee_plasma:  pairwise energy of the plasma creating the electric field
//   ene:        self energy of the plasma creating the electric field
//---------------------------------------------------------------------------------------------
extern "C" void gpu_get_self(double* ee_plasma, double* ene)
{
  if ((gpu->ee_plasma != *ee_plasma) || (*ene != gpu->self_energy)) {
    gpu->ee_plasma  = *ee_plasma;
    gpu->self_energy = *ene;
    gpuCopyConstants();
  }
}

//---------------------------------------------------------------------------------------------
// Electrostatic IPS parameters:
//---------------------------------------------------------------------------------------------
static const double aipse0 = -35.0 / 16.0;
static const double aipse1 =  35.0 / 16.0;
static const double aipse2 = -21.0 / 16.0;
static const double aipse3 =   5.0 / 16.0;
static const double pipsec =   1.0 + aipse0 + aipse1 + aipse2 + aipse3;
static const double pipse0 =   aipse0 - pipsec;
static const double bipse1 =   2.0 * aipse1;
static const double bipse2 =   4.0 * aipse2;
static const double bipse3 =   6.0 * aipse3;

//---------------------------------------------------------------------------------------------
// Dispersion IPS parameters:
//---------------------------------------------------------------------------------------------
static const double aipsvc0 =   7.0 / 16.0;
static const double aipsvc1 =   9.0 / 14.0;
static const double aipsvc2 =  -3.0 / 28.0;
static const double aipsvc3 =   6.0 / 7.0;
static const double pipsvcc =   1.0 + aipsvc0 + aipsvc1 + aipsvc2 + aipsvc3;
static const double pipsvc0 =   aipsvc0 - pipsvcc;
static const double bipsvc1 =   2.0 * aipsvc1;
static const double bipsvc2 =   4.0 * aipsvc2;
static const double bipsvc3 =   6.0 * aipsvc3;

//---------------------------------------------------------------------------------------------
// Repulsion IPS parameters:
//---------------------------------------------------------------------------------------------
static const double aipsva0 =  5.0 / 787.0;
static const double aipsva1 =  9.0 /  26.0;
static const double aipsva2 = -3.0 /  13.0;
static const double aipsva3 = 27.0 /  26.0;
static const double pipsvac =  1.0 + aipsva0 + aipsva1 + aipsva2 + aipsva3;
static const double pipsva0 =  aipsva0 - pipsvac;
static const double bipsva1 =  4.0 * aipsva1;
static const double bipsva2 =  8.0 * aipsva2;
static const double bipsva3 = 12.0 * aipsva3;

//---------------------------------------------------------------------------------------------
// gpu_ips_setup_: setup for Isotropic Periodic Sums
//
// Arguments:
//   rips:        the cutoff for IPS interactions
//
//---------------------------------------------------------------------------------------------
extern "C" void gpu_ips_setup_(double* rips, double atm_qterm[], int* ntypes, int iac[],
                               int ico[], double cn1[], double cn2[], double* eipssnb,
                               double* eipssel, double* virips)
{
  PRINTMETHOD("gpu_ips_setup");
  gpu->sim.rips = *rips;

  // Copied data
  gpu->sim.eipssnb = *eipssnb;
  gpu->sim.eipssel = *eipssel;
  gpu->sim.virips  = 0.5 * *virips;

  // Derived coefficients
  gpu->sim.rips2   = *rips * *rips;
  gpu->sim.ripsr   = 1.0 / (*rips);
  gpu->sim.rips2r  = 1.0 / (*rips * *rips);
  gpu->sim.rips6r  = pow(1.0 / *rips, 6);
  gpu->sim.rips12r = pow(1.0 / *rips, 12);

  // Premultiplied coefficients
  gpu->sim.aipse0     = aipse0 / *rips;
  gpu->sim.aipse1     = aipse1 / pow(*rips, 3);
  gpu->sim.aipse2     = aipse2 / pow(*rips, 5);
  gpu->sim.aipse3     = aipse3 / pow(*rips, 7);
  gpu->sim.bipse1     = bipse1 / pow(*rips, 3);
  gpu->sim.bipse2     = bipse2 / pow(*rips, 5);
  gpu->sim.bipse3     = bipse3 / pow(*rips, 7);
  gpu->sim.pipsec     = pipsec / *rips;

  gpu->sim.aipsvc0    = aipsvc0 / pow(*rips, 6);
  gpu->sim.aipsvc1    = aipsvc1 / pow(*rips, 8);
  gpu->sim.aipsvc2    = aipsvc2 / pow(*rips, 10);
  gpu->sim.aipsvc3    = aipsvc3 / pow(*rips, 12);
  gpu->sim.bipsvc1    = bipsvc1 / pow(*rips, 8);
  gpu->sim.bipsvc2    = bipsvc2 / pow(*rips, 10);
  gpu->sim.bipsvc3    = bipsvc3 / pow(*rips, 12);
  gpu->sim.pipsvcc    = pipsvcc / pow(*rips, 6);

  gpu->sim.aipsva0    = aipsva0 / pow(*rips, 12);
  gpu->sim.aipsva1    = aipsva1 / pow(*rips, 16);
  gpu->sim.aipsva2    = aipsva2 / pow(*rips, 20);
  gpu->sim.aipsva3    = aipsva3 / pow(*rips, 24);
  gpu->sim.bipsva1    = bipsva1 / pow(*rips, 16);
  gpu->sim.bipsva2    = bipsva2 / pow(*rips, 20);
  gpu->sim.bipsva3    = bipsva3 / pow(*rips, 24);
  gpu->sim.pipsvac    = pipsvac / pow(*rips, 12);

  // Calculate system energies
  double EEL    = 0.0;
  double ENB    = 0.0;
  double ripsr  = 1.0 / *rips;
  double rips6r = pow(ripsr, 6);
  for (int i = 0; i < gpu->sim.atoms; i++) {
    double qi   = atm_qterm[i];
    double qiqj = qi * qi;
    EEL         += 0.5 * qiqj * pipse0 * ripsr;
    int j       = iac[i] - 1;
    int nbtype  = ico[*ntypes * j + j] - 1;
    double A    = cn1[nbtype];
    double B    = cn2[nbtype];
    ENB         += 0.5 * (A * pipsva0 * rips6r - B * pipsvc0) * rips6r;
  }
  gpu->sim.EIPSNB = ENB;
  gpu->sim.EIPSEL = EEL;
  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gpu_ips_update_: update critical values (non-bonded sum, electrostatic sum, virial) for
//                  Isotropic Periodic sums
//
// Arguments:
//   eipssnb:    non-bonded IPS energy
//   eipssel:    electrostatic IPS energy
//   virips:     IPS virial
//---------------------------------------------------------------------------------------------
extern "C" void gpu_ips_update_(double* eipssnb, double* eipssel, double* virips)
{
  PRINTMETHOD("gpu_ips_update");
  gpu->sim.eipssnb = *eipssnb;
  gpu->sim.eipssel = *eipssel;
  gpu->sim.virips  = 0.5 * *virips;
}

//---------------------------------------------------------------------------------------------
// gpu_get_grid_weights_: gets the PME grid weights for each charge by launching an
//                        encapsulating function which in turn calls one of a list of more
//                        specific functions tuned for different system sizes.
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_get_grid_weights_()
{
  PRINTMETHOD("gpu_get_grid_weights");
  kPMEGetGridWeights(gpu);
}

//---------------------------------------------------------------------------------------------
// gpu_fill_charge_grid_: fill the charge grid on the GPU.  Again, not sure where this routine
//                        gets called, if at all.
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_fill_charge_grid_()
{
  PRINTMETHOD("gpu_fill_charge_grid");
  gpu_build_neighbor_list_();
  kPMEGetGridWeights(gpu);
  kPMEClearChargeGridBuffer(gpu);
  kPMEFillChargeGridBuffer(gpu);
  kPMEReduceChargeGridBuffer(gpu);
}

//---------------------------------------------------------------------------------------------
// gpu_fft3drc_forward_: perform the forward 3D FFT, taking the charge grid into reciprocal
//                       space.
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_fft3drc_forward_()
{
  PRINTMETHOD("gpu_fft3drc_forward");
  kPMEForwardFFT(gpu);
}

//---------------------------------------------------------------------------------------------
// gpu_scalar_sumrc_: get the scalar sum once the forward FFT is complete.
//
// Arguments:
//   vol:     box volume
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_scalar_sumrc_(double* vol)
{
  PRINTMETHOD("gpu_scalar_sumrc");
  kPMEScalarSumRC(gpu, *vol);
}

//---------------------------------------------------------------------------------------------
// gpu_fft3drc_back_: perform the backward 3D FFT, obtaining the electrostatic potential in
//                    real space.
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_fft3drc_back_()
{
  PRINTMETHOD("gpu_fft3drc_back");
  kPMEBackwardFFT(gpu);
}

//---------------------------------------------------------------------------------------------
// gpu_grad_sum_: compute the gradient sum to obtain forces on each charged particle based on
//                the electrostatic potential.
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_grad_sum_()
{
  PRINTMETHOD("gpu_grad_sum");
  kPMEGradSum(gpu);
}

//---------------------------------------------------------------------------------------------
// gpu_self_: compute the self energy of the charges for Ewald calculations.
//
// Arguments:
//   ee_plasma:   electrostatic energy of the infinite lattice of charges interacting at
//                infinite range
//   ene:         self energy of each charge
//   factor:      factor for computing self contribution to dph/dlambda if iphmd=3
//   vol:         unit cell volume for computing self contribution to dph/dlambda if iphmd=3
//   ewaldcofpi:  ewaldcof / sqrt_pi for iphmd=3  
//---------------------------------------------------------------------------------------------
extern "C" void gpu_self_(double* ee_plasma, double* ene, double* factor, double* vol, double* ewaldcof)
{
  double sqrt_pi;
  if (gpu->iphmd == 3) {
    sqrt_pi = sqrt(PI);
    gpu->ee_plasma = *factor * gpu->sumq * gpu->sumq / *vol;
    gpu->self_energy = - gpu->sumq2 * *ewaldcof / sqrt_pi + gpu->ee_plasma;
    *ee_plasma = gpu->ee_plasma;
    *ene = gpu->self_energy;
  }
  else {
    gpu->ee_plasma   = *ee_plasma;
    gpu->self_energy = *ene;
  }
}

//---------------------------------------------------------------------------------------------
// gpu_vdw_correction_: long-ranged dispersion energy corrections based on the homogeneity
//                      approximation.  This is computed on teh CPU and merely passed to the
//                      GPU once in the simulation.
//
// Arguments:
//   ene:       the long-ranged dispersion energy
//---------------------------------------------------------------------------------------------
extern "C" void gpu_vdw_correction_(double* ene)
{
    gpu->vdw_recip      = *ene;
}

//---------------------------------------------------------------------------------------------
// print_virial: print the components of the virial.
//
// Arguments:
//   message:    message to print as a qualifying statement, can help to identify the origin
//               of the call to this function
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void print_virial(char* message)
{
  static PMEDouble oldvirial[3] = {0.0, 0.0, 0.0};
  gpu->pbEnergyBuffer->Download();
  PMEDouble virial[EXTENDED_ENERGY_TERMS];
  for (int i = 0; i < 3; i++) {
    unsigned long long int val = gpu->pbEnergyBuffer->_pSysData[i + VIRIAL_OFFSET];
    if (val >= 0x8000000000000000ull) {
      virial[i] = -(PMEDouble)(val ^ 0xffffffffffffffffull) / ENERGYSCALE;
    }
    else {
      virial[i] = (PMEDouble)val / ENERGYSCALE;
    }
  }
  printf("%20s %20.10f %20.10f %20.10f\n", message, virial[0] - oldvirial[0],
         virial[1] - oldvirial[1], virial[2] - oldvirial[2]);
  oldvirial[0] = virial[0];
  oldvirial[1] = virial[1];
  oldvirial[2] = virial[2];
}

//---------------------------------------------------------------------------------------------
// gpu_force_new_neighborlist_: force the creation of a new neighborlist, to test the function
//                              of that machinery in various circumstances.  This is called by
//                              REMD applications when weighing exchanges.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_force_new_neighborlist_()
{
  PRINTMETHOD("gpu_force_new_neighborlist");
  gpu->bNeedNewNeighborList = true;
}

//---------------------------------------------------------------------------------------------
// gpu_allreduce: function for reducing a piece of data to the master thread.
//
// Arguments:
//   pbBuff:   array of data to reduce
//   size:     size of pbBuff (not each element, the whole thing)
//---------------------------------------------------------------------------------------------
extern "C" void gpu_allreduce(GpuBuffer<PMEAccumulator>* pbBuff, int size)
{
#ifdef MPI /* keep an empty body when MPI is not defined */
  PRINTMETHOD("gpu_allreduce");

  // Reduce nonbond and bonded forces if running conventional NTP
  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
    kReduceForces(gpu);
  }

  if (gpu->nGpus > 1) {
#  ifdef NCCL
    if (gpu->bNCCL) {
      ncclResult_t ret = ncclAllReduce(pbBuff->_pDevData,
                                       pbBuff->_pDevData,
                                       size,
                                       ncclInt64,
                                       ncclSum,
                                       gpu->ncclComm,
                                       cudaStreamDefault);
      if (ret != ncclSuccess) {
        printf("ncclAllreduce failed : %d\n", ret);
        exit(-1);
      }
    } else
#  endif
#  ifdef AMBER_PLATFORM_AMD
    if (gpu->bP2P && gpu->bPersistentP2P) {
      kPersistentAddAccumulators(gpu, pbBuff->_pDevData, size);
    }
    else
#  endif
    if (gpu->bP2P) {
      int cycle = 0;
      for (int i = 1; i < gpu->nGpus; i *= 2) {
        int j = gpu->gpuID ^ i;
        if (j < gpu->nGpus) {
          unsigned int tsize = size * sizeof(PMEAccumulator);
          unsigned int offset = 0;
          const int chunk = 992 * 1024;
          while (tsize > 0) {
            unsigned int csize = (tsize > chunk) ? chunk : tsize;
            cudaError_t status =
              cudaMemcpyAsync((char *)(gpu->pbPeerAccumulatorList->_pSysData[j] + cycle*size) + offset,
                              (char *)(pbBuff->_pDevData) + offset, csize, cudaMemcpyDefault);
            RTERROR(status, "gpu_allreduce Memcpy failed");
            tsize -= csize;
            offset += csize;
          }
          cudaError_t status = cudaDeviceSynchronize();
          RTERROR(status, "gpu_allreduce cudaDeviceSynchronize failed");
          MPI_Barrier(gpu->comm);
          kAddAccumulators(gpu, pbBuff->_pDevData,
                           gpu->pbPeerAccumulator->_pDevData + cycle * size, size);
          LAUNCHERROR("kAddAccumulators");
        }
        cycle++;
      }

      // Make sure there's at least 1 barrier between now and re-use of the first
      // accumulation buffer
      if (gpu->nGpus <= 2) {
        cudaDeviceSynchronize();
        MPI_Barrier(gpu->comm);
      }
    }
    else {
#  ifdef MVAPICH2GDR_GPU_DIRECT_COMM
      cudaDeviceSynchronize();
      MPI_Allreduce(MPI_IN_PLACE, pbBuff->_pDevData, size, MPI_PMEACCUMULATOR, MPI_SUM,
                    gpu->comm);
      pbBuff->Download();
#  else
      pbBuff->Download();
      MPI_Allreduce(MPI_IN_PLACE, pbBuff->_pSysData, size, MPI_PMEACCUMULATOR, MPI_SUM,
                    gpu->comm);
      pbBuff->Upload();
# endif
    }
  }
#endif
}

//---------------------------------------------------------------------------------------------
// gpu_pme_ene_: compute the energy of the system for a PME simulation
//
// Arguments:
//   ewaldcof:    the Ewald coefficient, 1 / (2*(Gaussian sigma))
//   vol:         simulation volume
//   pEnergy:     struct to hold various components of the energy function
//   enmr:        energy of NMR restraints
//   virial:      the trace of the virial tensor
//   ekcmt:       kinetic energy of the system center of mass
//   nstep:       the current step number
//   dt:          time step (i.e. 0.002ps)
//   atm_crd:     passed in to support a transfer of select coordinates to the host for
//                special force computations
//   atm_frc:     passed in to allow the results of host-side force computations to be
//                uploaded to the device
//---------------------------------------------------------------------------------------------
extern "C" void gpu_pme_ene_(double* ewaldcof, double* vol, pme_pot_ene_rec* pEnergy,
                             double enmr[3], double virial[3], double ekcmt[3], int* nstep,
                             double* dt, double atm_crd[][3], double atm_frc[][3],
                             int* ineb, int* nebfreq)
{
  PRINTMETHOD("gpu_pme_ene");
  // Rebuild neighbor list
  gpu_build_neighbor_list_();

  // Clear forces
  kClearForces(gpu, gpu->sim.NLNonbondEnergyWarps);
  if (gpu->ntf != 8) {

    // Download critical coordinates for force calculations to occur on the host side
    //if (gpu->sim.nShuttle > 0) {
    //  kRetreiveSimData(gpu, atm_crd);
    //}

#ifdef MPI
    if (gpu->bCalculateLocalForces) {
#endif
    // Local energy
    kExecuteBondWorkUnits(gpu, 1);
    kCalculateNMREnergy(gpu);
    if (gpu->iphmd == 3) {
      kCalculateSelfPHMD(*vol, *ewaldcof, gpu);
    }
#ifdef MPI
    }
    else {
      kPMEClearChargeGridBuffer(gpu);
    }
#endif
    // Non-bonded energy: every process will do either reciprocal space or direct sum
    // interactions and therefore must re-image coordinates.
    kPMEGetGridWeights(gpu);
#ifdef MPI
    if (gpu->bCalculateReciprocalSum) {
#endif
    // Reciprocal space component
    kPMEFillChargeGridBuffer(gpu);
    kPMEReduceChargeGridBuffer(gpu);
    kPMEForwardFFT(gpu);
    kPMEScalarSumRCEnergy(gpu, *vol);
    kPMEBackwardFFT(gpu);
    kPMEGradSum(gpu);
#ifdef MPI
    }
#endif

    // Direct sum (direct space) energy
#ifdef MPI
    if (gpu->bCalculateDirectSum) {
#endif
    kCalculatePMENonbondEnergy(gpu);

    // Run CPU (host)-side potential and force calculations while the non-bonded
    // force kernel is running.  This is the best opportunity to mask the cost.

#ifdef MPI
    }
    if (gpu->bCalculateLocalForces) {
#endif
    // Electric Field Energy
    if (gpu->sim.efx != 0 || gpu->sim.efy != 0 || gpu->sim.efz != 0) {
      kCalculateEFieldEnergy(gpu, *nstep, *dt);
    }
#ifdef MPI
    }
#endif

    // Commute forces on extra points
    if (gpu->sim.EPs > 0) {
      kOrientForces(gpu);
    }

    // Download forces on critical atoms for additional CPU calculations
    //if (gpu->sim.nShuttle > 0) {
    //  kRetrieveSimData(gpu, atm_frc, 2);
    //}

  }
  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
#ifdef MPI
    if (gpu->gpuID == 0) {
#endif
    kCalculateCOMKineticEnergy(gpu);
    kReduceCOMKineticEnergy(gpu);
#ifdef MPI
    }
#endif
    kCalculateMolecularVirial(gpu);
  }
#ifdef MPI
  gpu_allreduce(gpu->pbForceAccumulator, gpu->sim.stride3);
#endif
  gpu->pbEnergyBuffer->Download();
  PMEDouble energy[EXTENDED_ENERGY_TERMS];
  for (int i = 0; i < EXTENDED_ENERGY_TERMS; i++) {
    energy[i] = 0.0;
  }
  pEnergy->total = 0.0;
  for (int i = 0; i < VIRIAL_OFFSET; i++) {
    unsigned long long int val = gpu->pbEnergyBuffer->_pSysData[i];
    if (val >= 0x8000000000000000ull) {
      energy[i] = -(PMEDouble)(val ^ 0xffffffffffffffffull) / ENERGYSCALE;
    }
    else {
      energy[i] = (PMEDouble)val / ENERGYSCALE;
    }
    pEnergy->total += energy[i];
  }
  for (int i = VIRIAL_OFFSET; i < ENERGY_TERMS; i++) {
    unsigned long long int val = gpu->pbEnergyBuffer->_pSysData[i];
    if (val >= 0x8000000000000000ull) {
      energy[i] = -(PMEDouble)(val ^ 0xffffffffffffffffull) / FORCESCALE;
    }
    else {
      energy[i] = (PMEDouble)val / FORCESCALE;
    }
  }

#ifdef use_DPFP
  if (gpu->imin != 0)
  {
    for (int i = ENERGY_TERMS; i < EXTENDED_ENERGY_TERMS; i++)
    {
      unsigned long long int val = gpu->pbEnergyBuffer->_pSysData[i];
      energy[i] = (double)(*((long long int*)(&val)));
      pEnergy->total         += energy[i];
    }
    energy[1] += energy[EVDWE_OFFSET];
    energy[10] += energy[EEDE_OFFSET];
  }
#endif

  // Force and energy computations performed on the CPU retain a record in host RAM.
  // That is now drawn upon to add their contributions to the appropriate components
  // of the total energy.

#ifdef MPI
  if (gpu->gpuID == 0) {
#endif
  pEnergy->total                 += gpu->vdw_recip + gpu->self_energy;
#ifdef MPI
  }
//NEB
  if ((*ineb > 0) && (*nstep%*nebfreq == 0)) {
    MPI_Allgather(&pEnergy->total, 1, MPI_PMEDOUBLE, gpu->pbNEBEnergyAll->_pSysData,
                                   1, MPI_PMEDOUBLE, MPI_COMM_WORLD);
    gpu->pbNEBEnergyAll->Upload();
  }

  if (gpu->gpuID == 0) {
#endif
    pEnergy->vdw_tot = energy[1] + gpu->vdw_recip;
#ifdef MPI
  }
  else {
    pEnergy->vdw_tot = energy[1];
  }
#endif
  pEnergy->vdw_dir = energy[1];
  pEnergy->vdw_recip = gpu->vdw_recip;
#ifdef MPI
  if (gpu->gpuID == 0) {
#endif
  pEnergy->elec_tot = energy[9] + energy[10] + gpu->self_energy;
#ifdef MPI
  }
  else {
    pEnergy->elec_tot = energy[9] + energy[10];
  }
#endif
  pEnergy->elec_dir = energy[10];
  pEnergy->elec_recip = energy[9];
  pEnergy->elec_nb_adjust = 0.0;
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

  // Grab virial if needed
  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
    virial[0] = 0.5 * energy[VIRIAL_OFFSET + 0];
    virial[1] = 0.5 * energy[VIRIAL_OFFSET + 1];
    virial[2] = 0.5 * energy[VIRIAL_OFFSET + 2];
#ifdef MPI
    if (gpu->gpuID == 0) {
#endif
    virial[0] -= 0.5 * (gpu->ee_plasma + 2.0 * gpu->vdw_recip);
    virial[1] -= 0.5 * (gpu->ee_plasma + 2.0 * gpu->vdw_recip);
    virial[2] -= 0.5 * (gpu->ee_plasma + 2.0 * gpu->vdw_recip);
#ifdef MPI
    }
#endif
    ekcmt[0] = energy[VIRIAL_OFFSET + 3];
    ekcmt[1] = energy[VIRIAL_OFFSET + 4];
    ekcmt[2] = energy[VIRIAL_OFFSET + 5];
  }
}

//---------------------------------------------------------------------------------------------
// gpu_pme_force_: compute forces for (explicit solvent) PME simulations
//
// Arguments:
//   ewaldcof:    the Ewald coefficient, 1 / (2 * Gaussian sigma)
//   vol:         the simulation box volue
//   virial:      trace of the virial tensor
//   ekcmt:       kinetic energy of the center of mass of the sub-molecules (three elements
//                for three velocity components--it'll get summed into one number later)
//   nstep:       the step number
//   dt:          the time step, dt from the input file in units of ps
//   atm_crd:     passed in to support a transfer of select coordinates to the host for
//                special force computations
//   atm_frc:     passed in to allow the results of host-side force computations to be
//                uploaded to the device
//---------------------------------------------------------------------------------------------
extern "C" void gpu_pme_force_(double* ewaldcof, double* vol, double virial[3],
                               double ekcmt[3], int *nstep, double *dt, double atm_crd[][3],
                               double atm_frc[][3])
{
  PRINTMETHOD("gpu_pme_force");

  // Rebuild neighbor list
  gpu_build_neighbor_list_();

  // Clear forces
  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
    kClearForces(gpu, gpu->sim.NLNonbondEnergyWarps);
  }
  else {
    kClearForces(gpu, gpu->sim.NLNonbondForcesWarps);
  }

  if (gpu->ntf != 8) {

    // Download critical coordinates for force calculations to occur on the host side
    //if (gpu->sim.nShuttle > 0) {
    //  kRetreiveSimData(gpu, atm_crd);
    //}

    // Local forces
#ifdef MPI
    if (gpu->bCalculateLocalForces) {
#endif
    kExecuteBondWorkUnits(gpu, 0);
    kCalculateNMRForces(gpu);
    if (gpu->iphmd == 3) {
      kCalculateSelfPHMD(*vol, *ewaldcof, gpu);
    }
#ifdef MPI
    }
    else {
      kPMEClearChargeGridBuffer(gpu);
    }
#endif
    // Non-bonded energy: every process will do either reciprocal space or direct sum
    // interactions and therefore must re-image coordinates.
    kPMEGetGridWeights(gpu);
#ifdef MPI
    if (gpu->bCalculateReciprocalSum) {
#endif
    // Reciprocal space forces
    kPMEFillChargeGridBuffer(gpu);
    kPMEReduceChargeGridBuffer(gpu);
    kPMEForwardFFT(gpu);
    kPMEScalarSumRC(gpu, *vol);
    kPMEBackwardFFT(gpu);
    kPMEGradSum(gpu);
#ifdef MPI
    }
    if (gpu->bCalculateDirectSum) {
#endif

    // Direct space forces
    kCalculatePMENonbondForces(gpu);

    // Run CPU (host)-side potential and force calculations while the non-bonded
    // force kernel is running.  This is the best opportunity to mask the cost.

    // Electric Field Forces
#ifdef MPI
    }
    if (gpu->bCalculateLocalForces) {
#endif
    if (gpu->sim.efx != 0 || gpu->sim.efy != 0 || gpu->sim.efz != 0) {
      kCalculateEFieldForces(gpu, *nstep, *dt);
    }
#ifdef MPI
    }
#endif

    // Commute forces on extra points
    if (gpu->sim.EPs > 0) {
      kOrientForces(gpu);
    }

    // Download forces on critical atoms for additional CPU calculations
    //if (gpu->sim.nShuttle > 0) {
    //  kRetrieveSimData(gpu, atm_frc, 2);
    //}

  }
  // End contingency for ntf != 8 (if ntf == 8, all bonds, angles,
  // dihedrals, and non-bonded interactions would be omitted)

  // Grab virial and ekcmt if needed
  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
#ifdef MPI
    if (gpu->gpuID == 0) {
#endif
      kCalculateCOMKineticEnergy(gpu);
      kReduceCOMKineticEnergy(gpu);
#ifdef MPI
    }
#endif
    kCalculateMolecularVirial(gpu);
    gpu->pbEnergyBuffer->Download();
    for (int i = 0; i < 6; i++) {
      unsigned long long int val = gpu->pbEnergyBuffer->_pSysData[VIRIAL_OFFSET + i];
      PMEDouble dval;
      if (val >= 0x8000000000000000ull) {
        dval = -(PMEDouble)(val ^ 0xffffffffffffffffull) / FORCESCALE;
      }
      else {
        dval = (PMEDouble)val / FORCESCALE;
      }
      if (i < 3) {
        virial[i] = 0.5 * dval;
      }
      else {
        ekcmt[i - 3] = dval;
      }
    }
#ifdef MPI
    if (gpu->gpuID == 0) {
#endif
    virial[0] -= 0.5 * (gpu->ee_plasma + (2.0 * gpu->vdw_recip));
    virial[1] -= 0.5 * (gpu->ee_plasma + (2.0 * gpu->vdw_recip));
    virial[2] -= 0.5 * (gpu->ee_plasma + (2.0 * gpu->vdw_recip));
#ifdef MPI
    }
#endif
  }
  // End contingency for constant pressure simulations with Berendsen barostat

#ifdef MPI
  gpu_allreduce(gpu->pbForceAccumulator, gpu->sim.stride3);
#endif
}

//---------------------------------------------------------------------------------------------
// gpu_ti_pme_ene_: get the energy for a Thermodynamic Integration simulation involving
//                  Particle Mesh Ewald (implies peridoic boundary conditions and explicit
//                  solvent).
//
// Arguments:
//   ewaldcof:    the Ewald coefficient, 1 / (2 * Gaussian sigma)
//   vol:         the simulation box volue
//   pEnergy:     struct to hold various components of the energy function
//   pSCEnergy:   corresponding energy struct for the soft-core parts of the system
//   enmr:        energy of NMR restraints
//   virial:      trace of the virial tensor
//   ekcmt:       kinetic energy of the center of mass of the sub-molecules (three elements
//                for three velocity components--it'll get summed into one number later)
//   nstep:       the step number
//   dt:          the time step, dt from the input file in units of ps
//   atm_crd:     passed in to support a transfer of select coordinates to the host for
//                special force computations
//   atm_frc:     passed in to allow the results of host-side force computations to be
//                uploaded to the device
//---------------------------------------------------------------------------------------------
extern "C" void gpu_ti_pme_ene_(double* ewaldcof, double* vol, pme_pot_ene_rec* pEnergy,
                                afe_gpu_sc_ene_rec* pSCEnergy, double enmr[3],
                                double virial[3], double ekcmt[3], int* nstep, double* dt,
                                double atm_crd[][3], double atm_frc[][3], double bartot[])
{
  PRINTMETHOD("gpu_ti_pme_ene");

  // Rebuild neighbor list
  gpu_build_neighbor_list_();

  // Clear forces
  kClearForces(gpu, gpu->sim.NLNonbondEnergyWarps);

  // Contingency for any case where forces are calculated
  if (gpu->ntf != 8) {

    // Local energy
#ifdef MPI
    if (gpu->bCalculateLocalForces) {
#endif
    kExecuteBondWorkUnits(gpu, 1);
    kCalculateNMREnergy(gpu);
#ifdef MPI
    }
#endif

    // Reciprocal space energy
#ifdef MPI
    if (gpu->bCalculateReciprocalSum) {
#endif
    gpu->sim.AFE_recip_region = 1;
    kPMEGetGridWeights(gpu);
    kPMEClearChargeGridBuffer(gpu);
    kPMEFillChargeGridBuffer(gpu);
    kPMEReduceChargeGridBuffer(gpu);
    kPMEForwardFFT(gpu);
    kPMEScalarSumRCEnergy(gpu, *vol);
    kPMEBackwardFFT(gpu);
    kPMEGradSum(gpu);

    // Multiple calculations are only needed if we haven't removed charges
    if ((gpu->sim.ti_mode == 1) || (gpu->sim.ti_mode == 2)) {
      gpu->sim.AFE_recip_region = 2;
      kPMEGetGridWeights(gpu);
      kPMEClearChargeGridBuffer(gpu);
      kPMEFillChargeGridBuffer(gpu);
      kPMEReduceChargeGridBuffer(gpu);
      kPMEForwardFFT(gpu);
      kPMEScalarSumRCEnergy(gpu, *vol);
      kPMEBackwardFFT(gpu);
      kPMEGradSum(gpu);
    }
#ifdef MPI
  }
  else {
    kPMEGetGridWeights(gpu);
  }
#endif

  // Direct space energy
#ifdef MPI
  if (gpu->bCalculateDirectSum) {
#endif
  kCalculatePMENonbondEnergy(gpu);
#ifdef MPI
  }
#endif

  // Electric Field Energy
#ifdef MPI
  if (gpu->bCalculateLocalForces) {
#endif
  if (gpu->sim.efx != 0 || gpu->sim.efy != 0 || gpu->sim.efz != 0) {
    kCalculateEFieldEnergy(gpu, *nstep, *dt);
  }
#ifdef MPI
  }
#endif
  if (gpu->sim.EPs > 0)
    kOrientForces(gpu);
  }
  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
#ifdef MPI
    if (gpu->gpuID == 0) {
#endif
    kCalculateCOMKineticEnergy(gpu);
    kReduceCOMKineticEnergy(gpu);
#ifdef MPI
    }
#endif
    kCalculateMolecularVirial(gpu);
  }
  if (gpu->sim.TIlinearAtmCnt > 0) {
    kAFEExchangeFrc(gpu);
  }

#ifdef MPI
  gpu_allreduce(gpu->pbForceAccumulator, gpu->sim.stride3);
#endif
  gpu->pbEnergyBuffer->Download();
  PMEDouble energy[EXTENDED_ENERGY_TERMS];
  for (int i = 0; i < EXTENDED_ENERGY_TERMS; i++) {
    energy[i] = 0.0;
  }
  pEnergy->total = 0.0;
  for (int i = 0; i < VIRIAL_OFFSET; i++) {
    unsigned long long int val  = gpu->pbEnergyBuffer->_pSysData[i];
    if (val >= 0x8000000000000000ull) {
      energy[i] = -(PMEDouble)(val ^ 0xffffffffffffffffull) / ENERGYSCALE;
    }
    else {
      energy[i]               = (PMEDouble)val / ENERGYSCALE;
    }
    pEnergy->total         += energy[i];
  }


  for (int i = VIRIAL_OFFSET; i < ENERGY_TERMS; i++) {
    unsigned long long int val  = gpu->pbEnergyBuffer->_pSysData[i];
    if (val >= 0x8000000000000000ull) {
      energy[i] = -(PMEDouble)(val ^ 0xffffffffffffffffull) / FORCESCALE;
    }
    else {
      energy[i] = (PMEDouble)val / FORCESCALE;
    }
  }

#ifdef use_DPFP
  if (gpu->imin != 0)
  {
    for (int i = ENERGY_TERMS; i < EXTENDED_ENERGY_TERMS; i++)
    {
      unsigned long long int val = gpu->pbEnergyBuffer->_pSysData[i];
      energy[i] = (double)(*((long long int*)(&val)));
      pEnergy->total         += energy[i];
    }
    energy[1] += energy[EVDWE_OFFSET];
    energy[10] += energy[EEDE_OFFSET];
  }
#endif

  // Need to set DVDL to AFEBuffer[0]
  gpu->pbAFEBuffer->Download();
  PMEDouble AFEterm[AFE_TERMS];
  for (int i = 0; i < AFE_TERMS; i++) {
    AFEterm[i] = 0.0;
  }
  for (int i = 0; i < AFE_TERMS; i++) {
    unsigned long long int val  = gpu->pbAFEBuffer->_pSysData[i];
    if (val >= 0x8000000000000000ull) {
      AFEterm[i] = -(PMEDouble)(val ^ 0xffffffffffffffffull) / ENERGYSCALE;
    }
    else {
      AFEterm[i] = (PMEDouble)val / ENERGYSCALE;
    }
  }

#ifdef MPI
  if (gpu->gpuID == 0) {
#endif
  pEnergy->total += gpu->vdw_recip + gpu->self_energy;
#ifdef MPI
  }
  if (gpu->gpuID == 0) {
#endif
  pEnergy->vdw_tot = energy[1] + gpu->vdw_recip;
#ifdef MPI
  }
  else {
    pEnergy->vdw_tot = energy[1];
  }
#endif
  pEnergy->vdw_dir   = energy[1];
  pEnergy->vdw_recip = gpu->vdw_recip;
#ifdef MPI
  if (gpu->gpuID == 0) {
#endif
    pEnergy->elec_tot           = energy[9] + energy[10] + gpu->self_energy;
#ifdef MPI
  }
  else {
    pEnergy->elec_tot           = energy[9] + energy[10];
  }
#endif
  pEnergy->elec_dir               = energy[10];
  pEnergy->elec_recip             = energy[9];
  pEnergy->elec_nb_adjust         = 0.0;
  pEnergy->hbond                  = 0.0;
  pEnergy->bond                   = energy[3];
  pEnergy->angle                  = energy[4];
  pEnergy->dihedral               = energy[5];
  pEnergy->vdw_14                 = energy[7];
  pEnergy->elec_14                = energy[6];
  pEnergy->restraint              = energy[8] + energy[14] + energy[15] + energy[16];
  pEnergy->angle_ub               = energy[11];
  pEnergy->imp                    = energy[12];
  pEnergy->cmap                   = energy[13];
  enmr[0]                         = energy[14];
  enmr[1]                         = energy[15];
  enmr[2]                         = energy[16];
  pEnergy->efield                 = energy[17];
  pSCEnergy->dvdl                = AFEterm[0];
  pSCEnergy->bond_R1             = AFEterm[1];
  pSCEnergy->bond_R2             = AFEterm[2];
  pSCEnergy->angle_R1            = AFEterm[3];
  pSCEnergy->angle_R2            = AFEterm[4];
  pSCEnergy->dihedral_R1         = AFEterm[5];
  pSCEnergy->dihedral_R2         = AFEterm[6];
  pSCEnergy->sc_res_dist_R1      = AFEterm[7];
  pSCEnergy->sc_res_dist_R2      = AFEterm[8];
  pSCEnergy->sc_res_ang_R1       = AFEterm[9];
  pSCEnergy->sc_res_ang_R2       = AFEterm[10];
  pSCEnergy->sc_res_tors_R1      = AFEterm[11];
  pSCEnergy->sc_res_tors_R2      = AFEterm[12];
  pSCEnergy->vdw_dir_R1          = AFEterm[13];
  pSCEnergy->vdw_dir_R2          = AFEterm[14];
  pSCEnergy->elec_dir_R1         = AFEterm[15];
  pSCEnergy->elec_dir_R2         = AFEterm[16];
  pSCEnergy->vdw_14_R1           = AFEterm[17];
  pSCEnergy->vdw_14_R2           = AFEterm[18];
  pSCEnergy->elec_14_R1          = AFEterm[19];
  pSCEnergy->elec_14_R2          = AFEterm[20];
  pSCEnergy->vdw_der_R1          = AFEterm[21];
  pSCEnergy->vdw_der_R2          = AFEterm[22];
  pSCEnergy->elec_der_R1         = AFEterm[23];
  pSCEnergy->elec_der_R2         = AFEterm[24];

  if (gpu->sim.ifmbar > 0) {
    gpu->pbBarTot->Download(); //dies here. value is probably garbage on gpu
    PMEDouble mbarene;
    for (int i = 0; i < gpu->sim.bar_states; i++) {
      unsigned long long int val  = gpu->pbBarTot->_pSysData[i];
      if (val >= 0x8000000000000000ull) {
        mbarene = -(PMEDouble)(val ^ 0xffffffffffffffffull) / ENERGYSCALE;
      }
      else {
        mbarene = (PMEDouble)val / ENERGYSCALE;
      }
      bartot[i] += mbarene;
    }
  }

  // Grab virial if needed
  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
    virial[0] = 0.5 * energy[VIRIAL_OFFSET + 0];
    virial[1] = 0.5 * energy[VIRIAL_OFFSET + 1];
    virial[2] = 0.5 * energy[VIRIAL_OFFSET + 2];
#ifdef MPI
    if (gpu->gpuID == 0) {
#endif
      virial[0] -= 0.5*(gpu->ee_plasma + 2.0*gpu->vdw_recip);
      virial[1] -= 0.5*(gpu->ee_plasma + 2.0*gpu->vdw_recip);
      virial[2] -= 0.5*(gpu->ee_plasma + 2.0*gpu->vdw_recip);
#ifdef MPI
    }
#endif
    ekcmt[0] = energy[VIRIAL_OFFSET + 3];
    ekcmt[1] = energy[VIRIAL_OFFSET + 4];
    ekcmt[2] = energy[VIRIAL_OFFSET + 5];
  }
}

//---------------------------------------------------------------------------------------------
// gpu_ti_pme_force_: get forces for a TI simulation involving PME
//
// Arguments:
//   ewaldcof:    the Ewald coefficient, 1 / (2 * Gaussian sigma)
//   vol:         the simulation box volue
//   pEnergy:     struct to hold various components of the energy function
//   pSCEnergy:   corresponding energy struct for the soft-core parts of the system
//   enmr:        energy of NMR restraints
//   virial:      trace of the virial tensor
//   ekcmt:       kinetic energy of the center of mass of the sub-molecules (three elements
//                for three velocity components--it'll get summed into one number later)
//   nstep:       the step number
//   dt:          the time step, dt from the input file in units of ps
//   atm_crd:     passed in to support a transfer of select coordinates to the host for
//                special force computations
//   atm_frc:     passed in to allow the results of host-side force computations to be
//                uploaded to the device
//---------------------------------------------------------------------------------------------
extern "C" void gpu_ti_pme_force_(double* ewaldcof, double* vol, double virial[3],
                                  double ekcmt[3], int *nstep, double *dt, double atm_crd[][3],
                                  double atm_frc[][3])
{
  PRINTMETHOD("gpu_ti_pme_force");

  // Rebuild neighbor list
  gpu_build_neighbor_list_();

  // Clear forces
  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
    kClearForces(gpu, gpu->sim.NLNonbondEnergyWarps);
  }
  else {
    kClearForces(gpu, gpu->sim.NLNonbondForcesWarps);
  }
  if (gpu->ntf != 8) {

    // Local forces
#ifdef MPI
    if (gpu->bCalculateLocalForces) {
#endif
    kExecuteBondWorkUnits(gpu, 0);
    kCalculateNMRForces(gpu);
#ifdef MPI
    }
#endif

#ifdef MPI
    if (gpu->bCalculateReciprocalSum) {
#endif
    // Reciprocal forces
    gpu->sim.AFE_recip_region  = 1;
    kPMEGetGridWeights(gpu);
    kPMEClearChargeGridBuffer(gpu);
    kPMEFillChargeGridBuffer(gpu);
    kPMEReduceChargeGridBuffer(gpu);
    kPMEForwardFFT(gpu);
    kPMEScalarSumRC(gpu, *vol);
    kPMEBackwardFFT(gpu);
    kPMEGradSum(gpu);

    // The reciprocal portion needs to be done twice, but multiple
    // calculations are only needed if we haven't removed charges
    if ((gpu->sim.ti_mode == 1) || (gpu->sim.ti_mode == 2)) {
      gpu->sim.AFE_recip_region  = 2;
      kPMEGetGridWeights(gpu);
      kPMEClearChargeGridBuffer(gpu);
      kPMEFillChargeGridBuffer(gpu);
      kPMEReduceChargeGridBuffer(gpu);
      kPMEForwardFFT(gpu);
      kPMEScalarSumRC(gpu, *vol);
      kPMEBackwardFFT(gpu);
      kPMEGradSum(gpu);
    }
#ifdef MPI
    }
#endif

    // Direct energy
#ifdef MPI
    if (gpu->bCalculateDirectSum) {
#endif
    kCalculatePMENonbondForces(gpu);
#ifdef MPI
    }
#endif

    // Electric Field Forces
#ifdef MPI
    if (gpu->bCalculateLocalForces) {
#endif
    if (gpu->sim.efx != 0 || gpu->sim.efy != 0 || gpu->sim.efz != 0) {
      kCalculateEFieldForces(gpu, *nstep, *dt);
    }
#ifdef MPI
    }
#endif
    if (gpu->sim.EPs > 0) {
      kOrientForces(gpu);
    }
  }

  // Grab virial and ekcmt if needed
  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
#ifdef MPI
    if (gpu->gpuID == 0) {
#endif
      kCalculateCOMKineticEnergy(gpu);
      kReduceCOMKineticEnergy(gpu);
#ifdef MPI
    }
#endif
    kCalculateMolecularVirial(gpu);
    gpu->pbEnergyBuffer->Download();
    for (int i = 0; i < 6; i++) {
      unsigned long long int val = gpu->pbEnergyBuffer->_pSysData[VIRIAL_OFFSET + i];
      PMEDouble dval;
      if (val >= 0x8000000000000000ull) {
        dval = -(PMEDouble)(val ^ 0xffffffffffffffffull) / FORCESCALE;
      }
      else {
        dval = (PMEDouble)val / FORCESCALE;
      }
      if (i < 3) {
        virial[i] = 0.5 * dval;
      }
      else {
        ekcmt[i-3] = dval;
      }
    }
#ifdef MPI
    if (gpu->gpuID == 0) {
#endif
    virial[0]                  -= 0.5 * (gpu->ee_plasma + 2.0 * gpu->vdw_recip);
    virial[1]                  -= 0.5 * (gpu->ee_plasma + 2.0 * gpu->vdw_recip);
    virial[2]                  -= 0.5 * (gpu->ee_plasma + 2.0 * gpu->vdw_recip);
#ifdef MPI
    }
#endif
  }
  if (gpu->sim.TIlinearAtmCnt > 0) {

    // For now it doesn't matter if this is before or after the virial stuff,
    // because we don't allow barostat == 1.  But later, go back and figure
    // out if this should be before or after the virial stuff.
    kAFEExchangeFrc(gpu);
  }
#ifdef MPI
  gpu_allreduce(gpu->pbForceAccumulator, gpu->sim.stride3);
#endif
}

//---------------------------------------------------------------------------------------------
// gpu_ips_ene_: get the energy for the system in simulations with Isotropic Periodic Sums
//
// Arguments:
//   vol:       the simulation box volume
//   pEnergy:   potential energy accumulators for PME simulations (and the like... this is one
//              of the like)
//   enmr:      NMR restraint energy (see description in nmr_calls.F90 in the Fortran source)
//   virial:    trace of the virial tensor
//   ekcmt:     kinetic energy of the particles relative to the center of mass
//---------------------------------------------------------------------------------------------
extern "C" void gpu_ips_ene_(double* vol, pme_pot_ene_rec* pEnergy, double enmr[3],
                             double virial[3], double ekcmt[3])
{
  PRINTMETHOD("gpu_ips_ene");

  // Rebuild neighbor list
  gpu_build_neighbor_list_();

  // Clear forces, then do computations unless ntf == 8 (meaning that all
  // bonds, angles, dihedrals, and non-bonded interactions are omitted)
  kClearForces(gpu, gpu->sim.NLNonbondEnergyWarps);
  if (gpu->ntf != 8) {

    // Local energy
#ifdef MPI
    if (gpu->bCalculateLocalForces) {
#endif
    kExecuteBondWorkUnits(gpu, 1);
    kCalculateNMREnergy(gpu);
#ifdef MPI
    }
#endif
    kNLCalculateCellCoordinates(gpu);
    kCalculateIPSNonbondEnergy(gpu);
    if (gpu->sim.EPs > 0) {
      kOrientForces(gpu);
    }
  }
  // End contingency for ntf == 8

  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
#ifdef MPI
    if (gpu->gpuID == 0) {
#endif
      kCalculateCOMKineticEnergy(gpu);
      kReduceCOMKineticEnergy(gpu);
#ifdef MPI
    }
#endif
    kCalculateMolecularVirial(gpu);
  }

  gpu->pbEnergyBuffer->Download();
  PMEDouble energy[EXTENDED_ENERGY_TERMS];
  pEnergy->total = 0.0;
  for (int i = 0; i < VIRIAL_OFFSET; i++) {
    unsigned long long int val = gpu->pbEnergyBuffer->_pSysData[i];
    if (val >= 0x8000000000000000ull) {
      energy[i] = -(PMEDouble)(val ^ 0xffffffffffffffffull) / ENERGYSCALE;
    }
    else {
      energy[i] = (PMEDouble)val / ENERGYSCALE;
    }
    pEnergy->total += energy[i];
  }
  for (int i = VIRIAL_OFFSET; i < ENERGY_TERMS; i++) {
    unsigned long long int val = gpu->pbEnergyBuffer->_pSysData[i];
    if (val >= 0x8000000000000000ull) {
      energy[i] = -(PMEDouble)(val ^ 0xffffffffffffffffull) / FORCESCALE;
    }
    else {
      energy[i] = (PMEDouble)val / FORCESCALE;
    }
  }
  for (int i = ENERGY_TERMS; i < EXTENDED_ENERGY_TERMS; i++) {
    unsigned long long int val = gpu->pbEnergyBuffer->_pSysData[i];
    energy[i] = (double)(*((long long int*)(&val)));
  }

#ifdef MPI
  if (gpu->gpuID == 0) {
#endif
  pEnergy->total += gpu->sim.EIPSEL + gpu->sim.EIPSNB + gpu->sim.eipssnb + gpu->sim.eipssel;
#ifdef MPI
  }
  if (gpu->gpuID > 0) {
    pEnergy->vdw_tot = energy[1] + energy[EVDWE_OFFSET];
  }
  else {
#endif
    pEnergy->vdw_tot = energy[1] + energy[EVDWE_OFFSET] + gpu->sim.EIPSNB + gpu->sim.eipssnb;
#ifdef MPI
  }
#endif
  pEnergy->vdw_dir = pEnergy->vdw_tot;
#ifdef MPI
  if (gpu->gpuID > 0) {
    pEnergy->elec_tot = energy[10];
  }
  else {
#endif
  pEnergy->elec_tot = energy[10] + gpu->sim.EIPSEL + gpu->sim.eipssel;
#ifdef MPI
  }
#endif
  pEnergy->elec_dir = pEnergy->elec_tot;
  pEnergy->elec_recip = 0.0;
  pEnergy->elec_nb_adjust = 0.0;
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

  // Grab virial if needed
  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
    virial[0] = 0.5 * energy[VIRIAL_OFFSET + 0];
    virial[1] = 0.5 * energy[VIRIAL_OFFSET + 1];
    virial[2] = 0.5 * energy[VIRIAL_OFFSET + 2];
    ekcmt[0] = energy[VIRIAL_OFFSET + 3];
    ekcmt[1] = energy[VIRIAL_OFFSET + 4];
    ekcmt[2] = energy[VIRIAL_OFFSET + 5];
#ifdef MPI
    if (gpu->gpuID == 0) {
#endif
    virial[0] += gpu->sim.virips;
    virial[1] += gpu->sim.virips;
    virial[2] += gpu->sim.virips;
#ifdef MPI
    }
#endif
  }
  // End contingency for constant pressure simulations with Berendsen barostat

#ifdef MPI
  gpu_allreduce(gpu->pbForceAccumulator, gpu->sim.stride3);
#endif
}

//---------------------------------------------------------------------------------------------
// gpu_ips_force_: get the forces for simulations comptued with Isotropic Periodic Sums
//
// Arguments:
//   vol:      simulation box volume
//   virial:   trace of the virial tensor
//   ekcmt:    kinetic energy of the particles relative to the center of mass
//---------------------------------------------------------------------------------------------
extern "C" void gpu_ips_force_(double* vol, double virial[3], double ekcmt[3])
{
  PRINTMETHOD("gpu_ips_force");

  // Rebuild neighbor list
  gpu_build_neighbor_list_();

  // Clear forces
  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
    kClearForces(gpu, gpu->sim.NLNonbondEnergyWarps);
  }
  else {
    kClearForces(gpu, gpu->sim.NLNonbondForcesWarps);
  }

  // Contingency for calculating any bond, angle, dihedral, or non-bonded parameters
  if (gpu->ntf != 8) {

    // Local forces
#ifdef MPI
    if (gpu->bCalculateLocalForces) {
#endif
    kExecuteBondWorkUnits(gpu, 0);
    kCalculateNMRForces(gpu);
#ifdef MPI
    }
#endif

    // Direct energy and force reduction
#ifdef MPI
    if (gpu->bCalculateDirectSum) {
#endif
    kNLCalculateCellCoordinates(gpu);
    kCalculateIPSNonbondForces(gpu);
#ifdef MPI
    }
#endif
    if (gpu->sim.EPs > 0) {
      kOrientForces(gpu);
    }
  }

  // Grab virial and ekcmt if needed
  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
#ifdef MPI
    if (gpu->gpuID == 0) {
#endif
    kCalculateCOMKineticEnergy(gpu);
    kReduceCOMKineticEnergy(gpu);
#ifdef MPI
    }
#endif
    kCalculateMolecularVirial(gpu);
    gpu->pbEnergyBuffer->Download();
    for (int i = 0; i < 6; i++) {
      unsigned long long int val = gpu->pbEnergyBuffer->_pSysData[VIRIAL_OFFSET + i];
      PMEDouble dval;
      if (val >= 0x8000000000000000ull) {
        dval = -(PMEDouble)(val ^ 0xffffffffffffffffull) / FORCESCALE;
      }
      else {
        dval = (PMEDouble)val / FORCESCALE;
      }
      if (i < 3) {
        virial[i] = 0.5 * dval;
      }
      else {
        ekcmt[i - 3] = dval;
      }
    }
#ifdef MPI
    if (gpu->gpuID == 0) {
#endif
    virial[0] += gpu->sim.virips;
    virial[1] += gpu->sim.virips;
    virial[2] += gpu->sim.virips;
#ifdef MPI
    }
#endif
  }
  // End contingency for constant pressure simulations with Berendsen barostat

#ifdef MPI
  gpu_allreduce(gpu->pbForceAccumulator, gpu->sim.stride3);
#endif
}

//---------------------------------------------------------------------------------------------
// gpu_update_box_: update box information
//
// Arguments:
//   ucell:      transformation matrix from fractional coordinates to Cartesian space
//   recip:      inverse of ucell
//   uc_volume:  unit cell volume
//---------------------------------------------------------------------------------------------
extern "C" void gpu_update_box_(double ucell[3][3], double recip[3][3], double* uc_volume)
{
  PRINTMETHOD("gpu_update_box");

  // Recalculate nonbond skin
  double invU[9];
  invU[0] = ucell[0][0];
  invU[1] = ucell[1][0];
  invU[2] = ucell[2][0];
  invU[3] = ucell[0][1];
  invU[4] = ucell[1][1];
  invU[5] = ucell[2][1];
  invU[6] = ucell[0][2];
  invU[7] = ucell[1][2];
  invU[8] = ucell[2][2];
  double cdepth[3];
  HessianNorms(invU, cdepth);
  double skin = cdepth[0]/gpu->sim.xcells - gpu->sim.cut;
  double yskin = cdepth[1]/gpu->sim.ycells - gpu->sim.cut;
  double zskin = cdepth[2]/gpu->sim.zcells - gpu->sim.cut;

  if (yskin < skin) {
    skin = yskin;
  }
  if (zskin < skin) {
    skin = zskin;
  }

  skin *= 0.99;

  gpu->sim.xcell = (PMEFloat)(gpu->sim.a / (PMEDouble)gpu->sim.xcells);
  gpu->sim.ycell = (PMEFloat)(gpu->sim.b / (PMEDouble)gpu->sim.ycells);
  gpu->sim.zcell = (PMEFloat)(gpu->sim.c / (PMEDouble)gpu->sim.zcells);
  gpu->sim.pi_vol_inv = (PMEFloat)(1.0 / (PI * *uc_volume));
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      gpu->sim.ucellf[i][j] = ucell[j][i];
      gpu->sim.recipf[i][j] = recip[j][i];
      gpu->sim.ucell[i][j] = gpu->sim.ucellf[i][j];
      gpu->sim.recip[i][j] = gpu->sim.recipf[i][j];
    }
  }
  if (gpu->sim.is_orthog) {
    for (int i = 0; i < NEIGHBOR_CELLS; i++) {
      gpu->sim.cellOffset[i][0] = (PMEDouble)cellOffset[i][0] * gpu->sim.xcell;
      gpu->sim.cellOffset[i][1] = (PMEDouble)cellOffset[i][1] * gpu->sim.ycell;
      gpu->sim.cellOffset[i][2] = (PMEDouble)cellOffset[i][2] * gpu->sim.zcell;
    }
  }
  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gpu_pressure_scale_: carry out coordinate scaling to fulfill constant pressure conditions
//
// Arguments:
//   ucell:      transformation matrix from fractional coordinates to Cartesian space
//   recip:      inverse of ucell
//   uc_volume:  unit cell volume
//---------------------------------------------------------------------------------------------
extern "C" void gpu_pressure_scale_(double ucell[3][3], double recip[3][3], double* uc_volume)
{
  PRINTMETHOD("gpu_pressure_scale");
  NTPData* pNTPData = gpu->pbNTPData->_pSysData;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pNTPData->last_recip[j * 3 + i] = pNTPData->recip[j * 3 + i];
      pNTPData->ucell[j*3 + i] = ucell[i][j];
      pNTPData->ucellf[j*3 + i] = ucell[i][j];
      pNTPData->recip[j*3 + i] = recip[i][j];
      pNTPData->recipf[j*3 + i] = recip[i][j];
    }
  }

  // Recalculate a, b, c
  double a = pNTPData->ucell[0];
  double b = sqrt(pNTPData->ucell[1] * pNTPData->ucell[1] +
                  pNTPData->ucell[4] * pNTPData->ucell[4]);
  double c = sqrt(pNTPData->ucell[2] * pNTPData->ucell[2] +
                  pNTPData->ucell[5] * pNTPData->ucell[5] +
                  pNTPData->ucell[8] * pNTPData->ucell[8]);

  // Recalculate nonbond skin
  double invU[9];
  invU[0] = ucell[0][0];
  invU[1] = ucell[1][0];
  invU[2] = ucell[2][0];
  invU[3] = ucell[0][1];
  invU[4] = ucell[1][1];
  invU[5] = ucell[2][1];
  invU[6] = ucell[0][2];
  invU[7] = ucell[1][2];
  invU[8] = ucell[2][2];
  double cdepth[3];
  HessianNorms(invU, cdepth);
  double skin = cdepth[0]/gpu->sim.xcells - gpu->sim.cut;
  double yskin = cdepth[1]/gpu->sim.ycells - gpu->sim.cut;
  double zskin = cdepth[2]/gpu->sim.zcells - gpu->sim.cut;

  // The nonbond_skin attribute of the cudaSimulation struct is not updated because
  // it is not actually referenced during the pairlist margin check: the quantities
  // derived from it in gpu->sim.NTPData are.  The entire cudaSimulation is only
  // ported back up to the GPU in the MC barostat case.
  if (yskin < skin) {
    skin = yskin;
  }
  if (zskin < skin) {
    skin = zskin;
  }

  // Check if skin has fallen below acceptable threshold
  if (skin <= 0.5) {
    printf("ERROR: Calculation halted.");
    printf("  Periodic box dimensions have changed too much from their initial values.\n");
    printf("  Your system density has likely changed by a large amount, probably from\n");
    printf("  starting the simulation from a structure a long way from equilibrium.\n");
    printf("\n");
    printf("  [Although this error can also occur if the simulation has blown up for some "
           "reason]\n");
    printf("\n");
    printf("  The GPU code does not automatically reorganize grid cells and thus you\n");
    printf("  will need to restart the calculation from the previous restart file.\n");
    printf("  This will generate new grid cells and allow the calculation to continue.\n");
    printf("  It may be necessary to repeat this restarting multiple times if your system\n");
    printf("  is a long way from an equilibrated density.\n");
    printf("\n");
    printf("  Alternatively you can run with the CPU code until the density has converged\n");
    printf("  and then switch back to the GPU code.\n");
    printf("\n");
    exit(-1);
  }
  skin *= 0.99;
  pNTPData->one_half_nonbond_skin_squared  = (gpu->sim.skinPermit * skin);
  pNTPData->one_half_nonbond_skin_squared *= pNTPData->one_half_nonbond_skin_squared;
  pNTPData->cutPlusSkin2 = (gpu->sim.cut + skin) * (gpu->sim.cut + skin);
  gpu->pbNTPData->Upload();

  // Handle MC barostat
  if (gpu->sim.barostat != 1) {
    kClearSoluteCOM(gpu);
    gpu->sim.af    = (PMEFloat)a;
    gpu->sim.bf    = (PMEFloat)b;
    gpu->sim.cf    = (PMEFloat)c;
    gpu->sim.a     = (PMEDouble)gpu->sim.af;
    gpu->sim.b     = (PMEDouble)gpu->sim.bf;
    gpu->sim.c     = (PMEDouble)gpu->sim.cf;
    gpu->sim.xcell = (PMEFloat)(gpu->sim.a / (PMEDouble)gpu->sim.xcells);
    gpu->sim.ycell = (PMEFloat)(gpu->sim.b / (PMEDouble)gpu->sim.ycells);
    gpu->sim.zcell = (PMEFloat)(gpu->sim.c / (PMEDouble)gpu->sim.zcells);
    gpu->sim.pi_vol_inv = (PMEFloat)(1.0 / (PI * *uc_volume));
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        gpu->sim.ucellf[i][j] = ucell[j][i];
        gpu->sim.recipf[i][j] = recip[j][i];
        gpu->sim.ucell[i][j] = gpu->sim.ucellf[i][j];
        gpu->sim.recip[i][j] = gpu->sim.recipf[i][j];
      }
    }
    if (gpu->sim.is_orthog) {
      for (int i = 0; i < NEIGHBOR_CELLS; i++) {
        gpu->sim.cellOffset[i][0] = (PMEDouble)cellOffset[i][0] * gpu->sim.xcell;
        gpu->sim.cellOffset[i][1] = (PMEDouble)cellOffset[i][1] * gpu->sim.ycell;
        gpu->sim.cellOffset[i][2] = (PMEDouble)cellOffset[i][2] * gpu->sim.zcell;
      }
    }
    gpuCopyConstants();
  }
  kCalculateSoluteCOM(gpu);
  kPressureScaleCoordinates(gpu);
  if (gpu->sim.constraints > 0) {
    kPressureScaleConstraintCoordinates(gpu);
  }
}

//---------------------------------------------------------------------------------------------
// gpu_ucell_set_: Set new unit cell parameters and recalculate molecular
//                 centers of mass. This routine is currently only needed to
//                 update the unit cell after a successful HREMD exchange.
//                 TODO: This code should eventually be consolidated with the
//                       code in gpu_pressure_scale_. For now it is separate
//                       due to an abundance of caution.
//
// Arguments:
//   ucell:      transformation matrix from fractional coordinates to Cartesian space
//   recip:      inverse of ucell
//   uc_volume:  unit cell volume
//---------------------------------------------------------------------------------------------
extern "C" void gpu_ucell_set_(double ucell[3][3], double recip[3][3], double* uc_volume)
{
  PRINTMETHOD("gpu_ucell_set");
  NTPData* pNTPData = gpu->pbNTPData->_pSysData;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pNTPData->last_recip[j * 3 + i] = pNTPData->recip[j * 3 + i];
      pNTPData->ucell[j*3 + i] = ucell[i][j];
      pNTPData->ucellf[j*3 + i] = ucell[i][j];
      pNTPData->recip[j*3 + i] = recip[i][j];
      pNTPData->recipf[j*3 + i] = recip[i][j];
    }
  }

  // Recalculate a, b, c
  double a = pNTPData->ucell[0];
  double b = sqrt(pNTPData->ucell[1] * pNTPData->ucell[1] +
                  pNTPData->ucell[4] * pNTPData->ucell[4]);
  double c = sqrt(pNTPData->ucell[2] * pNTPData->ucell[2] +
                  pNTPData->ucell[5] * pNTPData->ucell[5] +
                  pNTPData->ucell[8] * pNTPData->ucell[8]);

  // Recalculate nonbond skin
  double invU[9];
  invU[0] = ucell[0][0];
  invU[1] = ucell[1][0];
  invU[2] = ucell[2][0];
  invU[3] = ucell[0][1];
  invU[4] = ucell[1][1];
  invU[5] = ucell[2][1];
  invU[6] = ucell[0][2];
  invU[7] = ucell[1][2];
  invU[8] = ucell[2][2];
  double cdepth[3];
  HessianNorms(invU, cdepth);
  double skin = cdepth[0]/gpu->sim.xcells - gpu->sim.cut;
  double yskin = cdepth[1]/gpu->sim.ycells - gpu->sim.cut;
  double zskin = cdepth[2]/gpu->sim.zcells - gpu->sim.cut;

  // The nonbond_skin attribute of the cudaSimulation struct is not updated because
  // it is not actually referenced during the pairlist margin check: the quantities
  // derived from it in gpu->sim.NTPData are.  The entire cudaSimulation is only
  // ported back up to the GPU in the MC barostat case.
  if (yskin < skin) {
    skin = yskin;
  }
  if (zskin < skin) {
    skin = zskin;
  }

  // Check if skin has fallen below acceptable threshold
  if (skin <= 0.5) {
    printf("ERROR: Calculation halted.");
    printf("  Periodic box dimensions have changed too much from their initial values.\n");
    printf("  Your system density has likely changed by a large amount, probably from\n");
    printf("  starting the simulation from a structure a long way from equilibrium.\n");
    printf("\n");
    printf("  [Although this error can also occur if the simulation has blown up for some "
           "reason]\n");
    printf("\n");
    printf("  The GPU code does not automatically reorganize grid cells and thus you\n");
    printf("  will need to restart the calculation from the previous restart file.\n");
    printf("  This will generate new grid cells and allow the calculation to continue.\n");
    printf("  It may be necessary to repeat this restarting multiple times if your system\n");
    printf("  is a long way from an equilibrated density.\n");
    printf("\n");
    printf("  Alternatively you can run with the CPU code until the density has converged\n");
    printf("  and then switch back to the GPU code.\n");
    printf("\n");
    exit(-1);
  }
  skin *= 0.99;
  pNTPData->one_half_nonbond_skin_squared  = (gpu->sim.skinPermit * skin);
  pNTPData->one_half_nonbond_skin_squared *= pNTPData->one_half_nonbond_skin_squared;
  pNTPData->cutPlusSkin2 = (gpu->sim.cut + skin) * (gpu->sim.cut + skin);
  gpu->pbNTPData->Upload();

  if (gpu->sim.ntp > 0) {
    // Handle anything related to pressure.
    // Handle MC barostat
    if (gpu->sim.barostat != 1) {
      kClearSoluteCOM(gpu);
      gpu->sim.af    = (PMEFloat)a;
      gpu->sim.bf    = (PMEFloat)b;
      gpu->sim.cf    = (PMEFloat)c;
      gpu->sim.a     = (PMEDouble)gpu->sim.af;
      gpu->sim.b     = (PMEDouble)gpu->sim.bf;
      gpu->sim.c     = (PMEDouble)gpu->sim.cf;
      gpu->sim.xcell = (PMEFloat)(gpu->sim.a / (PMEDouble)gpu->sim.xcells);
      gpu->sim.ycell = (PMEFloat)(gpu->sim.b / (PMEDouble)gpu->sim.ycells);
      gpu->sim.zcell = (PMEFloat)(gpu->sim.c / (PMEDouble)gpu->sim.zcells);
      gpu->sim.pi_vol_inv = (PMEFloat)(1.0 / (PI * *uc_volume));
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          gpu->sim.ucellf[i][j] = ucell[j][i];
          gpu->sim.recipf[i][j] = recip[j][i];
          gpu->sim.ucell[i][j] = gpu->sim.ucellf[i][j];
          gpu->sim.recip[i][j] = gpu->sim.recipf[i][j];
        }
      }
      if (gpu->sim.is_orthog) {
        for (int i = 0; i < NEIGHBOR_CELLS; i++) {
          gpu->sim.cellOffset[i][0] = (PMEDouble)cellOffset[i][0] * gpu->sim.xcell;
          gpu->sim.cellOffset[i][1] = (PMEDouble)cellOffset[i][1] * gpu->sim.ycell;
          gpu->sim.cellOffset[i][2] = (PMEDouble)cellOffset[i][2] * gpu->sim.zcell;
        }
      }
      gpuCopyConstants();
    }
    // Recalculate COM of everything to ensure correct virial.
    kCalculateCOM(gpu);
  }
}

//---------------------------------------------------------------------------------------------
// gpu_final_gb_setup_: finalize setup for Generalized Born simulations on the GPU
//
// Arguments:
//   igb:          setting for the type of Generalized Born (up to 10 so far)
//   alpb:
//   saltcon:      salt concentration
//   rgbmax:       cutoff for interactions that affect Generalized Born radius computations
//   gb_neckcut:   default value is 2.8, the distance between neck probe and atom centers at
//                 which the solvent-excluded volume will be defined
//   gb_fs_max:    maximum of the array atm_gx_fs
//   atm_gb_radii: the Generalized Born radii of all atoms
//   atm_gb_fs:    the atomic screening factors for GB calculations
//   neckMaxVal:   maximum value for the neck probe-to-atom distance (see the equation in
//                 gb_ene.F90 "Derivative of neck with respect to dij...")
//   neckMaxPos:   not sure, but related to neckMaxVal
//---------------------------------------------------------------------------------------------
extern "C" void gpu_final_gb_setup_(unsigned int* igb, unsigned int* alpb, double* saltcon,
                                    double* rgbmax, double *gb_neckcut, double* gb_fs_max,
                                    double atm_gb_radii[], double atm_gb_fs[],
                                    double neckMaxVal[][21], double neckMaxPos[][21]
                                    )
{
  PRINTMETHOD("gpu_final_gb_setup");
  static const PMEDouble ta = 1.0 / 3.0;
  static const PMEDouble tb = 2.0 / 5.0;
  static const PMEDouble tc = 3.0 / 7.0;
  static const PMEDouble td = 4.0 / 9.0;
  static const PMEDouble tdd = 5.0 / 11.0;
  static const PMEDouble te = 4.0 / 3.0;
  static const PMEDouble tf = 12.0 / 5.0;
  static const PMEDouble tg = 24.0 / 7.0;
  static const PMEDouble th = 40.0 / 9.0;
  static const PMEDouble thh = 60.0 / 11.0;
  static const PMEDouble alpb_alpha = 0.571412;

  // Count atom types
  unsigned int ntypes = 1;
  PMEDouble* gb_radii = new PMEDouble[gpu->sim.atoms];
  PMEDouble* gb_fs    = new PMEDouble[gpu->sim.atoms];
  unsigned int* gb_type = new unsigned int[gpu->sim.atoms];
  gpu->sim.igb = *igb;
  if (gpu->sim.igb != 6) {
    gpu->sim.alpb = *alpb;
    gpu->sim.saltcon = *saltcon;
    gpu->sim.rgbmax = *rgbmax;
  }

  // Set up simulation parameters
  if (gpu->sim.igb != 6) {
    gpu->sim.gb_alpha = 0.0;
    gpu->sim.gb_beta = 0.0;
    gpu->sim.gb_gamma = 0.0;
    gpu->sim.gb_fs_max = 0.0;
    gpu->sim.gb_kappa = 0.0;
    gpu->sim.gb_kappa_inv = 0.0;
    gpu->sim.extdiel_inv = 0.0;
  }
  gpu->sim.intdiel_inv = 0.0;
  gpu->sim.offset     = 0.09;
  if (gpu->sim.igb == 1) {
    gpu->sim.gb_alpha = 1.0;
    gpu->sim.gb_beta = 0.0;
    gpu->sim.gb_gamma = 0.0;
  }
  if (gpu->sim.igb == 2) {

    // Use our best guesses for Onufriev/Case GB  (GB^OBC I):
    gpu->sim.gb_alpha = 0.8;
    gpu->sim.gb_beta = 0.0;
    gpu->sim.gb_gamma = 2.909125;
  }
  else if (gpu->sim.igb == 5) {

    // Use our second best guesses for Onufriev/Case GB (GB^OBC II):
    gpu->sim.gb_alpha = 1.0;
    gpu->sim.gb_beta = 0.80;
    gpu->sim.gb_gamma = 4.850;
  }
  else if (gpu->sim.igb == 7) {
    gpu->sim.gb_alpha = 1.09511284;
    gpu->sim.gb_beta = 1.90792938;
    gpu->sim.gb_gamma = 2.50798245;
    gpu->sim.gb_neckscale = 0.361825;
    gpu->sim.gb_neckoffset = 1.0 - gpu->sim.offset;
  }
  else if (gpu->sim.igb == 8) {
    gpu->sim.gb_neckscale = 0.826836;
    gpu->sim.offset = 0.195141;
    gpu->sim.gb_neckoffset = 1.0 - gpu->sim.offset;
  }
  gpu->sim.gb_neckcut = *gb_neckcut + 2.0 * gpu->sim.offset;

  if ((gpu->sim.igb == 7) || (gpu->sim.igb == 8)) {
    delete gpu->pbNeckMaxValPos;
    gpu->pbNeckMaxValPos = new GpuBuffer<PMEFloat2>(441);
    for (int i = 0; i < 21; i++) {
      for (int j = 0; j < 21; j++) {
        gpu->pbNeckMaxValPos->_pSysData[i * 21 + j].x = neckMaxVal[j][i];
        gpu->pbNeckMaxValPos->_pSysData[i * 21 + j].y = neckMaxPos[j][i];
      }
    }
    gpu->pbNeckMaxValPos->Upload();
    gpu->sim.pNeckMaxValPos = gpu->pbNeckMaxValPos->_pDevData;
  }

  if (gpu->sim.saltcon >= 0.0) {
    gpu->sim.gb_kappa = 0.73 * sqrt(0.108060 * gpu->sim.saltcon);
  }

  if (gpu->sim.alpb == 0) {

    // Standard Still's GB
    if (gpu->sim.igb != 6) {
      gpu->sim.extdiel_inv = 1.0 / gpu->sim.extdiel;
    }
    gpu->sim.intdiel_inv = 1.0 / gpu->sim.intdiel;
  }
  else {

    // Sigalov Onufriev ALPB (epsilon-dependent GB):
    gpu->sim.alpb_beta = alpb_alpha * (gpu->sim.intdiel / gpu->sim.extdiel);
    gpu->sim.extdiel_inv = 1.0 / (gpu->sim.extdiel * (1.0 + gpu->sim.alpb_beta));
    gpu->sim.intdiel_inv = 1.0 / (gpu->sim.intdiel * (1.0 + gpu->sim.alpb_beta));
    gpu->sim.one_arad_beta = gpu->sim.alpb_beta / gpu->sim.arad;
    if (gpu->sim.gb_kappa != 0.0) {
      gpu->sim.gb_kappa_inv = 1.0 / gpu->sim.gb_kappa;
    }
  }
  gb_radii[0] = atm_gb_radii[0];
  gb_fs[0] = atm_gb_fs[0];
  for (int i = 0; i < gpu->sim.atoms; i++) {
    int j;
    for (j = 0; j < ntypes; j++) {
      if (((PMEDouble)atm_gb_radii[i] == gb_radii[j]) &&
          ((PMEDouble)atm_gb_fs[i] == gb_fs[j])) {
        break;
      }
    }
    gb_type[i] = j;

    // Add new type if discovered
    if (j == ntypes) {
      gb_radii[j] = atm_gb_radii[i];
      gb_fs[j] = atm_gb_fs[i];
      if (gpu->sim.gb_fs_max < gb_fs[j]) {
        gpu->sim.gb_fs_max = gb_fs[j];
      }
      ntypes++;
    }
  }

#ifdef GVERBOSE
  printf("%d GB atom types encountered.\n", ntypes);
  for (int j = 0; j < ntypes; j++) {
    printf("%5d %16.8f %16.8f\n", j, gb_radii[j], gb_fs[j]);
  }
#endif
  gpu->sim.rgbmax1i = 1.0 / gpu->sim.rgbmax;
  gpu->sim.rgbmax2i = gpu->sim.rgbmax1i * gpu->sim.rgbmax1i;
  gpu->sim.rgbmaxpsmax2 = (gpu->sim.rgbmax + gpu->sim.gb_fs_max) *
                          (gpu->sim.rgbmax + gpu->sim.gb_fs_max);
  gpuCopyConstants();

  // Lookup table is no longer generated.  See older revisions
  // in the repository for code to generated the table.

  // Free allocated memory
  delete[] gb_radii;
  delete[] gb_fs;
  delete[] gb_type;
}


//---------------------------------------------------------------------------------------------
// gpu_gb_gbion_setup_: setup for the special gbion!=0 case
//
// Arguments:
//   gi_coef, gi_coef_1, gi_coef_2, gb_neckscale_ion_1, gb_neckscale_ion_2,
//   intdiel_ion_1, intdiel_ion_2,
//   atom_ion_mask[],
//   gi_coef_1_p, gi_coef_1_n, gi_coef_2_pp, gi_coef_2_pn, gi_coef_2_nn,
//   gb_neckscale_ion_1_p, gb_neckscale_ion_1_n, &
//   gb_neckscale_ion_2_pp, gb_neckscale_ion_2_pn, gb_neckscale_ion_2_nn,
//   intdiel_ion_1_p, intdiel_ion_1_n,
//   intdiel_ion_2_pp, intdiel_ion_2_pn, intdiel_ion_2_nn
//---------------------------------------------------------------------------------------------
extern "C" void gpu_gb_gbion_setup_(unsigned int* gbion, double* gi_coef_1, double* gi_coef_2,
                                    double* intdiel_ion_1, double* intdiel_ion_2,
                                    double* gb_neckscale_ion_1, double* gb_neckscale_ion_2,
                                    unsigned int atom_ion_mask[],
                                    double* gi_coef_1_p, double* gi_coef_1_n, double* gi_coef_2_pp,
                                    double* gi_coef_2_pn, double* gi_coef_2_nn,
                                    double* gb_neckscale_ion_1_p, double* gb_neckscale_ion_1_n,
                                    double* gb_neckscale_ion_2_pp, double* gb_neckscale_ion_2_pn,
                                    double* gb_neckscale_ion_2_nn, double* intdiel_ion_1_p, double* intdiel_ion_1_n,
                                    double* intdiel_ion_2_pp, double* intdiel_ion_2_pn, double* intdiel_ion_2_nn)
{
  PRINTMETHOD("gpu_gb_gbion_setup");

  // GB w. explicit ions: gbion (0:off 2:on)
  gpu->sim.gbion = *gbion;
  if (gpu->sim.gbion==2) {
    gpu->sim.gi_coef_1 = *gi_coef_1;
    gpu->sim.gi_coef_2 = *gi_coef_2;
    gpu->sim.intdiel_ion_1 = *intdiel_ion_1;
    gpu->sim.intdiel_ion_2 = *intdiel_ion_2;
    gpu->sim.gb_neckscale_ion_1 = *gb_neckscale_ion_1;
    gpu->sim.gb_neckscale_ion_2 = *gb_neckscale_ion_2;
  }
  else if (gpu->sim.gbion==3)  {
    gpu->sim.gi_coef_1_p = *gi_coef_1_p;
    gpu->sim.gi_coef_1_n = *gi_coef_1_n;
    gpu->sim.gi_coef_2_pp = *gi_coef_2_pp;
    gpu->sim.gi_coef_2_pn = *gi_coef_2_pn;
    gpu->sim.gi_coef_2_nn = *gi_coef_2_nn;
    gpu->sim.intdiel_ion_1_p = *intdiel_ion_1_p;
    gpu->sim.intdiel_ion_1_n = *intdiel_ion_1_n;
    gpu->sim.intdiel_ion_2_pp = *intdiel_ion_2_pp;
    gpu->sim.intdiel_ion_2_pn = *intdiel_ion_2_pn;
    gpu->sim.intdiel_ion_2_nn = *intdiel_ion_2_nn;
    gpu->sim.gb_neckscale_ion_1_p = *gb_neckscale_ion_1_p;
    gpu->sim.gb_neckscale_ion_1_n = *gb_neckscale_ion_1_n;
    gpu->sim.gb_neckscale_ion_2_pp = *gb_neckscale_ion_2_pp;
    gpu->sim.gb_neckscale_ion_2_pn = *gb_neckscale_ion_2_pn;
    gpu->sim.gb_neckscale_ion_2_nn = *gb_neckscale_ion_2_nn;
  }

  // upload atom_ion_mask, based on "gpu_upload_charges_gb_cph"
  for (int i = 0; i < gpu->sim.atoms; i++) {
    gpu->pbAtomIonMask->_pSysData[i] = atom_ion_mask[i];
  }
  for (int i = gpu->sim.atoms; i < gpu->sim.paddedNumberOfAtoms; i++) {
    gpu->pbAtomIonMask->_pSysData[i] = 0;
  }
  gpu->pbAtomIonMask->Upload();

  gpuCopyConstants();
}


//---------------------------------------------------------------------------------------------
// gpu_gb_igb8_setup_: setup for the special igb == 8 case
//
// Arguments:
//   gb_alpha:   these are arrays when igb == 8, which is why we have this special setup
//   gb_beta:    function--each value for different atom types could be different.
//   gb_gamma:
//---------------------------------------------------------------------------------------------
extern "C" void gpu_gb_igb8_setup_(double gb_alpha[], double gb_beta[], double gb_gamma[])
{
  PRINTMETHOD("gpu_gb_igb8_setup");

  // Delete existing alpha, beta, and gamma arrays
  delete gpu->pbGBAlphaBetaGamma;

  // Allocate new arrays
  gpu->pbGBAlphaBetaGamma = new GpuBuffer<PMEFloat>(gpu->sim.stride3);
  gpu->sim.pgb_alpha = gpu->pbGBAlphaBetaGamma->_pDevData;
  gpu->sim.pgb_beta = gpu->sim.pgb_alpha + gpu->sim.stride;
  gpu->sim.pgb_gamma = gpu->sim.pgb_alpha + gpu->sim.stride2;

  // Fill arrays
  PMEFloat* pAlpha = gpu->pbGBAlphaBetaGamma->_pSysData;
  PMEFloat* pBeta = gpu->pbGBAlphaBetaGamma->_pSysData + gpu->sim.stride;
  PMEFloat* pGamma = gpu->pbGBAlphaBetaGamma->_pSysData + gpu->sim.stride2;
  for (int i = 0; i < gpu->sim.atoms; i++) {
    pAlpha[i] = gb_alpha[i];
    pBeta[i] = gb_beta[i];
    pGamma[i] = gb_gamma[i];
  }
  gpu->pbGBAlphaBetaGamma->Upload();
  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gpu_build_threadblock_work_list_: build a work plan for various thread blocks on the GPU.
//
// Arguments:
//   numex:
//   natex:
//---------------------------------------------------------------------------------------------
extern "C" void gpu_build_threadblock_work_list_(int numex[], int natex[])
{
  PRINTMETHOD("gpu_build_threadblock_work_list");
  const unsigned int dim = (gpu->sim.paddedNumberOfAtoms + (GRID - 1)) / GRID;
  const unsigned int cells = dim * (dim + 1) / 2;

  // Delete existing data
  delete gpu->pbWorkUnit;
  delete gpu->pbExclusion;
  delete gpu->pbGBPosition;

  // Determine number of work units with exclusions so we can place them at the head
  // of the list.  Generate symmetric exclusion list based on assymmetric -1 padded
  // list which is capped at 2x the previous list size, but is likely smaller due to
  // superfluous -1 entries
  unsigned int totalExclusions = 0;
  unsigned int* pExclCount = new unsigned int[gpu->sim.paddedNumberOfAtoms + 1];
  for (int i = 0; i < gpu->sim.atoms; i++) {
    totalExclusions += numex[i];
    pExclCount[i] = 0;
  }
  for (int i = gpu->sim.atoms; i < gpu->sim.paddedNumberOfAtoms; i++) {
    pExclCount[i] = 0;
  }

  // Count Exclusions and filter out -1 entries
  unsigned int* pExclList = new unsigned int[totalExclusions * 2];
  unsigned int* pExclOffset = new unsigned int[gpu->sim.paddedNumberOfAtoms + 1];
  unsigned int offset = 0;
  for (unsigned int x = 0; x < gpu->sim.atoms; x++) {
    for (unsigned int j = offset; j < offset + numex[x]; j++) {
      if (natex[j] > 0) {
        unsigned int y = natex[j] - 1;
        pExclCount[x]++;
        pExclCount[y]++;
      }
    }
    offset += numex[x];
  }

  // Calculate symmetric exclusion offsets
  unsigned int* pExclCounter = new unsigned int[gpu->sim.paddedNumberOfAtoms + 1];
  offset = 0;
  for (int i = 0; i < gpu->sim.atoms; i++) {
    pExclOffset[i]  = offset;
    pExclCounter[i] = offset;
    offset += pExclCount[i];
  }
  for (int i = gpu->sim.atoms; i <= gpu->sim.paddedNumberOfAtoms; i++) {
    pExclOffset[i]  = offset;
    pExclCounter[i] = offset;
  }

  // Now regenerate exclusions
  offset = 0;
  for (int x = 0; x < gpu->sim.atoms; x++) {
    for (int j = offset; j < offset + numex[x]; j++) {
      if (natex[j] > 0) {
        unsigned int y = natex[j] - 1;
        pExclList[pExclCounter[x]++] = y;
        pExclList[pExclCounter[y]++] = x;
      }
    }
    offset += numex[x];
  }

  // Partition nonbond calculation amongst nodes
#ifdef MPI
  unsigned int cellCounter = 0;
  unsigned int totalWorkUnits = (dim * (dim + 1)) / 2;
  unsigned int tileCounter = 0;
  unsigned int minTile = (gpu->gpuID * totalWorkUnits) / gpu->nGpus;
  unsigned int maxTile = ((gpu->gpuID + 1) * totalWorkUnits) / gpu->nGpus;
  if (gpu->gpuID == 0) {
    gpu->bCalculateLocalForces = true;
  }
  else {
    gpu->bCalculateLocalForces = false;
  }
#endif

  unsigned int excludedWorkUnits = 0;
  for (unsigned int y = 0; y < dim; y++) {
    for (unsigned int x = y; x < dim; x++) {
#ifdef MPI
      if ((tileCounter >= minTile) && (tileCounter < maxTile)) {
#endif
        unsigned int xstart = x * GRID;
        unsigned int ystart = y * GRID;
        unsigned int xend = xstart + GRID;
        unsigned int yend = ystart + GRID;
        bool excluded = false;
        for (int i = pExclOffset[xstart]; i < pExclOffset[xend]; i++) {
          if ((pExclList[i] >= ystart) && (pExclList[i] < yend)) {
            excluded = true;
            break;
          }
        }
        if (excluded) {
          excludedWorkUnits++;
        }
#ifdef MPI
        cellCounter++;
      }
      tileCounter++;
#endif
    }
  }
#ifdef MPI
  gpu->sim.workUnits = cellCounter;
  tileCounter = 0;
#else
  gpu->sim.workUnits = cells;
#endif

  // Decrease thread count for extra small molecules to spread computation
  // across entire chip
#ifdef AMBER_PLATFORM_AMD
  gpu->GBBornRadiiBlocks = (gpu->sim.workUnits /
                           (gpu->GBBornRadiiThreadsPerBlock / GRID)) + 1;
  gpu->GBNonbondEnergy1Blocks = (gpu->sim.workUnits /
                                (gpu->GBNonbondEnergy1ThreadsPerBlock / GRID)) + 1;
  gpu->GBNonbondEnergy2Blocks = (gpu->sim.workUnits /
                                (gpu->GBNonbondEnergy2ThreadsPerBlock / GRID)) + 1;
#endif
  int balancedWorkBlock = GRID * (1 + gpu->sim.workUnits / gpu->blocks);
  int activeWorkUnits = (gpu->blocks * gpu->GBBornRadiiThreadsPerBlock) / GRID;
  if (activeWorkUnits > gpu->sim.workUnits + gpu->blocks) {
#ifdef AMBER_PLATFORM_AMD
    gpu->GBBornRadiiBlocks *= (gpu->GBBornRadiiThreadsPerBlock / GRID);
#endif
    gpu->GBBornRadiiThreadsPerBlock = balancedWorkBlock;
  }
  activeWorkUnits = (gpu->blocks * gpu->GBNonbondEnergy1ThreadsPerBlock) / GRID;
  if (activeWorkUnits > gpu->sim.workUnits + gpu->blocks) {
#ifdef AMBER_PLATFORM_AMD
    gpu->GBNonbondEnergy1Blocks *= (gpu->GBNonbondEnergy1ThreadsPerBlock / GRID);
#endif
    gpu->GBNonbondEnergy1ThreadsPerBlock = balancedWorkBlock;
  }
  activeWorkUnits = (gpu->blocks * gpu->GBNonbondEnergy2ThreadsPerBlock) / GRID;
  if (activeWorkUnits > gpu->sim.workUnits + gpu->blocks) {
#ifdef AMBER_PLATFORM_AMD
    gpu->GBNonbondEnergy2Blocks *= (gpu->GBNonbondEnergy2ThreadsPerBlock / GRID);
#endif
    gpu->GBNonbondEnergy2ThreadsPerBlock = balancedWorkBlock;
  }

  // Build work unit list
  gpu->sim.excludedWorkUnits = excludedWorkUnits;
  unsigned int excludedOffset = 0;
  unsigned int unexcludedOffset = excludedWorkUnits;
  unsigned int exclusions = 0;
  GpuBuffer<PMEMask>* pbExclusion = new GpuBuffer<PMEMask>(GRID * excludedWorkUnits);
  PMEMask* pExclusion = pbExclusion->_pSysData;
  GpuBuffer<unsigned int>* pbWorkUnit = new GpuBuffer<unsigned int>(gpu->sim.workUnits);
  unsigned int* pWorkUnit = pbWorkUnit->_pSysData;
  gpu->pbGBPosition = new GpuBuffer<unsigned int>(3);
  for (unsigned int y = 0; y < dim; y++) {
    for (unsigned int x = y; x < dim; x++) {
#ifdef MPI
      if ((tileCounter >= minTile) && (tileCounter < maxTile)) {
#endif
        // Check for exclusions
        unsigned int xstart = x * GRID;
        unsigned int ystart = y * GRID;
        unsigned int xend = xstart + GRID;
        unsigned int yend = ystart + GRID;
        bool excluded = false;
        for (int i = pExclOffset[xstart]; i < pExclOffset[xend]; i++) {
          if ((pExclList[i] >= ystart) && (pExclList[i] < yend)) {
            excluded = true;
            break;
          }
        }

        // Add exclusions if present
        if (excluded) {
          // Create exclusion masks
          PMEMask excl[GRID];
          for (int i = 0; i < GRID; i++) {
            excl[i] = WARP_MASK;
          }
          for (int i = 0; i < GRID; i++) {
            unsigned int x = xstart + i;
            for (int j = pExclOffset[x]; j < pExclOffset[x + 1]; j++) {
              unsigned int y = pExclList[j];
              if ((y >= ystart) && (y < yend)) {
                excl[i] ^= ((PMEMask) 1) << (y - ystart);
              }
            }
          }

          // Skip padded atoms
          if (xend > gpu->sim.atoms) {
            for (int i = gpu->sim.atoms - xstart; i < GRID; i++) {
               excl[i] = 0;
            }
          }
          if (yend > gpu->sim.atoms) {
            unsigned int offset = yend - gpu->sim.atoms;
            for (int i = 0; i < GRID; i++) {
              excl[i] = (excl[i] << offset) >> offset;
            }
          }

          // Post-process exclusion masks
          for (int i = 0; i < GRID; i++) {
            unsigned int offset = i;
            if (xstart == ystart) {
              offset++;
            }
            excl[i] = (excl[i] >> offset) | (excl[i] << (GRID - offset));
            pExclusion[exclusions++]= excl[i];
          }
          pWorkUnit[excludedOffset] = (x << 17) | (y << 2);
          pWorkUnit[excludedOffset] |= 0x1;
          excludedOffset++;
        }
        else {
          pWorkUnit[unexcludedOffset] = (x << 17) | (y << 2);
          unexcludedOffset++;
        }
#ifdef MPI
      }
      tileCounter++;
#endif
    }
  }

  // Delete temporary data
  delete[] pExclOffset;
  delete[] pExclCount;
  delete[] pExclCounter;
  delete[] pExclList;

  // Set up GPU pointers and constants
  gpu->pbWorkUnit     = pbWorkUnit;
  gpu->sim.pWorkUnit  = pbWorkUnit->_pDevData;
  gpu->pbExclusion    = pbExclusion;
  gpu->sim.pExclusion = pbExclusion->_pDevData;
  gpu->sim.pGBBRPosition = gpu->pbGBPosition->_pDevData;
  gpu->sim.pGBNB1Position = gpu->pbGBPosition->_pDevData + 1;
  gpu->sim.pGBNB2Position = gpu->pbGBPosition->_pDevData + 2;
  gpu->sim.GBTotalWarps[0] = (gpu->GBBornRadiiThreadsPerBlock * gpu->GBBornRadiiBlocks) / GRID;
  gpu->sim.GBTotalWarps[1] = (gpu->GBNonbondEnergy1ThreadsPerBlock *
                              gpu->GBNonbondEnergy1Blocks) / GRID;
  gpu->sim.GBTotalWarps[2] = (gpu->GBNonbondEnergy2ThreadsPerBlock *
                              gpu->GBNonbondEnergy2Blocks) / GRID;
  pbWorkUnit->Upload();
  pbExclusion->Upload();
  gpuCopyConstants();

  return;
}

//---------------------------------------------------------------------------------------------
// gpu_final_setup_: just announce that there was a final setup?  This function appears to be
//                   a stub.
//
// This function is not called anywhere in the code at present.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_final_setup_()
{
  PRINTMETHOD("gpu_final_setup");
}

//---------------------------------------------------------------------------------------------
// gpu_update_phmd_: Langevin thermostatting of titration variables
//
// Arguments:
//   dtx:       time step
//---------------------------------------------------------------------------------------------
extern "C" void gpu_update_phmd_(double* dtx)
{
  kUpdatePHMD(*dtx, gpu);
  cudaDeviceSynchronize();
}
//---------------------------------------------------------------------------------------------
// gpu_update_: update critical parameters for Langevin thermostating
//
// Arguments:
//   dt:        time step
//   temp0:     target temperature
//   gamma_ln:  Langevin collision frequency
//---------------------------------------------------------------------------------------------
extern "C" void gpu_update_(double* dt, double* temp0, double* gamma_ln)
{
  PRINTMETHOD("gpu_update");
  if (gpu->ischeme == 0) {
      kUpdate(gpu, *dt, *temp0, *gamma_ln);
  }
  else {
      kUpdateMiddle(gpu, *dt, *temp0, *gamma_ln);
  }
  if (gpu->sim.ti_mode > 0 && gpu->sim.TIPaddedLinearAtmCnt > 0)
  {
    kAFEExchangeVel(gpu);
    kAFEExchangeCrd(gpu);
    kAFEExchangeFrc(gpu);
  }
}

//---------------------------------------------------------------------------------------------
// gpu_relaxmd_update_: update critical parameters for relaxed MD
//
// Arguments:
//   dt:        time step
//   temp0:     target temperature
//   gamma_ln:  Langevin collision frequency
//---------------------------------------------------------------------------------------------
extern "C" void gpu_relaxmd_update_(double* dt, double* temp0, double* gamma_ln)
{
  PRINTMETHOD("gpu_relaxmd_update");
  kRelaxMDUpdate(gpu, *dt, *temp0, *gamma_ln);
}



//---------------------------------------------------------------------------------------------
// gpu_sgld_update_: update critical parameters for SGLD
//
// Arguments:
//   dt:        time step
//   temp0:     target temperature
//   gamma_ln:  Langevin collision frequency
//---------------------------------------------------------------------------------------------
extern "C" void gpu_sgld_update_(double* dt, double* temp0, double* gamma_ln, double* gamma_sg, double com0sg[], double com1sg[], double com2sg[])
{
PRINTMETHOD("gpu_SGLD_update");

    kUpdateSGLD(gpu, *dt, *temp0, *gamma_ln, *gamma_sg, com0sg, com1sg, com2sg);

}


extern "C" void gpu_sgld_setup_(int* isgld, int* nsgatom, int sgatoms[], double* dt, double *tsgset, double* sgmsum, double* tsgavg,double* tsgavp, double* sgft, double* sgff, double* sgfg,double* fsgldg, double* tempsg, double sg_x0[][3], double sg_x1[][3], double sg_x2[][3], double sg_r[][3],double sg_fps[], double sg_pps[])
{
PRINTMETHOD("gpu_sgld_setup");
    // Determine what type of SGLD is being used
    gpu->sim.isgld                                   = *isgld;
    gpu->sim.nsgatom                                 = *nsgatom;
    // Allocate GPU data
    // Set up SGLD parameters
    gpu->sim.sgavg                                  = *dt/(*tsgavg);
    gpu->sim.sgavp                                  = *dt/(*tsgavp);
    gpu->sim.sgft                                   = *sgft;
    gpu->sim.sgff                                   = *sgff;
    gpu->sim.sgfg                                   = *sgfg;
    gpu->sim.fsgldg                                 = *fsgldg;
    gpu->sim.tempsg                                 = *tempsg;
    gpu->sim.tsgset                                 = *tsgset;
    gpu->sim.sgmsum                                 = *sgmsum;

  
    delete gpu->pbSGLDAverageBuffer;

     //  Allocate SGLD  buffer
    if ( gpu->sim.isgld > 0 )
    {
    if (gpu->bCanMapHostMemory)
        gpu->pbSGLDAverageBuffer                  = new GpuBuffer<SGLDAverage>(gpu->blocks, false, true);
    else
        gpu->pbSGLDAverageBuffer                  = new GpuBuffer<SGLDAverage>(gpu->blocks);

    gpu->sim.pSGLDAverage                         = gpu->pbSGLDAverageBuffer->_pDevData;
    }

    if (gpu->bNeighborList && (gpu->pbImageIndex != NULL))
    {
        delete  gpu->pbImageWsg;       // SGLD Image
        delete  gpu->pbImageX0sg;       // SGLD Image
        delete  gpu->pbImageX1sg;       // SGLD Image
        delete  gpu->pbImageX2sg;       // SGLD Image
        delete  gpu->pbImageRsg;       // SGLD Image
        delete  gpu->pbImageFPsg;       // SGLD Image
        delete  gpu->pbImagePPsg;       // SGLD Image

        gpu->pbImageWsg                         = new GpuBuffer<int>(gpu->sim.stride2);
        gpu->pbImageX0sg                         = new GpuBuffer<double>(2*gpu->sim.stride3);
        gpu->pbImageX1sg                         = new GpuBuffer<double>(2*gpu->sim.stride3);
        gpu->pbImageX2sg                         = new GpuBuffer<double>(2*gpu->sim.stride3);
        gpu->pbImageRsg                          = new GpuBuffer<double>(2*gpu->sim.stride3);
        gpu->pbImageFPsg                         = new GpuBuffer<double>(gpu->sim.stride2);
        gpu->pbImagePPsg                         = new GpuBuffer<double>(gpu->sim.stride2);

        gpu->sim.pImageWsg                      = gpu->pbImageWsg->_pDevData;
        gpu->sim.pImageWsg2                      = gpu->pbImageWsg->_pDevData+ gpu->sim.stride;
        gpu->sim.pImageX0sg                     = gpu->pbImageX0sg->_pDevData;
        gpu->sim.pImageY0sg                     = gpu->pbImageX0sg->_pDevData + gpu->sim.stride;
        gpu->sim.pImageZ0sg                     = gpu->pbImageX0sg->_pDevData + gpu->sim.stride2;
        gpu->sim.pImageX0sg2                     = gpu->pbImageX0sg->_pDevData+ gpu->sim.stride3;
        gpu->sim.pImageY0sg2                     = gpu->pbImageX0sg->_pDevData + gpu->sim.stride4;
        gpu->sim.pImageZ0sg2                     = gpu->pbImageX0sg->_pDevData + gpu->sim.stride2+ gpu->sim.stride3;
        gpu->sim.pImageX1sg                     = gpu->pbImageX1sg->_pDevData;
        gpu->sim.pImageY1sg                     = gpu->pbImageX1sg->_pDevData + gpu->sim.stride;
        gpu->sim.pImageZ1sg                     = gpu->pbImageX1sg->_pDevData + gpu->sim.stride2;
        gpu->sim.pImageX1sg2                     = gpu->pbImageX1sg->_pDevData + gpu->sim.stride3;
        gpu->sim.pImageY1sg2                     = gpu->pbImageX1sg->_pDevData + gpu->sim.stride4;
        gpu->sim.pImageZ1sg2                     = gpu->pbImageX1sg->_pDevData + gpu->sim.stride2 + gpu->sim.stride3;
        gpu->sim.pImageX2sg                     = gpu->pbImageX2sg->_pDevData;
        gpu->sim.pImageY2sg                     = gpu->pbImageX2sg->_pDevData + gpu->sim.stride;
        gpu->sim.pImageZ2sg                     = gpu->pbImageX2sg->_pDevData + gpu->sim.stride2;
        gpu->sim.pImageX2sg2                     = gpu->pbImageX2sg->_pDevData + gpu->sim.stride3;
        gpu->sim.pImageY2sg2                     = gpu->pbImageX2sg->_pDevData + gpu->sim.stride4;
        gpu->sim.pImageZ2sg2                     = gpu->pbImageX2sg->_pDevData + gpu->sim.stride2 + gpu->sim.stride3;
        gpu->sim.pImageRsgX                     = gpu->pbImageRsg->_pDevData;
        gpu->sim.pImageRsgY                     = gpu->pbImageRsg->_pDevData + gpu->sim.stride;
        gpu->sim.pImageRsgZ                     = gpu->pbImageRsg->_pDevData + gpu->sim.stride2;
        gpu->sim.pImageRsgX2                     = gpu->pbImageRsg->_pDevData + gpu->sim.stride3;
        gpu->sim.pImageRsgY2                     = gpu->pbImageRsg->_pDevData + gpu->sim.stride4;
        gpu->sim.pImageRsgZ2                     = gpu->pbImageRsg->_pDevData + gpu->sim.stride2 + gpu->sim.stride3;
        gpu->sim.pImageFPsg                      = gpu->pbImageFPsg->_pDevData;
        gpu->sim.pImageFPsg2                      = gpu->pbImageFPsg->_pDevData+ gpu->sim.stride;
        gpu->sim.pImagePPsg                      = gpu->pbImagePPsg->_pDevData;
        gpu->sim.pImagePPsg2                      = gpu->pbImagePPsg->_pDevData+ gpu->sim.stride;
        if (gpu->bNewNeighborList)
        {
            gpu->pbImageIndex->Download();
            gpu->bNewNeighborList                       = false;
        }
        unsigned int* pImageAtomLookup                  = &(gpu->pbImageIndex->_pSysData[gpu->sim.imageStride * 2]);
        int *pWsg                                       = gpu->pbImageWsg->_pSysData;
        double *pX0sg                                    = gpu->pbImageX0sg->_pSysData;
        double *pX1sg                                    = gpu->pbImageX1sg->_pSysData;
        double *pX2sg                                    = gpu->pbImageX2sg->_pSysData;
        double *pRsg                                    = gpu->pbImageRsg->_pSysData;
        double *pFPsg                                    = gpu->pbImageFPsg->_pSysData;
        double *pPPsg                                    = gpu->pbImagePPsg->_pSysData;
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            int i1                                      = pImageAtomLookup[i];
            pWsg[i1]                                    =  sgatoms[i];
            //pPsg[i1]                                    =  0.0;
            pX0sg[i1]                                    =  sg_x0[i][0];
            pX0sg[i1 + gpu->sim.stride]                  =  sg_x0[i][1];
            pX0sg[i1 + gpu->sim.stride2]                 =  sg_x0[i][2];
            pX1sg[i1]                                    =  sg_x1[i][0];
            pX1sg[i1 + gpu->sim.stride]                  =  sg_x1[i][1];
            pX1sg[i1 + gpu->sim.stride2]                 =  sg_x1[i][2];
            pX2sg[i1]                                    =  sg_x2[i][0];
            pX2sg[i1 + gpu->sim.stride]                  =  sg_x2[i][1];
            pX2sg[i1 + gpu->sim.stride2]                 =  sg_x2[i][2];
            pFPsg[i1]                                    = sg_fps[i];
            pPPsg[i1]                                    = sg_pps[i];
            pRsg[i1]                                     = sg_r[i][0];
            pRsg[i1 + gpu->sim.stride]                   = sg_r[i][1];
            pRsg[i1 + gpu->sim.stride2]                  = sg_r[i][2];
        }
        gpu->pbImageWsg->Upload();
        gpu->pbImageX0sg->Upload();
        gpu->pbImageX1sg->Upload();
        gpu->pbImageX2sg->Upload();
        gpu->pbImageRsg->Upload();
        gpu->pbImageFPsg->Upload();
        gpu->pbImagePPsg->Upload();
    }
    else
    {
        delete  gpu->pbWsg;       // SGLD 
        delete  gpu->pbX0sg;       // SGLD 
        delete  gpu->pbX1sg;       // SGLD
        delete  gpu->pbX2sg;       // SGLD 
        delete  gpu->pbRsg;       // SGLD
        delete  gpu->pbFPsg;       // SGLD 
        delete  gpu->pbPPsg;       // SGLD 

        gpu->pbWsg                          = new GpuBuffer<int>(gpu->sim.stride);
        gpu->pbX0sg                         = new GpuBuffer<double>(gpu->sim.stride3);
        gpu->pbX1sg                         = new GpuBuffer<double>(gpu->sim.stride3);
        gpu->pbX2sg                         = new GpuBuffer<double>(gpu->sim.stride3);
        gpu->pbRsg                          = new GpuBuffer<double>(gpu->sim.stride3);
        gpu->pbFPsg                         = new GpuBuffer<double>(gpu->sim.stride);
        gpu->pbPPsg                         = new GpuBuffer<double>(gpu->sim.stride);

        gpu->sim.pWsg                      = gpu->pbWsg->_pDevData;
        gpu->sim.pX0sg                     = gpu->pbX0sg->_pDevData;
        gpu->sim.pY0sg                     = gpu->pbX0sg->_pDevData + gpu->sim.stride;
        gpu->sim.pZ0sg                     = gpu->pbX0sg->_pDevData + gpu->sim.stride2;
        gpu->sim.pX1sg                     = gpu->pbX1sg->_pDevData;
        gpu->sim.pY1sg                     = gpu->pbX1sg->_pDevData + gpu->sim.stride;
        gpu->sim.pZ1sg                     = gpu->pbX1sg->_pDevData + gpu->sim.stride2;
        gpu->sim.pX2sg                     = gpu->pbX2sg->_pDevData;
        gpu->sim.pY2sg                     = gpu->pbX2sg->_pDevData + gpu->sim.stride;
        gpu->sim.pZ2sg                     = gpu->pbX2sg->_pDevData + gpu->sim.stride2;
        gpu->sim.pRsgX                     = gpu->pbRsg->_pDevData;
        gpu->sim.pRsgY                     = gpu->pbRsg->_pDevData + gpu->sim.stride;
        gpu->sim.pRsgZ                     = gpu->pbRsg->_pDevData + gpu->sim.stride2;
        gpu->sim.pFPsg                     = gpu->pbFPsg->_pDevData;
        gpu->sim.pPPsg                     = gpu->pbPPsg->_pDevData;
        for (int i = 0; i < gpu->sim.atoms; i++)
        {
            gpu->pbWsg->_pSysData[i]                     =  sgatoms[i];
            gpu->pbX0sg->_pSysData[i]                    = sg_x0[i][0];
            gpu->pbX0sg->_pSysData[i + gpu->sim.stride]  = sg_x0[i][1];
            gpu->pbX0sg->_pSysData[i + gpu->sim.stride2] = sg_x0[i][2];
            gpu->pbX1sg->_pSysData[i]                    = sg_x1[i][0];
            gpu->pbX1sg->_pSysData[i + gpu->sim.stride]  = sg_x1[i][1];
            gpu->pbX1sg->_pSysData[i + gpu->sim.stride2] = sg_x1[i][2];
            gpu->pbX2sg->_pSysData[i]                    = sg_x2[i][0];
            gpu->pbX2sg->_pSysData[i + gpu->sim.stride]  = sg_x2[i][1];
            gpu->pbX2sg->_pSysData[i + gpu->sim.stride2] = sg_x2[i][2];
            gpu->pbRsg->_pSysData[i]                    = sg_r[i][0];
            gpu->pbRsg->_pSysData[i + gpu->sim.stride]  = sg_r[i][1];
            gpu->pbRsg->_pSysData[i + gpu->sim.stride2] = sg_r[i][2];
            gpu->pbFPsg->_pSysData[i]                    = sg_fps[i];
            gpu->pbPPsg->_pSysData[i]                    = sg_pps[i];
        }
        gpu->pbWsg->Upload();
        gpu->pbX0sg->Upload();
        gpu->pbX1sg->Upload();
        gpu->pbX2sg->Upload();
        gpu->pbRsg->Upload();
        gpu->pbFPsg->Upload();
        gpu->pbPPsg->Upload();
    } 

    gpuCopyConstants();

    return;
}


extern "C" void gpu_sgld_rndf_(double* sgrndf)
{
PRINTMETHOD("gpu_sgld_rndf");

   gpu->sim.sgrndf = ((double) (*sgrndf) ) * 0.001987192 *0.5;
   gpuCopyConstants();

}


extern "C" void gpu_calculate_sgld_averages_(double* dt, double* temp0, double* gamma_ln, double* templf, double* temphf, double* sggamma, double com0[], double com1[], double com2[])
{
PRINTMETHOD("gpu_calculate_SGLD_averages");

     double tsgfac =gpu->sim.sgavg/(*dt)/20.455;
     kCalculateSGLDaverages(gpu, *dt);
     
     
    cudaDeviceSynchronize();

    
    if (!(gpu->bCanMapHostMemory))
        gpu->pbSGLDAverageBuffer->Download();


    
    PMEFloat E[12];
    for (int i = 0; i < 12; i++)E[i]       = 0.0f;     
    for (int i = 0; i < gpu->blocks; i++)
    {
        E[0]                                   += gpu->pbSGLDAverageBuffer->_pSysData[i].SGLDavg.EKT;
        E[1]                                   += gpu->pbSGLDAverageBuffer->_pSysData[i].SGLDavg.EKSG;
        E[2]                                   += gpu->pbSGLDAverageBuffer->_pSysData[i].SGLDavg.GAMSG;
        E[3]                                   += gpu->pbSGLDAverageBuffer->_pSysData[i].SGLDavg.COM0SG[0];
        E[4]                                   += gpu->pbSGLDAverageBuffer->_pSysData[i].SGLDavg.COM0SG[1];
        E[5]                                   += gpu->pbSGLDAverageBuffer->_pSysData[i].SGLDavg.COM0SG[2];
        E[6]                                   += gpu->pbSGLDAverageBuffer->_pSysData[i].SGLDavg.COM1SG[0];
        E[7]                                   += gpu->pbSGLDAverageBuffer->_pSysData[i].SGLDavg.COM1SG[1];
        E[8]                                   += gpu->pbSGLDAverageBuffer->_pSysData[i].SGLDavg.COM1SG[2];
        E[9]                                   += gpu->pbSGLDAverageBuffer->_pSysData[i].SGLDavg.COM2SG[0];
        E[10]                                   += gpu->pbSGLDAverageBuffer->_pSysData[i].SGLDavg.COM2SG[1];
        E[11]                                   += gpu->pbSGLDAverageBuffer->_pSysData[i].SGLDavg.COM2SG[2];
    }

    
    double sgavp0=1.0-gpu->sim.sgavp;
    double tempt = gpu->sim.tsgset;
    *templf=sgavp0*(*templf)+ gpu->sim.sgavp*(PMEDouble)E[1]*gpu->sim.tsgset/(PMEDouble)E[0];           // templfi=eksg/sgrndf/botz;
    *temphf=tempt-(*templf);           // templfi=tsgset-templf;
    *sggamma=E[2]/(double)(gpu->sim.nsgatom);                // sumgam;
    com0[0]=E[3]/(double)(gpu->sim.sgmsum);                // COM;
    com0[1]=E[4]/(double)(gpu->sim.sgmsum);                // COM;
    com0[2]=E[5]/(double)(gpu->sim.sgmsum);                // COM;
    com1[0]=E[6]/(double)(gpu->sim.sgmsum);                // COM;
    com1[1]=E[7]/(double)(gpu->sim.sgmsum);                // COM;
    com1[2]=E[8]/(double)(gpu->sim.sgmsum);                // COM;
    com2[0]=E[9]/(double)(gpu->sim.sgmsum);                // COM;
    com2[1]=E[10]/(double)(gpu->sim.sgmsum);                // COM;
    com2[2]=E[11]/(double)(gpu->sim.sgmsum);                // COM;
  
}

//---------------------------------------------------------------------------------------------
// gpu_download_vel_: download the velocities of all atoms from the device to the host.
//
// Arguments:
//   atm_vel:  the velocities of all atoms
//---------------------------------------------------------------------------------------------
extern "C" void gpu_download_sgld_(int sgatoms[],double avgx0[][3],double avgx1[][3],double avgx2[][3],double avgr[][3],double sgfps[],double sgpps[])
{
  PRINTMETHOD("gpu_download_sgld");
  if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
    if (gpu->bNewNeighborList) {
      gpu->pbImageIndex->Download();
      gpu->bNewNeighborList = false;
    }
    gpu->pbImageWsg->Download();
    gpu->pbImageX0sg->Download();
    gpu->pbImageX1sg->Download();
    gpu->pbImageX2sg->Download();
    gpu->pbImageRsg->Download();
    gpu->pbImageFPsg->Download();
    gpu->pbImagePPsg->Download();
    unsigned int* pImageAtomLookup = &(gpu->pbImageIndex->_pSysData[gpu->sim.imageStride * 2]);
        int *pWsg                                       = gpu->pbImageWsg->_pSysData;
        double *pX0sg                                    = gpu->pbImageX0sg->_pSysData;
        double *pX1sg                                    = gpu->pbImageX1sg->_pSysData;
        double *pX2sg                                    = gpu->pbImageX2sg->_pSysData;
        double *pRsg                                    = gpu->pbImageRsg->_pSysData;
        double *pFPsg                                    = gpu->pbImageFPsg->_pSysData;
        double *pPPsg                                    = gpu->pbImagePPsg->_pSysData;
    if (gpu->sim.pImageWsg != gpu->pbImageWsg->_pDevData) {
      pWsg = gpu->pbImageWsg->_pSysData + gpu->sim.stride;
    }
    if (gpu->sim.pImageX0sg != gpu->pbImageX0sg->_pDevData) {
      pX0sg = gpu->pbImageX0sg->_pSysData + gpu->sim.stride3;
    }
    if (gpu->sim.pImageX1sg != gpu->pbImageX1sg->_pDevData) {
      pX1sg = gpu->pbImageX1sg->_pSysData + gpu->sim.stride3;
    }
    if (gpu->sim.pImageX2sg != gpu->pbImageX2sg->_pDevData) {
      pX2sg = gpu->pbImageX2sg->_pSysData + gpu->sim.stride3;
    }
    if (gpu->sim.pImageRsgX != gpu->pbImageRsg->_pDevData) {
      pRsg = gpu->pbImageRsg->_pSysData + gpu->sim.stride3;
    }
    if (gpu->sim.pImageFPsg != gpu->pbImageFPsg->_pDevData) {
      pFPsg = gpu->pbImageFPsg->_pSysData + gpu->sim.stride;
    }
    if (gpu->sim.pImagePPsg != gpu->pbImagePPsg->_pDevData) {
      pPPsg = gpu->pbImagePPsg->_pSysData + gpu->sim.stride;
    }
    for (int i = 0; i < gpu->sim.atoms; i++) {
      int i1 = pImageAtomLookup[i];
      sgatoms[i] = pWsg[i1];
      avgx0[i][0] = pX0sg[i1];
      avgx0[i][1] = pX0sg[i1 + gpu->sim.stride];
      avgx0[i][2] = pX0sg[i1 + gpu->sim.stride2];
      avgx1[i][0] = pX1sg[i1];
      avgx1[i][1] = pX1sg[i1 + gpu->sim.stride];
      avgx1[i][2] = pX1sg[i1 + gpu->sim.stride2];
      avgx2[i][0] = pX2sg[i1];
      avgx2[i][1] = pX2sg[i1 + gpu->sim.stride];
      avgx2[i][2] = pX2sg[i1 + gpu->sim.stride2];
      avgr[i][0] = pRsg[i1];
      avgr[i][1] = pRsg[i1 + gpu->sim.stride];
      avgr[i][2] = pRsg[i1 + gpu->sim.stride2];
      sgfps[i] = pFPsg[i1];
      sgpps[i] = pPPsg[i1];
    }
  }
  else {
    gpu->pbWsg->Download();
    gpu->pbX0sg->Download();
    gpu->pbX1sg->Download();
    gpu->pbX2sg->Download();
    gpu->pbRsg->Download();
    gpu->pbFPsg->Download();
    gpu->pbPPsg->Download();
    for (int i = 0; i < gpu->sim.atoms; i++) {
      sgatoms[i] = gpu->pbWsg->_pSysData[i];
      avgx0[i][0] = gpu->pbX0sg->_pSysData[i];
      avgx0[i][1] = gpu->pbX0sg->_pSysData[i + gpu->sim.stride];
      avgx0[i][2] = gpu->pbX0sg->_pSysData[i + gpu->sim.stride2];
      avgx1[i][0] = gpu->pbX1sg->_pSysData[i];
      avgx1[i][1] = gpu->pbX1sg->_pSysData[i + gpu->sim.stride];
      avgx1[i][2] = gpu->pbX1sg->_pSysData[i + gpu->sim.stride2];
      avgx2[i][0] = gpu->pbX2sg->_pSysData[i];
      avgx2[i][1] = gpu->pbX2sg->_pSysData[i + gpu->sim.stride];
      avgx2[i][2] = gpu->pbX2sg->_pSysData[i + gpu->sim.stride2];
      avgr[i][0] = gpu->pbRsg->_pSysData[i];
      avgr[i][1] = gpu->pbRsg->_pSysData[i + gpu->sim.stride];
      avgr[i][2] = gpu->pbRsg->_pSysData[i + gpu->sim.stride2];
      sgfps[i] = gpu->pbFPsg->_pSysData[i];
      sgpps[i] = gpu->pbPPsg->_pSysData[i];
    }
  } 
}

//---------------------------------------------------------------------------------------------
// gpu_upload_last_vel_: upload the previous step's velocities.  This is called during system
//                       initialization.
//
// Arguments:
//   atm_last_vel:  the atomic velocities from a previous step
//---------------------------------------------------------------------------------------------
extern "C" void gpu_upload_sgld_(int sgatoms[],double avgx0[][3],double avgx1[][3],double avgx2[][3],double avgr[][3],double sgfps[],double sgpps[])
{
  PRINTMETHOD("gpu_upload_sgld");
  if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
    if (gpu->bNewNeighborList) {
      gpu->pbImageIndex->Download();
      gpu->bNewNeighborList = false;
    }
    unsigned int* pImageAtomLookup = &(gpu->pbImageIndex->_pSysData[gpu->sim.imageStride * 2]);
        int *pWsg                                       = gpu->pbImageWsg->_pSysData;
        double *pX0sg                                    = gpu->pbImageX0sg->_pSysData;
        double *pX1sg                                    = gpu->pbImageX1sg->_pSysData;
        double *pX2sg                                    = gpu->pbImageX2sg->_pSysData;
        double *pRsg                                    = gpu->pbImageRsg->_pSysData;
        double *pFPsg                                    = gpu->pbImageFPsg->_pSysData;
        double *pPPsg                                    = gpu->pbImagePPsg->_pSysData;
    if (gpu->sim.pImageWsg != gpu->pbImageWsg->_pDevData) {
      pWsg = gpu->pbImageWsg->_pSysData + gpu->sim.stride;
    }
    if (gpu->sim.pImageX0sg != gpu->pbImageX0sg->_pDevData) {
      pX0sg = gpu->pbImageX0sg->_pSysData + gpu->sim.stride3;
    }
    if (gpu->sim.pImageX1sg != gpu->pbImageX1sg->_pDevData) {
      pX1sg = gpu->pbImageX1sg->_pSysData + gpu->sim.stride3;
    }
    if (gpu->sim.pImageX2sg != gpu->pbImageX2sg->_pDevData) {
      pX2sg = gpu->pbImageX2sg->_pSysData + gpu->sim.stride3;
    }
    if (gpu->sim.pImageRsgX != gpu->pbImageRsg->_pDevData) {
      pRsg = gpu->pbImageRsg->_pSysData + gpu->sim.stride3;
    }
    if (gpu->sim.pImageFPsg != gpu->pbImageFPsg->_pDevData) {
      pFPsg = gpu->pbImageFPsg->_pSysData + gpu->sim.stride;
    }
    if (gpu->sim.pImagePPsg != gpu->pbImagePPsg->_pDevData) {
      pPPsg = gpu->pbImagePPsg->_pSysData + gpu->sim.stride;
    }
    for (int i = 0; i < gpu->sim.atoms; i++) {
      int i1 = pImageAtomLookup[i];
      pWsg[i1] = sgatoms[i];
      pX0sg[i1] = avgx0[i][0];
      pX0sg[i1 + gpu->sim.stride] = avgx0[i][1];
      pX0sg[i1 + gpu->sim.stride2] = avgx0[i][2];
      pX1sg[i1] = avgx1[i][0];
      pX1sg[i1 + gpu->sim.stride] = avgx1[i][1];
      pX1sg[i1 + gpu->sim.stride2] = avgx1[i][2];
      pX2sg[i1] = avgx2[i][0];
      pX2sg[i1 + gpu->sim.stride] = avgx2[i][1];
      pX2sg[i1 + gpu->sim.stride2] = avgx2[i][2];
      pRsg[i1] = avgr[i][0];
      pRsg[i1 + gpu->sim.stride] = avgr[i][1];
      pRsg[i1 + gpu->sim.stride2] = avgr[i][2];
      pFPsg[i1] = sgfps[i];
      pPPsg[i1] = sgpps[i];
    }
    gpu->pbImageWsg->Upload();
    gpu->pbImageX0sg->Upload();
    gpu->pbImageX1sg->Upload();
    gpu->pbImageX2sg->Upload();
    gpu->pbImageRsg->Upload();
    gpu->pbImageFPsg->Upload();
    gpu->pbImagePPsg->Upload();
  }
  else {
    for (int i = 0; i < gpu->sim.atoms; i++) {
      gpu->pbWsg->_pSysData[i] = sgatoms[i];
      gpu->pbX0sg->_pSysData[i] = avgx0[i][0];
      gpu->pbX0sg->_pSysData[i + gpu->sim.stride] = avgx0[i][1];
      gpu->pbX0sg->_pSysData[i + gpu->sim.stride2]= avgx0[i][2];
      gpu->pbX1sg->_pSysData[i] = avgx1[i][0];
      gpu->pbX1sg->_pSysData[i + gpu->sim.stride] = avgx1[i][1];
      gpu->pbX1sg->_pSysData[i + gpu->sim.stride2]= avgx1[i][2];
      gpu->pbX2sg->_pSysData[i] = avgx2[i][0];
      gpu->pbX2sg->_pSysData[i + gpu->sim.stride] = avgx2[i][1];
      gpu->pbX2sg->_pSysData[i + gpu->sim.stride2]= avgx2[i][2];
      gpu->pbRsg->_pSysData[i] = avgr[i][0];
      gpu->pbRsg->_pSysData[i + gpu->sim.stride] = avgr[i][1];
      gpu->pbRsg->_pSysData[i + gpu->sim.stride2]= avgr[i][2];
      gpu->pbFPsg->_pSysData[i] = sgfps[i];
      gpu->pbPPsg->_pSysData[i] = sgpps[i];
    }  
    gpu->pbWsg->Upload();
    gpu->pbX0sg->Upload();
    gpu->pbX1sg->Upload();
    gpu->pbX2sg->Upload();
    gpu->pbRsg->Upload();
    gpu->pbFPsg->Upload();
    gpu->pbPPsg->Upload();
  } 
}

//---------------------------------------------------------------------------------------------
// gpu_scale_sgld_: scale sgld arrays
//
// Arguments:
//   scale:   velocity scaling factor 
//   scalsg:   sgld scaling factor 
//---------------------------------------------------------------------------------------------
extern "C" void gpu_scale_sgld_(double* scale, double* scalsg,double* sgft, double* sgff, double* tsgset)
{
  PRINTMETHOD("gpu_scale_sgld");
  kScaleSGLD(gpu, *scale, *scalsg);

  gpu->sim.sgft                   = *sgft;
  gpu->sim.sgff                   = *sgff;
  gpu->sim.tsgset                 = *tsgset;
  gpuCopyConstants();

}

//---------------------------------------------------------------------------------------------
// gpu_shake_: execute SHAKE constraints on the GPU device
//---------------------------------------------------------------------------------------------
extern "C" void gpu_shake_()
{
  PRINTMETHOD("gpu_shake");
  kShake(gpu);
}

//---------------------------------------------------------------------------------------------
// gpu_rattle_: execute RATTLE constraints on the GPU device
//---------------------------------------------------------------------------------------------
extern "C" void gpu_rattle_(double* dt)
{
  PRINTMETHOD("gpu_rattle");
  kRattle(gpu, *dt);
}

//---------------------------------------------------------------------------------------------
// gpu_vrand_reset_velocities_: reset the velocities for initialization or application of the
//                              Andersen thermostat
//
// Arguments:
//   temp:      target temperature
//   half_dtx:  half of the scaled time step (half of dt * sqrt(418.4))
//---------------------------------------------------------------------------------------------
extern "C" void gpu_vrand_reset_velocities_(double* temp, double* half_dtx)
{
  PRINTMETHOD("gpu_vrand_reset_velocities");
  kResetVelocities(gpu, *temp, *half_dtx);
}

//---------------------------------------------------------------------------------------------
// gpu_recalculate_velocities_: recalculate velocities when applying geometry constraints
//
// Arguments:
//   dtx_inv:   inverse of the scaled time step (1.0 / [dt * sqt(418.4)])
//---------------------------------------------------------------------------------------------
extern "C" void gpu_recalculate_velocities_(double* dtx_inv)
{
  PRINTMETHOD("gpu_recalculate_velocities");
  if (gpu->ischeme == 0) {
      kRecalculateVelocities(gpu, *dtx_inv);
  }
  else {
      kMiddleRecalculateVelocities(gpu, *dtx_inv);
  }
}

//---------------------------------------------------------------------------------------------
// gpu_ti_calculate_kinetic_energy_:
//---------------------------------------------------------------------------------------------
extern "C" void gpu_ti_calculate_kinetic_energy_(double* c_ave, double* ti_eke_R1,
                                                 double* ti_eke_R2, double* ti_ekpbs_R1,
                                                 double* ti_ekpbs_R2, double* ti_ekph_R1,
                                                 double* ti_ekph_R2, double* ti_sc_eke_R1,
                                                 double* ti_sc_eke_R2)
{
  PRINTMETHOD("gpu_ti_calculate_kinetic_energy");
  kCalculateKineticEnergyAFE(gpu, *c_ave);
  cudaDeviceSynchronize();

  if (!(gpu->bCanMapHostMemory)) {
    gpu->pbAFEKineticEnergyBuffer->Download();
  }
  PMEFloat E[8];
  E[0] = 0.0f;     // ti_eke region 1
  E[1] = 0.0f;     // ti_eke region 2
  E[2] = 0.0f;     // ti_ekpbs region 1
  E[3] = 0.0f;     // ti_ekpbs region 2
  E[4] = 0.0f;     // ti_ekph region 1
  E[5] = 0.0f;     // ti_ekph region 2
  E[6] = 0.0f;     // ti_sc_eke region 1
  E[7] = 0.0f;     // ti_sc_eke region 2
  for (int i = 0; i < gpu->blocks; i++) {
    E[0] += gpu->pbAFEKineticEnergyBuffer->_pSysData[i].AFEKE.TI_EKER1;
    E[1] += gpu->pbAFEKineticEnergyBuffer->_pSysData[i].AFEKE.TI_EKER2;
    E[2] += gpu->pbAFEKineticEnergyBuffer->_pSysData[i].AFEKE.TI_EKPBSR1;
    E[3] += gpu->pbAFEKineticEnergyBuffer->_pSysData[i].AFEKE.TI_EKPBSR2;
    E[4] += gpu->pbAFEKineticEnergyBuffer->_pSysData[i].AFEKE.TI_EKPHR1;
    E[5] += gpu->pbAFEKineticEnergyBuffer->_pSysData[i].AFEKE.TI_EKPHR2;
    E[6] += gpu->pbAFEKineticEnergyBuffer->_pSysData[i].AFEKE.TI_SC_EKER1;
    E[7] += gpu->pbAFEKineticEnergyBuffer->_pSysData[i].AFEKE.TI_SC_EKER2;
  }
  *ti_eke_R1    = E[0];
  *ti_eke_R2    = E[1];
  *ti_ekpbs_R1  = E[2];
  *ti_ekpbs_R2  = E[3];
  *ti_ekph_R1   = E[4];
  *ti_ekph_R2   = E[5];
  *ti_sc_eke_R1 = E[6];
  *ti_sc_eke_R2 = E[7];
}

//---------------------------------------------------------------------------------------------
// gpu_scale_velocities_: scale velocities for Berendsen temperature coupling
//
// Arguments:
//   scale:   scaling factor to apply and make the system ever more like a flying ice cube
//---------------------------------------------------------------------------------------------
extern "C" void gpu_scale_velocities_(double* scale)
{
  PRINTMETHOD("gpu_scale_velocities");
  kScaleVelocities(gpu, *scale);
}

//---------------------------------------------------------------------------------------------
// gpu_calculate_kinetic_energy_: launch the kernel for calculating kinetic energy on the GPU.
//
// Arguments:
//   c_ave:    composite quantity referenced in Langevin dynamics
//   eke:      kinetic energy, seems to get calculated from last_vel
//   ekph:     additional kinetic energy energy accumulator calculated from vel
//   ekpbs:    additional kinetic energy energy accumulator calculated from vel
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_kinetic_energy_(double* c_ave, double* eke, double* ekph,
                                              double* ekpbs)
{
  PRINTMETHOD("gpu_calculate_kinetic_energy");
  kCalculateKineticEnergy(gpu, *c_ave);
  cudaDeviceSynchronize();

  if (!(gpu->bCanMapHostMemory)) {
    gpu->pbKineticEnergyBuffer->Download();
  }
  PMEFloat E[3];
  E[0] = 0.0f;   // EKE
  E[1] = 0.0f;   // EKPH
  E[2] = 0.0f;   // EKPBS
  for (int i = 0; i < gpu->blocks; i++) {
    E[0] += gpu->pbKineticEnergyBuffer->_pSysData[i].KE.EKE;
    E[1] += gpu->pbKineticEnergyBuffer->_pSysData[i].KE.EKPH;
    E[2] += gpu->pbKineticEnergyBuffer->_pSysData[i].KE.EKPBS;
  }
  *eke = E[0];
  *ekph = E[1];
  *ekpbs = E[2];
}

//---------------------------------------------------------------------------------------------
// gpu_recenter_molecule_: re-center a molecule.  This is used in Langevin thermostating.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_recenter_molecule_()
{
  PRINTMETHOD("gpu_recenter_molecule");
  kRecenter_Molecule(gpu);
}

//---------------------------------------------------------------------------------------------
// gpu_set_cpu_randoms_: set the choice of random number stream--GPU random numbers, or a C
//                       implementation of the Fortran CPU random number stream for GPU runs
//                       on all threads?
//
// Arguments:
//   buse_CPURandoms:  flag to make the choice above
//---------------------------------------------------------------------------------------------
extern "C" void gpu_set_cpu_randoms_(bool* buse_CPURandoms)
{
  gpu->bCPURandoms = *buse_CPURandoms;
}

//---------------------------------------------------------------------------------------------
// gpu_amrset_: set the random seed on the GPU
//
// Arguments:
//   seed:    the random seed to initialize the pseudo-random number sequence
//---------------------------------------------------------------------------------------------
extern "C" void gpu_amrset_(int* seed)
{
  PRINTMETHOD("gpu_amrset");
  if (gpu->bCPURandoms) {
    cpu_amrset(*seed);
    if (gpu->pbRandom->_pSysData == NULL) {
      delete gpu->pbRandom;
      gpu->pbRandom = new GpuBuffer<double>(gpu->sim.randomNumbers);
      gpu->sim.pRandom = gpu->pbRandom->_pDevData;
      gpu->sim.pRandomX = gpu->sim.pRandom;
      gpu->sim.pRandomY = gpu->sim.pRandomX +
                          gpu->sim.randomSteps * gpu->sim.paddedNumberOfAtoms;
      gpu->sim.pRandomZ = gpu->sim.pRandomY +
                          gpu->sim.randomSteps * gpu->sim.paddedNumberOfAtoms;
      gpuCopyConstants();
    }
    cpu_kRandom(gpu);
  }
  else {
    curandCreateGenerator(&gpu->RNG, CURAND_RNG_PSEUDO_DEFAULT);
    curandSetPseudoRandomGeneratorSeed(gpu->RNG, *seed);
    kRandom(gpu);
  }
  gpu->randomCounter = 0;
}

//---------------------------------------------------------------------------------------------
// gpu_dump_float_: dump an array of floating point descriptors, one for each atom / particle
//                  in the system.
//
// Arguments:
//   pFloat: the array of floating point numbers to dump (size must be the number of particles)
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_dump_float_(float* pFloat)
{
  FILE *fp = fopen("float.txt", "w");
  for (int i = 0; i < gpu->sim.atoms; i++) {
    fprintf(fp, "%32.15f\n", pFloat[i]);
    printf("%32.15f\n", pFloat[i]);
  }
  fclose(fp);
  exit(-1);
}

//---------------------------------------------------------------------------------------------
// gpu_dump_grid_weights_: dump the B-spline coefficients for all atoms in the system,
//                         numbered by their indices in the master topology.  This only works
//                         with 4th-order PME interpolation.
//
// Arguments:
//   atoms:     the number of atoms (particles) in the system
//   map:       the map of atom index in the B-spline list to atom index in the master topology
//   theta:     the B-spline coefficients X1, X2, X3, X4, Y1, Y2, ..., Z4
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_dump_grid_weights_(int* atoms, int* map, double* theta)
{
  for (int i = 1; i <= *atoms; i++) {
    for (int j = 0; j < *atoms; j++) {
      if (map[j] == i) {
        printf("%5d: %13.7f %13.7f %13.7f %13.7f\n", i - 1, theta[j * 12],
               theta[j * 12 + 1], theta[j * 12 + 2], theta[j * 12 + 3]);
        printf("%5d: %13.7f %13.7f %13.7f %13.7f\n", i - 1, theta[j * 12 + 4],
               theta[j * 12 + 5], theta[j * 12 + 6], theta[j * 12 + 7]);
        printf("%5d: %13.7f %13.7f %13.7f %13.7f\n", i - 1, theta[j * 12 + 8],
               theta[j * 12 + 9], theta[j * 12 + 10], theta[j * 12 + 11]);
      }
    }
  }
  exit(-1);
}

//---------------------------------------------------------------------------------------------
// gpu_dump_img_double_vector_: dump a vector of double-precision reals based on an image map.
//
// Arguments:
//   count:    the number of elements in the vector and map
//   img:      the map of which element to display in what order
//   pDouble:  the vector to display
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_dump_img_double_vector_(int* count, int* img, double* pDouble)
{
  for (int i = 1; i <= *count; i++) {
    for (int j = 0; j < *count; j++) {
      if (img[j] == i) {
        printf("%5d: %13.7f %13.7f %13.7f\n", i - 1, pDouble[j * 3],
               pDouble[j * 3 + 1], pDouble[j * 3 + 2]);
      }
    }
  }
  exit(-1);
}

//---------------------------------------------------------------------------------------------
// gpu_dump_img_int_vector_: dump a vector of integers based on an image map.
//
// Arguments:
//   count:    the number of elements in the vector and map
//   img:      the map of which element to display in what order
//   pInt:     the vector to display
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_dump_img_int_vector_(int* img, int* pInt)
{
  for (int i = 1; i <= gpu->sim.atoms; i++) {
    for (int j = 0; j < gpu->sim.atoms; j++) {
      if (img[j] == i) {
        printf("%5d: %9d %9d %9d\n", i - 1, pInt[j * 3], pInt[j * 3 + 1], pInt[j * 3 + 2]);
      }
    }
  }
}

//---------------------------------------------------------------------------------------------
// gpu_dump_complex_grid_: dump a grid of complex numbers, which is almost certainly the result
//                         of a forward FFT in PME.
//
// Arguments:
//   zxy_qt:           the grid of complex numbers in linear format (real imag real imag...)
//                     Probably also transposed in 3D.
//   fft_(x,y,z)_dim:  the dimensions of the grid
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_dump_complex_grid_(double zxy_qt[], int *fft_x_dim, int *fft_y_dim,
                                       int *fft_z_dim)
{
  for (int i = 0; i < *fft_x_dim; i++) {
    for (int j = 0; j < *fft_y_dim; j++) {
      for (int k = 0; k < *fft_z_dim; k++) {
        printf("%3d %3d %3d: %40.33f %40.33f\n", i, j, k,
               zxy_qt[j * (*fft_x_dim) * 2 * (*fft_z_dim) + i * 2 *(*fft_z_dim) + 2 * k],
               zxy_qt[j * (*fft_x_dim) * 2 * (*fft_z_dim) + i * 2 *(*fft_z_dim) + 2 * k + 1]);
      }
    }
  }
  exit(-1);
}

//---------------------------------------------------------------------------------------------
// gpu_load_complex_grid_: load a grid of complex numbers from disk, to substitute for an
//                         array in a calculation.  This is probably a way to intercept the
//                         PME reciprocal sum in flight and fix the outcome.
//
// Arguments:
//   zxy_qt:           the grid of complex numbers in linear format (real imag real imag...)
//                     Probably also transposed in 3D.
//   fft_(x,y,z)_dim:  the dimensions of the grid
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_load_complex_grid_(double zxy_qt[], int *fft_x_dim, int *fft_y_dim,
                                       int *fft_z_dim)
{
  gpu->pbXYZ_qt->Download();
  double maxerror = 0.0;
  for (int i = 0; i < *fft_x_dim; i++) {
    for (int j = 0; j < *fft_y_dim; j++) {
      for (int k = 0; k < *fft_z_dim; k++) {
        double error;
        error = fabs((zxy_qt[2*j*(*fft_x_dim)*(*fft_z_dim) + 2*i*(*fft_z_dim) + 2*k] -
                      gpu->pbXYZ_qt->_pSysData[(k*(*fft_y_dim) + j)*(*fft_x_dim) + i].x) /
                     gpu->pbXYZ_qt->_pSysData[(k*(*fft_y_dim) + j)*(*fft_x_dim) + i].x);
        if (error > maxerror) {
          printf("%3d %3d %3d: %40.30f %40.30f\n", i, j, k,
                 zxy_qt[j * (*fft_x_dim) * 2 * (*fft_z_dim) + i * 2 *(*fft_z_dim) + 2 * k],
                 gpu->pbXYZ_qt->_pSysData[(k * *fft_y_dim + j) * *fft_x_dim + i].x);
          maxerror = error;
        }
        error = fabs((zxy_qt[2*j*(*fft_x_dim)*(*fft_z_dim) + 2*i*(*fft_z_dim) + 2*k + 1] -
                      gpu->pbXYZ_qt->_pSysData[(k*(*fft_y_dim) + j)*(*fft_x_dim) + i].y) /
                     gpu->pbXYZ_qt->_pSysData[(k*(*fft_y_dim) + j)*(*fft_x_dim) + i].y);
        if (error > maxerror) {
          maxerror = error;
          printf("%3d %3d %3d: %40.30f %40.30f\n", i, j, k,
                 zxy_qt[j * (*fft_x_dim) * 2 * (*fft_z_dim) + i * 2 *(*fft_z_dim) + 2 * k + 1],
                 gpu->pbXYZ_qt->_pSysData[(k * *fft_y_dim + j) * *fft_x_dim + i].y);
          maxerror = error;
        }

        // This code was a bit messed up when I first saw it.  I think I've fixed it.
        gpu->pbXYZ_qt->_pSysData[(k * *fft_y_dim + j) * *fft_x_dim + i].x =
          zxy_qt[2*j*(*fft_x_dim)*(*fft_z_dim) + 2*i*(*fft_z_dim) + 2*k];
        gpu->pbXYZ_qt->_pSysData[(k * *fft_y_dim + j) * *fft_x_dim + i].y =
          zxy_qt[2*j*(*fft_x_dim)*(*fft_z_dim) + 2*i*(*fft_z_dim) + 2*k + 1];
      }
    }
  }
  printf("%40.20lf\n", maxerror);
  exit(-1);
}

//---------------------------------------------------------------------------------------------
// gpu_dump_grid_: dump a grid of real numbers.  Like gpu_dump_complex_grid_, but simpler to
//                 index.
//
// Arguments:
//   xyz_q:            the charge or electrostatic potential grid
//   fft_(x,y,z)_dim:  the dimensions of the grid
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_dump_grid_(double xyz_q[], int *fft_x_dim, int *fft_y_dim, int *fft_z_dim)
{
  double dmin = 99999.0;
  double dmax = -99999.0;
  int xstride = 2 * (*fft_x_dim / 2 + 1);
  for (int i = 0; i < *fft_x_dim; i++) {
    for (int j = 0; j < *fft_y_dim; j++) {
      for (int k = 0; k < *fft_z_dim; k++) {
        printf("%3d %3d %3d %32.15f\n", i, j, k, xyz_q[(k * *fft_y_dim + j) * xstride + i]);
        double v = xyz_q[(k * *fft_y_dim + j) * xstride + i];
        if (v > dmax) {
          dmax = v;
        }
        if (v < dmin) {
          dmin = v;
        }
      }
    }
  }
  printf("%f %f\n", dmin, dmax);
  exit(-1);
}

//---------------------------------------------------------------------------------------------
// gpu_dump_int_vector_: dump a vector of integers pertaining to every atom / particle.
//
// Arguments:
//   pInt:   the array of integers
//   atoms:  the number of particles
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_dump_int_vector_(int* pInt, int* atoms)
{
  FILE *fp = fopen("int.txt", "w");
  for (int i = 0; i < *atoms; i++) {
    fprintf(fp, "%5d: %9d %9d %9d\n", i, pInt[i*3], pInt[i*3+1], pInt[i*3+2]);
    printf("%5d: %9d %9d %9d\n", i, pInt[i*3], pInt[i*3+1], pInt[i*3+2]);
  }
  fclose(fp);
  exit(-1);
}

//---------------------------------------------------------------------------------------------
// gpu_dump_mapped_int_vector_: dump a mapped vector of integers pertaining to every atom /
//                              particle.
//
// Arguments:
//   pInt:   the array of integers
//   pMap:   the order in which to display each integer
//   atoms:  the number of particles (or just the number of integers to print)
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_dump_mapped_int_vector_(int* pInt, int* pMap, int* atoms)
{
  FILE *fp = fopen("int.txt", "w");
  for (int i = 1; i <= *atoms; i++) {
    for (int j = 0; j < *atoms; j++) {
      if (pMap[j] == i) {
        fprintf(fp, "%5d: %9d %9d %9d\n", i - 1, pInt[j*3], pInt[j*3+ 1], pInt[j*3+2]);
        printf("%5d: %9d %9d %9d\n", i - 1, pInt[j*3], pInt[j*3+1], pInt[j*3+2]);
      }
    }
  }
  fclose(fp);
  exit(-1);
}

//---------------------------------------------------------------------------------------------
// gpu_dump_mapped_double_vector_: dump a mapped vector of double precision reals pertaining
//                                 to every atom / particle.
//
// Arguments:
//   pDouble:   the array of reals
//   pMap:      the order in whcih to display each number
//   atoms:     the number of particles (or just the number of integers to print)
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_dump_mapped_double_vector_(double* pDouble, int* pMap, int* atoms)
{
  for (int i = 1; i <= *atoms; i++) {
    for (int j = 0; j < *atoms; j++) {
      if (pMap[j] == i) {
        printf("%5d: %20.15f %20.15f %20.15f\n", i - 1, pDouble[j*3], pDouble[j*3+1],
               pDouble[j*3+2]);
      }
    }
  }
  exit(-1);
}

//---------------------------------------------------------------------------------------------
// gpu_dump_int_: dump an integer.
//
// Arguments:
//   pInt:    the integer to dump
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_dump_int_(int* pInt)
{
  FILE *fp = fopen("int.txt", "w");
  for (int i = 0; i < gpu->sim.atoms; i++) {
    fprintf(fp, "%5d: %9d\n", i, pInt[i]);
    printf("%5d: %9d\n", i, pInt[i]);
  }
  fclose(fp);
  exit(-1);
}

//---------------------------------------------------------------------------------------------
// gpu_dump_double_: dump a double-precision real.
//
// Arguments:
//   pDouble:    the real number to dump
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_dump_double_(double* pDouble)
{
  FILE *fp = fopen("double.txt", "w");
  for (int i = 0; i < gpu->sim.atoms; i++) {
    fprintf(fp, "%32.15f\n", pDouble[i]);
    printf("%32.15f\n", pDouble[i]);
  }
  fclose(fp);
  exit(-1);
}

//---------------------------------------------------------------------------------------------
// gpu_dump_double_vector_: dump a double-precision real vector of triplets (i.e. coordinates,
//                          forces, velocities)
//
// Arguments:
//   atoms:    the number of real triplets to dump
//   pDouble:  the vector
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_dump_double_vector_(int* atoms, double* pDouble)
{
  for (int i = 0; i < *atoms; i++) {
    printf("%6d,%32.15f,%32.15f,%32.15f\n", i, pDouble[i*3], pDouble[i*3+1], pDouble[i*3+2]);
  }
}

//---------------------------------------------------------------------------------------------
// gpu_dump_dval_: how is this any different from gpu_dump_double_ above, except that the
//                 latter also writes to a file?
//
// Arguments:
//   pDouble:    the real number to dump
//
// This appears to be a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_dump_dval_(double* pDouble)
{
  printf("C %32.15f\n", *pDouble);
}

#ifdef MPI
//---------------------------------------------------------------------------------------------
// gpu_gather_gb_reff: gather the effective Born radii for all atoms
//---------------------------------------------------------------------------------------------
void gpu_gather_gb_reff()
{
  PRINTMETHOD("gpu_gather_gb_reff");
  gpu_allreduce(gpu->pbReffAccumulator, gpu->sim.stride);
  cudaDeviceSynchronize();
}

//---------------------------------------------------------------------------------------------
// gpu_gather_gb_temp7: gather parallel computations for GB calculations.
//---------------------------------------------------------------------------------------------
void gpu_gather_gb_temp7()
{
  PRINTMETHOD("gpu_gather_gb_temp7");
  gpu_allreduce(gpu->pbSumdeijdaAccumulator, gpu->sim.stride);
  cudaDeviceSynchronize();
}

//---------------------------------------------------------------------------------------------
// gpu_gather_gb_forces: gather forces on all particles computed across multiple GPUs.
//---------------------------------------------------------------------------------------------
void gpu_gather_gb_forces()
{
  PRINTMETHOD("gpu_gather_gb_forces");
  gpu_allreduce(gpu->pbForceAccumulator, gpu->sim.stride3);
  cudaDeviceSynchronize();
}
#endif

//---------------------------------------------------------------------------------------------
// gpu_gb_ene_: compute energies for Generalized Born simulations
//
// Arguments:
//   pEnergy:   record to store GB simulation energies
//   enmr:      NMR restraints energy (three components again)
//---------------------------------------------------------------------------------------------
extern "C" void gpu_gb_ene_(gb_pot_ene_rec* pEnergy, double enmr[3], int* ineb)
{
  PRINTMETHOD("gpu_gb_ene");
  kClearForces(gpu);
  kClearGBBuffers(gpu);
  if (gpu->ntf != 8) {
    if (gpu->sim.igb != 6) {
      kCalculateGBBornRadii(gpu);
#ifdef MPI
      gpu_gather_gb_reff();
#endif
      kReduceGBBornRadii(gpu);
    }
    kCalculateGBNonbondEnergy1(gpu);
    if (gpu->gbsa == 3) {
       // adding the slope and intercept maxSASA FOR energy only
#ifdef MPI
      if (gpu->gpuID == 0) {
#endif
        kReduceMaxsasaEsurf(gpu); //pwsasa only do when gbsa=3
#ifdef MPI
      }
#endif
    }
    if (gpu->sim.igb != 6) {
#ifdef MPI
      gpu_gather_gb_temp7();
      if (gpu->gpuID != 0) {
        kReduceGBTemp7(gpu);
      }
      else {
#endif
      kReduceGBTemp7Energy(gpu);
#ifdef MPI
      }
#endif
      kCalculateGBNonbondEnergy2(gpu);
    }
#ifdef MPI
    if (gpu->gpuID == 0) {
#endif
    kExecuteBondWorkUnits(gpu, 1);
    kCalculateNMREnergy(gpu);
#ifdef MPI
    }
#endif
  }
  // End contingency for ntf != 8

#ifdef MPI
  gpu_gather_gb_forces();
#endif
  gpu->pbEnergyBuffer->Download();
  PMEDouble energy[EXTENDED_ENERGY_TERMS];
  pEnergy->total = 0.0;
  for (int i = 0; i < VIRIAL_OFFSET; i++) {
    unsigned long long int val = gpu->pbEnergyBuffer->_pSysData[i];
    if (val >= 0x8000000000000000ull) {
      energy[i] = -(PMEDouble)(val ^ 0xffffffffffffffffull) / ENERGYSCALE;
    }
    else {
      energy[i] = (PMEDouble)val / ENERGYSCALE;
    }
    pEnergy->total += energy[i];
  }

#ifdef use_DPFP
  if (gpu->imin != 0)
  {
    for (int i = ENERGY_TERMS; i < EXTENDED_ENERGY_TERMS; i++)
    {
      unsigned long long int val = gpu->pbEnergyBuffer->_pSysData[i];
      energy[i] = (double)(*((long long int*)(&val)));
      pEnergy->total += energy[i];
    }
    energy[1] += energy[EVDWE_OFFSET];
    energy[0] += energy[EELTE_OFFSET];
  }
#endif

  // NEB all gather, followed by other result logging
#ifdef MPI
  if (*ineb > 0) {
      MPI_Allgather(&pEnergy->total, 1, MPI_PMEDOUBLE, gpu->pbNEBEnergyAll->_pSysData,
                    1, MPI_PMEDOUBLE, MPI_COMM_WORLD);
      gpu->pbNEBEnergyAll->Upload();
  }
#endif
  //printf("VVV %20.10f %20.10f %20.10f %20.10f\n", energy[0], energy[1], energy[EVDWE_OFFSET], energy[EELTE_OFFSET]);
  pEnergy->vdw_tot = energy[1];
  pEnergy->elec_tot = energy[0];
  pEnergy->gb = energy[2];
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
  if (gpu->gbsa == 3) {
    pEnergy->surf = energy[18]; // pwsasa and other sasa methods
  }
}

//---------------------------------------------------------------------------------------------
// gpu_clear_forces_: zero out all forces.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_clear_forces_()
{
  PRINTMETHOD("gpu_clear_forces");
  kClearForces(gpu);
}

//---------------------------------------------------------------------------------------------
// gpu_gb_forces_: compute forces for a Generalized Born simulation
//---------------------------------------------------------------------------------------------
extern "C" void gpu_gb_forces_()
{
  PRINTMETHOD("gpu_gb_forces");
  kClearForces(gpu);
  kClearGBBuffers(gpu);
  if (gpu->ntf != 8) {
    if (gpu->sim.igb != 6) {
      kCalculateGBBornRadii(gpu);
#ifdef MPI
      gpu_gather_gb_reff();
#endif
      kReduceGBBornRadii(gpu);
    }
    kCalculateGBNonbondForces1(gpu);
    if (gpu->sim.igb != 6) {
#ifdef MPI
      gpu_gather_gb_temp7();
#endif
      kReduceGBTemp7(gpu);
      kCalculateGBNonbondEnergy2(gpu);
    }
#ifdef MPI
    if (gpu->gpuID == 0) {
#endif
    kExecuteBondWorkUnits(gpu, 0);
    kCalculateNMRForces(gpu);
#ifdef MPI
    }
#endif
  }
#ifdef MPI
  gpu_gather_gb_forces();
#endif
}

//---------------------------------------------------------------------------------------------
// gpu_get_nb_energy_: get the non-bonded energy for a PME simulation.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_get_nb_energy_()
{
  PRINTMETHOD("gpu_get_nb_energy");
  kCalculatePMENonbondForces(gpu);
}

//---------------------------------------------------------------------------------------------
// gpu_gbsa3_setup_: setup for Generalized Born / Surface Area pwSASA calculations pwsasa
//---------------------------------------------------------------------------------------------
extern "C" void gpu_gbsa3_setup_(double sigma[], double epsilon[], double radius[],
                                 double maxsasa[])
{
  PRINTMETHOD("gpu_gbsa3_setup");

  int i;

  // Delete existing arrays
  delete gpu->pbGBSASigEpsiRadiMax;

  // Allocate new arrays
  gpu->pbGBSASigEpsiRadiMax        = new GpuBuffer<PMEFloat>(gpu->sim.stride4);
  gpu->sim.pgbsa_sigma             = gpu->pbGBSASigEpsiRadiMax->_pDevData;
  gpu->sim.pgbsa_epsilon           = gpu->sim.pgbsa_sigma + gpu->sim.stride;
  gpu->sim.pgbsa_radius            = gpu->sim.pgbsa_sigma + gpu->sim.stride2;
  gpu->sim.pgbsa_maxsasa           = gpu->sim.pgbsa_sigma + gpu->sim.stride3;

  // Fill arrays
  PMEFloat* pSigma                 = gpu->pbGBSASigEpsiRadiMax->_pSysData;
  PMEFloat* pEpsilon               = gpu->pbGBSASigEpsiRadiMax->_pSysData + gpu->sim.stride;
  PMEFloat* pRadius                = gpu->pbGBSASigEpsiRadiMax->_pSysData + gpu->sim.stride2;
  PMEFloat* pMaxsasa               = gpu->pbGBSASigEpsiRadiMax->_pSysData + gpu->sim.stride3;
  for (i = 0; i < gpu->sim.atoms; i++) {
    pSigma[i]   = sigma[i];
    pEpsilon[i] = epsilon[i];
    pRadius[i]  = radius[i];
    pMaxsasa[i] = maxsasa[i];
  }
  gpu->pbGBSASigEpsiRadiMax->Upload();
  gpuCopyConstants();
}


//---------------------------------------------------------------------------------------------
//

//---------------------------------------------------------------------------------------------
// gp
// gpu_gbsa_setup_: setup for Generalized Born / Surface Area calculations
//---------------------------------------------------------------------------------------------
extern "C" void gpu_gbsa_setup_()
{
  PRINTMETHOD("gpu_gbsa_setup");
  return;
}

//---------------------------------------------------------------------------------------------
// gpu_local_to_global_: commute local forces on extra points to the atoms with mass.  This is
//                       the GPU equivalent of the Fortran local_to_global subroutine.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_local_to_global_()
{
  PRINTMETHOD("gpu_local_to_global");
  kLocalToGlobal(gpu);
}

//---------------------------------------------------------------------------------------------
// gpu_amd_setup_: setup for accelerated Molecular Dynamics (aMD) on the GPU
//
// Arguments:
//   iamd:       flag to activate aMD
//   w_amd:      weight to use in aMD
//   iamdlag:    commence aMD after an equilibration period of this many steps
//   ntwx:       coordinate output frequency
//   EthreshP:
//   alphaP:     basic aMD parameters (see Donald Hamelberg's publications)
//   EthreshD:
//   alphaD:
//   EthreshP_w:
//   alphaP_w:   weighted versions of aMD parameters
//   EthreshD_w:
//   alphaD_w:
//   w_sign:     flag to switch from raising valleys to lowering barriers
//   temp0:      temperature at which to maintain the system
//---------------------------------------------------------------------------------------------
extern "C" void gpu_amd_setup_(int* iamd, int* w_amd, int* iamdlag, int* ntwx,
                               double* EthreshP, double* alphaP, double* EthreshD,
                               double* alphaD, double* EthreshP_w, double* alphaP_w,
                               double* EthreshD_w, double* alphaD_w, double* w_sign,
                               double* temp0)
{
  PRINTMETHOD("gpu_amd_setup");

  // Determine what type of AMD is being used
  gpu->sim.iamd    = *iamd;
  gpu->sim.w_amd   = *w_amd;
  gpu->sim.iamdlag = *iamdlag;

  // Allocate GPU data
  // Set up AMD parameters
  gpu->sim.amd_print_interval = *ntwx;
  gpu->sim.amd_EthreshP   = *EthreshP;
  gpu->sim.amd_alphaP     = *alphaP;
  gpu->sim.amd_EthreshD   = *EthreshD;
  gpu->sim.amd_alphaD     = *alphaD;
  gpu->sim.amd_EthreshP_w = *EthreshP_w;
  gpu->sim.amd_alphaP_w   = *alphaP_w;
  gpu->sim.amd_EthreshD_w = *EthreshD_w;
  gpu->sim.amd_alphaD_w   = *alphaD_w;
  gpu->sim.amd_temp0      = *temp0;
  gpu->sim.amd_w_sign     = *w_sign;
  gpu->pAmdWeightsAndEnergy = new PMEDouble[6];

  gpu->sim.AMDNumLag  = 0;
  gpu->sim.AMDNumRecs = 0;
  gpu->sim.AMDtboost  = 0.0;
  gpu->pbAMDfwgtd     = new GpuBuffer<PMEDouble>(1);
  gpu->sim.pAMDfwgtd  = gpu->pbAMDfwgtd->_pDevData;
  gpu->sim.AMDfwgt    = 1.0;
  gpuCopyConstants();

  return;
}

//---------------------------------------------------------------------------------------------
// gpu_download_amd_weights_: download all aMD weights (and energy)
//
// Arguments:
//   amd_weights_and_energy:  host-side array to store the critical values
//---------------------------------------------------------------------------------------------
extern "C" void gpu_download_amd_weights_(double amd_weights_and_energy[])
{
  PRINTMETHOD("gpu_download_amd_weights");
  amd_weights_and_energy[0] = gpu->pAmdWeightsAndEnergy[0];
  amd_weights_and_energy[1] = gpu->pAmdWeightsAndEnergy[1];
  amd_weights_and_energy[2] = gpu->pAmdWeightsAndEnergy[2];
  amd_weights_and_energy[3] = gpu->pAmdWeightsAndEnergy[3];
  amd_weights_and_energy[4] = gpu->pAmdWeightsAndEnergy[4];
  amd_weights_and_energy[5] = gpu->pAmdWeightsAndEnergy[5];
}

#ifdef MPI
//---------------------------------------------------------------------------------------------
// gpu_calculate_amd_dihedral_weight_: parallel version of the function to compute AMD
//                                     dihedral weight.  This applies to PME /
//                                     periodic simulations.
//
// Arguments:
//   totdih:   total dihedral potential
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_amd_dihedral_weight_(double* totdih)
{
  PRINTMETHOD("gpu_calculate_amd_dihedral_weight");

  // calculate AMD weight for dihedral boost (E_dih_boost)
  PMEDouble EthreshD    = gpu->sim.amd_EthreshD;
  PMEDouble alphaD      = gpu->sim.amd_alphaD;
  PMEDouble EthreshD_w  = gpu->sim.amd_EthreshD_w;
  PMEDouble alphaD_w    = gpu->sim.amd_alphaD_w;
  PMEDouble w_sign      = gpu->sim.amd_w_sign;
  PMEDouble E_dih_boost = 0.0;
  PMEDouble fwgtd       = 1.0;
  PMEDouble EV          = 0.0;
  PMEDouble windowed_factor, windowed_derivative;
  PMEDouble windowed_exp;

  EV = (EthreshD - (*totdih)) * w_sign;
  if (EV > 0.0) {
    if (gpu->sim.AMDNumLag == 0) {

      // Dih boost E in kcal/mol
      E_dih_boost = (EV * EV) / (alphaD + EV);
      fwgtd = (alphaD * alphaD) / ((alphaD + EV) * (alphaD + EV));
      if (gpu->sim.w_amd != 0 ) {
        windowed_exp = exp((EthreshD_w - *totdih) * w_sign/alphaD_w);
        windowed_factor = 1.0/(1.0 + windowed_exp);
        windowed_derivative = E_dih_boost * windowed_factor * windowed_factor * windowed_exp /
                              alphaD_w;
        E_dih_boost = E_dih_boost * windowed_factor;
        fwgtd = fwgtd * windowed_factor + windowed_derivative;
      }
    }
  }
  gpu->sim.AMDtboost = E_dih_boost;
  gpu->pbAMDfwgtd->_pSysData[0] = fwgtd;
  gpu->pbAMDfwgtd->Upload();
}

//---------------------------------------------------------------------------------------------
// gpu_calculate_amd_dihedral_energy_: parallel version of the function to compute dihedral
//                                     energy boost
//
// Arguments:
//   totdih:   total dihedral potential
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_amd_dihedral_energy_(double* totdih)
{
  PRINTMETHOD("gpu_calculate_amd_dihedral_energy");

  // Rebuild neighbor list
  gpu_build_neighbor_list_();

  // Local energy
  if (gpu->bCalculateLocalForces) {

    // Calculate AMD energy for dihedral boost (tboost)
    kExecuteBondWorkUnits(gpu, 2);
    gpu->pbEnergyBuffer->Download();
    PMEDouble totdih2 = (PMEDouble)gpu->pbEnergyBuffer->_pSysData[AMD_E_DIHEDRAL_OFFSET] /
                        ENERGYSCALE;
    *totdih = totdih2;
  }
  else {
    *totdih = 0.0;
  }
}
#else
//---------------------------------------------------------------------------------------------
// gpu_calculate_amd_dihedral_energy_weight_: combined function to compute both the energy and
//                                            weight of dihedral interactions in serial GPU
//                                            mode
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_amd_dihedral_energy_weight_()
{
  PRINTMETHOD("gpu_calculate_amd_dihedral_energy_weight");
  PMEDouble totdih        = 0.0;
  if (gpu->sim.AMDNumLag == 0) {

    // Rebuild neighbor list
    gpu_build_neighbor_list_();

    // Calculate AMD calculate energy for dihedral boost (tboost)
    kExecuteBondWorkUnits(gpu, 2);

    //kCalculateAMDDihedralEnergy(gpu);
    gpu->pbEnergyBuffer->Download();
    totdih = (PMEDouble)gpu->pbEnergyBuffer->_pSysData[AMD_E_DIHEDRAL_OFFSET] / ENERGYSCALE;
  }

  // Calculate AMD weight for dihedral boost (tboost)
  PMEDouble EthreshD    = gpu->sim.amd_EthreshD;
  PMEDouble alphaD      = gpu->sim.amd_alphaD;
  PMEDouble EthreshD_w  = gpu->sim.amd_EthreshD_w;
  PMEDouble alphaD_w    = gpu->sim.amd_alphaD_w;
  PMEDouble w_sign      = gpu->sim.amd_w_sign;
  PMEDouble E_dih_boost = 0.0;
  PMEDouble fwgtd       = 1.0;
  PMEDouble EV          = 0.0;
  PMEDouble windowed_factor, windowed_derivative;
  PMEDouble windowed_exp;

  EV = (EthreshD - (totdih)) * w_sign;
  if (EV > 0.0) {
    if (gpu->sim.AMDNumLag == 0) {
      E_dih_boost = (EV * EV) / (alphaD + EV);
      fwgtd = (alphaD * alphaD) / ((alphaD + EV) * (alphaD + EV));
      if (gpu->sim.w_amd != 0) {
        windowed_exp = exp((EthreshD_w - totdih) * w_sign/alphaD_w);
        windowed_factor = 1.0/(1.0 + windowed_exp);
        windowed_derivative = E_dih_boost * windowed_factor * windowed_factor * windowed_exp /
                              alphaD_w;
        E_dih_boost = E_dih_boost * windowed_factor;
        fwgtd = fwgtd * windowed_factor + windowed_derivative;
      }
    }
  }
  gpu->sim.AMDtboost = E_dih_boost;
  gpu->pbAMDfwgtd->_pSysData[0] = fwgtd;
  gpu->pbAMDfwgtd->Upload();
}
#endif

//---------------------------------------------------------------------------------------------
// gpu_calculate_and_apply_amd_weights_: calculate and apply the aMD weights.  Boosts for total
//                                       potential energy and dihedral energy are considered.
//
// Arguments:
//   pot_ene_tot:   total potential energy
//   dih_ene_tot:   total diehdral energy
//   amd_ene_tot:   total boost energy
//   num_amd_lag:   counter to track the number of steps through the lag phase, and eventually
//                  activate aMD when it gets large enough.  Mirrors sim.AMDNumLag on the GPU.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_and_apply_amd_weights_(double* pot_ene_tot, double* dih_ene_tot,
                                                     double* amd_ene_tot, double* num_amd_lag)
{
  PRINTMETHOD("gpu_calculate_and_apply_amd_weights");
  double E_dih_boost   = 0.0;
  double fwgtd         = 1.0;
  double fwgt          = 1.0;
  double E_total_boost = 0.0;
  double temp0         = gpu->sim.amd_temp0;
  double ONE_KB        = 1.0 / (temp0 * KB);
  PMEDouble windowed_factor, windowed_derivative;
  PMEDouble windowed_exp;

  if (gpu->sim.AMDNumLag == 0) {

    //calculate AMD weight, seting dihedral boost (E_dih_boost) to zero for now
    double EthreshP      = gpu->sim.amd_EthreshP;
    double alphaP        = gpu->sim.amd_alphaP;
    double EV            = 0.0;
    PMEDouble EthreshP_w = gpu->sim.amd_EthreshD_w;
    PMEDouble alphaP_w   = gpu->sim.amd_alphaD_w;
    PMEDouble w_sign     = gpu->sim.amd_w_sign;
    if ((gpu->sim.iamd == 2)||(gpu->sim.iamd == 3)) {
      E_dih_boost = gpu->sim.AMDtboost;
      fwgtd       = gpu->pbAMDfwgtd->_pSysData[0];
    }
    double totalenergy  = *pot_ene_tot + E_dih_boost;

    EV = (EthreshP - totalenergy) * w_sign;
    if (((gpu->sim.iamd == 1) || (gpu->sim.iamd == 3)) && (EV > 0.0)) {

      // PE boost in Kcal/mol
      E_total_boost = (EV*EV) / (alphaP + EV);
      fwgt = (alphaP * alphaP) / ((alphaP + EV) * (alphaP + EV));
      if (gpu->sim.w_amd != 0) {
        windowed_exp = exp((EthreshP_w - totalenergy) * w_sign/alphaP_w);
        windowed_factor = 1.0 / (1.0 + windowed_exp);
        windowed_derivative =  E_total_boost * windowed_factor * windowed_factor *
                               windowed_exp / alphaP_w;
        E_total_boost = E_total_boost * windowed_factor;
        fwgt = fwgt*windowed_factor + windowed_derivative;
      }
    }

    // Calculate AMD weight
    kCalculateAMDWeightAndScaleForces(gpu, *pot_ene_tot, *dih_ene_tot, fwgt);
  }

  int numrecs = gpu->sim.AMDNumRecs * 6;
  if ((gpu->sim.AMDNumRecs+1)  >= gpu->sim.amd_print_interval) {
    gpu->pAmdWeightsAndEnergy[0] = *pot_ene_tot;
    gpu->pAmdWeightsAndEnergy[1] = *dih_ene_tot;
    gpu->pAmdWeightsAndEnergy[2] = fwgt;
    gpu->pAmdWeightsAndEnergy[3] = fwgtd;
    gpu->pAmdWeightsAndEnergy[4] = E_total_boost;
    gpu->pAmdWeightsAndEnergy[5] = E_dih_boost;
  }

  if (gpu->sim.AMDNumLag == gpu->sim.iamdlag) {
    gpu->sim.AMDNumLag  = 0;
  }
  else {
    gpu->sim.AMDNumLag++;
  }

  gpu->sim.AMDNumRecs++;
  if (gpu->sim.AMDNumRecs  >= gpu->sim.amd_print_interval) {
    gpu->sim.AMDNumRecs = 0 ;
  }
  *num_amd_lag = gpu->sim.AMDNumLag;
  *amd_ene_tot = E_dih_boost + E_total_boost;
}

//---------------------------------------------------------------------------------------------
// AMD functions for gb, no need to do neighbor lists here
//---------------------------------------------------------------------------------------------
#ifdef MPI
//---------------------------------------------------------------------------------------------
// gpu_calculate_gb_amd_dihedral_energy_: GB equivalent of the function above for periodic
//                                        simulations.
//
// Arguments:
//   totdih:   total dihedral potential
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_gb_amd_dihedral_energy_(double* totdih)
{
  PRINTMETHOD("gpu_calculate_gb_amd_dihedral_energy");

  // Local energy
  if (gpu->bCalculateLocalForces) {

    // Calculate AMD energy for dihedral boost (tboost)
    kExecuteBondWorkUnits(gpu, 2);
    gpu->pbEnergyBuffer->Download();
    PMEDouble totdih2 = (PMEDouble)gpu->pbEnergyBuffer->_pSysData[AMD_E_DIHEDRAL_OFFSET] /
                        ENERGYSCALE;
    *totdih = totdih2;
  }
  else {
    *totdih             = 0.0;
  }
}

#else
//---------------------------------------------------------------------------------------------
// gpu_calculate_gb_amd_dihedral_energy_weight_: the name says it all.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_gb_amd_dihedral_energy_weight_()
{
  PRINTMETHOD("gpu_calculate_gb_amd_dihedral_energy_weight");
  PMEDouble totdih = 0.0;
  if (gpu->sim.AMDNumLag == 0) {

    // Calculate AMD energy for dihedral boost (tboost)
    kExecuteBondWorkUnits(gpu, 2);
    gpu->pbEnergyBuffer->Download();
    totdih = (PMEDouble)gpu->pbEnergyBuffer->_pSysData[AMD_E_DIHEDRAL_OFFSET] / ENERGYSCALE;
  }

  // Calculate AMD weight for dihedral boost (E_dih_boost)
  PMEDouble EthreshD      = gpu->sim.amd_EthreshD;
  PMEDouble alphaD        = gpu->sim.amd_alphaD;
  PMEDouble EthreshD_w    = gpu->sim.amd_EthreshD_w;
  PMEDouble alphaD_w      = gpu->sim.amd_alphaD_w;
  PMEDouble w_sign        = gpu->sim.amd_w_sign;
  PMEDouble E_dih_boost   = 0.0;
  PMEDouble fwgtd         = 1.0;
  PMEDouble EV            = 0.0;
  PMEDouble windowed_factor, windowed_derivative;
  PMEDouble windowed_exp;

  EV = (EthreshD - (totdih)) * w_sign;
  if (EV > 0.0) {
    if (gpu->sim.AMDNumLag == 0) {

      // Dih boost E in kcal/mol
      E_dih_boost = (EV * EV) / (alphaD + EV);
      fwgtd = (alphaD * alphaD)/((alphaD + EV) * (alphaD + EV));
      if (gpu->sim.w_amd != 0) {
        windowed_exp = exp((EthreshD_w - totdih) * w_sign/alphaD_w);
        windowed_factor = 1.0/(1.0 + windowed_exp);
        windowed_derivative = E_dih_boost * windowed_factor * windowed_factor * windowed_exp /
                              alphaD_w;
        E_dih_boost = E_dih_boost * windowed_factor;
        fwgtd = fwgtd * windowed_factor + windowed_derivative;
      }
    }
  }
  gpu->sim.AMDtboost = E_dih_boost;
  gpu->pbAMDfwgtd->_pSysData[0] = fwgtd;
  gpu->pbAMDfwgtd->Upload();
}
#endif

//---------------------------------------------------------------------------------------------
// gpu_gamd_setup_: setup for Gaussian Accelerated Molecular Dynamics simulations
//
// Arguments:
//   igamd:      flag to activate GaMD
//   igamdlag:   lag time for commencement of GaMD boosting (analogous to iamdlag in the
//               original aMD)
//   ntwx:       coordinate output frequency
//   EthreshP:   threshold total potential energy for adding total potential boost
//   kP:         harmonic constant for adding total potential boost
//   EtreshD:    threshold dihedral potential energy for adding dihedral boost
//   kD:         harmonic constant for adding dihedral boost
//   temp0:      target temperature for the simulation
//---------------------------------------------------------------------------------------------
extern "C" void gpu_gamd_setup_(int* igamd, int* igamdlag, int* ntwx, double* EthreshP,
                                double* kP, double* EthreshD, double* kD, double* temp0)
{
  PRINTMETHOD("gpu_gamd_setup");

  // Determine what type of GaMD is being used
  gpu->sim.igamd    = *igamd;
  gpu->sim.igamd0   = *igamd;
  gpu->sim.igamdlag = *igamdlag;

  // Allocate GPU data
  // Set up GaMD parameters
  gpu->sim.gamd_print_interval = *ntwx;
  gpu->sim.gamd_EthreshP       = *EthreshP;
  gpu->sim.gamd_kP             = *kP;
  gpu->sim.gamd_EthreshD       = *EthreshD;
  gpu->sim.gamd_kD             = *kD;
  gpu->sim.gamd_temp0          = *temp0;
  gpu->pGaMDWeightsAndEnergy   = new PMEDouble[6];
  gpu->sim.GaMDNumLag          = 0;
  gpu->sim.GaMDNumRecs         = 0;
  gpu->sim.GaMDtboost          = 0.0;
  gpu->pbGaMDfwgtd             = new GpuBuffer<PMEDouble>(1);
  gpu->sim.pGaMDfwgtd          = gpu->pbGaMDfwgtd->_pDevData;
  gpu->pbGaMDfwgtd->_pSysData[0] = 1.0;
  gpu->pbGaMDfwgtd->Upload();
  gpu->sim.GaMDfwgt            = 1.0;
  gpuCopyConstants();

  return;
}
//---------------------------------------------------------------------------------------------
// gpu_gamd_setup_: setup for Gaussian Accelerated Molecular Dynamics simulations
//
// Arguments:
//   igamd:      flag to activate GaMD
//   igamdlag:   lag time for commencement of GaMD boosting (analogous to iamdlag in the
//               original aMD)
//   ntwx:       coordinate output frequency
//   EthreshP:   threshold total potential energy for adding total potential boost
//   kP:         harmonic constant for adding total potential boost
//   EtreshD:    threshold dihedral potential energy for adding
//   kD:         harmonic constant for adding dihedral boost
//   temp0:      target temperature for the simulation
//---------------------------------------------------------------------------------------------
//
extern "C" void gpu_gamd_ppi_setup_(int* igamd, int* igamdlag, int* ntwx, double* EthreshP,
		                         double* kP, double* EthreshD, double* kD,double* EthreshB, double* k_B, double* temp0)
{
 PRINTMETHOD("gpu_gamd_setup");

  // Determine what type of GaMD is being use
  gpu->sim.igamd    = *igamd;
  gpu->sim.igamd0   = *igamd;
  gpu->sim.igamdlag = *igamdlag;
  // Allocate GPU data
  //   // Set up GaMD parameters
  gpu->sim.gamd_print_interval = *ntwx;
  gpu->sim.gamd_EthreshP       = *EthreshP;
  gpu->sim.gamd_kP             = *kP;
  gpu->sim.gamd_EthreshB       = *EthreshB;
  gpu->sim.gamd_kB             = *k_B;
  gpu->sim.gamd_EthreshD       = *EthreshD;
  gpu->sim.gamd_kD             = *kD;
  gpu->sim.gamd_temp0          = *temp0;
  gpu->pGaMDWeightsAndEnergy   = new PMEDouble[6];
  gpu->sim.GaMDNumLag          = 0;
  gpu->sim.GaMDNumRecs         = 0;
  gpu->sim.GaMDtboost          = 0.0;
  gpu->pbGaMDfwgtd             = new GpuBuffer<PMEDouble>(1);
  gpu->sim.pGaMDfwgtd          = gpu->pbGaMDfwgtd->_pDevData;
  gpu->pbGaMDfwgtd->_pSysData[0] = 1.0;
  gpu->pbGaMDfwgtd->Upload();
  gpu->sim.GaMDfwgt            = 1.0;
  gpu->sim.GaMDfwgtb            = 1.0;
  gpuCopyConstants();

  return;
}


//---------------------------------------------------------------------------------------------
// gpu_igamd_update_: update GaMD flag
//
// Arguments:
//   igamd: flag of GaMD
//---------------------------------------------------------------------------------------------
extern "C" void gpu_igamd_update_(int* igamd)
{
  PRINTMETHOD("gpu_igamd_update");
  gpu->sim.igamd    = *igamd;
  return;
}

//---------------------------------------------------------------------------------------------
// gpu_gamd_update_: update critical parameters for GaMD
//
// Arguments:
//   EthreshP:
//   kP:        GaMD parameters, with analogs in standard aMD.  kD and kP are boost factors
//   EthreshD:  for dihedral and total potential energy, respectively.
//   kD:
//---------------------------------------------------------------------------------------------
extern "C" void gpu_gamd_update_(double* EthreshP, double* kP, double* EthreshD, double* kD)
{
  PRINTMETHOD("gpu_gamd_update");
  gpu->sim.gamd_EthreshP = *EthreshP;
  gpu->sim.gamd_kP       = *kP;
  gpu->sim.gamd_EthreshD = *EthreshD;
  gpu->sim.gamd_kD       = *kD;
  gpuCopyConstants();

  return;
}

//---------------------------------------------------------------------------------------------
// gpu_gamd_triple_update_: update critical parameters for GaMD
//
// Arguments:
//   EthreshP:
//   kP:        GaMD parameters, with analogs in standard aMD.  kD and kP K_B are boost factors
//   EthreshD:  for first,second and third potential energy, respectively.
//   kD:
//   EthreshB:
//   k_B:
//
//--------------------------------------------------------------------------------------------
extern "C" void gpu_gamd_triple_update_(double* EthreshP, double* kP, double* EthreshD, double* kD,double* EthreshB, double* k_B)
{
  PRINTMETHOD("gpu_gamd_update");
  gpu->sim.gamd_EthreshP = *EthreshP;
  gpu->sim.gamd_kP       = *kP;
  gpu->sim.gamd_EthreshD = *EthreshD;
  gpu->sim.gamd_kD       = *kD;
  gpu->sim.gamd_EthreshB=*EthreshB;
  gpu->sim.gamd_kB=*k_B;

  gpuCopyConstants();
  return;
}
// gpu_download_gamd_weights_: download critical parameters and results from a GaMD simulation
//
// Arguments:
//   gamd_weights_and_energy:   obvious?
//---------------------------------------------------------------------------------------------
extern "C" void gpu_download_gamd_weights_(double gamd_weights_and_energy[])
{
  PRINTMETHOD("gpu_download_gamd_weights");
  gamd_weights_and_energy[0] = gpu->pGaMDWeightsAndEnergy[0];
  gamd_weights_and_energy[1] = gpu->pGaMDWeightsAndEnergy[1];
  gamd_weights_and_energy[2] = gpu->pGaMDWeightsAndEnergy[2];
  gamd_weights_and_energy[3] = gpu->pGaMDWeightsAndEnergy[3];
  gamd_weights_and_energy[4] = gpu->pGaMDWeightsAndEnergy[4];
  gamd_weights_and_energy[5] = gpu->pGaMDWeightsAndEnergy[5];
}

#ifdef MPI
//---------------------------------------------------------------------------------------------
// gpu_calculate_gamd_dihedral_weight_: parallel function for periodic / PME simulations with
//                                      GaMD.  Calculates the weight of the dihedral boost.
//
// Arguments:
//   totdih:   the total dihedral potential (this is a local variable in the serial GPU
//             variant of the function)
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_gamd_dihedral_weight_(double* totdih)
{
  PRINTMETHOD("gpu_calculate_gamd_dihedral_weight");

  // calculate GaMD weight for dihedral boost (tboost)
  PMEDouble EthreshD      = gpu->sim.gamd_EthreshD;
  PMEDouble kD            = gpu->sim.gamd_kD;
  PMEDouble tboost        = 0.0;
  PMEDouble fwgtd         = 1.0;
  PMEDouble EV            = 0.0;

  EV                      = (EthreshD - (*totdih));
  if (EV > 0.0) {
    if (gpu->sim.GaMDNumLag == 0) {
      tboost = 0.5 * kD * (EV * EV);
      fwgtd = 1.0 - kD * EV ;
    }
  }
  gpu->sim.GaMDtboost = tboost;
  gpu->pbGaMDfwgtd->_pSysData[0] = fwgtd;
  gpu->pbGaMDfwgtd->Upload();
}

//---------------------------------------------------------------------------------------------
// gpu_calculate_gamd_dihedral_energy_: the name says it all, a function that goes with
//                                      gpu_calculate_gamd_dihedral_weight_ above to carry out
//                                      GaMD with dihedral boosting.
//
// Arguments:
//   totdih:   the total dihedral potential (this is a local variable in the serial GPU
//             variant of the function)
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_gamd_dihedral_energy_(double* totdih)
{
  PRINTMETHOD("gpu_calculate_gamd_dihedral_energy");

  // Rebuild neighbor list
  gpu_build_neighbor_list_();

  // Local energy
  if (gpu->bCalculateLocalForces) {

    // Calculate GaMD energy for dihedral boost (tboost)
    kExecuteBondWorkUnits(gpu, 3);
    gpu->pbEnergyBuffer->Download();
    PMEDouble totdih2 = (PMEDouble)gpu->pbEnergyBuffer->_pSysData[GAMD_E_DIHEDRAL_OFFSET] /
                        ENERGYSCALE;
    *totdih = totdih2;
  }
  else {
    *totdih             = 0.0;
  }
}

#else
//---------------------------------------------------------------------------------------------
// gpu_calculate_gamd_dihedral_energy_weight_: as with standard GaMD, the serial GPU version
//                                             merges two functions in the parallel code.  What
//                                             was an argument to those functions (totdih) is a
//                                             local variable in this one.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_gamd_dihedral_energy_weight_(){
  PRINTMETHOD("gpu_calculate_gamd_dihedral_energy_weight");

  PMEDouble totdih = 0.0;
  if (gpu->sim.GaMDNumLag == 0) {

    // Rebuild neighbor list
    gpu_build_neighbor_list_();

    // calculate energy for dihedral boost (tboost)
    kExecuteBondWorkUnits(gpu, 3);
    gpu->pbEnergyBuffer->Download();
    totdih = (PMEDouble)gpu->pbEnergyBuffer->_pSysData[GAMD_E_DIHEDRAL_OFFSET] / ENERGYSCALE;
  }

  // Calculate GAMD weight for dihedral boost (tboost)
  PMEDouble EthreshD = gpu->sim.gamd_EthreshD;
  PMEDouble kD       = gpu->sim.gamd_kD;
  PMEDouble tboost   = 0.0;
  PMEDouble fwgtd    = 1.0;
  PMEDouble EV       = 0.0;
  EV = EthreshD - (totdih);

  // Calculate GaMD weight for dihedral boost (tboost)
  if (EV > 0.0) {
    if (gpu->sim.GaMDNumLag == 0) {
      tboost = 0.5 * kD * (EV * EV);
      fwgtd = 1.0 - kD * EV ;
    }
  }
  gpu->sim.GaMDtboost = tboost;
  gpu->pbGaMDfwgtd->_pSysData[0] = fwgtd;
  gpu->pbGaMDfwgtd->Upload();
}

//---------------------------------------------------------------------------------------------
// gpu_calculate_gamd_dihedral_energy_weight_ti_: as with standard GaMD, the serial GPU version
//                                             merges two functions in the parallel code.  What
//                                             was an argument to those functions (totdih) is a
//                                             local variable in this one.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_gamd_dihedral_energy_weight_ti_(double* totdih, double* dih_boost){
  PRINTMETHOD("gpu_calculate_gamd_dihedral_energy_weight_ti");

  // Calculate GAMD weight for dihedral boost (tboost)
  PMEDouble EthreshD = gpu->sim.gamd_EthreshD;
  PMEDouble kD       = gpu->sim.gamd_kD;
  PMEDouble tboost   = 0.0;
  PMEDouble fwgtd    = 1.0;
  PMEDouble EV       = 0.0;
  EV = EthreshD - (*totdih);

  // Calculate GaMD weight for dihedral boost (tboost)
  if (EV > 0.0) {
    if (gpu->sim.GaMDNumLag == 0) {
      tboost = 0.5 * kD * (EV * EV);
      fwgtd = 1.0 - kD * EV ;
    }
  }
  gpu->sim.GaMDtboost = tboost;
  gpu->pbGaMDfwgtd->_pSysData[0] = fwgtd;
  gpu->pbGaMDfwgtd->Upload();
  *dih_boost = tboost;
  // printf("Debug-GaMD-GPU) (EthreshD, fwgtd, totdih, tboost, totdih+tboost) = (%12.5f, %12.5f, %12.5f, %12.5f, %12.5f)\n", EthreshD, fwgtd, *totdih, tboost, (*totdih+tboost));

  if (gpu->sim.GaMDNumLag == 0) {

    // Rebuild neighbor list
    gpu_build_neighbor_list_();

    // calculate energy for dihedral boost (tboost)
    kExecuteBondWorkUnits(gpu, 3);
    gpu->pbEnergyBuffer->Download();
    PMEDouble totdih_new = (PMEDouble)gpu->pbEnergyBuffer->_pSysData[GAMD_E_DIHEDRAL_OFFSET] / ENERGYSCALE;
    // printf("Debug-GaMD-GPU) totdih_new = (%12.5f)\n", totdih_new);
  }

}
#endif

//---------------------------------------------------------------------------------------------
// gpu_calculate_and_apply_gamd_weights_: analogous to the function in standard aMD for
//                                        carrying out the weighting and boosting protocol.
//
// Arguments:
//   pot_ene_tot:    total potential energy
//   dih_ene_tot:    total dihedral energy
//   gamd_ene_tot:   total GaMD boost potential
//   num_gamd_lag:   counter to track when the gamdlag phase ends
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_and_apply_gamd_weights_(double* pot_ene_tot, double* dih_ene_tot,
                                                      double* gamd_ene_tot,
                                                      double* num_gamd_lag)
{
  PRINTMETHOD("gpu_calculate_and_apply_gamd_weights");
  double tboost    = 0.0;
  double fwgtd     = 1.0;
  double fwgt      = 1.0;
  double tboostall = 0.0;
  double temp0     = gpu->sim.gamd_temp0;
  double ONE_KB    = 1.0 / (temp0 * KB);

  if (gpu->sim.GaMDNumLag == 0) {

    // Calculate GaMD weight, seting dihedral boost (tboost) to zero for now
    double EthreshP = gpu->sim.gamd_EthreshP;
    double kP       = gpu->sim.gamd_kP;
    double EV       = 0.0;

    if ((gpu->sim.igamd == 2)||(gpu->sim.igamd == 3)||(gpu->sim.igamd == 5)) {
      tboost = gpu->sim.GaMDtboost;
      fwgtd = gpu->pbGaMDfwgtd->_pSysData[0];
    }
    double totalenergy = *pot_ene_tot + tboost;
    EV = (EthreshP - totalenergy);
    if (((gpu->sim.igamd == 1) || (gpu->sim.igamd == 3) || (gpu->sim.igamd == 4) ||
         (gpu->sim.igamd == 5)) && (totalenergy < EthreshP)) {
      tboostall = 0.5 * kP * (EV * EV);
      fwgt = 1.0 - kP * EV;
    }

    // Calculate GaMD weight
    kCalculateGAMDWeightAndScaleForces(gpu, *pot_ene_tot, *dih_ene_tot, fwgt);
  }

  int numrecs = gpu->sim.GaMDNumRecs * 6;
  if ((gpu->sim.GaMDNumRecs+1)  >= gpu->sim.gamd_print_interval) {
    gpu->pGaMDWeightsAndEnergy[0] = *pot_ene_tot;
    gpu->pGaMDWeightsAndEnergy[1] = *dih_ene_tot;
    gpu->pGaMDWeightsAndEnergy[2] = fwgt;
    gpu->pGaMDWeightsAndEnergy[3] = fwgtd;
    gpu->pGaMDWeightsAndEnergy[4] = tboostall;
    gpu->pGaMDWeightsAndEnergy[5] = tboost;
  }
  if (gpu->sim.GaMDNumLag == gpu->sim.igamdlag) {
    gpu->sim.GaMDNumLag  = 0;
  }
  else {
    gpu->sim.GaMDNumLag++;
  }
  gpu->sim.GaMDNumRecs++;
  if (gpu->sim.GaMDNumRecs  >= gpu->sim.gamd_print_interval) {
    gpu->sim.GaMDNumRecs = 0 ;
  }
  *num_gamd_lag = gpu->sim.GaMDNumLag;
  *gamd_ene_tot = tboost + tboostall;
}

//---------------------------------------------------------------------------------------------
// gpu_calculate_and_apply_gamd_weights_nb_: analogous to the function in standard aMD for
//                                        carrying out the weighting and boosting protocol.
//
// Arguments:
//   pot_ene_tot:    total potential energy
//   dih_ene_tot:    total dihedral energy
//   gamd_ene_tot:   total GaMD boost potential
//   num_gamd_lag:   counter to track when the gamdlag phase ends
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_and_apply_gamd_weights_nb_(double* pot_ene_nb,
                                                         double* dih_ene_tot,
                                                         double* gamd_ene_tot,
                                                         double* num_gamd_lag)
{
  PRINTMETHOD("gpu_calculate_and_apply_gamd_weights");
  double tboost    = 0.0;
  double fwgtd     = 1.0;
  double fwgt      = 1.0;
  double tboostall = 0.0;
  double temp0     = gpu->sim.gamd_temp0;
  double ONE_KB    = 1.0 / (temp0 * KB);

  if (gpu->sim.GaMDNumLag == 0) {

    // Calculate GaMD weight, seting dihedral boost (tboost) to zero for now
    double EthreshP = gpu->sim.gamd_EthreshP;
    double kP       = gpu->sim.gamd_kP;
    double EV       = 0.0;

    if ((gpu->sim.igamd == 2) || (gpu->sim.igamd == 3) || (gpu->sim.igamd == 5) ||
	(gpu->sim.igamd == 20)) {
      tboost = gpu->sim.GaMDtboost;
      fwgtd = gpu->pbGaMDfwgtd->_pSysData[0];
    }
    double totalenergy = *pot_ene_nb;
    EV = (EthreshP - totalenergy);
    if (((gpu->sim.igamd == 1) || (gpu->sim.igamd == 3) || (gpu->sim.igamd == 4) ||
	 (gpu->sim.igamd == 5) || (gpu->sim.igamd == 20)) && (totalenergy < EthreshP)) {
      tboostall = 0.5 * kP * (EV * EV);
      fwgt = 1.0 - kP * EV;
    }

    // Calculate GaMD weight
    kCalculateGAMDWeightAndScaleForces_nb(gpu, *pot_ene_nb, *dih_ene_tot, fwgt);
  }

  int numrecs = gpu->sim.GaMDNumRecs * 6;
  if ((gpu->sim.GaMDNumRecs+1)  >= gpu->sim.gamd_print_interval) {
    gpu->pGaMDWeightsAndEnergy[0] = *pot_ene_nb;
    gpu->pGaMDWeightsAndEnergy[1] = *dih_ene_tot;
    gpu->pGaMDWeightsAndEnergy[2] = fwgt;
    gpu->pGaMDWeightsAndEnergy[3] = fwgtd;
    gpu->pGaMDWeightsAndEnergy[4] = tboostall;
    gpu->pGaMDWeightsAndEnergy[5] = tboost;
  }
  if (gpu->sim.GaMDNumLag == gpu->sim.igamdlag) {
    gpu->sim.GaMDNumLag  = 0;
  }
  else {
    gpu->sim.GaMDNumLag++;
  }
  gpu->sim.GaMDNumRecs++;
  if (gpu->sim.GaMDNumRecs  >= gpu->sim.gamd_print_interval) {
    gpu->sim.GaMDNumRecs = 0 ;
  }
  *num_gamd_lag = gpu->sim.GaMDNumLag;
  *gamd_ene_tot = tboost + tboostall;

}

//---------------------------------------------------------------------------------------------
// gpu_calculate_and_apply_gamd_weights_gb_: analogous to the function in standard aMD for
//                                        carrying out the weighting and boosting protocol.
//
// Arguments:
//   pot_ene_tot:    total potential energy
//   dih_ene_tot:    total dihedral energy
//   gamd_ene_tot:   total GaMD boost potential
//   num_gamd_lag:   counter to track when the gamdlag phase ends
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_and_apply_gamd_weights_gb_(double* pot_ene_tot, double* dih_ene_tot,
                                                      double* gamd_ene_tot,
                                                      double* num_gamd_lag)
{
  PRINTMETHOD("gpu_calculate_and_apply_gamd_weights");
  double tboost    = 0.0;
  double fwgtd     = 1.0;
  double fwgt      = 1.0;
  double tboostall = 0.0;
  double temp0     = gpu->sim.gamd_temp0;
  double ONE_KB    = 1.0 / (temp0 * KB);

  if (gpu->sim.GaMDNumLag == 0) {

    // Calculate GaMD weight, seting dihedral boost (tboost) to zero for now
    double EthreshP = gpu->sim.gamd_EthreshP;
    double kP       = gpu->sim.gamd_kP;
    double EV       = 0.0;

    if ((gpu->sim.igamd == 2)||(gpu->sim.igamd == 3)||(gpu->sim.igamd == 5)) {
      tboost = gpu->sim.GaMDtboost;
      fwgtd = gpu->pbGaMDfwgtd->_pSysData[0];
    }
    double totalenergy = *pot_ene_tot + tboost;
    EV = (EthreshP - totalenergy);
    if (((gpu->sim.igamd == 1) || (gpu->sim.igamd == 3) || (gpu->sim.igamd == 4) ||
         (gpu->sim.igamd == 5)) && (totalenergy < EthreshP)) {
      tboostall = 0.5 * kP * (EV * EV);
      fwgt = 1.0 - kP * EV;
    }

    // Calculate GaMD weight
    kCalculateGAMDWeightAndScaleForces_gb(gpu, *pot_ene_tot, *dih_ene_tot, fwgt);
  }

  int numrecs = gpu->sim.GaMDNumRecs * 6;
  if ((gpu->sim.GaMDNumRecs+1)  >= gpu->sim.gamd_print_interval) {
    gpu->pGaMDWeightsAndEnergy[0] = *pot_ene_tot;
    gpu->pGaMDWeightsAndEnergy[1] = *dih_ene_tot;
    gpu->pGaMDWeightsAndEnergy[2] = fwgt;
    gpu->pGaMDWeightsAndEnergy[3] = fwgtd;
    gpu->pGaMDWeightsAndEnergy[4] = tboostall;
    gpu->pGaMDWeightsAndEnergy[5] = tboost;
  }
  if (gpu->sim.GaMDNumLag == gpu->sim.igamdlag) {
    gpu->sim.GaMDNumLag  = 0;
  }
  else {
    gpu->sim.GaMDNumLag++;
  }
  gpu->sim.GaMDNumRecs++;
  if (gpu->sim.GaMDNumRecs  >= gpu->sim.gamd_print_interval) {
    gpu->sim.GaMDNumRecs = 0 ;
  }
  *num_gamd_lag = gpu->sim.GaMDNumLag;
  *gamd_ene_tot = tboost + tboostall;
}

//---------------------------------------------------------------------------------------------
// NOTE: This funciton should not be used for the moment as pNBForce is not saved separately with GB.
//
// gpu_calculate_and_apply_gamd_weights_gb_nb_: analogous to the function in standard aMD for
//                                        carrying out the weighting and boosting protocol.
//
// Arguments:
//   pot_ene_tot:    total potential energy
//   dih_ene_tot:    total dihedral energy
//   gamd_ene_tot:   total GaMD boost potential
//   num_gamd_lag:   counter to track when the gamdlag phase ends
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_and_apply_gamd_weights_gb_nb_(double* pot_ene_nb,
                                                         double* dih_ene_tot,
                                                         double* gamd_ene_tot,
                                                         double* num_gamd_lag)
{
  PRINTMETHOD("gpu_calculate_and_apply_gamd_weights");
  double tboost    = 0.0;
  double fwgtd     = 1.0;
  double fwgt      = 1.0;
  double tboostall = 0.0;
  double temp0     = gpu->sim.gamd_temp0;
  double ONE_KB    = 1.0 / (temp0 * KB);

  if (gpu->sim.GaMDNumLag == 0) {

    // Calculate GaMD weight, seting dihedral boost (tboost) to zero for now
    double EthreshP = gpu->sim.gamd_EthreshP;
    double kP       = gpu->sim.gamd_kP;
    double EV       = 0.0;

    if ((gpu->sim.igamd == 2) || (gpu->sim.igamd == 3) || (gpu->sim.igamd == 5) ||
        (gpu->sim.igamd == 20)) {
      tboost = gpu->sim.GaMDtboost;
      fwgtd = gpu->pbGaMDfwgtd->_pSysData[0];
    }
    double totalenergy = *pot_ene_nb;
    EV = (EthreshP - totalenergy);
    if (((gpu->sim.igamd == 1) || (gpu->sim.igamd == 3) || (gpu->sim.igamd == 4) ||
         (gpu->sim.igamd == 5) || (gpu->sim.igamd == 20)) && (totalenergy < EthreshP)) {
      tboostall = 0.5 * kP * (EV * EV);
      fwgt = 1.0 - kP * EV;
    }

    // Calculate GaMD weight
    kCalculateGAMDWeightAndScaleForces_gb_nb(gpu, *pot_ene_nb, *dih_ene_tot, fwgt);
  }

  int numrecs = gpu->sim.GaMDNumRecs * 6;
  if ((gpu->sim.GaMDNumRecs+1)  >= gpu->sim.gamd_print_interval) {
    gpu->pGaMDWeightsAndEnergy[0] = *pot_ene_nb;
    gpu->pGaMDWeightsAndEnergy[1] = *dih_ene_tot;
    gpu->pGaMDWeightsAndEnergy[2] = fwgt;
    gpu->pGaMDWeightsAndEnergy[3] = fwgtd;
    gpu->pGaMDWeightsAndEnergy[4] = tboostall;
    gpu->pGaMDWeightsAndEnergy[5] = tboost;
  }
  if (gpu->sim.GaMDNumLag == gpu->sim.igamdlag) {
    gpu->sim.GaMDNumLag  = 0;
  }
  else {
    gpu->sim.GaMDNumLag++;
  }
  gpu->sim.GaMDNumRecs++;
  if (gpu->sim.GaMDNumRecs  >= gpu->sim.gamd_print_interval) {
    gpu->sim.GaMDNumRecs = 0 ;
  }
  *num_gamd_lag = gpu->sim.GaMDNumLag;
  *gamd_ene_tot = tboost + tboostall;

}

//---------------------------------------------------------------------------------------------
// gpu_calculate_and_apply_gamd_weights_bond_: analogous to the function in standard aMD for
//                                        carrying out the weighting and boosting protocol.
//
// Arguments:
//   pot_ene_tot:    total potential energy
//   dih_ene_tot:    total dihedral energy
//   gamd_ene_tot:   total GaMD boost potential
//   num_gamd_lag:   counter to track when the gamdlag phase ends
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_and_apply_gamd_weights_bond_(double* pot_ene_tot,
                                                         double* dih_ene_tot,
                                                         double* gamd_ene_tot,
                                                         double* num_gamd_lag)
{
  PRINTMETHOD("gpu_calculate_and_apply_gamd_weights_bond");
  double tboost    = 0.0;
  double fwgtd     = 1.0;
  double tboostall = 0.0;
  double fwgt      = 1.0;
  double temp0     = gpu->sim.gamd_temp0;
  double ONE_KB    = 1.0 / (temp0 * KB);

 if (gpu->sim.GaMDNumLag == 0) {
   // Calculate GaMD weight, seting dihedral boost (tboost) to zero for now
    double EthreshP = gpu->sim.gamd_EthreshD;
    double kP       = gpu->sim.gamd_kD;
    double EV       = 0.0;

    double totalenergy = *dih_ene_tot ;
    EV = (EthreshP - totalenergy);
    if (((gpu->sim.igamd == 20)) && (totalenergy < EthreshP)) {
      tboost = 0.5 * kP * (EV * EV);
      fwgtd = 1.0 - kP * EV;

      // Calculate GaMD weight
      kCalculateGAMDWeightAndScaleForces_bond(gpu, *pot_ene_tot, *dih_ene_tot, fwgtd);
    }
 }

  int numrecs = gpu->sim.GaMDNumRecs * 6;
  if ((gpu->sim.GaMDNumRecs+1)  >= gpu->sim.gamd_print_interval) {
    gpu->pGaMDWeightsAndEnergy[0] = *pot_ene_tot;
    gpu->pGaMDWeightsAndEnergy[1] = *dih_ene_tot;
    gpu->pGaMDWeightsAndEnergy[2] = fwgt;
    gpu->pGaMDWeightsAndEnergy[3] = fwgtd;
    gpu->pGaMDWeightsAndEnergy[4] = tboostall;
    gpu->pGaMDWeightsAndEnergy[5] = tboost;
  }
  if (gpu->sim.GaMDNumLag == gpu->sim.igamdlag) {
    gpu->sim.GaMDNumLag  = 0;
  }
  else {
    gpu->sim.GaMDNumLag++;
  }
  gpu->sim.GaMDNumRecs++;
  if (gpu->sim.GaMDNumRecs  >= gpu->sim.gamd_print_interval) {
    gpu->sim.GaMDNumRecs = 0 ;
  }
  gpu->sim.GaMDtboost = tboost;
  gpu->pbGaMDfwgtd->_pSysData[0] = fwgtd;
  gpu->pbGaMDfwgtd->Upload();
  *num_gamd_lag = gpu->sim.GaMDNumLag;
  *gamd_ene_tot = *gamd_ene_tot + tboost + tboostall;

}


//---------------------------------------------------------------------------------------------
// gpu_calculate_and_apply_gamd_weights_ti_region_: analogous to the function in standard aMD for
//                                        carrying out the weighting and boosting protocol.
//
// Arguments:
//   pot_ene_tot:    total potential energy
//   dih_ene_tot:    total dihedral energy
//   gamd_ene_tot:   total GaMD boost potential
//   num_gamd_lag:   counter to track when the gamdlag phase ends
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_and_apply_gamd_weights_ti_region_(double* pot_ene_ti_region,
                                                         double* dih_ene_tot,
                                                         double* gamd_ene_tot,
                                                         double* fwgt,
                                                         double* num_gamd_lag)
{
  PRINTMETHOD("gpu_calculate_and_apply_gamd_weights_ti_region");
  *fwgt     = 1.0;
  double fwgtd     = 1.0;
  double tboost    = 0.0;
  double tboostall = 0.0;
  double temp0     = gpu->sim.gamd_temp0;
  double ONE_KB    = 1.0 / (temp0 * KB);

    // Calculate GaMD weight, seting dihedral boost (tboost) to zero for now
    double EthreshP = gpu->sim.gamd_EthreshP;
    double kP       = gpu->sim.gamd_kP;
    double EV       = 0.0;

    if (gpu->sim.igamd == 7 || gpu->sim.igamd == 9 || gpu->sim.igamd == 11) {
      tboost = gpu->sim.GaMDtboost;
      fwgtd = gpu->pbGaMDfwgtd->_pSysData[0];
    }
    double totalenergy = *pot_ene_ti_region;
    EV = (EthreshP - totalenergy);

    if (((gpu->sim.igamd == 6) || (gpu->sim.igamd == 7)) && (totalenergy < EthreshP)) {
      tboostall = 0.5 * kP * (EV * EV);
      *fwgt = 1.0 - kP * EV;

      // Calculate GaMD force scaling factor
      kCalculateGAMDWeightAndScaleForces_ti(gpu, *pot_ene_ti_region, *dih_ene_tot, *fwgt);
    }

    if (((gpu->sim.igamd == 8) || (gpu->sim.igamd == 9)) && (totalenergy < EthreshP)) {
      tboostall = 0.5 * kP * (EV * EV);
      *fwgt = 1.0 - kP * EV;

      // Calculate GaMD force scaling factor
      kCalculateGAMDWeightAndScaleForces_ti_bonded(gpu, *pot_ene_ti_region, *dih_ene_tot, *fwgt);
    }

    if (((gpu->sim.igamd == 10) || (gpu->sim.igamd == 11) || (gpu->sim.igamd == 100) || (gpu->sim.igamd == 101) || (gpu->sim.igamd == 102) || (gpu->sim.igamd == 103) || (gpu->sim.igamd == 104)) && (totalenergy < EthreshP)) {
      tboostall = 0.5 * kP * (EV * EV);
      *fwgt = 1.0 - kP * EV;

      // Calculate GaMD force scaling factor
      kCalculateGAMDWeightAndScaleForces_ti_nonbonded(gpu, *pot_ene_ti_region, *dih_ene_tot, *fwgt);
    }
    if ((gpu->sim.igamd == 12 || gpu->sim.igamd == 13 || gpu->sim.igamd == 16|| gpu->sim.igamd == 17) && (totalenergy < EthreshP)) {

      tboostall = 0.5 * kP * (EV * EV);
      *fwgt = 1.0 - kP * EV;

      if((gpu->sim.igamd == 16|| gpu->sim.igamd == 17)&& fabs(totalenergy)<=0.0001){
      tboostall =0.0;
      *fwgt = 1.0;
      }

      // Calculate GaMD force scaling factor
      kCalculateGAMDWeightAndScaleForces_sc_nonbonded(gpu, *pot_ene_ti_region, *dih_ene_tot, *fwgt);
 //   printf("Debug-GaMD-GPU-ti) ep,ev,tboostall,kP =%12.5f, %12.5f, %12.5f, %12.5f\n",totalenergy,EV,tboostall, kP);
    }

   if (((gpu->sim.igamd == 14 || gpu->sim.igamd == 15||gpu->sim.igamd ==111||gpu->sim.igamd ==113||gpu->sim.igamd ==115)) && (totalenergy < EthreshP)) {
      tboostall = 0.5 * kP * (EV * EV);
      *fwgt = 1.0 - kP * EV;

      // Calculate GaMD force scaling factor
      kCalculateGAMDWeightAndScaleForces_sc_ppi(gpu, *pot_ene_ti_region, *dih_ene_tot, *fwgt);
    }
   if ((gpu->sim.igamd == 117 || gpu->sim.igamd == 118||gpu->sim.igamd == 119) && (totalenergy < EthreshP)) {
      tboostall = 0.5 * kP * (EV * EV);
      *fwgt = 1.0 - kP * EV;

      // Calculate GaMD force scaling factor for test
      kCalculateGAMDWeightAndScaleForces_sc_bonded_ppi(gpu, *pot_ene_ti_region, *dih_ene_tot, *fwgt);
    }
     if (gpu->sim.igamd == 120  && (totalenergy < EthreshP)) {
      tboostall = 0.5 * kP * (EV * EV);
      *fwgt = 1.0 - kP * EV;
      // Calculate GaMD force scaling factor for test

    kCalculateGAMDWeightAndScaleForces_sc_bonded_ppi(gpu, *pot_ene_ti_region, *dih_ene_tot, *fwgt);
    }

     if (gpu->sim.igamd == 18  && (totalenergy < EthreshP)) {
         tboostall = 0.5 * kP * (EV * EV);
        *fwgt = 1.0 - kP * EV;
      // Calculate GaMD force scaling factor
       kCalculateGAMDWeightAndScaleForces_sc_nonbonded(gpu, *pot_ene_ti_region, *dih_ene_tot, *fwgt);
     }

    if (gpu->sim.igamd == 19  && (totalenergy < EthreshP)) {
         tboostall = 0.5 * kP * (EV * EV);
        *fwgt = 1.0 - kP * EV;
      // Calculate GaMD force scaling factor
       kCalculateGAMDWeightAndScaleForces_sc_nonbonded(gpu, *pot_ene_ti_region, *dih_ene_tot, *fwgt);
     }

    if (gpu->sim.igamd == 21  && (totalenergy < EthreshP)) {
         tboostall = 0.5 * kP * (EV * EV);
        *fwgt = 1.0 - kP * EV;
      // Calculate GaMD force scaling factor
       kCalculateGAMDWeightAndScaleForces_sc_nonbonded(gpu, *pot_ene_ti_region, *dih_ene_tot, *fwgt);
     }

    if ((gpu->sim.igamd == 110 ||gpu->sim.igamd == 112 ||gpu->sim.igamd == 114 ||gpu->sim.igamd == 116) && (totalenergy < EthreshP)) {
      tboostall = 0.5 * kP * (EV * EV);
      *fwgt = 1.0 - kP * EV;
      // Calculate GaMD force scaling factor
      kCalculateGAMDWeightAndScaleForces_sc_nonbonded(gpu, *pot_ene_ti_region, *dih_ene_tot, *fwgt);
     }
    gpu->pGaMDWeightsAndEnergy[0] = *pot_ene_ti_region;
    gpu->pGaMDWeightsAndEnergy[1] = *dih_ene_tot;
    gpu->pGaMDWeightsAndEnergy[2] = *fwgt;
    gpu->pGaMDWeightsAndEnergy[3] = fwgtd;
    gpu->pGaMDWeightsAndEnergy[4] = tboostall;
    gpu->pGaMDWeightsAndEnergy[5] = tboost;
  *gamd_ene_tot = tboost + tboostall;
   // printf("Debug-GaMD-GPU-ti) ep,ev,tboostall,fwgt=%12.5f, %12.5f, %12.5f, %12.5f\n",totalenergy,EV,tboostall, gpu->pGaMDWeightsAndEnergy[3]);
  // printf("Debug-GaMD-GPU-ti) gpu->pGaMDWeightsAndEnergy = %12.5f, %12.5f, %12.5f, %12.5f, %12.5f, %12.5f\n", gpu->pGaMDWeightsAndEnergy[0],gpu->pGaMDWeightsAndEnergy[1],gpu->pGaMDWeightsAndEnergy[2],gpu->pGaMDWeightsAndEnergy[3],gpu->pGaMDWeightsAndEnergy[4],gpu->pGaMDWeightsAndEnergy[5]);

}

//---------------------------------------------------------------------------------------------
// gpu_calculate_and_apply_gamd_weights_ti_others_: analogous to the function in standard aMD for
//                                        carrying out the weighting and boosting protocol.
//
// Arguments:
//   pot_ene_tot:    total potential energy
//   dih_ene_tot:    total dihedral energy
//   gamd_ene_tot:   total GaMD boost potential
//   num_gamd_lag:   counter to track when the gamdlag phase ends
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_and_apply_gamd_weights_ti_others_(double* pot_ene_tot, double* dih_ene_tot,
                                                      double* gamd_ene_tot,
                                                      double* num_gamd_lag)
{
  PRINTMETHOD("gpu_calculate_and_apply_gamd_weights_ti_others");
  double tboost    = 0.0;
  double fwgtd     = 1.0;
  double tboostall = 0.0;
  double fwgt      = 1.0;
  double temp0     = gpu->sim.gamd_temp0;
  double ONE_KB    = 1.0 / (temp0 * KB);

    // Calculate GaMD weight, seting dihedral boost (tboost) to zero for now
    double EthreshP = gpu->sim.gamd_EthreshD;
    double kP       = gpu->sim.gamd_kD;
    double EV       = 0.0;

    double totalenergy = *dih_ene_tot ;
    EV = (EthreshP - totalenergy);
    if (((gpu->sim.igamd == 7)) && (totalenergy < EthreshP)) {
      tboost = 0.5 * kP * (EV * EV);
      fwgtd = 1.0 - kP * EV;

      // Calculate GaMD weight
      kCalculateGAMDWeightAndScaleForces_ti_others(gpu, *pot_ene_tot, *dih_ene_tot, fwgtd);
    }

    if (((gpu->sim.igamd == 9)) && (totalenergy < EthreshP)) {
      tboost = 0.5 * kP * (EV * EV);
      fwgtd = 1.0 - kP * EV;

      // Calculate GaMD weight
      kCalculateGAMDWeightAndScaleForces_ti_bonded_others(gpu, *pot_ene_tot, *dih_ene_tot,
							  fwgtd);
    }

    if (((gpu->sim.igamd == 11 || gpu->sim.igamd == 100)) && (totalenergy < EthreshP)) {
      tboost = 0.5 * kP * (EV * EV);
      fwgtd = 1.0 - kP * EV;

      // Calculate GaMD weight
      kCalculateGAMDWeightAndScaleForces_ti_nonbonded_others(gpu, *pot_ene_tot, *dih_ene_tot,
							     fwgtd);
    }

    if (((gpu->sim.igamd == 15)) && (totalenergy < EthreshP)) {
      tboost = 0.5 * kP * (EV * EV);
      fwgtd = 1.0 - kP * EV;

      // Calculate GaMD weight
      kCalculateGAMDWeightAndScaleForces_sc_others_ppi(gpu, *pot_ene_tot, *dih_ene_tot, fwgtd);
    }

    if (((gpu->sim.igamd == 18 ||gpu->sim.igamd == 21)) && (totalenergy < EthreshP)) {
      tboost = 0.5 * kP * (EV * EV);
      fwgtd = 1.0 - kP * EV;

      // Calculate GaMD weight
      kCalculateGAMDWeightAndScaleForces_sc_bonded_ppi(gpu, *pot_ene_tot, *dih_ene_tot, fwgtd);
    }


    gpu->pGaMDWeightsAndEnergy[0] = *pot_ene_tot;
    gpu->pGaMDWeightsAndEnergy[1] = *dih_ene_tot;
    // gpu->pGaMDWeightsAndEnergy[2] = fwgt;
    gpu->pGaMDWeightsAndEnergy[3] = fwgtd;
    // gpu->pGaMDWeightsAndEnergy[4] = tboostall;
    gpu->pGaMDWeightsAndEnergy[5] = tboost;

  *gamd_ene_tot = *gamd_ene_tot + tboost + tboostall;
  // printf("Debug-GaMD-GPU-ti-others) gpu->pGaMDWeightsAndEnergy = %12.5f, %12.5f, %12.5f, %12.5f, %12.5f, %12.5f\n", gpu->pGaMDWeightsAndEnergy[0],gpu->pGaMDWeightsAndEnergy[1],gpu->pGaMDWeightsAndEnergy[2],gpu->pGaMDWeightsAndEnergy[3],gpu->pGaMDWeightsAndEnergy[4],gpu->pGaMDWeightsAndEnergy[5]);
}

//---------------------------------------------------------------------------------------------
// gpu_calculate_and_apply_gamd_weights_sc_dual_nonbonded: analogous to the function in standard aMD for
//                                        carrying out the weighting and boosting protocol.
//
// Arguments:
//   pot_ene_tot:    total potential energy
//   dih_ene_tot:    total dihedral energy
//   gamd_ene_tot:   total GaMD boost potential
//   num_gamd_lag:   counter to track when the gamdlag phase ends
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_and_apply_gamd_weights_sc_dual_nonbonded_(double* pot_ene_tot, double* dih_ene_tot,
                                                      double* gamd_ene_tot,
                                                      double* num_gamd_lag)
{
  PRINTMETHOD("gpu_calculate_and_apply_gamd_weights_sc_dual");
  double tboost    = 0.0;
  double fwgtd     = 1.0;
  double tboostall = 0.0;
  double fwgt      = 1.0;
  double temp0     = gpu->sim.gamd_temp0;
  double ONE_KB    = 1.0 / (temp0 * KB);

    // Calculate GaMD weight, seting dihedral boost (tboost) to zero for now
    double EthreshP_D = gpu->sim.gamd_EthreshD;
    double kP_D       = gpu->sim.gamd_kD;
    double EV_D       = 0.0;
    double EthreshP_P = gpu->sim.gamd_EthreshP;
    double kP_P       = gpu->sim.gamd_kP;
    double EV_P       = 0.0;

    double totalenergy_D = *dih_ene_tot;
    double totalenergy_P = *pot_ene_tot;
    EV_D = (EthreshP_D - totalenergy_D);
    EV_P = (EthreshP_P - totalenergy_P);

   if (gpu->sim.igamd == 13||gpu->sim.igamd == 17) {
      if(totalenergy_P < EthreshP_P) {
      tboostall = 0.5 * kP_P * (EV_P * EV_P);
      fwgt = 1.0 - kP_P * EV_P;
      }
      if(fabs(totalenergy_P)<=0.01&&gpu->sim.igamd==17){
      tboostall =0;
      fwgt = 1.0;
      }
      if(totalenergy_D < EthreshP_D) {
      tboost = 0.5 * kP_D * (EV_D * EV_D);
      fwgtd = 1.0 - kP_D * EV_D;
      }

      // Calculate GaMD weight
      kCalculateGAMDWeightAndScaleForces_sc_dual_nonbonded_ppi(gpu, *pot_ene_tot, *dih_ene_tot, fwgt, fwgtd);
    }
    gpu->pGaMDWeightsAndEnergy[0] = *pot_ene_tot;
    gpu->pGaMDWeightsAndEnergy[1] = *dih_ene_tot;
    gpu->pGaMDWeightsAndEnergy[2] = fwgt;
    gpu->pGaMDWeightsAndEnergy[3] = fwgtd;
    gpu->pGaMDWeightsAndEnergy[4] = tboostall;
    gpu->pGaMDWeightsAndEnergy[5] = tboost;

    *gamd_ene_tot = tboost + tboostall;

  // printf("Debug-GaMD-GPU-ti-others) gpu->pGaMDWeightsAndEnergy = %12.5f, %12.5f, %12.5f, %12.5f, %12.5f, %12.5f\n", gpu->pGaMDWeightsAndEnergy[0],gpu->pGaMDWeightsAndEnergy[1],gpu->pGaMDWeightsAndEnergy[2],gpu->pGaMDWeightsAndEnergy[3],gpu->pGaMDWeightsAndEnergy[4],gpu->pGaMDWeightsAndEnergy[5]);
}

//---------------------------------------------------------------------------------------------
// gpu_calculate_and_apply_gamd_weights_sc_dual_bonded: analogous to the function in standard aMD for
//                                        carrying out the weighting and boosting protocol.
//
// Arguments:
//   pot_ene_tot:    total potential energy
//   dih_ene_tot:    total dihedral energy
//   gamd_ene_tot:   total GaMD boost potential
//   num_gamd_lag:   counter to track when the gamdlag phase ends
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_and_apply_gamd_weights_sc_dual_bonded_(double* pot_ene_tot, double* dih_ene_tot,
                                                      double* gamd_ene_tot,
                                                      double* num_gamd_lag)
{
  PRINTMETHOD("gpu_calculate_and_apply_gamd_weights_sc_dual");
  double tboost    = 0.0;
  double fwgtd     = 1.0;
  double tboostall = 0.0;
  double fwgt      = 1.0;
  double temp0     = gpu->sim.gamd_temp0;
  double ONE_KB    = 1.0 / (temp0 * KB);

    // Calculate GaMD weight, seting dihedral boost (tboost) to zero for now
    double EthreshP_D = gpu->sim.gamd_EthreshD;
    double kP_D       = gpu->sim.gamd_kD;
    double EV_D       = 0.0;
    double EthreshP_P = gpu->sim.gamd_EthreshP;
    double kP_P       = gpu->sim.gamd_kP;
    double EV_P       = 0.0;

    double totalenergy_D = *dih_ene_tot;
    double totalenergy_P = *pot_ene_tot;
    EV_D = (EthreshP_D - totalenergy_D);
    EV_P = (EthreshP_P - totalenergy_P);

   if (gpu->sim.igamd == 15||gpu->sim.igamd==18) {
      if(totalenergy_P < EthreshP_P) {
      tboostall = 0.5 * kP_P * (EV_P * EV_P);
      fwgt = 1.0 - kP_P * EV_P;
      }
      if(totalenergy_D < EthreshP_D) {
      tboost = 0.5 * kP_D * (EV_D * EV_D);
      fwgtd = 1.0 - kP_D * EV_D;
      }

      // Calculate GaMD weight
      kCalculateGAMDWeightAndScaleForces_sc_dual_bonded_ppi(gpu, *pot_ene_tot, *dih_ene_tot, fwgt, fwgtd);
    }
    gpu->pGaMDWeightsAndEnergy[0] = *pot_ene_tot;
    gpu->pGaMDWeightsAndEnergy[1] = *dih_ene_tot;
    gpu->pGaMDWeightsAndEnergy[2] = fwgt;
    gpu->pGaMDWeightsAndEnergy[3] = fwgtd;
    gpu->pGaMDWeightsAndEnergy[4] = tboostall;
    gpu->pGaMDWeightsAndEnergy[5] = tboost;

    *gamd_ene_tot = tboost + tboostall;

  // printf("Debug-GaMD-GPU-ti-others) gpu->pGaMDWeightsAndEnergy = %12.5f, %12.5f, %12.5f, %12.5f, %12.5f, %12.5f\n", gpu->pGaMDWeightsAndEnergy[0],gpu->pGaMDWeightsAndEnergy[1],gpu->pGaMDWeightsAndEnergy[2],gpu->pGaMDWeightsAndEnergy[3],gpu->pGaMDWeightsAndEnergy[4],gpu->pGaMDWeightsAndEnergy[5]);
}

//---------------------------------------------------------------------------------------------
// gpu_calculate_and_apply_gamd_weights_sc_dual_bonded: analogous to the function in standard aMD for
//                                        carrying out the weighting and boosting protocol.
//
// Arguments:
//   pot_ene_tot:    total potential energy
//   dih_ene_tot:    total dihedral energy
//  ppi_bond_tot:    total bonded energy
//   gamd_ene_tot:   total GaMD boost potential
//   num_gamd_lag:   counter to track when the gamdlag phase ends
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_and_apply_gamd_weights_sc_triple_(double* pot_ene_tot, double* dih_ene_tot,
                                                      double* ppi_bond_tot,double* gamd_ene_tot,
                                                      double* num_gamd_lag)
{
  PRINTMETHOD("gpu_calculate_and_apply_gamd_weights_sc_dual");
  double tboost    = 0.0;
  double fwgtd     = 1.0;
  double fwgtb     = 1.0;
  double tboostall = 0.0;
  double tboostb   = 0.0;
  double fwgt      = 1.0;
  double temp0     = gpu->sim.gamd_temp0;
  double ONE_KB    = 1.0 / (temp0 * KB);

    // Calculate GaMD weight, seting dihedral boost (tboost) to zero for now
    double EthreshP_D = gpu->sim.gamd_EthreshD;
    double kP_D       = gpu->sim.gamd_kD;
    double EV_D       = 0.0;
    double EthreshP_P = gpu->sim.gamd_EthreshP;
    double kP_P       = gpu->sim.gamd_kP;
    double EV_P       = 0.0;
    double EthreshP_B = gpu->sim.gamd_EthreshB;
    double kP_B       = gpu->sim.gamd_kB;
    double EV_B       = 0.0;
    double totalenergy_D = *dih_ene_tot;
    double totalenergy_P = *pot_ene_tot;
    double totalenergy_B = *ppi_bond_tot;
    EV_D = (EthreshP_D - totalenergy_D);
    EV_P = (EthreshP_P - totalenergy_P);
    EV_B = (EthreshP_B - totalenergy_B);


   if (gpu->sim.igamd == 22 || gpu->sim.igamd==21||gpu->sim.igamd==20||gpu->sim.igamd==23||gpu->sim.igamd==24||gpu->sim.igamd==25||gpu->sim.igamd==26||gpu->sim.igamd == 27||gpu->sim.igamd==28) {
      if(totalenergy_P < EthreshP_P) {
      tboostall = 0.5 * kP_P * (EV_P * EV_P);
      fwgt = 1.0 - kP_P * EV_P;
      }
      if(totalenergy_D < EthreshP_D) {
      tboost = 0.5 * kP_D * (EV_D * EV_D);
      fwgtd = 1.0 - kP_D * EV_D;
      }
      if(totalenergy_B < EthreshP_B) {
      tboostb = 0.5 * kP_B * (EV_B * EV_B);
      fwgtb = 1.0 - kP_B * EV_B;
      }
      if(fabs(totalenergy_P)<0.001&&fabs(totalenergy_B)<0.001&&gpu->sim.igamd==20){
      tboostall =0.0;
      fwgt = 1.0;
      tboostb = 0.0;
      fwgtb = 1.0;
      }
      if(fabs(totalenergy_P)<0.001&&gpu->sim.igamd==26){
      tboostall =0.0;
      fwgt = 1.0;
      }
      if(gpu->sim.igamd == 22){
      // Calculate GaMD weight
      kCalculateGAMDWeightAndScaleForces_sc_triple_ppi(gpu, *pot_ene_tot, *dih_ene_tot, *ppi_bond_tot,fwgt, fwgtd,fwgtb);
      }
      if (gpu->sim.igamd == 21||gpu->sim.igamd == 26||gpu->sim.igamd == 28){
      kCalculateGAMDWeightAndScaleForces_sc_triple2_ppi(gpu, *pot_ene_tot, *dih_ene_tot, *ppi_bond_tot,fwgt, fwgtd,fwgtb);
      }
      if (gpu->sim.igamd == 20||gpu->sim.igamd == 23){
      kCalculateGAMDWeightAndScaleForces_sc_triple3_ppi(gpu, *pot_ene_tot, *dih_ene_tot, *ppi_bond_tot,fwgt, fwgtd,fwgtb);
      }
      if (gpu->sim.igamd == 24||gpu->sim.igamd == 25){
      kCalculateGAMDWeightAndScaleForces_sc_triple4_ppi(gpu, *pot_ene_tot, *dih_ene_tot, *ppi_bond_tot,fwgt, fwgtd,fwgtb);
      }
      if (gpu->sim.igamd == 27){
      kCalculateGAMDWeightAndScaleForces_sc_triple5_ppi(gpu, *pot_ene_tot, *dih_ene_tot, *ppi_bond_tot,fwgt, fwgtd,fwgtb);
      }
    }
    gpu->pGaMDWeightsAndEnergy[0] = *pot_ene_tot;
    gpu->pGaMDWeightsAndEnergy[1] = *dih_ene_tot;
    gpu->pGaMDWeightsAndEnergy[2] = fwgt;
    gpu->pGaMDWeightsAndEnergy[3] = fwgtd;
    gpu->pGaMDWeightsAndEnergy[4] = tboostall+tboostb;
    gpu->pGaMDWeightsAndEnergy[5] = tboost;

    *gamd_ene_tot = tboost + tboostall+tboostb;

    //printf("Debug-GaMD-GPU-triple tboostb,fwgtb,kP_B,EthreshP_B %12.5f, %12.5f, %12.5f, %12.5f\n",tboostb,fwgtb,kP_B,EthreshP_B);
   //printf("Debug-GaMD-GPU-ti-others) gpu->pGaMDWeightsAndEnergy = %12.5f, %12.5f, %12.5f, %12.5f, %12.5f, %12.5f\n", gpu->pGaMDWeightsAndEnergy[0],gpu->pGaMDWeightsAndEnergy[1],gpu->pGaMDWeightsAndEnergy[2],gpu->pGaMDWeightsAndEnergy[3],gpu->pGaMDWeightsAndEnergy[4],gpu->pGaMDWeightsAndEnergy[5]);
}

//---------------------------------------------------------------------------------------------
// GaMD functions for gb, no need to do neighbor lists here
//---------------------------------------------------------------------------------------------
#ifdef MPI
//---------------------------------------------------------------------------------------------
// gpu_calculate_gb_gamd_dihedral_energy_: the name says it all.
//
// Arguments:
//   totdih:     total dihedral energy
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_gb_gamd_dihedral_energy_(double* totdih)
{
  PRINTMETHOD("gpu_calculate_gb_gamd_dihedral_energy");

  // Local energy
  if (gpu->bCalculateLocalForces) {

    // Calculate GaMD energy for dihedral boost (tboost)
    kExecuteBondWorkUnits(gpu, 3);
    gpu->pbEnergyBuffer->Download();
    PMEDouble totdih2 = (PMEDouble)gpu->pbEnergyBuffer->_pSysData[GAMD_E_DIHEDRAL_OFFSET] /
                        ENERGYSCALE;
    *totdih = totdih2;
  }
  else {
    *totdih = 0.0;
  }
}

#else
//---------------------------------------------------------------------------------------------
// gpu_calculate_gb_gamd_dihedral_energy_weight_: the name says it all
//---------------------------------------------------------------------------------------------
extern "C" void gpu_calculate_gb_gamd_dihedral_energy_weight_(){
  PRINTMETHOD("gpu_calculate_gb_gamd_dihedral_energy_weight");
  PMEDouble totdih = 0.0;
  if (gpu->sim.GaMDNumLag == 0) {

    // Calculate AMD calculate energy for dihedral boost (tboost)
    kExecuteBondWorkUnits(gpu, 3);
    gpu->pbEnergyBuffer->Download();
    totdih = (PMEDouble)gpu->pbEnergyBuffer->_pSysData[GAMD_E_DIHEDRAL_OFFSET] / ENERGYSCALE;
  }

  // Calculate GaMD weight for dihedral boost (tboost)
  PMEDouble EthreshD = gpu->sim.gamd_EthreshD;
  PMEDouble kD       = gpu->sim.gamd_kD;
  PMEDouble tboost   = 0.0;
  PMEDouble fwgtd    = 1.0;
  PMEDouble EV       = 0.0;
  EV = (EthreshD - (totdih));
  if (EV > 0.0) {
    if (gpu->sim.GaMDNumLag == 0) {
      tboost = 0.5 * kD * EV * EV;
      fwgtd = 1.0 - kD * EV;
    }
  }
  gpu->sim.GaMDtboost = tboost;
  gpu->pbGaMDfwgtd->_pSysData[0] = fwgtd;
  gpu->pbGaMDfwgtd->Upload();
}
#endif

//---------------------------------------------------------------------------------------------
// gpu_scaledmd_setup_: setup for scaled MD simulations
//
// Arguments:
//   scaledMD:         flag to activate scaled MD
//   scaledMD_lambda:  mixing factor
//---------------------------------------------------------------------------------------------
extern "C" void gpu_scaledmd_setup_(int* scaledMD, double* scaledMD_lambda)
{
  PRINTMETHOD("gpu_scaledmd_setup");

  // Set up scaledmd parameters
  gpu->sim.scaledMD = *scaledMD;
  gpu->sim.scaledMD_lambda = *scaledMD_lambda;
  gpuCopyConstants();
  return;
}

//---------------------------------------------------------------------------------------------
// gpu_scaledmd_scale_frc_: scale the forces in a simulation based on the lambda value
//
// Arguments:
//   pot_ene_tot:  the total (unscaled) potential energy of the system
//---------------------------------------------------------------------------------------------
extern "C" void gpu_scaledmd_scale_frc_(double* pot_ene_tot)
{
  PRINTMETHOD("gpu_scaledmd_scale_frc");
  double lambda             = 1.0;
  lambda = gpu->sim.scaledMD_lambda;
  gpu->sim.scaledMD_energy =  *pot_ene_tot * lambda;
  gpu->sim.scaledMD_weight = -( *pot_ene_tot * (1.0-lambda));
  gpu->sim.scaledMD_unscaled_energy =  *pot_ene_tot;
  kScaledMDScaleForces(gpu, *pot_ene_tot, lambda);
}

//---------------------------------------------------------------------------------------------
// gpu_download_scaledmd_weights_: get the scaled MD weights from the GPU.  This will be called
//                                 by the master thread so that the log file can be written.
//
// Arguments:
//   scaledMD_energy:          the scaled MD energy
//   scaledMD_weight:          the scaled MD weight
//   scaledMD_unscaled_energy: the system energy without scaling
//---------------------------------------------------------------------------------------------
extern "C" void gpu_download_scaledmd_weights_(double* scaledMD_energy,
                                               double* scaledMD_weight,
                                               double* scaledMD_unscaled_energy)
{
  PRINTMETHOD("gpu_download_scaledmd_weights");
  *scaledMD_energy = gpu->sim.scaledMD_energy;
  *scaledMD_weight = gpu->sim.scaledMD_weight;
  *scaledMD_unscaled_energy = gpu->sim.scaledMD_unscaled_energy;
}

//---------------------------------------------------------------------------------------------
// gpu_download_ucell: grabs the ucell dimension
//
// Arguments:
//   ucell: returns ucell
//---------------------------------------------------------------------------------------------
extern "C" void gpu_download_ucell_(double ucell[3][3])
{
  for (int i = 0; i < 3; i++)
  {
      for (int j = 0; j < 3; j++)
      {
          ucell[j][i]                     = (double)gpu->sim.ucell[i][j];
      }
  }
}

//---------------------------------------------------------------------------------------------
// gpu_update_natoms: update the number of atoms that will be active in constant pH
//                    simulations.
//
// Arguments:
//   natoms:       the number of active atoms
//   NeighborList: flag to indicate whether a neighbor list is available
//---------------------------------------------------------------------------------------------
extern "C" void gpu_update_natoms_(int* natoms, bool* NeighborList)
{
  PRINTMETHOD("gpu_update_natoms");

  // Allocate system based on atom count
  gpu->sim.atoms = *natoms;
  gpu->bNeighborList = *NeighborList;
  gpu->sim.paddedNumberOfAtoms =
    ((*natoms + gpu->sim.grid - 1) >> gpu->sim.gridBits) << gpu->sim.gridBits;

  // Calculate stride to ensure buffers begin on safe texture boundaries
  gpu->sim.stride = ((*natoms + 63) >> 6) << 6;
  gpu->sim.stride2 = 2 * gpu->sim.stride;
  gpu->sim.stride3 = 3 * gpu->sim.stride;
  gpu->sim.stride4 = 4 * gpu->sim.stride;
  gpu->sim.stride5 = 5 * gpu->sim.stride;
  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gpu_set_first_update_atom_: set the first atom that will be updated in constant pH
//                             simulations.
//
// Arguments:
//   first:    the first atom that will get updated
//---------------------------------------------------------------------------------------------
extern "C" void gpu_set_first_update_atom_(int* first)
{
  PRINTMETHOD("gpu_set_first_update_atom");
  gpu->sim.first_update_atom = *first;
  gpuCopyConstants();
}

//---------------------------------------------------------------------------------------------
// gpu_ti_exchange_$VECTOR: functions to synchronize $VECTOR between regions 1 and 2 for
//                          linear atoms
//---------------------------------------------------------------------------------------------
extern "C" void gpu_ti_exchange_frc_() {
  PRINTMETHOD("gpu_ti_exchange_frc");
  kAFEExchangeFrc(gpu);
}

extern "C" void gpu_ti_exchange_vel_() {
  PRINTMETHOD("gpu_ti_exchange_vel");
  kAFEExchangeVel(gpu);
}

extern "C" void gpu_ti_exchange_crd_() {
  PRINTMETHOD("gpu_ti_exchange_crd");
  kAFEExchangeCrd(gpu);
}

//---------------------------------------------------------------------------------------------
// gpu_check_forces_: function for checking all forces for anomalies.  This launches an
//                    eponymous kernel.
//
// Arguments:
//   frc_chk:   indication of the stage at which the force is being checked (for reporting
//              purposes)
//
// This is a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_check_forces_(int *frc_chk)
{
  kCheckGpuForces(gpu, *frc_chk);
}

//---------------------------------------------------------------------------------------------
// gpu_check_consistency_: function for checking dynamics propagation for anomalies.  This
//                         launches an eponymous kernel.
//
// Arguments:
//   iter:      the iteration counter
//
// This is a debugging function.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_check_consistency_(int *iter, int *cchk, int *nchk)
{
  kCheckGpuConsistency(gpu, *iter, *cchk, *nchk);
}

//---------------------------------------------------------------------------------------------
// gpu_sumq_: gets the sums of the charges and squares of the charges for self in iphmd=3
//
// Arguments:
//   sumq:      the sum of charges
//   sumq2:     the sum of the squares of the charges
//
//---------------------------------------------------------------------------------------------
extern "C" void gpu_sumq_(double* sumq, double* sumq2)
{
  *sumq = gpu->sumq;
  *sumq2 = gpu->sumq2;
}

//---------------------------------------------------------------------------------------------
// GetGpuContextPtr: get a pointer to the global gpuContext in use by this library.
//---------------------------------------------------------------------------------------------
gpuContext GetGpuContextPtr()
{
  return gpu;
}
