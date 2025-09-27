#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// kCalculateCOM: launch the appropriate kernel to calculate the center of mass.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kCalculateCOM(gpuContext gpu)
{
  if (gpu->sim.soluteMolecules <= gpu->maxSoluteMolecules) {
    if (gpu->sim.ti_mode == 0) {
      kPMECalculateCOM_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    }
    else {
      kPMECalculateCOMAFE_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    }
  }
  else {
    if (gpu->sim.ti_mode == 0) {
      kPMECalculateCOMLarge_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    }
    else {
      kPMECalculateCOMLargeAFE_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
    }
  }
  LAUNCHERROR("kCalculateCOM");
}

//---------------------------------------------------------------------------------------------
// kReduceSoluteCOM_kernel: reduce the center of mass initially computed by kCalculateCOM in
//                          gpu_ntp_setup_() (see gpu.cpp).  This kernel is only called in
//                          that one instance.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(THREADS_PER_BLOCK, 1)
kReduceSoluteCOM_kernel()
{
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  while (pos < cSim.soluteMolecules) {
    PMEDouble invMass         = cSim.pSoluteInvMass[pos];
    PMEUllInt ullCOMX         = cSim.pSoluteUllCOMX[pos];
    PMEUllInt ullCOMY         = cSim.pSoluteUllCOMY[pos];
    PMEUllInt ullCOMZ         = cSim.pSoluteUllCOMZ[pos];
    cSim.pSoluteUllCOMX[pos]  = 0;
    cSim.pSoluteUllCOMY[pos]  = 0;
    cSim.pSoluteUllCOMZ[pos]  = 0;
    invMass                  *= ONEOVERENERGYSCALE;
    PMEDouble CX, CY, CZ;
    if (ullCOMX >= 0x8000000000000000ull) {
      CX = -(PMEDouble)(ullCOMX ^ 0xffffffffffffffffull);
    }
    else {
      CX = (PMEDouble)ullCOMX;
    }
    cSim.pSoluteCOMX[pos] = invMass * CX;
    if (ullCOMY >= 0x8000000000000000ull) {
      CY = -(PMEDouble)(ullCOMY ^ 0xffffffffffffffffull);
    }
    else {
      CY = (PMEDouble)ullCOMY;
    }
    cSim.pSoluteCOMY[pos] = invMass * CY;
    if (ullCOMZ >= 0x8000000000000000ull) {
      CZ = -(PMEDouble)(ullCOMZ ^ 0xffffffffffffffffull);
    }
    else {
      CZ = (PMEDouble)ullCOMZ;
    }
    cSim.pSoluteCOMZ[pos] = invMass * CZ;
    pos += blockDim.x * gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kReduceSoluteCOM: host function to launch the kernel above.  It is called in one instance
//                   in gpu_ntp_setup_() as part of the constant pressure initialization.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kReduceSoluteCOM(gpuContext gpu)
{
  kReduceSoluteCOM_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
  LAUNCHERROR("kReduceSoluteCOM");
}

//---------------------------------------------------------------------------------------------
// kClearSoluteCOM_kernel: clear, or initialize, the positions of the centers of mass for all
///                        molecules in the system.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(THREADS_PER_BLOCK, 1)
kClearSoluteCOM_kernel()
{
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  while (pos < cSim.soluteMolecules) {
    cSim.pSoluteUllCOMX[pos] = 0;
    cSim.pSoluteUllCOMY[pos] = 0;
    cSim.pSoluteUllCOMZ[pos] = 0;
    pos += blockDim.x * gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kClearSoluteCOM: host function to launch the kernel above.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kClearSoluteCOM(gpuContext gpu)
{
  kClearSoluteCOM_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
  LAUNCHERROR("kClearSoluteCOM");
}

//---------------------------------------------------------------------------------------------
// kCalculateSoluteCOM: compute the solute center of mass.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kCalculateSoluteCOM(gpuContext gpu)
{
  if (gpu->sim.soluteMolecules <= gpu->maxSoluteMolecules) {
    kPMECalculateSoluteCOM_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
  }
  else {
    kPMECalculateSoluteCOMLarge_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
  }
  LAUNCHERROR("kCalculateSoluteCOM");
}

//---------------------------------------------------------------------------------------------
// kCalculateCOMKineticEnergy: compute the kinetic energy of the solute center of mass.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kCalculateCOMKineticEnergy(gpuContext gpu)
{
  if (gpu->sim.soluteMolecules <= gpu->maxSoluteMolecules) {
    kPMECalculateCOMKineticEnergy_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
  }
  else {
    kPMECalculateCOMKineticEnergyLarge_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
  }
  LAUNCHERROR("kCalculateCOMKineticEnergy");
}

//---------------------------------------------------------------------------------------------
// kReduceCOMKineticEnergy_kernel: kernel for computing the center of mass kinetic energy.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(THREADS_PER_BLOCK, 1)
kReduceCOMKineticEnergy_kernel()
{
  struct COMKineticEnergy {
    double EKCOMX;
    double EKCOMY;
    double EKCOMZ;
  };

  __shared__ volatile COMKineticEnergy sE[THREADS_PER_BLOCK];

  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  PMEDouble EKCOMX = (PMEDouble)0.0;
  PMEDouble EKCOMY = (PMEDouble)0.0;
  PMEDouble EKCOMZ = (PMEDouble)0.0;
  while (pos < cSim.soluteMolecules) {
    PMEDouble invMass   = cSim.pSoluteInvMass[pos];
    PMEUllInt ullEKCOMX = cSim.pSoluteUllEKCOMX[pos];
    PMEUllInt ullEKCOMY = cSim.pSoluteUllEKCOMY[pos];
    PMEUllInt ullEKCOMZ = cSim.pSoluteUllEKCOMZ[pos];
    cSim.pSoluteUllEKCOMX[pos] = 0;
    cSim.pSoluteUllEKCOMY[pos] = 0;
    cSim.pSoluteUllEKCOMZ[pos] = 0;
    invMass *= ONEOVERFORCESCALESQUARED;
    PMEDouble EX, EY, EZ;
    if (ullEKCOMX >= 0x8000000000000000ull) {
      EX = -(PMEDouble)(ullEKCOMX ^ 0xffffffffffffffffull);
    }
    else {
      EX = (PMEDouble)ullEKCOMX;
    }
    EKCOMX += invMass * EX * EX;
    if (ullEKCOMY >= 0x8000000000000000ull) {
      EY = -(PMEDouble)(ullEKCOMY ^ 0xffffffffffffffffull);
    }
    else {
      EY = (PMEDouble)ullEKCOMY;
    }
    EKCOMY += invMass * EY * EY;
    if (ullEKCOMZ >= 0x8000000000000000ull) {
      EZ = -(PMEDouble)(ullEKCOMZ ^ 0xffffffffffffffffull);
    }
    else {
      EZ = (PMEDouble)ullEKCOMZ;
    }
    EKCOMZ += invMass * EZ * EZ;
    pos    += blockDim.x * gridDim.x;
  }
  sE[threadIdx.x].EKCOMX = EKCOMX;
  sE[threadIdx.x].EKCOMY = EKCOMY;
  sE[threadIdx.x].EKCOMZ = EKCOMZ;
  __syncthreads();
  unsigned int m = 1;
  while (m < blockDim.x) {
    int p = threadIdx.x + m;

    // OPTIMIZE
    PMEDouble EX = ((p < blockDim.x) ? sE[p].EKCOMX : (PMEDouble)0.0);
    PMEDouble EY = ((p < blockDim.x) ? sE[p].EKCOMY : (PMEDouble)0.0);
    PMEDouble EZ = ((p < blockDim.x) ? sE[p].EKCOMZ : (PMEDouble)0.0);
    __syncthreads();
    sE[threadIdx.x].EKCOMX += EX;
    sE[threadIdx.x].EKCOMY += EY;
    sE[threadIdx.x].EKCOMZ += EZ;
    __syncthreads();
    m *= 2;
  }
  if (threadIdx.x == 0) {
    unsigned long long int val1 = llitoulli(lrint(sE[0].EKCOMX * FORCESCALE));
    unsigned long long int val2 = llitoulli(lrint(sE[0].EKCOMY * FORCESCALE));
    unsigned long long int val3 = llitoulli(lrint(sE[0].EKCOMZ * FORCESCALE));
    if (val1 != 0) {
      atomicAdd(cSim.pEKCOMX, val1);
    }
    if (val2 != 0) {
      atomicAdd(cSim.pEKCOMY, val2);
    }
    if (val3 != 0) {
      atomicAdd(cSim.pEKCOMZ, val3);
    }
  }
}

//---------------------------------------------------------------------------------------------
// kReduceCOMKineticEnergy: host function to launch the above kernel.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kReduceCOMKineticEnergy(gpuContext gpu)
{
  kReduceCOMKineticEnergy_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
  LAUNCHERROR("kReduceCOMKineticEnergy");
}

//---------------------------------------------------------------------------------------------
// kCalculateMolecularVirial: calculate the molecular virial, the standard used by AMBER in
//                            determining the system pressure (as opposed to an atomic virial)
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kCalculateMolecularVirial(gpuContext gpu)
{
  if (gpu->sim.soluteMolecules <= gpu->maxSoluteMolecules) {
    kCalculateMolecularVirial_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
  }
  else {
    kCalculateMolecularVirialLarge_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
  }
  LAUNCHERROR("kCalculateMolecularVirial");
}

//---------------------------------------------------------------------------------------------
// kPressureScaleCoordinates: scale coordinates to accommodate box size changes in constant
//                            pressure simulations.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kPressureScaleCoordinates(gpuContext gpu)
{
  if (gpu->sim.soluteMolecules <= gpu->maxPSSoluteMolecules) {
    kPressureScaleCoordinates_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
  }
  else {
    kPressureScaleCoordinatesLarge_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
  }
  LAUNCHERROR("kPressureScaleCoordinates");
}

//---------------------------------------------------------------------------------------------
// kPressureScaleConstraintCoordinates_kernel: scale the locations of positional restraints
//                                             to track the simulation box expansion and
//                                             contraction in constant pressure simulations.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(THREADS_PER_BLOCK, 1)
kPressureScaleConstraintCoordinates_kernel()
{
  __shared__ PMEDouble sUcell[9];

  // Read transformation matrices
  if (threadIdx.x < 9) {
    sUcell[threadIdx.x] = cSim.pNTPData->ucell[threadIdx.x];
  }
  __syncthreads();

  // Iterate over constraints
  unsigned int pos = blockIdx.x*blockDim.x + threadIdx.x;
  while (pos < cSim.constraints) {

    // Read the COM.  The COM is pre-scaled into fractional coordinates so as to always
    // start with fresh data and avoid heinous FPRE over long simulations in SPFP and
    // SPDP precision modes.  The original input coordinates of any restraints must
    // also always be scaled, in particular if the coordinates were initially centered,
    // so that the system does not gradually shrink towards the origin and subtly
    // increase pressure over 1M iterations or so.
    PMEDouble fx = cSim.pConstraintCOMX[pos];
    PMEDouble fy = cSim.pConstraintCOMY[pos];
    PMEDouble fz = cSim.pConstraintCOMZ[pos];
    PMEDouble x  = cSim.pConstraintAtomX[pos];
    PMEDouble y  = cSim.pConstraintAtomY[pos];
    PMEDouble z  = cSim.pConstraintAtomZ[pos];
    PMEDouble2 constraint1 = cSim.pConstraint1[pos];

    // Calculate COM displacement
    PMEDouble cx = fx*sUcell[0] + fy*sUcell[1] + fz*sUcell[2];
    PMEDouble cy =                fy*sUcell[4] + fz*sUcell[5];
    PMEDouble cz =                               fz*sUcell[8];
    x += cx;
    y += cy;
    z += cz;
    PMEDouble2 constraint2 = {y, z};
    constraint1.y          = x;
    cSim.pConstraint2[pos] = constraint2;
    cSim.pConstraint1[pos] = constraint1;

    // Update the bond work units data array
    unsigned int iupdate = cSim.pBwuCnstUpdateIdx[pos];
    cSim.pBwuCnst[iupdate       ] = constraint1;
    cSim.pBwuCnst[iupdate + GRID] = constraint2;

    // Increment the counter
    pos += blockDim.x * gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kPressureScaleConstraintCoordinates: host function to launch the above kernel.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void kPressureScaleConstraintCoordinates(gpuContext gpu)
{
  kPressureScaleConstraintCoordinates_kernel<<<gpu->blocks, gpu->threadsPerBlock>>>();
  LAUNCHERROR("kPressureScaleConstraintCoordinates");
}
