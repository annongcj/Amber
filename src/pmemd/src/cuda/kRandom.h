#include "copyright.i"

#include "hip_definitions.h" 

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------
static double u[98];
static double c = 0.0;
static double cd = 0.0;
static double cm = 0.0;
static int i97 = 0;
static int j97 = 0;

//---------------------------------------------------------------------------------------------
// cpu_amrset: set the random number generator according to the CPU
//
// Arguments:
//   iseed:    pseudo-random number generator seed
//---------------------------------------------------------------------------------------------
void cpu_amrset(int iseed)
{
  int is1;
  int is2;
  int is1max = 31328;
  int is2max = 30081;
  int i, ii, j, jj, k, l, m;
  double s, t;

  is1 = max((iseed / is2max) + 1, 1);
  is1 = min(is1, is1max);
  is2 = max(1, (iseed % is2max) + 1);
  is2 = min(is2, is2max);
  i = ((is1/177) % 177) + 2;
  j = ((is1)     % 177) + 2;
  k = ((is2/169) % 178) + 1;
  l = ((is2)     % 169);
  for (ii = 1; ii <= 97; ii++) {
    s = 0.00;
    t = 0.50;
    for (jj = 1; jj <= 24; jj++) {
      m = (((i * j) % 179) * k) % 179;
      i = j;
      j = k;
      k = m;
      l = (53 * l + 1) % 169;
      if (((l * m) % 64) >= 32) {
        s = s + t;
      }
      t = 0.50 * t;
    }
    u[ii] = s;
  }
  c  = 362436.0   / 16777216.0;
  cd = 7654321.0  / 16777216.0;
  cm = 16777213.0 / 16777216.0;
  i97 = 97;
  j97 = 33;

  return;
}

//---------------------------------------------------------------------------------------------
// cpu_gauss: CPU-based Gaussian random number generator
//
// Arguments:
//   am, sd: parameters for the distribution
//---------------------------------------------------------------------------------------------
double cpu_gauss(double am, double sd)
{
  double tmp1, tmp2;
  double uni;
  double zeta1, zeta2;

  while (1) {
    uni = u[i97] - u[j97];
    if (uni < 0.0) {
      uni = uni + 1.0;
    }
    u[i97] = uni;
    i97 = i97 - 1;
    if (i97 == 0) {
      i97 = 97;
    }
    j97 = j97 - 1;
    if (j97 == 0) {
      j97 = 97;
    }
    c = c - cd;
    if (c < 0.0) {
      c = c + cm;
    }
    uni = uni - c;
    if (uni < 0.0) {
      uni = uni + 1.0;
    }
    zeta1 = uni + uni - 1.;
    uni = u[i97] - u[j97];
    if (uni < 0.0) {
      uni = uni + 1.0;
    }
    u[i97] = uni;
    i97 = i97 - 1;
    if (i97 == 0) {
      i97 = 97;
    }
    j97 = j97 - 1;
    if (j97 == 0) {
      j97 = 97;
    }
    c = c - cd;
    if (c < 0.0) {
      c = c + cm;
    }
    uni = uni - c;
    if (uni < 0.0) {
      uni = uni + 1.0;
    }
    zeta2 = uni + uni - 1.0;
    tmp1 = zeta1*zeta1 + zeta2*zeta2;
    if ((tmp1 < 1.0) && (tmp1 != 0.0)) {
      tmp2 = sd * sqrt(-2.0 * log(tmp1) / tmp1);
      return zeta1*tmp2 + am;
    }
  }
}

//---------------------------------------------------------------------------------------------
// cpu_kRandomPHMD: generate Gaussian random numbers on th eCPU, then upload them to the GPU.  This
//              is NOT efficient but useful for debugging purposes.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void cpu_kRandomPHMD(gpuContext gpu)
{
  for (int i = 0; i < gpu->sim.randomNumbersPHMD; i++) {
    gpu->pbRandomPHMD->_pSysData[i] = cpu_gauss(0.0, 1.0);
  }
  gpu->pbRandomPHMD->Upload();
}

//---------------------------------------------------------------------------------------------
// kRandomPHMD: launch the kernel to generate normally distributed (Gaussian) random numbers on
//          the GPU.  This is the efficient approach for production MD.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kRandomPHMD(gpuContext gpu)
{
    curandGenerateNormalDouble(gpu->RNG, gpu->pbRandomPHMD->_pDevData,
                               gpu->sim.randomNumbersPHMD, 0.0, 1.0);
}

//---------------------------------------------------------------------------------------------
// cpu_kRandom: generate Gaussian random numbers on th eCPU, then upload them to the GPU.  This
//              is NOT efficient but useful for debugging purposes.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
void cpu_kRandom(gpuContext gpu)
{
  for (int i = 0; i < gpu->sim.randomNumbers; i++) {
    gpu->pbRandom->_pSysData[i] = cpu_gauss(0.0, 1.0);
  }
  gpu->pbRandom->Upload();
}

//---------------------------------------------------------------------------------------------
// kRandom: launch the kernel to generate normally distributed (Gaussian) random numbers on
//          the GPU.  This is the efficient approach for production MD.
//
// Arguments:
//   gpu: overarching type for storing all parameters, coordinates, and the energy function
//---------------------------------------------------------------------------------------------
extern "C" void kRandom(gpuContext gpu)
{
  if ((gpu->ntt == 3) || (gpu->ntt == 2) || (gpu->ischeme == 1)) {
    curandGenerateNormalDouble(gpu->RNG, gpu->pbRandom->_pDevData,
                               gpu->sim.randomNumbers, 0.0, 1.0);
    LAUNCHERROR("kRandom");
  }
}

