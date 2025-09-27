#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------
{
#ifdef AFE_NEIGHBORLIST
#  ifdef AFE_Vel
#    define VELX cSim.pImageVelX
#    define VELY cSim.pImageVelY
#    define VELZ cSim.pImageVelZ
#  endif

#  ifdef AFE_Crd
#    define ATOMX cSim.pImageX
#    define ATOMY cSim.pImageY
#    define ATOMZ cSim.pImageZ
#  endif

#else  // AFE_NEIGHBORLIST
#  ifdef AFE_Vel
#    define VELX cSim.pVelX
#    define VELY cSim.pVelY
#    define VELZ cSim.pVelZ
#  endif

#  ifdef AFE_Crd
#    define ATOMX cSim.pAtomX
#    define ATOMY cSim.pAtomY
#    define ATOMZ cSim.pAtomZ
#  endif
#endif // AFE_NEIGHBORLIST

  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

  // The number of linear atoms in region 1 will always
  // equal the number of linear atoms in region 2
  if (pos < cSim.TIlinearAtmCnt) {

    // Set indexes to be the index of the linear atoms w.r.t. the total atom list
    unsigned int index_atm1 = cSim.pImageTILinearAtmID[pos];
    unsigned int index_atm2 = cSim.pImageTILinearAtmID[pos+cSim.TIPaddedLinearAtmCnt];
#ifdef AFE_Frc
    PMEAccumulator fx_1 = cSim.pForceXAccumulator[index_atm1];
    PMEAccumulator fy_1 = cSim.pForceYAccumulator[index_atm1];
    PMEAccumulator fz_1 = cSim.pForceZAccumulator[index_atm1];
    PMEAccumulator fx_2 = cSim.pForceXAccumulator[index_atm2];
    PMEAccumulator fy_2 = cSim.pForceYAccumulator[index_atm2];
    PMEAccumulator fz_2 = cSim.pForceZAccumulator[index_atm2];

    // For forces, combine the forces from region 1 with forces from region 2.
    atomicAdd((unsigned long long int*)&cSim.pForceXAccumulator[index_atm1], fx_2);
    atomicAdd((unsigned long long int*)&cSim.pForceXAccumulator[index_atm2], fx_1);
    atomicAdd((unsigned long long int*)&cSim.pForceYAccumulator[index_atm1], fy_2);
    atomicAdd((unsigned long long int*)&cSim.pForceYAccumulator[index_atm2], fy_1);
    atomicAdd((unsigned long long int*)&cSim.pForceZAccumulator[index_atm1], fz_2);
    atomicAdd((unsigned long long int*)&cSim.pForceZAccumulator[index_atm2], fz_1);
#endif // AFE_Frc

#ifdef AFE_Vel
    VELX[index_atm2] = VELX[index_atm1];
    VELY[index_atm2] = VELY[index_atm1];
    VELZ[index_atm2] = VELZ[index_atm1];
#endif // AFE_Vel

#ifdef AFE_Crd
    ATOMX[index_atm2] = ATOMX[index_atm1];
    ATOMY[index_atm2] = ATOMY[index_atm1];
    ATOMZ[index_atm2] = ATOMZ[index_atm1];
#endif // AFE_Crd
  }
#ifdef AFE_Vel
#  undef VELX
#  undef VELY
#  undef VELZ
#endif
#ifdef AFE_Crd
#  undef ATOMX
#  undef ATOMY
#  undef ATOMZ
#endif
}
