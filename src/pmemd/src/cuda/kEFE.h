#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// This file is included by kCalculateEFieldEnergy.cu for computing PME energies in the context
// of an external electric field.
//
// #defines: PME_ENERGY, NEIGHBOR_LIST
//---------------------------------------------------------------------------------------------
{
  // Precompute Electric Field Constants
  PMEFloat phase = cos((2*PI*cSim.effreq/1000) * ((PMEFloat)dt * (PMEFloat)nstep) -
                       (PI / 180) * cSim.efphase);
  PMEFloat loc_efx = phase * (PMEFloat)cSim.efx;
  PMEFloat loc_efy = phase * (PMEFloat)cSim.efy;
  PMEFloat loc_efz = phase * (PMEFloat)cSim.efz;
#ifdef PME_ENERGY
  PMEForce sEEField = (PMEForce)0;
#endif
  unsigned int pos       = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int increment = gridDim.x * blockDim.x;
  unsigned int imgPos    = cSim.pImageAtomLookup[pos];
  if (cSim.efn == 1) {
    loc_efx *= (PMEFloat)cSim.recip[0][0];
    loc_efy *= (PMEFloat)cSim.recip[1][1];
    loc_efz *= (PMEFloat)cSim.recip[2][2];
  }
  if (pos < cSim.atoms) {

    //Convert internal charge to electron charge
    PMEDouble electron_charge = (PMEDouble)cSim.pImageCharge[imgPos] / (PMEDouble)18.2223;
    PMEDouble ef_frcx = electron_charge * loc_efx;
    PMEDouble ef_frcy = electron_charge * loc_efy;
    PMEDouble ef_frcz = electron_charge * loc_efz;

#ifdef use_SPFP
    atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[imgPos],
              llitoulli(ef_frcx * FORCESCALEF));
    atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[imgPos],
              llitoulli(ef_frcy * FORCESCALEF));
    atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[imgPos],
              llitoulli(ef_frcz * FORCESCALEF));
#else  // use_DPFP
    atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[imgPos],
              llitoulli(llrint((PMEForce)ef_frcx * FORCESCALE)));
    atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[imgPos],
              llitoulli(llrint((PMEForce)ef_frcy * FORCESCALE)));
    atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[imgPos],
              llitoulli(llrint((PMEForce)ef_frcz * FORCESCALE)));
#endif // End pre-processor branch over precision modes

#ifdef PME_ENERGY
    PMEDouble AtomX = (PMEDouble)cSim.pImageX[imgPos];
    PMEDouble AtomY = (PMEDouble)cSim.pImageY[imgPos];
    PMEDouble AtomZ = (PMEDouble)cSim.pImageZ[imgPos];
    PMEDouble ef_vx = AtomX - (PMEDouble)cSim.ucellf[0][0];
    PMEDouble ef_vy = AtomY - (PMEDouble)cSim.ucellf[1][1];
    PMEDouble ef_vz = AtomZ - (PMEDouble)cSim.ucellf[2][2];
#  ifndef use_DPFP
    sEEField       -= fast_llrintf(ENERGYSCALEF*(PMEFloat)(ef_vx*ef_frcx + ef_vy*ef_frcy +
                                                           ef_vz*ef_frcz));
#  else
    sEEField       -= (ef_vx*(PMEDouble)ef_frcx + ef_vy*(PMEDouble)ef_frcy +
                       ef_vz*(PMEDouble)ef_frcz);
#  endif
#endif
    pos            += increment;
#ifdef PME_ENERGY
#  ifndef use_DPFP
    atomicAdd(cSim.pEEField, llitoulli(sEEField));
#  else
    atomicAdd(cSim.pEEField, llitoulli(llrint(sEEField * ENERGYSCALE)));
#  endif
#endif
  }
}
