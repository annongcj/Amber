#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// This is included by kCalculateLocalForces.cu for different situations.
//---------------------------------------------------------------------------------------------
{
#if !defined(AMBER_PLATFORM_AMD)
#  define VOLATILE volatile
#else
#  define VOLATILE
#endif
  __shared__ volatile PMEDouble atmcrdx[BOND_WORK_UNIT_THREADS_PER_BLOCK];
  __shared__ volatile PMEDouble atmcrdy[BOND_WORK_UNIT_THREADS_PER_BLOCK];
  __shared__ volatile PMEDouble atmcrdz[BOND_WORK_UNIT_THREADS_PER_BLOCK];
#ifdef PHMD
  __shared__ volatile PMEFloat charges[BOND_WORK_UNIT_THREADS_PER_BLOCK];
  __shared__ volatile int indices[BOND_WORK_UNIT_THREADS_PER_BLOCK];
#endif
#ifdef use_DPFP
  __shared__ PMEAccumulator atmfrcx[BOND_WORK_UNIT_THREADS_PER_BLOCK];
  __shared__ PMEAccumulator atmfrcy[BOND_WORK_UNIT_THREADS_PER_BLOCK];
  __shared__ PMEAccumulator atmfrcz[BOND_WORK_UNIT_THREADS_PER_BLOCK];
  #ifdef PHMD
  __shared__ PMEAccumulator atmdudl[BOND_WORK_UNIT_THREADS_PER_BLOCK];
  __shared__ PMEAccumulator atmdudlplus[BOND_WORK_UNIT_THREADS_PER_BLOCK];
  #endif
#else
  __shared__ int atmfrcx[BOND_WORK_UNIT_THREADS_PER_BLOCK];
  __shared__ int atmfrcy[BOND_WORK_UNIT_THREADS_PER_BLOCK];
  __shared__ int atmfrcz[BOND_WORK_UNIT_THREADS_PER_BLOCK];
  #ifdef PHMD
  __shared__ int atmdudl[BOND_WORK_UNIT_THREADS_PER_BLOCK];
  __shared__ int atmdudlplus[BOND_WORK_UNIT_THREADS_PER_BLOCK];
  #endif
#endif
  __shared__ VOLATILE int insrIdx, wuidx, curr_wuidx;
  __shared__ VOLATILE int insrWarpIdx[BOND_WORK_UNIT_WARPS_PER_BLOCK];
  __shared__ VOLATILE unsigned int instructions[BOND_WORK_UNIT_THREADS_PER_BLOCK];
#ifdef LOCAL_ENERGY
  __shared__ VOLATILE PMEAccumulator nrgACC[EACC_TOTAL_SIZE];
#  if defined(LOCAL_AFE) && defined(LOCAL_MBAR)
  __shared__ PMEAccumulator sBarTot[BOND_WORK_UNIT_THREADS_PER_BLOCK];
#  endif
#endif
  // Set reference indices.  Initialize __shared__ memory accumulators.
  int warpIdx = threadIdx.x >> GRID_BITS;
 unsigned int tgx = threadIdx.x & GRID_BITS_MASK;
  if (threadIdx.x == 0) {
    wuidx = blockIdx.x;
  }
  __syncthreads();
  
  // The kernel launch will handle clusters of consecutive bond work units in each thread
  // block.  For most systems the cluster size will simply be 1, but for very large systems
  // with more than 65536 bond work units the cluster size will increase.
#ifdef AMBER_PLATFORM_AMD
  if (wuidx < cSim.bondWorkUnits) {
#else
  while (wuidx < cSim.bondWorkUnits) {
#endif

    // Read the instructions into __shared__ memory, where they will be kept resident
    // to capitalize on the fast operation of broadcasting one element of __shared__
    // to an entire warp.
    int readpos = (3 * BOND_WORK_UNIT_THREADS_PER_BLOCK * wuidx) + threadIdx.x;
    instructions[threadIdx.x] = cSim.pBwuInstructions[readpos];
    __syncthreads();

    // Read coordinates
    if (threadIdx.x < instructions[ATOM_IMPORT_COUNT_IDX]) {
#ifdef LOCAL_NEIGHBORLIST
      readpos = cSim.pBwuInstructions[readpos + 2*BOND_WORK_UNIT_THREADS_PER_BLOCK];
      atmcrdx[threadIdx.x] = cSim.pImageX[readpos];
      atmcrdy[threadIdx.x] = cSim.pImageY[readpos]; 
      atmcrdz[threadIdx.x] = cSim.pImageZ[readpos]; 
#ifdef PHMD
      indices[threadIdx.x] = readpos;
      charges[threadIdx.x] = cSim.pImageCharge[readpos];
#endif
#else
      readpos = cSim.pBwuInstructions[readpos + BOND_WORK_UNIT_THREADS_PER_BLOCK];
      atmcrdx[threadIdx.x] = cSim.pAtomX[readpos];
      atmcrdy[threadIdx.x] = cSim.pAtomY[readpos];
      atmcrdz[threadIdx.x] = cSim.pAtomZ[readpos];
#ifdef PHMD
      indices[threadIdx.x] = readpos;
      charges[threadIdx.x] = cSim.pAtomChargeSP[readpos];
#endif
#endif
    }

    // Initialize __shared__ memory force accumulators
#ifdef use_DPFP
    atmfrcx[threadIdx.x] = (PMEAccumulator)0;
    atmfrcy[threadIdx.x] = (PMEAccumulator)0;
    atmfrcz[threadIdx.x] = (PMEAccumulator)0;
  #ifdef PHMD
    atmdudl[threadIdx.x] = (PMEAccumulator)0;
    atmdudlplus[threadIdx.x] = (PMEAccumulator)0;
  #endif
#else
    atmfrcx[threadIdx.x] = 0;
    atmfrcy[threadIdx.x] = 0;
    atmfrcz[threadIdx.x] = 0;
  #ifdef PHMD
    atmdudl[threadIdx.x] = 0;
    atmdudlplus[threadIdx.x] = 0;
  #endif
#endif
#ifdef LOCAL_ENERGY
    if (threadIdx.x < EACC_TOTAL_SIZE) {
      nrgACC[threadIdx.x] = (PMEAccumulator)0;
    }
#  if defined(LOCAL_AFE) && defined(LOCAL_MBAR)
    sBarTot[threadIdx.x] = (PMEAccumulator)0;
#  endif
#endif

    // Set the master warp instruction counter for the whole block
    if (threadIdx.x == 0) {
      insrIdx = BOND_WORK_UNIT_THREADS_PER_BLOCK >> GRID_BITS;

      // Increment the counter to access the next work unit.  This can happen
      // between this pair of __syncthreads() calls as there is no other
      // dependency on wuidx.  The work in the following sections before the
      // big loop's final __syncthreads() call will rely on curr_wuidx.
      curr_wuidx = wuidx;
#ifndef AMBER_PLATFORM_AMD
      wuidx = atomicAdd(&cSim.pFrcBlkCounters[1], 1);
#endif
    }
    __syncthreads();

    // Start executing warp instructions.  Each warp starts on its natural index, but
    // thereafter takes the next available warp instruction.
    unsigned int iWarpIdx = warpIdx;
#ifdef AMBER_PLATFORM_AMD
    while (iWarpIdx < instructions[WARP_INSTRUCTION_COUNT_IDX]) {
#else
    if (tgx == 0) {
      insrWarpIdx[warpIdx] = warpIdx;
    }
    __SYNCWARP(WARP_MASK);
    while (insrWarpIdx[warpIdx] < instructions[WARP_INSTRUCTION_COUNT_IDX]) {
      iWarpIdx = insrWarpIdx[warpIdx];
#endif

      // Extract the directions and index into the atom index / parameter arrays.
      // In what follows, atmI, atmJ, ... atmM are all indices into the array of
      // imported atoms in the range [0, BOND_WORK_UNIT_THREADS_PER_BLOCK).
      unsigned int myinsr = instructions[WARP_INSTRUCTION_OFFSET + iWarpIdx];
      int action = myinsr >> 24;
      int startidx = (myinsr & 0xffffff) * GRID;

      // Lay out some variables that will be needed by most any interaction,
      // with pre-processor contingencies for TI and MBAR.
#ifdef LOCAL_AFE
      int TIregion, CVterm;
#  ifdef LOCAL_ENERGY
      int SCterm = 0;
      PMEAccumulator scEcomp[2], edvdl;
      scEcomp[0] = (PMEAccumulator)0;
      scEcomp[1] = (PMEAccumulator)0;
      edvdl      = (PMEAccumulator)0;
      TIregion   = 0;
      CVterm     = 0;
#    ifdef LOCAL_MBAR
      PMEDouble mbarTerm, mbarRefLambda;
#    endif
#  endif
#endif
#ifdef LOCAL_ENERGY
      PMEAccumulator eterm;
#endif

      // Case switch over all interactions
#if defined(AMD_SUM)
      if (action == DIHE_CODE) {
#include "kBWU_dihe.h"
      }
#elif defined(GAMD_SUM)
      if (action == DIHE_CODE) {
#include "kBWU_dihe.h"
      }
#else  // Pre-processor branch over aMD energy sums, GaMD energy sums, or full evaluations
      if (action == BOND_CODE) {
#include "kBWU_bond.h"
      }
      else if (action == ANGL_CODE) {
#include "kBWU_angl.h"
      }
      else if (action == DIHE_CODE) {
#include "kBWU_dihe.h"
      }
      else if (action == CMAP_CODE) {
#include "kBWU_cmap.h"
      }
      else if (action == NB14_CODE) {
#include "kBWU_nb14.h"
      }
      else if (action == NMR2_CODE) {
#include "kBWU_nmr2.h"
      }
      else if (action == NMR3_CODE) {
#include "kBWU_nmr3.h"
      }
      else if (action == NMR4_CODE) {
#include "kBWU_nmr4.h"
      }
      else if (action == UREY_CODE) {
#include "kBWU_urey.h"
      }
      else if (action == CIMP_CODE) {
#include "kBWU_cimp.h"
      }
      else if (action == CNST_CODE) {
#include "kBWU_cnst.h"
      }
#endif // Pre-processor branch over aMD energy sums, GaMD energy sums, or full evaluations
#ifdef AMBER_PLATFORM_AMD
      iWarpIdx += (BOND_WORK_UNIT_THREADS_PER_BLOCK >> GRID_BITS);
#else
      // Take the next available warp instruction
      __SYNCWARP(WARP_MASK);
      if (tgx == 0) {
        insrWarpIdx[warpIdx] = atomicAdd((int*)&insrIdx, 1);
      }
      __SYNCWARP(WARP_MASK);
#endif
    }
    __syncthreads();

    // Write forces back to global accumulators
    if (threadIdx.x < instructions[ATOM_IMPORT_COUNT_IDX]) {
#ifdef LOCAL_NEIGHBORLIST
      int atmX = cSim.pBwuInstructions[((3*curr_wuidx + 2)*BOND_WORK_UNIT_THREADS_PER_BLOCK) +
                                       threadIdx.x];
#else
      int atmX = cSim.pBwuInstructions[((3*curr_wuidx + 1)*BOND_WORK_UNIT_THREADS_PER_BLOCK) +
                                       threadIdx.x];
#endif
#ifdef use_DPFP
      PMEAccumulator ifx = atmfrcx[threadIdx.x];
      PMEAccumulator ify = atmfrcy[threadIdx.x];
      PMEAccumulator ifz = atmfrcz[threadIdx.x];
  #ifdef PHMD
      PMEAccumulator dudli = atmdudl[threadIdx.x];
      PMEAccumulator dudliplus = atmdudlplus[threadIdx.x];
  #endif
#else
      PMEAccumulator ifx = (PMEAccumulator)(atmfrcx[threadIdx.x]) * BWU_PROMOTION_FACTOR;
      PMEAccumulator ify = (PMEAccumulator)(atmfrcy[threadIdx.x]) * BWU_PROMOTION_FACTOR;
      PMEAccumulator ifz = (PMEAccumulator)(atmfrcz[threadIdx.x]) * BWU_PROMOTION_FACTOR;
  #ifdef PHMD
      PMEAccumulator dudli = (PMEAccumulator)(atmdudl[threadIdx.x]) * BWU_PROMOTION_FACTOR;
      PMEAccumulator dudliplus = (PMEAccumulator)(atmdudlplus[threadIdx.x]) * BWU_PROMOTION_FACTOR;
  #endif
#endif
      if (instructions[FORCE_ACCUMULATOR_IDX] == BOND_FORCE_ACCUMULATOR) {
        atomicAdd((unsigned long long int*)&cSim.pBondedForceXAccumulator[atmX],
                  llitoulli(ifx));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceYAccumulator[atmX],
                  llitoulli(ify));
        atomicAdd((unsigned long long int*)&cSim.pBondedForceZAccumulator[atmX],
                  llitoulli(ifz));
      }
      else {
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[atmX],
                  llitoulli(ifx));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[atmX],
                  llitoulli(ify));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[atmX],
                  llitoulli(ifz));
      }
  #ifdef PHMD
      int h = cSim.pImageGrplist[atmX] - 1;
      if(h >= 0) {
        atomicAdd((unsigned long long int*)&cSim.pdph_accumulator[atmX],
                   llitoulli(dudli));
        if(cSim.psp_grp[h] > 0) {
          atomicAdd((unsigned long long int*)&cSim.pdph_plus_accumulator[atmX],
                     llitoulli(dudliplus));
        }
      }
  #endif
    }
    __syncthreads();
#ifdef LOCAL_ENERGY
    // Reduce the energy accumulators.  Would be rather nice to have 'em all in one array
    // up in cudaSimulation, too, but that might be a long thread to pull on and not for
    // much optimization at all.  This kernel must be launched with at least 128 threads
    // due to the code below, however, spreading the global atomics over multiple warps.
    if ((threadIdx.x < EACC_TOTAL_SIZE) &&
        (threadIdx.x % BOND_WORK_UNIT_WARPS_PER_BLOCK == 0)) {
      for (int i = 1; i < BOND_WORK_UNIT_WARPS_PER_BLOCK; i++) {
        nrgACC[threadIdx.x] += nrgACC[threadIdx.x + i];
      }
    }
    __syncthreads();
#  if defined(AMD_SUM)
    if (threadIdx.x == 0) {
      atomicAdd(cSim.pAMDEDihedral, llitoulli(nrgACC[DIHE_EACC_OFFSET]));
    }
#  elif defined(GAMD_SUM)
    if (threadIdx.x == 0) {
      atomicAdd(cSim.pGaMDEDihedral, llitoulli(nrgACC[DIHE_EACC_OFFSET]));
    }
#  else
    if (threadIdx.x == 0) {
      atomicAdd(cSim.pEBond, llitoulli(nrgACC[BOND_EACC_OFFSET]));
      atomicAdd(cSim.pEAngle, llitoulli(nrgACC[ANGL_EACC_OFFSET]));
      atomicAdd(cSim.pEDihedral, llitoulli(nrgACC[DIHE_EACC_OFFSET]));
      atomicAdd(cSim.pEED, llitoulli(nrgACC[QQXC_EACC_OFFSET]));
    }
    else if (threadIdx.x == GRID) {
      atomicAdd(cSim.pENB14, llitoulli(nrgACC[SCNB_EACC_OFFSET]));
      atomicAdd(cSim.pEEL14, llitoulli(nrgACC[SCEE_EACC_OFFSET]));
      atomicAdd(cSim.pEConstraint, llitoulli(nrgACC[CNST_EACC_OFFSET]));
    }
    else if (threadIdx.x == 2*GRID) {
      atomicAdd(cSim.pEAngle_UB, llitoulli(nrgACC[UREY_EACC_OFFSET]));
      atomicAdd(cSim.pEImp, llitoulli(nrgACC[CIMP_EACC_OFFSET]));
      atomicAdd(cSim.pECmap, llitoulli(nrgACC[CMAP_EACC_OFFSET]));
    }
    else if (threadIdx.x == 3*GRID) {
      atomicAdd(cSim.pENMRDistance, llitoulli(nrgACC[NMR2_EACC_OFFSET]));
      atomicAdd(cSim.pENMRAngle, llitoulli(nrgACC[NMR3_EACC_OFFSET]));
      atomicAdd(cSim.pENMRTorsion, llitoulli(nrgACC[NMR4_EACC_OFFSET]));
    }
#  endif
#  if defined(LOCAL_AFE) && defined(LOCAL_MBAR)
    if (threadIdx.x < cSim.bar_states) {
      atomicAdd(&cSim.pBarTot[threadIdx.x], llitoulli(sBarTot[threadIdx.x]));
    }
#  endif
#endif
  }

  // Continuing on, clear the charge grid buffer.  This work is independent of
  // the bond work units and is therefore an ideal way to backfill idle threads.
  // However, it needn't be done twice for the Accelerated MD modes.
#ifndef AMBER_PLATFORM_AMD
#if !defined(AMD_SUM) && !defined(GAMD_SUM)
  while (wuidx < cSim.clearQBWorkUnits) {
    int pos = CHARGE_BUFFER_STRIDE*(wuidx - cSim.bondWorkUnits) + threadIdx.x;
    int i;
    for (i = 0; i < CHARGE_BUFFER_CYCLES; i++) {
      if (pos < cSim.XYZStride) {
#  ifdef use_DPFP
        cSim.plliXYZ_q[pos] = (long long int)0;
#  else
        cSim.plliXYZ_q[pos] = 0;
#  endif
      }
      pos += blockDim.x;
    }
    __syncthreads();
    if (threadIdx.x == 0) {
      wuidx = atomicAdd(&cSim.pFrcBlkCounters[1], 1);
    }
    __syncthreads();
  }
#endif
#endif
#undef VOLATILE
}
