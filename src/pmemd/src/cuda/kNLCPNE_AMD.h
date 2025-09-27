#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// This is included by kCalculatePMENonbondEnergy.cu, many times for different kernels
// calculating forces, energies, and virials, all with different numbers of atoms per warp.
// VDW force-switching code could use optimization.
//
// #defines: PME_ENERGY, PME_VIRIAL, PME_IS_ORTHOGONAL, PME_ATOMS_PER_WARP, PME_MINIMIZATION
//
// In this file, there is another added convention that nested pre-processor directives are
// also indented, as there are so many of them.
//---------------------------------------------------------------------------------------------
{
#define uint unsigned int
  struct NLAtom {
    PMEFloat x;
    PMEFloat y;
    PMEFloat z;
    PMEFloat q;
    unsigned int LJID;
    unsigned int ID;
  };

  struct NLWarp {
    NLEntry nlEntry;
    uint pos;
    bool bHomeCell;
  };

  const PMEFloat delta = 1.0e-5;
#ifdef PME_ENERGY
  const PMEFloat sixth = (PMEFloat)0.16666666666667;
  const PMEFloat twelvth = (PMEFloat)0.083333333333333;
#endif

#if defined(PME_VIRIAL) || defined(PME_ENERGY)
  const int THREADS_PER_BLOCK = PMENONBONDENERGY_THREADS_PER_BLOCK;
#else
  const int THREADS_PER_BLOCK = PMENONBONDFORCES_THREADS_PER_BLOCK;
#endif

  __shared__ unsigned int sNLEntries;
#if defined(PME_VIRIAL)
  __shared__ PMEFloat sUcellf[9];
#endif
  __shared__ NLWarp sWarp[THREADS_PER_BLOCK / GRID];

  // Read static data
  if (threadIdx.x == 0) {
    sNLEntries = *(cSim.pNLEntries);
  }

#ifdef PME_VIRIAL
  if (threadIdx.x < 9) {
    sUcellf[threadIdx.x] = cSim.pNTPData->ucellf[threadIdx.x];
  }
#endif
#ifdef PME_ENERGY
  PMEForceAccumulator eed, evdw;
  eed = (PMEForceAccumulator)0;
  evdw = (PMEForceAccumulator)0;
#  if defined(use_DPFP) && defined(PME_MINIMIZATION)
  long long int eede = 0;
  long long int evdwe = 0;
#  endif
#endif

#ifdef PME_VIRIAL
  PMEVirialAccumulator vir_11, vir_22, vir_33;
  vir_11 = (PMEVirialAccumulator)0;
  vir_22 = (PMEVirialAccumulator)0;
  vir_33 = (PMEVirialAccumulator)0;
#  if defined(use_DPFP) && defined(PME_MINIMIZATION)
  long long int vir_11E = 0;
  long long int vir_22E = 0;
  long long int vir_33E = 0;
#  endif
#endif

#if defined(use_DPFP) && defined(PME_MINIMIZATION)
  PMEFloat i;
#endif

  // Position in the warp and warp indexing into __shared__ arrays
  uint tgx = threadIdx.x & (GRID - 1);
  uint warp = THREADS_PER_BLOCK == GRID ? 0 : threadIdx.x / GRID;

  // The number of iterations required to compute all interactions in the first tile.
  // There are PME_ATOMS_PER_WARP / 2 * (PME_ATOMS_PER_WARP - 1) unique interactions: for every
  // two atoms i and j we compute i-j (upper triangle of a matrix), but do not compute
  // 1) j-i (lower triangle) since it's just a reversed i-j;
  // 2) i-i, j-j (diagonal), i.e. self-interactions.
  // * for 16 atoms per warp, the first tile has 120 ((1+15)/2*16) interactions to compute,
  //   which requires 4 iterations (32-lane warps) or 2 iterations (64-lane warps);
  // * for 32 atoms per warp, there are 496 ((1+31)/2*32) interactions,
  //   16 iterations (32-lane warps) or 8 iterations (64-lane warps).
  // * for 8 atoms per warp, there are only 28 ((1+7)/2*8) interactions, 1 iteration;
  uint totalInteractions = PME_ATOMS_PER_WARP / 2 * (PME_ATOMS_PER_WARP - 1);
  uint iftlim = (totalInteractions + GRID - 1) / GRID;

  uint jgroup = tgx >> cSim.NLAtomsPerWarpBits;
  uint joffset = tgx & ~cSim.NLAtomsPerWarpBitsMask;

  NLWarp* psWarp = &sWarp[warp];
  if (tgx == 0) {
    psWarp->pos = blockIdx.x * (THREADS_PER_BLOCK / GRID) + warp;
  }
  __syncthreads();

  // Massive loop over all neighbor list entries
  while (psWarp->pos < sNLEntries) {

    // Read Neighbor List entry
    if (tgx == 0) {
      psWarp->nlEntry.NL = cSim.pNLEntry[psWarp->pos].NL;
      psWarp->bHomeCell = psWarp->nlEntry.NL.ymax & NLENTRY_HOME_CELL_MASK;
      psWarp->nlEntry.NL.ymax >>= NLENTRY_YMAX_SHIFT;
    }
    __SYNCWARP(WARP_MASK);
    uint offset = psWarp->nlEntry.NL.offset;

    // Read y atoms into registers
    PMEFloat xi;
    PMEFloat yi;
    PMEFloat zi;
    PMEFloat qi;
    unsigned int LJIDi;
    unsigned int index = psWarp->nlEntry.NL.ypos + (tgx & cSim.NLAtomsPerWarpBitsMask);
    if (index < psWarp->nlEntry.NL.ymax) {
#ifndef use_DPFP
#  if defined(__CUDA_ARCH__) && ((__CUDA_ARCH__ == 700) || (__CUDA_ARCH__ >= 800))
      PMEFloat2 xy = cSim.pAtomXYSP[index];
      PMEFloat2 qljid = cSim.pAtomChargeSPLJID[index];
      zi = cSim.pAtomZSP[index];
#  else
      PMEFloat2 xy = tex1Dfetch<float2>(cSim.texAtomXYSP, index);
      PMEFloat2 qljid = tex1Dfetch<float2>(cSim.texAtomChargeSPLJID, index);
      zi = tex1Dfetch<float>(cSim.texAtomZSP, index);
#  endif
#else
      PMEFloat2 xy = cSim.pAtomXYSP[index];
      PMEFloat2 qljid = cSim.pAtomChargeSPLJID[index];
      zi = cSim.pAtomZSP[index];
#endif
      xi = xy.x;
      yi = xy.y;
      qi = qljid.x;
#ifdef use_DPFP
      LJIDi = __double_as_longlong(qljid.y);
#else
      LJIDi = __float_as_uint(qljid.y);
#endif
    }
    else {
      xi = (PMEFloat)10000.0 * index;
      yi = (PMEFloat)10000.0 * index;
      zi = (PMEFloat)10000.0 * index;
      qi = (PMEFloat)0.0;
      LJIDi = 0;
    }
#ifndef PME_IS_ORTHOGONAL
    // Transform into cartesian space
#  ifdef PME_VIRIAL
    xi = sUcellf[0]*xi + sUcellf[1]*yi + sUcellf[2]*zi;
    yi =                 sUcellf[4]*yi + sUcellf[5]*zi;
    zi =                                 sUcellf[8]*zi;
#  else
    xi = cSim.ucellf[0][0]*xi + cSim.ucellf[0][1]*yi + cSim.ucellf[0][2]*zi;
    yi =                        cSim.ucellf[1][1]*yi + cSim.ucellf[1][2]*zi;
    zi =                                               cSim.ucellf[2][2]*zi;
#  endif
#endif

    // Tile accumulators for forces
    PMEFloat TLx_i = (PMEFloat)0.0;
    PMEFloat TLy_i = (PMEFloat)0.0;
    PMEFloat TLz_i = (PMEFloat)0.0;

    // Special-case first tile: in the branch and loop that follows, each thread in a warp
    // will take the perspective of one "i" atom and loop over 8, 16, or 32 "j" atoms
    // depending on the compilation.
    if (psWarp->bHomeCell) {
      unsigned int exclusion = cSim.pNLAtomList[offset + (tgx & cSim.NLAtomsPerWarpBitsMask)];
      offset += cSim.NLAtomsPerWarp;
      if (PME_ATOMS_PER_WARP < GRID) {
        exclusion >>= iftlim * jgroup;
      }
      // In the first tile, atoms A, B, C, ... N are interacting with each other, not
      // some new group of 8, 16, or 32 atoms.  shAtom, a buffer for each thread to
      // store its j atom in addition to its i atom, is not needed in this context.
      // However, the Lennard Jones index of each i atom gets multiplied by the total
      // number of types to prepare for indexing into the table, so store the original
      // value for future __shfl() reference.
      unsigned int LJIDj = LJIDi;
      LJIDi *= cSim.LJTypes;

      // Tile accumulators for energy and virial components
#ifdef PME_ENERGY
      PMEFloat TLeed  = (PMEFloat)0.0;
      PMEFloat TLevdw = (PMEFloat)0.0;
#endif
#ifdef PME_VIRIAL
      PMEFloat TLvir_11 = (PMEFloat)0.0;
      PMEFloat TLvir_22 = (PMEFloat)0.0;
      PMEFloat TLvir_33 = (PMEFloat)0.0;
#endif
      // The first tile loop uniquely skips the
      // first interaction of particle held by each thread tgx, as it would
      // be each i atom with itself.  In subsequent tiles, the interaction
      // of the i atom in thread tgx with the j atom in that same thread
      // is necessary--they are different atoms.
      // Fast-forward iftlim iterations in the final PME_ATOMS_PER_WARP lanes for the first
      // tile computation.
      int j = ((tgx + iftlim * jgroup + 1) & cSim.NLAtomsPerWarpBitsMask) | joffset;
      int jrec = ((PME_ATOMS_PER_WARP + tgx - (iftlim * jgroup + 1))
                  & cSim.NLAtomsPerWarpBitsMask) | joffset;
      PMEFloat xj = __SHFL(WARP_MASK, xi, j);
      PMEFloat yj = __SHFL(WARP_MASK, yi, j);
      PMEFloat zj = __SHFL(WARP_MASK, zi, j);
      PMEFloat qj = __SHFL(WARP_MASK, qi, j);
      LJIDj       = __SHFL(WARP_MASK, LJIDj, j);
      #pragma unroll 2
      for (int ift = 0; ift < iftlim; ift++) {
        PMEFloat xij = xi - xj;
        PMEFloat yij = yi - yj;
        PMEFloat zij = zi - zj;
        xj = WarpRotateLeft<PME_ATOMS_PER_WARP>(xj);
        yj = WarpRotateLeft<PME_ATOMS_PER_WARP>(yj);
        zj = WarpRotateLeft<PME_ATOMS_PER_WARP>(zj);
        PMEFloat r2 = xij*xij + yij*yij + zij*zij;
        unsigned int index = LJIDi + LJIDj;
        LJIDj = WarpRotateLeft<PME_ATOMS_PER_WARP>(LJIDj);
        PMEFloat qiqj = qi * qj;
        qj = WarpRotateLeft<PME_ATOMS_PER_WARP>(qj);
        PMEFloat df = (PMEFloat)0.0;
        bool inrange = ((r2 < cSim.cut2) && r2 > delta);
        if (inrange) {
#if defined(use_SPFP) && !defined(PME_FSWITCH) && !defined(PME_ENERGY)
          PMEFloat r2inv = (PMEFloat)1.0 / r2;
          uint cidx = 2*(__float_as_uint(r2) >> 18) + (exclusion & 0x1);
#  if defined(__CUDA_ARCH__) && ((__CUDA_ARCH__ == 700) || (__CUDA_ARCH__ >= 800))
          PMEFloat4 coef = cSim.pErfcCoeffsTable[cidx];
#  else
          PMEFloat4 coef = tex1Dfetch<float4>(cSim.texErfcCoeffsTable, cidx);
#  endif
#else
          PMEFloat rinv  = rsqrt(r2);
          PMEFloat r     = r2 * rinv;
          PMEFloat r2inv = rinv * rinv;
#endif
#ifndef use_DPFP
#  if defined(__CUDA_ARCH__) && ((__CUDA_ARCH__ == 700) || (__CUDA_ARCH__ >= 800))
          PMEFloat2 term = cSim.pLJTerm[index];
#  else
          PMEFloat2 term = tex1Dfetch<float2>(cSim.texLJTerm, index);
#  endif
#else
          PMEFloat2 term = cSim.pLJTerm[index];
#endif
          PMEFloat r6inv = r2inv * r2inv * r2inv;
#if defined(use_SPFP) && !defined(PME_FSWITCH) && !defined(PME_ENERGY)
          PMEFloat d_swtch_dx = r2*coef.x + coef.y + r2inv*(coef.z + r2inv*coef.w);
#else
#  ifdef use_DPFP
          PMEFloat swtch = erfc(cSim.ew_coeffSP * r);
#  else
          PMEFloat swtch = fasterfc(r);
#  endif
          PMEFloat d_swtch_dx = cSim.negTwoEw_coeffRsqrtPI * exp(-cSim.ew_coeff2 * r2);
#endif
#ifdef PME_ENERGY
          PMEFloat fnrange = (PMEFloat)inrange;
          // Skip last "inactive" interactions
          if (ift * GRID + tgx >= totalInteractions) {
            fnrange = (PMEFloat)0.0;
          }
#endif
#ifdef PME_FSWITCH
          if (!(exclusion & 0x1)) {
            PMEFloat r3inv = rinv*r2inv;
            PMEFloat r12inv = r6inv*r6inv;
            PMEFloat df12f = (PMEFloat)(-1) * term.x * cSim.cut6invcut6minfswitch6 * r2inv *
                             r6inv * (r6inv - cSim.cut6inv);
            PMEFloat df6f = (PMEFloat)(-1) * term.y * cSim.cut3invcut3minfswitch3 * r3inv *
                            r2inv * (r3inv - cSim.cut3inv);
            PMEFloat df12 = (PMEFloat)(-1) * term.x * r2inv * r12inv;
            PMEFloat df6 = (PMEFloat)(-1) * term.y * r2inv * r6inv;
#  ifdef PME_ENERGY
            PMEFloat f12f = term.x * cSim.cut6invcut6minfswitch6 *
                            (r6inv - cSim.cut6inv) * (r6inv - cSim.cut6inv);
            PMEFloat f6f = term.y * cSim.cut3invcut3minfswitch3 *
                           (r3inv - cSim.cut3inv) * (r3inv - cSim.cut3inv);
            PMEFloat f12 = term.x*r12inv - term.x*cSim.invfswitch6cut6;
            PMEFloat f6 = term.y*r6inv - term.y*cSim.invfswitch3cut3;
#  endif
            if (r2 > cSim.fswitch2) {
              df12 = df12f;
              df6 = df6f;
#  ifdef PME_ENERGY
              f12 = f12f;
              f6 = f6f;
#  endif
            }
#  ifdef PME_ENERGY
            TLevdw += fnrange * (f12*twelvth - f6*sixth);
#  endif
            df += df6 - df12;
          }
#else  // PME_FSWITCH
          if (!(exclusion & 0x1)) {
            PMEFloat f6 = term.y * r6inv;
            PMEFloat f12 = term.x * r6inv * r6inv;
            df += (f12 - f6) * r2inv;
#  ifdef PME_ENERGY
            TLevdw += fnrange * (f12*twelvth - f6*sixth);
#  endif
          }
#endif // PME_FSWITCH
#if defined(use_DPFP) || defined(PME_FSWITCH) || defined(PME_ENERGY)
          else {
            swtch -= (PMEFloat)1.0;
          }
#endif
          // This ends a branch for "not an exclusion"--the non-bonded interaction is
          // to be counted.  0x1 is simply 1 in hexadecimal.

#ifdef PME_ENERGY
          PMEFloat b0 = qiqj * swtch * rinv;
          PMEFloat b1 = b0 - qiqj * d_swtch_dx;
          df += b1 * r2inv;
          TLeed += fnrange * b0;
#else  // PME_ENERGY
#  if defined(use_SPFP) && !defined(PME_FSWITCH)
          df += qiqj * d_swtch_dx;
#  else
          df += qiqj * (swtch * rinv - d_swtch_dx) * r2inv;
#  endif
#endif // PME_ENERGY
#if !defined(use_DPFP) && defined(PME_MINIMIZATION)
          df = max(-10000.0f, min(df, 10000.0f));
#endif
        } // inrange
        // Skip last "inactive" interactions
        if (ift * GRID + tgx >= totalInteractions) {
          df = (PMEFloat)0.0;
        }
        PMEFloat dfdx = df * xij;
        PMEFloat dfdy = df * yij;
        PMEFloat dfdz = df * zij;

        // Accumulate into registers for i atoms only, but accumulate
        // both the action and the equal and opposite reaction
        TLx_i += dfdx;
        TLy_i += dfdy;
        TLz_i += dfdz;
        TLx_i -= __SHFL(WARP_MASK, dfdx, jrec);
        TLy_i -= __SHFL(WARP_MASK, dfdy, jrec);
        TLz_i -= __SHFL(WARP_MASK, dfdz, jrec);
#ifdef PME_VIRIAL
        TLvir_11 -= xij * dfdx;
        TLvir_22 -= yij * dfdy;
        TLvir_33 -= zij * dfdz;
#endif

        // Shift bits one to the right in the exclusion tracker to move on to the next atom.
        exclusion >>= 1;
        jrec = ((jrec + PME_ATOMS_PER_WARP - 1) & cSim.NLAtomsPerWarpBitsMask) | joffset;
      }

      // Commit tile accumulators
#ifdef PME_ENERGY
      // The factor of 1/2 was folded into the Lennard-Jones calculation
      // above, but not the electrostatic calculation.
#  ifdef use_SPFP
      eed  += fast_llrintf(TLeed * ENERGYSCALEF);
      evdw += fast_llrintf(TLevdw * ENERGYSCALEF);
#  else
#    ifdef PME_MINIMIZATION
      TLevdw = modf(TLevdw, &i);
      evdwe += llrint(i);
      TLeed = modf(TLeed, &i);
      eede += llrint(i);
#    endif
      evdw += llrint(TLevdw * ENERGYSCALE);
      eed += llrint(TLeed * ENERGYSCALE);
#  endif
#endif

#ifdef PME_VIRIAL
#  ifdef use_SPFP
      vir_11 += fast_llrintf(TLvir_11 * FORCESCALEF);
      vir_22 += fast_llrintf(TLvir_22 * FORCESCALEF);
      vir_33 += fast_llrintf(TLvir_33 * FORCESCALEF);
#  else 
#    if defined(PME_MINIMIZATION)
      TLvir_11 = modf(TLvir_11, &i);
      vir_11E += llrint(i);
      TLvir_22 = modf(TLvir_22, &i);
      vir_22E += llrint(i);
      TLvir_33 = modf(TLvir_33, &i);
      vir_33E += llrint(i);
#    endif
      vir_11 += llrint(TLvir_11 * FORCESCALE);
      vir_22 += llrint(TLvir_22 * FORCESCALE);
      vir_33 += llrint(TLvir_33 * FORCESCALE);
#  endif
#endif
    }//if (psWarp->bHomeCell)
    else {
      LJIDi *= cSim.LJTypes;
    }
    // This ends the branch for special-casing the first tile to have
    // 32 atoms interact with each other.  LJIDi, the LJ index for the
    // ith atom, has to be set for future operations, no matter what.

    // Handle the remainder of the line
    int tx = 0;
    while (tx < psWarp->nlEntry.NL.xatoms) {

      // Read atom ID and exclusion data
      NLAtom shAtom;
      shAtom.ID = cSim.pNLAtomList[offset + tgx];
      offset += GRID;
      PMEMask fullExclusion =
        ((PMEMask*)&cSim.pNLAtomList[offset])[tgx & cSim.NLAtomsPerWarpBitsMask];
      offset += cSim.NLAtomsPerWarp * sizeof(PMEMask) / sizeof(unsigned int);
      if (PME_ATOMS_PER_WARP < GRID) {
        fullExclusion >>= joffset;
      }
      unsigned int exclusion = (unsigned int)fullExclusion & cSim.NLAtomsPerWarpMask;
      // Clear j atom forces
      PMEFloat shFx = (PMEFloat)0.0;
      PMEFloat shFy = (PMEFloat)0.0;
      PMEFloat shFz = (PMEFloat)0.0;

      // Read shared memory data
      if (tx + tgx < psWarp->nlEntry.NL.xatoms) {
        unsigned int atom = shAtom.ID >> NLATOM_CELL_SHIFT;
#ifndef use_DPFP
#  if defined(__CUDA_ARCH__) && ((__CUDA_ARCH__ == 700) || (__CUDA_ARCH__ >= 800))
        PMEFloat2 xy = cSim.pAtomXYSP[atom];
        PMEFloat2 qljid = cSim.pAtomChargeSPLJID[atom];
        shAtom.z = cSim.pAtomZSP[atom];
#  else
        PMEFloat2 xy = tex1Dfetch<float2>(cSim.texAtomXYSP, atom);
        PMEFloat2 qljid = tex1Dfetch<float2>(cSim.texAtomChargeSPLJID, atom);
        shAtom.z = tex1Dfetch<float>(cSim.texAtomZSP, atom);
#  endif
#else
        PMEFloat2 xy = cSim.pAtomXYSP[atom];
        PMEFloat2 qljid = cSim.pAtomChargeSPLJID[atom];
        shAtom.z = cSim.pAtomZSP[atom];
#endif
        shAtom.x = xy.x;
        shAtom.y = xy.y;
        shAtom.q = qljid.x;
#ifdef use_DPFP
        shAtom.LJID = __double_as_longlong(qljid.y);
#else
        shAtom.LJID = __float_as_uint(qljid.y);
#endif
      }
      else {
        shAtom.x = (PMEFloat)-10000.0 * tgx;
        shAtom.y = (PMEFloat)-10000.0 * tgx;
        shAtom.z = (PMEFloat)-10000.0 * tgx;
        shAtom.q = (PMEFloat)0.0;
        shAtom.LJID = 0;
      }

      // Translate all atoms into a local coordinate system within one unit
      // cell of the first atom read to avoid PBC handling within inner loops
      int cell = shAtom.ID & NLATOM_CELL_TYPE_MASK;
#if defined(PME_VIRIAL) && defined(PME_IS_ORTHOGONAL)
      shAtom.x += sUcellf[0] * cSim.cellOffset[cell][0];
      shAtom.y += sUcellf[4] * cSim.cellOffset[cell][1];
      shAtom.z += sUcellf[8] * cSim.cellOffset[cell][2];
#else
      shAtom.x += cSim.cellOffset[cell][0];
      shAtom.y += cSim.cellOffset[cell][1];
      shAtom.z += cSim.cellOffset[cell][2];
#endif
#ifndef PME_IS_ORTHOGONAL
#  ifdef PME_VIRIAL
      shAtom.x = sUcellf[0]*shAtom.x + sUcellf[1]*shAtom.y + sUcellf[2]*shAtom.z;
      shAtom.y = sUcellf[4]*shAtom.y + sUcellf[5]*shAtom.z;
      shAtom.z = sUcellf[8]*shAtom.z;
#  else
      shAtom.x = cSim.ucellf[0][0]*shAtom.x + cSim.ucellf[0][1]*shAtom.y +
                 cSim.ucellf[0][2]*shAtom.z;
      shAtom.y = cSim.ucellf[1][1]*shAtom.y + cSim.ucellf[1][2]*shAtom.z;
      shAtom.z = cSim.ucellf[2][2]*shAtom.z;
#  endif
#endif

#ifdef PME_ENERGY
      PMEFloat TLeed  = (PMEFloat)0.0;
      PMEFloat TLevdw = (PMEFloat)0.0;
#endif
#ifdef PME_VIRIAL
      PMEFloat TLvir_11 = (PMEFloat)0.0;
      PMEFloat TLvir_22 = (PMEFloat)0.0;
      PMEFloat TLvir_33 = (PMEFloat)0.0;
#endif

      if (__ANY(WARP_MASK, exclusion)) {
        #pragma unroll 2
        for (int it = 0; it < PME_ATOMS_PER_WARP; it++) {
          PMEFloat xij = xi - shAtom.x;
          PMEFloat yij = yi - shAtom.y;
          PMEFloat zij = zi - shAtom.z;
          shAtom.x = WarpRotateLeft<PME_ATOMS_PER_WARP>(shAtom.x);
          shAtom.y = WarpRotateLeft<PME_ATOMS_PER_WARP>(shAtom.y);
          shAtom.z = WarpRotateLeft<PME_ATOMS_PER_WARP>(shAtom.z);
          PMEFloat r2  = xij*xij + yij*yij + zij*zij;
          unsigned int index = LJIDi + shAtom.LJID;
          shAtom.LJID = WarpRotateLeft<PME_ATOMS_PER_WARP>(shAtom.LJID);
          PMEFloat qiqj = qi * shAtom.q;
          shAtom.q = WarpRotateLeft<PME_ATOMS_PER_WARP>(shAtom.q);
          PMEFloat df = (PMEFloat)0.0;
          bool inrange = ((r2 < cSim.cut2) && r2 > delta);
          if (inrange) {
#if defined(use_SPFP) && !defined(PME_FSWITCH) && !defined(PME_ENERGY)
            PMEFloat r2inv = (PMEFloat)1.0 / r2;
            uint cidx = 2*(__float_as_uint(r2) >> 18) + (exclusion & 0x1);
#  if defined(__CUDA_ARCH__) && ((__CUDA_ARCH__ == 700) || (__CUDA_ARCH__ >= 800))
            PMEFloat4 coef = cSim.pErfcCoeffsTable[cidx];
#  else
            PMEFloat4 coef = tex1Dfetch<float4>(cSim.texErfcCoeffsTable, cidx);
#  endif
#else
            PMEFloat rinv      = rsqrt(r2);
            PMEFloat r         = r2 * rinv;
            PMEFloat r2inv     = rinv * rinv;
#endif
#ifndef use_DPFP
#  if defined(__CUDA_ARCH__) && ((__CUDA_ARCH__ == 700) || (__CUDA_ARCH__ >= 800))
            PMEFloat2 term = cSim.pLJTerm[index];
#  else
            PMEFloat2 term = tex1Dfetch<float2>(cSim.texLJTerm, index);
#  endif
#else
            PMEFloat2 term = cSim.pLJTerm[index];
#endif
            PMEFloat r6inv     = r2inv * r2inv * r2inv;
#if defined(use_SPFP) && !defined(PME_FSWITCH) && !defined(PME_ENERGY)
            PMEFloat d_swtch_dx = r2*coef.x + coef.y + r2inv*(coef.z + r2inv*coef.w);
#else
#  ifdef use_DPFP
            PMEFloat swtch = erfc(cSim.ew_coeffSP * r);
#  else
            PMEFloat swtch = fasterfc(r);
#  endif
            PMEFloat d_swtch_dx = cSim.negTwoEw_coeffRsqrtPI * exp(-cSim.ew_coeff2 * r2);
#endif
#ifdef PME_ENERGY
            PMEFloat fnrange = (PMEFloat)inrange;
#endif
#ifdef PME_FSWITCH
            if (!(exclusion & 0x1)) {
              PMEFloat r3inv = rinv*r2inv;
              PMEFloat r12inv = r6inv*r6inv;
              PMEFloat df12f = (PMEFloat)(-1) * term.x * cSim.cut6invcut6minfswitch6 *
                               r2inv * r6inv * (r6inv - cSim.cut6inv);
              PMEFloat df6f = (PMEFloat)(-1) * term.y * cSim.cut3invcut3minfswitch3 *
                              r3inv * r2inv * (r3inv - cSim.cut3inv);
              PMEFloat df12 = (PMEFloat)(-1) * term.x * r2inv * r12inv;
              PMEFloat df6 = (PMEFloat)(-1) * term.y * r2inv * r6inv;
#  ifdef PME_ENERGY
              PMEFloat f12f = term.x * cSim.cut6invcut6minfswitch6 *
                              (r6inv - cSim.cut6inv)*(r6inv - cSim.cut6inv);
              PMEFloat f6f = term.y * cSim.cut3invcut3minfswitch3 *
                             (r3inv - cSim.cut3inv)*(r3inv - cSim.cut3inv);
              PMEFloat f12 = (term.x * r12inv) - (term.x * cSim.invfswitch6cut6);
              PMEFloat f6 = (term.y * r6inv) - (term.y * cSim.invfswitch3cut3);
#  endif
              if (r2 > cSim.fswitch2) {
                df12 = df12f;
                df6 = df6f;
#  ifdef PME_ENERGY
                f12 = f12f;
                f6 = f6f;
#  endif
              }
              df += df6 - df12;
#  ifdef PME_ENERGY
              TLevdw += fnrange * (f12*twelvth - f6*sixth);
#  endif
            }
#else  // PME_FSWITCH
            if (!(exclusion & 0x1)) {
              PMEFloat f6 = term.y * r6inv;
              PMEFloat f12 = term.x * r6inv * r6inv;
              df += (f12 - f6) * r2inv;
#  ifdef PME_ENERGY
              TLevdw += fnrange * (f12*twelvth - f6*sixth);
#  endif
            }
#endif // PME_FSWITCH
#if defined(use_DPFP) || defined(PME_FSWITCH) || defined(PME_ENERGY)
            else {
              swtch -= (PMEFloat)1.0;
            }
#endif
#ifdef PME_ENERGY
            PMEFloat b0 = qiqj * swtch * rinv;
            PMEFloat b1 = b0 - qiqj * d_swtch_dx;
            df += b1 * r2inv;
            TLeed += fnrange * b0;
#else  // PME_ENERGY
#  if defined(use_SPFP) && !defined(PME_FSWITCH)
            df += qiqj * d_swtch_dx;
#  else
            df += qiqj*(swtch*rinv - d_swtch_dx)*r2inv;
#  endif
#endif // PME_ENERGY
#if !defined(use_DPFP) && defined(PME_MINIMIZATION)
            df = max(-10000.0f, min(df, 10000.0f));
#endif
          }
          PMEFloat dfdx = df * xij;
          PMEFloat dfdy = df * yij;
          PMEFloat dfdz = df * zij;
          TLx_i += dfdx;
          TLy_i += dfdy;
          TLz_i += dfdz;
          shFx -= dfdx;
          shFy -= dfdy;
          shFz -= dfdz;
          shFx = WarpRotateLeft<PME_ATOMS_PER_WARP>(shFx);
          shFy = WarpRotateLeft<PME_ATOMS_PER_WARP>(shFy);
          shFz = WarpRotateLeft<PME_ATOMS_PER_WARP>(shFz);
#ifdef PME_VIRIAL
          TLvir_11 -= xij * dfdx;
          TLvir_22 -= yij * dfdy;
          TLvir_33 -= zij * dfdz;
#endif
          exclusion >>= 1;
        }
        // End for loop covering non-bonded computations when
        // there IS at least one exclusion somewhere in the pile
      }
      else {
        #pragma unroll 2
        for (int it = 0; it < PME_ATOMS_PER_WARP; it++) {
          // Read properties for the other atom
          PMEFloat xij = xi - shAtom.x;
          PMEFloat yij = yi - shAtom.y;
          PMEFloat zij = zi - shAtom.z;
          shAtom.x = WarpRotateLeft<PME_ATOMS_PER_WARP>(shAtom.x);
          shAtom.y = WarpRotateLeft<PME_ATOMS_PER_WARP>(shAtom.y);
          shAtom.z = WarpRotateLeft<PME_ATOMS_PER_WARP>(shAtom.z);

          // Perform the range test
          PMEFloat r2  = xij*xij + yij*yij + zij*zij;
          unsigned int index = LJIDi + shAtom.LJID;
          shAtom.LJID = WarpRotateLeft<PME_ATOMS_PER_WARP>(shAtom.LJID);
          PMEFloat qiqj = qi * shAtom.q;
          shAtom.q = WarpRotateLeft<PME_ATOMS_PER_WARP>(shAtom.q);
          bool inrange = ((r2 < cSim.cut2) && r2 > delta);
          PMEFloat df = (PMEFloat)0.0;
          if (inrange) {
#if defined(use_SPFP) && !defined(PME_FSWITCH) && !defined(PME_ENERGY)
            PMEFloat r2inv = (PMEFloat)1.0 / r2;
            uint cidx = 2*(__float_as_uint(r2) >> 18);
#  if defined(__CUDA_ARCH__) && ((__CUDA_ARCH__ == 700) || (__CUDA_ARCH__ >= 800))
            PMEFloat4 coef = cSim.pErfcCoeffsTable[cidx];
#  else
            PMEFloat4 coef = tex1Dfetch<float4>(cSim.texErfcCoeffsTable, cidx);
#  endif
#else
            PMEFloat rinv      = rsqrt(r2);
            PMEFloat r         = r2 * rinv;
            PMEFloat r2inv     = rinv * rinv;
#endif
#ifndef use_DPFP
#  if defined(__CUDA_ARCH__) && ((__CUDA_ARCH__ == 700) || (__CUDA_ARCH__ >= 800))
            PMEFloat2 term = cSim.pLJTerm[index];
#  else
            PMEFloat2 term = tex1Dfetch<float2>(cSim.texLJTerm, index);
#  endif
#else
            PMEFloat2 term = cSim.pLJTerm[index];
#endif
            PMEFloat r6inv     = r2inv * r2inv * r2inv;
#if defined(use_SPFP) && !defined(PME_FSWITCH) && !defined(PME_ENERGY)
            PMEFloat d_swtch_dx = r2*coef.x + coef.y + r2inv*(coef.z + r2inv*coef.w);
#else
#  ifdef use_DPFP
            PMEFloat swtch = erfc(cSim.ew_coeffSP * r) * rinv;
#  else
            PMEFloat swtch = fasterfc(r) * rinv;
#  endif
            PMEFloat d_swtch_dx = cSim.negTwoEw_coeffRsqrtPI * exp(-cSim.ew_coeff2 * r2);
#endif
#ifdef PME_ENERGY
            PMEFloat fnrange = (PMEFloat)inrange;
#endif
#ifdef PME_FSWITCH
            PMEFloat r3inv = rinv*r2inv;
            PMEFloat r12inv = r6inv*r6inv;
            PMEFloat df12f = (PMEFloat)(-1) * term.x * cSim.cut6invcut6minfswitch6 *
                             r2inv * r6inv * (r6inv - cSim.cut6inv);
            PMEFloat df6f = (PMEFloat)(-1) * term.y * cSim.cut3invcut3minfswitch3 *
                            r3inv * r2inv * (r3inv - cSim.cut3inv);
            PMEFloat df12 = (PMEFloat)(-1) * term.x * r2inv * r12inv;
            PMEFloat df6 = (PMEFloat)(-1) * term.y * r2inv * r6inv;
#  ifdef PME_ENERGY
            PMEFloat f12f = term.x * cSim.cut6invcut6minfswitch6 *
                            (r6inv - cSim.cut6inv)*(r6inv - cSim.cut6inv);
            PMEFloat f6f = term.y * cSim.cut3invcut3minfswitch3 *
                           (r3inv - cSim.cut3inv)*(r3inv - cSim.cut3inv);
            PMEFloat f12 = term.x * r12inv - term.x * cSim.invfswitch6cut6;
            PMEFloat f6 = term.y * r6inv - term.y * cSim.invfswitch3cut3;
#  endif
            if (r2 > cSim.fswitch2) {
              df12 = df12f;
              df6 = df6f;
#  ifdef PME_ENERGY
              f12 = f12f;
              f6 = f6f;
#  endif
            }
            df += df6 - df12;
#  ifdef PME_ENERGY
            TLevdw += fnrange * (f12*twelvth - f6*sixth);
            PMEFloat b0 = qiqj * swtch;
            PMEFloat b1 = b0 - qiqj * d_swtch_dx;
            df += b1 * r2inv;
            TLeed += fnrange * b0;
#  else  // PME_ENERGY
            df += qiqj * (swtch - d_swtch_dx) * r2inv;
#  endif // PME_ENERGY
#else // PME_FSWITCH
            PMEFloat f6 = term.y * r6inv;
            PMEFloat f12 = term.x * r6inv * r6inv;
            df += (f12-f6) * r2inv;
#  ifdef PME_ENERGY
            TLevdw += fnrange * (f12*twelvth - f6*sixth);
            PMEFloat b0 = qiqj * swtch;
            PMEFloat b1 = b0 - qiqj * d_swtch_dx;
            df += b1 * r2inv;
            TLeed += fnrange * b0;
#  else  // PME_ENERGY
#    if defined(use_SPFP) && !defined(PME_FSWITCH)
            df += qiqj * d_swtch_dx;
#    else
            df += qiqj * (swtch - d_swtch_dx) * r2inv;
#    endif
#  endif // PME_ENERGY
#endif //PME_FSWITCH
#if !defined(use_DPFP) && defined(PME_MINIMIZATION)
            df = max(-10000.0f, min(df, 10000.0f));
#endif
          }
          PMEFloat dfdx = df * xij;
          PMEFloat dfdy = df * yij;
          PMEFloat dfdz = df * zij;
          TLx_i += dfdx;
          TLy_i += dfdy;
          TLz_i += dfdz;
          shFx -= dfdx;
          shFy -= dfdy;
          shFz -= dfdz;
          shFx = WarpRotateLeft<PME_ATOMS_PER_WARP>(shFx);
          shFy = WarpRotateLeft<PME_ATOMS_PER_WARP>(shFy);
          shFz = WarpRotateLeft<PME_ATOMS_PER_WARP>(shFz);
#ifdef PME_VIRIAL
          TLvir_11 -= xij * dfdx;
          TLvir_22 -= yij * dfdy;
          TLvir_33 -= zij * dfdz;
#endif
        }
        // Ends for loop for processing a pile of non-bonded
        // interactions with no exclusions to worry about
      }
      // End branch dealing with the presence of exclusions
      // in the pile of non-bonded interactions.

#ifdef PME_ENERGY
#  ifdef use_SPFP
      eed  += fast_llrintf(TLeed * ENERGYSCALEF);
      evdw += fast_llrintf(TLevdw * ENERGYSCALEF);
#  else
#    ifdef PME_MINIMIZATION
      TLevdw = modf(TLevdw, &i);
      evdwe += llrint(i);
      TLeed = modf(TLeed, &i);
      eede += llrint(i);
#    endif
      eed  += llrint(TLeed * ENERGYSCALE);
      evdw += llrint(TLevdw * ENERGYSCALE);
#  endif
#endif

#ifdef PME_VIRIAL
#  ifdef use_SPFP
      vir_11 += fast_llrintf(TLvir_11 * FORCESCALEF);
      vir_22 += fast_llrintf(TLvir_22 * FORCESCALEF);
      vir_33 += fast_llrintf(TLvir_33 * FORCESCALEF);
#  else
#    if defined(PME_MINIMIZATION)
      TLvir_11 = modf(TLvir_11, &i);
      vir_11E += llrint(i);
      TLvir_22 = modf(TLvir_22, &i);
      vir_22E += llrint(i);
      TLvir_33 = modf(TLvir_33, &i);
      vir_33E += llrint(i);
#    endif
      vir_11 += llrint(TLvir_11 * FORCESCALE);
      vir_22 += llrint(TLvir_22 * FORCESCALE);
      vir_33 += llrint(TLvir_33 * FORCESCALE);
#  endif
#endif
      // Dump j atom forces
      if (tx + tgx < psWarp->nlEntry.NL.xatoms) {
        int offset = (shAtom.ID >> NLATOM_CELL_SHIFT);
#ifdef use_SPFP
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset],
                  llitoulli(fast_llrintf(shFx * FORCESCALEF)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset],
                  llitoulli(fast_llrintf(shFy * FORCESCALEF)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset],
                  llitoulli(fast_llrintf(shFz * FORCESCALEF)));
#else // use_DPFP
#  ifdef PME_MINIMIZATION
        shFx = modf(shFx, &i);
        atomicAdd((unsigned long long int*)&cSim.pIntForceXAccumulator[offset],
                  llitoulli(llrint(i)));
        shFy = modf(shFy, &i);
        atomicAdd((unsigned long long int*)&cSim.pIntForceYAccumulator[offset],
                  llitoulli(llrint(i)));
        shFz = modf(shFz, &i);
        atomicAdd((unsigned long long int*)&cSim.pIntForceZAccumulator[offset],
                  llitoulli(llrint(i)));
#  endif
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset],
                  llitoulli(llrint(shFx * FORCESCALE)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset],
                  llitoulli(llrint(shFy * FORCESCALE)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset],
                  llitoulli(llrint(shFz * FORCESCALE)));
#endif
      }
      // End contingency for dumping j atom forces to global

      // Advance to next x tile
      tx += GRID;
    }

    // Reduce register forces
    for (unsigned int rlev = GRID >> 1; rlev >= PME_ATOMS_PER_WARP; rlev /= 2) {
      TLx_i += __SHFL(WARP_MASK, TLx_i, tgx + rlev);
      TLy_i += __SHFL(WARP_MASK, TLy_i, tgx + rlev);
      TLz_i += __SHFL(WARP_MASK, TLz_i, tgx + rlev);
    }

    PMEForceAccumulator fx_i, fy_i, fz_i;
#ifdef use_SPFP
    fx_i = fast_llrintf(TLx_i * FORCESCALEF);
    fy_i = fast_llrintf(TLy_i * FORCESCALEF);
    fz_i = fast_llrintf(TLz_i * FORCESCALEF);
#else
#  ifdef PME_MINIMIZATION
    TLx_i = modf(TLx_i, &i);
    long long int fxe_i = llrint(i);
    TLy_i = modf(TLy_i, &i);
    long long int fye_i = llrint(i);
    TLz_i = modf(TLz_i, &i);
    long long int fze_i = llrint(i);
#  endif
    fx_i = llrint(TLx_i * FORCESCALE);
    fy_i = llrint(TLy_i * FORCESCALE);
    fz_i = llrint(TLz_i * FORCESCALE);
#endif

    // Dump register forces
    if (psWarp->nlEntry.NL.ypos + tgx < psWarp->nlEntry.NL.ymax) {
      int offset = psWarp->nlEntry.NL.ypos + tgx;
      atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset],
                llitoulli(fx_i));
      atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset],
                llitoulli(fy_i));
      atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset],
                llitoulli(fz_i));
#if defined(use_DPFP) && defined(PME_MINIMIZATION)
      atomicAdd((unsigned long long int*)&cSim.pIntForceXAccumulator[offset],
                llitoulli(fxe_i));
      atomicAdd((unsigned long long int*)&cSim.pIntForceYAccumulator[offset],
                llitoulli(fye_i));
      atomicAdd((unsigned long long int*)&cSim.pIntForceZAccumulator[offset],
                llitoulli(fze_i));
#endif
    }
    // End of contingency for dumping register forces.

    // Get next Neighbor List entry
    if (tgx == 0) {
      psWarp->pos = atomicAdd(&cSim.pFrcBlkCounters[0], 1);
    }
    __SYNCWARP(WARP_MASK);
  }
  // End of massive while loop iterating psWarp->pos up to sNLEntries

  // Reduce virials and/or energies
#if defined(PME_VIRIAL) || defined(PME_ENERGY)
  for (unsigned int rlev = GRID >> 1; rlev > 0; rlev /= 2) {
#  ifdef PME_ENERGY
#    if defined(use_DPFP) && defined(PME_MINIMIZATION)
    eede += __SHFL(WARP_MASK, eede, tgx + rlev);
    evdwe += __SHFL(WARP_MASK, evdwe, tgx + rlev);
#    endif
    eed += __SHFL(WARP_MASK, eed, tgx + rlev);
    evdw += __SHFL(WARP_MASK, evdw, tgx + rlev);
#  endif

#  ifdef PME_VIRIAL
#    if defined(use_DPFP) && defined(PME_MINIMIZATION)
    vir_11E += __SHFL(WARP_MASK, vir_11E, tgx + rlev);
    vir_22E += __SHFL(WARP_MASK, vir_22E, tgx + rlev);
    vir_33E += __SHFL(WARP_MASK, vir_33E, tgx + rlev);
#    endif
    vir_11 += __SHFL(WARP_MASK, vir_11, tgx + rlev);
    vir_22 += __SHFL(WARP_MASK, vir_22, tgx + rlev);
    vir_33 += __SHFL(WARP_MASK, vir_33, tgx + rlev);
#  endif
  }
#endif

#ifdef PME_VIRIAL
  if (tgx == 0) {
#  if defined(use_DPFP) && defined(PME_MINIMIZATION)
    atomicAdd(cSim.pVirial_11E, llitoulli(vir_11E));
    atomicAdd(cSim.pVirial_22E, llitoulli(vir_22E));
    atomicAdd(cSim.pVirial_33E, llitoulli(vir_33E));
#  endif
    atomicAdd(cSim.pVirial_11, llitoulli(vir_11));
    atomicAdd(cSim.pVirial_22, llitoulli(vir_22));
    atomicAdd(cSim.pVirial_33, llitoulli(vir_33));
  }
#endif

#ifdef PME_ENERGY
  if (tgx == 0) {
    atomicAdd(cSim.pEVDW, llitoulli(evdw));
    atomicAdd(cSim.pEED, llitoulli(eed));
#  if defined(use_DPFP) && defined(PME_MINIMIZATION)
    atomicAdd(cSim.pEVDWE, llitoulli(evdwe));
    atomicAdd(cSim.pEEDE, llitoulli(eede));
#  endif
  }
#endif // PME_ENERGY
}
