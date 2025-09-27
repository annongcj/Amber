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
#if !defined(AMBER_PLATFORM_AMD)
#  define VOLATILE volatile
#else
#  define VOLATILE
#endif
#define uint unsigned int
  struct NLAtom {
    PMEFloat x;
    PMEFloat y;
    PMEFloat z;
    PMEFloat q;
#ifdef PHMD
    PMEFloat qj_phmd;
    PMEFloat qg;
    PMEFloat qxg;
    PMEFloat lgtauto;
    PMEFloat igtaut;
    PMEFloat xg;
    PMEFloat factg;
    PMEFloat dxg;
    bool lgtitr;
    int g;
    int psp_grp_g;
    unsigned int indexj;
#endif
    unsigned int LJID;
    unsigned int ID;
  };

  struct NLWarp {
    NLEntry nlEntry;
    uint pos;
    bool bHomeCell;
  };

  const PMEFloat delta = 1.0e-5;
#if defined(PME_ENERGY) || defined(PHMD)
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
  __shared__ VOLATILE NLWarp sWarp[THREADS_PER_BLOCK / GRID];

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

  VOLATILE NLWarp* psWarp = &sWarp[warp];
  if (tgx == 0) {
    psWarp->pos = blockIdx.x * (THREADS_PER_BLOCK / GRID) + warp;
  }
  __syncthreads();

  // Massive loop over all neighbor list entries
  while (psWarp->pos < sNLEntries) {

    // Read Neighbor List entry
    if (tgx < 4) {
      psWarp->nlEntry.array[tgx] = cSim.pNLEntry[psWarp->pos].array[tgx];
    }
    __SYNCWARP(WARP_MASK);
    if (tgx == 0) {
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
    PMEForceAccumulator fx_i, fy_i, fz_i;
    fx_i = (PMEForceAccumulator)0;
    fy_i = (PMEForceAccumulator)0;
    fz_i = (PMEForceAccumulator)0;
#ifdef PHMD
    PMEFloat dudli, dudliplus;
    PMEFloat qh, qxh, lambda, qprot, qunprot, qi_phmd, qj_phmd;
    PMEFloat dudlj, dudljplus, factphmd;
    PMEFloat xh, ihtaut, facth, dxh, xg, factg, dxg;
    PMEFloat2 vstate1;
    PMEFloat lambdah, de, dfphmd;
    PMEFloat2 pqstate1, pqstate2;
    PMEFloat radh = 0;
    int h, psp_grp, psp_grp_g;
    bool lhtauto, lgtauto;
    bool lhtitr, lgtitr;
    PMEFloat x2 = 1.0;
#endif
#if defined(use_DPFP) && defined(PME_MINIMIZATION)
    long long int fxe_i = 0;
    long long int fye_i = 0;
    long long int fze_i = 0;
#endif
    unsigned int index = psWarp->nlEntry.NL.ypos + (tgx & cSim.NLAtomsPerWarpBitsMask);
#ifdef PHMD
      int indexi = index;
#endif
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
#ifdef PHMD
      h = cSim.pImageGrplist[index] - 1;
      qi_phmd = cSim.pImageCharge_phmd[index];
      qh = 0;
      qxh = 0;
      lhtitr = false;
      dudli = 0;
      dudliplus = 0;
      if (h >= 0) {
        pqstate1 = cSim.pImageQstate1[index];
        pqstate2 = cSim.pImageQstate2[index];
        psp_grp = cSim.psp_grp[h];
        vstate1 = cSim.pvstate1[cSim.pImageIndexPHMD[index]];
        radh = vstate1.x;
        lambda = sin(cSim.pph_theta[h]);
        lambda *= lambda;
        lambdah = lambda;
        facth = 1.0 - lambdah;
        xh = 1.0;
        if (psp_grp > 0) {
          lhtauto = true;
          x2 = sin(cSim.pph_theta[h+1]);
          x2 *= x2;
          radh += vstate1.y;
          if (vstate1.x > 0) {
            ihtaut = 1.0;
            xh = x2;
          }
          else if (vstate1.y > 0) {
            ihtaut = -1;
            xh = 1.0 - x2;
          }
          if (psp_grp == 1) {
            facth = 1.0 - lambda * xh;
            dxh = -ihtaut * lambda;
          }
          else if (psp_grp == 3) {
            facth = (1.0 - lambda) * xh;
            dxh = ihtaut * (1.0 - lambda);
          }
          qunprot = lambda * (pqstate2.x - pqstate2.y);
          qprot = (1 - lambda) * (pqstate1.x - pqstate1.y);
          qxh = qunprot + qprot;
        }
        lhtitr = (radh > 0);
        qunprot = x2 * pqstate2.x + (1.0 - x2) * pqstate2.y;
        qprot = x2 * pqstate1.x + (1.0 - x2) * pqstate1.y;
        qh = qunprot - qprot;
      }
#endif
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
#ifdef PHMD
      h = -1;
#endif
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

      // Tile accumulators for forces, energy, and virial components
      PMEFloat TLx_i  = (PMEFloat)0.0;
      PMEFloat TLy_i  = (PMEFloat)0.0;
      PMEFloat TLz_i  = (PMEFloat)0.0;
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
      #pragma unroll 2
      for (int ift = 0; ift < iftlim; ift++) {
        PMEFloat xij = xi - __SHFL(WARP_MASK, xi, j);
        PMEFloat yij = yi - __SHFL(WARP_MASK, yi, j);
        PMEFloat zij = zi - __SHFL(WARP_MASK, zi, j);
        PMEFloat r2  = xij*xij + yij*yij + zij*zij;
        unsigned int index = LJIDi + __SHFL(WARP_MASK, LJIDj, j);
        PMEFloat qiqj = qi * __SHFL(WARP_MASK, qi, j);
        PMEFloat df = (PMEFloat)0.0;
        int inrange = ((r2 < cSim.cut2) && r2 > delta);
#ifdef PHMD
        int g = __SHFL(0xFFFFFFFF, h, j);
        int indexj = __SHFL(0xFFFFFFFF, indexi, j);
        PMEFloat qg = __SHFL(0xFFFFFFFF, qh, j);
        PMEFloat qxg = __SHFL(0xFFFFFFFF, qxh, j);
        qj_phmd = __SHFL(0xFFFFFFFF, qi_phmd, j);
        psp_grp_g = __SHFL(0xFFFFFFFF, psp_grp, j);
        xg = __SHFL(0xFFFFFFFF, xh, j);
        dxg = __SHFL(0xFFFFFFFF, dxh, j);
        factg = __SHFL(0xFFFFFFFF, facth, j);
        lgtauto = __SHFL(0xFFFFFFFF, lhtauto, j);
        lgtitr = __SHFL(0xFFFFFFFF, lhtitr, j);
#endif
        if (inrange) {
#if defined(use_SPFP) && !defined(PME_FSWITCH) && !defined(PME_ENERGY)
          PMEFloat r2inv = (PMEFloat)1.0 / r2;
          uint cidx = 2*(__float_as_uint(r2) >> 18) + (exclusion & 0x1);
#  if defined(__CUDA_ARCH__) && ((__CUDA_ARCH__ == 700) || (__CUDA_ARCH__ >= 800))
          PMEFloat4 coef = cSim.pErfcCoeffsTable[cidx];
#  else
          PMEFloat4 coef = tex1Dfetch<float4>(cSim.texErfcCoeffsTable, cidx);
#  endif
#  ifdef PHMD
          PMEFloat r     = sqrt(r2);
          PMEFloat rinv  = (PMEFloat)1.0 / r;
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
#  ifdef PHMD
#    ifdef use_DPFP
          PMEFloat swtch = erfc(cSim.ew_coeffSP * r);
#    else
          PMEFloat swtch = fasterfc(r);
#    endif
#  endif
#else
#  ifdef use_DPFP
          PMEFloat swtch = erfc(cSim.ew_coeffSP * r);
#  else
          PMEFloat swtch = fasterfc(r);
#  endif
          PMEFloat d_swtch_dx = cSim.negTwoEw_coeffRsqrtPI * exp(-cSim.ew_coeff2 * r2);
#endif
#if defined(PME_ENERGY) || defined(PHMD)
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
#  if defined(PME_ENERGY) || defined(PHMD)
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
#  ifdef PHMD
            dfphmd = df6 - df12;
            de = fnrange * (f12 * twelvth - f6 * sixth);
#  endif
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
#  ifdef PHMD
            dfphmd = (f12 - f6) * r2inv;
            de = fnrange * (f12 * twelvth - f6 * sixth);
#  endif
#  ifdef PME_ENERGY
            TLevdw += fnrange * (f12*twelvth - f6*sixth);
#  endif
          }
#endif // PME_FSWITCH
#if defined(use_DPFP) || defined(PME_FSWITCH) || defined(PME_ENERGY) || defined(PHMD)
          else {
            swtch -= (PMEFloat)1.0;
          }
#endif
          // This ends a branch for "not an exclusion"--the non-bonded interaction is
          // to be counted.  0x1 is simply 1 in hexadecimal.

#ifdef PHMD
          //First
          dudlj = 0;
          dudljplus = 0;
          factphmd = 1.0;
          if (!(exclusion & 0x1)) {
            if (lhtitr && lgtitr) {
              if (h != g) {
                factphmd = facth * factg;
                dudli -= xh * factg * de;
                dudlj -= xg * facth * de;
                if (lhtauto) {
                  dudliplus += dxh * factg * de;
                }
                if (lgtauto) {
                  dudljplus += dxg * facth * de;
                }
              }
              else if (psp_grp == 1) {
                factphmd = 1.0 - lambdah;
                dudli -= de;
              }
            }
            else if (lhtitr) {
              factphmd = facth;
              dudli -= xh * de;
              if (lhtauto) {
                dudliplus += dxh * de;
              }
            }
            else if (lgtitr) {
              factphmd = factg;
              dudlj -= xg * de;
              if (lgtauto) {
                dudljplus += dxg * de;
              }
            }
            df += (factphmd - 1.0) * dfphmd;
#  ifdef PME_ENERGY
            TLevdw += (factphmd - 1.0) * de;
#   endif
          }
          if (h >= 0) {
            dudli += qj_phmd * qh * swtch * rinv * fnrange;
            if (psp_grp > 0) {
              dudliplus += qj_phmd * qxh * swtch * rinv * fnrange;
            }
          }
          if (g >= 0) {
            dudlj += qi_phmd * qg * swtch * rinv * fnrange;
            if (psp_grp_g > 0) {
              dudljplus += qi_phmd * qxg * swtch * rinv * fnrange;
          }
          atomicAdd((unsigned long long int*)&cSim.pdph_accumulator[indexj],
                    llitoulli(fast_llrintf(dudlj * FORCESCALEF)));
          atomicAdd((unsigned long long int*)&cSim.pdph_plus_accumulator[indexj],
                    llitoulli(fast_llrintf(dudljplus * FORCESCALEF)));
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
        j = ((j + 1) & cSim.NLAtomsPerWarpBitsMask) | joffset;
        jrec = ((jrec + PME_ATOMS_PER_WARP - 1) & cSim.NLAtomsPerWarpBitsMask) | joffset;
      }

      // Commit tile accumulators
#ifdef use_SPFP
      fx_i += fast_llrintf(TLx_i * FORCESCALEF);
      fy_i += fast_llrintf(TLy_i * FORCESCALEF);
      fz_i += fast_llrintf(TLz_i * FORCESCALEF);
#else
#  ifdef PME_MINIMIZATION
      PMEFloat i;
      TLx_i = modf(TLx_i, &i);
      fxe_i += llrint(i);
      TLy_i = modf(TLy_i, &i);
      fye_i += llrint(i);
      TLz_i = modf(TLz_i, &i);
      fze_i += llrint(i);
#  endif
      fx_i += llrint(TLx_i * FORCESCALE);
      fy_i += llrint(TLy_i * FORCESCALE);
      fz_i += llrint(TLz_i * FORCESCALE);
#endif

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
#ifdef PHMD
        shAtom.indexj = atom;
        shAtom.g = cSim.pImageGrplist[atom] - 1;
        shAtom.qg = 0;
        shAtom.qxg = 0;
        shAtom.lgtitr = false;
        shAtom.qj_phmd = cSim.pImageCharge_phmd[atom];
        if (shAtom.g >= 0) {
          pqstate1 = cSim.pImageQstate1[atom];
          pqstate2 = cSim.pImageQstate2[atom];
          shAtom.psp_grp_g = cSim.psp_grp[shAtom.g];
          vstate1 = cSim.pvstate1[cSim.pImageIndexPHMD[atom]];
          radh = vstate1.x;
          lambda = sin(cSim.pph_theta[shAtom.g]);
          lambda *= lambda;
          shAtom.factg = 1.0 - lambda;
          shAtom.xg = 1.0;
          if (shAtom.psp_grp_g > 0) {
            shAtom.lgtauto = true;
            x2 = sin(cSim.pph_theta[shAtom.g+1]);
            x2 *= x2;
            radh += vstate1.y;
            if (vstate1.x > 0) {
              shAtom.igtaut = 1.0;
              shAtom.xg = x2;
            }
            else if (vstate1.y > 0) {
              shAtom.igtaut = -1.0;
              shAtom.xg = 1.0 - x2;
            }
            if (shAtom.psp_grp_g == 1) {
              shAtom.factg = 1.0 - lambda * shAtom.xg;
              shAtom.dxg = -shAtom.igtaut * lambda;
            }
            else if (shAtom.psp_grp_g == 3) {
              shAtom.factg = (1.0 - lambda) * shAtom.xg;
              shAtom.dxg = shAtom.igtaut * (1.0 - lambda);
            }
            qunprot = lambda * (pqstate2.x - pqstate2.y);
            qprot = (1 - lambda) * (pqstate1.x - pqstate1.y);
            shAtom.qxg = qunprot + qprot;
          }
          shAtom.lgtitr = (radh > 0);
          qunprot = x2 * pqstate2.x + (1.0 - x2) * pqstate2.y;
          qprot = x2 * pqstate1.x + (1.0 - x2) * pqstate1.y;
          shAtom.qg = qunprot - qprot;
        }
#endif
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
      // Initialize tile-specific accumulators
      PMEFloat TLx_i  = (PMEFloat)0.0;
      PMEFloat TLy_i  = (PMEFloat)0.0;
      PMEFloat TLz_i  = (PMEFloat)0.0;

      // Initialize iterators: j indicates that this thread should seek out the
      // "j atom" atom owned by thread j, whereas jrec indicates that thread jrec
      // will have information to contribute to the "j atom" owned by this thread.
      int j = tgx;
      int jrec = tgx;
      if (__ANY(WARP_MASK, exclusion)) {
        #pragma unroll 2
        for (int it = 0; it < PME_ATOMS_PER_WARP; it++) {
          PMEFloat xij = xi - __SHFL(WARP_MASK, shAtom.x, j);
          PMEFloat yij = yi - __SHFL(WARP_MASK, shAtom.y, j);
          PMEFloat zij = zi - __SHFL(WARP_MASK, shAtom.z, j);
          PMEFloat r2  = xij*xij + yij*yij + zij*zij;
          unsigned int index = LJIDi + __SHFL(WARP_MASK, shAtom.LJID, j);
          PMEFloat qiqj = qi * __SHFL(WARP_MASK, shAtom.q, j);
          PMEFloat df = (PMEFloat)0.0;
          bool inrange = ((r2 < cSim.cut2) && r2 > delta);
#ifdef PHMD
          int g = __SHFL(WARP_MASK, shAtom.g, j);
          PMEFloat qg = __SHFL(WARP_MASK, shAtom.qg, j);
          PMEFloat qxg = __SHFL(WARP_MASK, shAtom.qxg, j);
          int indexj = __SHFL(WARP_MASK, shAtom.indexj, j);
          qj_phmd = __SHFL(WARP_MASK, shAtom.qj_phmd, j);
          psp_grp_g = __SHFL(WARP_MASK, shAtom.psp_grp_g, j);
          xg = __SHFL(WARP_MASK, shAtom.xg, j);
          dxg = __SHFL(WARP_MASK, shAtom.dxg, j);
          factg = __SHFL(WARP_MASK, shAtom.factg, j);
          lgtauto = __SHFL(WARP_MASK, shAtom.lgtauto, j);
          lgtitr = __SHFL(WARP_MASK, shAtom.lgtitr, j);
#endif
          if (inrange) {
#if defined(use_SPFP) && !defined(PME_FSWITCH) && !defined(PME_ENERGY)
            PMEFloat r2inv = (PMEFloat)1.0 / r2;
            uint cidx = 2*(__float_as_uint(r2) >> 18) + (exclusion & 0x1);
#  if defined(__CUDA_ARCH__) && ((__CUDA_ARCH__ == 700) || (__CUDA_ARCH__ >= 800))
            PMEFloat4 coef = cSim.pErfcCoeffsTable[cidx];
#  else
            PMEFloat4 coef = tex1Dfetch<float4>(cSim.texErfcCoeffsTable, cidx);
#  endif
#  ifdef PHMD
          PMEFloat r     = sqrt(r2);
          PMEFloat rinv  = (PMEFloat)1.0 / r;
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
#  ifdef PHMD
#    ifdef use_DPFP
            PMEFloat swtch = erfc(cSim.ew_coeffSP * r);
#    else
            PMEFloat swtch = fasterfc(r);
#    endif
#  endif
#else
#  ifdef use_DPFP
            PMEFloat swtch = erfc(cSim.ew_coeffSP * r);
#  else
            PMEFloat swtch = fasterfc(r);
#  endif
            PMEFloat d_swtch_dx = cSim.negTwoEw_coeffRsqrtPI * exp(-cSim.ew_coeff2 * r2);
#endif
#if defined(PME_ENERGY) || defined(PHMD)
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
#  if defined(PME_ENERGY) || defined(PHMD)
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
#ifdef PHMD
              dfphmd = df6 - df12;
              de = fnrange * (f12 * twelvth - f6 * sixth);
#endif
#  ifdef PME_ENERGY
              TLevdw += fnrange * (f12*twelvth - f6*sixth);
#  endif
            }
#else  // PME_FSWITCH
            if (!(exclusion & 0x1)) {
              PMEFloat f6 = term.y * r6inv;
              PMEFloat f12 = term.x * r6inv * r6inv;
              df += (f12 - f6) * r2inv;
#ifdef PHMD
              dfphmd = (f12 - f6) * r2inv;
              de = fnrange * (f12 * twelvth - f6 * sixth);
#endif
#  ifdef PME_ENERGY
              TLevdw += fnrange * (f12*twelvth - f6*sixth);
#  endif
            }
#endif // PME_FSWITCH
#if defined(use_DPFP) || defined(PME_FSWITCH) || defined(PME_ENERGY) || defined(PHMD)
            else {
              swtch -= (PMEFloat)1.0;
            }
#endif
#ifdef PHMD
            //Second
            dudlj = 0;
            dudljplus = 0;
            factphmd = 1.0;
            if (!(exclusion & 0x1)) {
              if (lhtitr && lgtitr) {
                if (h != g) {
                  factphmd = facth * factg;
                  dudli -= xh * factg * de;
                  dudlj -= xg * facth * de;
                  if (lhtauto) {
                    dudliplus += dxh * factg * de;
                  }
                  if (lgtauto) {
                    dudljplus += dxg * facth * de;
                  }
                }
                else if (psp_grp == 1) {
                  factphmd = 1.0 - lambdah;
                  dudli -= de;
                }
              }
              else if (lhtitr) {
                factphmd = facth;
                dudli -= xh * de;
                if (lhtauto) {
                  dudliplus += dxh * de;
                }
              }
              else if (lgtitr) {
                factphmd = factg;
                dudlj -= xg * de;
                if (lgtauto) {
                  dudljplus += dxg * de;
                }
              }
              df += (factphmd - 1.0) * dfphmd;
#  ifdef PME_ENERGY
              TLevdw += (factphmd - 1.0) * de;
#   endif
          }
          if (h >= 0) {
            dudli += qj_phmd * qh * swtch * rinv * fnrange;
            if (psp_grp > 0) {
              dudliplus += qj_phmd * qxh * swtch * rinv * fnrange;
            }
          }
          if (g >= 0) {
            dudlj += qi_phmd * qg * swtch * rinv * fnrange;
            if (psp_grp_g > 0) {
              dudljplus += fnrange * qi_phmd * qxg * swtch * rinv;
            }
            atomicAdd((unsigned long long int*)&cSim.pdph_accumulator[indexj],
                      llitoulli(fast_llrintf(dudlj * FORCESCALEF)));
            atomicAdd((unsigned long long int*)&cSim.pdph_plus_accumulator[indexj],
                      llitoulli(fast_llrintf(dudljplus * FORCESCALEF)));
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
          shFx -= __SHFL(WARP_MASK, dfdx, jrec);
          shFy -= __SHFL(WARP_MASK, dfdy, jrec);
          shFz -= __SHFL(WARP_MASK, dfdz, jrec);
#ifdef PME_VIRIAL
          TLvir_11 -= xij * dfdx;
          TLvir_22 -= yij * dfdy;
          TLvir_33 -= zij * dfdz;
#endif
          exclusion >>= 1;
          j = ((j + 1) & cSim.NLAtomsPerWarpBitsMask) | joffset;
          jrec = ((jrec + PME_ATOMS_PER_WARP - 1) & cSim.NLAtomsPerWarpBitsMask) | joffset;
        }
        // End for loop covering non-bonded computations when
        // there IS at least one exclusion somewhere in the pile
      }
      else {
        #pragma unroll 2
        for (int it = 0; it < PME_ATOMS_PER_WARP; it++) {
          // Read properties for the other atom
          PMEFloat xij = xi - __SHFL(WARP_MASK, shAtom.x, j);
          PMEFloat yij = yi - __SHFL(WARP_MASK, shAtom.y, j);
          PMEFloat zij = zi - __SHFL(WARP_MASK, shAtom.z, j);
#ifdef PHMD
          int g = __SHFL(WARP_MASK, shAtom.g, j);
          PMEFloat qg = __SHFL(WARP_MASK, shAtom.qg, j);
          PMEFloat qxg = __SHFL(WARP_MASK, shAtom.qxg, j);
          int indexj = __SHFL(WARP_MASK, shAtom.indexj, j);
          qj_phmd = __SHFL(WARP_MASK, shAtom.qj_phmd, j);
          psp_grp_g = __SHFL(WARP_MASK, shAtom.psp_grp_g, j);
          xg = __SHFL(WARP_MASK, shAtom.xg, j);
          dxg = __SHFL(WARP_MASK, shAtom.dxg, j);
          factg = __SHFL(WARP_MASK, shAtom.factg, j);
          lgtauto = __SHFL(WARP_MASK, shAtom.lgtauto, j);
          lgtitr = __SHFL(WARP_MASK, shAtom.lgtitr, j);
#endif
          // Perform the range test
          PMEFloat r2  = xij*xij + yij*yij + zij*zij;
          unsigned int index = LJIDi + __SHFL(WARP_MASK, shAtom.LJID, j);
          PMEFloat qiqj = qi * __SHFL(WARP_MASK, shAtom.q, j);
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
#  ifdef PHMD
          PMEFloat r     = sqrt(r2);
          PMEFloat rinv  = (PMEFloat)1.0 / r;
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
#  ifdef PHMD
#    ifdef use_DPFP
            PMEFloat swtch = erfc(cSim.ew_coeffSP * r) * rinv;
#    else
            PMEFloat swtch = fasterfc(r) * rinv;
#    endif
#  endif
#else
#  ifdef use_DPFP
            PMEFloat swtch = erfc(cSim.ew_coeffSP * r) * rinv;
#  else
            PMEFloat swtch = fasterfc(r) * rinv;
#  endif
            PMEFloat d_swtch_dx = cSim.negTwoEw_coeffRsqrtPI * exp(-cSim.ew_coeff2 * r2);
#endif
#if defined(PME_ENERGY) || defined(PHMD)
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
#  if defined(PME_ENERGY) || defined(PHMD)
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
#  ifdef PHMD
            //Third
          dudlj = 0;
          dudljplus = 0;
          de = fnrange * (f12 * twelvth - f6 * sixth);
          factphmd = 1.0;
          if (lhtitr && lgtitr) {
            if (h != g) {
              factphmd = facth * factg;
              dudli -= xh * factg * de;
              dudlj -= xg * facth * de;
              if (lhtauto) {
                dudliplus += dxh * factg * de;
              }
              if (lgtauto) {
                dudljplus += dxg * facth * de;
              }
            }
            else if (psp_grp == 1) {
              factphmd = 1.0 - lambdah;
              dudli -= de;
            }
          }
          else if (lhtitr) {
            factphmd = facth;
            dudli -= xh * de;
            if (lhtauto) {
              dudliplus += dxh * de;
            }
          }
          else if (lgtitr) {
            factphmd = factg;
            dudlj -= xg * de;
            if (lgtauto) {
              dudljplus += dxg * de;
            }
          }
          df += (factphmd - 1.0) * (df6 - df12);
#  ifdef PME_ENERGY
          TLevdw += (factphmd - 1.0) * de;
#   endif
          if (h >= 0) {
            dudli += qj_phmd * qh * swtch;
            if (psp_grp > 0) {
              dudliplus += qj_phmd * qxh * swtch;
            }
          }
          if (g >= 0) {
            dudlj += qi_phmd * qg * swtch;
            if (psp_grp_g > 0) {
              dudljplus += qi_phmd * qxg * swtch;
            }
            atomicAdd((unsigned long long int*)&cSim.pdph_accumulator[indexj],
                      llitoulli(fast_llrintf(dudlj * FORCESCALEF)));
            atomicAdd((unsigned long long int*)&cSim.pdph_plus_accumulator[indexj],
                      llitoulli(fast_llrintf(dudljplus * FORCESCALEF)));
          }
#  endif
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
#ifdef PHMD
            //Fourth
            dudlj = 0;
            dudljplus = 0;
            de = fnrange * (f12 * twelvth - f6 * sixth);
            factphmd = 1.0;
            if (lhtitr && lgtitr) {
              if (h != g) {
                factphmd = facth * factg;
                dudli -= xh * factg * de;
                dudlj -= xg * facth * de;
                if (lhtauto) {
                  dudliplus += dxh * factg * de;
                }
                if (lgtauto) {
                  dudljplus += dxg * facth * de;
                }
              }
              else if (psp_grp == 1) {
                factphmd = 1.0 - lambdah;
                dudli -= de;
              }
            }
              else if (lhtitr) {
                factphmd = facth;
                dudli -= xh * de;
                if (lhtauto) {
                  dudliplus += dxh * de;
                }
              }
              else if (lgtitr) {
                factphmd = factg;
                dudlj -= xg * de;
                if (lgtauto) {
                  dudljplus += dxg * de;
                }
              }
              df += (factphmd - 1.0) * (f12 - f6) * r2inv;
#  ifdef PME_ENERGY
              TLevdw += (factphmd - 1.0) * de;
#   endif
            if (h >= 0) {
              dudli += fnrange * qj_phmd * qh * swtch;
              if (psp_grp > 0) {
                dudliplus += fnrange * qj_phmd * qxh * swtch;
              }
            }
            if (g >= 0) {
              dudlj += fnrange * qi_phmd * qg * swtch;
              if (psp_grp_g > 0) {
                dudljplus += fnrange * qi_phmd * qxg * swtch;
              }
              atomicAdd((unsigned long long int*)&cSim.pdph_accumulator[indexj],
                        llitoulli(fast_llrintf(dudlj * FORCESCALEF)));
              atomicAdd((unsigned long long int*)&cSim.pdph_plus_accumulator[indexj],
                        llitoulli(fast_llrintf(dudljplus * FORCESCALEF)));
            }
#endif
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
          shFx -= __SHFL(WARP_MASK, dfdx, jrec);
          shFy -= __SHFL(WARP_MASK, dfdy, jrec);
          shFz -= __SHFL(WARP_MASK, dfdz, jrec);
#ifdef PME_VIRIAL
          TLvir_11 -= xij * dfdx;
          TLvir_22 -= yij * dfdy;
          TLvir_33 -= zij * dfdz;
#endif
          j = ((j + 1) & cSim.NLAtomsPerWarpBitsMask) | joffset;
          jrec = ((jrec + PME_ATOMS_PER_WARP - 1) & cSim.NLAtomsPerWarpBitsMask) | joffset;
        }
        // Ends for loop for processing a pile of non-bonded
        // interactions with no exclusions to worry about
      }
      // End branch dealing with the presence of exclusions
      // in the pile of non-bonded interactions.

#ifdef use_SPFP
      // Commit tile accumulators to registers
      fx_i += fast_llrintf(TLx_i * FORCESCALEF);
      fy_i += fast_llrintf(TLy_i * FORCESCALEF);
      fz_i += fast_llrintf(TLz_i * FORCESCALEF);
#else
#  ifdef PME_MINIMIZATION
      PMEFloat i;
      TLx_i = modf(TLx_i, &i);
      fxe_i += llrint(i);
      TLy_i = modf(TLy_i, &i);
      fye_i += llrint(i);
      TLz_i = modf(TLz_i, &i);
      fze_i += llrint(i);
#  endif
      fx_i += llrint(TLx_i * FORCESCALE);
      fy_i += llrint(TLy_i * FORCESCALE);
      fz_i += llrint(TLz_i * FORCESCALE);
#endif

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
        PMEFloat i;
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
      fx_i += __SHFL(WARP_MASK, fx_i, tgx + rlev);
      fy_i += __SHFL(WARP_MASK, fy_i, tgx + rlev);
      fz_i += __SHFL(WARP_MASK, fz_i, tgx + rlev);
#if defined(use_DPFP) && defined(PME_MINIMIZATION)
      fxe_i += __SHFL(WARP_MASK, fxe_i, tgx + rlev);
      fye_i += __SHFL(WARP_MASK, fye_i, tgx + rlev);
      fze_i += __SHFL(WARP_MASK, fze_i, tgx + rlev);
#endif
    }

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
#ifdef PHMD
      if (h >= 0) {
        atomicAdd((unsigned long long int*)&cSim.pdph_accumulator[indexi],
                   llitoulli(fast_llrintf(dudli * FORCESCALEF)));
        atomicAdd((unsigned long long int*)&cSim.pdph_plus_accumulator[indexi],
                   llitoulli(fast_llrintf(dudliplus * FORCESCALEF)));
      }
#endif

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
#undef VOLATILE
}
