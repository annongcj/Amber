#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// This file is included repeatedly in kCalculatePMENonbondEnergy.cu.
// #defines: IPS_ENERGY, IPS_VIRIAL, IPS_IS_ORTHOGONAL, IPS_ATOMS_PER_WARP, IPS_MINIMIZATION
//---------------------------------------------------------------------------------------------
{
#if !defined(AMBER_PLATFORM_AMD)
#  define VOLATILE volatile
#else
#  define VOLATILE
#endif
  // Neighbor list atom
  struct NLAtom {
    PMEFloat x;
    PMEFloat y;
    PMEFloat z;
    PMEFloat q;
    unsigned int LJID;
    unsigned int ID;
  };

  // Neighbor list force on any one atom in the neighbor list of another
  struct NLForce {
    PMEForce x;
    PMEForce y;
    PMEForce z;
  };

  // Neighbor list virial
  struct NLVirial {
    long long int vir_11;
    long long int vir_22;
    long long int vir_33;
  };

  // Neighbor list warp to group atoms for import and processing
  struct NLWarp {
    NLEntry nlEntry;
    uint pos;
    bool bHomeCell;
  };

#define PSATOMX(i) shAtom.x
#define PSATOMY(i) shAtom.y
#define PSATOMZ(i) shAtom.z
#define PSATOMQ(i) shAtom.q
#define PSATOMLJID(i) shAtom.LJID
#define PSATOMID(i) shAtom.ID

#if defined(IPS_VIRIAL) || defined(IPS_ENERGY)
  const int THREADS_PER_BLOCK = IPSNONBONDENERGY_THREADS_PER_BLOCK;
#else
  const int THREADS_PER_BLOCK = IPSNONBONDFORCES_THREADS_PER_BLOCK;
#endif

  __shared__ unsigned int sNLEntries;
#if defined(IPS_VIRIAL)
  __shared__ PMEFloat sUcellf[9];
#endif
  __shared__ unsigned int sNext[GRID];
  __shared__ VOLATILE NLWarp sWarp[THREADS_PER_BLOCK / GRID];
  NLAtom shAtom;
#ifdef IPS_VIRIAL
  __shared__ NLVirial sWarpVirial[THREADS_PER_BLOCK / GRID];
#endif

  // Read static data
  if (threadIdx.x == 0) {
    sNLEntries = *(cSim.pNLEntries);
  }
  if (threadIdx.x < GRID) {
    unsigned int offset = cSim.NLAtomsPerWarp * (threadIdx.x >> cSim.NLAtomsPerWarpBits);
    sNext[threadIdx.x] = ((threadIdx.x + 1) & cSim.NLAtomsPerWarpBitsMask) + offset;
  }
#if (IPS_ATOMS_PER_WARP <= 16)
  const int iftlim = (IPS_ATOMS_PER_WARP * IPS_ATOMS_PER_WARP + GRID - 1) / GRID;
#endif

#ifdef IPS_VIRIAL
  if (threadIdx.x < 9) {
    sUcellf[threadIdx.x] = cSim.pNTPData->ucellf[threadIdx.x];
  }
  if (threadIdx.x < THREADS_PER_BLOCK / GRID) {
    sWarpVirial[threadIdx.x].vir_11 = 0;
    sWarpVirial[threadIdx.x].vir_22 = 0;
    sWarpVirial[threadIdx.x].vir_33 = 0;
  }
#endif

#ifdef IPS_ENERGY
  PMEForce eed  = (PMEForce)0;
  PMEForce evdw = (PMEForce)0;
#endif
  unsigned int tgx = threadIdx.x & (GRID - 1);
  unsigned int warp = THREADS_PER_BLOCK == GRID ? 0 : threadIdx.x / GRID;

  unsigned int jgroup = tgx >> cSim.NLAtomsPerWarpBits;
  unsigned int joffset = cSim.NLAtomsPerWarp * jgroup;
#if (IPS_ATOMS_PER_WARP <= 16)
  unsigned int jj = iftlim * jgroup;
#endif

  VOLATILE NLWarp* psWarp = &sWarp[warp];
  if (tgx == 0)
    psWarp->pos = blockIdx.x * (THREADS_PER_BLOCK / GRID) + warp;
  __syncthreads();

  // Step over all neighbor list entries, advancing one tile at a time
  while (psWarp->pos < sNLEntries) {
#ifdef IPS_VIRIAL
    PMEVirial vir_11 = (PMEVirial)0;
    PMEVirial vir_22 = (PMEVirial)0;
    PMEVirial vir_33 = (PMEVirial)0;
#endif

    // Read Neighbor List entry
    NLForce psF;
    if (tgx < 4) {
      psWarp->nlEntry.array[tgx] = cSim.pNLEntry[psWarp->pos].array[tgx];
    }
    if (tgx == 0) {
      psWarp->bHomeCell = psWarp->nlEntry.NL.ymax & NLENTRY_HOME_CELL_MASK;
      psWarp->nlEntry.NL.ymax >>= NLENTRY_YMAX_SHIFT;
    }
    __SYNCWARP(WARP_MASK);
    unsigned int offset = psWarp->nlEntry.NL.offset;

    // Read y atoms into registers
    PMEFloat xi;
    PMEFloat yi;
    PMEFloat zi;
    PMEFloat qi;
    unsigned int LJIDi;
    PMEForce fx_i = (PMEForce)0;
    PMEForce fy_i = (PMEForce)0;
    PMEForce fz_i = (PMEForce)0;
    unsigned int index = psWarp->nlEntry.NL.ypos + (tgx & cSim.NLAtomsPerWarpBitsMask);
    if (index < psWarp->nlEntry.NL.ymax) {
      PMEFloat2 xy    = cSim.pAtomXYSP[index];
      PMEFloat2 qljid = cSim.pAtomChargeSPLJID[index];
      zi    = cSim.pAtomZSP[index];
      xi    = xy.x;
      yi    = xy.y;
      qi    = qljid.x;
#ifdef use_DPFP
      LJIDi = __double_as_longlong(qljid.y);
#else
      LJIDi = __float_as_uint(qljid.y);
#endif
    }
    else {
      xi    = (PMEFloat)10000.0 * index;
      yi    = (PMEFloat)10000.0 * index;
      zi    = (PMEFloat)10000.0 * index;
      qi    = (PMEFloat)0.0;
      LJIDi = 0;
    }

#ifndef IPS_IS_ORTHOGONAL
    // Transform into cartesian space
#ifdef IPS_VIRIAL
    xi = xi*sUcellf[0] + yi*sUcellf[1] + zi*sUcellf[2];
    yi =                 yi*sUcellf[4] + zi*sUcellf[5];
    zi =                                 zi*sUcellf[8];
#else
    xi = xi*cSim.ucellf[0][0] + yi*cSim.ucellf[0][1] + zi*cSim.ucellf[0][2];
    yi =                        yi*cSim.ucellf[1][1] + zi*cSim.ucellf[1][2];
    zi =                                               zi*cSim.ucellf[2][2];
#endif
#endif
    // Special-case first tile
    // Copy register data into shared memory
    if (psWarp->bHomeCell) {
      unsigned int exclusion = cSim.pNLAtomList[offset + (tgx & cSim.NLAtomsPerWarpBitsMask)];
      offset += cSim.NLAtomsPerWarp;
#if (IPS_ATOMS_PER_WARP <= 16)
      // SIZE DEPENDENT 2, 8
      exclusion >>= jj;
#endif
      PSATOMX(tgx)     = xi;
      PSATOMY(tgx)     = yi;
      PSATOMZ(tgx)     = zi;
      PSATOMQ(tgx)     = qi;
      PSATOMLJID(tgx)  = LJIDi;
      LJIDi           *= cSim.LJTypes;

      // Set up iteration counts
#if (IPS_ATOMS_PER_WARP == 32)
      unsigned int j = ((tgx + 1) & cSim.NLAtomsPerWarpBitsMask) + joffset;
      unsigned int shIdx = j;
      shAtom.x    = __SHFL(WARP_MASK, shAtom.x, j);
      shAtom.y    = __SHFL(WARP_MASK, shAtom.y, j);
      shAtom.z    = __SHFL(WARP_MASK, shAtom.z, j);
      shAtom.q    = __SHFL(WARP_MASK, shAtom.q, j);
      shAtom.LJID = __SHFL(WARP_MASK, shAtom.LJID, j);
      PMEMask mask1 = __BALLOT(WARP_MASK, j != tgx);
      while (j != tgx) {
#else
      unsigned int j = (tgx + 1 + jj) & cSim.NLAtomsPerWarpBitsMask;
      unsigned int jend = (tgx < GRID - cSim.NLAtomsPerWarp) ? 1 : 0;
      unsigned int end = (((tgx + iftlim + jend) & cSim.NLAtomsPerWarpBitsMask) +
                           jj) & cSim.NLAtomsPerWarpBitsMask;
      unsigned int shIdx = ((tgx + 1) & cSim.NLAtomsPerWarpBitsMask) + joffset;
      shAtom.x    = __SHFL(WARP_MASK, shAtom.x, j);
      shAtom.y    = __SHFL(WARP_MASK, shAtom.y, j);
      shAtom.z    = __SHFL(WARP_MASK, shAtom.z, j);
      shAtom.q    = __SHFL(WARP_MASK, shAtom.q, j);
      shAtom.LJID = __SHFL(WARP_MASK, shAtom.LJID, j);
      PMEMask mask1 = __BALLOT(WARP_MASK, j != end);
      while (j != end) {
#endif
        PMEFloat xij = xi - PSATOMX(j);
        PMEFloat yij = yi - PSATOMY(j);
        PMEFloat zij = zi - PSATOMZ(j);
        PMEFloat r2  = xij * xij + yij * yij + zij * zij;
        if (r2 < cSim.cut2) {
          PMEFloat rinv      = rsqrt(r2);
          PMEFloat r2inv     = rinv * rinv;
          unsigned int LJIDj = PSATOMLJID(j);
          unsigned int index = LJIDi + LJIDj;
#ifndef use_DPFP
          PMEFloat2 term = tex1Dfetch<float2>(cSim.texLJTerm, index);
#else
          PMEFloat2 term = cSim.pLJTerm[index];
#endif
          PMEFloat b0     = qi * PSATOMQ(j);
          PMEFloat r4     = r2 * r2;
          PMEFloat r6inv  = r2inv * r2inv * r2inv;
          PMEFloat r12inv = r6inv * r6inv;
          PMEFloat f6     = term.y;
          PMEFloat f12    = term.x;
          PMEFloat dpipse = -r2 * (cSim.bipse1 + r2*(cSim.bipse2 + r2*cSim.bipse3));
          PMEFloat dvcu   = -r2 * (cSim.bipsvc1 + r2*(cSim.bipsvc2 + r2*cSim.bipsvc3));
          PMEFloat dvau   = -r4 * (cSim.bipsva1 + r4*(cSim.bipsva2 + r4*cSim.bipsva3));
#ifdef IPS_ENERGY
          PMEFloat pipse  = cSim.aipse0 + r2*(cSim.aipse1 +
                                              r2*(cSim.aipse2 + r2*cSim.aipse3)) - cSim.pipsec;
          PMEFloat pvc    = cSim.aipsvc0 +
                            r2*(cSim.aipsvc1 + r2*(cSim.aipsvc2 + r2*cSim.aipsvc3)) -
                            cSim.pipsvcc;
          PMEFloat pva    = cSim.aipsva0 +
                            r4*(cSim.aipsva1 + r4*(cSim.aipsva2 + r4*cSim.aipsva3)) -
                            cSim.pipsvac;
#endif
          if (!(exclusion & 0x1)) {
            dpipse += rinv;
            dvcu   += (PMEFloat)6.0  * r6inv;
            dvau   += (PMEFloat)12.0 * r12inv;
#ifdef IPS_ENERGY
            pipse += rinv;
            pvc   += r6inv;
            pva   += r12inv;
#endif
          }
#ifdef IPS_ENERGY
#  ifndef use_DPFP
          evdw += fast_llrintf((PMEFloat)0.5 * ENERGYSCALEF * (f12 * pva - f6 * pvc));
          eed  += fast_llrintf((PMEFloat)0.5 * ENERGYSCALEF * b0 * pipse);
#  else
          evdw += (PMEFloat)0.5 * (f12 * pva - f6 * pvc);
          eed  += (PMEFloat)0.5 * b0 * pipse;
#  endif
#endif
          PMEFloat df = (f12 * dvau - f6 * dvcu + b0 * dpipse) * r2inv;
#ifndef use_DPFP
#  ifdef IPS_MINIMIZATION
          df = max(-10000.0f, min(df, 10000.0f));
#  endif
          df *= FORCESCALEF;
#  endif
          PMEFloat dfdx = df * xij;
          PMEFloat dfdy = df * yij;
          PMEFloat dfdz = df * zij;

          // Accumulate into registers only
#ifndef use_DPFP
          fx_i += fast_llrintf(dfdx);
          fy_i += fast_llrintf(dfdy);
          fz_i += fast_llrintf(dfdz);
#  ifdef IPS_VIRIAL
          vir_11 -= fast_llrintf((PMEFloat)0.5 * xij * dfdx);
          vir_22 -= fast_llrintf((PMEFloat)0.5 * yij * dfdy);
          vir_33 -= fast_llrintf((PMEFloat)0.5 * zij * dfdz);
#  endif
#else
          fx_i += (PMEForce)dfdx;
          fy_i += (PMEForce)dfdy;
          fz_i += (PMEForce)dfdz;
#  ifdef IPS_VIRIAL
          vir_11 -= (PMEVirial)((PMEFloat)0.5 * xij * dfdx);
          vir_22 -= (PMEVirial)((PMEFloat)0.5 * yij * dfdy);
          vir_33 -= (PMEVirial)((PMEFloat)0.5 * zij * dfdz);
#  endif
#endif
        }

        // Shift the exclusion counter for the next atom
        exclusion >>= 1;
        shAtom.x    = __SHFL(mask1, shAtom.x, shIdx);
        shAtom.y    = __SHFL(mask1, shAtom.y, shIdx);
        shAtom.z    = __SHFL(mask1, shAtom.z, shIdx);
        shAtom.q    = __SHFL(mask1, shAtom.q, shIdx);
        shAtom.LJID = __SHFL(mask1, shAtom.LJID, shIdx);
        j = sNext[j];
#if (IPS_ATOMS_PER_WARP == 32)
        mask1 = __BALLOT(mask1, j != tgx);
      }
#else
        mask1 = __BALLOT(mask1, j != end);
      }
#endif
      // Here ends the while loop incrementing j, bounded
      // according to the number of IPS atoms per warp.
    }
    else {
      LJIDi *= cSim.LJTypes;
    }
    // End of branch for whether the atoms in question are part of the home cell

    // Handle remainder of line
    int tx = 0;
    while (tx < psWarp->nlEntry.NL.xatoms) {

      // Read atom ID and exclusion data
      PSATOMID(tgx) = cSim.pNLAtomList[offset + tgx];
      offset += GRID;
      PMEMask fullExclusion =
        ((PMEMask*)&cSim.pNLAtomList[offset])[tgx & cSim.NLAtomsPerWarpBitsMask];
      offset += cSim.NLAtomsPerWarp * sizeof(PMEMask) / sizeof(unsigned int);
#if (IPS_ATOMS_PER_WARP < 32)
      fullExclusion >>= cSim.NLAtomsPerWarp * (tgx >> cSim.NLAtomsPerWarpBits);
#endif
      unsigned int exclusion = (unsigned int)fullExclusion & cSim.NLAtomsPerWarpMask;

      // Clear shared memory forces
      psF.x = (PMEForce)0;
      psF.y = (PMEForce)0;
      psF.z = (PMEForce)0;

      // Read shared memory data
      if (tx + tgx < psWarp->nlEntry.NL.xatoms) {
        unsigned int atom = PSATOMID(tgx) >> NLATOM_CELL_SHIFT;
#ifndef use_DPFP
        PMEFloat2 xy    = tex1Dfetch<float2>(cSim.texAtomXYSP, atom);
        PMEFloat2 qljid = tex1Dfetch<float2>(cSim.texAtomChargeSPLJID, atom);
        PSATOMZ(tgx)    = tex1Dfetch<float>(cSim.texAtomZSP, atom);
#else
        PMEFloat2 xy    = cSim.pAtomXYSP[atom];
        PMEFloat2 qljid = cSim.pAtomChargeSPLJID[atom];
        PSATOMZ(tgx)    = cSim.pAtomZSP[atom];
#endif
        PSATOMX(tgx)    = xy.x;
        PSATOMY(tgx)    = xy.y;
        PSATOMQ(tgx)    = qljid.x;
#ifdef use_DPFP
        PSATOMLJID(tgx) = __double_as_longlong(qljid.y);
#else
        PSATOMLJID(tgx) = __float_as_uint(qljid.y);
#endif
      }
      else {
        PSATOMX(tgx)    = (PMEFloat)-10000.0 * tgx;
        PSATOMY(tgx)    = (PMEFloat)-10000.0 * tgx;
        PSATOMZ(tgx)    = (PMEFloat)-10000.0 * tgx;
        PSATOMQ(tgx)    = (PMEFloat)0.0;
        PSATOMLJID(tgx) = 0;
      }

      // Translate all atoms into a local coordinate system within one unit
      // cell of the first atom read to avoid PBC handling within inner loops
      int cell = PSATOMID(tgx) & NLATOM_CELL_TYPE_MASK;
#if defined(IPS_VIRIAL) && defined(IPS_IS_ORTHOGONAL)
      PSATOMX(tgx) += sUcellf[0] * cSim.cellOffset[cell][0];
      PSATOMY(tgx) += sUcellf[4] * cSim.cellOffset[cell][1];
      PSATOMZ(tgx) += sUcellf[8] * cSim.cellOffset[cell][2];
#else
      PSATOMX(tgx) += cSim.cellOffset[cell][0];
      PSATOMY(tgx) += cSim.cellOffset[cell][1];
      PSATOMZ(tgx) += cSim.cellOffset[cell][2];
#endif

#ifndef IPS_IS_ORTHOGONAL
#  ifdef IPS_VIRIAL
      PSATOMX(tgx) = sUcellf[0]*PSATOMX(tgx) + sUcellf[1]*PSATOMY(tgx) +
                     sUcellf[2]*PSATOMZ(tgx);
      PSATOMY(tgx) = sUcellf[4]*PSATOMY(tgx) + sUcellf[5]*PSATOMZ(tgx);
      PSATOMZ(tgx) = sUcellf[8]*PSATOMZ(tgx);
#  else
      PSATOMX(tgx) = cSim.ucellf[0][0]*PSATOMX(tgx) + cSim.ucellf[0][1]*PSATOMY(tgx) +
                     cSim.ucellf[0][2]*PSATOMZ(tgx);
      PSATOMY(tgx) = cSim.ucellf[1][1]*PSATOMY(tgx) + cSim.ucellf[1][2]*PSATOMZ(tgx);
      PSATOMZ(tgx) = cSim.ucellf[2][2]*PSATOMZ(tgx);
#  endif
#endif
      int j = tgx;
      int shIdx = ((tgx + 1) & cSim.NLAtomsPerWarpBitsMask) + joffset;

      // Branch the inner loop based on whether exclusions are present
      // to avoid having to check each interaction unnecessarily.
      if (__ANY(WARP_MASK, exclusion != 0)) {
        PMEMask mask1 = WARP_MASK;
        #pragma unroll 2
        for (unsigned int ii = 0; ii < IPS_ATOMS_PER_WARP; ++ii) {
          PMEFloat xij = xi - PSATOMX(j);
          PMEFloat yij = yi - PSATOMY(j);
          PMEFloat zij = zi - PSATOMZ(j);
          PMEFloat r2  = xij*xij + yij*yij + zij*zij;
          if (r2 < cSim.cut2) {
            PMEFloat rinv      = rsqrt(r2);
            PMEFloat r2inv     = rinv * rinv;
            unsigned int LJIDj = PSATOMLJID(j);
            unsigned int index = LJIDi + LJIDj;
#ifndef use_DPFP
            PMEFloat2 term = tex1Dfetch<float2>(cSim.texLJTerm, index);
#else
            PMEFloat2 term = cSim.pLJTerm[index];
#endif
            PMEFloat f6     = term.y;
            PMEFloat f12    = term.x;
            PMEFloat b0     = qi * PSATOMQ(j);
            PMEFloat r4     = r2 * r2;
            PMEFloat r6inv  = r2inv * r2inv * r2inv;
            PMEFloat r12inv = r6inv * r6inv;
            PMEFloat dpipse = -r2 * (cSim.bipse1 + r2*(cSim.bipse2 + r2*cSim.bipse3));
            PMEFloat dvcu   = -r2 * (cSim.bipsvc1 + r2*(cSim.bipsvc2 + r2*cSim.bipsvc3));
            PMEFloat dvau   = -r4 * (cSim.bipsva1 + r4*(cSim.bipsva2 + r4*cSim.bipsva3));
#ifdef IPS_ENERGY
            PMEFloat pipse  = cSim.aipse0 +
                              r2*(cSim.aipse1 + r2*(cSim.aipse2 + r2*cSim.aipse3)) -
                              cSim.pipsec;
            PMEFloat pvc    = cSim.aipsvc0 +
                              r2*(cSim.aipsvc1 + r2*(cSim.aipsvc2 + r2*cSim.aipsvc3)) -
                              cSim.pipsvcc;
            PMEFloat pva    = cSim.aipsva0 +
                              r4*(cSim.aipsva1 + r4*(cSim.aipsva2 + r4*cSim.aipsva3)) -
                              cSim.pipsvac;
#endif
            if (!(exclusion & 0x1)) {
              dpipse += rinv;
              dvcu   += (PMEFloat)6.0  * r6inv;
              dvau   += (PMEFloat)12.0 * r12inv;
#ifdef IPS_ENERGY
              pipse += rinv;
              pvc   += r6inv;
              pva   += r12inv;
#endif
            }
#ifdef IPS_ENERGY
#  ifndef use_DPFP
            evdw += fast_llrintf(ENERGYSCALEF * (f12 * pva - f6 * pvc));
            eed  += fast_llrintf(ENERGYSCALEF * b0 * pipse);
#  else
            evdw += f12 * pva - f6 * pvc;
            eed  += b0 * pipse;
#  endif
#endif
            PMEFloat df = (f12*dvau - f6*dvcu + b0*dpipse) * r2inv;
#ifndef use_DPFP
#  ifdef IPS_MINIMIZATION
            df = max(-10000.0f, min(df, 10000.0f));
#  endif
            df *= FORCESCALEF;
#endif
            PMEFloat dfdx = df * xij;
            PMEFloat dfdy = df * yij;
            PMEFloat dfdz = df * zij;
#ifndef use_DPFP
            long long int dfdx1 = fast_llrintf(dfdx);
            long long int dfdy1 = fast_llrintf(dfdy);
            long long int dfdz1 = fast_llrintf(dfdz);
            fx_i += dfdx1;
            fy_i += dfdy1;
            fz_i += dfdz1;
            psF.x -= dfdx1;
            psF.y -= dfdy1;
            psF.z -= dfdz1;
#  ifdef IPS_VIRIAL
            vir_11 -= fast_llrintf(xij * dfdx);
            vir_22 -= fast_llrintf(yij * dfdy);
            vir_33 -= fast_llrintf(zij * dfdz);
#  endif
#else
            PMEForce dfdx1 = (PMEForce)dfdx;
            PMEForce dfdy1 = (PMEForce)dfdy;
            PMEForce dfdz1 = (PMEForce)dfdz;
            fx_i += dfdx1;
            fy_i += dfdy1;
            fz_i += dfdz1;
            psF.x -= dfdx1;
            psF.y -= dfdy1;
            psF.z -= dfdz1;
#  ifdef IPS_VIRIAL
            vir_11 -= (PMEForce)(xij * dfdx);
            vir_22 -= (PMEForce)(yij * dfdy);
            vir_33 -= (PMEForce)(zij * dfdz);
#  endif
#endif
          }

          // Shift the exclusion counter to prepare for the next atom
          exclusion >>= 1;
          shAtom.x    = __SHFL(mask1, shAtom.x, shIdx);
          shAtom.y    = __SHFL(mask1, shAtom.y, shIdx);
          shAtom.z    = __SHFL(mask1, shAtom.z, shIdx);
          shAtom.q    = __SHFL(mask1, shAtom.q, shIdx);
          shAtom.LJID = __SHFL(mask1, shAtom.LJID, shIdx);
          psF.x       = __SHFL(mask1, psF.x, shIdx);
          psF.y       = __SHFL(mask1, psF.y, shIdx);
          psF.z       = __SHFL(mask1, psF.z, shIdx);
          j = ((j + 1) & cSim.NLAtomsPerWarpBitsMask) + joffset;
      }
        // Here ends the inner loop for the case that exclusions ARE present.
      }
      else {
        PMEMask mask1 = WARP_MASK;
        #pragma unroll 2
        for (unsigned int ii = 0; ii < IPS_ATOMS_PER_WARP; ++ii) {
          PMEFloat xij = xi - PSATOMX(j);
          PMEFloat yij = yi - PSATOMY(j);
          PMEFloat zij = zi - PSATOMZ(j);
          PMEFloat r2  = xij*xij + yij*yij + zij*zij;
          if (r2 < cSim.cut2) {
            PMEFloat rinv  = rsqrt(r2);
            PMEFloat r2inv = rinv * rinv;
            unsigned int LJIDj = PSATOMLJID(j);
            unsigned int index = LJIDi + LJIDj;
#ifndef use_DPFP
            PMEFloat2 term = tex1Dfetch<float2>(cSim.texLJTerm, index);
#else
            PMEFloat2 term = cSim.pLJTerm[index];
#endif
            PMEFloat f6     = term.y;
            PMEFloat f12    = term.x;
            PMEFloat b0     = qi * PSATOMQ(j);
            PMEFloat r4     = r2 * r2;
            PMEFloat r6inv  = r2inv * r2inv * r2inv;
            PMEFloat r12inv = r6inv * r6inv;
            PMEFloat dpipse = rinv - r2*(cSim.bipse1 + r2*(cSim.bipse2 + r2*cSim.bipse3));
            PMEFloat dvcu = (PMEFloat)6.0 * r6inv -
                            r2*(cSim.bipsvc1 + r2*(cSim.bipsvc2 + r2 * cSim.bipsvc3));
            PMEFloat dvau = (PMEFloat)12.0 * r12inv -
                            r4*(cSim.bipsva1 + r4*(cSim.bipsva2 + r4*cSim.bipsva3));
            PMEFloat df = (f12*dvau - f6*dvcu + b0*dpipse) * r2inv;
#ifdef IPS_ENERGY
            PMEFloat pipse = rinv + cSim.aipse0 +
                             r2*(cSim.aipse1 + r2*(cSim.aipse2 + r2*cSim.aipse3)) -
                             cSim.pipsec;
            PMEFloat pvc   = r6inv + cSim.aipsvc0 +
                             r2*(cSim.aipsvc1 + r2*(cSim.aipsvc2 + r2*cSim.aipsvc3)) -
                             cSim.pipsvcc;
            PMEFloat pva   = r12inv + cSim.aipsva0 +
                             r4*(cSim.aipsva1 + r4*(cSim.aipsva2 + r4 * cSim.aipsva3)) -
                             cSim.pipsvac;
#  ifndef use_DPFP
            evdw += fast_llrintf(ENERGYSCALEF * (f12 * pva - f6 * pvc));
            eed  += fast_llrintf(ENERGYSCALEF * b0 * pipse);
#  else  // use_DPFP
            evdw += f12 * pva - f6 * pvc;
            eed  += b0 * pipse;
#  endif // use_DPFP
#endif
#ifndef use_DPFP
#  ifdef IPS_MINIMIZATION
            df = max(-10000.0f, min(df, 10000.0f));
#  endif
            df *= FORCESCALEF;
#endif
            PMEFloat dfdx = df * xij;
            PMEFloat dfdy = df * yij;
            PMEFloat dfdz = df * zij;
#ifndef use_DPFP
            long long int dfdx1 = fast_llrintf(dfdx);
            long long int dfdy1 = fast_llrintf(dfdy);
            long long int dfdz1 = fast_llrintf(dfdz);
            fx_i += dfdx1;
            fy_i += dfdy1;
            fz_i += dfdz1;
            psF.x -= dfdx1;
            psF.y -= dfdy1;
            psF.z -= dfdz1;
#  ifdef IPS_VIRIAL
            vir_11 -= fast_llrintf(xij * dfdx);
            vir_22 -= fast_llrintf(yij * dfdy);
            vir_33 -= fast_llrintf(zij * dfdz);
#  endif
#else
            PMEForce dfdx1 = (PMEForce)dfdx;
            PMEForce dfdy1 = (PMEForce)dfdy;
            PMEForce dfdz1 = (PMEForce)dfdz;
            fx_i += dfdx1;
            fy_i += dfdy1;
            fz_i += dfdz1;
            psF.x -= dfdx1;
            psF.y -= dfdy1;
            psF.z -= dfdz1;
#  ifdef IPS_VIRIAL
            vir_11 -= (PMEForce)(xij * dfdx);
            vir_22 -= (PMEForce)(yij * dfdy);
            vir_33 -= (PMEForce)(zij * dfdz);
#  endif
#endif
          }
          shAtom.x    = __SHFL(mask1, shAtom.x, shIdx);
          shAtom.y    = __SHFL(mask1, shAtom.y, shIdx);
          shAtom.z    = __SHFL(mask1, shAtom.z, shIdx);
          shAtom.q    = __SHFL(mask1, shAtom.q, shIdx);
          shAtom.LJID = __SHFL(mask1, shAtom.LJID, shIdx);
          psF.x       = __SHFL(mask1, psF.x, shIdx);
          psF.y       = __SHFL(mask1, psF.y, shIdx);
          psF.z       = __SHFL(mask1, psF.z, shIdx);
          j = ((j + 1) & cSim.NLAtomsPerWarpBitsMask) + joffset;
      }
        // Here ends the inner loop for the case that no exclusions will be encountered.
      }
      // End of branch for the innter loop based on the presence of non-bonded exclusions

      // Dump shared memory forces
      if (tx + tgx < psWarp->nlEntry.NL.xatoms) {
        int offset = (PSATOMID(tgx) >> NLATOM_CELL_SHIFT);
#ifndef use_DPFP
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset],
                  llitoulli(psF.x));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset],
                  llitoulli(psF.y));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset],
                  llitoulli(psF.z));
#else
#  ifdef PME_MINIMIZATION
        psF.x = max((PMEFloat)-10000.0, min(psF.x, (PMEFloat)10000.0));
        psF.y = max((PMEFloat)-10000.0, min(psF.y, (PMEFloat)10000.0));
        psF.z = max((PMEFloat)-10000.0, min(psF.z, (PMEFloat)10000.0));
#  endif
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset],
                  llitoulli(llrint(psF.x * FORCESCALE)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset],
                  llitoulli(llrint(psF.y * FORCESCALE)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset],
                  llitoulli(llrint(psF.z * FORCESCALE)));
#endif
      }

      // Advance to next x tile
      tx += GRID;
    }

    // Reduce register forces if necessary
#if (IPS_ATOMS_PER_WARP <= 16)
#  ifdef AMBER_PLATFORM_AMD_WARP64
    fx_i += __SHFL(WARP_MASK, fx_i, tgx + 32);
    fy_i += __SHFL(WARP_MASK, fy_i, tgx + 32);
    fz_i += __SHFL(WARP_MASK, fz_i, tgx + 32);
#  endif
    fx_i += __SHFL(WARP_MASK, fx_i, tgx + 16);
    fy_i += __SHFL(WARP_MASK, fy_i, tgx + 16);
    fz_i += __SHFL(WARP_MASK, fz_i, tgx + 16);
#  if (IPS_ATOMS_PER_WARP == 8)
    fx_i += __SHFL(WARP_MASK, fx_i, tgx + 8);
    fy_i += __SHFL(WARP_MASK, fy_i, tgx + 8);
    fz_i += __SHFL(WARP_MASK, fz_i, tgx + 8);
#  endif
#endif

    // Dump register forces
    if (psWarp->nlEntry.NL.ypos + tgx < psWarp->nlEntry.NL.ymax) {
      int offset = psWarp->nlEntry.NL.ypos + tgx;
#ifndef use_DPFP
      atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset], llitoulli(fx_i));
      atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset], llitoulli(fy_i));
      atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset], llitoulli(fz_i));
#else
#ifdef PME_MINIMIZATION
      fx_i = max((PMEFloat)-10000.0, min(fx_i, (PMEFloat)10000.0));
      fy_i = max((PMEFloat)-10000.0, min(fy_i, (PMEFloat)10000.0));
      fz_i = max((PMEFloat)-10000.0, min(fz_i, (PMEFloat)10000.0));
#endif
      atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset],
                llitoulli(llrint(fx_i * FORCESCALE)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset],
                llitoulli(llrint(fy_i * FORCESCALE)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset],
                llitoulli(llrint(fz_i * FORCESCALE)));
#endif
    }

#ifdef IPS_VIRIAL
    // Reduce virial per warp and convert to fixed point if necessary
    NLVirial* psV = &sWarpVirial[warp];
    for (unsigned int stride = warpSize >> 1; stride > 0; stride >>=1) {
      vir_11 += __SHFL_DOWN(WARP_MASK, vir_11, stride);
      vir_22 += __SHFL_DOWN(WARP_MASK, vir_22, stride);
      vir_33 += __SHFL_DOWN(WARP_MASK, vir_33, stride);
    }

    if (tgx == 0) {
#ifndef use_DPFP
      psV->vir_11 += vir_11;
      psV->vir_22 += vir_22;
      psV->vir_33 += vir_33;
#else
      psV->vir_11 += llrint(vir_11 * FORCESCALE);
      psV->vir_22 += llrint(vir_22 * FORCESCALE);
      psV->vir_33 += llrint(vir_33 * FORCESCALE);
#endif
    }
#endif

    // Get next Neighbor List entry
    if (tgx == 0) {
      psWarp->pos = atomicAdd(&cSim.pFrcBlkCounters[0], 1);
    }
    __SYNCWARP(WARP_MASK);
  }

#ifdef IPS_VIRIAL
  NLVirial* psV = &sWarpVirial[warp];
  if ((threadIdx.x & GRID_BITS_MASK) == 0) {
    unsigned long long int val1 = llitoulli(psV->vir_11);
    unsigned long long int val2 = llitoulli(psV->vir_22);
    unsigned long long int val3 = llitoulli(psV->vir_33);
    atomicAdd(cSim.pVirial_11, val1);
    atomicAdd(cSim.pVirial_22, val2);
    atomicAdd(cSim.pVirial_33, val3);
  }
#endif

#ifdef IPS_ENERGY
#define sEED(i)  eed
#define sEVDW(i) evdw
  for (unsigned int stride = warpSize >> 1; stride > 0; stride >>=1) {
    eed += __SHFL_DOWN(WARP_MASK, eed, stride);
    evdw += __SHFL_DOWN(WARP_MASK, evdw, stride);
  }
  if ((threadIdx.x & GRID_BITS_MASK) == 0) {
#ifndef use_DPFP
    atomicAdd(cSim.pEED, llitoulli(sEED(threadIdx.x)));
    atomicAdd(cSim.pEVDW, llitoulli(sEVDW(threadIdx.x)));
#else  // use_DPFP
    atomicAdd(cSim.pEED, llitoulli(llrint(ENERGYSCALE * sEED(threadIdx.x))));
    atomicAdd(cSim.pEVDW, llitoulli(llrint(ENERGYSCALE * sEVDW(threadIdx.x))));
#endif
  }
#undef sEED
#undef sEVDW
#endif
#undef PSATOMX
#undef PSATOMY
#undef PSATOMZ
#undef PSATOMQ
#undef PSATOMLJID
#undef PSATOMID
#undef VOLATILE
}
