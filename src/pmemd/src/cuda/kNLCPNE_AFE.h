#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// Inner loops for non-bonded interactions in the context of alchemical free energy problems.
// VDW force-switching code could use optimization.  This file is included in
// kCalculatePMENonbondEnergy.cu in the same way as the original kNLCPNE.h.
//
// #defines: PME_ENERGY, PME_VIRIAL, PME_IS_ORTHOGONAL, PME_ATOMS_PER_WARP, PME_MINIMIZATION,
//           PME_SCTI, PME_SCMBAR, AFE_VERBOSE
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
#ifdef PME_SCTI
    int TI;
#else
#endif
  };

  struct NLForce {
    PMEForce x;
    PMEForce y;
    PMEForce z;
#if defined(AFE_VERBOSE) && defined(PME_ENERGY)
    PMEForce vdwR1;
    PMEForce vdwR2;
    PMEForce eedR1;
    PMEForce eedR2;
#endif
  };

  struct NLVirial
  {
    long long int vir_11;
    long long int vir_22;
    long long int vir_33;
  };

  struct NLWarp
  {
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
#ifdef PME_SCTI
#  define PSATOMTI(i) shAtom.TI
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
  __shared__ unsigned int sNext[GRID];
#if (PME_ATOMS_PER_WARP <= 16)
  __shared__ unsigned int sStart[GRID];
  __shared__ unsigned int sEnd[GRID];
#endif
#ifdef PME_SCTI
  __shared__ PMEFloat sLambda[3];
#  ifdef PME_ENERGY
  __shared__ PMEFloat sTISigns[2];
#  endif
#endif

  __shared__ volatile NLWarp sWarp[THREADS_PER_BLOCK / GRID];
  __shared__ volatile NLForce sForce[THREADS_PER_BLOCK];
  NLAtom shAtom;
#ifdef PME_VIRIAL
  __shared__ volatile NLVirial sWarpVirial[THREADS_PER_BLOCK / GRID];
#endif

#ifdef PME_SCMBAR
  PMEForce bar_cont[16];
  for (int i = 0; i < 16; i++) {
    bar_cont[i] = 0;
  }
#endif

  // Read static data
  if (threadIdx.x == 0) {
    sNLEntries = *(cSim.pNLEntries);
  }
  if (threadIdx.x < GRID) {
    unsigned int offset = cSim.NLAtomsPerWarp * (threadIdx.x >> cSim.NLAtomsPerWarpBits);
    sNext[threadIdx.x] = ((threadIdx.x + 1) & cSim.NLAtomsPerWarpBitsMask) + offset;
  }

#if (PME_ATOMS_PER_WARP <= 16)
  const int iftlim = (PME_ATOMS_PER_WARP * PME_ATOMS_PER_WARP + GRID - 1) / GRID;
  if (threadIdx.x < GRID) {
    // SIZE DEPENDENT 2, 8
    sStart[threadIdx.x] = (threadIdx.x + 1 + iftlim*(threadIdx.x >> cSim.NLAtomsPerWarpBits)) &
                          cSim.NLAtomsPerWarpBitsMask;
    // SIZE DEPENDENT (3, 2, 2) (9, 8, 8)
    if (threadIdx.x < GRID - cSim.NLAtomsPerWarp) {
      sEnd[threadIdx.x] = (((threadIdx.x + iftlim + 1) & cSim.NLAtomsPerWarpBitsMask) +
                           iftlim*(threadIdx.x >> cSim.NLAtomsPerWarpBits)) &
                          cSim.NLAtomsPerWarpBitsMask;
    }
    else {
      sEnd[threadIdx.x] = (((threadIdx.x + iftlim) & cSim.NLAtomsPerWarpBitsMask) +
                           iftlim*(threadIdx.x >> cSim.NLAtomsPerWarpBits)) &
                          cSim.NLAtomsPerWarpBitsMask;
    }
  }
#endif

#ifdef PME_VIRIAL
  if (threadIdx.x < 9) {
    sUcellf[threadIdx.x] = cSim.pNTPData->ucellf[threadIdx.x];
  }
  if (threadIdx.x < THREADS_PER_BLOCK / GRID) {
    sWarpVirial[threadIdx.x].vir_11 = 0;
    sWarpVirial[threadIdx.x].vir_22 = 0;
    sWarpVirial[threadIdx.x].vir_33 = 0;
  }
#endif

#ifdef PME_SCTI
  //1, 1 - lambda, lambda
  sLambda[0] = (PMEFloat)1.0;
#  ifdef use_DPFP
  sLambda[1] = cSim.AFElambda[0];
  sLambda[2] = cSim.AFElambda[1];
#  else
  sLambda[1] = cSim.AFElambdaSP[0];
  sLambda[2] = cSim.AFElambdaSP[1];
#  endif
#  ifdef PME_ENERGY
  sTISigns[0] = cSim.TIsigns[0];
  sTISigns[1] = cSim.TIsigns[1];
#  endif
#endif

#ifdef PME_ENERGY
  PMEForce eed        = (PMEForce)0;
  PMEForce evdw       = (PMEForce)0;
#  ifdef PME_SCTI
  PMEForce sc_dvdl    = (PMEForce)0;
  PMEForce sc_evdw[2] = {(PMEForce)0, (PMEForce)0};
  PMEForce sc_eed[2]  = {(PMEForce)0, (PMEForce)0};
#  endif
#endif
  unsigned int tgx = threadIdx.x & (GRID - 1);
  volatile NLWarp* psWarp = &sWarp[threadIdx.x / GRID];
  if (tgx == 0)
    psWarp->pos = (blockIdx.x * blockDim.x + threadIdx.x) / GRID;
  __syncthreads();
  while (psWarp->pos < sNLEntries) {
#ifdef PME_VIRIAL
    PMEVirial vir_11 = (PMEVirial)0;
    PMEVirial vir_22 = (PMEVirial)0;
    PMEVirial vir_33 = (PMEVirial)0;
#endif

    // Read Neighbor List entry
    unsigned int tbx = threadIdx.x - tgx;
    volatile NLForce* psF = &sForce[tbx];
    if (tgx < 4) {
      psWarp->nlEntry.array[tgx] = cSim.pNLEntry[psWarp->pos].array[tgx];
    }
    if (tgx == 0) {
      psWarp->bHomeCell = psWarp->nlEntry.NL.ymax & NLENTRY_HOME_CELL_MASK;
      psWarp->nlEntry.NL.ymax >>= NLENTRY_YMAX_SHIFT;
    }
    __SYNCWARP(WARP_MASK);

    // Read y atoms into registers
    PMEFloat xi;
    PMEFloat yi;
    PMEFloat zi;
    PMEFloat qi;
#ifdef PME_SCTI
    unsigned int TIi;
#endif
    unsigned int LJIDi;
    PMEForce fx_i = (PMEForce)0;
    PMEForce fy_i = (PMEForce)0;
    PMEForce fz_i = (PMEForce)0;
    unsigned int index = psWarp->nlEntry.NL.ypos + (tgx & cSim.NLAtomsPerWarpBitsMask);
    if (index < psWarp->nlEntry.NL.ymax) {
      PMEFloat2 xy    = cSim.pAtomXYSP[index];
      PMEFloat2 qljid = cSim.pAtomChargeSPLJID[index];
      zi = cSim.pAtomZSP[index];
#ifdef PME_SCTI
      TIi = cSim.pImageTIRegion[index];
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
#ifdef PME_SCTI
      TIi = 0;
#endif
    }

#ifndef PME_IS_ORTHOGONAL
    // Transform into cartesian space
#  ifdef PME_VIRIAL
    xi     = xi*sUcellf[0] + yi*sUcellf[1] + zi*sUcellf[2];
    yi     =                 yi*sUcellf[4] + zi*sUcellf[5];
    zi     =                                 zi*sUcellf[8];
#  else
    xi     = xi*cSim.ucellf[0][0] + yi*cSim.ucellf[0][1] + zi*cSim.ucellf[0][2];
    yi     =                        yi*cSim.ucellf[1][1] + zi*cSim.ucellf[1][2];
    zi     =                                               zi*cSim.ucellf[2][2];
#  endif
#endif
    // Special-case first tile
    // Copy register data into shared memory
    if (psWarp->bHomeCell) {
      unsigned int exclusion = cSim.pNLAtomList[psWarp->nlEntry.NL.offset +
                                                (tgx & cSim.NLAtomsPerWarpBitsMask)];
      __SYNCWARP(WARP_MASK);
      if (tgx == 0)
        psWarp->nlEntry.NL.offset += cSim.NLAtomsPerWarp;
      __SYNCWARP(WARP_MASK);
#if (PME_ATOMS_PER_WARP <= 16)
      // SIZE DEPENDENT 2, 8
      exclusion >>= iftlim * (tgx >> cSim.NLAtomsPerWarpBits);
#endif
      PSATOMX(tgx)    = xi;
      PSATOMY(tgx)    = yi;
      PSATOMZ(tgx)    = zi;
      PSATOMQ(tgx)    = qi;
      PSATOMLJID(tgx) = LJIDi;
#ifdef PME_SCTI
      PSATOMTI(tgx)   = TIi;
#endif
      LJIDi          *= cSim.LJTypes;

      // Set up iteration counts
#if (PME_ATOMS_PER_WARP == 32)
      unsigned int j = sNext[tgx];
#else
      unsigned int j   = sStart[tgx];
      unsigned int end = sEnd[tgx];
#endif
      unsigned int shIdx = j;
      shAtom.x    = __SHFL(WARP_MASK, shAtom.x, j);
      shAtom.y    = __SHFL(WARP_MASK, shAtom.y, j);
      shAtom.z    = __SHFL(WARP_MASK, shAtom.z, j);
      shAtom.q    = __SHFL(WARP_MASK, shAtom.q, j);
      shAtom.LJID = __SHFL(WARP_MASK, shAtom.LJID, j);
#ifdef PME_SCTI
      shAtom.TI   = __SHFL(WARP_MASK, shAtom.TI, j);
#endif

#ifdef PME_SCTI
      if (__ALL(WARP_MASK, TIi == 0)) {
#endif
#if (PME_ATOMS_PER_WARP == 32)
        PMEMask mask1 = __BALLOT(WARP_MASK, j != tgx);
        while (j != tgx) {
#else
        PMEMask mask1 = __BALLOT(WARP_MASK, j != end);
        while (j != end) {
#endif
          PMEFloat xij = xi - PSATOMX(j);
          PMEFloat yij = yi - PSATOMY(j);
          PMEFloat zij = zi - PSATOMZ(j);
          PMEFloat r2  = xij * xij + yij * yij + zij * zij;
          if (r2 < cSim.cut2) {
            PMEFloat rinv      = rsqrt(r2);
            PMEFloat r         = r2 * rinv;
            PMEFloat r2inv     = rinv * rinv;
            PMEFloat r6inv     = r2inv * r2inv * r2inv;
            unsigned int LJIDj = PSATOMLJID(j);
            unsigned int index = LJIDi + LJIDj;
#ifndef use_DPFP
            PMEFloat2 term = tex1Dfetch<float2>(cSim.texLJTerm, index);
#else
            PMEFloat2 term = cSim.pLJTerm[index];
#endif
            PMEFloat df    = (PMEFloat)0.0;
            PMEFloat qiqj  = qi * PSATOMQ(j);
#ifdef use_DPFP
            PMEFloat swtch = erfc(cSim.ew_coeffSP * r);
#else
            PMEFloat swtch = fasterfc(r);
#endif
            PMEFloat d_swtch_dx = cSim.negTwoEw_coeffRsqrtPI * exp(-cSim.ew_coeff2 * r2);
            if (!(exclusion & 0x1)) {
#ifdef PME_FSWITCH
              PMEFloat r3inv  = rinv*r2inv;
              PMEFloat r12inv = r6inv*r6inv;
              PMEFloat df12f  = (PMEFloat)(-1) * term.x * cSim.cut6invcut6minfswitch6 * r2inv *
                                r6inv * (r6inv - cSim.cut6inv);
              PMEFloat df6f   = (PMEFloat)(-1) * term.y * cSim.cut3invcut3minfswitch3 * r3inv *
                                r2inv * (r3inv - cSim.cut3inv);
              PMEFloat df12   = (PMEFloat)(-1) * term.x * r2inv * r12inv;
              PMEFloat df6    = (PMEFloat)(-1) * term.y * r2inv * r6inv;
#  ifdef PME_ENERGY
              PMEFloat f12f   = term.x * cSim.cut6invcut6minfswitch6 * (r6inv - cSim.cut6inv) *
                                (r6inv - cSim.cut6inv);
              PMEFloat f6f    = term.y * cSim.cut3invcut3minfswitch3 * (r3inv - cSim.cut3inv) *
                                (r3inv - cSim.cut3inv);
              PMEFloat f12    = term.x*r12inv - term.x*cSim.invfswitch6cut6;
              PMEFloat f6     = term.y*r6inv - term.y*cSim.invfswitch3cut3;
#  endif
              if (r2 > cSim.fswitch2) {
                df12 = df12f;
                df6  = df6f;
#  ifdef PME_ENERGY
                f12  = f12f;
                f6   = f6f;
#  endif
              }
#  ifdef PME_ENERGY
#    ifndef use_DPFP
              evdw += fast_llrintf(ENERGYSCALEF * (f12*(PMEFloat)(0.5 / 12.0) -
                                                   f6*(PMEFloat)(0.5 / 6.0)));
#    else
              evdw += f12*(PMEFloat)(0.5 / 12.0) - f6*(PMEFloat)(0.5 / 6.0);
#    endif
#  endif
              df   += df6 - df12;
#else  //PME_FSWITCH
              PMEFloat f6  = term.y * r6inv;
              PMEFloat f12 = term.x * r6inv * r6inv;
              df += (f12 - f6) * r2inv;
#  ifdef PME_ENERGY
#    ifndef use_DPFP
              evdw += fast_llrintf(ENERGYSCALEF * (f12*(PMEFloat)(0.5 / 12.0) -
                                                   f6*(PMEFloat)(0.5 / 6.0)));
#    else
              evdw += (f12*(PMEFloat)(0.5 / 12.0) - f6*(PMEFloat)(0.5 / 6.0));
#    endif
#  endif

#endif //PME_FSWITCH
            }
            else {
              swtch -= (PMEFloat)1.0;
            }
#ifdef PME_ENERGY
            PMEFloat b0 = qiqj * swtch * rinv;
            PMEFloat b1 = b0 - qiqj * d_swtch_dx;
            df  += b1 * r2inv;
#  ifndef use_DPFP
            eed += fast_llrintf(b0 * (PMEFloat)0.5 * ENERGYSCALEF);
#  else
            eed += b0 * (PMEFloat)0.5;
#  endif
#else
            df  += qiqj * (swtch*rinv - d_swtch_dx) * r2inv;
#endif
#if !defined(use_DPFP) && defined(PME_MINIMIZATION)
            df   = max(-10000.0f, min(df, 10000.0f));
#endif
#ifdef use_SPFP
            df  *= FORCESCALEF;
#endif
            PMEFloat dfdx = df * xij;
            PMEFloat dfdy = df * yij;
            PMEFloat dfdz = df * zij;

            // Accumulate into registers only
#ifdef use_SPFP
            fx_i += fast_llrintf(dfdx);
            fy_i += fast_llrintf(dfdy);
            fz_i += fast_llrintf(dfdz);
#  ifdef PME_VIRIAL
            vir_11 -= fast_llrintf((PMEFloat)0.5 * xij * dfdx);
            vir_22 -= fast_llrintf((PMEFloat)0.5 * yij * dfdy);
            vir_33 -= fast_llrintf((PMEFloat)0.5 * zij * dfdz);
#  endif
#else  // use_DPFP
            fx_i += (PMEForce)dfdx;
            fy_i += (PMEForce)dfdy;
            fz_i += (PMEForce)dfdz;
#  ifdef PME_VIRIAL
            vir_11 -= (PMEVirial)((PMEFloat)0.5 * xij * dfdx);
            vir_22 -= (PMEVirial)((PMEFloat)0.5 * yij * dfdy);
            vir_33 -= (PMEVirial)((PMEFloat)0.5 * zij * dfdz);
#  endif
#endif
          }
          exclusion >>= 1;
          shAtom.x    = __SHFL(mask1, shAtom.x, shIdx);
          shAtom.y    = __SHFL(mask1, shAtom.y, shIdx);
          shAtom.z    = __SHFL(mask1, shAtom.z, shIdx);
          shAtom.q    = __SHFL(mask1, shAtom.q, shIdx);
          shAtom.LJID = __SHFL(mask1, shAtom.LJID, shIdx);
          j = sNext[j];
#if (PME_ATOMS_PER_WARP == 32)
          mask1 = __BALLOT(mask1, j != tgx);
        }
#else
          mask1 = __BALLOT(mask1, j != end);
        }
#endif
        // Here ends the while loop iterating over j, branched by the pre-processor
        // depending on whether PME_ATOMS_PER_WARP is 32 or not.  This all took
        // place under the condition that there are no TI interactions to consider
        // among this pile of non-bonded interactions.
#ifdef PME_SCTI
      }
      else {

        // In the alternative to the above contingency, there ARE TI interactions
        // to consider in this pile of non-bonded interactions.  The same
        // pre-processor differentiation of the inner loop occurs.
#if (PME_ATOMS_PER_WARP == 32)
        PMEMask mask1 = __BALLOT(WARP_MASK, j != tgx);
        while (j != tgx) {
#else
        PMEMask mask1 = __BALLOT(WARP_MASK, j != end);
        while (j != end) {
#endif
          unsigned int TIj = PSATOMTI(j);

          // Exclude all interactions across TI regions.  TIi can be 101, 011, 100,
          // 010, where anytime you're at 6, you're cross region (111 or 110)
          if ((TIi | TIj) < 6) {
            PMEFloat xij = xi - PSATOMX(j);
            PMEFloat yij = yi - PSATOMY(j);
            PMEFloat zij = zi - PSATOMZ(j);
            PMEFloat r2  = xij*xij + yij*yij + zij*zij;
            if (r2 < cSim.cut2) {
              PMEFloat rinv       = rsqrt(r2);
              PMEFloat r          = r2 * rinv;
              PMEFloat r2inv      = rinv * rinv;
              PMEFloat r6inv      = r2inv * r2inv * r2inv;
              unsigned int LJIDj  = PSATOMLJID(j);
              bool bSCtoCV        = ((TIi ^ TIj) & 1);
              bool bSCtoSC        = (TIi & TIj & 1);
              int region          = ((TIi | TIj) >> 1);
              PMEFloat lambda     = sLambda[region]; //1, lambda, 1-lambda
              unsigned int index  = LJIDi + LJIDj + (bSCtoCV * cSim.LJOffset);
#  ifndef use_DPFP
              PMEFloat2 term      = tex1Dfetch<float2>(cSim.texLJTerm, index);
#  else
              PMEFloat2 term      = cSim.pLJTerm[index];
#  endif
#  ifdef use_DPFP
              PMEFloat swtch      = erfc(cSim.ew_coeffSP * r);
#  else
              PMEFloat swtch      = fasterfc(r);
#  endif
              PMEFloat d_swtch_dx = cSim.negTwoEw_coeffRsqrtPI * exp(-cSim.ew_coeff2 * r2);
              PMEFloat qiqj       = qi * PSATOMQ(j);
              PMEFloat df         = (PMEFloat)0.0;
#  ifdef PME_ENERGY
              PMEFloat inv_swtch  = (PMEFloat)1.0 / swtch;

              // We need these for sc-c forces
              PMEFloat b0         = qiqj * swtch * rinv;
#  endif
              // Check for alpha/beta scaled SC to C/V interaction
              if (!bSCtoCV) {
                PMEFloat f6       = term.y * r6inv;
                PMEFloat f12      = term.x * r6inv * r6inv;
                if (!(exclusion & 0x1)) {
                  if (!bSCtoSC) {
                    df += lambda * (f12 - f6) * r2inv;
#  ifdef PME_SCMBAR
                    // MBAR VDW LINEAR
                    if (region != 0) {
                      for (int bar_i = 0; bar_i < cSim.bar_states; bar_i++) {
#    ifndef use_DPFP
                        bar_cont[bar_i] +=
                          fast_llrintf((PMEFloat)0.50 *
                                       (f12/(PMEFloat)12.0 - f6/(PMEFloat)6.0) *
                                       (cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i] -
                                        lambda) * ENERGYSCALEF);
#    else
                        bar_cont[bar_i] +=
                          (PMEFloat)0.50 * (f12/(PMEFloat)12.0 - f6/(PMEFloat)6.0) *
                          (cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i] - lambda);
#    endif
                      }
                    }
#  endif
                  }
                  if (bSCtoSC) {
                    df += (f12 - f6) * r2inv;
                  }
#  ifdef PME_ENERGY
#    ifndef use_DPFP
                  PMEFloat de = ENERGYSCALEF * (f12*(PMEFloat)(0.5 / 12.0) -
                                                f6*(PMEFloat)(0.5 / 6.0));

                  // c-c, need to count for both regions
                  if (region == 0) {
                    evdw += fast_llrintf(de);
                  }
                  if (region != 0) {
                    if (!bSCtoSC) {
                      sc_dvdl += fast_llrintf(sTISigns[region - 1] * de);
                      evdw    += fast_llrintf(lambda * de);
                    }
                    if (bSCtoSC) {
                      sc_evdw[region - 1] += fast_llrintf(de);
                    }
                  }

#    else  // use_DPFP
                  PMEFloat de = f12*(PMEFloat)(0.5 / 12.0) - f6*(PMEFloat)(0.5 / 6.0);
                  if (region == 0) {
                    evdw += de;
                  }
                  if (region != 0) {
                    if (!bSCtoSC) {
                      sc_dvdl += sTISigns[region - 1] * de;
                      evdw    += lambda * de;
                    }
                    if (bSCtoSC) {
                      sc_evdw[region - 1] += de;
                    }
                  }
#    endif // use_DPFP
#  endif // PME_ENERGY
                }
                else {
                  swtch -= (PMEFloat)1.0;
                }
                if (bSCtoSC) {
                  df += lambda * qiqj * r2inv * (swtch*rinv - d_swtch_dx);
                  if (swtch > 0) {
                    df += qiqj * r2inv * rinv * ((PMEFloat)1.0 - lambda);
#  ifdef PME_SCMBAR
                    // MBAR EELEC SC-SC
                    PMEFloat b0 = qiqj * swtch * rinv;
                    for (int bar_i = 0; bar_i < cSim.bar_states; bar_i++) {
#    ifndef use_DPFP
                      bar_cont[bar_i] +=
                        fast_llrintf((PMEFloat)0.5 * (b0 - rinv*qiqj) *
                                     (cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i] -
                                      lambda) * ENERGYSCALEF);
#    else
                      bar_cont[bar_i] +=
                        (PMEFloat)0.5 * (b0 - rinv*qiqj) *
                        (cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i] - lambda);
#    endif
                    }
#  endif // PME_SCMBAR
                  }
#  ifdef PME_SCMBAR
                  // MBAR EELEC EXCLUSIONS SC-SC
                  if (region != 0 && swtch < 0) {
                    for (int bar_i = 0; bar_i < cSim.bar_states; bar_i++) {
#    ifndef use_DPFP
                      bar_cont[bar_i] +=
                        fast_llrintf((PMEFloat)0.5 * qiqj * swtch * rinv *
                                     (cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i] -
                                      lambda) * ENERGYSCALEF);
#    else
                      bar_cont[bar_i] +=
                        (PMEFloat)0.5 * qiqj * swtch * rinv *
                        (cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i] - lambda);
#    endif
                    }
                  }
#  endif // PME_SCMBAR
                }
                if (!bSCtoSC) {
                  df += lambda * qiqj * r2inv * (swtch*rinv - d_swtch_dx);
#  ifdef PME_SCMBAR
                  // MBAR EELEC LINEAR
                  PMEFloat b0 = qiqj * swtch * rinv;
                  if (region != 0 && swtch > 0) {
                    for (int bar_i = 0; bar_i < cSim.bar_states; bar_i++) {
#    ifndef use_DPFP
                      bar_cont[bar_i] +=
                        fast_llrintf((PMEFloat)0.5 * b0 *
                                     (cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i] -
                                      lambda) * ENERGYSCALEF);
#    else
                      bar_cont[bar_i] +=
                        (PMEFloat)0.5 * b0 *
                        (cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i] - lambda);
#    endif
                    }
                  }

                  // MBAR EELEC EXCLUSIONS LINEAR
                  if (region !=0 && swtch < 0) {
                    for (int bar_i = 0; bar_i < cSim.bar_states; bar_i++) {
#    ifndef use_DPFP
                      bar_cont[bar_i] +=
                        fast_llrintf((PMEFloat)0.5 * qiqj * swtch * rinv *
                                     (cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i] -
                                      lambda) * ENERGYSCALEF);
#    else
                      bar_cont[bar_i] +=
                        (PMEFloat)0.5 * qiqj * rinv * swtch *
                        (cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i] - lambda);
#    endif
                    }
                  }
#  endif // PME_SCMBAR
                }
#  ifdef PME_ENERGY
#    ifndef use_DPFP
                PMEFloat de = ENERGYSCALEF * b0 * (PMEFloat)0.5;
                if (region == 0) {
                  if (swtch < 0) {
                    eed -= fast_llrintf(de * inv_swtch);
                  }
                  eed += fast_llrintf(de);
                }
                if (region != 0) {
                  if (bSCtoSC) {
                    eed -= fast_llrintf((PMEFloat)0.5 * ENERGYSCALEF * lambda *
                                        (qiqj*rinv - b0));
                    sc_dvdl -= fast_llrintf((PMEFloat)0.5 * ENERGYSCALEF *
                                            sTISigns[region - 1] * (qiqj * rinv - b0));
                    sc_eed[region - 1] += fast_llrintf((PMEFloat)0.5 * ENERGYSCALEF *
                                                       qiqj * rinv);
                  }
                  if (!bSCtoSC) {
                    if (swtch < 0) {
                      eed -= fast_llrintf(lambda * de * inv_swtch);
                      sc_dvdl -= fast_llrintf(sTISigns[region - 1] * de * inv_swtch);
                    }
                    eed     += fast_llrintf(lambda * de);
                    sc_dvdl += fast_llrintf(sTISigns[region - 1] * de);
                  }
                }
#    else  // use_DPFP
                PMEFloat de = b0 * (PMEFloat)0.5;
                if (region == 0) {
                  if (swtch < 0) {

                    // inv_swtch is stored from before we subtracted 1 from swtch
                    eed -= de * inv_swtch;
                  }
                  eed += de;
                }
                if (region != 0) {
                  if (bSCtoSC) {
                    eed -= (PMEFloat)0.5 * lambda * (qiqj * rinv - b0);
                    sc_dvdl -= (PMEFloat)0.5 * sTISigns[region - 1] * (qiqj * rinv - b0);
                    sc_eed[region - 1] += (PMEFloat)0.5 * qiqj * rinv;
                  }
                  if (!bSCtoSC) {
                    if (swtch < 0) {

                      // Have to take out the reciprocal which was scaled by lambda
                      eed -= lambda * de * inv_swtch;
                      sc_dvdl -= sTISigns[region - 1] * de * inv_swtch;
                    }
                    eed += lambda * de;
                    sc_dvdl += sTISigns[region - 1] * de;
                  }
                }
#    endif // use_DPFP
#  endif //PME_ENERGY
              }
              // Here ends a contingency for "not bSCtoCV"--interactions are not occuring
              // between the soft core and CV regions.

              // This next contingency takes care of interactions that DO occur between the
              // soft core and CV regions
              if (bSCtoCV) {

                // SC to C/V interaction
                PMEFloat r6 = r2 * r2 * r2;

                // Different potential for region 1 and region 2 to avoid numerical
                // instabilities--don't need to consider region == 0 because bSCtoCV is true
                PMEFloat SCtoCVDenom = (region & 0x1) ? sLambda[2] : sLambda[1];
                PMEFloat f6          = (PMEFloat)1.0 / (cSim.scalphaSP*SCtoCVDenom +
                                                        r6*term.x);
                PMEFloat f12         = f6 * f6;
                if (!(exclusion & 0x1)) {
                  df += lambda * term.y * r2 * r2 * f12 * term.x *
                        ((PMEFloat)12.0*f6 - (PMEFloat)6.0);
#  ifdef PME_SCMBAR
                  // MBAR VDW SC-C
                  for (int bar_i = 0; bar_i < cSim.bar_states; bar_i++) {
#    ifdef use_SPFP
                    PMEFloat mbarf6  =
                      (PMEFloat)1.0 /
                      (cSim.scalphaSP * ((PMEFloat)1.0 -
                                         cSim.pBarLambda[(region - 1)*cSim.bar_stride +
                                                         bar_i]) +
                       r6*term.x);
#    else
                    PMEFloat mbarf6  =
                      (PMEFloat)1.0 /
                      (cSim.scalpha * ((PMEFloat)1.0 -
                                       cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i]) +
                       r6*term.x);
#    endif
                    PMEFloat mbarf12 = mbarf6 * mbarf6;
#    ifndef use_DPFP
                    bar_cont[bar_i] -=
                      fast_llrintf((((PMEFloat)0.50 * term.y * (f12 - f6) * lambda) -
                                    (PMEFloat)0.5 * term.y * (mbarf12 - mbarf6) *
                                    cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i]) *
                                   ENERGYSCALEF);
#    else
                    bar_cont[bar_i] -=
                      (PMEFloat)0.50 * ((term.y  * (f12 - f6) * lambda) -
                                        (term.y * (mbarf12 - mbarf6)) *
                                        cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i]);
#    endif
                  }
#  endif
#  ifdef PME_ENERGY
#    ifndef use_DPFP
                  evdw += fast_llrintf(ENERGYSCALEF * (PMEFloat)0.5 * lambda * term.y *
                                       (f12 - f6));
                  sc_dvdl += fast_llrintf(ENERGYSCALEF * (PMEFloat)0.5 * sTISigns[region-1] *
                                          term.y * (f12 - f6));
#      ifdef use_SPFP
                  sc_dvdl += fast_llrintf(ENERGYSCALEF * (PMEFloat)0.5 * sTISigns[region-1] *
                                          term.y  * cSim.scalphaSP * lambda *
                                          (f12 * ((PMEFloat)2.0 * f6 - (PMEFloat)1.0)));
#      else
                  sc_dvdl += fast_llrintf(ENERGYSCALEF * (PMEFloat)0.5 * sTISigns[region-1] *
                                          term.y * cSim.scalpha * lambda *
                                          (f12 * ((PMEFloat)2.0 * f6 - (PMEFloat)1.0)));
#      endif // use_SPFP
#    else  // use_DPFP
                  evdw += (PMEFloat)0.5 * lambda * term.y * (f12 - f6);
                  sc_dvdl += (PMEFloat)0.5 * sTISigns[region-1] * term.y * (f12 - f6);
                  sc_dvdl += (PMEFloat)0.5 * sTISigns[region-1] * term.y * cSim.scalphaSP *
                             lambda * (f12 * ((PMEFloat)2.0 * f6 - (PMEFloat)1.0));
#    endif // use_DPFP
#  endif // PME_ENERGY
                }
                else {
                  swtch -= (PMEFloat)1.0;
                }
#  ifdef use_SPFP
                PMEFloat denom = (PMEFloat)1.0 / sqrt(r2 + cSim.scbetaSP * SCtoCVDenom);
#  else
                PMEFloat denom = (PMEFloat)1.0 / sqrt(r2 + cSim.scbeta * SCtoCVDenom);
#  endif
                PMEFloat denom_n = denom * denom * denom;
                if (swtch < 0) {
                  df += lambda * qiqj * r2inv * (swtch * rinv - d_swtch_dx);
#  ifdef PME_SCMBAR
                  // MBAR EXCLUSIONS SC-C
                  for (int bar_i = 0; bar_i < cSim.bar_states; bar_i++) {
#    ifndef use_DPFP
                    bar_cont[bar_i] +=
                      fast_llrintf((PMEFloat)0.5 * qiqj * swtch * rinv *
                                   (cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i] -
                                    lambda) * ENERGYSCALEF);
#    else
                    bar_cont[bar_i] +=
                      (PMEFloat)0.5 * qiqj * swtch * rinv *
                      (cSim.pBarLambda[(region - 1) * cSim.bar_stride + bar_i] - lambda);
#    endif
                  }
#  endif // PME_SCMBAR
                }
                if (swtch > 0) {
                  df += lambda * qiqj * (swtch * denom_n - d_swtch_dx * denom * rinv);
#  ifdef PME_SCMBAR
                  // MBAR EELEC SC-C
                  for (int bar_i = 0; bar_i < cSim.bar_states; bar_i++) {
#    ifdef use_SPFP
                    PMEFloat mbardenom =
                      (PMEFloat)1.0 /
                      sqrt(r2 + cSim.scbetaSP*((PMEFloat)1.0 -
                                               cSim.pBarLambda[(region - 1)*cSim.bar_stride +
                                                               bar_i]));
#    else
                    PMEFloat mbardenom =
                      (PMEFloat)1.0 /
                      sqrt(r2 + cSim.scbeta*((PMEFloat)1.0 -
                                             cSim.pBarLambda[(region - 1)*cSim.bar_stride +
                                                             bar_i]));
#    endif
#    ifndef use_DPFP
                    bar_cont[bar_i] -=
                      fast_llrintf((((PMEFloat)0.5 * qiqj * swtch * denom * lambda) -
                                    ((PMEFloat)0.5 * qiqj * swtch * mbardenom *
                                     cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i])) *
                                   ENERGYSCALEF);
#    else
                    bar_cont[bar_i] -=
                      ((PMEFloat)0.5 * qiqj * swtch * denom * lambda) -
                      ((PMEFloat)0.5 * qiqj * swtch * mbardenom *
                       cSim.pBarLambda[(region - 1) * cSim.bar_stride + bar_i]);
#    endif
                  }
#endif
                }
#  ifdef PME_ENERGY
#    ifndef use_DPFP
                if (swtch > 0) {
                  eed += fast_llrintf(ENERGYSCALEF * (PMEFloat)0.5 * qiqj * swtch * denom *
                                      lambda);
                  sc_dvdl += fast_llrintf(ENERGYSCALEF * (PMEFloat)0.5 * qiqj * swtch * denom *
                                          sTISigns[region - 1]);
                  sc_dvdl += fast_llrintf(ENERGYSCALEF * sTISigns[region - 1] * qiqj *
                                          (PMEFloat)0.25 * swtch * denom_n * cSim.scbeta *
                                          lambda);
                }
                if (swtch < 0) {

                  // Reciprocal contribution from excluded terms:
                  // swtch is really swtch - 1 since all terms were excluded
                  eed += fast_llrintf(ENERGYSCALEF * (PMEFloat)0.5 * qiqj * rinv * lambda *
                                      swtch);
                  sc_dvdl += fast_llrintf(ENERGYSCALEF * (PMEFloat)0.5 * qiqj * rinv *
                                          sTISigns[region - 1] * swtch);
                }
#    else  // use_DPFP
                if (swtch > 0) {
                  eed     += (PMEFloat)0.5 * qiqj * swtch * denom * lambda;
                  sc_dvdl += (PMEFloat)0.5 * qiqj * swtch * denom * sTISigns[region - 1];
                  sc_dvdl += sTISigns[region - 1] * qiqj * (PMEFloat)0.25 * swtch * denom_n *
                             cSim.scbeta * lambda;
                }
                if (swtch < 0) {

                  // Reciprocal contribution from excluded terms:
                  // swtch is really swtch - 1 since all terms were excluded
                  eed += (PMEFloat)0.5 * qiqj * rinv * lambda * swtch;
                  sc_dvdl += (PMEFloat)0.5 * qiqj * rinv * sTISigns[region - 1] * swtch;
                }
#    endif
#  endif // PME_ENERGY
              }
#  if !defined(use_DPFP) && defined(PME_MINIMIZATION)
              df = max(-10000.0f, min(df, 10000.0f));
#  endif
#  ifdef use_SPFP
              df *= FORCESCALEF;
#  endif
              PMEFloat dfdx = df * xij;
              PMEFloat dfdy = df * yij;
              PMEFloat dfdz = df * zij;

              // Accumulate into registers only
#  ifdef use_SPFP
              fx_i += fast_llrintf(dfdx);
              fy_i += fast_llrintf(dfdy);
              fz_i += fast_llrintf(dfdz);
#    ifdef PME_VIRIAL
              vir_11 -= fast_llrintf((PMEFloat)0.5 * xij * dfdx);
              vir_22 -= fast_llrintf((PMEFloat)0.5 * yij * dfdy);
              vir_33 -= fast_llrintf((PMEFloat)0.5 * zij * dfdz);
#    endif
#  else  // use_DPFP
              fx_i += (PMEForce)dfdx;
              fy_i += (PMEForce)dfdy;
              fz_i += (PMEForce)dfdz;
#    ifdef PME_VIRIAL
              vir_11 -= (PMEVirial)((PMEFloat)0.5 * xij * dfdx);
              vir_22 -= (PMEVirial)((PMEFloat)0.5 * yij * dfdy);
              vir_33 -= (PMEVirial)((PMEFloat)0.5 * zij * dfdz);
#    endif
#  endif
            }
          }
          exclusion >>= 1;
          shAtom.x    = __SHFL(mask1, shAtom.x, shIdx);
          shAtom.y    = __SHFL(mask1, shAtom.y, shIdx);
          shAtom.z    = __SHFL(mask1, shAtom.z, shIdx);
          shAtom.q    = __SHFL(mask1, shAtom.q, shIdx);
          shAtom.LJID = __SHFL(mask1, shAtom.LJID, shIdx);
          shAtom.TI   = __SHFL(mask1, shAtom.TI, shIdx);
          j = sNext[j];
#if (PME_ATOMS_PER_WARP == 32)
          mask1 = __BALLOT(mask1, j != tgx);
        }
#else
          mask1 = __BALLOT(mask1, j != end);
        }
#endif
        // Here ends the second of two renditions of the inner loop, this one
        // handling TI interactions.
      }
#endif // PME_SCTI
    }
    else {
      LJIDi *= cSim.LJTypes;
    }
    // Here ends a massive branch whoe purpose is to special case the
    // home cell::home cell interactions.

    // Handle remainder of line
    int tx = 0;
    while (tx < psWarp->nlEntry.NL.xatoms) {

      // Read atom ID and exclusion data
      PSATOMID(tgx) = cSim.pNLAtomList[psWarp->nlEntry.NL.offset + tgx];
      __SYNCWARP(WARP_MASK);
      if (tgx == 0)
        psWarp->nlEntry.NL.offset += GRID;
      __SYNCWARP(WARP_MASK);
      PMEMask fullExclusion = ((PMEMask*)&cSim.pNLAtomList[psWarp->nlEntry.NL.offset])
                              [tgx & cSim.NLAtomsPerWarpBitsMask];
      __SYNCWARP(WARP_MASK);
      if (tgx == 0) {
        psWarp->nlEntry.NL.offset += cSim.NLAtomsPerWarp *
                                     sizeof(PMEMask) / sizeof(unsigned int);
      }
      __SYNCWARP(WARP_MASK);
#if (PME_ATOMS_PER_WARP < 32)
      fullExclusion >>= cSim.NLAtomsPerWarp * (tgx >> cSim.NLAtomsPerWarpBits);
#endif
      unsigned int exclusion = (unsigned int)fullExclusion & cSim.NLAtomsPerWarpMask;

      // Clear shared memory forces
      psF[tgx].x = (PMEForce)0;
      psF[tgx].y = (PMEForce)0;
      psF[tgx].z = (PMEForce)0;

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
#ifdef PME_SCTI
        PSATOMTI(tgx)   = cSim.pImageTIRegion[atom];
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
#ifdef PME_SCTI
        PSATOMTI(tgx)   = 0;
#endif
      }

      // Translate all atoms into a local coordinate system within one unit
      // cell of the first atom read to avoid PBC handling within inner loops
      int cell = PSATOMID(tgx) & NLATOM_CELL_TYPE_MASK;
#if defined(PME_VIRIAL) && defined(PME_IS_ORTHOGONAL)
      PSATOMX(tgx) += sUcellf[0] * cSim.cellOffset[cell][0];
      PSATOMY(tgx) += sUcellf[4] * cSim.cellOffset[cell][1];
      PSATOMZ(tgx) += sUcellf[8] * cSim.cellOffset[cell][2];
#else
      PSATOMX(tgx) += cSim.cellOffset[cell][0];
      PSATOMY(tgx) += cSim.cellOffset[cell][1];
      PSATOMZ(tgx) += cSim.cellOffset[cell][2];
#endif

#ifndef PME_IS_ORTHOGONAL
#  ifdef PME_VIRIAL
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
      int shIdx = sNext[tgx];
#ifdef PME_SCTI
      if (__ANY(WARP_MASK, TIi) || __ANY(WARP_MASK, PSATOMTI(tgx))) {
        PMEMask mask1 = WARP_MASK;
#ifndef AMBER_PLATFORM_AMD
#pragma unroll 2
#endif
        do {
          unsigned int TIj = PSATOMTI(j);

          // If (TIi | TIj) evaluates to 111 or 110, the interaction cross two TI regions.
          // Mask those interactions out, we don't want them here.
          if ((TIi | TIj) < 6) {
            PMEFloat xij = xi - PSATOMX(j);
            PMEFloat yij = yi - PSATOMY(j);
            PMEFloat zij = zi - PSATOMZ(j);
            PMEFloat r2  = xij * xij + yij * yij + zij * zij;
            if (r2 < cSim.cut2) {
              PMEFloat rinv       = rsqrt(r2);
              PMEFloat r          = r2 * rinv;
              PMEFloat r2inv      = rinv * rinv;
              PMEFloat r6inv      = r2inv * r2inv * r2inv;
              unsigned int LJIDj  = PSATOMLJID(j);
              int bSCtoC          = ((TIi ^ TIj) & 0x1);
              int bSCtoSC         = (TIi & TIj & 1);
              int region          = ((TIi | TIj) >> 1);
              PMEFloat lambda     = sLambda[region];
              unsigned int index  = LJIDi + LJIDj + (bSCtoC * cSim.LJOffset);
#  ifndef use_DPFP
              PMEFloat2 term      = tex1Dfetch<float2>(cSim.texLJTerm, index);
#  else
              PMEFloat2 term      = cSim.pLJTerm[index];
#  endif
#  ifdef use_DPFP
              PMEFloat swtch      = erfc(cSim.ew_coeffSP * r);
#  else
              PMEFloat swtch      = fasterfc(r);
#  endif
              PMEFloat d_swtch_dx = cSim.negTwoEw_coeffRsqrtPI * exp(-cSim.ew_coeff2 * r2);
              PMEFloat qiqj       = qi * PSATOMQ(j);
              PMEFloat df         = (PMEFloat)0.0;
#  ifdef PME_ENERGY
              PMEFloat inv_swtch  = (PMEFloat)1.0 / swtch;

              // Need these for soft core - common (SC to C) forces
              PMEFloat b0         = qiqj * swtch * rinv;
#  endif

              // Check for alpha/beta scaled SC to C interaction
              if (!bSCtoC) {
                PMEFloat f6  = term.y * r6inv;
                PMEFloat f12 = term.x * r6inv * r6inv;

                // Don't scale van-der Waals for interactions between two soft core atoms
                PMEFloat lambdaVDW = (TIi & TIj & 0x1) ? (PMEFloat)1.0 : lambda;
                if (!(exclusion & 0x1)) {
                  df += lambdaVDW * (f12 - f6) * r2inv;
#  ifdef PME_SCMBAR
                  // MBAR VDW LINEAR
                  if (region != 0) {
                    for (int bar_i = 0; bar_i < cSim.bar_states; bar_i++) {
#    ifndef use_DPFP
                      bar_cont[bar_i] +=
                        fast_llrintf((PMEFloat)(f12/(PMEFloat)12.0 - f6/(PMEFloat)6.0) *
                                     (cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i] -
                                      lambdaVDW) * ENERGYSCALEF);
#    else
                      bar_cont[bar_i] +=
                        (f12/(PMEFloat)12.0 - f6/(PMEFloat)6.0) *
                        (cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i] - lambdaVDW);
#    endif
                    }
                  }
#  endif
#  ifdef PME_ENERGY
#    ifndef use_DPFP
                  if (region == 0) {
                    evdw += fast_llrintf(ENERGYSCALEF * (f12*(PMEFloat)(1.0 / 12.0) -
                                                         f6*(PMEFloat)(1.0 / 6.0)));
                  }
                  // Filter out C to C
                  if (region != 0) {
                    if (bSCtoSC) {
                      sc_evdw[region - 1] += fast_llrintf(ENERGYSCALEF *
                                                          (f12*(PMEFloat)(1.0 / 12.0) -
                                                           f6*(PMEFloat)(1.0 / 6.0)));
                    }
                    if (!bSCtoSC) {
                      evdw += fast_llrintf(ENERGYSCALEF * lambdaVDW *
                                           (f12*(PMEFloat)(1.0 / 12.0) -
                                            f6*(PMEFloat)(1.0 / 6.0)));
                      sc_dvdl += fast_llrintf(ENERGYSCALEF * sTISigns[region - 1] *
                                              (f12*(PMEFloat)(1.0 / 12.0) -
                                               f6*(PMEFloat)(1.0 / 6.0)));
                    }
                  }
#    else  // use_DPFP
                  if (region == 0) {
                    evdw += (f12 * (PMEFloat)(1.0 / 12.0) - f6 * (PMEFloat)(1.0 / 6.0));
                  }
                  if (region != 0) {
                    if (bSCtoSC == 1) {
                      sc_evdw[region - 1] += (f12*(PMEFloat)(1.0 / 12.0) -
                                              f6*(PMEFloat)(1.0 / 6.0));
                    }
                    if (!bSCtoSC) {
                      evdw += lambdaVDW * (f12*(PMEFloat)(1.0 / 12.0) -
                                           f6*(PMEFloat)(1.0 / 6.0));
                      sc_dvdl += sTISigns[region - 1] * (f12*(PMEFloat)(1.0 / 12.0) -
                                                         f6*(PMEFloat)(1.0 / 6.0));
                    }
                  }
#    endif // use_DPFP
#  endif // PME_ENERGY
                }
                else {
                  swtch -= (PMEFloat)1.0;
                }
                if (bSCtoSC) {
                  df += lambda * qiqj * r2inv * (swtch * rinv - d_swtch_dx);
                  if (swtch > 0) {
                    df += qiqj * r2inv * rinv * ((PMEFloat)1.0 - lambda);
#  ifdef PME_SCMBAR
                    // MBAR EELEC SC-SC
                    PMEFloat b0 = qiqj * swtch * rinv;
                    for(int bar_i = 0; bar_i < cSim.bar_states; bar_i++) {
#    ifndef use_DPFP
                      bar_cont[bar_i] +=
                        fast_llrintf((b0 - rinv * qiqj) *
                                     (cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i] -
                                      lambda) * ENERGYSCALEF);
#    else
                      bar_cont[bar_i] +=
                        (b0 - rinv*qiqj) *
                        (cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i] - lambda);
#    endif
                    }
#  endif
                  }
#  ifdef PME_SCMBAR
                  // MBAR EELEC EXCLUSIONS SC-SC
                  if (region !=0 && swtch < 0) {
                    for (int bar_i = 0; bar_i < cSim.bar_states; bar_i++) {
#    ifndef use_DPFP
                      bar_cont[bar_i] +=
                        fast_llrintf(qiqj * swtch * rinv *
                                     (cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i] -
                                      lambda) * ENERGYSCALEF);
#    else
                      bar_cont[bar_i] +=
                        qiqj * swtch * rinv *
                        (cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i] - lambda);
#    endif
                    }
                  }
#  endif
                }
                if (!bSCtoSC) {
                  df += lambdaVDW * qiqj * r2inv * (swtch * rinv - d_swtch_dx);
#  ifdef PME_SCMBAR
                  // MBAR EELEC LINEAR
                  PMEFloat b0 = qiqj * swtch * rinv;
                  if (region != 0 && swtch > 0) {
                    for (int bar_i = 0; bar_i < cSim.bar_states; bar_i++) {
#    ifndef use_DPFP
                      bar_cont[bar_i] +=
                        fast_llrintf(b0 * ENERGYSCALEF *
                                     (cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i] -
                                      lambda));
#    else
                      bar_cont[bar_i] +=
                        b0 * (cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i] - lambda);
#    endif
                    }
                  }

                  //MBAR EELEC EXCLUSIONS LINEAR
                  if (region !=0 && swtch < 0) {
                    for (int bar_i = 0; bar_i < cSim.bar_states; bar_i++) {
#    ifndef use_DPFP
                      bar_cont[bar_i] +=
                        fast_llrintf(qiqj * swtch * rinv * ENERGYSCALEF *
                                     (cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i] -
                                      lambda));
#    else
                      bar_cont[bar_i] +=
                        qiqj * rinv * swtch *
                        (cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i] - lambda);
#    endif
                    }
                  }
#  endif
                }
#  ifdef PME_ENERGY
#    ifndef use_DPFP
                if (region == 0) {
                  if (swtch < 0) {
                    eed -= fast_llrintf(b0 * inv_swtch * ENERGYSCALEF);
                  }
                  eed += fast_llrintf(b0 * ENERGYSCALEF);
                }
                if (region != 0) {
                  if (bSCtoSC) {
                    sc_dvdl -= fast_llrintf(sTISigns[region - 1] * ENERGYSCALEF *
                                            (qiqj * rinv - b0));
                    eed -= fast_llrintf(lambda * (qiqj * rinv - b0) * ENERGYSCALEF);
                    sc_eed[region - 1] += fast_llrintf(qiqj * rinv * ENERGYSCALEF);
                  }
                  if (!bSCtoSC) {

                    // No sc_dvdl contribution here.  We require a single topology for the
                    // linear atoms, so all C-linear terms are the same in V0 and V1
                    if (swtch < 0) {
                      eed -= fast_llrintf(b0 * inv_swtch * lambdaVDW * ENERGYSCALEF);
                      sc_dvdl -= fast_llrintf(b0 * sTISigns[region - 1] * ENERGYSCALEF *
                                              inv_swtch);
                    }
                    eed += fast_llrintf(b0 * lambdaVDW * ENERGYSCALEF);
                    sc_dvdl += fast_llrintf(b0 * sTISigns[region - 1] * ENERGYSCALEF);
                  }
                }
#    else  // use_DPFP
                if (region == 0) {

                  // Might not need this with lambdaVDW
                  if (swtch < 0) {
                    eed -= b0 * inv_swtch;
                  }
                  eed += b0;
                }
                if (region != 0) {
                  if (bSCtoSC) {
                    eed -= lambda * (qiqj * rinv - b0); //this becomes qiqj/r * swtch - 1
                    sc_dvdl -= sTISigns[region - 1] * (qiqj * rinv - b0);
                    sc_eed[region - 1] += qiqj * rinv;
                  }
                  if (!bSCtoSC) {
                    if (swtch < 0) {
                      eed     -= b0 * inv_swtch * lambdaVDW;
                      sc_dvdl -= sTISigns[region - 1] * b0 * inv_swtch;
                    }
                    eed                += b0 * lambdaVDW;
                    sc_dvdl            += sTISigns[region - 1] * b0;
                  }
                }
#    endif // use_DPFP
#  endif // PME_ENERGY
              }
              if (bSCtoC) {

                // Soft core to "common" (SC to C) interaction
                PMEFloat r6          = r2 * r2 * r2;
                PMEFloat SCtoCVDenom = (region & 0x1) ? sLambda[2] : sLambda[1];
                PMEFloat f6 = (PMEFloat)1.0 / (SCtoCVDenom*cSim.scalpha + r6*term.x);
                PMEFloat f12 = f6 * f6;
                if (!(exclusion & 0x1)) {
                  df += lambda * term.y * r2 * r2 * f12 * term.x *
                        ((PMEFloat)12.0*f6 - (PMEFloat)6.0);
#  ifdef PME_SCMBAR
                  //MBAR VDW SC-C
                  for (int bar_i = 0; bar_i < cSim.bar_states; bar_i++) {
#    ifdef use_SPFP
                    PMEFloat mbarf6  =
                      (PMEFloat)1.0 /
                      (cSim.scalphaSP * ((PMEFloat)1.0 -
                                         cSim.pBarLambda[(region - 1)*cSim.bar_stride +
                                                         bar_i]) +
                       r6*term.x);
#    else
                    PMEFloat mbarf6  =
                      (PMEFloat)1.0 /
                      (cSim.scalpha * ((PMEFloat)1.0 -
                                       cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i]) +
                       r6*term.x);
#    endif
                    PMEFloat mbarf12 = mbarf6 * mbarf6;
#    ifndef use_DPFP
                    bar_cont[bar_i] -=
                      fast_llrintf(((term.y * (f12 - f6) * lambda) -
                                    (term.y * (mbarf12 - mbarf6) *
                                     cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i])) *
                                   ENERGYSCALEF);
#    else
                    bar_cont[bar_i] -=
                      (term.y  * (f12 - f6) * lambda) -
                      (term.y * (mbarf12 - mbarf6) *
                       cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i]);
#    endif
                  }
#  endif // PME_SCMBAR
#  ifdef PME_ENERGY
#    ifndef use_DPFP
                  evdw += fast_llrintf(lambda * ENERGYSCALEF * term.y * (f12 - f6));
                  sc_dvdl += fast_llrintf(sTISigns[region - 1] * ENERGYSCALEF * term.y *
                                          (f12 - f6));
#      ifdef use_SPFP
                  sc_dvdl += fast_llrintf(sTISigns[region - 1] * ENERGYSCALEF * term.y *
                                          cSim.scalphaSP * lambda *
                                          (f12 * ((PMEFloat)2.0*f6 - (PMEFloat)1.0)));
#      else
                  sc_dvdl += fast_llrintf(sTISigns[region - 1] * ENERGYSCALEF * term.y *
                                          cSim.scalpha * lambda *
                                          (f12 * ((PMEFloat)2.0*f6 - (PMEFloat)1.0)));
#      endif //SPFP
#    else  // use_DPFP
                  evdw    += lambda * term.y * (f12 - f6);
                  sc_dvdl += sTISigns[region - 1] * term.y * (f12 - f6);
                  sc_dvdl += term.y * sTISigns[region-1] * cSim.scalpha * lambda *
                             (f12 * ((PMEFloat)2.0*f6 - (PMEFloat)1.0));
#    endif // use_DPFP
#  endif // PME_ENERGY
                }
                else {
                  swtch -= (PMEFloat)1.0;
                }
#  ifdef use_SPFP
                PMEFloat denom = (PMEFloat)1.0 / sqrt(r2 + cSim.scbetaSP * SCtoCVDenom);
#  else
                PMEFloat denom = (PMEFloat)1.0 / sqrt(r2 + cSim.scbeta * SCtoCVDenom);
#  endif
                PMEFloat denom_n = denom * denom * denom;
                if (swtch < 0) {
                  df += lambda * qiqj * r2inv * (swtch * rinv - d_swtch_dx);
#  ifdef PME_SCMBAR
                  //MBAR EXCLUSIONS SC-C
                  for (int bar_i = 0; bar_i < cSim.bar_states; bar_i++) {
#    ifndef use_DPFP
                    bar_cont[bar_i] +=
                      fast_llrintf(qiqj * swtch * rinv * ENERGYSCALEF *
                                   (cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i] -
                                    lambda));
#    else
                    bar_cont[bar_i] +=
                      qiqj * swtch * rinv *
                      (cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i] - lambda);
#    endif
                  }
#  endif
                }
                if (swtch > 0) {
                  df += lambda * qiqj * (swtch * denom_n - d_swtch_dx * denom * rinv);
#  ifdef PME_SCMBAR
                  //MBAR EELEC SC-C
                  for (int bar_i = 0; bar_i < cSim.bar_states; bar_i++) {
#    ifdef use_SPFP
                    PMEFloat mbardenom =
                      (PMEFloat)1.0 /
                      sqrt(r2 + cSim.scbetaSP*((PMEFloat)1.0 -
                                               cSim.pBarLambda[(region - 1)*cSim.bar_stride +
                                                               bar_i]));
#    else
                    PMEFloat mbardenom =
                      (PMEFloat)1.0 /
                      sqrt(r2 + cSim.scbeta*((PMEFloat)1.0 -
                                             cSim.pBarLambda[(region - 1)*cSim.bar_stride +
                                                             bar_i]));
#    endif
#    ifndef use_DPFP
                    bar_cont[bar_i] -=
                      fast_llrintf(((qiqj * swtch * denom * lambda) -
                                    (qiqj * swtch * mbardenom *
                                     cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i])) *
                                   ENERGYSCALEF);
#    else
                    bar_cont[bar_i] -=
                      (qiqj * swtch * denom * lambda) -
                      (qiqj * swtch * mbardenom *
                       cSim.pBarLambda[(region - 1)*cSim.bar_stride + bar_i]);
#    endif
                  }
#  endif
                }
#  ifdef PME_ENERGY
#    ifndef use_DPFP
                if (swtch > 0) {
                  eed += fast_llrintf(qiqj * lambda * swtch * denom * ENERGYSCALEF);
                  sc_dvdl += fast_llrintf(qiqj * swtch * denom * sTISigns[region - 1] *
                                          ENERGYSCALEF);
                  sc_dvdl += fast_llrintf(swtch * qiqj * (PMEFloat)0.5 * denom_n *
                                          cSim.scbeta * lambda * sTISigns[region - 1] *
                                          ENERGYSCALEF);
                }
                if (swtch < 0) {

                  // Reciprocal contribution from excluded terms:
                  // swtch is really swtch - 1 since all terms were excluded
                  eed += fast_llrintf(qiqj * rinv * lambda * swtch * ENERGYSCALEF);
                  sc_dvdl += fast_llrintf(qiqj * rinv * sTISigns[region - 1] * swtch *
                                          ENERGYSCALEF);
                }
#    else  // use_DPFP
                if (swtch > 0) {
                  eed += qiqj * swtch * denom * lambda;
                  sc_dvdl += qiqj * swtch * denom * sTISigns[region - 1];
                  sc_dvdl += swtch * qiqj * (PMEFloat)0.5 * denom_n * cSim.scbeta * lambda *
                             sTISigns[region - 1];
                }
                if (swtch < 0) {

                  // Reciprocal contribution from excluded terms
                  // swtch is really swtch - 1 since all terms were excluded
                  eed += qiqj * rinv * lambda * swtch;
                  sc_dvdl += qiqj * rinv * sTISigns[region - 1] * swtch;
                }
#    endif // use_DPFP
#  endif // PME_ENERGY
              }
              // Here ends the branch over whether interactions occur between the soft core
              // and common regions ("to be SC to C, or not to be SC to C").

#  if !defined(use_DPFP) && defined(PME_MINIMIZATION)
              df = max(-10000.0f, min(df, 10000.0f));
#  endif
#  ifdef use_SPFP
              df *= FORCESCALEF;
#  endif
              PMEFloat dfdx = df * xij;
              PMEFloat dfdy = df * yij;
              PMEFloat dfdz = df * zij;
#  ifdef use_SPFP
              long long int dfdx1 = fast_llrintf(dfdx);
              long long int dfdy1 = fast_llrintf(dfdy);
              long long int dfdz1 = fast_llrintf(dfdz);
              fx_i += dfdx1;
              fy_i += dfdy1;
              fz_i += dfdz1;
              psF[j].x -= dfdx1;
              psF[j].y -= dfdy1;
              psF[j].z -= dfdz1;
#    ifdef PME_VIRIAL
              vir_11 -= fast_llrintf(xij * dfdx);
              vir_22 -= fast_llrintf(yij * dfdy);
              vir_33 -= fast_llrintf(zij * dfdz);
#    endif
#  else  // use_DPFP
              PMEForce dfdx1 = (PMEForce)dfdx;
              PMEForce dfdy1 = (PMEForce)dfdy;
              PMEForce dfdz1 = (PMEForce)dfdz;
              fx_i += dfdx1;
              fy_i += dfdy1;
              fz_i += dfdz1;
              psF[j].x -= dfdx1;
              psF[j].y -= dfdy1;
              psF[j].z -= dfdz1;
#    ifdef PME_VIRIAL
              vir_11 -= (PMEForce)(xij * dfdx);
              vir_22 -= (PMEForce)(yij * dfdy);
              vir_33 -= (PMEForce)(zij * dfdz);
#    endif
#  endif // End pre-processor branch over different precision modes
            }
            // Here ends the contingency for interactions being within the cutoff
          }
          // Here ends the contingency for masking off TI interactions that cross regions
          // (111 or 110, see the possible three-bit outcomes way up above)

          exclusion >>= 1;
          shAtom.x    = __SHFL(mask1, shAtom.x, shIdx);
          shAtom.y    = __SHFL(mask1, shAtom.y, shIdx);
          shAtom.z    = __SHFL(mask1, shAtom.z, shIdx);
          shAtom.q    = __SHFL(mask1, shAtom.q, shIdx);
          shAtom.LJID = __SHFL(mask1, shAtom.LJID, shIdx);
          shAtom.TI   = __SHFL(mask1, shAtom.TI, shIdx);
          j = sNext[j];
          mask1 = __BALLOT(mask1, j != tgx);
        } while (j != tgx);
      }
      else
#endif // PME_SCTI
        if (__ANY(WARP_MASK, exclusion != 0)) {
          PMEMask mask1 = WARP_MASK;
#ifndef AMBER_PLATFORM_AMD
#pragma unroll 2
#endif
          do {
            PMEFloat xij = xi - PSATOMX(j);
            PMEFloat yij = yi - PSATOMY(j);
            PMEFloat zij = zi - PSATOMZ(j);
            PMEFloat r2  = xij*xij + yij*yij + zij*zij;
            if (r2 < cSim.cut2) {
              PMEFloat rinv       = rsqrt(r2);
              PMEFloat r          = r2 * rinv;
              PMEFloat r2inv      = rinv * rinv;
              PMEFloat r6inv      = r2inv * r2inv * r2inv;
              unsigned int LJIDj  = PSATOMLJID(j);
              unsigned int index  = LJIDi + LJIDj;
#ifndef use_DPFP
              PMEFloat2 term      = tex1Dfetch<float2>(cSim.texLJTerm, index);
#else
              PMEFloat2 term      = cSim.pLJTerm[index];
#endif
              PMEFloat df         = (PMEFloat)0.0;
              PMEFloat qiqj       = qi * PSATOMQ(j);
#ifdef use_DPFP
              PMEFloat swtch      = erfc(cSim.ew_coeffSP * r);
#else
              PMEFloat swtch      = fasterfc(r);
#endif
              PMEFloat d_swtch_dx = cSim.negTwoEw_coeffRsqrtPI * exp(-cSim.ew_coeff2 * r2);
              if (!(exclusion & 0x1)) {
#ifdef PME_FSWITCH
                PMEFloat r3inv  = rinv*r2inv;
                PMEFloat r12inv = r6inv*r6inv;
                PMEFloat df12f  = (PMEFloat)(-1) * term.x * cSim.cut6invcut6minfswitch6 *
                                  r2inv * r6inv * (r6inv - cSim.cut6inv);
                PMEFloat df6f   = (PMEFloat)(-1) * term.y * cSim.cut3invcut3minfswitch3 *
                                  r3inv * r2inv * (r3inv - cSim.cut3inv);
                PMEFloat df12   = (PMEFloat)(-1) * term.x * r2inv * r12inv;
                PMEFloat df6    = (PMEFloat)(-1) * term.y * r2inv * r6inv;
#  ifdef PME_ENERGY
                PMEFloat f12f   = term.x * cSim.cut6invcut6minfswitch6 *
                                  (r6inv - cSim.cut6inv) * (r6inv - cSim.cut6inv);
                PMEFloat f6f    = term.y * cSim.cut3invcut3minfswitch3 *
                                  (r3inv - cSim.cut3inv)*(r3inv - cSim.cut3inv);
                PMEFloat f12    = (term.x * r12inv) - (term.x * cSim.invfswitch6cut6);
                PMEFloat f6     = (term.y * r6inv) - (term.y * cSim.invfswitch3cut3);
#  endif
                if (r2 > cSim.fswitch2) {
                  df12 = df12f;
                  df6  = df6f;
#  ifdef PME_ENERGY
                  f12 = f12f;
                  f6  = f6f;
#  endif
                }
                df += df6 - df12;
#  ifdef PME_ENERGY
#    ifndef use_DPFP
                evdw += fast_llrintf(ENERGYSCALEF * (f12*(PMEFloat)(1.0 / 12.0) -
                                                     f6*(PMEFloat)(1.0 / 6.0)));
#    else
                evdw += f12 * (PMEFloat)(1.0 / 12.0) - f6 * (PMEFloat)(1.0 / 6.0);
#    endif
#  endif // PME_ENERGY
#else //PME_FSWITCH
                PMEFloat f6  = term.y * r6inv;
                PMEFloat f12 = term.x * r6inv * r6inv;
                df += (f12 - f6) * r2inv;
#ifdef PME_ENERGY
#  ifndef use_DPFP
                evdw += fast_llrintf(ENERGYSCALEF * (f12*(PMEFloat)(1.0 / 12.0) -
                                                     f6*(PMEFloat)(1.0 / 6.0)));
#  else
                evdw += f12 * (PMEFloat)(1.0 / 12.0) - f6 * (PMEFloat)(1.0 / 6.0);
#  endif
#endif
#endif //PME_FSWITCH
              }
              else {
                swtch                      -= (PMEFloat)1.0;
              }
#ifdef PME_ENERGY
              PMEFloat b0 = qiqj * swtch * rinv;
              PMEFloat b1 = b0 - qiqj * d_swtch_dx;
              df += b1 * r2inv;
#  ifndef use_DPFP
              eed += fast_llrintf(ENERGYSCALEF * b0);
#  else
              eed += b0;
#  endif
#else
              df += qiqj * (swtch * rinv - d_swtch_dx) * r2inv;
#endif
#if !defined(use_DPFP) && defined(PME_MINIMIZATION)
              df = max(-10000.0f, min(df, 10000.0f));
#endif
#ifdef use_SPFP
              df *= FORCESCALEF;
#endif
              PMEFloat dfdx = df * xij;
              PMEFloat dfdy = df * yij;
              PMEFloat dfdz = df * zij;
#ifdef use_SPFP
              long long int dfdx1 = fast_llrintf(dfdx);
              long long int dfdy1 = fast_llrintf(dfdy);
              long long int dfdz1 = fast_llrintf(dfdz);
              fx_i += dfdx1;
              fy_i += dfdy1;
              fz_i += dfdz1;
              psF[j].x -= dfdx1;
              psF[j].y -= dfdy1;
              psF[j].z -= dfdz1;
#  ifdef PME_VIRIAL
              vir_11 -= fast_llrintf(xij * dfdx);
              vir_22 -= fast_llrintf(yij * dfdy);
              vir_33 -= fast_llrintf(zij * dfdz);
#  endif
#else  // use_DPFP
              PMEForce dfdx1 = (PMEForce)dfdx;
              PMEForce dfdy1 = (PMEForce)dfdy;
              PMEForce dfdz1 = (PMEForce)dfdz;
              fx_i += dfdx1;
              fy_i += dfdy1;
              fz_i += dfdz1;
              psF[j].x -= dfdx1;
              psF[j].y -= dfdy1;
              psF[j].z -= dfdz1;
#  ifdef PME_VIRIAL
              vir_11 -= (PMEForce)(xij * dfdx);
              vir_22 -= (PMEForce)(yij * dfdy);
              vir_33 -= (PMEForce)(zij * dfdz);
#  endif
#endif // End pre-processor branch over different precision modes
            }
            exclusion >>= 1;
            shAtom.x    = __SHFL(mask1, shAtom.x, shIdx);
            shAtom.y    = __SHFL(mask1, shAtom.y, shIdx);
            shAtom.z    = __SHFL(mask1, shAtom.z, shIdx);
            shAtom.q    = __SHFL(mask1, shAtom.q, shIdx);
            shAtom.LJID = __SHFL(mask1, shAtom.LJID, shIdx);
            j = sNext[j];
            mask1 = __BALLOT(mask1, j != tgx);
          } while (j != tgx);
        }
        else {
          PMEMask mask1 = WARP_MASK;
#ifndef AMBER_PLATFORM_AMD
#pragma unroll 2
#endif
          do {
            PMEFloat xij = xi - PSATOMX(j);
            PMEFloat yij = yi - PSATOMY(j);
            PMEFloat zij = zi - PSATOMZ(j);
            PMEFloat r2  = xij*xij + yij*yij + zij*zij;
            if (r2 < cSim.cut2) {
              PMEFloat rinv      = rsqrt(r2);
              PMEFloat r         = r2 * rinv;
              PMEFloat r2inv     = rinv * rinv;
              PMEFloat r6inv     = r2inv * r2inv * r2inv;
              unsigned int LJIDj = PSATOMLJID(j);
              unsigned int index = LJIDi + LJIDj;
#ifndef use_DPFP
              PMEFloat2 term      = tex1Dfetch<float2>(cSim.texLJTerm, index);
#else
              PMEFloat2 term      = cSim.pLJTerm[index];
#endif
              PMEFloat qiqj       = qi * PSATOMQ(j);
#ifdef use_DPFP
              PMEFloat swtch      = erfc(cSim.ew_coeffSP * r) * rinv;
#else
              PMEFloat swtch      = fasterfc(r) * rinv;
#endif
              PMEFloat d_swtch_dx = cSim.negTwoEw_coeffRsqrtPI * exp(-cSim.ew_coeff2 * r2);

              PMEFloat df         = 0.0;
#ifdef PME_FSWITCH
              PMEFloat r3inv  = rinv*r2inv;
              PMEFloat r12inv = r6inv*r6inv;
              PMEFloat df12f  = (PMEFloat)(-1) * term.x * cSim.cut6invcut6minfswitch6 *
                                r2inv * r6inv * (r6inv - cSim.cut6inv);
              PMEFloat df6f   = (PMEFloat)(-1) * term.y * cSim.cut3invcut3minfswitch3 *
                                r3inv * r2inv * (r3inv - cSim.cut3inv);
              PMEFloat df12   = (PMEFloat)(-1) * term.x * r2inv * r12inv;
              PMEFloat df6    = (PMEFloat)(-1) * term.y * r2inv * r6inv;
#  ifdef PME_ENERGY
              PMEFloat f12f = term.x * cSim.cut6invcut6minfswitch6 *
                              (r6inv - cSim.cut6inv) * (r6inv - cSim.cut6inv);
              PMEFloat f6f  = term.y * cSim.cut3invcut3minfswitch3 *
                              (r3inv - cSim.cut3inv) * (r3inv - cSim.cut3inv);
              PMEFloat f12  = (term.x * r12inv) - (term.x * cSim.invfswitch6cut6);
              PMEFloat f6   = (term.y * r6inv) - (term.y * cSim.invfswitch3cut3);
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
#    ifndef use_DPFP
              evdw += fast_llrintf(ENERGYSCALEF * (f12*(PMEFloat)(1.0 / 12.0) -
                                                   f6*(PMEFloat)(1.0 / 6.0)));
#    else
              evdw += f12*(PMEFloat)(1.0 / 12.0) - f6*(PMEFloat)(1.0 / 6.0);
#    endif
              PMEFloat b0  = qiqj * swtch;
              PMEFloat b1  = b0 - qiqj * d_swtch_dx;
              df += b1 * r2inv;
#    ifndef use_DPFP
              eed += fast_llrintf(ENERGYSCALEF * b0);
#    else
              eed += b0;
#    endif
#  else
              df += qiqj * (swtch - d_swtch_dx) * r2inv;
#  endif // PME_ENERGY

#else //PME_FSWITCH
              PMEFloat f6  = term.y * r6inv;
              PMEFloat f12 = term.x * r6inv * r6inv;
              df += (f12-f6) * r2inv;
#  ifdef PME_ENERGY
#    ifndef use_DPFP
              evdw += fast_llrintf(ENERGYSCALEF * (f12 * (PMEFloat)(1.0 / 12.0) -
                                                    f6 * (PMEFloat)(1.0 /  6.0)));
#    else
              evdw += f12 * (PMEFloat)(1.0 / 12.0) - f6 * (PMEFloat)(1.0 / 6.0);
#    endif
              PMEFloat b0 = qiqj * swtch;
              PMEFloat b1 = b0 - qiqj * d_swtch_dx;
              df += b1 * r2inv;
#    ifndef use_DPFP
              eed += fast_llrintf(ENERGYSCALEF * b0);
#    else
              eed += b0;
#    endif
#  else
              df += qiqj * (swtch - d_swtch_dx) * r2inv;
#  endif // PME_ENERGY
#endif //PME_FSWITCH

#if !defined(use_DPFP) && defined(PME_MINIMIZATION)
              df = max(-10000.0f, min(df, 10000.0f));
#endif
#ifdef use_SPFP
              df *= FORCESCALEF;
#endif

              PMEFloat dfdx = df * xij;
              PMEFloat dfdy = df * yij;
              PMEFloat dfdz = df * zij;
#ifdef use_SPFP
              long long int dfdx1 = fast_llrintf(dfdx);
              long long int dfdy1 = fast_llrintf(dfdy);
              long long int dfdz1 = fast_llrintf(dfdz);
              fx_i += dfdx1;
              fy_i += dfdy1;
              fz_i += dfdz1;
              psF[j].x -= dfdx1;
              psF[j].y -= dfdy1;
              psF[j].z -= dfdz1;
#  ifdef PME_VIRIAL
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
              psF[j].x -= dfdx1;
              psF[j].y -= dfdy1;
              psF[j].z -= dfdz1;
#  ifdef PME_VIRIAL
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
            j = sNext[j];
            mask1 = __BALLOT(mask1, j != tgx);
          } while (j != tgx);
        }

      // Dump shared memory forces
      if (tx + tgx < psWarp->nlEntry.NL.xatoms) {
        int offset = (PSATOMID(tgx) >> NLATOM_CELL_SHIFT);
#ifdef use_SPFP
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset],
                  llitoulli(psF[tgx].x));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset],
                  llitoulli(psF[tgx].y));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset],
                  llitoulli(psF[tgx].z));
#else  // use_DPFP
#  ifdef PME_MINIMIZATION
        psF[tgx].x = max((PMEFloat)-10000.0, min(psF[tgx].x, (PMEFloat)10000.0));
        psF[tgx].y = max((PMEFloat)-10000.0, min(psF[tgx].y, (PMEFloat)10000.0));
        psF[tgx].z = max((PMEFloat)-10000.0, min(psF[tgx].z, (PMEFloat)10000.0));
#  endif
        atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset],
                  llitoulli(llrint(psF[tgx].x * FORCESCALE)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset],
                  llitoulli(llrint(psF[tgx].y * FORCESCALE)));
        atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset],
                  llitoulli(llrint(psF[tgx].z * FORCESCALE)));
#endif // End pre-processor branch over different precision modes
      }

      // Advance to next x tile
      tx += GRID;
    }

    // Reduce register forces if necessary
#if (PME_ATOMS_PER_WARP <= 16)
#  ifdef AMBER_PLATFORM_AMD_WARP64
    fx_i += __SHFL(WARP_MASK, fx_i, tgx + 32);
    fy_i += __SHFL(WARP_MASK, fy_i, tgx + 32);
    fz_i += __SHFL(WARP_MASK, fz_i, tgx + 32);
#  endif
    fx_i += __SHFL(WARP_MASK, fx_i, tgx + 16);
    fy_i += __SHFL(WARP_MASK, fy_i, tgx + 16);
    fz_i += __SHFL(WARP_MASK, fz_i, tgx + 16);
#  if (PME_ATOMS_PER_WARP == 8)
    fx_i += __SHFL(WARP_MASK, fx_i, tgx + 8);
    fy_i += __SHFL(WARP_MASK, fy_i, tgx + 8);
    fz_i += __SHFL(WARP_MASK, fz_i, tgx + 8);
#  endif
#endif

    // Dump register forces
    if (psWarp->nlEntry.NL.ypos + tgx < psWarp->nlEntry.NL.ymax) {
      int offset = psWarp->nlEntry.NL.ypos + tgx;
#ifdef use_SPFP
      atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset], llitoulli(fx_i));
      atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset], llitoulli(fy_i));
      atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset], llitoulli(fz_i));
#else  // use_DPFP
#  ifdef PME_MINIMIZATION
      fx_i = max((PMEFloat)-10000.0, min(fx_i, (PMEFloat)10000.0));
      fy_i = max((PMEFloat)-10000.0, min(fy_i, (PMEFloat)10000.0));
      fz_i = max((PMEFloat)-10000.0, min(fz_i, (PMEFloat)10000.0));
#  endif
      atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset],
                llitoulli(llrint(fx_i * FORCESCALE)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset],
                llitoulli(llrint(fy_i * FORCESCALE)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset],
                llitoulli(llrint(fz_i * FORCESCALE)));
#endif
    }


#ifdef PME_VIRIAL
    // Reduce virial per warp and convert to fixed point if necessary
    volatile NLVirial* psV = &sWarpVirial[threadIdx.x / GRID];
    psF[tgx].x  = vir_11;
    psF[tgx].y  = vir_22;
    psF[tgx].z  = vir_33;
    psF[tgx].x += psF[tgx ^ 1].x;
    psF[tgx].y += psF[tgx ^ 1].y;
    psF[tgx].z += psF[tgx ^ 1].z;
    psF[tgx].x += psF[tgx ^ 2].x;
    psF[tgx].y += psF[tgx ^ 2].y;
    psF[tgx].z += psF[tgx ^ 2].z;
    psF[tgx].x += psF[tgx ^ 4].x;
    psF[tgx].y += psF[tgx ^ 4].y;
    psF[tgx].z += psF[tgx ^ 4].z;
    psF[tgx].x += psF[tgx ^ 8].x;
    psF[tgx].y += psF[tgx ^ 8].y;
    psF[tgx].z += psF[tgx ^ 8].z;
    psF[tgx].x += psF[tgx ^ 16].x;
    psF[tgx].y += psF[tgx ^ 16].y;
    psF[tgx].z += psF[tgx ^ 16].z;
#ifdef AMBER_PLATFORM_AMD_WARP64
    psF[tgx].x += psF[tgx ^ 32].x;
    psF[tgx].y += psF[tgx ^ 32].y;
    psF[tgx].z += psF[tgx ^ 32].z;
#endif
    if (tgx == 0) {
#  ifndef use_DPFP
      psV->vir_11 += psF[tgx].x;
      psV->vir_22 += psF[tgx].y;
      psV->vir_33 += psF[tgx].z;
#  else  // use_DPFP
      psV->vir_11 += llrint(psF[tgx].x * FORCESCALE);
      psV->vir_22 += llrint(psF[tgx].y * FORCESCALE);
      psV->vir_33 += llrint(psF[tgx].z * FORCESCALE);
#  endif
    }
#endif // PME_VIRIAL

    // Get next Neighbor List entry
    if (tgx == 0) {
      psWarp->pos = atomicAdd(&cSim.pFrcBlkCounters[0], 1);
    }
    __SYNCWARP(WARP_MASK);
  }

#ifdef PME_VIRIAL
  volatile NLVirial* psV = &sWarpVirial[threadIdx.x / GRID];
  if ((threadIdx.x & GRID_BITS_MASK) == 0) {
    unsigned long long int val1 = llitoulli(psV->vir_11);
    unsigned long long int val2 = llitoulli(psV->vir_22);
    unsigned long long int val3 = llitoulli(psV->vir_33);
    atomicAdd(cSim.pVirial_11, val1);
    atomicAdd(cSim.pVirial_22, val2);
    atomicAdd(cSim.pVirial_33, val3);
  }
#endif

#ifdef PME_ENERGY
#  define sEED(i)  sForce[i].x
#  define sEVDW(i) sForce[i].y
#  define sDVDL(i) sForce[i].z
#  define sVDWR1(i) sForce[i].vdwR1
#  define sVDWR2(i) sForce[i].vdwR2
#  define sEEDR1(i) sForce[i].eedR1
#  define sEEDR2(i) sForce[i].eedR2
  sEED(threadIdx.x)   = eed;
  sEVDW(threadIdx.x)  = evdw;
  sEED(threadIdx.x)  += sEED(threadIdx.x ^ 1);
  sEVDW(threadIdx.x) += sEVDW(threadIdx.x ^ 1);
  sEED(threadIdx.x)  += sEED(threadIdx.x ^ 2);
  sEVDW(threadIdx.x) += sEVDW(threadIdx.x ^ 2);
  sEED(threadIdx.x)  += sEED(threadIdx.x ^ 4);
  sEVDW(threadIdx.x) += sEVDW(threadIdx.x ^ 4);
  sEED(threadIdx.x)  += sEED(threadIdx.x ^ 8);
  sEVDW(threadIdx.x) += sEVDW(threadIdx.x ^ 8);
  sEED(threadIdx.x)  += sEED(threadIdx.x ^ 16);
  sEVDW(threadIdx.x) += sEVDW(threadIdx.x ^ 16);
#ifdef AMBER_PLATFORM_AMD_WARP64
  sEED(threadIdx.x)  += sEED(threadIdx.x ^ 32);
  sEVDW(threadIdx.x) += sEVDW(threadIdx.x ^ 32);
#endif
#  ifdef PME_SCTI
  sDVDL(threadIdx.x)  = sc_dvdl;
  sDVDL(threadIdx.x) += sDVDL(threadIdx.x ^ 1);
  sDVDL(threadIdx.x) += sDVDL(threadIdx.x ^ 2);
  sDVDL(threadIdx.x) += sDVDL(threadIdx.x ^ 4);
  sDVDL(threadIdx.x) += sDVDL(threadIdx.x ^ 8);
  sDVDL(threadIdx.x) += sDVDL(threadIdx.x ^ 16);
#ifdef AMBER_PLATFORM_AMD_WARP64
  sDVDL(threadIdx.x) += sDVDL(threadIdx.x ^ 32);
#endif
#    ifdef AFE_VERBOSE
  sVDWR1(threadIdx.x) = sc_evdw[0];
  sVDWR2(threadIdx.x) = sc_evdw[1];
  sEEDR1(threadIdx.x) = sc_eed[0];
  sEEDR2(threadIdx.x) = sc_eed[1];
#    endif
#  endif
  if ((threadIdx.x & GRID_BITS_MASK) == 0) {
#  ifndef use_DPFP
    atomicAdd(cSim.pEED, llitoulli(sEED(threadIdx.x)));
    atomicAdd(cSim.pEVDW, llitoulli(sEVDW(threadIdx.x)));
#    ifdef PME_SCTI
    atomicAdd(cSim.pDVDL, llitoulli(sDVDL(threadIdx.x)));
#    endif
#  else  // use_DPFP
    atomicAdd(cSim.pEED, llitoulli(llrint(ENERGYSCALE * sEED(threadIdx.x))));
    atomicAdd(cSim.pEVDW, llitoulli(llrint(ENERGYSCALE * sEVDW(threadIdx.x))));
#    ifdef PME_SCTI
    atomicAdd(cSim.pDVDL, llitoulli(llrint(ENERGYSCALE * sDVDL(threadIdx.x))));
#    endif
#  endif
  }
#  undef sEED
#  undef sEVDW
#ifdef PME_SCMBAR
#define sMBAR(i)  sForce[i].x
  for (int i=0; i < cSim.bar_states; i++) {
    sMBAR(threadIdx.x)  = bar_cont[i];
    sMBAR(threadIdx.x) += sMBAR(threadIdx.x ^ 1);
    sMBAR(threadIdx.x) += sMBAR(threadIdx.x ^ 2);
    sMBAR(threadIdx.x) += sMBAR(threadIdx.x ^ 4);
    sMBAR(threadIdx.x) += sMBAR(threadIdx.x ^ 8);
    sMBAR(threadIdx.x) += sMBAR(threadIdx.x ^ 16);
#ifdef AMBER_PLATFORM_AMD_WARP64
    sMBAR(threadIdx.x) += sMBAR(threadIdx.x ^ 32);
#endif
    if ((threadIdx.x & GRID_BITS_MASK) == 0) {
#  ifndef use_DPFP
       atomicAdd(&cSim.pBarTot[i], llitoulli(sMBAR(threadIdx.x)));
#  else
       atomicAdd(&cSim.pBarTot[i], llitoulli(llrint(ENERGYSCALE * sMBAR(threadIdx.x))));
#  endif
    }
  }
#undef sMBAR
#endif

#  ifdef PME_SCTI
#    undef sDVDL
#    ifdef AFE_VERBOSE
#      undef sVDWR1
#      undef sVDWR2
#      undef sEEDR1
#      undef sEEDR2
#    endif
#  endif
#endif
#undef PSATOMX
#undef PSATOMY
#undef PSATOMZ
#undef PSATOMQ
#undef PSATOMLJID
#undef PSATOMID
}
