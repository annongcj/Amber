#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------
{
  struct Atom {
    PMEFloat x;
    PMEFloat y;
    PMEFloat z;
    PMEFloat q;
    unsigned int LJID;
    PMEFloat r;
    PMEFloat sig; //pwsasa sigma for a given atom
    PMEFloat eps; //pwsasa epsilon for a given atom
    PMEFloat rad; //pwsasa radius for a given atom
  };

  struct NLForce {
    PMEForce x;
    PMEForce y;
    PMEForce z;
  };

  struct NLFloat {
    PMEFloat x;
    PMEFloat y;
    PMEFloat z;
  };

#define PSATOMX(i)    shAtom.x
#define PSATOMY(i)    shAtom.y
#define PSATOMZ(i)    shAtom.z
#define PSATOMQ(i)    shAtom.q
#define PSATOMLJID(i) shAtom.LJID
#define PSATOMR(i)    shAtom.r
//Defining pwsasa variables
#define PSATOMSIG(i)  shAtom.sig //pwsasa shuffle sigma
#define PSATOMEPS(i)  shAtom.eps //pwsasa shuffle epsilon
#define PSATOMRAD(i)  shAtom.rad //pssasa shuffle radius
#ifdef AMBER_PLATFORM_AMD
#ifdef use_SPFP
#define PSFX(i)       psF.x
#define PSFY(i)       psF.y
#define PSFZ(i)       psF.z
#else
#define PSFX(i)       psF[i].x
#define PSFY(i)       psF[i].y
#define PSFZ(i)       psF[i].z
#endif
#define GBVOLATILE
#else
#define PSFX(i)       psF[i].x
#define PSFY(i)       psF[i].y
#define PSFZ(i)       psF[i].z
#define GBVOLATILE    volatile
#endif

  Atom shAtom;
#if defined(use_SPFP) && defined(AMBER_PLATFORM_AMD)
  NLFloat psF;
#else
  GBVOLATILE __shared__ NLForce sForce[GBNONBONDENERGY1_THREADS_PER_BLOCK];
#endif

#ifndef AMBER_PLATFORM_AMD
  volatile __shared__ unsigned int sPos[GBNONBONDENERGY1_THREADS_PER_BLOCK / GRID];
  volatile __shared__ unsigned int sNext[GRID];

  // Read static data
  if (threadIdx.x < GRID) {
    sNext[threadIdx.x] = (threadIdx.x + 1) & (GRID - 1);
  }
  __syncthreads();
#endif

#ifdef GB_ENERGY
  PMEAccumulator egb  = (PMEAccumulator)0;
  PMEAccumulator eelt = (PMEAccumulator)0;
  PMEAccumulator evdw = (PMEAccumulator)0;
#  if defined(use_DPFP) && defined(GB_MINIMIZATION)
  PMEAccumulator eelte = (PMEAccumulator)0;
  PMEAccumulator evdwe = (PMEAccumulator)0;
#  endif
#  ifdef GB_GBSA3
  PMEAccumulator esurf = (PMEAccumulator)0;
#  endif
#endif

  unsigned int ppos = 0;
#ifdef AMBER_PLATFORM_AMD
  ppos = ((blockIdx.x*blockDim.x + threadIdx.x) >> GRID_BITS);
  if (ppos < cSim.workUnits) {
#else
  // Initialize queue position
  GBVOLATILE unsigned int* psPos = &sPos[threadIdx.x >> GRID_BITS];
  *psPos                         = (blockIdx.x*blockDim.x + threadIdx.x) >> GRID_BITS;
  while (*psPos < cSim.workUnits) {
    ppos = *psPos;
#endif
#if defined(GB_ENERGY) && defined(use_DPFP) && defined(GB_MINIMIZATION)
    PMEFloat TLevdw = (PMEFloat)0;
    PMEFloat TLeelt = (PMEFloat)0;
#  ifdef GB_GBSA3
    PMEFloat TLesurf = (PMEFloat)0;
#  endif
    PMEFloat TLegb = (PMEFloat)0;
#endif

    // Extract cell coordinates from appropriate work unit
    unsigned int x   = cSim.pWorkUnit[ppos];
    unsigned int y   = ((x >> 2) & 0x7fff) << GRID_BITS;
    x                = (x >> 17) << GRID_BITS;
    unsigned int tgx = threadIdx.x & (GRID - 1);
    unsigned int i   = x + tgx;
    PMEFloat2 xyi    = cSim.pAtomXYSP[i];
    PMEFloat zi      = cSim.pAtomZSP[i];
    PMEFloat2 qljid  = cSim.pAtomChargeSPLJID[i];
#ifdef GB_GBSA3
    PMEFloat signpi  = cSim.pgbsa_sigma[i]; //pwsasa sigma
    PMEFloat radnpi  = cSim.pgbsa_radius[i];  //pwsasa radius
    PMEFloat epsnpi  = cSim.pgbsa_epsilon[i]; //pwsasa maxsasa
#endif
#ifndef GB_IGB6
    PMEFloat ri      = cSim.pReffSP[i];
#endif
    PMEFloat qi      = qljid.x;
    PMEMask excl = WARP_MASK;
    if (ppos < cSim.excludedWorkUnits) {
      excl = cSim.pExclusion[ppos * GRID + tgx];
    }

#if defined(use_SPFP) && defined(AMBER_PLATFORM_AMD)
    PMEFloat fx_i        = (PMEFloat)0;
    PMEFloat fy_i        = (PMEFloat)0;
    PMEFloat fz_i        = (PMEFloat)0;
#  ifndef GB_IGB6
    PMEFloat sumdeijda_i = (PMEFloat)0;
    PMEFloat sumdeijda_j = (PMEFloat)0;
#  endif
#else
    PMEForce fx_i        = (PMEForce)0;
    PMEForce fy_i        = (PMEForce)0;
    PMEForce fz_i        = (PMEForce)0;
#  ifndef GB_IGB6
    PMEForce sumdeijda_i = (PMEForce)0;
    PMEForce sumdeijda_j = (PMEForce)0;
#  endif //GB_IGB6
    unsigned int tbx = threadIdx.x - tgx;
    GBVOLATILE NLForce* psF = &sForce[tbx];
#endif
    PSFX(tgx) = 0;
    PSFY(tgx) = 0;
    PSFZ(tgx) = 0;
    __SYNCWARP(WARP_MASK);
#ifdef AMBER_PLATFORM_AMD
    unsigned int shIdx = ((tgx + 1) & (GRID - 1));
#else
    unsigned int shIdx = sNext[tgx];
#endif
    if (x == y) {
      PMEFloat xi        = xyi.x;
      PMEFloat yi        = xyi.y;
      PSATOMX(tgx)       = xi;
      PSATOMY(tgx)       = yi;
#ifdef use_DPFP
      unsigned int LJIDi = __double_as_longlong(qljid.y) * cSim.LJTypes;
#else
      unsigned int LJIDi = __float_as_uint(qljid.y) * cSim.LJTypes;
#endif
      PSATOMZ(tgx)       = zi;
      PSATOMQ(tgx)       = qi;
#ifdef use_DPFP
      PSATOMLJID(tgx)    = __double_as_longlong(qljid.y);
#else
      PSATOMLJID(tgx)    = __float_as_uint(qljid.y);
#endif
#ifdef GB_GBSA3
      PSATOMSIG(tgx)     = signpi; //pwsasa
      PSATOMEPS(tgx)     = epsnpi; //pwsasa
      PSATOMRAD(tgx)     = radnpi; //pwsasa
#endif
#ifndef GB_IGB6
      PSATOMR(tgx)       = ri;
#endif

#ifdef AMBER_PLATFORM_AMD
      shAtom.x = WarpRotateLeft<GRID>(shAtom.x);
      shAtom.y = WarpRotateLeft<GRID>(shAtom.y);
      shAtom.z = WarpRotateLeft<GRID>(shAtom.z);
      shAtom.q = WarpRotateLeft<GRID>(shAtom.q);
      shAtom.r = WarpRotateLeft<GRID>(shAtom.r);
      shAtom.LJID = WarpRotateLeft<GRID>(shAtom.LJID);
#ifdef GB_GBSA3
      shAtom.sig = WarpRotateLeft<GRID>(shAtom.sig); //pwsasa
      shAtom.eps = WarpRotateLeft<GRID>(shAtom.eps); //pwsasa
      shAtom.rad = WarpRotateLeft<GRID>(shAtom.rad); //pwsasa
#endif
      unsigned int j = ((tgx + 1) & (GRID - 1));
#else
      shAtom.x    = __SHFL(WARP_MASK, shAtom.x, shIdx);
      shAtom.y    = __SHFL(WARP_MASK, shAtom.y, shIdx);
      shAtom.z    = __SHFL(WARP_MASK, shAtom.z, shIdx);
      shAtom.q    = __SHFL(WARP_MASK, shAtom.q, shIdx);
      shAtom.r    = __SHFL(WARP_MASK, shAtom.r, shIdx);
      shAtom.LJID = __SHFL(WARP_MASK, shAtom.LJID, shIdx);
#ifdef GB_GBSA3
      shAtom.sig  = __SHFL(WARP_MASK, shAtom.sig, shIdx); //pwsasa
      shAtom.eps  = __SHFL(WARP_MASK, shAtom.eps, shIdx); //pwsasa
      shAtom.rad  = __SHFL(WARP_MASK, shAtom.rad, shIdx); //pwsasa
#endif
      unsigned int j = sNext[tgx];
#endif
      PMEMask mask1 = __BALLOT(WARP_MASK, j != tgx);
      while (j != tgx) {
        PMEFloat xij  = xi - PSATOMX(j);
        PMEFloat yij  = yi - PSATOMY(j);
        PMEFloat zij  = zi - PSATOMZ(j);
        PMEFloat r2   = xij*xij + yij*yij + zij*zij;
        PMEFloat qiqj = qi * PSATOMQ(j);
        PMEFloat v5   = rsqrt(r2);
#ifndef GB_IGB6
        PMEFloat rj     = PSATOMR(j);
        // GB with explicit ions gbion==2
        PMEFloat v1;
        if (cSim.gbion==2) {
                unsigned int j_glb = y + j;
                unsigned int ionmaski = cSim.pAtomIonMask[i];
                unsigned int ionmaskj = cSim.pAtomIonMask[j_glb];
                if (ionmaski>0 && ionmaskj>0){ // ion-ion pair
                        v1 = exp(-r2 / ((PMEFloat)4.0 * cSim.gi_coef_2 * ri * rj));
                }
                else if (ionmaski==0 && ionmaskj==0) { // solute-solute pair
                        v1 = exp(-r2 / ((PMEFloat)4.0 * ri * rj)); // no change
                }
                else { // ion-solute pair
                        v1 = exp(-r2 / ((PMEFloat)4.0 * cSim.gi_coef_1 * ri * rj));
                }
        }
        else if (cSim.gbion==3) { // separate cation and anion
                unsigned int j_glb = y + j;
                unsigned int ionmaski = cSim.pAtomIonMask[i];
                unsigned int ionmaskj = cSim.pAtomIonMask[j_glb];
                unsigned int ion_mask_pn_i;
                unsigned int ion_mask_pn_j;
                unsigned int ion_mask_pn_index;
                if (ionmaski == 0){ //atom i is solute
                        ion_mask_pn_i=0;
                }
                else if (ionmaski == 3){ // atom i is anion
                        ion_mask_pn_i=1;
                }
                else { // atom i is cation
                        ion_mask_pn_i=2;
                }
                if (ionmaskj == 0){ //atom j is solute
                        ion_mask_pn_j=0;
                }
                else if (ionmaskj == 3){ // atom j is anion
                        ion_mask_pn_j=1;
                }
                else { // atom j is cation
                        ion_mask_pn_j=2;
                }
                ion_mask_pn_index = ion_mask_pn_i*3+ion_mask_pn_j;
                if (ion_mask_pn_index == 0) { // sulote-solute pair
                        v1 = exp(-r2 / ((PMEFloat)4.0 * ri * rj)); // no change
                }
                else if (ion_mask_pn_index == 1 || ion_mask_pn_index == 3){ // solute-anion
                        v1 = exp(-r2 / ((PMEFloat)4.0 * cSim.gi_coef_1_n * ri * rj));
                }
                else if (ion_mask_pn_index == 2 || ion_mask_pn_index == 6){ // solute-cation
                        v1 = exp(-r2 / ((PMEFloat)4.0 * cSim.gi_coef_1_p * ri * rj));
                }
                else if (ion_mask_pn_index == 5 || ion_mask_pn_index == 7){ // cation-anion
                        v1 = exp(-r2 / ((PMEFloat)4.0 * cSim.gi_coef_2_pn * ri * rj));
                }
                else if (ion_mask_pn_index == 4){ // anion-anion
                        v1 = exp(-r2 / ((PMEFloat)4.0 * cSim.gi_coef_2_nn * ri * rj));
                }
                else { // (ion_mask_pn_index == 8)  cation-cation
                        v1 = exp(-r2 / ((PMEFloat)4.0 * cSim.gi_coef_2_pp * ri * rj));
                }
        }
        else { // gbion!=2 or 3, no explicit ions, normal GB
                v1     = exp(-r2 / ((PMEFloat)4.0 * ri * rj));
        }
        PMEFloat v3     = r2 + rj*ri*v1;
        PMEFloat v2     = rsqrt(v3);
        PMEFloat expmkf = cSim.extdiel_inv;
        PMEFloat fgbk   = (PMEFloat)0.0;
        PMEFloat fgbi   = v2;
#ifdef GB_ENERGY
        PMEFloat mul    = fgbi;
#endif
        if (cSim.gb_kappa != (PMEFloat)0.0) {
          v3           = -cSim.gb_kappa / v2;
          PMEFloat v4  = exp(v3);
          expmkf      *= v4;
          fgbk         = v3 * expmkf;
          if (cSim.alpb == 1) {
            fgbk += fgbk * cSim.one_arad_beta * (-v3 * cSim.gb_kappa_inv);
#ifdef GB_ENERGY
            mul += cSim.one_arad_beta;
#endif
          }
        }
        PMEFloat dl = cSim.intdiel_inv - expmkf;
        
        if (cSim.gbion == 2) { // GB with explict ions
                unsigned int j_glb = y + j;
                unsigned int ionmaski = cSim.pAtomIonMask[i];
                unsigned int ionmaskj = cSim.pAtomIonMask[j_glb];
                PMEFloat intdiel_ion_1_inv = 1.0 / cSim.intdiel_ion_1;
                PMEFloat intdiel_ion_2_inv = 1.0 / cSim.intdiel_ion_2;
                if (ionmaski>0 && ionmaskj>0){ // ion-ion pair
                        dl = cSim.intdiel_inv * intdiel_ion_2_inv - expmkf;
                }
                else if (ionmaski==0 && ionmaskj==0) { // solute-solute pair
                         // no change
                }
                else { // ion-solute pair
                        dl = cSim.intdiel_inv * intdiel_ion_1_inv - expmkf;
                }
        }
        else if (cSim.gbion == 3) { // GB with explict ions, separate cation and anion
                unsigned int j_glb = y + j;
                unsigned int ionmaski = cSim.pAtomIonMask[i];
                unsigned int ionmaskj = cSim.pAtomIonMask[j_glb];
                unsigned int ion_mask_pn_i;
                unsigned int ion_mask_pn_j;
                unsigned int ion_mask_pn_index;
                PMEFloat intdiel_ion_1_p_inv = 1.0 / cSim.intdiel_ion_1_p;
                PMEFloat intdiel_ion_1_n_inv = 1.0 / cSim.intdiel_ion_1_n;
                PMEFloat intdiel_ion_2_pp_inv = 1.0 / cSim.intdiel_ion_2_pp;
                PMEFloat intdiel_ion_2_pn_inv = 1.0 / cSim.intdiel_ion_2_pn;
                PMEFloat intdiel_ion_2_nn_inv = 1.0 / cSim.intdiel_ion_2_nn;
                if (ionmaski == 0){ //atom i is solute
                        ion_mask_pn_i=0;
                }
                else if (ionmaski == 3){ // atom i is anion
                        ion_mask_pn_i=1;
                }
                else { // atom i is cation
                        ion_mask_pn_i=2;
                }
                if (ionmaskj == 0){ //atom j is solute
                        ion_mask_pn_j=0;
                }
                else if (ionmaskj == 3){ // atom j is anion
                        ion_mask_pn_j=1;
                }
                else { // atom j is cation
                        ion_mask_pn_j=2;
                }
                ion_mask_pn_index = ion_mask_pn_i*3+ion_mask_pn_j;
                if (ion_mask_pn_index == 0) { // sulote-solute pair
                        // no change
                }
                else if (ion_mask_pn_index == 1 || ion_mask_pn_index == 3){ // solute-anion
                        dl = cSim.intdiel_inv * intdiel_ion_1_n_inv - expmkf;
                }
                else if (ion_mask_pn_index == 2 || ion_mask_pn_index == 6){ // solute-cation
                        dl = cSim.intdiel_inv * intdiel_ion_1_p_inv - expmkf;
                }
                else if (ion_mask_pn_index == 5 || ion_mask_pn_index == 7){ // cation-anion
                        dl = cSim.intdiel_inv * intdiel_ion_2_pn_inv - expmkf;
                }
                else if (ion_mask_pn_index == 4){ // anion-anion
                        dl = cSim.intdiel_inv * intdiel_ion_2_nn_inv - expmkf;
                }
                else { // (ion_mask_pn_index == 8)  cation-cation
                        dl = cSim.intdiel_inv * intdiel_ion_2_pp_inv - expmkf;
                }
        }

#ifdef GB_ENERGY
        PMEFloat e = -qiqj * dl * mul;
#  ifndef use_DPFP
        egb += fast_llrintf((PMEFloat)0.5 * ENERGYSCALEF * e);
#  elif defined(GB_MINIMIZATION)
        TLegb += (PMEForce)((PMEFloat)0.5 * e);
#  else
        egb += llrint((PMEFloat)0.5 * ENERGYSCALE * e);
#  endif
#endif //GB_ENERGY
        // 1.0 / fij^3
        PMEFloat temp4 = fgbi * fgbi * fgbi;

        // Here, and in the gas-phase part, "de" contains -(1/r)(dE/dr)
        PMEFloat temp6  = -qiqj * temp4 * (dl + fgbk);
        PMEFloat temp1  = v1;
#ifdef use_SPFP
        PMEFloat de     = temp6 * (FORCESCALEF - (PMEFloat)0.25 * FORCESCALEF * temp1);
        PMEFloat temp5  = (PMEFloat)0.50 * FORCESCALEF * temp1 * temp6 *
                          (ri*rj + (PMEFloat)0.25*r2);
#ifdef AMBER_PLATFORM_AMD
        sumdeijda_i    += (ri * temp5);
#else
        sumdeijda_i    += fast_llrintf(ri * temp5);
#endif
#else // use_DPFP
        PMEFloat de     = temp6 * ((PMEFloat)1.0 - (PMEFloat)0.25 * temp1);
        PMEFloat temp5  = (PMEFloat)0.50 * temp1 * temp6 * (ri*rj + (PMEFloat)0.25*r2);
        sumdeijda_i    += (PMEDouble)(ri * temp5);
#endif
#ifdef GB_GBSA3
        // pwsasa: calculate reflect-vdwlike
        if (epsnpi > (PMEFloat)1.0e-8 && PSATOMEPS(j) > (PMEFloat)1.0e-8) {
          PMEFloat dist     = (PMEFloat)1.0 / v5; // 1.0 / rij
          if (dist < (radnpi + PSATOMRAD(j))) {
            PMEFloat tempsum  = sqrt(epsnpi * PSATOMEPS(j)) / (PMEFloat)6.0;
            PMEFloat numer   = signpi + PSATOMSIG(j);
            PMEFloat numer2  = numer * numer;
            PMEFloat numer4  = numer2 * numer2;
            PMEFloat numer10 = numer4 * numer4 * numer2;
            PMEFloat Sij     = radnpi + PSATOMRAD(j) + numer;
            PMEFloat denom   = (PMEFloat)1.0 / (Sij - dist);
            PMEFloat denom2  = denom * denom;
            PMEFloat denom4  = denom2 * denom2;
            PMEFloat denom10 = denom4 * denom4 * denom2;
            PMEFloat reflectvdwA = (PMEFloat)4.0  * tempsum * numer10 * denom10;
            PMEFloat reflectvdwB = (PMEFloat)10.0 * tempsum * numer4 * denom4;
            PMEFloat AgNPe       = (reflectvdwB - reflectvdwA - sqrt(epsnpi * PSATOMEPS(j))) *
                                   cSim.surften * (PMEFloat)0.6;
            PMEFloat tempreflect = denom;
            PMEFloat AgNPde      = ((PMEFloat)10.0*reflectvdwA - (PMEFloat)4.0*reflectvdwB) *
                                   tempreflect * v5 * (PMEFloat)1.2 * cSim.surften;
#ifdef use_SPFP
            de    += FORCESCALEF * AgNPde;
#  ifdef GB_ENERGY
            esurf += fast_llrintf(ENERGYSCALEF * AgNPe);
#  endif // end if GB_ENERGY
#else // use_DPFP
            de    += AgNPde;
#  ifdef GB_ENERGY
#    ifdef GB_MINIMIZATION
            TLesurf += (double)AgNPe;
#    else
            esurf += llrint(ENERGYSCALE * AgNPe);
#    endif
#  endif //GB_ENERGY
#endif
          }
        }
// pwsasa calculations end
#endif //GB_GBSA3

#else  //GB_IGB6
        PMEFloat de     = 0.0;
#endif //GB_IGB6
        PMEFloat rinv  = v5; // 1.0 / rij
        PMEFloat r2inv = rinv * rinv;
        PMEFloat eel   = cSim.intdiel_inv * qiqj * rinv;

        if (cSim.gbion == 2) { // GB with explict ions
                unsigned int j_glb = y + j;
                unsigned int ionmaski = cSim.pAtomIonMask[i];
                unsigned int ionmaskj = cSim.pAtomIonMask[j_glb];
                PMEFloat intdiel_ion_1_inv = 1.0 / cSim.intdiel_ion_1;
                PMEFloat intdiel_ion_2_inv = 1.0 / cSim.intdiel_ion_2;
                if (ionmaski>0 && ionmaskj>0){ // ion-ion pair
                        eel = cSim.intdiel_inv * intdiel_ion_2_inv * qiqj * rinv;
                }
                else if (ionmaski==0 && ionmaskj==0) { // solute-solute pair
                         // no change
                }
                else { // ion-solute pair
                        eel = cSim.intdiel_inv * intdiel_ion_1_inv * qiqj * rinv;
                }
        }
        else if (cSim.gbion == 3) { // GB with explict ions, separate cation and anion
                unsigned int j_glb = y + j;
                unsigned int ionmaski = cSim.pAtomIonMask[i];
                unsigned int ionmaskj = cSim.pAtomIonMask[j_glb];
                unsigned int ion_mask_pn_i;
                unsigned int ion_mask_pn_j;
                unsigned int ion_mask_pn_index;
                PMEFloat intdiel_ion_1_p_inv = 1.0 / cSim.intdiel_ion_1_p;
                PMEFloat intdiel_ion_1_n_inv = 1.0 / cSim.intdiel_ion_1_n;
                PMEFloat intdiel_ion_2_pp_inv = 1.0 / cSim.intdiel_ion_2_pp;
                PMEFloat intdiel_ion_2_pn_inv = 1.0 / cSim.intdiel_ion_2_pn;
                PMEFloat intdiel_ion_2_nn_inv = 1.0 / cSim.intdiel_ion_2_nn;
                if (ionmaski == 0){ //atom i is solute
                        ion_mask_pn_i=0;
                }
                else if (ionmaski == 3){ // atom i is anion
                        ion_mask_pn_i=1;
                }
                else { // atom i is cation
                        ion_mask_pn_i=2;
                }
                if (ionmaskj == 0){ //atom j is solute
                        ion_mask_pn_j=0;
                }
                else if (ionmaskj == 3){ // atom j is anion
                        ion_mask_pn_j=1;
                }
                else { // atom j is cation
                        ion_mask_pn_j=2;
                }
                ion_mask_pn_index = ion_mask_pn_i*3+ion_mask_pn_j;
                if (ion_mask_pn_index == 0) { // sulote-solute pair
                        // no change
                }
                else if (ion_mask_pn_index == 1 || ion_mask_pn_index == 3){ // solute-anion
                        eel = cSim.intdiel_inv * intdiel_ion_1_n_inv * qiqj * rinv;
                }
                else if (ion_mask_pn_index == 2 || ion_mask_pn_index == 6){ // solute-cation
                        eel = cSim.intdiel_inv * intdiel_ion_1_p_inv * qiqj * rinv;
                }
                else if (ion_mask_pn_index == 5 || ion_mask_pn_index == 7){ // cation-anion
                        eel = cSim.intdiel_inv * intdiel_ion_2_pn_inv * qiqj * rinv;
                }
                else if (ion_mask_pn_index == 4){ // anion-anion
                        eel = cSim.intdiel_inv * intdiel_ion_2_nn_inv * qiqj * rinv;
                }
                else { // (ion_mask_pn_index == 8)  cation-cation
                        eel = cSim.intdiel_inv * intdiel_ion_2_pp_inv * qiqj * rinv;
                }
        }

        PMEFloat r6inv = r2inv * r2inv * r2inv;
        unsigned int LJIDj = PSATOMLJID(j);
        unsigned int index = LJIDi + LJIDj;
#ifndef use_DPFP
        PMEFloat2 term = tex1Dfetch<float2>(cSim.texLJTerm, index);
#else
        PMEFloat2 term = cSim.pLJTerm[index];
#endif
        PMEFloat f6    = term.y * r6inv;
        PMEFloat f12   = term.x * r6inv * r6inv;
#ifdef use_SPFP // Beginning of pre-processor branch over precision modes
        if (excl & 0x1) {
#  ifdef GB_MINIMIZATION
          de += FORCESCALEF * max(-10000.0f, min(((f12 - f6) + eel) * r2inv, 10000.0f));
#  else
          de += FORCESCALEF * ((f12 - f6) + eel) * r2inv;
#  endif
#  ifdef GB_ENERGY
          eelt += fast_llrintf((PMEFloat)0.5 * ENERGYSCALEF * eel);
          evdw += fast_llrintf(ENERGYSCALEF * ((PMEFloat)(0.5 / 12.0)*f12 -
                                               (PMEFloat)(0.5 / 6.0)*f6));
#  endif
        }
        PMEFloat dedx = de * xij;
        PMEFloat dedy = de * yij;
        PMEFloat dedz = de * zij;
#ifdef AMBER_PLATFORM_AMD
        fx_i += dedx;
        fy_i += dedy;
        fz_i += dedz;
#else
        fx_i += fast_llrintf(dedx);
        fy_i += fast_llrintf(dedy);
        fz_i += fast_llrintf(dedz);
#endif
#else  // use_DPFP
        if (excl & 0x1) {
          de += ((f12 - f6) + eel) * r2inv;
#  ifdef GB_ENERGY
#    ifdef GB_MINIMIZATION
          TLeelt += (double)((PMEFloat)0.5 * eel);
          TLevdw += (double)((PMEFloat)(0.5 / 12.0)*f12 - (0.5 / 6.0)*f6);
#    else
          eelt += llrint((PMEFloat)0.5 * ENERGYSCALE * eel);
          evdw += llrint(ENERGYSCALE * ((PMEFloat)(0.5 / 12.0)*f12 -
                                        (PMEFloat)(0.5 / 6.0)*f6));
#    endif
#  endif
        }
        PMEFloat dedx = de * xij;
        PMEFloat dedy = de * yij;
        PMEFloat dedz = de * zij;
        fx_i += (PMEDouble)dedx;
        fy_i += (PMEDouble)dedy;
        fz_i += (PMEDouble)dedz;
#endif // End of pre-processor branch over different precision modes
        excl >>= 1;
#ifdef AMBER_PLATFORM_AMD
        shAtom.x = WarpRotateLeft<GRID>(shAtom.x);
        shAtom.y = WarpRotateLeft<GRID>(shAtom.y);
        shAtom.z = WarpRotateLeft<GRID>(shAtom.z);
        shAtom.q = WarpRotateLeft<GRID>(shAtom.q);
        shAtom.r = WarpRotateLeft<GRID>(shAtom.r);
        shAtom.LJID = WarpRotateLeft<GRID>(shAtom.LJID);
#ifdef GB_GBSA3
        shAtom.sig = WarpRotateLeft<GRID>(shAtom.sig); //pwsasa
        shAtom.eps = WarpRotateLeft<GRID>(shAtom.eps); //pwsasa
        shAtom.rad = WarpRotateLeft<GRID>(shAtom.rad); //pwsasa
#endif
        j = ((j + 1) & (GRID - 1));
#else
        shAtom.x    = __SHFL(mask1, shAtom.x, shIdx);
        shAtom.y    = __SHFL(mask1, shAtom.y, shIdx);
        shAtom.z    = __SHFL(mask1, shAtom.z, shIdx);
        shAtom.q    = __SHFL(mask1, shAtom.q, shIdx);
        shAtom.r    = __SHFL(mask1, shAtom.r, shIdx);
        shAtom.LJID = __SHFL(mask1, shAtom.LJID, shIdx);
#ifdef GB_GBSA3
        shAtom.sig  = __SHFL(mask1, shAtom.sig, shIdx); //pwsasa
        shAtom.eps  = __SHFL(mask1, shAtom.eps, shIdx); //pwsasa
        shAtom.rad  = __SHFL(mask1, shAtom.rad, shIdx); //pwsasa
#endif
        j = sNext[j];
        mask1 = __BALLOT(mask1, j != tgx);
#endif
      }
      int offset = x + tgx;
#ifdef use_SPFP
#ifdef AMBER_PLATFORM_AMD
      atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset],
                llitoulli(fast_llrintf(fx_i)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset],
                llitoulli(fast_llrintf(fy_i)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset],
                llitoulli(fast_llrintf(fz_i)));
#  ifndef GB_IGB6
      atomicAdd((unsigned long long int*)&cSim.pSumdeijdaAccumulator[offset],
                llitoulli(fast_llrintf(sumdeijda_i)));
#  endif
#else
      atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset], llitoulli(fx_i));
      atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset], llitoulli(fy_i));
      atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset], llitoulli(fz_i));
#  ifndef GB_IGB6
      atomicAdd((unsigned long long int*)&cSim.pSumdeijdaAccumulator[offset],
                llitoulli(sumdeijda_i));
#  endif
#endif
#else  // use_DPFP
#ifdef GB_MINIMIZATION
      PMEFloat i;
      fx_i = modf(fx_i, &i);
      atomicAdd((unsigned long long int*)&cSim.pIntForceXAccumulator[offset], llitoulli(llrint(i)));
      fy_i = modf(fy_i, &i);
      atomicAdd((unsigned long long int*)&cSim.pIntForceYAccumulator[offset], llitoulli(llrint(i)));
      fz_i = modf(fz_i, &i);
      atomicAdd((unsigned long long int*)&cSim.pIntForceZAccumulator[offset], llitoulli(llrint(i)));
#endif
      atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset],
                llitoulli(llrint(fx_i * FORCESCALE)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset],
                llitoulli(llrint(fy_i * FORCESCALE)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset],
                llitoulli(llrint(fz_i * FORCESCALE)));
#  ifndef GB_IGB6
      //printf("PI %d %20.10f\n", offset, sumdeijda_i);
      atomicAdd((unsigned long long int*)&cSim.pSumdeijdaAccumulator[offset],
                llitoulli(llrint(sumdeijda_i * FORCESCALE)));
#  endif
#endif
    }
    else {
      unsigned int j   = y + tgx;
      PMEFloat2 xyj    = cSim.pAtomXYSP[j];
      PMEFloat2 qljidj = cSim.pAtomChargeSPLJID[j];
      PSATOMZ(tgx)     = cSim.pAtomZSP[j];
      PSATOMR(tgx)     = cSim.pReffSP[j];
      PSATOMQ(tgx)     = qljidj.x;
#ifdef GB_GBSA3
      PSATOMSIG(tgx)   = cSim.pgbsa_sigma[j]; //pwsasa
      PSATOMRAD(tgx)   = cSim.pgbsa_radius[j]; //pwsasa
      PSATOMEPS(tgx)   = cSim.pgbsa_epsilon[j]; //pwsasa
#endif
#ifdef use_DPFP
      PSATOMLJID(tgx)  = __double_as_longlong(qljidj.y);
#else
      PSATOMLJID(tgx)  = __float_as_uint(qljidj.y);
#endif
      PMEFloat xi        = xyi.x;
      PMEFloat yi        = xyi.y;
      PMEFloat qi        = qljid.x;
#ifdef use_DPFP
      unsigned int LJIDi = __double_as_longlong(qljid.y) * cSim.LJTypes;
#else
      unsigned int LJIDi = __float_as_uint(qljid.y) * cSim.LJTypes;
#endif
      PSATOMX(tgx)       = xyj.x;
      PSATOMY(tgx)       = xyj.y;
      j = tgx;
      PMEMask mask1 = WARP_MASK;
#ifdef AMBER_PLATFORM_AMD
      for (unsigned int xx = 0; xx < GRID; ++xx) {
#else
      do {
#endif
        PMEFloat xij  = xi - PSATOMX(j);
        PMEFloat yij  = yi - PSATOMY(j);
        PMEFloat zij  = zi - PSATOMZ(j);
        PMEFloat r2   = xij * xij + yij * yij + zij * zij;
        PMEFloat qiqj = qi * PSATOMQ(j);
        PMEFloat v5   = rsqrt(r2);
#ifndef GB_IGB6
        PMEFloat rj     = PSATOMR(j);
        // GB with explicit ions gbion==2
        PMEFloat v1;
        if (cSim.gbion==2) {
                unsigned int j_glb = y + j;
                unsigned int ionmaski = cSim.pAtomIonMask[i];
                unsigned int ionmaskj = cSim.pAtomIonMask[j_glb];
                if (ionmaski>0 && ionmaskj>0){ // ion-ion pair
                        v1 = exp(-r2 / ((PMEFloat)4.0 * cSim.gi_coef_2 * ri * rj));
                }
                else if (ionmaski==0 && ionmaskj==0) { // solute-solute pair
                        v1 = exp(-r2 / ((PMEFloat)4.0 * ri * rj)); // no change
                }
                else { // ion-solute pair
                        v1 = exp(-r2 / ((PMEFloat)4.0 * cSim.gi_coef_1 * ri * rj));
                }
        }
        else if (cSim.gbion==3) { // separate cation and anion
                unsigned int j_glb = y + j;
                unsigned int ionmaski = cSim.pAtomIonMask[i];
                unsigned int ionmaskj = cSim.pAtomIonMask[j_glb];
                unsigned int ion_mask_pn_i;
                unsigned int ion_mask_pn_j;
                unsigned int ion_mask_pn_index;
                if (ionmaski == 0){ //atom i is solute
                        ion_mask_pn_i=0;
                }
                else if (ionmaski == 3){ // atom i is anion
                        ion_mask_pn_i=1;
                }
                else { // atom i is cation
                        ion_mask_pn_i=2;
                }
                if (ionmaskj == 0){ //atom j is solute
                        ion_mask_pn_j=0;
                }
                else if (ionmaskj == 3){ // atom j is anion
                        ion_mask_pn_j=1;
                }
                else { // atom j is cation
                        ion_mask_pn_j=2;
                }
                ion_mask_pn_index = ion_mask_pn_i*3+ion_mask_pn_j;
                if (ion_mask_pn_index == 0) { // sulote-solute pair
                        v1 = exp(-r2 / ((PMEFloat)4.0 * ri * rj)); // no change
                }
                else if (ion_mask_pn_index == 1 || ion_mask_pn_index == 3){ // solute-anion
                        v1 = exp(-r2 / ((PMEFloat)4.0 * cSim.gi_coef_1_n * ri * rj));
                }
                else if (ion_mask_pn_index == 2 || ion_mask_pn_index == 6){ // solute-cation
                        v1 = exp(-r2 / ((PMEFloat)4.0 * cSim.gi_coef_1_p * ri * rj));
                }
                else if (ion_mask_pn_index == 5 || ion_mask_pn_index == 7){ // cation-anion
                        v1 = exp(-r2 / ((PMEFloat)4.0 * cSim.gi_coef_2_pn * ri * rj));
                }
                else if (ion_mask_pn_index == 4){ // anion-anion
                        v1 = exp(-r2 / ((PMEFloat)4.0 * cSim.gi_coef_2_nn * ri * rj));
                }
                else { // (ion_mask_pn_index == 8)  cation-cation
                        v1 = exp(-r2 / ((PMEFloat)4.0 * cSim.gi_coef_2_pp * ri * rj));
                }
        }
        else { // gbion!=2 or 3, no explicit ions, normal GB
                v1     = exp(-r2 / ((PMEFloat)4.0 * ri * rj));
        }
        // PMEFloat v1     = exp(-r2 / ((PMEFloat)4.0 * ri * rj));
        PMEFloat v3     = r2 + rj*ri*v1;
        PMEFloat v2     = rsqrt(v3);
        PMEFloat expmkf = cSim.extdiel_inv;
        PMEFloat fgbk   = (PMEFloat)0.0;
        PMEFloat fgbi   = v2;
#  ifdef GB_ENERGY
        PMEFloat mul    = fgbi;
#  endif
        if (cSim.gb_kappa != (PMEFloat)0.0) {
          v3           = -cSim.gb_kappa / v2;
          PMEFloat v4  = exp(v3);
          expmkf      *= v4;
          fgbk         = v3 * expmkf;
          if (cSim.alpb == 1) {
            fgbk += fgbk * cSim.one_arad_beta * (-v3 * cSim.gb_kappa_inv);
#  ifdef GB_ENERGY
            mul  += cSim.one_arad_beta;
#  endif
          }
        }
        PMEFloat dl = cSim.intdiel_inv - expmkf;

        if (cSim.gbion == 2) { // GB with explict ions
                unsigned int j_glb = y + j;
                unsigned int ionmaski = cSim.pAtomIonMask[i];
                unsigned int ionmaskj = cSim.pAtomIonMask[j_glb];
                PMEFloat intdiel_ion_1_inv = 1.0 / cSim.intdiel_ion_1;
                PMEFloat intdiel_ion_2_inv = 1.0 / cSim.intdiel_ion_2;
                if (ionmaski>0 && ionmaskj>0){ // ion-ion pair
                        dl = cSim.intdiel_inv * intdiel_ion_2_inv - expmkf;
                }
                else if (ionmaski==0 && ionmaskj==0) { // solute-solute pair
                         // no change
                }
                else { // ion-solute pair
                        dl = cSim.intdiel_inv * intdiel_ion_1_inv - expmkf;
                }
        }
        else if (cSim.gbion == 3) { // GB with explict ions, separate cation and anion
                unsigned int j_glb = y + j;
                unsigned int ionmaski = cSim.pAtomIonMask[i];
                unsigned int ionmaskj = cSim.pAtomIonMask[j_glb];
                unsigned int ion_mask_pn_i;
                unsigned int ion_mask_pn_j;
                unsigned int ion_mask_pn_index;
                PMEFloat intdiel_ion_1_p_inv = 1.0 / cSim.intdiel_ion_1_p;
                PMEFloat intdiel_ion_1_n_inv = 1.0 / cSim.intdiel_ion_1_n;
                PMEFloat intdiel_ion_2_pp_inv = 1.0 / cSim.intdiel_ion_2_pp;
                PMEFloat intdiel_ion_2_pn_inv = 1.0 / cSim.intdiel_ion_2_pn;
                PMEFloat intdiel_ion_2_nn_inv = 1.0 / cSim.intdiel_ion_2_nn;
                if (ionmaski == 0){ //atom i is solute
                        ion_mask_pn_i=0;
                }
                else if (ionmaski == 3){ // atom i is anion
                        ion_mask_pn_i=1;
                }
                else { // atom i is cation
                        ion_mask_pn_i=2;
                }
                if (ionmaskj == 0){ //atom j is solute
                        ion_mask_pn_j=0;
                }
                else if (ionmaskj == 3){ // atom j is anion
                        ion_mask_pn_j=1;
                }
                else { // atom j is cation
                        ion_mask_pn_j=2;
                }
                ion_mask_pn_index = ion_mask_pn_i*3+ion_mask_pn_j;
                if (ion_mask_pn_index == 0) { // sulote-solute pair
                        // no change
                }
                else if (ion_mask_pn_index == 1 || ion_mask_pn_index == 3){ // solute-anion
                        dl = cSim.intdiel_inv * intdiel_ion_1_n_inv - expmkf;
                }
                else if (ion_mask_pn_index == 2 || ion_mask_pn_index == 6){ // solute-cation
                        dl = cSim.intdiel_inv * intdiel_ion_1_p_inv - expmkf;
                }
                else if (ion_mask_pn_index == 5 || ion_mask_pn_index == 7){ // cation-anion
                        dl = cSim.intdiel_inv * intdiel_ion_2_pn_inv - expmkf;
                }
                else if (ion_mask_pn_index == 4){ // anion-anion
                        dl = cSim.intdiel_inv * intdiel_ion_2_nn_inv - expmkf;
                }
                else { // (ion_mask_pn_index == 8)  cation-cation
                        dl = cSim.intdiel_inv * intdiel_ion_2_pp_inv - expmkf;
                }
        }

#  ifdef GB_ENERGY
        PMEFloat e = -qiqj * dl * mul;
#    ifndef use_DPFP
        egb += fast_llrintf(ENERGYSCALEF * e);
#    elif defined(GB_MINIMIZATION)
        TLegb += (PMEForce)e;
#    else
        egb += llrint(ENERGYSCALE * e);
#    endif
#  endif
        // 1.0 / fij^3
        PMEFloat temp4 = fgbi * fgbi * fgbi;

        // Here, and in the gas-phase part, "de" contains -(1/r)(dE/dr)
        PMEFloat temp6 = -qiqj * temp4 * (dl + fgbk);
        PMEFloat temp1 = v1;
#  ifdef use_SPFP
        PMEFloat de    = temp6 * (FORCESCALEF - (PMEFloat)0.25 * FORCESCALEF * temp1);
        PMEFloat temp5 = (PMEFloat)0.50 * FORCESCALEF * temp1 * temp6 *
                         (ri * rj + (PMEFloat)0.25 * r2);
#ifdef AMBER_PLATFORM_AMD
        sumdeijda_i += (ri * temp5);
        sumdeijda_j += (rj * temp5);
#else
        sumdeijda_i += fast_llrintf(ri * temp5);
        sumdeijda_j += fast_llrintf(rj * temp5);
#endif
#  else // use_DPFP
        PMEFloat de     = temp6 * ((PMEFloat)1.0 - (PMEFloat)0.25 * temp1);
        PMEFloat temp5  = (PMEFloat)0.50 * temp1 * temp6 * (ri * rj + (PMEFloat)0.25 * r2);
        sumdeijda_i    += (PMEForce)(ri * temp5);
        sumdeijda_j    += (PMEForce)(rj * temp5);
#  endif
#ifdef GB_GBSA3
        // pwsasa: calculate reflect-vdwlike
        if (epsnpi > (PMEFloat)1.0e-8 && PSATOMEPS(j) > (PMEFloat)1.0e-8) {
          PMEFloat dist     = (PMEFloat)1.0 / v5; // 1.0 / rij
          if (dist < (radnpi + PSATOMRAD(j)) ) {
            PMEFloat tempsum = sqrt(epsnpi * PSATOMEPS(j)) / (PMEFloat)6.0;
            PMEFloat numer   = signpi + PSATOMSIG(j);
            PMEFloat numer2  = numer * numer;
            PMEFloat numer4  = numer2 * numer2;
            PMEFloat numer10 = numer4 * numer4 * numer2;
            PMEFloat Sij     = radnpi + PSATOMRAD(j) + numer;
            PMEFloat denom   = (PMEFloat)1.0 / (Sij - dist);
            PMEFloat denom2  = denom * denom;
            PMEFloat denom4  = denom2 * denom2;
            PMEFloat denom10 = denom4 * denom4 * denom2;
            PMEFloat reflectvdwA = (PMEFloat)4.0  * tempsum * numer10 * denom10;
            PMEFloat reflectvdwB = (PMEFloat)10.0 * tempsum * numer4 * denom4;
            PMEFloat AgNPe       = (reflectvdwB - reflectvdwA - sqrt(epsnpi * PSATOMEPS(j))) *
                                   cSim.surften * (PMEFloat)0.6;
            PMEFloat AgNPde      = ((PMEFloat)10.0*reflectvdwA - (PMEFloat)4.0*reflectvdwB) *
                                   denom * v5 * (PMEFloat)1.2 * cSim.surften;
#  ifdef use_SPFP
            de    += FORCESCALEF * AgNPde;
#    ifdef GB_ENERGY
            esurf += fast_llrintf((PMEFloat)2.0 * ENERGYSCALEF * AgNPe);
#    endif // end if GB_ENERGY
#  else // use_DPFP
            de    += AgNPde;
#    ifdef GB_ENERGY
#      ifdef GB_MINIMIZATION
            TLesurf += (PMEFloat)((PMEFloat)2.0 * AgNPe);
#      else
            esurf += llrint((PMEFloat)2.0 * ENERGYSCALE * AgNPe);
#      endif
#    endif //GB_ENERGY
#  endif // use_SPFP
          }
        }
#endif //GB_GBSA3
#else  // GB_IGB6
        PMEFloat de     = 0.0;
#endif // GB_IGB6
        // 1.0 / rij
        PMEFloat rinv      = v5;
        PMEFloat r2inv     = rinv * rinv;
        PMEFloat eel       = cSim.intdiel_inv * qiqj * rinv;

        if (cSim.gbion == 2) { // GB with explict ions
                unsigned int j_glb = y + j;
                unsigned int ionmaski = cSim.pAtomIonMask[i];
                unsigned int ionmaskj = cSim.pAtomIonMask[j_glb];
                PMEFloat intdiel_ion_1_inv = 1.0 / cSim.intdiel_ion_1;
                PMEFloat intdiel_ion_2_inv = 1.0 / cSim.intdiel_ion_2;
                if (ionmaski>0 && ionmaskj>0){ // ion-ion pair
                        eel = cSim.intdiel_inv * intdiel_ion_2_inv * qiqj * rinv;
                }
                else if (ionmaski==0 && ionmaskj==0) { // solute-solute pair
                         // no change
                }
                else { // ion-solute pair
                        eel = cSim.intdiel_inv * intdiel_ion_1_inv * qiqj * rinv;
                }
        }
        else if (cSim.gbion == 3) { // GB with explict ions, separate cation and anion
                unsigned int j_glb = y + j;
                unsigned int ionmaski = cSim.pAtomIonMask[i];
                unsigned int ionmaskj = cSim.pAtomIonMask[j_glb];
                unsigned int ion_mask_pn_i;
                unsigned int ion_mask_pn_j;
                unsigned int ion_mask_pn_index;
                PMEFloat intdiel_ion_1_p_inv = 1.0 / cSim.intdiel_ion_1_p;
                PMEFloat intdiel_ion_1_n_inv = 1.0 / cSim.intdiel_ion_1_n;
                PMEFloat intdiel_ion_2_pp_inv = 1.0 / cSim.intdiel_ion_2_pp;
                PMEFloat intdiel_ion_2_pn_inv = 1.0 / cSim.intdiel_ion_2_pn;
                PMEFloat intdiel_ion_2_nn_inv = 1.0 / cSim.intdiel_ion_2_nn;
                if (ionmaski == 0){ //atom i is solute
                        ion_mask_pn_i=0;
                }
                else if (ionmaski == 3){ // atom i is anion
                        ion_mask_pn_i=1;
                }
                else { // atom i is cation
                        ion_mask_pn_i=2;
                }
                if (ionmaskj == 0){ //atom j is solute
                        ion_mask_pn_j=0;
                }
                else if (ionmaskj == 3){ // atom j is anion
                        ion_mask_pn_j=1;
                }
                else { // atom j is cation
                        ion_mask_pn_j=2;
                }
                ion_mask_pn_index = ion_mask_pn_i*3+ion_mask_pn_j;
                if (ion_mask_pn_index == 0) { // sulote-solute pair
                        // no change
                }
                else if (ion_mask_pn_index == 1 || ion_mask_pn_index == 3){ // solute-anion
                        eel = cSim.intdiel_inv * intdiel_ion_1_n_inv * qiqj * rinv;
                }
                else if (ion_mask_pn_index == 2 || ion_mask_pn_index == 6){ // solute-cation
                        eel = cSim.intdiel_inv * intdiel_ion_1_p_inv * qiqj * rinv;
                }
                else if (ion_mask_pn_index == 5 || ion_mask_pn_index == 7){ // cation-anion
                        eel = cSim.intdiel_inv * intdiel_ion_2_pn_inv * qiqj * rinv;
                }
                else if (ion_mask_pn_index == 4){ // anion-anion
                        eel = cSim.intdiel_inv * intdiel_ion_2_nn_inv * qiqj * rinv;
                }
                else { // (ion_mask_pn_index == 8)  cation-cation
                        eel = cSim.intdiel_inv * intdiel_ion_2_pp_inv * qiqj * rinv;
                }
        }

        PMEFloat r6inv     = r2inv * r2inv * r2inv;
        unsigned int LJIDj = PSATOMLJID(j);
        unsigned int index = LJIDi + LJIDj;
#ifndef use_DPFP
        PMEFloat2 term = tex1Dfetch<float2>(cSim.texLJTerm, index);
#else
        PMEFloat2 term = cSim.pLJTerm[index];
#endif
        PMEFloat f6    = term.y * r6inv;
        PMEFloat f12   = term.x * r6inv * r6inv;
#ifdef use_SPFP // Beginning of pre-processor branch over precision modes
        if (excl & 0x1) {
#  ifdef GB_MINIMIZATION
          de += FORCESCALEF * max(-10000.0f, min(((f12 - f6) + eel) * r2inv, 10000.0f));
#  else
          de += FORCESCALEF * ((f12 - f6) + eel) * r2inv;
#  endif
#  ifdef GB_ENERGY
          eelt += fast_llrintf(ENERGYSCALEF * eel);
          evdw += fast_llrintf(ENERGYSCALEF * ((PMEFloat)(1.0 / 12.0) * f12 -
                                               (PMEFloat)(1.0 / 6.0) * f6));
#  endif
        }
#ifdef AMBER_PLATFORM_AMD
        PMEFloat dedx = de * xij;
        PMEFloat dedy = de * yij;
        PMEFloat dedz = de * zij;
#else
        long long int dedx = fast_llrintf(de * xij);
        long long int dedy = fast_llrintf(de * yij);
        long long int dedz = fast_llrintf(de * zij);
#endif
        fx_i += dedx;
        fy_i += dedy;
        fz_i += dedz;
        PSFX(j) -= dedx;
        PSFY(j) -= dedy;
        PSFZ(j) -= dedz;
#else  // use_DPFP
        if (excl & 0x1) {
#  ifdef GB_ENERGY
#    ifdef GB_MINIMIZATION
          TLeelt += (PMEForce)eel;
          TLevdw += (PMEForce)((PMEFloat)(1.0 / 12) * f12 - (PMEFloat)(1.0 / 6.0) * f6);
#    else
          eelt += llrint(ENERGYSCALE * eel);
          evdw += llrint(ENERGYSCALE * ((PMEFloat)(1.0 / 12.0) * f12 -
                                        (PMEFloat)(1.0 / 6.0) * f6));
#    endif
#  endif
          de += ((f12 - f6) + eel) * r2inv;
        }
        PMEForce dedx = (PMEForce)(de * xij);
        PMEForce dedy = (PMEForce)(de * yij);
        PMEForce dedz = (PMEForce)(de * zij);
        fx_i += dedx;
        fy_i += dedy;
        fz_i += dedz;
        PSFX(j) -= dedx;
        PSFY(j) -= dedy;
        PSFZ(j) -= dedz;
#endif // End of pre-processor branch over precision modes
        __SYNCWARP(WARP_MASK);
        excl >>= 1;
#ifdef AMBER_PLATFORM_AMD
        //
        shAtom.x = WarpRotateLeft<GRID>(shAtom.x);
        shAtom.y = WarpRotateLeft<GRID>(shAtom.y);
        shAtom.z = WarpRotateLeft<GRID>(shAtom.z);
        shAtom.q = WarpRotateLeft<GRID>(shAtom.q);
        shAtom.r = WarpRotateLeft<GRID>(shAtom.r);
        shAtom.LJID = WarpRotateLeft<GRID>(shAtom.LJID);
#ifdef GB_GBSA3
        shAtom.sig = WarpRotateLeft<GRID>(shAtom.sig); //pwsasa
        shAtom.eps = WarpRotateLeft<GRID>(shAtom.eps); //pwsasa
        shAtom.rad = WarpRotateLeft<GRID>(shAtom.rad); //pwsasa
#endif
#if defined(use_SPFP)
        PSFX(j) = WarpRotateLeft<GRID>(PSFX(j));
        PSFY(j) = WarpRotateLeft<GRID>(PSFY(j));
        PSFZ(j) = WarpRotateLeft<GRID>(PSFZ(j));
#endif
#ifndef GB_IGB6
        sumdeijda_j = WarpRotateLeft<GRID>(sumdeijda_j);
#endif //GB_IGB6
        j = ((j + 1) & (GRID - 1));
      }
#else
        shAtom.x    = __SHFL(mask1, shAtom.x, shIdx);
        shAtom.y    = __SHFL(mask1, shAtom.y, shIdx);
        shAtom.z    = __SHFL(mask1, shAtom.z, shIdx);
        shAtom.q    = __SHFL(mask1, shAtom.q, shIdx);
        shAtom.r    = __SHFL(mask1, shAtom.r, shIdx);
        shAtom.LJID = __SHFL(mask1, shAtom.LJID, shIdx);
#ifdef GB_GBSA3
        shAtom.sig  = __SHFL(mask1, shAtom.sig, shIdx); //pwsasa
        shAtom.eps  = __SHFL(mask1, shAtom.eps, shIdx); //pwsasa
        shAtom.rad  = __SHFL(mask1, shAtom.rad, shIdx); //pwsasa
#endif
#ifndef GB_IGB6
        sumdeijda_j = __SHFL(mask1, sumdeijda_j, shIdx);
#endif //GB_IGB6
        j = sNext[j];
        mask1 = __BALLOT(mask1, j != tgx);
      } while (j != tgx);
#endif

      // Write forces.
#ifdef use_SPFP // Here begins another pre-processor branch over precision modes
      int offset = x + tgx;
#ifdef AMBER_PLATFORM_AMD
      atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset],
                llitoulli(fast_llrintf(fx_i)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset],
                llitoulli(fast_llrintf(fy_i)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset],
                llitoulli(fast_llrintf(fz_i)));
#  ifndef GB_IGB6
      atomicAdd((unsigned long long int*)&cSim.pSumdeijdaAccumulator[offset],
                llitoulli(fast_llrintf(sumdeijda_i)));
#  endif
      offset = y + tgx;
      atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset],
                llitoulli(fast_llrintf(PSFX(tgx))));
      atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset],
                llitoulli(fast_llrintf(PSFY(tgx))));
      atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset],
                llitoulli(fast_llrintf(PSFZ(tgx))));
#  ifndef GB_IGB6
      atomicAdd((unsigned long long int*)&cSim.pSumdeijdaAccumulator[offset],
                llitoulli(fast_llrintf(sumdeijda_j)));
#  endif
#else
      atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset], llitoulli(fx_i));
      atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset], llitoulli(fy_i));
      atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset], llitoulli(fz_i));
#  ifndef GB_IGB6
      atomicAdd((unsigned long long int*)&cSim.pSumdeijdaAccumulator[offset],
                llitoulli(sumdeijda_i));
#  endif
      offset = y + tgx;
      atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset],
                llitoulli(PSFX(tgx)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset],
                llitoulli(PSFY(tgx)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset],
                llitoulli(PSFZ(tgx)));
#  ifndef GB_IGB6
      atomicAdd((unsigned long long int*)&cSim.pSumdeijdaAccumulator[offset],
                llitoulli(sumdeijda_j));
#  endif
#endif
#else  // use_DPFP
      int offset = x + tgx;
      //printf("PII %d %20.10f %20.10f %20.10f\n", offset, fx_i, fy_i, fz_i);
#  ifdef GB_MINIMIZATION
      PMEFloat i;
      fx_i = modf(fx_i, &i);
      atomicAdd((unsigned long long int*)&cSim.pIntForceXAccumulator[offset], llitoulli(llrint(i)));
      fy_i = modf(fy_i, &i);
      atomicAdd((unsigned long long int*)&cSim.pIntForceYAccumulator[offset], llitoulli(llrint(i)));
      fz_i = modf(fz_i, &i);
      atomicAdd((unsigned long long int*)&cSim.pIntForceZAccumulator[offset], llitoulli(llrint(i)));
#  endif
      atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset],
                llitoulli(llrint(fx_i * FORCESCALE)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset],
                llitoulli(llrint(fy_i * FORCESCALE)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset],
                llitoulli(llrint(fz_i * FORCESCALE)));
#  ifndef GB_IGB6
      //printf("PII %d %20.10f\n", offset, sumdeijda_i);
      atomicAdd((unsigned long long int*)&cSim.pSumdeijdaAccumulator[offset],
                llitoulli(llrint(sumdeijda_i * FORCESCALE)));
#  endif
      offset = y + tgx;
      //printf("PIJ %d %20.10f %20.10f %20.10f\n", offset, PSFX(tgx), PSFY(tgx), PSFZ(tgx));
#  ifdef GB_MINIMIZATION
      PSFX(tgx) = modf(PSFX(tgx), &i);
      atomicAdd((unsigned long long int*)&cSim.pIntForceXAccumulator[offset], llitoulli(llrint(i)));
      PSFY(tgx) = modf(PSFY(tgx), &i);
      atomicAdd((unsigned long long int*)&cSim.pIntForceYAccumulator[offset], llitoulli(llrint(i)));
      PSFZ(tgx) = modf(PSFZ(tgx), &i);
      atomicAdd((unsigned long long int*)&cSim.pIntForceZAccumulator[offset], llitoulli(llrint(i)));
#  endif
      atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset],
                llitoulli(llrint(PSFX(tgx) * FORCESCALE)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset],
                llitoulli(llrint(PSFY(tgx) * FORCESCALE)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset],
                llitoulli(llrint(PSFZ(tgx) * FORCESCALE)));
#  ifndef GB_IGB6
      //printf("PIJ %d %20.10f\n", offset, sumdeijda_j);
      atomicAdd((unsigned long long int*)&cSim.pSumdeijdaAccumulator[offset],
                llitoulli(llrint(sumdeijda_j * FORCESCALE)));
#  endif
#endif
    }

#if defined(use_DPFP) && defined(GB_ENERGY) && defined(GB_MINIMIZATION)
    PMEFloat f;
    TLevdw = modf(TLevdw, &f);
    evdwe += llrint(f);
    evdw += llrint(TLevdw * ENERGYSCALE);
    TLeelt = modf(TLeelt, &f);
    eelte += llrint(f);
    eelt += llrint(TLeelt * ENERGYSCALE);
    egb += llrint(TLegb * ENERGYSCALE);
# ifdef GB_GBSA3
    esurf += llrint(TLesurf * ENERGYSCALE);
# endif
#endif

#ifndef AMBER_PLATFORM_AMD
    if (tgx == 0) {
      *psPos = atomicAdd(cSim.pGBNB1Position, 1);
    }
#endif
  }
  // Here ends the while loop over work units

#ifdef GB_ENERGY
  // Reduce and write energies
  for (unsigned int stride = warpSize >> 1; stride > 0; stride >>=1) {
    egb += __SHFL_DOWN(WARP_MASK, egb, stride);
    evdw += __SHFL_DOWN(WARP_MASK, evdw, stride);
    eelt += __SHFL_DOWN(WARP_MASK, eelt, stride);
#  ifdef GB_GBSA3  //pwsasa
    esurf += __SHFL_DOWN(WARP_MASK, esurf, stride);
#  endif
#  if defined(use_DPFP) && defined(GB_MINIMIZATION)
    evdwe += __SHFL_DOWN(WARP_MASK, evdwe, stride);
    eelte += __SHFL_DOWN(WARP_MASK, eelte, stride);
#  endif
  }

  // Write out energies
  if ((threadIdx.x & GRID_BITS_MASK) == 0) {
    atomicAdd(cSim.pEGB, llitoulli(egb));
    atomicAdd(cSim.pEVDW, llitoulli(evdw));
    atomicAdd(cSim.pEELT, llitoulli(eelt));
#  ifdef GB_GBSA3
    atomicAdd(cSim.pESurf, llitoulli(esurf)); //pwsasa
#  endif
#  if defined(use_DPFP) && defined(GB_MINIMIZATION)
    atomicAdd(cSim.pEVDWE, llitoulli(evdwe));
    atomicAdd(cSim.pEELTE, llitoulli(eelte));
#  endif
  }

#endif // GB_ENERGY
}
