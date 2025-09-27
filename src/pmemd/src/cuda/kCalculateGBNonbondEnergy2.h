#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------
{
  struct NLForce {
    PMEForce x;
    PMEForce y;
    PMEForce z;
  };

  struct Atom {
    PMEFloat x;
    PMEFloat y;
    PMEFloat z;
    PMEFloat r;
    PMEFloat r1i;
    PMEFloat s;
    PMEFloat temp7;
  };

  struct NLFloat {
    PMEFloat x;
    PMEFloat y;
    PMEFloat z;
  };

#define PSATOMX(i)     shAtom.x
#define PSATOMY(i)     shAtom.y
#define PSATOMZ(i)     shAtom.z
#define PSATOMR(i)     shAtom.r
#define PSATOMR1I(i)   shAtom.r1i
#define PSATOMS(i)     shAtom.s
#define PSATOMTEMP7(i) shAtom.temp7
#define PSFX(i)        psF[i].x
#define PSFY(i)        psF[i].y
#define PSFZ(i)        psF[i].z
#ifdef AMBER_PLATFORM_AMD
#define GBVOLATILE
#else
#define GBVOLATILE     volatile
#endif

  Atom shAtom;
#if defined(use_SPFP) && defined(AMBER_PLATFORM_AMD)
  __shared__ NLFloat sFloat[GBNONBONDENERGY2_THREADS_PER_BLOCK];
#else
  GBVOLATILE __shared__ NLForce sForce[GBNONBONDENERGY2_THREADS_PER_BLOCK];
#endif

#ifdef GB_IGB78
  __shared__ PMEFloat2 sNeckMaxValPos[21 * 21];
#endif
#ifndef AMBER_PLATFORM_AMD
  volatile __shared__ unsigned int sPos[GBNONBONDENERGY2_THREADS_PER_BLOCK / GRID];
  volatile __shared__ unsigned int sNext[GRID];
#endif

#ifdef GB_IGB78
  // Read neck lookup table
  unsigned int pos = threadIdx.x;
  while (pos < 21 * 21) {
    sNeckMaxValPos[pos]  = cSim.pNeckMaxValPos[pos];
    pos                 += blockDim.x;
  }
#endif
#ifndef AMBER_PLATFORM_AMD
  if (threadIdx.x < GRID) {
    sNext[threadIdx.x] = (threadIdx.x + 1) & (GRID - 1);
  }
#endif
  __syncthreads();

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

    // Extract cell coordinates from appropriate work unit
    unsigned int x   = cSim.pWorkUnit[ppos];
    unsigned int y   = ((x >> 2) & 0x7fff) << GRID_BITS;
    x                = (x >> 17) << GRID_BITS;
    unsigned int tgx = threadIdx.x & (GRID - 1);
    unsigned int i   = x + tgx;
    PMEFloat2 xyi    = cSim.pAtomXYSP[i];
    PMEFloat zi      = cSim.pAtomZSP[i];
    PMEFloat ri      = cSim.pAtomRBorn[i];
    PMEFloat si      = cSim.pAtomS[i];
    PMEFloat temp7_i = cSim.pTemp7[i];
    unsigned int tbx  = threadIdx.x - tgx;
#if defined(use_SPFP) && defined(AMBER_PLATFORM_AMD)
    NLFloat* psF = &sFloat[tbx];
    PMEFloat fx_i    = (PMEFloat)0;
    PMEFloat fy_i    = (PMEFloat)0;
    PMEFloat fz_i    = (PMEFloat)0;
    PSFX(tgx) = (PMEFloat)0;
    PSFY(tgx) = (PMEFloat)0;
    PSFZ(tgx) = (PMEFloat)0;
#else
    GBVOLATILE NLForce* psF = &sForce[tbx];
    PMEForce fx_i    = (PMEForce)0;
    PMEForce fy_i    = (PMEForce)0;
    PMEForce fz_i    = (PMEForce)0;
    PSFX(tgx) = (PMEForce)0;
    PSFY(tgx) = (PMEForce)0;
    PSFZ(tgx) = (PMEForce)0;
#endif
    __SYNCWARP(WARP_MASK);
#ifdef AMBER_PLATFORM_AMD
    unsigned int shIdx = ((tgx + 1) & (GRID - 1));
#else
    unsigned int shIdx = sNext[tgx];
#endif

    // Handle diagonals uniquely at 50% efficiency, skipping i == j interactions
    if (x == y) {
      PMEFloat xi  = xyi.x;
      PMEFloat yi  = xyi.y;
      PSATOMX(tgx) = xi;
      PSATOMY(tgx) = yi;
      PSATOMZ(tgx) = zi;
      ri             -= cSim.offset;
      PMEFloat ri1i   = (PMEFloat)1.0 / ri;
      PSATOMR(tgx)    = ri;
      PSATOMR1I(tgx)  = ri1i;
      PSATOMS(tgx)    = si;

#ifdef AMBER_PLATFORM_AMD
      shAtom.x = WarpRotateLeft<GRID>(shAtom.x);
      shAtom.y = WarpRotateLeft<GRID>(shAtom.y);
      shAtom.z = WarpRotateLeft<GRID>(shAtom.z);
      shAtom.r = WarpRotateLeft<GRID>(shAtom.r);
      shAtom.r1i = WarpRotateLeft<GRID>(shAtom.r1i);
      shAtom.s = WarpRotateLeft<GRID>(shAtom.s);
      unsigned int j = ((tgx + 1) & (GRID - 1));
#else
      shAtom.x   = __SHFL(WARP_MASK, shAtom.x, shIdx);
      shAtom.y   = __SHFL(WARP_MASK, shAtom.y, shIdx);
      shAtom.z   = __SHFL(WARP_MASK, shAtom.z, shIdx);
      shAtom.r   = __SHFL(WARP_MASK, shAtom.r, shIdx);
      shAtom.r1i = __SHFL(WARP_MASK, shAtom.r1i, shIdx);
      shAtom.s   = __SHFL(WARP_MASK, shAtom.s, shIdx);
      unsigned int j = sNext[tgx];
#endif
      PMEMask mask1 = __BALLOT(WARP_MASK, j != tgx);
      while (j != tgx) {
        PMEFloat xij = xi - PSATOMX(j);
        PMEFloat yij = yi - PSATOMY(j);
        PMEFloat zij = zi - PSATOMZ(j);
        PMEFloat r2  = xij*xij + yij*yij + zij*zij;
        PMEFloat dij = sqrt(r2);
        PMEFloat sj  = PSATOMS(j);
        if (dij <= cSim.rgbmax + sj) {
          PMEFloat dij1i = (PMEFloat)1.0 / dij;
          PMEFloat dij2i = dij1i * dij1i;
          PMEFloat dij3i = dij2i * dij1i;
          PMEFloat datmp = (PMEFloat)0.0;
          PMEFloat v3    = (PMEFloat)1.0 / (dij + sj);
          if (dij > cSim.rgbmax - sj) {
            PMEFloat temp1 = (PMEFloat)1.0 / (dij - sj);
            datmp = (PMEFloat)0.125 * dij3i * ((r2 + sj*sj)*(temp1*temp1 - cSim.rgbmax2i) -
                                               (PMEFloat)2.0 * log(cSim.rgbmax * temp1));
          }
          else if (dij > (PMEFloat)4.0 * sj) {
            PMEFloat tmpsd = sj * sj * dij2i;
            PMEFloat dumbo = te + tmpsd*(tf + tmpsd*(tg + tmpsd*(th + tmpsd*thh)));
            datmp          = tmpsd * sj * dij2i * dij2i * dumbo;
          }
          else if (dij > ri + sj) {
            PMEFloat v2 = (PMEFloat)1.0 / (r2 - sj*sj);
            PMEFloat v4 = log(v3 * (dij - sj));
            datmp       = v2*sj*((PMEFloat)-0.5 * dij2i + v2) + (PMEFloat)0.25*dij3i*v4;
          }
          else if (dij > abs(ri - sj)) {
            PMEFloat v2 = (PMEFloat)1.0 / (dij + sj);
            PMEFloat v4 = log(v3 * ri);
            datmp       = (PMEFloat)-0.25 *
                          ((PMEFloat)-0.5*(r2 - ri*ri + sj*sj)*dij3i*ri1i*ri1i +
                           dij1i*v2*(v2 - dij1i) - dij3i*v4);
          }
          else if (ri < sj) {
            PMEFloat v2 = (PMEFloat)1.0 / (r2 - sj*sj);
            PMEFloat v4 = log(v3 * (sj - dij));
            datmp       = (PMEFloat)-0.5 * (sj*dij2i*v2 - (PMEFloat)2.0*sj*v2*v2 -
                                            (PMEFloat)0.5*dij3i*v4);
          }
#ifdef GB_IGB78
          if (dij < ri + PSATOMR(j) + cSim.gb_neckcut) {
            unsigned int ii      = round((ri - cSim.gb_neckoffset) * (PMEFloat)20.0);
            unsigned int jj      = round((PSATOMR(j) - cSim.gb_neckoffset) * (PMEFloat)20.0);
            PMEFloat2 neckValPos = sNeckMaxValPos[ii * 21 + jj];
            PMEFloat mdist       = dij - neckValPos.y;
            PMEFloat mdist2      = mdist * mdist;
            PMEFloat mdist3      = mdist2 * mdist;
            PMEFloat mdist5      = mdist2 * mdist3;
            PMEFloat mdist6      = mdist3 * mdist3;
            PMEFloat temp1       = (PMEFloat)1.0 + mdist2 + (PMEFloat)0.3 * mdist6;
            temp1               *= temp1 * dij;
            // GB with explicit ions gbion==2
            PMEFloat gb_neckscale_temp1 = cSim.gb_neckscale;
            if (cSim.gbion==2) {
              unsigned int j_glb = y + j;
              unsigned int ionmaski = cSim.pAtomIonMask[i];
              unsigned int ionmaskj = cSim.pAtomIonMask[j_glb];
              if (ionmaski>0 && ionmaskj>0){ // ion-ion pair
                gb_neckscale_temp1 = cSim.gb_neckscale * cSim.gb_neckscale_ion_2;
              }
              else if (ionmaski==0 && ionmaskj==0) { // solute-solute pair
                gb_neckscale_temp1 = cSim.gb_neckscale; // no change
              }
              else { // ion-solute pair
                gb_neckscale_temp1 = cSim.gb_neckscale * cSim.gb_neckscale_ion_1;
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
                      gb_neckscale_temp1 = cSim.gb_neckscale; // no change
              }
              else if (ion_mask_pn_index == 1 || ion_mask_pn_index == 3){ // solute-anion
                      gb_neckscale_temp1 = cSim.gb_neckscale * cSim.gb_neckscale_ion_1_n;
              }
              else if (ion_mask_pn_index == 2 || ion_mask_pn_index == 6){ // solute-cation
                      gb_neckscale_temp1 = cSim.gb_neckscale * cSim.gb_neckscale_ion_1_p;
              }
              else if (ion_mask_pn_index == 5 || ion_mask_pn_index == 7){ // cation-anion
                      gb_neckscale_temp1 = cSim.gb_neckscale * cSim.gb_neckscale_ion_2_pn;
              }
              else if (ion_mask_pn_index == 4){ // anion-anion
                      gb_neckscale_temp1 = cSim.gb_neckscale * cSim.gb_neckscale_ion_2_nn;
              }
              else { // (ion_mask_pn_index == 8)  cation-cation
                      gb_neckscale_temp1 = cSim.gb_neckscale * cSim.gb_neckscale_ion_2_pp;
              }
            }

            datmp               += (((PMEFloat)2.0*mdist + (PMEFloat)(9.0/5.0)*mdist5) *
                                    neckValPos.x * gb_neckscale_temp1) / temp1;
            
            // datmp               += (((PMEFloat)2.0*mdist + (PMEFloat)(9.0/5.0)*mdist5) *
            //                         neckValPos.x * cSim.gb_neckscale) / temp1;
          }
#endif
          datmp *= -temp7_i;
#ifdef use_SPFP
#ifdef AMBER_PLATFORM_AMD
          PMEFloat f_x = (xij * datmp);
          PMEFloat f_y = (yij * datmp);
          PMEFloat f_z = (zij * datmp);
#else
          long long int f_x  = fast_llrintf(FORCESCALEF * xij * datmp);
          long long int f_y  = fast_llrintf(FORCESCALEF * yij * datmp);
          long long int f_z  = fast_llrintf(FORCESCALEF * zij * datmp);
#endif
          fx_i -= f_x;
          fy_i -= f_y;
          fz_i -= f_z;
          PSFX(j) += f_x;
          PSFY(j) += f_y;
          PSFZ(j) += f_z;
#else // use_DPFP
          PMEForce f_x  = (PMEForce)(xij * datmp);
          PSFX(j)      += f_x;
          fx_i         -= f_x;
          PMEForce f_y  = (PMEForce)(yij * datmp);
          PSFY(j)      += f_y;
          fy_i         -= f_y;
          PMEForce f_z  = (PMEForce)(zij * datmp);
          PSFZ(j)      += f_z;
          fz_i         -= f_z;
#endif // End pre-processor branch over precision modes
        }
        __SYNCWARP(WARP_MASK);
#ifdef AMBER_PLATFORM_AMD
        shAtom.x = WarpRotateLeft<GRID>(shAtom.x);
        shAtom.y = WarpRotateLeft<GRID>(shAtom.y);
        shAtom.z = WarpRotateLeft<GRID>(shAtom.z);
        shAtom.r = WarpRotateLeft<GRID>(shAtom.r);
        shAtom.r1i = WarpRotateLeft<GRID>(shAtom.r1i);
        shAtom.s = WarpRotateLeft<GRID>(shAtom.s);
        j = ((j + 1) & (GRID - 1));
#else
        shAtom.x   = __SHFL(mask1, shAtom.x, shIdx);
        shAtom.y   = __SHFL(mask1, shAtom.y, shIdx);
        shAtom.z   = __SHFL(mask1, shAtom.z, shIdx);
        shAtom.r   = __SHFL(mask1, shAtom.r, shIdx);
        shAtom.r1i = __SHFL(mask1, shAtom.r1i, shIdx);
        shAtom.s   = __SHFL(mask1, shAtom.s, shIdx);
        j = sNext[j];
        mask1 = __BALLOT(mask1, j != tgx);
#endif
      }
      // Here ends a loop over j that has the non-bonded calculation hopping around
      // in a list of interactions

      int offset = x + tgx;
      fx_i += PSFX(tgx);
      fy_i += PSFY(tgx);
      fz_i += PSFZ(tgx);
#ifdef use_SPFP
#ifdef AMBER_PLATFORM_AMD
      atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset],
                llitoulli(fast_llrintf(FORCESCALEF * fx_i)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset],
                llitoulli(fast_llrintf(FORCESCALEF * fy_i)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset],
                llitoulli(fast_llrintf(FORCESCALEF * fz_i)));
#else
      atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset], llitoulli(fx_i));
      atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset], llitoulli(fy_i));
      atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset], llitoulli(fz_i));
#endif
#else  // use_DPFP
      atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset],
                llitoulli(llrint(fx_i * FORCESCALE)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset],
                llitoulli(llrint(fy_i * FORCESCALE)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset],
                llitoulli(llrint(fz_i * FORCESCALE)));
#endif // End pre-processor branch over precision modes
    }
    else {
      unsigned int j   = y + tgx;
      PMEFloat2 xyj    = cSim.pAtomXYSP[j];
      PSATOMZ(tgx)     = cSim.pAtomZSP[j];
      PSATOMR(tgx)     = cSim.pAtomRBorn[j];
      PSATOMS(tgx)     = cSim.pAtomS[j];
      PSATOMTEMP7(tgx) = cSim.pTemp7[j];
      PMEFloat xi = xyi.x;
      PMEFloat yi = xyi.y;
      ri             -= cSim.offset;
      PMEFloat ri1i   = (PMEFloat)1.0 / ri;
      PMEFloat si2    = si * si;
      PSATOMX(tgx)    = xyj.x;
      PSATOMY(tgx)    = xyj.y;
      PSATOMR(tgx)   -= cSim.offset;
      PSATOMR1I(tgx)  = (PMEFloat)1.0 / PSATOMR(tgx);
      j               = tgx;
      PMEMask mask1 = WARP_MASK;
#ifdef AMBER_PLATFORM_AMD
      for (unsigned int xx = 0; xx < GRID; ++xx) {
#else
      do {
#endif
        PMEFloat xij   = xi - PSATOMX(j);
        PMEFloat yij   = yi - PSATOMY(j);
        PMEFloat zij   = zi - PSATOMZ(j);
        PMEFloat r2    = xij*xij + yij*yij + zij*zij;
        PMEFloat dij   = sqrt(r2);
        PMEFloat dij1i = (PMEFloat)1.0 / dij;
        PMEFloat dij2i = dij1i * dij1i;
        PMEFloat dij3i = dij2i * dij1i;

        // Atom i forces
        PMEFloat sj = PSATOMS(j);
        if (dij <= cSim.rgbmax + sj) {
          PMEFloat datmp = (PMEFloat)0.0;
          PMEFloat v3    = (PMEFloat)1.0 / (dij + sj);
          if (dij > cSim.rgbmax - sj) {
            PMEFloat temp1 = (PMEFloat)1.0 / (dij - sj);
            datmp          = (PMEFloat)0.125 * dij3i *
                             ((r2 + sj*sj)*(temp1*temp1 - cSim.rgbmax2i) -
                              (PMEFloat)2.0 * log(cSim.rgbmax * temp1));
          }
          else if (dij > (PMEFloat)4.0 * sj) {
            PMEFloat tmpsd = sj * sj * dij2i;
            PMEFloat dumbo = te + tmpsd*(tf + tmpsd*(tg + tmpsd*(th + tmpsd*thh)));
            datmp          = tmpsd * sj * dij2i * dij2i * dumbo;
          }
          else if (dij > ri + sj) {
            PMEFloat v2         = (PMEFloat)1.0 / (r2 - sj*sj);
            PMEFloat v4         = log(v3 * (dij - sj));
            datmp               = v2*sj*((PMEFloat)-0.5*dij2i + v2) +
                                  (PMEFloat)0.25*dij3i*v4;
          }
          else if (dij > abs(ri - sj)) {
            PMEFloat v2 = (PMEFloat)1.0 / (dij + sj);
            PMEFloat v4 = log(v3 * ri);
            datmp    = (PMEFloat)-0.25 *
                       ((PMEFloat)-0.5*(r2 - ri*ri + sj*sj)*dij3i*ri1i*ri1i +
                        dij1i*v2*(v2 - dij1i) - dij3i*v4);
          }
          else if (ri < sj) {
            PMEFloat v2 = (PMEFloat)1.0 / (r2 - sj*sj);
            PMEFloat v4 = log(v3 * (sj - dij));
            datmp = (PMEFloat)-0.5 * (sj*dij2i*v2 - (PMEFloat)2.0*sj*v2*v2 -
                                      (PMEFloat)0.5*dij3i*v4);
          }
#ifdef GB_IGB78
          if (dij < ri + PSATOMR(j) + cSim.gb_neckcut) {
            unsigned int ii      = round((ri - cSim.gb_neckoffset) * (PMEFloat)20.0);
            unsigned int jj      = round((PSATOMR(j) - cSim.gb_neckoffset) * (PMEFloat)20.0);
            PMEFloat2 neckValPos = sNeckMaxValPos[ii*21 + jj];
            PMEFloat mdist       = dij - neckValPos.y;
            PMEFloat mdist2      = mdist * mdist;
            PMEFloat mdist3      = mdist2 * mdist;
            PMEFloat mdist5      = mdist2 * mdist3;
            PMEFloat mdist6      = mdist3 * mdist3;
            PMEFloat temp1  = (PMEFloat)1.0 + mdist2 + (PMEFloat)0.3*mdist6;
            temp1          *= temp1 * dij;

            // GB with explicit ions gbion==2
            PMEFloat gb_neckscale_temp1 = cSim.gb_neckscale;
            if (cSim.gbion==2) {
              unsigned int j_glb = y + j;
              unsigned int ionmaski = cSim.pAtomIonMask[i];
              unsigned int ionmaskj = cSim.pAtomIonMask[j_glb];
              if (ionmaski>0 && ionmaskj>0){ // ion-ion pair
                gb_neckscale_temp1 = cSim.gb_neckscale * cSim.gb_neckscale_ion_2;
              }
              else if (ionmaski==0 && ionmaskj==0) { // solute-solute pair
                gb_neckscale_temp1 = cSim.gb_neckscale; // no change
              }
              else { // ion-solute pair
                gb_neckscale_temp1 = cSim.gb_neckscale * cSim.gb_neckscale_ion_1;
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
                      gb_neckscale_temp1 = cSim.gb_neckscale; // no change
              }
              else if (ion_mask_pn_index == 1 || ion_mask_pn_index == 3){ // solute-anion
                      gb_neckscale_temp1 = cSim.gb_neckscale * cSim.gb_neckscale_ion_1_n;
              }
              else if (ion_mask_pn_index == 2 || ion_mask_pn_index == 6){ // solute-cation
                      gb_neckscale_temp1 = cSim.gb_neckscale * cSim.gb_neckscale_ion_1_p;
              }
              else if (ion_mask_pn_index == 5 || ion_mask_pn_index == 7){ // cation-anion
                      gb_neckscale_temp1 = cSim.gb_neckscale * cSim.gb_neckscale_ion_2_pn;
              }
              else if (ion_mask_pn_index == 4){ // anion-anion
                      gb_neckscale_temp1 = cSim.gb_neckscale * cSim.gb_neckscale_ion_2_nn;
              }
              else { // (ion_mask_pn_index == 8)  cation-cation
                      gb_neckscale_temp1 = cSim.gb_neckscale * cSim.gb_neckscale_ion_2_pp;
              }
            }
            datmp          += (((PMEFloat)2.0*mdist + (PMEFloat)(9.0/5.0)*mdist5) *
                               neckValPos.x * gb_neckscale_temp1) / temp1;
            // datmp          += (((PMEFloat)2.0*mdist + (PMEFloat)(9.0/5.0)*mdist5) *
            //                    neckValPos.x*cSim.gb_neckscale) / temp1;
          }
#endif
          datmp *= -temp7_i;
#ifdef use_SPFP
#ifdef AMBER_PLATFORM_AMD
          PMEFloat f_x = (xij * datmp);
          PMEFloat f_y = (yij * datmp);
          PMEFloat f_z = (zij * datmp);
#else
          long long int f_x  = fast_llrintf(FORCESCALEF * xij * datmp);
          long long int f_y  = fast_llrintf(FORCESCALEF * yij * datmp);
          long long int f_z  = fast_llrintf(FORCESCALEF * zij * datmp);
#endif
          fx_i -= f_x;
          fy_i -= f_y;
          fz_i -= f_z;
          PSFX(j) += f_x;
          PSFY(j) += f_y;
          PSFZ(j) += f_z;
#else  // use_DPFP
          PMEForce f_x  = (PMEForce)(xij * datmp);
          PSFX(j)      += f_x;
          fx_i         -= f_x;
          PMEForce f_y  = (PMEForce)(yij * datmp);
          PSFY(j)      += f_y;
          fy_i         -= f_y;
          PMEForce f_z  = (PMEForce)(zij * datmp);
          PSFZ(j)      += f_z;
          fz_i         -= f_z;
#endif // End pre-processor branch over precision modes
        }

        // Atom j forces
        if (dij <= cSim.rgbmax + si) {
          PMEFloat datmp = (PMEFloat)0.0;
          PMEFloat v3    = (PMEFloat)1.0 / (dij + si);
          if (dij > cSim.rgbmax - si) {
            PMEFloat temp1 = (PMEFloat)1.0 / (dij - si);
            datmp = (PMEFloat)0.125 * dij3i * ((r2 + si2)*(temp1*temp1 - cSim.rgbmax2i) -
                                               (PMEFloat)2.0*log(cSim.rgbmax * temp1));
          }
          else if (dij > (PMEFloat)4.0 * si) {
            PMEFloat tmpsd = si2 * dij2i;
            PMEFloat dumbo = te + tmpsd*(tf + tmpsd*(tg + tmpsd*(th + tmpsd*thh)));
            datmp          = tmpsd * si * dij2i * dij2i * dumbo;
          }
          else if (dij > PSATOMR(j) + si) {
            PMEFloat v2 = (PMEFloat)1.0 / (r2 - si2);
            PMEFloat v4 = log(v3 * (dij - si));
            datmp       = v2*si*((PMEFloat)-0.5*dij2i + v2) + (PMEFloat)0.25*dij3i*v4;
          }
          else if (dij > abs(PSATOMR(j) - si)) {
            PMEFloat v2 = (PMEFloat)1.0 / (dij + si);
            PMEFloat v4 = log(v3 * PSATOMR(j));
            datmp       = (PMEFloat)-0.25 *
                          ((PMEFloat)-0.5 * (r2 - PSATOMR(j)*PSATOMR(j) + si2) *
                           dij3i*PSATOMR1I(j)*PSATOMR1I(j) + dij1i*v2*(v2 - dij1i) - dij3i*v4);
          }
          else if (PSATOMR(j) < si) {
            PMEFloat v2 = (PMEFloat)1.0 / (r2 - si2);
            PMEFloat v4 = log(v3 * (si - dij));
            datmp       = (PMEFloat)-0.5 * (si*dij2i*v2 - (PMEFloat)2.0*si*v2*v2 -
                                            (PMEFloat)0.5*dij3i*v4);
          }
#ifdef GB_IGB78
          if (dij < ri + PSATOMR(j) + cSim.gb_neckcut) {
            unsigned int ii      = round((ri - cSim.gb_neckoffset) * (PMEFloat)20.0);
            unsigned int jj      = round((PSATOMR(j) - cSim.gb_neckoffset) * (PMEFloat)20.0);
            PMEFloat2 neckValPos = sNeckMaxValPos[jj*21 + ii];
            PMEFloat mdist       = dij - neckValPos.y;
            PMEFloat mdist2      = mdist * mdist;
            PMEFloat mdist3      = mdist2 * mdist;
            PMEFloat mdist5      = mdist2 * mdist3;
            PMEFloat mdist6      = mdist3 * mdist3;
            PMEFloat temp1       = (PMEFloat)1.0 + mdist2 + (PMEFloat)0.3*mdist6;
            temp1               *= temp1 * dij;
            // GB with explicit ions gbion==2
            PMEFloat gb_neckscale_temp1 = cSim.gb_neckscale;
            if (cSim.gbion==2) {
              unsigned int j_glb = y + j;
              unsigned int ionmaski = cSim.pAtomIonMask[i];
              unsigned int ionmaskj = cSim.pAtomIonMask[j_glb];
              if (ionmaski>0 && ionmaskj>0){ // ion-ion pair
                gb_neckscale_temp1 = cSim.gb_neckscale * cSim.gb_neckscale_ion_2;
              }
              else if (ionmaski==0 && ionmaskj==0) { // solute-solute pair
                gb_neckscale_temp1 = cSim.gb_neckscale; // no change
              }
              else { // ion-solute pair
                gb_neckscale_temp1 = cSim.gb_neckscale * cSim.gb_neckscale_ion_1;
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
                      gb_neckscale_temp1 = cSim.gb_neckscale; // no change
              }
              else if (ion_mask_pn_index == 1 || ion_mask_pn_index == 3){ // solute-anion
                      gb_neckscale_temp1 = cSim.gb_neckscale * cSim.gb_neckscale_ion_1_n;
              }
              else if (ion_mask_pn_index == 2 || ion_mask_pn_index == 6){ // solute-cation
                      gb_neckscale_temp1 = cSim.gb_neckscale * cSim.gb_neckscale_ion_1_p;
              }
              else if (ion_mask_pn_index == 5 || ion_mask_pn_index == 7){ // cation-anion
                      gb_neckscale_temp1 = cSim.gb_neckscale * cSim.gb_neckscale_ion_2_pn;
              }
              else if (ion_mask_pn_index == 4){ // anion-anion
                      gb_neckscale_temp1 = cSim.gb_neckscale * cSim.gb_neckscale_ion_2_nn;
              }
              else { // (ion_mask_pn_index == 8)  cation-cation
                      gb_neckscale_temp1 = cSim.gb_neckscale * cSim.gb_neckscale_ion_2_pp;
              }
            }
            datmp               += (((PMEFloat)2.0*mdist + (PMEFloat)(9.0/5.0)*mdist5) *
                                    neckValPos.x * gb_neckscale_temp1) / temp1;
            // datmp               += (((PMEFloat)2.0*mdist + (PMEFloat)(9.0/5.0)*mdist5) *
            //                         neckValPos.x * cSim.gb_neckscale) / temp1;
          }
#endif
          datmp *= PSATOMTEMP7(j);
#ifdef use_SPFP
#ifdef AMBER_PLATFORM_AMD
          PMEFloat f_x = (xij * datmp);
          PMEFloat f_y = (yij * datmp);
          PMEFloat f_z = (zij * datmp);
#else
          long long int f_x  = fast_llrintf(FORCESCALEF * xij * datmp);
          long long int f_y  = fast_llrintf(FORCESCALEF * yij * datmp);
          long long int f_z  = fast_llrintf(FORCESCALEF * zij * datmp);
#endif
          fx_i += f_x;
          fy_i += f_y;
          fz_i += f_z;
          PSFX(j) -= f_x;
          PSFY(j) -= f_y;
          PSFZ(j) -= f_z;
#else // use_DPFP
          PMEForce f_x  = (PMEForce)(xij * datmp);
          PSFX(j)      -= f_x;
          fx_i         += f_x;
          PMEForce f_y  = (PMEForce)(yij * datmp);
          PSFY(j)      -= f_y;
          fy_i         += f_y;
          PMEForce f_z  = (PMEForce)(zij * datmp);
          PSFZ(j)      -= f_z;
          fz_i         += f_z;
#endif
        }
        __SYNCWARP(WARP_MASK);

#ifdef AMBER_PLATFORM_AMD
        shAtom.x = WarpRotateLeft<GRID>(shAtom.x);
        shAtom.y = WarpRotateLeft<GRID>(shAtom.y);
        shAtom.z = WarpRotateLeft<GRID>(shAtom.z);
        shAtom.r = WarpRotateLeft<GRID>(shAtom.r);
        shAtom.r1i = WarpRotateLeft<GRID>(shAtom.r1i);
        shAtom.s = WarpRotateLeft<GRID>(shAtom.s);
        shAtom.temp7 = WarpRotateLeft<GRID>(shAtom.temp7);
        j = ((j + 1) & (GRID - 1));
      }
#else
        shAtom.x     = __SHFL(mask1, shAtom.x, shIdx);
        shAtom.y     = __SHFL(mask1, shAtom.y, shIdx);
        shAtom.z     = __SHFL(mask1, shAtom.z, shIdx);
        shAtom.r     = __SHFL(mask1, shAtom.r, shIdx);
        shAtom.r1i   = __SHFL(mask1, shAtom.r1i, shIdx);
        shAtom.s     = __SHFL(mask1, shAtom.s, shIdx);
        shAtom.temp7 = __SHFL(mask1, shAtom.temp7, shIdx);
        j = sNext[j];
        mask1 = __BALLOT(mask1, j != tgx);
      } while (j != tgx);
#endif
            // Write forces
#ifdef use_SPFP
      int offset = x + tgx;
#ifdef AMBER_PLATFORM_AMD
      atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset],
                llitoulli(fast_llrintf(FORCESCALEF * fx_i)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset],
                llitoulli(fast_llrintf(FORCESCALEF * fy_i)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset],
                llitoulli(fast_llrintf(FORCESCALEF * fz_i)));
      offset     = y + tgx;
      atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset],
                llitoulli(fast_llrintf(FORCESCALEF * PSFX(tgx))));
      atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset],
                llitoulli(fast_llrintf(FORCESCALEF * PSFY(tgx))));
      atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset],
                llitoulli(fast_llrintf(FORCESCALEF * PSFZ(tgx))));
#else
      atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset], llitoulli(fx_i));
      atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset], llitoulli(fy_i));
      atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset], llitoulli(fz_i));
      offset     = y + tgx;
      atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset],
                llitoulli(PSFX(tgx)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset],
                llitoulli(PSFY(tgx)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset],
                llitoulli(PSFZ(tgx)));
#endif
#else  // use_DPFP
      int offset = x + tgx;
      atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset],
                llitoulli(llrint(fx_i * FORCESCALE)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset],
                llitoulli(llrint(fy_i * FORCESCALE)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset],
                llitoulli(llrint(fz_i * FORCESCALE)));
      offset     = y + tgx;
      atomicAdd((unsigned long long int*)&cSim.pNBForceXAccumulator[offset],
                llitoulli(llrint(PSFX(tgx) * FORCESCALE)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceYAccumulator[offset],
                llitoulli(llrint(PSFY(tgx) * FORCESCALE)));
      atomicAdd((unsigned long long int*)&cSim.pNBForceZAccumulator[offset],
                llitoulli(llrint(PSFZ(tgx) * FORCESCALE)));
#endif
    }
#ifndef AMBER_PLATFORM_AMD
    if (tgx == 0) {
      *psPos = atomicAdd(cSim.pGBNB2Position, 1);
    }
#endif
  }
}
