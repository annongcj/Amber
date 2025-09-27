#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

// This file is meant to be included in kCalculateGBBornRadii.cu, twice
{
struct Atom {
  PMEFloat x;
  PMEFloat y;
  PMEFloat z;
  PMEFloat r;
  PMEFloat r1i;
  PMEFloat s;
  PMEFloat s2;
};

#define PSATOMX(i)   shAtom.x
#define PSATOMY(i)   shAtom.y
#define PSATOMZ(i)   shAtom.z
#define PSATOMR(i)   shAtom.r
#define PSATOMR1I(i) shAtom.r1i
#define PSATOMS(i)   shAtom.s
#define PSATOMS2(i)  shAtom.s2
#ifdef AMBER_PLATFORM_AMD
#ifdef use_SPFP
#define PSREFF(i)    psReff
#else
#define PSREFF(i)    psReff[i]
#endif
#define GBVOLATILE
#else
#define PSREFF(i)    psReff[i]
#define GBVOLATILE   volatile
#endif

Atom shAtom;
#if defined(use_SPFP) && defined(AMBER_PLATFORM_AMD)
PMEFloat psReff;
#else
GBVOLATILE __shared__ PMEForce sReff[GBBORNRADII_THREADS_PER_BLOCK];
#endif

#ifdef GB_IGB78
__shared__ PMEFloat2 sNeckMaxValPos[21 * 21];
#endif
#ifndef AMBER_PLATFORM_AMD
volatile __shared__ unsigned int sPos[GBBORNRADII_THREADS_PER_BLOCK / GRID];
volatile __shared__ unsigned int sNext[GRID];
#endif
#ifdef GB_IGB78
  // Read neck lookup table
  unsigned int pos = threadIdx.x;
  while (pos < 21 * 21) {
    sNeckMaxValPos[pos] = cSim.pNeckMaxValPos[pos];
    pos += blockDim.x;
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
    PMEFloat si2     = si * si;
    unsigned int tbx = threadIdx.x - tgx;
#ifdef AMBER_PLATFORM_AMD
    unsigned int shIdx = ((tgx + 1) & (GRID - 1));
#else
    unsigned int shIdx = sNext[tgx];
#endif
#if defined(use_SPFP) && defined(AMBER_PLATFORM_AMD)
    PMEFloat reff_i  = (PMEFloat)0;
    PSREFF(tgx)      = (PMEFloat)0;
#else
    GBVOLATILE PMEForce* psReff = &sReff[tbx];
    PMEForce reff_i  = (PMEForce)0;
    PSREFF(tgx)      = (PMEForce)0;
    __SYNCWARP(WARP_MASK);
#endif

    // Handle diagonals uniquely at 50% efficiency, skipping i == j interactions
    // x and y are always consistent within a warp, no diverge here
    if (x == y) {
      PMEFloat xi = xyi.x;
      PMEFloat yi = xyi.y;
      PSATOMX(tgx) = xi;
      PSATOMY(tgx) = yi;
      PSATOMZ(tgx) = zi;
      ri -= cSim.offset;
      PMEFloat ri1i = (PMEFloat)1.0 / ri;
      PSATOMR(tgx) = ri;
      PSATOMS(tgx) = si;
      PSATOMS2(tgx) = si2;
      PSATOMR1I(tgx) = ri1i;

#ifdef AMBER_PLATFORM_AMD
      shAtom.x = WarpRotateLeft<GRID>(shAtom.x);
      shAtom.y = WarpRotateLeft<GRID>(shAtom.y);
      shAtom.z = WarpRotateLeft<GRID>(shAtom.z);
      shAtom.r = WarpRotateLeft<GRID>(shAtom.r);
      shAtom.s = WarpRotateLeft<GRID>(shAtom.s);
      shAtom.s2 = WarpRotateLeft<GRID>(shAtom.s2);
      shAtom.r1i = WarpRotateLeft<GRID>(shAtom.r1i);
      unsigned int j = ((tgx + 1) & (GRID - 1));
#else
      shAtom.x   = __SHFL(WARP_MASK, shAtom.x, shIdx);
      shAtom.y   = __SHFL(WARP_MASK, shAtom.y, shIdx);
      shAtom.z   = __SHFL(WARP_MASK, shAtom.z, shIdx);
      shAtom.r   = __SHFL(WARP_MASK, shAtom.r, shIdx);
      shAtom.s   = __SHFL(WARP_MASK, shAtom.s, shIdx);
      shAtom.s2  = __SHFL(WARP_MASK, shAtom.s2, shIdx);
      shAtom.r1i = __SHFL(WARP_MASK, shAtom.r1i, shIdx);
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
        if (dij < cSim.rgbmax + sj) {
          PMEFloat dij1i = (PMEFloat)1.0 / dij;
          PMEFloat dij2i = dij1i * dij1i;
          PMEFloat dr;
          if (dij > cSim.rgbmax - sj) {
            PMEFloat uij = (PMEFloat)1.0 / (dij - sj);
            dr = (PMEFloat)0.125*dij1i*((PMEFloat)1.0 + (PMEFloat)2.0*dij*uij +
                 cSim.rgbmax2i*(r2 - (PMEFloat)4.0*cSim.rgbmax*dij - PSATOMS2(j)) +
                 (PMEFloat)2.0*log((dij - sj)*cSim.rgbmax1i));
          }
          else if (dij > (PMEFloat)4.0 * sj) {
            PMEFloat tmpsd = PSATOMS2(j) * dij2i;
            PMEFloat dumbo = ta + tmpsd*(tb + tmpsd*(tc + tmpsd*(td + tmpsd*tdd)));
            dr = tmpsd * sj * dij2i * dumbo;
          }
          else {
            PMEFloat v2 = (PMEFloat)1.0 / (dij + sj);
            if (dij > ri + sj) {
              PMEFloat v4 = log(v2 * (dij - sj));
              dr = (PMEFloat)0.5*(sj/(r2 - PSATOMS2(j)) + (PMEFloat)0.5*dij1i*v4);
            }
            else if (dij > fabs(ri - sj)) {
              PMEFloat v4 = log(v2 * ri);
              PMEFloat theta = (PMEFloat)0.5 * ri1i * dij1i * (r2 + ri*ri - PSATOMS2(j));
              dr = (PMEFloat)0.25 * (ri1i*((PMEFloat)2.0 - theta) - v2 + dij1i*v4);
            }
            else if (ri < sj) {
              PMEFloat v4 = log(v2*(sj - dij));
              dr = (PMEFloat)0.5*(sj/(r2 - PSATOMS2(j)) + (PMEFloat)2.0*ri1i +
                   (PMEFloat)0.5*dij1i*v4);
            }
          }
#ifdef GB_IGB78
          if (dij < cSim.gb_neckcut + ri + PSATOMR(j)) {
            unsigned int ii = round((ri - cSim.gb_neckoffset) * (PMEFloat)20.0);
            unsigned int jj = round((PSATOMR(j) - cSim.gb_neckoffset) * (PMEFloat)20.0);
            PMEFloat2 neckValPos = sNeckMaxValPos[ii * 21 + jj];
            PMEFloat mdist = dij - neckValPos.y;
            PMEFloat mdist2 = mdist * mdist;
            PMEFloat mdist6 = mdist2 * mdist;
            mdist6 = mdist6 * mdist6;
            PMEFloat neck = neckValPos.x / ((PMEFloat)1.0 + mdist2 + (PMEFloat)0.3*mdist6);
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
            dr += gb_neckscale_temp1 * neck;
            // dr += cSim.gb_neckscale * neck;
          }
#endif
#ifdef use_SPFP
#ifdef AMBER_PLATFORM_AMD
          reff_i -= dr;
#else
          reff_i -= fast_llrintf(FORCESCALEF * dr);
#endif
#else
          reff_i -= (PMEForce)dr;
#endif
        }
#ifdef AMBER_PLATFORM_AMD
        shAtom.x = WarpRotateLeft<GRID>(shAtom.x);
        shAtom.y = WarpRotateLeft<GRID>(shAtom.y);
        shAtom.z = WarpRotateLeft<GRID>(shAtom.z);
        shAtom.r = WarpRotateLeft<GRID>(shAtom.r);
        shAtom.s = WarpRotateLeft<GRID>(shAtom.s);
        shAtom.s2 = WarpRotateLeft<GRID>(shAtom.s2);
        shAtom.r1i = WarpRotateLeft<GRID>(shAtom.r1i);
        j = ((j + 1) & (GRID - 1));
#else
        shAtom.x   = __SHFL(mask1, shAtom.x, shIdx);
        shAtom.y   = __SHFL(mask1, shAtom.y, shIdx);
        shAtom.z   = __SHFL(mask1, shAtom.z, shIdx);
        shAtom.r   = __SHFL(mask1, shAtom.r, shIdx);
        shAtom.s   = __SHFL(mask1, shAtom.s, shIdx);
        shAtom.s2  = __SHFL(mask1, shAtom.s2, shIdx);
        shAtom.r1i = __SHFL(mask1, shAtom.r1i, shIdx);
        j = sNext[j];
        mask1 = __BALLOT(mask1, j != tgx);
#endif
      }
      int offset = x + tgx;
#ifdef use_SPFP
#ifdef AMBER_PLATFORM_AMD
      atomicAdd((unsigned long long int*)&cSim.pReffAccumulator[offset],
                llitoulli(fast_llrintf(FORCESCALEF * reff_i)));
#else
      atomicAdd((unsigned long long int*)&cSim.pReffAccumulator[offset], llitoulli(reff_i));
#endif
#else // use_DPFP
      atomicAdd((unsigned long long int*)&cSim.pReffAccumulator[offset],
                llitoulli(llrint(reff_i * FORCESCALE)));
#endif
    }
    else {
      unsigned int j  = y + tgx;
      PMEFloat2 xyj   = cSim.pAtomXYSP[j];
      PSATOMS(tgx)    = cSim.pAtomS[j];
      PSATOMR(tgx)    = cSim.pAtomRBorn[j];
      PSATOMZ(tgx)    = cSim.pAtomZSP[j];
      PMEFloat xi     = xyi.x;
      PMEFloat yi     = xyi.y;
      ri             -= cSim.offset;
      PMEFloat ri1i   = (PMEFloat)1.0 / ri;
      PSATOMX(tgx)    = xyj.x;
      PSATOMY(tgx)    = xyj.y;
      PSATOMS2(tgx)   = PSATOMS(tgx) * PSATOMS(tgx);
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
        PMEFloat r2    = xij * xij + yij * yij + zij * zij;
        PMEFloat dij   = sqrt(r2);
        PMEFloat dij1i = (PMEFloat)1.0 / dij;
        PMEFloat dij2i = dij1i * dij1i;
        PMEFloat sj    = PSATOMS(j);
        if (dij < cSim.rgbmax + sj) {
          PMEFloat dr;
          if (dij > cSim.rgbmax - sj) {
            PMEFloat uij = (PMEFloat)1.0 / (dij - sj);
            dr = (PMEFloat)0.125 * dij1i * ((PMEFloat)1.0 + (PMEFloat)2.0 * dij * uij +
                 cSim.rgbmax2i * (r2 - (PMEFloat)4.0 * cSim.rgbmax * dij - PSATOMS2(j)) +
                 (PMEFloat)2.0 * log((dij - sj) * cSim.rgbmax1i));
          }
          else if (dij > (PMEFloat)4.0 * sj) {
            PMEFloat tmpsd = PSATOMS2(j) * dij2i;
            PMEFloat dumbo = ta + tmpsd*(tb + tmpsd*(tc + tmpsd*(td + tmpsd*tdd)));
            dr = tmpsd * sj * dij2i * dumbo;
          }
          else {
            PMEFloat v2 = (PMEFloat)1.0 / (dij + sj);
            if (dij > ri + sj) {
              PMEFloat v4 = log(v2 * (dij - sj));
              dr = (PMEFloat)0.5 * (sj / (r2 - PSATOMS2(j)) + (PMEFloat)0.5 * dij1i * v4);
            }
            else if (dij > fabs(ri - sj)) {
              PMEFloat v4 = log(v2 * ri);
              PMEFloat theta = (PMEFloat)0.5 * ri1i * dij1i * (r2 + ri*ri - PSATOMS2(j));
              dr = (PMEFloat)0.25 * (ri1i*((PMEFloat)2.0 - theta) - v2 + dij1i*v4);
            }
            else if (ri < sj) {
              PMEFloat v4 = log(v2 * (sj - dij));
              dr = (PMEFloat)0.5*(sj / (r2 - PSATOMS2(j)) + (PMEFloat)2.0*ri1i +
                   (PMEFloat)0.5*dij1i*v4);
            }
          }
#ifdef GB_IGB78
          if (dij < cSim.gb_neckcut + ri + PSATOMR(j)) {
            unsigned int ii = round((ri - cSim.gb_neckoffset)*(PMEFloat)20.0);
            unsigned int jj = round((PSATOMR(j) - cSim.gb_neckoffset)*(PMEFloat)20.0);
            PMEFloat2 neckValPos = sNeckMaxValPos[ii*21 + jj];
            PMEFloat mdist = dij - neckValPos.y;
            PMEFloat mdist2 = mdist * mdist;
            PMEFloat mdist6 = mdist2 * mdist;
            mdist6 = mdist6 * mdist6;
            PMEFloat neck = neckValPos.x / ((PMEFloat)1.0 + mdist2 + (PMEFloat)0.3 * mdist6);
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
            dr += gb_neckscale_temp1 * neck;
            // dr += cSim.gb_neckscale * neck;
          }
#endif
#ifdef use_SPFP
#ifdef AMBER_PLATFORM_AMD
          reff_i -= dr;
#else
          reff_i -= fast_llrintf(FORCESCALEF * dr);
#endif
#else // use_DPFP
          reff_i -= (PMEForce)dr;
#endif
        }
        if (dij < cSim.rgbmax + si) {
          PMEFloat dr;
          if (dij > cSim.rgbmax - si) {
            PMEFloat uij = (PMEFloat)1.0 / (dij - si);
            dr = (PMEFloat)0.125 * dij1i * ((PMEFloat)1.0 + (PMEFloat)2.0 * dij * uij +
                 cSim.rgbmax2i * (r2 - (PMEFloat)4.0 * cSim.rgbmax * dij - si2) +
                 (PMEFloat)2.0 * log((dij - si) * cSim.rgbmax1i));
          }
          else if (dij > (PMEFloat)4.0 * si) {
            PMEFloat tmpsd = si2 * dij2i;
            PMEFloat dumbo = ta + tmpsd*(tb + tmpsd*(tc + tmpsd*(td + tmpsd*tdd)));
            dr = tmpsd * si * dij2i * dumbo;
          }
          else {
            PMEFloat v2 = (PMEFloat)1.0 / (dij + si);
            if (dij > PSATOMR(j) + si) {
              PMEFloat v4 = log(v2*(dij - si));
              dr = (PMEFloat)0.5 * (si / (r2 - si2) + (PMEFloat)0.5 * dij1i * v4);
            }
            else if (dij > fabs(PSATOMR(j) - si)) {
              PMEFloat v4 = log(v2 * PSATOMR(j));
              PMEFloat theta = (PMEFloat)0.5 * PSATOMR1I(j) * dij1i *
                               (r2 + PSATOMR(j)*PSATOMR(j) - si2);
              dr = (PMEFloat)0.25*(PSATOMR1I(j)*((PMEFloat)2.0 - theta) - v2 + dij1i*v4);
            }
            else if (PSATOMR(j) < si) {
              PMEFloat v4 = log(v2 * (si - dij));
              dr = (PMEFloat)0.5 * (si / (r2 - si2) + (PMEFloat)2.0 * PSATOMR1I(j) +
                   (PMEFloat)0.5 * dij1i * v4);
            }
          }
#ifdef GB_IGB78
          if (dij < cSim.gb_neckcut + ri + PSATOMR(j)) {
            unsigned int ii = round((ri - cSim.gb_neckoffset) * (PMEFloat)20.0);
            unsigned int jj = round((PSATOMR(j) - cSim.gb_neckoffset) * (PMEFloat)20.0);
            PMEFloat2 neckValPos = sNeckMaxValPos[jj * 21 + ii];
            PMEFloat mdist  = dij - neckValPos.y;
            PMEFloat mdist2 = mdist * mdist;
            PMEFloat mdist6 = mdist2 * mdist;
            mdist6          = mdist6 * mdist6;
            PMEFloat neck   = neckValPos.x / ((PMEFloat)1.0 + mdist2 + (PMEFloat)0.3 * mdist6);
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
            dr += gb_neckscale_temp1 * neck;
            // dr += cSim.gb_neckscale * neck;
          }
#endif
#ifdef use_SPFP
#ifdef AMBER_PLATFORM_AMD
          PSREFF(j) -= dr;
#else
          PSREFF(j) -= fast_llrintf(FORCESCALEF * dr);
#endif
#else  // use_DPFP
          PSREFF(j) -= (PMEForce)dr;
#endif // End case switch for different precision modes
        }

#if defined(use_SPFP) && defined(AMBER_PLATFORM_AMD)
        PSREFF(j)  = __SHFL(mask1, PSREFF(j), shIdx);
#endif
#ifdef AMBER_PLATFORM_AMD
        shAtom.x = WarpRotateLeft<GRID>(shAtom.x);
        shAtom.y = WarpRotateLeft<GRID>(shAtom.y);
        shAtom.z = WarpRotateLeft<GRID>(shAtom.z);
        shAtom.r = WarpRotateLeft<GRID>(shAtom.r);
        shAtom.s = WarpRotateLeft<GRID>(shAtom.s);
        shAtom.s2 = WarpRotateLeft<GRID>(shAtom.s2);
        shAtom.r1i = WarpRotateLeft<GRID>(shAtom.r1i);
        j = ((j + 1) & (GRID - 1));
      }
#else
        shAtom.x   = __SHFL(mask1, shAtom.x, shIdx);
        shAtom.y   = __SHFL(mask1, shAtom.y, shIdx);
        shAtom.z   = __SHFL(mask1, shAtom.z, shIdx);
        shAtom.r   = __SHFL(mask1, shAtom.r, shIdx);
        shAtom.s   = __SHFL(mask1, shAtom.s, shIdx);
        shAtom.s2  = __SHFL(mask1, shAtom.s2, shIdx);
        shAtom.r1i = __SHFL(mask1, shAtom.r1i, shIdx);
        j = sNext[j];
        mask1 = __BALLOT(mask1, j != tgx);
      } while (j != tgx);
#endif
      // End do ... while loop beginning roughly 130 lines ago

#ifdef use_SPFP
      int offset = x + tgx;
#ifdef AMBER_PLATFORM_AMD
      atomicAdd((unsigned long long int*)&cSim.pReffAccumulator[offset],
                llitoulli(fast_llrintf(FORCESCALEF * reff_i)));
      offset = y + tgx;
      atomicAdd((unsigned long long int*)&cSim.pReffAccumulator[offset],
                llitoulli(fast_llrintf(FORCESCALEF * PSREFF(tgx))));
#else
      atomicAdd((unsigned long long int*)&cSim.pReffAccumulator[offset], llitoulli(reff_i));
      offset = y + tgx;
      atomicAdd((unsigned long long int*)&cSim.pReffAccumulator[offset],
                llitoulli(PSREFF(tgx)));
#endif
#else  // use_DPFP
      int offset = x + tgx;
      atomicAdd((unsigned long long int*)&cSim.pReffAccumulator[offset],
                llitoulli(llrint(reff_i * FORCESCALE)));
      offset = y + tgx;
      atomicAdd((unsigned long long int*)&cSim.pReffAccumulator[offset],
                llitoulli(llrint(PSREFF(tgx) * FORCESCALE)));
#endif // End case switch for different precision modes
    }
#ifndef AMBER_PLATFORM_AMD
    if (tgx == 0) {
      *psPos = atomicAdd(cSim.pGBBRPosition, 1);
    }
#endif
  }
}
