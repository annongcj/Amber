#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------
{

#ifdef PME_ENERGY
  __shared__ PMEDouble sE[GENERAL_THREADS_PER_BLOCK];
#endif

#ifdef PME_VIRIAL
  struct ScalarSumVirial {
    PMEForce vir_11;
    PMEForce vir_22;
    PMEForce vir_33;
  };
  __shared__ ScalarSumVirial sV[GENERAL_THREADS_PER_BLOCK];
#endif

#if defined(MBAR) && defined(PME_ENERGY)
  PMEAccumulator bar_cont[16];
  for (int i = 0; i < 16; i++) {
    bar_cont[i] = 0;
  }
#endif
#if defined(use_DPFP)
#  define PREFACSIZE 2030
#else
#  define PREFACSIZE 4000
#endif
  __shared__ PMEFloat sPrefac[PREFACSIZE];
#ifdef PME_VIRIAL
  __shared__ PMEFloat sRecipf[9];
#endif
  int increment = blockDim.x * gridDim.x;

  // Insure q[0][0][0] == 0
  if ((blockIdx.x == 0) && (threadIdx.x == 0)) {
    PMEComplex cmplx = {(PMEFloat)0.0, (PMEFloat)0.0};
    cSim.pXYZ_qt[0]  = cmplx;
  }

#ifdef PME_ENERGY
  // Initialize energy parameters
  sE[threadIdx.x] = (PMEDouble)0.0;
#endif

#ifdef PME_VIRIAL
  // Initialize energy parameters
  PMEForce vir_11 = (PMEForce)0;
  PMEForce vir_22 = (PMEForce)0;
  PMEForce vir_33 = (PMEForce)0;
  if (threadIdx.x < 9) {
    sRecipf[threadIdx.x] = cSim.pNTPData->recipf[threadIdx.x];
  }
  __syncthreads();
#endif
  if (cSim.nSum <= PREFACSIZE) {

    // Load prefac data
    unsigned int pos = threadIdx.x;
    while (pos < cSim.nSum) {
      sPrefac[pos]  = cSim.pPrefac1[pos];
      pos          += blockDim.x;
    }
    __syncthreads();
    PMEFloat* psPrefac2 = &sPrefac[cSim.n2Offset];
    PMEFloat* psPrefac3 = &sPrefac[cSim.n3Offset];
    pos = blockIdx.x * blockDim.x + threadIdx.x;
    while (pos < cSim.fft_x_y_z_quarter_dim) {

      // Calculate position
      int k3 = pos / cSim.fft_y_dim_times_x_dim;
      int offset = pos - (k3 * cSim.fft_y_dim_times_x_dim);
      int k2 = offset / cSim.fft_x_dim;
      int k1 = offset - (k2 * cSim.fft_x_dim);

      // m{1,2,3} will now be the position on the complex grid, while k{1,2,3} get incremented
      // to access arrays imported from the Fortran code that retain the Fortran indexing.
      int m1 = k1;
      int m2 = k2;
      k1++;
      k2++;
      m1 -= (k1 > cSim.nf1) * cSim.nfft1;
      m2 -= (k2 > cSim.nf2) * cSim.nfft2;
      PMEFloat fm1 = m1;
      PMEFloat fm2 = m2;
#ifdef PME_VIRIAL
      PMEFloat mhat1 = fm1*sRecipf[0];
      PMEFloat mhat2 = fm1*sRecipf[3] + fm2*sRecipf[4];
      PMEFloat mhat3_base = fm1*sRecipf[6] + fm2*sRecipf[7];
      PMEFloat sP_base    = sPrefac[k1] * psPrefac2[k2] * pi_vol_inv;
#else
      PMEFloat mhat1 = fm1*cSim.recipf[0][0];
      PMEFloat mhat2 = fm1*cSim.recipf[1][0] + fm2*cSim.recipf[1][1];
      PMEFloat mhat3_base = fm1*cSim.recipf[2][0] + fm2*cSim.recipf[2][1];
      PMEFloat sP_base    = sPrefac[k1] * psPrefac2[k2] * cSim.pi_vol_inv;
#endif
      PMEFloat msq_base   = mhat1*mhat1 + mhat2*mhat2;

      // Loop four times to conserve pre-calculated values.  Increment k3 if the first
      // iteration is oging to get skipped (i.e. when k1 == k2 == 1).
      int i;
      int imin = (pos == 0);
      k3 += imin * (cSim.fft_z_dim >> 2);
      for (i = imin; i < 4; i++) {
        int posiz = pos + i*cSim.fft_x_y_z_quarter_dim;
        PMEComplex q = cSim.pXYZ_qt[posiz];
        int m3 = k3;
        k3++;
        m3 -= (k3 > cSim.nf3) * cSim.nfft3;
        PMEFloat fm3 = m3;
#ifdef PME_VIRIAL
        PMEFloat mhat3 = mhat3_base + fm3*sRecipf[8];
#else
        PMEFloat mhat3 = mhat3_base + fm3*cSim.recipf[2][2];
#endif
        PMEFloat msq = msq_base + mhat3*mhat3;
        PMEFloat msq_inv = (PMEFloat)1.0 / msq;
        PMEFloat eterm = exp(-cSim.fac * msq) * sP_base * psPrefac3[k3] * msq_inv;
#ifdef PME_VIRIAL
        PMEFloat vterm = cSim.fac2 + 2.0 * msq_inv;
#endif
#if defined(PME_ENERGY) || defined (PME_VIRIAL)
        PMEFloat eterm_struc2 = eterm * (q.x*q.x + q.y*q.y);
#endif
#ifdef PME_ENERGY
	sE[threadIdx.x] += (PMEDouble)eterm_struc2;
#endif
#ifdef PME_VIRIAL
#  ifndef use_DPFP
        eterm_struc2 *= FORCESCALEF;
        PMEFloat v11 = eterm_struc2 * (vterm*mhat1*mhat1 - (PMEFloat)1.0);
        PMEFloat v22 = eterm_struc2 * (vterm*mhat2*mhat2 - (PMEFloat)1.0);
        PMEFloat v33 = eterm_struc2 * (vterm*mhat3*mhat3 - (PMEFloat)1.0);
#  else
        vir_11 += eterm_struc2 * (vterm*mhat1*mhat1 - (PMEDouble)1.0);
        vir_22 += eterm_struc2 * (vterm*mhat2*mhat2 - (PMEDouble)1.0);
        vir_33 += eterm_struc2 * (vterm*mhat3*mhat3 - (PMEDouble)1.0);
#  endif
#endif
#if defined(PME_ENERGY) || defined(PME_VIRIAL)
        if ((k1 > 1) && (k1 <= cSim.nfft1)) {
          int k1s, k2s, k3s;
          int m1s, m2s, m3s;

          k1s = cSim.nfft1 - k1 + 2;
          k2s = ((cSim.nfft2 - k2 + 1) % cSim.nfft2) + 1;
          k3s = ((cSim.nfft3 - k3 + 1) % cSim.nfft3) + 1;
          m1s = k1s - 1;
          m2s = k2s - 1;
          m3s = k3s - 1;
          m1s -= (k1s > cSim.nf1) * cSim.nfft1;
          m2s -= (k2s > cSim.nf2) * cSim.nfft2;
          m3s -= (k3s > cSim.nf3) * cSim.nfft3;
          PMEFloat fm1s = m1s;
          PMEFloat fm2s = m2s;
          PMEFloat fm3s = m3s;
#  ifdef PME_VIRIAL
          PMEFloat mhat1s = fm1s*sRecipf[0];
          PMEFloat mhat2s = fm1s*sRecipf[3] + fm2s*sRecipf[4];
          PMEFloat mhat3s = fm1s*sRecipf[6] + fm2s*sRecipf[7] + fm3s*sRecipf[8];
#  else
          PMEFloat mhat1s = fm1s*cSim.recipf[0][0];
          PMEFloat mhat2s = fm1s*cSim.recipf[1][0] + fm2s*cSim.recipf[1][1];
          PMEFloat mhat3s = fm1s*cSim.recipf[2][0] + fm2s*cSim.recipf[2][1] +
                            fm3s*cSim.recipf[2][2];
#  endif
          PMEFloat msqs = mhat1s*mhat1s + mhat2s*mhat2s + mhat3s*mhat3s;
          PMEFloat msqs_inv = (PMEFloat)1.0 / msqs;
#  ifdef PME_VIRIAL
          PMEFloat eterms = exp(-cSim.fac * msqs) * sPrefac[k1s] * psPrefac2[k2s] *
                            psPrefac3[k3s] * pi_vol_inv * msqs_inv;
          PMEFloat vterms = cSim.fac2 + (PMEFloat)2.0 * msqs_inv;
#  else
          PMEFloat eterms = exp(-cSim.fac * msqs) * sPrefac[k1s] * psPrefac2[k2s] *
                            psPrefac3[k3s] * cSim.pi_vol_inv * msqs_inv;
#  endif
          PMEFloat eterms_struc2s = eterms * (q.x*q.x + q.y*q.y);
#  ifdef PME_ENERGY
          sE[threadIdx.x] += eterms_struc2s;
#  endif
#  ifdef PME_VIRIAL
#    ifndef use_DPFP
          eterms_struc2s *= FORCESCALEF;
          v11 += eterms_struc2s * (vterms*mhat1*mhat1 - (PMEFloat)1.0);
          v22 += eterms_struc2s * (vterms*mhat2*mhat2 - (PMEFloat)1.0);
          v33 += eterms_struc2s * (vterms*mhat3*mhat3 - (PMEFloat)1.0);
#    else
          vir_11 += eterms_struc2s * (vterms*mhat1s*mhat1s - (PMEDouble)1.0);
          vir_22 += eterms_struc2s * (vterms*mhat2s*mhat2s - (PMEDouble)1.0);
          vir_33 += eterms_struc2s * (vterms*mhat3s*mhat3s - (PMEDouble)1.0);
#    endif
#  endif
        }
#  if defined(PME_VIRIAL) && !defined(use_DPFP)
        vir_11 += fast_llrintf(v11);
        vir_22 += fast_llrintf(v22);
        vir_33 += fast_llrintf(v33);
#  endif
#endif // defined(PME_VIRIAL) || defined(PME_ENERGY)
        q.x *= eterm;
        q.y *= eterm;
        cSim.pXYZ_qt[posiz] = q;

        // Increment
        k3 += cSim.fft_quarter_z_dim_m1;
      }

      // Increment
      pos += increment;
    }
  }
  else {

    // Rework the incrementation--this uses the old scheme for iterating through the grid.
    // But, there should hardly ever be a case where the grid is so large that pre-factors
    // cannot be cached in __shared__.
    int zIncrement  = increment / (cSim.fft_y_dim_times_x_dim);
    increment      -= zIncrement * cSim.fft_y_dim_times_x_dim;
    int yIncrement  = increment / cSim.fft_x_dim;
    int xIncrement  = increment - yIncrement * cSim.fft_x_dim;
    xIncrement--;
    yIncrement--;
    zIncrement--;

    // Calculate initial position
    int offset  = blockIdx.x*blockDim.x + threadIdx.x + 1;
    int k3      = offset / (cSim.fft_y_dim_times_x_dim);
    offset     -= k3 * cSim.fft_y_dim_times_x_dim;
    int k2      = offset / cSim.fft_x_dim;
    int k1      = offset - k2*cSim.fft_x_dim;
    int m1, m2, m3;

    while (k3 < cSim.fft_z_dim) {

      // Read data
      unsigned int pos = (k3*cSim.fft_y_dim + k2) * cSim.fft_x_dim + k1;
      PMEComplex q = cSim.pXYZ_qt[pos];
      k1++;
      k2++;
      k3++;

      // Generate internal coordinates
      m1 = k1 - 1;
      m2 = k2 - 1;
      m3 = k3 - 1;
      if (k1 > cSim.nf1) {
        m1 -= cSim.nfft1;
      }
      if (k2 > cSim.nf2) {
        m2 -= cSim.nfft2;
      }
      if (k3 > cSim.nf3) {
        m3 -= cSim.nfft3;
      }
#ifdef PME_VIRIAL
      PMEFloat mhat1 = m1*sRecipf[0];
      PMEFloat mhat2 = m1*sRecipf[3] + m2*sRecipf[4];
      PMEFloat mhat3 = m1*sRecipf[6] + m2*sRecipf[7] + m3*sRecipf[8];
#else
      PMEFloat mhat1 = m1*cSim.recipf[0][0];
      PMEFloat mhat2 = m1*cSim.recipf[1][0] + m2*cSim.recipf[1][1];
      PMEFloat mhat3 = m1*cSim.recipf[2][0] + m2*cSim.recipf[2][1] + m3*cSim.recipf[2][2];
#endif
      PMEFloat msq     = mhat1*mhat1 + mhat2*mhat2 + mhat3*mhat3;
      PMEFloat msq_inv = (PMEFloat)1.0 / msq;
#ifdef PME_VIRIAL
      PMEFloat eterm = exp(-cSim.fac * msq) * cSim.pPrefac1[k1] * cSim.pPrefac2[k2] *
                       cSim.pPrefac3[k3] * pi_vol_inv * msq_inv;
      PMEFloat vterm = cSim.fac2 + (PMEFloat)2.0 * msq_inv;
#else
      PMEFloat eterm = exp(-cSim.fac * msq) * cSim.pPrefac1[k1] * cSim.pPrefac2[k2] *
                       cSim.pPrefac3[k3] * cSim.pi_vol_inv * msq_inv;
#endif
#if defined(PME_VIRIAL) || defined(PME_ENERGY)
      PMEFloat eterm_struc2 = eterm * (q.x*q.x + q.y*q.y);
#endif
#ifdef PME_ENERGY
      sE[threadIdx.x] += (PMEDouble)eterm_struc2;
#endif
#ifdef PME_VIRIAL
#  ifndef use_DPFP
      eterm_struc2 *= FORCESCALEF;
      vir_11 += fast_llrintf(eterm_struc2 * (vterm*mhat1*mhat1 - (PMEFloat)1.0));
      vir_22 += fast_llrintf(eterm_struc2 * (vterm*mhat2*mhat2 - (PMEFloat)1.0));
      vir_33 += fast_llrintf(eterm_struc2 * (vterm*mhat3*mhat3 - (PMEFloat)1.0));
#  else
      vir_11 += eterm_struc2 * (vterm*mhat1*mhat1 - (PMEDouble)1.0);
      vir_22 += eterm_struc2 * (vterm*mhat2*mhat2 - (PMEDouble)1.0);
      vir_33 += eterm_struc2 * (vterm*mhat3*mhat3 - (PMEDouble)1.0);
#  endif
#endif

#if defined(PME_VIRIAL) || defined(PME_ENERGY)
      if ((k1 > 1) && (k1 <= cSim.nfft1)) {
        int k1s, k2s, k3s;
        int m1s, m2s, m3s;
        k1s = cSim.nfft1 - k1 + 2;
        k2s = ((cSim.nfft2 - k2 + 1) % cSim.nfft2) + 1;
        k3s = ((cSim.nfft3 - k3 + 1) % cSim.nfft3) + 1;
        m1s = k1s - 1;
        m2s = k2s - 1;
        m3s = k3s - 1;
        if (k1s > cSim.nf1) {
          m1s -= cSim.nfft1;
        }
        if (k2s > cSim.nf2) {
          m2s -= cSim.nfft2;
        }
        if (k3s > cSim.nf3) {
          m3s -= cSim.nfft3;
        }
#  ifdef PME_VIRIAL
        PMEFloat mhat1s = m1s*sRecipf[0];
        PMEFloat mhat2s = m1s*sRecipf[3] + m2s*sRecipf[4];
        PMEFloat mhat3s = m1s*sRecipf[6] + m2s*sRecipf[7] + m3s*sRecipf[8];
#  else
        PMEFloat mhat1s = m1s*cSim.recipf[0][0];
        PMEFloat mhat2s = m1s*cSim.recipf[1][0] + m2s*cSim.recipf[1][1];
        PMEFloat mhat3s = m1s*cSim.recipf[2][0] + m2s*cSim.recipf[2][1] +
                          m3s*cSim.recipf[2][2];
#  endif
        PMEFloat msqs = mhat1s*mhat1s + mhat2s*mhat2s + mhat3s*mhat3s;
        PMEFloat msqs_inv = (PMEFloat)1.0 / msqs;
#  ifdef PME_VIRIAL
        PMEFloat eterms = exp(-cSim.fac * msqs) * cSim.pPrefac1[k1s] * cSim.pPrefac2[k2s] *
                          cSim.pPrefac3[k3s] * pi_vol_inv * msqs_inv;
        PMEFloat vterms = cSim.fac2 + (PMEFloat)2.0 * msqs_inv;
#  else
        PMEFloat eterms = exp(-cSim.fac * msqs) * cSim.pPrefac1[k1s] * cSim.pPrefac2[k2s] *
                          cSim.pPrefac3[k3s] * cSim.pi_vol_inv * msqs_inv;
#  endif
        PMEFloat eterms_struc2s = eterms * (q.x * q.x + q.y * q.y);
#  ifdef PME_ENERGY
        sE[threadIdx.x] += (PMEDouble)eterms_struc2s;
#  endif
#  ifdef PME_VIRIAL
#    ifndef use_DPFP
        eterms_struc2s *= FORCESCALEF;
        vir_11 += fast_llrintf(eterms_struc2s * (vterms*mhat1*mhat1 - (PMEFloat)1.0));
        vir_22 += fast_llrintf(eterms_struc2s * (vterms*mhat2*mhat2 - (PMEFloat)1.0));
        vir_33 += fast_llrintf(eterms_struc2s * (vterms*mhat3*mhat3 - (PMEFloat)1.0));
#    else
        vir_11 += eterms_struc2s * (vterms*mhat1s*mhat1s - (PMEDouble)1.0);
        vir_22 += eterms_struc2s * (vterms*mhat2s*mhat2s - (PMEDouble)1.0);
        vir_33 += eterms_struc2s * (vterms*mhat3s*mhat3s - (PMEDouble)1.0);
#    endif
#  endif
      }
#endif // defined(PME_VIRIAL) || defined(PME_ENERGY)
      q.x *= eterm;
      q.y *= eterm;
      cSim.pXYZ_qt[pos] = q;

      // Increment position
      k1 += xIncrement;
      if (k1 >= cSim.fft_x_dim) {
        k1 -= cSim.fft_x_dim;
        k2++;
      }
      k2 += yIncrement;
      if (k2 >= cSim.fft_y_dim) {
        k2 -= cSim.fft_y_dim;
        k3++;
      }
      k3 += zIncrement;
    }
  }

#if defined(PME_ENERGY) || defined(PME_VIRIAL)
#  ifdef PME_VIRIAL
  sV[threadIdx.x].vir_11 = vir_11;
  sV[threadIdx.x].vir_22 = vir_22;
  sV[threadIdx.x].vir_33 = vir_33;
#  endif

  // Reduce virial and energy
  __syncthreads();
  int m = 1;
  while (m < blockDim.x) {
    int p = threadIdx.x + m;
#  ifdef PME_ENERGY
    PMEDouble energy = ((p < blockDim.x) ? sE[p] : (PMEDouble)0.0);
#  endif
#  ifdef PME_VIRIAL
    PMEDouble vir_11 = ((p < blockDim.x) ? sV[p].vir_11 : (PMEDouble)0.0);
    PMEDouble vir_22 = ((p < blockDim.x) ? sV[p].vir_22 : (PMEDouble)0.0);
    PMEDouble vir_33 = ((p < blockDim.x) ? sV[p].vir_33 : (PMEDouble)0.0);
#  endif
    __syncthreads();
#  ifdef PME_ENERGY
    sE[threadIdx.x] += energy;
#  endif
#  ifdef PME_VIRIAL
    sV[threadIdx.x].vir_11 += vir_11;
    sV[threadIdx.x].vir_22 += vir_22;
    sV[threadIdx.x].vir_33 += vir_33;
#  endif
    __syncthreads();
    m *= 2;
  }
#  ifdef PME_ENERGY
#    if !defined(SC_REGION_1) && !defined(SC_REGION_2)
  unsigned long long int eer = llitoulli(llrint((PMEDouble)0.5 * sE[threadIdx.x] *
                                                ENERGYSCALE));
#    endif
#    ifdef SC_REGION_1
  // For energies, we don't want signed, we want a full potential.
  // For dvdl, we want the full unscaled potential, multiplied by the ti_sign.
  unsigned long long int eer = llitoulli(llrint((PMEDouble)0.5 * sE[threadIdx.x] *
                                                cSim.AFElambda[0] * ENERGYSCALE));
  unsigned long long int dvdl = llitoulli(llrint((PMEDouble)0.5 * sE[threadIdx.x] *
                                                 cSim.TIsigns[0] * ENERGYSCALE));
#      ifdef MBAR
  for (int i = 0; i < cSim.bar_states; i++) {
    bar_cont[i] = llrint(sE[threadIdx.x] * (PMEDouble)0.5 * ENERGYSCALE *
                         (cSim.pBarLambda[i] - cSim.AFElambda[0]));
  }
#      endif // MBAR
#    endif
#    ifdef SC_REGION_2
  // Don't need a sign for region 2
  unsigned long long int eer = llitoulli(llrint((PMEDouble)0.5 * sE[threadIdx.x] *
                                                cSim.AFElambda[1] * ENERGYSCALE));
  unsigned long long int dvdl = llitoulli(llrint((PMEDouble)0.5 * sE[threadIdx.x] *
                                                 ENERGYSCALE));
#      ifdef MBAR
  for (int i = 0; i < cSim.bar_states; i++) {
    bar_cont[i] = llrint(sE[threadIdx.x] * (PMEDouble)0.5 * ENERGYSCALE *
                         (cSim.pBarLambda[cSim.bar_stride + i] - cSim.AFElambda[1]));
  }
#      endif // MBAR
#    endif
#  endif
#  ifdef PME_VIRIAL
#    ifndef use_DPFP
  unsigned long long int vir11 = llitoulli(sV[threadIdx.x].vir_11 >> 1);
  unsigned long long int vir22 = llitoulli(sV[threadIdx.x].vir_22 >> 1);
  unsigned long long int vir33 = llitoulli(sV[threadIdx.x].vir_33 >> 1);
#    else
  unsigned long long int vir11 = llitoulli(llrint((PMEDouble)0.5 * sV[threadIdx.x].vir_11 *
                                                  FORCESCALE));
  unsigned long long int vir22 = llitoulli(llrint((PMEDouble)0.5 * sV[threadIdx.x].vir_22 *
                                                  FORCESCALE));
  unsigned long long int vir33 = llitoulli(llrint((PMEDouble)0.5 * sV[threadIdx.x].vir_33 *
                                                  FORCESCALE));
#    endif
#  endif

  // Write out energies
  if (threadIdx.x == 0) {
#  ifdef PME_ENERGY
    atomicAdd(cSim.pEER, eer);
#    if defined(SC_REGION_1) || defined(SC_REGION_2)
    atomicAdd(cSim.pDVDL, dvdl);
#      ifdef MBAR
      for(int i=0; i<cSim.bar_states;i++)
      {
        atomicAdd(&cSim.pBarTot[i], llitoulli(bar_cont[i]));
      }
#      endif
#    endif
#  endif
#  ifdef PME_VIRIAL
    atomicAdd(cSim.pVirial_11, vir11);
    atomicAdd(cSim.pVirial_22, vir22);
    atomicAdd(cSim.pVirial_33, vir33);
#  endif
  }
#endif // defined(PME_ENERGY) || defined(PME_VIRIAL)
#undef PREFACSIZE
}
