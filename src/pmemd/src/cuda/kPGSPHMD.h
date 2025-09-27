{
  const PMEFloat ONE = 1.0;
  const PMEFloat TWO = 2.0;
  const PMEFloat THREE = 3.0;
  const PMEFloat ONEHALF = 0.5;
  const PMEFloat ONETHIRD = (1.0 / 3.0);
  unsigned int titrindex = (blockIdx.x * blockDim.x + threadIdx.x) * cSim.phmd_stride;
  if (titrindex < cSim.ntitratoms) {
    for( int i = 0; i < cSim.phmd_stride; i++ ) {
      unsigned int atom = cSim.ptitratoms[titrindex + i];
      int h = cSim.pgrplist[atom] - 1;
      PMEFloat qunprot, qprot;
      PMEFloat qh = 0;
      PMEFloat qxh = 0;
      int psp_grp = cSim.psp_grp[h];
      PMEFloat pph_theta = cSim.pph_theta[h];
      PMEFloat pluspph_theta = cSim.pph_theta[h+1];
      PMEFloat lambda = sin(pph_theta) * sin(pph_theta);
      PMEFloat x2 = 1.0;
      int ImageIndex = cSim.pImageAtomLookup[atom];
      PMEFloat2 pqstate1 = cSim.pImageQstate1[ImageIndex];
      PMEFloat2 pqstate2 = cSim.pImageQstate2[ImageIndex];
      if (psp_grp > 0) {
        x2 = sin(pluspph_theta) * sin(pluspph_theta);
        qunprot = lambda * (pqstate2.x - pqstate2.y);
        qprot = (1 - lambda) * (pqstate1.x - pqstate1.y);
        qxh = qunprot + qprot;
      }
      qunprot = x2 * pqstate2.x + (1 - x2) * pqstate2.y;
      qprot = x2 * pqstate1.x + (1 - x2) * pqstate1.y;
      qh = qunprot - qprot;
      PMEFloat fx = cSim.pFractX[ImageIndex];
      PMEFloat fy = cSim.pFractY[ImageIndex];
      PMEFloat fz = cSim.pFractZ[ImageIndex];
      int ix = int(fx);
      int iy = int(fy);
      int iz = int(fz);
      fx -= ix;
      fy -= iy;
      fz -= iz;
      ix -= cSim.orderMinusOne;
      iy -= cSim.orderMinusOne;
      iz -= cSim.orderMinusOne;
      if (ix < 0) {
        ix += cSim.nfft1;
      }
      if (iy < 0) {
        iy += cSim.nfft2;
      }
      if (iz < 0) {
        iz += cSim.nfft3;
      }
      // Order 2 B-spline accumulation
      PMEFloat tx[4];
      tx[0] = ONE - fx;
      tx[1] = fx;
      // Order 3 B-spline accumulation
      tx[2] = ONEHALF * fx * tx[1];
      tx[1] = ONEHALF * ((fx + ONE) * tx[0] + (TWO - fx) * tx[1]);
      tx[0] = ONEHALF * (ONE - fx) * tx[0];
      // Order 4 B-spline accumulation
      tx[3] = ONETHIRD * fx * tx[2];
      tx[2] = ONETHIRD * ((fx + ONE) * tx[1] + (THREE - fx) * tx[2]);
      tx[1] = ONETHIRD *((fx + TWO) * tx[0] + (TWO - fx) * tx[1]);
      tx[0] = ONETHIRD * (ONE - fx) * tx[0];
      // Order 2 B-spline accumulation
      PMEFloat ty[4];
      ty[0] = ONE - fy;
      ty[1] = fy;
      // Order 3 B-spline accumulation
      ty[2] = ONEHALF * fy * ty[1];
      ty[1] = ONEHALF * ((fy + ONE) * ty[0] + (TWO - fy) * ty[1]);
      ty[0] = ONEHALF * (ONE - fy) * ty[0];
      // Order 4 B-spline accumulation
      ty[3] = ONETHIRD * fy * ty[2];
      ty[2] = ONETHIRD * ((fy + ONE) * ty[1] + (THREE - fy) * ty[2]);
      ty[1] = ONETHIRD * ((fy + TWO) * ty[0] + (TWO - fy) * ty[1]);
      ty[0] = ONETHIRD * (ONE - fy) * ty[0];
      //Order 2 B-spline accumulation
      PMEFloat tz[4];
      tz[0] = ONE - fz;
      tz[1] = fz;
      // Order 3
      tz[2] = ONEHALF * fz * tz[1];
      tz[1] = ONEHALF * ((fz + ONE) * tz[0] + (TWO - fz) * tz[1]);
      tz[0] = ONEHALF * (ONE - fz) * tz[0];
      // Order 4
      tz[3] = ONETHIRD * fz * tz[2];
      tz[2] = ONETHIRD * ((fz + ONE) * tz[1] + (THREE - fz) * tz[2]);
      tz[1] = ONETHIRD * ((fz + TWO) * tz[0] + (TWO - fz) * tz[1]);
      tz[0] = ONETHIRD * (ONE - fz) * tz[0];
      PMEFloat dudl = 0;
      PMEFloat dudlplus = 0;
      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          for (int k = 0; k < 4; k++) {
            int ix_local = ix + i;
            int iy_local = iy + j;
            int iz_local = iz + k;
            ix_local -= (ix_local >= cSim.nfft1) * cSim.nfft1;
            iy_local -= (iy_local >= cSim.nfft2) * cSim.nfft2;
            iz_local -= (iz_local >= cSim.nfft3) * cSim.nfft3;
#ifdef use_DPFP
            int2 iterm = tex1Dfetch<int2>(cSim.texXYZ_q, ((iz_local * cSim.nfft2 + iy_local) * cSim.nfft1 +
                                                           ix_local));
            PMEFloat qterm = __hiloint2double(iterm.y, iterm.x);
#else
            PMEFloat qterm = tex1Dfetch<float>(cSim.texXYZ_q, ((iz_local * cSim.nfft2 + iy_local) * cSim.nfft1 +
                                                                ix_local));
#endif
            dudl += tx[i] * ty[j] * tz[k] * qterm * qh;
            dudlplus += tx[i] * ty[j] * tz[k] * qterm * qxh;
          }
        }
      }
#ifdef use_DPFP
      atomicAdd((unsigned long long int*)&cSim.pdph_accumulator[ImageIndex],
                 llitoulli(llrint((PMEDouble)dudl * FORCESCALE)));
      if (psp_grp > 0) {
        atomicAdd((unsigned long long int*)&cSim.pdph_plus_accumulator[ImageIndex],
                  llitoulli(llrint((PMEDouble)dudlplus * FORCESCALE)));
      }
#else
      atomicAdd((unsigned long long int*)&cSim.pdph_accumulator[ImageIndex],
                llitoulli(fast_llrintf(FORCESCALEF * dudl)));
      if (psp_grp > 0) {
        atomicAdd((unsigned long long int*)&cSim.pdph_plus_accumulator[ImageIndex],
                  llitoulli(fast_llrintf(FORCESCALEF * dudlplus)));
      }
#endif
    }
  }
}
