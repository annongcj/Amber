#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------
{
#ifdef SHAKE_NEIGHBORLIST
#  define PATOMX(i) cSim.pImageX[i]
#  define PATOMY(i) cSim.pImageY[i]
#  define PATOMZ(i) cSim.pImageZ[i]
#else
#  define PATOMX(i) cSim.pAtomX[i]
#  define PATOMY(i) cSim.pAtomY[i]
#  define PATOMZ(i) cSim.pAtomZ[i]
#endif

#ifdef SHAKE_HMR
#  define INVMASSH invMassH
#else
#  define INVMASSH cSim.invMassH
#endif

  unsigned int pos                            = blockIdx.x * blockDim.x + threadIdx.x;

  if (pos < cSim.shakeOffset) {
    if (pos < cSim.shakeConstraints) {

      // Read SHAKE network data
#ifdef SHAKE_NEIGHBORLIST
      int4 shakeID = cSim.pImageShakeID[pos];
#ifdef TISHAKE2
      int4 refShakeID = cSim.pImageShakeID[pos];
#endif
#else
      int4 shakeID = cSim.pShakeID[pos];
#ifdef TISHAKE2
      int4 refShakeID = cSim.pShakeID[pos];
#endif
#endif
      double2 shakeParm = cSim.pShakeParm[pos];
#ifdef TISHAKE2
      double2 shakeParm2 = cSim.pShakeParm2[pos];
      if(shakeID.x < -1)
      {
        shakeID.x = abs(shakeID.x) - 2;
      }
      if(shakeID.y < -1)
      {
        shakeID.y = abs(shakeID.y) - 2;
      }
      if(shakeID.z < -1)
      {
        shakeID.z = abs(shakeID.z) - 2;
      }
      if(shakeID.w < -1)
      {
        shakeID.w = abs(shakeID.w) - 2;
      }
#endif
#ifdef SHAKE_HMR
      double invMassH = cSim.pShakeInvMassH[pos];
#endif

      // Read SHAKE network components
#ifdef NODPTEXTURE
        double xi = cSim.pOldAtomX[shakeID.x];
        double yi = cSim.pOldAtomY[shakeID.x];
        double zi = cSim.pOldAtomZ[shakeID.x];
        double xij = cSim.pOldAtomX[shakeID.y];
        double yij = cSim.pOldAtomY[shakeID.y];
        double zij = cSim.pOldAtomZ[shakeID.y];
#else
        int2 ixi = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID.x);
        int2 iyi = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID.x + cSim.stride);
        int2 izi = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID.x + cSim.stride2);
        int2 ixij = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID.y);
        int2 iyij = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID.y + cSim.stride);
        int2 izij = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID.y + cSim.stride2);
        double xi = __hiloint2double(ixi.y, ixi.x);
        double yi = __hiloint2double(iyi.y, iyi.x);
        double zi = __hiloint2double(izi.y, izi.x);
        double xij = __hiloint2double(ixij.y, ixij.x);
        double yij = __hiloint2double(iyij.y, iyij.x);
        double zij = __hiloint2double(izij.y, izij.x);
#endif
        double xpi = PATOMX(shakeID.x);
        double ypi = PATOMY(shakeID.x);
        double zpi = PATOMZ(shakeID.x);
        double xpj = PATOMX(shakeID.y);
        double ypj = PATOMY(shakeID.y);
        double zpj = PATOMZ(shakeID.y);
        double invMassI = shakeParm.x;
        double toler = shakeParm.y;

      // Optionally read 2nd hydrogen
      double xpk, ypk, zpk, xik, yik, zik;
      if (shakeID.z != -1) {
#ifdef NODPTEXTURE
        xik = cSim.pOldAtomX[shakeID.z];
        yik = cSim.pOldAtomY[shakeID.z];
        zik = cSim.pOldAtomZ[shakeID.z];
#else
        int2 ixik = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID.z);
        int2 iyik = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID.z + cSim.stride);
        int2 izik = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID.z + cSim.stride2);
        xik = __hiloint2double(ixik.y, ixik.x);
        yik = __hiloint2double(iyik.y, iyik.x);
        zik = __hiloint2double(izik.y, izik.x);
#endif
        xpk = PATOMX(shakeID.z);
        ypk = PATOMY(shakeID.z);
        zpk = PATOMZ(shakeID.z);
      }

      // Optionally read 3rd hydrogen into shared memory
      double xpl, ypl, zpl, xil, yil, zil;
      if (shakeID.w != -1) {
#ifdef NODPTEXTURE
        xil = cSim.pOldAtomX[shakeID.w];
        yil = cSim.pOldAtomY[shakeID.w];
        zil = cSim.pOldAtomZ[shakeID.w];
#else
        int2 ixil = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID.w);
        int2 iyil = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID.w + cSim.stride);
        int2 izil = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID.w + cSim.stride2);
        xil = __hiloint2double(ixil.y, ixil.x);
        yil = __hiloint2double(iyil.y, iyil.x);
        zil = __hiloint2double(izil.y, izil.x);
#endif
        xpl = PATOMX(shakeID.w);
        ypl = PATOMY(shakeID.w);
        zpl = PATOMZ(shakeID.w);
      }

      // Calculate unchanging quantities
      if (shakeID.y != -1) {
        xij = xi - xij;
        yij = yi - yij;
        zij = zi - zij;
        if (shakeID.z != -1) {
          xik = xi - xik;
          yik = yi - yik;
          zik = zi - zik;
        }
        if (shakeID.w != -1) {
          xil = xi - xil;
          yil = yi - yil;
          zil = zi - zil;
        }
      }
      bool done = false;

      for (int i = 0; i < 3000; i++) {
        done = true;
        if (shakeID.y == -1) break;

        // Calculate nominal distance squared
        double xpxx = xpi - xpj;
        double ypxx = ypi - ypj;
        double zpxx = zpi - zpj;
        double rpxx2 = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;

#ifdef TISHAKE2
        if(refShakeID.y < -1)
        {
          invMassI = shakeParm2.x;
          toler    = shakeParm2.y;
        }
        else
        {
          invMassI = shakeParm.x;
          toler    = shakeParm.y;
        }
#endif

        // Apply correction
        double diff = toler - rpxx2;
        if (abs(diff) >= toler * cSim.tol) {
          done = false;

          // Shake resetting of coordinate is done here
          double rrpr = xij * xpxx + yij * ypxx + zij * zpxx;
          if (rrpr >= toler * (double)1.0e-06) {
            double acor = diff / (rrpr * (double)2.0 * (invMassI + INVMASSH));
            double h = xij * acor;
            xpi += h * invMassI;
            xpj -= h * INVMASSH;
            h = yij * acor;
            ypi += h * invMassI;
            ypj -= h * INVMASSH;
            h = zij * acor;
            zpi += h * invMassI;
            zpj -= h * INVMASSH;
          }
        }
        

        // Second bond if present
        if (shakeID.z != -1) {
          xpxx  = xpi - xpk;
          ypxx  = ypi - ypk;
          zpxx  = zpi - zpk;
          rpxx2 = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;

#ifdef TISHAKE2
          if(refShakeID.z < -1)
          {
            invMassI = shakeParm2.x;
            toler    = shakeParm2.y;
          }
          else
          {
            invMassI = shakeParm.x;
            toler    = shakeParm.y;
          }
#endif

          // Apply correction
          diff = toler - rpxx2;
          if (abs(diff) >= toler * cSim.tol) {
            done = false;

            // Shake resetting of coordinate is done here
            double rrpr = xik * xpxx + yik * ypxx + zik * zpxx;
            if (rrpr >= toler * (double)1.0e-06) {
              double acor = diff / (rrpr * (double)2.0 * (invMassI + INVMASSH));
              double h    = xik * acor;
              xpi += h * invMassI;
              xpk -= h * INVMASSH;
              h    = yik * acor;
              ypi += h * invMassI;
              ypk -= h * INVMASSH;
              h    = zik * acor;
              zpi += h * invMassI;
              zpk -= h * INVMASSH;
            }
          }
        }

        // Third bond if present
        if (shakeID.w != -1) {
          xpxx  = xpi - xpl;
          ypxx  = ypi - ypl;
          zpxx  = zpi - zpl;
          rpxx2 = xpxx * xpxx + ypxx * ypxx + zpxx * zpxx;

#ifdef TISHAKE2
          if(refShakeID.w < -1)
          {
            invMassI = shakeParm2.x;
            toler    = shakeParm2.y;
          }
          else
          {
            invMassI = shakeParm.x;
            toler    = shakeParm.y;
          }
#endif

          // Apply correction
          diff = toler - rpxx2;
          if (abs(diff) >= toler * cSim.tol) {
            done = false;

            // Shake resetting of coordinate is done here
            double rrpr = xil * xpxx + yil * ypxx + zil * zpxx;
            if (rrpr >= toler * (double)1.0e-06) {
              double acor = diff / (rrpr * (double)2.0 * (invMassI + INVMASSH));
              double h    = xil * acor;
              xpi += h * invMassI;
              xpl -= h * INVMASSH;
              h    = yil * acor;
              ypi += h * invMassI;
              ypl -= h * INVMASSH;
              h    = zil * acor;
              zpi += h * invMassI;
              zpl -= h * INVMASSH;
            }
          }
        }

        // Check for convergence
        if (done) {
          break;
        }
      }

      // Write out results if converged, but there's no really good
      // way to indicate failure so we'll let the simulation heading
      // off to Neptune do that for us.  Wish there were a better way,
      // but until the CPU needs something from the GPU, those are the
      // the breaks.  I guess, technically, we could just set a flag to NOP
      // the simulation from here and then carry that result through upon
      // the next ntpr, ntwc, or ntwx update, but I leave that up to you
      // guys to implement that (or not).
      if (done) {
        if (shakeID.y != -1) {
          PATOMX(shakeID.x) = xpi;
          PATOMY(shakeID.x) = ypi;
          PATOMZ(shakeID.x) = zpi;
#ifndef SHAKE_NEIGHBORLIST
          PMEFloat2 xyi = { (PMEFloat)xpi, (PMEFloat)ypi };
          cSim.pAtomXYSP[shakeID.x] = xyi;
          cSim.pAtomZSP[shakeID.x] = zpi;
#endif
          PATOMX(shakeID.y) = xpj;
          PATOMY(shakeID.y) = ypj;
          PATOMZ(shakeID.y) = zpj;
#ifndef SHAKE_NEIGHBORLIST
          PMEFloat2 xyj = { (PMEFloat)xpj, (PMEFloat)ypj };
          cSim.pAtomXYSP[shakeID.y] = xyj;
          cSim.pAtomZSP[shakeID.y] = zpj;
#endif
        }
        if (shakeID.z != -1) {
          PATOMX(shakeID.z) = xpk;
          PATOMY(shakeID.z) = ypk;
          PATOMZ(shakeID.z) = zpk;
#ifndef SHAKE_NEIGHBORLIST
          PMEFloat2 xyk = {(PMEFloat)xpk, (PMEFloat)ypk};
          cSim.pAtomXYSP[shakeID.z] = xyk;
          cSim.pAtomZSP[shakeID.z]  = zpk;
#endif
        }
        if (shakeID.w != -1) {
          PATOMX(shakeID.w)           = xpl;
          PATOMY(shakeID.w)           = ypl;
          PATOMZ(shakeID.w)           = zpl;
#ifndef SHAKE_NEIGHBORLIST
          PMEFloat2 xyl               = {(PMEFloat)xpl, (PMEFloat)ypl};
          cSim.pAtomXYSP[shakeID.w]   = xyl;
          cSim.pAtomZSP[shakeID.w]    = zpl;
#endif
        }
      }
    }
  }
  else if (pos < cSim.fastShakeOffset) {
    pos -= cSim.shakeOffset;
    if (pos < cSim.fastShakeConstraints) {

      // Read atom data
#ifdef SHAKE_NEIGHBORLIST
      int4 shakeID = cSim.pImageFastShakeID[pos];
#else
      int4 shakeID = cSim.pFastShakeID[pos];
#endif
#ifdef NODPTEXTURE
      double x1 = cSim.pOldAtomX[shakeID.x];
      double y1 = cSim.pOldAtomY[shakeID.x];
      double z1 = cSim.pOldAtomZ[shakeID.x];
      double x2 = cSim.pOldAtomX[shakeID.y];
      double y2 = cSim.pOldAtomY[shakeID.y];
      double z2 = cSim.pOldAtomZ[shakeID.y];
      double x3 = cSim.pOldAtomX[shakeID.z];
      double y3 = cSim.pOldAtomY[shakeID.z];
      double z3 = cSim.pOldAtomZ[shakeID.z];
#else
      int2 ix1  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID.x);
      int2 iy1  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID.x + cSim.stride);
      int2 iz1  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID.x + cSim.stride2);
      int2 ix2  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID.y);
      int2 iy2  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID.y + cSim.stride);
      int2 iz2  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID.y + cSim.stride2);
      int2 ix3  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID.z);
      int2 iy3  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID.z + cSim.stride);
      int2 iz3  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID.z + cSim.stride2);
      double x1 = __hiloint2double(ix1.y, ix1.x);
      double y1 = __hiloint2double(iy1.y, iy1.x);
      double z1 = __hiloint2double(iz1.y, iz1.x);
      double x2 = __hiloint2double(ix2.y, ix2.x);
      double y2 = __hiloint2double(iy2.y, iy2.x);
      double z2 = __hiloint2double(iz2.y, iz2.x);
      double x3 = __hiloint2double(ix3.y, ix3.x);
      double y3 = __hiloint2double(iy3.y, iy3.x);
      double z3 = __hiloint2double(iz3.y, iz3.x);
#endif
      double xp1 = PATOMX(shakeID.x);
      double yp1 = PATOMY(shakeID.x);
      double zp1 = PATOMZ(shakeID.x);
      double xp2 = PATOMX(shakeID.y);
      double yp2 = PATOMY(shakeID.y);
      double zp2 = PATOMZ(shakeID.y);
      double xp3 = PATOMX(shakeID.z);
      double yp3 = PATOMY(shakeID.z);
      double zp3 = PATOMZ(shakeID.z);

      // Step1  A1_prime:
      double xb0  = x2 - x1;
      double yb0  = y2 - y1;
      double zb0  = z2 - z1;
      double xc0  = x3 - x1;
      double yc0  = y3 - y1;
      double zc0  = z3 - z1;
      double xcom = xp1*cSim.wo_div_wohh + (xp2 + xp3)*cSim.wh_div_wohh;
      double ycom = yp1*cSim.wo_div_wohh + (yp2 + yp3)*cSim.wh_div_wohh;
      double zcom = zp1*cSim.wo_div_wohh + (zp2 + zp3)*cSim.wh_div_wohh;
      double xa1 = xp1 - xcom;
      double ya1 = yp1 - ycom;
      double za1 = zp1 - zcom;
      double xb1 = xp2 - xcom;
      double yb1 = yp2 - ycom;
      double zb1 = zp2 - zcom;
      double xc1 = xp3 - xcom;
      double yc1 = yp3 - ycom;
      double zc1 = zp3 - zcom;
      double xakszd = yb0*zc0 - zb0*yc0;
      double yakszd = zb0*xc0 - xb0*zc0;
      double zakszd = xb0*yc0 - yb0*xc0;
      double xaksxd = ya1*zakszd - za1*yakszd;
      double yaksxd = za1*xakszd - xa1*zakszd;
      double zaksxd = xa1*yakszd - ya1*xakszd;
      double xaksyd = yakszd*zaksxd - zakszd*yaksxd;
      double yaksyd = zakszd*xaksxd - xakszd*zaksxd;
      double zaksyd = xakszd*yaksxd - yakszd*xaksxd;
      double axlng_inv = rsqrt(xaksxd*xaksxd + yaksxd*yaksxd + zaksxd*zaksxd);
      double aylng_inv = rsqrt(xaksyd*xaksyd + yaksyd*yaksyd + zaksyd*zaksyd);
      double azlng_inv = rsqrt(xakszd*xakszd + yakszd*yakszd + zakszd*zakszd);
      double trns11 = xaksxd * axlng_inv;
      double trns21 = yaksxd * axlng_inv;
      double trns31 = zaksxd * axlng_inv;
      double trns12 = xaksyd * aylng_inv;
      double trns22 = yaksyd * aylng_inv;
      double trns32 = zaksyd * aylng_inv;
      double trns13 = xakszd * azlng_inv;
      double trns23 = yakszd * azlng_inv;
      double trns33 = zakszd * azlng_inv;
      double xb0d = trns11*xb0 + trns21*yb0 + trns31*zb0;
      double yb0d = trns12*xb0 + trns22*yb0 + trns32*zb0;
      double xc0d = trns11*xc0 + trns21*yc0 + trns31*zc0;
      double yc0d = trns12*xc0 + trns22*yc0 + trns32*zc0;
      double za1d = trns13*xa1 + trns23*ya1 + trns33*za1;
      double xb1d = trns11*xb1 + trns21*yb1 + trns31*zb1;
      double yb1d = trns12*xb1 + trns22*yb1 + trns32*zb1;
      double zb1d = trns13*xb1 + trns23*yb1 + trns33*zb1;
      double xc1d = trns11*xc1 + trns21*yc1 + trns31*zc1;
      double yc1d = trns12*xc1 + trns22*yc1 + trns32*zc1;
      double zc1d = trns13*xc1 + trns23*yc1 + trns33*zc1;

      // Step2  A2_prime:
      double sinphi = za1d * cSim.ra_inv;
      double cosphi = sqrt(1.0 - sinphi*sinphi);
      double sinpsi = (zb1d - zc1d) / (cSim.rc2 * cosphi);
      double cospsi = sqrt(1.0 - sinpsi * sinpsi);
      double ya2d =  cSim.ra * cosphi;
      double xb2d = -cSim.rc * cospsi;
      double yb2d = -cSim.rb*cosphi - cSim.rc*sinpsi*sinphi;
      double yc2d = -cSim.rb*cosphi + cSim.rc*sinpsi*sinphi;
      xb2d = -0.5 * sqrt(cSim.hhhh - (yb2d - yc2d)*(yb2d - yc2d) -
                         (zb1d - zc1d)*(zb1d - zc1d));

      // Step3  al,be,ga:
      double alpa = (xb2d * (xb0d-xc0d) + yb0d * yb2d + yc0d * yc2d);
      double beta = (xb2d * (yc0d-yb0d) + xb0d * yb2d + xc0d * yc2d);
      double gama = xb0d * yb1d - xb1d * yb0d + xc0d * yc1d - xc1d * yc0d;

      double al2be2 =  alpa*alpa + beta*beta;
      double sinthe = (alpa*gama - beta*sqrt(al2be2 - gama*gama)) / al2be2;

      // Step4  A3_prime:
      double costhe = sqrt(1.0 - sinthe*sinthe);
      double xa3d   = -ya2d * sinthe;
      double ya3d   =  ya2d * costhe;
      double za3d   =  za1d;
      double xb3d   =  xb2d*costhe - yb2d*sinthe;
      double yb3d   =  xb2d*sinthe + yb2d*costhe;
      double zb3d   =  zb1d;
      double xc3d   = -xb2d*costhe - yc2d*sinthe;
      double yc3d   = -xb2d*sinthe + yc2d*costhe;
      double zc3d   =  zc1d;

      // Step5  A3:
      PATOMX(shakeID.x) = xcom + trns11*xa3d + trns12*ya3d + trns13*za3d;
      PATOMY(shakeID.x) = ycom + trns21*xa3d + trns22*ya3d + trns23*za3d;
      PATOMZ(shakeID.x) = zcom + trns31*xa3d + trns32*ya3d + trns33*za3d;
      PATOMX(shakeID.y) = xcom + trns11*xb3d + trns12*yb3d + trns13*zb3d;
      PATOMY(shakeID.y) = ycom + trns21*xb3d + trns22*yb3d + trns23*zb3d;
      PATOMZ(shakeID.y) = zcom + trns31*xb3d + trns32*yb3d + trns33*zb3d;
      PATOMX(shakeID.z) = xcom + trns11*xc3d + trns12*yc3d + trns13*zc3d;
      PATOMY(shakeID.z) = ycom + trns21*xc3d + trns22*yc3d + trns23*zc3d;
      PATOMZ(shakeID.z) = zcom + trns31*xc3d + trns32*yc3d + trns33*zc3d;
    }
  }
  else if (pos < cSim.slowShakeOffset) {
    pos -= cSim.fastShakeOffset;
    if (pos < cSim.slowShakeConstraints) {

      // Read SHAKE network data
#ifdef SHAKE_NEIGHBORLIST
      int  shakeID1 = cSim.pImageSlowShakeID1[pos];
      int4 shakeID2 = cSim.pImageSlowShakeID2[pos];
#ifdef TISHAKE2
      int4 refShakeID  = cSim.pImageSlowShakeID2[pos];
#endif
#else
      int  shakeID1 = cSim.pSlowShakeID1[pos];
      int4 shakeID2 = cSim.pSlowShakeID2[pos];
#ifdef TISHAKE2
      int4 refShakeID = cSim.pSlowShakeID2[pos];
#endif
#endif
      double2 shakeParm = cSim.pSlowShakeParm[pos];
#ifdef TISHAKE2
      double2 shakeParm2 = cSim.pSlowShakeParm2[pos];
      if(shakeID2.x < -1)
      {
          shakeID2.x = abs(shakeID2.x) - 2;
      }
      if(shakeID2.y < -1)
      {
          shakeID2.y = abs(shakeID2.y) - 2;
      }
      if(shakeID2.z < -1)
      {
          shakeID2.z = abs(shakeID2.z) - 2;
      }
      if(shakeID2.w < -1)
      {
          shakeID2.w = abs(shakeID2.w) - 2;
      }
#endif
#ifdef SHAKE_HMR
      double invMassH = cSim.pSlowShakeInvMassH[pos];
#endif
      // Read SHAKE network components
#ifdef NODPTEXTURE
      double xi  = cSim.pOldAtomX[shakeID1];
      double yi  = cSim.pOldAtomY[shakeID1];
      double zi  = cSim.pOldAtomZ[shakeID1];
      double xij = cSim.pOldAtomX[shakeID2.x];
      double yij = cSim.pOldAtomY[shakeID2.x];
      double zij = cSim.pOldAtomZ[shakeID2.x];
      double xik = cSim.pOldAtomX[shakeID2.y];
      double yik = cSim.pOldAtomY[shakeID2.y];
      double zik = cSim.pOldAtomZ[shakeID2.y];
      double xil = cSim.pOldAtomX[shakeID2.z];
      double yil = cSim.pOldAtomY[shakeID2.z];
      double zil = cSim.pOldAtomZ[shakeID2.z];
      double xim = cSim.pOldAtomX[shakeID2.w];
      double yim = cSim.pOldAtomY[shakeID2.w];
      double zim = cSim.pOldAtomZ[shakeID2.w];
#else
      int2 ixi   = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID1);
      int2 iyi   = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID1 + cSim.stride);
      int2 izi   = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID1 + cSim.stride2);
      int2 ixij  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.x);
      int2 iyij  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.x + cSim.stride);
      int2 izij  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.x + cSim.stride2);
      int2 ixik  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.y);
      int2 iyik  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.y + cSim.stride);
      int2 izik  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.y + cSim.stride2);
      int2 ixil  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.z);
      int2 iyil  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.z + cSim.stride);
      int2 izil  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.z + cSim.stride2);
      int2 ixim  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.w);
      int2 iyim  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.w + cSim.stride);
      int2 izim  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.w + cSim.stride2);
      double xi  = __hiloint2double(ixi.y, ixi.x);
      double yi  = __hiloint2double(iyi.y, iyi.x);
      double zi  = __hiloint2double(izi.y, izi.x);
      double xij = __hiloint2double(ixij.y, ixij.x);
      double yij = __hiloint2double(iyij.y, iyij.x);
      double zij = __hiloint2double(izij.y, izij.x);
      double xik = __hiloint2double(ixik.y, ixik.x);
      double yik = __hiloint2double(iyik.y, iyik.x);
      double zik = __hiloint2double(izik.y, izik.x);
      double xil = __hiloint2double(ixil.y, ixil.x);
      double yil = __hiloint2double(iyil.y, iyil.x);
      double zil = __hiloint2double(izil.y, izil.x);
      double xim = __hiloint2double(ixim.y, ixim.x);
      double yim = __hiloint2double(iyim.y, iyim.x);
      double zim = __hiloint2double(izim.y, izim.x);
#endif
      double xpi = PATOMX(shakeID1);
      double ypi = PATOMY(shakeID1);
      double zpi = PATOMZ(shakeID1);
      double xpj = PATOMX(shakeID2.x);
      double ypj = PATOMY(shakeID2.x);
      double zpj = PATOMZ(shakeID2.x);
      double xpk = PATOMX(shakeID2.y);
      double ypk = PATOMY(shakeID2.y);
      double zpk = PATOMZ(shakeID2.y);
      double xpl = PATOMX(shakeID2.z);
      double ypl = PATOMY(shakeID2.z);
      double zpl = PATOMZ(shakeID2.z);
      double xpm = PATOMX(shakeID2.w);
      double ypm = PATOMY(shakeID2.w);
      double zpm = PATOMZ(shakeID2.w);
      double invMassI = shakeParm.x;
      double toler    = shakeParm.y;

      // Calculate unchanging quantities
      xij = xi - xij;
      yij = yi - yij;
      zij = zi - zij;
      xik = xi - xik;
      yik = yi - yik;
      zik = zi - zik;
      xil = xi - xil;
      yil = yi - yil;
      zil = zi - zil;
      xim = xi - xim;
      yim = yi - yim;
      zim = zi - zim;

      bool done = false;
      for (int i = 0; i < 3000; i++) {
        done = true;

        // Calculate nominal distance squared
        double xpxx  = xpi - xpj;
        double ypxx  = ypi - ypj;
        double zpxx  = zpi - zpj;
        double rpxx2 = xpxx*xpxx + ypxx*ypxx + zpxx*zpxx;

#ifdef TISHAKE2
        if(refShakeID.x < -1)
        {
          invMassI = shakeParm2.x;
          toler    = shakeParm2.y;
        }
        else
        {
          invMassI = shakeParm.x;
          toler    = shakeParm.y;
        }

        // Apply correction to first hydrogen
        double diff = toler - rpxx2;
        if (abs(diff) >= toler * cSim.tol) {
          done = false;

          // Shake resetting of coordinate is done here
          double rrpr = xij * xpxx + yij * ypxx + zij * zpxx;
          if (rrpr >= toler * (double)1.0e-06) {
            double acor = diff / (rrpr * (double)2.0 * (invMassI + INVMASSH));
            double h    = xij * acor;
            xpi += h * invMassI;
            xpj -= h * INVMASSH;
            h    = yij * acor;
            ypi += h * invMassI;
            ypj -= h * INVMASSH;
            h    = zij * acor;
            zpi += h * invMassI;
            zpj -= h * INVMASSH;
          }
        }
        xpxx  = xpi - xpk;
        ypxx  = ypi - ypk;
        zpxx  = zpi - zpk;
        rpxx2 = xpxx*xpxx + ypxx*ypxx + zpxx*zpxx;

        if(refShakeID.y < -1)
        {
          invMassI = shakeParm2.x;
          toler    = shakeParm2.y;
        }
        else
        {
          invMassI = shakeParm.x;
          toler    = shakeParm.y;
        }

        // Apply correction to second hydrogen
        diff = toler - rpxx2;
        if (abs(diff) >= toler * cSim.tol) {
          done = false;

          // Shake resetting of coordinate is done here
          double rrpr = xik*xpxx + yik*ypxx + zik*zpxx;
          if (rrpr >= toler * (double)1.0e-06) {
            double acor             = diff / (rrpr * (double)2.0 * (invMassI + INVMASSH));
            double h                = xik * acor;
            xpi                    += h * invMassI;
            xpk                    -= h * INVMASSH;
            h                       = yik * acor;
            ypi                    += h * invMassI;
            ypk                    -= h * INVMASSH;
            h                       = zik * acor;
            zpi                    += h * invMassI;
            zpk                    -= h * INVMASSH;
          }
        }
        xpxx  = xpi - xpl;
        ypxx  = ypi - ypl;
        zpxx  = zpi - zpl;
        rpxx2 = xpxx*xpxx + ypxx*ypxx + zpxx*zpxx;

        if(refShakeID.z < -1)
        {
          invMassI = shakeParm2.x;
          toler    = shakeParm2.y;
        }
        else
        {
          invMassI = shakeParm.x;
          toler    = shakeParm.y;
        }

        // Apply correction to third hydrogen
        diff = toler - rpxx2;
        if (abs(diff) >= toler * cSim.tol) {
          done = false;

          // Shake resetting of coordinate is done here
          double rrpr = xil*xpxx + yil*ypxx + zil*zpxx;
          if (rrpr >= toler * (double)1.0e-06) {
            double acor = diff / (rrpr * (double)2.0 * (invMassI + INVMASSH));
            double h = xil * acor;
            xpi += h * invMassI;
            xpl -= h * INVMASSH;
            h    = yil * acor;
            ypi += h * invMassI;
            ypl -= h * INVMASSH;
            h    = zil * acor;
            zpi += h * invMassI;
            zpl -= h * INVMASSH;
          }
        }
        xpxx  = xpi - xpm;
        ypxx  = ypi - ypm;
        zpxx  = zpi - zpm;
        rpxx2 = xpxx*xpxx + ypxx*ypxx + zpxx*zpxx;

        if(refShakeID.w < -1)
        {
          invMassI = shakeParm2.x;
          toler    = shakeParm2.y;
        }
        else
        {
          invMassI = shakeParm.x;
          toler    = shakeParm.y;
        }

        // Apply correction to fourth hydrogen
        diff = toler - rpxx2;
        if (abs(diff) >= toler * cSim.tol) {
          done = false;

          // Shake resetting of coordinate is done here
          double rrpr = xim*xpxx + yim*ypxx + zim*zpxx;
          if (rrpr >= toler * (double)1.0e-06) {
            double acor = diff / (rrpr * (double)2.0 * (invMassI + INVMASSH));
            double h    = xim * acor;
            xpi += h * invMassI;
            xpm -= h * INVMASSH;
            h    = yim * acor;
            ypi += h * invMassI;
            ypm -= h * INVMASSH;
            h    = zim * acor;
            zpi += h * invMassI;
            zpm -= h * INVMASSH;
          }
        }

        // Check for convergence
        if (done) {
          break;
        }
      }

      // Write out results if converged, but there's no really good
      // way to indicate failure so we'll let the simulation heading
      // off to Neptune do that for us.  Wish there were a better way,
      // but until the CPU needs something from the GPU, those are the
      // the breaks.  I guess, technically, we could just set a flag to NOP
      // the simulation from here and then carry that result through upon
      // the next ntpr, ntwc, or ntwx update, but I leave that up to you
      // guys to implement that (or not).
      if (done) {
        PATOMX(shakeID1) = xpi;
        PATOMY(shakeID1) = ypi;
        PATOMZ(shakeID1) = zpi;
#ifndef SHAKE_NEIGHBORLIST
        PMEFloat2 xyi = {(PMEFloat)xpi, (PMEFloat)ypi};
        cSim.pAtomXYSP[shakeID1] = xyi;
        cSim.pAtomZSP[shakeID1]  = zpi;
#endif
        PATOMX(shakeID2.x) = xpj;
        PATOMY(shakeID2.x) = ypj;
        PATOMZ(shakeID2.x) = zpj;
#ifndef SHAKE_NEIGHBORLIST
        PMEFloat2 xyj = {(PMEFloat)xpj, (PMEFloat)ypj};
        cSim.pAtomXYSP[shakeID2.x] = xyj;
        cSim.pAtomZSP[shakeID2.x]  = zpj;
#endif
        PATOMX(shakeID2.y)              = xpk;
        PATOMY(shakeID2.y)              = ypk;
        PATOMZ(shakeID2.y)              = zpk;
#ifndef SHAKE_NEIGHBORLIST
        PMEFloat2 xyk = {(PMEFloat)xpk, (PMEFloat)ypk};
        cSim.pAtomXYSP[shakeID2.y] = xyk;
        cSim.pAtomZSP[shakeID2.y]  = zpk;
#endif
        PATOMX(shakeID2.z) = xpl;
        PATOMY(shakeID2.z) = ypl;
        PATOMZ(shakeID2.z) = zpl;
#ifndef SHAKE_NEIGHBORLIST
        PMEFloat2 xyl = {(PMEFloat)xpl, (PMEFloat)ypl};
        cSim.pAtomXYSP[shakeID2.z] = xyl;
        cSim.pAtomZSP[shakeID2.z]  = zpl;
#endif
        PATOMX(shakeID2.w) = xpm;
        PATOMY(shakeID2.w) = ypm;
        PATOMZ(shakeID2.w) = zpm;
#ifndef SHAKE_NEIGHBORLIST
        PMEFloat2 xym = {(PMEFloat)xpm, (PMEFloat)ypm};
        cSim.pAtomXYSP[shakeID2.w] = xym;
        cSim.pAtomZSP[shakeID2.w]  = zpm;
#endif
      }
    }
  }
  else if(pos < cSim.slowTIShakeOffset)
  {
    pos -= cSim.slowShakeOffset;
    if (pos < cSim.slowTIShakeConstraints) {

      // Read SHAKE network data
#ifdef SHAKE_NEIGHBORLIST
      int  shakeID1 = cSim.pImageSlowTIShakeID1[pos];
      int4 shakeID2 = cSim.pImageSlowTIShakeID2[pos];
      int4 shakeID3 = cSim.pImageSlowTIShakeID3[pos];
      int4 refShakeID2 = cSim.pImageSlowTIShakeID2[pos];
      int4 refShakeID3 = cSim.pImageSlowTIShakeID3[pos];
#else
      int  shakeID1 = cSim.pSlowTIShakeID1[pos];
      int4 shakeID2 = cSim.pSlowTIShakeID2[pos];
      int4 shakeID3 = cSim.pSlowTIShakeID2[pos];
      int4 refShakeID2 = cSim.pSlowTIShakeID2[pos];
      int4 refShakeID3 = cSim.pSlowTIShakeID3[pos];
#endif
      if(shakeID2.x < -1)
      {
          shakeID2.x = abs(shakeID2.x) - 2;
      }
      if(shakeID2.y < -1)
      {
          shakeID2.y = abs(shakeID2.y) - 2;
      }
      if(shakeID2.z < -1)
      {
          shakeID2.z = abs(shakeID2.z) - 2;
      }
      if(shakeID2.w < -1)
      {
          shakeID2.w = abs(shakeID2.w) - 2;
      }
      if(shakeID3.x < -1)
      {
          shakeID3.x = abs(shakeID3.x) - 2;
      }
      if(shakeID3.y < -1)
      {
          shakeID3.y = abs(shakeID3.y) - 2;
      }
      if(shakeID3.z < -1)
      {
          shakeID3.z = abs(shakeID3.z) - 2;
      }
      if(shakeID3.w < -1)
      {
          shakeID3.w = abs(shakeID3.w) - 2;
      }
      double2 shakeParm  = cSim.pSlowTIShakeParm[pos];
      double2 shakeParm2 = cSim.pSlowTIShakeParm2[pos];
#ifdef SHAKE_HMR
      double invMassH = cSim.pSlowTIShakeInvMassH[pos];
#endif

      // Read SHAKE network components
#ifdef NODPTEXTURE
      double xi  = cSim.pOldAtomX[shakeID1];
      double yi  = cSim.pOldAtomY[shakeID1];
      double zi  = cSim.pOldAtomZ[shakeID1];
      double xij = cSim.pOldAtomX[shakeID2.x];
      double yij = cSim.pOldAtomY[shakeID2.x];
      double zij = cSim.pOldAtomZ[shakeID2.x];
      double xik = cSim.pOldAtomX[shakeID2.y];
      double yik = cSim.pOldAtomY[shakeID2.y];
      double zik = cSim.pOldAtomZ[shakeID2.y];
      double xil = cSim.pOldAtomX[shakeID2.z];
      double yil = cSim.pOldAtomY[shakeID2.z];
      double zil = cSim.pOldAtomZ[shakeID2.z];
      double xim = cSim.pOldAtomX[shakeID2.w];
      double yim = cSim.pOldAtomY[shakeID2.w];
      double zim = cSim.pOldAtomZ[shakeID2.w];
      double xin = 0.0;
      double yin = 0.0;
      double zin = 0.0;
      double xio = 0.0;
      double yio = 0.0;
      double zio = 0.0;
      double xip = 0.0;
      double yip = 0.0;
      double zip = 0.0;
      double xiq = 0.0;
      double yiq = 0.0;
      double ziq = 0.0;
      if(shakeID3.x != -1) {
        xin = cSim.pOldAtomX[shakeID3.x];
        yin = cSim.pOldAtomY[shakeID3.x];
        zin = cSim.pOldAtomZ[shakeID3.x];
      }
      if(shakeID3.y != -1) {
        xio = cSim.pOldAtomX[shakeID3.y];
        yio = cSim.pOldAtomY[shakeID3.y];
        zio = cSim.pOldAtomZ[shakeID3.y];
      }
      if(shakeID3.z != -1) {
        xip = cSim.pOldAtomX[shakeID3.z];
        yip = cSim.pOldAtomY[shakeID3.z];
        zip = cSim.pOldAtomZ[shakeID3.z];
      }
      if(shakeID3.w != -1) {
        xiq = cSim.pOldAtomX[shakeID3.w];
        yiq = cSim.pOldAtomY[shakeID3.w];
        ziq = cSim.pOldAtomZ[shakeID3.w];
      }
#else
      int2 ixi   = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID1);
      int2 iyi   = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID1 + cSim.stride);
      int2 izi   = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID1 + cSim.stride2);
      int2 ixij  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.x);
      int2 iyij  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.x + cSim.stride);
      int2 izij  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.x + cSim.stride2);
      int2 ixik  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.y);
      int2 iyik  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.y + cSim.stride);
      int2 izik  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.y + cSim.stride2);
      int2 ixil  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.z);
      int2 iyil  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.z + cSim.stride);
      int2 izil  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.z + cSim.stride2);
      int2 ixim  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.w);
      int2 iyim  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.w + cSim.stride);
      int2 izim  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID2.w + cSim.stride2);
      double xi  = __hiloint2double(ixi.y, ixi.x);
      double yi  = __hiloint2double(iyi.y, iyi.x);
      double zi  = __hiloint2double(izi.y, izi.x);
      double xij = __hiloint2double(ixij.y, ixij.x);
      double yij = __hiloint2double(iyij.y, iyij.x);
      double zij = __hiloint2double(izij.y, izij.x);
      double xik = __hiloint2double(ixik.y, ixik.x);
      double yik = __hiloint2double(iyik.y, iyik.x);
      double zik = __hiloint2double(izik.y, izik.x);
      double xil = __hiloint2double(ixil.y, ixil.x);
      double yil = __hiloint2double(iyil.y, iyil.x);
      double zil = __hiloint2double(izil.y, izil.x);
      double xim = __hiloint2double(ixim.y, ixim.x);
      double yim = __hiloint2double(iyim.y, iyim.x);
      double zim = __hiloint2double(izim.y, izim.x);
      double xin = 0.0;
      double yin = 0.0;
      double zin = 0.0;
      double xio = 0.0;
      double yio = 0.0;
      double zio = 0.0;
      double xip = 0.0;
      double yip = 0.0;
      double zip = 0.0;
      double xiq = 0.0;
      double yiq = 0.0;
      double ziq = 0.0;
      if(shakeID3.x != -1) {
        int2 ixin  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID3.x);
        int2 iyin  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID3.x + cSim.stride);
        int2 izin  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID3.x + cSim.stride2);
        xin = __hiloint2double(ixin.y, ixin.x);
        yin = __hiloint2double(iyin.y, iyin.x);
        zin = __hiloint2double(izin.y, izin.x);
      }
      if(shakeID3.y != -1) {
        int2 ixio  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID3.y);
        int2 iyio  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID3.y + cSim.stride);
        int2 izio  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID3.y + cSim.stride2);
        xio = __hiloint2double(ixio.y, ixio.x);
        yio = __hiloint2double(iyio.y, iyio.x);
        zio = __hiloint2double(izio.y, izio.x);
      }
      if(shakeID3.z != -1) {
        int2 ixip  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID3.z);
        int2 iyip  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID3.z + cSim.stride);
        int2 izip  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID3.z + cSim.stride2);
        xip = __hiloint2double(ixip.y, ixip.x);
        yip = __hiloint2double(iyip.y, iyip.x);
        zip = __hiloint2double(izip.y, izip.x);
      }
      if(shakeID3.w != -1) {
        int2 ixiq  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID3.w);
        int2 iyiq  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID3.w + cSim.stride);
        int2 iziq  = tex1Dfetch<int2>(cSim.texOldAtomX, shakeID3.w + cSim.stride2);
        xiq = __hiloint2double(ixiq.y, ixiq.x);
        yiq = __hiloint2double(iyiq.y, iyiq.x);
        ziq = __hiloint2double(iziq.y, iziq.x);
      }
#endif
      double xpi = PATOMX(shakeID1);
      double ypi = PATOMY(shakeID1);
      double zpi = PATOMZ(shakeID1);
      double xpj = PATOMX(shakeID2.x);
      double ypj = PATOMY(shakeID2.x);
      double zpj = PATOMZ(shakeID2.x);
      double xpk = PATOMX(shakeID2.y);
      double ypk = PATOMY(shakeID2.y);
      double zpk = PATOMZ(shakeID2.y);
      double xpl = PATOMX(shakeID2.z);
      double ypl = PATOMY(shakeID2.z);
      double zpl = PATOMZ(shakeID2.z);
      double xpm = PATOMX(shakeID2.w);
      double ypm = PATOMY(shakeID2.w);
      double zpm = PATOMZ(shakeID2.w);
      double xpn = -1;
      double ypn = -1;
      double zpn = -1;
      double xpo = -1;
      double ypo = -1;
      double zpo = -1;
      double xpp = -1;
      double ypp = -1;
      double zpp = -1;
      double xpq = -1;
      double ypq = -1;
      double zpq = -1;
      if(shakeID3.x != -1) {
        xpn = PATOMX(shakeID3.x);
        ypn = PATOMY(shakeID3.x);
        zpn = PATOMZ(shakeID3.x);
      }
      if(shakeID3.y != -1) {
        xpo = PATOMX(shakeID3.y);
        ypo = PATOMY(shakeID3.y);
        zpo = PATOMZ(shakeID3.y);
      }
      if(shakeID3.z != -1) {
        xpp = PATOMX(shakeID3.z);
        ypp = PATOMY(shakeID3.z);
        zpp = PATOMZ(shakeID3.z);
      }
      if(shakeID3.w != -1) {
        xpq = PATOMX(shakeID3.w);
        ypq = PATOMY(shakeID3.w);
        zpq = PATOMZ(shakeID3.w);
      }
      double invMassI = shakeParm.x;
      double toler    = shakeParm.y;

      // Calculate unchanging quantities
      xij = xi - xij;
      yij = yi - yij;
      zij = zi - zij;
      xik = xi - xik;
      yik = yi - yik;
      zik = zi - zik;
      xil = xi - xil;
      yil = yi - yil;
      zil = zi - zil;
      xim = xi - xim;
      yim = yi - yim;
      zim = zi - zim;
      xin = xi - xin;
      yin = yi - yin;
      zin = zi - zin;
      if(shakeID3.y != -1) {
        xio = xi - xio;
        yio = yi - yio;
        zio = zi - zio;
      }
      if(shakeID3.z != -1) {
        xip = xi - xip;
        yip = yi - yip;
        zip = zi - zip;
      }
      if(shakeID3.w != -1) {
        xiq = xi - xiq;
        yiq = yi - yiq;
        ziq = zi - ziq;
      }
      
      bool done = false;
      for (int i = 0; i < 3000; i++) {
        done = true;

        // Calculate nominal distance squared
        double xpxx  = xpi - xpj;
        double ypxx  = ypi - ypj;
        double zpxx  = zpi - zpj;
        double rpxx2 = xpxx*xpxx + ypxx*ypxx + zpxx*zpxx;

        if(refShakeID2.x < -1)
        {
          invMassI = shakeParm2.x;
          toler    = shakeParm2.y;
        }
        else
        {
          invMassI = shakeParm.x;
          toler    = shakeParm.y;
        }
#endif //TISHAKE2

        // Apply correction to first hydrogen
        double diff = toler - rpxx2;
        if (abs(diff) >= toler * cSim.tol) {
          done = false;

          // Shake resetting of coordinate is done here
          double rrpr = xij * xpxx + yij * ypxx + zij * zpxx;
          if (rrpr >= toler * (double)1.0e-06) {
            double acor = diff / (rrpr * (double)2.0 * (invMassI + INVMASSH));
            double h    = xij * acor;
            xpi += h * invMassI;
            xpj -= h * INVMASSH;
            h    = yij * acor;
            ypi += h * invMassI;
            ypj -= h * INVMASSH;
            h    = zij * acor;
            zpi += h * invMassI;
            zpj -= h * INVMASSH;
          }
        }
        xpxx  = xpi - xpk;
        ypxx  = ypi - ypk;
        zpxx  = zpi - zpk;
        rpxx2 = xpxx*xpxx + ypxx*ypxx + zpxx*zpxx;

#ifdef TISHAKE2
        if(refShakeID2.y < -1)
        {
          invMassI = shakeParm2.x;
          toler    = shakeParm2.y;
        }
        else
        {
          invMassI = shakeParm.x;
          toler    = shakeParm.y;
        }
#endif

        // Apply correction to second hydrogen
        diff = toler - rpxx2;
        if (abs(diff) >= toler * cSim.tol) {
          done = false;

          // Shake resetting of coordinate is done here
          double rrpr = xik*xpxx + yik*ypxx + zik*zpxx;
          if (rrpr >= toler * (double)1.0e-06) {
            double acor             = diff / (rrpr * (double)2.0 * (invMassI + INVMASSH));
            double h                = xik * acor;
            xpi                    += h * invMassI;
            xpk                    -= h * INVMASSH;
            h                       = yik * acor;
            ypi                    += h * invMassI;
            ypk                    -= h * INVMASSH;
            h                       = zik * acor;
            zpi                    += h * invMassI;
            zpk                    -= h * INVMASSH;
          }
        }
        xpxx  = xpi - xpl;
        ypxx  = ypi - ypl;
        zpxx  = zpi - zpl;
        rpxx2 = xpxx*xpxx + ypxx*ypxx + zpxx*zpxx;

#ifdef TISHAKE2
        if(refShakeID2.z < -1)
        {
          invMassI = shakeParm2.x;
          toler    = shakeParm2.y;
        }
        else
        {
          invMassI = shakeParm.x;
          toler    = shakeParm.y;
        }
#endif

        // Apply correction to third hydrogen
        diff = toler - rpxx2;
        if (abs(diff) >= toler * cSim.tol) {
          done = false;

          // Shake resetting of coordinate is done here
          double rrpr = xil*xpxx + yil*ypxx + zil*zpxx;
          if (rrpr >= toler * (double)1.0e-06) {
            double acor = diff / (rrpr * (double)2.0 * (invMassI + INVMASSH));
            double h = xil * acor;
            xpi += h * invMassI;
            xpl -= h * INVMASSH;
            h    = yil * acor;
            ypi += h * invMassI;
            ypl -= h * INVMASSH;
            h    = zil * acor;
            zpi += h * invMassI;
            zpl -= h * INVMASSH;
          }
        }
        xpxx  = xpi - xpm;
        ypxx  = ypi - ypm;
        zpxx  = zpi - zpm;
        rpxx2 = xpxx*xpxx + ypxx*ypxx + zpxx*zpxx;

#ifdef TISHAKE2
        if(refShakeID2.w < -1)
        {
          invMassI = shakeParm2.x;
          toler    = shakeParm2.y;
        }
        else
        {
          invMassI = shakeParm.x;
          toler    = shakeParm.y;
        }
#endif

        // Apply correction to fourth hydrogen
        diff = toler - rpxx2;
        if (abs(diff) >= toler * cSim.tol) {
          done = false;

          // Shake resetting of coordinate is done here
          double rrpr = xim*xpxx + yim*ypxx + zim*zpxx;
          if (rrpr >= toler * (double)1.0e-06) {
            double acor = diff / (rrpr * (double)2.0 * (invMassI + INVMASSH));
            double h    = xim * acor;
            xpi += h * invMassI;
            xpm -= h * INVMASSH;
            h    = yim * acor;
            ypi += h * invMassI;
            ypm -= h * INVMASSH;
            h    = zim * acor;
            zpi += h * invMassI;
            zpm -= h * INVMASSH;
          }
        }
#ifdef TISHAKE2
        if(shakeID3.x != -1) {
          xpxx  = xpi - xpn;
          ypxx  = ypi - ypn;
          zpxx  = zpi - zpn;
          rpxx2 = xpxx*xpxx + ypxx*ypxx + zpxx*zpxx;

          if(refShakeID3.x < -1)
          {
            invMassI = shakeParm2.x;
            toler    = shakeParm2.y;
          }
          else
          {
            invMassI = shakeParm.x;
            toler    = shakeParm.y;
          }

          // Apply correction to fifth hydrogen
          diff = toler - rpxx2;
          if (abs(diff) >= toler * cSim.tol) {
            done = false;

            // Shake resetting of coordinate is done here
            double rrpr = xin*xpxx + yin*ypxx + zin*zpxx;
            if (rrpr >= toler * (double)1.0e-06) {
              double acor = diff / (rrpr * (double)2.0 * (invMassI + INVMASSH));
              double h    = xin * acor;
              xpi += h * invMassI;
              xpn -= h * INVMASSH;
              h    = yin * acor;
              ypi += h * invMassI;
              ypn -= h * INVMASSH;
              h    = zin * acor;
              zpi += h * invMassI;
              zpn -= h * INVMASSH;
            }
          }
        }
        if(shakeID3.y != -1) {
          xpxx  = xpi - xpo;
          ypxx  = ypi - ypo;
          zpxx  = zpi - zpo;
          rpxx2 = xpxx*xpxx + ypxx*ypxx + zpxx*zpxx;

          if(refShakeID3.y < -1)
          {
            invMassI = shakeParm2.x;
            toler    = shakeParm2.y;
          }
          else
          {
            invMassI = shakeParm.x;
            toler    = shakeParm.y;
          }

          // Apply correction to sixth hydrogen
          diff = toler - rpxx2;
          if (abs(diff) >= toler * cSim.tol) {
            done = false;

            // Shake resetting of coordinate is done here
            double rrpr = xio*xpxx + yio*ypxx + zio*zpxx;
            if (rrpr >= toler * (double)1.0e-06) {
              double acor = diff / (rrpr * (double)2.0 * (invMassI + INVMASSH));
              double h    = xio * acor;
              xpi += h * invMassI;
              xpo -= h * INVMASSH;
              h    = yio * acor;
              ypi += h * invMassI;
              ypo -= h * INVMASSH;
              h    = zio * acor;
              zpi += h * invMassI;
              zpo -= h * INVMASSH;
            }
          }
        }
        if(shakeID3.z != -1) {
          xpxx  = xpi - xpp;
          ypxx  = ypi - ypp;
          zpxx  = zpi - zpp;
          rpxx2 = xpxx*xpxx + ypxx*ypxx + zpxx*zpxx;

          if(refShakeID3.z < -1)
          {
            invMassI = shakeParm2.x;
            toler    = shakeParm2.y;
          }
          else
          {
            invMassI = shakeParm.x;
            toler    = shakeParm.y;
          }

          // Apply correction to seventh hydrogen
          diff = toler - rpxx2;
          if (abs(diff) >= toler * cSim.tol) {
            done = false;

            // Shake resetting of coordinate is done here
            double rrpr = xip*xpxx + yip*ypxx + zip*zpxx;
            if (rrpr >= toler * (double)1.0e-06) {
              double acor = diff / (rrpr * (double)2.0 * (invMassI + INVMASSH));
              double h    = xip * acor;
              xpi += h * invMassI;
              xpp -= h * INVMASSH;
              h    = yip * acor;
              ypi += h * invMassI;
              ypp -= h * INVMASSH;
              h    = zip * acor;
              zpi += h * invMassI;
              zpp -= h * INVMASSH;
            }
          }
        }
        if(shakeID3.w != -1) {
          xpxx  = xpi - xpq;
          ypxx  = ypi - ypq;
          zpxx  = zpi - zpq;
          rpxx2 = xpxx*xpxx + ypxx*ypxx + zpxx*zpxx;

          if(refShakeID3.w < -1)
          {
            invMassI = shakeParm2.x;
            toler    = shakeParm2.y;
          }
          else
          {
            invMassI = shakeParm.x;
            toler    = shakeParm.y;
          }

          // Apply correction to eighth hydrogen
          diff = toler - rpxx2;
          if (abs(diff) >= toler * cSim.tol) {
            done = false;

            // Shake resetting of coordinate is done here
            double rrpr = xiq*xpxx + yiq*ypxx + ziq*zpxx;
            if (rrpr >= toler * (double)1.0e-06) {
              double acor = diff / (rrpr * (double)2.0 * (invMassI + INVMASSH));
              double h    = xiq * acor;
              xpi += h * invMassI;
              xpq -= h * INVMASSH;
              h    = yiq * acor;
              ypi += h * invMassI;
              ypq -= h * INVMASSH;
              h    = ziq * acor;
              zpi += h * invMassI;
              zpq -= h * INVMASSH;
            }
          }
        }
#endif //TISHAKE2

        // Check for convergence
        if (done) {
          break;
        }
      }

      // Write out results if converged, but there's no really good
      // way to indicate failure so we'll let the simulation heading
      // off to Neptune do that for us.  Wish there were a better way,
      // but until the CPU needs something from the GPU, those are the
      // the breaks.  I guess, technically, we could just set a flag to NOP
      // the simulation from here and then carry that result through upon
      // the next ntpr, ntwc, or ntwx update, but I leave that up to you
      // guys to implement that (or not).
      if (done) {
        PATOMX(shakeID1) = xpi;
        PATOMY(shakeID1) = ypi;
        PATOMZ(shakeID1) = zpi;
#ifndef SHAKE_NEIGHBORLIST
        PMEFloat2 xyi = {(PMEFloat)xpi, (PMEFloat)ypi};
        cSim.pAtomXYSP[shakeID1] = xyi;
        cSim.pAtomZSP[shakeID1]  = zpi;
#endif
        PATOMX(shakeID2.x) = xpj;
        PATOMY(shakeID2.x) = ypj;
        PATOMZ(shakeID2.x) = zpj;
#ifndef SHAKE_NEIGHBORLIST
        PMEFloat2 xyj = {(PMEFloat)xpj, (PMEFloat)ypj};
        cSim.pAtomXYSP[shakeID2.x] = xyj;
        cSim.pAtomZSP[shakeID2.x]  = zpj;
#endif
        PATOMX(shakeID2.y)              = xpk;
        PATOMY(shakeID2.y)              = ypk;
        PATOMZ(shakeID2.y)              = zpk;
#ifndef SHAKE_NEIGHBORLIST
        PMEFloat2 xyk = {(PMEFloat)xpk, (PMEFloat)ypk};
        cSim.pAtomXYSP[shakeID2.y] = xyk;
        cSim.pAtomZSP[shakeID2.y]  = zpk;
#endif
        PATOMX(shakeID2.z) = xpl;
        PATOMY(shakeID2.z) = ypl;
        PATOMZ(shakeID2.z) = zpl;
#ifndef SHAKE_NEIGHBORLIST
        PMEFloat2 xyl = {(PMEFloat)xpl, (PMEFloat)ypl};
        cSim.pAtomXYSP[shakeID2.z] = xyl;
        cSim.pAtomZSP[shakeID2.z]  = zpl;
#endif
        PATOMX(shakeID2.w) = xpm;
        PATOMY(shakeID2.w) = ypm;
        PATOMZ(shakeID2.w) = zpm;
#ifndef SHAKE_NEIGHBORLIST
        PMEFloat2 xym = {(PMEFloat)xpm, (PMEFloat)ypm};
        cSim.pAtomXYSP[shakeID2.w] = xym;
        cSim.pAtomZSP[shakeID2.w]  = zpm;
#endif
#ifdef TISHAKE2
        PATOMX(shakeID3.x) = xpn;
        PATOMY(shakeID3.x) = ypn;
        PATOMZ(shakeID3.x) = zpn;
#ifndef SHAKE_NEIGHBORLIST
        PMEFloat2 xyn = {(PMEFloat)xpn, (PMEFloat)ypn};
        cSim.pAtomXYSP[shakeID2.x] = xyn;
        cSim.pAtomZSP[shakeID2.x]  = zpn;
#endif
        if(shakeID3.y != -1) {
          PATOMX(shakeID3.y) = xpo;
          PATOMY(shakeID3.y) = ypo;
          PATOMZ(shakeID3.y) = zpo;
#ifndef SHAKE_NEIGHBORLIST
          PMEFloat2 xyo = {(PMEFloat)xpo, (PMEFloat)ypo};
          cSim.pAtomXYSP[shakeID2.y] = xyo;
          cSim.pAtomZSP[shakeID2.y]  = zpo;
#endif
        }
        if(shakeID3.z != -1) {
          PATOMX(shakeID3.z) = xpp;
          PATOMY(shakeID3.z) = ypp;
          PATOMZ(shakeID3.z) = zpp;
#ifndef SHAKE_NEIGHBORLIST
          PMEFloat2 xyp = {(PMEFloat)xpp, (PMEFloat)ypp};
          cSim.pAtomXYSP[shakeID2.z] = xyp;
          cSim.pAtomZSP[shakeID2.z]  = zpp;
#endif
        }
        if(shakeID3.w != -1) {
          PATOMX(shakeID3.w) = xpq;
          PATOMY(shakeID3.w) = ypq;
          PATOMZ(shakeID3.w) = zpq;
#ifndef SHAKE_NEIGHBORLIST
          PMEFloat2 xyq = {(PMEFloat)xpq, (PMEFloat)ypq};
          cSim.pAtomXYSP[shakeID2.w] = xyq;
          cSim.pAtomZSP[shakeID2.w]  = zpq;
#endif
        }
#endif
      }
    }
  }
#undef PATOMX
#undef PATOMY
#undef PATOMZ
#undef INVMASSH
}

