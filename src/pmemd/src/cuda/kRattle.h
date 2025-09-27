#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// Oct 2021, by zhf
//---------------------------------------------------------------------------------------------
{
#ifdef RATTLE_NEIGHBORLIST
#  define PATOMX(i) cSim.pImageX[i]
#  define PATOMY(i) cSim.pImageY[i]
#  define PATOMZ(i) cSim.pImageZ[i]
#  define VELX(i) cSim.pImageVelX[i]
#  define VELY(i) cSim.pImageVelY[i]
#  define VELZ(i) cSim.pImageVelZ[i]
#else
#  define PATOMX(i) cSim.pAtomX[i]
#  define PATOMY(i) cSim.pAtomY[i]
#  define PATOMZ(i) cSim.pAtomZ[i]
#  define VELX(i) cSim.pVelX[i]
#  define VELY(i) cSim.pVelY[i]
#  define VELZ(i) cSim.pVelZ[i]
#endif

#ifdef RATTLE_HMR
#define INVMASSH invMassH
#else
#define INVMASSH cSim.invMassH
#endif

    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

    if (pos < cSim.shakeOffset) {
        if (pos < cSim.shakeConstraints) {

            //Read SHAKE network data
#ifdef RATTLE_NEIGHBORLIST
            int4 shakeID = cSim.pImageShakeID[pos];
#else
            int4 shakeID = cSim.pShakeID[pos];
#endif
            double2 shakeParm = cSim.pShakeParm[pos];
#ifdef RATTLE_HMR
            double invMassH = cSim.pShakeInvMassH[pos];
#endif

            double xpi = PATOMX(shakeID.x);
            double ypi = PATOMY(shakeID.x);
            double zpi = PATOMZ(shakeID.x);
            double xpj = PATOMX(shakeID.y);
            double ypj = PATOMY(shakeID.y);
            double zpj = PATOMZ(shakeID.y);
            double vxpi = VELX(shakeID.x);
            double vypi = VELY(shakeID.x);
            double vzpi = VELZ(shakeID.x);
            double vxpj = VELX(shakeID.y);
            double vypj = VELY(shakeID.y);
            double vzpj = VELZ(shakeID.y);

            double invMassI = shakeParm.x;
            double toler = shakeParm.y;

            double xpk, ypk, zpk, vxpk, vypk, vzpk;

            if (shakeID.z != -1) {
                xpk = PATOMX(shakeID.z);
                ypk = PATOMY(shakeID.z);
                zpk = PATOMZ(shakeID.z);
                vxpk = VELX(shakeID.z);
                vypk = VELY(shakeID.z);
                vzpk = VELZ(shakeID.z);
            }

            double xpl, ypl, zpl, vxpl, vypl, vzpl;
            if (shakeID.w != -1) {
                xpl = PATOMX(shakeID.w);
                ypl = PATOMY(shakeID.w);
                zpl = PATOMZ(shakeID.w);
                vxpl = VELX(shakeID.w);
                vypl = VELY(shakeID.w);
                vzpl = VELZ(shakeID.w);
            }

            bool done = false;
            double tol2 = cSim.tol / dt;
            for (int i = 0; i < 3000; ++i) {
                done = true;

                //calculate rvdot
                double xpxx = xpi - xpj;
                double ypxx = ypi - ypj;
                double zpxx = zpi - zpj;
                double vxpxx = vxpi - vxpj;
                double vypxx = vypi - vypj;
                double vzpxx = vzpi - vzpj;
                double rvdot = xpxx * vxpxx + ypxx * vypxx + zpxx * vzpxx;

                double acor = -rvdot / (toler * (invMassI + INVMASSH));

                if (abs(acor) >= tol2) {
                    done = false;
                    double h = xpxx * acor;
                    vxpi += h * invMassI;
                    vxpj -= h * INVMASSH;
                    h = ypxx * acor;
                    vypi += h * invMassI;
                    vypj -= h * INVMASSH;
                    h = zpxx * acor;
                    vzpi += h * invMassI;
                    vzpj -= h * INVMASSH;
                }

                // Second bond if present
                if (shakeID.z != -1) {
                    xpxx = xpi - xpk;
                    ypxx = ypi - ypk;
                    zpxx = zpi - zpk;
                    vxpxx = vxpi - vxpk;
                    vypxx = vypi - vypk;
                    vzpxx = vzpi - vzpk;
                    rvdot = xpxx * vxpxx + ypxx * vypxx + zpxx * vzpxx;

                    acor = -rvdot / (toler * (invMassI + INVMASSH));

                    if (abs(acor) >= tol2) {
                        done = false;

                        double h = xpxx * acor;
                        vxpi += h * invMassI;
                        vxpk -= h * INVMASSH;
                        h = ypxx * acor;
                        vypi += h * invMassI;
                        vypk -= h * INVMASSH;
                        h = zpxx * acor;
                        vzpi += h * invMassI;
                        vzpk -= h * INVMASSH;
                    }
                }

                // Third bond if present
                if (shakeID.w != -1) {
                    xpxx = xpi - xpl;
                    ypxx = ypi - ypl;
                    zpxx = zpi - zpl;
                    vxpxx = vxpi - vxpl;
                    vypxx = vypi - vypl;
                    vzpxx = vzpi - vzpl;
                    rvdot = xpxx * vxpxx + ypxx * vypxx + zpxx * vzpxx;

                    acor = -rvdot / (toler * (invMassI + INVMASSH));

                    if (abs(acor) >= tol2) {
                        done = false;

                        double h = xpxx * acor;
                        vxpi += h * invMassI;
                        vxpl -= h * INVMASSH;
                        h = ypxx * acor;
                        vypi += h * invMassI;
                        vypl -= h * INVMASSH;
                        h = zpxx * acor;
                        vzpi += h * invMassI;
                        vzpl -= h * INVMASSH;
                    }

                }

                if (done) {
                    break;
                }
            }


            if (done) {
                VELX(shakeID.x) = vxpi;
                VELY(shakeID.x) = vypi;
                VELZ(shakeID.x) = vzpi;
                VELX(shakeID.y) = vxpj;
                VELY(shakeID.y) = vypj;
                VELZ(shakeID.y) = vzpj;

                if (shakeID.z != -1) {
                    VELX(shakeID.z) = vxpk;
                    VELY(shakeID.z) = vypk;
                    VELZ(shakeID.z) = vzpk;
                }

                if (shakeID.w != -1) {
                    VELX(shakeID.w) = vxpl;
                    VELY(shakeID.w) = vypl;
                    VELZ(shakeID.w) = vzpl;
                }
            }
        }
    }
    else if (pos < cSim.fastShakeOffset) {
        pos -= cSim.shakeOffset;
        if (pos < cSim.fastShakeConstraints) {

            // Read atom data
#ifdef RATTLE_NEIGHBORLIST
            int4 shakeID = cSim.pImageFastShakeID[pos];
#else
            int4 shakeID = cSim.pFastShakeID[pos];
#endif

            //double2 shakeParm = cSim.pShakeParm[pos];
#ifdef RATTLE_HMR
            //double invMassH = cSim.pShakeInvMassH[pos];
#endif
            //double invMassI = shakeParm.x;

            //double wo = 1.0 / invMassI;
            //double wh = 1.0 / INVMASSH;
            double wo = cSim.wo;
            double wh = cSim.wh;
            double woh = wo + wh;
            double whh = wh + wh;
            double woh2 = woh * 2.0;
            double wowh2 = wo * whh;
            double whwh = wh * wh;
            double wohwoh = woh * woh;
            double wohwo = woh * wo;
            double whhwh = whh * wh;


            double xp1  = PATOMX(shakeID.x);
            double yp1  = PATOMY(shakeID.x);
            double zp1  = PATOMZ(shakeID.x);
            double vxp1 = VELX(shakeID.x);
            double vyp1 = VELY(shakeID.x);
            double vzp1 = VELZ(shakeID.x);
            double xp2  = PATOMX(shakeID.y);
            double yp2  = PATOMY(shakeID.y);
            double zp2  = PATOMZ(shakeID.y);
            double vxp2 = VELX(shakeID.y);
            double vyp2 = VELY(shakeID.y);
            double vzp2 = VELZ(shakeID.y);
            double xp3  = PATOMX(shakeID.z);
            double yp3  = PATOMY(shakeID.z);
            double zp3  = PATOMZ(shakeID.z);
            double vxp3 = VELX(shakeID.z);
            double vyp3 = VELY(shakeID.z);
            double vzp3 = VELZ(shakeID.z);

            // Step1 AB, VAB
            double xab = xp2 - xp1;
            double yab = yp2 - yp1;
            double zab = zp2 - zp1;
            double xbc = xp3 - xp2;
            double ybc = yp3 - yp2;
            double zbc = zp3 - zp2;
            double xca = xp1 - xp3;
            double yca = yp1 - yp3;
            double zca = zp1 - zp3;

            double xvab = vxp2 - vxp1;
            double yvab = vyp2 - vyp1;
            double zvab = vzp2 - vzp1;
            double xvbc = vxp3 - vxp2;
            double yvbc = vyp3 - vyp2;
            double zvbc = vzp3 - vzp2;
            double xvca = vxp1 - vxp3;
            double yvca = vyp1 - vyp3;
            double zvca = vzp1 - vzp3;

            // Step2 eab
            double ablng = sqrt(xab * xab + yab * yab + zab * zab);
            double bclng = sqrt(xbc * xbc + ybc * ybc + zbc * zbc);
            double calng = sqrt(xca * xca + yca * yca + zca * zca);

            double xeab = xab / ablng;
            double yeab = yab / ablng;
            double zeab = zab / ablng;
            double xebc = xbc / bclng;
            double yebc = ybc / bclng;
            double zebc = zbc / bclng;
            double xeca = xca / calng;
            double yeca = yca / calng;
            double zeca = zca / calng;

            // Step3 vabab
            double vabab = xvab * xeab + yvab * yeab + zvab * zeab;
            double vbcbc = xvbc * xebc + yvbc * yebc + zvbc * zebc;
            double vcaca = xvca * xeca + yvca * yeca + zvca * zeca;

            // Step4 tab
            double cosa = -xeab * xeca - yeab * yeca - zeab * zeca;
            double cosb = -xebc * xeab - yebc * yeab - zebc * zeab;
            double cosc = -xeca * xebc - yeca * yebc - zeca * zebc;
            double abmc = wh * cosa * cosb - woh * cosc;
            double bcma = wo * cosb * cosc - whh * cosa;
            double camb = wh * cosc * cosa - woh * cosb;

            double tabd = vabab * (woh2 - wo * cosc * cosc) 
                + vbcbc * camb + vcaca * bcma;
            double tbcd = vbcbc * (wohwoh - whwh * cosa * cosa) 
                + vcaca * abmc * wo + vabab * camb * wo;
            double tcad = vcaca * (woh2 - wo * cosb * cosb) 
                + vabab * bcma + vbcbc * abmc;
            double deno = 2.0 * wohwoh + wowh2 * cosa * cosb * cosc 
                - whhwh * cosa * cosa 
                - wohwo * (cosb * cosb + cosc * cosc);

            // Step5 v
            VELX(shakeID.x) = vxp1 + (xeab * tabd - xeca * tcad) * wh / deno;
            VELY(shakeID.x) = vyp1 + (yeab * tabd - yeca * tcad) * wh / deno;
            VELZ(shakeID.x) = vzp1 + (zeab * tabd - zeca * tcad) * wh / deno;
            VELX(shakeID.y) = vxp2 + (xebc * tbcd - xeab * tabd * wo) / deno;
            VELY(shakeID.y) = vyp2 + (yebc * tbcd - yeab * tabd * wo) / deno;
            VELZ(shakeID.y) = vzp2 + (zebc * tbcd - zeab * tabd * wo) / deno;
            VELX(shakeID.z) = vxp3 + (xeca * tcad * wo - xebc * tbcd) / deno;
            VELY(shakeID.z) = vyp3 + (yeca * tcad * wo - yebc * tbcd) / deno;
            VELZ(shakeID.z) = vzp3 + (zeca * tcad * wo - zebc * tbcd) / deno;
        }
    }
    else if (pos < cSim.slowShakeOffset) {
        pos -= cSim.fastShakeOffset;
        if (pos < cSim.slowShakeConstraints) {

            // Read SHAKE network data
#ifdef SHAKE_NEIGHBORLIST
            int  shakeID1 = cSim.pImageSlowShakeID1[pos];
            int4 shakeID2 = cSim.pImageSlowShakeID2[pos];
#else
            int  shakeID1 = cSim.pSlowShakeID1[pos];
            int4 shakeID2 = cSim.pSlowShakeID2[pos];
#endif
            double2 shakeParm = cSim.pSlowShakeParm[pos];
#ifdef RATTLE_HMR
            double invMassH = cSim.pSlowShakeInvMassH[pos];
#endif

            double  xpi = PATOMX(shakeID1);
            double  ypi = PATOMY(shakeID1);
            double  zpi = PATOMZ(shakeID1);
            double  xpj = PATOMX(shakeID2.x);
            double  ypj = PATOMY(shakeID2.x);
            double  zpj = PATOMZ(shakeID2.x);
            double  xpk = PATOMX(shakeID2.y);
            double  ypk = PATOMY(shakeID2.y);
            double  zpk = PATOMZ(shakeID2.y);
            double  xpl = PATOMX(shakeID2.z);
            double  ypl = PATOMY(shakeID2.z);
            double  zpl = PATOMZ(shakeID2.z);
            double  xpm = PATOMX(shakeID2.w);
            double  ypm = PATOMY(shakeID2.w);
            double  zpm = PATOMZ(shakeID2.w);
            double vxpi = VELX(shakeID1);
            double vypi = VELY(shakeID1);
            double vzpi = VELZ(shakeID1);
            double vxpj = VELX(shakeID2.x);
            double vypj = VELY(shakeID2.x);
            double vzpj = VELZ(shakeID2.x);
            double vxpk = VELX(shakeID2.y);
            double vypk = VELY(shakeID2.y);
            double vzpk = VELZ(shakeID2.y);
            double vxpl = VELX(shakeID2.z);
            double vypl = VELY(shakeID2.z);
            double vzpl = VELZ(shakeID2.z);
            double vxpm = VELX(shakeID2.w);
            double vypm = VELY(shakeID2.w);
            double vzpm = VELZ(shakeID2.w);

            double invMassI = shakeParm.x;
            double toler    = shakeParm.y;

            bool done = false;
            double tol2 = cSim.tol / dt;

            for (int i = 0; i < 3000; ++i) {
                done = true;

                //calculate rvdot;
                double  xpxx = xpi - xpj;
                double  ypxx = ypi - ypj;
                double  zpxx = zpi - zpj;
                double vxpxx = vxpi - vxpj;
                double vypxx = vypi - vypj;
                double vzpxx = vzpi - vzpj;
                double rvdot = xpxx * vxpxx + ypxx * vypxx + zpxx * vzpxx;

                // apply correction to first hydrogen
                double acor = -rvdot / (toler * (invMassI + INVMASSH));

                if (abs(acor) >= tol2) {
                    done = false;
                    double h = xpxx * acor;
                    vxpi += h * invMassI;
                    vxpj -= h * INVMASSH;
                    h = ypxx * acor;
                    vypi += h * invMassI;
                    vypj -= h * INVMASSH;
                    h = zpxx * acor;
                    vzpi += h * invMassI;
                    vzpj -= h * INVMASSH;
                }

                xpxx = xpi - xpk;
                ypxx = ypi - ypk;
                zpxx = zpi - zpk;
                vxpxx = vxpi - vxpk;
                vypxx = vypi - vypk;
                vzpxx = vzpi - vzpk;
                rvdot = xpxx * vxpxx + ypxx * vypxx + zpxx * vzpxx;

                // apply correction to second hydrogen
                acor = -rvdot / (toler * (invMassI + INVMASSH));

                if (abs(acor) >= tol2) {
                    done = false;

                    double h = xpxx * acor;
                    vxpi += h * invMassI;
                    vxpk -= h * INVMASSH;
                    h = ypxx * acor;
                    vypi += h * invMassI;
                    vypk -= h * INVMASSH;
                    h = zpxx * acor;
                    vzpi += h * invMassI;
                    vzpk -= h * INVMASSH;
                }

                xpxx = xpi - xpl;
                ypxx = ypi - ypl;
                zpxx = zpi - zpl;
                vxpxx = vxpi - vxpl;
                vypxx = vypi - vypl;
                vzpxx = vzpi - vzpl;
                rvdot = xpxx * vxpxx + ypxx * vypxx + zpxx * vzpxx;

                // apply correction to third hydrogen
                acor = -rvdot / (toler * (invMassI + INVMASSH));

                if (abs(acor) >= tol2) {
                    done = false;

                    double h = xpxx * acor;
                    vxpi += h * invMassI;
                    vxpl -= h * INVMASSH;
                    h = ypxx * acor;
                    vypi += h * invMassI;
                    vypl -= h * INVMASSH;
                    h = zpxx * acor;
                    vzpi += h * invMassI;
                    vzpl -= h * INVMASSH;
                }

                xpxx = xpi - xpm;
                ypxx = ypi - ypm;
                zpxx = zpi - zpm;
                vxpxx = vxpi - vxpm;
                vypxx = vypi - vypm;
                vzpxx = vzpi - vzpm;
                rvdot = xpxx * vxpxx + ypxx * vypxx + zpxx * vzpxx;

                // apply correction to fourth hydrogen
                acor = -rvdot / (toler * (invMassI + INVMASSH));

                if (abs(acor) >= tol2) {
                    done = false;

                    double h = xpxx * acor;
                    vxpi += h * invMassI;
                    vxpm -= h * INVMASSH;
                    h = ypxx * acor;
                    vypi += h * invMassI;
                    vypm -= h * INVMASSH;
                    h = zpxx * acor;
                    vzpi += h * invMassI;
                    vzpm -= h * INVMASSH;
                }

                if (done) {
                    break;
                }


            }

            if (done) {
                VELX(shakeID1) = vxpi;
                VELY(shakeID1) = vypi;
                VELZ(shakeID1) = vzpi;
                VELX(shakeID2.x) = vxpj;
                VELY(shakeID2.x) = vypj;
                VELZ(shakeID2.x) = vzpj;
                VELX(shakeID2.y) = vxpk;
                VELY(shakeID2.y) = vypk;
                VELZ(shakeID2.y) = vzpk;
                VELX(shakeID2.z) = vxpl;
                VELY(shakeID2.z) = vypl;
                VELZ(shakeID2.z) = vzpl;
                VELX(shakeID2.w) = vxpm;
                VELY(shakeID2.w) = vypm;
                VELZ(shakeID2.w) = vzpm;
            }
        }
    }

#undef PATOMX
#undef PATOMY
#undef PATOMZ
#undef VELX
#undef VELY
#undef VELZ
#undef INVMASSH
}


