#ifdef BUILD_PAIRS_CALC_EFV
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_efv_2cut
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file.
!*******************************************************************************

subroutine pairs_calc_efv_2cut(img_frc, img_crd, img_qterm, ef_tbl, eed_cub, &
                               ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)

#else

!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_efv
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file.
!*******************************************************************************

subroutine pairs_calc_efv(img_frc, img_crd, img_qterm, ef_tbl, eed_cub, &
                          ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)

#endif /* BUILD_PAIRS_CALC_NOVEC_2CUT */
#endif /* BUILD_PAIRS_CALC_EFV */

#ifdef BUILD_PAIRS_CALC_FV
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_fv_2cut
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file.
!*******************************************************************************

subroutine pairs_calc_fv_2cut(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
                              ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)

#else

!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_fv
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file.
!*******************************************************************************

subroutine pairs_calc_fv(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
                         ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)

#endif /* BUILD_PAIRS_CALC_NOVEC_2CUT */
#endif /* BUILD_PAIRS_CALC_FV */

#ifdef BUILD_PAIRS_CALC_F
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_f_2cut
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file.
!*******************************************************************************

subroutine pairs_calc_f_2cut(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
                             ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)

#else

!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_f
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file.
!*******************************************************************************

subroutine pairs_calc_f(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
                        ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)

#endif /* BUILD_PAIRS_CALC_NOVEC_2CUT */
#endif /* BUILD_PAIRS_CALC_F */

  implicit none

! Formal arguments:

  double precision              :: img_frc(3, *)
  double precision, intent(in)  :: img_crd(3, *)
  double precision, intent(in)  :: img_qterm(*)
#ifdef NEED_ENE
  double precision, intent(in)  :: ef_tbl(*)
#else
  double precision, intent(in)  :: f_tbl(*)
#endif
  double precision, intent(in)  :: eed_cub(*)
  integer, intent(in)           :: ico(*)
  integer                       :: ipairs_sublst(*)
  integer, intent(in)           :: img_iac(*)
  double precision, intent(in)  :: cn1(*), cn2(*)
  double precision, intent(in)  :: x_tran(1:3, 0:17)

! Local variables:

  integer, parameter            :: mask27 = int(Z"07FFFFFF")
  double precision, parameter   :: half = 1.d0/2.d0
  double precision, parameter   :: third = 1.d0/3.d0

  double precision      :: dumx, dumy, dumz

  double precision      :: cgi
  double precision      :: cgi_cgj
  double precision      :: b0, b1
  double precision      :: df
  double precision      :: dfx, dfy, dfz
  double precision      :: f6, r6, f12
#ifdef HAS_10_12
  double precision      :: f10, r10
#endif
  double precision      :: du, du2, du3 ! 'u' is another name for delr2
  double precision      :: del_efs
  double precision      :: dens_efs
  double precision      :: lowest_efs_u

  ! Variables used with erfc switch table; name are historical:

  double precision      :: switch
  double precision      :: d_switch_dx
  double precision      :: x, dx, e3dx, e4dx2
  double precision      :: delr, delrinv

  integer               :: iaci
  integer               :: ic
  integer               :: ind
  integer               :: nxt_img_j, img_j
  integer               :: itran
  integer               :: sublst_idx
  integer               :: saved_pairlist_val

  double precision      :: nxt_delx, nxt_dely, nxt_delz
  double precision      :: delx, dely, delz, delr2, delr2inv

  !IPS variables
  double precision      ::  uips,uips2,uips4,uips2r,uips6r,uips12r
  double precision      ::  pipse,eipse,dpipse,pvc,dvcu,pva,dvau
  double precision      ::  ecur
  integer               ::  icount

  dens_efs = efs_tbl_dens
  del_efs = 1.d0 / dens_efs
  lowest_efs_u = lowest_efs_delr2

! First loop over the ee evaluation-only pairs:

  dumx = 0.d0
  dumy = 0.d0
  dumz = 0.d0

  cgi = img_qterm(img_i)

  ! The pairlist must have one dummy end entry to cover reading past the
  ! end of the list...

  saved_pairlist_val = ipairs_sublst(ee_eval_cnt + full_eval_cnt + 1)

  ipairs_sublst(ee_eval_cnt + full_eval_cnt + 1) = &
    ipairs_sublst(ee_eval_cnt + full_eval_cnt)

  if (common_tran .eq. 1) then
    nxt_img_j = ipairs_sublst(1)
    itran = 13
  else
    nxt_img_j = iand(ipairs_sublst(1), mask27)
    itran = ishft(ipairs_sublst(1), -27)
  end if

  nxt_delx = img_crd(1, nxt_img_j) + x_tran(1, itran)
  nxt_dely = img_crd(2, nxt_img_j) + x_tran(2, itran)
  nxt_delz = img_crd(3, nxt_img_j) + x_tran(3, itran)

  icount = 0

  do sublst_idx = 2, ee_eval_cnt + 1

    img_j = nxt_img_j
    delx = nxt_delx
    dely = nxt_dely
    delz = nxt_delz

    if (common_tran .eq. 1) then
      nxt_img_j = ipairs_sublst(sublst_idx)
      itran = 13
    else
      nxt_img_j = iand(ipairs_sublst(sublst_idx), mask27)
      itran = ishft(ipairs_sublst(sublst_idx), -27)
    end if

    nxt_delx = img_crd(1, nxt_img_j) + x_tran(1, itran)
    nxt_dely = img_crd(2, nxt_img_j) + x_tran(2, itran)
    nxt_delz = img_crd(3, nxt_img_j) + x_tran(3, itran)

    delr2 = delx * delx + dely * dely + delz * delz

#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
    if (delr2 .lt. es_cut2) then
#else
    if (delr2 .lt. max_nb_cut2) then
#endif

      cgi_cgj = cgi * img_qterm(img_j)

        delr = sqrt(delr2)
        delrinv = 1.d0 / delr
!         delrinv = cache_df(im_new)
         delr2inv = delrinv*delrinv
!ifdef LES no
!         b0 = cgi*charge(j)*ripsr
         b0 = cgi_cgj*ripsr

         !  -- use the ips potential:
         !          
         uips=ripsr*delr2*delrinv
         uips2=rips2r*delr2
         pipse=1.d0/uips+AIPSE(0)+uips2*(AIPSE(1)+uips2*(AIPSE(2)+uips2*AIPSE(3)))
         dpipse=-1.d0/uips+uips2*(BIPSE(1)+uips2*(BIPSE(2)+uips2*BIPSE(3)))
         ecur=b0*(pipse-pipsec)
#ifdef NEED_ENE
!        eelt = eelt + ecur
         eed_stk = eed_stk + ecur
#endif /* NEED_ENE */
         df=-b0*dpipse*delr2inv

      dfx = delx * df
      dfy = dely * df
      dfz = delz * df

#ifdef NEED_VIR
      vxx = vxx - delx * dfx
      vxy = vxy - delx * dfy
      vxz = vxz - delx * dfz
      vyy = vyy - dely * dfy
      vyz = vyz - dely * dfz
      vzz = vzz - delz * dfz
#endif /* NEED_VIR */

      dumx = dumx + dfx
      dumy = dumy + dfy
      dumz = dumz + dfz

      img_frc(1, img_j) = img_frc(1, img_j) + dfx
      img_frc(2, img_j) = img_frc(2, img_j) + dfy
      img_frc(3, img_j) = img_frc(3, img_j) + dfz

    end if

  end do

  iaci = ntypes_stk * (img_iac(img_i) - 1)

  do sublst_idx = ee_eval_cnt + 2, ee_eval_cnt + full_eval_cnt + 1

    img_j = nxt_img_j
    delx = nxt_delx
    dely = nxt_dely
    delz = nxt_delz

    if (common_tran .eq. 1) then
      nxt_img_j = ipairs_sublst(sublst_idx)
      itran = 13
    else
      nxt_img_j = iand(ipairs_sublst(sublst_idx), mask27)
      itran = ishft(ipairs_sublst(sublst_idx), -27)
    end if

    nxt_delx = img_crd(1, nxt_img_j) + x_tran(1, itran)
    nxt_dely = img_crd(2, nxt_img_j) + x_tran(2, itran)
    nxt_delz = img_crd(3, nxt_img_j) + x_tran(3, itran)

    delr2 = delx * delx + dely * dely + delz * delz

    if (delr2 .lt. max_nb_cut2) then

      icount = icount + 1
      ic = ico(iaci + img_iac(img_j))

#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
      if (delr2 .lt. es_cut2) then
#endif

      cgi_cgj = cgi * img_qterm(img_j)

        delr = sqrt(delr2)
        delrinv = 1.d0 / delr
!         delrinv = cache_df(im_new)
         delr2inv = delrinv*delrinv
!ifdef LES no
!         b0 = cgi*charge(j)*ripsr
         b0 = cgi_cgj*ripsr

         !  -- use the ips potential:
         !          
         uips=ripsr*delr2*delrinv
         uips2=rips2r*delr2
         pipse=1.d0/uips+AIPSE(0)+uips2*(AIPSE(1)+uips2*(AIPSE(2)+uips2*AIPSE(3)))
         dpipse=-1.d0/uips+uips2*(BIPSE(1)+uips2*(BIPSE(2)+uips2*BIPSE(3)))
         ecur=b0*(pipse-pipsec)
#ifdef NEED_ENE
!        eelt = eelt + ecur
         eed_stk = eed_stk + ecur
#endif /* NEED_ENE */
         df=-b0*dpipse*delr2inv

#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
      else
        delr2inv = 1.d0 / delr2
        df = 0.d0
      end if
#endif

#ifdef HAS_10_12
      if (ic .gt. 0) then
#endif
   f6 = cn2(ic)*rips6r
   f12 = cn1(ic)*rips12r
!  L-J r6 term
!   etr6=1/r6+a0+r2*(a1+r2*(a2+a3*r2))
!   detr6/dr*r1=-6/r6+r2*(d1+r2*(d2+d3*r2))
!
        UIPS2R=(DELR2INV*RIPS2)
        UIPS2=1.d0/UIPS2R
        UIPS4=UIPS2*UIPS2
        UIPS6R=UIPS2R*UIPS2R*UIPS2R
        UIPS12R=UIPS6R*UIPS6R
            PVC=UIPS6R+AIPSVC(0)+UIPS2*(AIPSVC(1)+UIPS2*(AIPSVC(2)+UIPS2*AIPSVC(3)))
            DVCU=-6.d0*UIPS6R+UIPS2*(BIPSVC(1)+UIPS2*(BIPSVC(2)+UIPS2*BIPSVC(3)))
!  L-J r12 term 
!   etr12=1/r12+a0+r2*(a1+r2*(a2+a3*r2))
!   detr12/dr*r1=-12/r12+r4*(d1+r4*(d2+d3*r4))
!
            PVA=UIPS12R+AIPSVA(0)+UIPS4*(AIPSVA(1)+UIPS4*(AIPSVA(2)+UIPS4*AIPSVA(3)))
            DVAU=-12.d0*UIPS12R+UIPS4*(BIPSVA(1)+UIPS4*(BIPSVA(2)+UIPS4*BIPSVA(3) ))
#ifdef NEED_ENE
            evdw_stk = evdw_stk +F12*(PVA-PIPSVAC)-F6*(PVC-PIPSVCC)
#endif /* NEED_ENE */
            df = df -(F12*DVAU-F6*DVCU)*DELR2INV


#ifdef HAS_10_12
      else if (ic.lt.0) then
        ! This code allows 10-12 terms; in many (most?) (all?) cases, the
        ! only "nominal" 10-12 terms are on waters, where the asol and bsol
        ! parameters are always zero; hence we can skip this part.
        ic = - ic
        r10 = delr2inv * delr2inv * delr2inv * delr2inv * delr2inv
        f10 = gbl_bsol(ic) * r10
        f12 = gbl_asol(ic) * (r10 * delr2inv)
#ifdef NEED_ENE
        ehb_stk = ehb_stk + f12 - f10
#endif /* NEED_ENE */
        df = df + (12.d0 * f12 - 10.d0 * f10) * delr2inv
      end if
#endif /* HAS_10_12 */

      dfx = delx * df
      dfy = dely * df
      dfz = delz * df

#ifdef NEED_VIR
      vxx = vxx - delx * dfx
      vxy = vxy - delx * dfy
      vxz = vxz - delx * dfz
      vyy = vyy - dely * dfy
      vyz = vyz - dely * dfz
      vzz = vzz - delz * dfz
#endif /* NEED_VIR */

      dumx = dumx + dfx
      dumy = dumy + dfy
      dumz = dumz + dfz

      img_frc(1, img_j) = img_frc(1, img_j) + dfx
      img_frc(2, img_j) = img_frc(2, img_j) + dfy
      img_frc(3, img_j) = img_frc(3, img_j) + dfz

    end if

  end do

  img_frc(1, img_i) = img_frc(1, img_i) - dumx
  img_frc(2, img_i) = img_frc(2, img_i) - dumy
  img_frc(3, img_i) = img_frc(3, img_i) - dumz

  ipairs_sublst(ee_eval_cnt + full_eval_cnt + 1) = saved_pairlist_val

  nnbips=nnbips+icount*2

  return

#ifdef BUILD_PAIRS_CALC_EFV
#undef NEED_ENE
#undef NEED_VIR
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
end subroutine pairs_calc_efv_2cut
#else
end subroutine pairs_calc_efv
#endif
#endif /* BUILD_PAIRS_CALC_EFV */

#ifdef BUILD_PAIRS_CALC_FV
#undef NEED_VIR
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
end subroutine pairs_calc_fv_2cut
#else
end subroutine pairs_calc_fv
#endif
#endif /* BUILD_PAIRS_CALC_FV */

#ifdef BUILD_PAIRS_CALC_F
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
end subroutine pairs_calc_f_2cut
#else
end subroutine pairs_calc_f
#endif
#endif /* BUILD_PAIRS_CALC_F */

