#ifdef BUILD_PAIRS_CALC_EFV
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
#ifdef FSWITCH
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_efv_2cut_fs
!
! Description:  Direct force computation on one atom using fswitch.  
! Different implementations are covered by conditional compilation in this file.
!*******************************************************************************
subroutine pairs_calc_efv_2cut_fs(img_frc, img_crd, img_qterm, ef_tbl, eed_cub, &
                               ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)

#else
#ifdef VACUUM
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_efv_2cut_vac
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file.
!*******************************************************************************
subroutine pairs_calc_efv_2cut_vac(img_frc, img_crd, img_qterm, ef_tbl, eed_cub, &
                               ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)
#else
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_efv_2cut
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file.
!*******************************************************************************
subroutine pairs_calc_efv_2cut(img_frc, img_crd, img_qterm, ef_tbl, eed_cub, &
                               ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)
#endif /*VACUUM*/
#endif /*FSWITCH*/
#else
#ifdef FSWITCH
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_efv_fs
!
! Description:  Direct force computation on one atom using fswitch. 
! Different implementations are covered by conditional compilation in this file.
!*******************************************************************************

subroutine pairs_calc_efv_fs(img_frc, img_crd, img_qterm, ef_tbl, eed_cub, &
                          ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)

#else
#ifdef VACUUM
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_efv_vac
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file.
!*******************************************************************************

subroutine pairs_calc_efv_vac(img_frc, img_crd, img_qterm, ef_tbl, eed_cub, &
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
#endif /*VACUUM*/
#endif /*FSWITCH*/
#endif /* BUILD_PAIRS_CALC_NOVEC_2CUT */
#endif /* BUILD_PAIRS_CALC_EFV */

#ifdef BUILD_PAIRS_CALC_FV
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
#ifdef FSWITCH
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_fv_2cut_fs
!
! Description:  Direct force computation on one atom using fswitch. 
! Different implementations are covered by conditional compilation in this file.
!*******************************************************************************

subroutine pairs_calc_fv_2cut_fs(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
                              ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)

#else
#ifdef VACUUM
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_fv_2cut_vac
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file.
!*******************************************************************************
subroutine pairs_calc_fv_2cut_vac(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
                              ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)
#else
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_fv_2cut
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file.
!*******************************************************************************
subroutine pairs_calc_fv_2cut(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
                              ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)
#endif /*VACUUM*/
#endif /*FSWITCH*/
#else
#ifdef FSWITCH
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_fv_fs
!
! Description:  Direct force computation on one atom using fswitch. 
! Different implementations are covered by conditional compilation in this file.
!*******************************************************************************

subroutine pairs_calc_fv_fs(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
                         ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)

#else
#ifdef VACUUM
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_fv_vac
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file.
!*******************************************************************************
subroutine pairs_calc_fv_vac(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
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

#endif /* VACUUM */
#endif /* FSWITCH */
#endif /* BUILD_PAIRS_CALC_NOVEC_2CUT */
#endif /* BUILD_PAIRS_CALC_FV */

#ifdef BUILD_PAIRS_CALC_F
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
#ifdef FSWITCH
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_f_2cut_fs
!
! Description:  Direct force computation on one atom using fswitch. 
! Different implementations are covered by conditional compilation in this file.
!*******************************************************************************

subroutine pairs_calc_f_2cut_fs(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
                             ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)
#else
#ifdef VACUUM
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_f_2cut_vac
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file.
!*******************************************************************************

subroutine pairs_calc_f_2cut_vac(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
                             ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)
#else
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_f_2cut
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file.
!*******************************************************************************

subroutine pairs_calc_f_2cut(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
                             ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)
#endif /* VACUUM */
#endif /* FSWITCH */
#else
#ifdef FSWITCH
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_f_fs
!
! Description:  Direct force computation on one atom using fswitch. 
! Different implementations are covered by conditional compilation in this file.
!*******************************************************************************

subroutine pairs_calc_f_fs(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
                        ico, ipairs_sublst, img_iac, cn1, cn2, x_tran)

#else
#ifdef VACUUM
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_f_vac
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file.
!*******************************************************************************

subroutine pairs_calc_f_vac(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
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

#endif /* VACUUM */
#endif /* FSWITCH */
#endif /* BUILD_PAIRS_CALC_NOVEC_2CUT */
#endif /* BUILD_PAIRS_CALC_F */
  use mdin_ctrl_dat_mod, only : fswitch
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
#ifdef NEED_VIR
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

#ifdef FSWITCH
  double precision      :: delr3inv, delr12inv
  double precision      :: cut, cut6, cut3, cutinv, cut2inv, cut3inv, cut6inv
  double precision      :: fswitch2, fswitch3, fswitch6
  double precision      :: p12, p6, a_energy, b_energy, df12, df6

  double precision      :: invfswitch6cut6, invfswitch3cut3
#endif /*FSWITCH*/

  integer               :: iaci
  integer               :: ic
  integer               :: ind
  integer               :: nxt_img_j, img_j
  integer               :: itran
  integer               :: sublst_idx
  integer               :: saved_pairlist_val

  double precision      :: nxt_delx, nxt_dely, nxt_delz
  double precision      :: delx, dely, delz, delr2, delr2inv
  dens_efs = efs_tbl_dens
  del_efs = 1.d0 / dens_efs
  lowest_efs_u = lowest_efs_delr2

! Precompute some constants

#ifdef FSWITCH
  fswitch2 = fswitch * fswitch
  fswitch3 = fswitch2 * fswitch
  fswitch6 = fswitch3 * fswitch3
  cut = sqrt(max_nb_cut2)
  cut3 = cut*cut*cut
  cut6 = cut3*cut3
  cutinv = 1/cut
  cut2inv = cutinv*cutinv
  cut3inv = cut2inv*cutinv
  cut6inv = cut3inv*cut3inv

  invfswitch6cut6 = cut6inv / fswitch6
  invfswitch3cut3 = cut3inv / fswitch3
#endif /*FSWITCH*/

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
      if (delr2 .ge. lowest_efs_u) then

        ! Do the Coulomb part of the direct sum using efs: 

        ind = int(dens_efs * delr2)
        du = delr2 - dble(ind) * del_efs
        du2 = du * du
        du3 = du * du2
#ifdef NEED_ENE
        ind = ishft(ind, 3)             ! 8 * ind

        b0 = cgi_cgj * (ef_tbl(1 + ind) + du * ef_tbl(2 + ind) + &
             du2 * ef_tbl(3 + ind) + du3 * ef_tbl(4 + ind))

        df = cgi_cgj * (ef_tbl(5 + ind) + du * ef_tbl(6 + ind) + &
             du2 * ef_tbl(7 + ind) + du3 * ef_tbl(8 + ind))
#else
        ind = ishft(ind, 2)             ! 4 * ind

        df = cgi_cgj * (f_tbl(1 + ind) + du * f_tbl(2 + ind) + &
             du2 * f_tbl(3 + ind) + du3 * f_tbl(4 + ind))
#endif /* NEED_ENE */

#ifdef NEED_ENE
        eed_stk = eed_stk + b0
#endif /* NEED_ENE */

#ifdef NEED_VIR
        eedvir_stk = eedvir_stk - df * delr2 
#endif /* NEED_VIR */
      else
        ! Do the Coulomb part of the direct sum using erfc spline table: 

        delr = sqrt(delr2)
        delrinv = 1.0 / delr

        x = dxdr * delr
        ind = int(eedtbdns_stk * x)
        dx = x - dble(ind) * del
        ind = ishft(ind, 2)             ! 4 * ind

        e3dx  = dx * eed_cub(3 + ind)
        e4dx2 = dx * dx *  eed_cub(4 + ind)

#ifdef VACUUM
        switch = 1.0
        d_switch_dx = 1.0
#else
        switch = eed_cub(1 + ind) + &
                 dx * (eed_cub(2 + ind) + half * (e3dx + third * e4dx2))

        d_switch_dx = eed_cub(2 + ind) + e3dx + half * e4dx2
#endif
        b0 = cgi_cgj * delrinv * switch
        b1 = (b0 - cgi_cgj * d_switch_dx * dxdr)

#ifdef NEED_VIR
        eedvir_stk = eedvir_stk - b1
#endif /* NEED_VIR */
#ifdef NEED_ENE
        eed_stk = eed_stk + b0
#endif /* NEED_ENE */

        df = b1 * delrinv * delrinv

      end if

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

      ic = ico(iaci + img_iac(img_j))

#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
      if (delr2 .lt. es_cut2) then
#endif

      cgi_cgj = cgi * img_qterm(img_j)
!write(0,*)"Test",cgi_cgj,cgi,delr2
      if (delr2 .ge. lowest_efs_u) then

        ! Do the Coulomb part of the direct sum using efs: 

        ind = int(dens_efs * delr2)
        du = delr2 - dble(ind) * del_efs
        du2 = du * du
        du3 = du * du2

#ifdef NEED_ENE
        ind = ishft(ind, 3)             ! 8 * ind

        b0 = cgi_cgj * (ef_tbl(1 + ind) + du * ef_tbl(2 + ind) + &
             du2 * ef_tbl(3 + ind) + du3 * ef_tbl(4 + ind))

        df = cgi_cgj * (ef_tbl(5 + ind) + du * ef_tbl(6 + ind) + &
             du2 * ef_tbl(7 + ind) + du3 * ef_tbl(8 + ind))
#else
        ind = ishft(ind, 2)             ! 4 * ind

        df = cgi_cgj * (f_tbl(1 + ind) + du * f_tbl(2 + ind) + &
             du2 * f_tbl(3 + ind) + du3 * f_tbl(4 + ind))
#endif /* NEED_ENE */
#ifdef NEED_ENE
        eed_stk = eed_stk + b0
#endif /* NEED_ENE */

#ifdef NEED_VIR
        eedvir_stk = eedvir_stk - df * delr2
#endif /* NEED_VIR */
        delr2inv = 1.0 / delr2
      else    
        ! Do the Coulomb part of the direct sum using erfc spline table: 

        delr = sqrt(delr2)
        delrinv = 1.0 / delr

        x = dxdr * delr
        ind = int(eedtbdns_stk * x)
        dx = x - dble(ind) * del
        ind = ishft(ind, 2)             ! 4 * ind

        e3dx  = dx * eed_cub(3 + ind)
        e4dx2 = dx * dx *  eed_cub(4 + ind)
#ifdef VACUUM
        switch = 1.0
        d_switch_dx = 1.0
#else
        switch = eed_cub(1 + ind) + &
                 dx * (eed_cub(2 + ind) + half * (e3dx + third * e4dx2))

        d_switch_dx = eed_cub(2 + ind) + e3dx + half * e4dx2
#endif

        b0 = cgi_cgj * delrinv * switch
!write(0,*)b0,cgi_cgj,delrinv,switch
        b1 = (b0 - cgi_cgj * d_switch_dx * dxdr)

#ifdef NEED_ENE
        eed_stk = eed_stk + b0
#endif /* NEED_ENE */
#ifdef NEED_VIR
        eedvir_stk = eedvir_stk - b1
#endif /* NEED_VIR */

        delr2inv = delrinv * delrinv
        df = b1 * delrinv * delrinv

      end if

#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
      else
        delr2inv = 1.0 / delr2
        df = 0.0
      end if
#endif

#ifdef HAS_10_12
      if (ic .gt. 0) then
#endif
        r6 = delr2inv * delr2inv * delr2inv
#ifdef FSWITCH
          delr2inv=1.d0 / delr2
          delr3inv=delr2inv * sqrt(delr2inv)
          delr12inv=r6*r6
            
          if(delr2 .gt. fswitch2) then ! r < ron
            p12=cut6/(cut6-fswitch6)
            p6= cut3/(cut3-fswitch3)
            f12=cn1(ic)*p12*(r6-cut6inv)*(r6-cut6inv)
            f6=cn2(ic)*p6*(delr3inv-cut3inv)*(delr3inv-cut3inv)
            df12=-12.d0*p12*delr2inv*r6*(r6-cut6inv)
            df6=-6.d0*p6*(delr3inv-cut3inv)*delr3inv*delr2inv
          else  ! r<fswitch
            f12=cn1(ic)*delr12inv-cn1(ic)*invfswitch6cut6
            f6=cn2(ic)*r6-cn2(ic)*invfswitch3cut3
            df12=-12.d0*delr2inv*delr12inv
            df6=-6.d0*delr2inv*r6
          end if
#ifdef NEED_ENE
          evdw_stk = evdw_stk + f12 - f6
#endif /* NEED ENE */
          df = df - cn1(ic)*df12 + cn2(ic)*df6
#else
          f6 = cn2(ic) * r6
          f12 = cn1(ic) * (r6 * r6)
#ifdef NEED_ENE
          evdw_stk = evdw_stk + f12 - f6
#endif /* NEED_ENE */
          df = df + (12.0 * f12 - 6.0 * f6) * delr2inv
#endif /* FSWITCH */

#ifdef HAS_10_12
      else if( ic.lt.0 ) then
        ! This code allows 10-12 terms; in many (most?) (all?) cases, the
        ! only "nominal" 10-12 terms are on waters, where the asol and bsol
        ! parameters are always zero; hence we can skip this part.
        ic = - ic
        r10 = delr2inv * delr2inv * delr2inv * delr2inv * delr2inv
        f10 = gbl_bsol(ic) * r10
        f12 = gbl_asol(ic) * (r10 * delr2inv)
#ifdef NEED_ENE
        ehb_stk = ehb_stk + f12 - f10
      end if
#endif /* NEED_ENE */
        df = df + (12.0 * f12 - 10.0 * f10) * delr2inv
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
  return

#ifdef BUILD_PAIRS_CALC_EFV
#undef NEED_ENE
#undef NEED_VIR
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
#ifdef FSWITCH
end subroutine pairs_calc_efv_2cut_fs
#else
#ifdef VACUUM
end subroutine pairs_calc_efv_2cut_vac
#else
end subroutine pairs_calc_efv_2cut
#endif /* VACUUM */
#endif /* FSWITCH */
#else
#ifdef FSWITCH
end subroutine pairs_calc_efv_fs
#else
#ifdef VACUUM
end subroutine pairs_calc_efv_vac
#else
end subroutine pairs_calc_efv
#endif /* VACUUM */
#endif /* FSWITCH */
#endif
#endif /* BUILD_PAIRS_CALC_EFV */

#ifdef BUILD_PAIRS_CALC_FV
#undef NEED_VIR
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
#ifdef FSWITCH
end subroutine pairs_calc_fv_2cut_fs
#else
#ifdef VACUUM
end subroutine pairs_calc_fv_2cut_vac
#else
end subroutine pairs_calc_fv_2cut
#endif /* VACUUM */
#endif /* FSWITCH */
#else
#ifdef FSWITCH
end subroutine pairs_calc_fv_fs
#else
#ifdef VACUUM
end subroutine pairs_calc_fv_vac
#else
end subroutine pairs_calc_fv
#endif /*VACUUM*/
#endif /* FSWITCH */
#endif
#endif /* BUILD_PAIRS_CALC_FV */

#ifdef BUILD_PAIRS_CALC_F
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
#ifdef FSWITCH
end subroutine pairs_calc_f_2cut_fs
#else
#ifdef VACUUM
end subroutine pairs_calc_f_2cut_vac
#else
end subroutine pairs_calc_f_2cut
#endif /* VACUUM */
#endif /* FSWITCH */
#else
#ifdef FSWITCH
end subroutine pairs_calc_f_fs
#else
#ifdef VACUUM
end subroutine pairs_calc_f_vac
#else
end subroutine pairs_calc_f
#endif /* VACUUM */
#endif /* FSWITCH */
#endif
#endif /* BUILD_PAIRS_CALC_F */


