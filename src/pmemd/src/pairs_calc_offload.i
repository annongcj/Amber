#ifdef MIC_offload

#ifdef BUILD_PAIRS_CALC_EFV
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_efv_2cut_offload
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file.
!*******************************************************************************
!DIR$ ATTRIBUTES OFFLOAD : MIC :: pairs_calc_efv_2cut_offload
subroutine pairs_calc_efv_2cut_offload(img_frc, img_crd, img_qterm, ef_tbl, eed_cub, &
                          ico, ipairs_sublst, img_iac, cn1, cn2, x_tran, img_i,ee_eval_cnt, full_eval_cnt,common_tran,evdw_stk,eed_stk, eedvir_stk, ehb_stk, vxx, vyy, vzz,img_frc_thread,my_thread_id, atm_cnt, es_cut2, max_nb_cut2,eedtbdns_stk)

#else
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_efv_offload
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file, but all
!               versions here should only be invoked if DIRFRC_NOVEC is
!               defined.
!*******************************************************************************

!DIR$ ATTRIBUTES OFFLOAD : MIC :: pairs_calc_efv_offload
subroutine pairs_calc_efv_offload(img_frc, img_crd, img_qterm, ef_tbl, eed_cub, &
                          ico, ipairs_sublst, img_iac, cn1, cn2, x_tran, img_i,ee_eval_cnt, full_eval_cnt,common_tran,evdw_stk,eed_stk, eedvir_stk, ehb_stk, vxx, vyy, vzz,img_frc_thread,my_thread_id, atm_cnt, es_cut2, max_nb_cut2,eedtbdns_stk)

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

!DIR$ ATTRIBUTES OFFLOAD : MIC :: pairs_calc_fv_2cut_offload
subroutine pairs_calc_fv_2cut_offload(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
                         ico, ipairs_sublst, img_iac, cn1, cn2, x_tran, img_i,ee_eval_cnt, full_eval_cnt,common_tran, eedvir_stk,vxx, vyy, vzz,img_frc_thread,my_thread_id, atm_cnt, es_cut2, max_nb_cut2,eedtbdns_stk)
#else
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_fv
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file, but all
!               versions here should only be invoked if DIRFRC_NOVEC is
!               defined.
!*******************************************************************************

!DIR$ ATTRIBUTES OFFLOAD : MIC :: pairs_calc_fv_offload
subroutine pairs_calc_fv_offload(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
                         ico, ipairs_sublst, img_iac, cn1, cn2, x_tran, img_i,ee_eval_cnt, full_eval_cnt,common_tran,eedvir_stk, vxx, vyy, vzz,img_frc_thread,my_thread_id, atm_cnt, es_cut2, max_nb_cut2,eedtbdns_stk)
#endif
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
!DIR$ ATTRIBUTES OFFLOAD : MIC :: pairs_calc_f_2cut_offload
subroutine pairs_calc_f_2cut_offload(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
                        ico, ipairs_sublst, img_iac, cn1, cn2, x_tran, img_i,ee_eval_cnt, full_eval_cnt,common_tran, eedvir_stk, vxx, vyy, vzz,img_frc_thread,my_thread_id, atm_cnt, es_cut2, max_nb_cut2,eedtbdns_stk)
#else

!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_f
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file, but all
!               versions here should only be invoked if DIRFRC_NOVEC is
!               defined.
!*******************************************************************************

!DIR$ ATTRIBUTES OFFLOAD : MIC :: pairs_calc_f_offload
subroutine pairs_calc_f_offload(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
                        ico, ipairs_sublst, img_iac, cn1, cn2, x_tran, img_i,ee_eval_cnt, full_eval_cnt,common_tran, eedvir_stk, vxx, vyy, vzz,img_frc_thread,my_thread_id, atm_cnt, es_cut2, max_nb_cut2,eedtbdns_stk)
#endif
#endif /* BUILD_PAIRS_CALC_F */


  implicit none

! Formal arguments:

  double precision              :: img_frc(3, *)
/*openmp*/
  double precision              :: img_frc_thread(3, *)
/*done*/
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
  double precision, intent(in) :: es_cut2, max_nb_cut2,eedtbdns_stk

/*Ashraf: for Openmp the following variables are also passed*/
  integer, intent(in)           :: img_i, ee_eval_cnt, full_eval_cnt, common_tran,my_thread_id, atm_cnt
#ifdef NEED_ENE
    double precision 		:: evdw_stk,eed_stk,ehb_stk
#endif
    double precision 		:: eedvir_stk, vxx, vyy, vzz
/* done passing */
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
  double precision      :: eedvir_simd, eed_simd, evdw_simd, dumx_simd, dumy_simd, dumz_simd
  double precision      :: vxx_simd, vyy_simd, vzz_simd
  integer		:: offset
#ifdef Faster_erfc
  ! for error function
  double precision :: a1=0.254829592d0, a2=-0.284496736d0, a3=1.421413741d0, a4=-1.453152027d0, a5=1.061405429d0, p=0.3275911d0
  double precision	:: t, exp_term
#endif /*Faster_erfc*/
  double precision	:: neg2InvSqrtPi
  double precision,parameter    :: PI = 3.141592653589793d0
  neg2InvSqrtPi = -2.0 / sqrt(PI) 

  offset = my_thread_id*atm_cnt 
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


!following are simd related
	dumx_simd = 0.0
	dumy_simd = 0.0
	dumz_simd = 0.0
#ifdef NEED_VIR
	eedvir_simd = eedvir_stk
	vxx_simd = vxx
!	vxy_simd = vxy ! these values are not used later, so removed it
!	vxz_simd = vxz
	vyy_simd = vyy
!	vyz_simd = vyz
	vzz_simd = vzz
#endif
#ifdef NEED_ENE
	eed_simd = eed_stk
#endif
!$omp simd private(img_j, itran, delx, dely, delz, delr2, cgi_cgj, ind, du, du2, du3, b0, df, dfx, dfy, dfz, delr,delrinv, x, dx, e3dx, e4dx2, switch, d_switch_dx, b1 ) & 
#ifdef NEED_VIR
!$omp& Reduction(-:eedvir_simd) &
!$omp& Reduction(+:vxx_simd, vyy_simd,vzz_simd) &
#endif
#ifdef NEED_ENE
!$omp& Reduction(+:eed_simd) &
#endif
!$omp& Reduction(+:dumx_simd, dumy_simd, dumz_simd) 
  !do sublst_idx = 2, ee_eval_cnt + 1
  do sublst_idx = 1, ee_eval_cnt


!    img_j = nxt_img_j
!    delx = nxt_delx
!    dely = nxt_dely
!    delz = nxt_delz


    if (common_tran .eq. 1) then
      !nxt_img_j = ipairs_sublst(sublst_idx)
      img_j = ipairs_sublst(sublst_idx)
      itran = 13
    else
      !nxt_img_j = iand(ipairs_sublst(sublst_idx), mask27)
      img_j = iand(ipairs_sublst(sublst_idx), mask27)
      itran = ishft(ipairs_sublst(sublst_idx), -27)
    end if

    delx = img_crd(1, img_j) + x_tran(1, itran)
    dely = img_crd(2, img_j) + x_tran(2, itran)
    delz = img_crd(3, img_j) + x_tran(3, itran)
    !nxt_delx = img_crd(1, nxt_img_j) + x_tran(1, itran)
    !nxt_dely = img_crd(2, nxt_img_j) + x_tran(2, itran)
    !nxt_delz = img_crd(3, nxt_img_j) + x_tran(3, itran)

    delr2 = delx * delx + dely * dely + delz * delz

#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
    if (delr2 .lt. es_cut2) then
#else
    if (delr2 .lt. max_nb_cut2) then
#endif

      cgi_cgj = cgi * img_qterm(img_j)

#if 0
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
        !eed_stk = eed_stk + b0
        eed_simd = eed_simd + b0
#endif /* NEED_ENE */

#ifdef NEED_VIR
        !eedvir_stk = eedvir_stk - df * delr2
        eedvir_simd = eedvir_simd - df * delr2
#endif /* NEED_VIR */

      else
#endif /* if 0*/
        ! Do the Coulomb part of the direct sum using erfc spline table: 

        delr = sqrt(delr2)
        delrinv = 1.d0 / delr
        x = dxdr * delr
#if 0
        ind = int(eedtbdns_stk * x)
        dx = x - dble(ind) * del
        ind = ishft(ind, 2)             ! 4 * ind

        e3dx  = dx * eed_cub(3 + ind)
        e4dx2 = dx * dx *  eed_cub(4 + ind)

        switch = eed_cub(1 + ind) + &
                 dx * (eed_cub(2 + ind) + half * (e3dx + third * e4dx2))

        d_switch_dx = eed_cub(2 + ind) + e3dx + half * e4dx2
#endif /*if 0*/
#ifdef Faster_erfc
	exp_term = exp(-x * x)
        d_switch_dx = neg2InvSqrtPi * exp_term 
	t = 1.d0/(1.d0 + p * x )	
	switch =  t*(a1 + t * (a2 + t * (a3 + t * (a4 + t*a5)))) * exp_term 
#else
	switch = erfc(x) 
        d_switch_dx = (neg2InvSqrtPi) * exp(-x*x)
#endif

        b0 = cgi_cgj * delrinv * switch
        b1 = (b0 - cgi_cgj * d_switch_dx * dxdr)
#ifdef NEED_VIR
        !eedvir_stk = eedvir_stk + (- b1)
        eedvir_simd = eedvir_simd - b1
#endif /* NEED_VIR */
#ifdef NEED_ENE
        !eed_stk = eed_stk + b0
        eed_simd = eed_simd + b0
#endif /* NEED_ENE */

        df = b1 * delrinv * delrinv
#if 0
      end if
#endif

      dfx = delx * df
      dfy = dely * df
      dfz = delz * df

#ifdef NEED_VIR
#if 0
      vxx = vxx + (- delx * dfx)
!      vxy = vxy - delx * dfy
!      vxz = vxz - delx * dfz
      vyy = vyy + (- dely * dfy)
!      vyz = vyz - dely * dfz
      vzz = vzz + (- delz * dfz)
#endif
#if 1
      vxx_simd = vxx_simd - delx * dfx
!      vxy_simd = vxy_simd - delx * dfy
!      vxz_simd = vxz_simd - delx * dfz
      vyy_simd = vyy_simd - dely * dfy
!     vyz_simd = vyz_simd - dely * dfz
      vzz_simd = vzz_simd - delz * dfz
#endif
#endif /* NEED_VIR */
   

      !dumx = dumx + dfx
      !dumy = dumy + dfy
      !dumz = dumz + dfz
      dumx_simd = dumx_simd + dfx
      dumy_simd = dumy_simd + dfy
      dumz_simd = dumz_simd + dfz

!	!$OMP ATOMIC UPDATE
!      img_frc(1, img_j) = img_frc(1, img_j) + dfx
!	!$OMP ATOMIC UPDATE
!      img_frc(2, img_j) = img_frc(2, img_j) + dfy
!	!$OMP ATOMIC UPDATE
!      img_frc(3, img_j) = img_frc(3, img_j) + dfz

/* the following is for openmp */
      img_frc_thread(1, offset + img_j) = img_frc_thread(1, offset + img_j) + dfx
      img_frc_thread(2, offset + img_j) = img_frc_thread(2, offset + img_j) + dfy
      img_frc_thread(3, offset + img_j) = img_frc_thread(3, offset + img_j) + dfz

    end if
  end do

  iaci = ntypes_stk * (img_iac(img_i) - 1)

#ifdef NEED_ENE
	evdw_simd = evdw_stk
#endif
!$omp simd private(img_j, itran, delx, dely, delz, delr2, ic, cgi_cgj, ind, du, du2, du3, b0, df, delr2inv, dfx, dfy, dfz, delr,delrinv, x, dx, e3dx, e4dx2, switch, d_switch_dx, b1,r6, f6, f12 ) &
#ifdef NEED_VIR
!$omp& Reduction(-:eedvir_simd) &
!$omp& Reduction(+:vxx_simd, vyy_simd, vzz_simd) &
#endif
#ifdef NEED_ENE
!$omp& Reduction(+:eed_simd) &
!$omp& Reduction(+:evdw_simd) &
#endif
!$omp& Reduction(+:dumx_simd, dumy_simd, dumz_simd) 
  !do sublst_idx = ee_eval_cnt + 2, ee_eval_cnt + full_eval_cnt + 1
  do sublst_idx = ee_eval_cnt + 1, ee_eval_cnt + full_eval_cnt

!    img_j = nxt_img_j
  !  delx = nxt_delx
  !  dely = nxt_dely
  !  delz = nxt_delz

    if (common_tran .eq. 1) then
      img_j = ipairs_sublst(sublst_idx)
      !nxt_img_j = ipairs_sublst(sublst_idx)
      itran = 13
    else
      !nxt_img_j = iand(ipairs_sublst(sublst_idx), mask27)
      img_j = iand(ipairs_sublst(sublst_idx), mask27)
      itran = ishft(ipairs_sublst(sublst_idx), -27)
    end if

    !nxt_delx = img_crd(1, nxt_img_j) + x_tran(1, itran)
    !nxt_dely = img_crd(2, nxt_img_j) + x_tran(2, itran)
    !nxt_delz = img_crd(3, nxt_img_j) + x_tran(3, itran)

    delx = img_crd(1, img_j) + x_tran(1, itran)
    dely = img_crd(2, img_j) + x_tran(2, itran)
    delz = img_crd(3, img_j) + x_tran(3, itran)
    delr2 = delx * delx + dely * dely + delz * delz

    if (delr2 .lt. max_nb_cut2) then

      ic = ico(iaci + img_iac(img_j))

#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
      if (delr2 .lt. es_cut2) then
#endif

      cgi_cgj = cgi * img_qterm(img_j)
#if 0
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
        !eed_stk = eed_stk + b0
        eed_simd = eed_simd + b0
#endif /* NEED_ENE */

#ifdef NEED_VIR
        !eedvir_stk = eedvir_stk + (- df * delr2)
        eedvir_simd= eedvir_simd - df * delr2
#endif /* NEED_VIR */

        delr2inv = 1.d0 / delr2

      else
#endif   /*if 0*/   
        ! Do the Coulomb part of the direct sum using erfc spline table: 

        delr = sqrt(delr2)
        delrinv = 1.d0 / delr
        x = dxdr * delr
#if 0
        ind = int(eedtbdns_stk * x)
        dx = x - dble(ind) * del
        ind = ishft(ind, 2)             ! 4 * ind

        e3dx  = dx * eed_cub(3 + ind)
        e4dx2 = dx * dx *  eed_cub(4 + ind)

        switch = eed_cub(1 + ind) + &
                 dx * (eed_cub(2 + ind) + half * (e3dx + third * e4dx2))

        d_switch_dx = eed_cub(2 + ind) + e3dx + half * e4dx2
#endif /*if 0*/
#ifdef Faster_erfc
	exp_term = exp(-x * x)
        d_switch_dx = neg2InvSqrtPi * exp_term 
	t = 1.d0/(1.d0 + p * x )	
	switch =  t*(a1 + t * (a2 + t * (a3 + t * (a4 + t*a5)))) * exp_term 
#else
	switch = erfc(x) 
        d_switch_dx = (neg2InvSqrtPi) * exp(-x*x)
#endif /*Faster_erfc*/
        b0 = cgi_cgj * delrinv * switch
        b1 = b0 - cgi_cgj * d_switch_dx * dxdr
#ifdef NEED_ENE
        !eed_stk = eed_stk + b0
        eed_simd = eed_simd + b0
#endif /* NEED_ENE */
#ifdef NEED_VIR
        !eedvir_stk = eedvir_stk +( - b1)
        eedvir_simd = eedvir_simd - b1
#endif /* NEED_VIR */

        delr2inv = delrinv * delrinv
        df = b1 * delr2inv ! delrinv * delrinv
#if 0
      end if
#endif
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
      else
        delr2inv = 1.d0 / delr2
        df = 0.d0
      end if
#endif

#ifdef HAS_10_12
      if (ic .gt. 0) then
#endif
        r6 = delr2inv * delr2inv * delr2inv

        f6 = cn2(ic) * r6
        f12 = cn1(ic) * (r6 * r6)
#ifdef NEED_ENE
        !evdw_stk = evdw_stk + (f12 - f6)
        evdw_simd = evdw_simd + f12 - f6
#endif /* NEED_ENE */
        df = df + (12.d0 * f12 - 6.d0 * f6) * delr2inv

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
#endif /* NEED_ENE */
        df = df + (12.d0 * f12 - 10.d0 * f10) * delr2inv
      end if
#endif /* HAS_10_12 */

      dfx = delx * df
      dfy = dely * df
      dfz = delz * df

#ifdef NEED_VIR
#if 0
      vxx = vxx +(- delx * dfx)
!      vxy = vxy - delx * dfy
!      vxz = vxz - delx * dfz
      vyy = vyy +(- dely * dfy)
!      vyz = vyz - dely * dfz
      vzz = vzz +(- delz * dfz)
#endif
#if 1
      vxx_simd = vxx_simd - delx * dfx
   !   vxy_simd = vxy_simd - delx * dfy
   !   vxz_simd = vxz_simd - delx * dfz
      vyy_simd = vyy_simd - dely * dfy
   !   vyz_simd = vyz_simd - dely * dfz
      vzz_simd = vzz_simd - delz * dfz
#endif
#endif /* NEED_VIR */

#if 0
      dumx = dumx + dfx
      dumy = dumy + dfy
      dumz = dumz + dfz
#endif
#if 1
      dumx_simd = dumx_simd + dfx
      dumy_simd = dumy_simd + dfy
      dumz_simd = dumz_simd + dfz
#endif
!	!$OMP ATOMIC UPDATE
!      img_frc(1, img_j) = img_frc(1, img_j) + dfx
!	!$OMP ATOMIC UPDATE
!      img_frc(2, img_j) = img_frc(2, img_j) + dfy
!	!$OMP ATOMIC UPDATE
!      img_frc(3, img_j) = img_frc(3, img_j) + dfz

      img_frc_thread(1, offset + img_j) = img_frc_thread(1, offset +  img_j) + dfx
      img_frc_thread(2, offset + img_j) = img_frc_thread(2, offset + img_j) + dfy
      img_frc_thread(3, offset + img_j) = img_frc_thread(3, offset + img_j) + dfz

    end if

  end do

	dumx = dumx_simd
    	dumy = dumy_simd
      	dumz = dumz_simd
#ifdef NEED_VIR ! this is True
	eedvir_stk = eedvir_simd 
	vxx = vxx_simd 
   ! 	vxy = vxy_simd
   !   	vxz = vxz_simd
      	vyy = vyy_simd
   !   	vyz = vyz_simd
      	vzz = vzz_simd
#endif
#ifdef NEED_ENE 
	eed_stk= eed_simd
	evdw_stk= evdw_simd 
#endif
!	!$OMP ATOMIC UPDATE
  img_frc(1, img_i) = img_frc(1, img_i) - dumx
!	!$OMP ATOMIC UPDATE
  img_frc(2, img_i) = img_frc(2, img_i) - dumy
!	!$OMP ATOMIC UPDATE
   img_frc(3, img_i) = img_frc(3, img_i) - dumz
 ! ipairs_sublst(ee_eval_cnt + full_eval_cnt + 1) = saved_pairlist_val

  return


#ifdef BUILD_PAIRS_CALC_EFV
#undef NEED_ENE
#undef NEED_VIR
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
end subroutine pairs_calc_efv_2cut_offload
#else
end subroutine pairs_calc_efv_offload
#endif
#endif /* BUILD_PAIRS_CALC_EFV */

#ifdef BUILD_PAIRS_CALC_FV
#undef NEED_VIR
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
end subroutine pairs_calc_fv_2cut_offload
#else
end subroutine pairs_calc_fv_offload
#endif
#endif /* BUILD_PAIRS_CALC_FV */

#ifdef BUILD_PAIRS_CALC_F
#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
end subroutine pairs_calc_f_2cut_offload
#else
end subroutine pairs_calc_f_offload
#endif
#endif /* BUILD_PAIRS_CALC_F */

#endif
