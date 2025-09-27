#define BATCH_SIZE 192

#ifdef BUILD_PAIRS_CALC_EFV
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_midpoint_efv
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file.
!*******************************************************************************
#ifdef _OPENMP_

#ifdef pmemd_SPDP
subroutine pairs_calc_midpoint_efv(img_i,ee_eval_cnt, full_eval_cnt, &
                                evdw_stk,eed_stk, eedvir_stk, ehb_stk, &
                                vxx, vyy, vzz,vxy,vxz, vyz, &
                                img_frc, img_crd_q_sp, ef_tbl, &
                              eed_cub, ico, ipairs_sublst, &
                              img_iac, cn1, cn2, x_i, y_i, z_i)
#else
subroutine pairs_calc_midpoint_efv(img_i,ee_eval_cnt, full_eval_cnt, &
                                evdw_stk,eed_stk, eedvir_stk, ehb_stk, &
                                vxx, vyy, vzz,vxy,vxz, vyz, &
                                img_frc, img_crd, img_qterm, ef_tbl, &
                              eed_cub, ico, ipairs_sublst, &
                              img_iac, cn1, cn2, x_i, y_i, z_i)
#endif /* pmemd_SPDP */

#else
#ifdef pmemd_SPDP /* _OPENMP_ : SPDP w/o OMP */
 subroutine pairs_calc_midpoint_efv(img_frc, img_crd_q_sp, ef_tbl, eed_cub, &
                          ico, ipairs_sublst, img_iac, cn1, cn2, x_i, y_i, z_i)
#else /* DPDP w/o OMP */
subroutine pairs_calc_midpoint_efv(img_frc, img_crd, img_qterm, ef_tbl, eed_cub, &
                          ico, ipairs_sublst, img_iac, cn1, cn2, x_i, y_i, z_i)
#endif                          
#endif /* _OPENMP_ */

#endif /* BUILD_PAIRS_CALC_EFV */

#ifdef BUILD_PAIRS_CALC_FV
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_midpoint_fv
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file.
!*******************************************************************************
#ifdef _OPENMP_

#ifdef pmemd_SPDP
subroutine pairs_calc_midpoint_fv(img_i, ee_eval_cnt, full_eval_cnt, & 
                                vxx, vyy, vzz,vxy,vxz, vyz, &
                                img_frc, img_crd_q_sp, f_tbl, &
                               eed_cub, ico, ipairs_sublst, &
                              img_iac, cn1, cn2, x_i, y_i, z_i)
#else
subroutine pairs_calc_midpoint_fv(img_i, ee_eval_cnt, full_eval_cnt, & 
                                vxx, vyy, vzz,vxy,vxz, vyz, &
                                img_frc, img_crd, img_qterm, f_tbl, &
                               eed_cub, ico, ipairs_sublst, &
                              img_iac, cn1, cn2, x_i, y_i, z_i)
#endif /* pmemd_SPDP */

#else
#ifdef pmemd_SPDP
 subroutine pairs_calc_midpoint_fv(img_frc, img_crd_q_sp , f_tbl, eed_cub, &
                         ico, ipairs_sublst, img_iac, cn1, cn2, x_i, y_i, z_i)
#else
subroutine pairs_calc_midpoint_fv(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
                         ico, ipairs_sublst, img_iac, cn1, cn2, x_i, y_i, z_i)
#endif

#endif /* _OPENMP_ */
#endif /* BUILD_PAIRS_CALC_FV */

#ifdef BUILD_PAIRS_CALC_F
!*******************************************************************************
!
! Internal Subroutine:  pairs_calc_midpoint_f
!
! Description:  Direct force computation on one atom.  Different implementations
!               are covered by conditional compilation in this file.
!*******************************************************************************

#ifdef _OPENMP_

#ifdef pmemd_SPDP
subroutine pairs_calc_midpoint_f(img_i, ee_eval_cnt, full_eval_cnt, &
                                img_frc, img_crd_q_sp, f_tbl, eed_cub, &
                                ico, ipairs_sublst, img_iac,&
                                cn1, cn2, x_i, y_i, z_i)
#else
subroutine pairs_calc_midpoint_f(img_i, ee_eval_cnt, full_eval_cnt, &
                                img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
                                ico, ipairs_sublst, img_iac,&
                                cn1, cn2, x_i, y_i, z_i)
#endif /* pmemd_SPDP */

#else
#ifdef pmemd_SPDP
subroutine pairs_calc_midpoint_f(img_frc, img_crd_q_sp, f_tbl, eed_cub, &
                        ico, ipairs_sublst, img_iac, cn1, cn2, x_i, y_i, z_i)
#else                        
subroutine pairs_calc_midpoint_f(img_frc, img_crd, img_qterm, f_tbl, eed_cub, &
                        ico, ipairs_sublst, img_iac, cn1, cn2, x_i, y_i, z_i)
#endif                        
#endif /* _OPENMP_ */
#endif /* BUILD_PAIRS_CALC_F */
!  use mdin_ctrl_dat_mod, only : fswitch
   use processor_mod, only : proc_atm_to_full_list
#ifdef pmemd_SPDP   
   use processor_mod, only : neg2InvSqrtPi
#endif /* pmemd_SPDP */
   use parallel_dat_mod
  implicit none

! Formal arguments:

#ifdef pmemd_SPDP
  pme_float, intent(in)         :: img_crd_q_sp(4, *)
  pme_double                    :: img_frc(3, *)
  pme_float                    :: x_i, y_i, z_i
#else
  double precision              :: img_frc(3, *)
  double precision              :: x_i, y_i, z_i
  double precision, intent(in)  :: img_crd(3, *)
  double precision, intent(in)  :: img_qterm(*)
#endif  /* pmemd_SPDP */
#ifdef NEED_ENE
  double precision, intent(in)  :: ef_tbl(*)
#ifdef _OPENMP_
  double precision              :: evdw_stk,eed_stk, ehb_stk
#endif
#else
  double precision, intent(in)  :: f_tbl(*)
#endif

#ifdef NEED_VIR
#ifdef _OPENMP_
  double precision              :: eedvir_stk, vxx, vyy, vzz, vxy,vxz, vyz
#endif
#endif

#ifdef _OPENMP_
  integer, intent(in)           :: img_i, ee_eval_cnt, full_eval_cnt
#endif

  double precision, intent(in)  :: eed_cub(*)
  integer, intent(in)           :: ico(*)
  integer                       :: ipairs_sublst(*)
  integer, intent(in)           :: img_iac(*)

#ifdef pmemd_SPDP  
 pme_float , intent(in)  :: cn1(*), cn2(*)
! pme_float  ,intent(in)  :: x_tran(1:3, 0:17)
#else
 double precision , intent(in)  :: cn1(*), cn2(*)
 !double precision ,intent(in)  :: x_tran(1:3, 0:17)
#endif


! Local variables:

  integer, parameter            :: mask27 = int(Z"07FFFFFF")
  double precision      :: dumx, dumy, dumz
#ifndef pmemd_SPDP  
  double precision, parameter   :: half = 1.d0/2.d0
  double precision, parameter   :: third = 1.d0/3.d0
  double precision, parameter   :: one  = 1.d0
  double precision, parameter   :: zero = 0.d0
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
#else

 pme_float, parameter   :: half = 0.5
 pme_float, parameter   :: third = 0.3333
 pme_float, parameter   :: one = 1.0
 pme_float, parameter   :: zero = 0.0
 pme_float :: cgi
 pme_float :: cgi_cgj
 pme_float :: b0, b1
 pme_float :: df
 pme_float :: dfx, dfy, dfz
 pme_float :: f6, r6, f12
#ifdef HAS_10_12
 pme_float :: f10, r10
#endif
 pme_float   :: du, du2, du3 ! 'u' is another name for delr2
 pme_float   :: del_efs
 pme_float   :: dens_efs
 pme_float   :: lowest_efs_u
  ! Variables used with erfc switch table; name are historical:
 pme_float   :: switch
 pme_float   :: d_switch_dx
 pme_float   :: x, dx, e3dx, e4dx2
 pme_float   :: delr, delrinv

#endif /* pmemd_SPDP  */

  integer               :: iaci
  integer               :: ic
  integer               :: ind
  integer               :: nxt_img_j, img_j
  integer               :: itran
  integer               :: sublst_idx
  integer               :: saved_pairlist_val

#ifdef MP_VEC
#ifdef NEED_ENE
  !double precision      :: eed_simd=0.d0, evdw_simd=0.d0,ehb_simd=0 
  double precision      :: eed_simd, evdw_simd, ehb_simd 
#endif
#ifdef NEED_VIR
  !double precision      ::eedvir_simd,  vxx_simd, vyy_simd, vzz_simd, vxy_simd, vxz_simd, vyz_simd
  double precision      :: eedvir_simd, vxx_simd, vyy_simd, vzz_simd, vxy_simd, vxz_simd, vyz_simd
#endif
#endif /* MP_VEC */
  pme_float      :: t, exp_term
  pme_float      :: nxt_delx, nxt_dely, nxt_delz
  pme_float      :: delx, dely, delz, delr2, delr2inv
  pme_float      :: delx_arr(BATCH_SIZE), dely_arr(BATCH_SIZE), delz_arr(BATCH_SIZE)
  pme_float      :: j_arr(BATCH_SIZE)!, ic_arr(BATCH_SIZE)

#ifdef pmemd_SPDP
  pme_float      :: delq_arr(BATCH_SIZE)     
#else
  double precision :: cgi_cgj_arr(BATCH_SIZE) 
#endif
integer             :: sublst_idx1, jn, jmax, j_out, counter, loop_count, my_index
#ifdef pmemd_SPDP
 !double precision       , parameter  :: PI = 3.1415926535897930
 pme_float  :: a1=0.2548295920e0, a2=-0.2844967360e0, a3=1.4214137410e0, a4=-1.4531520270e0, &
                           a5=1.0614054290e0, p=0.32759110e0
#endif /* pmemd_SPDP */


  dens_efs = efs_tbl_dens
  del_efs = one / dens_efs
  lowest_efs_u = lowest_efs_delr2

! Precompute some constants


! First loop over the ee evaluation-only pairs:

  dumx = zero
  dumy = zero
  dumz = zero

#ifdef MP_VEC
#ifdef NEED_VIR
  eedvir_simd = eedvir_stk
  vxx_simd = vxx
  vxy_simd = vxy
  vxz_simd = vxz
  vyy_simd = vyy
  vyz_simd = vyz
  vzz_simd = vzz
#endif /*NEED_VIR*/
#ifdef NEED_ENE
  eed_simd = eed_stk
  evdw_simd = evdw_stk
  ehb_simd = ehb_stk
#endif/*NEED_ENE*/
#endif /* MP_VEC */

#ifdef pmemd_SPDP
  cgi = img_crd_q_sp(4,img_i)
#else
  cgi = img_qterm(img_i)
#endif /* pmemd_SPDP */


  ! The pairlist must have one dummy end entry to cover reading past the
  ! end of the list...
  
  !saved_pairlist_val = ipairs_sublst(ee_eval_cnt + full_eval_cnt + 1)

  !ipairs_sublst(ee_eval_cnt + full_eval_cnt + 1) = &
  !  ipairs_sublst(ee_eval_cnt + full_eval_cnt)
!  if (common_tran .eq. 1) then
!    nxt_img_j = ipairs_sublst(1)
!    itran = 13
!  else
    !nxt_img_j = iand(ipairs_sublst(1), mask27)
    !itran = ishft(ipairs_sublst(1), -27)
!  end if
!  nxt_delx = img_crd(1, nxt_img_j) + x_tran(1, itran)
!  nxt_dely = img_crd(2, nxt_img_j) + x_tran(2, itran)
!  nxt_delz = img_crd(3, nxt_img_j) + x_tran(3, itran)


!####################################################
!
!
!###################################################

 do j_out = 1, ee_eval_cnt, BATCH_SIZE
 
     jmax = j_out + (BATCH_SIZE - 1)
     if (jmax .gt. ee_eval_cnt) then 
       jmax = ee_eval_cnt
       loop_count = jmax - j_out + 1
     else
       loop_count = BATCH_SIZE
     end if
     counter = 0
  !do sublst_idx = j_out, jmax
  do sublst_idx = 1, loop_count 
     
     my_index = sublst_idx + j_out - 1
     img_j =ipairs_sublst(my_index)
#ifdef pmemd_SPDP
    delx = img_crd_q_sp(1, img_j) -x_i !+ x_tran(1, itran)
    dely = img_crd_q_sp(2, img_j) -y_i !+ x_tran(2, itran)
    delz = img_crd_q_sp(3, img_j) -z_i !+ x_tran(3, itran)
#else    
    delx = img_crd(1, img_j) - x_i ! + x_tran(1, itran)
    dely = img_crd(2, img_j) - y_i !+ x_tran(2, itran)
    delz = img_crd(3, img_j)  -z_i  !+ x_tran(3, itran)
#endif /* pmemd_SPDP */
    delr2 = delx * delx + dely * dely + delz * delz

#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
    if (delr2 .lt. es_cut2) then
#else
    if (delr2 .lt. max_nb_cut2) then
#endif
      counter = counter + 1
#ifdef pmemd_SPDP
      delx_arr(counter ) = delx 
      dely_arr(counter ) = dely
      delz_arr(counter ) = delz 
      delq_arr(counter ) = cgi * img_crd_q_sp(4,img_j)
#else
      delx_arr(counter ) = delx 
      dely_arr(counter ) = dely
      delz_arr(counter ) = delz 
      cgi_cgj_arr(counter ) = cgi * img_qterm(img_j)
#endif /* pmemd_SPDP */
      j_arr(counter ) = img_j
   end if 
  end do 

#ifdef MP_VEC
!$omp simd private(img_j, itran, delx, dely, delz, delr2, ic, cgi_cgj, ind, du, du2, du3, b0, df, &
!$omp& delr2inv, dfx, dfy, dfz, delr, delrinv, x, dx, e3dx, e4dx2, switch, d_switch_dx, b1, t, exp_term) & 
#ifdef NEED_VIR
!$omp& Reduction(-:eedvir_simd) &
!$omp& Reduction(-:vxx_simd, vyy_simd, vzz_simd,vxy_simd,vxz_simd, vyz_simd) &
#endif
#ifdef NEED_ENE
!$omp& Reduction(+:eed_simd) &
#endif
!$omp& Reduction(+:dumx, dumy, dumz) 
#endif /*MP_VEC*/
  do sublst_idx1 = 1, counter 
     
#ifdef pmemd_SPDP
      delx = delx_arr(sublst_idx1)
      dely = dely_arr(sublst_idx1) 
      delz = delz_arr(sublst_idx1) 
#else
      delx = delx_arr(sublst_idx1)  
      dely = dely_arr(sublst_idx1)  
      delz = delz_arr(sublst_idx1)  
#endif /* pmemd_SPDP */
      img_j = j_arr(sublst_idx1) 

    !if (common_tran .eq. 1) then
    !  img_j = ipairs_sublst(sublst_idx)
    !  itran = 13
    !else
    !==>  img_j =ipairs_sublst(sublst_idx)
     ! img_j = iand(ipairs_sublst(sublst_idx), mask27)
     ! itran = ishft(ipairs_sublst(sublst_idx), -27)
    !end if
  
!#ifdef pmemd_SPDP
!    delx = img_crd_q_sp(1, img_j) -x_i !+ x_tran(1, itran)
!    dely = img_crd_q_sp(2, img_j) -y_i !+ x_tran(2, itran)
!    delz = img_crd_q_sp(3, img_j) -z_i !+ x_tran(3, itran)
!#else    
!    delx = img_crd(1, img_j) - x_i ! + x_tran(1, itran)
!    dely = img_crd(2, img_j) - y_i !+ x_tran(2, itran)
!    delz = img_crd(3, img_j)  -z_i  !+ x_tran(3, itran)
!#endif /* pmemd_SPDP */

    delr2 = delx * delx + dely * dely + delz * delz

!#ifdef pmemd_SPDP
!      cgi_cgj = cgi * img_crd_q_sp(4,img_j)
!#else
!      cgi_cgj = cgi * img_qterm(img_j)
!#endif /* pmemd_SPDP */
  
#ifdef pmemd_SPDP
      cgi_cgj = delq_arr(sublst_idx1) 
#else
      cgi_cgj = cgi_cgj_arr(sublst_idx1)  
#endif /* pmemd_SPDP */
  
#ifndef pmemd_SPDP
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
#ifdef MP_VEC
        eed_simd = eed_simd + b0
#else
        eed_stk = eed_stk + b0
#endif
#endif

#ifdef NEED_VIR
#ifdef MP_VEC
        eedvir_simd = eedvir_simd - df * delr2 
#else
        eedvir_stk = eedvir_stk - df * delr2 
#endif
#endif 
      else
#endif /* pmemd_SPDP */
        ! Do the Coulomb part of the direct sum using erfc spline table: 

        delr = sqrt(delr2)
        delrinv = one / delr

        x = dxdr * delr
#ifndef pmemd_SPDP
        ind = int(eedtbdns_stk * x)
        dx = x - dble(ind) * del
        ind = ishft(ind, 2)             ! 4 * ind

        e3dx  = dx * eed_cub(3 + ind)
        e4dx2 = dx * dx *  eed_cub(4 + ind)

        switch = eed_cub(1 + ind) + &
                 dx * (eed_cub(2 + ind) + half * (e3dx + third * e4dx2))

        d_switch_dx = eed_cub(2 + ind) + e3dx + half * e4dx2
#else /* pmemd_SPDP */
#ifdef pmemd_SPDP
        exp_term = exp(-x * x)
        d_switch_dx = neg2InvSqrtPi * exp_term 
        t = one /(one + p * x )
        switch =  t*(a1 + t * (a2 + t * (a3 + t * (a4 + t*a5)))) * exp_term 
#else
        exp_term = exp(-x * x)
        d_switch_dx = neg2InvSqrtPi * exp_term 
        t = 1.d0/(1.d0 + p * x )
        switch =  t*(a1 + t * (a2 + t * (a3 + t * (a4 + t*a5)))) * exp_term 
#endif /* pmemd_SPDP */        
#endif /* pmemd_SPDP */
        b0 = cgi_cgj * delrinv * switch
        b1 = (b0 - cgi_cgj * d_switch_dx * dxdr)


#ifdef NEED_VIR
#ifdef MP_VEC
        eedvir_simd = eedvir_simd - b1
#else
        eedvir_stk = eedvir_stk - b1
#endif
#endif 

#ifdef NEED_ENE
#ifdef MP_VEC
        eed_simd = eed_simd + b0
#else
        eed_stk = eed_stk + b0
#endif
#endif

        df = b1 * delrinv * delrinv


#ifndef pmemd_SPDP
      end if
#endif /* pmemd_SPDP */

      dfx = delx * df
      dfy = dely * df
      dfz = delz * df

#ifdef NEED_VIR
#ifdef MP_VEC
      vxx_simd = vxx_simd - delx * dfx
      vxy_simd = vxy_simd - delx * dfy
      vxz_simd = vxz_simd - delx * dfz
      vyy_simd = vyy_simd - dely * dfy
      vyz_simd = vyz_simd - dely * dfz
      vzz_simd = vzz_simd - delz * dfz
#else
      vxx = vxx - delx * dfx
      vxy = vxy - delx * dfy
      vxz = vxz - delx * dfz
      vyy = vyy - dely * dfy
      vyz = vyz - dely * dfz
      vzz = vzz - delz * dfz
#endif
#endif

      dumx = dumx + dfx
      dumy = dumy + dfy
      dumz = dumz + dfz

      img_frc(1, img_j) = img_frc(1, img_j) + dfx
      img_frc(2, img_j) = img_frc(2, img_j) + dfy
      img_frc(3, img_j) = img_frc(3, img_j) + dfz

  end do /*sublst_idx = 1, (counter) (old: ee_eval_cnt) */
end do /** do j_out = 1, ee_eval_cnt, BATCH_SIZE */
  iaci = ntypes_stk * (img_iac(img_i) - 1)

!####################################################
!
!
!###################################################


 do j_out = ee_eval_cnt + 1, ee_eval_cnt + full_eval_cnt, BATCH_SIZE
 
     jmax = j_out + (BATCH_SIZE - 1)
     if (jmax .gt. ee_eval_cnt+full_eval_cnt) then 
       jmax = ee_eval_cnt + full_eval_cnt
       loop_count = jmax - j_out + 1
     else
       loop_count = BATCH_SIZE
     end if

     counter = 0
  !do sublst_idx = j_out, jmax
  do sublst_idx = 1, loop_count 
     
     my_index = sublst_idx + j_out - 1
     img_j =ipairs_sublst(my_index)
#ifdef pmemd_SPDP
    delx = img_crd_q_sp(1, img_j) -x_i !+ x_tran(1, itran)
    dely = img_crd_q_sp(2, img_j) -y_i !+ x_tran(2, itran)
    delz = img_crd_q_sp(3, img_j) -z_i !+ x_tran(3, itran)
#else    
    delx = img_crd(1, img_j) - x_i ! + x_tran(1, itran)
    dely = img_crd(2, img_j) - y_i !+ x_tran(2, itran)
    delz = img_crd(3, img_j)  -z_i  !+ x_tran(3, itran)
#endif /* pmemd_SPDP */
    delr2 = delx * delx + dely * dely + delz * delz

    if (delr2 .lt. max_nb_cut2) then
      counter = counter + 1
#ifdef pmemd_SPDP
      delx_arr(counter ) = delx 
      dely_arr(counter ) = dely
      delz_arr(counter ) = delz 
      delq_arr(counter ) = cgi * img_crd_q_sp(4,img_j)
#else
      delx_arr(counter ) = delx 
      dely_arr(counter ) = dely
      delz_arr(counter ) = delz 
      cgi_cgj_arr(counter ) = cgi * img_qterm(img_j)
#endif /* pmemd_SPDP */
      j_arr(counter ) = img_j
      !ic_arr(counter) = ico(iaci + img_iac(img_j))
   end if 
  end do 

#ifdef MP_VEC
!$omp simd private( img_j, itran, delx, dely, delz, delr2, ic, cgi_cgj, ind, du, du2, du3, b0, &
!$omp& df, delr2inv, dfx, dfy, dfz, delr, delrinv, x, dx, e3dx, e4dx2, switch, d_switch_dx, b1, t, exp_term) & 
#ifdef NEED_VIR
!$omp& Reduction(-:eedvir_simd) &
!$omp& Reduction(-:vxx_simd, vyy_simd, vzz_simd,vxy_simd, vyz_simd, vxz_simd) &
#endif
#ifdef NEED_ENE
!$omp& Reduction(+:eed_simd, evdw_simd ) &
#ifdef HAS_10_12
!$omp& Reduction(+:ehb_simd ) &
#endif
#endif
!$omp& Reduction(+:dumx, dumy, dumz) 
#endif /* MP_VEC */
  !do sublst_idx = ee_eval_cnt + 1, ee_eval_cnt + full_eval_cnt 
  do sublst_idx1 = 1 , counter 
   ! if (common_tran .eq. 1) then
   !   nxt_img_j = ipairs_sublst(sublst_idx)
   !   itran = 13
   ! else
      ! => img_j = ipairs_sublst(sublst_idx)
      !img_j = iand(ipairs_sublst(sublst_idx), mask27)
      !itran = ishft(ipairs_sublst(sublst_idx), -27)
   ! end if
!#ifdef pmemd_SPDP
!    delx = img_crd_q_sp(1, img_j) -x_i !+ x_tran(1, itran)
!    dely = img_crd_q_sp(2, img_j) -y_i !+ x_tran(2, itran)
!    delz = img_crd_q_sp(3, img_j) -z_i !+ x_tran(3, itran)
!#else
!    delx = img_crd(1, img_j) -x_i !+ x_tran(1, itran)
!    dely = img_crd(2, img_j) -y_i !+ x_tran(2, itran)
!    delz = img_crd(3, img_j) -z_i !+ x_tran(3, itran)
!#endif  /* pmemd_SPDP */

#ifdef pmemd_SPDP
      delx = delx_arr(sublst_idx1)
      dely = dely_arr(sublst_idx1) 
      delz = delz_arr(sublst_idx1) 
#else
      delx = delx_arr(sublst_idx1)  
      dely = dely_arr(sublst_idx1)  
      delz = delz_arr(sublst_idx1)  
#endif /* pmemd_SPDP */
      img_j = j_arr(sublst_idx1) 

    delr2 = delx * delx + dely * dely + delz * delz

    !if (delr2 .lt. max_nb_cut2) then
      ic = ico(iaci + img_iac(img_j))
      !ic = ic_arr(sublst_idx1)

#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
      if (delr2 .lt. es_cut2) then
#endif


!#ifdef pmemd_SPDP
!      cgi_cgj = cgi * img_crd_q_sp(4,img_j)
!#else
!      cgi_cgj = cgi * img_qterm(img_j)
!#endif /* pmemd_SPDP */

#ifdef pmemd_SPDP
      cgi_cgj = delq_arr(sublst_idx1) 
#else
      cgi_cgj = cgi_cgj_arr(sublst_idx1)  
#endif /* pmemd_SPDP */
  
#ifndef pmemd_SPDP
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
#ifdef MP_VEC
        eed_simd = eed_simd + b0
#else
        eed_stk = eed_stk + b0
#endif
#endif

#ifdef NEED_VIR
#ifdef MP_VEC
        eedvir_simd = eedvir_simd - df * delr2
#else
        eedvir_stk = eedvir_stk - df * delr2
#endif 
#endif /* NEED_VIR */
        delr2inv = one / delr2
      else    
#endif /* pmemd_SPDP */

        ! Do the Coulomb part of the direct sum using erfc spline table: 

        delr = sqrt(delr2)
        delrinv = one / delr

        x = dxdr * delr
#ifndef pmemd_SPDP
        ind = int(eedtbdns_stk * x)
        dx = x - dble(ind) * del
        ind = ishft(ind, 2)             ! 4 * ind

        e3dx  = dx * eed_cub(3 + ind)
        e4dx2 = dx * dx *  eed_cub(4 + ind)

        switch = eed_cub(1 + ind) + &
                 dx * (eed_cub(2 + ind) + half * (e3dx + third * e4dx2))

        d_switch_dx = eed_cub(2 + ind) + e3dx + half * e4dx2

#else /* pmemd_SPDP */
#ifdef pmemd_SPDP
        exp_term = exp(-x * x)
        d_switch_dx = neg2InvSqrtPi * exp_term 
        t =one / (one + p * x )
        switch =  t*(a1 + t * (a2 + t * (a3 + t * (a4 + t*a5)))) * exp_term 
#else        
        exp_term = exp(-x * x)
        d_switch_dx = neg2InvSqrtPi * exp_term 
        t = 1.d0 / (1.d0 + p * x )
        switch =  t*(a1 + t * (a2 + t * (a3 + t * (a4 + t*a5)))) * exp_term 
#endif /* pmemd_SPDP */
#endif /* pmemd_SPDP */

        b0 = cgi_cgj * delrinv * switch
        b1 = (b0 - cgi_cgj * d_switch_dx * dxdr)

#ifdef NEED_ENE
#ifdef MP_VEC
        eed_simd = eed_simd + b0
#else
        eed_stk = eed_stk + b0
#endif 
#endif /* NEED_ENE */

#ifdef NEED_VIR
#ifdef MP_VEC
        eedvir_simd = eedvir_simd - b1
#else
        eedvir_stk = eedvir_stk - b1
#endif 
#endif /* NEED_VIR */

        delr2inv = delrinv * delrinv
        df = b1 * delrinv * delrinv

#ifndef pmemd_SPDP
      end if
#endif /* pmemd_SPDP */

#ifdef BUILD_PAIRS_CALC_NOVEC_2CUT
      else
        delr2inv = one / delr2
        df = zero
      end if
#endif

#ifdef HAS_10_12
      if (ic .gt. 0) then
#endif
        r6 = delr2inv * delr2inv * delr2inv
        f6 = cn2(ic) * r6
        f12 = cn1(ic) * (r6 * r6)
#ifdef NEED_ENE
#ifdef MP_VEC
          evdw_simd = evdw_simd + f12 - f6
#else
          evdw_stk = evdw_stk + f12 - f6
#endif 
#endif /* NEED_ENE */
          df = df + (12.0 * f12 - 6.0 * f6) * delr2inv

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
        ehb_simd = ehb_simd + f12 - f10
        !ehb_stk = ehb_stk + f12 - f10
      end if
#endif /* NEED_ENE */
        df = df + (12.0 * f12 - 10.0 * f10) * delr2inv
      end if
#endif /* HAS_10_12 */
      dfx = delx * df
      dfy = dely * df
      dfz = delz * df

#ifdef NEED_VIR
#ifdef MP_VEC
      vxx_simd = vxx_simd - delx * dfx
      vxy_simd = vxy_simd - delx * dfy
      vxz_simd = vxz_simd - delx * dfz
      vyy_simd = vyy_simd - dely * dfy
      vyz_simd = vyz_simd - dely * dfz
      vzz_simd = vzz_simd - delz * dfz
#else
      vxx = vxx - delx * dfx
      vxy = vxy - delx * dfy
      vxz = vxz - delx * dfz
      vyy = vyy - dely * dfy
      vyz = vyz - dely * dfz
      vzz = vzz - delz * dfz
#endif
#endif

      dumx = dumx + dfx
      dumy = dumy + dfy
      dumz = dumz + dfz

      img_frc(1, img_j) = img_frc(1, img_j) + dfx
      img_frc(2, img_j) = img_frc(2, img_j) + dfy
      img_frc(3, img_j) = img_frc(3, img_j) + dfz
    !end if

  end do /* 1, counter*/
 end do /* ee_eval_cnt + 1, ee_eval_cnt + full_eval_cnt, BATCH_SIZE */

  img_frc(1, img_i) = img_frc(1, img_i) - dumx
  img_frc(2, img_i) = img_frc(2, img_i) - dumy
  img_frc(3, img_i) = img_frc(3, img_i) - dumz

#ifdef MP_VEC
#ifdef NEED_VIR 
  eedvir_stk = eedvir_simd 
  vxx =  vxx_simd 
  vxy =  vxy_simd 
  vxz =  vxz_simd 
  vyy =  vyy_simd
  vyz =  vyz_simd
  vzz =  vzz_simd
#endif
#ifdef NEED_ENE 
  eed_stk=  eed_simd
  evdw_stk= evdw_simd 
  ehb_stk= ehb_simd 
#endif
#endif /*MP_VEC*/

!  ipairs_sublst(ee_eval_cnt + full_eval_cnt + 1) = saved_pairlist_val

  return

#ifdef BUILD_PAIRS_CALC_EFV
#undef NEED_ENE
#undef NEED_VIR
end subroutine pairs_calc_midpoint_efv
#endif /* BUILD_PAIRS_CALC_EFV */

#ifdef BUILD_PAIRS_CALC_FV
#undef NEED_VIR
end subroutine pairs_calc_midpoint_fv
#endif /* BUILD_PAIRS_CALC_FV */

#ifdef BUILD_PAIRS_CALC_F
end subroutine pairs_calc_midpoint_f
#endif /* BUILD_PAIRS_CALC_F */


