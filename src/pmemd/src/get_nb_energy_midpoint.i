!*******************************************************************************
!
! Subroutine:  get_nb_energy
!
! Description:
!              
! The main routine for non bond energy (vdw and hbond) as well as direct part
! of ewald sum.  It is structured for parallelism.
!
!*******************************************************************************


#ifdef pmemd_SPDP
subroutine get_nb_energy_midpoint(img_frc, img_crd_q_sp , eed_cub, &
                         ipairs, need_pot_enes, need_virials, &
                         eed, evdw, ehb, eedvir, virial)
#else
subroutine get_nb_energy_midpoint(img_frc, img_crd, img_qterm, eed_cub, &
                         ipairs, need_pot_enes, need_virials, &
                         eed, evdw, ehb, eedvir, virial)
#endif /* pmemd_SPDP */

  use img_mod
  use timers_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use ene_frc_splines_mod
  use processor_mod, only: proc_num_atms, proc_ghost_num_atms, proc_iac, start_index &
      ,bkt_atm_cnt, bkt_atm_lst, int_bkt_dimx, int_bkt_dimy, int_bkt_dimz 
#ifdef pmemd_SPDP
  use processor_mod, only: proc_gbl_cn1_sp, proc_gbl_cn2_sp
#endif

#ifdef _OPENMP_
  use omp_lib
#endif

  implicit none

! Formal arguments:

  double precision, intent(in out) :: img_frc(3, *)
#ifdef pmemd_SPDP
  pme_float, intent(in)     :: img_crd_q_sp(4, *)
#else
  double precision, intent(in)     :: img_crd(3, *)
  double precision , intent(in)    :: img_qterm(*)
#endif /* pmemd_SPDP */
  double precision, intent(in)     :: eed_cub(*)
  integer                          :: ipairs(*)
  logical, intent(in)              :: need_pot_enes
  logical, intent(in)              :: need_virials
  double precision, intent(out)    :: eed
  double precision, intent(out)    :: evdw
  double precision, intent(out)    :: ehb
  double precision, intent(out)    :: eedvir
  double precision, intent(out)    :: virial(3, 3)

! Local variables and parameters:
#ifdef pmemd_SPDP
  pme_float                       del
  pme_float                       dxdr
  pme_float                       max_nb_cut2, es_cut2, es_cut
  pme_float                       x_i, y_i, z_i
!  pme_float                       x_tran(1:3, 0:17)
#else
  double precision                  del
  double precision                  dxdr
  double precision                  max_nb_cut2, es_cut2, es_cut
  double precision                  x_i, y_i, z_i
!  double precision                  x_tran(1:3, 0:17)
#endif /* pmemd_SPDP */   
  double precision      eedtbdns_stk
  double precision      eedvir_stk, eed_stk, evdw_stk, ehb_stk
  double precision      vxx, vxy, vxz, vyy, vyz, vzz
  integer               i
  integer               ipairs_idx
  integer               ntypes_stk
  integer               img_i
  integer               ee_eval_cnt
  integer               full_eval_cnt
  integer               common_tran    ! flag - 1 if translation not needed
  logical               cutoffs_equal
  integer               bkt, atm_i_idx, bkt_cnt

#ifdef _OPENMP_
   integer              :: tid , sd_num_atms, arr_ind
#endif
! Proc  dont have this check, do we need?
!#ifdef MPI
!  if (my_img_lo .gt. my_img_hi) return
!#endif /* MPI */

  ntypes_stk = ntypes 
  bkt_cnt = (int_bkt_dimx+2)*(int_bkt_dimy+2)*(int_bkt_dimz+2) 

  eedvir_stk = 0.d0
  eed_stk = 0.d0
  evdw_stk = 0.d0
  ehb_stk = 0.d0

  dxdr = ew_coeff
  eedtbdns_stk = eedtbdns
  del = 1.d0 / eedtbdns_stk
  max_nb_cut2 = vdw_cutoff * vdw_cutoff
  es_cut = es_cutoff
  es_cut2 = es_cut * es_cut
  cutoffs_equal = (vdw_cutoff .eq. es_cutoff)

  vxx = 0.d0
!  vxy = 0.d0
!  vxz = 0.d0
  vyy = 0.d0
!  vyz = 0.d0
  vzz = 0.d0
! ipairs_idx =1
#ifdef _OPENMP_
 sd_num_atms = proc_num_atms + proc_ghost_num_atms
#endif

  if (need_pot_enes) then
#ifdef _OPENMP_
!$omp parallel do default(shared) & 
!$omp& schedule (dynamic,16) &
!!$omp& shared(bkt_cnt, bkt_atm_cnt, bkt_atm_lst, proc_num_atms , proc_ghost_num_atms, start_index, ipairs, &
!!$omp& img_crd, img_qterm, img_frc, cutoffs_equal, efs_tbl, &
!!$omp& eed_cub, typ_ico, proc_iac, gbl_cn1, gbl_cn2, &
!!$omp& sd_num_atms) &
!$omp& private(bkt, atm_i_idx, img_i, ipairs_idx, ee_eval_cnt, &
!$omp&  full_eval_cnt, x_i, y_i, z_i, i, tid, arr_ind ) &
!$omp&  reduction(+: eedvir_stk) &
!$omp&  reduction(+: eed_stk, evdw_stk, ehb_stk) &
!$omp&  reduction(+: vxx, vyy, vzz,vxy,vxz, vyz)
#endif /* _OPENMP_ */
    do bkt = 1,bkt_cnt
    do atm_i_idx = 1, bkt_atm_cnt(bkt)
      img_i = bkt_atm_lst(atm_i_idx, bkt)
!   do img_i = 1 , proc_num_atms + proc_ghost_num_atms ! my_img_lo, my_img_hi

#ifdef _OPENMP_
      tid = omp_get_thread_num()
      arr_ind = (tid*sd_num_atms) + 1
#endif
      ! Common translation (ie. no translation) flag is packed at
      ! the front of each sublist followed by the count(s) of sublist
      ! image pair entries.

      !common_tran = ipairs(ipairs_idx)
      !if(ipairs_idx .ne. start_index(img_i)) print *, ipairs_idx, start_index(img_i), mytaskid
      !ipairs_idx = ipairs_idx + 1
      ipairs_idx = start_index(img_i) + 1
    
      ! Electrostatic evaluation-only count followed by
      ! full evaluation count packed at the front of each pair sublist.

      ee_eval_cnt = ipairs(ipairs_idx)
      full_eval_cnt = ipairs(ipairs_idx + 1)
      ipairs_idx = ipairs_idx + 2
      if (ee_eval_cnt + full_eval_cnt .gt. 0) then
#ifdef pmemd_SPDP
        x_i = img_crd_q_sp(1, img_i)
        y_i = img_crd_q_sp(2, img_i)
        z_i = img_crd_q_sp(3, img_i)
#else
        x_i = img_crd(1, img_i)
        y_i = img_crd(2, img_i)
        z_i = img_crd(3, img_i)
#endif /*  pmemd_SPDP */
        !if (common_tran .eq. 0) then
          ! We need all the translation vectors:
       ! else
       !   ! Just put the x,y,z values in the middle cell
       !   x_tran(1, 13) = - x_i
       !   x_tran(2, 13) = - y_i
       !   x_tran(3, 13) = - z_i
       ! end if
        ! We always need virials from this routine if we need energies,
        ! because the virials are used in estimating the pme error.

        if (cutoffs_equal) then
#ifdef _OPENMP_
#ifdef pmemd_SPDP
          call pairs_calc_midpoint_efv(img_i,ee_eval_cnt, full_eval_cnt, &
                                evdw_stk,eed_stk, eedvir_stk, ehb_stk, &
                                vxx, vyy, vzz,vxy,vxz, vyz, &
                                img_frc(:,arr_ind), img_crd_q_sp, efs_tbl, &
                              eed_cub, typ_ico, ipairs(ipairs_idx), &
                              proc_iac, proc_gbl_cn1_sp, proc_gbl_cn2_sp, x_i, y_i, z_i)
#else /* pmemd_SPDP */
          call pairs_calc_midpoint_efv(img_i,ee_eval_cnt, full_eval_cnt, &
                                evdw_stk,eed_stk, eedvir_stk, ehb_stk, &
                                vxx, vyy, vzz,vxy,vxz, vyz, &
                                img_frc(:,arr_ind), img_crd, img_qterm, efs_tbl, &
                              eed_cub, typ_ico, ipairs(ipairs_idx), &
                              proc_iac, gbl_cn1, gbl_cn2, x_i, y_i, z_i)
#endif /* pmemd_SPDP */
#else
#ifdef pmemd_SPDP
          call pairs_calc_midpoint_efv(img_frc, img_crd_q_sp, efs_tbl, &
                              eed_cub, typ_ico, ipairs(ipairs_idx), &
                              proc_iac, proc_gbl_cn1_sp, proc_gbl_cn2_sp, x_i, y_i, z_i)
#else
          call pairs_calc_midpoint_efv(img_frc, img_crd, img_qterm, efs_tbl, &
                              eed_cub, typ_ico, ipairs(ipairs_idx), &
                              proc_iac, gbl_cn1, gbl_cn2, x_i, y_i, z_i)
#endif                               
#endif /* _OPENMP_ */
         end if

       !ipairs_idx = ipairs_idx + ee_eval_cnt + full_eval_cnt
  
      end if/* ee_eval_cnt + full_eval_cnt .gt. 0  */

    end do
    end do
#ifdef _OPENMP_
!$omp end parallel do 
#endif /* _OPENMP_ */

  else ! do not need energies...
  !ipairs_idx = 1

#ifdef _OPENMP_
!$omp parallel do default(shared) & 
!$omp& schedule (dynamic,16) &
!!$omp& shared(bkt_cnt, bkt_atm_cnt, bkt_atm_lst, proc_num_atms , proc_ghost_num_atms, start_index, ipairs, &
!!$omp& img_crd, img_qterm, img_frc, cutoffs_equal, tranvec, fs_tbl, &
!!$omp& eed_cub, typ_ico, proc_iac, gbl_cn1, gbl_cn2, need_virials,  &
!!$omp& sd_num_atms) &
!$omp& private(bkt, atm_i_idx, img_i, ipairs_idx, ee_eval_cnt, &
!$omp&  full_eval_cnt,  x_i, y_i, z_i, i, tid, arr_ind) &
!$omp&  reduction(+: vxx, vyy, vzz,vxy,vxz, vyz)
#endif /* _OPENMP_ */
    do bkt = 1,bkt_cnt
    do atm_i_idx = 1, bkt_atm_cnt(bkt)
      img_i = bkt_atm_lst(atm_i_idx, bkt)
!   do img_i = 1 , proc_num_atms + proc_ghost_num_atms !my_img_lo, my_img_hi

#ifdef _OPENMP_
      tid = omp_get_thread_num()
      arr_ind = tid*sd_num_atms + 1
#endif

      ! Common translation (ie. no translation) flag is packed at
      ! the front of each sublist followed by the count(s) of sublist
      ! image pair entries.

      !common_tran = ipairs(ipairs_idx)
      !ipairs_idx = ipairs_idx + 1
      ipairs_idx = start_index(img_i) + 1
    
      ! Electrostatic evaluation-only count followed by
      ! full evaluation count packed at the front of each pair sublist.

      ee_eval_cnt = ipairs(ipairs_idx)
      full_eval_cnt = ipairs(ipairs_idx + 1)
      ipairs_idx = ipairs_idx + 2

      if (ee_eval_cnt + full_eval_cnt .gt. 0) then
#ifdef pmemd_SPDP
        x_i = img_crd_q_sp(1, img_i)
        y_i = img_crd_q_sp(2, img_i)
        z_i = img_crd_q_sp(3, img_i)
#else        
        x_i = img_crd(1, img_i)
        y_i = img_crd(2, img_i)
        z_i = img_crd(3, img_i)
#endif /* pmemd_SPDP */
       ! if (common_tran .eq. 0) then
          ! We need all the translation vectors:
#if 0
          do i = 0, 17
            x_tran(1, i) = tranvec(1, i) - x_i
            x_tran(2, i) = tranvec(2, i) - y_i
            x_tran(3, i) = tranvec(3, i) - z_i
          end do
#endif
      !  else
      !    ! Just put the x,y,z values in the middle cell
      !    x_tran(1, 13) = - x_i
      !    x_tran(2, 13) = - y_i
      !    x_tran(3, 13) = - z_i
      !  end if

        if (need_virials) then

          if (cutoffs_equal) then
#ifdef _OPENMP_
#ifdef pmemd_SPDP
            call pairs_calc_midpoint_fv(img_i,ee_eval_cnt, full_eval_cnt, vxx, vyy, vzz,vxy,vxz, vyz, &
                                img_frc(:,arr_ind), img_crd_q_sp, fs_tbl, &
                               eed_cub, typ_ico, ipairs(ipairs_idx), &
                              proc_iac, proc_gbl_cn1_sp, proc_gbl_cn2_sp,  x_i, y_i, z_i)
#else
            call pairs_calc_midpoint_fv(img_i,ee_eval_cnt, full_eval_cnt, vxx, vyy, vzz,vxy,vxz, vyz, &
                                img_frc(:,arr_ind), img_crd, img_qterm, fs_tbl, &
                               eed_cub, typ_ico, ipairs(ipairs_idx), &
                              proc_iac, gbl_cn1, gbl_cn2,  x_i, y_i, z_i)
#endif /* pmemd_SPDP */
#else
#ifdef pmemd_SPDP
            call pairs_calc_midpoint_fv(img_frc, img_crd_q_sp, fs_tbl, &
                               eed_cub, typ_ico, ipairs(ipairs_idx), &
                              proc_iac, proc_gbl_cn1_sp, proc_gbl_cn2_sp,  x_i, y_i, z_i)
#else
            call pairs_calc_midpoint_fv(img_frc, img_crd, img_qterm, fs_tbl, &
                               eed_cub, typ_ico, ipairs(ipairs_idx), &
                              proc_iac, gbl_cn1, gbl_cn2,  x_i, y_i, z_i)
#endif                              
#endif /* _OPENMP_ */
          end if

        else

          if (cutoffs_equal) then
#ifdef _OPENMP_
#ifdef pmemd_SPDP
            call pairs_calc_midpoint_f(img_i,ee_eval_cnt, full_eval_cnt, &
                               img_frc(:,arr_ind), img_crd_q_sp, fs_tbl, &
                              eed_cub, typ_ico, ipairs(ipairs_idx), &
                              proc_iac, proc_gbl_cn1_sp, proc_gbl_cn2_sp,  x_i, y_i, z_i)
#else
            call pairs_calc_midpoint_f(img_i,ee_eval_cnt, full_eval_cnt, &
                               img_frc(:,arr_ind), img_crd, img_qterm, fs_tbl, &
                              eed_cub, typ_ico, ipairs(ipairs_idx), &
                              proc_iac, gbl_cn1, gbl_cn2,  x_i, y_i, z_i)
#endif /* pmemd_SPDP */
#else
#ifdef pmemd_SPDP
            call pairs_calc_midpoint_f(img_frc, img_crd_q_sp, fs_tbl, &
                              eed_cub, typ_ico, ipairs(ipairs_idx), &
                              proc_iac, proc_gbl_cn1_sp, proc_gbl_cn2_sp,  x_i, y_i, z_i)
#else
            call pairs_calc_midpoint_f(img_frc, img_crd, img_qterm, fs_tbl, &
                              eed_cub, typ_ico, ipairs(ipairs_idx), &
                              proc_iac, gbl_cn1, gbl_cn2,  x_i, y_i, z_i)
#endif                              
#endif /* _OPENMP_ */
          end if
        end if/*need_virials*/

       !ipairs_idx = ipairs_idx + ee_eval_cnt + full_eval_cnt

      end if /*ee_eval_cnt + full_eval_cnt .gt. 0 */

    end do
    end do

#ifdef _OPENMP_
!$omp end parallel do 
#endif /* _OPENMP_ */
  end if /*need_pot_enes*/
  ! Save the energies:
                                                                                
  eedvir = eedvir_stk
  eed = eed_stk
  evdw = evdw_stk
  ehb = ehb_stk

  ! Save the virials.

  virial(1, 1) = vxx
!  virial(1, 2) = vxy
!  virial(2, 1) = vxy
!  virial(1, 3) = vxz
!  virial(3, 1) = vxz
  virial(2, 2) = vyy
!  virial(2, 3) = vyz
!  virial(3, 2) = vyz
  virial(3, 3) = vzz
  return

contains

#define BUILD_PAIRS_CALC_EFV
#define NEED_ENE
#define NEED_VIR
#include "pairs_calc_midpoint.i"
#undef NEED_VIR
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_EFV

#define BUILD_PAIRS_CALC_FV
#define NEED_VIR
#include "pairs_calc_midpoint.i"
#undef NEED_VIR
#undef BUILD_PAIRS_CALC_FV

#define BUILD_PAIRS_CALC_F
#include "pairs_calc_midpoint.i"
#undef BUILD_PAIRS_CALC_F

end subroutine get_nb_energy_midpoint
