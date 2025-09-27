!*******************************************************************************
!
! Subroutine:  get_nb_ips_energy
!
! Description:
!              
! The main routine for non bond energy (vdw and hbond) as well as direct part
! of ewald sum.  It is structured for parallelism.
!
!*******************************************************************************

subroutine get_nb_ips_energy(img_frc, img_crd, img_qterm, eed_cub, &
                         ipairs, tranvec, need_pot_enes, need_virials, &
                         eed, evdw, ehb, eedvir, virial)
  use img_mod
  use timers_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use ene_frc_splines_mod
  use nbips_mod

  implicit none

! Formal arguments:

  double precision, intent(in out) :: img_frc(3, *)
  double precision, intent(in)     :: img_crd(3, *)
  double precision , intent(in)    :: img_qterm(*)
  double precision, intent(in)     :: eed_cub(*)
  integer                          :: ipairs(*)
  double precision, intent(in)     :: tranvec(1:3, 0:17)
  logical, intent(in)              :: need_pot_enes
  logical, intent(in)              :: need_virials
  double precision, intent(out)    :: eed
  double precision, intent(out)    :: evdw
  double precision, intent(out)    :: ehb
  double precision, intent(out)    :: eedvir
  double precision, intent(out)    :: virial(3, 3)

! Local variables and parameters:

  double precision      del
  double precision      dxdr
  double precision      eedtbdns_stk
  double precision      eedvir_stk, eed_stk, evdw_stk, ehb_stk
  double precision      max_nb_cut2, es_cut2, es_cut
  double precision      x_i, y_i, z_i
  double precision      x_tran(1:3, 0:17)
  double precision      vxx, vxy, vxz, vyy, vyz, vzz
  integer               i
  integer               ipairs_idx
  integer               ntypes_stk
  integer               img_i
  integer               ee_eval_cnt
  integer               full_eval_cnt
  integer               common_tran    ! flag - 1 if translation not needed
  logical               cutoffs_equal

#ifdef MPI
  if (my_img_lo .gt. my_img_hi) return
#endif /* MPI */

#ifdef CUDA
  call gpu_get_nb_energy()
#else

  ntypes_stk = ntypes

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
  vxy = 0.d0
  vxz = 0.d0
  vyy = 0.d0
  vyz = 0.d0
  vzz = 0.d0

  ipairs_idx = 1

  if (need_pot_enes) then

    do img_i = my_img_lo, my_img_hi

      ! Common translation (ie. no translation) flag is packed at
      ! the front of each sublist followed by the count(s) of sublist
      ! image pair entries.

      common_tran = ipairs(ipairs_idx)
      ipairs_idx = ipairs_idx + 1
    
      ! Electrostatic evaluation-only count followed by
      ! full evaluation count packed at the front of each pair sublist.

      ee_eval_cnt = ipairs(ipairs_idx)
      full_eval_cnt = ipairs(ipairs_idx + 1)
      ipairs_idx = ipairs_idx + 2

      if (ee_eval_cnt + full_eval_cnt .gt. 0) then

        x_i = img_crd(1, img_i)
        y_i = img_crd(2, img_i)
        z_i = img_crd(3, img_i)

        if (common_tran .eq. 0) then
          ! We need all the translation vectors:
          do i = 0, 17
            x_tran(1, i) = tranvec(1, i) - x_i
            x_tran(2, i) = tranvec(2, i) - y_i
            x_tran(3, i) = tranvec(3, i) - z_i
          end do
        else
          ! Just put the x,y,z values in the middle cell
          x_tran(1, 13) = - x_i
          x_tran(2, 13) = - y_i
          x_tran(3, 13) = - z_i
        end if

        ! We always need virials from this routine if we need energies,
        ! because the virials are used in estimating the pme error.

        if (cutoffs_equal) then
          call pairs_calc_efv(img_frc, img_crd, img_qterm, efs_tbl, &
                              eed_cub, typ_ico, ipairs(ipairs_idx), &
                              gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
        else
          call pairs_calc_efv_2cut(img_frc, img_crd, img_qterm, efs_tbl, &
                                   eed_cub, typ_ico, ipairs(ipairs_idx), &
                                   gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
        end if
        ipairs_idx = ipairs_idx + ee_eval_cnt + full_eval_cnt
      end if
    end do

  else ! do not need energies...

    do img_i = my_img_lo, my_img_hi

      ! Common translation (ie. no translation) flag is packed at
      ! the front of each sublist followed by the count(s) of sublist
      ! image pair entries.

      common_tran = ipairs(ipairs_idx)
      ipairs_idx = ipairs_idx + 1
    
      ! Electrostatic evaluation-only count followed by
      ! full evaluation count packed at the front of each pair sublist.

      ee_eval_cnt = ipairs(ipairs_idx)
      full_eval_cnt = ipairs(ipairs_idx + 1)
      ipairs_idx = ipairs_idx + 2

      if (ee_eval_cnt + full_eval_cnt .gt. 0) then

        x_i = img_crd(1, img_i)
        y_i = img_crd(2, img_i)
        z_i = img_crd(3, img_i)

        if (common_tran .eq. 0) then
          ! We need all the translation vectors:
          do i = 0, 17
            x_tran(1, i) = tranvec(1, i) - x_i
            x_tran(2, i) = tranvec(2, i) - y_i
            x_tran(3, i) = tranvec(3, i) - z_i
          end do
        else
          ! Just put the x,y,z values in the middle cell
          x_tran(1, 13) = - x_i
          x_tran(2, 13) = - y_i
          x_tran(3, 13) = - z_i
        end if

        if (need_virials) then
          if (cutoffs_equal) then
            call pairs_calc_fv(img_frc, img_crd, img_qterm, fs_tbl, &
                               eed_cub, typ_ico, ipairs(ipairs_idx), &
                               gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
          else
            call pairs_calc_fv_2cut(img_frc, img_crd, img_qterm, fs_tbl, &
                                    eed_cub, typ_ico, ipairs(ipairs_idx), &
                                    gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
          end if
        else
          if (cutoffs_equal) then
            call pairs_calc_f(img_frc, img_crd, img_qterm, fs_tbl, &
                              eed_cub, typ_ico, ipairs(ipairs_idx), &
                              gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
          else
            call pairs_calc_f_2cut(img_frc, img_crd, img_qterm, fs_tbl, &
                                   eed_cub, typ_ico, ipairs(ipairs_idx), &
                                   gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
          end if
        end if
        ipairs_idx = ipairs_idx + ee_eval_cnt + full_eval_cnt
      end if
    end do

  end if

  ! Save the energies:
                                                                                
  eedvir = eedvir_stk
  eed = eed_stk
  evdw = evdw_stk
  ehb = ehb_stk

  ! Save the virials.

  virial(1, 1) = vxx
  virial(1, 2) = vxy
  virial(2, 1) = vxy
  virial(1, 3) = vxz
  virial(3, 1) = vxz
  virial(2, 2) = vyy
  virial(2, 3) = vyz
  virial(3, 2) = vyz
  virial(3, 3) = vzz
#endif /* CUDA */
  return

contains

!IPS CODE: GET_NB_ENE_IPS
#define BUILD_PAIRS_CALC_EFV
#define NEED_ENE
#define NEED_VIR
#include "pairs_calc_ips.i"
#undef NEED_VIR
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_EFV

#define BUILD_PAIRS_CALC_FV
#define NEED_VIR
#include "pairs_calc_ips.i"
#undef NEED_VIR
#undef BUILD_PAIRS_CALC_FV

#define BUILD_PAIRS_CALC_F
#include "pairs_calc_ips.i"
#undef BUILD_PAIRS_CALC_F

#define BUILD_PAIRS_CALC_EFV
#define NEED_ENE
#define NEED_VIR
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#include "pairs_calc_ips.i"
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef NEED_VIR
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_EFV

#define BUILD_PAIRS_CALC_FV
#define NEED_VIR
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#include "pairs_calc_ips.i"
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_FV

#define BUILD_PAIRS_CALC_F
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#include "pairs_calc_ips.i"
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef BUILD_PAIRS_CALC_F

end subroutine get_nb_ips_energy
