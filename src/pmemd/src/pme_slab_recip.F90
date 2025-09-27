#include "copyright.i"

!*******************************************************************************
!
! Module: pme_slab_recip_mod
!
! Description: <TBS>
!
! NOTE NOTE NOTE: This code assumes frc_int .eq. 0 and should only be used
!                 under these conditions!!!
!*******************************************************************************

module pme_slab_recip_mod

  use pme_recip_dat_mod

!precision setting for SPDP, DPDP, SPSP
#include "include_precision.i"
  implicit none

  private       ! Everything here is private unless specified public

#ifdef MPI
  public        :: pme_slab_recip_setup, &
                   do_slab_pmesh_kspace, &
                   claim_slab_recip_imgs, &
                   claim_slab_recip_imgs_nonorthog
#else
  public        :: pme_slab_recip_setup, &
                   do_slab_pmesh_kspace
#endif /* MPI */
contains

!*******************************************************************************
!
! Subroutine:  pme_slab_recip_setup
!
! Description:
!              
!*******************************************************************************

subroutine pme_slab_recip_setup(num_ints, num_reals)

  use gbl_constants_mod
  use loadbal_mod
  use pmemd_lib_mod
  use pme_slab_fft_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

#ifdef MPI
  if (nrespa .eq. 1 .and. imin .eq. 0 .and. numtasks * bspl_order * 2 .gt. nfft3) then
    fft_redist_enabled = .true.
    fft_redist_needed = .true.
  end if

  call slab_fft_setup(recip_numtasks, fft_redist_enabled, &
                      num_ints, num_reals)

#else
  call slab_fft_setup(num_ints, num_reals)
#endif

  return

end subroutine pme_slab_recip_setup

!*******************************************************************************
!
! Subroutine:  do_slab_pmesh_kspace
!
! Description:  <TBS>
!
! INPUT:
!
! OUTPUT:
!
! img_frc:      Forces incremented by k-space sum.
! virial:       Virial due to k-space sum (valid for atomic scaling;
!               rigid molecule virial needs a correction term not computed here.
!*******************************************************************************

subroutine do_slab_pmesh_kspace(img_frc, eer, virial, need_pot_enes, &
                                need_virials)

  use img_mod
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use pbc_mod
  use pme_slab_fft_mod
  use prmtop_dat_mod
  use timers_mod

  implicit none

! Formal arguments:
  double precision      :: img_frc(3, *)
  double precision      :: eer
  double precision      :: virial(3, 3)
  logical               :: need_pot_enes
  logical               :: need_virials

! Local variables:

! NOTE - For FFTW, *_q need to be 16 byte aligned if they are moved off
!        the stack (and the stack must be aligned if they are on the stack).

#ifdef MPI
  integer               :: my_chg_cnt
#endif

  double precision      :: theta(bspl_order * 3 * gbl_used_recip_img_cnt)
  double precision      :: dtheta(bspl_order * 3 * gbl_used_recip_img_cnt)
#ifdef MPI
  double precision      :: xyz_q(max(xy_slab_dbl_cnt * max_xy_slab_cnt, &
                                     zx_slab_dbl_cnt * max_zx_slab_cnt))
  double precision      :: zxy_q(max(xy_slab_dbl_cnt * max_xy_slab_cnt, &
                                     zx_slab_dbl_cnt * max_zx_slab_cnt))
  integer               :: my_chgs(gbl_used_recip_img_cnt)
#else /*MPI*/
  double precision      :: xyz_q(2 * fft_x_dim * fft_y_dim * fft_z_dim)
  double precision      :: zxy_q(2 * fft_x_dim * fft_y_dim * fft_z_dim)
#endif /*MPI*/
  integer               :: ifracts(3, gbl_used_recip_img_cnt)

  ! Exponential tables are only practical to use for orthogonal unit cells, and
  ! only need updating when constant pressure runs are being done.
 
  if (ntp .gt. 0 .and. is_orthog .ne. 0) call load_m_exp_tbls

#ifdef MPI

  if (is_orthog .ne. 0) then
    call get_grid_weights(natom, gbl_img_crd, ifracts, theta, dtheta, &
                          bspl_order, gbl_used_recip_img_cnt, &
                          gbl_used_recip_img_lst, my_chgs, my_chg_cnt)
  else
    call get_grid_weights_nonorthog(natom, atm_crd, gbl_img_atm_map, &
                                    ifracts, theta, dtheta, bspl_order, &
                                    gbl_used_recip_img_cnt, &
                                    gbl_used_recip_img_lst, &
                                    my_chgs, my_chg_cnt)
  end if
#else
  if (is_orthog .ne. 0) then
    call get_grid_weights(natom, gbl_img_crd, ifracts, theta, dtheta, &
                          bspl_order)
  else
    call get_grid_weights_nonorthog(natom, atm_crd, gbl_img_atm_map, ifracts, &
                                    theta, dtheta, bspl_order)
  end if
#endif

  call update_pme_time(bspline_timer)

! Fill Charge Grid.  Charges are approximated on an even grid.

#ifdef MPI
  call fill_charge_grid(xyz_q, theta, bspl_order, gbl_img_qterm, ifracts, &
                        my_chgs, my_chg_cnt)
#else /*MPI*/
  call fill_charge_grid(xyz_q, theta, bspl_order, gbl_img_qterm, ifracts)
#endif /*MPI*/

  call update_pme_time(grid_charges_timer)

  call slab_fft3drc_forward(xyz_q, zxy_q, fft_x_dim, fft_y_dim, fft_z_dim)

  call update_pme_time(fft_timer)

  if (is_orthog .ne. 0) then
    if (need_virials) then
      call scalar_sumrc_uv(zxy_q, ew_coeff, uc_volume, &
                           gbl_prefac1, gbl_prefac2, gbl_prefac3, &
                           nfft1, nfft2, nfft3, &
                           fft_x_dim, fft_y_dim, fft_z_dim, eer, virial)
    else if (need_pot_enes) then
!      call scalar_sumrc_u(zxy_q, ew_coeff, uc_volume, &
!                          gbl_prefac1, gbl_prefac2, gbl_prefac3, &
!                          nfft1, nfft2, nfft3, &
!                          fft_x_dim, fft_y_dim, fft_z_dim, eer)
! Ross Walker - 7/8/2014 - bug bug - The PME Error estimate requires the vir
! since it compares virial vs energy. 
! pme_err_est =   vir_vs_ene = vir%elec_recip(1, 1) + &
!               vir%elec_recip(2, 2) + &
!               vir%elec_recip(3, 3) + &
!               vir%eedvir + &
!               vir%elec_nb_adjust(1, 1) + &
!               vir%elec_nb_adjust(2, 2) + &
!               vir%elec_nb_adjust(3, 3)
! So we can't just naively skip it here unless we choose not to print the 
! Ewald error estimate anymore (like in the GPU code). For now we will calculate
! the reciprocal virial whenever need_pot_enes is true.
      call scalar_sumrc_uv(zxy_q, ew_coeff, uc_volume, &
                           gbl_prefac1, gbl_prefac2, gbl_prefac3, &
                           nfft1, nfft2, nfft3, &
                           fft_x_dim, fft_y_dim, fft_z_dim, eer, virial)
    else
      call scalar_sumrc(zxy_q, ew_coeff, uc_volume, &
                        gbl_prefac1, gbl_prefac2, gbl_prefac3, &
                        nfft1, nfft2, nfft3, &
                        fft_x_dim, fft_y_dim, fft_z_dim)
    end if
  else
    if (need_virials) then
      call scalar_sumrc_nonorthog_uv(zxy_q, ew_coeff, uc_volume, &
                                     gbl_prefac1, gbl_prefac2, gbl_prefac3, &
                                     nfft1, nfft2, nfft3, &
                                     fft_x_dim, fft_y_dim, fft_z_dim, &
                                     eer, virial)
    else if (need_pot_enes) then
!      call scalar_sumrc_nonorthog_u(zxy_q, ew_coeff, uc_volume, &
!                                    gbl_prefac1, gbl_prefac2, gbl_prefac3, &
!                                    nfft1, nfft2, nfft3, &
!                                    fft_x_dim, fft_y_dim, fft_z_dim, &
!                                    eer)
! Ross Walker - 7/8/2014 - bug bug - The PME Error estimate requires the vir
! since it compares virial vs energy. 
! pme_err_est =   vir_vs_ene = vir%elec_recip(1, 1) + &
!               vir%elec_recip(2, 2) + &
!               vir%elec_recip(3, 3) + &
!               vir%eedvir + &
!               vir%elec_nb_adjust(1, 1) + &
!               vir%elec_nb_adjust(2, 2) + &
!               vir%elec_nb_adjust(3, 3)
! So we can't just naively skip it here unless we choose not to print the 
! Ewald error estimate anymore (like in the GPU code). For now we will calculate
! the reciprocal virial whenever need_pot_enes is true.
      call scalar_sumrc_nonorthog_uv(zxy_q, ew_coeff, uc_volume, &
                                     gbl_prefac1, gbl_prefac2, gbl_prefac3, &
                                     nfft1, nfft2, nfft3, &
                                     fft_x_dim, fft_y_dim, fft_z_dim, &
                                     eer, virial)
    else
      call scalar_sumrc_nonorthog(zxy_q, ew_coeff, uc_volume, &
                                  gbl_prefac1, gbl_prefac2, gbl_prefac3, &
                                  nfft1, nfft2, nfft3, &
                                  fft_x_dim, fft_y_dim, fft_z_dim)
    end if
  end if

  call update_pme_time(scalar_sum_timer)

  call slab_fft3drc_back(zxy_q, xyz_q, fft_x_dim, fft_y_dim, fft_z_dim)

  call update_pme_time(fft_timer)

#ifdef MPI

  call grad_sum(xyz_q, theta, dtheta, bspl_order, gbl_img_qterm, &
                img_frc, ifracts, my_chgs, my_chg_cnt, &
                fft_x_dim, fft_y_dim, fft_z_dim)

#else /*MPI*/
  call grad_sum(xyz_q, theta, dtheta, bspl_order, gbl_img_qterm, &
                img_frc, ifracts, fft_x_dim, fft_y_dim, fft_z_dim)
#endif /* MPI */

  call update_pme_time(grad_sum_timer)

  return

end subroutine do_slab_pmesh_kspace

!*******************************************************************************
!
! Subroutine:  fill_charge_grid
!
! Description: <TBS>
!
! INPUT:
!
! theta1, theta2, theta3:       Spline coeff arrays.
! ifracts:                      int(scaled fractional coords).
! nfft1, nfft2, nfft3:          Logical charge grid dimensions.
!
! fft_x_dim, fft_y_dim, fft_z_dim: Physical charge grid dims.
!
! order:                    Order of spline interpolation.
!
! OUTPUT:
!
! q:                            Charge grid.
!              
!*******************************************************************************

#ifdef MPI
subroutine fill_charge_grid(q, theta, order, img_qterm, ifracts, &
                            my_chgs, my_chg_cnt)
#else
subroutine fill_charge_grid(q, theta, order, img_qterm, ifracts)
#endif

  use img_mod
  use mdin_ewald_dat_mod
  use pme_slab_fft_mod
  use prmtop_dat_mod
  use ti_mod
#ifdef TIDECOMP
  use ti_decomp_mod
  use mdin_ctrl_dat_mod
#endif

  implicit none

! Formal arguments:

  integer               :: order
  double precision      :: theta(order, 3, *)
  double precision      :: img_qterm(*)
  integer               :: ifracts(3, *)
  double precision      :: q(2 * fft_x_dim, fft_y_dim, my_xy_slab_cnt)

#ifdef MPI
  integer               :: my_chgs(*), my_chg_cnt
#endif

! Local variables:

  integer               :: i, j, k
  integer               :: i_base, j_base, k_base
  integer               :: ith1, ith2, ith3
  integer               :: my_chg_idx
#ifdef TIDECOMP
  integer               :: atmidx
#endif
  double precision      :: charge
  double precision      :: k_term, j_term
#ifdef MPI
! integer               :: kbot0
  integer               :: kbot, ktop
#endif
  integer               :: i_tbl(-order : nfft1)
  integer               :: j_tbl(-order : nfft2)
  integer               :: k_tbl(-order : nfft3)
  ! For TI:
  integer               :: iatom
#ifdef MPI
#ifdef RANGE_DBG
  integer               :: dbg_img_map(natom)
  integer               :: dbg_img_map_res(natom)

  dbg_img_map(:) = 0
#endif /* RANGE_DBG */
#endif

#ifdef MPI
! Old usages; kept for documentation during code conversion...
! kbot0 = my_xy_slab_start
! kbot = kbot0 + 1
! ktop = kbot0 + my_xy_slab_cnt
#endif

! Zero the Charge grids:

#ifdef MPI
  call zero_q_array(q, 2 * fft_x_dim * fft_y_dim * my_xy_slab_cnt)
#else
  call zero_q_array(q, 2 * fft_x_dim * fft_y_dim * fft_z_dim)
#endif

  ! Initialize the indexing tables.  It actually produces faster code to
  ! do this here rather than caching the info remotely.

  do i = 0, nfft1
    i_tbl(i) = i + 1
  end do
  do i = -order, -1
    i_tbl(i) = i + nfft1 + 1
  end do

  do j = 0, nfft2
    j_tbl(j) = j + 1
  end do
  do j = -order, -1
    j_tbl(j) = j + nfft2 + 1
  end do

#ifdef MPI
  kbot = 1
  ktop = my_xy_slab_cnt
  do k = 0, nfft3
    k_tbl(k) = k - my_xy_slab_start + 1
    if (k_tbl(k) .lt. kbot .or. k_tbl(k) .gt. ktop) k_tbl(k) = 0
  end do
  do k = -order, -1
    k_tbl(k) = k + nfft3 - my_xy_slab_start + 1
    if (k_tbl(k) .lt. kbot .or. k_tbl(k) .gt. ktop) k_tbl(k) = 0
  end do
#else
  do k = 0, nfft3
    k_tbl(k) = k + 1
  end do
  do k = -order, -1
    k_tbl(k) = k + nfft3 + 1
  end do
#endif

  ! We special-case order 4, the default.

  if (order .eq. 4) then
#ifdef MPI
    do my_chg_idx = 1, my_chg_cnt
      if (ti_mode .eq. 1 .or. ti_mode .eq. 2) then
        !Determine atom # that we are dealing with 
        iatom = gbl_img_atm_map(my_chgs(my_chg_idx)) ! get atm number from idx
        if(ti_lst(ti_mask_piece,iatom) .ne. 0) cycle
      end if
      charge = img_qterm(my_chgs(my_chg_idx))
#else
    do my_chg_idx = 1, natom
      if (ti_mode .eq. 1 .or. ti_mode .eq. 2) then
        !Determine atom # that we are dealing with 
        iatom = gbl_img_atm_map(my_chg_idx) ! get atm number from idx
        if(ti_lst(ti_mask_piece,iatom) .ne. 0) cycle
      end if
      charge = img_qterm(my_chg_idx)
#endif
#ifdef TIDECOMP
      atmidx =  gbl_img_atm_map(my_chg_idx)     
#endif

      i_base = ifracts(1, my_chg_idx)
      j_base = ifracts(2, my_chg_idx)
      k_base = ifracts(3, my_chg_idx)

      do ith3 = 1, 4

        k = k_tbl(k_base + ith3)

#ifdef MPI
        if (k .ne. 0) then

#ifdef RANGE_DBG
          dbg_img_map(my_chgs(my_chg_idx)) = dbg_img_map(my_chgs(my_chg_idx))+1
#endif /* RANGE_DBG */
#endif
          k_term =  theta(ith3, 3, my_chg_idx) * charge
          do ith2 = 1, 4
            j = j_tbl(j_base + ith2)
            j_term = k_term * theta(ith2, 2, my_chg_idx)
            q(i_tbl(i_base+1), j, k) = q(i_tbl(i_base+1), j, k) + &
                                       theta(1, 1, my_chg_idx) * j_term
            q(i_tbl(i_base+2), j, k) = q(i_tbl(i_base+2), j, k) + &
                                       theta(2, 1, my_chg_idx) * j_term
            q(i_tbl(i_base+3), j, k) = q(i_tbl(i_base+3), j, k) + &
                                       theta(3, 1, my_chg_idx) * j_term
            q(i_tbl(i_base+4), j, k) = q(i_tbl(i_base+4), j, k) + &
                                       theta(4, 1, my_chg_idx) * j_term
#ifdef TIDECOMP
            if(ntwd .gt. 0) then
                pme_cell_cnt(i_tbl(i_base+1),j,k,1) = pme_cell_cnt(i_tbl(i_base+1),j,k,1)+1
                pme_cell_cnt(i_tbl(i_base+2),j,k,1) = pme_cell_cnt(i_tbl(i_base+2),j,k,1)+1
                pme_cell_cnt(i_tbl(i_base+3),j,k,1) = pme_cell_cnt(i_tbl(i_base+3),j,k,1)+1
                pme_cell_cnt(i_tbl(i_base+4),j,k,1) = pme_cell_cnt(i_tbl(i_base+4),j,k,1)+1
                pme_cell_cnt(i_tbl(i_base+1),j,k,2) = pme_cell_cnt(i_tbl(i_base+1),j,k,2)+order*3-ith3-ith2
                pme_cell_cnt(i_tbl(i_base+2),j,k,2) = pme_cell_cnt(i_tbl(i_base+2),j,k,2)+order*3-ith3-ith2-1
                pme_cell_cnt(i_tbl(i_base+3),j,k,2) = pme_cell_cnt(i_tbl(i_base+3),j,k,2)+order*3-ith3-ith2-2
                pme_cell_cnt(i_tbl(i_base+4),j,k,2) = pme_cell_cnt(i_tbl(i_base+4),j,k,2)+order*3-ith3-ith2-3
                pme_cell_atom(i_tbl(i_base+1),j,k,pme_cell_cnt(i_tbl(i_base+1),j,k,1)) = atmidx
                pme_cell_atom(i_tbl(i_base+2),j,k,pme_cell_cnt(i_tbl(i_base+2),j,k,1)) = atmidx
                pme_cell_atom(i_tbl(i_base+3),j,k,pme_cell_cnt(i_tbl(i_base+3),j,k,1)) = atmidx
                pme_cell_atom(i_tbl(i_base+4),j,k,pme_cell_cnt(i_tbl(i_base+4),j,k,1)) = atmidx
                pme_cell_weight(i_tbl(i_base+1),j,k,pme_cell_cnt(i_tbl(i_base+1),j,k,1)) = order*3-ith3-ith2
                pme_cell_weight(i_tbl(i_base+2),j,k,pme_cell_cnt(i_tbl(i_base+2),j,k,1)) = order*3-ith3-ith2-1
                pme_cell_weight(i_tbl(i_base+3),j,k,pme_cell_cnt(i_tbl(i_base+3),j,k,1)) = order*3-ith3-ith2-2
                pme_cell_weight(i_tbl(i_base+4),j,k,pme_cell_cnt(i_tbl(i_base+4),j,k,1)) = order*3-ith3-ith2-3
            end if
#endif
          end do
#ifdef MPI
        end if
#endif
      end do
    end do
  else if (order .eq. 6) then
#ifdef MPI
    do my_chg_idx = 1, my_chg_cnt
      if (ti_mode .eq. 1 .or. ti_mode .eq. 2) then
        !Determine atom # that we are dealing with 
        iatom = gbl_img_atm_map(my_chgs(my_chg_idx)) ! get atm number from idx
        if(ti_lst(ti_mask_piece,iatom) .ne. 0) cycle
      end if
      charge = img_qterm(my_chgs(my_chg_idx))
#else
    do my_chg_idx = 1, natom
      if (ti_mode .eq. 1 .or. ti_mode .eq. 2) then
        !Determine atom # that we are dealing with 
        iatom = gbl_img_atm_map(my_chg_idx) ! get atm number from idx
        if(ti_lst(ti_mask_piece,iatom) .ne. 0) cycle
      end if
      charge = img_qterm(my_chg_idx)
#endif

      i_base = ifracts(1, my_chg_idx)
      j_base = ifracts(2, my_chg_idx)
      k_base = ifracts(3, my_chg_idx)

      do ith3 = 1, 6

        k = k_tbl(k_base + ith3)

#ifdef MPI
        if (k .ne. 0) then

#ifdef RANGE_DBG
          dbg_img_map(my_chgs(my_chg_idx)) = dbg_img_map(my_chgs(my_chg_idx))+1
#endif /* RANGE_DBG */
#endif
          k_term =  theta(ith3, 3, my_chg_idx) * charge
          do ith2 = 1, 6
            j = j_tbl(j_base + ith2)
            j_term = k_term * theta(ith2, 2, my_chg_idx)
            q(i_tbl(i_base+1), j, k) = q(i_tbl(i_base+1), j, k) + &
                                       theta(1, 1, my_chg_idx) * j_term
            q(i_tbl(i_base+2), j, k) = q(i_tbl(i_base+2), j, k) + &
                                       theta(2, 1, my_chg_idx) * j_term
            q(i_tbl(i_base+3), j, k) = q(i_tbl(i_base+3), j, k) + &
                                       theta(3, 1, my_chg_idx) * j_term
            q(i_tbl(i_base+4), j, k) = q(i_tbl(i_base+4), j, k) + &
                                       theta(4, 1, my_chg_idx) * j_term
            q(i_tbl(i_base+5), j, k) = q(i_tbl(i_base+5), j, k) + &
                                       theta(5, 1, my_chg_idx) * j_term
            q(i_tbl(i_base+6), j, k) = q(i_tbl(i_base+6), j, k) + &
                                       theta(6, 1, my_chg_idx) * j_term
          end do
#ifdef MPI
        end if
#endif
      end do
    end do
  else
#ifdef MPI
    do my_chg_idx = 1, my_chg_cnt
      if (ti_mode .eq. 1 .or. ti_mode .eq. 2) then
        !Determine atom # that we are dealing with 
        iatom = gbl_img_atm_map(my_chgs(my_chg_idx)) ! get atm number from idx
        if(ti_lst(ti_mask_piece,iatom) .ne. 0) cycle
      end if
      charge = img_qterm(my_chgs(my_chg_idx))
#else
    do my_chg_idx = 1, natom
      if (ti_mode .eq. 1 .or. ti_mode .eq. 2) then
        !Determine atom # that we are dealing with 
        iatom = gbl_img_atm_map(my_chg_idx) ! get atm number from idx
        if(ti_lst(ti_mask_piece,iatom) .ne. 0) cycle
      end if
      charge = img_qterm(my_chg_idx)
#endif

      i_base = ifracts(1, my_chg_idx)
      j_base = ifracts(2, my_chg_idx)
      k_base = ifracts(3, my_chg_idx)

      do ith3 = 1, order

        k = k_tbl(k_base + ith3)

#ifdef MPI
        if (k .ne. 0) then

#ifdef RANGE_DBG
          dbg_img_map(my_chgs(my_chg_idx)) = dbg_img_map(my_chgs(my_chg_idx))+1
#endif /* RANGE_DBG */
#endif
          k_term =  theta(ith3, 3, my_chg_idx) * charge
          do ith2 = 1, order
            j = j_tbl(j_base + ith2)
            j_term = k_term * theta(ith2, 2, my_chg_idx)
            do ith1 = 1, order
              i = i_tbl(i_base + ith1)
              q(i, j, k) = q(i, j, k) + theta(ith1, 1, my_chg_idx) * j_term
            end do
          end do
#ifdef MPI
        end if
#endif
      end do
    end do
  end if

#ifdef RANGE_DBG
  call mpi_reduce(dbg_img_map, dbg_img_map_res, size(dbg_img_map), &
                  mpi_integer, mpi_sum, 0, pmemd_comm, err_code_mpi)
  if (master) then
    do i = 1, natom
      if (dbg_img_map_res(i) .ne. order) then
        write(0,*)'DBG: IMAGE ', i, 'processed ', &
                  dbg_img_map_res(i), 'times in fill_charge_grid code!!!'
      end if
    end do
  end if
#endif /* RANGE_DBG */

  return

end subroutine fill_charge_grid

!*******************************************************************************
!
! Subroutine:  grad_sum
!
! Description: <TBS>
!              
!*******************************************************************************

#ifdef MPI

subroutine grad_sum(q, theta, dtheta, order, img_qterm, img_frc, ifracts, &
                    my_chgs, my_chg_cnt, x_dim, y_dim, z_dim)

#else
subroutine grad_sum(q, theta, dtheta, order, img_qterm, img_frc, ifracts, &
                    x_dim, y_dim, z_dim)
#endif

  use img_mod
  use mdin_ewald_dat_mod
  use prmtop_dat_mod
  use pbc_mod
  use pme_slab_fft_mod
  use ti_mod

  implicit none

! Formal arguments:

  integer               :: x_dim, y_dim, z_dim
  double precision      :: q(2 * x_dim, y_dim, z_dim)
  integer               :: order
  double precision      :: theta(order, 3, *)
  double precision      :: dtheta(order, 3, *)
  double precision      :: img_qterm(*)
  double precision      :: img_frc(3, *)
  integer               :: ifracts(3, *)

#ifdef MPI
  integer               :: my_chgs(*), my_chg_cnt
#endif

! Local variables:

  integer               :: i, j, k
  integer               :: i_base, j_base, k_base
  integer               :: ith1, ith2, ith3
  integer               :: my_chg_idx
  double precision      :: f1, f2, f3
  double precision      :: charge, qterm
  double precision      :: f1_term, f2_term, f3_term
  double precision      :: dfx, dfy, dfz
  double precision      :: recip_11, recip_22, recip_33
  double precision      :: dnfft1, dnfft2, dnfft3
#ifdef MPI
  integer               :: img_i
! integer               :: kbot0
  integer               :: kbot, ktop
#endif
  integer               :: i_tbl(-order : nfft1)
  integer               :: j_tbl(-order : nfft2)
  integer               :: k_tbl(-order : nfft3)
  ! For TI:
  integer               :: iatom
#ifdef MPI
#ifdef RANGE_DBG
  integer               :: dbg_img_map(natom)
  integer               :: dbg_img_map_res(natom)

  dbg_img_map(:) = 0
#endif /* RANGE_DBG */
#endif


#ifdef MPI
! Old usages; kept for documentation during code conversion...
! kbot0 = my_xy_slab_start
! kbot = kbot0 + 1
! ktop = kbot0 + my_xy_slab_cnt
#endif

  ! Initialize the indexing tables.  It actually produces faster code to
  ! do this here rather than caching the info remotely.

  do i = 0, nfft1
    i_tbl(i) = i + 1
  end do
  do i = -order, -1
    i_tbl(i) = i + nfft1 + 1
  end do

  do j = 0, nfft2
    j_tbl(j) = j + 1
  end do
  do j = -order, -1
    j_tbl(j) = j + nfft2 + 1
  end do

#ifdef MPI
  kbot = 1
  ktop = my_xy_slab_cnt
  do k = 0, nfft3
    k_tbl(k) = k - my_xy_slab_start + 1
    if (k_tbl(k) .lt. kbot .or. k_tbl(k) .gt. ktop) k_tbl(k) = 0
  end do
  do k = -order, -1
    k_tbl(k) = k + nfft3 - my_xy_slab_start + 1
    if (k_tbl(k) .lt. kbot .or. k_tbl(k) .gt. ktop) k_tbl(k) = 0
  end do
#else
  do k = 0, nfft3
    k_tbl(k) = k + 1
  end do
  do k = -order, -1
    k_tbl(k) = k + nfft3 + 1
  end do
#endif

  recip_11 = recip(1, 1)
  recip_22 = recip(2, 2)
  recip_33 = recip(3, 3)

  dnfft1 = dble(nfft1)
  dnfft2 = dble(nfft2)
  dnfft3 = dble(nfft3)

#ifdef MPI
  do my_chg_idx = 1, my_chg_cnt
    img_i = my_chgs(my_chg_idx)
    if (ti_mode .eq. 1 .or. ti_mode .eq. 2) then
      !Determine atom # that we are dealing with 
      iatom = gbl_img_atm_map(img_i) ! get atm number from idx
      if(ti_lst(ti_mask_piece,iatom) .ne. 0) cycle
    end if
    charge = img_qterm(img_i)
#else
  do my_chg_idx = 1, natom
    if (ti_mode .eq. 1 .or. ti_mode .eq. 2) then
      !Determine atom # that we are dealing with 
      iatom = gbl_img_atm_map(my_chg_idx) ! get atm number from idx
      if(ti_lst(ti_mask_piece,iatom) .ne. 0) cycle
    end if
    charge = img_qterm(my_chg_idx)
#endif

    i_base = ifracts(1, my_chg_idx)
    j_base = ifracts(2, my_chg_idx)
    k_base = ifracts(3, my_chg_idx)

    f1 = 0.d0
    f2 = 0.d0
    f3 = 0.d0

    ! We special-case order 4, the default.

    if (order .eq. 4) then
      do ith3 = 1, 4

        k = k_tbl(k_base + ith3)

#ifdef MPI
        if (k .ne. 0) then

#ifdef RANGE_DBG
          dbg_img_map(my_chgs(my_chg_idx)) = dbg_img_map(my_chgs(my_chg_idx))+1
#endif /* RANGE_DBG */

#endif
          do ith2 = 1, 4

            j = j_tbl(j_base + ith2)

            f1_term = theta(ith2, 2, my_chg_idx) * theta(ith3, 3, my_chg_idx)
            f2_term = dtheta(ith2, 2, my_chg_idx) * theta(ith3, 3, my_chg_idx)
            f3_term = theta(ith2, 2, my_chg_idx) * dtheta(ith3, 3, my_chg_idx)

! Force is negative of grad:

            qterm = q(i_tbl(i_base+1), j, k)
            f1 = f1 - qterm * dtheta(1, 1, my_chg_idx) * f1_term
            f2 = f2 - qterm * theta(1, 1, my_chg_idx) * f2_term
            f3 = f3 - qterm * theta(1, 1, my_chg_idx) * f3_term
            qterm = q(i_tbl(i_base+2), j, k)
            f1 = f1 - qterm * dtheta(2, 1, my_chg_idx) * f1_term
            f2 = f2 - qterm * theta(2, 1, my_chg_idx) * f2_term
            f3 = f3 - qterm * theta(2, 1, my_chg_idx) * f3_term
            qterm = q(i_tbl(i_base+3), j, k)
            f1 = f1 - qterm * dtheta(3, 1, my_chg_idx) * f1_term
            f2 = f2 - qterm * theta(3, 1, my_chg_idx) * f2_term
            f3 = f3 - qterm * theta(3, 1, my_chg_idx) * f3_term
            qterm = q(i_tbl(i_base+4), j, k)
            f1 = f1 - qterm * dtheta(4, 1, my_chg_idx) * f1_term
            f2 = f2 - qterm * theta(4, 1, my_chg_idx) * f2_term
            f3 = f3 - qterm * theta(4, 1, my_chg_idx) * f3_term

          end do

#ifdef MPI
        end if
#endif
      end do
    else if (order .eq. 6) then
      do ith3 = 1, 6

        k = k_tbl(k_base + ith3)

#ifdef MPI
        if (k .ne. 0) then

#ifdef RANGE_DBG
          dbg_img_map(my_chgs(my_chg_idx)) = dbg_img_map(my_chgs(my_chg_idx))+1
#endif /* RANGE_DBG */

#endif
          do ith2 = 1, 6

            j = j_tbl(j_base + ith2)

            f1_term = theta(ith2, 2, my_chg_idx) * theta(ith3, 3, my_chg_idx)
            f2_term = dtheta(ith2, 2, my_chg_idx) * theta(ith3, 3, my_chg_idx)
            f3_term = theta(ith2, 2, my_chg_idx) * dtheta(ith3, 3, my_chg_idx)

! Force is negative of grad:

            qterm = q(i_tbl(i_base+1), j, k)
            f1 = f1 - qterm * dtheta(1, 1, my_chg_idx) * f1_term
            f2 = f2 - qterm * theta(1, 1, my_chg_idx) * f2_term
            f3 = f3 - qterm * theta(1, 1, my_chg_idx) * f3_term
            qterm = q(i_tbl(i_base+2), j, k)
            f1 = f1 - qterm * dtheta(2, 1, my_chg_idx) * f1_term
            f2 = f2 - qterm * theta(2, 1, my_chg_idx) * f2_term
            f3 = f3 - qterm * theta(2, 1, my_chg_idx) * f3_term
            qterm = q(i_tbl(i_base+3), j, k)
            f1 = f1 - qterm * dtheta(3, 1, my_chg_idx) * f1_term
            f2 = f2 - qterm * theta(3, 1, my_chg_idx) * f2_term
            f3 = f3 - qterm * theta(3, 1, my_chg_idx) * f3_term
            qterm = q(i_tbl(i_base+4), j, k)
            f1 = f1 - qterm * dtheta(4, 1, my_chg_idx) * f1_term
            f2 = f2 - qterm * theta(4, 1, my_chg_idx) * f2_term
            f3 = f3 - qterm * theta(4, 1, my_chg_idx) * f3_term
            qterm = q(i_tbl(i_base+5), j, k)
            f1 = f1 - qterm * dtheta(5, 1, my_chg_idx) * f1_term
            f2 = f2 - qterm * theta(5, 1, my_chg_idx) * f2_term
            f3 = f3 - qterm * theta(5, 1, my_chg_idx) * f3_term
            qterm = q(i_tbl(i_base+6), j, k)
            f1 = f1 - qterm * dtheta(6, 1, my_chg_idx) * f1_term
            f2 = f2 - qterm * theta(6, 1, my_chg_idx) * f2_term
            f3 = f3 - qterm * theta(6, 1, my_chg_idx) * f3_term

          end do

#ifdef MPI
        end if
#endif
      end do
    else
      do ith3 = 1, order

        k = k_tbl(k_base + ith3)

#ifdef MPI
        if (k .ne. 0) then

#ifdef RANGE_DBG
          dbg_img_map(my_chgs(my_chg_idx)) = dbg_img_map(my_chgs(my_chg_idx))+1
#endif /* RANGE_DBG */

#endif
          do ith2 = 1, order

            j = j_tbl(j_base + ith2)

            f1_term = theta(ith2, 2, my_chg_idx) * theta(ith3, 3, my_chg_idx)
            f2_term = dtheta(ith2, 2, my_chg_idx) * theta(ith3, 3, my_chg_idx)
            f3_term = theta(ith2, 2, my_chg_idx) * dtheta(ith3, 3, my_chg_idx)

            do ith1 = 1, order

              qterm = q(i_tbl(i_base+ith1), j, k)

! Force is negative of grad:

              f1 = f1 - qterm * dtheta(ith1, 1, my_chg_idx) * f1_term
              f2 = f2 - qterm * theta(ith1, 1, my_chg_idx) * f2_term
              f3 = f3 - qterm * theta(ith1, 1, my_chg_idx) * f3_term
  
            end do

          end do

#ifdef MPI
        end if
#endif
      end do
    end if

    f1 = f1 * dnfft1 * charge
    f2 = f2 * dnfft2 * charge
    f3 = f3 * dnfft3 * charge

    if (is_orthog .ne. 0) then
      dfx = recip_11 * f1
      dfy = recip_22 * f2
      dfz = recip_33 * f3
    else
      dfx = recip(1, 1) * f1 + recip(1, 2) * f2 + recip(1, 3) * f3
      dfy = recip(2, 1) * f1 + recip(2, 2) * f2 + recip(2, 3) * f3
      dfz = recip(3, 1) * f1 + recip(3, 2) * f2 + recip(3, 3) * f3
    end if

#ifdef MPI
    img_frc(1, img_i) = img_frc(1, img_i) + dfx
    img_frc(2, img_i) = img_frc(2, img_i) + dfy
    img_frc(3, img_i) = img_frc(3, img_i) + dfz
#else
    img_frc(1, my_chg_idx) = img_frc(1, my_chg_idx) + dfx
    img_frc(2, my_chg_idx) = img_frc(2, my_chg_idx) + dfy
    img_frc(3, my_chg_idx) = img_frc(3, my_chg_idx) + dfz
#endif

  end do

#ifdef RANGE_DBG
  call mpi_reduce(dbg_img_map, dbg_img_map_res, size(dbg_img_map), &
                  mpi_integer, mpi_sum, 0, pmemd_comm, err_code_mpi)
  if (master) then
    do i = 1, natom
      if (dbg_img_map_res(i) .ne. order) then
        write(0,*)'DBG: IMAGE ', i, 'processed ', &
                  dbg_img_map_res(i), 'times in grad_sum code!!!'
      end if
    end do
  end if
#endif /* RANGE_DBG */

  return

end subroutine grad_sum

#define POT_ENES 1
#define VIRIAL 1
#include "scalar_sum_slab.i"
#undef VIRIAL
#define VIRIAL 0
#include "scalar_sum_slab.i"
#undef POT_ENES
#define POT_ENES 0
#include "scalar_sum_slab.i"
#undef POT_ENES
#undef VIRIAL

!*******************************************************************************
!
! Subroutine:  get_grid_weights
!
! Description: <TBS>
!              
!*******************************************************************************

#ifdef MPI
subroutine get_grid_weights(img_cnt, img_crd, ifracts, theta, dtheta, order, &
                            used_recip_img_cnt, used_recip_img_lst, &
                            my_chgs, my_chg_cnt)
#else
subroutine get_grid_weights(img_cnt, img_crd, ifracts, theta, dtheta, order)
#endif

  use bspline_mod
  use pme_slab_fft_mod
  use gbl_datatypes_mod
  use img_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: img_cnt
  double precision      :: img_crd(3, *)
  integer               :: ifracts(3, *)
  integer               :: order
  double precision      :: theta(order, 3, *)
  double precision      :: dtheta(order, 3, *)
#ifdef MPI
  integer               :: used_recip_img_cnt
  integer               :: used_recip_img_lst(*)
  integer               :: my_chgs(*), my_chg_cnt
#endif

! Local variables:

  double precision      :: box_x, box_y, box_z
  double precision      :: crd_x, crd_y, crd_z
  double precision      :: factor1, factor2, factor3
  double precision      :: fract(3)
  double precision      :: weight(3)
  integer               :: i, k

#ifdef MPI
  integer               :: kbot0, ktop1
  integer               :: lst_idx
  logical               :: kwraps
#endif /* MPI */

#ifdef MPI
  my_chg_cnt = 0

  kbot0 = xy_slab_start(mytaskid)
  ktop1 = kbot0 + my_xy_slab_cnt + order - 2
  kwraps = (ktop1 .ge. nfft3)
#endif /* MPI */

  box_x = ucell(1, 1)
  box_y = ucell(2, 2)
  box_z = ucell(3, 3)

  ! Scaling factors to get from cit table indexes to nfft indexes:

  factor1 = dble(nfft1) * recip(1, 1)
  factor2 = dble(nfft2) * recip(2, 2)
  factor3 = dble(nfft3) * recip(3, 3)

#ifdef MPI
  do lst_idx = 1, used_recip_img_cnt
     i = used_recip_img_lst(lst_idx)
#else /* MPI */
  do i = 1, img_cnt
#endif
    crd_x = img_crd(1, i)
    crd_y = img_crd(2, i)
    crd_z = img_crd(3, i)

    if (crd_x .ge. box_x) then
      crd_x = crd_x - box_x
    else if (crd_x .lt. 0.d0) then
      crd_x = crd_x + box_x
    end if

    if (crd_y .ge. box_y) then
      crd_y = crd_y - box_y
    else if (crd_y .lt. 0.d0) then
      crd_y = crd_y + box_y
    end if

    if (crd_z .ge. box_z) then
      crd_z = crd_z - box_z
    else if (crd_z .lt. 0.d0) then
      crd_z = crd_z + box_z
    end if

    fract(1) = factor1 * crd_x
    fract(2) = factor2 * crd_y
    fract(3) = factor3 * crd_z
#ifdef MPI
    k = int(fract(3))
#endif

#ifdef MPI

    if (kwraps) then
      if (k .lt. kbot0 .and. k .gt. ktop1 - nfft3) cycle
    else
      if (k .lt. kbot0 .or. k .gt. ktop1) cycle
    end if

#endif /* MPI */

#ifdef MPI
    my_chg_cnt = my_chg_cnt + 1
    my_chgs(my_chg_cnt) = i
#endif /* MPI */

    weight(:) = fract(:) - int(fract(:))

#ifdef MPI
    ifracts(:, my_chg_cnt) = int(fract(:)) - order

    call fill_bspline_1_3d(weight, order, &
                           theta(1, 1, my_chg_cnt), dtheta(1, 1, my_chg_cnt))
#else
    ifracts(:, i) = int(fract(:)) - order

    call fill_bspline_1_3d(weight, order, theta(1, 1, i), dtheta(1, 1, i))
#endif

  end do

  return

end subroutine get_grid_weights

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  claim_slab_recip_imgs
!
! Description: <TBS>
!              
!*******************************************************************************
subroutine claim_slab_recip_imgs(img_cnt, fraction, box, crd_idx_tbl,  &
#ifdef pmepack
                                 crd_q_pack, img_iac, mapped_img_cnt, &
#else
                                 img_crd, img_qterm, img_iac, mapped_img_cnt, &
#endif
                                 mapped_img_lst, img_atm_map, used_img_map, &
                                 used_img_cnt, used_img_lst)

  use cit_mod
  use pme_slab_fft_mod
  use gbl_datatypes_mod
  use img_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod

!precision setting for SPDP, DPDP, SPSP
#include "include_precision.i"
  implicit none

! Formal arguments:

  integer               :: img_cnt
  double precision      :: fraction(3, img_cnt)
  double precision      :: box(3)
  type(cit_tbl_rec)     :: crd_idx_tbl(0 : cit_tbl_x_dim - 1, &
                                       0 : cit_tbl_y_dim - 1, &
                                       0 : cit_tbl_z_dim - 1)
  double precision      :: img_crd(3, img_cnt)
  double precision      :: img_qterm(img_cnt)
  integer               :: img_iac(img_cnt)
  integer               :: mapped_img_cnt
  integer               :: mapped_img_lst(img_cnt)
  integer               :: img_atm_map(img_cnt)
  integer(byte)         :: used_img_map(img_cnt)
  integer               :: used_img_cnt
  integer               :: used_img_lst(*)

! Local variables:

  integer               :: k_bot0, k_top1
  integer               :: i, j, k
  integer               :: atm_id, img_id
  integer               :: z_bkt_lo, z_bkt_hi
  logical               :: z_wraps
  double precision      :: x_box, y_box, z_box
  double precision      :: z_lo, z_hi, z_cur

  gbl_used_recip_img_cnt = 0

  if (.not. i_do_recip) return    ! Map null range; not actually ever used...

  x_box = box(1)
  y_box = box(2)
  z_box = box(3)

  ! If there really is only one task doing reciprocal force calcs, we have to
  ! catch it here, because the code further below assumes that one process does
  ! not need to map the entire range of images.  Besides, this is more
  ! efficient for the special case...

  if (my_xy_slab_cnt .eq. fft_z_dim) then
    do img_id = 1, img_cnt
      if (used_img_map(img_id) .eq. 0) then
        used_img_map(img_id) = 1
        used_img_cnt = used_img_cnt + 1
        used_img_lst(used_img_cnt) = img_id
        if (img_iac(img_id) .eq. 0) then
          atm_id = img_atm_map(img_id)
          img_crd(1, img_id) = fraction(1, atm_id) * x_box
          img_crd(2, img_id) = fraction(2, atm_id) * y_box
          img_crd(3, img_id) = fraction(3, atm_id) * z_box
          img_qterm(img_id) = atm_qterm(atm_id)
          img_iac(img_id) = atm_iac(atm_id)
          mapped_img_cnt = mapped_img_cnt + 1
          mapped_img_lst(mapped_img_cnt) = img_id
        end if
      end if
      gbl_used_recip_img_cnt = gbl_used_recip_img_cnt + 1
      gbl_used_recip_img_lst(gbl_used_recip_img_cnt) = img_id
    end do
    return
  end if

  k_bot0 = xy_slab_start(mytaskid)
  k_top1 = k_bot0 + my_xy_slab_cnt + bspl_order - 2

  ! We have to allow for motion of skinnb, plus we have to allow for fact that
  ! atom z crds are rounded down to the grid point.  Hence, the + 1 for z_hi

  z_lo = dble(k_bot0) * z_box / dble(nfft3) - (0.5d0 * skinnb)
  z_hi = dble(k_top1 + 1) * z_box / dble(nfft3) + (0.5d0 * skinnb)

  if (z_lo .lt. 0.d0) then
    z_wraps = .true.
    z_lo = z_lo + z_box
  else if (z_hi .ge. z_box) then
    z_wraps = .true.
    z_hi = z_hi - z_box
  else
    z_wraps = .false.
  end if

  z_bkt_lo = int(z_lo * dble(cit_tbl_z_dim) / z_box)
  z_bkt_hi = int(z_hi * dble(cit_tbl_z_dim) / z_box)

  if (z_bkt_lo .eq. z_bkt_hi) z_wraps = .false.

  if (z_wraps) then

    do k = z_bkt_lo, cit_tbl_z_dim - 1
      do j = 0, cit_tbl_y_dim - 1
        do i = 0, cit_tbl_x_dim - 1
          do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
            atm_id = img_atm_map(img_id)
            z_cur = fraction(3, atm_id) * z_box
            if (z_cur .ge. z_lo) then
              if (used_img_map(img_id) .eq. 0) then
                used_img_map(img_id) = 1
                used_img_cnt = used_img_cnt + 1
                used_img_lst(used_img_cnt) = img_id
                if (img_iac(img_id) .eq. 0) then
                  img_crd(1, img_id) = fraction(1, atm_id) * x_box
                  img_crd(2, img_id) = fraction(2, atm_id) * y_box
                  img_crd(3, img_id) = z_cur
                  img_qterm(img_id) = atm_qterm(atm_id)
                  img_iac(img_id) = atm_iac(atm_id)
                  mapped_img_cnt = mapped_img_cnt + 1
                  mapped_img_lst(mapped_img_cnt) = img_id
                end if
              end if
              gbl_used_recip_img_cnt = gbl_used_recip_img_cnt + 1
              gbl_used_recip_img_lst(gbl_used_recip_img_cnt) = img_id
            end if
          end do
        end do
      end do
    end do

    do k = 0, z_bkt_hi
      do j = 0, cit_tbl_y_dim - 1
        do i = 0, cit_tbl_x_dim - 1
          do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
            atm_id = img_atm_map(img_id)
            z_cur = fraction(3, atm_id) * z_box
            if (z_cur .lt. z_hi) then
              if (used_img_map(img_id) .eq. 0) then
                used_img_map(img_id) = 1
                used_img_cnt = used_img_cnt + 1
                used_img_lst(used_img_cnt) = img_id
                if (img_iac(img_id) .eq. 0) then
                  img_crd(1, img_id) = fraction(1, atm_id) * x_box
                  img_crd(2, img_id) = fraction(2, atm_id) * y_box
                  img_crd(3, img_id) = z_cur
                  img_qterm(img_id) = atm_qterm(atm_id)
                  img_iac(img_id) = atm_iac(atm_id)
                  mapped_img_cnt = mapped_img_cnt + 1
                  mapped_img_lst(mapped_img_cnt) = img_id
                end if
              end if
              gbl_used_recip_img_cnt = gbl_used_recip_img_cnt + 1
              gbl_used_recip_img_lst(gbl_used_recip_img_cnt) = img_id
            end if
          end do
        end do
      end do
    end do

  else

    do k = z_bkt_lo, z_bkt_hi
      do j = 0, cit_tbl_y_dim - 1
        do i = 0, cit_tbl_x_dim - 1
          do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
            atm_id = img_atm_map(img_id)
            z_cur = fraction(3, atm_id) * z_box
            if (z_cur .ge. z_lo .and. z_cur .lt. z_hi) then
              if (used_img_map(img_id) .eq. 0) then
                used_img_map(img_id) = 1
                used_img_cnt = used_img_cnt + 1
                used_img_lst(used_img_cnt) = img_id
                if (img_iac(img_id) .eq. 0) then
                  img_crd(1, img_id) = fraction(1, atm_id) * x_box
                  img_crd(2, img_id) = fraction(2, atm_id) * y_box
                  img_crd(3, img_id) = z_cur
                  img_qterm(img_id) = atm_qterm(atm_id)
                  img_iac(img_id) = atm_iac(atm_id)
                  mapped_img_cnt = mapped_img_cnt + 1
                  mapped_img_lst(mapped_img_cnt) = img_id
                end if
              end if
              gbl_used_recip_img_cnt = gbl_used_recip_img_cnt + 1
              gbl_used_recip_img_lst(gbl_used_recip_img_cnt) = img_id
            end if
          end do
        end do
      end do
    end do

  end if

  return

end subroutine claim_slab_recip_imgs
#endif /* MPI */

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  claim_slab_recip_imgs_nonorthog
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine claim_slab_recip_imgs_nonorthog(img_cnt, fraction, crd_idx_tbl, &
                                           img_crd, img_qterm, img_iac, &
                                           mapped_img_cnt, mapped_img_lst, &
                                           img_atm_map, used_img_map, &
                                           used_img_cnt, used_img_lst)

  use cit_mod
  use pme_slab_fft_mod
  use gbl_datatypes_mod
  use img_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use pbc_mod

!precision setting for SPDP, DPDP, SPSP
#include "include_precision.i"

  implicit none

! Formal arguments:

  integer               :: img_cnt
  double precision      :: fraction(3, img_cnt)
  type(cit_tbl_rec)     :: crd_idx_tbl(0 : cit_tbl_x_dim - 1, &
                                       0 : cit_tbl_y_dim - 1, &
                                       0 : cit_tbl_z_dim - 1)
  double precision      :: img_crd(3, img_cnt)
  double precision      :: img_qterm(img_cnt)
  integer               :: img_iac(img_cnt)
  integer               :: mapped_img_cnt
  integer               :: mapped_img_lst(img_cnt)
  integer               :: img_atm_map(img_cnt)
  integer(byte)         :: used_img_map(img_cnt)
  integer               :: used_img_cnt
  integer               :: used_img_lst(*)

! Local variables:

  integer               :: k_bot0, k_top1
  integer               :: i, j, k
  integer               :: atm_id, img_id
  integer               :: z_bkt_lo, z_bkt_hi
  logical               :: z_wraps
  double precision      :: ucell_stk(3, 3)
  double precision      :: f1, f2, f3
  double precision      :: z_lo, z_hi

  ucell_stk(:,:) = ucell(:,:)

  gbl_used_recip_img_cnt = 0

  if (.not. i_do_recip) return    ! Map null range; not actually ever used...

  ! If there really is only one task doing reciprocal force calcs, we have to
  ! catch it here, because the code further below assumes that one process does
  ! not need to map the entire range of images.  Besides, this is more
  ! efficient for the special case...

  if (my_xy_slab_cnt .eq. fft_z_dim) then
    do img_id = 1, img_cnt
      if (used_img_map(img_id) .eq. 0) then
        used_img_map(img_id) = 1
        used_img_cnt = used_img_cnt + 1
        used_img_lst(used_img_cnt) = img_id
        if (img_iac(img_id) .eq. 0) then
          atm_id = img_atm_map(img_id)
          f1 = fraction(1, atm_id)
          f2 = fraction(2, atm_id)
          f3 = fraction(3, atm_id)
          ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0, so
          ! we can simplify the expression in this critical inner loop
          img_crd(1, img_id) = f1 * ucell_stk(1, 1) + &
                               f2 * ucell_stk(1, 2) + &
                               f3 * ucell_stk(1, 3)
          img_crd(2, img_id) = f2 * ucell_stk(2, 2) + &
                               f3 * ucell_stk(2, 3)
          img_crd(3, img_id) = f3 * ucell_stk(3, 3)
          img_qterm(img_id) = atm_qterm(atm_id)
          img_iac(img_id) = atm_iac(atm_id)
          mapped_img_cnt = mapped_img_cnt + 1
          mapped_img_lst(mapped_img_cnt) = img_id
        end if
      end if
      gbl_used_recip_img_cnt = gbl_used_recip_img_cnt + 1
      gbl_used_recip_img_lst(gbl_used_recip_img_cnt) = img_id
    end do
    return
  end if

  k_bot0 = xy_slab_start(mytaskid)
  k_top1 = k_bot0 + my_xy_slab_cnt + bspl_order - 2

  ! We have to allow for motion of skinnb, plus we have to allow for fact that
  ! atom z crds are rounded down to the grid point.  Hence, the + 1 for z_hi

  z_lo = dble(k_bot0) / dble(nfft3) - &
         (0.5d0 * cut_factor(3) * skinnb) / ucell_stk(3,3)      ! a fractional

  z_hi = dble(k_top1 + 1) / dble(nfft3) + &
         (0.5d0 * cut_factor(3) * skinnb) / ucell_stk(3,3)      ! a fractional

  if (z_lo .lt. 0.d0) then
    z_wraps = .true.
    z_lo = z_lo + 1.d0
  else if (z_hi .ge. 1.d0) then
    z_wraps = .true.
    z_hi = z_hi - 1.d0
  else
    z_wraps = .false.
  end if

  z_bkt_lo = int(z_lo * dble(cit_tbl_z_dim))
  z_bkt_hi = int(z_hi * dble(cit_tbl_z_dim))

  if (z_bkt_lo .eq. z_bkt_hi) z_wraps = .false.

  if (z_wraps) then

    do k = z_bkt_lo, cit_tbl_z_dim - 1
      do j = 0, cit_tbl_y_dim - 1
        do i = 0, cit_tbl_x_dim - 1
          do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
            atm_id = img_atm_map(img_id)
            f3 = fraction(3, atm_id)
            if (f3 .ge. z_lo) then
              if (used_img_map(img_id) .eq. 0) then
                used_img_map(img_id) = 1
                used_img_cnt = used_img_cnt + 1
                used_img_lst(used_img_cnt) = img_id
                if (img_iac(img_id) .eq. 0) then
                  f1 = fraction(1, atm_id)
                  f2 = fraction(2, atm_id)
                  ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0, so
                  ! we can simplify the expression in this critical inner loop
                  img_crd(1, img_id) = f1 * ucell_stk(1, 1) + &
                                       f2 * ucell_stk(1, 2) + &
                                       f3 * ucell_stk(1, 3)
                  img_crd(2, img_id) = f2 * ucell_stk(2, 2) + &
                                       f3 * ucell_stk(2, 3)
                  img_crd(3, img_id) = f3 * ucell_stk(3, 3)
                  img_qterm(img_id) = atm_qterm(atm_id)
                  img_iac(img_id) = atm_iac(atm_id)
                  mapped_img_cnt = mapped_img_cnt + 1
                  mapped_img_lst(mapped_img_cnt) = img_id
                end if
              end if
              gbl_used_recip_img_cnt = gbl_used_recip_img_cnt + 1
              gbl_used_recip_img_lst(gbl_used_recip_img_cnt) = img_id
            end if
          end do
        end do
      end do
    end do

    do k = 0, z_bkt_hi
      do j = 0, cit_tbl_y_dim - 1
        do i = 0, cit_tbl_x_dim - 1
          do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
            atm_id = img_atm_map(img_id)
            f3 = fraction(3, atm_id)
            if (f3 .lt. z_hi) then
              if (used_img_map(img_id) .eq. 0) then
                used_img_map(img_id) = 1
                used_img_cnt = used_img_cnt + 1
                used_img_lst(used_img_cnt) = img_id
                if (img_iac(img_id) .eq. 0) then
                  f1 = fraction(1, atm_id)
                  f2 = fraction(2, atm_id)
                  ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0, so
                  ! we can simplify the expression in this critical inner loop
                  img_crd(1, img_id) = f1 * ucell_stk(1, 1) + &
                                       f2 * ucell_stk(1, 2) + &
                                       f3 * ucell_stk(1, 3)
                  img_crd(2, img_id) = f2 * ucell_stk(2, 2) + &
                                       f3 * ucell_stk(2, 3)
                  img_crd(3, img_id) = f3 * ucell_stk(3, 3)
                  img_qterm(img_id) = atm_qterm(atm_id)
                  img_iac(img_id) = atm_iac(atm_id)
                  mapped_img_cnt = mapped_img_cnt + 1
                  mapped_img_lst(mapped_img_cnt) = img_id
                end if
              end if
              gbl_used_recip_img_cnt = gbl_used_recip_img_cnt + 1
              gbl_used_recip_img_lst(gbl_used_recip_img_cnt) = img_id
            end if
          end do
        end do
      end do
    end do

  else

    do k = z_bkt_lo, z_bkt_hi
      do j = 0, cit_tbl_y_dim - 1
        do i = 0, cit_tbl_x_dim - 1
          do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
            atm_id = img_atm_map(img_id)
            f3 = fraction(3, atm_id)
            if (f3 .ge. z_lo .and. f3 .lt. z_hi) then
              if (used_img_map(img_id) .eq. 0) then
                used_img_map(img_id) = 1
                used_img_cnt = used_img_cnt + 1
                used_img_lst(used_img_cnt) = img_id
                if (img_iac(img_id) .eq. 0) then
                  f1 = fraction(1, atm_id)
                  f2 = fraction(2, atm_id)
                  ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0, so
                  ! we can simplify the expression in this critical inner loop
                  img_crd(1, img_id) = f1 * ucell_stk(1, 1) + &
                                       f2 * ucell_stk(1, 2) + &
                                       f3 * ucell_stk(1, 3)
                  img_crd(2, img_id) = f2 * ucell_stk(2, 2) + &
                                       f3 * ucell_stk(2, 3)
                  img_crd(3, img_id) = f3 * ucell_stk(3, 3)
                  img_qterm(img_id) = atm_qterm(atm_id)
                  img_iac(img_id) = atm_iac(atm_id)
                  mapped_img_cnt = mapped_img_cnt + 1
                  mapped_img_lst(mapped_img_cnt) = img_id
                end if
              end if
              gbl_used_recip_img_cnt = gbl_used_recip_img_cnt + 1
              gbl_used_recip_img_lst(gbl_used_recip_img_cnt) = img_id
            end if
          end do
        end do
      end do
    end do

  end if

  return

end subroutine claim_slab_recip_imgs_nonorthog
#endif /* MPI */

!*******************************************************************************
!
! Subroutine:  get_grid_weights_nonorthog
!
! Description: <TBS>
!              
!*******************************************************************************

#ifdef MPI
subroutine get_grid_weights_nonorthog(img_cnt, crd, img_atm_map, &
                                      ifracts, theta, dtheta, order, &
                                      used_recip_img_cnt, used_recip_img_lst, &
                                      my_chgs, my_chg_cnt)
#else
subroutine get_grid_weights_nonorthog(img_cnt, crd, img_atm_map, ifracts, &
                                      theta, dtheta, order)
#endif

  use bspline_mod
  use gbl_datatypes_mod
  use pme_slab_fft_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: img_cnt
  double precision      :: crd(3, *)
  integer               :: img_atm_map(*)
  integer               :: ifracts(3, *)
  integer               :: order
  double precision      :: theta(order, 3, *)
  double precision      :: dtheta(order, 3, *)
#ifdef MPI
  integer               :: used_recip_img_cnt
  integer               :: used_recip_img_lst(*)
  integer               :: my_chgs(*), my_chg_cnt
#endif

! Local variables:

  double precision      :: crd_x, crd_y, crd_z
  double precision      :: fract(3)
  double precision      :: weight(3)
  double precision      :: recip_stk(3, 3)
  integer               :: atm_id
  integer               :: i, k

#ifdef MPI
  integer               :: kbot0, ktop1
  logical               :: kwraps
  integer               :: lst_idx
#endif /* MPI */

#ifdef MPI
  my_chg_cnt = 0

  kbot0 = xy_slab_start(mytaskid)
  ktop1 = kbot0 + my_xy_slab_cnt + order - 2
  kwraps = (ktop1 .ge. nfft3)
#endif /* MPI */

  recip_stk(:,:) = recip(:,:)


#ifdef MPI
  do lst_idx = 1, used_recip_img_cnt
    i = used_recip_img_lst(lst_idx)
#else
  do i = 1, img_cnt
#endif

! Unfortunately we need fractional coords that are not available without going
! back to atm_crd().  This is the case because the coordinates in crd
! were based on fractional coords at list build time, but as movement occurs
! you can't reconstruct the fractional.

    atm_id = img_atm_map(i)
    crd_x = crd(1, atm_id)
    crd_y = crd(2, atm_id)
    crd_z = crd(3, atm_id)

    fract(1) = crd_x * recip_stk(1,1) + &
               crd_y * recip_stk(2,1) +  &
               crd_z * recip_stk(3,1)

    fract(1) = dble(nfft1) * (fract(1) - dnint(fract(1)) + 0.5d0)

    fract(2) = crd_x * recip_stk(1,2) + &
               crd_y * recip_stk(2,2) + &
               crd_z * recip_stk(3,2)

    fract(2) = dble(nfft2) * (fract(2) - dnint(fract(2)) + 0.5d0)

    fract(3) = crd_x * recip_stk(1,3) + &
               crd_y * recip_stk(2,3) + &
               crd_z * recip_stk(3,3)

    fract(3) = dble(nfft3) * (fract(3) - dnint(fract(3)) + 0.5d0)

#ifdef MPI
    k = int(fract(3))

    if (kwraps) then
      if (k .lt. kbot0 .and. k .gt. ktop1 - nfft3) cycle
    else
      if (k .lt. kbot0 .or. k .gt. ktop1) cycle
    end if

    my_chg_cnt = my_chg_cnt + 1
    my_chgs(my_chg_cnt) = i
#endif

    weight(:) = fract(:) - int(fract(:))

#ifdef MPI
    ifracts(:, my_chg_cnt) = int(fract(:)) - order

    call fill_bspline_1_3d(weight, order, &
                           theta(1, 1, my_chg_cnt), dtheta(1, 1, my_chg_cnt))
#else
    ifracts(:, i) = int(fract(:)) - order

    call fill_bspline_1_3d(weight, order, theta(1, 1, i), dtheta(1, 1, i))
#endif

  end do

  return

end subroutine get_grid_weights_nonorthog

end module pme_slab_recip_mod
