#include "copyright.i"

!*******************************************************************************
!
! Module: pme_blk_recip_mod
!
! Description: <TBS>
!
! NOTE NOTE NOTE: This code assumes frc_int .eq. 0 and should only be used
!                 under these conditions!!!
!*******************************************************************************

module pme_blk_recip_mod

  use pme_recip_dat_mod
!precision setting for SPDP, DPDP, SPSP
#include "include_precision.i"

  implicit none

  private       ! Everything here is private unless specified public

  ! All the code is mpi-specific...

#ifdef MPI
  public        :: pme_blk_recip_setup, &
                   claim_blk_recip_imgs, &
                   claim_blk_recip_imgs_nonorthog, &
                   do_blk_pmesh_kspace

contains

!*******************************************************************************
!
! Subroutine:  pme_blk_recip_setup
!
! Description:
!              
!*******************************************************************************

subroutine pme_blk_recip_setup(frc_ene_numtasks, frc_ene_task_lst, &
                               num_ints, num_reals)

  use loadbal_mod
  use pmemd_lib_mod
  use pme_blk_fft_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: frc_ene_numtasks
  integer, intent(in)           :: frc_ene_task_lst(frc_ene_numtasks)

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: min_numtasks, max_numtasks

  call fft_dat_setup(num_ints, num_reals)

  max_numtasks = frc_ene_numtasks
  min_numtasks = int(blk_fft_workload_estimate * dble(max_numtasks) + 0.5d0)

  ! We make no special effort on handling respa or minimization yet; will
  ! have to check out later...
  call blk_fft_setup(min_numtasks, recip_numtasks, &
                     frc_ene_numtasks, frc_ene_task_lst, &
                     0, num_ints, num_reals)

  if (recip_numtasks .lt. frc_ene_numtasks) then
    fft_redist_enabled = .true.
    fft_redist_needed = .true.
  end if

  return

end subroutine pme_blk_recip_setup

!*******************************************************************************
!
! Subroutine:  do_blk_pmesh_kspace
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

subroutine do_blk_pmesh_kspace(img_frc, eer, virial, need_pot_enes, &
                               need_virials)

  use img_mod
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use pbc_mod
  use pme_blk_fft_mod
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

  integer               :: my_chg_cnt

  double precision      :: theta(bspl_order * 3 * gbl_used_recip_img_cnt)
  double precision      :: dtheta(bspl_order * 3 * gbl_used_recip_img_cnt)
  integer               :: my_chgs(gbl_used_recip_img_cnt)
  double precision      :: xyz_q(2 * fft_x_dim, &
                                 fft_y_cnts1(my_grid_idx1), &
                                 fft_z_cnts(my_grid_idx2))
  double precision      :: zxy_q(2 * fft_z_dim, &
                                 fft_x_cnts(my_grid_idx1), &
                                 fft_y_cnts2(my_grid_idx2))

  integer               :: ifracts(3, gbl_used_recip_img_cnt)

  ! Exponential tables are only practical to use for orthogonal unit cells, and
  ! only need updating when constant pressure runs are being done.

  if (ntp .gt. 0 .and. is_orthog .ne. 0) call load_m_exp_tbls

  if (is_orthog .ne. 0) then
    call get_grid_weights(natom, gbl_img_crd, ifracts, theta, dtheta, &
                          bspl_order, &
                          gbl_used_recip_img_cnt, gbl_used_recip_img_lst, &
                          my_chgs, my_chg_cnt)
  else
    call get_grid_weights_nonorthog(natom, atm_crd, gbl_img_atm_map, &
                                    ifracts, theta, dtheta, bspl_order, &
                                    gbl_used_recip_img_cnt, &
                                    gbl_used_recip_img_lst, &
                                    my_chgs, my_chg_cnt)
  end if

  call update_pme_time(bspline_timer)

! Fill Charge Grid.  Charges are approximated on an even grid.

  call fill_charge_grid(xyz_q, theta, bspl_order, gbl_img_qterm, ifracts, &
                        my_chgs, my_chg_cnt)

  call update_pme_time(grid_charges_timer)

  call blk_fft3drc_forward(xyz_q, zxy_q, fft_x_dim, fft_y_dim, fft_z_dim)

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

  call blk_fft3drc_back(zxy_q, xyz_q, fft_x_dim, fft_y_dim, fft_z_dim)

  call update_pme_time(fft_timer)

  call grad_sum(xyz_q, theta, dtheta, bspl_order, gbl_img_qterm, &
                img_frc, ifracts, my_chgs, my_chg_cnt)

  call update_pme_time(grad_sum_timer)

  return

end subroutine do_blk_pmesh_kspace

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

subroutine fill_charge_grid(q, theta, order, img_qterm, ifracts, &
                            my_chgs, my_chg_cnt)

  use img_mod
  use mdin_ewald_dat_mod
  use pme_blk_fft_mod
  use parallel_dat_mod
  use prmtop_dat_mod
! BEGIN DBG
  use pmemd_lib_mod
! END DBG
  use ti_mod

  implicit none

! Formal arguments:

  integer               :: order
  double precision      :: theta(order, 3, *)
  double precision      :: img_qterm(*)
  integer               :: ifracts(3, *)
  double precision      :: q(2 * fft_x_dim, &
                             fft_y_cnts1(my_grid_idx1), &
                             fft_z_cnts(my_grid_idx2))
  integer               :: my_chgs(*), my_chg_cnt

! Local variables:

  integer               :: i, j, k
  integer               :: i_base, j_base, k_base
  integer               :: ith1, ith2, ith3
  integer               :: my_chg_idx
  double precision      :: charge
  double precision      :: k_term, j_term
! integer               :: kbot0
  integer               :: y_off, z_off
  integer               :: y_cnt, z_cnt
  integer               :: jbot, jtop
  integer               :: kbot, ktop
  integer               :: i_tbl(-order : nfft1)
  integer               :: j_tbl(-order : nfft2)
  integer               :: k_tbl(-order : nfft3)
  ! For TI
  integer               :: iatom
#ifdef RANGE_DBG
  integer               :: dbg_img_map(natom)
  integer               :: dbg_img_map_res(natom)

  dbg_img_map(:) = 0
#endif /* RANGE_DBG */

! Old usages; kept for documentation during code conversion...
! kbot0 = my_xy_slab_start
! kbot = kbot0 + 1
! ktop = kbot0 + my_xy_slab_cnt

  y_off = fft_y_offs1(my_grid_idx1)
  y_cnt = fft_y_cnts1(my_grid_idx1)

  z_off = fft_z_offs(my_grid_idx2)
  z_cnt = fft_z_cnts(my_grid_idx2)

! Zero the Charge grids:

  call zero_q_array(q, 2 * fft_x_dim * y_cnt * z_cnt)

  ! Initialize the indexing tables.  It actually produces faster code to
  ! do this here rather than caching the info remotely.

  do i = 0, nfft1
    i_tbl(i) = i + 1
  end do
  do i = -order, -1
    i_tbl(i) = i + nfft1 + 1
  end do

  jbot = 1
  jtop = y_cnt
  do j = 0, nfft2
    j_tbl(j) = j - y_off + 1
    if (j_tbl(j) .lt. jbot .or. j_tbl(j) .gt. jtop) j_tbl(j) = 0
  end do
  do j = -order, -1
    j_tbl(j) = j + nfft2 - y_off + 1
    if (j_tbl(j) .lt. jbot .or. j_tbl(j) .gt. jtop) j_tbl(j) = 0
  end do

  kbot = 1
  ktop = z_cnt
  do k = 0, nfft3
    k_tbl(k) = k - z_off + 1
    if (k_tbl(k) .lt. kbot .or. k_tbl(k) .gt. ktop) k_tbl(k) = 0
  end do
  do k = -order, -1
    k_tbl(k) = k + nfft3 - z_off + 1
    if (k_tbl(k) .lt. kbot .or. k_tbl(k) .gt. ktop) k_tbl(k) = 0
  end do

  ! We special-case order 4, the default.

  if (order .eq. 4) then
    do my_chg_idx = 1, my_chg_cnt

      if (ti_mode .eq. 1 .or. ti_mode .eq. 2) then
        !Determine atom # that we are dealing with 
        iatom = gbl_img_atm_map(my_chgs(my_chg_idx)) ! get atm number from idx
        if(ti_lst(ti_mask_piece,iatom) .ne. 0) cycle
      end if

      charge = img_qterm(my_chgs(my_chg_idx))

      i_base = ifracts(1, my_chg_idx)
      j_base = ifracts(2, my_chg_idx)
      k_base = ifracts(3, my_chg_idx)

      do ith3 = 1, 4

        k = k_tbl(k_base + ith3)

        if (k .ne. 0) then

#ifdef RANGE_DBG
          dbg_img_map(my_chgs(my_chg_idx)) = dbg_img_map(my_chgs(my_chg_idx))+1
#endif /* RANGE_DBG */
          k_term =  theta(ith3, 3, my_chg_idx) * charge
          do ith2 = 1, 4
            j = j_tbl(j_base + ith2)
            if (j .ne. 0) then
              j_term = k_term * theta(ith2, 2, my_chg_idx)
              q(i_tbl(i_base+1), j, k) = q(i_tbl(i_base+1), j, k) + &
                                         theta(1, 1, my_chg_idx) * j_term
              q(i_tbl(i_base+2), j, k) = q(i_tbl(i_base+2), j, k) + &
                                         theta(2, 1, my_chg_idx) * j_term
              q(i_tbl(i_base+3), j, k) = q(i_tbl(i_base+3), j, k) + &
                                         theta(3, 1, my_chg_idx) * j_term
              q(i_tbl(i_base+4), j, k) = q(i_tbl(i_base+4), j, k) + &
                                         theta(4, 1, my_chg_idx) * j_term
            end if
          end do
        end if
      end do
    end do
  else
    do my_chg_idx = 1, my_chg_cnt

      if (ti_mode .eq. 1 .or. ti_mode .eq. 2) then
        !Determine atom # that we are dealing with 
        iatom = gbl_img_atm_map(my_chgs(my_chg_idx)) ! get atm number from idx
        if(ti_lst(ti_mask_piece,iatom) .ne. 0) cycle
      end if

      charge = img_qterm(my_chgs(my_chg_idx))

      i_base = ifracts(1, my_chg_idx)
      j_base = ifracts(2, my_chg_idx)
      k_base = ifracts(3, my_chg_idx)

      do ith3 = 1, order

        k = k_tbl(k_base + ith3)

        if (k .ne. 0) then

#ifdef RANGE_DBG
          dbg_img_map(my_chgs(my_chg_idx)) = dbg_img_map(my_chgs(my_chg_idx))+1
#endif /* RANGE_DBG */
          k_term =  theta(ith3, 3, my_chg_idx) * charge
          do ith2 = 1, order
            j = j_tbl(j_base + ith2)
            if (j .ne. 0) then
              j_term = k_term * theta(ith2, 2, my_chg_idx)
              do ith1 = 1, order
                i = i_tbl(i_base + ith1)
                q(i, j, k) = q(i, j, k) + theta(ith1, 1, my_chg_idx) * j_term
              end do
            end if
          end do
        end if
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

subroutine grad_sum(q, theta, dtheta, order, img_qterm, img_frc, ifracts, &
                    my_chgs, my_chg_cnt)

  use img_mod
  use mdin_ewald_dat_mod
  use prmtop_dat_mod
  use pbc_mod
  use pme_blk_fft_mod
  use ti_mod

  implicit none

! Formal arguments:

  double precision      :: q(2 * fft_x_dim, &
                             fft_y_cnts1(my_grid_idx1), &
                             fft_z_cnts(my_grid_idx2))
  integer               :: order
  double precision      :: theta(order, 3, *)
  double precision      :: dtheta(order, 3, *)
  double precision      :: img_qterm(*)
  double precision      :: img_frc(3, *)
  integer               :: ifracts(3, *)

  integer               :: my_chgs(*), my_chg_cnt

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
  integer               :: img_i
  integer               :: y_cnt, y_off
  integer               :: z_cnt, z_off
  integer               :: jbot, jtop
  integer               :: kbot, ktop
  integer               :: i_tbl(-order : nfft1)
  integer               :: j_tbl(-order : nfft2)
  integer               :: k_tbl(-order : nfft3)
  ! For TI:
  integer               :: iatom
#ifdef RANGE_DBG
  integer               :: dbg_img_map(natom)
  integer               :: dbg_img_map_res(natom)

  dbg_img_map(:) = 0
#endif /* RANGE_DBG */

  y_cnt = fft_y_cnts1(my_grid_idx1)
  y_off = fft_y_offs1(my_grid_idx1)
  z_cnt = fft_z_cnts(my_grid_idx2)
  z_off = fft_z_offs(my_grid_idx2)

  ! Initialize the indexing tables.  It actually produces faster code to
  ! do this here rather than caching the info remotely.

  do i = 0, nfft1
    i_tbl(i) = i + 1
  end do
  do i = -order, -1
    i_tbl(i) = i + nfft1 + 1
  end do

  jbot = 1
  jtop = y_cnt
  do j = 0, nfft2
    j_tbl(j) = j - y_off + 1
    if (j_tbl(j) .lt. jbot .or. j_tbl(j) .gt. jtop) j_tbl(j) = 0
  end do
  do j = -order, -1
    j_tbl(j) = j + nfft2 - y_off + 1
    if (j_tbl(j) .lt. jbot .or. j_tbl(j) .gt. jtop) j_tbl(j) = 0
  end do

  kbot = 1
  ktop = z_cnt
  do k = 0, nfft3
    k_tbl(k) = k - z_off + 1
    if (k_tbl(k) .lt. kbot .or. k_tbl(k) .gt. ktop) k_tbl(k) = 0
  end do
  do k = -order, -1
    k_tbl(k) = k + nfft3 - z_off + 1
    if (k_tbl(k) .lt. kbot .or. k_tbl(k) .gt. ktop) k_tbl(k) = 0
  end do

  recip_11 = recip(1, 1)
  recip_22 = recip(2, 2)
  recip_33 = recip(3, 3)

  dnfft1 = dble(nfft1)
  dnfft2 = dble(nfft2)
  dnfft3 = dble(nfft3)

  do my_chg_idx = 1, my_chg_cnt
    img_i = my_chgs(my_chg_idx)

    if (ti_mode .eq. 1 .or. ti_mode .eq. 2) then
      !Determine atom # that we are dealing with 
      iatom = gbl_img_atm_map(img_i) ! get atom number from index
      if(ti_lst(ti_mask_piece,iatom) .ne. 0) cycle
    end if

    charge = img_qterm(img_i)

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

        if (k .ne. 0) then

#ifdef RANGE_DBG
          dbg_img_map(my_chgs(my_chg_idx)) = dbg_img_map(my_chgs(my_chg_idx))+1
#endif /* RANGE_DBG */

          do ith2 = 1, 4

            j = j_tbl(j_base + ith2)

            if (j .ne. 0) then
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
            end if

          end do

        end if
      end do
    else
      do ith3 = 1, order

        k = k_tbl(k_base + ith3)

        if (k .ne. 0) then

#ifdef RANGE_DBG
          dbg_img_map(my_chgs(my_chg_idx)) = dbg_img_map(my_chgs(my_chg_idx))+1
#endif /* RANGE_DBG */

          do ith2 = 1, order

            j = j_tbl(j_base + ith2)

            if (j .ne. 0) then
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
            end if

          end do

        end if
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

    img_frc(1, img_i) = img_frc(1, img_i) + dfx
    img_frc(2, img_i) = img_frc(2, img_i) + dfy
    img_frc(3, img_i) = img_frc(3, img_i) + dfz

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
#include "scalar_sum_blk.i"
#undef VIRIAL
#define VIRIAL 0
#include "scalar_sum_blk.i"
#undef POT_ENES
#define POT_ENES 0
#include "scalar_sum_blk.i"
#undef POT_ENES
#undef VIRIAL

!*******************************************************************************
!
! Subroutine:  get_grid_weights
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine get_grid_weights(img_cnt, img_crd, ifracts, theta, dtheta, order, &
                            used_recip_img_cnt, used_recip_img_lst, &
                            my_chgs, my_chg_cnt)

  use bspline_mod
  use pme_blk_fft_mod
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
  integer               :: used_recip_img_cnt
  integer               :: used_recip_img_lst(*)
  integer               :: my_chgs(*), my_chg_cnt

! Local variables:

  double precision      :: box_x, box_y, box_z
  double precision      :: crd_x, crd_y, crd_z
  double precision      :: factor1, factor2, factor3
  double precision      :: fract(3)
  double precision      :: weight(3)
  integer               :: i, k

  integer               :: kbot0, ktop1
  integer               :: lst_idx
  logical               :: kwraps
! BEGIN DBG
  integer               :: i_base, j_base, k_base
! END DBG

  my_chg_cnt = 0

  kbot0 = fft_z_offs(my_grid_idx2)
  ktop1 = kbot0 + fft_z_cnts(my_grid_idx2) + order - 2
  kwraps = (ktop1 .ge. nfft3)

  box_x = ucell(1, 1)
  box_y = ucell(2, 2)
  box_z = ucell(3, 3)

  ! Scaling factors to get from cit table indexes to nfft indexes:

  factor1 = dble(nfft1) * recip(1, 1)
  factor2 = dble(nfft2) * recip(2, 2)
  factor3 = dble(nfft3) * recip(3, 3)


  do lst_idx = 1, used_recip_img_cnt
    i = used_recip_img_lst(lst_idx)

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

    k = int(fract(3))

    if (kwraps) then
      if (k .lt. kbot0 .and. k .gt. ktop1 - nfft3) cycle
    else
      if (k .lt. kbot0 .or. k .gt. ktop1) cycle
    end if

    my_chg_cnt = my_chg_cnt + 1
    my_chgs(my_chg_cnt) = i

    weight(:) = fract(:) - int(fract(:))

    ifracts(:, my_chg_cnt) = int(fract(:)) - order

    call fill_bspline_1_3d(weight, order, &
                           theta(1, 1, my_chg_cnt), dtheta(1, 1, my_chg_cnt))

  end do

  return

end subroutine get_grid_weights

!*******************************************************************************
!
! Subroutine:  claim_blk_recip_imgs
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine claim_blk_recip_imgs(img_cnt, fraction, box, crd_idx_tbl, img_crd, &
                                img_qterm, img_iac, mapped_img_cnt, &
                                mapped_img_lst, img_atm_map, used_img_map, &
                                used_img_cnt, used_img_lst)

  use cit_mod
  use pme_blk_fft_mod
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

  integer               :: j_bot0, j_top1
  integer               :: k_bot0, k_top1
  integer               :: i, j, k
  integer               :: atm_id, img_id
  integer               :: y_bkt_lo, y_bkt_hi
  integer               :: z_bkt_lo, z_bkt_hi
  logical               :: y_wraps
  logical               :: z_wraps
  double precision      :: x_box, y_box, z_box
  double precision      :: y_lo, y_hi, y_cur
  double precision      :: z_lo, z_hi, z_cur

  gbl_used_recip_img_cnt = 0

  if (.not. i_do_recip) return    ! Map null range; not actually ever used...

  x_box = box(1)
  y_box = box(2)
  z_box = box(3)

  ! There will always be 4 or more tasks for block fft-based reciprocal force
  ! calcs.

  k_bot0 = fft_z_offs(my_grid_idx2)
  k_top1 = k_bot0 + fft_z_cnts(my_grid_idx2) + bspl_order - 2
  j_bot0 = fft_y_offs1(my_grid_idx1)
  j_top1 = j_bot0 + fft_y_cnts1(my_grid_idx1) + bspl_order - 2

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

  ! We have to allow for motion of skinnb, plus we have to allow for fact that
  ! atom y crds are rounded down to the grid point.  Hence, the + 1 for y_hi

  y_lo = dble(j_bot0) * y_box / dble(nfft2) - (0.5d0 * skinnb)
  y_hi = dble(j_top1 + 1) * y_box / dble(nfft2) + (0.5d0 * skinnb)

  if (y_lo .lt. 0.d0) then
    y_wraps = .true.
    y_lo = y_lo + y_box
  else if (y_hi .ge. y_box) then
    y_wraps = .true.
    y_hi = y_hi - y_box
  else
    y_wraps = .false.
  end if

  y_bkt_lo = int(y_lo * dble(cit_tbl_y_dim) / y_box)
  y_bkt_hi = int(y_hi * dble(cit_tbl_y_dim) / y_box)

  z_bkt_lo = int(z_lo * dble(cit_tbl_z_dim) / z_box)
  z_bkt_hi = int(z_hi * dble(cit_tbl_z_dim) / z_box)

  if (y_bkt_lo .eq. y_bkt_hi) y_wraps = .false.
  if (z_bkt_lo .eq. z_bkt_hi) z_wraps = .false.

  if (z_wraps) then

    do k = z_bkt_lo, cit_tbl_z_dim - 1

      if (y_wraps) then

        do j = y_bkt_lo, cit_tbl_y_dim - 1
          do i = 0, cit_tbl_x_dim - 1
            do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
              atm_id = img_atm_map(img_id)
              z_cur = fraction(3, atm_id) * z_box
              if (z_cur .ge. z_lo) then
                y_cur = fraction(2, atm_id) * y_box
                if (y_cur .ge. y_lo) then
                  if (used_img_map(img_id) .eq. 0) then
                    used_img_map(img_id) = 1
                    used_img_cnt = used_img_cnt + 1
                    used_img_lst(used_img_cnt) = img_id
                    if (img_iac(img_id) .eq. 0) then
                      img_crd(1, img_id) = fraction(1, atm_id) * x_box
                      img_crd(2, img_id) = y_cur
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
              end if
            end do
          end do
        end do

        do j = 0, y_bkt_hi
          do i = 0, cit_tbl_x_dim - 1
            do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
              atm_id = img_atm_map(img_id)
              z_cur = fraction(3, atm_id) * z_box
              if (z_cur .ge. z_lo) then
                y_cur = fraction(2, atm_id) * y_box
                if (y_cur .lt. y_hi) then
                  if (used_img_map(img_id) .eq. 0) then
                    used_img_map(img_id) = 1
                    used_img_cnt = used_img_cnt + 1
                    used_img_lst(used_img_cnt) = img_id
                    if (img_iac(img_id) .eq. 0) then
                      img_crd(1, img_id) = fraction(1, atm_id) * x_box
                      img_crd(2, img_id) = y_cur
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
              end if
            end do
          end do
        end do

      else

        do j = y_bkt_lo, y_bkt_hi
          do i = 0, cit_tbl_x_dim - 1
            do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
              atm_id = img_atm_map(img_id)
              z_cur = fraction(3, atm_id) * z_box
              if (z_cur .ge. z_lo) then
                y_cur = fraction(2, atm_id) * y_box
                if (y_cur .ge. y_lo .and. y_cur .lt. y_hi) then
                  if (used_img_map(img_id) .eq. 0) then
                    used_img_map(img_id) = 1
                    used_img_cnt = used_img_cnt + 1
                    used_img_lst(used_img_cnt) = img_id
                    if (img_iac(img_id) .eq. 0) then
                      img_crd(1, img_id) = fraction(1, atm_id) * x_box
                      img_crd(2, img_id) = y_cur
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
              end if
            end do
          end do
        end do

      end if

    end do

    do k = 0, z_bkt_hi

      if (y_wraps) then

        do j = y_bkt_lo, cit_tbl_y_dim - 1
          do i = 0, cit_tbl_x_dim - 1
            do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
              atm_id = img_atm_map(img_id)
              z_cur = fraction(3, atm_id) * z_box
              if (z_cur .lt. z_hi) then
                y_cur = fraction(2, atm_id) * y_box
                if (y_cur .ge. y_lo) then
                  if (used_img_map(img_id) .eq. 0) then
                    used_img_map(img_id) = 1
                    used_img_cnt = used_img_cnt + 1
                    used_img_lst(used_img_cnt) = img_id
                    if (img_iac(img_id) .eq. 0) then
                      img_crd(1, img_id) = fraction(1, atm_id) * x_box
                      img_crd(2, img_id) = y_cur
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
              end if
            end do
          end do
        end do

        do j = 0, y_bkt_hi
          do i = 0, cit_tbl_x_dim - 1
            do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
              atm_id = img_atm_map(img_id)
              z_cur = fraction(3, atm_id) * z_box
              if (z_cur .lt. z_hi) then
                y_cur = fraction(2, atm_id) * y_box
                if (y_cur .lt. y_hi) then
                  if (used_img_map(img_id) .eq. 0) then
                    used_img_map(img_id) = 1
                    used_img_cnt = used_img_cnt + 1
                    used_img_lst(used_img_cnt) = img_id
                    if (img_iac(img_id) .eq. 0) then
                      img_crd(1, img_id) = fraction(1, atm_id) * x_box
                      img_crd(2, img_id) = y_cur
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
              end if
            end do
          end do
        end do

      else

        do j = y_bkt_lo, y_bkt_hi
          do i = 0, cit_tbl_x_dim - 1
            do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
              atm_id = img_atm_map(img_id)
              z_cur = fraction(3, atm_id) * z_box
              if (z_cur .lt. z_hi) then
                y_cur = fraction(2, atm_id) * y_box
                if (y_cur .ge. y_lo .and. y_cur .lt. y_hi) then
                  if (used_img_map(img_id) .eq. 0) then
                    used_img_map(img_id) = 1
                    used_img_cnt = used_img_cnt + 1
                    used_img_lst(used_img_cnt) = img_id
                    if (img_iac(img_id) .eq. 0) then
                      img_crd(1, img_id) = fraction(1, atm_id) * x_box
                      img_crd(2, img_id) = y_cur
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
              end if
            end do
          end do
        end do

      end if

    end do

  else

    do k = z_bkt_lo, z_bkt_hi

      if (y_wraps) then

        do j = y_bkt_lo, cit_tbl_y_dim - 1
          do i = 0, cit_tbl_x_dim - 1
            do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
              atm_id = img_atm_map(img_id)
              z_cur = fraction(3, atm_id) * z_box
              if (z_cur .ge. z_lo .and. z_cur .lt. z_hi) then
                y_cur = fraction(2, atm_id) * y_box
                if (y_cur .ge. y_lo) then
                  if (used_img_map(img_id) .eq. 0) then
                    used_img_map(img_id) = 1
                    used_img_cnt = used_img_cnt + 1
                    used_img_lst(used_img_cnt) = img_id
                    if (img_iac(img_id) .eq. 0) then
                      img_crd(1, img_id) = fraction(1, atm_id) * x_box
                      img_crd(2, img_id) = y_cur
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
              end if
            end do
          end do
        end do

        do j = 0, y_bkt_hi
          do i = 0, cit_tbl_x_dim - 1
            do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
              atm_id = img_atm_map(img_id)
              z_cur = fraction(3, atm_id) * z_box
              if (z_cur .ge. z_lo .and. z_cur .lt. z_hi) then
                y_cur = fraction(2, atm_id) * y_box
                if (y_cur .lt. y_hi) then
                  if (used_img_map(img_id) .eq. 0) then
                    used_img_map(img_id) = 1
                    used_img_cnt = used_img_cnt + 1
                    used_img_lst(used_img_cnt) = img_id
                    if (img_iac(img_id) .eq. 0) then
                      img_crd(1, img_id) = fraction(1, atm_id) * x_box
                      img_crd(2, img_id) = y_cur
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
              end if
            end do
          end do
        end do

      else

        do j = y_bkt_lo, y_bkt_hi
          do i = 0, cit_tbl_x_dim - 1
            do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
              atm_id = img_atm_map(img_id)
              z_cur = fraction(3, atm_id) * z_box
              if (z_cur .ge. z_lo .and. z_cur .lt. z_hi) then
                y_cur = fraction(2, atm_id) * y_box
                if (y_cur .ge. y_lo .and. y_cur .lt. y_hi) then
                  if (used_img_map(img_id) .eq. 0) then
                    used_img_map(img_id) = 1
                    used_img_cnt = used_img_cnt + 1
                    used_img_lst(used_img_cnt) = img_id
                    if (img_iac(img_id) .eq. 0) then
                      img_crd(1, img_id) = fraction(1, atm_id) * x_box
                      img_crd(2, img_id) = y_cur
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
              end if
            end do
          end do
        end do

      end if

    end do

  end if

  return

end subroutine claim_blk_recip_imgs

!*******************************************************************************
!
! Subroutine:  claim_blk_recip_imgs_nonorthog
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine claim_blk_recip_imgs_nonorthog(img_cnt, fraction, crd_idx_tbl, &
                                          img_crd, img_qterm, img_iac, &
                                          mapped_img_cnt, mapped_img_lst, &
                                          img_atm_map, used_img_map, &
                                          used_img_cnt, used_img_lst)

  use cit_mod
  use pme_blk_fft_mod
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

  integer               :: j_bot0, j_top1
  integer               :: k_bot0, k_top1
  integer               :: i, j, k
  integer               :: atm_id, img_id
  integer               :: y_bkt_lo, y_bkt_hi
  integer               :: z_bkt_lo, z_bkt_hi
  logical               :: y_wraps
  logical               :: z_wraps
  double precision      :: ucell_stk(3, 3)
  double precision      :: f1, f2, f3
  double precision      :: y_lo, y_hi
  double precision      :: z_lo, z_hi

  ucell_stk(:,:) = ucell(:,:)

  gbl_used_recip_img_cnt = 0

  if (.not. i_do_recip) return    ! Map null range; not actually ever used...
  
  ! There will always be 4 or more tasks for block fft-based reciprocal force
  ! calcs.

  j_bot0 = fft_y_offs1(my_grid_idx1)
  j_top1 = j_bot0 + fft_y_cnts1(my_grid_idx1) + bspl_order - 2
  k_bot0 = fft_z_offs(my_grid_idx2)
  k_top1 = k_bot0 + fft_z_cnts(my_grid_idx2) + bspl_order - 2

  ! We have to allow for motion of skinnb, plus we have to allow for fact that
  ! atom z crds are rounded down to the grid point.  Hence, the + 1 for z_hi

  y_lo = dble(j_bot0) / dble(nfft2) - &
         (0.5d0 * cut_factor(2) * skinnb) / ucell_stk(2,2)      ! a fractional

  y_hi = dble(j_top1 + 1) / dble(nfft2) + &
         (0.5d0 * cut_factor(2) * skinnb) / ucell_stk(2,2)      ! a fractional

  if (y_lo .lt. 0.d0) then
    y_wraps = .true.
    y_lo = y_lo + 1.d0
  else if (y_hi .ge. 1.d0) then
    y_wraps = .true.
    y_hi = y_hi - 1.d0
  else
    y_wraps = .false.
  end if

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

  y_bkt_lo = int(y_lo * dble(cit_tbl_y_dim))
  y_bkt_hi = int(y_hi * dble(cit_tbl_y_dim))

  z_bkt_lo = int(z_lo * dble(cit_tbl_z_dim))
  z_bkt_hi = int(z_hi * dble(cit_tbl_z_dim))

  if (y_bkt_lo .eq. y_bkt_hi) y_wraps = .false.
  if (z_bkt_lo .eq. z_bkt_hi) z_wraps = .false.

  if (z_wraps) then

    do k = z_bkt_lo, cit_tbl_z_dim - 1

      if (y_wraps) then

        do j = y_bkt_lo, cit_tbl_y_dim - 1
          do i = 0, cit_tbl_x_dim - 1
            do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
              atm_id = img_atm_map(img_id)
              f3 = fraction(3, atm_id)
              if (f3 .ge. z_lo) then
                f2 = fraction(2, atm_id)
                if (f2 .ge. y_lo) then
                  if (used_img_map(img_id) .eq. 0) then
                    used_img_map(img_id) = 1
                    used_img_cnt = used_img_cnt + 1
                    used_img_lst(used_img_cnt) = img_id
                    if (img_iac(img_id) .eq. 0) then
                      f1 = fraction(1, atm_id)
                      ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0,
                      ! so we simplify the expression in the critical inner loop
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
              end if
            end do
          end do
        end do

        do j = 0, y_bkt_hi
          do i = 0, cit_tbl_x_dim - 1
            do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
              atm_id = img_atm_map(img_id)
              f3 = fraction(3, atm_id)
              if (f3 .ge. z_lo) then
                f2 = fraction(2, atm_id)
                if (f2 .lt. y_hi) then
                  if (used_img_map(img_id) .eq. 0) then
                    used_img_map(img_id) = 1
                    used_img_cnt = used_img_cnt + 1
                    used_img_lst(used_img_cnt) = img_id
                    if (img_iac(img_id) .eq. 0) then
                      f1 = fraction(1, atm_id)
                      ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0,
                      ! so we simplify the expression in the critical inner loop
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
              end if
            end do
          end do
        end do

      else

        do j = y_bkt_lo, y_bkt_hi
          do i = 0, cit_tbl_x_dim - 1
            do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
              atm_id = img_atm_map(img_id)
              f3 = fraction(3, atm_id)
              if (f3 .ge. z_lo) then
                f2 = fraction(2, atm_id)
                if (f2 .ge. y_lo .and. f2 .lt. y_hi) then
                  if (used_img_map(img_id) .eq. 0) then
                    used_img_map(img_id) = 1
                    used_img_cnt = used_img_cnt + 1
                    used_img_lst(used_img_cnt) = img_id
                    if (img_iac(img_id) .eq. 0) then
                      f1 = fraction(1, atm_id)
                      ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0,
                      ! so we simplify the expression in the critical inner loop
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
              end if
            end do
          end do
        end do

      end if

    end do

    do k = 0, z_bkt_hi

      if (y_wraps) then

        do j = y_bkt_lo, cit_tbl_y_dim - 1
          do i = 0, cit_tbl_x_dim - 1
            do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
              atm_id = img_atm_map(img_id)
              f3 = fraction(3, atm_id)
              if (f3 .lt. z_hi) then
                f2 = fraction(2, atm_id)
                if (f2 .ge. y_lo) then
                  if (used_img_map(img_id) .eq. 0) then
                    used_img_map(img_id) = 1
                    used_img_cnt = used_img_cnt + 1
                    used_img_lst(used_img_cnt) = img_id
                    if (img_iac(img_id) .eq. 0) then
                      f1 = fraction(1, atm_id)
                      ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0,
                      ! so we simplify the expression in the critical inner loop
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
              end if
            end do
          end do
        end do

        do j = 0, y_bkt_hi
          do i = 0, cit_tbl_x_dim - 1
            do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
              atm_id = img_atm_map(img_id)
              f3 = fraction(3, atm_id)
              if (f3 .lt. z_hi) then
                f2 = fraction(2, atm_id)
                if (f2 .lt. y_hi) then
                  if (used_img_map(img_id) .eq. 0) then
                    used_img_map(img_id) = 1
                    used_img_cnt = used_img_cnt + 1
                    used_img_lst(used_img_cnt) = img_id
                    if (img_iac(img_id) .eq. 0) then
                      f1 = fraction(1, atm_id)
                      ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0,
                      ! so we simplify the expression in the critical inner loop
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
              end if
            end do
          end do
        end do

      else

        do j = y_bkt_lo, y_bkt_hi
          do i = 0, cit_tbl_x_dim - 1
            do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
              atm_id = img_atm_map(img_id)
              f3 = fraction(3, atm_id)
              if (f3 .lt. z_hi) then
                f2 = fraction(2, atm_id)
                if (f2 .ge. y_lo .and. f2 .lt. y_hi) then
                  if (used_img_map(img_id) .eq. 0) then
                    used_img_map(img_id) = 1
                    used_img_cnt = used_img_cnt + 1
                    used_img_lst(used_img_cnt) = img_id
                    if (img_iac(img_id) .eq. 0) then
                      f1 = fraction(1, atm_id)
                      ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0,
                      ! so we simplify the expression in the critical inner loop
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
              end if
            end do
          end do
        end do

      end if

    end do

  else

    do k = z_bkt_lo, z_bkt_hi

      if (y_wraps) then

        do j = y_bkt_lo, cit_tbl_y_dim - 1
          do i = 0, cit_tbl_x_dim - 1
            do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
              atm_id = img_atm_map(img_id)
              f3 = fraction(3, atm_id)
              if (f3 .ge. z_lo .and. f3 .lt. z_hi) then
                f2 = fraction(2, atm_id)
                if (f2 .ge. y_lo) then
                  if (used_img_map(img_id) .eq. 0) then
                    used_img_map(img_id) = 1
                    used_img_cnt = used_img_cnt + 1
                    used_img_lst(used_img_cnt) = img_id
                    if (img_iac(img_id) .eq. 0) then
                      f1 = fraction(1, atm_id)
                      ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0,
                      ! so we simplify the expression in the critical inner loop
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
              end if
            end do
          end do
        end do

        do j = 0, y_bkt_hi
          do i = 0, cit_tbl_x_dim - 1
            do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
              atm_id = img_atm_map(img_id)
              f3 = fraction(3, atm_id)
              if (f3 .ge. z_lo .and. f3 .lt. z_hi) then
                f2 = fraction(2, atm_id)
                if (f2 .lt. y_hi) then
                  if (used_img_map(img_id) .eq. 0) then
                    used_img_map(img_id) = 1
                    used_img_cnt = used_img_cnt + 1
                    used_img_lst(used_img_cnt) = img_id
                    if (img_iac(img_id) .eq. 0) then
                      f1 = fraction(1, atm_id)
                      ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0,
                      ! so we simplify the expression in the critical inner loop
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
              end if
            end do
          end do
        end do

      else

        do j = y_bkt_lo, y_bkt_hi
          do i = 0, cit_tbl_x_dim - 1
            do img_id = crd_idx_tbl(i,j,k)%img_lo, crd_idx_tbl(i,j,k)%img_hi
              atm_id = img_atm_map(img_id)
              f3 = fraction(3, atm_id)
              if (f3 .ge. z_lo .and. f3 .lt. z_hi) then
                f2 = fraction(2, atm_id)
                if (f2 .ge. y_lo .and. f2 .lt. y_hi) then
                  if (used_img_map(img_id) .eq. 0) then
                    used_img_map(img_id) = 1
                    used_img_cnt = used_img_cnt + 1
                    used_img_lst(used_img_cnt) = img_id
                    if (img_iac(img_id) .eq. 0) then
                      f1 = fraction(1, atm_id)
                      ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0,
                      ! so we simplify the expression in the critical inner loop
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
              end if
            end do
          end do
        end do

      end if

    end do

  end if

  return


end subroutine claim_blk_recip_imgs_nonorthog

!*******************************************************************************
!
! Subroutine:  get_grid_weights_nonorthog
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine get_grid_weights_nonorthog(img_cnt, crd, img_atm_map, &
                                      ifracts, theta, dtheta, order, &
                                      used_recip_img_cnt, used_recip_img_lst, &
                                      my_chgs, my_chg_cnt)

  use bspline_mod
  use gbl_datatypes_mod
  use pme_blk_fft_mod
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
  integer               :: used_recip_img_cnt
  integer               :: used_recip_img_lst(*)
  integer               :: my_chgs(*), my_chg_cnt

! Local variables:

  double precision      :: crd_x, crd_y, crd_z
  double precision      :: fract(3)
  double precision      :: weight(3)
  double precision      :: recip_stk(3, 3)
  integer               :: atm_id
  integer               :: i, k

  integer               :: kbot0, ktop1
  logical               :: kwraps
  integer               :: lst_idx

  my_chg_cnt = 0

  kbot0 = fft_z_offs(my_grid_idx2)
  ktop1 = kbot0 + fft_z_cnts(my_grid_idx2) + order - 2
  kwraps = (ktop1 .ge. nfft3)

  recip_stk(:,:) = recip(:,:)

  do lst_idx = 1, used_recip_img_cnt
    i = used_recip_img_lst(lst_idx)

! Unfortunately we need fractional coords that are not available without going
! back to atm_crd().  This is the case because the coordinates in img_crd
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

    k = int(fract(3))

    if (kwraps) then
      if (k .lt. kbot0 .and. k .gt. ktop1 - nfft3) cycle
    else
      if (k .lt. kbot0 .or. k .gt. ktop1) cycle
    end if

    my_chg_cnt = my_chg_cnt + 1
    my_chgs(my_chg_cnt) = i

    weight(:) = fract(:) - int(fract(:))

    ifracts(:, my_chg_cnt) = int(fract(:)) - order

    call fill_bspline_1_3d(weight, order, &
                           theta(1, 1, my_chg_cnt), dtheta(1, 1, my_chg_cnt))

  end do

  return

end subroutine get_grid_weights_nonorthog

#endif /* MPI */

end module pme_blk_recip_mod
