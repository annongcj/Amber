#include "copyright.i"

!*******************************************************************************
!
! Module: barostats_mod
!
! Description: Monte-Carlo barostat
!              
!*******************************************************************************

module barostats_mod

  use random_mod, only : random_state, amrset_gen, amrand_gen
#ifdef CUDA
  use pbc_mod
#endif

  implicit none

  private

  double precision, save :: dvmax = 0.02d0

  integer, save :: total_mcbar_attempts = 0
  integer, save :: total_mcbar_successes = 0

  integer, save      :: mcbar_attempts = 0
  integer, save      :: mcbar_successes = 0
  integer, parameter :: dvmax_interval = 10

  type(random_state), save :: mcbar_gen

! Variable explanations
!
! Public variables
!   dvmax           : Size of trial dV move
!   total_mcbar_attempts  : # of trial moves we've attempted with MC barostat
!   total_mcbar_successes : # of trial moves we've accepted with MC barostat
!
!   mcbar_attempts  : # of trial moves; used to adjust size of volume move
!   mcbar_successes : # of successes; used to adjust size of volume move
!   dvmax_interval  : # of exchange attempts to do before deciding whether or
!                     not to change the size of the volume move

  public mcbar_trial, mcbar_setup, mcbar_summary

contains

!*******************************************************************************
!
! Subroutine:  mcbar_setup
!
! Description: Sets up the random number generator for the MC barostat
!              
!*******************************************************************************

#ifdef MPI
! Make sure all threads on the same replica have the exact same random # stream
subroutine mcbar_setup(ig)

  use parallel_dat_mod, only : master, err_code_mpi, pmemd_comm

  implicit none
  include 'mpif.h'
  integer, intent(in) :: ig
  integer :: local_ig

!RCW: huh? - Where does the number 10 come from here? Was it just picked because
!it is a nice round number or is there solid reasoning behind this?
  local_ig = ig + 10

  call mpi_bcast(local_ig, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)

  call amrset_gen(mcbar_gen, local_ig)

  return

end subroutine mcbar_setup
#else
! Non-MPI case below
subroutine mcbar_setup(ig)

  implicit none

  integer, intent(in) :: ig

!RCW - same again - why was 10 chosen here?
  call amrset_gen(mcbar_gen, ig + 10)

  return

end subroutine mcbar_setup
#endif

!*******************************************************************************
!
! Subroutine:  mcbar_trial
!
! Description: Perform the trial move for the MC barostat
!              
!*******************************************************************************

#ifdef MPI
subroutine mcbar_trial(atm_cnt, crd, frc, mass, my_atm_lst, new_list, &
                       orig_pot_ene, verbose, pme_err_est, vel, numextra)
#else
subroutine mcbar_trial(atm_cnt, crd, frc, mass, new_list, orig_pot_ene, &
                       verbose, pme_err_est)
#endif

  use cit_mod, only : set_cit_tbl_dims
  use file_io_dat_mod, only : mdout
  use img_mod, only : gbl_img_atm_map, gbl_atm_img_map
  use mdin_ctrl_dat_mod, only : vdw_cutoff, ntp, temp0, pres0, nmropt, &
                                csurften, gamma_ten, ninterface, usemidpoint, infe, &
                                baroscalingdir
  use nmr_calls_mod, only : nmrdcp, skip_print
  use nfe_lib_mod, only : nfe_real_mdstep

#ifdef MPI
  use extra_pnts_nb14_mod, only : ep_frames, gbl_frame_cnt, zero_extra_pnts_vec
  use loadbal_mod, only : atm_redist_needed
  use parallel_mod, only : mpi_allgathervec
  use parallel_processor_mod
  use bonds_midpoint_mod
  use angles_midpoint_mod
  use dihedrals_midpoint_mod
  use pme_direct_mod
  use nb_exclusions_mod
#endif
  use parallel_dat_mod, only : master
  use pbc_mod, only : uc_volume, ucell, cut_factor, pbc_box, &
                      pressure_scale_pbc_data
#ifdef MPI
  use pme_force_mod, only : pme_pot_ene_rec, pme_force, pme_force_midpoint
#else
  use pme_force_mod, only : pme_pot_ene_rec, pme_force
#endif
  use mdin_ewald_dat_mod, only : skinnb
  use mdin_ctrl_dat_mod, only : ntf, ntc
  use mol_list_mod, only : gbl_mol_cnt
  use processor_mod
  use parallel_processor_mod

  implicit none

! Passed parameters
  
  integer, intent(in)     :: atm_cnt
  integer, intent(in)     :: verbose
  logical, intent(in out) :: new_list
#ifdef MPI
  integer, intent(in)     :: numextra
  integer, intent(in)     :: my_atm_lst(*)

  double precision, intent(in out)      :: vel(3, atm_cnt)
#endif
  double precision, intent(in out)      :: crd(3, atm_cnt)
  double precision, intent(in out)      :: frc(3, atm_cnt)
  double precision, intent(in)          :: mass(atm_cnt)
  double precision, intent(in out)      :: pme_err_est

  type(pme_pot_ene_rec), intent(in out) :: orig_pot_ene

! Local parameters

  double precision, dimension(3) :: dv
  double precision, dimension(3) :: rmu
  double precision, dimension(3) :: virial
  double precision, dimension(3) :: ekcmt
  double precision, dimension(3, atm_cnt) :: tmp_frc
#ifdef MPI
  double precision, dimension(3,proc_num_atms+proc_ghost_num_atms) :: frc_state
  double precision, dimension(3,proc_num_atms+proc_ghost_num_atms) :: crd_state
  double precision, dimension(3,proc_num_atms+proc_ghost_num_atms) :: vel_state
  double precision, dimension(3,proc_num_atms+proc_ghost_num_atms) :: last_vel_state
  double precision, dimension(proc_num_atms+proc_ghost_num_atms) :: mass_state
  double precision, dimension(proc_num_atms+proc_ghost_num_atms) :: qterm_state
  integer, dimension(proc_num_atms+proc_ghost_num_atms) :: full_list_state
  integer, dimension(proc_num_atms+proc_ghost_num_atms) :: iac_state
  integer, dimension(3,proc_num_atms+proc_ghost_num_atms) :: wrap_state
  integer :: proc_num_atms_old, proc_num_atms_ghost_old
  integer :: i, j, atmid
  type(bond_rec)     :: h_bond_state(cit_nbonh)
  type(bond_rec)     :: a_bond_state(cit_nbona)
  type(angle_rec)    :: h_angle_state(cit_ntheth)
  type(angle_rec)    :: a_angle_state(cit_ntheta)
  type(dihed_rec)    :: h_dihed_state(cit_nphih)
  type(dihed_rec)    :: a_dihed_state(cit_nphia)
  integer            :: cit_nbonh_old, cit_nbona_old, cit_ntheth_old, &
                        cit_ntheta_old, cit_nphih_old, cit_nphia_old
#endif

  double precision :: randval
  double precision :: pv_work
  double precision :: delta_area
  double precision :: orig_vol
  double precision :: local_pme_err_est
  double precision :: expfac
  double precision :: nbeta

  logical :: scaled_new_list

  integer :: aniso_dim

  type(pme_pot_ene_rec) :: new_pot_ene_rec

! Constants
  double precision, parameter :: ONE_THIRD = 1.d0 / 3.d0
  double precision, parameter :: ONE_HALF = 0.5d0
  ! converts dyne.A^2/cm to kcal/mol (2.390057e-27 * 6.022e23)
  double precision, parameter :: TENSION_CONV = 0.0014393264316443592

  ! This is another mcbar attempt
  total_mcbar_attempts = total_mcbar_attempts + 1
  mcbar_attempts = mcbar_attempts + 1

  ! Back up the original volume

  orig_vol = uc_volume

  ! Get the dV move we plan on doing

  if (ntp .eq. 1) then
    ! Isotropic
    call amrand_gen(mcbar_gen, randval)
    dv(1) = (randval - ONE_HALF) * dvmax
    dv(2) = dv(1)
    dv(3) = dv(1)
  else if (ntp .eq. 2) then
    ! Anisotropic -- pick one dimension to change(if baroscalingdir==0)
    if (baroscalingdir .eq. 0) then
    ! 'baroscalingdir' -- controlling the direction for anisotropic
    ! pressure scaling. Choose one direction(x, y or z), the box only
    ! scales along the chosen direction, its lengths along the other two
    ! directions are fixed. Default value 0, randomly pick x, y or
    ! z-direction to apply box size change. 'baroscalingdir = 1'-->
    ! only x-direction is allowed to scale; '= 2'--> y-direction;
    ! '= 3'--> z-direction. --Yeyue Xiong
      call amrand_gen(mcbar_gen, randval)
      aniso_dim = int(randval * 3.d0 * 0.99999999d0) + 1
      call amrand_gen(mcbar_gen, randval)
      dv(:) = 0.d0
      dv(aniso_dim) = (randval - ONE_HALF) * dvmax
    else
      select case (baroscalingdir)
        case(1)
          aniso_dim = 1
          call amrand_gen(mcbar_gen, randval)
          dv(:) = 0.d0
          dv(aniso_dim) = (randval - ONE_HALF) * dvmax
        case(2)
          aniso_dim = 2
          call amrand_gen(mcbar_gen, randval)
          dv(:) = 0.d0
          dv(aniso_dim) = (randval - ONE_HALF) * dvmax
        case(3)
          aniso_dim = 3
          call amrand_gen(mcbar_gen, randval)
          dv(:) = 0.d0
          dv(aniso_dim) = (randval - ONE_HALF) * dvmax
      end select
    end if
  else if (ntp .eq. 3) then
    ! Semi-isotropic -- pick one dimension to change (X,Y count as 1 dim)
    call amrand_gen(mcbar_gen, randval)
    aniso_dim = int(randval * 2.d0 * 0.99999999d0) + 1
    call amrand_gen(mcbar_gen, randval)
    if (aniso_dim .eq. 1) then
      select case (csurften)
        case(1)
          dv(2) = (randval - ONE_HALF) * dvmax
          dv(3) = dv(2)
          dv(1) = 0.d0
        case(2)
          dv(1) = (randval - ONE_HALF) * dvmax
          dv(3) = dv(1)
          dv(2) = 0.d0
        case(3)
          dv(1) = (randval - ONE_HALF) * dvmax
          dv(2) = dv(1)
          dv(3) = 0.d0
      end select
    else
      select case (csurften)
        case(1)
          dv(1) = (randval - ONE_HALF) * dvmax
          dv(2) = 0.d0
          dv(3) = 0.d0
        case(2)
          dv(2) = (randval - ONE_HALF) * dvmax
          dv(1) = 0.d0
          dv(3) = 0.d0
        case(3)
          dv(3) = (randval - ONE_HALF) * dvmax
          dv(1) = 0.d0
          dv(2) = 0.d0
      end select
    end if
  end if

  if (csurften .gt. 0) then
    select case (csurften)
      case(1)
        delta_area = ucell(2,2) * ucell(3,3)
      case(2)
        delta_area = ucell(1,1) * ucell(3,3)
      case(3)
        delta_area = ucell(1,1) * ucell(2,2)
    end select
  end if

  rmu(1) = (1.d0 + dv(1)) ** ONE_THIRD
  rmu(2) = (1.d0 + dv(2)) ** ONE_THIRD
  rmu(3) = (1.d0 + dv(3)) ** ONE_THIRD
#ifdef MPI
  if(usemidpoint) then
    !call comm_3dimensions_3dbls(proc_atm_vel)
    !call comm_3dimensions_3dbls(proc_atm_last_vel)
    ! Save current state
    crd_state(:,:)=proc_atm_crd(:,1:proc_num_atms+proc_ghost_num_atms)
    frc_state(:,:)=proc_atm_frc(:,1:proc_num_atms+proc_ghost_num_atms)
    vel_state(:,:)=proc_atm_vel(:,1:proc_num_atms+proc_ghost_num_atms)
    mass_state(:)=proc_atm_mass(1:proc_num_atms+proc_ghost_num_atms)
    qterm_state(:)=proc_atm_qterm(1:proc_num_atms+proc_ghost_num_atms)
    last_vel_state(:,:)=proc_atm_last_vel(:,1:proc_num_atms+proc_ghost_num_atms)
    full_list_state(:)=proc_atm_to_full_list(1:proc_num_atms+proc_ghost_num_atms)
    iac_state(:)=proc_iac(1:proc_num_atms+proc_ghost_num_atms)
    wrap_state(:,:)=proc_atm_wrap(:,1:proc_num_atms+proc_ghost_num_atms)
    proc_num_atms_old = proc_num_atms
    proc_num_atms_ghost_old = proc_ghost_num_atms
    do j = 1, cit_nbonh
      h_bond_state(j)%atm_i = cit_h_bond(j)%atm_i
      h_bond_state(j)%atm_j = cit_h_bond(j)%atm_j
      h_bond_state(j)%parm_idx = cit_h_bond(j)%parm_idx
    end do
    do j = 1, cit_nbona
      a_bond_state(j)%atm_i = cit_a_bond(j)%atm_i
      a_bond_state(j)%atm_j = cit_a_bond(j)%atm_j
      a_bond_state(j)%parm_idx = cit_a_bond(j)%parm_idx
    end do
    do j = 1, cit_ntheth
      h_angle_state(j)%atm_i = cit_h_angle(j)%atm_i
      h_angle_state(j)%atm_j = cit_h_angle(j)%atm_j
      h_angle_state(j)%atm_k = cit_h_angle(j)%atm_k
      h_angle_state(j)%parm_idx = cit_h_angle(j)%parm_idx
    end do
    do j = 1, cit_ntheta
      a_angle_state(j)%atm_i = cit_a_angle(j)%atm_i
      a_angle_state(j)%atm_j = cit_a_angle(j)%atm_j
      a_angle_state(j)%atm_k = cit_a_angle(j)%atm_k
      a_angle_state(j)%parm_idx = cit_a_angle(j)%parm_idx
    end do
    do j = 1, cit_nphih
      h_dihed_state(j)%atm_i = cit_h_dihed(j)%atm_i
      h_dihed_state(j)%atm_j = cit_h_dihed(j)%atm_j
      h_dihed_state(j)%atm_k = cit_h_dihed(j)%atm_k
      h_dihed_state(j)%atm_l = cit_h_dihed(j)%atm_l
      h_dihed_state(j)%parm_idx = cit_h_dihed(j)%parm_idx
    end do
    do j = 1, cit_nphia
      a_dihed_state(j)%atm_i = cit_a_dihed(j)%atm_i
      a_dihed_state(j)%atm_j = cit_a_dihed(j)%atm_j
      a_dihed_state(j)%atm_k = cit_a_dihed(j)%atm_k
      a_dihed_state(j)%atm_l = cit_a_dihed(j)%atm_l
      a_dihed_state(j)%parm_idx = cit_a_dihed(j)%parm_idx
    end do
    cit_nbonh_old = cit_nbonh
    cit_nbona_old = cit_nbona
    cit_ntheth_old = cit_ntheth
    cit_ntheta_old = cit_ntheta
    cit_nphih_old = cit_nphih
    cit_nphia_old = cit_nphia
  end if
#endif
if(usemidpoint) then
#ifdef MPI
!    call comm_3dimensions_3ints(proc_atm_wrap, .true.)
!    call comm_3dimensions_3dbls(proc_atm_crd, .true.)
#endif
endif
#ifdef MPI
  call scale_system_volume(rmu, verbose, atm_cnt, crd, mass, scaled_new_list, &
                           my_atm_lst)
#else
  call scale_system_volume(rmu, verbose, atm_cnt, crd, mass, scaled_new_list)
#endif
  new_list = new_list .or. scaled_new_list

#ifdef MPI
if(.not. usemidpoint) then
  ! redistribute coordinates if we're rebuilding the list. (always do this for
  ! now since something appears to be wrong in parallel...)
  ! TODO: Figure out how to avoid global redistributions every time we attempt a
  !       volume change
! if (new_list) then
    call mpi_allgathervec(atm_cnt, crd)
!   if (atm_redist_needed) then
      ! Also need velocities and forces, since we haven't propagated dynamics
      ! yet and the next force call will automagically redistribute atoms...
      call mpi_allgathervec(atm_cnt, vel)
      call mpi_allgathervec(atm_cnt, frc)
      if (numextra .gt. 0) &
        call zero_extra_pnts_vec(vel, ep_frames, gbl_frame_cnt)
!   end if
! end if
end if !usemidpoint
#endif

#ifdef CUDA
  call gpu_download_frc(frc)
#endif
  ! Now we need to calculate our next energy, but do not print any nmropt values

  if (nmropt .gt. 0) skip_print = .true.
  if (infe .gt. 0) nfe_real_mdstep = .false.
#ifdef MPI
!Last argument here is -1 since this is nstep and is not relevant for barostat
!call.
if(usemidpoint) then
  call pme_force_midpoint(atm_cnt, crd, tmp_frc, &
                 my_atm_lst, .true., .true., .false., &
                 new_pot_ene_rec, -1, virial, ekcmt, local_pme_err_est)
else
  call pme_force(atm_cnt, crd, tmp_frc, gbl_img_atm_map, gbl_atm_img_map, &
                 my_atm_lst, new_list, .true., .false., &
                 new_pot_ene_rec, -1, virial, ekcmt, local_pme_err_est)
end if
#else
  call pme_force(atm_cnt, crd, tmp_frc, gbl_img_atm_map, gbl_atm_img_map, &
                 new_list, .true., .false., &
                 new_pot_ene_rec, -1, virial, ekcmt, local_pme_err_est)
#endif /* MPI */

  ! Decrement the NMR opt counter since this force call should not 'count'

  if (nmropt .gt. 0) call nmrdcp

  ! p*dV (6.02204d-2 / 4184.0d0 converts from bar*A^3/particle to kcal/mol)
  pv_work = pres0 * (uc_volume - orig_vol) * 6.02204d-2 / 4184.0d0

  if (csurften .gt. 0) then
    select case (csurften)
      case(1)
        delta_area = ucell(2,2) * ucell(3,3) - delta_area
      case(2)
        delta_area = ucell(1,1) * ucell(3,3) - delta_area
      case(3)
        delta_area = ucell(1,1) * ucell(2,2) - delta_area
    end select
    pv_work = pv_work - gamma_ten * ninterface * delta_area * TENSION_CONV
  end if

  nbeta = -4.184d0 / (8.31441d-3 * temp0)

  expfac = new_pot_ene_rec%total - orig_pot_ene%total + pv_work + &
           gbl_mol_cnt * log(rmu(1) * rmu(2) * rmu(3)) / nbeta

  call amrand_gen(mcbar_gen, randval)

  ! Monte carlo decision

  ! DEBUG
  !if (master) &
  !  write(mdout, '(a)', advance='NO') '| Attempting MC barostat volume change: '
   if (randval .lt. exp(nbeta * expfac)) then
    ! Accept -- update force array
if(.not. usemidpoint) then
    frc(1,:) = tmp_frc(1,:)
    frc(2,:) = tmp_frc(2,:)
    frc(3,:) = tmp_frc(3,:)
end if
#ifdef MPI
    reorient_flag = .false.
#endif
    pme_err_est = local_pme_err_est
    total_mcbar_successes = total_mcbar_successes + 1
    mcbar_successes = mcbar_successes + 1
    orig_pot_ene = new_pot_ene_rec
    ! DEBUG
!   if(master) write(mdout, '(a)') 'Accepted'
  else
    ! Reject -- rescale everything and put back original forces
    rmu(1) = 1.d0 / rmu(1)
    rmu(2) = 1.d0 / rmu(2)
    rmu(3) = 1.d0 / rmu(3)
#ifdef MPI
    call scale_system_volume(rmu, verbose, atm_cnt, crd, mass, scaled_new_list, &
                             my_atm_lst)
    if(usemidpoint) then
        ! Load old state
        proc_num_atms = proc_num_atms_old
        proc_ghost_num_atms = proc_num_atms_ghost_old
        proc_atm_crd(:,1:proc_num_atms+proc_ghost_num_atms) = crd_state(:,:)
        proc_atm_frc(:,1:proc_num_atms+proc_ghost_num_atms) = frc_state(:,:)
        proc_atm_vel(:,1:proc_num_atms+proc_ghost_num_atms) = vel_state(:,:)
        proc_atm_mass(1:proc_num_atms+proc_ghost_num_atms) = mass_state(:)
        proc_atm_qterm(1:proc_num_atms+proc_ghost_num_atms) = qterm_state(:)
        proc_atm_last_vel(:,1:proc_num_atms+proc_ghost_num_atms) = last_vel_state(:,:)
        proc_atm_to_full_list(1:proc_num_atms+proc_ghost_num_atms) = full_list_state(:)
        proc_iac(1:proc_num_atms+proc_ghost_num_atms) = iac_state(:)
        proc_atm_wrap(:,1:proc_num_atms+proc_ghost_num_atms) = wrap_state(:,:)
        reorient_flag = .true.

        cit_nbonh = cit_nbonh_old
        cit_nbona = cit_nbona_old
        cit_ntheth = cit_ntheth_old
        cit_ntheta =cit_ntheta_old
        cit_nphih = cit_nphih_old
        cit_nphia = cit_nphia_old
        do j = 1, cit_nbonh
          cit_h_bond(j)%atm_i = h_bond_state(j)%atm_i
          cit_h_bond(j)%atm_j = h_bond_state(j)%atm_j
          cit_h_bond(j)%parm_idx = h_bond_state(j)%parm_idx
        end do
        do j = 1, cit_nbona
          cit_a_bond(j)%atm_i = a_bond_state(j)%atm_i
          cit_a_bond(j)%atm_j = a_bond_state(j)%atm_j
          cit_a_bond(j)%parm_idx = a_bond_state(j)%parm_idx
        end do
        do j = 1, cit_ntheth
          cit_h_angle(j)%atm_i = h_angle_state(j)%atm_i
          cit_h_angle(j)%atm_j = h_angle_state(j)%atm_j
          cit_h_angle(j)%atm_k = h_angle_state(j)%atm_k
          cit_h_angle(j)%parm_idx = h_angle_state(j)%parm_idx
        end do
        do j = 1, cit_ntheta
          cit_a_angle(j)%atm_i = a_angle_state(j)%atm_i
          cit_a_angle(j)%atm_j = a_angle_state(j)%atm_j
          cit_a_angle(j)%atm_k = a_angle_state(j)%atm_k
          cit_a_angle(j)%parm_idx = a_angle_state(j)%parm_idx
        end do
        do j = 1, cit_nphih
          cit_h_dihed(j)%atm_i = h_dihed_state(j)%atm_i 
          cit_h_dihed(j)%atm_j = h_dihed_state(j)%atm_j
          cit_h_dihed(j)%atm_k = h_dihed_state(j)%atm_k
          cit_h_dihed(j)%atm_l = h_dihed_state(j)%atm_l
          cit_h_dihed(j)%parm_idx = h_dihed_state(j)%parm_idx
        end do
        do j = 1, cit_nphia
          cit_a_dihed(j)%atm_i = a_dihed_state(j)%atm_i
          cit_a_dihed(j)%atm_j = a_dihed_state(j)%atm_j
          cit_a_dihed(j)%atm_k = a_dihed_state(j)%atm_k
          cit_a_dihed(j)%atm_l = a_dihed_state(j)%atm_l
          cit_a_dihed(j)%parm_idx = a_dihed_state(j)%parm_idx
        end do
        ! nstep is really only used for output so we set it to -5 so output is
        ! obvious where code breaks
        call pme_list_midpoint(atm_cnt, crd, atm_nb_maskdata, atm_nb_mask &
                      , proc_atm_nb_maskdata, proc_atm_nb_mask,5)
        
    end if ! end of load state if use midpoint
#else
    call scale_system_volume(rmu, verbose, atm_cnt, crd, mass, scaled_new_list)
#endif
    new_list = new_list .or. scaled_new_list
#ifdef CUDA
    ! For GPUs, reload the old forces back to the GPU
    call gpu_upload_frc(frc)
#endif
    ! DEBUG
!   if(master) write(mdout, '(a)') 'Rejected'
  end if

  if (mod(mcbar_attempts, dvmax_interval) .eq. 0) then

    ! If our success fraction is too large or too small, adjust dvmax

    if (mcbar_successes .ge. 0.75d0 * mcbar_attempts) then

      dvmax = dvmax * 1.1d0
      mcbar_attempts = 0
      mcbar_successes = 0
      if (master) &
        write(mdout, '(a)') '| MC Barostat: Increasing size of volume moves'

    else if (mcbar_successes .le. 0.25d0 * mcbar_attempts) then

      dvmax = dvmax / 1.1d0
      mcbar_attempts = 0
      mcbar_successes = 0
      if (master) &
        write(mdout, '(a)') '| MC Barostat: Decreasing size of volume moves'

    end if

  end if

end subroutine mcbar_trial

!*******************************************************************************
!
! Subroutine:  scale_system_volume
!
! Description: Scales the system volume and checks if a new pairlist is needed
!              
!*******************************************************************************

#ifdef MPI
subroutine scale_system_volume(rmu, verbose, atm_cnt, crd, mass, new_list, &
                               my_atm_lst)
#else
subroutine scale_system_volume(rmu, verbose, atm_cnt, crd, mass, new_list)
#endif

  use cit_mod, only : set_cit_tbl_dims, bkt_size
  use constraints_mod, only : natc, atm_xc
  use dynamics_dat_mod, only : gbl_mol_mass_inv, gbl_my_mol_lst, gbl_mol_com
  use mdin_ctrl_dat_mod, only : vdw_cutoff, ntp, ntr, usemidpoint
  use mdin_ewald_dat_mod, only : skinnb
#ifdef MPI
  use nb_pairlist_mod, only : gbl_atm_saved_crd, check_my_atom_movement, &
                              proc_check_my_atom_movement
  use processor_mod
  use parallel_processor_mod
#else
  use nb_pairlist_mod, only : gbl_atm_saved_crd, check_all_atom_movement
#endif
  use parallel_dat_mod ! for MPI functions and MPI-related variables
  use pbc_mod, only : pressure_scale_pbc_data, pressure_scale_crds, &
                      pressure_scale_restraint_crds, &
#ifdef MPI
                      pressure_scale_crds_midpoint, &
#endif
                      pbc_box, cut_factor
  use mdin_ctrl_dat_mod, only : usemidpoint
  implicit none

! Formal arguments
  
  integer, intent(in)  :: atm_cnt
  integer, intent(in)  :: verbose
#ifdef MPI
  integer, intent(in)  :: my_atm_lst(*)
#endif
  logical, intent(out) :: new_list

  double precision, intent(in out) :: crd(3, atm_cnt)
  double precision, intent(in)     :: mass(atm_cnt)
  double precision, intent(in)     :: rmu(3)

  integer, dimension(2) :: new_list_buf

  new_list = .false.
  new_list_buf(:) = 0

! Local variables
  
  call pressure_scale_pbc_data(rmu, vdw_cutoff + skinnb, verbose)
#ifdef CUDA
  call gpu_pressure_scale(ucell, recip, uc_volume)
#else
  call set_cit_tbl_dims(pbc_box, vdw_cutoff + skinnb, cut_factor)
#endif

#ifdef MPI
  if(usemidpoint) then
    call pressure_scale_crds_midpoint(proc_atm_crd, proc_atm_mass, gbl_mol_mass_inv, gbl_mol_com)
    reorient_flag = .true.
  else
    call pressure_scale_crds(crd, mass, gbl_mol_mass_inv, gbl_my_mol_lst, gbl_mol_com)
  end if
#else
  call pressure_scale_crds(crd, mass, gbl_mol_mass_inv, gbl_mol_com)
#endif

#ifndef CUDA
  ! If we have restraints, we need to scale restraint coordinates
  if (ntr .gt. 0 .and. natc .gt. 0) then
#ifdef MPI
    call pressure_scale_restraint_crds(atm_xc, mass, gbl_mol_mass_inv, &
                                       gbl_my_mol_lst)
#else
    call pressure_scale_restraint_crds(atm_xc, mass, gbl_mol_mass_inv)
#endif

if(usemidpoint) then
#ifdef MPI
    call comm_3dimensions_3ints(proc_atm_wrap, .true.)
    call comm_3dimensions_3dbls(proc_atm_crd, .true.)
#endif
endif

  end if
#endif

  ! Now we have to see if we need to rebuild our pairlist

#ifdef CUDA
  call gpu_skin_test()
#else

#ifdef MPI
  if(usemidpoint) then
    call proc_check_my_atom_movement(proc_atm_crd, proc_saved_atm_crd, &
                                     skinnb, ntp, new_list)
  else
    call check_my_atom_movement(crd, gbl_atm_saved_crd, my_atm_lst, &
                              skinnb, ntp, new_list)
  end if
  if (new_list) &
    new_list_buf(1) = 1
  
  ! Make sure everyone knows if we need to rebuild a list according to any
  ! worker thread
  call mpi_allreduce(new_list_buf(1), new_list_buf(2), 1, mpi_integer, &
                     mpi_sum, pmemd_comm, err_code_mpi)
  
  new_list = new_list_buf(2) .gt. 0
#else
  call check_all_atom_movement(atm_cnt, crd, gbl_atm_saved_crd, skinnb, ntp, &
                               new_list)
#endif /* MPI */

#endif /* CUDA */

end subroutine scale_system_volume

!*******************************************************************************
!
! Subroutine: mcbar_summary
!
! Description: Print out a summary of the MC barostat statistics
!
!*******************************************************************************

subroutine mcbar_summary()

  use parallel_dat_mod, only : master
  use file_io_dat_mod, only : mdout

  implicit none

  double precision :: success_percent

  if (.not. master) return

  success_percent = &
     dble(total_mcbar_successes) / dble(total_mcbar_attempts) * 100.d0

  write(mdout, '(("| MC Barostat: ",i10," volume changes attempted."))') &
     total_mcbar_attempts
  
  write(mdout, '(("| MC Barostat: ",i10," changes successful (",1f6.2,"%)"))') &
     total_mcbar_successes, success_percent
  write(mdout, '(t2,78("-"),/)')

  return

end subroutine mcbar_summary

end module barostats_mod
