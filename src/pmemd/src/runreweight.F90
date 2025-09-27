#include "copyright.i"
#include "dbg_arrays.i"

!*******************************************************************************
!
! Module: runreweight_mod
!
! Description: Sander 9-compatible MD.
!
!*******************************************************************************

module runreweight_mod

  implicit none

contains

!*******************************************************************************
!
! Subroutine:  runreweight
!
! Description: Driver routine for reweighting trajectories.  Stripped down
!              version of runmd.  Results slightly off due to wrapping /
!              precision.  Also, the reweight usually doesn't have the first
!              frame because the first energy printed is usually based off of
!              the first set of coordinates, which a trajectory does not reprint
!              so step counters are offset by 1 when comparing.
!
!*******************************************************************************

#ifdef MPI
subroutine runreweight(atm_cnt, crd, mass, frc, vel, last_vel, my_atm_lst, &
                 local_remd_method, local_numexchg)
#else
subroutine runreweight(atm_cnt, crd, mass, frc, vel, last_vel)
#endif

  use amd_mod, only : amd_weights_and_energy, write_amd_weights
  use barostats_mod, only : mcbar_trial, mcbar_summary
  use bintraj_mod, only : end_binary_frame, write_binary_cell_dat
  use constante_dat_mod, only : on_cestep
  use constante_mod, only : cnste_explicitmd, cnste_end_step, cnste_write_ceout, &
                 cnste_write_restart, cnste_begin_step, cnste_write_ceout
  use constantph_mod, only : cnstph_begin_step, cnstph_end_step, cnstph_explicitmd, &
                 cnstph_update_pairs, cnstph_write_cpout, cnstph_write_restart
  use constantph_dat_mod, only : on_cpstep
  use constraints_mod, only : atm_igroup, bellyf
  use degcnt_mod, only : degcnt
  use dynamics_dat_mod, only : tmass
  use emap_mod, only: emap_move
  use gamd_mod, only : gamd_weights_and_energy, write_gamd_weights,igamd0
  use gb_force_mod
  use gb_ene_mod
  use gbl_constants_mod, only : KB
  use img_mod, only : gbl_atm_img_map, gbl_img_atm_map
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod, only : skinnb, frameon
  use nmr_calls_mod, only : nmrptx, ndvptx
  use parallel_dat_mod, only : master
  use pbc_mod, only : bcast_pbc, init_pbc, pbc_box, recip, ucell, uc_volume
  use pme_force_mod  ! , only : irespa, pme_force, pme_force_midpoint
  use prmtop_dat_mod, only : atm_qterm, bonda_idx, gbl_bond, ifbox
  use runfiles_mod, only : corpac, mdeng, mdwrit, prntmd, wrapped_corpac, &
                 wrapped_mdwrit
  use scaledMD_mod, only : scaledmd_energy, scaledmd_unscaled_energy, &
                 scaledmd_weight, write_scaledmd_log
  use sgld_mod, only : trxsgld,sgld_avg,sgld_fluc
                 
  use shake_mod, only : num_noshake
  use timers_mod, only : fcve_dist_time, runmd_time, run_end_walltime, &
                 run_start_walltime, update_time, wall, zero_time
  use ti_mod
#ifdef TIDECOMP
  use ti_decomp_mod
#endif

#ifdef DBG_ARRAYS
  use dbg_arrays_mod
#endif
#ifdef EMIL
  use emil_mod
#endif
#ifdef GTI
  use gti_mod
#endif
#ifdef MPI
  use gb_parallel_mod, only : gb_mpi_gathervec
  use inpcrd_dat_mod, only : traj_frames, init_reweight_dat, bcast_inpcrd_dat
  use nebread_mod, only : beadid, nbona, neb_nbead
  use neb_mod
  use nfe_abmd_mod, only : abmd_mdstep => mdstep, abmd_selection_freq => selection_freq
  use nfe_bbmd_mod, only : bbmd_active => active, bbmd_mdstep => mdstep, &
                 bbmd_exchg_freq => exchange_freq
  use nfe_lib_mod, only : nfe_init, nfe_remember_initremd => remember_initremd
  use nfe_setup_mod, only : nfe_on_mdstep => on_mdstep, nfe_on_mdwrit => on_mdwrit
  use parallel_mod, only : mpi_gathervec
  use parallel_processor_mod, only : comm_3dimensions_3ints, &
                 comm_3dimensions_3dbls, mpi_output_transfer
  use pme_recip_midpoint_mod, only : proc_atm_crd, proc_atm_wrap, reorient_flag
  use remd_mod, only : remd_method, remd_types, replica_indexes
#else
  use inpcrd_dat_mod, only : traj_frames, init_reweight_dat
  use nfe_lib_mod, only : nfe_init
  use nfe_setup_mod, only : nfe_on_mdwrit => on_mdwrit
#endif


!---------------------------------------------------------------

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: crd(3, atm_cnt)
  double precision      :: mass(atm_cnt)
  double precision      :: frc(3, atm_cnt)
  double precision      :: vel(3, atm_cnt)
  double precision      :: last_vel(3, atm_cnt)
#ifdef MPI
  integer               :: my_atm_lst(*)
  integer               :: local_remd_method
  integer, intent(inout):: local_numexchg
#else
  integer               :: local_remd_method = 0
  integer               :: remd_types(1) = 0
  integer               :: replica_indexes(1) = 0
#endif

! Local variables:

  double precision      :: aamass
  double precision      :: boltz2
  double precision      :: dtcp
  double precision      :: dtx, dtx_inv, half_dtx
  type(gb_pot_ene_rec)  :: gb_pot_ene
  type(pme_pot_ene_rec) :: pme_pot_ene
  double precision      :: pme_err_est
  double precision      :: eke
  double precision      :: ekmh
  double precision      :: ekpbs
  double precision      :: ekph
  double precision      :: ekin0
  double precision      :: etot_save
  double precision      :: si(si_cnt)
  double precision      :: sit(si_cnt), sit_tmp(si_cnt)
  double precision      :: sit2(si_cnt), sit2_tmp(si_cnt)
#ifdef GBTimer
  integer               :: wall_s, wall_u
  double precision, save  :: strt_time_ms
#endif

#ifdef MPI
  double precision              :: remd_ptot

  double precision              :: proc_new_list_cnt
  double precision              :: new_list_cnt

#endif
  ! Always needed for TI even in the serial code
  double precision, save        :: reduce_buf_in(12)
  double precision, save        :: reduce_buf_out(12)

  double precision      :: fac(3)
  double precision      :: rmu(3)
  double precision      :: factt
  double precision      :: pconv
  double precision      :: rndf, rndfp, rndfs, oldrndfp, oldrndfs
  double precision      :: scaltp
  double precision      :: tmassinv
  double precision      :: ocm(3), xcm(3), vcm(3)       ! for COM velocities
  double precision      :: tspan
  double precision      :: velocity2
  double precision      :: vmax
  double precision      :: wfac
  double precision      :: winf
  double precision      :: sys_x, sys_y, sys_z
  double precision      :: sys_range(3, 2)

#ifdef MPI
  integer               :: atm_lst_idx
#endif
  ! Always needed for TI even in the serial code
  integer               :: buf_ctr

  integer               :: mdloop
  integer               :: i, j, m, k
  integer               :: nrx
  integer               :: nstep
  integer               :: steps_per_avg
  integer               :: nvalid
  integer               :: nvalidi
  integer               :: total_nstep   ! For REMD
  integer               :: total_nstlim  ! For REMD

  logical               :: belly
  logical               :: lout
  logical              :: proc_new_list
  logical              :: do_output_mpi
  logical               :: new_list
  logical               :: onstep
  logical               :: reset_velocities
  logical               :: use_vlimit
  logical               :: vlimit_exceeded
  logical               :: write_restrt
  logical               :: is_langevin       ! Is this a Langevin dynamics simu.
  logical               :: do_mcbar_trial    ! Do an MC barostat trial move?

  logical               :: need_pot_enes     ! refers to potential energy...
  logical               :: need_virials
  logical               :: dump_frcs
  logical               :: timlim_exceeded

#ifdef MPI
  logical               :: collect_crds
  logical               :: collect_vels
  logical               :: all_crds_valid
  logical               :: all_vels_valid
  logical               :: print_exch_data
#endif

  double precision              :: crd_bak(3,atm_cnt)
  double precision, parameter   :: one_third = 1.d0/3.d0
  double precision, parameter   :: pressure_constant = 6.85695d+4

! xyz components of pressure, ekcmt, virial:

  double precision      :: ekcmt(3)
  double precision      :: press(3)
  double precision      :: virial(3)

! Variables used for Langevin dynamics:

  double precision      :: c_ave
  double precision      :: gammai

! Variables and parameters for constant surface tension:
  double precision, parameter :: ten_conv = 100.0d0 !ten_conv - converts
                                                    !dyne/cm to bar angstroms
  double precision      :: pres0x
  double precision      :: pres0y
  double precision      :: pres0z
  double precision      :: gamma_ten_int
  double precision      :: press_tan_ave

!SGLD
  double precision      :: sgsta_rndfp, sgend_rndfp, ignore_solvent

! Variables used for GaMD
!  integer               :: igamd0
  logical               :: update_gamd
  double precision      :: VmaxDt,VminDt,VavgDt,sigmaVDt
  double precision      :: VmaxPt,VminPt,VavgPt,sigmaVPt
  integer,save          :: counts=0

! Reweight
  double precision      :: box_alpha, box_beta, box_gamma
  double precision      :: box(3)
  integer               :: inpcrd_natom
  character(80)         :: inpcrd_title
  double precision      :: input_time ! time read from INPCRD
  integer               :: num_ints, num_reals

! Runmd operates in kcal/mol units for energy, amu for masses,
! and angstoms for distances.  To convert the input time parameters
! from picoseconds to internal units, multiply by 20.455
! (which is 10.0 * sqrt(4.184)).

  call zero_time()

! Here we set total_nstlim. This will just be set to nstlim to start with, but
! will be set to numexchg * nstlim if we're doing REMD. Printing frequency,
! etc. will be based on total_nstep, which is incremented with nstep but not
! reset upon exchanges (as usual)
  if(nstlim .gt. traj_frames) then
    nstlim = traj_frames 
  end if
  total_nstlim = nstlim

! Initialize some variables:
  mdloop = 0
  timlim_exceeded = .false.
  do_mcbar_trial = .false.

  ! We init for case where they are not used but may be referenced in mpi.

  eke = 0.d0
  ekph = 0.d0
  ekpbs = 0.d0
  ekcmt(:) = 0.d0
  press(:) = 0.d0       ! for ntp 0, this prevents dumping of junk.
  virial(:) = 0.d0
  belly = ibelly .gt. 0
  lout = .true.
  use_vlimit = vlimit .gt. 1.0d-7
  ekmh = 0.d0

  tmassinv = 1.d0 / tmass       ! total mass inverse for all atoms.

! If ntwprt .ne. 0, only print the solute atoms in the coordinate/vel. archives.

  nrx  = natom * 3
  if (ntwprt .gt. 0) nrx = ntwprt * 3

! Get the degrees of freedom for the solute (rndfp) and solvent (rndfs).
! Then correct the solute value for the minimum degrees of freedom in the
! system (due to removal of translational or rotational motion, etc.).
  if (usemidpoint) then
#ifdef MPI
    if (.not. allocated(atm_igroup)) allocate(atm_igroup(atm_cnt/numtasks*2))
    if (nbona .eq. 0) then !such as in water
      call degcnt(ibelly, natom, atm_igroup, natom, gbl_bond, gbl_bond(1), &
                  ntc, rndfp, rndfs)
    !gbl_bond and gbl_bond(1) are same
    !we do this since degcnt() is called by others,so keep the format
    else
      call degcnt(ibelly, natom, atm_igroup, natom, gbl_bond, gbl_bond(bonda_idx), &
                  ntc, rndfp, rndfs)
    end if
#endif /*MPI*/
  else
    call degcnt(ibelly, natom, atm_igroup, natom, gbl_bond, gbl_bond(bonda_idx), &
                ntc, rndfp, rndfs)
  end if


  rndfp = rndfp - dble(ndfmin) + dble(num_noshake)

  rndf = rndfp + rndfs  ! total degrees of freedom in system

! Correct the total degrees of freedom for extra points (lone pairs):

  rndf = rndf - 3.d0 * dble(numextra)

! BUGBUG - NOTE that rndfp, rndfs are uncorrected in an extra points context!


! Runmd operates in kcal/mol units for energy, amu for masses,
! and angstoms for distances.  To convert the input time parameters
! from picoseconds to internal units, multiply by 20.455
! (which is 10.0 * sqrt(4.184)).

  boltz2 = KB * 0.5d0

  pconv = 1.6604345d+04
  pconv = pconv * 4.184d0 ! factor to convert pressure from kcal/mole to bar.

  dtx = dt * 20.455d+00
  dtx_inv = 1.0d0 / dtx
  half_dtx = dtx * 0.5d0

  fac(1) = boltz2 * rndf
  fac(2) = boltz2 * rndfp
  if (rndfp .lt. 0.1d0) fac(2) = 1.d-6
  fac(3) = boltz2 * rndfs
  if (rndfs .lt. 0.1d0) fac(3) = 1.d-6

  factt = rndf / (rndf + ndfmin)

  ekin0  = fac(1) * temp0

!I don't think we need these calls to exchange vel. the velocities shouldn't have changed
!at this point.
  if (ti_mode .ne. 0) then
    ! Calculate dof/temperature factors for TI runs
    call ti_degcnt(atm_cnt, num_noshake, boltz2)
    if (ti_latm_cnt(1) .gt. 0) then
      ! Make sure the initial velocities are the same, may be different due to
      ! random initialization
#ifdef CUDA
!sync initial velocities on GPU
#ifdef GTI
      call gti_sync_vector(1, gti_syn_mass+1) ! copy vel over: (linear-combination=1; copy-over =2)
#else
      call gpu_ti_exchange_vel()
#endif  /* GTI  */
#endif  /* CUDA  */
      call ti_exchange_vec(atm_cnt, vel, .false.)
    end if
  end if

! Langevin dynamics setup:

  is_langevin = (gamma_ln .gt. 0.d0)

  gammai = gamma_ln / 20.455d0
  c_ave    = 1.d0 + gammai * half_dtx

  if (is_langevin .and. ifbox .eq. 0) &
    call get_position(atm_cnt, crd, sys_x, sys_y, sys_z, sys_range)

  if (ntp .gt. 0 .and. barostat .eq. 1) dtcp = comp * 1.0d-06 * dt / taup

#ifdef CUDA
  if (ntp .gt. 0) call gpu_ntp_setup()
#endif

! Constant surface tension setup:

  if (csurften > 0) then

    ! Set pres0 in direction of surface tension.
    ! The reference pressure is held constant in on direction dependent
    ! on what the surface tension direction is set to.
    if (csurften .eq. 1) then           ! pres0 in the x direction
      pres0x = pres0

    else if (csurften .eq. 2) then      ! pres0 in the y direction
      pres0y = pres0

    !else if (csurften .eq. 3) then      ! pres0 in the z direction
    else
      pres0z = pres0

    end if

    ! Multiply surface tension by the number of interfaces
    gamma_ten_int = dble(ninterface) * gamma_ten

  end if

! GaMD Initialization
  igamd0 = igamd
  update_gamd = .false.
  VminDt = 1.0d99
  VmaxDt = -1.0d99
  VavgDt = 0.0d0
  sigmaVDt = 0.d0
  VminPt = 1.0d99
  VmaxPt = -1.0d99
  VavgPt = 0.0d0
  sigmaVPt = 0.d0

  nvalid = 0
  nvalidi = 0
  nstep = 0
  total_nstep = 0
  steps_per_avg = ene_avg_sampling

  si(:) = 0.d0
  sit(:) = 0.d0
  sit2(:) = 0.d0

  vel(:,:) = 0.d0

! The _tmp variables below are only used if ntave .gt. 0.  They are used for
! scratch space in the ntave calcs, and to hold the last sit* values between
! calls to the ntave code.

  sit_tmp(:) = 0.d0
  sit2_tmp(:) = 0.d0

#if defined(MPI) && defined(CUDA)
  if (ineb .gt. 0 .and. (beadid .eq. 1 .or. beadid .eq. neb_nbead)) then
    gamma_ln=0.d0
    call gpu_upload_gamma_ln(gamma_ln)
  end if
#endif

#ifdef MPI
  if (usemidpoint) then
    reorient_flag = .false.
  end if

! The following flags keep track of coordinate/velocity update status when
! using mpi

  all_crds_valid = .true.
  all_vels_valid = .true.
#endif

DBG_ARRAYS_TIME_STEP(nstep)

! Clean up the velocity if belly run:

  if (belly) call bellyf(atm_cnt, atm_igroup, vel)

! Make a first dynamics step:

  irespa = 1    ! PME respa step counter.
  if (usemidpoint) then
   proc_new_list = .true.
  else
    nfe_init = 0
    ! We must build a new pairlist the first time we run force:
    new_list = .true.
  end if

  need_pot_enes = .true.        ! maybe not, but most simple for step 0
  need_virials = (ntp .gt. 0 .and. barostat .ne. 2) .or. csurften .gt. 0

#ifdef EMIL
  if ( emil_do_calc .gt. 0 ) then
        call emil_init( atm_cnt, nstep, 1.0/(temp0 * 2 * boltz2 ), &
                     mass, crd, frc, vel)
  end if
#endif

#ifdef CUDA
  call gpu_upload_last_vel(last_vel)
#endif

  if (nstlim .eq. 0) return

  ! Print out the initial state information to cpout (step 0)
  if (.not. usemidpoint) then
    if (icnstph .eq. 1 .and. local_remd_method .eq. 0) &
      call cnstph_write_cpout(0, total_nstlim, t, local_remd_method, remd_types, replica_indexes)

    ! Print out the initial state information to cpout (step 0)
    if (icnste .eq. 1 .and. local_remd_method .eq. 0) &
      call cnste_write_ceout(0, total_nstlim, t, local_remd_method, remd_types, replica_indexes)
  end if

    !=======================================================================
    ! MAIN LOOP FOR PERFORMING THE DYNAMICS STEP:
    ! At this point, the coordinates are a half-step "ahead" of the velocities;
    ! the variable EKMH holds the kinetic energy at these "-1/2" velocities,
    ! which are stored in the array last_vel.
    !=======================================================================
 
    if (.not. usemidpoint) then
      nfe_init = 1
    end if
    do
      DBG_ARRAYS_TIME_STEP(nstep)
#ifdef MPI   
      if(master) then
#endif
        nstep = nstep + 1
        total_nstep = nstep
        ! Gets new coordinate frame

        call init_reweight_dat(num_ints, num_reals, inpcrd_natom, &
                           box_alpha, box_beta, box_gamma, &
                           box, input_time, inpcrd_title, nstep, .false.)
        if(ntb .ne. 0) then
           call init_pbc(box(1), box(2), box(3), box_alpha, box_beta, box_gamma, &
                         vdw_cutoff + skinnb)
        end if
#ifdef MPI
      end if

      !call bcast_inpcrd_dat(natom)
      if (ntb .ne. 0) then
          call bcast_pbc
      end if
#endif

      if(ti_mode .ne. 0) call ti_zero_arrays

#ifdef CUDA
      call gpu_upload_crd(crd)
      call gpu_update_box(ucell, recip, uc_volume)
      call gpu_force_new_neighborlist()
#endif 

      ! Update the igamd flag upon new or restarted GaMD simulation
      if ((igamd.gt.0).or.(igamd0.gt.0)) then
        if (ntcmd.gt.0) then
          if (total_nstep.eq.0) then
            igamd = 0
            if (master) write(mdout,'(/,a,i10,/)') &
              ' | GaMD: Run initial conventional MD with no boost; igamd = ', igamd
          else if ((total_nstep).eq.ntcmd) then
            igamd = igamd0
            if (master) write(mdout,'(/,a,i10,/)') &
              ' | GaMD: Apply boost potential after finishing conventional MD; igamd = ', igamd
          end if
        end if
        update_gamd = .false.
      end if

#ifdef TIDECOMP
      call zero_decomp_energies
#endif

      ! Calculate the force. This also does ekcmt if a regular pme run:

      ! Full energies are only calculated every nrespa steps, and may actually be
      ! calculated less frequently, depending on needs.

      onstep = mod(irespa, nrespa) .eq. 0
      do_mcbar_trial = ntp .ne. 0 .and. barostat .eq. 2 .and. &
        mod(total_nstep+1, mcbarint) .eq. 0

      if (.not. usemidpoint) then
        ! Deciding if in the current step we are gonna perform a constant pH titration
        ! attempt and/or a constant Redox potential titration attempt

        if (icnstph .gt. 0) then
#ifdef MPI
          on_cpstep = mod(irespa+(mdloop-1)*nstlim, ntcnstph) .eq. 0
#else
          on_cpstep = mod(irespa, ntcnstph) .eq. 0
#endif
        end if
        if (icnste .gt. 0) then
#ifdef MPI
          on_cestep = mod(irespa+(mdloop-1)*nstlim, ntcnste) .eq. 0
#else
          on_cestep = mod(irespa, ntcnste) .eq. 0
#endif
        end if

        ! If we are running constant pH simulation, we need to determine if we attempt
        ! a constant pH protonation change and if we need to update the cph pairlist

        if (icnstph .gt. 0) then
          ! Only update the cnstph pairlist with the main pairlist in PME to make
          ! sure that all of the coordinates are valid on every node (this
          ! prevents de-synchronization which destroys the simulation with extreme
          ! prejudice)
          if ( (using_pme_potential .and. new_list) .or. &
#ifdef MPI
          (using_gb_potential .and. mod(irespa+(mdloop-1)*nstlim, nsnb) .eq. 0) ) &
#else
          (using_gb_potential .and. mod(irespa, nsnb) .eq. 0) ) &
#endif
          call cnstph_update_pairs(crd)

          if (on_cpstep) then
            if (icnstph .eq. 1) then
              call cnstph_begin_step
            else
#ifdef MPI
              call cnstph_explicitmd(atm_cnt, crd, mass, frc, vel, last_vel, &
                total_nstep+1, total_nstlim, t, my_atm_lst)
#else
              call cnstph_explicitmd(atm_cnt, crd, mass, frc, vel, last_vel, &
                total_nstep+1, total_nstlim, t)
#endif
            end if
          end if

        end if

        ! If constant pH and constant Redox potential titration attempts are performed in the
        ! same step for implicit solvent, then we need to compute gb_pot_ene here and finish
        ! constant pH

        if (on_cestep .and. icnste .eq. 1 .and. on_cpstep .and. icnstph .eq. 1) then
          need_pot_enes = .true.
          if (using_gb_potential) then
#ifdef MPI
            if (local_remd_method .eq. 0 .or. nstep .gt. 0) &
#endif
            call gb_force(atm_cnt, crd, frc, gb_pot_ene, irespa, need_pot_enes)
          end if

          call cnstph_end_step(gb_pot_ene%dvdl, 0, 0)
          call cnstph_write_cpout(total_nstep+1, total_nstlim, t, local_remd_method, &
                                  remd_types, replica_indexes)
        end if

        ! If we are running constant Redox potential simulation, we need to determine if we attempt
        ! a constant Redox potential protonation change and if we need to update the ce pairlist

        if (icnste .gt. 0) then

          if (on_cestep) then
            if (icnste .eq. 1) then
              call cnste_begin_step
            else
#ifdef MPI
              call cnste_explicitmd(atm_cnt, crd, mass, frc, vel, last_vel, &
                total_nstep+1, total_nstlim, t, my_atm_lst)
#else
              call cnste_explicitmd(atm_cnt, crd, mass, frc, vel, last_vel, &
                total_nstep+1, total_nstlim, t)
#endif
            end if
          end if

        end if
      end if ! not usemidpoint

      if(using_pme_potential) then
#ifdef MPI
        ! If we are doing REMD, there's no need to run pme_force on the first
        ! step after a REMD exchange attempt since we already did it above.

        if (usemidpoint) then
          if (local_remd_method .eq. 0 .or. nstep .gt. 0) &
            call pme_force_midpoint(atm_cnt, crd, frc, &
            my_atm_lst, .true., need_pot_enes, need_virials, &
            pme_pot_ene, nstep, virial, ekcmt, pme_err_est)
          DBG_ARRAYS_DUMP_3DBLE("runmd_frc2", proc_atm_to_full_list,proc_atm_frc,proc_num_atms)
        else ! usemidpoint
          if (local_remd_method .eq. 0 .or. nstep .gt. 0) &
            call pme_force(atm_cnt, crd, frc, gbl_img_atm_map, gbl_atm_img_map, &
            my_atm_lst, .true., need_pot_enes, need_virials, &
            pme_pot_ene, nstep, virial, ekcmt, pme_err_est)
        end if ! usemidpoint
#else
        call pme_force(atm_cnt, crd, frc, gbl_img_atm_map, gbl_atm_img_map, &
          .true., need_pot_enes, need_virials, &
          pme_pot_ene, nstep, virial, ekcmt, pme_err_est)
#ifdef CUDA
        ! Need a sync thread so we download crd this fixes truncated octahedrons.
        call gpu_download_crd(crd_bak)
#endif
#endif /* MPI */
        counts = counts + 1

        ! Store energy terms in state info array for printout.

        si(si_pot_ene) = pme_pot_ene%total
        si(si_vdw_ene) = pme_pot_ene%vdw_tot
        si(si_elect_ene) = pme_pot_ene%elec_tot
        si(si_hbond_ene) = pme_pot_ene%hbond
        si(si_bond_ene) = pme_pot_ene%bond
        si(si_angle_ene) = pme_pot_ene%angle
        si(si_dihedral_ene) = pme_pot_ene%dihedral
        si(si_vdw_14_ene) = pme_pot_ene%vdw_14
        si(si_elect_14_ene) = pme_pot_ene%elec_14
        si(si_restraint_ene) = pme_pot_ene%restraint
        si(si_pme_err_est) = pme_err_est
        si(si_angle_ub_ene) = pme_pot_ene%angle_ub
        si(si_dihedral_imp_ene) = pme_pot_ene%imp
        si(si_cmap_ene) = pme_pot_ene%cmap
        si(si_amd_boost) = pme_pot_ene%amd_boost
        si(si_gamd_boost) = pme_pot_ene%gamd_boost
        si(si_emap_ene) = pme_pot_ene%emap
        si(si_efield_ene) = pme_pot_ene%efield
        ! Total virial will be printed in mden regardless of ntp value.
        si(si_tot_virial) = virial(1) + virial(2) + virial(3)

      else if (using_gb_potential) then
        if (.not. usemidpoint) then
#ifdef MPI
          if (local_remd_method .eq. 0 .or. nstep .gt. 0) &
#endif
          call gb_force(atm_cnt, crd, frc, gb_pot_ene, irespa, need_pot_enes)
        end if

        si(si_pot_ene) = gb_pot_ene%total
        si(si_vdw_ene) = gb_pot_ene%vdw_tot
        si(si_elect_ene) = gb_pot_ene%elec_tot
        si(si_hbond_ene) = gb_pot_ene%gb                ! temporary hack
        si(si_surf_ene) = gb_pot_ene%surf
        si(si_bond_ene) = gb_pot_ene%bond
        si(si_angle_ene) = gb_pot_ene%angle
        si(si_dihedral_ene) = gb_pot_ene%dihedral
        si(si_vdw_14_ene) = gb_pot_ene%vdw_14
        si(si_elect_14_ene) = gb_pot_ene%elec_14
        si(si_restraint_ene) = gb_pot_ene%restraint
        si(si_pme_err_est) = 0.d0
        si(si_angle_ub_ene) = gb_pot_ene%angle_ub
        si(si_dihedral_imp_ene) = gb_pot_ene%imp
        si(si_cmap_ene) = gb_pot_ene%cmap
        si(si_amd_boost) = gb_pot_ene%amd_boost
        si(si_gamd_boost) = gb_pot_ene%gamd_boost
        si(si_emap_ene) = gb_pot_ene%emap
        si(si_efield_ene) = 0.0d0
      end if

      !debug forces
      !#ifdef CUDA
      !      call gpu_download_frc(frc)
      !#endif

      ! Now that we've called force, evaluate the success (or failure) of the
      ! constant pH move and write the result to the output file, if not done yet

      if (.not. usemidpoint) then
        if (on_cpstep .and. icnstph .eq. 1 .and. .not. on_cestep) then
          call cnstph_end_step(gb_pot_ene%dvdl, 0, 0)
          call cnstph_write_cpout(total_nstep+1, total_nstlim, &
            t, local_remd_method, remd_types, replica_indexes)
        end if

        ! Now that we've called force, evaluate the success (or failure) of the
        ! constant Redox potential move and write the result to the output file

        if (on_cestep .and. icnste .eq. 1) then
          call cnste_end_step(gb_pot_ene%dvdl, 0, 0)
          call cnste_write_ceout(total_nstep+1, total_nstlim, &
            t, local_remd_method, remd_types, replica_indexes)
        end if
      end if

      ! PHMD
      if (iphmd /=0 ) then
        call updatephmd(dtx,atm_qterm,total_nstep+1,si(si_kin_ene))
      end if

      ! Figure out if we have to do the force dump here.

      dump_frcs = .false.
      if (ntwf .ne. 0) then
        if (ntwf .lt. 0) then
          if (mod(total_nstep + 1, ntwx) .eq. 0) &
            dump_frcs = .true.
        else
          if (mod(total_nstep + 1, ntwf) .eq. 0) &
            dump_frcs = .true.
        end if
      end if

      ! Call the EMIL absolute free energy calculation, which can modify
      ! forces.
#ifdef EMIL
      if ( emil_do_calc .gt. 0 ) then

        call emil_step(atm_cnt, nstep, 1.0 / (temp0 * 2 * boltz2),&
          mass, crd, frc, vel, &
          gb_pot_ene, pme_pot_ene, &
          ti_ene(1,si_dvdl), si)
      end if
#endif

      ! Reset quantities depending on temp0 and tautp (which may have been
      ! changed by modwt during force call):

      ekin0 = fac(1) * temp0

      if (ti_mode .ne. 0) then
        si(si_dvdl) = ti_ene(1, si_dvdl)
        if (ti_mode .gt. 1 .and. using_gb_potential) si(si_dvdl) = si(si_dvdl) + gb_pot_ene%dvdl
        if (ti_latm_cnt(1) .gt. 0) then
#ifndef CUDA
          !if EMIL is ever supported on GPU, we will need to rework when the forces are
          !syncd on the gpu side because emil will change the forces post pme_force (at
          !least, as currently implemented on the CPU) and we would need to re-upload the
          !forces here

          !we do the force sync in gpu_ti_pme_ene/force
          call ti_exchange_vec(atm_cnt, frc, .true.)
#endif

#if defined(GTI) && defined(CUDA)
          call gti_sync_vector(0, 0) ! copy force over
#endif
        end if
        ti_ekin0(1) = ti_fac(1,1) * temp0
        ti_ekin0(2) = ti_fac(2,1) * temp0
      end if
      ! Pressure coupling:

      if (ntp .gt. 0) then

        si(si_volume) = uc_volume
        if (ti_mode .eq. 0) then
          si(si_density) = tmass / (0.602204d0 * si(si_volume))
        else
          ti_ene(1,si_density) = ti_tmass(1) / (0.602204d0 * si(si_volume))
          ti_ene(2,si_density) = ti_tmass(2) / (0.602204d0 * si(si_volume))
        end if
        si(si_tot_ekcmt) = 0.d0
        si(si_tot_press) = 0.d0

        do m = 1, 3
          ekcmt(m) = ekcmt(m) * 0.5d0
          si(si_tot_ekcmt) = si(si_tot_ekcmt) + ekcmt(m)
          if (need_virials) then
            press(m) = (pconv + pconv) * (ekcmt(m) - virial(m)) / si(si_volume)
            si(si_tot_press) = si(si_tot_press) + press(m)
          end if
        end do

        si(si_tot_press) = si(si_tot_press) / 3.d0

        ! Constant surface tension output:

        if (csurften .gt. 0) then

          if (csurften .eq. 1) then          ! Surface tension in the x direction
            si(si_gamma_ten) = &
              pbc_box(1) * (press(1) - 0.5d0 * (press(2) + press(3))) / ( ninterface * ten_conv )
            !pbc_box(1) / ninterface / ten_conv * (press(1) - 0.5d0 * (press(2) + press(3)))

          else if (csurften .eq. 2) then     ! Surface tension in the y direction
            si(si_gamma_ten) = &
              pbc_box(2) * (press(2) - 0.5d0 * (press(1) + press(3))) / ( ninterface * ten_conv )
            !pbc_box(2) / ninterface / ten_conv * (press(2) - 0.5d0 * (press(1) + press(3)))

          else
            !else if (csurften .eq. 3) then     ! Surface tension in the z direction
            si(si_gamma_ten) = &
              pbc_box(3) * (press(3) - 0.5d0 * (press(1) + press(2))) / ( ninterface * ten_conv )
            ! pbc_box(3) / ninterface / ten_conv * (press(3) - 0.5d0 * (press(1) + press(2)))

          end if

        end if

      end if

#ifdef GBTimer
      call get_wall_time(wall_s, wall_u)
      print *, "reset vel =", dble(wall_s) * 1000.0d0 + &
        dble(wall_u) / 1000.0d0 - strt_time_ms

      call get_wall_time(wall_s, wall_u)
      strt_time_ms = dble(wall_s) * 1000.0d0 + dble(wall_u) / 1000.0d0
#endif

#ifdef CUDA
      ! Now we do the force dump, before the update since we overwrite the forces
      ! with the 'old' coordinates in the coordinate update.
      if (dump_frcs) then
        !for TI, we still haven't sync'd forces
        call gpu_download_frc(frc)
#ifdef MPI
        ! We need to synchronize the forces to the master here
        if (using_pme_potential) then
          call mpi_gathervec(atm_cnt, frc)
        else if (using_gb_potential) then
          call gb_mpi_gathervec(atm_cnt, frc)
        end if
#endif
        call corpac(nrx, frc, 1, mdfrc)
        !why do we do the update here? shouldn't it be post gpu_update?
        !I think we need to call this here because we're not doing the langevin setvel
        !stuff on the GPU but we are changing velocities.
      end if ! dump_frcs

      !we might need to do something with the velocities here
      !like do a vector exchange
      !    if (ti_mode .eq. 0) then

#ifdef GTI
      ! copy vel and crd over before updates (linear-combination=1; copy-over =2)
      if (ti_mode .ne.0)  then
        ! call gti_sync_vector(1, gti_syn_mass+1) ! vel
        ! call gti_sync_vector(2, 2) ! crd
        call gti_sync_vector(3, gti_syn_mass+3) ! crd
      end if

#endif

      !call gpu_update(dt, temp0, gamma_ln)

#endif /* ifdef CUDA */

#ifdef GBTimer
      call get_wall_time(wall_s, wall_u)
      print *, "wfac  or sgld =", dble(wall_s) * 1000.0d0 + &
        dble(wall_u) / 1000.0d0 - strt_time_ms
#endif

      if (.not. usemidpoint) then
        ! Update EMAP rigid domains
        if (iemap > 0) call emap_move()
      end if

      ! Consider vlimit:

#ifdef GBTimer
      call get_wall_time(wall_s, wall_u)
      strt_time_ms = dble(wall_s) * 1000.0d0 + dble(wall_u) / 1000.0d0
#endif

#ifdef GBTimer
      call get_wall_time(wall_s, wall_u)
      print *, "vlimit exceeded =", dble(wall_s) * 1000.0d0 + &
            dble(wall_u) / 1000.0d0 - strt_time_ms

      call get_wall_time(wall_s, wall_u)
      strt_time_ms = dble(wall_s) * 1000.0d0 + dble(wall_u) / 1000.0d0
#endif

      ! Update the positions, putting the "old" positions into frc():

#ifndef CUDA

      ! Now we do the force dump (we already did it on the GPU, so do it only on
      ! the CPU here)

      if (dump_frcs) then
#ifdef MPI
        ! We need to synchronize the forces to the master here
        if (using_pme_potential) then
          call mpi_gathervec(atm_cnt, frc)
        else if (.not. usemidpoint .and. using_gb_potential) then
          call gb_mpi_gathervec(atm_cnt, frc)
        end if
#endif
        if (master) call corpac(nrx, frc, 1, mdfrc)
      end if

#endif /* ifndef CUDA */
#ifdef GBTimer
      call get_wall_time(wall_s, wall_u)
      print *, "dump frcs =", dble(wall_s) * 1000.0d0 + &
         dble(wall_u) / 1000.0d0 - strt_time_ms

      call get_wall_time(wall_s, wall_u)
      strt_time_ms = dble(wall_s) * 1000.0d0 + dble(wall_u) / 1000.0d0
#endif

#ifdef MPI
      if (usemidpoint) then
        call update_time(runmd_time)
        ! Now distribute the coordinates:
        call comm_3dimensions_3dbls(proc_atm_crd, .true.)
        call comm_3dimensions_3ints(proc_atm_wrap, .true.)

        DBG_ARRAYS_DUMP_CRD("runmd_crd1", proc_atm_to_full_list, proc_atm_crd,proc_num_atms, proc_atm_wrap)

        call update_time(fcve_dist_time)
      end if ! usemidpoint

  ! Update the step counter and the integration time:
        nstep = nstep + 1
        total_nstep = total_nstep + 1
        t = t + dt


    ! Here we update running logging averages as appropriate:

        if (loadbal_verbose .gt. 1) then

          log_owned_img_cnt_avg = log_owned_img_cnt_avg + &
                                  (log_owned_img_cnt * log_recip_nstep)

          log_used_img_cnt_avg = log_used_img_cnt_avg + &
                                  (log_used_img_cnt * log_recip_nstep)

          log_owned_atm_cnt_avg = log_owned_atm_cnt_avg + &
                                  (log_owned_atm_cnt * log_recip_nstep)

          log_used_atm_cnt_avg = log_used_atm_cnt_avg + &
                                  (log_used_atm_cnt * log_recip_nstep)

          log_used_atm_source_cnt_avg = log_used_atm_source_cnt_avg + &
                                  (log_used_atm_source_cnt * log_recip_nstep)

          log_provided_atm_cnt_avg = log_provided_atm_cnt_avg + &
                                  (log_provided_atm_cnt * log_recip_nstep)

          log_provided_atm_sink_cnt_avg = log_provided_atm_sink_cnt_avg + &
                                  (log_provided_atm_sink_cnt * log_recip_nstep)

        end if
#endif

  ! nvalid is the number of steps where all energies are calculated; the onstep
  ! variable indicates that the kinetic energy is available (when using respa).

#if defined(CUDA) && defined (GTI)
        if ((total_nstep .gt. gti_eq_nstep) .and. (mod (total_nstep, steps_per_avg) .eq. 0)) then
#else
        if (mod(total_nstep, steps_per_avg) .eq. 0) then
#endif

          nvalid = nvalid + 1

          if (master) then
            do m = 1, si_cnt
              sit(m) = sit(m) + si(m)
              sit2(m) = sit2(m) + si(m) * si(m)
            end do
            if (ti_mode .ne. 0) call ti_update_avg_ene
          end if

        end if

        lout = mod(total_nstep, ntpr) .eq. 0 .and. onstep

        ! added for rbornstat
        if (mod(irespa,nrespai) == 0 .or. irespa < 2) nvalidi = nvalidi + 1

        irespa = irespa + 1

        ! Output from this step if required:

        ! Only the master needs to do the output:

        if (usemidpoint) then
#ifdef MPI
          ! Restrt:

          do_output_mpi = .false.

          write_restrt = .false.

          if (total_nstep .eq. total_nstlim) then
            write_restrt = .true.
          else
            if (mod(total_nstep, ntwr) .eq. 0 .and. onstep) write_restrt = .true.
          end if

          if (write_restrt) then
            call mpi_output_transfer(crd, vel)
            do_output_mpi = .true.
            if (master) then

              if (iwrap .eq. 0) then
                call mdwrit(total_nstep, atm_cnt, crd, vel, t)
              else
                call wrapped_mdwrit(total_nstep, atm_cnt, crd, vel, t)
              end if

              !          if (icnstph .gt. 0) call cnstph_write_restart
            end if !master
          end if

          ! Coordinate archive:

          if (ntwx .gt. 0) then

            if (mod(total_nstep, ntwx) .eq. 0) then
              if (.not. do_output_mpi) then
                call mpi_output_transfer(crd, vel)
                do_output_mpi = .true.
              end if

              if (local_remd_method .eq. 1 .and. master .and. ioutfm .eq. 0) then
                !            if (trxsgld)then
                !            write (mdcrd, '(a,3(1x,i8),1x,i8)') 'RXSGLD ', repnum, &
                !              mdloop, total_nstep, replica_indexes(1)
                !            else
                write (mdcrd, '(a,3(1x,i8),1x,f8.3)') 'REMD ', repnum, &
                  mdloop, total_nstep, temp0
                !            end if
              end if
              if (master) then
                if (iwrap .eq. 0) then
                  call corpac(nrx, crd, 1, mdcrd)
                else
                  call wrapped_corpac(atm_cnt, nrx, crd, 1, mdcrd)
                end if

                if (ntb .gt. 0) then
                  if (ioutfm .eq. 1) then
                    call write_binary_cell_dat
                  else
                    call corpac(3, pbc_box, 1, mdcrd)
                  end if
                end if

                ! For ioutfm == 1, coordinate archive may also contain vels.

                if (ioutfm .eq. 1 .and. ntwv .lt. 0) then
                  call corpac(nrx, vel, 1, mdvel)

                end if
              end if !master

              ! For ioutfm == 1, coordinate archive may also contain vels

            end if

          end if

          ! Velocity archive:

          if (ntwv .gt. 0) then

            if (mod(total_nstep, ntwv) .eq. 0) then
              if (.not. do_output_mpi) then
                call mpi_output_transfer(crd, vel)
                do_output_mpi = .true.
              end if
              if (master) then

                call corpac(nrx, vel, 1, mdvel)
              end if
            end if

          end if

          ! Flush binary netCDF file(s) if necessary...

          if (ioutfm .eq. 1) then

            if (ntwx .gt. 0) then
              if (mod(total_nstep, ntwx) .eq. 0) then
                if (.not. do_output_mpi) then
                  call mpi_output_transfer(crd, vel)
                  do_output_mpi = .true.
                end if
                if (master) then
                  call end_binary_frame(mdcrd)
                end if
              end if
            end if

            if (ntwv .gt. 0) then

              if (mod(total_nstep, ntwv) .eq. 0) then
                if (.not. do_output_mpi) then
                  call mpi_output_transfer(crd, vel)
                  do_output_mpi = .true.
                end if
                if (master) then

                  call end_binary_frame(mdvel)
                end if
              end if
            end if

            if (ntwf .gt. 0) then
              if (mod(total_nstep, ntwf) .eq. 0) then
                if (master) then
                  call end_binary_frame(mdfrc)
                end if
              end if
            end if

          end if

#endif /*MPI*/
        else ! usemidpoint

          if (master) then

            ! Restrt:
            write_restrt = .false.

            if (total_nstep .eq. total_nstlim) then
              write_restrt = .true.
            else
              if (mod(total_nstep, ntwr) .eq. 0 .and. onstep) write_restrt = .true.
            end if

            if (write_restrt) then
#ifdef CUDA
              call gpu_download_crd(crd)
              call gpu_download_vel(vel)
#endif
#ifdef MPI
#else
              if (numextra .gt. 0 .and. frameon .ne. 0) &
                call zero_extra_pnts_vec(vel, ep_frames, gbl_frame_cnt)
#endif

              if (iwrap .eq. 0) then
                call mdwrit(total_nstep, atm_cnt, crd, vel, t)
              else
                call wrapped_mdwrit(total_nstep, atm_cnt, crd, vel, t)
              end if

              if (icnstph .gt. 0) call cnstph_write_restart(total_nstep)

              if (icnste .gt. 0) call cnste_write_restart(total_nstep)

              !PHMD
              if (iphmd .gt. 0) call phmdwriterestart

              ! NFE
              if (infe.gt.0) call nfe_on_mdwrit()

            end if

            ! Coordinate archive:

            if (ntwx .gt. 0) then

              if (mod(total_nstep, ntwx) .eq. 0) then
#ifdef CUDA
                call gpu_download_crd(crd)
#endif

#ifdef MPI
                if (local_remd_method .eq. 1 .and. master .and. ioutfm .eq. 0) then
                  if (trxsgld)then
                    write (mdcrd, '(a,3(1x,i8),1x,i8)') 'RXSGLD ', repnum, &
                      mdloop, total_nstep, replica_indexes(1)
                  else
                    write (mdcrd, '(a,3(1x,i8),1x,f8.3)') 'REMD ', repnum, &
                      mdloop, total_nstep, temp0
                  end if
                end if
#endif
                if (iwrap .eq. 0) then
                  call corpac(nrx, crd, 1, mdcrd)
                else
                  call wrapped_corpac(atm_cnt, nrx, crd, 1, mdcrd)
                end if

                if (ntb .gt. 0) then
                  if (ioutfm .eq. 1) then
                    call write_binary_cell_dat
                  else
                    call corpac(3, pbc_box, 1, mdcrd)
                  end if
                end if

                ! For ioutfm == 1, coordinate archive may also contain vels.

                if (ioutfm .eq. 1 .and. ntwv .lt. 0) then
#ifdef CUDA
                  call gpu_download_vel(vel)
#endif
#ifdef MPI
#else
                  ! The zeroing of extra points may be redundant here.
                  ! Not a big deal overall though.
                  if (numextra .gt. 0 .and. frameon .ne. 0) &
                    call zero_extra_pnts_vec(vel, ep_frames, gbl_frame_cnt)
#endif /* MPI */
                  call corpac(nrx, vel, 1, mdvel)

                end if

                ! For ioutfm == 1, coordinate archive may also contain vels

                if (iamd.gt.0)then

                  !Flush amdlog file containing the list of weights and energies at each step
#ifdef CUDA
                  if (master) call gpu_download_amd_weights(amd_weights_and_energy)
#endif
                  if (master) call write_amd_weights(ntwx,total_nstep)

                end if

                if (igamd.gt.0)then

                  !Flush gamdlog file containing the list of weights and energies at each step
#ifdef CUDA
                  if (master) call gpu_download_gamd_weights(gamd_weights_and_energy)
#endif
                  if (master) call write_gamd_weights(ntwx,total_nstep)

                end if

                if (scaledMD.gt.0)then

                  !Flush scaledMDlog file containing the list of lambda and total energy at each step
#ifdef CUDA
                  if (master) call gpu_download_scaledmd_weights(scaledMD_energy,scaledMD_weight,scaledMD_unscaled_energy)
#endif
                  if (master) call write_scaledMD_log(ntwx,total_nstep)

                end if

              end if

            end if

            ! Velocity archive:

            if (ntwv .gt. 0) then

              if (mod(total_nstep, ntwv) .eq. 0) then
#ifdef CUDA
                call gpu_download_vel(vel)
#endif
#ifdef MPI
#else
               ! The zeroing of extra points may be redundant here.
               ! Not a big deal overall though.
               if (numextra .gt. 0 .and. frameon .ne. 0) &
                 call zero_extra_pnts_vec(vel, ep_frames, gbl_frame_cnt)
#endif /* MPI */

                call corpac(nrx, vel, 1, mdvel)

              end if

            end if

          end if ! Exit master for velocity/crd printing

        end if ! usemidpoint

        if (master) then

          ! Energy archive:

          if (ntwe .gt. 0) then
            if (mod(total_nstep, ntwe) .eq. 0 .and. onstep) then
              call mdeng(total_nstep, t, si, fac, press, virial, ekcmt)
            end if
          end if

          if (lout) then
            if (ti_mode .eq. 0) then
              call prntmd(total_nstep, total_nstlim, t, si, fac, 7, .false., &
                mdloop)
            else
              if (ifmbar .ne. 0 .and. do_mbar) &
                call ti_print_mbar_ene(si(si_pot_ene))
              do ti_region = 1, 2

                if ( emil_sc_lcl .ne. 1 ) write(mdout, 601) ti_region

                si(si_tot_ene) = ti_ene_aug(ti_region,ti_tot_tot_ene)
                si(si_kin_ene) = ti_ene_aug(ti_region,ti_tot_kin_ene)
                si(si_density) = ti_ene(ti_region,si_density)
                ti_temp = ti_ene_aug(ti_region,ti_tot_temp)

                if ( (emil_sc_lcl .ne. 1) .or. (ti_region .eq. 1) ) then
                  call prntmd(total_nstep, total_nstlim, t, si, fac, 7, .false., &
                    mdloop)
                  call ti_print_ene(ti_region, ti_ene(ti_region,:), &
                    ti_ene_aug(ti_region,:))
                end if
              end do
            end if
            if (.not. usemidpoint) then
              if (nmropt .ne. 0) call nmrptx(mdout)
            end if
          end if

          if (.not. usemidpoint) then
            ! Flush binary netCDF file(s) if necessary...

            if (ioutfm .eq. 1) then

              if (ntwx .gt. 0) then
                if (mod(total_nstep, ntwx) .eq. 0) then
                  call end_binary_frame(mdcrd)
                end if
              end if

              if (ntwv .gt. 0) then

                if (mod(total_nstep, ntwv) .eq. 0) then
#ifdef MPI
#else
                  ! The zeroing of extra points may be redundant here.
                  ! Not a big deal overall though.
                  if (numextra .gt. 0 .and. frameon .ne. 0) &
                    call zero_extra_pnts_vec(vel, ep_frames, gbl_frame_cnt)
#endif /* MPI */

                  call end_binary_frame(mdvel)

                end if
              end if

              if (ntwf .gt. 0) then
                if (mod(total_nstep, ntwf) .eq. 0) then
                  call end_binary_frame(mdfrc)
                end if
              end if

            end if
          end if ! not usemidpoint

          ! Output running average:

          if (ntave .gt. 0) then
            if (mod(total_nstep, ntave) .eq. 0 .and. onstep) then

              write(mdout, 542)

              tspan = ntave / nrespa

              ! Coming into this loop, the _tmp variables hold the values of
              ! sit, sit2 when this routine was last called (or 0.d0).  The _tmp
              ! vars are used here as scatch space and then updated with the
              ! current sit, sit2.

              do m = 1, si_cnt
                sit_tmp(m) = (sit(m) - sit_tmp(m)) / tspan
                sit2_tmp(m) = (sit2(m) - sit2_tmp(m)) / tspan - &
                  sit_tmp(m) * sit_tmp(m)
                if (sit2_tmp(m) .lt. 0.d0) sit2_tmp(m) = 0.d0
                sit2_tmp(m) = sqrt(sit2_tmp(m))
              end do
              if (ti_mode .eq. 0) then
                if (.not. usemidpoint) then
                  if (isgld>0)call sgld_avg(.false.)
                end if
                write(mdout, 540) ntave / nrespa
                call prntmd(total_nstep, total_nstlim, t, sit_tmp, fac, 0, &
                  .false., mdloop)
                if (.not. usemidpoint) then
                  if (isgld>0)call sgld_fluc(.false.)
                end if ! not usemidpoint
                write(mdout, 550)
                call prntmd(total_nstep, total_nstlim, t, sit2_tmp, fac, 0, &
                  .true., mdloop)
              else
                !also potentially bunk in GPU version because ti_update_ene_avg might be bunk
                call ti_calc_avg_ene(tspan, .false.)
                do ti_region = 1, 2
                  if ( emil_sc_lcl .ne. 1 ) write(mdout, 601) ti_region
                  if ( (emil_sc_lcl .ne. 1) .or. (ti_region .eq. 1) ) then
                    write(mdout, 540) ntave / nrespa
                  end if

                  sit_tmp(si_tot_ene) = ti_ene_aug_tmp(ti_region,ti_tot_tot_ene)
                  sit_tmp(si_kin_ene) = ti_ene_aug_tmp(ti_region,ti_tot_kin_ene)
                  sit_tmp(si_density) = ti_ene_tmp(ti_region,si_density)
                  ti_temp = ti_ene_aug_tmp(ti_region,ti_tot_temp)

                  if ( (emil_sc_lcl .ne. 1) .or. (ti_region .eq. 1) ) then
                    call prntmd(total_nstep, total_nstlim, t, sit_tmp, fac, 0, &
                      .false., mdloop)
                    call ti_print_ene(ti_region, ti_ene_tmp(ti_region,:), &
                      ti_ene_aug_tmp(ti_region,:))
                  end if

                  ! RMSD
                  if ( (emil_sc_lcl .ne. 1) .or. (ti_region .eq. 1) ) write(mdout, 550)

                  sit2_tmp(si_tot_ene) = ti_ene_aug_tmp2(ti_region,ti_tot_tot_ene)
                  sit2_tmp(si_kin_ene) = ti_ene_aug_tmp2(ti_region,ti_tot_kin_ene)
                  sit2_tmp(si_density) = ti_ene_tmp2(ti_region,si_density)
                  ti_temp = ti_ene_aug_tmp2(ti_region,ti_tot_temp)

                  if ( (emil_sc_lcl .ne. 1) .or. (ti_region .eq. 1) ) then
                    call prntmd(total_nstep, total_nstlim, t, sit2_tmp, fac, 0, &
                      .true., mdloop)
                    call ti_print_ene(ti_region, ti_ene_tmp2(ti_region,:), &
                      ti_ene_aug_tmp2(ti_region,:))

                  end if

                end do

                ! Following sander, we only output the average DV/DL contributions
                write(mdout, 541) ntave / nrespa
                call prntmd(total_nstep, total_nstlim, t,ti_ene_delta_tmp, fac, &
                  0, .true., mdloop)
              end if

              write(mdout, 542)

              ! GaMD: save the energy average and standard deviation
              if ( ((igamd.gt.0).or.(igamd0.gt.0)) .and. &
                ( (total_nstep.ge.ntcmdprep .and. total_nstep.le.ntcmd) .or. &
                (total_nstep.ge.(ntcmd+ntebprep) .and. total_nstep.le.(ntcmd+nteb)) ) ) then
                VavgDt = sit_tmp(si_dihedral_ene)
                sigmaVDt = sit2_tmp(si_dihedral_ene)
                if ((igamd.eq.4 .or. igamd.eq.5) .or. (igamd0.eq.4 .or. igamd0.eq.5)) then
                  VavgPt = sit_tmp(si_elect_ene) + sit_tmp(si_vdw_ene) + &
                    sit_tmp(si_vdw_14_ene) + sit_tmp(si_elect_14_ene) + sit_tmp(si_efield_ene)
                  sigmaVPt = sqrt(sit2_tmp(si_elect_ene)**2 + sit2_tmp(si_vdw_ene)**2 + &
                    sit2_tmp(si_vdw_14_ene)**2 + sit2_tmp(si_elect_14_ene)**2 + sit2_tmp(si_efield_ene)**2)
                else
                  VavgPt = sit_tmp(si_pot_ene)
                  sigmaVPt = sit2_tmp(si_pot_ene)
                end if
              end if

              sit_tmp(:) = sit(:)
              sit2_tmp(:) = sit2(:)
            end if
          end if

          if ((igamd.gt.0).or.(igamd0.gt.0)) then
            ! run cMD to calculate Vmax, Vmin, Vavg, sigmaV and then GaMD parameters
            if ((total_nstep.le.ntcmd).and.onstep) then
              if (total_nstep.gt.ntcmdprep) then
                if ((igamd.eq.4 .or. igamd.eq.5) .or. (igamd0.eq.4 .or. igamd0.eq.5)) then
                  VmaxPt = max(VmaxPt,si(si_elect_ene)+si(si_vdw_ene)+si(si_vdw_14_ene)+si(si_elect_14_ene)+si(si_efield_ene))
                  VminPt = min(VminPt,si(si_elect_ene)+si(si_vdw_ene)+si(si_vdw_14_ene)+si(si_elect_14_ene)+si(si_efield_ene))
                else
                  VmaxPt = max(VmaxPt,si(si_pot_ene))
                  VminPt = min(VminPt,si(si_pot_ene))
                end if

                VmaxDt = max(VmaxDt,si(si_dihedral_ene))
                VminDt = min(VminDt,si(si_dihedral_ene))

                if (mod(total_nstep, ntave).eq.0 .and. onstep) then
                  VmaxP = max(VmaxP,VmaxPt)
                  VminP = min(VminP,VminPt)
                  VavgP = VavgPt
                  sigmaVP = sigmaVPt

                  VmaxD = max(VmaxD,VmaxDt)
                  VminD = min(VminD,VminDt)
                  VavgD = VavgDt
                  sigmaVD = sigmaVDt

                  write(mdout,'(a,i10,4f14.4)') ' Energy statistics: step,VmaxP,VminP,VavgP,sigmaVP = ', &
                    total_nstep,VmaxP,VminP,VavgP,sigmaVP
                  write(mdout,'(a,i10,4f14.4)') ' Energy statistics: step,VmaxD,VminD,VavgD,sigmaVD = ', &
                    total_nstep,VmaxD,VminD,VavgD,sigmaVD
                end if
              end if
            end if

            ! update Vmax, Vmin, Vavg, sigmaV and then GaMD parameters after adding boost potential
            if ((total_nstep.le.(ntcmd+nteb)).and.(total_nstep.ge.ntcmd).and.(mod(total_nstep,nrespa).eq.0).and.onstep) then
              if ( (igamd.eq.1 .or. igamd.eq.3 .or. igamd.eq.4 .or. igamd.eq.5) ) then
                if (total_nstep.ge.(ntcmd+ntebprep)) then
                  if (igamd.eq.4 .or. igamd.eq.5) then
                    VmaxPt = max(VmaxPt,si(si_elect_ene)+si(si_vdw_ene)+si(si_vdw_14_ene)+si(si_elect_14_ene)+si(si_efield_ene))
                    VminPt = min(VminPt,si(si_elect_ene)+si(si_vdw_ene)+si(si_vdw_14_ene)+si(si_elect_14_ene)+si(si_efield_ene))
                  else
                    VmaxPt = max(VmaxPt,si(si_pot_ene))
                    VminPt = min(VminPt,si(si_pot_ene))
                  end if
                end if
                ! update GaMD parameters
                if (total_nstep.eq.ntcmd) then
                  update_gamd = .true.
                else
                  call gamd_calc_Vstat(VmaxPt,VminPt,VavgPt,sigmaVPt,cutoffP,iVm,ntave,total_nstep, &
                    update_gamd,VmaxP,VminP,VavgP,sigmaVP)
                end if
                if ( update_gamd ) then
                  call gamd_calc_E_k0(iE,sigma0P,VmaxP,VminP,VavgP,sigmaVP,EthreshP,k0P)
                  kP = k0P/(VmaxP - VminP)
                  write(mdout,'(a,i10,7f14.4)') '| GaMD updated parameters: step,VmaxP,VminP,VavgP,sigmaVP,k0P,kP,EthreshP = '
                  write(mdout,'(a,i10,7f14.4)') '| ', total_nstep,VmaxP,VminP,VavgP,sigmaVP,k0P,kP,EthreshP
                end if
              end if
              if ( igamd.eq.2 .or. igamd.eq.3 .or. igamd.eq.5) then
                if (total_nstep.gt.(ntcmd+ntebprep)) then
                  VmaxDt = max(VmaxDt,si(si_dihedral_ene))
                  VminDt = min(VminDt,si(si_dihedral_ene))
                end if
                ! update GaMD parameters
                if (total_nstep.eq.ntcmd) then
                  update_gamd = .true.
                else
                  call gamd_calc_Vstat(VmaxDt,VminDt,VavgDt,sigmaVDt,cutoffD,iVm,ntave,total_nstep, &
                    update_gamd,VmaxD,VminD,VavgD,sigmaVD)
                end if
                if ( update_gamd ) then
                  call gamd_calc_E_k0(iE,sigma0D,VmaxD,VminD,VavgD,sigmaVD,EthreshD,k0D)
                  kD = k0D/(VmaxD - VminD)
                  write(mdout,'(a,i10,7f14.4)') '| GaMD updated parameters: step,VmaxD,VminD,VavgD,sigmaVD,k0D,kD,EthreshD = '
                  write(mdout,'(a,i10,7f14.4)') '| ', total_nstep,VmaxD,VminD,VavgD,sigmaVD,k0D,kD,EthreshD
                end if
              end if
            end if

            ! save gamd restart file
            if ( (mod(total_nstep, ntwr).eq.0) .and. (total_nstep.le.(ntcmd+nteb)) .and. onstep) then
              open(unit=gamdres,file=gamdres_name,action='write')
              if (igamd.eq.1 .or. igamd.eq.3 .or. igamd.eq.4 .or. igamd.eq.5) then
                write(gamdres,*) VmaxP, VminP, VavgP, sigmaVP
                write(mdout,'(a,4f16.4)')'| Saved GaMD restart parameters: VmaxP,VminP,VavgP,sigmaVP = '
                write(mdout,'(a,4f16.4)')'| ', VmaxP,VminP,VavgP,sigmaVP
              end if
              if (igamd.eq.2 .or. igamd.eq.3 .or. igamd.eq.5) then
                write(gamdres,*) VmaxD, VminD, VavgD, sigmaVD
                write(mdout,'(a,4f16.4)')'| Saved GaMD restart parameters: VmaxD,VminD,VavgD,sigmaVD = '
                write(mdout,'(a,4f16.4)')'| ', VmaxD,VminD,VavgD,sigmaVD
              end if
              close(gamdres)
              write(mdout,*) ''
            end if
          end if ! end GaMD

        end if ! end master's output

        if ( igamd.gt.0 ) then
#ifdef MPI
          call mpi_bcast(update_gamd, 1, mpi_logical, 0, &
            pmemd_comm, err_code_mpi)
          call mpi_bcast(EthreshP, 1, mpi_double_precision, 0, &
            pmemd_comm, err_code_mpi)
          call mpi_bcast(kP, 1, mpi_double_precision, 0, &
            pmemd_comm, err_code_mpi)
          call mpi_bcast(EthreshD, 1, mpi_double_precision, 0, &
            pmemd_comm, err_code_mpi)
          call mpi_bcast(kD, 1, mpi_double_precision, 0, &
            pmemd_comm, err_code_mpi)
#endif /* MPI */
#ifdef CUDA
          if ( update_gamd ) call gpu_gamd_update(EthreshP, kP, EthreshD, kD)
#endif
        end if ! igamd.gt.0

        ! TI: Update dynamic lambda
        if (ti_mode .ne. 0) then
          if (ntave .gt. 0 .and. dynlmb .gt. 0.d0) then
            if (mod(nstep, ntave) .eq. 0 .and. onstep) then
              call ti_dynamic_lambda
            end if
          end if
        end if
#ifdef MPI
        if (infe.gt.0) then
#ifdef CUDA
          if (abmd_selection_freq.ne.0) then
            call gpu_download_crd(crd)
            if (abmd_mdstep.gt.1.and.mod(abmd_mdstep,abmd_selection_freq).eq.0) &
              call gpu_download_vel(vel)
          end if

          if (bbmd_active.eq.1.and.mod(bbmd_mdstep+1,bbmd_exchg_freq).eq.0) then
            call gpu_download_crd(crd)
            call gpu_download_vel(vel)
          end if
#endif
          if (.not. usemidpoint) then
            call nfe_remember_initremd(remd_method.gt.0)
            call nfe_on_mdstep(pme_pot_ene%total, crd, vel)
          end if ! not usemidpoint
#ifdef CUDA
          if (abmd_selection_freq.ne.0) then
            if (abmd_mdstep.gt.1.and.mod(abmd_mdstep,abmd_selection_freq).eq.0) then
              call gpu_upload_crd(crd)
              call gpu_upload_vel(vel)
            end if
          end if

          if (bbmd_active.eq.1.and.mod(bbmd_mdstep,bbmd_exchg_freq).eq.0) then
            call gpu_upload_crd(crd)
            call gpu_upload_vel(vel)
          end if
#endif
        end if
#endif /* MPI */
        !--------------------------------------------------------------------
        !--------------------------------------------------------------------

        if (nstep .ge. nstlim) then
#ifdef TIDECOMP
            if (ntwd .gt. 0) then
                call mddecomp(total_nstep)
            end if
            exit
#endif
        end if
        ! Stop the run if timlim is set and reached. One more step will be taken
        if (nstep .ge. nstlim) exit

        ! Stop the run if timlim is set and reached. One more step will be taken
        ! to ensure a final restart is written.
        if (timlim .gt. 0) then
#       ifdef MPI
          ! In parallel this check is only done on the master node and then
          ! broadcast everywhere else if necessary.
          if (worldrank .eq. 0) then
#       endif
            call wall(run_end_walltime)
            if ((run_end_walltime - run_start_walltime) .gt. timlim) &
              timlim_exceeded = .true.
#       ifdef MPI
          end if
          call mpi_bcast(timlim_exceeded,1,mpi_logical,0,mpi_comm_world,err_code_mpi)
#       endif
          if (timlim_exceeded) then
            ! Reset limits so next step is the final one.
            nstlim = nstep + 1
            total_nstlim = total_nstep + 1
            if (master) &
              write(mdout,'(2(a,i8))') '| Wall clock limit reached (timlim = ', &
              timlim, ' s), step ', total_nstep
          end if
        end if  ! end timlim check

        call update_time(runmd_time)

      end do ! Major cycle back to new step unless we have reached our limit:

    if (master) then

    tspan = nvalid

    if (nvalid .gt. 0) then

      do m = 1, si_cnt
        sit(m) = sit(m)/tspan
        sit2(m) = sit2(m)/tspan - sit(m) * sit(m)
        if (sit2(m) .lt. 0.d0) sit2(m) = 0.d0
        sit2(m) =  sqrt(sit2(m))
      end do

      if (ti_mode .eq. 0) then
        if (.not. usemidpoint) then
          if (isgld .gt. 0) call sgld_avg(.true.)
        end if
        write(mdout, 540) nvalid
        call prntmd(total_nstep, total_nstlim, t, sit, fac, 0, .false., mdloop)
        if (.not. usemidpoint) then
          if (nmropt .ne. 0) call nmrptx(mdout)
          if (isgld>0)call sgld_fluc(.true.)
        end if ! not usemidpoint
        write(mdout, 550)
        call prntmd(total_nstep, total_nstlim, t, sit2, fac, 0, .true., mdloop)
      else
        !still bunk for GPU
        call ti_calc_avg_ene(tspan, .true.)
        do ti_region = 1, 2
          if ( emil_sc_lcl .ne. 1) write(mdout, 601) ti_region

          ! Average energy
          if (( emil_sc_lcl .ne. 1) .or. ( ti_region .eq. 1 )) write(mdout, 540) nvalid

          sit(si_tot_ene) = ti_ene_aug_tmp(ti_region,ti_tot_tot_ene)
          sit(si_kin_ene) = ti_ene_aug_tmp(ti_region,ti_tot_kin_ene)
          sit(si_density) = ti_ene_tmp(ti_region,si_density)
          ti_temp = ti_ene_aug_tmp(ti_region,ti_tot_temp)

          if (( emil_sc_lcl .ne. 1) .or. ( ti_region .eq. 1 )) then
            call prntmd(total_nstep, total_nstlim, t, sit, fac, 0, .false., &
              mdloop)
            if (.not. usemidpoint) then
              if (nmropt .ne. 0) call nmrptx(mdout)
            end if ! not usemidpoint

            call ti_print_ene(ti_region, ti_ene_tmp(ti_region,:), &
              ti_ene_aug_tmp(ti_region,:))
            ! RMSD
            write(mdout, 550)
          end if

          sit2(si_tot_ene) = ti_ene_aug_tmp2(ti_region,ti_tot_tot_ene)
          sit2(si_kin_ene) = ti_ene_aug_tmp2(ti_region,ti_tot_kin_ene)
          sit2(si_density) = ti_ene_tmp2(ti_region,si_density)
          ti_temp = ti_ene_aug_tmp2(ti_region,ti_tot_temp)

          if ( ( emil_sc_lcl .ne. 1) .or. ( ti_region .eq. 1 ) ) then
            call prntmd(total_nstep, total_nstlim, t, sit2, &
              fac, 0, .true., mdloop)
            call ti_print_ene(ti_region, ti_ene_tmp2(ti_region,:), &
              ti_ene_aug_tmp2(ti_region,:))
          end if

        end do

        ! Following sander, we only output the average DV/DL contributions
        write(mdout, 541) nvalid
        call prntmd(total_nstep, total_nstlim, t, ti_ene_delta_tmp, fac, 0, &
          .true., mdloop)

        if (logdvdl .ne. 0) call ti_print_dvdl_values
      end if

      if (.not. usemidpoint) then
        if (nmropt .ne. 0) then
          write(mdout, '(/,a,/)') ' NMR restraints on final step:'
          call ndvptx(crd, frc, mdout)
        end if
      end if

      ! Print Born radii statistics:

      if (using_gb_potential) then
        if (rbornstat .eq. 1) then
          ! Born radii stats collected every nrespai step not nrespa step
          tspan = nvalidi
          write(mdout, 580) nstep
          write(mdout, 590)

          if (.not. usemidpoint) then
            do m = 1, atm_cnt
              gbl_rbave(m) = gbl_rbave(m) / tspan
              gbl_rbfluct(m) = gbl_rbfluct(m) / tspan - &
                gbl_rbave(m) * gbl_rbave(m)
              gbl_rbfluct(m) = sqrt(gbl_rbfluct(m))
              write(mdout, 600) m, gbl_rbmax(m), gbl_rbmin(m), &
                gbl_rbave(m), gbl_rbfluct(m)
            end do
          end if

        end if
      end if

    end if ! (nvalid .gt. 0)

    if (ntp .ne. 0 .and. barostat .eq. 2) call mcbar_summary

  end if ! master

  if (ti_mode .ne. 0) call ti_cleanup

  call update_time(runmd_time)

  return

  540 format(/5x, ' A V E R A G E S   O V E R ', i7, ' S T E P S', /)
  541 format(/5x,' DV/DL, AVERAGES OVER ',i7,' STEPS',/)
  542 format('|',79('='))
  550 format(/5x, ' R M S  F L U C T U A T I O N S', /)
  580 format('STATISTICS OF EFFECTIVE BORN RADII OVER ',i7,' STEPS')
  590 format('ATOMNUM     MAX RAD     MIN RAD     AVE RAD     FLUCT')
  600 format(i4, 2x, 4f12.4)
  601 format(/,'| TI region ', i2, /)

end subroutine runreweight

!*******************************************************************************
!
! Subroutine:   get_position
!
! Description:  Find the center of a set of atoms, and the extent.
!
! Input:  n  number of atoms
!         x  coordinates
! Output: xc, yc, zc  coordinates geometric center
!         e   xmin, ymin, zmin, xmax, ymax, zmax
!
!*******************************************************************************

subroutine get_position(n, x, xc, yc, zc, e)

  use file_io_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: n
  double precision      :: x(3,n)
  double precision      :: xc, yc, zc
  double precision      :: e(3,2)

! Local variables:

  integer               :: i, j

  xc = 0.d0
  yc = 0.d0
  zc = 0.d0

  do i = 1, 3
    e(i,1) = x(i,1)
    e(i,2) = x(i,1)
  end do

  do i = 1, n
    do j = 1, 3
      e(j,1) = min(e(j,1), x(j,i))
      e(j,2) = max(e(j,2), x(j,i))
    end do
  end do

  xc = (e(1,2) + e(1,1)) * 0.5d0
  yc = (e(2,2) + e(2,1)) * 0.5d0
  zc = (e(3,2) + e(3,1)) * 0.5d0

  return

end subroutine get_position

end module runreweight_mod
