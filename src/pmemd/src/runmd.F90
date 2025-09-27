#include "copyright.i"
#include "dbg_arrays.i"

!*******************************************************************************
!
! Module: runmd_mod
!
! Description: Sander 9-compatible MD.
!
!*******************************************************************************

module runmd_mod

  implicit none

contains

!*******************************************************************************
!
! Subroutine:  runmd
!
! Description: Driver routine for molecular dynamics.
!
!*******************************************************************************

#ifdef MPI
subroutine runmd(atm_cnt, crd, mass, frc, vel, last_vel, my_atm_lst, &
                 local_remd_method, local_numexchg)
#else
subroutine runmd(atm_cnt, crd, mass, frc, vel, last_vel)
#endif

  use barostats_mod, only : mcbar_trial, mcbar_summary
  use cit_mod
  use constantph_mod
  use constantph_dat_mod, only : on_cpstep
  use constante_mod
  use constante_dat_mod, only : on_cestep
  use get_cmdline_mod, only : cpein_specified
  use constraints_mod
  use degcnt_mod
  use dynamics_mod
  use dynamics_dat_mod
  use emap_mod, only: emap_move
#ifdef EMIL
  use emil_mod
#endif
  use extra_pnts_nb14_mod
  use gb_force_mod
  use gb_ene_mod
  use gbl_constants_mod, only : KB
  use pbc_mod
  use pme_force_mod
  use file_io_mod
  use file_io_dat_mod
  use img_mod
  use loadbal_mod
  use mcres_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use mol_list_mod
  use nb_pairlist_mod
  use nmr_calls_mod
  use parallel_dat_mod
  use parallel_mod
  use gb_parallel_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
#ifdef MPI
  use multipmemd_mod!, only : free_comms   test by DG
  use remd_mod
  use remd_exchg_mod
#endif
#ifndef MPI
#ifndef NOXRAY
  use xray_globals_module, only : real_kind, num_hkl, ntwsf
  use xray_interface2_module, only: xray_get_f_calc => get_f_calc
#endif
#endif
  use resamplekin_mod
  use ramd_mod
  use runfiles_mod
  use shake_mod
  use timers_mod
  use state_info_mod
  use bintraj_mod
  use amd_mod
  use gamd_mod
  use scaledMD_mod
  use ti_mod
#ifdef TIDECOMP
  use ti_decomp_mod
#endif
  use sams_mod  
#ifdef GTI
  use gti_mod
#endif
! Self-Guided molecular/Langevin Dynamics (SGLD)
#ifdef CUDA
  use sgld_mod, only : trxsgld,sgld_init_ene_rec,sgenergy,sgld_avg,sgld_fluc, &
      sgldw,sgmdw, sg_fix_degree_count,sggamma,com0sg,com1sg,com2sg,     &
           gpu_sgenergy
#else
  use sgld_mod, only : trxsgld,sgld_init_ene_rec,sgenergy,sgld_avg,sgld_fluc, &
      sgldw,sgmdw,sg_fix_degree_count,sggamma,com0sg,com1sg,com2sg,templf,temphf  
#endif
#ifdef MPI
  use sgld_mod, only :sg_allgather
#endif
! Modified by Feng Pan
#ifdef MPI
  use nfe_setup_mod, only : nfe_on_mdstep => on_mdstep, nfe_on_mdwrit => on_mdwrit
  use nfe_lib_mod,   only : nfe_prt, nfe_real_mdstep
  use remd_mod,      only : remd_types, replica_indexes, remd_volume, remd_pressure, use_pv
  use nfe_abmd_mod,  only : abmd_mdstep => mdstep, abmd_selection_freq => selection_freq
  use nfe_bbmd_mod,  only : bbmd_mdstep => mdstep, bbmd_exchg_freq => &
                            exchange_freq,bbmd_active => active
#else
  use nfe_setup_mod, only : nfe_on_mdwrit => on_mdwrit
  use nfe_lib_mod,   only : nfe_prt, nfe_real_mdstep
#endif

!Adde by DG
#ifdef MPI
  use nebread_mod
  use neb_mod
#endif
  use processor_mod
  use parallel_processor_mod
  use pme_recip_midpoint_mod

#ifdef DBG_ARRAYS
  use dbg_arrays_mod
#endif

! ASM
#ifdef MPI
  use asm_mod, only : asm_define, asm_calc
#endif
! END_ASM

! used by middle-scheme
  use md_scheme
!---------------------------------------------------------------

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: crd(3, atm_cnt)
  double precision      :: mass(atm_cnt)
  double precision      :: frc(3, atm_cnt)
  double precision      :: vel(3, atm_cnt)
  double precision      :: last_vel(3, atm_cnt)
  double precision      :: onstepvel(3, atm_cnt)
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
  double precision      :: target_ekin               !< Target EKIN, used in Bussi's thermostat
  integer               :: target_ekin_update_nstep  !< Target EKIN update rate, used in Bussi's thermostat
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

  double precision      :: temp

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
  logical               :: dump_crds
  logical               :: dump_rest
  logical               :: dump_vels
  logical               :: save_onstepvel
  logical               :: timlim_exceeded

  logical               :: update_bussi_target_kin_energy_on_current_step
  logical               :: update_kin_energy_on_current_step

#ifdef MPI
  logical               :: collect_crds
  logical               :: collect_vels
  logical               :: all_crds_valid
  logical               :: all_vels_valid
  logical               :: print_exch_data
#endif

  logical               :: mcresaccept

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


! Variables used for GaMD
  !integer               :: igamd0
  logical               :: update_gamd
  double precision      :: VmaxDt,VminDt,VavgDt,sigmaVDt,VPt
  double precision      :: VmaxPt,VminPt,VavgPt,sigmaVPt,VDt
  ! for triple-boost
  double precision      :: VmaxBt,VminBt,VavgBt,sigmaVBt
  ! LiGaMD
  double precision      :: pot_ene_ll, pot_ene_l_pe, pot_ene_notl       ! not used
  integer               :: blig,blig_min,blig0,atm_i,atm_j,nfmsd,ntwin,ift
  double precision      :: dlig,dmin,pos_tmp(3),fr(3)
  double precision, allocatable      :: lcrd0(:,:,:),lcrd1(:,:,:)

  integer,save          :: counts=0

! PLUMED
  double precision      :: plumed_box(3,3), plumed_virial(3,3), plumed_kbt
  integer               :: plumed_version, plumed_stopflag
  double precision      :: plumed_energyUnits, plumed_timeUnits, plumed_lengthUnits
  double precision      :: plumed_chargeUnits
  integer               :: plumed_need_pot_enes
  double precision      :: plumed_frc(3,atm_cnt)

#ifndef NOXRAY
! X-ray
#if defined(CUDA) && !defined(MPI)
  integer               :: alloc_status
  real(real_kind), allocatable :: r_Fcalc(:), i_Fcalc(:)
#endif
#endif

! ASM
  double precision      :: asm_e, asm_frc(3,atm_cnt)
! END_ASM

  !SGLD
  double precision      :: sgsta_rndfp, sgend_rndfp, ignore_solvent

! middle-scheme
  double precision      :: xold(3, atm_cnt)

! Runmd operates in kcal/mol units for energy, amu for masses,
! and angstoms for distances.  To convert the input time parameters
! from picoseconds to internal units, multiply by 20.455
! (which is 10.0 * sqrt(4.184)).

  call zero_time()

! Here we set total_nstlim. This will just be set to nstlim to start with, but
! will be set to numexchg * nstlim if we're doing REMD. Printing frequency,
! etc. will be based on total_nstep, which is incremented with nstep but not
! reset upon exchanges (as usual)
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
  mcresaccept = .false.
  belly = ibelly .gt. 0
  lout = .true.
  use_vlimit = vlimit .gt. 1.0d-7
  ekmh = 0.d0

  tmassinv = 1.d0 / tmass       ! total mass inverse for all atoms.

! If ntwprt .ne. 0, only print the solute atoms in the coordinate/vel. archives.

  nrx  = natom * 3
  if (ntwprt .gt. 0) nrx = ntwprt * 3

   rndfp=0.0
   rndfs=0.0
  
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
    if (.not. allocated(atm_igroup)) allocate(atm_igroup(atm_cnt))
    call degcnt(ibelly, natom, atm_igroup, natom, gbl_bond, gbl_bond(bonda_idx), &
                ntc, rndfp, rndfs)
  end if


  rndfp = rndfp - dble(ndfmin) + dble(num_noshake)

  rndf = rndfp + rndfs  ! total degrees of freedom in system

! Correct the total degrees of freedom for extra points (lone pairs):

  rndf = rndf - 3.d0 * dble(numextra)

! BUGBUG - NOTE that rndfp, rndfs are uncorrected in an extra points context!
  !  initialization
  if (isgld > 0) then
    ! number of degrees of freedom in the SGLD part
    if (isgsta == 1) then
      sgsta_rndfp = 0
    else
      call degcnt(ibelly, natom, atm_igroup, isgsta-1, gbl_bond, gbl_bond(bonda_idx), &
      ntc, sgsta_rndfp,ignore_solvent)
    end if
    if (isgend == natom) then
      sgend_rndfp = rndf
    else
      call degcnt(ibelly, natom, atm_igroup, isgend, gbl_bond, gbl_bond(bonda_idx), &
      ntc, sgend_rndfp,ignore_solvent)
    end if
    ! Warning - NOTE that the solute ndf outputs above from degcnt are uncorrected
    ! for qmmm_struct%noshake_overlap, num_noshake, and extra points;
    ! also ndfmin is not always being handled.
    call sg_fix_degree_count(sgsta_rndfp, sgend_rndfp, ndfmin, rndf)
  end if


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

  nvalid = 0
  nvalidi = 0
  nstep = 0
  total_nstep = 0
  steps_per_avg = ene_avg_sampling

  si(:) = 0.d0
  sit(:) = 0.d0
  sit2(:) = 0.d0

! The _tmp variables below are only used if ntave .gt. 0.  They are used for
! scratch space in the ntave calcs, and to hold the last sit* values between
! calls to the ntave code.

  sit_tmp(:) = 0.d0
  sit2_tmp(:) = 0.d0

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
  VminBt = 1.0d99
  VmaxBt = -1.0d99
  VavgBt = 0.0d0
  sigmaVBt = 0.d0

! GaMD settings
  if (igamd.gt.0 .and. ibblig.gt.0) then
     blig0 = 1
     blig = 1
     open(unit=gamdlig,file=gamdlig_name,action='write')
     write(gamdlig,'(a)')'#step     ligand     distance'
     ! allocate dynamic arrays to calculate ligand MSD
     if (ibblig.eq.2) then
        nfmsd = ntmsd/ntwx
        ntwin = ntmsd+nftau*ntwx
        allocate(lcrd0(3,nftau,nlig), lcrd1(3,nftau,nlig))
     end if
  end if
  if ((igamd.gt.0) .and. (ntcmd.gt.0)) then
     igamd = 0
     if (master) write(mdout,'(/,a,i10,/)') &
       '| GaMD: Run initial conventional MD with no boost; igamd = ', igamd
#ifdef CUDA
     call gpu_igamd_update(igamd)
#endif
  end if

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

! PLUMED
  if (plumed == 1) then
#   include "Plumed_init.inc"
  endif

! ASM
#ifdef MPI
  if (asm == 1) call asm_define(crd)
#endif
! END_ASM

! Make a first dynamics step:

  irespa = 1    ! PME respa step counter.
  if (usemidpoint) then
   proc_new_list = .true.
  else
    ! We must build a new pairlist the first time we run force:
    new_list = .true.
  end if

  need_pot_enes = .true.        ! maybe not, but most simple for step 0
  need_virials = (ntp .gt. 0 .and. barostat .ne. 2) .or. csurften .gt. 0

! irest = 0:  General startup if not continuing a previous run:
  if (irest .eq. 0 ) then

    irespa = 0

    if(.not.usemidpoint) then
       if (infe .ne. 0) nfe_real_mdstep = .False.
    end if

! PLUMED
  if (plumed == 1) then
    plumed_stopflag = 0
    plumed_need_pot_enes = 0

    call plumed_f_gcmd("setStep"//char(0), nstep)
    call plumed_f_gcmd("prepareDependencies"//char(0), 0)

    ! May be call to isEnergyNeeded isn't necessary here because of explicit assignment above.
    ! However, it could change in future so I'll let it be here
    call plumed_f_gcmd("isEnergyNeeded"//char(0), plumed_need_pot_enes)

    if (plumed_need_pot_enes > 0) then
      need_pot_enes = .true.
    end if
  end if

    ! Calculate the force.  This also does ekcmt if a regular pme run:

    call update_time(runmd_time)

    if (using_pme_potential) then

#ifdef CUDA
      call gpu_upload_crd(crd)
#endif

#ifdef MPI
      if (usemidpoint) then
        call pme_force_midpoint(atm_cnt, crd, frc, &
                       my_atm_lst, proc_new_list, need_pot_enes, need_virials, &
                       pme_pot_ene, nstep, virial, ekcmt, pme_err_est)
      else
        call pme_force(atm_cnt, crd, frc, gbl_img_atm_map, gbl_atm_img_map, &
                       my_atm_lst, new_list, need_pot_enes, need_virials, &
                       pme_pot_ene, nstep, virial, ekcmt, pme_err_est)
      end if
#else
      if (usemidpoint) then
#ifdef MPI
        call pme_force_midpoint(atm_cnt, crd, frc, &
                                proc_new_list, need_pot_enes, need_virials, &
                                pme_pot_ene, nstep, virial, ekcmt, pme_err_est)
        DBG_ARRAYS_DUMP_3DBLE("runmd_frc1", proc_atm_to_full_list,proc_atm_frc,proc_num_atms)
#endif
      else
        call pme_force(atm_cnt, crd, frc, gbl_img_atm_map, gbl_atm_img_map, &
                       new_list, need_pot_enes, need_virials, &
                       pme_pot_ene, nstep, virial, ekcmt, pme_err_est)
      end if
#endif /* MPI */
#ifdef CUDA
      call gpu_download_frc(frc)
#endif

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
      si(si_phmd_ene) = pme_pot_ene%phmd
      if((igamd.ge.12.and.igamd.le.19).or.(igamd0.ge.12.and.igamd0.le.19))then
      si(si_dvdl)=ppi_inter_ene
      si(si_gamd_ppi)=ppi_dihedral_ene
      endif
      if((igamd.ge.20.and.igamd.le.28).or.(igamd0.ge.20.and.igamd0.le.28))then
         si(si_dvdl)=ppi_inter_ene
         si(si_gamd_ppi)=ppi_dihedral_ene
         si(si_gamd_ppi2)=ppi_bond_ene
      endif
      if((igamd.ge.110.and.igamd.le.120).or.(igamd0.ge.110.and.igamd0.le.120))then
      si(si_dvdl)=ppi_inter_ene
      endif

      ! Total virial will be printed in mden regardless of ntp value.
      si(si_tot_virial) = virial(1) + virial(2) + virial(3)

    else if (using_gb_potential) then
#ifdef CUDA
      call gpu_upload_crd(crd)
#endif
      if (.not. usemidpoint) then
        call gb_force(atm_cnt, crd, frc, gb_pot_ene, irespa, need_pot_enes)
      end if
#ifdef CUDA
      call gpu_download_frc(frc)
#endif

      si(si_pot_ene) = gb_pot_ene%total
      si(si_vdw_ene) = gb_pot_ene%vdw_tot
      si(si_elect_ene) = gb_pot_ene%elec_tot
      si(si_hbond_ene) = gb_pot_ene%gb !hbonds are not supported in GB, so reuse the slot here
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
      si(si_efield_ene) = 0.d0
      si(si_phmd_ene) = gb_pot_ene%phmd
    end if

    if (ti_mode .ne. 0) si(si_dvdl) = ti_ene(1, si_dvdl)

#if defined(CUDA)

#if defined(GTI)

!    if (ti_mode .ne. 0 .and. (igamd.gt.0 .or. igamd0.gt.0)) &
    if (igamd.eq.0 .and. igamd0.eq.0) &
    call gti_sync_vector(0, 0) ! copy force over
    call gpu_calculate_kinetic_energy(c_ave, eke, ekph, ekpbs)
    eke=ekph
    if (ti_mode .ne. 0) then
      call gti_kinetic(c_ave)
      call gti_update_kin_energy_from_gpu(ti_ekph, ti_eke)
      ti_ene(1:2,si_kin_ene) = ti_kin_ene(1:2,ti_eke)
    end if
#endif  /* GTI */


    if (ti_mode .ne. 0) then
      si(si_dvdl) = ti_ene(1, si_dvdl)
!not ideal, but seems to be the best way to maintain encapsulation for dvdl
!while only added a minor conditional for gas phase ti
      if (ti_mode .gt. 1 .and. using_gb_potential) si(si_dvdl) = si(si_dvdl) + gb_pot_ene%dvdl
    end if
#else /* CUDA */
!we do the force sync in gpu_ti_pme_ene/force
    if (ti_mode .ne. 0) then
      if (ti_latm_cnt(1) .gt. 0) then
        call ti_exchange_vec(atm_cnt, frc, .true.)
      end if
    end if
#endif /* CUDA */

    ! This force call does not count as a "step". call nmrdcp to decrement
    ! local NMR step counter:
    if (.not. usemidpoint) then
      call nmrdcp
    end if

! PLUMED
    if (plumed == 1) then
#ifdef CUDA
      call gpu_download_crd(crd)

      if (plumed_need_pot_enes > 0) then
        call gpu_download_frc(frc)
      end if
#endif

#     include "Plumed_force.inc"

#ifdef CUDA
      if (plumed_need_pot_enes > 0) then
        call gpu_upload_frc(frc)
      else
        call gpu_upload_frc_add(plumed_frc)
      end if
#endif
    end if

! ASM
#ifdef MPI
    if (asm == 1) then
#ifdef CUDA
      call gpu_download_crd(crd)
#endif
      call asm_calc(crd, asm_frc, asm_e)
#ifdef CUDA
      call gpu_upload_frc_add(asm_frc)
#else
      frc = frc + asm_frc
#endif
    end if
#endif
! END_ASM

    irespa = 1

    ! The coordinates will not be changed between here and the next
    ! run of force, so we can just set new_list to .false.

    if (usemidpoint) then
      proc_new_list = .false.
    else
      new_list = .false.
    end if

    ! Reset quantities depending on temp0 and tautp (which may have been
    ! changed by modwt during force call):

    ekin0  = fac(1) * temp0
    target_ekin = ekin0
    target_ekin_update_nstep = tautp / dt

    if (ti_mode .ne. 0) then
      ti_ekin0(1) = ti_fac(1,1) * temp0
      ti_ekin0(2) = ti_fac(2,1) * temp0
    end if

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

    end if

    eke = 0.d0

#ifdef GBTimer
    call get_wall_time(wall_s, wall_u)
    strt_time_ms = dble(wall_s) * 1000.0d0 + dble(wall_u) / 1000.0d0
#endif

#if defined(MPI) && !defined(CUDA)
    if (usemidpoint) then
      do j = 1, proc_num_atms
        eke = eke + proc_atm_mass(j) * &
               (proc_atm_vel(1,j) * proc_atm_vel(1,j) + &
                proc_atm_vel(2,j) * proc_atm_vel(2,j) + &
                proc_atm_vel(3,j) * proc_atm_vel(3,j))
      end do
    else
      do atm_lst_idx = 1, my_atm_cnt
        j = my_atm_lst(atm_lst_idx)
        eke = eke + mass(j) * &
               (vel(1,j) * vel(1,j) + vel(2,j) * vel(2,j) + vel(3,j) * vel(3,j))
      end do
    end if
#else
    do j = 1, atm_cnt
      eke = eke + mass(j) * &
             (vel(1,j) * vel(1,j) + vel(2,j) * vel(2,j) + vel(3,j) * vel(3,j))
    end do
#endif

    eke = 0.5d0 * eke


#ifdef GBTimer
    call get_wall_time(wall_s, wall_u)
    print *, "eke +  mass =", dble(wall_s) * 1000.0d0 + &
       dble(wall_u) / 1000.0d0 - strt_time_ms
#endif

    ! New addition - is this right?
    if (ti_mode .eq. 0) then
#if defined(MPI) && !defined(CUDA)
      call update_time(runmd_time)
! Sum up the partial kinetic energies. The all_reduce is expensive, but only
! occurs once in an irest .eq. 0 run.

      reduce_buf_in(1) = eke
      call mpi_allreduce(reduce_buf_in, reduce_buf_out, 1, &
        mpi_double_precision, mpi_sum, pmemd_comm, err_code_mpi)

      eke = reduce_buf_out(1)
      call update_time(fcve_dist_time)
#endif
    else
#if defined(MPI) && !defined(CUDA)
      call ti_calc_kin_ene(atm_cnt, my_atm_cnt, my_atm_lst, mass, vel, vel, &
                           1.d0, ti_eke)
      call update_time(runmd_time)

      reduce_buf_in(1) = eke
      reduce_buf_in(2) = ti_kin_ene(1,ti_eke)
      reduce_buf_in(3) = ti_kin_ene(2,ti_eke)
      reduce_buf_in(4) = ti_kin_ene(1,ti_sc_eke)
      reduce_buf_in(5) = ti_kin_ene(2,ti_sc_eke)
      call mpi_allreduce(reduce_buf_in, reduce_buf_out, 5, &
        mpi_double_precision, mpi_sum, pmemd_comm, err_code_mpi)

      eke = reduce_buf_out(1)
      ti_ene(1,si_kin_ene) = reduce_buf_out(2)
      ti_ene(2,si_kin_ene) = reduce_buf_out(3)
      ti_kin_ene(1,ti_sc_eke) = reduce_buf_out(4)
      ti_kin_ene(2,ti_sc_eke) = reduce_buf_out(5)
      call update_time(fcve_dist_time)
#else
      call ti_calc_kin_ene(atm_cnt, mass, vel, vel, 1.d0, ti_eke)
      ti_ene(1,si_kin_ene) = ti_kin_ene(1,ti_eke)
      ti_ene(2,si_kin_ene) = ti_kin_ene(2,ti_eke)
#endif /* #if defined(GTI) && defined(CUDA)  */
    end if

    if (ti_mode .ne. 0) then
      ti_kin_ene(1,ti_eke) = eke - ti_ene(2,si_kin_ene)
      ti_kin_ene(2,ti_eke) = eke - ti_ene(1,si_kin_ene)

      ti_ene_aug(1,ti_tot_kin_ene) = ti_kin_ene(1,ti_eke)
      ti_ene_aug(2,ti_tot_kin_ene) = ti_kin_ene(2,ti_eke)
      ti_ene_aug(1,ti_tot_temp) = ti_kin_ene(1,ti_eke) / ti_fac(1,1)
      ti_ene_aug(2,ti_tot_temp) = ti_kin_ene(2,ti_eke) / ti_fac(2,1)
      ti_ene_aug(1,ti_tot_tot_ene) = ti_kin_ene(1,ti_eke) + si(si_pot_ene)
      ti_ene_aug(2,ti_tot_tot_ene) = ti_kin_ene(2,ti_eke) + si(si_pot_ene)

      ti_ene(1,si_kin_ene) = ti_kin_ene(1,ti_sc_eke)
      ti_ene(2,si_kin_ene) = ti_kin_ene(2,ti_sc_eke)

      ti_ene(1,si_tot_ene) = ti_ene(1,si_kin_ene) + ti_ene(1,si_pot_ene)
      ti_ene(2,si_tot_ene) = ti_ene(2,si_kin_ene) + ti_ene(2,si_pot_ene)
    end if

    si(si_solute_kin_ene) = eke
    si(si_kin_ene) = si(si_solute_kin_ene)
    si(si_tot_ene) = si(si_kin_ene) + si(si_pot_ene)
    eke = 0.d0

#ifdef CUDA
!this seems unnecessary, as the forces appear to be the same as at line 462.
    call gpu_download_frc(frc)
#endif
#ifdef GBTimer
    call get_wall_time(wall_s, wall_u)
    strt_time_ms = dble(wall_s) * 1000.0d0 + dble(wall_u) / 1000.0d0
#endif

    if (usemidpoint) then
#ifdef MPI
      do j = 1, proc_num_atms
        winf = 1.0d0/proc_atm_mass(j) * half_dtx
        do i = 1, 3
          proc_atm_vel(i,j) = proc_atm_vel(i,j) - proc_atm_frc(i,j)*winf
          if (use_vlimit) vel(i, j) = sign(min(abs(vel(i, j)), vlimit), vel(i, j))
        end do
      end do
#endif /*MPI*/
    else ! usemidpoint
#if defined(MPI) && !defined(CUDA)
      do atm_lst_idx = 1, my_atm_cnt
        j = my_atm_lst(atm_lst_idx)
#else
      do j = 1, atm_cnt
#endif
        winf = atm_mass_inv(j) * half_dtx
        do i = 1, 3
          vel(i, j) = vel(i, j) - frc(i, j) * winf
#ifndef CUDA
          !Vlimit not supported on GPUs
          if (use_vlimit) vel(i, j) = sign(min(abs(vel(i, j)), vlimit), vel(i, j))
#endif
        end do
      end do
    end if ! usemidpoint


#ifdef CUDA
#ifdef MPI
    if (ineb .gt. 0 .and. (beadid .eq. 1 .or. beadid .eq. neb_nbead)) then
      do i=1, atm_cnt
        vel(:,i)=0.d0
      end do
    end if
#endif
#endif

#ifdef CUDA
    call gpu_upload_vel(vel)
    if (ischeme == 1) call gpu_rattle(dt)
#endif
  if (ischeme == 1) then
            if (usemidpoint) then
#ifdef MPI
                  call comm_3dimensions_3dbls(proc_atm_vel)
                call rattle(proc_atm_crd, proc_atm_vel)
                call rattle_fastwater(proc_atm_crd, proc_atm_vel)
                  !call comm_3dimensions_3dbls(proc_atm_vel)

                  !--to be fixed
                  DBG_ARRAYS_DUMP_CRD("runmd_crd2", proc_atm_to_full_list, proc_atm_crd,proc_num_atms, proc_atm_wrap)
#endif
            else
                call rattle(crd, vel)
                call rattle_fastwater(crd, vel)
            end if
    end if


#ifdef GBTimer
    call get_wall_time(wall_s, wall_u)
    print *, "vel +  winf =", dble(wall_s) * 1000.0d0 + &
          dble(wall_u) / 1000.0d0 - strt_time_ms
#endif

    if (ntt == 1 .or. ntt == 11) then
      ekmh = max(si(si_solute_kin_ene), fac(1) * 10.d0)
      if (ti_mode .ne. 0) then
        ti_kin_ene(1,ti_ekmh) =  max(ti_kin_ene(1,ti_eke), ti_fac(1,1) * 10.d0)
        ti_kin_ene(2,ti_ekmh) =  max(ti_kin_ene(2,ti_eke), ti_fac(2,1) * 10.d0)
      end if
    end if

  end if !irest .eq. 0

#ifdef EMIL
  if ( emil_do_calc .gt. 0 ) then
        call emil_init( atm_cnt, nstep, 1.0/(temp0 * 2 * boltz2 ), &
                     mass, crd, frc, vel)
  end if
#endif

! irest = 1:  Continuation of a previous trajectory:
!
! Note: if the last printed energy from the previous trajectory was
!       at time "t", then the restrt file has velocities at time
!       t + 0.5dt, and coordinates at time t + dt

#ifdef GBTimer
  call get_wall_time(wall_s, wall_u)
  strt_time_ms = dble(wall_s) * 1000.0d0 + dble(wall_u) / 1000.0d0
#endif

  ekmh = 0.0d0
#if defined(GTI) && defined(CUDA)
  call gpu_calculate_kinetic_energy(c_ave, eke, ekph, ekpbs)
  ekmh=ekph
  if (ti_mode .ne. 0) then
    call gti_kinetic(c_ave)
    call gti_update_kin_energy_from_gpu(ti_ekph, ti_ekmh)
    ti_ene(1:2,si_kin_ene) = ti_kin_ene(1:2,ti_ekmh)
  end if
#else

  if (usemidpoint) then
#ifdef MPI
    do i = 1, proc_num_atms
      ekmh = ekmh + proc_atm_mass(i) * &
             (proc_atm_vel(1,i) * proc_atm_vel(1,i) + proc_atm_vel(2,i) &
              * proc_atm_vel(2,i) + proc_atm_vel(3,i) * proc_atm_vel(3,i))
    end do
#endif /*MPI*/
  else
#if defined(MPI) && !defined(CUDA)
    do atm_lst_idx = 1, my_atm_cnt
      i = my_atm_lst(atm_lst_idx)
#else
    do i = 1, atm_cnt
#endif
      ekmh = ekmh + mass(i) * &
             (vel(1,i) * vel(1,i) + vel(2,i) * vel(2,i) + vel(3,i) * vel(3,i))
    end do
  end if

#ifdef GBTimer
  call get_wall_time(wall_s, wall_u)
  print *, "ekmh +  mass =", dble(wall_s) * 1000.0d0 + &
       dble(wall_u) / 1000.0d0 - strt_time_ms
#endif

  if (ti_mode .eq. 0) then
#if defined(MPI) && !defined(CUDA)
! Sum up ekmh. The all_reduce is expensive, but only occurs once per run.

    call update_time(runmd_time)
    reduce_buf_in(1) = ekmh

    call mpi_allreduce(reduce_buf_in, reduce_buf_out, 1, mpi_double_precision, &
                       mpi_sum, pmemd_comm, err_code_mpi)

    ekmh = reduce_buf_out(1)
    call update_time(fcve_dist_time)
#endif
    ekmh = ekmh * 0.5d0
  else
#if defined(MPI) && !defined(CUDA)
! Sum up ekmh. The all_reduce is expensive, but only occurs once per run.
    call ti_calc_kin_ene(atm_cnt, my_atm_cnt, my_atm_lst, mass, vel, vel, &
                         1.d0, ti_ekmh)

    call update_time(runmd_time)

    reduce_buf_in(1) = ekmh
    reduce_buf_in(2) = ti_kin_ene(1,ti_ekmh)
    reduce_buf_in(3) = ti_kin_ene(2,ti_ekmh)
    reduce_buf_in(4) = ti_kin_ene(1,ti_sc_eke)
    reduce_buf_in(5) = ti_kin_ene(2,ti_sc_eke)

    call mpi_allreduce(reduce_buf_in, reduce_buf_out, 5, mpi_double_precision, &
                     mpi_sum, pmemd_comm, err_code_mpi)

    ekmh = reduce_buf_out(1)
    ti_ene(1,si_kin_ene) = reduce_buf_out(2)
    ti_ene(2,si_kin_ene) = reduce_buf_out(3)
    ti_kin_ene(1,ti_sc_eke) = reduce_buf_out(4)
    ti_kin_ene(2,ti_sc_eke) = reduce_buf_out(5)

    call update_time(fcve_dist_time)
#else
    call ti_calc_kin_ene(atm_cnt, mass, vel, vel, 1.d0, ti_ekmh)
    ti_ene(1,si_kin_ene) = ti_kin_ene(1,ti_ekmh)
    ti_ene(2,si_kin_ene) = ti_kin_ene(2,ti_ekmh)
    ekmh = ekmh * 0.5d0
#endif
  end if
#endif  /* #if defined(GTI) && defined(CUDA) */

  if (ti_mode .ne. 0) then
    ti_kin_ene(1,ti_ekmh) = ekmh - ti_ene(2,si_kin_ene)
    ti_kin_ene(2,ti_ekmh) = ekmh - ti_ene(1,si_kin_ene)

    ti_ene(1,si_kin_ene) = ti_kin_ene(1,ti_sc_eke)
    ti_ene(2,si_kin_ene) = ti_kin_ene(2,ti_sc_eke)
    ti_ene(1,si_tot_ene) = ti_ene(1,si_kin_ene) + ti_ene(1,si_pot_ene)
    ti_ene(2,si_tot_ene) = ti_ene(2,si_kin_ene) + ti_ene(2,si_pot_ene)
  end if
  if (usemidpoint) then
#ifdef MPI
    proc_atm_last_vel(:,:) = proc_atm_vel(:,:)
#endif /*MPI*/
  else ! usemidpoint
    last_vel(:,:) = vel(:,:)
  end if ! usemidpoint
#ifdef CUDA
  call gpu_upload_last_vel(last_vel)
#endif

  if (irest .eq. 0) then

    ! Print the initial energies and temperatures

    if (nstep .le. 0 .and. master) then
      if (.not. usemidpoint) then
        if (isgld > 0) call sgenergy(si(si_pot_ene))
      end if
      if (ti_mode .eq. 0) then
        call prntmd(nstep, total_nstlim, t, si, fac, 7, .false., mdloop)
      else
        if (ifmbar .ne. 0 .and. do_mbar) then
          call ti_print_mbar_ene(si(si_pot_ene))
        end if
        do ti_region = 1, 2

          if ( emil_sc_lcl .ne. 1 ) write(mdout, 601) ti_region

          si(si_tot_ene) = ti_ene_aug(ti_region,ti_tot_tot_ene)
          si(si_kin_ene) = ti_ene_aug(ti_region,ti_tot_kin_ene)
          si(si_density) = ti_ene(ti_region,si_density)
          ti_temp = ti_ene_aug(ti_region,ti_tot_temp)

          if ( (emil_sc_lcl .ne. 1) .or. (ti_region .eq. 1) ) then
            call prntmd(nstep, total_nstlim, t, si, fac, 7, .false., mdloop)
            call ti_print_ene(ti_region, ti_ene(ti_region,:), &
              ti_ene_aug(ti_region,:))
            if (igamd.gt.0 .or. igamd0.gt.0) &
              call ti_others_print_ene(ti_region, ti_others_ene(ti_region,:))
          end if
        end do
#ifdef GTI
        call gti_output_ti_result(nstep*t)        
#endif
      end if
      if (.not. usemidpoint) then
        if (nmropt .ne. 0) call nmrptx(mdout)
        if (infe .ne. 0) call nfe_prt(mdout)
      end if
    end if
    if (nstlim .eq. 0) return

  end if

  ! Print out the initial state information to cpout (step 0)
  if (.not. usemidpoint) then
    if ((icnstph .eq. 1 .or. (icnste .eq. 1 .and. cpein_specified)) .and. local_remd_method .eq. 0) &
      call cnstph_write_cpout(0, total_nstlim, t, local_remd_method, remd_types, replica_indexes)

    ! Print out the initial state information to cpout (step 0)
    if (icnste .eq. 1 .and. local_remd_method .eq. 0 .and. .not. cpein_specified) &
      call cnste_write_ceout(0, total_nstlim, t, local_remd_method, remd_types, replica_indexes)
  end if

#ifdef MPI
  !=============================================================================
  ! REMD LOOP.
  ! We loop here around the main loop for performing the dynamics. Before the
  ! main loop, we perform the exchange attempts.
  !=============================================================================

  ! Make sure we go through the runmd loop at least once

  if (local_numexchg .eq. 0) local_numexchg = local_numexchg + 1

  ! Some basic REMD setup that only has to be done once

  if (local_remd_method .ne. 0) then
    total_nstlim = nstlim * local_numexchg
  end if

  do mdloop = 1, local_numexchg

    ! Reset the step counters

    nstep = 0
    irespa = 1

!Init Energy record for SGLD
    if (.not. usemidpoint) then
      if (isgld.gt.0) call sgld_init_ene_rec
    end if

    if (local_remd_method .ne. 0) then
      ! We need to get the potential energies again, since the positions
      ! have been updated, so our energies no longer correspond to our
      ! coordinates

      print_exch_data = mod(mdloop, exch_ntpr) .eq. 0 .or. &
                        mdloop .eq. local_numexchg

      need_pot_enes = .true.

      if (using_pme_potential) then
        if (usemidpoint) then
#ifdef MPI
          call pme_force_midpoint(atm_cnt, crd, frc, &
                         my_atm_lst, proc_new_list, need_pot_enes, need_virials, &
                         pme_pot_ene, nstep, virial, ekcmt, pme_err_est)
          DBG_ARRAYS_DUMP_3DBLE("runmd_frc2", proc_atm_to_full_list,proc_atm_frc,proc_num_atms)
#endif
        else ! usemidpoint
          call pme_force(atm_cnt, crd, frc, gbl_img_atm_map, gbl_atm_img_map, &
                         my_atm_lst, new_list, need_pot_enes, need_virials, &
                         pme_pot_ene, nstep, virial, ekcmt, pme_err_est)
        end if ! usemidpoint
        remd_ptot = pme_pot_ene%total

      else if (using_gb_potential) then
        if (.not. usemidpoint) then
          call gb_force(atm_cnt, crd, frc, gb_pot_ene, irespa, need_pot_enes)
        end if

        remd_ptot = gb_pot_ene%total

      end if

      if (.not. usemidpoint) then
        ! Set pressure/volume in case PV correction is active.
        if (use_pv) then
          remd_pressure = pres0
          ! On first exchange no dynamics have been done, use uc_volume. Otherwise
          ! use the last recorded volume.
          if (mdloop.eq.1) then
          !  ! Use a fake scale to calculate the volume
          !  !rmu(1:3) = 1.d0
          !  !call pressure_scale_pbc_data(rmu, vdw_cutoff + skinnb, verbose)
            remd_volume = uc_volume
          else
            remd_volume = si(si_volume)
          endif
        endif
        select case (local_remd_method)
        case(-1)
          ! similar to what we do below because this may
          ! include hamiltonian_exchange or reservoir_exchange()
          if (using_pme_potential) then
            call mpi_allgathervec(atm_cnt, frc)
            call mpi_allgathervec(atm_cnt, vel)
          else if (using_gb_potential) then
            call gb_mpi_allgathervec(atm_cnt, frc)
            call gb_mpi_allgathervec(atm_cnt, vel)
          end if

          call multid_exchange(atm_cnt, crd, vel, frc, remd_ptot, &
            si(si_kin_ene)/fac(1), my_atm_lst, &
            gbl_img_atm_map, gbl_atm_img_map, mdloop)
          new_list = .true.
        case(1)
          if (isgld==1)then
            call rxsgld_exchange(atm_cnt, my_atm_lst, mass, crd, vel, remd_ptot, 1, numgroups, &
              si(si_kin_ene) / fac(1), print_exch_data, &
              mdloop)
              new_list = .true.
            else if (rremd_type .gt. 0 .and. hybridgb .le. 0) then
            ! We need to synchronize the forces and velocities since
            ! my_atm_lst may change, crd is OK because it will be done
            ! inside reservoir_exchange()
            if (using_pme_potential) then
              call mpi_allgathervec(atm_cnt, frc)
              call mpi_allgathervec(atm_cnt, vel)
            else if (using_gb_potential) then
              call gb_mpi_allgathervec(atm_cnt, frc)
              call gb_mpi_allgathervec(atm_cnt, vel)
            end if

            call reservoir_exchange(atm_cnt, crd, vel, remd_ptot, 1, numgroups, &
                                    si(si_kin_ene) / fac(1), print_exch_data, &
                                    mdloop)
            new_list = .true.
          else if (hybridgb .gt. 0) then
            ! We need to synchronize the coordinates, forces and velocities since
            ! my_atm_lst may change.
            ! Use mpi_allgathervec() since it hybridgb uses pme_potential.
            call mpi_allgathervec(atm_cnt, frc)
            call mpi_allgathervec(atm_cnt, vel)
            call mpi_allgathervec(atm_cnt, crd)

            call hybridsolvent_exchange(atm_cnt, crd, vel, frc, numwatkeep, &
                          remd_ptot, 1, numgroups, si(si_kin_ene) / fac(1), & 
                          print_exch_data, mdloop)
            new_list = .true.
          else
            call temperature_exchange(atm_cnt, vel, remd_ptot, 1, numgroups, &
                                  si(si_kin_ene) / fac(1), print_exch_data, mdloop)
          end if
        case(3)
          ! We need to synchronize the forces and velocities since
          ! my_atm_lst may change, crd is OK because it will be done
          ! inside hamiltonian_exchange()
          if (using_pme_potential) then
            call mpi_allgathervec(atm_cnt, frc)
            call mpi_allgathervec(atm_cnt, vel)
          else if (using_gb_potential) then
            call gb_mpi_allgathervec(atm_cnt, frc)
            call gb_mpi_allgathervec(atm_cnt, vel)
          end if

          call hamiltonian_exchange(atm_cnt, my_atm_lst, crd, vel, frc, &
            remd_ptot, gbl_img_atm_map, &
            gbl_atm_img_map, 1, numgroups, &
            print_exch_data, mdloop)
          new_list = .true.
        case(4)
          call ph_remd_exchange(1, numgroups, mdloop)
        case(5)
          call e_remd_exchange(1, numgroups, mdloop)
        end select
      end if

#ifndef CUDA
      ! We need to update last_vel (on GPUs this is done correctly already)
      if (usemidpoint) then
        proc_atm_last_vel(:, :) = proc_atm_vel(:, :)
      else ! usemidpoint
        last_vel(:, :) = vel(:, :)
      end if ! usemidpoint
#endif
    end if ! local_remd_method .ne. 0
#endif /* MPI */

    !=======================================================================
    ! MAIN LOOP FOR PERFORMING THE DYNAMICS STEP:
    ! At this point, the coordinates are a half-step "ahead" of the velocities;
    ! the variable EKMH holds the kinetic energy at these "-1/2" velocities,
    ! which are stored in the array last_vel.
    !=======================================================================

    do
      DBG_ARRAYS_TIME_STEP(nstep)

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

        if (icnstph .gt. 0 .or. (icnste .gt. 0 .and. cpein_specified)) then
#ifdef MPI
          on_cpstep = mod(irespa+(mdloop-1)*nstlim, ntcnstph) .eq. 0
#else
          on_cpstep = mod(irespa, ntcnstph) .eq. 0
#endif
        end if
        if (icnste .gt. 0 .and. .not. cpein_specified) then
#ifdef MPI
          on_cestep = mod(irespa+(mdloop-1)*nstlim, ntcnste) .eq. 0
#else
          on_cestep = mod(irespa, ntcnste) .eq. 0
#endif
        end if

        ! If we are running constant pH simulation, we need to determine if we attempt
        ! a constant pH protonation change and if we need to update the cph pairlist

        if (icnstph .gt. 0 .or. (icnste .gt. 0 .and. cpein_specified)) then
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
            if (icnstph .eq. 1 .or. (icnste .eq. 1 .and. cpein_specified)) then
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

        if (on_cestep .and. icnste .eq. 1 .and. on_cpstep .and. icnstph .eq. 1 &
            .and. .not. cpein_specified) then
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

        if (icnste .gt. 0 .and. .not. cpein_specified) then

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

      if (vrand .ne. 0 .and. ntt .eq. 2) then
        need_pot_enes = .true.
      else if (ntwe .gt. 0) then
        need_pot_enes = (mod(total_nstep+1, steps_per_avg) .eq. 0) .or. &
                        (mod(total_nstep+1, ntwe) .eq. 0 .and. onstep)
      else
        need_pot_enes = mod(total_nstep+1, steps_per_avg) .eq. 0
      end if

#ifdef EMIL
      if ( emil_do_calc .gt. 0 ) then
        need_pot_enes = .true.
      end if
#endif

#ifdef CUDA
      if ((mod(irespa, ntpr) .eq. 0) .or. &
          (iamd .gt. 0) .or. (igamd .gt. 0) .or. (icnstph .eq. 1 .and. on_cpstep) &
          .or. (icnste .eq. 1 .and. on_cestep) .or. (icnste .eq. 1 .and. on_cpstep .and. cpein_specified)) then
        need_pot_enes = .true.
      end if
#else
      if ((iamd .gt. 0) .or. (igamd .gt. 0)) then
        need_pot_enes = .true.
      end if

#endif

#ifdef MPI
      if ( ineb .gt. 0 .and. (mod(nstep, nebfreq)==0)) then
        need_pot_enes = .true.
      end if
#endif

#ifdef MPI
      ! nfe_bbmd may use pot
      if (infe .ne. 0 .and. bbmd_active .eq. 1) &
        need_pot_enes = .true.
#endif

      ! With TI we may need the energy more frequently.
      if (ti_mode .ne. 0) then
        if (ifmbar .ne. 0) then
          do_mbar = .false.
          if (mod(nstep+1,bar_intervall) == 0) then
            do_mbar = .true.
            need_pot_enes = .true.
          end if
        end if
        if (logdvdl .ne. 0) then
          need_pot_enes = .true. ! Always need the energy.
        end if
      end if


      ! Need the potential energy for MC barostat steps and SGLD reweighting
      ! calculation
      if (do_mcbar_trial .or. isgld .gt. 0) need_pot_enes = .true.

! PLUMED
      if (plumed == 1) then
        plumed_stopflag = 0
        plumed_need_pot_enes = 0

        call plumed_f_gcmd("setStep"//char(0), nstep)
        call plumed_f_gcmd("prepareDependencies"//char(0), 0)
        call plumed_f_gcmd("isEnergyNeeded"//char(0), plumed_need_pot_enes)

        if (plumed_need_pot_enes > 0) then
          need_pot_enes = .true.
        end if
      end if

      call update_time(runmd_time)

      if (using_pme_potential) then

        ! Monte Carlo water movement and exchange - mcint > 0
        ! Call mcres which will calculate energy with current coord set.
        ! attempts random residue move (solvent only) and then calculates
        ! new energy and applies metropolis criteria.
        if(mcwat .gt. 0) then

          mcint = nmd

          if(mod(nstep,mcint) .eq. 0) then
            call mcres(mcresaccept,mcresstr, atm_cnt, crd, frc, gbl_img_atm_map, &
                       gbl_atm_img_map, new_list, need_pot_enes, &
                       need_virials, pme_pot_ene, nstep, virial, ekcmt, &
                       pme_err_est)
          end if
        end if
#ifdef MPI
        ! If we are doing REMD in NVT ensemble, there's no need to run pme_force 
        ! on the first step after a REMD exchange attempt since we already did
        ! it above. However, if running in NPT ensemble and we exchanged, force
        ! should be called again to ensure the virial is properly updated.
        if (usemidpoint) then
          if (local_remd_method .eq. 0 .or. nstep .gt. 0 .or. use_pv .or. rremd_type .gt. 0) &
            call pme_force_midpoint(atm_cnt, crd, frc, &
            my_atm_lst, proc_new_list, need_pot_enes, need_virials, &
            pme_pot_ene, nstep, virial, ekcmt, pme_err_est)
          DBG_ARRAYS_DUMP_3DBLE("runmd_frc2", proc_atm_to_full_list,proc_atm_frc,proc_num_atms)
        else ! usemidpoint
          if (local_remd_method .eq. 0 .or. nstep .gt. 0 .or. use_pv .or. rremd_type .gt. 0) &
            call pme_force(atm_cnt, crd, frc, gbl_img_atm_map, gbl_atm_img_map, &
            my_atm_lst, new_list, need_pot_enes, need_virials, &
            pme_pot_ene, nstep, virial, ekcmt, pme_err_est)
        end if ! usemidpoint
#else
        call pme_force(atm_cnt, crd, frc, gbl_img_atm_map, gbl_atm_img_map, &
          new_list, need_pot_enes, need_virials, &
          pme_pot_ene, nstep, virial, ekcmt, pme_err_est)
#endif /* MPI */
        mcresaccept = .false.
        counts = counts + 1

        ! Do the MC barostat trial here (the Berendsen scaling is done after the
        ! coordinate propagation, which is just wrong -- this way the 'proper'
        ! energy terms will be printed (corresponding to the new volume if
        ! accepted, and the old volume otherwise)

        ! turn logdvdl off for hte mcbar cycle
        if (do_mcbar_trial) then
          k=logdvdl
          logdvdl=0
        endif

#ifdef MPI
        if (usemidpoint) then
          if (do_mcbar_trial) &
            call mcbar_trial(atm_cnt, crd, frc, mass, my_atm_lst, proc_new_list, &
            pme_pot_ene, verbose, pme_err_est, vel, numextra)
        else ! usemidpoint
          if (do_mcbar_trial) &
            call mcbar_trial(atm_cnt, crd, frc, mass, my_atm_lst, new_list, &
            pme_pot_ene, verbose, pme_err_est, vel, numextra)
        end if ! usemidpoint
#else
        if (do_mcbar_trial) &
          call mcbar_trial(atm_cnt, crd, frc, mass, new_list, pme_pot_ene, &
          verbose, pme_err_est)
#endif
        ! turn logdvdl back 
        if (do_mcbar_trial) logdvdl=k

        if (ramdint .gt. 0) then
          if (mod(nstep, ramdint) .eq. 0) then
            call run_ramd(atm_cnt, crd, frc, mass, new_list, nstep)
          end if
        end if

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
        si(si_phmd_ene) = pme_pot_ene%phmd
        if((igamd.ge.12.and.igamd.le.19).or.(igamd0.ge.12.and.igamd0.le.19))then
          si(si_dvdl)=ppi_inter_ene
          si(si_gamd_ppi)=ppi_dihedral_ene
        endif
        if((igamd.ge.110.and.igamd.le.120).or.(igamd0.ge.110.and.igamd0.le.120))then
         si(si_dvdl)=ppi_inter_ene
        endif
        if((igamd.ge.20.and.igamd.le.28).or.(igamd0.ge.20 .and.igamd0.le.28))then
           si(si_dvdl)=ppi_inter_ene
           si(si_gamd_ppi)=ppi_dihedral_ene
           si(si_gamd_ppi2)=ppi_bond_ene
        endif
        ! Total virial will be printed in mden regardless of ntp value.
        si(si_tot_virial) = virial(1) + virial(2) + virial(3)
        if (iphmd .eq. 3) then
          call updatephmd(dtx,atm_qterm,total_nstep+1,si(si_kin_ene))
        end if

      else if (using_gb_potential) then
        if (.not. usemidpoint) then
#ifdef MPI
          if (local_remd_method .eq. 0 .or. nstep .gt. 0 .or. rremd_type .gt. 0) &
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
        si(si_phmd_ene) = gb_pot_ene%phmd
      end if

      !debug forces
      !#ifdef CUDA
      !      call gpu_download_frc(frc)
      !#endif

      ! Now that we've called force, evaluate the success (or failure) of the
      ! constant pH move and write the result to the output file, if not done yet
      if (.not. usemidpoint) then
        if (on_cpstep .and. (icnstph .eq. 1 .or. (icnste .eq. 1 .and. cpein_specified)) .and. .not. on_cestep) then
          call cnstph_end_step(gb_pot_ene%dvdl, 0, 0)
          call cnstph_write_cpout(total_nstep+1, total_nstlim, &
            t, local_remd_method, remd_types, replica_indexes)
        end if

        ! Now that we've called force, evaluate the success (or failure) of the
        ! constant Redox potential move and write the result to the output file

        if (on_cestep .and. icnste .eq. 1 .and. .not. cpein_specified) then
          call cnste_end_step(gb_pot_ene%dvdl, 0, 0)
          call cnste_write_ceout(total_nstep+1, total_nstlim, &
            t, local_remd_method, remd_types, replica_indexes)
        end if
      end if

      ! PHMD
      if (iphmd .eq. 1 ) then
        call updatephmd(dtx,atm_qterm,total_nstep+1,si(si_kin_ene))
      end if

      ! Determine if coordinates will be written
      dump_crds = .false.
      if (ntwx .gt. 0) then
        if (mod(total_nstep + 1, ntwx) .eq. 0) dump_crds = .true.
      end if

      ! Determine if restart will be written
      dump_rest = .false.
      if (mod(total_nstep + 1, ntwr) .eq. 0) dump_rest = .true.
      ! Do we need a final restart?
      if (nstep + 1 .ge. nstlim) dump_rest = .true.

      ! Determine if velocities will be written
      dump_vels = .false.
      if (ntwv .gt. 0) then
        if (mod(total_nstep + 1, ntwv) .eq. 0) dump_vels = .true.
      else if (dump_crds .and. ntwv .lt. 0) then
        dump_vels = .true.
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
      target_ekin_update_nstep = tautp / dt
      
      update_bussi_target_kin_energy_on_current_step = (ntt==11 .and. target_ekin_update_nstep /=0 .and. &
          mod(total_nstep, target_ekin_update_nstep)==0)
      update_kin_energy_on_current_step = ntt==1 .or. onstep .or. update_bussi_target_kin_energy_on_current_step
      
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
          !if (ti_mode .ne. 0 .and. (igamd.gt.0 .or. igamd0.gt.0)) &
          if (igamd.eq.0 .and. igamd0.eq.0) &
          call gti_sync_vector(0, 0) ! copy force over
#endif
        end if
        ti_ekin0(1) = ti_fac(1,1) * temp0
        ti_ekin0(2) = ti_fac(2,1) * temp0
      end if

! PLUMED
      if (plumed == 1) then
#ifdef CUDA
        call gpu_download_crd(crd)

        if (plumed_need_pot_enes > 0) then
          call gpu_download_frc(frc)
        end if
#endif

#       include "Plumed_force.inc"

#ifdef CUDA
        if (plumed_need_pot_enes > 0) then
          call gpu_upload_frc(frc)
        else
          call gpu_upload_frc_add(plumed_frc)
        end if
#endif
      end if

! ASM
#ifdef MPI
    if (asm == 1) then
#ifdef CUDA
      call gpu_download_crd(crd)
#endif
      call asm_calc(crd, asm_frc, asm_e)
#ifdef CUDA
      call gpu_upload_frc_add(asm_frc)
#else
      frc = frc + asm_frc
#endif
    end if
#endif
! END_ASM

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

      ! Do randomization of velocities, if needed:
      reset_velocities = .false.
      
      ! Assign new random velocities every Vrand steps, if ntt .eq. 2:
      if (vrand .ne. 0 .and. ntt .eq. 2) then
        if (mod((total_nstep+1), vrand) .eq.  0) reset_velocities = .true.
      end if
      
#ifdef MPI
      ! enforce reset_velocities if set (by H-REMD) 
      if (enforce_reset_velocities) then 
        reset_velocities =.true.
        enforce_reset_velocities =.false.
      endif  
#endif

#ifdef GBTimer
      call get_wall_time(wall_s, wall_u)
      strt_time_ms = dble(wall_s) * 1000.0d0 + dble(wall_u) / 1000.0d0
#endif

      if (reset_velocities) then
        if (master) then
          write(mdout,'(a,i8)') 'Setting new random velocities at step ', &
            total_nstep + 1
        end if
#ifdef CUDA
        call gpu_vrand_reset_velocities(temp0 * factt, half_dtx)
#ifdef GTI
        if (ti_mode .ne. 0) then
          ! copy over vel: (gti_syn_mass=0: weighted-combination; gti_syn_mass=1: copy over)
          if (ti_latm_cnt(1) .gt. 0) call gti_sync_vector(1, gti_syn_mass+1)
        end if
#endif
#else
        if (ti_mode .eq. 0) then
          call vrand_set_velocities(atm_cnt, vel, atm_mass_inv, temp0 * factt)
        else
          call ti_vrand_set_velocities(atm_cnt, vel, atm_mass_inv, temp0)
          if (ti_latm_cnt(1) .gt. 0) call ti_exchange_vec(atm_cnt, vel, .false.)
        end if
        if (belly) call bellyf(atm_cnt, atm_igroup, vel)

        ! At this point in the code, the velocities lag the positions
        ! by half a timestep.  If we intend for the velocities to be drawn
        ! from a Maxwell distribution at the timepoint where the positions and
        ! velocities are synchronized, we have to correct these newly
        ! redrawn velocities by backing them up half a step using the
        ! current force.
        ! Note that this fix only works for Newtonian dynamics.

        if (.not. is_langevin) then
#if defined(MPI)
          if (usemidpoint) then
            do j=1,proc_num_atms
              wfac = 1.0d0/proc_atm_mass(j) * half_dtx
              proc_atm_vel(:, j) = proc_atm_vel(:, j) - proc_atm_frc(:, j) * wfac
            end do
          else ! usemidpoint
            do atm_lst_idx = 1, my_atm_cnt
              j = my_atm_lst(atm_lst_idx)
              wfac = atm_mass_inv(j) * half_dtx
              vel(:, j) = vel(:, j) - frc(:, j) * wfac
            end do
          end if ! usemidpoint
#else
          do j = 1, atm_cnt
            wfac = atm_mass_inv(j) * half_dtx
            vel(:, j) = vel(:, j) - frc(:, j) * wfac
          end do
#endif
        end if

#endif  /* ifdef CUDA */
      end if  ! (reset_velocities)
#ifdef GBTimer
      call get_wall_time(wall_s, wall_u)
      print *, "reset vel =", dble(wall_s) * 1000.0d0 + &
        dble(wall_u) / 1000.0d0 - strt_time_ms

      call get_wall_time(wall_s, wall_u)
      strt_time_ms = dble(wall_s) * 1000.0d0 + dble(wall_u) / 1000.0d0
#endif

      ! Do the velocity update:
#ifdef MPI
      if (vv==1) then
#ifdef CUDA
        call gpu_download_frc(frc)   !DG: check if the frc and vel download is necessary
        call gpu_download_vel(vel)
#endif
        call quench(frc,vel)     !added by DG for quench md for NEB.
#ifdef CUDA
        call gpu_upload_vel(vel)
#endif
      end if
#endif /* MPI */

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
        if (master) call corpac(nrx, frc, 1, mdfrc)
        !why do we do the update here? shouldn't it be post gpu_update?
        !I think we need to call this here because we're not doing the langevin setvel
        !stuff on the GPU but we are changing velocities.
      end if ! dump_frcs

      !we might need to do something with the velocities here
      !like do a vector exchange
      !    if (ti_mode .eq. 0) then

#ifdef GTI
      if (ti_mode .ne.0)  then 
       if (do_localheating) then
          if (mod(total_nstep, ntpr) .eq. 0) then 
                     
            if (total_nstep .le. gti_hsteps1) then
              write(mdout,'(a,f8.1,a)') 'TI-local heating at ', gti_tempi, 'K'
            else 
              temp=temp0 + (gti_tempi-temp0)*  &
                (gti_hsteps2*1.0-total_nstep*1.0)  &
               /(gti_hsteps2*1.0-gti_hsteps1*1.0)
              write(mdout,'(a,f8.1,a)') 'TI-local heating at ', temp, 'K'
              call gti_setup_localheating(temp)
            endif          

            if (total_nstep .gt. gti_hsteps2) then
              write(mdout,'(a)') 'Turn off TI-local heating '
              call gti_turnoff_localheating
              do_localheating=.false.
            endif
          
           end if
       endif
    endif

#endif

#if ( defined(_DEBUG) )
#ifdef CUDA
      call gpu_download_frc(frc)
#endif
      if (icfe.ne.0 .and. igamd.gt.0) then
      do i = 1, atm_cnt
            if ( ti_sc_lst(i) .gt. 0 ) then
                 write(*,'(a,2I10,3F15.5)') "GaMD runmd before gpu_update) nstep, atom, frc = ", nstep+1, i, frc(:,i)
            end if
      end do
    endif
#endif
    if (isgld .gt. 0) then
         ! update accumulators
        call gpu_sgenergy(si(si_pot_ene))
        ! SGLD time step
        call gpu_sgld_update(dt, temp0, gamma_ln,sggamma,com0sg,com1sg,com2sg)
    else
      if (ischeme == 1) then
        call gpu_update(dt, temp0, therm_par)
      else 
        call gpu_update(dt, temp0, gamma_ln)
      end if
    end if
#ifdef GTI
      if (icfe.ne.0 ) then 
        call gti_sync_vector(3, gti_syn_mass+1) ! sync crd & vel: (gti_syn_mass=0: weighted-combination; gti_syn_mass=1: copy over)
      endif
#endif
#else /* CUDA */
      if (.not. is_langevin .and. ischeme == 0) then

        ! ---Newtonian dynamics:
        ! Applying guiding force effect:
        if (isgld .gt. 0) then

          if (.not. usemidpoint) then
#ifdef MPI
            call sgmdw(atm_cnt, my_atm_cnt, dtx, si(si_pot_ene), mass, &
              atm_mass_inv, crd, frc, vel, my_atm_lst)
#else
            call sgmdw(natom, atm_cnt,  dtx, si(si_pot_ene), mass, &
              atm_mass_inv, crd, frc, vel)
#endif
          end if
        end if


        if (usemidpoint) then
#ifdef MPI
          DBG_ARRAYS_DUMP_3DBLE("runmd_vel1", proc_atm_to_full_list, proc_atm_vel, proc_num_atms)
          do j=1,proc_num_atms
            wfac = 1.0d0/proc_atm_mass(j) * dtx
            proc_atm_vel(:, j) = proc_atm_vel(:, j) + proc_atm_frc(:, j) * wfac
          end do
          DBG_ARRAYS_DUMP_3DBLE("runmd_vel2", proc_atm_to_full_list, proc_atm_vel,proc_num_atms)
#endif /*MPI*/
        else ! usemidpoint
#if defined(MPI) && !defined(CUDA)
#ifdef _OPENMP_
       !$omp parallel default(none) &
       !$omp& shared(my_atm_cnt,atm_mass_inv,vel,frc,dtx) &
       !$omp& private(wfac,j,atm_lst_idx)
       !$omp do
#endif /*_OPENMP_*/
          do atm_lst_idx = 1, my_atm_cnt
            j = my_atm_lst(atm_lst_idx)
#else
          do j = 1, atm_cnt
#endif /*  MPI */
            wfac = atm_mass_inv(j) * dtx
            vel(:, j) = vel(:, j) + frc(:, j) * wfac
          end do
#if defined(MPI) && !defined(CUDA)
#ifdef _OPENMP_
        !$omp end do
        !$omp end parallel
#endif /*_OPENMP_*/
#endif
        end if ! usemidpoint

      else if (isgld .gt. 0) then
        !  Using SGLD algorithm:
        if (.not. usemidpoint) then
#ifdef MPI
          call sgldw(atm_cnt,my_atm_cnt,  dtx, si(si_pot_ene), mass, &
            atm_mass_inv, crd, frc, vel, my_atm_lst)
#else
          call sgldw(natom,atm_cnt,  dtx,  si(si_pot_ene), mass, &
            atm_mass_inv, crd, frc, vel)
#endif
        end if ! not usemidpoint

      else if (ischeme == 1) then

        if (usemidpoint) then
#ifdef MPI
          DBG_ARRAYS_DUMP_3DBLE("runmd_vel1", proc_atm_to_full_list, proc_atm_vel, proc_num_atms)
          do j=1,proc_num_atms
            wfac = 1.0d0/proc_atm_mass(j) * dtx
            proc_atm_vel(:, j) = proc_atm_vel(:, j) + proc_atm_frc(:, j) * wfac
          end do
          DBG_ARRAYS_DUMP_3DBLE("runmd_vel2", proc_atm_to_full_list, proc_atm_vel,proc_num_atms)
#endif /*MPI*/
        else ! usemidpoint
#if defined(MPI) && !defined(CUDA)
#ifdef _OPENMP_
       !$omp parallel default(none) &
       !$omp& shared(my_atm_cnt,atm_mass_inv,vel,frc,dtx) &
       !$omp& private(wfac,j,atm_lst_idx)
       !$omp do
#endif /*_OPENMP_*/
          do atm_lst_idx = 1, my_atm_cnt
            j = my_atm_lst(atm_lst_idx)
#else
          do j = 1, atm_cnt
#endif /*  MPI */
            wfac = atm_mass_inv(j) * dtx
            vel(:, j) = vel(:, j) + frc(:, j) * wfac
          end do
#if defined(MPI) && !defined(CUDA)
#ifdef _OPENMP_
        !$omp end do
        !$omp end parallel
#endif /*_OPENMP_*/
#endif
        end if ! usemidpoint
        ! used by middle-scheme
        if (ntc .ne. 1) then
            if (usemidpoint) then
#ifdef MPI
                  call comm_3dimensions_3dbls(proc_atm_vel)
                call rattle(proc_atm_crd, proc_atm_vel)
                call rattle_fastwater(proc_atm_crd, proc_atm_vel)

                  DBG_ARRAYS_DUMP_CRD("runmd_crd2", proc_atm_to_full_list, proc_atm_crd,proc_num_atms, proc_atm_wrap)
#endif
            else
                call rattle(crd, vel)
                call rattle_fastwater(crd, vel)
            end if
        end if

      else  !  ntt .eq. 3, we are doing langevin dynamics

        ! Simple model for Langevin dynamics, basically taken from
        ! Loncharich, Brooks and Pastor, Biopolymers 32:523-535 (1992),
        ! Eq. 11.  (Note that the first term on the rhs of Eq. 11b
        ! should not be there.)

        if (usemidpoint) then
#ifdef MPI
          call langevin_setvel_midpoint(atm_cnt, proc_num_atms, proc_atm_vel, proc_atm_frc, &
            proc_atm_mass, dt, temp0, gamma_ln)
#endif /*MPI*/
        else ! usemidpoint
          call langevin_setvel(atm_cnt, vel, frc, mass, atm_mass_inv, &
            dt, temp0, gamma_ln)
        end if ! usemidpoint
        if (ti_latm_cnt(1) .gt. 0) then
          call ti_exchange_vec(atm_cnt, vel, .false.)
        end if

      end if      ! (is_langevin)
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

#ifndef CUDA
      !vlimit is not supported on GPU for speed reasons.

      if (use_vlimit) then

        ! We here provide code that is most efficient if vlimit is not exceeded.

        vlimit_exceeded = .false.

#if defined(MPI) && !defined(CUDA)
        if (usemidpoint) then
          do j = 1, proc_num_atms
            if (abs(proc_atm_vel(1, j)) .gt. vlimit .or. &
              abs(proc_atm_vel(2, j)) .gt. vlimit .or. &
              abs(proc_atm_vel(3, j)) .gt. vlimit) then
              vlimit_exceeded = .true.
              exit
            end if
          end do
        else !usemidpoint
#ifdef _OPENMP_
         !$omp parallel default(none) &
         !$omp& shared(vel,vlimit,vlimit_exceeded,my_atm_cnt) &
         !$omp& private(atm_lst_idx, j)
         !$omp do
#endif /*_OPENMP_*/
          do atm_lst_idx = 1, my_atm_cnt
            j = my_atm_lst(atm_lst_idx)
#else
          do j = 1, atm_cnt
#endif
            if (abs(vel(1, j)) .gt. vlimit .or. &
                abs(vel(2, j)) .gt. vlimit .or. &
                abs(vel(3, j)) .gt. vlimit) then
              vlimit_exceeded = .true.
#ifndef _OPENMP_
              exit
#endif
            end if
          end do
#ifdef _OPENMP_
        !$omp end do
        !$omp end parallel
#endif /*_OPENMP_*/

#ifdef MPI
        end if !usemidpoint
#endif
        if (vlimit_exceeded) then
          vmax = 0.d0
#if defined(MPI) && !defined(CUDA)
          if (usemidpoint) then
            do j = 1, proc_num_atms
              do i = 1,3
                vmax = max(vmax, abs(proc_atm_vel(i, j)))
                vel(i, j) = sign(min(abs(proc_atm_vel(i, j)), vlimit), proc_atm_vel(i, j))
              end do
            end do
          else !usemidpoint

            do atm_lst_idx = 1, my_atm_cnt
              j = my_atm_lst(atm_lst_idx)
#else
            do j = 1, atm_cnt
#endif
              do i = 1, 3
                vmax = max(vmax, abs(vel(i, j)))
                vel(i, j) = sign(min(abs(vel(i, j)), vlimit), vel(i, j))
              end do
            end do
            ! Only violations on the master node are actually reported
            ! to avoid both MPI communication and non-master writes.
            write(mdout, '(a,i6,a,f10.4)')  'vlimit exceeded for step ', &
                                        total_nstep, '; vmax = ', vmax
#ifdef MPI
          end if !usemidpoint
#endif
        end if
      end if

#endif /*ifndef CUDA*/
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

      if (ischeme == 1) then
          if (usemidpoint) then
#ifdef MPI
              do i = 1, proc_num_atms + proc_ghost_num_atms
              proc_atm_frc(:,i) = proc_atm_crd(:,i)
              proc_atm_crd(:,i) = proc_atm_crd(:,i) + proc_atm_vel(:,i) * dtx / 2.0
              end do

              ! used by middle-scheme
              call middle_langevin_thermostat_midpoint(atm_cnt, proc_num_atms, proc_atm_vel, &
                  proc_atm_mass, dt, temp0, therm_par)

              
              do i = 1, proc_num_atms + proc_ghost_num_atms
              !proc_atm_frc(:,i) = proc_atm_crd(:,i)
              proc_atm_crd(:,i) = proc_atm_crd(:,i) + proc_atm_vel(:,i) * dtx / 2.0
              end do

#endif /*MPI*/

          else
#ifdef MPI
#ifdef _OPENMP_
        !$omp parallel default(none) &
         !$omp& shared(vel,my_atm_cnt,frc,crd,dtx) &
         !$omp& private(atm_lst_idx, i)
         !$omp do
#endif /*_OPENMP_*/
                do atm_lst_idx = 1, my_atm_cnt
                  i = my_atm_lst(atm_lst_idx)
#else
                do i = 1, atm_cnt
#endif
                  frc(:, i) = crd(:, i)
                  crd(:, i) = crd(:, i) + vel(:, i) * dtx / 2.0
                end do
#ifdef _OPENMP_
       !$omp end do
       !$omp end parallel
#endif /*_OPENMP_*/
                
                call middle_langevin_thermostat(atm_cnt, vel, mass, atm_mass_inv, &
                    dt, temp0, therm_par)
#ifdef MPI
#ifdef _OPENMP_
        !$omp parallel default(none) &
         !$omp& shared(vel,my_atm_cnt,frc,crd,dtx) &
         !$omp& private(atm_lst_idx, i)
         !$omp do
#endif /*_OPENMP_*/
                do atm_lst_idx = 1, my_atm_cnt
                  i = my_atm_lst(atm_lst_idx)
#else
                do i = 1, atm_cnt
#endif
                  !frc(:, i) = crd(:, i)
                  crd(:, i) = crd(:, i) + vel(:, i) * dtx / 2.0
                end do
#ifdef _OPENMP_
       !$omp end do
       !$omp end parallel
#endif /*_OPENMP_*/
                
          end if
      else  ! non-middle-scheme
          if (usemidpoint) then
#ifdef MPI
            do i = 1, proc_num_atms + proc_ghost_num_atms
              proc_atm_frc(:,i) = proc_atm_crd(:,i)
              proc_atm_crd(:,i) = proc_atm_crd(:,i) + proc_atm_vel(:,i) * dtx
            end do
#endif /*MPI*/
          else ! usemidpoint
#ifdef MPI
#ifdef _OPENMP_
        !$omp parallel default(none) &
         !$omp& shared(vel,my_atm_cnt,frc,crd,dtx) &
         !$omp& private(atm_lst_idx, i)
         !$omp do
#endif /*_OPENMP_*/
            do atm_lst_idx = 1, my_atm_cnt
              i = my_atm_lst(atm_lst_idx)
#else
            do i = 1, atm_cnt
#endif
              frc(:, i) = crd(:, i)
              crd(:, i) = crd(:, i) + vel(:, i) * dtx
            end do
#ifdef _OPENMP_
       !$omp end do
       !$omp end parallel
#endif /*_OPENMP_*/
          end if ! usemidpoint
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
#endif /*MPI*/

      ! If shake is being used, update new positions to fix bond lengths:

      if (ntc .ne. 1) then

        if (ischeme .ne. 0) then
            if (usemidpoint) then
#ifdef MPI
                do i = 1, proc_num_atms
                    xold(:, i) = proc_atm_crd(:, i)
                end do
#endif
            else
                xold = crd;
            end if
        end if

        call update_time(runmd_time)
#ifdef CUDA
        call gpu_shake()
#ifdef GTI
        if(tishake .eq. 2) call gti_sync_vector(2, 2) ! Copy crds V0 to V1
#endif
#else
        if (usemidpoint) then
#ifdef MPI      
          ! We can't shake until all crds are updated
          call shake(proc_atm_frc, proc_atm_crd)
          call shake_fastwater(proc_atm_frc, proc_atm_crd)

          call update_time(shake_time)

          call comm_3dimensions_3dbls(proc_atm_crd, .true.)

          DBG_ARRAYS_DUMP_CRD("runmd_crd2", proc_atm_to_full_list, proc_atm_crd,proc_num_atms, proc_atm_wrap)

          call update_time(fcve_dist_time)
#endif /*MPI*/
        else ! usemidpoint
          call shake(frc, crd)
          call shake_fastwater(frc, crd)
          if(ti_mode .ne. 0 .and. ti_mode .ne. 1 .and. tishake .eq. 2) then
              call ti_copy_common_core(atm_cnt,crd)
          end if
        end if ! usemidpoint
#endif
        call update_time(shake_time)

        ! Must update extra point coordinates after the frame is moved:
        if (numextra .gt. 0 ) then
          if (frameon .ne. 0 .and. gbl_frame_cnt .ne. 0) &
#ifdef CUDA
          call gpu_local_to_global()
#else
          call local_to_global(crd, ep_frames, ep_lcl_crd, gbl_frame_cnt)
#endif
        end if

#ifdef GBTimer
        call get_wall_time(wall_s, wall_u)
        print *, "Shake wall time =", dble(wall_s) * 1000.0d0 + &
          dble(wall_u) / 1000.0d0 - strt_time_ms
#endif
        ! Need to synchronize coordinates for linearly scaled atoms after shake
        if (ti_mode .ne. 0) then
          if (ti_latm_cnt(1) .gt. 0) then
#ifdef CUDA
#ifdef GTI
            ! sync crd: (gti_syn_mass=0: weighted-combination; 
            ! gti_syn_mass=1: copy over, from V0 to V1
            ! gti_syn_mass=2: copy over, from V1 to V0
            ! gti_syn_mass=3: copy over, from SHAKE to non-SHAKE
            call gti_sync_vector(2, gti_syn_mass+1)
       
#else
            call gpu_ti_exchange_crd()
#endif

#else
            call ti_exchange_vec(atm_cnt, crd, .false.)
#endif
          end if
        end if

        ! Re-estimate velocities from differences in positions:
#ifdef GBTimer
        call get_wall_time(wall_s, wall_u)
        strt_time_ms = dble(wall_s) * 1000.0d0 + dble(wall_u) / 1000.0d0
#endif

#ifdef CUDA
        if (.not. reset_velocities) call gpu_recalculate_velocities(dtx_inv)


        if (ti_mode .ne. 0) then
          if (ti_latm_cnt(1) .gt. 0) then
#ifdef GTI
           ! sync vel: (gti_syn_mass=0: weighted-combination; gti_syn_mass=1: copy over)
            call gti_sync_vector(1, gti_syn_mass+1)
#else
            call gpu_ti_exchange_vel()
#endif
          end if
        end if

        ! rattle used by middle-scheme 
        if (ischeme .ne. 0) then
            call gpu_rattle(dt)
        end if
#else

        if (usemidpoint) then
#ifdef MPI
          DBG_ARRAYS_DUMP_3DBLE("runmd_vel3", proc_atm_to_full_list, proc_atm_vel,proc_num_atms)
          do j=1,proc_num_atms
            if (ischeme == 0) then
                proc_atm_vel(:,j)=(proc_atm_crd(:,j) - proc_atm_frc(:,j))*dtx_inv
            else
                ! used by middle-scheme
                proc_atm_vel(:,j) = proc_atm_vel(:,j) + (proc_atm_crd(:,j) - xold(:,j))*dtx_inv
            end if
          end do
          DBG_ARRAYS_DUMP_3DBLE("runmd_vel4", proc_atm_to_full_list, proc_atm_vel,proc_num_atms)
#endif /*MPI*/
        else ! usemidpoint
#ifdef MPI
#ifdef _OPENMP_
        !$omp parallel default(none) &
        !$omp& shared(vel,crd,dtx_inv,my_atm_cnt,frc) &
        !$omp& private(atm_lst_idx,j)
        !$omp do
#endif /*_OPENMP_*/
          do atm_lst_idx = 1, my_atm_cnt
            j = my_atm_lst(atm_lst_idx)
#else
          do j = 1, atm_cnt
#endif
            if (ischeme == 0) then
                vel(:, j) = (crd(:, j) - frc(:, j)) * dtx_inv
            else
                vel(:, j) = vel(:, j) + (crd(:, j) - xold(:, j)) * dtx_inv
            end if
          end do
#ifdef _OPENMP_
        !$omp end do
        !$omp end parallel
#endif /*_OPENMP_*/
        end if ! usemidpoint

        if (ischeme .ne. 0) then
            if (usemidpoint) then
#ifdef MPI
                  call comm_3dimensions_3dbls(proc_atm_vel)
                call rattle(proc_atm_crd, proc_atm_vel)
                call rattle_fastwater(proc_atm_crd, proc_atm_vel)

                  DBG_ARRAYS_DUMP_CRD("runmd_crd2", proc_atm_to_full_list, proc_atm_crd,proc_num_atms, proc_atm_wrap)
#endif
            else
                call rattle(crd, vel)
                call rattle_fastwater(crd, vel)
            end if

        end if


#endif /* ifdef CUDA */
      else !ntc .ne. 1
        ! Must update extra point coordinates after the frame is moved:
        if (numextra .gt. 0 ) then
          if (frameon .ne. 0 .and. gbl_frame_cnt .ne. 0) &
#ifdef CUDA
            call gpu_local_to_global()
#else
            call local_to_global(crd, ep_frames, ep_lcl_crd, gbl_frame_cnt)
#endif
        end if

      end if !ntc .ne. 1
#ifdef GBTimer
      call get_wall_time(wall_s, wall_u)
      print *, "Recalculate velocities =", dble(wall_s) * 1000.0d0 + &
         dble(wall_u) / 1000.0d0 - strt_time_ms
      call get_wall_time(wall_s, wall_u)
      strt_time_ms = dble(wall_s) * 1000.0d0 + dble(wall_u) / 1000.0d0
#endif
      ! Extra points velocities are known to be bogus at this point.  We zero
      ! them when they are output; ignore them otherwise.
      ! If running REMD, it is possible that ntpr > nstlim, so when determining
      ! whether GPU needs to calc. KE use total step count instead of irespa

#ifdef CUDA
#ifdef MPI
      if (ineb>0.and.(beadid==1.or.beadid==neb_nbead) ) then
        call gpu_clear_vel()
      end if
#endif
#endif

#ifndef CUDA
#ifdef MPI
      if (ineb>0.and.(beadid==1.or.beadid==neb_nbead) ) then
        do atm_lst_idx = 1, my_atm_cnt
          i = my_atm_lst(atm_lst_idx)
          crd(:, i)=frc(:, i)
          vel(:,i)=0.d0
        end do
      end if
#endif /* MPI */
#endif /* ndefCUDA */

      if (update_kin_energy_on_current_step) then
        
#ifdef CUDA
        if (update_kin_energy_on_current_step .or. (mod(irespa, ene_avg_sampling) .eq. 0) .or. &

            ((ntwe .gt. 0) .and. (mod(total_nstep, max(1,ntwe)) .eq. 0)) .or. &
#ifdef MPI
            (ntpr .gt. nstlim .and. mod(total_nstep, ntpr) .eq. 0) .or. &
#endif
            (mod(irespa, ntpr) .eq. 0)) then
            call gpu_calculate_kinetic_energy(c_ave, eke, ekph, ekpbs )
#ifdef GTI
            if (ntt.eq.1 .or. mod(total_nstep+1, ntpr) .eq. 0 ) then
              
              if (ti_mode .ne. 0 ) then
                  call gti_kinetic(c_ave)
                  call gti_update_kin_energy_from_gpu(ti_eke, ti_eke, .true. )
                  call gti_update_kin_energy_from_gpu(ti_ekpbs, ti_ekpbs, .false.)
                  call gti_update_kin_energy_from_gpu(ti_ekph, ti_ekph, .false. )
              end if ! ti_mode
            end if
#endif
        end if
#else

        ! Get the kinetic energy, either for printing or for Berendsen.
        ! The process is completed later under mpi, in order to avoid an
        ! extra all_reduce.

        eke = 0.d0
        ekph = 0.d0
        ekpbs = 0.d0

        if (.not. is_langevin) then
          if (usemidpoint) then
#ifdef MPI
            do j=1,proc_num_atms
              eke = eke + proc_atm_mass(j) * &
                ((proc_atm_vel(1,j) + proc_atm_last_vel(1,j))**2 + &
                (proc_atm_vel(2,j) + proc_atm_last_vel(2,j))**2 + &
                (proc_atm_vel(3,j) + proc_atm_last_vel(3,j))**2)

              ekpbs = ekpbs + proc_atm_mass(j) * &
                ((proc_atm_vel(1,j) * proc_atm_last_vel(1,j)) + &
                (proc_atm_vel(2,j) * proc_atm_last_vel(2,j)) + &
                (proc_atm_vel(3,j) * proc_atm_last_vel(3,j)))

              ekph  = ekph  + proc_atm_mass(j) * &
                ((proc_atm_vel(1,j) * proc_atm_vel(1,j)) + &
                (proc_atm_vel(2,j) * proc_atm_vel(2,j)) + &
                (proc_atm_vel(3,j) * proc_atm_vel(3,j)))
            end do
#endif /*MPI*/
          else ! usemidpoint
#ifdef MPI
#ifdef _OPENMP_
        !$omp parallel default(none) &
        !$omp& private(j,atm_lst_idx) &
        !$omp& shared(eke,my_atm_cnt,vel,last_vel,mass,ekpbs,ekph)
        !$omp do reduction(+: eke, ekpbs, ekph)
#endif /*_OPENMP_*/
            do atm_lst_idx = 1, my_atm_cnt
              j = my_atm_lst(atm_lst_idx)
#else
            do j = 1, atm_cnt
#endif
              eke = eke + mass(j) * &
                     ((vel(1,j) + last_vel(1,j))**2 + &
                      (vel(2,j) + last_vel(2,j))**2 + &
                      (vel(3,j) + last_vel(3,j))**2)

              ekpbs = ekpbs + mass(j) * &
                     ((vel(1,j) * last_vel(1,j)) + &
                      (vel(2,j) * last_vel(2,j)) + &
                      (vel(3,j) * last_vel(3,j)))

              ekph  = ekph  + mass(j) * &
                     ((vel(1,j) * vel(1,j)) + &
                      (vel(2,j) * vel(2,j)) + &
                      (vel(3,j) * vel(3,j)))
            end do
#ifdef _OPENMP_
      !$omp end do
      !$omp end parallel
#endif /*_OPENMP_*/
          end if ! usemidpoint

          eke = 0.25d0 * eke

        else

          if (usemidpoint) then
#ifdef MPI
            do j=1, proc_num_atms
              eke = eke + proc_atm_mass(j) * 0.25d0 * c_ave *&
                ((proc_atm_vel(1,j) + proc_atm_last_vel(1,j))**2 + &
                (proc_atm_vel(2,j) + proc_atm_last_vel(2,j))**2 + &
                (proc_atm_vel(3,j) + proc_atm_last_vel(3,j))**2)
            end do
#endif /*MPI*/
          else ! usemidpoint
#ifdef MPI
#ifdef _OPENMP_
            !$omp parallel default(none) &
            !$omp& private(j,atm_lst_idx) &
            !$omp& shared(eke,my_atm_cnt,vel,last_vel,c_ave,mass)
            !$omp do reduction(+: eke)
#endif /*_OPENMP_*/
            do atm_lst_idx = 1, my_atm_cnt
              j = my_atm_lst(atm_lst_idx)
#else
              do j = 1, atm_cnt
#endif
                eke = eke + mass(j) * 0.25d0 * c_ave *&
                  ((vel(1,j) + last_vel(1,j))**2 + &
                  (vel(2,j) + last_vel(2,j))**2 + &
                  (vel(3,j) + last_vel(3,j))**2)
              end do
#ifdef _OPENMP_
              !$omp end do
              !$omp end parallel
#endif /*_OPENMP_*/
            end if ! usemidpoint

          end if ! (is_langevin)

          eke = 0.5d0 * eke
          ekph = ekph * 0.5d0
          ekpbs = ekpbs * 0.5d0

#endif /* CUDA */

          ! NOTE - These ke terms are not yet summed under MPI!!!

        end if ! (update_kin_energy_on_current_step)
#ifdef GBTimer
        call get_wall_time(wall_s, wall_u)
        print *, "kinetic energy =", dble(wall_s) * 1000.0d0 + &
          dble(wall_u) / 1000.0d0 - strt_time_ms
#endif

        ! Scale coordinates if constant pressure run:

        ! ntp=1, isotropic pressure coupling
        if (ntp .eq. 1 .and. barostat .eq. 1) then

          rmu(1) = (1.d0 - dtcp * (pres0 - si(si_tot_press)))**one_third
          rmu(2) = rmu(1)
          rmu(3) = rmu(1)

          ! ntp=2, anisotropic pressure coupling
        else if (ntp .eq. 2 .and. barostat .eq. 1) then

          ! csurften>0, constant surface tension adjusts the tangential pressures
          ! See Zhang. J. Chem. Phys. 1995
          if (csurften .gt. 0) then

            if (csurften.eq.1) then  ! For surface tension in the x direction
              pres0y = pres0x - gamma_ten_int * ten_conv / pbc_box(1)
              pres0z = pres0y

            else if (csurften.eq.2) then !For surface tension in the y direction
              pres0x = pres0y - gamma_ten_int * ten_conv / pbc_box(2)
              pres0z = pres0x

              !else if (csurften.eq.3) then !For surface tension in the z direction
            else
              pres0x = pres0z - gamma_ten_int * ten_conv / pbc_box(3)
              pres0y = pres0x

            end if

            rmu(1) = (1.d0 - dtcp * (pres0x - press(1)))**one_third
            rmu(2) = (1.d0 - dtcp * (pres0y - press(2)))**one_third
            rmu(3) = (1.d0 - dtcp * (pres0z - press(3)))**one_third

          else  ! csurften = 0:

            rmu(1) = (1.d0 - dtcp * (pres0 - press(1)))**one_third
            rmu(2) = (1.d0 - dtcp * (pres0 - press(2)))**one_third
            rmu(3) = (1.d0 - dtcp * (pres0 - press(3)))**one_third

          end if

          ! ntp=3, semiisotropic pressure coupling
          ! (currently only for csurften>0, constant surface tension)

        else if (ntp .eq. 3 .and. barostat .eq. 1) then

          ! csurften>0, constant surface tension
          if (csurften .gt. 0) then

            if (csurften.eq.1) then    ! For surface tension in the x direction
              pres0y = pres0x - gamma_ten_int * ten_conv / pbc_box(1)
              pres0z = pres0y
              press_tan_ave = (press(2)+press(3))/2
              rmu(1) = (1.d0 - dtcp * (pres0x - press(1)))**one_third
              rmu(2) = (1.d0 - dtcp * (pres0y - press_tan_ave))**one_third
              rmu(3) = (1.d0 - dtcp * (pres0z - press_tan_ave))**one_third

            else if (csurften.eq.2) then ! For surface tension in the y direction
              pres0x = pres0y - gamma_ten_int * ten_conv / pbc_box(2)
              pres0z = pres0x
              press_tan_ave = (press(1)+press(3))/2
              rmu(1) = (1.d0 - dtcp * (pres0x - press_tan_ave))**one_third
              rmu(2) = (1.d0 - dtcp * (pres0y - press(2)))**one_third
              rmu(3) = (1.d0 - dtcp * (pres0z - press_tan_ave))**one_third

              !else if (csurften.eq.3) then ! For surface tension in the z direction
            else
              pres0x = pres0z - gamma_ten_int * ten_conv / pbc_box(3)
              pres0y = pres0x
              press_tan_ave = (press(1)+press(2))/2
              rmu(1) = (1.d0 - dtcp * (pres0x - press_tan_ave))**one_third
              rmu(2) = (1.d0 - dtcp * (pres0y - press_tan_ave))**one_third
              rmu(3) = (1.d0 - dtcp * (pres0z - press(3)))**one_third

            end if
          end if
          !Add semiisotropic pressure scaling in any direction with csurften=0 here
        else if (ntp .eq. 4) then
          rmu(:)=1.0
          if (mod(total_nstep+1, mcbarint) .eq. 0) then
             temp=target_n-(total_nstep*1.0)/(mcbarint*1.0)
             if (temp<1) temp=1.0
             rmu(1)=1.0+1.0/temp*(target_a-pbc_box(1))/pbc_box(1)
             rmu(2)=1.0+1.0/temp*(target_b-pbc_box(2))/pbc_box(2)
             rmu(3)=1.0+1.0/temp*(target_c-pbc_box(3))/pbc_box(3)
          endif
        end if

#ifdef GBTimer
        call get_wall_time(wall_s, wall_u)
        strt_time_ms = dble(wall_s) * 1000.0d0 + dble(wall_u) / 1000.0d0
#endif
        if (ntp .eq. 4) then

          call pressure_scale_pbc_data(rmu, vdw_cutoff + skinnb, verbose)

        else if (ntp .gt. 0 .and. barostat .eq. 1) then
           
          ! WARNING!!   This is not correct for non-orthogonal boxes if NTP > 1
          ! (i.e. non-isotropic scaling).  Currently general cell updates which
          ! allow cell angles to change are not implemented.  The virial tensor
          ! computed for ewald is the general Nose Klein; however the cell
          ! response needs a more general treatment.

          call pressure_scale_pbc_data(rmu, vdw_cutoff + skinnb, verbose)
#ifdef CUDA
          call gpu_pressure_scale(ucell, recip, uc_volume)
#else
          call set_cit_tbl_dims(pbc_box, vdw_cutoff + skinnb, cut_factor)

#ifdef MPI
          if (usemidpoint) then
            call pressure_scale_crds_midpoint(proc_atm_crd, proc_atm_mass, gbl_mol_mass_inv, &
              gbl_mol_com)
            reorient_flag = .true.
          else ! usemidpoint
            call pressure_scale_crds(crd, mass, gbl_mol_mass_inv, gbl_my_mol_lst, &
              gbl_mol_com)
          end if ! usemidpoint
#else
          call pressure_scale_crds(crd, mass, gbl_mol_mass_inv, gbl_mol_com)
#endif /* MPI */

          if (ntr .gt. 0 .and. natc .gt. 0) then
#ifdef MPI
            call pressure_scale_restraint_crds(atm_xc, mass, gbl_mol_mass_inv, &
              gbl_my_mol_lst)
#else
            call pressure_scale_restraint_crds(atm_xc, mass, gbl_mol_mass_inv)
#endif /* MPI */
          end if
#endif /* CUDA */
          DBG_ARRAYS_DUMP_CRD("runmd_crd3", proc_atm_to_full_list, proc_atm_crd,proc_num_atms, proc_atm_wrap)
        end if        ! ntp .gt. 0

#ifdef MPI
        if (usemidpoint) then
          call comm_3dimensions_3ints(proc_atm_wrap, .true.)
          call comm_3dimensions_3dbls(proc_atm_crd, .true.)
        end if
#endif

#ifdef GBTimer
        call get_wall_time(wall_s, wall_u)
        print *, "pressure scale =", dble(wall_s) * 1000.0d0 + &
          dble(wall_u) / 1000.0d0 - strt_time_ms
#endif

        if (using_pme_potential) then
#ifdef CUDA
          call gpu_skin_test()
#else
          ! Now we can do a skin check to see if we will have to rebuild the
          ! pairlist next time...

#ifdef MPI
          if (usemidpoint) then
            call proc_check_my_atom_movement(proc_atm_crd, proc_saved_atm_crd, &
              skinnb, ntp, proc_new_list)
          else ! usemidpoint
            call check_my_atom_movement(crd, gbl_atm_saved_crd, my_atm_lst, &
              skinnb, ntp, new_list)
          end if ! usemidpoint

#else /*serial version below*/
          call check_all_atom_movement(atm_cnt, crd, gbl_atm_saved_crd, skinnb, &
            ntp, new_list)
#endif /* MPI */
#endif /* CUDA */
        end if
        ! Nothing to do for Generalized Born...

#ifdef GBTimer
        call get_wall_time(wall_s, wall_u)
        strt_time_ms = dble(wall_s) * 1000.0d0 + dble(wall_u) / 1000.0d0
#endif

#if !defined(MPI)
        if (ti_mode .ne. 0) then
          !currently, ti on GPU only supports ntt .eq. 3
          if (update_kin_energy_on_current_step) then
            if (.not. is_langevin) then

#if defined(CUDA)

#if defined(GTI)
              if (ntt.eq.1 .or. mod(total_nstep+1, ntpr) .eq. 0) then
                  c_ave=1.0
                  call gti_kinetic(c_ave)
                  call gti_update_kin_energy_from_gpu(ti_eke, ti_eke, .true.)
              else  
                  call gti_update_kin_energy_from_gpu(ti_eke, ti_eke, .false.)
              end if
              call gti_update_kin_energy_from_gpu(ti_ekpbs, ti_ekpbs, .false.)
              call gti_update_kin_energy_from_gpu(ti_ekph, ti_ekph, .false. )
#else /* GTI */
              call gpu_ti_calculate_kinetic_energy(c_ave, ti_kin_ene(1,ti_eke), &
                ti_kin_ene(2,ti_eke), ti_kin_ene(1,ti_ekpbs), ti_kin_ene(2,ti_ekpbs), &
                ti_kin_ene(1,ti_ekph), ti_kin_ene(2,ti_ekph), ti_kin_ene(1,ti_sc_eke), &
                ti_kin_ene(2,ti_sc_eke))
#endif /* GTI  */

#else /* CUDA */
              call ti_calc_kin_ene(atm_cnt, mass, vel, last_vel, 1.d0, ti_eke)
              call ti_calc_kin_ene(atm_cnt, mass, vel, last_vel, 1.d0, ti_ekpbs)
              call ti_calc_kin_ene(atm_cnt, mass, vel, vel, 1.d0, ti_ekph)

#endif /* CUDA */

              buf_ctr = 4
              reduce_buf_out(buf_ctr + 1) = ti_kin_ene(1,ti_eke)
              reduce_buf_out(buf_ctr + 2) = ti_kin_ene(2,ti_eke)
              reduce_buf_out(buf_ctr + 3) = ti_kin_ene(1,ti_ekpbs)
              reduce_buf_out(buf_ctr + 4) = ti_kin_ene(2,ti_ekpbs)
              reduce_buf_out(buf_ctr + 5) = ti_kin_ene(1,ti_ekph)
              reduce_buf_out(buf_ctr + 6) = ti_kin_ene(2,ti_ekph)
              reduce_buf_out(buf_ctr + 7) = ti_kin_ene(1,ti_sc_eke)
              reduce_buf_out(buf_ctr + 8) = ti_kin_ene(2,ti_sc_eke)
              buf_ctr = buf_ctr + 8
            else
#if defined CUDA
#if defined(GTI)
              call gti_update_kin_energy_from_gpu(ti_eke, ti_eke, .false.)
#else
              call gpu_ti_calculate_kinetic_energy(c_ave, ti_kin_ene(1,ti_eke), &
                ti_kin_ene(2,ti_eke), ti_kin_ene(1,ti_ekpbs), ti_kin_ene(2,ti_ekpbs), &
                ti_kin_ene(1,ti_ekph), ti_kin_ene(2,ti_ekph), ti_kin_ene(1,ti_sc_eke), &
                ti_kin_ene(2,ti_sc_eke))
#endif /* GTI */

#else
              call ti_calc_kin_ene(atm_cnt, mass, vel, last_vel, c_ave, ti_eke)

#endif /* CUDA  */

              !save the ti_/ti_sc_ kinetic energies energies
              buf_ctr = 1
              reduce_buf_out(buf_ctr + 1) = ti_kin_ene(1,ti_eke)
              reduce_buf_out(buf_ctr + 2) = ti_kin_ene(2,ti_eke)
              reduce_buf_out(buf_ctr + 3) = ti_kin_ene(1,ti_sc_eke)
              reduce_buf_out(buf_ctr + 4) = ti_kin_ene(2,ti_sc_eke)
              buf_ctr = buf_ctr + 4
            end if
          end if
        end if
#endif

#ifdef GBTimer
        call get_wall_time(wall_s, wall_u)
        print *, "ti calc kinetic energy", dble(wall_s) * 1000.0d0 + &
          dble(wall_u) / 1000.0d0 - strt_time_ms
#endif

#if defined(MPI) && !defined(CUDA)
        if (usemidpoint) then
          if (proc_new_list) then
            proc_new_list_cnt = 1.d0
          else
            proc_new_list_cnt = 0.d0
          end if
        else ! usemidpoint
          if (new_list) then
            new_list_cnt = 1.d0
          else
            new_list_cnt = 0.d0
          end if
        end if ! usemidpoint
        ! Sum up the partial kinetic energies and complete the skin check.

        if (update_kin_energy_on_current_step) then
          reduce_buf_in(1) = eke
          reduce_buf_in(2) = ekph
          reduce_buf_in(3) = ekpbs
          if (usemidpoint) then
            reduce_buf_in(4) = proc_new_list_cnt
          else ! usemidpoint
            reduce_buf_in(4) = new_list_cnt
          end if ! usemidpoint
          buf_ctr = 4
        else
          if (usemidpoint) then
            reduce_buf_in(1) = proc_new_list_cnt
          else ! usemidpoint
            reduce_buf_in(1) = new_list_cnt
          end if ! usemidpoint
          buf_ctr = 1
        end if

        if (ti_mode .ne. 0) then
          if (update_kin_energy_on_current_step) then
            if (.not. is_langevin) then
              call ti_calc_kin_ene(atm_cnt, my_atm_cnt, my_atm_lst, mass, vel, &
                last_vel, 1.d0, ti_eke)
              call ti_calc_kin_ene(atm_cnt, my_atm_cnt, my_atm_lst, mass, vel, &
                last_vel, 1.d0, ti_ekpbs)
              call ti_calc_kin_ene(atm_cnt, my_atm_cnt, my_atm_lst, mass, vel, &
                vel, 1.d0, ti_ekph)
              reduce_buf_in(buf_ctr + 1) = ti_kin_ene(1,ti_eke)
              reduce_buf_in(buf_ctr + 2) = ti_kin_ene(2,ti_eke)
              reduce_buf_in(buf_ctr + 3) = ti_kin_ene(1,ti_ekpbs)
              reduce_buf_in(buf_ctr + 4) = ti_kin_ene(2,ti_ekpbs)
              reduce_buf_in(buf_ctr + 5) = ti_kin_ene(1,ti_ekph)
              reduce_buf_in(buf_ctr + 6) = ti_kin_ene(2,ti_ekph)
              reduce_buf_in(buf_ctr + 7) = ti_kin_ene(1,ti_sc_eke)
              reduce_buf_in(buf_ctr + 8) = ti_kin_ene(2,ti_sc_eke)
              buf_ctr = buf_ctr + 8
            else
              call ti_calc_kin_ene(atm_cnt, my_atm_cnt, my_atm_lst, mass, vel, &
                last_vel, c_ave, ti_eke)
              reduce_buf_in(buf_ctr + 1) = ti_kin_ene(1,ti_eke)
              reduce_buf_in(buf_ctr + 2) = ti_kin_ene(2,ti_eke)
              reduce_buf_in(buf_ctr + 3) = ti_kin_ene(1,ti_sc_eke)
              reduce_buf_in(buf_ctr + 4) = ti_kin_ene(2,ti_sc_eke)
              buf_ctr = buf_ctr + 4
            end if
          end if
        end if

        call update_time(runmd_time)
#ifdef COMM_TIME_TEST
        call start_test_timer(7, 'allreduce eke,list_cnt', 0)
#endif

        call mpi_allreduce(reduce_buf_in, reduce_buf_out, buf_ctr, &
          mpi_double_precision, &
          mpi_sum, pmemd_comm, err_code_mpi)

#ifdef COMM_TIME_TEST
        call stop_test_timer(7)
#endif

        call update_time(fcve_dist_time)

        if (update_kin_energy_on_current_step) then
          eke = reduce_buf_out(1)
          ekph = reduce_buf_out(2)
          ekpbs = reduce_buf_out(3)
          if (usemidpoint) then
            proc_new_list_cnt = reduce_buf_out(4)
          else ! usemidpoint
            new_list_cnt = reduce_buf_out(4)
          end if ! usemidpoint
        else
          if (usemidpoint) then
            proc_new_list_cnt = reduce_buf_out(1)
          else ! usemidpoint
            new_list_cnt = reduce_buf_out(1)
          end if ! usemidpoint
        end if

        ! Determine if any process saw an atom exceed the skin check.  We use the
        ! comparison to 0.5d0 to handle any rounding error issues; we use
        ! double precision for new_list_cnt in order to be able to piggyback the
        ! reduce.

        if (usemidpoint) then
          proc_new_list = (proc_new_list_cnt .ge. 0.5d0)
        else ! usemidpoint
          new_list = (new_list_cnt .ge. 0.5d0)
        end if ! usemidpoint

        if (.not.usemidpoint) then
          if (infe.ne.0.and.abmd_selection_freq.gt.0)then
             if(mod(abmd_mdstep,abmd_selection_freq).eq.0) &
             new_list = .true.
          endif
        end if

        if (using_pme_potential) then

          if (usemidpoint) then
            call check_new_list_limit(proc_new_list)
          else
            call check_new_list_limit(new_list)
          end if

        end if
        ! Nothing to do for Generalized Born...

        ! Now distribute the coordinates:

        call update_time(runmd_time)

        if(isgld>0.and.new_list)call sg_allgather(atm_cnt,using_pme_potential,using_gb_potential)

        if (usemidpoint) then
#ifdef EMIL
          if (emil_do_calc .gt. 0 .or. proc_new_list .or. nmropt .ne. 0) then
#else /*EMIL*/
          if (proc_new_list .or. nmropt .ne. 0) then
#endif /*EMIL*/
            all_crds_valid = .true.
          else
            ! It may be necessary for the master to have a complete copy of the
            ! coordinates for archiving and writing the restart file.

            collect_crds = .false.

            if (ntwx .gt. 0) then
              if (mod(total_nstep + 1, ntwx) .eq. 0) collect_crds = .true.
            end if

            if (mod(total_nstep + 1, ntwr) .eq. 0) collect_crds = .true.

            if (nstep + 1 .ge. nstlim) collect_crds = .true.

            all_crds_valid = .false.
          end if
        else ! usemidpoint
#ifdef EMIL
          if (using_gb_potential .or. emil_save_gb) then
#else
          if (using_gb_potential) then
#endif
            ! Generalized Born has it's own coordinate distribution scheme...

            call gb_mpi_allgathervec(atm_cnt, crd)

            all_crds_valid = .true.

! not sure if we need proc_new list or not
#ifdef EMIL
          else if (emil_do_calc .gt. 0 .or. new_list .or. &
                   nmropt .ne. 0 .or. infe .ne. 0) then
#else /*EMIL*/
          else if (new_list .or. nmropt .ne. 0 .or. infe .ne. 0) then
#endif /*EMIL*/
            ! This always gets done for a work redistribution cycle:
            call mpi_allgathervec(atm_cnt, crd)
            ! If this is a const pressure run and there are coordinate constraints,
            ! we will need to send the adjusted coordinate constraints around when
            ! we redistribute the work load.  Bummer.

            if (ntp .gt. 0) then
              if (ntr .gt. 0 .and. natc .gt. 0) then
                if (atm_redist_needed .or. fft_redist_needed) then
                  call mpi_allgathervec(atm_cnt, atm_xc)
                end if
              end if
            end if
            all_crds_valid = .true.

          else
            call distribute_crds(natom, crd)
            ! It may be necessary for the master to have a complete copy of the
            ! coordinates for archiving and writing the restart file.

            collect_crds = .false.

            if (ntwx .gt. 0) then
              if (mod(total_nstep + 1, ntwx) .eq. 0) collect_crds = .true.
            end if

            if (mod(total_nstep + 1, ntwr) .eq. 0) collect_crds = .true.

            if (nstep + 1 .ge. nstlim) collect_crds = .true.

            if (collect_crds) call mpi_gathervec(atm_cnt, crd)

            all_crds_valid = .false.
          end if
        end if ! usemidpoint

        call update_time(fcve_dist_time)

#endif /* MPI */
          ! Do for MPI and serial runs, see above ...
        if (ti_mode .ne. 0) then
#if !defined(CUDA) 
          !we don't need to subtract out any values for the GPU because
          ! we've already taken that into account in the kernel calculation
          if (update_kin_energy_on_current_step) then
            if (.not. is_langevin) then
              ! buf_ctr=12 Already calculated
              ti_kin_ene(1,ti_eke) = eke - reduce_buf_out(buf_ctr - 6)
              ti_kin_ene(2,ti_eke) = eke - reduce_buf_out(buf_ctr - 7)
              ti_kin_ene(1,ti_ekpbs) = ekpbs - reduce_buf_out(buf_ctr - 4)
              ti_kin_ene(2,ti_ekpbs) = ekpbs - reduce_buf_out(buf_ctr - 5)
              ti_kin_ene(1,ti_ekph) = ekph - reduce_buf_out(buf_ctr - 2)
              ti_kin_ene(2,ti_ekph) = ekph - reduce_buf_out(buf_ctr - 3)
              ti_kin_ene(1,ti_sc_eke) = reduce_buf_out(buf_ctr - 1)
              ti_kin_ene(2,ti_sc_eke) = reduce_buf_out(buf_ctr)
              ti_ene(1,si_kin_ene) = reduce_buf_out(buf_ctr - 7)
              ti_ene(2,si_kin_ene) = reduce_buf_out(buf_ctr - 6)
            else
              ! buf_ctr=5 Already calculated above
              ti_kin_ene(1,ti_eke) = eke - reduce_buf_out(buf_ctr - 2)
              ti_kin_ene(2,ti_eke) = eke - reduce_buf_out(buf_ctr - 3)
              ti_ene(1,si_kin_ene) = reduce_buf_out(buf_ctr - 3)
              ti_ene(2,si_kin_ene) = reduce_buf_out(buf_ctr - 2)
              ti_kin_ene(1,ti_sc_eke) = reduce_buf_out(buf_ctr - 1)
              ti_kin_ene(2,ti_sc_eke) = reduce_buf_out(buf_ctr)
            end if
          end if
#endif
#ifdef GTI
          if (update_kin_energy_on_current_step) then

            temp = ti_kin_ene(1,ti_eke)
            ti_kin_ene(1,ti_eke) = eke - ti_kin_ene(2,ti_eke)
            ti_kin_ene(2,ti_eke) = eke - temp

            if (.not. is_langevin) then
              temp = ti_kin_ene(1,ti_ekpbs)
              ti_kin_ene(1,ti_ekpbs) = ekpbs - ti_kin_ene(2,ti_ekpbs)
              ti_kin_ene(2,ti_ekpbs) = ekpbs - temp

              temp = ti_kin_ene(1,ti_ekph) 
              ti_kin_ene(1,ti_ekph) = ekph - ti_kin_ene(2,ti_ekph) 
              ti_kin_ene(2,ti_ekph) = ekph - temp
            end if
          end if
#endif

          !still need to store the total ti atoms kinetic energy
          !do we need to put ti_kin_ene(1/2, ti_sc_eke) somewhere?
          ti_ene_aug(1,ti_tot_kin_ene) = ti_kin_ene(1,ti_eke)
          ti_ene_aug(2,ti_tot_kin_ene) = ti_kin_ene(2,ti_eke)
          ti_ene_aug(1,ti_tot_temp) = ti_kin_ene(1,ti_eke) / ti_fac(1,1)
          ti_ene_aug(2,ti_tot_temp) = ti_kin_ene(2,ti_eke) / ti_fac(2,1)
        end if

        if (ntt .eq. 1) then

          !           --- following is from T.E. Cheatham, III and B.R. Brooks,
          !               Theor. Chem. Acc. 99:279, 1998.
          if (ti_mode .eq. 0) then
#ifdef CUDA
            scaltp = sqrt(1.d0+2.d0 * (dt / tautp)*(ekin0 - eke)/(ekmh + ekph))
            call gpu_scale_velocities(scaltp)
#else

            if (usemidpoint) then
#ifdef MPI
              if (proc_num_atms .gt. 0) &
                scaltp = sqrt(1.d0+2.d0 * (dt / tautp)*(ekin0 - eke)/(ekmh + ekph))
              do j = 1, proc_num_atms
                proc_atm_vel(:,j) = proc_atm_vel(:,j) * scaltp
              end do
#endif /*MPI*/
            else ! usemidpoint
#ifdef MPI
            if (my_atm_cnt .gt. 0) &
              scaltp = sqrt(1.d0+2.d0 * (dt / tautp)*(ekin0 - eke)/(ekmh + ekph))
            do atm_lst_idx = 1, my_atm_cnt
              j = my_atm_lst(atm_lst_idx)
#else
              scaltp = sqrt(1.d0+2.d0 * (dt / tautp)*(ekin0 - eke)/(ekmh + ekph))
              do j = 1, atm_cnt
#endif
                vel(:, j) = vel(:, j) * scaltp
              end do
            end if ! usemidpoint

#endif /* CUDA */
          else !TI

            ti_scaltp(1) = sqrt(1.d0 + 2.d0 * (dt / tautp) * (ti_ekin0(1) - &
              ti_kin_ene(1,ti_eke))/(ti_kin_ene(1,ti_ekmh) + &
              ti_kin_ene(1,ti_ekph)))
            ti_scaltp(2) = sqrt(1.d0 + 2.d0 * (dt / tautp) * (ti_ekin0(2) - &
              ti_kin_ene(2,ti_eke))/(ti_kin_ene(2,ti_ekmh) + &
              ti_kin_ene(2,ti_ekph)))
            scaltp = ti_scaltp(1) * ti_weights(1) + ti_scaltp(2) * ti_weights(2)

#ifdef CUDA
            call gpu_scale_velocities(scaltp)
#else

#ifdef MPI

            do atm_lst_idx = 1, my_atm_cnt
              j = my_atm_lst(atm_lst_idx)
#else
            do j = 1, atm_cnt
#endif
              if (ti_lst(3,j) .ne. 0) then
                vel(:, j) = vel(:, j) * scaltp
              else if (ti_lst(2,j) .ne. 0) then
                vel(:, j) = vel(:, j) * ti_scaltp(2)
              else
                vel(:, j) = vel(:, j) * ti_scaltp(1)
              end if
            end do
#endif /* CUDA */

          end if

        else if (update_bussi_target_kin_energy_on_current_step) then
          if (ti_mode == 0) then
            
            target_ekin  = resamplekin(eke, ekin0, int(rndf, 4), dt / tautp)
            scaltp = sqrt(target_ekin/eke)

#ifdef CUDA
            call gpu_scale_velocities(scaltp)
#else
            do j = 1, atm_cnt
              vel(:, j) = vel(:, j) * scaltp
            end do
#endif
          else
            ti_scaltp(1) = resamplekin(ti_kin_ene(1,ti_eke), ti_ekin0(1), int(rndf, 4), dt / tautp)
            ti_scaltp(2) = resamplekin(ti_kin_ene(2,ti_eke), ti_ekin0(2), int(rndf, 4), dt / tautp)
            ! Next 20 lines of code are copy-pasted from ntt=1 above
#ifdef CUDA
            call gpu_scale_velocities(scaltp)
#else

#ifdef MPI

            do atm_lst_idx = 1, my_atm_cnt
              j = my_atm_lst(atm_lst_idx)
#else
            do j = 1, atm_cnt
#endif
              if (ti_lst(3,j) /= 0) then
                vel(:, j) = vel(:, j) * scaltp
              else if (ti_lst(2,j) /= 0) then
                vel(:, j) = vel(:, j) * ti_scaltp(2)
              else
                vel(:, j) = vel(:, j) * ti_scaltp(1)
              end if
            end do
#endif /* CUDA */
          end if
        end if

        ! Pastor, Brooks, Szabo conserved quantity for harmonic oscillator:
        ! Eqn. 4.7b of Mol. Phys. 65:1409-1419, 1988:

        si(si_solvent_kin_ene) = ekpbs + si(si_pot_ene)

        si(si_solute_kin_ene) = eke
        si(si_kin_ene) = si(si_solute_kin_ene)

        if (ntt .eq. 1 .and. onstep) then
          ekmh = max(ekph, fac(1) * 10.d0)
        end if

        ! If velocities were reset, the KE is not accurate; fudge it here to keep
        ! the same total energy as on the previous step.  Note that this only
        ! affects printout and averages for Etot and KE -- it has no effect on the
        ! trajectory, or on any averages of potential energy terms.

        if (reset_velocities) si(si_kin_ene) = etot_save - si(si_pot_ene)

        ! --- total energy is sum of KE + PE:

        si(si_tot_ene) = si(si_kin_ene) + si(si_pot_ene)
        etot_save = si(si_tot_ene)

        ! Finish calculating KE for TI
        if (ti_mode .ne. 0) then
          ti_kin_ene(1,ti_solvent_kin_ene) = ti_kin_ene(1,ti_ekpbs)+si(si_pot_ene)
          ti_kin_ene(2,ti_solvent_kin_ene) = ti_kin_ene(2,ti_ekpbs)+si(si_pot_ene)

          ti_kin_ene(1,ti_solute_kin_ene) = ti_kin_ene(1,ti_eke)
          ti_kin_ene(2,ti_solute_kin_ene) = ti_kin_ene(2,ti_eke)

          if (ntt .eq. 1 .and. onstep) then
            ti_kin_ene(1,ti_ekmh) = max(ti_kin_ene(1,ti_ekph), ti_fac(1,1)*10.d0)
            ti_kin_ene(2,ti_ekmh) = max(ti_kin_ene(2,ti_ekph), ti_fac(2,1)*10.d0)
          end if

          if (reset_velocities) then
            ti_kin_ene(1,ti_eke) = ti_kin_ene(1,ti_etot_save) - si(si_pot_ene)
            ti_kin_ene(2,ti_eke) = ti_kin_ene(2,ti_etot_save) - si(si_pot_ene)
          end if

          ti_ene_aug(1,ti_tot_tot_ene) = ti_kin_ene(1,ti_eke) + si(si_pot_ene)
          ti_ene_aug(2,ti_tot_tot_ene) = ti_kin_ene(2,ti_eke) + si(si_pot_ene)

          ti_kin_ene(1,ti_etot_save) = ti_ene_aug(1,ti_tot_tot_ene)
          ti_kin_ene(2,ti_etot_save) = ti_ene_aug(2,ti_tot_tot_ene)

          ti_ene(1,si_kin_ene) = ti_kin_ene(1,ti_sc_eke)
          ti_ene(2,si_kin_ene) = ti_kin_ene(2,ti_sc_eke)
          ti_ene(1,si_tot_ene) = ti_ene(1,si_kin_ene) + ti_ene(1,si_pot_ene)
          ti_ene(2,si_tot_ene) = ti_ene(2,si_kin_ene) + ti_ene(2,si_pot_ene)
        end if


            ! For periodic PME, Zero COM velocity if requested; used for preventing ewald
            ! "block of ice flying thru space" phenomenon. We make this correction for pme
            ! explicitly before collecting velocities to the master...

        if (nscm .ne. 0) then
          if (mod(total_nstep + 1, nscm) .eq. 0) then
            if (ifbox .ne. 0) then
              if (.not. is_langevin .and. ischeme == 0) then !GPU TI only ntt=3 currently; unless configured with -gti
#ifdef CUDA
                call gpu_download_vel(vel)
#endif
                if (ti_mode .eq. 0) then
                  vcm(:) = 0.d0

                  if (usemidpoint) then
#ifdef MPI
                    do j=1,proc_num_atms
                      aamass = proc_atm_mass(j)
                      vcm(:) = vcm(:) + aamass * proc_atm_vel(:, j)
                    end do
#endif /*MPI*/
                  else ! usemidpoint
#if defined(MPI) && !defined(CUDA)
                    do atm_lst_idx = 1, my_atm_cnt
                    j = my_atm_lst(atm_lst_idx)
#else
                    do j = 1, atm_cnt
#endif
                      aamass = mass(j)
                      vcm(:) = vcm(:) + aamass * vel(:, j)
                    end do
                  end if ! usemidpoint

#if defined(MPI) && !defined(CUDA)
                  call update_time(runmd_time)
                  reduce_buf_in(1:3) = vcm(:)

#ifdef COMM_TIME_TEST
                  call start_test_timer(8, 'allreduce vcm', 0)
#endif

                  call mpi_allreduce(reduce_buf_in, reduce_buf_out, 3, &
                                     mpi_double_precision, &
                                     mpi_sum, pmemd_comm, err_code_mpi)

#ifdef COMM_TIME_TEST
                  call stop_test_timer(8)
#endif

                  vcm(:) = reduce_buf_out(1:3)
                  call update_time(fcve_dist_time)
#endif /* MPI */
                  vcm(:) = vcm(:) * tmassinv

                  if (master) then
                    velocity2 =vcm(1) * vcm(1) + vcm(2) * vcm(2) + vcm(3) * vcm(3)
                    write(mdout,'(a,f15.6,f9.2,a)') 'check COM velocity, temp: ',&
                                                sqrt(velocity2), 0.5d0 * tmass * &
                                                velocity2 / fac(1), '(Removed)'
                  end if

                  if (usemidpoint) then
#ifdef MPI
                    do j=1,proc_num_atms
                      proc_atm_vel(:,j)=proc_atm_vel(:,j) - vcm(:)
                    end do
#endif /*MPI*/
                  else ! usemidpoint
#if defined(MPI) && !defined(CUDA)
                    do atm_lst_idx = 1, my_atm_cnt
                      j = my_atm_lst(atm_lst_idx)
#else
                    do j = 1, atm_cnt
#endif
                      vel(:, j) = vel(:, j) - vcm(:)
                    end do
                  end if ! usemidpoint

                else ! if (ti_mode .ne. 0)
                  !for CUDA with ti, only langevin dynamics supported ; unless configured with -gti
                  !This is exactly what sander does with sc_mix_velocities
                  ti_vcm(:,:) = 0.d0

#if defined(MPI) && !defined(CUDA)
                  do atm_lst_idx = 1, my_atm_cnt
                    j = my_atm_lst(atm_lst_idx)
#else
                  do j = 1, atm_cnt
#endif
                    aamass = mass(j)
                    if (ti_lst(1,j) .ne. 0) then
                      ti_vcm(1,:) = ti_vcm(1,:) + aamass * vel(:, j)
                    else if (ti_lst(2,j) .ne. 0) then
                      ti_vcm(2,:) = ti_vcm(2,:) + aamass * vel(:, j)
                    else
                      ti_vcm(1,:) = ti_vcm(1,:) + aamass * vel(:, j)
                      ti_vcm(2,:) = ti_vcm(2,:) + aamass * vel(:, j)
                    end if
                  end do

#if defined(MPI) && !defined(CUDA)
                  call update_time(runmd_time)
                  reduce_buf_in(1:3) = ti_vcm(1,:)
                  reduce_buf_in(4:6) = ti_vcm(2,:)

#ifdef COMM_TIME_TEST
                  call start_test_timer(8, 'allreduce vcm', 0)
#endif

                  call mpi_allreduce(reduce_buf_in, reduce_buf_out, 2 * 3, &
                                   mpi_double_precision, &
                                   mpi_sum, pmemd_comm, err_code_mpi)

#ifdef COMM_TIME_TEST
                  call stop_test_timer(8)
#endif
                  ti_vcm(1,:) = reduce_buf_out(1:3)
                  ti_vcm(2,:) = reduce_buf_out(4:6)
                  call update_time(fcve_dist_time)
#endif /* MPI */
                  do j = 1, 2
                    ti_vcm(j,:) = ti_vcm(j,:) * ti_tmassinv(j)
                    if (master) then
                      velocity2 = ti_vcm(j,1) * ti_vcm(j,1) + &
                        ti_vcm(j,2) * ti_vcm(j,2) + &
                        ti_vcm(j,3) * ti_vcm(j,3)
                      write(mdout,'(a,i2,a,f15.6,f9.2,a)') &
                        'check COM velocity, temp for TI region ',j,':',&
                        sqrt(velocity2), 0.5d0 * ti_tmass(j) * &
                        velocity2 / ti_fac(j,1), '(Removed)'
                    end if
                  end do
                  vcm(:) = ti_vcm(1,:) * ti_weights(1) + ti_vcm(2,:) * ti_weights(2)
#if defined(MPI) && !defined(CUDA)
                  do atm_lst_idx = 1, my_atm_cnt
                    j = my_atm_lst(atm_lst_idx)
#else
                  do j = 1, atm_cnt
#endif
                    if (ti_lst(1,j) .ne. 0) then
                      vel(:, j) = vel(:, j) - ti_vcm(1,:)
                    else if (ti_lst(2,j) .ne. 0) then
                      vel(:, j) = vel(:, j) - ti_vcm(2,:)
                    else
                      vel(:, j) = vel(:, j) - vcm(:)
                    end if
                  end do
                end if  ! (ti_mode .eq. 0)
#ifdef CUDA
                call gpu_upload_vel(vel)
#endif
              end if  ! (.not. is_langevin)
            end if    ! (ifbox.ne.0)
          end if      ! (mod(total_nstep + 1, nscm) .eq. 0)
        end if        ! if nscm .ne. 0

#ifdef MPI
        ! Determine whether we want to save the on-step velocities.
        save_onstepvel = .false.
        if ((ionstepvelocities.ne.0).and.dump_vels) then
          save_onstepvel = .true.
          if (usemidpoint) then
            ! Midpoint method
            !if (all_vels_valid) then
            !  onstepvel(:,:) = (proc_atm_vel(:,:) + proc_atm_last_vel(:,:)) / 2.d0
            !else
              do j=1, proc_num_atms
                onstepvel(:,j) = (proc_atm_vel(:,j) + proc_atm_last_vel(:,j)) / 2.d0
              end do
            !end if
          else
            ! Not the midpoint method
            !if (all_vels_valid) then
            !  onstepvel(:,:) = (vel(:,:) + last_vel(:,:)) / 2.d0
            !else
              do atm_lst_idx = 1, my_atm_cnt
                j = my_atm_lst(atm_lst_idx)
                onstepvel(:,j) = (vel(:,j) + last_vel(:, j)) / 2.d0
              end do
            !end if
          end if
        endif ! Save on-step velocities
#endif /* MPI */

#if defined(MPI) && !defined(CUDA)
  ! It may be necessary for the master to have a complete copy of the velocities
  ! for archiving and writing the restart file. In any case, all processors must
  ! have a complete copy of velocities if atom redistribution is going to happen
  ! in the next call to force.

      ! For GB runs, always execute the "else"

        call update_time(runmd_time)
        ! not sure if we need proc_new list or not
        if (usemidpoint) then
          if (using_pme_potential .and. proc_new_list) then ! .and. &
            all_vels_valid = .true.
          else
            collect_vels = .false.

            if (ntwv .gt. 0) then
              if (mod(total_nstep + 1, ntwv) .eq. 0) collect_vels = .true.
            else if (collect_crds .and. ntwv .lt. 0) then
              collect_vels = .true.
            end if

            if (mod(total_nstep + 1, ntwr) .eq. 0) collect_vels = .true.

            if (nstep + 1 .ge. nstlim) collect_vels = .true.

            all_vels_valid = .false.
          end if
        else ! usemidpoint
          if (using_pme_potential .and. new_list .and. &
            (atm_redist_needed .or. fft_redist_needed)) then
            if (numextra .gt. 0 .and. frameon .ne. 0) &
              call zero_extra_pnts_vec(vel, ep_frames, gbl_frame_cnt)

            call mpi_allgathervec(atm_cnt, vel)
            if (save_onstepvel) call mpi_allgathervec(atm_cnt, onstepvel)
            all_vels_valid = .true.

          else

            collect_vels = .false.

            if (ntwv .gt. 0) then
              if (mod(total_nstep + 1, ntwv) .eq. 0) collect_vels = .true.
            else if (collect_crds .and. ntwv .lt. 0) then
              collect_vels = .true.
            end if

            if (mod(total_nstep + 1, ntwr) .eq. 0) collect_vels = .true.

            if (nstep + 1 .ge. nstlim) collect_vels = .true.

            if (collect_vels) then
              if (numextra .gt. 0 .and. frameon .ne. 0) &
                call zero_extra_pnts_vec(vel, ep_frames, gbl_frame_cnt)

              if (using_pme_potential) then
                call mpi_gathervec(atm_cnt, vel)
                if (save_onstepvel) call mpi_allgathervec(atm_cnt, onstepvel)
              else if (using_gb_potential) then
                call gb_mpi_gathervec(atm_cnt, vel)
                if (save_onstepvel) call gb_mpi_gathervec(atm_cnt, onstepvel)
              end if
            end if

            all_vels_valid = .false.

          end if
        end if ! usemidpoint

        ! Only need to collect forces for a force dump

        call update_time(fcve_dist_time)
#endif /* MPI */

        ! Zero COM velocity if requested; here we are doing the adjustment for a
        ! nonperiodic system, not pme.

        if (nscm .ne. 0) then
          if (mod(total_nstep + 1, nscm) .eq. 0) then
            if (ifbox .eq. 0) then
#if defined(MPI) && !defined(CUDA)
                ! WARNING - currently only GB code has ifbox .eq. 0, and currently
              ! all coordinates are always valid for GB.  We include the conditional
              ! below more-or-less as maintenance insurance...
              if (.not. usemidpoint) then
                if (.not. all_crds_valid) then ! Currently always false...
                  call update_time(runmd_time)
                  call gb_mpi_allgathervec(atm_cnt, crd)
                  all_crds_valid = .true.
                  call update_time(fcve_dist_time)
                end if
              end if
#endif /* MPI */
              if (.not. is_langevin .and. ischeme == 0) then !don't do this for GPU TI right now; unless configured with -gti
#ifdef CUDA
                call gpu_download_crd(crd)
                call gpu_download_vel(vel)
#endif

#if defined(MPI) && !defined(CUDA)
                  ! WARNING - currently GB code never updates all velocities unless
                !           forced by this scenario...
                if (.not. usemidpoint) then
                  if (.not. all_vels_valid) then ! Currently always true...

                    ! The following conditional is currently never .true., but we
                    ! include it if there is ever extra points support under GB, or
                    ! if other types of nonperiodic simulations with extra points
                    ! are eventually supported...

                    if (numextra .gt. 0 .and. frameon .ne. 0) &
                      call zero_extra_pnts_vec(vel, ep_frames, gbl_frame_cnt)

                    call update_time(runmd_time)
                    call gb_mpi_allgathervec(atm_cnt, vel)
                    all_vels_valid = .true.
                    call update_time(fcve_dist_time)
                  end if
                end if
#endif /* MPI */
                ! Nonperiodic simulation.  Remove both translation and rotation.
                ! Back the coords up a half step so they correspond to the velocities,
                ! temporarily storing them in frc(:,:).
                if (usemidpoint) then
#ifdef MPI
                  proc_atm_frc(:,:) = proc_atm_crd(:,:) - proc_atm_vel(:,:)*half_dtx
#endif /*MPI*/
                else ! usemidpoint
                  frc(:,:) = crd(:,:) - vel(:,:) * half_dtx
                end if ! usemidpoint
                ! Now compute COM motion and remove it; then recompute (sander
                ! compatibility...).

                ! NOTE - if mass can change, that has to be taken into account for
                !        tmass, tmassinv (say for TI).

                if (ti_mode .ne. 1) then
                  call cenmas(atm_cnt, frc, vel, tmass, tmassinv, mass, xcm, &
                    vcm, ocm)
                  call stopcm(atm_cnt, frc, vel, xcm, vcm, ocm)
                  call cenmas(atm_cnt, frc, vel, tmass, tmassinv, mass, xcm, &
                    vcm, ocm)
                else
                  call ti_cenmas(atm_cnt, frc, vel, tmass, tmassinv, mass, xcm, &
                    vcm, ocm)
                  call stopcm(atm_cnt, frc, vel, xcm, vcm, ocm)
                  call ti_cenmas(atm_cnt, frc, vel, tmass, tmassinv, mass, xcm, &
                    vcm, ocm)
                end if
#ifdef CUDA
                call gpu_upload_vel(vel)
#endif
              else  !  if (is_langevin):


#ifdef CUDA
                !this literally just adjusts the frame of reference. It doesn't actually change
                !coordinates
                call gpu_recenter_molecule()
#else
                ! Get current center of the system:
                call get_position(atm_cnt, crd, vcm(1), vcm(2), vcm(3), sys_range)

                ! Recenter system to the original center:
                ! All threads call so we do not need to send adjust coordinates or
                ! restraint coordinates around.

                call re_position(atm_cnt, ntr, crd, atm_xc, vcm(1), vcm(2), vcm(3), &
                  sys_x, sys_y, sys_z, sys_range, 0)

#endif
              end if  ! (is_langevin)
            end if    ! (ifbox .eq. 0)
          end if      ! (mod(total_nstep + 1, nscm) .eq. 0)
        end if        ! if nscm .ne. 0

        ! Zero out any non-moving velocities if a belly is active:

        if (belly) call bellyf(atm_cnt, atm_igroup, vel)

  ! Save old, or last velocities...
#ifndef CUDA
#ifdef MPI
        ! Save the plus-half-step-velocities
        if (usemidpoint) then
          if (all_vels_valid) then
            proc_atm_last_vel(:,:) = proc_atm_vel(:,:)
          else
            do j=1, proc_num_atms
              proc_atm_last_vel(:,j) = proc_atm_vel(:,j)
            end do
          end if
        else ! usemidpoint
          if (all_vels_valid) then
            last_vel(:,:) = vel(:,:)
          else
            do atm_lst_idx = 1, my_atm_cnt
              j = my_atm_lst(atm_lst_idx)
              last_vel(:, j) = vel(:, j)
            end do
          end if
        end if ! usemidpoint
#else
        ! Determine whether we want to save the on-step velocities.
        if ((ionstepvelocities.ne.0).and.dump_vels) then
          onstepvel(:,:) = (vel(:,:) + last_vel(:,:)) / 2.d0
        endif
        ! Save the plus-half-step velocities
        last_vel(:,:) = vel(:,:)
#endif /* MPI */
#endif /* CUDA */

  ! Update the step counter and the integration time:
        nstep = nstep + 1
        total_nstep = total_nstep + 1
        t = t + dt

#ifdef MPI
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

              !          if (icnstph .gt. 0 .or. (icnste .gt. 0 .and. cpein_specified)) call cnstph_write_restart
            end if !master
          end if

          ! Coordinate archive:

          if (ntwx .gt. 0) then

            if (mod(total_nstep, ntwx) .eq. 0) then
              if (.not. do_output_mpi) then
                if (ionstepvelocities.eq.0) then
                  call mpi_output_transfer(crd, vel)
                else
                  call mpi_output_transfer(crd, onstepvel)
                endif
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
                  if (ionstepvelocities.eq.0) then
                    call corpac(nrx, vel, 1, mdvel)
                  else
                    call corpac(nrx, onstepvel, 1, mdvel)
                  endif
                end if
              end if !master

              ! For ioutfm == 1, coordinate archive may also contain vels

              if (iamd.gt.0)then

                !Flush amdlog file containing the list of weights and energies at each step
                !            if (master) call write_amd_weights(ntwx,total_nstep)

              end if

              if (igamd.gt.0)then

                !Flush gamdlog file containing the list of weights and energies at each step
                !            if (master) call write_gamd_weights(ntwx,total_nstep)

              end if

              if (scaledMD.gt.0)then

                !Flush scaledMDlog file containing the list of lambda and total energy at each step
                !            if (master) call write_scaledMD_log(ntwx,total_nstep)

              end if

            end if

          end if

          ! Velocity archive:

          if (ntwv .gt. 0) then

            if (mod(total_nstep, ntwv) .eq. 0) then
              if (.not. do_output_mpi) then
                if (ionstepvelocities.eq.0) then
                  call mpi_output_transfer(crd, vel)
                else
                  call mpi_output_transfer(crd, onstepvel)
                endif
                do_output_mpi = .true.
              end if
              if (master) then

                if (ionstepvelocities.eq.0) then
                  call corpac(nrx, vel, 1, mdvel)
                else
                  call corpac(nrx, onstepvel, 1, mdvel)
                endif
              end if
            end if

          end if

          ! Flush binary netCDF file(s) if necessary...

          if (ioutfm .eq. 1) then

            if (ntwx .gt. 0) then
              if (mod(total_nstep, ntwx) .eq. 0) then
                if (.not. do_output_mpi) then
                  call mpi_output_transfer(crd, vel) ! DRR - I dont think this is necessary
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
                  call mpi_output_transfer(crd, vel) ! DRR - I dont think this is necessary
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

              if (icnstph .gt. 0 .or. (icnste .gt. 0 .and. cpein_specified)) call cnstph_write_restart(total_nstep)

              if (icnste .gt. 0 .and. .not. cpein_specified) call cnste_write_restart(total_nstep)

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
                  if (ionstepvelocities.eq.0) then
                    call corpac(nrx, vel, 1, mdvel)
                  else
                    call corpac(nrx, onstepvel, 1, mdvel)
                  endif
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

                if (ionstepvelocities.eq.0) then
                  call corpac(nrx, vel, 1, mdvel)
                else
                  call corpac(nrx, onstepvel, 1, mdvel)
                endif
              end if

            end if

            ! Structure factor archive
#ifndef MPI
#ifndef NOXRAY
            if (ntwsf .gt. 0) then
               if (mod(total_nstep, ntwsf) .eq. 0 .or. &
                   (nstlim .eq. 0 .and. total_nstep .eq. 0 .and. ntwsf .gt. 0)) then
                if (ioutfm .eq. 0) then
                  call formatted_strfac(xray_get_f_calc(), 1, num_hkl, sfout)
                else if (ioutfm .eq. 1) then
                  call write_binary_strfac(xray_get_f_calc(), num_hkl, sfout)
                end if
              end if
            end if
#endif
#endif            
          end if ! Exit master for velocity/crd printing

        end if ! usemidpoint

        ! Update the igamd flag upon new or restarted GaMD simulation
        if (((igamd.gt.0).or.(igamd0.gt.0)) .and. ntcmd.gt.0 .and. total_nstep.eq.ntcmd) then
              igamd = igamd0
#ifdef CUDA
              call gpu_igamd_update(igamd)
#endif
        end if

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
                  if (igamd.gt.0 .or. igamd0.gt.0) &
                    call ti_others_print_ene(ti_region, ti_others_ene(ti_region,:))
                end if
              end do
#ifdef GTI
              call gti_output_ti_result(nstep*t)    
#endif
            end if
            if (.not. usemidpoint) then
              if (nmropt .ne. 0) call nmrptx(mdout)
              if (infe .ne. 0) call nfe_prt(mdout)
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

                  if (isgld>0)then
                    call sgld_avg(.false.)
                  endif
                end if
                write(mdout, 540) ntave / nrespa
                call prntmd(total_nstep, total_nstlim, t, sit_tmp, fac, 0, &
                  .false., mdloop)
                if (.not. usemidpoint) then
                  if (isgld>0)then
                    call sgld_fluc(.false.)
                  endif
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
                    if (igamd.gt.0 .or. igamd0.gt.0) &
                      call ti_others_print_ene(ti_region, ti_others_ene_tmp(ti_region,:))
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
                    if (igamd.gt.0 .or. igamd0.gt.0) &
                      call ti_others_print_ene(ti_region, ti_others_ene_tmp2(ti_region,:))
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
                if ((igamd.eq.4 .or. igamd.eq.5) .or. (igamd0.eq.4 .or. igamd0.eq.5)) then
                VavgDt = sit_tmp(si_dihedral_ene)
                sigmaVDt = sit2_tmp(si_dihedral_ene)
                  VavgPt = sit_tmp(si_elect_ene) + sit_tmp(si_vdw_ene) + &
                    sit_tmp(si_vdw_14_ene) + sit_tmp(si_elect_14_ene) + sit_tmp(si_efield_ene)
                  sigmaVPt = sqrt(sit2_tmp(si_elect_ene)**2 + sit2_tmp(si_vdw_ene)**2 + &
                          sit2_tmp(si_vdw_14_ene)**2 + sit2_tmp(si_elect_14_ene)**2 + sit2_tmp(si_efield_ene)**2)
                else if ((igamd.eq.6 .or. igamd.eq.7) .or. (igamd0.eq.6 .or. igamd0.eq.7)) then
                  VavgPt = ti_ene_tmp(1,si_bond_ene) + ti_ene_tmp(1,si_angle_ene) + ti_ene_tmp(1,si_dihedral_ene) + &
                           ti_ene_tmp(1,si_vdw_ene) + ti_ene_tmp(1,si_elect_ene) + &
                           ti_others_ene_tmp(1,si_vdw_ene) + ti_others_ene_tmp(1,si_elect_ene)
                  sigmaVPt = sqrt(ti_ene_tmp2(1,si_bond_ene)**2 + ti_ene_tmp2(1,si_angle_ene)**2 + &
                             ti_ene_tmp2(1,si_dihedral_ene)**2 + &
                             ti_ene_tmp2(1,si_vdw_ene)**2 + ti_ene_tmp2(1,si_elect_ene)**2 + &
                             ti_others_ene_tmp2(1,si_vdw_ene)**2 + ti_others_ene_tmp2(1,si_elect_ene)**2)
                  VavgDt = sit_tmp(si_pot_ene) + ti_ene_tmp(1,si_vdw_14_ene) + ti_ene_tmp(1,si_elect_14_ene)
                  sigmaVDt = sqrt(sit2_tmp(si_pot_ene)**2 + &
                          ti_ene_tmp2(1, si_vdw_14_ene)**2 + ti_ene_tmp2(1, si_elect_14_ene)**2)
                else if ((igamd.eq.8 .or. igamd.eq.9) .or. (igamd0.eq.8 .or. igamd0.eq.9)) then
                  VavgPt = ti_ene_tmp(1,si_bond_ene) + ti_ene_tmp(1,si_angle_ene) + ti_ene_tmp(1,si_dihedral_ene)
                  sigmaVPt = sqrt(ti_ene_tmp2(1,si_bond_ene)**2 + ti_ene_tmp2(1,si_angle_ene)**2 + &
                          ti_ene_tmp2(1,si_dihedral_ene)**2 )
                  VavgDt = sit_tmp(si_pot_ene) + ti_ene_tmp(1,si_vdw_14_ene) + ti_ene_tmp(1,si_elect_14_ene) + &
                           ti_ene_tmp(1,si_vdw_ene) + ti_ene_tmp(1,si_elect_ene) + &
                           ti_others_ene_tmp(1,si_vdw_ene) + ti_others_ene_tmp(1,si_elect_ene)
                  sigmaVDt = sqrt(sit2_tmp(si_pot_ene)**2 + &
                          ti_ene_tmp2(1,si_vdw_14_ene)**2 + ti_ene_tmp2(1,si_elect_14_ene)**2 + &
                          ti_ene_tmp2(1,si_vdw_ene)**2 + ti_ene_tmp2(1,si_elect_ene)**2 + &
                          ti_others_ene_tmp2(1,si_vdw_ene)**2 + ti_others_ene_tmp2(1,si_elect_ene)**2)
                else if ((igamd.eq.10 .or. igamd.eq.11) .or. (igamd0.eq.10 .or. igamd0.eq.11)) then
                  VavgPt = ti_ene_tmp(1,si_vdw_ene) + ti_ene_tmp(1,si_elect_ene) + &
                           ti_others_ene_tmp(1,si_vdw_ene) + ti_others_ene_tmp(1,si_elect_ene)
                  sigmaVPt = sqrt(ti_ene_tmp2(1,si_vdw_ene)**2 + ti_ene_tmp2(1,si_elect_ene)**2 + &
                             ti_others_ene_tmp2(1,si_vdw_ene)**2 + ti_others_ene_tmp2(1,si_elect_ene)**2)
                  VavgDt = sit_tmp(si_pot_ene) + &
                           ti_ene_tmp(1,si_bond_ene) + ti_ene_tmp(1,si_angle_ene) + ti_ene_tmp(1,si_dihedral_ene) + &
                           ti_ene_tmp(1,si_elect_14_ene) + ti_ene_tmp(1,si_vdw_14_ene)
                  sigmaVDt = sqrt(sit2_tmp(si_pot_ene)**2 + &
                             ti_ene_tmp2(1,si_bond_ene)**2 + ti_ene_tmp2(1,si_angle_ene)**2 + &
                             ti_ene_tmp2(1,si_dihedral_ene)**2 + &
                             ti_ene_tmp2(1,si_elect_14_ene)**2 + ti_ene_tmp2(1,si_vdw_14_ene)**2)
                else if ((igamd.eq.100) .or. (igamd0.eq.100)) then
                  VavgDt = sit_tmp(si_pot_ene) + &
                           ti_ene_tmp(1,si_bond_ene) + ti_ene_tmp(1,si_angle_ene) + ti_ene_tmp(1,si_dihedral_ene) + &
                           ti_ene_tmp(1,si_elect_14_ene) + ti_ene_tmp(1,si_vdw_14_ene)
                  sigmaVDt = sqrt(sit2_tmp(si_pot_ene)**2 + &
                             ti_ene_tmp2(1,si_bond_ene)**2 + ti_ene_tmp2(1,si_angle_ene)**2 + &
                             ti_ene_tmp2(1,si_dihedral_ene)**2 + &
                             ti_ene_tmp2(1,si_elect_14_ene)**2 + ti_ene_tmp2(1,si_vdw_14_ene)**2)
                  VavgPt = ti_ene_tmp(1,si_vdw_ene) + ti_ene_tmp(1,si_elect_ene) + &
                           ti_others_ene_tmp(1,si_vdw_ene) + ti_others_ene_tmp(1,si_elect_ene)
                  sigmaVPt = sqrt(ti_ene_tmp2(1,si_vdw_ene)**2 + ti_ene_tmp2(1,si_elect_ene)**2 + &
                           ti_others_ene_tmp2(1,si_vdw_ene)**2 + ti_others_ene_tmp2(1,si_elect_ene)**2)
                else if ((igamd.eq.12 .or. igamd.eq.13) .or. (igamd0.eq.12 .or. igamd0.eq.13)) then
                  VavgDt = sit_tmp(si_gamd_ppi)
                  sigmaVDt = sit2_tmp(si_gamd_ppi)
                  VavgPt = sit_tmp(si_dvdl)
                  sigmaVPt =  sit2_tmp(si_dvdl)
                else if ((igamd.eq.14 .or. igamd.eq.15) .or. (igamd0.eq.14 .or. igamd0.eq.15)) then
                  VavgDt = sit_tmp(si_gamd_ppi)
                  sigmaVDt = sit2_tmp(si_gamd_ppi)
                  VavgPt = sit_tmp(si_dvdl)
                  sigmaVPt =  sit2_tmp(si_dvdl)
                else if ((igamd.eq.16 .or. igamd.eq.17) .or. (igamd0.eq.16 .or. igamd0.eq.17)) then
                  VavgDt = sit_tmp(si_gamd_ppi)
                  sigmaVDt = sit2_tmp(si_gamd_ppi)
                  VavgPt = sit_tmp(si_dvdl)
                  sigmaVPt =  sit2_tmp(si_dvdl)
                else if (igamd.eq.18  .or. igamd0.eq.18) then
                  VavgDt = sit_tmp(si_gamd_ppi)
                  sigmaVDt = sit2_tmp(si_gamd_ppi)
                  VavgPt = sit_tmp(si_dvdl)
                  sigmaVPt =  sit2_tmp(si_dvdl)
                else if (igamd.eq.19 .or. igamd0.eq.19) then
                  VavgDt = sit_tmp(si_dihedral_ene)
                  sigmaVDt =sit2_tmp(si_dihedral_ene)
                  VavgPt = sit_tmp(si_dvdl)
                  sigmaVPt =  sit2_tmp(si_dvdl)
                else if ((igamd.ge.20 .and. igamd.le.28).or.(igamd0.ge.20 .and. igamd0.le.28)) then
                   VavgDt = sit_tmp(si_gamd_ppi)
                   sigmaVDt = sit2_tmp(si_gamd_ppi)
                   VavgPt = sit_tmp(si_dvdl)
                   sigmaVPt =  sit2_tmp(si_dvdl)
                   VavgBt =sit_tmp(si_gamd_ppi2)
                   sigmaVBt = sit2_tmp(si_gamd_ppi2)
                else if (igamd.eq.110 .or. igamd0.eq.110) then
                  VavgDt = sit_tmp(si_dihedral_ene)
                  sigmaVDt = sit2_tmp(si_dihedral_ene)
                  VavgPt = sit_tmp(si_dvdl)
                  sigmaVPt =  sit2_tmp(si_dvdl)
                else if (igamd.eq.111 .or. igamd0.eq.111) then
                  VavgDt = sit_tmp(si_dihedral_ene)
                  sigmaVDt = sit2_tmp(si_dihedral_ene)
                  VavgPt = sit_tmp(si_dvdl)
                  sigmaVPt =  sit2_tmp(si_dvdl)
                else if (igamd.eq.112  .or. igamd0.eq.112) then
                  VavgDt = sit_tmp(si_dihedral_ene)
                  sigmaVDt = sit2_tmp(si_dihedral_ene)
                  VavgPt = sit_tmp(si_dvdl)
                  sigmaVPt =  sit2_tmp(si_dvdl)
                 else if (igamd.eq.113 .or. igamd0.eq.113) then
                  VavgDt = sit_tmp(si_dihedral_ene)
                  sigmaVDt = sit2_tmp(si_dihedral_ene)
                  VavgPt = sit_tmp(si_dvdl)
                  sigmaVPt =  sit2_tmp(si_dvdl)
                else if (igamd.eq.114 .or. igamd0.eq.114) then
                  VavgDt = sit_tmp(si_dihedral_ene)
                  sigmaVDt = sit2_tmp(si_dihedral_ene)
                  VavgPt = sit_tmp(si_dvdl)
                  sigmaVPt =  sit2_tmp(si_dvdl)
                else if (igamd.eq.115 .or. igamd0.eq.115) then
                  VavgDt = sit_tmp(si_dihedral_ene)
                  sigmaVDt = sit2_tmp(si_dihedral_ene)
                  VavgPt = sit_tmp(si_dvdl)
                  sigmaVPt =  sit2_tmp(si_dvdl)
                 else if (igamd.eq.116 .or. igamd0.eq.116) then
                  VavgDt = sit_tmp(si_dihedral_ene)
                  sigmaVDt = sit2_tmp(si_dihedral_ene)
                  VavgPt = sit_tmp(si_dvdl)
                  sigmaVPt =  sit2_tmp(si_dvdl)
                else if (igamd.eq.117 .or. igamd0.eq.117) then
                  VavgDt = sit_tmp(si_dihedral_ene)
                  sigmaVDt = sit2_tmp(si_dihedral_ene)
                  VavgPt = sit_tmp(si_dvdl)
                  sigmaVPt =  sit2_tmp(si_dvdl)
                else if (igamd.eq.118 .or. igamd0.eq.118) then
                  VavgDt = sit_tmp(si_dihedral_ene)
                  sigmaVDt = sit2_tmp(si_dihedral_ene)
                  VavgPt = sit_tmp(si_dvdl)
                  sigmaVPt =  sit2_tmp(si_dvdl)
                else if (igamd.eq.119 .or. igamd0.eq.119) then
                  VavgDt = sit_tmp(si_dihedral_ene)
                  sigmaVDt = sit2_tmp(si_dihedral_ene)
                  VavgPt = sit_tmp(si_dvdl)
                  sigmaVPt =  sit2_tmp(si_dvdl)
                else if (igamd.eq.120 .or. igamd0.eq.120) then
                  VavgDt = sit_tmp(si_dihedral_ene)
                  sigmaVDt = sit2_tmp(si_dihedral_ene)
                  VavgPt = sit_tmp(si_dvdl)
                  sigmaVPt =  sit2_tmp(si_dvdl)
                else if ((igamd.eq.101) .or. (igamd0.eq.101)) then
                  VavgDt = sit_tmp(si_dihedral_ene)
                  sigmaVDt = sit2_tmp(si_dihedral_ene)
                  VavgPt = ti_ene_tmp(1,si_elect_14_ene) + ti_others_ene_tmp(1,si_elect_14_ene)
                  sigmaVPt = sqrt(ti_ene_tmp2(1,si_elect_14_ene)**2 + ti_others_ene_tmp2(1,si_elect_14_ene)**2)
                else if ((igamd.eq.102) .or. (igamd0.eq.102)) then
                  VavgDt = sit_tmp(si_dihedral_ene)
                  sigmaVDt = sit2_tmp(si_dihedral_ene)
                  VavgPt = ti_ene_tmp(1,si_vdw_14_ene) + ti_others_ene_tmp(1,si_vdw_14_ene)
                  sigmaVPt = sqrt(ti_ene_tmp2(1,si_vdw_14_ene)**2 + ti_others_ene_tmp2(1,si_vdw_14_ene)**2)
                else if ((igamd.eq.103) .or. (igamd0.eq.103)) then
                  VavgDt = sit_tmp(si_dihedral_ene)
                  sigmaVDt = sit2_tmp(si_dihedral_ene)
                  VavgPt = ti_ene_tmp(1,si_elect_ene) + ti_others_ene_tmp(1,si_elect_ene)
                  sigmaVPt = sqrt(ti_ene_tmp2(1,si_elect_ene)**2 + ti_others_ene_tmp2(1,si_elect_ene)**2)
                else if ((igamd.eq.104) .or. (igamd0.eq.104)) then
                  VavgDt = sit_tmp(si_dihedral_ene)
                  sigmaVDt = sit2_tmp(si_dihedral_ene)
                  VavgPt = ti_ene_tmp(1,si_vdw_ene) + ti_others_ene_tmp(1,si_vdw_ene)
                  sigmaVPt = sqrt(ti_ene_tmp2(1,si_vdw_ene)**2 + ti_others_ene_tmp2(1,si_vdw_ene)**2)
                else
                  VavgDt = sit_tmp(si_dihedral_ene)
                  sigmaVDt = sit2_tmp(si_dihedral_ene)
                  VavgPt = sit_tmp(si_pot_ene)
                  sigmaVPt = sit2_tmp(si_pot_ene)
                end if
              end if

              sit_tmp(:) = sit(:)
              sit2_tmp(:) = sit2(:)
            end if
          end if

          ! GaMD
          if ((igamd.gt.0).or.(igamd0.gt.0)) then
            ! run cMD to calculate Vmax, Vmin, Vavg, sigmaV and then GaMD parameters
            if ((total_nstep.le.ntcmd).and.onstep) then
               if (total_nstep.ge.ntcmdprep) then
                 if ((igamd.eq.4 .or. igamd.eq.5) .or. (igamd0.eq.4 .or. igamd0.eq.5)) then
                   VmaxPt = max(VmaxPt,si(si_elect_ene)+si(si_vdw_ene)+si(si_vdw_14_ene)+si(si_elect_14_ene)+si(si_efield_ene))
                   VminPt = min(VminPt,si(si_elect_ene)+si(si_vdw_ene)+si(si_vdw_14_ene)+si(si_elect_14_ene)+si(si_efield_ene))
                 else if ((igamd.eq.6 .or. igamd.eq.7) .or. (igamd0.eq.6 .or. igamd0.eq.7)) then
                   VPt = ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                           ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                           ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
                   VmaxPt = max(VmaxPt,VPt)
                   VminPt = min(VminPt,VPt)
                 else if ((igamd.eq.8 .or. igamd.eq.9) .or. (igamd0.eq.8 .or. igamd0.eq.9)) then
                   VmaxPt = max(VmaxPt,ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene))
                   VminPt = min(VminPt,ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene))
                 else if ((igamd.eq.10 .or. igamd.eq.11) .or. (igamd0.eq.10 .or. igamd0.eq.11)) then
                   VPt = ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                            ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
                   VmaxPt = max(VmaxPt, VPt)
                   VminPt = min(VminPt, VPt)
                 else if ((igamd.eq.12 .or. igamd.eq.13) .or. (igamd0.eq.12 .or. igamd0.eq.13)) then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if ((igamd.eq.14 .or. igamd.eq.15) .or. (igamd0.eq.14 .or. igamd0.eq.15)) then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if ((igamd.eq.16 .or. igamd.eq.17) .or. (igamd0.eq.16 .or. igamd0.eq.17)) then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if (igamd.eq.18 .or. igamd0.eq.18)then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if (igamd.eq.19 .or. igamd0.eq.19)then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if (igamd.eq.20 .or. igamd0.eq.20)then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if (igamd.eq.21 .or. igamd0.eq.21)then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if (igamd.eq.22 .or. igamd0.eq.22)then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if (igamd.eq.23 .or. igamd0.eq.23)then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if (igamd.eq.24 .or. igamd0.eq.24)then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if (igamd.eq.25 .or. igamd0.eq.25)then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if (igamd.eq.26 .or. igamd0.eq.26)then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if (igamd.eq.27 .or. igamd0.eq.27)then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if (igamd.eq.28 .or. igamd0.eq.28)then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if ((igamd.eq.110 .or. igamd.eq.111) .or. (igamd0.eq.110 .or. igamd0.eq.111)) then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if ((igamd.eq.112 .or. igamd.eq.113) .or. (igamd0.eq.112 .or. igamd0.eq.113)) then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if ((igamd.eq.114 .or. igamd.eq.115) .or. (igamd0.eq.114 .or. igamd0.eq.115)) then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if ((igamd.eq.116 .or. igamd.eq.117) .or. (igamd0.eq.116 .or. igamd0.eq.117)) then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if ((igamd.eq.118 .or. igamd.eq.119) .or. (igamd0.eq.118 .or. igamd0.eq.119)) then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if (igamd.eq.120 .or. igamd0.eq.120) then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                else if ((igamd.eq.100) .or. (igamd0.eq.100)) then
                   VmaxPt = max(VmaxPt, ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                           ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene))
                   VminPt = min(VminPt, ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                           ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene))
                 else if ((igamd.eq.101) .or. (igamd0.eq.101)) then
                   VmaxPt = max(VmaxPt,ti_ene(1,si_elect_14_ene) + ti_others_ene(1,si_elect_14_ene))
                   VminPt = min(VminPt,ti_ene(1,si_elect_14_ene) + ti_others_ene(1,si_elect_14_ene))
                 else if ((igamd.eq.102) .or. (igamd0.eq.102)) then
                   VmaxPt = max(VmaxPt,ti_ene(1,si_vdw_14_ene) + ti_others_ene(1,si_vdw_14_ene))
                   VminPt = min(VminPt,ti_ene(1,si_vdw_14_ene) + ti_others_ene(1,si_vdw_14_ene))
                 else if ((igamd.eq.103) .or. (igamd0.eq.103)) then
                   VmaxPt = max(VmaxPt,ti_ene(1,si_elect_ene) + ti_others_ene(1,si_elect_ene))
                   VminPt = min(VminPt,ti_ene(1,si_elect_ene) + ti_others_ene(1,si_elect_ene))
                 else if ((igamd.eq.104) .or. (igamd0.eq.104)) then
                   VmaxPt = max(VmaxPt,ti_ene(1,si_vdw_ene) + ti_others_ene(1,si_vdw_ene))
                   VminPt = min(VminPt,ti_ene(1,si_vdw_ene) + ti_others_ene(1,si_vdw_ene))
                else
                  VmaxPt = max(VmaxPt,si(si_pot_ene))
                  VminPt = min(VminPt,si(si_pot_ene))
                end if
   
                 if ((igamd.eq.6 .or. igamd.eq.7) .or. (igamd0.eq.6 .or. igamd0.eq.7)) then
                   VDt = si(si_pot_ene) + ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene)
                   VmaxDt = max(VmaxDt,VDt)
                   VminDt = min(VminDt,VDt)
                 else if ((igamd.eq.8 .or. igamd.eq.9) .or. (igamd0.eq.8 .or. igamd0.eq.9)) then
                   VmaxDt = max(VmaxDt,si(si_pot_ene) + &
                      ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                      ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene) + &
                      ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene))
                   VminDt = min(VminDt,si(si_pot_ene) + &
                      ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                      ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene) + &
                      ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene))
                 else if ((igamd.eq.10 .or. igamd.eq.11) .or. (igamd0.eq.10 .or. igamd0.eq.11)) then
                   VDt = si(si_pot_ene) + &
                           ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                           ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene)
                   VmaxDt = max(VmaxDt, VDt)
                   VminDt = min(VminDt, VDt)
                  elseif((igamd.eq.13) .or. (igamd0.eq.13))then
                    VmaxDt = max(VmaxDt,si(si_gamd_ppi))
                    VminDt = min(VminDt,si(si_gamd_ppi))
                  elseif((igamd.eq.15) .or. (igamd0.eq.15))then
                   VmaxDt = max(VmaxDt,si(si_gamd_ppi))
                   VminDt = min(VminDt,si(si_gamd_ppi))
                   elseif((igamd.eq.17) .or. (igamd0.eq.17))then
                   VmaxDt = max(VmaxDt,si(si_gamd_ppi))
                   VminDt = min(VminDt,si(si_gamd_ppi))
                   elseif((igamd.eq.18) .or. (igamd0.eq.18))then
                   VmaxDt = max(VmaxDt,si(si_gamd_ppi))
                   VminDt = min(VminDt,si(si_gamd_ppi))
                 elseif(igamd.eq.20.or.igamd0.eq.20)then
                   VmaxDt = max(VmaxDt,si(si_gamd_ppi))
                   VminDt = min(VminDt,si(si_gamd_ppi))
                   VmaxBt = max(VmaxBt,si(si_gamd_ppi2))
                   VminBt = min(VminBt,si(si_gamd_ppi2))
                 elseif((igamd.eq.21) .or. (igamd0.eq.21))then
                   VmaxDt = max(VmaxDt,si(si_gamd_ppi))
                   VminDt = min(VminDt,si(si_gamd_ppi))
                   VmaxBt = max(VmaxBt,si(si_gamd_ppi2))
                   VminBt = min(VminBt,si(si_gamd_ppi2))
                 elseif((igamd.eq.22) .or. (igamd0.eq.22))then
                    VmaxDt = max(VmaxDt,si(si_gamd_ppi))
                    VminDt = min(VminDt,si(si_gamd_ppi))
                    VmaxBt = max(VmaxBt,si(si_gamd_ppi2))
                    VminBt = min(VminBt,si(si_gamd_ppi2))
                 elseif((igamd.ge.23.and.igamd.le.28) .or. (igamd0.ge.23.and.igamd0.le.28))then
                    VmaxDt = max(VmaxDt,si(si_gamd_ppi))
                    VminDt = min(VminDt,si(si_gamd_ppi))
                    VmaxBt = max(VmaxBt,si(si_gamd_ppi2))
                    VminBt = min(VminBt,si(si_gamd_ppi2))
                 else if ((igamd.eq.100) .or. (igamd0.eq.100)) then
                   VmaxDt = max(VmaxDt, si(si_pot_ene) + &
                           ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                           ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene))
                   VminDt = min(VminDt, si(si_pot_ene) + &
                           ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                           ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene))
                 else
                VmaxDt = max(VmaxDt,si(si_dihedral_ene))
                VminDt = min(VminDt,si(si_dihedral_ene))
                 end if

                if (mod(total_nstep, ntave).eq.0 .and. onstep) then
                  VmaxP = max(VmaxP,VmaxPt)
                  VminP = min(VminP,VminPt)
                  VavgP = VavgPt
                  sigmaVP = sigmaVPt
     
                  VmaxD = max(VmaxD,VmaxDt)
                  VminD = min(VminD,VminDt)
                  VavgD = VavgDt
                  sigmaVD = sigmaVDt

                   VmaxB = max(VmaxB,VmaxBt)
                   VminB = min(VminB,VminBt)
                   VavgB = VavgBt
                   sigmaVB = sigmaVBt

                  write(mdout,'(a,i10,4f14.4)') ' Energy statistics: step,VmaxP,VminP,VavgP,sigmaVP = ', &
                    total_nstep,VmaxP,VminP,VavgP,sigmaVP
                  write(mdout,'(a,i10,4f14.4)') ' Energy statistics: step,VmaxD,VminD,VavgD,sigmaVD = ', &
                    total_nstep,VmaxD,VminD,VavgD,sigmaVD
                  if((igamd.ge.20.and.igamd.le.28).or.(igamd0.ge.20.and.igamd0.le.28))then
                  write(mdout,'(a,i10,4f14.4)') ' Energy statistics: step,VmaxB,VminB,VavgB,sigmaVB = ', &
                     total_nstep,VmaxB,VminB,VavgB,sigmaVB
                  endif
                end if
              end if
            end if
     
            ! update Vmax, Vmin, Vavg, sigmaV and then GaMD parameters after adding boost potential
            if ((total_nstep.le.(ntcmd+nteb)).and.(total_nstep.ge.ntcmd).and.(mod(total_nstep,nrespa).eq.0).and.onstep) then
             if ( igamd.eq.1 .or. igamd.eq.3 .or. igamd.eq.4 .or. igamd.eq.5 .or. &
                   igamd.eq.6 .or. igamd.eq.7 .or. igamd.eq.8 .or. igamd.eq.9 .or. &
                   igamd.eq.10 .or. igamd.eq.11 .or.igamd.eq.12 .or. igamd.eq.13 .or. igamd.eq.14 .or. igamd.eq.100 .or. &
                   igamd.eq.15 .or.igamd.eq.16 .or. igamd.eq.17 .or. igamd.eq.18 .or.igamd.eq.19.or.&
                   igamd.eq.20 .or.igamd.eq.21 .or. igamd.eq.22.or.igamd.eq.23.or.igamd.eq.24.or.igamd.eq.25.or.&
                   igamd.eq.26 .or.igamd.eq.27 .or. igamd.eq.28.or. igamd.eq.101 .or. igamd.eq.102 .or. igamd.eq.103 .or. &
                   igamd.eq.104 .or.igamd.eq.110 .or. igamd.eq.111 .or. igamd.eq.112 .or. igamd.eq.113 .or. & 
                   igamd.eq.114 .or. igamd.eq.115 .or. igamd.eq.116 .or. igamd.eq.117.or. &
                   igamd.eq.118 .or. igamd.eq.119 .or. igamd.eq.120 .or. igamd.eq.20) then
               if (total_nstep.ge.(ntcmd+ntebprep)) then
                 if (igamd.eq.4 .or. igamd.eq.5) then
                   VmaxPt = max(VmaxPt,si(si_elect_ene)+si(si_vdw_ene)+si(si_vdw_14_ene)+si(si_elect_14_ene)+si(si_efield_ene))
                   VminPt = min(VminPt,si(si_elect_ene)+si(si_vdw_ene)+si(si_vdw_14_ene)+si(si_elect_14_ene)+si(si_efield_ene))
                 else if ((igamd.eq.6 .or. igamd.eq.7) .or. (igamd0.eq.6 .or. igamd0.eq.7)) then
                   VPt = ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene)+&
                           ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                           ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
                   VmaxPt = max(VmaxPt,VPt)
                   VminPt = min(VminPt,VPt)
                 else if ((igamd.eq.8 .or. igamd.eq.9) .or. (igamd0.eq.8 .or. igamd0.eq.9)) then
                   VmaxPt = max(VmaxPt,ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene))
                   VminPt = min(VminPt,ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene))
                 else if ((igamd.eq.10 .or. igamd.eq.11) .or. (igamd0.eq.10 .or. igamd0.eq.11)) then
                   VPt = ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                            ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
                   VmaxPt = max(VmaxPt, VPt)
                   VminPt = min(VminPt, VPt)
                 else if ((igamd.eq.12 .or. igamd.eq.13) .or. (igamd0.eq.12 .or. igamd0.eq.13)) then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if ((igamd.eq.14 .or. igamd.eq.15) .or. (igamd0.eq.14 .or. igamd0.eq.15)) then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if ((igamd.eq.16 .or. igamd.eq.17) .or. (igamd0.eq.16 .or. igamd0.eq.17)) then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if ((igamd.eq.18 .or. igamd0.eq.18))then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if ((igamd.eq.19 .or. igamd0.eq.19))then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if ((igamd.eq.20 .or. igamd.eq.21) .or. (igamd0.eq.20 .or. igamd0.eq.21)) then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if ((igamd.eq.22 .or. igamd0.eq.22))then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if ((igamd.ge.23 .and. igamd.le.28).or.(igamd0.ge.23.and.igamd0.le.28))then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if ((igamd.eq.110 .or. igamd.eq.111) .or. (igamd0.eq.110 .or. igamd0.eq.111)) then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if ((igamd.eq.112 .or. igamd.eq.113) .or. (igamd0.eq.112 .or. igamd0.eq.113)) then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if ((igamd.eq.114 .or. igamd.eq.115) .or. (igamd0.eq.114 .or. igamd0.eq.115)) then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                 else if ((igamd.eq.116 .or. igamd.eq.117) .or. (igamd0.eq.116 .or. igamd0.eq.117)) then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                  else if ((igamd.eq.118 .or. igamd.eq.119) .or. (igamd0.eq.118 .or. igamd0.eq.119)) then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                  else if (igamd.eq.120 .or. igamd0.eq.120) then
                   VmaxPt = max(VmaxPt, si(si_dvdl))
                   VminPt = min(VminPt, si(si_dvdl))
                else if ((igamd.eq.100) .or. (igamd0.eq.100)) then
                   VmaxPt = max(VmaxPt, ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                           ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene))
                   VminPt = min(VminPt, ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                           ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene))
                 else if ((igamd.eq.101) .or. (igamd0.eq.101)) then
                   VmaxPt = max(VmaxPt,ti_ene(1,si_elect_14_ene) + ti_others_ene(1,si_elect_14_ene))
                   VminPt = min(VminPt,ti_ene(1,si_elect_14_ene) + ti_others_ene(1,si_elect_14_ene))
                 else if ((igamd.eq.102) .or. (igamd0.eq.102)) then
                   VmaxPt = max(VmaxPt,ti_ene(1,si_vdw_14_ene) + ti_others_ene(1,si_vdw_14_ene))
                   VminPt = min(VminPt,ti_ene(1,si_vdw_14_ene) + ti_others_ene(1,si_vdw_14_ene))
                 else if ((igamd.eq.103) .or. (igamd0.eq.103)) then
                   VmaxPt = max(VmaxPt,ti_ene(1,si_elect_ene) + ti_others_ene(1,si_elect_ene))
                   VminPt = min(VminPt,ti_ene(1,si_elect_ene) + ti_others_ene(1,si_elect_ene))
                 else if ((igamd.eq.104) .or. (igamd0.eq.104)) then
                   VmaxPt = max(VmaxPt,ti_ene(1,si_vdw_ene) + ti_others_ene(1,si_vdw_ene))
                   VminPt = min(VminPt,ti_ene(1,si_vdw_ene) + ti_others_ene(1,si_vdw_ene))
                  else
                    VmaxPt = max(VmaxPt,si(si_pot_ene))
                    VminPt = min(VminPt,si(si_pot_ene))
                  end if
               endif
                ! update GaMD parameters
               update_gamd = .false.
               if (total_nstep.eq.ntcmd) update_gamd = .true.
                  call gamd_calc_Vstat(VmaxPt,VminPt,VavgPt,sigmaVPt,cutoffP,iVm,ntave,total_nstep, &
                    update_gamd,VmaxP,VminP,VavgP,sigmaVP)
               if( update_gamd ) then
                 call gamd_calc_E_k0(iEP,sigma0P,VmaxP,VminP,VavgP,sigmaVP,EthreshP,k0P)
                 kP = k0P/(VmaxP - VminP)
                 write(mdout,'(a,i10,7f14.4)') '| GaMD updated parameters: step,VmaxP,VminP,VavgP,sigmaVP,k0P,kP,EthreshP = ', &
                                               total_nstep,VmaxP,VminP,VavgP,sigmaVP,k0P,kP,EthreshP
               endif
                end if
             if ( igamd.eq.2 .or. igamd.eq.3 .or. igamd.eq.5 .or. igamd.eq.7 .or. igamd.eq.9 .or. &
                  igamd.eq.11 .or. igamd.eq.13 .or. igamd.eq.15 .or. igamd.eq.17.or.igamd.eq.18.or.&
                  igamd.eq.19.or. igamd.eq.20.or. igamd.eq.21.or. igamd.eq.22.or.igamd.eq.23.or.&
                  igamd.eq.24.or.igamd.eq.25.or.igamd.eq.26.or.igamd.eq.27.or.igamd.eq.28.or.igamd.eq.117.or.igamd.eq.100 ) then
               if (total_nstep.ge.(ntcmd+ntebprep)) then
                   if( igamd.eq.13 .or.  igamd0.eq.13)then
                   VmaxDt = max(VmaxDt,si(si_gamd_ppi))
                   VminDt = min(VminDt,si(si_gamd_ppi))
                   elseif(igamd.eq.15 .or. igamd0.eq.15)then
                   VmaxDt = max(VmaxDt,si(si_gamd_ppi))
                   VminDt = min(VminDt,si(si_gamd_ppi))
                   elseif(igamd.eq.17 .or. igamd0.eq.17)then
                   VmaxDt = max(VmaxDt,si(si_gamd_ppi))
                   VminDt = min(VminDt,si(si_gamd_ppi))
                   elseif(igamd.eq.18 .or. igamd0.eq.18)then
                   VmaxDt = max(VmaxDt,si(si_gamd_ppi))
                   VminDt = min(VminDt,si(si_gamd_ppi))
                  elseif(igamd.eq.19 .or. igamd0.eq.19)then
                   VmaxDt = max(VmaxDt,si(si_gamd_ppi))
                   VminDt = min(VminDt,si(si_gamd_ppi))
                  elseif((igamd.ge.20 .and.igamd.le.28) .or. (igamd0.ge.20.and.igamd0.le.28))then
                   VmaxDt = max(VmaxDt,si(si_gamd_ppi))
                   VminDt = min(VminDt,si(si_gamd_ppi))
                   VmaxBt = max(VmaxBt,si(si_gamd_ppi2))
                   VminBt = min(VminBt,si(si_gamd_ppi2))
                  elseif((igamd.eq.6 .or. igamd.eq.7) .or. (igamd0.eq.6 .or. igamd0.eq.7)) then
                   VDt = si(si_pot_ene) + ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene)
                   VmaxDt = max(VmaxDt,VDt)
                   VminDt = min(VminDt,VDt)
                 else if ((igamd.eq.8 .or. igamd.eq.9) .or. (igamd0.eq.8 .or. igamd0.eq.9)) then
                   VmaxDt = max(VmaxDt,si(si_pot_ene) + &
                      ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                      ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene) + &
                      ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene))
                   VminDt = min(VminDt,si(si_pot_ene) + &
                      ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                      ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene) + &
                      ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene))
                 else if ((igamd.eq.10 .or. igamd.eq.11) .or. (igamd0.eq.10 .or. igamd0.eq.11)) then
                   VDt = si(si_pot_ene) + &
                           ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                           ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene)
                   VmaxDt = max(VmaxDt, VDt)
                   VminDt = min(VminDt, VDt)
                 else if ((igamd.eq.100) .or. (igamd0.eq.100)) then
                   VmaxDt = max(VmaxDt, si(si_pot_ene) + &
                           ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                           ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene))
                   VminDt = min(VminDt, si(si_pot_ene) + &
                           ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                           ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene))
                 else
                  VmaxDt = max(VmaxDt,si(si_dihedral_ene))
                  VminDt = min(VminDt,si(si_dihedral_ene))
                 endif
               endif
                ! update GaMD parameters
               update_gamd = .false.
               if (total_nstep.eq.ntcmd) update_gamd = .true.
                  call gamd_calc_Vstat(VmaxDt,VminDt,VavgDt,sigmaVDt,cutoffD,iVm,ntave,total_nstep, &
                    update_gamd,VmaxD,VminD,VavgD,sigmaVD)
                   if((igamd.ge.20.and.igamd.le.28).or.(igamd0.ge.20.and.igamd0.le.28)) call gamd_calc_Vstat(VmaxBt,VminBt,VavgBt,&
                           sigmaVBt,cutoffB,iVm,ntave,total_nstep,update_gamd,VmaxB,VminB,VavgB,sigmaVB)
               if( update_gamd ) then
                  call gamd_calc_E_k0(iED,sigma0D,VmaxD,VminD,VavgD,sigmaVD,EthreshD,k0D)
                  kD = k0D/(VmaxD - VminD)
                 write(mdout,'(a,i10,7f14.4)') '| GaMD updated parameters: step,VmaxD,VminD,VavgD,sigmaVD,k0D,kD,EthreshD = ', &
                                               total_nstep,VmaxD,VminD,VavgD,sigmaVD,k0D,kD,EthreshD
                 if((igamd.ge.20.and.igamd.le.28).or.(igamd0.ge.20.and.igamd0.le.28))then
                   call gamd_calc_E_k0(iEB,sigma0B,VmaxB,VminB,VavgB,sigmaVB,EthreshB,k0B)
                   k_B = k0B/(VmaxB - VminB)
                   write(mdout,'(a,i10,7f14.4)') '| GaMD updated parameters: step,VmaxB,VminB,VavgB,sigmaVB,k0B,kB,EthreshB = ', &
                                               total_nstep,VmaxB,VminB,VavgB,sigmaVB,k0B,k_B,EthreshB
                 endif
               endif
              end if
            end if

            ! save gamd restart file
            if ( (mod(total_nstep, ntwr).eq.0) .and. (total_nstep.le.(ntcmd+nteb)) .and. onstep) then
              open(unit=gamdres,file=gamdres_name,action='write')
               if (igamd.eq.1 .or. igamd.eq.3 .or. igamd.eq.4 .or. igamd.eq.5 .or. &
                   igamd.eq.6 .or. igamd.eq.7 .or. igamd.eq.8 .or. igamd.eq.9 .or. &
                   igamd.eq.10 .or. igamd.eq.11 .or. igamd.eq.12 .or. igamd.eq.13 .or.igamd.eq.100 .or. &
                   igamd.eq.14 .or. igamd.eq.15 .or. igamd.eq.16 .or. igamd.eq.17 .or.igamd.eq.18 .or.&
                   igamd.eq.19.or.igamd.eq.20.or.igamd.eq.21.or.igamd.eq.22.or.igamd.eq.23.or.igamd.eq.24 .or.&
                   igamd.eq.25.or.igamd.eq.26.or.igamd.eq.27.or.igamd.eq.28.or.igamd.eq.101 .or. igamd.eq.102 .or. &
                   igamd.eq.103 .or.igamd.eq.104.or.igamd.eq.110 .or. igamd.eq.111 .or. igamd.eq.112 .or. igamd.eq.113 .or.&
                   igamd.eq.114 .or. igamd.eq.115 .or. igamd.eq.116 .or. igamd.eq.117.or. &
                   igamd.eq.118 .or. igamd.eq.119 .or. igamd.eq.120) then
                 write(gamdres,*) VmaxP, VminP, VavgP, sigmaVP
                 write(mdout,'(a,4f16.4)')'| Saved GaMD restart parameters: VmaxP,VminP,VavgP,sigmaVP = '
                 write(mdout,'(a,4f16.4)')'| ', VmaxP,VminP,VavgP,sigmaVP
               end if
             if ( igamd.eq.2 .or. igamd.eq.3 .or. igamd.eq.5 .or. igamd.eq.7 .or. igamd.eq.9 .or. igamd.eq.11 .or.&
                  igamd.eq.13 .or. igamd.eq.15 .or. igamd.eq.17.or.igamd.eq.18.or.igamd.eq.19.or.igamd.eq.20.or.igamd.eq.21&
                  .or.igamd.eq.22.or.igamd.eq.23.or.igamd.eq.24.or.igamd.eq.25.or.igamd.eq.26.or.igamd.eq.27.or.igamd.eq.28&
                  .or.igamd.eq.100) then
                 write(gamdres,*) VmaxD, VminD, VavgD, sigmaVD
                 write(mdout,'(a,4f16.4)')'| Saved GaMD restart parameters: VmaxD,VminD,VavgD,sigmaVD = '
                 write(mdout,'(a,4f16.4)')'| ', VmaxD,VminD,VavgD,sigmaVD
               end if
               if(igamd.ge.20.and.igamd.le.28)then
                 write(gamdres,*) VmaxB, VminB, VavgB, sigmaVB
                 write(mdout,'(a,4f16.4)')'| Saved GaMD restart parameters: VmaxB,VminB,VavgB,sigmaVB = '
                 write(mdout,'(a,4f16.4)')'| ', VmaxB,VminB,VavgB,sigmaVB
               endif
               close(gamdres)
               write(mdout,*) ''
            end if

            ! swap atomic coordinates of TI selected ligand residue with those of the bound ligand
            if (ibblig.gt.0 .and. onstep) then
            if ( (mod(total_nstep, ntwx).eq.0) ) then
#ifdef CUDA
                      call gpu_download_crd(crd)
#endif
                  if ( ibblig.eq.2 ) then
                    do i = 1, nlig
                       atm_i = ti_atm_lst(1,1) + ti_ti_atm_cnt(1)*(i-1) + atom_l-1
                       pos_tmp = crd(1:3,atm_i)
                       ! use PBC wrapped atomic coordinate
                       ! call PBC_r2f(crd(1:3,atm_i), fr)
                       ! call PBC_f2r(fr, pos_tmp, 0, 0, 0)

                       ift = mod(total_nstep, ntwin)/ntwx
                       if (ift.lt.nftau) lcrd0(1:3,ift+1,i) = pos_tmp
                       if (ift.ge.nfmsd) lcrd1(1:3,ift-nfmsd+1,i) = pos_tmp
                    end do
                  end if

                  ! calculate protein-ligand atom distance
                  dmin = 1.0d99
                  if (ibblig.eq.1) then
                    do i = 1, nlig
                       atm_i = ti_atm_lst(1,1) + ti_ti_atm_cnt(1)*(i-1) + atom_l-1
                       if (ntb.eq.0) then
                          dlig = distance(crd(1:3,atom_p), crd(1:3,atm_i))
                       else
                          dlig = PBC_distance(crd(1:3,atom_p), crd(1:3,atm_i))
                       end if
                       if (dmin.gt.dlig) then
                          dmin = dlig
                          blig_min = i
                       end if
                    end do
                    if (dmin.le.dblig) then
                       blig = blig_min
                       write(mdout,'(a,2I10,F14.6)')'| GaMD step, ligand, distance_min: ', total_nstep, blig, dmin
                       write(gamdlig,'(2I10,F14.6)') total_nstep, blig, dmin
                    end if
                  ! calculate ligand MSD
                  else if (ibblig.eq.2 .and. (mod(total_nstep, ntwin).eq.0)) then
                    do i = 1, nlig
                       dlig = msd(lcrd0(:,:,i), lcrd1(:,:,i))
                       if (dmin.gt.dlig) then
                          dmin = dlig
                          blig_min = i
                       end if
                    end do
                    if (dmin.le.dblig) then
                       blig = blig_min
                       write(mdout,'(a,2I10,F14.6)')'| GaMD step, ligand, MSD_min: ', total_nstep, blig, dmin
                       write(gamdlig,'(2I10,F14.6)') total_nstep, blig, dmin
                    end if
                  end if 

                  ! swap ligand coordinates
                  if (blig.ne.blig0) then
                      write(mdout,'(a,2I10)')'| GaMD swap coordinates, forces and velocities of ligands: ', blig0, blig
#ifdef CUDA
                      call gpu_download_vel(vel)
                      call gpu_download_frc(frc)
#endif
                      do i = 1, ti_ti_atm_cnt(1)
                         atm_i = ti_atm_lst(1,i)
                         atm_j = atm_i + ti_ti_atm_cnt(1)*(blig-1)
                         ! coordinates
                         pos_tmp(1:3) = crd(1:3,atm_i)
                         crd(1:3,atm_i) = crd(1:3,atm_j)
                         crd(1:3,atm_j) = pos_tmp(1:3)
                         ! forces
                         pos_tmp(1:3) = frc(1:3,atm_i)
                         frc(1:3,atm_i) = frc(1:3,atm_j)
                         frc(1:3,atm_j) = pos_tmp(1:3)
                         ! velocities
                         pos_tmp(1:3) = vel(1:3,atm_i)
                         vel(1:3,atm_i) = vel(1:3,atm_j)
                         vel(1:3,atm_j) = pos_tmp(1:3)
                      end do
#ifdef CUDA
                      call gpu_upload_crd(crd)
                      call gpu_upload_frc(frc)
                      call gpu_upload_vel(vel)
#endif
                      blig = blig0
                  end if
               end if
            end if

          end if ! end GaMD

          if (((igamd.gt.0).or.(igamd0.gt.0)) .and. ntcmd.gt.0 .and. total_nstep.eq.ntcmd) then
                write(mdout,'(/,a,i10,/)') &
                  '| GaMD: Apply boost potential after finishing conventional MD; igamd = ', igamd
          end if

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
          if ( update_gamd )then
            if(.not.(igamd0.ge.20.and.igamd0.le.28))then
            call gpu_gamd_update(EthreshP, kP, EthreshD, kD)
            else
            call  gpu_gamd_triple_update(EthreshP, kP, EthreshD, kD,EthreshB, k_B)
           endif
         endif
#endif
        end if ! igamd.gt.0

      ! TI: Update SAMS lambda or dynamic lambda
      if (ti_mode .ne. 0 .and. onstep) then
        if (do_mbar .and. ifsams .gt. 0 .and. total_nstep .ge. sams_init_start) then
          if (total_nstep .ge. sams_stat_start) sams_do_average=.true.
          call sams_update(si(si_dvdl), new_list)
          !TBD call sams_update(si(si_dvdl))
          new_list=.true.
          if (sams_onstep) call sams_output(total_nstep)
        else if (ntave .gt. 0) then
          if (mod(nstep, ntave) .eq. 0 .and. dynlmb .gt. 0.d0) &
              call ti_dynamic_lambda
        end if
      end if


#ifdef MPI
        if (.not. usemidpoint) then
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
            if (abmd_selection_freq.gt.0.and.mod(abmd_mdstep,abmd_selection_freq).eq.0) then
              if (using_pme_potential) then
                 call mpi_gathervec(atm_cnt, vel)
              else if (using_gb_potential) then
                 call gb_mpi_gathervec(atm_cnt, vel)
              end if
              all_vels_valid = .true.
            end if
            call nfe_on_mdstep(pme_pot_ene%total, crd, vel)

            if (abmd_selection_freq.ne.0) then
              if (abmd_mdstep.gt.1.and.mod(abmd_mdstep,abmd_selection_freq).eq.0) then
#ifdef CUDA
                if (using_pme_potential) &
                  call gpu_force_new_neighborlist()
                call gpu_upload_crd(crd)
                call gpu_upload_vel(vel)
#else
                last_vel(:,:) = vel(:,:)
#endif
              end if
            end if
            if (bbmd_active.eq.1.and.mod(bbmd_mdstep,bbmd_exchg_freq).eq.0) then
#ifdef CUDA
              call gpu_upload_vel(vel)
#else
              if (all_vels_valid) then
                 last_vel(:,:) = vel(:,:)
              else
                 do atm_lst_idx = 1, my_atm_cnt
                    j = my_atm_lst(atm_lst_idx)
                    last_vel(:, j) = vel(:, j)
                 end do
              end if
#endif
            end if
          end if  ! infe
        end if  ! not midpoint
#endif /* MPI */
        !--------------------------------------------------------------------
        !--------------------------------------------------------------------

        if (nstep .ge. nstlim) then
#ifdef TIDECOMP
                if (ntwd .gt. 0) then
               call mddecomp(total_nstep)
           end if
#endif
           exit
        end if

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

! PLUMED
        if (plumed .ne. 0 .and. plumed_stopflag .ne. 0) exit
      end do ! Major cycle back to new step unless we have reached our limit:
      
#ifdef MPI
      ! If timlim has been exceeded exit the REMD loop
      if (timlim_exceeded .and.local_remd_method .ne. 0) then
        local_numexchg = mdloop
        if (master) &
          write(mdout,'(2(a,i8))') '| Wall clock limit reached (timlim = ', &
          timlim, ' s), exchange ', mdloop
        exit ! Exit the REMD loop.
      end if

    end do ! REMD loop (mdloop)

    ! Print averages, but not for REMD:

    if (master .and. local_remd_method .eq. 0) then
#else
    if (master) then
#endif

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
          if (isgld>0)then
            call sgld_fluc(.true.)
          endif
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
            if (igamd.gt.0 .or. igamd0.gt.0) &
              call ti_others_print_ene(ti_region, ti_others_ene_tmp(ti_region,:))

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
            if (igamd.gt.0 .or. igamd0.gt.0) &
              call ti_others_print_ene(ti_region, ti_others_ene_tmp2(ti_region,:))
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

  if (igamd.gt.0 .and. ibblig.gt.0) close(gamdlig)
  if (igamd.gt.0 .and. ibblig.eq.2) then
     if (allocated(lcrd0)) deallocate(lcrd0)
     if (allocated(lcrd1)) deallocate(lcrd1)
  end if

  call update_time(runmd_time)

! PLUMED
  if (plumed .ne. 0) call plumed_f_gfinalize()

  return

  540 format(/5x, ' A V E R A G E S   O V E R ', i7, ' S T E P S', /)
  541 format(/5x,' DV/DL, AVERAGES OVER ',i7,' STEPS',/)
  542 format('|',79('='))
  550 format(/5x, ' R M S  F L U C T U A T I O N S', /)
  580 format('STATISTICS OF EFFECTIVE BORN RADII OVER ',i7,' STEPS')
  590 format('ATOMNUM     MAX RAD     MIN RAD     AVE RAD     FLUCT')
  600 format(i4, 2x, 4f12.4)
  601 format(/,'| TI region ', i2, /)
end subroutine runmd


!*******************************************************************************
!
! Subroutine:  cenmas
!
! Description:  Calculate the translational and rotational kinetic energies
!               and velocities.
!
!*******************************************************************************

subroutine cenmas(atm_cnt, x, v, tmass, tmassinv, amass, xcm, vcm, ocm)

  use file_io_dat_mod
  use gbl_constants_mod
  use parallel_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: x(*)
  double precision      :: v(*)
  double precision      :: tmass
  double precision      :: tmassinv
  double precision      :: amass(*)
  double precision      :: xcm(*)
  double precision      :: vcm(*)
  double precision      :: ocm(*)

! Local variables:

  double precision      :: aamass
  double precision      :: comvel
  double precision      :: det
  integer               :: i, j, m, n
  integer               :: i0
  integer               :: lh(3)
  integer               :: mh(3)
  double precision      :: acm(3)
  double precision      :: tcm(3, 3)
  double precision      :: xx, xy, xz, yy, yz, zz
  double precision      :: x1, x2, x3
  double precision      :: ekcm, ekrot

  double precision, parameter   :: crit = 1.d-06

  i0 = 3 * atm_cnt

  ! calculate the center of mass coordinates:

  xcm(1) = 0.d0
  xcm(2) = 0.d0
  xcm(3) = 0.d0

  i = 0
  do j = 1, atm_cnt
    aamass = amass(j)
    xcm(1) = xcm(1) + x(i+1) * aamass
    xcm(2) = xcm(2) + x(i+2) * aamass
    xcm(3) = xcm(3) + x(i+3) * aamass
    i = i + 3
  end do

  xcm(1) = xcm(1) * tmassinv
  xcm(2) = xcm(2) * tmassinv
  xcm(3) = xcm(3) * tmassinv

  ! calculate velocity and translational kinetic energy of the center of mass:

  ekcm = 0.d0
  vcm(1) = 0.0d0
  vcm(2) = 0.0d0
  vcm(3) = 0.0d0

  i = 0
  do j = 1, atm_cnt
    aamass = amass(j)
    vcm(1) = vcm(1) + v(i+1) * aamass
    vcm(2) = vcm(2) + v(i+2) * aamass
    vcm(3) = vcm(3) + v(i+3) * aamass
    i = i + 3
  end do

  do i = 1, 3
    vcm(i) = vcm(i) * tmassinv
    ekcm = ekcm + vcm(i) * vcm(i)
  end do
  ekcm = ekcm * tmass * 0.5d0
  comvel = sqrt(vcm(1) * vcm(1) + vcm(2) * vcm(2) + vcm(3) * vcm(3))

  ! calculate the angular momentum about the center of mass:

  acm(1) = 0.0d0
  acm(2) = 0.0d0
  acm(3) = 0.0d0

  i = 0
  do j = 1, atm_cnt
    aamass = amass(j)
    acm(1) = acm(1) + (x(i+2) * v(i+3)-x(i+3) * v(i+2)) * aamass
    acm(2) = acm(2) + (x(i+3) * v(i+1)-x(i+1) * v(i+3)) * aamass
    acm(3) = acm(3) + (x(i+1) * v(i+2)-x(i+2) * v(i+1)) * aamass
    i = i + 3
 end do

  acm(1) = acm(1) - (xcm(2) * vcm(3)-xcm(3) * vcm(2)) * tmass
  acm(2) = acm(2) - (xcm(3) * vcm(1)-xcm(1) * vcm(3)) * tmass
  acm(3) = acm(3) - (xcm(1) * vcm(2)-xcm(2) * vcm(1)) * tmass

  ! calculate the inertia tensor:

  xx = 0.d0
  xy = 0.d0
  xz = 0.d0
  yy = 0.d0
  yz = 0.d0
  zz = 0.d0

  i = 0
  do j = 1, atm_cnt
    x1 = x(i+1) - xcm(1)
    x2 = x(i+2) - xcm(2)
    x3 = x(i+3) - xcm(3)
    aamass = amass(j)
    xx = xx + x1 * x1 * aamass
    xy = xy + x1 * x2 * aamass
    xz = xz + x1 * x3 * aamass
    yy = yy + x2 * x2 * aamass
    yz = yz + x2 * x3 * aamass
    zz = zz + x3 * x3 * aamass
    i = i + 3
  end do

  tcm(1,1) = yy + zz
  tcm(2,1) = -xy
  tcm(3,1) = -xz
  tcm(1,2) = -xy
  tcm(2,2) = xx + zz
  tcm(3,2) = -yz
  tcm(1,3) = -xz
  tcm(2,3) = -yz
  tcm(3,3) = xx + yy

  ! invert the inertia tensor:

  call matinv(tcm, 3, det, lh, mh)

  if (abs(det) .le. crit) then
    write(mdout, '(a,a)') error_hdr, 'Zero determinant calculated in cenmas()'
    call mexit(6, 1)
  end if

  ! calculate the angular velocity about the center of mass and the rotational
  ! kinetic energy:

  ekrot = 0.d0
  do m = 1, 3
    ocm(m) = 0.d0
    do n = 1, 3
      ocm(m) = ocm(m) + tcm(m,n) * acm(n)
    end do
    ekrot = ekrot + ocm(m) * acm(m)
  end do
  ekrot = ekrot * 0.5d0

  if (master) then
    write(mdout, '(/3x, a, f11.4, 3x, a, f11.4, 3x, a, f12.6)') 'KE Trans =', &
          ekcm, 'KE Rot =', ekrot, 'C.O.M. Vel =', comvel
  end if

  return

end subroutine cenmas

!*******************************************************************************
!
! Subroutine:  ti_cenmas
!
! Description:  Calculate the translational and rotational kinetic energies
!               and velocities. Modified for TI so that we do not double count
!               velocities in mode 1.
!
!*******************************************************************************

subroutine ti_cenmas(atm_cnt, x, v, tmass, tmassinv, amass, xcm, vcm, ocm)

  use file_io_dat_mod
  use gbl_constants_mod
  use parallel_dat_mod
  use pmemd_lib_mod
  use ti_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: x(*)
  double precision      :: v(*)
  double precision      :: tmass
  double precision      :: tmassinv
  double precision      :: amass(*)
  double precision      :: xcm(*)
  double precision      :: vcm(*)
  double precision      :: ocm(*)

! Local variables:

  double precision      :: aamass
  double precision      :: comvel
  double precision      :: det
  integer               :: i, j, m, n
  integer               :: i0
  integer               :: lh(3)
  integer               :: mh(3)
  double precision      :: acm(3)
  double precision      :: tcm(3, 3)
  double precision      :: xx, xy, xz, yy, yz, zz
  double precision      :: x1, x2, x3
  double precision      :: ekcm, ekrot

  double precision, parameter   :: crit = 1.d-06

  i0 = 3 * atm_cnt

  ! calculate the center of mass coordinates:

  xcm(1) = 0.d0
  xcm(2) = 0.d0
  xcm(3) = 0.d0


  i = 0
  do j = 1, atm_cnt
    if (ti_lst(2,j) .ne. 0) then
      i = i + 3
      cycle
    end if
    aamass = amass(j)
    xcm(1) = xcm(1) + x(i+1) * aamass
    xcm(2) = xcm(2) + x(i+2) * aamass
    xcm(3) = xcm(3) + x(i+3) * aamass
    i = i + 3
  end do

  xcm(1) = xcm(1) * tmassinv
  xcm(2) = xcm(2) * tmassinv
  xcm(3) = xcm(3) * tmassinv

  ! calculate velocity and translational kinetic energy of the center of mass:

  ekcm = 0.d0
  vcm(1) = 0.0d0
  vcm(2) = 0.0d0
  vcm(3) = 0.0d0

  i = 0
  do j = 1, atm_cnt
    if (ti_lst(2,j) .ne. 0) then
      i = i + 3
      cycle
    end if
    aamass = amass(j)
    vcm(1) = vcm(1) + v(i+1) * aamass
    vcm(2) = vcm(2) + v(i+2) * aamass
    vcm(3) = vcm(3) + v(i+3) * aamass
    i = i + 3
  end do

  do i = 1, 3
    vcm(i) = vcm(i) * tmassinv
    ekcm = ekcm + vcm(i) * vcm(i)
  end do
  ekcm = ekcm * tmass * 0.5d0
  comvel = sqrt(vcm(1) * vcm(1) + vcm(2) * vcm(2) + vcm(3) * vcm(3))

  ! calculate the angular momentum about the center of mass:

  acm(1) = 0.0d0
  acm(2) = 0.0d0
  acm(3) = 0.0d0

  i = 0
  do j = 1, atm_cnt
    if (ti_lst(2,j) .ne. 0) then
      i = i + 3
      cycle
    end if
    aamass = amass(j)
    acm(1) = acm(1) + (x(i+2) * v(i+3)-x(i+3) * v(i+2)) * aamass
    acm(2) = acm(2) + (x(i+3) * v(i+1)-x(i+1) * v(i+3)) * aamass
    acm(3) = acm(3) + (x(i+1) * v(i+2)-x(i+2) * v(i+1)) * aamass
    i = i + 3
  end do

  acm(1) = acm(1) - (xcm(2) * vcm(3)-xcm(3) * vcm(2)) * tmass
  acm(2) = acm(2) - (xcm(3) * vcm(1)-xcm(1) * vcm(3)) * tmass
  acm(3) = acm(3) - (xcm(1) * vcm(2)-xcm(2) * vcm(1)) * tmass

  ! calculate the inertia tensor:

  xx = 0.d0
  xy = 0.d0
  xz = 0.d0
  yy = 0.d0
  yz = 0.d0
  zz = 0.d0

  i = 0
  do j = 1, atm_cnt
    if (ti_lst(2,j) .ne. 0) then
      i = i + 3
      cycle
    end if
    x1 = x(i+1) - xcm(1)
    x2 = x(i+2) - xcm(2)
    x3 = x(i+3) - xcm(3)
    aamass = amass(j)
    xx = xx + x1 * x1 * aamass
    xy = xy + x1 * x2 * aamass
    xz = xz + x1 * x3 * aamass
    yy = yy + x2 * x2 * aamass
    yz = yz + x2 * x3 * aamass
    zz = zz + x3 * x3 * aamass
    i = i + 3
  end do

  tcm(1,1) = yy + zz
  tcm(2,1) = -xy
  tcm(3,1) = -xz
  tcm(1,2) = -xy
  tcm(2,2) = xx + zz
  tcm(3,2) = -yz
  tcm(1,3) = -xz
  tcm(2,3) = -yz
  tcm(3,3) = xx + yy

  ! invert the inertia tensor:

  call matinv(tcm, 3, det, lh, mh)

  if (abs(det) .le. crit) then
    write(mdout, '(a,a)') error_hdr, 'Zero determinant calculated in cenmas()'
    call mexit(6, 1)
  end if

  ! calculate the angular velocity about the center of mass and the rotational
  ! kinetic energy:

  ekrot = 0.d0
  do m = 1, 3
    ocm(m) = 0.d0
    do n = 1, 3
      ocm(m) = ocm(m) + tcm(m,n) * acm(n)
    end do
    ekrot = ekrot + ocm(m) * acm(m)
  end do
  ekrot = ekrot * 0.5d0

  if (master) then
    write(mdout, '(/3x, a, f11.4, 3x, a, f11.4, 3x, a, f12.6)') 'KE Trans =', &
          ekcm, 'KE Rot =', ekrot, 'C.O.M. Vel =', comvel
  end if

  return

end subroutine ti_cenmas


!*******************************************************************************
!
! Subroutine:   stopcm
!
! Description:  Remove Center of Mass transrotational motion.
!
!*******************************************************************************

subroutine stopcm(nr, x, v, xcm, vcm, ocm)

  use file_io_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: nr
  double precision      :: x(*)
  double precision      :: v(*)
  double precision      :: xcm(*)
  double precision      :: vcm(*)
  double precision      :: ocm(*)

! Local variables:

  integer               :: i, j, m
  double precision      :: x1, x2, x3

  ! stop the center of mass translation:

  i = 0
  do j = 1, nr
    do m = 1, 3
      i = i + 1
      v(i) = v(i) - vcm(m)
    end do
  end do

  ! stop the rotation about the center of mass:

  i = 0
  do j = 1, nr
      x1 = x(i+1) - xcm(1)
      x2 = x(i+2) - xcm(2)
      x3 = x(i+3) - xcm(3)
      v(i+1) = v(i+1) - ocm(2) * x3 + ocm(3) * x2
      v(i+2) = v(i+2) - ocm(3) * x1 + ocm(1) * x3
      v(i+3) = v(i+3) - ocm(1) * x2 + ocm(2) * x1
      i = i + 3
  end do

  if (master) write(mdout, 9008)
9008 format(/3x, 'Translational and rotational motion removed')

  return

end subroutine stopcm

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

!*******************************************************************************
!
! Subroutine:   re_position
!
! Description:  move the center of a set of atoms
!
! Input:  n  number of atoms
!         x  coordinates
!         xr reference coordinates
!         xc, yc, zc current center
!         x0, y0, z0 new center
!
! Output: x, xr  moved coordinates
!
!*******************************************************************************

subroutine re_position(n, ntr, x, xr, xc, yc, zc, x0, y0, z0, vec, verbose)

  use file_io_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: n
  integer               :: ntr
  double precision      :: x(3,n)
  double precision      :: xr(3,n)
  double precision      :: xc, yc, zc
  double precision      :: x0, y0, z0
  double precision      :: vec(3)
  integer               :: verbose

! Local variables:

  integer               :: i, j
  double precision      :: xd, yd, zd

  xd = x0 - xc
  yd = y0 - yc
  zd = z0 - zc

  if (master) &
    write(mdout, '(a, 3f10.6)')"| RE_POSITION Moving by ", xd, yd, zd

  vec(1) = vec(1) + xd
  vec(2) = vec(2) + yd
  vec(3) = vec(3) + zd

  do i = 1, n
    x(1,i) = xd + x(1,i)
    x(2,i) = yd + x(2,i)
    x(3,i) = zd + x(3,i)
  end do

  if (ntr .gt. 0) then
    do i = 1, n
      xr(1,i) = xd + xr(1,i)
      xr(2,i) = yd + xr(2,i)
      xr(3,i) = zd + xr(3,i)
    end do
  end if
  if (verbose .gt. 0 .and. master) then
    write(mdout, *)"*********** Coords moved ****************"
    write(mdout, *) "delta x,y,z ", xd, yd, zd
  end if

  return

end subroutine re_position

!*******************************************************************************
!
! Subroutine:  matinv
!
! Description: Standard ibm matrix inversion routine.
!
! Arguments:
!   a:  square matrix of dimension nxn
!   d:  resultant determinant
!   l:  work vector of length n
!   m:  work vector of length n
!
!*******************************************************************************

subroutine matinv(a, n, d, l, m)

  implicit none
!  implicit double precision    :: (a-h,o-z)

! Formal arguments:

  double precision      :: a(*)
  integer               :: n
  double precision      :: d
  integer               :: l(*)
  integer               :: m(*)

! Local variables:

  double precision      :: biga
  double precision      :: hold
  integer               :: i, ij, ik, iz
  integer               :: j, ji, jk, jp, jq, jr
  integer               :: k, ki, kj, kk
  integer               :: nk

    ! Search for largest element:

    d = 1.0d0
    nk = -n

    do 80 k = 1, n
      nk = nk + n
      l(k) = k
      m(k) = k
      kk = nk + k
      biga = a(kk)
      do 20 j = k, n
        iz = n * (j - 1)
        do 20 i = k, n
          ij = iz + i
          if (abs(biga) - abs(a(ij))) 15, 20, 20
 15         biga = a(ij)
            l(k) = i
            m(k) = j
 20         continue

    ! Interchange rows:

      j = l(k)
      if (j - k) 35, 35, 25
 25     ki = k - n
        do 30 i = 1, n
          ki = ki + n
          hold = -a(ki)
          ji = ki - k + j
          a(ki) = a(ji)
 30       a(ji) = hold

    ! Interchange columns:

 35     i = m(k)

      if (i - k) 45, 45, 38
 38     jp = n * (i - 1)
        do 40 j = 1, n
          jk = nk + j
          ji = jp + j
          hold = -a(jk)
          a(jk) = a(ji)
 40       a(ji) = hold

    ! Divide column by minus pivot:

 45     if (biga) 48, 46, 48
 46       d = 0.0d0
          goto 150
 48       do 55 i = 1, n
            if (i - k) 50, 55, 50
 50           ik = nk + i
              a(ik) = a(ik)/(-biga)
 55           continue

     ! Reduce matrix:

      do 65 i = 1, n
        ik = nk + i
        hold = a(ik)
        ij = i - n
        do 65 j = 1, n
          ij = ij + n
          if (i - k) 60, 65, 60
 60         if (j - k) 62, 65, 62
 62           kj = ij - i + k
              a(ij) = hold * a(kj) + a(ij)
 65       continue

    ! Divide row by pivot:

      kj = k - n
      do 75 j = 1, n
        kj = kj + n
        if (j - k) 70, 75, 70
 70       a(kj) = a(kj)/biga
 75     continue

    ! Product of pivots:

      d = d * biga

    ! Replace pivot by reciprocal:

      a(kk) = 1.0d0/biga

 80 continue

    ! Final row and column interchange:

    k = n
100 k = (k - 1)
    if (k) 150, 150, 105
105   i = l(k)
      if (i - k) 120, 120, 108
108     jq = n * (k - 1)
        jr = n * (i - 1)
        do 110 j = 1, n
          jk = jq + j
          hold = a(jk)
          ji = jr + j
          a(jk) = -a(ji)
110       a(ji) = hold
120     j = m(k)
      if (j - k) 100, 100, 125
125     ki = k - n
        do 130 i = 1, n
          ki = ki + n
          hold = a(ki)
          ji = ki - k + j
          a(ki) = -a(ji)
130       a(ji) = hold
    goto 100

150 return

end subroutine matinv

end module runmd_mod
