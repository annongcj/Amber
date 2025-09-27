#include "copyright.i"

!*******************************************************************************
!
! Program:      PMEMD 22,   originally based on sander version 6
!
! Description:  The Molecular Dynamics/NMR Refinement/Modeling Module of the
!               AMBER Package, high performance version.
!
! Current Principal Author: Ross C. Walker
!
!               GPU Code developed and maintained by
!               Ross C. Walker & Scott Le Grand
!
!               GPU Implmentation is Copyright Ross C. Walker
!               All Rights Reserved. Licensed for use in PMEMD.
!
!               Modern CPU Optimizations and Xeon Phi Implementation by
!               Ross C. Walker and Ashraf Bhuiyan
!
!               Original development of PMEMD was by Robert E. Duke of the
!               University of North Carolina-Chapel Hill Chemistry Department
!               and NIEHS.
!
!               PMEMD was originally based on Amber 6.0 sander source code.
!               and intended to be to be a fast implementation of the Particle
!               Mesh Ewald method. Under development by the group of Ross C.
!               Walker at UC San Diego it has now diverged from the original
!               Sander support to be an MD program in its own right. It now
!               supports PME, gas phase, GB and IPS simulations.
!
!               Current funding is through gift and royalty funding to Prof.
!               Ross C. Walker from NVIDIA, Intel, Pfizer and Exxact Corp. as
!               well as NSF SI2-SSE funding to Ross C. Walker and Adrian E.
!               Roitberg.
!
!               The original implementation was funded under NIH and NSF
!               grants to Prof. Lee Pedersen and Dr. Tom Darden.
!
!               Extensions to support the CHARMM force field were made by
!               Mike Crowley, Mark Williamson, and Ross C. Walker.
!
!               NVIDIA GPU acceleration support was written by Scott Le Grand,
!               Duncan Poole and Ross C. Walker, is Copyright to Ross C. Walker
!               and is licensed royalty free for use in pmemd.
!
!               Gas Phase Support (IGB=6) was developed by Dan Mermelstein and Ross
!               C. Walker.
!
!               Thermodynamic Integration was developed by Joe Kauss and Ross C.
!               Walker.
!
!               Replica exchange support, igb=8, and GB/SA written by
!               Jason Swails, Adrian Roitberg, and Ross C. Walker,
!
!               IPS support was added by Romelia Salomon and Ross C. Walker
!
!               Constant pH support was adapted from Sander by Jason Swails.
!               GPU support for Constant pH is by Perri Needham and Ross C.
!               Walker.
!
!               Constant Redox Potential support was expanded from the Constant pH
!               code by Vinicius Wilian D. Cruzeiro.
!
!               Force switching was added by Charles Lin and Ross C. Walker
!
!               Xeon Phi support was developed by Perri Needham, Ashraf Bhuiyan
!               and Ross C. Walker
!
!               Continuous constant pH MD based on GBNeck2 model was
!               developed by Robert C. Harris, Yandong Huang and Jana Shen
!
!*******************************************************************************

program pmemd

  use gb_alltasks_setup_mod
  use pme_alltasks_setup_mod
  use constraints_mod
#ifdef MPI
  use constantph_mod, only : cnstph_bcast, cnstph_setup, cph_igb, cphfirst_sol, cnstph_deallocate
  use constante_mod, only : cnste_bcast, cnste_setup, ce_igb, cefirst_sol, cnste_deallocate
  use get_cmdline_mod, only : get_cmdline_bcast, cpein_specified
#else
  use constantph_mod, only : cnstph_setup, cph_igb, cphfirst_sol, cnstph_deallocate
  use constante_mod, only : cnste_setup, ce_igb, cefirst_sol, cnste_deallocate
  use get_cmdline_mod, only : cpein_specified
#endif
  use external_mod, only : external_init, external_cleanup
  use extra_pnts_nb14_mod
  use emap_mod,only: pemap,qemap
  use ene_frc_splines_mod
  use gb_force_mod
#ifdef _OPENMP_
  use gb_ene_hybrid_mod
#endif
  use gb_ene_mod
  use gbl_constants_mod
  use gbl_datatypes_mod
  use pme_setup_mod
  use img_mod
  use inpcrd_dat_mod
  use master_setup_mod
#ifndef MPI
#ifndef NOXRAY
  use xray_interface_module, only: xray_init => init, xray_fini => finalize, xray_init_clean_up => init_clean_up
  use xray_globals_module, only: xray_active, ntwsf
#endif
#endif
  use barostats_mod, only : mcbar_setup
  use mcres_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
#ifdef MPI
  use multipmemd_mod, only : setup_groups, free_comms
  use pme_recip_midpoint_mod, only: pme_recip_setup
  use nb_exclusions_mod, only: alloc_nb_mask_data,make_atm_excl_mask_list
  use dynamics_dat_mod, only: alloc_atm_rel_crd
#endif
  use nb_pairlist_mod
  use nmr_calls_mod
  use parallel_dat_mod
  use cit_mod, only: bkt_size
  use processor_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
  use ramd_mod
  use random_mod
  use runfiles_mod
  use runmd_mod
  use runmin_mod
  use runreweight_mod
  use shake_mod
  use timers_mod
# ifdef BINTRAJ
  use AmberNetcdf_mod, only : NC_setupAmberNetcdf
# endif
  use bintraj_mod
  use pbc_mod
  use nbips_mod
#ifdef MPI
  use pme_force_mod,  only : alloc_force_mem, dealloc_force_mem, &
                             proc_alloc_force_mem
  use remd_mod,       only : remd_method, bcast_remd_method, remd_setup, &
                             slave_remd_setup, remd_cleanup, setup_pv_correction, &
                             rremd_type
  use remd_exchg_mod, only : setup_remd_randgen
  use hybridsolvent_remd_mod, only : hybridsolvent_remd_bcast, hybridsolvent_remd_atm_cnt
#endif
  use amd_mod
  use gamd_mod
  use scaledMD_mod
  use file_io_mod,       only : amopen
  use ti_mod
! SGLD
   use sgld_mod, only : psgld

  use sams_mod
#ifdef GTI
  use gti_mod
  use reaf_mod
#endif


! Modified by Feng Pan
   use nfe_setup_mod, only : &
         nfe_on_pmemd_init => on_pmemd_init, &
         nfe_on_pmemd_exit => on_pmemd_exit
! ----------------------------------
! DG: for neb calculation
   use neb_mod
   use nebread_mod

!PHMD
  use phmd_mod, only : phmdzero

  implicit none

! Local variables:

  double precision      :: max_erfc_relerr
  integer               :: new_stack_limit    ! new stack limit
  integer               :: i
  integer               :: num_ints = 0
  integer               :: num_reals = 0
  logical               :: terminal_flag = .false.
  integer               :: natom_cnt_save
  integer, allocatable  :: use_atm_map(:)
  character(512)        :: char_tmp_512
  integer               :: inerr

#ifdef CUDA
  integer               :: errlen, abortsig
  character, allocatable :: cuda_errmsg(:)
#endif
#ifdef GTI
  double precision, allocatable:: atm_qterm_nonTI(:)
#endif

#ifdef MPI_DEBUGGER
  integer, volatile     :: release_debug

  call mpi_init(err_code_mpi)
  call mpi_comm_rank(mpi_comm_world, release_debug, err_code_mpi)


  ! Lock global master into an infinite loop while release_debug == 0.
  ! This must be changed in a debugger (gdb) via
  ! "set variable release_debug=1" or something
  do
     if (release_debug .ne. 0) exit
  end do

  ! Hold us at a barrier here. GDB (or IDBC, etc.) must be attached to at least
  ! the master process, which takes time. Release everybody once master has been
  ! broken out of the above infinite loop.
  call mpi_barrier(mpi_comm_world, err_code_mpi)
#endif
  call second(run_start_cputime)
  call wall(run_start_walltime)

#ifdef MPI
! Establish pmemd communicators and process ranks
  call setup_groups

  master = mytaskid .eq. 0 ! Make task 0 the master:

! CUDA version with GB can run on 1 GPU / task.
! MPI+OpenMP version can run on 1 MPI rank.
! All other combinations require .gt. 1
#ifndef _OPENMP_
#ifndef CUDA
  if (numtasks .lt. 2 .and. master) then
    write(mdout, *) &
      'MPI version of PMEMD must be used with 2 or more processors!'
    call mexit(6, 1)
  end if
#endif
#endif /*ndef _OPENMP_*/

#else
  master = .true.  ! In the single-threaded version, the 1 process is master
#endif

#ifdef _OPENMP_
if(usemidpoint) then
  inerr = 0
  call get_environment_variable("OMP_NUM_THREADS", char_tmp_512, inerr)
  if (inerr .eq. 0) then
    call omp_set_num_threads(1)
  end if
else ! usemidpoint
! Check whether every node or MPI rank has OMP_NUM_THREADS set otherwise exit.
! This is to make sure the node doesn't get thrashed with excessive default threads.
! TODO: Can process all MPI ranks for the variable and print just once.
! Note: get_environment_variable is part of the F2003 standard but seems
!       to be supported by GNU, Intel, IBM and Portland (2010+) compilers
  inerr = 0
  call get_environment_variable("OMP_NUM_THREADS", char_tmp_512, inerr)
  if (inerr .eq. 0) then
    call get_environment_variable("HOSTNAME", char_tmp_512, inerr)
    if (inerr .eq. 0) then
        write(mdout, *) 'OMP_NUM_THREADS not set on Hostname: Unknown, please refer manual.'
    else
        write(mdout, *) 'OMP_NUM_THREADS not set on ', trim(char_tmp_512), ', please refer manual.'
    end if
    call mexit(6, 1)
  end if
endif ! usemidpoint
#endif /*_OPENMP_*/


# ifdef BINTRAJ
  ! Activate NetCDF interface
  call NC_setupAmberNetcdf(mdout, "pmemd", VERSION_STRING)
# endif

! Reset the stack limits if you can:

  call unlimit_stack(new_stack_limit)

#ifdef TIME_TEST
  call init_test_timers   ! For mpi performance monitoring
  call enable_test_timers
#endif

! Create gpu context
#ifdef CUDA
#ifdef MPI
  call gpu_startup(mytaskid, numtasks, pmemd_comm_number)
#else
  call gpu_startup()
#endif
#endif

  if (master) then
    call master_setup(num_ints, num_reals, new_stack_limit, terminal_flag)
#ifdef CUDA
  else
    !RCW: Call to set device is mostly redundant now due to
    !     removal of -gpu command line argument. But leave for now.
    call gpu_set_device(-1)
    call gpu_init()
#ifdef MPI
    call gpu_send_slave_device_info()
#endif
#endif
  end if

#ifdef MPI
  call bcast_logical(terminal_flag, pmemd_comm)
  if (terminal_flag) call mexit(6,0)
  ! Do generic broadcasts of data used in pretty much all circumstances...
  call bcast_remd_method
  call bcast_mdin_ctrl_dat
#endif

! If a terminal flag was put on the command-line (--help, --version, etc.), then
! just exit here
if (terminal_flag) call mexit(6,0)

! General purpose set up should occur after terminal_flag has been checked.

  if (master) then
#ifndef MPI
#ifndef NOXRAY
    call xray_init()
    call xray_init_clean_up()
#endif
#endif
  end if

! Determine if we use the GPU random number stream or a C implementation
! of the Fortran  CPU random number stream for GPU runs on all threads.
  if ( ig < -1 ) then
#ifdef CUDA
    if (master) write (mdout, '(a)') "Note: ig <= -1 - using CPU random number generator in place of GPU version."
    call gpu_set_cpu_randoms(.true.)
  else
    call gpu_set_cpu_randoms(.false.)
#else
    if (master) write (mdout, '(a)') "Note: ig <= -1 is a debugging option for GPU code only. Ignoring and setting ig = abs(ig)."
#endif  /* CUDA */
  end if
!Either way whatever ig was set to we set it positive here since the random
!number generators expect a positive seed and there is lots of code later for
!remd etc that adds integer values to ig and therefore assumes it is positive.
  ig = abs(ig)

#ifdef MPI
  if ( no_ntt3_sync ) then
    !Here we are not synching the random number generator across threads for
    !NTT=3 but we don't want every thread to have the same random seed. So for
    !the moment just add mytask id to ig.
#ifndef CUDA
    !In the CUDA case this seems to cause massive problems with NTT3 in parallel
    !so
    !only do it for regular CPU runs. NO_NTT3_SYNC shouldn't help much with GPUs
    !anyway since they use their own random number system.
    ig = ig+mytaskid
#endif
  end if
#endif /* MPI */

    if (ifsams .eq. 1) then
      call sams_setup(irest, i, natom)
      !TBD call sams_setup(irest, i)
#ifdef CUDA
      call gpu_init_pbc(pbc_box(1), pbc_box(2), pbc_box(3), pbc_alpha, &
                        pbc_beta, pbc_gamma, uc_volume, uc_sphere, &
                        vdw_cutoff + skinnb, pbc_box, reclng, cut_factor, &
                        ucell, recip)
      call gpu_force_new_neighborlist()
#endif  
      if (i < 0) then
        write(mdout, '(/,a,/)') ' Error in setting up SAMS  Please check '
        call mexit(6, 1)
      endif
    end if
  
#ifdef MPI
  call bcast_amber_prmtop_dat
  call bcast_inpcrd_dat(natom)
  call bcast_constraints_dat(natom, ibelly, ntr)
  call bcast_extra_pnts_nb14_dat
  call bcast_shake(natom)
  call parallel_dat_setup(natom, num_ints, num_reals)
  call alloc_force_mem(natom,num_reals,ips)
  call get_cmdline_bcast
  ! bcast constant pH stuff
if(.not. usemidpoint) then
  if (icnstph .gt. 0 .or. (icnste .gt. 0 .and. cpein_specified)) call cnstph_bcast(num_ints, num_reals)
! PHMD
  call phmd_bcast(err_code_mpi)

  ! bcast constant Redox potential stuff
  if (icnste .gt. 0 .and. .not. cpein_specified) call cnste_bcast(num_ints, num_reals)

  ! bcast Hybrid Solvent REMD stuff
  if (hybridgb .gt. 0) call hybridsolvent_remd_bcast
endif
#ifdef CUDA
  if (ntb .ne. 0) then
    call bcast_pbc
    call bcast_mdin_ewald_dat
    if (.not. master) then
      call gpu_init_pbc(pbc_box(1), pbc_box(2), pbc_box(3), pbc_alpha, &
                        pbc_beta, pbc_gamma, uc_volume, uc_sphere, &
                        vdw_cutoff + skinnb, pbc_box, reclng, cut_factor, &
                        ucell, recip)
    end if
  end if
  if (nmropt .ne. 0) call bcast_nmr_dat
#endif /* CUDA */
  ! Set up TI
  if (icfe .ne. 0) then
    call ti_bcast_dat(natom, ntypes)
  end if
  ! Set up REMD if we're actually performing REMD
  if (remd_method .ne. 0) then
    call setup_pv_correction(mdout, ntp, master)
    if (master) then
      call remd_setup(numexchg)
      ! Now that REMD has been set up, perform any additional setup for
      ! netcdf files (e.g. for multi-REMD need to store indices)
      if (master) call setup_remd_indices
    else
      call slave_remd_setup
    end if
if(.not. usemidpoint) then
    call setup_remd_randgen
endif
  end if
#endif /* MPI */
  ! Set up MC barostat
  if (ntp .gt. 0 .and. barostat .eq. 2) call mcbar_setup(ig)
  ! Set up TI
  if (icfe .ne. 0) then
#ifdef GTI
    ti_lam_sch = gti_lam_sch
#endif
    call ti_change_weights(clambda)
  endif

  ! Set up REAF
#ifdef GTI
  if (ifreaf .ne. 0) then
    reaf_sch=gti_tau_sch
    call read_update_weights
  endif
#endif

if(.not. usemidpoint) then
  ! Set up AMD
  if (iamd .gt. 0) then
    call amd_setup(ntwx)
  endif
  ! Set up GaMD
  if (igamd .gt. 0) then
    call gamd_setup(ntwx)
  endif
  ! Set up scaledMD
  if (scaledMD .gt. 0) then
    call scaledMD_setup(ntwx)
  endif
endif

  ! Set up Nudged Elastic Band

! Initialize system

#ifdef CUDA

#ifdef GTI  /* cheat the gpu code as if this is a non-TI calc. */

    call gpu_setup_system(natom, imin, tol, ntf, ntb, ips, ntp, barostat, ntt, gbsa, &
        vrand, icnstph, icnste, hybridgb, 0, surften, ischeme, ithermostat, ntc)
    call gti_init_md_parameters(vlimit)
    
    if (ti_mode .ne. 0) then
        if (igamd .gt. 0) then
          if(igamd.eq.12.or.igamd.eq.13.or.igamd0.eq.12.or.igamd0.eq.13 .or.&
            (igamd.ge.110 .and. igamd.le.120).or.(igamd0.ge.110 .and. igamd0.le.120).or. &
            (igamd.eq.14 .or. igamd.eq.15).or.(igamd0.eq.14 .or. igamd0.eq.15).or.&
            (igamd.eq.16 .or. igamd.eq.17).or.(igamd0.eq.16 .or. igamd0.eq.17).or.&
            (igamd.eq.18 .or. igamd0.eq.18).or.(igamd.eq.19 .or.igamd0.eq.19).or.&
            (igamd.ge.20 .and. igamd.le.28).or.(igamd0.ge.20.and.igamd0.le.28))then
          call gti_init_ti_parameters(ti_latm_lst, ti_lst, ti_sc_lst, &
            ti_sc_bat_lst, ctrlVar)
           ti_mode=0
          else !for lig-gamd
            call gti_init_gamd_ti_parameters(ti_latm_lst, ti_lst, ti_sc_lst, &
                     ti_sc_bat_lst, clambda, klambda, scalpha, scbeta, scgamma, &
                    gti_cut, gti_add_sc, gti_ele_gauss, gti_ele_sc, &
                    gti_auto_alpha, gti_vdw_sc, gti_vdw_cap)
         endif

        else
          call gti_init_ti_parameters(ti_latm_lst, ti_lst, ti_sc_lst, &
            ti_sc_bat_lst, ctrlVar)
        endif
        i=-1
        if (ifsams .eq. 1 .and. sams_type .eq. 3) i=1
        if (ifmbar .ne. 0) call gti_init_mbar_parameters(bar_states, bar_lambda, i)

        if (gti_heating .ne. 0) call gti_setup_localheating(gti_tempi)
    endif

    if (reaf_mode .ge. 0) call gti_init_reaf_parameters(reaf_mode, reaf_tau, gti_add_re, reaf_atom_list)

#else
  call gpu_setup_system(natom, imin, tol, ntf, ntb, ips, ntp, barostat, ntt, gbsa, &
                        vrand, icnstph, icnste, hybridgb, ti_mode, surften, ischeme, ithermostat, ntc)
  if (icfe .ne. 0) then
    call gpu_setup_softcore_ti(ti_lst_repacked, ti_sc_lst, ti_latm_lst_repacked, &
                               ti_latm_cnt(1))
    call gpu_set_ti_constants(clambda, scalpha, scbeta, ifmbar, bar_intervall, bar_states, &
                              bar_lambda)
  end if
#endif /* GTI */

  call gpu_upload_iphmd(iphmd)
  call gpu_upload_crd(atm_crd)
  call gpu_upload_charges(atm_qterm)
  call gpu_upload_masses(atm_mass)
  call gpu_upload_frc(atm_frc)
  call gpu_upload_vel(atm_vel)
  call gpu_upload_last_vel(atm_last_vel)
  
  call gpu_init_extra_pnts_nb14(gbl_frame_cnt, ep_frames, ep_lcl_crd)
  call gpu_constraints_setup(natc, atm_jrc, atm_weight, atm_xc)

#endif /* CUDA */
! The following call does more uniprocessor setup and mpi master/slave setup.

  if (using_pme_potential) then
    call pme_alltasks_setup(num_ints, num_reals)
  else if (using_gb_potential) then
if(.not. usemidpoint) then
    call gb_alltasks_setup(num_ints, num_reals)
endif
  end if

  ! Set up mcres
  if (mcwat .gt. 0) then
    call setup_mcres(mcresstr)
  endif

  if(ramdint .gt. 0) then
    call init_ramd
  endif

if(.not. usemidpoint) then
! The following does constant pH setup stuff

  if (icnstph .gt. 0 .or. (icnste .gt. 0 .and. cpein_specified)) then
    call cnstph_setup(atm_crd)
  end if

!PHMD
  if (iphmd .gt. 0 .and. master) then
    call phmdzero()
  end if

! The following does constant Redox setup stuff

  if (icnste .gt. 0 .and. .not. cpein_specified) then
    call cnste_setup(atm_crd)
  end if
endif

! Call shake_setup, which will tag those bonds which are part of 3-point
! water molecules and also set up data structures for non-fastwater shake.
! Constraints will be done for waters using a fast analytic routine -- dap.
! Currently, ALL processes call shake_setup, as it is now very fast, and this
! is probably cheaper than broadcasting the results (also, there is now atom
! selection in shake setup under mpi).
  call shake_setup(num_ints, num_reals)

!  if (using_pme_potential) then
!    call final_pme_setup(num_ints, num_reals)
!    if (icnstph .eq. 2 .or. (icnste .eq. 2 .and. cpein_specified)) &
!      call final_gb_setup(natom, num_ints, num_reals, cph_igb)
!    if (icnste .eq. 2 .and. .not. cpein_specified) &
!      call final_gb_setup(natom, num_ints, num_reals, ce_igb)
!  else if (using_gb_potential) then
!    call final_gb_setup(natom, num_ints, num_reals, igb)
!  end if

  if (using_pme_potential) then
    if (icnstph .eq. 2 .or. (icnste .eq. 2 .and. cpein_specified)) then
#ifdef CUDA
      natom_cnt_save = natom
      natom = cphfirst_sol - 1
      call gpu_update_natoms(natom, .false.)
#endif
if(.not. usemidpoint) then
      call final_gb_setup(natom, num_ints, num_reals, cph_igb)
endif
#ifdef CUDA
      natom = natom_cnt_save
      call gpu_update_natoms(natom, .true.)
#endif
    else if (icnste .eq. 2 .and. .not. cpein_specified) then
#ifdef CUDA
      natom_cnt_save = natom
      natom = cefirst_sol - 1
      call gpu_update_natoms(natom, .false.)
#endif
if(.not. usemidpoint) then
      call final_gb_setup(natom, num_ints, num_reals, ce_igb)
endif
#ifdef CUDA
      natom = natom_cnt_save
      call gpu_update_natoms(natom, .true.)
#endif
    else if (hybridgb .gt. 0) then
#ifdef MPI
#ifdef CUDA
      natom_cnt_save = natom
      natom = hybridsolvent_remd_atm_cnt
      call gpu_update_natoms(natom, .false.)
#endif
if(.not. usemidpoint) then
      call final_gb_setup(hybridsolvent_remd_atm_cnt, num_ints, num_reals, hybridgb)
endif
#ifdef CUDA
      natom = natom_cnt_save
      call gpu_update_natoms(natom, .true.)
#endif
#endif /* End #ifdef MPI */
    endif

    call final_pme_setup(num_ints, num_reals)

#if defined(CUDA) && defined(GTI)
    if (ti_mode.ne.0 .or. reaf_mode .ge. 0 ) then
        call gti_setup_update(atm_qterm, gti_syn_mass, atm_mass)
    endif
    if(igamd.eq.12.or.igamd.eq.13.or.igamd0.eq.12.or.igamd0.eq.13 .or.&
            (igamd.ge.110 .and. igamd.le.120).or.(igamd0.ge.110 .and. igamd0.le.120).or. &
            (igamd.eq.14 .or. igamd.eq.15).or.(igamd0.eq.14 .or. igamd0.eq.15).or.&
            (igamd.eq.16 .or. igamd.eq.17).or.(igamd0.eq.16 .or. igamd0.eq.17).or.&
            (igamd.eq.18 .or. igamd0.eq.18).or.(igamd.eq.19 .or.igamd0.eq.19).or.&
            (igamd.ge.20 .and. igamd.le.28).or.(igamd0.ge.20 .and. igamd0.le.28))then
       call gti_setup_update(atm_qterm, gti_syn_mass, atm_mass)
   endif

    call gti_update_simulation_const;
#endif /* GTI */

  else if (using_gb_potential) then
if(.not. usemidpoint) then
#ifdef _OPENMP_
    call final_frc_setup()
    call final_gb_setup_hybrid(natom, atm_crd, num_ints, num_reals, igb)
#else
    call final_gb_setup(natom, num_ints, num_reals, igb)
#endif
endif
  end if

! Deallocate data that is not needed after setup:

  if (using_pme_potential .and. &
           icnstph .ne. 2 .and. icnste .ne. 2 .and. hybridgb .le. 0 .and. &
      emil_do_calc .le. 0) then
if(.not. usemidpoint) then
    if(isgld==0.or.sgtype==1)then
    num_ints = num_ints - size(atm_numex) - size(gbl_natex)
    deallocate(atm_numex, gbl_natex)
    endif
endif
  end if

  if (imin .ne. 0) then
    num_reals = num_reals - size(atm_vel)
    deallocate(atm_vel)
  end if

! Should be all done with setup. Report dynamic memory allocation totals:

  if (master) then
    write(mdout, '(a)')       '| Dynamic Memory, Types Used:'
    write(mdout, '(a,i14)')   '| Reals      ', num_reals
    write(mdout, '(a,i14,/)') '| Integers   ', num_ints
    if (using_pme_potential) then
      write(mdout, '(a,i12,/)') '| Nonbonded Pairs Initial Allocation:', &
                               ipairs_maxsize
    end if
  end if

#ifdef CUDA
#ifdef MPI
    if (master) then
#endif
      allocate(cuda_errmsg(32768))
      call gpu_final_case_inhibitors(cuda_errmsg, errlen, abortsig, atm_crd)
      do i = 1, errlen
        write(mdout, '(a)', ADVANCE="NO") cuda_errmsg(i)
      end do
      deallocate(cuda_errmsg)
      if (abortsig == 1) then
        call mexit(0, 0)
      end if
      call gpu_write_memory_info(mdout)
#ifdef MPI
    end if
#endif
#endif

#ifdef MPI
  if (master) then
    write(mdout, '(a,i4,a,/)') '| Running AMBER/MPI version on ', numtasks, &
                               ' MPI task'
#ifdef _OPENMP_
    write(mdout,'(a,I5)') &
    "| Number of openmp threads per MPI task: ", omp_get_max_threads()
#endif /*_OPENMP_*/
    write(mdout, '(a)') ' '
  end if
#endif
  
if(.not. usemidpoint) then
  ! Set up SGLD
#ifdef MPI
  if (isgld > 0) call psgld(natom,atm_mass,atm_crd,atm_vel, remd_method,gbl_my_atm_lst)
#else /* MPI */
  if (isgld > 0) call psgld(natom,atm_mass,atm_crd,atm_vel)
#endif /* MPI */
  ! Prepare for EMAP constraints
   if (iemap>0) call pemap(dt,temp0)
! Prepare for Isotropic periodic sum of nonbonded interaction
   if (ips .gt. 0)then
     call ipssys(atm_crd)
   endif
endif ! not usemidpoint

  if (master) write(mdout,'(80(''-'')/''   4.  RESULTS'',/80(''-'')/)')
#ifndef CUDA
  if (using_pme_potential) then
    call chk_switch(max_erfc_relerr)
    call chk_ef_spline_tables(max_erfc_relerr)
  end if
#endif

! We reset the random number generator here to be consistent with sander 10.
! I think this is pretty much unnecessary in pmemd, but we will only get
! sander-consistent results if we do it.

  call amrset(ig + 1)

#ifndef NO_F95_RANDOM_NUMBER_INTRINSIC
  if (irandom == 2) then
    !We will use the F95 RANDOM_NUMBER intrinsic for uniform random numbers
    call amrset_f95_intrinsic()
  end if
#endif

#ifdef MPI
if(usemidpoint) then
  call setup_processor(pbc_box, bkt_size)
  call alloc_full_list_to_proc(atm_crd, atm_frc, atm_vel, atm_last_vel, natom,pbc_box)
! call proc_data_dump(pbc_box,proc_num_atms)
 ! call bonds_midpoint_setup !called in pme_direct
  call alloc_nb_mask_data
  call alloc_atm_rel_crd
  call make_atm_excl_mask_list(natom, atm_numex, gbl_natex)
  call pme_recip_setup()
  allocate(use_atm_map(natom))
  use_atm_map(:) = 0
  call nb14_setup_midpoint(num_ints, num_reals, use_atm_map)
  call proc_alloc_force_mem
! call make_nb_adjust_pairlst_midpoint(gbl_nb_adjust_pairlst)
endif ! usemidpoint

! If we're running replica exchange, we need to reset our random number
! generator to match the random streams that sander produces so we can
! reproduce those results here

  if (remd_method .ne. 0) call amrset(ig + 17 * (repnum-1))

  if (irandom == 2) then
    !We will use the F95 RANDOM_NUMBER intrinsic for uniform random numbers
    call amrset_f95_intrinsic()
  end if

! In this implementation of pmemd, when running molecular dynamics,
! parallelization occurs at the level of subroutine runmd(), and when running
! minimizations, parallelization occurs at the level of subroutine force().

#ifndef CUDA
if(.not. usemidpoint) then
  if (nmropt .ne. 0) call bcast_nmr_dat
endif
#endif

  call second(run_setup_end_cputime)
  call wall(run_setup_end_walltime)

  ! Parallelization of minimization, nonmaster nodes:

  if (imin .ne. 0) then
    if (.not. master) then         ! All nodes do only force calc
if(.not. usemidpoint) then
      call runmin_slave(natom, atm_crd, atm_frc, gbl_img_atm_map)
endif
      call second(run_end_cputime)
#ifdef _OPENMP_
      call wall(run_end_walltime)
#endif
      call profile_cpu(imin, igb, loadbal_verbose)

#ifdef TIME_TEST
      call print_test_timers    ! Debugging output of performance timings
#endif
      call free_comms
      call mexit(0, 0)
    endif
  end if
#else
  call second(run_setup_end_cputime)
  call wall(run_setup_end_walltime)
#endif /* end MPI */

  ! Initialize the printing of ongoing time and performance summaries. We call this
  ! here after all the setup is done so we don't end up including all the startup
  ! time etc.
  if (master) call print_ongoing_time_summary(0,0,0.0d0,0)

! Now do the dynamics or minimization:

  if (imin .eq. 0) then         ! Do molecular dynamics:
if(.not. usemidpoint) then
!-------NFE MOD----------------------------------------------------------
    if (infe.gt.0) then
#ifdef MPI
        call nfe_on_pmemd_init(atm_igraph,atm_mass,atm_crd, remd_method)
#else
        call nfe_on_pmemd_init(atm_igraph,atm_mass,atm_crd, 0)
#endif /* MPI */
    endif
!-------NFE MOD--------- ------------------------------------------------
!DG
#ifdef MPI
    if(ineb>0) then
      call neb_init()
#ifdef CUDA
      call nebread()
      call gpu_setup_shuttle_info(pcnt, 0, partial_mask)
      call gpu_upload_neb_info(neb_nbead, beadid, nattgtrms, nattgtfit, &
           last_neb_atom, skmin, skmax, tmode, rmsmask, fitmask,idx_mask)
      call gpu_setup_neb()
#else
        if (master) then
          call nebread()
        end if
        call mpi_bcast(nattgtrms, 1, MPI_INTEGER, 0, mpi_comm_world, err_code_mpi)
        call mpi_bcast(rmsmask, nattgtrms, MPI_INTEGER, 0, mpi_comm_world, err_code_mpi)
#endif
    end if
#endif
!DG

endif ! not usemidpoint

    if(reweight .eq. 0) then
#ifdef MPI
    call runmd(natom, atm_crd, atm_mass, atm_frc, atm_vel, atm_last_vel, &
               gbl_my_atm_lst, remd_method, numexchg)
#else
    call runmd(natom, atm_crd, atm_mass, atm_frc, atm_vel, atm_last_vel)
#endif /* MPI */
    else
#ifdef MPI
        call runreweight(natom, atm_crd, atm_mass, atm_frc, atm_vel, atm_last_vel, &
                   gbl_my_atm_lst, remd_method, numexchg)
#else
        call runreweight(natom, atm_crd, atm_mass, atm_frc, atm_vel, atm_last_vel)
#endif /* MPI */
    end if
  else                          ! Do minimization:
if(.not. usemidpoint) then
    call runmin_master(natom, atm_crd, atm_frc, atm_igraph)
endif
    ! Write restart file. atm_vel is not actually used.
    call write_restart(restrt, natom, atm_crd, atm_vel, &
                       0.d0, .true., restrt_name)
    close (restrt)
  endif

if(.not. usemidpoint) then
  ! finish up EMAP
  if (iemap>0) call qemap()

  ! Finish up external interface
  if (master .and. iextpot .ne. 0) then
    call external_cleanup()
  end if

  ! Deallocate constant pH arrays
  if (icnstph .ne. 0 .or. (icnste .ne. 0 .and. cpein_specified)) then
    call cnstph_deallocate()
  end if
  ! Deallocate constant redox arrays
  if (icnste .ne. 0 .and. .not. cpein_specified) then
    call cnste_deallocate()
  end if

  ! finish up GB if it was used
  if (using_pme_potential .and. (icnstph .eq. 2 .or. icnste .eq. 2 .or. hybridgb .gt. 0)) then
    if (icnstph .eq. 2) then
      call gb_cleanup(num_ints, num_reals, cph_igb)
    else if (icnste .eq. 2) then
      call gb_cleanup(num_ints, num_reals, ce_igb)
    else if (hybridgb .gt. 0) then
      call gb_cleanup(num_ints, num_reals, hybridgb)
    end if
  else if (using_gb_potential) then
#ifdef _OPENMP_
    call gb_cleanup_hybrid(num_ints, num_reals, igb)
#else
    call gb_cleanup(num_ints, num_reals, igb)
#endif
  end if
endif ! not usemidpoint

#ifdef MPI
! If doing minimization, set and broadcast notdone to inform other nodes that
! we are finished calling force()

  if (imin .ne. 0) then
    notdone = 0
    call mpi_bcast(notdone, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)
  endif

  call dealloc_force_mem(ips)

  if (imin .eq. 0) then
    call remd_cleanup
  end if


  if (master) then
    if (ioutfm .eq. 1) then
      call close_binary_files
    else
      if (ntwx .gt. 0) close(mdcrd)
      if (ntwv .gt. 0) close(mdvel)
      if (ntwf .gt. 0) close(mdfrc)
    end if
    if (ntwe .gt. 0) close(mden)
  else

    call second(run_end_cputime)
#ifdef _OPENMP_
    call wall(run_end_walltime)
#endif
    call profile_cpu(imin, igb, loadbal_verbose)

#ifdef TIME_TEST
    call print_test_timers      ! Debugging output of performance timings
#endif
#ifdef CUDA
! Shut down GPU
  call gpu_shutdown()
#endif
    call free_comms
    call mexit(0, 0)
  endif
#endif /* end MPI */

if(.not. usemidpoint) then
  !--------------------- NFE_MOD----------------------
  if (infe.gt.0) then
    call nfe_on_pmemd_exit()
  endif
  !---------------------------------------------------
!DG
#ifdef MPI
  if(ineb>0) call neb_finalize()
#endif
!DG
endif ! not usemidpoint

  ! Master prints timings (averages for mpi) in mdout:

  call second(run_end_cputime)
  call wall(run_end_walltime)

  write(mdout,'(80(1H-)/,''   5.  TIMINGS'',/80(1H-)/)')
#ifdef MPI
  call profile_cpu(imin, igb, loadbal_verbose)
#else
  call profile_cpu(imin, igb, 0)
#endif

#ifdef MPI
  !Write final performance numbers to mdout
  if (imin .eq. 0 .and. remd_method .eq. 0) then
    call print_ongoing_time_summary(nstlim, nstlim, dt, mdout)
  else if (imin .eq. 0) then
    call print_ongoing_time_summary(nstlim*numexchg, nstlim*numexchg, dt, mdout)
  end if
#else
  !Write final performance numbers to mdout
  if (imin .eq. 0) call print_ongoing_time_summary(nstlim,nstlim,dt,mdout)
#endif

#ifdef MPI
  write(mdout, '(/, a, f11.2, a)') &
        '|  Master Setup CPU time:     ', &
        run_setup_end_cputime - run_start_cputime, ' seconds'
  write(mdout, '(   a, f11.2, a)') &
        '|  Master NonSetup CPU time:  ', &
        run_end_cputime - run_setup_end_cputime, ' seconds'
  write(mdout, '(   a, f11.2, a, f9.2, a)') &
        '|  Master Total CPU time:     ', &
        run_end_cputime - run_start_cputime, ' seconds', &
        (run_end_cputime - run_start_cputime) / 3600.d0, ' hours'
  write(mdout, '(/, a, i9, a)') &
        '|  Master Setup wall time:   ', &
        run_setup_end_walltime - run_start_walltime, '    seconds'
  write(mdout, '(   a, i9, a)') &
        '|  Master NonSetup wall time:', &
        run_end_walltime - run_setup_end_walltime, '    seconds'
  write(mdout, '(   a, i9, a, f9.2, a)') &
        '|  Master Total wall time:   ', &
        run_end_walltime - run_start_walltime, '    seconds', &
        dble(run_end_walltime - run_start_walltime) / 3600.d0, ' hours'
#else
  write(mdout, '(/, a, f11.2, a)') &
        '|  Setup CPU time:     ', &
        run_setup_end_cputime - run_start_cputime, ' seconds'
  write(mdout, '(   a, f11.2, a)') &
        '|  NonSetup CPU time:  ', &
        run_end_cputime - run_setup_end_cputime, ' seconds'
  write(mdout, '(   a, f11.2, a, f9.2, a)') &
        '|  Total CPU time:     ', &
        run_end_cputime - run_start_cputime, ' seconds', &
        (run_end_cputime - run_start_cputime) / 3600.d0, ' hours'
  write(mdout, '(/, a, i9, a)') &
        '|  Setup wall time:   ', &
        run_setup_end_walltime - run_start_walltime, '    seconds'
  write(mdout, '(   a, i9, a)') &
        '|  NonSetup wall time:', &
        run_end_walltime - run_setup_end_walltime, '    seconds'
  write(mdout, '(   a, i9, a, f9.2, a)') &
        '|  Total wall time:   ', &
        run_end_walltime - run_start_walltime, '    seconds', &
        dble(run_end_walltime - run_start_walltime) / 3600.d0, ' hours'
#endif

  !close(mdout)  DGDG
#ifndef MPI
#ifndef NOXRAY
  call xray_fini()
#endif
#endif
#ifdef CUDA
! Shut down GPU
  call gpu_shutdown()
#endif

#ifdef TIME_TEST
  call print_test_timers        ! Debugging output of performance timings
#endif
  close(mdout) !DGDG
#ifndef MPI
#ifndef NOXRAY
  if (ntwsf .gt. 0) then
    close(sfout)
  end if
#endif
#endif
#ifdef MPI
  call free_comms               ! free MPI communicators
#endif
  call mexit(6, 0)

end program pmemd
