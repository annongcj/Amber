#include "copyright.i"
#include "dbg_arrays.i"

!*******************************************************************************
!
! Module: pme_force_mod
!
! Description: <TBS>
!
!*******************************************************************************

module pme_force_mod

  use gbl_datatypes_mod
  use prmtop_dat_mod, only : numextra
  use energy_records_mod, only : pme_pot_ene_rec, pme_pot_ene_rec_size, &
             null_pme_pot_ene_rec, pme_virial_rec, pme_virial_rec_size, &
             null_pme_virial_rec
  use extra_pnts_nb14_mod
#ifdef TIME_TEST
  use timers_mod
#endif
#ifdef GTI
  use gti_mod
  use reaf_mod
#endif
  use external_mod

  implicit none

  double precision, allocatable :: img_frc(:,:)
  double precision, allocatable :: nb_frc(:,:)
  double precision, allocatable :: ips_frc(:,:)
  integer, allocatable, save            :: gbl_nvdwcls(:)

  ! The following variables don't need to be broadcast:

  integer, save         :: irespa = 0

  ! Molecular virial correction factor; used internally.

  double precision, save, private       :: molvir_netfrc_corr(3, 3)

contains

!*******************************************************************************
!
! Subroutine:  alloc_pme_force_mem
!
! Description: <TBS>
!
!*******************************************************************************

subroutine alloc_pme_force_mem(ntypes, num_ints, num_reals)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer                       :: ntypes

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed

  allocate(gbl_nvdwcls(ntypes), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_ints = num_ints + size(gbl_nvdwcls)

  gbl_nvdwcls(:) = 0

  return

end subroutine alloc_pme_force_mem

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  bcast_pme_force_dat
!
! Description: <TBS>
!
!*******************************************************************************

subroutine bcast_pme_force_dat(ntypes)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer       :: ntypes

! Local variables:

  integer       :: num_ints, num_reals  ! returned values discarded

  ! Nothing to broadcast.  We just allocate storage in the non-master nodes.

  if (.not. master) then
    num_ints = 0
    num_reals = 0
    call alloc_pme_force_mem(ntypes, num_ints, num_reals)
  end if

  ! The allocated data is not initialized from the master node.

  return

end subroutine bcast_pme_force_dat

!*******************************************************************************
!
! Subroutine:  alloc_force_mem
!
! Description: <TBS>
!
!*******************************************************************************

subroutine alloc_force_mem(atm_cnt,num_reals,ips)

  use pmemd_lib_mod
  use processor_mod
  implicit none

! Formal arguments:

  integer                       :: atm_cnt

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_reals

  integer                       :: ips

! Local variables:

  integer               :: alloc_failed

  if (ips .gt. 0) then
   allocate(img_frc(3, atm_cnt), &
            nb_frc(3, atm_cnt), &
            ips_frc(3, atm_cnt), &
            stat = alloc_failed)

   if (alloc_failed .ne. 0) call setup_alloc_error

   num_reals = num_reals + size(img_frc) + size(nb_frc) + size(ips_frc)
  else
   allocate(img_frc(3, atm_cnt), &
            nb_frc(3, atm_cnt), &
            stat = alloc_failed)

   if (alloc_failed .ne. 0) call setup_alloc_error

   num_reals = num_reals + size(img_frc) + size(nb_frc)
  end if

  return

end subroutine alloc_force_mem

!*******************************************************************************
!
! Subroutine:  dealloc_force_mem
!
! Description: <TBS>
!
!*******************************************************************************

subroutine dealloc_force_mem(ips)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer                       :: ips

  if (ips .gt. 0) then
    deallocate(img_frc, nb_frc, ips_frc)
  else
    deallocate(img_frc, nb_frc)
  end if

  return

end subroutine dealloc_force_mem
#endif

!*******************************************************************************
!
! Subroutine:  pme_force
!
! Description: <TBS>
!
!*******************************************************************************

#ifdef MPI
subroutine pme_force(atm_cnt, crd, frc, img_atm_map, atm_img_map, &
                     my_atm_lst, new_list, need_pot_enes, need_virials, &
                     pot_ene, nstep, virial, ekcmt, pme_err_est)

  use bonds_midpoint_mod
  use mol_list_mod, only : gbl_mol_cnt
  use angles_mod
  use angles_ub_mod
  use bonds_mod
  use constraints_mod
  use dihedrals_mod
  use dihedrals_imp_mod
  use extra_pnts_nb14_mod, only : ep_lcl_crd
  use cmap_mod

  use dynamics_dat_mod
  use dynamics_mod
  use ene_frc_splines_mod
  use emap_mod,only:emapforce
  use nb_exclusions_mod
  use pme_direct_mod
  use pme_recip_dat_mod
  use pme_blk_recip_mod
  use pme_slab_recip_mod
  use loadbal_mod
  use img_mod
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use nb_pairlist_mod
  use nmr_calls_mod
#ifndef MPI
#ifndef NOXRAY
  use xray_interface_module, only: xray_get_derivative
  use xray_globals_module, only: xray_active
#endif
#endif
  use parallel_dat_mod
  use parallel_mod
  use pbc_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
  use runfiles_mod
  use timers_mod
  use charmm_mod, only : charmm_active
  use mdin_debugf_dat_mod, only : do_charmm_dump_gold
  use nbips_mod
  use amd_mod
  use gamd_mod
  use scaledMD_mod
  use ti_mod
#ifdef TIDECOMP
  use ti_decomp_mod
#endif
  use phmd_mod
  ! This #ifdef section can be merged with the section on top
  use processor_mod !, only : proc_atm_crd, proc_atm_frc, proc_atm_qterm, &
  !       proc_num_atms
  use parallel_processor_mod
  use pme_recip_midpoint_mod
  use ensure_alloc_mod
#ifdef _OPENMP_
  use omp_lib
#endif /* _OPENMP_ */
  use nfe_setup_mod, only : nfe_on_force => on_force
  use nfe_lib_mod
  use neb_mod
  use nebread_mod
  use multipmemd_mod

#ifdef DBG_ARRAYS
  use dbg_arrays_mod
#endif


  implicit none

  ! Formal arguments:

  integer                       :: atm_cnt
  double precision              :: crd(3, atm_cnt)
  double precision              :: frc(3, atm_cnt)
  integer                       :: img_atm_map(atm_cnt)
  integer                       :: atm_img_map(atm_cnt)
  integer                       :: my_atm_lst(*)
  logical                       :: new_list
  logical                       :: need_pot_enes
  logical                       :: need_virials
  type(pme_pot_ene_rec)         :: pot_ene
  double precision              :: virial(3)            ! Only used for MD
  double precision              :: ekcmt(3)             ! Only used for MD
  double precision              :: pme_err_est          ! Only used for MD
  integer, intent(in)           :: nstep

  ! Local variables:

  type(pme_virial_rec)          :: proc_vir
  integer                       :: offset_loc
  logical                       :: send_data
  integer                       :: iatm, jatm, task_id
  type(pme_virial_rec)          :: vir
  double precision              :: enmr(3), entr
#ifndef MPI
#ifndef NOXRAY
  double precision              :: xray_e
#endif
#endif
  double precision              :: vir_vs_ene
  integer                       :: atm_lst_idx
  integer                       :: alloc_failed
  integer                       :: i, j, k
  logical                       :: params_may_change
  logical                       :: onstep
  double precision              :: net_frcs(3)

  double precision              :: evdwex, eelex
  double precision              :: enemap
  double precision              :: temp_tot_dih
  double precision              :: totdih
  integer                       :: buf_size

! GaMD Local Variables
  double precision              :: pot_ene_nb, pot_ene_lig, pot_ene_dih, fwgt
  double precision, allocatable :: nb_frc_gamd(:,:)
  double precision, allocatable :: nb_frc_tmp(:,:)
  double precision              :: net_frc_ti(3)

  double precision              :: bias_frc(3,atm_cnt), temp_frc(3,atm_cnt)

  double precision              :: temp_holding

  call zero_time()
  call zero_pme_time()

  if (.not. usemidpoint) then
    pot_ene = null_pme_pot_ene_rec
    if (ti_mode .ne. 0) call ti_zero_arrays
  end if

#ifdef CUDA
  if (nmropt .ne. 0) then
    call nmr_weight(atm_cnt, crd, 6)
  end if
  params_may_change = .false.
  if (infe .ne. 0) params_may_change = .true.


if (iskip_ff .ne. 1) then

  if (use_pme .ne. 0) then
    call self(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
              params_may_change)
    if (vdwmeth .eq. 1) then
      call vdw_correction(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                          params_may_change)
    end if
    call update_time(pme_misc_timer)
    if (need_pot_enes) then ! DO NOT TOUCH THIS, SEE RUNMD.F90 TO SET need_pot_enes
      if ((iamd .eq. 2) .or. (iamd .eq. 3)) then
        totdih = 0.0
        if (num_amd_lag .eq. 0) then
          temp_tot_dih = 0.0
          call gpu_calculate_amd_dihedral_energy(totdih)
          buf_size = 1
          call mpi_allreduce(totdih, temp_tot_dih, &
                             buf_size, mpi_double_precision, &
                             mpi_sum, pmemd_comm, err_code_mpi)
          totdih = temp_tot_dih
        end if
        call gpu_calculate_amd_dihedral_weight(totdih)
      end if
      if ( igamd.eq.2 .or. igamd.eq.3 .or. igamd.eq.5 .or. igamd.eq.7 .or. igamd.eq.9 .or. &
                  igamd.eq.11) then
        totdih = 0.0
        if (num_gamd_lag .eq. 0) then
          temp_tot_dih = 0.0
          call gpu_calculate_gamd_dihedral_energy(totdih)
          buf_size = 1
          call mpi_allreduce(totdih, temp_tot_dih, &
                             buf_size, mpi_double_precision, &
                             mpi_sum, pmemd_comm, err_code_mpi)
          totdih = temp_tot_dih
        end if
        call gpu_calculate_gamd_dihedral_weight(totdih)
      end if

    end if

#ifdef GTI

    !! MPI part of GTI routines
    call gti_build_nl_list(ti_mode, lj1264, plj1264, reaf_mode) ! C4PairwiseCUDA
    call gti_clear(need_pot_enes, ti_mode, crd)
    call gti_ele_recip(need_pot_enes, need_virials, uc_volume, ti_mode, reaf_mode)
    call gti_nb(need_pot_enes, need_virials, ti_mode, lj1264, plj1264, reaf_mode) ! C4PairwiseCUDA
    call gti_bonded(need_pot_enes, need_virials, ti_mode, reaf_mode)
    call gti_others(need_pot_enes, nstep, dt, uc_volume, ew_coeff)
    call gti_finalize_force_virial(numextra, need_virials, ti_mode, netfrc, frc)
    if (ineb>0) then
#ifdef TIME_TEST
        call start_test_timer(14, 'neb_exchange_crd', 0)
#endif
      call transfer_fit_neb_crd(nstep)
#ifdef TIME_TEST
        call stop_test_timer(14)
#endif
    endif

    if (need_pot_enes .or. need_virials) call gti_update_md_ene(pot_ene, enmr, virial, ekcmt, ineb, nebfreq, nstep)

    if (ti_mode .gt. 0) then
      if (need_pot_enes) then
        call gti_update_ti_pot_energy_from_gpu(pot_ene, enmr, ti_mode, gti_cpu_output)
        if (ifmbar .ne. 0) call gti_update_mbar_from_gpu
      end if
      if (need_virials) call gti_get_virial(ti_weights(1:2), virial, ekcmt, ti_mode)
    end if

#else /* GTI  */

    if (need_pot_enes) then
      if (ineb>0) then
        call transfer_fit_neb_crd(nstep)
      endif
      if (ti_mode .eq. 0) then
        call gpu_pme_ene(ew_coeff, uc_volume, pot_ene, enmr, virial, &
                         ekcmt, nstep, dt, crd, frc, ineb, nebfreq)
      else
        call gpu_ti_pme_ene(ew_coeff, uc_volume, pot_ene, enmr, virial, &
                            ekcmt, nstep, dt, crd, frc, bar_cont)
      end if
    else
      call gpu_pme_force(ew_coeff, uc_volume, virial, ekcmt, nstep, dt, crd, frc)
    end if

#endif /* GTI  */

    call update_time(nonbond_time)
    if (need_virials) then
      vir%molecular(1,1) = virial(1)
      vir%molecular(2,2) = virial(2)
      vir%molecular(3,3) = virial(3)
    end if
    call dist_enes_virs_netfrcs(pot_ene, vir, ekcmt, net_frcs, &
                                need_pot_enes, need_virials)
    if (need_virials) then
      virial(1) = vir%molecular(1,1)
      virial(2) = vir%molecular(2,2)
      virial(3) = vir%molecular(3,3)
    end if

  else !use_pme .eq. 0 below
    call update_time(pme_misc_timer)
    call ipsupdate(ntb)
    if (need_pot_enes) then ! DO NOT TOUCH THIS, SEE RUNMD.F90 TO SET need_pot_enes
      if ((iamd .eq. 2) .or. (iamd .eq. 3)) then
        totdih = 0.0
        if (num_amd_lag .eq. 0) then
          temp_tot_dih = 0.0
          call gpu_calculate_amd_dihedral_energy(totdih)
          buf_size = 1
          call mpi_allreduce(totdih, temp_tot_dih, &
                             buf_size, mpi_double_precision, &
                             mpi_sum, pmemd_comm, err_code_mpi)
          totdih = temp_tot_dih
        end if
        call gpu_calculate_amd_dihedral_weight(totdih)
      end if
      if ( igamd.eq.2 .or. igamd.eq.3 .or. igamd.eq.5 .or. igamd.eq.7 .or. igamd.eq.9 .or. &
                  igamd.eq.11) then
        totdih = 0.0
        if (num_gamd_lag .eq. 0) then
          temp_tot_dih = 0.0
          call gpu_calculate_gamd_dihedral_energy(totdih)
          buf_size = 1
          call mpi_allreduce(totdih, temp_tot_dih, &
                             buf_size, mpi_double_precision, &
                             mpi_sum, pmemd_comm, err_code_mpi)
          totdih = temp_tot_dih
        end if
        call gpu_calculate_gamd_dihedral_weight(totdih)
      end if
      call gpu_ips_ene(uc_volume, pot_ene, enmr, virial, ekcmt)
      call update_time(nonbond_time)
      if (need_virials) then
        vir%molecular(1,1) = virial(1)
        vir%molecular(2,2) = virial(2)
        vir%molecular(3,3) = virial(3)
      end if
      call dist_enes_virs_netfrcs(pot_ene, vir, ekcmt, net_frcs, &
                                  need_pot_enes, need_virials)
      if (need_virials) then
        virial(1) = vir%molecular(1,1)
        virial(2) = vir%molecular(2,2)
        virial(3) = vir%molecular(3,3)
      end if
    else
      call gpu_ips_force(uc_volume, virial, ekcmt)
      call update_time(nonbond_time)

      if (need_virials) then
        vir%molecular(1,1) = virial(1)
        vir%molecular(2,2) = virial(2)
        vir%molecular(3,3) = virial(3)
        call dist_enes_virs_netfrcs(pot_ene, vir, ekcmt, net_frcs, &
                                    need_pot_enes, need_virials)
        virial(1) = vir%molecular(1,1)
        virial(2) = vir%molecular(2,2)
        virial(3) = vir%molecular(3,3)
      end if
    end if
  end if

  if (nmropt .ne. 0) then
    call nmr_calc(crd, frc, enmr, 6)
  end if

  !AMD calculate weight and scale forces
  if (iamd .gt. 0) then
    call gpu_calculate_and_apply_amd_weights(pot_ene%total, pot_ene%dihedral, &
                                             pot_ene%amd_boost,num_amd_lag)
    ! Update total energy
    pot_ene%total = pot_ene%total + pot_ene%amd_boost
  endif
  !GaMD calculate weight and scale forces
  if(igamd.gt.0)then
    if(igamd.eq.1 .or. igamd.eq.2 .or. igamd.eq.3) then
      call gpu_calculate_and_apply_gamd_weights(pot_ene%total, pot_ene%dihedral, &
      pot_ene%gamd_boost,num_gamd_lag)
       ! Update total energy
       pot_ene%total = pot_ene%total + pot_ene%gamd_boost
     else if(igamd.eq.4 .or. igamd.eq.5 ) then
      pot_ene_nb = pot_ene%elec_tot + pot_ene%vdw_tot + pot_ene%vdw_14 + pot_ene%elec_14 + pot_ene%efield
      call gpu_calculate_and_apply_gamd_weights_nb(pot_ene_nb, pot_ene%dihedral, &
                                                   pot_ene%gamd_boost,num_gamd_lag)
    ! Update total energy
    pot_ene%total = pot_ene%total + pot_ene%gamd_boost
     else if( (igamd.eq.6 .or. igamd.eq.7) .and. (ti_mode.ne.0)) then
       pot_ene_lig = ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                     ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                     ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
       pot_ene_dih = pot_ene%total + &
                      ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene)
       call gpu_calculate_and_apply_gamd_weights_ti_region(pot_ene_lig, pot_ene_dih, &
                                               pot_ene%gamd_boost,fwgt,num_gamd_lag)
       ! Update total energy
       ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
     else if( (igamd.eq.8 .or. igamd.eq.9) .and. (ti_mode.ne.0)) then
       pot_ene_lig = ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene)
       pot_ene_dih = pot_ene%total + &
                     ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene) + &
                     ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                     ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
       call gpu_calculate_and_apply_gamd_weights_ti_region(pot_ene_lig, pot_ene_dih, &
                                               pot_ene%gamd_boost,fwgt,num_gamd_lag)
       ! Update total energy
       ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
     else if( (igamd.eq.10 .or. igamd.eq.11) .and. (ti_mode.ne.0)) then
        pot_ene_lig = ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                      ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
        pot_ene_dih = pot_ene%total + &
                      ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                      ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene)
       call gpu_calculate_and_apply_gamd_weights_ti_region(pot_ene_lig, pot_ene_dih, &
                                               pot_ene%gamd_boost,fwgt,num_gamd_lag)
       ! Update total energy
       ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
     end if
  end if

  !scaledMD calculate scale forces
  if (scaledMD .gt. 0) then
    call gpu_scaledmd_scale_frc(pot_ene%total)
    pot_ene%total = pot_ene%total * scaledMD_lambda
  end if

  call update_time(fcve_dist_time)

  if (infe .gt. 0 .or. iextpot .gt. 0) then
      if (infe .gt.0 .and. nfe_first) then
         call mpi_bcast(nfe_atm_cnt, 1, MPI_INTEGER, 0, pmemd_comm, err_code_mpi)
         if (mytaskid .ne. 0) allocate(nfe_atm_lst(nfe_atm_cnt))
         call mpi_bcast(nfe_atm_lst, size(nfe_atm_lst), MPI_INTEGER, 0, pmemd_comm, err_code_mpi)
         call gpu_setup_shuttle_info(nfe_atm_cnt, 0, nfe_atm_lst)
         nfe_first = .false.
      end if

      call gpu_shuttle_retrieve_data(crd, 0)
      call gpu_shuttle_retrieve_data(frc, 1)

      bias_frc(:,:) = 0.d0

      if (infe .gt. 0) then
         call nfe_on_force(crd, bias_frc, pot_ene%nfe)
      end if

      ! If calling an External library: may have 
      ! nothing to do with nfe, but can take advantage of same 
      ! logic/bandwidth for moving data on/off the card(s).
      if (iextpot .gt. 0) then
         temp_holding = 0.0
         call pme_external(crd, bias_frc, temp_holding)
      endif


      if (mytaskid .eq. 0) then
         do i=1,nfe_atm_cnt
            j=nfe_atm_lst(i)

            frc(1,j) = frc(1,j)+bias_frc(1,j)
            frc(2,j) = frc(2,j)+bias_frc(2,j)
            frc(3,j) = frc(3,j)+bias_frc(3,j)
         enddo
      end if

      call mpi_bcast(frc, 3*atm_cnt, MPI_DOUBLE_PRECISION, 0, pmemd_comm, err_code_mpi)

      if (infe .gt. 0) then
         pot_ene%total     = pot_ene%total + pot_ene%nfe
         pot_ene%restraint = pot_ene%restraint + pot_ene%nfe
      endif 
      if (iextpot .gt. 0) then
         pot_ene%total     = pot_ene%total + temp_holding
         pot_ene%restraint = pot_ene%restraint + temp_holding
      endif

      call gpu_shuttle_post_data(frc, 1)
  end if

  if (ineb .gt. 0) then
#ifdef TIME_TEST
        call start_test_timer(15, 'neb_force', 0)
#endif
    call full_neb_forces(nstep)
#ifdef TIME_TEST
        call stop_test_timer(15)
#endif

    if (beadid==1 .or. beadid==neb_nbead) then
      call gpu_clear_forces()
    endif
  end if

! Ending if of iskip_ff
end if

#else /*Not cuda version below*/

  if (.not. usemidpoint) then
    call do_load_balancing(new_list, atm_cnt)

    ! Do weight changes, if requested.

    if (nmropt .ne. 0) call nmr_weight(atm_cnt, crd, mdout)
  end if

  ! If no ene / force calcs are to be done, clear everything and bag out now.
  if (ntf .eq. 8) then
    pot_ene = null_pme_pot_ene_rec
    virial(:) = 0.d0
    ekcmt(:) = 0.d0
    pme_err_est = 0.d0
    frc(:,:) = 0.d0
    return
  end if

  if (.not. usemidpoint) then
    ! GaMD
    if (igamd .eq. 4 .or. igamd .eq. 5) then
      allocate(nb_frc_gamd(3, atm_cnt), &
               nb_frc_tmp(3, atm_cnt), &
               stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
    end if
  end if ! not usemidpoint

  ! Calculate the non-bonded contributions:

  ! Direct part of ewald plus vdw, hbond, pairlist setup and image claiming:

  if (iphmd .eq. 3) then
    call runphmd(crd,frc,atm_iac,typ_ico,atm_numex,gbl_natex,belly_atm_cnt)
  end if
  if (ntp .gt. 0) call fill_tranvec(gbl_tranvec)
  ! The following encapsulation (save_imgcrds) seems to be necessary to
  ! prevent an optimization bug with the SGI f90 compiler.  Sigh...

  if (new_list) then

    if (usemidpoint) then
      call pme_list_midpoint(atm_cnt, crd, atm_nb_maskdata, atm_nb_mask, &
                             proc_atm_nb_maskdata, proc_atm_nb_mask,nstep)
    else ! usemidpoint
      call pme_list(atm_cnt, crd, atm_nb_maskdata, atm_nb_mask)
      call start_loadbal_timer

      call save_imgcrds(gbl_img_crd, gbl_saved_imgcrd, &
                        gbl_used_img_cnt, gbl_used_img_lst)
    end if ! usemidpoint
  else
    if (.not. usemidpoint) then
      call start_loadbal_timer
      call adjust_imgcrds(atm_cnt, gbl_img_crd, img_atm_map, &
                          gbl_used_img_cnt, gbl_used_img_lst, &
                          gbl_saved_imgcrd, crd, gbl_atm_saved_crd, ntp)
    end if ! not usemidpoint
  end if

  ! Zero energies that are stack or call parameters:

  pot_ene = null_pme_pot_ene_rec
  vir = null_pme_virial_rec

  virial(:) = 0.d0
  ekcmt(:) = 0.d0
  pme_err_est = 0.d0

  ! Zero internal energies, virials, etc.

  enmr(:) = 0.d0
  vir_vs_ene = 0.d0

  net_frcs(:) = 0.d0
  molvir_netfrc_corr(:,:) = 0.d0

! If calling an External library
if (iextpot .gt. 0) then
  do j = 1, gbl_used_img_cnt
    i = gbl_used_img_lst(j)
    img_frc(:, i) = 0.d0
  end do
  do i = 1, atm_cnt
    frc(:,i) = 0.d0
  end do
  virial(:) = 0.d0
  
#ifdef MPI
  if (mytaskid.eq.0) then
#endif /* MPI */
 
!!NEED TO REWORK LOGIC TO ALLOW RESTRAINTS TO BE CALLED EVEN IF PMEMD IS NOT CALLED
!!With below code commented-in, then restraints are called twice which is not good.
#if 0
      
      if (nmropt .ne. 0) then
         call nmr_calc(crd, temp_frc, enmr, 6)
         pot_ene%restraint = pot_ene%restraint + enmr(1) + enmr(2) + enmr(3)
         pot_ene%total = pot_ene%total + enmr(1) + enmr(2) + enmr(3)
      end if
      if (natc .gt. 0) then
          call get_crd_constraint_energy(natc, entr, atm_jrc, &
                                         crd, temp_frc, atm_xc, atm_weight)
          pot_ene%restraint = pot_ene%restraint + entr
          pot_ene%total = pot_ene%total + entr
      end if
#endif
      
#ifdef MPI
  endif
#endif /* MPI */

  call mpi_bcast(temp_frc, 3*atm_cnt, MPI_DOUBLE_PRECISION, 0, pmemd_comm, err_code_mpi)
  if (new_list) call get_img_frc_distribution(atm_cnt)
  do atm_lst_idx = 1, my_atm_cnt
    i = my_atm_lst(atm_lst_idx)
    frc(1,i) = temp_frc(1,i)
    frc(2,i) = temp_frc(2,i)
    frc(3,i) = temp_frc(3,i)
  end do
! If not calling an External library
endif

! If doing forcefield calculation
if (iskip_ff .ne. 1) then

  if (usemidpoint) then
    !Resetting the forces
    call ensure_alloc_dble3(nb_frc, proc_num_atms_min_bound+proc_ghost_num_atms)
    do j = 1, proc_num_atms+proc_ghost_num_atms
      proc_atm_frc(:, j) = 0.d0
      nb_frc(:,j)=0.0d0
    end do
  else
    do j = 1, gbl_used_img_cnt
      i = gbl_used_img_lst(j)
      img_frc(:, i) = 0.d0
    end do

    if (ti_mode .ne. 0) then
      call ti_zero_arrays
      do j = 1, gbl_used_img_cnt
        i = gbl_used_img_lst(j)
        ti_img_frc(:, :, i) = 0.d0
      end do
      call zero_extra_used_atm_ti_img_frcs(ti_img_frc)
    end if
  end if

  ! We also need to zero anything on the extra used atom list since nonbonded
  ! forces for it will get sent to the atom owner.  We don't calc any such
  ! forces for these atoms, but they are on the send atom list in order to
  ! get their coordinates updated.

  if (.not. usemidpoint) then
    call zero_extra_used_atm_img_frcs(img_frc)
  end if

  ! Don't do recip if PME is not invoked. Don't do it this step unless
  ! mod(irespa,nrepsa) = 0

  onstep = mod(irespa, nrespa) .eq. 0

  params_may_change = (nmropt .ne. 0) .or. (infe .ne. 0)

  ! Calculate the exclusions correction for ips if using it
  ! Clear the force array and the atom-ordered nonbonded force array.
  ! We delay copying the nonbonded forces into the force
  ! array in order to be able to schedule i/o later, and batch it up.

  evdwex = 0.0
  eelex = 0.0

  if (.not. usemidpoint) then
    if ( ips .gt. 0 ) then

      do atm_lst_idx = 1, my_atm_cnt
        i = my_atm_lst(atm_lst_idx)
        frc(:,i) = 0.d0
        nb_frc(:,i) = 0.d0
        ips_frc(:,i) = 0.d0
      end do
      call eexips(evdwex,eelex,ips_frc,crd, &
                  gbl_img_qterm,gbl_ipairs,atm_nb_maskdata, &
                  atm_nb_mask,img_atm_map,my_atm_cnt,gbl_tranvec,my_atm_lst)

      if (new_list) call get_img_frc_ips_distribution(atm_cnt)
      call distribute_img_ips_frcs(atm_cnt, ips_frc)
    else
      do atm_lst_idx = 1, my_atm_cnt
        i = my_atm_lst(atm_lst_idx)
        frc(:,i) = 0.d0
        nb_frc(:,i) = 0.d0
      end do
      if (ti_mode .ne. 0) then
        do atm_lst_idx = 1, my_atm_cnt
          i = my_atm_lst(atm_lst_idx)
          ti_nb_frc(:, :, i) = 0.d0
        end do
      end if
    end if
  end if

  if (use_pme .ne. 0 .and. onstep) then

    ! Self energy:

    ! The small amount of time used here gets lumped with the recip stuff...

    if (usemidpoint) then
      !!!?? Need to verify if the following version is ok or not
      if (master) then
        call self(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
                  params_may_change)
      end if
    else
      if (gbl_frc_ene_task_lst(1) .eq. mytaskid) then
        call self(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
                  params_may_change)
      end if
    end if


    ! Reciprocal energy:

    ! We intentionally keep the load balance counter running through the
    ! fft dist transposes in the recip code; synchronization will mess up the
    ! times a bit, but when not all nodes are doing recip calcs, it will
    ! really help with load balancing.
    !shall we use min max here?

    if (usemidpoint) then
      call update_pme_time(pme_misc_timer)
#if defined(pmemd_SPDP)
      do i = 1 , proc_num_atms + proc_ghost_num_atms ! my_img_lo, my_img_hi
        proc_crd_q_sp(1:3,i) = proc_atm_crd(1:3,i)
        proc_crd_q_sp(4,i) = proc_atm_qterm(i)
      end do
#ifdef _OPENMP_
      call update_pme_time(get_nb_pack_timer)
#endif
#endif /* pmemd_SPDP*/
      call do_kspace_midpoint(nb_frc,pot_ene%elec_recip, vir%elec_recip, need_pot_enes, need_virials)
      !DBG_ARRAYS_DUMP_3DBLE("frc_fft", proc_atm_to_full_list, proc_atm_frc, proc_num_atms+proc_ghost_num_atms)

    else ! not usemidpoint below
      if (i_do_recip) then
        if (ti_mode .eq. 0) then
          call update_pme_time(pme_misc_timer)
          call update_loadbal_timer(elapsed_100usec_atmuser)
          if (block_fft .eq. 0) then
            call do_slab_pmesh_kspace(img_frc, pot_ene%elec_recip, &
                                      vir%elec_recip, need_pot_enes, need_virials)
          else
            call do_blk_pmesh_kspace(img_frc, pot_ene%elec_recip, &
                                     vir%elec_recip, need_pot_enes, need_virials)
          end if
          call update_loadbal_timer(elapsed_100usec_recipfrc)
          if (nrespa .gt. 1) call respa_scale(atm_cnt, img_frc, nrespa)
        else
          call update_pme_time(pme_misc_timer)
          call update_loadbal_timer(elapsed_100usec_atmuser)
          ti_mask_piece = 2
          if (block_fft .eq. 0) then
            call do_slab_pmesh_kspace(img_frc, pot_ene%elec_recip, ti_vir(1,:,:), need_pot_enes, need_virials)

          else
            call do_blk_pmesh_kspace(img_frc, pot_ene%elec_recip, &
                                     ti_vir(1,:,:), need_pot_enes, need_virials)
          end if
          ! It is much more efficient to do it this way, rather
          ! than passing ti_img_frc(1,:,:) above ...
          do j = 1, gbl_used_img_cnt
            i = gbl_used_img_lst(j)
            ti_img_frc(1, :, i) = img_frc(:, i)
            img_frc(:,i) = 0.d0
          end do

          ti_pot(1) = pot_ene%elec_recip

          pot_ene%elec_recip = 0.d0
          ti_mask_piece = 1
          if (ti_mode .ne. 3) then
            if (block_fft .eq. 0) then
              call do_slab_pmesh_kspace(img_frc, pot_ene%elec_recip, ti_vir(2,:,:), need_pot_enes, need_virials)

            else
              call do_blk_pmesh_kspace(img_frc, pot_ene%elec_recip, ti_vir(2,:,:), need_pot_enes, need_virials)
            end if
            do j = 1, gbl_used_img_cnt
              i = gbl_used_img_lst(j)
              ti_img_frc(2, :, i) = img_frc(:, i)
              img_frc(:,i) = 0.d0
            end do
            ti_pot(2) = pot_ene%elec_recip
          else
            ti_pot(2) = ti_pot(1)
            do j = 1, gbl_used_img_cnt
              i = gbl_used_img_lst(j)
              ti_img_frc(2, :, i) = ti_img_frc(1, :, i)
            end do
            ti_vir(2,:,:) = ti_vir(1,:,:)
          end if
          call update_loadbal_timer(elapsed_100usec_recipfrc)

          do j = 1, gbl_used_img_cnt
            i = gbl_used_img_lst(j)
            ti_img_frc(1, :, i) = ti_img_frc(1, :, i) * ti_weights(1)
            ti_img_frc(2, :, i) = ti_img_frc(2, :, i) * ti_weights(2)
          end do

          vir%elec_recip = ti_vir(1,:,:) * ti_weights(1) + &
            ti_vir(2,:,:) * ti_weights(2)

          ti_vve(ti_vir0) = (ti_vir(1,1,1) + &
                             ti_vir(1,2,2) + &
                             ti_vir(1,3,3)) * ti_weights(1)
          ti_vve(ti_vir1) = (ti_vir(2,1,1) + &
                             ti_vir(2,2,2) + &
                             ti_vir(2,3,3)) * ti_weights(2)

          call ti_update_ene_all(ti_pot, si_elect_ene, pot_ene%elec_recip,3)

          if (nrespa .gt. 1) then
            call ti_respa_scale(atm_cnt, ti_img_frc, nrespa)
          end if
        end if
      end if
    end if ! usemidpoint


    ! Long range dispersion contributions:

    ! Continuum method:

    if (vdwmeth .eq. 1) then
      if ((.not. usemidpoint .and. gbl_frc_ene_task_lst(1) .eq. mytaskid) .or. &
          (usemidpoint .and. master)) then
        call vdw_correction(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                            params_may_change)
      end if
    end if
  end if      ! respa .and. use_pme

  if (.not. usemidpoint) then
    call update_loadbal_timer(elapsed_100usec_atmuser)
  end if ! not usemidpoint

  ! Direct part of ewald plus vdw, hbond, force and energy calculations:

  call update_pme_time(pme_misc_timer)
  if (ips .gt. 0) then
    call get_nb_ips_energy(img_frc, gbl_img_crd, gbl_img_qterm, gbl_eed_cub, &
                           gbl_ipairs, gbl_tranvec, need_pot_enes, need_virials, &
                           pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                           vir%eedvir, vir%elec_direct)
  else if (ti_mode .ne. 0) then
    ! This is the only part of the TI code that uses img_frc, everything else
    ! writes to ti_img_frc
    call get_nb_energy_ti(img_frc, gbl_img_crd, gbl_img_qterm, gbl_eed_cub, &
                          gbl_ipairs, gbl_tranvec, need_pot_enes, need_virials, &
                          pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                          vir%eedvir, vir%elec_direct, gbl_img_atm_map)
  else if ((lj1264 .eq. 1) .and. (plj1264 .eq. 0)) then ! New2022
    call get_nb_energy_1264(img_frc, gbl_img_crd, gbl_img_qterm, gbl_eed_cub, &
                            gbl_ipairs, gbl_tranvec, need_pot_enes, need_virials, &
                            pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                            vir%eedvir, vir%elec_direct)
  else if ((lj1264 .eq. 0) .and. (plj1264 .eq. 1)) then
    call get_nb_energy_p1264(img_frc, gbl_img_crd, gbl_img_qterm, gbl_eed_cub, &
                            gbl_ipairs, gbl_tranvec, need_pot_enes, need_virials, &
                            pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                            vir%eedvir, vir%elec_direct)
  else if ((lj1264 .eq. 1) .and. (plj1264 .eq. 1)) then
    call get_nb_energy_1264p1264(img_frc, gbl_img_crd, gbl_img_qterm, gbl_eed_cub, &
                            gbl_ipairs, gbl_tranvec, need_pot_enes, need_virials, &
                            pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                            vir%eedvir, vir%elec_direct) 
  else

    if (usemidpoint) then
#ifdef _OPENMP_
      !print *,"Before: force_thrd:",frc_thread(:,1)
#ifdef pmemd_SPDP
      call get_nb_energy_midpoint(frc_thread, proc_crd_q_sp, gbl_eed_cub, &
                                  proc_ipairs, need_pot_enes, need_virials, &
                                  pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                                  vir%eedvir, vir%elec_direct)
#else
      call get_nb_energy_midpoint(frc_thread, proc_atm_crd, proc_atm_qterm, gbl_eed_cub, &
                                  proc_ipairs, need_pot_enes, need_virials, &
                                  pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                                  vir%eedvir, vir%elec_direct)
#endif /* pmemd_SPDP */
#else
#ifdef pmemd_SPDP
      call get_nb_energy_midpoint(nb_frc, proc_crd_q_sp, gbl_eed_cub, &
                                  proc_ipairs, need_pot_enes, need_virials, &
                                  pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                                  vir%eedvir, vir%elec_direct)
#else
      call get_nb_energy_midpoint(nb_frc, proc_atm_crd, proc_atm_qterm, gbl_eed_cub, &
                                  proc_ipairs, need_pot_enes, need_virials, &
                                  pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                                  vir%eedvir, vir%elec_direct)
#endif
#endif /* _OPENMP_ */
      !DBG_ARRAYS_DUMP_3DBLE("frc_nbe", proc_atm_to_full_list, proc_atm_frc, proc_num_atms+proc_ghost_num_atms)
    else
      call get_nb_energy(img_frc, gbl_img_crd, gbl_img_qterm, gbl_eed_cub, &
                         gbl_ipairs, gbl_tranvec, need_pot_enes, need_virials, &
                         pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                         vir%eedvir, vir%elec_direct)
    end if
  end if

  if (usemidpoint) then
#ifdef _OPENMP_
    do j = 0, nthreads-1
      offset_loc = j*(proc_num_atms + proc_ghost_num_atms)
      do k = 1,  proc_num_atms + proc_ghost_num_atms
        proc_atm_frc(1:3,k) = proc_atm_frc(1:3,k) + frc_thread(1:3,k+offset_loc)
        frc_thread(1:3,k+offset_loc) = 0.d0
      end do  ! k = 1,  proc_num_atms + proc_ghost_num_atms
    end do  ! do j=0,nthreads-1
#endif /* _OPENMP_ */
  end if

  call update_pme_time(dir_frc_sum_timer)

  if (.not. usemidpoint) then
    if (efx .ne. 0 .or. efy .ne. 0 .or. efz .ne. 0) then
      call get_efield_energy(img_frc, crd, gbl_img_qterm, gbl_img_atm_map, need_pot_enes, &
                             pot_ene%efield, atm_cnt, nstep)
    end if
  end if

  call update_pme_time(dir_frc_sum_timer)
  if (.not. usemidpoint) then
    call update_loadbal_timer(elapsed_100usec_dirfrc)
  end if

  call update_pme_time(pme_misc_timer)

  ! Calculate 1-4 electrostatic energies, forces:
  if (usemidpoint) then
    if (need_virials) then
      call get_nb14_energy_midpoint(proc_atm_qterm, proc_atm_crd, nb_frc, proc_iac, typ_ico, &
                                    gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt, &
                                    pot_ene%elec_14, pot_ene%vdw_14, &
                                    vir%elec_14)
    else
      call get_nb14_energy_midpoint(proc_atm_qterm, proc_atm_crd, nb_frc, proc_iac, typ_ico, &
                                    gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt, &
                                    pot_ene%elec_14, pot_ene%vdw_14)
    end if
  else
    if (charmm_active) then
      if (need_virials) then
        call get_nb14_energy(atm_qterm, crd, nb_frc, atm_iac, typ_ico, &
                             gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt, &
                             pot_ene%elec_14, pot_ene%vdw_14, &
                             vir%elec_14)
      else
        call get_nb14_energy(atm_qterm, crd, nb_frc, atm_iac, typ_ico, &
                             gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt, &
                             pot_ene%elec_14, pot_ene%vdw_14)
      end if
    else
      if (need_virials) then
        call get_nb14_energy(atm_qterm, crd, nb_frc, atm_iac, typ_ico, &
                             gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt, &
                             pot_ene%elec_14, pot_ene%vdw_14, &
                             vir%elec_14)
      else
        call get_nb14_energy(atm_qterm, crd, nb_frc, atm_iac, typ_ico, &
                             gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt, &
                             pot_ene%elec_14, pot_ene%vdw_14)
      end if
    end if
  end if

  call update_pme_time(dir_frc_sum_timer)

  ! Adjust energies, forces for masked out pairs:

  if (use_pme .ne. 0) then
    if (usemidpoint) then
      call nb_adjust_midpoint(proc_num_atms+proc_ghost_num_atms, proc_atm_qterm, &
                              proc_atm_crd, gbl_nb_adjust_pairlst, gbl_eed_cub, &
                              nb_frc, pot_ene%elec_nb_adjust, vir%elec_nb_adjust)
    else
      call nb_adjust(my_atm_cnt, my_atm_lst, &
                     atm_qterm, crd, gbl_nb_adjust_pairlst, gbl_eed_cub, &
                     nb_frc, pot_ene%elec_nb_adjust, vir%elec_nb_adjust)
    end if
  end if

  call update_pme_time(adjust_masked_timer)

  ! Get COM-relative coordinates.
  ! Calculate total nonbonded energy components.

  pot_ene%vdw_tot = pot_ene%vdw_dir + &
                    pot_ene%vdw_recip

  pot_ene%elec_tot = pot_ene%elec_dir + &
                     pot_ene%elec_recip + &
                     pot_ene%elec_nb_adjust + &
                     pot_ene%elec_self

  if (ips .gt. 0) then
    pot_ene%vdw_tot = pot_ene%vdw_tot + evdwex
    pot_ene%elec_tot = pot_ene%elec_tot + eelex
  end if

  if (need_virials) then
    if (usemidpoint) then
      call ensure_alloc_dble3(atm_rel_crd, proc_num_atms_min_bound+proc_ghost_num_atms)
      call get_atm_rel_crd_midpoint(gbl_mol_cnt, gbl_mol_com, proc_atm_crd, atm_rel_crd, pbc_box)
    else
      call get_atm_rel_crd(my_mol_cnt, gbl_mol_com, crd, atm_rel_crd, gbl_my_mol_lst)
    end if
  end if

  call update_pme_time(pme_misc_timer)
  call update_time(nonbond_time)

  if (usemidpoint) then
    if (igamd .eq. 4 .or. igamd .eq. 5) then
      nb_frc_gamd = frc
    end if
  end if

  ! Calculate the other contributions:

  ! Global to Local atomic index
  ! In order to get local atomic index from global(i.e.as it is written in inpcrd)
  ! index
  !!!NOT sure if we need here since we already have it in build_atm_space
  !routine
  if (usemidpoint) then
    call pme_bonded_force_midpoint(pot_ene, new_list)
    ! DBG_ARRAYS_DUMP_3DBLE("frc_bond", proc_atm_to_full_list, proc_atm_frc, proc_num_atms+proc_ghost_num_atms)

    do i=1,proc_num_atms + proc_ghost_num_atms
      proc_atm_frc(:,i) = proc_atm_frc(:,i) + nb_frc(:,i)
    end do
  else
    call pme_bonded_force(crd, frc, pot_ene)
  end if

  ! The above stuff gets lumped as "other time"...

  ! --- calculate the EMAP constraint energy ---

  if (.not. usemidpoint) then
    if (iemap .gt. 0) then
      call emapforce(natom,enemap,atm_mass,crd,frc )
      pot_ene%emap = enemap
      pot_ene%restraint = pot_ene%restraint + enemap
    end if
  end if

  ! Sum up total potential energy for this task:

  pot_ene%total = pot_ene%vdw_tot + &
                  pot_ene%elec_tot + &
                  pot_ene%hbond + &
                  pot_ene%bond + &
                  pot_ene%angle + &
                  pot_ene%dihedral + &
                  pot_ene%vdw_14 + &
                  pot_ene%elec_14 + &
                  pot_ene%restraint + &
                  pot_ene%imp + &
                  pot_ene%angle_ub + &
                  pot_ene%efield + &
                  pot_ene%cmap

  if (.not. usemidpoint) then
    call update_loadbal_timer(elapsed_100usec_atmowner)

    !add the ips non bonded contribution
    if ( ips .gt. 0 ) then
      do atm_lst_idx = 1, my_atm_cnt
        i = my_atm_lst(atm_lst_idx)
        nb_frc(:,i) = nb_frc(:,i) + ips_frc(:,i)
      end do
    end if
  end if

  if (usemidpoint) then
    ! Return forces to local mpi for collection

    call comm_frc(proc_atm_frc)
    call comm_frc(nb_frc)

  else
    if (new_list) call get_img_frc_distribution(atm_cnt)

    if (ti_mode .eq. 0) then
      call distribute_img_frcs(atm_cnt, img_frc, nb_frc)
    else
      call distribute_img_ti_frcs(atm_cnt, ti_img_frc, ti_nb_frc)
    end if
  end if

  call update_time(fcve_dist_time)

  call zero_pme_time()

  ! If using extra points and a frame (checked internal to subroutine),
  ! transfer force and torque from the extra points to the parent atom:

  if (.not. usemidpoint) then
    if (numextra .gt. 0 .and. frameon .ne. 0) then
      if (ti_mode .eq. 0) then
        call orient_frc(crd, nb_frc, vir%ep_frame, ep_frames, ep_lcl_crd, gbl_frame_cnt)
      else
        ti_vir(:,:,:) = 0.d0
        ! Use nb_frc as a temporary array, so we don't need to splice
        ! the entire force array (ti_nb_frc)
        do atm_lst_idx = 1, my_atm_cnt
          i = my_atm_lst(atm_lst_idx)
          nb_frc(:,i) = ti_nb_frc(1, :, i)
        end do
        call orient_frc(crd, nb_frc, ti_vir(1,:,:), ep_frames, ep_lcl_crd,  gbl_frame_cnt)
        do atm_lst_idx = 1, my_atm_cnt
          i = my_atm_lst(atm_lst_idx)
          ti_nb_frc(1, :, i) = nb_frc(:,i)
          nb_frc(:,i) = ti_nb_frc(2, :, i)
        end do
        call orient_frc(crd, nb_frc, ti_vir(2,:,:), ep_frames, ep_lcl_crd, gbl_frame_cnt)
        do atm_lst_idx = 1, my_atm_cnt
          i = my_atm_lst(atm_lst_idx)
          ti_nb_frc(2, :, i) = nb_frc(:,i)
          nb_frc(:,i) = 0.d0
        end do
        vir%ep_frame = ti_vir(1,:,:) + ti_vir(2,:,:)
      end if
    end if
  end if
  ! If the net force correction is in use, here we determine the net forces by
  ! looking at the sum of all nonbonded forces.  This should give the same
  ! result as just looking at the reciprocal forces, but it is more
  ! computationally convenient, especially for extra points, to do it this way.

  if (netfrc .gt. 0 .and. onstep) then
    if (ti_mode .eq. 0) then
      if (usemidpoint) then
        do i = 1, proc_num_atms
          net_frcs(:) = net_frcs(:) + nb_frc(:, i)
        end do
      else
        do atm_lst_idx = 1, my_atm_cnt
          i = my_atm_lst(atm_lst_idx)
          net_frcs(:) = net_frcs(:) + nb_frc(:, i)
        end do
      end if
    else
      do atm_lst_idx = 1, my_atm_cnt
        i = my_atm_lst(atm_lst_idx)
        ti_net_frcs(1, :) = ti_net_frcs(1, :) + ti_nb_frc(1, :, i)
        ti_net_frcs(2, :) = ti_net_frcs(2, :) + ti_nb_frc(2, :, i)
      end do
    end if
  end if

  ! First phase of virial work.  We need just the nonbonded forces at this
  ! stage.

  if (need_virials) then

    ! The relative coordinates for the extra points will mess up the molvir
    ! netfrc correction here unless we zero them out.  They are of no import
    ! in ekcom calc because extra points have a mass of 0.

    if (usemidpoint) then
      do atm_lst_idx = 1, proc_num_atms
        vir%molecular(:,1) = vir%molecular(:,1) + nb_frc(:,atm_lst_idx) * atm_rel_crd(1,atm_lst_idx)
        vir%molecular(:,2) = vir%molecular(:,2) + nb_frc(:,atm_lst_idx) * atm_rel_crd(2,atm_lst_idx)
        vir%molecular(:,3) = vir%molecular(:,3) + nb_frc(:,atm_lst_idx) * atm_rel_crd(3,atm_lst_idx)
        molvir_netfrc_corr(:, 1) = molvir_netfrc_corr(:, 1) + atm_rel_crd(1,atm_lst_idx)
        molvir_netfrc_corr(:, 2) = molvir_netfrc_corr(:, 2) + atm_rel_crd(2,atm_lst_idx)
        molvir_netfrc_corr(:, 3) = molvir_netfrc_corr(:, 3) + atm_rel_crd(3,atm_lst_idx)
      end do
    else ! not usemidpoint below
      if (numextra .gt. 0 .and. frameon .ne. 0) &
        call zero_extra_pnts_vec(atm_rel_crd, ep_frames, gbl_frame_cnt)

      if (ti_mode .eq. 0) then
        do atm_lst_idx = 1, my_atm_cnt
          i = my_atm_lst(atm_lst_idx)
          vir%molecular(:,1) = vir%molecular(:,1) + nb_frc(:,i) * atm_rel_crd(1,i)
          vir%molecular(:,2) = vir%molecular(:,2) + nb_frc(:,i) * atm_rel_crd(2,i)
          vir%molecular(:,3) = vir%molecular(:,3) + nb_frc(:,i) * atm_rel_crd(3,i)
          molvir_netfrc_corr(:, 1) = molvir_netfrc_corr(:, 1) + atm_rel_crd(1, i)
          molvir_netfrc_corr(:, 2) = molvir_netfrc_corr(:, 2) + atm_rel_crd(2, i)
          molvir_netfrc_corr(:, 3) = molvir_netfrc_corr(:, 3) + atm_rel_crd(3, i)
        end do
      else
        do atm_lst_idx = 1, my_atm_cnt
          i = my_atm_lst(atm_lst_idx)

          vir%molecular(:,1) = vir%molecular(:,1) + &
            ti_nb_frc(1,:,i) * atm_rel_crd_sc(1,1,i) + &
            ti_nb_frc(2,:,i) * atm_rel_crd_sc(2,1,i)
          vir%molecular(:,2) = vir%molecular(:,2) + &
            ti_nb_frc(1,:,i) * atm_rel_crd_sc(1,2,i) + &
            ti_nb_frc(2,:,i) * atm_rel_crd_sc(2,2,i)
          vir%molecular(:,3) = vir%molecular(:,3) + &
            ti_nb_frc(1,:,i) * atm_rel_crd_sc(1,3,i) + &
            ti_nb_frc(2,:,i) * atm_rel_crd_sc(2,3,i)

          ti_molvir_netfrc_corr(1,:,1) = ti_molvir_netfrc_corr(1,:,1) + &
            atm_rel_crd_sc(1,1,i)
          ti_molvir_netfrc_corr(1,:,2) = ti_molvir_netfrc_corr(1,:,2) + &
            atm_rel_crd_sc(1,2,i)
          ti_molvir_netfrc_corr(1,:,3) = ti_molvir_netfrc_corr(1,:,3) + &
            atm_rel_crd_sc(1,3,i)

          ti_molvir_netfrc_corr(2,:,1) = ti_molvir_netfrc_corr(2,:,1) + &
            atm_rel_crd_sc(2,1,i)
          ti_molvir_netfrc_corr(2,:,2) = ti_molvir_netfrc_corr(2,:,2) + &
            atm_rel_crd_sc(2,2,i)
          ti_molvir_netfrc_corr(2,:,3) = ti_molvir_netfrc_corr(2,:,3) + &
            atm_rel_crd_sc(2,3,i)
        end do
      end if
    end if ! usemidpoint
    ! Finish up virial work; Timing is inconsequential...

    vir%atomic(:,:) = vir%elec_recip(:,:) + &
      vir%elec_direct(:,:) + &
      vir%elec_nb_adjust(:,:) + &
      vir%elec_recip_vdw_corr(:,:) + &
      vir%elec_recip_self(:,:) + &
      vir%elec_14(:,:) + &
      vir%ep_frame(:,:)

    !add the ips contribution to virials
    if ( ips .gt. 0 ) then
      vir%atomic(:,:) = vir%atomic(:,:) + VIREXIPS(:,:)
    end if

    vir%molecular(:,:) = vir%molecular(:,:) + vir%atomic(:,:)

    if (usemidpoint) then
      call get_ekcom_midpoint(gbl_mol_cnt, gbl_mol_mass_inv, ekcmt, proc_atm_vel, &
        proc_atm_mass)
    else
      call get_ekcom(my_mol_cnt, gbl_mol_mass_inv, ekcmt, atm_vel, &
        atm_mass, gbl_my_mol_lst)
    end if
  end if

  call update_time(nonbond_time)
  call update_pme_time(pme_misc_timer)

  ! Add a bunch of values (energies, virials, etc.) together as needed.

  call dist_enes_virs_netfrcs(pot_ene, vir, ekcmt, net_frcs, &
    need_pot_enes, need_virials)

  if (ti_mode .ne. 0 ) then
    call ti_dist_enes_virs_netfrcs(onstep, need_virials)
    call ti_combine_frcs(atm_cnt, my_atm_cnt, my_atm_lst, nb_frc, ti_nb_frc)
  end if

  call update_time(fcve_dist_time)
  call zero_pme_time()

  if (.not. usemidpoint) then
    if (iamd.gt.1) then
      call calculate_amd_dih_weights(atm_cnt,pot_ene%dihedral,frc,crd)
    end if

  if((igamd.eq.2).or.(igamd.eq.3).or.(igamd.eq.5).or.(igamd.eq.7).or.(igamd.eq.9))then
      call calculate_gamd_dih_weights(atm_cnt,pot_ene%dihedral,frc,crd)
  endif
endif ! not usemidpoint

  ! Copy image forces to atom forces, correcting for netfrc as you go, if
  ! appropriate.  Do not remove net force if netfrc = 0; e.g. in minimization.

  if (netfrc .gt. 0 .and. onstep) then

    if (ti_mode .eq. 0) then
      if (usemidpoint) then
        net_frcs(:) = net_frcs(:) / dble(atm_cnt - numextra)

        do i = 1, proc_num_atms
          proc_atm_frc(:, i) = proc_atm_frc(:, i) - net_frcs(:)
        end do
      else
        net_frcs(:) = net_frcs(:) / dble(atm_cnt - numextra)

        do atm_lst_idx = 1, my_atm_cnt
          i = my_atm_lst(atm_lst_idx)
          frc(:, i) = frc(:, i) + nb_frc(:, i) - net_frcs(:)
        end do
      end if
    else
      ti_net_frcs(1,:) = ti_net_frcs(1,:)/dble(ti_atm_cnt(1)-ti_numextra_pts(1))
      ti_net_frcs(2,:) = ti_net_frcs(2,:)/dble(ti_atm_cnt(2)-ti_numextra_pts(2))
      net_frcs(:) = ti_net_frcs(1,:) + ti_net_frcs(2,:)

      do atm_lst_idx = 1, my_atm_cnt
        i = my_atm_lst(atm_lst_idx)
        ! This matches how sander removes netfrcs in TI runs
        if (ti_lst(1,i) .ne. 0) then
          frc(:, i) = frc(:, i) + nb_frc(:, i) - ti_net_frcs(1,:)
        else if (ti_lst(2,i) .ne. 0) then
          frc(:, i) = frc(:, i) + nb_frc(:, i) - ti_net_frcs(2,:)
        else
          frc(:, i) = frc(:, i) + nb_frc(:, i) - net_frcs(:)
        end if
      end do
    end if

    if (numextra .gt. 0 .and. frameon .ne. 0) &
      call zero_extra_pnts_vec(frc, ep_frames, gbl_frame_cnt)

    ! Correct the molecular virial for netfrc:

    if (need_virials) then
      if (ti_mode .eq. 0) then
        vir%molecular(1,:) = vir%molecular(1,:) - &
          molvir_netfrc_corr(1,:) * net_frcs(1)
        vir%molecular(2,:) = vir%molecular(2,:) - &
          molvir_netfrc_corr(2,:) * net_frcs(2)
        vir%molecular(3,:) = vir%molecular(3,:) - &
          molvir_netfrc_corr(3,:) * net_frcs(3)
      else
        vir%molecular(1,:) = vir%molecular(1,:) - &
          ti_molvir_netfrc_corr(1,1,:) * ti_net_frcs(1,1) - &
          ti_molvir_netfrc_corr(2,1,:) * ti_net_frcs(2,1)
        vir%molecular(2,:) = vir%molecular(2,:) - &
          ti_molvir_netfrc_corr(1,2,:) * ti_net_frcs(1,2) - &
          ti_molvir_netfrc_corr(2,2,:) * ti_net_frcs(2,2)
        vir%molecular(3,:) = vir%molecular(3,:) - &
          ti_molvir_netfrc_corr(1,3,:) * ti_net_frcs(1,3) - &
          ti_molvir_netfrc_corr(2,3,:) * ti_net_frcs(2,3)
      end if
    end if

  else

    if (.not. usemidpoint) then
      do atm_lst_idx = 1, my_atm_cnt
        i = my_atm_lst(atm_lst_idx)
        frc(1, i) = frc(1, i) + nb_frc(1, i)
        frc(2, i) = frc(2, i) + nb_frc(2, i)
        frc(3, i) = frc(3, i) + nb_frc(3, i)
      end do
    end if

  end if

  ! Calculate vir_vs_ene in master.

  if (master) then

    if (ti_mode .eq. 0) then
      vir_vs_ene = vir%elec_recip(1, 1) + &
        vir%elec_recip(2, 2) + &
        vir%elec_recip(3, 3) + &
        vir%eedvir + &
        vir%elec_nb_adjust(1, 1) + &
        vir%elec_nb_adjust(2, 2) + &
        vir%elec_nb_adjust(3, 3)

      ! Avoid divide-by-zero for pure neutral systems (l-j spheres), or if
      ! energies are not being calculated this step...

      if (pot_ene%elec_tot .ne. 0.0d0) then
        vir_vs_ene = abs(vir_vs_ene + pot_ene%elec_tot)/abs(pot_ene%elec_tot)
      else
        vir_vs_ene = 0.0d0
      end if
    else
      ! Avoid divide-by-zero for pure neutral systems (l-j spheres), or if
      ! energies are not being calculated this step...
      vir_vs_ene = 0.d0
      if (ti_vve(ti_ene0) .ne. 0.0d0) then
        vir_vs_ene = abs(ti_vve(ti_vir0) + &
          ti_vve(ti_ene0) * ti_weights(1))/abs(ti_vve(ti_ene0))
      end if
      if (ti_vve(ti_ene1) .ne. 0.d0) then
        vir_vs_ene = vir_vs_ene + abs(ti_vve(ti_vir1) + &
          ti_vve(ti_ene1) * ti_weights(2))/abs(ti_vve(ti_ene1))
      end if
    end if
    pme_err_est = vir_vs_ene
  end if

  ! Save virials in form used in runmd:

  if (need_virials) then
    virial(1) = 0.5d0 * vir%molecular(1, 1)
    virial(2) = 0.5d0 * vir%molecular(2, 2)
    virial(3) = 0.5d0 * vir%molecular(3, 3)
  end if

! End If doing forcefield calculation (iskip_ff)
endif

    ! Built-in X-ray target function and gradient
#ifndef MPI
#ifndef NOXRAY
    if (xray_active) then
       call xray_get_derivative(crd,frc,nstep,xray_e)
       pot_ene%restraint = pot_ene%restraint + xray_e
       pot_ene%total = pot_ene%total + xray_e
    end if
#endif
#endif

! Calculate the NMR restraint energy contributions, if requested.
if (.not. usemidpoint) then
  if (nmropt .ne. 0) then
    call nmr_calc(crd, frc, enmr, 6)
    pot_ene%restraint = pot_ene%restraint + enmr(1) + enmr(2) + enmr(3)
    pot_ene%total     = pot_ene%total     + enmr(1) + enmr(2) + enmr(3)
  end if
end if

if (ti_mode .ne. 0 .and. need_pot_enes) then
  call ti_calc_dvdl
end if

if (master .and. verbose .gt. 0) then
  call write_netfrc(net_frcs)
  call pme_verbose_print(pot_ene, vir, vir_vs_ene)
end if
  
!restart skipping of forcefield
if (iskip_ff .ne. 1) then

!AMD DUAL BOOST CALC START
   if(iamd.gt.0)then
!calculate totboost and apply weight to frc. frc=frc*fwgt
if(.not. usemidpoint) then
      call calculate_amd_total_weights(atm_cnt, pot_ene%total, pot_ene%dihedral,&
        pot_ene%amd_boost, frc, crd, my_atm_lst)
endif ! not usemidpoint
    ! Update total energy
    pot_ene%total = pot_ene%total + pot_ene%amd_boost
  end if

if(.not. usemidpoint) then
!GaMD DUAL BOOST CALC START
   if(igamd.gt.0)then
!calculate totboost and apply weight to frc. frc=frc*fwgt
     if(igamd.eq.1 .or. igamd.eq.2 .or. igamd.eq.3) then
       call calculate_gamd_total_weights(atm_cnt, &
        pot_ene%total,pot_ene%dihedral,pot_ene%gamd_boost,frc,crd,my_atm_lst)
       ! Update total energy
       pot_ene%total = pot_ene%total + pot_ene%gamd_boost
     else if(igamd.eq.4 .or. igamd.eq.5 ) then
        pot_ene_nb = pot_ene%elec_tot + pot_ene%vdw_tot + pot_ene%vdw_14 + pot_ene%elec_14 + pot_ene%efield
        nb_frc_tmp = nb_frc_gamd
        call calculate_gamd_nb_weights(atm_cnt,pot_ene_nb,pot_ene%dihedral, &
          pot_ene%gamd_boost,nb_frc_gamd,crd,my_atm_lst)
       ! Update total force
        frc = frc - nb_frc_tmp + nb_frc_gamd
       ! Update total energy
       pot_ene%total = pot_ene%total + pot_ene%gamd_boost
     else if( (igamd.eq.6 .or. igamd.eq.7) .and. (ti_mode.ne.0)) then
       pot_ene_lig = ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                     ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                     ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
       pot_ene_dih = pot_ene%total + &
                      ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene)
       call calculate_gamd_ti_region_weights(atm_cnt,pot_ene_lig, pot_ene_dih, &
                       pot_ene%gamd_boost,fwgt,frc,crd,ti_sc_lst,my_atm_lst)
       ! Calculate net force on TI region
       if (mod(nstep+1,ntave).eq.0) then
         net_frc_ti    = 0.0
         do i = 1, atm_cnt
               if ( ti_sc_lst(i) .gt. 0 ) then
                       net_frc_ti(:) = net_frc_ti(:) + frc(:,i)
               end if
         end do
         write(*,'(a,i10,4f15.5)') "| GaMD: step, frc_weight, frc_lig(1:3) = ", nstep+1, fwgt, net_frc_ti
       end if
       ! Update total energy
       ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
     else if( (igamd.eq.8 .or. igamd.eq.9) .and. (ti_mode.ne.0)) then
       pot_ene_lig = ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene)
       pot_ene_dih = pot_ene%total + &
                     ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene) + &
                     ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                     ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
       call calculate_gamd_ti_region_weights(atm_cnt,pot_ene_lig, pot_ene_dih, &
                       pot_ene%gamd_boost,fwgt,frc,crd,ti_sc_lst,my_atm_lst)
       ! Update total energy
       ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
     else if( (igamd.eq.10 .or. igamd.eq.11) .and. (ti_mode.ne.0)) then
        pot_ene_lig = ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                      ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
        pot_ene_dih = pot_ene%total + &
                      ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                      ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene)
       call calculate_gamd_ti_region_weights(atm_cnt,pot_ene_lig, pot_ene_dih, &
                       pot_ene%gamd_boost,fwgt,frc,crd,ti_sc_lst,my_atm_lst)
       ! Update total energy
       ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
      end if
    end if

!scaledMD scale forces
   if(scaledMD.gt.0)then
!apply scale to frc. frc=frc*scaledMD_lambda
      call scaledMD_scale_frc(atm_cnt, pot_ene%total, &
                              frc, crd, my_atm_lst)
      pot_ene%total = pot_ene%total * scaledMD_lambda
    end if

 endif ! not usemidpoint

! Endff of iskip_ff (skipping forcefield)
endif

  ! If belly is on then set the belly atom forces to zero:

  if (ibelly .gt. 0) call bellyf(atm_cnt, atm_igroup, frc)

  call update_time(nonbond_time)
  call update_pme_time(pme_misc_timer)

  if (.not. usemidpoint) then
    if (charmm_active .and. do_charmm_dump_gold .eq. 1) then
      call mpi_gathervec(atm_cnt, frc)
      if (master) then
        call pme_charmm_dump_gold(atm_cnt,frc,pot_ene)
        write(mdout, '(a)') 'charmm_gold() completed. Exiting'
      end if
      call mexit(6, 0)
    end if

    if (igamd .eq. 4 .or. igamd .eq. 5) then
      deallocate(nb_frc_gamd,nb_frc_tmp)
    end if

    if (infe.gt.0 .or. iextpot.gt.0) then
        if (infe.gt.0) then
           bias_frc(:,:) = 0.d0
           call nfe_on_force(crd,bias_frc, pot_ene%nfe)
           do atm_lst_idx = 1, my_atm_cnt
              i = my_atm_lst(atm_lst_idx)
              frc(1,i)=frc(1,i)+bias_frc(1,i)
              frc(2,i)=frc(2,i)+bias_frc(2,i)
              frc(3,i)=frc(3,i)+bias_frc(3,i)
           end do
        endif

        if (iextpot.gt.0) then
           temp_holding = 0.0
           call pme_external(crd, frc, temp_holding)
        endif

        if (infe.gt.0) then
           pot_ene%total     = pot_ene%total + pot_ene%nfe
           pot_ene%restraint = pot_ene%restraint + pot_ene%nfe
        endif

        if (iextpot.gt.0) then
           pot_ene%total     = pot_ene%total + temp_holding
           pot_ene%restraint = pot_ene%restraint + temp_holding
        endif
    end if




    if (ineb .gt. 0) then
      call mpi_gathervec(atm_cnt, frc)
      call mpi_gathervec(atm_cnt, crd)

      if (mytaskid .eq. 0) then
        call full_neb_forces(atm_mass, crd, frc, pot_ene%total, fitmask, rmsmask, nstep)
      end if
      call mpi_bcast(neb_force, 3*nattgtrms, MPI_DOUBLE_PRECISION, 0, pmemd_comm, err_code_mpi)

      do jatm=1,nattgtrms
        iatm=rmsmask(jatm)
        if (iatm .eq. 0) cycle
        if (gbl_atm_owner_map(iatm) .eq. mytaskid) then
          frc(1,iatm)= frc(1,iatm)+neb_force(1,jatm)
          frc(2,iatm)= frc(2,iatm)+neb_force(2,jatm)
          frc(3,iatm)= frc(3,iatm)+neb_force(3,jatm)
        end if
      enddo

      if (beadid.eq.1 .or. beadid.eq.neb_nbead) then
        do i=1, my_atm_cnt
          iatm = my_atm_lst(i)
          frc(:,iatm)=0.d0
        enddo
      end if

    end if

  end if ! not usemidpoint

#endif /* CUDA */

  return

end subroutine pme_force

subroutine pme_force_midpoint(atm_cnt, crd, frc, &
                              my_atm_lst, new_list, need_pot_enes, need_virials, &
                              pot_ene, nstep, virial, ekcmt, pme_err_est)

  use bonds_midpoint_mod
  use mol_list_mod, only : gbl_mol_cnt
  use angles_mod
  use angles_ub_mod
  use bonds_mod
  use constraints_mod
  use dihedrals_mod
  use dihedrals_imp_mod
  use cmap_mod
#ifndef MPI
#ifndef NOXRAY
  use xray_interface_module, only: xray_get_derivative
#endif
#endif
  use dynamics_dat_mod
  use dynamics_mod
  use ene_frc_splines_mod
  use emap_mod,only:emapforce
  use nb_exclusions_mod
  use pme_direct_mod
  use pme_recip_dat_mod
  use pme_blk_recip_mod
  use pme_slab_recip_mod
  use loadbal_mod
  use img_mod
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use nb_pairlist_mod
  use nmr_calls_mod
  use parallel_dat_mod
  use parallel_mod
  use pbc_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
  use runfiles_mod
  use timers_mod
  use charmm_mod, only : charmm_active
  use mdin_debugf_dat_mod, only : do_charmm_dump_gold
  use nbips_mod
  use amd_mod
  use gamd_mod
  use scaledMD_mod
  use ti_mod
  use phmd_mod
  ! This #ifdef section can be merged with the section on top
  use processor_mod !, only : proc_atm_crd, proc_atm_frc, proc_atm_qterm, &
  !       proc_num_atms
  use parallel_processor_mod
  use pme_recip_midpoint_mod
  use ensure_alloc_mod
#ifdef _OPENMP_
  use omp_lib
#endif /* _OPENMP_ */
  use nfe_setup_mod, only : nfe_on_force => on_force
  use nfe_lib_mod
  use neb_mod
  use nebread_mod
  use multipmemd_mod

#ifdef DBG_ARRAYS
  use dbg_arrays_mod
#endif

  implicit none

  ! Formal arguments:

  integer                       :: atm_cnt
  double precision              :: crd(3, atm_cnt)
  double precision              :: frc(3, atm_cnt)
  integer                       :: my_atm_lst(*)
  logical                       :: new_list
  logical                       :: need_pot_enes
  logical                       :: need_virials
  type(pme_pot_ene_rec)         :: pot_ene
  double precision              :: virial(3)            ! Only used for MD
  double precision              :: ekcmt(3)             ! Only used for MD
  double precision              :: pme_err_est          ! Only used for MD
  integer, intent(in)           :: nstep

  ! Local variables:

  type(pme_virial_rec)          :: proc_vir
  integer                       :: offset_loc
  logical                       :: firsttime=.true.
  logical                       :: send_data
  integer                       :: iatm, jatm, task_id
  type(pme_virial_rec)          :: vir
  double precision              :: enmr(3), entr
#ifndef MPI
#ifndef NOXRAY
  double precision              :: xray_e
#endif
#endif
  double precision              :: vir_vs_ene
  integer                       :: atm_lst_idx
  integer                       :: alloc_failed
  integer                       :: i, j, k
  logical                       :: params_may_change
  logical                       :: onstep
  double precision              :: net_frcs(3)

  double precision              :: evdwex, eelex
  double precision              :: enemap
  double precision              :: temp_tot_dih
  double precision              :: totdih
  integer                       :: buf_size

! GaMD Local Variables
  double precision              :: pot_ene_nb, pot_ene_lig, pot_ene_dih, fwgt
  double precision, allocatable :: nb_frc_gamd(:,:)
  double precision, allocatable :: nb_frc_tmp(:,:)
  double precision              :: net_frc_ti(3)

  double precision              :: bias_frc(3,atm_cnt)

  call zero_time()
  call zero_pme_time()
  if (.not. usemidpoint) then
    pot_ene = null_pme_pot_ene_rec
    if (ti_mode .ne. 0) call ti_zero_arrays
  end if
#ifdef CUDA
  if (nmropt .ne. 0) then
    call nmr_weight(atm_cnt, crd, mdout)
    params_may_change = .true.
  else
    params_may_change = .false.
  end if

  if (infe .ne. 0) params_may_change = .true.

  if (use_pme .ne. 0) then
    call self(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
      params_may_change)
    if (vdwmeth .eq. 1) then
      call vdw_correction(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
        params_may_change)
    end if
    call update_time(pme_misc_timer)
    if (need_pot_enes) then ! DO NOT TOUCH THIS, SEE RUNMD.F90 TO SET need_pot_enes
      if ((iamd .eq. 2) .or. (iamd .eq. 3)) then
        totdih = 0.0
        if (num_amd_lag .eq. 0) then
          temp_tot_dih = 0.0
          call gpu_calculate_amd_dihedral_energy(totdih)
          buf_size = 1
          call mpi_allreduce(totdih, temp_tot_dih, &
            buf_size, mpi_double_precision, &
            mpi_sum, pmemd_comm, err_code_mpi)
          totdih = temp_tot_dih
        end if
        call gpu_calculate_amd_dihedral_weight(totdih)
      end if
      if((igamd.eq.2).or.(igamd.eq.3).or.(igamd.eq.5).or.(igamd.eq.7).or.(igamd.eq.9))then
        totdih = 0.0
        if (num_gamd_lag .eq. 0) then
          temp_tot_dih = 0.0
          call gpu_calculate_gamd_dihedral_energy(totdih)
          buf_size = 1
          call mpi_allreduce(totdih, temp_tot_dih, &
                             buf_size, mpi_double_precision, &
                             mpi_sum, pmemd_comm, err_code_mpi)
          totdih = temp_tot_dih
        end if
        call gpu_calculate_gamd_dihedral_weight(totdih)
      end if

    end if

#ifdef GTI

    call gti_build_nl_list(ti_mode, lj1264, plj1264, reaf_mode) ! C4PairwiseCUDA
    call gti_clear(need_pot_enes, ti_mode, crd)

    call gti_ele_recip(need_pot_enes, need_virials, uc_volume, ti_mode, reaf_mode)
    call gti_nb(need_pot_enes, need_virials, ti_mode, lj1264, plj1264, reaf_mode) ! C4PairwiseCUDA

    call gti_bonded(need_pot_enes, need_virials, ti_mode, reaf_mode)
    call gti_others(need_pot_enes, nstep, dt, uc_volume, ew_coeff)

    call gti_finalize_force_virial(numextra, need_virials, ti_mode, netfrc, frc)
    if (ineb>0) then
      call transfer_fit_neb_crd(nstep)
    endif

    if (need_pot_enes .or. need_virials) call gti_update_md_ene(pot_ene, enmr, virial, ekcmt, ineb, nebfreq, nstep)

    if (ti_mode .gt. 0) then
      if (need_pot_enes) then
        call gti_update_ti_pot_energy_from_gpu(pot_ene, enmr, ti_mode, gti_cpu_output)
        if (ifmbar .ne. 0) call gti_update_mbar_from_gpu
      end if
      if (need_virials) call gti_get_virial(ti_weights(1:2), virial, ekcmt, ti_mode)
    end if

#else /* GTI */

    if (need_pot_enes) then
      if (ti_mode .eq. 0) then
        call gpu_pme_ene(ew_coeff, uc_volume, pot_ene, enmr, virial, &
                         ekcmt, nstep, dt, crd, frc)
      else
        call gpu_ti_pme_ene(ew_coeff, uc_volume, pot_ene, enmr, virial, &
                            ekcmt, nstep, dt, crd, frc, bar_cont)
      end if
    else
      call gpu_pme_force(ew_coeff, uc_volume, virial, ekcmt, nstep, dt, crd, frc)
    end if
    if (ineb>0) then
      call transfer_fit_neb_crd(nstep)
    endif
  else !(use_pme /= 0)
    call ipsupdate(ntb)
    if (need_pot_enes) then ! DO NOT TOUCH THIS, SEE RUNMD.F90 TO SET need_pot_enes
      if((iamd.eq.2).or.(iamd.eq.3))then
        call gpu_calculate_amd_dihedral_energy_weight()
      endif
      if((igamd.eq.2).or.(igamd.eq.3).or.(igamd.eq.5).or.(igamd.eq.7).or.(igamd.eq.9))then
        call gpu_calculate_gamd_dihedral_energy_weight()
      endif
      call gpu_ips_ene(uc_volume, pot_ene, enmr, virial, ekcmt)
  end if  !(use_pme /= 0)

  if (nmropt .ne. 0) then
    call nmr_calc(crd, frc, enmr, 6)
  end if
#ifndef GTI
  if (ti_mode .ne. 0) then
  !set these values after calling nmr_calc
    ti_ene_aug(1,ti_rest_dist_ene) = sc_pot_ene%sc_res_dist_R1
    ti_ene_aug(2,ti_rest_dist_ene) = sc_pot_ene%sc_res_dist_R2
    ti_ene_aug(1,ti_rest_ang_ene) = sc_pot_ene%sc_res_ang_R1
    ti_ene_aug(2,ti_rest_ang_ene) = sc_pot_ene%sc_res_ang_R2
    ti_ene_aug(1,ti_rest_tor_ene) = sc_pot_ene%sc_res_tors_R1
    ti_ene_aug(2,ti_rest_tor_ene) = sc_pot_ene%sc_res_tors_R2
  end if
#endif /* !GTI */

!same call as CPU
  if (ti_mode .ne. 0 .and. need_pot_enes .and. nstep .ge. 0) then
    call ti_calc_dvdl
  end if

  if(iamd.gt.0)then
!check if the energy here is updated
!AMD calculate weight and scale forces
    call gpu_calculate_and_apply_amd_weights(pot_ene%total, pot_ene%dihedral,&
                                             pot_ene%amd_boost,num_amd_lag)
     ! Update total energy
     pot_ene%total = pot_ene%total + pot_ene%amd_boost
  endif
  if(igamd.gt.0)then
!GAMD calculate weight and scale forces
    if(igamd.eq.1 .or. igamd.eq.2 .or. igamd.eq.3) then
       call gpu_calculate_and_apply_gamd_weights(pot_ene%total, pot_ene%dihedral, &
                                               pot_ene%gamd_boost,num_gamd_lag)
       ! Update total energy
       pot_ene%total = pot_ene%total + pot_ene%gamd_boost
     else if(igamd.eq.4 .or. igamd.eq.5 ) then
       pot_ene_nb = pot_ene%elec_tot + pot_ene%vdw_tot + pot_ene%vdw_14 + pot_ene%elec_14 + pot_ene%efield
       call gpu_calculate_and_apply_gamd_weights_nb(pot_ene_nb, pot_ene%dihedral, &
                                               pot_ene%gamd_boost,num_gamd_lag)
       ! Update total energy
       pot_ene%total = pot_ene%total + pot_ene%gamd_boost
     else if( (igamd.eq.6 .or. igamd.eq.7) .and. (ti_mode.ne.0)) then
       pot_ene_lig = ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                     ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                     ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
       pot_ene_dih = pot_ene%total + &
                      ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene)
       call gpu_calculate_and_apply_gamd_weights_ti_region(pot_ene_lig, pot_ene_dih, &
                                               pot_ene%gamd_boost,fwgt,num_gamd_lag)
       ! Update total energy
       ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
     else if( (igamd.eq.8 .or. igamd.eq.9) .and. (ti_mode.ne.0)) then
       pot_ene_lig = ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene)
       pot_ene_dih = pot_ene%total + &
                     ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene) + &
                     ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                     ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
       call gpu_calculate_and_apply_gamd_weights_ti_region(pot_ene_lig, pot_ene_dih, &
                                               pot_ene%gamd_boost,fwgt,num_gamd_lag)
       ! Update total energy
       ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
     else if( (igamd.eq.10 .or. igamd.eq.11) .and. (ti_mode.ne.0)) then
        pot_ene_lig = ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                      ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
        pot_ene_dih = pot_ene%total + &
                      ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                      ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene)
       call gpu_calculate_and_apply_gamd_weights_ti_region(pot_ene_lig, pot_ene_dih, &
                                               pot_ene%gamd_boost,fwgt,num_gamd_lag)
       ! Update total energy
       ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
     end if
  endif
  if(scaledMD.gt.0)then
!scaledMD scale forces
    call gpu_scaledmd_scale_frc(pot_ene%total)
    pot_ene%total = pot_ene%total * scaledMD_lambda
  endif

  call update_time(nonbond_time)

  if (infe .gt. 0) then
      if (nfe_first) then
         call mpi_bcast(nfe_atm_cnt, 1, MPI_INTEGER, 0, pmemd_comm, err_code_mpi)
         if (mytaskid .ne. 0) allocate(nfe_atm_lst(nfe_atm_cnt))
         call mpi_bcast(nfe_atm_lst, size(nfe_atm_lst), MPI_INTEGER, 0, pmemd_comm, err_code_mpi)
         call gpu_setup_shuttle_info(nfe_atm_cnt, 0, nfe_atm_lst)
         nfe_first = .false.
      end if

      call gpu_shuttle_retrieve_data(crd, 0)
      call gpu_shuttle_retrieve_data(frc, 1)

      bias_frc(:,:) = 0.d0
      call nfe_on_force(crd,bias_frc,pot_ene%nfe)

      if (mytaskid .eq. 0) then
         do i=1,nfe_atm_cnt
            j=nfe_atm_lst(i)

            frc(1,j) = frc(1,j)+bias_frc(1,j)
            frc(2,j) = frc(2,j)+bias_frc(2,j)
            frc(3,j) = frc(3,j)+bias_frc(3,j)
         enddo
      end if

      call mpi_bcast(frc, 3*atm_cnt, MPI_DOUBLE_PRECISION, 0, pmemd_comm, err_code_mpi)
      pot_ene%total = pot_ene%total + pot_ene%nfe
      pot_ene%restraint = pot_ene%restraint + pot_ene%nfe
      call gpu_shuttle_post_data(frc, 1)
  end if

#endif /* GTI  */

    call update_time(nonbond_time)
    if (need_virials) then
      vir%molecular(1,1) = virial(1)
      vir%molecular(2,2) = virial(2)
      vir%molecular(3,3) = virial(3)
    end if
    call dist_enes_virs_netfrcs(pot_ene, vir, ekcmt, net_frcs, &
                                need_pot_enes, need_virials)
    if (need_virials) then
      virial(1) = vir%molecular(1,1)
      virial(2) = vir%molecular(2,2)
      virial(3) = vir%molecular(3,3)
    end if

  else ! use_pme .ne. 0 below
    call update_time(pme_misc_timer)
    call ipsupdate(ntb)
    if (need_pot_enes) then ! DO NOT TOUCH THIS, SEE RUNMD.F90 TO SET need_pot_enes
      if ((iamd .eq. 2).or.(iamd .eq. 3)) then
        totdih = 0.0
        if (num_amd_lag .eq. 0) then
          temp_tot_dih = 0.0
          call gpu_calculate_amd_dihedral_energy(totdih)
          buf_size = 1
          call mpi_allreduce(totdih, temp_tot_dih, &
            buf_size, mpi_double_precision, &
            mpi_sum, pmemd_comm, err_code_mpi)
          totdih = temp_tot_dih
        end if
        call gpu_calculate_amd_dihedral_weight(totdih)
      end if
      if((igamd.eq.2).or.(igamd.eq.3).or.(igamd.eq.5).or.(igamd.eq.7).or.(igamd.eq.9))then
        totdih = 0.0
        if (num_gamd_lag .eq. 0) then
          temp_tot_dih = 0.0
          call gpu_calculate_gamd_dihedral_energy(totdih)
          buf_size = 1
          call mpi_allreduce(totdih, temp_tot_dih, &
            buf_size, mpi_double_precision, &
            mpi_sum, pmemd_comm, err_code_mpi)
          totdih = temp_tot_dih
        end if
        call gpu_calculate_gamd_dihedral_weight(totdih)
      end if
      call gpu_ips_ene(uc_volume, pot_ene, enmr, virial, ekcmt)
      call update_time(nonbond_time)
      if (need_virials) then
        vir%molecular(1,1) = virial(1)
        vir%molecular(2,2) = virial(2)
        vir%molecular(3,3) = virial(3)
      end if
      call dist_enes_virs_netfrcs(pot_ene, vir, ekcmt, net_frcs, &
        need_pot_enes, need_virials)
      if (need_virials) then
        virial(1) = vir%molecular(1,1)
        virial(2) = vir%molecular(2,2)
        virial(3) = vir%molecular(3,3)
      end if
    else
      call gpu_ips_force(uc_volume, virial, ekcmt)
      call update_time(nonbond_time)

      if (need_virials) then
        vir%molecular(1,1) = virial(1)
        vir%molecular(2,2) = virial(2)
        vir%molecular(3,3) = virial(3)
        call dist_enes_virs_netfrcs(pot_ene, vir, ekcmt, net_frcs, &
          need_pot_enes, need_virials)
        virial(1) = vir%molecular(1,1)
        virial(2) = vir%molecular(2,2)
        virial(3) = vir%molecular(3,3)
      end if
    end if
  end if
  if (nmropt .ne. 0) then
    call nmr_calc(crd, frc, enmr, 6)
  end if
  !AMD calculate weight and scale forces
  if (iamd .gt. 0) then
    call gpu_calculate_and_apply_amd_weights(pot_ene%total, pot_ene%dihedral, &
      pot_ene%amd_boost,num_amd_lag)
    ! Update total energy
    pot_ene%total = pot_ene%total + pot_ene%amd_boost
  end if
  !GaMD calculate weight and scale forces
  if (igamd .gt. 0) then
     if (igamd .eq. 1 .or. igamd .eq. 2 .or. igamd .eq. 3) then
      call gpu_calculate_and_apply_gamd_weights(pot_ene%total, pot_ene%dihedral, &
      pot_ene%gamd_boost,num_gamd_lag)
       ! Update total energy
       pot_ene%total = pot_ene%total + pot_ene%gamd_boost
     else if (igamd .eq. 4 .or. igamd .eq. 5 ) then
      pot_ene_nb = pot_ene%elec_tot + pot_ene%vdw_tot + pot_ene%vdw_14 + pot_ene%elec_14 + pot_ene%efield
      call gpu_calculate_and_apply_gamd_weights_nb(pot_ene_nb, pot_ene%dihedral, &
        pot_ene%gamd_boost,num_gamd_lag)
       ! Update total energy
       pot_ene%total = pot_ene%total + pot_ene%gamd_boost
     else if( (igamd.eq.6 .or. igamd.eq.7) .and. (ti_mode.ne.0)) then
       pot_ene_lig = ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                     ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                     ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
       pot_ene_dih = pot_ene%total + &
                      ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene)
       call calculate_gamd_ti_region_weights(atm_cnt,pot_ene_lig, pot_ene_dih, &
                       pot_ene%gamd_boost,fwgt,frc,crd,ti_sc_lst,my_atm_lst)
       ! Calculate net force on TI region
       if (mod(nstep+1,ntave).eq.0) then
         net_frc_ti    = 0.0
         do i = 1, atm_cnt
               if ( ti_sc_lst(i) .gt. 0 ) then
                       net_frc_ti(:) = net_frc_ti(:) + frc(:,i)
               end if
         end do
         write(*,'(a,i10,4f15.5)') "| GaMD: step, frc_weight, frc_lig(1:3) = ", nstep+1, fwgt, net_frc_ti
       end if
       ! Update total energy
       ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
     else if( (igamd.eq.8 .or. igamd.eq.9) .and. (ti_mode.ne.0)) then
       pot_ene_lig = ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene)
       pot_ene_dih = pot_ene%total + &
                     ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene) + &
                     ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                     ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
       call calculate_gamd_ti_region_weights(atm_cnt,pot_ene_lig, pot_ene_dih, &
                       pot_ene%gamd_boost,fwgt,frc,crd,ti_sc_lst,my_atm_lst)
       ! Update total energy
       ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
     else if( (igamd.eq.10 .or. igamd.eq.11) .and. (ti_mode.ne.0)) then
       pot_ene_lig = ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                     ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
       pot_ene_dih = pot_ene%total + &
                     ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                     ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene)
       call calculate_gamd_ti_region_weights(atm_cnt,pot_ene_lig, pot_ene_dih, &
                       pot_ene%gamd_boost,fwgt,frc,crd,ti_sc_lst,my_atm_lst)
       ! Update total energy
       ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
    end if
  end if
  !scaledMD calculate scale forces
  if (scaledMD .gt. 0) then
    call gpu_scaledmd_scale_frc(pot_ene%total)
    pot_ene%total = pot_ene%total * scaledMD_lambda
  end if

  call update_time(fcve_dist_time)

  if (infe .gt. 0) then
      if (nfe_first) then
         call mpi_bcast(nfe_atm_cnt, 1, MPI_INTEGER, 0, pmemd_comm, err_code_mpi)
         if (mytaskid .ne. 0) allocate(nfe_atm_lst(nfe_atm_cnt))
         call mpi_bcast(nfe_atm_lst, size(nfe_atm_lst), MPI_INTEGER, 0, pmemd_comm, err_code_mpi)
         call gpu_setup_shuttle_info(nfe_atm_cnt, 0, nfe_atm_lst)
         nfe_first = .false.
      end if

      call gpu_shuttle_retrieve_data(crd, 0)
      call gpu_shuttle_retrieve_data(frc, 1)

      bias_frc(:,:) = 0.d0
      call nfe_on_force(crd,bias_frc,pot_ene%nfe)

      if (mytaskid .eq. 0) then
         do i=1,nfe_atm_cnt
            j=nfe_atm_lst(i)

            frc(1,j) = frc(1,j)+bias_frc(1,j)
            frc(2,j) = frc(2,j)+bias_frc(2,j)
            frc(3,j) = frc(3,j)+bias_frc(3,j)
         enddo
      end if

      call mpi_bcast(frc, 3*atm_cnt, MPI_DOUBLE_PRECISION, 0, pmemd_comm, err_code_mpi)
      pot_ene%total = pot_ene%total + pot_ene%nfe
      pot_ene%restraint = pot_ene%restraint + pot_ene%nfe
      call gpu_shuttle_post_data(frc, 1)
  end if

  if (ineb .gt. 0) then
    call full_neb_forces(nstep)
    if (beadid==1 .or. beadid==neb_nbead) then
      call gpu_clear_forces()
    endif
    end if

#else /*Not cuda version below*/

  if (.not. usemidpoint) then
    call do_load_balancing(new_list, atm_cnt)

    ! Do weight changes, if requested.

    if (nmropt .ne. 0) call nmr_weight(atm_cnt, crd, mdout)
  end if

  ! If no ene / force calcs are to be done, clear everything and bag out now.
  if (ntf .eq. 8) then
    pot_ene = null_pme_pot_ene_rec
    virial(:) = 0.d0
    ekcmt(:) = 0.d0
    pme_err_est = 0.d0
    frc(:,:) = 0.d0
    return
  end if

  if (.not. usemidpoint) then
    if (igamd.eq.4 .or. igamd.eq.5) then
      allocate( nb_frc_gamd(3, atm_cnt), &
        nb_frc_tmp(3, atm_cnt), &
        stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
    end if
  end if ! not usemidpoint

  ! Calculate the non-bonded contributions:

  ! Direct part of ewald plus vdw, hbond, pairlist setup and image claiming:

  if (ntp .gt. 0) call fill_tranvec(gbl_tranvec)
  ! The following encapsulation (save_imgcrds) seems to be necessary to
  ! prevent an optimization bug with the SGI f90 compiler.  Sigh...

  if (new_list) then
    call pme_list_midpoint(atm_cnt, crd, atm_nb_maskdata, atm_nb_mask &
      , proc_atm_nb_maskdata, proc_atm_nb_mask,nstep)
  end if

  ! Zero energies that are stack or call parameters:

  pot_ene = null_pme_pot_ene_rec
  vir = null_pme_virial_rec

  virial(:) = 0.d0
  ekcmt(:) = 0.d0
  pme_err_est = 0.d0

  ! Zero internal energies, virials, etc.

  enmr(:) = 0.d0
  vir_vs_ene = 0.d0

  net_frcs(:) = 0.d0
  molvir_netfrc_corr(:,:) = 0.d0

  if (usemidpoint) then
    !Resetting the forces
    call ensure_alloc_dble3(nb_frc, proc_num_atms_min_bound+proc_ghost_num_atms)
    do j = 1, proc_num_atms+proc_ghost_num_atms
      proc_atm_frc(:, j) = 0.d0
      nb_frc(:,j)=0.0d0
    end do
  else
    do j = 1, gbl_used_img_cnt
      i = gbl_used_img_lst(j)
      img_frc(:, i) = 0.d0
    end do

    if (ti_mode .ne. 0) then
      call ti_zero_arrays
      do j = 1, gbl_used_img_cnt
        i = gbl_used_img_lst(j)
        ti_img_frc(:, :, i) = 0.d0
      end do
      call zero_extra_used_atm_ti_img_frcs(ti_img_frc)
    end if
  end if

  ! We also need to zero anything on the extra used atom list since nonbonded
  ! forces for it will get sent to the atom owner.  We don't calc any such
  ! forces for these atoms, but they are on the send atom list in order to
  ! get their coordinates updated.

  if (.not. usemidpoint) then
    call zero_extra_used_atm_img_frcs(img_frc)
  end if

  ! Don't do recip if PME is not invoked. Don't do it this step unless
  ! mod(irespa,nrepsa) = 0

  onstep = mod(irespa, nrespa) .eq. 0

  params_may_change = (nmropt .ne. 0) .or. (infe .ne. 0)

  ! Calculate the exclusions correction for ips if using it
  ! Clear the force array and the atom-ordered nonbonded force array.
  ! We delay copying the nonbonded forces into the force
  ! array in order to be able to schedule i/o later, and batch it up.

  evdwex = 0.0
  eelex = 0.0

  if (use_pme .ne. 0 .and. onstep) then

    ! Self energy:

    ! The small amount of time used here gets lumped with the recip stuff...

    if (iphmd .eq. 3) then
      call runphmd(crd,frc,atm_iac,typ_ico,atm_numex,gbl_natex,belly_atm_cnt)
    end if
    if (usemidpoint) then
      !!!?? Need to verify if the following version is ok or not
      if (master) then
        call self(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
          params_may_change)
      end if
    else
      if (gbl_frc_ene_task_lst(1) .eq. mytaskid) then
        call self(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
          params_may_change)
      end if
    end if


    ! Reciprocal energy:

    ! We intentionally keep the load balance counter running through the
    ! fft dist transposes in the recip code; synchronization will mess up the
    ! times a bit, but when not all nodes are doing recip calcs, it will
    ! really help with load balancing.
    !shall we use min max here?

    if (usemidpoint) then
      call update_pme_time(pme_misc_timer)
#if defined(pmemd_SPDP)
      do i = 1 , proc_num_atms + proc_ghost_num_atms ! my_img_lo, my_img_hi
        proc_crd_q_sp(1:3,i) = proc_atm_crd(1:3,i)
        proc_crd_q_sp(4,i) = proc_atm_qterm(i)
      end do
#ifdef _OPENMP_
      call update_pme_time(get_nb_pack_timer)
#endif
#endif /* pmemd_SPDP*/
      ! Midpoint just started  supporting
      !print *, "before Recip new energy",pot_ene%elec_recip,vir%elec_recip,mytaskid
      call do_kspace_midpoint(nb_frc,pot_ene%elec_recip, vir%elec_recip, need_pot_enes, need_virials)
      !DBG_ARRAYS_DUMP_3DBLE("frc_fft", proc_atm_to_full_list, proc_atm_frc, proc_num_atms+proc_ghost_num_atms)
      !print *, "after Recip new energy",pot_ene%elec_recip,vir%elec_recip,mytaskid

    else ! not usemidpoint below
      if (i_do_recip) then
        if (ti_mode .eq. 0) then
          call update_pme_time(pme_misc_timer)
          call update_loadbal_timer(elapsed_100usec_atmuser)
          if (block_fft .eq. 0) then
            call do_slab_pmesh_kspace(img_frc, pot_ene%elec_recip, &
              vir%elec_recip, need_pot_enes, need_virials)
          else
            call do_blk_pmesh_kspace(img_frc, pot_ene%elec_recip, &
              vir%elec_recip, need_pot_enes, need_virials)
          end if
          call update_loadbal_timer(elapsed_100usec_recipfrc)
          if (nrespa .gt. 1) call respa_scale(atm_cnt, img_frc, nrespa)
        else
          call update_pme_time(pme_misc_timer)
          call update_loadbal_timer(elapsed_100usec_atmuser)
          ti_mask_piece = 2
          if (block_fft .eq. 0) then
            call do_slab_pmesh_kspace(img_frc, pot_ene%elec_recip, ti_vir(1,:,:), need_pot_enes, need_virials)

          else
            call do_blk_pmesh_kspace(img_frc, pot_ene%elec_recip, &
              ti_vir(1,:,:), need_pot_enes, need_virials)
          end if
          ! It is much more efficient to do it this way, rather
          ! than passing ti_img_frc(1,:,:) above ...
          do j = 1, gbl_used_img_cnt
            i = gbl_used_img_lst(j)
            ti_img_frc(1, :, i) = img_frc(:, i)
            img_frc(:,i) = 0.d0
          end do

          ti_pot(1) = pot_ene%elec_recip

          pot_ene%elec_recip = 0.d0
          ti_mask_piece = 1
          if (ti_mode .ne. 3) then
            if (block_fft .eq. 0) then
              call do_slab_pmesh_kspace(img_frc, pot_ene%elec_recip, ti_vir(2,:,:), need_pot_enes, need_virials)
            else
              call do_blk_pmesh_kspace(img_frc, pot_ene%elec_recip, &
                ti_vir(2,:,:), need_pot_enes, need_virials)
            end if
            do j = 1, gbl_used_img_cnt
              i = gbl_used_img_lst(j)
              ti_img_frc(2, :, i) = img_frc(:, i)
              img_frc(:,i) = 0.d0
            end do
            ti_pot(2) = pot_ene%elec_recip
          else
            ti_pot(2) = ti_pot(1)
            do j = 1, gbl_used_img_cnt
              i = gbl_used_img_lst(j)
              ti_img_frc(2, :, i) = ti_img_frc(1, :, i)
            end do
            ti_vir(2,:,:) = ti_vir(1,:,:)
          end if
          call update_loadbal_timer(elapsed_100usec_recipfrc)

          do j = 1, gbl_used_img_cnt
            i = gbl_used_img_lst(j)
            ti_img_frc(1, :, i) = ti_img_frc(1, :, i) * ti_weights(1)
            ti_img_frc(2, :, i) = ti_img_frc(2, :, i) * ti_weights(2)
          end do

          vir%elec_recip = ti_vir(1,:,:) * ti_weights(1) + &
            ti_vir(2,:,:) * ti_weights(2)

          ti_vve(ti_vir0) = (ti_vir(1,1,1) + &
            ti_vir(1,2,2) + &
            ti_vir(1,3,3)) * ti_weights(1)
          ti_vve(ti_vir1) = (ti_vir(2,1,1) + &
            ti_vir(2,2,2) + &
            ti_vir(2,3,3)) * ti_weights(2)

          call ti_update_ene_all(ti_pot, si_elect_ene, pot_ene%elec_recip,3)

          if (nrespa .gt. 1) then
            call ti_respa_scale(atm_cnt, ti_img_frc, nrespa)
          end if
        end if
      end if
    end if ! usemidpoint


    ! Long range dispersion contributions:

    ! Continuum method:

    if (usemidpoint) then
      if (vdwmeth .eq. 1 .and. master) then
        call vdw_correction(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
          params_may_change)
      end if
    else
      if (vdwmeth .eq. 1 .and. gbl_frc_ene_task_lst(1) .eq. mytaskid) then
        call vdw_correction(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
          params_may_change)
      end if
    end if
  end if      ! respa .and. use_pme


  if (.not. usemidpoint) then
    call update_loadbal_timer(elapsed_100usec_atmuser)
  end if

  ! Direct part of ewald plus vdw, hbond, force and energy calculations:

  call update_pme_time(pme_misc_timer)
  if (ips .gt. 0) then
    call get_nb_ips_energy(img_frc, gbl_img_crd, gbl_img_qterm, gbl_eed_cub, &
                           gbl_ipairs, gbl_tranvec, need_pot_enes, need_virials, &
                           pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                           vir%eedvir, vir%elec_direct)
  else if (ti_mode .ne. 0) then
    ! This is the only part of the TI code that uses img_frc, everything else
    ! writes to ti_img_frc
    call get_nb_energy_ti(img_frc, gbl_img_crd, gbl_img_qterm, gbl_eed_cub, &
                          gbl_ipairs, gbl_tranvec, need_pot_enes, need_virials, &
                          pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                          vir%eedvir, vir%elec_direct, gbl_img_atm_map)
  else if ((lj1264 .eq. 1) .and. (plj1264 .eq. 0)) then
    call get_nb_energy_1264(img_frc, gbl_img_crd, gbl_img_qterm, gbl_eed_cub, &
                            gbl_ipairs, gbl_tranvec, need_pot_enes, need_virials, &
                            pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                            vir%eedvir, vir%elec_direct)
  else if ((lj1264 .eq. 0) .and. (plj1264 .eq. 1)) then ! New2022
    call get_nb_energy_p1264(img_frc, gbl_img_crd, gbl_img_qterm, gbl_eed_cub, &
                              gbl_ipairs, gbl_tranvec, need_pot_enes, need_virials, &
                              pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                              vir%eedvir, vir%elec_direct)
  else if ((lj1264 .eq. 1) .and. (plj1264 .eq. 1)) then ! New2022
    call get_nb_energy_1264p1264(img_frc, gbl_img_crd, gbl_img_qterm, gbl_eed_cub, &
                              gbl_ipairs, gbl_tranvec, need_pot_enes, need_virials, &
                              pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                              vir%eedvir, vir%elec_direct)
  else

    if (usemidpoint) then
#ifdef _OPENMP_
#ifdef pmemd_SPDP
      call get_nb_energy_midpoint(frc_thread, proc_crd_q_sp, gbl_eed_cub, &
                                  proc_ipairs, need_pot_enes, need_virials, &
                                  pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                                  vir%eedvir, vir%elec_direct)
#else
      call get_nb_energy_midpoint(frc_thread, proc_atm_crd, proc_atm_qterm, gbl_eed_cub, &
                                  proc_ipairs, need_pot_enes, need_virials, &
                                  pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                                  vir%eedvir, vir%elec_direct)
#endif /* pmemd_SPDP */
#else
#ifdef pmemd_SPDP
      call get_nb_energy_midpoint(nb_frc, proc_crd_q_sp, gbl_eed_cub, &
                                  proc_ipairs, need_pot_enes, need_virials, &
                                  pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                                  vir%eedvir, vir%elec_direct)
#else
      call get_nb_energy_midpoint(nb_frc, proc_atm_crd, proc_atm_qterm, gbl_eed_cub, &
                                  proc_ipairs, need_pot_enes, need_virials, &
                                  pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                                  vir%eedvir, vir%elec_direct)
#endif
#endif /* _OPENMP_ */
      !DBG_ARRAYS_DUMP_3DBLE("frc_nbe", proc_atm_to_full_list, proc_atm_frc, proc_num_atms+proc_ghost_num_atms)
    else
      call get_nb_energy(img_frc, gbl_img_crd, gbl_img_qterm, gbl_eed_cub, &
                         gbl_ipairs, gbl_tranvec, need_pot_enes, need_virials, &
                         pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                         vir%eedvir, vir%elec_direct)
    end if
  end if

  if (usemidpoint) then
#ifdef _OPENMP_
    do j = 0, nthreads-1
      offset_loc = j*(proc_num_atms + proc_ghost_num_atms)
      do k = 1,  proc_num_atms + proc_ghost_num_atms
        nb_frc(1:3,k) = nb_frc(1:3,k) + frc_thread(1:3,k+offset_loc)
        frc_thread(1:3,k+offset_loc) = 0.d0
      end do  ! k = 1,  proc_num_atms + proc_ghost_num_atms
    end do  ! do j=0,nthreads-1
#endif /* _OPENMP_ */
  end if

  call update_pme_time(dir_frc_sum_timer)
  if (.not. usemidpoint) then
    if (efx .ne. 0 .or. efy .ne. 0 .or. efz .ne. 0) then
      call get_efield_energy(img_frc, crd, gbl_img_qterm, gbl_img_atm_map, need_pot_enes, &
        pot_ene%efield, atm_cnt, nstep)
    end if
  end if

  call update_pme_time(dir_frc_sum_timer)
  if (.not. usemidpoint) then
    call update_loadbal_timer(elapsed_100usec_dirfrc)
  end if

  call update_pme_time(pme_misc_timer)

  ! Calculate 1-4 electrostatic energies, forces:
  if (usemidpoint) then
    if (need_virials) then
      call get_nb14_energy_midpoint(proc_atm_qterm, proc_atm_crd, nb_frc, proc_iac, typ_ico, &
                                    gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt, &
                                    pot_ene%elec_14, pot_ene%vdw_14, &
                                    vir%elec_14)
    else
      call get_nb14_energy_midpoint(proc_atm_qterm, proc_atm_crd, nb_frc, proc_iac, typ_ico, &
                                    gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt, &
                                    pot_ene%elec_14, pot_ene%vdw_14)
    end if
  else
    if (charmm_active) then
      if (need_virials) then
        call get_nb14_energy(atm_qterm, crd, nb_frc, atm_iac, typ_ico, &
                             gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt, &
                             pot_ene%elec_14, pot_ene%vdw_14, &
                             vir%elec_14)
      else
        call get_nb14_energy(atm_qterm, crd, nb_frc, atm_iac, typ_ico, &
                             gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt, &
                             pot_ene%elec_14, pot_ene%vdw_14)
      end if
    else
      if (need_virials) then
        call get_nb14_energy(atm_qterm, crd, nb_frc, atm_iac, typ_ico, &
                             gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt, &
                             pot_ene%elec_14, pot_ene%vdw_14, &
                             vir%elec_14)
      else
        call get_nb14_energy(atm_qterm, crd, nb_frc, atm_iac, typ_ico, &
                             gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt, &
                             pot_ene%elec_14, pot_ene%vdw_14)
      end if
    end if
  end if

  call update_pme_time(dir_frc_sum_timer)

  ! Adjust energies, forces for masked out pairs:

  if (use_pme .ne. 0) then
    if (usemidpoint) then
      call nb_adjust_midpoint(proc_num_atms+proc_ghost_num_atms, proc_atm_qterm, &
                              proc_atm_crd, gbl_nb_adjust_pairlst, gbl_eed_cub, &
                              nb_frc, pot_ene%elec_nb_adjust, vir%elec_nb_adjust)
    else
      call nb_adjust(my_atm_cnt, my_atm_lst, &
                     atm_qterm, crd, gbl_nb_adjust_pairlst, gbl_eed_cub, &
                     nb_frc, pot_ene%elec_nb_adjust, vir%elec_nb_adjust)
    end if
  end if

  call update_pme_time(adjust_masked_timer)

  ! Get COM-relative coordinates.
  ! Calculate total nonbonded energy components.

  pot_ene%vdw_tot = pot_ene%vdw_dir + &
                    pot_ene%vdw_recip

  pot_ene%elec_tot = pot_ene%elec_dir + &
                     pot_ene%elec_recip + &
                     pot_ene%elec_nb_adjust + &
                     pot_ene%elec_self

  if ( ips .gt. 0 ) then
    pot_ene%vdw_tot = pot_ene%vdw_tot + evdwex
    pot_ene%elec_tot = pot_ene%elec_tot + eelex
  end if

  if (need_virials) then
    if (usemidpoint) then
      call ensure_alloc_dble3(atm_rel_crd, proc_num_atms_min_bound+proc_ghost_num_atms)
      call get_atm_rel_crd_midpoint(gbl_mol_cnt, gbl_mol_com, proc_atm_crd, atm_rel_crd, pbc_box)
    else
      call get_atm_rel_crd(my_mol_cnt, gbl_mol_com, crd, atm_rel_crd, gbl_my_mol_lst)
    end if
  end if

  call update_pme_time(pme_misc_timer)
  call update_time(nonbond_time)

  if (usemidpoint) then
    if (igamd.eq.4 .or. igamd.eq.5) then
      nb_frc_gamd = frc
    end if
  end if
  ! Calculate the other contributions:

  ! Global to Local atomic index
  ! In order to get local atomic index from global(i.e.as it is written in inpcrd)
  ! index
  !!!NOT sure if we need here since we already have it in build_atm_space
  !routine
  if (usemidpoint) then
    call pme_bonded_force_midpoint(pot_ene, new_list)
    ! DBG_ARRAYS_DUMP_3DBLE("frc_bond", proc_atm_to_full_list, proc_atm_frc, proc_num_atms+proc_ghost_num_atms)

    do i=1,proc_num_atms + proc_ghost_num_atms
      proc_atm_frc(:,i) = proc_atm_frc(:,i) + nb_frc(:,i)
    end do
  else
    call pme_bonded_force(crd, frc, pot_ene)
  end if

  ! The above stuff gets lumped as "other time"...

  ! --- calculate the EMAP constraint energy ---

  if (.not. usemidpoint) then
    if (iemap .gt. 0) then
      call emapforce(natom,enemap,atm_mass,crd,frc )
      pot_ene%emap = enemap
      pot_ene%restraint = pot_ene%restraint + enemap
    end if
  end if

  ! Sum up total potential energy for this task:

  pot_ene%total = pot_ene%vdw_tot + &
                  pot_ene%elec_tot + &
                  pot_ene%hbond + &
                  pot_ene%bond + &
                  pot_ene%angle + &
                  pot_ene%dihedral + &
                  pot_ene%vdw_14 + &
                  pot_ene%elec_14 + &
                  pot_ene%restraint + &
                  pot_ene%imp + &
                  pot_ene%angle_ub + &
                  pot_ene%efield + &
                  pot_ene%cmap

  if (.not. usemidpoint) then
    call update_loadbal_timer(elapsed_100usec_atmowner)

    !add the ips non bonded contribution
    if ( ips .gt. 0 ) then
      do atm_lst_idx = 1, my_atm_cnt
        i = my_atm_lst(atm_lst_idx)
        nb_frc(:,i) = nb_frc(:,i) + ips_frc(:,i)
      end do
    end if
  end if

  if (usemidpoint) then
    ! Return forces to local mpi for collection

    call comm_frc(proc_atm_frc)
    call comm_frc(nb_frc)

  else
    if (new_list) call get_img_frc_distribution(atm_cnt)

    if (ti_mode .eq. 0) then
      call distribute_img_frcs(atm_cnt, img_frc, nb_frc)
    else
      call distribute_img_ti_frcs(atm_cnt, ti_img_frc, ti_nb_frc)
    end if
  end if

  call update_time(fcve_dist_time)

  call zero_pme_time()

  ! If using extra points and a frame (checked internal to subroutine),
  ! transfer force and torque from the extra points to the parent atom:

  if (.not. usemidpoint) then
    if (numextra .gt. 0 .and. frameon .ne. 0) then
      if (ti_mode .eq. 0) then
        call orient_frc(crd, nb_frc, vir%ep_frame, ep_frames, ep_lcl_crd, gbl_frame_cnt)
      else
        ti_vir(:,:,:) = 0.d0
        ! Use nb_frc as a temporary array, so we don't need to splice
        ! the entire force array (ti_nb_frc)
        do atm_lst_idx = 1, my_atm_cnt
          i = my_atm_lst(atm_lst_idx)
          nb_frc(:,i) = ti_nb_frc(1, :, i)
        end do
        call orient_frc(crd, nb_frc, ti_vir(1,:,:), ep_frames, ep_lcl_crd, gbl_frame_cnt)
        do atm_lst_idx = 1, my_atm_cnt
          i = my_atm_lst(atm_lst_idx)
          ti_nb_frc(1, :, i) = nb_frc(:,i)
          nb_frc(:,i) = ti_nb_frc(2, :, i)
        end do
        call orient_frc(crd, nb_frc, ti_vir(2,:,:), ep_frames, ep_lcl_crd, gbl_frame_cnt)
        do atm_lst_idx = 1, my_atm_cnt
          i = my_atm_lst(atm_lst_idx)
          ti_nb_frc(2, :, i) = nb_frc(:,i)
          nb_frc(:,i) = 0.d0
        end do
        vir%ep_frame = ti_vir(1,:,:) + ti_vir(2,:,:)
      end if
    end if
  end if ! not usemidpoint
  ! If the net force correction is in use, here we determine the net forces by
  ! looking at the sum of all nonbonded forces.  This should give the same
  ! result as just looking at the reciprocal forces, but it is more
  ! computationally convenient, especially for extra points, to do it this way.

  if (netfrc .gt. 0 .and. onstep) then
    if (ti_mode .eq. 0) then
      if (usemidpoint) then
        do i = 1, proc_num_atms
          net_frcs(:) = net_frcs(:) + nb_frc(:, i)
        end do
      else
        do atm_lst_idx = 1, my_atm_cnt
          i = my_atm_lst(atm_lst_idx)
          net_frcs(:) = net_frcs(:) + nb_frc(:, i)
        end do
      end if
    else
      do atm_lst_idx = 1, my_atm_cnt
        i = my_atm_lst(atm_lst_idx)
        ti_net_frcs(1, :) = ti_net_frcs(1, :) + ti_nb_frc(1, :, i)
        ti_net_frcs(2, :) = ti_net_frcs(2, :) + ti_nb_frc(2, :, i)
      end do
    end if
  end if

  ! First phase of virial work.  We need just the nonbonded forces at this
  ! stage.

  if (need_virials) then

    ! The relative coordinates for the extra points will mess up the molvir
    ! netfrc correction here unless we zero them out.  They are of no import
    ! in ekcom calc because extra points have a mass of 0.

    if (usemidpoint) then
      do atm_lst_idx = 1, proc_num_atms
        vir%molecular(:,1) = vir%molecular(:,1) + nb_frc(:,atm_lst_idx) * atm_rel_crd(1,atm_lst_idx)
        vir%molecular(:,2) = vir%molecular(:,2) + nb_frc(:,atm_lst_idx) * atm_rel_crd(2,atm_lst_idx)
        vir%molecular(:,3) = vir%molecular(:,3) + nb_frc(:,atm_lst_idx) * atm_rel_crd(3,atm_lst_idx)
        molvir_netfrc_corr(:, 1) = molvir_netfrc_corr(:, 1) + atm_rel_crd(1,atm_lst_idx)
        molvir_netfrc_corr(:, 2) = molvir_netfrc_corr(:, 2) + atm_rel_crd(2,atm_lst_idx)
        molvir_netfrc_corr(:, 3) = molvir_netfrc_corr(:, 3) + atm_rel_crd(3,atm_lst_idx)
      end do
    else
      if (numextra .gt. 0 .and. frameon .ne. 0) &
        call zero_extra_pnts_vec(atm_rel_crd, ep_frames, gbl_frame_cnt)

      if (ti_mode .eq. 0) then
        do atm_lst_idx = 1, my_atm_cnt
          i = my_atm_lst(atm_lst_idx)
          vir%molecular(:,1) = vir%molecular(:,1) + nb_frc(:,i) * atm_rel_crd(1,i)
          vir%molecular(:,2) = vir%molecular(:,2) + nb_frc(:,i) * atm_rel_crd(2,i)
          vir%molecular(:,3) = vir%molecular(:,3) + nb_frc(:,i) * atm_rel_crd(3,i)
          molvir_netfrc_corr(:, 1) = molvir_netfrc_corr(:, 1) + atm_rel_crd(1, i)
          molvir_netfrc_corr(:, 2) = molvir_netfrc_corr(:, 2) + atm_rel_crd(2, i)
          molvir_netfrc_corr(:, 3) = molvir_netfrc_corr(:, 3) + atm_rel_crd(3, i)
        end do
      else
        do atm_lst_idx = 1, my_atm_cnt
          i = my_atm_lst(atm_lst_idx)

          vir%molecular(:,1) = vir%molecular(:,1) + &
                               ti_nb_frc(1,:,i) * atm_rel_crd_sc(1,1,i) + &
                               ti_nb_frc(2,:,i) * atm_rel_crd_sc(2,1,i)
          vir%molecular(:,2) = vir%molecular(:,2) + &
                               ti_nb_frc(1,:,i) * atm_rel_crd_sc(1,2,i) + &
                               ti_nb_frc(2,:,i) * atm_rel_crd_sc(2,2,i)
          vir%molecular(:,3) = vir%molecular(:,3) + &
                               ti_nb_frc(1,:,i) * atm_rel_crd_sc(1,3,i) + &
                               ti_nb_frc(2,:,i) * atm_rel_crd_sc(2,3,i)

          ti_molvir_netfrc_corr(1,:,1) = ti_molvir_netfrc_corr(1,:,1) + &
                                         atm_rel_crd_sc(1,1,i)
          ti_molvir_netfrc_corr(1,:,2) = ti_molvir_netfrc_corr(1,:,2) + &
                                         atm_rel_crd_sc(1,2,i)
          ti_molvir_netfrc_corr(1,:,3) = ti_molvir_netfrc_corr(1,:,3) + &
                                         atm_rel_crd_sc(1,3,i)

          ti_molvir_netfrc_corr(2,:,1) = ti_molvir_netfrc_corr(2,:,1) + &
                                         atm_rel_crd_sc(2,1,i)
          ti_molvir_netfrc_corr(2,:,2) = ti_molvir_netfrc_corr(2,:,2) + &
                                         atm_rel_crd_sc(2,2,i)
          ti_molvir_netfrc_corr(2,:,3) = ti_molvir_netfrc_corr(2,:,3) + &
                                         atm_rel_crd_sc(2,3,i)
        end do
      end if
    end if
    ! Finish up virial work; Timing is inconsequential...

    vir%atomic(:,:) = vir%elec_recip(:,:) + &
      vir%elec_direct(:,:) + &
      vir%elec_nb_adjust(:,:) + &
      vir%elec_recip_vdw_corr(:,:) + &
      vir%elec_recip_self(:,:) + &
      vir%elec_14(:,:) + &
      vir%ep_frame(:,:)

    !add the ips contribution to virials
    if ( ips .gt. 0 ) then
      vir%atomic(:,:) = vir%atomic(:,:) + VIREXIPS(:,:)
    end if

    vir%molecular(:,:) = vir%molecular(:,:) + vir%atomic(:,:)

    if (usemidpoint) then
      call get_ekcom_midpoint(gbl_mol_cnt, gbl_mol_mass_inv, ekcmt, proc_atm_vel, &
                              proc_atm_mass)
    else
      call get_ekcom(my_mol_cnt, gbl_mol_mass_inv, ekcmt, atm_vel, &
                     atm_mass, gbl_my_mol_lst)
    end if
  end if

  call update_time(nonbond_time)
  call update_pme_time(pme_misc_timer)

  ! Add a bunch of values (energies, virials, etc.) together as needed.

  call dist_enes_virs_netfrcs(pot_ene, vir, ekcmt, net_frcs, &
    need_pot_enes, need_virials)

  if (ti_mode .ne. 0 ) then
    call ti_dist_enes_virs_netfrcs(onstep, need_virials)
    call ti_combine_frcs(atm_cnt, my_atm_cnt, my_atm_lst, nb_frc, ti_nb_frc)
  end if

  call update_time(fcve_dist_time)
  call zero_pme_time()

  if (.not. usemidpoint) then
    if (iamd.gt.1) then
      call calculate_amd_dih_weights(atm_cnt,pot_ene%dihedral,frc,crd)
    end if

    if((igamd.eq.2).or.(igamd.eq.3).or.(igamd.eq.5).or.(igamd.eq.7).or.(igamd.eq.9))then
      call calculate_gamd_dih_weights(atm_cnt,pot_ene%dihedral,frc,crd)
    end if
  end if

  ! Copy image forces to atom forces, correcting for netfrc as you go, if
  ! appropriate.  Do not remove net force if netfrc = 0; e.g. in minimization.

  if (netfrc .gt. 0 .and. onstep) then

    if (ti_mode .eq. 0) then
      if (usemidpoint) then
        net_frcs(:) = net_frcs(:) / dble(atm_cnt - numextra)

        do i = 1, proc_num_atms
          proc_atm_frc(:, i) = proc_atm_frc(:, i) - net_frcs(:)
        end do
      else
        net_frcs(:) = net_frcs(:) / dble(atm_cnt - numextra)

        do atm_lst_idx = 1, my_atm_cnt
          i = my_atm_lst(atm_lst_idx)
          frc(:, i) = frc(:, i) + nb_frc(:, i) - net_frcs(:)
        end do
      end if
    else
      ti_net_frcs(1,:) = ti_net_frcs(1,:)/dble(ti_atm_cnt(1)-ti_numextra_pts(1))
      ti_net_frcs(2,:) = ti_net_frcs(2,:)/dble(ti_atm_cnt(2)-ti_numextra_pts(2))
      net_frcs(:) = ti_net_frcs(1,:) + ti_net_frcs(2,:)

      do atm_lst_idx = 1, my_atm_cnt
        i = my_atm_lst(atm_lst_idx)
        ! This matches how sander removes netfrcs in TI runs
        if (ti_lst(1,i) .ne. 0) then
          frc(:, i) = frc(:, i) + nb_frc(:, i) - ti_net_frcs(1,:)
        else if (ti_lst(2,i) .ne. 0) then
          frc(:, i) = frc(:, i) + nb_frc(:, i) - ti_net_frcs(2,:)
        else
          frc(:, i) = frc(:, i) + nb_frc(:, i) - net_frcs(:)
        end if
      end do
    end if

    if (numextra .gt. 0 .and. frameon .ne. 0) &
      call zero_extra_pnts_vec(frc, ep_frames, gbl_frame_cnt)

    ! Correct the molecular virial for netfrc:

    if (need_virials) then
      if (ti_mode .eq. 0) then
        vir%molecular(1,:) = vir%molecular(1,:) - &
                             molvir_netfrc_corr(1,:) * net_frcs(1)
        vir%molecular(2,:) = vir%molecular(2,:) - &
                             molvir_netfrc_corr(2,:) * net_frcs(2)
        vir%molecular(3,:) = vir%molecular(3,:) - &
                             molvir_netfrc_corr(3,:) * net_frcs(3)
      else
        vir%molecular(1,:) = vir%molecular(1,:) - &
                             ti_molvir_netfrc_corr(1,1,:) * ti_net_frcs(1,1) - &
                             ti_molvir_netfrc_corr(2,1,:) * ti_net_frcs(2,1)
        vir%molecular(2,:) = vir%molecular(2,:) - &
                             ti_molvir_netfrc_corr(1,2,:) * ti_net_frcs(1,2) - &
                             ti_molvir_netfrc_corr(2,2,:) * ti_net_frcs(2,2)
        vir%molecular(3,:) = vir%molecular(3,:) - &
                             ti_molvir_netfrc_corr(1,3,:) * ti_net_frcs(1,3) - &
                             ti_molvir_netfrc_corr(2,3,:) * ti_net_frcs(2,3)
      end if
    end if

  else

    if (.not. usemidpoint) then
      do atm_lst_idx = 1, my_atm_cnt
        i = my_atm_lst(atm_lst_idx)
        frc(1, i) = frc(1, i) + nb_frc(1, i)
        frc(2, i) = frc(2, i) + nb_frc(2, i)
        frc(3, i) = frc(3, i) + nb_frc(3, i)
      end do
    end if

  end if

  ! Calculate vir_vs_ene in master.

  if (master) then

    if (ti_mode .eq. 0) then
      vir_vs_ene = vir%elec_recip(1, 1) + &
                   vir%elec_recip(2, 2) + &
                   vir%elec_recip(3, 3) + &
                   vir%eedvir + &
                   vir%elec_nb_adjust(1, 1) + &
                   vir%elec_nb_adjust(2, 2) + &
                   vir%elec_nb_adjust(3, 3)

      ! Avoid divide-by-zero for pure neutral systems (l-j spheres), or if
      ! energies are not being calculated this step...

      if (pot_ene%elec_tot .ne. 0.0d0) then
        vir_vs_ene = abs(vir_vs_ene + pot_ene%elec_tot)/abs(pot_ene%elec_tot)
      else
        vir_vs_ene = 0.0d0
      end if
    else
      ! Avoid divide-by-zero for pure neutral systems (l-j spheres), or if
      ! energies are not being calculated this step...
      vir_vs_ene = 0.d0
      if (ti_vve(ti_ene0) .ne. 0.0d0) then
        vir_vs_ene = abs(ti_vve(ti_vir0) + &
          ti_vve(ti_ene0) * ti_weights(1))/abs(ti_vve(ti_ene0))
      end if
      if (ti_vve(ti_ene1) .ne. 0.d0) then
        vir_vs_ene = vir_vs_ene + abs(ti_vve(ti_vir1) + &
          ti_vve(ti_ene1) * ti_weights(2))/abs(ti_vve(ti_ene1))
      end if
    end if
    pme_err_est = vir_vs_ene
  end if

  ! Save virials in form used in runmd:

  if (need_virials) then
    virial(1) = 0.5d0 * vir%molecular(1, 1)
    virial(2) = 0.5d0 * vir%molecular(2, 2)
    virial(3) = 0.5d0 * vir%molecular(3, 3)
  end if

  ! Calculate the NMR restraint energy contributions, if requested.
  if (.not. usemidpoint) then
    if (nmropt .ne. 0) then
      call nmr_calc(crd, frc, enmr, 6)
      write(mdout,*) '3280 Adding to restraint_ene', pot_ene%restraint, enmr(1)
      pot_ene%restraint = pot_ene%restraint + enmr(1) + enmr(2) + enmr(3)
      pot_ene%total = pot_ene%total + enmr(1) + enmr(2) + enmr(3)
    end if

    ! Built-in X-ray target function and gradient
#ifndef MPI
#ifndef NOXRAY
   call xray_get_derivative(crd,frc,nstep,xray_e)
       pot_ene%restraint = pot_ene%restraint + xray_e
       pot_ene%total = pot_ene%total + xray_e
#endif
#endif
  end if

  if (ti_mode .ne. 0 .and. need_pot_enes) then
    call ti_calc_dvdl
  end if

  if (master .and. verbose .gt. 0) then
    call write_netfrc(net_frcs)
    call pme_verbose_print(pot_ene, vir, vir_vs_ene)
  end if

  !AMD DUAL BOOST CALC START
  if (iamd.gt.0) then
    !calculate totboost and apply weight to frc. frc=frc*fwgt
    if (.not. usemidpoint) then
      call calculate_amd_total_weights(atm_cnt, pot_ene%total, pot_ene%dihedral,&
        pot_ene%amd_boost, frc, crd, my_atm_lst)
    end if
    ! Update total energy
    pot_ene%total = pot_ene%total + pot_ene%amd_boost
  end if

  if (.not. usemidpoint) then
    !GaMD DUAL BOOST CALC START
    if (igamd .gt. 0) then
      !calculate totboost and apply weight to frc. frc=frc*fwgt
      if (igamd .eq. 1 .or. igamd .eq. 2 .or. igamd .eq. 3) then
        call calculate_gamd_total_weights(atm_cnt, pot_ene%total, pot_ene%dihedral, &
                                          pot_ene%gamd_boost, frc, crd, my_atm_lst)
        ! Update total energy
        pot_ene%total = pot_ene%total + pot_ene%gamd_boost
      else if (igamd .eq. 4 .or. igamd .eq. 5 ) then
        pot_ene_nb = pot_ene%elec_tot + pot_ene%vdw_tot + pot_ene%vdw_14 + pot_ene%elec_14 + pot_ene%efield
        nb_frc_tmp = nb_frc_gamd
        call calculate_gamd_nb_weights(atm_cnt,pot_ene_nb,pot_ene%dihedral, &
                                       pot_ene%gamd_boost,nb_frc_gamd,crd,my_atm_lst)
        frc = frc - nb_frc_tmp + nb_frc_gamd
      ! Update total energy
      pot_ene%total = pot_ene%total + pot_ene%gamd_boost
      else if( (igamd.eq.6 .or. igamd.eq.7) .and. (ti_mode.ne.0)) then
        pot_ene_lig = ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                      ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                      ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
        pot_ene_dih = pot_ene%total + &
                      ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene)
        call calculate_gamd_ti_region_weights(atm_cnt,pot_ene_lig, pot_ene_dih, &
                        pot_ene%gamd_boost,fwgt,frc,crd,ti_sc_lst,my_atm_lst)
        ! Update total energy
        ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
      else if( (igamd.eq.8 .or. igamd.eq.9) .and. (ti_mode.ne.0)) then
        pot_ene_lig = ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene)
        pot_ene_dih = pot_ene%total + &
                      ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene) + &
                      ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                      ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
        call calculate_gamd_ti_region_weights(atm_cnt,pot_ene_lig, pot_ene_dih, &
                        pot_ene%gamd_boost,fwgt,frc,crd,ti_sc_lst,my_atm_lst)
        ! Update total energy
        ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
      else if( (igamd.eq.10 .or. igamd.eq.11) .and. (ti_mode.ne.0)) then
        pot_ene_lig = ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                      ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
        pot_ene_dih = pot_ene%total + &
                      ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                      ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene)
        call calculate_gamd_ti_region_weights(atm_cnt,pot_ene_lig, pot_ene_dih, &
                        pot_ene%gamd_boost,fwgt,frc,crd,ti_sc_lst,my_atm_lst)
        ! Update total energy
        ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
      end if
    end if

    !scaledMD scale forces
    if (scaledMD.gt.0) then
      !apply scale to frc. frc=frc*scaledMD_lambda
      call scaledMD_scale_frc(atm_cnt, pot_ene%total, &
                              frc, crd, my_atm_lst)
      pot_ene%total = pot_ene%total * scaledMD_lambda
    end if

  end if

  ! If belly is on then set the belly atom forces to zero:
  if (ibelly .gt. 0) call bellyf(atm_cnt, atm_igroup, frc)

  call update_time(nonbond_time)
  call update_pme_time(pme_misc_timer)

  if (.not. usemidpoint) then
    if (charmm_active .and. do_charmm_dump_gold .eq. 1) then
      call mpi_gathervec(atm_cnt, frc)
      if (master) then
        call pme_charmm_dump_gold(atm_cnt,frc,pot_ene)
        write(mdout, '(a)') 'charmm_gold() completed. Exiting'
      end if
      call mexit(6, 0)
    end if

    if (igamd.eq.4 .or. igamd.eq.5) then
      deallocate(nb_frc_gamd,nb_frc_tmp)
    end if

    if (infe.gt.0) then
        bias_frc(:,:) = 0.d0
        call nfe_on_force(crd,bias_frc,pot_ene%nfe)

        do atm_lst_idx = 1, my_atm_cnt
          i = my_atm_lst(atm_lst_idx)
          frc(1,i)=frc(1,i)+bias_frc(1,i)
          frc(2,i)=frc(2,i)+bias_frc(2,i)
          frc(3,i)=frc(3,i)+bias_frc(3,i)
        end do

        pot_ene%total = pot_ene%total + pot_ene%nfe
        pot_ene%restraint = pot_ene%restraint + pot_ene%nfe
    end if

    if (ineb .gt. 0) then
      call mpi_gathervec(atm_cnt, frc)
      call mpi_gathervec(atm_cnt, crd)

      if (mytaskid .eq. 0) then
        call full_neb_forces(atm_mass, crd, frc, pot_ene%total, fitmask, rmsmask, nstep)
      end if
      call mpi_bcast(neb_force, 3*nattgtrms, MPI_DOUBLE_PRECISION, 0, pmemd_comm, err_code_mpi)

      do jatm=1,nattgtrms
        iatm=rmsmask(jatm)
        if (iatm .eq. 0) cycle
        if (gbl_atm_owner_map(iatm) .eq. mytaskid) then
          frc(1,iatm)= frc(1,iatm)+neb_force(1,jatm)
          frc(2,iatm)= frc(2,iatm)+neb_force(2,jatm)
          frc(3,iatm)= frc(3,iatm)+neb_force(3,jatm)
        end if
      enddo

      if (beadid .eq. 1 .or. beadid .eq. neb_nbead) then
        do i=1, my_atm_cnt
          iatm = my_atm_lst(i)
          frc(:, iatm) = 0.d0
        enddo
      end if

    end if

  end if ! not usemidpoint
#endif /* CUDA */

  return

end subroutine pme_force_midpoint

#else /* uniprocessor version below: */

subroutine pme_force(atm_cnt, crd, frc, img_atm_map, atm_img_map, &
                     new_list, need_pot_enes, need_virials, &
                     pot_ene, nstep, virial, ekcmt, pme_err_est)

  use angles_mod
  use bonds_mod
  use constraints_mod
  use dihedrals_mod
  use dynamics_dat_mod
  use dynamics_mod
  use emap_mod, only : emapforce
  use ene_frc_splines_mod
  use nb_exclusions_mod
  use pme_direct_mod
  use pme_slab_recip_mod
  use loadbal_mod
  use img_mod
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use nb_pairlist_mod
  use nmr_calls_mod
#ifndef NOXRAY
  use xray_interface_module, only: xray_get_derivative
#endif
  use parallel_dat_mod
  use parallel_mod
  use pbc_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
  use runfiles_mod
  use timers_mod
  use charmm_mod, only : charmm_active
  use mdin_debugf_dat_mod, only : do_charmm_dump_gold
  use nbips_mod
  use amd_mod
  use gamd_mod
  use scaledMD_mod
  use ti_mod
#ifdef TIDECOMP
  use ti_decomp_mod
#endif
  use nfe_setup_mod, only : nfe_on_force => on_force
  use nfe_lib_mod
#if defined(CUDA) && !defined(GTI)
    use energy_records_mod, only : afe_gpu_sc_ene_rec, null_afe_gpu_sc_ene_rec
#endif

  implicit none

  ! Formal arguments:

  integer                       :: atm_cnt
  double precision              :: crd(3, atm_cnt)
  double precision              :: frc(3, atm_cnt)
  integer                       :: img_atm_map(atm_cnt)
  integer                       :: atm_img_map(atm_cnt)
  logical                       :: new_list
  logical                       :: need_pot_enes
  logical                       :: need_virials
  type(pme_pot_ene_rec)         :: pot_ene
  integer, intent(in)           :: nstep
  double precision, optional    :: virial(3)            ! Only used for MD
  double precision, optional    :: ekcmt(3)             ! Only used for MD
  double precision, optional    :: pme_err_est          ! Only used for MD

  ! Local variables:

  type(pme_virial_rec)          :: vir
  double precision              :: enmr(3), entr
#ifndef MPI
#ifndef NOXRAY
  double precision              :: xray_e, residual, Fcalc_scale, norm_scale
#endif
#endif
  double precision              :: enemap
  double precision              :: vir_vs_ene
  integer                       :: alloc_failed
  integer                       :: i, j
  logical                       :: params_may_change
  logical                       :: onstep
  double precision              :: net_frcs(3)
  double precision, allocatable :: img_frc(:,:)
  double precision              :: evdwex, eelex
  double precision              :: bias_frc(3,atm_cnt)
#if defined(CUDA) && !defined(GTI)
    type(afe_gpu_sc_ene_rec)      :: sc_pot_ene
#endif

  ! GaMD Local Variables
  logical,save                  :: need_virials_bk
  type(pme_pot_ene_rec)         :: pot_ene2
  double precision              :: pot_ene_nb, pot_ene_lig, pot_ene_dih, fwgt, dih_boost
  double precision, allocatable :: nb_frc_gamd(:,:)
  double precision, allocatable :: nb_frc_tmp(:,:)
  double precision              :: net_frc_ti(3)
  double precision              :: enmr2(3)
  double precision              :: virial2(3)            ! Only used for MD
  double precision              :: ekcmt2(3)             ! Only used for MD

  double precision              :: temp_holding

  double precision              :: frc2(3,atm_cnt)
  integer ,save                 ::ti_mode_bk

  call zero_time()
  call zero_pme_time()
  onstep = mod(irespa, nrespa) .eq. 0

  need_virials_bk=need_virials

  if (ti_mode .ne. 0) then
    call ti_zero_arrays

#if defined(CUDA) && !defined(GTI)
      if (ti_mode .gt. 1 .and. need_pot_enes) then
        sc_pot_ene = null_afe_gpu_sc_ene_rec
      end if
#endif
    end if

    pot_ene = null_pme_pot_ene_rec
    vir = null_pme_virial_rec
    virial(:) = 0.d0
    ekcmt(:) = 0.d0
    pme_err_est = 0.d0

    ! Zero internal energies, virials, etc.

    enmr(:) = 0.d0
    vir_vs_ene = 0.d0
    net_frcs(:) = 0.d0
    molvir_netfrc_corr(:,:) = 0.d0

    if (nmropt .ne. 0) then
      call nmr_weight(atm_cnt, crd, mdout)
      pot_ene = null_pme_pot_ene_rec
      enmr(:) = 0.d0
    end if
    params_may_change = .false.

    params_may_change = .false.

    if (infe .ne. 0) params_may_change = .true.

    call update_time(pme_misc_timer)

#ifdef CUDA
    
    if (use_pme .ne. 0) then
      call self(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
        params_may_change)
      if (vdwmeth .eq. 1) then
        call vdw_correction(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
          params_may_change)
      end if
      if (need_pot_enes) then ! DO NOT TOUCH THIS, SEE RUNMD.F90 TO SET need_pot_enes
        if (iamd .eq. 2 .or. iamd .eq. 3) then
          call gpu_calculate_amd_dihedral_energy_weight()
        end if
        if (igamd.eq.2 .or. igamd.eq.3 .or. igamd.eq.5) then
          call gpu_calculate_gamd_dihedral_energy_weight()
        end if
      end if

#ifdef GTI
      if (igamd.gt.0 .or. ntcmd.gt.0 .or. nteb.gt.0) then
          if(igamd.eq.13.or.igamd.eq.12.or.igamd.eq.19.or.igamd0.eq.13.or.igamd0.eq.12.or.igamd0.eq.19)then
          ti_mode=2
          ! need_virials=.false.
          call gti_build_nl_list_gamd(ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_clear_gamd(need_pot_enes, ti_mode, crd) ! # _gamd
          call gti_nb_ppi_gamd(need_pot_enes, .false., ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_update_ti_pot_energy_from_gpu(pot_ene2, enmr, ti_mode, gti_cpu_output)
          ppi_inter_ene = ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene)
          ti_mode=0
          !need_virials=need_virials_bk
          call self_gamd(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
                  params_may_change)
          if (vdwmeth .eq. 1) then
           call vdw_correction_gamd(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                                  params_may_change)
          end if
          if (igamd.eq.19) then
             call gpu_calculate_gamd_dihedral_energy_weight()
          end if
          elseif(igamd.eq.23.or.igamd0.eq.23)then
          ti_mode=2
          ! need_virials=.false.
          call gti_build_nl_list_gamd(ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_clear_gamd(need_pot_enes, ti_mode, crd) ! # _gamd
          call gti_nb_ppi_gamd(need_pot_enes, .false., ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_update_ti_pot_energy_from_gpu(pot_ene2, enmr, ti_mode, gti_cpu_output)
          ppi_inter_ene = ti_ene(1,si_elect_ene)
          ppi_bond_ene = ti_ene(1,si_vdw_ene)
          ti_mode=0
          !need_virials=need_virials_bk
          call self_gamd(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
                                    params_may_change)
          if (vdwmeth .eq. 1) then
                  call vdw_correction_gamd(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                                           params_may_change)
          endif
          elseif(igamd.eq.14.or.igamd0.eq.14.or.igamd.eq.15.or.igamd0.eq.15.or.igamd.eq.24.or.&
                 igamd0.eq.24.or.igamd.eq.25.or.igamd0.eq.25.or.igamd.eq.27.or.igamd0.eq.27)then
          ti_mode=2
          ! need_virials=.false.
          call gti_build_nl_list_gamd(ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_clear_gamd(need_pot_enes, ti_mode, crd) ! # _gamd
          call gti_nb_ppi3_gamd(need_pot_enes, .false., ti_mode, lj1264, plj1264) ! C4PairwiseCUDA #ppi contained 14NB
          call gti_bonded_ppi(need_pot_enes, .false., ti_mode)
          call gti_update_ti_pot_energy_from_gpu(pot_ene2, enmr, ti_mode, gti_cpu_output)
          if(igamd.eq.14.or.igamd0.eq.14.or.igamd.eq.15.or.igamd0.eq.15)then
           ppi_inter_ene = ti_ene(1,si_bond_ene)+ ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                          ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene) + ti_ene(1,si_vdw_ene) + &
                          ti_ene(1,si_elect_ene)
          endif
          if(igamd.eq.24.or.igamd0.eq.24)then
             ppi_inter_ene = ti_ene(1,si_bond_ene)+ ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                             ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene) + ti_ene(1,si_elect_ene)
             ppi_bond_ene  = ti_ene(1,si_vdw_ene)
          endif
          if(igamd.eq.25.or.igamd0.eq.25)then
             ppi_inter_ene = ti_ene(1,si_bond_ene)+ ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                             ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene) + ti_ene(1,si_vdw_ene)
             ppi_bond_ene  =ti_ene(1,si_elect_ene)
          endif
          if(igamd.eq.27.or.igamd0.eq.27)then
             ppi_inter_ene = ti_ene(1,si_bond_ene)+ ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                             ti_ene(1,si_vdw_ene) +ti_ene(1,si_elect_ene)
             ppi_bond_ene  = ti_ene(1,si_bond_ene)+ ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene)
           endif
          ti_mode=0
          !need_virials=need_virials_bk
          call self_gamd(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
                  params_may_change)
          if (vdwmeth .eq. 1) then
           call vdw_correction_gamd(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                                  params_may_change)
          end if
          elseif(igamd.eq.16.or.igamd.eq.17.or.igamd0.eq.16.or.igamd0.eq.17.or.igamd.eq.20.or.igamd0.eq.20&
                    .or.igamd.eq.26.or.igamd0.eq.26.or.igamd.eq.28.or.igamd0.eq.28)then
          ti_mode=2
          ! need_virials=.false.
          call gti_build_nl_list_gamd(ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_clear_gamd(need_pot_enes, ti_mode, crd) ! # _gamd
          call gti_nb_ppi2_gamd(need_pot_enes, .false.,ti_mode,lj1264,plj1264,bgpro2atm,edpro2atm) ! C4PairwiseCUDA
          call gti_update_ti_pot_energy_from_gpu(pot_ene2, enmr, ti_mode, gti_cpu_output)
          ppi_inter_ene = ti_ene(1,si_elect_ene) + ti_ene(1,si_vdw_ene)
          if(igamd.eq.20.or.igamd0.eq.20)then
            ppi_inter_ene = ti_ene(1,si_elect_ene)
            ppi_bond_ene  = ti_ene(1,si_vdw_ene)
          elseif(igamd.eq.26.or.igamd0.eq.26.or.igamd.eq.28.or.igamd0.eq.28)then
            ppi_inter_ene = ti_ene(1,si_elect_ene)+ti_ene(1,si_vdw_ene)
          endif

          if(abs(ti_ene(1,si_elect_ene)).lt.0.001.and.abs(ti_ene(1,si_vdw_ene)).le.0.001)then
            ppi_inter_ene=0.0  !!not to scale zero non-bond interaction
          endif
          ti_mode=0
          !need_virials=need_virials_bk
          call self_gamd(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
                  params_may_change)
          if (vdwmeth .eq. 1) then
           call vdw_correction_gamd(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                                  params_may_change)
          end if
          elseif(igamd.eq.18.or.igamd0.eq.18)then
          ti_mode=2
          ! need_virials=.false.
          call gti_build_nl_list_gamd(ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_clear_gamd(need_pot_enes, ti_mode, crd) ! # _gamd
          call gti_nb_ppi_gamd(need_pot_enes, .false., ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_bonded_ppi(need_pot_enes, .false., ti_mode)
          call gti_update_ti_pot_energy_from_gpu(pot_ene2, enmr, ti_mode, gti_cpu_output)
          ppi_inter_ene= ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene)+ti_ene(1,si_dihedral_ene)
!          ppi_dihedral_ene=ti_ene(1,si_dihedral_ene)
          ti_mode=0
          !need_virials=need_virials_bk
          call self_gamd(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
                  params_may_change)
          if (vdwmeth .eq. 1) then
           call vdw_correction_gamd(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                                  params_may_change)
          end if
          elseif(igamd.eq.21.or.igamd0.eq.21)then
          ti_mode=2
          ! need_virials=.false.
          call gti_build_nl_list_gamd(ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_clear_gamd(need_pot_enes, ti_mode, crd) ! # _gamd
          call gti_nb_ppi_gamd(need_pot_enes, .false., ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
!          call gti_bonded_ppi(need_pot_enes, need_virials, ti_mode)
          call gti_update_ti_pot_energy_from_gpu(pot_ene2, enmr, ti_mode, gti_cpu_output)
          ppi_inter_ene= ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene)
!          ppi_dihedral_ene=ti_ene(1,si_bond_ene)+ti_ene(1,si_angle_ene)+ti_ene(1,si_dihedral_ene)
          ti_mode=0
          !need_virials=need_virials_bk
          call self_gamd(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
                  params_may_change)
          if (vdwmeth .eq. 1) then
           call vdw_correction_gamd(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                                  params_may_change)
          end if
          elseif(igamd.eq.22.or.igamd0.eq.22)then
           ti_mode=2
           ! need_virials=.false.
           call gti_build_nl_list_gamd(ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
           call gti_clear_gamd(need_pot_enes, ti_mode, crd) ! # _gamd
           call gti_nb_ppi3_gamd(need_pot_enes, .false., ti_mode, lj1264, plj1264) ! C4PairwiseCUDA #ppi contained 14NB
           call gti_bonded_ppi(need_pot_enes, .false., ti_mode)
           call gti_update_ti_pot_energy_from_gpu(pot_ene2,enmr, ti_mode, gti_cpu_output)
           ppi_bond_ene = ti_ene(1,si_bond_ene)+ti_ene(1,si_angle_ene)+ti_ene(1,si_dihedral_ene)+ &
                                                ti_ene(1,si_elect_14_ene)+ti_ene(1,si_vdw_14_ene) !!include 14NB
!           ppi_bond_ene = ti_ene(1,si_bond_ene)+ti_ene(1,si_angle_ene)+ti_ene(1,si_dihedral_ene) !!excluded 14NB
           ppi_inter_ene= ti_ene(1,si_vdw_ene)+ti_ene(1,si_elect_ene)
           ti_mode=0
           !need_virials=need_virials_bk
           call self_gamd(pot_ene%elec_self,ew_coeff,uc_volume,vir%elec_recip_self,params_may_change)
           if(vdwmeth .eq. 1) then
           call vdw_correction_gamd(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, params_may_change)
           endif
          elseif(igamd.eq.111.or.igamd0.eq.111)then
          ti_mode=2
!          need_virials=.false.
          call gti_build_nl_list_gamd(ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_clear_gamd(need_pot_enes, ti_mode, crd) ! # _gamd
          call gti_nb_ppi_gamd(need_pot_enes, .false., ti_mode, lj1264, plj1264) ! C4PairwiseCUDA #ppi contained 14NB
          call gti_bonded_ppi(need_pot_enes, .false., ti_mode)
          call gti_update_ti_pot_energy_from_gpu(pot_ene2, enmr, ti_mode, gti_cpu_output)
          ppi_inter_ene = ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) +ti_ene(1,si_dihedral_ene) + ti_ene(1,si_vdw_ene) + &
                          ti_ene(1,si_elect_ene)
          ti_mode=0
          !need_virials=need_virials_bk
          call self_gamd(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
                  params_may_change)
          if (vdwmeth .eq. 1) then
           call vdw_correction_gamd(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                                  params_may_change)
          endif
          elseif(igamd.eq.113.or.igamd0.eq.113)then
          ti_mode=2
          need_virials=.false.
          call gti_build_nl_list_gamd(ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_clear_gamd(need_pot_enes, ti_mode, crd) ! # _gamd
          call gti_nb_ppi_gamd(need_pot_enes, need_virials, ti_mode, lj1264, plj1264) ! C4PairwiseCUDA #ppi contained 14NB
          call gti_bonded_ppi(need_pot_enes, need_virials, ti_mode)
          call gti_update_ti_pot_energy_from_gpu(pot_ene2, enmr, ti_mode, gti_cpu_output)
          ppi_inter_ene = ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene)
          ti_mode=0
          need_virials=need_virials_bk
          call self_gamd(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
                  params_may_change)
          if (vdwmeth .eq. 1) then
           call vdw_correction_gamd(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                                  params_may_change)
          endif
          elseif(igamd.eq.115.or.igamd0.eq.115)then
          ti_mode=2
          need_virials=.false.
          call gti_build_nl_list_gamd(ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_clear_gamd(need_pot_enes, ti_mode, crd) ! # _gamd
          call gti_nb_ppi_gamd(need_pot_enes, need_virials, ti_mode, lj1264, plj1264) ! C4PairwiseCUDA #ppi contained 14NB
          call gti_bonded_ppi(need_pot_enes, need_virials, ti_mode)
          call gti_update_ti_pot_energy_from_gpu(pot_ene2, enmr, ti_mode, gti_cpu_output)
          ppi_inter_ene = ti_ene(1,si_dihedral_ene) + ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene)
          ti_mode=0
          need_virials=need_virials_bk
          call self_gamd(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
                  params_may_change)
          if (vdwmeth .eq. 1) then
           call vdw_correction_gamd(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                                  params_may_change)
          endif
          elseif(igamd.eq.117.or.igamd0.eq.117)then
          ti_mode=2
          need_virials=.false.
          call gti_build_nl_list_gamd(ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_clear_gamd(need_pot_enes, ti_mode, crd) ! # _gamd
          call gti_bonded_ppi(need_pot_enes, need_virials, ti_mode)
          call gti_update_ti_pot_energy_from_gpu(pot_ene2, enmr, ti_mode, gti_cpu_output)
          ppi_inter_ene = ti_ene(1,si_dihedral_ene)
          ti_mode=0
          need_virials=need_virials_bk
          call self_gamd(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
                  params_may_change)
          if (vdwmeth .eq. 1) then
           call vdw_correction_gamd(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                                  params_may_change)
          endif
          elseif(igamd.eq.118.or.igamd0.eq.118)then
          ti_mode=2
          need_virials=.false.
          call gti_build_nl_list_gamd(ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_clear_gamd(need_pot_enes, ti_mode, crd) ! # _gamd
          call gti_bonded_ppi(need_pot_enes, need_virials, ti_mode)
          call gti_update_ti_pot_energy_from_gpu(pot_ene2, enmr, ti_mode, gti_cpu_output)
          ppi_inter_ene = ti_ene(1,si_bond_ene)
          ti_mode=0
          need_virials=need_virials_bk
          call self_gamd(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
                  params_may_change)
          if (vdwmeth .eq. 1) then
           call vdw_correction_gamd(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                                  params_may_change)
          endif
          elseif(igamd.eq.119.or.igamd0.eq.119)then
          ti_mode=2
          need_virials=.false.
          call gti_build_nl_list_gamd(ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_clear_gamd(need_pot_enes, ti_mode, crd) ! # _gamd
          call gti_bonded_ppi(need_pot_enes, need_virials, ti_mode)
          call gti_update_ti_pot_energy_from_gpu(pot_ene2, enmr, ti_mode, gti_cpu_output)
          ppi_inter_ene = ti_ene(1,si_angle_ene)
          ti_mode=0
          need_virials=need_virials_bk
          call self_gamd(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
                  params_may_change)
          if (vdwmeth .eq. 1) then
           call vdw_correction_gamd(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                                  params_may_change)
          endif
          elseif(igamd.eq.120.or.igamd0.eq.120)then
          ti_mode=2
          need_virials=.false.
          call gti_build_nl_list_gamd(ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_clear_gamd(need_pot_enes, ti_mode, crd) ! # _gamd
          call gti_bonded_ppi(need_pot_enes, need_virials, ti_mode)
          call gti_update_ti_pot_energy_from_gpu(pot_ene2, enmr, ti_mode, gti_cpu_output)
          ppi_inter_ene = ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene)
          ti_mode=0
          need_virials=need_virials_bk
          call self_gamd(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
                  params_may_change)
          if (vdwmeth .eq. 1) then
           call vdw_correction_gamd(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                                  params_may_change)
          endif
          elseif(igamd.eq.110.or.igamd0.eq.110)then
          ti_mode=2
          need_virials=.false.
          call gti_build_nl_list_gamd(ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_clear_gamd(need_pot_enes, ti_mode, crd) ! # _gamd
          call gti_nb_ppi_gamd(need_pot_enes, need_virials, ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_update_ti_pot_energy_from_gpu(pot_ene2, enmr, ti_mode, gti_cpu_output)
          ppi_inter_ene = ti_ene(1,si_elect_ene)
          ti_mode=0
          need_virials=need_virials_bk
          call self_gamd(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
                  params_may_change)
          if (vdwmeth .eq. 1) then
           call vdw_correction_gamd(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                                  params_may_change)
          end if
          elseif(igamd.eq.112.or.igamd0.eq.112)then
          ti_mode=2
          need_virials=.false.
          call gti_build_nl_list_gamd(ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_clear_gamd(need_pot_enes, ti_mode, crd) ! # _gamd
          call gti_nb_ppi_gamd(need_pot_enes, need_virials, ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_update_ti_pot_energy_from_gpu(pot_ene2, enmr, ti_mode, gti_cpu_output)
          ppi_inter_ene = ti_ene(1,si_elect_ene)
          ti_mode=0
          need_virials=need_virials_bk
          call self_gamd(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
                  params_may_change)
          if (vdwmeth .eq. 1) then
           call vdw_correction_gamd(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                                  params_may_change)
          end if
          elseif(igamd.eq.114.or.igamd0.eq.114)then
          ti_mode=2
          need_virials=.false.
          call gti_build_nl_list_gamd(ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_clear_gamd(need_pot_enes, ti_mode, crd) ! # _gamd
          call gti_nb_ppi_gamd(need_pot_enes, need_virials, ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_update_ti_pot_energy_from_gpu(pot_ene2, enmr, ti_mode, gti_cpu_output)
          ppi_inter_ene = ti_ene(1,si_vdw_ene)
          ti_mode=0
          need_virials=need_virials_bk
          call self_gamd(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
                  params_may_change)
          if (vdwmeth .eq. 1) then
           call vdw_correction_gamd(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                                  params_may_change)
          end if
          elseif(igamd.eq.116.or.igamd0.eq.116)then
          ti_mode=2
          need_virials=.false.
          call gti_build_nl_list_gamd(ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_clear_gamd(need_pot_enes, ti_mode, crd) ! # _gamd
          call gti_nb_ppi_gamd(need_pot_enes, need_virials, ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
          call gti_update_ti_pot_energy_from_gpu(pot_ene2, enmr, ti_mode, gti_cpu_output)
          ppi_inter_ene = ti_ene(1,si_vdw_ene)
          ti_mode=0
          need_virials=need_virials_bk
          call self_gamd(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
                  params_may_change)
          if (vdwmeth .eq. 1) then
           call vdw_correction_gamd(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                                  params_may_change)
          end if

          endif
        if(.not.(igamd.eq.12.or.igamd.eq.13.or.igamd0.eq.12.or.igamd0.eq.13 .or.&
                (igamd.ge.110 .and. igamd.le.120).or.(igamd0.ge.110 .and. igamd0.le.120).or. &
                (igamd.eq.14 .or. igamd.eq.15).or.(igamd0.eq.14 .or. igamd0.eq.15).or.&
                (igamd.eq.16 .or. igamd.eq.17).or.(igamd0.eq.16 .or. igamd0.eq.17).or.&
                (igamd.eq.18 .or. igamd0.eq.18).or.(igamd.eq.19 .or. igamd0.eq.19).or.&
                (igamd.eq.20 .or. igamd0.eq.20).or.(igamd.eq.21 .or. igamd0.eq.21).or.&
                (igamd.eq.22 .or. igamd0.eq.22).or.(igamd.eq.23 .or. igamd0.eq.23).or.&
                (igamd.eq.24 .or. igamd0.eq.24).or.(igamd.eq.25 .or. igamd0.eq.25).or.&
                (igamd.eq.26 .or. igamd0.eq.26).or.(igamd.eq.27 .or. igamd0.eq.27).or.&
                (igamd.eq.28 .or. igamd0.eq.28)))then
           call gti_build_nl_list_gamd(ti_mode, lj1264, plj1264) ! C4PairwiseCUDA
           call gti_clear_gamd(need_pot_enes, ti_mode, crd) ! # _gamd
           call gti_ele_recip(need_pot_enes, need_virials, uc_volume, ti_mode, reaf_mode)
           call gti_nb_gamd(need_pot_enes, need_virials, ti_mode, lj1264, plj1264)    ! C4PairwiseCUDA # _gamd
           call gti_bonded(need_pot_enes, need_virials, ti_mode, reaf_mode)
           call gti_others(need_pot_enes, nstep, dt, uc_volume, ew_coeff)
           call gti_finalize_force_virial_gamd(numextra, need_virials, ti_mode, netfrc, frc) ! # _gamd
           ! call gti_finalize_force_gamd(numextra, need_virials, ti_mode, netfrc, frc) ! # _gamd
       else
           call gti_clear(need_pot_enes, ti_mode, crd)
           call gti_ele_recip(need_pot_enes, need_virials, uc_volume, ti_mode, reaf_mode)
           call gti_nb(need_pot_enes, need_virials, ti_mode, lj1264, plj1264, reaf_mode) ! C4PairwiseCUDA
           call gti_bonded(need_pot_enes, need_virials, ti_mode, reaf_mode)
           call gti_others(need_pot_enes, nstep, dt, uc_volume, ew_coeff)
           call gti_finalize_force_virial(numextra, need_virials, ti_mode, netfrc, frc)
      endif

        if (need_pot_enes .or. need_virials) call gti_update_md_ene(pot_ene, enmr, virial, ekcmt)

        if (ti_mode .gt. 0) then
          if (need_pot_enes) then
            if (ifmbar .ne. 0) call gti_update_mbar_from_gpu
            call gti_update_ti_pot_energy_from_gpu_gamd(pot_ene, enmr, ti_mode, gti_cpu_output)   ! # _gamd
          end if
          if (need_virials) call gti_get_virial_gamd(ti_weights(1:2), virial, ekcmt, ti_mode)     ! # _gamd
        end if

      else ! gti with igamd .gt. 0

      !! Regualr serial version of GTI routines
      call gti_build_nl_list(ti_mode, lj1264, plj1264, reaf_mode) ! C4PairwiseCUDA
      call gti_clear(need_pot_enes, ti_mode, crd)
      call gti_ele_recip(need_pot_enes, need_virials, uc_volume, ti_mode, reaf_mode)
      call gti_nb(need_pot_enes, need_virials, ti_mode, lj1264, plj1264, reaf_mode) ! C4PairwiseCUDA
      call gti_bonded(need_pot_enes, need_virials, ti_mode, reaf_mode)
      call gti_others(need_pot_enes, nstep, dt, uc_volume, ew_coeff)
      call gti_finalize_force_virial(numextra, need_virials, ti_mode, netfrc, frc)

      if (need_pot_enes .or. need_virials) call gti_update_md_ene(pot_ene, enmr, virial, ekcmt, ineb, nebfreq, nstep)
      if (ti_mode .gt. 0) then
        if (need_pot_enes) then
          call gti_update_ti_pot_energy_from_gpu(pot_ene, enmr, ti_mode, gti_cpu_output)
          if (ifmbar .ne. 0) call gti_update_mbar_from_gpu
        end if
        if (need_virials) call gti_get_virial(ti_weights(1:2), virial, ekcmt, ti_mode)
      end if
      end if ! gti only

#else /* GTI */

      if (need_pot_enes) then
        if (ti_mode .eq. 0) then
          call gpu_pme_ene(ew_coeff, uc_volume, pot_ene, enmr, virial, &
                           ekcmt, nstep, dt, crd, frc)
        else
          call gpu_ti_pme_ene(ew_coeff, uc_volume, pot_ene, sc_pot_ene, enmr, virial, &
                              ekcmt, nstep, dt, crd, frc, bar_cont)

          !for now need to write our sc energies here because our array will empty when we
          !return to runmd
          !we've already added self terms and vdw correction to ti_ene(1,si_dvdl)
          ti_ene(1,si_dvdl) = ti_ene(1,si_dvdl) + sc_pot_ene%dvdl
          ti_ene(1,si_bond_ene) = sc_pot_ene%bond_R1
          ti_ene(2,si_bond_ene) = sc_pot_ene%bond_R2
          ti_ene(1,si_angle_ene) = sc_pot_ene%angle_R1
          ti_ene(2,si_angle_ene) = sc_pot_ene%angle_R2
          ti_ene(1,si_dihedral_ene) = sc_pot_ene%dihedral_R1
          ti_ene(2,si_dihedral_ene) = sc_pot_ene%dihedral_R2
          ti_ene(1,si_vdw_ene) = sc_pot_ene%vdw_dir_R1
          ti_ene(2,si_vdw_ene) = sc_pot_ene%vdw_dir_R2
          ti_ene(1,si_elect_ene) = sc_pot_ene%elec_dir_R1
          ti_ene(2,si_elect_ene) = sc_pot_ene%elec_dir_R2
          ti_ene(1,si_vdw_14_ene) = sc_pot_ene%vdw_14_R1
          ti_ene(2,si_vdw_14_ene) = sc_pot_ene%vdw_14_R2
          ti_ene(1,si_elect_14_ene) = sc_pot_ene%elec_14_R1
          ti_ene(2,si_elect_14_ene) = sc_pot_ene%elec_14_R2
          !these are output as well, but need to deal with them differently
          !ti_ene(1,si_vdw_der_ene) = sc_pot_ene%vdw_der_R1
          !ti_ene(2,si_vdw_der_ene) = sc_pot_ene%vdw_der_R2
          !ti_ene(1,si_elect_der_ene) = sc_pot_ene%elec_der_R1
          !ti_ene(2,si_elect_der_ene) = sc_pot_ene%elec_der_R2
        end if
      else
        ! Potential energies are not needed below
        if (ti_mode .eq. 0) then
          call gpu_pme_force(ew_coeff, uc_volume, virial, ekcmt, nstep, dt, crd, frc)
        else
          call gpu_ti_pme_force(ew_coeff, uc_volume, virial, ekcmt, nstep, dt, crd, frc)
        end if
      end if

#endif /* GTI  */

    else !(use_pme .ne. 0)
      call ipsupdate(ntb)
      if (need_pot_enes) then ! DO NOT TOUCH THIS, SEE RUNMD.F90 TO SET need_pot_enes
        if (iamd .eq. 2 .or. iamd .eq. 3) then
          call gpu_calculate_amd_dihedral_energy_weight()
        end if
        if (igamd.eq.2 .or. igamd.eq.3 .or. igamd.eq.5.or.igamd.eq.19) then
          call gpu_calculate_gamd_dihedral_energy_weight()
        end if
        call gpu_ips_ene(uc_volume, pot_ene, enmr, virial, ekcmt)
      else
        call gpu_ips_force(uc_volume, virial, ekcmt)
      end if
    end if  !(use_pme .ne. 0)

    ! Built-in X-ray target function and gradient
#ifndef NOXRAY
    call xray_get_derivative(crd, frc, nstep, xray_e)
       pot_ene%restraint = pot_ene%restraint + xray_e
       pot_ene%total = pot_ene%total + xray_e
#endif

    if (nmropt .ne. 0) then
      call nmr_calc(crd, frc, enmr, 6)
    end if
#ifndef GTI
    if (ti_mode .ne. 0) then
      !set these values after calling nmr_calc
      ti_ene_aug(1,ti_rest_dist_ene) = sc_pot_ene%sc_res_dist_R1
      ti_ene_aug(2,ti_rest_dist_ene) = sc_pot_ene%sc_res_dist_R2
      ti_ene_aug(1,ti_rest_ang_ene) = sc_pot_ene%sc_res_ang_R1
      ti_ene_aug(2,ti_rest_ang_ene) = sc_pot_ene%sc_res_ang_R2
      ti_ene_aug(1,ti_rest_tor_ene) = sc_pot_ene%sc_res_tors_R1
      ti_ene_aug(2,ti_rest_tor_ene) = sc_pot_ene%sc_res_tors_R2
    end if
#endif /* !GTI */

    if (iamd.gt.0) then
      !check if the energy here is updated
      !AMD calculate weight and scale forces
      call gpu_calculate_and_apply_amd_weights(pot_ene%total, pot_ene%dihedral,&
                                               pot_ene%amd_boost,num_amd_lag)
      ! Update total energy
      pot_ene%total = pot_ene%total + pot_ene%amd_boost
    end if
    if(igamd0.eq.13.or.igamd0.eq.15.or.igamd0.eq.17 .or.igamd0.eq.18) ppi_dihedral_ene=pot_ene%total-ppi_inter_ene
    if(igamd0.eq.21.or.igamd0.eq.26.or.igamd0.eq.28)then
       if(charmm_active)then
       ppi_bond_ene = pot_ene%bond+pot_ene%angle+pot_ene%dihedral+pot_ene%angle_ub+pot_ene%imp+pot_ene%cmap
       else
       ppi_bond_ene = pot_ene%bond+pot_ene%angle+pot_ene%dihedral
       endif
       ppi_dihedral_ene=pot_ene%total-ppi_inter_ene-ppi_bond_ene
    endif
    if(igamd0.eq.20.or.igamd0.eq.23.or.igamd0.eq.24.or.igamd0.eq.25)then
       ppi_dihedral_ene=pot_ene%total-ppi_inter_ene-ppi_bond_ene
    endif
    if(igamd0.eq.22)then
       ppi_dihedral_ene=pot_ene%total-ppi_inter_ene-ppi_bond_ene
    endif
    if(igamd0.eq.27)then
       ppi_bond_ene = pot_ene%bond+pot_ene%angle+pot_ene%dihedral-ppi_bond_ene
       ppi_dihedral_ene=pot_ene%total-ppi_inter_ene-ppi_bond_ene
    endif

    if (igamd.gt.0) then
      !GAMD calculate weight and scale forces
      if (igamd.eq.1 .or. igamd.eq.2 .or. igamd.eq.3) then
        call gpu_calculate_and_apply_gamd_weights(pot_ene%total, pot_ene%dihedral, &
        pot_ene%gamd_boost,num_gamd_lag)
        ! Update total energy
        pot_ene%total = pot_ene%total + pot_ene%gamd_boost
      else if (igamd.eq.4 .or. igamd.eq.5) then
        pot_ene_nb = pot_ene%elec_tot + pot_ene%vdw_tot + pot_ene%vdw_14 + pot_ene%elec_14 + pot_ene%efield
        call gpu_calculate_and_apply_gamd_weights_nb(pot_ene_nb, pot_ene%dihedral, &
                                                     pot_ene%gamd_boost,num_gamd_lag)
        ! Update total energy
        pot_ene%total = pot_ene%total + pot_ene%gamd_boost
      else if( (igamd.eq.6 .or. igamd.eq.7) .and. (ti_mode.ne.0)) then
        ! Calculate net force on TI region
        if (mod(nstep+1,ntave).eq.0) then
          call gpu_download_frc(frc)
          net_frc_ti    = 0.0
          do i = 1, atm_cnt
                if ( ti_sc_lst(i) .gt. 0 ) then
                        net_frc_ti(:) = net_frc_ti(:) + frc(:,i)
                        ! write(*,'(a,2I10,3F15.5)') "| GaMD pme_force before scaling) step, atom, frc = ", nstep+1, i, frc(:,i)
                end if
          end do
          ! write(*,'(a,i10,3f15.5)') "| GaMD pme_force before scaling) step, frc_lig(1:3) = ", nstep+1, net_frc_ti

        end if
        pot_ene_lig = ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                      ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                      ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
        pot_ene_dih = pot_ene%total + &
                      ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene)
        call gpu_calculate_and_apply_gamd_weights_ti_region(pot_ene_lig, pot_ene_dih, &
                                                pot_ene%gamd_boost,fwgt,num_gamd_lag)
        if (igamd.eq.7) call gpu_calculate_and_apply_gamd_weights_ti_others(pot_ene_lig, pot_ene_dih, &
                                                pot_ene%gamd_boost,num_gamd_lag)
        ! Calculate net force on TI region
        if (mod(nstep+1,ntave).eq.0) then
          call gpu_download_frc(frc)
          net_frc_ti    = 0.0
          do i = 1, atm_cnt
                if ( ti_sc_lst(i) .gt. 0 ) then
                        net_frc_ti(:) = net_frc_ti(:) + frc(:,i)
                        ! write(*,'(a,2I10,3F15.5)') "| GaMD pme_force after scaling) step, atom, frc = ", nstep+1, i, frc(:,i)
                end if
          end do
          ! write(*,'(a,i10,4f15.5)') "| GaMD pme_force after scaling) step, frc_weight, frc_lig(1:3) = ", nstep+1, fwgt, net_frc_ti
        end if
        ! Update total energy
        ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
      else if( (igamd.eq.8 .or. igamd.eq.9) .and. (ti_mode.ne.0)) then
        pot_ene_lig = ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene)
        pot_ene_dih = pot_ene%total + &
                      ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                      ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene) + &
                      ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene)
        call gpu_calculate_and_apply_gamd_weights_ti_region(pot_ene_lig, pot_ene_dih, &
                                                pot_ene%gamd_boost,fwgt,num_gamd_lag)
        if (igamd.eq.9) call gpu_calculate_and_apply_gamd_weights_ti_others(pot_ene_lig, pot_ene_dih, &
                                                pot_ene%gamd_boost,num_gamd_lag)
        ! Update total energy
        ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
      else if( (igamd.eq.10 .or. igamd.eq.11) .and. (ti_mode.ne.0)) then
        ! write(*,'(a,i10,4f15.5)') "| GaMD: step, ti_ene, ti_others_ene = ", ti_ene, ti_others_ene
        pot_ene_lig = ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                      ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
        pot_ene_dih = pot_ene%total + &
                      ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                      ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene)
        ! write(*,'(a,i10,2f15.5)') "| GaMD: step, pot_ene_lig, pot_ene_dih = ", nstep+1, pot_ene_lig, pot_ene_dih
        call gpu_calculate_and_apply_gamd_weights_ti_region(pot_ene_lig, pot_ene_dih, &
                                                pot_ene%gamd_boost,fwgt,num_gamd_lag)
        if (igamd.eq.11) call gpu_calculate_and_apply_gamd_weights_ti_others(pot_ene_lig, pot_ene_dih, &
                                                pot_ene%gamd_boost,num_gamd_lag)
        ! update virial
        ! if (need_virials) call gti_virial_gamd(numextra, need_virials, ti_mode, netfrc, frc) ! # _gamd
      ! Update total energy
        ! ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
        pot_ene%total = pot_ene%total + pot_ene%gamd_boost
        else if (igamd.eq.12 .or. igamd.eq.13 ) then
!        pot_ene_dih = pot_ene%total-pot_ene_ppi
        ppi_dihedral_ene=pot_ene%total-ppi_inter_ene
        if(igamd.eq.12)then
           call gpu_calculate_and_apply_gamd_weights_ti_region(ppi_inter_ene, ppi_dihedral_ene, &
                                                pot_ene%gamd_boost,fwgt,num_gamd_lag)
        elseif(igamd.eq.13)then
           call gpu_calculate_and_apply_gamd_weights_sc_dual_nonbonded(ppi_inter_ene, ppi_dihedral_ene, &
                                                pot_ene%gamd_boost,num_gamd_lag)
        endif

        !call gpu_download_gamd_weights(gamd_weights_and_energy)

        ! Update energy
        !ppi_inter_ene=ppi_inter_ene+gamd_weights_and_energy(5)
        !ppi_dihedral_ene=ppi_dihedral_ene+gamd_weights_and_energy(6)
        pot_ene%total=pot_ene%total+pot_ene%gamd_boost
        !print*,"pme_force ppi)dihedral",ppi_dihedral_ene
    else if (igamd.eq.14 .or. igamd.eq.15) then
        ppi_dihedral_ene = pot_ene%total-ppi_inter_ene
        if(igamd.eq.14)then
         call gpu_calculate_and_apply_gamd_weights_ti_region(ppi_inter_ene, ppi_dihedral_ene, &
                                                pot_ene%gamd_boost,fwgt,num_gamd_lag)
        elseif (igamd.eq.15) then
        call gpu_calculate_and_apply_gamd_weights_sc_dual_bonded(ppi_inter_ene, ppi_dihedral_ene, &
                                                pot_ene%gamd_boost,num_gamd_lag)
        endif
        call gpu_download_gamd_weights(gamd_weights_and_energy)

        ! Update energy
        !ppi_inter_ene=ppi_inter_ene+gamd_weights_and_energy(5)
        !ppi_dihedral_ene=ppi_dihedral_ene+gamd_weights_and_energy(6)
        pot_ene%total = pot_ene%total + pot_ene%gamd_boost

      else if (igamd.eq.16 .or. igamd.eq.17 ) then
        ppi_dihedral_ene=pot_ene%total-ppi_inter_ene
!        pot_ene_dih = pot_ene%total-pot_ene_ppi
        if(igamd.eq.16) then
        call gpu_calculate_and_apply_gamd_weights_ti_region(ppi_inter_ene, ppi_dihedral_ene, &
                                                pot_ene%gamd_boost,fwgt,num_gamd_lag)
        elseif (igamd.eq.17) then
             call gpu_calculate_and_apply_gamd_weights_sc_dual_nonbonded(ppi_inter_ene, ppi_dihedral_ene, &
                                                pot_ene%gamd_boost,num_gamd_lag)
        endif
        !call gpu_download_gamd_weights(gamd_weights_and_energy)

        ! Update energy
        !ppi_inter_ene=ppi_inter_ene+gamd_weights_and_energy(5)
        !ppi_dihedral_ene=ppi_dihedral_ene+gamd_weights_and_energy(6)
        pot_ene%total=pot_ene%total+pot_ene%gamd_boost
      else if (igamd.eq.18) then
!         pot_ene_dih = pot_ene%dihedral
!         call gpu_calculate_and_apply_gamd_weights_ti_region(ppi_inter_ene, ppi_dihedral_ene, &
!                                             pot_ene%gamd_boost,fwgt,num_gamd_lag)
!         call gpu_calculate_and_apply_gamd_weights_ti_others(ppi_inter_ene, ppi_dihedral_ene, &
!                                                pot_ene%gamd_boost,num_gamd_lag)
       call gpu_calculate_and_apply_gamd_weights_sc_dual_bonded(ppi_inter_ene, ppi_dihedral_ene, &
                                                               pot_ene%gamd_boost,num_gamd_lag)
!       call gpu_download_gamd_weights(gamd_weights_and_energy)

        ! Update energy
       ! ppi_inter_ene=ppi_inter_ene+gamd_weights_and_energy(5)
       ! ppi_dihedral_ene=ppi_dihedral_ene+gamd_weights_and_energy(6)
        pot_ene%total=pot_ene%total+pot_ene%gamd_boost

      else if (igamd.eq.19) then
         pot_ene_dih = pot_ene%dihedral
         call gpu_calculate_and_apply_gamd_weights_ti_region(ppi_inter_ene, pot_ene_dih, &
                                             pot_ene%gamd_boost,fwgt,num_gamd_lag)

        ! call gpu_download_gamd_weights(gamd_weights_and_energy)
        ! Update energy
        !ppi_inter_ene=ppi_inter_ene+gamd_weights_and_energy(5)
        !ppi_dihedral_ene=pot_ene%dihedral
        pot_ene%total=pot_ene%total+pot_ene%gamd_boost

      else if (igamd.ge.20.and.igamd.le.28) then
        call gpu_calculate_and_apply_gamd_weights_sc_triple(ppi_inter_ene, ppi_dihedral_ene, &
                                              ppi_bond_ene, pot_ene%gamd_boost,num_gamd_lag)
      ! Update energy
      pot_ene%total=pot_ene%total+pot_ene%gamd_boost
      else if (igamd.eq.110 .or. igamd.eq.111) then
        pot_ene_dih = pot_ene%dihedral
        call gpu_calculate_and_apply_gamd_weights_ti_region(ppi_inter_ene, pot_ene_dih, &
                                                pot_ene%gamd_boost,fwgt,num_gamd_lag)
        ! Update total energy
        ppi_inter_ene=ppi_inter_ene+pot_ene%gamd_boost
        pot_ene%total = pot_ene%total + pot_ene%gamd_boost
      else if (igamd.eq.112 .or. igamd.eq.113) then
        pot_ene_dih = pot_ene%dihedral
        call gpu_calculate_and_apply_gamd_weights_ti_region(ppi_inter_ene, pot_ene_dih, &
                                                pot_ene%gamd_boost,fwgt,num_gamd_lag)
        ! Update total energy
        ppi_inter_ene=ppi_inter_ene+pot_ene%gamd_boost
        pot_ene%total = pot_ene%total + pot_ene%gamd_boost
      else if (igamd.eq.114 .or. igamd.eq.115) then
        pot_ene_dih = pot_ene%dihedral
        call gpu_calculate_and_apply_gamd_weights_ti_region(ppi_inter_ene, pot_ene_dih, &
                                                pot_ene%gamd_boost,fwgt,num_gamd_lag)
        ! Update total energy
        ppi_inter_ene=ppi_inter_ene+pot_ene%gamd_boost
        pot_ene%total = pot_ene%total + pot_ene%gamd_boost
      else if (igamd.eq.116 .or. igamd.eq.117) then
        pot_ene_dih = pot_ene%dihedral
        call gpu_calculate_and_apply_gamd_weights_ti_region(ppi_inter_ene, pot_ene_dih, &
                                                pot_ene%gamd_boost,fwgt,num_gamd_lag)
        ! Update total energy
        ppi_inter_ene=ppi_inter_ene+pot_ene%gamd_boost
        pot_ene%total = pot_ene%total + pot_ene%gamd_boost
      else if (igamd.eq.118 .or. igamd.eq.119 .or. igamd.eq.120) then
        pot_ene_dih = pot_ene%dihedral
        call gpu_calculate_and_apply_gamd_weights_ti_region(ppi_inter_ene, pot_ene_dih, &
                                                pot_ene%gamd_boost,fwgt,num_gamd_lag)
        ! Update total energy
        ppi_inter_ene=ppi_inter_ene+pot_ene%gamd_boost
        pot_ene%total = pot_ene%total + pot_ene%gamd_boost

      else if( (igamd.eq.100) .and. (ti_mode.ne.0)) then
        pot_ene_lig = ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
        pot_ene_dih = pot_ene%total + &
                      ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                      ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene)
        call gpu_calculate_and_apply_gamd_weights_ti_region(pot_ene_lig, pot_ene_dih, &
                                                pot_ene%gamd_boost,fwgt,num_gamd_lag)
        call gpu_calculate_and_apply_gamd_weights_ti_others(pot_ene_lig, pot_ene_dih, &
                                                pot_ene%gamd_boost,num_gamd_lag)
        ! update virial
        ! if (need_virials) call gti_virial_gamd(numextra, need_virials, ti_mode, netfrc, frc) ! # _gamd
        ! Update total energy
        ! ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
        pot_ene%total = pot_ene%total + pot_ene%gamd_boost
      else if( (igamd.eq.101) .and. (ti_mode.ne.0)) then
        pot_ene_lig = ti_ene(1,si_elect_14_ene) + ti_others_ene(1,si_elect_14_ene)
        pot_ene_dih = pot_ene%dihedral
        call gpu_calculate_and_apply_gamd_weights_ti_region(pot_ene_lig, pot_ene_dih, &
                                                pot_ene%gamd_boost,fwgt,num_gamd_lag)
        ! update virial
        ! if (need_virials) call gti_virial_gamd(numextra, need_virials, ti_mode, netfrc, frc) ! # _gamd
        ! Update total energy
        ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
      else if( (igamd.eq.102) .and. (ti_mode.ne.0)) then
        pot_ene_lig = ti_ene(1,si_vdw_14_ene) + ti_others_ene(1,si_vdw_14_ene)
        pot_ene_dih = pot_ene%dihedral
        call gpu_calculate_and_apply_gamd_weights_ti_region(pot_ene_lig, pot_ene_dih, &
                                                pot_ene%gamd_boost,fwgt,num_gamd_lag)
        ! update virial
        ! if (need_virials) call gti_virial_gamd(numextra, need_virials, ti_mode, netfrc, frc) ! # _gamd
        ! Update total energy
        ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
      else if( (igamd.eq.103) .and. (ti_mode.ne.0)) then
        pot_ene_lig = ti_ene(1,si_elect_ene) + ti_others_ene(1,si_elect_ene)
        pot_ene_dih = pot_ene%dihedral
        call gpu_calculate_and_apply_gamd_weights_ti_region(pot_ene_lig, pot_ene_dih, &
                                                pot_ene%gamd_boost,fwgt,num_gamd_lag)
        ! update virial
        ! if (need_virials) call gti_virial_gamd(numextra, need_virials, ti_mode, netfrc, frc) ! # _gamd
        ! Update total energy
        ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
      else if( (igamd.eq.104) .and. (ti_mode.ne.0)) then
        pot_ene_lig = ti_ene(1,si_vdw_ene) + ti_others_ene(1,si_vdw_ene)
        pot_ene_dih = pot_ene%dihedral
        call gpu_calculate_and_apply_gamd_weights_ti_region(pot_ene_lig, pot_ene_dih, &
                                                pot_ene%gamd_boost,fwgt,num_gamd_lag)
        ! update virial
        ! if (need_virials) call gti_virial_gamd(numextra, need_virials, ti_mode, netfrc, frc) ! # _gamd
        ! Update total energy
        ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
      end if
    end if
    if (scaledMD.gt.0) then
      !scaledMD scale forces
      call gpu_scaledmd_scale_frc(pot_ene%total)
      pot_ene%total = pot_ene%total * scaledMD_lambda
    end if

    !same call as CPU
    if (ti_mode .ne. 0 .and. need_pot_enes .and. nstep .ge. 0) then
      call ti_calc_dvdl
    end if

    call update_time(nonbond_time)


    if (infe .gt. 0 .or. iextpot .gt. 0) then
        if (nfe_first) then
           call gpu_setup_shuttle_info(nfe_atm_cnt, 0, nfe_atm_lst)
           nfe_first = .false.
        end if

        call gpu_shuttle_retrieve_data(crd, 0)
        call gpu_shuttle_retrieve_data(frc, 1)

        if (infe .gt. 0) then
           bias_frc(:,:) = 0.d0
           call nfe_on_force(crd,bias_frc,pot_ene%nfe)
           do i=1,nfe_atm_cnt
              j=nfe_atm_lst(i)

              frc(1,j) = frc(1,j)+bias_frc(1,j)
              frc(2,j) = frc(2,j)+bias_frc(2,j)
              frc(3,j) = frc(3,j)+bias_frc(3,j)
           enddo

           pot_ene%total = pot_ene%total + pot_ene%nfe
           pot_ene%restraint = pot_ene%restraint + pot_ene%nfe
        end if

        !take advantage of the fact that the nfe code 
        !already gives us coordinate and force access,
        !use that data for any external potential modifiers,
        !even if infe and iextpot aren't supported together.
        if (iextpot .gt. 0) then
           temp_holding = 0.0
           
           call pme_external(crd, frc, temp_holding)
           pot_ene%restraint = pot_ene%restraint + temp_holding
           pot_ene%total     = pot_ene%total     + temp_holding
        endif


        !send force modifications back up.
        call gpu_shuttle_post_data(frc, 1)
        
    end if !end if (infe .or. iextpot)
#else /*  CUDA */

    ! Zero energies that are stack or call parameters:

    pot_ene = null_pme_pot_ene_rec
    vir = null_pme_virial_rec

    virial(:) = 0.d0
    ekcmt(:) = 0.d0
    pme_err_est = 0.d0

    ! Zero internal energies, virials, etc.

    enmr(:) = 0.d0
    vir_vs_ene = 0.d0

    net_frcs(:) = 0.d0
    molvir_netfrc_corr(:,:) = 0.d0

!!Avoid double-add of NMR energy    
#if 0
    if (nmropt .ne. 0) then
       call nmr_calc(crd, frc, enmr, 6)
       pot_ene%restraint = pot_ene%restraint + enmr(1) + enmr(2) + enmr(3)
       pot_ene%total = pot_ene%total + enmr(1) + enmr(2) + enmr(3)
    end if
    if (natc .gt. 0) then
        call get_crd_constraint_energy(natc, entr, atm_jrc, &
                                       crd, frc, atm_xc, atm_weight)
        pot_ene%restraint = pot_ene%restraint + entr
        pot_ene%total = pot_ene%total + entr
    end if
#endif

    ! Do weight changes, if requested.

    if (nmropt .ne. 0) call nmr_weight(atm_cnt, crd, mdout)

    ! If no force calcs are to be done, clear the frc array and bag out now.

    if (ntf .eq. 8) then
      frc(:,:) = 0.d0
      return
    end if

    allocate(img_frc(3, atm_cnt), &
      stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    if (igamd .eq. 4 .or. igamd .eq. 5) then
      allocate(nb_frc_gamd(3, atm_cnt), &
               nb_frc_tmp(3, atm_cnt), &
               stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
    end if

    ! Calculate the non-bonded contributions:

    ! Direct part of ewald plus vdw, hbond, pairlist setup and image claiming:

    params_may_change = (nmropt .ne. 0) .or. (infe .ne. 0)

    if (ntp .gt. 0) call fill_tranvec(gbl_tranvec)

    ! The following encapsulation (save_imgcrds) seems to be necessary to
    ! prevent an optimization bug with the SGI f90 compiler.  Sigh...

    if (new_list) then
      call pme_list(atm_cnt, crd, atm_nb_maskdata, atm_nb_mask)
      call save_imgcrds(atm_cnt, gbl_img_crd, gbl_saved_imgcrd)
    else
      call adjust_imgcrds(atm_cnt, gbl_img_crd, img_atm_map, &
                          gbl_saved_imgcrd, crd, gbl_atm_saved_crd, ntp)
    end if

    do i = 1, atm_cnt
      img_frc(1, i) = 0.d0
      img_frc(2, i) = 0.d0
      img_frc(3, i) = 0.d0
    end do
    if (ti_mode .ne. 0) then
      call ti_zero_arrays
      do i = 1, atm_cnt
        ti_img_frc(:, :, i) = 0.d0
        ti_nb_frc(:, :, i) = 0.d0
      end do
    end if

    evdwex = 0.0
    eelex = 0.0

    if ( ips .gt. 0 ) then
      frc(:,:) = 0.d0
      call eexips(evdwex,eelex,frc,crd, &
                  gbl_img_qterm,gbl_ipairs,atm_nb_maskdata, &
                  atm_nb_mask,img_atm_map,atm_cnt,gbl_tranvec)
    end if

    ! Don't do recip if PME is not invoked. Don't do it this step unless
    ! mod(irespa,nrepsa) = 0

    if (use_pme .ne. 0 .and. onstep) then

      ! Self energy:

      ! The small amount of time used here gets lumped with the recip stuff...

      call self(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
                params_may_change)

      ! Reciprocal energy:

      call update_pme_time(pme_misc_timer)
      ! Only the slab fft implementation is used in uniprocessor runs...
      if (ti_mode .eq. 0) then
        call do_slab_pmesh_kspace(img_frc, pot_ene%elec_recip, vir%elec_recip, &
                                  need_pot_enes, need_virials)
        if (nrespa .gt. 1) call respa_scale(atm_cnt, img_frc, nrespa)
      else
        ti_mask_piece = 2
#ifdef TIDECOMP
        decomp_region = 1
        if(ti_mode .eq. 3) decomp_region = 0
#endif
        call do_slab_pmesh_kspace(img_frc, pot_ene%elec_recip, &
                                  ti_vir(1,:,:), need_pot_enes, need_virials)
        do i = 1, atm_cnt
          ti_img_frc(1,:,i) = img_frc(:, i)
          img_frc(:, i) = 0.d0
        end do
        ti_pot(1) = pot_ene%elec_recip

        pot_ene%elec_recip = 0.d0
        ti_mask_piece = 1
        if (ti_mode .ne. 3) then
#ifdef TIDECOMP
          decomp_region =2
#endif
          call do_slab_pmesh_kspace(img_frc, pot_ene%elec_recip, &
                                    ti_vir(2,:,:), need_pot_enes, need_virials)
          do i = 1, atm_cnt
            ti_img_frc(2,:,i) = img_frc(:, i)
            img_frc(:, i) = 0.d0
          end do
          ti_pot(2) = pot_ene%elec_recip
        else
          ti_pot(2) = ti_pot(1)
          do i = 1, atm_cnt
            ti_img_frc(2, :, i) = ti_img_frc(1, :, i)
          end do
          ti_vir(2,:,:) = ti_vir(1,:,:)
        end if

        do i = 1, atm_cnt
          ti_img_frc(1, :, i) = ti_img_frc(1, :, i) * ti_item_weights(3,1)
          ti_img_frc(2, :, i) = ti_img_frc(2, :, i) * ti_item_weights(3,2)
        end do

        vir%elec_recip = ti_vir(1,:,:) * ti_weights(1) + &
          ti_vir(2,:,:) * ti_weights(2)

        ti_vve(ti_vir0) = (ti_vir(1,1,1) + &
          ti_vir(1,2,2) + &
          ti_vir(1,3,3)) * ti_weights(1)
        ti_vve(ti_vir1) = (ti_vir(2,1,1) + &
          ti_vir(2,2,2) + &
          ti_vir(2,3,3)) * ti_weights(2)

        call ti_update_ene_all(ti_pot, si_elect_ene, pot_ene%elec_recip,3)

        if (nrespa .gt. 1) then
          call ti_respa_scale(atm_cnt, ti_img_frc, nrespa)
        end if
      end if
      ! Long range dispersion contributions:

      ! Continuum method:

      if (vdwmeth .eq. 1) then
        call vdw_correction(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                            params_may_change)
      end if

    end if      ! respa

    ! Direct part of ewald plus vdw, hbond, force and energy calculations:

    call update_pme_time(pme_misc_timer)

    if (ips .ne. 0) then
      call get_nb_ips_energy(img_frc, gbl_img_crd, gbl_img_qterm,  gbl_eed_cub, &
                             gbl_ipairs, gbl_tranvec, need_pot_enes, need_virials, &
                             pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                             vir%eedvir, vir%elec_direct)

    else if (ti_mode .ne. 0) then
      call get_nb_energy_ti(img_frc, gbl_img_crd, gbl_img_qterm, gbl_eed_cub, &
                            gbl_ipairs, gbl_tranvec, need_pot_enes, need_virials, &
                            pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                            vir%eedvir, vir%elec_direct, gbl_img_atm_map)
    else if ((lj1264 .eq. 1) .and. (plj1264 .eq. 0)) then ! New2022
      call get_nb_energy_1264(img_frc, gbl_img_crd, gbl_img_qterm, gbl_eed_cub, &
                              gbl_ipairs, gbl_tranvec, need_pot_enes, need_virials, &
                              pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                              vir%eedvir, vir%elec_direct)
    else if ((lj1264 .eq. 0) .and. (plj1264 .eq. 1)) then ! New2022
      call get_nb_energy_p1264(img_frc, gbl_img_crd, gbl_img_qterm, gbl_eed_cub, &
                              gbl_ipairs, gbl_tranvec, need_pot_enes, need_virials, &
                              pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                              vir%eedvir, vir%elec_direct)
    else if ((lj1264 .eq. 1) .and. (plj1264 .eq. 1)) then ! New2022
      call get_nb_energy_1264p1264(img_frc, gbl_img_crd, gbl_img_qterm, gbl_eed_cub, &
                              gbl_ipairs, gbl_tranvec, need_pot_enes, need_virials, &
                              pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                              vir%eedvir, vir%elec_direct)
    else
      call get_nb_energy(img_frc, gbl_img_crd, gbl_img_qterm, gbl_eed_cub, &
                         gbl_ipairs, gbl_tranvec, need_pot_enes, need_virials, &
                         pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                         vir%eedvir, vir%elec_direct)
    end if

    if (efx .ne. 0 .or. efy .ne. 0 .or. efz .ne. 0) then
      call get_efield_energy(img_frc, crd, gbl_img_qterm, gbl_img_atm_map, need_pot_enes, &
                             pot_ene%efield, atm_cnt, nstep)
    end if

    call update_pme_time(dir_frc_sum_timer)

    ! Transfer image forces to the force array.


    if ( ips .gt. 0 ) then
      do i = 1, atm_cnt
        j = atm_img_map(i)
        frc(1, i) = frc(1, i) + img_frc(1, j)
        frc(2, i) = frc(2, i) + img_frc(2, j)
        frc(3, i) = frc(3, i) + img_frc(3, j)
      end do
    else
      do i = 1, atm_cnt
        j = atm_img_map(i)
        frc(1, i) = img_frc(1, j)
        frc(2, i) = img_frc(2, j)
        frc(3, i) = img_frc(3, j)
      end do
    end if

    call update_pme_time(pme_misc_timer)

    ! Calculate 1-4 electrostatic energies, forces:

    if (charmm_active) then
      if (need_virials) then
        call get_nb14_energy(atm_qterm, crd, frc, atm_iac, typ_ico, &
                             gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt, &
                             pot_ene%elec_14, pot_ene%vdw_14, &
                             vir%elec_14)
      else
        call get_nb14_energy(atm_qterm, crd, frc, atm_iac, typ_ico, &
                             gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt, &
                             pot_ene%elec_14, pot_ene%vdw_14)
      end if
    else
      if (need_virials) then
        call get_nb14_energy(atm_qterm, crd, frc, atm_iac, typ_ico, &
                             gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt, &
                             pot_ene%elec_14, pot_ene%vdw_14, &
                             vir%elec_14)
      else
        call get_nb14_energy(atm_qterm, crd, frc, atm_iac, typ_ico, &
                             gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt, &
                             pot_ene%elec_14, pot_ene%vdw_14)
      end if
    end if

    call update_pme_time(dir_frc_sum_timer)

    ! Adjust energies, forces for masked out pairs:

    if (use_pme .ne. 0) then
      call nb_adjust(atm_cnt, atm_qterm, crd, gbl_nb_adjust_pairlst, &
                     gbl_eed_cub, frc, pot_ene%elec_nb_adjust, vir%elec_nb_adjust)
    end if

    call update_pme_time(adjust_masked_timer)

    if (ti_mode .ne. 0) then
      do i = 1, atm_cnt
        j = atm_img_map(i)
        ti_nb_frc(:, 1, i) = ti_nb_frc(:, 1, i) + ti_img_frc(:, 1, j)
        ti_nb_frc(:, 2, i) = ti_nb_frc(:, 2, i) + ti_img_frc(:, 2, j)
        ti_nb_frc(:, 3, i) = ti_nb_frc(:, 3, i) + ti_img_frc(:, 3, j)
      end do
    end if
    ! If using extra points and a frame (checked internal to subroutine),
    ! transfer force and torque from the extra points to the parent atom:

    if (numextra .gt. 0 .and. frameon .ne. 0) then
      if (ti_mode .eq. 0) then
        call orient_frc(crd, frc, vir%ep_frame, ep_frames, ep_lcl_crd, gbl_frame_cnt)
      else
        ti_vir(:,:,:) = 0.d0
        do i = 1, atm_cnt
          frc(:, i) = ti_nb_frc(1,:,i)
        end do
        call orient_frc(crd, frc, ti_vir(1,:,:), ep_frames, ep_lcl_crd, gbl_frame_cnt)
        do i = 1, atm_cnt
          ti_nb_frc(1,:,i) = frc(:, i)
          frc(:, i) = ti_nb_frc(2,:,i)
        end do
        call orient_frc(crd, frc, ti_vir(2,:,:), ep_frames, ep_lcl_crd, gbl_frame_cnt)
        do i = 1, atm_cnt
          ti_nb_frc(2,:,i) = frc(:, i)
          frc(:, i) = 0.d0
        end do
        vir%ep_frame = ti_vir(1,:,:) + ti_vir(2,:,:)
      end if
    end if
    call update_pme_time(dir_frc_sum_timer)

    ! Calculate total nonbonded energy components.

    pot_ene%vdw_tot = pot_ene%vdw_dir + &
                      pot_ene%vdw_recip

    pot_ene%elec_tot = pot_ene%elec_dir + &
                       pot_ene%elec_recip + &
                       pot_ene%elec_nb_adjust + &
                       pot_ene%elec_self

    if ( ips .gt. 0 ) then
      pot_ene%vdw_tot = pot_ene%vdw_tot + evdwex
      pot_ene%elec_tot = pot_ene%elec_tot + eelex
    end if

    if (ti_mode .ne. 0 ) then
      call ti_combine_frcs(atm_cnt, frc, ti_nb_frc)
    end if

    ! If the net force correction is in use, here we determine the net forces by
    ! looking at the sum of all nonbonded forces.  This should give the same
    ! result as just looking at the reciprocal forces, but it is more
    ! computationally convenient, especially for extra points, to do it this way.

    ! NOTE: this is not currently done for the GPU

#ifndef CUDA
    if (netfrc .gt. 0 .and. onstep) then

      if (ti_mode .eq. 0) then
        do i = 1, atm_cnt
          net_frcs(:) = net_frcs(:) + frc(:, i)
        end do

        ! Now do the correction:

        net_frcs(:) = net_frcs(:) / dble(atm_cnt - numextra)

        do i = 1, atm_cnt
          frc(:, i) = frc(:, i) - net_frcs(:)
        end do
      else
        do i = 1, atm_cnt
          ti_net_frcs(1, :) = ti_net_frcs(1, :) + ti_nb_frc(1, :, i)
          ti_net_frcs(2, :) = ti_net_frcs(2, :) + ti_nb_frc(2, :, i)
        end do

        ti_net_frcs(1,:) = ti_net_frcs(1,:)/dble(ti_atm_cnt(1)-ti_numextra_pts(1))
        ti_net_frcs(2,:) = ti_net_frcs(2,:)/dble(ti_atm_cnt(2)-ti_numextra_pts(2))
        net_frcs(:) = ti_net_frcs(1,:) + ti_net_frcs(2,:)

        do i = 1, atm_cnt
          ! This matches how sander removes netfrcs in TI runs
          if (ti_lst(1,i) .ne. 0) then
            frc(:, i) = frc(:, i) - ti_net_frcs(1,:)
          else if (ti_lst(2,i) .ne. 0) then
            frc(:, i) = frc(:, i) - ti_net_frcs(2,:)
          else
            frc(:, i) = frc(:, i) - net_frcs(:)
          end if
        end do
      end if
      ! Any extra points must have their 0.d0 forces reset...

      if (numextra .gt. 0 .and. frameon .ne. 0) &
        call zero_extra_pnts_vec(frc, ep_frames, gbl_frame_cnt)

    end if
#endif

    call update_pme_time(pme_misc_timer)

    ! First phase of virial work.  We need just the nonbonded forces at this
    ! stage, though there will be a further correction in get_dihed.

    if (need_virials) then
      call get_atm_rel_crd(my_mol_cnt, gbl_mol_com, crd, atm_rel_crd)
      if (ti_mode .eq. 0) then
        do i = 1, atm_cnt
          vir%molecular(:,1) = vir%molecular(:,1) + frc(:,i) * atm_rel_crd(1,i)
          vir%molecular(:,2) = vir%molecular(:,2) + frc(:,i) * atm_rel_crd(2,i)
          vir%molecular(:,3) = vir%molecular(:,3) + frc(:,i) * atm_rel_crd(3,i)
        end do
      else
        do i = 1, atm_cnt
          vir%molecular(:,1) = vir%molecular(:,1) + &
                               (ti_nb_frc(1,:,i) - ti_net_frcs(1,:)) * atm_rel_crd_sc(1,1,i) + &
                               (ti_nb_frc(2,:,i) - ti_net_frcs(2,:)) * atm_rel_crd_sc(2,1,i)
          vir%molecular(:,2) = vir%molecular(:,2) + &
                               (ti_nb_frc(1,:,i) - ti_net_frcs(1,:)) * atm_rel_crd_sc(1,2,i) + &
                               (ti_nb_frc(2,:,i) - ti_net_frcs(2,:)) * atm_rel_crd_sc(2,2,i)
          vir%molecular(:,3) = vir%molecular(:,3) + &
                               (ti_nb_frc(1,:,i) - ti_net_frcs(1,:)) * atm_rel_crd_sc(1,3,i) + &
                               (ti_nb_frc(2,:,i) - ti_net_frcs(2,:)) * atm_rel_crd_sc(2,3,i)
        end do
      end if
      call get_ekcom(my_mol_cnt, gbl_mol_mass_inv, ekcmt, atm_vel, atm_mass)
    end if

    call update_time(nonbond_time)
    call update_pme_time(pme_misc_timer)

    ! GaMD
    if (igamd.eq.4 .or. igamd.eq.5) then
      nb_frc_gamd = frc
    end if

    ! Calculate the other contributions:

    call pme_bonded_force(crd, frc, pot_ene)

    ! --- calculate the EMAP constraint energy ---

    if (iemap .gt. 0) then
      call emapforce(natom, enemap, atm_mass, crd, frc)
      pot_ene%emap = enemap
      pot_ene%restraint = pot_ene%restraint + enemap
    end if

    ! Sum up total potential energy for this task:

    pot_ene%total = pot_ene%vdw_tot + &
                    pot_ene%elec_tot + &
                    pot_ene%hbond + &
                    pot_ene%bond + &
                    pot_ene%angle + &
                    pot_ene%dihedral + &
                    pot_ene%vdw_14 + &
                    pot_ene%elec_14 + &
                    pot_ene%restraint + &
                    pot_ene%imp + &
                    pot_ene%angle_ub + &
                    pot_ene%efield + &
                    pot_ene%cmap

    if (iamd .gt. 1) then

      ! Calculate the boosting weight for amd
      call calculate_amd_dih_weights(atm_cnt, pot_ene%dihedral, frc, crd)
    end if

    if ( igamd.eq.2 .or. igamd.eq.3 .or. igamd.eq.5 .or. igamd.eq.7 .or. igamd.eq.9 .or. &
                  igamd.eq.11 ) then
      ! Calculate the boosting weight for gamd
      call calculate_gamd_dih_weights(atm_cnt, pot_ene%dihedral, frc, crd)
    end if

    ! Adjustment of total energy for constraint energies does not seem
    ! consistent, but it matches sander...

    call zero_pme_time()

    ! Finish up virial work; Timing is inconsequential...

    if (need_virials) then

      vir%atomic(:,:) = vir%elec_recip(:,:) + &
                        vir%elec_direct(:,:) + &
                        vir%elec_nb_adjust(:,:) + &
                        vir%elec_recip_vdw_corr(:,:) + &
                        vir%elec_recip_self(:,:) + &
                        vir%elec_14(:,:) + &
                        vir%ep_frame(:,:)

      if ( ips .gt. 0 ) then
        vir%atomic(:,:) = vir%atomic(:,:) + VIREXIPS(:,:)
      end if

      vir%molecular(:,:) = vir%molecular(:,:) + vir%atomic(:,:)

      ! Save virials in form used in runmd:

      virial(1) = 0.5d0 * vir%molecular(1, 1)
      virial(2) = 0.5d0 * vir%molecular(2, 2)
      virial(3) = 0.5d0 * vir%molecular(3, 3)

    end if

    vir_vs_ene = vir%elec_recip(1, 1) + &
                 vir%elec_recip(2, 2) + &
                 vir%elec_recip(3, 3) + &
                 vir%eedvir + &
                 vir%elec_nb_adjust(1, 1) + &
                 vir%elec_nb_adjust(2, 2) + &
                 vir%elec_nb_adjust(3, 3)

    ! Avoid divide-by-zero for pure neutral systems (l-j spheres), or if
    ! energies were not calculated this step...

    if (ti_mode .eq. 0) then
      if (pot_ene%elec_tot .ne. 0.0d0) then
        vir_vs_ene = abs(vir_vs_ene + pot_ene%elec_tot)/abs(pot_ene%elec_tot)
      else
        vir_vs_ene = 0.0d0
      end if
    else
      vir_vs_ene = 0.d0
      if (ti_vve(ti_ene0) .ne. 0.0d0) then
        vir_vs_ene = abs(ti_vve(ti_vir0) + &
          ti_vve(ti_ene0) * ti_weights(1))/abs(ti_vve(ti_ene0))
      end if
      if (ti_vve(ti_ene1) .ne. 0.d0) then
        vir_vs_ene = vir_vs_ene + abs(ti_vve(ti_vir1) + &
          ti_vve(ti_ene1) * ti_weights(2))/abs(ti_vve(ti_ene1))
      end if
    end if
    pme_err_est = vir_vs_ene

    ! Calculate the NMR restraint energy contributions, if requested.

    if (nmropt .ne. 0) then
      call nmr_calc(crd, frc, enmr, 6)
      pot_ene%restraint = pot_ene%restraint + enmr(1) + enmr(2) + enmr(3)
      pot_ene%total = pot_ene%total + enmr(1) + enmr(2) + enmr(3)
    end if

#ifndef NOXRAY
    ! Built-in X-ray target function and gradient
   call xray_get_derivative(crd, frc, nstep, xray_e)
       pot_ene%restraint = pot_ene%restraint + xray_e
       pot_ene%total = pot_ene%total + xray_e
#endif

    if (ti_mode .ne. 0 .and. need_pot_enes .and. nstep .ge. 0) then
      call ti_calc_dvdl
    end if

    if (verbose .gt. 0) then
      call write_netfrc(net_frcs)
      call pme_verbose_print(pot_ene, vir, vir_vs_ene)
    end if

    !AMD DUAL BOOST CALC START
    if (iamd .gt. 0) then
      call calculate_amd_total_weights(atm_cnt,pot_ene%total,pot_ene%dihedral,&
                                       pot_ene%amd_boost,frc,crd)
      ! Update total energy
      pot_ene%total = pot_ene%total + pot_ene%amd_boost
    end if

    !GAMD DUAL BOOST CALC START
    if (igamd .gt. 0) then
      if (igamd .eq. 1 .or. igamd .eq. 2 .or. igamd .eq. 3) then
        call calculate_gamd_total_weights(atm_cnt, pot_ene%total, pot_ene%dihedral, &
                                          pot_ene%gamd_boost, frc, crd)
        ! Update total energy
        pot_ene%total = pot_ene%total + pot_ene%gamd_boost
      else if (igamd .eq. 4 .or. igamd .eq. 5) then
        pot_ene_nb = pot_ene%elec_tot + pot_ene%vdw_tot + pot_ene%vdw_14 + pot_ene%elec_14 + pot_ene%efield
        nb_frc_tmp = nb_frc_gamd
        call calculate_gamd_nb_weights(atm_cnt, pot_ene_nb,pot_ene%dihedral,pot_ene%gamd_boost,nb_frc_gamd,crd)
        ! Update total force
        frc = frc - nb_frc_tmp + nb_frc_gamd
        ! Update total energy
        pot_ene%total = pot_ene%total + pot_ene%gamd_boost
      else if( (igamd.eq.6 .or. igamd.eq.7) .and. (ti_mode.ne.0)) then
        pot_ene_lig = ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                      ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                      ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
        pot_ene_dih = pot_ene%total + &
                      ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene)
        call calculate_gamd_ti_region_weights(atm_cnt,pot_ene_lig,pot_ene_dih,pot_ene%gamd_boost,fwgt,frc,crd,ti_sc_lst)
        ! Calculate net force on TI region
        if (mod(nstep+1,ntave).eq.0) then
          net_frc_ti    = 0.0
          do i = 1, atm_cnt
                if ( ti_sc_lst(i) .gt. 0 ) then
                        net_frc_ti(:) = net_frc_ti(:) + frc(:,i)
      end if
          end do
          write(*,'(a,i10,4f15.5)') "| GaMD: step, frc_weight, frc_lig(1:3) = ", nstep+1, fwgt, net_frc_ti
        end if
        ! Update total energy
        ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
      else if( (igamd.eq.8 .or. igamd.eq.9) .and. (ti_mode.ne.0)) then
        pot_ene_lig = ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene)
        pot_ene_dih = pot_ene%total + &
                      ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene) + &
                      ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                      ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
        call calculate_gamd_ti_region_weights(atm_cnt,pot_ene_lig,pot_ene_dih,pot_ene%gamd_boost,fwgt,frc,crd,ti_sc_lst)
      ! Update total energy
        ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
      else if( (igamd.eq.10 .or. igamd.eq.11) .and. (ti_mode.ne.0)) then
        pot_ene_lig = ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) + &
                      ti_others_ene(1,si_vdw_ene) + ti_others_ene(1,si_elect_ene)
        pot_ene_dih = pot_ene%total + &
                      ti_ene(1,si_bond_ene) + ti_ene(1,si_angle_ene) + ti_ene(1,si_dihedral_ene) + &
                      ti_ene(1,si_elect_14_ene) + ti_ene(1,si_vdw_14_ene)
        call calculate_gamd_ti_region_weights(atm_cnt,pot_ene_lig,pot_ene_dih,pot_ene%gamd_boost,fwgt,frc,crd,ti_sc_lst)
        ! Update total energy
        ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
      end if
    end if

    !scaledMD scale forces
    if (scaledMD .gt. 0) then
      call scaledMD_scale_frc(atm_cnt, pot_ene%total, frc, crd)
      pot_ene%total = pot_ene%total * scaledMD_lambda
    end if

    ! If belly is on then set the belly atom forces to zero:

    if (ibelly .gt. 0) call bellyf(atm_cnt, atm_igroup, frc)

    call update_time(nonbond_time)
    call update_pme_time(pme_misc_timer)
    if (charmm_active .and. do_charmm_dump_gold .eq. 1) then
      call pme_charmm_dump_gold(atm_cnt, frc, pot_ene)
      write(mdout, '(a)') 'charmm_gold() completed. Exiting'
      call mexit(6, 0)
    end if

    if (igamd .eq. 4 .or. igamd .eq. 5) then
      deallocate(nb_frc_gamd,nb_frc_tmp)
    end if

    if (infe.gt.0 .or. iextpot.gt.0) then
        bias_frc(:,:) = 0.d0
 
        if (infe.gt.0) then
           call nfe_on_force(crd, bias_frc, pot_ene%nfe)
        endif
        if (iextpot.gt.0) then
           temp_holding = 0.d0
           call pme_external(crd, bias_frc, temp_holding)
        endif

        do i = 1, atm_cnt
          frc(1,i) = frc(1,i) + bias_frc(1,i)
          frc(2,i) = frc(2,i) + bias_frc(2,i)
          frc(3,i) = frc(3,i) + bias_frc(3,i)
        end do

        if (infe.gt.0) then
           pot_ene%total     = pot_ene%total + pot_ene%nfe
           pot_ene%restraint = pot_ene%restraint + pot_ene%nfe
        endif
        if (iextpot.gt.0) then
           pot_ene%total     = pot_ene%total     + temp_holding
           pot_ene%restraint = pot_ene%restraint + temp_holding
        endif
    end if

    deallocate(img_frc)
! Ending if of external library
!end if
#endif /* CUDA */
    return


  end subroutine pme_force

#endif /* MPI */

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  dist_enes_virs_netfrcs
!
! Description: We reduce the appropriate subset of values in the ene array,
!              the ekcmt array, and the pme_ene_vir common block.
!*******************************************************************************

subroutine dist_enes_virs_netfrcs(pot_ene, vir, ekcmt, net_frcs, &
                                  need_pot_enes, need_virials)

  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  type(pme_pot_ene_rec) :: pot_ene
  type(pme_virial_rec)  :: vir
  double precision      :: ekcmt(3)
  double precision      :: net_frcs(3)
  logical, intent(in)   :: need_pot_enes
  logical, intent(in)   :: need_virials

! Local variables:

  type pme_dat
    sequence
    type(pme_pot_ene_rec)       :: pot_ene
    type(pme_virial_rec)        :: vir
    double precision            :: ekcmt(3)
    double precision            :: molvir_netfrc_corr(3,3)
    double precision            :: net_frcs(3)
  end type pme_dat

  integer, parameter            :: pme_dat_size = &
                                   pme_pot_ene_rec_size + &
                                   pme_virial_rec_size + 3 + 9 + 3

  integer, parameter            :: pme_dat_no_pot_ene_size = &
                                   pme_virial_rec_size + 3 + 9 + 3

  type(pme_dat), save           :: dat_in               ! used by all tasks
  type(pme_dat), save           :: dat_out              ! used by master

  integer                       :: buf_size

#ifdef COMM_TIME_TEST
  call start_test_timer(6, 'dist_enes_virs_netfrcs', 0)
#endif

  if (need_pot_enes .or. nmropt .ne. 0 .or. verbose .ne. 0 .or. infe .ne. 0) then

    ! When you need potential energies, you need parts of the virial info in
    ! the master and adding in net frc info is then a minor impact, so send
    ! it all.

    dat_in%pot_ene                 = pot_ene
    dat_in%vir                     = vir
    dat_in%ekcmt(:)                = ekcmt(:)
    dat_in%molvir_netfrc_corr(:,:) = molvir_netfrc_corr(:,:)
    dat_in%net_frcs(:)             = net_frcs(:)

    buf_size = pme_dat_size

    call mpi_allreduce(dat_in%pot_ene%total, dat_out%pot_ene%total, &
                       buf_size, mpi_double_precision, &
                       mpi_sum, pmemd_comm, err_code_mpi)

    pot_ene = dat_out%pot_ene
    vir     = dat_out%vir
    if (.not. usemidpoint) then
      ekcmt(:) = dat_out%ekcmt(:)
    end if
    molvir_netfrc_corr(:,:) = dat_out%molvir_netfrc_corr(:,:)
    net_frcs(:)             = dat_out%net_frcs(:)

  else if (need_virials) then

    ! The virials are distributed, with net frc info added as it has
    ! a minor impact.

    dat_in%vir                     = vir
    dat_in%ekcmt(:)                = ekcmt(:)
    dat_in%molvir_netfrc_corr(:,:) = molvir_netfrc_corr(:,:)
    dat_in%net_frcs(:)             = net_frcs(:)

    buf_size = pme_dat_no_pot_ene_size

    call mpi_allreduce(dat_in%vir%molecular(1, 1), &
                       dat_out%vir%molecular(1, 1), &
                       buf_size, mpi_double_precision, &
                       mpi_sum, pmemd_comm, err_code_mpi)

    vir                     = dat_out%vir
    ekcmt(:)                = dat_out%ekcmt(:)
    molvir_netfrc_corr(:,:) = dat_out%molvir_netfrc_corr(:,:)
    net_frcs(:)             = dat_out%net_frcs(:)

  else if (netfrc .ne. 0) then

    ! We just need net frc info!

    dat_in%net_frcs(:) = net_frcs(:)
    buf_size = 3

    call mpi_allreduce(dat_in%net_frcs(1), &
                       dat_out%net_frcs(1), &
                       buf_size, mpi_double_precision, &
                       mpi_sum, pmemd_comm, err_code_mpi)

    net_frcs(:) = dat_out%net_frcs(:)

  end if

#ifdef COMM_TIME_TEST
  call stop_test_timer(6)
#endif

  return

end subroutine dist_enes_virs_netfrcs
#endif

!*******************************************************************************
!
! Subroutine:  self
!
! Description: <TBS>
!
!*******************************************************************************

subroutine self(ene, ewaldcof, vol, vir, params_may_change)

  use gbl_constants_mod
  use prmtop_dat_mod
  use ti_mod
#ifdef TIDECOMP
  use ti_decomp_mod
#endif
#ifdef GTI
  use reaf_mod
#endif

  implicit none

! Formal arguments:

  double precision      :: ene, ewaldcof, vol, vir(3, 3)
  logical               :: params_may_change

! Local variables:

  integer                       :: i
#ifdef TIDECOMP
  integer                       :: j
#endif
  double precision              :: ee_plasma
  logical, save                 :: setup_not_done = .true.
  double precision, save        :: factor
  double precision, save        :: sqrt_pi
  double precision, save        :: sumq
  double precision, save        :: sumq2
  ! Only used for TI
  double precision, save        :: sumql(2)
  double precision, save        :: sumq2l(2)
  double precision, save        :: sumq2_sc(2)
  double precision              :: ti_ee_plasma(2)

! Only compute sumq and sumq2 at beginning. They don't change. This code is
! only executed by the master, so we precalc anything we can...

  if (setup_not_done .or. params_may_change) then

    factor = -0.5d0 * PI / (ewaldcof * ewaldcof)

    sqrt_pi = sqrt(PI)

    sumq = 0.d0
    sumq2 = 0.d0

    if (ti_mode .eq. 0) then
      do i = 1, natom
        sumq = sumq + atm_qterm(i)
        sumq2 = sumq2 + atm_qterm(i) * atm_qterm(i)
      end do
#ifdef GTI
      if (reaf_mode .ge. 0) then
        sumq = 0.d0
        sumq2 = 0.d0
        do i = 1, natom
          if (reaf_atom_list(reaf_mode+1,i)>0) then
            sumq = sumq + atm_qterm(i)*read_item_weights(5,1)  !TypeElecSC; internal
            sumq2 = sumq2 + atm_qterm(i) * atm_qterm(i) &
                * read_item_weights(5,1) * read_item_weights(5,1) 
          else
            sumq = sumq + atm_qterm(i)
            sumq2 = sumq2 + atm_qterm(i) * atm_qterm(i)
          end if
        end do
      endif
#endif
    else
#ifdef GTI
      if (reaf_mode .ge. 0) then
        sumq = 0.d0
        sumq2 = 0.d0
        do i = 1, natom
          if (reaf_atom_list(reaf_mode+1,i)>0) then
            if (ti_lst(2-reaf_mode,i) .eq. 0) then
              sumq = sumq + atm_qterm(i)*read_item_weights(5,1)  !TypeElecSC; internal
              sumq2 = sumq2 + atm_qterm(i) * atm_qterm(i) &
                * read_item_weights(5,1) * read_item_weights(5,1) 
            end if
          else
            if (ti_lst(2-reaf_mode,i) .eq. 0) then
              sumq = sumq + atm_qterm(i)
              sumq2 = sumq2 + atm_qterm(i) * atm_qterm(i)
            end if
          endif
        end do
      endif
#endif
      sumql(:) = 0.d0
      sumq2l(:) = 0.d0
      sumq2_sc(:)=0.0
      do i = 1, natom
        if (ti_lst(1,i) .ne. 0) then
          sumql(1) = sumql(1) + atm_qterm(i)
          sumq2l(1) = sumq2l(1) + atm_qterm(i) * atm_qterm(i)
          if (ti_sc_lst(i) .ne. 0) then
            sumq2_sc(1)=sumq2_sc(1)+atm_qterm(i) * atm_qterm(i)
          end if
        else if (ti_lst(2,i) .ne. 0) then
          sumql(2) = sumql(2) + atm_qterm(i)
          sumq2l(2) = sumq2l(2) + atm_qterm(i) * atm_qterm(i)
          if (ti_sc_lst(i) .ne. 0) then
            sumq2_sc(2)=sumq2_sc(2)+atm_qterm(i) * atm_qterm(i)
          end if
        else
          sumql(1) = sumql(1) + atm_qterm(i)
          sumq2l(1) = sumq2l(1) + atm_qterm(i) * atm_qterm(i)
          sumql(2) = sumql(2) + atm_qterm(i)
          sumq2l(2) = sumq2l(2) + atm_qterm(i) * atm_qterm(i)
        end if
      end do
    end if
    setup_not_done = .false.

  end if

#ifdef GTI
  if (ti_mode .eq. 0 .or. reaf_mode .ge. 0) then
#else
  if (ti_mode .eq. 0) then
#endif
    ee_plasma = factor * sumq * sumq / vol
    ene = - sumq2 * ewaldcof / sqrt_pi + ee_plasma
  else
    ti_pot(:) = 0.d0
    ee_plasma = 0.d0
    ti_ee_plasma(:) = 0.d0
    do i = 1, 2 !this is TI region to use
      ti_ee_plasma(i) = (factor * sumql(i) * sumql(i) / vol)
      ti_pot(i) = ti_ee_plasma(i) - sumq2l(i) * ewaldcof / sqrt_pi
#ifdef TIDECOMP
      ! This block of code is for decomp most likely incorrect
      do j = 1, natom
        if(ti_lst(i,j) .ne. 0 .or. ti_lst(3,j) .ne. 0) then
            lig_dvdl_decomp(id_atom_region(j)) = &
                            lig_dvdl_decomp(id_atom_region(j)) - &
                            atm_qterm(j) * atm_qterm(j) * ewaldcof / sqrt_pi * &
                            ti_item_dweights(3,i) + ti_ee_plasma(i)/natom * ti_signs(i)

            lig_ti_decomp(id_atom_region(j), j) = &
                            lig_ti_decomp(id_atom_region(j),j) - &
                            atm_qterm(j) * atm_qterm(j) * ewaldcof / sqrt_pi * &
                            ti_item_dweights(3,i) + ti_ee_plasma(i)/natom * ti_signs(i)

            ti_decomp(j) = ti_decomp(j) - atm_qterm(j) * atm_qterm(j) * ewaldcof / sqrt_pi * &
                            ti_item_dweights(3,i) + ti_ee_plasma(i)/natom * ti_signs(i)
        end if
      end do
#endif
      ee_plasma = ee_plasma + ti_ee_plasma(i) * ti_weights(i)
    end do
    call ti_update_ene_all(ti_pot, si_elect_ene, ene, 3)

    ti_ene_tmp(:,5) = - sumq2_sc(:) * ewaldcof / sqrt_pi
    ti_ene_tmp(:,4) = - ( sumq2l(:) - sumq2_sc(:) )* ewaldcof / sqrt_pi

    ti_ene_tmp(:,3) = ti_pot(:) -  ti_ene_tmp(:,4) - ti_ene_tmp(:,5)

    ti_ene_tmp(2,4) = ti_ene_tmp(2,4)-ti_ene_tmp(1,4)
    ti_ene_tmp(1,4) = 0.0

  end if

#ifdef CUDA
  call gpu_self(ee_plasma, ene, factor, vol, ewaldcof)
#endif

  ! The off-diagonal elements are already zero.

  vir(1,1) = -ee_plasma
  vir(2,2) = -ee_plasma
  vir(3,3) = -ee_plasma

  return

end subroutine self

!*******************************************************************************
!
! Subroutine:  self_gamd
!
! Description: <TBS>
!
!*******************************************************************************

subroutine self_gamd(ene, ewaldcof, vol, vir, params_may_change)

  use gbl_constants_mod
  use prmtop_dat_mod
  use ti_mod

  implicit none

! Formal arguments:

  double precision      :: ene, ewaldcof, vol, vir(3, 3)
  logical               :: params_may_change

! Local variables:

  integer                       :: i
  integer                       ::ti_mode_bk
  double precision              :: ee_plasma
  logical, save                 :: setup_not_done = .true.
  double precision, save        :: factor
  double precision, save        :: sqrt_pi
  double precision, save        :: sumq
  double precision, save        :: sumq2
  ! Only used for TI
  double precision, save        :: sumql(2)
  double precision, save        :: sumq2l(2)
  double precision              :: ti_ee_plasma(2)

! Only compute sumq and sumq2 at beginning. They don't change. This code is
! only executed by the master, so we precalc anything we can...

!! to get the cMD self energy set ti_mode=0
  ti_mode_bk=ti_mode
  ti_mode=0

  if (setup_not_done .or. params_may_change) then

    factor = -0.5d0 * PI / (ewaldcof * ewaldcof)

    sqrt_pi = sqrt(PI)

    sumq = 0.d0
    sumq2 = 0.d0

    if (ti_mode .eq. 0) then
      do i = 1, natom
        sumq = sumq + atm_qterm(i)
        sumq2 = sumq2 + atm_qterm(i) * atm_qterm(i)
      end do
    else
      sumql(:) = 0.d0
      sumq2l(:) = 0.d0
      do i = 1, natom
        if (ti_lst(1,i) .ne. 0) then
          sumql(1) = sumql(1) + atm_qterm(i)
          sumq2l(1) = sumq2l(1) + atm_qterm(i) * atm_qterm(i)
        else if (ti_lst(2,i) .ne. 0) then
          sumql(2) = sumql(2) + atm_qterm(i)
          sumq2l(2) = sumq2l(2) + atm_qterm(i) * atm_qterm(i)
        else
          sumql(1) = sumql(1) + atm_qterm(i)
          sumq2l(1) = sumq2l(1) + atm_qterm(i) * atm_qterm(i)
          sumql(2) = sumql(2) + atm_qterm(i)
          sumq2l(2) = sumq2l(2) + atm_qterm(i) * atm_qterm(i)
        end if
      end do
    end if
    setup_not_done = .false.

  end if

  if (ti_mode .eq. 0) then
    ee_plasma = factor * sumq * sumq / vol
    ene = - sumq2 * ewaldcof / sqrt_pi + ee_plasma
  else
    ti_pot(:) = 0.d0
    ee_plasma = 0.d0
    ti_ee_plasma(:) = 0.d0
    do i = 1, 2 !this is TI region to use
      ti_ee_plasma(i) = (factor * sumql(i) * sumql(i) / vol)
      ti_pot(i) = ti_ee_plasma(i) - sumq2l(i) * ewaldcof / sqrt_pi
      ee_plasma = ee_plasma + ti_ee_plasma(i) * ti_weights(i)
    end do
    call ti_update_ene_all(ti_pot, si_elect_ene, ene)
  end if

#ifdef CUDA
  call gpu_self(ee_plasma, ene, factor, vol, ewaldcof)
#endif

  ! The off-diagonal elements are already zero.

  vir(1,1) = -ee_plasma
  vir(2,2) = -ee_plasma
  vir(3,3) = -ee_plasma

  ti_mode = ti_mode_bk

  return

end subroutine self_gamd

!*******************************************************************************
!
! Subroutine:  vdw_correction
!
! Description:  Get analytic estimate of energy and virial corrections due to
!               dispersion interactions beyond the cutoff.
!*******************************************************************************

subroutine vdw_correction(ene, virial, params_may_change)

  use mdin_ctrl_dat_mod
  use pbc_mod
  use prmtop_dat_mod
  use ti_mod

  implicit none

! Formal arguments:

  double precision      :: ene, virial(3, 3)
  logical               :: params_may_change

! Local variables:

  integer                       :: i, j, ic, iaci
  logical, save                 :: setup_not_done = .true.
  double precision, save        :: ene_factor            ! Result of precalc.
  double precision              :: prefac, term
  ! Only for TI
  double precision, save        :: ti_ene_factor(2)
  integer                       :: types(ntypes + 1) !+1 for scratch space
  integer                       :: loc_nvdwcls(ntypes + 1) !+1 for scratch space
  integer                       :: itypes
  integer                       :: k
  double precision, save        :: reaf_ene_factor           
  double precision              :: term2, temp

! Only compute ene_factor at the beginning. It doesn't change. This code is
! only executed by the master, so we precalc anything we can...

  if (setup_not_done .or. params_may_change) then

    term = 0.d0

    ! Will later divide by volume, which is all that could change:

    prefac = 2.d0 * PI / (3.d0 * vdw_cutoff**3)

    if (ti_mode .eq. 0) then
      do i = 1, ntypes
        iaci = ntypes * (i - 1)
        do j = 1, ntypes
          ic = typ_ico(iaci + j)
          if (ic .gt. 0) term = term + &
                                gbl_nvdwcls(i) * gbl_nvdwcls(j) * gbl_cn2(ic)
        end do
      end do
      ene_factor = -prefac * term
    else
      do k = 1, 2 !this is the TI region to mask
        itypes = 1
        types(:) = 0
        loc_nvdwcls(:) = 0
        term = 0.d0
        ! All this does is create a sublist of vdw types used by atoms in the
        ! defined groups so the vdw correction term is for the subset of atoms
        do i = 1, natom
          if (ti_lst(k,i) .ne. 0) cycle !mask contribution from TI region k

          j = atm_iac(i)
          types(itypes) = j

          ! Also redo the atom count, since we are effectively changing
          ! the number of atoms in the system
          loc_nvdwcls(j) = loc_nvdwcls(j) + 1

          if (itypes .ge. 2) then
            do j = 1, itypes - 1
              if (types(j) .eq. types(itypes)) then
                itypes = itypes - 1 !don't advance the list
                exit
              end if
            end do
          end if
          itypes = itypes + 1
        end do

        itypes = itypes - 1 !the last entry is just scratch space

        do i = 1, itypes
          iaci = ntypes * (types(i) - 1)
          do j = 1, itypes
            ic = typ_ico(iaci + types(j))
            if (ic .gt. 0) term = term + loc_nvdwcls(types(i)) * &
                                         loc_nvdwcls(types(j)) * gbl_cn2(ic)
          end do
        end do
        if (usemidpoint) then
          ti_ene_factor(2-k+1) = -prefac * term
        else
          ti_ene_factor(3-k) = -prefac * term
        end if
      end do
    end if

#ifdef GTI
    if(reaf_mode .ge. 0) then

      loc_nvdwcls(1:ntypes) = 0
      do i = 1, natom
        if (reaf_atom_list(reaf_mode+1,i)>0) then
          j = atm_iac(i)
          loc_nvdwcls(j)=loc_nvdwcls(j)+1
        endif
      enddo
      
      !! re-do everything without REAF atoms
      term=0
      term2=0
      do i = 1, ntypes
        iaci = ntypes * (i - 1)
        do j = 1, ntypes
          ic = typ_ico(iaci + j)
          if (ic .gt. 0) then
            term = term + &
             (gbl_nvdwcls(i)-loc_nvdwcls(i)) *  &
             (gbl_nvdwcls(j)-loc_nvdwcls(j)) * gbl_cn2(ic)
            term2 = term2 + &
             loc_nvdwcls(i) * loc_nvdwcls(j) &
             * gbl_cn2(ic)            
          endif
        end do
      end do
      term = - prefac * term
      term2 = - prefac * term2      
      
      if (ti_mode .eq. 0) then
        temp=ene_factor
      else
        temp=ti_ene_factor(reaf_mode+1)
      endif
              
      temp = term + (temp -term) * read_item_weights(6,1) &
         + term2*(read_item_weights(6,2)-read_item_weights(6,1))
      
      if (ti_mode .eq. 0) then
        ene_factor=temp
      else
        ti_ene_factor(reaf_mode+1)=temp
      endif
  
    endif
#endif
    
    setup_not_done = .false.
  end if

  ene = ene_factor / uc_volume

  if (ti_mode .ne. 0) then
    ti_pot(:) = ti_ene_factor(:) / uc_volume
#ifdef GTI
    if (reaf_mode .ge. 0) then 
      temp = 0
    else 
      temp = ti_pot(1)
    endif
    ti_ene_tmp(:,6) =  ti_pot(:)-temp
#else
    ti_ene_tmp(1,6)=0.0
    ti_ene_tmp(2,6)= ti_pot(2) * ti_item_dweights(6,2) + ti_pot(1) * ti_item_dweights(6,1)
#endif    
    call ti_update_ene_all(ti_pot, si_vdw_ene, ene, 6) !write out to ene
  end if
#ifdef CUDA
  call gpu_vdw_correction(ene)
#endif


  ! The off-diagonal elements are already zero.

  virial(1, 1) = - 2.d0 * ene
  virial(2, 2) = - 2.d0 * ene
  virial(3, 3) = - 2.d0 * ene

  return

end subroutine vdw_correction

!*******************************************************************************
!
! Subroutine:  vdw_correction_gamd
!
! Description:  Get analytic estimate of energy and virial corrections due to
!               dispersion interactions beyond the cutoff.
!*******************************************************************************

subroutine vdw_correction_gamd(ene, virial, params_may_change)

  use mdin_ctrl_dat_mod
  use pbc_mod
  use prmtop_dat_mod
  use ti_mod

  implicit none

! Formal arguments:

  double precision      :: ene, virial(3, 3)
  logical               :: params_may_change

! Local variables:

  integer                       :: i, j, ic, iaci
  logical, save                 :: setup_not_done = .true.
  double precision, save        :: ene_factor            ! Result of precalc.
  double precision              :: prefac, term
  ! Only for TI
  double precision, save        :: ti_ene_factor(2)
  integer                       :: types(ntypes + 1) !+1 for scratch space
  integer                       :: loc_nvdwcls(ntypes + 1) !+1 for scratch space
  integer                       :: itypes
  integer                       :: k

! Only compute ene_factor at the beginning. It doesn't change. This code is
! only executed by the master, so we precalc anything we can...

  ti_mode=0
  setup_not_done = .true.
  if (setup_not_done .or. params_may_change) then

    term = 0.d0

    ! Will later divide by volume, which is all that could change:

    prefac = 2.d0 * PI / (3.d0 * vdw_cutoff**3)

    if (ti_mode .eq. 0) then
      do i = 1, ntypes
        iaci = ntypes * (i - 1)
        do j = 1, ntypes
          ic = typ_ico(iaci + j)
          if (ic .gt. 0) term = term + &
                                gbl_nvdwcls(i) * gbl_nvdwcls(j) * gbl_cn2(ic)
        end do
      end do
      ene_factor = -prefac * term
    else
      do k = 1, 2 !this is the TI region to mask
        itypes = 1
        types(:) = 0
        loc_nvdwcls(:) = 0
        term = 0.d0
        ! All this does is create a sublist of vdw types used by atoms in the
        ! defined groups so the vdw correction term is for the subset of atoms
        do i = 1, natom
          if (ti_lst(k,i) .ne. 0) cycle !mask contribution from TI region k

          j = atm_iac(i)
          types(itypes) = j


          ! Also redo the atom count, since we are effectively changing
          ! the number of atoms in the system
          loc_nvdwcls(j) = loc_nvdwcls(j) + 1

          if (itypes .ge. 2) then
            do j = 1, itypes - 1
              if (types(j) .eq. types(itypes)) then
                itypes = itypes - 1 !don't advance the list
                exit
              end if
            end do
          end if
          itypes = itypes + 1
        end do

        itypes = itypes - 1 !the last entry is just scratch space

        do i = 1, itypes
          iaci = ntypes * (types(i) - 1)
          do j = 1, itypes
            ic = typ_ico(iaci + types(j))
            if (ic .gt. 0) term = term + loc_nvdwcls(types(i)) * &
                                         loc_nvdwcls(types(j)) * gbl_cn2(ic)
          end do
        end do
        if (usemidpoint) then
          ti_ene_factor(2-k+1) = -prefac * term
        else
          ti_ene_factor(3-k) = -prefac * term
        end if
      end do
    end if
    setup_not_done = .false.
  end if

  ene = ene_factor / uc_volume

  if (ti_mode .ne. 0) then
    ti_pot(:) = ti_ene_factor(:) / uc_volume
    call ti_update_ene_all(ti_pot, si_vdw_ene, ene) !write out to ene
  end if
#ifdef CUDA
  call gpu_vdw_correction(ene)
#endif


  ! The off-diagonal elements are already zero.

  virial(1, 1) = - 2.d0 * ene
  virial(2, 2) = - 2.d0 * ene
  virial(3, 3) = - 2.d0 * ene

  return

end subroutine vdw_correction_gamd



!*******************************************************************************
!
! Subroutine:  respa_scale
!
! Description: <TBS>
!
!*******************************************************************************

subroutine respa_scale(atm_cnt, img_frc, nrespa)

  use img_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: img_frc(3, atm_cnt)
  integer               :: nrespa

! Local variables:

  integer               :: img_idx

#ifdef MPI
  integer               :: lst_idx
#endif /* MPI */

#ifdef MPI
  do lst_idx = 1, gbl_used_img_cnt
    img_idx = gbl_used_img_lst(lst_idx)
#else
  do img_idx = 1, atm_cnt
#endif
    img_frc(1, img_idx) = nrespa * img_frc(1, img_idx)
    img_frc(2, img_idx) = nrespa * img_frc(2, img_idx)
    img_frc(3, img_idx) = nrespa * img_frc(3, img_idx)
  end do

  return

end subroutine respa_scale

!*******************************************************************************
!
! Subroutine:  ti_respa_scale
!
! Description: <TBS>
!
!*******************************************************************************

subroutine ti_respa_scale(atm_cnt, img_frc, nrespa)

  use img_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: img_frc(2, 3, atm_cnt)
  integer               :: nrespa

! Local variables:

  integer               :: img_idx

#ifdef MPI
  integer               :: lst_idx
#endif /* MPI */

#ifdef MPI
  do lst_idx = 1, gbl_used_img_cnt
    img_idx = gbl_used_img_lst(lst_idx)
#else
  do img_idx = 1, atm_cnt
#endif
    img_frc(1, 1, img_idx) = nrespa * img_frc(1, 1, img_idx)
    img_frc(1, 2, img_idx) = nrespa * img_frc(1, 2, img_idx)
    img_frc(1, 3, img_idx) = nrespa * img_frc(1, 3, img_idx)
    img_frc(2, 1, img_idx) = nrespa * img_frc(2, 1, img_idx)
    img_frc(2, 2, img_idx) = nrespa * img_frc(2, 2, img_idx)
    img_frc(2, 3, img_idx) = nrespa * img_frc(2, 3, img_idx)
  end do

  return

end subroutine ti_respa_scale

#ifdef MPI

!*******************************************************************************
!
! Subroutine:  pme_bonded_force_midpoint
!
! Description: Bonded energy and force calcualtions
!
!*******************************************************************************

subroutine pme_bonded_force_midpoint(pot_ene, new_list)

  use bonds_midpoint_mod
  use angles_midpoint_mod
  use dihedrals_midpoint_mod
  use mdin_ctrl_dat_mod
  use nmr_calls_mod
  use parallel_dat_mod
  use timers_mod
  use processor_mod !, only : proc_atm_crd, proc_atm_frc, proc_atm_qterm, &
  !                         proc_num_atms
  use parallel_processor_mod

  implicit none

! Formal arguments:

  !double precision              :: proc_atm_crd(3, *)
  !double precision              :: proc_atm_frc(3, *)
  logical                       :: new_list
  type(pme_pot_ene_rec)         :: pot_ene

  ! These energy variables are temporaries, for summing. DON'T use otherwise!

  double precision              :: bond_ene
  double precision              :: ub_ene
  double precision              :: angle_ene
  double precision              :: dihedral_ene
  double precision              :: dihedral_imp_ene
  double precision              :: cmap_ene

#ifdef MPI
  !if (my_atm_cnt .eq. 0) return
  if (proc_num_atms .eq. 0) return
#endif /* MPI */
  bond_ene  = 0.d0
  angle_ene = 0.d0
  dihedral_ene = 0.d0
! Bond energy contribution:

! The ebdev/eadev stuff currently only is output under nmr_calls for non-mpi
! code, so we basically drop it here under mpi.
!if (mytaskid .eq. 0)
!print *, "Rank,cit_nbonh,cit_nbona",mytaskid,cit_nbonh,cit_nbona
#ifdef MPI
#else
  ebdev = 0.d0
#endif

#if 1

  if (ntf .le. 1) then
    if (cit_nbonh .gt. 0) then
      call get_bond_energy_midpoint(cit_nbonh, cit_h_bond, proc_atm_crd, proc_atm_frc, bond_ene,new_list,mult_vech_bond)
      pot_ene%bond = bond_ene
      !print *, "pme_force nbonh",pot_ene%bond, bond_ene, cit_nbonh, mytaskid
    end if
  end if

  if (ntf .le. 2) then
    if (cit_nbona .gt. 0) then
      call get_bond_energy_midpoint(cit_nbona, cit_a_bond, proc_atm_crd, proc_atm_frc, bond_ene,new_list,mult_veca_bond)
      pot_ene%bond = pot_ene%bond + bond_ene
      !print *, "pme_force nbona",pot_ene%bond, bond_ene, cit_nbona, mytaskid
    end if
  end if
#endif

#ifdef MPI
#else
    if (cit_nbonh + cit_nbona .gt. 0) &
      ebdev = sqrt(ebdev / (cit_nbonh + cit_nbona))
#endif
  call update_time(bond_time)

! Angle energy contribution:

#ifdef MPI
#else
  eadev = 0.d0
#endif
  if (ntf .le. 3) then
    if (cit_ntheth .gt. 0) then
      call get_angle_energy_midpoint(cit_ntheth, cit_h_angle, proc_atm_crd, proc_atm_frc, angle_ene,mult_vech_angle)
      pot_ene%angle = angle_ene
    end if
  end if

  if (ntf .le. 4) then
    if (cit_ntheta .gt. 0) then
      call get_angle_energy_midpoint(cit_ntheta, cit_a_angle, proc_atm_crd, proc_atm_frc, angle_ene,mult_veca_angle)
      pot_ene%angle = pot_ene%angle + angle_ene
    end if
  end if
#ifdef MPI
#else
  if (cit_ntheth + cit_ntheta .gt. 0) &
    eadev = 57.296 * sqrt(eadev / (cit_ntheth + cit_ntheta))
#endif

  call update_time(angle_time)

! Dihedral energy contribution:
  !print *, "PMEFORCE:dihedwih. w/o H:",cit_nphih, cit_nphia

  if (ntf .le. 5) then
    if (cit_nphih .gt. 0) then
      call get_dihed_midpoint_energy(cit_nphih, cit_h_dihed, proc_atm_crd, proc_atm_frc, dihedral_ene,mult_vech_dihed)
      pot_ene%dihedral = dihedral_ene
    end if
  end if

  if (ntf .le. 6) then
    if (cit_nphia .gt. 0) then
      call get_dihed_midpoint_energy(cit_nphia, cit_a_dihed, proc_atm_crd, proc_atm_frc, &
                            dihedral_ene,mult_veca_dihed)
      pot_ene%dihedral = pot_ene%dihedral + dihedral_ene
    end if
  end if

  call update_time(dihedral_time)


  return

end subroutine pme_bonded_force_midpoint

!*******************************************************************************
!
! Subroutine:  proc_alloc_force_mem
!
! Description: <TBS>
!
!*******************************************************************************

subroutine proc_alloc_force_mem
  use processor_mod
  implicit none

! Formal arguments:

  if (allocated(nb_frc))deallocate(nb_frc)
  allocate(nb_frc(3,proc_atm_alloc_size))
  return

end subroutine proc_alloc_force_mem

#endif

!*******************************************************************************
!
! Subroutine:  pme_bonded_force
!
! Description: <TBS>
!
!*******************************************************************************

subroutine pme_bonded_force(crd, frc, pot_ene)

  use angles_mod
  use angles_ub_mod
  use bonds_mod
  use constraints_mod
  use dihedrals_mod
  use dihedrals_imp_mod
  use cmap_mod
  use dynamics_mod
  use dynamics_dat_mod
  use mdin_ctrl_dat_mod
  use nmr_calls_mod
  use parallel_dat_mod
  use timers_mod

  implicit none

  ! Formal arguments:

  double precision              :: crd(3, *)
  double precision              :: frc(3, *)
  type(pme_pot_ene_rec)         :: pot_ene

  ! These energy variables are temporaries, for summing. DON'T use otherwise!

  double precision              :: bond_ene
  double precision              :: ub_ene
  double precision              :: angle_ene
  double precision              :: dihedral_ene
  double precision              :: dihedral_imp_ene
  double precision              :: cmap_ene

#ifdef MPI
  if (my_atm_cnt .eq. 0) return
#endif /* MPI */

  ! Bond energy contribution:

  ! The ebdev/eadev stuff currently only is output under nmr_calls for non-mpi
  ! code, so we basically drop it here under mpi.

#ifdef MPI
#else
  ebdev = 0.d0
#endif
  if (ntf .le. 1) then
    if (cit_nbonh .gt. 0) then
      call get_bond_energy(cit_nbonh, cit_h_bond, crd, frc, bond_ene)
      pot_ene%bond = bond_ene
    end if
  end if

  if (ntf .le. 2) then
    if (cit_nbona .gt. 0) then
      call get_bond_energy(cit_nbona, cit_a_bond, crd, frc, bond_ene)
      pot_ene%bond = pot_ene%bond + bond_ene
    end if
  end if
#ifdef MPI
#else
  if (cit_nbonh + cit_nbona .gt. 0) &
    ebdev = sqrt(ebdev / (cit_nbonh + cit_nbona))
#endif

  call update_time(bond_time)

  ! UB Angle energy contribution
  if (cit_nthet_ub .gt. 0) then
    call get_angle_ub_energy(cit_nthet_ub, cit_angle_ub, crd, frc, ub_ene)
    pot_ene%angle_ub = ub_ene
  end if


  ! Angle energy contribution:

#ifdef MPI
#else
  eadev = 0.d0
#endif
  if (ntf .le. 3) then
    if (cit_ntheth .gt. 0) then
      call get_angle_energy(cit_ntheth, cit_angle, crd, frc, angle_ene)
      pot_ene%angle = angle_ene
    end if
  end if

  if (ntf .le. 4) then
    if (cit_ntheta .gt. 0) then
      call get_angle_energy(cit_ntheta, cit_angle(cit_ntheth+1), &
        crd, frc, angle_ene)
      pot_ene%angle = pot_ene%angle + angle_ene
    end if
  end if
#ifdef MPI
#else
  if (cit_ntheth + cit_ntheta .gt. 0) &
    eadev = 57.296 * sqrt(eadev / (cit_ntheth + cit_ntheta))
#endif

  call update_time(angle_time)

  ! Improper dihedral energy contributions:

  if (cit_nimphi .gt. 0) then
    call get_dihed_imp_energy(cit_nimphi, cit_dihed_imp, crd, frc, dihedral_imp_ene)
    pot_ene%imp = dihedral_imp_ene
  end if

  if (cit_cmap_term_count .gt. 0) then
    call get_cmap_energy(cit_cmap_term_count, cit_cmap, crd, frc, cmap_ene)
    pot_ene%cmap = cmap_ene
  end if



  ! Dihedral energy contribution:

  if (iamd.gt.1) then
    if (ntf .le. 5) then
      if (cit_nphih .gt. 0) then
        call get_dihed_energy_amd(cit_nphih, cit_dihed, crd, dihedral_ene)
        pot_ene%dihedral = dihedral_ene
      end if
    end if

    if (ntf .le. 6) then
      if (cit_nphia .gt. 0) then
        call get_dihed_energy_amd(cit_nphia, cit_dihed(cit_nphih + 1), crd, &
                                  dihedral_ene)
        pot_ene%dihedral = pot_ene%dihedral + dihedral_ene
      end if
    end if
 else if((igamd.eq.2).or.(igamd.eq.3).or.(igamd.eq.5).or.(igamd.eq.7).or.(igamd.eq.9))then
    if (ntf .le. 5) then
      if (cit_nphih .gt. 0) then
        call get_dihed_energy_gamd(cit_nphih, cit_dihed, crd, dihedral_ene)
        pot_ene%dihedral = dihedral_ene
      end if
    end if

    if (ntf .le. 6) then
      if (cit_nphia .gt. 0) then
        call get_dihed_energy_gamd(cit_nphia, cit_dihed(cit_nphih + 1), crd, &
                                   dihedral_ene)
        pot_ene%dihedral = pot_ene%dihedral + dihedral_ene
      end if
    end if
  else
    if (ntf .le. 5) then
      if (cit_nphih .gt. 0) then
        call get_dihed_energy(cit_nphih, cit_dihed, crd, frc, dihedral_ene)
        pot_ene%dihedral = dihedral_ene
      end if
    end if

    if (ntf .le. 6) then
      if (cit_nphia .gt. 0) then
        call get_dihed_energy(cit_nphia, cit_dihed(cit_nphih + 1), crd, frc, &
                              dihedral_ene)
        pot_ene%dihedral = pot_ene%dihedral + dihedral_ene
      end if
    end if
  end if

  call update_time(dihedral_time)

  ! Calculate the position constraint energy:

  if (natc .gt. 0) then
    call get_crd_constraint_energy(natc, pot_ene%restraint, atm_jrc, &
                                   crd, frc, atm_xc, atm_weight)
  end if

  return

end subroutine pme_bonded_force


!*******************************************************************************
!
! Subroutine:  write_netfrc
!
! Description:  Get the netfrc's back into the external axis order and print
!               them out.  We do this all in a separate subroutine just to
!               keep from cluttering up critical code.
!
!*******************************************************************************

subroutine write_netfrc(net_frcs)

  use axis_optimize_mod
  use file_io_dat_mod

  implicit none

! Formal arguments:

  double precision      :: net_frcs(3)

! Local variables:

  integer               :: ord1, ord2, ord3

  ord1 = axis_flipback_ords(1)
  ord2 = axis_flipback_ords(2)
  ord3 = axis_flipback_ords(3)

  write(mdout, 33) net_frcs(ord1), net_frcs(ord2), net_frcs(ord3)

  return

33     format(1x, 'NET FORCE PER ATOM: ', 3(1x, e12.4))

end subroutine write_netfrc

!*******************************************************************************
!
! Subroutine:  pme_verbose_print
!
! Description: <TBS>
!
!*******************************************************************************

subroutine pme_verbose_print(pot_ene, vir, vir_vs_ene)

  use axis_optimize_mod
  use file_io_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  type(pme_pot_ene_rec) :: pot_ene
  type(pme_virial_rec)  :: vir
  double precision      :: vir_vs_ene

! Local variables:

  integer       :: ord1, ord2, ord3

  if (.not. master) return

  ord1 = axis_flipback_ords(1)
  ord2 = axis_flipback_ords(2)
  ord3 = axis_flipback_ords(3)

  if (verbose .ge. 1) then
    write(mdout, '(3(/,5x,a,f24.12))') &
          'Evdw                   = ', pot_ene%vdw_tot, &
          'Ehbond                 = ', pot_ene%hbond, &
          'Ecoulomb               = ', pot_ene%elec_tot
    write(mdout, '(2(/,5x,a,f24.12))') &
          'Iso virial             = ',  &
          vir%molecular(1, 1) + vir%molecular(2, 2) + vir%molecular(3, 3), &
          'Eevir vs. Ecoulomb     = ', vir_vs_ene
  end if

  if (verbose .ge. 2) then
    write(mdout, '(4(/,5x,a,f24.12),/)') &
          'E electrostatic (self) = ', pot_ene%elec_self, &
          '                (rec)  = ', pot_ene%elec_recip, &
          '                (dir)  = ', pot_ene%elec_dir, &
          '                (adj)  = ', pot_ene%elec_nb_adjust
    write(mdout, 30) vir%molecular(ord1, ord1), &
                     vir%molecular(ord1, ord2), &
                     vir%molecular(ord1, ord3)
    write(mdout, 30) vir%molecular(ord2, ord1), &
                     vir%molecular(ord2, ord2), &
                     vir%molecular(ord2, ord3)
    write(mdout, 30) vir%molecular(ord3, ord1), &
                     vir%molecular(ord3, ord2), &
                     vir%molecular(ord3, ord3)
30     format(5x, 'MOLECULAR VIRIAL: ', 3(1x, e14.7))
  end if

  if (verbose .eq. 3) then
    write(mdout, *) '--------------------------------------------'
    write(mdout, 31) vir%elec_recip(ord1, ord1), &
                     vir%elec_recip(ord1, ord2), &
                     vir%elec_recip(ord1, ord3)
    write(mdout, 31) vir%elec_recip(ord2, ord1), &
                     vir%elec_recip(ord2, ord2), &
                     vir%elec_recip(ord2, ord3)
    write(mdout, 31) vir%elec_recip(ord3, ord1), &
                     vir%elec_recip(ord3, ord2), &
                     vir%elec_recip(ord3, ord3)
    write(mdout, *) '..................'
31     format(5x, 'Reciprocal VIRIAL: ', 3(1x, e14.7))
    write(mdout, 32) vir%elec_direct(ord1, ord1), &
                     vir%elec_direct(ord1, ord2), &
                     vir%elec_direct(ord1, ord3)
    write(mdout, 32) vir%elec_direct(ord2, ord1), &
                     vir%elec_direct(ord2, ord2), &
                     vir%elec_direct(ord2, ord3)
    write(mdout, 32) vir%elec_direct(ord3, ord1), &
                     vir%elec_direct(ord3, ord2), &
                     vir%elec_direct(ord3, ord3)
    write(mdout, *) '..................'
32     format(5x, 'Direct VIRIAL: ', 3(1x, e14.7))
    write(mdout, 38) vir%eedvir
    write(mdout, *) '..................'
38     format(5x, 'Dir Sum EE vir trace: ', e14.8)
    write(mdout, 33) vir%elec_nb_adjust(ord1, ord1), &
                     vir%elec_nb_adjust(ord1, ord2), &
                     vir%elec_nb_adjust(ord1, ord3)
    write(mdout, 33) vir%elec_nb_adjust(ord2, ord1), &
                     vir%elec_nb_adjust(ord2, ord2), &
                     vir%elec_nb_adjust(ord2, ord3)
    write(mdout, 33) vir%elec_nb_adjust(ord3, ord1), &
                     vir%elec_nb_adjust(ord3, ord2), &
                     vir%elec_nb_adjust(ord3, ord3)
    write(mdout, *) '..................'
33     format(5x, 'Adjust VIRIAL: ', 3(1x, e14.7))
    write(mdout, 34) vir%elec_recip_vdw_corr(ord1, ord1), &
                     vir%elec_recip_vdw_corr(ord1, ord2), &
                     vir%elec_recip_vdw_corr(ord1, ord3)
    write(mdout, 34) vir%elec_recip_vdw_corr(ord2, ord1), &
                     vir%elec_recip_vdw_corr(ord2, ord2), &
                     vir%elec_recip_vdw_corr(ord2, ord3)
    write(mdout, 34) vir%elec_recip_vdw_corr(ord3, ord1), &
                     vir%elec_recip_vdw_corr(ord3, ord2), &
                     vir%elec_recip_vdw_corr(ord3, ord3)
    write(mdout, *) '..................'
34     format(5x, 'Recip Disp. VIRIAL: ', 3(1x, e14.7))
    write(mdout, 35) vir%elec_recip_self(ord1, ord1), &
                     vir%elec_recip_self(ord1, ord2), &
                     vir%elec_recip_self(ord1, ord3)
    write(mdout, 35) vir%elec_recip_self(ord2, ord1), &
                     vir%elec_recip_self(ord2, ord2), &
                     vir%elec_recip_self(ord2, ord3)
    write(mdout, 35) vir%elec_recip_self(ord3, ord1), &
                     vir%elec_recip_self(ord3, ord2), &
                     vir%elec_recip_self(ord3, ord3)
    write(mdout, *) '..................'
35     format(5x, 'Self VIRIAL: ', 3(1x, e14.7))
    write(mdout, 36) vir%elec_14(ord1, ord1), &
                     vir%elec_14(ord1, ord2), &
                     vir%elec_14(ord1, ord3)
    write(mdout, 36) vir%elec_14(ord2, ord1), &
                     vir%elec_14(ord2, ord2), &
                     vir%elec_14(ord2, ord3)
    write(mdout, 36) vir%elec_14(ord3, ord1), &
                     vir%elec_14(ord3, ord2), &
                     vir%elec_14(ord3, ord3)
    write(mdout, *) '..................'
36     format(5x, 'E14 VIRIAL: ', 3(1x, e14.7))
    write(mdout, 37) vir%atomic(ord1, ord1), &
                     vir%atomic(ord1, ord2), &
                     vir%atomic(ord1, ord3)
    write(mdout, 37) vir%atomic(ord2, ord1), &
                     vir%atomic(ord2, ord2), &
                     vir%atomic(ord2, ord3)
    write(mdout, 37) vir%atomic(ord3, ord1), &
                     vir%atomic(ord3, ord2), &
                     vir%atomic(ord3, ord3)
37     format(5x, 'Atomic VIRIAL: ', 3(1x, e14.7))
    write(mdout, *)'--------------------------------------------'
  end if

  return

end subroutine pme_verbose_print

end module pme_force_mod
