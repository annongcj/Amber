#include "copyright.i"

!*******************************************************************************
!
! Module: gb_force_mod
!
! Description: <TBS>
!
!*******************************************************************************

module gb_force_mod

  use energy_records_mod, only : gb_pot_ene_rec, gb_pot_ene_rec_size, &
                                 null_gb_pot_ene_rec
  use gbl_datatypes_mod
  use external_mod, only : gb_external

! PHMD
  use phmd_mod

  implicit none

#ifdef _OPENMP_

#ifdef __INTEL_COMPILER
#define LENGTH atm_cnt
#else
#define LENGTH *
#endif /* INTEL_COMPILER */

  abstract interface
      subroutine func(crd, frc, rborn, fs, charge, iac, ico, numex, &
                 natex, atm_cnt, natbel, egb, eelt, evdw, esurf, irespa, &
                 skip_radii_)
        integer               :: atm_cnt
        double precision      :: rborn(LENGTH)
        double precision      :: fs(LENGTH)
        double precision      :: charge(LENGTH)
        integer               :: iac(LENGTH)
        integer               :: ico(LENGTH)
        integer               :: natbel
        integer, intent(in)   :: irespa
        logical               :: skip_radii ! for converting the optional argument
        double precision      :: egb, eelt, evdw, esurf
        logical, optional, intent(in)         :: skip_radii_
        integer               :: numex(*)
        integer               :: natex(*)
#ifdef __INTEL_COMPILER
        double precision      :: crd(atm_cnt*3)
        double precision,  intent(inout)      :: frc(atm_cnt*3)
#else
        double precision      :: crd(*)
        double precision      :: frc(*)
#endif
      end subroutine func
   end interface

   procedure (func), pointer :: force_ptr => null ()
   procedure (func), pointer :: energy_ptr => null ()
   procedure (func), pointer :: gbene_ptr => null ()


  ! Potential energies, with breakdown, from GB.  This is intended to be the
  ! external interface to potential energy data produced by this module, in
  ! particular by subroutine gb_force().

contains

subroutine final_frc_setup()
  use ti_mod
  use mdin_ctrl_dat_mod
  use remd_mod
  use gb_ene_hybrid_mod
  use gb_ene_mod
  gbene_ptr => gb_ene

  if (using_pme_potential .and. (icnstph .eq. 2 .or. icnste .eq. 2 .or. hybridgb .gt. 0)) then
    energy_ptr => gb_ene
    force_ptr => gb_ene
  else
    if (ifsc .gt. 0) then
      energy_ptr => gb_ene_hyb_energy_ifsc
      force_ptr => gb_ene_hyb_energy_ifsc
    else
      if (ti_mode .ne. 0) then
        energy_ptr => gb_ene_hyb_energy_timode
        if ((numexchg .ne. 0 .and. remd_method .ne. 0) .or. calc_emil .eqv. .true.) then
          force_ptr => gb_ene_hyb_energy_timode
        else
          force_ptr => gb_ene_hyb_force_timode
        end if
      else
        energy_ptr => gb_ene_hyb_energy
        if ((numexchg .ne. 0 .and. remd_method .ne. 0) .or. calc_emil .eqv. .true.) then
          force_ptr => gb_ene_hyb_energy
        else
          force_ptr => gb_ene_hyb_force
        end if
      end if
    end if
  end if ! pme_potential .and. icnstph .eq. 2
end subroutine
#else
contains
#endif /*_OPENMP_*/

!*******************************************************************************
!
! Subroutine:  gb_force
!
! Description: <TBS>
!
!*******************************************************************************

subroutine gb_force(atm_cnt, crd, frc, pot_ene, irespa, need_pot_ene)

  use angles_mod
  use angles_ub_mod
  use bonds_mod
  use charmm_mod, only : charmm_active
  use constantph_dat_mod, only : on_cpstep, proposed_qterm, cnstph_frc
  use constante_dat_mod, only : on_cestep, cnste_frc
  use constraints_mod
  use cmap_mod
  use dihedrals_mod
  use dihedrals_imp_mod
  use emap_mod,only:emapforce
  use extra_pnts_nb14_mod
#ifdef _OPENMP_
  use gb_ene_hybrid_mod
#endif
  use gb_ene_mod
  use gb_parallel_mod
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod, only : frameon
  use mdin_debugf_dat_mod, only : do_charmm_dump_gold
  use nmr_calls_mod
  use parallel_dat_mod
  use parallel_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
  use runfiles_mod
  use timers_mod
  use charmm_mod, only : charmm_active
  use mdin_debugf_dat_mod, only : do_charmm_dump_gold
  use amd_mod
  use gamd_mod
  use scaledMD_mod
  use gbsa_mod
  use ti_mod
  use nfe_setup_mod, only : nfe_on_force => on_force

#ifdef MPI
  use neb_mod
  use nebread_mod
  use multipmemd_mod
#endif /* MPI */
  !----------------------------------------------------

  implicit none

! Formal arguments:

  integer                       :: atm_cnt
  double precision              :: crd(3, atm_cnt)
  double precision              :: frc(3, atm_cnt)
  type(gb_pot_ene_rec)          :: pot_ene
  type(gb_pot_ene_rec)          :: pot_ene_prop ! CpH or CE proposed state energy rec.
  integer, intent(in)           :: irespa
  logical, intent(in)           :: need_pot_ene

! Local variables:

#ifdef MPI
  integer, dimension(MPI_STATUS_SIZE) :: st           !DG
#endif /* MPI */
  integer                       :: iatm, jatm   !added by DG
  double precision              :: enmr(3)
  double precision              :: enemap
  integer                       :: i, j
  double precision              :: temp_tot_dih
  double precision              :: totdih
  integer                       :: buf_size

  ! Constant pH energy terms for proposed state

  double precision              :: cph_egb
  double precision              :: cph_eel
  double precision              :: cph_eel14
  double precision              :: cph_vdw
  double precision              :: cph_vdw14
  double precision              :: cph_surf

  ! End constant pH energy terms for proposed state

  ! Constant Redox potential energy terms for proposed state

  double precision              :: ce_egb
  double precision              :: ce_eel
  double precision              :: ce_eel14
  double precision              :: ce_vdw
  double precision              :: ce_vdw14
  double precision              :: ce_surf

  ! End constant Redox potential energy terms for proposed state

! biased forces, added by FENG PAN
  double precision              :: bias_frc(3,atm_cnt)

  ! Placeholder for virial computations so that orient_frc (intended for PME simulations)
  ! can be called
  double precision              :: vir_placeholder(3, 3)
  
! GaMD Local Variables
  integer                       :: alloc_failed
  double precision              :: pot_ene_nb, pot_ene_lig, pot_ene_dih, fwgt
  double precision, allocatable :: nb_frc_gamd(:,:)
  double precision, allocatable :: nb_frc_tmp(:,:)

! temp holder for sc non bonded terms calculated in gb_ene_sc
  double precision              :: gb_sc_eel
  double precision              :: gb_sc_evdw

! temp holder for ene from external force modifiers
  double precision              :: temp_holder

#ifdef CUDA


  ! zero-out the potential energy record: (before adding any extra in external)
  if (nmropt .ne. 0) then
    pot_ene = null_gb_pot_ene_rec !Zero the entire structure
    enmr(:) = 0
  end if

  ! Update any weight modifications
  if (nmropt .ne. 0) then
    call nmr_weight(atm_cnt, crd, mdout)
  end if

  if (need_pot_ene) then !DO NOT CHANGE THIS OR ADD ADDITIONAL CLAUSES, THIS IS SET ONCE IN RUNMD.F90.
#ifdef MPI
    if((iamd.eq.2).or.(iamd.eq.3))then
      totdih = 0.0
      if(num_amd_lag .eq. 0)then
        temp_tot_dih = 0.0
        !Note this double calculates the dihedral energy since we do it
        !here for AMD and then repeat it later - consider merging these
        !at some point for efficiency.
        call gpu_calculate_gb_amd_dihedral_energy(totdih)
        buf_size = 1
        call mpi_allreduce(totdih, temp_tot_dih, &
                       buf_size, mpi_double_precision, &
                       mpi_sum, pmemd_comm, err_code_mpi)
        totdih = temp_tot_dih
      end if
      call gpu_calculate_amd_dihedral_weight(totdih)
    end if
#else
    if((iamd.eq.2).or.(iamd.eq.3))then
      call gpu_calculate_gb_amd_dihedral_energy_weight()
    end if
#endif

#ifdef MPI
    if((igamd.eq.2).or.(igamd.eq.3).or.(igamd.eq.5))then
      totdih = 0.0
      if(num_gamd_lag .eq. 0)then
        temp_tot_dih = 0.0
        !Note this double calculates the dihedral energy since we do it
        !here for GaMD and then repeat it later - consider merging these
        !at some point for efficiency.
        call gpu_calculate_gb_gamd_dihedral_energy(totdih)
        buf_size = 1
        call mpi_allreduce(totdih, temp_tot_dih, &
                       buf_size, mpi_double_precision, &
                       mpi_sum, pmemd_comm, err_code_mpi)
        totdih = temp_tot_dih
      end if
      call gpu_calculate_gamd_dihedral_weight(totdih)
    end if
#else
    if((igamd.eq.2).or.(igamd.eq.3).or.(igamd.eq.5))then
      call gpu_calculate_gb_gamd_dihedral_energy_weight()
    end if
#endif

! BEGIN GBSA implementation to work with the GPU code.
! It runs on the CPU for now, needs to write a kernel for this.
! Regular CPU call is located in gb_ene
! pmemd.cuda only takes nrespa=1 for now, so GBSA contribution must ALWAYS
! be calculated (regardless of need_pot_ene)
    if (gbsa .eq. 1) then
      frc=0.d0
      gbsafrc=0.d0
      call gpu_download_crd(crd)
      call gbsa_ene(crd, gbsafrc, pot_ene%surf ,atm_cnt, jj, r2x, belly_atm_cnt)
    end if
! END GBSA implementation to work with the GPU code.
    if (on_cpstep) then
      call gpu_refresh_charges(cit_nb14, gbl_one_scee, proposed_qterm)
      call gpu_gb_ene(pot_ene_prop, enmr, ineb)
      call gpu_refresh_charges(cit_nb14, gbl_one_scee, atm_qterm)
      call gpu_gb_ene(pot_ene, enmr, ineb)
      pot_ene%dvdl = pot_ene_prop%total - pot_ene%total
    else if (on_cestep) then
      call gpu_refresh_charges(cit_nb14, gbl_one_scee, proposed_qterm)
      call gpu_gb_ene(pot_ene_prop, enmr, ineb)
      call gpu_refresh_charges(cit_nb14, gbl_one_scee, atm_qterm)
      call gpu_gb_ene(pot_ene, enmr, ineb)
      pot_ene%dvdl = pot_ene_prop%total - pot_ene%total
    else
      call gpu_gb_ene(pot_ene, enmr, ineb)
    end if
! BEGIN GBSA implementation to work with the GPU code.
! force elements sent to the GPU to be added to the other force contributions
    if (gbsa .eq. 1) then
#ifdef MPI
      buf_size = atm_cnt *3
        call mpi_allreduce(gbsafrc, frc, &
                       buf_size, mpi_double_precision, &
                       mpi_sum, pmemd_comm, err_code_mpi)
      call gpu_gbsa_frc_add(frc)
#else
      call gpu_gbsa_frc_add(gbsafrc)
#endif
      pot_ene%total =  pot_ene%total +  pot_ene%surf
    end if

! END GBSA implementation to work with the GPU code.
#ifdef MPI
    call gb_distribute_enes(pot_ene)
#endif
  else !need_pot_ene

! BEGIN GBSA implementation to work with the GPU code.
! It runs on the CPU for now, needs to write a kernel for this.
! Regular CPU call is located in gb_ene
! pmemd.cuda only takes nrespa=1 for now, so GBSA contribution must ALWAYS
! be calculated (regardless of need_pot_ene)
    if (gbsa .eq. 1) then
      frc=0.d0
      gbsafrc=0.d0
      call gpu_download_crd(crd)
      call gbsa_ene(crd, gbsafrc, pot_ene%surf, atm_cnt, jj, r2x, belly_atm_cnt)
    end if
! END GBSA implementation to work with the GPU code.
    call gpu_gb_forces()
! BEGIN GBSA implementation to work with the GPU code.
! force elements sent to the GPU to be added to the other force contributions
    if (gbsa .eq. 1) then
#ifdef MPI
      buf_size = atm_cnt *3
        call mpi_allreduce(gbsafrc, frc, &
                       buf_size, mpi_double_precision, &
                       mpi_sum, pmemd_comm, err_code_mpi)
     call gpu_gbsa_frc_add(frc)
#else
     call gpu_gbsa_frc_add(gbsafrc)
#endif
    end if
! END GBSA implementation to work with the GPU code.
  end if
  if (nmropt .ne. 0) then
    call nmr_calc(crd, frc, enmr, 6)
  end if
  if(iamd.gt.0)then
!AMD calculate weight and scale forces
    call gpu_calculate_and_apply_amd_weights(pot_ene%total, pot_ene%dihedral, &
                                             pot_ene%amd_boost, num_amd_lag)
     ! Update total energy
     pot_ene%total = pot_ene%total + pot_ene%amd_boost
  end if
  if(igamd.gt.0)then
!GaMD calculate weight and scale forces
    if(igamd.eq.1 .or. igamd.eq.2 .or. igamd.eq.3) &
       call gpu_calculate_and_apply_gamd_weights_gb(pot_ene%total, pot_ene%dihedral, &
                                               pot_ene%gamd_boost,num_gamd_lag)
     if(igamd.eq.4 .or. igamd.eq.5 ) then       ! THIS DOES NOT work as pNBForce is not saved separately with GB.
       pot_ene_nb = pot_ene%elec_tot + pot_ene%vdw_tot + pot_ene%gb + pot_ene%vdw_14 + pot_ene%elec_14
       call gpu_calculate_and_apply_gamd_weights_gb_nb(pot_ene_nb, pot_ene%dihedral, &
                                               pot_ene%gamd_boost,num_gamd_lag)
     end if
     ! Update total energy
     pot_ene%total = pot_ene%total + pot_ene%gamd_boost
  end if

  if(scaledMD.gt.0)then
!scaledMD scale forces
    call gpu_scaledmd_scale_frc(pot_ene%total)
    pot_ene%total = pot_ene%total * scaledMD_lambda
  end if

  call update_time(nonbond_time)

  if (infe.gt.0) then
      call gpu_download_crd(crd)

      bias_frc(:,:) = 0.d0
      call nfe_on_force(crd,bias_frc,pot_ene%nfe)

      pot_ene%total = pot_ene%total + pot_ene%nfe
      pot_ene%restraint = pot_ene%restraint + pot_ene%nfe

      call gpu_upload_frc_add(bias_frc)
  end if

  
  ! If calling an External library
  if (iextpot .gt. 0) then
#ifdef MPI
      if (mytaskid.eq.0) then
#endif /* MPI */
         call gpu_download_crd(crd) !avoidable performance hit, but simpler code
         temp_holder = 0.0
         bias_frc(:,:) = 0.d0 !!reuse the force holder from nfe.
         call gb_external(crd, bias_frc, temp_holder) !Save "external" energyinto the restraint field.
         pot_ene%total      = pot_ene%total     + temp_holder
         pot_ene%restraint  = pot_ene%restraint + temp_holder

         call gpu_upload_frc_add(bias_frc)
#ifdef MPI
      endif
#endif /* MPI */
  ! End if not calling an External library
  endif




#ifdef MPI
  if (ineb>0) then
   call transfer_fit_neb_crd(irespa)
   call full_neb_forces(irespa)

   if (beadid==1 .or. beadid==neb_nbead) then
      call gpu_clear_forces()
   end if

  end if
#endif /* MPI */
!END_NEB

  return
#else /* below here, not CUDA */

! Zero energies that are stack or call parameters:

  pot_ene = null_gb_pot_ene_rec
  pot_ene_prop = null_gb_pot_ene_rec

! Zero internal energies:

  enmr(:) = 0.d0

  if (ti_mode .ne. 0) then
    call ti_zero_arrays
    ti_nb_frc(:,:,:) = 0.d0
  end if

#ifdef _OPENMP_
  if  (mod(irespa, ntpr) .ne. 0) then
      gbene_ptr => force_ptr
  else
      gbene_ptr => energy_ptr
  end if
#endif /*_OPENMP_*/

! Do weight changes, if requested.

  if (nmropt .ne. 0) call nmr_weight(atm_cnt, crd, mdout)

! If no force calcs are to be done, clear the frc array and bag out now.

  if (ntf .eq. 8) then
    frc(:,:) = 0.d0
    return
  end if

  call zero_time()
  call zero_gb_time()

! GaMD
  if (igamd.eq.4 .or. igamd.eq.5) then
   allocate( nb_frc_gamd(3, atm_cnt), &
             nb_frc_tmp(3, atm_cnt), &
            stat = alloc_failed)
   if (alloc_failed .ne. 0) call setup_alloc_error
  end if

! Calculate the non-bonded contributions:

  frc(:,:) = 0.d0

  ! If calling an External library
  if (iextpot .gt. 0) then
       !Save "external" energy into the restraint field.
       !Even if we eventually get comfortable running restraints+external, it should
       !still be possible to separate the two by setting extern as everything not 
       !already defined as a bond, angle, dihedral etc NMR restraint.
       call gb_external(crd, frc, pot_ene%restraint)

!!NEED TO REWORK LOGIC TO ALLOW RESTRAINTS TO BE CALLED EVEN IF PMEMD IS NOT CALLED
!!With below code commented-in, then restraints are called twice which is not good.
#if 0
    if (nmropt .ne. 0) then
       call nmr_calc(crd, frc, enmr, 6)
       pot_ene%restraint = pot_ene%restraint + enmr(1) + enmr(2) + enmr(3)
       pot_ene%total = pot_ene%total + enmr(1) + enmr(2) + enmr(3)
    end if
#endif    

! If not calling an External library
  endif

!#if 0
  if (ti_mode .eq. 0) then
#ifdef _OPENMP_
    call gbene_ptr(crd, frc, atm_gb_radii, atm_gb_fs, atm_qterm, atm_iac, &
                typ_ico, atm_numex, gbl_natex, atm_cnt, belly_atm_cnt, &
                pot_ene%gb, pot_ene%elec_tot, pot_ene%vdw_tot,         &
                pot_ene%surf, irespa)
#else
    call gb_ene(crd, frc, atm_gb_radii, atm_gb_fs, atm_qterm, atm_iac, &
                typ_ico, atm_numex, gbl_natex, atm_cnt, belly_atm_cnt, &
                pot_ene%gb, pot_ene%elec_tot, pot_ene%vdw_tot,         &
                pot_ene%surf, irespa)
#endif /*_OPENMP_*/
  else
    ti_mask_piece = 2
    if (ti_mode .eq. 1) then !have to call twice
    ! NOTE: By calling gb_ene twice, this gives us at most 50% efficiency.
    ! However, if we try to split the function (as with the
    ! other energy calculations), we'll have problems if the rborn of all of
    ! the atoms in the TI regions are not the same. This causes the reff
    ! of common atoms to be different in V0/V1, so we lose most of the speed
    ! advantage in the general case.

    ! we might be able to take advantage of this in the gas phase case since we
    ! don't use rborn in gas phase
#ifdef _OPENMP_
      call gbene_ptr(crd, frc, atm_gb_radii, atm_gb_fs, atm_qterm, atm_iac, &
                  typ_ico, atm_numex, gbl_natex, atm_cnt, belly_atm_cnt, &
                  ti_pot_gb(1, 1), ti_pot_gb(1, 2), ti_pot_gb(1, 3),     &
                  ti_pot_gb(1, 4), irespa)
#else
      call gb_ene(crd, frc, atm_gb_radii, atm_gb_fs, atm_qterm, atm_iac, &
                  typ_ico, atm_numex, gbl_natex, atm_cnt, belly_atm_cnt, &
                  ti_pot_gb(1, 1), ti_pot_gb(1, 2), ti_pot_gb(1, 3),     &
                  ti_pot_gb(1, 4), irespa)
#endif /*_OPENMP_*/
    else !ti_mode > 1 i.e. scti
      gb_sc_eel = 0.d0
      gb_sc_evdw = 0.d0
!for now, OPENMP, ifsc gt 0, and igb gt 0 is not optimized
!not a big deal, since gas phase is on so few atoms anyways
!it gets like microseconds/day
      call gb_ene_sc(crd, frc, atm_gb_radii, atm_gb_fs, atm_qterm, atm_iac, &
                  typ_ico, atm_numex, gbl_natex, atm_cnt, belly_atm_cnt, &
                  ti_pot_gb(1, 1), ti_pot_gb(1, 2), ti_pot_gb(1, 3),     &
                  ti_pot_gb(1, 4), pot_ene%dvdl, gb_sc_eel, gb_sc_evdw, &
                  irespa)
!easier and hopefully more transparent to add sc energies here
!rather than pass than back out to runmd and write them there
!since gb other than igb6 can't really ever support softcore anyways
      ti_ene(1,si_vdw_ene) = ti_ene(1,si_vdw_ene) + gb_sc_evdw
      ti_ene(1,si_elect_ene) = ti_ene(1,si_elect_ene) + gb_sc_eel
      gb_sc_eel = 0.d0
      gb_sc_evdw = 0.d0
    end if
    ti_nb_frc(1,:,:) = frc(:,:) * ti_weights(1)
    frc(:,:) = 0.d0
    ti_mask_piece = 1
    if (ti_mode .eq. 1) then !have to call twice
#ifdef _OPENMP_
      call gbene_ptr(crd, frc, atm_gb_radii, atm_gb_fs, atm_qterm, atm_iac, &
                  typ_ico, atm_numex, gbl_natex, atm_cnt, belly_atm_cnt, &
                  ti_pot_gb(2, 1), ti_pot_gb(2, 2), ti_pot_gb(2, 3),     &
                  ti_pot_gb(2, 4), irespa)
#else
      call gb_ene(crd, frc, atm_gb_radii, atm_gb_fs, atm_qterm, atm_iac, &
                  typ_ico, atm_numex, gbl_natex, atm_cnt, belly_atm_cnt, &
                  ti_pot_gb(2, 1), ti_pot_gb(2, 2), ti_pot_gb(2, 3),     &
                  ti_pot_gb(2, 4), irespa)
#endif /*_OPENMP_*/
    else !ti_mode > 1 i.e. scti
      call gb_ene_sc(crd, frc, atm_gb_radii, atm_gb_fs, atm_qterm, atm_iac, &
                  typ_ico, atm_numex, gbl_natex, atm_cnt, belly_atm_cnt, &
                  ti_pot_gb(2, 1), ti_pot_gb(2, 2), ti_pot_gb(2, 3),     &
                  ti_pot_gb(2, 4), pot_ene%dvdl, gb_sc_eel, gb_sc_evdw, &
                  irespa)
      ti_ene(2,si_vdw_ene) = ti_ene(2,si_vdw_ene) + gb_sc_evdw
      ti_ene(2,si_elect_ene) = ti_ene(2,si_elect_ene) + gb_sc_eel
    end if
    frc(:,:) = ti_nb_frc(1,:,:) + frc(:,:) * ti_weights(2)
    call ti_update_ene_all(ti_pot_gb(:,1), si_hbond_ene, pot_ene%gb)
    call ti_update_ene_all(ti_pot_gb(:,2), si_elect_ene, pot_ene%elec_tot)
    call ti_update_ene_all(ti_pot_gb(:,3), si_vdw_ene,   pot_ene%vdw_tot)
    call ti_update_ene_all(ti_pot_gb(:,4), si_surf_ene,  pot_ene%surf)
  end if
       
  !PHMD
  if (iphmd .eq. 1) then
    call runphmd(crd,frc,atm_iac,typ_ico,atm_numex,gbl_natex,belly_atm_cnt)
  end if

! Calculate the 1-4 vdw and electrostatics contributions:

  if (charmm_active) then
#ifdef _OPENMP_
    call get_nb14_energy_gb(atm_cnt,atm_qterm, crd, frc, atm_iac, typ_ico, &
                         gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt, &
                         pot_ene%elec_14, pot_ene%vdw_14)
#else
    call get_nb14_energy(atm_qterm, crd, frc, atm_iac, typ_ico, &
                         gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt, &
                         pot_ene%elec_14, pot_ene%vdw_14)
#endif

  else
#ifdef _OPENMP_
    call get_nb14_energy_gb(atm_cnt,atm_qterm, crd, frc, atm_iac, typ_ico, &
                         gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt, &
                         pot_ene%elec_14, pot_ene%vdw_14)
#else
    call get_nb14_energy(atm_qterm, crd, frc, atm_iac, typ_ico, &
                         gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt, &
                         pot_ene%elec_14, pot_ene%vdw_14)
#endif
  end if

! For constant pH, get the EEL/VDW energies for the proposed state. This allows
! for easy extension to VDW modifications

  if (on_cpstep) then

    cnstph_frc(:,:) = 0.d0
    cph_egb         = 0.d0
    cph_eel         = 0.d0
    cph_eel14       = 0.d0
    cph_vdw         = 0.d0
    cph_vdw14       = 0.d0
    cph_surf        = 0.d0

    ! The last .true. argument forces gb_ene to skip calculating the effective
    ! radii, thereby saving time (since they will be unchanged from the previous
    ! gb_ene call)

#ifdef _OPENMP_
    call gbene_ptr(crd, cnstph_frc, atm_gb_radii, atm_gb_fs, proposed_qterm, &
                atm_iac, typ_ico, atm_numex, gbl_natex, atm_cnt,          &
                belly_atm_cnt, cph_egb, cph_eel, cph_vdw, cph_surf, irespa, &
                .true.)
#else
    call gb_ene(crd, cnstph_frc, atm_gb_radii, atm_gb_fs, proposed_qterm, &
                atm_iac, typ_ico, atm_numex, gbl_natex, atm_cnt,          &
                belly_atm_cnt, cph_egb, cph_eel, cph_vdw, cph_surf,       &
                irespa, .true.)
#endif /*_OPENMP_*/

    if (charmm_active) then
#ifdef _OPENMP_
      call get_nb14_energy_gb(atm_cnt,proposed_qterm, crd, cnstph_frc, atm_iac, typ_ico, &
                           gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt,      &
                           cph_eel14, cph_vdw14)
#else
      call get_nb14_energy(proposed_qterm, crd, cnstph_frc, atm_iac, typ_ico, &
                           gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt,      &
                           cph_eel14, cph_vdw14)
#endif /*_OPENMP_*/
    else
#ifdef _OPENMP_
      call get_nb14_energy_gb(atm_cnt,proposed_qterm, crd, cnstph_frc, atm_iac, typ_ico, &
                           gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt,          &
                           cph_eel14, cph_vdw14)
#else
      call get_nb14_energy(proposed_qterm, crd, cnstph_frc, atm_iac, typ_ico, &
                           gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt,          &
                           cph_eel14, cph_vdw14)
#endif /*_OPENMP_*/
    end if

    ! Now compute dvdl = proposed state - current state

    pot_ene%dvdl = cph_egb + cph_eel + cph_vdw + cph_surf + &
                   cph_eel14 + cph_vdw14 - &
                   pot_ene%gb - pot_ene%elec_tot - pot_ene%vdw_tot - &
                   pot_ene%surf - pot_ene%elec_14 - pot_ene%vdw_14

  end if ! on_cpstep

! For constant Redox potential, get the EEL/VDW energies for the proposed state. This allows
! for easy extension to VDW modifications

  if (on_cestep) then

    cnste_frc(:,:) = 0.d0
    ce_egb         = 0.d0
    ce_eel         = 0.d0
    ce_eel14       = 0.d0
    ce_vdw         = 0.d0
    ce_vdw14       = 0.d0
    ce_surf        = 0.d0

    ! The last .true. argument forces gb_ene to skip calculating the effective
    ! radii, thereby saving time (since they will be unchanged from the previous
    ! gb_ene call)

#ifdef _OPENMP_
    call gbene_ptr(crd, cnste_frc, atm_gb_radii, atm_gb_fs, proposed_qterm, &
                atm_iac, typ_ico, atm_numex, gbl_natex, atm_cnt,          &
                belly_atm_cnt, ce_egb, ce_eel, ce_vdw, ce_surf, irespa, &
                .true.)
#else
    call gb_ene(crd, cnste_frc, atm_gb_radii, atm_gb_fs, proposed_qterm, &
                atm_iac, typ_ico, atm_numex, gbl_natex, atm_cnt,          &
                belly_atm_cnt, ce_egb, ce_eel, ce_vdw, ce_surf, irespa, &
                .true.)
#endif /*_OPENMP_*/

    if (charmm_active) then
#ifdef _OPENMP_
      call get_nb14_energy_gb(atm_cnt,proposed_qterm, crd, cnste_frc, atm_iac, typ_ico, &
                           gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt,      &
                           ce_eel14, ce_vdw14)
#else
      call get_nb14_energy(proposed_qterm, crd, cnste_frc, atm_iac, typ_ico, &
                           gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt,      &
                           ce_eel14, ce_vdw14)
#endif /*_OPENMP_*/
    else
#ifdef _OPENMP_
      call get_nb14_energy_gb(atm_cnt,proposed_qterm, crd, cnste_frc, atm_iac, typ_ico, &
                           gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt,          &
                           ce_eel14, ce_vdw14)
#else
      call get_nb14_energy(proposed_qterm, crd, cnste_frc, atm_iac, typ_ico, &
                           gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt,          &
                           ce_eel14, ce_vdw14)
#endif /*_OPENMP_*/
    end if

    ! Now compute dvdl = proposed state - current state

    pot_ene%dvdl = ce_egb + ce_eel + ce_vdw + ce_surf + &
                   ce_eel14 + ce_vdw14 - &
                   pot_ene%gb - pot_ene%elec_tot - pot_ene%vdw_tot - &
                   pot_ene%surf - pot_ene%elec_14 - pot_ene%vdw_14

  end if ! on_cestep

  call update_time(nonbond_time)

! GaMD
  if (igamd.eq.4 .or. igamd.eq.5) then
    nb_frc_gamd = frc
  end if

  ! Calculate the other contributions:

  call gb_bonded_force(crd, frc, pot_ene)

  ! --- calculate the EMAP constraint energy ---

  if(iemap>0) then   ! (EMAP restraints)
    call emapforce(natom,enemap,atm_mass,crd,frc )
    pot_ene%emap = enemap
    pot_ene%restraint = pot_ene%restraint + enemap
  end if

  ! Sum up total potential energy for this task:

  pot_ene%total = pot_ene%vdw_tot + &
                  pot_ene%elec_tot + &
                  pot_ene%gb + &
                  pot_ene%surf + &
                  pot_ene%bond + &
                  pot_ene%angle + &
                  pot_ene%dihedral + &
                  pot_ene%vdw_14 + &
                  pot_ene%elec_14 + &
                  pot_ene%restraint + &
                  pot_ene%angle_ub + &
                  pot_ene%imp + &
                  pot_ene%cmap


#ifdef MPI
! Distribute energies to all processes.

  call gb_distribute_enes(pot_ene)

  if (ti_mode .ne. 0) then
    call ti_dist_enes_virs_netfrcs(.true.,.false.)
  end if

#endif /* MPI */

  ! Orient the extra points forces
  if (numextra .gt. 0 .and. frameon .ne. 0) then
    call orient_frc(crd, frc, vir_placeholder, ep_frames, ep_lcl_crd, gbl_frame_cnt)
  end if
  
 if(iamd.gt.1)then
 ! Calculate the boosting weight for dihedrals amd
   call calculate_amd_dih_weights(atm_cnt,pot_ene%dihedral,frc,crd)
 end if

 if((igamd.eq.2).or.(igamd.eq.3).or.(igamd.eq.5))then
 ! Calculate the boosting weight for dihedrals gamd
   call calculate_gamd_dih_weights(atm_cnt,pot_ene%dihedral,frc,crd)
 end if

#ifdef MPI
! Distribute forces to atom owners; this is a reduce sum operation:

  call gb_frcs_distrib(atm_cnt, frc, dbl_mpi_recv_buf)

  ! Clear the unused forces
  if (ti_mode .ne. 0) then
    do i = 1, atm_cnt
      if (gbl_atm_owner_map(i) .ne. mytaskid) frc(:,i) = 0.d0
    end do
  end if
#endif /* MPI */

  if (charmm_active .and. do_charmm_dump_gold == 1) then
#ifdef MPI
    call gb_mpi_gathervec(atm_cnt, frc)
    if (master) then
      call gb_charmm_dump_gold(atm_cnt,frc,pot_ene)
      write(mdout, '(a)') 'charmm_gold() completed. Exiting'
    end if
#else
    call gb_charmm_dump_gold(atm_cnt,frc,pot_ene)
    write(mdout, '(a)') 'charmm_gold() completed. Exiting'
#endif
    call mexit(6, 0)
  end if

! Calculate the NMR restraint energy contributions, if requested.

  if (nmropt .ne. 0) then
    call nmr_calc(crd, frc, enmr, 6)
    pot_ene%restraint = pot_ene%restraint + enmr(1) + enmr(2) + enmr(3)
    pot_ene%total = pot_ene%total + enmr(1) + enmr(2) + enmr(3)
  end if

  ! Built-in X-ray target function and gradient
  if (ti_mode .ne. 0) then
    call ti_calc_dvdl
  end if

!AMD DUAL BOOST CALC START
  if(iamd.gt.0)then
!calculate totboost and apply weight to frc. frc=frc*fwgt
    call calculate_amd_total_weights_gb(atm_cnt,pot_ene%total,pot_ene%dihedral,&
                                        pot_ene%amd_boost,frc,crd)
    ! Update total energy
    pot_ene%total = pot_ene%total + pot_ene%amd_boost
  end if

!GaMD DUAL BOOST CALC START
  if(igamd.gt.0)then
!calculate totboost and apply weight to frc. frc=frc*fwgt
     if(igamd.eq.1 .or. igamd.eq.2 .or. igamd.eq.3) then
       call calculate_gamd_total_weights_gb(atm_cnt, &
                       pot_ene%total,pot_ene%dihedral,pot_ene%gamd_boost,frc,crd)
       ! Update total energy
       pot_ene%total = pot_ene%total + pot_ene%gamd_boost
     else if(igamd.eq.4 .or. igamd.eq.5 ) then
       pot_ene_nb = pot_ene%elec_tot + pot_ene%vdw_tot + pot_ene%gb + pot_ene%vdw_14 + pot_ene%elec_14
       nb_frc_tmp = nb_frc_gamd
       call calculate_gamd_nb_weights_gb(atm_cnt, pot_ene_nb,pot_ene%dihedral,pot_ene%gamd_boost,nb_frc_gamd,crd)
       ! Update total force
       frc = frc - nb_frc_tmp + nb_frc_gamd
       ! Update total energy
       pot_ene%total = pot_ene%total + pot_ene%gamd_boost
     else if( (igamd.eq.6 .or. igamd.eq.7) .and. (ti_mode.ne.0)) then
       pot_ene_lig = ti_ene(1,si_pot_ene)
       pot_ene_dih = pot_ene%dihedral - ti_ene(1,si_dihedral_ene)
       call calculate_gamd_ti_region_weights_gb(atm_cnt,pot_ene_lig, pot_ene_dih,pot_ene%gamd_boost,fwgt,frc,crd,ti_sc_lst)
       ! Update total energy
       ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
     else if( (igamd.eq.8 .or. igamd.eq.9) .and. (ti_mode.ne.0)) then
       pot_ene_lig = ti_ene(1,si_vdw_ene) + ti_ene(1,si_elect_ene) - ti_ene(1,si_dvdl)
       pot_ene_dih = pot_ene%dihedral - ti_ene(1,si_dihedral_ene)
       call calculate_gamd_ti_region_weights_gb(atm_cnt,pot_ene_lig, pot_ene_dih,pot_ene%gamd_boost,fwgt,frc,crd,ti_sc_lst)
       ! Update total energy
       ti_ene(1,si_pot_ene) = ti_ene(1,si_pot_ene) + pot_ene%gamd_boost
     end if
   end if

!scaledMD
  if(scaledMD.gt.0)then
!apply scale to frc. frc=frc*scaledMD_lambda
    call scaledMD_scale_frc_gb(atm_cnt,pot_ene%total,frc,crd)
  end if

  ! If belly is on then set the belly atom forces to zero:

  if (ibelly .gt. 0) call bellyf(atm_cnt, atm_igroup, frc)

  if (infe.gt.0) then
      bias_frc(:,:) = 0.d0
      call nfe_on_force(crd,bias_frc,pot_ene%nfe)

      do i = 1, atm_cnt
#ifdef MPI
        if (gbl_atm_owner_map(i) .eq. mytaskid) then
#endif
          frc(1, i) = frc(1,i)+bias_frc(1, i)
          frc(2, i) = frc(2,i)+bias_frc(2, i)
          frc(3, i) = frc(3,i)+bias_frc(3, i)
#ifdef MPI
        end if
#endif
      end do

      pot_ene%total = pot_ene%total + pot_ene%nfe
      pot_ene%restraint = pot_ene%restraint + pot_ene%nfe
  end if

  if (igamd.eq.4 .or. igamd.eq.5) then
    deallocate(nb_frc_gamd,nb_frc_tmp)
  end if

#ifdef MPI
  if (ineb>0) then
    call gb_mpi_gathervec(atm_cnt, crd)
    call gb_mpi_gathervec(atm_cnt, frc)

    if(mytaskid.eq.0) then
      call full_neb_forces(atm_mass, crd, frc, pot_ene%total, fitmask, rmsmask, irespa)
    end if
    call mpi_bcast(neb_force, 3*nattgtrms, MPI_DOUBLE_PRECISION, 0, pmemd_comm, err_code_mpi)

    do jatm=1,nattgtrms
      iatm=rmsmask(jatm)
      if(iatm .eq. 0) cycle
      if (gbl_atm_owner_map(iatm) .eq. mytaskid) then
        frc(1,iatm)= frc(1,iatm)+neb_force(1,jatm)
        frc(2,iatm)= frc(2,iatm)+neb_force(2,jatm)
        frc(3,iatm)= frc(3,iatm)+neb_force(3,jatm)
      end if
    end do

    if (beadid==1 .or. beadid==neb_nbead) then
      do i=1, atm_cnt
        if (gbl_atm_owner_map(i) .eq. mytaskid) then
          frc(:,i)=0.d0
        endif
      enddo
    end if
  end if
#endif /* MPI */

! Ending if of External library
!end if

#endif /* CUDA */
  return

end subroutine gb_force

!*******************************************************************************
!
! Subroutine:  gb_bonded_force
!
! Description: <TBS>
!
!*******************************************************************************

subroutine gb_bonded_force(crd, frc, pot_ene)

  use angles_mod
  use angles_ub_mod
  use bonds_mod
  use constraints_mod
  use dihedrals_mod
  use dihedrals_imp_mod
  use cmap_mod
  use dynamics_dat_mod
  use mdin_ctrl_dat_mod
  use nmr_calls_mod
  use timers_mod

  implicit none

! Formal arguments:

  double precision              :: crd(3, *)
  double precision              :: frc(3, *)
  type(gb_pot_ene_rec)          :: pot_ene

! Local variables:

  ! These energy variables are temporaries, for summing. DON'T use otherwise!

  double precision              :: bond_ene
  double precision              :: angle_ene
  double precision              :: dihedral_ene
  double precision              :: ub_ene
  double precision              :: dihedral_imp_ene
  double precision              :: cmap_ene

  double precision              :: molvir(3, 3)         ! dummy argument
  double precision              :: e14_vir(3, 3)        ! dummy argument

! Bond energy contribution:

! The ebdev/eadev stuff currently only is output under nmr_calls for non-mpi
! code, so we basically drop it here under mpi.

#ifdef MPI
#else
  ebdev = 0.d0
#endif
  if (ntf .le. 1) then
    if (cit_nbonh .gt. 0) then
#ifdef _OPENMP_
      call get_bond_energy_gb(cit_nbonh, cit_h_bond, crd, frc, bond_ene)
#else
      call get_bond_energy(cit_nbonh, cit_h_bond, crd, frc, bond_ene)
#endif /*_OPENMP_*/
      pot_ene%bond = bond_ene
    end if
  end if

  if (ntf .le. 2) then
    if (cit_nbona .gt. 0) then
#ifdef _OPENMP_
      call get_bond_energy_gb(cit_nbona, cit_a_bond, crd, frc, bond_ene)
#else
      call get_bond_energy(cit_nbona, cit_a_bond, crd, frc, bond_ene)
#endif /*_OPENMP_*/
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
#ifdef _OPENMP_
       call get_angle_energy_gb(cit_ntheth, cit_angle, crd, frc, angle_ene)
#else
        call get_angle_energy(cit_ntheth, cit_angle, crd, frc, angle_ene)
#endif /*_OPENMP_*/
      pot_ene%angle = angle_ene
    end if
  end if

  if (ntf .le. 4) then
    if (cit_ntheta .gt. 0) then
#ifdef _OPENMP_
         call get_angle_energy_gb(cit_ntheta, cit_angle(cit_ntheth+1), &
                            crd, frc, angle_ene)
#else
         call get_angle_energy(cit_ntheta, cit_angle(cit_ntheth+1), &
                            crd, frc, angle_ene)
#endif /*_OPENMP_*/
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

!CMAP
  if (cit_cmap_term_count .gt. 0) then
    call get_cmap_energy(cit_cmap_term_count, cit_cmap, crd, frc, cmap_ene)
    pot_ene%cmap = cmap_ene
  end if

  ! Dihedral energy contribution:

  if(iamd.gt.1)then
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
  else if((igamd.eq.2).or.(igamd.eq.3).or.(igamd.eq.5))then
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
#ifdef _OPENMP_
        call get_dihed_energy_gb(cit_nphih, cit_dihed, crd, frc, dihedral_ene)
#else
        call get_dihed_energy(cit_nphih, cit_dihed, crd, frc, dihedral_ene)
#endif /*_OPENMP_*/
        pot_ene%dihedral = dihedral_ene
      end if
    end if

    if (ntf .le. 6) then
      if (cit_nphia .gt. 0) then
#ifdef _OPENMP_
        call get_dihed_energy_gb(cit_nphia, cit_dihed(cit_nphih + 1), crd, frc, &
          dihedral_ene)
#else
        call get_dihed_energy(cit_nphia, cit_dihed(cit_nphih + 1), crd, frc, &
          dihedral_ene)
#endif /*_OPENMP_*/
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

end subroutine gb_bonded_force

!*******************************************************************************
!
! Subroutine:  gb_cph_ene
!
! Description: Gets CpHMD dvdl term
!
!*******************************************************************************

subroutine gb_cph_ene(atm_cnt, crd, qterm, dvdl, skipradii)

  use charmm_mod, only          : charmm_active
  use constantph_dat_mod, only  : cnstph_frc
  use constraints_mod, only     : belly_atm_cnt
  use gb_ene_mod, only          : gb_ene
  use extra_pnts_nb14_mod, only : cit_nb14, cit_nb14_cnt, get_nb14_energy
  use parallel_dat_mod
  use prmtop_dat_mod, only      : atm_gb_radii, atm_gb_fs, atm_iac, typ_ico, &
                                  atm_numex, gbl_natex, gbl_cn114, gbl_cn214, &
                                  gbl_cn1, gbl_cn2
  use mdin_ctrl_dat_mod, only   : ineb
  
  implicit none

  ! Passed arguments

  integer, intent(in)             :: atm_cnt
  double precision, intent(in)    :: crd(3, atm_cnt)
  double precision, intent(in)    :: qterm(atm_cnt)
  double precision, intent(out)   :: dvdl
  logical, intent(in)             :: skipradii

  ! Local arguments

  double precision :: egb = 0.d0
  double precision :: eel = 0.d0
  double precision :: dvdl_buf = 0.d0
  double precision :: surf = 0.d0
  double precision :: vdw = 0.d0
  double precision :: eel14 = 0.d0
  double precision :: vdw14 = 0.d0
  type(gb_pot_ene_rec)          :: pot_ene
  double precision              :: enmr(3)

  ! Compute the energy of the current state
#ifdef CUDA
  enmr(:) = 0.d0
  pot_ene = null_gb_pot_ene_rec
  call gpu_gb_ene(pot_ene, enmr, ineb)
  egb = pot_ene%gb
  eel = pot_ene%elec_tot
  vdw = pot_ene%vdw_tot
  surf = pot_ene%surf
  eel14 = pot_ene%elec_14
  vdw14 = pot_ene%vdw_14
#else
!remove dvdl_buf or else make constant pH accessible with this.
  call gb_ene(crd, cnstph_frc, atm_gb_radii, atm_gb_fs, qterm,   &
              atm_iac, typ_ico, atm_numex, gbl_natex, atm_cnt,   &
              belly_atm_cnt, egb, eel, vdw, surf, 0, skipradii)
#endif

  if (charmm_active) then
    call get_nb14_energy(qterm, crd, cnstph_frc, atm_iac, typ_ico, &
                         gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt, &
                         eel14, vdw14)
  else
    call get_nb14_energy(qterm, crd, cnstph_frc, atm_iac, typ_ico, &
                         gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt,     &
                         eel14, vdw14)
  end if

  ! Store the constant pH-relevant portion of the current state energy

  dvdl = egb + eel + vdw + surf + eel14 + vdw14

#ifdef MPI
  ! We need to reduce our total energy in parallel
  call mpi_allreduce(dvdl, dvdl_buf, 1, mpi_double_precision, mpi_sum, &
                     pmemd_comm, err_code_mpi)

  dvdl = dvdl_buf
#endif

  return

end subroutine gb_cph_ene

!*******************************************************************************
!
! Subroutine:  gb_ce_ene
!
! Description: Gets CEMD dvdl term
!
!*******************************************************************************

subroutine gb_ce_ene(atm_cnt, crd, qterm, dvdl, skipradii)

  use charmm_mod, only          : charmm_active
  use constante_dat_mod, only  : cnste_frc
  use constraints_mod, only     : belly_atm_cnt
  use gb_ene_mod, only          : gb_ene
  use extra_pnts_nb14_mod, only : cit_nb14, cit_nb14_cnt, get_nb14_energy
  use parallel_dat_mod
  use prmtop_dat_mod, only      : atm_gb_radii, atm_gb_fs, atm_iac, typ_ico, &
                                  atm_numex, gbl_natex, gbl_cn114, gbl_cn214, &
                                  gbl_cn1, gbl_cn2
  use mdin_ctrl_dat_mod, only   : ineb
  
  implicit none

  ! Passed arguments

  integer, intent(in)             :: atm_cnt
  double precision, intent(in)    :: crd(3, atm_cnt)
  double precision, intent(in)    :: qterm(atm_cnt)
  double precision, intent(out)   :: dvdl
  logical, intent(in)             :: skipradii

  ! Local arguments

  double precision :: egb = 0.d0
  double precision :: eel = 0.d0
  double precision :: dvdl_buf = 0.d0
  double precision :: surf = 0.d0
  double precision :: vdw = 0.d0
  double precision :: eel14 = 0.d0
  double precision :: vdw14 = 0.d0
  type(gb_pot_ene_rec)          :: pot_ene
  double precision              :: enmr(3)

  ! Compute the energy of the current state
#ifdef CUDA
  enmr(:) = 0.d0
  pot_ene = null_gb_pot_ene_rec
  call gpu_gb_ene(pot_ene, enmr, ineb)
  egb = pot_ene%gb
  eel = pot_ene%elec_tot
  vdw = pot_ene%vdw_tot
  surf = pot_ene%surf
  eel14 = pot_ene%elec_14
  vdw14 = pot_ene%vdw_14
#else
  call gb_ene(crd, cnste_frc, atm_gb_radii, atm_gb_fs, qterm,   &
              atm_iac, typ_ico, atm_numex, gbl_natex, atm_cnt,   &
              belly_atm_cnt, egb, eel, vdw, surf, 0, skipradii)
#endif

  if (charmm_active) then
    call get_nb14_energy(qterm, crd, cnste_frc, atm_iac, typ_ico, &
                         gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt, &
                         eel14, vdw14)
  else
    call get_nb14_energy(qterm, crd, cnste_frc, atm_iac, typ_ico, &
                         gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt,     &
                         eel14, vdw14)
  end if

  ! Store the constant E-relevant portion of the current state energy

  dvdl = egb + eel + vdw + surf + eel14 + vdw14

#ifdef MPI
  ! We need to reduce our total energy in parallel
  call mpi_allreduce(dvdl, dvdl_buf, 1, mpi_double_precision, mpi_sum, &
                     pmemd_comm, err_code_mpi)

  dvdl = dvdl_buf
#endif

!  write (0, '(5(a,f15.8,x))') &
!     'EGB = ', egb, 'EEL = ', eel, 'VDW = ', vdw, 'EEL 1-4 = ', eel14, &
!     'VDW 1-4 = ', vdw14

  return

end subroutine gb_ce_ene

!*******************************************************************************
!
! Subroutine:  gb_hybridsolvent_remd_ene
!
! Description: <TBS>
!
!*******************************************************************************

subroutine gb_hybridsolvent_remd_ene(atm_cnt, crd, frc, pot_ene, irespa, need_pot_ene)

  use charmm_mod, only          : charmm_active
  use constraints_mod, only     : belly_atm_cnt
  use gb_ene_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use extra_pnts_nb14_mod
  use gbsa_mod
  use mdin_ctrl_dat_mod
  use timers_mod

  implicit none

! Formal arguments:

  integer                       :: atm_cnt
  double precision              :: crd(3, atm_cnt)
  double precision              :: frc(3, atm_cnt)
  type(gb_pot_ene_rec)          :: pot_ene
  integer, intent(in)           :: irespa
  logical, intent(in)           :: need_pot_ene

! Local variables:

  double precision              :: enmr(3)
!  double precision              :: egb = 0.d0
!  double precision              :: eel = 0.d0
!  double precision              :: vdw = 0.d0
!  double precision              :: surf = 0.d0
!  double precision              :: eel14 = 0.d0
!  double precision              :: vdw14 = 0.d0
!  double precision              :: bond = 0.d0
!  double precision              :: angle = 0.d0
!  double precision              :: dihedral = 0.d0
!  double precision              :: restraint = 0.d0
!  double precision              :: angle_ub = 0.d0
!  double precision              :: imp = 0.d0
!  double precision              :: cmap = 0.d0

!#ifdef CUDA
!  ! Update any weight modifications and zero-out the potential energy record
!  enmr(:) = 0
!  pot_ene = null_gb_pot_ene_rec !Zero's entire structure
!  call gpu_upload_charges_gb_cph(atm_qterm)
!  call gpu_gb_ene(pot_ene, enmr, ineb)
!  write(mdout,'(a,f13.4)') 'pot_ene%total=', pot_ene%total
!  write(mdout,'(a,f13.4)') 'pot_ene%egb=', pot_ene%gb
!  write(mdout,'(a,f13.4)') 'pot_ene%eel=', pot_ene%elec_tot
!  write(mdout,'(a,f13.4)') 'pot_ene%vdw=', pot_ene%vdw_tot
!  write(mdout,'(a,f13.4)') 'pot_ene%surf=', pot_ene%surf
!  write(mdout,'(a,f13.4)') 'pot_ene%eel14=', pot_ene%elec_14
!  write(mdout,'(a,f13.4)') 'pot_ene%vdw14=', pot_ene%vdw_14
!  write(mdout,'(a,f13.4)') 'pot_ene%bond=', pot_ene%bond
!  write(mdout,'(a,f13.4)') 'pot_ene%angle=', pot_ene%angle
!  write(mdout,'(a,f13.4)') 'pot_ene%dihedral=', pot_ene%dihedral
!
!! BEGIN GBSA implementation to work with the GPU code.
!! It runs on the CPU for now, needs to write a kernel for this.
!! Regular CPU call is located in gb_ene
!! pmemd.cuda only takes nrespa=1 for now, so GBSA contribution must ALWAYS
!! be calculated (regardless of need_pot_ene)
!  if (gbsa .eq. 1) then
!    frc=0.d0
!    gbsafrc=0.d0
!    call gpu_download_crd(crd)
!    call gbsa_ene(crd, gbsafrc, pot_ene%surf ,atm_cnt, jj, r2x, belly_atm_cnt)
!    pot_ene%total =  pot_ene%total +  pot_ene%surf
!  end if
!! END GBSA implementation to work with the GPU code.
!
!#else /* CUDA */

! Zero energies that are stack or call parameters:

  pot_ene = null_gb_pot_ene_rec

! Zero internal energies:

  enmr(:) = 0.d0

  call zero_time()
  call zero_gb_time()

  frc(:,:) = 0.d0

  call gb_ene(crd, frc, atm_gb_radii, atm_gb_fs, atm_qterm, atm_iac, &
              typ_ico, atm_numex, gbl_natex, atm_cnt, belly_atm_cnt, &
              pot_ene%gb, pot_ene%elec_tot, pot_ene%vdw_tot,         &
              pot_ene%surf, irespa, .false.)

  if (charmm_active) then
#ifdef _OPENMP_
     call get_nb14_energy_gb(atm_cnt,atm_qterm, crd, frc, atm_iac, typ_ico, &
                         gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt, &
                         pot_ene%elec_14, pot_ene%vdw_14)
#else
     call get_nb14_energy(atm_qterm, crd, frc, atm_iac, typ_ico, &
                         gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt, &
                         pot_ene%elec_14, pot_ene%vdw_14)
#endif

  else
#ifdef _OPENMP_
     call get_nb14_energy_gb(atm_cnt,atm_qterm, crd, frc, atm_iac, typ_ico, &
                         gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt, &
                         pot_ene%elec_14, pot_ene%vdw_14)
#else
     call get_nb14_energy(atm_qterm, crd, frc, atm_iac, typ_ico, &
                         gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt, &
                         pot_ene%elec_14, pot_ene%vdw_14)
#endif
  end if

  ! Calculate the other contributions:

  call gb_bonded_force(crd, frc, pot_ene)

  call update_time(nonbond_time)

  ! Sum up total potential energy for this task:

  pot_ene%total = pot_ene%vdw_tot + &
                  pot_ene%elec_tot + &
                  pot_ene%gb + &
                  pot_ene%surf + &
                  pot_ene%bond + &
                  pot_ene%angle + &
                  pot_ene%dihedral + &
                  pot_ene%vdw_14 + &
                  pot_ene%elec_14 + &
                  pot_ene%restraint + &
                  pot_ene%angle_ub + &
                  pot_ene%imp + &
                  pot_ene%cmap

#ifdef MPI
! Distribute energies to all processes.

  call gb_distribute_enes(pot_ene)

#endif /* MPI */

!#endif /* CUDA */
  return

end subroutine gb_hybridsolvent_remd_ene

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  gb_distribute_enes
!
! Description: <TBS>
!
!*******************************************************************************

subroutine gb_distribute_enes(pot_ene)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  type(gb_pot_ene_rec) :: pot_ene

! Local variables:

  type(gb_pot_ene_rec), save    :: dat_in, dat_out

  dat_in = pot_ene

  call mpi_allreduce(dat_in%total, dat_out%total, &
                     gb_pot_ene_rec_size, mpi_double_precision, &
                     mpi_sum, pmemd_comm, err_code_mpi)

  pot_ene = dat_out

  return

end subroutine gb_distribute_enes
#endif

end module gb_force_mod
