
!************************************************************************
!                              AMBER                                   **
!                                                                      **
!                        Copyright (c) 2008                            **
!                Regents of the University of California               **
!                       All Rights Reserved.                           **
!                                                                      **
!  This software provided pursuant to a license agreement containing   **
!  restrictions on its disclosure, duplication, and use. This software **
!  contains confidential and proprietary information, and may not be   **
!  extracted or distributed, in whole or in part, for any purpose      **
!  whatsoever, without the express written permission of the authors.  **
!  This notice, and the associated author list, must be attached to    **
!  all copies, or extracts, of this software. Any additional           **
!  restrictions set forth in the license agreement also apply to this  **
!  software.                                                           **
!************************************************************************

#ifdef EMIL

!+++++++++++++++++++++++++++++++++++++++++++++
!This module contains an interface and some state for an EMIL free energy calc.
!The meaningful work is done outside pmemd, in the EMIL library code  
!(kept in ../../../AmberTools/src/emil at time of writing).
!+++++++++++++++++++++++++++++++++++++++++++++
module emil_mod
   implicit none
   private

   !state variables for this module
   double precision, public, save    :: emil_lambdamol, emil_lambdarest
   double precision, public, save    :: emil_softcore_scale

   !!Stuff for Monte Carlo move state
    !using int and not bool for cross-platform stability
   integer, public, save             :: emil_mc_done
   
   !!flag to indicate if a new forces calc is needed.
   integer, public, save             :: emil_calc_AMBER
   logical, public, save             :: emil_save_pme, emil_save_gb

   !temporary storage for energies 
   double precision, public          :: emil_base_epot, emil_base_evdw, emil_base_eelec 

   !subroutines
   public  emil_init
   public  emil_step


 contains

subroutine emil_init (atm_cnt, nstep, beta, mass, crd, frc, vel)
  
  !import some state variables
  use mdin_ctrl_dat_mod,   only : nstlim, using_gb_potential, vdw_cutoff
  use prmtop_dat_mod,      only : atm_iac, atm_igraph, atm_nsp, atm_qterm, gbl_cn1, gbl_cn2, &
                                  gbl_labres, gbl_res_atms, nres, ntypes, typ_ico, &
                                  atm_numex, gbl_natex
  use mol_list_mod,        only : gbl_mol_cnt 
  use pbc_mod,             only : pbc_box
  use mdin_emil_dat_mod,   only : emil_paramfile, emil_logfile, emil_model_infile, &
                                  emil_model_outfile
  use ti_mod,              only : emil_sc_lcl
#ifdef MPI
  use parallel_dat_mod,    only : mytaskid, numtasks, pmemd_comm
#endif

  implicit none

  !Arguments
  integer, intent(in)               :: atm_cnt
  integer, intent(in)               :: nstep
  double precision, intent(in)      :: beta 
  double precision, intent(in)      :: mass(atm_cnt)
  double precision, intent(inout)   :: crd(3, atm_cnt), frc(3, atm_cnt), vel(3, atm_cnt)


  emil_lambdamol      = 0.0
  emil_lambdarest     = 0.0
  emil_softcore_scale = 0.0

  emil_save_pme       = .false.
  emil_save_gb        = .false.


  ! Setup emil - position wells, open outfiles, all that
#ifdef MPI
  call c_emil_init( atm_cnt, nres, gbl_mol_cnt, gbl_res_atms, atm_nsp, mass, atm_igraph, &
             gbl_labres, ntypes, vdw_cutoff, atm_iac, typ_ico, atm_qterm, &
             crd, frc, vel, pbc_box, gbl_cn1,  gbl_cn2, &
             emil_lambdamol, emil_lambdarest,  emil_softcore_scale, &
             beta,  mytaskid, numtasks, pmemd_comm, nstlim, &
             atm_numex, gbl_natex, & 
             trim(emil_paramfile)//char(0), &   !C strings need a null-terminator
             trim(emil_logfile)  //char(0), &
             trim(emil_model_infile)//char(0), &
             trim(emil_model_outfile)//char(0), emil_sc_lcl )
#else
  call c_emil_init( atm_cnt, nres, gbl_mol_cnt, gbl_res_atms, atm_nsp, mass, atm_igraph, &
             gbl_labres, ntypes, vdw_cutoff, atm_iac, typ_ico, atm_qterm, &
             crd, frc, vel, pbc_box, gbl_cn1,  gbl_cn2, &
             emil_lambdamol, emil_lambdarest,  emil_softcore_scale, beta,  0, 1, nstlim, &
             atm_numex, gbl_natex, &
             trim(emil_paramfile)//char(0), &   !C strings need a null-terminator
             trim(emil_logfile)  //char(0), &
             trim(emil_model_infile)//char(0) , &
             trim(emil_model_outfile)//char(0), emil_sc_lcl )
#endif


  return !End subroutine emil_init
  end subroutine emil_init

subroutine emil_step (atm_cnt,        &      !total atoms
                      nstep,          &      !step count
                      beta,           &      !1/kBT in AMBER units
                      mass,           &      !atomic masses
                      crd,            &      !coordinates (whole system on each node)
                      frc,            &      !forces      (whole system on each node)
                      vel,            &      !velocities     "
                      gb_pot_ene,     &      !potential energies datatype
                      pme_pot_ene,    &      ! "
                      softcore_dvdl,  &      !dvdl only if amber-side softcoring
                      si)                    !general state info
 
  !import derived datatypes
  use gb_force_mod,        only : gb_pot_ene_rec, null_gb_pot_ene_rec
  use pme_force_mod,       only : pme_pot_ene_rec, null_pme_pot_ene_rec
  use state_info_mod

  !import state
  use mdin_ctrl_dat_mod,   only : using_gb_potential, using_pme_potential
  use pbc_mod,             only : pbc_box
  use ti_mod,              only : emil_sc_lcl

#ifdef MPI
  use parallel_dat_mod,    only : mytaskid, numtasks, pmemd_comm, &
                                  my_atm_cnt, gbl_my_atm_lst
#endif

  implicit none

  !Arguments
  integer, intent(in)                  :: atm_cnt
  double precision, intent(in)         :: beta 
  integer, intent(in)                  :: nstep
  double precision, intent(in)         :: softcore_dvdl
  double precision, intent(in)         :: mass(atm_cnt)
  double precision, intent(inout)      :: crd(3, atm_cnt), frc(3, atm_cnt), vel(3, atm_cnt)
  double precision, intent(inout)      :: si(si_cnt)
  type(gb_pot_ene_rec),  intent(inout) :: gb_pot_ene
  type(pme_pot_ene_rec), intent(inout) :: pme_pot_ene
  

  !Local variables
  double precision      :: potential_energy, temp_dvdl

  !Where do we find the total energy of the system?
  if ( using_gb_potential )  then
         potential_energy = gb_pot_ene%total
  else if ( using_pme_potential ) then
         potential_energy = pme_pot_ene%total
  else 
         potential_energy = 0.0
  endif

  !make sure that the C-code does not get hold of a pointer
  !to AMBER's dvdl value.
  if ( emil_sc_lcl .ne. 0 ) then
         temp_dvdl = softcore_dvdl
  else
         temp_dvdl = 0.0
  endif

  !Try some MC moves
  call c_emil_mcattempt( crd, frc, vel, pbc_box, emil_base_epot, emil_mc_done, &
                           beta, emil_base_evdw, emil_base_eelec )

#ifdef MPI      
  !determine forces and energies due to the fancy emil restraining potentials
  call c_emil_forces_pme( crd, frc, pbc_box, temp_dvdl, potential_energy, &
                        0.0, beta, nstep, gbl_my_atm_lst, my_atm_cnt, &
                        emil_calc_AMBER )
#else
  !determine forces and energies due to the fancy emil restraining potentials
  call c_emil_forces_pme( crd, frc, pbc_box, temp_dvdl, potential_energy, &
                        0.0, beta, nstep, 0, atm_cnt, &
                        emil_calc_AMBER )
#endif

  !Is the AMBER Hamiltonian in play at all?  If not, can save time:
  ! This block should switch OFF whichever of pme or gb is ON
  if ( emil_calc_AMBER .eq. 0 ) then
     if ( using_pme_potential ) then
        using_pme_potential = .false.
        emil_save_pme       = .true.
        
        !reset the logging
        si(si_pot_ene)       = 0.d0
        si(si_vdw_ene)       = 0.d0
        si(si_elect_ene)     = 0.d0
        si(si_hbond_ene)     = 0.d0
        si(si_bond_ene)      = 0.d0
        si(si_angle_ene)     = 0.d0
        si(si_dihedral_ene)  = 0.d0
        si(si_vdw_14_ene)    = 0.d0
        si(si_elect_14_ene)  = 0.d0
        si(si_restraint_ene) = 0.d0
        si(si_pme_err_est)   = 0.d0
        si(si_angle_ub_ene)  = 0.d0
        si(si_dihedral_imp_ene) = 0.d0
        si(si_cmap_ene)      = 0.d0

     else if ( using_gb_potential ) then
        using_gb_potential  = .false.
        emil_save_gb        = .true.
   
       !reset the logging
        si(si_pot_ene)       = 0.d0
        si(si_vdw_ene)       = 0.d0
        si(si_elect_ene)     = 0.d0
        si(si_hbond_ene)     = 0.d0
        si(si_surf_ene)      = 0.d0
        si(si_bond_ene)      = 0.d0
        si(si_angle_ene)     = 0.d0
        si(si_dihedral_ene)  = 0.d0
        si(si_vdw_14_ene)    = 0.d0
        si(si_elect_14_ene)  = 0.d0
        si(si_restraint_ene) = 0.d0
        si(si_angle_ub_ene)  = 0.d0
        si(si_dihedral_imp_ene) = 0.d0
        si(si_cmap_ene)      = 0.d0
     end if
  else 
     !This block should switch ON whichever of them was previously OFF
     if ( emil_save_pme ) then
        using_pme_potential = .true.
        emil_save_pme       = .false.
     else if ( emil_save_gb ) then
        using_gb_potential  = .true.
        emil_save_gb        = .false.
     end if
  end if

  return !End subroutine emil_step
  end subroutine emil_step


end module emil_mod
#endif 
!===================================================================================================
