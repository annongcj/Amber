#include "copyright.i"

!*******************************************************************************
!
! Module: gb_alltasks_setup_mod
!
! Description:  Setup of data structures for uniprocessor code as well as
!               mpi master and slave processors.
!              
!*******************************************************************************

module gb_alltasks_setup_mod

use file_io_dat_mod

  implicit none

contains

!*******************************************************************************
!
! Subroutine:  gb_alltasks_setup
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine gb_alltasks_setup(num_ints, num_reals)

  use angles_mod
  use angles_ub_mod
  use cmap_mod
  use bonds_mod
  use constraints_mod
  use dihedrals_mod
  use dihedrals_imp_mod
  use dynamics_mod
  use dynamics_dat_mod
  use extra_pnts_nb14_mod
  use mol_list_mod
  use prfs_mod
  use img_mod
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use gb_parallel_mod
  use parallel_mod
  use prfs_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
  use random_mod
  use charmm_mod, only : charmm_active
  use mdin_debugf_dat_mod
#ifdef CUDA
  use nmr_calls_mod
#endif

  implicit none

! Formal arguments:

! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed

#ifdef MPI
  integer               :: i
  integer               :: node
  integer               :: taskmap_idx
#endif
  integer               :: use_atm_map(natom)

  use_atm_map(:) = 0

  ! Molecule lists are not needed for GB:
  if  ( emil_do_calc .ne. 0 ) then !...but still useful for emil calc.
    call setup_molecule_lists(num_ints, num_reals, natom, nspm, atm_nsp, &
                              nbona, gbl_bond(bonda_idx), &
                              no_intermolecular_bonds)
  else
    num_ints = num_ints - size(atm_nsp)
    deallocate(atm_nsp)
    nspm = 0
  end if

  ! Initialize dynamics data.

  call dynamics_dat_init(natom, ntp, imin, atm_mass, atm_crd, &
                         num_ints, num_reals)

#ifdef MPI
  call setup_prfs(num_ints, num_reals, natom, nres, gbl_res_atms, &
                  nbonh, gbl_bond, gbl_frame_cnt, ep_frames)

  ! For the first GB implementation, we attempt to evenly divide the
  ! atom workload on prf boundaries.  We will keep track of atom
  ! ownership, but will keep all coordinates updated in all processes,
  ! assuming that large cutoff size and small atom count make it somewhere
  ! between unnecessary and perhaps even detrimental to do otherwise
  ! (in other words, attaining processor spatial locality is not practical
  ! with a small irregularly shaped population of atoms, given the large
  ! nonbonded cutoffs typical of GB).

  call gb_parallel_setup(num_ints, num_reals)

  call bcast_mdin_debugf_dat()

  ! Bond-angle-dihedral ownership needs to be established.  We use the
  ! use_atm_map here; in reality the atom usage info will not be kept for
  ! GB though.

  call bonds_setup(num_ints, num_reals, use_atm_map)
  call angles_setup(num_ints, num_reals, use_atm_map)
  call dihedrals_setup(num_ints, num_reals, use_atm_map)
#ifdef CUDA 
  if (nmropt .ne. 0) then
    call cuda_nmr_setup()         
  end if
#endif  
  if (charmm_active) then
    call angles_ub_setup(num_ints, num_reals, use_atm_map)
    call dihedrals_imp_setup(num_ints, num_reals, use_atm_map)
  endif
  ! if CHARMM or CMAP is active
  if (charmm_active .or. cmap_term_count > 0) then
    call cmap_setup(num_ints, num_reals, use_atm_map)
  endif
  call nb14_setup(num_ints, num_reals, use_atm_map)

#else /* begin non-MPI code */

  ! Uniprocessor values:

  my_mol_cnt = gbl_mol_cnt

  call bonds_setup(num_ints, num_reals, use_atm_map)
  call angles_setup(num_ints, num_reals, use_atm_map)
  call dihedrals_setup(num_ints, num_reals, use_atm_map)
#ifdef CUDA 
  if (nmropt .ne. 0) then
    call cuda_nmr_setup()         
  end if
#endif  
  if (charmm_active) then
    call angles_ub_setup(num_ints, num_reals, use_atm_map)
    call dihedrals_imp_setup(num_ints, num_reals, use_atm_map)
  endif
  ! if CHARMM or CMAP is active
  if (charmm_active .or. cmap_term_count > 0) then
    call cmap_setup(num_ints, num_reals, use_atm_map)
  endif
  call nb14_setup(num_ints, num_reals, use_atm_map)

  ! gbl_bond is still needed for shake setup and resetup, but gbl_angle and
  ! gbl_dihed can be deallocated .

  num_ints = num_ints - size(gbl_angle) * angle_rec_ints
  num_ints = num_ints - size(gbl_dihed) * dihed_rec_ints
  if (charmm_active) then
    num_ints = num_ints - size(gbl_angle_ub) * angle_ub_rec_ints
    num_ints = num_ints - size(gbl_dihed_imp) * dihed_imp_rec_ints
    num_ints = num_ints - size(gbl_cmap) * cmap_rec_ints

    deallocate(gbl_angle_ub, gbl_dihed_imp, gbl_cmap)
  ! if CMAP is active
  else if (cmap_term_count > 0) then
    num_ints = num_ints - size(gbl_cmap) * cmap_rec_ints

    deallocate(gbl_cmap)
  endif

  deallocate(gbl_angle, gbl_dihed)

#endif /* not MPI */

! Initialize random number generator at same point in all processors. Then, if
! random initial velocities are needed, generate them in all processors. 
! In general, we must be careful to generate the same sequence of random
! numbers in all processors.

  call amrset(ig)

  if (ntx .eq. 1 .or. ntx .eq. 2) then
    call all_atom_setvel(natom, atm_vel, atm_mass_inv, tempi)
#ifdef CUDA
    call gpu_upload_vel(atm_vel)
#endif
  end if

  if (ibelly .gt. 0) then
    call all_atom_belly(natom, atm_igroup, atm_vel)
  end if

#ifdef CUDA
  call gpu_remap_bonded_setup(typ_ico, gbl_cn1, gbl_cn2, atm_numex, gbl_natex, ntypes)
#endif

  return

end subroutine gb_alltasks_setup

end module gb_alltasks_setup_mod
