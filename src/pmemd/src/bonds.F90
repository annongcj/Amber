#include "copyright.i"
#include "hybrid_datatypes.i"

!*******************************************************************************
!
! Module:  bonds_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module bonds_mod

  use gbl_datatypes_mod
#ifdef _OPENMP_
  use omp_lib_kinds
#endif

  implicit none

! The following are derived from prmtop bond info:

  integer, save                         :: cit_nbonh, cit_nbona

  type(bond_rec), allocatable, save     :: cit_h_bond(:)
  type(bond_rec), allocatable, save     :: cit_a_bond(:)

#ifdef GBTimer
  integer                :: wall_sec, wall_usec
  double precision, save :: start_time_ms
#endif

#ifdef _OPENMP_
  double precision, allocatable                      :: lbe_arr(:)
  integer(kind=omp_lock_kind), allocatable, private  :: omplk(:)
#endif

contains

!*******************************************************************************
!
! Subroutine:  bonds_setup
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine bonds_setup(num_ints, num_reals, use_atm_map)

  use parallel_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
#if defined(_OPENMP_) || defined(CUDA)
  use mdin_ctrl_dat_mod
#endif

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals
  integer                       :: use_atm_map(natom)

! Local variables:

  integer               :: alloc_failed
  type(bond_rec)        :: bonds_copy(nbonh + nbona)
#ifdef _OPENMP_
  integer               :: i 
#endif

  ! This routine can handle reallocation, and thus can be called multiple
  ! times.

  call find_my_bonds(nbonh, gbl_bond, cit_nbonh, bonds_copy, use_atm_map)

  if (cit_nbonh .gt. 0) then
    if (allocated(cit_h_bond)) then
      if (size(cit_h_bond) .lt. cit_nbonh) then
        num_ints = num_ints - size(cit_h_bond) * bond_rec_ints
        deallocate(cit_h_bond)
        allocate(cit_h_bond(cit_nbonh), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
        num_ints = num_ints + size(cit_h_bond) * bond_rec_ints
      end if
    else
      allocate(cit_h_bond(cit_nbonh), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(cit_h_bond) * bond_rec_ints
    end if
    cit_h_bond(1:cit_nbonh) = bonds_copy(1:cit_nbonh)
  end if

  if  (cit_nbona+nbona .gt. 0) &
  call find_my_bonds(nbona, gbl_bond(bonda_idx), cit_nbona, bonds_copy, &
                         use_atm_map)

  if (cit_nbona .gt. 0) then
    if (allocated(cit_a_bond)) then
      if (size(cit_a_bond) .lt. cit_nbona) then
        num_ints = num_ints - size(cit_a_bond) * bond_rec_ints
        deallocate(cit_a_bond)
        allocate(cit_a_bond(cit_nbona), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
        num_ints = num_ints + size(cit_a_bond) * bond_rec_ints
      end if
    else
      allocate(cit_a_bond(cit_nbona), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(cit_a_bond) * bond_rec_ints
    end if
    cit_a_bond(1:cit_nbona) = bonds_copy(1:cit_nbona)
  end if

#ifdef CUDA
    call gpu_bonds_setup(cit_nbona, cit_a_bond, cit_nbonh, cit_h_bond, gbl_req, gbl_rk);    
#endif

#ifdef _OPENMP_
  if(using_gb_potential) then
    if (flag_bond .eqv. .true.) then
     flag_bond = .false. ! 
     allocate(lbe_arr(natom), &
              omplk(natom), &
           stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error

      do i = 1, natom
        call omp_init_lock(omplk(i))
      end do
    end if
  end if
#endif /*_OPENMP_*/

  return

end subroutine bonds_setup

!*******************************************************************************
!
! Subroutine:  find_my_bonds
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine find_my_bonds(bond_cnt, bonds, my_bond_cnt, my_bonds, &
                             use_atm_map)

  use parallel_dat_mod
#ifdef _OPENMP_
  use mdin_ctrl_dat_mod
#endif

  implicit none

! Formal arguments:

  integer               :: bond_cnt
  type(bond_rec)        :: bonds(bond_cnt)
  integer               :: my_bond_cnt
  type(bond_rec)        :: my_bonds(*)
  integer               :: use_atm_map(*)

! Local variables:

  integer               :: atm_i, atm_j, bonds_idx

! Find all bonds for which this process owns either atom:

  my_bond_cnt = 0

  do bonds_idx = 1, bond_cnt

    atm_i = bonds(bonds_idx)%atm_i
    atm_j = bonds(bonds_idx)%atm_j

#if defined(MPI) && !defined(CUDA)
#ifdef _OPENMP_
  if(using_gb_potential) then
   if(master) then
      my_bond_cnt = my_bond_cnt + 1
      my_bonds(my_bond_cnt) = bonds(bonds_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
   end if
  else
#endif
    if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
      my_bond_cnt = my_bond_cnt + 1
      my_bonds(my_bond_cnt) = bonds(bonds_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
    else if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
      my_bond_cnt = my_bond_cnt + 1
      my_bonds(my_bond_cnt) = bonds(bonds_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
    end if
#ifdef _OPENMP_
  end if
#endif
#else
    my_bond_cnt = my_bond_cnt + 1
    my_bonds(my_bond_cnt) = bonds(bonds_idx)
    use_atm_map(atm_i) = 1
    use_atm_map(atm_j) = 1
#endif

  end do

  return

end subroutine find_my_bonds

!*******************************************************************************
!
! Subroutine:  get_bond_energy
!
! Description:
!              
! Routine to get bond energy and forces for the potential of cb*(b-b0)**2.
!
!*******************************************************************************

subroutine get_bond_energy(bond_cnt, bond, x, frc, bond_energy)
#include "bonds.i"
end subroutine get_bond_energy

#ifdef _OPENMP_
subroutine get_bond_energy_gb(bond_cnt, bond, x, frc, bond_energy)
#define GBorn
#include "bonds.i"
#undef GBorn
end subroutine get_bond_energy_gb
#endif /*_OPENMP_*/

end module bonds_mod
