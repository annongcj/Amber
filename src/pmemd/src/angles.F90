#include "copyright.i"
#include "hybrid_datatypes.i"

!*******************************************************************************
!
! Module:  angles_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module angles_mod

  use gbl_datatypes_mod
#ifdef _OPENMP_
  use omp_lib_kinds
#endif

  implicit none

! The following are derived from prmtop angle info:

  integer, save                         :: cit_ntheth, cit_ntheta

  type(angle_rec), allocatable, save    :: cit_angle(:)
#ifdef GBTimer
  integer                :: wall_sec, wall_usec
  double precision, save :: start_time_ms
#endif
 
#ifdef _OPENMP_
  double precision, allocatable, save   :: eadev_arr(:)
  double precision, allocatable, save   :: angle_energy_arr(:)
  integer(kind=omp_lock_kind), allocatable, private  :: omplk(:)
#endif

contains

!*******************************************************************************
!
! Subroutine:  angles_setup
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine angles_setup(num_ints, num_reals, use_atm_map)

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
#ifdef _OPENMP_
  integer               :: i
#endif
  type(angle_rec)       :: angles_copy(ntheth + ntheta)
  integer               :: atm_i, atm_j, atm_k, angles_idx
  integer               :: my_angle_cnt

! This routine can handle reallocation, and thus can be called multiple
! times.

! Find all angles for which this process owns either atom:

#ifdef _OPENMP_
if(using_gb_potential) then
  if (flag_angles .eqv. .true.) then
    flag_angles = .false. 
    allocate(eadev_arr(natom), &  
           angle_energy_arr(natom), &         
           omplk(natom), &
           stat = alloc_failed)
    if (alloc_failed .ne. 0) call setup_alloc_error
    do i = 1, natom
     call omp_init_lock(omplk(i))
    end do
  end if   
end if
#endif
  my_angle_cnt = 0

  do angles_idx = 1, ntheth

    atm_i = gbl_angle(angles_idx)%atm_i
    atm_j = gbl_angle(angles_idx)%atm_j
    atm_k = gbl_angle(angles_idx)%atm_k

#if defined(MPI) && !defined(CUDA)
#ifdef _OPENMP_
  if(using_gb_potential) then
    if(master) then
      my_angle_cnt = my_angle_cnt + 1
      angles_copy(my_angle_cnt) = gbl_angle(angles_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
    end if
  else
#endif
    if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
      my_angle_cnt = my_angle_cnt + 1
      angles_copy(my_angle_cnt) = gbl_angle(angles_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
    else if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
      my_angle_cnt = my_angle_cnt + 1
      angles_copy(my_angle_cnt) = gbl_angle(angles_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
    else if (gbl_atm_owner_map(atm_k) .eq. mytaskid) then
      my_angle_cnt = my_angle_cnt + 1
      angles_copy(my_angle_cnt) = gbl_angle(angles_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
    end if
#ifdef _OPENMP_
  end if
#endif
#else
    my_angle_cnt = my_angle_cnt + 1
    angles_copy(my_angle_cnt) = gbl_angle(angles_idx)
    use_atm_map(atm_i) = 1
    use_atm_map(atm_j) = 1
    use_atm_map(atm_k) = 1
#endif

  end do

  cit_ntheth = my_angle_cnt

  do angles_idx = anglea_idx, anglea_idx + ntheta - 1

    atm_i = gbl_angle(angles_idx)%atm_i
    atm_j = gbl_angle(angles_idx)%atm_j
    atm_k = gbl_angle(angles_idx)%atm_k

#if defined(MPI) && !defined(CUDA) 
#ifdef _OPENMP_
  if(using_gb_potential) then
    if(master) then
      my_angle_cnt = my_angle_cnt + 1
      angles_copy(my_angle_cnt) = gbl_angle(angles_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
    end if
  else
#endif
    if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
      my_angle_cnt = my_angle_cnt + 1
      angles_copy(my_angle_cnt) = gbl_angle(angles_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
    else if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
      my_angle_cnt = my_angle_cnt + 1
      angles_copy(my_angle_cnt) = gbl_angle(angles_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
    else if (gbl_atm_owner_map(atm_k) .eq. mytaskid) then
      my_angle_cnt = my_angle_cnt + 1
      angles_copy(my_angle_cnt) = gbl_angle(angles_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
    end if
#ifdef _OPENMP_
  end if
#endif
#else
    my_angle_cnt = my_angle_cnt + 1
    angles_copy(my_angle_cnt) = gbl_angle(angles_idx)
    use_atm_map(atm_i) = 1
    use_atm_map(atm_j) = 1
    use_atm_map(atm_k) = 1
#endif

  end do

  cit_ntheta = my_angle_cnt - cit_ntheth

  if (my_angle_cnt .gt. 0) then
    if (allocated(cit_angle)) then
      if (size(cit_angle) .lt. my_angle_cnt) then
        num_ints = num_ints - size(cit_angle) * angle_rec_ints
        deallocate(cit_angle)
        allocate(cit_angle(my_angle_cnt), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
        num_ints = num_ints + size(cit_angle) * angle_rec_ints
      end if
    else
      allocate(cit_angle(my_angle_cnt), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(cit_angle) * angle_rec_ints
    end if
    cit_angle(1:my_angle_cnt) = angles_copy(1:my_angle_cnt)
  end if

#ifdef CUDA
    call gpu_angles_setup(my_angle_cnt, cit_ntheth, cit_angle, gbl_teq, gbl_tk)
#endif

  return

end subroutine angles_setup

!*******************************************************************************
!
! Subroutine:  get_angle_energy
!
! Description:  Routine to get the bond energies and forces for potentials of
!               the type ct*(t-t0)**2.
!
!*******************************************************************************

subroutine get_angle_energy(angle_cnt, angle, x, frc, angle_energy)
#include  "angles.i"
end subroutine get_angle_energy

#ifdef _OPENMP_
!*******************************************************************************
!
! Subroutine:  get_angle_energy_gb
!
! Description:  Routine to get the bond energies and forces for potentials of
!               the type ct*(t-t0)**2. This is an hybrid implementation with
!               OpenMP+MPI.
!
!*******************************************************************************

subroutine get_angle_energy_gb(angle_cnt, angle, x, frc, angle_energy)
#define GBorn
#include  "angles.i"
#undef GBorn
end subroutine get_angle_energy_gb
#endif /*_OPENMP_*/

end module angles_mod
