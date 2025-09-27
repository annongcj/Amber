#include "copyright.i"

!*******************************************************************************
!
! Module: gb_parallel_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module gb_parallel_mod

  use parallel_dat_mod

  implicit none

  private       ! Everything here is private unless specified public

#if defined(MPI)

 integer, allocatable, save    :: atm_offsets(:)                       
 integer, allocatable, save    :: vec_offsets(:)                       
 integer, allocatable, save    :: vec_rcvcnts(:)                       

  public        :: gb_parallel_setup, &
                   gb_frcs_distrib, &
                   gb_mpi_allgathervec, &
                   gb_mpi_gathervec,  &
                   gb_amd_apply_weight, &
                   gb_gamd_apply_weight

contains

!*******************************************************************************
!
! Subroutine:  gb_parallel_setup
!
! Description:  Set up gb-specific data structures and do initial workload
!               division for generalized Born method.  Currently, there is
!               no loadbalancing after the initial setup.
!*******************************************************************************

subroutine gb_parallel_setup(num_ints, num_reals)

  use gbl_constants_mod
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use pmemd_lib_mod
  use prfs_mod
  use prmtop_dat_mod
  use ti_mod

  implicit none

! Formal arguments:

! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed
  integer               :: i
  integer               :: task_id 

! The following checks for minimum atoms and residues per
! processor are intended to avoid hitting conditions we may not
! have considered.  The use of more processors than atoms would be a bit
! ridiculous; with GB though, we may actually get close to 1 processor per
! residue.  We thus allow a lower limit there of 1.01 residue per processor.

! Check for sufficient atoms for the number of processors:

#ifdef CUDA
  if (natom .lt. 32 * numtasks) then
    write(mdout, '(a,a)') error_hdr, 'Must have 32x more atoms than processors!'
    call mexit(6, 1)
  end if
#else
  if (natom .lt. 10 * numtasks) then
    write(mdout, '(a,a)') error_hdr, 'Must have 10x more atoms than processors!'
    call mexit(6, 1)
  end if
#endif

! Check for sufficient prf's for the number of processors:

  if (dble(gbl_prf_cnt) .lt. dble(numtasks) * 1.01d0) then
    write(mdout, '(a,a)') error_hdr, &
      'Must have 1.01x more pseudo-residue fragments than processors!'
    call mexit(6, 1)
  end if

! Allocate storage for various structures associated with keeping track of the
! division of the atom workload.

  allocate(vec_offsets(0:numtasks), &
           atm_offsets(0:numtasks), &
           vec_rcvcnts(0:numtasks), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_ints = num_ints + size(vec_offsets) + &
                        size(atm_offsets) + &
                        size(vec_rcvcnts)

  ! Do prf-based atom division.

  call gb_atom_division

  ! All the other data structures can be set up using vec_rcvcnts.

  my_atm_cnt = vec_rcvcnts(mytaskid)

  atm_offsets(0) = 0

  do task_id = 0, numtasks - 1
    atm_offsets(task_id + 1) = atm_offsets(task_id) + vec_rcvcnts(task_id)
  end do

  ! The atom list is not strictly necessary in a GB context, as atoms
  ! are owned in contiguous blocks.  However, the guts of pmemd is set up
  ! without this assumption, so...
  
  do i = 1, my_atm_cnt
    gbl_my_atm_lst(i) = atm_offsets(mytaskid) + i
  end do

  do task_id = 0, numtasks - 1
    vec_rcvcnts(task_id) = 3 * vec_rcvcnts(task_id)
  end do

  vec_rcvcnts(numtasks) = 0

  vec_offsets(0) = 0

  do task_id = 0, numtasks - 1
    vec_offsets(task_id + 1) = vec_offsets(task_id) + vec_rcvcnts(task_id)
  end do

  ! Allocate buffers for mpi i/o that will remain allocated throughout the
  ! run.  This is done because some mpi implementations (myrinet in particular)
  ! mmap the buffers, and using stack buffers could have negative performance
  ! impacts in result (ifort also mmap/munmaps the dynamic stack space).
  ! These buffers are not used when other static allocations are available, say
  ! when broadcasting initialization data from the master.  All i/o using
  ! these buffers must be completed before routine exit (ie., waits MUST be
  ! done on nonblocking i/o).  We must be sure that the buffer size variables
  ! have been updated before this call...
  call set_minimum_mpi_bufs_size(3 * natom, num_reals)

  return

contains

!*******************************************************************************
!
! Internal Subroutine:  gb_atom_division
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine gb_atom_division
#ifdef _OPENMP_
  use omp_lib
#endif 

  implicit none

  integer               :: task_id
  integer               :: prf_id, first_prf_id, last_prf_id
  integer               :: offset, listcnt, atm_id, i
  integer               :: prfs_assigned
  integer               :: task_prf_cnt
  integer               :: task_prf_cnts(0:numtasks - 1)
#ifdef _OPENMP_
  integer               :: numthreads 
  integer               :: tid_prf_cnt
  integer,allocatable   :: tid_prf_cnts(:)
#endif
  double precision      :: per_task_target
  double precision      :: total_target

  ! BUGBUG - We should further improve on this method, actually targetting an
  !          atom count, and rounding up or down at the prf to hit the target...

  per_task_target = dble(gbl_prf_cnt) / dble(numtasks)
  total_target = 0.d0
  prfs_assigned = 0

#ifdef _OPENMP_
  task_prf_cnts(:) = 0 

#ifdef EMIL
  if ( emil_do_calc .gt. 0 ) then
     calc_emil = .true.
  end if
#endif
#endif
  do task_id = 0, numtasks - 2
    total_target = total_target + per_task_target
    task_prf_cnt = int(total_target) - prfs_assigned
    prfs_assigned = prfs_assigned + task_prf_cnt
    task_prf_cnts(task_id) = task_prf_cnt
  end do

  task_prf_cnts(numtasks - 1) = gbl_prf_cnt - prfs_assigned

  ! Initialize the various data structures that will reflect the atom division:

  vec_rcvcnts(:) = 0
  gbl_atm_owner_map(:) = -1

  ! Use the task_prf_cnts() array to do the division...

  first_prf_id = 1

  do task_id = 0, numtasks - 1
    vec_rcvcnts(task_id) = 0
    last_prf_id = first_prf_id + task_prf_cnts(task_id) - 1
    do prf_id = first_prf_id, last_prf_id
      offset = gbl_prf_listdata(prf_id)%offset
      listcnt = gbl_prf_listdata(prf_id)%cnt
      vec_rcvcnts(task_id) = vec_rcvcnts(task_id) + listcnt
      do i = offset + 1, offset + listcnt
        gbl_atm_owner_map(gbl_prf_lists(i)) = task_id
      end do
    end do
    first_prf_id = last_prf_id + 1
  end do

#ifdef _OPENMP_
!$omp parallel 
     numthreads = omp_get_num_threads()
!$omp end parallel 

  allocate(tid_prf_cnts(0:numthreads-1), &
           stat = alloc_failed)
  if (alloc_failed .ne. 0) call setup_alloc_error

  per_task_target = dble(gbl_prf_cnt) / dble(numthreads)
  total_target = 0.d0
  prfs_assigned = 0

  do task_id = 0, numthreads - 2
    total_target = total_target + per_task_target
    tid_prf_cnt = int(total_target) - prfs_assigned
    prfs_assigned = prfs_assigned + tid_prf_cnt
    tid_prf_cnts(task_id) = tid_prf_cnt
 end do

  tid_prf_cnts(numthreads - 1) = gbl_prf_cnt - prfs_assigned
  first_prf_id = 1

  do task_id = 0, numthreads - 1
    last_prf_id = first_prf_id + tid_prf_cnts(task_id) - 1
    do prf_id = first_prf_id, last_prf_id
      offset = gbl_prf_listdata(prf_id)%offset
      listcnt = gbl_prf_listdata(prf_id)%cnt
      do i = offset + 1, offset + listcnt
        gbl_atm_tid_map(gbl_prf_lists(i)) = task_id
      end do
    end do
    first_prf_id = last_prf_id + 1
  end do
#endif /*_OPENMP_*/
  
  return

end subroutine gb_atom_division

end subroutine gb_parallel_setup

!*******************************************************************************
!
! Subroutine:  gb_frcs_distrib
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine gb_frcs_distrib(atm_cnt, frc, recv_buf)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: frc(3 * atm_cnt)
  double precision      :: recv_buf(*)

! Local variables:

  integer               :: my_recvoff
  integer               :: my_recvcnt

  my_recvoff = vec_offsets(mytaskid)
  my_recvcnt = vec_rcvcnts(mytaskid)

  call mpi_reduce_scatter(frc, recv_buf, vec_rcvcnts, &
                          mpi_double_precision, mpi_sum, pmemd_comm, &
                          err_code_mpi)

  frc(my_recvoff + 1 : my_recvoff + my_recvcnt) = &
    recv_buf(1 : my_recvcnt)

  return

end subroutine gb_frcs_distrib

!*******************************************************************************
!
! Subroutine:  gb_mpi_allgathervec
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine gb_mpi_allgathervec(atm_cnt, vec)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: vec(3, atm_cnt)

! Local variables:

  double precision      :: recv_buf(3, atm_cnt)

  call mpi_allgatherv(vec(1, atm_offsets(mytaskid) + 1), &
                      vec_rcvcnts(mytaskid), mpi_double_precision, &
                      recv_buf, vec_rcvcnts, vec_offsets, &
                      mpi_double_precision, pmemd_comm, err_code_mpi)

  vec = recv_buf

  return

end subroutine gb_mpi_allgathervec

!*******************************************************************************
!
! Subroutine:  gb_mpi_gathervec
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine gb_mpi_gathervec(atm_cnt, vec)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: vec(3, atm_cnt)

! Local variables:

  double precision      :: recv_buf(3, atm_cnt)

  call mpi_gatherv(vec(1, atm_offsets(mytaskid) + 1), &
                   vec_rcvcnts(mytaskid), mpi_double_precision, &
                   recv_buf, vec_rcvcnts, vec_offsets, &
                   mpi_double_precision, 0, pmemd_comm, err_code_mpi)

  if (mytaskid == 0) then
    vec = recv_buf
  end if


  return

end subroutine gb_mpi_gathervec

!*******************************************************************************
!
! Subroutine:  gb_amd_apply_weight
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine gb_amd_apply_weight(atm_cnt,frc,fwgt)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: frc(3 * atm_cnt)
  double precision      :: fwgt

! Local variables:

  integer               :: my_recvoff
  integer               :: my_recvcnt

  my_recvoff = vec_offsets(mytaskid)
  my_recvcnt = vec_rcvcnts(mytaskid)

  frc(my_recvoff + 1 : my_recvoff + my_recvcnt) = &
  frc(my_recvoff + 1 : my_recvoff + my_recvcnt) * fwgt

  return

end subroutine gb_amd_apply_weight

!*******************************************************************************
!
! Subroutine:  gb_gamd_apply_weight
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine gb_gamd_apply_weight(atm_cnt,frc,fwgt)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: frc(3 * atm_cnt)
  double precision      :: fwgt

! Local variables:

  integer               :: my_recvoff
  integer               :: my_recvcnt

  my_recvoff = vec_offsets(mytaskid)
  my_recvcnt = vec_rcvcnts(mytaskid)

  frc(my_recvoff + 1 : my_recvoff + my_recvcnt) = &
  frc(my_recvoff + 1 : my_recvoff + my_recvcnt) * fwgt

  return

end subroutine gb_gamd_apply_weight

#endif /* MPI */

end module gb_parallel_mod
