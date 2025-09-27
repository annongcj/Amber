#include "copyright.i"
#include "dbg_arrays.i"

!*******************************************************************************
!
! Module: dbg_arrays_mod
!
! Description: Arrays debugging by dumping selected arrays contents into files
!              to compare the program status at different stages
!              
!*******************************************************************************
#ifdef DBG_ARRAYS
#ifdef MPI
module dbg_arrays_mod

  use parallel_dat_mod, only: mytaskid
!use ...

  implicit none


  logical, private, save                :: is_initialized=.false.
  logical                               :: dbg_arrays_enabled
#if defined(DBG_ARRAYS_STEPS) & defined(DBG_ARRAYS_STEPS_CNT)
  integer                               :: dbg_arrays_steps(DBG_ARRAYS_STEPS_CNT)=(/ DBG_ARRAYS_STEPS /) 
#else
  integer                               :: dbg_arrays_steps(3)=(/ 0, 1, 2 /) 
#endif

!write(0,*) __FILE__,__LINE__
contains

!*******************************************************************************
!
! Subroutine:  init_dbg_arrays
!
! Description: one time initialization 
!*******************************************************************************
subroutine init_dbg_arrays()
  implicit none

  if(mytaskid==0) print*, "initializing the debugging arrays"
#if defined(DBG_ARRAYS_STEPS) & defined(DBG_ARRAYS_STEPS_CNT)
  if(mytaskid==0) print*, "Printing info from time steps:", DBG_ARRAYS_STEPS
#else
  if(mytaskid==0) print*, "Printing info from time steps:", dbg_arrays_steps
#endif
end subroutine init_dbg_arrays

!*******************************************************************************
!
! Subroutine:  timetamp_dbg_arrays
!
! Description: write the timestep at all the debugging files
!*******************************************************************************
subroutine timestep_dbg_arrays(nstep)
  implicit none
  include 'mpif.h'

  Integer nstep 

  !local variables
  Integer err_code_mpi

  if(.not. is_initialized) then
    is_initialized=.true.
    call init_dbg_arrays
  end if

  call mpi_barrier(pmemd_comm, err_code_mpi)
  dbg_arrays_enabled=.false.
  if(ANY(dbg_arrays_steps==nstep)) dbg_arrays_enabled=.true.
  if(dbg_arrays_enabled) write(DBG_ARRAYS+mytaskid,'(A,I0.10)') "timestep", nstep


end subroutine timestep_dbg_arrays

!*******************************************************************************
!
! Subroutine:  DBG_ARRAYS_DUMP_CRD_FUN
!
! Description:   dumps the contents of three components double array
!*******************************************************************************
subroutine DBG_ARRAYS_DUMP_CRD_FUN(tag, id, dble3, N, wrap)
  implicit none
  Character(len=*)  :: tag
  Integer            :: N
  Integer            :: id(N)
  Double precision   :: dble3(3,N)
  Double precision   :: wrap(3,N)

  ! local variables
  integer i
  Double precision   :: tmp(3,N)
  if(.not. dbg_arrays_enabled) return

  tmp(:,:) = dble3(:,1:N)+wrap(:,1:N)
  call DBG_ARRAYS_DUMP_3DBLE_FUN(tag, id, tmp, N)

end subroutine DBG_ARRAYS_DUMP_CRD_FUN

!*******************************************************************************
!
! Subroutine:  DBG_ARRAYS_DUMP_3DBLE_FUN
!
! Description:   dumps the contents of three components double array
!*******************************************************************************
subroutine DBG_ARRAYS_DUMP_3DBLE_FUN(tag, id, dble3, N)
  implicit none
  Character(len=*)  :: tag
  Integer            :: N
  Integer            :: id(N)
  Double precision   :: dble3(3,N)

  ! local variables
  integer i
  if(.not. dbg_arrays_enabled) return

  do i=1, N
    write(DBG_ARRAYS+mytaskid,'(A," ",I0.10,":",3F)') tag, id(i), dble3(:,i)
  end do

end subroutine DBG_ARRAYS_DUMP_3DBLE_FUN

end module dbg_arrays_mod

#endif /*DBG_ARRAYS*/
#endif /*MPI*/
