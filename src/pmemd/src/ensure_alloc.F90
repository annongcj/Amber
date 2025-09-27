module ensure_alloc_mod

 implicit none
real  :: alloc_factor=1.2

contains

!*******************************************************************************
!
! Subroutines: ensure_alloc_*
!
! Description: tests whether the array size is sufficient and reallocates if
!              needed. The data is copied to the new array when retain_data=true
!
!*******************************************************************************

subroutine ensure_alloc_float4(ar, new_size, retain_data_)
  implicit none
  real, allocatable :: ar(:,:)
  integer                       :: new_size
  logical,optional              :: retain_data_

  ! local variables
  real, allocatable :: tmp(:,:)
  integer          :: old_size
  logical          :: retain_data=.false.

  old_size = size(ar,2)
  if(present(retain_data_)) retain_data = retain_data_

  if(.not. allocated(ar)) then
    allocate(ar(4,int(new_size*alloc_factor)))
  else if(new_size .gt. old_size) then ! Resize array
    if(retain_data) then
      allocate(tmp(4,old_size))
      tmp = ar
      deallocate(ar)
      allocate(ar(4,int(new_size*alloc_factor)))
      ar(:,1:old_size) = tmp
      deallocate(tmp)
    else
      deallocate(ar)
      allocate(ar(4,int(new_size*alloc_factor)))
    endif
  endif
end subroutine ensure_alloc_float4

subroutine ensure_alloc_float3(ar, new_size, retain_data_)
  implicit none
  real, allocatable :: ar(:,:)
  integer                       :: new_size
  logical,optional              :: retain_data_

  ! local variables
  real, allocatable :: tmp(:,:)
  integer          :: old_size
  logical          :: retain_data=.false.

  old_size = size(ar,2)
  if(present(retain_data_)) retain_data = retain_data_

  if(.not. allocated(ar)) then
    allocate(ar(3,int(new_size*alloc_factor)))
  else if(new_size .gt. old_size) then ! Resize array
    if(retain_data) then
      allocate(tmp(3,old_size))
      tmp = ar
      deallocate(ar)
      allocate(ar(3,int(new_size*alloc_factor)))
      ar(:,1:old_size) = tmp
      deallocate(tmp)
    else
      deallocate(ar)
      allocate(ar(3,int(new_size*alloc_factor)))
    endif
  endif
end subroutine ensure_alloc_float3

subroutine ensure_alloc_dble3(ar, new_size, retain_data_)
  implicit none
  double precision, allocatable :: ar(:,:)
  integer                       :: new_size
  logical,optional              :: retain_data_

  ! local variables
  double precision, allocatable :: tmp(:,:)
  integer          :: old_size
  logical          :: retain_data=.false.

  old_size = size(ar,2)
  if(present(retain_data_)) retain_data = retain_data_

  if(.not. allocated(ar)) then
    allocate(ar(3,int(new_size*alloc_factor)))
  else if(new_size .gt. old_size) then ! Resize array
    if(retain_data) then
      allocate(tmp(3,old_size))
      tmp(:,1:old_size) = ar(:,1:old_size)
      deallocate(ar)
      allocate(ar(3,int(new_size*alloc_factor)))
      ar(1:3,1:old_size) = tmp(1:3,1:old_size)
      deallocate(tmp)
    else
      deallocate(ar)
      allocate(ar(3,int(new_size*alloc_factor)))
    endif
  endif
end subroutine ensure_alloc_dble3

subroutine ensure_alloc_dble2d(ar, new_size, size2, retain_data_)
  implicit none
  double precision, allocatable :: ar(:,:)
  integer                       :: new_size, size2
  logical,optional              :: retain_data_

  ! local variables
  double precision, allocatable :: tmp(:,:)
  integer          :: old_size
  logical          :: retain_data=.false.

  old_size = size(ar,1)
  if(present(retain_data_)) retain_data = retain_data_

  if(.not. allocated(ar)) then
    allocate(ar(int(new_size*alloc_factor), size2))
  else if(new_size .gt. old_size) then ! Resize array
    if(retain_data) then
      allocate(tmp(old_size, size2))
      tmp = ar
      deallocate(ar)
      allocate(ar(int(new_size*alloc_factor), size2))
      ar(1:old_size, :) = tmp
      deallocate(tmp)
    else
      deallocate(ar)
      allocate(ar(int(new_size*alloc_factor), size2))
    endif
  endif
end subroutine ensure_alloc_dble2d

subroutine ensure_alloc_dble(ar, new_size, retain_data_)
  implicit none
  double precision, allocatable :: ar(:)
  integer                       :: new_size
  logical,optional              :: retain_data_

  ! local variables
  double precision, allocatable :: tmp(:)
  integer          :: old_size
  logical          :: retain_data=.false.

  old_size = size(ar)
  if(present(retain_data_)) retain_data = retain_data_

  if(.not. allocated(ar)) then
    allocate(ar(int(new_size*alloc_factor)))
  else if(new_size .gt. old_size) then ! Resize array
    if(retain_data) then
      allocate(tmp(old_size))
      tmp = ar
      deallocate(ar)
      allocate(ar(int(new_size*alloc_factor)))
      ar(1:old_size) = tmp
      deallocate(tmp)
    else
      deallocate(ar)
      allocate(ar(int(new_size*alloc_factor)))
    endif
  endif
end subroutine ensure_alloc_dble

subroutine ensure_alloc_int3(ar, new_size, retain_data_)
  implicit none
  integer, allocatable :: ar(:,:)
  integer                       :: new_size
  logical,optional              :: retain_data_

  ! local variables
  integer, allocatable :: tmp(:,:)
  integer          :: old_size
  logical          :: retain_data=.false.

  old_size = size(ar,2)
  if(present(retain_data_)) retain_data = retain_data_

  if(.not. allocated(ar)) then
    allocate(ar(3,int(new_size*alloc_factor)))
  else if(new_size .gt. old_size) then ! Resize array
    if(retain_data) then
      allocate(tmp(3,old_size))
      tmp = ar
      deallocate(ar)
      allocate(ar(3,int(new_size*alloc_factor)))
      ar(:,1:old_size) = tmp
      deallocate(tmp)
    else
      deallocate(ar)
      allocate(ar(3,int(new_size*alloc_factor)))
    endif
  endif
end subroutine ensure_alloc_int3

subroutine ensure_alloc_int2d(ar, new_size, size2, retain_data_)
  implicit none
  integer, allocatable :: ar(:,:)
  integer                       :: new_size, size2
  logical,optional              :: retain_data_

  ! local variables
  integer, allocatable :: tmp(:,:)
  integer          :: old_size
  logical          :: retain_data=.false.

  old_size = size(ar,1)
  if(present(retain_data_)) retain_data = retain_data_

  if(.not. allocated(ar)) then
    allocate(ar(int(new_size*alloc_factor), size2))
  else if(new_size .gt. old_size) then ! Resize array
    if(retain_data) then
      allocate(tmp(old_size, size2))
      tmp = ar
      deallocate(ar)
      allocate(ar(int(new_size*alloc_factor), size2))
      ar(1:old_size, :) = tmp
      deallocate(tmp)
    else
      deallocate(ar)
      allocate(ar(int(new_size*alloc_factor), size2))
    endif
  endif
end subroutine ensure_alloc_int2d

subroutine ensure_alloc_int(ar, new_size, retain_data_)
  implicit none
  integer, allocatable :: ar(:)
  integer                       :: new_size
  logical,optional              :: retain_data_

  ! local variables
  integer, allocatable :: tmp(:)
  integer          :: old_size
  logical          :: retain_data=.false.

  old_size = size(ar)
  if(present(retain_data_)) retain_data = retain_data_

  if(.not. allocated(ar)) then
    allocate(ar(int(new_size*alloc_factor)))
  else if(new_size .gt. old_size) then ! Resize array
    if(retain_data) then
      allocate(tmp(old_size))
      tmp = ar
      deallocate(ar)
      allocate(ar(int(new_size*alloc_factor)))
      ar(1:old_size) = tmp
      deallocate(tmp)
    else
      deallocate(ar)
      allocate(ar(int(new_size*alloc_factor)))
    endif
  endif
end subroutine ensure_alloc_int

subroutine ensure_alloc_char4(ar, new_size, retain_data_)
  implicit none
  character(4), allocatable :: ar(:)
  integer                       :: new_size
  logical,optional              :: retain_data_

  ! local variables
  character(4), allocatable :: tmp(:)
  integer          :: old_size
  logical          :: retain_data=.false.

  old_size = size(ar)
  if(present(retain_data_)) retain_data = retain_data_

  if(.not. allocated(ar)) then
    allocate(ar(int(new_size*alloc_factor)))
  else if(new_size .gt. old_size) then ! Resize array
    if(retain_data) then
      allocate(tmp(old_size))
      tmp = ar
      deallocate(ar)
      allocate(ar(int(new_size*alloc_factor)))
      ar(1:old_size) = tmp
      deallocate(tmp)
    else
      deallocate(ar)
      allocate(ar(int(new_size*alloc_factor)))
    endif
  endif
end subroutine ensure_alloc_char4

end module ensure_alloc_mod
