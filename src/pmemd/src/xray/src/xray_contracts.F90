!
! This module should be updated to use the existing error reporting mechanism
! of pmemd; see the commented code below.  Consistency counts.  In addition,
! error stop is a F2008 feature that is not supported by all the compilers
! that Amber intends to support, e.g., intel 17.  
!
module xray_contracts_module

  ! Commented but preferred manner of error reporting:
  !use file_io_dat_mod, only : mdout
  !use gbl_constants_mod, only : error_hdr
  !use pmemd_lib_mod,   only : mexit

  implicit none
  
contains

subroutine check_requirement(condition, message)
  implicit none
  logical, intent(in) :: condition
  character(len = *), optional, intent(in) :: message

  if (.not. condition) then
    if (present(message)) then
      ! Commented but preferred manner of error reporting:
      !write(mdout, '(a,a)') error_hdr, 'Unmet requirement: ' // message
      !call mexit(6, 1)
      print message
      error stop 'Unmet requirement.'
    else
      stop "Unmet requirement."
    end if
  endif

end subroutine check_requirement

subroutine check_assertion(condition, message)
  implicit none
  logical, intent(in) :: condition
  character(len = *), optional, intent(in) :: message

#ifndef NDEBUG
  if (.not. condition) then
    if (present(message)) then
      print message
      error stop "Assertion failed: "
    else
      error stop "Assertion failed."
    end if
  endif
#endif

end subroutine check_assertion

subroutine check_precondition(condition, message)
  implicit none
  logical, intent(in) :: condition
  character(len=*), optional, intent(in) :: message

#ifndef NDEBUG
  if (.not. condition) then
    if (present(message)) then
      print message
      error stop "Precondition failed: "
    else
      error stop "Precondition failed."
    end if
  endif
#endif

end subroutine check_precondition

subroutine check_postcondition(condition, message)
  implicit none
  logical, intent(in) :: condition
  character(len=*), optional, intent(in) :: message

#ifndef NDEBUG
  if (.not. condition) then
    if (present(message)) then
      print message
      error stop "Postcondition failed: "
    else
      error stop "Postcondition failed."
    end if
  endif
#endif

end subroutine check_postcondition

end module xray_contracts_module
