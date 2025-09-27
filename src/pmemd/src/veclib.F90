#include "copyright.i"

#ifndef MKL
!  vectorized function calls

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ vectorized exponential

subroutine vdexp(n, x, y)
   
  implicit none

  integer               :: n
  double precision      :: x(n), y(n)
   
  y(1:n) = exp(x(1:n))

  return

end subroutine vdexp 
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ vectorized logarithm

subroutine vdln(n, x, y)
   
  implicit none

  integer               :: n
  double precision      :: x(n), y(n)
   
  y(1:n) = log(x(1:n))
   
  return

end subroutine vdln 
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ vectorized inverse square root

subroutine vdinvsqrt(n, x, y)
   
  implicit none

  integer               :: n
  double precision      :: x(n), y(n)
   
  y(1:n) = 1.d0 / sqrt(x(1:n))

  return

end subroutine vdinvsqrt 
!--------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ vectorized inverse

subroutine vdinv(n, x, y)
   
  implicit none

  integer               :: n
  double precision      :: x(n), y(n)
   
  y(1:n) = 1.d0 / x(1:n)

  return

end subroutine vdinv 
#else

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ dummy subroutine to avoid compiling an empty file

subroutine vd_dummy_to_avoid_empty_file()
   
  implicit none

  return

end subroutine vd_dummy_to_avoid_empty_file 
!--------------------------------------------------------------
#endif
