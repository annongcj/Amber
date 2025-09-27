#include "copyright.i"

!*******************************************************************************
!
! Module:  random_mod
!
! Description: 
!
!*******************************************************************************

module random_mod

use file_io_dat_mod

  implicit none

#include "random.i"

  type(random_state), private, save :: default_rand_gen

contains

!*******************************************************************************
!
! Subroutine: amrset
!
! Description: Wrapper to call amrset_gen with default_rand_gen. This should
!              only be called once (to set the default random generator), so
!              we move gpu_amrset() call in here.
!
!*******************************************************************************

subroutine amrset(iseed)
  
  implicit none

  integer iseed

  call amrset_gen(default_rand_gen, iseed)

#ifdef CUDA
    call gpu_amrset(iseed)
#endif  

  return

end subroutine amrset

!*******************************************************************************
!
! Subroutine:  amrset_gen
!
! Description: 
!
! Initialization routine for Marsaglias random number generator
! as implemented in Amber 3.0 Rev A by George Seibel.  See doc in amrand.
!
! Testing:  Call amrset with iseed = 54185253.  This should result
! in is1 = 1802 and is2 = 9373.  Call amrand 20000 times, then six
! more times, printing the six random numbers * 2**24 (ie, 4096*4096)
! They should be: (6f12.1)
! 6533892.0  14220222.0  7275067.0  6172232.0  8354498.0  10633180.0
!              
! INPUT:        iseed
!
! OUTPUT:       to module private variables
!
!*******************************************************************************

subroutine amrset_gen(gen, iseed)

  implicit none

! Formal arguments:

  type(random_state), intent(inout) :: gen

  integer       iseed   ! integer seed (absolute value used)

! Local variables:

! Two internal seeds used in Marsaglia algorithm:

  integer           is1
  integer           is2

! Max value of first seed (is1), 31328:

  integer           is1max

! Max value of second seed (is2), 30081

  integer           is2max

  integer           i, ii, j, jj, k, l, m

  double precision  s, t

  data is1max, is2max /31328, 30081/

! Construct two internal seeds from single unbound Amber seed: 
!
! is1 and is2 are quotient and remainder of iseed/IS2MAX.  We add
! one to keep zero and one results from both mapping to one.
! max and min functions keep is1 and is2 in required bounds.
  iseed = abs(iseed)

  is1 = max((iseed / is2max) + 1, 1)
  is1 = min(is1, is1max)

  is2 = max(1, mod(iseed, is2max) + 1)
  is2 = min(is2, is2max)

  i = mod(is1/177, 177) + 2
  j = mod(is1    , 177) + 2
  k = mod(is2/169, 178) + 1
  l = mod(is2    , 169)

  do ii = 1, 97
    s = 0.0d0
    t = 0.5d0
    do jj = 1, 24
      m = mod(mod(i*j, 179)*k, 179)
      i = j
      j = k
      k = m
      l = mod(53*l + 1, 169)
      if (mod(l*m, 64) .ge. 32) s = s + t
      t = 0.5d0 * t
    end do
    gen%u(ii) = s
  end do

  gen%c  = 362436.0d0   / 16777216.0d0
  gen%cd = 7654321.0d0  / 16777216.0d0
  gen%cm = 16777213.0d0 / 16777216.0d0

  gen%i97 = 97
  gen%j97 = 33

  gen%set = .true.
  
  return

end subroutine amrset_gen

!*******************************************************************************
!
! Subroutine: amrand
!
! Description: Wrapper call for amrand_gen, using the default generator state
!
!*******************************************************************************

subroutine amrand(y)

  implicit none

  double precision, intent(out) :: y

  call amrand_gen(default_rand_gen, y)

  return

end subroutine amrand

!*******************************************************************************
!
! Subroutine:  amrand_gen
!
! Description: 
!
! Portable Random number generator by George Marsaglia
! Amber 3.0 Rev A implementation by George Seibel
!
! This random number generator originally appeared in *Toward a Universal
! Random Number Generator* by George Marsaglia and Arif Zaman.  Florida
! State University Report: FSU-SCRI-87-50 (1987)
!
! It was later modified by F. James and published in *A Review of Pseudo-
! random Number Generators*
!
! This is claimed to be the best known random number generator available.
! It passes ALL of the tests for random number generators and has a
! period of 2^144, is completely portable (gives bit identical results on
! all machines with at least 24-bit mantissas in the floating point
! representation).
!
! The algorithm is a combination of a Fibonacci sequence (with lags of 97
! and 33, and operation "subtraction plus one, modulo one") and an
! "arithmetic sequence" (using subtraction).
!
! INPUT:        from module private variables
!
! OUTPUT:       y:  A random number between 0.0 and 1.0
!
!*******************************************************************************

subroutine amrand_gen(gen, y)

  use pmemd_lib_mod
#ifndef NO_F95_RANDOM_NUMBER_INTRINSIC
  use mdin_ctrl_dat_mod, only : irandom
#endif

  implicit none

! Formal arguments:
  
  type(random_state), intent(inout) :: gen

  double precision, intent(out)     :: y

! Local variables:

  double precision  uni

#ifndef NO_F95_RANDOM_NUMBER_INTRINSIC
  ! if we are using the intrinsic random number generator there is only one
  ! stream that s advanced - multiple random number generators cannot be generated
  ! so gen is ignored here.

  if (irandom == 2) then
     call amrand_f95_intrinsic(y)
     return
  endif
#endif

  if (.not. gen%set) then
    write(mdout, '(a)') 'amrand not initialized'
    call mexit(6, 1)
  end if

  uni = gen%u(gen%i97) - gen%u(gen%j97)
  if (uni .lt. 0.0d0) uni = uni + 1.0d0
  gen%u(gen%i97) = uni
  gen%i97 = gen%i97 - 1
  if (gen%i97 .eq. 0) gen%i97 = 97
  gen%j97 = gen%j97 - 1
  if (gen%j97 .eq. 0) gen%j97 = 97
  gen%c = gen%c - gen%cd
  if (gen%c .lt. 0.0d0) gen%c = gen%c + gen%cm
  uni = uni - gen%c
  if (uni .lt. 0.0d0) uni = uni + 1.0d0
  y = uni

  return

end subroutine amrand_gen

!*******************************************************************************
!
! Subroutine: gauss
!
! Description: Wrapper for gauss_gen using default random generator state
!
!*******************************************************************************

subroutine gauss(am, sd, v)

  implicit none

  double precision, intent(in)                  :: am
  double precision, intent(in)                  :: sd
  double precision, optional, intent(out)       :: v

  call gauss_gen(default_rand_gen, am, sd, v)

  return

end subroutine gauss

!*******************************************************************************
!
! Subroutine:   gauss_gen
!
! Description:  Generate a pseudo-random Gaussian sequence, with mean am and
!               std. dev. sd.  This is a version of amrand() that adds the
!               constraint of a gaussian distribution, with mean "AM" and
!               standard deviation "SD".  Output is to variable "V". It also
!               requires amrset to have been called first, and "uses up" the
!               same sequence that amrand() does.
!
!*******************************************************************************

subroutine gauss_gen(gen, am, sd, v)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  type(random_state), intent(inout)             :: gen
  double precision, intent(in)                  :: am
  double precision, intent(in)                  :: sd
  double precision, optional, intent(out)       :: v

! Local variables:

  double precision              :: tmp1, tmp2
  double precision              :: uni
  double precision              :: zeta1, zeta2

  if (.not. gen%set) then
    write(mdout, '(a)') 'amrand not initialized!'
    call mexit(6, 1)
  end if

  ! Use the method of Box and Muller

  ! for some applications, one could use both "v" and "veven" in random
  ! sequence; but this won't work for most things we need (e.g. Langevin
  ! dynamics,) since the two adjacent variables are themselves highly
  ! correlated.  Hence we will just use the first ("v") variable.

  ! get two random numbers, even on (-1,1):

  do

    uni = gen%u(gen%i97) - gen%u(gen%j97)
    if (uni .lt. 0.0d0) uni = uni + 1.0d0
    gen%u(gen%i97) = uni
    gen%i97 = gen%i97 - 1
    if (gen%i97 .eq. 0) gen%i97 = 97
    gen%j97 = gen%j97 - 1
    if (gen%j97 .eq. 0) gen%j97 = 97
    gen%c = gen%c - gen%cd
    if (gen%c .lt. 0.0d0) gen%c = gen%c + gen%cm
    uni = uni - gen%c
    if (uni .lt. 0.0d0) uni = uni + 1.0d0
    zeta1 = uni + uni - 1.d0

    uni = gen%u(gen%i97) - gen%u(gen%j97)
    if (uni .lt. 0.0d0) uni = uni + 1.0d0
    gen%u(gen%i97) = uni
    gen%i97 = gen%i97 - 1
    if (gen%i97 .eq. 0) gen%i97 = 97
    gen%j97 = gen%j97 - 1
    if (gen%j97 .eq. 0) gen%j97 = 97
    gen%c = gen%c - gen%cd
    if (gen%c .lt. 0.0d0) gen%c = gen%c + gen%cm
    uni = uni - gen%c
    if (uni .lt. 0.0d0) uni = uni + 1.0d0
    zeta2 = uni + uni - 1.d0

    tmp1 = zeta1 * zeta1 + zeta2 * zeta2

    if (tmp1 .lt. 1.d0 .and. tmp1 .ne. 0.d0) then

      if (present(v)) then
        tmp2 = sd * sqrt(-2.d0 * log(tmp1)/tmp1)
        v = zeta1 * tmp2 + am
      end if

      return

    end if

  end do

end subroutine gauss_gen

#ifndef NO_F95_RANDOM_NUMBER_INTRINSIC
!*******************************************************************************
!Routines for using F95 intrinsic random number generator
!*******************************************************************************

!------
!This routine sets up the random seeds needed for the F95 intrinsic
!random number generator. The intrinsic is implementation specific
!and so needs different numbers of seeds on different hardware and
!different compilers. It is also not possible to have multiple instances
!of this random number generator. As such using it (irandom=2) necessarily
!breaks the test cases.
!
!The advantage is that this doesn't suffer from the issues with AMBER's own
!random number generator that generates too many whole fractions (including
!zero).
!
!We seed it by calling the AMBER random number generator with the ig set in
!&cntrl to get the number of random sees that we need. This at least makes it
!reproducible across runs on the same hardware / OS / Compiler.
!------
subroutine amrset_f95_intrinsic()

  use pmemd_lib_mod, only : alloc_error
  use mdin_ctrl_dat_mod, only : irandom

  implicit none

  integer, allocatable :: seed(:)
  integer :: alloc_failed
  integer :: n, i, x
  double precision :: y  

  !First figure out how many random seeds we need.

  call random_seed(size = n)
  allocate(seed(n), stat=alloc_failed)
  if (alloc_failed .ne. 0) call alloc_error("amrset_f95_intrinsic", "seed(n)")

  !We use AMBER random number generator to get the seeds - assumes it has
  !be initialized already.

  !set irandom=1 so amrand does not try to use the intrinsic

  irandom = 1
  do i = 1, n    
    call amrand(y)
    !disallow zero since it is not valid as a seed.
    do 
      x = floor(y * huge(x))
      if (x /= 0) exit
    end do

    seed(i) = x

  end do
  !reset irandom back to 2
  irandom = 2

  call random_seed(put=seed)

  return

end subroutine amrset_f95_intrinsic

subroutine amrand_f95_intrinsic(y)

  implicit none

  double precision, intent(out) :: y

  call random_number(y)

end subroutine amrand_f95_intrinsic

#else
!Dummy routines to allow compilation.
subroutine amrset_f95_intrinsic()

  implicit none

end subroutine amrset_f95_intrinsic

subroutine amrand_f95_intrinsic(y)

  implicit none

  double precision, intent(out) :: y

end subroutine amrand_f95_intrinsic

#endif
!*******************************************************************************
!End F95 Random number intrinsic routines.
!*******************************************************************************

end module random_mod
