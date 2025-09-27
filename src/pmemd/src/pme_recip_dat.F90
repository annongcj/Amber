#include "copyright.i"

!*******************************************************************************
!
! Module: pme_recip_dat_mod
!
! Description: <TBS>
!
! NOTE NOTE NOTE: This code assumes frc_int .eq. 0 and should only be used
!                 under these conditions!!!
!*******************************************************************************

module pme_recip_dat_mod

!precision setting for SPDP, DPDP, SPSP
#include "include_precision.i"
  implicit none

  double precision, allocatable, save           :: m1_exp_tbl(:)
  double precision, allocatable, save           :: m2_exp_tbl(:)
  double precision, allocatable, save           :: m3_exp_tbl(:)

  double precision, allocatable, save           :: gbl_prefac1(:)
  double precision, allocatable, save           :: gbl_prefac2(:)
  double precision, allocatable, save           :: gbl_prefac3(:)

#ifdef MPI

  logical, allocatable, save                    :: is_recip_task(:)
  logical, save                                 :: i_do_recip = .false.
  integer, save                                 :: recip_numtasks = 0
  integer, allocatable, save                    :: gbl_used_recip_img_lst(:)

#endif /* MPI */
  integer, save                                 :: gbl_used_recip_img_cnt = 0

contains

!*******************************************************************************
!
! Subroutine:  pme_recip_dat_setup
!
! Description:
!              
!*******************************************************************************

subroutine pme_recip_dat_setup(num_ints, num_reals)

  use pmemd_lib_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed

  allocate(gbl_prefac1(nfft1), &
           gbl_prefac2(nfft2), &
           gbl_prefac3(nfft3), &
#ifdef MPI
           gbl_used_recip_img_lst(natom), &
           is_recip_task(0: numtasks-1), &
#endif /* MPI */
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_reals = num_reals + size(gbl_prefac1) + size(gbl_prefac2) + &
                          size(gbl_prefac3)

  gbl_prefac1(:) = 0.d0
  gbl_prefac2(:) = 0.d0
  gbl_prefac3(:) = 0.d0

  call load_prefacs(gbl_prefac1, gbl_prefac2, gbl_prefac3)

#ifdef MPI
  num_ints = num_ints + size(gbl_used_recip_img_lst) + &
                        size(is_recip_task)

  gbl_used_recip_img_lst(:) = 0
  is_recip_task(:) = .false.
#else
  gbl_used_recip_img_cnt = natom
#endif /* MPI */

  if (is_orthog .ne. 0) then
    call allocate_m_exp_tbls
    call load_m_exp_tbls
  end if

  return

contains

!*******************************************************************************
!
! Internal Subroutine:  load_prefacs
!
! Description:  This routine loads the moduli of the inverse DFT of the
!               B splines.  bsp_mod1-3 hold these values, and nfft1-3 are the
!               grid dimensions.  Order is the order of the B spline approx.
!
!*******************************************************************************

subroutine load_prefacs(prefac1, prefac2, prefac3)

  use bspline_mod
  use mdin_ewald_dat_mod

  implicit none

  ! Formal arguments:

  double precision      :: prefac1(*), prefac2(*), prefac3(*)

  ! Local variables:

  double precision      :: array(bspl_order)
  double precision      :: w
  double precision      :: bsp_arr(max(nfft1,nfft2,nfft3))
  double precision      :: bsp_mod(max(nfft1,nfft2,nfft3))
  integer               :: i

  w = 0.d0
  call fill_bspline_0(w, bspl_order, array)

  bsp_arr(:) = 0.d0

  do i = 1, bspl_order - 1
   bsp_arr(i) = array(bspl_order - i)
  end do

  call dftmod(bsp_mod, bsp_arr, nfft1)

  call factor_lambda(bsp_mod, nfft1, prefac1)

  call dftmod(bsp_mod, bsp_arr, nfft2)

  call factor_lambda(bsp_mod, nfft2, prefac2)

  call dftmod(bsp_mod, bsp_arr, nfft3)

  call factor_lambda(bsp_mod, nfft3, prefac3)

  return

end subroutine load_prefacs

!*******************************************************************************
!
! Internal Subroutine:  dftmod
!
! Description:  Computes the modulus of the discrete fourier transform of
!               bsp_arr, storing it into bsp_mod.
!
!*******************************************************************************

subroutine dftmod(bsp_mod, bsp_arr, nfft)

  implicit none

  integer               nfft

  double precision      bsp_mod(nfft), bsp_arr(nfft)

  integer               j, k

  double precision      sum1, sum2, twopi, arg, tiny

  twopi = 2.d0 * PI

  tiny = 1.d-7

  do k = 1, nfft
    sum1 = 0.d0
    sum2 = 0.d0

    do j = 1, nfft
      arg = twopi * (k - 1) * (j - 1)/nfft
      sum1 = sum1 + bsp_arr(j) * cos(arg)
      sum2 = sum2 + bsp_arr(j) * sin(arg)
    end do

    bsp_mod(k) = sum1**2 + sum2**2

  end do

! Fix the ONE case where exponential Euler spline interpolation fails.
! This arbitrary assignment to avoid division by zero is okay
! since it happens with highest frequency reciprocal vectors out in tail
! of reciprocal sum.

  do k = 1, nfft
    if (bsp_mod(k) .lt. tiny) bsp_mod(k) = 0.5d0 * (bsp_mod(k-1) + bsp_mod(k+1))
  end do

  return

end subroutine dftmod

!*******************************************************************************
!
! Internal Subroutine:  factor_lambda
!
! Description:  Factor in optimal lambda coefficient for bspline approximant
!               of complex exponentials, thus modify influence function.
!
!*******************************************************************************

subroutine factor_lambda(bsp_mod, nfft, prefac)

  use mdin_ewald_dat_mod

  implicit none

  integer               nfft

  double precision      bsp_mod(nfft), prefac(nfft)

  integer               k, nf, m, order2

  double precision      lambda, gsum, gsum2

  integer, parameter    :: kcut = 50

  nf = nfft / 2
  order2 = 2 * bspl_order

  do k = 1, nfft
    m = k - 1
    if (k .gt. nf)m = k - 1 - nfft
    order2 = 2 * bspl_order
    if (m .eq. 0) then
      lambda = 1.d0
    else
      call gamma_sum(gsum, m, nfft, bspl_order, kcut)
      call gamma_sum(gsum2, m, nfft, order2, kcut)
      lambda = gsum/gsum2
    end if
    prefac(k) = lambda**2/bsp_mod(k)
  end do

  return

end subroutine factor_lambda

!*******************************************************************************
!
! Internal Subroutine:  gamma_sum
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine gamma_sum(gsum, m, nfft, order, kcut)

  implicit none

  double precision      gsum

  integer               m, nfft, order, kcut

  double precision      frac, x

  integer               k

  if (m .eq. 0) then
    gsum = 1.d0
    return
  end if

  frac = dble(m)/nfft 
  x = PI * frac
  gsum = 1.d0

  do k = 1, kcut
    gsum = gsum + (x/(x + PI * k))**order
  end do 

  do k = 1, kcut
    gsum = gsum + (x/(x - PI * k))**order
  end do 

  return

end subroutine gamma_sum

end subroutine pme_recip_dat_setup

!*******************************************************************************
!
! Subroutine:  allocate_m_exp_tbls
!
! Description: <TBS>
!
!*******************************************************************************

subroutine allocate_m_exp_tbls

  use mdin_ewald_dat_mod
  use pmemd_lib_mod

  implicit none

! Local variables:

  integer               :: alloc_failed

  allocate(m1_exp_tbl(-(nfft1/2 + 1) : nfft1/2 + 1), &
           m2_exp_tbl(-(nfft2/2 + 1) : nfft2/2 + 1), &
           m3_exp_tbl(-(nfft3/2 + 1) : nfft3/2 + 1), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  return

end subroutine allocate_m_exp_tbls

!*******************************************************************************
!
! Subroutine:  load_m_exp_tbls
!
! Description: <TBS>
!
!*******************************************************************************

subroutine load_m_exp_tbls

  use pmemd_lib_mod
  use mdin_ewald_dat_mod
  use pbc_mod

  implicit none

! Local variables:

  double precision      :: fac
  double precision      :: recip_11_2, recip_22_2, recip_33_2
  integer               :: i

  fac = -(PI * PI)/(ew_coeff * ew_coeff)

  recip_11_2 = recip(1,1) * recip(1,1)
  recip_22_2 = recip(2,2) * recip(2,2)
  recip_33_2 = recip(3,3) * recip(3,3)

  do i = -(nfft1/2 + 1), nfft1/2 + 1
    m1_exp_tbl(i) = exp(fac * dble(i) * dble(i) * recip_11_2)
  end do

  do i = -(nfft2/2 + 1), nfft2/2 + 1
    m2_exp_tbl(i) = exp(fac * dble(i) * dble(i) * recip_22_2)
  end do

  do i = -(nfft3/2 + 1), nfft3/2 + 1
    m3_exp_tbl(i) = exp(fac * dble(i) * dble(i) * recip_33_2)
  end do

  return

end subroutine load_m_exp_tbls

!*******************************************************************************!
! Subroutine:  zero_q_array
!
! Description: <TBS>
!
!*******************************************************************************
subroutine zero_q_array(array, num)

  implicit none
  double precision      array(*)
  integer               num

  integer               i

  do i = 1, num
    array(i) = 0.d0
  end do

  return

end subroutine zero_q_array

!*******************************************************************************!
! Subroutine:  zero_q_array_midpoint
!
! Description: <TBS>
!
!*******************************************************************************
subroutine zero_q_array_midpoint(array, num)

  implicit none
#ifdef pmemd_SPDP
  pme_float      array(*)
#else /*pmemd_SPDP*/
  double precision      array(*)
#endif /*pmemd_SPDP2*/
  integer               num

  integer               i

  do i = 1, num
    array(i) = 0.d0
  end do

  return

end subroutine zero_q_array_midpoint

end module pme_recip_dat_mod
