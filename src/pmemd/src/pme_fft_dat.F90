#include "copyright.i"

!*******************************************************************************
!
! Module: pme_fft_dat_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module pme_fft_dat_mod

  use fft1d_mod

  implicit none

! None of this stuff is broadcast; fft setup occurs in all tasks.
! The data here is common between different pme fft implementations -
! either slab or block.

  ! alpha/beta coefficients
  double precision, allocatable, save   :: fftrc_coefs(:,:)

  ! Handles to the 1d fft's data:

  type(fft1d_rec), pointer, save        :: fft_x_hdl => null()
  type(fft1d_rec), pointer, save        :: fft_y_hdl => null()
  type(fft1d_rec), pointer, save        :: fft_z_hdl => null()

  ! Physical dimensions of fft data.  The data is considered complex, but
  ! stored as double precision in the form:
  ! double precision data(2, fft_x_dim, fft_y_dim, fft_z_dim)
  ! The above index ordering is the default input/output form.
  ! Intermediate data  may be reordered due to transpositions.

  integer, save         :: fft_x_dim
  integer, save         :: fft_y_dim
  integer, save         :: fft_z_dim

#ifdef MPI
  integer, save                 :: siz_fft_mpi_buf
#endif

contains

!*******************************************************************************
!
! Subroutine:  fft_dat_setup
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine fft_dat_setup(num_ints, num_reals)

  use gbl_constants_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed
  integer               :: i
  double precision      :: pi2n
  double precision      :: theta

  ! Set fft dimensions - these are physical array dimensions:

  if (mod(nfft1, 2) .ne. 0) then
    write(mdout, *)'RCFFT requires that nfft1 is even!'
    call mexit(6, 1)
  end if

  fft_x_dim = (nfft1 + 2) / 2
  fft_y_dim = nfft2
  fft_z_dim = nfft3

  ! Create 1 dimensional fft handles, used to generically access all
  ! 1d fft implementations.

  call fft1d_create(fft_x_hdl, nfft1 / 2, num_ints, num_reals)
  call fft1d_create(fft_y_hdl, nfft2, num_ints, num_reals)
  call fft1d_create(fft_z_hdl, nfft3, num_ints, num_reals)

  ! Allocate and initialize arrays used in rc-cr fft conversions:

  allocate(fftrc_coefs(2, nfft1 / 2), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_reals = num_reals + size(fftrc_coefs)

  pi2n = 2 * PI / nfft1

  theta = 0.d0
  do i = 1, nfft1 / 2
    fftrc_coefs(1, i) = cos(theta)
    fftrc_coefs(2, i) = sin(theta)
    theta = pi2n * i
  end do

  return

end subroutine fft_dat_setup

end module pme_fft_dat_mod
