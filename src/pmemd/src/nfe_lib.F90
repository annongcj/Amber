! Pre-processing
#ifndef NFE_UTILS_H
#define NFE_UTILS_H

#ifndef NFE_DISABLE_ASSERT
#  define nfe_assert(stmt) if (.not.(stmt)) call afailed(__FILE__, __LINE__)
#  define nfe_assert_not_reached() call afailed(__FILE__, __LINE__)
#  define NFE_PURE_EXCEPT_ASSERT
#  define NFE_USE_AFAILED 
#else
#  define nfe_assert(s) 
#  define nfe_assert_not_reached() 
#  define NFE_PURE_EXCEPT_ASSERT pure
#  define NFE_USE_AFAILED
#endif /* NFE_DISABLE_ASSERT */

#define NFE_OUT_OF_MEMORY call out_of_memory(__FILE__, __LINE__)

#ifdef MPI
#  define NFE_MASTER_ONLY_BEGIN if (mytaskid.eq.0) then
#  define NFE_MASTER_ONLY_END end if
#else
#  define NFE_MASTER_ONLY_BEGIN
#  define NFE_MASTER_ONLY_END
#endif /* MPI */

#define NFE_ERROR   ' ** NFE-Error ** : '
#define NFE_WARNING ' ** NFE-Warning ** : '
#define NFE_INFO    ' NFE : '

#endif /* NFE_UTILS_H */
!----------------------------------------------------------------

module nfe_lib_mod

use file_io_dat_mod
use mdin_ctrl_dat_mod
use pmemd_lib_mod, only : mexit
use prmtop_dat_mod
use parallel_dat_mod
#ifdef BINTRAJ
use netcdf, only : nf90_double, nf90_float
#endif /* BINTRAJ */

implicit none

private

! Constant part--------------------------------------------------
double precision, public, parameter :: ZERO  = 0.d0
double precision, public, parameter :: ONE   = 1.d0
double precision, public, parameter :: TWO   = 2.d0
double precision, public, parameter :: THREE = 3.d0
double precision, public, parameter :: FOUR  = 4.d0
double precision, public, parameter :: kB = 1.9872041d-3  !Boltzmann's constant in (kcal/mol)/K
double precision, public, parameter :: TINY = 0.00000100000000000000D0 ! dble(0.000001)

integer, public, parameter :: STRING_LENGTH = 256

integer, public, parameter :: ERR_UNIT = mdout
integer, public, parameter :: OUT_UNIT = mdout
integer, public, parameter :: LAST_UNIT = 77

integer, public, parameter :: ABMD_MONITOR_UNIT = LAST_UNIT + 1
integer, public, parameter :: SMD_OUTPUT_UNIT = LAST_UNIT + 2
integer, public, parameter :: PMD_OUTPUT_UNIT = LAST_UNIT + 3
integer, public, parameter :: STSM_OUTPUT_UNIT = LAST_UNIT + 4

#ifdef MPI
integer, public, parameter :: REM_MDIN_UNIT = LAST_UNIT + 5
integer, public, parameter :: PMD_REMLOG_UNIT = LAST_UNIT + 6
integer, public, parameter :: ABMD_REMLOG_UNIT = LAST_UNIT + 7
integer, public, parameter :: BBMD_MONITOR_UNIT = LAST_UNIT + 8
integer, public, parameter :: BBMD_LOG_UNIT = LAST_UNIT + 9
integer, public, parameter :: STSM_REMLOG_UNIT = LAST_UNIT + 10
#endif /* MPI */

integer, public, parameter :: EVEC_UNIT1 = LAST_UNIT + 11
integer, public, parameter :: CRD_UNIT1  = LAST_UNIT + 12
integer, public, parameter :: REF_UNIT1 = LAST_UNIT + 13
integer, public, parameter :: IDX_UNIT1 = LAST_UNIT + 14

integer, public, parameter :: ABMD_CV_UNIT = LAST_UNIT + 15
integer, public, parameter :: SMD_CV_UNIT = LAST_UNIT + 16
integer, public, parameter :: PMD_CV_UNIT = LAST_UNIT + 17
integer, public, parameter :: BBMD_CV_UNIT = LAST_UNIT + 18
integer, public, parameter :: STSM_CV_UNIT = LAST_UNIT + 19
! ---------------------------------------------------------------

! PMEMD-proxy part-----------------------------------------------
public :: pmemd_mdin_name
public :: pmemd_mdin_unit
public :: pmemd_mdout_name

public :: multipmemd_rem
public :: multipmemd_initremd
public :: multipmemd_numgroup

public :: terminate
public :: is_master

public :: proxy_finalize

public :: pmemd_imin
public :: pmemd_natoms
public :: pmemd_mdtime
public :: pmemd_sgft
public :: pmemd_sgff
public :: pmemd_temp0
public :: pmemd_timestep
public :: pmemd_init
public :: pmemd_nstlim
public :: pmemd_ntp
public :: pmemd_ntb
public :: pmemd_nsolut 

#ifdef MPI
public :: set_pmemd_temp0
#endif

public :: pmemd_atom_name
public :: remember_atom_names

character(len = 4), private, pointer, save :: atom_names(:) => null()

public :: flush_UNIT

public :: remember_rem
public :: remember_initremd

public :: nfe_prt

integer, private, save :: saved_rem = -3212341
#ifdef MPI
logical, private, save :: saved_initremd = .true.
#else
logical, private, save :: saved_initremd = .false.
#endif /* MPI */
integer, public, save :: nfe_init
logical, public, save :: nfe_real_mdstep = .true.

! NFE potential energy contribution, added by F. Pan
type nfe_pot_ene_rec
    sequence
    double precision    :: smd
    double precision    :: pmd
    double precision    :: abmd
    double precision    :: bbmd
    double precision    :: stsm
    double precision    :: total  ! restraint energy in NFE module
end type nfe_pot_ene_rec

integer, parameter, public    :: nfe_pot_ene_rec_size = 6

type(nfe_pot_ene_rec), parameter, public      :: null_nfe_pot_ene_rec = &
    nfe_pot_ene_rec(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)

type(nfe_pot_ene_rec), public, save :: nfe_pot_ene = null_nfe_pot_ene_rec

! atom mask
integer, public, save :: nfe_atm_cnt = 0
integer, public, allocatable, save :: nfe_atm_lst(:)
logical, public, save :: nfe_first = .true.

!---------------------- U T I L S ------------------------------- 
public :: fatal
public :: out_of_memory

public :: close_UNIT

#ifdef MPI
public :: cpus_enter
public :: cpus_leave
#endif /* MPI */

#ifndef NFE_DISABLE_ASSERT
public :: afailed
#endif /* NFE_DISABLE_ASSERT */

interface swap
   module procedure swap_i, swap_r4, swap_r8, swap_cSL
end interface swap

private :: swap_i, swap_r4, swap_r8, swap_cSL

public :: swap

interface pfmt
   module procedure i0_format_1, f0_format_1
end interface

private :: i0_format_1, f0_format_1, n_digits_i, n_digits_r

public :: pfmt ! this is needed because "I0" and "F0.X" editing is not portable

public :: AddToList

!============================================================================+
!             + x +   U M B R E L L A    I N T E R F A C E   + x +           !
!============================================================================+

integer, public, parameter :: UMBRELLA_MIN_EXTENT   = 5
integer, public, parameter :: UMBRELLA_MAX_NEXTENTS = 4

type, public :: umbrella_t

   private

   integer :: nextents, extents(UMBRELLA_MAX_NEXTENTS)
   logical :: periodicity(UMBRELLA_MAX_NEXTENTS)

   double precision :: origin(UMBRELLA_MAX_NEXTENTS)
   double precision :: spacing(UMBRELLA_MAX_NEXTENTS)

   double precision, pointer :: coeffs(:) ! row-major

#ifndef NFE_DISABLE_ASSERT
   logical :: inited = .false.
#endif /* NFE_DISABLE_ASSERT */
end type umbrella_t

public :: umbrella_init
public :: umbrella_fini

public :: umbrella_nextents
public :: umbrella_extent
public :: umbrella_origin
public :: umbrella_spacing
public :: umbrella_periodicity

public :: umbrella_hill
public :: umbrella_swap
public :: umbrella_transfer

public :: umbrella_eval_v
public :: umbrella_eval_vdv
! Modified by M Moradi
public :: umbrella_eval_laplacian
! Moradi end

! By F Pan
public :: umbrella_wt_mod
public :: umbrella_copy
! Pan end

#ifdef MPI
public :: umbrella_bcast
public :: umbrella_send_coeffs
public :: umbrella_recv_coeffs
#endif /* MPI */

#ifdef BINTRAJ
public :: umbrella_load
public :: umbrella_save
integer, private, parameter :: coeffs_type = nf90_double
#endif /* BINTRAJ */

!============================================================================+
!                      + x +   D E T A I L S   + x +                         !
!============================================================================+

private :: m4_v
private :: m4_vdv

double precision, private, parameter :: MINUS_ONE = -ONE
double precision, private, parameter :: ONE_THIRD = ONE/THREE
double precision, private, parameter :: MINUS_ONE_THIRD = -ONE_THIRD
double precision, private, parameter :: ONE_SIXTH = ONE/6
double precision, private, parameter :: TWO_THIRD = TWO/THREE
double precision, private, parameter :: HALF = ONE/TWO
double precision, private, parameter :: MINUS_HALF = -HALF
! ---------------------------------------------------------------

contains

! ---------------------------------------------------------------

subroutine terminate()
   implicit none
   call mexit(6, 1)
end subroutine terminate

!-----------------------------------------------------------------------------

character(len=STRING_LENGTH) function pmemd_mdin_name()
   implicit none
   pmemd_mdin_name = mdin_name
end function pmemd_mdin_name

!-----------------------------------------------------------------------------
character(len=STRING_LENGTH) function pmemd_mdout_name()
   implicit none
   pmemd_mdout_name = mdout_name
end function pmemd_mdout_name

!-----------------------------------------------------------------------------

pure integer function pmemd_mdin_unit()
   implicit none
   pmemd_mdin_unit = 5
end function pmemd_mdin_unit

!-----------------------------------------------------------------------------

pure integer function multipmemd_numgroup()
   implicit none
   multipmemd_numgroup = numgroups
end function multipmemd_numgroup

!-----------------------------------------------------------------------------

NFE_PURE_EXCEPT_ASSERT integer function multipmemd_rem()
   implicit none
   multipmemd_rem = saved_rem
end function multipmemd_rem

!-----------------------------------------------------------------------------

subroutine remember_rem(r)
   implicit none
   integer, intent(in) :: r
   saved_rem = r
end subroutine remember_rem

!-----------------------------------------------------------------------------

subroutine remember_initremd(i)
   implicit none
   logical, intent(in) :: i
   saved_initremd = i
end subroutine remember_initremd

!-----------------------------------------------------------------------------

pure logical function multipmemd_initremd()
   implicit none
   multipmemd_initremd = saved_initremd
end function multipmemd_initremd

!-----------------------------------------------------------------------------

subroutine nfe_prt(i)
  implicit none
  integer, intent(in) :: i

  if (infe .eq. 0) return

  write(i, 100) nfe_pot_ene%smd, nfe_pot_ene%pmd, nfe_pot_ene%abmd
  write(i, 101) nfe_pot_ene%bbmd, nfe_pot_ene%stsm
  write(i, 102)
  return

100 format(' NFE restraints:    SMD  :',f9.3,4x,'PMD  : ',f9.3,4x, &
          'ABMD : ',f9.3)
101 format(20x,'BBMD :',f9.3,4x,'STSM : ',f9.3)
102 format(79('='))  
end subroutine nfe_prt

!-----------------------------------------------------------------------------
!
! for use in nfe_assert() & co [where unneeded indirection is acceptable]
!

logical function is_master()

   implicit none

#ifdef MPI
   nfe_assert(pmemd_comm /= mpi_comm_null)
   is_master = (mytaskid == 0)
#else
   is_master = .true.
#endif /* MPI */

end function is_master

!-----------------------------------------------------------------------------

subroutine proxy_finalize()
   implicit none
   if (associated(atom_names)) &
      deallocate(atom_names)
end subroutine proxy_finalize

!-----------------------------------------------------------------------------

pure integer function pmemd_imin()
   implicit none
   pmemd_imin = imin
end function pmemd_imin

!-----------------------------------------------------------------------------
pure integer function pmemd_nsolut()
   use shake_mod,only : ibgwat
   implicit none
   pmemd_nsolut = natom-(nres-ibgwat+1)*4
end function pmemd_nsolut
!-----------------------------------------------------------------------------

pure integer function pmemd_natoms()
   implicit none
   pmemd_natoms = natom
end function pmemd_natoms

!-----------------------------------------------------------------------------

pure double precision function pmemd_mdtime()
   implicit none
   pmemd_mdtime = t
end function pmemd_mdtime

!-----------------------------------------------------------------------------

pure double precision function pmemd_sgft()
   implicit none
   pmemd_sgft = sgft
end function pmemd_sgft

!-----------------------------------------------------------------------------

pure double precision function pmemd_sgff()
   implicit none
   pmemd_sgff = sgff
end function pmemd_sgff

!-----------------------------------------------------------------------------

pure double precision function pmemd_temp0()
   implicit none
   pmemd_temp0 = temp0
end function pmemd_temp0

!-----------------------------------------------------------------------------

#ifdef MPI
subroutine set_pmemd_temp0(new_temp0)
   implicit none
   double precision, intent(in) :: new_temp0
   temp0 = new_temp0
end subroutine set_pmemd_temp0
#endif /* MPI */

!-----------------------------------------------------------------------------

pure double precision function pmemd_timestep()
   implicit none
   pmemd_timestep = dt
end function pmemd_timestep

!-----------------------------------------------------------------------------

pure integer function pmemd_init()
   implicit none
   pmemd_init = nfe_init
end function pmemd_init

!-----------------------------------------------------------------------------

pure integer function pmemd_nstlim()
   implicit none
   pmemd_nstlim = nstlim
end function pmemd_nstlim

!-----------------------------------------------------------------------------

pure integer function pmemd_ntp()
   implicit none
   pmemd_ntp = ntp
end function pmemd_ntp

!-----------------------------------------------------------------------------

pure integer function pmemd_ntb()
   implicit none
   pmemd_ntb = ntb
end function pmemd_ntb

!-----------------------------------------------------------------------------

character(len = 4) function pmemd_atom_name(n)

   implicit none

   integer, intent(in) :: n

   nfe_assert(n > 0)
   nfe_assert(n <= pmemd_natoms())
   nfe_assert(associated(atom_names))

   pmemd_atom_name = atom_names(n)

end function pmemd_atom_name

!-----------------------------------------------------------------------------

subroutine remember_atom_names(ih)

   implicit none

   character(len = 4), intent(in) :: ih(*)

   integer :: n, error

   if (associated(atom_names)) &
      deallocate(atom_names)

   nfe_assert(natom > 0)

   allocate(atom_names(natom), stat = error)
   if (error /= 0) then
      write (unit = ERR_UNIT, fmt = '(a,a)') &
         NFE_ERROR, 'out of memory in remember_atom_names()'
      call terminate()
   end if

   do n = 1, natom
      atom_names(n) = ih(n)
   end do

end subroutine remember_atom_names

!-----------------------------------------------------------------------------

subroutine flush_UNIT(lun,filename)

   implicit none

   integer, intent(in) :: lun
   character(len = *), intent(in) :: filename

! from runfiles.F90  srb july 2016:

! Flushing actually does not work particularly reliably for a number of
! machines and compilers, and in the more benign cases simply fails, but in
! the more malign cases can actually corrupt the stack (due to a compiler-
! dependent flush() call interface change).  We therefore no longer do
! flushes of anything in PMEMD; if it needs to go out, we close it and reopen
! it.

! thus, this is commented out until the above close/reopen is implemented in
! order to fix the broken compilation of pmemd on mic's  srb july 2016:
!   call amflsh(lun)
   close(lun)
   open(unit=lun, file=filename, status='OLD', position='APPEND')

end subroutine flush_UNIT

!-----------------------------------------------------------------------------

subroutine fatal(message)

   implicit none

   character(len = *), intent (in) :: message

   if (is_master()) &
      write(unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR, message

   call terminate()

end subroutine fatal

!-----------------------------------------------------------------------------

subroutine out_of_memory(filename, lineno)

   implicit none

   character(len = *), intent(in) :: filename
   integer,            intent(in) :: lineno

   write (unit = ERR_UNIT, fmt = '(a,a,a,a,'//pfmt(lineno)//')') &
      NFE_ERROR, 'memory allocation failed at ', filename, ':', lineno

   call terminate()

end subroutine out_of_memory

!-----------------------------------------------------------------------------

subroutine close_UNIT(u)

   implicit none

   integer, intent(in) :: u
   logical :: o

   inquire (unit = u, opened = o)
   if (o) close (unit = u)

end subroutine close_UNIT

!-----------------------------------------------------------------------------

#ifndef NFE_DISABLE_ASSERT
subroutine afailed(filename, lineno)

   implicit none

   character(len = *), intent(in) :: filename
   integer,            intent(in) :: lineno

   write(unit = ERR_UNIT, fmt = '(/a,a,a,'//pfmt(lineno)//',a/)') &
         NFE_ERROR, filename, ':', lineno, ': nfe_assert() failed'
   call flush_UNIT(ERR_UNIT,filename)
   call terminate()

end subroutine afailed
#endif /* NFE_DISABLE_ASSERT */

!-----------------------------------------------------------------------------

#ifdef MPI

!
! stolen from mpb-1.4.2
!

subroutine cpus_enter(comm, tag)

   implicit none

   integer, intent(in) :: comm, tag

   integer :: commrank, commsize
   integer :: recv_tag, recv_status(MPI_STATUS_SIZE), error

   call mpi_comm_rank(comm, commrank, error)
   nfe_assert(error.eq.0)

   call mpi_comm_size(comm, commsize, error)
   nfe_assert(error.eq.0)

   nfe_assert(commrank.ge.0)
   nfe_assert(commsize.gt.0)
   nfe_assert(commrank.lt.commsize)

   if (commrank.gt.0) then
      recv_tag = tag - 1
      call mpi_recv(recv_tag, 1, MPI_INTEGER, &
         commrank - 1, tag, comm, recv_status, error)
      nfe_assert(error.eq.0)
      nfe_assert(recv_tag.eq.tag)
   end if

end subroutine cpus_enter

!-----------------------------------------------------------------------------

subroutine cpus_leave(comm, tag)

   implicit none

   integer, intent(in) :: comm, tag

   integer :: error, commsize, commrank

   call mpi_comm_rank(comm, commrank, error)
   nfe_assert(error.eq.0)

   call mpi_comm_size(comm, commsize, error)
   nfe_assert(error.eq.0)

   nfe_assert(commrank.ge.0)
   nfe_assert(commsize.gt.0)
   nfe_assert(commrank.lt.commsize)

   if (commrank.ne.(commsize - 1)) then
      call mpi_send(tag, 1, MPI_INTEGER, &
         commrank + 1, tag, comm, error)
      nfe_assert(error.eq.0)
   end if

end subroutine cpus_leave
#endif /* MPI */

!-----------------------------------------------------------------------------

subroutine swap_i(a, b)

   implicit none

   integer, intent(inout) :: a, b

   integer :: tmp

   tmp = a
   a = b
   b = tmp

end subroutine swap_i

!-----------------------------------------------------------------------------

subroutine swap_r4(a, b)

   implicit none

   real(4), intent(inout) :: a, b

   real(4) :: tmp

   tmp = a
   a = b
   b = tmp

end subroutine swap_r4

!-----------------------------------------------------------------------------

subroutine swap_r8(a, b)

   implicit none

   real(8), intent(inout) :: a, b

   real(8) :: tmp

   tmp = a
   a = b
   b = tmp

end subroutine swap_r8

!-----------------------------------------------------------------------------

subroutine swap_cSL(a, b)

   implicit none

   character(STRING_LENGTH), intent(inout) :: a, b
   
   character(STRING_LENGTH) :: tmp

   tmp = a
   a = b
   b = tmp

end subroutine swap_cSL

!-----------------------------------------------------------------------------

function i0_format_1(x) result(f)

   implicit none

   character(8) :: f
   integer, intent(in) :: x

   integer :: n

   n = n_digits_i(x)

   if (n.le.9) then
      write (unit = f, fmt = '(a,i1)') 'i', n
   else
      write (unit = f, fmt = '(a,i2)') 'i', n
   end if

end function i0_format_1

!-----------------------------------------------------------------------------

function f0_format_1(x, y) result(f)

   implicit none

   character(8) :: f
   real(8), intent(in) :: x
   integer, intent(in) :: y

   integer :: n

   n = n_digits_r(x)

   if ((n + y + 1).le.9) then
      write (unit = f, fmt = '(a,i1,a,i1)') 'f', (n + y + 1), '.', y
   else
      write (unit = f, fmt = '(a,i2,a,i1)') 'f', (n + y + 1), '.', y
   end if

end function f0_format_1

!-----------------------------------------------------------------------------

pure function n_digits_i(i) result(n)

   implicit none

   integer :: n

   integer, intent(in)  :: i

   n = n_digits_r(dble(i))

end function n_digits_i

!-----------------------------------------------------------------------------

pure function n_digits_r(r) result(n)

   implicit none

   integer :: n

   real(8), intent(in)  :: r

   if (r.gt.1.0D0) then
      n = 1 + int(floor(log10(r)))
   else if (r.ge.0.0D0) then
      n = 1
   else if (r.ge.-1.0D0) then
      n = 2
   else
      n = 2 + int(floor(log10(-r)))
   end if

end function n_digits_r

!-----------------------------------------------------------------------------
subroutine AddToList(list, element)

   implicit none

   integer :: i, isize
   integer, intent(in) :: element
   integer, dimension(:), allocatable, intent(inout) :: list
   integer, dimension(:), allocatable :: clist


   if(allocated(list)) then
        isize = size(list)
        allocate(clist(isize+1))
        do i=1,isize
        clist(i) = list(i)
        end do
        clist(isize+1) = element

        deallocate(list)
        call move_alloc(clist, list)

   else
        allocate(list(1))
        list(1) = element
   end if

end subroutine AddToList

!========================================================================
subroutine umbrella_init(this, nextents, extents, origin, spacing, periodicity)

   implicit none

   type(umbrella_t), intent(inout) :: this

   integer, intent(in) :: nextents
   integer, intent(in) :: extents(*)

   double precision, intent(in) :: origin(*)
   double precision, intent(in) :: spacing(*)

   logical, intent(in) :: periodicity(*)

   integer :: n, ncoeffs, ierr, tmp

   nfe_assert(.not.this%inited)

   nfe_assert(nextents.ge.1)
   nfe_assert(nextents.le.UMBRELLA_MAX_NEXTENTS)

   this%nextents = nextents

   ncoeffs = 1
   do n = 1, nextents
      nfe_assert(extents(n).ge.UMBRELLA_MIN_EXTENT)
      nfe_assert(spacing(n).gt.ZERO)

      this%origin(n) = origin(n)
      this%spacing(n) = spacing(n)

      this%extents(n) = extents(n)
      this%periodicity(n) = periodicity(n)

      tmp = ncoeffs
      ncoeffs = ncoeffs*extents(n)
      if (ncoeffs/extents(n).ne.tmp) &
         call fatal('overflow in umbrella_init()')
   end do

   allocate(this%coeffs(ncoeffs), stat = ierr)
   if (ierr /= 0) &
      NFE_OUT_OF_MEMORY

   this%coeffs(1:ncoeffs) = ZERO

#ifndef NFE_DISABLE_ASSERT
   this%inited = .true.
#endif /* NFE_DISABLE_ASSERT */

end subroutine umbrella_init

!=============================================================================

subroutine umbrella_fini(this)

   NFE_USE_AFAILED

   implicit none

   type(umbrella_t), intent(inout) :: this

   nfe_assert(this%inited)
   nfe_assert(associated(this%coeffs))

   deallocate(this%coeffs)

#ifndef NFE_DISABLE_ASSERT
   this%inited = .false.
#endif /* NFE_DISABLE_ASSERT */

end subroutine umbrella_fini

!=============================================================================

integer function umbrella_nextents(this)

   NFE_USE_AFAILED

   implicit none

   type(umbrella_t), intent(in) :: this

   nfe_assert(this%inited)
   umbrella_nextents = this%nextents

end function umbrella_nextents

!=============================================================================

integer function umbrella_extent(this, n)

   NFE_USE_AFAILED

   implicit none

   type(umbrella_t), intent(in) :: this
   integer, intent(in):: n

   nfe_assert(this%inited)
   nfe_assert(n.ge.1)
   nfe_assert(n.le.this%nextents)

   umbrella_extent = this%extents(n)

end function umbrella_extent

!=============================================================================

double precision function umbrella_origin(this, n)

   NFE_USE_AFAILED

   implicit none

   type(umbrella_t), intent(in) :: this
   integer, intent(in):: n

   nfe_assert(this%inited)
   nfe_assert(n.ge.1)
   nfe_assert(n.le.this%nextents)

   umbrella_origin = this%origin(n)

end function umbrella_origin

!=============================================================================

double precision function umbrella_spacing(this, n)

   NFE_USE_AFAILED

   implicit none

   type(umbrella_t), intent(in) :: this
   integer, intent(in):: n

   nfe_assert(this%inited)
   nfe_assert(n.ge.1)
   nfe_assert(n.le.this%nextents)

   umbrella_spacing = this%spacing(n)

end function umbrella_spacing

!=============================================================================

logical function umbrella_periodicity(this, n)

   NFE_USE_AFAILED

   implicit none

   type(umbrella_t), intent(in) :: this
   integer, intent(in):: n

   nfe_assert(this%inited)
   nfe_assert(n.ge.1)
   nfe_assert(n.le.this%nextents)

   umbrella_periodicity = this%periodicity(n)

end function umbrella_periodicity

!=============================================================================

subroutine umbrella_hill(this, x, alt)

   NFE_USE_AFAILED

   implicit none

   type(umbrella_t), intent(inout) :: this

   double precision, intent(in) :: x(*)
   double precision, intent(in) :: alt

   double precision :: xs, accum, hill_v(4, UMBRELLA_MAX_NEXTENTS)
   integer :: gc(UMBRELLA_MAX_NEXTENTS), i, j, n, p, o

   nfe_assert(this%inited)
   nfe_assert(this%nextents.ge.1)
   nfe_assert(this%nextents.le.UMBRELLA_MAX_NEXTENTS)
   nfe_assert(associated(this%coeffs))

   do i = 1, this%nextents
      xs = (x(i) - this%origin(i))/this%spacing(i)
      gc(i) = int(floor(xs))

      do j = 1, 4
         hill_v(j, i) = hill(xs - dble(j + gc(i) - 2))
      end do
   end do

   ! loop over 4 x 4 x ... x 4

   outer: do n = 1, ishft(1, ishft(this%nextents, 1))

      p = n
      o = 0

      accum = alt

      do i = 1, this%nextents

         j = gc(i) + iand(p, 3) - 1

         if (j.lt.0.or.j.ge.this%extents(i)) then
            if (this%periodicity(i)) then
               if (j.lt.0) then
                  j = this%extents(i) - 1 + mod(j + 1, this%extents(i))
               else
                  j = mod(j, this%extents(i))
               end if
            else
               cycle outer
            end if
         end if

         if (i.eq.this%nextents) then
            o = o + j + 1
         else
            o = o + this%extents(i + 1)*(o + j)
         end if

         accum = accum*hill_v(iand(p, 3) + 1, i)

         p = ishft(p, -2)

      end do ! loop over i

      this%coeffs(o) = this%coeffs(o) + accum

   end do outer

!-----------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------------

! "hill" centered at 0 (-2 < x < 2)
pure double precision function hill(x)

   implicit none

   double precision, intent(in) :: x

   double precision, parameter :: SCALE = 48.000000000000000000000D0 / 41.000000000000000000000D0!dble(48)/dble(41)

   double precision :: x2

   x2 = (HALF*x)**2

!   if (x2.lt.ONE) then
    hill = SCALE*(x2 - ONE)**2
!   else
!      hill = ZERO
!   end if

end function hill

!-----------------------------------------------------------------------------

end subroutine umbrella_hill

!=============================================================================

! sets 'this' from 'other' [does *NOT* interpolate rigorously]
subroutine umbrella_transfer(this, other)

   NFE_USE_AFAILED

   implicit none

   type(umbrella_t), intent(inout) :: this
   type(umbrella_t), intent(in) :: other

   double precision :: pos(UMBRELLA_MAX_NEXTENTS)

   integer :: i, j, p, n, o, ncoeffs

   nfe_assert(this%inited)
   nfe_assert(other%inited)

   nfe_assert(this%nextents.ge.1)
   nfe_assert(this%nextents.le.UMBRELLA_MAX_NEXTENTS)

   nfe_assert(other%nextents.ge.1)
   nfe_assert(other%nextents.le.UMBRELLA_MAX_NEXTENTS)

   nfe_assert(this%nextents.eq.other%nextents)

   ncoeffs = 1
   do n = 1, this%nextents
      ncoeffs = ncoeffs*this%extents(n)
   end do

   ! this is *not* interpolation
   do n = 1, ncoeffs
      p = n
      o = 0
      do i = 1, this%nextents
         j = mod(p, this%extents(i))
         pos(i) = this%origin(i) + j*this%spacing(i)
         if (i.eq.this%nextents) then
            o = o + j + 1
         else
            o = o + this%extents(i + 1)*(o + j)
         end if
         p = p/this%extents(i)
      end do
      this%coeffs(o) = umbrella_eval_v(other, pos)
   end do

end subroutine umbrella_transfer

!=============================================================================

subroutine umbrella_copy(this, other)
   NFE_USE_AFAILED

   implicit none

   type(umbrella_t), intent(inout) :: this
   type(umbrella_t), intent(in) :: other

   integer :: n, ncoeffs
   double precision :: tmp

   nfe_assert(this%inited)
   nfe_assert(other%inited)

   nfe_assert(this%nextents.ge.1)
   nfe_assert(this%nextents.le.UMBRELLA_MAX_NEXTENTS)

   nfe_assert(other%nextents.ge.1)
   nfe_assert(other%nextents.le.UMBRELLA_MAX_NEXTENTS)

   nfe_assert(this%nextents.eq.other%nextents)

   ncoeffs = 1
   do n = 1, this%nextents
      ncoeffs = ncoeffs*this%extents(n)
   end do

   do n = 1, ncoeffs
      tmp = other%coeffs(n)
      this%coeffs(n) = tmp
   end do
end subroutine umbrella_copy

!=============================================================================

! swaps without any checks
subroutine umbrella_swap(this, other)

   implicit none

   type(umbrella_t), intent(inout) :: this, other

   type(umbrella_t) :: tmp

   tmp%nextents = this%nextents
   this%nextents = other%nextents
   other%nextents = tmp%nextents

   tmp%extents = this%extents
   this%extents = other%extents
   other%extents = tmp%extents

   tmp%periodicity = this%periodicity
   this%periodicity = other%periodicity
   other%periodicity = tmp%periodicity

   tmp%origin = this%origin
   this%origin = other%origin
   other%origin = tmp%origin

   tmp%spacing = this%spacing
   this%spacing = other%spacing
   other%spacing = tmp%spacing

#ifndef NFE_DISABLE_ASSERT
   tmp%inited = this%inited
   this%inited = other%inited
   other%inited = tmp%inited
#endif /* NFE_DISABLE_ASSERT */

   tmp%coeffs => this%coeffs
   this%coeffs => other%coeffs
   other%coeffs => tmp%coeffs

end subroutine umbrella_swap

!=============================================================================

double precision function umbrella_eval_v(this, x) result(v)

   NFE_USE_AFAILED

   implicit none

   type(umbrella_t), intent(in) :: this
   double precision, intent(in) :: x(*)

   double precision :: xs, m4v_prod, m4v(4, UMBRELLA_MAX_NEXTENTS)
   integer :: gc(UMBRELLA_MAX_NEXTENTS), i, j, n, p, o

   nfe_assert(this%inited)
   nfe_assert(this%nextents.ge.1)
   nfe_assert(this%nextents.le.UMBRELLA_MAX_NEXTENTS)
   nfe_assert(associated(this%coeffs))

   v = ZERO

   do i = 1, this%nextents
      xs = (x(i) - this%origin(i))/this%spacing(i)
      gc(i) = int(floor(xs))

      do j = 1, 4
         m4v(j, i) = m4_v(xs - dble(j + gc(i) - 2))
      end do
   end do

   ! loop over 4 x 4 x ... x 4

   outer: do n = 1, ishft(1, ishft(this%nextents, 1))

      p = n
      o = 0

      m4v_prod = ONE

      do i = 1, this%nextents

         j = gc(i) + iand(p, 3) - 1

         if (j.lt.0.or.j.ge.this%extents(i)) then
            if (this%periodicity(i)) then
               if (j.lt.0) then
                  j = this%extents(i) - 1 + mod(j + 1, this%extents(i))
               else
                  j = mod(j, this%extents(i))
               end if
            else
               cycle outer
            end if
         end if

         if (i.eq.this%nextents) then
            o = o + j + 1
         else
            o = o + this%extents(i + 1)*(o + j)
         end if

         m4v_prod = m4v_prod*m4v(iand(p, 3) + 1, i)

         p = ishft(p, -2)

      end do ! loop over i

      v = v + this%coeffs(o)*m4v_prod

   end do outer

end function umbrella_eval_v

!=============================================================================

subroutine umbrella_eval_vdv(this, x, v, dv)

   NFE_USE_AFAILED

   implicit none

   type(umbrella_t), intent(in) :: this

   double precision, intent(in)  :: x(*)
   double precision, intent(out) :: v
   double precision, intent(out) :: dv(*)

   double precision :: xs, m4v_prod, m4dv_prod(UMBRELLA_MAX_NEXTENTS)
   double precision :: m4v(4, UMBRELLA_MAX_NEXTENTS), m4dv(4, UMBRELLA_MAX_NEXTENTS)

   integer :: gc(UMBRELLA_MAX_NEXTENTS), i, j, n, p, o

   nfe_assert(this%inited)
   nfe_assert(this%nextents.ge.1)
   nfe_assert(this%nextents.le.UMBRELLA_MAX_NEXTENTS)
   nfe_assert(associated(this%coeffs))

   v = ZERO

   do i = 1, this%nextents
      dv(i) = ZERO

      xs = (x(i) - this%origin(i))/this%spacing(i)
      gc(i) = int(floor(xs))

      do j = 1, 4
         call m4_vdv(xs - dble(j + gc(i) - 2), m4v(j, i), m4dv(j, i))
      end do
   end do

   ! loop over 4 x 4 x ... x 4

   outer: do n = 1, ishft(1, ishft(this%nextents, 1))

      p = n
      o = 0

      m4v_prod = ONE

      do i = 1, this%nextents
         m4dv_prod(i) = ONE
      end do

      do i = 1, this%nextents

         j = gc(i) + iand(p, 3) - 1

         if (j.lt.0.or.j.ge.this%extents(i)) then
            if (this%periodicity(i)) then
               if (j.lt.0) then
                  j = this%extents(i) - 1 + mod(j + 1, this%extents(i))
               else
                  j = mod(j, this%extents(i))
               end if
            else
               cycle outer
            end if
         end if

         if (i.eq.this%nextents) then
            o = o + j + 1
         else
            o = o + this%extents(i + 1)*(o + j)
         end if

         m4v_prod = m4v_prod*m4v(iand(p, 3) + 1, i)

         do j = 1, i - 1
            m4dv_prod(j) = m4dv_prod(j)*m4v(iand(p, 3) + 1, i)
         end do

         m4dv_prod(i) = m4dv_prod(i)*m4dv(iand(p, 3) + 1, i)

         do j = i + 1, this%nextents
            m4dv_prod(j) = m4dv_prod(j)*m4v(iand(p, 3) + 1, i)
         end do

         p = ishft(p, -2)

      end do ! loop over i

      v = v + this%coeffs(o)*m4v_prod

      do i = 1, this%nextents
         dv(i) = dv(i) + this%coeffs(o)*m4dv_prod(i)
      end do

   end do outer

   do i = 1, this%nextents
      dv(i) = dv(i)/this%spacing(i)
   end do

end subroutine umbrella_eval_vdv

!=============================================================================
! Added by M Moradi
! Laplacian of umbrella potential

subroutine umbrella_eval_laplacian(this, x, v, d2v)

   NFE_USE_AFAILED

   implicit none

   type(umbrella_t), intent(in) :: this

   double precision, intent(in)  :: x(*)
   double precision, intent(out) :: v
   double precision, intent(out) :: d2v

   double precision :: xs, m4v_prod, m4d2v_prod(UMBRELLA_MAX_NEXTENTS)
   double precision :: m4v(4, UMBRELLA_MAX_NEXTENTS), m4d2v(4, UMBRELLA_MAX_NEXTENTS)

   integer :: gc(UMBRELLA_MAX_NEXTENTS), i, j, n, p, o

   nfe_assert(this%inited)
   nfe_assert(this%nextents.ge.1)
   nfe_assert(this%nextents.le.UMBRELLA_MAX_NEXTENTS)
   nfe_assert(associated(this%coeffs))

   v = ZERO
   d2v = ZERO

   do i = 1, this%nextents
      xs = (x(i) - this%origin(i))/this%spacing(i)
      gc(i) = int(floor(xs))
      do j = 1, 4
         call m4_vd2v(xs - dble(j + gc(i) - 2), m4v(j, i), m4d2v(j, i))
      end do
   end do

   ! loop over 4 x 4 x ... x 4

   outer: do n = 1, ishft(1, ishft(this%nextents, 1))

      p = n
      o = 0

      m4v_prod = ONE

      do i = 1, this%nextents
         m4d2v_prod(i) = ONE
      end do

      do i = 1, this%nextents

         j = gc(i) + iand(p, 3) - 1

         if (j.lt.0.or.j.ge.this%extents(i)) then
            if (this%periodicity(i)) then
               if (j.lt.0) then
                  j = this%extents(i) - 1 + mod(j + 1, this%extents(i))
               else
                  j = mod(j, this%extents(i))
               end if
            else
               cycle outer
            end if
         end if

         if (i.eq.this%nextents) then
            o = o + j + 1
         else
            o = o + this%extents(i + 1)*(o + j)
         end if

         m4v_prod = m4v_prod*m4v(iand(p, 3) + 1, i)

         do j = 1, i - 1
            m4d2v_prod(j) = m4d2v_prod(j)*m4v(iand(p, 3) + 1, i)
         end do

         m4d2v_prod(i) = m4d2v_prod(i)*m4d2v(iand(p, 3) + 1, i)

         do j = i + 1, this%nextents
            m4d2v_prod(j) = m4d2v_prod(j)*m4v(iand(p, 3) + 1, i)
         end do

         p = ishft(p, -2)

      end do ! loop over i

      v = v + this%coeffs(o)*m4v_prod

      do i = 1, this%nextents
         d2v = d2v + this%coeffs(o)*m4d2v_prod(i)/(this%spacing(i)**2)
      end do

   end do outer

end subroutine umbrella_eval_laplacian

!=============================================================================

#ifdef MPI

subroutine umbrella_bcast(this, comm, root) ! inefficient && memory-leak prone

   implicit none

   type(umbrella_t), intent(inout) :: this

   integer, intent(in) :: comm, root

   integer :: p(UMBRELLA_MAX_NEXTENTS)

   integer :: n, error, commrank, commsize, ncoeffs

   nfe_assert(comm.ne.mpi_comm_null)

   call mpi_comm_rank(comm, commrank, error)
   nfe_assert(error.eq.0)

   call mpi_comm_size(comm, commsize, error)
   nfe_assert(error.eq.0)

   nfe_assert(root.ge.0)
   nfe_assert(root.lt.commsize)

#  ifndef NFE_DISABLE_ASSERT
   if (commrank.eq.root) then
      nfe_assert(this%inited)
   else
      nfe_assert(.not.this%inited)
   end if
#  endif /* NFE_DISABLE_ASSERT */

   call mpi_bcast(this%nextents, 1, MPI_INTEGER, root, comm, error)
   nfe_assert(error.eq.0)

   call mpi_bcast(this%extents, this%nextents, MPI_INTEGER, root, comm, error)
   nfe_assert(error.eq.0)

   do n = 1, this%nextents
      if (this%periodicity(n)) then
         p(n) = 1
      else
         p(n) = 0
      end if
   end do

   call mpi_bcast(p, this%nextents, MPI_INTEGER, root, comm, error)
   nfe_assert(error.eq.0)

   do n = 1, this%nextents
      this%periodicity(n) = (p(n).eq.1)
   end do

   call mpi_bcast(this%origin, this%nextents, MPI_DOUBLE_PRECISION, &
      root, comm, error)
   nfe_assert(error.eq.0)

   call mpi_bcast(this%spacing, this%nextents, MPI_DOUBLE_PRECISION, &
      root, comm, error)
   nfe_assert(error.eq.0)

   ncoeffs = 1
   do n = 1, this%nextents
      error = ncoeffs
      ncoeffs = ncoeffs*this%extents(n)
      if (ncoeffs/this%extents(n).ne.error) &
         call fatal('overflow in umbrella_bcast()')
   end do

   if (commrank.ne.root) then
      allocate(this%coeffs(ncoeffs), stat = error)
      if (error.ne.0) &
         NFE_OUT_OF_MEMORY
   end if

   call mpi_bcast(this%coeffs, ncoeffs, MPI_DOUBLE_PRECISION, root, comm, error)
   nfe_assert(error.eq.0)

#  ifndef NFE_DISABLE_ASSERT
   this%inited = .true.
#  endif /* NFE_DISABLE_ASSERT */

end subroutine umbrella_bcast

!=============================================================================

subroutine umbrella_send_coeffs(this, dst, comm) ! unsafe

   implicit none

   type(umbrella_t), intent(inout) :: this

   integer, intent(in) :: dst, comm

   integer :: n, ncoeffs, error

   nfe_assert(this%inited)

   ncoeffs = 1
   do n = 1, this%nextents
      error = ncoeffs
      ncoeffs = ncoeffs*this%extents(n)
      if (ncoeffs/this%extents(n).ne.error) &
         call fatal('overflow in umbrella_send_coeffs()')
   end do

   call mpi_send(this%coeffs, ncoeffs, MPI_DOUBLE_PRECISION, &
                 dst, this%nextents, comm, error)

   nfe_assert(error.eq.0)

end subroutine umbrella_send_coeffs

!=============================================================================

subroutine umbrella_recv_coeffs(this, src, comm) ! unsafe

   implicit none

   type(umbrella_t), intent(inout) :: this

   integer, intent(in) :: src, comm

   integer :: n, ncoeffs, error

   nfe_assert(this%inited)

   ncoeffs = 1
   do n = 1, this%nextents
      error = ncoeffs
      ncoeffs = ncoeffs*this%extents(n)
      if (ncoeffs/this%extents(n).ne.error) &
         call fatal('overflow in umbrella_recv_coeffs()')
   end do

   call mpi_recv(this%coeffs, ncoeffs, MPI_DOUBLE_PRECISION, &
                 src, this%nextents, comm, MPI_STATUS_IGNORE, error)

   nfe_assert(error.eq.0)

end subroutine umbrella_recv_coeffs

#endif /* MPI */

!=============================================================================

#ifdef BINTRAJ
subroutine umbrella_load(this, filename)

   use netcdf

   implicit none

   type(umbrella_t), intent(inout) :: this
   character(*), intent(in) :: filename

   integer :: n, rc, setid, dimid, varid, ncoeffs, tmp, nextents, coeffs_len
   integer :: periodicity(UMBRELLA_MAX_NEXTENTS) ! netCDF doesn't like logical

   nfe_assert(.not.this%inited)

   rc = nf90_open(filename, nf90_nowrite, setid)
   call check_rc()

   rc = nf90_get_att(setid, nf90_global, 'nextents', nextents)
   call check_rc()

   if (nextents.lt.1.or.nextents.gt.UMBRELLA_MAX_NEXTENTS) then
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a,'//pfmt(nextents)//',a/)') &
         NFE_ERROR, 'umbrella_load(filename=''', trim(filename), &
         ''') : nextents is out of range (', nextents, ')'
      call terminate()
   end if

   this%nextents = nextents

   rc = nf90_get_att(setid, nf90_global, 'extents', this%extents(1:nextents))
   call check_rc()

   rc = nf90_get_att(setid, nf90_global, 'periodicity', periodicity(1:nextents))
   call check_rc()

   rc = nf90_get_att(setid, nf90_global, 'origin', this%origin(1:nextents))
   call check_rc()

   rc = nf90_get_att(setid, nf90_global, 'spacing', this%spacing(1:nextents))
   call check_rc()

   ncoeffs = 1
   do n = 1, nextents
      tmp = ncoeffs
      ncoeffs = ncoeffs*this%extents(n)
      if (ncoeffs/this%extents(n).ne.tmp) &
         call fatal('overflow in umbrella_load()')
      if (this%extents(n).lt.UMBRELLA_MIN_EXTENT) then
         write (unit = ERR_UNIT, fmt = '(/a,a,a,a,i1,a,i1/)') NFE_ERROR, &
            'umbrella_load(filename=''', trim(filename), &
            ''') : extents(', n, ').lt.', UMBRELLA_MIN_EXTENT
         call terminate()
      end if
      if (this%spacing(n).le.ZERO) then
         write (unit = ERR_UNIT, fmt = '(/a,a,a,a,i1,a/)') NFE_ERROR, &
            'umbrella_load(filename=''', trim(filename), &
            ''') : spacing(', n, ').le.0'
         call terminate()
      end if
      if (periodicity(n).ne.0) then
         this%periodicity(n) = .true.
      else
         this%periodicity(n) = .false.
      end if
   end do

   rc = nf90_inq_dimid(setid, 'row-major', dimid)
   call check_rc()

   rc = nf90_inquire_dimension(setid, dimid, len = coeffs_len)
   call check_rc()

   if (coeffs_len.ne.ncoeffs) then
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') NFE_ERROR, &
         'umbrella_load(filename=''', trim(filename), &
         ''') : number of ''coeffs'' is wrong'
      call terminate()
   end if

   rc = nf90_inq_varid(setid, 'coeffs', varid)
   call check_rc()

   allocate(this%coeffs(ncoeffs), stat = rc)
   if (rc.ne.0) &
      NFE_OUT_OF_MEMORY

   rc = nf90_get_var(setid, varid, this%coeffs)
   call check_rc()

   rc = nf90_close(setid)
   call check_rc()

#ifndef NFE_DISABLE_ASSERT
   this%inited = .true.
#endif /* NFE_DISABLE_ASSERT */

!-----------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------------

subroutine check_rc()

   implicit none

   character(len = 80) :: errmsg

   if (rc.ne.nf90_noerr) then
      errmsg = nf90_strerror(rc)
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a,a/)') NFE_ERROR, &
         'umbrella_load(filename=''', trim(filename), ''') : ', trim(errmsg)
      call terminate()
   end if

end subroutine check_rc

!-----------------------------------------------------------------------------

end subroutine umbrella_load

!=============================================================================

subroutine umbrella_save(this, filename)

   use netcdf

   implicit none

   type(umbrella_t), intent(in) :: this
   character(*), intent(in) :: filename

   integer :: n, rc, setid, dimid, varid, nextents, ncoeffs
   integer :: periodicity(UMBRELLA_MAX_NEXTENTS) ! netCDF doesn't like logical

   nfe_assert(this%inited)
   nfe_assert(this%nextents.ge.1)
   nfe_assert(this%nextents.le.UMBRELLA_MAX_NEXTENTS)
   nfe_assert(associated(this%coeffs))

   rc = nf90_create(filename, cmode = nf90_clobber, ncid = setid)
   call check_rc()

   nextents = this%nextents

   ncoeffs = 1
   do n = 1, nextents
      if (this%periodicity(n)) then
         periodicity(n) = 1
      else
         periodicity(n) = 0
      end if
      ncoeffs = ncoeffs*this%extents(n)
   end do

   rc = nf90_put_att(setid, nf90_global, 'nextents', nextents)
   call check_rc()

   rc = nf90_put_att(setid, nf90_global, 'extents', this%extents(1:nextents))
   call check_rc()

   rc = nf90_put_att(setid, nf90_global, 'periodicity', periodicity(1:nextents))
   call check_rc()

   rc = nf90_put_att(setid, nf90_global, 'origin', this%origin(1:nextents))
   call check_rc()

   rc = nf90_put_att(setid, nf90_global, 'spacing', this%spacing(1:nextents))
   call check_rc()

   rc = nf90_def_dim(setid, 'row-major', ncoeffs, dimid)
   call check_rc()

   rc = nf90_def_var(setid, 'coeffs', coeffs_type, dimid, varid)
   call check_rc()

   rc = nf90_enddef(setid)
   call check_rc()

   rc = nf90_put_var(setid, varid, this%coeffs)
   call check_rc()

   rc = nf90_close(setid)
   call check_rc()

!-----------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------------

subroutine check_rc()

   implicit none

   character(len = 80) :: errmsg

   if (rc.ne.nf90_noerr) then
      errmsg = nf90_strerror(rc)
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a,a/)') NFE_ERROR, &
         'umbrella_save(filename=''', trim(filename), ''') : ', trim(errmsg)
      call terminate()
   end if

end subroutine check_rc

!-----------------------------------------------------------------------------

end subroutine umbrella_save
#endif /* BINTRAJ */

!=============================================================================
subroutine umbrella_wt_mod(this, temp0, wt_temp)

   NFE_USE_AFAILED

   implicit none

   type(umbrella_t), intent(inout) :: this
   double precision, intent(in)    :: temp0, wt_temp

   integer :: n, ncoeffs
   
   nfe_assert(this%inited)

   nfe_assert(this%nextents.ge.1)
   nfe_assert(this%nextents.le.UMBRELLA_MAX_NEXTENTS)

   ncoeffs = 1
   do n = 1, this%nextents
      ncoeffs = ncoeffs*this%extents(n)
   end do

   do n = 1, ncoeffs
      this%coeffs(n) = this%coeffs(n) * (1 + temp0/wt_temp)
   end do

end subroutine umbrella_wt_mod

!=============================================================================
!
! evaluates value(f) of cubic B-spline centered
! at 0.0D0 (works correctly only for -2 < x < 2)
!

pure double precision function m4_v(x)

   implicit none

   double precision, intent(in)  :: x

   if (x.lt.MINUS_ONE) then
      m4_v = ONE_SIXTH*(x + TWO)**3
   else if (x.lt.ZERO) then
      m4_v = MINUS_HALF*x*x*(x + TWO) + TWO_THIRD
   else if (x.lt.ONE) then
      m4_v = MINUS_HALF*x*x*(TWO - x) + TWO_THIRD
   else
      m4_v = ONE_SIXTH*(TWO - x)**3
   end if

end function m4_v

!=============================================================================

!
! evaluates value(f) and derivative(df)
! of cubic B-spline centered at 0.0D0
! (works correctly only for -2 < x < 2)
!

subroutine m4_vdv(x, f, df)

   implicit none

   double precision, intent(in)  :: x
   double precision, intent(out) :: f, df

   double precision :: x2

   if (x.lt.MINUS_ONE) then
      x2 = x + TWO
      df = HALF*x2*x2
       f = ONE_THIRD*x2*df
   else if (x.lt.ZERO) then
      x2 = MINUS_HALF*x
       f = x2*x*(x + TWO) + TWO_THIRD
      df = x2*(THREE*x + FOUR)
   else if (x.lt.ONE) then
      x2 = HALF*x
       f = x2*x*(x - TWO) + TWO_THIRD
      df = x2*(THREE*x - FOUR)
   else
      x2 = TWO - x
      df = MINUS_HALF*x2*x2
       f = MINUS_ONE_THIRD*x2*df
   end if

end subroutine m4_vdv

!=============================================================================

! Added by M Moradi
!
! evaluates value(f) and second derivative(d2f)
! of cubic B-spline centered at 0.0D0
! (works correctly only for -2 < x < 2)
!

subroutine m4_vd2v(x, f, d2f)

   implicit none

   double precision, intent(in)  :: x
   double precision, intent(out) :: f, d2f

   double precision :: x2, df

   if (x.lt.MINUS_ONE) then
      x2 = x + TWO
      df = HALF*x2*x2
       f = ONE_THIRD*x2*df
     d2f = x2
   else if (x.lt.ZERO) then
      x2 = MINUS_HALF*x
       f = x2*x*(x + TWO) + TWO_THIRD
      df = x2*(THREE*x + FOUR)
     d2f = -(THREE*x + TWO)
   else if (x.lt.ONE) then
      x2 = HALF*x
       f = x2*x*(x - TWO) + TWO_THIRD
      df = x2*(THREE*x - FOUR)
     d2f = THREE*x - TWO
   else
      x2 = TWO - x
      df = MINUS_HALF*x2*x2
       f = MINUS_ONE_THIRD*x2*df
     d2f = x2
   end if

end subroutine m4_vd2v

! Moradi end
!=============================================================================
end module nfe_lib_mod
