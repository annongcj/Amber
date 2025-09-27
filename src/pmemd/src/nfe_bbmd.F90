!---------------------!
! Bizarrely Biased MD !
!---------------------!

!------------------------------------------------------------------
! Updated in Sep. 2018, since the compatible H-REMD is implemented by
! F. Pan (with -rem 3), this module seems redundant, but we still
! keep it for records.
!------------------------------------------------------------------

! Pre-processing
#ifndef NFE_UTILS_H
#define NFE_UTILS_H

#ifndef NFE_DISABLE_ASSERT
#  define nfe_assert(stmt) if (.not.(stmt)) call afailed(__FILE__, __LINE__)
#  define nfe_assert_not_reached() call afailed(__FILE__, __LINE__)
#  define NFE_PURE_EXCEPT_ASSERT
#  define NFE_USE_AFAILED use nfe_lib_mod, only : afailed
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

#define ASSUME_GFORTRAN yes
!-------------------------------------------------------------------------

module nfe_bbmd_mod

#ifdef MPI

use nfe_lib_mod, only : SL => STRING_LENGTH, LOG_UNIT => BBMD_LOG_UNIT, BBMD_MONITOR_UNIT, &
                        umbrella_t, MAX_NUMBER_OF_COLVARS => UMBRELLA_MAX_NEXTENTS, &
                        BBMD_CV_UNIT
use nfe_colvar_mod, only : colvar_t

implicit none

public :: on_pmemd_init
public :: on_pmemd_exit

public :: on_force
public :: on_mdstep
public :: on_mdwrit

private :: ctxt_init
private :: ctxt_fini

private :: ctxt_print
private :: ctxt_bcast

private :: ctxt_on_force
private :: ctxt_on_mdwrit

private :: ctxt_Um
private :: ctxt_Uo

private :: ctxt_send
private :: ctxt_recv

private :: ctxt_close_units
private :: ctxt_open_units
!<>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><>
! mt19937 seed part
!<>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><>
intrinsic :: bit_size

private

integer,  parameter  :: intg = selected_int_kind( 9 )
integer,  parameter  :: long = selected_int_kind( 18 )
integer,  parameter  :: flot = selected_real_kind( 6, 37 )
integer,  parameter  :: dobl = selected_real_kind( 15, 307 )

integer,  parameter :: wi = intg
integer,  parameter :: wl = long
integer,  parameter :: wr = dobl

! Period parameters
integer( kind = wi ), parameter :: n = 624_wi
integer( kind = wi ), parameter :: m = 397_wi
integer( kind = wi ), parameter :: hbs = bit_size( n ) / 2_wi
integer( kind = wi ), parameter :: qbs = hbs / 2_wi
integer( kind = wi ), parameter :: tbs = 3_wi * qbs

type :: mt19937_t
  private
  integer(kind = wi) :: mt(n) ! the array for the state vector
  logical(kind = wi) :: mtinit = .false._wi ! means mt[N] is not initialized
  integer(kind = wi) :: mti = n + 1_wi ! mti==N+1 means mt[N] is not initialized
end type mt19937_t

private :: init_by_seed
private :: init_by_array

private :: random_int32
private :: random_int31

private :: random_real1
private :: random_real2
private :: random_real3

private :: random_res53

private :: check_ncrc
private :: mt19937_save
private :: mt19937_load
private :: mt19937_bcast
!<>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><>

integer, private, parameter :: MONITOR_UNIT = BBMD_MONITOR_UNIT
integer, private, parameter :: CV_UNIT = BBMD_CV_UNIT

character(*), private, parameter :: &
   DEFAULT_MONITOR_FILE = 'nfe-bbmd-monitor', &
   DEFAULT_UMBRELLA_FILE = 'nfe-bbmd-umbrella', &
   DEFAULT_WT_UMBRELLA_FILE = 'nfe-bbmd-wt-umbrella', &
   DEFAULT_SNAPSHOTS_BASENAME = 'nfe-bbmd-umbrella-snapshot'

integer, private, parameter :: MODE_NONE = 5432
integer, private, parameter :: MODE_ANALYSIS = 1234
integer, private, parameter :: MODE_UMBRELLA = 2345
integer, private, parameter :: MODE_FLOODING = 3456

integer, public, save :: active = 0 ! integer for MPI

type, private :: bbmd_ctxt_t

   character(SL) :: mdout ! master only
   character(SL) :: restrt ! master only
   character(SL) :: mdvel ! master only
   character(SL) :: mden ! master only
   character(SL) :: mdcrd ! master only
   character(SL) :: mdinfo ! master only

   integer :: ioutfm
   integer :: ntpr
   integer :: ntwr
   integer :: ntwx

   character(SL) :: monitor_file ! master only
   character(SL) :: umbrella_file ! master only
   character(SL) :: wt_umbrella_file ! master only
   character(SL) :: snapshots_basename ! master only

   character(SL) :: monitor_fmt ! master only

   integer :: monitor_freq
   integer :: snapshots_freq ! master only

   double precision :: timescale ! master only

   integer :: imode
   integer :: ncolvars

   type(colvar_t)   :: colvars(MAX_NUMBER_OF_COLVARS)
   type(umbrella_t) :: umbrella ! master only (mytaskid.eq.0)
   type(umbrella_t) :: wt_umbrella ! master only (mytaskid.eq.0)

   logical :: umbrella_file_existed ! home-master only
   logical :: umbrella_discretization_changed  ! home-master only

! Added by M Moradi
! Well-tempered ABMD
   double precision :: pseudo
! Driven ABMD
   integer   :: drivenw
   integer   :: drivenu
   double precision :: driven_cutoff
! Moradi end

#ifndef NFE_DISABLE_ASERT
   logical :: initialized = .false.
#endif /* NFE_DISABLE_ASSERT */

end type bbmd_ctxt_t

type(bbmd_ctxt_t), private, allocatable, save :: contexts(:)

integer, public,  save :: mdstep = 0
integer, private, save :: exchno

integer, private, save :: current
integer, private, save :: partner
integer, private, save :: partner_masterrank

double precision, private, save :: U_mm, U_mo, U_om, U_oo

character(SL),    private, save :: mode = 'NONE'
character(SL),    private, save :: monitor_file ! master only
integer,          private, save :: monitor_freq = 50
double precision, private, save :: timescale = 1 ! master only
character(SL),    private, save :: umbrella_file ! master only
character(SL),    private, save :: snapshots_basename ! master only
integer,          private, save :: snapshots_freq = -1 ! master only
double precision, private, save :: wt_temperature = 0.0
character(SL),    private, save :: wt_umbrella_file ! master only
character(SL),    private, save :: cv_file = 'nfe-bbmd-cv'
character(SL),    private, save :: driven_weight = 'NONE'
double precision, private, save :: driven_cutoff = 0.0


integer, public, save :: exchange_freq = 1
character(SL), private, save :: exchange_log_file
integer, private, save :: exchange_log_freq

integer, private, save :: mt19937_seed
character(SL), private, save :: mt19937_file

type(mt19937_t), private, save :: mersenne_twister
integer, private, allocatable, save :: permutation(:)

! symmetric matrices with 0s on the diagonal
! [row major : (2,1), (3,1), (3,2), (4,1), ... ]
integer, private, allocatable, save :: local_exchg_attempts(:)
integer, private, allocatable, save :: local_exchg_accepted(:)

! on masterrank.eq.0 (layout as above)
integer, private, allocatable, save :: exchg_attempts(:)
integer, private, allocatable, save :: exchg_accepted(:)

namelist / bbmd /    mode, monitor_file, monitor_freq, timescale,&
                     umbrella_file, snapshots_basename, snapshots_freq,&
                     wt_temperature, wt_umbrella_file, cv_file, &
                     driven_weight, driven_cutoff, &
                     exchange_freq, exchange_log_file, exchange_log_freq, &
                     mt19937_seed, mt19937_file


contains

!<>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><>

subroutine on_pmemd_init(mdin_unit, amass)

   use nfe_lib_mod
   use parallel_dat_mod
   use file_io_mod
   use mdin_ctrl_dat_mod, only : ioutfm, ntwv, ntwf

   implicit none

   integer, intent(in) :: mdin_unit

   double precision, intent(in) :: amass(*)

   integer :: n, error, mt19937_loaded, ifind
   integer, allocatable :: active_all(:)

#  ifndef NFE_DISABLE_ASSERT
   logical, save :: called = .false.
#  endif /* NFE_DISABLE_ASSERT */

   logical :: found

   nfe_assert(active.eq.0)
   nfe_assert(multipmemd_rem().eq.0)
   nfe_assert(multipmemd_numgroup().gt.1)

   nfe_assert(.not.called)
   nfe_assert(.not.allocated(contexts))

   if (ioutfm.ne.0) then
       if ((ntwv.ne.0).or.(ntwf.ne.0)) then
         write (unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR, &
              'please set ntwv and ntwf to zero for binary NETCDF output'
         call terminate()
       end if
   end if

#  ifndef NFE_DISABLE_ASSERT
   called = .true.
#  endif /* NFE_DISABLE_ASSERT */

   NFE_MASTER_ONLY_BEGIN
      rewind(mdin_unit)
      call nmlsrc('bbmd', mdin_unit, ifind)
      if (ifind.ne.0) then
         active = 1
      else
         active = 0
      end if

      nfe_assert(pmemd_master_comm.ne.MPI_COMM_NULL)

      allocate(active_all(master_size), stat = error)
      if (error.ne.0) &
         NFE_OUT_OF_MEMORY

      call mpi_gather(active, 1, MPI_INTEGER, active_all, &
         1, MPI_INTEGER, 0, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      if (master_rank.eq.0) then
         do n = 1, master_size
            if (active.ne.active_all(n)) &
               call fatal('either all or none of MDIN files &
                          &should contain the &bbmd namelist')
         end do
      end if ! master_rank.eq.0

      deallocate(active_all)
   NFE_MASTER_ONLY_END

   call mpi_bcast(active, 1, MPI_INTEGER, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)

   if (active.eq.0) &
      return

   ! <<>> BBMD is requested in all replicas <<>>

   call mpi_bcast(master_size, 1, MPI_INTEGER, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)

   call mpi_bcast(master_rank, 1, MPI_INTEGER, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)

   allocate(contexts(master_size), &
      permutation(master_size + mod(master_size, 2)), stat = error)
   if (error.ne.0) &
      NFE_OUT_OF_MEMORY


   NFE_MASTER_ONLY_BEGIN
   ! setup defaults

      write (unit = monitor_file, fmt = '(a,a,i3.3)') &
        DEFAULT_MONITOR_FILE, '-', (master_rank + 1)
      write (unit = umbrella_file, fmt = '(a,a,i3.3,a)') &
        DEFAULT_UMBRELLA_FILE, '-', (master_rank + 1), '.nc'
      write (unit = snapshots_basename, fmt = '(a,a,i3.3)') &
        DEFAULT_SNAPSHOTS_BASENAME, '-', (master_rank + 1)
      write (unit = wt_umbrella_file, fmt = '(a,a,i3.3,a)') &
        DEFAULT_WT_UMBRELLA_FILE, '-', (master_rank + 1), '.nc'

      mt19937_file = 'nfe-mt19937.nc'
      mt19937_seed = 5489
      exchange_log_file = 'nfe-bbmd-exchange-log'
      exchange_log_freq = 100

      rewind(mdin_unit)
      read(mdin_unit,nml=bbmd,err=666)
   NFE_MASTER_ONLY_END

   call ctxt_init(contexts(master_rank + 1), amass)

   do n = 1, master_size
      call ctxt_bcast(contexts(n), n - 1, amass)
   end do

!   call mpi_barrier(commworld, error) ! unneeded actually
!   nfe_assert(error.eq.0)

   ! exchange frequency, logfile, etc

   NFE_MASTER_ONLY_BEGIN

      ! try load the twister first
      mt19937_loaded = 0

#     ifdef BINTRAJ
      if (master_rank.eq.0) then
         inquire (file = mt19937_file, exist = found)
         if (found) then
            call mt19937_load(mersenne_twister, mt19937_file)
            mt19937_loaded = 1
         end if ! found
      end if ! master_rank.eq.0
#     else
      write (unit = ERR_UNIT, fmt = '(a,a)') NFE_WARNING, &
         'mt19937 : netCDF is not available (try ''-bintraj'' configure option)'
#     endif /* BINTRAJ */

      if (mt19937_loaded.eq.0) then
         call init_by_seed(mersenne_twister, mt19937_seed)
      end if ! .not.mt19937_loaded

      call mt19937_bcast(mersenne_twister, pmemd_master_comm, 0)

      call mpi_bcast(mt19937_file, len(mt19937_file), &
                     MPI_CHARACTER, 0, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(mt19937_seed, 1, MPI_INTEGER, 0, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(mt19937_loaded, 1, MPI_INTEGER, 0, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(exchange_freq, 1, MPI_INTEGER, 0, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(exchange_log_file, len(exchange_log_file), &
                     MPI_CHARACTER, 0, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(exchange_log_freq, 1, MPI_INTEGER, 0, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      allocate(local_exchg_attempts(master_size*(master_size - 1)/2), &
         local_exchg_accepted(master_size*(master_size - 1)/2), stat = error)
      if (error.ne.0) &
         NFE_OUT_OF_MEMORY

      local_exchg_attempts(:) = 0
      local_exchg_accepted(:) = 0

      if (master_rank.eq.0) then
         allocate(exchg_attempts(master_size*(master_size - 1)/2), &
            exchg_accepted(master_size*(master_size - 1)/2), stat = error)
         if (error.ne.0) &
            NFE_OUT_OF_MEMORY

         open (unit = LOG_UNIT, file = exchange_log_file, iostat = error, &
               form = 'FORMATTED', action = 'WRITE', status = 'REPLACE')
         if (error.ne.0) then
            write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
               NFE_ERROR, 'failed to open ''', &
               trim(exchange_log_file), ''' for writing'
            call terminate()
         end if
      end if ! master_rank.eq.0
   NFE_MASTER_ONLY_END

   mdstep = 0
   exchno = 0

   current = master_rank + 1
   partner = -1

   call mpi_bcast(exchange_freq, 1, MPI_INTEGER, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)

   if (mytaskid.gt.0) &
      return

   write (unit = OUT_UNIT, fmt = '(/a,a)') NFE_INFO, &
      '/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/ B. B. M. D. /=/=/=/=/=/=/=/=/=/=/=/=/=/'
   write (unit = OUT_UNIT, &
      fmt = '(a,/a,a,'//pfmt(current)//',a,'//pfmt(master_size)//')') &
      NFE_INFO, NFE_INFO, &
      'this is replica ', current, ' of ', master_size
   write (unit = OUT_UNIT, fmt = '(a,a,'//pfmt(exchange_freq)//')') &
      NFE_INFO, 'number of MD steps between exchange attempts = ', exchange_freq
#  ifndef BINTRAJ
   write (unit = OUT_UNIT, fmt = '(a,a,'//pfmt(mt19937_seed)//',a)') &
      NFE_INFO, 'mersenne twister seed = ', mt19937_seed, ' (no netCDF)'
#  else
   write (unit = OUT_UNIT, fmt = '(a,a,a,a)') NFE_INFO, &
      'mersenne twister file = ''', trim(mt19937_file), ''''
   if (mt19937_loaded.eq.1) then
      write (unit = OUT_UNIT, fmt = '(a,a)') NFE_INFO, &
         'mersenne twister state loaded'
   else
      write (unit = OUT_UNIT, fmt = '(a,a,'//pfmt(mt19937_seed)//')') &
         NFE_INFO, 'mersenne twister seed = ', mt19937_seed
   end if
#  endif /* BINTRAJ */
   write (unit = OUT_UNIT, fmt = '(a,a,a,a)') NFE_INFO, &
      'exchange log file = ''', trim(exchange_log_file), ''''
   write (unit = OUT_UNIT, fmt = '(a,a,'//pfmt(exchange_log_freq)//',/a)') &
      NFE_INFO, 'exchange log update frequency = ', &
      exchange_log_freq, NFE_INFO

   call ctxt_print(contexts(master_rank + 1), OUT_UNIT)

   write (unit = OUT_UNIT, fmt = '(a,/a,a/)') NFE_INFO, NFE_INFO, &
      '/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/'

   return
666 write(unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR,'Cannot read &bbmd namelist!'
    call terminate()
end subroutine on_pmemd_init

!<>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><>

subroutine on_pmemd_exit()

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   integer :: n

   if (active.eq.0) then
      nfe_assert(.not.allocated(contexts))
      nfe_assert(.not.allocated(permutation))
      return
   end if ! active.eq.0

   nfe_assert(master_size.gt.1)

   do n = 1, master_size
      call ctxt_fini(contexts(n))
   end do

   deallocate(contexts, permutation)

   NFE_MASTER_ONLY_BEGIN
      deallocate(local_exchg_attempts, local_exchg_accepted)

      if (master_rank.eq.0) then
         deallocate(exchg_attempts, exchg_accepted)
         close (LOG_UNIT)
      end if ! master_rank.eq.0
   NFE_MASTER_ONLY_END

   active = 0

end subroutine on_pmemd_exit

!<>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><>

! Modified by M Moradi
! Driven ABMD
subroutine on_force(x, f, work, udr, pot)
! Moradi end

   NFE_USE_AFAILED

   use parallel_dat_mod

   implicit none

   double precision, intent(in) :: x(3,*)

   double precision, intent(inout) :: f(3,*)

! Modified by M Moradi
! Driven ABMD
   double precision, intent(in) :: work
   double precision, intent(in) :: udr
! Moradi end
   double precision, intent(inout) :: pot

   integer :: n, k, error

   if (active.eq.0) &
      return

   nfe_assert(master_size.gt.1)

   nfe_assert(current.gt.0)
   nfe_assert(current.le.master_size)
   nfe_assert(allocated(contexts))
   nfe_assert(allocated(permutation))

   call ctxt_on_force(contexts(current), x, f, mdstep, work, udr, pot)

   nfe_assert(partner.lt.0)

   ! prepare for an exchange attempt if needed
   if (mod(mdstep + 1, exchange_freq).ne.0) &
      return

   ! generate a random permutation
   NFE_MASTER_ONLY_BEGIN
      do n = 1, master_size + mod(master_size, 2)
         permutation(n) = n
      end do
      do n = 1, master_size + mod(master_size, 2) - 1
         k = 1 + n + mod(random_int31(mersenne_twister), &
            master_size + mod(master_size, 2) - n)
         error = permutation(k)
         permutation(k) = permutation(n)
         permutation(n) = error
      end do
      partner_masterrank = -1
      do n = 1, master_size + mod(master_size + 1, 2), 2
         if (permutation(n).le.master_size.and. &
            permutation(n + 1).le.master_size) then
            if (permutation(n).eq.(master_rank + 1)) then
               partner_masterrank = permutation(n + 1) - 1
               exit
            end if
            if (permutation(n + 1).eq.(master_rank + 1)) then
               partner_masterrank  = permutation(n) - 1
               exit
            end if
         end if
      end do
   NFE_MASTER_ONLY_END

   call mpi_bcast(partner_masterrank, 1, MPI_INTEGER, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)

   nfe_assert(master_rank.ne.partner_masterrank)

   if (partner_masterrank.lt.0) & ! for mod(master_size, 2).ne.0
      return

   ! get the partner's context number
   NFE_MASTER_ONLY_BEGIN
      call mpi_sendrecv(current, 1, MPI_INTEGER, partner_masterrank, mdstep, &
                        partner, 1, MPI_INTEGER, partner_masterrank, mdstep, &
                        pmemd_master_comm, MPI_STATUS_IGNORE, error)
      nfe_assert(error.eq.0)
      nfe_assert(partner.ge.1.and.partner.le.master_size)
   NFE_MASTER_ONLY_END

   call mpi_bcast(partner, 1, MPI_INTEGER, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)

   ! U_?? are needed only in the replica with smaller masterrank
   if (master_rank.lt.partner_masterrank) then
      call ctxt_Um(contexts(current), partner_masterrank, x, U_mm, U_mo)
      call ctxt_Uo(contexts(partner), partner_masterrank, x, U_om, U_oo)
   else
      call ctxt_Um(contexts(partner), partner_masterrank, x, U_mm, U_mo)
      call ctxt_Uo(contexts(current), partner_masterrank, x, U_om, U_oo)
   end if ! master_rank.lt.partner_masterrank

end subroutine on_force

!<>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><>

subroutine on_mdstep(E_p, v)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   double precision, intent(in) :: E_p
   double precision, intent(inout) :: v(*)

   double precision :: v_scale
   double precision :: E_px, partner_E_px, temp0, partner_temp0
   double precision :: beta_m, beta_o, delta, random_number

   integer :: error, exchange, n

   if (active.eq.0) &
      return

   nfe_assert(master_size.gt.1)

   nfe_assert(current.gt.0)
   nfe_assert(current.le.master_size)
   nfe_assert(allocated(contexts))

   mdstep = mdstep + 1

   if (mod(mdstep, exchange_freq).ne.0) then
      if (mdstep.eq.pmemd_nstlim()) &
         goto 2 ! update_log()
      return
   end if ! mod(mdstep, exchange_freq).ne.0

   random_number = random_res53(mersenne_twister)

   exchno = exchno + 1

   if (partner_masterrank.lt.0) & ! this replica does not participate
      goto 2 ! update_log()

   nfe_assert(master_rank.ne.partner_masterrank)
   nfe_assert(partner.ge.1.and.partner.le.master_size)
   nfe_assert(partner.ne.current)

   temp0 = pmemd_temp0()

   if (mytaskid.ne.0) &
      goto 1 ! bcast(exchange)

   ! added by F Pan, since we already include the potential contribution
   ! from bbmd in E_p, we need to deduct it
   E_px = E_p - nfe_pot_ene%bbmd

   ! get partner's E_p
   if (master_rank.lt.partner_masterrank) then
      call mpi_recv(partner_E_px, 1, MPI_DOUBLE_PRECISION, &
                    partner_masterrank, mdstep, pmemd_master_comm, &
                    MPI_STATUS_IGNORE, error)
      nfe_assert(error.eq.0)
   else
      call mpi_send(E_px, 1, MPI_DOUBLE_PRECISION, &
                    partner_masterrank, mdstep, pmemd_master_comm, error)
      nfe_assert(error.eq.0)
   end if

   call mpi_sendrecv &
      (temp0, 1, MPI_DOUBLE_PRECISION, partner_masterrank, mdstep, &
       partner_temp0, 1, MPI_DOUBLE_PRECISION, partner_masterrank, mdstep, &
       pmemd_master_comm, MPI_STATUS_IGNORE, error)
   nfe_assert(error.eq.0)

   if (master_rank.lt.partner_masterrank) then

      nfe_assert(temp0.gt.ZERO)
      nfe_assert(partner_temp0.gt.ZERO)

      beta_m = 503.01D0/temp0
      beta_o = 503.01D0/partner_temp0

      delta = (beta_o - beta_m)*(E_px - partner_E_px) &
         + beta_m*(U_mo - U_mm) - beta_o*(U_oo - U_om)

      if (delta.lt.ZERO) then
         exchange = 1
      else if (exp(-delta).gt.random_number) then
         exchange = 1
      else
         exchange = 0
      end if

      call mpi_send(exchange, 1, MPI_INTEGER, partner_masterrank, &
                    mdstep, pmemd_master_comm, error)
      nfe_assert(error.eq.0)
   else
      call mpi_recv(exchange, 1, MPI_INTEGER, partner_masterrank, &
                    mdstep, pmemd_master_comm, MPI_STATUS_IGNORE, error)
      nfe_assert(error.eq.0)
   end if ! master_rank.lt.partner_masterrank

   if (current.lt.partner) then
      error = (partner - 1)*(partner - 2)/2 + current
      local_exchg_attempts(error) = local_exchg_attempts(error) + 1
      if (exchange.eq.1) &
         local_exchg_accepted(error) = local_exchg_accepted(error) + 1
   end if ! current.lt.partner

1  call mpi_bcast(exchange, 1, MPI_INTEGER, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)

   if (exchange.eq.1) then

      NFE_MASTER_ONLY_BEGIN
         write (unit = OUT_UNIT, &
         fmt ='(a,a,'//pfmt(partner)//',a,'//pfmt(pmemd_mdtime(), 3)//')') &
            NFE_INFO, 'BBMD : exchanged coordinates/velocities with ', &
            partner, ' at t = ', pmemd_mdtime()
         call ctxt_close_units(contexts(current))
         if (master_rank.lt.partner_masterrank) then
            call ctxt_recv(contexts(partner), partner_masterrank)
            call ctxt_send(contexts(current), partner_masterrank)
         else
            call ctxt_send(contexts(current), partner_masterrank)
            call ctxt_recv(contexts(partner), partner_masterrank)
         end if ! master_rank.lt.partner_masterrank
         call ctxt_open_units(contexts(partner))
      NFE_MASTER_ONLY_END

      current = partner

      call mpi_bcast(partner_temp0, 1, MPI_DOUBLE_PRECISION, &
                     0, pmemd_comm, error)
      nfe_assert(error.eq.0)

      if (abs(temp0 - partner_temp0).gt.TINY) then
         v_scale = sqrt(partner_temp0/temp0)
         do n = 1, 3*pmemd_natoms()
            v(n) = v_scale*v(n)
         end do
!         ekmh = ekmh*v_scale*v_scale
      end if

      call set_pmemd_temp0(partner_temp0)
   end if ! exchange.eq.1

2  NFE_MASTER_ONLY_BEGIN
      if (mod(exchno, exchange_log_freq).eq.0.or.mdstep.eq.pmemd_nstlim()) &
         call update_log()
   NFE_MASTER_ONLY_END

#  ifndef NFE_DISABLE_ASSERT
   partner = -1
#  endif /* NFE_DISABLE_ASSERT */

contains

subroutine update_log()

   implicit none

   integer :: n1, n2, idx
   character(32) :: log_fmt

   nfe_assert(mytaskid.eq.0)

   if (master_rank.eq.0) then
      call mpi_reduce(local_exchg_attempts, exchg_attempts, &
         master_size*(master_size - 1)/2, MPI_INTEGER, MPI_SUM, &
         0, pmemd_master_comm, error)
   else
      call mpi_reduce(local_exchg_attempts, 0, &
         master_size*(master_size - 1)/2, MPI_INTEGER, MPI_SUM, &
         0, pmemd_master_comm, error)
   end if
   nfe_assert(error.eq.0)

   if (master_rank.eq.0) then
      call mpi_reduce(local_exchg_accepted, exchg_accepted, &
         master_size*(master_size - 1)/2, MPI_INTEGER, MPI_SUM, &
         0, pmemd_master_comm, error)
   else
      call mpi_reduce(local_exchg_accepted, 0, &
         master_size*(master_size - 1)/2, MPI_INTEGER, MPI_SUM, &
         0, pmemd_master_comm, error)
   end if
   nfe_assert(error.eq.0)

   if (master_rank.ne.0) &
      return

   nfe_assert(allocated(exchg_attempts))
   nfe_assert(allocated(exchg_accepted))

   if (exchno.gt.exchange_log_freq) &
      write (unit = LOG_UNIT, fmt = '(80(''=''))')

   write (unit = LOG_UNIT, &
      fmt = '(/a,'//pfmt((master_size - mod(master_size, 2)/2)/2)//',a,'//pfmt &
      (exchno)//',a/)') ' <> exchange statistics over ', &
      ((master_size - mod(master_size, 2)/2)/2),' x ', exchno, ' attempts <>'

   idx = 2 + int(floor(log10(ONE + dble(maxval(exchg_attempts)))))
   write (unit = log_fmt, fmt = '(a,'//pfmt(idx)//',a)') '(i', idx, ')'

   write (unit = LOG_UNIT, fmt = '(a/)') ' * number of trials:'
   write (unit = LOG_UNIT, fmt = '(6x)', advance = 'NO')

   do n1 = 1, master_size - 1
      write (unit = LOG_UNIT, fmt = log_fmt, advance = 'NO') n1
   end do

   do n1 = 2, master_size
      write (unit = LOG_UNIT, fmt = '(/i3,1x,'':'',1x)', advance = 'NO') n1
      do n2 = 1, n1 - 1
         idx = (n1 - 1)*(n1 - 2)/2 + n2
         write (unit = LOG_UNIT, fmt = log_fmt, advance = 'NO') &
            exchg_attempts(idx)
      end do
   end do

   write (unit = LOG_UNIT, fmt = '(//a/)') ' * number of exchanges:'
   write (unit = LOG_UNIT, fmt = '(6x)', advance = 'NO')

   do n1 = 1, master_size - 1
      write (unit = LOG_UNIT, fmt = log_fmt, advance = 'NO') n1
   end do

   do n1 = 2, master_size
      write (unit = LOG_UNIT, fmt = '(/i3,1x,'':'',1x)', advance = 'NO') n1
      do n2 = 1, n1 - 1
         idx = (n1 - 1)*(n1 - 2)/2 + n2
         write (unit = LOG_UNIT, fmt = log_fmt, advance = 'NO') &
            exchg_accepted(idx)
      end do
   end do

   write (unit = LOG_UNIT, fmt = '(//a/)') ' * acceptance rates:'
   write (unit = LOG_UNIT, fmt = '(6x)', advance = 'NO')

   do n1 = 1, master_size - 1
      write (unit = LOG_UNIT, fmt = '(i6)', advance = 'NO') n1
   end do

   do n1 = 2, master_size
      write (unit = LOG_UNIT, fmt = '(/i3,1x,'':'',1x)', advance = 'NO') n1
      do n2 = 1, n1 - 1
         idx = (n1 - 1)*(n1 - 2)/2 + n2
         if (exchg_attempts(idx).gt.0) then
            write (unit = LOG_UNIT, fmt = '(f6.1)', advance = 'NO') &
              ((dble(100)*exchg_accepted(idx))/exchg_attempts(idx))
         else
            write (unit = LOG_UNIT, fmt = '(a)', advance = 'NO') '  XX.X'
         end if ! exchg_attempts(idx).gt.0
      end do
   end do

   write (unit = LOG_UNIT, fmt = '(/a)') ''
   call flush_UNIT(LOG_UNIT, exchange_log_file)

end subroutine update_log

end subroutine on_mdstep

!<>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><>

subroutine on_mdwrit()

   NFE_USE_AFAILED
   use parallel_dat_mod

   implicit none

   if (active.eq.0) &
      return

   nfe_assert(master_size.gt.1)
   nfe_assert(current.gt.0)
   nfe_assert(current.le.master_size)
   nfe_assert(allocated(contexts))

   call ctxt_on_mdwrit(contexts(current))

#  ifdef BINTRAJ
   if (mytaskid.eq.0.and.master_rank.eq.0) &
      call mt19937_save(mersenne_twister, mt19937_file)
#  endif /* BINTRAJ */

end subroutine on_mdwrit
!-------------------------------------------------------------------------------

subroutine ctxt_init(self, amass)

   use nfe_lib_mod
   use nfe_colvar_mod
   use file_io_dat_mod
   use file_io_mod
   use parallel_dat_mod
   use mdin_ctrl_dat_mod

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self
   double precision, intent(in) :: amass(*)

   logical :: umbrella_file_exists
   type(umbrella_t) :: umbrella_from_file

   logical :: do_transfer
   integer :: n, error, ifind, i

   integer :: cv_extents(UMBRELLA_MAX_NEXTENTS)
   logical :: cv_periodicity(UMBRELLA_MAX_NEXTENTS)

   double precision :: cv_origin(UMBRELLA_MAX_NEXTENTS)
   double precision :: cv_spacing(UMBRELLA_MAX_NEXTENTS)

   double precision :: tmp
   character(len = 80) :: buf

   nfe_assert(.not.self%initialized)

   if (mytaskid.gt.0) &
      goto 1

   ! store PMEMD's filenames

   self%mdout  = mdout_name
   self%restrt = restrt_name
   self%mdvel  = mdvel_name
   self%mden   = mden_name
   self%mdcrd  = mdcrd_name
   self%mdinfo = mdinfo_name

   self%ioutfm = ioutfm
   self%ntpr = ntpr
   self%ntwr = ntwr
   self%ntwx = ntwx

   self%ncolvars = 0

   ! discover the run-mode

   if (mode == 'NONE') then
      self%imode = MODE_NONE
      goto 1
   else if (mode == 'ANALYSIS') then
      self%imode = MODE_ANALYSIS
   else if (mode == 'UMBRELLA') then
      self%imode = MODE_UMBRELLA
   else if (mode == 'FLOODING') then
      self%imode = MODE_FLOODING
   else
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
         NFE_ERROR, 'unknown mode ''', trim(mode), ''''
      call terminate()
   end if

! Added by M Moradi
! Well-tempered ABMD
   self%pseudo = ZERO
   if (wt_temperature .gt. ZERO) then
      self%pseudo = ONE/wt_temperature
   end if
   self%wt_umbrella_file = wt_umbrella_file
!  Driven ABMD
   if (driven_weight == 'NONE') then
      self%drivenw = 0
      self%drivenu = 0
   else if (driven_weight == 'CONSTANT') then
      self%drivenw = 1
      self%drivenu = 1
   else if (driven_weight == 'PULLING') then
      self%drivenw = 1
      self%drivenu = 0
   else
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
         NFE_ERROR, 'unknown driven weight scheme ''', trim(driven_weight), ''''
      call terminate()
   end if

   self%driven_cutoff = driven_cutoff
! Moradi end

   self%umbrella_file = umbrella_file

#ifndef BINTRAJ
   umbrella_file_exists = .false.
   write (unit = ERR_UNIT, fmt = '(a,a)') NFE_WARNING, &
      'netCDF is not available (try ''-bintraj'' configure option)'
#else
   inquire (file = self%umbrella_file, exist = umbrella_file_exists)
#endif /* BINTRAJ */

   if (.not.umbrella_file_exists.and.self%imode.eq.MODE_UMBRELLA) then
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') NFE_ERROR, '''', &
      trim(self%umbrella_file), ''' does not exist (required for UMBRELLA mode)'
      call terminate()
   end if

   if (self%imode.eq.MODE_ANALYSIS) &
      umbrella_file_exists = .false.

#ifdef BINTRAJ
   if (umbrella_file_exists) &
      call umbrella_load(umbrella_from_file, self%umbrella_file)
#endif /* BINTRAJ */

   ! collective variables
   nfe_assert(self%ncolvars.eq.0)

   call amopen(CV_UNIT, cv_file, 'O', 'F', 'R')

   do
     call nmlsrc('colvar', CV_UNIT, ifind)
     if (ifind.eq.0) exit
     read(CV_UNIT,'(a80)') buf
     self%ncolvars = self%ncolvars + 1
   end do

   if (self%ncolvars.eq.0) &
      call fatal('no variable(s) in the CV file')

   if (self%ncolvars.gt.MAX_NUMBER_OF_COLVARS) &
      call fatal('too many variables in the CV file')

   if (umbrella_file_exists) then
      if(umbrella_nextents(umbrella_from_file).ne.self%ncolvars) &
         call fatal('number of variables in the CV file does not &
                 &match with the number of extents found in the umbrella_file')
   end if ! umbrella_file_exists

   n = 1

   do while (n.le.self%ncolvars)

         call colvar_nlread(CV_UNIT, self%colvars(n))

         if (self%imode.eq.MODE_FLOODING) then
            cv_spacing(n) = resolution
            cv_spacing(n) = cv_spacing(n)/4
            if ((resolution.eq.ZERO).and..not.umbrella_file_exists) then
               write (unit = ERR_UNIT, fmt = '(/a,a,i1/)') NFE_ERROR, &
                  'could not determine ''resolution'' for CV #', n
               call terminate()
            end if

            if (resolution.eq.ZERO) &
               cv_spacing(n) = umbrella_spacing(umbrella_from_file, n)

            cv_periodicity(n) = colvar_is_periodic(self%colvars(n))

            if (cv_periodicity(n)) then

               nfe_assert(colvar_has_min(self%colvars(n)))
               nfe_assert(colvar_has_max(self%colvars(n)))

               cv_origin(n) = colvar_min(self%colvars(n))

               nfe_assert(cv_spacing(n).gt.ZERO)
               nfe_assert(colvar_max(self%colvars(n)).gt.cv_origin(n))

               cv_extents(n) = &
               int((colvar_max(self%colvars(n)) - cv_origin(n))/cv_spacing(n))

               if (cv_extents(n).lt.UMBRELLA_MIN_EXTENT) then
                  write (unit = ERR_UNIT, fmt = '(/a,a,i1,a/)') NFE_ERROR, &
                     'CV #', n, ' : ''resolution'' is too big'
                  call terminate()
               end if

               cv_spacing(n) = &
                  (colvar_max(self%colvars(n)) - cv_origin(n))/cv_extents(n)

            else ! .not.periodic

               cv_origin(n) = cv_min
               tmp = cv_max

               if (cv_origin(n).ge.tmp) then
                  write (unit = ERR_UNIT, fmt = '(/a,a,i1/)') NFE_ERROR, &
                     'min.ge.max for CV #', n
                  call terminate()
               end if

               nfe_assert(cv_spacing(n).gt.ZERO)
               cv_extents(n) = 1 &
                  + int((tmp - cv_origin(n))/cv_spacing(n))

               if (cv_extents(n).lt.UMBRELLA_MIN_EXTENT) then
                  write (unit = ERR_UNIT, fmt = '(/a,a,i1,a/)') NFE_ERROR, &
                     'CV #', n, ' : the ''resolution'' is too big'
                  call terminate()
               end if

               cv_spacing(n) = (tmp - cv_origin(n))/(cv_extents(n) - 1)

            end if ! cv_periodicity(n)
         end if ! imode.eq.MODE_FLOODING
         if (colvar_has_refcrd(self%colvars(n))) then
            refcrd_file = trim(refcrd_file)
            refcrd_len = len_trim(refcrd_file)
         end if
         if (colvar_is_quaternion(self%colvars(n))) then 
          allocate(self%colvars(n)%q_index, stat = error)
           if (error.ne.0) &
            NFE_OUT_OF_MEMORY
            self%colvars(n)%q_index = q_index
         end if 
         if (colvar_has_axis(self%colvars(n))) then
           allocate(self%colvars(n)%axis(3), stat = error)
             if (error.ne.0) &
                NFE_OUT_OF_MEMORY
           i = 1
           do while (i.le.3)
             self%colvars(n)%axis(i) = axis(i)
             i = i + 1
           end do
         end if
         n = n + 1

         nfe_atm_cnt = nfe_atm_cnt + cv_ni
         do i=1,cv_ni
            call AddToList(nfe_atm_lst,cv_i(i))
         end do
   end do

   ! monitor
   self%monitor_file = monitor_file

   self%monitor_freq = monitor_freq

   self%monitor_freq = min(self%monitor_freq, pmemd_nstlim())
   self%monitor_freq = max(1, self%monitor_freq)

   ! umbrella snapshots
   self%snapshots_basename = snapshots_basename

   self%snapshots_freq = snapshots_freq

   if (self%imode.eq.MODE_FLOODING) then
      self%timescale = timescale
      if (self%timescale.eq.ZERO) &
         call fatal('timescale cannot be zero !')

   end if ! imode.eq.MODE_FLOODING

1  continue ! mytaskid.gt.0 jumps here

   call mpi_bcast(self%imode, 1, MPI_INTEGER, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)

   if (self%imode.eq.MODE_NONE) &
      goto 2

   nfe_assert(.not.is_master().or.self%ncolvars.gt.0)

   call mpi_bcast(self%ncolvars, 1, MPI_INTEGER, 0, pmemd_comm, error)
   nfe_assert(error.eq.0) 

   call mpi_bcast(refcrd_len, 1, MPI_INTEGER, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
   call mpi_bcast(refcrd_file, refcrd_len, MPI_CHARACTER, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)


   call mpi_bcast(self%monitor_freq, 1, MPI_INTEGER, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)

   nfe_assert(self%ncolvars.gt.0)
   nfe_assert(self%ncolvars.le.MAX_NUMBER_OF_COLVARS)

   if (pmemd_imin().ne.0) &
      call fatal('imin.ne.0 is not supported')

   do n = 1, self%ncolvars
      call colvar_bootstrap(self%colvars(n), n, amass)
   end do

   if (mytaskid.gt.0) &
      goto 2

   do_transfer = .false.

   if (self%imode.eq.MODE_UMBRELLA) then
      nfe_assert(umbrella_file_exists)
      do_transfer = .false.
      call umbrella_swap(self%umbrella, umbrella_from_file)
   else if (self%imode.eq.MODE_FLOODING) then
      if (umbrella_file_exists) then
         do_transfer = .false.
         do n = 1, self%ncolvars
            do_transfer = do_transfer &
               .or.(cv_extents(n).ne.umbrella_extent(umbrella_from_file, n))
            do_transfer = do_transfer &
               .or.(cv_periodicity(n).neqv.&
                  umbrella_periodicity(umbrella_from_file, n))
            do_transfer = do_transfer &
               .or.(abs(cv_origin(n) - umbrella_origin(umbrella_from_file, n)) &
                  .gt.TINY)
            do_transfer = do_transfer &
               .or.(abs(cv_spacing(n) - umbrella_spacing(umbrella_from_file, &
                  n)).gt.TINY)
            if (do_transfer) &
               exit
         end do
         if (do_transfer) then
            call umbrella_init(self%umbrella, self%ncolvars, cv_extents, &
                               cv_origin, cv_spacing, cv_periodicity)
            call umbrella_transfer(self%umbrella, umbrella_from_file)
            call umbrella_fini(umbrella_from_file)
         else
            call umbrella_swap(self%umbrella, umbrella_from_file)
         end if ! do_transfer
      else
         call umbrella_init(self%umbrella, self%ncolvars, cv_extents, &
                            cv_origin, cv_spacing, cv_periodicity)
      end if ! umbrella_file_exits

      call umbrella_init(self%wt_umbrella, self%ncolvars, cv_extents, &
                            cv_origin, cv_spacing, cv_periodicity)
   end if ! self%imode.eq.MODE_FLOODING

   self%umbrella_file_existed = umbrella_file_exists
   self%umbrella_discretization_changed = do_transfer

   ! prepare monitor_fmt & open MONITOR_UNIT

   open (unit = MONITOR_UNIT, file = self%monitor_file, iostat = error, &
         form = 'FORMATTED', action = 'WRITE', status = 'REPLACE')

   if (error.ne.0) then
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
         NFE_ERROR, 'failed to open ''', trim(self%monitor_file), &
         ''' for writing'
      call terminate()
   end if

   write (unit = MONITOR_UNIT, fmt = '(a,/a)', advance = 'NO') &
      '#', '# MD time (ps), '
   do n = 1, self%ncolvars - 1
      write (unit = MONITOR_UNIT, fmt = '(a,i1,a)', advance = 'NO') &
         'CV #', n, ', '
   end do

   if (self%imode == MODE_FLOODING) then
      write (unit = MONITOR_UNIT, fmt = '(a,i1,a,/a)') &
         'CV #', self%ncolvars, ', E_{bias} (kcal/mol)', '#'
      write (unit = self%monitor_fmt, fmt = '(a,i1,a)') &
         '(f12.4,', self%ncolvars, '(1x,f16.10),1x,f16.10)'
   else
      write (unit = MONITOR_UNIT, fmt = '(a,i1,/a)') &
         'CV #', self%ncolvars, '#'
      write (unit = self%monitor_fmt, fmt = '(a,i1,a)') &
         '(f12.4,', self%ncolvars, '(1x,f16.10))'
   end if

   call flush_UNIT(MONITOR_UNIT, self%monitor_file)

2  continue ! mytaskid.gt.0 jump here (or imode.eq.MODE_NONE)

#  ifndef NFE_DISABLE_ASSERT
   self%initialized = .true.
#  endif /* NFE_DISABLE_ASSERT */

end subroutine ctxt_init

!-------------------------------------------------------------------------------

subroutine ctxt_fini(self)

   use nfe_colvar_mod
   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self

   integer :: n

   nfe_assert(self%initialized)

   if (self%imode.ne.MODE_NONE) then
      nfe_assert(self%ncolvars.gt.0)
      do n = 1, self%ncolvars
         call colvar_cleanup(self%colvars(n))
      end do
   end if

   NFE_MASTER_ONLY_BEGIN
      if (self%imode.eq.MODE_FLOODING.or.self%imode.eq.MODE_UMBRELLA) &
         call umbrella_fini(self%umbrella)
   NFE_MASTER_ONLY_END

   self%imode = MODE_NONE

#  ifndef NFE_DISABLE_ASSERT
   self%initialized = .false.
#  endif /* NFE_DISABLE_ASSERT */

end subroutine ctxt_fini

!-------------------------------------------------------------------------------

subroutine ctxt_print(self, lun)

   use nfe_lib_mod
   use nfe_colvar_mod

   implicit none

   type(bbmd_ctxt_t), intent(in) :: self
   integer, intent(in) :: lun

   integer :: n
   double precision :: tmp

   nfe_assert(self%initialized)

   write (unit = lun, fmt = '(a,a)', advance = 'NO') NFE_INFO, 'mode = '

   select case(self%imode)
      case(MODE_NONE)
         write (unit = lun, fmt = '(a)') 'NONE'
         goto 1
      case(MODE_ANALYSIS)
         write (unit = lun, fmt = '(a)') 'ANALYSIS'
      case(MODE_UMBRELLA)
         write (unit = lun, fmt = '(a)') 'UMBRELLA'
      case(MODE_FLOODING)
         write (unit = lun, fmt = '(a)') 'FLOODING'
      case default
         nfe_assert_not_reached()
         continue
   end select

   write (unit = lun, fmt = '(a)') NFE_INFO

   do n = 1, self%ncolvars
      write (unit = lun, fmt = '(a,a,i1)') NFE_INFO, 'CV #', n
      call colvar_print(self%colvars(n), lun)
      write (unit = lun, fmt = '(a)') NFE_INFO
      if (colvar_is_quaternion(self%colvars(n))) then
         write (unit = OUT_UNIT, fmt = '(a,a,I3)') NFE_INFO, &
         ' q_index = ', self%colvars(n)%q_index
      end if
      if (colvar_has_axis(self%colvars(n))) then
         write (unit = OUT_UNIT, fmt = '(a,a,f8.4,a,f8.4,a,f8.4,a)') NFE_INFO, &
         ' axis = [',self%colvars(n)%axis(1),', ', self%colvars(n)%axis(2), &
                ', ', self%colvars(n)%axis(3),']'
      end if

   end do

   write (unit = lun, fmt = '(a,a,a)') NFE_INFO, &
      'monitor_file = ', trim(self%monitor_file)
   write (unit = lun, fmt = '(a,a,'//pfmt &
      (self%monitor_freq)//',a,'//pfmt &
      (self%monitor_freq*pmemd_timestep(), 4)//',a)') NFE_INFO, &
      'monitor_freq = ', self%monitor_freq, ' (', &
      self%monitor_freq*pmemd_timestep(), ' ps)'

   if (self%imode.eq.MODE_ANALYSIS) &
      goto 1

   write (unit = lun, fmt = '(a,a,a,a)', advance = 'NO') NFE_INFO, &
      'umbrella_file = ', trim(self%umbrella_file), ' ('

   if (self%umbrella_file_existed) then
      write (unit = lun, fmt = '(a)') 'loaded)'
   else
      write (unit = lun, fmt = '(a)') 'not found)'
   end if

   write (unit = lun, fmt = '(a)') NFE_INFO
   write (unit = lun, fmt = '(a,a)', advance = 'NO') NFE_INFO, &
      'umbrella discretization '

   if (self%umbrella_file_existed) then
      if (self%umbrella_discretization_changed) then
         write (unit = lun, fmt = '(a)') '(modified) :'
      else
         write (unit = lun, fmt = '(a)') '(unchanged) :'
      end if
   else
      write (unit = lun, fmt = '(a)') '(new) :'
   end if

   do n = 1, self%ncolvars
      write (unit = lun, fmt = '(a,a,i1)', advance = 'NO') &
         NFE_INFO, 'CV #', n
      if (umbrella_periodicity(self%umbrella, n)) then
         write (unit = lun, fmt = '(a)', advance = 'NO') ' periodic, '
         tmp = umbrella_origin(self%umbrella, n) &
         + umbrella_spacing(self%umbrella, n)*umbrella_extent(self%umbrella, n)
      else
         write (unit = lun, fmt = '(a)', advance = 'NO') ' not periodic, '
         tmp = umbrella_origin(self%umbrella, n) &
         + umbrella_spacing(self%umbrella, n)*(umbrella_extent(self%umbrella, n) - 1)
      end if

      write (unit = lun, &
         fmt = '('//pfmt(umbrella_extent(self%umbrella, n))//',a,'//pfmt &
         (umbrella_origin(self%umbrella, n), 6)//',a,'//pfmt(tmp, 6)//')') &
         umbrella_extent(self%umbrella, n), ' points, min/max = ', &
         umbrella_origin(self%umbrella, n), '/', tmp
   end do

   if (self%imode.eq.MODE_UMBRELLA) &
      goto 1

   write (unit = lun, fmt = '(a/,a,a,'//pfmt(self%timescale, 3)//',a)') &
      NFE_INFO, NFE_INFO, 'flooding timescale = ', self%timescale, ' ps'

   if (self%snapshots_freq.gt.0) then
      write (unit = lun, fmt = '(a,a,a)') NFE_INFO, &
         'snapshots_basename = ', trim(self%snapshots_basename)
      write (unit = lun, &
         fmt = '(a,a,'//pfmt(self%snapshots_freq)//',a,'//pfmt &
         (self%snapshots_freq*pmemd_timestep(), 4)//',a)') NFE_INFO, &
         'snapshots_freq = ', self%snapshots_freq, ' (', &
         self%snapshots_freq*pmemd_timestep(), ' ps)'
   end if

! Modified by M Moradi
! Well-tempered ABMD
   if (self%pseudo.gt.ZERO) then
      write (unit = lun, &
      fmt = '(a/,a,a)')&
          NFE_INFO, NFE_INFO,'well-tempered ABMD:'
      write (unit = lun, &
      fmt = '(a,a,'//pfmt(ONE/self%pseudo,6)//')')&
          NFE_INFO, 'pseudo-temperature = ', ONE/self%pseudo
      write (unit = lun, fmt = '(a,a,a)') NFE_INFO, &
         'wt_umbrella_file = ', trim(self%wt_umbrella_file)
   end if
! Driven ABMD
   if (self%drivenw.gt.ZERO) then
      write (unit = lun, &
      fmt = '(a/,a,a)')&
          NFE_INFO, NFE_INFO,'driven ABMD:'
      if (self%drivenu.gt.ZERO) then
         write (unit = lun, &
         fmt = '(a,a)')&
             NFE_INFO,'CONSTANT weighting scheme (use delta work)'
      else
         write (unit = lun, &
         fmt = '(a,a)')&
             NFE_INFO,'PULLING weighting scheme (use work only)'
      end if
      write (unit = lun, &
      fmt = '(a,a,'//pfmt(self%driven_cutoff,6)//')')&
          NFE_INFO, 'driven (delta)work cutoff = ', self%driven_cutoff
   end if
! Moradi end

1  call flush_UNIT(lun, self%mdout)

end subroutine ctxt_print

!-------------------------------------------------------------------------------

subroutine ctxt_bcast(self, masterroot, amass)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self
   integer, intent(in) :: masterroot
   double precision, intent(in) :: amass(*)

   integer :: n, error

   nfe_assert(masterroot.ge.0.and.masterroot.lt.master_size)

#  ifndef NFE_DISABLE_ASSERT
   if (masterroot.eq.master_rank) then
      nfe_assert(self%initialized)
   else
      self%initialized = .true.
   end if ! masterroot.eq.master_rank
#  endif /* NFE_DISABLE_ASSERT */

   ! PMEMD's files (no matter what the mode is)

   NFE_MASTER_ONLY_BEGIN
      nfe_assert(pmemd_master_comm.ne.MPI_COMM_NULL)

      call mpi_bcast(self%mdout, len(self%mdout), MPI_CHARACTER, &
                     masterroot, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%restrt, len(self%restrt), MPI_CHARACTER, &
                     masterroot, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%mdvel, len(self%mdvel), MPI_CHARACTER, &
                     masterroot, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%mden, len(self%mden), MPI_CHARACTER, &
                     masterroot, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%mdcrd, len(self%mdcrd), MPI_CHARACTER, &
                     masterroot, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%mdinfo, len(self%mdinfo), MPI_CHARACTER, &
                     masterroot, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%ioutfm, 1, MPI_INTEGER, masterroot, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%ntpr, 1, MPI_INTEGER, masterroot, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%ntwr, 1, MPI_INTEGER, masterroot, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%ntwx, 1, MPI_INTEGER, masterroot, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%imode, 1, MPI_INTEGER, masterroot, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

   NFE_MASTER_ONLY_END

   if (master_rank.ne.masterroot) then
      call mpi_bcast(self%imode, 1, MPI_INTEGER, 0, pmemd_comm, error)
      nfe_assert(error.eq.0)
   end if ! master_rank.ne.masterroot

   if (self%imode.eq.MODE_NONE) &
      return

   ! CVs

   nfe_assert(self%ncolvars.gt.0.or.master_rank.ne.masterroot)

   NFE_MASTER_ONLY_BEGIN
      call mpi_bcast(self%ncolvars, 1, &
                     MPI_INTEGER, masterroot, pmemd_master_comm, error)
      nfe_assert(error.eq.0)
   NFE_MASTER_ONLY_END

   if (master_rank.ne.masterroot) then
      call mpi_bcast(self%ncolvars, 1, MPI_INTEGER, 0, pmemd_comm, error)
      nfe_assert(error.eq.0)
   end if ! master_rank.ne.masterroot

   nfe_assert(self%ncolvars.gt.0)

   do n = 1, self%ncolvars
      call bcast_colvar(self%colvars(n), n + 10*master_rank)
   end do

   ! mode = ANALYSIS

   NFE_MASTER_ONLY_BEGIN
      call mpi_bcast(self%monitor_file, len(self%monitor_file), &
                     MPI_CHARACTER, masterroot, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%monitor_fmt, len(self%monitor_fmt), &
                     MPI_CHARACTER, masterroot, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%monitor_freq, 1, &
                     MPI_INTEGER, masterroot, pmemd_master_comm, error)
      nfe_assert(error.eq.0)
   NFE_MASTER_ONLY_END

   if (master_rank.ne.masterroot) then
      call mpi_bcast(self%monitor_freq, 1, MPI_INTEGER, 0, pmemd_comm, error)
      nfe_assert(error.eq.0)
   end if ! master_rank.ne.masterroot

   if (self%imode.eq.MODE_ANALYSIS) &
      return

   ! mode = UMBRELLA | FLOODING (these are on masters only)

   NFE_MASTER_ONLY_BEGIN
      call mpi_bcast(self%umbrella_file, len(self%umbrella_file), &
                     MPI_CHARACTER, masterroot, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%snapshots_basename, len(self%snapshots_basename), &
                     MPI_CHARACTER, masterroot, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%snapshots_freq, 1, &
                     MPI_INTEGER, masterroot, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      call mpi_bcast(self%timescale, 1, &
                     MPI_DOUBLE_PRECISION, masterroot, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      call umbrella_bcast(self%umbrella, pmemd_master_comm, masterroot)

      if (self%imode.eq.MODE_FLOODING) then
         call mpi_bcast(self%wt_umbrella_file, len(self%wt_umbrella_file), &
                     MPI_CHARACTER, masterroot, pmemd_master_comm, error)
         nfe_assert(error.eq.0)

         call mpi_bcast(self%pseudo, 1, &
                     MPI_DOUBLE_PRECISION, masterroot, pmemd_master_comm, error)
         nfe_assert(error.eq.0)

         call mpi_bcast(self%drivenw, 1, &
                     MPI_INTEGER, masterroot, pmemd_master_comm, error)
         nfe_assert(error.eq.0)

         call mpi_bcast(self%drivenu, 1, &
                     MPI_INTEGER, masterroot, pmemd_master_comm, error)
         nfe_assert(error.eq.0)

         call mpi_bcast(self%driven_cutoff, 1, &
                     MPI_DOUBLE_PRECISION, masterroot, pmemd_master_comm, error)
         nfe_assert(error.eq.0)

         call umbrella_bcast(self%wt_umbrella, pmemd_master_comm, masterroot)
      end if

   NFE_MASTER_ONLY_END

contains

subroutine bcast_colvar(cv, cvno)

   use nfe_colvar_mod, only : colvar_bootstrap

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer, intent(in) :: cvno

   integer :: bcastdata(3)

   NFE_MASTER_ONLY_BEGIN

      if (master_rank.eq.masterroot) then

         nfe_assert(cv%type.gt.0)

         bcastdata(1) = cv%type

         bcastdata(2) = 0
         if (associated(cv%i)) &
            bcastdata(2) = size(cv%i)

         bcastdata(3) = 0
         if (associated(cv%r)) &
            bcastdata(3) = size(cv%r)

      end if ! master_rank.eq.masterroot

      call mpi_bcast(bcastdata, 3, MPI_INTEGER, masterroot, pmemd_master_comm, error)
      nfe_assert(error.eq.0)

      if (master_rank.ne.masterroot) then
         cv%type = bcastdata(1)

         if (bcastdata(2).gt.0) then
            allocate(cv%i(bcastdata(2)), stat = error)
            if (error.ne.0) &
               NFE_OUT_OF_MEMORY
         end if ! bcastdata(2).gt.0

         if (bcastdata(3).gt.0) then
            allocate(cv%r(bcastdata(3)), stat = error)
            if (error.ne.0) &
               NFE_OUT_OF_MEMORY
         end if ! bcastdata(3).gt.0

      end if ! master_rank.ne.masterroot

      if (bcastdata(2).gt.0) then
         call mpi_bcast(cv%i, bcastdata(2), MPI_INTEGER, &
                        masterroot, pmemd_master_comm, error)
         nfe_assert(error.eq.0)
      end if ! bcastdata(2).gt.0

      if (bcastdata(3).gt.0) then
         call mpi_bcast(cv%r, bcastdata(3), MPI_DOUBLE_PRECISION, &
                        masterroot, pmemd_master_comm, error)
         nfe_assert(error.eq.0)
      end if ! bcastdata(3).gt.0

   NFE_MASTER_ONLY_END

   if (master_rank.ne.masterroot) &
      call colvar_bootstrap(cv, cvno, amass)

end subroutine bcast_colvar

end subroutine ctxt_bcast

!-------------------------------------------------------------------------------

! Modified by M Moradi
! Driven ABMD
subroutine ctxt_on_force(self, x, f, mdstep, wdriven, udriven, pot)
! Moradi end

   use nfe_lib_mod
   use nfe_colvar_mod
   use parallel_dat_mod

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self

   double precision, intent(in) :: x(3,*)

   double precision, intent(inout) :: f(3,*)

   integer, intent(in) :: mdstep

! Modified by M Moradi
! Driven ABMD
   double precision, intent(in) :: wdriven
   double precision, intent(in) :: udriven
! Moradi end
   double precision, intent(inout) :: pot

   double precision :: u_derivative(UMBRELLA_MAX_NEXTENTS)
   double precision :: instantaneous(UMBRELLA_MAX_NEXTENTS)
   double precision :: alt, u_value

   character(len = SL + 16) :: snapshot

   integer :: n, error, m 
   logical :: real_mdstep
   integer, DIMENSION(4) :: cv_q = (/COLVAR_QUATERNION0, COLVAR_QUATERNION1, &
                                     COLVAR_QUATERNION2, COLVAR_QUATERNION3/)
   double precision :: norm4(100), cv_N(4, 100)

   nfe_assert(self%initialized)

   if (self%imode.eq.MODE_NONE) &
      return

   real_mdstep = (pmemd_init().eq.1)  !feng

   nfe_assert(self%ncolvars.gt.0)

   if (self%imode.eq.MODE_ANALYSIS) then
      if (nfe_real_mdstep.and.mod(mdstep, self%monitor_freq).eq.0) then
         do n = 1, self%ncolvars
            instantaneous(n) = colvar_value(self%colvars(n), x)
         end do

         NFE_MASTER_ONLY_BEGIN
         do n = 1, self%ncolvars
          if (colvar_is_quaternion(self%colvars(n))) then
           do m = 1, 4
            if (self%colvars(n)%type == cv_q(m)) then
              cv_N(m, self%colvars(n)%q_index) = instantaneous(n)
            end if
           end do
          end if
         end do
         do n = 1, self%ncolvars
          if (colvar_is_quaternion(self%colvars(n))) then
            norm4(self%colvars(n)%q_index) =sqrt(cv_N(1,self%colvars(n)%q_index)**2 &
                                               + cv_N(2,self%colvars(n)%q_index)**2 &
                                               + cv_N(3,self%colvars(n)%q_index)**2 &
                                               + cv_N(4,self%colvars(n)%q_index)**2)
          end if
         end do
         do n = 1, self%ncolvars
          if (colvar_is_quaternion(self%colvars(n))) then
             instantaneous(n) = instantaneous(n) / norm4(self%colvars(n)%q_index)
          else
             instantaneous(n) = instantaneous(n)
          end if
         end do

         write (unit = MONITOR_UNIT, fmt = self%monitor_fmt) &
               pmemd_mdtime(), instantaneous(1:self%ncolvars)
            call flush_UNIT(MONITOR_UNIT, self%monitor_file)
         NFE_MASTER_ONLY_END
      end if

      return
   end if ! self%imode.eq.MODE_ANALYSIS

   !
   ! either UMBRELLA or FLOODING
   !

   do n = 1, self%ncolvars
      instantaneous(n) = colvar_value(self%colvars(n), x)
   end do

   NFE_MASTER_ONLY_BEGIN
   do n = 1, self%ncolvars
    if (colvar_is_quaternion(self%colvars(n))) then
      do m = 1, 4
        if (self%colvars(n)%type == cv_q(m)) then
          cv_N(m, self%colvars(n)%q_index) = instantaneous(n)
        end if
      end do
    end if
   end do
   do n = 1, self%ncolvars
    if (colvar_is_quaternion(self%colvars(n))) then
       norm4(self%colvars(n)%q_index) =sqrt(cv_N(1,self%colvars(n)%q_index)**2 &
                                          + cv_N(2,self%colvars(n)%q_index)**2 &
                                          + cv_N(3,self%colvars(n)%q_index)**2 &
                                          + cv_N(4,self%colvars(n)%q_index)**2)
    end if
   end do
   do n = 1, self%ncolvars
    if (colvar_is_quaternion(self%colvars(n))) then
      instantaneous(n) = instantaneous(n) / norm4(self%colvars(n)%q_index)
    else
      instantaneous(n) = instantaneous(n)
    end if
   end do

      call umbrella_eval_vdv(self%umbrella, instantaneous, &
                             u_value, u_derivative)
      pot = umbrella_eval_v(self%umbrella, instantaneous)
   NFE_MASTER_ONLY_END

   call mpi_bcast(u_derivative, self%ncolvars, &
      MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
   call mpi_bcast(pot, 1, &
      MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)

   ! FIXME: virial
   do n = 1, self%ncolvars
      call colvar_force(self%colvars(n), x, -u_derivative(n), f)
   end do

   if (.not.nfe_real_mdstep.or.mytaskid.ne.0) &
      return

   nfe_assert(self%imode.eq.MODE_UMBRELLA.or.self%imode.eq.MODE_FLOODING)

   if (mod(mdstep, self%monitor_freq).eq.0) then
      if (self%imode.eq.MODE_FLOODING) then
         write (unit = MONITOR_UNIT, fmt = self%monitor_fmt) &
            pmemd_mdtime(), instantaneous(1:self%ncolvars), u_value
      else
         write (unit = MONITOR_UNIT, fmt = self%monitor_fmt) &
            pmemd_mdtime(), instantaneous(1:self%ncolvars)
      end if
      call flush_UNIT(MONITOR_UNIT, self%monitor_file)
   end if

   if (self%imode.eq.MODE_FLOODING) then
#  ifdef BINTRAJ
      if (self%snapshots_freq.gt.0 &
          .and.mod(mdstep, self%snapshots_freq).eq.0) then
         write (unit = snapshot, fmt = '(a,a,i10.10,a)') &
            trim(self%snapshots_basename), '.', mdstep, '.nc'
         call umbrella_save(self%umbrella, snapshot)
         write (unit = OUT_UNIT, fmt = '(/a,a,'//pfmt &
            (pmemd_mdtime(), 4)//',a,/a,a,a,a)') &
            NFE_INFO, 'biasing potential snapshot at t = ', &
            pmemd_mdtime(), ' ps', NFE_INFO, 'saved as ''', &
            trim(snapshot), ''''
      end if
#  endif /* BINTRAJ */
! Modified by M Moradi
       alt = pmemd_timestep()/self%timescale
! Well-tempered ABMD (based on Barducci, Bussi, and Parrinello, PRL(2008) 100:020603)
      if (self%pseudo.gt.ZERO) &
         alt = alt * &
         exp(-self%pseudo*umbrella_eval_v(self%umbrella,instantaneous)/kB)
! Driven ABMD (based on Moradi and Tajkhorshid, JPCL(2013) 4:1882)
      if (self%drivenw.ne.0) then
         if (self%drivenw*wdriven-self%drivenu*udriven.gt.self%driven_cutoff) then
            alt = alt * &
             exp(-(self%drivenw*wdriven-self%drivenu*udriven)/(kB*pmemd_temp0()))
         else
            alt = alt * &
             exp(-self%driven_cutoff/(kB*pmemd_temp0()))
         end if
      end if ! self%drivenw.ne.0
! Moradi end
      call umbrella_hill(self%umbrella, instantaneous, alt)
   end if ! self%imode.eq.MODE_FLOODING

end subroutine ctxt_on_force

!-------------------------------------------------------------------------------

subroutine ctxt_on_mdwrit(self)

   use nfe_lib_mod

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self
   double precision :: t

   nfe_assert(self%initialized)

#ifdef BINTRAJ
   if (self%imode.eq.MODE_FLOODING) then
      nfe_assert(is_master())
      call umbrella_save(self%umbrella, self%umbrella_file)
! Modified by F Pan
      if (abs(self%pseudo-ZERO).gt.TINY) then
         t = pmemd_temp0()
         call umbrella_copy(self%wt_umbrella, self%umbrella)
         call umbrella_wt_mod(self%wt_umbrella,t,1/self%pseudo)
         call umbrella_save(self%wt_umbrella, self%wt_umbrella_file)
      end if
! Pan end
   end if
#endif /* BINTRAJ*/

end subroutine ctxt_on_mdwrit

!-------------------------------------------------------------------------------

!
! U_m? are valid for master_rank.lt.r_master_rank
!

! assumes that local self is up to date
subroutine ctxt_Um(self, r_master_rank, x, U_mm, U_mo)

   use nfe_colvar_mod
   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self

   integer, intent(in) :: r_master_rank
   double precision, intent(in) :: x(*)

   double precision, intent(out) :: U_mm, U_mo

   double precision :: local_inst(UMBRELLA_MAX_NEXTENTS)
   double precision :: remote_inst(UMBRELLA_MAX_NEXTENTS)

   integer :: n, error
   nfe_assert(self%initialized)

   U_mm = ZERO
   U_mo = ZERO

   if (self%imode.eq.MODE_NONE.or.self%imode.eq.MODE_ANALYSIS) &
      return

   nfe_assert(self%ncolvars.gt.0)
   do n = 1, self%ncolvars
      local_inst(n) = colvar_value(self%colvars(n), x)
   end do

   if (mytaskid.gt.0) &
      return

   nfe_assert(self%imode.eq.MODE_UMBRELLA.or.self%imode.eq.MODE_FLOODING)

   if (master_rank.lt.r_master_rank) then
      call mpi_recv(remote_inst, self%ncolvars, MPI_DOUBLE_PRECISION, &
                    r_master_rank, 0, pmemd_master_comm, MPI_STATUS_IGNORE, error)
      nfe_assert(error.eq.0)
      U_mm = umbrella_eval_v(self%umbrella, local_inst)
      U_mo = umbrella_eval_v(self%umbrella, remote_inst)
   else
      call mpi_send(local_inst, self%ncolvars, MPI_DOUBLE_PRECISION, &
                    r_master_rank, 0, pmemd_master_comm, error)
      nfe_assert(error.eq.0)
   end if ! master_rank.lt.r_master_rank

end subroutine ctxt_Um

! assumes that remote self is up to date (if mode.eq.MODE_FLOODING)
subroutine ctxt_Uo(self, r_master_rank, x, U_om, U_oo)

   use nfe_colvar_mod
   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self

   integer, intent(in) :: r_master_rank
   double precision, intent(in) :: x(*)

   double precision, intent(out) :: U_om, U_oo

   double precision :: local_inst(UMBRELLA_MAX_NEXTENTS)
   double precision :: remote_inst(UMBRELLA_MAX_NEXTENTS)

   double precision :: tmp(2)

   integer :: n, error

   nfe_assert(self%initialized)

   U_om = ZERO
   U_oo = ZERO

   if (self%imode.eq.MODE_NONE.or.self%imode.eq.MODE_ANALYSIS) &
      return

   nfe_assert(self%ncolvars.gt.0)
   do n = 1, self%ncolvars
      local_inst(n) = colvar_value(self%colvars(n), x)
   end do

   if (mytaskid.gt.0) &
      return

   if (self%imode.eq.MODE_UMBRELLA) then
      ! remote umbrella is same as local
      if (master_rank.lt.r_master_rank) then
         call mpi_recv(remote_inst, self%ncolvars, MPI_DOUBLE_PRECISION, &
                       r_master_rank, 0, pmemd_master_comm, MPI_STATUS_IGNORE, error)
         nfe_assert(error.eq.0)
         U_om = umbrella_eval_v(self%umbrella, local_inst)
         U_oo = umbrella_eval_v(self%umbrella, remote_inst)
      else
         call mpi_send(local_inst, self%ncolvars, MPI_DOUBLE_PRECISION, &
                       r_master_rank, 0, pmemd_master_comm, error)
         nfe_assert(error.eq.0)
      end if ! master_rank.lt.r_master_rank
   else
      nfe_assert(self%imode.eq.MODE_FLOODING)
      if (master_rank.gt.r_master_rank) then
         call mpi_recv(remote_inst, self%ncolvars, MPI_DOUBLE_PRECISION, &
                       r_master_rank, 0, pmemd_master_comm, MPI_STATUS_IGNORE, error)
         nfe_assert(error.eq.0)
         tmp(1) = umbrella_eval_v(self%umbrella, local_inst)
         tmp(2) = umbrella_eval_v(self%umbrella, remote_inst)
         call mpi_send(tmp, 2, MPI_DOUBLE_PRECISION, &
                       r_master_rank, 1, pmemd_master_comm, error)
         nfe_assert(error.eq.0)
      else
         call mpi_send(local_inst, self%ncolvars, MPI_DOUBLE_PRECISION, &
                       r_master_rank, 0, pmemd_master_comm, error)
         nfe_assert(error.eq.0)
         call mpi_recv(tmp, 2, MPI_DOUBLE_PRECISION, &
                       r_master_rank, 1, pmemd_master_comm, MPI_STATUS_IGNORE, error)
         nfe_assert(error.eq.0)
         U_om = tmp(2)
         U_oo = tmp(1)
      end if ! master_rank.gt.r_master_rank
   end if ! self%imode.eq.MODE_UMBRELLA

end subroutine ctxt_Uo

!-------------------------------------------------------------------------------

subroutine ctxt_send(self, dst_master_rank)

   NFE_USE_AFAILED

   use nfe_lib_mod, only : umbrella_send_coeffs
   use parallel_dat_mod

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self
   integer, intent(in) :: dst_master_rank

   integer :: error

   nfe_assert(mytaskid.eq.0)
   nfe_assert(self%initialized)

   if (self%imode.eq.MODE_FLOODING) then
      call umbrella_send_coeffs(self%umbrella, dst_master_rank, pmemd_master_comm)
   else
      ! for synchronization
      call mpi_send(master_rank, 1, MPI_INTEGER, &
                    dst_master_rank, 8, pmemd_master_comm, error)
      nfe_assert(error.eq.0)
   end if

end subroutine ctxt_send

!-------------------------------------------------------------------------------

subroutine ctxt_recv(self, src_master_rank)

   NFE_USE_AFAILED

   use nfe_lib_mod, only : umbrella_recv_coeffs
   use parallel_dat_mod

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self
   integer, intent(in) :: src_master_rank

   integer :: error, itmp

   nfe_assert(mytaskid.eq.0)
   nfe_assert(self%initialized)

   if (self%imode.eq.MODE_FLOODING) then
      call umbrella_recv_coeffs(self%umbrella, src_master_rank, pmemd_master_comm)
   else
      ! for synchronization
      call mpi_recv(itmp, 1, MPI_INTEGER, src_master_rank, 8, &
                    pmemd_master_comm, MPI_STATUS_IGNORE, error)
      nfe_assert(error.eq.0)
      nfe_assert(itmp.eq.src_master_rank)
   end if

end subroutine ctxt_recv

!-------------------------------------------------------------------------------

subroutine ctxt_close_units(self)

   use nfe_lib_mod
   use parallel_dat_mod
   use file_io_dat_mod
   use mdin_ctrl_dat_mod
   use bintraj_mod, only : close_binary_files
   use get_cmdline_mod, only : cpein_specified

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self

   nfe_assert(self%initialized)
   nfe_assert(mytaskid.eq.0)

   if (self%imode.ne.MODE_NONE) &
      close (MONITOR_UNIT)

   if (self%mdout.ne.'stdout') &
      close (OUT_UNIT)

!   call close_dump_files()
   if (master) then
       if ( ioutfm.eq.0) then
          if ( ntwx > 0 ) close( mdcrd )
          if ( ntwv > 0 ) close( mdvel )
          if ( ntwf > 0 ) close( mdfrc )
       else
          call close_binary_files
       end if
       if ( ntwe > 0 ) close( mden )
       if ( ntpr > 0 ) close( mdinfo )
       if ( icnstph /= 0 .and. .not. cpein_specified) close ( cpout )
       if ( icnste /= 0  .and. .not. cpein_specified) close ( ceout )
       if ( (icnstph /= 0 .or. icnste /= 0)  .and. cpein_specified) close ( cpeout )
   end if

end subroutine ctxt_close_units

!-------------------------------------------------------------------------------

subroutine ctxt_open_units(self)

   use nfe_lib_mod
   use file_io_dat_mod
   use mdin_ctrl_dat_mod
   use parallel_dat_mod
   use prmtop_dat_mod, only : prmtop_ititl
   use file_io_mod, only : amopen
   use bintraj_mod, only : open_binary_files_append
   use get_cmdline_mod, only : cpein_specified

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self

   integer :: error

   nfe_assert(self%initialized)
   nfe_assert(mytaskid.eq.0)

   if (self%imode.ne.MODE_NONE) then
      open (unit = MONITOR_UNIT, file = self%monitor_file, iostat = error, &
        form = 'FORMATTED', action = 'WRITE', status = 'OLD', position = 'APPEND')
      if (error.ne.0) then
         write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
            NFE_ERROR, 'failed to open ''', trim(self%monitor_file), &
            ''' for appending'
         call terminate()
      end if
   end if ! self%imode.ne.MODE_NONE

   mdout_name = self%mdout
   mdinfo_name = self%mdinfo
   restrt_name = self%restrt
   mdvel_name = self%mdvel
   mden_name = self%mden
   mdcrd_name = self%mdcrd

   ioutfm = self%ioutfm

   ntpr = self%ntpr
   ntwr = self%ntwr
   ntwx = self%ntwx

   if (self%mdout.ne.'stdout') &
      call amopen(OUT_UNIT, mdout_name, 'O', 'F', 'A')

   call amopen(mdinfo, mdinfo_name, 'U', 'F', 'W')

   owrite = 'U'
!   call open_dump_files()
   if (master) then
       if (ioutfm.eq.0) then
           if (ntwx > 0) &
              call amopen(mdcrd,mdcrd_name,'U','F','A')
           if (ntwv > 0) &
              call amopen(mdvel,mdvel_name,owrite,'F','A')
           if (ntwf > 0) &
              call amopen(mdfrc,mdfrc_name,owrite,'F','A')
       else
           call open_binary_files_append
       end if

       if (icnstph /= 0 .and. .not. cpein_specified) &
           call amopen(cpout, cpout_name, owrite, 'F', 'A')
       if (icnste /= 0 .and. .not. cpein_specified) &
           call amopen(ceout, ceout_name, owrite, 'F', 'A')
       if ((icnstph /= 0 .and. icnste /= 0) .and. cpein_specified) &
           call amopen(cpeout, ceout_name, owrite, 'F', 'A')
       if (ntwe > 0) &
           call amopen(mden,mden_name,owrite,'F','A')
   end if
   1000 format(a80)

end subroutine ctxt_open_units

!-------------------------------------------------------------------------------
!<>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><>
! mt19937 subroutines
!<>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><>
  elemental function uiadd( a, b ) result( c )

    implicit none

    intrinsic :: ibits, ior, ishft

    integer( kind = wi ), intent( in )  :: a, b

    integer( kind = wi )  :: c

    integer( kind = wi )  :: a1, a2, b1, b2, s1, s2

    a1 = ibits( a, 0, hbs )
    a2 = ibits( a, hbs, hbs )
    b1 = ibits( b, 0, hbs )
    b2 = ibits( b, hbs, hbs )
    s1 = a1 + b1
    s2 = a2 + b2 + ibits( s1, hbs, hbs )
    c  = ior( ishft( s2, hbs ), ibits( s1, 0, hbs ) )

  end function uiadd

  elemental function uisub( a, b ) result( c )

    implicit none

    intrinsic :: ibits, ior, ishft

    integer( kind = wi ), intent( in )  :: a, b

    integer( kind = wi )  :: c

    integer( kind = wi )  :: a1, a2, b1, b2, s1, s2

    a1 = ibits( a, 0, hbs )
    a2 = ibits( a, hbs, hbs )
    b1 = ibits( b, 0, hbs )
    b2 = ibits( b, hbs, hbs )
    s1 = a1 - b1
    s2 = a2 - b2 + ibits( s1, hbs, hbs )
    c  = ior( ishft( s2, hbs ), ibits( s1, 0, hbs ) )

  end function uisub

  elemental function uimlt( a, b ) result( c )

    implicit none

    intrinsic :: ibits, ior, ishft

    integer( kind = wi ), intent( in )  :: a, b

    integer( kind = wi )  :: c

    integer( kind = wi )  :: a0, a1, a2, a3
    integer( kind = wi )  :: b0, b1, b2, b3
    integer( kind = wi )  :: p0, p1, p2, p3

    a0 = ibits( a, 0, qbs )
    a1 = ibits( a, qbs, qbs )
    a2 = ibits( a, hbs, qbs )
    a3 = ibits( a, tbs, qbs )
    b0 = ibits( b, 0, qbs )
    b1 = ibits( b, qbs, qbs )
    b2 = ibits( b, hbs, qbs )
    b3 = ibits( b, tbs, qbs )
    p0 = a0 * b0
    p1 = a1 * b0 + a0 * b1 + ibits( p0, qbs, tbs )
    p2 = a2 * b0 + a1 * b1 + a0 * b2 + ibits( p1, qbs, tbs )
    p3 = a3 * b0 + a2 * b1 + a1 * b2 + a0 * b3 + ibits( p2, qbs, tbs )
    c  = ior( ishft( p1, qbs ), ibits( p0, 0, qbs ) )
    c  = ior( ishft( p2, hbs ), ibits( c, 0, hbs ) )
    c  = ior( ishft( p3, tbs ), ibits( c, 0, tbs ) )

  end function uimlt

  ! initializes mt[N] with a seed
  subroutine init_by_seed(self, s)

    implicit none

    intrinsic :: iand, ishft, ieor, ibits

    type(mt19937_t), intent(inout) :: self
    integer( kind = wi ), intent( in )  :: s

    integer( kind = wi )  :: i, mult_a

#ifndef ASSUME_GFORTRAN
    data mult_a /z'6C078965'/ ! gfortran does not like this
#else
    mult_a = ieor(ishft(int(z'6C07'), 16), z'8965') ! but this is okay
#endif /* ASSUME_GFORTRAN */

    self%mtinit = .true._wi
    self%mt(1) = ibits( s, 0, 32 )
    do i = 2, n, 1
      self%mt(i) = ieor( self%mt(i-1), ishft(self%mt(i-1), -30 ) )
      self%mt(i) = uimlt( self%mt(i), mult_a)
      self%mt(i) = uiadd( self%mt(i), uisub( i, 1_wi ) )
      ! See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.
      ! In the previous versions, MSBs of the seed affect
      ! only MSBs of the array mt[].
      ! 2002/01/09 modified by Makoto Matsumoto
      self%mt(i) = ibits( self%mt(i), 0, 32 )
      ! for >32 bit machines
    end do
    self%mti = n + 1_wi

  end subroutine init_by_seed

  ! initialize by an array with array-length
  ! init_key is the array for initializing keys
  ! key_length is its length
  subroutine init_by_array(self, init_key)

    implicit none

    intrinsic :: iand, ishft, ieor

    type(mt19937_t), intent(inout) :: self
    integer( kind = wi ), intent( in )  :: init_key(:)

    integer( kind = wi )  :: i, j, k, tp, key_length
    integer( kind = wi )  :: seed_d, mult_a, mult_b, msb1_d

    data seed_d /z'12BD6AA'/
    data mult_a /z'19660D'/

#ifndef ASSUME_GFORTRAN
    data mult_b /z'5D588B65'/
    data msb1_d /z'80000000'/
#else
    mult_b = ieor(ishft(int(z'5D58'), 16), z'8B65')
    msb1_d = 1
    msb1_d = ishft(msb1_d, 31)
#endif /* ASSUME_GFORTRAN */

    key_length = size( init_key, dim = 1 )
    call init_by_seed(self, seed_d)
    i = 2_wi
    j = 1_wi
    do k = max( n, key_length ), 1, -1
      tp = ieor( self%mt(i-1), ishft( self%mt(i-1), -30 ) )
      tp = uimlt( tp, mult_a )
      self%mt(i) = ieor( self%mt(i), tp )
      self%mt(i) = uiadd( self%mt(i), uiadd( init_key(j), uisub( j, 1_wi ) ) ) ! non linear
      self%mt(i) = ibits( self%mt(i), 0, 32 ) ! for WORDSIZE > 32 machines
      i = i + 1_wi
      j = j + 1_wi
      if ( i > n ) then
        self%mt(1) = self%mt(n)
        i = 2_wi
      end if
      if ( j > key_length) j = 1_wi
    end do
    do k = n-1, 1, -1
      tp = ieor( self%mt(i-1), ishft( self%mt(i-1), -30 ) )
      tp = uimlt( tp, mult_b )
      self%mt(i) = ieor( self%mt(i), tp )
      self%mt(i) = uisub( self%mt(i), uisub( i, 1_wi ) ) ! non linear
      self%mt(i) = ibits( self%mt(i), 0, 32 ) ! for WORDSIZE > 32 machines
      i = i + 1_wi
      if ( i > n ) then
        self%mt(1) = self%mt(n)
        i = 2_wi
      end if
    end do
    self%mt(1) = msb1_d ! MSB is 1; assuring non-zero initial array
  end subroutine init_by_array

  ! generates a random number on [0,0xffffffff]-interval
  function random_int32(self) result( y )

    implicit none

    type(mt19937_t), intent(inout) :: self

    intrinsic :: iand, ishft, ior, ieor, btest, ibset, mvbits

    integer( kind = wi )  :: y

    integer( kind = wi )  :: kk
    integer( kind = wi )  :: seed_d
    data seed_d   /z'5489'/

#ifndef ASSUME_GFORTRAN
    integer( kind = wi)   :: matrix_a, matrix_b, temper_a, temper_b

    data matrix_a /z'9908B0DF'/
    data matrix_b /z'0'/
    data temper_a /z'9D2C5680'/
    data temper_b /z'EFC60000'/
#else
#   define matrix_a z'9908B0DF'
#   define matrix_b z'0'
#   define temper_a z'9D2C5680'
#   define temper_b z'EFC60000'
#endif /* ASSUME_GFORTRAN */

    if ( self%mti > n ) then ! generate N words at one time
      if ( .not. self%mtinit ) call init_by_seed(self, seed_d) ! if init_genrand() has not been called, a default initial seed is used
      do kk = 1, n-m, 1
        y = ibits( self%mt(kk+1), 0, 31 )
        call mvbits( self%mt(kk), 31, 1, y, 31 )
        if ( btest( y, 0 ) ) then
          self%mt(kk) = int(ieor( ieor( self%mt(kk+m), ishft( y, -1 ) ), matrix_a ))
        else
          self%mt(kk) = int(ieor( ieor( self%mt(kk+m), ishft( y, -1 ) ), matrix_b ))
        end if
      end do
      do kk = n-m+1, n-1, 1
        y = ibits( self%mt(kk+1), 0, 31 )
        call mvbits( self%mt(kk), 31, 1, y, 31 )
        if ( btest( y, 0 ) ) then
          self%mt(kk) = int(ieor( ieor( self%mt(kk+m-n), ishft( y, -1 ) ), matrix_a ))
        else
          self%mt(kk) = int(ieor( ieor( self%mt(kk+m-n), ishft( y, -1 ) ), matrix_b ))
        end if
      end do
      y = ibits( self%mt(1), 0, 31 )
      call mvbits( self%mt(n), 31, 1, y, 31 )
      if ( btest( y, 0 ) ) then
        self%mt(kk) = int(ieor( ieor( self%mt(m), ishft( y, -1 ) ), matrix_a ))
      else
        self%mt(kk) = int(ieor( ieor( self%mt(m), ishft( y, -1 ) ), matrix_b ))
      end if
      self%mti = 1_wi
    end if
    y = self%mt(self%mti)
    self%mti = self%mti + 1_wi
    ! Tempering
    y = ieor( y, ishft( y, -11) )
    y = int(ieor( y, iand( ishft( y, 7 ), temper_a ) ))
    y = int(ieor( y, iand( ishft( y, 15 ), temper_b ) ))
    y = ieor( y, ishft( y, -18 ) )

  end function random_int32

  ! generates a random number on [0,0x7fffffff]-interval
  function random_int31(self) result( i )

    implicit none

    type(mt19937_t), intent(inout) :: self

    intrinsic :: ishft

    integer( kind = wi )  :: i

    i = ishft(random_int32(self), -1)

  end function random_int31

  ! generates a random number on [0,1]-real-interval
  function random_real1(self) result( r )

    implicit none

    type(mt19937_t), intent(inout) :: self

    real( kind = wr )  :: r

    integer( kind = wi )  :: a, a1, a0

    a = random_int32(self)
    a0 = ibits( a, 0, hbs )
    a1 = ibits( a, hbs, hbs )
    r = real( a0, kind = wr ) / 4294967295.0_wr
    r = real( a1, kind = wr ) * ( 65536.0_wr / 4294967295.0_wr ) + r
    ! divided by 2^32-1

  end function random_real1

  ! generates a random number on [0,1)-real-interval
  function random_real2(self) result( r )

    implicit none

    intrinsic :: ibits

    type(mt19937_t), intent(inout) :: self

    real( kind = wr )  :: r

    integer( kind = wi )  :: a, a1, a0

    a = random_int32(self)
    a0 = ibits( a, 0, hbs )
    a1 = ibits( a, hbs, hbs )
    r = real( a0, kind = wr ) / 4294967296.0_wr
    r = real( a1, kind = wr ) / 65536.0_wr + r
    ! divided by 2^32

  end function random_real2

  ! generates a random number on (0,1)-real-interval
  function random_real3(self) result( r )

    implicit none

    type(mt19937_t), intent(inout) :: self

    real( kind = wr )  :: r

    integer( kind = wi )  :: a, a1, a0

    a = random_int32(self)
    a0 = ibits( a, 0, hbs )
    a1 = ibits( a, hbs, hbs )
    r = ( real( a0, kind = wr ) + 0.5_wr ) / 4294967296.0_wr
    r = real( a1, kind = wr ) / 65536.0_wr + r
    ! divided by 2^32

  end function random_real3

  ! generates a random number on [0,1) with 53-bit resolution
  function random_res53(self)  result( r )

    implicit none

    intrinsic :: ishft

    type(mt19937_t), intent(inout) :: self

    real( kind = wr )  :: r

    integer( kind = wi )  :: a, a0, a1
    integer( kind = wi )  :: b, b0, b1

    a = ishft(random_int32(self), -5 )
    a0 = ibits( a, 0, hbs )
    a1 = ibits( a, hbs, hbs )
    b = ishft(random_int32(self), -6 )
    b0 = ibits( b, 0, hbs )
    b1 = ibits( b, hbs, hbs )
    r = real( a1, kind = wr ) / 2048.0_wr
    r = real( a0, kind = wr ) / 134217728.0_wr + r
    r = real( b1, kind = wr ) / 137438953472.0_wr + r
    r = real( b0, kind = wr ) / 9007199254740992.0_wr + r

  end function random_res53
  ! These real versions are due to Isaku Wada, 2002/01/09 added

subroutine check_ncrc(rc, sbrtn, filename)

   use netcdf
   use nfe_lib_mod, only : ERR_UNIT, terminate

   implicit none

   integer, intent(in) :: rc
   character(*), intent(in) :: sbrtn ! this is FORTRAN, after all :-)
   character(*), intent(in) :: filename

   character(len = 80) :: errmsg

   if (rc.ne.nf90_noerr) then
      errmsg = nf90_strerror(rc)
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a,a,a/)') NFE_ERROR, &
         trim(sbrtn), '(filename=''', trim(filename), ''') : ', trim(errmsg)
      call terminate()
   end if

end subroutine check_ncrc

subroutine mt19937_save(self, filename)

   use netcdf

   implicit none

   type(mt19937_t), intent(in) :: self
   character(*), intent(in) :: filename

   character(*), parameter :: sbrtn = 'mt19937%save'

   integer :: rc, setid, imtinit

   rc = nf90_create(filename, cmode = nf90_clobber, ncid = setid)
   call check_ncrc(rc, sbrtn, filename)

   rc = nf90_put_att(setid, nf90_global, 'mt', self%mt)
   call check_ncrc(rc, sbrtn, filename)

   if (self%mtinit) then
      imtinit = 1
   else
      imtinit = 0
   end if

   rc = nf90_put_att(setid, nf90_global, 'mtinit', imtinit)
   call check_ncrc(rc, sbrtn, filename)

   rc = nf90_put_att(setid, nf90_global, 'mti', self%mti)
   call check_ncrc(rc, sbrtn, filename)

   rc = nf90_close(setid)
   call check_ncrc(rc, sbrtn, filename)

end subroutine mt19937_save

subroutine mt19937_load(self, filename)

   NFE_USE_AFAILED

   use netcdf

   implicit none

   type(mt19937_t), intent(inout) :: self
   character(*), intent(in) :: filename

   character(*), parameter :: sbrtn = 'mt19937%load'

   integer :: rc, setid, imtinit

   rc = nf90_open(filename, nf90_nowrite, setid)
   call check_ncrc(rc, sbrtn, filename)

   rc = nf90_get_att(setid, nf90_global, 'mt', self%mt)
   call check_ncrc(rc, sbrtn, filename)

   rc = nf90_get_att(setid, nf90_global, 'mtinit', imtinit)
   call check_ncrc(rc, sbrtn, filename)

   select case(imtinit)
      case (0)
         self%mtinit = .false.
      case (1)
         self%mtinit = .true.
      case default
         nfe_assert_not_reached()
         continue
   end select

   rc = nf90_get_att(setid, nf90_global, 'mti', self%mti)
   call check_ncrc(rc, sbrtn, filename)

   rc = nf90_close(setid)
   call check_ncrc(rc, sbrtn, filename)

end subroutine mt19937_load

subroutine mt19937_bcast(self, comm, root)

   use nfe_lib_mod

   implicit none

   type(mt19937_t), intent(inout) :: self
   integer, intent(in) :: comm, root

   include 'mpif.h'

   integer :: error

   nfe_assert(comm.ne.MPI_COMM_NULL)

   call mpi_bcast(self%mt, size(self%mt), MPI_INTEGER, root, comm, error)
   nfe_assert(error.eq.0)

   call mpi_bcast(self%mtinit, 1, MPI_INTEGER, root, comm, error)
   nfe_assert(error.eq.0)

   call mpi_bcast(self%mti, 1, MPI_INTEGER, root, comm, error)
   nfe_assert(error.eq.0)

end subroutine mt19937_bcast
!<>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><>

#endif /* MPI */
end module nfe_bbmd_mod
