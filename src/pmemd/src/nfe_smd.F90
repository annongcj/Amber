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
!-------------------------------------------------------------------------

module nfe_smd_mod

use nfe_lib_mod, only : SL => STRING_LENGTH, SMD_OUTPUT_UNIT, SMD_CV_UNIT
use nfe_colvar_mod

implicit none

private

public :: on_pmemd_init
public :: on_pmemd_exit

public :: on_force

!- - - - - - - - - - - - - - - - P R I V A T E - - - - - - - - - - - - - - - -

integer, private, parameter :: smd_UNIT = SMD_OUTPUT_UNIT
integer, private, parameter :: CV_UNIT = SMD_CV_UNIT

character(*), private, parameter :: DEFAULT_OUTPUT_FILE = 'nfe-smd.txt'
character(*), private, parameter :: DEFAULT_CV_FILE = 'nfe-smd-cv'

integer, private, parameter :: DEFAULT_OUTPUT_FREQ = 50

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

type, private :: colvar_t_smd

   type(colvar_t) :: parent

   double precision, pointer :: cvar_path(:) => null() ! master only
   double precision, pointer :: harm_path(:) => null() ! master only

   double precision :: inst, curr, harm
   
   integer :: harm_mode
   integer :: path_mode

end type colvar_t_smd

integer, private, save :: ncolvars = 0 ! .gt.0 means "active"
type(colvar_t_smd), private, allocatable, save :: cv(:)

double precision, private, allocatable, save :: fcv_curr(:)

double precision, private, save :: work

integer, private, save :: output_freq = 50
integer, private, save :: smd_mdstep ! = runmd.f::nstep + 1
character(len = SL), private, save :: output_file = DEFAULT_OUTPUT_FILE
character(len = SL), private, save :: cv_file = DEFAULT_CV_FILE

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

integer, private, parameter :: MODE_LINES = 4123
integer, private, parameter :: MODE_SPLINE = 3112

private :: mode_eval
private :: mode_write
private :: mode_from_string

private :: lines, spline

namelist / smd /     output_file, output_freq, cv_file

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

contains

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine on_pmemd_init(mdin_unit, amass, acrds)

   use nfe_lib_mod
   use nfe_colvar_mod
   use file_io_mod
   use parallel_dat_mod
   use pmemd_lib_mod, only : upper

   implicit none

   integer, intent(in) :: mdin_unit

   double precision, intent(in) :: amass(*)
   double precision, intent(in) :: acrds(*)

   logical :: found
   integer :: i, n, error, ifind
   character(SL) :: modestr
   character(80) :: buf

#  ifdef MPI
   nfe_assert(multipmemd_rem().eq.0)
#  endif /* MPI */

   nfe_assert(ncolvars.eq.0)

   NFE_MASTER_ONLY_BEGIN

   rewind(mdin_unit)
   call nmlsrc('smd', mdin_unit, ifind)
   
   ! no smd section
   if (ifind.eq.0) goto 1

   if (pmemd_imin().ne.0) &
      call fatal('imin.ne.0 is not supported')

   rewind(mdin_unit)
   read(mdin_unit,nml=smd,err=666)
   ncolvars = 0

   ! collective variables
   call amopen(CV_UNIT, cv_file, 'O', 'F', 'R')

   do
     call nmlsrc('colvar', CV_UNIT, ifind)
     if (ifind.eq.0) exit
     read(CV_UNIT,'(a80)') buf
     ncolvars = ncolvars + 1
   end do

   if (ncolvars.eq.0) &
      call fatal('no variable(s) in the CV file')

   allocate(cv(ncolvars), fcv_curr(ncolvars), stat = error)
   if (error /= 0) &
      NFE_OUT_OF_MEMORY

   n = 1

   do while (n.le.ncolvars)
      call colvar_nlread(CV_UNIT,cv(n)%parent)
      if (colvar_has_refcrd(cv(n)%parent)) then 
            refcrd_file = trim(refcrd_file)
            refcrd_len = len_trim(refcrd_file)
      end if
      if (colvar_is_quaternion(cv(n)%parent)) then
        allocate(cv(n)%parent%q_index, stat = error)
          if (error.ne.0) &
            NFE_OUT_OF_MEMORY
            cv(n)%parent%q_index = q_index
      end if
      if (colvar_has_axis(cv(n)%parent)) then
           allocate(cv(n)%parent%axis(3), stat = error)
             if (error.ne.0) &
                NFE_OUT_OF_MEMORY
           i = 1
           do while (i.le.3)
             cv(n)%parent%axis(i) = axis(i)
             i = i + 1
           end do
      end if 
      n = n + 1
   end do

   output_freq = min(output_freq, pmemd_nstlim())
   output_freq = max(1, output_freq)

1  continue

   NFE_MASTER_ONLY_END

#  ifdef MPI
   call mpi_bcast(ncolvars, 1, MPI_INTEGER, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)

   call mpi_bcast(refcrd_len, 1, MPI_INTEGER, 0, pmemd_comm, error) 
   nfe_assert(error.eq.0)
   call mpi_bcast(refcrd_file, refcrd_len, MPI_CHARACTER, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

   if (ncolvars.eq.0) &
      return

#ifdef MPI
   if (mytaskid.ne.0) then
      allocate(cv(ncolvars), fcv_curr(ncolvars), stat = error)
      if (error.ne.0) &
         NFE_OUT_OF_MEMORY
   end if
#endif /* MPI */

   do n = 1, ncolvars
      call colvar_bootstrap(cv(n)%parent, n, amass)
      cv(n)%inst = colvar_value(cv(n)%parent, acrds)
      cv(n)%curr = ZERO
      cv(n)%harm = ZERO
      fcv_curr(n) = ZERO
   end do

   smd_mdstep = 0
   work = ZERO

#ifdef MPI
   if (mytaskid.ne.0) &
      return
#endif /* MPI */

   n = 1
   rewind(CV_UNIT)
   
   do while (n.le.ncolvars)
         path(:) = 0.0
         harm(:) = 0.0
         npath = 0
         nharm = 0
         path_mode = 'SPLINE'
         harm_mode = 'SPLINE'
         
         read(CV_UNIT,nml=colvar,err=667)
         
         ! path 

         if (npath.lt.2) then
            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(n)//',a/)') &
              NFE_ERROR, '''path'' for CV #', n, ' contains less than 2 values'
            call terminate()
         end if

         allocate(cv(n)%cvar_path(npath), stat = error)
         if (error.ne.0) &
            NFE_OUT_OF_MEMORY

         i = 1
         do while (i.le.npath)
            cv(n)%cvar_path(i) = path(i)
            i = i + 1
         end do

         ! harm

         if (nharm.lt.1) then
            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(n)//',a/)') &
               NFE_ERROR, '''harm'' for CV #', n, ' contains less than 1 value'
            call terminate()
         end if

         allocate(cv(n)%harm_path(nharm), stat = error)
         if (error.ne.0) &
            NFE_OUT_OF_MEMORY

         i = 1
         do while (i.le.nharm)
            cv(n)%harm_path(i) = harm(i)
            if (cv(n)%harm_path(i).lt.ZERO) then
               write (unit = ERR_UNIT, &
                 fmt = '(/a,a,'//pfmt(n)//',a,'//pfmt(i)//',a/)') NFE_ERROR, &
                  'CV #', n, ' : harm(', i, ') is negative'
               call terminate()
            end if
           i = i + 1
         end do

         call upper(harm_mode)
         ! harm_mode
         modestr = harm_mode
         cv(n)%harm_mode = mode_from_string(modestr)
         if (cv(n)%harm_mode.lt.0) then
               write (unit = ERR_UNIT, &
                 fmt = '(/a,a,'//pfmt(n)//',a,a,a/)') NFE_ERROR, &
                  'CV #', n, ' : unknown harm_mode = ''', trim(modestr), ''''
               call terminate()
         end if

         call upper(path_mode)
         ! path_mode
         modestr = path_mode
         cv(n)%path_mode = mode_from_string(modestr)
         if (cv(n)%path_mode.lt.0) then
               write (unit = ERR_UNIT, &
                 fmt = '(/a,a,'//pfmt(n)//',a,a,a/)') NFE_ERROR, &
                  'CV #', n, ' : unknown path_mode = ''', trim(modestr), ''''
               call terminate()
         end if

         nfe_atm_cnt = nfe_atm_cnt + cv_ni
         do i=1,cv_ni
            call AddToList(nfe_atm_lst,cv_i(i))
         end do

         n = n + 1
   end do
   close (CV_UNIT)

   ! done with parsing

   open (unit = smd_UNIT, file = output_file, iostat = error, &
         form = 'FORMATTED', action = 'WRITE', status = 'REPLACE')

   if (error.ne.0) then
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
         NFE_ERROR, 'failed to open ''', trim(output_file), ''' for writing'
      call terminate()
   end if

   write (unit = smd_UNIT, fmt = '(a,/a,/a)') '#', &
      '# MD time (ps), CV, handle_position, spring_constant, work', '#'

   call flush_UNIT(smd_UNIT,output_file)

   ! print summary & we'r done

   write (unit = OUT_UNIT, fmt = '(a,a)') NFE_INFO, &
      ' *  *  *  *  *  *  *  *  S T E E R E D  M.D.  *  *  *  *  *  *  *  *  *'

   write (unit = OUT_UNIT, fmt = '(a,/a,a,a)') NFE_INFO, NFE_INFO, &
      'output_file = ', trim(output_file)
   write (unit = OUT_UNIT, &
      fmt = '(a,a,'//pfmt(output_freq)//',a,'//pfmt(output_freq*pmemd_timestep(), 4)//',a)') &
      NFE_INFO, &
      'output_freq = ', output_freq, ' (', &
      output_freq*pmemd_timestep(), ' ps)'

   write (unit = OUT_UNIT, fmt = '(a)') NFE_INFO
   do n = 1, ncolvars
      write (unit = OUT_UNIT, fmt = '(a,a,'//pfmt(n)//',/a)') &
         NFE_INFO, 'CV #', n, NFE_INFO
      write (unit = OUT_UNIT, fmt = '(a,a)', advance = 'NO') &
         NFE_INFO, ' <> path = ('
      do i = 1, size(cv(n)%cvar_path)
         write (unit = OUT_UNIT, &
            fmt = '('//pfmt(cv(n)%cvar_path(i), 4)//')', advance = 'NO') &
            cv(n)%cvar_path(i)
         if (i.eq.size(cv(n)%cvar_path)) then
            write (unit = OUT_UNIT, fmt = '(a)') ')'
         else if (mod(i + 1, 5).eq.0) then
            write (unit = OUT_UNIT, fmt = '(a,/a,a)', advance = 'NO') &
               ',', NFE_INFO, '     '
         else
            write (unit = OUT_UNIT, fmt = '(a)', advance = 'NO') ', '
         end if
      end do
      write (unit = OUT_UNIT, fmt = '(a,a)', advance = 'NO') &
         NFE_INFO, ' <> path_mode = '
      call mode_write(cv(n)%path_mode, OUT_UNIT)

      write (unit = OUT_UNIT, fmt = '(a,a)', advance = 'NO') &
         NFE_INFO, ' <> harm = ('
      do i = 1, size(cv(n)%harm_path)
         write (unit = OUT_UNIT, &
            fmt = '('//pfmt(cv(n)%harm_path(i), 4)//')', advance = 'NO') &
            cv(n)%harm_path(i)
         if (i.eq.size(cv(n)%harm_path)) then
            write (unit = OUT_UNIT, fmt = '(a)') ')'
         else if (mod(i + 1, 5).eq.0) then
            write (unit = OUT_UNIT, fmt = '(a,/a,a)', advance = 'NO') &
               ',', NFE_INFO, '     '
         else
            write (unit = OUT_UNIT, fmt = '(a)', advance = 'NO') ', '
         end if
      end do
      write (unit = OUT_UNIT, fmt = '(a,a)', advance = 'NO') &
         NFE_INFO, ' <> harm_mode = '
      call mode_write(cv(n)%harm_mode, OUT_UNIT)
      if (colvar_is_quaternion(cv(n)%parent)) then
      !if (cv(n)%parent%q_index /= 0) then
          write (unit = OUT_UNIT, fmt = '(a,a,I3)') NFE_INFO, &
          ' <> q_index = ', cv(n)%parent%q_index
      end if
      if (colvar_has_axis(cv(n)%parent)) then
          write (unit = OUT_UNIT, fmt = '(a,a,f8.4,a,f8.4,a,f8.4,a)') NFE_INFO, &
           ' <> axis = [',cv(n)%parent%axis(1),', ', cv(n)%parent%axis(2),', ', cv(n)%parent%axis(3), ']'
      end if

      write (unit = OUT_UNIT, fmt = '(a)') NFE_INFO
      call colvar_print(cv(n)%parent, OUT_UNIT)
      write (unit = OUT_UNIT, fmt = '(a)') NFE_INFO
   end do

   write (unit = OUT_UNIT, fmt = '(a,a/)') NFE_INFO, &
      ' *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *'
   return
666 write(unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR,'Cannot read &smd namelist!'
    call terminate()

667 write(unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR,'Cannot read &colvar namelist!'
    call terminate()
end subroutine on_pmemd_init

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine on_pmemd_exit()

   use nfe_colvar_mod, only : colvar_cleanup
   use parallel_dat_mod

   implicit none

   integer :: n

   if (ncolvars.gt.0) then
      do n = 1, ncolvars
         call colvar_cleanup(cv(n)%parent)
         NFE_MASTER_ONLY_BEGIN
         deallocate(cv(n)%cvar_path)
         NFE_MASTER_ONLY_END
      end do

      deallocate(cv, fcv_curr)

      NFE_MASTER_ONLY_BEGIN
      write (unit = smd_UNIT, fmt = '(a/,a,f16.10/,a)') &
         '#', '# <> total work done: ', work, '#'
      close (smd_UNIT)
      NFE_MASTER_ONLY_END
   end if

   ncolvars = 0

end subroutine on_pmemd_exit

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Modified by M Moradi
! for driven ABMD
subroutine on_force(x, f, wdriven, udriven, pot)

   use nfe_colvar_mod
   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   double precision, intent(in) :: x(3,*)

   double precision, intent(inout) :: f(3,*)
   double precision, intent(out) :: wdriven
   double precision, intent(out) :: udriven
   double precision, intent(inout) :: pot

! Moradi end

   integer :: n, m

   double precision :: position, dwork, f2, dcv, dharm 
   double precision :: norm4_path(100), cv_N_path(4, 100), norm4(100), cv_N(4, 100)
   integer, DIMENSION(4) :: cv_q = (/COLVAR_QUATERNION0, COLVAR_QUATERNION1, &
                                   COLVAR_QUATERNION2, COLVAR_QUATERNION3/)


#  ifdef MPI
   integer :: error
#  endif /* MPI */

   if (ncolvars.eq.0) &
      return

   nfe_assert(multipmemd_rem().eq.0)

   position = dble(smd_mdstep)/dble(pmemd_nstlim())

   do n = 1, ncolvars
      cv(n)%inst = colvar_value(cv(n)%parent, x)
   end do
! Modified by M Moradi
! for driven ABMD
   NFE_MASTER_ONLY_BEGIN
   udriven = ZERO
   dwork = ZERO
!  cv_N to normalize quaternions
   do n = 1, ncolvars
      if (colvar_is_quaternion(cv(n)%parent)) then
        do m = 1, 4
         if (cv(n)%parent%type == cv_q(m)) then
            cv_N_path(m, cv(n)%parent%q_index) = mode_eval(cv(n)%path_mode, cv(n)%cvar_path, position)
         end if
        end do
      end if
   end do
   do n = 1, ncolvars
      if (colvar_is_quaternion(cv(n)%parent)) then
        norm4_path(cv(n)%parent%q_index) = sqrt(cv_N_path(1,cv(n)%parent%q_index)**2 &
                                         + cv_N_path(2,cv(n)%parent%q_index)**2 &
                                         + cv_N_path(3,cv(n)%parent%q_index)**2 &
                                         + cv_N_path(4,cv(n)%parent%q_index)**2)
      end if
   end do
! normalize quaternions (inst) 
   do n = 1, ncolvars
      if (colvar_is_quaternion(cv(n)%parent)) then
        do m = 1, 4
         if (cv(n)%parent%type == cv_q(m)) then
            cv_N(m, cv(n)%parent%q_index) = cv(n)%inst
         end if
        end do
      end if
   end do
   do n = 1, ncolvars
      if (colvar_is_quaternion(cv(n)%parent)) then
        norm4(cv(n)%parent%q_index) = sqrt(cv_N(1,cv(n)%parent%q_index)**2 &
                                         + cv_N(2,cv(n)%parent%q_index)**2 &
                                         + cv_N(3,cv(n)%parent%q_index)**2 &
                                         + cv_N(4,cv(n)%parent%q_index)**2)
      end if
   end do

   do n = 1, ncolvars
      dharm = mode_eval(cv(n)%harm_mode, cv(n)%harm_path, position) - &
            cv(n)%harm
      if (nfe_real_mdstep) &
            cv(n)%harm = mode_eval(cv(n)%harm_mode, cv(n)%harm_path, position)
      dcv = mode_eval(cv(n)%path_mode, cv(n)%cvar_path, position) - &
            cv(n)%curr
      if (nfe_real_mdstep) &
            cv(n)%curr = mode_eval(cv(n)%path_mode, cv(n)%cvar_path, position)
      if (colvar_is_quaternion(cv(n)%parent)) then
          cv(n)%curr = cv(n)%curr / norm4_path(cv(n)%parent%q_index)
      else
          cv(n)%curr = cv(n)%curr
      end if
      ! Fix the quaternions themself	
      if (colvar_is_quaternion(cv(n)%parent)) then
          cv(n)%inst = cv(n)%inst / norm4(cv(n)%parent%q_index)
      else
          cv(n)%inst = cv(n)%inst
      end if

      f2 = fcv_curr(n)
      fcv_curr(n) = cv(n)%harm*colvar_difference(cv(n)%parent,cv(n)%curr,cv(n)%inst)
      pot = pot + cv(n)%harm*colvar_difference(cv(n)%parent,cv(n)%curr,cv(n)%inst)**2/2
      f2 = (f2 + fcv_curr(n))/2
      dwork = dwork + f2*dcv + dharm*&
              colvar_difference(cv(n)%parent,cv(n)%curr,cv(n)%inst)**2/2
      udriven = udriven + cv(n)%harm*colvar_difference(cv(n)%parent,cv(n)%curr,cv(n)%inst)**2/2
   end do
   wdriven = work + dwork
   NFE_MASTER_ONLY_END
! Moradi end

#  ifdef MPI
   call mpi_bcast(fcv_curr, ncolvars, &
      MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
   call mpi_bcast(pot, 1, &
      MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#  endif /* MPI */

   ! FIXME: virial
   do n = 1, ncolvars
      call colvar_force(cv(n)%parent, x, fcv_curr(n), f)
   end do

   NFE_MASTER_ONLY_BEGIN
   if (nfe_real_mdstep) then

      if (smd_mdstep.gt.0) &
         work = work + dwork

      if (mod(smd_mdstep, output_freq).eq.0) then
         write (unit = smd_UNIT, fmt = '(f12.4,1x)', advance = 'NO') &
            pmemd_mdtime()
         do n = 1, ncolvars
            write (unit = smd_UNIT, fmt = '(f16.8,1x)', advance = 'NO') &
               cv(n)%inst
         end do
         do n = 1, ncolvars
            write (unit = smd_UNIT, fmt = '(f16.8,1x)', advance = 'NO') &
               cv(n)%curr
         end do
         do n = 1, ncolvars
            write (unit = smd_UNIT, fmt = '(f16.8,1x)', advance = 'NO') &
               cv(n)%harm
         end do
         write (unit = smd_UNIT, fmt = '(f16.8)') work
         call flush_UNIT(smd_UNIT, output_file)
      end if

      smd_mdstep = smd_mdstep + 1

   end if ! nfe_real_mdstep
   NFE_MASTER_ONLY_END

end subroutine on_force

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function mode_eval(mode, path, t) result(y)

   NFE_USE_AFAILED

   implicit none

   double precision :: y
   integer, intent(in) :: mode
   double precision, intent(in) :: path(:), t

   if (mode.eq.MODE_LINES) then
      y = lines(path, t)
   else if (mode.eq.MODE_SPLINE) then
      y = spline(path, t)
   else
      nfe_assert_not_reached()
      y = 77.0D0
   end if

end function mode_eval

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function lines(path, t) result(y)

   use nfe_lib_mod
 
   implicit none

   double precision :: y
   double precision, intent(in) :: path(:), t

   integer :: npoints, n

   double precision :: s

   npoints = size(path)

   if (npoints.eq.1) then
      y = path(1)
      return
   end if

   if (t.le.ZERO) then
      y = path(1)
   else if (t.ge.ONE) then
      y = path(npoints)
   else
      n = 1 + int(floor((npoints - 1)*t))
      nfe_assert(n.ge.1)
      nfe_assert(n.lt.npoints)

      ! 0 < s < 1 between path(n) and path(n + 1)
      s = (npoints - 1)*t - n + 1
      y = (ONE - s)*path(n) + s*path(n + 1)
   end if

end function lines

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function spline(path, t) result(y)

   use nfe_lib_mod

   implicit none

   double precision :: y
   double precision, intent(in) :: path(:), t

   integer :: npoints, n

   double precision :: m1, m2
   double precision :: s, s2, s3

   npoints = size(path)

   if (npoints.eq.1) then
      y = path(1)
      return
   end if

   if (t.le.ZERO) then
      y = path(1)
   else if (t.ge.ONE) then
      y = path(npoints)
   else
      n = 1 + int(floor((npoints - 1)*t))
      nfe_assert(n.ge.1)
      nfe_assert(n.lt.npoints)

      if (npoints.eq.2) then
         m1 = ZERO
         m2 = ZERO
      else if (n.eq.1) then
         m1 = ZERO
         m2 = (path(3) - path(1))/2
      else if ((n + 1).eq.npoints) then
         m1 = (path(n + 1) - path(n - 1))/2
         m2 = ZERO
      else
         m1 = (path(n + 1) - path(n - 1))/2
         m2 = (path(n + 2) - path(n))/2
      end if

      ! compute the value
      s = (npoints - 1)*t - n + 1

      s2 = s*s
      s3 = 2*s - dble(3)

      y = (s2*s3 + ONE)*path(n) &
        + s*(s*(s - TWO) + ONE)*m1 &
        - s2*s3*path(n + 1) &
        + s2*(s - ONE)*m2
   end if

end function spline

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine mode_write(m, lun)

   implicit none

   integer, intent(in) :: m
   integer, intent(in) :: lun

   if (m.eq.MODE_LINES) then
      write (unit = lun, fmt = '(a)') 'LINES'
   else if (m.eq.MODE_SPLINE) then
      write (unit = lun, fmt = '(a)') 'SPLINE'
   else
      write (unit = lun, fmt = '(a)') 'UNKNOWN'
   end if

end subroutine mode_write

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function mode_from_string(str) result(mode)

    implicit none

    integer :: mode
    character(*), intent(in) :: str

    if (str.eq.'LINES') then
        mode = MODE_LINES
    else if (str.eq.'SPLINE') then
        mode = MODE_SPLINE
    else
        mode = -1
    end if

end function mode_from_string

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module nfe_smd_mod
