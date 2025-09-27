#include "copyright.i"
!*******************************************************************************
!
! Module: get_cmdline_mod
!
! Description: Module controlling getting command line arguments from either
!              a line in a group file or from the command-line itself
!
!*******************************************************************************

module get_cmdline_mod

  implicit none

  private       add_suffix
  logical, save, public :: cpein_specified = .false.

#ifdef MPI
  public get_cmdline_bcast
#endif /* MPI */

contains

#ifdef USE_PXFGETARG
!*******************************************************************************
!
! Subroutine:  pmemd_getarg
!
! Description: wrapper for intrinsic getarg to make it work with nonstandard
!              systems
!
!*******************************************************************************

subroutine pmemd_getarg(i, c)

  implicit none

! Formal arguments:

  character(*)  :: c
  integer       :: i

! Local variables:

  integer       :: ilen
  integer       :: ierror

! Note that the pxfgetarg subroutine is broken in ifc 7.1; use getarg instead
! for any compiler that supports it!

  call pxfgetarg(i, c, ilen, ierror)

  if (ierror .ne. 0) then
    write(mdout, '(a)') 'ERROR: Command line parsing failed in getarg!'
    call mexit(6, 1)
  end if

  return

end subroutine pmemd_getarg
#else
! The preferred command line subroutine is getarg.
#define pmemd_getarg    getarg
#endif /* USE_PXFGETARG */

!*******************************************************************************
!
! Subroutine:  get_cmdline
!
! Description: Gets command-line arguments either from CL or from groupfile
!
! OUTPUT: File names and command-line arguments
!
! Author: George Seibel. Modifications by JMS to handle strings for multipmemd
!
!*******************************************************************************

subroutine get_cmdline(terminal_flag)

  use file_io_dat_mod
  use gbl_constants_mod, only : error_hdr, VERSION_STRING
  use pmemd_lib_mod
  use parallel_dat_mod
#ifdef MPI
  use remd_mod,         only  : remd_method, rremd_type, remd_random_partner
#endif

  implicit none

! Passed variables:

  logical, intent(in out)  :: terminal_flag

! Local variables:

  character(max_fn_len) :: arg  ! temp buffer for each of the whitespace-
                        ! delimited command line words.
  integer       :: iarg ! arg pointer, final number of arguments.
  integer       :: indx

#ifdef MPI
  logical       :: in_groupfile ! if we're reading this from the groupfile
                                ! or _really_ on the command-line
#endif

! Flags that determine if/when the -suffix gets used:

  logical       :: mdout_specified
  logical       :: mdinfo_specified
  logical       :: mdcrd_specified
  logical       :: mdvel_specified
  logical       :: mdfrc_specified
  logical       :: mden_specified
  logical       :: restrt_specified
  logical       :: logfile_specified
  logical       :: cpin_specified
  logical       :: cpout_specified
  logical       :: cprestrt_specified
  logical       :: cein_specified
  logical       :: ceout_specified
  logical       :: cerestrt_specified
  logical       :: cpeout_specified
  logical       :: cperestrt_specified
  logical       :: outfile_suffix_specified
  logical       :: amdlog_specified
  logical       :: gamdlog_specified
  !!for gamd-remd
  logical       :: gamdres_name_specified
  !! end
  logical       :: scaledMDlog_specified
  logical       :: always_apply_suffix
  
  integer       :: suffix_len

  mdout_specified = .false.
  mdinfo_specified = .false.
  mdcrd_specified = .false.
  mdvel_specified = .false.
  mdfrc_specified = .false.
  mden_specified = .false.
  restrt_specified = .false.
  logfile_specified = .false.
  cpin_specified = .false.
  cpout_specified = .false.
  cprestrt_specified = .false.
  cein_specified = .false.
  ceout_specified = .false.
  cerestrt_specified = .false.
  cpein_specified = .false.
  cpeout_specified = .false.
  cperestrt_specified = .false.
  outfile_suffix_specified = .false.
  amdlog_specified = .false.
  gamdlog_specified = .false.
  gamdres_name_specified = .false.
  scaledMDlog_specified = .false.
  always_apply_suffix = .false.

! Default file names:

  mdin_name   = 'mdin'
  mdout_name  = 'mdout'
  inpcrd_name = 'inpcrd'
  prmtop_name = 'prmtop'
  restrt_name = 'restrt'
  refc_name   = 'refc'
  mdvel_name  = 'mdvel'
  mdfrc_name  = 'mdfrc'
  mden_name   = 'mden'
  mdcrd_name  = 'mdcrd'
  mdinfo_name = 'mdinfo'
  cpin_name   = 'cpin'
  cpout_name  = 'cpout'
  cprestrt_name = 'cprestrt'
  cein_name   = 'cein'
  ceout_name  = 'ceout'
  cerestrt_name = 'cerestrt'
  cpein_name   = 'cpein'
  cpeout_name  = 'cpeout'
  cperestrt_name = 'cperestrt'
  logfile_name = 'logfile'
! PHMD
  phmdin_name = 'phmdin'
  phmdout_name = 'phmdout'
  phmdstrt_name = 'phmdstrt'
  phmdparm_name = 'phmdparm'
  phmdrestrt_name = 'phmdrestrt'
#ifdef TIDECOMP
  decomp_name = 'decomp.log'
#endif
#ifdef MPI
  proc_map_name = ' '

! To make sure we have unique output file names in all independent replicas,
! activate the outfile_suffix_specified, and build the suffix from the replica
! number (i.e. .001, .002, etc.). Activate the outfile_suffix_specified here,
! and if -suffix is provided, take that as a key to apply it to all files, even
! the ones we DID specify in the groupfile.

  if ( numgroups .gt. 1) then
    outfile_suffix_specified = .true.
    call assign_group_suffix(outfile_suffix)
  end if

  in_groupfile = len_trim(groupline_buffer) .ne. 0
#endif /* MPI */


! Default status of output: new

  owrite = 'N'

! Get com line arguments:

  iarg = 0
  indx = iargc_wrap()

  if (indx .eq. 0) goto 20

  10 continue

  iarg = iarg + 1

  call getarg_wrap(iarg, arg)

  ! Look for a --version or -V flag and dump out the version of pmemd
  if (arg .eq. '-V' .or. arg .eq. '--version') then
#ifdef MPI
    if (in_groupfile) then
      write(0, '(2a)') trim(arg), &
                     ' Bad flag in groupfile. Must be on command-line!'
      call mexit(mdout, 0)
    end if
#endif
    call getarg_wrap(0, arg)
    call basename(arg)
    write(6, '(3a)') trim(arg), ': ', trim(VERSION_STRING)
    terminal_flag = .true.
    return ! no need to do anything more here.
  else if (arg .eq. '-O') then
    owrite = 'U'
  else if (arg .eq. '-i') then
    iarg = iarg + 1
    call getarg_wrap(iarg, mdin_name)
  else if (arg .eq. '-o') then
    iarg = iarg + 1
    call getarg_wrap(iarg, mdout_name)
    mdout_specified = .true.
  else if (arg .eq. '-p') then
    iarg = iarg + 1
    call getarg_wrap(iarg, prmtop_name)
  else if (arg .eq. '-c') then
    iarg = iarg + 1
    call getarg_wrap(iarg, inpcrd_name)
  else if (arg .eq. '-r') then
    iarg = iarg + 1
    call getarg_wrap(iarg, restrt_name)
    restrt_specified = .true.
  else if (arg .eq. '-ref' .or. arg .eq.'-z') then
    iarg = iarg + 1
    call getarg_wrap(iarg, refc_name)
  else if (arg .eq. '-e') then
    iarg = iarg + 1
    call getarg_wrap(iarg, mden_name)
    mden_specified = .true.
  else if (arg .eq. '-v') then
    iarg = iarg + 1
    call getarg_wrap(iarg, mdvel_name)
    mdvel_specified = .true.
  else if (arg .eq. '-frc') then
    iarg = iarg + 1
    call getarg_wrap(iarg, mdfrc_name)
    mdfrc_specified = .true.
  else if (arg .eq. '-x' .or. arg .eq.'-t') then
    iarg = iarg + 1
    call getarg_wrap(iarg, mdcrd_name)
    mdcrd_specified = .true.
  else if (arg .eq. '-inf') then
    iarg = iarg + 1
    call getarg_wrap(iarg, mdinfo_name)
    mdinfo_specified = .true.
  else if (arg .eq. '-l') then
    iarg = iarg + 1
    call getarg_wrap(iarg, logfile_name)
    logfile_specified = .true.
  else if (arg .eq. '-suffix') then
    iarg = iarg + 1
    call getarg_wrap(iarg, outfile_suffix)
    if (numgroups .gt. 1) always_apply_suffix = .true.
    outfile_suffix_specified = .true.
  else if (arg .eq. '-AllowSmallBox') then
    smbx_permit = 1
  ! PHMD
  else if (arg .eq. '-phmdin') then
    iarg = iarg + 1
    call getarg_wrap(iarg, phmdin_name)
  else if (arg .eq. '-phmdout') then
    iarg = iarg + 1
    call getarg_wrap(iarg, phmdout_name)
  else if (arg .eq. '-phmdrestrt') then
    iarg = iarg + 1
    call getarg_wrap(iarg, phmdrestrt_name)
  else if (arg .eq. '-phmdstrt') then
    iarg = iarg + 1
    call getarg_wrap(iarg, phmdstrt_name)
  else if (arg .eq. '-phmdparm') then
    iarg = iarg + 1
    call getarg_wrap(iarg, phmdparm_name)
#ifdef TIDECOMP
  else if (arg .eq. '-decomp') then
    iarg = iarg + 1
    call getarg_wrap(iarg, decomp_name)
#endif
  else if (arg .eq. '-cpin') then
    iarg = iarg + 1
    call getarg_wrap(iarg, cpin_name)
    cpin_specified = .true.
  else if (arg .eq. '-cpout') then
    iarg = iarg + 1
    call getarg_wrap(iarg, cpout_name)
    cpout_specified = .true.
  else if (arg .eq. '-cprestrt') then
    iarg = iarg + 1
    call getarg_wrap(iarg, cprestrt_name)
    cprestrt_specified = .true.
  else if (arg .eq. '-cein') then
    iarg = iarg + 1
    cein_specified = .true.
    call getarg_wrap(iarg, cein_name)
  else if (arg .eq. '-ceout') then
    iarg = iarg + 1
    call getarg_wrap(iarg, ceout_name)
    ceout_specified = .true.
  else if (arg .eq. '-cerestrt') then
    iarg = iarg + 1
    call getarg_wrap(iarg, cerestrt_name)
    cerestrt_specified = .true.
  else if (arg .eq. '-cpein') then
    iarg = iarg + 1
    call getarg_wrap(iarg, cpein_name)
    cpein_specified = .true.
  else if (arg .eq. '-cpeout') then
    iarg = iarg + 1
    call getarg_wrap(iarg, cpeout_name)
    cpeout_specified = .true.
  else if (arg .eq. '-cperestrt') then
    iarg = iarg + 1
    call getarg_wrap(iarg, cperestrt_name)
    cperestrt_specified = .true.
  else if (arg .eq. '-amd') then
    iarg = iarg + 1
    call getarg_wrap(iarg, amdlog_name)
    amdlog_specified = .true.
  else if (arg .eq. '-gamd') then
    iarg = iarg + 1
    call getarg_wrap(iarg, gamdlog_name)
    gamdlog_specified = .true.
  else if (arg .eq. '-scaledMD') then
    iarg = iarg + 1
    call getarg_wrap(iarg, scaledMDlog_name)
    scaledMDlog_specified = .true.
!! SAMS  
  else if (arg .eq. '-sams_init') then
    iarg = iarg + 1
    call getarg_wrap(iarg, sams_init_name)
  else if (arg .eq. '-sams_log') then
    iarg = iarg + 1
    call getarg_wrap(iarg, sams_log_name)
  else if (arg .eq. '-sams_rest') then
    iarg = iarg + 1
    call getarg_wrap(iarg, sams_restart_name)
!!lambda-sceduling
  else if (arg .eq. '-lambda_sch') then
    iarg = iarg + 1
    call getarg_wrap(iarg, lambda_scheduling_name)
!!tau-sceduling
  else if (arg .eq. '-tau_sch') then
    iarg = iarg + 1
    call getarg_wrap(iarg, tau_scheduling_name)
  else if (arg .eq. '-reservoir') then
    iarg = iarg + 1
    call getarg_wrap(iarg, reservoir_name(1))
  else if (arg .eq. '-reservoir2') then
    iarg = iarg + 1
    call getarg_wrap(iarg, reservoir_name(2))  

  else if (arg .eq. '-help' .or. arg .eq. '--help' .or. arg .eq. '-H') then
#ifdef MPI
    if (in_groupfile) then
      write(0, '(2a)') trim(arg), &
                     ' Bad flag in groupfile. Must be on command-line!'
      call mexit(mdout, 0)
    end if
#endif
    write(mdout, 9000)
    terminal_flag = .true.
    return

#ifdef MPI
  else if (arg .eq. '-ng') then
    iarg = iarg + 1
    call getarg_wrap(iarg, arg)
    read(arg, '(1i5)') numgroups
  else if (arg .eq. '-rem') then
    ! Make sure we didn't call
    iarg = iarg + 1
    call getarg_wrap(iarg, arg)
    read(arg, '(1i5)', err=666) remd_method
    if (remd_method .ne. -1 .and. len_trim(remd_dimension_name) .ne. 0) goto 667
  else if (arg .eq. '-remd-file') then
    iarg = iarg + 1
    call getarg_wrap(iarg, remd_dimension_name)
    if (remd_method .ne. 0) goto 667
    ! Make remd_method -1 here, since that's multi-D REMD indicator
    remd_method = -1
  else if (arg .eq. '-remlog') then
    iarg = iarg + 1
    call getarg_wrap(iarg, remlog_name)
  else if (arg .eq. '-remtype') then
    iarg = iarg + 1
    call getarg_wrap(iarg, remtype_name)
  else if (arg .eq. '-remrandompartner') then
    iarg = iarg + 1
    call getarg_wrap(iarg, arg)
    read(arg, '(1i5)', err=668) remd_random_partner

  else if (arg .eq. '-rremd') then
    iarg = iarg + 1
    call getarg_wrap(iarg, arg)
    read(arg, '(1i5)', err=666) rremd_type
  else if (arg .eq. '-hybridsolvent_remd_traj') then
    iarg = iarg + 1
    call getarg_wrap(iarg, hybridsolvent_remd_traj_name)
  else if (arg .eq. '-groupfile') then
    iarg = iarg + 1
    call getarg_wrap(iarg, groupfile_name)
  else if (arg .eq. '-ng-nonsequential') then
    ng_nonsequential = .true.
  else if (arg .eq. '-gpes') then
    iarg = iarg + 1
    call getarg_wrap(iarg, proc_map_name)
  else if (arg .eq. '-p4pg') then
    iarg = iarg + 1
  else if (arg .eq. '-p4wd') then
    iarg = iarg + 1
  else if (arg .eq. '-np') then
    iarg = iarg + 1
  else if (arg .eq. '-mpedbg') then
    continue
  else if (arg .eq. '-dbx') then
    continue
  else if (arg .eq. '-gdb') then
    continue
#endif
  else
    if (arg .eq. ' ') goto 20
    write(mdout, '(/,2x,a,a)') 'unknown flag: ', arg
    write(mdout, 9000)
    call mexit(6, 1)
  end if

  if (cpein_specified) then
    if (cpin_specified) then
      write(mdout, '(a)') 'ERROR: The cpein file cannot be specified together with the cpin file!'
      call mexit(6, 1)
    end if
    if (cein_specified) then
      write(mdout, '(a)') 'ERROR: The cpein file cannot be specified together with the cein file!'
      call mexit(6, 1)
    end if
  end if

  if (iarg .lt. indx) goto 10

  20 continue

  if (outfile_suffix_specified) then

    suffix_len = len_trim(outfile_suffix)

    if (.not. mdout_specified .or. always_apply_suffix) &
      call add_suffix(mdout_name, outfile_suffix, suffix_len)

    if (.not. mdinfo_specified .or. always_apply_suffix) &
      call add_suffix(mdinfo_name, outfile_suffix, suffix_len)

    if (.not. mdcrd_specified .or. always_apply_suffix) &
      call add_suffix(mdcrd_name, outfile_suffix, suffix_len)

    if (.not. mdvel_specified .or. always_apply_suffix) &
      call add_suffix(mdvel_name, outfile_suffix, suffix_len)

    if (.not. mdfrc_specified .or. always_apply_suffix) &
      call add_suffix(mdfrc_name, outfile_suffix, suffix_len)

    if (.not. mden_specified .or. always_apply_suffix) &
      call add_suffix(mden_name, outfile_suffix, suffix_len)

    if (.not. restrt_specified .or. always_apply_suffix) &
      call add_suffix(restrt_name, outfile_suffix, suffix_len)

    if (.not. logfile_specified .or. always_apply_suffix) &
      call add_suffix(logfile_name, outfile_suffix, suffix_len)

    if (.not. cpout_specified .or. always_apply_suffix) &
      call add_suffix(cpout_name, outfile_suffix, suffix_len)

    if (.not. cprestrt_specified .or. always_apply_suffix) &
      call add_suffix(cprestrt_name, outfile_suffix, suffix_len)

    if (.not. ceout_specified .or. always_apply_suffix) &
      call add_suffix(ceout_name, outfile_suffix, suffix_len)

    if (.not. cerestrt_specified .or. always_apply_suffix) &
      call add_suffix(cerestrt_name, outfile_suffix, suffix_len)

    if (.not. cpeout_specified .or. always_apply_suffix) &
      call add_suffix(cpeout_name, outfile_suffix, suffix_len)

    if (.not. cperestrt_specified .or. always_apply_suffix) &
      call add_suffix(cperestrt_name, outfile_suffix, suffix_len)

    if (.not. amdlog_specified .or. always_apply_suffix) &
      call add_suffix(amdlog_name, outfile_suffix, suffix_len)

    if (.not. gamdlog_specified .or. always_apply_suffix) &
      call add_suffix(gamdlog_name, outfile_suffix, suffix_len)
    if (.not. gamdres_name_specified .or. always_apply_suffix) &
      call add_suffix(gamdres_name, outfile_suffix, suffix_len)

    if (.not. scaledMDlog_specified .or. always_apply_suffix) &
      call add_suffix(scaledMDlog_name, outfile_suffix, suffix_len)
  end if

  return
#ifdef MPI
9000 format(/, 2x, &
      'usage: pmemd  [-O] -i mdin -o mdout -p prmtop -c inpcrd -r restrt', &
      /14x, '[-ref refc -x mdcrd -v mdvel -frc mdfrc -e mden &
             &-inf mdinfo -l logfile]', &
      /14x, '[-ng numgroups -groupfile groupfile -rem remd_method]', &
      /14x, '[-rremd rremd_type -reservoir reservoir.nc]', &
      /14x, '[-hybridsolvent_remd_traj hybridsolvent_remd_traj_name]', &
      /14x, '[-amd amdlog_name -gamd gamdlog_name -scaledMD scaledMDlog_name &
             &-suffix output_files_suffix]', &
      /14x, '[-sf sfout_name]'/)

666 write(mdout, '(2a)') error_hdr, 'bad integer read on command-line (-rem)'
    call mexit(mdout, 1)

667 write(mdout, '(2a)') error_hdr, &
      'do not specify -rem # when you supply a -remdim <dimension file>'
    call mexit(mdout, 1)

668 write(mdout, '(2a)') error_hdr, 'bad integer read on command-line (-remrandompartner)'
    call mexit(mdout, 1)

#else
9000 format(/, 2x, &
      'usage: pmemd  [-O] -i mdin -o mdout -p prmtop -c inpcrd -r restrt', &
      /14x, '[-ref refc -x mdcrd -v mdvel -frc mdfrc -e mden &
             &-inf mdinfo -l logfile]', &
      /14x, '[-amd amdlog_name -gamd gamdlog_name -scaledMD scaledMDlog_name &
             &-suffix output_files_suffix]',/)
#endif /* MPI */

end subroutine get_cmdline

!*******************************************************************************
!
! Subroutine:  add_suffix
!
! Description: adds a suffix to each file
!
!*******************************************************************************

subroutine add_suffix(file_name, suffix, suffix_len)

  use file_io_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:

  character(max_fn_len), intent(in out) :: file_name
  character(max_fn_len), intent(in)     :: suffix
  integer, intent(in)                   :: suffix_len

  if (suffix_len + len_trim(file_name) .lt. max_fn_len) then
    file_name = trim(file_name) // '.' // trim(suffix)
  else
    write(mdout, '(a)') 'ERROR: Filename with suffix too long!'
    call mexit(6, 1)
  end if

  return

end subroutine add_suffix

!*******************************************************************************
!
! Function: iargc_wrap
!
! Description: wrapper for iargc() intrinsic so we can read arguments from a
!              groupfile. Counts number of arguments in groupfile line or on
!              command-line
!
!*******************************************************************************

integer function iargc_wrap()

#ifdef MPI
  use file_io_dat_mod,  only : groupline_buffer
  use pmemd_lib_mod,    only : get_num_tokens
#endif

  implicit none

#ifdef USE_PXFGETARG
  integer  :: ipxfargc
#else
  integer  :: iargc  ! intrinsic iargc if we're on command-line
#endif

#ifdef MPI
  ! If numgroups > 1 and groupfile has been found, parse groupline

  if (len_trim(groupline_buffer) .ne. 0) then
    call get_num_tokens(groupline_buffer, iargc_wrap)
  else
# ifdef USE_PXFGETARG
    iargc_wrap = ipxfargc()
# else
    iargc_wrap = iargc()
# endif
  end if

#else

# ifdef USE_PXFGETARG
  iargc_wrap = ipxfargc()
# else
  iargc_wrap = iargc()
# endif /* USE_PXFGETARG */

#endif /* MPI */

end function iargc_wrap

!*******************************************************************************
!
! Subroutine: getarg_wrap
!
! Description: Wrapper subroutine for intrinsic getarg to support parsing
!              arguments from a groupfile
!
!*******************************************************************************

subroutine getarg_wrap(iarg, arg)

#ifdef MPI
  use file_io_dat_mod, only : groupline_buffer
  use pmemd_lib_mod, only : get_token
#endif

  implicit none

! Passed variables

  integer, intent(in)       :: iarg  ! argument number

  character(*), intent(out) :: arg   ! argument string

#ifdef MPI
  if (len_trim(groupline_buffer) .gt. 0) then
    call get_token(groupline_buffer, iarg, arg)
  else
    call pmemd_getarg(iarg, arg)
  end if
#else
  call pmemd_getarg(iarg, arg)
#endif /* MPI */

end subroutine getarg_wrap


#ifdef MPI
!*******************************************************************************
!
! Subroutine: get_cmdline_bcast
!
! Description: Broadcasts get_cmdline data to the whole MPI universe
!
!*******************************************************************************

subroutine get_cmdline_bcast()

  use parallel_dat_mod

  implicit none

  call mpi_bcast(cpein_specified, 1, mpi_logical, 0, pmemd_comm, err_code_mpi)

end subroutine get_cmdline_bcast
#endif /* MPI */

end module get_cmdline_mod
