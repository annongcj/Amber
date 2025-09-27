!*******************************************************************************
!
! Module: file_io_dat_mod
!
! Description: <TBS>
!
!*******************************************************************************

module file_io_dat_mod

  implicit none

! File-related data.  No need to broadcast this stuff, since only the master
! does i/o:

  integer, parameter            :: max_fn_len = 256
#ifdef MPI
  integer, parameter            :: grpfl_line_len = 4096
#endif

  character(max_fn_len), save   :: mdin_name
  character(max_fn_len), save   :: mdout_name
  character(max_fn_len), save   :: mdinfo_name
  character(max_fn_len), save   :: prmtop_name
  character(max_fn_len), save   :: inpcrd_name
  character(max_fn_len), save   :: refc_name
  character(max_fn_len), save   :: mdcrd_name
  character(max_fn_len), save   :: mdvel_name
  character(max_fn_len), save   :: mdfrc_name
  character(max_fn_len), save   :: mden_name
  character(max_fn_len), save   :: restrt_name
  character(max_fn_len), save   :: logfile_name
  character(max_fn_len), save   :: infile_suffix
  character(max_fn_len), save   :: outfile_suffix
  character(max_fn_len), save   :: cpin_name
  character(max_fn_len), save   :: cpout_name
  character(max_fn_len), save   :: cprestrt_name
  character(max_fn_len), save   :: cein_name
  character(max_fn_len), save   :: ceout_name
  character(max_fn_len), save   :: cerestrt_name
  character(max_fn_len), save   :: cpein_name
  character(max_fn_len), save   :: cpeout_name
  character(max_fn_len), save   :: cperestrt_name
!PHMD
  character(max_fn_len), save   :: phmdin_name
  character(max_fn_len), save   :: phmdout_name
  character(max_fn_len), save   :: phmdrestrt_name
  character(max_fn_len), save   :: phmdstrt_name
  character(max_fn_len), save   :: phmdparm_name
  character(max_fn_len), save   :: decomp_name

!AMD file for reweighting values
  character(max_fn_len), save   :: amdlog_name = 'amd.log'
!GaMD file for reweighting values
  character(max_fn_len), save   :: gamdlog_name = 'gamd.log'
  character(max_fn_len), save   :: gamdres_name = 'gamd-restart.dat'
  character(max_fn_len), save   :: gamdlig_name = 'gamd-ligB.dat'
!scaledMD file for reweighting values
  character(max_fn_len), save   :: scaledMDlog_name = 'scaledMD.log'

!SAMS
  character(max_fn_len), save   :: sams_restart_name = 'sams.rest'
  character(max_fn_len), save   :: sams_log_name = 'sams.log'
  character(max_fn_len), save   :: sams_init_name = 'sams.init'

! lambda_scheduling
  character(max_fn_len), save   :: lambda_scheduling_name = 'lambda.sch'
! tau_scheduling
  character(max_fn_len), save   :: tau_scheduling_name = 'tau.sch'

#ifdef MPI
! Set file names for multipmemd and REMD. proc_map refers
! to the file describing how many processors to assign to
! each communicator
  character(max_fn_len), save   :: groupfile_name = ' '
  character(max_fn_len), save   :: remlog_name = 'rem.log'
  character(max_fn_len), save   :: remtype_name = 'rem.type'
  !!character(max_fn_len), save   :: reservoir_name(2) = (/ 'reservoir.nc ',  'reservoir2.nc' /)
  character(max_fn_len), save   :: hybridsolvent_remd_traj_name = ' '
  character(max_fn_len), save   :: proc_map_name
  character(max_fn_len), save   :: remd_dimension_name = ' '
  character(grpfl_line_len), save :: groupline_buffer = ' '
  character(grpfl_line_len), save :: proc_map_buffer = ' '
#endif
  character(max_fn_len), save   :: reservoir_name(2) = (/ 'reservoir.nc ',  'reservoir2.nc' /)
  character, save       :: owrite

! Logical unit numbers for pmemd files.

  integer, parameter    :: mdin        = 115
  integer, parameter    :: mdout       = 116
  integer, parameter    :: mdinfo      =  7
  integer, parameter    :: prmtop      =  8
  integer, parameter    :: inpcrd      =  9
  integer, parameter    :: refc        = 10
  ! who was 11?
  integer, parameter    :: mdcrd       = 12
  integer, parameter    :: mdvel       = 13
  integer, parameter    :: mdfrc       = 14
  integer, parameter    :: mden        = 15
  integer, parameter    :: restrt      = 16
  integer, parameter    :: restrt2     = 17
  integer, parameter    :: logfile     = 18
#ifdef MPI
  integer, parameter    :: groupfile     = 19
  integer, parameter    :: remlog        = 20
  integer, parameter    :: proc_map      = 21
  integer, parameter    :: remd_file     = 22
#endif
  integer, parameter    :: cpin          = 23
  integer, parameter    :: cpout         = 24
  integer, parameter    :: cprestrt      = 25

#ifdef MPI
! Maybe not quite the right place to put this, but it can't be in
! multipmemd_mod, since we'll get a cyclic dependency.

  logical, save         :: ng_nonsequential = .false.
#endif
  integer, parameter    :: amdlog        = 26
  integer, parameter    :: scaledMDlog   = 27
  integer, parameter    :: gamdlog       = 28
  integer, parameter    :: gamdres       = 29
  integer, parameter    :: cein          = 30
  integer, parameter    :: ceout         = 31
  integer, parameter    :: cerestrt      = 32
!PHMD
  integer, parameter    :: phmd_unit     = 33
  integer, parameter    :: phmdout_unit  = 34
  integer, parameter    :: phmdstrt_unit = 35
  integer, parameter    :: phmdparm_unit = 36

  integer, parameter    :: cpein         = 37
  integer, parameter    :: cpeout        = 38
  integer, parameter    :: cperestrt     = 39
  integer, parameter    :: remtype       = 40
!SAMS
  integer, parameter    ::  sams_log      = 41
  integer, parameter    ::  sams_init     = 42
  integer, parameter    ::  sams_rest  = 43
!lambda-scheduling
  integer, parameter    ::  lambda_sch    = 44
  integer, parameter    :: sfout         = 45

! File for writing stripped coordinates during Hybrid Solvent REMD
  integer, parameter    :: hybridsolvent_remd_traj = 45
! Lig-GaMD  
  integer, parameter    :: gamdlig       = 47
  integer, parameter    :: decomp_unit   = 48
! Flags to allow certain risky simulations
  integer, save         :: smbx_permit = 0

end module file_io_dat_mod
