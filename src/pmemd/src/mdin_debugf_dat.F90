!*******************************************************************************
!
! Module: mdin_debugf_dat_mod
!
! Description: Module supporting the debugf namelist. Currently only
!              a very small subset of options is supported.
!              
!*******************************************************************************

module mdin_debugf_dat_mod

use file_io_dat_mod

  implicit none

! WARNING: Only a small subset of the debugf namelist is supported.
!          Current options supported =
!
!          do_charmm_dump_gold = 1  Turns on dumping of CHARMM force
!                                   field energies and forces.

  integer, parameter    :: mdin_debugf_int_cnt = 1

  integer                       do_charmm_dump_gold

  common / mdin_debugf_int /    do_charmm_dump_gold
                                

  save  :: / mdin_debugf_int /

  private       :: debugf

  namelist /debugf/     do_charmm_dump_gold

contains

!*******************************************************************************
!
! Subroutine:  init_mdin_debugf_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine init_mdin_debugf_dat()

  use file_io_mod

  implicit none

! Local variables:

  integer               :: ifind

  do_charmm_dump_gold = 0

  rewind(mdin)                  ! Insurance against maintenance mods.

  call nmlsrc('debugf', mdin, ifind)

  if (ifind .ne. 0) read(mdin, nml = debugf)   ! Namelist found. Read it:

  return

end subroutine init_mdin_debugf_dat

subroutine validate_mdin_debugf_dat()

  use charmm_mod, only : charmm_active
  use mdin_ewald_dat_mod, only : int_legal_range
  use gbl_constants_mod, only : error_hdr
  use pmemd_lib_mod, only : mexit

  implicit none

! Local variables

  integer       :: inerr

! Check on bogus data and unsupported options:

  inerr = 0

  call int_legal_range('do_charmm_dump_gold', do_charmm_dump_gold, 0, 1, inerr)

  if (do_charmm_dump_gold /=0 .and. .not. charmm_active) then
    write(mdout, '(a,a)') error_hdr, 'do_charmm_dump_gold requires a CHARMM enabled prmtop.'
    inerr = 1
  end if

! Field any errors and bag out.

  if (inerr .eq. 1) then
    write(mdout, '(/,a)') ' Input errors occurred. Terminating execution.'
    call mexit(6, 1)
  else
    write(mdout, '(a)') ' '
  end if

  return

end subroutine validate_mdin_debugf_dat

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  bcast_mdin_debugf_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine bcast_mdin_debugf_dat

  use parallel_dat_mod

  implicit none

  call mpi_bcast(do_charmm_dump_gold, mdin_debugf_int_cnt, mpi_integer, &
                 0, pmemd_comm, err_code_mpi)

  return

end subroutine bcast_mdin_debugf_dat
#endif

end module mdin_debugf_dat_mod
