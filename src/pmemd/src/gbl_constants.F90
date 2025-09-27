#include "copyright.i"

!*******************************************************************************
!
! Module:  gbl_constants_mod
!
! Description: A central repository of constants that otherwise create
!              circular dependency headaches.
!
!*******************************************************************************

module gbl_constants_mod

  implicit none

! Global constants:

  double precision, parameter   :: PI = 3.1415926535897932384626433832795d0
  double precision, parameter   :: DEG_TO_RAD = PI / 180.d0
  double precision, parameter   :: RAD_TO_DEG = 180.d0 / PI
  double precision, parameter   :: AMBER_ELECTROSTATIC = 18.2223d0
  double precision, parameter   :: ONE_AMBER_ELECTROSTATIC = 1.d0 / &
                                                             AMBER_ELECTROSTATIC
  double precision, parameter   :: LN_TO_LOG = log(10.d0)

  double precision, parameter   :: BOLTZMANN = 8.31441d0
  double precision, parameter   :: JOULES_PER_KCAL = 4.184d3
  double precision, parameter   :: KB = BOLTZMANN / JOULES_PER_KCAL

  double precision, parameter   :: FARADAY = 23.06054801 ! (kcal/mol)/V

  double precision, parameter :: BOHRS_TO_A = 0.529177249D0   ! Bohrs * this = angstroms - Same constants as used in dynamo v2.
  double precision, parameter :: AU_TO_EV = 27.21d0
  double precision, parameter :: EV_TO_KCAL = 23.061d0  !Dynamo's conversion

  ! big_int is largest int that fits in an i8 field, for max nscm, etc.
  integer, parameter    :: big_int = 99999999

  integer, parameter    :: RETIRED_INPUT_OPTION = -10301        ! from sander8
  integer, parameter    :: UNSUPPORTED_INPUT_OPTION = -10302
  integer, parameter    :: NO_INPUT_VALUE = 12344321

  character(11), parameter      :: info_hdr =       '| INFO:    '
  character(11), parameter      :: warn_hdr =       '| WARNING: '
  character(11), parameter      :: error_hdr =      '| ERROR:   '
  character(11), parameter      :: extra_line_hdr = '|          '
! Keep this line so grep finds these version stamps:  Amber24 Amber 24
  character(8), parameter       :: prog_name =      'PMEMD 24'
  character(12), parameter      :: VERSION_STRING = 'Version 24.0'
  character(41), parameter      :: use_sander = &
                                   '|           Please use sander instead.'

end module gbl_constants_mod
