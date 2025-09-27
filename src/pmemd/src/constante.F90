!*******************************************************************************
!
! Module: constante_mod
!
! Description: This module houses the constant Redox potential functionality adapted
!              from the constant pH functionality implemented on PMEMD by Jason Swails
!              (GPU Implementation by Ross C. Walker and Perri Needham) and on Sander
!              by John Mongan.
!
!              Adapted here by Vinicius Wilian D. Cruzeiro
!
!*******************************************************************************

module constante_mod

use constante_dat_mod, only : on_cestep
use constantph_dat_mod, only : on_cpstep, proposed_qterm
use constantph_mod, only : cph_success
use file_io_mod,     only : amopen, nmlsrc
use file_io_dat_mod, only : cein, cein_name, ceout, ceout_name, &
                            cerestrt, cerestrt_name, mdout, owrite, &
                            max_fn_len
use random_mod,      only : random_state, amrand_gen, amrset_gen

implicit none


private ! everything here is private

! Limits on the sizes of the data structures
integer, save :: MAX_TITR_RES, MAX_TITR_STATES, MAX_ATOM_CHRG

! Constant Redox Potential variables

type :: cnste_info
  sequence
  integer :: num_states
  integer :: first_atom
  integer :: num_atoms
  integer :: first_state
  integer :: first_charge
end type cnste_info
integer, parameter :: SIZE_CNSTE_INFO = 5

type(cnste_info), parameter :: NULL_CE_INFO = cnste_info(0,0,0,0,0)
type(random_state), save     :: cnste_randgen

double precision, save, allocatable :: gbl_chrgdat(:)
double precision, save, allocatable :: gbl_statene(:)
double precision, save, allocatable :: gbl_eo_corr(:)

character (len=40), save, allocatable :: gbl_resname(:)

integer, save, allocatable :: gbl_resstate(:)
integer, save, allocatable :: gbl_eleccnt(:)
integer, save, allocatable :: mobile_atoms(:)
integer, save              :: trescnt

integer, save, allocatable :: iselres(:)
integer, save, allocatable :: iselstat(:)

integer, save, public :: ce_igb

double precision, save :: ce_intdiel ! not implemented

logical, save              :: first_ceout ! first pass through ceout file
logical, save              :: first_cerestrt ! first pass through cerestrt file
logical, save              :: ce_success ! to keep track of whether or not a
                                          ! reduction state changed

type(cnste_info), save, allocatable :: gbl_stateinf(:)

#ifdef MPI
public cnste_setup, cnste_begin_step, cnste_end_step, &
       cnste_write_restart, cnste_write_ceout, cnste_read, cnste_bcast, &
       total_reduction, explicit_cnste_begin_step, cnste_explicitmd, &
       cnste_read_limits, cnste_deallocate
#else
public cnste_setup, cnste_begin_step, cnste_end_step, &
       cnste_write_restart, cnste_write_ceout, cnste_read, &
       total_reduction, explicit_cnste_begin_step, cnste_explicitmd, &
       cnste_read_limits, cnste_deallocate
#endif
integer, public, save :: cefirst_sol


contains

!*******************************************************************************
!
! Subroutine: cnstph_read_limits
!
! Description: Reads the limits from the constant pH input file
!
!*******************************************************************************

subroutine cnste_read_limits()

  implicit none

  integer :: ntres, ntstates, natchrg, maxh, ifind

  namelist /cnstphe_limits/ ntres, ntstates, natchrg, maxh

  ! Set default values
  ntres = 50
  ntstates = 200
  natchrg = 1000

  ! Open the cein file and read the cnstphe_limits namelist

  call amopen(cein, cein_name, 'O', 'F', 'R')

  call nmlsrc('cnstphe_limits', cein, ifind)

  if (ifind .ne. 0) then        ! Namelist found. Read it:
    read (cein, nml=cnstphe_limits)
  end if

  ! Setting limits
  MAX_TITR_RES = ntres
  MAX_TITR_STATES = ntstates
  MAX_ATOM_CHRG = natchrg

end subroutine cnste_read_limits

!*******************************************************************************
!
! Subroutine: cnste_read
!
! Description: Reads the constant Redox Potential input file and initializes the data
!
!*******************************************************************************

subroutine cnste_read(num_ints, num_reals)

  use gbl_constants_mod, only : AMBER_ELECTROSTATIC
  use mdin_ctrl_dat_mod
  use pmemd_lib_mod,     only : mexit, get_atomic_number
  use prmtop_dat_mod, only    : natom, atm_gb_fs, atm_atomicnumber, &
            loaded_atm_atomicnumber, atm_igraph, atm_mass, atm_gb_radii

  implicit none

! Passed variables

  integer, intent(in out) :: num_ints
  integer, intent(in out) :: num_reals

! Local variables

  ! Namelist variables -- same as gbl_ counterparts above

  double precision  :: chrgdat(0:MAX_ATOM_CHRG-1)
  double precision  :: statene(0:MAX_TITR_STATES-1)
  double precision  :: eo_corr(0:MAX_TITR_STATES-1)
  type(cnste_info) :: stateinf(0:MAX_TITR_RES-1)

  integer       :: eleccnt(0:MAX_TITR_STATES-1)
  integer       :: resstate(0:MAX_TITR_RES-1)
  character(40) :: resname(0:MAX_TITR_RES)

  integer       :: itres        ! residue iterator
  integer       :: itatm        ! atom iterator
  integer       :: last_res     ! last residue (test for overflow)
  integer       :: last_charge  ! last charge  (test for overflow)
  integer       :: last_state   ! last state   (test for overflow)
  integer       :: atomicnumber ! holder for atomic numbers
  integer       :: i            ! counter

  namelist /cnste/ stateinf, resstate, eleccnt, chrgdat, statene, eo_corr, &
                    trescnt, resname, cefirst_sol, ce_igb, ce_intdiel

  ! Allocate the global constant E data arrays

  call cnste_allocate(num_ints, num_reals)

  ! Zero out our temporary arrays

  call cnste_zero(stateinf, eleccnt, resname, resstate, chrgdat, statene, eo_corr)

  ! Open the cein file and read the cnste namelist

  write(mdout, '(a,a)') '|reading charge increments from file: ', cein_name

  call amopen(cein, cein_name, 'O', 'F', 'R')

  ! Initialize the internal dielectric in case it's not specified

  ce_intdiel = 1.0d0

  read (cein, nml=cnste)

  ! Scale the charges to amber internal units

  do itatm = 0, MAX_ATOM_CHRG - 1
    chrgdat(itatm) = chrgdat(itatm) * AMBER_ELECTROSTATIC
  end do

  ! Check to see if any of our fields are overflowing...

  if (trescnt .gt. MAX_TITR_RES) then
    write(mdout, '(a)') 'Too many titrating residues. Alter ntres in the &
                        &cnstphe_limits namelist'
    call mexit(mdout, 1)
  end if

  do itres = 0, trescnt - 1

    if (stateinf(itres)%first_state + stateinf(itres)%num_states .gt. &
        MAX_TITR_STATES) then
      write(mdout, '(a)') 'Too many titrating states. Alter ntstates in the &
                          &cnstphe_limits namelist'
      call mexit(mdout, 1)
    end if

    if (stateinf(itres)%first_charge + stateinf(itres)%num_atoms .gt. &
        MAX_ATOM_CHRG) then
      write(mdout, '(a)') 'Too much charge data. Alter natchrg in the &
                          &cnstphe_limits namelist'
      call mexit(mdout, 1)
    end if

  end do

  ! Copy the data into the main arrays

  gbl_stateinf(:) = stateinf(:)
  gbl_resstate(:) = resstate(:)
  gbl_eleccnt(:)  = eleccnt(:)
  gbl_chrgdat(:)  = chrgdat(:)
  gbl_statene(:)  = statene(:)
  gbl_eo_corr(:)  = eo_corr(:)
  gbl_resname(:)  = resname(:)

  ! See if we need to set up GB parameters

  if (icnste .eq. 2) then
    gb_cutoff = 1000.d0  ! GB cutoff is infinite here
    if (ce_igb .eq. 1) then
      gb_alpha = 1.d0
      gb_beta = 0.d0
      gb_gamma = 0.d0
    else if (ce_igb .eq. 2) then
      ! Use our best guesses for Onufriev/Case GB  (GB^OBC I):
      gb_alpha = 0.8d0
      gb_beta = 0.d0
      gb_gamma = 2.909125d0
    else if (ce_igb .eq. 5) then
      ! Use our second best guesses for Onufriev/Case GB (GB^OBC II):
      gb_alpha = 1.d0
      gb_beta = 0.8d0
      gb_gamma = 4.85d0
    else if (ce_igb .eq. 7) then
      ! Use parameters for Mongan et al. CFA GBNECK:
      gb_alpha = 1.09511284d0
      gb_beta = 1.90792938d0
      gb_gamma = 2.50798245d0
      gb_neckscale = 0.361825d0
    else if (ce_igb .eq. 8) then
      gb_alpha_h   = 0.788440d0
      gb_beta_h    = 0.798699d0
      gb_gamma_h   = 0.437334d0
      gb_alpha_c   = 0.733756d0
      gb_beta_c    = 0.506378d0
      gb_gamma_c   = 0.205844d0
      gb_alpha_n   = 0.503364d0
      gb_beta_n    = 0.316828d0
      gb_gamma_n   = 0.192915d0
      gb_alpha_os  = 0.867814d0
      gb_beta_os   = 0.876635d0
      gb_gamma_os  = 0.387882d0
      gb_alpha_p   = 1.0d0
      gb_beta_p    = 0.8d0
      gb_gamma_p   = 4.85d0
      gb_neckscale = 0.826836d0
      screen_h = 1.425952d0
      screen_c = 1.058554d0
      screen_n = 0.733599d0
      screen_o = 1.061039d0
      screen_s = -0.703469d0
      screen_p = 0.5d0
      offset  = 0.195141d0
    end if

    ! Set gb_kappa as long as saltcon is .ge. 0.d0 (and catch bug below if not).

    if (saltcon .ge. 0.d0) then

      ! Get Debye-Huckel kappa (A**-1) from salt concentration (M), assuming:
      !   T = 298.15, epsext=78.5,
      gb_kappa = sqrt(0.10806d0 * saltcon)

      ! Scale kappa by 0.73 to account(?) for lack of ion exclusions:

      gb_kappa = 0.73d0 * gb_kappa

    end if

    ! We need to set up the data structures for igb == 7 or igb == 8 since this
    ! is done in subroutine init_prmtop_dat by default, and that was already
    ! called so we would know the value of 'natom'

    if (ce_igb .eq. 7) then

      write(mdout,'(a)') &
        ' Replacing prmtop screening parameters with GBn (igb=7) values'

      do i = 1, natom

        if (loaded_atm_atomicnumber) then
          atomicnumber = atm_atomicnumber(i)
        else
          call get_atomic_number(atm_igraph(i), atm_mass(i), atomicnumber)
        end if

        if (atomicnumber .eq. 6) then
          atm_gb_fs(i) = 4.84353823306d-1
        else if (atomicnumber .eq. 1) then
          atm_gb_fs(i) = 1.09085413633d0
        else if (atomicnumber .eq. 7) then
          atm_gb_fs(i) = 7.00147318409d-1
        else if (atomicnumber .eq. 8) then
          atm_gb_fs(i) = 1.06557401132d0
        else if (atomicnumber .eq. 16) then
          atm_gb_fs(i) = 6.02256336067d-1
        else
          atm_gb_fs(i) = 5.d-1 ! not optimized
        end if

      end do

    else if (ce_igb .eq. 8) then

      write(mdout, '(a)') &
        ' Replacing prmtop screening parameters with GBn2 (igb=8) values'

      do i = 1, natom

        if(loaded_atm_atomicnumber) then
          atomicnumber = atm_atomicnumber(i)
        else
          call get_atomic_number(atm_igraph(i), atm_mass(i), atomicnumber)
        end if

        ! The screen_ variables are found in mdin_ctrl_dat_mod

        if (atomicnumber .eq. 6) then
          atm_gb_fs(i) = screen_c
        else if (atomicnumber .eq. 1) then
          atm_gb_fs(i) = screen_h
        else if (atomicnumber .eq. 7) then
          atm_gb_fs(i) = screen_n
        else if (atomicnumber .eq. 8) then
          atm_gb_fs(i) = screen_o
        else if (atomicnumber .eq. 16) then
          atm_gb_fs(i) = screen_s
        else if (atomicnumber .eq. 15) then
          atm_gb_fs(i) = screen_p
        else
          atm_gb_fs(i) = 5.d-1 ! not optimized
        end if

      end do

    end if ! ce_igb .eq. 7

    if (ce_igb .eq. 7 .or. ce_igb .eq. 8) then

      ! In this case, we need to re-compute atm_gb_fs(i), since we changed its
      ! value above

      gb_fs_max = 0.d0

      do i = 1, natom
        atm_gb_fs(i) = atm_gb_fs(i) * (atm_gb_radii(i) - offset)
        gb_fs_max = max(gb_fs_max, atm_gb_fs(i))
      end do

    end if

    mobile_atoms(1:cefirst_sol-1) = 0
    mobile_atoms(cefirst_sol:natom) = 1
  end if ! (icnste .eq. 2)

  return

end subroutine cnste_read

!*******************************************************************************
!
! Subroutine: cnste_allocate
!
! Description: Allocates constant Redox potential data structures
!
!*******************************************************************************

subroutine cnste_allocate(num_ints, num_reals)

  use constante_dat_mod, only: allocate_cnste_dat
  use mdin_ctrl_dat_mod, only : icnste
  use pmemd_lib_mod,  only    : mexit
  use prmtop_dat_mod, only    : natom

  implicit none

! Passed variables

  integer, intent(in out) :: num_ints
  integer, intent(in out) :: num_reals

! Local variables

  integer :: alloc_failed

  allocate(gbl_chrgdat(0:MAX_ATOM_CHRG-1),    &
           gbl_statene(0:MAX_TITR_STATES-1),  &
           gbl_eo_corr(0:MAX_TITR_STATES-1),  &
           gbl_stateinf(0:MAX_TITR_RES-1),    &
           gbl_eleccnt(0:MAX_TITR_STATES-1),  &
           gbl_resname(0:MAX_TITR_RES),       &
           gbl_resstate(0:MAX_TITR_RES-1),    &
           iselres(0:MAX_TITR_RES),           &
           iselstat(2),                       &
           stat=alloc_failed )

  if (alloc_failed .ne. 0) then
    write(6,'(a)') 'Error in constant Redox potential allocation!'
    call mexit(6, 1)
  end if

  call allocate_cnste_dat(natom, num_reals, alloc_failed)
  if (alloc_failed .ne. 0) then
    write(6,'(a)') 'Error in constant Redox potential allocation!'
    call mexit(6, 1)
  end if

  ! Allocate our mobile_atoms array

  if (icnste .eq. 2) then
    allocate(mobile_atoms(natom), stat=alloc_failed)
    if (alloc_failed .ne. 0) then
      write(6,'(a)') 'Error in constant Redox potential allocation!'
      call mexit(6, 1)
    end if
    num_ints = num_ints + size(mobile_atoms)
  end if

  num_ints = num_ints + size(gbl_stateinf)   &
                      + size(gbl_eleccnt)    &
                      + size(gbl_resname)    &
                      + size(gbl_resstate)   &
                      + size(iselres)        &
                      + size(iselstat)

  num_reals = num_reals + size(gbl_chrgdat)  &
                        + size(gbl_statene)  &
                        + size(gbl_eo_corr)

end subroutine cnste_allocate

!*******************************************************************************
!
! Subroutine: cnste_deallocate
!
! Description: Deallocates constant Redox potential data structures
!
!*******************************************************************************

subroutine cnste_deallocate()

  use constante_dat_mod, only: cleanup_cnste_dat
  use mdin_ctrl_dat_mod, only: icnste
  use pmemd_lib_mod,  only   : mexit

  implicit none

! Local variables

  integer :: dealloc_failed

  deallocate(gbl_chrgdat,  &
             gbl_statene,  &
             gbl_eo_corr,  &
             gbl_stateinf, &
             gbl_eleccnt,  &
             gbl_resname,  &
             gbl_resstate, &
             iselres,      &
             iselstat,     &
             stat=dealloc_failed )

  if (dealloc_failed .ne. 0) then
    write(6,'(a)') 'Error in constant Redox potential deallocation!'
    call mexit(6, 1)
  end if

  call cleanup_cnste_dat()

  ! Deallocate our mobile_atoms array

  if (icnste .eq. 2) then
    deallocate(mobile_atoms, stat=dealloc_failed)
    if (dealloc_failed .ne. 0) then
      write(6,'(a)') 'Error in constant Redox potential deallocation!'
      call mexit(6, 1)
    end if
  end if

end subroutine cnste_deallocate

!*********************************************************************************************
!
! Subroutine: cnste_zero
!
! Description: Zeroes out all of the main data arrays/types used for constant Redox potential
!
!*********************************************************************************************

subroutine cnste_zero(stateinf, eleccnt, resname, resstate, chrgdat, statene, eo_corr)

  use mdin_ctrl_dat_mod, only : icnste

  implicit none

! Passed Variables

  type(cnste_info), intent(out)  :: stateinf(0:MAX_TITR_RES-1)
  integer, intent(out)            :: eleccnt(0:MAX_TITR_STATES-1)
  integer, intent(out)            :: resstate(0:MAX_TITR_RES-1)
  character (len=40), intent(out) :: resname(0:MAX_TITR_RES)
  double precision, intent(out)   :: chrgdat(0:MAX_ATOM_CHRG-1)
  double precision, intent(out)   :: statene(0:MAX_TITR_STATES-1)
  double precision, intent(out)   :: eo_corr(0:MAX_TITR_STATES-1)

  ! cnste_info type
  stateinf(:) = NULL_CE_INFO

  ! Integers
  eleccnt(:)  = 0
  resstate(:) = 0
  iselres(:)  = 0
  iselstat(:) = 0

  ! Chars
  resname(:)  = ' '

  ! Reals
  chrgdat(:)    = 0.d0
  statene(:)    = 0.d0
  eo_corr(:)    = 0.d0

  return

end subroutine cnste_zero

#ifdef MPI
!*******************************************************************************
!
! Subroutine: cnste_bcast
!
! Description: Broadcasts all constant Redox potential data to the whole MPI universe
!
!*******************************************************************************

subroutine cnste_bcast(num_ints, num_reals)

  use mdin_ctrl_dat_mod, only : icnste
  use prmtop_dat_mod, only : atm_qterm, natom
  use parallel_dat_mod  ! Provides MPI routines/constants, pmemd_comm and master

  implicit none

! Passed variables

  integer, intent(in out) :: num_ints
  integer, intent(in out) :: num_reals

  ! The slave nodes haven't been allocated yet...

  call mpi_bcast(MAX_TITR_RES, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)
  call mpi_bcast(MAX_TITR_STATES, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)
  call mpi_bcast(MAX_ATOM_CHRG, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)
  if (.not. master) then
    call cnste_allocate(num_ints, num_reals)
  end if

  call mpi_bcast(gbl_stateinf, MAX_TITR_RES * SIZE_CNSTE_INFO, mpi_integer, &
                 0, pmemd_comm, err_code_mpi)

  call mpi_bcast(gbl_resstate, MAX_TITR_RES, mpi_integer, 0, pmemd_comm, &
                 err_code_mpi)

  call mpi_bcast(gbl_eleccnt, MAX_TITR_STATES, mpi_integer, 0, pmemd_comm, &
                 err_code_mpi)

  call mpi_bcast(trescnt, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)

  call mpi_bcast(gbl_resname, 40 * (MAX_TITR_RES + 1), mpi_character, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(gbl_chrgdat, MAX_ATOM_CHRG, mpi_double_precision, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(gbl_statene, MAX_TITR_STATES, mpi_double_precision, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(gbl_eo_corr, MAX_TITR_STATES, mpi_double_precision, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(ce_igb, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)
  call mpi_bcast(cefirst_sol, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)
  call mpi_bcast(ce_intdiel, 1, mpi_double_precision, 0, &
                 pmemd_comm, err_code_mpi)
  if (icnste .eq. 2) &
    call mpi_bcast(mobile_atoms, natom, mpi_integer, 0, &
                   pmemd_comm, err_code_mpi)

end subroutine cnste_bcast
#endif /* MPI */

!*******************************************************************************
!
! Subroutine: cnste_setup
!
! Description:
!
! Sets up constant Redox potential calculation and launching the random number generator.
! This should be called by all threads.
!
!*******************************************************************************

subroutine cnste_setup(crd)

  use mdin_ctrl_dat_mod, only : ig
  use parallel_dat_mod
  use prmtop_dat_mod,    only : natom, atm_qterm, gbl_one_scee
  use extra_pnts_nb14_mod

  implicit none

! Passed variables

  double precision, intent(in) :: crd(3, natom)

! Local variables

  integer  :: itres   ! residue counter
  integer  :: itstate ! state counter
  integer  :: itatm   ! atom counter
  integer  :: i       ! loop counter
  integer  :: randseed! Random number seed

  ! Open up the ceout file if you're the master

  if (master) then
    ! For some reason, owrite is set to 'U' given the -O flag, but this just
    ! appends to the ceout file instead of overwriting it. Passing 'R' to amopen
    ! overwrites the file (this is the behavior in sander). Since this has been
    ! around since Bob first wrote the code, I don't want to remove it. Hence
    ! the kludge here.
    if (owrite .eq. 'U') then
      call amopen(ceout, ceout_name, 'R', 'F', 'W')
    else
      call amopen(ceout, ceout_name, owrite, 'F', 'W')
    end if
  end if

  ! Set initial charges based on initial states

  do itres = 0, trescnt - 1
    call cnste_update_charge(atm_qterm, itres, gbl_resstate(itres))
  end do

#ifdef CUDA
    call gpu_refresh_charges(cit_nb14, gbl_one_scee, atm_qterm)
#endif

  ! Now copy this set of charges to the proposed charge set

  proposed_qterm(:) = atm_qterm(:)

  ! Set up the random number generator

  randseed = ig

#ifdef MPI
  call mpi_bcast(randseed, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)
#endif

  call amrset_gen(cnste_randgen, randseed)

  ! Set up our first time through

  first_ceout = .true.
  first_cerestrt = .true.

  ! We don't do an exchange attempt on our first step

  on_cestep = .false.

end subroutine cnste_setup

!*******************************************************************************
!
! Subroutine: cnste_begin_step
!
! Description:
!
! Begins the constant Redox potential step -- selects the titratable residue(s) that will be
! titrated as well as the state that will be attempted.
!
! We randomly select a residue to titrate. Update the iselstat/iselres arrays to list which
! residues we're changing to which states.
!
!*******************************************************************************

subroutine cnste_begin_step

  use prmtop_dat_mod, only : natom, gbl_one_scee
  use extra_pnts_nb14_mod
  use mdin_ctrl_dat_mod, only : icnste
  implicit none

  integer          :: i
  double precision :: random_value

  ! No multisite

  iselres(0) = 1

  ! Select a random residue

  call amrand_gen(cnste_randgen, random_value)
  iselres(1) = int( (random_value * 0.99999999d0) * trescnt)

  ! Select a random, unique state

  call amrand_gen(cnste_randgen, random_value)
  iselstat(1) = int( (random_value * 0.99999999d0) * &
                     (gbl_stateinf(iselres(1))%num_states - 1) )

  if (iselstat(1) .ge. gbl_resstate(iselres(1))) &
    iselstat(1) = iselstat(1) + 1

  ! Update our proposed state charge array

  call cnste_update_charge(proposed_qterm, iselres(1), iselstat(1))
#ifdef CUDA
  if (icnste .eq. 2)then
    call gpu_refresh_charges(cit_nb14, gbl_one_scee, proposed_qterm)
  endif
#endif

end subroutine cnste_begin_step

!*******************************************************************************
!
! Subroutine: cnste_end_step
!
! Description: Ends the constant Redox potential step and accepts/rejects the MC move
!
!*******************************************************************************

subroutine cnste_end_step(dvdl, tresi, holder_natom)

  use gbl_constants_mod, only : KB, FARADAY
  use mdin_ctrl_dat_mod, only : solve, temp0, icnste
  use prmtop_dat_mod
  use charmm_mod
  use extra_pnts_nb14_mod

  implicit none

! Passed variables

  double precision, intent(in) :: dvdl

! Local variables

  double precision :: random_value
  double precision :: deltae
  double precision :: metrop

  integer          :: statebase
  integer          :: i, tresi, holder_natom

  deltae = 0.d0

  do i = 1, iselres(0)

    statebase = gbl_stateinf(iselres(i))%first_state

    ! deltae = delta G ref (proposed state) - delta G ref (current state)

    deltae = deltae + gbl_statene(iselstat(i) + statebase) - &
                      gbl_statene(gbl_resstate(iselres(i)) + statebase)

    ! Adjust for the Eo term (Eoref * FARADAY)

    deltae = deltae + (gbl_eo_corr(iselstat(i) + statebase) - &
                      gbl_eo_corr(gbl_resstate(iselres(i)) + statebase)) * &
                      FARADAY

    ! Adjust for E (delta electrons * E * FARADAY)

    deltae = deltae - (gbl_eleccnt(iselstat(i) + statebase) - &
                       gbl_eleccnt(gbl_resstate(iselres(i)) + statebase)) * &
                       solve * FARADAY

  end do

  call amrand_gen(cnste_randgen, random_value)

  metrop = exp( (deltae - dvdl) / (KB * temp0) )

! write(0, '(a,f16.5)') 'deltae  = ', deltae
! write(0, '(a,f16.5)') 'dvdl    = ', dvdl
! write(0, '(a,f16.5)') 'randval = ', random_value
! write(0, '(a,f16.5)') 'metrop  = ', metrop
! write(0, '(72("="))')

! do i = 1, natom
!   write(0, '(a,i3,a,f16.5,a,i3,a,f16.5,a,f16.5)') ' proposed_qterm(',i,') = ', proposed_qterm(i), &
!       ' atm_qterm(',i,') = ', atm_qterm(i), ' diff = ',(proposed_qterm(i)-atm_qterm(i))
! end do
! write(0, '(72("="))')

  if (random_value .le. metrop) then
    do i = 1, iselres(0)
      gbl_resstate(iselres(i)) = iselstat(i)
      call cnste_update_charge(atm_qterm, iselres(i), iselstat(i))
    end do

    ce_success = .true.

#ifdef CUDA
    if (icnste .eq. 2) then
       if (trescnt .eq. tresi) then
          call gpu_update_natoms(holder_natom, .true.)
          call gpu_upload_charges_pme_cph(proposed_qterm)
          call gpu_refresh_charges(cit_nb14, gbl_one_scee, proposed_qterm)
       endif
    else
       call gpu_refresh_charges(cit_nb14, gbl_one_scee, proposed_qterm)
    endif
#endif

  else
    ! Revert charges if we rejected the change

    do i = 1, iselres(0)
      call cnste_update_charge(proposed_qterm, iselres(i), &
                                gbl_resstate(iselres(i)) )
    end do
#ifdef CUDA
    if (icnste .eq. 2) then
       if (trescnt .eq. tresi) then
          call gpu_update_natoms(holder_natom, .true.)
          call gpu_upload_charges_pme_cph(atm_qterm)
          call gpu_refresh_charges(cit_nb14, gbl_one_scee, atm_qterm)
       endif
    endif
#endif

  end if

end subroutine cnste_end_step

!*******************************************************************************
!
! Subroutine: cnste_update_charge
!
! Description: Update the given charge array with proposed states
!
!*******************************************************************************

subroutine cnste_update_charge(chg_arry, selres, selstat)

  use prmtop_dat_mod, only : natom

  implicit none

! Passed variables

  double precision, intent(in out) :: chg_arry(natom) ! charge array to change
  integer, intent(in)              :: selres          ! selected residue
  integer, intent(in)              :: selstat         ! selected state

! Local variables

  integer :: itatm ! atom iterator

  do itatm = 0, gbl_stateinf(selres)%num_atoms - 1

    chg_arry(gbl_stateinf(selres)%first_atom + itatm) = &
      gbl_chrgdat(gbl_stateinf(selres)%first_charge + &
      selstat * gbl_stateinf(selres)%num_atoms + itatm)
  end do

end subroutine cnste_update_charge

!*******************************************************************************
!
! Subroutine: cnste_write_ceout
!
! Description: Writes constant Redox potential titration information to ceout file
!
!*******************************************************************************

subroutine cnste_write_ceout(nstep, nstlim, time, remd_method, remd_types, replica_indexes)

  use mdin_ctrl_dat_mod, only : ntwx, solve, solvph, temp0, ntcnste
  use parallel_dat_mod, only  : master

  implicit none

! Passed variables

  integer, intent(in) :: nstep
  integer, intent(in) :: nstlim
  integer, intent(in) :: remd_method
  integer, dimension(:), intent(in) :: remd_types, replica_indexes

  double precision, intent(in) :: time

! Local variables

  logical :: full
  integer :: i
  character(len=300) :: txt

  ! Only the master writes ceout info

  if (.not. master) return

  if (ntwx .gt. 0) then
    full = (mod(nstep, ntwx) .eq. 0 .or. nstep .eq. nstlim .or. first_ceout)
  else
    full = (nstep .eq. nstlim .or. first_ceout)
  end if

   write(txt, '()')
   if (remd_method>0) then
      if (remd_method == 1) then
         write(txt, '(a,f7.2,a)') ' T: ', temp0, ' K'
      else if (remd_method == 3) then
         write(txt, '(a,i3)') ' H: ', replica_indexes(1)
      else if (remd_method == 4) then
         write(txt, '(a,f7.3)') ' pH: ', solvph
      else if (remd_method == 5) then
         write(txt, '(a,f12.7,a)') ' E: ', solve, ' V'
      end if
   else if (remd_method == -1) then
      do i = 1, size(remd_types)
         if (remd_types(i) == 1) then
            write(txt, '(a,a,f7.2,a)') trim(txt), ' T: ', temp0, ' K'
         else if (remd_types(i) == 3) then
            write(txt, '(a,a,i3)') trim(txt), ' H: ', replica_indexes(i)
         else if (remd_types(i) == 4) then
            write(txt, '(a,a,f7.3)') trim(txt), ' pH: ', solvph
         else if (remd_types(i) == 5) then
            write(txt, '(a,a,f12.7,a)') trim(txt), ' E: ', solve, ' V'
         end if
      end do
   end if
  if (full) then

    write(ceout, '(a,f12.7,a,f7.2,a)') 'Redox potential: ', solve, ' V Temperature: ', temp0, ' K'
    write(ceout, '(a,i8)')   'Monte Carlo step size: ', ntcnste
    write(ceout, '(a,i8)')   'Time step: ', nstep
    write(ceout, '(a,f14.3)') 'Time: ', time

    do i = 0, trescnt - 1
      write(ceout, '(a,i4,a,i2,a)') &
          'Residue ', i, ' State: ', gbl_resstate(i), trim(txt)
    end do

  else

    do i = 1, iselres(0)
      write(ceout, '(a,i4,a,i2,a)') 'Residue ', iselres(i), ' State: ', &
            gbl_resstate(iselres(i)), trim(txt)
    end do

  end if ! full

  write(ceout, '()')

  ! No longer first time around

  first_ceout = .false.

end subroutine cnste_write_ceout

!*******************************************************************************
!
! Subroutine: cnste_write_restart
!
! Description: writes the cerestrt file
!
!*******************************************************************************

subroutine cnste_write_restart(nstep)

  use gbl_constants_mod, only : ONE_AMBER_ELECTROSTATIC
  use parallel_dat_mod, only  : master
  use pmemd_lib_mod,     only : mexit
  use mdin_ctrl_dat_mod, only : ntwr

  implicit none

  integer, intent(in) :: nstep

  ! Namelist variables -- same as gbl_ counterparts above

  double precision  :: chrgdat(0:MAX_ATOM_CHRG-1)
  double precision  :: statene(0:MAX_TITR_STATES-1)
  double precision  :: eo_corr(0:MAX_TITR_STATES-1)
  type(cnste_info) :: stateinf(0:MAX_TITR_RES-1)

  integer       :: eleccnt(0:MAX_TITR_STATES-1)
  integer       :: resstate(0:MAX_TITR_RES-1)
  character(40) :: resname(0:MAX_TITR_RES)

  character(7)  :: stat

  integer        :: i, istart, iend
  character(12)  :: num
  character(256) :: cerestrt2_name

  integer :: ntres, ntstates, natchrg

  namelist /cnstphe_limits/ ntres, ntstates, natchrg

  namelist /cnste/ stateinf, resstate, eleccnt, chrgdat, statene, eo_corr, &
                    trescnt, resname, cefirst_sol, ce_igb, ce_intdiel

  if (.not. master) return

  ntres = MAX_TITR_RES
  ntstates = MAX_TITR_STATES
  natchrg = MAX_ATOM_CHRG

  chrgdat(:) = gbl_chrgdat(:) * ONE_AMBER_ELECTROSTATIC
  statene(:) = gbl_statene(:)
  eo_corr(:) = gbl_eo_corr(:)
  stateinf(:) = gbl_stateinf(:)
  eleccnt(:) = gbl_eleccnt(:)
  resstate(:) = gbl_resstate(:)
  resname(:) = gbl_resname(:)

  if (first_cerestrt) then
    if (owrite .eq. 'N') then
      stat = 'NEW'
    else if (owrite .eq. 'O') then
      stat = 'OLD'
    else if (owrite .eq. 'R') then
      stat = 'REPLACE'
    else if (owrite .eq. 'U') then
      stat = 'UNKNOWN'
    end if
    open(unit=cerestrt, file=cerestrt_name, status=stat, form='FORMATTED', &
         delim='APOSTROPHE', err=666)
    first_cerestrt = .false.
  else
    open(unit=cerestrt, file=cerestrt_name, status='OLD', form='FORMATTED', &
         delim='APOSTROPHE')
  end if

  write(cerestrt, nml=cnstphe_limits)
  write(cerestrt, nml=cnste)

  close(cerestrt)

  ! Consider whether to save 2ndary restrt:

  if (ntwr .ge. 0) return

  do iend = 1, max_fn_len
    if (cerestrt_name(iend:iend) .le. ' ') exit
  end do

  iend = iend - 1

  write(num, '(i12)') nstep

  do istart = 1, 12
    if (num(istart:istart) .ne. ' ') exit
  end do

  write(cerestrt2_name, '(a,a,a)') cerestrt_name(1:iend), '_', num(istart:12)

  if (owrite .eq. 'N') then
    stat = 'NEW'
  else if (owrite .eq. 'O') then
    stat = 'OLD'
  else if (owrite .eq. 'R') then
    stat = 'REPLACE'
  else if (owrite .eq. 'U') then
    stat = 'UNKNOWN'
  end if
  open(unit=cerestrt, file=cerestrt2_name, status=stat, form='FORMATTED', &
       delim='APOSTROPHE', err=667)

  write(cerestrt, nml=cnste)

  close(cerestrt)

  return

666 write(mdout, '(a)') 'Error opening ', trim(cerestrt_name)
    call mexit(mdout, 1)

667 write(mdout, '(a)') 'Error opening ', trim(cerestrt2_name)
    call mexit(mdout, 1)

end subroutine cnste_write_restart

!*******************************************************************************
!
! Function: total_reduction
!
! Description: Finds the total number of 'active' electrons
!
!*******************************************************************************

integer function total_reduction()

  implicit none

  integer :: i

  total_reduction = 0

  do i = 0, trescnt - 1
    total_reduction = total_reduction + &
             gbl_eleccnt(gbl_stateinf(i)%first_state + gbl_resstate(i))
  end do

  return

end function total_reduction

!*******************************************************************************
!
! Subroutine: explicit_cnste_begin_step
!
! Description: The beginning of a constant Redox potential step in explicit solvent
!
!*******************************************************************************

subroutine explicit_cnste_begin_step(dcharge, done_res, num_done)

  implicit none

! Passed variables

  double precision, intent(in out) :: dcharge(1:*)
  integer, intent (in out)         :: done_res(1:MAX_TITR_RES)
  integer, intent(in)              :: num_done

! Local variables

  integer :: i
  integer :: j
  integer :: holder_res(1:MAX_TITR_RES)

  double precision :: randval

! This subroutine will randomly select what to titrate, as long as it has not
! yet been titrated. That is the only difference from cnste_begin_step

  iselres(0) = 1

  ! Select a random residue
  call amrand_gen(cnste_randgen, randval)
  iselres(1) = int((randval * 0.9999999d0) * (trescnt - num_done))

  if (num_done .eq. 0) then
    done_res(1) = iselres(1)
  else

    ! Now move the selected residue off of those that have been selected already
    ! THIS ARRAY SHOULD ALREADY BE SORTED (see below)

    do i = 1, num_done
      if (iselres(1) .ge. done_res(i)) iselres(1) = iselres(1) + 1
    end do

    ! If the residue we chose was larger than our last number, add it to the end
    if (iselres(1) .gt. done_res(num_done)) &
      done_res(num_done + 1) = iselres(1)

  end if ! (num_done .eq. 0)

  do i = 1, num_done
    if (iselres(1) .lt. done_res(i)) then
      holder_res(i:num_done) = done_res(i:num_done)
      done_res(i) = iselres(1)
      done_res(i+1:num_done+1) = holder_res(i:num_done)
      exit
    end if
  end do

  ! Select a random, unique state

  call amrand_gen(cnste_randgen, randval)
  iselstat(1) = int( (randval * 0.99999999d0) * &
                     (gbl_stateinf(iselres(1))%num_states - 1) )

  if (iselstat(1) .ge. gbl_resstate(iselres(1))) &
    iselstat(1) = iselstat(1) + 1

  ! Update our proposed state charge array

  call cnste_update_charge(proposed_qterm, iselres(1), iselstat(1))

  return

end subroutine explicit_cnste_begin_step

!*******************************************************************************
!
! Subroutine: cnste_explicitmd
!
! Description: Do MC cycles in explicit solvent
!
!*******************************************************************************

#ifdef MPI
subroutine cnste_explicitmd(atm_cnt, crd, mass, frc, vel, last_vel, &
                             nstep, nstlim, time, my_atm_lst)
#else
subroutine cnste_explicitmd(atm_cnt, crd, mass, frc, vel, last_vel, &
                             nstep, nstlim, time)
#endif

  use gb_force_mod, only       : gb_ce_ene
  use mdin_ctrl_dat_mod, only  : ntp, ntb, ips, igb, nmropt, nscm, ntrelax, ntrelaxe, iamd, igamd
  use nmr_calls_mod, only      : nmrdcp
#ifdef MPI
  use parallel_mod, only       : mpi_allgathervec
#endif
  use parallel_dat_mod, only   : master
  use prmtop_dat_mod, only     : natom, atm_qterm, atm_gb_radii, gbl_one_scee
  use remd_mod, only           : remd_method, remd_types, replica_indexes
  use relaxmd_mod, only        : relaxmd
  use extra_pnts_nb14_mod

  implicit none

  ! Passed variables

#ifdef MPI
  integer, intent(in)     :: my_atm_lst(*)
#endif
  integer, intent(in)     :: nstep
  integer, intent(in)     :: nstlim
  integer, intent(in out) :: atm_cnt
  double precision, intent(in)     :: time
  double precision, intent(in out) :: crd(3, atm_cnt)
  double precision, intent(in)     :: mass(atm_cnt)
  double precision, intent(in out) :: frc(3, atm_cnt)
  double precision, intent(in out) :: vel(3, atm_cnt)
  double precision, intent(in out) :: last_vel(3, atm_cnt)

  ! Local variables

  integer :: holder_ntb, holder_natom, holder_ntp, holder_ips, j, i, k, &
             natom3, selres_holder(1:MAX_TITR_RES), holder_nscm, holder_iamd, holder_igamd

#ifdef CUDA
  integer :: fat
#endif

  double precision, dimension(3, atm_cnt) :: vtemp, last_vtemp

  double precision :: dvdl_current
  double precision :: dvdl_proposed
  double precision :: dvdl

  logical          :: any_success
  ! This subroutine is the main driver for running constant Redox potential MD in explicit
  ! solvent. The first thing it does is call being_step to set up the MC.
  ! Then it removes the periodic boundary conditions (PBC), selects a GB model,
  ! then calls force with nstep = 0 which will force oncpstep to be .true., so
  ! we'll get dvdl back. This is then passed to cnsteendstep to evaluate the
  ! transition. If it fails, then we restore the PBC, turn off the GB model,
  ! restore the CUToff, and return to the calling routine. If it's successful,
  ! we still return the above variables, but then we also turn on belly and fix
  ! the protein while we call runmd to relax the solvent for ntrelaxe steps.
  ! Some important things we do are:
  !  o  turn off trajectory/energy printing
  !  o  zero-out the v and vold arrays for the solute
  !  o  turn on belly and assign belly arrays

! Collect all of the coordinates

#ifdef MPI
  ! Collect all coordinates
  call mpi_allgathervec(atm_cnt, crd)
#endif

  natom3          = atm_cnt * 3
  holder_natom    = atm_cnt
  holder_ntb      = ntb
  holder_ntp      = ntp
  holder_ips      = ips

  natom   = cefirst_sol - 1
  atm_cnt = cefirst_sol - 1
  ntb     = 0
  igb     = ce_igb

  dvdl_current = 0.d0
  dvdl_proposed = 0.d0

  ! Zero-out the holder for the selected residues

  selres_holder(1:trescnt) = 0

  ! Initiate trial moves for reduction states
#ifdef CUDA
  call gpu_download_crd(crd)
  call gpu_update_natoms(natom, .false.)
  call gpu_upload_crd_gb_cph(crd)
#endif

  any_success = .false.

  ! Get the initial energy
#ifdef CUDA
    call gpu_upload_charges_gb_cph(atm_qterm)
#endif
  call gb_ce_ene(atm_cnt, crd, atm_qterm, dvdl_current, .false.)

  do i = 1, trescnt

    ce_success = .false.

    call explicit_cnste_begin_step(proposed_qterm, selres_holder, i-1)

#ifdef CUDA
! These routines can be combined in the next GPU-exCE version
    call gpu_upload_charges_gb_cph(proposed_qterm)
#endif

    ! Call gb force to get energies
    call gb_ce_ene(atm_cnt, crd, proposed_qterm, dvdl_proposed, .true.)

    dvdl = dvdl_proposed - dvdl_current

    ! End the step
    call cnste_end_step(dvdl, i, holder_natom)

    ! If we succeeded, update our 'current' state energy
    if (ce_success) &
      dvdl_current = dvdl_proposed

    any_success = any_success .or. ce_success

    ! Decrement the nmropt counter

    if (nmropt .ne. 0) call nmrdcp

  end do

  ! Restore PBC variables
  natom   = holder_natom
  atm_cnt = holder_natom
  ntb     = holder_ntb
  igb     = 0
  ntp     = holder_ntp
  ips     = holder_ips

  ! Use iselres to stipulate that every residue is printed out

  iselres(0) = trescnt
  do i = 1, trescnt
    iselres(i) = i - 1
  end do

  if (master) call cnste_write_ceout(nstep, nstlim, time, remd_method, remd_types, replica_indexes)

#ifdef CUDA
  call gpu_upload_crd(crd)
#endif

  ! If we did not succeed, bail out
  if (.not. ((on_cpstep .and. cph_success) .or. any_success)) return

  holder_nscm     = nscm
  holder_iamd     = iamd
  holder_igamd    = igamd
  nscm    = 0
  iamd    = 0
  igamd    = 0

  ! If we succeeded, adjust the variables for relaxation dynamics. Don't print
  ! during relaxation steps. Also zero the velocity array, and don't remove COM
  ! motion

#ifdef CUDA
  fat = cefirst_sol - 1
  call gpu_download_vel(vel)
  call gpu_set_first_update_atom(fat)
#endif

#ifdef MPI
  call mpi_allgathervec(atm_cnt, vel)
  call mpi_allgathervec(atm_cnt, last_vel)
#endif

  do i = 1, cefirst_sol - 1
    vtemp(:, i) = vel(:, i)
    last_vtemp(:, i) = last_vel(:, i)
    vel(:, i) = 0.d0
  end do
#ifdef CUDA
  call gpu_upload_vel(vel)
#endif
  ! Call routine to do relaxation dynamics
  if (on_cpstep .and. cph_success) then
    i = MAX(ntrelax,ntrelaxe)
  else
    i = ntrelaxe
  end if
#ifdef MPI
  call relaxmd(atm_cnt, crd, mass, frc, vel, last_vel, my_atm_lst, &
               mobile_atoms, i)
  call mpi_allgathervec(atm_cnt, crd)
  call mpi_allgathervec(atm_cnt, vel)
  call mpi_allgathervec(atm_cnt, last_vel)
#else
  call relaxmd(atm_cnt, crd, mass, frc, vel, last_vel, &
               mobile_atoms, i)
#endif
  ! Restore the original velocities
#ifdef CUDA
  call gpu_download_vel(vel)
#endif

  do i = 1, cefirst_sol - 1
    vel(:,i) = vtemp(:,i)
    last_vel(:,i) = last_vtemp(:,i)
  end do

#ifdef CUDA
  fat = 0
  call gpu_upload_vel(vel)
  call gpu_set_first_update_atom(fat)
#endif

  nscm = holder_nscm
  iamd = holder_iamd
  igamd = holder_igamd

  return

end subroutine cnste_explicitmd

end module constante_mod
