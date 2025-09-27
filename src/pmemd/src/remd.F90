#include "copyright.i"

!*******************************************************************************
!
! Module: remd_mod
!
! Description: 
!
! Module for controlling replica exchange functionality. Many parts were adapted
! from sander's remd module for pmemd. The sander implementation was done by
! Daniel Roe, based on the original implementation by Guanglei Cui and Carlos
! Simmerling. 
!
! This is being designed such that exchanges can be done in an arbitrary number
! of dimensions. This is accomplished through the use of several data 
! structures that track information for each dimension. All exchanges will be
! performed with the replica *higher* in the replica ladder, with the replicas
! performing the exchanges alternating between "even" and "odd" replicas (to
! get every possible exchange represented).
!
! If you want to add a new exchange type that is visible by the multi-D REMD,
! you'll have to provide a new integer index for it (see remd_types description)
! and add it into all of the select case() blocks. You will also have to add a
! standalone setup routine for it if you want it accessible via -rem 4 for
! instance. If you are effectively swapping Hamiltonians, then you end up moving
! around in replica-space, so you will have to swap group assignments and 
! indexes with your neighbor (see how temperature exchanges are handled). If you
! instead swap coordinates, you stand still in replica-space, so this additional
! bookkeeping is unnecessary (see how hamiltonian exchanges are handled).
! Careful, though, since when coupled with dimensions that swap state indices,
! the replica order could change from what you expect.
!
! CPU Implementation by Jason Swails
! GPU Implementation by Ross C. Walker
!              
!*******************************************************************************

module remd_mod

  use mdin_ctrl_dat_mod
  use file_io_dat_mod
  use parallel_dat_mod
  use random_mod, only : random_state, amrset_gen, amrand_gen

  implicit none

  integer, save :: remd_method = 0 ! REMD method to use
  integer, save :: rremd_type = 0 ! RREMD type to use
  integer, save :: remd_random_partner = 0 ! Decide use random partners (0 for regular, and any other value for random)

  double precision, parameter :: BIG = 9999999.d0

#ifndef MPI
! remd_types and replica_indexes are required by constantph.F90 and constante.F90
  integer, save    :: remd_types(1) = 0
  integer, save    :: replica_indexes(1) = 0
#endif

#ifdef MPI
  
  ! Everything is public by default

! Variables
!
!  exchange_successes   : used to track acceptance ratio
!  partners             : which replicas this one can attempt exchanges with
!  remd_dimension       : how many dimensions we will attempt to exchange in
!  control_exchange     : Determines whether the replica controls the exchange
!  even_exchanges       : Does the even replica exchange or not?
!  crd_temp             : temporary coordinate array for hamiltonian exchanges
!  frc_temp             : temporary force array for hamiltonian exchanges
!  temperatures         : array with temperatures for every replica
!                         exchange or let our partner do it
!  phs                  : array with phs for every replica
!                         exchange or let our partner do it
!  redoxs               : array with redox potentials for every replica
!                         exchange or let our partner do it
!  total_left_fe        : calculated FEP free energy for exchanging to the left
!  total_right_fe       : calculated FEP free energy for exchanging to the right
!  num_right_exchg      : # of times we've exchanged to the right
!  num_left_exchg       : # of times we've exchanged to the left
!  remd_modwt           : Forces modwt to re-read the temperatures, etc., or the
!                         NMR restraint facility will overwrite temp0
!  use_pv               : If true enable press/vol correction to exch calc.
!  exchange_ucell       : If true, cell info needs to be exchanged during HREMD.
!  group_num            : Stores *my* group number in each dimension
!  remd_repidx          : Overall replica index for this replica.
!  remd_crdidx          : Overall coordinate index for this replica.
!  replica_indexes      : A collection of which ranks this replica is in each of
!                         its REMD dimensions, beginning from 1
!  index_list           : The list of replica_indexes for every replica in the
!                         dimension we are exchanging in. Set in set_partners
!  remd_types           : Which type of REMD each dimension is. More types can
!                         be added easily, just add to the table below
!                         TEMPERATURE : 1
!                         HAMILTONIAN : 3
!                         PH :          4
!                         REDOX :       5

! General for all REMD methods

  integer, allocatable, save    :: exchange_successes(:,:)
  
  integer, save                 :: partners(2)
  integer, save                 :: remd_dimension
  integer, save                 :: remd_repidx
  integer, save                 :: remd_crdidx

  integer, allocatable, save    :: group_num(:)
  integer, allocatable, save    :: replica_indexes(:)
  integer, allocatable, save    :: index_list(:)
  integer, allocatable, save    :: remd_types(:)

  logical, allocatable, save    :: control_exchange(:)
  logical, allocatable, save    :: even_exchanges(:)

  logical, save                 :: remd_modwt
  logical, save                 :: use_pv = .false.
  logical, save                 :: exchange_ucell = .false.
  logical, save                 :: enforce_reset_velocities = .false.
  type(random_state), save      :: remd_randgen
  
! Variables for pressure/volume correction to exchange calc
  double precision, save        :: remd_pressure = 0.d0
  double precision, save        :: remd_volume   = 0.d0
  
! Specific to T-REMD

  double precision, allocatable :: temperatures(:)

! Specific to pH-REMD

  double precision, allocatable :: phs(:)

! Specific to E-REMD

  double precision, allocatable :: redoxs(:)

! Specific to RXSGLD 

  double precision, allocatable :: tempsgtable(:),sgfttable(:),sgfftable(:),sgfgtable(:)

! Specific to H-REMD

  double precision, allocatable :: crd_temp(:,:)
  double precision, allocatable :: frc_temp(:,:)

  double precision, allocatable, save :: total_left_fe(:,:,:)
  double precision, allocatable, save :: total_right_fe(:,:,:)

  integer, allocatable, save    :: num_right_exchg(:,:,:)
  integer, allocatable, save    :: num_left_exchg(:,:,:)

! Specific to Reservoir REMD
  double precision, allocatable :: rremd_crd(:,:)
  double precision, allocatable :: rremd_vel(:,:)

! Specific to Hybrid Solvent REMD
  double precision, allocatable :: hybridsolvent_remd_crd(:,:)
  double precision, allocatable :: hybridsolvent_remd_vel(:,:)
  double precision, allocatable :: hybridsolvent_remd_frc(:,:)

! Bookkeeping for multi-D REMD
  
  type :: remlog_data
    sequence
    double precision :: scaling         ! T-REMD
    double precision :: real_temp       ! T-REMD
    double precision :: new_temp0       ! T-REMD
    double precision :: struct_num      ! T-REMD
    double precision :: pot_ene_tot     ! T-REMD & H-REMD
    double precision :: temp0           ! T-REMD & H-REMD
    double precision :: success_ratio   ! T-REMD & H-REMD & pH-REMD & E-REMD
    double precision :: num_rep         ! T-REMD & H-REMD & pH-REMD & E-REMD
    double precision :: group_num       ! T-REMD & H-REMD & pH-REMD & E-REMD
    double precision :: left_fe         !          H-REMD
    double precision :: right_fe        !          H-REMD
    double precision :: neighbor_rep    !          H-REMD
    double precision :: nei_pot_ene     !          H-REMD
    double precision :: success         !          H-REMD & pH-REMD & E-REMD
    double precision :: repnum          !          H-REMD
    double precision :: left_exchg      !          H-REMD
    double precision :: right_exchg     !          H-REMD
    double precision :: my_ph           !                   pH-REMD
    double precision :: nei_ph          !                   pH-REMD
    double precision :: nprot           !                   pH-REMD
    double precision :: my_e            !                             E-REMD
    double precision :: nei_e           !                             E-REMD
    double precision :: nelec           !                             E-REMD
  end type remlog_data

  integer, parameter :: SIZE_REMLOG_DATA = 23

  type(remlog_data), allocatable, save :: multid_print_data(:)
  type(remlog_data), allocatable, save :: multid_print_data_buf(:,:)
  type(remlog_data), parameter         :: NULL_REMLOG_DATA = &
    remlog_data(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                0.d0,0.d0)

                

contains

!*******************************************************************************
!
! Subroutine: setup_pv_correction 
!
! Description: Determine if the pressure-volume correction should be used. 
!
!*******************************************************************************
subroutine setup_pv_correction(outu, ntpin, master)
  ! USE STATEMENTS
  ! ARGUMENTS
  implicit none
  integer, intent(in) :: outu, ntpin
  logical, intent(in) :: master

  if (ntpin > 0) then
    exchange_ucell = .true. ! Only matters for hremd_exchange
    use_pv = .true.
    if (master) &
      write(outu,'(a)') '| REMD: Pressure/volume correction to exchange calc active for TREMD/HREMD.'
  else
    exchange_ucell = .false.
    use_pv = .false.
  endif
end subroutine setup_pv_correction

!*******************************************************************************
!
! Subroutine: remd_setup
!
! Description: Sets up the REMD run. Opens the log file, sets up exchange tables
!              etc.
!
!*******************************************************************************

subroutine remd_setup(numexchg)

  use file_io_mod,       only : amopen
  use gbl_constants_mod
  use mdin_ctrl_dat_mod, only : ig, temp0, isgld
  use pmemd_lib_mod,     only : mexit
  use reservoir_mod,     only : load_reservoir_files 

  implicit none

  ! Passed variables

  integer, intent(in) :: numexchg

  ! Local variables

  integer             :: alloc_failed
  integer, dimension(mpi_status_size) :: stat_array ! needed for mpi_recv

  ! Force an initial re-reading of the variables (doesn't hurt)

  remd_modwt = .true.

  ! Nullify communicator
  
  remd_comm = mpi_comm_null

  ! Set replica index and coordinate index.
  remd_repidx = master_rank
  remd_crdidx = master_rank ! TODO should also be set from restart 

  ! Set random number generator
  call amrset_gen(remd_randgen, ig)

  ! Open remlog file and write some initial info

  if (master_master .and. remd_method .ne. -1) then

    call amopen(remlog, remlog_name, owrite, 'F', 'W')
    write(remlog, '(a)')     '# Replica Exchange log file'
    write(remlog, '(a,i10)') '# numexchg is ', numexchg

    if (isgld == 1) then
     ! Write out RXSGLD filenames
      write(remlog,'(a)') "# RXSGLD(RXSGMD) filenames:"
      write(remlog,'(a,a)') "#   rxsgldlog= ",trim(remlog_name)
      write(remlog,'(a,a)') "#   rxsgldtype= ",trim(remtype_name)
    else
      write(remlog, '(a)')     '# REMD filenames:'
      write(remlog, '(a,a)')   '#   remlog= ', trim(remlog_name)
      write(remlog, '(a,a)')   '#   remtype= ', trim(remtype_name)
      if (remd_method .ne. 3) then 
      if (rremd_type .gt. 0) then
        write(remlog, '(a, i5)') '#   RREMD type = ', rremd_type
        write(remlog, '(a, a)') '#   reservoir = ', trim(reservoir_name(1))
        if (rremd_type.eq.4 .or. rremd_type.eq.6) & 
          write(remlog, '(a, a)') '#   reservoir = ', trim(reservoir_name(1))
        if (rremd_type.eq.5 .or. rremd_type.eq.6) & 
          write(remlog, '(a, a)') '#   reservoir2 = ', trim(reservoir_name(2))
      endif
      
        ! Not Implemented yet
      !  if (rremd_type .eq. 3) then
      !    if (cluster_info_file_name .ne. 'cluster.info') then
      !      write(remlog, '(a,a)') '#   clusterinfo = ', trim(cluster_info_file_name)
      !    end if
      !  end if
      end if
    endif
  end if

  ! If a REMD exchange definition file was specified, open that and set up
  ! REMD vars. Otherwise, define them the way they'd be defined in sander.

  if (remd_method .eq. -1) then
    call set_dimensions_from_file(numexchg)
  else

    ! Allocate remd_types with dimension 1
    allocate(remd_types(1), &
             stat=alloc_failed)

    if (alloc_failed .ne. 0) &
      call alloc_error('remd_setup', 'remd_types not allocated')
      
    remd_types(1) = remd_method
    
    if (remd_method .gt. 0) then
      if (master_rank.eq.0) then
        call mpi_send(reservoir_name, max_fn_len*2, mpi_char, &
            (numgroups-1), remd_tag, pmemd_master_comm, err_code_mpi)
      elseif (master_rank.eq.(numgroups-1)) then
        call mpi_recv(reservoir_name, max_fn_len*2, mpi_char, &
            0, remd_tag, pmemd_master_comm, stat_array, err_code_mpi)
      endif
      endif
    
    if (remd_method .eq. 1) then
       call temp_remd_setup
       if(rremd_type .gt. 0) then
          call load_reservoir_files(rremd_type)
       end if
    else if (remd_method .eq. 3) then
      call h_remd_setup
       if (rremd_type .ge. 4 .and. rremd_type .le. 6) then
          call load_reservoir_files(rremd_type)
       end if      
    else if (remd_method .eq. 4) then
      call ph_remd_setup
    else if (remd_method .eq. 5) then
      call e_remd_setup
    end if

  end if

end subroutine remd_setup

!*******************************************************************************
!
! Subroutine: slave_remd_setup
!
! Description: Performs necessary setup for slave nodes
!
!*******************************************************************************

subroutine slave_remd_setup
  
  use prmtop_dat_mod, only : natom

  implicit none

  integer :: alloc_failed

  ! Slaves do not participate in inter-replica communications

  remd_comm = mpi_comm_null
  group_master_comm = mpi_comm_null

  if (remd_method .eq. 3) then

    allocate(crd_temp(3, natom), frc_temp(3, natom), stat=alloc_failed)

    if (alloc_failed .ne. 0) &
      call alloc_error('slave_remd_setup', 'Error allocating temporary arrays')

  else if (remd_method .eq. 1 .and. hybridgb .gt. 0) then

    allocate(hybridsolvent_remd_crd(3, natom), hybridsolvent_remd_vel(3, natom), &
             hybridsolvent_remd_frc(3, natom), stat=alloc_failed)

    if (alloc_failed .ne. 0) &
      call alloc_error('slave_remd_setup', 'Error allocating temporary arrays')

  else if (remd_method .eq. -1) then

    ! Receive the remd_dimension
    call mpi_bcast(remd_dimension, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)

    allocate(crd_temp(3, natom), &
             frc_temp(3, natom), &
             remd_types(remd_dimension), &
             stat=alloc_failed)
    
    if (alloc_failed .ne. 0) &
      call alloc_error('slave_remd_setup', 'Error allocating temporary arrays')

    ! Receive the remd_types
    call mpi_bcast(remd_types, remd_dimension, mpi_integer, 0, &
                   pmemd_comm, err_code_mpi)

  end if

end subroutine slave_remd_setup

!*******************************************************************************
!
! Subroutine: temp_remd_setup
!
! Description: This sets up the REMD data structures for remd_method .eq. 1
!
!*******************************************************************************

subroutine temp_remd_setup

  use gbl_constants_mod, only : error_hdr
  use pmemd_lib_mod,     only : mexit
  use mdin_ctrl_dat_mod, only : isgld

  implicit none

! Local variables

  double precision   :: my_temp       ! temperature of my replica

  integer            :: alloc_failed
  integer            :: i


  ! Because we are only exchanging in 1 dimension, our remd_comm and remd_rank
  ! are just our pmemd_master_comm and master_rank

  remd_dimension = 1
  call mpi_comm_dup(pmemd_master_comm, remd_comm, err_code_mpi)
  call mpi_comm_rank(remd_comm, remd_rank, err_code_mpi)
  remd_master = remd_rank .eq. 0
  group_master_comm = mpi_comm_null
  group_master_rank = 0

  ! Allocate temperature table and exchange successes, and initialize them

  allocate( exchange_successes(remd_dimension, numgroups), &
            temperatures(numgroups),                       &
            control_exchange(remd_dimension),                  &
            even_exchanges(remd_dimension),                &
            group_num(remd_dimension),                     &
            replica_indexes(remd_dimension),               &
            index_list(numgroups),                         &
            stat = alloc_failed)

  if (alloc_failed .ne. 0) &
    call alloc_error('temp_remd_setup', 'allocation error')

! allocate memory for self-guiding temperature table: sgfttable
   if (isgld .eq. 1) then
      allocate( tempsgtable(numgroups), sgfttable(numgroups), sgfftable(numgroups), sgfgtable(numgroups), stat=alloc_failed)
      if (alloc_failed .ne. 0) &
       call alloc_error('temp_remd_setup', 'SGLD allocation error')
   end if

  exchange_successes(:,:) = 0
  even_exchanges(:) = .false.
  group_num(1) = 1
  replica_indexes(1) = 1 ! reset in collect_t_statevars

  ! Set up the temperature table

  call collect_t_statevars(1, numgroups)

  ! Write some T-REMD specific rem.log info

  if (master_master) then
    if (isgld==1) then
       write(remlog,'(a)')"# RXSGLD setup: stagid  temp0     sgft      sgff "
       do i=1,numgroups
          write(remlog,'(a,i4,f10.2,3f10.4)')"#      replica: ", &
                  i, temperatures(i),  sgfttable(i),  sgfftable(i)
       enddo
       write(remlog,'(a)') &
         "# Rep Stagid Vscale SGscale Temp    Templf Ep Acceptance(i,i+1)"
    else
       write(remlog, '(a)') '# Rep#, Velocity Scaling, T, Eptot, Temp0, NewTemp0, &
                         &Success rate (i,i+1), ResStruct#'
    endif
  endif

  my_temp = temperatures(master_rank + 1)

    if (isgld/=1) then
  ! Check to make sure we have no duplicate temperatures. Then set the partners
  ! array via a call to set_temp_partners

    do i = 1, numgroups
      if (temperatures(i) .eq. my_temp .and. i .ne. master_rank+1) then
        write(mdout, '(2a)') error_hdr, &
            'two temperatures are identical in temp_remd_setup'
        call mexit(mdout, 1)
      end if
    end do
  endif

  call set_partners(1, numgroups, remd_random_partner)

  return

end subroutine temp_remd_setup

!*******************************************************************************
!
! Subroutine: h_remd_setup
!
! Description: Sets up Hamiltonian-REMD in which coordinates are traded between
!              replicas.
!
!*******************************************************************************

subroutine h_remd_setup

  use prmtop_dat_mod,    only : natom

  implicit none

  integer :: alloc_failed

  ! Because we are only exchanging in H-space, our remd_comm and remd_rank are
  ! just our pmemd_master_comm and master_rank

  remd_dimension = 1
  call mpi_comm_dup(pmemd_master_comm, remd_comm, err_code_mpi)
  remd_rank = master_rank
  remd_master = master_master
  group_master_comm = mpi_comm_null
  group_master_rank = 0

  ! First allocate necessary data structures and initialize them

  allocate( exchange_successes(remd_dimension, numgroups), &
            control_exchange(remd_dimension),                  &
            even_exchanges(remd_dimension),                &
            crd_temp(3, natom),                            &
            frc_temp(3, natom),                            &
            group_num(1),                                  &
            replica_indexes(remd_dimension),               &
            total_left_fe(remd_dimension, numgroups, numgroups),  &
            total_right_fe(remd_dimension, numgroups, numgroups), &
            num_left_exchg(remd_dimension, numgroups, numgroups), &
            num_right_exchg(remd_dimension, numgroups, numgroups),&
            index_list(numgroups), &
            stat = alloc_failed)

  if (alloc_failed .ne. 0) call alloc_error('h_remd_setup', 'allocation error')

  ! Initialize counters

  num_right_exchg(:,:,:)  = 0
  num_left_exchg(:,:,:)   = 0
  total_left_fe(:,:,:)    = 0.d0
  total_right_fe(:,:,:)   = 0.d0
  exchange_successes(:,:) = 0
  even_exchanges(:)       = .false.
  group_num(1)            = 1
  replica_indexes(1)      = master_rank + 1 ! My replica index is rank + 1

  ! Write the header line for the rem.log

  if (master_master) &
    write(remlog, '(a)') '# Rep#, Neibr#, Temp0, PotE(x_1), PotE(x_2), left_fe,&
                        & right_fe, Success, Success rate (i,i+1)'

  ! Set partners here

  call set_partners(1, numgroups, remd_random_partner)

end subroutine h_remd_setup

!*******************************************************************************
!
! Subroutine: ph_remd_setup
!
! Description: This sets up the REMD data structures for remd_method .eq. 4
!
!*******************************************************************************

subroutine ph_remd_setup

  use gbl_constants_mod, only : error_hdr
  use pmemd_lib_mod,     only : mexit

  implicit none

! Local variables

  double precision   :: my_ph       ! ph of my replica

  integer            :: alloc_failed
  integer            :: i


  ! Because we are only exchanging in 1 dimension, our remd_comm and remd_rank
  ! are just our pmemd_master_comm and master_rank

  remd_dimension = 1
  call mpi_comm_dup(pmemd_master_comm, remd_comm, err_code_mpi)
  call mpi_comm_rank(remd_comm, remd_rank, err_code_mpi)
  remd_master = remd_rank .eq. 0
  group_master_comm = mpi_comm_null
  group_master_rank = 0

  ! Allocate ph table and exchange successes, and initialize them

  allocate( exchange_successes(remd_dimension, numgroups), &
            phs(numgroups),                                &
            control_exchange(remd_dimension),                  &
            even_exchanges(remd_dimension),                &
            group_num(remd_dimension),                     &
            replica_indexes(remd_dimension),               &
            index_list(numgroups),                         &
            stat = alloc_failed)

  if (alloc_failed .ne. 0) &
    call alloc_error('ph_remd_setup', 'allocation error')

  exchange_successes(:,:) = 0
  even_exchanges(:) = .false.
  group_num(1) = 1
  replica_indexes(1) = 1 ! reset in collect_ph_statevars

  ! Set up the ph table

  call collect_ph_statevars(1, numgroups)

  ! Write some pH-REMD specific rem.log info

  if (master_master) &
    write(remlog, '(a)') "# Rep#, N_prot, old_pH, new_pH, Success rate (i,i+1)"

  my_ph = phs(master_rank + 1)

  ! Check to make sure we have no duplicate phs. Then set the partners
  ! array via a call to set_temp_partners

  do i = 1, numgroups
    if (phs(i) .eq. my_ph .and. i .ne. master_rank+1) then
      write(mdout, '(2a)') error_hdr, &
            'two solution pHs are identical in ph_remd_setup'
      call mexit(mdout, 1)
    end if
  end do

  call set_partners(1, numgroups, remd_random_partner)

  return

end subroutine ph_remd_setup

!*******************************************************************************
!
! Subroutine: e_remd_setup
!
! Description: This sets up the REMD data structures for remd_method .eq. 5
!
!*******************************************************************************

subroutine e_remd_setup

  use gbl_constants_mod, only : error_hdr
  use pmemd_lib_mod,     only : mexit

  implicit none

! Local variables

  double precision   :: my_e       ! redox potential of my replica

  integer            :: alloc_failed
  integer            :: i


  ! Because we are only exchanging in 1 dimension, our remd_comm and remd_rank
  ! are just our pmemd_master_comm and master_rank

  remd_dimension = 1
  call mpi_comm_dup(pmemd_master_comm, remd_comm, err_code_mpi)
  call mpi_comm_rank(remd_comm, remd_rank, err_code_mpi)
  remd_master = remd_rank .eq. 0
  group_master_comm = mpi_comm_null
  group_master_rank = 0

  ! Allocate redox potential table and exchange successes, and initialize them

  allocate( exchange_successes(remd_dimension, numgroups), &
            redoxs(numgroups),                             &
            control_exchange(remd_dimension),                  &
            even_exchanges(remd_dimension),                &
            group_num(remd_dimension),                     &
            replica_indexes(remd_dimension),               &
            index_list(numgroups),                         &
            stat = alloc_failed)

  if (alloc_failed .ne. 0) &
    call alloc_error('e_remd_setup', 'allocation error')

  exchange_successes(:,:) = 0
  even_exchanges(:) = .false.
  group_num(1) = 1
  replica_indexes(1) = 1 ! reset in collect_e_statevars

  ! Set up the redox table

  call collect_e_statevars(1, numgroups)

  ! Write some E-REMD specific rem.log info

  if (master_master) &
    write(remlog, '(a)') "# Rep#, N_elec, old_E, new_E, Success rate (i,i+1)"

  my_e = redoxs(master_rank + 1)

  ! Check to make sure we have no duplicate redox potentials. Then set the partners
  ! array via a call to set_temp_partners

  do i = 1, numgroups
    if (redoxs(i) .eq. my_e .and. i .ne. master_rank+1) then
      write(mdout, '(2a)') error_hdr, &
            'two solution redox potentials are identical in e_remd_setup'
      call mexit(mdout, 1)
    end if
  end do

  call set_partners(1, numgroups, remd_random_partner)

  return

end subroutine e_remd_setup

!*******************************************************************************
!
! Subroutine: alloc_error
!
! Description: Called in the case of an allocation error
!
!*******************************************************************************

subroutine alloc_error(routine, message)

  use gbl_constants_mod, only : error_hdr, extra_line_hdr
  use pmemd_lib_mod,     only : mexit

  implicit none

  ! Passed variables

  character(*), intent(in) :: routine, message

  write(mdout, '(a,a,a)') error_hdr, 'Error in ', routine
  write(mdout, '(a,a)') extra_line_hdr, message

  call mexit(mdout, 1)

end subroutine alloc_error

!*******************************************************************************
!
! Subroutine: set_dimensions_from_file
!
! Description: Parses remd_dimension file to define the remd dimensions/ranks
!              and performs other basic setup required for each dimension
!              NOTE: Every GOTO in this function jumps to a specific error
!              printout followed by a call mexit(mdout, 1), which matches to
!              what happens following a read() error
!
!*******************************************************************************

subroutine set_dimensions_from_file(numexchg)

  use file_io_mod, only       : amopen, nmlsrc
  use gbl_constants_mod, only : error_hdr, extra_line_hdr
  use pmemd_lib_mod, only     : mexit, upper, strip
  use prmtop_dat_mod, only    : natom
  use AmberNetcdf_mod, only   : NC_readRestartIndices
  use mdin_ctrl_dat_mod, only : irest, temp0, solvph, solve

  implicit none

! Parameter -- group size. This is because we can't have allocatables
! in namelists
  integer, parameter            :: GRPS = 1000

! Passed variables

  integer, intent(in) :: numexchg

! Local variables
  
  ! &multirem namelist
  integer, dimension(GRPS,GRPS) :: group
  character(80)                 :: exch_type, desc
  integer                       :: repidxIn, crdidxIn

  ! Utility variables (counters, error markers, etc.)
  integer                   :: i, j, idx
  integer                   :: ifind
  integer                   :: alloc_failed
  integer                   :: group_counter
  integer, allocatable      :: replica_assignments(:,:)
  integer, allocatable      :: replica_indexes_buf(:,:)
  character(len=5)          :: extension
  character(len=max_fn_len) :: filename
  character(len=80)         :: buf
  character(len=80),allocatable :: replica_desc(:)

  namelist / multirem /        group, exch_type, desc

! Variable descriptions:
!
!  group               : This is a set of integers that denote the replica #s
!                        involved in that set of communicating replicas
!  exch_type           : The type of exchange (TEMPERATURE)
!  temptype            : Token to declare an exchange in Temperature space
!  hamtype             : Token to declare an exchange in Hamiltonian space
!  i,j, group_counter  : counters
!  multirem            : namelist to extract information from
!  ifind               : for nmlsrc -- did we find our namelist? 0 = no, 1 = yes
!  alloc_failed        : did memory allocation fail?
!  replica_assignments : buffer for master_master when reading which group
!                        each replica belongs to in each dimension
!  replica_indexes_buf : buffer for reading replica ranks by master_master
!  extension           : file name extension (dimension #) for rem.log
!  filename            : full rem.log file name with ".extension" (above)
!  replica_desc        : The description string for each replica dimension


  if (master_master) then

    call amopen(remd_file, remd_dimension_name, 'O', 'F', 'R')
  
    ! Peek at the remd_file to see how many &multiterm namelists there are, and
    ! that will be how many dimensions we're exchanging in. We need this to
    ! allocate some of the necessary data structures
    remd_dimension = 0
    do
      call nmlsrc('multirem', remd_file, ifind)
      if (ifind .eq. 0) exit
      ! nmlsrc backs up to the start of the namelist, so we have to eat this
      ! line if we want nmlsrc to go to the next &multirem namelist
      read(remd_file, '(a80)') buf
      remd_dimension = remd_dimension + 1
    end do

    if (remd_dimension .eq. 0) go to 665 ! we need at least 1 dimension

    ! Allocate our necessary data structures
    allocate(remd_types(remd_dimension), &
             replica_assignments(remd_dimension, numgroups), &
             replica_indexes_buf(remd_dimension, numgroups), &
             replica_desc(remd_dimension), &
             stat=alloc_failed)

    if (alloc_failed .ne. 0) &
      call alloc_error('set_dimensions_from_file', 'remd_types not allocated')

    replica_assignments(:,:) = 0
    replica_indexes_buf(:,:) = 0

    ! Loop through all of the remd dimensions. At this point, remd_file has been
    ! rewound by subroutine nmlsrc
    do i = 1, remd_dimension

      call nmlsrc('multirem', remd_file, ifind)
     
      ! This should never happen because of the preliminary &multirem counting
      if (ifind .eq. 0) then
        write(mdout, '(2a)') error_hdr, 'set_dimensions_from_file: &
                        &Should not be here...'
        call mexit(mdout, 1)
      end if

      ! Initialize our namelist variables
      group(:,:) = 0
      exch_type = ' '
      desc = ' '

      ! Read that namelist
      read(remd_file, nml=multirem, err=666)

      ! Make the exch_type upper-case to make it case-insensitive
      call upper(exch_type)

      ! Store our replica description. It's just cosmetic. A form of self-
      ! documentation for the user which will also be dumped to the rem.log
      ! file for users' benefit. It has no impact on the simulation
      replica_desc(i) = desc

      ! Assign the remd_type
      select case(trim(exch_type))
        case ('TEMPERATURE')
          remd_types(i) = 1
        case ('TEMP')
          remd_types(i) = 1
        case ('HAMILTONIAN')
          remd_types(i) = 3
        case ('HREMD')
          remd_types(i) = 3
        case ('PH')
          remd_types(i) = 4
        case ('REDOX')
          remd_types(i) = 5
        case default
          write(mdout, '(3a)') error_hdr, 'Unrecognized EXCH_TYPE ', &
                               trim(exch_type)
          call mexit(mdout, 1)
      end select
      
      do group_counter = 1, GRPS

        ! jump out of the loop if this group is not assigned. That's
        ! the last group in this dimension
        if (group(group_counter, 1) .eq. 0) exit

        ! Assign all of the replicas in this group
        do j = 1, GRPS

          idx = group(group_counter, j)

          ! bail out if this is 0 -- that means we've reached the last replica
          ! in this group (make sure we have an even number of replicas
          if (idx .eq. 0) then
            if (mod(j, 2) .eq. 0) go to 669
            exit
          end if

          ! Catch bad replica assignment
          if (idx .le. 0 .or. idx .gt. numgroups) go to 671

          ! Assign the replica, but make sure each replica is only assigned once
          if (replica_assignments(i, idx) .eq. 0) then
            replica_assignments(i, idx) = group_counter
            replica_indexes_buf(i, idx) = j
          else
            go to 670 ! duplicate assignment
          end if

        end do

      end do ! group_counter = 1, GRPS

      ! Make sure everyone was assigned here
      do j = 1, numgroups
        if (replica_assignments(i,j) .eq. 0) go to 670
      end do

    end do ! i = 1, num_dim

    close(remd_file)

  end if ! (master_master)

  group_master_comm = mpi_comm_null

  ! Broadcast the remd_dimension along pmemd_master_comm
  call mpi_bcast(remd_dimension, 1, mpi_integer, 0, &
                 pmemd_master_comm, err_code_mpi)

  call mpi_barrier(pmemd_master_comm, err_code_mpi)

  ! Everybody really needs remd_dimension to allocate all of the necessary
  ! arrays. We don't use MPI_COMM_WORLD above since doing so would make the
  ! (possibly erroneous) assumption that the master_master thread is also the
  ! worldmaster thread...

  call mpi_bcast(remd_dimension, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)

  ! Now the non-master_master's need to allocate remd_types

  if (.not. master_master) then
    allocate(remd_types(remd_dimension), &
             replica_assignments(remd_dimension, numgroups), &
             replica_indexes_buf(remd_dimension, numgroups), &
             stat=alloc_failed)
    if (alloc_failed .ne. 0) &
      call alloc_error('set_dimensions_from_file', 'remd_types not allocated')
  end if

  ! Now allocate all data structures that everybody here needs
  
  allocate(group_num(remd_dimension),       &
           replica_indexes(remd_dimension), &
           exchange_successes(remd_dimension, numgroups), &
           control_exchange(remd_dimension),    &
           even_exchanges(remd_dimension),  &
           index_list(numgroups),           &
           temperatures(numgroups),         &
           phs(numgroups),                  &
           redoxs(numgroups),               &
           multid_print_data(numgroups),    &
           multid_print_data_buf(numgroups, numgroups), &
           crd_temp(3, natom), &
           frc_temp(3, natom), &
           num_right_exchg(remd_dimension, numgroups, numgroups),&
           num_left_exchg(remd_dimension, numgroups, numgroups), &
           total_left_fe(remd_dimension, numgroups, numgroups),  &
           total_right_fe(remd_dimension, numgroups, numgroups), &
           stat=alloc_failed)

  if (alloc_failed .ne. 0) &
    call alloc_error('set_dimensions_from_file', 'main data allocation failed')

  ! Broadcast the remd_types

  call mpi_bcast(remd_types, remd_dimension, mpi_integer, 0, &
                 pmemd_master_comm, err_code_mpi)

  ! Broadcast the remd_types to the rest of the pmemd_comm, which has to be
  ! received in the slave_setup routine just like for remd_dimension

  call mpi_bcast(remd_types, remd_dimension, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)
  
  ! Zero out a bunch of our data arrays so they aren't filled with junk

  multid_print_data(:)       = NULL_REMLOG_DATA
  multid_print_data_buf(:,:) = NULL_REMLOG_DATA
  exchange_successes(:,:)    = 0
  even_exchanges(:)          = .true.
  num_right_exchg(:,:,:)     = 0
  num_left_exchg(:,:,:)      = 0
  total_left_fe(:,:,:)       = 0.d0
  total_right_fe(:,:,:)      = 0.d0
  
  ! Now set up all data structures for each dimension

  do i = 1, remd_dimension

    ! First scatter our placement in the group/replica ladders in each dimension

    call mpi_scatter(replica_assignments(i,:), 1, mpi_integer, &
                     group_num(i), 1, mpi_integer,             &
                     0, pmemd_master_comm, err_code_mpi)
    call mpi_scatter(replica_indexes_buf(i,:), 1, mpi_integer, &
                     replica_indexes(i), 1, mpi_integer,       &
                     0, pmemd_master_comm, err_code_mpi)

  end do ! i = 1, remd_dimension

  if ( irest .eq. 1 ) then
    alloc_failed = NC_readRestartIndices( inpcrd_name, replica_indexes, &
                                          group_num, repidxIn, crdidxIn,&
                                          remd_dimension )
    if (repidxIn .gt. -1 .and. crdidxIn .gt. -1) then
      remd_repidx = repidxIn
      remd_crdidx = crdidxIn
      write(mdout,'(2(a,i6))') '| Overall replica indices from restart: RepIdx=',&
                               remd_repidx, ' CrdIdx=', remd_crdidx
    else
      write(6,'(2(a,i6))') '| Initial overall replica indices: RepIdx=',&
                           remd_repidx, ' CrdIdx=', remd_crdidx
    endif
    if ( alloc_failed .ne. 0 ) then ! TODO: mexit if -1?
      write(mdout,'(a)') '| Warning: Replica indices will NOT be used to &
                         &restart MultiD-REMD run.'
    else
      write(mdout,'(a)') '| Restarting MultiD-REMD run. This replica will use indices:'
      write(mdout,'(a,13i6)') '| ', (replica_indexes(i),i=1,remd_dimension)
    endif

  endif

  do i = 1, remd_dimension
    
    ! Now we set up our remd_comm and get our respective sizes

    if (remd_comm .ne. mpi_comm_null) &
      call mpi_comm_free(remd_comm, err_code_mpi)
    remd_comm = mpi_comm_null
    call mpi_barrier(pmemd_master_comm, err_code_mpi)
    call mpi_comm_split(pmemd_master_comm, group_num(i), &
                        master_rank, remd_comm, err_code_mpi)
    call mpi_comm_size(remd_comm, remd_size, err_code_mpi)
    call mpi_comm_rank(remd_comm, remd_rank, err_code_mpi)

    ! Here we do any additional setup that's needed for any of the exchange
    ! types. For instance, we have to collect_t_statevars for T-REMD in order to set
    ! up the temperature ladder. H-REMD doesn't require any additional setup.

    select case(remd_types(i))
      case(1) ! TEMPERATURE
        call collect_t_statevars(i, remd_size)
      case(4) ! PH
        call collect_ph_statevars(i, remd_size)
      case(5) ! REDOX
        call collect_e_statevars(i, remd_size)
    end select

    ! Now we set up the rem.log files. Each dimension will get its own rem.log
    ! taking the -remlog <filename> and appending .1, .2, .3, ... etc. to it
    ! based on the dimension
    
    if (master_master) then
      
      write(extension, '(i5)') i
      call strip(extension)

      if (len_trim(remlog_name) + len_trim(extension) + 1 .gt. max_fn_len) then
        write(mdout, '(2a)') error_hdr, 'rem.log file name overflow. Increase &
          &max_fn_len in file_io_dat.F90 and recompile'
        call mexit(mdout, 1)
      end if

      filename = trim(remlog_name) // '.' // trim(extension)
      call amopen(remlog, filename, owrite, 'F', 'W')

      write(remlog,'(a)')       '# Replica Exchange log file'
      write(remlog,'(a,i10)')   '# numexchg is ', numexchg
      write(remlog,'(2(a,i4))') '# Dimension ', i, ' of ', remd_dimension
      if (len_trim(replica_desc(i)) .gt. 0) &
         write(remlog,'(2a)')      '# Description: ', trim(replica_desc(i))

      ! Write the dimension type here

      select case (remd_types(i))
        case(1)
          write(remlog, '(a)') '# exchange_type = TEMPERATURE'
          write(remlog, '(a)') '# REMD filenames:'
          write(remlog, '(2a)') '# remlog= ', trim(filename)
          write(remlog, '(2a)') '# remd dimension file= ', &
                                trim(remd_dimension_name)
          write(remlog, '(a)') '# Rep#, Velocity Scaling, T, Eptot, Temp0, &
                               &NewTemp0, Success rate (i,i+1), ResStruct#'
        case(3)
          write(remlog, '(a)') '# exchange_type = HAMILTONIAN'
          write(remlog, '(a)') '# REMD filenames:'
          write(remlog, '(2a)') '# remlog= ', trim(filename)
          write(remlog, '(2a)') '# remd dimension file= ', &
                                trim(remd_dimension_name)
          write(remlog, '(a)') '# Rep#, Neibr#, Temp0, PotE(x_1), PotE(x_2), &
                               &left_fe, right_fe, Success, Success rate (i,i+1)'
        case(4)
          write(remlog, '(a)') '# exchange_type = PH'
          write(remlog, '(a)') '# REMD filenames:'
          write(remlog, '(2a)') '# remlog= ', trim(filename)
          write(remlog, '(2a)') '# remd dimension file= ', &
                                trim(remd_dimension_name)
          write(remlog, '(a)') '# Rep#, N_prot, old_pH, new_pH, &
                                &Success rate (i,i+1)'
        case(5)
          write(remlog, '(a)') '# exchange_type = REDOX'
          write(remlog, '(a)') '# REMD filenames:'
          write(remlog, '(2a)') '# remlog= ', trim(filename)
          write(remlog, '(2a)') '# remd dimension file= ', &
                                trim(remd_dimension_name)
          write(remlog, '(a)') '# Rep#, N_elec, old_E, new_E, &
                                &Success rate (i,i+1)'
      end select

      close(remlog)

    end if ! master_master

  end do ! i = 1, remd_dimension

  ! Deallocate the temporary allocatable arrays
  if (allocated(replica_assignments)) deallocate(replica_assignments)
  if (allocated(replica_indexes_buf)) deallocate(replica_indexes_buf)
  if (allocated(replica_desc)) deallocate(replica_desc)

  return

! Reading errors (differentiate to be helpful)

665 write(mdout, '(4a)') error_hdr, 'Could not find &multirem in ', &
                         trim(remd_dimension_name), '!'
    call mexit(mdout, 1)

666 write(mdout, '(2a,i3,a)') error_hdr, 'Could not read the ', i, &
                              'th &multirem namelist!'
    call mexit(mdout, 1)

667 write(mdout, '(2a)') error_hdr, 'Every &multirem namelist needs EXCH_TYPE!'
    call mexit(mdout, 1)

668 write(mdout, '(4a)') error_hdr, &
            'Bad or incomplete replica assignment! Every replica needs ', &
            extra_line_hdr, 'to be assigned in each dimension!'
    call mexit(mdout, 1)

669 write(mdout, '(2a)') error_hdr, &
      'All replica groups need an even number of replicas.'
    call mexit(mdout, 1)

670 write(mdout, '(2a)') error_hdr, &
     'Each replica must be assigned to one and only one group in each dimension'
    call mexit(mdout, 1)

671 write(mdout, '(2a,i4)') error_hdr, &
     'Bad replica assignment. Replicas must be between 1 and ', numgroups
    call mexit(mdout, 1)

    return ! this will never be reached

end subroutine set_dimensions_from_file

!*******************************************************************************
!
! Subroutine: set_partners
!
! Description: Sets the partners array in my dimension based on each replica's
!              rank in this given dimension.
!
!*******************************************************************************

subroutine set_partners(rem_dim, num_replicas, rem_random_partner)

  use pmemd_lib_mod, only : mexit

  implicit none

! Passed Variables

  integer, intent(in) :: rem_dim      ! temperature dimension
  integer, intent(in) :: num_replicas ! # of replicas in this dimension
  integer, intent(in) :: rem_random_partner ! # Decide to pick a random partner or not. (Any value different than zero will pick random)

! Local variables

  integer             :: lower_neibr  ! what the lower neighbor rank will be
  integer             :: higher_neibr ! what the higher neighbor rank will be
  integer             :: my_idx       ! my index in this dimension
  integer             :: partner_list(num_replicas) ! list of partners for each replica, if picked randomly

  integer             :: i ! counter
  integer             :: ip ! partner index
  integer             :: irandval
  double precision    :: random_value ! random number generated

  my_idx = replica_indexes(rem_dim)

  if (my_idx .eq. 1) then
    lower_neibr = num_replicas
    higher_neibr = 2
  else if (my_idx .eq. num_replicas) then
    lower_neibr = my_idx - 1
    higher_neibr = 1
  else
    lower_neibr = my_idx - 1
    higher_neibr = my_idx + 1
  end if

  ! Find out everyone's rank in this REMD dimension, then search for our
  ! neighbors. NOTE that when we find our partners(1/2), they will be indexed
  ! starting from 1, not 0

  call mpi_allgather(my_idx, 1, mpi_integer, index_list, 1, &
                     mpi_integer, remd_comm, err_code_mpi)

  ! Find random partners
  if (master_master .and. rem_random_partner .ne. 0) then
    partner_list(:) = -1
    do i = 1, num_replicas
      if (count(partner_list==-1) .eq. 1) cycle
      if (partner_list(i) .ne. -1) cycle
      irandval = -1
      do while (irandval .eq. -1 .or. irandval .eq. index_list(i) .or. find_location(partner_list, num_replicas, irandval) .ne. 0)
        call amrand_gen(remd_randgen, random_value)
        irandval = nint(random_value*(num_replicas-1)+1)
      end do
      partner_list(i) = irandval
      ! Change partner info to match what was decided here
      ip = find_location(index_list, num_replicas, irandval)
      partner_list(ip) = index_list(i)
    end do
  end if

  call mpi_bcast(partner_list, num_replicas, mpi_integer, 0, remd_comm, err_code_mpi)

  partners(:) = -1
  if (rem_random_partner .eq. 0) then
    do i = 1, num_replicas
      if (index_list(i) .eq. lower_neibr)  partners(1) = i
      if (index_list(i) .eq. higher_neibr) partners(2) = i
    end do
  else
    partners(1) = find_location(index_list, num_replicas, partner_list(find_location(index_list, num_replicas, my_idx)))
    partners(2) = partners(1)
  end if

  if (partners(1) .eq. -1 .or. partners(2) .eq. -1) then
    write(0, '(a)') 'Creation of partners array failed! Could not find partners'
    write(0, '(a)') 'Check that the numbering of replicas in the replica dimension file is'
    write(0, '(a)') 'consistent with the total number of replicas (and if present, that'
    write(0, '(a)') 'the replica indices in the input coordinates are as well).'
    do i = 1, num_replicas
      write(6, '(2(a,i6))') '    index_list(', i, ')= ', index_list(i)
    end do
    call mexit(mdout, 1)
  end if

  if (gremd_acyc .gt. 0) then
    if (my_idx .eq. 1) partners(1)=-1
    if (my_idx .eq. num_replicas) partners(2)=-1
  end if
  ! Determine if we are an even replica here. Note that an even replica means we
  ! number beginning from zero. However, my_idx begins from 1, so our "even
  ! replicas" here are actually odd in my_idx
  ! VWDC: control_exchange is also used to decide who will control the replica exchange attempt,
  !       so in case of random partners, the decision of who controls is based on which replica
  !       number is smaller
  if (rem_random_partner .eq. 0) then
    control_exchange(rem_dim) = mod(my_idx, 2) .eq. 1
  else
    control_exchange(rem_dim) = my_idx .lt. partners(1)
  end if

  return

end subroutine set_partners

!*******************************************************************************
!
! Subroutine: rescale_velocities
!
! Description: rescale the velocities for the new temperature
!
!*******************************************************************************

subroutine rescale_velocities(atm_cnt, vel, vel_scale_factor)

  use mdin_ctrl_dat_mod, only : temp0

  implicit none

! Passed variables

  integer, intent(in)             :: atm_cnt
  double precision, intent(inout) :: vel(3, atm_cnt)
  double precision, intent(in)    :: vel_scale_factor

#ifdef VERBOSE_REMD
  if (master) &
    write(mdout, '(a,f8.3,a,f8.3)') '| REMD: scaling velocities by ', &
      vel_scale_factor, ' to match new bath T ', temp0
#endif

#ifdef CUDA
  call gpu_scale_velocities(vel_scale_factor)
#else
  vel(:, :) = vel(:, :) * vel_scale_factor
#endif

end subroutine rescale_velocities

!*******************************************************************************
!
! Subroutine: bcast_remd_method
!
! Description: Broadcasts the replica exchange method we're using so everyone in
!              the world knows what it is
!
!*******************************************************************************

subroutine bcast_remd_method

  implicit none

  call mpi_bcast(remd_method, 1, mpi_integer, 0, mpi_comm_world, err_code_mpi)
  call mpi_bcast(rremd_type, 1, mpi_integer, 0, mpi_comm_world, err_code_mpi)
  call mpi_bcast(remd_random_partner, 1, mpi_integer, 0, mpi_comm_world, err_code_mpi)

  return

end subroutine bcast_remd_method

!*******************************************************************************
!
! Subroutine: collect_t_statevars
!
! Description: Collects new temperatures in case they were changed by NMR
!              restraints
!
!*******************************************************************************

subroutine collect_t_statevars(t_dim, num_replicas)

  use gbl_constants_mod, only : error_hdr
  use mdin_ctrl_dat_mod, only : temp0,isgld,tempsg,sgft,sgff,sgfg
  use pmemd_lib_mod, only     : mexit
  !use sgld_mod, only: tsgld,sgfti,sgffi,sgfgi

  implicit none

! Passed variables

  integer, intent(in) :: t_dim        ! temperature dimension
  integer, intent(in) :: num_replicas ! number of replicas in this dimension

! Local variables

  integer             :: i ! counter
  double precision    :: my_temp
  double precision    :: param ! thermodynamic parameter
  double precision    :: my_tempsg,my_sgft,my_sgff,my_sgfg

  if (.not. master) return

  param = temp0

  temperatures(:) = 0.d0

  call mpi_allgather(   param,     1, mpi_double_precision, &
                     temperatures, 1, mpi_double_precision, &
                     remd_comm, err_code_mpi)
  if(isgld==1)then
    tempsgtable(:) = 0.d0
    sgfttable(:) = 0.d0
    sgfftable(:) = 0.d0
    sgfgtable(:) = 0.d0
    call mpi_allgather(tempsg,     1, mpi_double_precision, &
                     tempsgtable, 1, mpi_double_precision, &
                     remd_comm, err_code_mpi)
    call mpi_allgather(sgft,     1, mpi_double_precision, &
                     sgfttable, 1, mpi_double_precision, &
                     remd_comm, err_code_mpi)
    call mpi_allgather(sgff,     1, mpi_double_precision, &
                     sgfftable, 1, mpi_double_precision, &
                     remd_comm, err_code_mpi)
    call mpi_allgather(sgfg,     1, mpi_double_precision, &
                     sgfgtable, 1, mpi_double_precision, &
                     remd_comm, err_code_mpi)
    my_tempsg = tempsgtable(remd_rank + 1)
    my_sgft = sgfttable(remd_rank + 1)
    my_sgff = sgfftable(remd_rank + 1)
    my_sgfg = sgfgtable(remd_rank + 1)
  endif
  ! Now determine your overall rank in the temperature ladder

  replica_indexes(t_dim) = 0
  my_temp = temperatures(remd_rank + 1)
  
  do i = 1, num_replicas
    if (i .eq. remd_rank + 1) continue
    if (temperatures(i) .lt. my_temp) &
        replica_indexes(t_dim) = replica_indexes(t_dim) + 1
    if (temperatures(i) .eq. my_temp .and. i .ne. remd_rank + 1) then
    if(isgld==1)then
      if (tempsgtable(i) .lt. my_tempsg) &
        replica_indexes(t_dim) = replica_indexes(t_dim) + 1
      if (tempsgtable(i) .eq. my_tempsg.and.sgfttable(i) .lt. my_sgft) &
        replica_indexes(t_dim) = replica_indexes(t_dim) + 1
      if(tempsgtable(i) .eq. my_tempsg.and.sgfttable(i) .eq. my_sgft .and. sgfftable(i).gt. my_sgff)&
        replica_indexes(t_dim) = replica_indexes(t_dim) + 1
      if(tempsgtable(i) .eq. my_tempsg.and.sgfttable(i) .eq. my_sgft .and. sgfftable(i).eq. my_sgff &
      .and. i.le.remd_rank)replica_indexes(t_dim) = replica_indexes(t_dim) + 1
    else
      if(remd_method .eq. 1) then
         if (rremd_type .gt. 0 .and. hybridgb .le. 0) then
            write(mdout, '(2a)') error_hdr, 'Duplicate temperatures in Reservoir REMD'
         else if (hybridgb .gt. 0) then
            write(mdout, '(2a)') error_hdr, 'Duplicate temperatures in Hybrid Solvent REMD'
         else
            write(mdout, '(2a)') error_hdr, 'Duplicate temperatures in T-REMD'
         end if
      end if
      call mexit(mdout, 1)
    end if
    endif
  end do

  replica_indexes(t_dim) = replica_indexes(t_dim) + 1

  return

end subroutine collect_t_statevars

!*******************************************************************************
!
! Subroutine: collect_ph_statevars
!
! Description: Collects new pHs in case they were changed by NMR
!              restraints
!
!*******************************************************************************

subroutine collect_ph_statevars(ph_dim, num_replicas)

  use gbl_constants_mod, only : error_hdr
  use mdin_ctrl_dat_mod, only : solvph
  use pmemd_lib_mod, only     : mexit

  implicit none

! Passed variables

  integer, intent(in) :: ph_dim        ! ph dimension
  integer, intent(in) :: num_replicas ! number of replicas in this dimension

! Local variables

  integer             :: i ! counter
  double precision    :: my_ph
  double precision    :: param ! thermodynamic parameter

  if (.not. master) return

  param = solvph

  phs(:) = 0.d0

  call mpi_allgather(   param,     1, mpi_double_precision, &
                     phs, 1, mpi_double_precision, &
                     remd_comm, err_code_mpi)
  ! Now determine your overall rank in the ph ladder

  replica_indexes(ph_dim) = 0
  my_ph = phs(remd_rank + 1)
  
  do i = 1, num_replicas
    if (i .eq. remd_rank + 1) continue
    if (phs(i) .lt. my_ph) &
        replica_indexes(ph_dim) = replica_indexes(ph_dim) + 1
    if (phs(i) .eq. my_ph .and. i .ne. remd_rank + 1) then
      write(mdout, '(2a)') error_hdr, 'Duplicate pHs in pH-REMD'
      call mexit(mdout, 1)
    end if
  end do

  replica_indexes(ph_dim) = replica_indexes(ph_dim) + 1

  return

end subroutine collect_ph_statevars

!*******************************************************************************
!
! Subroutine: collect_e_statevars
!
! Description: Collects new redoxs in case they were changed by NMR
!              restraints
!
!*******************************************************************************

subroutine collect_e_statevars(e_dim, num_replicas)

  use gbl_constants_mod, only : error_hdr
  use mdin_ctrl_dat_mod, only : solve
  use pmemd_lib_mod, only     : mexit

  implicit none

! Passed variables

  integer, intent(in) :: e_dim        ! redox potential dimension
  integer, intent(in) :: num_replicas ! number of replicas in this dimension

! Local variables

  integer             :: i ! counter
  double precision    :: my_e
  double precision    :: param ! thermodynamic parameter

  if (.not. master) return

  param = solve

  redoxs(:) = 0.d0

  call mpi_allgather(   param,     1, mpi_double_precision, &
                     redoxs, 1, mpi_double_precision, &
                     remd_comm, err_code_mpi)
  ! Now determine your overall rank in the redox potential ladder

  replica_indexes(e_dim) = 0
  my_e = redoxs(remd_rank + 1)
  
  do i = 1, num_replicas
    if (i .eq. remd_rank + 1) continue
    if (redoxs(i) .lt. my_e) &
        replica_indexes(e_dim) = replica_indexes(e_dim) + 1
    if (redoxs(i) .eq. my_e .and. i .ne. remd_rank + 1) then
      write(mdout, '(2a)') error_hdr, 'Duplicate redox potentials in E-REMD'
      call mexit(mdout, 1)
    end if
  end do

  replica_indexes(e_dim) = replica_indexes(e_dim) + 1

  return

end subroutine collect_e_statevars

!*******************************************************************************
!
! Subroutine: remd_cleanup
!
! Description: Cleans up after a REMD run. It closes the REMD files and 
!              deallocates arrays
!
!*******************************************************************************

subroutine remd_cleanup
#ifdef BINTRAJ
  use AmberNetcdf_mod,    only: NC_close 
#endif
use reservoir_mod

implicit none

  integer::i


! Return if no REMD was run

  if (remd_method .eq. 0) return

! If REMD was run, deallocate any of the allocated data structures

  if (allocated(exchange_successes)   ) deallocate(exchange_successes)
  if (allocated(crd_temp)             ) deallocate(crd_temp)
  if (allocated(frc_temp)             ) deallocate(frc_temp)
  if (allocated(temperatures)         ) deallocate(temperatures)
  if (allocated(phs)                  ) deallocate(phs)
  if (allocated(redoxs)               ) deallocate(redoxs)
  if (allocated(control_exchange)         ) deallocate(control_exchange)
  if (allocated(even_exchanges)       ) deallocate(even_exchanges)
  if (allocated(group_num)            ) deallocate(group_num)
  if (allocated(replica_indexes)      ) deallocate(replica_indexes)
  if (allocated(multid_print_data)    ) deallocate(multid_print_data)
  if (allocated(multid_print_data_buf)) deallocate(multid_print_data_buf)
  if (allocated(total_left_fe)        ) deallocate(total_left_fe)
  if (allocated(total_right_fe)       ) deallocate(total_right_fe)
  if (allocated(num_left_exchg)       ) deallocate(num_left_exchg)
  if (allocated(num_right_exchg)      ) deallocate(num_right_exchg)
  if (allocated(index_list)           ) deallocate(index_list)
  if (allocated(remd_types)           ) deallocate(remd_types)
  if (allocated(rremd_crd)            ) deallocate(rremd_crd)
  if (allocated(rremd_vel)            ) deallocate(rremd_vel)
  if (allocated(reservoir_structure_energies_array)   ) &
                deallocate(reservoir_structure_energies_array)
  if (allocated(hybridsolvent_remd_crd) ) deallocate(hybridsolvent_remd_crd)
  if (allocated(hybridsolvent_remd_vel) ) deallocate(hybridsolvent_remd_vel)
  if (allocated(hybridsolvent_remd_frc) ) deallocate(hybridsolvent_remd_frc)

  if (rremd_type .gt. 0 .and. master_master) then
    do i=1,2
     if(reservoir_ncid(i) .ne. -1) then
#ifdef BINTRAJ
        call NC_close(reservoir_ncid(i))
#endif
     end if
    enddo
  
  end if

! Close any files that were opened
  
  if (remd_method .ne. -1 .and. master_master) close(remlog)

end subroutine remd_cleanup

!*******************************************************************************
!
! Function: find_location
!
! Description: Find the location of element in integer array
!
!*******************************************************************************
function find_location(iarray, iardim, ival)
  implicit none
  integer find_location
  integer, intent(in) :: iardim, ival
  integer, intent(in) :: iarray(iardim)
  integer i

  find_location = 0

  do i = 1, iardim
    if (iarray(i) .eq. ival) then
      find_location = i
      exit
    end if
end do
end function find_location

#endif /* MPI */

end module remd_mod
