#include "copyright.i"

!*******************************************************************************
!
! Module: multipmemd_mod
!
! Description: Contains the routines and such used to build the additional
!              communicators necessary for REMD and other multisander-like
!              applications
!              
!*******************************************************************************

module multipmemd_mod

#ifdef MPI

private

public setup_groups, &
       free_comms

!*******************************************************************************
!
! Subroutine: setup_groups
!
! Description: This subroutine splits up the global communicator into individual
!              replicas governed by pmemd_comm. For replicas that have to share
!              information, it sets up a communicator between the various pmemd
!              masters. It also does the parallel setup (moved from mstartup)
!
!*******************************************************************************

contains

subroutine setup_groups

  use file_io_dat_mod, only : grpfl_line_len, max_fn_len, &
                              groupfile_name, groupfile, &
                              groupline_buffer, proc_map_name, &
                              ng_nonsequential
  use file_io_mod,     only : amopen
  use get_cmdline_mod, only : get_cmdline
  use parallel_dat_mod
  use pmemd_lib_mod,   only : mexit
  use remd_mod,        only : remd_method, rremd_type, remd_random_partner

  implicit none
  
! Local variables

  ! store the lines of the groupfile

  character (len=grpfl_line_len), allocatable :: grouplines(:)
  character (len=grpfl_line_len), allocatable :: proc_map_lines(:)
  character (len=grpfl_line_len)              :: buf   ! buffer
  character (len=10)                          :: formt ! format

  integer :: i, j
  integer :: alloc_failed
  
  integer :: masterid     ! used to create pmemd_master_comm
!  integer :: master_size  ! size of comm_master, modified by FENG PAN

  logical :: version_requested ! see if we asked for a version number on the CL

! Get the world statistics here
  
#ifndef MPI_DEBUGGER
  call mpi_init(err_code_mpi)
#endif
  call mpi_comm_size(mpi_comm_world, worldsize, err_code_mpi)
  call mpi_comm_rank(mpi_comm_world, worldrank, err_code_mpi)

! By default, nullify all of the REMD-specific communicators so we don't choke
! trying to free them at the end

  remd_comm = mpi_comm_null
  group_master_comm = mpi_comm_null
  pmemd_master_comm = mpi_comm_null

! Get command line args here, to see if we need to set up groups. numgroups is
! initialized in get_cmdline.

  version_requested = .false.
  if (worldrank .eq. 0) &
    call get_cmdline(version_requested)

  ! See if we asked for the version number. If we did, exit quietly.
  call mpi_bcast(version_requested, 1, mpi_logical, 0, mpi_comm_world, &
                 err_code_mpi)
  if (version_requested) call mexit(6, 0)

  ! Broadcast multipmemd variables taken from command-line to everyone else so
  ! they'll know what to do

  call mpi_bcast(numgroups, 1, mpi_integer, 0, mpi_comm_world, err_code_mpi)
  call mpi_bcast(groupfile_name, max_fn_len, mpi_character, 0, &
                 mpi_comm_world, err_code_mpi)
  call mpi_bcast(proc_map_name, max_fn_len, mpi_character, 0, &
                 mpi_comm_world, err_code_mpi)
  call mpi_bcast(remd_method, 1, mpi_integer, 0, mpi_comm_world, err_code_mpi)
  call mpi_bcast(remd_random_partner, 1, mpi_integer, 0, mpi_comm_world, err_code_mpi)
  ! Broadcast RREMD too for Reservoir REMD 
  call mpi_bcast(rremd_type, 1, mpi_integer, 0, mpi_comm_world, err_code_mpi)
  call mpi_bcast(ng_nonsequential, 1, mpi_logical, 0, mpi_comm_world, &
                 err_code_mpi)

  ! Set up pmemd_comm communicators. Start with mpi_comm_world and split up for
  ! additional replicas if necessary. Set up a comm between the rank 0 threads
  ! of each separate pmemd comm

  pmemd_comm = mpi_comm_world
  numtasks = worldsize
  mytaskid = worldrank

  if (numgroups .gt. 1) then

    ! If no processor map file was given on the command-line, then use either
    ! a sequential (default) or non-sequential (via -ng-nonsequential) processor
    ! distribution, and require worldsize be a multiple of numgroups. If we were
    ! given a processor map file, then use that to allocate processors. We allow
    ! for numgroups - 1 lines, and all remaining, unassigned threads will be
    ! assigned to the last pmemd job

    if (len_trim(proc_map_name) .eq. 0) then

      if ( mod(worldsize, numgroups) .ne. 0 ) &
        call comm_error('setup_groups: MPI size is not a multiple of -ng')

      if (ng_nonsequential) then
        pmemd_comm_number = mod(worldrank, numgroups)
      else
        pmemd_comm_number = worldrank / (worldsize / numgroups)
      end if

    else

      allocate(proc_map_lines(numgroups), stat=alloc_failed)

      if (worldrank .eq. 0) call get_proc_map(proc_map_lines)

      call mpi_bcast(proc_map_lines, grpfl_line_len * numgroups, &
                     mpi_character, 0, mpi_comm_world, err_code_mpi)

      call assign_communicator(proc_map_lines)

      deallocate(proc_map_lines)

    end if

    pmemd_comm = mpi_comm_null

    call mpi_comm_split(mpi_comm_world, pmemd_comm_number, worldrank, &
                        pmemd_comm, err_code_mpi)

    if ( pmemd_comm .eq. mpi_comm_null ) &
      call comm_error('setup_groups: pmemd_comm not set up properly')
    
    call mpi_comm_size(pmemd_comm, numtasks, err_code_mpi)
    call mpi_comm_rank(pmemd_comm, mytaskid, err_code_mpi)

    pmemd_master_comm = mpi_comm_null
    masterid = MPI_UNDEFINED
    master_rank = MPI_UNDEFINED

    if ( mytaskid .eq. 0 ) &
      masterid = 0

    call mpi_comm_split(mpi_comm_world, masterid, worldrank, &
                        pmemd_master_comm, err_code_mpi)

    if (err_code_mpi .ne. MPI_SUCCESS) &
      call comm_error('setup_groups: Error creating master communicator')

    ! Create the master communicator and have master_master parse the groupfile
    ! and scatter the results

    if ( mytaskid .eq. 0 ) then
      call mpi_comm_size(pmemd_master_comm, master_size, err_code_mpi)
      call mpi_comm_rank(pmemd_master_comm, master_rank, err_code_mpi)
      master_master = master_rank .eq. 0
      

      allocate(grouplines(numgroups), stat = alloc_failed)
      
      if (master_master) then

        write(formt, '(a,i5,a)') '(a', grpfl_line_len, ')'

        call amopen(groupfile, groupfile_name, 'O', 'F', 'R')
      
        
        i = 1
        do while (i .le. numgroups)
          read(unit=groupfile, fmt=formt, end=1000) buf
  
          if (len_trim(buf) .eq. 0) cycle

          do j = 1, len_trim(buf)
             if (buf(j:j) .le. ' ') cycle
             exit
          end do

          if (buf(j:j) .eq. '#' .or. buf(j:j) .eq. '!') cycle
         
          grouplines(i) = buf
          i = i + 1

        end do

        close(groupfile)
  
      end if ! master_master

     ! Scatter the grouplines to the buffers if we have multiple groups

      call mpi_scatter(grouplines, grpfl_line_len, MPI_CHARACTER, &
                       groupline_buffer, grpfl_line_len, MPI_CHARACTER, &
                       0, pmemd_master_comm, err_code_mpi)

    end if ! mytaskid .eq. 0
    
    ! Deallocate grouplines

    if (allocated(grouplines)) deallocate(grouplines)

    ! Broadcast our replica number

    if (mytaskid .eq. 0) repnum = master_rank + 1

    call mpi_bcast(repnum, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)

    ! Write multipmemd details to stdout a la multisander

    if (master_master) then
      write(6,'(a)') ''
      write(6,'(a)') ' Running multipmemd version of pmemd Amber24'
      write(6,'(a,i5)') '    Total processors = ', worldsize
      write(6,'(a,i5)') '    Number of groups = ', numgroups
      write(6, '(a)') ''
    end if

  else
    ! with only 1 group, make pmemd_comm a copy of our mpi_comm_world
    call mpi_comm_dup(mpi_comm_world, pmemd_comm, err_code_mpi)
  end if ! numgroups .gt. 1

  ! create the pmemd_group

  call mpi_comm_group(pmemd_comm, pmemd_group, err_code_mpi)

  call mpi_barrier(mpi_comm_world, err_code_mpi)

  numtasks_recip = numtasks
  mytaskid_recip = mytaskid

  return

1000 call comm_error('setup_groups: unexpected end of file in groupfile')

end subroutine setup_groups

!*******************************************************************************
!
! Subroutine: free_comms
!
! Description: Frees the replica communicators. Only call this before an 
!              MPI_Finalize on mpi_comm_world. For MPI_Abort, no need to call
!              this. Needs to be here because in pmemd_lib we get a cyclic
!              dependency that prevents compilation
!
!*******************************************************************************

subroutine free_comms

  use parallel_dat_mod

  implicit none

  ! Free the comms and groups
  call mpi_group_free(pmemd_group, err_code_mpi)
  call mpi_comm_free(pmemd_comm, err_code_mpi)
  if (pmemd_master_comm .ne. mpi_comm_null) &
    call mpi_comm_free(pmemd_master_comm, err_code_mpi)
  if (remd_comm .ne. mpi_comm_null) &
    call mpi_comm_free(remd_comm, err_code_mpi)
  if (group_master_comm .ne. mpi_comm_null) &
    call mpi_comm_free(group_master_comm, err_code_mpi)

end subroutine free_comms

!*******************************************************************************
!
! Subroutine: comm_error
!
! Description: Called in the case of an error creating communicators
!
!*******************************************************************************

subroutine comm_error(error_message)

  use pmemd_lib_mod, only : mexit

  implicit none

  character(*)    :: error_message

  write(0,'(a)') error_message

  call mexit(6,1)

end subroutine comm_error

!*******************************************************************************
!
! Subroutine: get_proc_map
!
! Description: Parses the processor mapping file and fills the array of lines
!              so each processor can figure out where it belongs.
!
!*******************************************************************************

subroutine get_proc_map(proc_map_lines)

  use file_io_mod,      only : amopen
  use file_io_dat_mod,  only : grpfl_line_len, proc_map_name, proc_map
  use parallel_dat_mod, only : numgroups

  implicit none

! Passed variables

  character (len=grpfl_line_len), intent(out) :: proc_map_lines(numgroups)

! Local variables

  character (len=grpfl_line_len) :: buf
  character (len=10)             :: formt ! format
  integer i, j

  proc_map_lines(:) = ' '

  call amopen(proc_map, proc_map_name, 'O', 'F', 'R')

  write(formt, '(a,i5,a)') '(a', grpfl_line_len, ')'

  i = 1
  do while (i .le. numgroups)

    read(unit=proc_map, fmt=formt, end=1000) buf

    ! Skip over comments and blank lines

    if (len_trim(buf) .eq. 0) cycle

    do j = 1, len_trim(buf)
       if (buf(j:j) .le. ' ') cycle
       exit
    end do

    if (buf(j:j) .eq. '#' .or. buf(j:j) .eq. '!') cycle

    proc_map_lines(i) = buf
    i = i + 1
    
  end do

  return

  close(proc_map)

! Now we allow for one less line than group, and assign remaining procs to
! the last thread

  1000 continue
  
  if (i .lt. numgroups - 1) call comm_error('unexpected end of file for ' // &
                                             trim(proc_map_name))
  
  return

end subroutine get_proc_map

!*******************************************************************************
!
! Subroutine: assign_communicator
!
! Description: Assigns each processor to a communicator based on file input
!
!*******************************************************************************

subroutine assign_communicator(proc_map_lines)

  use file_io_dat_mod,  only : grpfl_line_len
  use parallel_dat_mod, only : numgroups, worldsize, worldrank, pmemd_comm_number
  use pmemd_lib_mod,    only : get_num_tokens, get_token

  implicit none

! Passed variables

  character (len=grpfl_line_len), intent(in) :: proc_map_lines(numgroups)

! Local variables

  character (len=10)    :: token    ! holder for token parsed from string

  integer               :: proc_num ! holder for integer value of above token
  integer               :: num_tokens
  integer               :: i 
  integer               :: j

  pmemd_comm_number = -1

  ! Now parse through each line of proc_map_lines. The last line may be blank,
  ! other than that they should all be full. Those errors will have been caught
  ! in get_proc_map. The errors we have to catch here are specifying a non-
  ! existent processor ( < 0 or > worldsize), double-specifying a processor, or
  ! NOT specifying a processor, and the last line is not blank.

  do i = 1, numgroups

    call get_num_tokens(proc_map_lines(i), num_tokens)

    do j = 1, num_tokens

      call get_token(proc_map_lines(i), j, token)

      read(token, '(i5)', err=666) proc_num
      
      ! Check for errors here

      if (proc_num .lt. 0) then
        call comm_error('negative processor assigned in processor map file!')
      else if (proc_num .ge. worldsize) then
        call comm_error('Only processors 0 to ng-1 can be assigned in &
                        &processor map file!')
      else if (proc_num .eq. worldrank .and. pmemd_comm_number .ne. -1) then
        call comm_error('Processor assignments must be unique!')
      else if (proc_num .eq. worldrank) then
        pmemd_comm_number = i
      end if
      
    end do

  end do

  ! Check out the last line here, to see if it's blank, or if any procs have
  ! been left over.

  if (len_trim(proc_map_lines(numgroups)) .eq. 0) then
    if (pmemd_comm_number .eq. -1 ) &
      pmemd_comm_number = numgroups
  else if (pmemd_comm_number .eq. -1) then
    call comm_error('not all processors were assigned to threads!')
  end if

  return

  666 continue

  call comm_error('error reading integer from processor map file!')

end subroutine assign_communicator

#endif /* MPI */

end module multipmemd_mod
