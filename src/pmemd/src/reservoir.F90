#include "copyright.i"

!*******************************************************************************
!
! Module: reservoir_mod
!
! Description: 
!
! Module for holding the necessary subroutines required for Reservoir REMD.
! The sander implementation was done by Daniel Roe, based on the original
! implementation by Guanglei Cui and Carlos Simmerling. Adapted to PMEMD by Koushik K.
!
! The main exchange subroutines for this REMD calculation are still
! in remd.F90 and remd_exchg.F90. Here, only the non-exchange related subroutines
! are included, which are listed below.
!
! Reservoir REMD subroutines:
! load_reservoir_files() - Loading the Reservoir files
! load_reservoir_structure() - Load a reservoir structure if exchange with the 
!                              reservoir was successful
!
!*******************************************************************************

module reservoir_mod

! ... Constants:
!  rremd_type            : The type of reservoir replica exchange
!                           0 - no reservoir (default)
!                           1 - Boltzmann weighted reservoir
!                           2 - 1/N weighted reservoir
!                           3 - reservoir with weights defined by
!                                dihedral clustering (not yet implemented)
!  reservoir_velocity    : If the reservoir restart files have velocities
!  reservoir_structure_energies_array() : Energy of each structure in the reservoir
!  reservoir_temperature : temperature of reservoir
!  rremd_idx             : index of random structure from reservoir
!  reservoir_size        : total number of structures in reservoir
!  reservoir_ncid        : NCID of netcdf reservoir; -1 when not netcdf

  logical, save                          :: reservoir_velocity(2)
  double precision, allocatable, save    :: reservoir_structure_energies_array(:)
  double precision, save                 :: reservoir_temperature(2)
  integer, save                          :: rremd_idx(2) = -1
  integer, save                          :: reservoir_size(2) = -1
  integer, save                          :: reservoir_ncid(2) = -1
  integer, save                          :: reservoir_iseed

! required for netcdf reservoir
#ifdef BINTRAJ
   integer, save :: coordVID, velocityVID, cellAngleVID, cellLengthVID
#endif

! NOT IMPLEMENTED YET
! Reservoir REMD Dihedral Clustering

! Specific for Reservoir REMD Dihedral Clustering:
! clusternum()    : cluster that reservoir structure belongs to
! clusterid(,)    : 2-D array of bin values for each dihedral, 1 set per cluster
! clustersize()   : number of structures from reservoir in this cluster (note 0
!                   is unknown cluster, currently set to 0)
! dihclustat(,)   : 2-D, 4 atoms that define each of the dihedrals used for
!                    clustering
!     TODO: dihclustnbin is only used to hold values for output - could be axed.
! dihclust_nbin() : number of bins used for each of the dihedrals (each is
!                     dihclust_nbin/360 in width)
! dihclustmin()   : For each dihedral, torsion value that should get bin 0
! dihcluststep()  : Step size for each dihedral (calculated during read of
!                     clusterinfo as 360/Nbins).
! currdihid()     : diheddral bin values for the current MD structure
!                     (to assign it to a cluster)
! incluster       : cluster for the current MD structure
! nclust          : number of clusters
! nclustdih       : number of dihedrals used for clustering
  integer, allocatable, save :: clusternum(:)
!  integer, allocatable, save :: clusterid(:,:)
!  integer, allocatable, save :: clustersize(:)
!  integer, allocatable, save :: currdihid(:)
!  integer, allocatable, save :: dihclustnbin(:)
!  integer, allocatable, save :: dihclustat(:,:)
!  double precision, allocatable, save :: dihclustmin(:)
!  double precision, allocatable, save :: dihcluststep(:)
!  integer, save :: incluster
!  integer, save :: nclust
!  integer, save :: nclustdih
! Reservoir REMD Dihedral Clustering Not Implemented yet

! End Reservoir REMD

contains

!*******************************************************************************
!
! Subroutine: load_reservoir_files
!
! Description: Load reservoir structure files for Reservoir REMD.
!              Implemented with remd_method .eq. 5 or remd_method .eq. 6
!
!*******************************************************************************

subroutine load_reservoir_files(rremd_type)

#ifdef BINTRAJ
  use AmberNetcdf_mod, only : NC_checkTraj, NC_openRead, NC_setupReservoir, &
                              NC_readReservoir, NC_setupBox
#endif
  use file_io_mod,     only : amopen
  use pmemd_lib_mod,   only : mexit, alloc_error
  use parallel_dat_mod
  use file_io_dat_mod
  use mdin_ctrl_dat_mod

  implicit none

! Formal variable
  
  integer, intent(in)    :: rremd_type

! Local variables

  integer                :: my_atom_count(2)
  integer                :: alloc_failed
  integer                :: i, load_index
  integer                :: ios
  integer                :: local_reservoir_velocity
  
#ifdef MPI
  integer, dimension(mpi_status_size) :: stat_array ! needed for mpi_recv
#endif  

#ifdef BINTRAJ
  integer       :: eptotVID
  integer       :: binsVID
#endif
  
  load_index=0
  
  if(rremd_type.le.3) &  ! non-TI
     load_index=1
#ifdef MPI
  if(remd_rank.eq.(numgroups-1) .and. (rremd_type.eq.5 .or. rremd_type.eq.6)) &  ! TI, 2nd resrvoir
     load_index=2
  if (remd_rank.eq.0 .and. (rremd_type.eq.4 .or. rremd_type.eq.6)) &  ! TI, 1st reservoir
     load_index=1
#endif

  if (load_index.le.0) return

#ifndef BINTRAJ
      write (remtype,'(a)') &
          'Reservoir: only netcdf format is supported now, Please re-compile with the BINTRAJ flag'
      call mexit(mdout,1)
   return !only netcdf format is supported now
#endif

#ifdef BINTRAJ
    ! Check if "reservoir_name" is a Netcdf Reservoir
  if (NC_checkTraj(reservoir_name(load_index) )) then
      if (NC_openRead(reservoir_name(load_index), reservoir_ncid(load_index) )) &
        call mexit(mdout, 1)
      if (NC_setupReservoir(reservoir_ncid(load_index), &
          reservoir_size(load_index), reservoir_temperature(load_index) , &
          my_atom_count(load_index), coordVID, velocityVID, eptotVID, &
          binsVID, reservoir_iseed)) call mexit(mdout, 1)
      reservoir_velocity(load_index) = (velocityVID .ne. -1)
  else
      write(mdout,'(a)') 'Reservoir is not in Netcdf format. Please use a &
                       Netcdf reservoir or use Sander'
  end if

  ! Allocate memory for reservoir structure energies
  if (load_index.gt.0) then
    if (.not.allocated(reservoir_structure_energies_array)) &
      allocate(reservoir_structure_energies_array(reservoir_size(load_index) ), stat=alloc_failed)

    if (alloc_failed .ne. 0) &
      call alloc_error('load_reservoir_files','error allocating array for&
                  reservoir structure energies')
  endif
    
  if (rremd_type .lt. 4 .or. rremd_type .gt. 6) then
    ! Read energy for each structure. All we really need to load are the energies.
    ! if the exchange is successful we can grab the coords (and velocity if necessary
    ! from the disk during the run. See load_reservoir_structure() in remd_exchg.F90
    ! for details on how the structures are loaded.

    ! Netcdf reservoir read binsVID and clusternum are only for rremd_type = 3 which is
    ! not implemented yet
    if (.not.allocated(clusternum)) allocate(clusternum(1)) !! avoid runtime error if 
    if (NC_readReservoir(reservoir_ncid(load_index), &
        reservoir_size(load_index), eptotVID, &
        binsVID, reservoir_structure_energies_array, clusternum))&
        call mexit(mdout,1)

  endif
  
  call NC_setupBox(reservoir_ncid(load_index), cellLengthVID, cellAngleVID)
  
#endif  
  
#ifdef MPI
  !! pass the second reservoir info to the first replica
  if (rremd_type.gt.4) then
    load_index=2
    if (remd_rank.eq.(numgroups-1)) then
      call mpi_send(my_atom_count(load_index), 1, mpi_integer, &
          0, 1231, pmemd_master_comm, err_code_mpi)
      call mpi_send(reservoir_size(load_index), 1, mpi_integer, &
          0, 1232, pmemd_master_comm, err_code_mpi)
      call mpi_send(reservoir_temperature(load_index), 1, mpi_integer, &
          0, 1233, pmemd_master_comm, err_code_mpi)
    endif
    if (remd_rank.eq.0) then
      call mpi_recv(my_atom_count(load_index), 1, mpi_integer, &
          (numgroups-1), 1231, pmemd_master_comm, MPI_STATUS_IGNORE, err_code_mpi)
      call mpi_recv(reservoir_size(load_index), 1, mpi_integer, &
          (numgroups-1), 1232, pmemd_master_comm, MPI_STATUS_IGNORE, err_code_mpi)   
      call mpi_recv(reservoir_temperature(load_index), 1, mpi_integer, &
          (numgroups-1), 1233, pmemd_master_comm, MPI_STATUS_IGNORE, err_code_mpi)   
    endif
  endif

  if (master_master) then
    call amopen(remtype, remtype_name, owrite, 'F', 'W')
    
  ! Print rremd_type information
    write (remtype,'(a)') '================'
    write (remtype,'(a,i5)') 'RREMD: Reservoir type ', rremd_type
    if (rremd_type .eq. 1) then
      write (remtype,'(a)') &
          ' Boltzmann weighted reservoir exchange uses delta beta'
    else if (rremd_type .eq. 2) then
      write (remtype,'(a)') &
          ' Non-Boltzmann 1/N weighted reservoir, exchange uses beta'
    else if (rremd_type .eq. 4) then
      write (remtype,'(a)') &
          ' reservoir for TI lambd=0 real state only'
    else if (rremd_type .eq. 5) then
      write (remtype,'(a)') &
          ' reservoir for TI lambd=1 real state only'
    else if (rremd_type .eq. 6) then
      write (remtype,'(a)') &
          ' reservoirs for both TI lambd=0 and 1 real states'
    else
      write (remtype,'(a,i5)') &
          'Unknown reservoir type: rremd_type=', rremd_type
      call mexit(mdout,1)
    end if
      
    do load_index=1,2
       if(load_index.eq.1.and.rremd_type.eq.5) cycle
       if(load_index.eq.2.and.rremd_type.eq.4) cycle
        
        write (remtype,'(a)') '================'
        write (remtype,'(a,i2,a)') 'Reservoir #: ',load_index
        write (remtype,'(a,a)')    '  Info from file ', reservoir_name(load_index)
        write (remtype,'(a,i5)')   '  NumAtoms= ', my_atom_count(load_index)
        write (remtype,'(a,i5)')   '  ReservoirSize= ', reservoir_size(load_index)
        write (remtype,'(a,f6.2)') '  ReservoirTemp(K)= ', reservoir_temperature(load_index)
        write (remtype,'(a,i10)')  '  ReservoirRandomSeed= ', reservoir_iseed

        if (reservoir_velocity(load_index)) then
            write (remtype,'(a)') ' Velocities will be read from reservoir'
        else
            write (remtype,'(a)') ' Velocities will be assigned to &
                                    structure after exchange'
        end if
    enddo ! load_index

    ! Write energy in the remtype file just to make sure everything has been read properly
    if (rremd_type .lt. 4 .or. rremd_type .gt. 6) then
      do i = 1, reservoir_size(1)
        write (remtype,'(a,i10,f12.4)') "frame, energy ", i, reservoir_structure_energies_array(i)
      enddo
    end if
        
    write (remtype,'(a)') "Done reading reservoirs"
    close(remtype)
           
  end if

  ! Set random # generator based on the reservoir_iseed from reservoir_name file.

  ! In sander, the reservoir_iseed and reservoir_velocity were broadcasted to all the threads. In pmemd,
  ! all masters are aware of reservoir_iseed and reservoir_velocity. No need to broadcast them here.
  ! We also do not call amrset here. See setup_remd_randgen() on how the reservoir_iseed effects
  ! the random generator for REMD
  ! All threads do this.
  
  ! Comment taken from Sander code, made by Daniel Roe
  ! This reservoir_iseed can be used to do independent runs with different sequence of structures
  ! from the reservoir. Here we have no concept of the ig seed in md.in and we don't want to
  ! hardcode. So, the reservoir_iseed in reservoir related files is a convenient location.

#endif  // MPI  
  return

end subroutine load_reservoir_files

!*******************************************************************************
!
! Subroutine: load_reservoir_structure
!
! Description: Load reservoir structure coordinates if exchange is successful.
!
!*******************************************************************************

subroutine load_reservoir_structure(rremd_crd, rremd_vel, rremd_atm_cnt, &
                           a,b,c, load_index)

#  ifdef BINTRAJ
   use netcdf
   use AmberNetcdf_mod, only: NC_error, NC_readRestartBox
#  endif

  use file_io_mod,     only : amopen
  use pmemd_lib_mod,   only : mexit
  use mdin_ctrl_dat_mod ! For ntp check
  use prmtop_dat_mod,  only : ifbox
  use parallel_dat_mod
  use file_io_dat_mod

  implicit none

! Passed variables
  integer, intent(in)             :: rremd_atm_cnt
  double precision, intent(inout) :: rremd_crd(3,rremd_atm_cnt)
  double precision, intent(inout) :: rremd_vel(3,rremd_atm_cnt)
  double precision, intent(out), optional :: a,b,c
  integer, intent(in), optional:: load_index
  
! Local variables

  integer                :: rremd_crd_cnt
  integer                :: i, j, local_index
  integer                :: alloc_failed
  double precision, dimension(3)::box
  double precision::alpha,beta,gamma

  character(80)          :: line
   
  if(present(load_index)) then
    local_index=load_index
  else
    local_index=1
  endif
       
  rremd_crd_cnt = rremd_atm_cnt * 3

   ! If we exchanged with the reservoir, swap the coordinates
   ! here (after the inpcrd have been read), master process only.

   ! Read the coordinates and velocities for frame rremd_idx from
   ! cpptraj generated NetCDF reservoir.

   if (master) then
#ifdef VERBOSE_REMD
      write(mdout,'(21("="),a, 21("="))') 'Reservoir Read'
      write(mdout, '(a,a)') 'reservoir_name= ',trim(reservoir_name(local_index))
      write(mdout, '(a,i8)') 'rremd_idx =', rremd_idx(local_index)
#endif

#ifdef BINTRAJ
      ! Read coordinates
      if (NC_error(nf90_get_var(reservoir_ncid(local_index), coordVID, &
                             rremd_crd(1:3, 1:rremd_atm_cnt ), &
                             start = (/ 1, 1, rremd_idx(local_index) /), &
                             count = (/ 3, rremd_atm_cnt , 1 /)), &
                    'reading reservoir coordinates') ) call mexit(mdout,1)
      ! Read velocities
      if (reservoir_velocity(local_index)) then
         if (NC_error(nf90_get_var(reservoir_ncid(local_index), velocityVID, &
                                   rremd_vel(1:3, 1:rremd_atm_cnt), &
                                   start = (/ 1, 1, rremd_idx(local_index) /), &
                                   count = (/ 3, rremd_atm_cnt, 1 /)), &
                      'reading reservoir velocities')) call mexit(mdout,1)
      end if

      if (ifbox .gt. 0) then
      ! Read box information if present, not necessary for NetCDF
      ! trajectories because the indexing is already done
        
        if (present(a) .and. present(b) .and. present(c))  then 
          if (.not. NC_error(nf90_get_var(reservoir_ncid(local_index), cellLengthVID, box(1:3), &
                                  start = (/1, rremd_idx(local_index) /),  count = (/ 3, 1/)),&
                     'Getting box lengths')) then
            a=box(1)
            b=box(2)
            c=box(3)
          else
            if (NC_readRestartBox(reservoir_ncid(rremd_idx(local_index) ),a,b,c,alpha,beta,gamma)) then
              call mexit(mdout,1)
            end if   
          endif
        endif
      end if
      

#endif

#ifdef VERBOSE_REMD
      write(mdout, '(a,i8)') 'RREMD: coords of first 5 atoms read for frame ', rremd_idx
      ! Print coordinates for reservoir structure
      do j = 1 , 5
        write(mdout,*) (rremd_crd(i,j), i=1, 3)
      end do

      if (reservoir_velocity) then
           write(mdout,'(a,i8)') "RREMD: velocities of first 5 atoms read for frame ", rremd_idx
      ! Print velocities for reservoir structure
         do j = 1, 5
            write(mdout,*) (rremd_vel(i,j), i=1, 3)
         end do
      end if
#endif

      ! If we have const P, need box change. Currently unsupported
      ! TI is an exception (for TI-reservoir))
      if (ntp .gt. 0 .and. icfe .eq. 0) then
         write (mdout,'(a)') 'const P not allowed for Reservoir REMD'
         call mexit (mdout,1)
      end if
   end if ! master 

#ifdef VERBOSE_REMD
   if (master) then
      write (mdout,'(21("="),a,21("="))') 'End Reservoir Structure Read'
   end if
#endif

   return

end subroutine load_reservoir_structure

end module reservoir_mod
