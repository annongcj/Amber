#include "copyright.i"

!*******************************************************************************
! Module:  binrestart_mod
!
! Description:  Module for generating Amber NetCDF restart files. Originally 
!               developed by Dan Roe, August 2011. Based on the Amber NetCDF 
!               trajectory format originally developed by John Mongan.
!*******************************************************************************

module binrestart_mod
   private

   integer, save :: coordVID, velocityVID, cellAngleVID, cellLengthVID
   integer, save :: timeVID, repidx_var_id, crdidx_var_id
   integer, save :: remd_indices_var_id, remd_groups_var_id
   integer, save :: remd_values_var_id, remd_types_var_id
   integer, save :: remd_dimension_var_id

   ! If WRITE_NC_RESTART is called for the main restart file, this variable
   ! indicates whether the file has been set up.
   logical, save :: mainRestartSetup = .false.

   public write_nc_restart, &
          read_nc_restart, &
          read_nc_traj, &
          read_nc_restart_atoms, &
          read_nc_frame_cnt, &
          read_nc_refc, &
          read_nc_remd_dimension, &
          read_nc_remd_types
contains

#ifndef BINTRAJ
subroutine no_nc_error
  use AmberNetcdf_mod, only : NC_NoNetcdfError
  use file_io_dat_mod, only : mdout
  use pmemd_lib_mod,   only : mexit
  call NC_NoNetcdfError(mdout)
  call mexit(mdout, 1)
end subroutine no_nc_error
#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Write Netcdf restart file.
!-------------------------------------------------------------------
!     --- WRITE_NC_RESTART ---
!-------------------------------------------------------------------
!     Write Netcdf Restart file with given filename and title. 
!     owrite indicates overwrite status (N is no overwrite), natom 
!     is the # atoms, ntb>0 indicates presence of box coords, isMain
!     indicates this is the main restart file. The main restart file
!     can only be created and set-up once, all others will always
!     be created and set-up. 
!     Coords and Velo are the coordinates and velocities, temp0 is 
!     the current  temperature and Time is the current simulation 
!     time. If imin is 0, this is for MD so velocities will be 
!     written, otherwise no velocity info will be saved.
subroutine write_nc_restart(filename,title,owrite,natom,ntb,isMain,Coords,Velo,temp0,Time,imin,&
                            box,alpha,beta,gamma,&
#ifdef MPI
                            remd_method,remd_dimension,remd_types,group_num,replica_indexes,&
                            remd_repidx, remd_crdidx, &
#endif
                           solvph,solve)
#ifdef BINTRAJ
  use netcdf
  use AmberNetcdf_mod, only : NC_create, NC_defineRemdIndices, checkNCerror,&
                              NC_openWrite, NC_close
  use axis_optimize_mod
  use pmemd_lib_mod,   only : mexit
  use file_io_dat_mod, only : mdout
#endif
  use mdin_ctrl_dat_mod, only : usemidpoint
  implicit none
  ! Formal arguments:
  character(len=*), intent(in)                 :: filename
  character(len=*), intent(in)                 :: title
  character, intent(in)                        :: owrite
  integer, intent(in)                          :: natom,ntb
  logical, intent(in)                          :: isMain 
  double precision, dimension(:,:), intent(in) :: Coords, Velo
  double precision, intent(in)                 :: temp0, Time, solvph
  double precision, intent(in)                 :: solve
  integer, intent(in)                          :: imin
  double precision, dimension(3), intent(in)   :: box
  double precision, intent(in)                 :: alpha, beta, gamma
#ifdef MPI
  integer, intent(in)               :: remd_method, remd_dimension
  integer, intent(in), dimension(:) :: remd_types, group_num, replica_indexes
  integer, intent(in)               :: remd_repidx, remd_crdidx
#endif
#ifdef BINTRAJ
  ! Local variables:
  integer :: ncid, natom3, ord1, ord2, ord3
  logical :: has1DRemdValues = .false.
  integer :: frcVID ! dummy variable
# ifdef MPI
  double precision, dimension(remd_dimension) :: remd_values
  integer :: i
  has1DRemdValues = (remd_method.gt.0)
# endif
  ! ------------------------- File Setup --------------------------
  if ( (isMain .eqv. .false.) .or. (mainRestartSetup .eqv. .false.) ) then
    ! If first call, create the file and set up all dimensions and vars
    if (NC_create(filename, owrite, .true., natom, .true., &
                  (imin.eq.0), (ntb.gt.0), has1DRemdValues, .true., .false., &
                  title, ncid, timeVID, coordVID, velocityVID, frcVID, &
                  cellLengthVID, cellAngleVID, remd_values_var_id)) call mexit(mdout,1)
#   ifdef MPI
    ! REMD indices
    if (NC_defineRemdIndices(ncid, remd_dimension, remd_indices_var_id,&
                               repidx_var_id, crdidx_var_id,&
                               remd_types, .true., (remd_method .ne. 0),&
                               (remd_method .eq. -1),&
                               remd_groupsVID=remd_groups_var_id,&
                               remd_valuesVID=remd_values_var_id,&
                               remd_typesVID=remd_types_var_id ))&
      call mexit(mdout,1)
#   endif
    ! If this is the main restart, indicate it has been setup
    if (isMain) mainRestartSetup=.true.
  else
    ! If not the first call, just reopen the existing file
    if (NC_openWrite(filename, ncid)) then 
      write (mdout,'(a)') 'write_nc_restart(): Could not open restart'
      call mexit(mdout,1)
    endif
  endif

  ! ------------------------- File Write --------------------------
  ord1 = axis_flipback_ords(1)
  ord2 = axis_flipback_ords(2)
  ord3 = axis_flipback_ords(3)

  ! Write time
  call checkNCerror(nf90_put_var(ncid, timeVID, Time), 'write time')
  ! Write coords
  call write_nc_coords(ncid,coordVID,natom,Coords,ord1,ord2,ord3)
  ! Write velocities
  if (imin .eq. 0) then
    call write_nc_coords(ncid,velocityVID,natom,Velo,ord1,ord2,ord3)
  endif
  ! Write box information
  if (ntb > 0) then
    call checkNCerror(nf90_put_var(ncid,cellLengthVID, &
                                   (/ box(ord1),box(ord2),box(ord3) /), &
                                   start = (/ 1 /), count = (/ 3 /) ), &
                      'write cell lengths')
    call checkNCerror(nf90_put_var(ncid,cellAngleVID, &
                                   (/ alpha,beta,gamma /), &
                                   start = (/ 1 /), count = (/ 3 /) ), &
                      'write cell angles')
  endif
# ifdef MPI
  ! Write replica temperature, solvent pH, redox potential, indices
  if (remd_method.ne.0) then
    ! Write overall coordinate and replica index
     call checkNCerror(nf90_put_var(ncid, repidx_var_id, remd_repidx), &
                       'write overall replica index')
     call checkNCerror(nf90_put_var(ncid, crdidx_var_id, remd_crdidx), &
                       'write overall coordinate index')
     call checkNCerror(nf90_put_var(ncid, remd_indices_var_id, replica_indexes(:), &
                       start = (/ 1 /), count = (/ remd_dimension /)), &
                       'write replica index for each dimension')
    ! multi-D remd: Store indices of this replica in each dimension
    if (remd_method .eq. -1) then
      call checkNCerror(nf90_put_var(ncid, remd_groups_var_id, group_num(:), &
                        start = (/ 1 /), count = (/ remd_dimension /)), &
                        'write replica group for each dimension')
      ! Preparing remd_values vector
      do i = 1, remd_dimension
        if (remd_types(i) == 1) then
          remd_values(i) = temp0
        else if (remd_types(i) == 3) then
          remd_values(i) = replica_indexes(i)
        else if (remd_types(i) == 4) then
          remd_values(i) = solvph
        else if (remd_types(i) == 5) then
          remd_values(i) = solve
        end if
      end do
      call checkNCerror(nf90_put_var(ncid, remd_values_var_id, remd_values(:), &
                        start = (/ 1 /), count = (/ remd_dimension /)), &
                        'write replica values to be restarted for each dimension')
      call checkNCerror(nf90_put_var(ncid, remd_types_var_id, remd_types(:), &
                        start = (/ 1 /), count = (/ remd_dimension /)), &
                        'write remd_types')
    endif
    if (remd_method .eq. 1) then
      call checkNCerror(nf90_put_var(ncid, remd_values_var_id, temp0), 'write temp0')
    else if (remd_method .eq. 4) then
      call checkNCerror(nf90_put_var(ncid, remd_values_var_id, solvph), 'write solvph')
    else if (remd_method .eq. 5) then
      call checkNCerror(nf90_put_var(ncid, remd_values_var_id, solve), 'write solve')
    end if
  endif
# endif
  ! Close restart file
  call NC_close(ncid)
#else
  call no_nc_error
#endif
end subroutine write_nc_restart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read number of frames from a Netcdf trajectory file.
!-------------------------------------------------------------------
!     --- READ_NC_FRAME_CNT ---
!-------------------------------------------------------------------
!     This routine is necessary so PMEMD can allocate memory for 
!     input coordinates before reading them in.
!     It is assumed file has already been checked for conventions.
subroutine read_nc_frame_cnt(filename, trajframes)

#ifdef BINTRAJ
  use AmberNetcdf_mod, only : NC_openRead, GetDimInfo, NC_close
  use file_io_dat_mod, only : mdout
  use pmemd_lib_mod,   only : mexit
#endif
  implicit none
  ! Formal arguments
  character(len=*), intent(in) :: filename
  integer, intent(out)         :: trajframes
#ifdef BINTRAJ
  ! Local variables
  integer :: ncid, err

  ! Open file
  if (NC_openRead(filename, ncid)) then
    write(mdout,'(a)') "read_nc_restart_atoms(): Could not open traj file."
    call mexit(mdout,1)
  endif
  ! Get number of atoms. err is used as a dummy arg for coordVID and velocityVID
  if (GetDimInfo(ncid, "frame", trajframes) .eq. 0) call mexit(mdout,1)
  ! Close file
  call NC_close(ncid)
#else
  call no_nc_error
#endif

end subroutine read_nc_frame_cnt

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read number of atoms from a Netcdf restart file.
!-------------------------------------------------------------------
!     --- READ_NC_RESTART_ATOMS ---
!-------------------------------------------------------------------
!     This routine is necessary so PMEMD can allocate memory for 
!     input coordinates before reading them in.
!     It is assumed file has already been checked for conventions.
subroutine read_nc_restart_atoms(filename,natom)
#ifdef BINTRAJ
  use AmberNetcdf_mod, only : NC_openRead, NC_setupCoordsVelo, NC_close
  use file_io_dat_mod, only : mdout
  use pmemd_lib_mod,   only : mexit
#endif
  implicit none
  ! Formal arguments
  character(len=*), intent(in) :: filename
  integer, intent(out)         :: natom
#ifdef BINTRAJ
  ! Local variables
  integer :: ncid, err

  ! Open file
  if (NC_openRead(filename, ncid)) then
    write(mdout,'(a)') "read_nc_restart_atoms(): Could not open coordinate file."
    call mexit(mdout,1)
  endif
  ! Get number of atoms. err is used as a dummy arg for coordVID and velocityVID
  if (NC_setupCoordsVelo(ncid, natom, err, err)) call mexit(mdout,1)
  ! Close file
  call NC_close(ncid) 
#else
  call no_nc_error
#endif
end subroutine read_nc_restart_atoms

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read coord/velocity information from a Netcdf restart file.
!-------------------------------------------------------------------
!     --- READ_NC_RESTART ---
!-------------------------------------------------------------------
!     Read coordinates and velocities from the Netcdf Restart file 
!     with specified filename. Title will be read in and set. 
!     natom is the number of atoms in the restart file.
!     Coords and Velo are the coordinates and velocities, temp0 is 
!     the temperature (if present) and Time is the time.
!     box_found is set to true if restart contains box coords.
!     velocities_found set to true if restart contains velocities.
subroutine read_nc_restart(filename,title,natom,Coords,Velo,remd_values,remd_values_dim,Time,&
                           box,alpha,beta,gamma,box_found,velocities_found)
#ifdef BINTRAJ
  use netcdf
  use AmberNetcdf_mod, only : NC_openRead, NC_error, NC_setupRestart, &
                              NC_readRestartBox, NC_close
  use file_io_dat_mod, only : mdout
  use pmemd_lib_mod,   only : mexit
#endif
  implicit none
  ! Formal Arguments
  character(len=*), intent(in)                  :: filename
  character(len=80), intent(out)                :: title
  integer, intent(in)                           :: natom, remd_values_dim
  double precision, dimension(:,:), intent(out) :: Coords, Velo
  double precision, intent(out)                 :: Time
  double precision, dimension(remd_values_dim), intent(out) :: remd_values
  double precision, dimension(3), intent(out)   :: box
  double precision, intent(out)                 :: alpha, beta, gamma
  logical, intent(out)                          :: box_found, velocities_found
#ifdef BINTRAJ
  ! Local variables
  integer :: ncid, ncatom

  ! ---=== Open file
  if (NC_openRead(filename, ncid)) then
    write(mdout,'(a)') "read_nc_restart(): Could not open coordinate file."
    call mexit(mdout,1)
  endif
  ! Setup restart: title, coordVID, velocityVID, time, remd_values_var_id
  if (NC_setupRestart(ncid, title, ncatom, coordVID, velocityVID,&
                      remd_values_var_id, Time)) call mexit(mdout,1)
  ! Check number of atoms; this may be redundant since read_nc_restart
  ! is called after read_nc_restart_atoms, but since the design of this
  ! module and the code flow doesn't enforce that calling sequence, better
  ! safe than sorry; srb aug 2013.
  if (natom .ne. ncatom) then
    write(mdout, '(/2x,a)') 'FATAL: NATOM mismatch in NetCDF constraint coord and prmtop files'
    call mexit(mdout, 1)
  end if
  ! Get coords
  if (NC_error(nf90_get_var(ncid, coordVID, Coords(1:3,1:ncatom), &
                            start = (/ 1, 1 /), count = (/ 3, ncatom /)),&
               'reading restart coordinates')) call mexit(mdout,1)
  ! Get velocities
  if (velocityVID.ne.-1) then
    velocities_found=.true.
    if (NC_error(nf90_get_var(ncid, velocityVID, Velo(1:3,1:ncatom), &
                              start = (/ 1, 1 /), count = (/ 3, ncatom /)),&
                 "read_nc_restart(): Getting velocities")) call mexit(mdout,1)
  endif
  ! Get box information
  if (NC_readRestartBox(ncid,box(1),box(2),box(3),alpha,beta,gamma)) then
    box_found = .false.
  else
    box_found = .true.
  endif
  ! ---=== Replica remd_values
  if (remd_values_var_id .ne. -1) then
    if (NC_error(nf90_get_var(ncid, remd_values_var_id, remd_values, &
                              start = (/ 1 /), count = (/ remd_values_dim /)),&
                 "read_nc_restart(): Getting restart remd_values")) &
      call mexit(mdout,1)
  else
    remd_values(:)=1234321
  endif

  ! NOTE: TO BE ADDED
  !labelDID;
  !int cell_spatialDID, cell_angularDID;
  !int spatialVID, cell_spatialVID, cell_angularVID;
  
  ! ---=== Close file
  call NC_close(ncid) 
#else
  call no_nc_error
#endif
end subroutine read_nc_restart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read coord information from a Netcdf trajectory file.
!-------------------------------------------------------------------
!     --- READ_NC_TRAJ ---
!-------------------------------------------------------------------
!     Read coordinates and velocities from the Netcdf trajectory file
!     with specified filename. Title will be read in and set. 
!     natom is the number of atoms in the restart file.
!     Coords and Velo are the coordinates and velocities, temp0 is 
!     the temperature (if present) and Time is the time.
!     box_found is set to true if restart contains box coords.
!     velocities_found set to true if restart contains velocities.
subroutine read_nc_traj(filename,title,natom,Coords,Velo,remd_values,remd_values_dim,Time,&
                           box,alpha,beta,gamma,box_found,velocities_found, cur_frame)
#ifdef BINTRAJ
  use netcdf
  use AmberNetcdf_mod, only : NC_openRead, NC_error, NC_setupTraj, &
                              NC_readTrajBox, NC_close
  use file_io_dat_mod, only : mdout
  use pmemd_lib_mod,   only : mexit
#endif
  implicit none
  ! Formal Arguments
  character(len=*), intent(in)                  :: filename
  character(len=80), intent(out)                :: title
  integer, intent(in)                           :: natom, remd_values_dim
  double precision, dimension(:,:), intent(out) :: Coords, Velo
  double precision, intent(out)                 :: Time
  double precision, dimension(remd_values_dim), intent(out) :: remd_values
  double precision, dimension(3), intent(out)   :: box
  double precision, intent(out)                 :: alpha, beta, gamma
  logical, intent(out)                          :: box_found, velocities_found
  integer, intent(in)                           :: cur_frame
#ifdef BINTRAJ
  ! Local variables
  integer :: ncid, ncatom

  ! ---=== Open file
  if (NC_openRead(filename, ncid)) then
    write(mdout,'(a)') "read_nc_restart(): Could not open coordinate file."
    call mexit(mdout,1)
  endif
  ! Setup restart: title, coordVID, velocityVID, time, remd_values_var_id
  if (NC_setupTraj(ncid, title, ncatom, coordVID, velocityVID,&
                      remd_values_var_id, Time)) call mexit(mdout,1)
  ! Check number of atoms; this may be redundant since read_nc_restart
  ! is called after read_nc_restart_atoms, but since the design of this
  ! module and the code flow doesn't enforce that calling sequence, better
  ! safe than sorry; srb aug 2013.
  if (natom .ne. ncatom) then
    write(mdout, '(/2x,a)') 'FATAL: NATOM mismatch in NetCDF constraint coord and prmtop files'
    call mexit(mdout, 1)
  end if
  ! Get coords
  if (NC_error(nf90_get_var(ncid, coordVID, Coords(1:3,1:ncatom), &
                            start = (/ 1,1,cur_frame /), count = (/ 3, ncatom,1 /)),&
               'reading restart coordinates')) call mexit(mdout,1)
  ! Get box information
  if (NC_readTrajBox(ncid,box(1),box(2),box(3),alpha,beta,gamma, cur_frame)) then
    box_found = .false.
  else
    box_found = .true.
  endif
  ! ---=== Replica remd_values
  if (remd_values_var_id .ne. -1) then
    if (NC_error(nf90_get_var(ncid, remd_values_var_id, remd_values, &
                              start = (/ 1 /), count = (/ remd_values_dim /)),&
                 "read_nc_restart(): Getting restart remd_values")) &
      call mexit(mdout,1)
  else
    remd_values(:)=1234321
  endif

  ! NOTE: TO BE ADDED
  !labelDID;
  !int cell_spatialDID, cell_angularDID;
  !int spatialVID, cell_spatialVID, cell_angularVID;

  ! ---=== Close file
  call NC_close(ncid)
#else
  call no_nc_error
#endif
end subroutine read_nc_traj

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read coord information from a Netcdf restart file.
!-------------------------------------------------------------------
!     --- READ_NC_REFC ---
!-------------------------------------------------------------------
!     Read coordinates from the Netcdf Restart file with specified 
!     filename. Title will be read in and set. 
!     natom is the expected number of atoms in the restart file.
!     Coords are the coordinates. 
subroutine read_nc_refc(filename,title,natom,Coords)
#ifdef BINTRAJ
  use netcdf
  use AmberNetcdf_mod, only: NC_openRead, NC_setupRestart, NC_error, NC_close 
  use file_io_dat_mod, only: mdout
  use pmemd_lib_mod, only: mexit
#endif
  implicit none
  ! Formal Arguments
  character(len=*), intent(in)                :: filename
  character(len=80), intent(out)              :: title
  integer, intent(in)                         :: natom
  double precision, dimension(*), intent(out) :: Coords
#ifdef BINTRAJ
  ! Local variables
  integer :: ncid, ncatom, nr3
  double precision Time

  ! ---=== Open file
  if (NC_openRead(filename, ncid)) then
    write(mdout,'(a)') "read_nc_refc(): Could not open coordinate file."
    call mexit(mdout,1)
  endif
  ! Setup restart: title, coordVID, velocityVID, time, remd_values_var_id
  if (NC_setupRestart(ncid, title, ncatom, coordVID, velocityVID,&
                      remd_values_var_id, Time)) call mexit(mdout,1)
  ! Check number of atoms
  if (natom .ne. ncatom) then
    write(mdout, '(/2x,a)') 'FATAL: NATOM mismatch in NetCDF constraint coord and prmtop files'
    call mexit(mdout, 1)
  end if
  nr3 = ncatom * 3
  ! Get coords
  if (NC_error(nf90_get_var(ncid, coordVID, Coords(1:nr3), &
                            start = (/ 1, 1 /), count = (/ 3, ncatom /)),&
               "read_nc_refc(): Getting coordinates")) call mexit(mdout,1)
  ! ---=== Close file
  call NC_close(ncid)  
#else
  call no_nc_error
#endif
end subroutine read_nc_refc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read remd_dimension information from a Netcdf restart file.
!-------------------------------------------------------------------
!     --- READ_NC_REMD_DIMENSION ---
!-------------------------------------------------------------------
!     Read remd_dimension from the Netcdf Restart file with specified 
!     filename. Title will be read in and set.
subroutine read_nc_remd_dimension(filename,title,remd_dimension)
#ifdef BINTRAJ
  use netcdf
  use AmberNetcdf_mod, only: NC_openRead, NC_setupRemdDimension, NC_error, NC_close 
  use file_io_dat_mod, only: mdout
  use pmemd_lib_mod, only: mexit
#endif
  implicit none
  ! Formal Arguments
  character(len=*), intent(in)                :: filename
  character(len=80), intent(out)              :: title
  integer, intent(out) :: remd_dimension
#ifdef BINTRAJ
  ! Local variables
  integer :: ncid

  ! ---=== Open file
  if (NC_openRead(filename, ncid)) then
    write(mdout,'(a)') "read_nc_remd_dimension(): Could not open restart file."
    call mexit(mdout,1)
  endif
  ! Get remd_dimension
  remd_dimension_var_id = NC_setupRemdDimension(ncid,remd_dimension)
  if (remd_dimension_var_id.eq.-1) remd_dimension = -1
  ! ---=== Close file
  call NC_close(ncid)
#else
  call no_nc_error
#endif
end subroutine read_nc_remd_dimension

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read remd_types information from a Netcdf restart file.
!-------------------------------------------------------------------
!     --- READ_NC_REMD_TYPES ---
!-------------------------------------------------------------------
!     Read remd_types from the Netcdf Restart file with specified 
!     filename. Title will be read in and set.
subroutine read_nc_remd_types(filename,title,remd_types,remd_dimension)
#ifdef BINTRAJ
  use netcdf
  use AmberNetcdf_mod, only: NC_openRead, NC_setupRemdTypes, NC_error, NC_close 
  use file_io_dat_mod, only: mdout
  use pmemd_lib_mod, only: mexit
#endif
  implicit none
  ! Formal Arguments
  character(len=*), intent(in)                    :: filename
  character(len=80), intent(out)                  :: title
  integer, intent(in)                             :: remd_dimension
  integer, dimension(remd_dimension), intent(out) :: remd_types
#ifdef BINTRAJ
  ! Local variables
  integer :: ncid

  ! ---=== Open file
  if (NC_openRead(filename, ncid)) then
    write(mdout,'(a)') "read_nc_remd_types(): Could not open restart file."
    call mexit(mdout,1)
  endif
  ! Setup restart: remd_types_var_id
  remd_types_var_id = NC_setupRemdTypes(ncid)
  if (remd_types_var_id.eq.-1) then
    remd_types(1) = -1
  else
    ! Get remd_types
    if (NC_error(nf90_get_var(ncid, remd_types_var_id, remd_types, &
                              start=(/ 1 /), count=(/ remd_dimension /)),&
                 "read_nc_remd_types(): Getting remd_types")) call mexit(mdout,1)
  end if
  ! ---=== Close file
  call NC_close(ncid)  
#else
  call no_nc_error
#endif
end subroutine read_nc_remd_types

! ======================== PRIVATE SUBROUTINES =========================
#ifdef BINTRAJ
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Write given Netcdf coordinates.
!-------------------------------------------------------------------
!     --- WRITE_NC_COORDS ---
!-------------------------------------------------------------------
!     Write given coordinates to VID in given ncid. Based on the 
!     axis, write flipped if necessary.
subroutine write_nc_coords(ncid,VID,natom,arrayIn,ord1,ord2,ord3)
  use AmberNetcdf_mod, only: checkNCerror
  use netcdf
  implicit none
  integer, intent(in)                               :: ncid, VID, natom
  double precision, dimension(3, natom), intent(in) :: arrayIn
  integer, intent(in)                               :: ord1, ord2, ord3
  ! Local variables for dealing with flipped axis
  integer :: i
  real    :: buf(3, natom)

  ! Normal axis
  if (ord1 .eq. 1 .and. ord2 .eq. 2) then

    call checkNCerror(nf90_put_var(ncid, VID, arrayIn(:,:), &
                                   start=(/ 1, 1 /), &
                                   count=(/ 3, natom /)), &
                      'NetCDF write coords')
  ! Flipped axis
  else
    do i = 1, natom
      buf(1, i) = arrayIn(ord1, i)
      buf(2, i) = arrayIn(ord2, i)
      buf(3, i) = arrayIn(ord3, i)
    end do

    call checkNCerror(nf90_put_var(ncid, VID, buf(:,:), &
                                   start=(/ 1, 1 /), &
                                   count=(/ 3, natom /)), &
                      'NetCDF write flipped coords')

  end if
end subroutine write_nc_coords

#endif

end module binrestart_mod
