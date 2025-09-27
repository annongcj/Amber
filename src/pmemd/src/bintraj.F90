#include "copyright.i"

!*******************************************************************************!
! Module:  bintraj_mod
!
! Description:  Module for generating binary trajectory-type output in NetCDF
!               format.  Originally developed by John Mongan, November 2005
!*******************************************************************************                                                                                
module bintraj_mod

  use file_io_dat_mod, only: mdout

  implicit none

  private

  ! Variables used for mdcrd:
  integer, save         :: mdcrd_ncid
  integer, save         :: coord_var_id
  integer, save         :: cell_length_var_id
  integer, save         :: cell_angle_var_id
  integer, save         :: mdcrd_time_var_id
  integer, save         :: mdcrd_veloc_var_id
  integer, save         :: mdcrd_frc_var_id
  integer, save         :: mdcrd_frame
  integer, save         :: remd_indices_var_id
  integer, save         :: remd_values_var_id
  integer, save         :: remd_types_var_id
  integer, save         :: repidx_var_id
  integer, save         :: crdidx_var_id
  ! frame_dim_id must be global since it is used by open_binary_files and 
  ! setup_remd_indices
  integer, save         :: frame_dim_id

  ! Variables used for mdvel:
  integer, save         :: mdvel_ncid
  integer, save         :: mdvel_time_var_id
  integer, save         :: mdvel_veloc_var_id
  integer, save         :: mdvel_frame

  ! Variables used for mdfrc
  integer, save         :: mdfrc_ncid
  integer, save         :: mdfrc_time_var_id
  integer, save         :: mdfrc_frc_var_id
  integer, save         :: mdfrc_frame

  ! Variables used for X-ray structure factors
#ifndef MPI
  integer, save         :: mdsf_ncid
  integer, save         :: mdsf_time_var_id
  integer, save         :: mdsf_sf_var_id
  integer, save         :: mdsf_frame
#endif
  
  ! In this context, count of atoms dumped to file:
  integer, save         :: atom_wrt_cnt

  public        open_binary_files, &
                setup_remd_indices, &
                close_binary_files, &
                write_binary_crds, &
                write_binary_vels, &
                write_binary_frcs, &
#ifndef NOXRAY
                write_binary_strfac, &
#endif
                write_binary_cell_dat, &
                end_binary_frame, &
                open_binary_files_append ! added by Feng Pan, 2015
contains

#ifndef BINTRAJ
subroutine no_nc_error
  use AmberNetcdf_mod, only: NC_NoNetcdfError
  use pmemd_lib_mod, only  : mexit
  call NC_NoNetcdfError(mdout)
  call mexit(mdout, 1)
end subroutine no_nc_error
#endif

!*******************************************************************************!
! Subroutine:  open_binary_files
!
! Description:  Open the coordinate, velocity, and energy files.
!
!*******************************************************************************
subroutine open_binary_files
#ifndef MPI
  use file_io_mod
#ifndef NOXRAY
  use xray_globals_module, only: ntwsf, num_hkl, sf_outfile
#endif
#endif
#ifdef BINTRAJ
  use AmberNetcdf_mod, only: NC_create
  use file_io_dat_mod
  use pmemd_lib_mod
  use mdin_ctrl_dat_mod
  use prmtop_dat_mod
# ifdef MPI
    use remd_mod, only : remd_method
# endif
  implicit none
  ! Local variables:
  character(80) :: vel_title
  character(7)  :: exp_stat
  logical       :: fexists
  logical       :: crd_file_exists = .false.
  logical       :: vel_file_exists = .false.
  logical       :: frc_file_exists = .false.
  logical       :: has1DRemdValues = .false.
  integer       :: ierr

  mdcrd_ncid = -1
  mdvel_ncid = -1
  mdfrc_ncid = -1
#ifndef MPI
  mdsf_ncid  = -1
#ifndef NOXRAY
#ifndef NOXRAY
  if (ntwv .eq. 0 .and. ntwx .eq. 0 .and. ntwf .eq. 0 .and. ntwsf .eq. 0) return
#else
  if (ntwv .eq. 0 .and. ntwx .eq. 0 .and. ntwf .eq. 0 .and. ntwsf .eq. 0) return
#endif
#else
  if (ntwv .eq. 0 .and. ntwx .eq. 0 .and. ntwf .eq. 0 ) return
#endif
#else
  if (ntwv .eq. 0 .and. ntwx .eq. 0 .and. ntwf .eq. 0) return
#endif
  
  if (ntwprt .gt. 0) then
    atom_wrt_cnt = ntwprt
  else
    atom_wrt_cnt = natom
  end if
# ifdef MPI
  ! Is this a run requiring replica values?
  has1DRemdValues = (remd_method .gt. 0)
# endif
  ! Set up output coordinates
  if (ntwx .gt. 0) then
    ! NOTE: No append in pmemd?
    ! Create NetCDF trajectory file
    if ( NC_create( mdcrd_name, owrite, .false., atom_wrt_cnt, &
                    .true., ntwv.lt.0, ntb.gt.0, has1DRemdValues, .true., &
                    ntwf.lt.0, prmtop_ititl, mdcrd_ncid, mdcrd_time_var_id, &
                    coord_var_id, mdcrd_veloc_var_id, mdcrd_frc_var_id, &
                    cell_length_var_id, cell_angle_var_id, remd_values_var_id, &
                    frameDID=frame_dim_id ) ) then
      write(mdout, '(a)') "Creation of NetCDF trajectory file failed."
#     ifndef DUMB
      write (0,*) 'Error on opening ', mdcrd_name
#     endif
      call mexit(mdout,1)
    endif
    mdcrd_frame = 1
  endif
  ! Set up output velocities
  if (ntwv .gt. 0) then
    ! Create velocity traj - ierr is a dummy var
    if ( NC_create( mdvel_name, owrite, .false., atom_wrt_cnt, &
                    .false., .true., .false., .false., .true., .false., &
                    prmtop_ititl, mdvel_ncid, mdvel_time_var_id, ierr, &
                    mdvel_veloc_var_id, ierr, ierr, ierr, ierr, ierr ) ) then
      write(mdout, '(a)') "Creation of NetCDF velocity file failed."
#     ifndef DUMB
      write (0,*) 'Error on opening ', mdvel_name
#     endif
      call mexit(mdout,1)
    endif
    mdvel_frame = 1
  else if (ntwv .lt. 0) then
    if (ntwx.eq.0 .or. mdcrd_ncid.eq.-1) then
      write(mdout, '(a)') "ntwx = 0 and ntwv < 0 not allowed."
      call mexit(mdout,1)
    endif
    ! Combined mdcrd/mdvel
    mdvel_ncid = mdcrd_ncid
    mdvel_frame = mdcrd_frame
  endif

  ! Set up output forces
  if (ntwf .gt. 0) then
    ! Create frc traj - ierr is a dummy var
    if ( NC_create( mdfrc_name, owrite, .false., atom_wrt_cnt, &
                    .false., .false., .false., .false., .true., .true., &
                    prmtop_ititl, mdfrc_ncid, mdfrc_time_var_id, ierr, ierr, &
                    mdfrc_frc_var_id, ierr, ierr, ierr, ierr) ) then
      write(mdout, '(a)') "Creation of NetCDF force file failed."
#     ifndef DUMB
      write(0, *) 'Error on opening ', mdfrc_name
#     endif
      call mexit(mdout, 1)
    end if
    mdfrc_frame = 1
  else if (ntwf .lt. 0) then
    if (ntwx .eq. 0 .or. mdcrd_ncid .eq. -1) then
      write(mdout, '(a)') "ntwx = 0 and ntwf < 0 not allowed."
      call mexit(mdout, 1)
    end if
    ! Combined mdcrd/mdfrc
    mdfrc_ncid = mdcrd_ncid
    mdfrc_frame = mdcrd_frame
  end if
#  ifndef MPI
#  ifndef NOXRAY
  ! Set up output structure factors
  if (ntwsf .gt. 0) then
    if (owrite .eq. 'N') then
      exp_stat = 'NEW'
    else if (owrite .eq. 'O') then
      exp_stat = 'OLD'
    else if (owrite .eq. 'U') then
      exp_stat = 'UNKNOWN'
    end if
    open(sfout, file=sf_outfile, form='unformatted', status=exp_stat, iostat=ierr)
    if (ierr .ne. 0) then
      write(mdout, '(/,2x,a,i4,a,a)') 'Unit ', sfout, ' Error on OPEN: ', sf_outfile
    end if
#if 0
    if ( NC_create( sf_outfile, owrite, .false., num_hkl, &
                    .false., .false., .false., .false., .false., .false., &
                    'This is the title', mdsf_ncid, mdsf_time_var_id, ierr, ierr, &
                    mdsf_sf_var_id, ierr, ierr, ierr, ierr) ) then
      write(mdout, '(a)') "Creation of NetCDF structure factors file failed."
#     ifndef DUMB
      write(0, *) 'Error on opening ', sf_outfile
#     endif
      call mexit(mdout, 1)
    end if
    mdsf_frame = 1
#endif
  end if
#  endif
#  endif

#else
  call no_nc_error 
#endif
  return
end subroutine open_binary_files

!*******************************************************************************!
! Subroutine:  setup_remd_indices
!
! Description:  Set up dimension information for multi-D REMD. This is in a
!               subroutine separate from open_binary_files since REMD setup
!               occurs after the call to open_binary_files.
!
!******************************************************************************
subroutine setup_remd_indices
#ifdef BINTRAJ
  use AmberNetcdf_mod, only   : NC_defineRemdIndices
  use file_io_dat_mod, only   : mdout
  use pmemd_lib_mod, only     : mexit
# ifdef MPI
  use remd_mod, only : remd_method, remd_dimension, remd_types, group_num
  use mdin_ctrl_dat_mod, only: ioutfm
# endif
  implicit none
# ifdef MPI
  ! Only setup remd indices if ioutfm and multid remd
  if ( ioutfm .ne. 0 ) then
    if ( NC_defineRemdIndices(mdcrd_ncid, remd_dimension, remd_indices_var_id, &
                              repidx_var_id, crdidx_var_id, &
                              remd_types, .false., (remd_method .ne. 0), &
                              (remd_method .eq. -1), frameDID=frame_dim_id, &
                              remd_valuesVID=remd_values_var_id, &
                              remd_typesVID=remd_types_var_id ) ) &
      call mexit(mdout,1)
  endif
# endif
#else
  call no_nc_error 
#endif
  return
end subroutine setup_remd_indices

!*******************************************************************************!
! Subroutine:  close_binary_files
!
! Description:  Close the coordinate, velocity, and energy files.
!
!*******************************************************************************
subroutine close_binary_files
#ifdef BINTRAJ
  use AmberNetcdf_mod,   only: NC_close
  use mdin_ctrl_dat_mod, only: ntwx, ntwv, ntwf
#ifndef MPI
  use file_io_dat_mod, only: sfout
#ifndef NOXRAY
  use xray_globals_module, only: ntwsf
#endif
#endif
  implicit none
  if (ntwx  .gt. 0) call NC_close(mdcrd_ncid) 
  if (ntwv  .gt. 0) call NC_close(mdvel_ncid)
  if (ntwf  .gt. 0) call NC_close(mdfrc_ncid)
#ifndef MPI
#ifndef NOXRAY
  if (ntwsf .gt. 0) close(sfout)
#endif
#if 0
  if (ntwsf .gt. 0) call NC_close(mdsf_ncid)
#endif
#endif
#else
  call no_nc_error
#endif
  return
end subroutine close_binary_files

!*******************************************************************************!
! Subroutine:   write_binary_crds
!
! Description:  Emit coordinates to trajectory file mdcrd
!
!*******************************************************************************

subroutine write_binary_crds(crd)
#ifdef BINTRAJ
  use netcdf
  use AmberNetcdf_mod,   only : checkNCerror
  use axis_optimize_mod
# ifdef MPI
  use mdin_ctrl_dat_mod, only : temp0, solvph, solve
  use remd_mod,          only : remd_method, replica_indexes, remd_dimension, &
                                remd_repidx, remd_crdidx, remd_types
  use sgld_mod,          only : trxsgld
# endif
#endif
  implicit none
  ! Formal arguments:
  double precision :: crd(3, atom_wrt_cnt)
#ifdef BINTRAJ
  ! Local variables:
  integer :: ord1, ord2, ord3
# ifdef MPI
  double precision, dimension(remd_dimension) :: remd_values
  integer :: i
# endif

  ord1 = axis_flipback_ords(1)
  ord2 = axis_flipback_ords(2)
  ord3 = axis_flipback_ords(3)

  if (ord1 .eq. 1 .and. ord2 .eq. 2) then
    call checkNCerror(nf90_put_var(mdcrd_ncid, coord_var_id, crd(:,:), &
                                   start=(/ 1, 1, mdcrd_frame /), &
                                   count=(/ 3, atom_wrt_cnt, 1 /)), &
                      'write atom coords')
  else
    call write_binary_crds_axis_flipped(crd, ord1, ord2, ord3)
  end if
# ifdef MPI
  if (remd_method .ne. 0) then
    ! Store overall replica and coordinate indices
    call checkNCerror(nf90_put_var(mdcrd_ncid, repidx_var_id, remd_repidx, &
                                   start = (/ mdcrd_frame /)), &
                      'write overall replica index')
    call checkNCerror(nf90_put_var(mdcrd_ncid, crdidx_var_id, remd_crdidx, &
                                   start = (/ mdcrd_frame /)), &
                      'write overall coordinate index')
    call checkNCerror(nf90_put_var(mdcrd_ncid, remd_indices_var_id, &
                                   replica_indexes(:), &
                                   start = (/ 1, mdcrd_frame /), &
                                   count = (/ remd_dimension, 1 /)), &
                      'write replica index for each dimension')
    ! multi-D remd: Store indices of this replica in each dimension
    if (remd_method .eq. -1) then
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
      call checkNCerror(nf90_put_var(mdcrd_ncid, remd_values_var_id, remd_values(:), &
                        start = (/ 1, mdcrd_frame /), count = (/ remd_dimension, 1 /)), &
                        'write remd_values')
    else if (remd_method .eq. 4) then
      call checkNCerror(nf90_put_var(mdcrd_ncid, remd_values_var_id, solvph, &
                                     start = (/ mdcrd_frame /)), &
                        'write replica pH')
    else if (remd_method .eq. 5) then
      call checkNCerror(nf90_put_var(mdcrd_ncid, remd_values_var_id, solve, &
                                     start = (/ mdcrd_frame /)), &
                        'write replica Redox potential')
    else if (trxsgld) then
      call checkNCerror(nf90_put_var(mdcrd_ncid, remd_values_var_id, &
                                     REAL(replica_indexes(1)), &
                                     start = (/ mdcrd_frame /)), &
                        'write SGLD replica index')
    else
      call checkNCerror(nf90_put_var(mdcrd_ncid, remd_values_var_id, temp0, &
                                     start = (/ mdcrd_frame /)), &
                        'write replica temperature')
    end if
  end if
# endif
#else
  call no_nc_error
#endif
  return
end subroutine write_binary_crds

#ifdef BINTRAJ
!*******************************************************************************!
! Subroutine:   write_binary_crds_axis_flipped
!
! Description:  Emit coordinates to trajectory file mdcrd
!
!*******************************************************************************
subroutine write_binary_crds_axis_flipped(crd, ord1, ord2, ord3)
  use netcdf
  use AmberNetcdf_mod, only: checkNCerror
  implicit none
  ! Formal arguments:
  double precision      :: crd(3, atom_wrt_cnt)
  integer               :: ord1, ord2, ord3
  ! Local variables:
  integer       :: i
  real          :: buf(3, atom_wrt_cnt)

  do i = 1, atom_wrt_cnt
    buf(1, i) = crd(ord1, i)
    buf(2, i) = crd(ord2, i)
    buf(3, i) = crd(ord3, i)
  end do

  call checkNCerror(nf90_put_var(mdcrd_ncid, coord_var_id, buf(:,:), &
                                 start=(/ 1, 1, mdcrd_frame /), &
                                 count=(/ 3, atom_wrt_cnt, 1 /)), &
                    'write atom coords')
  return
end subroutine write_binary_crds_axis_flipped
#endif /* BINTRAJ */

!*******************************************************************************!
! Subroutine:   write_binary_vels
!
! Description:  Emit velocities to trajectory or velocities file (mdcrd or
!               mdvel, dependent on value of ntwv).
!
!*******************************************************************************
subroutine write_binary_vels(velocities)
#ifdef BINTRAJ
  use netcdf
  use AmberNetcdf_mod,   only : checkNCerror
  use axis_optimize_mod, only : axis_flipback_ords
  use mdin_ctrl_dat_mod, only : ntwv
#endif
  implicit none
  ! Formal arguments:
  double precision:: velocities(3, atom_wrt_cnt)
  ! Local variables:
  integer       :: ncid
  integer       :: frame
  integer       :: var_id
  integer       :: ord1, ord2, ord3
#ifdef BINTRAJ
  if (ntwv .lt. 0) then
    ncid = mdcrd_ncid
    frame = mdcrd_frame
    var_id = mdcrd_veloc_var_id
  else
    ncid = mdvel_ncid
    frame = mdvel_frame
    var_id = mdvel_veloc_var_id
  end if

  ord1 = axis_flipback_ords(1)
  ord2 = axis_flipback_ords(2)
  ord3 = axis_flipback_ords(3)

  if (ord1 .eq. 1 .and. ord2 .eq. 2) then
    call checkNCerror(nf90_put_var(ncid, var_id, velocities(:,:), &
                                   start=(/ 1, 1, frame /), &
                                   count=(/ 3, atom_wrt_cnt, 1 /)), &
                      'write velocities')
  else
    call write_binary_vels_axis_flipped(ncid, frame, var_id, ord1, ord2, ord3, velocities)
  end if
#else
  call no_nc_error
#endif
  return
end subroutine write_binary_vels

#ifdef BINTRAJ
!*******************************************************************************!
! Subroutine:   write_binary_vels_axis_flipped
!
! Description:  Emit velocities to trajectory or velocities file (mdcrd or
!               mdvel, dependent on value of ntwv).
!
!*******************************************************************************
! FIXME: Why is this separate from write_binary_crds_axis_flipped?
subroutine write_binary_vels_axis_flipped(ncid, frame, var_id, ord1, ord2, ord3, velocities)
  use netcdf
  use AmberNetcdf_mod, only: checkNCerror
  implicit none
  ! Formal arguments:
  integer       :: ncid
  integer       :: frame
  integer       :: var_id
  integer       :: ord1, ord2, ord3
  double precision :: velocities(3, atom_wrt_cnt)
  ! Local variables:
  integer       :: i
  real          :: buf(3, atom_wrt_cnt)

  do i = 1, atom_wrt_cnt
    buf(1, i) = velocities(ord1, i)
    buf(2, i) = velocities(ord2, i)
    buf(3, i) = velocities(ord3, i)
  end do

  call checkNCerror(nf90_put_var(ncid, var_id, buf(:,:), &
                                 start=(/ 1, 1, frame /), &
                                 count=(/ 3, atom_wrt_cnt, 1 /)), &
                    'write velocities')
  return
end subroutine write_binary_vels_axis_flipped
#endif /* BINTRAJ */

!*******************************************************************************!
! Subroutine:   write_binary_frcs
!
! Description:  Emit forces to trajectory or velocities file (mdcrd or
!               mdfrc, dependent on value of ntwf).
!
!*******************************************************************************
subroutine write_binary_frcs
#ifdef BINTRAJ
  use netcdf
  use AmberNetcdf_mod, only: checkNCerror
  use axis_optimize_mod
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  implicit none
  ! Local variables:
  integer       :: ncid
  integer       :: frame
  integer       :: var_id
  integer       :: ord1, ord2, ord3

  if (ntwf .lt. 0) then
    ncid = mdcrd_ncid
    frame = mdcrd_frame
    var_id = mdcrd_frc_var_id
  else
    ncid = mdfrc_ncid
    frame = mdfrc_frame
    var_id = mdfrc_frc_var_id
  end if
# ifdef _NETCDF_DEBUG
  write(mdout,'(3(a,i6))') &
    '| DEBUG: force write ncid= ', ncid, ' frame= ', frame, ' var_id= ', var_id
# endif
  ord1 = axis_flipback_ords(1)
  ord2 = axis_flipback_ords(2)
  ord3 = axis_flipback_ords(3)

  if (ord1 .eq. 1 .and. ord2 .eq. 2) then
    call checkNCerror(nf90_put_var(ncid, var_id, atm_frc(:,:), &
                                   start=(/ 1, 1, frame /), &
                                   count=(/ 3, atom_wrt_cnt, 1 /)), &
                      'write forces')
  else
    call write_binary_frcs_axis_flipped(ncid, frame, var_id, ord1, ord2, ord3)
  end if
#else
  call no_nc_error
#endif
  return
end subroutine write_binary_frcs

#ifdef BINTRAJ
!*******************************************************************************!
! Subroutine:   write_binary_frcs_axis_flipped
!
! Description:  Emit forces to trajectory or velocities file (mdcrd or
!               mdfrc, dependent on value of ntwf).
!
!*******************************************************************************
! FIXME: Why is this separate from write_binary_crds_axis_flipped?
subroutine write_binary_frcs_axis_flipped(ncid, frame, var_id, ord1, ord2, ord3)
  use netcdf
  use AmberNetcdf_mod, only: checkNCerror
  use inpcrd_dat_mod
  implicit none
  ! Formal arguments:
  integer       :: ncid
  integer       :: frame
  integer       :: var_id
  integer       :: ord1, ord2, ord3
  ! Local variables:
  integer       :: i
  real          :: buf(3, atom_wrt_cnt)

  do i = 1, atom_wrt_cnt
    buf(1, i) = atm_frc(ord1, i)
    buf(2, i) = atm_frc(ord2, i)
    buf(3, i) = atm_frc(ord3, i)
  end do

  call checkNCerror(nf90_put_var(ncid, var_id, buf(:,:), &
                                 start=(/ 1, 1, frame /), &
                                 count=(/ 3, atom_wrt_cnt, 1 /)), &
                    'write forces (flipped)')
  return
end subroutine write_binary_frcs_axis_flipped

#endif /* BINTRAJ */

#ifndef NOXRAY
!----------------------------------------------------------------------------------------------
! write_binary_strfac: write structure factors to a binary output file
!
! Arguments:
!   Fcalc:     array of calculated strcuture factors
!   num_hkl:   the number of reflections (each reflection has a structure factor to write)
!   sfout:     output file unit
!----------------------------------------------------------------------------------------------
subroutine write_binary_strfac(Fcalc, num_hkl, sfout)

  use xray_globals_module, only : real_kind
  ! Formal arguments
  complex(real_kind) :: Fcalc(*)
  integer :: num_hkl, sfout

#ifdef BINTRAJ
  ! Local variables
  double precision :: iobuf(2, 120)
  integer :: i, j

  i = 1
  do while (i .le. num_hkl + 120)
    j = 0
    do while (i + j .le. num_hkl .and. j .lt. 120)
      iobuf(1, j + 1) =  real(Fcalc(i + j))
      iobuf(2, j + 1) = aimag(Fcalc(i + j))
      j = j + 1
    end do
    i = i + 120
    if (j .gt. 0) then
      write(sfout) iobuf(:,1:j)
    end if
  end do 
#endif
  return
end subroutine write_binary_strfac
#endif

!*******************************************************************************!
! Subroutine:   write_binary_cell_dat
!
! Description:  Emit cell length and angle data to mdcrd.
!
!*******************************************************************************

subroutine write_binary_cell_dat
#ifdef BINTRAJ
  use netcdf
  use AmberNetcdf_mod, only: checkNCerror
  use axis_optimize_mod
  use file_io_dat_mod
  use pbc_mod
  implicit none
  ! Local variables:
  integer       :: ord1, ord2, ord3

  ord1 = axis_flipback_ords(1)
  ord2 = axis_flipback_ords(2)
  ord3 = axis_flipback_ords(3)

  call checkNCerror(nf90_put_var(mdcrd_ncid, cell_length_var_id, &
                                 (/ pbc_box(ord1), pbc_box(ord2), pbc_box(ord3) /), &
                                 start=(/ 1, mdcrd_frame /), count=(/ 3, 1 /)), &
                    'write cell lengths')

  call checkNCerror(nf90_put_var(mdcrd_ncid, cell_angle_var_id, &
                                 (/ pbc_alpha, pbc_beta, pbc_gamma /), &
                                 start=(/ 1, mdcrd_frame /), count=(/ 3, 1 /)), &
                    'write cell angles')
#else
  call no_nc_error
#endif
  return
end subroutine write_binary_cell_dat

!*******************************************************************************!
! Subroutine:   end_binary_frame
!
! Description:  Write scalar data and increment frame counter.
!
!*******************************************************************************
subroutine end_binary_frame(unit)
#ifdef BINTRAJ
  use netcdf
  use AmberNetcdf_mod, only: checkNCerror
  use file_io_dat_mod
  use mdin_ctrl_dat_mod
#endif
  implicit none
  integer, intent(in) :: unit
#ifdef BINTRAJ
  if (unit .eq. mdcrd) then
    call checkNCerror(nf90_put_var(mdcrd_ncid, mdcrd_time_var_id, (/ t /), &
                      start=(/ mdcrd_frame /), count=(/ 1 /)), 'write time')
   
    call checkNCerror(nf90_sync(mdcrd_ncid))

    mdcrd_frame = mdcrd_frame + 1
    ! Sync frame counters if combined trajectory
    if (ntwv .lt. 0) mdvel_frame = mdvel_frame + 1
    if (ntwf .lt. 0) mdfrc_frame = mdfrc_frame + 1

  else if (unit .eq. mdvel) then
    call checkNCerror(nf90_put_var(mdvel_ncid, mdvel_time_var_id, (/ t /), &
                      start=(/ mdvel_frame /), count=(/ 1 /)), 'write time')
   
    call checkNCerror(nf90_sync(mdvel_ncid))

    mdvel_frame = mdvel_frame + 1
  else if (unit .eq. mdfrc) then
    call checkNCerror(nf90_put_var(mdfrc_ncid, mdfrc_time_var_id, (/ t /), &
                      start=(/ mdfrc_frame /), count=(/ 1 /)), 'write time')
    
    call checkNCerror(nf90_sync(mdfrc_ncid))

    mdfrc_frame = mdfrc_frame + 1
#ifndef MPI
  else if (unit .eq. sfout) then
    call checkNCerror(nf90_put_var(mdsf_ncid, mdsf_time_var_id, (/ t /), &
                      start=(/ mdsf_frame /), count=(/ 1 /)), 'write time')
    
    call checkNCerror(nf90_sync(mdsf_ncid))

    mdsf_frame = mdsf_frame + 1
#endif
  else
    write (mdout, *) 'Error: unhandled unit ', unit, &
                 ' selected for flush in bintraj'
  end if
#else
  call no_nc_error
#endif
  return
end subroutine end_binary_frame
!*******************************************************************************!
! Subroutine:  open_binary_files_append
!
! Description:  Open the existing coordinate, velocity, and energy files to append
!               written by Feng Pan, to use in nfe_bbmd
!
!*******************************************************************************
subroutine open_binary_files_append
#ifdef BINTRAJ
  use AmberNetcdf_mod, only: NC_openWrite, NC_setupMdcrd
  use file_io_dat_mod
  use pmemd_lib_mod
  use mdin_ctrl_dat_mod
  use prmtop_dat_mod
# ifdef MPI
  use remd_mod, only : remd_method
# endif
#ifndef MPI
#ifndef NOXRAY
  use xray_globals_module, only : ntwsf
#endif
#endif
  implicit none
  ! Local variables:
  character(80) :: vel_title
  logical       :: crd_file_exists = .false.
  logical       :: vel_file_exists = .false.
  logical       :: frc_file_exists = .false.
  logical       :: hasRemdValues = .false.
  integer       :: ierr, ncatom

  mdcrd_ncid = -1
  mdvel_ncid = -1
  mdfrc_ncid = -1
#ifndef MPI
#ifndef NOXRAY
  if (ntwv .eq. 0 .and. ntwx .eq. 0 .and. ntwf .eq. 0 .and. ntwsf .eq. 0) return
#else
  if (ntwv .eq. 0 .and. ntwx .eq. 0 .and. ntwf .eq. 0 ) return
#endif
#else
  if (ntwv .eq. 0 .and. ntwx .eq. 0 .and. ntwf .eq. 0) return
#endif
  
  if (ntwprt .gt. 0) then
    atom_wrt_cnt = ntwprt
  else
    atom_wrt_cnt = natom
  end if
# ifdef MPI
  ! Is this a run requiring replica remd_values?
  hasRemdValues = (remd_method .ne. 0)
# endif
  ! Set up output coordinates
  if (ntwx .gt. 0) then 
    if (NC_openWrite(mdcrd_name,mdcrd_ncid)) call mexit(6,1)
    if (NC_setupMdcrd(mdcrd_ncid, prmtop_ititl, mdcrd_frame, ncatom, &
                      coord_var_id, mdcrd_veloc_var_id, mdcrd_time_var_id, &
                      cell_length_var_id, cell_angle_var_id, remd_values_var_id)) call mexit(6,1)
    mdcrd_frame=mdcrd_frame+1
    ! Check for natom mismatch
    if (ncatom .ne. atom_wrt_cnt) then
        write(mdout,'(2(a,i8))') 'Error: NetCDF traj file has ', ncatom, ' atoms, but&
                         & current topology has ',  atom_wrt_cnt
        call mexit(6,1)
    endif                  
    ! Check for box mismatch
    if ( (ntb.gt.0).neqv.(cell_length_var_id.ne.-1) ) then
        write(mdout,'(a)') 'Error: Cannot append to netcdf traj, box info mismatch.'
        call mexit(6,1)
    endif
    ! Check for remd_values mismatch 
    if ( (remd_values_var_id.ne.-1).neqv.hasRemdValues ) then
        write(mdout,'(a)') 'Error: Cannot append to netcdf traj, remd_values info mismatch.'
        call mexit(6,1)
    endif

  endif

  ! Set up output velocities
  if (ntwv .ne. 0) then
      write(mdout, *) 'Error: Cannot append netcdf velocities. Please set ntwv=0 for nfe_bbmd'
      call mexit(6, 1)
  end if

  ! Set up output forces
  if (ntwf .ne. 0) then
      write(mdout, *) 'Error: Cannot append netcdf forces. Please set ntwf=0 for nfe_bbmd'
      call mexit(6, 1)
  end if

  ! Set up output forces
#  ifndef MPI
#  ifndef NOXRAY
  if (ntwsf .ne. 0) then
      write(mdout, *) 'Error: Cannot append netcdf structure factors. Please set ntwsf=0 for &
                  &nfe_bbmd'
      call mexit(6, 1)
  end if
#  endif
#  endif
#else
  call no_nc_error 
#endif
  return
end subroutine open_binary_files_append

end module bintraj_mod
