#include "copyright.i"
#include "mpi_err_check.i"

!*******************************************************************************
!
! Module: pme_fft_midpoint_mod
!
! Description: <TBS>
!
! NOTE NOTE NOTE: This code assumes frc_int .eq. 0 and should only be used
!                 under these conditions!!!
!*******************************************************************************

module pme_fft_midpoint_mod

#ifdef MPI

  use fft1d_mod
  use pme_fft_dat_mod
#include "include_precision.i"
  implicit none

! Data describing fft blocks and ownership:

! These are arrays describing how the fft grid is divided between tasks;
! these will all be of dimension fft_blks_dim.

#ifdef pmemd_SPDP
  real, allocatable       :: proc_send_buf(:)
  real, allocatable       :: proc_recv_buf(:)
#else
  double precision, allocatable       :: proc_send_buf(:)
  double precision, allocatable       :: proc_recv_buf(:)
#endif
  integer, allocatable, save            :: fft_x_offs(:)
  integer, allocatable, save            :: fft_x_cnts(:)
  integer, allocatable, save            :: fft_y_offs1(:)
  integer, allocatable, save            :: fft_y_offs2(:)
  integer, allocatable, save            :: fft_y_cnts1(:)
  integer, allocatable, save            :: fft_y_cnts2(:)
  integer, allocatable, save            :: fft_z_offs(:)
  integer, allocatable, save            :: fft_z_cnts(:)

  integer, save                         :: max_fft_x_cnts = 0
  integer, save                         :: max_fft_y_cnts1 = 0
  integer, save                         :: max_fft_y_cnts2 = 0
  integer, save                         :: max_fft_z_cnts = 0

! The "task grid" is essentially used to assign tasks to fft runs and control
! how transposes are done.  Transposes are done on groupings with 1 constant
! task index; there are notes in the code as to how it is done, but basically
! for an xy transpose the 2nd idx is constant within a group of tasks doing
! the transpose, and for a yz transpose the 1st idx is constant within a group
! of tasks doing the transpose.

  integer, allocatable, save            :: fft_task_grid(:,:)

  integer, allocatable, save            :: fft_xy_idxlst(:)
  integer, allocatable, save            :: fft_yz_idxlst(:)

  integer, save                         :: my_grid_idx1 = 0
  integer, save                         :: my_grid_idx2 = 0

! This is the common count of blks in each dimension, currently in use,
! initialized to the minimum value allowed:

  integer, save                         :: fft_blks_dim1 = 2
  integer, save                         :: fft_blks_dim2 = 2
  integer, save                         :: max_fft_blks_dim1
  integer, save                         :: max_fft_blks_dim2

  integer, save                         :: siz_fft_mpi_buf1
  integer, save                         :: siz_fft_mpi_buf2

#ifdef pmemd_SPDP
  real, parameter           :: blk_fft_workload_estimate = 0.15d0
#else
  double precision, parameter           :: blk_fft_workload_estimate = 0.15d0
#endif
contains

!*******************************************************************************
!
! Subroutine:  blk_fft_setup_midpoint
!
! Description:
!              
!*******************************************************************************

subroutine blk_fft_setup_midpoint(fft_x_dim, fft_y_dim, fft_z_dim, new_fft_x_cnts, &
                              new_fft_x_offs, new_fft_y_cnts, new_fft_y_offs) 

  use gbl_constants_mod
  use pmemd_lib_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use pbc_mod
  use processor_mod, only:send_xfft_cnts, send_yfft_cnts, & 
                          send_zfft_cnts, my_xfft_cnt, my_yfft_cnt, &
                          my_zfft_cnt, x_offset, y_offset, z_offset, &
                          proc_dimx_recip, proc_dimy_recip, proc_dimz_recip, &
                          proc_x_comm_lst, &
                          blockunflatten_recip, &
                          proc_dimx_recip, proc_dimy_recip, proc_dimz_recip, &
                          proc_min_nfft1, proc_min_nfft2, proc_min_nfft3
  implicit none


  !Formal arguements 
  integer                 :: fft_x_dim, fft_y_dim, fft_z_dim
  integer                 :: new_fft_x_cnts(0:proc_dimx_recip*proc_dimy_recip-1)
  integer                 :: new_fft_x_offs(0:proc_dimx_recip*proc_dimy_recip-1)
  integer                 :: new_fft_y_cnts(0:proc_dimx_recip*proc_dimy_recip-1)
  integer                 :: new_fft_y_offs(0:proc_dimx_recip*proc_dimy_recip-1)

  !local arguements 
  integer                 :: i,j,k
  integer                 :: max_buf_size
  integer                 :: quot, rem, cnts_idx
  integer                 :: max_fft_y_cnts1, max_fft_y_cnts2 
  integer                 :: max_fft_x_cnts, max_fft_z_cnts 
  integer                 :: my_x, my_y, my_z, grid_idx1, grid_idx2 
  integer                 :: taskmap_idx
  integer                 :: ofset

!  print *, " pme_fft_midpoint proc_dim*",proc_dimx_recip, proc_dimy_recip,proc_dimz_recip,mytaskid_recip

  fft_blks_dim1 = proc_dimx_recip*proc_dimy_recip 
  fft_blks_dim2 = proc_dimz_recip

  if (allocated(fft_task_grid)) then

    deallocate(fft_x_offs, &
               fft_x_cnts, &
               fft_y_offs1, &
               fft_y_offs2, &
               fft_y_cnts1, &
               fft_y_cnts2, &
               fft_z_offs, &
               fft_z_cnts, &
               fft_task_grid, &
               fft_xy_idxlst, &
               fft_yz_idxlst)

  end if

  allocate(fft_x_offs(fft_blks_dim1), &
           fft_x_cnts(fft_blks_dim1), &
           fft_y_offs1(fft_blks_dim1), &
           fft_y_cnts1(fft_blks_dim1), &
           fft_y_offs2(fft_blks_dim2), &
           fft_y_cnts2(fft_blks_dim2), &
           fft_z_offs(fft_blks_dim2), &
           fft_z_cnts(fft_blks_dim2), &
           fft_task_grid(fft_blks_dim1, fft_blks_dim2), &
           fft_xy_idxlst(fft_blks_dim1), &
           fft_yz_idxlst(fft_blks_dim2))
               
!  if (alloc_failed .ne. 0) call setup_alloc_error

  fft_x_offs(:) = 0
  fft_x_cnts(:) = 0
  fft_y_offs1(:) = 0
  fft_y_cnts1(:) = 0
  fft_y_offs2(:) = 0
  fft_y_cnts2(:) = 0
  fft_z_offs(:) = 0
  fft_z_cnts(:) = 0
  fft_task_grid(:,:) = 0
  fft_xy_idxlst(:) = 0
  fft_yz_idxlst(:) = 0

!if(mytaskid_recip .eq. 0) print *, "fft_x_dim",fft_x_dim, "fft_y_dim",fft_y_dim,"fft_z_dim",fft_z_dim

  ! -- first x...
!divide the fft_x_dim into the x dim MPI ranks = fft_blk_dim1
  quot = fft_x_dim / fft_blks_dim1
  rem = fft_x_dim - quot * fft_blks_dim1

  fft_x_cnts(:) = quot

  if (rem .ne. 0) then
    do i = 1, rem
      cnts_idx = fft_blks_dim1 * i / rem
      fft_x_cnts(cnts_idx) = quot + 1
    end do
    max_fft_x_cnts = quot + 1
  else
    max_fft_x_cnts = quot
  end if

  offset = 0
  do i = 1, fft_blks_dim1
    fft_x_offs(i) = offset
    offset = offset + fft_x_cnts(i)
  end do

  ! -- then y in dim 1...

!divide the fft_y_dim into the y dim MPI ranks = x dim MPI rank = fft_blk_dim1
 !since already we developed it in pme_recip_midpoint setup code we just copy
 !here the way old code divide y_dim in ffy_y_cnts1 and new code divides in
 !new_fft_y_cnts are little bit different, so we copy the new to the old data
 !structure
  ! we can actually build new_fft_y_cnts() here instead of pme_recip_midpoint
  ! code, then we can avoid copies  
  do i = 0,proc_dimx_recip*proc_dimy_recip - 1   ! proc_dimx_recip*proc_dimy_recip = fft_blks_dim1
    fft_y_cnts1(i+1) = new_fft_y_cnts(i)
    !fft_y_offs1(i+1) = new_fft_y_offs(i)
  end do
  offset = 0
  do i = 0,proc_dimx_recip*proc_dimy_recip - 1   ! proc_dimx_recip*proc_dimy_recip = fft_blks_dim1
    fft_y_offs1(i+1) = offset !new_fft_y_offs(i)
    offset = offset +  new_fft_y_cnts(i)
  end do
  max_fft_y_cnts1 = maxval(new_fft_y_cnts(0:proc_dimx_recip*proc_dimy_recip-1))  
  ! -- then y in dim 2...
!divide the fft_y_dim into the other y dim MPI ranks = z dim MPI rank = fft_blk_dim2

  quot = fft_y_dim / fft_blks_dim2
  rem = fft_y_dim - quot * fft_blks_dim2

  fft_y_cnts2(:) = quot

  if (rem .ne. 0) then
    do i = 1, rem
      cnts_idx = fft_blks_dim2 * i / rem
      fft_y_cnts2(cnts_idx) = quot + 1
    end do
    max_fft_y_cnts2 = quot + 1
  else
    max_fft_y_cnts2 = quot
  end if

  offset = 0
  do i = 1, fft_blks_dim2
    fft_y_offs2(i) = offset
    offset = offset + fft_y_cnts2(i)
  end do

  ! -- then z...
!divide the fft_z_dim into the z dim MPI ranks = fft_blk_dim2
!since we already divided it in the new code in processor.F90 , we
!just copy it here

  do i = 0,proc_dimz_recip - 1 ! proc_dimz_recip = fft_blks_dim2
    fft_z_cnts(i+1) = send_zfft_cnts(i)
    fft_z_offs(i+1) = z_offset(i)
  end do
  max_fft_z_cnts = maxval(fft_z_cnts(1:proc_dimz_recip))  
!print *,"max_fft_z_cnts",max_fft_z_cnts,fft_z_cnts(1:proc_dimz_recip),mytaskid_recip
!my_grid_idx1 and my_grid_idx2
  call blockunflatten_recip(mytaskid_recip, my_x,my_y,my_z) !input is mytaskid_recip, output is 3d co-ord
  my_grid_idx1 = my_x + my_y * proc_dimx_recip + 1 !2d (x,y) converted to 1d; +1 is
                                             !to move the range, 1 to the right
  my_grid_idx2 = my_z + 1 ! +1 is to move the range, 1 to the right

!creating the fft_task_grid for 2d grid, 1 dim is proc_dimx_recip*proc_dimy_recip long
! 2nd dim is proc_dimz_recip long
  do k = 1, proc_dimz_recip
    do j = 1, proc_dimy_recip
      do i = 1, proc_dimx_recip
        grid_idx1 = i + (j-1) * proc_dimx_recip !x and y flattened to 1d
        grid_idx2 = k
        fft_task_grid(grid_idx1, grid_idx2) = i + (j-1)*proc_dimx_recip + (k-1)*proc_dimx_recip*proc_dimy_recip -1 ! -1 to move range, 1 to the left
      end do
    end do
  end do
  !print *, "from pme_fft_midpoint",fft_task_grid, my_grid_idx1, my_grid_idx2,mytaskid_recip
!my_grid_idx1 and 2 start from 1, but fft_task_grid contains value from 0

!fft_xy_idxlst() stores the other MPI's (x,y) domain flattened to 1d
!,for the same z, basically it is list of 1:fft_blks_dim1 except my_grid_idx1
    taskmap_idx = my_grid_idx1 + 1
    do i = 1, fft_blks_dim1 - 1
      if (taskmap_idx .gt. fft_blks_dim1) taskmap_idx = 1
      fft_xy_idxlst(i) = taskmap_idx
      taskmap_idx = taskmap_idx + 1
    end do

!fft_yz_idxlst() stores the other MPI's (y,z) domain flattened to 1d
!,for the same x, basically it is list of 1:fft_blks_dim2 except my_grid_idx2
    taskmap_idx = my_grid_idx2 + 1
    do i = 1, fft_blks_dim2 - 1
      if (taskmap_idx .gt. fft_blks_dim2) taskmap_idx = 1
      fft_yz_idxlst(i) = taskmap_idx
      taskmap_idx = taskmap_idx + 1
    end do

    ! Allocate space for multidimensional fft transposes.
  
    siz_fft_mpi_buf1 = 2 * max_fft_x_cnts * max_fft_y_cnts1 * &
                            max_fft_z_cnts
    siz_fft_mpi_buf2 = 2 * max_fft_x_cnts * max_fft_y_cnts2 * &
                            max_fft_z_cnts

    max_buf_size = max0(siz_fft_mpi_buf1*fft_blks_dim1,siz_fft_mpi_buf2*fft_blks_dim2)
!later we can have just one set of buffer for both non-bonded code and
!reciprocal code MPi trasfer
! I used another set of buffers for x dim data transfer before the fft start
!we can combine all of the to only one set of buffers
!print *, "blk fft",max_fft_x_cnts,max_fft_y_cnts1,max_fft_z_cnts, mytaskid_recip
  !print *, "max_fft_y_cnt1",max_fft_y_cnts1,max_fft_y_cnts2,max_fft_z_cnts,max_buf_size,mytaskid_recip
  !print *, "proc_dimx_recip", proc_dimx_recip, "proc_dimy_recip", proc_dimy_recip, "proc_dimz_recip", proc_dimz_recip,"fft_blks_dim1" , fft_blks_dim1," fft_blks_dim2", fft_blks_dim2, mytaskid_recip
  if(allocated(proc_send_buf)) deallocate(proc_send_buf)
  if(allocated(proc_recv_buf)) deallocate(proc_recv_buf)
  allocate (proc_send_buf(max_buf_size))
  allocate (proc_recv_buf(max_buf_size))

!print *, "max buf size",siz_fft_mpi_buf1,siz_fft_mpi_buf2,fft_blks_dim1-1,fft_blks_dim2-1,max_buf_size,mytaskid_recip

  !call set_minimum_mpi_bufs_size(siz_fft_mpi_buf1 * (fft_blks_dim1 - 1))
  !call set_minimum_mpi_bufs_size(siz_fft_mpi_buf2 * (fft_blks_dim2 - 1))
#if 0
  if(mytaskid_recip .eq. 0) print *, "fft xim", fft_x_dim, fft_y_dim, fft_z_dim
  if(mytaskid_recip .eq. 0) print *, "latest cnts",fft_x_cnts,"y1 cnts",fft_y_cnts1,"y2 cnts",fft_y_cnts2,"z cnts",fft_z_cnts 
  if(mytaskid_recip .eq. 0) print *, "laters offs",fft_x_offs,"y1 off",fft_y_offs1,"y2 offs",fft_y_offs2,"z off",fft_z_offs 
 ! print *, "sum cnt", sum(fft_x_cnts), sum(fft_y_cnts1), sum(fft_y_cnts2),sum(fft_z_cnts),mytaskid_recip, "buf size ", max_buf_size
#endif
  return

end subroutine blk_fft_setup_midpoint 

!*******************************************************************************
!
! Subroutine:  fft3d_forward_midpoint
!
! Description:  <TBS>
!
! INPUT:
!
! OUTPUT:
!
! frc:      Forces incremented by k-space sum.
! virial:       Virial due to k-space sum (valid for atomic scaling;
!               rigid molecule virial needs a correction term not computed here.
!*******************************************************************************

subroutine fft3d_forward_midpoint(xyz_data, yxz_data, zxy_data,new_fft_y_cnts,new_fft_y_offs,&
                                x_dim, y_dim, z_dim,new_fft_x_cnts, & 
                                    new_fft_x_offs)

  use processor_mod, only: fft_halo_size, send_xfft_cnts, send_yfft_cnts, & 
                          send_zfft_cnts, my_xfft_cnt, my_yfft_cnt, &
                          my_zfft_cnt, x_offset, y_offset, z_offset, &
                          proc_dimx_recip, proc_dimy_recip, proc_dimz_recip, &
                          proc_x_comm_lst, &
                          blockunflatten_recip, &
                          proc_dimx_recip, proc_dimy_recip, proc_dimz_recip, &
                          proc_min_nfft1, proc_min_nfft2, proc_min_nfft3, &
                          proc_max_nfft1, proc_max_nfft2, proc_max_nfft3
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use pbc_mod
  use prmtop_dat_mod
  use timers_mod

  implicit none
  !Formal arguments
  integer               :: x_dim, y_dim, z_dim
#ifdef pmemd_SPDP
  real      :: xyz_data(-1*fft_halo_size:2 * x_dim+fft_halo_size+1, &
                                -1*fft_halo_size:proc_max_nfft2-proc_min_nfft2+fft_halo_size+1, &
                                -1*fft_halo_size:proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)
  real       :: yxz_data(2 * y_dim, &
                                fft_x_cnts(my_grid_idx1), &
                                fft_z_cnts(my_grid_idx2))
                                !proc_max_nfft1-proc_min_nfft1+1, &
                                !proc_max_nfft3-proc_min_nfft3+1)
  real      :: zxy_data(2 * z_dim, &
                                fft_x_cnts(my_grid_idx1), &
                                fft_y_cnts2(my_grid_idx2))
                                ! proc_max_nfft1-proc_min_nfft1+1, &
                                ! proc_max_nfft2-proc_min_nfft2+1)
#else
  double precision      :: xyz_data(-1*fft_halo_size:2 * x_dim+fft_halo_size+1, &
                                -1*fft_halo_size:proc_max_nfft2-proc_min_nfft2+fft_halo_size+1, &
                                -1*fft_halo_size:proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)

  double precision      :: yxz_data(2 * y_dim, &
                                fft_x_cnts(my_grid_idx1), &
                                fft_z_cnts(my_grid_idx2))
                                !proc_max_nfft1-proc_min_nfft1+1, &
                                !proc_max_nfft3-proc_min_nfft3+1)
  double precision      :: zxy_data(2 * z_dim, &
                                fft_x_cnts(my_grid_idx1), &
                                fft_y_cnts2(my_grid_idx2))
                                ! proc_max_nfft1-proc_min_nfft1+1, &
                                ! proc_max_nfft2-proc_min_nfft2+1)
#endif
  !double precision      :: xyz_data(2 * (nfft1+1)/2, my_yfft_cnt, my_zfft_cnt) !change it from charlie
  integer               :: new_fft_y_cnts(0:proc_dimx_recip*proc_dimy_recip-1)
  integer               :: new_fft_y_offs(0:proc_dimx_recip*proc_dimy_recip-1)
  integer               :: new_fft_x_cnts(0:proc_dimx_recip*proc_dimy_recip-1) !division of x fft based on fft_x_dim
  integer               :: new_fft_x_offs(0:proc_dimx_recip*proc_dimy_recip-1)!offset of x fft based on fft_x_dim
!  double precision      :: send_bufy(*), recv_bufy(*)
!  double precision      :: send_bufz(*), recv_bufz(*)

  !local arguments
  integer               :: send_x_offset(proc_dimx_recip) 
  integer               :: send_y_offset(proc_dimy_recip) 
  integer               :: seond_z_offset(proc_dimz_recip) 

  integer               :: x_idx
  integer               :: other_offset, other_mpi_x, other_cnt_x,&
                          other_cnt_y, other_cnt_z, other_cnt
  integer               :: my_cnt_x, my_cnt_y, my_cnt_z
  integer               :: x_tag=13
  integer               :: x,y,z !3d co-ordinate of MPI ranks
  integer               :: my_x,my_y,my_z,my_xy_index !3d co-ordinate of my MPI ranks
  integer               :: my_y_offset, other_y_offset
  
  !this is temporary, this one should come from pme_recip_midpoint setup routine
  integer               :: send_cnt, recv_cnt
  integer               :: i,j,k,buf_base, buf_off 
  integer               :: other_x,other_y, other_z 
  !MPI data structure
  integer               :: irecv_stat(mpi_status_size)
  integer               :: isend_stat(mpi_status_size, proc_dimx_recip - 1)
  integer               :: recv_req(proc_dimx_recip - 1)
  integer               :: send_req(proc_dimx_recip - 1)
  integer               :: cnt_idx

  !offset setting
  send_x_offset(:) = 1

  !alltoall in x direction
 ! proc_min_nfft1, proc_min_nfft2, and proc_min_nfft3 are the offset that
 ! my MPI's "block of fft grid" start

  !post the non-blocking recv
  !I will receive the chunk of other_x_cnt * my_y_cnt * my_z_cnt
  !my_y_cnt is smaller than my_yfft_cnt, and should come from new_fft_y_cnts
  !at least two dim of fft chunk(here my_y_cnt and  my_z_cnt) that I receive should be same as mine, 
  !otherwise I cannot align them and cannot have a full pencil
  !since we are communicating in x direction, x_dim can be different
  call blockunflatten_recip(mytaskid_recip,my_x,my_y,my_z)
  my_cnt_y = new_fft_y_cnts(my_x + proc_dimx_recip*my_y) ! this is how the array was
                                                   ! created, I will also
                                                   ! receieve this amount 
  my_cnt_z = my_zfft_cnt        !Calcuated in processor, I will receieve this amount 
!************Pack Data starts*********
#ifdef separate
  my_cnt_x = my_xfft_cnt !also =  send_xfft_cnt(my_x) 
  my_xy_index = my_y * proc_dimx_recip !this is the index offset from where my_yfft_cnt
  buf_base = 1
  do x_idx = 0, proc_dimx_recip - 1 ! make sure it is right.. only send to MPIs in x direction, to get a full x_dim   
  
    !pack the data
    !other_mpi_x  = proc_x_comm_lst(x_idx)  ! convert to real MPI rank
    if(x_idx .eq. my_x) cycle !avoid sending to me
    !send_cnt = 0 ! proc_min_nfft1 !offset to start my block in x
    other_cnt_y =    new_fft_y_cnts(x_idx+my_xy_index)
    other_y_offset = new_fft_y_offs(x_idx+my_xy_index)
    buf_off = buf_base - 1
    do k  = 1, my_cnt_z !my_zfft_cnt
      do j = other_y_offset + 1, other_y_offset + other_cnt_y !my_yfft_cnt/proc_dimx_recip
                          !the offset is local only to divide my_yfft_cnt, and get the proper chunk 
        proc_send_buf(buf_off+1:buf_off+my_cnt_x) = xyz_data(proc_min_nfft1:proc_min_nfft1+my_cnt_x-1,j,k) 
                                 !also x_offset can be used
                                 !proc_min_fft1 is my offset of the total nfft1
                                  !since for x, full dim is allocated, and my
                                  !chunk is properly spaced with proc_min_nfft1
                                  !offset, that is why we need proc_min_nfft1
        buf_off = buf_off + my_cnt_x 
      end do
    end do
    buf_base = buf_base + my_cnt_x * other_cnt_y * my_cnt_z ! these are sent in one MPI_isend 
  end do
#endif
!************Pack Data Done *********
#ifdef fft_time
  call update_pme_time(data_pack_timer)
#endif

  other_offset = 0
  cnt_idx = 1
  do x_idx = 0, proc_dimx_recip - 1 ! only send to MPIs in x direction   

    !post the asynchronous receives first
    other_mpi_x  = proc_x_comm_lst(x_idx)  ! real MPI rank in x direction
    if(other_mpi_x .eq. mytaskid_recip) cycle !avoid receiving from me
    !call blockunflatten_recip(other_mpi_x,other_x,other_y,other_z)
    other_cnt_x = send_xfft_cnts(x_idx) !other_x) 
    other_cnt   = other_cnt_x * my_cnt_y * my_cnt_z ! we dont need to use *2 here 

#ifdef pmemd_SPDP
    call mpi_irecv(proc_recv_buf(other_offset+1), other_cnt, mpi_real, &
#else
    call mpi_irecv(proc_recv_buf(other_offset+1), other_cnt, mpi_double_precision, &
#endif
                   other_mpi_x, x_tag, pmemd_comm, &
                   recv_req(cnt_idx), err_code_mpi); MPI_ERR_CHK("in fft ")
    other_offset = other_offset + other_cnt 
    cnt_idx = cnt_idx+1
  end do

  !pack data and  post the non-blocking send
  !I will send the chunk of my_x_cnt * other_y_cnt * (my=other)z_cnt
  !other_y_cnt is smaller than my_yfft_cnt, and should come from new_fft_y_cnts
  my_cnt_x = my_xfft_cnt !also =  send_xfft_cnt(my_x) 
  my_xy_index = my_y * proc_dimx_recip !this is the index offset from where my_yfft_cnt
                                 !starts in full arrays new_fft_y_cnts() and new_fft_y_offs()
  buf_base = 1
  cnt_idx  = 1
  do x_idx = 0, proc_dimx_recip - 1 ! make sure it is right.. only send to MPIs in x direction, to get a full x_dim   
  
    !pack the data
    other_mpi_x  = proc_x_comm_lst(x_idx)  ! convert to real MPI rank
    if(other_mpi_x .eq. mytaskid_recip) cycle !avoid sending to me
    !send_cnt = 0 ! proc_min_nfft1 !offset to start my block in x
    other_cnt_y =    new_fft_y_cnts(x_idx+my_xy_index)!x_idx=other_x
    other_y_offset = new_fft_y_offs(x_idx+my_xy_index)
    buf_off = buf_base - 1
    !if(buf_base+my_cnt_x * other_cnt_y * my_cnt_z .gt. size(proc_send_buf)) 
!   if (mytaskid_recip .eq. 1 .or. mytaskid_recip .eq. 0) print *, "mpi send",my_cnt_x * other_cnt_y *my_cnt_z,size(proc_send_buf),other_mpi_x,mytaskid_recip
#ifndef separate
    do k  = 1, my_cnt_z !my_zfft_cnt
      do j = other_y_offset + 1, other_y_offset + other_cnt_y !my_yfft_cnt/proc_dimx_recip
                          !the offset is local only to divide my_yfft_cnt, and get the proper chunk 
        proc_send_buf(buf_off+1:buf_off+my_cnt_x) = xyz_data(proc_min_nfft1:proc_min_nfft1+my_cnt_x-1,j,k) 
                                 !also x_offset can be used
                                 !proc_min_fft1 is my offset of the total nfft1
                                  !since for x, full dim is allocated, and my
                                  !chunk is properly spaced with proc_min_nfft1
                                  !offset, that is why we need proc_min_nfft1
        buf_off = buf_off + my_cnt_x !send_cnt is offset for send_bufx
      end do
    end do
#endif
    !post the non-blocking send
#ifdef separate
    send_cnt = my_cnt_x * other_cnt_y * my_cnt_z
    call mpi_isend(proc_send_buf(buf_base), send_cnt, & !buf_off - buf_base + 1, &
#else
    call mpi_isend(proc_send_buf(buf_base), buf_off - buf_base + 1, &
#endif
#ifdef pmemd_SPDP
                   mpi_real, other_mpi_x, x_tag, pmemd_comm, &
#else
                   mpi_double_precision, other_mpi_x, x_tag, pmemd_comm, &
#endif
                   send_req(cnt_idx), err_code_mpi)
    buf_base = buf_base + my_cnt_x * other_cnt_y * my_cnt_z ! these many are sent in one MPI_isend 
    cnt_idx = cnt_idx + 1
  end do
  ! Wait for all the MPI receive to complete
  do x_idx = 1, proc_dimx_recip - 1 ! receiving from proc_dimx_recip-1 ranks
    call mpi_waitany(proc_dimx_recip - 1, recv_req, cnt_idx, irecv_stat, &
                     err_code_mpi)!instead of cnt_idx any other variable can be used, but dont use x_idx, that will change the index of the loop,a big bug showed up when x_idx was used here in palce of cnt_idx 
  end do

  ! Wait for all MPI sends to complete:
  call mpi_waitall(proc_dimx_recip - 1, send_req, isend_stat, err_code_mpi)

#ifdef fft_time
  call update_pme_time(mpi_forward_timer)
#endif

  !unpack the recv data
  recv_cnt = 0
  my_y_offset = new_fft_y_offs(my_x+my_xy_index)
                                !To place y into proper place in xyz_data
                                !since we divided my_yfft-cnt into more
                                !proc_dimx_recip divisions
                                !we could build new_fft_y_offset() as
                                !0:proc_dimx_recip-1 long, then my_xy_index would
                                !not be required as index for new_fft_y_offset
  do x_idx = 0, proc_dimx_recip - 1
    !other_mpi_x = proc_x_comm_lst(x_idx) 
    if(x_idx .eq. my_x) cycle !avoid receiving from  me
    !call blockunflatten_recip(other_mpi_x,other_x,other_y,other_z)!getting other x,y,z
    other_offset = x_offset(x_idx) !other_x) !x_idx can be used too; at this location start to store
    other_cnt_x = send_xfft_cnts(x_idx) !other_x) 
    do k = 1, my_cnt_z
      do j = my_y_offset + 1, my_y_offset + my_cnt_y ! though it was previously in another place of y, now we all moved it to start from 1 
      !do j = 1, my_cnt_y ! though it was previously in another place of y, now we all moved it to start from 1 
        !do i=1, other_cnt_x
          xyz_data(other_offset+1:other_offset+other_cnt_x,j,k) = proc_recv_buf(recv_cnt+1:recv_cnt+other_cnt_x)
        !end do
        recv_cnt = recv_cnt + other_cnt_x
      end do
    end do
  end do
  !if my y does not start from 1, then move it to 1, we are doing it so that it
  !will match the position of y with old and new code
#if 0
  if(my_y_offset .gt. 0) then ! it means it starts not from 1, so copying and pasting down at 1 
    do k = 1, my_cnt_z
      do j = 1, my_cnt_y 
        do i = proc_min_nfft1, proc_min_nfft1+my_cnt_x-1
          xyz_data(i,j,k) = xyz_data(i,my_y_offset+j,k) 
        end do
      end do
    end do

  end if
#endif 
#if 0
  do k=1, fft_z_cnts(my_grid_idx2)
    do j=1, fft_y_cnts1(my_grid_idx1)
      do i=1, 2*x_dim
          write(15+mytaskid_recip,*) "new bef x",i, fft_y_offs1(my_grid_idx1)+j,fft_z_offs(my_grid_idx2)+k,xyz_data(i,j,k) 
      end do
     end do
   end do
#endif
!put ifdef fft_time for each timer
#ifdef fft_time
  call update_pme_time(data_unpack_timer)
#endif

  call do_forward_x_fft(x_dim, xyz_data,my_y_offset)
!  call do_forward_x_fft(x_dim, xyz_data)
#ifdef fft_time
  call update_pme_time(forward_x_timer)
#endif

  call do_xy_transpose(x_dim, y_dim, xyz_data, yxz_data, &
                       proc_send_buf, proc_recv_buf,my_y_offset)
#ifdef fft_time
  call update_pme_time(xy_transpose_timer)
#endif

  call do_forward_y_fft(y_dim, yxz_data)
#ifdef fft_time
  call update_pme_time(forward_y_timer)
#endif

  call do_yz_transpose(y_dim, z_dim, yxz_data, zxy_data, &
                       proc_send_buf, proc_recv_buf)
#ifdef fft_time
  call update_pme_time(yz_transpose_timer)
#endif

  call do_forward_z_fft(z_dim, zxy_data)
#ifdef fft_time
  call update_pme_time(forward_z_timer)
#endif

  return

end subroutine fft3d_forward_midpoint
!*******************************************************************************
!
! Subroutine:  do_forward_x_fft
!
! Description:  <TBS>
!
!*******************************************************************************

!subroutine do_forward_x_fft(x_dim, x_runs)
subroutine do_forward_x_fft(x_dim, x_runs, my_y_offset)
  use processor_mod, only: fft_halo_size, myblockid, proc_max_nfft2,proc_min_nfft2, &
                                 proc_max_nfft3,proc_min_nfft3
  implicit none

! Formal arguments:

  integer               :: x_dim
#ifdef pmemd_SPDP
  real      :: x_runs(2, -1*fft_halo_size/2:x_dim+fft_halo_size/2+1, & !1st time it is divided into 2
                                 ! fft_y_cnts1(my_grid_idx1), &
                                 ! fft_z_cnts(my_grid_idx2))
                                -1*fft_halo_size: proc_max_nfft2-proc_min_nfft2+fft_halo_size+1, &
                                -1*fft_halo_size: proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)
! Local variables:

  real      :: a, b, c, d
  real      :: buf(2, x_dim - 1)
#else
  double precision      :: x_runs(2, -1*fft_halo_size/2:x_dim+fft_halo_size/2+1, & !1st time it is divided into 2
                                 ! fft_y_cnts1(my_grid_idx1), &
                                 ! fft_z_cnts(my_grid_idx2))
                                -1*fft_halo_size: proc_max_nfft2-proc_min_nfft2+fft_halo_size+1, &
                                -1*fft_halo_size: proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)
! Local variables:

  double precision      :: a, b, c, d
  double precision      :: buf(2, x_dim - 1)
#endif
  integer               :: i, j, k
  integer               :: x_lim        ! indexing limit
  integer               :: m
  integer               :: my_y_offset
! real to complex:
#if 0
  do k=1, fft_z_cnts(my_grid_idx2)
    do j=1, fft_y_cnts1(my_grid_idx1)
      do i=1, x_dim
       do m = 1,2
          write(15+myblockid,*) "new",2*(i-1)+m,fft_y_offs1(my_grid_idx1)+j,fft_z_offs(my_grid_idx2)+k,x_runs(m,i,j,k)                                            
       end do
      end do
     end do
   end do
#endif

  x_lim = x_dim + 1

  do k = 1, fft_z_cnts(my_grid_idx2)

    do j = my_y_offset+1, my_y_offset + fft_y_cnts1(my_grid_idx1)
    !do j = 1, fft_y_cnts1(my_grid_idx1)

      buf(:,:) = x_runs(:, 1:x_dim - 1, j, k)
      call fft1d_forward(fft_x_hdl, buf)
      do i = 2, x_dim - 1
        a =  (buf(1, i) + buf(1, x_lim - i)) ! Real F even * 2
        b =  (buf(2, i) - buf(2, x_lim - i)) ! Imag F even * 2
        c =  (buf(2, i) + buf(2, x_lim - i)) ! Real F odd * 2
        d = -(buf(1, i) - buf(1, x_lim - i)) ! Imag F odd * 2
        
        x_runs(1, i, j, k) = .5d0 * (a + fftrc_coefs(1, i) * c + &
                                     fftrc_coefs(2, i) * d)
        x_runs(2, i, j, k) = .5d0 * (b + fftrc_coefs(1, i) * d - &
                                     fftrc_coefs(2, i) * c)
      end do

      ! DC and nyquist:

      x_runs(1, 1, j, k) = buf(1, 1) + buf(2, 1)
      x_runs(2, 1, j, k) = 0.d0
      x_runs(1, x_dim, j, k) = buf(1, 1) - buf(2, 1)
      x_runs(2, x_dim, j, k) = 0.d0

    end do

  end do
#if 0
  do k=1, fft_z_cnts(my_grid_idx2)
    do j=1, fft_y_cnts1(my_grid_idx1)
      do i=1, x_dim
       do m = 1,2
          write(15+myblockid,*) "new aft x",2*(i-1)+m,fft_y_offs1(my_grid_idx1)+j,fft_z_offs(my_grid_idx2)+k,x_runs(m,i,j,k)                                            
       end do
      end do
     end do
   end do
#endif

  return
      
end subroutine do_forward_x_fft

!*******************************************************************************
!
! Subroutine:  do_xy_transpose
!
! Description: Fully asynchronous implementation!
!              
!*******************************************************************************

!subroutine do_xy_transpose(x_dim, y_dim, xyz_data, yxz_data, send_buf, recv_buf)
subroutine do_xy_transpose(x_dim, y_dim, xyz_data, yxz_data, send_buf,recv_buf,my_y_offset)

  use parallel_dat_mod
  use processor_mod, only: fft_halo_size, myblockid, proc_max_nfft2,proc_min_nfft2, &
                                 proc_max_nfft3,proc_min_nfft3

  implicit none

! Formal arguments:

  integer               :: x_dim
  integer               :: y_dim
#ifdef pmemd_SPDP
  real      :: xyz_data(-1*fft_halo_size:2 * x_dim+fft_halo_size+1, &
                                -1*fft_halo_size: proc_max_nfft2-proc_min_nfft2+fft_halo_size+1, &
                                -1*fft_halo_size: proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)
  real      :: yxz_data(2 * y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
  real      :: send_buf(siz_fft_mpi_buf1 * fft_blks_dim1)
  real      :: recv_buf(siz_fft_mpi_buf1, fft_blks_dim1)
#else
  double precision      :: xyz_data(-1*fft_halo_size:2 * x_dim+fft_halo_size+1, &
                                -1*fft_halo_size: proc_max_nfft2-proc_min_nfft2+fft_halo_size+1, &
                                -1*fft_halo_size: proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)
  double precision      :: yxz_data(2 * y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
  double precision      :: send_buf(siz_fft_mpi_buf1 * fft_blks_dim1)
  double precision      :: recv_buf(siz_fft_mpi_buf1, fft_blks_dim1)
#endif

! Local variables:

  integer               :: i, j, k
  integer               :: numval, buf_base, buf_off
  integer               :: xy_idxlst_idx, recv_task, send_task
  integer               :: other_grid_idx1
  integer               :: my_x_cnt
  integer               :: my_y_cnt
  integer               :: other_x_off, other_x_cnt
  integer               :: other_y_cnt
  integer               :: z_cnt
  integer               :: irecv_stat(mpi_status_size)
  integer               :: isend_stat(mpi_status_size, fft_blks_dim1 - 1)
  integer               :: recv_req(fft_blks_dim1)
  integer               :: send_req(fft_blks_dim1)
  integer               :: xyyxt_tags=199
  integer               :: my_y_offset
#ifdef COMM_TIME_TEST
  call start_test_timer(9, 'do_xy_transpose', 0)
#endif
  !print *, "xy transpose starts",mytaskid_recip 
!the goal of this routine is to line up and create a pencil in y direction
!To do that in y direction,  my and "other MPI sending me" should have same z and x cnt
!but y cnt could be different
!The following send and recv could be confusing a bit since currently all the
!MPI sending and receiving has full x but partial y and z
!For example there are 4 MPI which have full x and partial y and z
!Need to MPI communicate the data, so that in future 4 MPI will have all y and
!partial x and z
!so,each MPI will send other_x, my_y, and my_z to other 3 MPI
!each MPI will receive my_x, other_y, and my_z = other_z from other MPI
!Thus, I am getting 3 chunks, each of my_x * my_z * other_y. Thus adding with my
!current chunck my_x*my_z*my_y with the 3 chunks, I get the full pencil in full y, 
!since given same x and z, I get all 4 chuncks of y, and after glueing them
!together in y dierction, I get a pencit with full y, so there will be 4 pencils
!in y direction owned by 4 MPIs

!on the contrary, if I send my x, my y and my z, then we only get multiple data
!chunk, we dont get all chunk of y, so that dont work (proved in paper also)


! my_grid_idx1 kind of my MPI rank position in x and in x directions,
! using my_grid_idx1, I find the fft cnt I own in both x and y direction
  my_x_cnt = fft_x_cnts(my_grid_idx1)
  my_y_cnt = fft_y_cnts1(my_grid_idx1)

  ! The z data is constant for the xy transpose:

!my_grid_idx1 is my MPI position in x when x_fft is divided by the MPI ranks 
!same my_grid_idx1 is my MPI position in y when y_fft is divided by the MPI ranks
!my_grid_idx2 is my position wnen z_fft is divided by the MPI ranks

  z_cnt = fft_z_cnts(my_grid_idx2) !should be same for me and all the senders
                                   !MPI to me


!fft_blks_dim1 is the total # of divisions in x direction
!also fft_blks_dim1 is the total # of divisions in y direction
!fft_blks_dim2 is the division in z direction
!fft_xy_idxlst() contains list of MPI ranks in x other than me..
!for xy case, meaning list of fft_blks_dim1-1 MPI ranks

  ! Post the asynchronous receives first:
  do xy_idxlst_idx = 1, fft_blks_dim1 - 1
    other_grid_idx1 = fft_xy_idxlst(xy_idxlst_idx)
    other_y_cnt = fft_y_cnts1(other_grid_idx1)
    !Actual MPI rank can be found from fft_task_grid
    ! my_grid_idx2 is used for z, since I will receive from same fft z index
    ! my_x_cnt is also same for other MPI
    ! only other_y_cnt is different
    recv_task = fft_task_grid(other_grid_idx1, my_grid_idx2)

    numval = 2 * my_x_cnt * other_y_cnt * z_cnt ! my_x_cnt and z_cnt should be
!    if(recv_task .eq. 1)print *, numval, my_x_cnt, other_y_cnt, z_cnt,other_grid_idx1 ,other_y_cnt  , recv_task ,xy_idxlst_idx  ,mytaskid_recip
                        !same for me and the mpi I receive data, since only y_cnt is different
! since it is column major, for fixed xy_idxlst_idx, it receives data in the row
! starting from 1
! alwasy MPI receives one dimension data, so 3d is straighten to 1d
! so, all numval values is recieved in 1 dimension
 ! print *, "Receving in x : I am ",mytaskid_recip,"receiving from",recv_task,"cnt ",numval
#ifdef pmemd_SPDP
    call mpi_irecv(recv_buf(1, xy_idxlst_idx), numval, mpi_real, &
#else
    call mpi_irecv(recv_buf(1, xy_idxlst_idx), numval, mpi_double_precision, &
#endif
                   recv_task, xyyxt_tags, pmemd_comm, &
                   recv_req(xy_idxlst_idx), err_code_mpi); MPI_ERR_CHK("in fft ")

  end do
  ! Now set up and post the asynchronous sends:

  buf_base = 1
  do xy_idxlst_idx = 1, fft_blks_dim1 - 1
    other_grid_idx1 = fft_xy_idxlst(xy_idxlst_idx)
    other_x_off = fft_x_offs(other_grid_idx1) * 2 !other_x_off for other mpi and
                                            !my mpi are same,since we both are
                                            !in the same line. same for z. only y are
                                            !different
                      !since my MPI has full x
                      ! dim, I need to know from which offset I should send
                      ! 2 is for complex
    other_x_cnt = fft_x_cnts(other_grid_idx1) * 2 !other x_cnt and my x cnt are
                                                    !same.

    send_task = fft_task_grid(other_grid_idx1, my_grid_idx2)
    buf_off = buf_base - 1

    do k = 1, z_cnt
     ! do j = 1, my_y_cnt
      do j = my_y_offset+1, my_y_offset+my_y_cnt
        ! copy "1:other_x_cnts" data with proper offset from my big array to send_buf for specific y and z
        ! send_buf is single dim.
        ! the way send_buf is filled is 
        ! all x for y,z, then next all x for z,y+1 then next all x for z, y+2 and
        ! so on
        send_buf(buf_off+1:buf_off+other_x_cnt) = &
          xyz_data(other_x_off+1:other_x_off+other_x_cnt, j, k)
        buf_off = buf_off + other_x_cnt !move the offset to accomodate all of my y and z
      end do
    end do

 ! print *, "Sending in x : I am ",mytaskid_recip,"Sending to",send_task,"cnt ",buf_off - buf_base + 1
    call mpi_isend(send_buf(buf_base), buf_off - buf_base + 1, &
#ifdef pmemd_SPDP
                   mpi_real, send_task, xyyxt_tags, pmemd_comm, &
#else
                   mpi_double_precision, send_task, xyyxt_tags, pmemd_comm, &
#endif
                   send_req(xy_idxlst_idx), err_code_mpi)

    buf_base = buf_base + siz_fft_mpi_buf1 !send_buf is used for sending to all
                                            ! fft_blks_dim1-1. Thus, it needs
                                            ! the buf_base since send_buf keeps
                                            ! all value as it is asynchronous
                                            ! and out of order send

  end do

  !print *, "xy transpose ends",mytaskid_recip 
  ! Do your own transposes, overlapped with i/o.
  ! my specific part whch I did not MPI send or receive, transpose x <-> y
  !call xy_transpose_my_data(x_dim, y_dim, xyz_data, yxz_data)
  call xy_transpose_my_data(x_dim, y_dim, xyz_data, yxz_data,my_y_offset)

  ! Wait on and process the pending receive requests:

  do j = 1, fft_blks_dim1 - 1

    call mpi_waitany(fft_blks_dim1 - 1, recv_req, xy_idxlst_idx, irecv_stat, &
                     err_code_mpi)
    other_grid_idx1 = fft_xy_idxlst(xy_idxlst_idx)
    recv_task = fft_task_grid(other_grid_idx1, my_grid_idx2)
    !transposing after receive, x->y and y->x
    call do_xy_transpose_recv(y_dim, other_grid_idx1, &
                              recv_buf(1, xy_idxlst_idx), yxz_data)

  end do

  ! Wait for all sends to complete:

  call mpi_waitall(fft_blks_dim1 - 1, send_req, isend_stat, err_code_mpi)

#if 0
  do k=1, fft_z_cnts(my_grid_idx2)
    do j=1, fft_x_cnts(my_grid_idx1)
      do i=1, 2*y_dim
          write(15+myblockid,*) "new aft xy",fft_x_offs(my_grid_idx1)+j,i,fft_z_offs(my_grid_idx2)+k,yxz_data(i,j,k) 
      end do
     end do
   end do
#endif

#ifdef COMM_TIME_TEST
  call stop_test_timer(9)
#endif

  return

contains

!*******************************************************************************
!
! Internal Subroutine:  do_xy_transpose_recv
!
! Description: Fully asynchronous implementation!
!             !nothing to be asynchronouse 
!*******************************************************************************

subroutine do_xy_transpose_recv(y_dim, other_grid_idx1, recv_buf, yxz_data)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: y_dim
  integer               :: other_grid_idx1
#ifdef pmemd_SPDP
  real      :: recv_buf(fft_x_cnts(my_grid_idx1) * 2 *&
                                    fft_y_cnts1(other_grid_idx1) * &
                                    fft_z_cnts(my_grid_idx2))
  real      :: yxz_data(2, y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
#else
#ifdef fft_opt
  double precision      :: recv_buf(2,fft_x_cnts(my_grid_idx1) *&
                                    fft_y_cnts1(other_grid_idx1) * &
                                    fft_z_cnts(my_grid_idx2))
#else
  double precision      :: recv_buf(fft_x_cnts(my_grid_idx1) * 2 *&
                                    fft_y_cnts1(other_grid_idx1) * &
                                    fft_z_cnts(my_grid_idx2))
#endif
  double precision      :: yxz_data(2, y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
#endif

! Local variables:

  integer               :: i, j, k
  integer               :: buf_idx
  integer               :: my_x_cnt
  integer               :: other_y_off, other_y_cnt
  integer               :: z_cnt

  my_x_cnt = fft_x_cnts(my_grid_idx1)
  other_y_cnt = fft_y_cnts1(other_grid_idx1)
  other_y_off = fft_y_offs1(other_grid_idx1)
  z_cnt = fft_z_cnts(my_grid_idx2)

  buf_idx = 0
  !converting 1d to 3d and also transposing x<->y at the same time
#ifdef fft_opt
!probably we can do some optimization by interchanging i and j and collapsing
!1,2 to 1:2 in the 1st dimension
  do k = 1, z_cnt
    do j = 1, other_y_cnt
      do i = 1, my_x_cnt
        yxz_data(1:2, other_y_off + j, i, k) = recv_buf(1:2,buf_idx) !transposing x to y for real
        buf_idx = buf_idx + 2
      end do
    end do
  end do
#else
  do k = 1, z_cnt
    do j = 1, other_y_cnt
      do i = 1, my_x_cnt
        buf_idx = buf_idx + 1
        yxz_data(1, other_y_off + j, i, k) = recv_buf(buf_idx) !transposing x to y for real
        buf_idx = buf_idx + 1
        yxz_data(2, other_y_off + j, i, k) = recv_buf(buf_idx) !transposing x to y for imaginary
      end do
    end do
  end do
#endif
 
  return

end subroutine do_xy_transpose_recv

!*******************************************************************************
!
! Internal Subroutine:  xy_transpose_my_data
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine xy_transpose_my_data(x_dim, y_dim, xyz_data, yxz_data,my_y_offset)
!subroutine xy_transpose_my_data(x_dim, y_dim, xyz_data, yxz_data)

  use parallel_dat_mod
  use processor_mod, only: fft_halo_size, myblockid, proc_max_nfft2,proc_min_nfft2, &
                                 proc_max_nfft3,proc_min_nfft3

  implicit none

! Formal arguments:

  integer               :: x_dim
  integer               :: y_dim
#ifdef pmemd_SPDP
  real      :: xyz_data(2, -1*fft_halo_size/2:x_dim+fft_halo_size/2+1, &
                                -1*fft_halo_size: proc_max_nfft2-proc_min_nfft2+fft_halo_size+1, &
                                -1*fft_halo_size: proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)
                                 !   fft_y_cnts1(my_grid_idx1), &
                                 !   fft_z_cnts(my_grid_idx2))
  real      :: yxz_data(2, y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
#else
  double precision      :: xyz_data(2, -1*fft_halo_size/2:x_dim+fft_halo_size/2+1, &
                                -1*fft_halo_size: proc_max_nfft2-proc_min_nfft2+fft_halo_size+1, &
                                -1*fft_halo_size: proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)
                                 !   fft_y_cnts1(my_grid_idx1), &
                                 !   fft_z_cnts(my_grid_idx2))
  double precision      :: yxz_data(2, y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
#endif

! Local variables:

  integer               :: i, j, k
  integer               :: my_x_off, my_x_cnt
  integer               :: my_y_off, my_y_cnt
  integer               :: z_cnt
  integer               :: my_y_offset

  my_x_off = fft_x_offs(my_grid_idx1)
  my_x_cnt = fft_x_cnts(my_grid_idx1)
  my_y_off = fft_y_offs1(my_grid_idx1)
  my_y_cnt = fft_y_cnts1(my_grid_idx1)
  z_cnt = fft_z_cnts(my_grid_idx2)

  do k = 1, z_cnt
    do j = 1, my_y_cnt
      do i = 1, my_x_cnt
        !yxz_data(:, my_y_off + j, i, k) = xyz_data(:, my_x_off + i, j, k)
        yxz_data(:, my_y_off + j, i, k) = xyz_data(:, my_x_off + i, my_y_offset+j, k)
      end do
    end do
  end do
 
  return

end subroutine xy_transpose_my_data
end subroutine do_xy_transpose

!*******************************************************************************
!
! Subroutine:  do_forward_y_fft
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine do_forward_y_fft(y_dim, y_runs)

  use parallel_dat_mod
  use processor_mod, only: fft_halo_size, myblockid

  implicit none

! Formal arguments:

  integer               :: y_dim
#ifdef pmemd_SPDP
  real      :: y_runs(2, y_dim, &
                                  fft_x_cnts(my_grid_idx1), &
                                  fft_z_cnts(my_grid_idx2))
#else
  double precision      :: y_runs(2, y_dim, &
                                  fft_x_cnts(my_grid_idx1), &
                                  fft_z_cnts(my_grid_idx2))
#endif

! Local variables:

  integer               :: j, k, i, m

#if 0
  do k=1, fft_z_cnts(my_grid_idx2)
    do j=1, fft_x_cnts(my_grid_idx1)
      do i=1, y_dim
        do m=1,2
          write(15+myblockid,*) "new bef y",2*(i-1)+m,fft_x_offs(my_grid_idx1)+j,fft_z_offs(my_grid_idx2)+k,y_runs(m,i,j,k) 
        end do
      end do
     end do
   end do
#endif
  do k = 1, fft_z_cnts(my_grid_idx2)
    do j = 1, fft_x_cnts(my_grid_idx1)
      !that is why we needed the transpose, when we send the array to fft
      !routine, we need to have contiguous data in memory to do the fft. For
      !xyz(), x data are contiguous, but y data is skipped every fft_x_dim
      !thus in xyz() format we cannot do fft in y-dim. Thus, xyz() is transposed
      !to yxz() so that in yxz(), y data are now contigous. Thus we can do fft
      !in y-dim
      !major and all of 1st dimension can be sent and received
      call fft1d_forward(fft_y_hdl, y_runs(1, 1, j, k))


    end do
  end do
#if 0
  do k=1, fft_z_cnts(my_grid_idx2)
    do j=1, fft_x_cnts(my_grid_idx1)
      do i=1, y_dim
        do m=1,2
          write(15+myblockid,*) "new aft y",2*(i-1)+m,fft_x_offs(my_grid_idx1)+j,fft_z_offs(my_grid_idx2)+k,y_runs(m,i,j,k) 
        end do
      end do
     end do
   end do
#endif

  return
      
end subroutine do_forward_y_fft

!*******************************************************************************
!
! Subroutine:  do_yz_transpose
!
! Description: Fully asynchronous implementation!
!              
!*******************************************************************************

subroutine do_yz_transpose(y_dim, z_dim, yxz_data, zxy_data, send_buf, recv_buf)

  use parallel_dat_mod
  use processor_mod, only: fft_halo_size, myblockid

  implicit none

! Formal arguments:

  integer               :: y_dim
  integer               :: z_dim
#ifdef pmemd_SPDP
  real      :: yxz_data(2 * y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
  real      :: zxy_data(2 * z_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_y_cnts2(my_grid_idx2))
  real      :: send_buf(siz_fft_mpi_buf2 * fft_blks_dim2)
  real      :: recv_buf(siz_fft_mpi_buf2, fft_blks_dim2)
#else
  double precision      :: yxz_data(2 * y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
  double precision      :: zxy_data(2 * z_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_y_cnts2(my_grid_idx2))
  double precision      :: send_buf(siz_fft_mpi_buf2 * fft_blks_dim2)
  double precision      :: recv_buf(siz_fft_mpi_buf2, fft_blks_dim2)
#endif

! Local variables:

  integer               :: i, j, k
  integer               :: numval, buf_base, buf_off
  integer               :: yz_idxlst_idx, recv_task, send_task
  integer               :: other_grid_idx2
  integer               :: my_y_cnt
  integer               :: my_z_cnt
  integer               :: other_y_off, other_y_cnt
  integer               :: other_z_cnt
  integer               :: x_cnt
  integer               :: irecv_stat(mpi_status_size)
  integer               :: isend_stat(mpi_status_size, fft_blks_dim2 - 1)
  integer               :: recv_req(fft_blks_dim2)
  integer               :: send_req(fft_blks_dim2)

#ifdef COMM_TIME_TEST
  call start_test_timer(10, 'do_yz_transpose', 0)
#endif
#if 0
  do k=1, fft_z_cnts(my_grid_idx2)
    do j=1, fft_x_cnts(my_grid_idx1)
      do i=1, 2*y_dim
          write(15+myblockid,*) "new bef yz",fft_x_offs(my_grid_idx1)+j,i,fft_z_offs(my_grid_idx2)+k,yxz_data(i,j,k) 
      end do
     end do
   end do
#endif

  my_y_cnt = fft_y_cnts2(my_grid_idx2)
  my_z_cnt = fft_z_cnts(my_grid_idx2)

  ! The x data is constant for the yz transpose:

  x_cnt = fft_x_cnts(my_grid_idx1)

  ! Post the asynchronous receives first:

  do yz_idxlst_idx = 1, fft_blks_dim2 - 1
    other_grid_idx2 = fft_yz_idxlst(yz_idxlst_idx)
    other_z_cnt = fft_z_cnts(other_grid_idx2)

    recv_task = fft_task_grid(my_grid_idx1, other_grid_idx2)

    numval = 2 * my_y_cnt * x_cnt * other_z_cnt

#ifdef pmemd_SPDP
    call mpi_irecv(recv_buf(1, yz_idxlst_idx), numval, mpi_real, &
#else
    call mpi_irecv(recv_buf(1, yz_idxlst_idx), numval, mpi_double_precision, &
#endif
                   recv_task, yzzyt_tag, pmemd_comm, &
                   recv_req(yz_idxlst_idx), err_code_mpi); MPI_ERR_CHK("in fft ")
  end do

  ! Now set up and post the asynchronous sends:

  buf_base = 1

  do yz_idxlst_idx = 1, fft_blks_dim2 - 1
    other_grid_idx2 = fft_yz_idxlst(yz_idxlst_idx)
    other_y_off = fft_y_offs2(other_grid_idx2) * 2
    other_y_cnt = fft_y_cnts2(other_grid_idx2) * 2

    send_task = fft_task_grid(my_grid_idx1, other_grid_idx2)
    buf_off = buf_base - 1

    do k = 1, my_z_cnt
      do j = 1, x_cnt
        send_buf(buf_off+1:buf_off+other_y_cnt) = &
          yxz_data(other_y_off+1:other_y_off+other_y_cnt, j, k)
        buf_off = buf_off + other_y_cnt
      end do
    end do

    call mpi_isend(send_buf(buf_base), buf_off - buf_base + 1, &
#ifdef pmemd_SPDP
                   mpi_real, send_task, yzzyt_tag, pmemd_comm, &
#else
                   mpi_double_precision, send_task, yzzyt_tag, pmemd_comm, &
#endif
                   send_req(yz_idxlst_idx), err_code_mpi)

    buf_base = buf_base + siz_fft_mpi_buf2

  end do

  ! Do your own transposes, overlapped with i/o.

  call yz_transpose_my_data(y_dim, z_dim, yxz_data, zxy_data)

  ! Wait on and process the pending receive requests:

  do j = 1, fft_blks_dim2 - 1

    call mpi_waitany(fft_blks_dim2 - 1, recv_req, yz_idxlst_idx, irecv_stat, &
                     err_code_mpi)
    other_grid_idx2 = fft_yz_idxlst(yz_idxlst_idx)
    recv_task = fft_task_grid(my_grid_idx1, other_grid_idx2)
    call do_yz_transpose_recv(z_dim, other_grid_idx2, &
                              recv_buf(1, yz_idxlst_idx), zxy_data)

  end do

  ! Wait for all sends to complete:

  call mpi_waitall(fft_blks_dim2 - 1, send_req, isend_stat, err_code_mpi)

#if 0
  do k=1, fft_y_cnts2(my_grid_idx2)
    do j=1, fft_x_cnts(my_grid_idx1)
      do i=1, 2*z_dim
          write(15+myblockid,*) "new aft yz",fft_x_offs(my_grid_idx1)+j,fft_y_offs2(my_grid_idx2)+k,i,zxy_data(i,j,k) 
      end do
     end do
   end do
#endif
#ifdef COMM_TIME_TEST
  call stop_test_timer(10)
#endif

  return

contains

!*******************************************************************************
!
! Internal Subroutine:  do_yz_transpose_recv
!
! Description: Fully asynchronous implementation!
!              
!*******************************************************************************

subroutine do_yz_transpose_recv(z_dim, other_grid_idx2, recv_buf, zxy_data)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: z_dim
  integer               :: other_grid_idx2
#ifdef pmemd_SPDP
  real      :: recv_buf(2 * fft_y_cnts2(my_grid_idx2) * &
                                    fft_x_cnts(my_grid_idx1) * &
                                    fft_z_cnts(other_grid_idx2))
  real      :: zxy_data(2, z_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_y_cnts2(my_grid_idx2))
#else
  double precision      :: recv_buf(2 * fft_y_cnts2(my_grid_idx2) * &
                                    fft_x_cnts(my_grid_idx1) * &
                                    fft_z_cnts(other_grid_idx2))
  double precision      :: zxy_data(2, z_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_y_cnts2(my_grid_idx2))
#endif

! Local variables:

  integer               :: i, j, k
  integer               :: buf_idx
  integer               :: x_cnt
  integer               :: my_y_cnt
  integer               :: other_z_cnt
  integer               :: other_z_off

  x_cnt = fft_x_cnts(my_grid_idx1)
  my_y_cnt = fft_y_cnts2(my_grid_idx2)
  other_z_cnt = fft_z_cnts(other_grid_idx2)
  other_z_off = fft_z_offs(other_grid_idx2)

  buf_idx = 0
  do k = 1, other_z_cnt
    do j = 1, x_cnt
      do i = 1, my_y_cnt
        buf_idx = buf_idx + 1
        zxy_data(1, other_z_off + k, j, i) = recv_buf(buf_idx)
        buf_idx = buf_idx + 1
        zxy_data(2, other_z_off + k, j, i) = recv_buf(buf_idx)
      end do
    end do
  end do
 
  return

end subroutine do_yz_transpose_recv

!*******************************************************************************
!
! Internal Subroutine:  yz_transpose_my_data
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine yz_transpose_my_data(y_dim, z_dim, yxz_data, zxy_data)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: y_dim
  integer               :: z_dim
#ifdef pmemd_SPDP
  real      :: yxz_data(2, y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
  real      :: zxy_data(2, z_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_y_cnts2(my_grid_idx2))
#else
  double precision      :: yxz_data(2, y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
  double precision      :: zxy_data(2, z_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_y_cnts2(my_grid_idx2))
#endif

! Local variables:

  integer               :: i, j, k
  integer               :: x_cnt
  integer               :: my_y_off, my_y_cnt
  integer               :: my_z_off, my_z_cnt

  x_cnt = fft_x_cnts(my_grid_idx1)
  my_y_off = fft_y_offs2(my_grid_idx2)
  my_y_cnt = fft_y_cnts2(my_grid_idx2)
  my_z_off = fft_z_offs(my_grid_idx2)
  my_z_cnt = fft_z_cnts(my_grid_idx2)

  do k = 1, my_z_cnt
    do j = 1, x_cnt
      do i = 1, my_y_cnt
        zxy_data(:, my_z_off + k, j, i) = yxz_data(:, my_y_off + i, j, k)
      end do
    end do
  end do
 
  return

end subroutine yz_transpose_my_data

end subroutine do_yz_transpose

!*******************************************************************************
!
! Subroutine:  do_forward_z_fft
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine do_forward_z_fft(z_dim, z_runs)

  use parallel_dat_mod
  use processor_mod, only: fft_halo_size, myblockid
  implicit none

! Formal arguments:

  integer               :: z_dim
#ifdef pmemd_SPDP
  real                  :: z_runs(2, z_dim, &
                                  fft_x_cnts(my_grid_idx1), &
                                  fft_y_cnts2(my_grid_idx2))
#else /*pmemd_SPDP*/
  double precision      :: z_runs(2, z_dim, &
                                  fft_x_cnts(my_grid_idx1), &
                                  fft_y_cnts2(my_grid_idx2))
#endif

! Local variables:

  integer               :: j, k, i, m

  do k = 1, fft_y_cnts2(my_grid_idx2)
    do j = 1, fft_x_cnts(my_grid_idx1)

      call fft1d_forward(fft_z_hdl, z_runs(1, 1, j, k))

    end do
  end do
#if 0
print *, "new after for z",z_dim,fft_x_offs(my_grid_idx1),fft_x_cnts(my_grid_idx1), fft_y_offs2(my_grid_idx2),fft_y_cnts2(my_grid_idx2), myblockid
  do k=1, fft_y_cnts2(my_grid_idx2)
    do j=1, fft_x_cnts(my_grid_idx1)
      do i=1, z_dim
        do m=1,2
          write(15+myblockid,*) "new aft forw z",2*(i-1)+m,fft_x_offs(my_grid_idx1)+j,fft_y_offs2(my_grid_idx2)+k,z_runs(2,i,j,k) 
        end do
      end do
     end do
   end do
#endif

  return
      
end subroutine do_forward_z_fft

!*******************************************************************************
!
! Subroutine:  fft3d_back_midpoint
!
! Description:  <TBS>
!              
!*******************************************************************************

subroutine fft3d_back_midpoint(zxy_data, xyz_data, x_dim, y_dim, z_dim, &
                            new_fft_y_cnts, new_fft_y_offs)

  use parallel_dat_mod
  use processor_mod, only: fft_halo_size, myblockid, proc_max_nfft2,proc_min_nfft2, &
                                 proc_max_nfft3,proc_min_nfft3, &
                                 proc_dimx_recip,proc_dimy_recip, blockunflatten_recip 
  use timers_mod

  implicit none

! Formal arguments:

  integer                               :: x_dim, y_dim, z_dim

#ifdef pmemd_SPDP
  real, intent(in out)      :: zxy_data(2, z_dim, &
                                                    fft_x_cnts(my_grid_idx1), &
                                                    fft_y_cnts2(my_grid_idx2))

  real, intent(out)         :: xyz_data(2, -1*fft_halo_size/2:x_dim+fft_halo_size/2+1, &
                                -1*fft_halo_size: proc_max_nfft2-proc_min_nfft2+fft_halo_size+1, &
                                -1*fft_halo_size: proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)
#else
  double precision, intent(in out)      :: zxy_data(2, z_dim, &
                                                    fft_x_cnts(my_grid_idx1), &
                                                    fft_y_cnts2(my_grid_idx2))
  double precision, intent(out)         :: xyz_data(2, -1*fft_halo_size/2:x_dim+fft_halo_size/2+1, &
                                -1*fft_halo_size: proc_max_nfft2-proc_min_nfft2+fft_halo_size+1, &
                                -1*fft_halo_size: proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)
#endif
  integer                               ::new_fft_y_cnts(0:proc_dimx_recip*proc_dimy_recip-1)
  integer                               ::new_fft_y_offs(0:proc_dimx_recip*proc_dimy_recip-1)
  !integer                               :: new_fft_y_cnts(*)
  !integer                               :: new_fft_y_offs(*)

! Local variables:

#ifdef pmemd_SPDP
  real     :: yxz_data(2, y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
#else
  double precision      :: yxz_data(2, y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
#endif
  integer               :: my_x, my_y, my_z
  integer               :: my_y_offset 
  call blockunflatten_recip(mytaskid_recip,my_x,my_y,my_z) ! we can move it to setup
  my_y_offset = new_fft_y_offs(my_x+my_y*proc_dimx_recip)! we can move it to setup
  !my_y_offset helps keep the segment of data in y hanging to the offset, not
  !starting from 1. That way we can save moving to its place again in the last
  !all2all mpi call. It was used in the 1st all2all mpi call

  call do_backward_z_fft(z_dim, zxy_data)
#ifdef fft_time  
  call update_pme_time(back_z_timer)
#endif

  call do_zy_transpose(y_dim, z_dim, zxy_data, yxz_data, &
                       proc_send_buf, proc_recv_buf)
#ifdef fft_time  
  call update_pme_time(zy_transpose_timer)
#endif

  call do_backward_y_fft(y_dim, yxz_data)
#ifdef fft_time  
  call update_pme_time(back_y_timer)
#endif

  call do_yx_transpose(x_dim, y_dim, yxz_data, xyz_data, &
                       proc_send_buf, proc_recv_buf,my_y_offset)
#ifdef fft_time  
  call update_pme_time(yx_transpose_timer)
#endif

  call do_backward_x_fft(x_dim, xyz_data,my_y_offset)
#ifdef fft_time  
  call update_pme_time(back_x_timer)
#endif

  call  alltoall_last_midpoint(xyz_data, new_fft_y_cnts,new_fft_y_offs,&
                                x_dim, y_dim, z_dim,my_y_offset)
#ifdef fft_time  
  call update_pme_time(mpi_back_timer)
#endif
  return

end subroutine fft3d_back_midpoint

!*******************************************************************************
!
! Subroutine:  do_backward_z_fft
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine do_backward_z_fft(z_dim, z_runs)

  use parallel_dat_mod
  use processor_mod, only: fft_halo_size, myblockid
  implicit none

! Formal arguments:

  integer               :: z_dim
#ifdef pmemd_SPDP
  real      :: z_runs(2, z_dim, &
                                  fft_x_cnts(my_grid_idx1), &
                                  fft_y_cnts2(my_grid_idx2))
#else
  double precision      :: z_runs(2, z_dim, &
                                  fft_x_cnts(my_grid_idx1), &
                                  fft_y_cnts2(my_grid_idx2))
#endif
! Local variables:

  integer               :: j, k, i, m
  integer,save          :: counts= 0
  
  counts = counts + 1

  do k = 1, fft_y_cnts2(my_grid_idx2)
    do j = 1, fft_x_cnts(my_grid_idx1)

      call fft1d_back(fft_z_hdl, z_runs(1, 1, j, k))

    end do
  end do
#if 0
if(counts .eq. 2) then
print *, "new z",z_dim,fft_y_offs2(my_grid_idx2),fft_y_cnts2(my_grid_idx2), fft_x_offs(my_grid_idx1),fft_x_cnts(my_grid_idx1), myblockid
  do k=1, fft_x_cnts(my_grid_idx1)
    do j=1, fft_y_cnts2(my_grid_idx2)
      do i=1, z_dim
       do m = 1,2
          write(15+myblockid,*) "new back z",fft_x_offs(my_grid_idx1)+k,fft_y_offs2(my_grid_idx2)+j,2*(i-1)+m,z_runs(m,i,k,j)                                            
       end do
      end do
     end do
   end do
end if
#endif

  return
      
end subroutine do_backward_z_fft

!*******************************************************************************
!
! Subroutine:  do_zy_transpose
!
! Description: Fully asynchronous implementation!
!              
!*******************************************************************************

subroutine do_zy_transpose(y_dim, z_dim, zxy_data, yxz_data, send_buf, recv_buf)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: y_dim
  integer               :: z_dim
#ifdef pmemd_SPDP
  real      :: zxy_data(2 * z_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_y_cnts2(my_grid_idx2))
  real      :: yxz_data(2 * y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
  real      :: send_buf(siz_fft_mpi_buf2 * fft_blks_dim2)
  real      :: recv_buf(siz_fft_mpi_buf2, fft_blks_dim2)
#else
  double precision      :: zxy_data(2 * z_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_y_cnts2(my_grid_idx2))
  double precision      :: yxz_data(2 * y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
  double precision      :: send_buf(siz_fft_mpi_buf2 * fft_blks_dim2)
  double precision      :: recv_buf(siz_fft_mpi_buf2, fft_blks_dim2)
#endif

! Local variables:

  integer               :: i, j, k
  integer               :: numval, buf_base, buf_off
  integer               :: yz_idxlst_idx, recv_task, send_task
  integer               :: other_grid_idx2
  integer               :: my_y_cnt
  integer               :: my_z_cnt
  integer               :: other_y_cnt
  integer               :: other_z_cnt, other_z_off
  integer               :: x_cnt
  integer               :: irecv_stat(mpi_status_size)
  integer               :: isend_stat(mpi_status_size, fft_blks_dim2 - 1)
  integer               :: recv_req(fft_blks_dim2)
  integer               :: send_req(fft_blks_dim2)

#ifdef COMM_TIME_TEST
  call start_test_timer(11, 'do_zy_transpose', 0)
#endif

  my_y_cnt = fft_y_cnts2(my_grid_idx2)
  my_z_cnt = fft_z_cnts(my_grid_idx2)

  ! The x data is constant for the yz transpose:

  x_cnt = fft_x_cnts(my_grid_idx1)

  ! Post the asynchronous receives first:

  do yz_idxlst_idx = 1, fft_blks_dim2 - 1
    other_grid_idx2 = fft_yz_idxlst(yz_idxlst_idx)
    other_y_cnt = fft_y_cnts2(other_grid_idx2)

    recv_task = fft_task_grid(my_grid_idx1, other_grid_idx2)

    numval = 2 * my_z_cnt * x_cnt * other_y_cnt

#ifdef pmemd_SPDP
    call mpi_irecv(recv_buf(1, yz_idxlst_idx), numval, mpi_real, &
#else
    call mpi_irecv(recv_buf(1, yz_idxlst_idx), numval, mpi_double_precision, &
#endif
                   recv_task, zyyzt_tag, pmemd_comm, &
                   recv_req(yz_idxlst_idx), err_code_mpi); MPI_ERR_CHK("in fft ")
  end do


  ! Now set up and post the asynchronous sends:

  buf_base = 1

  do yz_idxlst_idx = 1, fft_blks_dim2 - 1
    other_grid_idx2 = fft_yz_idxlst(yz_idxlst_idx)
    other_z_off = fft_z_offs(other_grid_idx2) * 2
    other_z_cnt = fft_z_cnts(other_grid_idx2) * 2

    send_task = fft_task_grid(my_grid_idx1, other_grid_idx2)
    buf_off = buf_base - 1

    do k = 1, my_y_cnt
      do j = 1, x_cnt
        send_buf(buf_off+1:buf_off+other_z_cnt) = &
          zxy_data(other_z_off+1:other_z_off+other_z_cnt, j, k)
        buf_off = buf_off + other_z_cnt
      end do
    end do

    call mpi_isend(send_buf(buf_base), buf_off - buf_base + 1, &
#ifdef pmemd_SPDP
                   mpi_real, send_task, zyyzt_tag, pmemd_comm, &
#else
                   mpi_double_precision, send_task, zyyzt_tag, pmemd_comm, &
#endif
                   send_req(yz_idxlst_idx), err_code_mpi)

    buf_base = buf_base + siz_fft_mpi_buf2

  end do


  ! Do your own transposes, overlapped with i/o.

  call zy_transpose_my_data(y_dim, z_dim, zxy_data, yxz_data)

  ! Wait on and process the pending receive requests:

  do j = 1, fft_blks_dim2 - 1

    call mpi_waitany(fft_blks_dim2 - 1, recv_req, yz_idxlst_idx, irecv_stat, &
                     err_code_mpi)
    other_grid_idx2 = fft_yz_idxlst(yz_idxlst_idx)
    recv_task = fft_task_grid(my_grid_idx1, other_grid_idx2)
    call do_zy_transpose_recv(y_dim, other_grid_idx2, &
                              recv_buf(1, yz_idxlst_idx), yxz_data)

  end do

  ! Wait for all sends to complete:

  call mpi_waitall(fft_blks_dim2 - 1, send_req, isend_stat, err_code_mpi)

#ifdef COMM_TIME_TEST
  call stop_test_timer(11)
#endif

  return

contains

!*******************************************************************************
!
! Internal Subroutine:  do_zy_transpose_recv
!
! Description: Fully asynchronous implementation!
!              
!*******************************************************************************

subroutine do_zy_transpose_recv(y_dim, other_grid_idx2, recv_buf, yxz_data)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: y_dim
  integer               :: other_grid_idx2
#ifdef pmemd_SPDP
  real      :: recv_buf(2 * fft_z_cnts(my_grid_idx2) * &
                                    fft_x_cnts(my_grid_idx1) * &
                                    fft_y_cnts2(other_grid_idx2))
  real      :: yxz_data(2, y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
#else
  double precision      :: recv_buf(2 * fft_z_cnts(my_grid_idx2) * &
                                    fft_x_cnts(my_grid_idx1) * &
                                    fft_y_cnts2(other_grid_idx2))
  double precision      :: yxz_data(2, y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
#endif
! Local variables:

  integer               :: i, j, k
  integer               :: buf_idx
  integer               :: x_cnt
  integer               :: my_z_cnt
  integer               :: other_y_cnt
  integer               :: other_y_off

  x_cnt = fft_x_cnts(my_grid_idx1)
  my_z_cnt = fft_z_cnts(my_grid_idx2)
  other_y_cnt = fft_y_cnts2(other_grid_idx2)
  other_y_off = fft_y_offs2(other_grid_idx2)

  buf_idx = 0
  do k = 1, other_y_cnt
    do j = 1, x_cnt
      do i = 1, my_z_cnt
        buf_idx = buf_idx + 1
        yxz_data(1, other_y_off + k, j, i) = recv_buf(buf_idx)
        buf_idx = buf_idx + 1
        yxz_data(2, other_y_off + k, j, i) = recv_buf(buf_idx)
      end do
    end do
  end do
 
  return

end subroutine do_zy_transpose_recv

!*******************************************************************************
!
! Internal Subroutine:  zy_transpose_my_data
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine zy_transpose_my_data(y_dim, z_dim, zxy_data, yxz_data)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: y_dim
  integer               :: z_dim
#ifdef pmemd_SPDP
  real      :: zxy_data(2, z_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_y_cnts2(my_grid_idx2))
  real      :: yxz_data(2, y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
#else
  double precision      :: zxy_data(2, z_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_y_cnts2(my_grid_idx2))
  double precision      :: yxz_data(2, y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
#endif
! Local variables:

  integer               :: i, j, k
  integer               :: x_cnt
  integer               :: my_y_off, my_y_cnt
  integer               :: my_z_off, my_z_cnt

  x_cnt = fft_x_cnts(my_grid_idx1)
  my_y_off = fft_y_offs2(my_grid_idx2)
  my_y_cnt = fft_y_cnts2(my_grid_idx2)
  my_z_off = fft_z_offs(my_grid_idx2)
  my_z_cnt = fft_z_cnts(my_grid_idx2)

  do k = 1, my_z_cnt
    do j = 1, x_cnt
      do i = 1, my_y_cnt
        yxz_data(:, my_y_off + i, j, k) = zxy_data(:, my_z_off + k, j, i)
      end do
    end do
  end do
 
  return

end subroutine zy_transpose_my_data

end subroutine do_zy_transpose

!*******************************************************************************
!
! Subroutine:  do_backward_y_fft
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine do_backward_y_fft(y_dim, y_runs)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: y_dim
#ifdef pmemd_SPDP
  real      :: y_runs(2, y_dim, &
                                  fft_x_cnts(my_grid_idx1), &
                                  fft_z_cnts(my_grid_idx2))
#else
  double precision      :: y_runs(2, y_dim, &
                                  fft_x_cnts(my_grid_idx1), &
                                  fft_z_cnts(my_grid_idx2))
#endif
! Local variables:

  integer               :: j, k

  do k = 1, fft_z_cnts(my_grid_idx2)
    do j = 1, fft_x_cnts(my_grid_idx1)

      call fft1d_back(fft_y_hdl, y_runs(1, 1, j, k))
    end do
  end do

  return
      
end subroutine do_backward_y_fft

!*******************************************************************************
!
! Subroutine:  do_yx_transpose
!
! Description: Fully asynchronous implementation!
!              
!*******************************************************************************

subroutine do_yx_transpose(x_dim, y_dim, yxz_data, xyz_data, send_buf, &
recv_buf,my_y_offset)

  use parallel_dat_mod
  use processor_mod, only: fft_halo_size, myblockid, proc_max_nfft2,proc_min_nfft2, &
                                 proc_max_nfft3,proc_min_nfft3

  implicit none

! Formal arguments:

  integer               :: x_dim
  integer               :: y_dim
#ifdef pmemd_SPDP
  real      :: yxz_data(2 * y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
  real      :: xyz_data(-1*fft_halo_size:2 * x_dim+fft_halo_size+1, &
                                -1*fft_halo_size: proc_max_nfft2-proc_min_nfft2+fft_halo_size+1, &
                                -1*fft_halo_size: proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)
                                 !   fft_y_cnts1(my_grid_idx1), &
                                 !   fft_z_cnts(my_grid_idx2))
  real      :: send_buf(siz_fft_mpi_buf1 * fft_blks_dim1)
  real      :: recv_buf(siz_fft_mpi_buf1, fft_blks_dim1)
#else
  double precision      :: yxz_data(2 * y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
  double precision      :: xyz_data(-1*fft_halo_size:2 * x_dim+fft_halo_size+1, &
                                -1*fft_halo_size: proc_max_nfft2-proc_min_nfft2+fft_halo_size+1, &
                                -1*fft_halo_size: proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)
  double precision      :: send_buf(siz_fft_mpi_buf1 * fft_blks_dim1)
  double precision      :: recv_buf(siz_fft_mpi_buf1, fft_blks_dim1)
#endif
  integer               :: my_y_offset

! Local variables:

  integer               :: i, j, k
  integer               :: numval, buf_base, buf_off
  integer               :: xy_idxlst_idx, recv_task, send_task
  integer               :: other_grid_idx1
  integer               :: my_x_cnt
  integer               :: my_y_cnt
  integer               :: other_x_cnt
  integer               :: other_y_off, other_y_cnt
  integer               :: z_cnt
  integer               :: irecv_stat(mpi_status_size)
  integer               :: isend_stat(mpi_status_size, fft_blks_dim1 - 1)
  integer               :: recv_req(fft_blks_dim1)
  integer               :: send_req(fft_blks_dim1)

#ifdef COMM_TIME_TEST
  call start_test_timer(12, 'do_yx_transpose', 0)
#endif

  my_x_cnt = fft_x_cnts(my_grid_idx1)
  my_y_cnt = fft_y_cnts1(my_grid_idx1)

  ! The z data is constant for the xy transpose:

  z_cnt = fft_z_cnts(my_grid_idx2)

  ! Post the asynchronous receives first:

  do xy_idxlst_idx = 1, fft_blks_dim1 - 1
    other_grid_idx1 = fft_xy_idxlst(xy_idxlst_idx)
    other_x_cnt = fft_x_cnts(other_grid_idx1)

    recv_task = fft_task_grid(other_grid_idx1, my_grid_idx2)

    numval = 2 * my_y_cnt * other_x_cnt * z_cnt

#ifdef pmemd_SPDP
    call mpi_irecv(recv_buf(1, xy_idxlst_idx), numval, mpi_real, &
#else
    call mpi_irecv(recv_buf(1, xy_idxlst_idx), numval, mpi_double_precision, &
#endif
                   recv_task, yxxyt_tag, pmemd_comm, &
                   recv_req(xy_idxlst_idx), err_code_mpi); MPI_ERR_CHK("in fft ")
  end do

  ! Now set up and post the asynchronous sends:

  buf_base = 1

  do xy_idxlst_idx = 1, fft_blks_dim1 - 1
    other_grid_idx1 = fft_xy_idxlst(xy_idxlst_idx)
    other_y_off = fft_y_offs1(other_grid_idx1) * 2
    other_y_cnt = fft_y_cnts1(other_grid_idx1) * 2

    send_task = fft_task_grid(other_grid_idx1, my_grid_idx2)
    buf_off = buf_base - 1

    do k = 1, z_cnt
      do j = 1, my_x_cnt
        send_buf(buf_off+1:buf_off+other_y_cnt) = &
          yxz_data(other_y_off+1:other_y_off+other_y_cnt, j, k)
        buf_off = buf_off + other_y_cnt
      end do
    end do

    call mpi_isend(send_buf(buf_base), buf_off - buf_base + 1, &
#ifdef pmemd_SPDP
                   mpi_real, send_task, yxxyt_tag, pmemd_comm, &
#else
                   mpi_double_precision, send_task, yxxyt_tag, pmemd_comm, &
#endif
                   send_req(xy_idxlst_idx), err_code_mpi)

    buf_base = buf_base + siz_fft_mpi_buf1

  end do

  ! Do your own transposes, overlapped with i/o.

  call yx_transpose_my_data(x_dim, y_dim, yxz_data, xyz_data,my_y_offset)

  ! Wait on and process the pending receive requests:

  do j = 1, fft_blks_dim1 - 1

    call mpi_waitany(fft_blks_dim1 - 1, recv_req, xy_idxlst_idx, irecv_stat, &
                     err_code_mpi)
    other_grid_idx1 = fft_xy_idxlst(xy_idxlst_idx)
    recv_task = fft_task_grid(other_grid_idx1, my_grid_idx2)
    call do_yx_transpose_recv(x_dim, other_grid_idx1, &
                              recv_buf(1, xy_idxlst_idx), xyz_data,my_y_offset)

  end do

  ! Wait for all sends to complete:

  call mpi_waitall(fft_blks_dim1 - 1, send_req, isend_stat, err_code_mpi)

#ifdef COMM_TIME_TEST
  call stop_test_timer(12)
#endif

  return

contains

!*******************************************************************************
!
! Internal Subroutine:  do_yx_transpose_recv
!
! Description: Fully asynchronous implementation!
!              
!*******************************************************************************

subroutine do_yx_transpose_recv(x_dim, other_grid_idx1, recv_buf, &
xyz_data,my_y_offset)

  use parallel_dat_mod
  use processor_mod, only: fft_halo_size, myblockid, proc_max_nfft2,proc_min_nfft2, &
                                 proc_max_nfft3,proc_min_nfft3

  implicit none

! Formal arguments:

  integer               :: x_dim
  integer               :: other_grid_idx1
#ifdef pmemd_SPDP
  real      :: recv_buf(fft_y_cnts1(my_grid_idx1) * 2 *&
                                    fft_x_cnts(other_grid_idx1) * &
                                    fft_z_cnts(my_grid_idx2))
  real      :: xyz_data(2, -1*fft_halo_size/2:x_dim+fft_halo_size/2+1, &
                                -1*fft_halo_size:proc_max_nfft2-proc_min_nfft2+fft_halo_size+1, &
                                -1*fft_halo_size:proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)
#else
  double precision      :: recv_buf(fft_y_cnts1(my_grid_idx1) * 2 *&
                                    fft_x_cnts(other_grid_idx1) * &
                                    fft_z_cnts(my_grid_idx2))
  double precision      :: xyz_data(2, -1*fft_halo_size/2:x_dim+fft_halo_size/2+1, &
                                -1*fft_halo_size:proc_max_nfft2-proc_min_nfft2+fft_halo_size+1, &
                                -1*fft_halo_size:proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)
#endif
  integer               :: my_y_offset

! Local variables:

  integer               :: i, j, k
  integer               :: buf_idx
  integer               :: my_y_cnt
  integer               :: other_x_off, other_x_cnt
  integer               :: z_cnt

  my_y_cnt = fft_y_cnts1(my_grid_idx1)
  other_x_cnt = fft_x_cnts(other_grid_idx1)
  other_x_off = fft_x_offs(other_grid_idx1)
  z_cnt = fft_z_cnts(my_grid_idx2)

  buf_idx = 0
  do k = 1, z_cnt
    do j = 1, other_x_cnt
      do i = my_y_offset + 1, my_y_offset + my_y_cnt
        buf_idx = buf_idx + 1
        xyz_data(1, other_x_off + j, i, k) = recv_buf(buf_idx)
        buf_idx = buf_idx + 1
        xyz_data(2, other_x_off + j, i, k) = recv_buf(buf_idx)
      end do
    end do
  end do
 
  return

end subroutine do_yx_transpose_recv

!*******************************************************************************
!
! Internal Subroutine:  yx_transpose_my_data
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine yx_transpose_my_data(x_dim, y_dim, yxz_data, xyz_data,my_y_offset)

  use parallel_dat_mod
  use processor_mod, only: fft_halo_size, myblockid, proc_max_nfft2,proc_min_nfft2, &
                                 proc_max_nfft3,proc_min_nfft3

  implicit none

! Formal arguments:

  integer               :: x_dim
  integer               :: y_dim
#ifdef pmemd_SPDP
  real      :: yxz_data(2, y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
  real      :: xyz_data(2, -1*fft_halo_size/2:x_dim+fft_halo_size/2+1, &
                                -1*fft_halo_size:proc_max_nfft2-proc_min_nfft2+fft_halo_size+1, &
                                -1*fft_halo_size:proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)
#else
  double precision      :: yxz_data(2, y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
  double precision      :: xyz_data(2, -1*fft_halo_size/2:x_dim+fft_halo_size/2+1, &
                                -1*fft_halo_size:proc_max_nfft2-proc_min_nfft2+fft_halo_size+1, &
                                -1*fft_halo_size:proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)
#endif
  integer               :: my_y_offset

! Local variables:

  integer               :: i, j, k
  integer               :: my_x_off, my_x_cnt
  integer               :: my_y_off, my_y_cnt
  integer               :: z_cnt

  my_x_off = fft_x_offs(my_grid_idx1)
  my_x_cnt = fft_x_cnts(my_grid_idx1)
  my_y_off = fft_y_offs1(my_grid_idx1)
  my_y_cnt = fft_y_cnts1(my_grid_idx1)
  z_cnt = fft_z_cnts(my_grid_idx2)

  do k = 1, z_cnt
    do j = 1, my_y_cnt
      do i = 1, my_x_cnt
        xyz_data(:, my_x_off + i, my_y_offset+j, k) = yxz_data(:, my_y_off + j, i, k)
      end do
    end do
  end do
 
  return

end subroutine yx_transpose_my_data

end subroutine do_yx_transpose

!*******************************************************************************
!
! Subroutine:  do_backward_x_fft
!
! Description:
!
!*******************************************************************************

subroutine do_backward_x_fft(x_dim, x_runs,my_y_offset)

  use processor_mod, only: fft_halo_size, myblockid, proc_max_nfft2,proc_min_nfft2, &
                                 proc_max_nfft3,proc_min_nfft3
  implicit none

! Formal arguments:

  integer               :: x_dim
  integer               :: my_y_offset
#ifdef pmemd_SPDP
  real      :: x_runs(2, -1*fft_halo_size/2:x_dim+fft_halo_size/2+1, &
                                -1*fft_halo_size: proc_max_nfft2-proc_min_nfft2+fft_halo_size+1, &
                               -1*fft_halo_size:proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)
! Local variables:

  real      :: a, b, c, d
  real      :: buf(2, x_dim - 1)
#else
  double precision      :: x_runs(2, -1*fft_halo_size/2:x_dim+fft_halo_size/2+1, &
                                -1*fft_halo_size: proc_max_nfft2-proc_min_nfft2+fft_halo_size+1, &
                               -1*fft_halo_size:proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)
! Local variables:

  double precision      :: a, b, c, d
  double precision      :: buf(2, x_dim - 1)
#endif
  integer               :: i, j, k, m
  integer               :: x_lim                ! indexing limit
  integer,save          :: counts=0

  counts = counts + 1
! complex to real:

  x_lim = x_dim + 1

  do k = 1, fft_z_cnts(my_grid_idx2)

    do j = my_y_offset + 1, my_y_offset + fft_y_cnts1(my_grid_idx1)

      ! We have to put results in a temporary run to prevent overwriting of
      ! values needed later.

      do i = 2, x_dim - 1
        a = (x_runs(1, i, j, k) + x_runs(1, x_lim - i, j, k)) ! Real F even
        b = (x_runs(2, i, j, k) - x_runs(2, x_lim - i, j, k)) ! Imag F even
        c = (x_runs(2, i, j, k) + x_runs(2, x_lim - i, j, k)) ! F odd contrib
        d = (x_runs(1, i, j, k) - x_runs(1, x_lim - i, j, k)) ! F odd contrib
        buf(1, i) = a - fftrc_coefs(1, i) * c - fftrc_coefs(2, i) * d
        buf(2, i) = b + fftrc_coefs(1, i) * d - fftrc_coefs(2, i) * c
      end do

      buf(1, 1) = x_runs(1, 1, j, k) + x_runs(1, x_dim, j, k)
      buf(2, 1) = x_runs(1, 1, j, k) - x_runs(1, x_dim, j, k)

      call fft1d_back(fft_x_hdl, buf(1, 1))

      do i = 1, x_dim - 1
        x_runs(:, i, j, k) = buf(:, i)
      end do

    end do

  end do
#if 0
  do k=1, fft_z_cnts(my_grid_idx2)
    do j=1, fft_y_cnts1(my_grid_idx1)
      do i=1, x_dim
       do m = 1,2
          write(15+myblockid,*) "new back x",2*(i-1)+m,fft_y_offs1(my_grid_idx1)+j,fft_z_offs(my_grid_idx2)+k ,x_runs(m,i,j,k)                                            
       end do
      end do
     end do
   end do
#endif

  return

end subroutine do_backward_x_fft

!*******************************************************************************
!
! Subroutine: alltoall_last_midpoint 
!
! Description:  MPI communication in y direction using proc_dimx_recip division
!               to convert pencil fomat xyz_data() to cubic format xyz_data()
!               the pencil format dim was: full x_dim * partial y * partial z
!               The 3d MPI divison, it was 1 * (proc_dimx_recip*proc_dimy_recipi) * proc_dimz_recip  
!               The cubic output of xyz_data is xyz is part x * part y * part z 
!               The 3d MPI divison, it was proc_dimx_recip * proc_dimy_recipi * proc_dimz_recip  
!
! INPUT: xyz_data (in pencil fromat with full x dim)
! 
! OUTPUT: xyz_data (in cubic format with no full dim)
!
!*******************************************************************************
subroutine alltoall_last_midpoint(xyz_data, new_fft_y_cnts,new_fft_y_offs,&
                                x_dim, y_dim, z_dim,my_y_offset)

  use processor_mod, only: fft_halo_size, send_xfft_cnts, send_yfft_cnts, & 
                          send_zfft_cnts, my_xfft_cnt, my_yfft_cnt, &
                          my_zfft_cnt, x_offset, y_offset, z_offset, &
                          proc_dimx_recip, proc_dimy_recip, proc_dimz_recip, &
                          proc_x_comm_lst, &
                          blockunflatten_recip, &
                          proc_dimx_recip, proc_dimy_recip, proc_dimz_recip, &
                          proc_min_nfft1, proc_min_nfft2, proc_min_nfft3, &
                          proc_max_nfft1, proc_max_nfft2, proc_max_nfft3, &
                          myblockid
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use pbc_mod
  use prmtop_dat_mod
  use timers_mod

  implicit none
  !Formal arguments
  integer               :: x_dim, y_dim, z_dim
#ifdef pmemd_SPDP
  real      :: xyz_data(-1*fft_halo_size:2 * x_dim+fft_halo_size+1, &
                                -1*fft_halo_size:proc_max_nfft2-proc_min_nfft2+fft_halo_size+1, &
                                -1*fft_halo_size:proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)
#else
  double precision      :: xyz_data(-1*fft_halo_size:2 * x_dim+fft_halo_size+1, &
                                -1*fft_halo_size:proc_max_nfft2-proc_min_nfft2+fft_halo_size+1, &
                                -1*fft_halo_size:proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)
#endif
  integer               :: new_fft_y_cnts(0:proc_dimx_recip*proc_dimy_recip-1)
  integer               :: new_fft_y_offs(0:proc_dimx_recip*proc_dimy_recip-1)
  integer               :: my_y_offset

  !local arguments
  integer               :: x_idx
  integer               :: other_offset, other_mpi_x, other_cnt_x,&
                          other_cnt_y, other_cnt_z, other_cnt
  integer               :: my_cnt_x, my_cnt_y, my_cnt_z
  integer               :: x_tag=13
  integer               :: x,y,z !3d co-ordinate of MPI ranks
  integer               :: my_x,my_y,my_z,my_xy_index !3d co-ordinate of my MPI ranks
  integer               :: other_y_offset,my_x_offset, other_x_offset
  
  !this is temporary, this one should come from pme_recip_midpoint setup routine
  integer               :: send_cnt, recv_cnt
  integer               :: i,j,k,buf_base, buf_off 
  integer               :: other_x,other_y, other_z 
  !MPI data structure
  integer               :: irecv_stat(mpi_status_size)
  integer               :: isend_stat(mpi_status_size, proc_dimx_recip - 1)
  integer               :: recv_req(proc_dimx_recip - 1)
  integer               :: send_req(proc_dimx_recip - 1)
  integer               :: cnt_idx
  !alltoall in x direction

  !post the non-blocking recv
  !I will receive the chunk of my_x_cnt * other_y_cnt * my_z_cnt
  !at least two dim of fft chunk(here my_x_cnt and  my_z_cnt) that I receive should be same as mine, 
  !since we are communicating in y direction, y_dim can be different
  !But the communicating MPIs are originally in x direction

  call blockunflatten_recip(mytaskid_recip,my_x,my_y,my_z)
  my_cnt_x = send_xfft_cnts(my_x) 
  !my_cnt_y = new_fft_y_cnts(my_x + proc_dimx_recip*my_y) ! this is how the array was
                                                   ! created, I will also
                                                   ! receieve this amount 
  my_cnt_z = my_zfft_cnt ! = fft_z_cnts(my_grid_idx2)      
  other_offset = 0
  cnt_idx = 1
  do x_idx = 0, proc_dimx_recip - 1 ! only send to MPIs in x direction   

    !post the asynchronous receives first
    other_mpi_x  = proc_x_comm_lst(x_idx)  ! real MPI rank now in y direction
                                          !but orginally was x direction
    if(other_mpi_x .eq. mytaskid_recip) cycle !avoid receiving from me
    call blockunflatten_recip(other_mpi_x,other_x,other_y,other_z)
    other_cnt_y = new_fft_y_cnts(other_x + proc_dimx_recip*other_y) ! = fft_y_cnts(my_grid_idx1)      
    other_cnt   = my_cnt_x * other_cnt_y * my_cnt_z ! we dont need to use *2 here 

#ifdef pmemd_SPDP
    call mpi_irecv(proc_recv_buf(other_offset+1), other_cnt, mpi_real, &
#else
    call mpi_irecv(proc_recv_buf(other_offset+1), other_cnt, mpi_double_precision, &
#endif
                   other_mpi_x, x_tag, pmemd_comm, &
                   recv_req(cnt_idx), err_code_mpi); MPI_ERR_CHK("in fft ")
    other_offset = other_offset + other_cnt 
    cnt_idx = cnt_idx+1
    !print *, "Receiving cnt", other_cnt, "other_cnt_y", other_cnt_y, ",other_mpi_x ",other_mpi_x, mytaskid_recip
  end do

  !pack data and  post the non-blocking send
  !I will send the chunk of other_cnt_x * my_cnt_y * (my=other)cnt_z
  my_cnt_y = new_fft_y_cnts(my_x + proc_dimx_recip*my_y) ! = fft_y_cnts(my_grid_idx1)      
  my_xy_index = my_y * proc_dimx_recip !this is the index offset from where my_yfft_cnt
                                 !starts in full arrays new_fft_y_cnts() and new_fft_y_offs()
  buf_base = 1
  cnt_idx  = 1
  do x_idx = 0, proc_dimx_recip - 1 ! make sure it is right.. only send to MPIs in x direction, to get a full x_dim   
  
    !pack the data
    other_mpi_x  = proc_x_comm_lst(x_idx)  ! convert to real MPI rank
    if(other_mpi_x .eq. mytaskid_recip) cycle !avoid sending to me
    send_cnt = 0 ! proc_min_nfft1 !offset to start my block in x
    other_cnt_x =    send_xfft_cnts(x_idx)
    other_x_offset = x_offset(x_idx) 
    buf_off = buf_base - 1
    do k  = 1, my_cnt_z !my_zfft_cnt
      do j = my_y_offset + 1, my_y_offset + my_cnt_y !other_y_offset + 1, other_y_offset + other_cnt_y !my_yfft_cnt/proc_dimx_recip
                          !the offset is local only to divide my_yfft_cnt, and get the proper chunk 
        proc_send_buf(buf_off+1:buf_off+other_cnt_x) = xyz_data(other_x_offset+1:other_x_offset+other_cnt_x,j,k) 
                                  !chunk is properly spaced with other_x_offset
        buf_off = buf_off + other_cnt_x !send_cnt is offset for send_bufx
      end do
    end do

    !post the non-blocking send
    call mpi_isend(proc_send_buf(buf_base), buf_off - buf_base + 1, &
#ifdef pmemd_SPDP
                   mpi_real, other_mpi_x, x_tag, pmemd_comm, &
#else
                   mpi_double_precision, other_mpi_x, x_tag, pmemd_comm, &
#endif
                   send_req(cnt_idx), err_code_mpi)
    buf_base = buf_base + other_cnt_x * my_cnt_y * my_cnt_z ! these many are sent in one MPI_isend 
    cnt_idx = cnt_idx + 1
!    print *, "sending _cnt", other_cnt_x * my_cnt_y * my_cnt_z, "my_cnt_y", my_cnt_y, ",other_mpi_x ",other_mpi_x, mytaskid_recip
  end do
  
  ! Wait for all the MPI receive to complete
  do x_idx = 1, proc_dimx_recip - 1 ! receiving from proc_dimx_recip-1 ranks
    call mpi_waitany(proc_dimx_recip - 1, recv_req, cnt_idx, irecv_stat, &
                     err_code_mpi) !instead of cnt_idx any other variable can be used, but dont use x_idx, that will change the index of the loop,a big bug showed up when x_idx was used here in palce of cnt_idx 
  end do

  ! Wait for all MPI sends to complete:
  call mpi_waitall(proc_dimx_recip - 1, send_req, isend_stat, err_code_mpi)
  !unpack the recv data
  recv_cnt = 0
 ! my_y_offset = new_fft_y_offs(my_x+my_xy_index)
                                !To place y into proper place in xyz_data
                                !since we divided my_yfft-cnt into more
                                !proc_dimx_recip divisions
                                !we could build new_fft_y_offset() as
                                !0:proc_dimx_recip-1 long, then my_xy_index would
                                !not be required as index for new_fft_y_offset
  !if my y does not start from 1, but current data starts at 1,
  !so move it to my_y_offset position, so that incoming data is moved to
  !their original position(origina = before the fft started) 
#if 0
  if(my_y_offset .gt. 0) then ! it means it starts not from 1, so copying and pasting down at 1 
   !if(myblockid .eq. 1) print *, "for 1",proc_min_nfft1,my_cnt_x,my_y_offset,my_cnt_y, my_cnt_z
    do k = 1, my_cnt_z
      do j =  my_cnt_y,1,-1 !not 1,my_cnt_y, to resovle a bug, otherwise it will override its own data 
        do i = proc_min_nfft1, proc_min_nfft1+my_cnt_x-1
          xyz_data(i,my_y_offset+j,k) = xyz_data(i,j,k) 
        end do
      end do
    end do

  end if
#endif
  !unpack it and position it proper x and y offset 
  my_x_offset = x_offset(my_x)
  do x_idx = 0, proc_dimx_recip - 1
    other_mpi_x = proc_x_comm_lst(x_idx) 
    if(other_mpi_x .eq. mytaskid_recip) cycle !avoid receiving from  me
    call blockunflatten_recip(other_mpi_x,other_x,other_y,other_z)!getting other x,y,z
    other_y_offset = new_fft_y_offs(other_x + proc_dimx_recip*other_y) 
    other_cnt_y = new_fft_y_cnts(other_x + proc_dimx_recip*other_y) 
    do k = 1, my_cnt_z
      do j = other_y_offset+1, other_y_offset+other_cnt_y !position in proper place
        do i=1, my_cnt_x
          xyz_data(my_x_offset+i,j,k) = proc_recv_buf(recv_cnt+i)
        end do
        recv_cnt = recv_cnt + my_cnt_x
      end do
    end do
  end do
  
end subroutine alltoall_last_midpoint

#endif /*MPI*/

end module pme_fft_midpoint_mod
