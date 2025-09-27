#include "copyright.i"

!*******************************************************************************
!
! Module: pme_blk_fft_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module pme_blk_fft_mod

  use fft1d_mod
  use pme_recip_dat_mod
  use pme_fft_dat_mod
#ifdef TIME_TEST
  use timers_mod
#endif

!precision setting for SPDP, DPDP, SPSP
#include "include_precision.i"
  implicit none

! None of this stuff is broadcast; fft setup occurs in all tasks; note that
! the stuff here is only set up and used if block (blk), as opposed to slab 
! fft's are in effect (only used for mpi parallel code). A "block" is
! essentially a collection of 1D runs over which 1d fft's are to be performed.

#ifdef MPI
! Data describing fft blocks and ownership:

! These are arrays describing how the fft grid is divided between tasks;
! these will all be of dimension fft_blks_dim.

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

  double precision, parameter           :: blk_fft_workload_estimate = 0.15d0

contains

!*******************************************************************************
!
! Subroutine:  blk_fft_setup
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine blk_fft_setup(min_numtasks, recip_numtasks, &
                         frc_ene_numtasks, frc_ene_task_lst, &
                         step_ctr, num_ints, num_reals)

  use file_io_mod
  use gbl_constants_mod
  use mdin_ctrl_dat_mod, only : loadbal_verbose
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:

  ! recip_numtasks returns the number of tasks doing recip force calcs.

  integer, intent(in)           :: min_numtasks   ! suggested, not firm.
  integer, intent(out)          :: recip_numtasks ! output is new value.
  integer, intent(in)           :: frc_ene_numtasks
  integer, intent(in)           :: frc_ene_task_lst(frc_ene_numtasks)
  integer, intent(in)           :: step_ctr

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed
  integer               :: i, j, k
  integer               :: numtasks_candidate
  integer               :: new_fft_blks_dim2
  integer               :: quot, rem            ! quotient, remainder for
                                                ! grid division
  integer               :: cnts_idx, offs_idx
  integer               :: lst_idx
  integer               :: offset
  integer               :: taskid
  integer               :: taskmap_idx
  integer, save         :: distrib_cnt = 0

  distrib_cnt = distrib_cnt + 1 ! count of calls to this routine

  if (numtasks .lt. 4) then
    if (master) then
      write(mdout, '(a,a)') error_hdr, &
         'numtasks must be >= 4 in parallel_blk_setup()!'
    end if
    call mexit(mdout, 1)
  end if

  if (frc_ene_numtasks .lt. 4) then
    if (master) then
      write(mdout, '(a,a)') error_hdr, &
         'frc_ene_numtasks must be >= 4 in parallel_blk_setup()!'
    end if
    call mexit(mdout, 1)
  end if

  if (frc_ene_numtasks .lt. min_numtasks) then
    if (master) then
      write(mdout, '(a,a)') error_hdr, &
         'frc_ene_numtasks must be >= min_numtasks in parallel_blk_setup()!'
    end if
    call mexit(mdout, 1)
  end if

  max_fft_blks_dim1 = min(fft_x_dim, fft_y_dim)
  max_fft_blks_dim2 = min(fft_y_dim, fft_z_dim)

  fft_blks_dim1 = fft_blk_y_divisor

  ! Note that this algorithm is dependent on frc_ene_numtasks .gt. 4 and
  ! an initial value for fft_blks_dim2 of 2!

  new_fft_blks_dim2 = fft_blks_dim2

  do i = fft_blks_dim2, max_fft_blks_dim2
    numtasks_candidate = fft_blks_dim1 * i
    if (numtasks_candidate .le. frc_ene_numtasks) then
      new_fft_blks_dim2 = i
      recip_numtasks = numtasks_candidate
      if (numtasks_candidate .ge. min_numtasks) exit
    else
      exit
    end if
  end do

  fft_blks_dim2 = new_fft_blks_dim2

  ! It is possible that the size of the data structures has not changed, but
  ! it is probably better to avoid unnecessary calls to this routine instead
  ! of avoiding deallocation/identical reallocation.  We use the allocation
  ! state of fft_task_grid to indicate the allocation state of all dynamically
  ! allocated data in this module.

  if (allocated(fft_task_grid)) then

    num_ints = num_ints - size(fft_x_offs) - &
                          size(fft_x_cnts) - &
                          size(fft_y_offs1) - &
                          size(fft_y_cnts1) - &
                          size(fft_y_offs2) - &
                          size(fft_y_cnts2) - &
                          size(fft_z_offs) - &
                          size(fft_z_cnts) - &
                          size(fft_task_grid) - &
                          size(fft_xy_idxlst) - &
                          size(fft_yz_idxlst)

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
           fft_xy_idxlst(fft_blks_dim1 - 1), &
           fft_yz_idxlst(fft_blks_dim2 - 1), &
           stat = alloc_failed)
               
  if (alloc_failed .ne. 0) call setup_alloc_error

  num_ints = num_ints + size(fft_x_offs) + &
                        size(fft_x_cnts) + &
                        size(fft_y_offs1) + &
                        size(fft_y_cnts1) + &
                        size(fft_y_offs2) + &
                        size(fft_y_cnts2) + &
                        size(fft_z_offs) + &
                        size(fft_z_cnts) + &
                        size(fft_task_grid) + &
                        size(fft_xy_idxlst) + &
                        size(fft_yz_idxlst)

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

  ! Now setup the offsets and counts tables.

  ! -- first x...

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

  quot = fft_y_dim / fft_blks_dim1
  rem = fft_y_dim - quot * fft_blks_dim1

  fft_y_cnts1(:) = quot

  if (rem .ne. 0) then
    do i = 1, rem
      cnts_idx = fft_blks_dim1 * i / rem
      fft_y_cnts1(cnts_idx) = quot + 1
    end do
    max_fft_y_cnts1 = quot + 1
  else
    max_fft_y_cnts1 = quot
  end if

  offset = 0
  do i = 1, fft_blks_dim1
    fft_y_offs1(i) = offset
    offset = offset + fft_y_cnts1(i)
  end do

  ! -- then y in dim 2...

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

  quot = fft_z_dim / fft_blks_dim2
  rem = fft_z_dim - quot * fft_blks_dim2

  fft_z_cnts(:) = quot

  if (rem .ne. 0) then
    do i = 1, rem
      cnts_idx = fft_blks_dim2 * i / rem
      fft_z_cnts(cnts_idx) = quot + 1
    end do
    max_fft_z_cnts = quot + 1
  else
    max_fft_z_cnts = quot
  end if

  offset = 0
  do i = 1, fft_blks_dim2
    fft_z_offs(i) = offset
    offset = offset + fft_z_cnts(i)
  end do

  ! Now set up the tasks table. We choose tasks equally spaced through the
  ! id space, thereby hopefully not overloading one smp...  We also do this
  ! from top to bottom, which prevents use of the master task.  We actually
  ! have to use indirect references through the frc_ene_task_lst() because
  ! the task id space is not dense if there are any comm tasks or the master
  ! has been idled (not implemented yet).
  
  is_recip_task(:) = .false.
  i_do_recip = .false.

  my_grid_idx1 = 0
  my_grid_idx2 = 0
  i = fft_blks_dim1
  j = fft_blks_dim2
  do k = 0, recip_numtasks - 1
    lst_idx = frc_ene_numtasks - &
              int(dble(k) * dble(frc_ene_numtasks) / dble(recip_numtasks))
    taskid = frc_ene_task_lst(lst_idx)
    is_recip_task(taskid) = .true.
    fft_task_grid(i, j) = taskid
    if (taskid .eq. mytaskid) then
      my_grid_idx1 = i      
      my_grid_idx2 = j      
      i_do_recip = .true.
    end if
    i = i - 1
    if (i .eq. 0) then
      i = fft_blks_dim1
      j = j - 1
    end if
  end do

  ! Now set up the fft xy index list and mpi fft buffers in reciprocal tasks.
  ! The index list is used to control order of access to taskid's and data
  ! offsets and counts in the xy transpose.

  if (i_do_recip) then

    taskmap_idx = my_grid_idx1 + 1
    do i = 1, fft_blks_dim1 - 1
      if (taskmap_idx .gt. fft_blks_dim1) taskmap_idx = 1
      fft_xy_idxlst(i) = taskmap_idx
      taskmap_idx = taskmap_idx + 1
    end do

    taskmap_idx = my_grid_idx2 + 1
    do i = 1, fft_blks_dim2 - 1
      if (taskmap_idx .gt. fft_blks_dim2) taskmap_idx = 1
      fft_yz_idxlst(i) = taskmap_idx
      taskmap_idx = taskmap_idx + 1
    end do

    ! Allocate space for multidimensional fft transposes.
  
    siz_fft_mpi_buf1 = 2 * max_fft_x_cnts * max_fft_y_cnts1 * max_fft_z_cnts
    siz_fft_mpi_buf2 = 2 * max_fft_x_cnts * max_fft_y_cnts2 * max_fft_z_cnts

    call set_minimum_mpi_bufs_size(siz_fft_mpi_buf1 * (fft_blks_dim1 - 1), &
                                   num_reals)

    call set_minimum_mpi_bufs_size(siz_fft_mpi_buf2 * (fft_blks_dim2 - 1), &
                                   num_reals)
  end if

  if (master) then

! BEGIN DBG
! write(0, *)'fft_x,y,z_dim=', fft_x_dim, fft_y_dim, fft_z_dim
! write(0, *)'siz_fft_mpi_buf1=', &
!   2 * max_fft_x_cnts * max_fft_y_cnts1 * max_fft_z_cnts
! write(0, *)'siz_fft_mpi_buf2=', &
!   2 * max_fft_x_cnts * max_fft_y_cnts2 * max_fft_z_cnts
! write(0, *)'is_recip_task()=', is_recip_task(:)

! write(0, *)'fft_x_offs()=', fft_x_offs(:)
! write(0, *)'fft_x_cnts()=', fft_x_cnts(:)
! write(0, *)'fft_y_offs1()=', fft_y_offs1(:)
! write(0, *)'fft_y_cnts1()=', fft_y_cnts1(:)
! write(0, *)'fft_y_offs2()=', fft_y_offs2(:)
! write(0, *)'fft_y_cnts2()=', fft_y_cnts2(:)
! write(0, *)'fft_z_offs()=', fft_z_offs(:)
! write(0, *)'fft_z_cnts()=', fft_z_cnts(:)

! write(0, *)'max_fft_x_cnts=', max_fft_x_cnts
! write(0, *)'max_fft_y_cnts1=', max_fft_y_cnts1
! write(0, *)'max_fft_y_cnts2=', max_fft_y_cnts2
! write(0, *)'max_fft_z_cnts=', max_fft_z_cnts

! write(0, *)'fft_blks_dim1=', fft_blks_dim1
! write(0, *)'fft_blks_dim2=', fft_blks_dim2
! write(0, *)'max_fft_blks_dim1=', max_fft_blks_dim1
! write(0, *)'max_fft_blks_dim2=', max_fft_blks_dim2

! write(0, *)'fft_task_grid(:,:)=', fft_task_grid(:,:)
! END DBG

10 format(/, a, /)
12 format(/, a, i5, a, i7, a /)
20 format('  ', a, i4, a)
30 format('  ', a)
40 format('    ', 16i4)

    if (loadbal_verbose .le. 0) then

      if (distrib_cnt .le. 2) then

        ! Initial distribution info for run where dynamic slab distribution
        ! is turned on; data is for first two dynamic distributions based on 
        ! the initial workload guess followed by the first actual workload
        ! data.

        if (distrib_cnt .eq. 1) then
          write(logfile, 10) &
            'Initial FFT Block Distribution Based on Workload Estimate:'
        else
          write(logfile, 10) &
            'First FFT Block Distribution Based on Actual Workload:'
        end if

        write(logfile, 20) 'FFT blocks assigned to ', recip_numtasks, ' tasks'

      end if

    else

      ! Verbose mode.  Report fft slab redistribution for every call...

        write(logfile, 12) &
          'FFT Block Distribution No. ', distrib_cnt, ' at run step ', &
          step_ctr, ':'
        write(logfile, 20) 'FFT blocks assigned to ', recip_numtasks, ' tasks'

    end if

  end if

! BEGIN DBG
! if (i_do_recip) then
!   write(0, *)'task=', mytaskid, ' fft_xy_idxlst(:)=', fft_xy_idxlst(:)
!   write(0, *)'task=', mytaskid, ' fft_yz_idxlst(:)=', fft_yz_idxlst(:)
! end if
! END DBG

  return

end subroutine blk_fft_setup

!*******************************************************************************
!
! Subroutine:  blk_fft3drc_forward
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine blk_fft3drc_forward(xyz_data, zxy_data, x_dim, y_dim, z_dim)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer                               :: x_dim, y_dim, z_dim

  double precision, intent(in out)      :: xyz_data(2, x_dim, &
                                                    fft_y_cnts1(my_grid_idx1), &
                                                    fft_z_cnts(my_grid_idx2))
  
  double precision, intent(out)         :: zxy_data(2, z_dim, &
                                                    fft_x_cnts(my_grid_idx1), &
                                                    fft_y_cnts2(my_grid_idx2))
! Local variables:

  double precision      :: yxz_data(2, y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))

  ! Do all the forward x fft's in your block, leaving results in xyz_data;
  ! this data will not be used as such by the caller, but is left here to
  ! reduce the total data footprint and cache thrashing.
  
  call do_forward_x_fft(x_dim, xyz_data)

  call do_xy_transpose(x_dim, y_dim, xyz_data, yxz_data, &
                       dbl_mpi_send_buf, dbl_mpi_recv_buf)

  call do_forward_y_fft(y_dim, yxz_data)

  call do_yz_transpose(y_dim, z_dim, yxz_data, zxy_data, &
                       dbl_mpi_send_buf, dbl_mpi_recv_buf)

  call do_forward_z_fft(z_dim, zxy_data)

  return

end subroutine blk_fft3drc_forward

!*******************************************************************************
!
! Subroutine:  do_forward_x_fft
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine do_forward_x_fft(x_dim, x_runs)

  implicit none

! Formal arguments:

  integer               :: x_dim
  double precision      :: x_runs(2, x_dim, &
                                  fft_y_cnts1(my_grid_idx1), &
                                  fft_z_cnts(my_grid_idx2))

! Local variables:

  double precision      :: a, b, c, d
  double precision      :: buf(2, x_dim - 1)
  integer               :: i, j, k
  integer               :: x_lim        ! indexing limit

! real to complex:

  x_lim = x_dim + 1

  do k = 1, fft_z_cnts(my_grid_idx2)

    do j = 1, fft_y_cnts1(my_grid_idx1)

      buf(:,:) = x_runs(:, 1:x_dim - 1, j, k)

#ifndef MIDPOINT_SPDP /*this will not be used if MIDPOINT_SPDP is defined*/
      call fft1d_forward(fft_x_hdl, buf)
#endif

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

  return
      
end subroutine do_forward_x_fft

!*******************************************************************************
!
! Subroutine:  do_xy_transpose
!
! Description: Fully asynchronous implementation!
!              
!*******************************************************************************

subroutine do_xy_transpose(x_dim, y_dim, xyz_data, yxz_data, send_buf, recv_buf)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: x_dim
  integer               :: y_dim
  double precision      :: xyz_data(2 * x_dim, &
                                    fft_y_cnts1(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
  double precision      :: yxz_data(2 * y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
  double precision      :: send_buf(siz_fft_mpi_buf1 * (fft_blks_dim1 - 1))
  double precision      :: recv_buf(siz_fft_mpi_buf1, (fft_blks_dim1 - 1))

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
  integer               :: recv_req(fft_blks_dim1 - 1)
  integer               :: send_req(fft_blks_dim1 - 1)

#ifdef COMM_TIME_TEST
  call start_test_timer(9, 'do_xy_transpose', 0)
#endif

  my_x_cnt = fft_x_cnts(my_grid_idx1)
  my_y_cnt = fft_y_cnts1(my_grid_idx1)

  ! The z data is constant for the xy transpose:

  z_cnt = fft_z_cnts(my_grid_idx2)

  ! Post the asynchronous receives first:

  do xy_idxlst_idx = 1, fft_blks_dim1 - 1
    other_grid_idx1 = fft_xy_idxlst(xy_idxlst_idx)
    other_y_cnt = fft_y_cnts1(other_grid_idx1)

    recv_task = fft_task_grid(other_grid_idx1, my_grid_idx2)

    numval = 2 * my_x_cnt * other_y_cnt * z_cnt

    call mpi_irecv(recv_buf(1, xy_idxlst_idx), numval, mpi_double_precision, &
                   recv_task, xyyxt_tag, pmemd_comm, &
                   recv_req(xy_idxlst_idx), err_code_mpi) 
  end do

  ! Now set up and post the asynchronous sends:

  buf_base = 1

  do xy_idxlst_idx = 1, fft_blks_dim1 - 1
    other_grid_idx1 = fft_xy_idxlst(xy_idxlst_idx)
    other_x_off = fft_x_offs(other_grid_idx1) * 2
    other_x_cnt = fft_x_cnts(other_grid_idx1) * 2

    send_task = fft_task_grid(other_grid_idx1, my_grid_idx2)
    buf_off = buf_base - 1

    do k = 1, z_cnt
      do j = 1, my_y_cnt
        send_buf(buf_off+1:buf_off+other_x_cnt) = &
          xyz_data(other_x_off+1:other_x_off+other_x_cnt, j, k)
        buf_off = buf_off + other_x_cnt
      end do
    end do

    call mpi_isend(send_buf(buf_base), buf_off - buf_base + 1, &
                   mpi_double_precision, send_task, xyyxt_tag, pmemd_comm, &
                   send_req(xy_idxlst_idx), err_code_mpi)

    buf_base = buf_base + siz_fft_mpi_buf1

  end do

  ! Do your own transposes, overlapped with i/o.

  call xy_transpose_my_data(x_dim, y_dim, xyz_data, yxz_data)

  ! Wait on and process the pending receive requests:

  do j = 1, fft_blks_dim1 - 1

    call mpi_waitany(fft_blks_dim1 - 1, recv_req, xy_idxlst_idx, irecv_stat, &
                     err_code_mpi)
    other_grid_idx1 = fft_xy_idxlst(xy_idxlst_idx)
    recv_task = fft_task_grid(other_grid_idx1, my_grid_idx2)
    call do_xy_transpose_recv(y_dim, other_grid_idx1, &
                              recv_buf(1, xy_idxlst_idx), yxz_data)

  end do

  ! Wait for all sends to complete:

  call mpi_waitall(fft_blks_dim1 - 1, send_req, isend_stat, err_code_mpi)

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
!              
!*******************************************************************************

subroutine do_xy_transpose_recv(y_dim, other_grid_idx1, recv_buf, yxz_data)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: y_dim
  integer               :: other_grid_idx1
  double precision      :: recv_buf(fft_x_cnts(my_grid_idx1) * 2 *&
                                    fft_y_cnts1(other_grid_idx1) * &
                                    fft_z_cnts(my_grid_idx2))
  double precision      :: yxz_data(2, y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))

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
  do k = 1, z_cnt
    do j = 1, other_y_cnt
      do i = 1, my_x_cnt
        buf_idx = buf_idx + 1
        yxz_data(1, other_y_off + j, i, k) = recv_buf(buf_idx)
        buf_idx = buf_idx + 1
        yxz_data(2, other_y_off + j, i, k) = recv_buf(buf_idx)
      end do
    end do
  end do
 
  return

end subroutine do_xy_transpose_recv

!*******************************************************************************
!
! Internal Subroutine:  xy_transpose_my_data
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine xy_transpose_my_data(x_dim, y_dim, xyz_data, yxz_data)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: x_dim
  integer               :: y_dim
  double precision      :: xyz_data(2, x_dim, &
                                    fft_y_cnts1(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
  double precision      :: yxz_data(2, y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))

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
        yxz_data(:, my_y_off + j, i, k) = xyz_data(:, my_x_off + i, j, k)
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

  implicit none

! Formal arguments:

  integer               :: y_dim
  double precision      :: y_runs(2, y_dim, &
                                  fft_x_cnts(my_grid_idx1), &
                                  fft_z_cnts(my_grid_idx2))

! Local variables:

  integer               :: j, k

  do k = 1, fft_z_cnts(my_grid_idx2)
    do j = 1, fft_x_cnts(my_grid_idx1)

#ifndef MIDPOINT_SPDP /*this will not be used if MIDPOINT_SPDP is defined*/
      call fft1d_forward(fft_y_hdl, y_runs(1, 1, j, k))
#endif

    end do
  end do

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

  implicit none

! Formal arguments:

  integer               :: y_dim
  integer               :: z_dim
  double precision      :: yxz_data(2 * y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
  double precision      :: zxy_data(2 * z_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_y_cnts2(my_grid_idx2))
  double precision      :: send_buf(siz_fft_mpi_buf2 * (fft_blks_dim2 - 1))
  double precision      :: recv_buf(siz_fft_mpi_buf2, (fft_blks_dim2 - 1))

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
  integer               :: recv_req(fft_blks_dim2 - 1)
  integer               :: send_req(fft_blks_dim2 - 1)

#ifdef COMM_TIME_TEST
  call start_test_timer(10, 'do_yz_transpose', 0)
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

    call mpi_irecv(recv_buf(1, yz_idxlst_idx), numval, mpi_double_precision, &
                   recv_task, yzzyt_tag, pmemd_comm, &
                   recv_req(yz_idxlst_idx), err_code_mpi) 
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
                   mpi_double_precision, send_task, yzzyt_tag, pmemd_comm, &
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
  double precision      :: recv_buf(2 * fft_y_cnts2(my_grid_idx2) * &
                                    fft_x_cnts(my_grid_idx1) * &
                                    fft_z_cnts(other_grid_idx2))
  double precision      :: zxy_data(2, z_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_y_cnts2(my_grid_idx2))

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
  double precision      :: yxz_data(2, y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
  double precision      :: zxy_data(2, z_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_y_cnts2(my_grid_idx2))

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

  implicit none

! Formal arguments:

  integer               :: z_dim
  double precision      :: z_runs(2, z_dim, &
                                  fft_x_cnts(my_grid_idx1), &
                                  fft_y_cnts2(my_grid_idx2))

! Local variables:

  integer               :: j, k

  do k = 1, fft_y_cnts2(my_grid_idx2)
    do j = 1, fft_x_cnts(my_grid_idx1)

#ifndef MIDPOINT_SPDP /*this will not be used if MIDPOINT_SPDP is defined*/
      call fft1d_forward(fft_z_hdl, z_runs(1, 1, j, k))
#endif

    end do
  end do

  return
      
end subroutine do_forward_z_fft

!*******************************************************************************
!
! Subroutine:  blk_fft3drc_back
!
! Description:  <TBS>
!              
!*******************************************************************************

subroutine blk_fft3drc_back(zxy_data, xyz_data, x_dim, y_dim, z_dim)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer                               :: x_dim, y_dim, z_dim

  double precision, intent(in out)      :: zxy_data(2, z_dim, &
                                                    fft_x_cnts(my_grid_idx1), &
                                                    fft_y_cnts2(my_grid_idx2))

  double precision, intent(out)         :: xyz_data(2, x_dim, &
                                                    fft_y_cnts1(my_grid_idx1), &
                                                    fft_z_cnts(my_grid_idx2))

! Local variables:

  double precision      :: yxz_data(2, y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))


  call do_backward_z_fft(z_dim, zxy_data)

  call do_zy_transpose(y_dim, z_dim, zxy_data, yxz_data, &
                       dbl_mpi_send_buf, dbl_mpi_recv_buf)

  call do_backward_y_fft(y_dim, yxz_data)

  call do_yx_transpose(x_dim, y_dim, yxz_data, xyz_data, &
                       dbl_mpi_send_buf, dbl_mpi_recv_buf)

  call do_backward_x_fft(x_dim, xyz_data)

  return

end subroutine blk_fft3drc_back

!*******************************************************************************
!
! Subroutine:  do_backward_z_fft
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine do_backward_z_fft(z_dim, z_runs)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: z_dim
  double precision      :: z_runs(2, z_dim, &
                                  fft_x_cnts(my_grid_idx1), &
                                  fft_y_cnts2(my_grid_idx2))

! Local variables:

  integer               :: j, k

  do k = 1, fft_y_cnts2(my_grid_idx2)
    do j = 1, fft_x_cnts(my_grid_idx1)

#ifndef MIDPOINT_SPDP /*this will not be used if MIDPOINT_SPDP is defined*/
      call fft1d_back(fft_z_hdl, z_runs(1, 1, j, k))
#endif

    end do
  end do

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
  double precision      :: zxy_data(2 * z_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_y_cnts2(my_grid_idx2))
  double precision      :: yxz_data(2 * y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
  double precision      :: send_buf(siz_fft_mpi_buf2 * (fft_blks_dim2 - 1))
  double precision      :: recv_buf(siz_fft_mpi_buf2, (fft_blks_dim2 - 1))

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
  integer               :: recv_req(fft_blks_dim2 - 1)
  integer               :: send_req(fft_blks_dim2 - 1)

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

    call mpi_irecv(recv_buf(1, yz_idxlst_idx), numval, mpi_double_precision, &
                   recv_task, zyyzt_tag, pmemd_comm, &
                   recv_req(yz_idxlst_idx), err_code_mpi) 
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
                   mpi_double_precision, send_task, zyyzt_tag, pmemd_comm, &
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
  double precision      :: recv_buf(2 * fft_z_cnts(my_grid_idx2) * &
                                    fft_x_cnts(my_grid_idx1) * &
                                    fft_y_cnts2(other_grid_idx2))
  double precision      :: yxz_data(2, y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))

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
  double precision      :: zxy_data(2, z_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_y_cnts2(my_grid_idx2))
  double precision      :: yxz_data(2, y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))

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
  double precision      :: y_runs(2, y_dim, &
                                  fft_x_cnts(my_grid_idx1), &
                                  fft_z_cnts(my_grid_idx2))

! Local variables:

  integer               :: j, k

  do k = 1, fft_z_cnts(my_grid_idx2)
    do j = 1, fft_x_cnts(my_grid_idx1)

#ifndef MIDPOINT_SPDP /*this will not be used if MIDPOINT_SPDP is defined*/
      call fft1d_back(fft_y_hdl, y_runs(1, 1, j, k))
#endif

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

subroutine do_yx_transpose(x_dim, y_dim, yxz_data, xyz_data, send_buf, recv_buf)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: x_dim
  integer               :: y_dim
  double precision      :: yxz_data(2 * y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
  double precision      :: xyz_data(2 * x_dim, &
                                    fft_y_cnts1(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
  double precision      :: send_buf(siz_fft_mpi_buf1 * (fft_blks_dim1 - 1))
  double precision      :: recv_buf(siz_fft_mpi_buf1, (fft_blks_dim1 - 1))

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
  integer               :: recv_req(fft_blks_dim1 - 1)
  integer               :: send_req(fft_blks_dim1 - 1)

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
    other_x_cnt = fft_y_cnts1(other_grid_idx1)

    recv_task = fft_task_grid(other_grid_idx1, my_grid_idx2)

    numval = 2 * my_y_cnt * other_x_cnt * z_cnt

    call mpi_irecv(recv_buf(1, xy_idxlst_idx), numval, mpi_double_precision, &
                   recv_task, yxxyt_tag, pmemd_comm, &
                   recv_req(xy_idxlst_idx), err_code_mpi) 
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
                   mpi_double_precision, send_task, yxxyt_tag, pmemd_comm, &
                   send_req(xy_idxlst_idx), err_code_mpi)

    buf_base = buf_base + siz_fft_mpi_buf1

  end do

  ! Do your own transposes, overlapped with i/o.

  call yx_transpose_my_data(x_dim, y_dim, yxz_data, xyz_data)

  ! Wait on and process the pending receive requests:

  do j = 1, fft_blks_dim1 - 1

    call mpi_waitany(fft_blks_dim1 - 1, recv_req, xy_idxlst_idx, irecv_stat, &
                     err_code_mpi)
    other_grid_idx1 = fft_xy_idxlst(xy_idxlst_idx)
    recv_task = fft_task_grid(other_grid_idx1, my_grid_idx2)
    call do_yx_transpose_recv(x_dim, other_grid_idx1, &
                              recv_buf(1, xy_idxlst_idx), xyz_data)

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

subroutine do_yx_transpose_recv(x_dim, other_grid_idx1, recv_buf, xyz_data)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: x_dim
  integer               :: other_grid_idx1
  double precision      :: recv_buf(fft_y_cnts1(my_grid_idx1) * 2 *&
                                    fft_x_cnts(other_grid_idx1) * &
                                    fft_z_cnts(my_grid_idx2))
  double precision      :: xyz_data(2, x_dim, &
                                    fft_y_cnts1(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))

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
      do i = 1, my_y_cnt
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

subroutine yx_transpose_my_data(x_dim, y_dim, yxz_data, xyz_data)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: x_dim
  integer               :: y_dim
  double precision      :: yxz_data(2, y_dim, &
                                    fft_x_cnts(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))
  double precision      :: xyz_data(2, x_dim, &
                                    fft_y_cnts1(my_grid_idx1), &
                                    fft_z_cnts(my_grid_idx2))

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
        xyz_data(:, my_x_off + i, j, k) = yxz_data(:, my_y_off + j, i, k)
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

subroutine do_backward_x_fft(x_dim, x_runs)

  implicit none

! Formal arguments:

  integer               :: x_dim
  double precision      :: x_runs(2, x_dim, &
                                  fft_y_cnts1(my_grid_idx1), &
                                  fft_z_cnts(my_grid_idx2))
! Local variables:

  double precision      :: a, b, c, d
  double precision      :: buf(2, x_dim - 1)
  integer               :: i, j, k
  integer               :: x_lim                ! indexing limit

! complex to real:

  x_lim = x_dim + 1

  do k = 1, fft_z_cnts(my_grid_idx2)

    do j = 1, fft_y_cnts1(my_grid_idx1)

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

#ifndef MIDPOINT_SPDP /*this will not be used if MIDPOINT_SPDP is defined*/
      call fft1d_back(fft_x_hdl, buf(1, 1))
#endif

      do i = 1, x_dim - 1
        x_runs(:, i, j, k) = buf(:, i)
      end do

    end do

  end do

  return

end subroutine do_backward_x_fft

#endif /* MPI */

end module pme_blk_fft_mod
