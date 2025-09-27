#include "copyright.i"

!*******************************************************************************
!
! Module:  img_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module img_mod

  use gbl_datatypes_mod

!precision setting for SPDP, DPDP, SPSP
#include "include_precision.i"

  implicit none

  double precision, allocatable, save           :: gbl_img_crd(:,:)
  ! qterm is atom charge under pme or gb:

  double precision, allocatable, save           :: gbl_img_qterm(:)
  integer, allocatable, save                    :: gbl_img_iac(:)

  integer, allocatable, save                    :: gbl_atm_img_map(:)
  integer, allocatable, save                    :: gbl_img_atm_map(:)

#ifdef MPI
  ! used images are "1" in the map; unused images are "0"
  integer(byte), allocatable, save              :: gbl_used_img_map(:)
  integer(byte), allocatable, save              :: gbl_used_img_ips_map(:)
  integer, allocatable, save                    :: gbl_mapped_img_lst(:)
  integer, allocatable, save                    :: gbl_used_img_lst(:)
  integer, allocatable, save                    :: gbl_used_img_ips_lst(:)
  integer, allocatable, save                    :: gbl_owned_imgs_tbl(:,:)
#endif /* MPI */
  integer, allocatable, save                    :: gbl_excl_img_flags(:)
  double precision, save                        :: gbl_tranvec(1:3,0:17)
  ! my_img_* covers the range of images you "own"; ie. forces are accumulated
  ! in the local process for those images.  These are assigned WITHOUT wrapping.
  ! used_img_* covers the range of images you "use"; ie. you may need
  ! coordinate information and will report nonbonded force information to the
  ! owner for some images in this range. This range may wrap through natom to 1.

  integer, save         :: my_img_lo, my_img_hi

#ifdef MPI
  integer, save         :: gbl_mapped_img_cnt = 0
  integer, save         :: gbl_used_img_cnt = 0
  integer, save         :: gbl_used_img_ips_cnt = 0
#endif

contains

!*******************************************************************************
!
! Subroutine:  alloc_img_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_img_mem(atm_cnt, num_ints, num_reals)

  use parallel_dat_mod
  use pmemd_lib_mod
  use mdin_ctrl_dat_mod, only : ips

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer                       :: atm_cnt
  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed

#ifndef MPI
  ! In an MPI version, this stuff will be set elsewhere.

  my_img_lo = 1
  my_img_hi = atm_cnt
#endif

  ! The gbl_owned_imgs_tbl has one extra element at the end that is not
  ! currently used.

   if (ips .gt. 0)then
  allocate(gbl_img_crd(3, atm_cnt), &
           gbl_img_qterm(atm_cnt), &
           gbl_img_iac(atm_cnt), &
           gbl_atm_img_map(atm_cnt), &
           gbl_img_atm_map(atm_cnt), &
#ifdef MPI
           gbl_used_img_map(atm_cnt), &
           gbl_used_img_ips_map(atm_cnt), &
           gbl_mapped_img_lst(atm_cnt), &
           gbl_used_img_lst(atm_cnt), &
           gbl_used_img_ips_lst(atm_cnt), &
           gbl_owned_imgs_tbl(2, 0:numtasks), &
#endif /* MPI */
           gbl_excl_img_flags(atm_cnt), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_reals = num_reals + size(gbl_img_crd) + &
                          size(gbl_img_qterm)

  num_ints = num_ints + size(gbl_img_iac) + &
                        size(gbl_atm_img_map) + &
                        size(gbl_img_atm_map) + &
#ifdef MPI
                        size(gbl_used_img_map) / 4 + &
                        size(gbl_used_img_ips_map) / 4 + &
                        size(gbl_mapped_img_lst) + &
                        size(gbl_used_img_lst) + &
                        size(gbl_used_img_ips_lst) + &
                        size(gbl_owned_imgs_tbl) + &
#endif
                        size(gbl_excl_img_flags)

  gbl_img_crd(:,:) = 0.d0
  gbl_img_qterm(:) = 0.d0
  gbl_img_iac(:) = 0

  gbl_atm_img_map(:) = 0
  gbl_img_atm_map(:) = 0
#ifdef MPI
  gbl_used_img_map(:) = 0
  gbl_used_img_ips_map(:) = 0
  gbl_mapped_img_lst(:) = 0
  gbl_used_img_lst(:) = 0
  gbl_used_img_ips_lst(:) = 0
  gbl_owned_imgs_tbl(:,:) = 0
#endif
  gbl_excl_img_flags(:) = 0
  else
  allocate(gbl_img_crd(3, atm_cnt), &
           gbl_img_qterm(atm_cnt), &
           gbl_img_iac(atm_cnt), &
           gbl_atm_img_map(atm_cnt), &
           gbl_img_atm_map(atm_cnt), &
#ifdef MPI
           gbl_used_img_map(atm_cnt), &
           gbl_mapped_img_lst(atm_cnt), &
           gbl_used_img_lst(atm_cnt), &
           gbl_owned_imgs_tbl(2, 0:numtasks), &
#endif /* MPI */
           gbl_excl_img_flags(atm_cnt), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_reals = num_reals + size(gbl_img_crd) + &
                          size(gbl_img_qterm)

  num_ints = num_ints + size(gbl_img_iac) + &
                        size(gbl_atm_img_map) + &
                        size(gbl_img_atm_map) + &
#ifdef MPI
                        size(gbl_used_img_map) / 4 + &
                        size(gbl_mapped_img_lst) + &
                        size(gbl_used_img_lst) + &
                        size(gbl_owned_imgs_tbl) + &
#endif
                        size(gbl_excl_img_flags)

  gbl_img_crd(:,:) = 0.d0
  gbl_img_qterm(:) = 0.d0
  gbl_img_iac(:) = 0

  gbl_atm_img_map(:) = 0
  gbl_img_atm_map(:) = 0
#ifdef MPI
  gbl_used_img_map(:) = 0
  gbl_mapped_img_lst(:) = 0
  gbl_used_img_lst(:) = 0
  gbl_owned_imgs_tbl(:,:) = 0
#endif
  gbl_excl_img_flags(:) = 0
  endif

  return

end subroutine alloc_img_mem

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  bcast_img_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine bcast_img_dat(atm_cnt)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer       :: atm_cnt

! Local variables:

  integer       :: num_ints, num_reals  ! returned values discarded

  ! Nothing to broadcast.  We just allocate storage in the non-master nodes.

  if (.not. master) then
    num_ints = 0
    num_reals = 0
    call alloc_img_mem(atm_cnt, num_ints, num_reals)
  end if

  ! The allocated data is not initialized from the master node.

  return

end subroutine bcast_img_dat
#endif

!*******************************************************************************
!
! Subroutine:  fill_tranvec
!
! Description: <TBS>
!              
!*******************************************************************************
subroutine fill_tranvec(tranvec)
  use pbc_mod
  implicit none

  double precision      :: tranvec(3, 18) ! this is (1:3,0:17) externally.
  integer               :: iv, i0, i1, i2, i3

! This works for both orthogonal and nonorthogonal unit cells.

  iv = 0
  do i3 = -1, 0
    do i2 = -1, 1
      do i1 = -1, 1
        iv = iv + 1
        do i0 = 1, 3
          tranvec(i0, iv) = i1 * ucell(i0, 1) + &
                            i2 * ucell(i0, 2) + &
                            i3 * ucell(i0, 3)
        end do
      end do
    end do
  end do
  return
end subroutine fill_tranvec

#ifdef MPI
!*******************************************************************************
!
! Subroutine:   init_used_img_map
!
! Description:  Per name; initialization includes setting currently owned
!               images as "used"
!
!*******************************************************************************

subroutine init_used_img_map(used_img_map, used_img_cnt, used_img_lst)

  implicit none

! Formal arguments:

  integer(byte)         :: used_img_map(*)
  integer               :: used_img_cnt
  integer               :: used_img_lst(*)

! Local variables:

  integer               :: lst_idx, img_idx

  ! First clear the last used_img_map data structures.

  do lst_idx = 1, used_img_cnt
    img_idx = used_img_lst(lst_idx)
    used_img_map(img_idx) = 0
  end do

  used_img_cnt = 0

  ! Now go ahead and mark all owned images as used.

  do img_idx = my_img_lo, my_img_hi
    used_img_map(img_idx) = 1
    used_img_cnt = used_img_cnt + 1
    used_img_lst(used_img_cnt) = img_idx
  end do
    
  return

end subroutine init_used_img_map

!*******************************************************************************
!
! Subroutine:   init_used_img_ips_map
!
! Description:  Per name; initialization includes setting currently owned
!               images as "used"
!
!*******************************************************************************

subroutine init_used_img_ips_map(used_img_ips_map, used_img_ips_cnt, used_img_ips_lst)

  implicit none

! Formal arguments:

  integer(byte)         :: used_img_ips_map(*)
  integer               :: used_img_ips_cnt
  integer               :: used_img_ips_lst(*)

! Local variables:

  integer               :: lst_idx, img_idx

  ! First clear the last used_img_map data structures.

  do lst_idx = 1, used_img_ips_cnt
    img_idx = used_img_ips_lst(lst_idx)
    used_img_ips_map(img_idx) = 0
  end do

  used_img_ips_cnt = 0

  ! Now go ahead and mark all owned images as used.

  do img_idx = my_img_lo, my_img_hi
    used_img_ips_map(img_idx) = 1
    used_img_ips_cnt = used_img_ips_cnt + 1
    used_img_ips_lst(used_img_ips_cnt) = img_idx
  end do
    
  return

end subroutine init_used_img_ips_map
#endif /* MPI */

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  divide_images_evenly
!
! Description:  Set up image boundaries for parallel processing.
!
!*******************************************************************************

subroutine divide_images_evenly(img_cnt, owned_imgs_tbl, &
                                frc_ene_numtasks, frc_ene_task_lst, &
                                write_log, step_ctr)

  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use pmemd_lib_mod
  use mdin_ewald_dat_mod

  implicit none

! Formal arguments:

  integer       :: img_cnt
  integer       :: owned_imgs_tbl(2, 0:numtasks)
  integer       :: frc_ene_numtasks
  integer       :: frc_ene_task_lst(frc_ene_numtasks)
  logical       :: write_log
  integer       :: step_ctr

! Local variables:

  double precision      :: fraction
  double precision      :: next_div
  integer               :: lst_idx
  integer               :: taskid
  integer               :: next_used_img_cnt
  integer               :: used_img_cnt

  ! The fraction is guaranteed .gt. 0 by an earlier check...

  fraction = dble(img_cnt) / dble(frc_ene_numtasks)
  used_img_cnt = 0
  next_div = 0.d0

  ! For highscaling setup, there may be some tasks without frc ene workload:

  if (frc_ene_numtasks .ne. numtasks) owned_imgs_tbl(2, 0:numtasks - 1) = 0

  do lst_idx = 1, frc_ene_numtasks - 1
    taskid = frc_ene_task_lst(lst_idx)
    next_div = next_div + fraction
    next_used_img_cnt = int(next_div)
    owned_imgs_tbl(2, taskid) = next_used_img_cnt - used_img_cnt
    used_img_cnt = next_used_img_cnt
  end do

  taskid = frc_ene_task_lst(frc_ene_numtasks)
  owned_imgs_tbl(2, taskid) = img_cnt - used_img_cnt

  used_img_cnt = 0

  do taskid = 0, numtasks - 1
    owned_imgs_tbl(1, taskid) = used_img_cnt
    used_img_cnt = used_img_cnt + owned_imgs_tbl(2, taskid)
  end do

  my_img_lo = owned_imgs_tbl(1, mytaskid) + 1
  my_img_hi = my_img_lo + owned_imgs_tbl(2, mytaskid) - 1

  log_owned_img_cnt = dble(my_img_hi - my_img_lo + 1)

  if (master) then
    if (loadbal_verbose .gt. 0 .or. write_log) then
        write(logfile, 10) &
          'Image Distribution at run step ', step_ctr, ':'
        write(logfile, 20) 'Count of images assigned to each task:'
        write(logfile, 30) (owned_imgs_tbl(2, taskid), taskid = 0, numtasks - 1)
    end if
  end if

10 format(/, a, i7, a /)
20 format('  ', a)
30 format('    ', 8i9)

  return

end subroutine divide_images_evenly
#endif /* MPI */

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  divide_images_excl_recip
!
! Description:  Set up image boundaries for parallel processing.
!
!*******************************************************************************

subroutine divide_images_excl_recip(img_cnt, owned_imgs_tbl, &
                                    frc_ene_numtasks, frc_ene_task_lst, &
                                    is_recip_task, recip_numtasks, &
                                    write_log, step_ctr)

  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use pmemd_lib_mod
  use mdin_ewald_dat_mod

  implicit none

! Formal arguments:

  integer       :: img_cnt
  integer       :: owned_imgs_tbl(2, 0:numtasks)
  integer       :: frc_ene_numtasks
  integer       :: frc_ene_task_lst(frc_ene_numtasks)
  logical       :: is_recip_task(0:)
  integer       :: recip_numtasks
  logical       :: write_log
  integer       :: step_ctr

! Local variables:

  double precision      :: fraction
  double precision      :: next_div
  integer               :: lst_idx, dirlst_idx
  integer               :: taskid
  integer               :: next_used_img_cnt
  integer               :: used_img_cnt
  integer               :: dir_numtasks
  logical               :: saved_master_recip_stat
  integer               :: dir_task_lst(frc_ene_numtasks - recip_numtasks)

  saved_master_recip_stat = is_recip_task(0)

  dir_numtasks = frc_ene_numtasks - recip_numtasks

  if (excl_master .ne. 0) then
    is_recip_task(0) = .true. ! done to exclude master...
    dir_numtasks = dir_numtasks - 1
  end if

  ! Make the direct task list.  These tasks DO NOT do recip work...

  dirlst_idx = 0
  do lst_idx = 1, frc_ene_numtasks
    taskid = frc_ene_task_lst(lst_idx)
    if (is_recip_task(taskid)) cycle
    dirlst_idx = dirlst_idx + 1
    dir_task_lst(dirlst_idx) = taskid
  end do

  ! The fraction is guaranteed .gt. 0 by an earlier check...

  fraction = dble(img_cnt) / dble(dir_numtasks)
  used_img_cnt = 0
  next_div = 0.d0

  ! All the reciprocal tasks will not have a workload...

  owned_imgs_tbl(2, 0:numtasks - 1) = 0

  do dirlst_idx = 1, dir_numtasks - 1
    taskid = dir_task_lst(dirlst_idx)
    next_div = next_div + fraction
    next_used_img_cnt = int(next_div)
    owned_imgs_tbl(2, taskid) = next_used_img_cnt - used_img_cnt
    used_img_cnt = next_used_img_cnt
  end do

  taskid = dir_task_lst(dir_numtasks)
  owned_imgs_tbl(2, taskid) = img_cnt - used_img_cnt

  used_img_cnt = 0

  do taskid = 0, numtasks - 1
    owned_imgs_tbl(1, taskid) = used_img_cnt
    used_img_cnt = used_img_cnt + owned_imgs_tbl(2, taskid)
  end do

  my_img_lo = owned_imgs_tbl(1, mytaskid) + 1
  my_img_hi = my_img_lo + owned_imgs_tbl(2, mytaskid) - 1

  log_owned_img_cnt = dble(my_img_hi - my_img_lo + 1)

  if (master) then
    if (loadbal_verbose .gt. 0 .or. write_log) then
        write(logfile, 10) &
          'Image Distribution at run step ', step_ctr, ':'
        write(logfile, 20) 'Count of images assigned to each task:'
        write(logfile, 30) (owned_imgs_tbl(2, taskid), taskid = 0, numtasks - 1)
    end if
  end if

  is_recip_task(0) = saved_master_recip_stat

10 format(/, a, i7, a /)
20 format('  ', a)
30 format('    ', 8i9)

  return

end subroutine divide_images_excl_recip
#endif /* MPI */

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  divide_images_recip_biased
!
! Description:  Set up image boundaries for parallel processing.
!
!*******************************************************************************

subroutine divide_images_recip_biased(img_cnt, recip_workload, owned_imgs_tbl, &
                                      is_recip_task, recip_numtasks, &
                                      frc_ene_numtasks, frc_ene_task_lst, &
                                      write_log, step_ctr)

  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use pmemd_lib_mod
  use mdin_ewald_dat_mod

  implicit none

! Formal arguments:

  integer               :: img_cnt
  double precision      :: recip_workload
  integer               :: owned_imgs_tbl(2, 0:numtasks)
  logical               :: is_recip_task(0:)
  integer               :: recip_numtasks
  integer               :: frc_ene_numtasks
  integer               :: frc_ene_task_lst(frc_ene_numtasks)
  logical               :: write_log
  integer               :: step_ctr

! Local variables:

  double precision              :: nonrecip_frac
  double precision              :: recip_frac
  double precision              :: next_div
  integer                       :: lst_idx
  integer                       :: taskid
  integer                       :: next_used_img_cnt
  integer                       :: used_img_cnt


  ! If all tasks have a reciprocal workload, there is no bias in the
  ! assignment of workload.

  if (recip_numtasks .eq. frc_ene_numtasks) then
    call divide_images_evenly(img_cnt, owned_imgs_tbl, &
                              frc_ene_numtasks, frc_ene_task_lst, &
                              write_log, step_ctr)
    return
  end if

  ! The direct workload is 1.d0.  We first apportion the direct workload to
  ! give all tasks an equal workload, based on the assumption that the input
  ! recip workload, expressed as a fraction of the direct workload, is correct.
  ! Then, if the direct workload assigned to the tasks handling a reciprocal
  ! workload is less than 10% of the total direct workload, we readjust the
  ! direct workload, distributing 10% of it to the tasks that also handle a
  ! reciprocal workload.  This adjustment will place a higher total workload
  ! on the reciprocal tasks, but it will be quickly adjusted down as
  ! appropriate.  This algorithm is used because direct force load balancing
  ! does not adapt quickly if the initially assigned direct force load is
  ! near 0.

  nonrecip_frac = (1.d0 + recip_workload) / dble(frc_ene_numtasks)
  recip_frac = nonrecip_frac - (recip_workload / dble(recip_numtasks))

  ! Division by 0 insured to not occur by above check for
  ! recip_numtasks .eq. frc_ene_numtasks.

  if (recip_frac * dble(recip_numtasks) / &
      (nonrecip_frac * dble(frc_ene_numtasks - recip_numtasks)) .lt. 0.1d0) then
    recip_frac = 0.1d0 / dble(recip_numtasks)
    nonrecip_frac = 0.9d0 / dble(frc_ene_numtasks - recip_numtasks)
  end if

  recip_frac = recip_frac * dble(img_cnt)
  nonrecip_frac = nonrecip_frac * dble(img_cnt)

  ! For highscaling setup, there may be some tasks without frc ene workload:

  if (frc_ene_numtasks .ne. numtasks) owned_imgs_tbl(2, 0:numtasks - 1) = 0
  
  used_img_cnt = 0
  next_div = 0.d0

  do lst_idx = 1, frc_ene_numtasks - 1
    taskid = frc_ene_task_lst(lst_idx)
    if (.not. is_recip_task(taskid)) then
      next_div = next_div + nonrecip_frac
    else
      next_div = next_div + recip_frac
    end if
    next_used_img_cnt = int(next_div)
    owned_imgs_tbl(2, taskid) = next_used_img_cnt - used_img_cnt
    used_img_cnt = next_used_img_cnt
  end do

  taskid = frc_ene_task_lst(frc_ene_numtasks)
  owned_imgs_tbl(2, taskid) = img_cnt - used_img_cnt

  used_img_cnt = 0

  do taskid = 0, numtasks - 1
    owned_imgs_tbl(1, taskid) = used_img_cnt
    used_img_cnt = used_img_cnt + owned_imgs_tbl(2, taskid)
  end do

  my_img_lo = owned_imgs_tbl(1, mytaskid) + 1
  my_img_hi = my_img_lo + owned_imgs_tbl(2, mytaskid) - 1

  log_owned_img_cnt = dble(my_img_hi - my_img_lo + 1)

10 format(/, a, i7, a /)
20 format('  ', a, i4, a)
22 format('  ', a, f5.1, a)
30 format('  ', a)
40 format('    ', 8i9)

  if (master) then

    if (loadbal_verbose .gt. 2) then
      write(logfile, 10) &
        'Image Distribution Bias Data at run step ', step_ctr, ':'
      write(logfile, 20) 'Recip force workload assigned to ', recip_numtasks, &
                         ' tasks'
      write(logfile, 22) 'Recip force workload estimate is ', &
            recip_workload * 100.d0, ' percent'
    end if

    if (loadbal_verbose .gt. 0 .or. write_log) then
      write(logfile, 10) &
        'Image Distribution at run step ', step_ctr, ':'
      write(logfile, 30) 'Count of images assigned to each task:'
      write(logfile, 40) (owned_imgs_tbl(2, taskid), taskid = 0, numtasks - 1)
    end if

  end if

  return

end subroutine divide_images_recip_biased
#endif /* MPI */

end module img_mod
