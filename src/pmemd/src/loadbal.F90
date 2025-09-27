#include "copyright.i"

!*******************************************************************************
!
! Module:  loadbal_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module loadbal_mod

  implicit none

#ifdef MPI

  ! Load balancing timers and other data structures:

  integer, private, save        :: start_sec, start_usec
  integer, save                 :: elapsed_100usec_dirfrc = 0
  integer, save                 :: elapsed_100usec_recipfrc = 0
  integer, save                 :: elapsed_100usec_atmowner = 0 ! bnd,ang,dihed,
                                                                ! other atm work
  integer, save                 :: elapsed_100usec_atmuser = 0  ! used atm and
                                                                ! img work

  ! Some counters start at -1 to inhibit 1st pass redistribution.

  integer, private, save        :: img_redist_ctr = -1  
  integer, private, save        :: img_redist_trigger = 1 ! 1,2,4,8,16,32,64
  integer, parameter            :: img_redist_trigger_max = 64
  integer, parameter            :: img_redist_trigger_divisor = 4
  double precision, parameter   :: img_redist_time_tol = 0.01d0

  logical, save                 :: atm_redist_needed = .false.
  integer, private, save        :: atm_redist_ctr = -1

  integer, private, save        :: loadbal_step_ctr = 0
  integer, private, save        :: last_loadbal_step_ctr = 0

  integer, allocatable, save    :: gbl_loadbal_node_dat(:,:)

  ! All these are used for slab or block fft's:

  logical, save                 :: fft_redist_enabled = .false.
  logical, save                 :: fft_redist_needed = .false.
  integer, private, save        :: fft_redist_ctr = -1
  integer, private, save        :: fft_redist_trigger = 1 ! 1,2,3,...

contains

!*******************************************************************************
!
! Subroutine:  start_loadbal_timer
!
! Description: Used by cit code.
!              
!*******************************************************************************

subroutine start_loadbal_timer

  implicit none

  call get_wall_time(start_sec, start_usec)

  return

end subroutine start_loadbal_timer

!*******************************************************************************
!
! Subroutine:  update_loadbal_timer
!
! Description: Used by cit code. Note - Only can run one loadbal timer at a
!              time (ie., the adjustable and nonadjustable timer start/stop
!              cycles must not overlap (but they would never need to)).
!              
!*******************************************************************************

subroutine update_loadbal_timer(elapsed_100usec)

  implicit none

! Formal arguments:

  integer               :: elapsed_100usec

! Local variables:

  integer               :: stop_sec, stop_usec
  integer               :: next_start_sec, next_start_usec
  integer               :: this_exec_time

  call get_wall_time(stop_sec, stop_usec)

  next_start_sec = stop_sec
  next_start_usec = stop_usec
    
  ! For next round:

  if (stop_usec .lt. start_usec) then
    stop_usec = stop_usec + 1000000
    stop_sec = stop_sec - 1
  end if

  this_exec_time = (stop_sec - start_sec) * 10000 + &
                   (stop_usec - start_usec) / 100

  elapsed_100usec = elapsed_100usec + this_exec_time

  start_sec = next_start_sec
  start_usec = next_start_usec

  return

end subroutine update_loadbal_timer

!*******************************************************************************
!
! Subroutine:  do_load_balancing
!
! Description: The main entry into load balancing code; should be called at
!              the beginning of each force evaluation; it actually only
!              executes if there is a new pairlist.
!*******************************************************************************

subroutine do_load_balancing(new_list, atm_cnt)

  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_mod
  use timers_mod

  implicit none

! Formal arguments:

  logical               :: new_list
  integer               :: atm_cnt

  loadbal_step_ctr = loadbal_step_ctr + 1

  if (.not. new_list) return

  atm_redist_ctr = atm_redist_ctr + 1

  if (fft_redist_needed) then

    fft_redist_ctr = fft_redist_ctr + 1

    if (fft_redist_ctr .ge. fft_redist_trigger) then
      call do_fft_redistribution(atm_cnt)
      fft_redist_ctr = 0
      fft_redist_trigger = fft_redist_trigger + 1
      atm_redist_ctr = 0
      img_redist_ctr = 0
      img_redist_trigger = 1
      call update_pme_time(fft_reassign_timer)
    else
      img_redist_ctr = img_redist_ctr + 1
      if (img_redist_ctr .ge. 1 .and. img_redist_ctr .ge. &
          img_redist_trigger * nrespa / img_redist_trigger_divisor) then
        call do_img_redistribution
        img_redist_ctr = 0
        img_redist_trigger = min(img_redist_trigger * 2, img_redist_trigger_max)
        call update_pme_time(img_reassign_timer)
      end if
    end if

  else if (atm_redist_needed) then

    call do_atm_redistribution(atm_cnt, .true., loadbal_step_ctr)
    atm_redist_needed = .false.
    img_redist_ctr = 0
    img_redist_trigger = 1
    call update_pme_time(atm_reassign_timer)

  else

    img_redist_ctr = img_redist_ctr + 1
    if (img_redist_ctr .ge. 1 .and. img_redist_ctr .ge. &
        img_redist_trigger * nrespa / img_redist_trigger_divisor) then
      call do_img_redistribution
      img_redist_ctr = 0
      img_redist_trigger = min(img_redist_trigger * 2, img_redist_trigger_max)
      call update_pme_time(img_reassign_timer)
    end if

  end if

  if (atm_redist_ctr .ge. atm_redist_freq) then
    atm_redist_needed = .true.
    atm_redist_ctr = 0
  end if

  return

end subroutine do_load_balancing

!*******************************************************************************
!
! Subroutine:  check_new_list_limit
!
! Description: Used by pme code. Note - Only can run one loadbal timer at a
!              time (ie., the adjustable and nonadjustable timer start/stop
!              cycles must not overlap (but they would never need to)).
!              
!*******************************************************************************

subroutine check_new_list_limit(new_list)

  implicit none

! Formal arguments:

  logical               :: new_list

! Local variables:
  integer, save         :: last_new_list_cnt = 0
  integer, save         :: last_new_list_limit = 16
  integer, parameter    :: max_last_new_list_limit = 32

  if (new_list) then
    last_new_list_cnt = 0
  else
    last_new_list_cnt = last_new_list_cnt + 1
    if (last_new_list_cnt .ge. last_new_list_limit) then
      last_new_list_cnt = 0
      last_new_list_limit = last_new_list_limit + 1
      last_new_list_limit = min(last_new_list_limit, max_last_new_list_limit)
      new_list = .true.
    end if
  end if

  return

end subroutine check_new_list_limit

!*******************************************************************************
!
! Subroutine:  do_fft_redistribution
!
! Description:  <TBS>
!              
!*******************************************************************************

subroutine do_fft_redistribution(atm_cnt)

  use img_mod
  use pme_recip_dat_mod
  use pme_blk_fft_mod
  use pme_slab_fft_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt

! Local variables:

  double precision              :: avg_total_time
  double precision              :: recip_nodes_time
  double precision              :: recip_time_per_slab
  double precision              :: recip_workload
  double precision              :: total_time
  logical                       :: not_distributed
  integer                       :: i
  integer                       :: node
  integer                       :: num_ints, num_reals  ! values are not kept.
  integer                       :: use_recip_workload
  integer                       :: min_numtasks

  ! We place the following data in static storage to keep the Myrinet MPI
  ! implementation from going ape over mpi buffers on the stack (causes
  ! problems between mpi and ifc, resulting in huge stalls in processing).

  integer, save                 :: my_node_dat(4)
  integer, save                 :: fft_redist_dat(3)

  num_ints = 0
  num_reals = 0
  use_recip_workload = 0

  ! Gather times from the other processes to the master process.  The
  ! reference to gbl_loadbal_node_dat is only used in the master process.
  ! The dat collected from each node are the direct, reciprocal, atm owner,
  ! and atm (and image) user categories.

  if (i_do_frc_ene) then
    my_node_dat(1) = elapsed_100usec_dirfrc
    my_node_dat(2) = elapsed_100usec_recipfrc
    my_node_dat(3) = elapsed_100usec_atmowner
    my_node_dat(4) = elapsed_100usec_atmuser
  else
    my_node_dat(1:4) = 0
  end if

  call mpi_gather(my_node_dat, 4, mpi_integer, gbl_loadbal_node_dat, 4, &
                  mpi_integer, 0, pmemd_comm, err_code_mpi)

  if (master) then

    total_time = 0.d0
    recip_nodes_time = 0.d0

    do i = 0, numtasks - 1
      ! We include atom user time as reciprocal time because in reciprocal
      ! tasks it is largely a reciprocal cost.
      if (is_recip_task(i)) then
        recip_nodes_time = recip_nodes_time + &
                           dble(gbl_loadbal_node_dat(2, i) + &
                                gbl_loadbal_node_dat(4, i))
      end if
      ! We intentionally omit list building time from the total time because
      ! it does not occur every cycle, is not completely variable with image
      ! count, and basically adds an undesirable element of randomness.
      total_time = total_time + &
                   dble(gbl_loadbal_node_dat(1, i) + &
                        gbl_loadbal_node_dat(2, i) + &
                        gbl_loadbal_node_dat(3, i) + &
                        gbl_loadbal_node_dat(4, i))
    end do

    if (block_fft .eq. 0) &
      recip_time_per_slab = recip_nodes_time / dble(fft_z_dim)
    
    avg_total_time = total_time / dble(frc_ene_task_cnt)
    recip_workload = recip_nodes_time / total_time

! BEGIN DBG
!   if (block_fft .eq. 0) &
!     write(logfile,*)'DBG: recip_time_per_slab*max_xy_slabs=',&
!                     recip_time_per_slab*dble(max_xy_slab_cnt)
!   write(logfile,*)'DBG: avg_total_time=',avg_total_time
!   write(logfile,*)'DBG: elapsed steps=',&
!                   loadbal_step_ctr - last_loadbal_step_ctr + 1
! END DBG

    if (block_fft .eq. 0) then
      if (recip_time_per_slab * dble(max_xy_slab_cnt) .gt. avg_total_time) then
        if (fft_redist_ctr .eq. 1) then
            not_distributed = .true.
            use_recip_workload = 1
        else
          if (max_xy_slab_cnt .eq. 1 .or. &
              max_zx_slab_cnt .eq. 1 .or. &
              recip_numtasks .eq. frc_ene_task_cnt) then

            not_distributed = .false.
            
            ! We have split the slabs up as far as we can.
          else
            not_distributed = .true.
          end if
        end if
      else
        not_distributed = .false.
      end if
    else
      if (recip_numtasks .eq. frc_ene_task_cnt) then
        not_distributed = .false.
      else if (recip_workload .ge. &
               dble(recip_numtasks) / dble(frc_ene_task_cnt)) then
        if (fft_redist_ctr .eq. 1) then
            not_distributed = .true.
            use_recip_workload = 1
        else
          if ((fft_blks_dim2+1) * (fft_blks_dim1) .gt. frc_ene_task_cnt .or. &
              fft_blks_dim2 .eq. max_fft_blks_dim2) then
            ! We cannot assign the reciprocal space workload to more tasks.
            not_distributed = .false.
          else
            not_distributed = .true.
          end if
        end if
      else
        not_distributed = .false.
      end if
    end if

    if (.not. not_distributed) then
      ! Negative value indicates balancing done; don't need to redistrib
      ! slabs, but DO adjust images and atoms using workload.
      recip_workload = -recip_workload
    end if

    ! Integerize the reciprocal workload, and send the other fft
    ! redistribution control variables to the slaves.

    fft_redist_dat(1) = int(recip_workload * 1000000.d0)
    fft_redist_dat(2) = fft_redist_ctr
    fft_redist_dat(3) = use_recip_workload

    call mpi_bcast(fft_redist_dat, 3, mpi_integer, 0, pmemd_comm, &
                   err_code_mpi)

    if (loadbal_verbose .gt. 1 .and. recip_workload .gt. 0.d0) then
      write(logfile, '(/,2x,a,f5.1,a)') 'Recip force workload estimate is ', &
                         recip_workload * 100.d0, ' percent'
    end if

  else

    call mpi_bcast(fft_redist_dat, 3, mpi_integer, 0, pmemd_comm, &
                   err_code_mpi)

    fft_redist_ctr = fft_redist_dat(2)
    use_recip_workload = fft_redist_dat(3)

  end if

  ! Recip workload recalc'd in master also to insure exact same result as
  ! in slaves.

  recip_workload = dble(fft_redist_dat(1))/1000000.d0


  if (recip_workload .gt. 0.d0) then
    if (block_fft .eq. 0) then
      if (use_recip_workload .ne. 0) then
        recip_numtasks = ceiling(dble(frc_ene_task_cnt) * recip_workload)
      else
        ! We must only execute this code if max_xy_slab_cnt .gt. 1, currently
        ! insured by .not. not_distributed --> recip_workload .lt. 0.d0
        recip_numtasks = fft_z_dim / (max_xy_slab_cnt - 1)
        if (mod(fft_z_dim, max_xy_slab_cnt - 1) .ne. 0) then
          recip_numtasks = recip_numtasks + 1
        end if
      end if

      call distribute_slab_fft_slabs(recip_numtasks, fft_redist_enabled, &
                                     loadbal_step_ctr, num_ints, num_reals)
    else
      if (use_recip_workload .ne. 0) then
        min_numtasks = ceiling(dble(frc_ene_task_cnt) * recip_workload)
      else
        min_numtasks = (fft_blks_dim2 + 1) * (fft_blks_dim1)
      end if

      call blk_fft_setup(min_numtasks, recip_numtasks, &
                         frc_ene_task_cnt, gbl_frc_ene_task_lst, &
                         loadbal_step_ctr, num_ints, num_reals)
    end if
  else
    recip_workload = -recip_workload
    fft_redist_needed = .false.
  end if

  if (recip_workload .gt. 0.6d0) excl_recip = 0 ! just in case...

  if (excl_recip .eq. 0) then
    if (recip_numtasks .eq. frc_ene_task_cnt) then
      call divide_images_evenly(atm_cnt, gbl_owned_imgs_tbl, &
                                frc_ene_task_cnt, gbl_frc_ene_task_lst, &
                                .false., loadbal_step_ctr)
    else
      call divide_images_recip_biased(atm_cnt, recip_workload, &
                                      gbl_owned_imgs_tbl, is_recip_task, &
                                      recip_numtasks, &
                                      frc_ene_task_cnt, gbl_frc_ene_task_lst, &
                                      .false., loadbal_step_ctr)
    end if
  else
    call divide_images_excl_recip(atm_cnt, gbl_owned_imgs_tbl, &
                                  frc_ene_task_cnt, gbl_frc_ene_task_lst, &
                                  is_recip_task, recip_numtasks, &
                                  .false., loadbal_step_ctr)
  end if
  
  call do_atm_redistribution(atm_cnt, .false., loadbal_step_ctr)

  elapsed_100usec_dirfrc = 0
  elapsed_100usec_recipfrc = 0
  elapsed_100usec_atmowner = 0
  elapsed_100usec_atmuser = 0
  last_loadbal_step_ctr = loadbal_step_ctr

  return

end subroutine do_fft_redistribution

!*******************************************************************************
!
! Subroutine:  do_img_redistribution
!
! Description:  <TBS>
!              
!*******************************************************************************

subroutine do_img_redistribution

  use img_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use parallel_mod
  use mdin_ewald_dat_mod
  use pmemd_lib_mod

  implicit none

! Local variables:

  integer               :: i
  integer               :: node
  double precision      :: total_time
  double precision      :: recip_time
  logical, save         :: initial_img_div_logged = .false.
  logical               :: write_log

  ! We place the following data in static storage to keep the Myrinet MPI
  ! implementation from going ape over mpi buffers on the stack (causes
  ! problems between mpi and ifc, resulting in huge stalls in processing).

  integer, save :: my_node_dat(4)

  ! Gather times from the other processes to the master process.  The
  ! reference to gbl_loadbal_node_dat is only used in the master process.
  ! The dat collected from each node are the "direct force time" due to image
  ! nonbonded calcs, the "reciprocal force time" due to pme nonbonded
  ! reciprocal force calcs, the "atom owner time" due to computations required
  ! of atom owners, and the "atom user time" due to managing used atom and
  ! image data structures.  The reciprocal force component is artificially
  ! inflated by a factor of nrespa to keep pathological loadbalancing from
  ! occurring if respa is in effect.  For non-respa runs, this has no effect,
  ! as nrespa is 1.  For respa runs where all tasks have a slice of the
  ! reciprocal force workload, this has no effect.  For respa runs where some
  ! tasks don't have a reciprocal force workload, this keeps the reciprocal
  ! force tasks from being given a load that they can't handle when they are
  !actually doing reciprocal force calcs.
 
  if (i_do_frc_ene) then
    my_node_dat(1) = elapsed_100usec_dirfrc
    my_node_dat(2) = elapsed_100usec_recipfrc * nrespa
    my_node_dat(3) = elapsed_100usec_atmowner
    my_node_dat(4) = elapsed_100usec_atmuser
  else
    my_node_dat(1:4) = 0
  end if

  call mpi_gather(my_node_dat, 4, mpi_integer, gbl_loadbal_node_dat, 4, &
                  mpi_integer, 0, pmemd_comm, err_code_mpi)

  if (master) then

10 format(/, a, i7, a /)
20 format('  ', a, i4)
30 format(/, '  ', a)
40 format('    ', 8i9)

    if (loadbal_verbose .gt. 2) then

      write(logfile, 10) &
        'Image Distribution Calc Data at run step ', loadbal_step_ctr, ':'
      write(logfile, 20) 'Pairlist builds done between calls: ', &
                         log_listbuild_call_ctr
      write(logfile, 30) &
        'Direct force time (in units of 100 usec) for each task:'
      write(logfile, 40) (gbl_loadbal_node_dat(1, i), i = 0, numtasks - 1)
      write(logfile, 30) &
        'Reciprocal force time (in units of 100 usec) for each task:'
      write(logfile, 40) (gbl_loadbal_node_dat(2, i), i = 0, numtasks - 1)
      write(logfile, 30) &
        'Atom ownership associated time (in units of 100 usec) for each task:'
      write(logfile, 40) (gbl_loadbal_node_dat(3, i), i = 0, numtasks - 1)
      write(logfile, 30) &
        'Atom usage associated time (in units of 100 usec) for each task:'
      write(logfile, 40) (gbl_loadbal_node_dat(4, i), i = 0, numtasks - 1)
      write(logfile, 30) &
        'Total force time (in units of 100 usec) for each task:'
      write(logfile, 40) &
       (gbl_loadbal_node_dat(1, i) + &
        gbl_loadbal_node_dat(2, i) + &
        gbl_loadbal_node_dat(3, i) + &
        gbl_loadbal_node_dat(4, i), i = 0, numtasks - 1)

    end if

    if (excl_recip .eq. 0) then
      call calc_new_img_distribution(gbl_loadbal_node_dat, total_time)
    else
      call calc_new_img_distribution_excl_recip(gbl_loadbal_node_dat, &
                                                total_time)
    end if

    call mpi_bcast(gbl_owned_imgs_tbl, size(gbl_owned_imgs_tbl), mpi_integer, &
                   0, pmemd_comm, err_code_mpi)

    write_log = .false.

    if (.not. initial_img_div_logged) then
      if (img_redist_ctr .eq. &
          img_redist_trigger_max / img_redist_trigger_divisor) then
        if (fft_redist_enabled) then
          write_log = .true.
          initial_img_div_logged = .true.
        end if
      end if
    end if
    if (loadbal_verbose .gt. 0) then
      write_log = .true.
    end if

    if (write_log) then
      write(logfile, 110) &
        'Image Distribution at run step ', loadbal_step_ctr, ':'
      write(logfile, 120) 'Count of images assigned to each task:'
      write(logfile, 130) &
        (gbl_owned_imgs_tbl(2, i), i = 0, numtasks - 1)
    end if

110 format(/, a, i7, a /)
120 format('  ', a)
130 format('    ', 8i9)

  else

    call mpi_bcast(gbl_owned_imgs_tbl, size(gbl_owned_imgs_tbl), mpi_integer, &
                   0, pmemd_comm, err_code_mpi)

  end if

  my_img_lo = gbl_owned_imgs_tbl(1, mytaskid) + 1
  my_img_hi = my_img_lo + gbl_owned_imgs_tbl(2, mytaskid) - 1

  log_owned_img_cnt = dble(my_img_hi - my_img_lo + 1)

  elapsed_100usec_dirfrc = 0
  elapsed_100usec_recipfrc = 0
  elapsed_100usec_atmowner = 0
  elapsed_100usec_atmuser = 0
  log_listbuild_call_ctr = 0
  last_loadbal_step_ctr = loadbal_step_ctr

  return

end subroutine do_img_redistribution

!*******************************************************************************
!
! Subroutine:  calc_new_img_distribution
!
! Description:  <TBS>
!              
!*******************************************************************************

subroutine calc_new_img_distribution(node_dat, total_time)

  use img_mod
  use parallel_mod
  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer               :: node_dat(4, 0:numtasks - 1) ! used by master
  double precision      :: total_time

! Local variables:

  integer               :: img_cnt
  integer               :: node
  integer               :: images_per_node(0:numtasks - 1) ! final result
  double precision      :: avg_time_per_img
  double precision      :: factor
  double precision      :: total_images
  double precision      :: target_time
  double precision      :: target_node_adj_time
  double precision      :: dirfrc_avg_time
  double precision      :: dirfrc_time_per_image
  double precision      :: dirfrc_total_time
  integer               :: max_total_time
  integer               :: min_total_time
  integer               :: dirfrc_limit         ! below this, use different
                                                ! algorithm for assigning
                                                ! images to node
  integer               :: images_limit         ! ditto...
  double precision      :: max_time_variation
  double precision      :: old_images_per_node(0:numtasks - 1)
  double precision      :: new_images_per_node(0:numtasks - 1)
  integer               :: total_times(0:numtasks - 1)

  ! Get total time in last run and average (target) time you want each 
  ! node to use.

  total_times(:) = node_dat(1, :) + &
                   node_dat(2, :) + &
                   node_dat(3, :) + &
                   node_dat(4, :)

  max_total_time = maxval(total_times)
  min_total_time = minval(total_times)

  total_time = 0.d0

  do node = 0, numtasks - 1
    total_time = total_time + dble(total_times(node))
  end do

  ! target_time is passed back for use in the caller as an average task
  ! time (for recip workload balancing)

  target_time = total_time / dble(frc_ene_task_cnt)

  max_time_variation = dble(max_total_time - min_total_time) * &
                       dble(frc_ene_task_cnt) / total_time

  ! We only readjust the image division if it is worthwhile.

! write(0,*)'| DBG: max img redist time variation = ', max_time_variation

  if (max_time_variation .gt. img_redist_time_tol) then

    ! Load balancing is handled differently on nodes with less than 5% the
    ! average workload because the standard algorithm does not work well at
    ! very low times / image counts (to say nothing of the possibility of
    ! dividing by 0...)

    dirfrc_total_time = 0.d0
    do node = 0, numtasks - 1
      dirfrc_total_time = dirfrc_total_time + dble(node_dat(1, node))
    end do

    dirfrc_avg_time = dirfrc_total_time / dble(frc_ene_task_cnt)
    dirfrc_limit = ceiling(0.05d0 * dirfrc_avg_time)
    dirfrc_time_per_image = dirfrc_total_time / dble(natom)

    images_limit = ceiling(0.05d0 * dble(natom) / dble(frc_ene_task_cnt))

    ! Fill in old images per node array:

    do node = 0, numtasks - 1
      old_images_per_node(node) = dble(gbl_owned_imgs_tbl(2, node))
    end do

    ! Calculate new image distributions. We correct for the nonadjustable 
    ! (everything except direct) times, make an estimate of processing time
    ! required per image based on the last cycle, and interpolated between old
    ! and new calculated values in order to damp the fluctuations, unless we
    ! are below limits for either the direct force time or image count for the
    ! node.  In that case, we use the average time required per image and
    ! the available time for direct force processing to estimate the correct
    ! number of images.

    total_images = 0.d0

    do node = 0, numtasks - 1

      if (is_frc_ene_task(node)) then
        target_node_adj_time = target_time - &
                               dble(node_dat(2, node) + &
                                    node_dat(3, node) + &
                                    node_dat(4, node))

        if (old_images_per_node(node) .ge. images_limit .and. &
            node_dat(1, node) .ge. dirfrc_limit) then
          factor = target_node_adj_time / &
                   dble(node_dat(1, node))
          if (factor .gt. 0.d0) then
            new_images_per_node(node) = factor * old_images_per_node(node)
          else
            new_images_per_node(node) = 0.d0
          end if
        else
          if (target_node_adj_time .gt. 0.d0) then
            new_images_per_node(node) = &
              target_node_adj_time / dirfrc_time_per_image
          else
            new_images_per_node(node) = 0.d0
          end if
        end if

        new_images_per_node(node) = (new_images_per_node(node) + &
                                     old_images_per_node(node)) / 2.d0
        total_images = total_images + new_images_per_node(node)
      else
        new_images_per_node(node) = 0.d0
      end if
    end do

    ! Now basically scale the results to get the correct image total as an
    ! integer.  This can still be incorrect after all this scaling and
    ! adjustment, so we have to make a final adjustment pass on the integers,
    ! making sure the total is correct and none are 0.

    factor = dble(natom) / dble(total_images)
  
    img_cnt = 0

    do node = 0, numtasks - 1
      if (is_frc_ene_task(node)) then
        images_per_node(node) = nint(new_images_per_node(node) * factor)
        img_cnt = img_cnt + images_per_node(node)
      else
        images_per_node(node) = 0
      end if
    end do

    if (img_cnt .lt. natom) then

!     write(0,*)'DBG: Img redist - correcting low img_cnt of ', img_cnt
      node = 0

      do while (img_cnt .lt. natom)

        if (is_frc_ene_task(node)) then
          images_per_node(node) = images_per_node(node) + 1
          img_cnt = img_cnt + 1
        end if

        node = node + 1

        if (node .ge. numtasks) node = 0

      end do

    else if (img_cnt .gt. natom) then

!     write(0,*)'DBG: Img redist - correcting high img_cnt of ', img_cnt
      node = 0

      do while (img_cnt .gt. natom)

        if (images_per_node(node) .gt. 1) then
          images_per_node(node) = images_per_node(node) - 1
          img_cnt = img_cnt - 1
        end if

        node = node + 1

        if (node .ge. numtasks) node = 0

      end do

    end if

    gbl_owned_imgs_tbl(2, 0:numtasks-1) = images_per_node(0:numtasks-1)

    img_cnt = 0

    do node = 0, numtasks - 1
      gbl_owned_imgs_tbl(1, node) = img_cnt
      img_cnt = img_cnt + gbl_owned_imgs_tbl(2, node)
    end do

  end if

  return

end subroutine calc_new_img_distribution

!*******************************************************************************
!
! Subroutine:  calc_new_img_distribution_excl_recip
!
! Description:  <TBS>
!              
!*******************************************************************************

subroutine calc_new_img_distribution_excl_recip(node_dat, total_time)

  use img_mod
  use mdin_ewald_dat_mod
  use parallel_mod
  use parallel_dat_mod
  use pme_recip_dat_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer               :: node_dat(4, 0:numtasks - 1) ! used by master
  double precision      :: total_time

! Local variables:

  integer               :: img_cnt
  integer               :: node
  integer               :: images_per_node(0:numtasks - 1) ! final result
  double precision      :: avg_time_per_img
  double precision      :: factor
  double precision      :: total_images
  double precision      :: target_time
  double precision      :: target_node_adj_time
  double precision      :: dirfrc_avg_time
  double precision      :: dirfrc_time_per_image
  double precision      :: dirfrc_total_time
  integer               :: max_total_time
  integer               :: min_total_time
  integer               :: dirfrc_limit         ! below this, use different
                                                ! algorithm for assigning
                                                ! images to node
  integer               :: images_limit         ! ditto...
  integer               :: dir_numtasks
  logical               :: saved_master_recip_stat
  double precision      :: max_time_variation
  double precision      :: old_images_per_node(0:numtasks - 1)
  double precision      :: new_images_per_node(0:numtasks - 1)
  integer               :: total_times(0:numtasks - 1)
  logical               :: is_direct_task(0:numtasks - 1)

  ! Here we just load balance direct force computation, allowing time for
  ! bond-angle-dihedral computation.  The reciprocal tasks are out of the
  ! picture in this scenario.

  saved_master_recip_stat = is_recip_task(0)

  dir_numtasks = frc_ene_task_cnt - recip_numtasks

  if (excl_master .ne. 0) then
    is_recip_task(0) = .true.           ! done to exclude master...
    dir_numtasks = dir_numtasks - 1
  end if

  is_direct_task(:) = .not. is_recip_task(:)

  do node = 0, numtasks - 1
    if (is_direct_task(node)) then
      total_times(node) = node_dat(1, node) + &
                          node_dat(3, node) + &
                          node_dat(4, node)
    else
      total_times(node) = 0
    end if
  end do

  max_total_time = maxval(total_times)
  min_total_time = minval(total_times)

  total_time = 0.d0

  do node = 0, numtasks - 1
    total_time = total_time + dble(total_times(node))
  end do

  ! target_time is passed back for use in the caller as an average task
  ! time (for recip workload balancing)

  target_time = total_time / dble(dir_numtasks)

  max_time_variation = dble(max_total_time - min_total_time) * &
                       dble(dir_numtasks) / total_time

  ! We only readjust the image division if it is worthwhile.

! write(0,*)'| DBG: max img redist time variation = ', max_time_variation

  if (max_time_variation .gt. img_redist_time_tol) then

    ! Load balancing is handled differently on nodes with less than 5% the
    ! average workload because the standard algorithm does not work well at
    ! very low times / image counts (to say nothing of the possibility of
    ! dividing by 0...)

    dirfrc_total_time = 0.d0
    do node = 0, numtasks - 1
      dirfrc_total_time = dirfrc_total_time + dble(node_dat(1, node))
    end do

    dirfrc_avg_time = dirfrc_total_time / dble(dir_numtasks)
    dirfrc_limit = ceiling(0.05d0 * dirfrc_avg_time)
    dirfrc_time_per_image = dirfrc_total_time / dble(natom)

    images_limit = ceiling(0.05d0 * dble(natom) / dble(dir_numtasks))

    ! Fill in old images per node array:

    do node = 0, numtasks - 1
      old_images_per_node(node) = dble(gbl_owned_imgs_tbl(2, node))
    end do

    ! Calculate new image distributions. We correct for the nonadjustable 
    ! (recip and bnd-angle-dihed) times, make an estimate of processing time
    ! required per image based on the last cycle, and interpolated between old
    ! and new calculated values in order to damp the fluctuations, unless we
    ! are below limits for either the direct force time or image count for the
    ! node.  In that case, we use the average time required per image and
    ! the available time for direct force processing to estimate the correct
    ! number of images.

    total_images = 0.d0

    do node = 0, numtasks - 1

      if (is_direct_task(node)) then
        target_node_adj_time = target_time - dble(node_dat(3, node) + &
                                                  node_dat(4, node))

        if (old_images_per_node(node) .ge. images_limit .and. &
            node_dat(1, node) .ge. dirfrc_limit) then
          factor = target_node_adj_time / dble(node_dat(1, node))
          if (factor .gt. 0.d0) then
            new_images_per_node(node) = factor * old_images_per_node(node)
          else
            new_images_per_node(node) = 0.d0
          end if
        else
          if (target_node_adj_time .gt. 0.d0) then
            new_images_per_node(node) = &
              target_node_adj_time / dirfrc_time_per_image
          else
            new_images_per_node(node) = 0.d0
          end if
        end if

        new_images_per_node(node) = (new_images_per_node(node) + &
                                     old_images_per_node(node)) / 2.d0
        total_images = total_images + new_images_per_node(node)
      else
        new_images_per_node(node) = 0.d0
      end if
    end do

    ! Now basically scale the results to get the correct image total as an
    ! integer.  This can still be incorrect after all this scaling and
    ! adjustment, so we have to make a final adjustment pass on the integers,
    ! making sure the total is correct and none are 0.

    factor = dble(natom) / dble(total_images)
  
    img_cnt = 0

    do node = 0, numtasks - 1
      if (is_direct_task(node)) then
        images_per_node(node) = nint(new_images_per_node(node) * factor)
        img_cnt = img_cnt + images_per_node(node)
      else
        images_per_node(node) = 0
      end if
    end do

    if (img_cnt .lt. natom) then

!     write(0,*)'DBG: Img redist - correcting low img_cnt of ', img_cnt
      node = 0

      do while (img_cnt .lt. natom)

        if (is_direct_task(node)) then
          images_per_node(node) = images_per_node(node) + 1
          img_cnt = img_cnt + 1
        end if

        node = node + 1

        if (node .ge. numtasks) node = 0

      end do

    else if (img_cnt .gt. natom) then

!     write(0,*)'DBG: Img redist - correcting high img_cnt of ', img_cnt
      node = 0

      do while (img_cnt .gt. natom)

        if (images_per_node(node) .gt. 1) then
          images_per_node(node) = images_per_node(node) - 1
          img_cnt = img_cnt - 1
        end if

        node = node + 1

        if (node .ge. numtasks) node = 0

      end do

    end if

    gbl_owned_imgs_tbl(2, 0:numtasks-1) = images_per_node(0:numtasks-1)

    img_cnt = 0

    do node = 0, numtasks - 1
      gbl_owned_imgs_tbl(1, node) = img_cnt
      img_cnt = img_cnt + gbl_owned_imgs_tbl(2, node)
    end do

  end if

  is_recip_task(0) = saved_master_recip_stat

  return

end subroutine calc_new_img_distribution_excl_recip
#endif /* MPI */

end module loadbal_mod
