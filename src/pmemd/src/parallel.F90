#include "copyright.i"
!precision setting for SPDP, DPDP, SPSP
#include "include_precision.i"


!*******************************************************************************
!
! Module: parallel_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module parallel_mod

  use gbl_datatypes_mod
  use parallel_dat_mod
#ifdef TIME_TEST
  use timers_mod
#endif

  implicit none


  private       ! Everything below here is private unless specified public

  ! Public variables:

  logical, save, public :: i_do_frc_ene = .true.

  integer, save, public :: frc_ene_task_cnt = 0

  ! Under high scaling, the master is not assigned frc-ene work but is
  ! also not listed as a "comm task".  The master has a unique role, effectively
  ! as a master of the communications tasks.
  
  logical, allocatable, save, public    :: is_frc_ene_task(:)
  integer, allocatable, save, public    :: gbl_frc_ene_task_lst(:)

#if defined(MPI)

  ! Public subroutines:

  public        :: parallel_setup, &
                   setup_frc_ene_task_dat, &
                   get_send_atm_lst, &
                   get_send_ips_atm_lst, &
                   save_used_atom_crds, &
                   get_img_frc_distribution, &
                   get_img_frc_ips_distribution, &
                   distribute_img_frcs, &
                   distribute_img_ips_frcs, &
                   distribute_img_ti_frcs, &
                   distribute_crds, &
                   zero_extra_used_atm_img_frcs, &
                   zero_extra_used_atm_ti_img_frcs, &
                   mpi_allgathervec, &
                   do_atm_redistribution, &
                   mpi_gathervec

  ! Private data:

  ! PRIVATE DATA USED FOR THE LOWSCALING IMPLEMENTATION:

  integer, allocatable, save    :: pvt_taskmap(:)
  integer, allocatable, save    :: pvt_inv_taskmap(:)
  integer, allocatable, save    :: pvt_send_atm_lst(:)
  integer, allocatable, save    :: pvt_send_atm_cnts(:)
  integer, allocatable, save    :: pvt_recv_atm_cnts(:)

  integer, allocatable, save    :: pvt_send_atm_ips_lst(:)
  integer, allocatable, save    :: pvt_send_atm_ips_cnts(:)
  integer, allocatable, save    :: pvt_recv_atm_ips_cnts(:)

  ! This is allocated/reallocated only in do_atm_distribution() because it
  ! can need to change with atom workload reassignment.

  integer, allocatable, save    :: pvt_recv_atm_lsts(:,:)
  integer, allocatable, save    :: pvt_recv_atm_ips_lsts(:,:)

  ! PRIVATE DATA:

  integer, allocatable, save    :: pvt_owned_atm_cnts(:)
  integer, allocatable, save    :: pvt_atm_offsets(:)
  integer, allocatable, save    :: pvt_vec_offsets(:)
  integer, allocatable, save    :: pvt_vec_rcvcnts(:)

  integer, save                 :: extra_used_atm_cnt = 0
  integer, allocatable, save    :: pvt_extra_used_atms(:)

contains

!*******************************************************************************
!
! Subroutine:  setup_frc_ene_task_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine setup_frc_ene_task_dat(num_ints, num_reals)

  use file_io_dat_mod
  use gbl_constants_mod
  use parallel_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:

! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed
  integer               :: lst_ctr
  integer               :: taskid

  i_do_frc_ene = .true.
  frc_ene_task_cnt = numtasks
  allocate(is_frc_ene_task(0:numtasks - 1), &
           gbl_frc_ene_task_lst(numtasks), &
           stat = alloc_failed)
  if (alloc_failed .ne. 0) call setup_alloc_error
  num_ints = num_ints + size(is_frc_ene_task) + &
                        size(gbl_frc_ene_task_lst)
  is_frc_ene_task(:) = .true.

  lst_ctr = 0
  do taskid = 0, numtasks - 1
    if (is_frc_ene_task(taskid)) then
      lst_ctr = lst_ctr + 1
      gbl_frc_ene_task_lst(lst_ctr) = taskid
    end if
  end do

  return

end subroutine setup_frc_ene_task_dat

!*******************************************************************************
!
! Subroutine:  parallel_setup
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine parallel_setup(num_ints, num_reals)

  use gbl_constants_mod
  use img_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use pme_recip_dat_mod
  use pme_blk_fft_mod
  use pme_slab_fft_mod
  use pmemd_lib_mod
  use prfs_mod
  use prmtop_dat_mod
  use remd_mod, only : remd_method

  implicit none

! Formal arguments:

! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed
  integer               :: i
  integer               :: node
  double precision      :: recip_workload_estimate
  integer               :: taskmap_idx
  logical               :: write_log

! Used for logging workload distribution statistics:

  log_recip_nstep = 1.d0 / dble(nstlim)

! The pvt_taskmap is an ordering of tasks other than this task ordered from
! mytaskid + 1 to mytaskid - 1, with wraparound at numtasks.  Also
! Allocate storage for various structures associated with keeping track of the
! division of the atom workload.  Some of these data structures are needed
! prior to the initial atom division call.

  if (ips .gt. 0) then
    allocate(pvt_taskmap(numtasks - 1), &
             pvt_inv_taskmap(numtasks - 1), &
             pvt_owned_atm_cnts(0:numtasks), &
             pvt_atm_offsets(0:numtasks), &
             pvt_vec_offsets(0:numtasks), &
             pvt_vec_rcvcnts(0:numtasks), &
             pvt_send_atm_lst(natom), &
             pvt_send_atm_cnts(0:numtasks - 1), &
             pvt_recv_atm_cnts(0:numtasks - 1), &
             pvt_send_atm_ips_lst(natom), &
             pvt_send_atm_ips_cnts(0:numtasks - 1), &
             pvt_recv_atm_ips_cnts(0:numtasks - 1), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    num_ints = num_ints + size(pvt_taskmap) + &
                          size(pvt_inv_taskmap) + &
                          size(pvt_owned_atm_cnts) + &
                          size(pvt_atm_offsets) + &
                          size(pvt_vec_offsets) + &
                          size(pvt_vec_rcvcnts) + &
                          size(pvt_send_atm_lst) + &
                          size(pvt_send_atm_cnts) + &
                          size(pvt_recv_atm_cnts) + &
                          size(pvt_send_atm_ips_lst) + &
                          size(pvt_send_atm_ips_cnts) + &
                          size(pvt_recv_atm_ips_cnts)
  else
    allocate(pvt_taskmap(numtasks - 1), &
             pvt_inv_taskmap(numtasks - 1), &
             pvt_owned_atm_cnts(0:numtasks), &
             pvt_atm_offsets(0:numtasks), &
             pvt_vec_offsets(0:numtasks), &
             pvt_vec_rcvcnts(0:numtasks), &
             pvt_send_atm_lst(natom), &
             pvt_send_atm_cnts(0:numtasks - 1), &
             pvt_recv_atm_cnts(0:numtasks - 1), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    num_ints = num_ints + size(pvt_taskmap) + &
                          size(pvt_inv_taskmap) + &
                          size(pvt_owned_atm_cnts) + &
                          size(pvt_atm_offsets) + &
                          size(pvt_vec_offsets) + &
                          size(pvt_vec_rcvcnts) + &
                          size(pvt_send_atm_lst) + &
                          size(pvt_send_atm_cnts) + &
                          size(pvt_recv_atm_cnts) 
  endif ! (ips .gt. 0)

  ! Set up data structures used in n to n mpi exchanges.  The taskmaps set up
  ! the order for what is essentially a "synchronous shuffle exchange"
  ! (Tam and Wang) which we discovered independently.  We use it for our
  ! asynchronous comm.

  taskmap_idx = 0
  node = mytaskid + 1

  do
    if (node .ge. numtasks) node = 0
    if (node .eq. mytaskid) exit
    taskmap_idx = taskmap_idx + 1
    pvt_taskmap(taskmap_idx) = node
    node = node + 1
  end do

  ! The order of pvt_inv_taskmap is inverted.  If you receive and send
  ! sequentially, it works best to recv via pvt_taskmap while you
  ! simultaneously send via pvt_inv_taskmap (or vice versa).

  do taskmap_idx = 1, numtasks - 1
    pvt_inv_taskmap(taskmap_idx) = pvt_taskmap(numtasks - taskmap_idx)
  end do

  ! The following checks for minimum atoms, prf's, and molecules per
  ! processor are intended to avoid hitting conditions we may not
  ! have considered.  The use of more processors than atoms, prf's or
  ! molecules would be a bit ridiculous anyway...  These limits are arbitrary,
  ! and subject to change if we determine that problems don't occur.

  ! Check for sufficient atoms for the number of processors:
  ! (only for non-REMD)
  if (remd_method.le.0) then

  if (natom .lt. 10 * numtasks) then
    write(mdout, '(a,a)') error_hdr, 'Must have 10x more atoms than processors!'
    call mexit(6, 1)
  end if

  ! Check for sufficient prf's for the number of processors:

  if (gbl_prf_cnt .lt. 4 * numtasks) then
    write(mdout, '(a,a)') error_hdr, &
      'Must have 4x more pseudo-residue fragments than processors!'
    call mexit(6, 1)
  end if
  end if

  ! We no longer require a minimum number of molecules under NTP, since
  ! molecule fragmentation should take care of any issues...

  ! Set up image division for force calcs.  Images ARE assigned in contiguous
  ! blocks, without wraparound.  For respa runs and minimizations, asymmetric
  ! fft slab load balancing (basically, trying to use as few tasks doing
  ! fft's/recip force) will not be done for various technical reasons.

  if (block_fft .eq. 0) then
    recip_workload_estimate = slab_fft_workload_estimate
  else
    recip_workload_estimate = blk_fft_workload_estimate
  end if

if (.not. usemidpoint) then
  if (excl_recip .eq. 0 .and. .not. usemidpoint) then
    if (nrespa .eq. 1 .and. imin .eq. 0) then
      call divide_images_recip_biased(natom, recip_workload_estimate, &
                                      gbl_owned_imgs_tbl, is_recip_task, &
                                      recip_numtasks, &
                                      frc_ene_task_cnt, gbl_frc_ene_task_lst, &
                                      .false., 0)

    else
      call divide_images_evenly(natom, gbl_owned_imgs_tbl, &
                                frc_ene_task_cnt, gbl_frc_ene_task_lst, &
                                .true., 0)
    end if
  else
      call divide_images_excl_recip(natom, gbl_owned_imgs_tbl, &
                                    frc_ene_task_cnt, gbl_frc_ene_task_lst, &
                                    is_recip_task, recip_numtasks, &
                                    .true., 0)
  end if

  ! All nodes divide the work in runmd. It is a residue-based division
  ! for constant volume simulations and a molecule-based division for
  ! constant pressure simulations.  The division is map- and list-based, in
  ! order to improve locality.  The division will be redone periodically to
  ! maintain reasonable locality as atoms move.
  
  write_log = .not. (nrespa .eq. 1 .and. imin .eq. 0)

  call do_atm_distribution(natom, num_ints, num_reals, write_log, 0)
end if ! not usemidpoint
  return

end subroutine parallel_setup

!*******************************************************************************
!
! Subroutine:  get_send_atm_lst
!
! Description: <TBS>
!
!*******************************************************************************

subroutine get_send_atm_lst(atm_cnt)

  use img_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer, intent(in)   :: atm_cnt

! Local variables

  integer               :: i
  integer               :: send_atm_cnt
  integer               :: send_atm_cnts_sum
  integer               :: senders_cnt

#ifdef COMM_TIME_TEST
  call start_test_timer(13, 'get_send_atm_lst', 0)
#endif

  call pvt_get_send_atm_lst(atm_cnt, gbl_img_atm_map, gbl_used_img_map, &
                            gbl_used_img_cnt, gbl_used_img_lst, &
                            gbl_atm_img_map, gbl_atm_owner_map, &
                            pvt_atm_offsets, &
                            pvt_send_atm_lst, pvt_send_atm_cnts)

  ! Sum up the atoms for which this process sends data elsewhere.  This is a
  ! measure of locality sent to the master and used in load balancing.  This
  ! will be adjusted to an average between load balancings.  The loadbalancing
  ! code controls when this info is actually used.

  send_atm_cnts_sum = 0
  senders_cnt = 0       ! used for logging...

  do i = 0, numtasks - 1
    send_atm_cnt = pvt_send_atm_cnts(i)
    if (send_atm_cnt .gt. 0) then
      send_atm_cnts_sum = send_atm_cnts_sum + send_atm_cnt
      senders_cnt = senders_cnt + 1
    end if
  end do

  if (pvt_send_atm_cnts(mytaskid) .ne. 0) then
    send_atm_cnts_sum = send_atm_cnts_sum - pvt_send_atm_cnts(mytaskid)
    senders_cnt = senders_cnt - 1
  end if

  log_used_atm_cnt = dble(send_atm_cnts_sum)
  log_used_atm_source_cnt = dble(senders_cnt)
  log_listbuild_call_ctr = log_listbuild_call_ctr + 1

! write(0,*)'DBG: task, pvt_send_atm_cnts()=', mytaskid, pvt_send_atm_cnts(:)

#ifdef COMM_TIME_TEST
  call stop_test_timer(13)
#endif
  return

end subroutine get_send_atm_lst

!*******************************************************************************
!
! Subroutine:  pvt_get_send_atm_lst
!
! Description: <TBS>
!
!*******************************************************************************

subroutine pvt_get_send_atm_lst(atm_cnt, img_atm_map, used_img_map, &
                                used_img_cnt, used_img_lst, &
                                atm_img_map, atm_owner_map, off_tbl, &
                                send_atm_lst, send_atm_cnts)

  use gbl_datatypes_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  integer               :: img_atm_map(atm_cnt)
  integer(byte)         :: used_img_map(atm_cnt)
  integer               :: used_img_cnt
  integer               :: used_img_lst(*)
  integer               :: atm_img_map(atm_cnt)
  integer               :: atm_owner_map(atm_cnt)
  integer               :: off_tbl(0:numtasks)
  integer               :: send_atm_lst(atm_cnt)
  integer               :: send_atm_cnts(0 : numtasks - 1)

! Local variables

  integer               :: atm_id, img_id, node
  integer               :: lst_idx
  integer               :: i

  send_atm_cnts(0: numtasks - 1) = 0

  do lst_idx = 1, used_img_cnt
    img_id = used_img_lst(lst_idx)
    atm_id = img_atm_map(img_id)
    node = atm_owner_map(atm_id)
    send_atm_cnts(node) = send_atm_cnts(node) + 1
    send_atm_lst(off_tbl(node) + send_atm_cnts(node)) = atm_id
  end do

  ! Okay, now add the extra used atoms, which are atoms we don't own but
  ! that we have to keep coordinates updated on because they are used in
  ! bond, angle, and dihedral calculations.  A little gotcha...
  ! There are very few of these (like 5-10 at most typically), but they
  ! really screw things up if you don't update their coordinates!
  ! We only need the coordinates, but because the send atoms mechanism is
  ! shared, forces will also get propagated, so it is important to zero the
  ! forces if they are not in use.

  do i = 1, extra_used_atm_cnt
    atm_id = pvt_extra_used_atms(i)
    if (used_img_map(atm_img_map(atm_id)) .eq. 0) then
      ! It is not already marked as used, so we include it in send atms list.
      node = atm_owner_map(atm_id)
      send_atm_cnts(node) = send_atm_cnts(node) + 1
      send_atm_lst(off_tbl(node) + send_atm_cnts(node)) = atm_id
    end if
  end do

  return

end subroutine pvt_get_send_atm_lst

!*******************************************************************************
!
! Subroutine:  get_send_ips_atm_lst
!
! Description: <TBS>
!
!*******************************************************************************

subroutine get_send_ips_atm_lst(atm_cnt)

  use img_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer, intent(in)   :: atm_cnt

! Local variables

  integer               :: i
  integer               :: send_atm_cnt
  integer               :: send_atm_cnts_sum
  integer               :: senders_cnt

#ifdef COMM_TIME_TEST
  call start_test_timer(13, 'get_send_atm_lst', 0)
#endif


  call pvt_get_send_ips_atm_lst(atm_cnt, gbl_img_atm_map, gbl_used_img_ips_map, &
                            gbl_used_img_ips_cnt, gbl_used_img_ips_lst, &
                            gbl_atm_img_map, gbl_atm_owner_map, &
                            pvt_atm_offsets, &
                            pvt_send_atm_ips_lst, pvt_send_atm_ips_cnts)

  ! Sum up the atoms for which this process sends data elsewhere.  This is a
  ! measure of locality sent to the master and used in load balancing.  This
  ! will be adjusted to an average between load balancings.  The loadbalancing
  ! code controls when this info is actually used.

  send_atm_cnts_sum = 0
  senders_cnt = 0       ! used for logging...

  do i = 0, numtasks - 1
    send_atm_cnt = pvt_send_atm_ips_cnts(i)
    if (send_atm_cnt .gt. 0) then
      send_atm_cnts_sum = send_atm_cnts_sum + send_atm_cnt
      senders_cnt = senders_cnt + 1
    end if
  end do

  if (pvt_send_atm_ips_cnts(mytaskid) .ne. 0) then
    send_atm_cnts_sum = send_atm_cnts_sum - pvt_send_atm_ips_cnts(mytaskid)
    senders_cnt = senders_cnt - 1
  end if

  log_used_atm_cnt = dble(send_atm_cnts_sum)
  log_used_atm_source_cnt = dble(senders_cnt)
  log_listbuild_call_ctr = log_listbuild_call_ctr + 1



! write(0,*)'DBG: task, pvt_send_atm_ips_cnts()=', mytaskid, pvt_send_atm_ips_cnts(:)

#ifdef COMM_TIME_TEST
  call stop_test_timer(13)
#endif
  return

end subroutine get_send_ips_atm_lst

!*******************************************************************************
!
! Subroutine:  pvt_get_send_ips_atm_lst
!
! Description: <TBS>
!
!*******************************************************************************

subroutine pvt_get_send_ips_atm_lst(atm_cnt, img_atm_map, used_img_ips_map, &
                                used_img_ips_cnt, used_img_ips_lst, &
                                atm_img_map, atm_owner_map, off_tbl, &
                                send_atm_lst, send_atm_cnts)

  use gbl_datatypes_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  integer               :: img_atm_map(atm_cnt)
  integer(byte)         :: used_img_ips_map(atm_cnt)
  integer               :: used_img_ips_cnt
  integer               :: used_img_ips_lst(*)
  integer               :: atm_img_map(atm_cnt)
  integer               :: atm_owner_map(atm_cnt)
  integer               :: off_tbl(0:numtasks)
  integer               :: send_atm_lst(atm_cnt)
  integer               :: send_atm_cnts(0 : numtasks - 1)

! Local variables

  integer               :: atm_id, img_id, node
  integer               :: lst_idx
  integer               :: i

  send_atm_cnts(0: numtasks - 1) = 0

  do lst_idx = 1, used_img_ips_cnt
    img_id = used_img_ips_lst(lst_idx)
    atm_id = img_atm_map(img_id)
    node = atm_owner_map(atm_id)
    send_atm_cnts(node) = send_atm_cnts(node) + 1
    send_atm_lst(off_tbl(node) + send_atm_cnts(node)) = atm_id
  end do

  ! Okay, now add the extra used atoms, which are atoms we don't own but
  ! that we have to keep coordinates updated on because they are used in
  ! bond, angle, and dihedral calculations.  A little gotcha...
  ! There are very few of these (like 5-10 at most typically), but they
  ! really screw things up if you don't update their coordinates!
  ! We only need the coordinates, but because the send atoms mechanism is
  ! shared, forces will also get propagated, so it is important to zero the
  ! forces if they are not in use.

  do i = 1, extra_used_atm_cnt
    atm_id = pvt_extra_used_atms(i)
    if (used_img_ips_map(atm_img_map(atm_id)) .eq. 0) then
      ! It is not already marked as used, so we include it in send atms list.
      node = atm_owner_map(atm_id)
      send_atm_cnts(node) = send_atm_cnts(node) + 1
      send_atm_lst(off_tbl(node) + send_atm_cnts(node)) = atm_id
    end if
  end do

  return

end subroutine pvt_get_send_ips_atm_lst

!*******************************************************************************
!
! Subroutine:  save_used_atom_crds
!
! Description:  This is needed for the skin test for buffered pairlists, but
!               is also used in adjusting image coordinates between list
!               builds, hence the need to look at the send_atm_lst.
!*******************************************************************************

subroutine save_used_atom_crds(atm_cnt, crd, saved_crd)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: crd(3, atm_cnt)
  double precision      :: saved_crd(3, atm_cnt)


  call pvt_save_used_atom_crds(atm_cnt, crd, saved_crd, pvt_send_atm_lst, &
                               pvt_send_atm_cnts, pvt_atm_offsets)

  return

end subroutine save_used_atom_crds

!*******************************************************************************
!
! Subroutine:  pvt_save_used_atom_crds
!
! Description:  <TBS>
!*******************************************************************************

subroutine pvt_save_used_atom_crds(atm_cnt, crd, saved_crd, send_atm_lst, &
                                   send_atm_cnts, off_tbl)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: crd(3, atm_cnt)
  double precision      :: saved_crd(3, atm_cnt)
  integer               :: send_atm_lst(atm_cnt)
  integer               :: send_atm_cnts(0: numtasks - 1)
  integer               :: off_tbl(0:numtasks)

! Local variables:

  integer               :: atm_lst_idx, atm_id
  integer               :: i
  integer               :: node, node_offset
  integer               :: send_cnt

    do atm_lst_idx = 1, my_atm_cnt
      atm_id = gbl_my_atm_lst(atm_lst_idx)
      saved_crd(1, atm_id) = crd(1, atm_id)
      saved_crd(2, atm_id) = crd(2, atm_id)
      saved_crd(3, atm_id) = crd(3, atm_id)
    end do

    do node = 0, numtasks - 1

      if (node .ne. mytaskid) then
        node_offset = off_tbl(node)
        send_cnt = send_atm_cnts(node)

        ! This also picks up stuff on the extra atoms list, which is probably
        ! not necessary but not harmful.

        do i = 1, send_cnt
          atm_id = send_atm_lst(node_offset + i)
          saved_crd(1, atm_id) = crd(1, atm_id)
          saved_crd(2, atm_id) = crd(2, atm_id)
          saved_crd(3, atm_id) = crd(3, atm_id)
        end do
      end if

    end do

  return

end subroutine pvt_save_used_atom_crds

!*******************************************************************************
!
! Subroutine:  get_img_frc_distribution
!
! Description: <TBS>
!
!*******************************************************************************

subroutine get_img_frc_distribution(atm_cnt)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer       :: atm_cnt

! Local variables:

  integer       :: i
  integer       :: recv_atm_cnt
  integer       :: recv_atm_cnts_sum
  integer       :: receivers_cnt

#ifdef COMM_TIME_TEST
  call start_test_timer(1, 'get_img_frc_distribution', 0)
#endif

  call pvt_get_img_frc_distribution(atm_cnt, pvt_atm_offsets, pvt_taskmap, &
                                    pvt_inv_taskmap, pvt_send_atm_lst, &
                                    pvt_send_atm_cnts, pvt_recv_atm_lsts, &
                                    pvt_recv_atm_cnts, pvt_owned_atm_cnts)
  
  recv_atm_cnts_sum = 0
  receivers_cnt = 0     ! used for logging...

  if (my_atm_cnt .gt. 0) then
    do i = 0, numtasks - 1
      recv_atm_cnt = pvt_recv_atm_cnts(i)
      if (recv_atm_cnt .gt. 0) then
        recv_atm_cnts_sum = recv_atm_cnts_sum + recv_atm_cnt
        receivers_cnt = receivers_cnt + 1
      end if
    end do
  end if

  log_provided_atm_cnt = dble(recv_atm_cnts_sum)
  log_provided_atm_sink_cnt = dble(receivers_cnt)

#ifdef COMM_TIME_TEST
  call stop_test_timer(1)
#endif

  return

end subroutine get_img_frc_distribution

! This is the default async implementation...
!*******************************************************************************
!
! Subroutine:  pvt_get_img_frc_distribution
!
! Description: <TBS>
!
!*******************************************************************************

subroutine pvt_get_img_frc_distribution(atm_cnt, off_tbl, recv_taskmap, &
                                        send_taskmap, send_atm_lst, &
                                        send_atm_cnts, recv_atm_lsts, &
                                        recv_atm_cnts, owned_atm_cnts)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  integer               :: off_tbl(0:numtasks)
  integer               :: recv_taskmap(numtasks - 1)
  integer               :: send_taskmap(numtasks - 1)
  integer               :: send_atm_lst(atm_cnt)
  integer               :: send_atm_cnts(0 : numtasks - 1)
  integer               :: recv_atm_lsts(my_atm_cnt, numtasks - 1)
  integer               :: recv_atm_cnts(0 : numtasks - 1)
  integer               :: owned_atm_cnts(0:numtasks)

! Local variables

  integer               :: wait_call
  integer               :: taskmap_idx
  integer               :: recv_task, send_task
  integer               :: recv_cnt, send_cnt, send_offset
  logical, save         :: recv_atm_cnts_zeroed = .false.
  integer               :: irecv_stat(mpi_status_size)
  integer               :: isend_stat(mpi_status_size, numtasks - 1)
  integer               :: recv_req(numtasks - 1)
  integer               :: send_req(numtasks - 1)

  ! No need to receive from yourself, so there are numtasks - 1 bufs:

  ! Exchange atom lists...

  ! Post the asynchronous receives first:

  if (my_atm_cnt .ne. 0) then

    do taskmap_idx = 1, numtasks - 1

      recv_task = recv_taskmap(taskmap_idx)

      call mpi_irecv(recv_atm_lsts(1, taskmap_idx), my_atm_cnt, mpi_integer, &
                     recv_task, gifd_tag, pmemd_comm, &
                     recv_req(taskmap_idx), err_code_mpi) 
    end do

  end if

  ! Now set up and post the asynchronous sends:

  do taskmap_idx = 1, numtasks - 1

    send_task = send_taskmap(taskmap_idx)

    if (owned_atm_cnts(send_task) .ne. 0) then

      send_offset = off_tbl(send_task) + 1
      send_cnt = send_atm_cnts(send_task)

      call mpi_isend(send_atm_lst(send_offset), send_cnt, &
                     mpi_integer, send_task, gifd_tag, pmemd_comm, &
                     send_req(taskmap_idx), err_code_mpi)

    else

      send_req(taskmap_idx) = MPI_REQUEST_NULL

    end if

  end do

  ! Wait on and process the pending receive requests:

  if (my_atm_cnt .ne. 0) then

    do wait_call = 1, numtasks - 1

      call mpi_waitany(numtasks - 1, recv_req, taskmap_idx, irecv_stat, &
                       err_code_mpi)

      recv_task = recv_taskmap(taskmap_idx)

      call mpi_get_count(irecv_stat, mpi_integer, recv_cnt, err_code_mpi)

      recv_atm_cnts(recv_task) = recv_cnt

    end do

    recv_atm_cnts_zeroed = .false.

  else

    if (.not. recv_atm_cnts_zeroed) then
      recv_atm_cnts(:) = 0
      recv_atm_cnts_zeroed = .true.
    end if

  end if

  ! Wait for all sends to complete:

  call mpi_waitall(numtasks - 1, send_req, isend_stat, err_code_mpi)

  recv_atm_cnts(mytaskid) = 0

  return

end subroutine pvt_get_img_frc_distribution

!*******************************************************************************
!
! Subroutine:  get_img_frc_ips_distribution
!
! Description: <TBS>
!
!*******************************************************************************

subroutine get_img_frc_ips_distribution(atm_cnt)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer       :: atm_cnt

! Local variables:

  integer       :: i
  integer       :: recv_atm_cnt
  integer       :: recv_atm_cnts_sum
  integer       :: receivers_cnt

#ifdef COMM_TIME_TEST
  call start_test_timer(1, 'get_img_frc_distribution', 0)
#endif

  call pvt_get_img_frc_ips_distribution(atm_cnt, pvt_atm_offsets, pvt_taskmap, &
                                    pvt_inv_taskmap, pvt_send_atm_ips_lst, &
                                    pvt_send_atm_ips_cnts, pvt_recv_atm_ips_lsts, &
                                    pvt_recv_atm_ips_cnts, pvt_owned_atm_cnts)
  
  recv_atm_cnts_sum = 0
  receivers_cnt = 0     ! used for logging...

  if (my_atm_cnt .gt. 0) then
    do i = 0, numtasks - 1
      recv_atm_cnt = pvt_recv_atm_ips_cnts(i)
      if (recv_atm_cnt .gt. 0) then
        recv_atm_cnts_sum = recv_atm_cnts_sum + recv_atm_cnt
        receivers_cnt = receivers_cnt + 1
      end if
    end do
  end if

  log_provided_atm_cnt = dble(recv_atm_cnts_sum)
  log_provided_atm_sink_cnt = dble(receivers_cnt)

#ifdef COMM_TIME_TEST
  call stop_test_timer(1)
#endif

  return

end subroutine get_img_frc_ips_distribution

! This is the default async implementation...
!*******************************************************************************
!
! Subroutine:  pvt_get_img_frc_ips_distribution
!
! Description: <TBS>
!
!*******************************************************************************

subroutine pvt_get_img_frc_ips_distribution(atm_cnt, off_tbl, recv_taskmap, &
                                        send_taskmap, send_atm_lst, &
                                        send_atm_cnts, recv_atm_lsts, &
                                        recv_atm_cnts, owned_atm_cnts)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  integer               :: off_tbl(0:numtasks)
  integer               :: recv_taskmap(numtasks - 1)
  integer               :: send_taskmap(numtasks - 1)
  integer               :: send_atm_lst(atm_cnt)
  integer               :: send_atm_cnts(0 : numtasks - 1)
  integer               :: recv_atm_lsts(my_atm_cnt, numtasks - 1)
  integer               :: recv_atm_cnts(0 : numtasks - 1)
  integer               :: owned_atm_cnts(0:numtasks)

! Local variables

  integer               :: wait_call
  integer               :: taskmap_idx
  integer               :: recv_task, send_task
  integer               :: recv_cnt, send_cnt, send_offset
  logical, save         :: recv_atm_cnts_zeroed = .false.
  integer               :: irecv_stat(mpi_status_size)
  integer               :: isend_stat(mpi_status_size, numtasks - 1)
  integer               :: recv_req(numtasks - 1)
  integer               :: send_req(numtasks - 1)

  ! No need to receive from yourself, so there are numtasks - 1 bufs:

  ! Exchange atom lists...

  ! Post the asynchronous receives first:

  if (my_atm_cnt .ne. 0) then

    do taskmap_idx = 1, numtasks - 1

      recv_task = recv_taskmap(taskmap_idx)

      call mpi_irecv(recv_atm_lsts(1, taskmap_idx), my_atm_cnt, mpi_integer, &
                     recv_task, gifd_tag, pmemd_comm, &
                     recv_req(taskmap_idx), err_code_mpi) 
    end do

  end if

  ! Now set up and post the asynchronous sends:

  do taskmap_idx = 1, numtasks - 1

    send_task = send_taskmap(taskmap_idx)

    if (owned_atm_cnts(send_task) .ne. 0) then

      send_offset = off_tbl(send_task) + 1
      send_cnt = send_atm_cnts(send_task)

      call mpi_isend(send_atm_lst(send_offset), send_cnt, &
                     mpi_integer, send_task, gifd_tag, pmemd_comm, &
                     send_req(taskmap_idx), err_code_mpi)

    else

      send_req(taskmap_idx) = MPI_REQUEST_NULL

    end if

  end do

  ! Wait on and process the pending receive requests:

  if (my_atm_cnt .ne. 0) then

    do wait_call = 1, numtasks - 1

      call mpi_waitany(numtasks - 1, recv_req, taskmap_idx, irecv_stat, &
                       err_code_mpi)

      recv_task = recv_taskmap(taskmap_idx)

      call mpi_get_count(irecv_stat, mpi_integer, recv_cnt, err_code_mpi)

      recv_atm_cnts(recv_task) = recv_cnt

    end do

    recv_atm_cnts_zeroed = .false.

  else

    if (.not. recv_atm_cnts_zeroed) then
      recv_atm_cnts(:) = 0
      recv_atm_cnts_zeroed = .true.
    end if

  end if

  ! Wait for all sends to complete:

  call mpi_waitall(numtasks - 1, send_req, isend_stat, err_code_mpi)

  recv_atm_cnts(mytaskid) = 0 

  return

end subroutine pvt_get_img_frc_ips_distribution

!*******************************************************************************
!
! Subroutine:  distribute_img_frcs
!
! Description: <TBS>
!
!*******************************************************************************

subroutine distribute_img_frcs(atm_cnt, img_frc, frc)

  use img_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: img_frc(3, atm_cnt)
  double precision      :: frc(3, atm_cnt)

#ifdef COMM_TIME_TEST
  call start_test_timer(2, 'distribute_img_frcs', 0)
#endif

  call pvt_distribute_img_frcs(atm_cnt, img_frc, frc, gbl_atm_img_map, &
                               pvt_atm_offsets, pvt_taskmap, &
                               pvt_inv_taskmap, &
                               pvt_send_atm_lst, pvt_send_atm_cnts, &
                               pvt_recv_atm_lsts, pvt_recv_atm_cnts, &
                               dbl_mpi_send_buf, dbl_mpi_recv_buf)

#ifdef COMM_TIME_TEST
  call stop_test_timer(2)
#endif

  return

end subroutine distribute_img_frcs

!*******************************************************************************
!
! Subroutine:  pvt_distribute_img_frcs
!
! Description: <TBS>
!
!*******************************************************************************

subroutine pvt_distribute_img_frcs(atm_cnt, img_frc, frc, atm_img_map, &
                                   off_tbl, recv_taskmap, send_taskmap, &
                                   send_atm_lst, send_atm_cnts, &
                                   recv_atm_lsts, recv_atm_cnts, &
                                   send_frc_lst, recv_frc_lsts)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: img_frc(3, atm_cnt)
  double precision      :: frc(3, atm_cnt)
  integer               :: atm_img_map(atm_cnt)
  integer               :: off_tbl(0:numtasks)
  integer               :: recv_taskmap(numtasks - 1)
  integer               :: send_taskmap(numtasks - 1)
  integer               :: send_atm_lst(atm_cnt)
  integer               :: send_atm_cnts(0 : numtasks - 1)
  integer               :: recv_atm_lsts(my_atm_cnt, numtasks - 1)
  integer               :: recv_atm_cnts(0 : numtasks - 1)
  double precision      :: send_frc_lst(3, atm_cnt)
  double precision      :: recv_frc_lsts(3, my_atm_cnt, numtasks - 1)

  ! No need to receive from yourself, so there are numtasks - 1 bufs:

! Local variables:

  integer               :: atm_id, img_id, node, wait_call
  integer               :: atm_lst_idx
  integer               :: i, j
  integer               :: taskmap_idx
  integer               :: recvs_posted, sends_posted
  integer               :: recv_task, send_task
  integer               :: recv_cnt, send_cnt, send_offset
  integer               :: irecv_stat(mpi_status_size)
  integer               :: isend_stat(mpi_status_size, numtasks - 1)
  integer               :: recv_req(numtasks - 1)
  integer               :: send_req(numtasks - 1)

  ! The nonbonded force array should be initialized prior to entry for all
  ! owned atoms.

  recvs_posted = 0

  do taskmap_idx = 1, numtasks - 1

    recv_task = recv_taskmap(taskmap_idx)
    recv_cnt = recv_atm_cnts(recv_task) * 3

    if (recv_cnt .gt. 0) then
      call mpi_irecv(recv_frc_lsts(1, 1, taskmap_idx), recv_cnt, &
                     mpi_double_precision, recv_task, dif_tag, &
                     pmemd_comm, recv_req(taskmap_idx), err_code_mpi)
      recvs_posted = recvs_posted + 1
    else
      recv_req(taskmap_idx) = MPI_REQUEST_NULL
    end if

  end do

  sends_posted = 0

  do taskmap_idx = 1, numtasks - 1

    send_task = send_taskmap(taskmap_idx)
    send_cnt = send_atm_cnts(send_task) * 3
  
    if (send_cnt .gt. 0) then

!BEGIN DBG
!       if (off_tbl(send_task + 1) - off_tbl(send_task) .lt. &
!           send_atm_cnts(send_task)) then
!         write(0,*)'WARNING: send_cnt too big!!!'
!       end if
!END DBG

      send_offset = off_tbl(send_task)
      do i = send_offset + 1, send_offset + send_atm_cnts(send_task)
        atm_id = send_atm_lst(i)
        img_id = atm_img_map(atm_id)
        send_frc_lst(:, i) = img_frc(:, img_id)
      end do

      call mpi_isend(send_frc_lst(1, send_offset + 1), send_cnt, &
                     mpi_double_precision, send_task, dif_tag, pmemd_comm, &
                     send_req(taskmap_idx), err_code_mpi)

      sends_posted = sends_posted + 1
    else
      send_req(taskmap_idx) = MPI_REQUEST_NULL
    end if

  end do

  ! Copy your own img_frc info into the frc buffer:

  send_offset = off_tbl(mytaskid)
  do i = 1, send_atm_cnts(mytaskid)
    atm_id = send_atm_lst(send_offset + i)
    img_id = atm_img_map(atm_id)
    frc(:, atm_id) = frc(:, atm_id) + img_frc(:, img_id)
  end do

  ! Wait on and process any pending receive requests:

  do wait_call = 1, recvs_posted
    call mpi_waitany(numtasks - 1, recv_req, taskmap_idx, irecv_stat, &
                     err_code_mpi)

    recv_task = recv_taskmap(taskmap_idx)

! BEGIN DBG
!   call mpi_get_count(irecv_stat, mpi_double_precision, i, err_code_mpi)
!
!   if (i .ne. 3 * recv_atm_cnts(recv_task)) then
!     write(0,*)'WARNING: recv_atm_cnts value wrong!'
!   end if
!
!   if (i .gt. my_atm_cnt * 3) then
!     write(0,*)'WARNING: buffer overflow on recv!'
!   end if
! END DBG
    
    do i = 1, recv_atm_cnts(recv_task)
      j = recv_atm_lsts(i, taskmap_idx)
      frc(:, j) = frc(:, j) + recv_frc_lsts(:, i, taskmap_idx)
    end do
  end do

  ! Wait for all sends to complete:

  if (sends_posted .gt. 0) then
    call mpi_waitall(numtasks - 1, send_req, isend_stat, err_code_mpi)
  end if

  return

end subroutine pvt_distribute_img_frcs

!*******************************************************************************
!
! Subroutine:  distribute_img_ips_frcs
!
! Description: <TBS>
!
!*******************************************************************************

subroutine distribute_img_ips_frcs(atm_cnt,frc)

  use img_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: frc(3, atm_cnt)

#ifdef COMM_TIME_TEST
  call start_test_timer(2, 'distribute_img_frcs', 0)
#endif

  call pvt_distribute_img_ips_frcs(atm_cnt, frc, gbl_atm_img_map, &
                               pvt_atm_offsets, pvt_taskmap, &
                               pvt_inv_taskmap, &
                               pvt_send_atm_ips_lst, pvt_send_atm_ips_cnts, &
                               pvt_recv_atm_ips_lsts, pvt_recv_atm_ips_cnts, &
                               dbl_mpi_send_buf, dbl_mpi_recv_buf)

#ifdef COMM_TIME_TEST
  call stop_test_timer(2)
#endif

  return

end subroutine distribute_img_ips_frcs

!*******************************************************************************
!
! Subroutine:  pvt_distribute_img_ips_frcs
!
! Description: <TBS>
!
!*******************************************************************************

subroutine pvt_distribute_img_ips_frcs(atm_cnt, frc, atm_img_map, &
                                   off_tbl, recv_taskmap, send_taskmap, &
                                   send_atm_lst, send_atm_cnts, &
                                   recv_atm_lsts, recv_atm_cnts, &
                                   send_frc_lst, recv_frc_lsts)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: frc(3, atm_cnt)
  integer               :: atm_img_map(atm_cnt)
  integer               :: off_tbl(0:numtasks)
  integer               :: recv_taskmap(numtasks - 1)
  integer               :: send_taskmap(numtasks - 1)
  integer               :: send_atm_lst(atm_cnt)
  integer               :: send_atm_cnts(0 : numtasks - 1)
  integer               :: recv_atm_lsts(my_atm_cnt, numtasks - 1)
  integer               :: recv_atm_cnts(0 : numtasks - 1)
  double precision      :: send_frc_lst(3, atm_cnt)
  double precision      :: recv_frc_lsts(3, my_atm_cnt, numtasks - 1)

  ! No need to receive from yourself, so there are numtasks - 1 bufs:

! Local variables:

  integer               :: atm_id, img_id, node, wait_call
  integer               :: atm_lst_idx
  integer               :: i, j
  integer               :: taskmap_idx
  integer               :: recvs_posted, sends_posted
  integer               :: recv_task, send_task
  integer               :: recv_cnt, send_cnt, send_offset
  integer               :: irecv_stat(mpi_status_size)
  integer               :: isend_stat(mpi_status_size, numtasks - 1)
  integer               :: recv_req(numtasks - 1)
  integer               :: send_req(numtasks - 1)

  ! The nonbonded force array should be initialized prior to entry for all
  ! owned atoms.

  recvs_posted = 0

  do taskmap_idx = 1, numtasks - 1

    recv_task = recv_taskmap(taskmap_idx)
    recv_cnt = recv_atm_cnts(recv_task) * 3

    if (recv_cnt .gt. 0) then
      call mpi_irecv(recv_frc_lsts(1, 1, taskmap_idx), recv_cnt, &
                     mpi_double_precision, recv_task, dif_tag, &
                     pmemd_comm, recv_req(taskmap_idx), err_code_mpi)
      recvs_posted = recvs_posted + 1
    else
      recv_req(taskmap_idx) = MPI_REQUEST_NULL
    end if

  end do

  sends_posted = 0

  do taskmap_idx = 1, numtasks - 1

    send_task = send_taskmap(taskmap_idx)
    send_cnt = send_atm_cnts(send_task) * 3
  
    if (send_cnt .gt. 0) then

!BEGIN DBG
!       if (off_tbl(send_task + 1) - off_tbl(send_task) .lt. &
!           send_atm_cnts(send_task)) then
!         write(0,*)'WARNING: send_cnt too big!!!'
!       end if
!END DBG

      send_offset = off_tbl(send_task)
      do i = send_offset + 1, send_offset + send_atm_cnts(send_task)
        atm_id = send_atm_lst(i)
        img_id = atm_img_map(atm_id)
        send_frc_lst(:, i) = frc(:, atm_id)
        frc(:,atm_id) = 0.d0
      end do

      call mpi_isend(send_frc_lst(1, send_offset + 1), send_cnt, &
                     mpi_double_precision, send_task, dif_tag, pmemd_comm, &
                     send_req(taskmap_idx), err_code_mpi)

      sends_posted = sends_posted + 1
    else
      send_req(taskmap_idx) = MPI_REQUEST_NULL
    end if

  end do

!not  ! Copy your own img_frc info into the frc buffer:

!  send_offset = off_tbl(mytaskid)
!  do i = 1, send_atm_cnts(mytaskid)
!    atm_id = send_atm_lst(send_offset + i)
!    img_id = atm_img_map(atm_id)
!    frc(:, atm_id) = frc(:, atm_id) + img_frc(:, img_id)
!  end do

  ! Wait on and process any pending receive requests:

  do wait_call = 1, recvs_posted
    call mpi_waitany(numtasks - 1, recv_req, taskmap_idx, irecv_stat, &
                     err_code_mpi)

    recv_task = recv_taskmap(taskmap_idx)

! BEGIN DBG
!   call mpi_get_count(irecv_stat, mpi_double_precision, i, err_code_mpi)
!
!   if (i .ne. 3 * recv_atm_cnts(recv_task)) then
!     write(0,*)'WARNING: recv_atm_cnts value wrong!'
!   end if
!
!   if (i .gt. my_atm_cnt * 3) then
!     write(0,*)'WARNING: buffer overflow on recv!'
!   end if
! END DBG
    
    do i = 1, recv_atm_cnts(recv_task)
      j = recv_atm_lsts(i, taskmap_idx)
      frc(:, j) = frc(:, j) + recv_frc_lsts(:, i, taskmap_idx)
    end do
  end do

  ! Wait for all sends to complete:

  if (sends_posted .gt. 0) then
    call mpi_waitall(numtasks - 1, send_req, isend_stat, err_code_mpi)
  end if

  return

end subroutine pvt_distribute_img_ips_frcs


!*******************************************************************************
!
! Subroutine:  distribute_img_ti_frcs
!
! Description: The only difference from the original distribute_img_frcs is
! the force arrays are split in two to match ti_nb_frc and ti_img_frc. Doing
! this here rather than splicing the ti_*_frc arrays improves the effeciency
! at the cost of some duplicate code.
!
!*******************************************************************************

subroutine distribute_img_ti_frcs(atm_cnt, img_frc, frc)

  use img_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: img_frc(2, 3, atm_cnt)
  double precision      :: frc(2, 3, atm_cnt)
#ifdef COMM_TIME_TEST
  call start_test_timer(2, 'distribute_img_frcs', 0)
#endif
  call pvt_distribute_img_ti_frcs(atm_cnt, img_frc, frc, gbl_atm_img_map, &
                               pvt_atm_offsets, pvt_taskmap, &
                               pvt_inv_taskmap, &
                               pvt_send_atm_lst, pvt_send_atm_cnts, &
                               pvt_recv_atm_lsts, pvt_recv_atm_cnts, &
                               dbl_mpi_send_buf, dbl_mpi_recv_buf)

#ifdef COMM_TIME_TEST
  call stop_test_timer(2)
#endif

  return

end subroutine distribute_img_ti_frcs

!*******************************************************************************
!
! Subroutine:  pvt_distribute_img_ti_frcs
!
! Description: <TBS>
!
!*******************************************************************************

subroutine pvt_distribute_img_ti_frcs(atm_cnt, img_frc, frc, atm_img_map, &
                                   off_tbl, recv_taskmap, send_taskmap, &
                                   send_atm_lst, send_atm_cnts, &
                                   recv_atm_lsts, recv_atm_cnts, &
                                   send_frc_lst, recv_frc_lsts)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: img_frc(2, 3, atm_cnt)
  double precision      :: frc(2, 3, atm_cnt)
  integer               :: atm_img_map(atm_cnt)
  integer               :: off_tbl(0:numtasks)
  integer               :: recv_taskmap(numtasks - 1)
  integer               :: send_taskmap(numtasks - 1)
  integer               :: send_atm_lst(atm_cnt)
  integer               :: send_atm_cnts(0 : numtasks - 1)
  integer               :: recv_atm_lsts(my_atm_cnt, numtasks - 1)
  integer               :: recv_atm_cnts(0 : numtasks - 1)
  double precision      :: send_frc_lst(2, 3, atm_cnt)
  double precision      :: recv_frc_lsts(2, 3, my_atm_cnt, numtasks - 1)

  ! No need to receive from yourself, so there are numtasks - 1 bufs:

! Local variables:

  integer               :: atm_id, img_id, node, wait_call
  integer               :: atm_lst_idx
  integer               :: i, j
  integer               :: taskmap_idx
  integer               :: recvs_posted, sends_posted
  integer               :: recv_task, send_task
  integer               :: recv_cnt, send_cnt, send_offset
  integer               :: irecv_stat(mpi_status_size)
  integer               :: isend_stat(mpi_status_size, numtasks - 1)
  integer               :: recv_req(numtasks - 1)
  integer               :: send_req(numtasks - 1)

  ! The nonbonded force array should be initialized prior to entry for all
  ! owned atoms.

  recvs_posted = 0

  do taskmap_idx = 1, numtasks - 1

    recv_task = recv_taskmap(taskmap_idx)
    recv_cnt = recv_atm_cnts(recv_task) * 3 * 2

    if (recv_cnt .gt. 0) then
      call mpi_irecv(recv_frc_lsts(1, 1, 1, taskmap_idx), recv_cnt, &
                     mpi_double_precision, recv_task, dif_tag, &
                     pmemd_comm, recv_req(taskmap_idx), err_code_mpi)
      recvs_posted = recvs_posted + 1
    else
      recv_req(taskmap_idx) = MPI_REQUEST_NULL
    end if

  end do

  sends_posted = 0

  do taskmap_idx = 1, numtasks - 1

    send_task = send_taskmap(taskmap_idx)
    send_cnt = send_atm_cnts(send_task) * 3 * 2
  
    if (send_cnt .gt. 0) then

!BEGIN DBG
!       if (off_tbl(send_task + 1) - off_tbl(send_task) .lt. &
!           send_atm_cnts(send_task)) then
!         write(0,*)'WARNING: send_cnt too big!!!'
!       end if
!END DBG

      send_offset = off_tbl(send_task)
      do i = send_offset + 1, send_offset + send_atm_cnts(send_task)
        atm_id = send_atm_lst(i)
        img_id = atm_img_map(atm_id)
        send_frc_lst(:, :, i) = img_frc(:, :, img_id)
      end do

      call mpi_isend(send_frc_lst(1, 1, send_offset + 1), send_cnt, &
                     mpi_double_precision, send_task, dif_tag, pmemd_comm, &
                     send_req(taskmap_idx), err_code_mpi)

      sends_posted = sends_posted + 1
    else
      send_req(taskmap_idx) = MPI_REQUEST_NULL
    end if

  end do

  ! Copy your own img_frc info into the frc buffer:

  send_offset = off_tbl(mytaskid)
  do i = 1, send_atm_cnts(mytaskid)
    atm_id = send_atm_lst(send_offset + i)
    img_id = atm_img_map(atm_id)
    frc(:, :, atm_id) = frc(:, :, atm_id) + img_frc(:, :, img_id)
  end do

  ! Wait on and process any pending receive requests:

  do wait_call = 1, recvs_posted
    call mpi_waitany(numtasks - 1, recv_req, taskmap_idx, irecv_stat, &
                     err_code_mpi)

    recv_task = recv_taskmap(taskmap_idx)

! BEGIN DBG
!   call mpi_get_count(irecv_stat, mpi_double_precision, i, err_code_mpi)
!
!   if (i .ne. 3 * recv_atm_cnts(recv_task)) then
!     write(0,*)'WARNING: recv_atm_cnts value wrong!'
!   end if
!
!   if (i .gt. my_atm_cnt * 3) then
!     write(0,*)'WARNING: buffer overflow on recv!'
!   end if
! END DBG
    
    do i = 1, recv_atm_cnts(recv_task)
      j = recv_atm_lsts(i, taskmap_idx)
      frc(:, :, j) = frc(:, :, j) + recv_frc_lsts(:, :, i, taskmap_idx)
    end do
  end do

  ! Wait for all sends to complete:

  if (sends_posted .gt. 0) then
    call mpi_waitall(numtasks - 1, send_req, isend_stat, err_code_mpi)
  end if

  return

end subroutine pvt_distribute_img_ti_frcs

!*******************************************************************************
!
! Subroutine:  distribute_crds
!
! Description: <TBS>
!
!*******************************************************************************

subroutine distribute_crds(atm_cnt, crd)

  use parallel_dat_mod

  implicit none

  integer               :: atm_cnt
  double precision      :: crd(3, atm_cnt)

#ifdef COMM_TIME_TEST
  call start_test_timer(3, 'distribute_crds', 0)
#endif

  call pvt_distribute_crds(atm_cnt, crd, pvt_atm_offsets, pvt_inv_taskmap, &
                           pvt_taskmap, pvt_send_atm_lst, pvt_send_atm_cnts, &
                           pvt_recv_atm_lsts, pvt_recv_atm_cnts, &
                           dbl_mpi_send_buf, dbl_mpi_recv_buf)

#ifdef COMM_TIME_TEST
  call stop_test_timer(3)
#endif

  return

end subroutine distribute_crds

!*******************************************************************************
!
! Subroutine:  pvt_distribute_crds
!
! Description: <TBS>
!
!*******************************************************************************

subroutine pvt_distribute_crds(atm_cnt, crd, off_tbl, recv_taskmap, &
                               send_taskmap, recv_atm_lst, recv_atm_cnts, &
                               send_atm_lsts, send_atm_cnts, &
                               send_crd_lsts, recv_crd_lst)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: crd(3, atm_cnt)
  integer               :: off_tbl(0:numtasks)
  integer               :: recv_taskmap(numtasks - 1)
  integer               :: send_taskmap(numtasks - 1)
  integer               :: recv_atm_lst(atm_cnt)
  integer               :: recv_atm_cnts(0 : numtasks - 1)
  integer               :: send_atm_lsts(my_atm_cnt, numtasks - 1)
  integer               :: send_atm_cnts(0 : numtasks - 1)
  double precision      :: send_crd_lsts(3, my_atm_cnt, numtasks - 1)
  double precision      :: recv_crd_lst(3, atm_cnt)


! Local variables

  integer               :: wait_call
  integer               :: i
  integer               :: taskmap_idx
  integer               :: recvs_posted, sends_posted
  integer               :: recv_task, send_task
  integer               :: recv_cnt, send_cnt, recv_offset
  integer               :: irecv_stat(mpi_status_size)
  integer               :: isend_stat(mpi_status_size, numtasks - 1)
  integer               :: recv_req(numtasks - 1)
  integer               :: send_req(numtasks - 1)

  ! No need to send to yourself, so there are numtasks - 1 bufs.

  ! The global names pvt_send_atm_lst, cit_sent_atm_cnts, cit_recv_atm_lst,
  ! pvt_recv_atm_cnts all refer to the process required to distribute image
  ! forces to their owners.  Here we are gathering coordinates to be used with
  ! images from the atom owners, so the process is inverted (ie., we actually
  ! use send_atm_lst and send_atm_cnts to receive coordinates, and vice versa). 
  ! We therefore invert the names in the call list to make the names in this
  ! routine consistent.

  ! Post the asynchronous receives first:

  recvs_posted = 0

  do taskmap_idx = 1, numtasks - 1

    recv_task = recv_taskmap(taskmap_idx)
    recv_cnt = recv_atm_cnts(recv_task) * 3

    if (recv_cnt .gt. 0) then

      call mpi_irecv(recv_crd_lst(1, off_tbl(recv_task) + 1), recv_cnt, &
                     mpi_double_precision, recv_task, dc_tag, &
                     pmemd_comm, recv_req(taskmap_idx), err_code_mpi)
      recvs_posted = recvs_posted + 1
    else
      recv_req(taskmap_idx) = MPI_REQUEST_NULL
    end if

  end do

    ! Now set up and post the asynchronous sends:

  sends_posted = 0

  do taskmap_idx = 1, numtasks - 1

    send_task = send_taskmap(taskmap_idx)

    do i = 1, send_atm_cnts(send_task)
      send_crd_lsts(:, i, taskmap_idx) = crd(:, send_atm_lsts(i, taskmap_idx))
    end do
  
    send_cnt = send_atm_cnts(send_task) * 3

! BEGIN DBG
!       if (my_atm_cnt .lt. send_atm_cnts(send_task)) then
!         write(0,*)'WARNING: send_cnt too big, distribute_crds!!!'
!       end if
! END DBG

    if (send_cnt .gt. 0) then
      call mpi_isend(send_crd_lsts(1, 1, taskmap_idx), send_cnt, &
                     mpi_double_precision, send_task, dc_tag, pmemd_comm, &
                     send_req(taskmap_idx), err_code_mpi)
      sends_posted = sends_posted + 1
    else
      send_req(taskmap_idx) = MPI_REQUEST_NULL
    end if

  end do

  ! Wait on and process the pending receive requests:

  do wait_call = 1, recvs_posted

    call mpi_waitany(numtasks - 1, recv_req, taskmap_idx, irecv_stat, &
                     err_code_mpi)

    recv_task = recv_taskmap(taskmap_idx)
    recv_offset = off_tbl(recv_task)

! BEGIN DBG
!   call mpi_get_count(irecv_stat, mpi_double_precision, i, err_code_mpi)
!
!   if (i .ne. 3 * recv_atm_cnts(recv_task)) then
!     write(0,*)'WARNING: recv_atm_cnts value wrong, distribute_crds!'
!   end if
!
!   if (recv_atm_cnts(recv_task) .gt. &
!       off_tbl(recv_task + 1) - off_tbl(recv_task)) then
!     write(0,*)'WARNING: buffer overflow on recv, distribute_crds!'
!   end if
! END DBG
    
    do i = recv_offset + 1, recv_offset + recv_atm_cnts(recv_task)
      crd(:, recv_atm_lst(i)) = recv_crd_lst(:, i)
    end do

  end do

  ! Wait for all sends to complete:

  if (sends_posted .gt. 0) then
    call mpi_waitall(numtasks - 1, send_req, isend_stat, err_code_mpi)
  end if

  return

end subroutine pvt_distribute_crds

!*******************************************************************************
!
! Subroutine:  zero_extra_used_atm_img_frcs
!
! Description: <TBS>
!              
!*******************************************************************************
subroutine zero_extra_used_atm_img_frcs(img_frc)

  use img_mod

  implicit none

! Formal arguments:

  double precision      :: img_frc(3, *)
    
! Local variables:

  integer               :: i, j

  do j = 1, extra_used_atm_cnt
    i = gbl_atm_img_map(pvt_extra_used_atms(j))
    img_frc(:, i) = 0.d0
  end do

  return

end subroutine zero_extra_used_atm_img_frcs

!*******************************************************************************
!
! Subroutine:  zero_extra_used_atm_ti_img_frcs
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine zero_extra_used_atm_ti_img_frcs(img_frc)

  use img_mod

  implicit none

! Formal arguments:

  double precision      :: img_frc(2, 3, *)
    
! Local variables:

  integer               :: i, j

  do j = 1, extra_used_atm_cnt
    i = gbl_atm_img_map(pvt_extra_used_atms(j))
    img_frc(:, :, i) = 0.d0
  end do

  return

end subroutine zero_extra_used_atm_ti_img_frcs

!*******************************************************************************
!
! Subroutine:  mpi_allgathervec
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine mpi_allgathervec(atm_cnt, vec)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: vec(3, atm_cnt)

#ifdef COMM_TIME_TEST
  call start_test_timer(4, 'mpi_allgathervec', 0)
#endif

  call pvt_mpi_allgathervec(atm_cnt, vec, dbl_mpi_send_buf, dbl_mpi_recv_buf)

#ifdef COMM_TIME_TEST
  call stop_test_timer(4)
#endif

  return

end subroutine mpi_allgathervec

!*******************************************************************************
!
! Subroutine:  pvt_mpi_allgathervec
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine pvt_mpi_allgathervec(atm_cnt, vec, send_buf, recv_buf)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: vec(3, atm_cnt)
  double precision      :: send_buf(3, *)
  double precision      :: recv_buf(3, *)

! Local variables:

  integer               :: atm_lst_idx, atm_idx
  integer               :: my_sendoff
  integer               :: node
  integer               :: cur_node_off(0: numtasks - 1)

  ! We marshall, send/recv, and unmarshall from all processes including the
  ! current process in order to keep the current process from having to do a
  ! conditional across all the data in the unmarshalling.

  my_sendoff = pvt_atm_offsets(mytaskid)
  
  do atm_lst_idx = 1, my_atm_cnt
    atm_idx = gbl_my_atm_lst(atm_lst_idx)
    send_buf(:, my_sendoff + atm_lst_idx) = vec(:, atm_idx)
  end do

  call mpi_allgatherv(send_buf(1, my_sendoff + 1), &
                      pvt_vec_rcvcnts(mytaskid), &
                      mpi_double_precision, &
                      recv_buf, pvt_vec_rcvcnts, pvt_vec_offsets, &
                      mpi_double_precision, pmemd_comm, err_code_mpi)

  cur_node_off(0:numtasks - 1) = pvt_atm_offsets(0:numtasks - 1) + 1
  
  do atm_idx = 1, atm_cnt
    node = gbl_atm_owner_map(atm_idx)
    vec(:, atm_idx) = recv_buf(:, cur_node_off(node))
    cur_node_off(node) = cur_node_off(node) + 1
  end do

  return

end subroutine pvt_mpi_allgathervec

!*******************************************************************************
!
! Subroutine:  mpi_gathervec
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine mpi_gathervec(atm_cnt, vec)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: vec(3, atm_cnt)

#ifdef COMM_TIME_TEST
  call start_test_timer(5, 'mpi_gathervec', 0)
#endif

  call pvt_mpi_gathervec(atm_cnt, vec, &
                         dbl_mpi_send_buf, dbl_mpi_recv_buf)

#ifdef COMM_TIME_TEST
  call stop_test_timer(5)
#endif

  return

end subroutine mpi_gathervec

!*******************************************************************************
!
! Subroutine:  pvt_mpi_gathervec
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine pvt_mpi_gathervec(atm_cnt, vec, send_buf, recv_buf)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: vec(3, atm_cnt)
  double precision      :: send_buf(3, *)
  double precision      :: recv_buf(3, *)

! Local variables:

  integer               :: atm_lst_idx, atm_idx
  integer               :: my_sendoff
  integer               :: node
  integer               :: cur_node_off(0: numtasks - 1)

  ! We marshall, send/recv, and unmarshall from all processes including root, in
  ! order to keep the root from having to do a conditional across all the
  ! data in the unmarshalling.

  my_sendoff = pvt_atm_offsets(mytaskid)
  
  do atm_lst_idx = 1, my_atm_cnt
    atm_idx = gbl_my_atm_lst(atm_lst_idx)
    send_buf(:, my_sendoff + atm_lst_idx) = vec(:, atm_idx)
  end do

  call mpi_gatherv(send_buf(1, my_sendoff + 1), &
                   pvt_vec_rcvcnts(mytaskid), &
                   mpi_double_precision, &
                   recv_buf, pvt_vec_rcvcnts, pvt_vec_offsets, &
                   mpi_double_precision, &
                   0, pmemd_comm, err_code_mpi)

  if (mytaskid .eq. 0) then
    cur_node_off(0:numtasks - 1) = pvt_atm_offsets(0:numtasks - 1) + 1
    do atm_idx = 1, atm_cnt
      node = gbl_atm_owner_map(atm_idx)
      vec(:, atm_idx) = recv_buf(:, cur_node_off(node))
      cur_node_off(node) = cur_node_off(node) + 1
    end do
  end if

  return

end subroutine pvt_mpi_gathervec

!*******************************************************************************
!
! Subroutine:  do_atm_redistribution
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine do_atm_redistribution(atm_cnt, write_log, step_ctr)

  use dynamics_mod
  use dynamics_dat_mod
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
  use shake_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  logical               :: write_log
  integer               :: step_ctr
    
! Local variables:

  integer               :: alloc_failed
  integer               :: i
  integer               :: num_ints, num_reals  ! values are not kept.

! if (master) write(0,*)'| DBG: Doing atom redistribution.'

  num_ints = 0
  num_reals = 0

  call do_atm_distribution(atm_cnt, num_ints, num_reals, write_log, step_ctr)

! Redo the shake per-process setup:

#ifndef CUDA
#ifdef _OPENMP_
if(.not. usemidpoint) then
 if(using_gb_potential) then
  call claim_my_fastwat_residues_gb(num_ints, num_reals)
  call claim_my_nonfastwat_bonds_gb(num_ints, num_reals)
 else 
  call claim_my_fastwat_residues(num_ints, num_reals)
  call claim_my_nonfastwat_bonds(num_ints, num_reals)
 end if
else
  call claim_my_fastwat_residues(num_ints, num_reals)
  call claim_my_nonfastwat_bonds(num_ints, num_reals)
endif ! not usemidpoint
#else
  call claim_my_fastwat_residues(num_ints, num_reals)
  call claim_my_nonfastwat_bonds(num_ints, num_reals)
#endif /*_OPENMP_*/
#endif

! Recreate the entire molecule COM array.  You don't actually need the whole
! thing under mpi, but you have all the crds from which it is derived, and
! this code should execute very infrequently.

  if (ntp .gt. 0 .and. imin .eq. 0) then
    call get_all_mol_com(atm_crd, atm_mass, gbl_mol_mass_inv, gbl_mol_com)
  end if

! Clear the force array.  This avoids potential problems with non-zero data
! for atoms we don't own in the nmr routines and possibly elsewhere.

  atm_frc(:,:) = 0.d0

  return

end subroutine do_atm_redistribution

!*******************************************************************************
!
! Subroutine:  do_atm_distribution
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine do_atm_distribution(atm_cnt, num_ints, num_reals, &
                               write_log, step_ctr)

  use angles_mod
  use angles_ub_mod,     only : angles_ub_setup
  use cmap_mod,          only : cmap_setup
  use bonds_mod
  use cit_mod
  use dihedrals_mod
  use dihedrals_imp_mod, only : dihedrals_imp_setup
  use extra_pnts_nb14_mod
  use inpcrd_dat_mod
  use nb_exclusions_mod
  use pbc_mod
  use pmemd_lib_mod
  use parallel_dat_mod
  use nmr_calls_mod
  use charmm_mod,        only : charmm_active
  use mdin_ctrl_dat_mod, only : ips, nmropt
  use ti_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer                       :: atm_cnt
  integer, intent(in out)       :: num_ints, num_reals
  logical                       :: write_log
  integer                       :: step_ctr
    
! Local variables:

  integer               :: alloc_failed
  integer               :: i
  integer               :: extra_used_atms(atm_cnt)
  integer               :: use_atm_map(atm_cnt)

  double precision      :: fraction(3, atm_cnt) ! in range 0.0 - +0.999...
  integer               :: crd_idx_lst_tbl(0 : cit_tbl_x_dim - 1, &
                                           0 : cit_tbl_y_dim - 1, &
                                           0 : cit_tbl_z_dim - 1)
  type(atm_lst_rec)     :: atm_lst(atm_cnt)

! Divide atoms up among the processors.  The atom division is redone
! periodically under cit, and is prf-based, with locality.

! if (master) write(0,*)'| DBG: Doing atom distribution.'

  call get_fract_crds(atm_cnt, atm_crd, fraction)

  call setup_crd_idx_tbl(atm_cnt, fraction, crd_idx_lst_tbl, atm_lst)

  call divide_atoms(atm_cnt, fraction, crd_idx_lst_tbl, atm_lst, write_log, &
                    step_ctr)

  log_owned_atm_cnt = dble(my_atm_cnt)

  use_atm_map(:) = 0

  call bonds_setup(num_ints, num_reals, use_atm_map)
  call angles_setup(num_ints, num_reals, use_atm_map)
  call dihedrals_setup(num_ints, num_reals, use_atm_map)
#ifdef CUDA 
  if (nmropt .ne. 0) then
    call cuda_nmr_setup()         
  end if
#endif
  if (charmm_active) then
    call angles_ub_setup(num_ints, num_reals, use_atm_map)
    call dihedrals_imp_setup(num_ints, num_reals, use_atm_map)
  end if
  ! if CMAP is active
  if (cmap_term_count > 0) then
    call cmap_setup(num_ints, num_reals, use_atm_map)
  endif
  call nb14_setup(num_ints, num_reals, use_atm_map)

  call make_nb_adjust_pairlst(my_atm_cnt, gbl_my_atm_lst, &
                              gbl_atm_owner_map, use_atm_map, &
                              gbl_nb_adjust_pairlst, &
                              atm_nb_maskdata, atm_nb_mask)

  extra_used_atms(:) = 0
  extra_used_atm_cnt = 0

  do i = 1, atm_cnt
    if (use_atm_map(i) .ne. 0) then
      if (gbl_atm_owner_map(i) .ne. mytaskid) then
        extra_used_atm_cnt = extra_used_atm_cnt + 1
        extra_used_atms(extra_used_atm_cnt) = i
      end if
    end if
  end do

  if (extra_used_atm_cnt .gt. 0) then
    if (allocated(pvt_extra_used_atms)) then
      if (size(pvt_extra_used_atms) .lt. extra_used_atm_cnt) then
        deallocate(pvt_extra_used_atms)
        allocate(pvt_extra_used_atms(extra_used_atm_cnt), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
      end if
    else
      allocate(pvt_extra_used_atms(extra_used_atm_cnt), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
    end if

    pvt_extra_used_atms(1:extra_used_atm_cnt) = &
      extra_used_atms(1:extra_used_atm_cnt)
  end if

! write(0,*)'DBG: Task ', mytaskid, ' has ', extra_used_atm_cnt, &
!           'extra used atoms'

  if (allocated(pvt_recv_atm_lsts)) then
    if (size(pvt_recv_atm_lsts, 1) .lt. my_atm_cnt) then
      deallocate(pvt_recv_atm_lsts)
      allocate(pvt_recv_atm_lsts(my_atm_cnt, numtasks - 1), &
               stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
    end if
  else
    allocate(pvt_recv_atm_lsts(my_atm_cnt, numtasks - 1), &
             stat = alloc_failed)
    if (alloc_failed .ne. 0) call setup_alloc_error
  end if

  pvt_recv_atm_lsts(:,:) = 0

  if (ips .gt. 0)then
  if (allocated(pvt_recv_atm_ips_lsts)) then
    if (size(pvt_recv_atm_ips_lsts, 1) .lt. my_atm_cnt) then
      deallocate(pvt_recv_atm_ips_lsts)
      allocate(pvt_recv_atm_ips_lsts(my_atm_cnt, numtasks - 1), &
               stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
    end if
  else
    allocate(pvt_recv_atm_ips_lsts(my_atm_cnt, numtasks - 1), &
             stat = alloc_failed)
    if (alloc_failed .ne. 0) call setup_alloc_error
  end if

  pvt_recv_atm_ips_lsts(:,:) = 0
  endif

  if (ti_mode .eq. 0) then
     call set_minimum_mpi_bufs_size(max(3 * atm_cnt, &
                                    3 * my_atm_cnt * (numtasks-1)),&
                                    num_reals)
  else
     call set_minimum_mpi_bufs_size(max(2 * 3 * atm_cnt, &
                                    2 * 3 * my_atm_cnt * (numtasks-1)),&
                                    num_reals)
  end if

  return

end subroutine do_atm_distribution

!*******************************************************************************
!
! Subroutine:  divide_atoms
!
! Description:  Set up atom workload division for parallel processing.
!              
!*******************************************************************************

subroutine divide_atoms(atm_cnt, fraction, crd_idx_lst_tbl, atm_lst, &
                        write_log, step_ctr)

  use cit_mod
  use img_mod
  use dynamics_dat_mod
  use extra_pnts_nb14_mod
  use gbl_constants_mod
  use pme_recip_dat_mod
  use pme_slab_fft_mod
  use inpcrd_dat_mod
  use parallel_dat_mod
  use prfs_mod
  use mol_list_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
  use mdin_ctrl_dat_mod, only : ntp, imin, loadbal_verbose

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: fraction(3, atm_cnt) ! in range 0.0 - +0.999...
  integer               :: crd_idx_lst_tbl(0 : cit_tbl_x_dim - 1, &
                                           0 : cit_tbl_y_dim - 1, &
                                           0 : cit_tbl_z_dim - 1)
  type(atm_lst_rec)     :: atm_lst(atm_cnt)
  logical               :: write_log
  integer               :: step_ctr


! Local variables:

  integer               :: alloc_failed
  integer               :: atm_id
  integer               :: first_atm_id
  integer               :: i, j, k
  integer               :: molfrag_idx
  integer               :: frag_mol_idx
  integer               :: img_id, img_lo, img_hi
  integer               :: nxt_idx
  integer               :: mol_id
  integer               :: mol_offset, mol_listcnt
  integer               :: prf_offset, prf_listcnt
  integer               :: first_prf_id, prf_id 
  integer               :: ori_taskid 
  integer               :: srvr_taskid, srvr_off 
  integer               :: taskid 
  integer               :: added_atms_cnt
  integer               :: mol_atms_cnt
  integer, save         :: distrib_cnt = 0
  integer               :: srvr_my_owned_atm_cnt
  integer               :: target_cnt(0:numtasks - 1)

  integer               :: img_owner_map(atm_cnt)

  integer               :: atm_img_map(atm_cnt) ! NOTE local copy...

  type(cit_tbl_rec)     :: crd_idx_tbl(0 : cit_tbl_x_dim - 1, &
                                       0 : cit_tbl_y_dim - 1, &
                                       0 : cit_tbl_z_dim - 1)

  distrib_cnt = distrib_cnt + 1 ! count of calls to this routine

  ! Initialize the atom owner map to values that will catch bugs:

  gbl_atm_owner_map(:) = -1

  ! Make the image owner map, as well as the target cnts:

  do taskid = 0, numtasks - 1
    img_lo = gbl_owned_imgs_tbl(1, taskid) + 1
    img_hi = img_lo + gbl_owned_imgs_tbl(2, taskid) - 1
    if (img_lo .le. img_hi) img_owner_map(img_lo:img_hi) = taskid
    target_cnt(taskid) = img_hi - img_lo + 1
  end do

  img_hi = 0

  do k = 0, cit_tbl_z_dim - 1
    do j = 0, cit_tbl_y_dim - 1
      do i = 0, cit_tbl_x_dim - 1

        nxt_idx = crd_idx_lst_tbl(i, j, k)

        if (nxt_idx .ne. 0) then

          img_hi = img_hi + 1
          img_lo = img_hi

          do
            atm_img_map(atm_lst(nxt_idx)%idx) = img_hi
            nxt_idx = atm_lst(nxt_idx)%nxt
            if (nxt_idx .ne. 0) then
              img_hi = img_hi + 1
            else
              exit
            end if
          end do

          crd_idx_tbl(i, j, k)%img_lo = img_lo
          crd_idx_tbl(i, j, k)%img_hi = img_hi

          nxt_idx = crd_idx_lst_tbl(i, j, k)

        else
          crd_idx_tbl(i, j, k)%img_lo = 0
          crd_idx_tbl(i, j, k)%img_hi = -1
        end if

      end do
    end do
  end do

  ! For NVE/NVT:
  ! Find the first atom in each prf, and map this atom to its image.  Then find
  ! the image owner, and claim all the atoms in the prf for the owner.  We use
  ! a simpler mechanism than in the past, not worrying about targets for owned
  ! atom counts.

  my_atm_cnt = 0
  my_mol_cnt = 0
  my_frag_mol_cnt = 0
  frag_mol_idx = 0
  pvt_owned_atm_cnts(:) = 0
  pvt_vec_rcvcnts(:) = 0

  if (ntp .gt. 0 .and. imin .eq. 0) then ! NTP ensemble.  Need molecule data.

    ! Constant pressure is a bear for a variety of reasons.  One big problem
    ! is uneven atom division caused by big molecules.  To counteract this
    ! effect, we have implemented a molecule fragment model, and we try to not
    ! exceed a target count by spilling over into other nodes when necessary.
    ! The target count is set to match the image count for a task. We do not
    ! overflow to tasks with a recip force workload if we have a choice.

!BEGIN DBG
!   if (master) write(0,*)'| DBG: Target node atm count =', target_cnt(:)
!END DBG

    if (recip_numtasks .eq. numtasks) then

      do mol_id = 1, gbl_mol_cnt
  
        added_atms_cnt = gbl_mol_atms_listdata(mol_id)%cnt
        mol_offset = gbl_mol_prfs_listdata(mol_id)%offset
        mol_listcnt = gbl_mol_prfs_listdata(mol_id)%cnt

        first_prf_id = gbl_mol_prfs_lists(mol_offset + 1)

        atm_id = gbl_prf_lists(gbl_prf_listdata(first_prf_id)%offset + 1)

        if (added_atms_cnt .le. max_unfrag_mol_atms) then

          img_id = atm_img_map(atm_id)
          taskid = img_owner_map(img_id)
        
          if (pvt_owned_atm_cnts(taskid) + added_atms_cnt .gt. &
              target_cnt(taskid)) then
  
            ! If there is no previous assignment to this task, go ahead and
            ! assign this; otherwise look for a more appropriate task.
            ! To support high_scaling algorithms, we save the original task
            ! assignment so we can fall back on it should the algorithm below
            ! pick a non-frc-ene task.
  
            if (pvt_owned_atm_cnts(taskid) .ne. 0) then
              ori_taskid = taskid
              do i = 1, numtasks
                taskid = taskid + 1
                if (taskid .ge. numtasks) taskid = 0
                if (is_frc_ene_task(taskid)) then
                  if (pvt_owned_atm_cnts(taskid) .eq. 0 .or. &
                      pvt_owned_atm_cnts(taskid) + added_atms_cnt .le. &
                      target_cnt(taskid)) exit
                end if
              end do !i = 1, numtasks
              if (.not. is_frc_ene_task(taskid)) taskid = ori_taskid
            end if

          end if
            
          pvt_owned_atm_cnts(taskid) = &
            pvt_owned_atm_cnts(taskid) + added_atms_cnt
  
          if (taskid .ne. mytaskid) then
            do i = mol_offset + 1, mol_offset + mol_listcnt
              prf_id = gbl_mol_prfs_lists(i)
              prf_offset = gbl_prf_listdata(prf_id)%offset
              prf_listcnt = gbl_prf_listdata(prf_id)%cnt
              do j = prf_offset + 1, prf_offset + prf_listcnt
                atm_id = gbl_prf_lists(j)
                gbl_atm_owner_map(atm_id) = taskid
              end do !j = prf_offset + 1, prf_offset + prf_listcnt
            end do !i = mol_offset + 1, mol_offset + mol_listcnt
          else
            my_mol_cnt = my_mol_cnt + 1
            gbl_my_mol_lst(my_mol_cnt) = mol_id
            do i = mol_offset + 1, mol_offset + mol_listcnt
              prf_id = gbl_mol_prfs_lists(i)
              prf_offset = gbl_prf_listdata(prf_id)%offset
              prf_listcnt = gbl_prf_listdata(prf_id)%cnt
              do j = prf_offset + 1, prf_offset + prf_listcnt
                atm_id = gbl_prf_lists(j)
                gbl_atm_owner_map(atm_id) = taskid
                my_atm_cnt = my_atm_cnt + 1
                gbl_my_atm_lst(my_atm_cnt) = atm_id
              end do !j = prf_offset + 1, prf_offset + prf_listcnt
            end do !i = mol_offset + 1, mol_offset + mol_listcnt
          end if !(taskid .ne. mytaskid)

        else

          frag_mol_idx = frag_mol_idx + 1

          do molfrag_idx = gbl_frag_mols(frag_mol_idx)%first_molfrag_idx, &
                           gbl_frag_mols(frag_mol_idx)%first_molfrag_idx + &
                           gbl_frag_mols(frag_mol_idx)%molfrag_cnt - 1

            mol_offset = gbl_molfrags(molfrag_idx)%offset
            mol_listcnt = gbl_molfrags(molfrag_idx)%cnt

            added_atms_cnt = 0

            do i = mol_offset + 1, mol_offset + mol_listcnt
              prf_id = gbl_mol_prfs_lists(i)
              added_atms_cnt = added_atms_cnt + gbl_prf_listdata(prf_id)%cnt
            end do !i = mol_offset + 1, mol_offset + mol_listcnt

            first_prf_id = gbl_mol_prfs_lists(mol_offset + mol_listcnt)
            first_atm_id = &
              gbl_prf_lists(gbl_prf_listdata(first_prf_id)%offset+1)

            img_id = atm_img_map(first_atm_id)
            taskid = img_owner_map(img_id)

            if (pvt_owned_atm_cnts(taskid) + added_atms_cnt .gt. &
                target_cnt(taskid)) then
  
              ! If there is no previous assignment to this task, go ahead and
              ! assign this; otherwise look for a more appropriate task.
              ! To support high_scaling algorithms, we save the original task
              ! assignment so we can fall back on it should the algorithm below
              ! pick a non-frc-ene task.
  
              if (pvt_owned_atm_cnts(taskid) .ne. 0) then
                ori_taskid = taskid
                do i = 1, numtasks
                  taskid = taskid + 1
                  if (taskid .ge. numtasks) taskid = 0
                  if (is_frc_ene_task(taskid)) then
                    if (pvt_owned_atm_cnts(taskid) .eq. 0 .or. &
                        pvt_owned_atm_cnts(taskid) + added_atms_cnt .le. &
                        target_cnt(taskid)) exit
                  end if
                end do
                if (.not. is_frc_ene_task(taskid)) taskid = ori_taskid
              end if

            end if

            pvt_owned_atm_cnts(taskid) = &
              pvt_owned_atm_cnts(taskid) + added_atms_cnt

            if (taskid .ne. mytaskid) then
              do i = mol_offset + 1, &
                mol_offset + mol_listcnt
                prf_id = gbl_mol_prfs_lists(i)
                prf_offset = gbl_prf_listdata(prf_id)%offset
                prf_listcnt = gbl_prf_listdata(prf_id)%cnt
                do j = prf_offset + 1, prf_offset + prf_listcnt
                  atm_id = gbl_prf_lists(j)
                  gbl_atm_owner_map(atm_id) = taskid
                end do
              end do
            else
              if (my_frag_mol_cnt .ne. 0) then
                if (gbl_my_frag_mol_lst(my_frag_mol_cnt) .ne. frag_mol_idx) then
                  my_frag_mol_cnt = my_frag_mol_cnt + 1
                  gbl_my_frag_mol_lst(my_frag_mol_cnt) = frag_mol_idx
                end if
              else
                my_frag_mol_cnt = my_frag_mol_cnt + 1
                gbl_my_frag_mol_lst(my_frag_mol_cnt) = frag_mol_idx
              end if
              do i = mol_offset + 1, mol_offset + mol_listcnt
                prf_id = gbl_mol_prfs_lists(i)
                prf_offset = gbl_prf_listdata(prf_id)%offset
                prf_listcnt = gbl_prf_listdata(prf_id)%cnt
                do j = prf_offset + 1, prf_offset + prf_listcnt
                  atm_id = gbl_prf_lists(j)
                  gbl_atm_owner_map(atm_id) = taskid
                  my_atm_cnt = my_atm_cnt + 1
                  gbl_my_atm_lst(my_atm_cnt) = atm_id
                end do
              end do
            end if
            gbl_molfrags(molfrag_idx)%owner = taskid
          end do

        end if
  
      end do

    else

      do mol_id = 1, gbl_mol_cnt
  
        added_atms_cnt = gbl_mol_atms_listdata(mol_id)%cnt
        mol_offset = gbl_mol_prfs_listdata(mol_id)%offset
        mol_listcnt = gbl_mol_prfs_listdata(mol_id)%cnt

        first_prf_id = gbl_mol_prfs_lists(mol_offset + 1)

        atm_id = gbl_prf_lists(gbl_prf_listdata(first_prf_id)%offset + 1)

        if (added_atms_cnt .le. max_unfrag_mol_atms) then

          img_id = atm_img_map(atm_id)
          taskid = img_owner_map(img_id)

          if (pvt_owned_atm_cnts(taskid) + added_atms_cnt .gt. &
              target_cnt(taskid)) then
  
            ! Just go ahead and assign to the current task unless there is
            ! already some assignment or the current task has a recip workload.
            ! Otherwise, look for a task that would not be overloaded.
            ! To support high_scaling algorithms, we save the original task
            ! assignment so we can fall back on it should the algorithm below
            ! pick a non-frc-ene task.
  
            if (pvt_owned_atm_cnts(taskid) .ne. 0 .or. &
                is_recip_task(taskid)) then
              
              ori_taskid = taskid ! original task assignment saved for furture
              do i = 1, numtasks - 1
                taskid = taskid + 1
                if (taskid .ge. numtasks) taskid = 0
                if (is_frc_ene_task(taskid)) then
                  if (pvt_owned_atm_cnts(taskid) + added_atms_cnt .le. &
                      target_cnt(taskid)) exit
                end if !(is_frc_ene_task(taskid))
              end do
              if (.not. is_frc_ene_task(taskid)) taskid = ori_taskid
  
              ! Okay, did we end up on a node where we overflowed anyway, and it
              ! also has a recip force workload?  If so, if there are any tasks
              ! without a recip force workload, overflow to them regardless of
              ! their target count.

              if (pvt_owned_atm_cnts(taskid) + added_atms_cnt .gt. &
                  target_cnt(taskid)) then
                if (is_recip_task(taskid)) then ! Ashraf
                  taskid = img_owner_map(img_id)
                  ori_taskid = taskid
                  do i = 1, numtasks
                    taskid = taskid + 1
                    if (taskid .ge. numtasks) taskid = 0
                    if (is_frc_ene_task(taskid)) then
                      if (.not. is_recip_task(taskid)) exit ! Ashraf
                    end if !( .not is_frc_ene_task(taskid))
                  end do
                  if (.not. is_frc_ene_task(taskid)) taskid = ori_taskid
                end if !(is_recip_task(taskid))
              end if !(pvt_owned_atm_cnts(taskid) + added_atms_cnt .gt. & ....
  
            end if !(pvt_owned_atm_cnts(taskid) .ne. 0 .or. & ....
          end if ! (pvt_owned_atm_cnts(taskid) + added_atms_cnt .gt. & ....
            
          pvt_owned_atm_cnts(taskid) = &
            pvt_owned_atm_cnts(taskid) + added_atms_cnt
  
          if (taskid .ne. mytaskid) then
            do i = mol_offset + 1, mol_offset + mol_listcnt
              prf_id = gbl_mol_prfs_lists(i)
              prf_offset = gbl_prf_listdata(prf_id)%offset
              prf_listcnt = gbl_prf_listdata(prf_id)%cnt
              do j = prf_offset + 1, prf_offset + prf_listcnt
                atm_id = gbl_prf_lists(j)
                gbl_atm_owner_map(atm_id) = taskid
              end do
            end do
          else
            my_mol_cnt = my_mol_cnt + 1
            gbl_my_mol_lst(my_mol_cnt) = mol_id
            do i = mol_offset + 1, mol_offset + mol_listcnt
              prf_id = gbl_mol_prfs_lists(i)
              prf_offset = gbl_prf_listdata(prf_id)%offset
              prf_listcnt = gbl_prf_listdata(prf_id)%cnt
              do j = prf_offset + 1, prf_offset + prf_listcnt
                atm_id = gbl_prf_lists(j)
                gbl_atm_owner_map(atm_id) = taskid
                my_atm_cnt = my_atm_cnt + 1
                gbl_my_atm_lst(my_atm_cnt) = atm_id
              end do
            end do
          end if

        else !(added_atms_cnt .gt. max_unfrag_mol_atms)

          frag_mol_idx = frag_mol_idx + 1

          do molfrag_idx = gbl_frag_mols(frag_mol_idx)%first_molfrag_idx, &
                           gbl_frag_mols(frag_mol_idx)%first_molfrag_idx + &
                           gbl_frag_mols(frag_mol_idx)%molfrag_cnt - 1

            mol_offset = gbl_molfrags(molfrag_idx)%offset
            mol_listcnt = gbl_molfrags(molfrag_idx)%cnt

            added_atms_cnt = 0

            do i = mol_offset + 1, mol_offset + mol_listcnt
              prf_id = gbl_mol_prfs_lists(i)
              added_atms_cnt = added_atms_cnt + gbl_prf_listdata(prf_id)%cnt
            end do

            first_prf_id = gbl_mol_prfs_lists(mol_offset + mol_listcnt)
            first_atm_id = &
              gbl_prf_lists(gbl_prf_listdata(first_prf_id)%offset+1)

            img_id = atm_img_map(first_atm_id)
            taskid = img_owner_map(img_id)

            if (pvt_owned_atm_cnts(taskid) + added_atms_cnt .gt. &
              target_cnt(taskid)) then
  
              ! Just go ahead and assign to the current task unless there is
              ! already some assignment or the current task has a recip
              ! workload.  Otherwise, look for a task that would not be
              ! overloaded.
              ! To support high_scaling algorithms, we save the original task
              ! assignment so we can fall back on it should the algorithm below
              ! pick a non-frc-ene task.
  

              if (pvt_owned_atm_cnts(taskid) .ne. 0 .or. &
                  is_recip_task(taskid)) then

                ori_taskid = taskid
                do i = 1, numtasks - 1
                  taskid = taskid + 1
                  if (taskid .ge. numtasks) taskid = 0
                  if (is_frc_ene_task(taskid)) then
                    if (pvt_owned_atm_cnts(taskid) + added_atms_cnt .le. &
                      target_cnt(taskid)) exit

                  end if
                end do
                if (.not. is_frc_ene_task(taskid)) taskid = ori_taskid

                ! Okay, did we end up on a node where we overflowed anyway,
                ! and it also has a recip force workload?  If so, if there are
                ! any tasks without a recip force workload, overflow to them
                ! regardless of their target count.

                if (pvt_owned_atm_cnts(taskid) + added_atms_cnt .gt. &
                    target_cnt(taskid)) then
                  if (is_recip_task(taskid)) then
                    taskid = img_owner_map(img_id)
                    ori_taskid = taskid
                    do i = 1, numtasks
                      taskid = taskid + 1
                      if (taskid .ge. numtasks) taskid = 0
                      if (is_frc_ene_task(taskid)) then
                        if (.not. is_recip_task(taskid)) exit
                      end if
                    end do
                    if (.not. is_frc_ene_task(taskid)) taskid = ori_taskid
                  end if
                end if

              end if
            end if

            pvt_owned_atm_cnts(taskid) = &
              pvt_owned_atm_cnts(taskid) + added_atms_cnt

            if (taskid .ne. mytaskid) then
              do i = mol_offset + 1, &
                mol_offset + mol_listcnt
                prf_id = gbl_mol_prfs_lists(i)
                prf_offset = gbl_prf_listdata(prf_id)%offset
                prf_listcnt = gbl_prf_listdata(prf_id)%cnt
                do j = prf_offset + 1, prf_offset + prf_listcnt
                  atm_id = gbl_prf_lists(j)
                  gbl_atm_owner_map(atm_id) = taskid
                end do
              end do
            else
              if (my_frag_mol_cnt .ne. 0) then
                if (gbl_my_frag_mol_lst(my_frag_mol_cnt) .ne. frag_mol_idx) then
                  my_frag_mol_cnt = my_frag_mol_cnt + 1
                  gbl_my_frag_mol_lst(my_frag_mol_cnt) = frag_mol_idx
                end if
              else
                my_frag_mol_cnt = my_frag_mol_cnt + 1
                gbl_my_frag_mol_lst(my_frag_mol_cnt) = frag_mol_idx
              end if
              do i = mol_offset + 1, mol_offset + mol_listcnt
                prf_id = gbl_mol_prfs_lists(i)
                prf_offset = gbl_prf_listdata(prf_id)%offset
                prf_listcnt = gbl_prf_listdata(prf_id)%cnt
                do j = prf_offset + 1, prf_offset + prf_listcnt
                  atm_id = gbl_prf_lists(j)
                  gbl_atm_owner_map(atm_id) = taskid
                  my_atm_cnt = my_atm_cnt + 1
                  gbl_my_atm_lst(my_atm_cnt) = atm_id
                end do
              end do
            end if
            gbl_molfrags(molfrag_idx)%owner = taskid
          end do

        end if
  
      end do !mol_id = 1, gbl_mol_cnt

    end if ! recip_numtasks .ne. numtasks

  else  ! not NTP ensemble...  Molecules not used...

    do prf_id = 1, gbl_prf_cnt
  
      prf_offset = gbl_prf_listdata(prf_id)%offset
      prf_listcnt = gbl_prf_listdata(prf_id)%cnt

      atm_id = gbl_prf_lists(prf_offset + 1)
      img_id = atm_img_map(atm_id)
      taskid = img_owner_map(img_id)
#if 0 /*the following did not help for NVE*/
     !the following is added by Ashraf for better offload performance
      if (offload_tasks(taskid)) then
         ! look for a task that does not offload.
         do i = 1, numtasks - 1
            taskid = taskid + 1
            if (taskid .ge. numtasks) taskid = 0
              if (is_frc_ene_task(taskid)) then
                if (.not. offload_tasks(taskid)) exit
              end if !(is_frc_ene_task(taskid))
       end do
     end if
#endif
      if (taskid .ne. mytaskid) then

        do i = prf_offset + 1, prf_offset + prf_listcnt
          atm_id = gbl_prf_lists(i)
          gbl_atm_owner_map(atm_id) = taskid
        end do

      else

        do i = prf_offset + 1, prf_offset + prf_listcnt
          atm_id = gbl_prf_lists(i)
          gbl_atm_owner_map(atm_id) = taskid
          my_atm_cnt = my_atm_cnt + 1
          gbl_my_atm_lst(my_atm_cnt) = atm_id
        end do

      end if

      pvt_owned_atm_cnts(taskid) = pvt_owned_atm_cnts(taskid) + prf_listcnt

    end do

  end if !assignment of direct force and bond atoms to tasks are done 

  ! Set up for distributed molecule handling if needed. ALL tasks have to call
  ! this stuff, even if they don't participate in processing a particular, or
  ! any molecule.

  if (frag_mol_cnt .gt. 0) then
    call destroy_communicators()        ! only does something if comm's exist...
    call create_communicators()
  end if

  pvt_atm_offsets(0) = 0

  do taskid = 0, numtasks - 1
    pvt_atm_offsets(taskid + 1) = pvt_atm_offsets(taskid) + &
                                   pvt_owned_atm_cnts(taskid)
  end do

!BEGIN DBG
! if (master) write(0,*)'| DBG: Node atom cnts =', pvt_owned_atm_cnts(:)

! do i = 1, atm_cnt
!   if (gbl_atm_owner_map(i) .lt. 0 .or. &
!       gbl_atm_owner_map(i) .ge. numtasks) &
!       write(0,*)'| DBG_ERR: atom ', i, 'assigned to bad task ', &
!                 gbl_atm_owner_map(i)
! end do
!END DBG

  do taskid = 0, numtasks - 1
    pvt_vec_rcvcnts(taskid) = 3 * pvt_owned_atm_cnts(taskid)
  end do

  pvt_vec_rcvcnts(numtasks) = 0

  pvt_vec_offsets(0) = 0

  do taskid = 0, numtasks - 1
    pvt_vec_offsets(taskid + 1) = pvt_vec_offsets(taskid) + &
                                   pvt_vec_rcvcnts(taskid)
  end do

  if (master) then
    if (loadbal_verbose .gt. 0 .or. write_log) then
        write(logfile, 10) &
          'Atom Distribution No. ', distrib_cnt, ' at run step ', &
          step_ctr, ':'
        write(logfile, 20) 'Count of atoms assigned to each task:'
        write(logfile, 30) (pvt_owned_atm_cnts(i), i = 0, numtasks - 1)
    end if
  end if

  ! Set up extra points per-process frame list if needed:

  if (numextra .gt. 0) then

    ! Claim frames for which you own the first ep:

    my_ep_frame_cnt = 0

    do i = 1, gbl_frame_cnt
      if (gbl_atm_owner_map(ep_frames(i)%extra_pnt(1)) .eq. mytaskid) then
        my_ep_frame_cnt = my_ep_frame_cnt + 1
        gbl_my_ep_frame_lst(my_ep_frame_cnt) = i
      end if
    end do

  end if

  ! There is parallel code that assumes that my_atm_lst() contains
  ! ascending atom id's.

  call sort_my_atm_lst(my_atm_cnt, gbl_my_atm_lst)

10 format(/, a, i5, a, i7, a /)
20 format('  ', a)
30 format('    ', 8i9)

  return

end subroutine divide_atoms
!*******************************************************************************
!
! Subroutine:  sort_my_atm_lst
!
! Description:  Sort the current task's atom list, leaving it in ascending
!               order.  This is an efficient linked list sort, with efficiency
!               based on prior partial ordering.  This is only intended for
!               use by divide_atoms().
!              
!*******************************************************************************

subroutine sort_my_atm_lst(my_atm_cnt, my_atm_lst)

  use gbl_datatypes_mod
  use parallel_dat_mod, only : master

  implicit none

! Formal arguments:

  integer               :: my_atm_cnt
  integer               :: my_atm_lst(my_atm_cnt)

! Local variables:

  type(atm_lst_rec)     :: atm_lst(my_atm_cnt)
  integer               :: atm_id
  integer               :: head_node
  integer               :: nxt_node
  integer               :: cur_node
  integer               :: i

  if (my_atm_cnt .le. 0) return

! BEGIN DBG
! if (master) then
! write(0, *) 'DBG: BEFORE sorting, my_atm_lst for master: ='
! write(0, '(10i8)') my_atm_lst(1:my_atm_cnt)
! end if
! END DBG

  head_node = 1
  atm_lst(1)%idx = my_atm_lst(1)
  atm_lst(1)%nxt = 0

  ! Sort atm_id's in descending order into a linked list.

  do i = 2, my_atm_cnt

    atm_id = my_atm_lst(i)
    atm_lst(i)%idx = atm_id
    if (atm_id .gt. atm_lst(head_node)%idx) then
      atm_lst(i)%nxt = head_node
      head_node = i
    else
      cur_node = head_node
      nxt_node = atm_lst(cur_node)%nxt
      do while (nxt_node .ne. 0)
        if (atm_id .gt. atm_lst(nxt_node)%idx) then
          atm_lst(cur_node)%nxt = i
          atm_lst(i)%nxt = nxt_node
          exit
        end if
        cur_node = nxt_node
        nxt_node = atm_lst(cur_node)%nxt
      end do
      if (nxt_node .eq. 0) then
        atm_lst(i)%nxt = 0
        atm_lst(cur_node)%nxt = i
      end if
    end if

  end do

  ! Now store them back in my_atm_lst reversed.

  cur_node = head_node
  do i = my_atm_cnt, 1, -1
    my_atm_lst(i) = atm_lst(cur_node)%idx
    cur_node = atm_lst(cur_node)%nxt
  end do

! BEGIN DBG
! if (master) then
! write(0, *) 'DBG: AFTER sorting, my_atm_lst for master: ='
! write(0, '(10i8)') my_atm_lst(1:my_atm_cnt)
! end if
! END DBG

  return

end subroutine sort_my_atm_lst
#endif /* MPI */

end module parallel_mod
