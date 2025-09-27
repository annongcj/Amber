#include "copyright.i"

#include "include_precision.i"
!*******************************************************************************
!
! Module:  nb_pairlist_mod
!
! Description: <TBS>
!
!*******************************************************************************

module nb_pairlist_mod


  implicit none

! Data that is not broadcast:

  double precision, save, allocatable   :: gbl_atm_saved_crd(:,:)
  double precision, save, allocatable   :: gbl_saved_imgcrd(:,:)
  integer, allocatable, save            :: gbl_ipairs(:)
  integer, allocatable, save            :: proc_ipairs(:)
  integer, save                         :: proc_ipairs_maxsize
  double precision, save                :: gbl_saved_box(3)
  integer, save                         :: ipairs_maxsize
  ! Pairlist definitions for TI
  integer, parameter    :: nb_list_null = 0
  integer, parameter    :: nb_list_common = 1
  integer, parameter    :: nb_list_linear_V0 = 2
  integer, parameter    :: nb_list_linear_V1 = 3
  integer, parameter    :: nb_list_sc_common_V0 = 4
  integer, parameter    :: nb_list_sc_common_V1 = 5
  integer, parameter    :: nb_list_sc_sc_V0 = 6
  integer, parameter    :: nb_list_sc_sc_V1 = 7
  integer, parameter    :: nb_list_cnt = 7
contains

!*******************************************************************************
!
! Subroutine:  alloc_nb_pairlist_mem
!
! Description: <TBS>
!
!*******************************************************************************

subroutine alloc_nb_pairlist_mem(atm_cnt, cutlist, num_ints, num_reals)

  use parallel_dat_mod
  use pmemd_lib_mod
  use processor_mod
  use mdin_ctrl_dat_mod, only : usemidpoint
  implicit none

! Formal arguments:

  integer, intent(in)           :: atm_cnt
  double precision, intent(in)  :: cutlist

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer                       :: alloc_failed
  integer, parameter            :: ipairs_size_min = 1000
  double precision, parameter   :: ipairs_size_coef = 0.225d0

  if (.not. usemidpoint) then

    allocate(gbl_atm_saved_crd(3, atm_cnt), &
             gbl_saved_imgcrd(3, atm_cnt), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    num_reals = num_reals + size(gbl_atm_saved_crd) + &
                size(gbl_saved_imgcrd)

    ! Allocate pairs list:

    ! The following ipairs_maxsize calc is a crude heuristic assuming that the
    ! number of nonbonded pairs roughly scales with the cutoff volume and
    ! the number of atoms in the system.  If MPI is running, the total is
    ! divided up among the processes.  The 2 or 3 is for the counters and flags at
    ! the head of each atom sublist.

    ipairs_maxsize = int((ipairs_size_coef * cutlist**3 + 3.d0) * dble(atm_cnt))

#ifdef MPI
    if (numtasks .gt. 0) ipairs_maxsize = ipairs_maxsize / numtasks
#endif /* MPI */

    if (ipairs_maxsize .lt. ipairs_size_min) ipairs_maxsize = ipairs_size_min

    allocate(gbl_ipairs(ipairs_maxsize), stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

  else
#ifdef MPI
    proc_ipairs_maxsize = int((ipairs_size_coef * cutlist**3 + 3.d0) * dble(atm_cnt))
    !proc_num_atms cannot be used, since the routine with proc_num_atms is called later
    if (numtasks .gt. 0) proc_ipairs_maxsize = proc_ipairs_maxsize / numtasks*2 !2 was added to be safe
    if (proc_ipairs_maxsize .lt. ipairs_size_min) proc_ipairs_maxsize = ipairs_size_min
    allocate(proc_ipairs(0:proc_ipairs_maxsize))
#endif /* MPI */
  endif

  return

end subroutine alloc_nb_pairlist_mem

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  bcast_nb_pairlist_dat
!
! Description: <TBS>
!
!*******************************************************************************

subroutine bcast_nb_pairlist_dat(atm_cnt, cutlist)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: atm_cnt
  double precision, intent(in)  :: cutlist

! Local variables:

  integer               :: num_ints, num_reals  ! returned values discarded

  if (.not. master) then
    num_ints = 0
    num_reals = 0
    call alloc_nb_pairlist_mem(atm_cnt, cutlist, num_ints, num_reals)
  end if

  ! The allocated data is not initialized from the master node.

  return

end subroutine bcast_nb_pairlist_dat
#endif

!*******************************************************************************
!
! Subroutine:   map_pairlist_imgs
!
! Description:  The img_crd() array were not numbered before. In this routing the
! img_crd() is numbered following the flat_cit()%lo and %hi (that is
! crd_idx_tbl()%lo and hi in cit.F90). The indexing of img_crd() is monotonous
! and same numbering work for all MPIs
!
!*******************************************************************************

#ifdef MPI
subroutine map_pairlist_imgs(atm_cnt, flat_cit, fraction, charge, iac, &
                             img_crd, img_qterm, img_iac, mapped_img_cnt, &
                             mapped_img_lst, img_atm_map)
#else
subroutine map_pairlist_imgs(atm_cnt, fraction, charge, iac, &
                             img_crd, img_qterm, img_iac, img_atm_map)
#endif

  use cit_mod
  use gbl_datatypes_mod
  use img_mod
  use parallel_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
#ifdef MPI
  type(cit_tbl_rec)     :: flat_cit(0 : cit_tbl_x_dim * &
                                        cit_tbl_y_dim * &
                                        cit_tbl_z_dim - 1)
  !flat_cit is the flattened version of crd_idx_tbl()%lo and hi, which has the
  !lo and hi of each bucket
#endif /* MPI */
  double precision      :: fraction(3, atm_cnt)
  double precision      :: charge(atm_cnt)
  integer               :: iac(atm_cnt)
#ifdef MPI
  integer               :: mapped_img_cnt
  integer               :: mapped_img_lst(*)
#endif /* MPI */
  double precision      :: img_crd(3, atm_cnt)
  double precision      :: img_qterm(atm_cnt)
  integer               :: img_iac(atm_cnt)
  integer               :: img_atm_map(atm_cnt)

! Local variables:

! "bkt" refers to a bucket or cell in the flat cit table.
! "bkts" refers to the bucket mapping tables.

  double precision      :: f1, f2, f3
  double precision      :: x_box, y_box, z_box
  double precision      :: ucell_stk(3, 3)
#ifdef MPI
  double precision      :: scale_fac_x, scale_fac_y, scale_fac_z
#endif /* MPI */
  integer               :: atm_idx, img_idx
  integer               :: is_orthog_stk
#ifdef MPI
  integer               :: lst_idx
  integer               :: i_bkt_x_idx, i_bkt_y_idx, i_bkt_z_idx
  integer               :: x_bkts_lo, x_bkts_hi
  integer               :: y_bkts_lo, y_bkts_hi
  integer               :: z_bkts_lo, z_bkts_hi
  integer               :: img_i_lo, atm_i_lo
  integer               :: i_bkt_lo, i_bkt_hi
  integer               :: i_bkt, cur_bkt, yz_bkt, z_bkt
  integer               :: x_bkts_idx, y_bkts_idx, z_bkts_idx
  integer               :: x_bkts(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: y_bkts(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: z_bkts(0 : cit_tbl_z_dim * 2 - 1)
#endif /* MPI */
#ifdef bkt_print
  integer               :: x_bkts_lo_min = 1000000, x_bkts_hi_max = 0
  integer               :: y_bkts_lo_min = 1000000, y_bkts_hi_max = 0
  integer               :: z_bkts_lo_min = 1000000, z_bkts_hi_max = 0
#endif

  x_box = pbc_box(1)
  y_box = pbc_box(2)
  z_box = pbc_box(3)

  ucell_stk(:,:) = ucell(:,:)

  is_orthog_stk = is_orthog

#ifdef MPI

  do lst_idx = 1, mapped_img_cnt
    img_idx = mapped_img_lst(lst_idx)
    img_iac(img_idx) = 0
  end do
  mapped_img_cnt = 0

  if (my_img_lo .gt. my_img_hi) return

  scale_fac_x = dble(cit_tbl_x_dim)
  scale_fac_y = dble(cit_tbl_y_dim)
  scale_fac_z = dble(cit_tbl_z_dim)

! Set up the bucket mapping arrays.  Redoing this is only necessary if the
! box dimensions change, but it is cheap and there is a locality benefit
! to having it on the stack.

  call setup_cit_tbl_bkts(x_bkts, y_bkts, z_bkts)

! Get the low and high flat cit bucket indexes for this task:

  call get_flat_cit_idx(my_img_lo, img_atm_map, fraction, i_bkt_lo)
  call get_flat_cit_idx(my_img_hi, img_atm_map, fraction, i_bkt_hi)

! Map all the images you will reference in building the pairlist.
! Note that "mapping" is different than claiming the images as "used".
! When mapped, this just means all the image's data structures are valid.


  do i_bkt = i_bkt_lo, i_bkt_hi

    img_i_lo = flat_cit(i_bkt)%img_lo
    if (img_i_lo .eq. 0) cycle

    ! Get the x,y,z 3d cit indexes to be able to get the bucket ranges
    ! to search:

    atm_i_lo = img_atm_map(img_i_lo)
    i_bkt_x_idx = int(fraction(1, atm_i_lo) * scale_fac_x)
    i_bkt_y_idx = int(fraction(2, atm_i_lo) * scale_fac_y)
    i_bkt_z_idx = int(fraction(3, atm_i_lo) * scale_fac_z)

    ! Determine the bucket ranges that need to be searched:

    x_bkts_lo = i_bkt_x_idx + cit_tbl_x_dim - cit_bkt_delta
    x_bkts_hi = i_bkt_x_idx + cit_tbl_x_dim + cit_bkt_delta
    y_bkts_lo = i_bkt_y_idx + cit_tbl_y_dim - cit_bkt_delta
    y_bkts_hi = i_bkt_y_idx + cit_tbl_y_dim + cit_bkt_delta
    z_bkts_lo = i_bkt_z_idx + cit_tbl_z_dim - cit_bkt_delta
    z_bkts_hi = i_bkt_z_idx + cit_tbl_z_dim

    z_loop: &
    do z_bkts_idx = z_bkts_lo, z_bkts_hi
      z_bkt = z_bkts(z_bkts_idx)
      y_loop: &
      do y_bkts_idx = y_bkts_lo, y_bkts_hi
        yz_bkt = z_bkt + y_bkts(y_bkts_idx)
        x_loop: &
        do x_bkts_idx = x_bkts_lo, x_bkts_hi
          cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
          img_i_lo = flat_cit(cur_bkt)%img_lo
          if (img_i_lo .eq. 0) cycle
          if (img_iac(img_i_lo) .ne. 0) cycle
          if (is_orthog_stk .ne. 0) then
            do img_idx = img_i_lo, flat_cit(cur_bkt)%img_hi
              atm_idx = img_atm_map(img_idx)
              img_iac(img_idx) = iac(atm_idx)
              img_crd(1, img_idx) = fraction(1, atm_idx) * x_box
              img_crd(2, img_idx) = fraction(2, atm_idx) * y_box
              img_crd(3, img_idx) = fraction(3, atm_idx) * z_box
              img_qterm(img_idx) = charge(atm_idx)
              mapped_img_cnt = mapped_img_cnt + 1
              mapped_img_lst(mapped_img_cnt) = img_idx
            end do
          else
            do img_idx = img_i_lo, flat_cit(cur_bkt)%img_hi
              atm_idx = img_atm_map(img_idx)
              f1 = fraction(1, atm_idx)
              f2 = fraction(2, atm_idx)
              f3 = fraction(3, atm_idx)
              ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0, so
              ! we can simplify the expression in this critical inner loop
              img_crd(1, img_idx) = f1 * ucell_stk(1, 1) + &
                                    f2 * ucell_stk(1, 2) + &
                                    f3 * ucell_stk(1, 3)
              img_crd(2, img_idx) = f2 * ucell_stk(2, 2) + &
                                    f3 * ucell_stk(2, 3)
              img_crd(3, img_idx) = f3 * ucell_stk(3, 3)
              img_qterm(img_idx) = charge(atm_idx)
              img_iac(img_idx) = iac(atm_idx)
              mapped_img_cnt = mapped_img_cnt + 1
              mapped_img_lst(mapped_img_cnt) = img_idx
            end do
          end if
          if (cur_bkt .eq. i_bkt) exit z_loop ! Done with cit buckets
                                             ! image i claims pairs from.
        end do x_loop
      end do y_loop
    end do z_loop

  end do

#else

! For the uniprocessor, just map it all...

  if (is_orthog_stk .ne. 0) then
    do img_idx = 1, atm_cnt
      atm_idx = img_atm_map(img_idx)
      img_crd(1, img_idx) = fraction(1, atm_idx) * x_box
      img_crd(2, img_idx) = fraction(2, atm_idx) * y_box
      img_crd(3, img_idx) = fraction(3, atm_idx) * z_box
      img_qterm(img_idx) = charge(atm_idx)
      img_iac(img_idx) = iac(atm_idx)
    end do
  else
    do img_idx = 1, atm_cnt
      atm_idx = img_atm_map(img_idx)
      f1 = fraction(1, atm_idx)
      f2 = fraction(2, atm_idx)
      f3 = fraction(3, atm_idx)
      ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0, so
      ! we can simplify the expression in this critical inner loop
      img_crd(1, img_idx) = f1 * ucell_stk(1, 1) + &
                            f2 * ucell_stk(1, 2) + &
                            f3 * ucell_stk(1, 3)
      img_crd(2, img_idx) = f2 * ucell_stk(2, 2) + &
                            f3 * ucell_stk(2, 3)
      img_crd(3, img_idx) = f3 * ucell_stk(3, 3)
      img_qterm(img_idx) = charge(atm_idx)
      img_iac(img_idx) = iac(atm_idx)
    end do
  end if

#endif /* MPI */

  return

end subroutine map_pairlist_imgs


!*******************************************************************************
!
! Subroutine:   get_nb_list_gp
!
! Description:  Main routine for pairlist building.
!
!*******************************************************************************
subroutine get_nb_list_gp(atm_cnt, flat_cit, img_crd, img_atm_map, &
#ifdef MPI
                       used_img_map, used_img_cnt, used_img_lst, &
#endif
                       ico, fraction, tranvec, atm_maskdata, atm_mask, &
                       atm_img_map, excl_img_flags, &
                       img_iac, igroup, ntypes, ibelly, &
                       es_cutoff, vdw_cutoff, skinnb, verbose, ifail)

  use cit_mod
  use gbl_datatypes_mod
  use img_mod
  use parallel_dat_mod
  use pbc_mod
  use ti_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  type(cit_tbl_rec)     :: flat_cit(0 : cit_tbl_x_dim * &
                                        cit_tbl_y_dim * &
                                        cit_tbl_z_dim - 1)
  double precision      :: img_crd(3, atm_cnt)

  integer               :: img_atm_map(atm_cnt)
#ifdef MPI
  integer(byte)         :: used_img_map(atm_cnt)
  integer               :: used_img_cnt
  integer               :: used_img_lst(*)
#endif
  integer               :: ico(*)
  double precision      :: fraction(3, atm_cnt)
  double precision      :: tranvec(1:3, 0:17)
  type(listdata_rec)    :: atm_maskdata(*)
  integer               :: atm_mask(*)
  integer               :: atm_img_map(atm_cnt)
  integer               :: excl_img_flags(atm_cnt)
  integer               :: img_iac(atm_cnt)
  integer               :: igroup(atm_cnt)
  integer               :: ntypes
  integer               :: ibelly
!the following double precison stays before the loop
  double precision      :: es_cutoff
  double precision      :: vdw_cutoff
  double precision      :: skinnb
  integer               :: verbose
  integer               :: ifail

! Local variables:
  double precision      :: scale_fac_x, scale_fac_y, scale_fac_z
  double precision      :: cutlist_sq
  double precision      :: x_i, y_i, z_i
  integer               :: atm_i, img_i
  integer               :: i, j, k
  integer               :: atm_mask_idx
  integer               :: num_packed
  integer               :: i_bkt_x_idx, i_bkt_y_idx, i_bkt_z_idx
  integer               :: x_bkts_lo, x_bkts_hi
  integer               :: y_bkts_lo, y_bkts_hi
  integer               :: z_bkts_lo, z_bkts_hi
  integer               :: i_bkt, saved_i_bkt_img_hi
  integer               :: img_i_lo, img_i_hi
  integer               :: i_bkt_lo, i_bkt_hi
  logical               :: common_tran
  logical               :: dont_skip_belly_pairs
  integer               :: iaci, ntypes_stk
  integer               :: ee_eval_cnt
  integer               :: full_eval_cnt
  integer               :: is_orthog_stk

! Variables use in the included files.  We would put this stuff in
! contained subroutines, but it really whacks performance for itanium/ifort.

  double precision      :: dx, dy, dz
  double precision      :: r_sq
  double precision      :: i_tranvec(1:3, 0:17)
  integer               :: cur_bkt, yz_bkt, z_bkt
  integer               :: cur_tran, yz_tran, z_tran
  integer               :: img_j
  integer               :: x_bkts_idx, y_bkts_idx, z_bkts_idx
  integer               :: enc_img

! Variables used for double cutoffs:

  double precision      :: es_cut_sq
  logical               :: cutoffs_equal

  integer               :: x_bkts(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: x_trans(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: y_bkts(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: y_trans(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: z_bkts(0 : cit_tbl_z_dim * 2 - 1)
  integer               :: z_trans(0 : cit_tbl_z_dim * 2 - 1)

  integer               :: total_pairs  ! DBG
  integer               :: img_j_ee_eval(atm_cnt)
  integer               :: img_j_full_eval(atm_cnt)

  num_packed = 0        ! For "regular" pairlist
  total_pairs = 0       ! DBG

  dont_skip_belly_pairs = ibelly .eq. 0 .and. ti_mode .eq. 0
  ntypes_stk = ntypes
  is_orthog_stk = is_orthog

  scale_fac_x = dble(cit_tbl_x_dim)
  scale_fac_y = dble(cit_tbl_y_dim)
  scale_fac_z = dble(cit_tbl_z_dim)

  es_cut_sq = (es_cutoff + skinnb) * (es_cutoff + skinnb)
  cutoffs_equal = (es_cutoff .eq. vdw_cutoff)
  cutlist_sq = (vdw_cutoff + skinnb) * (vdw_cutoff + skinnb)

  do i = 1, atm_cnt

    atm_i = i
    img_i = i
      ee_eval_cnt = 0
      full_eval_cnt = 0
      iaci = ntypes_stk * (img_iac(img_i) - 1)

      atm_i = img_atm_map(img_i)
 
      ! Mark excluded j images:

      atm_mask_idx = atm_maskdata(atm_i)%offset
      do k = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%cnt
        excl_img_flags(atm_img_map(atm_mask(k))) = 1
      end do

      ! These are imaged coordinates:
      x_i = img_crd(1, img_i)
      y_i = img_crd(2, img_i)
      z_i = img_crd(3, img_i)
      do img_j = i+1, atm_cnt
         dx = img_crd(1, img_j) - x_i
         dy = img_crd(2, img_j) - y_i
         dz = img_crd(3, img_j) - z_i
         if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
            if (excl_img_flags(img_j) .eq. 0) then
               if (ico(iaci + img_iac(img_j)) .eq. 0) then
                 ee_eval_cnt = ee_eval_cnt + 1
                 img_j_ee_eval(ee_eval_cnt) = img_j
               else
                 full_eval_cnt = full_eval_cnt + 1
                 img_j_full_eval(full_eval_cnt) = img_j
               end if
            end if /* excl_img_flag*/
         end if
      end do
      total_pairs = total_pairs + full_eval_cnt + ee_eval_cnt       ! DBG
      if (dont_skip_belly_pairs) then
        call pack_nb_list(ee_eval_cnt, img_j_ee_eval, &
                          full_eval_cnt, img_j_full_eval, &
                          gbl_ipairs, num_packed)
      else if (ti_mode .ne. 0) then
        call pack_nb_list_skip_ti_pairs(ee_eval_cnt, img_j_ee_eval, &
                                        full_eval_cnt, img_j_full_eval, &
                                        gbl_ipairs, img_atm_map, &
                                        num_packed)
      else
        call pack_nb_list_skip_belly_pairs(ee_eval_cnt, img_j_ee_eval, &
                                           full_eval_cnt, img_j_full_eval, &
                                           gbl_ipairs, igroup, img_atm_map, &
                                            num_packed)
      end if

      ! Clear excluded j images flags:

      do k = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%cnt
        excl_img_flags(atm_img_map(atm_mask(k))) = 0
      end do

  end do ! end of all cit buckets with i images owned by this task

! write(0, *) 'DBG: total_pairs, listtot = ', &
!             total_pairs, num_packed
! write(0, *) 'DBG: avg packing = ', &
!             dble(total_pairs)/dble(num_packed)

  if (verbose .gt. 0) then
    if (master) write(mdout, *) 'listtot = ', num_packed
  end if
  return

contains

#include "pack_nb_list_dflt.i"

end subroutine get_nb_list_gp

!*******************************************************************************
!
! Subroutine:   get_nb_list
!
! Description:  Main routine for pairlist building.
!
!*******************************************************************************
subroutine get_nb_list(atm_cnt, flat_cit, img_crd, img_atm_map, &
#ifdef MPI
                       used_img_map, used_img_cnt, used_img_lst, &
#endif
                       ico, fraction, tranvec, atm_maskdata, atm_mask, &
                       atm_img_map, excl_img_flags, &
                       img_iac, igroup, ntypes, ibelly, &
                       es_cutoff, vdw_cutoff, skinnb, verbose, ifail)

  use cit_mod
  use gbl_datatypes_mod
  use img_mod
  use parallel_dat_mod
  use pbc_mod
  use ti_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  type(cit_tbl_rec)     :: flat_cit(0 : cit_tbl_x_dim * &
                                        cit_tbl_y_dim * &
                                        cit_tbl_z_dim - 1)
  double precision      :: img_crd(3, atm_cnt)

  integer               :: img_atm_map(atm_cnt)
#ifdef MPI
  integer(byte)         :: used_img_map(atm_cnt)
  integer               :: used_img_cnt
  integer               :: used_img_lst(*)
#endif
  integer               :: ico(*)
  double precision      :: fraction(3, atm_cnt)
  double precision      :: tranvec(1:3, 0:17)
  type(listdata_rec)    :: atm_maskdata(*)
  integer               :: atm_mask(*)
  integer               :: atm_img_map(atm_cnt)
  integer               :: excl_img_flags(atm_cnt)
  integer               :: img_iac(atm_cnt)
  integer               :: igroup(atm_cnt)
  integer               :: ntypes
  integer               :: ibelly
!the following double precison stays before the loop
  double precision      :: es_cutoff
  double precision      :: vdw_cutoff
  double precision      :: skinnb
  integer               :: verbose
  integer               :: ifail

! Local variables:
  double precision      :: scale_fac_x, scale_fac_y, scale_fac_z
  double precision      :: cutlist_sq
  double precision      :: x_i, y_i, z_i
  integer               :: atm_i, img_i
  integer               :: i
  integer               :: atm_mask_idx
  integer               :: num_packed
  integer               :: i_bkt_x_idx, i_bkt_y_idx, i_bkt_z_idx
  integer               :: x_bkts_lo, x_bkts_hi
  integer               :: y_bkts_lo, y_bkts_hi
  integer               :: z_bkts_lo, z_bkts_hi
  integer               :: i_bkt, saved_i_bkt_img_hi
  integer               :: img_i_lo, img_i_hi
  integer               :: i_bkt_lo, i_bkt_hi
  logical               :: common_tran
  logical               :: dont_skip_belly_pairs
  integer               :: iaci, ntypes_stk
  integer               :: ee_eval_cnt
  integer               :: full_eval_cnt
  integer               :: is_orthog_stk

! Variables use in the included files.  We would put this stuff in
! contained subroutines, but it really whacks performance for itanium/ifort.

  double precision      :: dx, dy, dz
  double precision      :: r_sq
  double precision      :: i_tranvec(1:3, 0:17)
  integer               :: cur_bkt, yz_bkt, z_bkt
  integer               :: cur_tran, yz_tran, z_tran
  integer               :: img_j
  integer               :: x_bkts_idx, y_bkts_idx, z_bkts_idx
  integer               :: enc_img

! Variables used for double cutoffs:

  double precision      :: es_cut_sq
  logical               :: cutoffs_equal

  integer               :: x_bkts(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: x_trans(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: y_bkts(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: y_trans(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: z_bkts(0 : cit_tbl_z_dim * 2 - 1)
  integer               :: z_trans(0 : cit_tbl_z_dim * 2 - 1)

  integer               :: total_pairs  ! DBG
  integer               :: img_j_ee_eval(atm_cnt)
  integer               :: img_j_full_eval(atm_cnt)

! "bkt" refers to a bucket or cell in the flat cit table.
! "bkts" refers to the bucket mapping tables.
! "trans" refers to the translation mapping tables.

#ifdef MPI
  if (my_img_lo .gt. my_img_hi) return
#endif /* MPI */
  num_packed = 0        ! For "regular" pairlist
  total_pairs = 0       ! DBG

  dont_skip_belly_pairs = ibelly .eq. 0 .and. ti_mode .eq. 0
  ntypes_stk = ntypes
  is_orthog_stk = is_orthog

  scale_fac_x = dble(cit_tbl_x_dim)
  scale_fac_y = dble(cit_tbl_y_dim)
  scale_fac_z = dble(cit_tbl_z_dim)

  es_cut_sq = (es_cutoff + skinnb) * (es_cutoff + skinnb)
  cutoffs_equal = (es_cutoff .eq. vdw_cutoff)
  cutlist_sq = (vdw_cutoff + skinnb) * (vdw_cutoff + skinnb)

  ! Set up the bucket and translation index mapping arrays.  We really only need
  ! to do this each time if the cit table dimensions are changing, but there
  ! may be advantages to having only local variables anyway.

  call setup_cit_tbl_bkts(x_bkts, y_bkts, z_bkts, x_trans, y_trans, z_trans)

  ! Get the low and high flat cit bucket indexes for this task:

#ifdef MPI
  call get_flat_cit_idx(my_img_lo, img_atm_map, fraction, i_bkt_lo)
  call get_flat_cit_idx(my_img_hi, img_atm_map, fraction, i_bkt_hi)
#else
  i_bkt_lo = 0
  i_bkt_hi = cit_tbl_x_dim * cit_tbl_y_dim * cit_tbl_z_dim - 1
#endif /* MPI */

  ! Loop over the atoms you own, finding their nonbonded partners...
  ! img_i is the atom for which you are finding the partner atoms, img_j

  do i_bkt = i_bkt_lo, i_bkt_hi

    img_i_lo = flat_cit(i_bkt)%img_lo
    img_i_hi = flat_cit(i_bkt)%img_hi

    if (img_i_lo .eq. 0) cycle
    saved_i_bkt_img_hi = img_i_hi
#ifdef MPI
    if (img_i_lo .lt. my_img_lo) img_i_lo = my_img_lo
    if (img_i_hi .gt. my_img_hi) img_i_hi = my_img_hi
#endif /* MPI */

    ! Get the x,y,z 3d cit indexes to be able to get the bucket ranges
    ! to search:

    atm_i = img_atm_map(img_i_lo)
    i_bkt_x_idx = int(fraction(1, atm_i) * scale_fac_x)
    i_bkt_y_idx = int(fraction(2, atm_i) * scale_fac_y)
    i_bkt_z_idx = int(fraction(3, atm_i) * scale_fac_z)

    ! Determine the bucket ranges that need to be searched:

    x_bkts_lo = i_bkt_x_idx + cit_tbl_x_dim - cit_bkt_delta
    x_bkts_hi = i_bkt_x_idx + cit_tbl_x_dim + cit_bkt_delta
    y_bkts_lo = i_bkt_y_idx + cit_tbl_y_dim - cit_bkt_delta
    y_bkts_hi = i_bkt_y_idx + cit_tbl_y_dim + cit_bkt_delta
    z_bkts_lo = i_bkt_z_idx + cit_tbl_z_dim - cit_bkt_delta
    z_bkts_hi = i_bkt_z_idx + cit_tbl_z_dim

    common_tran = x_trans(x_bkts_lo) .eq. 1 .and. &
                  x_trans(x_bkts_hi) .eq. 1 .and. &
                  y_trans(y_bkts_lo) .eq. 3 .and. &
                  y_trans(y_bkts_hi) .eq. 3 .and. &
                  z_trans(z_bkts_lo) .eq. 9 .and. &
                  z_trans(z_bkts_hi) .eq. 9

    do img_i = img_i_lo, img_i_hi
        flat_cit(i_bkt)%img_hi = img_i - 1

      ee_eval_cnt = 0
      full_eval_cnt = 0
      iaci = ntypes_stk * (img_iac(img_i) - 1)

      atm_i = img_atm_map(img_i)

      ! Mark excluded j images:

      atm_mask_idx = atm_maskdata(atm_i)%offset
      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%cnt
        excl_img_flags(atm_img_map(atm_mask(i))) = 1
      end do

      ! These are imaged coordinates:
      x_i = img_crd(1, img_i)
      y_i = img_crd(2, img_i)
      z_i = img_crd(3, img_i)
      if (cutoffs_equal) then

        if (common_tran) then

          z_loop1: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            y_loop1: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              x_loop1: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) - x_i
                  dy = img_crd(2, img_j) - y_i
                  dz = img_crd(3, img_j) - z_i
                  if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval(ee_eval_cnt) = img_j
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = img_j
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
#endif /* MPI */
                    end if  ! excl_img_flag
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop1 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop1
            end do y_loop1
          end do z_loop1

        else    ! not common_tran

          do i = 0, 17
            i_tranvec(1, i) = tranvec(1, i) - x_i
            i_tranvec(2, i) = tranvec(2, i) - y_i
            i_tranvec(3, i) = tranvec(3, i) - z_i
          end do

          z_loop2: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            z_tran = z_trans(z_bkts_idx)
            y_loop2: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              yz_tran = z_tran + y_trans(y_bkts_idx)
              x_loop2: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                cur_tran = x_trans(x_bkts_idx) + yz_tran
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                  dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                  dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                  if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      enc_img = ior(img_j, ishft(cur_tran, 27))
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval(ee_eval_cnt) = enc_img
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = enc_img
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
#endif /* MPI */
                    end if  ! excl_img_flags
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop2 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop2
            end do y_loop2
          end do z_loop2

        end if

      else  ! cutoffs not equal
        if (common_tran) then

          z_loop3: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            y_loop3: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              x_loop3: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) - x_i
                  dy = img_crd(2, img_j) - y_i
                  dz = img_crd(3, img_j) - z_i
                  r_sq = dx * dx + dy * dy + dz * dz
                  if (r_sq .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        if (r_sq .lt. es_cut_sq) then
                          ee_eval_cnt = ee_eval_cnt + 1
                          img_j_ee_eval(ee_eval_cnt) = img_j
                        end if
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = img_j
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop3 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop3
            end do y_loop3
          end do z_loop3

        else  ! not common_tran

          do i = 0, 17
            i_tranvec(1, i) = tranvec(1, i) - x_i
            i_tranvec(2, i) = tranvec(2, i) - y_i
            i_tranvec(3, i) = tranvec(3, i) - z_i
          end do

          z_loop4: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            z_tran = z_trans(z_bkts_idx)
            y_loop4: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              yz_tran = z_tran + y_trans(y_bkts_idx)
              x_loop4: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                cur_tran = x_trans(x_bkts_idx) + yz_tran
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi

                  dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                  dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                  dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                  r_sq = dx * dx + dy * dy + dz * dz
                  if (r_sq .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        if (r_sq .lt. es_cut_sq) then
                          ee_eval_cnt = ee_eval_cnt + 1
                          img_j_ee_eval(ee_eval_cnt) = &
                            ior(img_j,ishft(cur_tran, 27))
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
#endif /* MPI */
                        end if
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = &
                          ior(img_j,ishft(cur_tran,27))
#ifdef MPI
                        if (used_img_map(img_j) .eq. 0) then
                          used_img_map(img_j) = 1
                          used_img_cnt = used_img_cnt + 1
                          used_img_lst(used_img_cnt) = img_j
                        end if
#endif /* MPI */
                      end if
                    end if
                  end if
                end do

                if (cur_bkt .eq. i_bkt) exit z_loop4 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop4
            end do y_loop4
          end do z_loop4
        end if
      end if

      total_pairs = total_pairs + full_eval_cnt + ee_eval_cnt       ! DBG
      if (dont_skip_belly_pairs) then
        call pack_nb_list(ee_eval_cnt, img_j_ee_eval, &
                          full_eval_cnt, img_j_full_eval, &
                          gbl_ipairs, num_packed)
      else if (ti_mode .ne. 0) then
        call pack_nb_list_skip_ti_pairs(ee_eval_cnt, img_j_ee_eval, &
                                        full_eval_cnt, img_j_full_eval, &
                                        gbl_ipairs, img_atm_map, &
                                        num_packed)
      else
        call pack_nb_list_skip_belly_pairs(ee_eval_cnt, img_j_ee_eval, &
                                           full_eval_cnt, img_j_full_eval, &
                                           gbl_ipairs, igroup, img_atm_map, &
                                            num_packed)
      end if

      ! Clear excluded j images flags:

      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%cnt
        excl_img_flags(atm_img_map(atm_mask(i))) = 0
      end do

      if (ifail .eq. 1) then
        flat_cit(i_bkt)%img_hi = saved_i_bkt_img_hi
        return
      end if
    end do ! end of images in i_bkt

    flat_cit(i_bkt)%img_hi = saved_i_bkt_img_hi

  end do ! end of all cit buckets with i images owned by this task

! write(0, *) 'DBG: total_pairs, listtot = ', &
!             total_pairs, num_packed
! write(0, *) 'DBG: avg packing = ', &
!             dble(total_pairs)/dble(num_packed)

  if (verbose .gt. 0) then
    if (master) write(mdout, *) 'listtot = ', num_packed
  end if
  return

contains

#include "pack_nb_list_dflt.i"

end subroutine get_nb_list

!*******************************************************************************
!
! Subroutine:   get_nb_ips_list
!
! Description:  Main routine for pairlist building.
!
!*******************************************************************************

subroutine get_nb_ips_list(atm_cnt, flat_cit, img_crd, img_atm_map, &
#ifdef MPI
                       used_img_map, used_img_cnt, used_img_lst, &
                       used_img_ips_map, used_img_ips_cnt, used_img_ips_lst, &
#endif
                       ico, fraction, tranvec, atm_maskdata, atm_mask, &
                       atm_img_map, excl_img_flags, &
                       img_iac, igroup, ntypes, ibelly, &
                       es_cutoff, vdw_cutoff, skinnb, verbose, ifail)

  use cit_mod
  use gbl_datatypes_mod
  use img_mod
  use parallel_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  type(cit_tbl_rec)     :: flat_cit(0 : cit_tbl_x_dim * &
                                        cit_tbl_y_dim * &
                                        cit_tbl_z_dim - 1)

  double precision      :: img_crd(3, atm_cnt)
  integer               :: img_atm_map(atm_cnt)
#ifdef MPI
  integer(byte)         :: used_img_map(atm_cnt)
  integer               :: used_img_cnt
  integer               :: used_img_lst(*)
  integer(byte)         :: used_img_ips_map(atm_cnt)
  integer               :: used_img_ips_cnt
  integer               :: used_img_ips_lst(*)
#endif
  integer               :: ico(*)
  double precision      :: fraction(3, atm_cnt)
  double precision      :: tranvec(1:3, 0:17)
  type(listdata_rec)    :: atm_maskdata(*)
  integer               :: atm_mask(*)
  integer               :: atm_img_map(atm_cnt)
  integer               :: excl_img_flags(atm_cnt)
  integer               :: img_iac(atm_cnt)
  integer               :: igroup(atm_cnt)
  integer               :: ntypes
  integer               :: ibelly
  double precision      :: es_cutoff
  double precision      :: vdw_cutoff
  double precision      :: skinnb
  integer               :: verbose
  integer               :: ifail

! Local variables:

  double precision      :: scale_fac_x, scale_fac_y, scale_fac_z
  double precision      :: cutlist_sq
  double precision      :: x_i, y_i, z_i
  integer               :: atm_i, img_i
  integer               :: i
  integer               :: atm_mask_idx
  integer               :: num_packed
  integer               :: i_bkt_x_idx, i_bkt_y_idx, i_bkt_z_idx
  integer               :: x_bkts_lo, x_bkts_hi
  integer               :: y_bkts_lo, y_bkts_hi
  integer               :: z_bkts_lo, z_bkts_hi
  integer               :: i_bkt, saved_i_bkt_img_hi
  integer               :: img_i_lo, img_i_hi
  integer               :: i_bkt_lo, i_bkt_hi
  logical               :: common_tran
  logical               :: dont_skip_belly_pairs
  integer               :: iaci, ntypes_stk
  integer               :: ee_eval_cnt
  integer               :: full_eval_cnt
  integer               :: is_orthog_stk

! Variables use in the included files.  We would put this stuff in
! contained subroutines, but it really whacks performance for itanium/ifort.

  double precision      :: dx, dy, dz
  integer               :: cur_bkt, yz_bkt, z_bkt
  integer               :: cur_tran, yz_tran, z_tran
  integer               :: img_j
  integer               :: x_bkts_idx, y_bkts_idx, z_bkts_idx
  double precision      :: r_sq
  integer               :: enc_img
  double precision      :: i_tranvec(1:3, 0:17)

! Variables used for double cutoffs:

  double precision      :: es_cut_sq
  logical               :: cutoffs_equal

  integer               :: x_bkts(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: x_trans(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: y_bkts(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: y_trans(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: z_bkts(0 : cit_tbl_z_dim * 2 - 1)
  integer               :: z_trans(0 : cit_tbl_z_dim * 2 - 1)

!  integer               :: total_pairs  ! DBG

  integer               :: img_j_ee_eval(atm_cnt)
  integer               :: img_j_full_eval(atm_cnt)

! "bkt" refers to a bucket or cell in the flat cit table.
! "bkts" refers to the bucket mapping tables.
! "trans" refers to the translation mapping tables.

#ifdef MPI
  if (my_img_lo .gt. my_img_hi) return
#endif /* MPI */
  num_packed = 0        ! For "regular" pairlist
  !total_pairs = 0       ! DBG

  dont_skip_belly_pairs = ibelly .eq. 0
  ntypes_stk = ntypes
  is_orthog_stk = is_orthog

  scale_fac_x = dble(cit_tbl_x_dim)
  scale_fac_y = dble(cit_tbl_y_dim)
  scale_fac_z = dble(cit_tbl_z_dim)

  es_cut_sq = (es_cutoff + skinnb) * (es_cutoff + skinnb)
  cutoffs_equal = (es_cutoff .eq. vdw_cutoff)
  cutlist_sq = (vdw_cutoff + skinnb) * (vdw_cutoff + skinnb)

  ! Set up the bucket and translation index mapping arrays.  We really only need
  ! to do this each time if the cit table dimensions are changing, but there
  ! may be advantages to having only local variables anyway.

  call setup_cit_tbl_bkts(x_bkts, y_bkts, z_bkts, x_trans, y_trans, z_trans)

  ! Get the low and high flat cit bucket indexes for this task:

#ifdef MPI
  call get_flat_cit_idx(my_img_lo, img_atm_map, fraction, i_bkt_lo)
  call get_flat_cit_idx(my_img_hi, img_atm_map, fraction, i_bkt_hi)
#else
  i_bkt_lo = 0
  i_bkt_hi = cit_tbl_x_dim * cit_tbl_y_dim * cit_tbl_z_dim - 1
#endif /* MPI */

  ! Loop over the atoms you own, finding their nonbonded partners...
  ! img_i is the atom for which you are finding the partner atoms, img_j

  do i_bkt = i_bkt_lo, i_bkt_hi

    img_i_lo = flat_cit(i_bkt)%img_lo
    img_i_hi = flat_cit(i_bkt)%img_hi

    if (img_i_lo .eq. 0) cycle
    saved_i_bkt_img_hi = img_i_hi
#ifdef MPI
    if (img_i_lo .lt. my_img_lo) img_i_lo = my_img_lo
    if (img_i_hi .gt. my_img_hi) img_i_hi = my_img_hi
#endif /* MPI */

    ! Get the x,y,z 3d cit indexes to be able to get the bucket ranges
    ! to search:

    atm_i = img_atm_map(img_i_lo)
    i_bkt_x_idx = int(fraction(1, atm_i) * scale_fac_x)
    i_bkt_y_idx = int(fraction(2, atm_i) * scale_fac_y)
    i_bkt_z_idx = int(fraction(3, atm_i) * scale_fac_z)

    ! Determine the bucket ranges that need to be searched:

    x_bkts_lo = i_bkt_x_idx + cit_tbl_x_dim - cit_bkt_delta
    x_bkts_hi = i_bkt_x_idx + cit_tbl_x_dim + cit_bkt_delta
    y_bkts_lo = i_bkt_y_idx + cit_tbl_y_dim - cit_bkt_delta
    y_bkts_hi = i_bkt_y_idx + cit_tbl_y_dim + cit_bkt_delta
    z_bkts_lo = i_bkt_z_idx + cit_tbl_z_dim - cit_bkt_delta
    z_bkts_hi = i_bkt_z_idx + cit_tbl_z_dim

    common_tran = x_trans(x_bkts_lo) .eq. 1 .and. &
                  x_trans(x_bkts_hi) .eq. 1 .and. &
                  y_trans(y_bkts_lo) .eq. 3 .and. &
                  y_trans(y_bkts_hi) .eq. 3 .and. &
                  z_trans(z_bkts_lo) .eq. 9 .and. &
                  z_trans(z_bkts_hi) .eq. 9

    do img_i = img_i_lo, img_i_hi

      flat_cit(i_bkt)%img_hi = img_i - 1

      ee_eval_cnt = 0
      full_eval_cnt = 0
      iaci = ntypes_stk * (img_iac(img_i) - 1)

      atm_i = img_atm_map(img_i)

      ! Mark excluded j images:

      atm_mask_idx = atm_maskdata(atm_i)%offset
      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%cnt
        excl_img_flags(atm_img_map(atm_mask(i))) = 1
      end do

      ! These are imaged coordinates:

      x_i = img_crd(1, img_i)
      y_i = img_crd(2, img_i)
      z_i = img_crd(3, img_i)

      if (cutoffs_equal) then

        if (common_tran) then

          z_loop1: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            y_loop1: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              x_loop1: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) - x_i
                  dy = img_crd(2, img_j) - y_i
                  dz = img_crd(3, img_j) - z_i
                  if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval(ee_eval_cnt) = img_j
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = img_j
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
                    else
                      if (used_img_ips_map(img_j) .eq. 0) then
                        used_img_ips_map(img_j) = 1
                        used_img_ips_cnt = used_img_ips_cnt + 1
                        used_img_ips_lst(used_img_ips_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop1 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop1
            end do y_loop1
          end do z_loop1

        else    ! not common_tran

          do i = 0, 17
            i_tranvec(1, i) = tranvec(1, i) - x_i
            i_tranvec(2, i) = tranvec(2, i) - y_i
            i_tranvec(3, i) = tranvec(3, i) - z_i
          end do

          z_loop2: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            z_tran = z_trans(z_bkts_idx)
            y_loop2: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              yz_tran = z_tran + y_trans(y_bkts_idx)
              x_loop2: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                cur_tran = x_trans(x_bkts_idx) + yz_tran
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                  dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                  dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                  if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      enc_img = ior(img_j, ishft(cur_tran, 27))
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval(ee_eval_cnt) = enc_img
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = enc_img
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
                    else
                      if (used_img_ips_map(img_j) .eq. 0) then
                        used_img_ips_map(img_j) = 1
                        used_img_ips_cnt = used_img_ips_cnt + 1
                        used_img_ips_lst(used_img_ips_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop2 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop2
            end do y_loop2
          end do z_loop2

        end if

      else  ! cutoffs not equal

        if (common_tran) then

          z_loop3: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            y_loop3: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              x_loop3: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) - x_i
                  dy = img_crd(2, img_j) - y_i
                  dz = img_crd(3, img_j) - z_i
                  r_sq = dx * dx + dy * dy + dz * dz
                  if (r_sq .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        if (r_sq .lt. es_cut_sq) then
                          ee_eval_cnt = ee_eval_cnt + 1
                          img_j_ee_eval(ee_eval_cnt) = img_j
                        end if
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = img_j
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
                    else
                      if (used_img_ips_map(img_j) .eq. 0) then
                        used_img_ips_map(img_j) = 1
                        used_img_ips_cnt = used_img_ips_cnt + 1
                        used_img_ips_lst(used_img_ips_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop3 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop3
            end do y_loop3
          end do z_loop3

        else  ! not common_tran

          do i = 0, 17
            i_tranvec(1, i) = tranvec(1, i) - x_i
            i_tranvec(2, i) = tranvec(2, i) - y_i
            i_tranvec(3, i) = tranvec(3, i) - z_i
          end do

          z_loop4: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            z_tran = z_trans(z_bkts_idx)
            y_loop4: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              yz_tran = z_tran + y_trans(y_bkts_idx)
              x_loop4: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                cur_tran = x_trans(x_bkts_idx) + yz_tran
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                  dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                  dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                  r_sq = dx * dx + dy * dy + dz * dz
                  if (r_sq .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        if (r_sq .lt. es_cut_sq) then
                          ee_eval_cnt = ee_eval_cnt + 1
                          img_j_ee_eval(ee_eval_cnt) = &
                            ior(img_j,ishft(cur_tran, 27))
#ifdef MPI
                          if (used_img_map(img_j) .eq. 0) then
                            used_img_map(img_j) = 1
                            used_img_cnt = used_img_cnt + 1
                            used_img_lst(used_img_cnt) = img_j
                          end if
#endif /* MPI */
                        end if
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = &
                          ior(img_j,ishft(cur_tran,27))
#ifdef MPI
                        if (used_img_map(img_j) .eq. 0) then
                          used_img_map(img_j) = 1
                          used_img_cnt = used_img_cnt + 1
                          used_img_lst(used_img_cnt) = img_j
                        end if
#endif /* MPI */
                      end if
#ifdef MPI
                    else
                      if (used_img_ips_map(img_j) .eq. 0) then
                        used_img_ips_map(img_j) = 1
                        used_img_ips_cnt = used_img_ips_cnt + 1
                        used_img_ips_lst(used_img_ips_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                end do

                if (cur_bkt .eq. i_bkt) exit z_loop4 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop4
            end do y_loop4
          end do z_loop4
        end if
      end if

      !total_pairs = total_pairs + full_eval_cnt + ee_eval_cnt       ! DBG

      if (dont_skip_belly_pairs) then
        call pack_nb_list(ee_eval_cnt, img_j_ee_eval, &
                          full_eval_cnt, img_j_full_eval, &
                          gbl_ipairs, num_packed)
      else
        call pack_nb_list_skip_belly_pairs(ee_eval_cnt, img_j_ee_eval, &
                                           full_eval_cnt, img_j_full_eval, &
                                           gbl_ipairs, igroup, img_atm_map, &
                                            num_packed)
      end if

      ! Clear excluded j images flags:

      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%cnt
        excl_img_flags(atm_img_map(atm_mask(i))) = 0
      end do

      if (ifail .eq. 1) then
        flat_cit(i_bkt)%img_hi = saved_i_bkt_img_hi
        return
      end if

    end do ! end of images in i_bkt

    flat_cit(i_bkt)%img_hi = saved_i_bkt_img_hi

  end do ! end of all cit buckets with i images owned by this task

! write(0, *) 'DBG: total_pairs, listtot = ', &
!             total_pairs, num_packed
! write(0, *) 'DBG: avg packing = ', &
!             dble(total_pairs)/dble(num_packed)

  if (verbose .gt. 0) then
    if (master) write(mdout, *) 'listtot = ', num_packed
  end if

  return

contains

#include "pack_nb_list_dflt.i"

end subroutine get_nb_ips_list

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  check_my_atom_movement
!
! Description:  This routine checks if any atom on this task's atom list has
!               moved more than half the skin distance.  This is intended for
!               use only with parallel MD.
!
!*******************************************************************************

subroutine check_my_atom_movement(crd, saved_crd, my_atm_lst, skinnb, ntp, &
                                  new_list)

  use parallel_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  double precision      :: crd(3, *)
  double precision      :: saved_crd(3, *)
  integer               :: my_atm_lst(*)
  double precision      :: skinnb
  integer               :: ntp
  logical               :: new_list

! Local variables:

  integer               :: atm_lst_idx
  double precision      :: dx, dy, dz
  double precision      :: box_dx, box_dy, box_dz
  double precision      :: distance_limit
  integer               :: n

! If anybody moves more than half the nbskin added to cutoff in generating the
! pairlist, update the pairlist.

  distance_limit = 0.25d0 * skinnb * skinnb

  if (ntp .eq. 0) then  ! Constant volume

    do atm_lst_idx = 1, my_atm_cnt
      n = my_atm_lst(atm_lst_idx)

      dx = crd(1, n) - saved_crd(1, n)
      dy = crd(2, n) - saved_crd(2, n)
      dz = crd(3, n) - saved_crd(3, n)

      if (dx * dx + dy * dy + dz * dz .ge. distance_limit) then
        new_list = .true.
        return
      end if

    end do

  else  ! Constant pressure scaling.

    box_dx = pbc_box(1) / gbl_saved_box(1)
    box_dy = pbc_box(2) / gbl_saved_box(2)
    box_dz = pbc_box(3) / gbl_saved_box(3)

    do atm_lst_idx = 1, my_atm_cnt
      n = my_atm_lst(atm_lst_idx)

      dx = crd(1, n) - saved_crd(1, n) * box_dx
      dy = crd(2, n) - saved_crd(2, n) * box_dy
      dz = crd(3, n) - saved_crd(3, n) * box_dz

      if (dx * dx + dy * dy + dz * dz .ge. distance_limit) then
        new_list = .true.
        return
      end if

    end do

  end if

  new_list = .false.

  return

end subroutine check_my_atom_movement
#endif /* MPI */

!*******************************************************************************
!
! Subroutine:  check_all_atom_movement
!
! Description:  This routine checks if any atom has moved more than half the
!               skin distance.  This is intended for all minimizations and
!               for single processor MD.
!*******************************************************************************

subroutine check_all_atom_movement(atm_cnt, crd, saved_crd, skinnb, ntp, &
                                   new_list)

  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: crd(3, *)
  double precision      :: saved_crd(3, *)
  double precision      :: skinnb
  integer               :: ntp
  logical               :: new_list

! Local variables:

  double precision      :: dx, dy, dz
  double precision      :: box_dx, box_dy, box_dz
  double precision      :: distance_limit
  integer               :: n

! If anybody moves more than half the nbskin added to cutoff in generating the
! pairlist, update the pairlist.

  distance_limit = 0.25d0 * skinnb * skinnb

  if (ntp .eq. 0) then  ! Constant volume

    do n = 1, atm_cnt

      dx = crd(1, n) - saved_crd(1, n)
      dy = crd(2, n) - saved_crd(2, n)
      dz = crd(3, n) - saved_crd(3, n)

      if (dx * dx + dy * dy + dz * dz .ge. distance_limit) then
        new_list = .true.
        return
      end if

    end do

  else  ! Constant pressure scaling.

    box_dx = pbc_box(1) / gbl_saved_box(1)
    box_dy = pbc_box(2) / gbl_saved_box(2)
    box_dz = pbc_box(3) / gbl_saved_box(3)

    do n = 1, atm_cnt

      dx = crd(1, n) - saved_crd(1, n) * box_dx
      dy = crd(2, n) - saved_crd(2, n) * box_dy
      dz = crd(3, n) - saved_crd(3, n) * box_dz

      if (dx * dx + dy * dy + dz * dz .ge. distance_limit) then
        new_list = .true.
        return
      end if

    end do

  end if

  new_list = .false.
  return

end subroutine check_all_atom_movement

!*******************************************************************************
!
! Subroutine:  save_all_atom_crds
!
! Description:  This is needed for the skin test for buffered pairlists.
!*******************************************************************************

subroutine save_all_atom_crds(atm_cnt, crd, saved_crd)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: crd(3, atm_cnt)
  double precision      :: saved_crd(3, atm_cnt)

! Local variables:

  integer               :: i

  do i = 1, atm_cnt
    saved_crd(1, i) = crd(1, i)
    saved_crd(2, i) = crd(2, i)
    saved_crd(3, i) = crd(3, i)
  end do

  return

end subroutine save_all_atom_crds

!*******************************************************************************
!
! Subroutine:  save_imgcrds
!
! Description:  <TBS>
!
!*******************************************************************************

#ifdef MPI
subroutine save_imgcrds(img_crd, saved_imgcrd, used_img_cnt, used_img_lst)
#else /*MPI*/
subroutine save_imgcrds(img_cnt, img_crd, saved_imgcrd)
#endif

  use img_mod

  implicit none

! Formal arguments:

#ifdef MPI
  double precision      :: img_crd(3, *)
  double precision      :: saved_imgcrd(3, *)
  integer               :: used_img_cnt
  integer               :: used_img_lst(*)
#else
  integer               :: img_cnt
  double precision      :: img_crd(3, *)
  double precision      :: saved_imgcrd(3, *)
#endif /* MPI */

! Local variables:

#ifdef MPI
  integer               :: lst_idx
#endif /* MPI */
  integer               :: img_idx

#ifdef MPI
  do lst_idx = 1, used_img_cnt
    img_idx = used_img_lst(lst_idx)
    saved_imgcrd(:, img_idx) = img_crd(:, img_idx)
  end do
#else
  do img_idx = 1, img_cnt
    saved_imgcrd(:,img_idx) = img_crd(:,img_idx)
  end do
#endif
  return

end subroutine save_imgcrds

!*******************************************************************************
!
! Subroutine:  adjust_imgcrds
!
! Description:  <TBS>
!
!*******************************************************************************

#ifdef MPI
subroutine adjust_imgcrds(img_cnt, img_crd, img_atm_map, &
                          used_img_cnt, used_img_lst, &
                          saved_imgcrd, crd, saved_crd, ntp)
#else /*MPI*/
subroutine adjust_imgcrds(img_cnt, img_crd, img_atm_map, &
                          saved_imgcrd, crd, saved_crd, ntp)
#endif /* MPI */

  use gbl_datatypes_mod
  use img_mod
  use pbc_mod

  implicit none

! Formal arguments:

#ifdef MPI
  integer               :: img_cnt
  double precision      :: img_crd(3, *)
  integer               :: img_atm_map(*)
  integer               :: used_img_cnt
  integer               :: used_img_lst(*)
  double precision      :: saved_imgcrd(3, *)
  double precision      :: crd(3, *)
  double precision      :: saved_crd(3, *)
  integer               :: ntp
#else
  integer               :: img_cnt
  double precision      :: img_crd(3, *)
  integer               :: img_atm_map(*)
  double precision      :: saved_imgcrd(3, *)
  double precision      :: crd(3, *)
  double precision      :: saved_crd(3, *)
  integer               :: ntp
#endif /* MPI */

! Local variables:

  integer               :: atm_idx
  integer               :: img_idx
  double precision      :: box_del(3)
#ifdef MPI
  integer               :: lst_idx
#endif /* MPI */

#ifdef MPI
  if (ntp .eq. 0) then  ! Constant volume

      do lst_idx = 1, used_img_cnt
        img_idx = used_img_lst(lst_idx)
        atm_idx = img_atm_map(img_idx)
        img_crd(:,img_idx) = crd(:,atm_idx) + saved_imgcrd(:,img_idx) - &
                             saved_crd(:,atm_idx)
      end do

  else  ! Constant pressure scaling.

    box_del(:) = pbc_box(:) / gbl_saved_box(:)

      do lst_idx = 1, used_img_cnt
        img_idx = used_img_lst(lst_idx)
        atm_idx = img_atm_map(img_idx)
        img_crd(:,img_idx) = crd(:,atm_idx) + (saved_imgcrd(:,img_idx) - &
                             saved_crd(:,atm_idx)) * box_del(:)
      end do

  end if

#else  /*MPI*/

  if (ntp .eq. 0) then  ! Constant volume

    do img_idx = 1, img_cnt
        atm_idx = img_atm_map(img_idx)
        img_crd(:,img_idx) = crd(:,atm_idx) + saved_imgcrd(:,img_idx) - &
                             saved_crd(:,atm_idx)
    end do

  else  ! Constant pressure scaling.

    box_del(:) = pbc_box(:) / gbl_saved_box(:)

    do img_idx = 1, img_cnt
      atm_idx = img_atm_map(img_idx)
      img_crd(:,img_idx) = crd(:,atm_idx) + (saved_imgcrd(:,img_idx) - &
                           saved_crd(:,atm_idx)) * box_del(:)
    end do

  end if

#endif /* MPI */

  return

end subroutine adjust_imgcrds

#ifdef MPI
!*******************************************************************************
!
! Subroutine:   get_nb_list_midpoint
!
! Description:  Main routine for pairlist building using midpoint method
!
!*******************************************************************************
subroutine get_nb_list_midpoint( atm_cnt, atm_crd, &
                        ico, ntypes, ibelly, &
                       es_cutoff, vdw_cutoff, skinnb, verbose, ifail, proc_neighborbuckets, atm_maskdata,atm_mask)
!need excluding list (excl_img_flags) by using proc_atm_mask etc
!es cut-off and vdw cut-off are assumed same, do we need to different cutoff?
!how do we handle ntypes and ibelly -- shall we reuse?
! these works for single MPI, not tested for multi-MPI

  use cit_mod
  use gbl_datatypes_mod
  use img_mod
  use parallel_dat_mod
  use pbc_mod
  use ti_mod
  use nb_exclusions_mod
! use processor_mod , only: int_bkt_minx, int_bkt_miny, int_bkt_minz, &
!                          int_bkt_maxx, int_bkt_maxy, int_bkt_maxz,proc_iac,  &
!                         unflatten, int_bkt_dimx, &
!                         int_bkt_dimy, int_bkt_dimz, int_unflatten, &
!                         int_mpi_map, proc_bkt_minx, proc_bkt_miny, proc_bkt_minz
  use processor_mod
  use mdin_ctrl_dat_mod, only : ntp
  use bonds_midpoint_mod


  implicit none
!Formal parameter
  integer               :: atm_cnt
  double precision      :: atm_crd(:,:) ! proc version
  integer               :: ico(:) ! need to use proc version
  integer               :: ntypes
  integer               :: ibelly
!the following double precison stays before the loop
  double precision      :: es_cutoff
  double precision      :: vdw_cutoff
  double precision      :: skinnb
  integer               :: verbose
  integer               :: ifail
  integer               :: proc_neighborbuckets(:,:) ! global bkt index
  type(listdata_rec)    :: atm_maskdata(*)
  integer               :: atm_mask(*)
!local variables

 !DIR$ ATTRIBUTES ALIGN : 64 :: atm_j_ee_eval, atm_j_full_eval,x_crd_j,y_crd_j,z_crd_j, atm_j_lst, proc_iac_j_lst
!x_crd_j(1024),y_crd_j(1024),z_crd_j(1024)
  integer               :: atm_j_eval(8192)! only es neighbors, probably 1000 is enough
  integer               :: atm_j_ee_eval(4096)! only es neighbors, probably 1000 is enough
  integer               :: atm_j_full_eval(4096)! es and vdw neighbors,  probably 1000 is enough
  integer               :: numElements ! total neighboring buckets for each central bucket
  integer               :: bkt,i,  first_ind, ind, numAtoms
  integer               :: cur_bkt, nl_bkt, atm_i, atm_j, nxt_atm_j !me and neighbor
  integer               :: x, y, z, xj, yj, zj !3d co-ordinate of bucket
  integer               :: enc_atm_j ! wrapped version of atm_j
  integer               :: iaci, ntypes_stk
  integer               :: neighbor_bkt
  integer               :: num_packed,total_pairs=0
  integer               :: ext_bkt(3), int_bkt(3)
  integer               :: intflatid
  integer               :: bkt_i_list(500), bkt_i_size, bkt_i_idx
  !debug starts
  integer               :: count1,total_count=0, count2=0, x_bkt, y_bkt,z_bkt,flat_bkt
  !debug ends
  integer               :: atm_mask_idx
 ! integer               :: excl_img_flags(atm_cnt*numtasks)

  integer               :: bondi_start, bondi_end ! atm_i star/end idx in atm_mask
  integer               :: atm_i_idx, atm_j_idx, bkt_cnt
  integer               :: dummy_common_tran, ee_eval_cnt, full_eval_cnt, eval_cnt
  integer               :: bins(50),nbins
  integer               :: atm_j_cnt, atm_j_lst(8192), proc_iac_j_lst(8192)
 ! HARDCODED to the ICO range of Cellulose benchmark
  integer               :: j, bkt_i_enabled
  integer               :: ico_imask
  integer               :: swap, cnt, idx


  !DIR$ ATTRIBUTES ALIGN : 64 :: atm_ee_eval, atm_full_eval, bkt_crd_lst
  integer,allocatable,save              :: excl_img_flags2(:) !it can be byte size

  integer,allocatable,save              :: atm_ee_eval(:,:)
  integer,allocatable,save              :: atm_full_eval(:,:)
  integer,allocatable,save              :: atms_ee_eval_cnt(:)
  integer,allocatable,save              :: atms_full_eval_cnt(:)
  integer                               :: proc_atm_cnt_maxsize=0 !to dynamically expand the sparse neigbor list as needed

#if 0
! #ifdef pmemd_SPDP
  pme_float      :: two=2.0
  pme_float      :: x_i, y_i, z_i!3d atm co-ordinate
  pme_float      :: dx, dy, dz !distance
  pme_float      :: cutlist_sq
  pme_float      :: midpoint(3)
  pme_float                       :: x_crd_j(8192),y_crd_j(8192),z_crd_j(8192)
  pme_float,allocatable,save     :: bkt_crd_lst(:,:,:) ! SoA crd atoms in each bucket
#else
  double precision      :: two=2.d0
  double precision      :: x_i, y_i, z_i!3d atm co-ordinate
  double precision      :: dx, dy, dz !distance
  double precision      :: cutlist_sq
  double precision      :: midpoint(3)
  double precision                      :: x_crd_j(8192),y_crd_j(8192),z_crd_j(8192)
  double precision,allocatable,save     :: bkt_crd_lst(:,:,:) ! SoA crd atoms in each bucket
#endif
  integer,allocatable                   :: ico_mask(:)
  integer                               :: all_num_atms=0
  integer                               :: tot_packed, all_pairs, mytaskid2
  integer                               :: max_atms, min_atms, mean_atms

  bkt_cnt = (int_bkt_dimx+2)*(int_bkt_dimy+2)*(int_bkt_dimz+2)

  if(.not. allocated(excl_img_flags2)) then
    call MPI_AllReduce(proc_num_atms, all_num_atms, 1, MPI_INTEGER, MPI_SUM,pmemd_comm, err_code_mpi)

    allocate(excl_img_flags2(all_num_atms))
    excl_img_flags2 = 0
  endif

  if(.not. allocated(bkt_crd_lst)) allocate(bkt_crd_lst(256, 3, (int_bkt_dimx+2)*(int_bkt_dimy+2)*(int_bkt_dimz+2)))

  if(.not. allocated(ico_mask)) then
    allocate(ico_mask(ntypes))  ! assuming ico range is ntypes_stk
    ! test if ntypes ls larger than the bits stize of the mask
    if(bit_size(ico_mask(1)) .lt. ntypes) then
      print*, "ERROR: number of types is larger than the used data type"
      call exit(-1)
    endif
    ico_mask=0
    do i = 0,ntypes-1
      do j = 0,ntypes-1
        ind = i*ntypes + j + 1
        ico_imask=0
        if(ico(ind) .eq. 0) ico_imask = ibset(ico_imask,j)
        ico_mask(i+1) = ior(ico_mask(i+1), ico_imask)
      enddo
!if(mytaskid ==0) then
!  print*, i, ico(i*ntypes+1:i*ntypes+ntypes)
!  write(*,'(b32.32)') ico_mask(i+1)
!endif
    enddo
  end if


  dummy_common_tran = 1 ! in case we use common_tran, we modify
  num_packed = 0        ! For "regular" pairlist
  total_pairs=0
  !proc_corner_neighborbuckets() contains the list of neighbor buckets for each bucket
  !row: list of all buckets in the sub-domain; col:list of neighbor buckets
  !proc_corner_neighborbuckets(i,:)=list of neighbor buckets for i bucket
  ntypes_stk=ntypes

  i = 2 !bucket is half the cut_off+skin, if full cut-off+skin then i = 1
  numElements = (i*((2*i+1)**2))+(i*(2*i+1))+i+1 !63 for bucket size = half cut-off. we can move the calculatino
  cutlist_sq = (vdw_cutoff + skinnb) * (vdw_cutoff + skinnb)

  do bkt = 1,bkt_cnt
    do atm_i_idx = 1,bkt_atm_cnt(bkt)
      atm_i = bkt_atm_lst(atm_i_idx,bkt)
      bkt_atm_lst(atm_i_idx,bkt) = atm_i
      bkt_crd_lst(atm_i_idx,1, bkt) = atm_crd(1,atm_i)
      bkt_crd_lst(atm_i_idx,2, bkt) = atm_crd(2,atm_i)
      bkt_crd_lst(atm_i_idx,3, bkt) = atm_crd(3,atm_i)
    enddo
  enddo



  do bkt = 1,bkt_cnt
!  do x=0,int_bkt_dimx+1
!  do y=0,int_bkt_dimy+1
!  do z=0,int_bkt_dimz+1
!    bkt=int_flatten(x,y,z)

    ! Copy the neighbors coordinates from AoS to SoA
    bkt_i_enabled=1
    atm_j_cnt=0
    do neighbor_bkt = 1,  proc_corner_neighborbuckets_cnt(bkt)
      cur_bkt = proc_neighborbuckets(neighbor_bkt,bkt)

!      !$omp parallel do private(atm_j, idx)
      do atm_j_idx = 1, bkt_atm_cnt(cur_bkt)
        idx = atm_j_cnt+atm_j_idx
        atm_j = bkt_atm_lst(atm_j_idx, cur_bkt)
        atm_j_lst(idx) = atm_j
        x_crd_j(idx) = bkt_crd_lst(atm_j_idx,1, cur_bkt)
        y_crd_j(idx) = bkt_crd_lst(atm_j_idx,2, cur_bkt)
        z_crd_j(idx) = bkt_crd_lst(atm_j_idx,3, cur_bkt)
        proc_iac_j_lst(idx) = proc_iac(atm_j)-1
        excl_img_flags2(proc_atm_to_full_list(atm_j)) = idx
      end do
      atm_j_cnt =atm_j_cnt + bkt_atm_cnt(cur_bkt)
    end do
    if (bkt_i_enabled .eq. 1)  atm_j_cnt = atm_j_cnt-bkt_atm_cnt(bkt) ! to allow adding the current bkt atms dynamically

    ! build the neighbor list of the atoms in the current bucket
    do atm_i_idx = 1, bkt_atm_cnt(bkt)
      atm_i = bkt_atm_lst(atm_i_idx,bkt)
      !get the co-ordinates
      x_i = bkt_crd_lst(atm_i_idx,1, bkt)
      y_i = bkt_crd_lst(atm_i_idx,2, bkt)
      z_i = bkt_crd_lst(atm_i_idx,3, bkt)

      iaci = ntypes_stk * (proc_iac(atm_i) - 1)
      ico_imask = ico_mask(proc_iac(atm_i))

!if(atm_i .eq. 1) then
!  write(15+mytaskid,*)atm_i,bkt
!end if
      ! get the current atm bond mask
      atm_mask_idx = atm_maskdata(proc_atm_to_full_list(atm_i))%offset
      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(proc_atm_to_full_list(atm_i))%cnt
        if(excl_img_flags2(atm_mask(i)) .ne. 0) then
          atm_j_excl_list(excl_img_flags2(atm_mask(i))) = 1
        end if
      end do

      ! reset neigbor list counter
      eval_cnt = 0

      ! loop over neighbor buckets
#if defined(MPI) && !defined(CUDA)
 !DIR$ ASSUME_ALIGNED x_crd_j:64,y_crd_j:64,z_crd_j:64,atm_j_lst:64
 !DIR$ VECTOR ALWAYS
! !DIR$ UNROLL=4
! !DIR$ LOOP COUNT AVG=1000
#endif
      do atm_j_idx = 1, atm_j_cnt + (atm_i_idx-1)*bkt_i_enabled
        atm_j = atm_j_lst(atm_j_idx)
        !both real space (inside box) j and j beyond box, both are treated same way
        dx = x_crd_j(atm_j_idx) - x_i !translated version of atm_j -x_i
        dy = y_crd_j(atm_j_idx) - y_i !translated version of atm_j -y_i
        dz = z_crd_j(atm_j_idx) - z_i !translated version of atm_j -z_i
        !Grabs atom midpoint
        midpoint(1)=(x_i+x_crd_j(atm_j_idx))/two !dble(2)
        midpoint(2)=(y_i+y_crd_j(atm_j_idx))/two !dble(2)
        midpoint(3)=(z_i+z_crd_j(atm_j_idx))/two !dble(2)
        !enc_atm_j = ior(atm_j, ishft(13, 27)) ! kind of wrap, will unwrap in pairs_cal_f

        !!!we can get rid of enc_atm_j scalar and use atm_j all through

        enc_atm_j = atm_j
        if( btest(ico_imask, proc_iac_j_lst(atm_j_idx)) ) enc_atm_j=-enc_atm_j
        if(atm_j_excl_list(atm_j_idx) .ne. 0) enc_atm_j=0
        if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
          if(     midpoint(1) .gt. proc_min_x_crd .and. midpoint(1) .lt. proc_max_x_crd &
            .and. midpoint(2) .gt. proc_min_y_crd .and. midpoint(2) .lt. proc_max_y_crd &
            .and. midpoint(3) .gt. proc_min_z_crd .and. midpoint(3) .lt. proc_max_z_crd) then
!            if(atm_i .eq. atm_j) then
!              write(15+mytaskid,*)"Match",atm_i,atm_j
!              cycle
!            end if
            eval_cnt = eval_cnt + 1
            atm_j_eval(eval_cnt) = enc_atm_j
          end if  ! int_mpi_map
        end if  ! cutoff check
      end do  ! end atm_j loop

      ! Temporary transformation:
      ee_eval_cnt = 0
      full_eval_cnt = 0
      do i = 1,eval_cnt
        if(atm_j_eval(i) < 0) then
          ee_eval_cnt = ee_eval_cnt + 1
          atm_j_ee_eval(ee_eval_cnt) = -atm_j_eval(i)
        end if
        if(atm_j_eval(i) > 0) then
          full_eval_cnt = full_eval_cnt + 1
          atm_j_full_eval(full_eval_cnt) = atm_j_eval(i)
        endif

      end do

!#ifdef DBG_ARRAYS
!  if(dbg_arrays_enabled) then
!    write(DBG_ARRAYS+mytaskid,'(I0.2," ",I0.10,":",I)') 2, proc_atm_to_full_list(atm_i), sum(atm_j_full_eval(:full_eval_cnt))
!    write(DBG_ARRAYS+mytaskid,'(I0.2," ",I0.10,":",I)') 3, proc_atm_to_full_list(atm_i), sum(atm_j_ee_eval(:ee_eval_cnt))
!  endif
!#endif

      ! pack the neighbors in the master array
      start_index(atm_i) = num_packed+1
      total_pairs = total_pairs + full_eval_cnt + ee_eval_cnt       ! DBG
      call pack_nb_list_midpoint(ee_eval_cnt, atm_j_ee_eval, &
                    full_eval_cnt, atm_j_full_eval, &
                    proc_ipairs, num_packed)
      if (ifail .eq. 1) then
        return !go back and increase the proc_ipaires size
      end if

      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(proc_atm_to_full_list(atm_i))%cnt
        if(excl_img_flags2(atm_mask(i)) .ne. 0) then
          atm_j_excl_list(excl_img_flags2(atm_mask(i))) = 0
        end if
      end do

    enddo ! end atm_i loop

    ! reset the auxulary bond list indirect indices
    do neighbor_bkt = 1,  proc_corner_neighborbuckets_cnt(bkt)
      cur_bkt = proc_neighborbuckets(neighbor_bkt,bkt)
      do atm_j_idx = 1, bkt_atm_cnt(cur_bkt)
        atm_j = bkt_atm_lst(atm_j_idx, cur_bkt)
        excl_img_flags2(proc_atm_to_full_list(atm_j)) = 0
      end do
    end do

!  enddo ! loop over all bucket
!  enddo ! loop over all bucket
  enddo ! loop over all bucket

!write(0,"('0-',I0.3,' BKT atoms min: ',I4,' max: ',I4, ' mean: ', I4, ' sum: ',I6)") mytaskid,minval(bkt_atm_cnt), maxval(bkt_atm_cnt), sum(bkt_atm_cnt)/size(bkt_atm_cnt), sum(bkt_atm_cnt)

#if 0
call MPI_Reduce(atm_cnt, max_atms, 1, MPI_INTEGER, MPI_MAX, 0, pmemd_comm, err_code_mpi)
call MPI_Reduce(atm_cnt, min_atms, 1, MPI_INTEGER, MPI_MIN, 0, pmemd_comm, err_code_mpi)
call MPI_Reduce(atm_cnt, mean_atms, 1, MPI_INTEGER, MPI_SUM, 0, pmemd_comm, err_code_mpi)
mean_atms = mean_atms/numtasks
if(mytaskid == 0) write(0,"('8- MPI atoms count max: ',I,' min: ' I,' mean: ' I, ' max-mean range % :' F7.2)") &
                max_atms, min_atms , mean_atms, (max_atms-mean_atms)*100.d0/mean_atms
#endif


!print*, mytaskid, mp_cnt, no_mp_cnt, 1.0*no_mp_cnt/(no_mp_cnt+mp_cnt)
!  call MPI_Reduce(num_packed, tot_packed, 1, MPI_INTEGER, MPI_SUM, 0, pmemd_comm, err_code_mpi)
!  call MPI_Reduce(total_pairs, all_pairs, 1, MPI_INTEGER, MPI_SUM, 0, pmemd_comm, err_code_mpi)

!  if(mytaskid .eq. 0) write(0, *) "TIME STEP *******************"
!  call MPI_Reduce(num_packed, tot_packed, 1, MPI_INTEGER, MPI_SUM, 0, pmemd_comm, err_code_mpi)
!  call MPI_Reduce(total_pairs, all_pairs, 1, MPI_INTEGER, MPI_SUM, 0, pmemd_comm, err_code_mpi)
!  if(mytaskid .eq. 0) write(mdout, *) "DBG", all_pairs, tot_packed
  return

contains
#include "pack_nb_list_dflt_midpoint.i"

end subroutine get_nb_list_midpoint

!*******************************************************************************
! Subroutine:  proc_check_my_atom_movement
!
! Description:  This routine checks if any atom on this task's atom list has
!               moved more than half the skin distance.  This is intended for
!               use only with parallel MD.
!
!*******************************************************************************

subroutine proc_check_my_atom_movement(crd, saved_crd, skinnb, ntp, &
                                  new_list)

  use parallel_dat_mod
  use pbc_mod
  use processor_mod, only: proc_num_atms, proc_atm_to_full_list

  implicit none

! Formal arguments:

  double precision      :: crd(3, *)
  double precision      :: saved_crd(3, *)
  double precision      :: skinnb
  integer               :: ntp
  logical               :: new_list

! Local variables:

  integer               :: idx
  double precision      :: dx, dy, dz
  double precision      :: box_dx, box_dy, box_dz
  double precision      :: distance_limit
  integer               :: n

! If anybody moves more than half the nbskin added to cutoff in generating the
! pairlist, update the pairlist.

  distance_limit = 0.25d0 * skinnb * skinnb

  if (ntp .eq. 0) then  ! Constant volume

    do idx = 1, proc_num_atms
      dx = crd(1, idx) - saved_crd(1, idx)
      dy = crd(2, idx) - saved_crd(2, idx)
      dz = crd(3, idx) - saved_crd(3, idx)

      if (dx * dx + dy * dy + dz * dz .ge. distance_limit) then
        !print *, crd(1, idx), saved_crd(1, idx), idx, mytaskid
        !print *, crd(2, idx), saved_crd(2, idx), idx, mytaskid
        !print *, crd(3, idx), saved_crd(3, idx), idx, mytaskid
        new_list = .true.
        return
      end if

    end do

  else  ! Constant pressure scaling.

    box_dx = pbc_box(1) / gbl_saved_box(1)
    box_dy = pbc_box(2) / gbl_saved_box(2)
    box_dz = pbc_box(3) / gbl_saved_box(3)

    do idx = 1, proc_num_atms

      dx = crd(1, idx) - saved_crd(1, idx) * box_dx
      dy = crd(2, idx) - saved_crd(2, idx) * box_dy
      dz = crd(3, idx) - saved_crd(3, idx) * box_dz

      if (dx * dx + dy * dy + dz * dz .ge. distance_limit) then
        new_list = .true.
        return
      end if

    end do

  end if

  new_list = .false.

  return

end subroutine proc_check_my_atom_movement
!*******************************************************************************
!
! Subroutine:  proc_save_all_atom_crds
!
! Description:  This is needed for the skin test for buffered pairlists.
!*******************************************************************************

subroutine proc_save_all_atom_crds(crd, saved_crd)

  use processor_mod, only: proc_num_atms, proc_atm_alloc_size,&
             proc_num_atms_min_bound,proc_ghost_num_atms
  implicit none

! Formal arguments:

  double precision      :: crd(3, proc_atm_alloc_size)
  double precision      :: saved_crd(3, proc_atm_alloc_size)

! Local variables:

  integer               :: i

  do i = 1, proc_num_atms+proc_ghost_num_atms
    saved_crd(1, i) = crd(1, i)
    saved_crd(2, i) = crd(2, i)
    saved_crd(3, i) = crd(3, i)
  end do

  return

end subroutine proc_save_all_atom_crds
#endif /* MPI */

end module nb_pairlist_mod
