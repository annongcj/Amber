#include "copyright.i"
!precision setting for SPDP, DPDP, SPSP
#include "include_precision.i"

!*******************************************************************************
!
! Module:  cit_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module cit_mod

  implicit none

  type cit_tbl_rec
    integer             :: img_lo
    integer             :: img_hi
  end type cit_tbl_rec
  type proc_cit_tbl_rec
    integer             :: atm_lo
    integer             :: atm_hi
  end type proc_cit_tbl_rec
  Double Precision, save         :: bkt_size(3) 
  integer, save         :: cit_tbl_x_dim
  integer, save         :: cit_tbl_y_dim
  integer, save         :: cit_tbl_z_dim

  integer, parameter    :: cit_bkt_delta = 2 ! used in pairlist code.

contains

!*******************************************************************************
!
! Subroutine:  setup_crd_idx_tbl
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine setup_crd_idx_tbl(atm_cnt, fraction, crd_idx_lst_tbl, atm_lst)

  use gbl_datatypes_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: fraction(3, atm_cnt)
  integer               :: crd_idx_lst_tbl(0 : cit_tbl_x_dim - 1, &
                                           0 : cit_tbl_y_dim - 1, &
                                           0 : cit_tbl_z_dim - 1)

  type(atm_lst_rec)     :: atm_lst(atm_cnt)

! Local variables:

  integer               :: atm_id
  integer               :: nxt_atm_lst_idx
  integer               :: x_idx, y_idx, z_idx
  double precision      :: scale_fac_x, scale_fac_y, scale_fac_z

  integer               :: tail_idx_tbl(0 : cit_tbl_x_dim - 1, &
                                        0 : cit_tbl_y_dim - 1, &
                                        0 : cit_tbl_z_dim - 1)

! Pre-initialize as needed:

  nxt_atm_lst_idx = 1

  scale_fac_x = dble(cit_tbl_x_dim)
  scale_fac_y = dble(cit_tbl_y_dim)
  scale_fac_z = dble(cit_tbl_z_dim)

  crd_idx_lst_tbl(:,:,:) = 0                ! Marks empty entries 

! Load the atom id's:

  do atm_id = 1, atm_cnt

    x_idx = int(fraction(1, atm_id) * scale_fac_x)
    y_idx = int(fraction(2, atm_id) * scale_fac_y)
    z_idx = int(fraction(3, atm_id) * scale_fac_z)

    ! crd_idx_lst_tbl takes bucket x, y, z coordinates and gives you the first
    ! atom in the bucket.

    if (crd_idx_lst_tbl(x_idx, y_idx, z_idx) .eq. 0) then   ! New list.

      crd_idx_lst_tbl(x_idx, y_idx, z_idx) = nxt_atm_lst_idx

    else        ! List already started.  Follow the chain and add a node:

    ! atm_lst is a linked list that gives you the next atom in the linked list
    ! when you enter in the index of a certain atom.  Each linked list is sorted
    ! by bucket.  So say 1 -> 3 -> 8 -> 16.  All four of these atoms are in the
    ! same bucket.  2 -> 5->9->23.  All 4 of these are in the same bucket.  You
    ! use something like crd_idx_lst_tbl(1,3,5) to find the first atom in the bucket then
    ! iterate through the linked list to get the rest of the atoms. %nxt=0 means
    ! you have iterated through all atoms in the bucket.

      atm_lst(tail_idx_tbl(x_idx, y_idx, z_idx))%nxt = nxt_atm_lst_idx

    end if

    ! Add the node:

    atm_lst(nxt_atm_lst_idx)%idx = atm_id
    atm_lst(nxt_atm_lst_idx)%nxt = 0
    tail_idx_tbl(x_idx, y_idx, z_idx) = nxt_atm_lst_idx
    nxt_atm_lst_idx = nxt_atm_lst_idx + 1

  end do

  return

end subroutine setup_crd_idx_tbl

!*******************************************************************************
!
! Subroutine:  setup_cit_tbl_dims
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine set_cit_tbl_dims(pbc_box, list_cutoff, cut_factor)

  implicit none

! Formal arguments:

  double precision, intent(in)  :: pbc_box(3)
  double precision, intent(in)  :: list_cutoff
  double precision, intent(in)  :: cut_factor(3)

! Local variables:

  ! Used to avoid rounding problems in gridding:
  double precision, parameter   :: rnd_fudge = 0.00001d0

  bkt_size(1:3) = cut_factor(1:3) * list_cutoff * 0.5d0 + rnd_fudge

#ifdef larger_bucket /*have tried larger bucket but did not help*/
  cit_tbl_x_dim = floor(pbc_box(1) / &
                  (cut_factor(1) * list_cutoff * 1.d0 + rnd_fudge))
  cit_tbl_y_dim = floor(pbc_box(2) / &
                  (cut_factor(2) * list_cutoff * 1.d0 + rnd_fudge))
  cit_tbl_z_dim = floor(pbc_box(3) / &
                  (cut_factor(3) * list_cutoff * 1.d0 + rnd_fudge))
#else
  cit_tbl_x_dim = floor(pbc_box(1) / &
                  (cut_factor(1) * list_cutoff * 0.5d0 + rnd_fudge))
  cit_tbl_y_dim = floor(pbc_box(2) / &
                  (cut_factor(2) * list_cutoff * 0.5d0 + rnd_fudge))
  cit_tbl_z_dim = floor(pbc_box(3) / &
                  (cut_factor(3) * list_cutoff * 0.5d0 + rnd_fudge))
#endif
  return

end subroutine set_cit_tbl_dims

!*******************************************************************************
!
! Subroutine:  setup_cit
!
! Description: Create three arrays: crd_idx_tbl, atm_img_map, img_atm_map
!              
!*******************************************************************************

subroutine setup_cit(atm_cnt, fraction, crd_idx_tbl, &
                     atm_img_map, img_atm_map)

  use gbl_datatypes_mod
  use img_mod
  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: fraction(3, atm_cnt)
  type(cit_tbl_rec)     :: crd_idx_tbl(0 : cit_tbl_x_dim - 1, &
                                       0 : cit_tbl_y_dim - 1, &
                                       0 : cit_tbl_z_dim - 1)

  integer               :: atm_img_map(atm_cnt)
  integer               :: img_atm_map(atm_cnt)

! Local variables:

  integer               :: img_id
  integer               :: i, j, k
  integer               :: img_lo, img_hi
  integer               :: nxt_idx
  integer               :: crd_idx_lst_tbl(0 : cit_tbl_x_dim - 1, &
                                           0 : cit_tbl_y_dim - 1, &
                                           0 : cit_tbl_z_dim - 1)

  type(atm_lst_rec)     :: atm_lst(atm_cnt)

  call setup_crd_idx_tbl(atm_cnt, fraction, crd_idx_lst_tbl, atm_lst)

  img_hi = 0

  do k = 0, cit_tbl_z_dim - 1
    do j = 0, cit_tbl_y_dim - 1
      do i = 0, cit_tbl_x_dim - 1

        nxt_idx = crd_idx_lst_tbl(i, j, k)

        if (nxt_idx .ne. 0) then

          img_hi = img_hi + 1
          img_lo = img_hi

          do ! while do loop create the atm_img_map() array using atm_lst array
          ! BEGIN DBG
          ! if (atm_lst(nxt_idx)%idx .lt. 1 .or. &
          !     atm_lst(nxt_idx)%idx .gt. natom) then
          !   write(0,*)'DBG: Found bad value ', atm_lst(nxt_idx)%idx, &
          !             'in atm_lst!!!'
          ! end if
          ! END DBG
            atm_img_map(atm_lst(nxt_idx)%idx) = img_hi
            nxt_idx = atm_lst(nxt_idx)%nxt
            if (nxt_idx .ne. 0) then
              img_hi = img_hi + 1
            else
              exit
            end if
          end do
          !Here crd_idx_tbl()%lo and %hi images are stored for each bucket
          crd_idx_tbl(i, j, k)%img_lo = img_lo
          crd_idx_tbl(i, j, k)%img_hi = img_hi
          !since the orginal value of nxt_idx is lost, get the value one
          nxt_idx = crd_idx_lst_tbl(i, j, k)

          do img_id = img_lo, img_hi! this do loop create the img_atm_map() array using atm_lst array
            img_atm_map(img_id) = atm_lst(nxt_idx)%idx
            nxt_idx = atm_lst(nxt_idx)%nxt
          end do

        else ! no atom in this box 
          crd_idx_tbl(i, j, k)%img_lo = 0
          crd_idx_tbl(i, j, k)%img_hi = -1
        end if

      end do
    end do
  end do

  return

end subroutine setup_cit

!*******************************************************************************
!
! Subroutine:  setup_cit_tbl_bkts
!
! Description:  This subroutine sets up tables used to handle pbc wrapping
!               and translation vectors in working with a flat cit.  The
!               vector arguments are optional.
!              
!*******************************************************************************

subroutine setup_cit_tbl_bkts(x_bkts, y_bkts, z_bkts, x_trans, y_trans, z_trans)

  implicit none

! Formal arguments:

  integer, intent(out)  :: x_bkts(0 : cit_tbl_x_dim * 3 - 1)
  integer, intent(out)  :: y_bkts(0 : cit_tbl_y_dim * 3 - 1)
  integer, intent(out)  :: z_bkts(0 : cit_tbl_z_dim * 2 - 1)
  integer, optional, intent(out)  :: x_trans(0 : cit_tbl_x_dim * 3 - 1)
  integer, optional, intent(out)  :: y_trans(0 : cit_tbl_y_dim * 3 - 1)
  integer, optional, intent(out)  :: z_trans(0 : cit_tbl_z_dim * 2 - 1)

! Local variables:

  integer                       :: i

  if (present(x_trans)) then
    do i = 0, cit_tbl_x_dim - 1
      x_bkts(i) = i
      x_bkts(i + cit_tbl_x_dim) = i
      x_bkts(i + cit_tbl_x_dim + cit_tbl_x_dim) = i
      x_trans(i) = 0
      x_trans(i + cit_tbl_x_dim) = 1
      x_trans(i + cit_tbl_x_dim + cit_tbl_x_dim) = 2
    end do
  else
    do i = 0, cit_tbl_x_dim - 1
      x_bkts(i) = i
      x_bkts(i + cit_tbl_x_dim) = i
      x_bkts(i + cit_tbl_x_dim + cit_tbl_x_dim) = i
    end do
  end if

  if (present(y_trans)) then
    do i = 0, cit_tbl_y_dim - 1
      y_bkts(i) = i * cit_tbl_x_dim
      y_bkts(i + cit_tbl_y_dim) = i * cit_tbl_x_dim
      y_bkts(i + cit_tbl_y_dim + cit_tbl_y_dim) = i * cit_tbl_x_dim
      y_trans(i) = 0
      y_trans(i + cit_tbl_y_dim) = 3
      y_trans(i + cit_tbl_y_dim + cit_tbl_y_dim) = 6
    end do
  else
    do i = 0, cit_tbl_y_dim - 1
      y_bkts(i) = i * cit_tbl_x_dim
      y_bkts(i + cit_tbl_y_dim) = i * cit_tbl_x_dim
      y_bkts(i + cit_tbl_y_dim + cit_tbl_y_dim) = i * cit_tbl_x_dim
    end do
  end if

  if (present(z_trans)) then
    do i = 0, cit_tbl_z_dim - 1
      z_bkts(i) = i * cit_tbl_x_dim * cit_tbl_y_dim
      z_bkts(i + cit_tbl_z_dim) = i * cit_tbl_x_dim * cit_tbl_y_dim
      z_trans(i) = 0
      z_trans(i + cit_tbl_z_dim) = 9
    end do
  else
    do i = 0, cit_tbl_z_dim - 1
      z_bkts(i) = i * cit_tbl_x_dim * cit_tbl_y_dim
      z_bkts(i + cit_tbl_z_dim) = i * cit_tbl_x_dim * cit_tbl_y_dim
    end do
  end if

  return

end subroutine setup_cit_tbl_bkts

!*******************************************************************************
!
! Subroutine:  get_flat_cit_idx
!
! Description:  Find the flat cit "bucket" a given image is in.  Only valid
!               during a list build step.
!              
!*******************************************************************************

subroutine get_flat_cit_idx(img_id, img_atm_map, fraction, flat_cit_idx)

  implicit none

! Formal arguments:

  integer, intent(in)           :: img_id
  integer, intent(in)           :: img_atm_map(*)
  double precision, intent(in)  :: fraction(3, *)
  integer, intent(out)          :: flat_cit_idx

! Local variables:

  integer                       :: atm_id
  integer                       :: x_idx, y_idx, z_idx

  atm_id = img_atm_map(img_id)

  x_idx = int(fraction(1, atm_id) * dble(cit_tbl_x_dim))
  y_idx = int(fraction(2, atm_id) * dble(cit_tbl_y_dim))
  z_idx = int(fraction(3, atm_id) * dble(cit_tbl_z_dim))

  flat_cit_idx = x_idx + cit_tbl_x_dim * (y_idx + z_idx * cit_tbl_y_dim)

  return

end subroutine get_flat_cit_idx

#ifdef MPI

!*******************************************************************************
!
! Subroutine:  proc_get_flat_cit_idx
!
! Description:  Find the flat cit "bucket" a given image is in.  Only valid
!               during a list build step.
!              
!*******************************************************************************

subroutine proc_get_flat_cit_idx(atm_id, fraction, flat_cit_idx)

  implicit none

! Formal arguments:

  integer, intent(in)           :: atm_id
  double precision, intent(in)  :: fraction(3, *)
  integer, intent(out)          :: flat_cit_idx

! Local variables:

  integer                       :: x_idx, y_idx, z_idx

  x_idx = int(fraction(1, atm_id) * dble(cit_tbl_x_dim))
  y_idx = int(fraction(2, atm_id) * dble(cit_tbl_y_dim))
  z_idx = int(fraction(3, atm_id) * dble(cit_tbl_z_dim))

  flat_cit_idx = x_idx + cit_tbl_x_dim * (y_idx + z_idx * cit_tbl_y_dim)

  return

end subroutine proc_get_flat_cit_idx

!*******************************************************************************
!
! Subroutine:  proc_setup_crd_idx_tbl
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine proc_setup_crd_idx_tbl(atm_cnt, fraction, crd_idx_lst_tbl, tail_idx_tbl, atm_lst,int_x_dim, int_y_dim,int_z_dim)

  use gbl_datatypes_mod
  use parallel_dat_mod
  use processor_mod, only: proc_atm_alloc_size
  implicit none

! Formal arguments:

  integer               :: atm_cnt
  integer               :: int_x_dim
  integer               :: int_y_dim
  integer               :: int_z_dim
#if 0
!#ifdef pmemd_SPDP
  pme_float             :: fraction(:,:)
#else
  double precision      :: fraction(:,:)
#endif
  !originally it was -1, but now it is +1, since two buckets are for ghosts
  integer               :: crd_idx_lst_tbl(0 : int_x_dim + 1, &
                                           0 : int_y_dim + 1, &
                                           0 : int_z_dim + 1)
  !originally it was -1, but now it is +1, since two buckets are for ghosts
  integer               :: tail_idx_tbl(0 : int_x_dim + 1, &
                                        0 : int_y_dim + 1, &
                                        0 : int_z_dim + 1)


  type(atm_lst_rec)     :: atm_lst(proc_atm_alloc_size)

! Local variables:

#if 0
!#ifdef pmemd_SPDP
  pme_float             :: scale_fac_x, scale_fac_y, scale_fac_z
#else
  double precision      :: scale_fac_x, scale_fac_y, scale_fac_z
#endif
  integer               :: atm_id
  integer               :: nxt_atm_lst_idx
  integer               :: x_idx, y_idx, z_idx

! Pre-initialize as needed:

  nxt_atm_lst_idx = 1

  scale_fac_x = dble(int_x_dim+2)!2 for ghost buckets each side
  scale_fac_y = dble(int_y_dim+2)!2 for ghost buckets each side
  scale_fac_z = dble(int_z_dim+2)!2 for ghost buckets each side

  crd_idx_lst_tbl(:,:,:) = 0                ! Marks empty entries 
  atm_lst(:)%idx = 0 ! may be a bit expensive, we can reduce the size
  atm_lst(:)%nxt = 0 ! may be a bit expensive, we can reduce the size

! Load the atom id's:

  do atm_id = 1, atm_cnt

    x_idx = int(fraction(1, atm_id) * scale_fac_x)
    y_idx = int(fraction(2, atm_id) * scale_fac_y)
    z_idx = int(fraction(3, atm_id) * scale_fac_z)

    ! crd_idx_lst_tbl takes bucket x, y, z coordinates and gives you the first
    ! atom in the bucket.

    if (crd_idx_lst_tbl(x_idx, y_idx, z_idx) .eq. 0) then   ! New list.

      crd_idx_lst_tbl(x_idx, y_idx, z_idx) = nxt_atm_lst_idx

    else        ! List already started.  Follow the chain and add a node:

    ! atm_lst is a linked list that gives you the next atom in the linked list
    ! when you enter in the index of a certain atom.  Each linked list is sorted
    ! by bucket.  So say 1 -> 3 -> 8 -> 16.  All four of these atoms are in the
    ! same bucket.  2 -> 5->9->23.  All 4 of these are in the same bucket.  You
    ! use something like crd_idx_lst_tbl(1,3,5) to find the first atom in the bucket then
    ! iterate through the linked list to get the rest of the atoms. %nxt=0 means
    ! you have iterated through all atoms in the bucket.

      atm_lst(tail_idx_tbl(x_idx, y_idx, z_idx))%nxt = nxt_atm_lst_idx

    end if

    ! Add the node:

    atm_lst(nxt_atm_lst_idx)%idx = atm_id
    atm_lst(nxt_atm_lst_idx)%nxt = 0
    tail_idx_tbl(x_idx, y_idx, z_idx) = nxt_atm_lst_idx
    nxt_atm_lst_idx = nxt_atm_lst_idx + 1
  end do
  nxt_atm_lst_idx=0
!Debug
! do x_idx=0,int_x_dim+1
!   do y_idx=0, int_y_dim+1
!     do z_idx=0, int_z_dim+1
!       write(15+mytaskid,*)"CIT: ", x_idx, y_idx, z_idx, &
!                    crd_idx_lst_tbl(x_idx,y_idx,z_idx)
!       if(crd_idx_lst_tbl(x_idx,y_idx,z_idx) .gt. 0) then
!         nxt_atm_lst_idx=nxt_atm_lst_idx+1
!       end if
!     end do
!   end do
! end do
!
! write(15+mytaskid,*)"CIT Count: ",nxt_atm_lst_idx

  return

end subroutine proc_setup_crd_idx_tbl

!*******************************************************************************
!
! Subroutine:  proc_setup_crd_idx_tbl_ghost
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine proc_setup_crd_idx_tbl_ghost(atm_cnt, ghost_atm_cnt, fraction, &
           crd_idx_lst_tbl, tail_idx_tbl, atm_lst, int_x_dim, int_y_dim, int_z_dim)

  use gbl_datatypes_mod
  use processor_mod, only: int_bkt_minx, int_bkt_miny, int_bkt_minz, &
                           proc_bkt_minx, proc_bkt_miny, proc_bkt_minz, &
                           proc_old_local_id, proc_atm_crd, &
                           proc_atm_to_full_list
  use parallel_dat_mod
  implicit none

! Formal arguments:

  integer               :: atm_cnt
  integer               :: ghost_atm_cnt
  integer               :: int_x_dim
  integer               :: int_y_dim
  integer               :: int_z_dim
#if 0
!#ifdef pmemd_SPDP
  pme_float             :: fraction(:, :)
#else
  double precision      :: fraction(:, :)
#endif
  !originally it was -1, but now it is +1, since two buckets are for ghosts
  integer               :: crd_idx_lst_tbl(0 : int_x_dim + 1, &
                                           0 : int_y_dim + 1, &
                                           0 : int_z_dim + 1)
  integer              :: tail_idx_tbl(0 : int_x_dim + 1, &
                                        0 : int_y_dim + 1, &
                                        0 : int_z_dim + 1)
 
  type(atm_lst_rec)     :: atm_lst(:)

! Local variables:
#if 0
!#ifdef pmemd_SPDP
  pme_float             :: scale_fac_x, scale_fac_y, scale_fac_z
#else
  double precision      :: scale_fac_x, scale_fac_y, scale_fac_z
#endif
  integer               :: atm_id
  integer               :: nxt_atm_lst_idx
  integer               :: x_idx, y_idx, z_idx


! Pre-initialize as needed:

  nxt_atm_lst_idx = atm_cnt+1

! We work in external buckets first for ghosts

  scale_fac_x = dble(int_x_dim+2)!2 for ghost buckets each side
  scale_fac_y = dble(int_y_dim+2)!2 for ghost buckets each side
  scale_fac_z = dble(int_z_dim+2)!2 for ghost buckets each side

! Load the ghost atom id's:

  do atm_id = atm_cnt+1, atm_cnt + ghost_atm_cnt

    x_idx = int(fraction(1, atm_id) * scale_fac_x)
    y_idx = int(fraction(2, atm_id) * scale_fac_y)
    z_idx = int(fraction(3, atm_id) * scale_fac_z)

!if(x_idx .gt. int_x_dim+1) print*, "x g overflow", x_idx, int_x_dim+1, fraction(1, atm_id), scale_fac_x 
!if(y_idx .gt. int_y_dim+1) print*, "y g overflow", y_idx, int_y_dim+1, fraction(2, atm_id), scale_fac_y 
!if(z_idx .gt. int_z_dim+1) print*, "z g overflow", z_idx, int_z_dim+1, fraction(3, atm_id), scale_fac_z

!   if(proc_old_local_id(atm_id) .eq. 180) then
!     write(15+mytaskid,*)"atm 10",x_idx,y_idx,z_idx
!   end if
!   write(15+mytaskid,*)"cit: ",x_idx,y_idx,z_idx

    ! crd_idx_lst_tbl takes bucket x, y, z coordinates and gives you the first
    ! atom in the bucket.

    if (crd_idx_lst_tbl(x_idx, y_idx, z_idx) .eq. 0) then   ! New list.

      crd_idx_lst_tbl(x_idx, y_idx, z_idx) = nxt_atm_lst_idx
!write(0,"(I)") nxt_atm_lst_idx

    else        ! List already started.  Follow the chain and add a node:

    ! atm_lst is a linked list that gives you the next atom in the linked list
    ! when you enter in the index of a certain atom.  Each linked list is sorted
    ! by bucket.  So say 1 -> 3 -> 8 -> 16.  All four of these atoms are in the
    ! same bucket.  2 -> 5->9->23.  All 4 of these are in the same bucket.  You
    ! use something like crd_idx_lst_tbl(1,3,5) to find the first atom in the bucket then
    ! iterate through the linked list to get the rest of the atoms. %nxt=0 means
    ! you have iterated through all atoms in the bucket.
if (tail_idx_tbl(x_idx, y_idx, z_idx) == -1) print*, x_idx, y_idx, z_idx, atm_cnt+1, atm_id
      atm_lst(tail_idx_tbl(x_idx, y_idx, z_idx))%nxt = nxt_atm_lst_idx

    end if

    ! Add the node:

    atm_lst(nxt_atm_lst_idx)%idx = atm_id
    atm_lst(nxt_atm_lst_idx)%nxt = 0
    tail_idx_tbl(x_idx, y_idx, z_idx) = nxt_atm_lst_idx
    nxt_atm_lst_idx = nxt_atm_lst_idx + 1
  end do

!Debug
! do x_idx=0,int_x_dim+1
!   do y_idx=0, int_y_dim+1
!     do z_idx=0, int_z_dim+1
!       write(15+mytaskid,*)"CIT: ", x_idx, y_idx, z_idx, &
!                    crd_idx_lst_tbl(x_idx,y_idx,z_idx)
!       if(crd_idx_lst_tbl(x_idx,y_idx,z_idx) .gt. 0) then
!         nxt_atm_lst_idx=nxt_atm_lst_idx+1
!       end if
!     end do
!   end do
! end do

  return

end subroutine proc_setup_crd_idx_tbl_ghost
#endif /*MPI*/

end module cit_mod
