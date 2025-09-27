#include "copyright.i"

!*******************************************************************************
!
! Module: shake_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module shake_mod

  implicit none

  type shake_bond_rec
    integer             :: atm_i
    integer             :: atm_j
    double precision    :: parm
    double precision    :: mass_i
  end type shake_bond_rec

  ! Note: fastwat_res_lst and my_fastwat_res_lst are actually lists of the
  ! first atm_id in each water residue and each water residue processed by
  ! this processor, respectively.

#if defined(MPI) && !defined(CUDA)
  integer, allocatable, save, private   :: fastwat_res_lst(:)
  integer, allocatable, save, private   :: nonfastwat_bond_lst(:)
#endif

  integer, save, private        :: fastwat_res_cnt = 0
  integer, save, private        :: nonfastwat_bond_cnt = 0

  integer, allocatable, save, private   :: my_fastwat_res_lst(:)
  integer, save, private        :: my_fastwat_res_cnt = 0
  integer, save, private        :: select_my_fastwat_res_cnt = 0

  type(shake_bond_rec), allocatable, save, private :: my_nonfastwat_bond_dat(:)
  integer, save, private                           :: my_nonfastwat_bond_cnt = 0

  integer, allocatable, private         :: my_fastwat_res_lst_midpoint(:)
  integer                               :: my_fastwat_res_cnt_midpoint 

  type(shake_bond_rec), allocatable, private       :: my_nonfastwat_bond_dat_midpoint(:)
  integer                                          :: my_nonfastwat_bond_cnt_midpoint
  integer                                          :: select_bond_cnt

  integer, save, private                           :: atm_cnt, res_cnt

  integer, save, private                           :: iorwat

  double precision, save, private                  :: rbtarg(8)
  double precision, save, private                  :: box_half(3)

  integer, save                                    :: noshakemask_cnt
  integer, save                                    :: num_noshake
  integer, allocatable, save, private              :: noshakemask_lst(:)
  integer, save                                    :: ibgwat

  integer, allocatable, save            :: mult_vec_wat(:,:)
  integer, allocatable, save            :: mult_vec_shake_bonds(:,:)
  private       get_water_distances, setlep_setup, wrap_crds

contains

#ifdef MPI

!*******************************************************************************
!
! Subroutine:   shake_claim_midpoint
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine shake_claim_midpoint()

  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
  use processor_mod
  use pbc_mod
  use parallel_dat_mod
  use ensure_alloc_mod

  implicit none

! Formal arguments:
  integer i, atm_i, atm_j, atm_k, lead_atom
  integer ind1, ind2, ind3, new_size
  double precision x1i, x2i, x3i, x1j, x2j, x3j, x1k, x2k, x3k
  double precision half_box_x, half_box_y, half_box_z
  logical atm1bool, atm2bool

! Allocate by 2/numtasks to ensure enough space
!  if(.not. allocated(my_nonfastwat_bond_dat_midpoint)) then
!    allocate(my_nonfastwat_bond_dat_midpoint(my_nonfastwat_bond_cnt*2*ghost__mult/numtasks))
!  end if
!  new_size = my_nonfastwat_bond_cnt*2*ghost__mult/numtasks
  new_size = proc_num_atms_min_bound+proc_ghost_num_atms
  if(new_size .gt. size(my_nonfastwat_bond_dat_midpoint) .and. allocated(my_nonfastwat_bond_dat_midpoint)) &
                                 deallocate(my_nonfastwat_bond_dat_midpoint)
  if(.not. allocated(my_nonfastwat_bond_dat_midpoint)) allocate(my_nonfastwat_bond_dat_midpoint(new_size))
  
!  if(.not. allocated(my_fastwat_res_lst_midpoint)) then
!    allocate(my_fastwat_res_lst_midpoint(my_fastwat_res_cnt*2*ghost__mult/numtasks))
!  end if
!  if(.not. allocated(mult_vec_wat)) then
!    allocate(mult_vec_wat(3,my_fastwat_res_cnt*2*3*ghost__mult/numtasks))
!  end if
!  if(.not. allocated(mult_vec_shake_bonds)) then
!    allocate(mult_vec_shake_bonds(3,my_nonfastwat_bond_cnt*2*ghost__mult/numtasks))
!  end if


!  call ensure_alloc_int(my_fastwat_res_lst_midpoint, my_fastwat_res_cnt*2*ghost__mult/numtasks)
!  call ensure_alloc_int3(mult_vec_wat, my_fastwat_res_cnt*2*3*ghost__mult/numtasks)
!  call ensure_alloc_int3(mult_vec_shake_bonds, my_nonfastwat_bond_cnt*2*ghost__mult/numtasks)
  call ensure_alloc_int(my_fastwat_res_lst_midpoint, new_size)
  call ensure_alloc_int3(mult_vec_wat, 2*3*new_size)
  call ensure_alloc_int3(mult_vec_shake_bonds, 2*3*new_size)

  if (iorwat .eq. 1) then
    ind1 = 0
    ind2 = 1
    ind3 = 2
  else if (iorwat .eq. 2) then
    ind1 = 1
    ind2 = 2
    ind3 = 0
  else
    ind1 = 2
    ind2 = 0
    ind3 = 1
  end if

  half_box_x = pbc_box(1)/2.0
  half_box_y = pbc_box(2)/2.0
  half_box_z = pbc_box(3)/2.0

  my_fastwat_res_cnt_midpoint = 0
  my_nonfastwat_bond_cnt_midpoint = 0

! Claim fast. Just checking bounds using atm_i should be good because
! theoretically they're all waters.
  do i = 1, my_fastwat_res_cnt
    lead_atom = my_fastwat_res_lst(i)
    atm_i = proc_shake_space(lead_atom + ind1)
    atm_j = proc_shake_space(lead_atom + ind2)
    atm_k = proc_shake_space(lead_atom + ind3)
    if(atm_i .eq. 0 .or. atm_j .eq. 0 .or. atm_k .eq. 0) cycle
! Assumption is if one atom is in our owned space they probably all are.
    if(atm_i .le. proc_num_atms .or. atm_j .le. proc_num_atms &
       .or. atm_k .le. proc_num_atms) then
      my_fastwat_res_cnt_midpoint = my_fastwat_res_cnt_midpoint + 1
      my_fastwat_res_lst_midpoint(my_fastwat_res_cnt_midpoint) = lead_atom
      ! Calculation of the bond vector:
      x1i = proc_atm_crd(1,atm_i)
      x2i = proc_atm_crd(2,atm_i)
      x3i = proc_atm_crd(3,atm_i)
      x1j = proc_atm_crd(1,atm_j)
      x2j = proc_atm_crd(2,atm_j)
      x3j = proc_atm_crd(3,atm_j)
      x1k = proc_atm_crd(1,atm_k)
      x2k = proc_atm_crd(2,atm_k)
      x3k = proc_atm_crd(3,atm_k)

      mult_vec_wat(1,3*my_fastwat_res_cnt_midpoint+ind1)=0
      mult_vec_wat(2,3*my_fastwat_res_cnt_midpoint+ind1)=0
      mult_vec_wat(3,3*my_fastwat_res_cnt_midpoint+ind1)=0

      !---xk-xj setup, used for rattle added by zhf
      if(abs(x1k - x1j) .lt. half_box_x) then
        mult_vec_wat(1,3*my_fastwat_res_cnt_midpoint+ind1) = 0
      else if(abs(x1k - (x1j + pbc_box(1)) ) .lt. half_box_x) then
        mult_vec_wat(1,3*my_fastwat_res_cnt_midpoint+ind1) = 1
      else if (abs(x1k - (x1j - pbc_box(1)) ) .lt. half_box_x) then
        mult_vec_wat(1,3*my_fastwat_res_cnt_midpoint+ind1) = -1
      end if
      if(abs(x2k - x2j) .lt. half_box_y) then
        mult_vec_wat(2,3*my_fastwat_res_cnt_midpoint+ind1) = 0
      else if(abs(x2k - (x2j + pbc_box(2)) ) .le. half_box_y) then
        mult_vec_wat(2,3*my_fastwat_res_cnt_midpoint+ind1) = 1
      else if (abs(x2k - (x2j - pbc_box(2)) ) .le. half_box_y) then
        mult_vec_wat(2,3*my_fastwat_res_cnt_midpoint+ind1) = -1
      end if
      if(abs(x3k - x3j) .lt. half_box_z) then
        mult_vec_wat(3,3*my_fastwat_res_cnt_midpoint+ind1) = 0
      else if(abs(x3k - (x3j + pbc_box(3)) ) .le. half_box_z) then
        mult_vec_wat(3,3*my_fastwat_res_cnt_midpoint+ind1) = 1
      else if (abs(x3k - (x3j - pbc_box(3)) ) .le. half_box_z) then
        mult_vec_wat(3,3*my_fastwat_res_cnt_midpoint+ind1) = -1
      end if


      if(abs(x1i - x1j) .lt. half_box_x) then
        mult_vec_wat(1,3*my_fastwat_res_cnt_midpoint+ind2) = 0
      else if(abs(x1i - (x1j + pbc_box(1)) ) .lt. half_box_x) then
        mult_vec_wat(1,3*my_fastwat_res_cnt_midpoint+ind2) = 1
      else if (abs(x1i - (x1j - pbc_box(1)) ) .lt. half_box_x) then
        mult_vec_wat(1,3*my_fastwat_res_cnt_midpoint+ind2) = -1
      end if
      if(abs(x2i - x2j) .lt. half_box_y) then
        mult_vec_wat(2,3*my_fastwat_res_cnt_midpoint+ind2) = 0
      else if(abs(x2i - (x2j + pbc_box(2)) ) .le. half_box_y) then
        mult_vec_wat(2,3*my_fastwat_res_cnt_midpoint+ind2) = 1
      else if (abs(x2i - (x2j - pbc_box(2)) ) .le. half_box_y) then
        mult_vec_wat(2,3*my_fastwat_res_cnt_midpoint+ind2) = -1
      end if
      if(abs(x3i - x3j) .lt. half_box_z) then
        mult_vec_wat(3,3*my_fastwat_res_cnt_midpoint+ind2) = 0
      else if(abs(x3i - (x3j + pbc_box(3)) ) .le. half_box_z) then
        mult_vec_wat(3,3*my_fastwat_res_cnt_midpoint+ind2) = 1
      else if (abs(x3i - (x3j - pbc_box(3)) ) .le. half_box_z) then
        mult_vec_wat(3,3*my_fastwat_res_cnt_midpoint+ind2) = -1
      end if

      if(abs(x1i - x1k) .lt. half_box_x) then
        mult_vec_wat(1,3*my_fastwat_res_cnt_midpoint+ind3) = 0
      else if(abs(x1i - (x1k + pbc_box(1)) ) .lt. half_box_x) then
        mult_vec_wat(1,3*my_fastwat_res_cnt_midpoint+ind3) = 1
      else if (abs(x1i - (x1k - pbc_box(1)) ) .lt. half_box_x) then
        mult_vec_wat(1,3*my_fastwat_res_cnt_midpoint+ind3) = -1
      end if
      if(abs(x2i - x2k) .lt. half_box_y) then
        mult_vec_wat(2,3*my_fastwat_res_cnt_midpoint+ind3) = 0
      else if(abs(x2i - (x2k + pbc_box(2)) ) .le. half_box_y) then
        mult_vec_wat(2,3*my_fastwat_res_cnt_midpoint+ind3) = 1
      else if (abs(x2i - (x2k - pbc_box(2)) ) .le. half_box_y) then
        mult_vec_wat(2,3*my_fastwat_res_cnt_midpoint+ind3) = -1
      end if
      if(abs(x3i - x3k) .lt. half_box_z) then
        mult_vec_wat(3,3*my_fastwat_res_cnt_midpoint+ind3) = 0
      else if(abs(x3i - (x3k + pbc_box(3)) ) .le. half_box_z) then
        mult_vec_wat(3,3*my_fastwat_res_cnt_midpoint+ind3) = 1
      else if (abs(x3i - (x3k - pbc_box(3)) ) .le. half_box_z) then
        mult_vec_wat(3,3*my_fastwat_res_cnt_midpoint+ind3) = -1
      end if
    end if
  end do

! Claim nonfast checking just i should be good because they're all bonds.
  do i = 1, my_nonfastwat_bond_cnt
    atm_i = my_nonfastwat_bond_dat(i)%atm_i
    atm_j = my_nonfastwat_bond_dat(i)%atm_j
!   if (noshakemask_cnt .gt. 0) then
!     if (noshakemask_lst(atm_i)+noshakemask_lst(atm_j) .ge. 1) then
!        cycle
!     end if
!   end if
    atm_i = proc_shake_space(atm_i)
    atm_j = proc_shake_space(atm_j)

    if(atm_i .eq. 0 .or. atm_j .eq. 0) cycle
!   if(atm_i .le. proc_num_atms .or. atm_j .le. proc_num_atms) then
      my_nonfastwat_bond_cnt_midpoint = my_nonfastwat_bond_cnt_midpoint+1
      my_nonfastwat_bond_dat_midpoint(my_nonfastwat_bond_cnt_midpoint)%atm_i= &
                                                 atm_i
      my_nonfastwat_bond_dat_midpoint(my_nonfastwat_bond_cnt_midpoint)%atm_j= &
                                                 atm_j
      my_nonfastwat_bond_dat_midpoint(my_nonfastwat_bond_cnt_midpoint)%parm= &
                                                 my_nonfastwat_bond_dat(i)%parm
      x1i = proc_atm_crd(1,atm_i)
      x2i = proc_atm_crd(2,atm_i)
      x3i = proc_atm_crd(3,atm_i)
      x1j = proc_atm_crd(1,atm_j)
      x2j = proc_atm_crd(2,atm_j)
      x3j = proc_atm_crd(3,atm_j)
      !unwrapped
      !atm_i will be unchanged
      !atm_j  have 3 variations
      if(abs(x1i - x1j) .le. half_box_x) then
        mult_vec_shake_bonds(1,my_nonfastwat_bond_cnt_midpoint) = 0
      else if(abs(x1i - (x1j + pbc_box(1)) ) .le. half_box_x) then
        mult_vec_shake_bonds(1,my_nonfastwat_bond_cnt_midpoint) = 1
      else if (abs(x1i - (x1j - pbc_box(1)) ) .le. half_box_x) then
        mult_vec_shake_bonds(1,my_nonfastwat_bond_cnt_midpoint) = -1
      end if
       if(abs(x2i - x2j) .le. half_box_y) then
         mult_vec_shake_bonds(2,my_nonfastwat_bond_cnt_midpoint) = 0
       else if(abs(x2i - (x2j + pbc_box(2)) ) .le. half_box_y) then
         mult_vec_shake_bonds(2,my_nonfastwat_bond_cnt_midpoint) = 1
       else if (abs(x2i - (x2j - pbc_box(2)) ) .le. half_box_y) then
         mult_vec_shake_bonds(2,my_nonfastwat_bond_cnt_midpoint) = -1
       end if
       if(abs(x3i - x3j) .le.half_box_z) then
         mult_vec_shake_bonds(3,my_nonfastwat_bond_cnt_midpoint) = 0
       else if(abs(x3i - (x3j + pbc_box(3)) ) .le. half_box_z) then
         mult_vec_shake_bonds(3,my_nonfastwat_bond_cnt_midpoint) = 1
       else if (abs(x3i - (x3j - pbc_box(3)) ) .le. half_box_z) then
         mult_vec_shake_bonds(3,my_nonfastwat_bond_cnt_midpoint) = -1
       end if
!   end if
  end do

end subroutine shake_claim_midpoint

#endif /*MPI*/

#if defined(MPI) && !defined(CUDA)
!*******************************************************************************
!
! Subroutine:   claim_my_fastwat_residues
!
! Description:  <TBS>
!              
!*******************************************************************************

subroutine claim_my_fastwat_residues(num_ints, num_reals)

  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed
  integer               :: atm_id
  integer               :: i
  integer               :: my_fastwat_res(nres)

  if (imin .ne. 0) return       ! There is no reallocation for minimization.

  if (ntc .eq. 1) return        ! No shake...

  ! Find residues that you own.  Since we always divide on residue
  ! boundaries, we only need check the first residue.

  my_fastwat_res_cnt = 0

  do i = 1, fastwat_res_cnt
    atm_id = fastwat_res_lst(i)
    if (gbl_atm_owner_map(atm_id) .eq. mytaskid) then
      my_fastwat_res_cnt = my_fastwat_res_cnt + 1
      my_fastwat_res(my_fastwat_res_cnt) = atm_id
    end if
  end do

  ! We first constructed the fastwater residue list for this process on
  ! the stack because we did not know how big it would be.  Now we can
  ! allocate space and make a copy to static storage.  Note that this
  ! code allows for reallocation (ie., multiple calls).

  if (allocated(my_fastwat_res_lst)) then

    if (size(my_fastwat_res_lst) .lt. my_fastwat_res_cnt) then

      num_ints = num_ints - size(my_fastwat_res_lst)
      deallocate(my_fastwat_res_lst)

      allocate(my_fastwat_res_lst(my_fastwat_res_cnt), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(my_fastwat_res_lst)

    end if

  else

    allocate(my_fastwat_res_lst(my_fastwat_res_cnt), stat = alloc_failed)
    if (alloc_failed .ne. 0) call setup_alloc_error
    num_ints = num_ints + size(my_fastwat_res_lst)

  end if

  do i = 1, my_fastwat_res_cnt
    my_fastwat_res_lst(i) = my_fastwat_res(i)
  end do

  return

end subroutine claim_my_fastwat_residues

#ifdef _OPENMP_
subroutine claim_my_fastwat_residues_gb(num_ints, num_reals)

  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed
  integer               :: atm_id
  integer               :: i
  integer               :: my_fastwat_res(nres)

  if (imin .ne. 0) return       ! There is no reallocation for minimization.

  if (ntc .eq. 1) return        ! No shake...

  ! Find residues that you own.  Since we always divide on residue
  ! boundaries, we only need check the first residue.

  my_fastwat_res_cnt = 0

  if(master) then
      my_fastwat_res(1:fastwat_res_cnt) = &
         fastwat_res_lst(1: fastwat_res_cnt) 
      my_fastwat_res_cnt = fastwat_res_cnt
  end if

  ! We first constructed the fastwater residue list for this process on
  ! the stack because we did not know how big it would be.  Now we can
  ! allocate space and make a copy to static storage.  Note that this
  ! code allows for reallocation (ie., multiple calls).

  if (allocated(my_fastwat_res_lst)) then
    if (size(my_fastwat_res_lst) .lt. my_fastwat_res_cnt) then
      num_ints = num_ints - size(my_fastwat_res_lst)
      deallocate(my_fastwat_res_lst)
      allocate(my_fastwat_res_lst(my_fastwat_res_cnt), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(my_fastwat_res_lst)
    end if
  else

    allocate(my_fastwat_res_lst(my_fastwat_res_cnt), stat = alloc_failed)
    if (alloc_failed .ne. 0) call setup_alloc_error
    num_ints = num_ints + size(my_fastwat_res_lst)

  end if

  do i = 1, my_fastwat_res_cnt
    my_fastwat_res_lst(i) = my_fastwat_res(i)
  end do

  return

end subroutine claim_my_fastwat_residues_gb
#endif /*_OPENMP_*/

!*******************************************************************************
!
!*******************************************************************************
!
! Subroutine:   claim_my_nonfastwat_bonds
!
! Description:  <TBS>
!              
!*******************************************************************************

subroutine claim_my_nonfastwat_bonds(num_ints, num_reals)

  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed
  integer               :: atm_i, atm_j
  integer               :: bond_idx        
  integer               :: bonddat_idx        
  integer               :: bondlst_idx

  if (imin .ne. 0) return       ! There is no reallocation for minimization.

  if (ntc .eq. 1) return        ! No shake...

  ! Determine how many nonfastwater bonds you actually own.

  my_nonfastwat_bond_cnt = 0

  do bondlst_idx = 1, nonfastwat_bond_cnt

    bond_idx = nonfastwat_bond_lst(bondlst_idx)

    atm_i = gbl_bond(bond_idx)%atm_i
    atm_j = gbl_bond(bond_idx)%atm_j

    if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then

      my_nonfastwat_bond_cnt = my_nonfastwat_bond_cnt + 1
          
      if (gbl_atm_owner_map(atm_j) .ne. mytaskid) then
        write(mdout, *) 'Partition error in shake, task ', mytaskid
        call mexit(6, 1)
      end if

    else

      if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
        write(mdout, *) 'Partition error in shake, task ', mytaskid
        call mexit(6, 1)
      end if

    end if

  end do

  ! Allocate space. Note it is reallocatable.

  if (allocated(my_nonfastwat_bond_dat)) then

    if (size(my_nonfastwat_bond_dat) .lt. my_nonfastwat_bond_cnt) then

      num_ints = num_ints - 2 * size(my_nonfastwat_bond_dat)
      num_reals = num_reals - size(my_nonfastwat_bond_dat)

      deallocate(my_nonfastwat_bond_dat)

      allocate(my_nonfastwat_bond_dat(my_nonfastwat_bond_cnt), &
               stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error

      num_ints = num_ints + 2 * size(my_nonfastwat_bond_dat)
      num_reals = num_reals + size(my_nonfastwat_bond_dat)

    end if

  else

    allocate(my_nonfastwat_bond_dat(my_nonfastwat_bond_cnt), &
             stat = alloc_failed)
    if (alloc_failed .ne. 0) call setup_alloc_error

    num_ints = num_ints + 2 * size(my_nonfastwat_bond_dat)
    num_reals = num_reals + size(my_nonfastwat_bond_dat)

  end if

  bonddat_idx = 0

  ! Set up your bond data.

  do bondlst_idx = 1, nonfastwat_bond_cnt

    bond_idx = nonfastwat_bond_lst(bondlst_idx)

    atm_i = gbl_bond(bond_idx)%atm_i

    if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then

      bonddat_idx = bonddat_idx + 1

      my_nonfastwat_bond_dat(bonddat_idx)%atm_i = gbl_bond(bond_idx)%atm_i
      my_nonfastwat_bond_dat(bonddat_idx)%atm_j = gbl_bond(bond_idx)%atm_j
      my_nonfastwat_bond_dat(bonddat_idx)%parm = &
        gbl_req(gbl_bond(bond_idx)%parm_idx)**2 
      my_nonfastwat_bond_dat(bonddat_idx)%mass_i = &
        atm_mass(gbl_bond(bond_idx)%atm_i) 

    end if

  end do

  return

end subroutine claim_my_nonfastwat_bonds

#ifdef _OPENMP_
subroutine claim_my_nonfastwat_bonds_gb(num_ints, num_reals)

  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed
  integer               :: atm_i, atm_j
  integer               :: bond_idx        
  integer               :: bonddat_idx        
  integer               :: bondlst_idx

  if (imin .ne. 0) return       ! There is no reallocation for minimization.

  if (ntc .eq. 1) return        ! No shake...

  ! Determine how many nonfastwater bonds you actually own.

  my_nonfastwat_bond_cnt = 0

  if(master) then
    my_nonfastwat_bond_cnt = nonfastwat_bond_cnt
  end if


  bonddat_idx = 0

  if(master) then
      bonddat_idx = nonfastwat_bond_cnt
     do bondlst_idx = 1, nonfastwat_bond_cnt
      bond_idx = nonfastwat_bond_lst(bondlst_idx)
      my_nonfastwat_bond_dat(bondlst_idx)%atm_i = & 
            gbl_bond(bond_idx)%atm_i
      my_nonfastwat_bond_dat(bondlst_idx)%atm_j = & 
            gbl_bond(bond_idx)%atm_j
      my_nonfastwat_bond_dat(bondlst_idx)%parm = &
           gbl_req(gbl_bond(bond_idx)%parm_idx)**2 
     end do
  end if 

  return
end subroutine claim_my_nonfastwat_bonds_gb
#endif /*_OPENMP_*/

#endif /* MPI */

!*******************************************************************************
!
! Subroutine:   shake_setup
!
! Description:
!
! This routine determines which shaken bonds correspond to 3-point solvent
! solvent molecules (typically TIP3P or SPC waters). For these bonds,
! fastwat_bonds(i) is set to 1. For all other bonds, fastwat_bonds(i) is set
! to 0.  The code now also supports TIP4P and TIP5P water, basically by not
! insisting on a residue with 3 atoms.  A more thorough check would insure that
! the extra "atoms" are actually extra points, but since sander does not do
! this, we really don't need to do it either...
!
! Waters are defined as residues that meet the following criteria:
!   1) residue name = watnam
!   2) atom names are owatnm, hwatnm(1) and hwatnm(2) 
!
! Author: David A. Pearlman
! Date: 12/92
!
! Input:
!
! natom: Number of atoms
! atm_igraph(i): The name of atom i (i4)
! nres: The number of residues
! gbl_res_atms(i): gbl_res_atms(i)->gbl_res_atms(i+1)-1 are atoms of residue i
! gbl_labres(i): The name of residue i
! nbonh: The number of bonds to hydrogens
! nbona: The number of bonds to non-hydrogens
! gbl_bond(i)%atm_i
! gbl_bond(i)%atm_j: The two atoms of bond i.
! ibelly: > 0 if belly run is being performed.
! atm_igroup(i): > 0 if atom i is part of moving belly (only if ibelly > 0).
! iwtnm: The name of residues to be considered water (i)
! iowtnm: The name of oxygen atom in the water residues (i)
! ihwtnm(2): The name of 2 hydrogen atoms in the water residues (i)
! jfastw: = 0: Use fast water routine (possibly with modified names)
!              Use partitioning for fast tip3p-tip3p interactions in
!              nonbon routine (sander only).
!           4: Do not use fast waters anywhere.
!
! Output:
!
! fastwat_bonds(i): = 0 if bond should be constrained by standard routine.
!                   = 1 if bond should be constrained by fast 3-point const.
!                      routine.
! fastwat_res(i): = 0 if bonds of residue are constrained by standard routine.
!                 = 1 if bonds of residue are constrained by fast 3-point
!                   const. routine.
!
! ibgwat: The first water residue found.
! iorwat: The position of the oxygen in each residue (1, 2, or 3). 
!              
!*******************************************************************************

subroutine shake_setup(num_ints, num_reals)

  use pmemd_lib_mod
  use constraints_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use ti_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Record types and parameters:

  type fw_atm_rec
    integer             :: res_id
    integer             :: atm_type     ! unknown, fw_oxy, fw_h1, or fw_h2
  end type fw_atm_rec

  integer, parameter    :: unknown = 0
  integer, parameter    :: fw_oxy  = 1
  integer, parameter    :: fw_h1   = 2
  integer, parameter    :: fw_h2   = 3

  integer, parameter    :: h1_h2_bond  = 1       ! TIP3
  integer, parameter    :: oxy_h1_bond = 2
  integer, parameter    :: oxy_h2_bond = 3

  type fw_res_rec
    logical             :: h1_h2_bond   ! .true. if h1_h2 bond found for res.
    logical             :: oxy_h1_bond  ! .true. if oxy_h1 bond found for res.
    logical             :: oxy_h2_bond  ! .true. if oxy_h2 bond found for res.
  end type fw_res_rec

! Local variables:

  integer               :: alloc_failed
  integer               :: atm_i, atm_j
  integer               :: atm_i_type, atm_j_type
  integer               :: atm_id
  integer               :: bond_cnt
  integer               :: bond_type
  integer               :: i, j
!  integer               :: ibgwat  ! modified by FENG PAN
  integer               :: h1_id, h2_id, oxy_id
  integer               :: res_id

  type(fw_atm_rec)      :: fw_atm(natom)
  type(fw_res_rec)      :: fw_res(nres)

  integer               :: fastwat_res(nres)
  integer               :: fastwat_bonds(nbonh + nbona)

! If this is a minimization run, only the master runs shake, so we can bag
! out of everything that is shake-related if this is not the master process.
#if defined(MPI) && !defined(CUDA)
if(.not. usemidpoint) then
  if (imin .ne. 0 .and. .not. master) return
endif
#endif

  atm_cnt = natom
  bond_cnt = nbonh + nbona
  res_cnt = nres

! Initialize the various data arrays and scalars as required:

  fastwat_bonds(:) = 0
  fastwat_res(:) = 0

  ibgwat = 0

! If jfastw = 4 user has requested we not use fast water routine in any case.

  if (jfastw .ne. 4) then

    fw_atm(:) = fw_atm_rec(0, unknown)
    fw_res(:) = fw_res_rec(.false., .false., .false.)

! Search the list of residues for those 1) with the appropriate residue name;
! 2) with the appropriate atom names:

! If using belly, skip any residue where all atoms of residue
! cannot move (where all atoms do not have atm_igroup(i) > 0).

    fastwat_res_cnt = 0

    outer: do i = 1, res_cnt

! Screen out residue if wrong resname:

      if (gbl_labres(i) .ne. iwtnm) cycle outer

! Screen out residue if an atom is fixed (belly) or if an atom name does not
! match or is duplicated.  Also note atom #s.

      oxy_id = 0
      h1_id = 0
      h2_id = 0
      do j = 1, 3
        atm_id = gbl_res_atms(i) + j - 1
        if (ibelly .gt. 0) then
          if (atm_igroup(atm_id) .le. 0) cycle outer
        end if
        if (atm_igraph(atm_id) .eq. iowtnm) then
          if (oxy_id .ne. 0) cycle outer
          oxy_id = atm_id
        else if (atm_igraph(atm_id) .eq. ihwtnm(1)) then
          if (h1_id .ne. 0) cycle outer
          h1_id = atm_id
        else if (atm_igraph(atm_id) .eq. ihwtnm(2)) then
          if (h2_id .ne. 0) cycle outer
          h2_id = atm_id
        else
          cycle outer
        end if
      end do

! If we get here, this is a water molecule. Store atom info.

      fw_atm(oxy_id)%res_id = i
      fw_atm(oxy_id)%atm_type = fw_oxy

      fw_atm(h1_id)%res_id = i
      fw_atm(h1_id)%atm_type = fw_h1

      fw_atm(h2_id)%res_id = i
      fw_atm(h2_id)%atm_type = fw_h2

! ibgwat is now just used as a flag to know when to set iorwat.

      if (ibgwat .eq. 0) then
        ibgwat = i
        if (oxy_id .eq. gbl_res_atms(i)) then
          iorwat = 1
        else if (oxy_id .eq. gbl_res_atms(i) + 1) then
          iorwat = 2
        else
          iorwat = 3
        end if
      end if

      fastwat_res(i) = 1

      fastwat_res_cnt = fastwat_res_cnt + 1

    end do outer

  ! Now construct the list of bonds in residues using gbl_bond and the list
  ! just constructed.  This does a number of consistency checks.

    do i = 1, bond_cnt

        atm_i = gbl_bond(i)%atm_i
        atm_j = gbl_bond(i)%atm_j

        res_id = fw_atm(atm_i)%res_id

        if (res_id .ne. 0 .and. res_id .eq. fw_atm(atm_j)%res_id) then
  
          atm_i_type = fw_atm(atm_i)%atm_type
          atm_j_type = fw_atm(atm_j)%atm_type

          if (atm_i_type .eq. fw_oxy .and. atm_j_type .eq. fw_h1) then
            bond_type = oxy_h1_bond
          else if (atm_i_type .eq. fw_h1 .and. atm_j_type .eq. fw_oxy) then
            bond_type = oxy_h1_bond
          else if (atm_i_type .eq. fw_oxy .and. atm_j_type .eq. fw_h2) then
            bond_type = oxy_h2_bond
          else if (atm_i_type .eq. fw_h2 .and. atm_j_type .eq. fw_oxy) then
            bond_type = oxy_h2_bond
          else if (atm_i_type .eq. fw_h1 .and. atm_j_type .eq. fw_h2) then
            bond_type = h1_h2_bond
          else if (atm_i_type .eq. fw_h2 .and. atm_j_type .eq. fw_h1) then
            bond_type = h1_h2_bond
          else
            if (master) write(mdout, 9000)
            call mexit(6, 1)
          end if

          if (bond_type .eq. oxy_h1_bond) then

            if (.not. fw_res(res_id)%oxy_h1_bond) then
              fw_res(res_id)%oxy_h1_bond = .true.
            else
              if (master) write(mdout, 9000)
              call mexit(6, 1)
            end if

          else if (bond_type .eq. oxy_h2_bond) then

            if (.not. fw_res(res_id)%oxy_h2_bond) then
              fw_res(res_id)%oxy_h2_bond = .true.
            else
              if (master) write(mdout, 9000)
              call mexit(6, 1)
            end if

          else

            if (.not. fw_res(res_id)%h1_h2_bond) then
              fw_res(res_id)%h1_h2_bond = .true.
            else
              if (master) write(mdout, 9000)
              call mexit(6, 1)
            end if

          end if

          fastwat_bonds(i) = 1 ! Mark bond as fast water:

        end if

    end do

  ! Final consistency checks:

    do i = 1, res_cnt

      if (fastwat_res(i) .eq. 1 ) then

        if (.not. fw_res(i)%oxy_h1_bond .or. &
            .not. fw_res(i)%oxy_h2_bond .or. &
            .not. fw_res(i)%h1_h2_bond) then
          if (master) write(mdout, 9000)
          call mexit(6, 1)
        end if

      else

        if (fw_res(i)%oxy_h1_bond .or. &
            fw_res(i)%oxy_h2_bond .or. &
            fw_res(i)%h1_h2_bond) then
          if (master) write(mdout, 9000)
          call mexit(6, 1)
        end if

      end if

    end do

    call get_water_distances
  
  else

    fastwat_res_cnt = 0
    my_fastwat_res_cnt = 0

  end if

  if (master) write(mdout, 9001) fastwat_res_cnt

#if defined(MPI) && !defined(CUDA)
  if (fastwat_res_cnt .ne. 0) then

if(.not. usemidpoint) then
    if (imin .eq. 0) then
      ! Create the list that includes ALL fastwater residues.  We will use
      ! this to construct a list of fastwater residues we own at various
      ! points in time.

      allocate(fastwat_res_lst(fastwat_res_cnt), stat = alloc_failed)

      if (alloc_failed .ne. 0) call setup_alloc_error

      num_ints = num_ints + size(fastwat_res_lst)

      fastwat_res_cnt = 0

      do i = 1, res_cnt
        if (fastwat_res(i) .ne. 0) then
          atm_id = gbl_res_atms(i)
          fastwat_res_cnt = fastwat_res_cnt + 1
          fastwat_res_lst(fastwat_res_cnt) = atm_id
        end if
      end do

      ! Find residues that you own.  Since we always divide on residue
      ! boundaries, we only need check the first residue.

      call claim_my_fastwat_residues(num_ints, num_reals)

    else
      ! Only here if in master process, or in uniprocessor code:
      allocate(my_fastwat_res_lst(fastwat_res_cnt), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(my_fastwat_res_lst)
      my_fastwat_res_cnt = 0
      do i = 1, res_cnt
        if (fastwat_res(i) .ne. 0) then
          atm_id = gbl_res_atms(i)
          my_fastwat_res_cnt = my_fastwat_res_cnt + 1
          my_fastwat_res_lst(my_fastwat_res_cnt) = atm_id
        end if
      end do
    end if
else ! usemidpoint 
      ! Only here if in master process, or in uniprocessor code:
      allocate(my_fastwat_res_lst(fastwat_res_cnt), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(my_fastwat_res_lst)
      my_fastwat_res_cnt = 0
      do i = 1, res_cnt
        if (fastwat_res(i) .ne. 0) then
          atm_id = gbl_res_atms(i)
          my_fastwat_res_cnt = my_fastwat_res_cnt + 1
          my_fastwat_res_lst(my_fastwat_res_cnt) = atm_id
        end if
      end do
endif ! not usemidpoint
  end if
#else /*Serial code*/
  if (fastwat_res_cnt .ne. 0) then

      ! Only here if in master process, or in uniprocessor code:
      allocate(my_fastwat_res_lst(fastwat_res_cnt), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(my_fastwat_res_lst)
      my_fastwat_res_cnt = 0
      do i = 1, res_cnt
        if (fastwat_res(i) .ne. 0) then
          atm_id = gbl_res_atms(i)
          my_fastwat_res_cnt = my_fastwat_res_cnt + 1
          my_fastwat_res_lst(my_fastwat_res_cnt) = atm_id
        end if
      end do
  end if
#endif

  call nonfastwat_shake_setup
  
#ifdef CUDA
  call gpu_shake_setup(atm_mass, my_nonfastwat_bond_cnt, my_nonfastwat_bond_dat, my_fastwat_res_cnt, iorwat, &
                       my_fastwat_res_lst, tishake)
#ifdef GTI
  if (ti_mode .ne. 0 .and. tishake .le. 1) call gpu_shake_ti_setup(atm_mass, ti_lst, ti_sc_lst, ti_weights)
#endif
#endif  

  if (master) then
    if (ti_mode .ne. 0 .and. ti_mode .ne. 1) then
      write(mdout, 9003) 1, num_noshake + ti_num_noshake(1)
      write(mdout, 9003) 2, num_noshake + ti_num_noshake(2)
    else
      if (noshakemask_cnt .gt. 0) then
        write(mdout, 9002) num_noshake
      end if   
    end if
  end if

  return

! Format statements:

9000 format(' Error: Fast 3-point water residue, name and bond data incorrect!')
9001 format(' Number of triangulated 3-point waters found: ', i8)
9002 format(' Number of shake restraints removed: ', i8)
9003 format(' Number of shake restraints removed in TI region ',i2,' : ', i8)
contains

!*******************************************************************************
!
! Internal Subroutine:  nonfastwat_shake_setup
!
! Description: <TBS>
!
!*******************************************************************************

subroutine nonfastwat_shake_setup

  use pmemd_lib_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use ti_mod

  implicit none

! Local variables:

  integer               :: alloc_failed
  integer               :: bond_cnt
  integer               :: bond_idx        
  integer               :: bonddat_idx        
  integer               :: bondlst_idx
  integer               :: resi, resj
  integer               :: i
  integer               :: j
  integer               :: new_atm
  logical               :: skip
  ! For TI
  integer               :: ipartner, jpartner, ipartner_atm, jpartner_atm
  integer               :: partner_bond_idx
  logical               :: found_partner, write_msg
  double precision      :: req, partner_req
  double precision      :: mass_i
  
  if (ntc .eq. 2) then

    bond_cnt = nbonh

  else if (ntc .eq. 3) then

#if defined(MPI)
if(.not. usemidpoint) then
    if (master) write(mdout, *) 'this parallel version only works for ntc < 3'
    call mexit(6, 1)
endif
#endif

#ifdef CUDA
    write(mdout, *) 'CUDA version only works for ntc < 3'
    call mexit(6, 1)
#endif

    bond_cnt = nbonh + nbona

  else

    ! shake not done on anything.

    nonfastwat_bond_cnt = 0
    my_nonfastwat_bond_cnt = 0

    return

  end if

  ! Determine how many nonfastwater bonds there are.

  nonfastwat_bond_cnt = 0

  do bond_idx = 1, bond_cnt
    skip = .false.
    if (fastwat_bonds(bond_idx) .eq. 0) then
      ! Remove shake from bonds to softcore atoms.
      i = gbl_bond(bond_idx)%atm_i
      j = gbl_bond(bond_idx)%atm_j
      if (ti_mode .ne. 0 .and. ti_mode .ne. 1) then
        ! if this bond includes a softcore and non softcore atom
        if (ti_sc_lst(i)+ti_sc_lst(j) .eq. 2 .and. tishake .eq. 1) then
          if (ti_lst(1,i)+ti_lst(1,j) .ne. 0) then
            ti_num_noshake(1) = ti_num_noshake(1) + 1
            sc_num_noshake(1) = sc_num_noshake(1) + 1
          else
            ti_num_noshake(2) = ti_num_noshake(2) + 1          
            sc_num_noshake(2) = sc_num_noshake(2) + 1
          end if
          skip = .true.
        end if
      end if
      if (.not. skip .and. noshakemask_cnt .gt. 0) then
        if (noshakemask_lst(i)+noshakemask_lst(j) .ge. 1) then
          if (ti_mode .eq. 0 .or. ti_mode .eq. 1) then
            num_noshake = num_noshake + 1
            skip = .true.
          else            
            if (ti_lst(1,i)+ti_lst(1,j) .ne. 0) then
              ti_num_noshake(1) = ti_num_noshake(1) + 1
              if (ti_sc_lst(i)+ti_sc_lst(j) .ne. 0) then
                sc_num_noshake(1) = sc_num_noshake(1) + 1              
              end if
            else if (ti_lst(2,i)+ti_lst(2,j) .ne. 0) then
              ti_num_noshake(2) = ti_num_noshake(2) + 1
              if (ti_sc_lst(i)+ti_sc_lst(j) .ne. 0) then         
                sc_num_noshake(2) = sc_num_noshake(2) + 1
              end if
            else                       
              num_noshake = num_noshake + 1
            end if
            skip = .true.
          end if
        end if
      end if
      if (skip) then
        if (master) then
          do resi = 1, nres - 1
            if (i .ge. gbl_res_atms(resi) .and. &
                i .lt. gbl_res_atms(resi+1)) exit
          end do

          do resj = 1, nres - 1
            if (j .ge. gbl_res_atms(resj) .and. &
                j .lt. gbl_res_atms(resj+1)) exit
          end do

          write(mdout,'(a,a,a,a,i3,a,a,a,a,i3)') &
                  '   Removing shake constraints from ', &
                   atm_igraph(i),'  ',gbl_labres(resi),resi,' -- ',&
                   atm_igraph(j),'  ',gbl_labres(resj),resj
        end if
        cycle
      end if
      nonfastwat_bond_cnt = nonfastwat_bond_cnt + 1
    end if
  end do

#if defined(MPI) && !defined(CUDA)
if(.not. usemidpoint) then
  ! If you need it, set up the nonfastwater bonds list.

  if (imin .eq. 0) then         ! Molecular Dynamics

    ! Set up the nonfastwater bonds list.  This is only used for MPI
    ! molecular dynamics, not minimizations.

    allocate(nonfastwat_bond_lst(nonfastwat_bond_cnt), stat = alloc_failed)
    if (alloc_failed .ne. 0) call setup_alloc_error

    num_ints = num_ints + size(nonfastwat_bond_lst)

    bondlst_idx = 0

    do bond_idx = 1, bond_cnt
      if (fastwat_bonds(bond_idx) .eq. 0) then
        ! Remove shake from bonds to softcore atoms.
        i = gbl_bond(bond_idx)%atm_i
        j = gbl_bond(bond_idx)%atm_j
        if (ti_mode .ne. 0 .and. ti_mode .ne. 1) then
          ! if this bond includes a softcore and non softcore atom
          if (ti_sc_lst(i)+ti_sc_lst(j) .eq. 2 .and. tishake .eq. 1) cycle
        end if
        if (noshakemask_cnt .gt. 0) then
          if (noshakemask_lst(i)+noshakemask_lst(j) .ge. 1) then
            cycle
          end if
        end if
        bondlst_idx = bondlst_idx + 1
        nonfastwat_bond_lst(bondlst_idx) = bond_idx
      end if
    end do

  end if
end if ! not usemidpoint
#endif

  ! Claim all the nonfastwater bonds.  This is not true for mpi molecular
  ! dynamics, and will be corrected in claim_nonfastwater_bonds.

  my_nonfastwat_bond_cnt = nonfastwat_bond_cnt


#if defined(MPI) && !defined(CUDA)
  if (imin .eq. 0 .and. .not. usemidpoint) then         ! Molecular Dynamics
    call  claim_my_nonfastwat_bonds(num_ints, num_reals)
  else
#endif

    ! Allocate space.
    allocate(my_nonfastwat_bond_dat(nonfastwat_bond_cnt), stat = alloc_failed)
    if (alloc_failed .ne. 0) call setup_alloc_error

    num_ints = num_ints + 2 * size(my_nonfastwat_bond_dat)
    num_reals = num_reals + size(my_nonfastwat_bond_dat)

    bonddat_idx = 0

    do bond_idx = 1, bond_cnt
      if (fastwat_bonds(bond_idx) .eq. 0) then
        i = gbl_bond(bond_idx)%atm_i
        j = gbl_bond(bond_idx)%atm_j
        if (ti_mode .ne. 0 .and. ti_mode .ne. 1) then
          ! if this bond includes a softcore and non softcore atom
          if (ti_sc_lst(i)+ti_sc_lst(j) .eq. 2 .and. tishake .eq. 1) cycle
        end if
        if (noshakemask_cnt .gt. 0) then
          if (noshakemask_lst(i)+noshakemask_lst(j) .ge. 1) then
            cycle
          end if
        end if
        mass_i = 0.0
        if (ti_mode .ne. 0 .and. ti_mode .ne. 1 .and. tishake .eq. 2) then
          if(ti_sc_lst(i)+ti_sc_lst(j) .eq. 2) then
            if(ti_sc_lst(i) .eq. 2) then
              call ti_get_partner(j, new_atm)
              if(new_atm .ne. -1 .and. ti_lst(1,new_atm) .eq. 1) then
#ifdef CUDA
                if(atm_mass(j) .gt. atm_mass(i)) then
                    mass_i = atm_mass(j)
                    j=new_atm
                    i=i*(-1)
                else
                    mass_i = atm_mass(i)
                    j=new_atm*(-1)
                end if
#else
                j=new_atm
#endif
              end if
            else
              call ti_get_partner(i, new_atm)
              if(new_atm .ne. -1 .and. ti_lst(1,new_atm) .eq. 1) then
#ifdef CUDA
                if(atm_mass(j) .gt. atm_mass(i)) then
                    mass_i = atm_mass(j)
                    i=new_atm*(-1)
                else
                    mass_i = atm_mass(i)
                    j=j*(-1)
                    i=new_atm
                end if
#else
                i=new_atm
#endif
              end if
            end if
          end if
        end if
        bonddat_idx = bonddat_idx + 1
        if(tishake .eq. 2) then
          my_nonfastwat_bond_dat(bonddat_idx)%atm_i = i
          my_nonfastwat_bond_dat(bonddat_idx)%atm_j = j
        else
          my_nonfastwat_bond_dat(bonddat_idx)%atm_i = gbl_bond(bond_idx)%atm_i
          my_nonfastwat_bond_dat(bonddat_idx)%atm_j = gbl_bond(bond_idx)%atm_j
        end if
        my_nonfastwat_bond_dat(bonddat_idx)%parm = &
        gbl_req(gbl_bond(bond_idx)%parm_idx)**2
        if(mass_i .eq. 0.0 .and. tishake .eq. 2) then
          if(atm_mass(gbl_bond(bond_idx)%atm_i) .gt. atm_mass(gbl_bond(bond_idx)%atm_j)) then
            mass_i = atm_mass(gbl_bond(bond_idx)%atm_i)
          else
            mass_i = atm_mass(gbl_bond(bond_idx)%atm_j)
          end if
        end if
        my_nonfastwat_bond_dat(bonddat_idx)%mass_i = mass_i
      end if
    end do

#if defined(MPI) && !defined(CUDA)
  end if 
#endif

if(.not. usemidpoint) then
  ! Additional checks for SHAKE and TI
  ! This print a message if the equilibrium bond length between
  ! common atoms that are shaken is large. This can happen if 
  ! a hydrogen atom is changed into a heavy atom. 
  if (ti_mode .ne. 0) then
    write_msg = .false.    
    do bond_idx = 1, nbonh + nbona
      found_partner = .false.
      i = gbl_bond(bond_idx)%atm_i
      j = gbl_bond(bond_idx)%atm_j  
      req = gbl_req(gbl_bond(bond_idx)%parm_idx)
      
      ! No need to check non-TI atoms
      if (ti_lst(1,i)+ti_lst(1,j) .eq. 0 .and. &
          ti_lst(2,i)+ti_lst(2,j) .eq. 0) cycle
          
      ! Check for SHAKEn bond crossing into softcore region
      if (ti_sc_lst(i)+ti_sc_lst(j) .eq. 2 .and. &
          bond_idx .le. nbonh .and. tishake .ne. 1) then
        if (master) write(mdout,1113) i,j
        write_msg = .true.
        cycle
      end if
      
      ! Only check non-softcore TI atoms
      if (ti_sc_lst(i)+ti_sc_lst(j) .ne. 0) cycle
            
      call ti_get_partner(i, ipartner)
      call ti_get_partner(j, jpartner)
        
      ! If partner not found skip
      if (ipartner .lt. 0 .or. jpartner .lt. 0) cycle
        
      ! Find partner bond and get req
      do partner_bond_idx = 1, nbonh + nbona 
        ipartner_atm = gbl_bond(partner_bond_idx)%atm_i
        jpartner_atm = gbl_bond(partner_bond_idx)%atm_j
        if ((ipartner .eq. ipartner_atm .and. jpartner .eq. jpartner_atm) .or. &
            (jpartner .eq. ipartner_atm .and. ipartner .eq. jpartner_atm)) then
          partner_req = gbl_req(gbl_bond(partner_bond_idx)%parm_idx)
          found_partner = .true.
          exit
        end if
      end do
        
      ! If both bonds are not shaken then no need to check
      if (bond_idx .gt. nbonh .and. partner_bond_idx .gt. nbonh) then
        found_partner = .false.
      end if
        
      ! If both bonds are in noshake mask no need to check
      if (noshakemask_cnt .gt. 0) then
        if (noshakemask_lst(i)+noshakemask_lst(j) .eq. 0 .and. &
           noshakemask_lst(ipartner)+noshakemask_lst(jpartner) .eq. 0) then
        else
          found_partner = .false. ! skip if in noshake mask
        end if 
      end if 
        
      if (found_partner) then              
        if (abs(partner_req - req) .gt. 0.1d0) then             
           if (master) write (mdout, 1111) abs(partner_req - req), i,j,ipartner,jpartner
           write_msg = .true.
        end if                    
      end if 
    end do
    if (master .and. write_msg) write(mdout, 1112)
  end if
endif ! not usemidpoint
  
  1111 format(1x,'WARNING: Large deviation of ',f4.3, & 
      ' detected in equilibrium bond length for',/5x,'SHAKEn TI atoms', &
      1x,i8,2x,i8,2x,' and ',i8,2x,i8)
  1112 format(1x,'Coordinates for common atoms will be synchronized from V0. ',&
      /1x,'Consider using noshake mask or disabling SHAKE.',/)
  1113 format(1x,'WARNING: SHAKEn bond for atoms ', i8,1x,i8, ' is partially softcore.')
  return

end subroutine nonfastwat_shake_setup

end subroutine shake_setup

!*******************************************************************************
!
! Subroutine: shake_noshakemask
!
! Description: Init the noshakemask to prevent certain atoms from being shaken.
! 
!*******************************************************************************

subroutine shake_noshakemask(atm_cnt, nres, igraph, isymbl, res_atms, labres, &
                             crd, maskstr, num_ints, num_reals)
  use findmask_mod
  use file_io_dat_mod
  use pmemd_lib_mod
  use parallel_dat_mod
  use mdin_ctrl_dat_mod, only : usemidpoint

  implicit none

! Formal arguments:

  integer, intent(in)              :: atm_cnt
  integer, intent(in)              :: nres
  integer, intent(in)              :: res_atms(nres)
  character(len=4), intent(in)     :: igraph(atm_cnt)
  character(len=4), intent(in)     :: isymbl(atm_cnt)
  character(len=4), intent(in)     :: labres(nres)
  double precision, intent(in)     :: crd(3*atm_cnt)
  character(len=*), intent(inout)  :: maskstr
  integer                          :: num_ints, num_reals

! Local variables:
  integer               :: alloc_failed

  noshakemask_cnt = 0
  num_noshake = 0

  if (len_trim(maskstr) .eq. 0) return

  allocate(noshakemask_lst(atm_cnt), stat = alloc_failed)
  if (alloc_failed .ne. 0) call setup_alloc_error 
  num_ints = num_ints + size(noshakemask_lst)

  noshakemask_lst(:) = 0

if(.not. usemidpoint) &
  call atommask(atm_cnt, nres, 0, igraph, isymbl, res_atms, labres, crd, &
                maskstr, noshakemask_lst)

  noshakemask_cnt = sum(noshakemask_lst)
  if (master) then
     write (mdout, '(a,a,a,i7,a)') ' Noshake mask ',trim(maskstr),&
              ' matches ', noshakemask_cnt, ' atoms.'
  end if
  return

end subroutine shake_noshakemask


#ifdef MPI
!*******************************************************************************
!
! Subroutine:  shake_bcast
!
! Description: Broadcast noshakemask to all nodes
!              
!*******************************************************************************

subroutine bcast_shake(atm_cnt)
  use parallel_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:
  integer, intent(in)           :: atm_cnt
! Local variables:
  integer                       :: num_ints, num_reals !returned value discarded
  integer                       :: alloc_failed


  call mpi_bcast(noshakemask_cnt, 1, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  num_noshake = 0

  if (noshakemask_cnt .eq. 0) return

  if (.not. master) then
    num_ints = 0
    num_reals = 0
    allocate(noshakemask_lst(atm_cnt), stat = alloc_failed)
    if (alloc_failed .ne. 0) call setup_alloc_error 
  end if

  call mpi_bcast(noshakemask_lst, atm_cnt, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  return

end subroutine bcast_shake
#endif

!*******************************************************************************
!
! Subroutine:   get_water_distances
!
! Description:
!
! This routine determines the target lengths for the water molecules and
! derives from them the paramters required by fast 3-point water routine
! shake_fastwater. 
!
! Author: David A. Pearlman
! Date: 10/93
!
! Input:
!
! atm_igraph(i): The name of atom i (i4)
! nres: The number of residues
! gbl_labres(i): The name of residue i
! nbonh: The number of bonds to hydrogens
! nbona: The number of bonds to non-hydrogens
! gbl_bond(i)%atm_i
! gbl_bond(i)%atm_j: The two atoms of bond i.
! iwtnm: The name of residues to be considered water (i)
! iowtnm: The name of oxygen atom in the water residues (i)
! ihwtnm(2): The name of 2 hydrogen atoms in the water residues (i)
! jfastw: = 0: Use fast water routine (and default water names)
!              Use partitioning for fast TIP3P-TIP3P interactions in
!              nonbon routine (SANDER only).
!           4: Do not use fast waters anywhere.
! gbl_bond(i)%parm_idx: Parameter pointer for bond number i.
! req(i): req(icb(i)) is the target bond length of bond i.
! winv(i): The inverse mass of atom i.
!
! Input variables in global memory:
!
! gbl_res_atms(i): gbl_res_atms(i)->gbl_res_atms(i+1)-1 are atoms of residue i
!
! Output:
!
! rbtarg(8): ra, rb, rc, rc2, hhhh mass(o), mass(h), and mass (hoh)
! as required by the shake_fastwater routine.
!              
!*******************************************************************************

subroutine get_water_distances

  use constraints_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use pmemd_lib_mod

  implicit none

! Local variables:

  integer               :: atm_i, atm_j
  integer               :: i, j
  integer               :: atm_id
  integer               :: ihf1, ihf2, iof
  integer               :: nofast
  integer               :: parm_idx
  double precision      :: rhh, roh, roh1, roh2

  double precision, parameter   :: small = 1.d-4

! Search through the list for a water molecule. If one is found,
! determine the O-H and H-H distances ascribed to it.

  nofast = 1
  ! Following initialization guards against situation where this code gets
  ! called, but there is no recognized "water" in the system:
  iof = 0
  ihf1 = 0
  ihf2 = 0

  
outer: &
  do i = 1, res_cnt

    if (gbl_labres(i) .ne. iwtnm) cycle outer
    ! This initialization insures all values used come from the same water:
    iof = 0
    ihf1 = 0
    ihf2 = 0

    do j = 1, 3
      atm_id = gbl_res_atms(i) + j - 1

! Any waters not part of moving belly (atm_igroup(i)=0) will have had their bond
! parameters removed from list already and these must be skipped here:

      if (ibelly .gt. 0) then
        if (atm_igroup(atm_id) .le. 0) cycle outer
      end if

      if (atm_igraph(atm_id) .eq. iowtnm) then
        iof = atm_id
        rbtarg(6) = atm_mass(atm_id)
      else if (atm_igraph(atm_id) .eq. ihwtnm(1)) then
        ihf1 = atm_id
        rbtarg(7) = atm_mass(atm_id)
      else if (atm_igraph(atm_id) .eq. ihwtnm(2)) then
        ihf2 = atm_id
        rbtarg(7) = atm_mass(atm_id)
      else
        cycle outer
      end if
    end do

    if (iof .gt. 0 .and. ihf1 .gt. 0 .and. ihf2 .gt. 0) then
      nofast = 0
      exit outer
    end if

  end do outer

! We have found the three atoms. Now search the bonds list for the
! corresponding bonds. At that point, we can determine the target distances
! from gbl_req:

  roh1 = -10.0d0
  roh2 = -10.0d0
  rhh  = -10.0d0

  do j = 1, nbonh + nbona

    atm_i = gbl_bond(j)%atm_i
    atm_j = gbl_bond(j)%atm_j
    parm_idx = gbl_bond(j)%parm_idx

    if (atm_i .eq. iof .and. atm_j .eq. ihf1) then
      roh1 = gbl_req(parm_idx)
    else if (atm_j .eq. iof .and. atm_i .eq. ihf1) then
      roh1 = gbl_req(parm_idx)
    else if (atm_i .eq. iof .and. atm_j .eq. ihf2) then
      roh2 = gbl_req(parm_idx)
    else if (atm_j .eq. iof .and. atm_i .eq. ihf2) then
      roh2 = gbl_req(parm_idx)
    else if (atm_i .eq. ihf1.and. atm_j .eq. ihf2) then
      rhh  = gbl_req(parm_idx)
    else if (atm_j .eq. ihf1.and. atm_i .eq. ihf2) then
      rhh  = gbl_req(parm_idx)
    end if

  end do

! If all three bond lengths were not assigned in list, or if roh1 and
! roh2 are not the same length, assume something is a bit cockeyed,
! and exit with a message.

  if (roh1 .lt. 0.0d0 .or. roh2 .lt. 0.0d0 .or. rhh .lt. 0.0d0) then
    if (nofast .eq. 0 .and. master) then
       write(mdout, 1001)
       call mexit(6, 1)
    endif
  else if (abs(roh1 - roh2) .gt. small) then
    if (master) write(mdout, 1002)
    call mexit(6, 1)
  else
    roh = roh1
  end if

  rbtarg(8) = rbtarg(6) + 2.0d0 * rbtarg(7)

! setlep_setup sets the appropriate constants based on the bond lengths:

  call setlep_setup(rhh, roh)

#ifdef CUDA
  call gpu_get_water_distances(rbtarg)
#endif

  return

1001 format('ERROR: Bond lengths params not found in PARM file', &
            ' for fast water model')
1002 format('ERROR: R(O-H) bond lengths params found in PARM file', &
            ' for fast water model', /, t10, 'not identical')

end subroutine get_water_distances

!*******************************************************************************
!
! Internal Subroutine:  setlep_setup  (SETLep setup) 
!
! Description:
!
! This routine determines the values of ra, rb, rc, rc2, and hhhh
! required by the setlep program. These are calculated from the
! target rigid water distances rhh (h...h) and roh (o...h).
!
! The calulated values are returned in rbtarg(1-5):
!
!    rbtarg(1) ... ra
!    rbtarg(2) ... rb
!    rbtarg(3) ... rc
!    rbtarg(4) ... rc2
!    rbtarg(5) ... hhhh
!    (rbtarg(6->8) are masses already set in calling routine)
!
! Author: David A. Pearlman
! Date: 10/93
!              
!*******************************************************************************

subroutine setlep_setup(rhh, roh)

  implicit none

! Formal arguments:

  double precision  rhh, roh

! Local variables:

  double precision  comx, comy
  double precision  dis
  double precision  height
  double precision  wh, wo, wohh
  double precision  x(2, 3)

  double precision, parameter   :: zero = 0.d0, two = 2.d0

  wo = rbtarg(6)
  wh = rbtarg(7)
  wohh = rbtarg(8)

! Create triangle in x,y plane with base parallel to x:

  height = sqrt(roh**2 - (rhh/two)**2)
  x(1, 1) = -rhh/two
  x(2, 1) = -height
  x(1, 2) =  rhh/two
  x(2, 2) = -height
  x(1, 3) = zero
  x(2, 3) = zero

! Calculate the center of mass of the triangle:

  comx = (x(1, 1)*wh + x(1, 2)*wh + x(1, 3)*wo)/wohh
  comy = (x(2, 1)*wh + x(2, 2)*wh + x(2, 3)*wo)/wohh

! The distance between the center of mass and the apex is ra. rb is
! the height - ra.

  dis = sqrt(comx * comx + comy * comy)
  rbtarg(1) = dis
  rbtarg(2) = height - rbtarg(1)
  rbtarg(3) = rhh / two
  rbtarg(4) = rhh
  rbtarg(5) = rbtarg(4)*rbtarg(4)

  return

end subroutine setlep_setup

!*******************************************************************************
!
! Subroutine:  shake_fastwater
!
! Description:
!
! This routine makes the calls to carry out the fast, analytic 3-point
! water constraints.
!
! Author: David A. Pearlman
! Date: 12/92
!
! INPUT:
! x0(i): Coordinate array corresponding to time t-dt/2. Coordinates
!        are packed such that x(atom i), y(atom i) and z(atom i) are
!        found at x0(3*(i-1)+1), x0(3*(i-1)+2) and x0(3*(i-1)+3), respectively.
!
! xh(i): Coordinate array corresponding to time t+dt/2, but not corrected
!        (on input) for internal constraints, i.e. an unconstrained move.
!        Packing is the same as for x0.
!
! fastwat_res(i): If = 0, the bonds of this residue are constrained by shake.
!                 If fastwat_res(i)=1, the bonds of the residue are constrained
!                 by the fast, analytic 3-point routine.
!
! gbl_res_atms(i): gbl_res_atms(i): gbl_res_atms(i+1)-1 are atoms in residue i.
!
! iorwat: The position of the oxygen atom in each water (1,2, or 3).
!
! rbtarg(8): The five geometry parameters used by the fast water routine
!         to impose constraints, the masses of the atoms in the water,
!         and the total mass of the water.
!              
!*******************************************************************************
!*******************************************************************************
!
! Note - This subroutine has been subsumed into shake_fastwater, but the
!        header is retained for the documentation value.
!
! Internal Subroutine : setlep - reset positions of TIP3P waters
!
! Description:
!                                                               
!    Author : Shuichi Miyamoto                                  
!    Date of last update : Dec. 12, 1992                        
!                                                               
!    Reference for the SETTLE algorithm                         
!      S. Miyamoto et al., J. Comp. Chem.,  13, 952 (1992)      
!                                                               
!    Revisions:                                                 
!      12/92: Add iorwat to call list. Allows O atom to be at   
!             any position in the water. -- David Pearlman      
!                                                               
!      12/92: Modify routine to be consistent with the Amber    
!             data structures; merge x0,y0,z0 arrays into a     
!             single x0(3,mxatm) array and the x1,y1,z1 arrays  
!             into a single x1(3,mxatm) array.                  
!                                                               
!      10/93: Add rbtarg(8) to call list. This contains the     
!             values of the geometric parameters ra, rb, rc,    
!             rc2 and hhhh, and the masses mh, mo, and mhoh     
!             required for this routine. Previously these were  
!             hard-wired.                                       
!                                       -- David Pearlman       
!
! Most of the variable names correspond to those in the reference:
!
!     x0(3,i)  : position at previous time step (t0) (input)
!     x1(3,i)  : position at present  time step (t0 + dt)
!              : position before applying constraints (input)
!              : position after  applying constraints (output)
!     mxatm    : maximum No. of atoms in system
!     iorwat   : position of oxygen in water (1, 2, or 3)
!
!     x,y,zcom       : center of mass
!     x,y,zaksX(YZ)d : axis vectors of the alternative orthogonal
!                      coordinate system Xprime Yprime Zprime
!     trns..         : matrix of orthogonal transformation
!
!     wo,wh : mass of oxygen and hydrogen
!     wohh  : mass of water (H2O)
!     hhhh  : rc2*rc2 (square of HH distance)
!
!*******************************************************************************

subroutine shake_fastwater(x0, x1)

  use pbc_mod
  use processor_mod
  use mdin_ctrl_dat_mod, only : usemidpoint
  implicit none

! Formal arguments:

  double precision      :: x0(3, atm_cnt)
  double precision      :: x1(3, atm_cnt)

! Local variables:

  double precision      :: alpa, beta, gama
  double precision      :: al2be2
  double precision      :: axlng_inv, aylng_inv, azlng_inv
  double precision      :: cosphi, cospsi, costhe
! double precision      :: deltx
! double precision      :: hh2, hhhh
  double precision      :: hhhh
  double precision      :: ra, ra_inv, rb, rc
  double precision      :: rc2
  double precision      :: sinphi, sinpsi, sinthe
  double precision      :: trns11, trns12, trns13
  double precision      :: trns21, trns22, trns23
  double precision      :: trns31, trns32, trns33
  double precision      :: wh_div_wohh, wo_div_wohh
  double precision      :: xa1, ya1, za1
  double precision      :: xa3d, ya3d, za3d
  double precision      :: xaksxd, yaksxd, zaksxd
  double precision      :: xaksyd, yaksyd, zaksyd
  double precision      :: xakszd, yakszd, zakszd
  double precision      :: xb0d, yb0d
  double precision      :: xb1, yb1, zb1
  double precision      :: xb1d, yb1d, zb1d
  double precision      :: xb2d, yb2d
  double precision      :: xb3d, yb3d, zb3d
  double precision      :: xb2d2
  double precision      :: xc0d, yc0d
  double precision      :: xc1, yc1, zc1
  double precision      :: xc1d, yc1d, zc1d
  double precision      :: xc3d, yc3d, zc3d
  double precision      :: xb0, yb0, zb0
  double precision      :: xc0, yc0, zc0
  double precision      :: xcom, ycom, zcom
  double precision      :: ya2d
  double precision      :: yc2d
  double precision      :: za1d
  double precision      :: x1j, x2j, x3j, x1k, x2k, x3k

  integer               :: atm_1, atm_2, atm_3
  integer               :: first_res_atm
  integer               :: ind1, ind2, ind3
  integer               :: res_idx

#if defined(MPI) && !defined(CUDA)
  if (my_fastwat_res_cnt .eq. 0) return
#endif

  ra = rbtarg(1)
  ra_inv = 1.d0 / ra
  rb = rbtarg(2)
  rc = rbtarg(3)
  rc2 = rbtarg(4)
  hhhh = rbtarg(5)
  wo_div_wohh = rbtarg(6) / rbtarg(8)
  wh_div_wohh = rbtarg(7) / rbtarg(8)

  if (iorwat .eq. 1) then
    ind1 = 0
    ind2 = 1
    ind3 = 2
  else if (iorwat .eq. 2) then
    ind1 = 1
    ind2 = 2
    ind3 = 0
  else
    ind1 = 2
    ind2 = 0
    ind3 = 1
  end if

if(usemidpoint) then
  select_my_fastwat_res_cnt=my_fastwat_res_cnt_midpoint 
else
  select_my_fastwat_res_cnt = my_fastwat_res_cnt
endif

  do res_idx = 1, select_my_fastwat_res_cnt
if(usemidpoint) then
#ifdef MPI
    first_res_atm = my_fastwat_res_lst_midpoint(res_idx)

    atm_1 = proc_shake_space(first_res_atm + ind1)
    atm_2 = proc_shake_space(first_res_atm + ind2)
    atm_3 = proc_shake_space(first_res_atm + ind3)

    xb0 = x0(1, atm_2) + mult_vec_wat(1,3*res_idx+ind2) * pbc_box(1) - x0(1, atm_1)
    yb0 = x0(2, atm_2) + mult_vec_wat(2,3*res_idx+ind2) * pbc_box(2) - x0(2, atm_1)
    zb0 = x0(3, atm_2) + mult_vec_wat(3,3*res_idx+ind2) * pbc_box(3) - x0(3, atm_1)
    xc0 = x0(1, atm_3) + mult_vec_wat(1,3*res_idx+ind3) * pbc_box(1) - x0(1, atm_1)
    yc0 = x0(2, atm_3) + mult_vec_wat(2,3*res_idx+ind3) * pbc_box(2) - x0(2, atm_1)
    zc0 = x0(3, atm_3) + mult_vec_wat(3,3*res_idx+ind3) * pbc_box(3) - x0(3, atm_1)
    
    x1j = x1(1, atm_2) + mult_vec_wat(1,3*res_idx+ind2) * pbc_box(1)
    x2j = x1(2, atm_2) + mult_vec_wat(2,3*res_idx+ind2) * pbc_box(2)
    x3j = x1(3, atm_2) + mult_vec_wat(3,3*res_idx+ind2) * pbc_box(3)
    x1k = x1(1, atm_3) + mult_vec_wat(1,3*res_idx+ind3) * pbc_box(1)
    x2k = x1(2, atm_3) + mult_vec_wat(2,3*res_idx+ind3) * pbc_box(2)
    x3k = x1(3, atm_3) + mult_vec_wat(3,3*res_idx+ind3) * pbc_box(3)

    xcom = x1(1, atm_1) * wo_div_wohh + (x1j + &
           x1k) * wh_div_wohh
    ycom = x1(2, atm_1) * wo_div_wohh + (x2j + &
           x2k) * wh_div_wohh
    zcom = x1(3, atm_1) * wo_div_wohh + (x3j + &
           x3k) * wh_div_wohh
    xa1 = x1(1, atm_1) - xcom
    ya1 = x1(2, atm_1) - ycom
    za1 = x1(3, atm_1) - zcom
    xb1 = x1j - xcom
    yb1 = x2j - ycom
    zb1 = x3j - zcom
    xc1 = x1k - xcom
    yc1 = x2k - ycom
    zc1 = x3k - zcom
#endif
else ! usemidpoint
    first_res_atm = my_fastwat_res_lst(res_idx)

    atm_1 = first_res_atm + ind1
    atm_2 = first_res_atm + ind2
    atm_3 = first_res_atm + ind3

    xb0 = x0(1, atm_2) - x0(1, atm_1)
    yb0 = x0(2, atm_2) - x0(2, atm_1)
    zb0 = x0(3, atm_2) - x0(3, atm_1)
    xc0 = x0(1, atm_3) - x0(1, atm_1)
    yc0 = x0(2, atm_3) - x0(2, atm_1)
    zc0 = x0(3, atm_3) - x0(3, atm_1)

    xcom = x1(1, atm_1) * wo_div_wohh + (x1(1, atm_2) + &
           x1(1, atm_3)) * wh_div_wohh
    ycom = x1(2, atm_1) * wo_div_wohh + (x1(2, atm_2) + &
           x1(2, atm_3)) * wh_div_wohh
    zcom = x1(3, atm_1) * wo_div_wohh + (x1(3, atm_2) + &
           x1(3, atm_3)) * wh_div_wohh

    xa1 = x1(1, atm_1) - xcom
    ya1 = x1(2, atm_1) - ycom
    za1 = x1(3, atm_1) - zcom
    xb1 = x1(1, atm_2) - xcom
    yb1 = x1(2, atm_2) - ycom
    zb1 = x1(3, atm_2) - zcom
    xc1 = x1(1, atm_3) - xcom
    yc1 = x1(2, atm_3) - ycom
    zc1 = x1(3, atm_3) - zcom
endif ! usemidpoint

! Step1  A1_prime:

    xakszd = yb0 * zc0 - zb0 * yc0
    yakszd = zb0 * xc0 - xb0 * zc0
    zakszd = xb0 * yc0 - yb0 * xc0
    xaksxd = ya1 * zakszd - za1 * yakszd
    yaksxd = za1 * xakszd - xa1 * zakszd
    zaksxd = xa1 * yakszd - ya1 * xakszd
    xaksyd = yakszd * zaksxd - zakszd * yaksxd
    yaksyd = zakszd * xaksxd - xakszd * zaksxd
    zaksyd = xakszd * yaksxd - yakszd * xaksxd

    axlng_inv = 1.d0 / sqrt(xaksxd * xaksxd + yaksxd * yaksxd + zaksxd * zaksxd)
    aylng_inv = 1.d0 / sqrt(xaksyd * xaksyd + yaksyd * yaksyd + zaksyd * zaksyd)
    azlng_inv = 1.d0 / sqrt(xakszd * xakszd + yakszd * yakszd + zakszd * zakszd)

    trns11 = xaksxd * axlng_inv
    trns21 = yaksxd * axlng_inv
    trns31 = zaksxd * axlng_inv
    trns12 = xaksyd * aylng_inv
    trns22 = yaksyd * aylng_inv
    trns32 = zaksyd * aylng_inv
    trns13 = xakszd * azlng_inv
    trns23 = yakszd * azlng_inv
    trns33 = zakszd * azlng_inv

    xb0d = trns11 * xb0 + trns21 * yb0 + trns31 * zb0
    yb0d = trns12 * xb0 + trns22 * yb0 + trns32 * zb0
    xc0d = trns11 * xc0 + trns21 * yc0 + trns31 * zc0
    yc0d = trns12 * xc0 + trns22 * yc0 + trns32 * zc0
    za1d = trns13 * xa1 + trns23 * ya1 + trns33 * za1
    xb1d = trns11 * xb1 + trns21 * yb1 + trns31 * zb1
    yb1d = trns12 * xb1 + trns22 * yb1 + trns32 * zb1
    zb1d = trns13 * xb1 + trns23 * yb1 + trns33 * zb1
    xc1d = trns11 * xc1 + trns21 * yc1 + trns31 * zc1
    yc1d = trns12 * xc1 + trns22 * yc1 + trns32 * zc1
    zc1d = trns13 * xc1 + trns23 * yc1 + trns33 * zc1

! Step2  A2_prime:

    sinphi = za1d * ra_inv
    cosphi = sqrt(1.d0 - sinphi*sinphi)
    sinpsi = (zb1d - zc1d) / (rc2 * cosphi)
    cospsi = sqrt(1.d0 - sinpsi*sinpsi)
 
    ya2d =   ra * cosphi
    xb2d = - rc * cospsi
    yb2d = - rb * cosphi - rc *sinpsi * sinphi
    yc2d = - rb * cosphi + rc *sinpsi * sinphi
    xb2d2 = xb2d * xb2d
!   hh2 = 4.d0 * xb2d2 + (yb2d-yc2d) * (yb2d-yc2d) + (zb1d-zc1d) * (zb1d-zc1d)
!   deltx = 2.d0 * xb2d + sqrt(4.d0 * xb2d2 - hh2 + hhhh)
!   xb2d = xb2d - deltx * 0.5d0

    xb2d =  -0.5d0 * sqrt(hhhh - &
                          (yb2d-yc2d) * (yb2d-yc2d) - (zb1d-zc1d) * (zb1d-zc1d))

! Step3  al,be,ga:

    alpa = (xb2d * (xb0d-xc0d) + yb0d * yb2d + yc0d * yc2d)
    beta = (xb2d * (yc0d-yb0d) + xb0d * yb2d + xc0d * yc2d)
    gama = xb0d * yb1d - xb1d * yb0d + xc0d * yc1d - xc1d * yc0d

    al2be2 = alpa * alpa + beta * beta
    sinthe = (alpa*gama - beta * sqrt(al2be2 - gama * gama)) / al2be2

! Step4  A3_prime:

    costhe = sqrt(1.d0 - sinthe * sinthe)
    xa3d = -ya2d * sinthe
    ya3d = ya2d * costhe
    za3d = za1d
    xb3d = xb2d * costhe - yb2d * sinthe
    yb3d = xb2d * sinthe + yb2d * costhe
    zb3d = zb1d
    xc3d = -xb2d * costhe - yc2d * sinthe
    yc3d = -xb2d * sinthe + yc2d * costhe
    zc3d = zc1d

! Step5  A3:

    x1(1, atm_1) = xcom + trns11 * xa3d + trns12 * ya3d + trns13 * za3d
    x1(2, atm_1) = ycom + trns21 * xa3d + trns22 * ya3d + trns23 * za3d
    x1(3, atm_1) = zcom + trns31 * xa3d + trns32 * ya3d + trns33 * za3d
if(usemidpoint) then
    x1(1, atm_2) = xcom + trns11 * xb3d + trns12 * yb3d + trns13 * zb3d - &
                   mult_vec_wat(1,3*res_idx+ind2) * pbc_box(1)
    x1(2, atm_2) = ycom + trns21 * xb3d + trns22 * yb3d + trns23 * zb3d - &
                   mult_vec_wat(2,3*res_idx+ind2) * pbc_box(2)
    x1(3, atm_2) = zcom + trns31 * xb3d + trns32 * yb3d + trns33 * zb3d - &
                   mult_vec_wat(3,3*res_idx+ind2) * pbc_box(3)
    x1(1, atm_3) = xcom + trns11 * xc3d + trns12 * yc3d + trns13 * zc3d - &
                   mult_vec_wat(1,3*res_idx+ind3) * pbc_box(1)
    x1(2, atm_3) = ycom + trns21 * xc3d + trns22 * yc3d + trns23 * zc3d - &
                   mult_vec_wat(2,3*res_idx+ind3) * pbc_box(2)
    x1(3, atm_3) = zcom + trns31 * xc3d + trns32 * yc3d + trns33 * zc3d - &
                   mult_vec_wat(3,3*res_idx+ind3) * pbc_box(3)
else ! usemidpoint
    x1(1, atm_2) = xcom + trns11 * xb3d + trns12 * yb3d + trns13 * zb3d
    x1(2, atm_2) = ycom + trns21 * xb3d + trns22 * yb3d + trns23 * zb3d
    x1(3, atm_2) = zcom + trns31 * xb3d + trns32 * yb3d + trns33 * zb3d
    x1(1, atm_3) = xcom + trns11 * xc3d + trns12 * yc3d + trns13 * zc3d
    x1(2, atm_3) = ycom + trns21 * xc3d + trns22 * yc3d + trns23 * zc3d
    x1(3, atm_3) = zcom + trns31 * xc3d + trns32 * yc3d + trns33 * zc3d
endif ! usemidpoint

  end do

  return

end subroutine shake_fastwater

!******************************************************************************
!
! Subroutine : rattle_fastwater
!
! Description:
!
!
!******************************************************************************
subroutine rattle_fastwater(x1, v1)

    use pbc_mod
    use processor_mod
    use mdin_ctrl_dat_mod, only : usemidpoint
    implicit none

! Formal arguments:
    
    double precision      :: x1(3, atm_cnt)
    double precision      :: v1(3, atm_cnt)

! Local variables:
    
    double precision      :: wh, wo, woh, whh, woh2, wowh2, whwh
    double precision      :: wohwoh, wohwo, whhwh

    double precision      :: ablng, bclng, calng
    double precision      :: abmc, bcma, camb
    double precision      :: cosa, cosb, cosc
    double precision      :: deno

    double precision      :: tabd, tbcd, tcad
    double precision      :: vabab, vbcbc, vcaca
    double precision      :: xab, yab, zab, xbc, ybc, zbc, xca, yca, zca
    double precision      :: xeab, yeab, zeab
    double precision      :: xebc, yebc, zebc
    double precision      :: xeca, yeca, zeca
    double precision      :: xvab, yvab, zvab
    double precision      :: xvbc, yvbc, zvbc
    double precision      :: xvca, yvca, zvca

    integer               :: atm_1, atm_2, atm_3
    integer               :: first_res_atm
    integer               :: ind1, ind2, ind3
    integer               :: res_idx


#if defined(MPI) && !defined(CUDA)
    if (my_fastwat_res_cnt .eq. 0) return
#endif



    wo = rbtarg(6)
    wh = rbtarg(7)


    woh = wo + wh
    whh = wh + wh
    woh2 = woh * 2.0d0
    wowh2 = wo * whh
    whwh = wh * wh
    wohwoh = woh * woh
    wohwo = woh * wo
    whhwh = whh * wh

    if (iorwat .eq. 1) then
        ind1 = 0
        ind2 = 1
        ind3 = 2
    else if (iorwat .eq. 2) then 
        ind1 = 1
        ind2 = 2
        ind3 = 0
    else 
        ind1 = 2
        ind2 = 0
        ind3 = 1
    end if

    if (usemidpoint) then
        select_my_fastwat_res_cnt = my_fastwat_res_cnt_midpoint
    else 
        select_my_fastwat_res_cnt = my_fastwat_res_cnt
    end if

    do res_idx = 1, select_my_fastwat_res_cnt

! Step1 AB, VAB
        if (usemidpoint) then
#ifdef MPI
            first_res_atm = my_fastwat_res_lst_midpoint(res_idx)

            atm_1 = proc_shake_space(first_res_atm + ind1)
            atm_2 = proc_shake_space(first_res_atm + ind2)
            atm_3 = proc_shake_space(first_res_atm + ind3)

            xab = x1(1, atm_2) + mult_vec_wat(1, 3 * res_idx + ind2) * pbc_box(1) - x1(1, atm_1)
            yab = x1(2, atm_2) + mult_vec_wat(2, 3 * res_idx + ind2) * pbc_box(2) - x1(2, atm_1)
            zab = x1(3, atm_2) + mult_vec_wat(3, 3 * res_idx + ind2) * pbc_box(3) - x1(3, atm_1)
            xbc = x1(1, atm_2) + mult_vec_wat(1, 3 * res_idx + ind1) * pbc_box(1) - x1(1, atm_3)
            ybc = x1(2, atm_2) + mult_vec_wat(2, 3 * res_idx + ind1) * pbc_box(2) - x1(2, atm_3)
            zbc = x1(3, atm_2) + mult_vec_wat(3, 3 * res_idx + ind1) * pbc_box(3) - x1(3, atm_3)
            xca = x1(1, atm_3) + mult_vec_wat(1, 3 * res_idx + ind3) * pbc_box(1) - x1(1, atm_1)
            yca = x1(2, atm_3) + mult_vec_wat(2, 3 * res_idx + ind3) * pbc_box(2) - x1(2, atm_1)
            zca = x1(3, atm_3) + mult_vec_wat(3, 3 * res_idx + ind3) * pbc_box(3) - x1(3, atm_1)
#endif
        else
            first_res_atm = my_fastwat_res_lst(res_idx)

            atm_1 = first_res_atm + ind1
            atm_2 = first_res_atm + ind2
            atm_3 = first_res_atm + ind3

            xab = x1(1, atm_2) - x1(1, atm_1)
            yab = x1(2, atm_2) - x1(2, atm_1)
            zab = x1(3, atm_2) - x1(3, atm_1)
            xbc = x1(1, atm_3) - x1(1, atm_2)
            ybc = x1(2, atm_3) - x1(2, atm_2)
            zbc = x1(3, atm_3) - x1(3, atm_2)
            xca = x1(1, atm_1) - x1(1, atm_3)
            yca = x1(2, atm_1) - x1(2, atm_3)
            zca = x1(3, atm_1) - x1(3, atm_3)
        end if

        xvab = v1(1, atm_2) - v1(1, atm_1)
        yvab = v1(2, atm_2) - v1(2, atm_1)
        zvab = v1(3, atm_2) - v1(3, atm_1)
        xvbc = v1(1, atm_3) - v1(1, atm_2)
        yvbc = v1(2, atm_3) - v1(2, atm_2)
        zvbc = v1(3, atm_3) - v1(3, atm_2)
        xvca = v1(1, atm_1) - v1(1, atm_3)
        yvca = v1(2, atm_1) - v1(2, atm_3)
        zvca = v1(3, atm_1) - v1(3, atm_3)

! Step2 eab
        ablng = sqrt(xab * xab + yab * yab + zab * zab)
        bclng = sqrt(xbc * xbc + ybc * ybc + zbc * zbc)
        calng = sqrt(xca * xca + yca * yca + zca * zca)

        xeab = xab / ablng
        yeab = yab / ablng
        zeab = zab / ablng
        xebc = xbc / bclng
        yebc = ybc / bclng
        zebc = zbc / bclng
        xeca = xca / calng
        yeca = yca / calng
        zeca = zca / calng

! Step3 vabab
        vabab = xvab * xeab + yvab * yeab + zvab * zeab
        vbcbc = xvbc * xebc + yvbc * yebc + zvbc * zebc
        vcaca = xvca * xeca + yvca * yeca + zvca * zeca

! Step4 tab
        cosa = - xeab * xeca - yeab * yeca - zeab * zeca
        cosb = - xebc * xeab - yebc * yeab - zebc * zeab
        cosc = - xeca * xebc - yeca * yebc - zeca * zebc
        abmc = wh * cosa * cosb - woh * cosc
        bcma = wo * cosb * cosc - whh * cosa
        camb = wh * cosc * cosa - woh * cosb
        tabd = vabab * (woh2 - wo * cosc * cosc) &
            + vbcbc * camb + vcaca * bcma
        tbcd = vbcbc * (wohwoh - whwh * cosa * cosa) &
            + vcaca * abmc * wo + vabab * camb * wo
        tcad = vcaca * (woh2 - wo * cosb * cosb) &
            + vabab * bcma + vbcbc * abmc
        deno = 2.d0 * wohwoh + wowh2 * cosa * cosb * cosc &
            - whhwh * cosa * cosa &
            - wohwo * (cosb * cosb + cosc * cosc)

! Step5 V

        v1(1, atm_1) = v1(1, atm_1) + (xeab * tabd - xeca * tcad) * wh / deno
        v1(2, atm_1) = v1(2, atm_1) + (yeab * tabd - yeca * tcad) * wh / deno
        v1(3, atm_1) = v1(3, atm_1) + (zeab * tabd - zeca * tcad) * wh / deno
        v1(1, atm_2) = v1(1, atm_2) + (xebc * tbcd - xeab * tabd * wo) / deno
        v1(2, atm_2) = v1(2, atm_2) + (yebc * tbcd - yeab * tabd * wo) / deno
        v1(3, atm_2) = v1(3, atm_2) + (zebc * tbcd - zeab * tabd * wo) / deno
        v1(1, atm_3) = v1(1, atm_3) + (xeca * tcad * wo - xebc * tbcd) / deno
        v1(2, atm_3) = v1(2, atm_3) + (yeca * tcad * wo - yebc * tbcd) / deno
        v1(3, atm_3) = v1(3, atm_3) + (zeca * tcad * wo - zebc * tbcd) / deno

    end do


    return

end subroutine rattle_fastwater


!*******************************************************************************
!
! Subroutine:  shake
!
! Description:
!
! All the bonds involving hydrogen atoms are loaded first in struct gbl_bond
! followed by those involving non-hydrogen atoms and the perturbed atoms.
!
! Mods for version 4.1:
!   - Add fastwat_bonds(i), so that waters shaken by fast analytic 3-point
!     shake will not be shake-n here -- dap
!
#if defined(MPI) && !defined(CUDA)
!
! The only MPI code specific to this routine is wrappers for
! all the I/O to allow only the master to write.
!
#endif
!
!*******************************************************************************

subroutine shake(x, xp)
#include "shake.i"
end subroutine shake

#ifdef _OPENMP_
subroutine shake_gb(x, xp)
#define GBorn
#include "shake.i"
#undef GBorn
end subroutine shake_gb
#endif /*_OPENMP_*/

subroutine rattle(x, v)
#include "rattle.i"
end subroutine rattle

#ifdef _OPENMP_
subroutine rattle_gb(x, v)
#define GBorn
#include "rattle.i"
#undef GBorn
end subroutine rattle_gb
#endif /*_OPENMP_*/

!*******************************************************************************
!
! Internal Subroutine:  wrap_crds
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine wrap_crds(xpij)

  use pbc_mod

  implicit none

! Formal arguments:

  double precision  xpij(3)

! Local variables:

  integer           m

  do m = 1, 3

    if (xpij(m) .ge. box_half(m)) then

      xpij(m) = xpij(m) - pbc_box(m)

    else if (xpij(m) .lt. -box_half(m)) then

      xpij(m) = xpij(m) + pbc_box(m)

    end if

  end do

  return

end subroutine wrap_crds

end module shake_mod
