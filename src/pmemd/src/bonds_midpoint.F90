#include "copyright.i"
#include "hybrid_datatypes.i"

!*******************************************************************************
!
! Module:  bonds_midpoint_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module bonds_midpoint_mod
#ifdef MPI

  use gbl_datatypes_mod
  use processor_mod

  implicit none

! The following are derived from prmtop bond info:

  integer, save                         :: cit_nbonh, cit_nbona

  type(bond_rec), allocatable, save     :: cit_h_bond(:) !array holding hydrogen bond
  type(bond_rec), allocatable, save     :: cit_a_bond(:) !array holding non-hydrogen bond
  type(bond_rec), allocatable, save     :: old_cit_h_bond(:) !array holding hydrogen bond
  type(bond_rec), allocatable, save     :: old_cit_a_bond(:) !array holding non-hydrogen bond
  integer, allocatable, save            :: my_hbonds_leads(:)
  integer, allocatable, save            :: my_abonds_leads(:)
  type(bond_trans), allocatable, save   :: b_recv_buf(:,:), b_send_buf(:,:)
  integer, save            :: mult_bonded 


contains

!*******************************************************************************
!
! Subroutine:  bonds_setup
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine build_atm_space()

  use parallel_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
  use processor_mod, only: proc_atm_space, proc_atm_space_ghosts, &
                           proc_num_atms, proc_ghost_num_atms, &
                           proc_atm_to_full_list

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  !integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: i, atmid 
  !Can we optimize so that we can call it only one time, and next on, just
  !update as new atoms comes and goes, and with sorting 
 proc_atm_space(:) = 0 !very expensive, need to only reset the previous nb steps
 proc_atm_space_ghosts(:) = 0!very expensive, need to only reset the previous nb steps
  do i = 1, proc_num_atms
         atmid = proc_atm_to_full_list(i)
         proc_atm_space(atmid) = i
  end do
  do i = proc_num_atms+1, proc_num_atms+proc_ghost_num_atms
         atmid = proc_atm_to_full_list(i)
         proc_atm_space(atmid) = i
  end do
end subroutine build_atm_space
subroutine bonds_midpoint_setup()

  use parallel_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
  use processor_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  !integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: i, atmid 
  integer               :: alloc_failed
  integer               :: alloc_len
  type(bond_rec)        :: bonds_copy(nbonh + nbona)
  mult_bonded = 10 ! Multiplier to allocate space for bonded interactions
!  print *,"Rank,min,max",mytaskid, proc_min_x_crd, proc_max_x_crd, &
!proc_min_y_crd, proc_max_y_crd,proc_min_z_crd, proc_max_z_crd
 
  ! Fill array : local to node global
  !it is already done inthe above routine, so we can skip it
#if 1
  proc_atm_space(:) =0
  proc_atm_space_ghosts(:) =0
  do i = 1, proc_num_atms
         atmid = proc_atm_to_full_list(i)
         proc_atm_space(atmid) = i
  end do
  do i = proc_num_atms+1, proc_num_atms+proc_ghost_num_atms
         atmid = proc_atm_to_full_list(i)
         proc_atm_space_ghosts(atmid) = i
  end do
#endif
!    if (mytaskid .eq. 0) then 
!print *, &
!"Bonded:proc_atm_space:",proc_atm_space(1000:1010)
!print *, "proc_num_atms:",proc_num_atms
!end if
  ! Scan Bonds for the sub-domain 
  if (nbonh .gt. 0) then
    call find_mysd_bonds(nbonh, gbl_bond, cit_nbonh, bonds_copy)
  end if
    !print *, "Rank,cit_nbonh:",mytaskid,cit_nbonh
  ! This routine can handle reallocation, and thus can be called multiple
  ! times.

  !call find_my_bonds(nbonh, gbl_bond, cit_nbonh, bonds_copy, use_atm_map)

  alloc_len = max(cit_nbonh,proc_num_atms_min_bound) * mult_bonded
  if (.not. allocated(cit_h_bond)) allocate(cit_h_bond(alloc_len), stat =alloc_failed)
  if (.not. allocated(old_cit_h_bond)) allocate(old_cit_h_bond(alloc_len), stat =alloc_failed)
  if (alloc_failed .ne. 0) call setup_alloc_error
  if (cit_nbonh .gt. 0) then
    cit_h_bond(1:cit_nbonh) = bonds_copy(1:cit_nbonh)
  end if
    !print *,"Rank,size cit_h_bond",mytaskid, size(cit_h_bond)
  if(.not. allocated(my_hbonds_leads))  allocate(my_hbonds_leads(alloc_len), stat = alloc_failed)
  if (alloc_failed .ne. 0) call setup_alloc_error
  if (cit_nbonh .gt. 0) then
    do i = 1, cit_nbonh
      my_hbonds_leads(i) = cit_h_bond(i)%atm_i
    end do 
  end if 
  ! Create an array (my_hbonds_leads) of lead atoms that is useful after neighbor list 
  ! building to check the bond(with lead atom) which are to be transfred.
 
 
  ! Bonds without Hydorgen atoms: 
  if (nbona .gt. 0) then  
    call find_mysd_bonds(nbona, gbl_bond(bonda_idx), cit_nbona, bonds_copy)
  end if
  !call find_my_bonds(nbona, gbl_bond(bonda_idx), cit_nbona, bonds_copy, &
  !                       use_atm_map)
  alloc_len = max(cit_nbona,proc_num_atms_min_bound) * mult_bonded
  if (.not. allocated(cit_a_bond)) allocate(cit_a_bond(alloc_len), stat =alloc_failed)
  if (.not. allocated(old_cit_a_bond)) allocate(old_cit_a_bond(alloc_len), stat =alloc_failed)
  if (alloc_failed .ne. 0) call setup_alloc_error

  if (cit_nbona .gt. 0) then
    cit_a_bond(1:cit_nbona) = bonds_copy(1:cit_nbona)
  end if

  if(.not. allocated(my_abonds_leads))  allocate(my_abonds_leads(alloc_len), stat = alloc_failed)
  if (alloc_failed .ne. 0) call setup_alloc_error
  if (cit_nbona .gt. 0) then
    do i = 1, cit_nbona
      my_abonds_leads(i) = cit_a_bond(i)%atm_i
    end do 
  end if 
  alloc_len = max(cit_nbona+cit_nbonh,proc_num_atms_min_bound) * mult_bonded
  if(.not. allocated(b_recv_buf)) allocate(b_recv_buf(alloc_len, neighbor_mpi_cnt), stat = alloc_failed)
  if(.not. allocated(b_send_buf)) allocate(b_send_buf(alloc_len, neighbor_mpi_cnt), stat = alloc_failed)
  if (alloc_failed .ne. 0) call setup_alloc_error
  ! Create an array (my_abonds_leads) of lead atoms that is useful after neighbor list 
  ! building to check the bond(with lead atom) which are to be transfred.
  !print *, "BondSetup,Rank,withH,w/o H:",mytaskid,cit_nbonh, cit_nbona
   
  return

end subroutine bonds_midpoint_setup     


!*******************************************************************************
!
! Subroutine:  find_mysd_bonds
!
! Description:  Finding sub-doamin bonds
!
!*******************************************************************************

subroutine find_mysd_bonds(bond_cnt, bonds, my_bond_cnt, my_bonds)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: bond_cnt
  type(bond_rec)        :: bonds(bond_cnt)
  integer               :: my_bond_cnt
  type(bond_rec)        :: my_bonds(*)

! Local variables:

  integer               :: atm_i, atm_j, bonds_idx
  integer               :: sub_atm
  integer               :: cnt_1, cnt_2

! Find all bonds for which this process owns either atom:
  my_bond_cnt = 0
  cnt_1 = 0
  cnt_2 = 0


  do bonds_idx = 1, bond_cnt

    atm_i = bonds(bonds_idx)%atm_i
    atm_j = bonds(bonds_idx)%atm_j
    
    ! Scan/search atoms from sub-domain : 
    if ((proc_atm_space(atm_i) /= 0) .AND.(proc_atm_space(atm_j) /= 0) ) then
        my_bond_cnt = my_bond_cnt + 1
        cnt_1 = cnt_1 + 1
        my_bonds(my_bond_cnt)%atm_i = proc_atm_space(atm_i)
        my_bonds(my_bond_cnt)%atm_j = proc_atm_space(atm_j)
        my_bonds(my_bond_cnt)%parm_idx    = bonds(bonds_idx)%parm_idx
    else if ((proc_atm_space(atm_i) /= 0) .AND.(proc_atm_space_ghosts(atm_j) /= 0) ) then
        my_bond_cnt = my_bond_cnt + 1
        cnt_2 = cnt_2 + 1
        my_bonds(my_bond_cnt)%atm_i = proc_atm_space(atm_i)
        my_bonds(my_bond_cnt)%atm_j = proc_atm_space_ghosts(atm_j)
        my_bonds(my_bond_cnt)%parm_idx    = bonds(bonds_idx)%parm_idx
    end if 
    
  end do

end subroutine find_mysd_bonds 

!*******************************************************************************
!
! Subroutine:  get_bond_energy_midpoint
!
! Description:
!              
! Routine to get bond energy and forces for the potential of cb*(b-b0)**2.
!
!*******************************************************************************

subroutine get_bond_energy_midpoint(bond_cnt, bond, x, frc, bond_energy,new_list,mult_vec)

  use nmr_calls_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use processor_mod
  use pbc_mod, only: pbc_box

  implicit none

! Formal arguments:

  integer               :: bond_cnt
  type(bond_rec)        :: bond(*)
  double precision      :: x(3, *)
  double precision      :: frc(3, *)
  double precision      :: bond_energy
  logical               :: new_list
  integer               :: mult_vec(3, *)

! Local variables:

  double precision      :: da, da_pull, da_press
  double precision      :: df
  double precision      :: dfw
  integer               :: i, j, ic, jn
  double precision      :: lcl_bond_energy
  double precision      :: xa, ya, za
  double precision      :: rij
  double precision      :: xij, yij, zij
  double precision      :: x1i,x2i,x3i
  double precision      :: x1j,x2j,x3j
  !integer               :: mult_vec(3,proc_num_atms+proc_ghost_num_atms*2)


  lcl_bond_energy = 0.0d0

! Grand loop for the bond stuff:

  do jn = 1, bond_cnt
    i = bond(jn)%atm_i
    j = bond(jn)%atm_j
    
! Calculation of the bond vector:
    x1i = x(1,i) 
    x2i = x(2,i) 
    x3i = x(3,i) 
    x1j = x(1,j) 
    x2j = x(2,j) 
    x3j = x(3,j) 

    x1j = x1j + mult_vec(1,jn) * pbc_box(1)
    x2j = x2j + mult_vec(2,jn) * pbc_box(2)
    x3j = x3j + mult_vec(3,jn) * pbc_box(3)
    xij = x1i - x1j
    yij = x2i - x2j
    zij = x3i - x3j
    rij = sqrt(xij * xij + yij * yij + zij * zij)
! Calculation of the energy and deriv:

    ic = bond(jn)%parm_idx
    da = rij - gbl_req(ic)
    if (rij < gbl_rpresseq(ic)) then
       da_press = gbl_rpresseq(ic) - rij
    else
       da_press = 0.0
    end if
    if (rij > gbl_rpulleq(ic)) then
       da_pull = gbl_rpulleq(ic) - rij
    else
       da_pull = 0.0
    end if
    
#ifdef MPI
! If ebdev is ever supported under mpi, you will need to treat it like the
! energy...
#else
    ! For rms deviation from ideal bonds:
    ebdev = ebdev + da * da + da_press * da_press + da_pull * da_pull 
#endif
    
    df = gbl_rk(ic) * da + gbl_rpressk(ic) * da_press + gbl_rpullk(ic) * da_pull
    dfw = (df + df) / rij

! Calculation of the force:

    xa = dfw * xij
    ya = dfw * yij
    za = dfw * zij

#ifdef MPI

      lcl_bond_energy = lcl_bond_energy + df * da
      frc(1, i) = frc(1, i) - xa
      frc(2, i) = frc(2, i) - ya
      frc(3, i) = frc(3, i) - za

      frc(1, j) = frc(1, j) + xa
      frc(2, j) = frc(2, j) + ya
      frc(3, j) = frc(3, j) + za
#else
      lcl_bond_energy = lcl_bond_energy + df * da
      frc(1, i) = frc(1, i) - xa
      frc(2, i) = frc(2, i) - ya
      frc(3, i) = frc(3, i) - za
      frc(1, j) = frc(1, j) + xa
      frc(2, j) = frc(2, j) + ya
      frc(3, j) = frc(3, j) + za
#endif

  end do

  bond_energy = lcl_bond_energy
  return

end subroutine get_bond_energy_midpoint

!*******************************************************************************
!
! Subroutine:  build_mult_vector_bond
!
! Description:
!              
! Routine to build multiplier vector to which is called only during the nb list 
! and used in every time step in get_bond_energy_midpoint
!
!*******************************************************************************

subroutine build_mult_vector_bond(bond_cnt,bond,x,mult_vec)

  use processor_mod, only: proc_dimx, proc_dimy, proc_dimz
  use parallel_dat_mod
  use prmtop_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: bond_cnt
  type(bond_rec)        :: bond(:)
  double precision      :: x(:, :)
  integer               :: mult_vec(:,:)

! Local variables:

  integer               :: i, j, ic, jn
  double precision      :: xij, yij, zij
  double precision      :: x1i,x2i,x3i
  double precision      :: x1j,x2j,x3j
  double precision      :: half_box_x, half_box_y, half_box_z

!we can move the following conditionals in the call site
!we probably dont need the conditional, since the wrapping can happen 
!when the simulation progress
!if(proc_dimx .eq. 1 .or. proc_dimy .eq. 1 .or. proc_dimz .eq. 1) then
  half_box_x = pbc_box(1)/2
  half_box_y = pbc_box(2)/2
  half_box_z = pbc_box(3)/2

  do jn = 1, bond_cnt

    i = bond(jn)%atm_i
    j = bond(jn)%atm_j
    
! Calculation of the bond vector:
    x1i = x(1,i) 
    x2i = x(2,i) 
    x3i = x(3,i) 
    x1j = x(1,j) 
    x2j = x(2,j) 
    x3j = x(3,j) 
!the following happens only during nb list is built, the following
!takes care of the fact that, atom_j can be wrapped, then here they are 
!unwrapped
      !atm_i will be unchanged
      !atm_j  have 3 variations
      if(abs(x1i - x1j) .le.half_box_x) then 
        mult_vec(1,jn) = 0
      else if(abs(x1i - (x1j + pbc_box(1)) ) .le.half_box_x) then 
        mult_vec(1,jn) = 1
      else if (abs(x1i - (x1j - pbc_box(1)) ) .le.half_box_x) then 
        mult_vec(1,jn) = -1
      end if
      if(abs(x2i - x2j) .le.half_box_y) then 
        mult_vec(2,jn) = 0
      else if(abs(x2i - (x2j + pbc_box(2)) ) .le.half_box_y) then 
        mult_vec(2,jn) = 1
      else if (abs(x2i - (x2j - pbc_box(2)) ) .le.half_box_y) then 
        mult_vec(2,jn) = -1
      end if
      if(abs(x3i - x3j) .le.half_box_z) then 
        mult_vec(3,jn) = 0
      else if(abs(x3i - (x3j + pbc_box(3)) ) .le.half_box_z) then 
        mult_vec(3,jn) = 1
      else if (abs(x3i - (x3j - pbc_box(3)) ) .le.half_box_z) then 
        mult_vec(3,jn) = -1
      end if
    end do
!end if
  return
end subroutine build_mult_vector_bond
#endif
end module bonds_midpoint_mod
