#include "copyright.i"
#include "hybrid_datatypes.i"

!*******************************************************************************
!
! Module:  angles_midpoint_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module angles_midpoint_mod

#ifdef MPI

  use gbl_datatypes_mod
  use processor_mod
  use bonds_midpoint_mod

  implicit none

! The following are derived from prmtop angle info:

  integer, save                         :: cit_ntheth, cit_ntheta
  type(angle_rec), allocatable, save    :: cit_h_angle(:)
  type(angle_rec), allocatable, save    :: cit_a_angle(:)
  type(angle_rec), allocatable, save    :: old_cit_h_angle(:)
  type(angle_rec), allocatable, save    :: old_cit_a_angle(:)
  integer, allocatable, save            :: my_hangles_leads(:)
  integer, allocatable, save            :: my_aangles_leads(:)
  type(angle_trans),allocatable, save   :: a_recv_buf(:,:), a_send_buf(:,:)

contains

!*******************************************************************************
!
! Subroutine:  angles_midpoint_setup
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine angles_midpoint_setup()

  use parallel_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
  use processor_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  ! integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: i, atmid 
  integer               :: alloc_failed
  integer               :: mult_bond = 1.5
  type(angle_rec)       :: angles_copy(ntheth + ntheta)
  integer               :: atm_i, atm_j, atm_k, angles_idx
  integer               :: my_angle_cnt
  integer               :: alloc_len

  mult_bonded = 10 ! Multiplier to allocate space for bonded interactions
  ! Scan Angles for the sub-domain 
  ! This routine can handle reallocation, and thus can be called multiple
  ! times.
  ! Scanning angle list :
  
  !call find_my_bonds(nbonh, gbl_bond, cit_nbonh, bonds_copy, use_atm_map)
  call find_mysd_angles(ntheth, gbl_angle, cit_ntheth, angles_copy)
  alloc_len = max(cit_ntheth,proc_num_atms_min_bound) * mult_bonded
  if (.not. allocated(cit_h_angle)) allocate(cit_h_angle(alloc_len), stat =alloc_failed)
  if(.not. allocated(old_cit_h_angle)) allocate(old_cit_h_angle(alloc_len),stat = alloc_failed)
  if (alloc_failed .ne. 0) call setup_alloc_error
  if (cit_ntheth .gt. 0) then
    cit_h_angle(1:cit_ntheth) = angles_copy(1:cit_ntheth)
  end if
  
  ! Create an array for lead/captain atoms for angles (middle atom,atm_j)
  
  if(.not. allocated(my_hangles_leads)) allocate(my_hangles_leads(alloc_len), stat = alloc_failed)
  if (alloc_failed .ne. 0) call setup_alloc_error
  
  do i = 1, cit_ntheth
    my_hangles_leads(i) = cit_h_angle(i)%atm_j
  end do
 
  ! Angles without Hydorgen atoms: 

  call find_mysd_angles(ntheta, gbl_angle(anglea_idx), cit_ntheta, angles_copy)
 
  alloc_len = max(cit_ntheta,proc_num_atms_min_bound) * mult_bonded
  if (.not. allocated(cit_a_angle)) allocate(cit_a_angle(alloc_len), stat =alloc_failed)
  if(.not. allocated(old_cit_a_angle)) allocate(old_cit_a_angle(alloc_len),stat = alloc_failed)
  if (alloc_failed .ne. 0) call setup_alloc_error
  if (cit_ntheta .gt. 0) then
    cit_a_angle(1:cit_ntheta) = angles_copy(1:cit_ntheta)
  end if
  
  
  ! Create an array for lead/captain atoms for angles (middle atom,atm_j)
  
   if(.not. allocated(my_aangles_leads)) allocate(my_aangles_leads(alloc_len), stat = alloc_failed)
  if (alloc_failed .ne. 0) call setup_alloc_error
  
  do i = 1, cit_ntheta
    my_aangles_leads(i) = cit_a_angle(i)%atm_j
  end do
  alloc_len = max(cit_ntheta+cit_ntheth,proc_num_atms_min_bound) * mult_bonded
  if(.not. allocated(a_send_buf)) allocate(a_send_buf(alloc_len, neighbor_mpi_cnt), stat = alloc_failed)
  !if(.not. allocated(a_send_buf)) allocate(a_send_buf(cit_ntheth+cit_ntheta, neighbor_mpi_cnt), stat = alloc_failed)
  if(.not. allocated(a_recv_buf)) allocate(a_recv_buf(alloc_len, neighbor_mpi_cnt), stat = alloc_failed)
  !if(.not. allocated(a_recv_buf)) allocate(a_recv_buf(cit_ntheth+cit_ntheta, neighbor_mpi_cnt), stat = alloc_failed)
  return

end subroutine angles_midpoint_setup     

!*******************************************************************************
!
! Subroutine:  find_mysd_angles
!
! Description:  Finding sub-doamin angles
!
!*******************************************************************************

subroutine find_mysd_angles(angle_cnt, angles, my_angle_cnt, my_angles)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer                :: angle_cnt
  type(angle_rec)        :: angles(angle_cnt)
  integer                :: my_angle_cnt
  type(angle_rec)        :: my_angles(*)


! Local variables:

  integer               :: atm_i, atm_j, atm_k, angles_idx
  integer               :: sub_atm
  integer               :: cnt_1, cnt_2, cnt_3, cnt_4


! Find all angles for which this process owns either atom:
! The angle is marked to the sub-domain where the middle atom(atm_j) is marked
  my_angle_cnt = 0
  cnt_1 = 0
  cnt_2 = 0
  cnt_3 = 0
  cnt_4 = 0


  do angles_idx = 1, angle_cnt

    atm_i = angles(angles_idx)%atm_i
    atm_j = angles(angles_idx)%atm_j
    atm_k = angles(angles_idx)%atm_k
    
    ! Scan/search atoms from sub-domain : 
    if ((proc_atm_space(atm_i) /= 0) .AND.(proc_atm_space(atm_j) /= 0) .AND. &
            (proc_atm_space(atm_k) /= 0))              then
        my_angle_cnt = my_angle_cnt + 1
        cnt_1 = cnt_1 + 1
        my_angles(my_angle_cnt)%atm_i = proc_atm_space(atm_i)
        my_angles(my_angle_cnt)%atm_j = proc_atm_space(atm_j)
        my_angles(my_angle_cnt)%atm_k = proc_atm_space(atm_k)
        my_angles(my_angle_cnt)%parm_idx    = angles(angles_idx)%parm_idx
    else if ((proc_atm_space(atm_j) /= 0) .AND.(proc_atm_space_ghosts(atm_i) /= 0) .AND. &
            (proc_atm_space_ghosts(atm_k) /= 0))              then
        my_angle_cnt = my_angle_cnt + 1
        cnt_2 = cnt_2 + 1
        my_angles(my_angle_cnt)%atm_i = proc_atm_space_ghosts(atm_i)
        my_angles(my_angle_cnt)%atm_j = proc_atm_space(atm_j)
        my_angles(my_angle_cnt)%atm_k = proc_atm_space_ghosts(atm_k)
        my_angles(my_angle_cnt)%parm_idx    = angles(angles_idx)%parm_idx
    else if ((proc_atm_space(atm_j) /= 0) .AND.(proc_atm_space_ghosts(atm_i) /= 0) .AND. &
            (proc_atm_space(atm_k) /= 0))              then
        my_angle_cnt = my_angle_cnt + 1
        cnt_3 = cnt_3 + 1
        my_angles(my_angle_cnt)%atm_i = proc_atm_space_ghosts(atm_i)
        my_angles(my_angle_cnt)%atm_j = proc_atm_space(atm_j)
        my_angles(my_angle_cnt)%atm_k = proc_atm_space(atm_k)
        my_angles(my_angle_cnt)%parm_idx    = angles(angles_idx)%parm_idx
    else if ((proc_atm_space(atm_j) /= 0) .AND.(proc_atm_space(atm_i) /= 0) .AND. &
            (proc_atm_space_ghosts(atm_k) /= 0))              then
        my_angle_cnt = my_angle_cnt + 1
        cnt_4 = cnt_4 + 1
        my_angles(my_angle_cnt)%atm_i = proc_atm_space(atm_i)
        my_angles(my_angle_cnt)%atm_j = proc_atm_space(atm_j)
        my_angles(my_angle_cnt)%atm_k = proc_atm_space_ghosts(atm_k)
        my_angles(my_angle_cnt)%parm_idx    = angles(angles_idx)%parm_idx
    end if 
    
  end do
   ! print *, "##Angles,Rank,cnt_1. cnt_2,3,4:",mytaskid,cnt_1,cnt_2,cnt_3,cnt_4

end subroutine find_mysd_angles

!*******************************************************************************
!
! Subroutine:  get_angles_energy_midpoint
!
!
! Description:  Routine to get the bond energies and forces for potentials of
!               the type ct*(t-t0)**2.
!
!*******************************************************************************

subroutine get_angle_energy_midpoint(angle_cnt, angle, x, frc,angle_energy,mult_vec)

  use nmr_calls_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use processor_mod
  use pbc_mod, only: pbc_box

  implicit none

! Formal arguments:

  integer               :: angle_cnt
  type(angle_rec)       :: angle(*)
  double precision      :: x(3, *)
  double precision      :: frc(3, *)
  double precision      :: angle_energy
  integer               :: mult_vec(3,*)


! Local variables:

  double precision, parameter   :: pt999 = 0.9990d0

  double precision      :: ant
  double precision      :: cst
  double precision      :: dfw
  double precision      :: eaw
  double precision      :: rij, rik, rkj
  double precision      :: xij, yij, zij
  double precision      :: xkj, ykj, zkj
  double precision      :: cii, cik, ckk
  double precision      :: da
  double precision      :: df
  double precision      :: dt1, dt2, dt3, dt4, dt5, dt6
  double precision      :: sth
  integer               :: i, j, k, ic, jn, kn

  angle_energy = 0.0d0


! Grand loop for the angle stuff:

  do jn = 1, angle_cnt

    i = angle(jn)%atm_i
    j = angle(jn)%atm_j
    k = angle(jn)%atm_k
    ic = angle(jn)%parm_idx
    kn = (jn - 1) *2 + 1

! Calculation of the angle:

    xij = x(1, i) - (x(1, j) + mult_vec(1,kn) * pbc_box(1))
    xkj = x(1, k) + mult_vec(1,kn+1) * pbc_box(1) - (x(1, j) + mult_vec(1,kn) * pbc_box(1))

    yij = x(2, i) - (x(2, j) + mult_vec(2,kn) * pbc_box(2))
    ykj = x(2, k) + mult_vec(2,kn+1) * pbc_box(2) - (x(2, j) + mult_vec(2,kn) * pbc_box(2))

    zij = x(3, i) - (x(3, j) + mult_vec(3,kn) * pbc_box(3))
    zkj = x(3, k) + mult_vec(3,kn+1) * pbc_box(3) - (x(3, j) + mult_vec(3,kn) * pbc_box(3))

    rij = xij * xij + yij * yij + zij * zij
    rkj = xkj * xkj + ykj * ykj + zkj * zkj

    rik = sqrt(rij * rkj)

    cst = min(pt999, max(-pt999, (xij * xkj + yij * ykj + zij * zkj) / rik))

    ant = acos(cst)

! Calculation of the energy and deriv:

    da = ant - gbl_teq(ic)

#ifdef MPI
! If eadev is ever supported under mpi, you will need to split the contribution.
#else
    eadev = eadev + da * da         ! For rms deviation from ideal angles.
#endif

    df = gbl_tk(ic) * da
    eaw = df * da
    dfw = -(df + df) / sin(ant)

! Calculation of the force: 

    cik = dfw / rik
    sth = dfw * cst
    cii = sth / rij
    ckk = sth / rkj

    dt1 = cik * xkj - cii * xij
    dt2 = cik * ykj - cii * yij
    dt3 = cik * zkj - cii * zij
    dt4 = cik * xij - ckk * xkj
    dt5 = cik * yij - ckk * ykj
    dt6 = cik * zij - ckk * zkj
    

      angle_energy = angle_energy + eaw

      frc(1, i) = frc(1, i) - dt1
      frc(2, i) = frc(2, i) - dt2
      frc(3, i) = frc(3, i) - dt3

      frc(1, j) = frc(1, j) + dt1 + dt4
      frc(2, j) = frc(2, j) + dt2 + dt5
      frc(3, j) = frc(3, j) + dt3 + dt6


      frc(1, k) = frc(1, k) - dt4
      frc(2, k) = frc(2, k) - dt5
      frc(3, k) = frc(3, k) - dt6

        
  end do

  return

end subroutine get_angle_energy_midpoint

!*******************************************************************************
!
! Subroutine:  build_mult_vector_angle
!
! Description:
!              
! Routine to build multiplier vector to which is called only during the nb list 
! and used in every time step in get_bond_energy_midpoint
!
!*******************************************************************************

subroutine build_mult_vector_angle(angle_cnt,angle,x,mult_vec)

  use processor_mod, only: proc_num_atms, proc_num_atms_min_bound
  use parallel_dat_mod
  use prmtop_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: angle_cnt
  type(angle_rec)       :: angle(:)
  double precision      :: x(:,:)
  integer               :: mult_vec(:,:)

! Local variables:

  integer               :: i, j, k, ic, jn, kn
  double precision      :: xij, yij, zij
  double precision      :: x1i,x2i,x3i
  double precision      :: x1j,x2j,x3j
  double precision      :: x1k,x2k,x3k
  double precision      :: half_box_x, half_box_y, half_box_z
  

!we can move the following conditionals in the call site
!we probably dont need the conditional, since the wrapping can happen 
!when the simulation progress
!if(proc_dimx .eq. 1 .or. proc_dimy .eq. 1 .or. proc_dimz .eq. 1) then
  half_box_x = pbc_box(1)/2.0
  half_box_y = pbc_box(2)/2.0
  half_box_z = pbc_box(3)/2.0
  do jn = 1, angle_cnt

    i = angle(jn)%atm_i
    j = angle(jn)%atm_j
    k = angle(jn)%atm_k
    
! Calculation of the bond vector:
    x1i = x(1,i) 
    x2i = x(2,i) 
    x3i = x(3,i) 
    x1j = x(1,j) 
    x2j = x(2,j) 
    x3j = x(3,j) 
    x1k = x(1,k) 
    x2k = x(2,k) 
    x3k = x(3,k) 
  
    kn = (jn -1)*2 + 1  !index for mult_vec for j,k in respect of i

!the following happens only during nb list is built, the following
!takes care of the fact that, atom_j can be wrapped, then here they are 
!unwrapped
      !atm_i will be unchanged
      !atm_j  have 3 variations
      if(abs(x1i - x1j) .lt. half_box_x) then 
        mult_vec(1,kn) = 0
      else if(abs(x1i - (x1j + pbc_box(1)) ) .lt. half_box_x) then 
        mult_vec(1,kn) = 1
      else if (abs(x1i - (x1j - pbc_box(1)) ) .lt. half_box_x) then 
        mult_vec(1,kn) = -1
      end if
      if(abs(x2i - x2j) .lt. half_box_y) then 
        mult_vec(2,kn) = 0
      else if(abs(x2i - (x2j + pbc_box(2)) ) .le. half_box_y) then 
        mult_vec(2,kn) = 1
      else if (abs(x2i - (x2j - pbc_box(2)) ) .le. half_box_y) then 
        mult_vec(2,kn) = -1
      end if
      if(abs(x3i - x3j) .lt. half_box_z) then 
        mult_vec(3,kn) = 0
      else if(abs(x3i - (x3j + pbc_box(3)) ) .le. half_box_z) then 
        mult_vec(3,kn) = 1
      else if (abs(x3i - (x3j - pbc_box(3)) ) .le. half_box_z) then 
        mult_vec(3,kn) = -1
      end if

      if(abs(x1i - x1k) .lt. half_box_x) then 
        mult_vec(1,kn+1) = 0
      else if(abs(x1i - (x1k + pbc_box(1)) ) .lt. half_box_x) then 
        mult_vec(1,kn+1) = 1
      else if (abs(x1i - (x1k - pbc_box(1)) ) .lt. half_box_x) then 
        mult_vec(1,kn+1) = -1
      end if
      if(abs(x2i - x2k) .lt. half_box_y) then 
        mult_vec(2,kn+1) = 0
      else if(abs(x2i - (x2k + pbc_box(2)) ) .le. half_box_y) then 
        mult_vec(2,kn+1) = 1
      else if (abs(x2i - (x2k - pbc_box(2)) ) .le. half_box_y) then 
        mult_vec(2,kn+1) = -1
      end if
      if(abs(x3i - x3k) .lt. half_box_z) then 
        mult_vec(3,kn+1) = 0
      else if(abs(x3i - (x3k + pbc_box(3)) ) .le. half_box_z) then 
        mult_vec(3,kn+1) = 1
      else if (abs(x3i - (x3k - pbc_box(3)) ) .le. half_box_z) then 
        mult_vec(3,kn+1) = -1
      end if

    end do
  return
end subroutine build_mult_vector_angle
#endif /*MPI*/

end module angles_midpoint_mod

