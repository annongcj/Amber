#include "copyright.i"
#include "hybrid_datatypes.i"

!*******************************************************************************
!
! Module:  dihedrals_midpoint_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module dihedrals_midpoint_mod
#ifdef MPI

  use gbl_datatypes_mod
  use processor_mod
  use bonds_midpoint_mod

  implicit none

! The following are derived from prmtop dihedral info:

  integer, save                         :: cit_nphih, cit_nphia
  type(dihed_rec), allocatable, save    :: cit_h_dihed(:)
  type(dihed_rec), allocatable, save    :: cit_a_dihed(:)
  type(dihed_rec), allocatable, save    :: old_cit_h_dihed(:)
  type(dihed_rec), allocatable, save    :: old_cit_a_dihed(:)
  integer, allocatable, save            :: my_hdiheds_leads(:)
  integer, allocatable, save            :: my_adiheds_leads(:)

!AMD
  double precision, save                :: fwgtd

  type(dihed_rec), allocatable, save    :: cit_dihed(:)
  type(dihed_trans),allocatable, save   :: d_recv_buf(:,:), d_send_buf(:,:)

contains

!*******************************************************************************
!
! Subroutine:  dihedrals_midpoint_setup
!
! Description:  <TBS>
!
!*******************************************************************************

!subroutine dihedrals_midpoint_setup(num_ints, num_reals, use_atm_map)
subroutine dihedrals_midpoint_setup()

  use parallel_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
  use processor_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  !integer, intent(in out)       :: num_ints, num_reals
  !integer                       :: use_atm_map(natom)

! Local variables:
  integer               :: i, atmid
  integer               :: alloc_failed
  type(dihed_rec)       :: diheds_copy(nphih + nphia)
  integer               :: atm_i, atm_j, atm_k, atm_l, diheds_idx
  integer               :: my_dihed_cnt
  integer               :: alloc_len

  !print *, "DIHEDsetup:nphih,nphia:",nphih, nphia
  mult_bonded = 10 ! Multiplier to allocate space for bonded interactions

! This routine can handle reallocation, and thus can be called multiple
! times.

  ! Dihedrals with Hydorgen atoms: 

  call find_mysd_dihedrals(nphih, gbl_dihed, cit_nphih, diheds_copy)
  alloc_len = max(cit_nphih,proc_num_atms_min_bound)
  if(.not. allocated(cit_h_dihed)) allocate(cit_h_dihed(alloc_len*mult_bonded), stat = alloc_failed)
  if(.not. allocated(old_cit_h_dihed)) allocate(old_cit_h_dihed(alloc_len*mult_bonded), stat = alloc_failed)
  if (alloc_failed .ne. 0) call setup_alloc_error
  if (cit_nphih .gt. 0) then
    cit_h_dihed(1:cit_nphih) = diheds_copy(1:cit_nphih)
  end if
  
  ! Create an array for lead/captain atoms for dihedrals with Hydrogen atoms (2nd atom,atm_j)
    if(.not. allocated(my_hdiheds_leads)) allocate(my_hdiheds_leads(alloc_len*mult_bonded), stat = alloc_failed)
  
  do i = 1, cit_nphih
    my_hdiheds_leads(i) = cit_h_dihed(i)%atm_j
  end do
 
  ! Dihedrals without Hydorgen atoms: 

  call find_mysd_dihedrals(nphia, gbl_dihed(diheda_idx), cit_nphia, diheds_copy)
  alloc_len = max(cit_nphia,proc_num_atms_min_bound)
  if(.not. allocated(cit_a_dihed)) allocate(cit_a_dihed(alloc_len*mult_bonded), stat = alloc_failed)
  if(.not. allocated(old_cit_a_dihed)) allocate(old_cit_a_dihed(alloc_len*mult_bonded), stat = alloc_failed)
  if (alloc_failed .ne. 0) call setup_alloc_error

  if (cit_nphia .gt. 0) then
    cit_a_dihed(1:cit_nphia) = diheds_copy(1:cit_nphia)
  end if

  
  ! Create an array for lead/captain atoms for dihedrals without Hydrogen atoms (2nd atom,atm_j)
  
  if(.not. allocated(my_adiheds_leads)) allocate(my_adiheds_leads(alloc_len*mult_bonded), stat = alloc_failed)
    if (alloc_failed .ne. 0) call setup_alloc_error
  
  do i = 1, cit_nphia
    my_adiheds_leads(i) = cit_a_dihed(i)%atm_j
  end do
  !if(.not. allocated(d_send_buf)) allocate(d_send_buf(cit_nphih+cit_nphia, neighbor_mpi_cnt), stat = alloc_failed)
  !if(.not. allocated(d_recv_buf)) allocate(d_recv_buf(cit_nphih+cit_nphia, neighbor_mpi_cnt), stat = alloc_failed)
  alloc_len = max(cit_nphih+cit_nphia,proc_num_atms_min_bound) * mult_bonded
  if(.not. allocated(d_send_buf)) allocate(d_send_buf(alloc_len, neighbor_mpi_cnt), stat = alloc_failed)
  if(.not. allocated(d_recv_buf)) allocate(d_recv_buf(alloc_len, neighbor_mpi_cnt), stat = alloc_failed)

  return

end subroutine dihedrals_midpoint_setup

!*******************************************************************************
!
! Subroutine:  find_mysd_dihedrals
!
! Description:  Finding sub-doamin dihedrals
!
!*******************************************************************************

!subroutine find_mysd_angles(angle_cnt, angles, my_angle_cnt, my_angles)
subroutine find_mysd_dihedrals(dihed_cnt, dihedrals, my_dihed_cnt, my_dihedrals)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer                :: dihed_cnt
  type(dihed_rec)        :: dihedrals(dihed_cnt)
  integer                :: my_dihed_cnt
  type(dihed_rec)        :: my_dihedrals(*)


! Local variables:

  integer               :: atm_i, atm_j, atm_k, atm_l, diheds_idx
  integer               :: sub_atm
  integer               :: cnt_1, cnt_2, cnt_3, cnt_4, cnt_5, cnt_6, cnt_7, &
                           cnt_8
  integer               :: sign_k, sign_l

! Find all dihedrals for which this process owns either atom:
! The dihedal is marked to the sub-domain where the 2nd atom(atm_j) is marked
  my_dihed_cnt = 0
  cnt_1 = 0
  cnt_2 = 0
  cnt_3 = 0
  cnt_4 = 0
  cnt_5 = 0
  cnt_6 = 0
  cnt_7 = 0
  cnt_8 = 0


  do diheds_idx = 1, dihed_cnt

    atm_i = dihedrals(diheds_idx)%atm_i
    atm_j = dihedrals(diheds_idx)%atm_j
    atm_k = iabs(dihedrals(diheds_idx)%atm_k)
    atm_l = iabs(dihedrals(diheds_idx)%atm_l)
    sign_k = sign(1,dihedrals(diheds_idx)%atm_k)
    sign_l = sign(1,dihedrals(diheds_idx)%atm_l)

    ! Scan/search atoms from sub-domain : 
    if ((proc_atm_space(atm_i) /= 0) .AND.(proc_atm_space(atm_j) /= 0) .AND. &
            (proc_atm_space(atm_k) /= 0) .AND.(proc_atm_space(atm_l) /= 0) )  then
        cnt_1 = cnt_1 + 1
        my_dihed_cnt = my_dihed_cnt + 1
        my_dihedrals(my_dihed_cnt)%atm_i = proc_atm_space(atm_i)
        my_dihedrals(my_dihed_cnt)%atm_j = proc_atm_space(atm_j)
        my_dihedrals(my_dihed_cnt)%atm_k = proc_atm_space(atm_k) * sign_k
        my_dihedrals(my_dihed_cnt)%atm_l = proc_atm_space(atm_l) * sign_l
        my_dihedrals(my_dihed_cnt)%parm_idx    = dihedrals(diheds_idx)%parm_idx

    else if ((proc_atm_space(atm_i) /= 0) .AND.(proc_atm_space(atm_j) /= 0) .AND. &
            (proc_atm_space(atm_k) /= 0) .AND.(proc_atm_space_ghosts(atm_l) /= 0) )  then
        cnt_2 = cnt_2 + 1
        my_dihed_cnt = my_dihed_cnt + 1
        my_dihedrals(my_dihed_cnt)%atm_i = proc_atm_space(atm_i)
        my_dihedrals(my_dihed_cnt)%atm_j = proc_atm_space(atm_j)
        my_dihedrals(my_dihed_cnt)%atm_k = proc_atm_space(atm_k) * sign_k
        my_dihedrals(my_dihed_cnt)%atm_l = proc_atm_space_ghosts(atm_l) * sign_l
        my_dihedrals(my_dihed_cnt)%parm_idx    = dihedrals(diheds_idx)%parm_idx


    else if ((proc_atm_space(atm_i) /= 0) .AND.(proc_atm_space(atm_j) /= 0) .AND. &
            (proc_atm_space_ghosts(atm_k) /= 0) .AND.(proc_atm_space(atm_l) /= 0) )  then
        cnt_3 = cnt_3 + 1
        my_dihed_cnt = my_dihed_cnt + 1
        my_dihedrals(my_dihed_cnt)%atm_i = proc_atm_space(atm_i)
        my_dihedrals(my_dihed_cnt)%atm_j = proc_atm_space(atm_j)
        my_dihedrals(my_dihed_cnt)%atm_k = proc_atm_space_ghosts(atm_k) * sign_k
        my_dihedrals(my_dihed_cnt)%atm_l = proc_atm_space(atm_l) * sign_l
        my_dihedrals(my_dihed_cnt)%parm_idx    = dihedrals(diheds_idx)%parm_idx

    else if ((proc_atm_space_ghosts(atm_i) /= 0) .AND.(proc_atm_space(atm_j) /= 0) .AND. &
            (proc_atm_space(atm_k) /= 0) .AND.(proc_atm_space(atm_l) /= 0) )  then
        cnt_4 = cnt_4 + 1
        my_dihed_cnt = my_dihed_cnt + 1
        my_dihedrals(my_dihed_cnt)%atm_i = proc_atm_space_ghosts(atm_i)
        my_dihedrals(my_dihed_cnt)%atm_j = proc_atm_space(atm_j)
        my_dihedrals(my_dihed_cnt)%atm_k = proc_atm_space(atm_k) * sign_k
        my_dihedrals(my_dihed_cnt)%atm_l = proc_atm_space(atm_l) * sign_l
        my_dihedrals(my_dihed_cnt)%parm_idx    = dihedrals(diheds_idx)%parm_idx


    else if ((proc_atm_space(atm_i) /= 0) .AND.(proc_atm_space(atm_j) /= 0) .AND. &
            (proc_atm_space_ghosts(atm_k) /= 0) .AND.(proc_atm_space_ghosts(atm_l) /= 0) )  then
        cnt_5 = cnt_5 + 1
        my_dihed_cnt = my_dihed_cnt + 1
        my_dihedrals(my_dihed_cnt)%atm_i = proc_atm_space(atm_i)
        my_dihedrals(my_dihed_cnt)%atm_j = proc_atm_space(atm_j)
        my_dihedrals(my_dihed_cnt)%atm_k = proc_atm_space_ghosts(atm_k) * sign_k
        my_dihedrals(my_dihed_cnt)%atm_l = proc_atm_space_ghosts(atm_l) * sign_l
        my_dihedrals(my_dihed_cnt)%parm_idx    = dihedrals(diheds_idx)%parm_idx


    else if ((proc_atm_space_ghosts(atm_i) /= 0) .AND.(proc_atm_space(atm_j) /= 0) .AND. &
            (proc_atm_space(atm_k) /= 0) .AND.(proc_atm_space_ghosts(atm_l) /= 0) )  then
        cnt_6 = cnt_6 + 1
        my_dihed_cnt = my_dihed_cnt + 1
        my_dihedrals(my_dihed_cnt)%atm_i = proc_atm_space_ghosts(atm_i)
        my_dihedrals(my_dihed_cnt)%atm_j = proc_atm_space(atm_j)
        my_dihedrals(my_dihed_cnt)%atm_k = proc_atm_space(atm_k) * sign_k
        my_dihedrals(my_dihed_cnt)%atm_l = proc_atm_space_ghosts(atm_l) * sign_l
        my_dihedrals(my_dihed_cnt)%parm_idx    = dihedrals(diheds_idx)%parm_idx

    else if ((proc_atm_space_ghosts(atm_i) /= 0) .AND.(proc_atm_space(atm_j) /= 0) .AND. &
            (proc_atm_space_ghosts(atm_k) /= 0) .AND.(proc_atm_space_ghosts(atm_l) /= 0) )  then
        cnt_7 = cnt_7 + 1
        my_dihed_cnt = my_dihed_cnt + 1
        my_dihedrals(my_dihed_cnt)%atm_i = proc_atm_space_ghosts(atm_i)
        my_dihedrals(my_dihed_cnt)%atm_j = proc_atm_space(atm_j)
        my_dihedrals(my_dihed_cnt)%atm_k = proc_atm_space_ghosts(atm_k) * sign_k
        my_dihedrals(my_dihed_cnt)%atm_l = proc_atm_space_ghosts(atm_l) * sign_l
        my_dihedrals(my_dihed_cnt)%parm_idx    = dihedrals(diheds_idx)%parm_idx


    else if ((proc_atm_space_ghosts(atm_i) /= 0) .AND.(proc_atm_space(atm_j) /= 0) .AND. &
            (proc_atm_space_ghosts(atm_k) /= 0) .AND.(proc_atm_space(atm_l) /= 0) )  then
        cnt_8 = cnt_8 + 1
        my_dihed_cnt = my_dihed_cnt + 1
        my_dihedrals(my_dihed_cnt)%atm_i = proc_atm_space_ghosts(atm_i)
        my_dihedrals(my_dihed_cnt)%atm_j = proc_atm_space(atm_j)
        my_dihedrals(my_dihed_cnt)%atm_k = proc_atm_space_ghosts(atm_k) * sign_k
        my_dihedrals(my_dihed_cnt)%atm_l = proc_atm_space(atm_l) * sign_l
        my_dihedrals(my_dihed_cnt)%parm_idx    = dihedrals(diheds_idx)%parm_idx

   end if  

 end do
 ! print *,"Rank,dihedcnt,1-3:",mytaskid,my_dihed_cnt,cnt_1,cnt_2,cnt_3
 ! print *,"Rank,dihedcnt,sum:",mytaskid,cnt_4,cnt_5,cnt_6,cnt_7
 
end subroutine find_mysd_dihedrals

!*******************************************************************************
!
! Subroutine:  get_dihed_midpoint_energy
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine get_dihed_midpoint_energy(dihed_cnt, dihed, x, frc,ep,mult_vec)

  use gbl_constants_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use pbc_mod, only: pbc_box

  implicit none

! Formal arguments:

  integer                       :: dihed_cnt
  type(dihed_rec)               :: dihed(*)
  double precision              :: x(3, *)
  double precision              :: frc(3, *)
  double precision              :: ep
  integer                       :: mult_vec(3,*)

! Local variables:

  double precision      :: ap
  double precision      :: cosnp
  double precision      :: cphi
  double precision      :: ct, ct0
  double precision      :: dc1, dc2, dc3, dc4, dc5, dc6
  double precision      :: df
  double precision      :: dr1, dr2, dr3, dr4, dr5, dr6
  double precision      :: drx, dry, drz
  double precision      :: dums
  double precision      :: dx, dy, dz
  double precision      :: epl
  double precision      :: epw
  double precision      :: f1, f2
  double precision      :: fxi, fyi, fzi
  double precision      :: fxj, fyj, fzj
  double precision      :: fxk, fyk, fzk
  double precision      :: fxl, fyl, fzl
  double precision      :: g
  double precision      :: gmul(10)
  double precision      :: gx, gy, gz
  integer               :: ic, ic0
  integer               :: inc
  integer               :: i, j, k, kt, l, lt
  integer               :: jn
  double precision      :: one
  double precision      :: s
  double precision      :: sinnp
  double precision      :: sphi
  double precision      :: tenm3
  double precision      :: tm06, tm24
  double precision      :: xa, ya, za
  double precision      :: xij, yij, zij
  double precision      :: xkj, ykj, zkj
  double precision      :: xkl, ykl, zkl
  double precision      :: z1, z2
  double precision      :: z11, z12, z22
  double precision      :: zero
  integer               :: kn

  data gmul /0.d0, 2.d0, 0.d0, 4.d0, 0.d0, 6.d0, 0.d0, 8.d0, 0.d0, 10.d0/

  data tm24, tm06, tenm3/1.d-18, 1.d-06, 1.d-03/

  data zero, one /0.d0, 1.d0/

! Arrays gbl_gamc = gbl_pk * cos(phase) and gbl_gams = gbl_pk * sin(phase)

  epl = zero

! Grand loop for the dihedral stuff:

  do jn = 1, dihed_cnt

    i = dihed(jn)%atm_i
    j = dihed(jn)%atm_j
    kt = dihed(jn)%atm_k
    lt = dihed(jn)%atm_l
    k = iabs(kt)
    l = iabs(lt)
    kn = (jn - 1) *3 + 1
! Calculation of ij, kj, kl vectors:

    xij = x(1, i) - (x(1, j) + mult_vec(1,kn) * pbc_box(1))
    yij = x(2, i) - (x(2, j) + mult_vec(2,kn) * pbc_box(2))
    zij = x(3, i) - (x(3, j) + mult_vec(3,kn) * pbc_box(3))
    xkj = x(1, k) + mult_vec(1,kn+1) * pbc_box(1) - (x(1, j) + mult_vec(1,kn) * pbc_box(1))
    ykj = x(2, k) + mult_vec(2,kn+1) * pbc_box(2) - (x(2, j) + mult_vec(2,kn) * pbc_box(2))
    zkj = x(3, k) + mult_vec(3,kn+1) * pbc_box(3) - (x(3, j) + mult_vec(3,kn) * pbc_box(3))
    xkl = x(1, k) + mult_vec(1,kn+1) * pbc_box(1) - (x(1, l) + mult_vec(1,kn+2) *pbc_box(1))
    ykl = x(2, k) + mult_vec(2,kn+1) * pbc_box(2) - (x(2, l) + mult_vec(2,kn+2) *pbc_box(2))
    zkl = x(3, k) + mult_vec(3,kn+1) * pbc_box(3) - (x(3, l) + mult_vec(3,kn+2) *pbc_box(3))  

! Get the normal vector:

    dx = yij * zkj - zij * ykj
    dy = zij * xkj - xij * zkj
    dz = xij * ykj - yij * xkj
    gx = zkj * ykl - ykj * zkl
    gy = xkj * zkl - zkj * xkl
    gz = ykj * xkl - xkj * ykl

    fxi = sqrt(dx * dx + dy * dy + dz * dz + tm24)
    fyi = sqrt(gx * gx + gy * gy + gz * gz + tm24)
    ct = dx * gx + dy * gy + dz * gz

! Branch if linear dihedral:

    if (tenm3 .le. fxi) then
      z1 = one / fxi
    else
      z1 = zero
    endif

    if (tenm3 .le. fyi) then
      z2 = one / fyi
    else
      z2 = zero
    endif

    z12 = z1 * z2

    if (z12 .ne. zero) then
      fzi = one
    else
      fzi = zero
    end if

    s = xkj * (dz * gy - dy * gz) + ykj * (dx * gz - dz * gx) + &
        zkj * (dy * gx - dx * gy)

    ap = PI - sign(acos(max(-one, min(one, ct * z12))), s)

    cphi = cos(ap)
    sphi = sin(ap)

! Calculate the energy and the derivatives with respect to cosphi:

    ic = dihed(jn)%parm_idx
    inc = gbl_ipn(ic)
    ct0 = gbl_pn(ic) * ap
    cosnp = cos(ct0)
    sinnp = sin(ct0)
    epw = (gbl_pk(ic) + cosnp * gbl_gamc(ic) + sinnp * gbl_gams(ic)) * fzi
    
    dums = sphi + sign(tm24, sphi)

    if (tm06 .gt. abs(dums)) then
      df = fzi * gbl_gamc(ic) * (gbl_pn(ic) - gmul(inc) + gmul(inc) * cphi)
    else
      df = fzi * gbl_pn(ic) * (gbl_gamc(ic) * sinnp - gbl_gams(ic) * cosnp)/dums
    end if

! Now do torsional first derivatives:

! Now, set up array dc = 1st der. of cosphi w/respect to cartesian differences:

    z11 = z1 * z1
    z12 = z1 * z2
    z22 = z2 * z2
    dc1 = -gx * z12 - cphi * dx * z11
    dc2 = -gy * z12 - cphi * dy * z11
    dc3 = -gz * z12 - cphi * dz * z11
    dc4 =  dx * z12 + cphi * gx * z22
    dc5 =  dy * z12 + cphi * gy * z22
    dc6 =  dz * z12 + cphi * gz * z22

! Update the first derivative array:

    dr1 = df * ( dc3 * ykj - dc2 * zkj)
    dr2 = df * ( dc1 * zkj - dc3 * xkj)
    dr3 = df * ( dc2 * xkj - dc1 * ykj)
    dr4 = df * ( dc6 * ykj - dc5 * zkj)
    dr5 = df * ( dc4 * zkj - dc6 * xkj)
    dr6 = df * ( dc5 * xkj - dc4 * ykj)
    drx = df * (-dc2 * zij + dc3 * yij + dc5 * zkl - dc6 * ykl)
    dry = df * ( dc1 * zij - dc3 * xij - dc4 * zkl + dc6 * xkl)
    drz = df * (-dc1 * yij + dc2 * xij + dc4 * ykl - dc5 * xkl)
    fxi = -dr1
    fyi = -dr2
    fzi = -dr3
    fxj = -drx + dr1
    fyj = -dry + dr2
    fzj = -drz + dr3
    fxk = +drx + dr4
    fyk = +dry + dr5
    fzk = +drz + dr6
    fxl = -dr4
    fyl = -dr5
    fzl = -dr6

#ifdef MPI

      epl = epl  + epw
      frc(1, i) = frc(1, i) + fxi
      frc(2, i) = frc(2, i) + fyi
      frc(3, i) = frc(3, i) + fzi

      frc(1, j) = frc(1, j) + fxj 
      frc(2, j) = frc(2, j) + fyj
      frc(3, j) = frc(3, j) + fzj

      frc(1, k) = frc(1, k) + fxk
      frc(2, k) = frc(2, k) + fyk
      frc(3, k) = frc(3, k) + fzk

      frc(1, l) = frc(1, l) + fxl
      frc(2, l) = frc(2, l) + fyl
      frc(3, l) = frc(3, l) + fzl
#else
    epl = epl + epw

    frc(1, i) = frc(1, i) + fxi
    frc(2, i) = frc(2, i) + fyi
    frc(3, i) = frc(3, i) + fzi
    frc(1, j) = frc(1, j) + fxj 
    frc(2, j) = frc(2, j) + fyj
    frc(3, j) = frc(3, j) + fzj
    frc(1, k) = frc(1, k) + fxk
    frc(2, k) = frc(2, k) + fyk
    frc(3, k) = frc(3, k) + fzk
    frc(1, l) = frc(1, l) + fxl
    frc(2, l) = frc(2, l) + fyl
    frc(3, l) = frc(3, l) + fzl
#endif
    
  end do

  ep = epl

  return

end subroutine get_dihed_midpoint_energy
!*******************************************************************************
!
! Subroutine:  build_mult_vector_dihed
!
! Description:
!              
! Routine to build multiplier vector to which is called only during the nb list 
! and used in every time step in get_dihd_energy_midpoint
!
!*******************************************************************************

subroutine build_mult_vector_dihed(dihed_cnt,dihed,x,mult_vec)

  use processor_mod, only: proc_dimx, proc_dimy, proc_dimz
  use parallel_dat_mod
  use prmtop_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: dihed_cnt
  type(dihed_rec)       :: dihed(*)
  double precision      :: x(3, *)
  integer               :: mult_vec(3,*)

! Local variables:

  integer               :: i, j, k, l, ic, jn, kn
  double precision      :: xij, yij, zij
  double precision      :: x1i,x2i,x3i
  double precision      :: x1j,x2j,x3j
  double precision      :: x1k,x2k,x3k
  double precision      :: x1l,x2l,x3l
  double precision      :: half_box_x, half_box_y, half_box_z
  

!we can move the following conditionals in the call site
!we probably dont need the conditional, since the wrapping can happen 
!when the simulation progress
!if(proc_dimx .eq. 1 .or. proc_dimy .eq. 1 .or. proc_dimz .eq. 1) then
  half_box_x = pbc_box(1)/2.0
  half_box_y = pbc_box(2)/2.0
  half_box_z = pbc_box(3)/2.0

  do jn = 1, dihed_cnt

    i = dihed(jn)%atm_i
    j = dihed(jn)%atm_j
    k = iabs (dihed(jn)%atm_k)
    l = iabs (dihed(jn)%atm_l)
    
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
    x1l = x(1,l) 
    x2l = x(2,l) 
    x3l = x(3,l) 
    kn = (jn -1) * 3 + 1 !index for mult_vec for j,k,l in respect of i
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

      if(abs(x1i - x1l) .lt. half_box_x) then 
        mult_vec(1,kn+2) = 0
      else if(abs(x1i - (x1l + pbc_box(1)) ) .lt. half_box_x) then 
        mult_vec(1,kn+2) = 1
      else if (abs(x1i - (x1l - pbc_box(1)) ) .lt. half_box_x) then 
        mult_vec(1,kn+2) = -1
      end if
      if(abs(x2i - x2l) .lt. half_box_y) then 
        mult_vec(2,kn+2) = 0
      else if(abs(x2i - (x2l + pbc_box(2)) ) .le. half_box_y) then 
        mult_vec(2,kn+2) = 1
      else if (abs(x2i - (x2l - pbc_box(2)) ) .le. half_box_y) then 
        mult_vec(2,kn+2) = -1
      end if
      if(abs(x3i - x3l) .lt. half_box_z) then 
        mult_vec(3,kn+2) = 0
      else if(abs(x3i - (x3l + pbc_box(3)) ) .le. half_box_z) then 
        mult_vec(3,kn+2) = 1
      else if (abs(x3i - (x3l - pbc_box(3)) ) .le. half_box_z) then 
        mult_vec(3,kn+2) = -1
      end if
    end do
  return
end subroutine build_mult_vector_dihed
#endif

end module dihedrals_midpoint_mod
