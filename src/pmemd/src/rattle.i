  use constraints_mod
  use dynamics_dat_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use pbc_mod
  use prmtop_dat_mod
  use pmemd_lib_mod
#ifdef _OPENMP_
  use omp_lib
#endif
  use processor_mod

  implicit none

! Formal arguments:

  double precision       :: x(3, *)
  double precision       :: v(3, *)

! Local variables:

  double precision       :: acor
  double precision       :: box_half_min
  double precision       :: box_half_min_sq
  double precision       :: diff
  integer                :: i, j, m
  integer                :: iter_cnt
  integer                :: my_bond_idx
  integer                :: ns
  integer                :: ibelly_lcl
  logical                :: done

  double precision       :: rpij2
  double precision       :: rvdot
  double precision       :: rrpr
  double precision       :: tol_lcl
  double precision       :: toler
  double precision       :: winvi, winvj
  double precision       :: xh
  double precision       :: xij(3)
  double precision       :: vij(3)

  double precision, save :: zero = 0.0d0
#ifdef _OPENMP_
  integer                :: my_thread_id 
#endif

  integer,save,allocatable :: skips(:)
  if(allocated(skips)) deallocate(skips)
if(usemidpoint) then
#ifdef MPI
  allocate(skips(proc_num_atms_min_bound + proc_ghost_num_atms))
#endif
else
  allocate(skips(atm_cnt))
endif

if(usemidpoint) then
  if (my_nonfastwat_bond_cnt_midpoint .eq. 0) return
else
  if (my_nonfastwat_bond_cnt .eq. 0) return
endif

  ! BUGBUG - Initialization of the skips() array does not scale...

  ! Calc size of max vector that is guaranteed to fit in box w/out imaging:

  box_half(:) = pbc_box(:) * 0.5d0
  box_half_min = min(box_half(1), box_half(2), box_half(3))

  box_half_min_sq = box_half_min * box_half_min

  ibelly_lcl = ibelly
  tol_lcl = tol / dt

#ifdef GBorn
  skips(:) = 1
#else

if(usemidpoint) then
  do my_bond_idx = 1, my_nonfastwat_bond_cnt_midpoint
      i = my_nonfastwat_bond_dat_midpoint(my_bond_idx)%atm_i
      j = my_nonfastwat_bond_dat_midpoint(my_bond_idx)%atm_j
     skips(i) = 1
     skips(j) = 1
  end do
else
  do my_bond_idx = 1, my_nonfastwat_bond_cnt
      i = my_nonfastwat_bond_dat(my_bond_idx)%atm_i
      j = my_nonfastwat_bond_dat(my_bond_idx)%atm_j
     skips(i) = 1
     skips(j) = 1
  end do
endif

#endif

  do iter_cnt = 1, 3000

    done = .true.       ! until proven otherwise.

! Loop over all the bonds that are not 3-point waters:
#ifdef MPI
#ifdef GBorn
if(master) then
    !$omp parallel default (none) &
    !$omp& private(my_bond_idx,i,j,rpij2,m,xpij,ns,toler,diff,xij,rrpr, &
    !$omp& select_bond_cnt, &
    !$omp& winvi,winvj,acor,xh,my_thread_id) &
    !$omp& shared(my_nonfastwat_bond_cnt,skips,zero,xp,ntb,box_half_min_sq, &
    !$omp& atm_igroup, my_nonfastwat_bond_dat, ibelly_lcl,gbl_atm_tid_map,&
    !$omp& box_half,pbc_box,tol_lcl,x,master,imin,iter_cnt,atm_mass_inv,done, &
    !$omp& usemidpoint,my_nonfastwat_bond_cnt_midpoint, my_nonfastwat_bond_dat_midpoint, &
    !$omp& mult_vec_shake_bonds, proc_atm_mass)
    my_thread_id = omp_get_thread_num()
#endif
#endif

if(usemidpoint) then
    select_bond_cnt = my_nonfastwat_bond_cnt_midpoint
else
    select_bond_cnt = my_nonfastwat_bond_cnt
endif
    do my_bond_idx = 1, select_bond_cnt

if(usemidpoint) then
      i = my_nonfastwat_bond_dat_midpoint(my_bond_idx)%atm_i
      j = my_nonfastwat_bond_dat_midpoint(my_bond_idx)%atm_j
else
      i = my_nonfastwat_bond_dat(my_bond_idx)%atm_i
      j = my_nonfastwat_bond_dat(my_bond_idx)%atm_j
endif

#ifdef MPI
#ifdef GBorn
      if(gbl_atm_tid_map(i) .ne.  my_thread_id) cycle
#endif
#endif


if(.not. usemidpoint) then
      if (skips(i) .lt. iter_cnt .and. skips(j) .lt. iter_cnt) cycle
endif

! Calc nominal distance squared:

      rvdot = zero

if(usemidpoint) then
      do  m = 1, 3
        xij(m) = x(m, i) - x(m, j) - mult_vec_shake_bonds(m,my_bond_idx) * pbc_box(m) 
        vij(m) = v(m, i) - v(m, j)
        rvdot = rvdot + xij(m) * vij(m)
      end do
else
      do  m = 1, 3
        xij(m) = x(m, i) - x(m, j)
        vij(m) = v(m, i) - v(m, j)
        rvdot = rvdot + xij(m) * vij(m)
      end do
endif

! BUGBUG - The following boundary check is no longer used in sander 8+, and
!          is probably unnecessary since we are dealing with atom coords (not
!          image crds) here and this stuff is bonded...  We leave it in for
!          now since it is probably not doing any harm in terms of results and
!          minimal harm in terms of performance...

! If boundary condition is not present skip it:

      ns = 0 

      if (ntb .gt. 0) then

        if (rpij2 .ge. box_half_min_sq) then

! Apply the boundary & recalc the distance squared:

          ns = 1
          call wrap_crds(xij)
          rvdot = zero

          do m = 1, 3
            rpij2 = rpij2 + xij(m) * vij(m)
          end do

        end if

      end if

! Apply the correction:

if(usemidpoint) then
      toler = my_nonfastwat_bond_dat_midpoint(my_bond_idx)%parm
else
      toler = my_nonfastwat_bond_dat(my_bond_idx)%parm
endif

      diff = - rvdot

!      if (abs(diff) .lt. toler * tol_lcl) cycle
!
!      do  m = 1, 3
!if(usemidpoint) then
!        xij(m) = x(m, i) - x(m, j) - mult_vec_shake_bonds(m,my_bond_idx) * pbc_box(m)
!else
!        xij(m) = x(m, i) - x(m, j)
!endif
!      end do
!
!      if (ns .ne. 0) then
!        call wrap_crds(xij)
!      end if

! Shake resetting of coordinate is done here:

!      rrpr = zero
!
!      do  m = 1, 3
!        rrpr = rrpr + xij(m) * xpij(m)
!      end do
!
!      if (rrpr .lt. toler * 1.0d-06) then ! Deviation too large.  Kill PMEMD.
!        if (master) then
!          write(mdout, 321) iter_cnt, my_bond_idx, i, j
!          if (imin .eq. 1) write(mdout, 331)
!        end if
!        call mexit(6, 1)
!      end if

if(usemidpoint) then
#ifdef MPI
      winvi = 1.d0/proc_atm_mass(i)
      winvj = 1.d0/proc_atm_mass(j)
#endif
else
      winvi = atm_mass_inv(i)
      winvj = atm_mass_inv(j)
endif

! If belly option is on then resetting of the frozen atom is to be prevented:

      if (ibelly_lcl .gt. 0) then
        if (atm_igroup(i) .le. 0) winvi = zero
        if (atm_igroup(j) .le. 0) winvj = zero
      end if

!if(usemidpoint) then
!#ifdef MPI
!      acor = diff / (rrpr * (1.d0/proc_atm_mass(i) + 1.d0/proc_atm_mass(j) + &
!             1.d0/proc_atm_mass(i) + 1.d0/proc_atm_mass(j)))
!#endif
!else
!      acor = diff / (rrpr * (atm_mass_inv(i) + atm_mass_inv(j) + &
!             atm_mass_inv(i) + atm_mass_inv(j)))
!endif

      acor = diff / (toler * (winvi + winvj))
      if (abs(acor) < tol_lcl) cycle
    

      do m = 1, 3
        xh = xij(m) * acor
        v(m, i) = v(m, i) + xh * winvi
        v(m, j) = v(m, j) - xh * winvj
!        xp(m, i) = xp(m, i) + xh * winvi
!if(usemidpoint) then
!        xp(m, j) = xp(m, j) - xh * winvj ! - mult_vec_shake_bonds(m,my_bond_idx) * pbc_box(m)
!else
!        xp(m, j) = xp(m, j) - xh * winvj
!endif
      end do

      skips(i) = iter_cnt + 1 ! Process in this and the next iteration.
      skips(j) = iter_cnt + 1 ! Process in this and the next iteration.

      done = .false.

    end do
#ifdef MPI
#ifdef GBorn
    !$omp end parallel
end if
#endif
#endif

    if (done) return    ! Normal exit

  end do

  ! We failed to accomplish coordinate resetting.  Kill PMEMD.

  if (master) write(mdout, 311)

  if (imin .eq. 1 .and. master) write(mdout, 331)

  call mexit(6, 1)

311 format(/5x, 'Coordinate resetting (shake) was not accomplished', &
           /5x, 'within 3000 iterations')
321 format(/5x, 'Coordinate resetting cannot be accomplished,', &
           /5x, 'deviation is too large', &
           /5x, 'iter_cnt, my_bond_idx, i and j are :', 4i8)
331 format(1x, ' *** Especially for minimization, try ntc=1 (no shake)')
