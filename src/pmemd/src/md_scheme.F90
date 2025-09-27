module md_scheme
    implicit none

contains

!*******************************************************************************
!
! Subroutine:  middle_anderson_setvel
!
! Description: <TBS>
!              
!*******************************************************************************
!subroutine middle_anderson_setvel(atm_cnt, vel, mass, mass_inv, &
!                           dt, temp0, gamma_ln)
!  
!  use gbl_constants_mod, only : KB
!  use parallel_dat_mod
!  use random_mod
!  use mdin_ctrl_dat_mod, only : no_ntt3_sync
!  implicit none
!
!! Formal arguments:
!
!  integer               :: atm_cnt
!  double precision      :: vel(3, atm_cnt)
!  double precision      :: mass(atm_cnt)
!  double precision      :: mass_inv(atm_cnt)
!  double precision      :: dt
!  double precision      :: temp0
!  double precision      :: gamma_ln
!
!! Local variables:
!
!  double precision      :: dtx
!  double precision      :: fln1, fln2, fln3
!  double precision      :: half_dtx
!  integer               :: j,i
!
!  double precision      :: lgv_c1, lgv_c2, stdvel, rtkT
!
!  boltz2 = 8.31441d-3 * 0.5d0 / 4.184d0
!  dtx = dt * 20.455d+00
!  half_dtx = dtx * 0.5d0
!
!  lgv_c1 = exp(-gamma_ln * dt)
!  lgv_c2 = sqrt(1.0d0 - lgv_c1 * lgv_c1)
!  rtkT = sqrt(kB * temp0)
!
!! Split here depending on whether we are synching the random
!! number stream across MPI tasks for ntt=3. If ig=-1 then we
!! do not sync. This gives better scaling. We duplicate code here
!! to avoid an if statement in the inner loop.
!#if defined(MPI) && !defined(CUDA)
!  if (no_ntt3_sync) then
!    do j = 1, atm_cnt
!      if (gbl_atm_owner_map(j) .eq. mytaskid) then
!        stdvel = rtkT * sqrt(mass_inv(j))
!        call gauss(0.d0, stdvel, fln1)
!        call gauss(0.d0, stdvel, fln2)
!        call gauss(0.d0, stdvel, fln3)
!        vel(1, j) = lgv_c1 * vel(1, j) + lgv_c2 * fln1
!        vel(2, j) = lgv_c1 * vel(2, j) + lgv_c2 * fln2
!        vel(3, j) = lgv_c1 * vel(3, j) + lgv_c2 * fln3
!      end if
!    end do
!  else
!    do j = 1, atm_cnt
!      if (gbl_atm_owner_map(j) .eq. mytaskid) then
!        stdvel = rtkT * sqrt(mass_inv(j))
!        call gauss(0.d0, stdvel, fln1)
!        call gauss(0.d0, stdvel, fln2)
!        call gauss(0.d0, stdvel, fln3)
!        vel(1, j) = lgv_c1 * vel(1, j) + lgv_c2 * fln1
!        vel(2, j) = lgv_c1 * vel(2, j) + lgv_c2 * fln2
!        vel(3, j) = lgv_c1 * vel(3, j) + lgv_c2 * fln3
!      else
!        call gauss(0.d0, 1.d0)
!        call gauss(0.d0, 1.d0)
!        call gauss(0.d0, 1.d0)
!      end if
!    end do
!  end if
!#else
!  do j = 1, atm_cnt
!    stdvel = rtkT * sqrt(mass_inv(j))
!    call gauss(0.d0, stdvel, fln1)
!    call gauss(0.d0, stdvel, fln2)
!    call gauss(0.d0, stdvel, fln3)
!    vel(1, j) = lgv_c1 * vel(1, j) + lgv_c2 * fln1
!    vel(2, j) = lgv_c1 * vel(2, j) + lgv_c2 * fln2
!    vel(3, j) = lgv_c1 * vel(3, j) + lgv_c2 * fln3
!  end do
!#endif
!  return
!
!end subroutine
!*******************************************************************************
!
! Subroutine:  middle_langevin_thermostat
!
! Description: <TBS>
!              
!*******************************************************************************
subroutine middle_langevin_thermostat(atm_cnt, vel, mass, mass_inv, &
                           dt, temp0, gamma_ln)
  
  use gbl_constants_mod, only : KB
  use parallel_dat_mod
  use random_mod
  use mdin_ctrl_dat_mod, only : no_ntt3_sync
  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: vel(3, atm_cnt)
  double precision      :: mass(atm_cnt)
  double precision      :: mass_inv(atm_cnt)
  double precision      :: dt
  double precision      :: temp0
  double precision      :: gamma_ln

! Local variables:

  double precision      :: dtx
  double precision      :: fln1, fln2, fln3
  double precision      :: half_dtx
  integer               :: j,i

  double precision      :: lgv_c1, lgv_c2, stdvel, rtkT

  !boltz2 = 8.31441d-3 * 0.5d0 / 4.184d0
  dtx = dt * 20.455d+00
  half_dtx = dtx * 0.5d0

  lgv_c1 = exp(-gamma_ln * dt)
  lgv_c2 = sqrt(1.0d0 - lgv_c1 * lgv_c1)
  !kbp = 1.380658 * 6.0221367 / 4.184d3
  !rtkT = sqrt(kbp * temp0)
  rtkT = sqrt(kB * temp0)

! Split here depending on whether we are synching the random
! number stream across MPI tasks for ntt=3. If ig=-1 then we
! do not sync. This gives better scaling. We duplicate code here
! to avoid an if statement in the inner loop.
#if defined(MPI) && !defined(CUDA)
  if (no_ntt3_sync) then
    do j = 1, atm_cnt
      if (gbl_atm_owner_map(j) .eq. mytaskid) then
        stdvel = rtkT * sqrt(mass_inv(j))
        call gauss(0.d0, stdvel, fln1)
        call gauss(0.d0, stdvel, fln2)
        call gauss(0.d0, stdvel, fln3)
        vel(1, j) = lgv_c1 * vel(1, j) + lgv_c2 * fln1
        vel(2, j) = lgv_c1 * vel(2, j) + lgv_c2 * fln2
        vel(3, j) = lgv_c1 * vel(3, j) + lgv_c2 * fln3
      end if
    end do
  else
    do j = 1, atm_cnt
      if (gbl_atm_owner_map(j) .eq. mytaskid) then
        stdvel = rtkT * sqrt(mass_inv(j))
        call gauss(0.d0, stdvel, fln1)
        call gauss(0.d0, stdvel, fln2)
        call gauss(0.d0, stdvel, fln3)
        vel(1, j) = lgv_c1 * vel(1, j) + lgv_c2 * fln1
        vel(2, j) = lgv_c1 * vel(2, j) + lgv_c2 * fln2
        vel(3, j) = lgv_c1 * vel(3, j) + lgv_c2 * fln3
      else
        call gauss(0.d0, 1.d0)
        call gauss(0.d0, 1.d0)
        call gauss(0.d0, 1.d0)
      end if
    end do
  end if
#else

  do j = 1, atm_cnt
    stdvel = rtkT * sqrt(mass_inv(j))
    call gauss(0.d0, stdvel, fln1)
    call gauss(0.d0, stdvel, fln2)
    call gauss(0.d0, stdvel, fln3)
    vel(1, j) = lgv_c1 * vel(1, j) + lgv_c2 * fln1
    vel(2, j) = lgv_c1 * vel(2, j) + lgv_c2 * fln2
    vel(3, j) = lgv_c1 * vel(3, j) + lgv_c2 * fln3
  end do
#endif
  return

end subroutine

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  middle_langevin_thermostat_midpoint
!
! Description: <TBS>
!              
!*******************************************************************************
subroutine middle_langevin_thermostat_midpoint(gbl_atm_cnt, atm_cnt, vel, mass, &
                           dt, temp0, gamma_ln)
  
  use gbl_constants_mod, only : KB
  use parallel_dat_mod
  use random_mod
  use mdin_ctrl_dat_mod, only : no_ntt3_sync
  use processor_mod, only : proc_atm_to_full_list, proc_atm_space, mytaskid
  implicit none

! Formal arguments:
  integer               :: gbl_atm_cnt
  integer               :: atm_cnt
  double precision      :: vel(3, atm_cnt)
  double precision      :: mass(atm_cnt)
  double precision      :: dt
  double precision      :: temp0
  double precision      :: gamma_ln

! Local variables:

  double precision      :: dtx
  double precision      :: fln1, fln2, fln3
  double precision      :: half_dtx
  integer               :: j,i

  double precision      :: lgv_c1, lgv_c2, stdvel, rtkT

  !boltz2 = 8.31441d-3 * 0.5d0 / 4.184d0
  dtx = dt * 20.455d+00
  half_dtx = dtx * 0.5d0

  lgv_c1 = exp(-gamma_ln * dt)
  lgv_c2 = sqrt(1.0d0 - lgv_c1 * lgv_c1)
  rtkT = sqrt(kB * temp0)



! Split here depending on whether we are synching the random
! number stream across MPI tasks for ntt=3. If ig=-1 then we
! do not sync. This gives better scaling. We duplicate code here
! to avoid an if statement in the inner loop.
  if (no_ntt3_sync) then
    do j = 1, atm_cnt
        stdvel = rtkT * sqrt(1.0 / mass(j))
        call gauss(0.d0, stdvel, fln1)
        call gauss(0.d0, stdvel, fln2)
        call gauss(0.d0, stdvel, fln3)
        vel(1, j) = lgv_c1 * vel(1, j) + lgv_c2 * fln1
        vel(2, j) = lgv_c1 * vel(2, j) + lgv_c2 * fln2
        vel(3, j) = lgv_c1 * vel(3, j) + lgv_c2 * fln3
    end do
  else
    do j = 1, gbl_atm_cnt
      i = proc_atm_space(j)
      if (i .ne. 0 .and. i .le. atm_cnt) then
        stdvel = rtkT * sqrt(1.0 / mass(i))
        call gauss(0.d0, stdvel, fln1)
        call gauss(0.d0, stdvel, fln2)
        call gauss(0.d0, stdvel, fln3)

        vel(1, i) = lgv_c1 * vel(1, i) + lgv_c2 * fln1
        vel(2, i) = lgv_c1 * vel(2, i) + lgv_c2 * fln2
        vel(3, i) = lgv_c1 * vel(3, i) + lgv_c2 * fln3
      else
        call gauss(0.d0, 1.d0)
        call gauss(0.d0, 1.d0)
        call gauss(0.d0, 1.d0)
      end if
    end do
  end if


  return

end subroutine
#endif /*MPI*/

end module


