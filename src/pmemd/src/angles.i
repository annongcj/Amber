  use nmr_calls_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use ti_mod
#ifdef _OPENMP_
  use omp_lib
#endif
  use mdin_ctrl_dat_mod
#ifdef TIDECOMP
  use ti_decomp_mod
#endif
  implicit none

! Formal arguments:

  integer               :: angle_cnt
  type(angle_rec)       :: angle(*)
  double precision      :: x(3, *)
  double precision      :: frc(3, *)
  double precision      :: angle_energy

! Local variables:

  GBFloat, parameter   :: pt999 = Point999 

  GBFloat      :: ant
  GBFloat      :: cst
  GBFloat      :: dfw
  double precision      :: eaw
  GBFloat      :: rij, rik, rkj
  GBFloat      :: xij, yij, zij
  GBFloat      :: xkj, ykj, zkj
  GBFloat      :: cii, cik, ckk
  GBFloat      :: da
  GBFloat      :: df
  GBFloat      :: dt1, dt2, dt3, dt4, dt5, dt6
  GBFloat      :: sth
  integer               :: i, j, k, ic, jn, jmax,j_out
  
  integer               :: atm_cnt
#ifdef _OPENMP_
  double precision      :: eadev
#endif
  
#ifdef GBTimer
  call get_wall_time(wall_sec, wall_usec)
  start_time_ms = dble(wall_sec) * 1000.0d0 + dble(wall_usec) / 1000.0d0
#endif

  angle_energy = 0.d0

#ifdef MPI
#ifdef GBorn
  atm_cnt = size(gbl_atm_owner_map)
  eadev_arr(:) = 0.d0
  angle_energy_arr(:) = 0.d0
if(master) then
  ! Grand loop for the angle stuff:
  !$omp parallel default(private) &
  !$omp& firstprivate(ti_region)&
  !$omp& shared (angle_energy,x,angle,frc,omplk, &
  !$omp& gbl_teq,gbl_tk,ti_mode,ti_ene,&
  !$omp& ti_lst,ti_sc_lst,ti_weights,angle_cnt, &
  !$omp& eadev_arr, angle_energy_arr) 
  !$omp do
#endif
#endif
  do j_out = 1, angle_cnt, 64
  jmax = j_out + 63
  if(jmax .gt. angle_cnt) jmax = angle_cnt
#ifdef __INTEL_COMPILER
 !$omp simd 
#endif
  do jn = j_out, jmax 
    i = angle(jn)%atm_i
    j = angle(jn)%atm_j
    k = angle(jn)%atm_k
    ic = angle(jn)%parm_idx

! Calculation of the angle:

    xij = x(1, i) - x(1, j)
    xkj = x(1, k) - x(1, j)

    yij = x(2, i) - x(2, j)
    ykj = x(2, k) - x(2, j)

    zij = x(3, i) - x(3, j)
    zkj = x(3, k) - x(3, j)

    rij = xij * xij + yij * yij + zij * zij
    rkj = xkj * xkj + ykj * ykj + zkj * zkj

    rik = sqrt(rij * rkj)

    cst = min(pt999, max(-pt999, (xij * xkj + yij * ykj + zij * zkj) / rik))

    ant = acos(cst)

! Calculation of the energy and deriv:

    da = ant - gbl_teq(ic)

! If eadev is ever supported under mpi, you will need to split the contribution.
#ifdef MPI
#ifdef GBorn
    eadev_arr(i) = eadev_arr(i) + da * da         ! For rms deviation from ideal angles.
#endif
#else
    eadev = eadev + da * da         ! For rms deviation from ideal angles.
#endif

    df = gbl_tk(ic) * da
    eaw = df * da
    if (ti_mode .ne. 0) then
      ti_region = 0
      if (ti_lst(1,i) .eq. 1 .or. ti_lst(1,j) .eq. 1 &
          .or. ti_lst(1,k) .eq. 1) then
        ti_region = 1
      else if (ti_lst(2,i) .eq. 1 .or. ti_lst(2,j) .eq. 1 &
               .or. ti_lst(2,k) .eq. 1) then
        ti_region = 2
      end if
      if (ti_region .gt. 0) then 
        if (ti_mode .ne. 1) then
          if(gti_bat_sc .gt. 0) then
              ti_lscale = (angle(jn)%bat_sc .eq. 1)
          else
              ti_lscale = (ti_sc_lst(i) .lt. 2 .and. &
                           ti_sc_lst(j) .lt. 2 .and. &
                           ti_sc_lst(k) .lt. 2)

          end if
        end if
#ifdef MPI
#ifndef GBorn
        if (gbl_atm_owner_map(i) .eq. mytaskid) then
#endif 
#endif
          if (ti_mode .eq. 1 .or. ti_lscale) then
#ifdef MPI
#ifdef GBorn
            !$OMP CRITICAL
#endif
#endif
            call ti_update_ene(eaw, si_angle_ene, ti_region, 2)

#ifdef TIDECOMP
            lig_dvdl_decomp(id_atom_region(i)+id_atom_region(j)) = &
               lig_dvdl_decomp(id_atom_region(i)+id_atom_region(j)) + &
               eaw * ti_item_dweights(2, ti_region)
 
            lig_ti_decomp(id_atom_region(i)+id_atom_region(j),i) = &
               lig_ti_decomp(id_atom_region(i)+id_atom_region(j),i) + &
               eaw * ti_item_dweights(2, ti_region) / 3.0
            lig_ti_decomp(id_atom_region(i)+id_atom_region(j),j) = &
               lig_ti_decomp(id_atom_region(i)+id_atom_region(j),j) + &
               eaw * ti_item_dweights(2, ti_region) / 3.0
            lig_ti_decomp(id_atom_region(i)+id_atom_region(j),k) = &
               lig_ti_decomp(id_atom_region(i)+id_atom_region(j),k) + &
               eaw * ti_item_dweights(2, ti_region) / 3.0

            ti_decomp(i) = ti_decomp(i) + eaw * ti_item_dweights(2, ti_region) / 3.0d0
            ti_decomp(j) = ti_decomp(j) + eaw * ti_item_dweights(2, ti_region) / 3.0d0
            ti_decomp(k) = ti_decomp(k) + eaw * ti_item_dweights(2, ti_region) / 3.0d0
#endif
#ifdef MPI
#ifdef GBorn
            !$OMP END CRITICAL
#endif
#endif
            eaw = eaw * ti_item_weights(2,ti_region)
          else
#ifdef MPI
#ifdef GBorn
            !$OMP CRITICAL
#endif
#endif
            ti_ene(ti_region,si_angle_ene) = &
            ti_ene(ti_region,si_angle_ene) + eaw
#ifdef MPI
#ifdef GBorn
            !$OMP END CRITICAL
#endif
#endif
            eaw = 0.0
          end if
#ifdef MPI
#ifndef GBorn
      end if
#endif
#endif
        if (ti_mode .eq. 1 .or. ti_lscale) df = df * ti_item_weights(2,ti_region)
      end if
    end if
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
    
    ! We use atm_i to determine who sums up the energy...
#ifdef MPI
#endif/*MPI*/

#ifdef MPI
#ifdef GBorn
     call omp_set_lock(omplk(i))
      angle_energy_arr(i) = angle_energy_arr(i) + eaw
#else
     if (gbl_atm_owner_map(i) .eq. mytaskid) then
       angle_energy = angle_energy + eaw
#endif
       frc(1, i) = frc(1, i) - dt1
       frc(2, i) = frc(2, i) - dt2
       frc(3, i) = frc(3, i) - dt3
#ifdef GBorn
     call omp_unset_lock(omplk(i))
#else
   end if
#endif


#ifdef GBorn
     call omp_set_lock(omplk(j))
#else
     if (gbl_atm_owner_map(j) .eq. mytaskid) then
#endif
      frc(1, j) = frc(1, j) + dt1 + dt4
      frc(2, j) = frc(2, j) + dt2 + dt5
      frc(3, j) = frc(3, j) + dt3 + dt6
#ifdef GBorn
     call omp_unset_lock(omplk(j))
#else
   end if
#endif


#ifdef GBorn
     call omp_set_lock(omplk(k))
#else
     if (gbl_atm_owner_map(k) .eq. mytaskid) then
#endif

      frc(1, k) = frc(1, k) - dt4
      frc(2, k) = frc(2, k) - dt5
      frc(3, k) = frc(3, k) - dt6
#ifdef GBorn
     call omp_unset_lock(omplk(k))
#else
   end if
#endif
#else
    frc(1, i) = frc(1, i) - dt1
    frc(2, i) = frc(2, i) - dt2
    frc(3, i) = frc(3, i) - dt3
    frc(1, j) = frc(1, j) + dt1 + dt4
    frc(2, j) = frc(2, j) + dt2 + dt5
    frc(3, j) = frc(3, j) + dt3 + dt6
    frc(1, k) = frc(1, k) - dt4
    frc(2, k) = frc(2, k) - dt5
    frc(3, k) = frc(3, k) - dt6

    angle_energy = angle_energy + eaw
#endif
  end do
  end do

#ifdef MPI 
#ifdef GBorn
!$omp end do
!$omp end parallel
    eadev = eadev + sum(eadev_arr(1:atm_cnt))
    angle_energy = angle_energy + sum(angle_energy_arr(1:atm_cnt))
 end if
#endif
#endif

#ifdef GBTimer
  call get_wall_time(wall_sec, wall_usec)
  print *, "Angle wall time =", dble(wall_sec) * 1000.0d0 + &
         dble(wall_usec) / 1000.0d0 - start_time_ms
#endif
  return
