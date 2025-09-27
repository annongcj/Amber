#ifdef TIDECOMP
  use ti_decomp_mod
#endif
  use nmr_calls_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use ti_mod
  use mdin_ctrl_dat_mod

  implicit none

! Formal arguments:

  integer               :: bond_cnt
  type(bond_rec)        :: bond(*)
  double precision      :: x(3, *)
  double precision      :: frc(3, *)
  double precision      :: bond_energy

! Local variables:

  double precision      :: da, da_pull, da_press
  double precision      :: df, df_pull, df_press
  double precision      :: dfw
  integer               :: i, j, ic, jn,jmax,j_out
  double precision      :: lcl_bond_energy
  double precision      :: xa, ya, za
  double precision      :: rij
  double precision      :: xij, yij, zij
#ifdef GBTimer
  call get_wall_time(wall_sec, wall_usec)
  start_time_ms = dble(wall_sec) * 1000.0d0 + dble(wall_usec) / 1000.0d0
#endif


  lcl_bond_energy = 0.0d0

! Grand loop for the bond stuff:

#ifdef MPI
#  ifdef GBorn
 lbe_arr(:) = 0.0d0
 if(master) then
!$omp parallel default(none) &
  !$omp& private (jn,jmax,j_out,i,j,xij,yij,zij,rij,ic,da,df,ti_lscale, &
  !$omp& ti_region,dfw,xa,ya,za,da_press,da_pull,df_press,df_pull) &
  !$omp& shared (bond_cnt, bond,gbl_req, gbl_rk, ti_mode, ti_lst, &
  !$omp& ti_sc_lst, ti_ene, ti_weights,frc,omplk,lbe_arr,gbl_rpresseq,&
  !$omp& gbl_rpulleq,gbl_rpressk,gbl_rpullk)
!$omp do
#  endif
#endif
  do j_out = 1, bond_cnt, 64
  jmax = j_out + 63
  if(jmax .gt. bond_cnt) jmax = bond_cnt
!$dir simd
  do jn = j_out, jmax
    i = bond(jn)%atm_i
    j = bond(jn)%atm_j

! Calculation of the bond vector:

    xij = x(1, i) - x(1, j)
    yij = x(2, i) - x(2, j)
    zij = x(3, i) - x(3, j)

    rij = sqrt(xij * xij + yij * yij + zij * zij)

! Calculation of the energy and deriv:

    ic = bond(jn)%parm_idx
    da = rij - gbl_req(ic)
    if (rij < gbl_rpresseq(ic)) then
      da_press = rij - gbl_rpresseq(ic)
    else
      da_press = 0.0
    end if
    if (rij > gbl_rpulleq(ic)) then
      da_pull = rij - gbl_rpulleq(ic)
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
    
    df       = gbl_rk(ic) * da
    df_press = gbl_rpressk(ic) * da_press
    df_pull  = gbl_rpullk(ic) * da_pull

    if (ti_mode .ne. 0) then
      ti_region = 0
      if (ti_lst(1,i) .eq. 1 .or. ti_lst(1,j) .eq. 1) then
        ti_region = 1
      else if (ti_lst(2,i) .eq. 1 .or. ti_lst(2,j) .eq. 1) then
        ti_region = 2
      end if
      if (ti_region .gt. 0) then
        if (ti_mode .ne. 1) then
          if(gti_bat_sc .gt. 0) then
              ti_lscale = (bond(jn)%bat_sc .eq. 1)
          else
              ti_lscale = (ti_sc_lst(i) .lt. 2 .and. &
                           ti_sc_lst(j) .lt. 2)
          end if
        end if
#ifdef MPI
#  ifdef GBorn
        !$OMP CRITICAL
#  else
        if (gbl_atm_owner_map(i) .eq. mytaskid) then
#  endif    
#endif
          if (ti_mode .eq. 1 .or. ti_lscale) then 
      
            ! Linear scaling
            call ti_update_ene((df * da) + (df_press * da_press) + (df_pull * da_pull), &
                               si_bond_ene, ti_region, 2)

#ifdef TIDECOMP
            lig_ti_decomp(id_atom_region(i)+id_atom_region(j),i) = &
               lig_ti_decomp(id_atom_region(i)+id_atom_region(j),i) + &
               df * da * ti_item_dweights(2, ti_region) / 2.0

            lig_ti_decomp(id_atom_region(i)+id_atom_region(j),j) = &
               lig_ti_decomp(id_atom_region(i)+id_atom_region(j),j) + &
               df * da * ti_item_dweights(2, ti_region) / 2.0

            lig_dvdl_decomp(id_atom_region(i)+id_atom_region(j)) = &
               lig_dvdl_decomp(id_atom_region(i)+id_atom_region(j)) + &
               df * da * ti_item_dweights(2, ti_region)

            ti_decomp(i) = ti_decomp(i) + df * da * ti_item_dweights(2, ti_region) / 2.0d0
            ti_decomp(j) = ti_decomp(j) + df * da * ti_item_dweights(2, ti_region) / 2.0d0
#endif
          else

            ! No need to scale da as we always scale df
            ti_ene(ti_region,si_bond_ene) = ti_ene(ti_region,si_bond_ene) + &
                                            (df * da) + (df_press * da_press) + &
                                            (df_pull * da_pull)
            da = 0.d0
            da_press = 0.d0
            da_pull = 0.d0
          end if
#ifdef MPI
#  ifdef GBorn
       !$OMP END CRITICAL
#  else
        end if
#  endif
#endif
        if (ti_mode .eq. 1 .or. ti_lscale) then
           df       = df * ti_item_weights(2,ti_region)
           df_press = df_press * ti_item_weights(2,ti_region)
           df_pull  = df_pull * ti_item_weights(2,ti_region)
        end if
      end if
    end if
    
    dfw = 2.0 * (df + df_press + df_pull) / rij

! Calculation of the force:

    xa = dfw * xij
    ya = dfw * yij
    za = dfw * zij

#ifdef MPI
    ! We use atm_i to determine who sums up the energy...
#  ifdef GBorn
     call omp_set_lock(omplk(i))
     lbe_arr(i) = lbe_arr(i) + (df * da) + (df_press * da_press) + (df_pull * da_pull)
#  else
    if (gbl_atm_owner_map(i) .eq. mytaskid) then
      lcl_bond_energy = lcl_bond_energy + (df * da) + (df_press * da_press) + &
                        (df_pull * da_pull)
#  endif
      frc(1, i) = frc(1, i) - xa
      frc(2, i) = frc(2, i) - ya
      frc(3, i) = frc(3, i) - za
#  ifdef GBorn
    call omp_unset_lock(omplk(i))
#  else
    end if
#  endif

#  ifdef GBorn
    call omp_set_lock(omplk(j))
#  else
    if (gbl_atm_owner_map(j) .eq. mytaskid) then
#  endif
      frc(1, j) = frc(1, j) + xa
      frc(2, j) = frc(2, j) + ya
      frc(3, j) = frc(3, j) + za
#  ifdef GBorn
      call omp_unset_lock(omplk(j))
#  else
    end if
#  endif
#else
    lcl_bond_energy = lcl_bond_energy + (df * da) + (df_press * da_press) + (df_pull * da_pull)
    frc(1, i) = frc(1, i) - xa
    frc(2, i) = frc(2, i) - ya
    frc(3, i) = frc(3, i) - za
    frc(1, j) = frc(1, j) + xa
    frc(2, j) = frc(2, j) + ya
    frc(3, j) = frc(3, j) + za
#endif

  end do
  end do

#ifdef MPI
#ifdef GBorn
!$omp end do
!$omp end parallel
    bond_energy = sum(lbe_arr) 
 end if
#else
  bond_energy = lcl_bond_energy

#endif
#else
  bond_energy = lcl_bond_energy
#endif


#ifdef GBTimer
  call get_wall_time(wall_sec, wall_usec)
  print *, "Bond wall time =", dble(wall_sec) * 1000.0d0 + &
      dble(wall_usec) / 1000.0d0 - start_time_ms
#endif

  return
