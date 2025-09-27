#ifdef TIDECOMP
  use ti_decomp_mod
#endif
  use prmtop_dat_mod, only : ntypes, gbl_one_scee, gbl_one_scnb
  use parallel_dat_mod
  use ti_mod
  use phmd_mod
!PHMD
  use mdin_ctrl_dat_mod, only : iphmd, gti_add_sc
!PHMD
  implicit none

! Formal arguments:
#ifdef GBorn
  integer, intent(in)                           :: atm_cnt
#endif
  double precision, intent(in)                  :: charge(*)
  double precision, intent(in)                  :: crd(3, *)
#ifdef GBorn
  double precision, intent(in out)              :: frc(3, atm_cnt)
#else
  double precision, intent(in out)              :: frc(3, *)
#endif
  integer, intent(in)                           :: iac(*)
  integer, intent(in)                           :: ico(*)
  double precision, intent(in)                  :: cn1(*)
  double precision, intent(in)                  :: cn2(*)
  integer, intent(in)                           :: nb14(3, *)
  integer, intent(in)                           :: nb14_cnt
  double precision, intent(out)                 :: ee14
  double precision, intent(out)                 :: enb14
  double precision, optional, intent(out)       :: e14vir(3, 3)
   
! Local variables:

  integer               :: n, i, j, ic, parm_idx, ia1, ia2, ibig, isml
  integer               :: ntypes_lcl
  integer               :: nmax, n_out
  logical               :: do_virial
  double precision      :: dx, dy, dz, r2, r2inv, rinv
  double precision      :: scnb0, g, f6, f12, r6, df
#ifdef GBTimer
  integer                :: wall_s, wall_u
  double precision, save :: strt_time_ms


  call get_wall_time(wall_s, wall_u)
  strt_time_ms = dble(wall_s) * 1000.00 + dble(wall_u) / 1000.00
#endif

#ifdef MPI
#ifdef GBorn 
  ee14_arr(:) = 0.d0
  enb14_arr(:) = 0.d0
  e14vir11(:) = 0.d0
  e14vir12(:) = 0.d0
  e14vir13(:) = 0.d0
  e14vir22(:) = 0.d0
  e14vir23(:) = 0.d0
  e14vir33(:) = 0.d0
#endif
#endif


  ee14 = 0.d0
  enb14 = 0.d0

  do_virial = present(e14vir)

  if (do_virial) e14vir(:, :) = 0.d0

  if (nb14_cnt .eq. 0) return

  ntypes_lcl = ntypes
#ifdef MPI
#ifdef GBorn
 if(master) then
  !$omp parallel default(none) &
  !$omp& private (i,j,n,parm_idx,scnb0,dx,dy,dz,r2,rinv,&
  !$omp& r2inv,r6,ia1,ia2, ibig, isml,f6,f12,df,ti_region,ti_lscale,g, &
  !$omp& ic,nmax,n_out) &
  !$omp& shared (iphmd, nb14_cnt,nb14,gbl_one_scnb,crd,charge,ti_ene,&
  !$omp& gbl_one_scee,iac,cn2,ti_mode, ti_lst, ti_sc_lst, frc,omplk,& 
  !$omp& ti_weights, cn1,do_virial,ee14_arr,enb14_arr,gti_add_sc,&
  !$omp& e14vir11,e14vir12,e14vir13,e14vir22,e14vir23,e14vir33) 
  !$omp do 
#endif
#endif
  do n_out = 1, nb14_cnt, 64
  nmax = n_out + 63
  if(nmax .gt. nb14_cnt) nmax = nb14_cnt
  do n = n_out, nmax
    i = nb14(1, n)
    j = nb14(2, n)
    parm_idx = nb14(3, n)
    scnb0 = gbl_one_scnb(parm_idx)
    dx = crd(1, j) - crd(1, i)
    dy = crd(2, j) - crd(2, i)
    dz = crd(3, j) - crd(3, i)
    r2 = dx * dx + dy * dy + dz * dz
    rinv = sqrt(ONE / r2)
    r2inv = rinv * rinv
    r6 = r2inv * r2inv * r2inv
    g = charge(i) * charge(j) * rinv * gbl_one_scee(parm_idx)
    !  always use the 6-12 parameters, even if 10-12 are available:
    ia1 = iac(i)
    ia2 = iac(j)
    ibig = max0(ia1,ia2)
    isml = min0(ia1,ia2)
    ic = ibig*(ibig-1)/2+isml
    f6 = cn2(ic) * r6
    f12 = cn1(ic) * (r6 * r6)
    df = (g + scnb0 * (TWELVE * f12 - SIX * f6)) * r2inv

    if (ti_mode .ne. 0) then
      ti_region = 0
      if (ti_lst(1,i) .eq. 1 .or. ti_lst(1,j) .eq. 1) then
        ti_region = 1
      else if(ti_lst(2,i) .eq. 1 .or. ti_lst(2,j) .eq. 1) then
        ti_region = 2
      end if

      if (ti_region .gt. 0) then
        if (ti_mode .ne. 1) then
          if(gti_add_sc .eq. 1) then
            ti_lscale = (ti_sc_lst(i) .lt. 2 .and. &
                         ti_sc_lst(j) .lt. 2) .or. &
                        .not. (ti_sc_lst(i) .eq. 2 .and. &
                         ti_sc_lst(j) .eq. 2)
          else
            ti_lscale = (ti_sc_lst(i) .lt. 2 .and. &
                         ti_sc_lst(j) .lt. 2)
          end if
        end if
#ifdef MPI
#ifdef GBorn
      !$OMP CRITICAL
#else
        if (gbl_atm_owner_map(i) .eq. mytaskid) then
#endif
#endif
          if(ti_mode .eq. 1 .or. ti_lscale) then !linear scaling
            call ti_update_ene(g, si_elect_14_ene, ti_region, 4)
            call ti_update_ene((f12 - f6) * scnb0, si_vdw_14_ene, ti_region, 6)

#ifdef TIDECOMP
            lig_dvdl_decomp(id_atom_region(i)+id_atom_region(j)) = &
               lig_dvdl_decomp(id_atom_region(i)+id_atom_region(j)) + &
               g*ti_item_dweights(5, ti_region) + (f12-f6)*scnb0*ti_item_dweights(6, ti_region)

            lig_ti_decomp(id_atom_region(i)+id_atom_region(j),i) = &
               lig_ti_decomp(id_atom_region(i)+id_atom_region(j),i) + &
               g*ti_item_dweights(5, ti_region)/2.0 + (f12-f6)*scnb0*ti_item_dweights(6, ti_region) / 2.0
            lig_ti_decomp(id_atom_region(i)+id_atom_region(j),j) = &
               lig_ti_decomp(id_atom_region(i)+id_atom_region(j),j) + &
               g*ti_item_dweights(5, ti_region)/2.0 + (f12-f6)*scnb0*ti_item_dweights(6, ti_region) / 2.0

            ti_decomp(i) = ti_decomp(i) &
               + g/2.0d0*ti_item_dweights(5, ti_region) + (f12-f6)*scnb0/2.0d0*ti_item_dweights(6, ti_region)
            ti_decomp(j) = ti_decomp(j) &
               + g/2.0d0*ti_item_dweights(5, ti_region) + (f12-f6)*scnb0/2.0d0*ti_item_dweights(6, ti_region)
#endif

            g = g * ti_item_weights(5,ti_region)    
            f12 = f12 * ti_item_weights(6,ti_region)
            f6 = f6 * ti_item_weights(6,ti_region)       
          else
            ti_ene(ti_region,si_vdw_14_ene) = &
               ti_ene(ti_region,si_vdw_14_ene) + (f12 - f6) * scnb0
            ti_ene(ti_region,si_elect_14_ene) = &
               ti_ene(ti_region,si_elect_14_ene) + g
            f6 = 0.d0
            f12 = 0.d0
            g = 0.d0
          end if 
#ifdef MPI
#ifdef GBorn
      !$OMP END CRITICAL
#else
        else if (gbl_atm_owner_map(j) .eq. mytaskid) then
          if(ti_mode .eq. 1 .or. ti_lscale) then !linear scaling
            g = g * ti_item_weights(5,ti_region)    
            f12 = f12 * ti_item_weights(6,ti_region)
            f6 = f6 * ti_item_weights(6,ti_region)       
          else
            f6 = 0.d0
            f12 = 0.d0
            g = 0.d0
          end if
        end if
#endif
#endif
        if (ti_mode .eq. 1 .or. ti_lscale) then
            df = g * r2inv + scnb0 * (TWELVE*f12-SIX*f6) * r2inv
        end if
      end if
    end if

#ifdef MPI
! We use i to determine who sums up the energy...

#ifdef GBorn
      call omp_set_lock(omplk(i))
      ee14_arr(i) = ee14_arr(i) + g
      enb14_arr(i) = enb14_arr(i)  + (f12 - f6) * scnb0
#else
    if (gbl_atm_owner_map(i) .eq. mytaskid) then
      ee14 = ee14 + g
      enb14 = enb14 + (f12 - f6) * scnb0
#endif
!PHMD
        if (iphmd /=0) then
           call phmd14nb(crd,frc,i,j,rinv*gbl_one_scee(parm_idx),f12,f6,scnb0,r2inv)
        end if
!PHMD

      frc(1, i) = frc(1, i) - df * dx
      frc(2, i) = frc(2, i) - df * dy
      frc(3, i) = frc(3, i) - df * dz

      if (do_virial) then
#ifdef GBorn
        e14vir11(i) = e14vir11(i) - df * dx * dx
        e14vir12(i) = e14vir12(i) - df * dx * dy
        e14vir13(i) = e14vir13(i) - df * dx * dz
        e14vir22(i) = e14vir22(i) - df * dy * dy
        e14vir23(i) = e14vir23(i) - df * dy * dz
        e14vir33(i) = e14vir33(i) - df * dz * dz
#else
        e14vir(1, 1) = e14vir(1, 1) - df * dx * dx
        e14vir(1, 2) = e14vir(1, 2) - df * dx * dy
        e14vir(1, 3) = e14vir(1, 3) - df * dx * dz
        e14vir(2, 2) = e14vir(2, 2) - df * dy * dy
        e14vir(2, 3) = e14vir(2, 3) - df * dy * dz
        e14vir(3, 3) = e14vir(3, 3) - df * dz * dz
#endif
      end if

      if (ti_mode .ne. 0) then
        call ti_update_nb_frc(-df * dx, -df * dy ,-df * dz, i, ti_region,4)
      end if

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

      frc(1, j) = frc(1, j) + df * dx
      frc(2, j) = frc(2, j) + df * dy
      frc(3, j) = frc(3, j) + df * dz

      if (ti_mode .ne. 0) then
        call ti_update_nb_frc(df * dx, df * dy, df * dz, j, ti_region,4)
      end if
#ifdef GBorn
   call omp_unset_lock(omplk(j))
#else
    end if
#endif
     
#else
!PHMD
        if (iphmd /=0) then
           call phmd14nb(crd,frc,i,j,rinv*gbl_one_scee(parm_idx),f12,f6,scnb0,r2inv)
        end if
!PHMD

    ee14 = ee14 + g
    enb14 = enb14 + (f12 - f6) * scnb0
    frc(1, i) = frc(1, i) - df * dx
    frc(2, i) = frc(2, i) - df * dy
    frc(3, i) = frc(3, i) - df * dz
    frc(1, j) = frc(1, j) + df * dx
    frc(2, j) = frc(2, j) + df * dy
    frc(3, j) = frc(3, j) + df * dz

    if (do_virial) then
      e14vir(1, 1) = e14vir(1, 1) - df * dx * dx
      e14vir(1, 2) = e14vir(1, 2) - df * dx * dy
      e14vir(1, 3) = e14vir(1, 3) - df * dx * dz
      e14vir(2, 2) = e14vir(2, 2) - df * dy * dy
      e14vir(2, 3) = e14vir(2, 3) - df * dy * dz
      e14vir(3, 3) = e14vir(3, 3) - df * dz * dz
    end if

    if (ti_mode .ne. 0) then
      call ti_update_nb_frc(-df * dx, -df * dy, -df * dz, i, ti_region,4)
      call ti_update_nb_frc(df * dx, df * dy, df * dz, j, ti_region,4)
    end if

#endif /* MPI */
  end do
  end do

#ifdef MPI
#ifdef GBorn
  !$omp end do
  !$omp end parallel
  ee14 = ee14 + sum(ee14_arr)
  enb14 = enb14 + sum(enb14_arr)
  if (do_virial) then
    e14vir(1,1) = e14vir(1,1) + sum(e14vir11)
    e14vir(1,2) = e14vir(1,2) + sum(e14vir12)
    e14vir(1,3) = e14vir(1,3) + sum(e14vir13)
    e14vir(2,2) = e14vir(2,2) + sum(e14vir22)
    e14vir(2,3) = e14vir(2,3) + sum(e14vir23)
    e14vir(3,3) = e14vir(3,3) + sum(e14vir33)
  end if
end if
#endif
#endif

  if (do_virial) then
    e14vir(2, 1) = e14vir(1, 2)
    e14vir(3, 1) = e14vir(1, 3)
    e14vir(3, 2) = e14vir(2, 3)
  end if

#ifdef GBTimer
  call get_wall_time(wall_s, wall_u)
  print *, "nb14 reg ene =", dble(wall_s) * 1000.00 + &
           dble(wall_u) / 1000.00 - strt_time_ms
#endif
  return

