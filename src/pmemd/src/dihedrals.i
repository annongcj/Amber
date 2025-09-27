#ifdef TIDECOMP
  use ti_decomp_mod
#endif
  use gbl_constants_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use ti_mod
#ifdef _OPENMP_
  use omp_lib
#endif

  implicit none

! Formal arguments:

  integer                       :: dihed_cnt
  type(dihed_rec)               :: dihed(*)
  double precision              :: x(3, *)
  double precision              :: frc(3, *)
  double precision              :: ep

! Local variables:

  GBFloat      :: ap
  GBFloat      :: cosnp
  GBFloat      :: cphi
  GBFloat      :: ct, ct0
  GBFloat      :: dc1, dc2, dc3, dc4, dc5, dc6
  GBFloat      :: df
  GBFloat      :: dr1, dr2, dr3, dr4, dr5, dr6
  GBFloat      :: drx, dry, drz
  GBFloat      :: dums
  GBFloat      :: dx, dy, dz
  GBFloat      :: epl
  double precision      :: epw
  GBFloat      :: f1, f2
  GBFloat      :: fxi, fyi, fzi
  GBFloat      :: fxj, fyj, fzj
  GBFloat      :: fxk, fyk, fzk
  GBFloat      :: fxl, fyl, fzl
  GBFloat      :: g
  GBFloat      :: gmul(10)
  GBFloat      :: gx, gy, gz
  integer               :: ic, ic0
  integer               :: inc
  integer               :: i, j, k, kt, l, lt
  integer               :: jn, j_out,jmax
  GBFloat      :: one
  GBFloat      :: s
  GBFloat      :: sinnp
  GBFloat      :: sphi
  GBFloat      :: tenm3
  GBFloat      :: tm06, tm24
  GBFloat      :: xa, ya, za
  GBFloat      :: xij, yij, zij
  GBFloat      :: xkj, ykj, zkj
  GBFloat      :: xkl, ykl, zkl
  GBFloat      :: z1, z2
  GBFloat      :: z11, z12, z22
  GBFloat    ::   zero

 
#ifdef SPDP
  data gmul /0.0, 2.0, 0.0, 4.0, 0.0, 6.0, 0.0, 8.0, 0.0, 10.0/
  data tm24, tm06, tenm3/1.e-18, 1.e-06, 1.e-03/
  data zero, one /0.0, 1.0/
#else
  data gmul /0.d0, 2.d0, 0.d0, 4.d0, 0.d0, 6.d0, 0.d0, 8.d0, 0.d0, 10.d0/
  data tm24, tm06, tenm3/1.d-18, 1.d-06, 1.d-03/
  data zero, one /0.d0, 1.d0/
#endif

! Arrays gbl_gamc = gbl_pk * cos(phase) and gbl_gams = gbl_pk * sin(phase)

  epl = zero

! Grand loop for the dihedral stuff:

#ifdef GBTimer
  call get_wall_time(wall_sec, wall_usec)
  start_time_ms = dble(wall_sec) * 1000.0d0 + dble(wall_usec) / 1000.0d0
#endif

#ifdef MPI
#ifdef GBorn
  epl_arr(:) = 0.d0
if(master) then
  !$omp parallel default (none) &
  !$omp& private(jn,j_out,jmax,xij,yij,zij,xkj,ykj,zkj,xkl,ykl,zkl, &
  !$omp& i,j,kt,lt,k,l,dx,dy,dz,gx,gy,gz,ti_lscale,&
  !$omp& dc1,dc2,dc3,dc4,dc5,dc6,dr1,dr2,dr3, &
  !$omp& dr4,dr5,dr6,drx,dry,drz,fxj,fyj,fzj, &
  !$omp& z11,z22,fxk,fyk,fzk,fxl,fyl,fzl,dums,df,cosnp,sinnp, &
  !$omp& fxi,fyi,z12,fzi,s,ap,cphi,sphi,ic,inc,ct,ct0,z1,z2,epw) &
  !$omp& firstprivate(ti_region)&
  !$omp& shared (dihed_cnt,dihed,x,gbl_ipn,gbl_pn,frc,epl,omplk,&
  !$omp& gmul,tm24,tm06,tenm3,zero,one,gbl_gamc,gbl_pk,gbl_gams, &
  !$omp& ti_mode,ti_lst,ti_sc_lst,ti_weights,ti_ene, &
  !$omp& epl_arr, mytaskid, gbl_atm_owner_map) 
  !$omp do 
#endif
#endif

  do j_out = 1, dihed_cnt,64
  jmax = j_out + 63 
  if (jmax .gt. dihed_cnt) jmax = dihed_cnt
#ifdef __INTEL_COMPILER
 !$omp simd
#endif
  do jn = j_out, jmax 
    i = dihed(jn)%atm_i
    j = dihed(jn)%atm_j
    kt = dihed(jn)%atm_k
    lt = dihed(jn)%atm_l
    k = iabs(kt)
    l = iabs(lt)

! Calculation of ij, kj, kl vectors:

    xij = ToGBFloat(x(1, i) - x(1, j))
    yij = ToGBFloat(x(2, i) - x(2, j))
    zij = ToGBFloat(x(3, i) - x(3, j))
    xkj = ToGBFloat(x(1, k) - x(1, j))
    ykj = ToGBFloat(x(2, k) - x(2, j))
    zkj = ToGBFloat(x(3, k) - x(3, j))
    xkl = ToGBFloat(x(1, k) - x(1, l))
    ykl = ToGBFloat(x(2, k) - x(2, l))
    zkl = ToGBFloat(x(3, k) - x(3, l)) 

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
    if (ti_mode .ne. 0) then
      ti_region = 0
      if (ti_lst(1,i) .eq. 1 .or. ti_lst(1,j) .eq. 1 &
          .or. ti_lst(1,k) .eq. 1 .or. ti_lst(1,l) .eq. 1) then
        ti_region = 1
      else if (ti_lst(2,i) .eq. 1 .or. ti_lst(2,j) .eq. 1 &
               .or. ti_lst(2,k) .eq. 1 .or. ti_lst(2,l) .eq. 1) then
        ti_region = 2
      end if

      if (ti_region .gt. 0) then
        if (ti_mode .ne. 1) then
          if(gti_bat_sc .gt. 0) then
              ti_lscale = (dihed(jn)%bat_sc .eq. 1)
          else
              ti_lscale = (ti_sc_lst(i) .lt. 2 .and. &
                           ti_sc_lst(j) .lt. 2 .and. &
                           ti_sc_lst(k) .lt. 2 .and. &
                           ti_sc_lst(l) .lt. 2) 
          end if
        end if

#ifdef MPI
#ifdef GBorn
         !$OMP CRITICAL
#else
         if (gbl_atm_owner_map(i) .eq. mytaskid) then
#endif
#endif
          if (ti_mode .eq. 1 .or. ti_lscale) then !linear scaling 
            call ti_update_ene(epw, si_dihedral_ene, ti_region, 2)

#ifdef TIDECOMP
            lig_dvdl_decomp(id_atom_region(i)+id_atom_region(j)) = &
              lig_dvdl_decomp(id_atom_region(i)+id_atom_region(j)) + &
              epw * ti_item_dweights(2, ti_region)

            lig_ti_decomp(id_atom_region(i)+id_atom_region(j),i) = &
              lig_ti_decomp(id_atom_region(i)+id_atom_region(j),i) + &
              epw * ti_item_dweights(2, ti_region)/4.0
            lig_ti_decomp(id_atom_region(i)+id_atom_region(j),j) = &
              lig_ti_decomp(id_atom_region(i)+id_atom_region(j),j) + &
              epw * ti_item_dweights(2, ti_region)/4.0
            lig_ti_decomp(id_atom_region(i)+id_atom_region(j),k) = &
              lig_ti_decomp(id_atom_region(i)+id_atom_region(j),k) + &
              epw * ti_item_dweights(2, ti_region)/4.0
            lig_ti_decomp(id_atom_region(i)+id_atom_region(j),l) = &
              lig_ti_decomp(id_atom_region(i)+id_atom_region(j),l) + &
              epw * ti_item_dweights(2, ti_region)/4.0

            ti_decomp(i) = ti_decomp(i) + epw * ti_item_dweights(2, ti_region) / 4.0d0
            ti_decomp(j) = ti_decomp(j) + epw * ti_item_dweights(2, ti_region) / 4.0d0
            ti_decomp(k) = ti_decomp(k) + epw * ti_item_dweights(2, ti_region) / 4.0d0
            ti_decomp(l) = ti_decomp(l) + epw * ti_item_dweights(2, ti_region) / 4.0d0
#endif

            epw = epw * ti_item_weights(2,ti_region)
          else
            ti_ene(ti_region,si_dihedral_ene) = &
               ti_ene(ti_region,si_dihedral_ene) + epw
            epw = 0.d0    
          end if
#ifdef MPI
#ifdef GBorn
          !$OMP END CRITICAL
#else
  end if
#endif
#endif

        if (ti_mode .eq. 1 .or. ti_lscale) df = df * ti_item_weights(2,ti_region)
      end if
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


#ifdef GBorn
      call omp_set_lock(omplk(i))
      epl_arr(i) = epl_arr(i) + epw
#else
      if (gbl_atm_owner_map(i) .eq. mytaskid) then
         epl = epl  + epw
#endif
         frc(1, i) = frc(1, i) + fxi
         frc(2, i) = frc(2, i) + fyi
         frc(3, i) = frc(3, i) + fzi
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
         frc(1, j) = frc(1, j) + fxj 
         frc(2, j) = frc(2, j) + fyj
         frc(3, j) = frc(3, j) + fzj
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
         frc(1, k) = frc(1, k) + fxk
         frc(2, k) = frc(2, k) + fyk
         frc(3, k) = frc(3, k) + fzk
#ifdef GBorn 
     call omp_unset_lock(omplk(k))
#else
     end if
#endif


#ifdef GBorn 
     call omp_set_lock(omplk(l))
#else
     if (gbl_atm_owner_map(l) .eq. mytaskid) then
#endif
        frc(1, l) = frc(1, l) + fxl
        frc(2, l) = frc(2, l) + fyl
        frc(3, l) = frc(3, l) + fzl
#ifdef GBorn 
     call omp_unset_lock(omplk(l))
#else
    end if
#endif

#else /*MPI*/
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
#endif /*MPI*/
  end do
  end do
#ifdef MPI
#ifdef GBorn
!$omp end do
!$omp end parallel  
    epl = epl + sum(epl_arr)
 end if
#endif
#endif

#ifdef GBTimer
  call get_wall_time(wall_sec, wall_usec)
  print *, "Dihedral wall time =", dble(wall_sec) * 1000.0d0 + &
        dble(wall_usec) / 1000.0d0 - start_time_ms
#endif

  ep = epl
  
  return

