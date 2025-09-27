  use prmtop_dat_mod , only : ntypes, gbl_one_scee, gbl_one_scnb
  use parallel_dat_mod
  use processor_mod,only:index_14, box_len

  implicit none

! Formal arguments:
  double precision, intent(in)                  :: charge(*)
  double precision, intent(in)                  :: crd(3, *)
  double precision, intent(in out)              :: frc(3, *)
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

  ee14 = 0.d0
  enb14 = 0.d0

  do_virial = present(e14vir)

  if (do_virial) e14vir(:, :) = 0.d0

  if (nb14_cnt .eq. 0) return

  ntypes_lcl = ntypes

  do n_out = 1, nb14_cnt, 64
  nmax = n_out + 63
  if(nmax .gt. nb14_cnt) nmax = nb14_cnt
  do n = n_out, nmax
    i = nb14(1, n)
    j = nb14(2, n)
    parm_idx = nb14(3, n)
  !  if(parm_idx .gt. 33) write(0,*)"WARNING:",n,parm_idx
    scnb0 = gbl_one_scnb(parm_idx)
    dx = crd(1, j) + index_14(1,n) * box_len(1) - crd(1, i)
    dy = crd(2, j) + index_14(2,n) * box_len(2) - crd(2, i)
    dz = crd(3, j) + index_14(3,n) * box_len(3) - crd(3, i)
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

  end do
  end do

  if (do_virial) then
    e14vir(2, 1) = e14vir(1, 2)
    e14vir(3, 1) = e14vir(1, 3)
    e14vir(3, 2) = e14vir(2, 3)
  end if

  return

