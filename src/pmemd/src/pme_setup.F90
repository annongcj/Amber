#include "copyright.i"

!*******************************************************************************
!
! Module: pme_setup_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module pme_setup_mod

  use gbl_datatypes_mod

  implicit none

contains

!*******************************************************************************
!
! Subroutine:  final_pme_setup
!
! Description:  Called after dynamic memory allocations to fill in some arrays
!               and perform other initial chores.  Called by ALL processes.
!*******************************************************************************

subroutine final_pme_setup(num_ints, num_reals)

  use constraints_mod
  use pme_recip_dat_mod
  use pmemd_lib_mod
  use img_mod
  use pbc_mod
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use nb_pairlist_mod
  use parallel_dat_mod
  use ene_frc_splines_mod
  use prmtop_dat_mod
  use timers_mod
  use ti_mod
  use file_io_dat_mod
  
  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  double precision              :: tim1, tim2

#ifdef MPI
  integer                       :: alloc_failed
#endif

  ! check_neutral() is currently executed by all processes.  This is not
  ! really necessary but we leave it here for now to keep output consistent
  ! with sander.

  if (ti_mode .eq. 0) then
    call check_neutral
  else
    call ti_check_neutral
  end if

#if defined(CUDA) 
  if (ti_mode .eq. 0) then
    call gpu_upload_charges(atm_qterm)
  else
#if defined(GTI)
    ! TL: Upload charges to make sure consistency when compared w/ CPU results  
    call gpu_upload_charges(atm_qterm)
#endif
  endif  
#endif

  call second(tim1)

  call fill_eed_table

  if(igb .ne. 0 .or. use_pme .ne. 0 .or. ips .ne. 0) call fill_ef_spline_table

  call vdw_correct_setup
  call fill_tranvec(gbl_tranvec)

#ifdef MPI
  ! If all tasks are not both reciprocal and direct tasks, then the direct
  ! workload will be assymetric.  Oversize the pairlist by a factor of 2 so
  ! that we are not slowly growing it instead.  We do not do this for the
  ! special case of fft loadbalancing on 2 cpu's - the potential memory cost
  ! could cause some simulations to not run.

  if (recip_numtasks .ne. numtasks .and. numtasks .gt. 2) then
  if (allocated(gbl_ipairs))  deallocate(gbl_ipairs)
    ipairs_maxsize = 2 * ipairs_maxsize
    allocate(gbl_ipairs(ipairs_maxsize), stat = alloc_failed)
    if (alloc_failed .ne. 0) then
      call alloc_error('final_pme_setup', 'ipairs array reallocation failed!');
    end if
  end if
#endif

  call second(tim2)

#ifdef CUDA
  call gpu_neighbor_list_setup(atm_numex, gbl_natex, vdw_cutoff, &
                               skinnb, skin_permit, smbx_permit, tishake)
#endif

  return

end subroutine final_pme_setup

!*******************************************************************************
!
! Subroutine:  check_neutral
!
! Description:
!
! In general, the Ewald method is only truly applicable under conditions of
! charge neutrality.  When the system is not net neutral, the direct sum and
! reciprocal sums are not "beta" independent.  Regardless, the Ewald method can
! be applied with the ficticious assumption that there is a "uniform net
! neutralizing plasma".  
!
! This routine will remove any net charge resulting from conversion of the low
! precision parm topology charges in the case that the system is supposed to
! be net neutral.
!              
!*******************************************************************************

subroutine check_neutral

  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none

  integer               i
  double precision      sum

  sum = 0.d0

  do i = 1, natom
    sum = sum + atm_qterm(i)
  end do

  if (master) write(mdout, '(/,5x,a,f12.8)') &
    'Sum of charges from parm topology file = ',  sum / 18.2223d0
    
  if (abs(sum/18.2223d0) .gt. 0.01) then
    if (master) &
      write(mdout, '(5x,a,/)') 'Assuming uniform neutralizing plasma'
  else
    if (master) &
      write(mdout, '(5x,a,/)') 'Forcing neutrality...'
    sum = sum / natom
    do i = 1, natom
      atm_qterm(i) = atm_qterm(i) - sum
    end do

  end if

  return

end subroutine check_neutral

!*******************************************************************************
!
! Subroutine:  fill_eed_table
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine fill_eed_table

  use bspline_mod
  use pme_direct_mod
  use mdin_ewald_dat_mod

  implicit none

  double precision      del, x, switch, d_switch_dx

  integer               i, ibcbeg, ibcend

  double precision, dimension(mxeedtab) :: tau

  del = 1.d0 / eedtbdns

  x = 0.d0
  call get_ee_func(x, switch, d_switch_dx)
  gbl_eed_cub(2) = d_switch_dx
  x = (mxeedtab - 1) * del
  call get_ee_func(x, switch, d_switch_dx)
  gbl_eed_cub(2 + 4 * (mxeedtab - 1)) = d_switch_dx

  do i = 1, mxeedtab
    x = del * (i - 1)
    call get_ee_func(x, switch, d_switch_dx)
    tau(i) = x
    gbl_eed_cub(1 + 4 * (i - 1)) = switch
  end do

  ibcbeg = 1
  ibcend = 1
  call cubspl(tau, gbl_eed_cub, mxeedtab, ibcbeg, ibcend)

  return
  
end subroutine fill_eed_table

!*******************************************************************************
!
! Subroutine:  chk_switch
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine chk_switch(max_erfc_relerr)

  use pme_direct_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  double precision, intent(out) :: max_erfc_relerr
  
! Local variables:

  double precision      switch
  double precision      aswitch
  double precision      d_switch_dx
  double precision      ad_switch_dx
  double precision      err
  double precision      max_err
  double precision      max_err_at
  double precision      derr
  double precision      max_derr
  double precision      max_derr_at
  double precision      df_err
  double precision      max_df_err
  double precision      max_df_err_at

  double precision      del
  double precision      dx
  double precision      dxdr
  double precision      half
  integer               i, j
  integer               ind
  integer               max, min
  double precision      small
  double precision      third
  double precision      x
  double precision      df_term
  double precision      adf_term

  half = 0.5d0
  third = 1.d0/3.d0
  dxdr = ew_coeff

  if (master) then
    write(mdout, *) &
     '---------------------------------------------------'
  end if

  small = 1.d-4
  max = dxdr * es_cutoff * eedtbdns
  min = dxdr * eedtbdns
  min = 1
  max_err = 0.d0
  max_derr = 0.d0
  max_df_err = 0.d0

  max_err_at = - 1.d0
  max_derr_at = - 1.d0
  max_df_err_at = - 1.d0

  del = 1.d0 / eedtbdns

  if (master) then
    if (ips.gt.0)then
      write(mdout,'(a)')'     eedmeth=6: Using IPS method for electrostatic energy'
    else
      write(mdout, 1000)
      write(mdout, 1001) eedtbdns
      write(mdout, 1002)
    end if
  end if

1000  format(1x, 'APPROXIMATING switch and d/dx switch', &
                 ' using CUBIC SPLINE INTERPOLATION')

1001  format(1x, 'using ', f8.1, ' points per unit in tabled values')

1002  format(1x, 'TESTING RELATIVE ERROR over r ranging from ', '0.0 to cutoff')

  do i = min, max
    do j = 1, 10
      x = del * (i - 1) + j * 0.1d0 * del
      call get_ee_func(x, switch, d_switch_dx)

      ! cubic spline on switch (eedmeth 1) is all we do:

      ind = eedtbdns * x

      dx = x - ind * del

      ind = 4 * ind

      aswitch = gbl_eed_cub(1 + ind) + dx * (gbl_eed_cub(2 + ind) + &
                dx * (gbl_eed_cub(3 + ind) + &
                dx * gbl_eed_cub(4 + ind) * third)*0.5d0)

      ad_switch_dx = gbl_eed_cub(2 + ind) + dx * (gbl_eed_cub(3 + ind) + &
                     dx * gbl_eed_cub(4 + ind) * half)


      err = abs(aswitch - switch) / (abs(switch) + small)

      ! The relative energy error is the same as err above.

      derr = abs(d_switch_dx - ad_switch_dx) /  (abs(d_switch_dx) + small)

      ! To get the relative force error, you have to consider how df is
      ! calculated using the two terms that have error (aswitch, ad_switch_dx).
      ! You can cancel terms across the difference, so there is no need to
      ! actually calculate df per se.  Note x is dxdr * delr.

      df_term = switch - d_switch_dx * x
      adf_term = aswitch - ad_switch_dx * x
      df_err = abs(df_term - adf_term) / (abs(df_term) + small)

      if (err .gt. max_err) then
        max_err = err
        max_err_at = x
      end if

      if (derr .gt. max_derr) then
        max_derr = derr
        max_derr_at = x
      end if

      if (df_err .gt. max_df_err) then
        max_df_err = df_err
        max_df_err_at = x
      end if

    end do
  end do
  
  if (master) then
    write(mdout, 1003)max_err, max_err_at
    write(mdout, 1004)max_derr, max_derr_at
    write(mdout, *) '---------------------------------------------------'
  end if

  if (max_df_err .gt. max_err) then
    max_erfc_relerr = max_df_err
  else
    max_erfc_relerr = max_err
  end if

! write(mdout, *)'DBG: max df rel err =', max_df_err, ' at', max_df_err_at

1003  format('| CHECK switch(x): max rel err = ', e12.4, '   at ', f10.6)
1004  format('| CHECK d/dx switch(x): max rel err = ', e12.4, '   at ', f10.6)

  return

end subroutine chk_switch

!*******************************************************************************
!
! Subroutine:  get_ee_func
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine get_ee_func(x, switch, d_switch_dx)

  use gbl_constants_mod

  implicit none

  double precision      x, switch, d_switch_dx

! Get switch function multiplying the Coulomb interaction 1/r.
! r has been converted to x by x = dxdr*r for convenience (erfc for ewald):

  call derfcfun(x, switch)

  d_switch_dx = - 2.d0 * exp(-x * x)/sqrt(PI)

  return

end subroutine get_ee_func

!*******************************************************************************
!
! Subroutine:  vdw_correct_setup
!
! Description: <TBS>
!              
! Sets up the numbers of atoms in each vdw type.  Used for
! analytic pressure, energy, correction to vdw dispersion.
!
!*******************************************************************************

subroutine vdw_correct_setup

  use pme_force_mod, only   : gbl_nvdwcls
  use prmtop_dat_mod, only  : natom, atm_iac

  implicit none

  integer       n, j

  gbl_nvdwcls(:) = 0

  do n = 1, natom
    j = atm_iac(n)
    gbl_nvdwcls(j) = gbl_nvdwcls(j) + 1
  end do

  return

end subroutine vdw_correct_setup

end module pme_setup_mod
