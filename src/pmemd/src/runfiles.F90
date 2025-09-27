#include "copyright.i"
!*******************************************************************************
!
! Module: runfiles_mod
!
! Description: <TBS>
!
!*******************************************************************************

module runfiles_mod

use file_io_dat_mod

  implicit none

contains

!*******************************************************************************
!
! Subroutine:  wrapped_corpac
!
! Description:
!
! We are using a wrapper so that we only allocate the stack storage in crd_copy
! when we need it.
!
!*******************************************************************************

subroutine wrapped_corpac(atm_cnt, crd_cnt, crd, crd_start, nf)

  use mdin_ctrl_dat_mod
  use pbc_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  integer               :: crd_cnt
  double precision      :: crd(atm_cnt * 3)
  integer               :: crd_start
  integer               :: nf

! Local variables:

  double precision      :: crd_copy(atm_cnt * 3)

  crd_copy(:) = crd(:)

  if (ntb .gt. 0) then
    call wrap_molecules(crd_copy)

    if (ifbox .eq. 2 .and. iwrap .eq. 1)  &
      call wrap_to(crd_copy)
  end if

  call corpac(crd_cnt, crd_copy, crd_start, nf)

end subroutine wrapped_corpac

!*******************************************************************************
!
! Subroutine:  corpac
!
! Description: <TBS>
!
!*******************************************************************************

subroutine corpac(iend, crd, istart, nf)

  use mdin_ctrl_dat_mod
  use bintraj_mod

  implicit none

! Formal arguments:

  integer               :: iend
  double precision      :: crd(*)
  integer               :: istart
  integer               :: nf

  if (istart .gt. iend) return

  ! We further wrapper this mess so we can stack-allocate space in the event
  ! of needing to do an unformatted write.  Not sure it is ever even used...

  if (ioutfm .eq. 0) then
    call formatted_corpac(iend, crd, istart, nf)
  else if (ioutfm .eq. 1) then
    if (nf .eq. mdcrd) then
      call write_binary_crds(crd)
    else if (nf .eq. mdvel) then
      call write_binary_vels(crd)
    else if (nf .eq. mdfrc) then
      call write_binary_frcs
    end if
  end if

  return

end subroutine corpac

!*******************************************************************************
!
! Subroutine:  formatted_corpac
!
! Description: <TBS>
!
!*******************************************************************************

subroutine formatted_corpac(iend, crd, istart, nf)

  use axis_optimize_mod

  implicit none

! Formal arguments:

  integer               :: iend
  double precision      :: crd(*)
  integer               :: istart
  integer               :: nf

! Local variables:

  integer, save         :: imax = 0
  integer               :: i
  integer               :: iobuf_cnt
  integer               :: ord1, ord2, ord3
  double precision      :: iobuf(3, 150)

  double precision, parameter   :: rmin = -999.99d0
  double precision, parameter   :: rmax = 9999.99d0
  integer, parameter            :: iobuf_siz = size(iobuf, 2)

  ! It is ASSUMED that iend is a multiple of 3; anything else would be
  ! invalid!  Also istart is a multiple of 3 (including 0) + 1.

! NOTE - corpac input may actually be crds, vels, or box data.

  ord1 = axis_flipback_ords(1) - 1
  ord2 = axis_flipback_ords(2) - 1
  ord3 = axis_flipback_ords(3) - 1

  ! We only check for out-of-range values as long as they have not occurred
  ! before.  We assume that once they have occurred, it is unlikely that
  ! we will get reverse diffusion.

  ! BUGBUG - This check should be unnecessary for vels, as any system with
  !          vel velocities past f8.3 format is well past exploding.

  if (imax .eq. 0) then
    do i = istart, iend
      if (crd(i) .gt. rmax .or. crd(i) .lt. rmin) then
        imax = 1
        exit
      end if
    end do
  end if

  i = istart

  do

    iobuf_cnt = 0

    do

      if (iobuf_cnt .ge. iobuf_siz .or. i .gt. iend) exit

      iobuf_cnt = iobuf_cnt + 1

      iobuf(1, iobuf_cnt) = crd(i + ord1)
      iobuf(2, iobuf_cnt) = crd(i + ord2)
      iobuf(3, iobuf_cnt) = crd(i + ord3)

      i = i + 3

    end do

    if (imax .eq. 0) then
      write(nf, 1000) iobuf(:, 1:iobuf_cnt)
    else
      write(nf, 1001) iobuf(:, 1:iobuf_cnt)
    end if

    if (i .gt. iend) exit

  end do

 1000 format(10f8.3)
 1001 format(10f8.2)

  return

end subroutine formatted_corpac

#ifndef NOXRAY
!----------------------------------------------------------------------------------------------
! formatted_strfac: write formatted structure factors to a file, using 6e13.6 format.
!
! Arguments:
!   sfarray:   complex real_kind array of structure factors
!   istart:    starting point for printing structure factors (not assumed to be a multiple of
!              three as above, and in all likelihood will be 1)
!   iend:      ending point for printing structure factors (range is [istart, iend] inclusive).
!              This is also not assumed to be a multiple of 3.
!   sfout:     output file (this is not like mdout, more like a trajectory)
!----------------------------------------------------------------------------------------------
subroutine formatted_strfac(sfarray, istart, iend, sfout)
  use xray_globals_module, only : real_kind
  implicit none

  ! Formal arguments
  complex(real_kind) :: sfarray(*)
  integer            :: istart, iend, sfout
  
  ! Local variables
  double precision :: iobuf(2, 120)
  integer          :: i, j, jlim

  i = istart
  do while (i .le. iend + 120)
    j = 0
    do while (i + j .le. iend .and. j .lt. 120)
      iobuf(1, j + 1) =  real(sfarray(i + j))
      iobuf(2, j + 1) = aimag(sfarray(i + j))       
      j = j + 1
    end do
    i = i + 120
    if (j .gt. 0) then
      write(sfout, 1002) iobuf(:,1:j)
    end if
  end do

 1002 format(6e13.5)
  
  return
end subroutine formatted_strfac
#endif

!*******************************************************************************
!
! Subroutine:  wrapped_mdwrit
!
! Description:
!
! We are using a wrapper so that we only allocate the stack storage in crd_copy
! when we need it.
!
!*******************************************************************************

subroutine wrapped_mdwrit(nstep, atm_cnt, crd, vel, tt)

  use mdin_ctrl_dat_mod
  use pbc_mod
  use prmtop_dat_mod

  implicit none

  integer               :: nstep
  integer               :: atm_cnt
  double precision      :: crd(3, atm_cnt)
  double precision      :: vel(*)
  double precision      :: tt

! Use temp. array to hold coords. so that the master's values
! are always identical to those on all other nodes:

  double precision      :: crd_copy(3, atm_cnt)

  crd_copy(:,:) = crd(:,:)

  if (ntb .gt. 0) then
    call wrap_molecules(crd_copy)
    if (ifbox .eq. 2) call wrap_to(crd_copy)
  end if

  call mdwrit(nstep, atm_cnt, crd_copy, vel, tt)

  return

end subroutine wrapped_mdwrit

!*******************************************************************************
!
! Internal Subroutine:  mdwrit
!
! Description: <TBS>
!
!*******************************************************************************

subroutine mdwrit(nstep, atm_cnt, crd, vel, tt)

  use file_io_mod
  use mdin_ctrl_dat_mod

  implicit none

! Formal arguments:

  integer               :: nstep
  integer               :: atm_cnt
  double precision      :: crd(3, atm_cnt)
  double precision      :: vel(3, atm_cnt)
  double precision      :: tt

! Local variables:

  integer                 :: istart, iend
  character(12)           :: num
  character(max_fn_len+9) :: restrt2_name

! Write/rewind the restrt:

!We actually open the restart file here and close it afterwards rather than
!using the original method of just rewinding in order to force newer 'broken?'
!linuxes to flush the write buffer. Always open with "Unknown" state so we
!will overwrite the file -- it's already been checked for existence

  if (ntxo .le. 0) then
    call amopen(restrt, restrt_name, 'U', 'U', 'W')
  else if (ntxo .eq. 1) then
    call amopen(restrt, restrt_name, 'U', 'F', 'W')
  end if

  call write_restart(restrt, atm_cnt, crd, vel, tt, .true.,restrt_name)

  if (ntxo .ne. 2) close(restrt)

! Consider whether to save 2ndary restrt:

  if (ntwr .ge. 0) return

  do iend = 1, max_fn_len
    if (restrt_name(iend:iend) .le. ' ') exit
  end do

  iend = iend - 1

  write(num, '(i12)') nstep

  do istart = 1, 12
    if (num(istart:istart) .ne. ' ') exit
  end do

  write(restrt2_name, '(a,a,a)') restrt_name(1:iend), '_', num(istart:12)

  write(mdout, '(a,a)') ' writing ', restrt2_name

  if (ntxo .eq. 0) then
    call amopen(restrt2, restrt2_name, owrite, 'U', 'W')
  else if (ntxo .eq. 1) then
    call amopen(restrt2, restrt2_name, owrite, 'F', 'W')
  end if

  call write_restart(restrt2, atm_cnt, crd, vel, tt, .false.,restrt2_name)

  if (ntxo .ne. 2) close(restrt2)

  return

end subroutine mdwrit

!*******************************************************************************
!
! Subroutine:  mdeng
!
! Description: <TBS>
!
!*******************************************************************************

subroutine mdeng(nstep, time, si, fac, press, virial, ekcmt)

  use axis_optimize_mod
  use mdin_ctrl_dat_mod
  use pbc_mod
  use state_info_mod

  implicit none

! Formal arguments:

  integer               :: nstep
  double precision      :: time
  double precision      :: si(*)        ! State information.
  double precision      :: fac(*)
  double precision      :: press(3)
  double precision      :: virial(3)
  double precision      :: ekcmt(3)

! Local variables:

  integer               :: i
  integer               :: ord1, ord2, ord3
  logical, save         :: first = .true.
  character(16), save   :: labs(42)

  data labs/'Nsteps  ', 'time(ps)  ', 'Etot  ', 'EKinetic  ', &
            'Temp  ', 'T_solute ', 'T_solv  ', 'Pres_scal_solu ', &
            'Pres_scal_solv ', 'BoxX  ', 'BoxY  ', 'BoxZ  ', &
            'volume  ', 'pres_X  ', 'pres_Y  ', 'pres_Z  ', &
            'Pressure ', 'EKCoM_x ', 'EKCoM_y ', 'EKCoM_z', &
            'EKComTot ', 'VIRIAL_x ', 'VIRIAL_y ', 'VIRIAL_z ', &
            'VIRIAL_tot ', 'E_pot  ', 'E_vdw  ', 'E_el  ', &
            'E_hbon  ', 'E_bon  ', 'E_angle  ', 'E_dih  ', &
            'E_14vdw  ', 'E_14el  ', 'E_const  ', 'E_pol  ', &
            'AV_permMoment ', 'AV_indMoment ', 'AV_totMoment ', &
            'Density', 'dV/dlambda', 'surften'/

  ord1 = axis_flipback_ords(1)
  ord2 = axis_flipback_ords(2)
  ord3 = axis_flipback_ords(3)

! Define various terms:

  if (first) then
    ! up to Ekinetic:
    write(mden, 1) 'L0 ', (labs(i), i = 1, 4)
    ! up to Pres_scal_solu:
    write(mden, 1) 'L1 ', (labs(i), i = 5, 8)
    ! up to boxZ:
    write(mden, 1) 'L2 ', (labs(i), i = 9, 12)
    ! up to pres_Z:
    write(mden, 1) 'L3 ', (labs(i), i = 13, 16)
    ! up to EKCoM_z:
    write(mden, 1) 'L4 ', (labs(i), i = 17, 20)
    ! up to VIRIAL_z:
    write(mden, 1) 'L5 ', (labs(i), i = 21, 24)
    ! up to E_el:
    write(mden, 1) 'L6 ', (labs(i), i = 25, 28)
    ! up to E_dih:
    write(mden, 1) 'L7 ', (labs(i), i = 29, 32)
    ! up to E_pol:
    write(mden, 1) 'L8 ', (labs(i), i = 33, 36)
    ! up to Density or dV/dlambda:
    write(mden, 1) 'L9 ', (labs(i), i = 37, 41)
    ! surface tension info if constant surface tension in use.
    if (csurften > 0) &
       write(mden, 1) 'L10 ', labs(42)
1   format(a, 10(1x, a))
    first = .false.
  end if

! Write values for this step:

  ! Pres_scal_solu and Pres_scal_solv are not supported anymore; output values
  ! are fixed at 1.d0

  ! up to Ekinetic:
  write(mden, 2) 'L0 ', nstep, time, si(si_tot_ene), si(si_kin_ene)

  ! up to Pres_scal_solu:
  write(mden, 3) 'L1 ', si(si_kin_ene) / fac(1), &
                        si(si_solute_kin_ene) / fac(2), &
                        si(si_solvent_kin_ene) / fac(3), 1.d0

  ! up to boxZ:
  write(mden, 3) 'L2 ', 1.d0, pbc_box(ord1), pbc_box(ord2), pbc_box(ord3)

  ! up to pres_Z:
  write(mden, 3) 'L3 ', si(si_volume), press(ord1), press(ord2), press(ord3)

  ! up to EKCoM_z:
  write(mden, 3) 'L4 ', si(si_tot_press), ekcmt(ord1), ekcmt(ord2), ekcmt(ord3)

  ! up to VIRIAL_z:
  write(mden, 3) 'L5 ', si(si_tot_ekcmt), &
                        virial(ord1), virial(ord2), virial(ord3)

  ! up to E_el:
  write(mden, 3) 'L6 ', si(si_tot_virial), &
                        si(si_pot_ene), si(si_vdw_ene), si(si_elect_ene)

  ! up to E_dih:
  write(mden, 3) 'L7 ', si(si_hbond_ene), &
                        si(si_bond_ene), si(si_angle_ene), si(si_dihedral_ene)

  ! up to E_pol (currently 0.d0):
  write(mden, 3) 'L8 ', si(si_vdw_14_ene), si(si_elect_14_ene), &
                        si(si_restraint_ene), 0.d0

  ! up to dV/dlambda, includes 3 0.d0 fields:
  write(mden, 3) 'L9 ', 0.d0, 0.d0, 0.d0, si(si_density), si(si_dvdl)

  ! Constant surface tension info if running with constant surface tension
  if (csurften > 0) &
    write(mden, 3) 'L10 ', si(si_gamma_ten)

2 format(a, i8, 20(2x, e16.10))
3 format(a, 20(e16.10, 2x))

  return

end subroutine mdeng

!*******************************************************************************
!
! Subroutine:  write_restart
!
! Description: Routine to write final coordinates and velocities.
!
! EWALD: dump ewald specific box information to the restrt files
!        (may only be necessary with ntx=7).
!
! isMain and ncrst_filename are needed for multiple restart writes with netcdf.
! The isMain argument (which is true if this subroutine is being called for
! the main restart file, false otherwise) is required for netcdf restart
! writes to trigger setup of secondary nc restart files that would not
! occur otherwise (setup for the main restart occurs only once).
! The ncrst_filename argument is needed since netcdf routines cannot make use
! of fortran logical units; otherwise the filename for the 2ndary restart
! would not be known.
!*******************************************************************************

subroutine write_restart(nf, atm_cnt, crd, vel, tt, isMain, ncrst_filename)

  use axis_optimize_mod
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use pbc_mod
  use prmtop_dat_mod
  use binrestart_mod, only: write_nc_restart
#ifdef MPI
  use remd_mod, only : remd_method, remd_dimension, remd_types, &
                       group_num, replica_indexes, remd_repidx, remd_crdidx
#endif

  implicit none

! Formal arguments:

  integer               :: i
  integer               :: ord1, ord2, ord3
  integer               :: nf
  integer               :: atm_cnt
  double precision      :: crd(3, atm_cnt)
  double precision      :: vel(3, atm_cnt)
  double precision      :: tt           ! Not used if minimization
  double precision      :: param
  logical               :: isMain       ! true if main restart, false if 2ndary
  character(len=*)      :: ncrst_filename ! Restart filename, only for netcdf restart
  logical               :: odd_atm_cnt

  ord1 = axis_flipback_ords(1)
  ord2 = axis_flipback_ords(2)
  ord3 = axis_flipback_ords(3)

  odd_atm_cnt = mod(atm_cnt, 2) .ne. 0

  if (ntxo .eq. 2) then                 ! Netcdf write:

#ifdef MPI
    call write_nc_restart(ncrst_filename,prmtop_ititl,owrite,atm_cnt,ntb,isMain,crd,vel,temp0,tt,&
                          imin,pbc_box,pbc_alpha,pbc_beta,pbc_gamma,&
                          remd_method,remd_dimension,remd_types,group_num,replica_indexes,&
                          remd_repidx, remd_crdidx, solvph,solve)
#else
    call write_nc_restart(ncrst_filename,prmtop_ititl,owrite,atm_cnt,ntb,isMain,crd,vel,temp0,tt,&
                          imin,pbc_box,pbc_alpha,pbc_beta,pbc_gamma,solvph,solve)
#endif /*MPI*/

  else if (ntxo .eq. 1) then            ! Formatted writing:
    write(nf, 9008) prmtop_ititl
    if (imin .eq. 0) then
#ifdef MPI
      if (remd_method .eq. 4) then
        param = solvph
      else if (remd_method .eq. 5) then
        param = solve
      else
        param = temp0
      end if
      if (remd_method .ne. 0) then
        if (atm_cnt .lt. 100000) then
          write(nf, 9018) atm_cnt, tt, param
        else if (atm_cnt .lt. 1000000) then
          write(nf, 9019) atm_cnt, tt, param
        else if (atm_cnt .lt. 10000000) then
          write(nf, 9020) atm_cnt, tt, param
        else
          write(nf, 9021) atm_cnt, tt, param
        end if
      else
#endif
        if (atm_cnt .lt. 100000) then
          write(nf, 9018) atm_cnt, tt
        else if (atm_cnt .lt. 1000000) then
          write(nf, 9019) atm_cnt, tt
        else if (atm_cnt .lt. 10000000) then
          write(nf, 9020) atm_cnt, tt
        else
          write(nf, 9021) atm_cnt, tt
        end if
#ifdef MPI
      end if
#endif
      do i = 1, atm_cnt - 1, 2
        write(nf, 9028) crd(ord1, i), crd(ord2, i), crd(ord3, i), &
                        crd(ord1, i+1), crd(ord2, i+1), crd(ord3, i+1)
      end do
      if (odd_atm_cnt) then
        write(nf, 9028) crd(ord1, atm_cnt), &
                        crd(ord2, atm_cnt), &
                        crd(ord3, atm_cnt)
      end if
      do i = 1, atm_cnt - 1, 2
        write(nf, 9028) vel(ord1, i), vel(ord2, i), vel(ord3, i), &
                        vel(ord1, i+1), vel(ord2, i+1), vel(ord3, i+1)
      end do
      if (odd_atm_cnt) then
        write(nf, 9028) vel(ord1, atm_cnt), &
                        vel(ord2, atm_cnt), &
                        vel(ord3, atm_cnt)
      end if
    else
      if (atm_cnt .lt. 100000) then
        write(nf, 9018) atm_cnt
      else if (atm_cnt .lt. 1000000) then
        write(nf, 9019) atm_cnt
      else if (atm_cnt .lt. 10000000) then
        write(nf, 9020) atm_cnt
      else
        write(nf, 9021) atm_cnt
      end if

      do i = 1, atm_cnt - 1, 2
        write(nf, 9028) crd(ord1, i), crd(ord2, i), crd(ord3, i), &
                        crd(ord1, i+1), crd(ord2, i+1), crd(ord3, i+1)
      end do
      if (odd_atm_cnt) then
        write(nf, 9028) crd(ord1, atm_cnt), &
                        crd(ord2, atm_cnt), &
                        crd(ord3, atm_cnt)
      end if
    end if
    if (ntb .ne. 0) write(nf, 9028) pbc_box(ord1),pbc_box(ord2),pbc_box(ord3), &
                                    pbc_alpha,pbc_beta,pbc_gamma

  else                                  ! Binary writing:

    write(nf) prmtop_ititl
    if (imin .eq. 0) then
      write(nf) atm_cnt, tt
      write(nf) (crd(ord1, i), crd(ord2, i), crd(ord3, i), i = 1, atm_cnt)
      write(nf) (vel(ord1, i), vel(ord2, i), vel(ord3, i), i = 1, atm_cnt)
    else
      write(nf) atm_cnt, 0.d0
      write(nf) (crd(ord1, i), crd(ord2, i), crd(ord3, i), i = 1, atm_cnt)
    end if

    ! Sander does not provide box/angle info here, which strikes me as a bug.

    if (ntb .ne. 0) write(nf) pbc_box(ord1), pbc_box(ord2), pbc_box(ord3), &
                              pbc_alpha, pbc_beta, pbc_gamma
  end if

 9008 format(a80)
 9018 format(i5, 5e15.7)
 9019 format(i6, 5e15.7)
 9020 format(i7, 5e15.7)
 9021 format(i8, 5e15.7)
 9028 format(6f12.7)

  return

end subroutine write_restart

!*******************************************************************************
!
! Subroutine:  prntmd
!
! Description: <TBS>
!
!*******************************************************************************

subroutine prntmd(nstep, total_nstlim, time, si, fac, iout7, rms, mdloop)

  use charmm_mod, only : charmm_active
  use file_io_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use nmr_calls_mod
  use parallel_dat_mod
  use prmtop_dat_mod, only : cmap_term_count
#ifdef MPI
  use remd_mod, only   : remd_method, remd_dimension, group_num, replica_indexes
#endif
! Self-Guided molecular/Langevin Dynamics (SGLD)
  use sgld_mod, only : sgld_print,sgfti,sgffi
#ifndef MPI
#ifndef NOXRAY
  use xray_interface_module, only: xray_write_md_state
  use xray_globals_module, only: xray_active
#endif
#endif
  use state_info_mod
  use timers_mod, only : print_ongoing_time_summary
  use ti_mod
  use get_cmdline_mod, only : cpein_specified

#ifdef MPI
  use neb_mod !added by DG
#endif
  use sams_mod

#ifdef GTI  
  use gti_mod
#endif

  implicit none

! Formal arguments:

  integer               :: nstep
  integer               :: total_nstlim
  double precision      :: time
  double precision      :: si(*)
  double precision      :: fac(*)
  integer               :: iout7
  integer               :: mdloop
  logical               :: rms

! Local variables:

  double precision      :: epol
  double precision      :: aveper
  double precision      :: aveind
  double precision      :: avetot
  double precision      :: temp

  logical, save         :: first_6call = .true.
  integer, save         :: next_6flush_sec
  logical, save         :: first_7call = .true.
  integer, save         :: next_7flush_sec
  integer               :: current_sec
  integer               :: current_usec           ! Dummy, not used.
  integer               :: i                      ! loop counter

! Define various terms; most of these are placeholders we don't currently
! support, but that we would have to support in a polarizable ff.

  temp        = si(si_kin_ene) / fac(1)

  if (ti_mode .ne. 0) then
    temp = ti_temp
  end if

  epol        = 0.d0
  aveper      = 0.d0
  aveind      = 0.d0
  avetot      = 0.d0

  write(mdout, 9021) nstep, time, temp, si(si_tot_press)
  write(mdout, 9029) si(si_tot_ene), si(si_kin_ene), si(si_pot_ene)
  write(mdout, 9039) si(si_bond_ene), si(si_angle_ene), si(si_dihedral_ene)
  if (charmm_active) then
    write(mdout, 9040) si(si_angle_ub_ene), si(si_dihedral_imp_ene), &
                       si(si_cmap_ene)
  ! if CMAP is active
  else if (cmap_term_count > 0) then
    write(mdout, 9040) 0.d0,0.d0, si(si_cmap_ene)
  endif
  write(mdout, 9049) si(si_vdw_14_ene), si(si_elect_14_ene), si(si_vdw_ene)

  if (using_pme_potential) then
    write(mdout, 9061) si(si_elect_ene), si(si_hbond_ene), si(si_restraint_ene)
  else
    write(mdout, 9062) si(si_elect_ene), si(si_hbond_ene), si(si_restraint_ene)
  end if

  if(efx .ne. 0 .or. efy .ne. 0 .or. efz .ne. 0) then
    write(mdout, 9091) si(si_efield_ene)
  end if
#ifndef MPI
#ifndef NOXRAY
  if (xray_active) call xray_write_md_state(6)
#endif
#endif
  if (si(si_restraint_ene) .ne. 0.0) &
    write(mdout, 9079) si(si_pot_ene) - si(si_restraint_ene)

  if (iamd.gt.0) &
    write(mdout, 9080) si(si_amd_boost)

  if (rms .and. ntt .eq. 0) write(mdout,9075) si(si_solvent_kin_ene)

  if (ntp .gt. 0 .or. si(si_volume) .ne. 0.0) &
    write(mdout, 9081) si(si_tot_ekcmt), si(si_tot_virial), si(si_volume)

  if (csurften .gt. 0) &
    write(mdout, 9082) si(si_gamma_ten)

  if (epol .ne. 0.0) write(mdout, 9070) epol

  if (ntp .gt. 0 .or. si(si_volume) .ne. 0.0) write(mdout, 9083) si(si_density)

  if (gbsa .eq. 1 .or. gbsa .eq. 3) &
    write(mdout, 9063) si(si_surf_ene)

  if (ti_mode.ne.0) then
    write(mdout, 9089) si(si_dvdl)
  else if (si(si_dvdl).ne.0 ) then
    write(mdout, 9089) si(si_dvdl)
  endif

  if (igamd.gt.0) then
    if(.not.((igamd.ge.12.and.igamd.le.27).or. &
            (igamd.ge.110 .and. igamd.le.120)))then
      write(mdout, 9090) si(si_gamd_boost)
    else
      write(mdout, 9092) si(si_gamd_boost),si(si_gamd_ppi)
    endif
  endif

#ifndef CUDA
  !Ewald error estimate is not calculated when running on GPUs.
  !Do not print error estimate when running IPS.
  if (using_pme_potential .and. ips == 0) write(mdout, 9188) si(si_pme_err_est)
#endif

!  Printout SGLD guiding information
if(.not. usemidpoint) then
   if (isgld > 0) call sgld_print(mdout)
endif

#ifdef MPI
  if (remd_method .ne. 0 .and. iout7 .gt. 0 .and. remd_method .ne. -1) then
    if (remd_method .eq. 4) then
      write (mdout, 9066) solvph, repnum, mdloop
    else if (remd_method .eq. 5) then
      write (mdout, 9071) solve, repnum, mdloop
    else
      if (isgld > 0)then
        write(mdout,9064) temp0,sgfti,sgffi, replica_indexes(1),repnum,mdloop
      else
      write (mdout, 9065) temp0, repnum, mdloop
      endif
    end if
  else if (remd_method .eq. -1 .and. iout7 .gt. 0) then
    write (mdout, 967) remd_dimension, mdloop
    write (mdout, 968) (group_num(i), i=1, remd_dimension)
    write (mdout, 969) (replica_indexes(i), i=1, remd_dimension)
  end if
#endif

!Added by DG for NEB energy output.
#ifdef MPI
  if (ineb>0) call neb_energy_report(mdout)
#endif

!  wxw: EMAP energy
   if (iemap>0 ) then
      write (mdout,9078) si(si_emap_ene)
   endif

  write(mdout, 8088)

! Check if we need to force a flush of mdout. Barring changes in the unix
! system call, the clock is going to wrap in 2038, and flushing won't be
! strictly correct for a flush_interval...

! If we flush the mdout, we also flush the cpout and/or ceout here if we're doing CpHMD

  call get_wall_time(current_sec, current_usec)

  if (first_6call) then
    first_6call = .false.
    next_6flush_sec = current_sec + mdout_flush_interval
    close(mdout)
    open(unit=mdout, file=mdout_name, status='OLD', position='APPEND')
    if (icnstph .gt. 0 .and. .not. cpein_specified) then
      close(cpout)
      open(unit=cpout, file=cpout_name, status='OLD', position='APPEND')
    end if
    if (icnste .gt. 0 .and. .not. cpein_specified) then
      close(ceout)
      open(unit=ceout, file=ceout_name, status='OLD', position='APPEND')
    end if
    if ((icnstph .gt. 0 .or. icnste .gt. 0) .and. cpein_specified) then
      close(cpeout)
      open(unit=cpeout, file=cpeout_name, status='OLD', position='APPEND')
    end if
! PHMD
    if (iphmd .gt. 0) then
      close(phmdout_unit)
      open(unit=phmdout_unit, file=phmdout_name, status='OLD',position='APPEND')
    end if
#ifdef TIDECOMP
    if (ntwd .gt. 0) then
      close(decomp_unit)
      open(unit=decomp_unit, file=decomp_name, status='OLD',position='APPEND')
    end if
#endif
  else
    if (current_sec .ge. next_6flush_sec) then
      next_6flush_sec = current_sec + mdout_flush_interval
      close(mdout)
      open(unit=mdout, file=mdout_name, status='OLD', position='APPEND')
      if (icnstph .gt. 0 .and. .not. cpein_specified) then
        close(cpout)
        open(unit=cpout, file=cpout_name, status='OLD', position='APPEND')
      end if
      if (icnste .gt. 0 .and. .not. cpein_specified) then
        close(ceout)
        open(unit=ceout, file=ceout_name, status='OLD', position='APPEND')
      end if
      if ((icnstph .gt. 0 .or. icnste .gt. 0) .and. cpein_specified) then
        close(cpeout)
        open(unit=cpeout, file=cpeout_name, status='OLD', position='APPEND')
      end if
! PHMD
      if (iphmd .gt. 0) then
        close(phmdout_unit)
        open(unit=phmdout_unit, file=phmdout_name, status='OLD',position='APPEND')
      end if
    end if
  end if

  if (iout7 .eq. 0) return

! Flush i/o buffer:

! Flushing actually does not work particularly reliably for a number of
! machines and compilers, and in the more benign cases simply fails, but in
! the more malign cases can actually corrupt the stack (due to a compiler-
! dependent flush() call interface change).  We therefore no longer do
! flushes of anything in PMEMD; if it needs to go out, we close it and reopen
! it.

! Output the mdinfo file if requested, and if the flush interval has elapsed.

  if (first_7call) then
    first_7call = .false.
    next_7flush_sec = current_sec + mdinfo_flush_interval
  else
    if (current_sec .lt. next_7flush_sec) return
    next_7flush_sec = current_sec + mdinfo_flush_interval
  end if

  call amopen(mdinfo, mdinfo_name, 'U', 'F', 'W')

  write(mdinfo, 9021) nstep, time, temp, si(si_tot_press)
  write(mdinfo, 9029) si(si_tot_ene), si(si_kin_ene), si(si_pot_ene)
  write(mdinfo, 9039) si(si_bond_ene), si(si_angle_ene), si(si_dihedral_ene)
  if (charmm_active) then
    write(mdinfo, 9040) si(si_angle_ub_ene), si(si_dihedral_imp_ene), &
                        si(si_cmap_ene)
  ! if CMAP is active
  else if (cmap_term_count > 0) then
    write(mdinfo, 9040) 0.d0,0.d0,si(si_cmap_ene)
  endif
  write(mdinfo, 9049) si(si_vdw_14_ene), si(si_elect_14_ene), si(si_vdw_ene)

  if (using_pme_potential) then
    write(mdinfo, 9061) si(si_elect_ene), si(si_hbond_ene), si(si_restraint_ene)
  else
    write(mdinfo, 9062) si(si_elect_ene), si(si_hbond_ene), si(si_restraint_ene)
  end if


  if(efx .ne. 0 .or. efy .ne. 0 .or. efz .ne. 0) then
    write(mdinfo, 9091) si(si_efield_ene)
  end if
#ifndef MPI
#ifndef NOXRAY
  if (xray_active) call xray_write_md_state(7)
#endif
#endif
  if (si(si_restraint_ene) .ne. 0.0) &
    write(mdinfo, 9079) si(si_pot_ene) - si(si_restraint_ene)

  if (iamd.gt.0) &
    write(mdinfo, 9080) si(si_amd_boost)

  if (rms .and. ntt .eq. 0) write(mdinfo,9075) si(si_solvent_kin_ene)

  if (ntp .gt. 0 .or. si(si_volume) .ne. 0.0) &
    write(mdinfo, 9081) si(si_tot_ekcmt), si(si_tot_virial), si(si_volume)

  if (csurften .gt. 0) &
    write(mdinfo, 9082) si(si_gamma_ten)

  if (epol .ne. 0.0)  write(mdinfo, 9070) epol

  if (ntp .gt. 0 .or. si(si_volume) .ne. 0.0) write(mdinfo, 9083) si(si_density)

  if (gbsa .eq. 1 .or. gbsa .eq. 3) &
    write(mdinfo, 9063) si(si_surf_ene)

  if (si(si_dvdl) .ne. 0.d0) write(mdinfo, 9089) si(si_dvdl)

  if (igamd.gt.0) then
    if(.not.((igamd.ge.12.and.igamd.le.27).or.&
            (igamd.ge.110 .and. igamd.le.120)))then
        write(mdout, 9090) si(si_gamd_boost)
    else
        write(mdout, 9092) si(si_gamd_boost),si(si_gamd_ppi)
    endif
  endif

#ifndef CUDA
  !Ewald error estimate is not calculated when running on GPUs.
  !Do not print error estimate when running IPS.
  if (using_pme_potential .and. ips == 0) write(mdinfo, 9188) si(si_pme_err_est)
#endif

!  Printout SGLD guiding information
if(.not. usemidpoint) then
   if (isgld > 0) call sgld_print(mdinfo)
endif

#ifdef MPI
  if (remd_method .ne. 0 .and. iout7 .gt. 0 .and. remd_method .ne. -1) then
    if (remd_method .eq. 4) then
      write (mdinfo, 9066) solvph, repnum, mdloop
    else if (remd_method .eq. 5) then
      write (mdinfo, 9071) solve, repnum, mdloop
    else
      if (isgld > 0)then
        write(mdinfo,9064) temp0,sgfti,sgffi, replica_indexes(1),repnum,mdloop
      else
      write (mdinfo, 9065) temp0, repnum, mdloop
      endif
    end if
  else if (remd_method .eq. -1 .and. iout7 .gt. 0) then
    write (mdinfo, 967) remd_dimension, mdloop
    write (mdinfo, 968) (group_num(i), i=1, remd_dimension)
    write (mdinfo, 969) (replica_indexes(i), i=1, remd_dimension)
  end if
#endif

if(.not. usemidpoint) then
  if (nmropt .ne. 0) call nmrptx(mdinfo)
endif
!  wxw: EMAP energy
   if (iemap>0 ) then
      write (mdinfo,9078) si(si_emap_ene)
   endif

  !Print Timing estimates to mdinfo.
  if (nstep /= 0) call print_ongoing_time_summary(total_nstlim,nstep,dt,mdinfo)

  close(mdinfo)

  return

8088 format(t2, 78('-'), /)
9021 format(/1x, 'NSTEP =', i9, 3x, 'TIME(PS) =', f12.3,2X, &
            'TEMP(K) =', f9.2, 2x, 'PRESS =', f8.1)
!! two more digits for GTI on Windows
            !!TBL
#if ( defined(_WIN32) && defined(_DEBUG) )
9029 format(1x, 'Etot   = ', f14.6, 2x, 'EKtot   = ', f14.6, 2x, &
            'EPtot      = ', f14.6)
9039 format(1x, 'BOND   = ', f14.6, 2x, 'ANGLE   = ', f14.6, 2x, &
            'DIHED      = ', f14.6)
9040 format(1x, 'UB     = ', f14.6, 2x, 'IMP     = ', f14.6, 2x, &
            'CMAP       = ', f14.6)
9049 format(1x, '1-4 NB = ', f14.6, 2x, '1-4 EEL = ', f14.6, 2x, &
            'VDWAALS    = ', f14.6)
9061 format(1x, 'EELEC  = ', f14.6, 2x, 'EHBOND  = ', f14.6, 2x, &
            'RESTRAINT  = ', f14.6)
9062 format(1x, 'EELEC  = ', f14.6, 2x, 'EGB     = ', f14.6, 2x, &
            'RESTRAINT  = ', f14.6)
9063 format(1x, 'ESURF= ',f14.6)
#else
9029 format(1x, 'Etot   = ', f14.4, 2x, 'EKtot   = ', f14.4, 2x, &
            'EPtot      = ', f14.4)
9039 format(1x, 'BOND   = ', f14.4, 2x, 'ANGLE   = ', f14.4, 2x, &
            'DIHED      = ', f14.4)
9040 format(1x, 'UB     = ', f14.4, 2x, 'IMP     = ', f14.4, 2x, &
            'CMAP       = ', f14.4)
9049 format(1x, '1-4 NB = ', f14.4, 2x, '1-4 EEL = ', f14.4, 2x, &
            'VDWAALS    = ', f14.4)
9061 format(1x, 'EELEC  = ', f14.4, 2x, 'EHBOND  = ', f14.4, 2x, &
            'RESTRAINT  = ', f14.4)
9062 format(1x, 'EELEC  = ', f14.4, 2x, 'EGB     = ', f14.4, 2x, &
            'RESTRAINT  = ', f14.4)
9063 format(1x, 'ESURF= ',f14.4)
#endif
! Xiongwu: add sgft, sgff for RXSGLD
9064 format (1x,'TEMP= ',f6.1,1x,'SGFT= ',f7.4, 1X,'SGFF= ',f7.4, 1X, &
        'STAGE= ',i4,1x,'REPNUM= ',i4,1x,'EXCH= ',i6)
9065 format(1x,'TEMP0  = ',f14.4,2x,'REPNUM  = ',i14,2x,'EXCHANGE#  = ',i14)
9066 format(1x,'SOLVPH = ',f14.4,2x,'REPNUM  = ',i14,2x,'EXCHANGE#  = ',i14)
9071 format(1x,'SOLVE  = ',f14.6,2x,'REPNUM  = ',i14,2x,'EXCHANGE#  = ',i14)
967  format(1x,'REMD_DIMENSION = ',i14,2x,'EXCHANGE#  = ',i14)
968  format(1x,'MULTI-D REMD. GROUP INDEXES:   ',12i5)
969  format(1x,'MULTI-D REMD. REPLICA INDEXES: ', 12i5)
9070 format(1x, 'EPOLZ  = ', f14.4)
9075 format('|E(PBS) = ',f14.4)
9078 format (1x,'EMAP   = ',f14.4)
9079 format(1x, 'EAMBER (non-restraint)  = ', f14.4)
9080 format(1x, 'EAMD_BOOST  = ', f14.4)
9081 format(1x, 'EKCMT  = ', f14.4, 2x, 'VIRIAL  = ', f14.4, 2x, &
            'VOLUME     = ', f14.4)
9082 format(52x, 'SURFTEN    = ', f14.4)
9083 format(52x, 'Density    = ', f14.4)
#if ( defined(_WIN32) && defined(_DEBUG) )  
9089 format(1x, 'DV/DL  = ',f14.6)
#else
9089 format(1x, 'DV/DL  = ',f14.4)
#endif
9188 format(1x, 'Ewald error estimate: ', e12.4)
9090 format(1x, 'GAMD   = ', f14.4)
9091 format(1x, 'ELECTRIC_FIELD = ', f14.4)
9092 format(1x, 'GAMD   = ', f14.4, 1X,'GAMD_environment',f14.4)

end subroutine prntmd

end module runfiles_mod
