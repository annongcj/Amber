#include "copyright.i"

!*******************************************************************************
!
! Module:  constraints_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module constraints_mod

  implicit none

! Global data definitions.

  ! The following storage is per-process common; ie., it SHOULD be
  ! broadcast from the master to the other processes!

  integer, parameter    :: constraints_dat_int_cnt = 2

  integer                        natc, belly_atm_cnt

  common / constraints_dat_int / natc, belly_atm_cnt

  save  :: / constraints_dat_int /

  ! atm_igroup = belly group atom selection array
  ! atm_jrc = restraint group atom selection array
  ! atm_weight = atom weights for position constraints array.
  ! atm_xc = atom position coordinates for constraints array.

  integer,              allocatable, save       :: atm_igroup(:)

  ! For weights and constraints:

  integer,              allocatable, save       :: atm_jrc(:)
  double precision,     allocatable, save       :: atm_weight(:)
  double precision,     allocatable, save       :: atm_xc(:,:)

! Hide internal routines:

  private       alloc_constraints_mem

contains

!*******************************************************************************
!
! Subroutine:  init_constraints_dat
!
! Description: <TBS>
!
!*******************************************************************************

subroutine init_constraints_dat(atm_cnt, ibelly, ntr, num_ints, num_reals)

  implicit none

! Formal arguments:

  integer, intent(in)           :: atm_cnt
  integer, intent(in)           :: ibelly
  integer, intent(in)           :: ntr

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals     

  call alloc_constraints_mem(atm_cnt, ibelly, ntr, num_ints, num_reals)

  return

end subroutine init_constraints_dat

!*******************************************************************************
!
! Subroutine:  alloc_constraints_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_constraints_mem(atm_cnt, ibelly, ntr, num_ints, num_reals)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: atm_cnt
  integer, intent(in)           :: ibelly
  integer, intent(in)           :: ntr

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals     

! Local variables:

  integer                       :: alloc_failed

  ! If belly dynamics are used, allocate memory for the belly atoms list:

  if (ibelly .gt. 0) then

    if(allocated(atm_igroup)) deallocate(atm_igroup)
    allocate(atm_igroup(atm_cnt), stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    num_ints = num_ints + size(atm_igroup)

    atm_igroup(:) = 0
  end if

  ! If position restraints are used, allocate memory for them:

  !!if (ntr .gt. 0) then

    allocate(atm_jrc(atm_cnt), &
             atm_weight(atm_cnt), &
             atm_xc(3, atm_cnt), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    num_reals = num_reals + size(atm_weight) + size(atm_xc)

    num_ints = num_ints + size(atm_jrc)

    atm_jrc(:) = 0
    atm_weight(:) = 0.d0
    atm_xc(:,:) = 0.d0

  !!end if

  return

end subroutine alloc_constraints_mem

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  bcast_constraints_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine bcast_constraints_dat(atm_cnt, ibelly, ntr)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: atm_cnt
  integer, intent(in)           :: ibelly
  integer, intent(in)           :: ntr

! Local variables:

  integer               :: num_ints, num_reals  ! returned values discarded

  call mpi_bcast(natc, constraints_dat_int_cnt, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  if (.not. master) then
    num_ints = 0
    num_reals = 0
    call alloc_constraints_mem(atm_cnt, ibelly, ntr, num_ints, num_reals)
  end if

  if (ibelly .gt. 0) then
    call mpi_bcast(atm_igroup, atm_cnt, mpi_integer, 0, &
                   pmemd_comm, err_code_mpi)
  end if

  ! If position restraints are used, broadcast them:

  if (ntr .gt. 0) then
    call mpi_bcast(atm_jrc, atm_cnt, mpi_integer, 0, &
                   pmemd_comm, err_code_mpi)

    call mpi_bcast(atm_weight, atm_cnt, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)

    call mpi_bcast(atm_xc, 3 * atm_cnt, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)

  end if

  return

end subroutine bcast_constraints_dat
#endif

!*******************************************************************************
!
! Subroutine:  get_crd_constraint_energy
!
! Description: Routine to put harmonic constraints for position.
!
! Mods for Rev A by GLS.
!
!*******************************************************************************

subroutine get_crd_constraint_energy(natc, econ, jrc, x, frc, xc, weit)

  use parallel_dat_mod
  use ti_mod

  implicit none

! Formal arguments:

  integer               :: natc
  double precision      :: econ
  integer               :: jrc(*)
  double precision      :: x(3, *)
  double precision      :: frc(3, *)
  double precision      :: xc(3, *)
  double precision      :: weit(*)

! Local variables:

  double precision      :: ax, ay, az
  double precision      :: eadd
  integer               :: i, j
  double precision      :: wt
  double precision      :: wx, wy, wz

  econ = 0.0d+00

! BUGBUG - A more efficient implementation would modify jrc to contain
!          only atoms owned by this processor, but this constraint stuff
!          is probably not that much used...

  do j = 1, natc
    i = jrc(j)
#ifdef MPI
    if (gbl_atm_owner_map(i) .eq. mytaskid) then
#endif
      wt = weit(j)
      ax = x(1, i) - xc(1, i)
      ay = x(2, i) - xc(2, i)
      az = x(3, i) - xc(3, i)
      wx = wt * ax
      wy = wt * ay
      wz = wt * az
      eadd = wx * ax + wy * ay + wz * az

      if (ti_mode .ne. 0) then
        if (ti_mode .ne. 1) then
          ti_lscale = ti_sc_lst(i) .lt. 2
        end if
        if (ti_mode .ne. 1 .and. .not. ti_lscale) then
          if (ti_lst(1,i) .ne. 0) then
            ti_ene_aug(1,ti_rest_dist_ene) = &
              ti_ene_aug(1,ti_rest_dist_ene) + eadd 
            eadd = 0.d0
          else if(ti_lst(2,i) .ne. 0) then
            ti_ene_aug(2,ti_rest_dist_ene) = &
              ti_ene_aug(2,ti_rest_dist_ene) + eadd 
            eadd = 0.d0
          end if
        else if (ti_mode .eq. 1 .or. ti_lscale) then
          ti_region = 0
          if (ti_lst(1,i) .ne. 0) then
            ti_region = 1
          else if(ti_lst(2,i) .ne. 0) then
            ti_region = 2
          end if
          if (ti_region .gt. 0) then
            call ti_update_ene(eadd, si_restraint_ene, ti_region,2)
            eadd = eadd * ti_weights(ti_region)       
            wx = wx * ti_weights(ti_region)
            wy = wy * ti_weights(ti_region)
            wz = wz * ti_weights(ti_region)
          end if
        end if
      end if

      econ = econ + eadd
      frc(1, i) = frc(1, i) - (wx + wx)
      frc(2, i) = frc(2, i) - (wy + wy)
      frc(3, i) = frc(3, i) - (wz + wz)
#ifdef MPI
    end if
#endif
  end do

  return

end subroutine get_crd_constraint_energy

!*******************************************************************************
!
! Subroutine:  bellyf
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine bellyf(atm_cnt, igrp, frc)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  integer               :: igrp(*)
  double precision      :: frc(3, *)

! Local variables:

#ifdef MPI
  integer               :: atm_lst_idx
#endif
  integer               :: i

#ifdef MPI
  do atm_lst_idx = 1, my_atm_cnt
    i = gbl_my_atm_lst(atm_lst_idx)
#else
  do i = 1, atm_cnt
#endif
    if (igrp(i) .gt. 0) cycle
    frc(:, i) = 0.d0
  end do

  return

end subroutine bellyf

!*******************************************************************************
!
! Subroutine:  all_atom_belly
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine all_atom_belly(atm_cnt, igroup, vec3d)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  integer               :: igroup(*)
  double precision      :: vec3d(*)

! Local variables:

  integer                       :: i
  integer                       :: i3
  double precision, parameter   :: zero = 0.d0

  do i = 1, atm_cnt
    if (igroup(i) .le. 0) then
      i3 = 3 * i - 3
      vec3d(i3 + 1) = zero
      vec3d(i3 + 2) = zero
      vec3d(i3 + 3) = zero
    end if
  end do

  return

end subroutine all_atom_belly

!*******************************************************************************
!
! Subroutine:   read_restraints
!
! Description:  Routine to read the reference positions for restraining.
!               DRR - Obsolete 
!*******************************************************************************

!subroutine read_restraints(natom, ntrx, xc)
!
!  use file_io_mod
!  use file_io_dat_mod
!  use parallel_dat_mod
!
!  implicit none
!
!! Formal arguments:
!
!  integer               :: natom
!  integer               :: ntrx
!  double precision      :: xc(*)
!
!! Local variables:
!
!  ! DRR - Open and close in pvt_read_restraints
!  !if (ntrx .le. 0) then
!  !  call amopen(refc, refc_name, 'O', 'U', 'R')
!  !else if (ntrx .eq. 1) then
!  !  call amopen(refc, refc_name, 'O', 'F', 'R')
!  !end if
!
!  call pvt_read_restraints(natom, ntrx, xc)
!
!  !if (ntrx .ne. 2) close(refc)
!
!  return
!
!end subroutine read_restraints

!*******************************************************************************
!
! Subroutine:   pvt_read_restraints
!
! Description:  Routine to read the reference positions for restraining.
!               DRR - Making this read_restraints - no need for the extra 
!               function call which separated reading from file opening.
!              
!*******************************************************************************

subroutine read_restraints(natom, ntrx, xc)

  use file_io_mod
  use file_io_dat_mod
  use pmemd_lib_mod
  use AmberNetcdf_mod, only: NC_checkRestart
  use binrestart_mod, only: read_nc_refc

  implicit none

! Formal arguments:

  integer               :: natom
  integer               :: ntrx
  double precision      :: xc(*)


! Local variables:

  character(80)         :: title
  logical               :: formatted_input
  integer               :: i
  integer               :: nr3
  integer               :: refc_natom
  character             :: read_buf(80)            ! For format checking
  integer               :: inpcrd_version          ! For pmemd 2.01 format
  integer               :: zero_natom              ! For pmemd 2.01 format

  nr3 = 3 * natom
  write(mdout, 9100)

  formatted_input =  (ntrx .ne. 0)

  ! ---=== NetCDF Restart File
  if (NC_checkRestart(refc_name)) then
    call read_nc_refc(refc_name,title,natom, xc)
  
  ! ---=== Formatted (ASCII) Restart File
  else if (ntrx .eq. 1) then
    call amopen(refc, refc_name, 'O', 'F', 'R')

    ! Sander 7 uses a mechanism of checking if the 6th char is blank to see if
    ! this is an old i5 format sander 6 file or a new i6 format sander 7 file.
    ! We also check for the versioned pmemd 2.01 format, with an initial natom
    ! of 0.  This was really the best approach, but unfortunately sander 7
    ! did not implement something like it.

    read(refc, 9008) title
    read(refc, '(80a1)') read_buf
    rewind(refc)
    read(refc, 9008) title

    if (read_buf(6) .eq. ' ') then

      ! It is either sander 6 or pmemd 2.01 format...

       read(refc, '(i5)') refc_natom
       rewind(refc)
       read(refc, 9008) title

       if (refc_natom .ne. 0) then
         read(refc, 9018) refc_natom
       else
         read(refc, '(2i5, i10)') &
              zero_natom, inpcrd_version, refc_natom
         if (inpcrd_version .ne. 2) then
           write(mdout, 9128)
           call mexit(6, 1)
         end if
       end if
    else if (read_buf(7) .eq. ' ') then ! Sander 7/8/9/10 large system format...
      read(refc, 9019) refc_natom
    else if (read_buf(8) .eq. ' ') then ! Sander 11 - 1 mil+ format
      read(refc, 9020) refc_natom
    else        ! assume amber 11 VERY large system format. 10 mil+
      read(refc, 9021) refc_natom
    end if

    if (refc_natom .ne. natom) then
      write(mdout, 9118)
      call mexit(6, 1)
    end if

    read(refc, 9028, end = 1000, err = 1000) xc(1 : nr3)
    close(refc)

  ! ---=== Unformatted (Binary) Restart File
  else 
    call amopen(refc, refc_name, 'O', 'U', 'R')

    read(refc) title
    read(refc) refc_natom
    if (refc_natom .ne. natom) then
      write(mdout, 9118)
      call mexit(6, 1)
    end if

    read(refc, end = 1000, err = 1000) xc(1 : nr3)
    close(refc)
  endif

  write(mdout, 9009) title

  return

1000 continue

  write(mdout, '(a,a)') 'FATAL: Could not read constraint coords from ', &
                        refc_name
  call mexit(6, 1)

9008 format(a80)
9009 format(2x, a80)
9018 format(i5)
9019 format(i6)
9020 format(i7)
9021 format(i8)
9028 format(6f12.7)
9100 format(/, '   5.  REFERENCE ATOM COORDINATES', /)
9118 format(/2x, 'FATAL: NATOM mismatch in constraint coord and prmtop files')
9128 format(/2x, 'FATAL: Unrecognized format for constraint coord file')

end subroutine read_restraints

!*******************************************************************************
!
! Subroutine:   remove_nonbelly_bnd_ang_dihed
!
! Description:  Routine to do the necessary accomodations for protein belly
!               minimisations.  Only call if ibelly .gt. 0!
!*******************************************************************************

subroutine remove_nonbelly_bnd_ang_dihed

  use prmtop_dat_mod

  implicit none

! Local variables:

  integer           i

! Delete bonds which are in the belly alone: 

  call remove_nonbelly_bonds(nbonh, gbl_bond, atm_igroup)
  call remove_nonbelly_bonds(nbona, gbl_bond(bonda_idx), atm_igroup)

! Delete the angles which are in the belly alone:

  call remove_nonbelly_angles(ntheth, gbl_angle, atm_igroup)
  call remove_nonbelly_angles(ntheta, gbl_angle(anglea_idx), atm_igroup)

! Delete the dihedrals:

  call remove_nonbelly_dihedrals(nphih, gbl_dihed, atm_igroup)
  call remove_nonbelly_dihedrals(nphia, gbl_dihed(diheda_idx), atm_igroup)

  return

end subroutine remove_nonbelly_bnd_ang_dihed

!*******************************************************************************
!
! Subroutine:  remove_nonbelly_bonds
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine remove_nonbelly_bonds(nb, bond, igrp)

  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer               :: nb
  type(bond_rec)        :: bond(*)
  integer               :: igrp(*)

! Local variables:

  integer               :: i
  integer               :: nba

  nba = 0

  do i = 1, nb
    if (igrp(bond(i)%atm_i) .gt. 0 .or. igrp(bond(i)%atm_j) .gt. 0) then
      nba = nba + 1
      bond(nba) = bond(i)
    end if
  end do

  nb = nba

  return

end subroutine remove_nonbelly_bonds

!*******************************************************************************
!
! Subroutine:  remove_nonbelly_angles
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine remove_nonbelly_angles(nt, angle, igrp)

  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer               :: nt
  type(angle_rec)       :: angle(*)
  integer               :: igrp(*)

! Local variables:

  integer               :: i
  integer               :: nta

  nta = 0

  do i = 1, nt

    if (igrp(angle(i)%atm_i) .gt. 0 .or. &
        igrp(angle(i)%atm_j) .gt. 0 .or. &
        igrp(angle(i)%atm_k) .gt. 0) then
      nta = nta + 1
      angle(nta) = angle(i)
    end if

  end do

  nt = nta

  return

end subroutine remove_nonbelly_angles

!*******************************************************************************
!
! Subroutine:  remove_nonbelly_dihedrals
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine remove_nonbelly_dihedrals(np, dihed, igrp)

  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer               :: np
  type(dihed_rec)       :: dihed(*)
  integer               :: igrp(*)

! Local variables:

  integer               :: i
  integer               :: npa

  npa = 0

  do i = 1, np

    if (igrp(dihed(i)%atm_i) .gt. 0 .or. &
        igrp(dihed(i)%atm_j) .gt. 0 .or. &
        igrp(iabs(dihed(i)%atm_k)) .gt. 0 .or. &
        igrp(iabs(dihed(i)%atm_l)) .gt. 0) then
      npa = npa + 1
      dihed(npa) = dihed(i)
    end if

  end do

  np = npa

  return

end subroutine remove_nonbelly_dihedrals

!*******************************************************************************
!
! Subroutine:  count_belly_atoms
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine count_belly_atoms(atm_cnt, belly_atm_cnt, igrp)

  implicit none

  integer       atm_cnt
  integer       belly_atm_cnt
  integer       igrp(*)

  integer       i

  belly_atm_cnt = 0

  do i = 1, atm_cnt
    if (igrp(i) .gt. 0) belly_atm_cnt = belly_atm_cnt + 1
  end do

  return

end subroutine count_belly_atoms

end module constraints_mod
