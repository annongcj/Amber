#include "copyright.i"

!*******************************************************************************
!
! Module: hybridsolvent_remd_mod
!
! Description: 
!
! Module for holding the necessary subroutines required for Hybrid Solvent REMD.
! The sander implementation was done by Daniel Roe, based on the original
! implementation by Guanglei Cui and Carlos Simmerling. Adapted to PMEMD by Koushik K.
!
! The main exchange subroutines for this REMD calculation are still
! in remd.F90 and remd_exchg.F90. Here, only the non-exchange related subroutines
! are included, which are listed below.
!
! Hybrid Solvent REMD subroutines:
! hybridsolvent_remd_initial_setup - Populate the necessary variables
! hybridsolvent_remd_bcast - Broadcast some variables
! hybridsolvent_remd_strip_water - Strips the molecule to keep desired number of atoms
! open_hybridsolvent_remd_traj_binary_file - Open file for writing stripped coordinates
! hybridsolvent_remd_traj_setup_remd_indices - Setup indices for each replica
! write_hybridsolvent_remd_traj_binary_crds - Write the stripped coords to the file.
!
!*******************************************************************************

module hybridsolvent_remd_mod

#ifdef MPI
  ! For writing stripped coordinates in Netcdf format
  integer, save                 :: hybridsolvent_remd_traj_ncid
  integer, save                 :: hybridsolvent_remd_traj_coord_var_id
  integer, save                 :: hybridsolvent_remd_traj_time_var_id
  integer, save                 :: hybridsolvent_remd_traj_frame
  integer, save                 :: hybridsolvent_remd_traj_repidx_var_id
  integer, save                 :: hybridsolvent_remd_traj_crdidx_var_id
  integer, save                 :: hybridsolvent_remd_traj_remd_indices_var_id
  integer, save                 :: hybridsolvent_remd_traj_temp_var_id
  integer, save                 :: hybridsolvent_remd_traj_frame_dim_id
  integer, save                 :: hybridsolvent_remd_traj_remd_values_var_id
  integer, save                 :: hybridsolvent_remd_traj_remd_types_var_id
  integer, public, save         :: hybridsolvent_remd_atm_cnt

public hybridsolvent_remd_initial_setup, hybridsolvent_remd_bcast, &
       hybridsolvent_remd_strip_water, &
       open_hybridsolvent_remd_traj_binary_file, &
       hybridsolvent_remd_traj_setup_remd_indices, &
       write_hybridsolvent_remd_traj_binary_crds

contains

#ifndef BINTRAJ
subroutine no_nc_error
  use AmberNetcdf_mod, only : NC_NoNetcdfError
  use file_io_dat_mod, only : mdout
  use pmemd_lib_mod,   only : mexit
  call NC_NoNetcdfError(mdout)
  call mexit(mdout, 1)
end subroutine no_nc_error
#endif


!*******************************************************************************
!
! Subroutine: hybridsolvent_remd_initial_setup
!
! Description: Initiatizes Hybrid Solvent REMD data
!
!*******************************************************************************

subroutine hybridsolvent_remd_initial_setup(num_ints, num_reals)

  use mdin_ctrl_dat_mod
  use prmtop_dat_mod
  use pmemd_lib_mod,        only: get_atomic_number

  implicit none

! Passed variables

  integer, intent(in out) :: num_ints
  integer, intent(in out) :: num_reals

! Local variables

  integer       :: atomicnumber ! holder for atomic numbers
  integer       :: i,j,k            ! counter
  integer       :: first_water
  integer       :: num_ions
  integer       :: total_solute_res
  integer       :: total_solute_atm_cnt
  integer       :: total_water_atm_cnt

  if (hybridgb .gt. 0) then ! Should be true. Otherwise, what are we doing here?
    gb_cutoff = 9999.d0  ! GB cutoff is infinite here
    if (hybridgb .eq. 1) then
      gb_alpha = 1.d0
      gb_beta = 0.d0
      gb_gamma = 0.d0
    else if (hybridgb .eq. 2) then
      ! Use our best guesses for Onufriev/Case GB  (GB^OBC I):
      gb_alpha = 0.8d0
      gb_beta = 0.d0
      gb_gamma = 2.909125d0
    else if (hybridgb .eq. 5) then
      ! Use our second best guesses for Onufriev/Case GB (GB^OBC II):
      gb_alpha = 1.d0
      gb_beta = 0.8d0
      gb_gamma = 4.85d0
    else if (hybridgb .eq. 7) then
      ! Use parameters for Mongan et al. CFA GBNECK:
      gb_alpha = 1.09511284d0
      gb_beta = 1.90792938d0
      gb_gamma = 2.50798245d0
      gb_neckscale = 0.361825d0
    else if (hybridgb .eq. 8) then
      gb_alpha_h   = 0.788440d0
      gb_beta_h    = 0.798699d0
      gb_gamma_h   = 0.437334d0
      gb_alpha_c   = 0.733756d0
      gb_beta_c    = 0.506378d0
      gb_gamma_c   = 0.205844d0
      gb_alpha_n   = 0.503364d0
      gb_beta_n    = 0.316828d0
      gb_gamma_n   = 0.192915d0
      gb_alpha_os  = 0.867814d0
      gb_beta_os   = 0.876635d0
      gb_gamma_os  = 0.387882d0
      gb_alpha_p   = 0.418365d0  ! dac, 4/24: for now use gb_alpha_pnu values
      gb_beta_p    = 0.290054d0
      gb_gamma_p   = 0.1064245d0
      gb_neckscale = 0.826836d0
      screen_h = 1.425952d0
      screen_c = 1.058554d0
      screen_n = 0.733599d0
      screen_o = 1.061039d0
      if (iphmd .gt. 0) then
        screen_s = 1.061039d0
      else
        screen_s = -0.703469d0
      end if
      screen_p = 0.5d0
      offset  = 0.195141d0

      ! update gbneck2nu pars
      ! using offset and gb_neckscale parameters from GB8-protein for GBNeck2nu
      ! Scaling factors
      screen_hnu = 1.696538d0
      screen_cnu = 1.268902d0
      screen_nnu = 1.4259728d0
      screen_onu = 0.1840098d0
      screen_pnu = 1.5450597d0
      !alpha, beta, gamma for each atome element
      gb_alpha_hnu = 0.537050d0
      gb_beta_hnu = 0.362861d0
      gb_gamma_hnu = 0.116704d0
      gb_alpha_cnu = 0.331670d0
      gb_beta_cnu = 0.196842d0
      gb_gamma_cnu = 0.093422d0
      gb_alpha_nnu = 0.686311d0
      gb_beta_nnu = 0.463189d0
      gb_gamma_nnu = 0.138722d0
      gb_alpha_osnu = 0.606344d0
      gb_beta_osnu = 0.463006d0
      gb_gamma_osnu = 0.142262d0
      gb_alpha_pnu = 0.418365d0
      gb_beta_pnu = 0.290054d0
      gb_gamma_pnu = 0.1064245d0
    
    end if

    ! Set gb_kappa as long as saltcon is .ge. 0.d0 (and catch bug below if not).

    if (saltcon .ge. 0.d0) then

      ! Get Debye-Huckel kappa (A**-1) from salt concentration (M), assuming:
      !   T = 298.15, epsext=78.5,
      gb_kappa = sqrt(0.10806d0 * saltcon)

      ! Scale kappa by 0.73 to account(?) for lack of ion exclusions:

      gb_kappa = 0.73d0 * gb_kappa

    end if

   ! If hybridgb .eq. 7 use special S_x screening params; here we overwrite
   ! the tinker values read from the prmtop.

    if (hybridgb .eq. 7) then

      write(mdout,'(a)') &
        ' Replacing prmtop screening parameters with GBn (igb=7) values'

      do i = 1, natom

        if (loaded_atm_atomicnumber) then
          atomicnumber = atm_atomicnumber(i)
        else
          call get_atomic_number(atm_igraph(i), atm_mass(i), atomicnumber)
        end if

        if (atomicnumber .eq. 6) then
          atm_gb_fs(i) = 4.84353823306d-1
        else if (atomicnumber .eq. 1) then
          atm_gb_fs(i) = 1.09085413633d0
        else if (atomicnumber .eq. 7) then
          atm_gb_fs(i) = 7.00147318409d-1
        else if (atomicnumber .eq. 8) then
          atm_gb_fs(i) = 1.06557401132d0
        else if (atomicnumber .eq. 16) then
          atm_gb_fs(i) = 6.02256336067d-1
        else
          atm_gb_fs(i) = 5.d-1 ! not optimized
        end if

      end do

    else if (hybridgb .eq. 8) then

      write(mdout, '(a)') &
        ' Replacing prmtop screening parameters with GBn2 (igb=8) values'

      do i = 1, natom

        if(loaded_atm_atomicnumber) then
          atomicnumber = atm_atomicnumber(i)
        else
          call get_atomic_number(atm_igraph(i), atm_mass(i), atomicnumber)
        end if

        ! The screen_ variables are found in mdin_ctrl_dat_mod

        call isnucat(nucat,i,nres,60,gbl_res_atms(1:nres),gbl_labres(1:nres))
        ! check if atom belong to nucleic or protein residue
        ! 37 = the total number of different nucleic acid residue type in AMBER
        ! Need to update this number if adding new residue type
                ! update isnucat function in gb_ene.F90 too

        ! RCW - HACKY CODE ALERT!
        !       This (and the isnucat subroutine) is incredibly fragile and people will forget
        !       to update this. Plus one can in principal use any residue name they want. This
        !       is a total HACK and should NOT be here. Make it robust and move it to leap where
        !       it belogs please.

        if (nucat == 1) then
            ! If atom belongs to nuc, use GBNeck2nu pars
            if (atomicnumber .eq. 6) then
              atm_gb_fs(i) = screen_cnu
            else if (atomicnumber .eq. 1) then
              atm_gb_fs(i) = screen_hnu
            else if (atomicnumber .eq. 7) then
              atm_gb_fs(i) = screen_nnu
            else if (atomicnumber .eq. 8) then
              atm_gb_fs(i) = screen_onu
            else if (atomicnumber .eq. 15) then
              atm_gb_fs(i) = screen_pnu
            else
              atm_gb_fs(i) = 5.d-1 ! not optimized
            end if
        else
            !if atom does not belong to nuc, use protein pars
            if (atomicnumber .eq. 6) then
              atm_gb_fs(i) = screen_c
            else if (atomicnumber .eq. 1) then
              atm_gb_fs(i) = screen_h
            else if (atomicnumber .eq. 7) then
              atm_gb_fs(i) = screen_n
            else if (atomicnumber .eq. 8) then
              atm_gb_fs(i) = screen_o
            else if (atomicnumber .eq. 16) then
              atm_gb_fs(i) = screen_s
            else if (atomicnumber .eq. 15) then
              atm_gb_fs(i) = screen_p
            else
              atm_gb_fs(i) = 5.d-1 ! not optimized
            end if
         end if

      end do

    end if

  ! Recompute atm_gb_fs for GBNeck and GBNeck2 models
    if (hybridgb .eq. 7 .or. hybridgb .eq. 8) then
       gb_fs_max = 0.d0

       do i = 1, natom
          atm_gb_fs(i) = atm_gb_fs(i) * (atm_gb_radii(i) - offset)
          gb_fs_max = max(gb_fs_max, atm_gb_fs(i))
       end do
    end if

  ! Find number of atoms in stripped system
    first_water = -1
    num_ions = 0

  ! Find first water residue
    do i = 1, nres
      if (gbl_labres(i) .eq. "WAT") then
         first_water = i
         exit
      end if
    end do

    ! Count number of ions if present  
    do i = 1, nres
      ! Add checks for more ions here if necessary
      if (gbl_labres(i) .eq. "Na+" .or. &
          gbl_labres(i) .eq. "K+" .or. &
          gbl_labres(i) .eq. "Cl-") then
         num_ions = num_ions + 1
      end if
    end do

    ! Count number of solute residues. -1 is necessary because first_water is after solute ends    
    total_solute_res = first_water-num_ions-1
    
    ! Count number of solute atoms
    total_solute_atm_cnt = 0
    do j = 1, total_solute_res
      do k=gbl_res_atms(j), gbl_res_atms(j+1)-1
        total_solute_atm_cnt = total_solute_atm_cnt + 1
      end do
    end do

    ! Count number of water atoms. GB only supports 3 point water models.
    ! So number of water atoms will be 3*numwatkeep.
    total_water_atm_cnt = numwatkeep * 3

    hybridsolvent_remd_atm_cnt = total_solute_atm_cnt + total_water_atm_cnt
  end if

  return

end subroutine hybridsolvent_remd_initial_setup

!*******************************************************************************
!
! Subroutine: hybridsolvent_remd_bcast
!
! Description: Broadcasts all Hybrid Solvent REMD data
!
!*******************************************************************************

subroutine hybridsolvent_remd_bcast

  use mdin_ctrl_dat_mod, only : hybridgb, numwatkeep, igb
  use parallel_dat_mod  ! Provides MPI routines/constants, pmemd_comm and master

  implicit none

  call mpi_bcast(hybridgb, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)

  ! !!!!!!!A RISKY WAY TO DO IT!!!!!!!
  ! Assigning hybridgb to igb so that calc_born_radii() in gb_ene() can be used 
  ! as it is. Otherwise, for every 'if igb' condition in calc_born_radii(), 
  ! 'if hybridgb' conditions will have to be added. Setting igb to hybridgb here 
  ! seems to be a simpler solution. Setting it in here instead of in 
  ! mdin_ctrl_dat.F90 will bypass the checks in validate_mdin_ctrl_dat().
  ! Ideally, it is better to set 'if hybridgb' conditions in calc_born_radii() but
  ! I am not sure how much the speed will be affected with the additional 
  ! if conditions. Since Hybrid Solvent REMD is a very specific method, checking 
  ! if conditions for every step of regular MD may not be ideal.
  ! This is a simpler fix. BE CAREFUL THOUGH!!! 
  igb = hybridgb

  call mpi_bcast(igb, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)
  call mpi_bcast(numwatkeep, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)
  call mpi_bcast(hybridsolvent_remd_atm_cnt, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)

end subroutine hybridsolvent_remd_bcast

!*******************************************************************************
!
! Subroutine: hybridsolvent_remd_strip_water
!
! Description: Strips water from structure for replica exchange with mixed
!              solvent models
!
!*******************************************************************************

! This routine requires that the residues are organized so that we have the solute
! atoms first, followed by ions (currently only K+, Na+, Cl- ions are handled) if any,
! and then water molecules.

! It returns:
!    hybridsolvent_remd_atm_cnt --> New number of atoms (stripped system)
!    hybridsolvent_remd_crd --> New coordinates (solute + closest waters)
! The way it works:

! 1. Find where the water molecules start, puts that residue number into 'firstwat'. 
!    If there is any residue other than water after this point, the calculation stops
!    with an error. K+, Na+, Cl- ions can be present in the system which are going to  be
!    removed from the stripped system so that they are not involved in GB energy calculation.

! 2. Reimage the system so that waters surround the solute. Make sure
!    to image the waters correctly to the minimum distance.
!    (Remember that the hybrid-solvent REMD code uses GB, which is not periodic.)

! 3. For every water, calculate the distance from it's oxygen to EVERY
!    non-water atom in the system, keeping only the shortest one.
!    * After this point, we have a list of all the waters, their 
!      residue numbers and the distance to the closest non-water atom.

! 4. Sort the waters by distance. 
!    * At this point, the routine uses a "selection sort" algorithm, 
!      which scales as N**2, N being the number of waters. If this 
!      becomes a bottleneck, it can be changed later to "mergesort" 
!      or a "binary tree" sort, which are N*logN.

! 5. Set a new logical array "keep_water(nres)", which is a logical mask
!    indicating whether to keep this water or not. The first "numwatkeep"
!    waters are set to .true., all the rest are set to .false.
!    * After this sort we have three arrays ordered by increasing distance
!      to the closest non-solute atom:
!      water_resnum(nres) --> The residue number of the water
!      closedistance(nres) --> The distance to the closest non-water residue
!      keep_water(nres)   --> Logical. ".true." if we want to keep this water.

! 6. Now, sort those 3 arrays in the order by residue number instead

! 7. Get the coordinates for the atoms in the solute plus *only* the
!    waters we want to keep. That means looping through all residues and
!    copying the coordiantes only of the non-water atoms, plus the atoms
!    from water residues marked with "keep_water = .true."

! 8. Copy those coordinates to the main coordinate array, update the
! "hybridsolvent_remd_atm_cnt" variable with the number of atoms in this reduced system.

! 9. Return.
!=====================================================================

subroutine hybridsolvent_remd_strip_water(hybridsolvent_remd_atm_cnt, &
                                      hybridsolvent_remd_crd, numwatkeep)

  use prmtop_dat_mod, only : gbl_labres, gbl_res_atms, ifbox
  use prmtop_dat_mod, only : nres, atm_mass
  use pbc_mod
  use pmemd_lib_mod,  only : mexit
  use parallel_dat_mod
  use file_io_dat_mod

  implicit none

! Passed variables

  integer, intent(inout)             :: hybridsolvent_remd_atm_cnt
  double precision, intent(inout)    :: hybridsolvent_remd_crd(3, hybridsolvent_remd_atm_cnt)
  integer, intent(in)                :: numwatkeep

! Local variables

  integer                            :: i, j, m, n, k
  integer                            :: water_pointer ! water pointer
  integer                            :: total_water ! total water in system
  integer                            :: first_water ! first water in system
  integer                            :: num_ions ! number of ions in system
  integer                            :: total_solute_res ! solute residues in system

! Needed for finding closest water molecules and sort them
  double precision                   :: hybridsolvent_remd_crd_oxygen_x
  double precision                   :: hybridsolvent_remd_crd_oxygen_y
  double precision                   :: hybridsolvent_remd_crd_oxygen_z
  double precision                   :: xij, yij, zij
  double precision                   :: r2

  double precision                   :: closedistance(nres)
  integer                            :: water_resnum(nres)
  logical                            :: keep_water(nres)
  double precision                   :: temp_closedistance
  double precision                   :: temp_water_resnum
  logical                            :: temp_keep_water

! Temporary atom count and temporary coordinate array
  integer                            :: hybridsolvent_remd_new_atm_cnt
  double precision                   :: hybridsolvent_remd_tempcrd(3, hybridsolvent_remd_atm_cnt)

! To re-center the system around solute
  double precision                   :: aamass
  double precision                   :: tmassinv
  double precision                   :: xcm(3)


  ! Find out where water starts

  first_water = -1
  num_ions = 0

  do i = 1, nres
    if (gbl_labres(i) .eq. "WAT") then
       first_water = i
       exit
    end if
  end do

  do i = 1, nres
    ! Add checking for more ions here if necessary
    if (gbl_labres(i) .eq. "Na+" .or. &
        gbl_labres(i) .eq. "K+" .or. &
        gbl_labres(i) .eq. "Cl-") then
       num_ions = num_ions + 1
    end if
  end do

  if (first_water .le. 0) then
    if (master) then
      write(mdout,'(a)') 'There are no water molecules in the system. This is a &
                         hybrid-solvent REMD calculation. This requires water to be present'
    end if
  end if

  ! next we need to reimage the system so that the waters surround the solute
  ! we could use minimum image type calculation to get the
  ! closest distance but we still need the water imaged properly
  ! for the GB calculation (which is not periodic)
  ! follow the code for imaging from runmd's iwrap=1 code


  ! First, center the system on the CM of the solute

  xcm(1) = 0.d0
  xcm(2) = 0.d0
  xcm(3) = 0.d0

  ! Here, tmassinv is only for non-water and non-ions
  total_solute_res = first_water-num_ions-1
  tmassinv = 0.d0
  i = 0
  do j = 1, total_solute_res
    do k=gbl_res_atms(j), gbl_res_atms(j+1)-1
      aamass = atm_mass(k)
      xcm(1) = xcm(1) + hybridsolvent_remd_crd(i+1,k) * aamass
      xcm(2) = xcm(2) + hybridsolvent_remd_crd(i+2,k) * aamass
      xcm(3) = xcm(3) + hybridsolvent_remd_crd(i+3,k) * aamass
      i = i + 3
      tmassinv = tmassinv + aamass
    end do
  end do

  tmassinv = 1.d0/tmassinv

  xcm(1) = xcm(1) * tmassinv
  xcm(2) = xcm(2) * tmassinv
  xcm(3) = xcm(3) * tmassinv

  ! Center all atoms, not just solute

  do i = 1, hybridsolvent_remd_atm_cnt
      hybridsolvent_remd_crd(1,i) = hybridsolvent_remd_crd(1,i) - xcm(1)
      hybridsolvent_remd_crd(2,i) = hybridsolvent_remd_crd(2,i) - xcm(2)
      hybridsolvent_remd_crd(3,i) = hybridsolvent_remd_crd(3,i) - xcm(3)
  end do

  ! Now re-image the box

  call wrap_molecules(hybridsolvent_remd_crd)
  if (ifbox .eq. 2) call wrap_to(hybridsolvent_remd_crd)

  ! Now, start setting closest distance between water and solute.
  ! In fact we do not need the distance, the distance squared is
  ! sufficient because we only want to sort them

  water_pointer=0

  ! Check to make sure it is water

  do i = first_water, nres
    if (gbl_labres(i) .ne. "WAT") then
      if (master) then
        write (mdout,'(a,a)') "solvent molecule is not water: ", gbl_labres(i)
        write (mdout,'(a)') "Stopping water search. Please check your topology file"
      end if
      call mexit(mdout, 1)
    end if

    water_pointer = water_pointer + 1

    ! closedistance(i) is the distance from the water to the closest solute atom.
    ! water_resnum(i) is the corresponding residue number for that water in nres.

    closedistance(water_pointer) = 99999
    water_resnum(water_pointer) = i

    ! For water, just take distance to oxygen (first atom in water)

    m = gbl_res_atms(i)
    hybridsolvent_remd_crd_oxygen_x = hybridsolvent_remd_crd(1,m)
    hybridsolvent_remd_crd_oxygen_y = hybridsolvent_remd_crd(2,m)
    hybridsolvent_remd_crd_oxygen_z = hybridsolvent_remd_crd(3,m)

    ! loop over non-water residues/atoms

    do j = 1, total_solute_res
       do n = gbl_res_atms(j), gbl_res_atms(j+1)-1
          xij = hybridsolvent_remd_crd_oxygen_x - hybridsolvent_remd_crd(1,n)
          yij = hybridsolvent_remd_crd_oxygen_y - hybridsolvent_remd_crd(2,n)
          zij = hybridsolvent_remd_crd_oxygen_z - hybridsolvent_remd_crd(3,n)
          r2  = xij*xij + yij*yij + zij*zij

          if (r2 .lt. closedistance(water_pointer)) then
            closedistance(water_pointer) = r2
          end if
       end do
     end do
  end do

  total_water = water_pointer

  ! now we have a list of all waters - their residue number and distance
  ! to solute

  ! sort them by distance

  do i = 1, total_water-1
    do j = i+1, total_water
      if (closedistance(i) .gt. closedistance(j)) then
        temp_closedistance = closedistance(i)
        temp_water_resnum = water_resnum(i)
        closedistance(i)  = closedistance(j)
        water_resnum(i)  = water_resnum(j)
        closedistance(j)  = temp_closedistance
        water_resnum(j)  = temp_water_resnum
      end if
    end do
  end do

  ! now set save flags for closest numwatkeep

  do i = 1, numwatkeep
     keep_water(i) = .true.
  end do
  do i = numwatkeep+1, total_water
     keep_water(i) = .false.
  end do

  ! now sort them back into the order by water_resnum

  do i = first_water, nres
  ! i is the residue number we are looking for

  ! m is the current water at this residue number
  ! this is in the "1 to total_water" sense so we can pull those
  ! indices for swapping waters in list

  m = i-first_water+1

  ! look to see where this water is, and put it in place

    do j = 1, total_water
      if (water_resnum(j) .eq. i) then
        ! found it, so swap them
        temp_closedistance = closedistance(m)
        temp_water_resnum = water_resnum(m)
        temp_keep_water   = keep_water(m)
        closedistance(m)  = closedistance(j)
        water_resnum(m)  = water_resnum(j)
        keep_water(m)    = keep_water(j)
        closedistance(j)  = temp_closedistance
        water_resnum(j)  = temp_water_resnum
        keep_water(j)    = temp_keep_water
  !      write(mdout, '(a,i3,a,i3,a,i3)') 'i ', i, ' j ', j, ' m ', m
        exit
      end if
    enddo
  enddo

! now go through and write the new coordinates. we need to write to temp
! array because we can't remove atoms while writing the file itself since
! there is more than 1 molecule per line

  hybridsolvent_remd_new_atm_cnt=0
  do i = 1, nres
! i .lt. total_solute_res is for solute and keepwat(i-firstwat+1) is for closest waters
     if (i .le. total_solute_res) then
      !  write(mdout,'(a,i4)') 'ires =', i
        do j = gbl_res_atms(i), gbl_res_atms(i+1)-1 ! Copy hydrogens too
           hybridsolvent_remd_new_atm_cnt = hybridsolvent_remd_new_atm_cnt+1
           hybridsolvent_remd_tempcrd(1,hybridsolvent_remd_new_atm_cnt) = hybridsolvent_remd_crd(1,j)
           hybridsolvent_remd_tempcrd(2,hybridsolvent_remd_new_atm_cnt) = hybridsolvent_remd_crd(2,j)
           hybridsolvent_remd_tempcrd(3,hybridsolvent_remd_new_atm_cnt) = hybridsolvent_remd_crd(3,j)
        enddo
     end if
     if (i .ge. first_water) then
     !   write(mdout,'(a,i4,L6)') 'ires =', i, keep_water(i-first_water+1)
        if (keep_water(i-first_water+1) .eqv. .true.) then
     !      write(mdout,'(a,i4,L6)') 'ires =', i, keepwat(i-firstwat+1)
           do j = gbl_res_atms(i), gbl_res_atms(i+1)-1 ! Copy hydrogens too
              hybridsolvent_remd_new_atm_cnt = hybridsolvent_remd_new_atm_cnt+1
              hybridsolvent_remd_tempcrd(1,hybridsolvent_remd_new_atm_cnt) = hybridsolvent_remd_crd(1,j)
              hybridsolvent_remd_tempcrd(2,hybridsolvent_remd_new_atm_cnt) = hybridsolvent_remd_crd(2,j)
              hybridsolvent_remd_tempcrd(3,hybridsolvent_remd_new_atm_cnt) = hybridsolvent_remd_crd(3,j)
           enddo
        end if
     end if
  enddo

  ! now copy this array to old one, reset hybrid_atm_cnt

  hybridsolvent_remd_atm_cnt = hybridsolvent_remd_new_atm_cnt
  do i = 1, hybridsolvent_remd_atm_cnt
    hybridsolvent_remd_crd(1,i) = hybridsolvent_remd_tempcrd(1,i)
    hybridsolvent_remd_crd(2,i) = hybridsolvent_remd_tempcrd(2,i)
    hybridsolvent_remd_crd(3,i) = hybridsolvent_remd_tempcrd(3,i)
  enddo

  return

end subroutine hybridsolvent_remd_strip_water

!!*******************************************************************************
!!
!! Subroutine: hybridsolvent_remd_writetraj
!!
!! Description: Writes the stripped coordinates for Hybrid Solvent REMD
!!
!!*******************************************************************************
!
!subroutine hybridsolvent_remd_writetraj(hybridsolvent_remd_atm_cnt, &
!                   hybridsolvent_remd_crd, repnum, my_mdloop, local_total_nstep, &
!                   temperature)
!
!  use runfiles_mod,      only : prntmd, corpac
!  use mdin_ctrl_dat_mod, only : ioutfm
!  use file_io_dat_mod
!  use parallel_dat_mod
!  use remd_mod, only : remd_method, remd_dimension, remd_types, group_num
!
!  implicit none
!
!! Passed variables
!
!  integer, intent(in)             :: hybridsolvent_remd_atm_cnt
!  integer, intent(in)             :: repnum
!  integer, intent(in)             :: my_mdloop
!  integer, intent(in)             :: local_total_nstep
!
!  double precision, intent(in)    :: temperature
!  double precision, intent(in)    :: hybridsolvent_remd_crd(3, hybridsolvent_remd_atm_cnt)
!
!  integer                         :: hybridsolvent_remd_crd_cnt
!
!
!  if(master) then
!
!     hybridsolvent_remd_crd_cnt = hybridsolvent_remd_atm_cnt * 3
!   
!     if (ioutfm .eq. 0) then
!        write(hybridsolvent_remd_traj, '(a,3(1x,i8),1x,f8.3)') 'REMD ', repnum, &
!              my_mdloop, local_total_nstep, temp0
!        call corpac(hybridsolvent_remd_crd_cnt, hybridsolvent_remd_crd, 1, &
!                    hybridsolvent_remd_traj)
!     else
!       !! write hybrid traj in netcdf format
!       if(my_mdloop .eq. 0) then
!          call open_hybridsolvent_remd_traj_binary_file(hybridsolvent_remd_atm_cnt)
!          ! Need to pass remd arguments to avoid cyclic dependencies
!          call hybridsolvent_remd_traj_setup_remd_indices(remd_method, remd_dimension, &
!                                                    remd_types, group_num)
!          ! Need to pass remd arguments to avoid cyclic dependencies
!          call write_hybridsolvent_remd_traj_binary_crds(hybridsolvent_remd_crd, &
!                                 hybridsolvent_remd_atm_cnt, remd_method, &
!                                 replica_indexes, remd_dimension, &
!                                 remd_repidx, remd_crdidx)
!       else
!          ! Need to pass remd arguments to avoid cyclic dependencies
!          call write_hybridsolvent_remd_traj_binary_crds(hybridsolvent_remd_crd, &
!                                 hybridsolvent_remd_atm_cnt, &
!                                 remd_method, replica_indexes, remd_dimension, &
!                                 remd_repidx, remd_crdidx)
!       end if
!     end if
!  
!  end if
!
!end subroutine hybridsolvent_remd_writetraj

!*******************************************************************************
!
! Subroutine: open_hybridsolvent_remd_traj_binary_file
!
! Description: Opens files for writing the stripped coordinates in Netcdf format
!              for hybrid-solvent REMD simulations
!
!*******************************************************************************

subroutine open_hybridsolvent_remd_traj_binary_file(hybridsolvent_remd_atm_cnt)
#ifdef BINTRAJ
  use AmberNetcdf_mod, only: NC_create, NC_defineRemdIndices
  use file_io_dat_mod
  use pmemd_lib_mod
  use mdin_ctrl_dat_mod
  use prmtop_dat_mod

  implicit none

  integer, intent(in)       :: hybridsolvent_remd_atm_cnt

  ! Local variables:
  integer                   :: ierr

  hybridsolvent_remd_traj_ncid = -1

  ! Create NetCDF trajectory file for writing stripped coordinates
  if ( NC_create( hybridsolvent_remd_traj_name, owrite, .false., hybridsolvent_remd_atm_cnt, &
                  .true., .false., .false., .true. , .true., .false. , prmtop_ititl, &
                  hybridsolvent_remd_traj_ncid, hybridsolvent_remd_traj_time_var_id, &
                  hybridsolvent_remd_traj_coord_var_id, ierr, ierr, ierr, &
                  ierr, hybridsolvent_remd_traj_temp_var_id, &
                  frameDID=hybridsolvent_remd_traj_frame_dim_id ) ) then
    write(mdout, '(a)') "Creation of Hybrid Solvent REMD NetCDF trajectory file failed."
#   ifndef DUMB
    write (0,*) 'Error on opening ', hybridsolvent_remd_traj_name
#   endif
    call mexit(mdout,1)
  endif

  hybridsolvent_remd_traj_frame = 1

#else
  integer, intent(in)       :: hybridsolvent_remd_atm_cnt
  call no_nc_error
#endif
  return

end subroutine open_hybridsolvent_remd_traj_binary_file

!*******************************************************************************!
! Subroutine:  hybridsolvent_remd_traj_setup_remd_indices
!
! Description:  Set up dimension information for multi-D REMD. This is in a
!               subroutine separate from open_binary_files since REMD setup
!               occurs after the call to open_binary_files.
!
!******************************************************************************
subroutine hybridsolvent_remd_traj_setup_remd_indices(my_remd_method, my_remd_dimension, &
                                    my_remd_types, my_group_num)
#ifdef BINTRAJ
  use AmberNetcdf_mod, only   : NC_defineRemdIndices
  use file_io_dat_mod, only   : mdout
  use pmemd_lib_mod, only     : mexit
# ifdef MPI
!  use remd_mod, only : remd_method, remd_dimension, remd_types, group_num
  use mdin_ctrl_dat_mod, only: ioutfm
# endif
  implicit none
  ! Formal arguments
  integer, intent(in)      :: my_remd_method
  integer, intent(in)      :: my_remd_dimension
  integer, intent(in)      :: my_remd_types(my_remd_dimension)
  integer, intent(in)      :: my_group_num(my_remd_dimension)
# ifdef MPI
  ! Only setup remd indices if ioutfm and multid remd
    if ( NC_defineRemdIndices(hybridsolvent_remd_traj_ncid, my_remd_dimension, &
                            hybridsolvent_remd_traj_remd_indices_var_id, &
                            hybridsolvent_remd_traj_repidx_var_id, &
                            hybridsolvent_remd_traj_crdidx_var_id, &
                            my_remd_types, .false., (my_remd_method .ne. 0), (my_remd_method .eq. -1), &
                            frameDID=hybridsolvent_remd_traj_frame_dim_id, &
                            remd_valuesVID=hybridsolvent_remd_traj_remd_values_var_id, &
                            remd_typesVID=hybridsolvent_remd_traj_remd_types_var_id) ) call mexit(mdout,1)
# endif
#else
  integer, intent(in)      :: my_remd_method
  integer, intent(in)      :: my_remd_dimension
  integer, intent(in)      :: my_remd_types(my_remd_dimension)
  integer, intent(in)      :: my_group_num(my_remd_dimension)
  call no_nc_error
#endif
  return

end subroutine hybridsolvent_remd_traj_setup_remd_indices

!*******************************************************************************
!
! Subroutine: write_hybridsolvent_remd_traj_binary_crds
!
! Description: Writes stripped coordinates in Netcdf format
!              for Hybrid Solvent REMD simulations and update frame counters
!
!*******************************************************************************

subroutine write_hybridsolvent_remd_traj_binary_crds(hybridsolvent_remd_crd, &
                        hybridsolvent_remd_atm_cnt, &
                        my_remd_method, my_replica_indexes, my_remd_dimension, &
                        my_remd_repidx, my_remd_crdidx)

#ifdef BINTRAJ
  use netcdf
  use AmberNetcdf_mod,   only : checkNCerror
  use axis_optimize_mod
# ifdef MPI
  use mdin_ctrl_dat_mod, only : temp0, solvph, t, solve
!  use remd_mod,          only : remd_method, replica_indexes, remd_dimension, &
!                                remd_repidx, remd_crdidx
  use sgld_mod,          only : trxsgld
# endif
#endif
  implicit none
  ! Formal arguments:
  integer, intent(in)          :: hybridsolvent_remd_atm_cnt
  double precision, intent(in) :: hybridsolvent_remd_crd(3, hybridsolvent_remd_atm_cnt)
  integer, intent(in)          :: my_remd_method
  integer, intent(in)          :: my_remd_dimension
  integer, intent(in)          :: my_replica_indexes(my_remd_dimension)
  integer, intent(in)          :: my_remd_repidx
  integer, intent(in)          :: my_remd_crdidx
!  double precision, intent(in) :: temp0

  ! Local variables
  double precision :: buf(3, hybridsolvent_remd_atm_cnt)
  integer          :: i
#ifdef BINTRAJ
  ! Local variables:
  integer :: ord1, ord2, ord3

  ord1 = axis_flipback_ords(1)
  ord2 = axis_flipback_ords(2)
  ord3 = axis_flipback_ords(3)

  if (ord1 .eq. 1 .and. ord2 .eq. 2) then
    call checkNCerror(nf90_put_var(hybridsolvent_remd_traj_ncid, &
                                   hybridsolvent_remd_traj_coord_var_id, &
                                   hybridsolvent_remd_crd(:,:), &
                                   start=(/ 1, 1, hybridsolvent_remd_traj_frame /), &
                                   count=(/ 3, hybridsolvent_remd_atm_cnt, 1 /)), &
                      'write stripped atom coords')
  else
    do i = 1, hybridsolvent_remd_atm_cnt
      buf(1, i) = hybridsolvent_remd_crd(ord1, i)
      buf(2, i) = hybridsolvent_remd_crd(ord2, i)
      buf(3, i) = hybridsolvent_remd_crd(ord3, i)
    end do

    call checkNCerror(nf90_put_var(hybridsolvent_remd_traj_ncid, &
                                   hybridsolvent_remd_traj_coord_var_id, &
                                   buf(:,:), start=(/ 1, 1, hybridsolvent_remd_traj_frame /), &
                                   count=(/ 3, hybridsolvent_remd_atm_cnt, 1 /)), &
                      'write stripped atom coords')
  end if

# ifdef MPI
  if (my_remd_method .ne. 0) then
    ! Store overall replica and coordinate indices

    call checkNCerror(nf90_put_var(hybridsolvent_remd_traj_ncid, hybridsolvent_remd_traj_repidx_var_id, my_remd_repidx, &
                                   start = (/ hybridsolvent_remd_traj_frame /)), &
                      'write overall replica index')
    call checkNCerror(nf90_put_var(hybridsolvent_remd_traj_ncid, hybridsolvent_remd_traj_crdidx_var_id, my_remd_crdidx, &
                                   start = (/ hybridsolvent_remd_traj_frame /)), &
                      'write overall coordinate index')
    ! multi-D remd: Store indices of this replica in each dimension
    if (my_remd_method .eq. -1) then
      call checkNCerror(nf90_put_var(hybridsolvent_remd_traj_ncid, hybridsolvent_remd_traj_remd_indices_var_id, &
                                     my_replica_indexes(:), &
                                     start = (/ 1, hybridsolvent_remd_traj_frame /), &
                                     count = (/ my_remd_dimension, 1 /)), &
                        'write replica index for each dimension')
    endif
    if (my_remd_method .eq. 4) then
      call checkNCerror(nf90_put_var(hybridsolvent_remd_traj_ncid, hybridsolvent_remd_traj_temp_var_id, solvph, &
                                     start = (/ hybridsolvent_remd_traj_frame /)), &
                        'write replica pH')
    else if (my_remd_method .eq. 5) then
      call checkNCerror(nf90_put_var(hybridsolvent_remd_traj_ncid, hybridsolvent_remd_traj_temp_var_id, solve, &
                                     start = (/ hybridsolvent_remd_traj_frame /)), &
                        'write replica redox potential')
    else if (trxsgld) then
      call checkNCerror(nf90_put_var(hybridsolvent_remd_traj_ncid, hybridsolvent_remd_traj_temp_var_id, &
                                     REAL(my_replica_indexes(1)), &
                                     start = (/ hybridsolvent_remd_traj_frame /)), &
                        'write SGLD replica index')
    else
      call checkNCerror(nf90_put_var(hybridsolvent_remd_traj_ncid, hybridsolvent_remd_traj_temp_var_id, temp0, &
                                     start = (/ hybridsolvent_remd_traj_frame /)), &
                        'write replica temperature')
    end if
  end if
# endif
  call checkNCerror(nf90_put_var(hybridsolvent_remd_traj_ncid, hybridsolvent_remd_traj_time_var_id, (/ t /), &
                    start=(/ hybridsolvent_remd_traj_frame /), count=(/ 1 /)), 'write time')

  call checkNCerror(nf90_sync(hybridsolvent_remd_traj_ncid))

  hybridsolvent_remd_traj_frame = hybridsolvent_remd_traj_frame + 1

#else
  call no_nc_error
#endif

  return

end subroutine write_hybridsolvent_remd_traj_binary_crds
#endif /* End #ifdef MPI */

end module hybridsolvent_remd_mod
