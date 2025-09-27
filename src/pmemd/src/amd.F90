#include "copyright.i"

!*******************************************************************************
!
! Module: amd_mod
!
! Description: 
!
! Module for controlling accelerated molecular dynamics calculations
!
! Written by Romelia Salomon, 5/2011
!              
!*******************************************************************************

module amd_mod

  use file_io_dat_mod
  use gbl_constants_mod, only : KB
  use parallel_dat_mod
  use mdin_ctrl_dat_mod
  use dihedrals_mod

  implicit none

! Everything is private by default

! Variables
!

!AMD
  integer, save                       :: num_amd_recs, num_amd_lag, amd_ntwx
  double precision, allocatable, save :: amd_weights_and_energy(:)
  double precision, save              :: E_dih_boost ! Dihedral boost energy in kcal/mol
  double precision, save              :: E_total_boost ! Total Potential boost energy in kcal/mol
  double precision, save              :: w_sign ! Sign to switch from raising valleys to lowering barriers

contains

!*******************************************************************************
!
! Subroutine: amd_setup
!
! Description: Sets up the AMD run. Opens the log file, sets up exchange tables
!              etc.
!
!*******************************************************************************

subroutine amd_setup(ntwx)
  
  use pmemd_lib_mod

  implicit none

! Formal arguments:
  integer               :: ntwx

! Local variables:
  integer :: alloc_failed

  if (master) then
    allocate(amd_weights_and_energy(6), stat = alloc_failed)
    if (alloc_failed .ne. 0) call setup_alloc_error
  end if

  fwgtd = 1.0d0
  E_dih_boost = 0.0d0
  num_amd_recs = 1
  num_amd_lag = 0
  amd_ntwx = ntwx

  w_sign = 1.0d0
  if(w_amd .ne. 0) w_sign = -1.0d0

#ifdef CUDA
  call gpu_amd_setup(iamd, w_amd, iamdlag, ntwx, EthreshP, alphaP, EthreshD, alphaD, EthreshP_w, alphaP_w, EthreshD_w, &
  alphaD_w, w_sign, temp0)
#endif

  if (master) then
    write(amdlog,'(2x,a)')"# Accelerated Molecular Dynamics log file"
    write(amdlog,'(2x,a)')"# All energy terms stored in units of kcal/mol"
    write(amdlog,'(2x,a)')"# ntwx,total_nstep,Unboosted-Potential-Energy, &
              &Unboosted-Dihedral-Energy,Total-Force-Weight,Dihedral-Force-Weight,Boost-Energy-Potential,Boost-Energy-Dihedral"
  end if

  return
end subroutine amd_setup

!*******************************************************************************
!
! Subroutine:  calculate_amd_dih_weights
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine calculate_amd_dih_weights(atm_cnt,totdih_ene,frc,crd)

  implicit none

! Formal arguments:
  integer                       :: atm_cnt
  double precision              :: totdih_ene
  double precision              :: frc(3, atm_cnt)
  double precision              :: crd(3, atm_cnt)
! Local variables
  double precision              :: EV
  double precision              :: windowed_factor, windowed_derivative
  double precision              :: windowed_exp
! Calculate the boosting weight for amd

  fwgtd   = 1.0d0
  E_dih_boost = 0.0d0

  EV = (EthreshD - totdih_ene) * w_sign 
  if( EV .gt. 0.0d0 )then
!  if(totdih_ene.le.EthreshD)then
    if(num_amd_lag .eq. 0)then
      !EV = EthreshD - totdih_ene
      E_dih_boost = (EV*EV) / (alphaD + EV) ! Dih boost E in kcal/mol
      fwgtd = (alphaD**2)/((alphaD + EV)**2)
      if(w_amd .ne.0 )then
        windowed_exp = exp((EthreshD_w - totdih_ene) * w_sign/alphaD_w)
        windowed_factor = 1.0d0/(1.0d0 + windowed_exp)
        windowed_derivative =  E_dih_boost * windowed_factor * windowed_factor * windowed_exp / alphaD_w
        E_dih_boost = E_dih_boost * windowed_factor
        fwgtd = fwgtd * windowed_factor + windowed_derivative
      endif
    end if
  end if
  if (ntf .le. 5) then
    if (cit_nphih .gt. 0) then
      call get_dihed_forces_amd(cit_nphih, cit_dihed, crd, frc)
    end if
  end if

  if (ntf .le. 6) then
    if (cit_nphia .gt. 0) then
      call get_dihed_forces_amd(cit_nphia, cit_dihed(cit_nphih + 1), crd, frc)
    end if
  end if

  return

end subroutine calculate_amd_dih_weights

!*******************************************************************************
!
! Subroutine:  calculate_amd_total_weights
!
! Description: <TBS>
!              
!*******************************************************************************

#ifdef MPI
subroutine calculate_amd_total_weights(atm_cnt,tot_potenergy,totdih_ene,&
                                       E_boost,frc,crd,my_atm_lst)
#else
subroutine calculate_amd_total_weights(atm_cnt,tot_potenergy,totdih_ene,&
                                       E_boost,frc,crd)
#endif
  implicit none

! Formal arguments:
  integer                       :: atm_cnt
  double precision              :: tot_potenergy
  double precision              :: totdih_ene
  double precision, intent(out) :: E_boost
  double precision              :: frc(3, atm_cnt)
  double precision              :: crd(3, atm_cnt)
#ifdef MPI
  integer                       :: my_atm_lst(*)
#endif

! Local variables:
  integer              :: atm_lst_idx,i
  double precision     :: EV, ONE_KB, totalenergy, fwgt 
  double precision     :: windowed_factor, windowed_derivative
  double precision     :: windowed_exp


!AMD DUAL BOOST CALC START
   if(iamd.gt.0)then
     E_total_boost = 0.0d0
     fwgt   = 1.0d0
     totalenergy = tot_potenergy + E_dih_boost
     EV = (EthreshP - totalenergy ) * w_sign
     if (((iamd == 1).or.(iamd == 3)) .and. ( EV .gt. 0.0d0 )) then
!     if (((iamd == 1).or.(iamd == 3)) .and. (totalenergy.le.EthreshP)) then
       if (num_amd_lag .eq. 0) then
         !EV = EthreshP - totalenergy 
         E_total_boost = ((EV*EV) / (alphaP + EV)) ! PE boost in Kcal/mol 
         fwgt = (alphaP**2)/((alphaP + EV)**2)
         if(w_amd .ne.0 )then
           windowed_exp = exp((EthreshP_w - totalenergy) * w_sign/alphaP_w)
           windowed_factor = 1.0d0/(1.0d0 + windowed_exp)
           windowed_derivative =  E_total_boost * windowed_factor * windowed_factor * windowed_exp / alphaP_w
           E_total_boost = E_total_boost * windowed_factor
           fwgt = fwgt * windowed_factor + windowed_derivative
         endif
#ifdef MPI
         do atm_lst_idx = 1, my_atm_cnt
           i = my_atm_lst(atm_lst_idx)
#else
         do i = 1, atm_cnt
#endif
           frc(:,i) = frc(:,i)*fwgt
         enddo 
       end if
     end if
     if (master.and. (num_amd_recs.eq.amd_ntwx)) then
       ONE_KB = 1.0d0 / (temp0 * KB)
       amd_weights_and_energy(1) = tot_potenergy  
       amd_weights_and_energy(2) = totdih_ene
       amd_weights_and_energy(3) = fwgt
       amd_weights_and_energy(4) = fwgtd
       amd_weights_and_energy(5) = E_total_boost
       amd_weights_and_energy(6) = E_dih_boost
!       amd_weights_and_energy(5) = E_total_boost * ONE_KB
!       amd_weights_and_energy(6) = E_dih_boost * ONE_KB
     endif
     E_boost = E_total_boost + E_dih_boost
     num_amd_recs = num_amd_recs +1
     if(num_amd_lag .eq. iamdlag)then
       num_amd_lag = 0
     else
       num_amd_lag = num_amd_lag +1
     endif 
   end if

  return

end subroutine calculate_amd_total_weights


!*******************************************************************************
!
! Subroutine:  calculate_amd_total_weights_gb
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine calculate_amd_total_weights_gb(atm_cnt,tot_potenergy,totdih_ene,&
                                          E_boost,frc,crd)

  use gb_parallel_mod

  implicit none

! Formal arguments:
  integer                       :: atm_cnt
  double precision              :: tot_potenergy
  double precision              :: totdih_ene
  double precision, intent(out) :: E_boost
  double precision              :: frc(3, atm_cnt)
  double precision              :: crd(3, atm_cnt)

! Local variables:
  integer              :: atm_lst_idx,i
  double precision     :: EV, ONE_KB, totalenergy, fwgt
  double precision     :: windowed_factor, windowed_derivative
  double precision     :: windowed_exp

!AMD DUAL BOOST CALC START
   if(iamd.gt.0)then
     E_total_boost = 0.0d0
     fwgt   = 1.0d0
     totalenergy = tot_potenergy + E_dih_boost 
     EV = (EthreshP - totalenergy ) * w_sign
     if (((iamd == 1).or.(iamd == 3)) .and. ( EV .gt. 0.0d0 )) then
!     if (((iamd == 1).or.(iamd == 3)) .and. (totalenergy.le.EthreshP)) then
       if (num_amd_lag .eq. 0) then
         !EV = EthreshP - totalenergy 
         E_total_boost = ((EV*EV) / (alphaP + EV)) ! PE boost in kcal/mol 
         fwgt = (alphaP**2)/((alphaP + EV)**2)
         if(w_amd .ne.0 )then
           windowed_exp = exp((EthreshP_w - totalenergy) * w_sign/alphaP_w)
           windowed_factor = 1.0d0/(1.0d0 + windowed_exp)
           windowed_derivative =  E_total_boost * windowed_factor * windowed_factor * windowed_exp / alphaP_w
           E_total_boost = E_total_boost * windowed_factor
           fwgt = fwgt * windowed_factor + windowed_derivative
         endif
#ifdef MPI
         call gb_amd_apply_weight(atm_cnt,frc,fwgt)
#else
         frc(:,:) = frc(:,:)*fwgt
#endif
       end if
     end if
     if (master.and. (num_amd_recs.eq.amd_ntwx)) then
       ONE_KB = 1.0d0 / (temp0 * KB)
       amd_weights_and_energy(1) = tot_potenergy 
       amd_weights_and_energy(2) = totdih_ene
       amd_weights_and_energy(3) = fwgt
       amd_weights_and_energy(4) = fwgtd
       amd_weights_and_energy(5) = E_total_boost
       amd_weights_and_energy(6) = E_dih_boost
!       amd_weights_and_energy(5) = E_total_boost * ONE_KB
!       amd_weights_and_energy(6) = E_dih_boost * ONE_KB
     endif
     E_boost = E_total_boost + E_dih_boost
     num_amd_recs = num_amd_recs +1
     if(num_amd_lag .eq. iamdlag)then
       num_amd_lag = 0
     else
       num_amd_lag = num_amd_lag +1
     endif 
   end if

  return

end subroutine calculate_amd_total_weights_gb


!*******************************************************************************
!
! Subroutine:  write_amd_weights
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine write_amd_weights(ntwx,total_nstep)

  use file_io_mod
!  use mdin_ctrl_dat_mod

  implicit none

! Formal arguments:

  integer               :: ntwx
  integer               :: total_nstep

! Local variables:

    write(amdlog,'(2x,2i10,6f22.12)') ntwx,total_nstep, &
     amd_weights_and_energy(1), amd_weights_and_energy(2), &
     amd_weights_and_energy(3), amd_weights_and_energy(4), &
     amd_weights_and_energy(5), amd_weights_and_energy(6)
  
  num_amd_recs = 1

  return

end subroutine write_amd_weights

end module amd_mod
