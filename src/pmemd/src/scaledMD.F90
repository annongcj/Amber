#include "copyright.i"

!*******************************************************************************
!
! Module: scaledMD_mod
!
! Description: 
!
! Module for controlling accelerated molecular dynamics calculations
! Specifically for Scaled MD.
!
! Written by Romelia Salomon, 5/2013
!              
!*******************************************************************************

module scaledMD_mod

  use file_io_dat_mod
  use parallel_dat_mod
  use mdin_ctrl_dat_mod
  
  implicit none

! Everything is private by default

! scaledMD Variables
!
  double precision, save              ::     scaledMD_energy, scaledMD_unscaled_energy, scaledMD_weight

contains

!*******************************************************************************
!
! Subroutine: scaledMD_setup
!
! Description: Sets up the scaledMD run. Opens the log file, sets up exchange tables
!              etc.
!
!*******************************************************************************

subroutine scaledMD_setup(ntwx)
  
  use pmemd_lib_mod

  implicit none

! Formal arguments:
  integer               :: ntwx

#ifdef CUDA
  call gpu_scaledmd_setup(scaledMD, scaledMD_lambda)
#endif
  return
end subroutine scaledMD_setup


!*******************************************************************************
!
! Subroutine:  scaledMD_scale_frc
!
! Description: <TBS>
!              
!*******************************************************************************

#ifdef MPI
subroutine scaledMD_scale_frc(atm_cnt,tot_potenergy,frc,crd,my_atm_lst)
#else
subroutine scaledMD_scale_frc(atm_cnt,tot_potenergy,frc,crd)
#endif
  implicit none

! Formal arguments:
  integer                       :: atm_cnt
  double precision              :: tot_potenergy
  double precision              :: frc(3, atm_cnt)
  double precision              :: crd(3, atm_cnt)
#ifdef MPI
  integer                       :: my_atm_lst(*)
#endif

! Local variables:
  integer              :: atm_lst_idx,i

!scaledMD 
! Scaled energy
  scaledMD_energy = tot_potenergy * scaledMD_lambda
  scaledMD_weight = -(tot_potenergy * (1.0-scaledMD_lambda))
  scaledMD_unscaled_energy = tot_potenergy 

#ifdef MPI
  do atm_lst_idx = 1, my_atm_cnt
     i = my_atm_lst(atm_lst_idx)
#else
  do i = 1, atm_cnt
#endif
     frc(:,i) = frc(:,i) * scaledMD_lambda
  enddo 

  return

end subroutine scaledMD_scale_frc


!*******************************************************************************
!
! Subroutine:  write_amd_weights
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine write_scaledMD_log(ntwx,total_nstep)

  use file_io_mod

  implicit none

! Formal arguments:

  integer               :: ntwx
  integer               :: total_nstep

! Local variables:

    write(scaledMDlog,'(2x,2i10,4f22.12)') ntwx,total_nstep, &
     scaledMD_unscaled_energy,scaledMD_lambda,scaledMD_energy,scaledMD_weight
  
  return

end subroutine write_scaledMD_log


!*******************************************************************************
!
! Subroutine:  scaledMD_scale_frc_gb
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine scaledMD_scale_frc_gb(atm_cnt,tot_potenergy,frc,crd)

  use gb_parallel_mod

  implicit none

! Formal arguments:
  integer                       :: atm_cnt
  double precision              :: tot_potenergy
  double precision              :: frc(3, atm_cnt)
  double precision              :: crd(3, atm_cnt)

! Local variables:
  integer              :: atm_lst_idx,i

!scaledMD 
! Scaled energy
  scaledMD_energy = tot_potenergy * scaledMD_lambda
  scaledMD_weight = -(tot_potenergy * (1.0-scaledMD_lambda))
  scaledMD_unscaled_energy = tot_potenergy 

#ifdef MPI
  call gb_amd_apply_weight(atm_cnt,frc,scaledMD_lambda)
#else
  frc(:,:) = frc(:,:) * scaledMD_lambda
#endif

  return

end subroutine scaledMD_scale_frc_gb

end module scaledMD_mod
