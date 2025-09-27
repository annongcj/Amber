#include "copyright.i"

!*******************************************************************************
!
! Module: gamd_mod
!
! Description: 
!
! Module for controlling Gaussian Accelerated Molecular Dynamics calculations
!
! Written by Yinglong Miao, yinglong.miao@gmail.com, 2014 - 2015
!              
!*******************************************************************************

module gamd_mod

  use file_io_dat_mod
  use gbl_constants_mod, only : KB
  use parallel_dat_mod
  use mdin_ctrl_dat_mod
  use dihedrals_mod

  implicit none

! Everything is private by default

! Variables
!

!GaMD
  integer, save                       :: num_gamd_recs, num_gamd_lag, gamd_ntwx,igamd0
  double precision, allocatable, save :: gamd_weights_and_energy(:)
  double precision, save              :: E_dih_boost ! Dihedral boost energy in kcal/mol
  double precision, save              :: E_total_boost ! Total Potential boost energy in kcal/mol
  double precision, save              :: E_nb_boost ! Non-bonded Potential boost energy in kcal/mol
  double precision, save              :: E_ti_region_boost ! Non-bonded Potential boost energy in kcal/mol
  double precision, save              :: ppi_inter_ene       !ppi_interaction between 2 regions
  double precision, save              :: ppi_dihedral_ene    ! saving the second potential
  double precision, save              :: ppi_bond_ene    ! saving the third bonded part potential

contains

!*******************************************************************************
!
! Subroutine: gamd_setup
!
! Description: Sets up the GaMD run. Opens the log file, sets up exchange tables
!              etc.
!
!*******************************************************************************

subroutine gamd_setup(ntwx)
  
  use pmemd_lib_mod

  implicit none

! Formal arguments:
  integer               :: ntwx

! Local variables:
  integer :: alloc_failed

  if (master) then
    allocate(gamd_weights_and_energy(6), stat = alloc_failed)
    if (alloc_failed .ne. 0) call setup_alloc_error
    gamd_weights_and_energy(1) = 0.0d0 
    gamd_weights_and_energy(2) = 0.0d0
    gamd_weights_and_energy(3) = 1.0d0
    gamd_weights_and_energy(4) = 1.0d0
    gamd_weights_and_energy(5) = 0.0d0
    gamd_weights_and_energy(6) = 0.0d0
  end if

  fwgtd = 1.0d0
  E_dih_boost = 0.0d0
  E_total_boost = 0.0d0
  E_nb_boost = 0.0d0
  E_ti_region_boost = 0.0d0
  num_gamd_recs = 1
  num_gamd_lag = 0
  gamd_ntwx = ntwx

#ifdef CUDA
  if(.not.(igamd.ge.20.and.igamd.le.28))then
  call gpu_gamd_setup(igamd, igamdlag, gamd_ntwx, EthreshP, kP, EthreshD, kD, temp0)
  else
  call gpu_gamd_ppi_setup(igamd, igamdlag, gamd_ntwx, EthreshP,kP,EthreshD,kD,EthreshB,k_B,temp0)
  endif
#endif

  if (master) then
    write(gamdlog,'(2x,a)')"# Gaussian Accelerated Molecular Dynamics log file"
    write(gamdlog,'(2x,a)')"# All energy terms stored in units of kcal/mol"
    write(gamdlog,'(2x,a)')"# ntwx,total_nstep,Unboosted-Potential-Energy, &
              &Unboosted-Dihedral-Energy,Total-Force-Weight,Dihedral-Force-Weight,Boost-Energy-Potential,Boost-Energy-Dihedral"
  end if

  return
end subroutine gamd_setup

!*******************************************************************************
!
! Subroutine:  calculate_gamd_dih_weights
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine calculate_gamd_dih_weights(atm_cnt,totdih_ene,frc,crd)

  implicit none

! Formal arguments:
  integer                       :: atm_cnt
  double precision              :: totdih_ene
  double precision              :: frc(3, atm_cnt)
  double precision              :: crd(3, atm_cnt)
! Local variables
  double precision              :: EV

! Calculate the boosting weight for gamd
  fwgtd   = 1.0d0
  E_dih_boost = 0.0d0

  EV = (EthreshD - totdih_ene)
  if( EV .gt. 0.0d0 )then
    if(num_gamd_lag .eq. 0)then
      E_dih_boost = 0.5 * kD * (EV * EV); ! Dih boost E in kcal/mol
      fwgtd  = 1.0 - kD * EV ;
      ! write(*,'(a,2f15.5)') "GaMD) E_dih_boost, fwgtd = ", E_dih_boost, fwgtd
    end if
  end if
  if (ntf .le. 5) then
    if (cit_nphih .gt. 0) then
      call get_dihed_forces_gamd(cit_nphih, cit_dihed, crd, frc)
    end if
  end if

  if (ntf .le. 6) then
    if (cit_nphia .gt. 0) then
      call get_dihed_forces_gamd(cit_nphia, cit_dihed(cit_nphih + 1), crd, frc)
    end if
  end if

  return

end subroutine calculate_gamd_dih_weights

!*******************************************************************************
!
! Subroutine:  calculate_gamd_total_weights
!
! Description: <TBS>
!              
!*******************************************************************************

#ifdef MPI
subroutine calculate_gamd_total_weights(atm_cnt,tot_potenergy,totdih_ene,&
                                       E_boost,frc,crd,my_atm_lst)
#else
subroutine calculate_gamd_total_weights(atm_cnt,tot_potenergy,totdih_ene,&
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

!GaMD DUAL BOOST CALC START
   if(igamd.gt.0)then
     E_total_boost = 0.0d0
     fwgt   = 1.0d0
     if ((igamd == 1).or.(igamd == 3)) totalenergy = tot_potenergy + E_dih_boost
     EV = (EthreshP - totalenergy )
     if (((igamd == 1).or.(igamd == 3)) .and. (totalenergy.lt.EthreshP)) then
       if (num_gamd_lag .eq. 0) then
       E_total_boost = 0.5 * kP * (EV*EV); ! PE boost in Kcal/mol
       fwgt = 1.0 - kP * EV ;
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
     if (master.and. (num_gamd_recs.eq.gamd_ntwx)) then
       gamd_weights_and_energy(1) = tot_potenergy 
       gamd_weights_and_energy(2) = totdih_ene
       gamd_weights_and_energy(3) = fwgt
       gamd_weights_and_energy(4) = fwgtd
       gamd_weights_and_energy(5) = E_total_boost
       gamd_weights_and_energy(6) = E_dih_boost
     endif
     E_boost = E_total_boost + E_dih_boost
     num_gamd_recs = num_gamd_recs +1
     if(num_gamd_lag .eq. igamdlag)then
       num_gamd_lag = 0
     else
       num_gamd_lag = num_gamd_lag +1
     endif 
   end if

  return

end subroutine calculate_gamd_total_weights


!*******************************************************************************
!
! Subroutine:  calculate_gamd_nb_weights
!
! Description: <TBS>
!              
!*******************************************************************************

#ifdef MPI
subroutine calculate_gamd_nb_weights(atm_cnt,nb_potenergy,totdih_ene,&
                                       E_boost,nb_frc,crd,my_atm_lst)
#else
subroutine calculate_gamd_nb_weights(atm_cnt,nb_potenergy,totdih_ene,&
                                       E_boost,nb_frc,crd)
#endif
  implicit none

! Formal arguments:
  integer                       :: atm_cnt
  double precision              :: E_nb_boost
  double precision              :: nb_potenergy
  double precision              :: totdih_ene
  double precision, intent(out) :: E_boost
  double precision              :: nb_frc(3, atm_cnt)
  double precision              :: crd(3, atm_cnt)
#ifdef MPI
  integer                       :: my_atm_lst(*)
#endif

! Local variables:
  integer              :: atm_lst_idx,i
  double precision     :: EV, ONE_KB, fwgt

!GaMD DUAL BOOST CALC START
   if(igamd.gt.0)then
     E_nb_boost = 0.0d0
     fwgt   = 1.0d0
     EV = (EthreshP - nb_potenergy)
     if (((igamd == 4).or.(igamd == 5)) .and. (nb_potenergy.lt.EthreshP)) then
       if (num_gamd_lag .eq. 0) then
       E_nb_boost = 0.5 * kP * (EV*EV); ! PE boost in Kcal/mol
       fwgt = 1.0 - kP * EV ;
#ifdef MPI
       do atm_lst_idx = 1, my_atm_cnt
         i = my_atm_lst(atm_lst_idx)
#else
       do i = 1, atm_cnt
#endif
         nb_frc(:,i) = nb_frc(:,i)*fwgt
       enddo 
       end if
     end if
     if (master.and. (num_gamd_recs.eq.gamd_ntwx)) then
       gamd_weights_and_energy(1) = nb_potenergy 
       gamd_weights_and_energy(2) = totdih_ene
       gamd_weights_and_energy(3) = fwgt
       gamd_weights_and_energy(4) = fwgtd
       gamd_weights_and_energy(5) = E_nb_boost
       gamd_weights_and_energy(6) = E_dih_boost
     endif
     E_boost = E_nb_boost + E_dih_boost
     num_gamd_recs = num_gamd_recs +1
     if(num_gamd_lag .eq. igamdlag)then
       num_gamd_lag = 0
     else
       num_gamd_lag = num_gamd_lag +1
     endif 
   end if

  return

end subroutine calculate_gamd_nb_weights


!*******************************************************************************
!
! Subroutine:  calculate_gamd_nb_weights_gb
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine calculate_gamd_nb_weights_gb(atm_cnt,nb_potenergy,totdih_ene,&
                                          E_boost,nb_frc,crd)

  use gb_parallel_mod

  implicit none

! Formal arguments:
  integer                       :: atm_cnt
  double precision              :: E_nb_boost
  double precision              :: nb_potenergy
  double precision              :: totdih_ene
  double precision, intent(out) :: E_boost
  double precision              :: nb_frc(3, atm_cnt)
  double precision              :: crd(3, atm_cnt)

! Local variables:
  integer              :: atm_lst_idx,i
  double precision     :: EV, ONE_KB, fwgt

!GaMD DUAL BOOST CALC START
   if(igamd.gt.0)then
     E_nb_boost = 0.0d0
     fwgt   = 1.0d0
     EV = (EthreshP - nb_potenergy)
     if (((igamd == 4).or.(igamd == 5)) .and. (nb_potenergy.lt.EthreshP)) then
       if (num_gamd_lag .eq. 0) then
       E_nb_boost = 0.5 * kP * (EV*EV); ! PE boost in kcal/mol
       fwgt = 1.0 - kP * EV ;
#ifdef MPI
       call gb_gamd_apply_weight(atm_cnt,nb_frc,fwgt)
#else
       nb_frc(:,:) = nb_frc(:,:)*fwgt
#endif
       end if
     end if
     if (master.and. (num_gamd_recs.eq.gamd_ntwx)) then
       gamd_weights_and_energy(1) = nb_potenergy 
       gamd_weights_and_energy(2) = totdih_ene
       gamd_weights_and_energy(3) = fwgt
       gamd_weights_and_energy(4) = fwgtd
       gamd_weights_and_energy(5) = E_nb_boost
       gamd_weights_and_energy(6) = E_dih_boost
     endif
     E_boost = E_nb_boost + E_dih_boost
     num_gamd_recs = num_gamd_recs +1
     if(num_gamd_lag .eq. igamdlag)then
       num_gamd_lag = 0
     else
       num_gamd_lag = num_gamd_lag +1
     endif 
   end if

  return

end subroutine calculate_gamd_nb_weights_gb


!*******************************************************************************
!
! Subroutine:  calculate_gamd_ti_region_weights
!
! Description: <TBS>
!              
!*******************************************************************************

#ifdef MPI
subroutine calculate_gamd_ti_region_weights(atm_cnt,ti_region_potenergy,totdih_ene,&
                                       E_boost,fwgt,frc,crd,ti_sc_lst,my_atm_lst)
#else
subroutine calculate_gamd_ti_region_weights(atm_cnt,ti_region_potenergy,totdih_ene,&
                                       E_boost,fwgt,frc,crd,ti_sc_lst)
#endif
  implicit none

! Formal arguments:
  integer                       :: atm_cnt
  double precision              :: E_ti_region_boost
  double precision              :: ti_region_potenergy, fwgt
  double precision              :: totdih_ene
  double precision, intent(out) :: E_boost
  double precision              :: crd(3, atm_cnt)
  double precision              :: frc(3, atm_cnt)
  integer                       :: ti_sc_lst(atm_cnt)
#ifdef MPI
  integer                       :: my_atm_lst(*)
#endif

! Local variables:
  integer              :: atm_lst_idx,i
  double precision     :: EV, ONE_KB

!GaMD DUAL BOOST CALC START
   if(igamd.gt.0)then
     E_ti_region_boost = 0.0d0
     fwgt   = 1.0d0
     EV = (EthreshP - ti_region_potenergy)
     if (((igamd == 6).or.(igamd == 7)) .and. (ti_region_potenergy.lt.EthreshP)) then
       if (num_gamd_lag .eq. 0) then
       E_ti_region_boost = 0.5 * kP * (EV*EV); ! PE boost in Kcal/mol
       fwgt = 1.0 - kP * EV ;
       ! write(*,'(a,2f15.5)') "GaMD) E_ti_region_boost, fwgt = ", E_ti_region_boost, fwgt
#ifdef MPI
       do atm_lst_idx = 1, my_atm_cnt
         i = my_atm_lst(atm_lst_idx)
#else
       do i = 1, atm_cnt
#endif
         if (ti_sc_lst(i).ne.0) frc(:,i) = frc(:,i)*fwgt
       enddo 
       end if
     end if
     if (master.and. (num_gamd_recs.eq.gamd_ntwx)) then
       gamd_weights_and_energy(1) = ti_region_potenergy 
       gamd_weights_and_energy(2) = totdih_ene
       gamd_weights_and_energy(3) = fwgt
       gamd_weights_and_energy(4) = fwgtd
       gamd_weights_and_energy(5) = E_ti_region_boost
       gamd_weights_and_energy(6) = E_dih_boost
     endif
     E_boost = E_ti_region_boost + E_dih_boost
     num_gamd_recs = num_gamd_recs +1
     if(num_gamd_lag .eq. igamdlag)then
       num_gamd_lag = 0
     else
       num_gamd_lag = num_gamd_lag +1
     endif 
   end if

!   write(*,'(a,6f22.12)') "GaMD) gamd_weights: ", &
!    gamd_weights_and_energy(1), gamd_weights_and_energy(2), &
!    gamd_weights_and_energy(3), gamd_weights_and_energy(4), &
!    gamd_weights_and_energy(5), gamd_weights_and_energy(6)

  return

end subroutine calculate_gamd_ti_region_weights


!*******************************************************************************
!
! Subroutine:  calculate_gamd_ti_region_weights_gb
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine calculate_gamd_ti_region_weights_gb(atm_cnt,ti_region_potenergy,totdih_ene,&
                                          E_boost,fwgt,frc,crd,ti_sc_lst)

  use gb_parallel_mod

  implicit none

! Formal arguments:
  integer                       :: atm_cnt
  double precision              :: E_ti_region_boost
  double precision              :: ti_region_potenergy, fwgt
  double precision              :: totdih_ene
  double precision, intent(out) :: E_boost
  double precision              :: frc(3, atm_cnt)
  double precision              :: crd(3, atm_cnt)
  integer                       :: ti_sc_lst(atm_cnt)

! Local variables:
  integer              :: atm_lst_idx,i
  double precision     :: EV, ONE_KB

!GaMD DUAL BOOST CALC START
   if(igamd.gt.0)then
     E_ti_region_boost = 0.0d0
     fwgt   = 1.0d0
     EV = (EthreshP - ti_region_potenergy)
     if (((igamd == 4).or.(igamd == 5)) .and. (ti_region_potenergy.lt.EthreshP)) then
       if (num_gamd_lag .eq. 0) then
       E_ti_region_boost = 0.5 * kP * (EV*EV); ! PE boost in kcal/mol
       fwgt = 1.0 - kP * EV ;
#ifdef MPI
       call gb_gamd_apply_weight(atm_cnt,frc,fwgt)
#else
       if (ti_sc_lst(i).ne.0) frc(:,i) = frc(:,i)*fwgt
#endif
       end if
     end if
     if (master.and. (num_gamd_recs.eq.gamd_ntwx)) then
       gamd_weights_and_energy(1) = ti_region_potenergy 
       gamd_weights_and_energy(2) = totdih_ene
       gamd_weights_and_energy(3) = fwgt
       gamd_weights_and_energy(4) = fwgtd
       gamd_weights_and_energy(5) = E_ti_region_boost
       gamd_weights_and_energy(6) = E_dih_boost
     endif
     E_boost = E_ti_region_boost + E_dih_boost
     num_gamd_recs = num_gamd_recs +1
     if(num_gamd_lag .eq. igamdlag)then
       num_gamd_lag = 0
     else
       num_gamd_lag = num_gamd_lag +1
     endif 
   end if

  return

end subroutine calculate_gamd_ti_region_weights_gb


!*******************************************************************************
!
! Subroutine:  calculate_gamd_total_weights_gb
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine calculate_gamd_total_weights_gb(atm_cnt,tot_potenergy,totdih_ene,&
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

!GaMD DUAL BOOST CALC START
   if(igamd.gt.0)then
     E_total_boost = 0.0d0
     fwgt   = 1.0d0
     if ((igamd == 1).or.(igamd == 3)) totalenergy = tot_potenergy + E_dih_boost 
     EV = (EthreshP - totalenergy)
     if (((igamd == 1).or.(igamd == 3)) .and. (totalenergy.lt.EthreshP)) then
       if (num_gamd_lag .eq. 0) then
       E_total_boost = 0.5 * kP * (EV*EV); ! PE boost in kcal/mol
       fwgt = 1.0 - kP * EV ;
#ifdef MPI
       call gb_gamd_apply_weight(atm_cnt,frc,fwgt)
#else
       frc(:,:) = frc(:,:)*fwgt
#endif
       end if
     end if
     if (master.and. (num_gamd_recs.eq.gamd_ntwx)) then
       gamd_weights_and_energy(1) = tot_potenergy 
       gamd_weights_and_energy(2) = totdih_ene
       gamd_weights_and_energy(3) = fwgt
       gamd_weights_and_energy(4) = fwgtd
       gamd_weights_and_energy(5) = E_total_boost
       gamd_weights_and_energy(6) = E_dih_boost
     endif
     E_boost = E_total_boost + E_dih_boost
     num_gamd_recs = num_gamd_recs +1
     if(num_gamd_lag .eq. igamdlag)then
       num_gamd_lag = 0
     else
       num_gamd_lag = num_gamd_lag +1
     endif 
   end if

  return

end subroutine calculate_gamd_total_weights_gb


!*******************************************************************************
!
! Subroutine:  write_gamd_weights
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine write_gamd_weights(ntwx,total_nstep)

  use file_io_mod
!  use mdin_ctrl_dat_mod

  implicit none

! Formal arguments:

  integer               :: ntwx
  integer               :: total_nstep

! Local variables:

    write(gamdlog,'(2x,2i10,6f22.12)') ntwx,total_nstep, &
     gamd_weights_and_energy(1), gamd_weights_and_energy(2), &
     gamd_weights_and_energy(3), gamd_weights_and_energy(4), &
     gamd_weights_and_energy(5), gamd_weights_and_energy(6)

  num_gamd_recs = 1

  return

end subroutine write_gamd_weights

pure double precision function r2(v)

   implicit none

   double precision, intent(in) :: v(3)

   r2 = v(1)**2 + v(2)**2 + v(3)**2

end function r2

pure double precision function norm3(v)

   implicit none

   double precision, intent(in) :: v(3)

   norm3 = sqrt(v(1)**2 + v(2)**2 + v(3)**2)

end function norm3

pure double precision function distance(r1, r2)

   implicit none

   double precision, intent(in) :: r1(3), r2(3)

   distance = norm3(r1 - r2)

end function distance

subroutine PBC_r2f(r, f)

   use pbc_mod, only : recip 

   implicit none

   double precision, intent(in) :: r(3)
   double precision, intent(out) :: f(3)

   f(1) = r(1)*recip(1,1) + r(2)*recip(2,1) + r(3)*recip(3,1)
   f(2) = r(1)*recip(1,2) + r(2)*recip(2,2) + r(3)*recip(3,2)
   f(3) = r(1)*recip(1,3) + r(2)*recip(2,3) + r(3)*recip(3,3)

   f(1) = f(1) - floor(f(1))
   f(2) = f(2) - floor(f(2))
   f(3) = f(3) - floor(f(3))

end subroutine PBC_r2f

!=============================================================================

subroutine PBC_f2r(f, r, t1, t2, t3)

   use pbc_mod, only : ucell

   implicit none

   double precision, intent(in) :: f(3)
   double precision, intent(out) :: r(3)

   integer, intent(in) :: t1, t2, t3

   r(1) = &
      (f(1) + t1)*ucell(1,1) + (f(2) + t2)*ucell(1,2) + (f(3) + t3)*ucell(1,3)
   r(2) = &
      (f(1) + t1)*ucell(2,1) + (f(2) + t2)*ucell(2,2) + (f(3) + t3)*ucell(2,3)
   r(3) = &
      (f(1) + t1)*ucell(3,1) + (f(2) + t2)*ucell(3,2) + (f(3) + t3)*ucell(3,3)

end subroutine PBC_f2r

function PBC_distance(r1, r2) result(value)

   implicit none

   double precision :: value
   double precision, intent(in) :: r1(*), r2(*)

   double precision :: f1(3), f2(3), d(-1:1,-1:1,-1:1)
   double precision :: x1(3), x2(3), d_min

   integer :: i, j, k

   ! get the fractionals

   call PBC_r2f(r1(1:3), f1)
   call PBC_r2f(r2(1:3), f2)

   ! wrap r1 to the primary cell

   call PBC_f2r(f1, x1, 0, 0, 0)
   d_min = 1.0d99

   ! wrap r2 to the primary/neighboring cells

   do i = -1, 1
      do j = -1, 1
         do k = -1, 1
            call PBC_f2r(f2, x2, i, j, k)
            d(i,j,k) = distance(x1, x2)
            ! write(*,*) "i,j,k,d = ", i,j,k,d(i,j,k)
            if (d(i,j,k).lt.d_min) &
               d_min = d(i,j,k)
         end do
      end do
   end do

   value = d_min

end function PBC_distance

double precision function msd(lcrd0, lcrd1)

   implicit none

   integer :: i
   double precision, intent(in) :: lcrd0(:,:), lcrd1(:,:)

   msd = 0.0d0
   do i = 1, size(lcrd0,2)
     msd = msd + r2(lcrd1(:,i) - lcrd0(:,i))
   end do
   msd = msd/size(lcrd0,2)

end function msd

end module gamd_mod
