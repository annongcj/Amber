!*******************************************************************************
!
! Subroutine:  get_efield_energy
!
! Description:
!              
! The main routine for electric field energy. 
!
!*******************************************************************************

subroutine get_efield_energy(img_frc, crd, img_qterm, img_atm_map, &
                         need_pot_enes, efield, atm_cnt, nstep)

  use gbl_constants_mod
  use timers_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use ene_frc_splines_mod
  use img_mod
  use pbc_mod

  implicit none

! Formal arguments:

  double precision, intent(in out) :: img_frc(3, *)
  double precision, intent(in)     :: crd(3, *)
  double precision , intent(in)    :: img_qterm(*)
  logical, intent(in)              :: need_pot_enes
  double precision, intent(out)    :: efield
  integer, intent(in)              :: atm_cnt
  integer, intent(in)              :: nstep

! Local variables and parameters:

  double precision      charge
  double precision      vx, vy, vz
  integer               img_i
  double precision      loc_efx, loc_efy, loc_efz
  double precision      efrcx,efrcy,efrcz
  double precision      phase 
  integer               img_atm_map(*)
  integer               crd_i
#ifdef MPI
  integer               simple_load_lo
  integer               simple_load_hi
#endif

  loc_efx=0.d0   !Stores time based electric field charge
  loc_efy=0.d0   
  loc_efz=0.d0
  
  efrcx=0.d0  ! Stores force on atom.
  efrcy=0.d0
  efrcz=0.d0

  efield = 0.d0
 
  vx=0.d0
  vy=0.d0
  vz=0.d0

#ifdef MPI
  simple_load_lo = int(mytaskid*float(atm_cnt/numtasks)) + 1
  simple_load_hi = int((mytaskid+1) * float(atm_cnt/numtasks))
#endif

  phase= cos((2*PI*effreq/1000)*(dt*float(nstep))-(pi/180*efphase))

  loc_efx = phase * efx
  loc_efy = phase * efy
  loc_efz = phase * efz

  !Normalize efield only works in a box.  Note: add trap

  if(efn .eq. 1) then
     loc_efx=loc_efx*recip(1,1) 
     loc_efy=loc_efy*recip(2,2)
     loc_efz=loc_efz*recip(3,3)
  end if  

  if (need_pot_enes) then

#ifdef MPI
    do img_i = simple_load_lo, simple_load_hi 
#else
    do img_i = 1, atm_cnt
#endif

        crd_i = img_atm_map(img_i)    !convert image array to crd array
        
        charge = img_qterm(img_i) * ONE_AMBER_ELECTROSTATIC

        efrcx = charge*loc_efx
        efrcy = charge*loc_efy
        efrcz = charge*loc_efz
   
        img_frc(1, img_i) = efrcx + img_frc(1, img_i)
        img_frc(2, img_i) = efrcy + img_frc(2, img_i)
        img_frc(3, img_i) = efrcz + img_frc(3, img_i)
     
        vx=crd(1,crd_i)-ucell(1,1)
        vy=crd(2,crd_i)-ucell(2,2)
        vz=crd(3,crd_i)-ucell(3,3)

        efield = efield - (vx)*efrcx - (vy)*efrcy - (vz)*efrcz 

    end do

  else ! do not need energies...

#ifdef MPI
    do img_i = simple_load_lo, simple_load_hi 
#else
    do img_i = 1, atm_cnt
#endif

        crd_i = img_atm_map(img_i)    !convert image array to crd array

        charge = img_qterm(img_i) / AMBER_ELECTROSTATIC

        efrcx = charge*loc_efx
        efrcy = charge*loc_efy
        efrcz = charge*loc_efz

        img_frc(1, img_i) = efrcx + img_frc(1, img_i)
        img_frc(2, img_i) = efrcy + img_frc(2, img_i)
        img_frc(3, img_i) = efrcz + img_frc(3, img_i)

    end do

  end if

  return

end subroutine get_efield_energy
