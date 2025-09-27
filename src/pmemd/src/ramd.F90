#include "copyright.i"

!*******************************************************************************
!
! Module:  ramd_mod
!
! Description: Random Acceleration Molecular Dynamics
!              
!*******************************************************************************

module ramd_mod

  use random_mod, only : random_state, amrand_gen, amrset_gen

  implicit none

  type(random_state),save          :: ramd_randgen
  integer, allocatable             :: ligregion(:)
  integer, allocatable             :: protregion(:)

  double precision  :: save_dist

  double precision  :: max_dist
  
  integer  :: forward_time

  integer :: max_time

  double precision :: cur_boost

  double precision :: forward_boost

  logical :: forward

  logical :: first_time

  contains

!*******************************************************************************
!
! Module:  init_ramd_mod
!
! Description: Initializes thing(s)
!              
!*******************************************************************************

  subroutine init_ramd

    use mdin_ctrl_dat_mod

    implicit none

    call amrset_gen(ramd_randgen, ig)

    ! Set save dist to something absurd so we actually run the code first step
    save_dist =9999999.0d0

    cur_boost = ramdboost

    first_time = .true.
 
    forward= .true.

  end subroutine init_ramd

!*******************************************************************************
!
! Module:  init_ramd_mask
!
! Description: Initializes ligand mask
!             
!*******************************************************************************

  subroutine init_ramd_mask(atm_cnt, nres, igraph, isymbl, res_atms, labres, crd)

      use mdin_ctrl_dat_mod
      use file_io_dat_mod
      use findmask_mod
#if defined(CUDA)
      use pmemd_lib_mod !I think we just need setup_alloc_error
#endif

      implicit none

      ! Formal arguments:
      integer, intent(in)             :: atm_cnt, nres
      integer, intent(in)             :: res_atms(nres)
      character(len=4), intent(in)    :: igraph(atm_cnt)
      character(len=4), intent(in)    :: isymbl(atm_cnt)
      character(len=4), intent(in)    :: labres(nres)
      double precision, intent(in)    :: crd(3 * atm_cnt)
#if defined(CUDA)
      integer                         :: alloc_failed
#endif

      if(.not. allocated(ligregion)) allocate(ligregion(atm_cnt))
      if(.not. allocated(protregion)) allocate(protregion(atm_cnt))

      if(ramdligmask .ne. '') then
        call atommask(atm_cnt, nres, 0, igraph, isymbl, res_atms, labres, &
               crd, ramdligmask, ligregion)
      end if

      if(ramdproteinmask .ne. '') then
        call atommask(atm_cnt, nres, 0, igraph, isymbl, res_atms, labres, &
               crd, ramdproteinmask, protregion)
      end if

  end subroutine init_ramd_mask

!*******************************************************************************
!
! Module:  set_ramd_random_accel
!
! Description: Sets random velocitiesk
!              
!*******************************************************************************

  subroutine set_ramd_random_accel(atm_cnt, ramdaccel)

    use mdin_ctrl_dat_mod, only : ramdboost

    implicit none

    integer          :: atm_cnt
    double precision, intent(out) :: ramdaccel(3)

    double precision :: rand1, rand2
    double precision :: theta
    double precision :: psi

    call amrand_gen(ramd_randgen, rand1)
    call amrand_gen(ramd_randgen, rand2)
    
    theta = rand1 * 2.0 * 3.1415926535897932
    psi = acos(1.0 -2*rand2)
    ramdaccel(1) = cos(theta)*sin(psi)*cur_boost
    ramdaccel(2) = sin(theta)*sin(psi)*cur_boost
    ramdaccel(3) = cos(psi)*cur_boost

  end subroutine set_ramd_random_accel

!*******************************************************************************
!
! Module:  set_ramd_com
!
! Description: TBD
!              
!*******************************************************************************

  subroutine set_ramd_com(atm_cnt, crd, mass, lig_com, prot_com)

    integer :: atm_cnt
    double precision :: crd(3, atm_cnt)
    double precision :: mass(atm_Cnt)
    double precision, intent(out) :: lig_com(3)
    double precision, intent(out) :: prot_com(3)

    double precision :: mass_prot, mass_lig
    integer :: i

    mass_prot = 0.0d0
    mass_lig = 0.0d0
    lig_com(:) = 0.d0
    prot_com(:) = 0.0d0

    do i=1, atm_cnt
        lig_com(1)  = lig_com(1) + crd(1, i) * mass(i) * ligregion(i)
        lig_com(2)  = lig_com(2) + crd(2, i) * mass(i) * ligregion(i)
        lig_com(3)  = lig_com(3) + crd(3, i) * mass(i) * ligregion(i)
        mass_lig = mass_lig + mass(i) * ligregion(i)
        prot_com(1)  = prot_com(1) + crd(1, i) * mass(i) * protregion(i)
        prot_com(2)  = prot_com(2) + crd(2, i) * mass(i) * protregion(i)
        prot_com(3)  = prot_com(3) + crd(3, i) * mass(i) * protregion(i)
        mass_prot = mass_prot + mass(i) * protregion(i)
    end do

    lig_com(:) = lig_com(:) / mass_lig    
    prot_com(:) = prot_com(:) / mass_prot

  end subroutine set_ramd_com

!*******************************************************************************
!
! Module:  run_ramd
!
! Description: runs ramd workflow
!              
!*******************************************************************************
  
  subroutine run_ramd(atm_cnt, crd, frc, mass, new_list, nstep)

    use file_io_dat_mod
    use pmemd_lib_mod
    use mdin_ctrl_dat_mod, only : ramdmaxdist, ramdboost, &
                ramdboostrate,ramdboostfreq

    implicit none

    integer :: atm_cnt
    double precision :: crd(3,atm_cnt)
    double precision :: mass(atm_cnt)
    double precision :: frc(3,atm_cnt)
    logical, intent(inout) :: new_list
    integer, intent(in) :: nstep

    double precision :: ligcom(3), protcom(3), ramdaccel(3)
    double precision :: dist
    integer i

#ifdef CUDA
    call gpu_download_crd(crd)    
#endif

    call set_ramd_com(atm_cnt, crd, mass, ligcom, protcom)
    
    dist = (ligcom(3) - protcom(3))**2 + (ligcom(2) - protcom(2))**2 + (ligcom(1) - protcom(1))**2
    
    if(first_time) then
        max_dist = (sqrt(dist)+ramdmaxdist)*(sqrt(dist)+ramdmaxdist)
        first_time = .false.
    end if
    
    write(mdout, '(a,f25.15,/,a,f25.15)' ) " RAMD Distance^2 Old: ", save_dist, &
                                           " RAMD Distance^2 New: ", dist
    
    ! Update cur_boost
    if(mod(nstep,ramdboostfreq) .eq. 0 .and. nstep .gt. 0) cur_boost = cur_boost + ramdboostrate
    
    if(dist .le. save_dist) then
#ifdef CUDA
      call gpu_download_frc(frc)
#endif
      call set_ramd_random_accel(atm_cnt, ramdaccel)

      do i=1, atm_cnt
        if(ligregion(i) .eq. 1) then
          frc(:,i)=frc(:,i)+ramdaccel(:)*mass(i)
        end if
      end do
  
      write(mdout, '(a,/,3g25.15,/,a,f34.15)' ) " RAMD Random Accel (X,Y,Z):", &
         ramdaccel(1),ramdaccel(2),ramdaccel(3), " RAMD Boost: ", cur_boost

#ifdef CUDA
      call gpu_upload_frc(frc)    
      call gpu_force_new_neighborlist()
#endif
      new_list = .true.
    end if
    
    save_dist = dist
 
    if(dist .gt. max_dist) then
      write(mdout,'(a)') "RAMD Forward: Reached maxed distance."
      forward=.false.
      forward_time=nstep
      forward_boost = cur_boost
      cur_boost=ramdboost
      write(mdout,'(a)') "RAMD Forward: Reached maxed distance."
      write(mdout,*)"Forward boost ended at:",forward_boost
      write(mdout,*) "Forward time:", forward_time
      call mexit(6,1)
    end if
       
  end subroutine run_ramd

end module ramd_mod
