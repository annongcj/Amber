
!*******************************************************************************
!
! Module:  reaf_mod
!
! Description: Module for Replica-Exchange with Arbitrary Degrees if freedom (REAF) support. 
!              Written by Taisung Lee (Rutgers, 2020)
!              
!*******************************************************************************

module reaf_mod

#ifdef GTI

  use file_io_dat_mod
  use state_info_mod
  implicit none
  
  logical, save :: initialized=.false.

  integer, save                            :: reaf_mode = -1  ! 0: on Lambda=0 ; 1: on Lambda=1
  integer, dimension(:,:),allocatable,save :: reaf_atom_list 

  !  Current tau value
  double precision, save :: current_tau  
  ! General REAF weights 
  double precision, save :: reaf_weights(2) 

  ! Item-specific weights
  integer, save:: reaf_sch
  integer, parameter:: ItemTypeTotal=8 !! this must be consistent w/ gti.F90 and gti_simulationConst.h 
  double precision, save :: read_item_weights(ItemTypeTotal, 2) 

contains

!*******************************************************************************
!
! Subroutine:  read_alloc_data
!
!*******************************************************************************
subroutine read_alloc_data(atm_cnt, num_ints, num_reals)
  use pmemd_lib_mod
  use file_io_dat_mod
  use parallel_dat_mod
  use mdin_ctrl_dat_mod
  use ti_mod
  implicit none

! Formal arguments:
  integer, intent(in)           :: atm_cnt

! num_ints and num_reals are used to return allocation counts. 
  integer, intent(inout)        :: num_ints, num_reals     

! Local variables:
  integer                       :: alloc_failed

  ! allocate memory for the atom mask lists:
  allocate(reaf_atom_list(2, atm_cnt), &
           stat = alloc_failed)
  if (alloc_failed .ne. 0) call setup_alloc_error

  write (mdout, '(a)') '|--------------------------------------------------------------------------------------------' 
  write (mdout, '(a)') '|Info for the REAF module:' 

  if (reaf_tau .lt. 0) then
    if (reaf_temp .gt. 0) then
       if (reaf_temp .lt. temp0) then
         write (mdout, '(a)') '|  The input REAF temp must be higher than temp0 '
       end if
       current_tau=1.0-sqrt(temp0/reaf_temp)
       write (mdout, '(a, f10.6)') '| The input REAF temp is ', reaf_temp
    else
       write (mdout, '(a)') '|  Both reaf_temp and reaf_tau are not set.  REAF will not run.'
       current_tau=-999.0
    endif
  else
    current_tau=reaf_tau
  endif    
  
  write (mdout, '(a, f10.6)') '| Current REAF tau value:', current_tau
  write (mdout, '(a, i4)') '| Current gti_add_re value:', gti_add_re

  num_ints = num_ints + size(reaf_atom_list) 

end subroutine


!*******************************************************************************
!
! Subroutine:  read_init_mask
!
!*******************************************************************************

subroutine read_init_mask(atm_cnt, nres, igraph, isymbl, res_atms, labres, &
                            crd, icfe, clambda, reaf_mask1, reaf_mask2)
  !!use mdin_ctrl_dat_mod
  use file_io_dat_mod
  use findmask_mod
  use pmemd_lib_mod 
  use ti_mod

  implicit none

! Formal arguments:
  integer, intent(in)             :: atm_cnt, nres
  integer, intent(in)             :: res_atms(nres)
  character(len=4), intent(in)    :: igraph(atm_cnt)
  character(len=4), intent(in)    :: isymbl(atm_cnt)
  character(len=4), intent(in)    :: labres(nres)
  double precision, intent(in)    :: crd(3 * atm_cnt)
  integer, intent(in)             :: icfe
  double precision, intent(in)    :: clambda
  character(len=*), intent(inout) :: reaf_mask1
  character(len=*), intent(inout) :: reaf_mask2

! Local variables:
  integer                         :: i, j, k
  character(len=4), parameter     :: ep_symbl = 'EP  '
  integer                         :: mask(atm_cnt, 2)
  integer                         :: alloc_failed

  double precision,parameter      :: delta=1e-12 

  mask(:, :) = 0
  call atommask(atm_cnt, nres, 0, igraph, isymbl, res_atms, labres, &
       crd, reaf_mask1, mask(:, 1))
  call atommask(atm_cnt, nres, 0, igraph, isymbl, res_atms, labres, &
       crd, reaf_mask2, mask(:, 2))

  reaf_atom_list(1,:) = mask(:, 1)
  reaf_atom_list(2,:) = mask(:, 2)

  if (current_tau<0) then
    write (mdout, '(a,a)') '| REAF is not running since the scaling' &
                            ,' factor is negative or not set.'
    reaf_mode=-1
  else if (current_tau<delta) then 
    write (mdout, '(a,a)') '| REAF is not running since the scaling' &
                            ,' factor is unity (reaf_tau=0).'
    reaf_mode=-1
 
  else 
    write (mdout, '(a,a,a,i7,a)') '| REAF Mask 1 ',trim(reaf_mask1),&
        ' matches', sum(reaf_atom_list(1,:)),' atoms'
    write (mdout, '(a,a,a,i7,a)') '| REAF Mask 2 ',trim(reaf_mask2),&
        ' matches', sum(reaf_atom_list(2,:)),' atoms'

    if ( sum(reaf_atom_list(1,:))+sum(reaf_atom_list(2,:)) .eq. 0) then
      write (mdout, '(a,a)') '| REAF is not running since the REAF' &
                             ,' region is not defined.'
      reaf_mode=-1
    else 
      if (icfe .eq. 0) then
        reaf_mode=0
        write (mdout, '(a)') '| REAF is running without TI' 
      else 
        !! REAF at region #1
        if (abs(clambda)<delta) then 
          reaf_mode=0
          write (mdout, '(a)') '| REAF is running at TI/Lambda=0 real state.' 
          if ( sum(reaf_atom_list(1,:)) .eq. 0) then
            write (mdout, '(a,a)') '| REAF will not run since the REAF' &
                                   ,' region #1 is not defined.'
            reaf_mode=-1
          else if  ( sum(reaf_atom_list(1,:)*ti_lst(2,:)) .gt. 0 ) then
            write (mdout, '(a,a,a)') '| Error! REAF region can only be' &
              ,' defined in the non-TI or the same TI region,' & 
              ,'e.g., reaf_mask1 cannot overlap with timask2'
            call mexit(mdout, 1)
          endif

        !! REAF at region #2
        else if (abs(clambda-1.0)<delta) then 
          reaf_mode=1
          write (mdout, '(a)') '| REAF is running at TI/Lambda=1 real state' 
          if ( sum(reaf_atom_list(2,:)) .eq. 0) then
            write (mdout, '(a,a)') '| REAF will not run since the REAF' &
                                   ,' region #2 is not defined.'
            reaf_mode=-1
          else if  ( sum(reaf_atom_list(2,:)*ti_lst(1,:)) .gt. 0 ) then
            write (mdout, '(a,a,a)') '| Error! REAF region can only be' &
              ,' defined in the non-TI or the same TI region,' & 
              ,'e.g., reaf_mask1 cannot overlap with timask2'
            call mexit(mdout, 1)
          endif
        else 
          reaf_mode=-1
          write (mdout, '(a,a)') '| Currently REAF only runs with TI real states' &
                 ,'(lambda=0,1).  REAF will be disabled now.' 
        end if
      end if
    endif
  endif 

  write (mdout, '(a)') '|--------------------------------------------------------------------------------------------' 

end subroutine read_init_mask


subroutine read_update_weights

  if (reaf_sch.ne.0) call gti_reaf_tau_schedule(tau_scheduling_name)
 
  call gti_update_tau(current_tau, reaf_weights, read_item_weights, reaf_sch.ne.0)

end subroutine read_update_weights

#endif

end module
