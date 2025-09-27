
!*******************************************************************************
!
! Module:  rmsd_mod
!
! Description: Module for Replica-Exchange with Arbitrary Degrees if freedom (rmsd) support. 
!              Written by Taisung Lee (Rutgers, 2020)
!              
!*******************************************************************************

module rmsd_mod

#ifdef GTI

  use file_io_dat_mod
  use state_info_mod
  implicit none
  
  logical, save :: initialized=.false.

  integer, save :: rmsd_mode = 0
  integer, save :: number_rmsd_set = 0
  integer, save :: rmsd_energy_type = 0
  integer, dimension(:),allocatable,save :: rmsd_atom_count 
  integer, dimension(:),allocatable,save :: rmsd_ti_region 
  integer, dimension(:,:),allocatable,save :: rmsd_atom_list 
  double precision, dimension(:,:,:),allocatable,save :: rmsd_ref_crd 
  double precision, dimension(:,:),allocatable,save :: rmsd_ref_com
  double precision, save :: rmsd_weights(5) 

contains

!*******************************************************************************
!
! Subroutine:  rmsd_init_mask
!
!*******************************************************************************

subroutine rmsd_init_data(atm_cnt, nres, igraph, isymbl, res_atms, labres, &
                           crd, ref_crd, icfe, &
             rmsd_mask, rmsd_strength, rmsd_ti, rmsd_type, num_ints, num_reals)
  !!use mdin_ctrl_dat_mod
  use file_io_dat_mod
  use findmask_mod
  use ti_mod
  use pmemd_lib_mod,   only : setup_alloc_error, mexit

  implicit none

! Formal arguments:
  integer, intent(in)             :: atm_cnt, nres
  integer, intent(in)             :: res_atms(nres)
  integer, intent(in)             :: rmsd_ti(5)
  integer, intent(in)             :: rmsd_type
  character(len=4), intent(in)    :: igraph(atm_cnt)
  character(len=4), intent(in)    :: isymbl(atm_cnt)
  character(len=4), intent(in)    :: labres(nres)
  double precision, intent(in)    :: crd(3 * atm_cnt), ref_crd(3*atm_cnt)
  integer, intent(in)             :: icfe
  character(len=*), intent(inout), dimension(5) :: rmsd_mask
  double precision, intent(inout), dimension(5) :: rmsd_strength 
  ! num_ints and num_reals are used to return allocation counts. 
  integer, intent(inout)        :: num_ints, num_reals    

! Local variables:
  integer                         :: i, j, k, m, n,  count(5), mask(atm_cnt), &
      temp_list(atm_cnt, 5), temp_region(5)
  integer                         :: alloc_failed

  rmsd_energy_type=rmsd_type
  number_rmsd_set=0
  count=0
  do i=1,5
    mask = 0
    temp_list(:,i)=0
    temp_region(i)=3
    call atommask(atm_cnt, nres, 0, igraph, isymbl, res_atms, labres, &
       crd, rmsd_mask(i), mask)
    
    if (icfe .ne. 0) then
      if (sum(ti_lst(2,:)*mask(:)).eq.0 .and. sum(ti_lst(3,:)*mask(:)).eq.0) then
        temp_region(i)=1
      elseif  (sum(ti_lst(1,:)*mask(:)).eq.0 .and. sum(ti_lst(3,:)*mask(:)).eq.0) then
        temp_region(i)=2
      endif
    endif
    
    count(i)=0
    do j=1, atm_cnt
       if(mask(j) .ne. 0) then
          count(i)=count(i)+1
          temp_list(count(i), i)=j
       end if
    end do
    if (count(i)>m) m=count(i)
    if (count(i)>0) number_rmsd_set=number_rmsd_set+1
  end do
  
  if (number_rmsd_set>0) then
     allocate(rmsd_atom_list(m, number_rmsd_set),  & 
        rmsd_atom_count(number_rmsd_set), &
        rmsd_ti_region(number_rmsd_set), &
        rmsd_ref_crd(3, m, number_rmsd_set), &
        rmsd_ref_com(3, number_rmsd_set), stat = alloc_failed)
     if (alloc_failed .ne. 0) then
       call setup_alloc_error
       call mexit(mdout,1)
     endif

    num_ints = num_ints + size(rmsd_atom_count) &
      + size(rmsd_ti_region) + size(rmsd_atom_list) + size(rmsd_ref_crd) &
      + size(rmsd_ref_com)
       
     k=0
     do i=1, 5
        if (count(i)>0) then
           k=k+1  
           rmsd_atom_count(k)=count(i)
           
           if (rmsd_ti(i)<0 .or. rmsd_ti(i)>3) then
             rmsd_ti_region(k)=temp_region(i)
           else
             rmsd_ti_region(k)=rmsd_ti(i)
           end if
           
           rmsd_atom_list(1:count(i),k)=temp_list(1:count(i),i)
           rmsd_weights(k) = rmsd_strength(i)
           do j=1, count(i)
              m=(rmsd_atom_list(j,k)-1)*3
              rmsd_ref_crd(1:3,j, k)=ref_crd(m+1:m+3)
           end do
           rmsd_ref_com(1,k)=sum(rmsd_ref_crd(1,1:count(i),k))/(count(i)*1.d0)
           rmsd_ref_com(2,k)=sum(rmsd_ref_crd(2,1:count(i),k))/(count(i)*1.d0)
           rmsd_ref_com(3,k)=sum(rmsd_ref_crd(3,1:count(i),k))/(count(i)*1.d0) 
           rmsd_ref_crd(1,1:count(i), k)=rmsd_ref_crd(1,1:count(i), k)-rmsd_ref_com(1,k)
           rmsd_ref_crd(2,1:count(i), k)=rmsd_ref_crd(2,1:count(i), k)-rmsd_ref_com(2,k)
           rmsd_ref_crd(3,1:count(i), k)=rmsd_ref_crd(3,1:count(i), k)-rmsd_ref_com(3,k)           
        end if
     end do   
     
     rmsd_mode=1
     write (mdout, '(a)') '|--------------------------------------------------------------------------------------------' 
     write (mdout, '(a)') '|Info for the RMSD module: ' 
     if (rmsd_energy_type .eq. 0) then
        write (mdout, '(a, i2, a)') '| RMSD Type : ', rmsd_energy_type, '    E=w*RMSD'
     else
        write (mdout, '(a, i4, a)') '| RMSD Type : ', rmsd_energy_type, ' E=w*Sum(distance deviation)'
     endif
    
     write (mdout, '(a, i4)') '| Number of RMSD sets: ', number_rmsd_set
     do i=1, number_rmsd_set
       if (icfe.eq.0 .or. rmsd_ti_region(i).eq.3 ) then
          write (mdout, '(a, i2 , a, i4, a, f8.3)') '|  Set #: ', i, ' has ', rmsd_atom_count(i), &
              ' atoms and a weight of ', rmsd_weights(i)
       else
          write (mdout, '(a, i2 , a, i4, a, f8.3, a, i2)') '|  Set #: ', i, ' has ', rmsd_atom_count(i), &
              ' atoms and a weight of ', rmsd_weights(i), ' and is in TI region #', rmsd_ti_region(i)         
       endif  
     end do  
     write (mdout, '(a)') '|--------------------------------------------------------------------------------------------' 
     
  endif

  call gti_init_rmsd(number_rmsd_set, rmsd_energy_type, rmsd_atom_count, &
    rmsd_ti_region, rmsd_atom_list, rmsd_ref_crd, rmsd_ref_com, rmsd_weights)
  
end subroutine rmsd_init_data

#endif

end module
