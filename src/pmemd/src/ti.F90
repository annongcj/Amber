
#include "copyright.i"

!*******************************************************************************
!
! Module:  ti_mod
!
! Description: Module for TI/MBAR including softcore support. 
!              Written by Joe Kaus (UCSD, 2012)
!              
!*******************************************************************************

module ti_mod
  use file_io_dat_mod
  use state_info_mod
  implicit none
  
  ! The following are broadcast at startup:
  ! Count of atoms
  integer, save                            :: ti_atm_cnt(2) !V0/V1
  integer, save                            :: ti_ti_atm_cnt(2) !Just TI atoms
  integer, save                            :: ti_numextra_pts(2) !Extra points
  ! Use ti_latm_cnt/ti_latm_lst for all linear scaling terms
  integer, save                            :: ti_latm_cnt(2) !Just linear atoms

  ! So we can just use this module instead of mdin_ctrl_dat_mod
  integer, save                            :: ifmbar_lcl
   
  ! Internally we use ti_mode rather than icfe/ifsc
  ! ti_mode == 1 is icfe = 1, ifsc = 0
  ! ti_mode == 2 is icfe = 1, ifsc = 1
  ! ti_mode == 3 is the same as 2, but the charges on all of the TI atoms are
  ! 0 so we can do the recip sum once giving a performance increase
  integer, save                            :: ti_mode
  integer, save                            :: klambda_lcl

  ! ti_region = 1/2 for V0/V1
  ! (ti_region+1, atm_id) -> 1 for TI atm, 0 for common atm 
  ! (3, atm_id) gives a list of common atoms (0 for TI, 1 for common atm)
  integer, dimension(:,:),allocatable,save :: ti_lst 
  ! (atm_id) gives 1 for sc atom 0 for linearly scaled atom
  integer, dimension(:),allocatable,save   :: ti_sc_lst

  ! (ti_region+1, ti_atm_id) -> list of TI atms
  ! (3, nonti_atm_id) -> list of non TI atoms
  integer, dimension(:,:),allocatable,save :: ti_atm_lst  

  ! (ti_region, ti_atm_id) -> list of TI atms to be scaled linearly
  ! (3, nonti_atm_id) -> list of non TI atoms
  integer, dimension(:,:),allocatable,save :: ti_latm_lst  

  ! (atm_id for TI region 1): 2 for sc_bond; 4 for sc_angle; 8 for sc_torsion
  ! (atm_id for TI region 2): 16 for sc_bond; 32 for sc_angle; 64 for sc_torsion
  integer, dimension(:),allocatable,save   :: ti_sc_bat_lst

  ! TI atoms where vdw being zero-ed
  integer, dimension(:),allocatable,save   :: ti_vdw_mask_list

#if defined(CUDA)
  ! reformatted ti_lst for passing to the C code
  integer, dimension(:),allocatable,save :: ti_lst_repacked

  ! reformatted ti_latm_lst for passing to the C code
  integer, dimension(:),allocatable,save :: ti_latm_lst_repacked

#endif

  ! VDW parameters
  double precision, dimension(:),allocatable,save :: ti_foureps
  double precision, dimension(:),allocatable,save :: ti_sigma6

  ! These are broadcast in ti_change_weights

  double precision, save :: ti_weights(2)
  double precision, save :: ti_dweights(2)

  integer, save::ti_lam_sch
  integer, parameter:: ItemTypeTotal=9 !! this must be consistent w/ gti.F90 and gti_simulationConst.h 
  !!integer, parameter::MaxNumberTIAtom=500 !! this must be consistent w/ gti_controlVariable.i 
  double precision, save :: ti_item_weights(ItemTypeTotal, 2) 
  double precision, save :: ti_item_dweights(ItemTypeTotal, 2)

  double precision, save :: ti_signs(2)
  double precision, save :: ti_klambda_factor

  ! The following are not broadcast:
  ! For the recip sum, mask the TI region given in ti_mask_piece
  integer, save                            :: ti_mask_piece
  ! For restraints, keep track of the type of atoms involved
  integer, save                            :: res_ti_region

  ! Split the nonbonded force arrays so we can properly 
  ! calculate the virial (if needed) (ti_region, 3, atm_cnt)
  ! NOTE: For performance reasons, it is best to avoid slicing
  ! these arrays -- ti_nb_frc(1,:,:) -- when passing to subroutines
  double precision, dimension(:,:,:),allocatable,save :: ti_nb_frc
  double precision, dimension(:,:,:),allocatable,save :: ti_img_frc

  ! Needed to calculate the virial
  double precision, save                 :: ti_net_frcs(2, 3)
  double precision, save                 :: ti_molvir_netfrc_corr(2, 3, 3)

  ! (ti_region,3,3) 
  double precision, dimension(:,:,:),allocatable,save :: ti_vir 

  integer, parameter    :: ti_vir0 = 1 ! (1) Elec Virial for V0
  integer, parameter    :: ti_vir1 = 2 ! (2) Elec Virial for V1
  integer, parameter    :: ti_ene0 = 3 ! (3) Elec Energy for V0
  integer, parameter    :: ti_ene1 = 4 ! (4) Elec Energy for V1     
  integer, parameter    :: ti_vir_ene_cnt = 4 ! update as needed
  double precision, save:: ti_vve(ti_vir_ene_cnt)

  ! Temporary variables:
  ! Marks the current TI region for the given set of atms (1/2/3 = V0/V1/Common)
  integer, save          :: ti_region
  logical, save          :: ti_lscale ! scale terms linearly when ifsc > 0
  double precision, save :: ti_pot(2)
  double precision, save :: ti_pot_gb(2, 4) ! For GB 

  ! Note: The full dvdl value is stored at ti_ene(1,si_dvdl) not here:
  ! This holds dvdl temporarily during each part of the calculation
  double precision       :: ti_dvdl

  ! Energy arrays
  ! TI energy arrays
  double precision,save :: ti_ene(2,si_cnt)   !(ti_region, ene)
  double precision,save :: ti_ene_ave(2,si_cnt)
  double precision,save :: ti_ene_rms(2,si_cnt)
  double precision,save :: ti_ene_tmp(2,si_cnt)
  double precision,save :: ti_ene_tmp2(2,si_cnt)
  double precision,save :: ti_ene_old(2,si_cnt)
  double precision,save :: ti_ene_old2(2,si_cnt)

  ! energy arrays for interactions between TI and other atoms
  double precision,save :: ti_others_ene(2,si_cnt)   !(ti_region, ene)
  double precision,save :: ti_others_ene_ave(2,si_cnt)
  double precision,save :: ti_others_ene_rms(2,si_cnt)
  double precision,save :: ti_others_ene_tmp(2,si_cnt)
  double precision,save :: ti_others_ene_tmp2(2,si_cnt)
  double precision,save :: ti_others_ene_old(2,si_cnt)
  double precision,save :: ti_others_ene_old2(2,si_cnt)

  ! Augmented energy arrays - for terms not found in the si array
  integer, parameter    :: ti_der_term = 1 ! (1) ADDITIONAL DERIVATIVE TERM
  integer, parameter    :: ti_elect_der_dvdl = 2! (2) DVDL of ELECTROSTATICS
  integer, parameter    :: ti_vdw_der_dvdl = 3   ! (3) DVDL of VDW
  integer, parameter    :: ti_rest_dist_ene = 4  ! (4) RESTRAINT_DIST
  integer, parameter    :: ti_rest_ang_ene = 5   ! (5) RESTRAINT_ANG
  integer, parameter    :: ti_rest_tor_ene = 6   ! (6) RESTRAINT_TORSIONAL
  !These are for the entire system, since they differ between V0/V1
  integer, parameter    :: ti_tot_kin_ene = 7    ! (7) TOTAL KINETIC ENERGY
  integer, parameter    :: ti_tot_temp = 8       ! (8) TOTAL TEMPERATURE
  integer, parameter    :: ti_tot_tot_ene = 9   !  (9) TOTAL ENERGY-All atms
  
  ! Lambda-dependent restraint (only on GPU)
  integer, parameter    :: ti_rest_der_term = 11  ! Total Lambda-depdent restraint contribution
  integer, parameter    :: ti_rest_dist_der_dvdl = 12  ! RESTRAINT_DIST
  integer, parameter    :: ti_rest_ang_der_dvdl = 13  ! RESTRAINT_ANG
  integer, parameter    :: ti_rest_tor_der_dvdl = 14   ! RESTRAINT_TORSIONAL 
  integer, parameter    :: ti_rest_rmsd_der_dvdl = 15   ! RESTRAINT_RMSD
  
   integer, parameter    :: ti_ene_aug_cnt = 16           ! update as needed

  double precision,save :: ti_ene_aug(2,ti_ene_aug_cnt)   !(ti_region, ene)
  double precision,save :: ti_ene_aug_ave(2,ti_ene_aug_cnt)
  double precision,save :: ti_ene_aug_rms(2,ti_ene_aug_cnt)
  double precision,save :: ti_ene_aug_tmp(2,ti_ene_aug_cnt)
  double precision,save :: ti_ene_aug_tmp2(2,ti_ene_aug_cnt)
  double precision,save :: ti_ene_aug_old(2,ti_ene_aug_cnt)
  double precision,save :: ti_ene_aug_old2(2,ti_ene_aug_cnt)

  ! DVDL Contributions by energy
  double precision,save :: ti_ene_delta(si_cnt)   !(ti_region, ene)
  double precision,save :: ti_ene_delta_ave(si_cnt)
  double precision,save :: ti_ene_delta_rms(si_cnt)
  double precision,save :: ti_ene_delta_tmp(si_cnt)
  double precision,save :: ti_ene_delta_tmp2(si_cnt)
  double precision,save :: ti_ene_delta_old(si_cnt)
  double precision,save :: ti_ene_delta_old2(si_cnt)

  ! Kinetic Energies
  integer, parameter    :: ti_eke = 1
  integer, parameter    :: ti_sc_eke = 2
  integer, parameter    :: ti_ekph = 3
  integer, parameter    :: ti_ekpbs = 4
  integer, parameter    :: ti_ekmh = 5
  integer, parameter    :: ti_solvent_kin_ene = 6
  integer, parameter    :: ti_solute_kin_ene = 7
  integer, parameter    :: ti_etot_save = 8
  integer, parameter    :: ti_ekin_cnt = 8
  double precision,save :: ti_kin_ene(2,ti_ekin_cnt)

  ! For com vel removal
  double precision,save :: ti_vcm(2, 3)
 
  ! Temperature factors
  double precision,save :: ti_fac(2, 3)
  double precision,save :: ti_factt(2)
  double precision,save :: ti_ekin0(2)
  double precision,save :: ti_scaltp(2)
  double precision,save :: ti_sc_fac(2)
  double precision      :: ti_temp

  ! Degrees of freedom
  double precision,save :: ti_rndf(2)
  double precision,save :: ti_rndfp(2)
  double precision,save :: ti_rndfs(2)
  integer,save          :: ti_ndfmin(2)
  integer,save          :: ti_dof_shaked(2)
  integer,save          :: ti_num_noshake(2)  
  integer,save          :: sc_num_noshake(2)

  ! For BAR
  integer, save                                     :: bar_values, bar_states
  integer                                           :: bar_i
  double precision, dimension(:)  ,allocatable,save :: bar_cont
  double precision, dimension(:,:,:),allocatable,save :: bar_lambda
  !only collect every bar_intervall steps
  logical, save                                     :: do_mbar
  
  ! Local heating 
  logical, save                                     :: do_localheating
  
  

  ! Total mass
  double precision, save                            :: ti_tmass(2)
  double precision, save                            :: ti_tmassinv(2)

  ! log DVDL
  double precision, dimension(:),  allocatable,save :: ti_dvdl_values
  integer, save                                     :: ti_dvdl_pnt

  ! Scale bonding SC terms by lambda, use SC interactions for all atoms
  integer, save                                     :: emil_sc_lcl
contains

!*******************************************************************************
!
! Subroutine:  ti_get_partner
!
! Description: Given an atom in the linearly scaled region return the
!              corresponding partner atom.
!              
!*******************************************************************************

subroutine ti_get_partner(iatm, partner_atm)
  use mdin_ctrl_dat_mod
  implicit none

! Formal arguments:
  integer, intent(in)           :: iatm
  integer, intent(out)          :: partner_atm     

! Local variables:
  integer                       :: i
  
  partner_atm = -1
  do i = 1, ti_latm_cnt(1)
    if (ti_latm_lst(1,i) .eq. iatm) then
      partner_atm = ti_latm_lst(2,i)
      return
    else if (ti_latm_lst(2,i) .eq. iatm) then
      partner_atm = ti_latm_lst(1,i)       
      return
    end if
  end do

end subroutine ti_get_partner

!*******************************************************************************
!
! Subroutine:  ti_alloc_mem
!
! Description: Allocates memory for SC arrays.
!              
!*******************************************************************************

subroutine ti_alloc_mem(atm_cnt, ntypes, num_ints, num_reals)
  use pmemd_lib_mod
  use file_io_dat_mod
  use parallel_dat_mod
  use mdin_ctrl_dat_mod
  implicit none

! Formal arguments:
  integer, intent(in)           :: atm_cnt, ntypes
! num_ints and num_reals are used to return allocation counts. Don't zero.
  integer, intent(inout)        :: num_ints, num_reals     

! Local variables:
  integer                       :: alloc_failed
  integer                       :: i
  integer                       :: check_num_mbar_states
  double precision              :: cur_lambda

  ! allocate memory for the atom mask lists:
  allocate(ti_lst(3, atm_cnt), &
           ti_sc_bat_lst(atm_cnt), &
           ti_vdw_mask_list(atm_cnt), &
           ti_sc_lst(atm_cnt), &
           ti_atm_lst(3, atm_cnt), &
           ti_latm_lst(3, atm_cnt), &
           ti_nb_frc(2, 3, atm_cnt), &
           ti_img_frc(2, 3, atm_cnt), &
           ti_vir(2, 3, 3), &
           ti_foureps(ntypes*ntypes), &
           ti_sigma6(ntypes*ntypes), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_ints = num_ints + size(ti_lst) + size(ti_sc_lst) + size(ti_atm_lst) + &
             size(ti_latm_lst)
  num_reals = num_reals + size(ti_foureps) + size(ti_sigma6) + &
              size(ti_nb_frc) + size(ti_img_frc) + size(ti_vir) 

  if (logdvdl .ne. 0) then
     if(irest .eq. 0) then
!the extra 1 is needed if we have irest .eq. 0 because we call pme_force an
!additional time. Could change this, since it may lead to double counting, but
!would require a very extensive refactor to allow us to still get the necessary
!energies/values for step 0
       allocate(ti_dvdl_values(nstlim+1), stat = alloc_failed)
     else
       allocate(ti_dvdl_values(nstlim), stat = alloc_failed)
     end if
     if (alloc_failed .ne. 0) call setup_alloc_error
     num_reals = num_reals + size(ti_dvdl_values)

     ti_dvdl_values(:) = 0.d0
     ti_dvdl_pnt = 0
  end if

  ! Save commonly used switches so other functions don't need to
  ! use mdin_ctrl_dat_mod

  ti_mode = 0 ! This is set in ti_init_dat
  klambda_lcl = klambda
  emil_sc_lcl = emil_sc
  
  ! Init data

  ! Atom/molecule data
  ti_lst(:,:) = 0
  ti_sc_lst(:) = 0
  ti_sc_bat_lst(:)=0
  ti_vdw_mask_list(:)=0
  ti_latm_cnt(:) = 0
  ti_atm_lst(:,:) = 0
  ti_latm_lst(:,:) = 0

  ti_tmass(:) = 0.d0
  ti_tmassinv(:) = 0.d0

  ti_atm_cnt(:) = 0
  ti_numextra_pts(:) = 0

  ! Energy terms
  ti_ene(:,:) = 0.d0
  ti_others_ene(:,:) = 0.d0
  ti_ene_aug(:,:) = 0.d0
  ti_ene_delta(:) = 0.d0
  ti_pot(:) = 0.d0
  ! Force/Virial
  ti_mask_piece = 0
  res_ti_region = 0
  ti_foureps(:) = 0.d0
  ti_sigma6(:) = 0.d0
  ti_nb_frc(:,:,:) = 0.d0
  ti_img_frc(:,:,:) = 0.d0
  ti_net_frcs(:,:) = 0.d0
  ti_vir(:,:,:) = 0.d0
  ti_vve(:) = 0.d0  

  ! DOF/Temperature factors
  ti_num_noshake(:) = 0  
  sc_num_noshake(:) = 0
  ti_fac(:,:) = 0.d0
  ti_factt(:) = 0.d0
  ti_ekin0(:) = 0.d0
  ti_scaltp(:) = 0.d0
  ti_rndf(:) = 0.d0
  ti_rndfp(:) = 0.d0
  ti_rndfs(:) = 0.d0
  ti_ndfmin(:) = 0

  ! For running averages
  ti_ene_ave(:,:) = 0.d0
  ti_ene_rms(:,:) = 0.d0
  ti_ene_tmp(:,:) = 0.d0
  ti_ene_tmp2(:,:) = 0.d0
  ti_ene_old(:,:) = 0.d0
  ti_ene_old2(:,:) = 0.d0

  ti_others_ene_ave(:,:) = 0.d0
  ti_others_ene_rms(:,:) = 0.d0
  ti_others_ene_tmp(:,:) = 0.d0
  ti_others_ene_tmp2(:,:) = 0.d0
  ti_others_ene_old(:,:) = 0.d0
  ti_others_ene_old2(:,:) = 0.d0

  ti_ene_aug_ave(:,:) = 0.d0
  ti_ene_aug_rms(:,:) = 0.d0
  ti_ene_aug_tmp(:,:) = 0.d0
  ti_ene_aug_tmp2(:,:) = 0.d0
  ti_ene_aug_old(:,:) = 0.d0
  ti_ene_aug_old2(:,:) = 0.d0

  ti_ene_delta(:) = 0.d0
  ti_ene_delta_rms(:) = 0.d0
  ti_ene_delta_tmp(:) = 0.d0
  ti_ene_delta_tmp2(:) = 0.d0
  ti_ene_delta_old(:) = 0.d0
  ti_ene_delta_old2(:) = 0.d0

  ! MBAR - Do init here
  ifmbar_lcl = 0
  do_mbar = .false.

  do_localheating=.false.
  if (gti_heating .gt. 0) do_localheating=.true.

  if (ifmbar .ne. 0) then
     ifmbar_lcl = ifmbar
!deprecated mbar input from AMBER 16
!     bar_states = 1 + nint( (bar_l_max-bar_l_min) / bar_l_incr )
     bar_states = mbar_states
     bar_values = nstlim / bar_intervall

     allocate(bar_lambda(2,bar_states,7), &
              bar_cont(bar_states), &
              stat = alloc_failed)     

     if (alloc_failed .ne. 0) call setup_alloc_error  
     num_reals = num_reals + size(bar_lambda) + size(bar_cont)

!this will error out if mbar_states mismatches the number of values in
!mbar_lambda. have to do it twice or we get memory errors when we try to assign
!values to bar_lambda that haven't had space allocated when
!check_num_mbar_states > mbar_states

       check_num_mbar_states = 0
       do i = 1, size(mbar_lambda)
         if (mbar_lambda(i) .le. 1.0) then
           if (check_num_mbar_states .eq. mbar_states) then
             write (mdout,'(a,/,a)') '     ERROR: the number of mbar lambdas &
               &specified in mbar_lambda is greater than mbar_states'
             call mexit(mdout, 1)
           end if
           if(gti_lam_sch .eq. 0) then
               bar_lambda(1,i,:) = (1.d0 - mbar_lambda(i)) ** klambda
               bar_lambda(2,i,:)= mbar_lambda(i)
           else
               bar_lambda(1,i,4:6) = 6.d0*((1.d0-mbar_lambda(i))**5) - &
                                   15.d0*((1.d0-mbar_lambda(i))**4) + &
                                   10.d0*((1.d0-mbar_lambda(i))**3)
               bar_lambda(2,i,4:6) = 6.d0*(mbar_lambda(i)**5) - &
                                   15.d0*(mbar_lambda(i)**4) + &
                                   10.d0*(mbar_lambda(i)**3)
               bar_lambda(1,i,1:3)= 1.d0 - mbar_lambda(i)
               bar_lambda(2,i,1:3)= mbar_lambda(i)
               bar_lambda(1,i,7)= 1.d0 - mbar_lambda(i)
               bar_lambda(2,i,7)= mbar_lambda(i)
           end if
           check_num_mbar_states = check_num_mbar_states + 1
         else
           EXIT
         end if
       end do
 
       if (check_num_mbar_states .lt. mbar_states) then
         write (mdout,'(a,/,a)') '     ERROR: the number of mbar lambdas &
           &specified in mbar_lambda is less than mbar_states'
         call mexit(mdout, 1)
       end if

       bar_cont(:) = 0.d0

  end if

  return

end subroutine ti_alloc_mem

!*******************************************************************************
!
! Subroutine:  ti_init_mask
!
! Description: Fill in the ti_lst array, based on the amber mask.
!
!*******************************************************************************

subroutine ti_init_mask(atm_cnt, nres, igraph, isymbl, res_atms, labres, &
                            crd, maskstr, maskstr2, scmaskstr1, scmaskstr2, &
                            tivdwstr,  &
                            sc_bondstr1, sc_anglestr1, sc_torsionstr1, &
                            sc_bondstr2, sc_anglestr2, sc_torsionstr2)
  use file_io_dat_mod
  use findmask_mod
#if defined(CUDA)
  use pmemd_lib_mod !I think we just need setup_alloc_error
#endif
  use mdin_ctrl_dat_mod, only : usemidpoint
  implicit none

! Formal arguments:
  integer, intent(in)             :: atm_cnt, nres
  integer, intent(in)             :: res_atms(nres)
  character(len=4), intent(in)    :: igraph(atm_cnt)
  character(len=4), intent(in)    :: isymbl(atm_cnt)
  character(len=4), intent(in)    :: labres(nres)
  character(len=*), intent(inout) :: maskstr
  character(len=*), intent(inout) :: maskstr2
  character(len=*), intent(inout) :: tivdwstr
  character(len=*), intent(inout) :: scmaskstr1
  character(len=*), intent(inout) :: scmaskstr2
  character(len=*), intent(inout) :: sc_bondstr1
  character(len=*), intent(inout) :: sc_anglestr1
  character(len=*), intent(inout) :: sc_torsionstr1
  character(len=*), intent(inout) :: sc_bondstr2
  character(len=*), intent(inout) :: sc_anglestr2
  character(len=*), intent(inout) :: sc_torsionstr2

  double precision, intent(in)    :: crd(3 * atm_cnt)

! Local variables:
  integer                         :: i, j, k
  character(len=4), parameter     :: ep_symbl = 'EP  '
  integer                         :: mask(atm_cnt), mask2(atm_cnt)
#if defined(CUDA)
  integer                         :: alloc_failed
#endif

  if( emil_sc_lcl .ge. 1 ) then
       mask(:)  = 1
       mask2(:) = 0
  else
if(.not. usemidpoint) then
       mask(:) = 0
       call atommask(atm_cnt, nres, 0, igraph, isymbl, res_atms, labres, &
                crd, maskstr, mask)
       mask2(:) = 0
       call atommask(atm_cnt, nres, 0, igraph, isymbl, &
                res_atms, labres, &
                crd, maskstr2, mask2)
endif
  end if

  ti_lst(1,:) = mask(:)
  ti_lst(2,:) = mask2(:)
  ti_lst(3,:) = 1 - (ti_lst(1,:)+ti_lst(2,:)) !mask for common atoms

  ti_ti_atm_cnt(1) = sum(ti_lst(1,:))
  ti_ti_atm_cnt(2) = sum(ti_lst(2,:))

  ! Fill in list of TI atoms
  do j = 1, 3
     k = 1
     do i = 1, atm_cnt
        if( ti_lst(j,i) .ne. 0) then
           ti_atm_lst(j,k) = i
           k = k + 1
        end if
     end do
  end do

  ! Count extra points in each TI region
  ! How is this done before init_extra_pnts_nb14 ???
  do i = 1, atm_cnt
    if (isymbl(i) .eq. ep_symbl) then
       ti_numextra_pts(1) = ti_numextra_pts(1) + ti_lst(1,i) + ti_lst(3,i)
       ti_numextra_pts(2) = ti_numextra_pts(2) + ti_lst(2,i) + ti_lst(3,i)
    end if
  end do

  ti_atm_cnt(1) = sum(ti_lst(3,:)) + ti_ti_atm_cnt(1)
  ti_atm_cnt(2) = sum(ti_lst(3,:)) + ti_ti_atm_cnt(2)

  if( emil_sc_lcl .eq. 0 ) then
     write (mdout, '(a,a,a,i7,a)') '     TI Mask 1 ',trim(maskstr),&
         ' matches', ti_ti_atm_cnt(1),' atoms'
     write (mdout, '(a,a,a,i7,a)') '     TI Mask 2 ',trim(maskstr2),&
         ' matches', ti_ti_atm_cnt(2),' atoms'
     write (mdout, '(a,i7,a)') '     TI region 1: ',ti_atm_cnt(1),&
         ' atoms'
     write (mdout, '(a,i7,a)') '     TI region 2: ',ti_atm_cnt(2),&
         ' atoms'
  end if

  ! Determine zero-vdw list
  ti_vdw_mask_list(:)=0

  if (tivdwstr .ne. '') then
    write (mdout, '(a,a,a,i7,a)') '    TI VDW MASK ',trim(tivdwstr),&
           ' matches', ti_vdw_mask_list(mask),' atoms'
    call atommask(atm_cnt, nres, 0, igraph, isymbl, res_atms, labres, &
              crd, tivdwstr, ti_vdw_mask_list)
  endif

  ! Determine sc atom list if present
  if (scmaskstr1 .ne. '' .or. (emil_sc_lcl .ge. 1) ) then
     if( emil_sc_lcl .ge. 1 ) then
        mask(:) = 1
     end if
if(.not. usemidpoint) then
     if( emil_sc_lcl .lt. 1 ) then
        mask(:) = 0
        call atommask(atm_cnt, nres, 0, igraph, isymbl, res_atms, labres, &
                  crd, scmaskstr1, mask)
     end if
endif
     ti_sc_lst(:) = mask(:)
     ti_latm_cnt(1) = ti_ti_atm_cnt(1) - sum(mask)
     write (mdout, '(a,a,a,i7,a)') '     SC Mask 1 ',trim(scmaskstr1),&
           ' matches', sum(mask),' atoms'

     if(sc_bondstr1.ne. '') then
        mask(:) = 0
        call atommask(atm_cnt, nres, 0, igraph, isymbl, res_atms, labres, &
                  crd, sc_bondstr1, mask)
        ti_sc_bat_lst(:) = ti_sc_bat_lst(:) + 2*mask(:)
     endif
     if(sc_anglestr1.ne. '') then
        mask(:) = 0
        call atommask(atm_cnt, nres, 0, igraph, isymbl, res_atms, labres, &
                  crd, sc_anglestr1, mask)
        ti_sc_bat_lst(:) = ti_sc_bat_lst(:)+ 4*mask(:)
     endif
     if(sc_torsionstr1.ne. '') then
        mask(:) = 0
        call atommask(atm_cnt, nres, 0, igraph, isymbl, res_atms, labres, &
                  crd, sc_torsionstr1, mask)
        ti_sc_bat_lst(:) = ti_sc_bat_lst(:)+ 8*mask(:)
     endif

  else
    ti_latm_cnt(1) = ti_ti_atm_cnt(1)
  end if

  if ( scmaskstr2 .ne. '') then
if(.not. usemidpoint) then
    mask(:) = 0
    call atommask(atm_cnt, nres, 0, igraph, isymbl, res_atms, labres, &
                  crd, scmaskstr2, mask)
    ti_sc_lst(:) = ti_sc_lst(:) + mask(:)
    ti_latm_cnt(2) = ti_ti_atm_cnt(2) - sum(mask)
    write (mdout, '(a,a,a,i7,a)') '     SC Mask 2 ',trim(scmaskstr2),&
           ' matches', sum(mask),' atoms'
endif
     if(sc_bondstr2.ne. '') then
        mask(:) = 0
        call atommask(atm_cnt, nres, 0, igraph, isymbl, res_atms, labres, &
                  crd, sc_bondstr2, mask)
        ti_sc_bat_lst(:) = ti_sc_bat_lst(:) + 16*mask(:)
     endif
     if(sc_anglestr2.ne. '') then
        mask(:) = 0
        call atommask(atm_cnt, nres, 0, igraph, isymbl, res_atms, labres, &
                  crd, sc_anglestr2, mask)
        ti_sc_bat_lst(:) = ti_sc_bat_lst(:)+ 32*mask(:)
     endif
     if(sc_torsionstr2.ne. '') then
        mask(:) = 0
        call atommask(atm_cnt, nres, 0, igraph, isymbl, res_atms, labres, &
                  crd, sc_torsionstr2, mask)
        ti_sc_bat_lst(:) = ti_sc_bat_lst(:)+ 64*mask(:)
     endif
  else
    ti_latm_cnt(2) = ti_ti_atm_cnt(2)
  end if

  ! Fill in list of TI atoms
  do j = 1, 3
     k = 1
     do i = 1, atm_cnt
        if( ti_lst(j,i) .ne. 0 .and. ti_sc_lst(i) .eq. 0) then
           ti_latm_lst(j,k) = i
           k = k + 1
        end if
     end do
  end do

  ! Setup integer flags for type of softcoring:
  ! 0,2 ==> softcore interactions to atoms with opposite value
  ! 1   ==> softcore interactions to atoms with same value
  ! i.e. if( ti_sc_lst(i) + ti_sc_lst(j) .eq. 2 ) then softcore. 
  if( emil_sc_lcl .eq. 0 ) then
      ti_sc_lst(:) = ti_sc_lst(:) * 2
  end if

  return

end subroutine ti_init_mask

!*******************************************************************************
!
! Subroutine:  ti_init_dat
!
! Description: Init data for TI.
!
!*******************************************************************************

subroutine ti_init_dat(atm_cnt, crd, icfe, ifsc, qterm)

  use file_io_dat_mod
  use parallel_dat_mod
  use pmemd_lib_mod
  use mdin_ctrl_dat_mod, only : gti_lam_sch, mbar_lambda
  use charmm_mod, only : charmm_active
  use mdin_ctrl_dat_mod, only :igamd

  implicit none

! Formal arguments:
  integer, intent(in)              :: atm_cnt
  double precision, intent(inout)  :: crd(3, atm_cnt)
  integer, intent(in)              :: icfe
  integer, intent(in)              :: ifsc
  double precision, intent(inout)  :: qterm(atm_cnt)

! Local variables:
  integer                       :: i
  integer                       :: zero_cnt
  integer                       :: net_charge
  integer                       :: tot_ti_atm_cnt

  call ti_calc_parameters
  call ti_zero_arrays

  zero_cnt = 0

  if (icfe .ne. 0 .and. ifsc .eq. 0) then
    ti_mode = 1
  else if (icfe .ne. 0 .and. ifsc .eq. 1) then
    ti_mode = 2
  end if

  ! Now check the charges on the TI atms. If they are all 0.d0, then we can
  ! use ti_mode = 3, which only does the recip calculation once.
  if (ti_mode .eq. 2) then
     net_charge = 0
     do i = 1, atm_cnt     
        if((ti_lst(1,i) + ti_lst(2,i)) .ne. 0) then
           if(qterm(i) .eq. 0.d0) then
              net_charge = net_charge + 1
           end if
        end if
     end do

     tot_ti_atm_cnt = ti_ti_atm_cnt(1) + ti_ti_atm_cnt(2)
     if (net_charge .eq. tot_ti_atm_cnt) then
        ti_mode = 3
         write (mdout, '(a)') '    No charge on TI atoms. &
                              &Skipping extra recip sum.'
     end if
  end if

  if (ifmbar_lcl .ne. 0) then
     write (mdout, *)
     write (mdout,'(a)') '    MBAR - lambda values considered:'
     if(gti_lam_sch .eq. 1) then
     write (mdout,'(a,i4,a,20(1x,f6.4))') '    ', bar_states, ' total: ', &
             bar_lambda(2,1:bar_states,1)
     else
         write (mdout,'(a,i4,a,20(1x,f6.4))') '    ', bar_states, ' total: ', &
             mbar_lambda(1:bar_states)
     end if
     write (mdout,'(a,i6,a)')      '    Extra energies will be computed ',&
             bar_values, ' times.'
  end if

  ! Additional input checks
  if (charmm_active) then
     if(igamd.eq.6.or.igamd.eq.7.or.igamd.eq.8.or.igamd.eq.9.or.igamd.eq.10.or.igamd.eq.11.or. &
             igamd.eq.12.or.igamd.eq.13.or.igamd.eq.16.or.igamd.eq.17.or.igamd.eq.18.or.igamd.eq.19.or.&
             igamd.eq.20.or.igamd.eq.21.or.igamd.eq.23.or.igamd.eq.26.or.igamd.eq.28) then
       write (mdout,'(a)') '   WARING: TI is compatible with CHARMM in testing with igamd = 6-13,16-28'
     else
     write (mdout,'(a)') '    ERROR: TI is incompatible with CHARMM for now.'
     call mexit(mdout, 1)
     endif
  end if

  if (ti_mode .eq. 1) then
     ! Extra input checks for icfe = 1, ifsc = 0
     write (mdout,'(a)') '     Checking for mismatched coordinates.'
     if (ti_ti_atm_cnt(1) .eq. 0) then
        write (mdout,'(a,/,a)') '     ERROR: timask1/2 must match at least &
           &one atom for non-softcore run'
        call mexit(mdout, 1)
     end if
     if (ti_ti_atm_cnt(1) .ne. ti_ti_atm_cnt(2)) then
        write (mdout,'(a,/,a)') '     ERROR: timask1/2 must match the same &
           &number of atoms for non-softcore run'
       call mexit(mdout, 1)
     end if
     ti_latm_cnt(:) = ti_ti_atm_cnt(:)
     ti_latm_lst(:,:) = ti_atm_lst(:,:)
     call ti_check_and_adjust_crd(atm_cnt, crd)          
  else
    if (ti_latm_cnt(1) .ne. ti_latm_cnt(2)) then
        write (mdout,'(a,/,a)') '     ERROR: The number of linearly scaled &
           &atoms must be the same, check timask1/2 and scmask1/2'
       call mexit(mdout, 1)
    end if
    call ti_check_and_adjust_crd(atm_cnt, crd)  
  end if

  return

end subroutine ti_init_dat

#if defined(CUDA)
!*******************************************************************************
!
! Subroutine:  ti_alloc_gpu_array
!
! Description: c can't take a 3xn array as a passed argument, so we have to
! make a flipped array. I think we still need ti_lst and ti_latm_lst or we 
! could just make a !different array in ti_init_mask
! we need ti_lst for just generally figuring out the regions of each ti atom
! and we need ti_latm_lst for matching up region 1 atoms with region 2 atoms
! for vector exchange
!
!*******************************************************************************
subroutine ti_alloc_gpu_array(atm_cnt, linear_atm_cnt)
  use pmemd_lib_mod !just setup_alloc_error?

  implicit none

! Formal arguments:
  integer, intent(inout)          :: atm_cnt
  integer, intent(inout)          :: linear_atm_cnt

! Local variables:
  integer                         :: alloc_failed
  integer                         :: j,k


    allocate(ti_lst_repacked(3 * atm_cnt + 1), &
             ti_latm_lst_repacked(2 * linear_atm_cnt), &
             stat = alloc_failed)
    ti_lst_repacked(:) = 0
    ti_latm_lst_repacked(:) = 0
    do j = 1, linear_atm_cnt
      do k = 1, 2
        ti_latm_lst_repacked(j + (k - 1) * linear_atm_cnt) = ti_latm_lst(k,j)
      end do
    end do
!try to add the region 1, then region 2, and then common atoms individually
  do j = 1, atm_cnt
    do k = 1, 3
      ti_lst_repacked(j + atm_cnt * (k - 1)) = ti_lst(k,j)
    end do
  end do

  if (alloc_failed .ne. 0) call setup_alloc_error  
end subroutine
#endif

!*******************************************************************************
!
! Subroutine:  ti_check_and_adjust_crd
!
! Description: Check for mismatched crds and adjust for small deviations. 
!
!*******************************************************************************

subroutine ti_check_and_adjust_crd(atm_cnt, crd)

  use file_io_dat_mod
  use parallel_dat_mod
  use pmemd_lib_mod
#ifdef GTI
  use mdin_ctrl_dat_mod, only: gti_syn_mass, ntc, tishake
  use prmtop_dat_mod, only: atm_mass
#endif

  implicit none

! Formal arguments:
  integer, intent(in)             :: atm_cnt
  double precision, intent(inout) :: crd(3, atm_cnt)

! Local variables:
  double precision                :: crd_diff
  integer                         :: i, j, k
  integer                         :: m
  integer                         :: nadj
  integer                         :: atm_i
  integer                         :: atm_j
  integer                         :: ti_latm_lst_sort(atm_cnt)
  logical                         :: picked(atm_cnt)
  logical                         :: mismatched
  integer                         :: total_idx

  nadj = 0
  crd_diff = 0.d0

  !reorder atoms so that if linear atoms are out of order in the prmtop
  !the simulation will still run 
  ti_latm_lst_sort(:) = 0
  picked(:)=.false.
  do k = 1,20
  do i = 1, ti_latm_cnt(1)
      if (ti_latm_lst_sort(i) .eq. 0) then
    do j = 1, ti_latm_cnt(1)
      atm_i = ti_latm_lst(1,i)
      atm_j = ti_latm_lst(2,j)
        if (.not. picked(atm_j)) then
          crd_diff=0.0
      do m = 1, 3
            crd_diff = crd_diff+ (crd(m, atm_i) - crd(m, atm_j))**2
          end do
          crd_diff = sqrt(crd_diff)
          if (crd_diff .lt. (0.1d0*k) ) then  
              ti_latm_lst_sort(i) = atm_j
              picked(atm_j)=.true.
              exit
          end if
        end if
      end do
      end if
    end do
  end do

  ti_latm_lst(2,:) = ti_latm_lst_sort(:)

  do i = 1, ti_latm_cnt(1)
    atm_i = ti_latm_lst(1,i)
    atm_j = ti_latm_lst(2,i)
    if (atm_j .eq. 0) then
      write (mdout,'(a,i7,a)') '     Error : Atom ', &
             atm_i,' does not have match in V1 !'
      call mexit(mdout, 1)
    end if
  end do

#ifdef GTI
  k=0
  if(gti_syn_mass<0) k=1 ! no user-defined value is found
  mismatched=.false.
#endif

  do i = 1, ti_latm_cnt(1)
     atm_i = ti_latm_lst(1,i)
     atm_j = ti_latm_lst(2,i)

#ifdef GTI
    if (abs(atm_mass(atm_i)-atm_mass(atm_j))>0.01) then
      mismatched=.true.
      write(mdout, '(a, i4, f10.4, i4, f10.4)') &
        '| mismatched mass: atom #1, mass #1, atom #2, mass #2' &
        , atm_i, atm_mass(atm_i), atm_j, atm_mass(atm_j)
    endif
#endif

    crd_diff=0.0
     do m = 1, 3
      crd_diff = crd_diff+ (crd(m, atm_i) - crd(m, atm_j))**2
    end do
    if (crd_diff .gt. 1e-6) then
           if (crd_diff .gt. 0.1d0) then
              write (mdout,'(a,i7,a,i7,a)') '     WARNING: Local coordinate ', &
                   atm_i,' differs from partner coordinate ', atm_j,' !'
              write (mdout,'(a)') &
                 '     Atom coordinate disagreement, check input files.'
              !call mexit(mdout, 1)
           else
              nadj = nadj + 1
              if (nadj .lt. 11) then
                 if (nadj .lt. 10) then
                    write (mdout,'(a,i7,a,i7,a)') &
                       '     WARNING: Local coordinate ', &
                       atm_i,' differs from partner coordinate ', atm_j,' !'
                    write (mdout,'(a)') &
                       '     Deviation is small, changing partner coordinate.'
                 else
                    write (mdout,'(a)') &
                       '     ... making more adjustments ...'
                 end if
              end if
           end if
        end if
     end do

#ifdef GTI
  if (ntc .eq. 2  .and. tishake .eq. 0) then 
    k=3
    write(mdout, '(a)') &
    '| SHAKE is turn on and tishake is 0. gti_syn_mass will set to 3 if not explicitly specified.' 
  else if (mismatched) then
    k=0
  endif 

  if(gti_syn_mass<0) then 
    gti_syn_mass=k
    write (mdout,'(a, i4)') &
          '| gti_syn_mass has been set to',gti_syn_mass 
  end if
#endif

  if (nadj .gt. 9) then
     write (mdout,'(a,i7,a)') '     A total of ', nadj, &
        ' small coordinate adjustments were made, check results carefully.'
  end if

  return

end subroutine ti_check_and_adjust_crd

!*******************************************************************************
!
! Subroutine:  ti_calc_parameters
!
! Description: Calculate VDW parameters. 
!
!*******************************************************************************

subroutine ti_calc_parameters

  use prmtop_dat_mod

  implicit none

! Local variables:
  integer             :: i

  do i=1,((ntypes+1)*ntypes)/2
     if (gbl_cn1(i) .ne. 0.0d0) then ! catch zero vdW hydrogens
        ti_sigma6(i) = gbl_cn2(i) / gbl_cn1(i)
        ti_foureps(i) = gbl_cn2(i) * ti_sigma6(i)
     else
        ti_sigma6(i) = 1.0d0 ! will always be multiplied by zero
        ti_foureps(i) = 0.0d0
     end if
  end do

  return

end subroutine ti_calc_parameters

#if defined(MPI)
!*******************************************************************************
!
! Subroutine:  ti_bcast_dat
!
! Description: Broadcast TI data to all nodes. The weights are done seperately.
!              
!*******************************************************************************

subroutine ti_bcast_dat(atm_cnt, ntypes)
  use parallel_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:
  integer, intent(in)           :: atm_cnt, ntypes
! Local variables:
  integer                       :: num_ints, num_reals !returned value discarded
  integer                       :: alloc_failed

  if (.not. master) then
    num_ints = 0
    num_reals = 0
    call ti_alloc_mem(atm_cnt, ntypes, num_ints, num_reals)
    call ti_zero_arrays
  end if

  call mpi_bcast(ti_ti_atm_cnt, 2, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(ti_lst, 3 * atm_cnt, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(ti_sc_lst, atm_cnt, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(ti_latm_cnt, 2, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(ti_atm_lst, 3 * atm_cnt, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(ti_latm_lst, 3 * atm_cnt, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(ti_mode, 1, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(klambda_lcl, 1, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(ti_foureps, ntypes*ntypes, mpi_double_precision, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(ti_sigma6, ntypes*ntypes, mpi_double_precision, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(ti_atm_cnt, 2, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(ti_numextra_pts, 2, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  return

end subroutine ti_bcast_dat
#endif

!*******************************************************************************
!
! Subroutine:  ti_crgmask
!
! Description: Zero partial charges listed in maskstr. 
!
!*******************************************************************************

subroutine ti_crgmask(atm_cnt, nres, igraph, isymbl, res_atms, labres, crd, &
                      maskstr, qterm)

  use findmask_mod
  use file_io_dat_mod
  use mdin_ctrl_dat_mod, only : usemidpoint
  implicit none

! Formal arguments:

  integer, intent(in)             :: atm_cnt
  integer, intent(in)             :: nres
  integer, intent(in)             :: res_atms(nres)
  character(len=4), intent(in)    :: igraph(atm_cnt)
  character(len=4), intent(in)    :: isymbl(atm_cnt)
  character(len=4), intent(in)    :: labres(nres)
  character(len=*), intent(inout) :: maskstr
  double precision, intent(in)    :: crd(3 * atm_cnt)
  double precision, intent(inout) :: qterm(atm_cnt)

! Local variables:
  integer                         :: qmask(atm_cnt)
  integer                         :: i
  double precision                :: charge_removed

  if(len_trim(maskstr) .eq. 0) return

  charge_removed = 0.d0
if(.not. usemidpoint) &
  call atommask(atm_cnt, nres, 0, igraph, isymbl, res_atms, labres, crd, &
                maskstr, qmask)

  do i = 1, atm_cnt
     if(qmask(i) .eq. 1) then
        charge_removed = charge_removed + qterm(i) / 18.2223d0 
        write (mdout, '(a,f12.4,a,i7)') 'Removing charge of ', &
               qterm(i) / 18.2223d0,' from atom ',i
        qterm(i) = 0.0d0
     end if
  end do

  write (mdout, '(a,f12.8,a,i7,a)') 'Total charge of ', charge_removed, &
         ' removed from ', sum(qmask), ' atoms'

  return

end subroutine ti_crgmask

!*******************************************************************************
!
! Subroutine:  ti_check_neutral
!
! Description: Check a softcore system for neutrality and force neutrality if
!   we are below the cutoff. From check_neutral in pme_setup.fpp with 
!   modifications. These changes are needed to neutralize both V0/V1 while
!   keeping the common atoms syncronized.
!
!*******************************************************************************

subroutine ti_check_neutral
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none

! Local variables:
  integer            :: i, ib
  integer            :: j, jb
  integer            :: nsum(3)
  integer            :: nsum_net(2)
  integer            :: nchange(2)
  integer            :: ichange
  double precision   :: sum_val(3)
  double precision   :: sum_val_net(3)
  integer            :: charge_lst(3, natom)
  logical            :: neutralize

  neutralize =.true.
  if (gti_chg_keep>0) neutralize = .false.

  nsum(:) = 0
  sum_val(:) = 0.d0
  nchange(:) = 0
  ichange = 1

  if (ti_mode .eq. 1) then ! Non-softcore simulation
     charge_lst(:,:) = 0
     do j = 1, ti_ti_atm_cnt(1)
        ib = ti_atm_lst(1,j)
        jb = ti_atm_lst(2,j)
        if (abs((atm_qterm(ib)-atm_qterm(jb))/18.2223d0) > 1.d-10) then
           charge_lst(1,ib) = 1
           charge_lst(2,jb) = 1
        end if
     end do
     charge_lst(3,:) = 1 - (charge_lst(1,:) + charge_lst(2,:))
  else
     charge_lst(:,:) = ti_lst(:,:)
  end if

  do j = 1, natom
     if (charge_lst(1,j) .eq. 1) then
        sum_val(1) = sum_val(1) + atm_qterm(j)
        nsum(1) = nsum(1) + 1
     else if(charge_lst(2,j) .eq. 1) then
        sum_val(2) = sum_val(2) + atm_qterm(j)
        nsum(2) = nsum(2) + 1
     else
        ! So we don't double count common atoms with mode 1
        if (ti_mode .ne. 1 .or. ti_lst(2,j) .ne. 1) then
          sum_val(3) = sum_val(3) + atm_qterm(j)
          nsum(3) = nsum(3) + 1
        end if
     end if
  end do

  do i = 1, 2
     sum_val_net(i) = sum_val(i) + sum_val(3)
     nsum_net(i) = nsum(i) + nsum(3)
     if (master) write(mdout, '(/,5x,a,i2,a,f12.8)') &
       'Sum of charges for TI region ', i, ' = ', sum_val_net(i) / 18.2223d0       
     
     if (neutralize) then
       if (abs(sum_val_net(i)/18.2223d0) .gt. 0.01) then
          if (master) &
            write(mdout, '(5x,a,/)') 'Assuming uniform neutralizing plasma'
       else
          if (master) &
            write(mdout, '(5x,a,/)') 'Forcing neutrality...'
          nchange(i) = 1
          ichange = i
       end if
     else
       write(mdout, '(5x,a,/)') 'Skip neutralizing charges...'
     end if
  end do

  if (.not. neutralize) return

  if (nchange(1)+nchange(2) .eq. 1) then
     if (nsum_net(ichange) .gt. 0) then
        sum_val_net(ichange) = sum_val_net(ichange) / nsum_net(ichange)
     else
        sum_val_net(ichange) = 0.d0
     end if

     do j = 1, natom
        if (charge_lst(ichange,j) .eq. 1 .or. charge_lst(3,j) .eq. 1) then
           atm_qterm(j) = atm_qterm(j) - sum_val_net(ichange)
        end if
     end do  
  else if (nchange(1)+nchange(2) .eq. 2) then
     ! We need to remove the net charge on two effective simulation regions
     ! We first remove the average charge, then remove the remaining amount 
     ! from the softcore atoms. 
     
     ! Check for regions that have no charge
     if (nsum_net(1) .ne. 0 .and. nsum_net(2) .ne. 0) then
        sum_val_net(3) = 0.5d0 * (sum_val_net(1) / nsum_net(1) + &
                         sum_val_net(2) / nsum_net(2))
     else if (nsum_net(1) .ne. 0) then
        sum_val_net(3) = sum_val_net(1) / nsum_net(1)
     else if (nsum_net(2) .ne. 0) then
        sum_val_net(3) = sum_val_net(2) / nsum_net(2)
     else
        sum_val_net(3) = 0.d0
     end if

     ! For completely VDW runs/perturb to nothing runs these can be zero
     if (nsum(1) .gt. 0) then
        sum_val_net(1) = (sum_val_net(1) - sum_val_net(3) * nsum(3)) / nsum(1)
     else
        sum_val_net(1) = 0.d0
     end if

     if (nsum(2) .gt. 0) then
        sum_val_net(2) = (sum_val_net(2) - sum_val_net(3) * nsum(3)) / nsum(2)
     else
        sum_val_net(2) = 0.d0
     end if

     do j = 1, natom
        do i = 1, 3
           if(charge_lst(i,j) .eq. 1) atm_qterm(j) = atm_qterm(j) - sum_val_net(i)
        end do
     end do         
  end if

  return

end subroutine ti_check_neutral


!*******************************************************************************
!
! Subroutine:  ti_change_weights
!
! Description: Change and broadcast weights for TI.
!              
!*******************************************************************************

subroutine ti_change_weights(clambda, lambda_index)
  use parallel_dat_mod
  use prmtop_dat_mod
  use mdin_ctrl_dat_mod, only : gti_ele_sc, gti_vdw_sc, gti_lam_sch

  implicit none

! Formal arguments:
  double precision      :: clambda
  integer, optional :: lambda_index
! Local variables:
  double precision      :: w0, w1
  double precision      :: lam_sch_weights(2), lam_sch_dweights(2)
  integer               :: i, j, tt

#ifdef GTI
  !! Using c++ routine to ensure consistency

  if(ti_lam_sch.ne.0) call gti_reaf_lambda_schedule(lambda_scheduling_name)

  i=-1
  if (present(lambda_index))  i = lambda_index
  call gti_update_lambda(clambda, i-1, klambda_lcl, ti_weights, ti_dweights, &
     ti_item_weights, ti_item_dweights, ti_lam_sch.ne.0)

#else
  if(gti_lam_sch .eq. 1) then
    w0 = 1.d0 - clambda
    w1 = 1.d0 - w0
    ti_weights(1) = w0
    ti_weights(2) = w1
    ti_dweights(1) = -1
    ti_dweights(2) = 1
    lam_sch_weights(1) = 6.d0*(w0**5)-15.d0*(w0**4)+10.d0*(w0**3)
    lam_sch_weights(2) = 6.d0*(w1**5)-15.d0*(w1**4)+10.d0*(w1**3)
    lam_sch_dweights(1) = -1.d0*(30*(w0**4)-60.d0*(w0**3)+30.d0*(w0**2))
    lam_sch_dweights(2) = 30.d0*(w0**4)-60.d0*(w0**3)+30.d0*(w0**2)
  else
    w0 = (1.d0 - clambda) ** klambda_lcl
    w1 = 1.d0 - w0
    ti_weights(1) = w0
    ti_weights(2) = w1
    ti_dweights(1) = -1
    ti_dweights(2) = 1
  end if

  ! TypeGen = 1, TypeBAT = 2, TypeEleRec = 3, TypeEleCC = 4, TypeEleSC = 5, TypeVDW = 6, TypeRestBA = 7, TypeEleSSC = 8, TypeRMSD=9

  do i=1,ItemTypeTotal
    ti_item_weights(i,:) = ti_weights(:)
    ti_item_dweights(i,:) = ti_dweights(:)
  end do
  
  if(gti_lam_sch .eq. 1) then
    do i=4,6
      ti_item_weights(i,:) = lam_sch_weights(:)
      ti_item_dweights(i,:) = lam_sch_dweights(:)
    end do
  end if
#endif

  ti_signs(1) = -1.d0
  ti_signs(2) = 1.d0

  ! klambda scaling
  if (ti_mode .eq. 1) then
     if (klambda_lcl .ne. 1) then
        ti_klambda_factor = klambda_lcl * &
           ((1.d0 - clambda) ** (klambda_lcl - 1))
     else
        ti_klambda_factor = 1.d0
     end if
  else
     ti_klambda_factor = 1.d0
  end if

#if defined(MPI)
  call mpi_bcast(ti_weights, 2, mpi_double_precision, 0, &
                 pmemd_comm, err_code_mpi)
  call mpi_bcast(ti_dweights, 2, mpi_double_precision, 0, &
                 pmemd_comm, err_code_mpi)
  call mpi_bcast(ti_item_weights, SIZE(ti_item_weights), mpi_double_precision, 0, &
                 pmemd_comm, err_code_mpi)
  call mpi_bcast(ti_item_dweights, SIZE(ti_item_dweights), mpi_double_precision, 0, &
                 pmemd_comm, err_code_mpi)
  call mpi_bcast(ti_klambda_factor, 1, mpi_double_precision, 0, & 
                 pmemd_comm, err_code_mpi)                         
  call mpi_bcast(ti_signs, 2, mpi_double_precision, 0, &
                 pmemd_comm, err_code_mpi)
#endif

  return

end subroutine ti_change_weights

!*******************************************************************************
!
! Subroutine:  ti_degcnt
!
! Description: Calculates degrees of freedom for TI regions. Also calculates 
! temperature scaling factors.
!
!*******************************************************************************

subroutine ti_degcnt(atm_cnt, num_noshake, boltz2)
  use file_io_dat_mod
  use parallel_dat_mod
  use mdin_ctrl_dat_mod
  implicit none

! Formal arguments:
  integer, intent(in)           :: atm_cnt
  integer, intent(in)           :: num_noshake
  double precision, intent(in)  :: boltz2 

! Local variables:
  integer                       :: i
  integer                       :: sc_deg_lcl
  integer                       :: nadj
  integer                       :: sc_dof_cor
  integer                       :: sc_atm_cnt

  do i = 1, 2
     if (nscm .gt. 0 .and. ntb .eq. 0) then
        ti_ndfmin(i) = 6   ! both translation and rotation com motion removed
        nadj = atm_cnt - ti_ti_atm_cnt(2-i+1)
        if (nadj == 1) ti_ndfmin(i) = 3
        if (nadj == 2) ti_ndfmin(i) = 5
     else if (nscm .gt. 0) then
        ti_ndfmin(i) = 3    ! just translation com will be removed
     else
        ti_ndfmin(i) = 0
     end if  

     ! No COM motion removal for LD simulation
     if (gamma_ln .gt. 0.0d0) ti_ndfmin(i) = 0 
           
     ! Total dof
     ti_rndfp(i) = ti_rndfp(i) - dble(ti_ndfmin(i)) + dble(num_noshake) + &
                   dble(sc_num_noshake(i))
     ti_rndf(i) = ti_rndfp(i) + ti_rndfs(i) - 3.d0 * dble(ti_numextra_pts(i))

     ! Temperature scaling factors
     ti_fac(i,1) = boltz2 * ti_rndf(i)
     ti_fac(i,2) = boltz2 * ti_rndfp(i)
     if (ti_rndfp(i) .lt. 0.1d0) ti_fac(i,2) = 1.d-6
     ti_fac(i,3) = boltz2 * ti_rndfs(i)
     if (ti_rndfs(i) .lt. 0.1d0) ti_fac(i,3) = 1.d-6

     ti_factt(i) = ti_rndf(i) / (ti_rndf(i) + ti_ndfmin(i))

     ti_ekin0(i) = ti_fac(i,1) * temp0
     sc_atm_cnt = ti_ti_atm_cnt(i) - ti_latm_cnt(i)
     if (sc_atm_cnt .gt. 0 .and. ti_mode .ne. 1) then
        ! For now just assume that each SC atom has its full 3 dof
        ! w/ no additional internal constraints
        sc_deg_lcl = 3.0d0 * sc_atm_cnt

        ! Adjust for SHAKE constraints
        sc_dof_cor = ti_dof_shaked(i) - sc_num_noshake(i) 
        sc_deg_lcl = sc_deg_lcl - sc_dof_cor

        if (sc_dof_cor .eq. 0) then
           if (master) write(mdout,'(a,i2,a,i4)')  '   DOF for the SC part ',i,&
                             ' of the system: ', sc_deg_lcl
        else
           if (master) write(mdout,'(a,i2,a,i4,/,a,i4)')  &
              '   DOF for the SC part ',i,' of the system: ', &
              sc_deg_lcl,'   SHAKE constraints in the SC region: ', &
              sc_dof_cor
        end if
        ti_sc_fac(i) = 2 * 4.184 / ( 8.31441d-3 * sc_deg_lcl )
     else
        ti_sc_fac(i) = 0.0d0
     end if

  end do

  return

end subroutine ti_degcnt

!*******************************************************************************
!
! Subroutine:  ti_zero_arrays
!
! Description: Zero energy/force/virial arrays. Called by pme_force every frame.
!
!*******************************************************************************

subroutine ti_zero_arrays
  
  implicit none

  ! Energies
  ti_ene(:,:) = 0.0d0  
  ti_ene_aug(:,:) = 0.d0
  ti_ene_delta(:) = 0.d0
  if (ifmbar_lcl .ne. 0) then
    bar_cont(:) = 0.d0 
  end if
  ti_pot_gb(:,:) = 0.d0

  ! Forces
  ti_net_frcs(:,:) = 0.d0

  ! Virials
  ti_vir(:,:,:) = 0.d0 
  ti_vve(:) = 0.d0
  ti_molvir_netfrc_corr(:,:,:) = 0.d0

  return

end subroutine ti_zero_arrays

!*******************************************************************************
!
! Subroutine:  ti_update_ene
!
! Description: Update dvdl, delta V and mbar for a single potential (V0 or V1).
!              
!*******************************************************************************

subroutine ti_update_ene(pot_ene, ene_type, region, lambdaType)

  implicit none

! Formal arguments:
  double precision, intent(in) :: pot_ene
  integer, intent(in)          :: ene_type
  integer, intent(in)          :: region
  integer, intent(in)          :: lambdaType

! Local variables:
  integer                      :: i
  double precision             :: dvdl_temp
  
  dvdl_temp = pot_ene * ti_item_dweights(lambdaType, region)
  ti_ene(1,si_dvdl) = ti_ene(1,si_dvdl) + dvdl_temp
  ti_ene_delta(ene_type) = ti_ene_delta(ene_type) + dvdl_temp
  
  if (ifmbar_lcl .ne. 0 .and. do_mbar) then
    do i = 1, bar_states
      bar_cont(i) = bar_cont(i) + &
         pot_ene * (bar_lambda(region,i, lambdaType) - ti_item_weights(lambdaType,region))
    end do
  end if

  return

end subroutine ti_update_ene

!*******************************************************************************
!
! Subroutine:  ti_update_ene_all
!
! Description:  Update dvdl, delta V and mbar for both potentials (V0 and V1).
!              
!*******************************************************************************

subroutine ti_update_ene_all(pot_ene, ene_type, ene_combined, lambdaType)

  implicit none

! Formal arguments:
  double precision, intent(in) :: pot_ene(2)
  integer, intent(in)          :: ene_type
  double precision, intent(out):: ene_combined
  integer, optional :: lambdaType

! Local variables:
  integer                      :: i, l_type
  double precision             :: dvdl_temp
  
  if(present(lambdaType))then
    l_type = lambdaType
  else
    l_type = 1
  endif

  dvdl_temp = sum(pot_ene * ti_item_dweights(l_type,:))

  ti_ene(1,si_dvdl) = ti_ene(1,si_dvdl) + dvdl_temp
  ti_ene_delta(ene_type) = ti_ene_delta(ene_type) + dvdl_temp

  ene_combined = sum(pot_ene * ti_item_weights(l_type,:))

  if (ifmbar_lcl .ne. 0 .and. do_mbar) then
    do i = 1, bar_states
      bar_cont(i) = bar_cont(i) + &
       sum(pot_ene(:)*( bar_lambda(:,i,lambdaType)-ti_item_weights(l_type, :)) )
    end do
  end if

  if (ene_type .eq. si_elect_ene) then
    !only used for pme_err_est
    ti_vve(ti_ene0) = ti_vve(ti_ene0) + ti_pot(1)
    ti_vve(ti_ene1) = ti_vve(ti_ene1) + ti_pot(2)
    !store the self terms
    !ti_ene(:,ene_type)=pot_ene(:)
  end if

  if (ene_type .eq. si_vdw_ene) then
    !store the vdw-corr terms
    !ti_ene(:,ene_type)=pot_ene(:)
  end if

  return

end subroutine ti_update_ene_all

!*******************************************************************************
!
! Subroutine:  ti_update_ene_all_gamd
!
! Description:  Update dvdl, delta V and mbar for both potentials (V0 and V1).
!              
!*******************************************************************************

subroutine ti_update_ene_all_gamd(pot_ene, ene_type, ene_combined)

  implicit none

! Formal arguments:
  double precision, intent(in) :: pot_ene(2)
  integer, intent(in)          :: ene_type
  double precision, intent(out):: ene_combined
! Local variables:
  integer                      :: i
  double precision             :: dvdl_temp
  
  dvdl_temp = pot_ene(2) - pot_ene(1)

  ti_ene(1,si_dvdl) = ti_ene(1,si_dvdl) + dvdl_temp
  ti_ene_delta(ene_type) = ti_ene_delta(ene_type) + dvdl_temp

  ene_combined = pot_ene(1) + pot_ene(2)

  if (ifmbar_lcl .ne. 0 .and. do_mbar) then
    do i = 1, bar_states
      bar_cont(i) = bar_cont(i) + &
         pot_ene(1) * (bar_lambda(1,i,1) - ti_weights(1)) + &
         pot_ene(2) * (bar_lambda(2,i,1) - ti_weights(2))
    end do
  end if

  if (ene_type .eq. si_elect_ene) then
    !only used for pme_err_est
    ti_vve(ti_ene0) = ti_vve(ti_ene0) + ti_pot(1)
    ti_vve(ti_ene1) = ti_vve(ti_ene1) + ti_pot(2)
  end if

  return

end subroutine ti_update_ene_all_gamd

!*******************************************************************************
!
! Subroutine:  ti_update_vve
!
! Description: Update data for virial vs ene calculation. 
!              Only used in nb_adjust in nb_exclusions.
!              
!*******************************************************************************

subroutine ti_update_vve(pot_ene, vir, region)

  implicit none

! Formal arguments:
  double precision, intent(in) :: pot_ene
  double precision, intent(in) :: vir
  integer, intent(in)          :: region ! 0 means common, add to both

! Local variables:
  integer                      :: i
  double precision             :: vir_temp

  if (region .ne. 0) then
     ti_vve(region + ti_vir0 - 1) = ti_vve(region + ti_vir0 - 1) + vir
     ti_vve(region + ti_ene0 - 1) = ti_vve(region + ti_ene0 - 1) + pot_ene
  else
     ti_vve(ti_vir0) = ti_vve(ti_vir0) + vir * ti_weights(1)
     ti_vve(ti_vir1) = ti_vve(ti_vir1) + vir * ti_weights(2)

     ti_vve(ti_ene0) = ti_vve(ti_ene0) + pot_ene
     ti_vve(ti_ene1) = ti_vve(ti_ene1) + pot_ene
  end if
  
  return

end subroutine ti_update_vve

!*******************************************************************************
!
! Subroutine:  ti_update_nb_frc
!
! Description: Update ti_nb_frc. Used by get_nb14_energy and nb_adjust.
!              
!*******************************************************************************

subroutine ti_update_nb_frc(dfx, dfy, dfz, atm_i, region, itemType)

  implicit none

! Formal arguments:
 
  double precision, intent(in) :: dfx, dfy, dfz
  integer, intent(in)          :: atm_i
  integer, intent(in)          :: region
  integer, intent(in)          :: itemType

! Local variables:

  if (region .eq. 0) then
     ti_nb_frc(1, 1, atm_i) = ti_nb_frc(1, 1, atm_i) + dfx * ti_item_weights(itemType,1)
     ti_nb_frc(1, 2, atm_i) = ti_nb_frc(1, 2, atm_i) + dfy * ti_item_weights(itemType,1)
     ti_nb_frc(1, 3, atm_i) = ti_nb_frc(1, 3, atm_i) + dfz * ti_item_weights(itemType,1)
     ti_nb_frc(2, 1, atm_i) = ti_nb_frc(2, 1, atm_i) + dfx * ti_item_weights(itemType,2)
     ti_nb_frc(2, 2, atm_i) = ti_nb_frc(2, 2, atm_i) + dfy * ti_item_weights(itemType,2)
     ti_nb_frc(2, 3, atm_i) = ti_nb_frc(2, 3, atm_i) + dfz * ti_item_weights(itemType,2)
  else
     ti_nb_frc(region, 1, atm_i) = ti_nb_frc(region, 1, atm_i) + dfx 
     ti_nb_frc(region, 2, atm_i) = ti_nb_frc(region, 2, atm_i) + dfy 
     ti_nb_frc(region, 3, atm_i) = ti_nb_frc(region, 3, atm_i) + dfz
  end if

  return

end subroutine ti_update_nb_frc

!*******************************************************************************
!
! Subroutine:  ti_check_res_int
!
! Description: Parse the nmrat/nmrcom arrays to make sure no restraints have
!              been defined between the non-interacting TI regions.
!              Only used during initialization, see ti_check_res below for
!              the code that runs during a simulation.
!
!*******************************************************************************

subroutine ti_check_res_int(nmrnum, nmrat, nmrcom, iout)

  use file_io_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:
  ! The nmr code will read past what is defined here as the bounds for the array
  integer, intent(in) :: nmrnum
  integer, intent(in) :: nmrat(16, nmrnum)
  integer, intent(in) :: nmrcom(2, nmrnum)
  integer, intent(in) :: iout
! Local variables:
  integer             :: iat
  integer             :: jat
  integer             :: ip
  integer             :: i
  integer             :: j
  integer             :: max_iat
  integer             :: iat_num
  integer             :: res_region(2)

  do i = 1, nmrnum
     res_region(:) = 0
     if (nmrat(9, i) .lt. 0) then
        max_iat = 2
     else if(nmrat(10, i) .lt. 0) then
        max_iat = 3
     else
        max_iat = 4
     end if

     do iat = 1, max_iat
        jat = nmrat(iat, i)
        if (jat .lt. 0) then
           do ip = -nmrat(iat, i), -nmrat(8 + iat, i)
              do j = nmrcom(1, ip), nmrcom(2, ip)
                 iat_num = j
                 if(ti_lst(1,iat_num) .ne. 0) then
                    res_region(1) = 1
                 else if(ti_lst(2,iat_num) .ne. 0) then
                    res_region(2) = 1
                 end if
              end do
           end do
        else
           iat_num = (jat/3) + 1
           if(ti_lst(1,iat_num) .ne. 0) then
              res_region(1) = 1
           else if(ti_lst(2,iat_num) .ne. 0) then
              res_region(2) = 1
           end if      
        end if
     end do

     if (res_region(1) .ne. 0 .and. res_region(2) .ne. 0) then
     !   write (iout,'(a)') '    ERROR: Restraints can not be defined between &
     !      &non-interacting TI regions.'
     !   call mexit(iout, 1)
     end if
  end do

  return

end subroutine ti_check_res_int

!*******************************************************************************
!
! Subroutine:  ti_check_res
!
! Description: Parse the nmrat/nmrcom arrays to determinine if the restraint ene
!              should be counted as "common" or put into the sc ene arrays.
!              res_ti_region is 0 for common, or 1,2 for ti_region 1,2.
!              See nmrcms in nmr_lib.fpp. NOTE: Restraints can not be defined
!              between ti_regions as they are non interacting by definition.
!
!*******************************************************************************

subroutine ti_check_res(nmrat, nmrcom, i)

  implicit none

! Formal arguments:
  ! The nmr code will read past what is defined here as the bounds for the array
  integer, intent(in) :: nmrat(16, *)
  integer, intent(in) :: nmrcom(2, *)
  integer, intent(in) :: i  !index into nmrat array

! Local variables:
  integer             :: iat
  integer             :: jat
  integer             :: ip
  integer             :: j
  integer             :: max_iat
  integer             :: iat_num

  res_ti_region = 0

  if (nmrat(9, i) .lt. 0) then
     max_iat = 2
  else if(nmrat(10, i) .lt. 0) then
     max_iat = 3
  else
     max_iat = 4
  end if

  do iat = 1, max_iat

     jat = nmrat(iat, i)
     if (jat .lt. 0) then
        do ip = -nmrat(iat, i), -nmrat(8 + iat, i)
           do j = nmrcom(1, ip), nmrcom(2, ip)
              iat_num = j
              if(ti_lst(1,iat_num) .ne. 0) then
                 res_ti_region = 1
                 return
              else if(ti_lst(2,iat_num) .ne. 0) then
                 res_ti_region = 2
                 return
              end if
           end do
        end do
     else
        iat_num = (jat/3) + 1
        if(ti_lst(1,iat_num) .ne. 0) then
           res_ti_region = 1
           return
        else if(ti_lst(2,iat_num) .ne. 0) then
           res_ti_region = 2
           return
        end if      
     end if
  end do

  return

end subroutine ti_check_res

!*******************************************************************************
!
! Subroutine:  ti_check_res_lscale
!
! Description: Parse the nmrat/nmrcom arrays to determinine if the restraint ene
!              should be considered linear or not for ifsc > 0.
!
!*******************************************************************************

subroutine ti_check_res_lscale(nmrat, nmrcom, i)

  implicit none

! Formal arguments:
  ! The nmr code will read past what is defined here as the bounds for the array
  integer, intent(in) :: nmrat(16, *)
  integer, intent(in) :: nmrcom(2, *)
  integer, intent(in) :: i  !index into nmrat array

! Local variables:
  integer             :: iat
  integer             :: jat
  integer             :: ip
  integer             :: j
  integer             :: max_iat
  integer             :: iat_num

  ti_lscale = .false.

  if (nmrat(9, i) .lt. 0) then
     max_iat = 2
  else if(nmrat(10, i) .lt. 0) then
     max_iat = 3
  else
     max_iat = 4
  end if

  do iat = 1, max_iat
     jat = nmrat(iat, i)
     if (jat .lt. 0) then
        do ip = -nmrat(iat, i), -nmrat(8 + iat, i)
           do j = nmrcom(1, ip), nmrcom(2, ip)
              iat_num = j
              if(ti_sc_lst(iat_num) .ne. 0) then
                 return              
              end if
           end do
        end do
     else
        iat_num = (jat/3) + 1
        if(ti_sc_lst(iat_num) .ne. 0) then
           return
        end if      
     end if
  end do

  ti_lscale = .true.
  return

end subroutine ti_check_res_lscale

!*******************************************************************************
!
! Subroutine:  ti_dist_enes_virs_netfrcs
!
! Description: Combine the ene in the TI arrays across nodes. 
!              
!*******************************************************************************

subroutine ti_dist_enes_virs_netfrcs(need_pot_enes, need_virials)

  use parallel_dat_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod

  implicit none

! Formal Arguments:
  logical, intent(in)           :: need_pot_enes
  logical, intent(in)           :: need_virials

! Local variables:

  type ti_pme_dat
    sequence
    double precision              :: ti_ene(2, si_cnt)
    double precision              :: ti_ene_aug(2, ti_ene_aug_cnt)
    double precision              :: ti_ene_delta(si_cnt)
    double precision              :: ti_molvir_netfrc_corr(2, 3, 3)
    double precision              :: ti_vve(ti_vir_ene_cnt)
    double precision              :: ti_net_frcs(2, 3)
  end type ti_pme_dat
 
  integer, parameter              :: ti_pme_dat_size = &
                                     2 * si_cnt + 2 * ti_ene_aug_cnt + &
                                     si_cnt + 2 * 3 * 3 + ti_vir_ene_cnt + &
                                     2 * 3

  integer, parameter              :: ti_pme_dat_no_pot_ene_size = &
                                     2 * 3 * 3 + ti_vir_ene_cnt + 2 * 3

  type(ti_pme_dat), save          :: dat_in
  type(ti_pme_dat), save          :: dat_out

 ! Sent seperately as bar_states is not a constant
  double precision                :: bar_cont_tmp(bar_states)

  integer                         :: i
  integer                         :: buf_size

  if (need_pot_enes .or. nmropt .ne. 0 .or. verbose .ne. 0 .or. infe .ne. 0) then

#if defined(MPI)

     dat_in%ti_ene(:,:)                  = ti_ene(:,:)
     dat_in%ti_ene_aug(:,:)              = ti_ene_aug(:,:)
     dat_in%ti_ene_delta(:)              = ti_ene_delta(:)
     dat_in%ti_molvir_netfrc_corr(:,:,:) = ti_molvir_netfrc_corr(:,:,:)
     dat_in%ti_vve(:)                    = ti_vve(:)
     dat_in%ti_net_frcs(:,:)             = ti_net_frcs(:,:)

     buf_size = ti_pme_dat_size

     call mpi_allreduce(dat_in%ti_ene(1,1), dat_out%ti_ene(1,1), buf_size, &
                        mpi_double_precision, mpi_sum, pmemd_comm, err_code_mpi)

     ti_ene(:,:)                  = dat_out%ti_ene(:,:)
     ti_ene_aug(:,:)              = dat_out%ti_ene_aug(:,:)
     ti_ene_delta(:)              = dat_out%ti_ene_delta(:)
     ti_molvir_netfrc_corr(:,:,:) = dat_out%ti_molvir_netfrc_corr(:,:,:)
     ti_vve(:)                    = dat_out%ti_vve(:)
     ti_net_frcs(:,:)             = dat_out%ti_net_frcs(:,:)

     if (ifmbar .ne. 0 .and. do_mbar) then
        call mpi_allreduce(bar_cont, bar_cont_tmp, bar_states, &
           mpi_double_precision, mpi_sum, pmemd_comm, err_code_mpi)
        bar_cont(:) = bar_cont_tmp(:)
     end if

  else if(need_virials) then

     dat_in%ti_molvir_netfrc_corr(:,:,:) = ti_molvir_netfrc_corr(:,:,:)
     dat_in%ti_vve(:)                    = ti_vve(:)
     dat_in%ti_net_frcs(:,:)             = ti_net_frcs(:,:)

     buf_size = ti_pme_dat_no_pot_ene_size

     call mpi_allreduce(dat_in%ti_molvir_netfrc_corr(1,1,1), &
                        dat_out%ti_molvir_netfrc_corr(1,1,1), buf_size, &
                        mpi_double_precision, mpi_sum, pmemd_comm, err_code_mpi)

     ti_molvir_netfrc_corr(:,:,:) = dat_out%ti_molvir_netfrc_corr(:,:,:)
     ti_vve(:)                    = dat_out%ti_vve(:)
     ti_net_frcs(:,:)             = dat_out%ti_net_frcs(:,:)

  else if(netfrc .ne. 0) then

     dat_in%ti_net_frcs(:,:)             = ti_net_frcs(:,:)

     buf_size = 2 * 3

     call mpi_allreduce(dat_in%ti_net_frcs(1,1), dat_out%ti_net_frcs(1,1), &
                        buf_size, mpi_double_precision, mpi_sum, pmemd_comm, &
                        err_code_mpi)

     ti_net_frcs(:,:)             = dat_out%ti_net_frcs(:,:)

#endif
  end if

  return

end subroutine ti_dist_enes_virs_netfrcs

!*******************************************************************************
!
! Subroutine:  ti_calc_dvdl
!
! Description: Calculate total energy and scale dvdl by klambda.
!              
!*******************************************************************************

subroutine ti_calc_dvdl

  use parallel_dat_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod

  implicit none

! Local variables:

  integer                         :: i

   ! Apply klambda scaling
   if (klambda_lcl .ne. 1) then
      ti_ene(1, si_dvdl) = ti_ene(1, si_dvdl) * ti_klambda_factor
      ti_ene_delta(:) = ti_ene_delta(:) * ti_klambda_factor
   end if
        
   ! Region 1 is used to store dvdl    
   ti_ene(2, si_dvdl) = ti_ene(1, si_dvdl) 
  
   if (logdvdl .ne. 0 .and. master) then
      ti_dvdl_pnt = ti_dvdl_pnt + 1
      ti_dvdl_values(ti_dvdl_pnt) = ti_ene(1, si_dvdl)
   end if

   ! Sum up potential energy terms
   do i = 1, 2
      ti_ene(i,si_restraint_ene) = &
         ti_ene_aug(i,ti_rest_dist_ene) + &
         ti_ene_aug(i,ti_rest_ang_ene) + &
         ti_ene_aug(i,ti_rest_tor_ene)

      ! NOTE: This is the pot ene without restraint ene
      ti_ene(i,si_pot_ene) = &
         ti_ene(i,si_vdw_ene) + &
         ti_ene(i,si_elect_ene) + &
         ti_ene(i,si_hbond_ene) + &
         ti_ene(i,si_bond_ene) + &
         ti_ene(i,si_angle_ene) + &
         ti_ene(i,si_dihedral_ene) + &
         ti_ene(i,si_vdw_14_ene) + &
         ti_ene(i,si_elect_14_ene) + &
         ti_ene(i,si_angle_ub_ene) + &
         ti_ene(i,si_dihedral_imp_ene) + &
         ti_ene(i,si_cmap_ene) + &
         ti_ene(i,si_surf_ene)
   end do

   ! dvdl should be the same as the sum of the difference
   ! of potential energy terms
   ti_ene_delta(si_dvdl) = ti_ene(1, si_dvdl)      
   ti_ene_delta(si_pot_ene) = ti_ene(1, si_dvdl) 

  return

end subroutine ti_calc_dvdl

!*******************************************************************************
!
! Subroutine:  ti_combine_frcs
!
! Description: Combines the TI forces into the main frc array.
!              
!*******************************************************************************

#if defined(MPI) && !defined(CUDA)
subroutine ti_combine_frcs(atm_cnt, my_atm_cnt, my_atm_lst, frc, ti_frc)
#else
subroutine ti_combine_frcs(atm_cnt, frc, ti_frc)
#endif
  implicit none

! Formal arguments:
  integer, intent(in)              :: atm_cnt
#if defined(MPI) && !defined(CUDA)
  integer, intent(in)              :: my_atm_cnt
  integer, intent(in)              :: my_atm_lst(my_atm_cnt)
#endif
  double precision, intent(inout)  :: frc(3, atm_cnt)
  double precision, intent(in)     :: ti_frc(2, 3, atm_cnt)

! Local variables:
  integer                          :: j
  integer                          :: m
  integer                          :: iatom
  integer                          :: z

#if defined(MPI) && !defined(CUDA)
  do iatom = 1, my_atm_cnt
     j = my_atm_lst(iatom)
#else
  do iatom = 1, atm_cnt
     j = iatom
#endif
     frc(:,j) = 0.d0
     do z = 1, 2
        do m = 1, 3
           frc(m,j) = frc(m,j) + ti_frc(z,m,j)
        end do
     end do
  end do     

  return

end subroutine ti_combine_frcs

!*******************************************************************************
!
! Subroutine:  ti_calc_kin_ene
!
! Description: Calculates kinetic energy for TI atoms.
!              
!*******************************************************************************

#if defined(MPI) && !defined(CUDA)
subroutine ti_calc_kin_ene(atm_cnt, my_atm_cnt, my_atm_lst, amass, vel, vold,&
                            c_ave, mode)
#else
subroutine ti_calc_kin_ene(atm_cnt, amass, vel, vold, c_ave, mode)
#endif
  implicit none

! Formal arguments:
  integer, intent(in)           :: atm_cnt
#if defined(MPI) && !defined(CUDA)
  integer, intent(in)           :: my_atm_cnt
  integer, intent(in)           :: my_atm_lst(my_atm_cnt)
#endif
  double precision, intent(in)  :: amass(atm_cnt)
  double precision, intent(in)  :: vel(3, atm_cnt)
  double precision, intent(in)  :: vold(3, atm_cnt)
  double precision, intent(in)  :: c_ave
  integer, intent(in)           :: mode ! 1 for eke, 2 for ekpbs, 3 for ekph

! Local variables:
  double precision              :: aamass
  integer                       :: j
  integer                       :: m
  integer                       :: iatom
  integer                       :: z
  double precision              :: val
  logical                       :: calc_sc_kin

  calc_sc_kin = mode .eq. ti_eke .or. mode .eq. ti_ekmh !for GPU with ntt=3,
!always ti_eke
  ti_kin_ene(:,mode) = 0.d0 !ti_eke = 1
  if (calc_sc_kin) ti_kin_ene(:,ti_sc_eke) = 0.d0 !ti_sc_eke = 2
  do z = 1, 2
#if defined(MPI) && !defined(CUDA)
     do iatom = 1, my_atm_cnt
        j = my_atm_lst(iatom)
#else
     do iatom = 1, atm_cnt
        j = iatom
#endif
        if (ti_lst(z,j) .eq. 1) then
!for ti atoms
           aamass = amass(j) !get mass, store temporarily
           do m = 1, 3     !over x,y,z
              if(calc_sc_kin) then 
!c_ave = 1 + gammai * half_dtx
!gammai = gamma_ln/20.455
!gamma_ln = 1
                 val = aamass * 0.25d0 * c_ave * ( vel(m,j) + vold(m,j) ) ** 2
                 ti_kin_ene(z,mode) = ti_kin_ene(z,mode) + val                    
                 if (ti_mode .ne. 1 .and. ti_sc_lst(j) .ge. 1) then
!also increment sc kinetic energy
                   ti_kin_ene(z,ti_sc_eke) = ti_kin_ene(z,ti_sc_eke) + val
                 end if 
!ignore the else for ntt=3
              else !ekpbs/ekph
                 ti_kin_ene(z,mode) = ti_kin_ene(z,mode) + &
                    aamass * ( vel(m,j) * vold(m,j) )
              end if
           end do
        end if
     end do
  end do

  ti_kin_ene(:,mode) = ti_kin_ene(:,mode) * 0.5d0

  if (calc_sc_kin) ti_kin_ene(:,ti_sc_eke) = ti_kin_ene(:,ti_sc_eke)*0.5d0

  return

end subroutine ti_calc_kin_ene

!*******************************************************************************
!
! Subroutine:  ti_copy_common_core
!
! Description: Copy TI Mol 1 common core to TI Mol 2 for post shake
!              
!*******************************************************************************

subroutine ti_copy_common_core(atm_cnt, crd)

  implicit none

  integer, intent(in)       :: atm_cnt
  double precision          :: crd(3, atm_cnt)

  integer                   :: i, atm_i, atm_j

  do i=1, ti_latm_cnt(1)
      atm_i = ti_latm_lst(1,i)
      call ti_get_partner(ti_latm_lst(1,i), atm_j)
      crd(1,atm_j)=crd(1,atm_i)
      crd(2,atm_j)=crd(2,atm_i)
      crd(3,atm_j)=crd(3,atm_i)
  end do

end subroutine ti_copy_common_core

!*******************************************************************************
!
! Subroutine:  ti_exchange_vec
!
! Description: For linear scaling, combine/bcast the contributions for each TI
!              region. Combine determines whether we add contributions from both
!              regions (.true.) or just push from V0 to V1 (ie for vel).
!
!*******************************************************************************

subroutine ti_exchange_vec(atm_cnt, vec, combine)

  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  implicit none

! Formal arguments:
  integer, intent(in)             :: atm_cnt
  double precision, intent(inout) :: vec(3, atm_cnt) 
  logical, intent(in)             :: combine

#if defined(MPI) && !defined(CUDA)
! Local variables:
  double precision                :: vec_tmp(3, ti_latm_cnt(1))
  double precision                :: vec_tmp_red(3, ti_latm_cnt(1))
#endif
  integer                         :: i, j, m
  integer                         :: atm_i, atm_j

#if defined(MPI) && !defined(CUDA)
  if (imin .eq. 0) then

     !we need to all reduce the vector across nodes
     vec_tmp(:,:) = 0.d0
     !first put all atoms that we own into vec_tmp
     do i = 1, ti_latm_cnt(1)
        atm_i = ti_latm_lst(1,i)
        atm_j = ti_latm_lst(2,i)
        if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
           do m = 1, 3
              vec_tmp(m, i) = vec_tmp(m, i) + vec(m, atm_i)                 
           end do
        end if
        if ((gbl_atm_owner_map(atm_j) .eq. mytaskid) .and. combine) then
           do m = 1, 3
              vec_tmp(m, i) = vec_tmp(m, i) + vec(m, atm_j)
           end do
        end if
     end do

     ! Now all reduce the tmp array
     ! Ideally, we could use a more complex redistribution scheme
     ! but this is fine for now ...
     call mpi_allreduce(vec_tmp, vec_tmp_red, 3 * ti_latm_cnt(1), &
                        mpi_double_precision, mpi_sum, pmemd_comm, err_code_mpi)

     ! Finally restore the tmp array
     do i = 1, ti_latm_cnt(1)
        atm_i = ti_latm_lst(1,i)
        atm_j = ti_latm_lst(2,i)
        if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
           do m = 1, 3
              vec(m, atm_i) = vec_tmp_red(m, i)
           end do
        end if
        if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
           do m = 1, 3
              vec(m, atm_j) = vec_tmp_red(m, i)
           end do
        end if
     end do

  else
#endif
     !only one node or minimization
     do i = 1, ti_latm_cnt(1)
        atm_i = ti_latm_lst(1,i)
        atm_j = ti_latm_lst(2,i)
        do m = 1, 3
           if (combine) then
              vec(m, atm_i) = vec(m, atm_i) + vec(m, atm_j)
           end if
           vec(m, atm_j) = vec(m, atm_i)
        end do      
     end do

#if defined(MPI) && !defined(CUDA)
  end if
#endif
  return

end subroutine ti_exchange_vec

!*******************************************************************************
!
! Subroutine:   ti_vrand_set_velocities
!
! Description:  Assign velocities from a Maxwellian distribution. From 
!               dynamics.fpp, with modifications for TI. In sander, the 
!               velocities are pushed from V0 to V1 for the common atoms, 
!               do the same here.
!              
!*******************************************************************************

subroutine ti_vrand_set_velocities(atm_cnt, vel, mass_inv, temp0)
   
  use parallel_dat_mod
  use random_mod

  implicit none

! Formal arguments:
  integer, intent(in)           :: atm_cnt
  double precision, intent(out) :: vel(3, atm_cnt)
  double precision, intent(in)  :: mass_inv(atm_cnt)
  double precision, intent(in)  :: temp0

! Local variables:
  double precision      :: boltz
  double precision      :: sd
  double precision      :: temp(3)
  integer               :: j
  integer               :: i 
  integer               :: atm_id
  integer               :: nreg(3)

  temp(1) = temp0 * ti_factt(1)    
  temp(2) = temp0 * ti_factt(2)
  temp(3) = temp0 * ti_factt(1) !common atoms use V0

  nreg(1) = ti_ti_atm_cnt(1)
  nreg(2) = ti_ti_atm_cnt(2)
  nreg(3) = ti_atm_cnt(1) - ti_ti_atm_cnt(1)

  do i = 1, 3 !V0/V1/common
     if (temp(i) .lt. 1.d-6) then
        do atm_id = 1, nreg(i)
           j = ti_atm_lst(i,atm_id)
#if defined(MPI) && !defined(CUDA)
           if (gbl_atm_owner_map(j) .eq. mytaskid) then
#endif         
              vel(:,j) = 0.d0
#if defined(MPI) && !defined(CUDA)
           end if
#endif
        end do      
     else

        boltz = 8.31441d-3 * temp(i) / 4.184d0

        do atm_id = 1, nreg(i)
           j = ti_atm_lst(i,atm_id)

#if defined(MPI) && !defined(CUDA)
  ! In order to generate the same sequence of pseudorandom numbers that you
  ! would using a single processor or any other combo of multiple processors,
  ! you have to go through the atoms in order.  The unused results are not
  ! returned.

           if (gbl_atm_owner_map(j) .eq. mytaskid) then
              sd =  sqrt(boltz * mass_inv(j))
              call gauss(0.d0, sd, vel(1, j))
              call gauss(0.d0, sd, vel(2, j))
              call gauss(0.d0, sd, vel(3, j))
           else
              call gauss(0.d0, 1.d0)
              call gauss(0.d0, 1.d0)
              call gauss(0.d0, 1.d0)
           end if
#else
           sd =  sqrt(boltz * mass_inv(j))
           call gauss(0.d0, sd, vel(1, j))
           call gauss(0.d0, sd, vel(2, j))
           call gauss(0.d0, sd, vel(3, j))
#endif

        end do
     end if
  end do

  return

end subroutine ti_vrand_set_velocities

!*******************************************************************************
!
! Subroutine:  ti_update_avg_ene
!
! Description: Update averages for the softcore arrays.
!              
!*******************************************************************************

subroutine ti_update_avg_ene

  implicit none

! Local variables:
  integer      i

  do i = 1, 2
    ti_ene_ave(i,:) = ti_ene_ave(i,:) + ti_ene(i,:)
    ti_ene_rms(i,:) = ti_ene_rms(i,:) + ti_ene(i,:) ** 2

    ti_others_ene_ave(i,:) = ti_others_ene_ave(i,:) + ti_others_ene(i,:)
    ti_others_ene_rms(i,:) = ti_others_ene_rms(i,:) + ti_others_ene(i,:) ** 2

    ti_ene_aug_ave(i,:) = ti_ene_aug_ave(i,:) + &
                               ti_ene_aug(i,:)
    ti_ene_aug_rms(i,:) = ti_ene_aug_rms(i,:) + &
                               ti_ene_aug(i,:) ** 2
  end do

  ti_ene_delta_ave(:) = ti_ene_delta_ave(:) + ti_ene_delta(:)
  ti_ene_delta_rms(:) = ti_ene_delta_rms(:) + &
     ti_ene_delta(:) ** 2

  return

end subroutine ti_update_avg_ene

!*******************************************************************************
!
! Subroutine:  ti_calc_avg_ene
!
! Description: Determine the avg/rms energies and store the results in
!              the tmp/tmp2 arrays.
!              
!*******************************************************************************

subroutine ti_calc_avg_ene(tspan, final)
  use file_io_dat_mod
  implicit none

! Formal arguments:
  double precision, intent(in) :: tspan
  logical, intent(in)          :: final !true for final output

! Local variables:

! Work around what seems to be a recently-introduced GNU compiler optimization
! bug. We want to prevent m and i from being optimized away, due to potential
! branching in the loops where they're used. But GCC 4.1.2 (actually <= 4.2)
! doesn't support the volatile keyword, (and the bug doesn't seem to exist
! there), so use this ugly workaround to continue GCC 4.1.2 support
#if !defined(__GNUC__) || (__GNUC__ >= 4 && __GNUC_MINOR__ > 2)
  integer, volatile            :: m
  integer, volatile            :: i
#else
  integer                      :: m
  integer                      :: i
#endif

  if (final) then
     ti_ene_old(:,:) = 0.0d0
     ti_ene_old2(:,:) = 0.0d0

     ti_others_ene_old(:,:) = 0.0d0
     ti_others_ene_old2(:,:) = 0.0d0

     ti_ene_aug_old(:,:) = 0.0d0
     ti_ene_aug_old2(:,:) = 0.0d0

     ti_ene_delta_old(:) = 0.0d0
     ti_ene_delta_old2(:) = 0.0d0
  end if

  do m = 1, 2
     do i = 1, si_cnt
        ti_ene_tmp(m,i) = ti_ene_ave(m,i)-ti_ene_old(m,i)
        ti_ene_tmp2(m,i) = ti_ene_rms(m,i)-ti_ene_old2(m,i)
        ti_ene_old(m,i) = ti_ene_ave(m,i)
        ti_ene_old2(m,i) = ti_ene_rms(m,i)
        ti_ene_tmp(m,i) = ti_ene_tmp(m,i)/tspan
        ti_ene_tmp2(m,i) = &
           ti_ene_tmp2(m,i)/tspan - ti_ene_tmp(m,i)**2
        if (ti_ene_tmp2(m,i) < 0.0d0) ti_ene_tmp2(m,i) = 0.0d0
        ti_ene_tmp2(m,i) = sqrt(ti_ene_tmp2(m,i))
     end do
  end do

  do m = 1, 2
     do i = 1, si_cnt
        ti_others_ene_tmp(m,i) = ti_others_ene_ave(m,i)-ti_others_ene_old(m,i)
        ti_others_ene_tmp2(m,i) = ti_others_ene_rms(m,i)-ti_others_ene_old2(m,i)
        ti_others_ene_old(m,i) = ti_others_ene_ave(m,i)
        ti_others_ene_old2(m,i) = ti_others_ene_rms(m,i)
        ti_others_ene_tmp(m,i) = ti_others_ene_tmp(m,i)/tspan
        ti_others_ene_tmp2(m,i) = &
           ti_others_ene_tmp2(m,i)/tspan - ti_others_ene_tmp(m,i)**2
        if (ti_others_ene_tmp2(m,i) < 0.0d0) ti_others_ene_tmp2(m,i) = 0.0d0
        ti_others_ene_tmp2(m,i) = sqrt(ti_others_ene_tmp2(m,i))
     end do
  end do

  do m = 1, 2
     do i = 1, ti_ene_aug_cnt
        ti_ene_aug_tmp(m,i) = &
           ti_ene_aug_ave(m,i)-ti_ene_aug_old(m,i)
        ti_ene_aug_tmp2(m,i) = &
           ti_ene_aug_rms(m,i)-ti_ene_aug_old2(m,i)
        ti_ene_aug_old(m,i) = ti_ene_aug_ave(m,i)
        ti_ene_aug_old2(m,i) = ti_ene_aug_rms(m,i)
        ti_ene_aug_tmp(m,i) = ti_ene_aug_tmp(m,i)/tspan
        ti_ene_aug_tmp2(m,i) = &
           ti_ene_aug_tmp2(m,i)/tspan - ti_ene_aug_tmp(m,i)**2
        if (ti_ene_aug_tmp2(m,i) < 0.0d0) ti_ene_aug_tmp2(m,i) = 0.0d0
        ti_ene_aug_tmp2(m,i) = sqrt(ti_ene_aug_tmp2(m,i))
     end do
  end do

  do i = 1, si_cnt
     ti_ene_delta_tmp(i) = &
        ti_ene_delta_ave(i)-ti_ene_delta_old(i)
     ti_ene_delta_tmp2(i) = &
        ti_ene_delta_rms(i)-ti_ene_delta_old2(i)
     ti_ene_delta_old(i) = ti_ene_delta_ave(i)
     ti_ene_delta_old2(i) = ti_ene_delta_rms(i)
     ti_ene_delta_tmp(i) = ti_ene_delta_tmp(i)/tspan
     ti_ene_delta_tmp2(i) = &
        ti_ene_delta_tmp2(i)/tspan - ti_ene_delta_tmp(i)**2
     if (ti_ene_delta_tmp2(i) < 0.0d0) ti_ene_delta_tmp2(i) = 0.0d0
     ti_ene_delta_tmp2(i) = sqrt(ti_ene_delta_tmp2(i))
  end do
      
  return

end subroutine ti_calc_avg_ene

!*******************************************************************************
!
! Subroutine:  ti_dynamic_lambda
!
! Description: Update lambda for dynlmb run.
!              
!*******************************************************************************

subroutine ti_dynamic_lambda

  use mdin_ctrl_dat_mod
  use parallel_dat_mod
    
  call ti_update_lambda(dynlmb+clambda)

  !clear all average arrays
  if (master) then
    write (mdout,'(/,a,f12.4,a,f12.4,/)') &
       'Dynamically changing lambda: Increased clambda by ', &
       dynlmb, ' to ', clambda
  end if

  return
end subroutine ti_dynamic_lambda


subroutine ti_update_lambda(lambda, lambdaIndex)

  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  
  double precision, intent(in):: lambda 
  integer,optional::lambdaIndex
  
  clambda=max(min(lambda,1.0),0.0)

  call ti_change_weights(clambda, lambdaIndex)

#ifdef CUDA
#ifndef GTI
  call gpu_ti_dynamic_lambda(clambda, lambdaIndex)
#endif  
#endif

  !clear all average arrays
  if (master) then
    ti_ene(:,:) = 0.d0   
    ti_ene_ave(:,:) = 0.d0
    ti_ene_rms(:,:) = 0.d0
    ti_ene_tmp(:,:) = 0.d0
    ti_ene_tmp2(:,:) = 0.d0
    ti_ene_old(:,:) = 0.d0
    ti_ene_old2(:,:) = 0.d0

    ti_others_ene(:,:) = 0.d0   
    ti_others_ene_ave(:,:) = 0.d0
    ti_others_ene_rms(:,:) = 0.d0
    ti_others_ene_tmp(:,:) = 0.d0
    ti_others_ene_tmp2(:,:) = 0.d0
    ti_others_ene_old(:,:) = 0.d0
    ti_others_ene_old2(:,:) = 0.d0

    ti_ene_aug(:,:) = 0.d0   
    ti_ene_aug_ave(:,:) = 0.d0
    ti_ene_aug_rms(:,:) = 0.d0
    ti_ene_aug_tmp(:,:) = 0.d0
    ti_ene_aug_tmp2(:,:) = 0.d0
    ti_ene_aug_old(:,:) = 0.d0
    ti_ene_aug_old2(:,:) = 0.d0

    ti_ene_delta(:) = 0.d0  
    ti_ene_delta_ave(:) = 0.d0
    ti_ene_delta_rms(:) = 0.d0
    ti_ene_delta_tmp(:) = 0.d0
    ti_ene_delta_tmp2(:) = 0.d0
    ti_ene_delta_old(:) = 0.d0
    ti_ene_delta_old2(:) = 0.d0   


  end if

  return

end subroutine ti_update_lambda

!*******************************************************************************
!
! Subroutine:  ti_print_ene
!
! Description: Write out sc ene.
!              
!*******************************************************************************

subroutine ti_print_ene(iti, sc_ener_in, sc_ener_aug_in)

  use file_io_dat_mod
  use charmm_mod, only : charmm_active  
  use prmtop_dat_mod, only: cmap_term_count

  implicit none

! Formal arguments:
  integer, intent(in)             :: iti
  ! We input the arrays so that this routine can be used to print avg/rmsd
  double precision, intent(in) :: sc_ener_in(si_cnt)
  double precision, intent(in) :: sc_ener_aug_in(ti_ene_aug_cnt)

! Local variables:
  double precision     :: sc_temp_val
  double precision     :: temp
  double precision     :: etot
  integer              :: sc_atm_cnt

  sc_atm_cnt = ti_ti_atm_cnt(iti) - ti_latm_cnt(iti)    
  if (   sc_atm_cnt  .ne. 0 &
   .and. ti_mode     .ne. 1 &
   .and. emil_sc_lcl .ne. 1 ) then 
    
     sc_temp_val = sc_ener_in(si_kin_ene) * ti_sc_fac(iti) 
    
     write (mdout,900)  sc_atm_cnt, sc_temp_val

     write (mdout,950)  sc_ener_in(si_tot_ene), &
                        sc_ener_in(si_kin_ene), &
                        sc_ener_in(si_pot_ene)
     write (mdout,1000) sc_ener_in(si_bond_ene), &
                        sc_ener_in(si_angle_ene), &
                        sc_ener_in(si_dihedral_ene)
     write (mdout,1100) sc_ener_in(si_vdw_14_ene), &
                        sc_ener_in(si_elect_14_ene), &
                        sc_ener_in(si_vdw_ene)
     write (mdout,1200) sc_ener_in(si_elect_ene)

     if(charmm_active) then
       write (mdout,1400) sc_ener_in(si_angle_ub_ene), &
                          sc_ener_in(si_dihedral_imp_ene), &
                          sc_ener_in(si_cmap_ene) !charmm
     ! if CMAP is active
     else if(cmap_term_count > 0) then
       write(mdout,1400) 0., 0., &
                         sc_ener_in(si_cmap_ene) 
     end if

     write (mdout,1300) sc_ener_aug_in(ti_rest_dist_ene), &
                        sc_ener_aug_in(ti_rest_ang_ene), &
                        sc_ener_aug_in(ti_rest_tor_ene)

     write (mdout,1600) sc_ener_aug_in(ti_elect_der_dvdl), &
                        sc_ener_aug_in(ti_vdw_der_dvdl), &
                        sc_ener_aug_in(ti_der_term)
     write (mdout,2000)
  end if
   
  if (abs(sc_ener_aug_in(ti_rest_der_term))>1e-10) then
    write (mdout,1700) sc_ener_aug_in(ti_rest_der_term)  
    write (mdout,1750) sc_ener_aug_in(ti_rest_dist_der_dvdl), &
                        sc_ener_aug_in(ti_rest_ang_der_dvdl), &
                        sc_ener_aug_in(ti_rest_tor_der_dvdl)
    write (mdout,1751) sc_ener_aug_in(ti_rest_rmsd_der_dvdl)
  endif
   

900 format (2x,'Softcore part of the system: ',i7,' atoms,',7x,&
               'TEMP(K)    = ',f14.2)
    
!! two more digits for GTI on Windows
#if ( defined(_WIN32) && defined(_DEBUG) )  
   950 format (1x,'SC_Etot=    ',f14.6,2x,'SC_EKtot=    ',f14.6,2x,&
               'SC_EPtot   =    ',f14.6)
  1000 format (1x,'SC_BOND=    ',f14.6,2x,'SC_ANGLE=    ',f14.6,2x,&
               'SC_DIHED   =    ',f14.6)
  1100 format (1x,'SC_14NB=    ',f14.6,2x,'SC_14EEL=    ',f14.6,2x,&
               'SC_VDW     =    ',f14.6)
  1200 format (1x,'SC_EEL =    ',f14.6)
  1400 format (1x,'SC_UBANG=   ',f14.6,2x,'SC_DIH_IMP=  ',f14.6,2x,&
               'SC_CMAP=        ',f14.6)
  1300 format (1x,'SC_RES_DIST=',f14.6,2x,'SC_RES_ANG=  ',f14.6,2x,&
               'SC_RES_TORS=    ',f14.6)
  1600 format (1x,'SC_EEL_DER= ',f14.6,2x,'SC_VDW_DER=  ',f14.6,2x,&
               'SC_DERIV   =    ',f14.6)

  1700 format (1x,'Lambda-dependent restraint contribution to dV/dL: Total = ' &
      ,f14.6)
  1750  format (3x,'Dist = ', f14.6, 2x,'Angl =  ',f14.6,2x,&
      'Tors = ',f14.6, 2x)
  1751 format (3x,'RMSD = ' ,f14.6)   
         
#else
   950 format (1x,'SC_Etot=    ',f11.4,2x,'SC_EKtot=    ',f11.4,2x,&
               'SC_EPtot   =    ',f11.4)
  1000 format (1x,'SC_BOND=    ',f11.4,2x,'SC_ANGLE=    ',f11.4,2x,&
               'SC_DIHED   =    ',f11.4)
  1100 format (1x,'SC_14NB=    ',f11.4,2x,'SC_14EEL=    ',f11.4,2x,&
               'SC_VDW     =    ',f11.4)
  1200 format (1x,'SC_EEL =    ',f11.4)
  1400 format (1x,'SC_UBANG=   ',f11.4,2x,'SC_DIH_IMP=  ',f11.4,2x,&
               'SC_CMAP=        ',f11.4)
  1300 format (1x,'SC_RES_DIST=',f11.4,2x,'SC_RES_ANG=  ',f11.4,2x,&
               'SC_RES_TORS=    ',f11.4)
  1600 format (1x,'SC_EEL_DER= ',f11.4,2x,'SC_VDW_DER=  ',f11.4,2x,&
               'SC_DERIV   =    ',f11.4)
               
  1700 format (1x,'Lambda-dependent restraint contribution to dV/Dl: Total = ' &
      ,f11.4)
  1750 format (3x,'Dist = ',f11.4, 2x,'Angl =  ',f11.4,2x,&
      'Tors = ',f11.4, 2x)
  1751 format (3x,'RMSD = ' ,f11.4)   
#endif
  2000 format (t2,78('-'),/)
      
  return

end subroutine ti_print_ene

!*******************************************************************************
!
! Subroutine:  ti_others_print_ene
!
! Deti_othersription: Write out ti_others ene.
!              
!*******************************************************************************

subroutine ti_others_print_ene(iti, ti_others_ener_in)

  use file_io_dat_mod
  use charmm_mod, only : charmm_active  
  use prmtop_dat_mod, only: cmap_term_count

  implicit none

! Formal arguments:
  integer, intent(in)             :: iti
  ! We input the arrays so that this routine can be used to print avg/rmsd
  double precision, intent(in) :: ti_others_ener_in(si_cnt)

! Local variables:
  integer              :: sc_atm_cnt

  sc_atm_cnt = ti_ti_atm_cnt(iti) - ti_latm_cnt(iti)    
  if (   sc_atm_cnt  .ne. 0 &
   .and. ti_mode     .ne. 1 &
   .and. emil_sc_lcl .ne. 1 ) then 
    
     write (mdout,950) 0., 0., &
                        ti_others_ener_in(si_pot_ene)
     write (mdout,1000) ti_others_ener_in(si_bond_ene), &
                        ti_others_ener_in(si_angle_ene), &
                        ti_others_ener_in(si_dihedral_ene)
     write (mdout,1100) ti_others_ener_in(si_vdw_14_ene), &
                        ti_others_ener_in(si_elect_14_ene), &
                        ti_others_ener_in(si_vdw_ene)
     write (mdout,1200) ti_others_ener_in(si_elect_ene)

     if(charmm_active) then
       write (mdout,1400) ti_others_ener_in(si_angle_ub_ene), &
                          ti_others_ener_in(si_dihedral_imp_ene), &
                          ti_others_ener_in(si_cmap_ene) !charmm
     ! if CMAP is active
     else if(cmap_term_count > 0) then
       write(mdout,1400) 0., 0., &
                         ti_others_ener_in(si_cmap_ene) 
     end if

     write (mdout,2000)
  end if

!! two more digits for GTI on Windows
#if ( defined(_WIN32) && defined(_DEBUG) )  
   950 format (1x,'TI_Etot=    ',f14.6,2x,'TI_EKtot=    ',f14.6,2x,&
               'TI_EPtot   =    ',f14.6)
  1000 format (1x,'TI_BOND=    ',f14.6,2x,'TI_ANGLE=    ',f14.6,2x,&
               'TI_DIHED   =    ',f14.6)
  1100 format (1x,'TI_14NB=    ',f14.6,2x,'TI_14EEL=    ',f14.6,2x,&
               'TI_VDW     =    ',f14.6)
  1200 format (1x,'TI_EEL =    ',f14.6)
  1400 format (1x,'TI_UBANG=   ',f14.6,2x,'TI_DIH_IMP=  ',f14.6,2x,&
               'TI_CMAP=        ',f14.6)
#else
   950 format (1x,'TI_Etot=    ',f11.4,2x,'TI_EKtot=    ',f11.4,2x,&
               'TI_EPtot   =    ',f11.4)
  1000 format (1x,'TI_BOND=    ',f11.4,2x,'TI_ANGLE=    ',f11.4,2x,&
               'TI_DIHED   =    ',f11.4)
  1100 format (1x,'TI_14NB=    ',f11.4,2x,'TI_14EEL=    ',f11.4,2x,&
               'TI_VDW     =    ',f11.4)
  1200 format (1x,'TI_EEL =    ',f11.4)
  1400 format (1x,'TI_UBANG=   ',f11.4,2x,'TI_DIH_IMP=  ',f11.4,2x,&
               'TI_CMAP=        ',f11.4)
#endif
  2000 format (t2,78('-'),/)
      
  return

end subroutine ti_others_print_ene

!*******************************************************************************
!
! Subroutine:  ti_print_mbar_ene
!
! Description: Write out mbar energies.
!              
!*******************************************************************************

subroutine ti_print_mbar_ene(etot)
  use file_io_dat_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:
  double precision, intent(in)     :: etot

! Local variables:
  integer                          :: i, j, k, l
  double precision                 :: energy   

  
  if (ifmbar .ne. 0 .and. do_mbar .and. master .and. &
  (sams_verbose.ge.1 .or. ifsams.eq.0) ) then    
     j=1
     k=bar_states
     if (sams_type.eq.3) then
       l=int(clambda*(bar_states-1)+0.0001)+1
       j=max(1,l-1)
       k=min(bar_states, l+1)
     endif
     write (mdout,'(/,a)') 'MBAR Energy analysis:'
     do i = j, k
        if (ifmbar_net .ne. 0) then
            energy = bar_cont(i)
        else
        energy = etot + bar_cont(i)
        end if
        if(gti_lam_sch .eq. 1) then
            write (mdout,'(a,f6.4,a,f18.8)') 'Energy at ', &
                mbar_lambda(i), ' = ', energy
        else
            write (mdout,'(a,f6.4,a,f18.8)') 'Energy at ', &
                bar_lambda(2,i,1), ' = ', energy
        end if
     end do
     write (mdout,2000)
  end if
 
  2000 format (t2,78('-'),/) 
  return

end subroutine ti_print_mbar_ene


!*******************************************************************************
!
! Subroutine:  ti_print_dvdl_values
!
! Description: Write out mbar energies.
!              
!*******************************************************************************

subroutine ti_print_dvdl_values
  use file_io_dat_mod

  implicit none

! Local variables:
  integer                          :: i

  write (mdout,'(a,i8,a)') 'Summary of dvdl values over ',ti_dvdl_pnt,' steps:'
  do i = 1, ti_dvdl_pnt
     write (mdout,'(f11.4)') ti_dvdl_values(i)
  end do
  write (mdout,'(a)') 'End of dvdl summary'
  write (mdout,2000)
 
  2000 format (t2,78('-'),/)
  return

end subroutine ti_print_dvdl_values

!*******************************************************************************
!
! Subroutine:  ti_cleanup
!
! Description: Clean up allocated arrays
!
!*******************************************************************************

subroutine ti_cleanup

  implicit none

! Local variables:
  integer ier_alloc
  
  if(allocated(ti_lst)) deallocate (ti_lst,stat=ier_alloc)
  if(allocated(ti_atm_lst)) deallocate (ti_atm_lst,stat=ier_alloc)
  if(allocated(ti_foureps)) deallocate (ti_foureps,stat=ier_alloc)
  if(allocated(ti_sigma6)) deallocate (ti_sigma6,stat=ier_alloc)
  if(allocated(ti_nb_frc)) deallocate (ti_nb_frc,stat=ier_alloc)
  if(allocated(ti_img_frc)) deallocate (ti_img_frc,stat=ier_alloc)
  if(allocated(ti_vir)) deallocate (ti_vir,stat=ier_alloc)
  if(allocated(bar_lambda)) deallocate (bar_lambda,stat=ier_alloc)
  if(allocated(bar_cont)) deallocate (bar_cont,stat=ier_alloc)
#if defined(CUDA)
  if(allocated(ti_lst_repacked)) deallocate (ti_lst_repacked,stat=ier_alloc)
  if(allocated(ti_latm_lst_repacked)) deallocate (ti_latm_lst_repacked,stat=ier_alloc)
#endif

  return

  end subroutine ti_cleanup
  
!*******************************************************************************
!
! Subroutine:  ti_mergeReservoir
!
! Description: Add the missing SC atoms to the recervoir structure of the real state  
!
!*******************************************************************************

subroutine ti_mergeReservoir(atm_cnt, crd, vel, box, random_value, tiRegion)


  use reservoir_mod
  use prmtop_dat_mod, only: atm_mass
  use mdin_ctrl_dat_mod, only: gti_resv_fit, gti_resv_anch, &
  gti_resv_mapmode, gti_resv_velmode, vel_scale_t
  use gbl_constants_mod, only : KB
   
  implicit none

  integer, intent(in) ::atm_cnt
  integer, intent(in) ::tiRegion
  double precision, intent(inout) :: crd(3,atm_cnt), vel(3,atm_cnt), box(3)
  double precision, intent(in) :: random_value
  
#ifdef MPI  
  
  logical,save::firstTime(2)=.true.
  integer,save::nearestRealIndex(3, 2)
  integer,save::nearestDummyIndex(3, 2)
  integer,save::zmatIndex(6)
  integer,save::softCoreList(500, 2)  
  integer,save::softCoreCount(2), rremd_atm_cnt(2)

  double precision, allocatable, save :: softCoreCrd(:,:), softCoreVel(:,:)
  integer, allocatable, save :: indexShift(:,:)
  
  integer::i,j,k,l,m,foundCount
  integer::tempIndex(3)
  double precision:: nearestRealCrd(3,3)
  double precision:: nearestRealReservoirCrd(3,3)
  double precision:: t0, tt, ss, pp(2, 3), qq(2, 3), yy(2,3), minDist(3)
  double precision:: h(3,3), hp(3,3), u(3,3), s(3), is(3), vt(3,3), work(1200)
  double precision:: nearestDummyCrd(3,3)
  double precision:: tmat(3,3), vrot(3,3)
  double precision:: zmat(3,6), zmat_crd(6,3), zmat_copy(6,3)

  integer::info, iwork(200)
  
  double precision:: translation(3)
  double precision:: rotation(3,3)  

  double precision, allocatable, save::ti_vel(:,:)

  if (.not.allocated(indexShift))  allocate(indexShift(2,atm_cnt))
  
  if ( firstTime(tiRegion) ) then
    
    allocate(softCoreCrd(3,500), softCoreVel(3,500))

    if (gti_resv_velmode .eq. 3 ) allocate(ti_vel(3, count(ti_atm_lst(3-tiRegion,:)/=0)))
           
    softCoreList(:, tiRegion )=0
    softCoreCount(tiRegion)=0
    indexShift(tiRegion,:)=0

    j=0
    do i=1, atm_cnt
      if (ti_lst(3-tiRegion,i)>0) then
        j=j+1
        if(ti_sc_lst(i)>0) then
             softCoreCount(tiRegion)=softCoreCount(tiRegion)+1
             softCoreList(softCoreCount(tiRegion), tiRegion)=i
             softCoreCrd(1:3,softCoreCount(tiRegion))=crd(1:3,i)
        endif
      endif
      indexShift(tiRegion,i)=j
    end do
    softCoreList(softCoreCount(tiRegion)+1, tiRegion)=atm_cnt
    rremd_atm_cnt(tiRegion) = atm_cnt-softCoreCount(tiRegion)-ti_latm_cnt(tiRegion)
    
    nearestRealIndex(:, tiRegion)=0
    nearestDummyIndex(:, tiRegion)=0
    
    if(gti_resv_anch(1).gt.0 .and. gti_resv_anch(1).le. atm_cnt &
      .and. gti_resv_anch(2).gt.0 .and. gti_resv_anch(2).le. atm_cnt &
      .and. gti_resv_anch(3).gt.0 .and. gti_resv_anch(3).le. atm_cnt & 
      .and. gti_resv_anch(1).ne.gti_resv_anch(2) &
      .and. gti_resv_anch(2).ne.gti_resv_anch(3) &
      .and. gti_resv_anch(1).ne.gti_resv_anch(3) ) then      
      
        nearestRealIndex(1:3, tiRegion)= gti_resv_anch(1:3) 
    
    else    
        
    !!find out the nearest real state atoms to the dummy atoms
   
      minDist(1:3)=1.0E9
      foundCount=1
      do i=1, atm_cnt  
        if (ti_lst(3-tiRegion, i).gt.0 .and. ti_sc_lst(i).eq.0 .and. atm_mass(i)>3.5 ) then
      
          t0=1.0E9
          do j=1, softCoreCount(tiRegion)
              k=softCoreList(j,tiRegion)
              tt=(crd(1, i)-crd(1, k))**2 + (crd(2, i)-crd(2, k))**2 +(crd(3, i)-crd(3, k))**2
              if (tt<t0) t0=tt
          end do
        
            
          if ( t0<minDist(1) ) then
            nearestRealIndex(3, tiRegion)=nearestRealIndex(2, tiRegion)
            nearestRealIndex(2, tiRegion)=nearestRealIndex(1, tiRegion)
            nearestRealIndex(1, tiRegion)=i
            minDist(1)=t0
          else if  ( t0<minDist(2) ) then
            nearestRealIndex(3, tiRegion)=nearestRealIndex(2, tiRegion)
            nearestRealIndex(2, tiRegion)=i
            minDist(2)=t0
          else  if  ( t0<minDist(3) ) then
            nearestRealIndex(3, tiRegion)=i
            minDist(3)=t0
          endif
           
        endif
      end do

      minDist(1:3)=1.0E9
    
      if (gti_resv_fit .gt. 0 .and. sum(ti_sc_bat_lst(1:atm_cnt)) .ne. 0) then
         !!Find the bonded atom
        tempIndex=0
         do i=1, atm_cnt  
           if (atm_mass(i)>3.5 .and. ti_lst(3-tiRegion, i).gt.0 .and. ti_sc_lst(i).eq.0 ) then
              if ( iand( ti_sc_bat_lst(i) , (2**(7-3*tiRegion)) ) .gt. 0)  then
                 tempIndex(1)=i
                 cycle
              endif
           endif  
         enddo
        do i=1, atm_cnt  
           if (atm_mass(i)>3.5 .and.  ti_lst(3-tiRegion, i).gt.0 .and. ti_sc_lst(i).eq.0 .and. &
             tempIndex(1) .ne. i) then
              if ( iand( ti_sc_bat_lst(i) , (2**(8-3*tiRegion)) ) .gt. 0)  then
                 tempIndex(2)=i
                 cycle
              endif
           endif  
        enddo
        do i=1, atm_cnt  
           if (atm_mass(i)>3.5 .and.  ti_lst(3-tiRegion, i).gt.0 .and. ti_sc_lst(i).eq.0 .and. &
             tempIndex(1) .ne. i .and. tempIndex(2) .ne. i) then
              if ( iand( ti_sc_bat_lst(i) , (2**(9-3*tiRegion)) ) .gt. 0)  then
                 tempIndex(3)=i
                 cycle
              endif
           endif  
        enddo     
      
        do i=1, 3
           if(tempIndex(1).eq.0) then
             if (nearestRealIndex(i, tiRegion) .ne. tempIndex(2) .and. &
                  nearestRealIndex(i, tiRegion) .ne. tempIndex(3) ) then
                tempIndex(1)=nearestRealIndex(i, tiRegion) 
               endif
           endif
           if(tempIndex(2).eq.0) then
             if (nearestRealIndex(i, tiRegion) .ne. tempIndex(1) .and. &
                  nearestRealIndex(i, tiRegion) .ne. tempIndex(3) ) then
                tempIndex(2)=nearestRealIndex(i, tiRegion) 
               endif
           endif
           if(tempIndex(3).eq.0) then
             if (nearestRealIndex(i, tiRegion) .ne. tempIndex(1) .and. &
                  nearestRealIndex(i, tiRegion) .ne. tempIndex(2) ) then
                tempIndex(3)=nearestRealIndex(i, tiRegion) 
               endif
           endif
        enddo
        nearestRealIndex(1:3, tiRegion)=tempIndex(1:3)

      endif

      minDist(1:3)=1.0E9 

      ! Find Dummy atoms nearest to real atoms
      do j=1, softCoreCount(tiRegion)
        k=softCoreList(j,tiRegion)
        t0=1.0E9
        l=nearestRealIndex(1, tiRegion)
        tt=(crd(1, l)-crd(1, k))**2 + (crd(2, l)-crd(2, k))**2 +(crd(3, l)-crd(3, k))**2
        if (tt<t0) t0=tt
      
        if ( t0<minDist(1) ) then
          nearestDummyIndex(3, tiRegion)=nearestDummyIndex(2, tiRegion)
          nearestDummyIndex(2, tiRegion)=nearestDummyIndex(1, tiRegion)
          nearestDummyIndex(1, tiRegion)=k
          minDist(1)=t0
          zmatIndex(6) = zmatIndex(5)
          zmatIndex(5) = zmatIndex(4)
          zmatIndex(4) = j
        else if  ( t0<minDist(2) ) then
          nearestDummyIndex(3, tiRegion)=nearestDummyIndex(2, tiRegion)
          nearestDummyIndex(2, tiRegion)=k
          minDist(2)=t0
          zmatIndex(6) = zmatIndex(5)
          zmatIndex(5) = j
        else  if  ( t0<minDist(3) ) then
          nearestDummyIndex(3, tiRegion)=k
          minDist(3)=t0
          zmatIndex(6) = j
        endif
      end do
      do j=1,3
        zmatIndex(j) = nearestRealIndex(4-j,tiRegion)
      end do

    endif !! gti_resv_anch check
   
    firstTime(tiRegion)=.false.      
    
  endif !! first time check

  do i=1, 3
    nearestRealCrd(1:3,i)=crd(1:3, nearestRealIndex(i, tiRegion) )
    nearestDummyCrd(1:3,i)=crd(1:3, nearestDummyIndex(i, tiRegion) )
  end do
 
  if ( gti_resv_mapmode .eq. 1 ) then

    ! Make Z matrix
    zmat = 0
    zmat_crd = 0

    do i=1,3
      zmat_crd(i,1:3)=nearestRealCrd(1:3, 4-i)
      zmat_crd(3+i,1:3)=nearestDummyCrd(1:3, i)
    end do

    do i=2,6
      zmat(1,i) = sqrt(sum((zmat_crd(i,1:3)-zmat_crd((i-1),1:3))**2))
    if ( i .ge. 3 ) then 
      zmat(2,i) = CalculateAngle(zmat_crd(i,1:3), zmat_crd((i-1),1:3), zmat_crd((i-2),1:3))
    end if
    if ( i .ge. 4 ) then
      zmat(3,i) = CalculateDihedralAngle(zmat_crd(i,1:3), zmat_crd((i-1),1:3), &
        zmat_crd((i-2),1:3), zmat_crd((i-3),1:3))
    end if
    end do
  end if  !! gti_resv_mapmode .eq. 1

  do i=1, softCoreCount(tiRegion)
    j=softCoreList(i, (tiRegion))
    softCoreCrd(1:3,i)=crd(1:3,j)
    softCoreVel(1:3,i)=vel(1:3,j)
  enddo

  if (gti_resv_velmode .eq. 3 ) then 
    do i=1,count(ti_atm_lst(3-tiRegion,:)/=0)
      j = ti_atm_lst(3-tiRegion,i)
      ti_vel(1:3,i)=vel(1:3,j)
    end do 
  end if !! gti_resv_velmode .eq. 3 
   
  rremd_idx(tiRegion) = min(reservoir_size(tiRegion),  &
    int(random_value * reservoir_size(tiRegion)) + 1)

  call load_reservoir_structure(crd, vel, rremd_atm_cnt(tiRegion), box(1), box(2), box(3), tiRegion) 
  
  do i=atm_cnt-1, 1, -1
    j=indexShift(tiRegion,i)
    if (j .gt. 0) then
      crd(1:3,i+1)=crd(1:3,i-j+1)
      crd(1:3,i-j+1)=0.0  
      if (reservoir_velocity(tiRegion)) then
        vel(1:3,i+1)=vel(1:3,i-j+1)
        vel(1:3,i-j+1)=0.0  
      end if  
    end if
  end do
  
  do i=1, ti_latm_cnt(tiRegion)
     j = ti_latm_lst(3-tiRegion,i)
     k = ti_latm_lst(tiRegion,i)
     crd(1:3,j)=crd(1:3,k)
     if (reservoir_velocity(tiRegion)) then
        vel(1:3,j)=vel(1:3,k) 
     end if  
  end do
  
  do i=1, 3
    nearestRealReservoirCrd(1:3,i)=crd(1:3, nearestRealIndex(i, tiRegion) )
  end do

  if ( gti_resv_mapmode .eq. 1 ) then
    ! Use Z-Matrix to get location of dummy atoms
    ! Also velocities
    zmat_copy = zmat_crd
    do i=1,3
      zmat_crd(i,1:3)=nearestRealReservoirCrd(1:3, 4-i)
    end do

    do i=4,MIN(6,(3+softCoreCount(tiRegion)))
      j = nearestDummyIndex(i-3, tiRegion)
      ! angle = 180 degrees is a special case 
      if ( zmat(2,i) .eq. acos(-1.0) ) then
        zmat_crd(i,1:3) = CalcCrdLinear(zmat_crd, zmat, i)
      else
        zmat_crd(i,1:3) = CalcCrd(zmat_crd, zmat, i)
      end if
      crd(1:3,j) = zmat_crd(i,1:3) 
      ! rotate velocities
      vrot = velocity_rotation(nearestDummyCrd, zmat_crd, zmat_copy, i)
      !vel(1:3, j) = MATMUL(vrot,softCoreVel(1:3,i-3)) 
      vel(1:3, j) = MATMUL(vrot,softCoreVel(1:3,zmatIndex(i)))

    end do

    ! If Softcore region is greater than 3 atoms, use rotation matrix 
    ! to find location of remaining softcore atoms
    if ( softCoreCount(tiRegion) .gt. 3 ) then
      !nearestDummyCrd(1:3,i)=crd(1:3, nearestDummyIndex(i, tiRegion) )
      ! reset pp, pp=qq since zmatrix is used
      pp = 0
      ! reset other variables used in svd
      s=0
      u=0
      vt=0
      do i=1, 2
        pp(i, 1:3)=nearestDummyCrd(1:3,i+1)-nearestDummyCrd(1:3,1)
      end do
      qq=pp
      ! Unsure if rodrigues rotation needed; TBD
      h=matmul(transpose(pp),qq)
      info=0
      call dgesvd('A', 'A', 3, 3, h, 3, s, u, 3, vt, 3 , work, 100, info)

      tmat=0
      tmat(1,1)=1
      tmat(2,2)=1
      tmat(3,3)=MatrixDeterminant(matmul(u,vt))
    
      rotation=matmul(u,matmul(tmat,vt))

      ! apply these transformations to the additional dummy atoms 
      do i=1, softCoreCount(tiRegion)
        j=softCoreList(i, tiRegion)
        if (.not.( ANY(nearestDummyIndex .eq. j) ) ) then  
          s(1:3)=softCoreCrd(1:3,i)-nearestDummyCrd(1:3,1)
          do k=1, 3
            crd(k, j)=s(1)*rotation(1,k)+s(2)*rotation(2,k)+s(3)*rotation(3,k)
          end do
          crd(1:3, j)=crd(1:3, j)+zmat_crd(1:3,4)
        end if  
      enddo

      !! Copy the velocities of the dummy atoms over and then 
      !! rotate and rescale them. 
      if(reservoir_velocity(tiRegion)) then
        tt=0.0
        do i=1, softCoreCount(tiRegion)
          j=softCoreList(i, tiRegion)
          if (.not.( ANY(nearestDummyIndex .eq. j) ) ) then
            is(1:3)=softCoreVel(1:3,i)
            do k=1, 3
              vel(k,j)=is(1)*rotation(1,k)+is(2)*rotation(2,k)+is(3)*rotation(3,k)
            end do
          end if ! not in nearestDummyIndex
        end do ! is number of softcore atoms > 3
      end if ! reservoir velocities
    end if ! softCoreCount > 3

  else ! ( gti_resv_mapmode .ne. 1 )

    !! SVD rotation matrix based on 3 CC atoms 

    !! Get the translation/rotation matrice needed to map the current 3-atom 
    !! to the corresponding recervoir atoms
    do i=1, 2
      pp(i, 1:3)=nearestRealCrd(1:3,i+1)-nearestRealCrd(1:3,1)
    end do
    
    do i=1, 2
      qq(i, 1:3)=nearestRealReservoirCrd(1:3,i+1)-nearestRealReservoirCrd(1:3,1)
    end do
      
    !!if (gti_resv_fit .gt. 0 .and. sum(ti_sc_bat_lst(1:atm_cnt)) .ne. 0) then
    if (gti_resv_fit .gt. 0) then
      !! scale the connection bond
      do i=1,2 
        s(i)=sqrt(sum(pp(i,1:3)**2))
        is(i)=sqrt(sum(qq(i,1:3)**2))
        pp(i,1:3)=is(i)/s(i)*pp(i,1:3)
      end do
      
      !!Rodrigues' Rotation 
      t0=acos(sum(qq(1,1:3)*qq(2,1:3))/is(1)/is(2))    
      tt=acos(sum(pp(1,1:3)*pp(2,1:3))/is(1)/is(2))
      tt=t0-tt
      
      s(1)=pp(1,2)*pp(2,3)-pp(1,3)*pp(2,2)
      s(2)=pp(1,3)*pp(2,1)-pp(1,1)*pp(2,3)
      s(3)=pp(1,1)*pp(2,2)-pp(1,2)*pp(2,1)
      ss=sqrt(sum(s(1:3)**2))
      s(1:3)=s(1:3)/ss
      
      u=0
      u(1,2)=-s(3)
      u(1,3)=s(2)
      u(2,1)=s(3)
      u(2,3)=-s(1)
      u(3,1)=-s(2)
      u(3,2)=s(1)
      vt=0
      vt(1,1)=1
      vt(2,2)=1
      vt(3,3)=1
      vt=vt+sin(-tt)*u+(1.0-cos(-tt))*matmul(u,u)
      pp(2,1:3)=matmul(pp(2,1:3),vt)

      !!debug:  tt=acos(sum(pp(1,1:3)*pp(2,1:3))/sqrt(sum(pp(1,1:3)**2))/sqrt(sum(pp(2,1:3)**2)))
    endif
   
    h=matmul(transpose(pp),qq)
    info=0
    call dgesvd('A', 'A', 3, 3, h, 3, s, u, 3, vt, 3 , work, 100, info)

    tmat=0
    tmat(1,1)=1
    tmat(2,2)=1
    tmat(3,3)=MatrixDeterminant(matmul(u,vt))

    rotation=matmul(u,matmul(tmat,vt)) 
    !!debug: yy=matmul(pp,rotation)

    !! apply these transformations to the additional dummy atoms 
    do i=1, softCoreCount(tiRegion)
      j=softCoreList(i, tiRegion)

      s(1:3)=softCoreCrd(1:3,i)-nearestRealCrd(1:3,1)
      do k=1, 3
        crd(k, j)=s(1)*rotation(1,k)+s(2)*rotation(2,k)+s(3)*rotation(3,k)
      end do
     
      crd(1:3, j)=crd(1:3, j)+nearestRealReservoirCrd(1:3,1)
    enddo
    !! Copy the velocities of the dummy atoms over and then 
    !! rotate them. 
    do i=1, softCoreCount(tiRegion)
      j=softCoreList(i, tiRegion)
      is(1:3)=softCoreVel(1:3,i)
      do k=1, 3
        vel(k,j)=is(1)*rotation(1,k)+is(2)*rotation(2,k)+is(3)*rotation(3,k)
      end do
    end do

  end if ! gti_resv_mapmode
  
  !! Copy the velocities of the dummy atoms over and then 
  !! rescale them. 

  if(reservoir_velocity(tiRegion)) then
    tt=0.0
    do i=1, softCoreCount(tiRegion)
      j=softCoreList(i, tiRegion)
      is(1:3)=softCoreVel(1:3,i)
      ss=is(1)*is(1)+ is(2)*is(2) + is(3)*is(3)
      tt=tt+ss*atm_mass(j)
    end do
    if (gti_resv_velmode .eq. 3) then
      t0=reservoir_temperature(tiRegion)
      if (t0<10.0 .or. t0>1000) t0=300.0 !! in case that the temperature was set wrong
      tt = 0.0 
      do i=1,count(ti_atm_lst(3-tiRegion,:)/=0) 
        j = ti_atm_lst(3-tiRegion,i)
        is(1:3) = ti_vel(1:3,i)
        ss=is(1)*is(1)+ is(2)*is(2) + is(3)*is(3)
        tt=tt+ss*atm_mass(j)
      end do
      tt=tt/(count(ti_atm_lst(3-tiRegion,:)/=0)*3.0)/KB
      ss = sqrt(tt/t0)
      do i=1,count(ti_atm_lst(tiRegion,:)/=0) 
        j = ti_atm_lst(tiRegion,i)
        vel(1:3,j)=vel(1:3,j)*ss
      end do
    else if (gti_resv_velmode .ne. 1) then ! remain modes scale dummy atom velocities 
      tt=tt/(softCoreCount(tiRegion)*3.0)/KB
      if (gti_resv_velmode .eq. 2) then 
        t0 = vel_scale_t
      else 
        t0=reservoir_temperature(tiRegion)
      end if                  
      if (t0<10.0 .or. t0>1000) t0=300.0 !! in case that the temperature was set wrong
        ss=sqrt(t0/tt)
        do i=1, softCoreCount(tiRegion)
          j=softCoreList(i, tiRegion)
          vel(1:3,j)=vel(1:3,j)*ss
        end do
      end if
  end if

#endif  
end subroutine ti_mergeReservoir

function CalculateAngle(point1, point2, point3) result(angle)
  implicit none

  ! Declare variables
  real(8) :: point1(3), point2(3), point3(3)
  real(8) :: vector1(3), vector2(3)
  real(8) :: dotProduct, magnitude1, magnitude2, angle

  ! Calculate vectors
  vector1 = point2 - point1
  vector2 = point3 - point2

  ! Calculate dot product
  dotProduct = dot_product(vector1, vector2)

  ! Calculate magnitudes
  magnitude1 = sqrt(dot_product(vector1, vector1))
  magnitude2 = sqrt(dot_product(vector2, vector2))

  ! Calculate angle in radians
  angle = acos(dotProduct / (magnitude1 * magnitude2))

  ! Convert angle to degrees
  !angle = angle * 180.0 / acos(-1.0)
  angle = acos(-1.0)-angle

end function CalculateAngle

FUNCTION CalculateDihedralAngle(a, b, c, d) result(angle)
  implicit none
  REAL(8), INTENT(IN) :: a(3), b(3), c(3), d(3)
  real(8) :: angle
  real(8) :: b0(3), b1(3), b2(3), v(3), w(3)
  real(8) :: x, y

  ! Extract Cartesian points
  b0 = a - b 
  b1 = c - b
  b2 = d - c 

  ! Normalize b1
  b1 = b1 / sqrt(dot_product(b1, b1))

  ! Vector rejections
  v = b0 - dot_product(b0, b1) * b1
  w = b2 - dot_product(b2, b1) * b1

  ! Compute the torsion angle
  x = dot_product(v, w)
  y = dot_product(cross_product(b1, v), w)
  angle = atan2(y, x)

  ! Convert angle to degrees
  !angle = angle !* 180.0 / acos(-1.0)

END FUNCTION CalculateDihedralAngle

! Helper function for cross product
FUNCTION cross_product(a, b)
  implicit none
  REAL(8), INTENT(IN) :: a(3), b(3)
  REAL(8) :: result(3)
  REAL(8) :: cross_product(3)
  result(1) = a(2) * b(3) - a(3) * b(2)
  result(2) = a(3) * b(1) - a(1) * b(3)
  result(3) = a(1) * b(2) - a(2) * b(1)
  cross_product = result
END FUNCTION cross_product

! Function to calculate determinant for a 3x3 matrix
function MatrixDeterminant(mat) result(determinant)
  real(8), intent(in) :: mat(3, 3)
  real(8) :: determinant

  determinant = mat(1,1)*(mat(2,2)*mat(3,3) - mat(3,2)*mat(2,3)) - &
                mat(1,2)*(mat(2,1)*mat(3,3) - mat(3,1)*mat(2,3)) + &
                mat(1,3)*(mat(2,1)*mat(3,2) - mat(3,1)*mat(2,2))

end function MatrixDeterminant

function CalcCrdLinear(zmat_crd, zmat, i) result(new_crd)
  integer, intent(in) :: i
  real(8), intent(in) :: zmat_crd(6,3), zmat(3,6)
  real(8) :: new_crd(3), v32(3), m32
  v32 = zmat_crd(i-1,1:3) - zmat_crd(i-2,1:3) 
  m32 = SQRT(SUM(v32**2))
  new_crd = zmat_crd(i-1,1:3) + zmat(1,i)*(v32/m32)
end function CalcCrdLinear

function CalcCrd(zmat_crd, zmat, i) result(new_crd)
  integer, intent(in) :: i
  real(8), intent(in) :: zmat_crd(6,3), zmat(3,6)
  real(8) :: new_crd(3) 
  real(8) :: bond_vector(1,3), axes(3,3), disp_vector(1,3)
  new_crd = 0

  !new_crd
  axes=CalcAxes(zmat_crd(i-1,1:3),zmat_crd(i-2,1:3),zmat_crd(i-3,1:3)) 
  bond_vector(1,1:3)=CalcBondVector(zmat(1:3,i)) 
  disp_vector = matmul(bond_vector, axes) !(1,1:3)
  new_crd = zmat_crd(i-1,1:3) + disp_vector(1,1:3)
end function CalcCrd

function CalcBondVector(zmat) result(bv)
  real(8), intent(in) :: zmat(3)
  real(8) :: bv(3)
  bv(1) = zmat(1) * SIN(zmat(2)) * SIN(zmat(3))
  bv(2) = zmat(1) * SIN(zmat(2)) * COS(zmat(3))
  bv(3) = zmat(1) * COS(zmat(2))
end function CalcBondVector

function CalcAxes(p1, p2, p3) result(axes)
  real(8), intent(in) :: p1(3), p2(3), p3(3)
  real(8) :: axes(3,3), m21, m32, u21(3), u32(3)
  real(8) :: u23c21(3), u21c23c21(3)
  real(8) :: v32(3), v21(3)

  ! Calculate vectors  
  v21 = p2-p1
  v32 = p3-p2

  ! Calculate magnitude of vectors
  m21 = SQRT(SUM(v21**2)) ! magnitude of vector 1
  m32 = SQRT(SUM(v32**2)) ! magnitude of vector 2

  ! Calculate unit vector (v_hat) 
  u21 = v21/m21 
  u32 = v32/m32

  ! Unit Cross Product between Unit vectors   
  u23c21 = CalcUnitCrossProduct(u32, u21)
  u21c23c21 = CalcUnitCrossProduct(u21, u23c21)

  axes(3,1:3) = u21
  axes(2,1:3) = u21c23c21
  axes(1,1:3) = CalcUnitCrossProduct(axes(2,1:3), axes(3,1:3))
end function CalcAxes

function CalcUnitCrossProduct(u1, u2) result(ucp)
  real(8), intent(in) :: u1(3), u2(3)
  real(8) :: ucp(3)
  real(8) :: cos_u1u2, sin_u1u2
  ucp = 0
  cos_u1u2 = MAX(MIN(SUM(u1*u2), 1.0), -1.0)
  sin_u1u2 = SQRT(1 - cos_u1u2**2)
  
  ucp(1) = (u1(2) * u2(3) - u1(3) * u2(2))/ sin_u1u2
  ucp(2) = (u1(3) * u2(1) - u1(1) * u2(3))/ sin_u1u2
  ucp(3) = (u1(1) * u2(2) - u1(2) * u2(1))/ sin_u1u2
end function CalcUnitCrossProduct

function velocity_rotation(nearestDummyCrd, zmat_crd, zmat_copy, i) result(rotation_matrix)
  integer, intent(in) :: i
  real(8), intent(in) :: nearestDummyCrd(3,3), zmat_crd(6,3), zmat_copy(6,3)
  real(8) :: rotation_matrix(3,3), rx(3,3), ry(3,3), rz(3,3)
  real(8) :: vo(3), vn(3), cross(3), mo, mn, mcross, ucross(3)
  real(8) :: angle, sina, cosa, u(3,3), vt(3,3), t0, tt, sn, so
  vn = zmat_crd(i,1:3) - zmat_crd(i-1,1:3)
  vo = nearestDummyCrd(1:3,i-3) - zmat_copy(i-1,1:3)

  cross = cross_product(vn,vo)

  cross = cross_product(vn,vo)
  ucross = cross / sqrt(sum(cross**2))

  mn = SQRT(SUM(vn**2))
  mo = SQRT(SUM(vo**2))
  mcross = SQRT(SUM(cross**2))
  angle = ASIN(mcross/(mo*mn))
  sina = SIN(angle)
  cosa = COS(angle)

  sn=sqrt(sum(vn**2))
  so=sqrt(sum(vo**2))

  !!Rodrigues' Rotation 
  tt=acos(sum(vo*vn)/mn/mo)
 
  cross = cross_product(vn,vo)
  ucross = cross / sqrt(sum(cross**2))
  u=0
  u(1,2)=-ucross(3)
  u(1,3)=ucross(2)
  u(2,1)=ucross(3)
  u(2,3)=-ucross(1)
  u(3,1)=-ucross(2)
  u(3,2)=ucross(1)
  vt=0
  vt(1,1)=1
  vt(2,2)=1
  vt(3,3)=1

  vt=vt+sin(-tt)*u+(1.0-cos(-tt))*matmul(u,u)
  !pp(2,1:3)=matmul(pp(2,1:3),vt)
  rotation_matrix=vt

end function velocity_rotation

end module ti_mod
