module gti_mod
#if defined(GTI)

  use gbl_constants_mod
  use mdin_ctrl_dat_mod
  use energy_records_mod
  use state_info_mod
  use ti_mod
  use sams_mod
  use iso_c_binding

!! TI control variables for both C and FORTRAN are put in gti_controlVariable.i
!! to enforce consistency
#undef  C_COMPILER
#include "gti_controlVariable.i"

  !! This section must be consistent w/ gti_simulationConst.h and ti.F90

  integer, parameter:: TypeGen = 1, TypeBAT = 2, TypeEleRec = 3, &
         TypeEleCC = 4, TypeEleSC = 5, TypeVDW = 6,  TypeRestBA = 7, &
         TypeEleSSC = 8, TypeRMSD=9, &
         TypeTotal = 9 !! FORTRAN-index


  type(TypeCtrlVar)::ctrlVar


  double precision, save, dimension(2)::specicalEstimates=0.0, SE2=0.0

  character*80, dimension(GPUPotEnergyTerms), save:: items
  integer, dimension(GPUPotEnergyTerms), save:: lambdaType

  integer, dimension(GPUPotEnergyTerms), save:: si_index, sc_index, ti_aug_index
  type energyPtr
    real(8), pointer :: ptr
  end type energyPtr
  
  type(energyPtr), dimension(GPUPotEnergyTerms), save :: pEnergy
  double precision, save, dimension(GPUPotEnergyTerms,TIEnergyBufferMultiplier):: gti_potEnergy
  double precision, dimension(2) :: self_Rec, self_CC,self_SC, vdwCorr 

  logical, save::firsttime=.true.

  private::firsttime
  private::gti_init

  contains

subroutine gti_load_control_variable()

  ctrlVar%scalpha	=	scalpha
  ctrlVar%scbeta 	=	scbeta
  ctrlVar%scgamma	=	scgamma
  ctrlVar%klambda	=	klambda
		
  ctrlVar%addSC	=	gti_add_sc
  ctrlVar%batSC	=	gti_bat_sc
  ctrlVar%autoAlpha	=	gti_auto_alpha
  ctrlVar%scaleBeta	=	gti_scale_beta
		
  ctrlVar%eleExp	 = gti_ele_exp
  ctrlVar%eleGauss	 = gti_ele_gauss
  ctrlVar%eleSC	=	gti_ele_sc
		
  ctrlVar%vdwCap	 = 	gti_vdw_cap
  ctrlVar%vdwExp	 =	 gti_vdw_exp
  ctrlVar%vdwSC	 = 	gti_vdw_sc
		
  ctrlVar%tiCut	 = 	gti_cut

  ctrlVar%cut_sc	 = 	gti_cut_sc

  if (gti_cut_sc .gt. 0) then

    if (gti_cut_sc_off .gt. 0) then
      ctrlVar%cut_sc1	=	min(gti_cut_sc_off, vdw_cutoff)
    else
      ctrlVar%cut_sc1	=	vdw_cutoff
    end if

    if (gti_cut_sc_on .gt. 0) then
      if (fswitch .gt. 0) then
        ctrlVar%cut_sc0	=	min(gti_cut_sc_on, fswitch)
      else
        ctrlVar%cut_sc0	=	gti_cut_sc_on
      endif 
    else
      ctrlVar%cut_sc0=vdw_cutoff-2.0
    endif

    if (ctrlVar%cut_sc1 .le. ctrlVar%cut_sc0) &
        ctrlVar%cut_sc0 = ctrlVar%cut_sc1-2.0

  endif
      
  write(mdout,'(a)') ''
  write(mdout,'(a)') '|--------------------------------------------------------------------------------------------'
  write(mdout,'(a)') '| Extra TI control variables'

  write(mdout,'(a, 5x,4(a,i4))') '|', &
       'gti_add_sc     =', ctrlVar%addSC, &
       ', gti_ele_gauss  =', ctrlVar%eleGauss, &
	     ', gti_auto_alpha =', ctrlVar%autoAlpha, &
       ', gti_scale_beta =', ctrlVar%scaleBeta

  write(mdout,'(a, 5x,4(a,i4))') '|', &
         'gti_ele_exp    =', ctrlVar%eleExp	,  &
       ', gti_vdw_exp    =', ctrlVar%vdwExp, &
       ', gti_ele_sc     =', ctrlVar%eleSC	,  &
       ', gti_vdw_sc     =', ctrlVar%vdwSC

  write(mdout,'(a, 5x,2(a,i4))') '|', &
         'gti_cut        =', ctrlVar%tiCut, &
       ', gti_cut_sc     =', ctrlVar%cut_sc

  write(mdout,'(a, 5x,2(a,f8.4))') '|', &
       'gti_cut_sc_on    =', ctrlVar%cut_sc0, &
       ', gti_cut_sc_off    =', ctrlVar%cut_sc1

  write(mdout,'(a)') '|--------------------------------------------------------------------------------------------'
  write(mdout,'(a)') ''


end subroutine

subroutine gti_init()

  implicit none
  integer i

  if (firsttime) then

    firsttime=.false.

!! Names
    items(:)=''
    items(2)="vDW"

    items(4)="Bond"
    items(5)="Angle"
    items(6)="Torsion"

    items(7)="EE14-CC"
    items(8)="VDW14"  

    items(9)="Constraint"  
    items(10)="Elec-Rec"
    items(11)="Elec-CC"
    
    items(12)="UB-Angle"
    items(13)="Impprop"
    items(14)="CMap"

    items(15)="Res_Bond"
    items(16)="Res_Angl"
    items(17)="Res_Tors"  
    items(18)="Res_RMSD"  
      
    items(TIExtraShift+1)="Elec-SC"
    items(TIExtraShift+2)="EE14-SC"
    items(TIExtraShift+3)="Self-Rec"
    items(TIExtraShift+4)="Self-CC"
    items(TIExtraShift+5)="Self-SC"
    items(TIExtraShift+6)="vDW-Corr"

    items(TISpecialShift+1)="Orig BAT"
    items(TISpecialShift+2)="Current BAT"

!Lambda Types: 

    lambdaType(:)=TypeGen !

    lambdaType(2)=TypeVDW !"vDW"

    lambdaType(4)=TypeBAT !"Bond"
    lambdaType(5)=TypeBAT !"Angle"
    lambdaType(6)=TypeBAT !"Torsion"

    lambdaType(7)=TypeEleCC !"1-4 Ele"
    lambdaType(8)=TypeVDW !"1-4 VDW"  

    lambdaType(9)=TypeGen !"Constraint"  
    lambdaType(10)=TypeEleRec !"Elec-Rec"
    lambdaType(11)=TypeEleCC !"Elec-Dir"

    lambdaType(12)=TypeGen !"UB-Angle"
    lambdaType(13)=TypeGen !"Impprop"
    lambdaType(14)=TypeGen !"CMap"

    lambdaType(15)=TypeRestBA !"Rest_Bond"
    lambdaType(16)=TypeRestBA !"Rest_Angle"
    lambdaType(17)=TypeRestBA !"Rest_Tors"  
    
    lambdaType(18)=TypeRMSD !"Rest_RMSD"      

    lambdaType(TIExtraShift+1)=TypeEleSC !"Elec-SC"
    lambdaType(TIExtraShift+2)=TypeEleSC !"EE14-SC"
    lambdaType(TIExtraShift+3)=TypeEleRec !"Self-Rec"
    lambdaType(TIExtraShift+4)=TypeEleRec !"Self-CC"
    lambdaType(TIExtraShift+5)=TypeEleRec !"Self-SC"
    lambdaType(TIExtraShift+6)=TypeVDW !"VDW-corr"    
    
!! Map to the linear SI index (in state_info.F90)
    si_index=-1
    si_index(2)=si_vdw_ene    
    si_index(3)=si_vdw_ene 
    si_index(4)=si_bond_ene 
    si_index(5)=si_angle_ene 
    si_index(6)=si_dihedral_ene 
    si_index(7)=si_elect_14_ene  
    si_index(8)=si_vdw_14_ene  

    si_index(9)=si_restraint_ene   
  
    si_index(10)=si_elect_ene        
    si_index(11)=si_elect_ene    
    
    si_index(12)=si_angle_ub_ene
    si_index(13)=si_dihedral_imp_ene    
    si_index(14)=si_cmap_ene

    si_index(15)=si_restraint_ene   
    si_index(16)=si_restraint_ene   
    si_index(17)=si_restraint_ene   
    si_index(18)=si_restraint_ene    
    
    si_index(TIExtraShift+1)=si_elect_ene
    si_index(TIExtraShift+2)=si_elect_14_ene
    si_index(TIExtraShift+3)=si_elect_ene
    si_index(TIExtraShift+4)=si_elect_ene
    si_index(TIExtraShift+5)=si_elect_ene
    si_index(TIExtraShift+6)=si_vdw_ene

!! Map to the SC-only terms 
    sc_index=-1
    !! SC NB terms
    sc_index(2)=si_vdw_ene
    sc_index(10)=si_elect_ene
    sc_index(11)=si_elect_ene
    
    !! SC bonded terms
    sc_index(4)=si_bond_ene
    sc_index(5)=si_angle_ene
    sc_index(6)=si_dihedral_ene    
    sc_index(7)=si_elect_14_ene    
    sc_index(8)=si_vdw_14_ene   

!! Map to the lambda-dependent terms
    ti_aug_index=-1
    !! lambda-dependent vdw & ele for SC atoms 
    ti_aug_index(2)=ti_vdw_der_dvdl      
    ti_aug_index(8)=ti_vdw_der_dvdl !"1-4 VDW"       

    ti_aug_index(7)=ti_elect_der_dvdl !"1-4 Ele"
    ti_aug_index(11)=ti_elect_der_dvdl !!EE-CC
    ti_aug_index(TIExtraShift+1)=ti_elect_der_dvdl !!EE-SC
    ti_aug_index(TIExtraShift+2)=ti_elect_der_dvdl !!EE14-SC

    !! lambda-dependent restraints          
    ti_aug_index(15)=ti_rest_dist_der_dvdl 
    ti_aug_index(16)=ti_rest_ang_der_dvdl 
    ti_aug_index(17)=ti_rest_tor_der_dvdl 
    ti_aug_index(18)=ti_rest_rmsd_der_dvdl     

    do i =1,GPUPotEnergyTerms
     pEnergy(i)%ptr=>null()
    enddo

  endif

end subroutine

subroutine gti_update_ti_pot_energy_from_gpu(pot_ene, enmr, ti_mode, gti_cpu_output)
    
implicit none
    
    type(pme_pot_ene_rec), intent (inout), target :: pot_ene
    real(8), dimension(3), intent (inout), target :: enmr
    integer, intent (in) :: ti_mode, gti_cpu_output

    integer :: i,j
    double precision :: t0, tp(2)

    if (firsttime)  call gti_init

    pEnergy(2)%ptr=>pot_ene%vdw_dir  
    pEnergy(4)%ptr=>pot_ene%bond
    pEnergy(5)%ptr=>pot_ene%angle
    pEnergy(6)%ptr=>pot_ene%dihedral
    pEnergy(8)%ptr=>pot_ene%vdw_14
    pEnergy(7)%ptr=>pot_ene%elec_14   
    pEnergy(9)%ptr=>pot_ene%restraint

    pEnergy(10)%ptr=>pot_ene%elec_recip
    pEnergy(11)%ptr=>pot_ene%elec_dir
    pEnergy(12)%ptr=>pot_ene%angle_ub   
    pEnergy(13)%ptr=>pot_ene%imp   
    pEnergy(14)%ptr=>pot_ene%cmap

    pEnergy(TIExtraShift+1)%ptr=>pot_ene%elec_dir
    pEnergy(TIExtraShift+2)%ptr=>pot_ene%elec_14
 
    gti_potEnergy=0.0
    ti_ene(1,si_dvdl)=0.0 ! Zero out all linear results so far; everything will be recalculated
    call gti_get_pot_energy(gti_potEnergy(1:GPUPotEnergyTerms,  &
        1:TIEnergyBufferMultiplier), ti_mode)

    do i=3, 6
      gti_potEnergy(TIExtraShift+i, 1:2)=ti_ene_tmp(1:2,i)
    end do

    gti_potEnergy(9,3) = gti_potEnergy(9,3) + gti_potEnergy(15, 3) + &
        gti_potEnergy(16, 3) + gti_potEnergy(17, 3) + gti_potEnergy(18, 3)
   
    enmr(1:3)=enmr(1:3)+gti_potEnergy(15:17,3)

    ! Only updates those terms calculated by GPU-TI
    do i=1, TIExtraShift + TINumberExtraTerms
      if (associated(pEnergy(i)%ptr)) then
        if (si_index(i)>0) then
            t0 = 0.0
            if (ti_mode>0) then
              tp(1:2)=gti_potEnergy(i,1:2)
              call ti_update_ene_all(tp, si_index(i), &
                  t0, lambdaType(i))  
            end if
            pEnergy(i)%ptr = pEnergy(i)%ptr + gti_potEnergy(i,3) + t0
        else
            pEnergy(i)%ptr = pEnergy(i)%ptr + gti_potEnergy(i,3)   + &
                sum(gti_potEnergy(i, 1:2) * ti_item_weights( lambdaType(i), :))
        endif
      endif
    end do    

    !Self terms
    if (ti_mode>0) then

      !! SC-only terms
      do i=1, NumberTIEnergyTerms
        if (sc_index(i)>0) then
          ti_ene(1:2,sc_index(i)) = gti_potEnergy(i,TIEnergySCShift+1:TIEnergySCShift+2)
          if (gti_cpu_output .eq. 0) then
            pEnergy(i)%ptr = pEnergy(i)%ptr    + &
                sum(gti_potEnergy(i, TIEnergySCShift+1:TIEnergySCShift+2))
          end if
          ti_ene(1:2,si_pot_ene) = ti_ene(1:2,si_pot_ene) + &
                gti_potEnergy(i,TIEnergySCShift+1:TIEnergySCShift+2)
        endif
      end do    
       
      !! lambda-dependent terms for SC atoms 
      ti_ene_aug(1:2,ti_der_term)= 0.0;
      ti_ene_aug(1:2,ti_rest_der_term) = 0.0; 
      do i=1, TIExtraShift + TINumberExtraTerms
        if (ti_aug_index(i)>0) then
          if (i>=15 .and. i<=18) then  !! 15-17 are for NMR restraints in PI way, 18: RMSD restraint
            t0=gti_potEnergy(i,TIEnergyDLShift+3) 
            ti_ene_aug(1:2,ti_aug_index(i)) = t0
            ti_ene_aug(1:2,ti_rest_der_term) = ti_ene_aug(1:2,ti_rest_der_term) + t0
          else
            t0 = sum( gti_potEnergy(i,TIEnergyDLShift+1:TIEnergyDLShift+2) &
                     * ti_item_weights(lambdaType(i), 1:2))
            if (gti_cpu_output > 0) then
              ti_ene_aug(1:2,ti_aug_index(i))=ti_ene_aug(1:2,ti_aug_index(i)) + t0
            else 
              ti_ene_aug(1:2,ti_aug_index(i)) = &
                 ti_ene_aug(1:2,ti_aug_index(i)) +  &
                 gti_potEnergy(i,TIEnergyDLShift+1:TIEnergyDLShift+2)
            endif
            ti_ene_aug(1:2,ti_der_term) = ti_ene_aug(1:2,ti_der_term) + t0
          endif
        endif
      end do       

      !! Other self terms
      do i=TIExtraShift+3, TIExtraShift+6
          ti_ene(1,si_dvdl) = ti_ene(1,si_dvdl) + &
          sum(gti_potEnergy(i,1:2) &
             *ti_item_dweights(lambdaType(i), 1:2))
      end do
        
      !! final sums
      ti_ene(1:2,si_dvdl) = ti_ene(1:2,si_dvdl) + ti_ene_aug(1:2,ti_der_term) &
            + ti_ene_aug(1:2,ti_rest_der_term)

    endif
        
    pot_ene%vdw_tot = pot_ene%vdw_dir + pot_ene%vdw_recip
    pot_ene%elec_tot = pot_ene%elec_dir +  &
        pot_ene%elec_recip + &
        pot_ene%elec_self + &
        pot_ene%elec_nb_adjust
                   
    pot_ene%total = pot_ene%vdw_tot + &
                    pot_ene%elec_tot + &
        
                    pot_ene%hbond + &
                    pot_ene%bond + &
                    pot_ene%angle + &
                    pot_ene%dihedral + &
                    pot_ene%vdw_14 + &
                    pot_ene%elec_14 + &
                    pot_ene%restraint + &
                    pot_ene%imp + &
                    pot_ene%angle_ub + &
                    pot_ene%efield + &
                    pot_ene%cmap     
             
end subroutine

subroutine gti_update_ti_pot_energy_from_gpu_gamd(pot_ene, enmr, ti_mode, gti_cpu_output)
    
implicit none
    
    type(pme_pot_ene_rec), intent (inout), target :: pot_ene
    real(8), dimension(3), intent (inout), target :: enmr
    integer, intent (in) :: ti_mode, gti_cpu_output

    integer :: i,j
    double precision :: t0 

    if (firsttime)  call gti_init

    pEnergy(2)%ptr=>pot_ene%vdw_dir  
    pEnergy(4)%ptr=>pot_ene%bond
    pEnergy(5)%ptr=>pot_ene%angle
    pEnergy(6)%ptr=>pot_ene%dihedral
    pEnergy(8)%ptr=>pot_ene%vdw_14
    pEnergy(7)%ptr=>pot_ene%elec_14   
    pEnergy(9)%ptr=>pot_ene%restraint

    pEnergy(10)%ptr=>pot_ene%elec_recip
    pEnergy(11)%ptr=>pot_ene%elec_dir
    pEnergy(12)%ptr=>pot_ene%angle_ub   
    pEnergy(13)%ptr=>pot_ene%imp   
    pEnergy(14)%ptr=>pot_ene%cmap

    pEnergy(TIExtraShift+1)%ptr=>pot_ene%elec_dir
    pEnergy(TIExtraShift+2)%ptr=>pot_ene%elec_14
 
    gti_potEnergy=0.0
    ti_ene(1,si_dvdl)=0.0 ! Zero out all linear results so far; everything will be recalculated
    call gti_get_pot_energy(gti_potEnergy(1:GPUPotEnergyTerms,  &
        1:TIEnergyBufferMultiplier), ti_mode)

    do i=3, 6
      gti_potEnergy(TIExtraShift+i, 1:2)=ti_ene_tmp(1:2,i)
    end do

    gti_potEnergy(9,3) = gti_potEnergy(9,3) + gti_potEnergy(15, 3) + &
        gti_potEnergy(16, 3) + gti_potEnergy(17, 3)
   
    enmr(1:3)=enmr(1:3)+gti_potEnergy(15:17,3)

    ! Only updates those terms calculated by GPU-TI
    ti_others_ene(:,:) = 0.d0
    do i=1, NumberTIEnergyTerms
      ! collect ti_others interaction energies between SC and other atoms (TI+I)
      ! si_index(2)=si_vdw_ene
      ! si_index(11)=si_elect_ene
      ! si_index(7)=si_elect_14_ene
      ! si_index(8)=si_vdw_14_ene
      if (si_index(i)>0) then
          if(i.eq.2 .or. i.eq.11 .or. i.eq.7 .or. i.eq.8) then 
             ti_others_ene(1:2,si_index(i)) = ti_others_ene(1:2,si_index(i)) + gti_potEnergy(i,1:2)
             ti_others_ene(1:2,si_pot_ene) = ti_others_ene(1:2,si_pot_ene) + gti_potEnergy(i,1:2)
          endif   
      endif      
    end do

    do i=1, TIExtraShift + TINumberExtraTerms
      if (associated(pEnergy(i)%ptr)) then
        if (si_index(i)>0) then
            t0 = 0.0
            if (ti_mode>0) then
              call ti_update_ene_all(gti_potEnergy(i,1:2), si_index(i), &
                  t0, lambdaType(i))  
            end if
            ! if ((i.le.NumberTIEnergyTerms) then
            if(i.eq.2 .or. i.eq.11 .or. i.eq.7 .or. i.eq.8) then
              ! corrections with using ti_others_ene
              pEnergy(i)%ptr = pEnergy(i)%ptr + gti_potEnergy(i,3)
            else
              pEnergy(i)%ptr = pEnergy(i)%ptr + gti_potEnergy(i,3) + t0
            endif
        else
            pEnergy(i)%ptr = pEnergy(i)%ptr + gti_potEnergy(i,3)   + &
                sum(gti_potEnergy(i, 1:2) * ti_item_weights( lambdaType(i), :))
        endif
      endif
    end do    

    !Self terms
    if (ti_mode>0) then

      !! SC-only terms
      do i=1, NumberTIEnergyTerms
        if (sc_index(i)>0) then
          ti_ene(1:2,sc_index(i)) = gti_potEnergy(i,TIEnergySCShift+1:TIEnergySCShift+2)
          if (gti_cpu_output .eq. 0) then
            pEnergy(i)%ptr = pEnergy(i)%ptr    + &
                sum(gti_potEnergy(i, TIEnergySCShift+1:TIEnergySCShift+2))
          end if
          ti_ene(1:2,si_pot_ene) = ti_ene(1:2,si_pot_ene) + &
                gti_potEnergy(i,TIEnergySCShift+1:TIEnergySCShift+2)
        endif
      end do    

      !! lambda-dependent terms for SC atoms 
      ti_ene_aug(1:2,ti_der_term)= 0.0;
      ti_ene_aug(1:2,ti_rest_der_term) = 0.0; 
      do i=1, TIExtraShift + TINumberExtraTerms
        if (ti_aug_index(i)>0) then
          if (i>=15 .and. i<=17) then  !! 15-17 are for NMR restraints in PI way
            t0=gti_potEnergy(i,TIEnergyDLShift+3) 
            ti_ene_aug(1:2,ti_aug_index(i)) = t0
            ti_ene_aug(1:2,ti_rest_der_term) = ti_ene_aug(1:2,ti_rest_der_term) + t0
          else
            t0 = sum( gti_potEnergy(i,TIEnergyDLShift+1:TIEnergyDLShift+2) &
                     * ti_item_weights(lambdaType(i), 1:2))
            if (gti_cpu_output > 0) then
              ti_ene_aug(1:2,ti_aug_index(i))=ti_ene_aug(1:2,ti_aug_index(i)) + t0
            else 
              ti_ene_aug(1:2,ti_aug_index(i)) = &
                 ti_ene_aug(1:2,ti_aug_index(i)) +  &
                 gti_potEnergy(i,TIEnergyDLShift+1:TIEnergyDLShift+2)
            endif
            ti_ene_aug(1:2,ti_der_term) = ti_ene_aug(1:2,ti_der_term) + t0
          endif
        endif
      end do       

      !! Other self terms
      do i=TIExtraShift+3, TIExtraShift+6
          ti_ene(1,si_dvdl) = ti_ene(1,si_dvdl) + &
          sum(gti_potEnergy(i,1:2) &
             *ti_item_dweights(lambdaType(i), 1:2))
      end do
        
      !! final sums
      ti_ene(1:2,si_dvdl) = ti_ene(1:2,si_dvdl) + ti_ene_aug(1:2,ti_der_term) &
            + ti_ene_aug(1:2,ti_rest_der_term)

    endif
        
    pot_ene%vdw_tot = pot_ene%vdw_dir + pot_ene%vdw_recip
    pot_ene%elec_tot = pot_ene%elec_dir +  &
        pot_ene%elec_recip + &
        pot_ene%elec_self + &
        pot_ene%elec_nb_adjust
                   
    pot_ene%total = pot_ene%vdw_tot + &
                    pot_ene%elec_tot + &
        
                    pot_ene%hbond + &
                    pot_ene%bond + &
                    pot_ene%angle + &
                    pot_ene%dihedral + &
                    pot_ene%vdw_14 + &
                    pot_ene%elec_14 + &
                    pot_ene%restraint + &
                    pot_ene%imp + &
                    pot_ene%angle_ub + &
                    pot_ene%efield + &
                    pot_ene%cmap     

end subroutine

subroutine gti_update_mbar_from_gpu

    implicit none
 
    double precision, dimension(:,:,:), allocatable, save::  weight, dWeight
    double precision, dimension(:), allocatable, save::  lambda1, lambda2
    double precision, dimension(TypeTotal, 2), save:: currentWeight, currentdWeight
    double precision, save::clambda_old=200
    integer, save:: l_current=1
    logical, save:: init=.false., keep=.false.
    
    double precision, dimension(:), allocatable::  MBAREnergy
    double precision, dimension(TypeTotal, 3) :: currentPE 

    integer :: i,j
    double precision:: t0

    if (.not. init) then
      if (bar_states>0 .and. do_mbar)  then
          allocate(&
            lambda1(bar_states), lambda2(bar_states), &
            weight(bar_states, 2, TypeTotal), &
            dWeight(bar_states, 2, TypeTotal) )
      else
          return
      endif
      init=.true.
    endif 
    
    allocate(MBAREnergy(bar_states*TypeTotal *3))
    MBAREnergy=0.0

    keep = (dabs(clambda-clambda_old)<1e-10)

    currentPE(:,:)=0.0
    do i=1, TIExtraShift + TINumberExtraTerms
        currentPE(lambdaType(i),1:3)= currentPE(lambdaType(i),1:3) &
          +gti_potEnergy(i,1:3)
    end do
  
    !! MBAREnergy stores only the lambda-dependent contributions
    call gti_get_mbar_energy(MBAREnergy(1:bar_states*TypeTotal*3))
    
    if (.not.keep) then
      lambda1(1)=clambda
      lambda2(1)=1.0-clambda
      do i=0, TypeTotal-1 !! running over all types
        call gti_get_lambda_weights(i, 1, lambda1, &
             currentWeight(i+1, 1), currentdWeight(i+1, 1))
        call gti_get_lambda_weights(i, 1, lambda2, & 
            currentWeight(i+1, 2), currentdWeight(i+1, 2))
      end do

      t0=999.0
      do i=1, bar_states
         bar_cont(i)=0.0
         lambda1(i)=1.0-bar_lambda(1, i, 1)
         lambda2(i)=1.0-bar_lambda(2, i, 1)
         if (t0 .gt. abs(lambda1(i)-clambda)) then
            t0=abs(lambda1(i)-clambda)
            l_current=i
         endif
      enddo

      do i=0, TypeTotal-1 !! running over all types; C-index
        call gti_get_lambda_weights(i, bar_states, lambda1, weight(:,1,i+1), dWeight(:,1,i+1))
        call gti_get_lambda_weights(i, bar_states, lambda2, weight(:,2,i+1), dWeight(:,2,i+1))
      end do
    end if

    do i=1, bar_states
      bar_cont(i)=0.0
    enddo

    do i=1, TypeTotal !! running over all types; FORTRAN-index
      do j=1, bar_states
        bar_cont(j) = bar_cont(j) + & 
            (weight(j, 1, i) - currentWeight(i, 1) )* currentPE(i, 1)
        bar_cont(j) = bar_cont(j) + &
            MBAREnergy(j + bar_states * 3*(i-1) )* weight(j, 1, i)
      end do

      do j=1, bar_states
        bar_cont(j) = bar_cont(j) + & 
            (weight(j, 2, i) - currentWeight(i, 2) )* currentPE(i, 2)
        bar_cont(j) = bar_cont(j) + &
            MBAREnergy(j + bar_states * (1+3*(i-1)  ) )* weight(j, 2, i)
      end do

      do j=1, bar_states
        bar_cont(j) = bar_cont(j) + &
            MBAREnergy(j + bar_states * (2+3*(i-1)  ) )
      end do
    end do

    t0=minval(bar_cont(1:bar_states))
    t0=bar_cont(l_current)
    bar_cont(1:bar_states)=bar_cont(1:bar_states)-t0

    clambda_old=clambda

    if (allocated(MBAREnergy)) deallocate(MBAREnergy)
   
end subroutine


subroutine gti_update_kin_energy_from_gpu(gpuIndex, outIndex, update)
    
    implicit none
    integer, intent(in):: gpuIndex, outIndex
    logical, intent(in), optional:: update
    
!local variables
 
    integer, parameter:: numberKinEnergyTerms = 10
    real(8), dimension(numberKinEnergyTerms,3),save :: gti_kinEnergy
    
    if (.not. present(update) ) then
        call gti_get_kin_energy(gti_kinEnergy(1:numberKinEnergyTerms, 1:3))    
    else if (update) then
        call gti_get_kin_energy(gti_kinEnergy(1:numberKinEnergyTerms, 1:3))    
    end if   

    ti_kin_ene(1:2, outIndex) = gti_kinEnergy(outIndex,1:2)
    if (outIndex==gpuIndex) then
        ti_kin_ene(1:2, ti_sc_eke) = gti_kinEnergy(ti_sc_eke,1:2)
    else
        ti_kin_ene(1:2, ti_sc_eke) = gti_kinEnergy(gpuIndex,1:2)
    endif
              
end subroutine

subroutine gti_update_force_from_gpu(atm_cnt, frc, needVirial, index )

    implicit none

    integer, intent(in) ::atm_cnt   
    logical, intent(in) ::needVirial
    integer, optional, intent(in) :: index
    double precision, intent (out)  :: frc(3, atm_cnt)

!local

    double precision:: frc0(3, atm_cnt), frc1(3, atm_cnt), & 
        frc2(3, atm_cnt), frcX(3, atm_cnt)      
  
    call gti_download_frc(0, frc0, needVirial, .true.)   
    call gti_download_frc(1, frc1, needVirial, .false.)     
    call gti_download_frc(2, frc2, needVirial, .false.)  
    frcX=frc0*ti_weights(1)+frc1*ti_weights(2)+frc2;    
    
    if (present(index)) then
        if (index==0) frc=frc0
        if (index==1) frc=frc1
        if (index==2) frc=frc2  
        if (index==3) frc=frcX    
    else
        call gti_download_frc(-1, frc, needVirial, .true.)   
    endif

end subroutine

subroutine gti_add_ti_to_exclusion(atm_cnt, numex, n, natex)

  implicit none
  integer, intent(in)               :: atm_cnt, n
  integer, intent(inout)               :: numex(atm_cnt)
  integer, dimension(:), intent(inout)               :: natex
 
  integer:: oldIndex, newIndex
  integer:: i1,i2
  integer:: i,j,k
  integer, dimension(:), allocatable :: tempEx, region1List
 
  allocate(tempEx(n))
  tempEx=natex

  i1=ti_ti_atm_cnt(1)
  i2=ti_ti_atm_cnt(2)
  
  allocate(region1List(i1))
  k=1
  do i=1, atm_cnt
      if (ti_atm_lst(1,i)>0) then
          region1List(k)=i
          k=k+1
      endif
  enddo
   
  oldIndex=1
  newIndex=1
  do i=1, atm_cnt
     
    do j=0, numex(i)-1
        natex(newIndex+j)=tempEx(oldIndex+j)
    end do
          
    newIndex=newIndex+numex(i)
    oldIndex=oldIndex+numex(i)
      
    if (ti_lst(2, i)>0) then
        numex(i)=numex(i)+i1
        do j=0, i1-1
            natex(newIndex+j)=region1List(j+1)
        enddo
        newIndex=newIndex+i1
    endif
      
  end do
  
  deallocate(tempEx)
  deallocate(region1List)

end subroutine

subroutine gti_output_ti_result(time)

  implicit none

  double precision,intent(in)::time
  
  integer:: i,j,k
  integer,save::counter=0, counter_cut=200
  double precision, save::average(2)=0.0, average_sq(2)=0.0
  double precision, dimension(2,2)::tt
  double precision:: t1, t2, t3, t11, t22, t33

  if (ti_mode.eq.0) return

!! Detailed TI output
  if (gti_output.ne.0) then
    t11=0.0
    t22=0.0
    t33=0.0
    write(mdout, *) "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    write(mdout, *) "Detailed TI info at lambda=", clambda

    write(mdout,'(A6, A7, A14, 3A12)') "Region", "  ", &
      "H", "W",  &  
      "dH/dl", "dW/dl"
    do i=1, TIExtraShift + TINumberExtraTerms
      if (si_index(i) .gt. 0) then

        ! tt(:,1) holds H and tt(:,2) holds dH/dl
        tt(1:2,1)=gti_potEnergy(i,1:2)  !! linear terms
        tt(1:2,2)=gti_potEnergy(i,TIEnergyDLShift+1:TIEnergyDLShift+2) !! lambda-dep terms
        if (norm2(tt).ne.0) then
          do j=1,2
              write(mdout,'(A2, I2, A1, A11, F14.5, 3F12.5)') "TI", j, " ", items(i), &
                tt(j, 1), ti_item_weights(lambdaType(i),j),   &
                tt(j, 2), ti_item_dweights(lambdaType(i),j)
          end do
        end if

        t3=gti_potEnergy(i,TIEnergyDLShift+3)
        if (t3.ne.0) then
            write(mdout,'(A3, A11, A13, F14.5)') "TI ", items(i), "PI-TI cont: ", t3
        endif
        t1=sum(tt(1:2,1)*ti_item_dweights(lambdaType(i), 1:2)) ! dW(0)*H(0) + dW(1)*H(1)
        t2=sum(tt(1:2,2)*ti_item_weights(lambdaType(i), 1:2))  ! W(0)*dH(0) + W(1)*dH(1)

        if (norm2(tt).ne.0) then
          write(mdout,'("lambda = ", F5.3," : ", A9, "H=", F13.4, &
            " dU/dL: L=", F10.4 , " NL=", F10.4 ," Tot=", F11.5)' )  &
              clambda, items(i), &
              sum(tt(1:2,1)*ti_item_weights(lambdaType(i), 1:2)), & ! W(0)*H(0) + W(1)*H(1)
              t1, t2, t1+t2
          t11=t11+t1
          t22=t22+t2
          write(mdout, *) "------------------------------------------------------------------------"
        endif
        if (t3.ne.0) then
          write(mdout,'("lambda = ", F5.3," : ", A9, " PI=", F10.4)' )  &
                  clambda, items(i), t3
          t33=t33+t3
          write(mdout, *) "------------------------------------------------------------------------"
        endif

      endif

    end do      
    write(mdout,'("lambda = ", F5.3," : Total dU/dl:", F12.6, &
         "  L:", F11.5, "  NL:", F11.5, "  PI:", F11.5)') &
         clambda, t11+t22+t33, t11, t22, t33

    counter=counter+1
    if (counter .gt. counter_cut .and. norm2(gti_potEnergy(TISpecialShift+1,1:2)).ne.0) then
      t11=KB*temp0
      t22=1.0/t11
      t33=((counter-counter_cut)*1.0)
      specicalEstimates(1:2)= specicalEstimates(1:2) + exp( t22* &
          (gti_potEnergy(TISpecialShift+1,1:2) - gti_potEnergy(TISpecialShift+2,1:2)))
      write(mdout,'(A11, " #1: ", F10.5, " #2: ", F10.5, &
       "      ", A11, " #1: ", F10.5, " #2: ", F10.5)') &
      items(TISpecialShift+1), gti_potEnergy(TISpecialShift+1,1:2), &
      items(TISpecialShift+2), gti_potEnergy(TISpecialShift+2,1:2)
      t1=-log(specicalEstimates(1)/t33)*t11
      t2=-log(specicalEstimates(2)/t33)*t11
      
      average(1)= average(1)+t1
      average_sq(1) = average_sq(1) + t1*t1
      average(2)= average(2)+t2
      average_sq(2) = average_sq(2) + t2*t2
      write(mdout,'("Curr to Orig correction : #1: ", F10.4, " +-", F7.4, &
                                             " #2: ", F10.4, " +-", F7.4 )') & 
         t1, sqrt(average_sq(1)/t33 - (average(1)/t33)*(average(1)/t33)  ),  &
         t2, sqrt(average_sq(2)/t33 - (average(2)/t33)*(average(2)/t33)  )
    end if

    write(mdout, *) "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

  endif

!! SAMS output
  if (ifsams >0 .and. sams_onstep .and. sams_type .le. 2 ) then
    write(mdout, *) "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    write(mdout, *) "SAMS info"
    write(mdout, 9088) currentLambdaIndex, clambda, ti_ene(1,si_dvdl)
    write(mdout, *) "Iteration #", currentIteration, " Stage ", burnin_stage
    if (burnin_stage .eq. 0) then
      write(mdout, *) " State histogram:"
      write(mdout, '(20i4)') lambda_counter(1:bar_states)
      write(mdout, *) "logZ:"
      write(mdout, '(10f10.4)') logZ_report(1:bar_states)
    else
      write(mdout, &
        '(" Time: ",F10.2, " deltaF: ",F10.6, " dvdl_ave:", F10.6)') &
                time, -logZ_report(bar_states)*temp0*KB, dvdl_ave
    end if
    write(mdout, *) "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
  endif
   
9088 format(1x, 'Current Index:', i4, ' lambda = ', f14.6, 2x, 'DV/DL  = ',f14.6) 

end subroutine

#endif  /* GTI */

end module
