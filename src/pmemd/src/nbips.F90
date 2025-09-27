#include "copyright.i"
! <compile=optimized>

module nbips_mod
  
  !
  !      variables for Isotropic Periodic Sum (IPS) calculation
  !
  !      IPS         IPS options: 1--for both ele and l-j using 3D IPS
  !                               2--for ele only using 3D IPS
  !                               3--for l-j only using 3D IPS
  !                               4--for both ele and l-j using 3D IPS/DFFT
  !                               5--for ele only using 3D IPS/DFFT
  !                               6--for l-j only using 3D IPS/DFFT
  !      NNBIPS      Number of nonbonded atom pairs
  !      NNBIPST     Provious Number of nonbonded atom pairs
  !
  !integer,save:: NNBIPST,NNBIPS,IPSSIZ
  integer,save:: NNBIPS,IPSSIZ
  double precision,save:: NNBIPST
  
  
  !    IPS parameters
  !    AIPSE(0:3)    Electrostatic IPS coefficients
  !    BIPSE(3)    Electrostatic IPS derivative coefficients
  !    AIPSVA(0:3)   Lennard-Jones repulsion IPS coefficients
  !    BIPSVA(3)   Lennard-Jones repulsion IPS derivative coefficients
  !    AIPSVC(0:3)   Lennard-Jones dispersion IPS coefficients
  !    BIPSVC(3)   Lennard-Jones dispersion IPS derivative coefficients
  !    RIPS      Radius of IPS local region 
  !    RAIPS     Radius of IPS extensive local region 
  !    PIPS*0    Self IPS pair energies 
  !    PIPS*C    IPS boundary energies 
  !    EIPSS*C   IPS boundary energies 
  !    VIRIPS*C  IPS boundary virials 
  !    VBOXIPS   IPS local region volume
  !
  
  double precision,save::  AIPSE(0:3),BIPSE(3), &
       AIPSVA(0:3),BIPSVA(3), &
       AIPSVC(0:3),BIPSVC(3), &
       RIPS,RIPS2,RIPSR,RIPS2R,RIPS6R,RIPS12R, &
       RAIPS,GRIDIPS,VBOXIPS, &
       PIPSE0,PIPSVA0,PIPSVC0,PIPSEC,PIPSVAC,PIPSVCC,  &
       EIPSSNB,EIPSSEL,VIRIPS,EIPSANB,EIPSAEL,VIRAIPS
  
  double precision,save::  CGSUM,CIJSUM,AIJSUM,VIREXIPS(3,3)
  
  
  INTEGER,save::  MIPSX,MIPSY,MIPSZ,MIPSO
  ! To enhance performance in pmemd we have decided to allow only 
  ! IPS calculations with teips and tvips set to TRUE so we won't check for them
  !     TEIPS   ! Perform IPS for electrostatic interaction
  !     TVIPS   ! Perform IPS for Lennard-Jones interaction
  !
  
  integer,save :: ierr_allocate,alloc_err
  
  
  double precision,allocatable,dimension(:),save :: elearray,vdwarray
  double precision,allocatable,dimension(:),save :: elexx,elexy, &
       elexz,eleyy,eleyz,elezz,vdwxx,vdwxy,vdwxz,vdwyy,vdwyz,vdwzz
  
  double precision,allocatable,dimension(:),save ::prefac1,prefac2,prefac3
  !
  integer,save :: forder
  !
  INTEGER,save :: NFFT1,NFFT2,NFFT3,nfftdim1,nfftdim2,nfftdim3, &
       nfftable,nffwork,SIZ_Q,SIZFFTAB,SIZFFWRK
  double precision,allocatable,dimension(:),save :: fftable,ffwork
  
  !===========================================================================
contains
  !===========================================================================
  
  !call ipssys(natom,ntypes,ntb,x(l15), &
  !            cut,cn1,cn2,ix(i04),ix(i06),x(lcrd))
  SUBROUTINE IPSSYS(X)
    !-----------------------------------------------------------------------
    !    Calculate system energy and force in 3D periodic systems
    !
    use img_mod
    use timers_mod
    use mdin_ctrl_dat_mod
    use mdin_ewald_dat_mod
    use parallel_mod
    use parallel_dat_mod
    use prmtop_dat_mod!contains natom
    use ene_frc_splines_mod
    use gbl_datatypes_mod
    use gbl_constants_mod
    
    implicit none
    
    !#include "extra.h"
    
    !Formal arguments
    double precision, intent(in)     ::  X(3,*)
    !    double precision atm_qterm(*),gbl_cn1(*), gbl_cn2(*)
    !    INTEGER gbl_img_iac(*),typ_ico(*)
    
    !Internal variables
    double precision  CUT,sqrtcut
    INTEGER I,J
    INTEGER ITI,ITJ,ITMAX,IACI,IC!,NITI,NITJ
    double precision  NITI,NITJ
    double precision  CISUMI
    double precision  XI,YI,ZI,XIJ,YIJ,ZIJ,R2
    double precision  AIJ,CIJ,ANIJ
    double precision  CGIJSUM
    
    !INTEGER,allocatable,dimension(:) :: NSUMIT
    double precision,allocatable,dimension(:) :: NSUMIT
    double precision,allocatable,dimension(:) :: CISUM
#ifdef MPI

  ! To enhance performance in pmemd we have decided to allow only 
  ! IPS calculations with teips and tvips set to TRUE, i.e. ips=1

    teips=.true.
    tvips=.true.
#endif
    
    !  IPS Radius:
    sqrtcut = ips_cut
    cut = sqrtcut * sqrtcut

    RIPS2=cut
    RIPS2R=1.d0/cut
    RIPS6R=RIPS2R*RIPS2R*RIPS2R
    RIPS12R=RIPS6R*RIPS6R
    RIPS=SQRT(RIPS2)
    RIPSR=1.d0/RIPS
    
    !  Ele IPS parameters:
    AIPSE(0)=-35.0/16.0
    AIPSE(1)=35.0/16.0
    AIPSE(2)=-21.0/16.0D0
    AIPSE(3)=5.0D0/16.0D0
    !          PIPSEC=0.D0
    PIPSEC=1.D0+AIPSE(0)+AIPSE(1)+AIPSE(2)+AIPSE(3)
    PIPSE0=AIPSE(0)-PIPSEC
    BIPSE(1)=2.D0*AIPSE(1)
    BIPSE(2)=4.D0*AIPSE(2)
    BIPSE(3)=6.D0*AIPSE(3)
    
    !  Dispersion IPS parameters:
    AIPSVC(0)=7.0D0/16.0D0
    AIPSVC(1)=9.0D0/14.0D0
    AIPSVC(2)=-3.0D0/28.0D0
    AIPSVC(3)=6.0D0/7.0D0
    !          PIPSVCC=67.0D0/28.0D0+7.0D0/16.0D0
    PIPSVCC=1.D0+AIPSVC(0)+AIPSVC(1)+AIPSVC(2)+AIPSVC(3)
    PIPSVC0=AIPSVC(0)-PIPSVCC
    BIPSVC(1)=2.D0*AIPSVC(1)
    BIPSVC(2)=4.D0*AIPSVC(2)
    BIPSVC(3)=6.D0*AIPSVC(3)
    !  Repulsion IPS parameters:
    AIPSVA(0)=5.0D0/787.0D0
    AIPSVA(1)=9.0D0/26.0D0
    AIPSVA(2)=-3.0D0/13.0D0
    AIPSVA(3)=27.0D0/26.0D0
    !        PIPSVAC=28.0D0/13.0D0+5.0D0/787.0D0
    PIPSVAC=1.D0+AIPSVA(0)+AIPSVA(1)+AIPSVA(2)+AIPSVA(3)
    PIPSVA0=AIPSVA(0)-PIPSVAC
    BIPSVA(1)=4.D0*AIPSVA(1)
    BIPSVA(2)=8.D0*AIPSVA(2)
    BIPSVA(3)=2.D0*6.D0*AIPSVA(3)
    
    !
    EIPSSNB=0.0D0
    EIPSSEL=0.0D0
    
    !=======================================================================
    !   Main loop begin
    !=======================================================================

    if(allocated(NSUMIT))deallocate(NSUMIT)
    allocate(NSUMIT(ntypes),stat=alloc_err)
    if(alloc_err /= 0 ) write(0,*)"unable to allocate NSUMIT"
    if(allocated(CISUM))deallocate(CISUM)
    allocate(CISUM(ntypes),stat=alloc_err)
    if(alloc_err /= 0 ) write(0,*)"unable to allocate CISUM"
    
    DO ITI=1,ntypes
       NSUMIT(ITI)=0
       CISUM(ITI)=0
    end do
    CGSUM=0.0D0
    ITMAX=0
    DO I=1,NATOM
       CGSUM=CGSUM + atm_qterm(I)
       ITI=atm_iac(I)
       IF(ITI>ntypes)STOP "problem in IAC!"
       NSUMIT(ITI)=NSUMIT(ITI)+1
    end do
    CIJSUM=0.0D0
    AIJSUM=0.0D0
    NNBIPST=0
    
    !  system energy is calculated based on all pairs:
    DO ITI=1,ntypes
       iaci = ntypes * (ITI - 1)
       ic = typ_ico(iaci+ITI)
       !        IC=ITI*(ITI-1)/2+ITI
       if (ic .gt. 0) then
          AIJ=gbl_cn1(IC)
          CIJ=gbl_cn2(IC)
       else
          AIJ=0.0d0
          CIJ=0.0d0
       end if
       NITI=NSUMIT(ITI)
       ANIJ=NITI*NITI*0.5
       CISUMI=CISUM(ITI)+CIJ*NITI
       CIJSUM=CIJSUM+CIJ*ANIJ
       AIJSUM=AIJSUM+AIJ*ANIJ
       NNBIPST=NNBIPST+NITI*NITI
       DO ITJ=ITI+1,ntypes
          NITJ=NSUMIT(ITJ)
          !   iaci = ntypes * (gbl_img_iac(i) - 1)
          !   ic = typ_ico(iaci+gbl_img_iac(j))
          ic = typ_ico(iaci+ITJ)
          !          IC=ITJ*(ITJ-1)/2+ITI
          if (ic .gt. 0) then
            AIJ=gbl_cn1(IC)
            CIJ=gbl_cn2(IC)
          else
            AIJ=0.0d0
            CIJ=0.0d0
          endif
          ANIJ=NITI*NITJ !r
          CIJSUM=CIJSUM+CIJ*ANIJ !r
          AIJSUM=AIJSUM+AIJ*ANIJ !r
          CISUMI=CISUMI+CIJ*NITJ !r
          CISUM(ITJ)=CISUM(ITJ)+CIJ*NITI !r
          NNBIPST=NNBIPST+2*NITI*NITJ
       end do
       CISUM(ITI)=CISUMI
    end do
    CGIJSUM=CGSUM*CGSUM*0.5
    if( cgijsum < 1.D-15 ) cgijsum = 0.d0
    IF(TVIPS.and.master)THEN
       EIPSSNB=AIJSUM*PIPSVAC*RIPS12R
       EIPSSNB=EIPSSNB-CIJSUM*PIPSVCC*RIPS6R
    ELSE
       EIPSSNB=0.0D0
    end if
    IF(TEIPS.and.master)THEN
       EIPSSEL=0.0D0
       EIPSSEL=CGIJSUM*PIPSEC*RIPSR
    ELSE
       EIPSSEL=0.0D0
    end if
    deallocate(NSUMIT,CISUM)
    
    !=======================================================================
    !   Main loop end
    !=======================================================================
    
    ! Calculate volume virials:
    
    VIRIPS=-(EIPSSNB+EIPSSEL)
    VBOXIPS=4.d0*PI*RIPS2*RIPS*(1.d0/3.d0)
    IF(NTB.EQ.0)THEN
       ! For non-periodic system, system energy is calculated based on cutoff
       NNBIPS=0
       DO I=1,NATOM
          ! Atom i long-range reference and self-interaction:
          XI=X(1,I)
          YI=X(2,I)
          ZI=X(3,I)
          NNBIPS=NNBIPS+1
          DO J=I+1,NATOM
             XIJ=XI-X(1,J)
             YIJ=YI-X(2,J)
             ZIJ=ZI-X(3,J)
             R2=XIJ*XIJ+YIJ*YIJ+ZIJ*ZIJ
             !  RIPS2 is infinite except for non periodic conditions
             IF(R2<RIPS2)NNBIPS=NNBIPS+2
          end do
       end do
    end if
    
    if( master ) then
       WRITE(mdout,'(" ----------------------------------")')
       WRITE(mdout,'(" Using 3D-IPS algorithm")')
       WRITE(mdout,'("   IPS Radius: ",F6.2," A")')RIPS
       WRITE(mdout,'("   Using IPS for electrostatic energy")')
       WRITE(mdout,'("   Using IPS for L-J energy")')
       WRITE(mdout,'(" ----------------------------------")')
    else 
       EIPSSNB=0.0D0
       EIPSSEL=0.0D0
       VIRIPS=0.0D0
       NNBIPS=0
    end if

#ifdef CUDA
    call gpu_ips_setup(RIPS, atm_qterm, ntypes, atm_iac, typ_ico, gbl_cn1, gbl_cn2, EIPSSNB, EIPSSEL, VIRIPS)
#endif       
    
    RETURN 
    
  end subroutine ipssys
  
  
  SUBROUTINE IPSUPDATE(NTB)
    
    !-----------------------------------------------------------------------
    !     Update parameters once IPS radius or the box size changed
    !-----------------------------------------------------------------------
    
    use mdin_ctrl_dat_mod, only : dvbips 
#ifdef MPI
    use parallel_mod
    use parallel_dat_mod
#endif
    use pbc_mod, only : uc_volume
    implicit none
#ifdef MPI
    INTEGER NNBTMP
#endif
    INTEGER NTB
    double precision FIPS,CHANGE 
    
    IF(NTB>0)THEN
       FIPS=VBOXIPS/uc_volume
    ELSE
#ifdef MPI
       call mpi_allreduce(NNBIPS,NNBTMP,1,MPI_INTEGER, &
            mpi_sum,pmemd_comm,err_code_mpi)
       NNBIPS=NNBTMP 
#endif
       FIPS=NNBIPS*1.0D0/NNBIPST
    end if
    CHANGE=ABS(FIPS-1.0D0)
    
    ! Update system energies and forces:
    IF(CHANGE>DVBIPS)THEN
       VIRIPS=VIRIPS*FIPS
       EIPSSNB=EIPSSNB*FIPS
       EIPSSEL=EIPSSEL*FIPS
       VBOXIPS=uc_volume
       NNBIPST=NNBIPS
    end if
#ifdef CUDA
    call gpu_ips_update(EIPSSNB, EIPSSEL, VIRIPS)
#endif    
    
    RETURN
  end subroutine ipsupdate
  
  
  
#ifdef MPI
  SUBROUTINE EEXIPS(ENB,EEL,frc_vec,crd_vec,img_qterm,ipairs,atm_maskdata, atm_mask,img_atm_map,atm_cnt,tranvec,my_atm_lst)
#else
  SUBROUTINE EEXIPS(ENB,EEL,frc_vec,crd_vec,img_qterm,ipairs,atm_maskdata, atm_mask,img_atm_map,atm_cnt,tranvec)
#endif    
    !-----------------------------------------------------------------------
    !   3D IPS interaction between excluded atom pairs
    !   This routine must be called first to update IPS parameters when needed
    !       and to initalize electrostatic and vdw energies
    !
    !   by Romelia 
    !
    !-----------------------------------------------------------------------
    
    use img_mod
    use timers_mod
    use mdin_ctrl_dat_mod
    use mdin_ewald_dat_mod
    use parallel_mod
    use parallel_dat_mod
    use prmtop_dat_mod
    use ene_frc_splines_mod
    use nb_pairlist_mod

    implicit none
    ! Formal arguments:
    
    double precision, intent(in out) :: frc_vec(3, *)
    double precision, intent(in)     :: crd_vec(3, *)
    double precision , intent(in)    :: img_qterm(*)
    integer, intent(in)              :: ipairs(*)
    double precision, intent(in out) :: ENB, EEL
    type(listdata_rec)    :: atm_maskdata(*)
    integer               :: atm_mask(*)
    integer               :: img_atm_map(*)
    integer               :: atm_cnt
    double precision      :: tranvec(1:3, 0:17)
#ifdef MPI
  integer                       :: my_atm_lst(*)
#endif

    ! local:
    INTEGER I,J,K,L,IC,IACI,ITI,ITJ,Iim,Jim
    INTEGER  ILAST, IFIRST
    INTEGER ee_eval_cnt, full_eval_cnt,ipairs_idx
    double precision DXI, DYI, DZI,DIJ,DIJX,DIJY,DIJZ
    double precision ENBIJ,EELIJ
    double precision  XI,YI,ZI,XIJ,YIJ,ZIJ
    double precision  CGI,CGIJ,AIJ,CIJ
    double precision  R2,R6,R12,U1,U2,U4
    double precision  PE,DEU,PVA,DVAU,PVC,DVCU
    INTEGER common_tran
    INTEGER atm_mask_idx, atm_i
    double precision      x_tran(1:3, 0:17)
    INTEGER itran
    integer atm_lst_idx

    ! check to see if volume or atom pair changed 
    CALL IPSUPDATE(NTB)
    
    I=0
    J=0
    k=0
    Iim=0
    Jim=0

    ! Setup constants for use in inner loops
    ENB=EIPSSNB
    EEL=EIPSSEL
    DEU=0.d0
    DVAU=0.d0
    DVCU=0.d0
    CGIJ=0.d0
    AIJ=0.d0
    CIJ=0.d0
    VIREXIPS = 0.d0
    VIREXIPS(1,1) = virips
    VIREXIPS(2,2) = virips
    VIREXIPS(3,3) = virips
    
    !=======================================================================
    !   Main loop begin
    !=======================================================================
    
    NNBIPS=0
    ipairs_idx = 1
    
#ifdef MPI
      do atm_lst_idx = 1, atm_cnt
        i = my_atm_lst(atm_lst_idx)
         Iim=i
#else
      do i = 1, atm_cnt
         Iim=i
#endif /*  MPI */

      ! Mark excluded j images:

      atm_mask_idx = atm_maskdata(I)%offset
       NNBIPS=NNBIPS+1
   
       CGI=atm_qterm(I)
       CGIJ=CGI*CGI
       EELIJ=0.5*CGIJ*PIPSE0*RIPSR
       EEL=EEL+EELIJ
       ITI=atm_iac(I)
       iaci = ntypes * (iti - 1)
       ic = typ_ico(iaci+ITI)
       if (ic .gt. 0) then
          AIJ=gbl_cn1(IC)
          CIJ=gbl_cn2(IC)
       else
          AIJ=0.d0
          CIJ=0.d0
       end if
       ! Atom i long-range reference and self-interaction
       ENBIJ=0.5*(AIJ*PIPSVA0*RIPS6R-CIJ*PIPSVC0)*RIPS6R
       ENB=ENB+ENBIJ

       DXI=frc_vec(1,Iim)
       DYI=frc_vec(2,Iim)
       DZI=frc_vec(3,Iim)
       XI=crd_vec(1,I)
       YI=crd_vec(2,I)
       ZI=crd_vec(3,I)
       
       ! Electrostatic evaluation-only count followed by
       ! full evaluation count packed at the front of each pair sublist.
       
         do K = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(I)%cnt
          J = atm_mask(K)
#ifdef MPI
         Jim=J
#else
         Jim=J
#endif /*  MPI */
          IF(J<=I) cycle

          XIJ=crd_vec(1,J)-XI
          YIJ=crd_vec(2,J)-YI
          ZIJ=crd_vec(3,J)-ZI
          R2=XIJ*XIJ+YIJ*YIJ+ZIJ*ZIJ
          NNBIPS=NNBIPS+2
          R6=R2*R2*R2
          R12=R6*R6
          U2=R2*RIPS2R
          U1=SQRT(U2)

          !Electrostatics  
          !This routine adds the contribution to the electrostatic 
          !nonbonding energy for the images of atoms participating 
          !in bonds, angles etc it does not contain the 1/r term
          !(see pairs_ips_(no)vec.i 

          CGIJ=CGI*atm_qterm(J)/RIPS
          PE=AIPSE(0)+U2*(AIPSE(1)+U2*(AIPSE(2)+U2*(AIPSE(3))))
          DEU=U2*(BIPSE(1)+U2*(BIPSE(2)+U2*(BIPSE(3))))
          EELIJ=CGIJ*(PE-PIPSEC)
          EEL=EEL+EELIJ
          DIJ=-CGIJ*DEU/R2
          ITJ=atm_iac(J)
          ic = typ_ico(iaci+ITJ)
          if (ic .gt. 0) then
             AIJ=gbl_cn1(IC)*RIPS12R
             CIJ=gbl_cn2(IC)*RIPS6R
          else
             AIJ=0.d0
             CIJ=0.d0
          end if
          U4=U2*U2

          ! Lennard Jones
          !This routine adds the contribution to the L-J nonbonding energy 
          !for the images of atoms participating in bonds, angles etc 
          !it does not contain the 1/r6 for L-J r6 and 1/r12 for L-J r12 terms
          !(see pairs_ips_(no)vec.i 

          !  L-J r6 term
          !   etr6=a0+r2*(a1+r2*(a2+a3*r2))
          !   detr6/dr*r1=r2*(d1+r2*(d2+d3*r2))
          !
          PVC=CIJ*(AIPSVC(0)+U2*(AIPSVC(1)+U2*(AIPSVC(2)+U2*(AIPSVC(3))))-PIPSVCC)
          DVCU=CIJ*(U2*(BIPSVC(1)+U2*(BIPSVC(2)+U2*BIPSVC(3))))
          !  L-J r12 term 
          !   etr12=a0+r2*(a1+r2*(a2+a3*r2))
          !   detr12/dr*r1=r2*(d1+r2*(d2+d3*r2))
          !
          PVA=AIJ*(AIPSVA(0)+U4*(AIPSVA(1)+U4*(AIPSVA(2) +U4*(AIPSVA(3))))-PIPSVAC)
          DVAU=AIJ*(U4*(BIPSVA(1)+U4*(BIPSVA(2)+U4*BIPSVA(3))))
          ENBIJ=PVA-PVC
          ENB=ENB+ENBIJ
          DIJ=DIJ-(DVAU-DVCU)/R2
          
          DIJX=DIJ*XIJ
          DIJY=DIJ*YIJ
          DIJZ=DIJ*ZIJ
          DXI=DXI-DIJX
          DYI=DYI-DIJY
          DZI=DZI-DIJZ
          frc_vec(1,Jim)=frc_vec(1,Jim)+DIJX
          frc_vec(2,Jim)=frc_vec(2,Jim)+DIJY
          frc_vec(3,Jim)=frc_vec(3,Jim)+DIJZ
          VIREXIPS(1,1) = VIREXIPS(1,1) - DIJX*XIJ
          VIREXIPS(1,2) = VIREXIPS(1,2) - DIJX*YIJ
          VIREXIPS(1,3) = VIREXIPS(1,3) - DIJX*ZIJ
          VIREXIPS(2,2) = VIREXIPS(2,2) - DIJY*YIJ
          VIREXIPS(2,3) = VIREXIPS(2,3) - DIJY*ZIJ
          VIREXIPS(3,3) = VIREXIPS(3,3) - DIJZ*ZIJ
       end do
       frc_vec(1,Iim) =  DXI
       frc_vec(2,Iim) =  DYI
       frc_vec(3,Iim) =  DZI
    end do
    
    !=======================================================================
    !   Main loop end
    !=======================================================================
    VIREXIPS(2,1)=VIREXIPS(1,2)
    VIREXIPS(3,1)=VIREXIPS(1,3)
    VIREXIPS(3,2)=VIREXIPS(2,3)


   RETURN
    
  end subroutine EEXIPS
  
end module nbips_mod
