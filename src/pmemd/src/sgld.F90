#include "copyright.i"
!!TODO connect sgld-remd, check sandersize==numtasks? and print sgld energy to mdout
!not doing allreduce in place

!*******************************************************************************
! module: sgld_mod
! ================== self-guided molecular/langevin dynamics ====================
! romelia salomon, 2013
!*******************************************************************************

module sgld_mod

  use file_io_dat_mod
  use parallel_dat_mod
  use mdin_ctrl_dat_mod
  use gbl_datatypes_mod


  implicit none

  type sgld_ene_rec
    sequence
    double precision    :: templf
    double precision    :: temphf
    double precision    :: epotlf
    double precision    :: epothf
    double precision    :: epotllf
    double precision    :: sggamma
    double precision    :: sgwt
  end type sgld_ene_rec


  integer, parameter    :: sgld_ene_rec_size = 7

  type(sgld_ene_rec), parameter      :: null_sgld_ene_rec = &
       sgld_ene_rec(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)


!*******************************************************************************
!     head file for the self-guided langevin dynamics simulation  
!
!      variables for sgld simulation
!
! ... integers:
!
!
!    sgmd/sgld applying range
!      isgsta      begining atom index applying sgld
!      isgend      ending atom index applying sgld
!
!
!
! ... double precision:
!
!
!    sgmd/sgld variables
!     sgft    !  momentum guiding factor 
!     sgff    !  force guiding guiding factor 
!     sgfg    !  SGLD-GLE momentum guiding factor 
!     tsgavg  !  local average time, ps
!     tsgavp  !  convergence time, ps
!
!
!     sgavg0  !  local average remains
!     sgavg1  !  local average factor, sgavg1=1-sgavg0
!     sgavp0  !  convergence average remains
!     sgavp1  !  convergency average factor, sgavp1=1-sgavp0
!     templf  !  low frequency temperature, k
!     temphf  !  high frequency temperature, k
!     epotlf  !  low frequency potential energy
!     epothf  !  high frequency potential energy
!     epotllf !  low-low frequency potential energy
!*******************************************************************************



      double precision,save::  gammas, sggamma, &
          sgavg0,sgavg1,sgavp0,sgavp1, &
           sgwt, sgscale,tsgfac,sgrndf,  &
           sgmsum,com0sg(3),com1sg(3),com2sg(3)

      double precision  tsgset,&
          sgfti,sgffi,sgfgi, &
          templf,temphf,epotlf,epothf,epotllf,psgldg,fsgldg
         

!  common block for parallel broadcast
      common/sgldr/tsgset,&
          sgfti,sgffi,sgfgi,psgldg,fsgldg, &
          templf,temphf,epotlf,epothf,epotllf
          

!  number of broadcasting variables
integer, parameter :: nsgld_real=11
       
!*******************************************************************************
!
! ... flags:
!
!
!     tsgld      !  run sgmd/sgld  (isgld>0)
!     tsgldgle   !  sgld in generalized ensemble (isgld=2)
!     trxsgld   !  replica exchange sgld (isgld>0 & rem>0)
!     tsgbond   !  local average over bonded atoms
!     tsgmap    !  local average over atoms within sgsize distance
!
!*******************************************************************************

      logical, save :: tsgld,tsgldgle,trxsgld,tsgbond,tsgmap,tsggamma
      
!*******************************************************************************
!
! ...allocatable arrays:
!
!     avgx1     ! local averages of position
!     avgx2     ! local averages of local averages of position
!     avgp     ! local averages of momentum
!     avgr     ! local average of random forces
!     sgfps   !  averages of force-momentum product
!     sgpps   ! averages of momentum-momentum product
!*******************************************************************************

    double precision, dimension(:,:), allocatable, save :: avgx0,avgx1,avgx2,avgr
    double precision, dimension(:), allocatable, private,save :: sgfps,sgpps,sgmass
    integer, dimension(:), allocatable, private,save :: sgatoms
    type(sgld_ene_rec),private,save:: sgld_ene,sgld_sum,sgld_sum2
    type(sgld_ene_rec),private,save:: sgld_tmp,sgld_tmp2
    integer,private,save::  sgld_sumn,sgld_tmpn,nsgatom
    character(8), save:: sglabel

!*******************************************************************************
!
! spatial average variables:
!
!     atm_sg_maskdata     ! SG local structure of each atom (count, offset)
!     atm_sg_mask         ! atom id in SG local structures
!     NGRIDX,NGRIDY,NGRIDZ     ! SG map dimensions in x, y, z directions
!     NGRIDXY,NGRIDYZ,NGRIDZX     ! SG map dimensions in x-y, y-z, z-x plan
!     NGRIDXYZ     ! SG map dimensions in x-y-z
!     GXMIN,GXMAX,GYMIN,GYMAX,GZMIN,GZMAX     ! SG map boundary
!     RHOM   !  mass map
!     RHOVX,RHOVY,RHOVZ   ! velocity vx, vy, vz maps
!     RHOAX,RHOAY,RHOAZ   ! acceleration ax, ay, az maps
!*******************************************************************************
   type(listdata_rec), allocatable, save :: atm_sg_maskdata(:)
   integer, allocatable, save            :: atm_sg_mask(:)

   integer,private,save:: NGRIDX,NGRIDY,NGRIDZ,NGRIDXY,NGRIDYZ,NGRIDZX,NGRIDXYZ
   double precision,private,save::  GXMIN,GXMAX,GYMIN,GYMAX,GZMIN,GZMAX
   double precision, dimension(:), allocatable, private,save :: RHOM,RHOVX,RHOVY,RHOVZ,RHOAX,RHOAY,RHOAZ
 
contains

    subroutine sgld_init_ene_rec
!*********************************************************************
!               subroutine sgld_init_ene_rec
!*********************************************************************
!-----------------------------------------------------------------------
!     this routine performs initiation of energy record for the self-guided        
!       langevin dynamcs (sgld) simulaiton                 
!
!*********************************************************************
     implicit none

     sgld_ene = null_sgld_ene_rec
     sgld_sum = null_sgld_ene_rec
     sgld_sum2 = null_sgld_ene_rec
     sgld_sumn = 0
     sgld_tmp = null_sgld_ene_rec
     sgld_tmp2 = null_sgld_ene_rec
     sgld_tmpn = 0
     return
    end subroutine sgld_init_ene_rec

#ifdef MPI
    subroutine psgld(natom,amass,crd,vel,rem,my_atm_lst)
#else /* MPI */
    subroutine psgld(natom,amass,crd,vel)
#endif /* MPI */
!*********************************************************************
!               subroutine psgld
!*********************************************************************
!-----------------------------------------------------------------------
!     this routine performs initiation for the self-guided        
!       langevin dynamcs (sgld) simulaiton                 
!
!*********************************************************************
  use mdin_ctrl_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod, only: nres, &
           gbl_res_atms,gbl_labres,atm_igraph,atm_mass,atm_isymbl,  &
           bonda_idx,gbl_bond,atm_numex,gbl_natex,next,next_mult_fac
  use constraints_mod, only : atm_igroup
  use degcnt_mod
  use findmask_mod,only:atommask

      implicit none

      integer natom
      double precision      :: amass(natom)
      double precision      :: crd(3, natom)
      double precision      :: vel(3, natom)
      integer i,m,ierror,j,nsgsubi,idx_nbex,jatm
      double precision amassi,xi3,vi3,gamm,ekin,ekinsg,fact1,fact2
      logical is_langevin  ! is this a langevin dynamics simulation
#ifdef MPI
      integer               :: my_atm_lst(*)
      integer rem,atm_lst_idx,ierr
!
#ifdef VSCODE
!     mpi debug
      i=1
      do while (i==0)
        call sleep(5)
      enddo
#endif /* VSCODE */
#endif /* MPI */
!
      is_langevin = gamma_ln > 0.0d0
      tsgld = (isgld > 0)
      tsggamma = (isgld == 1)
      tsgbond = (sgtype > 1).and.(sgtype < 4)
      tsgmap = (sgtype > 3)
      tsgldgle = (abs(sgfg) > 1.0d-6).and.is_langevin
#ifdef MPI
      trxsgld = rem > 0 .and. tsgld
#endif /* MPI */
!  check for invalid sgld setting
      if(isgsta < 1)isgsta=1
      if(isgend > natom .or. isgend < 1)isgend=natom
      if(tsgavg.lt.dt)tsgavg=dt
      if(tsgavp.lt.dt)tsgavp=10.0d0*tsgavg
      sgavg1=dt/tsgavg
      sgavg0=1.0d0-sgavg1
      sgavp1=dt/tsgavp
      sgavp0=1.0d0-sgavp1
      is_langevin = gamma_ln > 0.0d0
      gammas=gamma_ln/20.455d0
      tsgfac=1.0d0/20.455d0/tsgavg
      sglabel="SGMD: "
      if(is_langevin)sglabel="SGLD: "
      tsgset=temp0
      if(tsgset<1.0d-6)tsgset=300.0d0
      if(isgld==3)then
        if(sgft<=1.0d0)then
          if(ABS(sgft)>ABS(sgff))then
            sgfti=sgft
            FACT1=9.0d0-SQRT(81.0d0-12.0d0*SGFTI*SGFTI*SGFTI)
            FACT2=(ABS(FACT1)*1.5d0)**(1.0/3.0d0)
            PSGLDG=SIGN(1.0d0,FACT1)*(FACT2/3.0d0+SGFTI/FACT2)-1.0d0
            sgffi=psgldg
          else
            sgffi=sgff
            sgfti=(1.0d0+SGFFI)*(1.0d0+SGFFI)-1.0d0/(1.0d0+SGFFI)
          endif
        else
          sgffi=sgff
          sgfti=(1.0d0+SGFFI)*(1.0d0+SGFFI)-1.0d0/(1.0d0+SGFFI)
        endif
      else 
        sgfti=sgft 
        sgffi=sgff
      endif
      if(ABS(SGFTI)>1.0d-8.and. sgfti<=1.0d0)then
        FACT1=9.0d0-SQRT(81.0d0-12.0d0*SGFTI*SGFTI*SGFTI)
        FACT2=(ABS(FACT1)*1.5d0)**(1.0/3.0d0)
        PSGLDG=SIGN(1.0d0,FACT1)*(FACT2/3.0d0+SGFTI/FACT2)-1.0d0
      else if(sgfti>-1.0d0)then 
        PSGLDG=SQRT(1.0d0+SGFTI)-1.0d0
      else
        PSGLDG=-1.0d0
      endif
      if(tsgldgle)then
        sgfgi=sgfg
        if(sgfgi<=1.0d0)then
          FSGLDG=SQRT(1.0d0-SGFGI)-1.0d0
        else 
          fsgldg=-1.0d0
        endif
        ! substract sgffi which contains random forces
        !fsgldg=fsgldg-sgffi
      else
        sgfgi=0.0d0
        fsgldg=0.0d0
      endif
      IF(TEMPSG>1.0d-6)THEN
        ! when tempsg is set
        FACT1=1.0d0-TSGSET/TEMPSG
        IF(SGFTI*SGFTI>1.0d-8)THEN
          SGFFI=PSGLDG-FACT1
        ELSE 
          PSGLDG=SGFFI+FACT1
          SGFTI=(1.0d0+PSGLDG)*(1.0d0+PSGLDG)-1.0d0/(1.0d0+PSGLDG)
        ENDIF
      ELSE
        TEMPSG=TSGSET/(1.0d0-PSGLDG+SGFFI)
      ENDIF
      ! identify guided atoms
#ifdef MPI
      if(numtasks>1)then
        call mpi_bcast(sgmask,256,MPI_CHARACTER,0,pmemd_comm,ierror)
      endif
#endif
      if(allocated(sgatoms))deallocate(sgatoms)
      allocate(sgatoms(NATOM),stat=ierror)
      if (ierror .ne. 0) call setup_alloc_error
      call atommask( natom, nres, 0, atm_igraph, atm_isymbl, &
      gbl_res_atms, gbl_labres, crd, sgmask, sgatoms )
    !     allocate working arrays
      allocate( avgx0(3,natom),avgx1(3,natom),avgx2(3,natom),  &
      avgr(3,natom ),sgfps(natom ),sgpps(natom ),&
      sgmass(natom),stat=ierror)
     if (ierror .ne. 0) call setup_alloc_error
     ! build bidireectional  exclusion lists
     if(tsgbond)then 
      allocate(atm_sg_maskdata(natom), &
               atm_sg_mask(next * next_mult_fac), &
               stat = ierror)
      if (ierror .ne. 0) call setup_alloc_error
      call make_sgavg_mask_list(natom, atm_numex, gbl_natex)
     endif
     if(tsgmap)then
       gxmin=1.0d8
       gxmax=-1.0d8
       gymin=1.0d8
       gymax=-1.0d8
       gzmin=1.0d8
       gzmax=-1.0d8
     endif
!    initialize arrays
     gamm=sqrt(dt/tsgavg)
     ekin=0.0d0
     ekinsg=0.0d0
     sgmsum=0.0d0
     nsgatom=0
      do i=1,natom
        do m=1,3
          xi3=crd(m,i)
          vi3=vel(m,i)
          amassi = amass(i)
          ekin=ekin+amassi*vi3*vi3
          if(tsgbond)then 
            nsgsubi=atm_sg_maskdata(i)%cnt
            xi3=amassi*xi3
            vi3=amassi*vi3
            do j=1,nsgsubi
              idx_nbex=atm_sg_maskdata(i)%offset + j 
              jatm=atm_sg_mask(idx_nbex)
              xi3=xi3+amass(jatm)*crd(m,jatm)
              vi3=vi3+amass(jatm)*vel(m,jatm)
              amassi=amassi+amass(jatm)
            enddo
            xi3=xi3/amassi
            vi3=vi3/amassi
          endif
          sgmass(i)=amassi
          amassi = amass(i)
          avgx0(m,i)=xi3
          avgx1(m,i)=xi3-vi3*sqrt(sgavg1)/tsgfac
          avgx2(m,i)=avgx1(m,i)-vi3*sgavg1/tsgfac
          avgr(m,i)=0.0d0
          ekinsg=ekinsg+amassi*vi3*vi3
        end do        
        sgpps(i)=amassi*0.001987*tsgset*sgavg1
        sgfps(i)=-0.01*sgpps(i)/20.455d0
        if(i>=isgsta.and.i<=isgend.and.sgatoms(i)>0)then
          nsgatom=nsgatom+1
          sgmsum=sgmsum+amassi
          sgatoms(i)=i
        else
          sgatoms(i)=0
        endif
      end do
      epotlf=2.0d10
      epotllf=2.0d10
      templf=tsgset*sgavg1
      temphf=tsgset-templf
      sgwt=0.0d0
      sggamma=0.01d0/20.455d0
      com0sg=0.0d0
      com1sg=0.0d0
      com2sg=0.0d0
#ifdef CUDA
     call gpu_sgld_setup(isgld,nsgatom,sgatoms, dt,tsgset,   &
       sgmsum,tsgavg,tsgavp,sgfti,sgffi,sgfgi,fsgldg,tempsg,  &
       avgx0,avgx1,avgx2,avgr,sgfps,sgpps)
#endif

#ifdef MPI
      if(master)then
#endif
        write(mdout,910)isgsta,isgend,nsgatom
        write(mdout,915)tsgavg,tsgavp
        if(is_langevin)then
          if(tsgldgle)then
            write(mdout,928)
            write(mdout,927)sgfgi,fsgldg
          else
            write(mdout,940)
          endif
          write(mdout,930)gamma_ln
        else
            write(mdout,941)
        endif
        write(mdout,925)sgfti,psgldg
        write(mdout,926)sgffi
        write(mdout,932)tempsg
        if(tsgbond)then
          if(sgtype==2)then
            write(mdout,942)
          else
            write(mdout,943)
          endif
        endif      
        if(tsgmap)then
          write(mdout,948)sgsize
        endif      
        write(mdout,935)
#ifdef MPI
      endif
#endif

910   format("  _________________ SGMD/SGLD parameters _________________"/  &
      "  Parameters for self-guided Molecular/Langevin dynamics (SGMD/SGLD) simulation"//  &
          "  Guiding range from ",i5,"  to ",i8, " with ",i8," guiding atoms")
915   format("  Local averaging time: tsgavg: ",f10.4," ps,  tsgavp: ",f10.4," ps")
925   format("  sgfti: ",f8.4," psgldg: ",f8.4)
926   format("  sgffi: ",f8.4)
927   format("  momentum factor sgfgi= ",f8.4," random force factor fsgldg=",f8.4)
928   format("  SGLD-GLE method is used to mantain a canonical distribution. ")
932   format("  Guided sampling effective temperature (TEMPSG): ",f8.2)
940   format("  SGLDg  method is used to enhance conformational search. ")
941   format("  SGMDg  method is used to enhance conformational search. ")
942   format("  SGTYPE=2, Guiding forces are averaged over 1-2,1-3 bonded structures" )
943   format("  SGTYPE=3, Guiding forces are averaged over 1-2,1-3,1-4 bonded structures" )
948   format("  SGTYPE=4, Guiding forces are averaged over a cutoff: ",F8.4 )
930   format("  Collision frequency:",f8.2," /ps" )
935   format("  Output properties:"    /  &
             "  SGMD/SGLD:  SGGAMMA TEMPLF  TEMPHF  EPOTLF EPOTHF EPOTLLF SGWT" /  &
             "         SGMD/SGLD weighting factor =exp(SGWT)"/  &
              " _______________________________________________________"/)
      return
    end subroutine psgld


subroutine make_sgavg_mask_list(atm_cnt, numex, natex)

  use gbl_constants_mod
  use file_io_dat_mod, only : mdout
  use constraints_mod, only : atm_igroup
  use mdin_ctrl_dat_mod, only : ibelly
  use pmemd_lib_mod
  use prmtop_dat_mod
  use extra_pnts_nb14_mod, only : gbl_nb14_cnt,gbl_nb14
  implicit none

! Formal arguments:

  integer               :: atm_cnt
  ! Excluded atom count for each atom.
  integer               :: numex(atm_cnt)
  ! Excluded atom concatenated list:
  integer               :: natex(next * next_mult_fac)

! Local variables:

  integer               :: atm_i, atm_j
  integer               :: lst_idx, sublst_idx, num_sublst
  integer               :: mask_idx
  integer               :: offset
  integer               :: total_excl
  ! nb14 variables
  integer               :: j,idx14,cnt14,list14(20)
  logical               :: addj

! Double the mask to deal with our list generator

! Pass 1: get pointers, check size

  lst_idx = 0

  atm_sg_maskdata(:)%cnt = 0       ! array assignment
  do atm_i = 1, atm_cnt - 1         ! last atom never has any...
    cnt14=0
    if(sgtype==2)then 
      do idx14=1,gbl_nb14_cnt
        if(gbl_nb14(1,idx14)==atm_i)then 
          cnt14=cnt14+1
          list14(cnt14)=gbl_nb14(2,idx14)
        endif
      end do
    endif
    num_sublst = numex(atm_i)
    do sublst_idx = 1, num_sublst
      atm_j = natex(lst_idx + sublst_idx)
      if (atm_j .gt. 0 ) then
        addj=.true.
        if(sgtype==2)then 
          do j=1,cnt14
            if(list14(j)==atm_j)addj=.false.
          enddo
        endif 
        if(addj)then
          atm_sg_maskdata(atm_i)%cnt = atm_sg_maskdata(atm_i)%cnt + 1
          atm_sg_maskdata(atm_j)%cnt = atm_sg_maskdata(atm_j)%cnt + 1
        endif 
      end if
    end do
    lst_idx = lst_idx + num_sublst
  end do

  total_excl = 0

  do atm_i = 1, atm_cnt
    total_excl = total_excl + atm_sg_maskdata(atm_i)%cnt
  end do

  if (total_excl .gt. next * next_mult_fac) then
    write(mdout, '(a,a)') error_hdr, &
         'The total number of sg substructure exceeds that stipulated by the'
    write(mdout, '(a,a)') error_hdr, &
         'prmtop.  This is likely due to a very high density of added extra points.'
    write(mdout, '(a,a)') error_hdr, &
         'Scale back the model detail, or contact the developers for a workaround.'
    call mexit(6, 1)
  end if

  offset = 0

  do atm_i = 1, atm_cnt
    atm_sg_maskdata(atm_i)%offset = offset
    offset = offset + atm_sg_maskdata(atm_i)%cnt
  end do

! Pass 2: fill mask array

  lst_idx = 0

  atm_sg_maskdata(:)%cnt = 0       ! array assignment
  
  if (ibelly .eq. 0) then
    do atm_i = 1, atm_cnt - 1
      if(sgtype==2)then 
        cnt14=0
        do idx14=1,gbl_nb14_cnt
          if(gbl_nb14(1,idx14)==atm_i)then 
            cnt14=cnt14+1
            list14(cnt14)=gbl_nb14(2,idx14)
          endif
        end do
      endif
      num_sublst = numex(atm_i)
      do sublst_idx = 1, num_sublst
        atm_j = natex(lst_idx + sublst_idx)
        if (atm_j .gt. 0 ) then
          addj=.true.
          if(sgtype==2)then 
            do j=1,cnt14
              if(list14(j)==atm_j)addj=.false.
            enddo
          endif 
          if(addj)then
            atm_sg_maskdata(atm_j)%cnt = atm_sg_maskdata(atm_j)%cnt + 1
            mask_idx = atm_sg_maskdata(atm_j)%offset + &
                     atm_sg_maskdata(atm_j)%cnt
            atm_sg_mask(mask_idx) = atm_i

            atm_sg_maskdata(atm_i)%cnt = atm_sg_maskdata(atm_i)%cnt + 1
            mask_idx = atm_sg_maskdata(atm_i)%offset + &
                     atm_sg_maskdata(atm_i)%cnt
            atm_sg_mask(mask_idx) = atm_j
          endif
        end if
      end do
      lst_idx = lst_idx + num_sublst
    end do
  else
    do atm_i = 1, atm_cnt - 1
      if(sgtype==2)then 
        cnt14=0
        do idx14=1,gbl_nb14_cnt
          if(gbl_nb14(1,idx14)==atm_i)then 
            cnt14=cnt14+1
            list14(cnt14)=gbl_nb14(2,idx14)
          endif
        end do
      endif
      num_sublst = numex(atm_i)
      do sublst_idx = 1, num_sublst
        atm_j = natex(lst_idx + sublst_idx)
        if (atm_j .gt. 0 ) then
          addj=.true.
          if(sgtype==2)then 
            do j=1,cnt14
              if(list14(j)==atm_j)addj=.false.
            enddo
          endif 
          if(addj)then
            if (atm_igroup(atm_i) .ne. 0 .and. atm_igroup(atm_j) .ne. 0) then

            atm_sg_maskdata(atm_j)%cnt = atm_sg_maskdata(atm_j)%cnt + 1
            mask_idx = atm_sg_maskdata(atm_j)%offset + &
                       atm_sg_maskdata(atm_j)%cnt
            atm_sg_mask(mask_idx) = atm_i

            atm_sg_maskdata(atm_i)%cnt = atm_sg_maskdata(atm_i)%cnt + 1
            mask_idx = atm_sg_maskdata(atm_i)%offset + &
            atm_sg_maskdata(atm_i)%cnt
            atm_sg_mask(mask_idx) = atm_j

            end if
          end if
        end if
      end do
      lst_idx = lst_idx + num_sublst
    end do
  end if
  return

end subroutine make_sgavg_mask_list


    subroutine sg_fix_degree_count(sgsta_rndfp, sgend_rndfp, ndfmin, rndf)
      !*********************************************************************
      !               subroutine sg_fix_degree_count
      !*********************************************************************
      !-----------------------------------------------------------------------
      !     correct the total number of degrees of freedom for a translatable COM,
      !       and compute the number of degrees of freedom in the sgld part.
      !       the latter is mostly done by the caller to avoid passing the long
      !       argument list needed by routine degcnt which also differs between
      !       sander and pmemd.
      !
      !*********************************************************************
            implicit none
            double precision, intent(in)    :: sgsta_rndfp, sgend_rndfp
            integer, intent(in)   :: ndfmin
            double precision, intent(inout) :: rndf
      
            sgrndf = sgend_rndfp - sgsta_rndfp
#ifdef CUDA
            call gpu_sgld_rndf(sgrndf)
#endif
            return
          end subroutine sg_fix_degree_count
      
#ifdef MPI
    subroutine sgldw(atm_cnt,my_atm_cnt, &
             dtx,pot_ene,amass,winv,crd,frc,vel, my_atm_lst)
#else
    subroutine sgldw(natom,atm_cnt, &
             dtx,pot_ene,amass,winv,crd,frc,vel)
#endif
!*********************************************************************
!               subroutine sgldw
!*********************************************************************
!-----------------------------------------------------------------------
!     perform sgld integration        
!
!*********************************************************************
  use random_mod
      implicit none
#ifdef MPI
      integer my_atm_cnt, ierr
#endif
      integer natom,atm_cnt
      double precision dtx,pot_ene
      double precision      :: amass(atm_cnt)
      double precision      :: winv(atm_cnt)
      double precision      :: crd(3, atm_cnt)
      double precision      :: vel(3, atm_cnt)
      double precision      :: frc(3, atm_cnt)
!
      integer i,m,jsta,jend,j,nsgsubi,idx_nbex,jatm
      double precision boltz,amassi,amassj,temp_t
      double precision fact,wfac,rsd,fln,gam,sgbeta
      double precision ekin,ekinsg
      double precision dcom1(3),dcom2(3),com0(3),com1(3),com2(3)
      double precision sumgam,sumfp,sumpp,sumgv,sumpv
      double precision sggammai,avgpi3,pi3t,avgdfi3,avgri3,fsgpi,fsgfi,fsgi3,frici
      double precision xi3,x1i3,x2i3,vi3t,vi3,fi3,x0i3,dx0i3,psgi3
      double precision avgpi(3),avgfi(3),temp1(12),temp2(12)
      parameter (boltz = 1.987192d-3)
#ifdef MPI
      integer               :: my_atm_lst(*)
      integer atm_lst_idx
#endif
!
    !
        if(isgld.eq.2)then
          dcom1=0.0d0
          dcom2=0.0d0
        else
          dcom1=com0sg-com1sg
          dcom2=com0sg-com2sg
        endif
    ! build spatial average maps
    if(tsgmap)then
#ifdef MPI
      call mapbuild(my_atm_cnt,dtx,amass,crd,vel,dcom1,dcom2,my_atm_lst)
#else
      call mapbuild(atm_cnt,dtx,amass,crd,vel,dcom1,dcom2)
#endif
    endif
    !
        gam=gammas*dtx
        jsta=isgsta
        jend=isgend
        sumgam=0.0d0
        ekin=0.0d0
        ekinsg=0.0d0
        com0=0.0d0
        com1=0.0d0
        com2=0.0d0
#ifdef MPI
        !if(new_list)call sg_allgather(atm_cnt)
        do i = 1, atm_cnt
          if (gbl_atm_owner_map(i) .ne. mytaskid) then
            if(no_ntt3_sync)cycle
            call gauss( 0.d0, 1.0d0, fln )
            call gauss( 0.d0, 1.0d0, fln )
            call gauss( 0.d0, 1.0d0, fln )
            cycle
          endif
#else
        do i=1,atm_cnt
#endif
          amassi = amass(i)
          wfac =  2.0d0*dtx*winv(i)
          rsd = sqrt(2.d0*gammas*boltz*tsgset*amassi/dtx)
          if(sgatoms(i)>0)then
            ! sggamma
            sggammai=-sgfps(i)/sgpps(i)
            sumgam=sumgam+sggammai
            if(tsggamma)sggammai=sggamma
            if(tsgmap)then
              ! spatial average using map
              call mapvalues(amassi,crd(:,i),avgpi,avgfi)
            endif
            sumfp=0.0d0
            sumpp=0.0d0
            sumgv=0.0d0
            sumpv=0.0d0
            do  m = 1,3
              !   generate random number 
              call gauss( 0.d0, rsd, fln )
              if(tsgmap)then
              ! spatial average using map
                avgpi3=avgpi(m)
                avgdfi3=avgfi(m)
              else
                ! avg(x)
                xi3=crd(m,i)
                if(tsgbond)then 
                  ! average over SG structures
                  nsgsubi=atm_sg_maskdata(i)%cnt
                  xi3=amassi*xi3
                  do j=1,nsgsubi
                    idx_nbex=atm_sg_maskdata(i)%offset + j 
                    jatm=atm_sg_mask(idx_nbex)
                    xi3=xi3+amass(jatm)*crd(m,jatm)
                  enddo
                  xi3=xi3/sgmass(i)
                endif
                ! adjust previous averages
                x0i3=xi3-vel(m,i)*dtx
                dx0i3=x0i3-avgx0(m,i)
                x1i3=avgx1(m,i)+dx0i3
                x2i3=avgx2(m,i)+dx0i3
                psgi3=tsgfac*amassi*(x0i3-x1i3)
                avgx0(m,i)=xi3
                ! Obtain local averages
                x1i3=sgavg0*(x1i3+dcom1(m))+sgavg1*xi3
                avgx1(m,i)=x1i3
                ! avgavg(x)
                x2i3=sgavg0*(x2i3+dcom2(m))+sgavg1*x1i3
                avgx2(m,i)=x2i3
                ! avg(p)
                avgpi3=tsgfac*amassi*(xi3-x1i3)
                pi3t=(avgpi3-sgavg0*psgi3)/sgavg1
                ! avg(f-avg(f))
                avgdfi3=tsgfac*(pi3t-2.0d0*avgpi3+tsgfac*amassi*(x1i3-x2i3))
                com0(m)=com0(m)+amassi*xi3
                com1(m)=com1(m)+amassi*x1i3
                com2(m)=com2(m)+amassi*x2i3
              endif /* tsgmap */
              ! sum(avg(f-avg(f))avg(p))
              sumfp=sumfp+avgdfi3*avgpi3
              ! sum(avg(p)avg(p))
              sumpp=sumpp+avgpi3*avgpi3
              ! average random forces
              avgri3=sgavg0*avgr(m,i)+sgavg1*fln
              avgr(m,i)=avgri3
              ! guiding forces
              fsgpi=sgfti*sggammai*avgpi3
              fsgfi=sgffi*avgdfi3
              fsgi3=sgfgi*gammas*avgpi3+(fsgldg-sgffi)*avgri3+fsgpi+fsgfi
              fi3=frc(m,i)+fln+fsgi3
              frc(m,i)=fi3
              ! estimate velocity at t+dt/2
              ! Using volocities at t avoid SHAKE complication
              vi3t=vel(m,i)
              ! sum(g*v)
              sumgv=sumgv+fsgi3*vi3t
              ! sum(p*v)
              sumpv=sumpv+amassi*vi3t*vi3t
              ekin=ekin+amassi*vi3t*vi3t
              ekinsg=ekinsg+avgpi3*avgpi3/amassi
            end do
            ! <(avg(f-avg(f))avg(v))>
            sgfps(i)=sgavp0*sgfps(i)+sgavp1*sumfp
            ! <(avg(p)avg(v))>
            sgpps(i)=sgavp0*sgpps(i)+sgavp1*sumpp
            ! energy conservation friction constant
            if(sumpv<1.0d-8)then
              sgbeta=0.0d0
            else
              sgbeta=sumgv/(2.0d0*sumpv/(2.0d0+gam)-sumgv*dtx)
            endif
            fact=dtx*(gammas+sgbeta)
            !write(*,'("DEBUG: ",i6,10f10.4)')i,crd(1,i),avgx1(1,i),avgx2(1,i),vel(1,i),crd(2,i),vel(3,i),sgbeta,fact,sumfp
            do  m = 1,3
              fi3=frc(m,i)
              vi3t=((2.0d0-fact)*vel(m,i)+fi3*wfac)/(2.0d0+fact)
              vel(m,i)=vi3t
           end do
          !write(*,'("DEBUG: ",i6,10f10.4)')i,crd(1,i),avgx1(1,i),avgx2(1,i),avgp(1,i)/amassi,avgp(2,i)/amassi,avgp(3,i)/amassi,sgbeta,fact,sumfp
           !if(i<20)write(*,'("DEBUG: ",i6,10f10.4)')i,crd(1,i),avgx0(1,i),avgx1(1,i),avgx2(1,i),sgbeta,fact,sumfp,sggammai
          else 
              ! without guiding forces
            do  m = 1,3
                !   generate random number 
                call gauss( 0.d0, rsd, fln )
                fi3=frc(m,i)+fln
                frc(m,i)=fi3
                vi3t=((2.0d0-gam)*vel(m,i)+fi3*wfac)/(2.0d0+gam)
                vel(m,i)=vi3t
              end do
          endif
        end do
#ifdef MPI
        if(numtasks > 1)then
!  combining all node results
!
          temp1(1)=sumgam
          temp1(2)=ekin
          temp1(3)=ekinsg
          temp1(4:6)=com0
          temp1(7:9)=com1
          temp1(10:12)=com2

          call mpi_allreduce(temp1,temp2,12, &
          mpi_double_precision,mpi_sum,pmemd_comm,ierr)

          sumgam=temp2(1)
          ekin=temp2(2)
          ekinsg=temp2(3)
          com0=temp2(4:6)
          com1=temp2(7:9)
          com2=temp2(10:12)

        endif
#endif
    ! Estimate low frequency temperatures
        !temp_t=ekin/sgrndf/boltz
        temp_t=tsgset
        !templf=sgavp0*templf+sgavp1*ekinsg/sgrndf/boltz
        templf=sgavp0*templf+sgavp1*ekinsg*tsgset/ekin
        temphf=temp_t-templf
        sggamma=sumgam/nsgatom
        sgscale=20.455d0*sggamma
        com0sg=com0/sgmsum
        com1sg=com1/sgmsum
        com2sg=com2/sgmsum
        !write(*,'("DEBUG: ",i6,f10.2,11f10.4)')i,sgmsum,sggamma,com0sg,com1sg,com2sg
        ! update accumulators
        call sgenergy(pot_ene)
        return
    end subroutine sgldw

#ifdef MPI
        subroutine sgmdw(atm_cnt,my_atm_cnt, &
             dtx,pot_ene,amass,winv,crd,frc,vel,my_atm_lst)
#else
        subroutine sgmdw(natom,atm_cnt, &
             dtx,pot_ene,amass,winv,crd,frc,vel)
#endif
!*********************************************************************
!               subroutine sgmdw
!*********************************************************************
!-----------------------------------------------------------------------
!     calculate guiding force using sgld method for md simulation        
!
!*********************************************************************
  use random_mod
      implicit none
#ifdef MPI
      integer my_atm_cnt,ierr
#endif
      integer natom,atm_cnt
      double precision dtx,pot_ene
      double precision      :: amass(atm_cnt)
      double precision      :: winv(atm_cnt)
      double precision      :: crd(3, atm_cnt)
      double precision      :: vel(3, atm_cnt)
      double precision      :: frc(3, atm_cnt)
!
#ifdef MPI
      integer               :: my_atm_lst(*)
      integer atm_lst_idx
#endif
      integer jsta,jend,i,m,j,nsgsubi,idx_nbex,jatm
      double precision boltz,amassi,temp_t
      double precision fact,wfac,ekin,ekinsg,sgbeta
      double precision dcom1(3),dcom2(3),com0(3),com1(3),com2(3)
      double precision sumgam,sumfp,sumpp,sumgv,sumpv
      double precision sggammai,avgpi3,pi3t,avgdfi3,avgri3,fsgpi,fsgfi,fsgi3,frici
      double precision xi3,x1i3,x2i3,vi3t,vi3,fi3,x0i3,dx0i3,psgi3
      double precision avgpi(3),avgfi(3),temp1(12),temp2(12)
      parameter (boltz = 1.987192d-3)
!
    !
        if(isgld.eq.2)then
          dcom1=0.0d0
          dcom2=0.0d0
        else
          dcom1=com0sg-com1sg
          dcom2=com0sg-com2sg
        endif
    ! build spatial average maps
      if(tsgmap)then
#ifdef MPI
        call mapbuild(my_atm_cnt,dtx,amass,crd,vel,dcom1,dcom2,my_atm_lst)
#else
        call mapbuild(atm_cnt,dtx,amass,crd,vel,dcom1,dcom2)
#endif
      endif
            !
        jsta=isgsta
        jend=isgend
        sumgam=0.0d0
        ekin=0.0d0
        ekinsg=0.0d0
        com0=0.0d0
        com1=0.0d0
        com2=0.0d0
#ifdef MPI
        do atm_lst_idx = 1, my_atm_cnt
        i = my_atm_lst(atm_lst_idx)
        !if (gbl_atm_owner_map(i) .ne. mytaskid) cycle
#else
        do i = 1, atm_cnt
#endif
          if(sgatoms(i)>0)then
            wfac = 0.5d0*winv(i)*dtx
            amassi=amass(i)
            ! sggamma
            sggammai=-sgfps(i)/sgpps(i)
            sumgam=sumgam+sggammai
            if(tsggamma)sggammai=sggamma
            if(tsgmap)then
              ! spatial average using map
              call mapvalues(amassi,crd(:,i),avgpi,avgfi)
            endif

            sumfp=0.0d0
            sumpp=0.0d0
            sumgv=0.0d0
            sumpv=0.0d0
            do  m = 1,3
              if(tsgmap)then
                ! spatial average using map
                  avgpi3=avgpi(m)
                  avgdfi3=avgfi(m)
              else
                ! avg(x)
                xi3=crd(m,i)
                if(tsgbond)then 
                  nsgsubi=atm_sg_maskdata(i)%cnt
                  xi3=amassi*xi3
                  do j=1,nsgsubi
                    idx_nbex=atm_sg_maskdata(i)%offset + j 
                    jatm=atm_sg_mask(idx_nbex)
                    xi3=xi3+amass(jatm)*crd(m,jatm)
                  enddo
                  xi3=xi3/sgmass(i)
                endif
                ! adjust previous averages
                x0i3=xi3-vel(m,i)*dtx
                dx0i3=x0i3-avgx0(m,i)
                x1i3=avgx1(m,i)+dx0i3
                x2i3=avgx2(m,i)+dx0i3
                psgi3=tsgfac*amassi*(x0i3-x1i3)
                avgx0(m,i)=xi3
                ! Obtain local averages
                x1i3=sgavg0*(x1i3+dcom1(m))+sgavg1*xi3
                avgx1(m,i)=x1i3
                ! avgavg(x)
                x2i3=sgavg0*(x2i3+dcom2(m))+sgavg1*x1i3
                avgx2(m,i)=x2i3
                ! avg(p)
                avgpi3=tsgfac*amassi*(xi3-x1i3)
                pi3t=(avgpi3-sgavg0*psgi3)/sgavg1
                ! avg(f-avg(f))
                avgdfi3=tsgfac*(pi3t-2.0d0*avgpi3+tsgfac*amassi*(x1i3-x2i3))
                com0(m)=com0(m)+amassi*xi3
                com1(m)=com1(m)+amassi*x1i3
                com2(m)=com2(m)+amassi*x2i3
              endif /* tsgmap */
            ! sum(avg(f-avg(f))avg(p))
              sumfp=sumfp+avgdfi3*avgpi3
              ! sum(avg(p)avg(p))
              sumpp=sumpp+avgpi3*avgpi3
              ! guiding forces
              fsgpi=sgfti*sggammai*avgpi3
              fsgfi=sgffi*avgdfi3
              fsgi3=fsgpi+fsgfi
              fi3=frc(m,i)+fsgi3
              frc(m,i)=fi3
             ! estimate velocity at t+dt/2
             !vi3t=vel(m,i)+fi3*wfac
             ! Using volocities at t avoid SHAKE complication
              vi3t=vel(m,i)
              ! sum(g*v)
              sumgv=sumgv+fsgi3*vi3t
              ! sum(p*v)
              sumpv=sumpv+amassi*vi3t*vi3t
              ekin=ekin+amassi*vi3t*vi3t
              ekinsg=ekinsg+avgpi3*avgpi3/amassi
            end do
            ! <(avg(f-avg(f))avg(v))>
            sgfps(i)=sgavp0*sgfps(i)+sgavp1*sumfp
            ! <(avg(p)avg(v))>
            sgpps(i)=sgavp0*sgpps(i)+sgavp1*sumpp
            ! energy conservation friction constant
            if(sumpv<1.0d-8)then
              sgbeta=0.0d0
            else
              sgbeta=2.0d0*sumgv/(2.0d0*sumpv-sumgv*dtx)
            endif
            fact=sgbeta/(1.0d0+0.5d0*sgbeta*dtx)
            do  m = 1,3
              fi3=frc(m,i)
              vi3t = vel(m,i) + fi3*wfac
              frici=fact*amassi*vi3t
              frc(m,i)=fi3-frici
            end do
          endif
          !if(i>natom-21)write(*,'("DEBUG: ",i6,10f10.4)')i,sggammai,crd(1,i),avgx0(1,i),avgx1(1,i),avgx2(1,i),sgfps(i),sgpps(i)
        end do
#ifdef MPI
        if(numtasks > 1)then
!  combining all node results
!
          temp1(1)=sumgam
          temp1(2)=ekin
          temp1(3)=ekinsg
          temp1(4:6)=com0
          temp1(7:9)=com1
          temp1(10:12)=com2

          call mpi_allreduce(temp1,temp2,12, &
          mpi_double_precision,mpi_sum,pmemd_comm,ierr)

          sumgam=temp2(1)
          ekin=temp2(2)
          ekinsg=temp2(3)
          com0=temp2(4:6)
          com1=temp2(7:9)
          com2=temp2(10:12)
        endif
#endif
    ! Estimate low frequency temperatures
        !tempi=ekin/sgrndf/boltz
        tempi=tsgset
        !templf=sgavp0*templf+sgavp1*ekinsg/sgrndf/boltz
        templf=sgavp0*templf+sgavp1*ekinsg*tsgset/ekin
        temphf=tempi-templf
        sggamma=sumgam/nsgatom
        sgscale=20.455d0*sggamma
        com0sg=com0/sgmsum
        com1sg=com1/sgmsum
        com2sg=com2/sgmsum
        !write(*,'("DEBUG: ",i6,f10.2,12f10.4)')i,sgmsum,tempi,templf,sggamma,com0sg,com1sg,com2sg
       ! update accumulators
        call sgenergy(pot_ene)
        return
      end subroutine sgmdw
     
      subroutine sgenergy(pot_ene)
!*********************************************************************
!               subroutine sgenergy
!*********************************************************************
!-----------------------------------------------------------------------
!     set the sgld_ene fields of sgld variables
!
!*********************************************************************
      implicit none
      double precision pot_ene
      double precision boltz
      parameter (boltz = 1.987192d-3)
      double precision epoti
    ! weighting accumulators
        epoti=pot_ene
        if(epotlf>1.0d10)then
          epotlf=epoti
          epotllf=epotlf
        else
          epotlf=sgavg0*epotlf+sgavg1*epoti
          epotllf=sgavg0*epotllf+sgavg1*epotlf
        endif
        epothf=epoti-epotlf
        sgwt=(psgldg-sgffi)*(epotlf-epotllf)/(boltz*tsgset)
    ! update ener structure
        sgld_ene%sggamma=sgscale
        sgld_ene%templf=templf
        sgld_ene%temphf=temphf
        sgld_ene%epotlf=epotlf
        sgld_ene%epothf=epothf
        sgld_ene%epotllf=epotllf
        sgld_ene%sgwt=sgwt
!     statistics of sgld_ene
        sgld_sumn=sgld_sumn+1
   ! averages
        sgld_sum%sggamma=sgld_sum%sggamma+sgld_ene%sggamma
        sgld_sum%templf=sgld_sum%templf+sgld_ene%templf
        sgld_sum%temphf=sgld_sum%temphf+sgld_ene%temphf
        sgld_sum%epotlf=sgld_sum%epotlf+sgld_ene%epotlf
        sgld_sum%epothf=sgld_sum%epothf+sgld_ene%epothf
        sgld_sum%epotllf=sgld_sum%epotllf+sgld_ene%epotllf
        sgld_sum%sgwt=sgld_sum%sgwt+sgld_ene%sgwt
   ! fluctuations
        sgld_sum2%sggamma=sgld_sum2%sggamma+sgld_ene%sggamma**2
        sgld_sum2%templf=sgld_sum2%templf+sgld_ene%templf**2
        sgld_sum2%temphf=sgld_sum2%temphf+sgld_ene%temphf**2
        sgld_sum2%epotlf=sgld_sum2%epotlf+sgld_ene%epotlf**2
        sgld_sum2%epothf=sgld_sum2%epothf+sgld_ene%epothf**2
        sgld_sum2%epotllf=sgld_sum2%epotllf+sgld_ene%epotllf**2
        sgld_sum2%sgwt=sgld_sum2%sgwt+sgld_ene%sgwt**2

      return
      end subroutine sgenergy

#ifdef CUDA
      subroutine gpu_sgenergy(pot_ene)
      implicit none
      double precision pot_ene
      ! Calculate SGLD quantities
      call gpu_calculate_sgld_averages(dt, temp0, gamma_ln, templf, temphf, sggamma,com0sg,com1sg,com2sg)
      sgscale=20.455d0*sggamma
      call sgenergy(pot_ene)
      return
      end subroutine gpu_sgenergy
#endif

      subroutine sgld_avg(total)
!*********************************************************************
!               subroutine sgld_avg
!*********************************************************************
!-----------------------------------------------------------------------
!     average the sgld data fields for output.
!
!*********************************************************************
      implicit none
      logical total

      if(total)then
        sgld_tmp=sgld_sum
        sgld_tmpn=sgld_sumn
      else
        sgld_tmp%sggamma=sgld_sum%sggamma-sgld_tmp%sggamma
        sgld_tmp%templf=sgld_sum%templf-sgld_tmp%templf
        sgld_tmp%temphf=sgld_sum%temphf-sgld_tmp%temphf
        sgld_tmp%epotlf=sgld_sum%epotlf-sgld_tmp%epotlf
        sgld_tmp%epothf=sgld_sum%epothf-sgld_tmp%epothf
        sgld_tmp%epotllf=sgld_sum%epotllf-sgld_tmp%epotllf
        sgld_tmp%sgwt=sgld_sum%sgwt-sgld_tmp%sgwt
        sgld_tmpn=sgld_sumn-sgld_tmpn
      endif
        sgld_ene%sggamma=sgld_tmp%sggamma/sgld_tmpn
        sgld_ene%templf=sgld_tmp%templf/sgld_tmpn
        sgld_ene%temphf=sgld_tmp%temphf/sgld_tmpn
        sgld_ene%epotlf=sgld_tmp%epotlf/sgld_tmpn
        sgld_ene%epothf=sgld_tmp%epothf/sgld_tmpn
        sgld_ene%epotllf=sgld_tmp%epotllf/sgld_tmpn
        sgld_ene%sgwt=sgld_tmp%sgwt/sgld_tmpn
      return
      end subroutine sgld_avg

      subroutine sgld_fluc(total)
!*********************************************************************
!               subroutine sgld_fluc
!*********************************************************************
!-----------------------------------------------------------------------
!     fluctuation of the sgld data fields for output.
!
!*********************************************************************
      implicit none
      logical total

      if(total)then
        sgld_tmp2=sgld_sum2
      else
        sgld_tmp2%sggamma=sgld_sum2%sggamma-sgld_tmp2%sggamma
        sgld_tmp2%templf=sgld_sum2%templf-sgld_tmp2%templf
        sgld_tmp2%temphf=sgld_sum2%temphf-sgld_tmp2%temphf
        sgld_tmp2%epotlf=sgld_sum2%epotlf-sgld_tmp2%epotlf
        sgld_tmp2%epothf=sgld_sum2%epothf-sgld_tmp2%epothf
        sgld_tmp2%epotllf=sgld_sum2%epotllf-sgld_tmp2%epotllf
        sgld_tmp2%sgwt=sgld_sum2%sgwt-sgld_tmp2%sgwt
      endif

      sgld_ene%sggamma=(sgld_tmp2%sggamma-sgld_tmp%sggamma**2/sgld_tmpn)/sgld_tmpn
      sgld_ene%templf=(sgld_tmp2%templf-sgld_tmp%templf**2/sgld_tmpn)/sgld_tmpn
      sgld_ene%temphf=(sgld_tmp2%temphf-sgld_tmp%temphf**2/sgld_tmpn)/sgld_tmpn
      sgld_ene%epotlf=(sgld_tmp2%epotlf-sgld_tmp%epotlf**2/sgld_tmpn)/sgld_tmpn
      sgld_ene%epothf=(sgld_tmp2%epothf-sgld_tmp%epothf**2/sgld_tmpn)/sgld_tmpn
      sgld_ene%epotllf=(sgld_tmp2%epotllf-sgld_tmp%epotllf**2/sgld_tmpn)/sgld_tmpn
      sgld_ene%sgwt=(sgld_tmp2%sgwt-sgld_tmp%sgwt**2/sgld_tmpn)/sgld_tmpn
      if(sgld_ene%sggamma>0.0d0)sgld_ene%sggamma=sqrt(sgld_ene%sggamma)
      if(sgld_ene%templf>0.0d0)sgld_ene%templf=sqrt(sgld_ene%templf)
      if(sgld_ene%temphf>0.0d0)sgld_ene%temphf=sqrt(sgld_ene%temphf)
      if(sgld_ene%epotlf>0.0d0)sgld_ene%epotlf=sqrt(sgld_ene%epotlf)
      if(sgld_ene%epothf>0.0d0)sgld_ene%epothf=sqrt(sgld_ene%epothf)
      if(sgld_ene%epotllf>0.0d0)sgld_ene%epotllf=sqrt(sgld_ene%epotllf)
      if(sgld_ene%sgwt>0.0d0)sgld_ene%sgwt=sqrt(sgld_ene%sgwt)
      sgld_tmp=sgld_sum
      sgld_tmp2=sgld_sum2
      sgld_tmpn=sgld_sumn
      return
      end subroutine sgld_fluc

      subroutine sgld_print( iout)
!*********************************************************************
!               subroutine sgld_print
!*********************************************************************
!-----------------------------------------------------------------------
!     print the sgld data fields during the per nstep output.
!
!*********************************************************************
      implicit none
      integer, intent(in) :: iout      ! The output unit.

      write(iout,1005)sglabel,sgld_ene%sggamma,sgld_ene%templf,sgld_ene%temphf,  &
      sgld_ene%epotlf,sgld_ene%epothf,sgld_ene%epotllf,sgld_ene%sgwt

   1005 format(1X,A6,F9.4,F8.2,F8.2,X,F12.2,F12.2,F12.2,F10.4)

      return
      end subroutine sgld_print

#ifdef MPI

!*********************************************************************
!               subroutine sg_allgather
!*********************************************************************
!  Gather  all SG vectors when new list is set
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine sg_allgather(atm_cnt,using_pme_potential,using_gb_potential)

  use parallel_mod
  use gb_parallel_mod

  implicit none

   integer, intent(in) :: atm_cnt
   logical using_pme_potential,using_gb_potential
   integer  i
   double precision:: buffer(3,atm_cnt)
   do i=1,atm_cnt
    buffer(1,i)=sgfps(i)
    buffer(2,i)=sgpps(i)
    buffer(3,i)=sgmass(i)
   enddo
   if (using_pme_potential) then
    call mpi_allgathervec(atm_cnt, avgx0)
    call mpi_allgathervec(atm_cnt, avgx1)
    call mpi_allgathervec(atm_cnt, avgx2)
    call mpi_allgathervec(atm_cnt, avgr)
    call mpi_allgathervec(atm_cnt, buffer)
   else if (using_gb_potential) then
    call gb_mpi_allgathervec(atm_cnt, avgx0)
    call gb_mpi_allgathervec(atm_cnt, avgx1)
    call gb_mpi_allgathervec(atm_cnt, avgx2)
    call gb_mpi_allgathervec(atm_cnt, avgr)
    call gb_mpi_allgathervec(atm_cnt, buffer)
   end if
   do i=1,atm_cnt
    sgfps(i)=buffer(1,i)
    sgpps(i)=buffer(2,i)
    sgmass(i)=buffer(3,i)
   enddo
   return
   end subroutine sg_allgather

!*********************************************************************
!               subroutine rxsgld_scale
!*********************************************************************
! scale sgld quantites  after exchange
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine rxsgld_scale(atm_cnt,my_atm_lst,myscaling,myscalsg,amass,crd,vel)

   implicit none

   integer, intent(in) :: atm_cnt
   integer               :: my_atm_lst(*)
   double precision, intent(in) :: myscaling,myscalsg
   double precision, intent(inout) :: crd(3,atm_cnt)
   double precision, intent(inout) :: vel(3,atm_cnt)
   double precision, intent(in)    :: amass(atm_cnt)
   integer ierror
   integer i,j,jsta,jend,atm_lst_idx
   double precision xi,x0i,x1i,x2i,amassi
!--------------------
         ! all processes scale velocities, avgx1, avgx2, vsg.
         ! dan roe: this could potentially be divided up as in runmd
         !  since when there are mutiple threads per group each thread 
         !  only ever knows about its own subset of velocities anyway.
#ifdef VERBOSE_REMD
         if (master) then
            write (mdout,'(a,f8.4,a,f8.4,a,3f8.4)') &
               "rxsgld: scaling velocities by ", myscaling, &
               " guiding properties by ",myscalsg,&
               " to match temp0, sgft, sgff: ",temp0,sgfti,sgffi
            !write (mdout,*) &
            !   "atm_cnt: ",atm_cnt,isgsta,isgend,numtasks
         endif
#endif
! ---=== broadcast rxsgld guiding effect ===---
      if(numtasks > 1)call mpi_bcast(tsgset,nsgld_real,mpi_double_precision,&
                              0,pmemd_comm,ierror)
#ifdef CUDA
      if(myscaling>0.0d0 .and. myscalsg>0.0d0)then
        call gpu_scale_sgld(myscaling,myscalsg,sgfti,sgffi,tsgset)
        templf=myscalsg*myscalsg*templf
      endif                         
#else
      if (myscaling > 0.0d0) then
            vel(:,:)=myscaling*vel(:,:)
      endif                         
      if (myscalsg > 0.0d0) then
        templf=myscalsg*myscalsg*templf
        avgr(:,:)=myscalsg*avgr(:,:)
        sgfps(:)=myscalsg*sgfps(:)
        sgpps(:)=myscalsg*myscalsg*sgpps(:)
  !ROMEEE substitute iparpt by the right variable and removelines and uncomment right ones.
!        jsta = iparpt(mytaskid) + 1
!        jend = iparpt(mytaskid+1)
        jsta = 1
        jend = atm_cnt
        !if(jsta < isgsta)jsta=isgsta
        !if(jend > isgend)jend=isgend
#if defined(MPI) && !defined(CUDA)
         do atm_lst_idx = 1, atm_cnt
           i = my_atm_lst(atm_lst_idx)
#else
         do i = 1, atm_cnt
#endif
            if(i<isgsta.or.i>isgend)cycle
            amassi=amass(i)
            do j=1,3
               xi=crd(j,i)
               x0i=avgx0(j,i)
               x1i=avgx1(j,i)
               x2i=avgx2(j,i)
               avgx1(j,i)=xi+myscalsg*(x1i-x0i)
               avgx2(j,i)=xi+myscalsg*(x2i-x0i)
           enddo
         enddo
      endif
#endif
   return

end subroutine rxsgld_scale


#endif /* MPI */

#ifdef MPI
  SUBROUTINE MAPBUILD(my_atm_cnt,DELTA,AMASS,CRD,VEL,DCOM1,DCOM2, my_atm_lst)
#else
  SUBROUTINE MAPBUILD(atm_cnt,DELTA,AMASS,CRD,VEL,DCOM1,DCOM2)
#endif
  !-----------------------------------------------------------------------
  !     This routine build velocity and acceleration maps
  !
  !-----------------------------------------------------------------------
  use mdin_ctrl_dat_mod,only:sgsize
  
  implicit none
  INTEGER atm_cnt
  double precision CRD(3,*),VEL(3,*),AMASS(*),DCOM1(3),DCOM2(3)
  double precision DELTA,XI,YI,ZI,AMASSI
  INTEGER X0,Y0,Z0,X1,Y1,Z1
  INTEGER I000,I100,I010,I001,I110,I101,I011,I111
  double precision GX,GY,GZ,A0,B0,C0,A1,B1,C1
  double precision ABC000,ABC100,ABC010,ABC001,ABC110,ABC101,ABC011,ABC111
  INTEGER IA,I
    double precision X1I,Y1I,Z1I,X2I,Y2I,Z2I,XT,YT,ZT
    double precision GX0,GY0,GZ0,GXI,GYI,GZI,PXT,PYT,PZT,FXI,FYI,FZI
    double precision FACT1,FACT2,RHOMI
#ifdef MPI
      integer               :: my_atm_lst(*)
      integer my_atm_cnt,atm_lst_idx,ierr
#endif
  !
  ! check mapsize
#ifdef MPI
  CALL MAPINIT(my_atm_cnt,CRD, my_atm_lst)
#else
  CALL MAPINIT(atm_cnt,CRD)
#endif
  ! interpolate structure to protein map
  !WXW Calculate constraint factor to eliminate energy input from guiding force
#ifdef MPI
        do atm_lst_idx = 1, my_atm_cnt
        i = my_atm_lst(atm_lst_idx)
        !if (gbl_atm_owner_map(i) .ne. mytaskid) cycle
#else
        do i = 1, atm_cnt
#endif
          if(i>=isgsta.and.i<=isgend.and.sgatoms(i)>0)then

            AMASSI=AMASS(I)
            XI=CRD(1,I)
            YI=CRD(2,I)
            ZI=CRD(3,I)
            XT=CRD(1,I)-VEL(1,I)*DELTA-AVGX0(1,I)
            YT=CRD(2,I)-VEL(2,I)*DELTA-AVGX0(2,I)
            ZT=CRD(3,I)-VEL(3,I)*DELTA-AVGX0(3,I)
            X1I=AVGX1(1,I)+XT
            Y1I=AVGX1(2,I)+YT
            Z1I=AVGX1(3,I)+ZT
            X2I=AVGX2(1,I)+XT
            Y2I=AVGX2(2,I)+YT
            Z2I=AVGX2(3,I)+ZT
            FACT1=TSGFAC*AMASSI
            GX0=FACT1*(XI-X1I)
            GY0=FACT1*(YI-Y1I)
            GZ0=FACT1*(ZI-Z1I)
            X1I=SGAVG0*(X1I+DCOM1(1))+SGAVG1*XI
            Y1I=SGAVG0*(Y1I+DCOM1(2))+SGAVG1*YI
            Z1I=SGAVG0*(Z1I+DCOM1(3))+SGAVG1*ZI
            X2I=SGAVG0*(X2I+DCOM2(1))+SGAVG1*X1I
            Y2I=SGAVG0*(Y2I+DCOM2(2))+SGAVG1*Y1I
            Z2I=SGAVG0*(Z2I+DCOM2(3))+SGAVG1*Z1I
            AVGX0(1,I)=XI
            AVGX0(2,I)=YI
            AVGX0(3,I)=ZI
            AVGX1(1,I)=X1I
            AVGX1(2,I)=Y1I
            AVGX1(3,I)=Z1I
            AVGX2(1,I)=X2I
            AVGX2(2,I)=Y2I
            AVGX2(3,I)=Z2I
            ! avg(p)
            GXI=FACT1*(XI-X1I)
            GYI=FACT1*(YI-Y1I)
            GZI=FACT1*(ZI-Z1I)
            PXT=(GXI-SGAVG0*GX0)/SGAVG1
            PYT=(GYI-SGAVG0*GY0)/SGAVG1
            PZT=(GZI-SGAVG0*GZ0)/SGAVG1
            ! average force deviation 
            FXI=tsgfac*(pxt-2.0d0*gxi+fact1*(x1i-x2i))
            FYI=tsgfac*(pyt-2.0d0*gyi+fact1*(y1i-y2i))
            FZI=tsgfac*(pzt-2.0d0*gzi+fact1*(z1i-z2i))
            ! grid distribution
            GX=(XI-GXMIN)/SGSIZE+1.0d0
            GY=(YI-GYMIN)/SGSIZE+1.0d0
            GZ=(ZI-GZMIN)/SGSIZE+1.0d0
            X0=INT(GX)
            Y0=INT(GY)
            Z0=INT(GZ)
            IF(X0<1.OR.Y0<1.OR.Z0<1)THEN
              STOP 'Position outside of lower boundary'
            ENDIF
            X1=X0+1
            Y1=Y0+1
            Z1=Z0+1
            IF(X1>NGRIDX.OR.Y1>NGRIDY.OR.Z1>NGRIDZ)THEN
              !write(mdout,*)'I= ',i,xi,yi,zi,x1,y1,z1
              STOP 'Position outside of higher boundary'
            ENDIF
            A0=X1-GX
            B0=Y1-GY
            C0=Z1-GZ
            A1=1.0d0-A0
            B1=1.0d0-B0
            C1=1.0d0-C0
            ABC000=A0*B0*C0
            ABC100=A1*B0*C0
            ABC010=A0*B1*C0
            ABC001=A0*B0*C1
            ABC110=A1*B1*C0
            ABC011=A0*B1*C1
            ABC101=A1*B0*C1
            ABC111=A1*B1*C1
            I000=X0+NGRIDX*(Y0-1+NGRIDY*(Z0-1))
            I100=I000+1
            I010=I000+NGRIDX
            I001=I000+NGRIDXY
            I110=I010+1
            I101=I001+1
            I011=I001+NGRIDX
            I111=I011+1
            RHOM(I000)=RHOM(I000)+ABC111*AMASSI
            RHOM(I100)=RHOM(I100)+ABC011*AMASSI
            RHOM(I010)=RHOM(I010)+ABC101*AMASSI
            RHOM(I001)=RHOM(I001)+ABC110*AMASSI
            RHOM(I110)=RHOM(I110)+ABC001*AMASSI
            RHOM(I101)=RHOM(I101)+ABC010*AMASSI
            RHOM(I011)=RHOM(I011)+ABC100*AMASSI
            RHOM(I111)=RHOM(I111)+ABC000*AMASSI
            RHOVX(I000)=RHOVX(I000)+ABC111*GXI
            RHOVY(I000)=RHOVY(I000)+ABC111*GYI
            RHOVZ(I000)=RHOVZ(I000)+ABC111*GZI
            RHOVX(I100)=RHOVX(I100)+ABC011*GXI
            RHOVY(I100)=RHOVY(I100)+ABC011*GYI
            RHOVZ(I100)=RHOVZ(I100)+ABC011*GZI
            RHOVX(I010)=RHOVX(I010)+ABC101*GXI
            RHOVY(I010)=RHOVY(I010)+ABC101*GYI
            RHOVZ(I010)=RHOVZ(I010)+ABC101*GZI
            RHOVX(I001)=RHOVX(I001)+ABC110*GXI
            RHOVY(I001)=RHOVY(I001)+ABC110*GYI
            RHOVZ(I001)=RHOVZ(I001)+ABC110*GZI
            RHOVX(I110)=RHOVX(I110)+ABC001*GXI
            RHOVY(I110)=RHOVY(I110)+ABC001*GYI
            RHOVZ(I110)=RHOVZ(I110)+ABC001*GZI
            RHOVX(I101)=RHOVX(I101)+ABC010*GXI
            RHOVY(I101)=RHOVY(I101)+ABC010*GYI
            RHOVZ(I101)=RHOVZ(I101)+ABC010*GZI
            RHOVX(I011)=RHOVX(I011)+ABC100*GXI
            RHOVY(I011)=RHOVY(I011)+ABC100*GYI
            RHOVZ(I011)=RHOVZ(I011)+ABC100*GZI
            RHOVX(I111)=RHOVX(I111)+ABC000*GXI
            RHOVY(I111)=RHOVY(I111)+ABC000*GYI
            RHOVZ(I111)=RHOVZ(I111)+ABC000*GZI

            RHOAX(I000)=RHOAX(I000)+ABC111*FXI
            RHOAY(I000)=RHOAY(I000)+ABC111*FYI
            RHOAZ(I000)=RHOAZ(I000)+ABC111*FZI
            RHOAX(I100)=RHOAX(I100)+ABC011*FXI
            RHOAY(I100)=RHOAY(I100)+ABC011*FYI
            RHOAZ(I100)=RHOAZ(I100)+ABC011*FZI
            RHOAX(I010)=RHOAX(I010)+ABC101*FXI
            RHOAY(I010)=RHOAY(I010)+ABC101*FYI
            RHOAZ(I010)=RHOAZ(I010)+ABC101*FZI
            RHOAX(I001)=RHOAX(I001)+ABC110*FXI
            RHOAY(I001)=RHOAY(I001)+ABC110*FYI
            RHOAZ(I001)=RHOAZ(I001)+ABC110*FZI
            RHOAX(I110)=RHOAX(I110)+ABC001*FXI
            RHOAY(I110)=RHOAY(I110)+ABC001*FYI
            RHOAZ(I110)=RHOAZ(I110)+ABC001*FZI
            RHOAX(I101)=RHOAX(I101)+ABC010*FXI
            RHOAY(I101)=RHOAY(I101)+ABC010*FYI
            RHOAZ(I101)=RHOAZ(I101)+ABC010*FZI
            RHOAX(I011)=RHOAX(I011)+ABC100*FXI
            RHOAY(I011)=RHOAY(I011)+ABC100*FYI
            RHOAZ(I011)=RHOAZ(I011)+ABC100*FZI
            RHOAX(I111)=RHOAX(I111)+ABC000*FXI
            RHOAY(I111)=RHOAY(I111)+ABC000*FYI
            RHOAZ(I111)=RHOAZ(I111)+ABC000*FZI
    ENDIF  /* ISGSTA */
  enddo 
    !
#ifdef MPI
        if(numtasks > 1)then
!  combining all node results
!
          call mpi_allreduce(MPI_IN_PLACE,rhom,ngridxyz, &
          mpi_double_precision,mpi_sum,pmemd_comm,ierr)
          call mpi_allreduce(MPI_IN_PLACE,rhovx,ngridxyz, &
          mpi_double_precision,mpi_sum,pmemd_comm,ierr)
          call mpi_allreduce(MPI_IN_PLACE,rhovy,ngridxyz, &
          mpi_double_precision,mpi_sum,pmemd_comm,ierr)
          call mpi_allreduce(MPI_IN_PLACE,rhovz,ngridxyz, &
          mpi_double_precision,mpi_sum,pmemd_comm,ierr)
          call mpi_allreduce(MPI_IN_PLACE,rhoax,ngridxyz, &
          mpi_double_precision,mpi_sum,pmemd_comm,ierr)
          call mpi_allreduce(MPI_IN_PLACE,rhoay,ngridxyz, &
          mpi_double_precision,mpi_sum,pmemd_comm,ierr)
          call mpi_allreduce(MPI_IN_PLACE,rhoaz,ngridxyz, &
          mpi_double_precision,mpi_sum,pmemd_comm,ierr)
        endif
#endif
  ! remove net velocity and acceleration
    AMASSI=SUM(RHOM)
    GXI=SUM(RHOVX)/AMASSI
    GYI=SUM(RHOVY)/AMASSI
    GZI=SUM(RHOVZ)/AMASSI
    FXI=SUM(RHOAX)/AMASSI
    FYI=SUM(RHOAY)/AMASSI
    FZI=SUM(RHOAZ)/AMASSI
    DO I=1,NGRIDXYZ
      RHOMI=RHOM(I)
      IF(RHOMI>1.0d-6)THEN
        RHOVX(I)=RHOVX(I)/RHOMI
        RHOVY(I)=RHOVY(I)/RHOMI
        RHOVZ(I)=RHOVZ(I)/RHOMI
        RHOAX(I)=RHOAX(I)/RHOMI
        RHOAY(I)=RHOAY(I)/RHOMI
        RHOAZ(I)=RHOAZ(I)/RHOMI
      ENDIF
    ENDDO
    RHOVX=RHOVX-GXI
    RHOVY=RHOVY-GYI
    RHOVZ=RHOVZ-GZI
    RHOAX=RHOAX-FXI
    RHOAY=RHOAY-FYI
    RHOAZ=RHOAZ-FZI
    RETURN
  END SUBROUTINE MAPBUILD

  SUBROUTINE mapvalues(amassi,crdi,avgpi,avgfi)

  !-----------------------------------------------------------------------
  !     This routine to find map values at input position
  !     throuth linear interpolation
  !-----------------------------------------------------------------------
  use mdin_ctrl_dat_mod,only:sgsize

  implicit none
  double precision AMASSI,TMASS
  double precision crdi(3),avgpi(3),avgfi(3)
  INTEGER X0,Y0,Z0,X1,Y1,Z1
  INTEGER I000,I100,I010,I001,I110,I101,I011,I111
  double precision XI,YI,ZI,VXI,VYI,VZI,AXI,AYI,AZI
  double precision GX,GY,GZ,A0,B0,C0,A1,B1,C1
  double precision ABC000,ABC100,ABC010,ABC001,ABC110,ABC101,ABC011,ABC111
     GX=(crdi(1)-GXMIN)/SGSIZE+1.0d0
     GY=(crdi(2)-GYMIN)/SGSIZE+1.0d0
     GZ=(crdi(3)-GZMIN)/SGSIZE+1.0d0
     X0=INT(GX)
     Y0=INT(GY)
     Z0=INT(GZ)
     IF(X0<1.OR.Y0<1.OR.Z0<1)THEN
       STOP 'Position outside of lower boundary'
     ENDIF
     X1=X0+1
     Y1=Y0+1
     Z1=Z0+1
     IF(X1>NGRIDX.OR.Y1>NGRIDY.OR.Z1>NGRIDZ)THEN
       STOP 'Position outside of higher boundary'
     ENDIF
     A0=X1-GX
     B0=Y1-GY
     C0=Z1-GZ
     A1=1.0d0-A0
     B1=1.0d0-B0
     C1=1.0d0-C0
     ABC000=A0*B0*C0
     ABC100=A1*B0*C0
     ABC010=A0*B1*C0
     ABC001=A0*B0*C1
     ABC110=A1*B1*C0
     ABC011=A0*B1*C1
     ABC101=A1*B0*C1
     ABC111=A1*B1*C1
     I000=X0+NGRIDX*(Y0-1+NGRIDY*(Z0-1))
     I100=I000+1
     I010=I000+NGRIDX
     I001=I000+NGRIDXY
     I110=I010+1
     I101=I001+1
     I011=I001+NGRIDX
     I111=I011+1
     TMASS=ABC111*RHOM(I000)+ABC110*RHOM(I001)+ABC101*RHOM(I010)+ABC011*RHOM(I100)+     &
     ABC100*RHOM(I011)+ABC010*RHOM(I101)+ABC001*RHOM(I110)+ABC000*RHOM(I111)
     VXI=ABC111*RHOVX(I000)+ABC110*RHOVX(I001)+ABC101*RHOVX(I010)+ABC011*RHOVX(I100)+     &
     ABC100*RHOVX(I011)+ABC010*RHOVX(I101)+ABC001*RHOVX(I110)+ABC000*RHOVX(I111)
     VYI=ABC111*RHOVY(I000)+ABC110*RHOVY(I001)+ABC101*RHOVY(I010)+ABC011*RHOVY(I100)+     &
     ABC100*RHOVY(I011)+ABC010*RHOVY(I101)+ABC001*RHOVY(I110)+ABC000*RHOVY(I111)
     VZI=ABC111*RHOVZ(I000)+ABC110*RHOVZ(I001)+ABC101*RHOVZ(I010)+ABC011*RHOVZ(I100)+     &
     ABC100*RHOVZ(I011)+ABC010*RHOVZ(I101)+ABC001*RHOVZ(I110)+ABC000*RHOVZ(I111)
     AXI=ABC111*RHOAX(I000)+ABC110*RHOAX(I001)+ABC101*RHOAX(I010)+ABC011*RHOAX(I100)+     &
     ABC100*RHOAX(I011)+ABC010*RHOAX(I101)+ABC001*RHOAX(I110)+ABC000*RHOAX(I111)
     AYI=ABC111*RHOAY(I000)+ABC110*RHOAY(I001)+ABC101*RHOAY(I010)+ABC011*RHOAY(I100)+     &
     ABC100*RHOAY(I011)+ABC010*RHOAY(I101)+ABC001*RHOAY(I110)+ABC000*RHOAY(I111)
     AZI=ABC111*RHOAZ(I000)+ABC110*RHOAZ(I001)+ABC101*RHOAZ(I010)+ABC011*RHOAZ(I100)+     &
     ABC100*RHOAZ(I011)+ABC010*RHOAZ(I101)+ABC001*RHOAZ(I110)+ABC000*RHOAZ(I111)
     avgpi(1)=AMASSI*VXI
     avgpi(2)=AMASSI*VYI
     avgpi(3)=AMASSI*VZI
     avgfi(1)=AMASSI*AXI
     avgfi(2)=AMASSI*AYI
     avgfi(3)=AMASSI*AZI
     RETURN
  END SUBROUTINE MAPVALUES
  

#ifdef MPI
  SUBROUTINE MAPINIT(my_atm_cnt,CRD, my_atm_lst)
#else
  SUBROUTINE MAPINIT(atm_cnt,CRD)
#endif
  !-----------------------------------------------------------------------
  !     This routine build and initialize maps
  !
  !-----------------------------------------------------------------------
  use mdin_ctrl_dat_mod,only:sgsize

  implicit none
  integer atm_cnt
  double precision crd(3,*)

  integer I,IA,alloc_err,dealloc_err
  double precision XI,YI,ZI
  double precision XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN,TEMP1(6),TEMP2(6)
  logical rmap
#ifdef MPI
      integer               :: my_atm_lst(*)
      integer my_atm_cnt,atm_lst_idx,ierr
#endif

  XMIN=1.0D8
  XMAX=-1.0D8
  YMIN=1.0D8
  YMAX=-1.0D8
  ZMIN=1.0D8
  ZMAX=-1.0D8
#ifdef MPI
        do atm_lst_idx = 1, my_atm_cnt
        i = my_atm_lst(atm_lst_idx)
        !if (gbl_atm_owner_map(i) .ne. mytaskid) cycle
#else
        do i = 1, atm_cnt
#endif

       XI=CRD(1,I)
       YI=CRD(2,I)
       ZI=CRD(3,I)
       IF(XI>XMAX)XMAX=XI
       IF(XI<XMIN)XMIN=XI
       IF(YI>YMAX)YMAX=YI
       IF(YI<YMIN)YMIN=YI
       IF(ZI>ZMAX)ZMAX=ZI
       IF(ZI<ZMIN)ZMIN=ZI
  ENDDO
#ifdef MPI
        if(numtasks > 1)then
!  combining all node results
!
          temp1(1)=XMAX
          temp1(2)=-XMIN
          temp1(3)=YMAX
          temp1(4)=-YMIN
          temp1(5)=ZMAX
          temp1(6)=-ZMIN

          call mpi_allreduce(temp1,temp2,6, &
          mpi_double_precision,mpi_max,pmemd_comm,ierr)

          xmax=temp2(1)
          xmin=-temp2(2)
          ymax=temp2(3)
          ymin=-temp2(4)
          zmax=temp2(5)
          zmin=-temp2(6)
        endif
#endif
  RMAP=(XMAX>GXMAX.OR.XMIN<GXMIN.OR.YMAX>GYMAX.OR.YMIN<GYMIN.OR.ZMAX>GZMAX.OR.ZMIN<GZMIN)
  IF(RMAP)THEN

    GXMAX=SGSIZE*(AINT(XMAX/SGSIZE)+1.0d0)
    GXMIN=SGSIZE*(AINT(XMIN/SGSIZE)-1.0d0)
    GYMAX=SGSIZE*(AINT(YMAX/SGSIZE)+1.0d0)
    GYMIN=SGSIZE*(AINT(YMIN/SGSIZE)-1.0d0)
    GZMAX=SGSIZE*(AINT(ZMAX/SGSIZE)+1.0d0)
    GZMIN=SGSIZE*(AINT(ZMIN/SGSIZE)-1.0d0)
    NGRIDX=INT((GXMAX-GXMIN)/SGSIZE)+1
    NGRIDY=INT((GYMAX-GYMIN)/SGSIZE)+1
    NGRIDZ=INT((GZMAX-GZMIN)/SGSIZE)+1
    NGRIDXY=NGRIDX*NGRIDY
    NGRIDYZ=NGRIDY*NGRIDZ
    NGRIDZX=NGRIDZ*NGRIDX
    NGRIDXYZ=NGRIDX*NGRIDY*NGRIDZ
    ! deallocate maps
    if(allocated(rhom))deallocate(rhom)  ! deallocate sgmaps
    if(allocated(rhovx))deallocate(rhovx)  ! deallocate sgmaps
    if(allocated(rhovy))deallocate(rhovy)  ! deallocate sgmaps
    if(allocated(rhovz))deallocate(rhovz)  ! deallocate sgmaps
    if(allocated(rhoax))deallocate(rhoax)  ! deallocate sgmaps
    if(allocated(rhoay))deallocate(rhoay)  ! deallocate sgmaps
    if(allocated(rhoaz))deallocate(rhoaz)  ! deallocate sgmaps
    !if(allocated(rhoaz))deallocate(rhoaz,stat=dealloc_err)  ! deallocate sgmaps
    !if(dealloc_err /= 0 ) then
    !  stop "unable to deallocate SG maps "
    !endif
    allocate( rhom(ngridxyz),rhovx(ngridxyz),rhovy(ngridxyz),rhovz(ngridxyz),rhoax(ngridxyz),rhoay(ngridxyz),rhoaz(ngridxyz),stat=alloc_err)  ! allocate sgmaps
    !if (alloc_err .ne. 0) then
    !  stop "unable to allocate SG maps "
    !endif
  ENDIF
    RHOM=0.0d0
    RHOVX=0.0d0
    RHOVY=0.0d0
    RHOVZ=0.0d0
    RHOAX=0.0d0
    RHOAY=0.0d0
    RHOAZ=0.0d0
    RETURN
  END SUBROUTINE MAPINIT
  
end module sgld_mod

