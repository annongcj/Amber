!****************************************************************************
! Continuous Constant pH Molecular Dynamics
! Both conformational and protonation state sampling with GBneck2
! Robert C. Harris, Yandong Huang, and Jana Shen
! University of Maryland School of Pharmacy, Feb 2018 
!******************************************************************************
module phmd_mod

use random_mod

use gbl_constants_mod, only: pi,kb,AMBER_ELECTROSTATIC
use prmtop_dat_mod, only : natom,ntypes,nres

! Dimension setting up, like MAX_TITR_STATES
! Limits on the sizes of the data structures
#define MAX_TITR_STATES 512
#define MAX_TYPES 128
#define MAX_CHARGES 128 

   ! double parameters
   double precision, parameter :: zero      = 0.0d0
   double precision, parameter :: one       = 1.0d0
   double precision, parameter :: two       = 2.0d0
   double precision, parameter :: three     = 3.0d0
   double precision, parameter :: four      = 4.0d0
   double precision, parameter :: five      = 5.0d0
   double precision, parameter :: seven     = 7.0d0
   double precision, parameter :: ten       = 10.0d0

   double precision, parameter :: half      = 0.5d0
   double precision, parameter :: fourth    = 0.25d0

   ! PHMD energy and force parameters
   double precision,save,allocatable,dimension(:) :: dph_theta, dph_avg, &
        para, parb, park, local_model_pka
   double precision,save,allocatable,dimension(:,:) :: qstate1, qstate2, &
        vstate1, vstate2, parmod, sp_par

   ! Separate charge from MD
   double precision,save,allocatable,dimension(:,:) :: qstate1_md, qstate2_md
   double precision,save,allocatable,dimension(:) :: charge_phmd

   double precision, save :: vph_theta(MAX_TITR_STATES), &
                             ph_theta(MAX_TITR_STATES), &
                             thetaold(MAX_TITR_STATES), &
                             dphold(MAX_TITR_STATES), &
                             vphold(MAX_TITR_STATES), &
                             barr(MAX_TITR_STATES)

    integer, save          :: link_res_list(MAX_TITR_STATES)

   ! PHMD topology arrays
   integer,allocatable,dimension(:) :: sp_grp, grplist, titrres, resboundsa, &
                                       resboundsb, titratoms, &
                                       linkgrplist
   logical,allocatable,dimension(:) :: baseres

   ! PHMD I/O parameters
   integer :: ntitr,nprint_phmd,phmd_igb, nlinktitr, &
        nsolute,iphfrq,max_numch, ntitratoms

   integer, save :: iupdate


  logical prlam,prderiv,qsetlam,qphmdstart,qbarread, add_energies, prke, &
          use_linked, use_forces, prdudls

   ! Thermostat parameters in PHMD 
   double precision, save :: qmass_phmd, temp_phmd, &
        kaytee, &
        phbeta, &
        phmdcut

   double precision, save :: phkin

   double precision, save :: timfac  
   double precision, save :: kbt  ! Boltzmann constant times bath temperature
   double precision, save :: gam  ! the product of gamma, timestep, and timfac

   integer outu  ! Unit number for all standard CHARMM output

   ! Input parameters original in LoadPHMD function
   character(len=6) :: res_name(MAX_TYPES)
   character(len=4) :: atom_name(MAX_TYPES,MAX_CHARGES)
   character(len=10) :: fmt00
   character(len=20) :: fmt120
   double precision, save :: ch(MAX_TYPES,MAX_CHARGES)
   double precision, save :: ch_md(MAX_TYPES,MAX_CHARGES)
   double precision, save :: rad(MAX_TYPES,MAX_CHARGES)
   double precision, save :: model_pka(MAX_TYPES,2),parameters(MAX_TYPES,12), &
            bar(MAX_TYPES,2)       
   
   integer prnlev ! print lambda values for the last step and print level.
   ! NCT: number of titratable residue type,
   ! NUMCH: number of atoms in one residue
   integer ngt, res_type(MAX_TYPES),numch(MAX_TYPES)
   character(len=4) :: masktitrres(MAX_TYPES) ! phmd selection by residue name
   integer masktitrrestypes ! Number of residue type for phmd

   ! Input for restraint cphmd by DoPHMDTest
   integer phtest

contains

subroutine phmd_zero()
  !***********************************
  ! Initialize arrays and variables.
  !***********************************
   use file_io_dat_mod, only : mdout

   timfac = 4.88882E-002
   qmass_phmd = ten
   temp_phmd = 300
   phbeta = five
   iphfrq = 1
   nprint_phmd = 500
   prlam = .false.
   prderiv = .false.
   prdudls = .false.
   use_forces = .true.
   res_name(:) = ''
   atom_name(:,:) = ''
   ch(:,:) = 0
   ch_md(:,:) = 0
   rad(:,:) = 0
   model_pka(:,:) = 0
   parameters(:,:) = 0
   bar(:,:) = 0
   outu = mdout 
   prnlev = 6
   iupdate = 0
   ngt = 0
   numch(:) = zero
   phtest = 0
   masktitrres(:) = '' 
   masktitrrestypes = zero
   ph_theta(:) = pi*fourth
   vph_theta(:) = zero
   phmdcut = 1000.0
   phmd_igb = 8
   qphmdstart = .true.
   qbarread = .false.
   add_energies = .true.
   prke = .false.
   barr(:) = zero
   use_linked = .false.

   return

end subroutine phmd_zero

SUBROUTINE phmdzero()
  !***********************************
  ! Start PHMD dynamics
  ! If a test is being run, the
  ! velocities of the lambda particles
  ! are set to zero.
  ! Otherwise, if this is a new run
  ! (QPHMDStart == true) set the
  ! initial velocities to a 
  ! Gaussian distribution.
  ! Otherwise, read velocities from
  ! phmd restart file.
  !***********************************

   use random_mod, only: gauss
#ifdef MPI
   use parallel_dat_mod, only : master
#endif
   use mdin_ctrl_dat_mod, only : iphmd, solvph, temp0
   use remd_mod, only : remd_method
   use prmtop_dat_mod, only : atm_numex, gbl_natex
   use gb_ene_mod,       only : excludeatoms

   implicit none

   integer i
   double precision r1, r2

#ifdef MPI
   if (master) then
#endif
     if (remd_method .eq. 1) then
       kaytee  = half * kb * temp0
     else
       kaytee  = half * kb * temp_phmd
     end if
     r1 = dsqrt(two * kaytee / qmass_phmd)
     phkin = zero
     if (qsetlam) then
        do i = 1, ntitr
           vph_theta(i) = zero
        end do
     else if (qphmdstart)  then
        do i=1,ntitr
           call gauss(0.d0,r1,r2)
           vphold(i) = r2
        end do
     else
        do I=1, ntitr
           vphold(i) = vph_theta(i)
        end do
     endif
#ifdef MPI
   end if
#endif

#ifdef CUDA
   call gpu_upload_phmd(qsetlam, ntitr, qmass_phmd, temp_phmd, &
   phbeta, para, parb, park, barr, ph_theta, dph_theta, vph_theta, &
   charge_phmd, grplist, sp_grp, sp_par, qstate1, qstate2, vstate1, &
   vstate2, qstate1_md, qstate2_md, parmod, resboundsa, resboundsb, baseres, &
   max_numch, excludeatoms, iphmd, solvph, &
   titratoms, ntitratoms, atm_numex, gbl_natex, linkgrplist, use_linked, &
   use_forces);
#endif
   return

end subroutine phmdzero

subroutine startphmd(charge,lbres,igraph,ipres)
   !***********************************************************
   ! read input parameters and set up parameter, charge, and
   ! vdw arrays.
   !***********************************************************

   ! As in constantph module, the following modules are used.
   use file_io_dat_mod, only : phmd_unit, phmdout_unit, phmdin_name, owrite
   use mdin_ctrl_dat_mod, only : iphmd, solvph, temp0, gb_kappa
   use parallel_dat_mod, only : master
   use file_io_mod,     only : amopen
   use prmtop_dat_mod,   only : atm_igraph, atm_isymbl
   use gb_ene_mod,       only : excludeatoms

   implicit none
   integer i, j, k, ier 

   ! Variables defined in LoadPHMD()
   character(len=10) fmt00
   character(len=20) fmt120
   integer ires, numfound
   double precision x1, x2, qunprot, qprot

   logical match
   logical iftitrres

   double precision, intent(inout) :: charge(*)
   character(len=4), intent(in) :: lbres(nres) ! residue name,(lbres(i),i=1,nres)
   integer, intent(in) :: ipres(*) ! pointer of residue, (ipres(i),i=1,nres), (j=ipres(i),ipres(i+1)-1)
   character(len=4), intent(in) :: igraph(*) ! atom type, (isymbl(i),i=1,natom)

   ! Read in namelists from three files, phmdin, phmdstrt, and phmdparm
   ! Control input structure phmdin
   namelist /phmdin/ nsolute, phmdcut, &
                 qmass_phmd, temp_phmd, &
                 phbeta, iphfrq, prdudls, &
                 qphmdstart, qbarread, add_energies, & 
                 nprint_phmd, prlam, prderiv, & 
                 prnlev, outu, link_res_list, &
                 phtest, barr, prke, &
                 masktitrres, masktitrrestypes, &
                 use_forces

   ! Parameter input structure phmdparm
   namelist /phmdparm/ ngt, numch, res_name, res_type, atom_name, ch, ch_md, &
                rad, model_pka, parameters, bar

   ! Restart input structure phmdstrt
   namelist /phmdstrt/ ph_theta, vph_theta ! status of phmd

   call amopen(phmd_unit, phmdin_name,'O', 'F', 'R')
   read(phmd_unit, nml=phmdin)

   call amopen(phmdparm_unit, phmdparm_name, 'O', 'F', 'R')
   read(phmdparm_unit, nml=phmdparm)

   if (.not. qphmdstart) then
      call amopen(phmdstrt_unit, phmdstrt_name, 'O', 'F', 'R')
      read(phmdstrt_unit, nml=phmdstrt)
   endif
   call amopen(PHMDOUT_UNIT,phmdout_name,owrite, 'F', 'w')
   iupdate = 0

   qsetlam = .FALSE.

   allocate(qstate1(2,natom), qstate2(2,natom), vstate1(2,natom), & 
            vstate2(2,natom), qstate1_md(2,natom), qstate2_md(2,natom), &
            dph_theta(MAX_TITR_STATES), dph_avg(MAX_TITR_STATES), &
            charge_phmd(natom), para(MAX_TITR_STATES), parb(MAX_TITR_STATES), &
            parmod(6,MAX_TITR_STATES), park(MAX_TITR_STATES), &
            local_model_pka(MAX_TITR_STATES), sp_par(2,MAX_TITR_STATES), &
            sp_grp(MAX_TITR_STATES), &
            titrres(MAX_TITR_STATES*2), grplist(natom), linkgrplist(natom), &
            resboundsa(MAX_TITR_STATES), resboundsb(MAX_TITR_STATES), &
            baseres(MAX_TITR_STATES), excludeatoms(natom), titratoms(natom), &
            stat=ier)


   ! Step 2: Assign QState, Par to atoms
   !************************************
   ! First set QState1/2 <- Charge
   ! 

   VState1(1:2,1:natom) = zero  ! vdw state: 1 for protonated H
   VState2(1:2,1:natom) = zero

   ntitratoms = 0
   ntitr = 0
   nlinktitr = 0
   do i = 1, natom
      grplist(i) = 0         ! titr group number
      linkgrplist(i) = 0
      qstate1(1,i) = charge(i)
      qstate1(2,i) = charge(i)
      qstate2(1,i) = charge(i)
      qstate2(2,i) = charge(i)
      qstate1_MD(1,i) = charge(i)
      qstate1_MD(2,i) = charge(i)
      qstate2_MD(1,i) = charge(i)
      qstate2_MD(2,i) = charge(i)

      charge_phmd(i) = charge(i)
      if ( ( atm_igraph(i)(1:2) .eq. 'HD' &
          .or. atm_igraph(i)(1:2) .eq. 'HE') &
          .and. atm_isymbl(i) .eq. 'HO') then
         excludeatoms(i) = 1
      else
         excludeatoms(i) = 0
      end if 
   end do

   max_numch = -1
   do ires = 1, nres ! input residues
      do j = 1, ngt  ! number of titratable residue types

         match=(lbres(ires).eq. res_name(j))

         ! Check if the residue is selected to be titrated.
         IfTitrRes = .false.
         do k = 1,masktitrrestypes
            if(masktitrres(k) .eq. lbres(ires)) then
               iftitrres = .true.
            end if
         end do

      
         if (match .and. iftitrres) then
            ! Check if ATOM_TYPE in phmd param file agree with IGRAPH.
            if (numch(j) .ge. max_numch) then
               max_numch = numch(j)
            end if
            numfound = 0
            do i = ipres(ires), ipres(ires+1)-1
               do k = 1,numch(j)
                  if (atom_name(j,k).eq.igraph(i)) then
                     numfound = numfound + 1 ! atom match
                  end if
               end do
            end do
            if (numfound .lt. numch(j)) THEN
               if (prnlev .gt. 5)  &
                 write(outu,'(a,1X,a,1X,a,1X,I4,1X,I4)') &
                 'insufficient atom match', &
                 lbres(ires), res_name(j), numfound, numch(j)
            else  ! all atoms match
               ! Single site
               if (res_type(j) .eq. 0) then
                  ntitr = ntitr + 1
                  resboundsa(ntitr) = ipres(ires)
                  resboundsb(ntitr) = ipres(ires+1)-1
                  baseres(ntitr) = .true.
                  titrres(ntitr) = ires
                  sp_grp(ntitr) = 0
                  sp_par(1,ntitr) = zero
                  sp_par(2,NTitr) = zero
                  para(ntitr) = parameters(j,1)
                  parb(ntitr) = parameters(j,2)
                  if(.not. qbarread) then
                    barr(ntitr) = bar(j,1)
                  end if
                  park(ntitr) = log(10.0d0) * kb * temp_phmd * &
                                (model_pka(j,1) - solvph)
                  local_model_pka(ntitr) = model_pka(j,1)
                  dph_theta(ntitr) = zero
                  dph_avg(ntitr) = zero
                  if (qphmdstart) then
                     if (solvph .ge. model_pka(j,1)) then
                        ph_theta(ntitr) = half * pi
                     else
                        ph_theta(ntitr) = 0
                     end if
                  end if
                  do i = ipres(ires), ipres(ires+1)-1
                     do k = 1,numch(j)
                        if (atom_name(j,k).eq.igraph(i)) then
                           ntitratoms = ntitratoms + 1
                           titratoms(ntitratoms) = i - 1
                           qstate1(1,i) = ch(j,k)*AMBER_ELECTROSTATIC
                           qstate2(1,i) = ch(j,k+numch(j))*AMBER_ELECTROSTATIC

                           qstate1(2,i) = qstate1(1,i)
                           qstate2(2,i) = qstate2(1,i)

                           qstate1_md(1,i) = ch_md(j,k)*AMBER_ELECTROSTATIC
                           qstate2_md(1,i) = ch_md(j,k+numch(j)) * &
                                             AMBER_ELECTROSTATIC

                           qstate1_md(2,i) = qstate1_md(1,i)
                           qstate2_md(2,i) = qstate2_md(1,i)

                           vstate1(1,i) = rad(j,k)
                           vstate2(1,i) = rad(j,k+numch(j))

                           grplist(i) = ntitr
                        end if
                     end do
                  end do
               ! histidine
               else if (res_type(j) .eq. 2) then
                  ntitr = ntitr + 2
                  titrres(ntitr) = ires
                  titrres(ntitr-1) = ires
                  resboundsa(ntitr-1) = ipres(ires)
                  resboundsb(ntitr-1) = ipres(ires+1)-1
                  baseres(ntitr-1) = .true.
                  baseres(ntitr) = .false.
                  sp_grp(ntitr) = 2
                  sp_grp(ntitr-1) = 1
                  para(ntitr) = parameters(j,3)
                  parb(ntitr) = parameters(j,4)
                  para(ntitr-1) = parameters(j,1)
                  parb(ntitr-1) = parameters(j,2)
                  sp_par(1,ntitr) = parameters(j,5)
                  sp_par(2,ntitr) = parameters(j,6)
                  if(.not. qbarread) then
                    barr(ntitr) = bar(j,2)
                    barr(ntitr-1) = bar(j,1)
                  end if
                  park(ntitr) = log(10.0d0) * kb * temp_phmd * &
                                (model_pka(j,2) - solvph)
                  local_model_pka(ntitr) = model_pka(j,2)
                  park(ntitr-1) = log(10.0d0) * kb * temp_phmd * &
                                  (model_pka(j,1) - solvph)
                  local_model_pka(ntitr-1) = model_pka(j,1)
                  dph_theta(ntitr) = zero
                  dph_avg(ntitr) = zero
                  dph_theta(ntitr-1) = zero
                  dph_avg(ntitr-1) = zero
                  if (qphmdstart) then
                     if (solvph .ge. model_pka(j,1)) then
                        ph_theta(ntitr-1) = half * pi
                        ph_theta(ntitr) = half * pi
                     else
                        ph_theta(ntitr-1) = 0
                        ph_theta(ntitr) = half * pi
                     end if
                  end if
                  do i = ipres(ires), ipres(ires+1)-1
                     do k = 1,numch(j)
                        if (atom_name(j,k).eq.igraph(i)) then
                            ntitratoms = ntitratoms + 1
                            titratoms(ntitratoms) = i - 1
                            qstate1(1,i) = ch(j,k)*AMBER_ELECTROSTATIC
                            qstate2(1,i) = ch(j,k+numch(j))*AMBER_ELECTROSTATIC
  
                            qstate1(2,i) = ch(j,k)*AMBER_ELECTROSTATIC
                            qstate2(2,i) = ch(j,k+2*numch(j)) * &
                                           AMBER_ELECTROSTATIC
  
                            qstate1_md(1,i) = ch_md(j,k)*AMBER_ELECTROSTATIC
                            qstate2_md(1,I) = ch_md(j,k+numch(j)) * &
                                              AMBER_ELECTROSTATIC
  
                            qstate1_md(2,i) = ch_md(j,k)*AMBER_ELECTROSTATIC
                            qstate2_MD(2,i) = ch_md(j,k+2*numch(j)) * &
                                              AMBER_ELECTROSTATIC

                            vstate1(1,i) = rad(j,k)
                            vstate2(1,i) = rad(j,k+2*numch(j))

                            vstate1(2,i) = rad(j,k+numch(j))
                            vstate2(2,i) = rad(j,k+2*numch(j))

                            grplist(i) = ntitr-1
                        end if
                     end do
                  end do
               !glu/asp
               else if (res_type(j) .eq. 4) then 
                  ntitr = ntitr + 2
                  titrres(ntitr) = ires
                  titrres(ntitr-1) = ires
                  resboundsa(ntitr-1) = ipres(ires)
                  resboundsb(ntitr-1) = ipres(ires+1)-1
                  baseres(ntitr-1) = .true.
                  baseres(ntitr) = .false.
                  sp_grp(ntitr) = 4
                  sp_grp(ntitr-1) = 3
                  do i=1,6
                     parmod(i,ntitr) = parameters(j,i)
                  end do
                  if(.not. qbarread) then
                    barr(ntitr) = bar(j,2)
                    barr(ntitr-1) = bar(j,1)
                  end if
                  park(ntitr) = log(10.0d0) * kb * temp_phmd * &
                                (model_pka(j,1) - solvph)
                  local_model_pka(ntitr) = model_pka(j,1)
                  park(ntitr-1) = log(10.0d0) * kb * temp_phmd * &
                                  (model_pka(j,1) - solvph)
                  local_model_pka(ntitr-1) = model_pka(j,1)
                  dph_theta(ntitr) = zero
                  dph_avg(ntitr) = zero
                  dph_theta(ntitr-1) = zero
                  dph_avg(ntitr-1) = zero
                  if (qphmdstart) then
                     if (solvph .ge. model_pka(j,1)) then
                        ph_theta(ntitr-1) = half * pi
                        ph_theta(ntitr) = half * pi
                     else
                        ph_theta(ntitr-1) = 0
                        ph_theta(ntitr) = half * pi
                     end if
                  end if
                  do i = ipres(ires), ipres(ires+1)-1
                     do k = 1,numch(j)
                        if (atom_name(j,k).eq.igraph(i)) then
                           ntitratoms = ntitratoms + 1
                           titratoms(ntitratoms) = i - 1
                           qstate1(1,i) = ch(j,k)*AMBER_ELECTROSTATIC
                           qstate2(1,i) = ch(j,k+2*numch(j)) * &
                                          AMBER_ELECTROSTATIC
  
                           qstate1(2,i) = ch(j,k+numch(j))*AMBER_ELECTROSTATIC
                           qstate2(2,i) = ch(j,k+2*numch(j)) * &
                                          AMBER_ELECTROSTATIC
                     
                           qstate1_md(1,I) = ch_md(j,k)*AMBER_ELECTROSTATIC
                           qstate2_md(1,i) = ch_md(j,k+2*numch(j)) * &
                                             AMBER_ELECTROSTATIC

                           qstate1_md(2,i) = ch_md(j,k+numch(j)) * &
                                             AMBER_ELECTROSTATIC
                           qstate2_md(2,i) = ch_md(j,k+2*numch(j)) * &
                                             AMBER_ELECTROSTATIC

                           vstate1(1,i) = rad(j,k)
                           vstate2(1,i) = rad(j,k+2*numch(j))

                           vstate1(2,i) = rad(j,k+numch(j))
                           vstate2(2,i) = rad(j,k+2*numch(j))

                           grplist(i) = ntitr-1
                        end if
                     end do
                  end do
               else if (res_type(j) .eq. -2) then
                  nlinktitr = nlinktitr + 1
                  do i = ipres(ires), ipres(ires+1)-1
                     do k = 1,numch(j)
                        if (atom_name(j,k) .eq. igraph(i)) then
                           qstate1(1,i) = ch(j,k)*AMBER_ELECTROSTATIC
                           qstate2(1,i) = ch(j,k+numch(j))*AMBER_ELECTROSTATIC

                           qstate1(2,i) = qstate1(1,i)
                           qstate2(2,i) = qstate2(1,i)

                           qstate1_md(1,i) = ch_md(j,k)*AMBER_ELECTROSTATIC
                           qstate2_md(1,i) = ch_md(j,k+numch(j)) * &
                                             AMBER_ELECTROSTATIC

                           qstate1_md(2,i) = qstate1_md(1,i)
                           qstate2_md(2,i) = qstate2_md(1,i)

                           vstate1(1,i) = rad(j,k)
                           vstate2(1,i) = rad(j,k+numch(j))
                           linkgrplist(i) = grplist(ipres(link_res_list(nlinktitr)))
                           use_linked = .true.
                        end if
                     end do
                  end do
                write(outu, '(a,i4,a,i4,a,i4)') ' PHMD> residue ', ires, &
                      ' linked to residue ', link_res_list(nlinktitr), &
                      ' titr ', linkgrplist(ipres(ires))
               end if ! res_type switch
            end if               ! numfound > numch
         end if                  ! match
      end do                       ! j
   end do                       ! ires

   if (prnlev .ge. 5) then    
      write(outu,'(a,F7.2)') ' PHMD> simulation pH = ',solvph
      write(outu,'(a,I4)')    ' PHMD> titr grps     = ',ntitr
   end if

   if ( phtest .eq. 1 ) then
      qsetlam = .true.
   end if

   !Assign initial charges for titratable groups
   do I=1,natom
      j = grplist(i)
      if (j .ne. 0) then
         x1 = dsin(ph_theta(j))**2
         if (sp_grp(j) .le. 0) then
            x2 = one
         else
            x2 = dsin(ph_theta(j+1))**2
         end if
         qprot = x2*qstate1(1,i) + (one-x2)*qstate1(2,i)
         qunprot = x2*qstate2(1,i) + (one-x2)*qstate2(2,i)
         charge_phmd(i) = (one-x1)* qprot + x1 * qunprot

         ! Initial charge for MD
         qprot = x2*qstate1_md(1,i) + (one-x2)*qstate1_md(2,i)
         qunprot = x2*qstate2_md(1,i) + (one-x2)*qstate2_md(2,i)
         charge(i) = (one-x1)* qprot + x1 * qunprot
      end if
      j = linkgrplist(i)
      ! first
      if (j .ne. 0) then
         x1 = dsin(ph_theta(j))**2
         charge_phmd(i) = x1 * qstate1(1, i) + (one - x1) * qstate2(1,i)
         charge(i) = x1 * qstate1_md(1, i) + (one - x1) * qstate2_md(1,i)
      end if
   end do

   if (iphmd .eq. 3) then
     write(outu,'(a)') 'Doing PME CpHMD. Titration charges set equal to MD charges.'
     do i=1,natom
       charge_phmd(i) = charge(i)
       qstate1(1,i) = qstate1_md(1,i)
       qstate1(2,i) = qstate1_md(2,i)
       qstate2(1,i) = qstate2_md(1,i)
       qstate2(2,i) = qstate2_md(2,i)
     end do
   end if

   do i=1,ntitr
      if (prnlev .gt. 5) write(outu,'(a,i4,2f6.2)')  &
          ' PHMD> lambda,barrier= ',i, dsin(ph_theta(i))**2,barr(i)
   end do


   if (prlam .or. prderiv .or. prke .and. prnlev .ge. 0 ) then
      ! Print RES,TITR,PKA INFO
      write(fmt00,'("(a,",i4,"i5)")') ntitr
      write(fmt120,'("(a, ",i4,"(f8.3))")') ntitr
      write(phmdout_unit,fmt00) '# ititr ',(i,i=1,ntitr)
      write(phmdout_unit,fmt00) '#  ires ',(titrres(i),i=1,ntitr)
      write(phmdout_unit,fmt00) '# itauto',(sp_grp(i),i=1,ntitr)
      write(phmdout_unit,fmt120) '# ParK',(park(i),i=1,ntitr)
   end if

 
   return

end subroutine startphmd

subroutine updatephmd(dtx, charge, istp, kin_ene)
   !***********************************************************
   ! Use Langevin integrator to update lambdas and velocities,
   ! update charges, and output to .lambda file
   !***********************************************************

   use mdin_ctrl_dat_mod, only : temp0
   use file_io_dat_mod, only : phmdout_unit
   use remd_mod, only : remd_method
#ifdef MPI
   use parallel_dat_mod, only : mpi_double_precision, pmemd_comm, master, mpi_sum
#endif

   implicit none

   ! Variables
   !
   ! Passed:
   !  dtx            : time step with the unit of AKMA
   !  charge         : atomic charge
   !  istp           : total time
   !  kin_ene        : kinetic energy
   ! Internal:
   !  full           : do we write a full record or not
   !  i,j            : loop counter

   double precision, intent(in)    :: dtx
   double precision, intent(inout) :: charge(*)
   double precision, intent(inout) :: kin_ene
   integer, intent(in)   :: istp
   double precision temp_dph_theta(MAX_TITR_STATES)
   double precision dudls(MAX_TITR_STATES)
   integer :: i,j,ierr

   double precision x1,x2,qunprot,qprot
   double precision rfd, frand ! variables for Langevin dynamics

   double precision thetanew(ntitr)
   character(len=10) fmt00
   character(len=20) fmt120
   character(len=20) fmt100, fmt110
   character(len=25) fmt115
#ifdef CUDA
   if ((.not. qsetlam) .and. (istp .ne. 0) ) then
      iupdate = iupdate + 1
      if(mod(iupdate,iphfrq) == 0) then
         call gpu_update_phmd(dtx)
      end if
   end if
   if(mod(istp,nprint_phmd) .eq. 0) then
      call gpu_download_phmd_theta(ph_theta, dph_theta, vph_theta, dudls)
   end if
   if (qsetlam) then
     call convphmd()
   end if
#else
   temp_dph_theta = 0
#ifdef MPI
   call mpi_allreduce(dpH_Theta, temp_dph_theta, MAX_TITR_STATES, &
                      mpi_double_precision, mpi_sum, pmemd_comm, ierr)
   dpH_Theta(:) = temp_dph_theta(:)
   if( master ) then
#endif
   call CONVPHMD()


   if ((.not. qsetlam) .and. (istp .ne. 0) ) then
      iupdate = iupdate + 1
      if(mod(iupdate,iphfrq) == 0)then
         ! KBT - Boltzmann constant times bath temperature.
         ! GAM - The product of gamma, timestep, and timfac.
         if (remd_method .eq. 1) then
           kbt=kb*temp0
         else
           kbt=kb*temp_phmd
         end if
         gam=timfac*phbeta*dtx
         rfd=sqrt(two*qmass_phmd*gam*kbt)/dtx
         do i=1,ntitr
            call gauss(0d0, rfd, frand)
            dph_theta(i)=dph_theta(i)+frand
            vph_theta(i)=vphold(i)-(dtx/two)*(dph_theta(i)/qmass_phmd) &
                  -(dtx/two)*(dphold(i)/qmass_phmd)

            vph_theta(i)=(one-gam)*vph_theta(i)
            thetanew(i) = ph_theta(i)+vph_theta(i)*dtx &
                  -dph_theta(i)*dtx*dtx/(two*qmass_phmd)
            dphold(i)=dph_theta(i)
            vphold(i)=vph_theta(i)
         end do

         do i=1,ntitr
            thetaold(i) = ph_theta(i)
            ph_theta(i) = thetanew(i)
         end do

      end if ! iupdate

! Update charges

      do i=1,natom
         j = grplist(i)
         if (j .ne. 0) then
            x1 = dsin(ph_theta(j))**2
            if (sp_grp(j) .le. 0) then
               x2 = one
            else
               x2 = dsin(ph_theta(j+1))**2
            end if
            qprot = x2*qstate1(1,i) + (one-x2)*qstate1(2,i)
            qunprot = x2*qstate2(1,i) + (one-x2)*qstate2(2,i)
            charge_phmd(i) = (one-x1)*qprot + x1 * qunprot

            ! Charge update for MD
            qprot = x2*qstate1_md(1,i) + (one-x2)*qstate1_md(2,i)
            qunprot = x2*qstate2_md(1,i) + (one-x2)*qstate2_md(2,i)
            charge(i) = (one-x1)*qprot + x1 * qunprot
         end if
         j = linkgrplist(i)
         ! second
         if (j .ne. 0) then
            x1 = dsin(ph_theta(j))**2
            charge_phmd(i) = x1 * qstate1(1, i) + (one - x1) * qstate2(1,i)
            charge(i) = x1 * qstate1_md(1, i) + (one - x1) * qstate2_md(1,i)
         end if
      end do
      if (add_energies) then
        do i = 1, ntitr
          kin_ene = kin_ene + 0.5 * qmass_phmd * vph_theta(i) ** 2
        end do
      end if

   end if ! QSetLam
#ifdef MPI
   endif
#endif
#endif
   !-------------------------------------------------
   ! Write data to the phmdout file
#ifdef MPI
   if(master) then
#endif
   if (prlam .and. ( prnlev .ge. 0 ) .and. ( istp .eq. 0 )) then
      ! Print RES,TITR,PKA INFO
      write(fmt00,'("(a,",i4,"i5)")') ntitr
      write(fmt120,'("(a, ",i4,"(f8.3))")') ntitr
      write(phmdout_unit,fmt00) '# ititr ',(i,i=1,ntitr)
      write(phmdout_unit,fmt00) '#  ires ',(titrres(i),i=1,ntitr)
      write(phmdout_unit,fmt00) '# itauto',(sp_grp(i),i=1,ntitr)
      write(phmdout_unit,fmt120) '# ParK',(park(i),i=1,ntitr)
   end if

   if ((istp .ne. 0) .and. (mod(istp,nprint_phmd) .eq. 0)) then
      write(fmt100,'("(i8, ",i4,"(f5.2))")') ntitr
      write(fmt110,'("(i8, ",i4,"(f10.4))")') ntitr
      write(fmt115,'("(i8, ",i4,"(f10.4,f10.4))")') ntitr
      if (prdudls) then
        write(phmdout_unit, fmt115) istp, (dsin(ph_theta(i))**2, &
              dudls(i),i=1,ntitr)
      else if (prderiv .and. prnlev .gt. 2) then
         write(phmdout_unit,fmt115) istp, (ph_theta(i), &
           dph_theta(i),i=1,ntitr)
      else if (prlam .and. prnlev .gt. 2) then
         write(phmdout_unit,fmt100) istp, (dsin(ph_theta(i))**2,i=1,ntitr)
      else if (prke .and. prnlev .gt. 2) then
         write(phmdout_unit, fmt115) istp, &
             (0.5 * qmass_phmd * vph_theta(i) ** 2, i = 1, ntitr)
      end if
   end if

#ifdef MPI
   end if !master
#endif

#ifndef CUDA
#ifdef MPI
   call mpi_bcast(charge, natom, mpi_double_precision, 0, pmemd_comm, ierr )
   call mpi_bcast(charge_phmd, natom, mpi_double_precision, 0, &
                 pmemd_comm, ierr )
   call mpi_bcast(ph_theta, MAX_TITR_STATES, mpi_double_precision, 0, &
                  pmemd_comm, ierr )
#endif
#endif

   return

end subroutine updatephmd

subroutine runphmd(x,f,iac,ico,numex,natex,natbel)
   !***************************************************************
   ! Calculates energy and derivatives for PHMD
   ! Charges on titratable group are functions of lambda
   ! para(ntitr) : First parameter in the 1-d model potential
   ! parb(ntitr) : Second parameter in the 1-d model potential
   ! r1,r2,r3,r4,r5,r6 are the parameters in the 2-d model potential
   ! Notice the negative sign in front of r4 in the model potential
   ! qstate1 : Charges of system in state 1
   ! qstate2 : Charges of system in state 2
   ! ph_energy : Energy of pH constraints
   ! ph_theta : States of titratable groups
   ! dph_theta : Derivatives
   !***************************************************************
   use mdin_ctrl_dat_mod, only : iphmd, solvph, temp0, gb_kappa
   use remd_mod, only: remd_method
   use prmtop_dat_mod, only : atm_gb_fs,atm_gb_radii
   use gb_ene_mod, only : calc_born_radii,reff
#ifdef MPI
   use parallel_dat_mod, only : master,numtasks,mytaskid
#endif

   implicit none

   ! Passed parameters
   double precision x,f
   integer iac,ico,numex,natex,natbel
  
   dimension x(3,*),f(3,*), &
             iac(*),ico(*),numex(*),natex(*)

   integer i

   double precision r1,r2,r3,r4,r5,r6,rr
   double precision x1,x2,lambda
   double precision umod,uph,ubarr,dumod,duph,dubarr
   double precision :: ph_energy

   ! Compute energy associated with pH Function
   !********************************************
   ! initialize
   dph_theta(1:ntitr) = zero

   ! Model,PH and Barrier potentials
   model_ene: do i=1,ntitr
#ifdef MPI
      if (.not. master) then
         exit model_ene
      end if
#endif
      if (sp_grp(i) .eq. 0) then ! single site
         lambda = dsin(ph_theta(i))**2
         ubarr = four*barr(i)*(lambda-half)**2
         dubarr = four*two*barr(i)*(lambda-half)
         uph = park(i)*lambda
         duph = park(i)
         if (remd_method .eq. 1) then
             uph = uph * temp0 / temp_phmd
             duph = duph * temp0 / temp_phmd
         end if
         umod = para(i)*(lambda-parb(i))**2
         dumod = two*para(i)*(lambda-parb(i))
         dph_theta(i) = -dumod - dubarr + duph

      else if (sp_grp(i) .eq. 1 .or. sp_grp(i) .eq. 3) then
         umod = zero
         ubarr = zero
         uph = zero

      else if (sp_grp(i) .eq. 2) then   ! double site pKa1 .ne. pka2
         lambda = dsin(ph_theta(i-1))**2
         x1 = dsin(ph_theta(i))**2
         x2 = ONE - x1
         r2 = -two*para(i)*parb(i)
         r1 = -two*para(i-1)*parb(i-1) -r2
         r3 = -two*sp_par(1,i)*sp_par(2,i) -r1
         ubarr = four*barr(i)*(x1-half)**2 + four*barr(i-1)*(lambda-half)**2
         uph = lambda*(park(i-1)*x1 + park(i)*x2)
         umod = sp_par(1,i)*lambda**2*x1**2 +r1*lambda*x1 &
              + r2*lambda + r3*lambda**2*x1 +para(i)*lambda**2
         if (remd_method .eq. 1) then
           uph = uph * temp0 / temp_phmd
         end if

         ! dU/dLamb
         !Barr
         dubarr = two*four*barr(i-1)*(lambda-half)
         !pH 
         duph = park(i-1)*x1 + park(i)*x2
         if (remd_method .eq. 1) then
           uph = uph * temp0 / temp_phmd
           duph = duph * temp0 / temp_phmd
         end if
         !mod
         dumod = two*sp_par(1,i)*lambda*x1**2 +r1*x1 &
              + r2 + two*r3*lambda*x1 +two*para(i)*lambda

         dph_theta(i-1) = -dumod - dubarr + duph

         ! dU/dx
         !Barr
         dubarr = two*four*barr(i)*(x1-half)
         !pH
         duph = (park(i-1)-park(i))*lambda
         if (remd_method .eq. 1) then
           duph = duph * temp0 / temp_phmd
         end if
         !mod
         dumod = two*sp_par(1,i)*lambda**2*x1 +r1*lambda &
              + r3*lambda**2

         dph_theta(I)= - dumod - dubarr + duph

      else if (sp_grp(i) .eq. 4) then  ! double site pKa1 .eq. pka2
         lambda = dsin(ph_theta(i-1))**2
         x1 = dsin(ph_theta(i))**2
         x2 = one - x1

         r1 = parmod(1,i)
         r2 = parmod(2,i)
         r3 = parmod(3,i)
         r4 = parmod(4,i)
         r5 = parmod(5,i) -r1*r4**2
         r6 = -two*parmod(5,i)*parmod(6,i) -r2*r4**2

         ubarr = four*barr(i)*(x1-half)**2 + four*barr(i-1)*(lambda-half)**2
         uph = lambda*park(i)
         if (remd_method .eq. 1) then
           uph = uph * temp0 / temp_phmd
         end if
         umod = (r1*lambda**2 + r2*lambda + r3)*(x1-r4)**2 &
              + r5*lambda**2 + r6*lambda

         !dU/dLamb
         !Barr
         dubarr = four*two*barr(i-1)*(lambda-half)
         !pH
         dupH = park(i-1)
         if (remd_method .eq. 1) then
            uph = uph * temp0 / temp_phmd
            duph = duph * temp0 / temp_phmd
         end if
         !mod
         dumod = (two*r1*lambda + r2)*(x1-r4)**2 + two*r5*lambda + r6
         dph_theta(i-1) = -dumod - dubarr + duph

         !dU/dx
         !Barr
         dubarr = four*two*barr(i)*(x1-half)
         !pH
         dupH = zero
         !mod
         dumod = two*(r1*lambda**2 + r2*lambda + r3)*(x1-r4)
         dph_theta(I) = - dumod - dubarr + dupH

      end if

      ph_energy = ph_energy - umod - ubarr + upH

   end do model_ene ! NTitr

   if (qsetlam) dph_theta(1:ntitr) = zero


   !  Compute the energy derivatives wrt lambda of the
   !  Non 1-4 nonbonded electrostatics and vdw energies
   !  and GB solvation free energies
   !*******************************************************
   if ( iphmd .eq. 1 ) then
      call gb_phmd(x,f,reff,iac,ico,numex,natex,natbel,natom)
   end if
   return
end subroutine runphmd

subroutine gb_phmd(x,f,reff,iac,ico,numex,natex,natbel,natom)
  !***********************************************************
  ! Compute the energy derivatives wrt lambda for the non-1-4
  ! electrostatic energy, vdw energy, and GB solvation free
  ! energy. Adjust the vdw forces in GB model.
  !***********************************************************

   use mdin_ctrl_dat_mod, only: alpb,iphmd,arad,extdiel,intdiel,gb_kappa
   use prmtop_dat_mod, only : gbl_cn1,gbl_cn2,ntypes
#ifdef MPI
   use parallel_dat_mod, only : numtasks,mytaskid,master
#endif

   implicit none

   ! passed parameters
   double precision x,f,reff
   integer iac,ico,numex,natex,natbel,natom
   double precision, parameter   :: alpb_alpha = 0.571412d0
 
   ! local parameters

   double precision extdieli,intdieli, &
          xi,yi,zi,ri,four_ri,xij,yij,zij, &
          xj,yj,zj,rij,qi,r2,reff_j, &
          onekappa,dl,e,expmkf,fgbi,fgbk,qiqj, &
          rinv,temp1,f6,f12,qi2h,qid2h,r6inv,r2inv

#ifdef HAS_10_12
   double precision :: r10inv,f10
#endif

   double precision :: alpb_beta, one_arad_beta

   integer i,j,k,iexcl,iaci,ic,maxi,jexcl,jexcl_last,jjv, &
           ier, icount 

   logical, dimension(:), allocatable :: skipv
   integer, dimension(:), allocatable :: temp_jj
   double precision, dimension(:), allocatable :: r2x,rjx,vectmp1,vectmp2, &
                         vectmp3,vectmp4,vectmp5               

   double precision x1,x2
   double precision lambdah, lambdag, lambda, lambdah2, lambdag2
   double precision qunprot,qprot,qg,qh,qxh,qxg,fact,xh,xg
   double precision dxh,dxg,facth,factg
   double precision radg, radh
   logical lhtitr, lgtitr, lhtauto, lgtauto
   integer h
   integer, dimension(:) :: g(natom)
   integer ihtaut,igtaut
   
   double precision :: tmpx,tmpy,tmpz
   double precision :: tx,ty,tz
   double precision :: fvdw,cut2

   dimension x(3,*),f(3,*),reff(natom), &
             iac(*),ico(*),numex(*),natex(*)


   allocate( r2x(natom),rjx(natom),vectmp1(natom), vectmp2(natom), &
             vectmp3(natom),vectmp4(natom),vectmp5(natom), &
             skipv(0:natom), temp_jj(natom), stat = ier )

   onekappa = zero

   if(alpb == 1) then
      ! Sigalov Onufriev ALPB (epsilon-dependent GB):
      alpb_beta = alpb_alpha*(intdiel/extdiel)
      extdieli = one/(extdiel*(one + alpb_beta))
      intdieli = one/(intdiel*(one + alpb_beta))
      one_arad_beta = alpb_beta/arad
      if (gb_kappa/=zero) onekappa = one/gb_kappa
   else
   !  Standard Still's GB - alpb=0
      extdieli = one/extdiel
      intdieli = one/intdiel
      one_arad_beta = zero
   end if

   !---------------------------------------------------------------------------
   !
   ! The effective Born radii are calculated during dynamics
   !
   !---------------------------------------------------------------------------

   iexcl = 1

   cut2 = phmdcut * phmdcut
#ifdef MPI
   do i=1,mytaskid
      iexcl = iexcl + numex(i)
   end do
#endif

   maxi = natom
   if(natbel > 0) maxi = natbel

   !--------------------------------------------------------------------------
   !
   !     Compute dU/dlambda from GB solvation, electrostatics, and vdw, as
   !     well as correcting the vdw forces on atom pairs i,j
   !
   !--------------------------------------------------------------------------

#ifdef MPI
    do i = mytaskid + 1, maxi, numtasks
#else
    do i=1,maxi
#endif

      tmpx = zero
      tmpz = zero
      tmpy = zero

      h = grplist(i)

      xi = x(1,i)
      yi = x(2,i)
      zi = x(3,i)
      qi = charge_phmd(i)
      ri = reff(i)
      four_ri = four*reff(i)
      iaci = ntypes * (iac(i) - 1)
      jexcl = iexcl
      jexcl_last = iexcl + numex(i) -1
   
     !         -- check the exclusion list for eel and vdw:

      do k=i+1,natom
         skipv(k) = .false.
      end do
      do jjv=jexcl,jexcl_last
         skipv(natex(jjv))=.true.
      end do

      icount = 0
      do j=i+1,natom
         if( h .gt. 0 .or. grplist(j) .gt. 0 ) then 
            xij = xi - x(1,j)
            yij = yi - x(2,j)
            zij = zi - x(3,j)
            r2 = xij*xij + yij*yij + zij*zij
          
            if( r2 <= cut2 ) then

              reff_j = reff(j)

              icount = icount + 1
              g(icount) = grplist(j)
              temp_jj(icount) = j
              r2x(icount) = r2

              rjx(icount) = reff_j

            end if !r2 <= cut
         end if !atom is titrating
      end do

      vectmp1(1:icount) = four_ri*rjx(1:icount)

      call vdinv( icount, vectmp1, vectmp1 ) !Invert things
         
      vectmp1(1:icount) = -r2x(1:icount)*vectmp1(1:icount)
         
      call vdexp( icount, vectmp1, vectmp1 )
      ! [ends up with Exp(-rij^2/[4*ai*aj])]

      vectmp3(1:icount) = r2x(1:icount) + rjx(1:icount)*ri*vectmp1(1:icount)
      ! [ends up with fij]

      call vdinvsqrt( icount, vectmp3, vectmp2 ) !vectmp2 = 1/fij

         
      if( gb_kappa /= zero ) then
         call vdinv( icount, vectmp2, vectmp3 )
         vectmp3(1:icount) = -gb_kappa*vectmp3(1:icount)
         call vdexp( icount, vectmp3, vectmp4 ) !exp(-kappa*fij)
      end if

      call vdinvsqrt( icount, r2x, vectmp5 ) !1/rij

      ! vectmp1 = Exp(-rij^2/[4*ai*aj])
      ! vectmp2 = 1/fij
      ! vectmp3 = -kappa*fij - if kappa/=zero, otherwise = fij
      ! vectmp4 = exp(-kappa*fij)
      ! vectmp5 = 1/rij

      !---- Start first outer loop ----

      do k=1,icount

         if ( ( h .gt. 0 ) .or. ( g(k) .gt. 0 ) ) then
            j = temp_jj(k)

            if (h .gt. 0) then
               if (sp_grp(h) .le. 0) then
                  x2 = one
               else
                  x2 = dsin(ph_theta(h+1))**2
                  lambda = dsin(ph_theta(h))**2

                  qunprot = lambda*(qstate2(1,i)-qstate2(2,i))
                  qprot = (1-lambda)*(qstate1(1,i)-qstate1(2,i))
                  qxh = charge_phmd(j)*(qunprot+qprot)
               end if
               qunprot = x2*qstate2(1,i) +(1-x2)*qstate2(2,i)
               qprot = x2*qstate1(1,i) + (1-x2)*qstate1(2,i)
               qh = charge_phmd(j)*(qunprot - qprot)

            end if

            if (g(k) .gt. 0) then
               if (sp_grp(g(k)) .le. 0) then
                  x2 = one
               else
                  x2 = dsin(ph_theta(g(k)+1))**2
                  lambda = dsin(ph_theta(g(k)))**2
                  qunprot = lambda*(qstate2(1,j)-qstate2(2,j))
                  qprot =(1-lambda)*(qstate1(1,j)-qstate1(2,j))

                  qxg = qi*(qunprot + qprot)
               endif
               qunprot = x2*qstate2(1,j) + (1-x2)*qstate2(2,j)
               qprot = x2*qstate1(1,j) + (1-X2)*qstate1(2,j)
               qg = qi*(qunprot - qprot)
            endif
            xij = xi - x(1,j)
            yij = yi - x(2,j)
            zij = zi - x(3,j)
            r2 = r2x(k)
            qiqj = qi * charge_phmd(j)

            if( gb_kappa == zero ) then
               fgbk = zero
               expmkf = extdieli
            else
               expmkf = vectmp4(k)*extdieli


               fgbk = vectmp3(k)*expmkf !-kappa*fij*exp(-kappa*fij)/Eout
               if(alpb == 1) then ! Sigalov Onufriev ALPB:
                  fgbk = fgbk+(fgbk*one_arad_beta*(-vectmp3(k)*onekappa))

                    ! (-kappa*fij*exp(-kappa*fij)(1 + fij*ab/A)/Eout)*(1/fij+ab/A)
                    ! Note: -vectmp3(k)*onekappa = fij
               end if
            end if
            dl = intdieli - expmkf

            fgbi = vectmp2(k)  !1/fij

            temp1 = -dl*(fgbi + one_arad_beta)

            if (g(k) .gt. 0) then
               dph_theta(g(k)) = dph_theta(g(k)) + qg*temp1

               if (sp_grp(g(k)) .gt. 0) then
                  dpH_Theta(G(k)+1) = dpH_Theta(G(k)+1) +  &
                                      QXG*temp1
               endif
            endif

            if (h .gt. 0) then
               dph_theta(h) = dph_theta(h) + qh*temp1
               if (sp_grp(h) .gt. 0) then
                  dph_theta(h+1) = dph_theta(h+1) +  &
                                   qxh*temp1
               endif
            endif

            if( .not. skipv(j) ) then

               !   -- gas-phase Coulomb energy:

               !we can use the cached values.
               rinv = vectmp5(k) !1/rij

               temp1 = intdieli*rinv

              if (g(k) .gt. 0) then
                 dph_theta(g(k)) = dph_theta(g(k)) + qg*temp1

                 if (sp_grp(g(k)) .gt. 0) then
                    dph_theta(g(k)+1) = dph_theta(g(k)+1) +  &
                        qxg*temp1
                 endif
              endif

              if (h .gt. 0) then
                 dph_theta(h) = dph_theta(h) + qh*temp1
                 if (sp_grp(h) .gt. 0) then
                    dph_theta(h+1) = dph_theta(h+1) +  &
                                     qxh*temp1
                 end if
              end if
              !    -- van der Waals energy:
              ! Set up some variables

              lhtitr = .false.
              lhtauto = .false.
              xh = one
              ihtaut = 0
              if (h .gt. 0) then
                 radh = vstate1(1,i)
                 lambdah = dsin(ph_theta(h))**2
                 lambdah2 = one - lambdah
                 facth = lambdah2
                 if (sp_grp(h) .eq. 1 .or. sp_grp(h) .eq. 3) then
                    lhtauto = .true.
                    x2 = dsin(ph_theta(h+1))**2
                    radh = radh + vstate1(2,i)
                    if (vstate1(1,i) .gt. zero) then
                       ihtaut = 1
                       xh = x2
                    else if (vstate1(2,i) .gt. zero) then
                       ihtaut = -1
                       xh = one - x2
                    end if
                    if (sp_grp(h) .eq. 1) then
                       facth = one - lambdah*xh
                       dxh = -ihtaut*lambdah
                    else if (sp_grp(h) .eq. 3) then
                       facth = lambdah2*xh
                       dxh = ihtaut*lambdah2
                    end if
                 end if
                 if (radh .gt. zero) lhtitr = .true.
              end if


              lgtitr = .false.
              lgtauto = .false.
              xg = one
              igtaut = 0
              if (g(k) .gt. 0) then
                 radg = vstate1(1,j)
                 lambdag = dsin(ph_theta(g(k)))**2
                 lambdag2 = one - lambdag
                 factg = lambdag2
                 if (sp_grp(g(k)).eq.1 .or. sp_grp(g(k)).eq.3) then
                    lgtauto = .true.
                    radg =  radg + vstate1(2,j)
                    x2 = dsin(ph_theta(g(k)+1))**2

                    if (vstate1(1,j) .gt. zero) then
                       igtaut = 1
                       xg = x2
                    else if (vstate1(2,j) .gt. zero) then
                       igtaut = -1
                       xg = one - x2
                    end if
                    if (sp_grp(g(k)).eq.1) then
                       factg = one - lambdag*xg
                       dxg = -igtaut*lambdag
                    else if (sp_grp(g(k)).eq.3) then
                       factg = lambdag2*xg
                       dxg = igtaut*lambdag2
                    end if
                 end if
                 if (radg .gt. zero) lgtitr = .true.
              end if

            
              if (lhtitr .or. lgtitr) then
 
                 ic = ico( iaci + iac(j) )

                 if( ic .gt. 0 ) then
                    ! 6-12 potential:
                    r2inv = rinv*rinv
                    r6inv = r2inv*r2inv*r2inv
                    f6 = gbl_cn2(ic)*r6inv
                    f12 = gbl_cn1(ic)*(r6inv*r6inv)

                    e = f12 - f6     ! regular LJ vdw energy
                    fvdw = 6*(2*f12-f6)*r2inv ! vdw force 

#ifdef HAS_10_12
                 !
                 ! The following could be commented out if the Cornell
                 ! et al. force field was always used, since then all hbond
                 ! terms are zero.

                 else if(ic .lt. 0)
                  ! 10-12 potential:
                    r10inv = r2inv*r2inv*r2inv*r2inv*r2inv
                    f10 = bsol(-ic)*r10inv
                    f12 = asol(-ic)*r10inv*r2inv
                    e = f12 -f10    ! regular 12-10 vdw energy
                    fvdw = (12*f12-10*f10)*r2inv ! vdw force
#endif
                 end if
                 ! Compute lambda energy and force of vdw for two and one titrating H's
                 if( ic .ne. 0 ) then
                    fact = one
                    if (lhtitr .and. lgtitr) then
                       if (h .ne. g(k)) then
                          fact = facth*factg

                          dph_theta(h) = dph_theta(h) - xh*factg*e
                          dph_theta(g(k)) = dph_theta(g(k)) - xg*facth*e

                          if (lhtauto) then
                             dph_theta(h+1)=dph_theta(h+1) + dxh*factg*e
                          end if
                          if (lgtauto) then
                             dph_theta(g(k)+1)=dph_theta(g(k)+1) + dxg*facth*e
                          end if

                       else if (sp_grp(h).eq.1) then
                          fact = one - lambdah
                          dph_theta(h) = dph_theta(h) - e
                       end if
                    else if (lhtitr) then
                       fact = facth
                       dph_theta(h) = dph_theta(h) - xh*e
                       if (lhtauto) then
                          dph_theta(h+1) = dph_theta(h+1) + dxh*e
                       end if
                    else if (lgtitr) then
                       fact = factg
                       dph_theta(g(k)) = dph_theta(g(k)) - xg*e
                       if (lgtauto) then
                          dph_theta(g(k)+1) = dph_theta(g(k)+1) + dxg*e
                       end if
                    end if
                    fact = fact-one

                    if( iphmd .ne. 2 ) then
                       ! Force on spatial coordinates
                       tx = Fact*fvdw*xij
                       ty = Fact*fvdw*yij
                       tz = Fact*fvdw*zij
                    
                       ! vdw force on atom j
                       f(1,j) = f(1,j) + tx
                       f(2,j) = f(2,j) + ty
                       f(3,j) = f(3,j) + tz

                       tmpx = tmpx + tx
                       tmpy = tmpy + ty
                       tmpz = tmpz + tz
                    
                    endif
                 endif ! if vdw
              endif  ! titratable group       
            endif ! ( .not. skipv(j) ) for gas-phase vdw and elec
         endif ! ( H .gt. 0 ) .or. ( G .gt. 0 )
      enddo ! k=1,icount, for epol and gas-phase vdw and elec

       ! vdw force on atom i     
      if( iphmd .ne. 2 ) then
         f(1,i) = f(1,i) - tmpx
         f(2,i) = f(2,i) - tmpy
         f(3,i) = f(3,i) - tmpz
      endif 

#ifdef MPI
      do k = i, min(i + numtasks - 1, natom)
         iexcl = iexcl + numex(k)
      end do
#else
      iexcl=iexcl+numex(i)
#endif

    enddo ! i=1,maxi

  !---- End first outer loop ----

      !  -- diagonal epol term

#ifdef MPI

   do i=mytaskid + 1, natom, numtasks
#else
   do i=1,natom
#endif
      qi = charge_PHMD(i)
      h = grplist(i)

      if (h .gt. 0) then
         if (sp_grp(h) .le. 0) then
            x2 = one
         else
            x2 = dsin(ph_theta(h+1))**2
            lambda = dsin(ph_theta(h))**2

            qunprot = lambda*(qstate2(1,i)-qstate2(2,i))
            qprot = (1-lambda)*(qstate1(1,i)-qstate1(2,i))
            qxh = qi*(qunprot+qprot)
         endif
         qunprot = x2*qstate2(1,i) +(1-x2)*qstate2(2,i)
         qprot = x2*qstate1(1,i) + (1-x2)*qstate1(2,i)
         qh = qi*(qunprot - qprot)
         expmkf = exp( -gb_kappa * reff(i) )*extdieli
         dl = intdieli - expmkf
         qi2h = half*qi*qi
         qid2h = qi2h * dl

         temp1 = (one / reff(i) + one_arad_beta)

         dph_theta(H) = dph_theta(h) - qh*dl*temp1

         if (sp_grp(h) .gt. 0) then
            dph_theta(h+1) = dph_theta(h+1) -  &
                             qxh*dl*temp1
         endif
      endif

   enddo

   if( allocated( r2x ) ) then
      deallocate(skipv, temp_jj, r2x,rjx,vectmp1,vectmp2, & 
                 vectmp3,vectmp4,vectmp5, stat = ier )
   endif
   
   return
end subroutine gb_phmd

subroutine phmd14nb(x,f,i,j,crfac,f12,f6,scnb0,r2inv)
  !***********************************************************
  ! Compute the derivatives wrt to lambda for the 1-4 terms
  ! and the corrections to the 1-4 vdw force
  !***********************************************************

   use mdin_ctrl_dat_mod, only : iphmd
   implicit none

   integer, intent(in) :: i,j
   double precision, intent(in) :: crfac,f12,f6,scnb0,r2inv
   double precision f,x

   integer h,g
   double precision qg,qh,qprot,qunprot,qxg,qxh,x2,lambda

   double precision lambdah, lambdag, lambdah2, lambdag2
   double precision xh,xg
   double precision dxh,dxg,facth,factg
   double precision fact
   double precision radg, radh
   logical lhtitr, lgtitr, lhtauto, lgtauto
   integer ihtaut,igtaut

   double precision :: tx,ty,tz
   double precision :: fvdw
   double precision :: e

   dimension x(3,*),f(3,*)

   h = GrpList(i)
   g = GrpList(j)

   e    = 0

   ! 1-4 ELE
   if (h .gt. 0) then
      if (sp_grp(h) .le. 0) then
         x2 = one
      else
         x2 = dsin(ph_theta(h+1))**2
         lambda = dsin(ph_theta(h))**2

         qunprot = lambda*(qstate2(1,i)-qstate2(2,i))
         qprot = (1-lambda)*(qstate1(1,i)-qstate1(2,i))
         qxh = charge_phmd(j)*(qunprot+qprot)
      endif
      qunprot = x2*qstate2(1,i) +(1-x2)*qstate2(2,i)
      qprot = x2*qstate1(1,i) + (1-x2)*qstate1(2,i)
      qh = charge_phmd(j)*(qunprot - qprot)
      dph_theta(h) = dph_theta(h) + qh*crfac
      if (sp_grp(h) .gt. 0) then
         dph_theta(h+1) = dph_theta(h+1) +  &
                          qxh*crfac
      endif
   endif

   if (g .gt. 0) then
      if (sp_grp(g) .le. 0) then
         x2 = one
      else
         x2 = dsin(ph_theta(g+1))**2
         lambda = dsin(ph_theta(g))**2
         qunprot = lambda*(qstate2(1,j)-qstate2(2,j))
         qprot =(1-lambda)*(qstate1(1,j)-qstate1(2,j))
         qxg = charge_phmd(i)*(qunprot + qprot)
      endif
         qunprot = x2*qstate2(1,j) + (1-x2)*qstate2(2,j)
         qprot = x2*qstate1(1,j) + (1-x2)*qstate1(2,j)
         qg = charge_phmd(i)*(qunprot - qprot)
         dph_theta(g) = dph_theta(g) + qg*crfac
         if (sp_grp(g) .gt. 0) then
            dph_theta(g+1) = dph_theta(g+1) +  &
                             QXG*crfac
         endif
   endif

   ! 1-4 VDW
   lhtitr = .false.
   lhtauto = .false.
   xh = one
   ihtaut = 0
   if (h .gt. 0) then
      radh = vstate1(1,i)
      lambdah = dsin(ph_theta(h))**2
      lambdah2 = one - lambdah
      facth = lambdah2
      if (sp_grp(h) .eq. 1 .or. sp_grp(h) .eq. 3) then
         lhtauto = .true.
         x2 = dsin(ph_theta(h+1))**2
         radh = radh + vstate1(2,i)
         if (vstate1(1,i) .gt. zero) then
            ihtaut = 1
            xh = x2
         else if (vstate1(2,i) .gt. zero) then
            ihtaut = -1
            xh = one - x2
         endif
         if (sp_grp(h) .eq. 1) then
            facth = one - lambdah*xh
            dxh = -ihtaut*lambdah
         else if (sp_grp(h) .eq. 3) then
            facth = lambdah2*xh
            dxh = ihtaut*lambdah2
         endif
      end if
      if (radh .gt. zero) lhtitr = .true.
   end if

   lgtitr = .false.
   lgtauto = .false.
   xg = one
   igtaut = 0
   if (g .gt. 0) then
      radg = vstate1(1,j)
      lambdag = dsin(ph_theta(g))**2
      lambdag2 = one - lambdag
      factg = lambdag2
      if (sp_grp(g).eq.1 .or. sp_grp(g).eq.3) then
          lgtauto = .TRUE.
          radg = radg + vstate1(2,j)
          x2 = dsin(ph_theta(g+1))**2

          if (vstate1(1,j) .gt. zero) then
             igtaut = 1
             xg = x2
          else if (vstate1(2,j) .gt. zero) then
             igtaut = -1
             xg = one - x2
          end if
          if (sp_grp(g).eq.1) then
             factg = one - lambdag*xg
             dxg = -igtaut*lambdag
          else if (sp_grp(g).eq.3) then
             factg = lambdag2*xg
             dxg = igtaut*lambdag2
          end if
      end if
      if (radg .gt. zero) lgtitr = .true.
   end if

   if (lhtitr .or. lgtitr) then

      e = (f12 - f6)*scnb0     ! regular LJ vdw energy
      fvdw = 6*(2*f12-f6)*r2inv*scnb0 ! vdw force 

      ! Compute lambda energy and force of 1-4 vdw for titrating H's

      fact = one
      if (lhtitr .and. lgtitr) then
         if (h .ne. g) then
            fact = facth*factg

            dph_theta(h) = dph_theta(h) - xh*factg*e
            dph_theta(g) = dph_theta(g) - xg*facth*e

            if (lhtauto) then
               dph_theta(h+1)=dph_theta(h+1) + dxh*factg*e
            end if
            if (lgtauto) then
               dph_theta(g+1)=dph_theta(g+1) + dxg*facth*e
            end if

         else if (sp_grp(h).eq.1) then
            fact = one - lambdah
            dph_theta(h) = dph_theta(h) - e
         endif
      else if (lhtitr) then
         fact = facth
         dph_theta(h) = dph_theta(h) - xh*e
         if (lhtauto) then
            dph_theta(h+1) = dph_theta(h+1) + dxh*e
         end if
      else if (lgtitr) then
         fact = factg
         dph_theta(g) = dph_theta(g) - xg*e
         if (lgtauto) then
            dph_theta(g+1) = dph_theta(g+1) + dxg*e
         end if
      end if
      fact = fact-one

      ! Force on spatial coordinates
      tx = fact*fvdw*(x(1,i)-x(1,j))
      ty = fact*fvdw*(x(2,i)-x(2,j))
      tz = fact*fvdw*(x(3,i)-x(3,j))
    
      ! vdw force on atom j
      f(1,j) = f(1,j) + tx
      f(2,j) = f(2,j) + ty
      f(3,j)   = f(3,j)   + tz

      ! vdw force on atom i     
      f(1,i) = f(1,i) - tx
      f(2,i) = f(2,i) - ty
      f(3,i)   = f(3,i)   - tz

   end if
     
   return
end subroutine phmd14nb

subroutine convphmd()
  !***********************************************************
  ! Convert from dU/dlambda to dU/dtheta
  !***********************************************************

   implicit none
   integer i

   do i = 1, ntitr
      dph_theta(i) = dph_theta(i)*dsin(two*ph_theta(i))
   end do

   return
end subroutine convphmd


subroutine phmdwriterestart()
  !***********************************************************
  ! Write the phmd restart file
  !***********************************************************
   use file_io_dat_mod, only: PHMD_UNIT, phmdrestrt_name, owrite, phmdout_unit
   implicit none

   ! Variable descriptions
   !

   logical, save     :: first_cprestrt = .true.
   character(len=7)  :: stat

   namelist /phmdrst/ pH_Theta,vpH_Theta ! status of phmd

   if (first_cprestrt) then
      if (owrite == 'N') then
         stat = 'NEW'
      else if (owrite == 'O') then
         stat = 'OLD'
      else if (owrite == 'R') then
         stat = 'REPLACE'
      else if (owrite == 'U') then
         stat = 'UNKNOWN'
      end if
      open(unit=phmd_unit, file=phmdrestrt_name, status=stat, form='FORMATTED', &
           delim='APOSTROPHE')
      first_cprestrt = .false.
   else
      open(unit=phmd_unit, file=phmdrestrt_name, status='OLD', form='FORMATTED', &
           delim='APOSTROPHE')
   end if

   write(phmd_unit, nml=phmdrst)
   close(phmd_unit)

   ! flush all cpout data
   close(phmdout_unit)
   open(unit=phmdout_unit, file=phmdout_name, status='OLD',position='APPEND')
  
   return
end subroutine phmdwriterestart

#ifdef MPI
subroutine phmd_bcast(ierr)
   !***********************************************************
   ! Broadcast master information to all threads
   !***********************************************************
   use mdin_ctrl_dat_mod
   use parallel_dat_mod, only : mpi_integer, mpi_double_precision, &
                                pmemd_comm, master, mpi_logical
   use gb_ene_mod, only : excludeatoms
   implicit none
   integer, intent(out) :: ierr
   if( ( iphmd .gt. 0 ) .and. ( .not. master ) ) then
      allocate(qstate1(2,natom), qstate2(2,natom), vstate1(2,natom), & 
            vstate2(2,natom), dph_theta(MAX_TITR_STATES), &
            dph_avg(MAX_TITR_STATES), para(MAX_TITR_STATES), &
            parb(MAX_TITR_STATES), parmod(6,MAX_TITR_STATES), &
            park(MAX_TITR_STATES), sp_par(2,MAX_TITR_STATES), &
            sp_grp(MAX_TITR_STATES), &
            titrres(MAX_TITR_STATES*2),grplist(NAtom), &
            local_model_pka(MAX_TITR_STATES), qstate1_md(2,NAtom), &
            qstate2_md(2,NAtom), charge_phmd(natom), excludeatoms(natom), & 
            stat=ierr) ! if successfully allocated, ier=0
   endif
   if( iphmd .gt. 0 ) then
      call mpi_bcast(qstate1, 2*natom, mpi_double_precision, 0, pmemd_comm, &
                     ierr )
      call mpi_bcast(qstate2, 2*natom, mpi_double_precision, 0, pmemd_comm, &
                     ierr )
      call mpi_bcast(qstate1_md, 2*natom, mpi_double_precision, 0, &
                     pmemd_comm, ierr )
      call mpi_bcast(qstate2_md, 2*natom, mpi_double_precision, 0, &
                     pmemd_comm, ierr )
      call mpi_bcast(charge_phmd, natom, mpi_double_precision, 0, &
                     pmemd_comm, ierr )
      call mpi_bcast(excludeatoms, natom, mpi_integer, 0, &
                     pmemd_comm, ierr )
      call mpi_bcast(vstate1, 2*natom, mpi_double_precision, 0, &
                     pmemd_comm, ierr )
      call mpi_bcast(vstate2, 2*natom, mpi_double_precision, 0, &
                     pmemd_comm, ierr )
      call mpi_bcast(dph_theta, MAX_TITR_STATES, mpi_double_precision, 0, &
                     pmemd_comm, ierr )
      call mpi_bcast(vph_theta, MAX_TITR_STATES, mpi_double_precision, 0, &
                     pmemd_comm, ierr )
      call mpi_bcast(ph_theta, MAX_TITR_STATES, mpi_double_precision, 0, &
                     pmemd_comm, ierr )
      call mpi_bcast(dph_avg, MAX_TITR_STATES, mpi_double_precision, 0, &
                     pmemd_comm, ierr )
      call mpi_bcast(para, MAX_TITR_STATES, mpi_double_precision, 0, &
                     pmemd_comm, ierr )
      call mpi_bcast(parb, MAX_TITR_STATES, mpi_double_precision, 0, &
                     pmemd_comm, ierr )
      call mpi_bcast(parmod, 6*MAX_TITR_STATES, mpi_double_precision, 0, &
                     pmemd_comm, ierr )
      call mpi_bcast(park, MAX_TITR_STATES, mpi_double_precision, 0, &
                     pmemd_comm, ierr )
      call mpi_bcast(sp_par, MAX_TITR_STATES, mpi_double_precision, 0, &
                     pmemd_comm, ierr )
      call mpi_bcast(barr, MAX_TITR_STATES, mpi_double_precision, 0, &
                     pmemd_comm, ierr )
      call mpi_bcast(sp_grp, MAX_TITR_STATES, mpi_integer, 0, pmemd_comm, &
                     ierr )
      call mpi_bcast(titrres, 2*MAX_TITR_STATES, mpi_integer, 0, pmemd_comm, &
                     ierr )
      call mpi_bcast(grplist, natom, mpi_integer, 0, pmemd_comm, ierr )
      call mpi_bcast(ntitr, 1, mpi_integer, 0, pmemd_comm, ierr )
      call mpi_bcast(nsolute, 1, mpi_integer, 0, pmemd_comm, ierr )
      call mpi_bcast(phmdcut, 1, mpi_double_precision, 0, pmemd_comm, ierr)
      call mpi_bcast(phmd_igb, 1, mpi_integer, 0, pmemd_comm, ierr)
      call mpi_bcast(local_model_pka, MAX_TITR_STATES, mpi_double_precision, &
                     0, pmemd_comm, ierr)
      call mpi_bcast(temp_phmd, 1, mpi_double_precision, 0, pmemd_comm, ierr)
      call mpi_bcast(qsetlam, 1, mpi_logical, 0, pmemd_comm, ierr)
      call mpi_bcast(ntitr, 1, mpi_integer, 0, pmemd_comm, ierr)
      call mpi_bcast(qmass_phmd, 1, mpi_double_precision, 0, pmemd_comm, ierr)
      call mpi_bcast(temp_phmd, 1, mpi_double_precision, 0, pmemd_comm, ierr)
      call mpi_bcast(phbeta, 1, mpi_double_precision, 0, pmemd_comm, ierr)
   endif
end subroutine phmd_bcast
#endif /* MPI */
double precision function double_protonation()
  !***********************************************************
  ! Compute the number of protons on the system for remd.
  ! Note that in this implementation the number of protons 
  ! is generally not an integer.
  !***********************************************************

   implicit none
   integer i
   do i = 1, ntitr
      if(sp_grp(i) .EQ. 0) then
         double_protonation = double_protonation + (1-dsin(ph_theta(i))**2)
      elseif(sp_grp(i) .EQ. 2 .OR. sp_grp(i) .EQ. 4) then
         double_protonation = double_protonation + (1-dsin(ph_theta(i-1))**2)
      endif
   enddo

   return

end function double_protonation

subroutine reset_park()
  !***********************************************************
  ! Reset the park array after remd switches.
  !***********************************************************

   use mdin_ctrl_dat_mod, only : solvph
   implicit none
   integer i
   do i = 1, ntitr
      park(i) = log(10.0D0) * kb * temp_phmd * &
                 (local_model_pka(i) - solvph)
   end do

#ifdef CUDA
   call gpu_reset_park(park)
#endif
   return
end subroutine reset_park

end module phmd_mod

