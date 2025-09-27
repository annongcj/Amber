#include "copyright.i"
!*******************************************************************************
!
! Module: mcres_mod
!
! Description: This code implements Monte Carlo water moves in molecular
! dynamics simulations. It uses a grid based method in order to block out
! sterically occupied regions and thus reduce the number of energy evaluations
! required.  It is intended to speed convergence of simulations of, for example,
! proteins with buried cavities.
! 
! The program was written by Ido Ben-Shalom, Charlie Lin, and Mike Gilson 
!
! Main Variables Description: 
! X, Y, Z: the atomic coordinates 
! xVxl, yVxl, zVxl: the voxel "coordinates" 
! watXVxl, watXVxl, watXVxl: the voxel coordinates of the water molecule
!
!*******************************************************************************

module mcres_mod

  use random_mod, only : random_state, amrand_gen, amrset_gen
  use file_io_dat_mod, only: mdout

  implicit none

  type(random_state),save      :: mcres_randgen   

  integer :: mcres_numatms

! contains
  double precision, parameter      :: GSP = 0.2  ! GRID SPACING
  double precision, parameter      :: WATERRADIUS = 0.9

  integer              :: error,stat,nmEmpVxls

  integer              :: rndmVxl(3)

  integer           :: maxCoarseXVxl, maxCoarseYVxl, maxCoarseZVxl
  double precision  :: coarseGSP
  double precision, allocatable :: steric_radius(:)
  integer, allocatable :: mcwatregionmask(:)

  character(4)                  :: mol_str

contains

!*******************************************************************************
!
! Subroutine:  setup_mcres
!
! Description: Checks if string for water residue exists.  Sets up RNG
!
!*******************************************************************************

subroutine setup_mcres(mol_str)

    use mdin_ctrl_dat_mod
    use mdin_ewald_dat_mod
    use file_io_dat_mod
    use prmtop_dat_mod
    use pmemd_lib_mod
    use mol_list_mod
    use pbc_mod

    implicit none
    integer i
    logical mc_err
    character(4)                  :: mol_str

    call amrset_gen(mcres_randgen, ig)

    ! Checks if you can spell.  Error checking
    ! Can be optimized later by only being used on first call.
    mc_err = .true.    
    mcres_numatms = 0

    mcrescyc = nmc

    do i= 1, nres
        if(gbl_labres(i) .eq. mol_str)  then
            mcres_numatms = gbl_res_atms(i+1) - gbl_res_atms(i)
            mc_err = .false.
        end if
    end do

#ifdef CUDA
    call gpu_coarsegrid_setup(mcres_numatms, ew_coeff, mcrescyc)
#endif

    if(mc_err .eqv. .true.) then
       write(mdout,'(a,a,a)') "Error in definition of mcres - residue with name ", mol_str, " does not exist."
       call mexit(mdout, 1)
    end if

end subroutine setup_mcres

!*******************************************************************************
!
! Module:  init_mcwat_mask
!
! Description: Initializes ligand mask
!             
!*******************************************************************************

subroutine init_mcwat_mask(atm_cnt, nres, igraph, isymbl, res_atms, labres, crd)

    use mdin_ctrl_dat_mod
    use file_io_dat_mod
    use findmask_mod

    implicit none

    ! Formal arguments:
    integer, intent(in)             :: atm_cnt, nres
    integer, intent(in)             :: res_atms(nres)
    character(len=4), intent(in)    :: igraph(atm_cnt)
    character(len=4), intent(in)    :: isymbl(atm_cnt)
    character(len=4), intent(in)    :: labres(nres)
    double precision, intent(in)    :: crd(3 * atm_cnt)

    if(.not. allocated(mcwatregionmask)) allocate(mcwatregionmask(atm_cnt))

    if(mcwatmask .ne. '') then
      call atommask(atm_cnt, nres, 0, igraph, isymbl, res_atms, labres, &
             crd, mcwatmask, mcwatregionmask)
    end if

end subroutine init_mcwat_mask


!*******************************************************************************
!
! Subroutine:  get_atom_fract_crd
!
! Description: Gets the fractional crd for 1 atom 
!
!*******************************************************************************

subroutine get_atom_fract_crd(X,Y,Z,fracx,fracy,fracz)

  use pbc_mod
  implicit none

  double precision, intent(in)          :: X,Y,Z
  double precision, intent(out)         :: fracx,fracy,fracz
  double precision                      :: f1,f2,f3
  if (is_orthog .ne. 0) then
    f1 = recip(1,1) * X 
    f2 = recip(2,2) * Y
    f3 = recip(3,3) * Z
    fracx = f1 - dnint(f1)
    fracy = f2 - dnint(f2)
    fracz = f3 - dnint(f3)
  else
    f1 = X * recip(1, 1) + Y * recip(2, 1) + &
         Z * recip(3, 1)
    f2 = X * recip(1, 2) + Y * recip(2, 2) + &
         Z * recip(3, 2)
    f3 = X * recip(1, 3) + Y * recip(2, 3) + &
         Z * recip(3, 3)
    fracx = f1 - dnint(f1)
    fracy = f2 - dnint(f2)
    fracz = f3 - dnint(f3)
  end if

  ! We must have fractional coordinates in the range of 0.0 - 0.999...
  ! The algorithm used above will produce fractionals in the range of 0.0 -
  ! 1.0, with some possibility of slight underflow (neg value) due to
  ! rounding error (confirmed). SO we force fractionals to be nonredundant
  ! and within the anticipated range here.
  if (fracx .lt. 0.d0) fracx = fracx + 1.d0
  if (fracx .ge. 1.d0) fracx = fracx - 1.d0
  if (fracy .lt. 0.d0) fracy = fracy + 1.d0
  if (fracy .ge. 1.d0) fracy = fracy - 1.d0
  if (fracz .lt. 0.d0) fracz = fracz + 1.d0
  if (fracz .ge. 1.d0) fracz = fracz - 1.d0

end subroutine get_atom_fract_crd

!*******************************************************************************
!
! Subroutine:  determineAtomRadius
!
! Description: This function determines the atom radius for each atom type,
! as default takes the smallest radius. It then increments the radius by the 
! effective water radius, to define the size of the spherical region this atom is
! considered to block on the steric grid.
! 
!*******************************************************************************

subroutine determineAtomRadius(aid,totalRadius)

    use prmtop_dat_mod
    use mdin_ctrl_dat_mod, only : clambda
    use ti_mod

    implicit none

    integer, intent(in)           :: aid
    double precision, intent(out) :: totalRadius 

! local variables 
    integer                       :: round  
    double precision              :: r  

    r = 0 
    round = nint(atm_mass(aid))
    if (round > 3) then 
        select case (round)
            case (7)                ! Li
                r = 1.82     
            case (12)               ! C 
                r = 1.70
            case (14)               ! N  
                r = 1.55
            case (16)               ! O 
                r = 1.52
            case (19)               ! F 
                r = 1.47
            case (23)               ! Na
                r = 2.27
            case (24)               ! Mg
                r = 1.73     
            case (31)               ! P
                r = 1.95
            case (32)               ! S
                r = 1.80
            case (35)               ! Cl            
                r = 1.75
            case (36)               ! Cl           
                r = 1.75
            case (39)               ! K 
                r = 2.75    
            case (40)               ! Ca
                r = 2.31 
            case (56)               ! Fe              
                r = 1.26    
            case (59)               ! Ni / Co 
                r = 1.63 
            case (63)               ! Cu               
                r = 1.40     
            case (64)               ! Cu                
                r = 1.40     
            case (65)               ! Zn  
                r = 1.39
            case (66)               ! Zn  
                r = 1.39
            case (79)               ! Br 
                r = 1.85
            case (106)              ! Pd   
                r = 1.63
            case (107)              ! Pd   
                r = 1.63
            case (127)              ! I                
                r = 1.98    
            case (195)              ! Pt
                r = 1.75
            case (200)              ! Hg
                r = 1.55    
            case (201)              ! Hg               
                r = 1.55    
            case (207)              ! Pb
                r = 2.02
            case default
                r = 1.39 ! set to be the lowest radius
            end select
        end if 

        if (r .gt. 0) then 
            totalRadius = r + WATERRADIUS
        else 
            totalRadius = 0 
        endif

        if (ti_mode .ge. 1) then
          if (ti_lst(2,aid) .eq. 1 .and. clambda .lt. 0.5) then
            totalRadius=0
          endif
          if (ti_lst(1,aid) .eq. 1 .and. clambda .ge. 0.5) then
            totalRadius=0
          endif
        endif
!       print *,'r, totalRadius = ', r, totalRadius

end subroutine determineAtomRadius


!*******************************************************************************
!
! Subroutine:  buildGrid
!
! Description: This function builds the entire grid for the first time using the function 
! "placeInGrid". For each atom the entire volume of the atom is incremented. 
! Nearby atoms can overlap in the grid volumes, resulting in values higher than 1.  
! This information is used when updating the grid as waters are moved.
!
!*******************************************************************************

subroutine buildGrid(fraction,stericGrid,atm_cnt,mxX,mxY,mxZ,mxXVxl,mxYVxl,mxZVxl)

    use prmtop_dat_mod

    implicit none

    integer, intent(in)           :: atm_cnt, mxXVxl,mxYVxl,mxZVxl
    double precision, intent(in)  :: fraction(3, atm_cnt)
    integer, intent(inout)        :: stericGrid(:,:,:)
    double precision, intent(in)  :: mxX, mxY, mxZ

    integer       :: counter,mass, l
    integer       :: i,ic
    double precision       :: X,Y,Z

    do i = 1, atm_cnt
        mass = nint(atm_mass(i))
        if (mass > 3) then
            X = fraction(1, i) - 0.5 ! atomic units
            Y = fraction(2, i) - 0.5
            Z = fraction(3, i) - 0.5
            if(X .lt. 0.0d0) then
                X=X+1.0d0
            end if
            if(Y .lt. 0.0d0) then
                Y=Y+1.0d0
            end if
            if(Z .lt. 0.0d0) then
                Z=Z+1.0d0
            end if
            call placeInGrid(stericGrid,mxX,mxY,mxZ,mxXVxl,mxYVxl,mxZVxl,X,Y,Z,i)
        endif
    end do

end subroutine buildGrid

!*******************************************************************************
!
! Subroutine:  placeInGrid
!
! Description: This function marks the volume of an atom on the grid as sterically
! blocked or unblocked. If the variable "del" (delta) is 1, the grid is incremented. if
! "del" is -1, it is decremented. A voxel is unblocked only when it is set to 0
!
!*******************************************************************************

subroutine placeInGrid(stericGrid,mxX,mxY,mxZ,mxXVxl,mxYVxl,mxZVxl,X,Y,Z,aid)
    implicit none

    integer, intent(inout)        :: stericGrid(:,:,:)
    double precision, intent(in)  :: mxX,mxY,mxZ
    double precision, intent(inout)  :: X,Y,Z
    integer, intent(in)           :: mxXVxl, mxYVxl, mxZVxl, aid

! local variables : 
    integer                       :: p, round, l, m, n  
    double precision              :: r, totalRadius  
    integer              :: xVxl, yVxl, zVxl
    integer              :: atmMnXVxl,atmMxXVxl,atmMnYVxl, atmMxYVxl,atmMnZVxl,atmMxZVxl
    double precision              :: xRange,yRange,Zrange

    totalRadius=steric_radius(aid)
    xVxl=X*mxXVxl+1
    yVxl=Y*mxYVxl+1
    zVxl=Z*mxZVxl+1

    zRange = totalRadius/GSP
    atmMnZVxl = int(zVxl - zRange) + 1  
    atmMxZVxl = int(zVxl + zRange)  
    do n = max(atmMnZVxl,1), min(atmMxZVxl, mxZVxl)
        yRange = sqrt(zRange*zRange - (zVxl - n)**2)

        atmMnYVxl = int(yVxl - yRange) +1 
        atmMxYVxl = int(yVxl + yRange)  

        do m = max(atmMnYVxl,1), min(atmMxYVxl, mxYVxl)

            xRange = sqrt(yRange*yRange - (yVxl - m)**2)

            if (xRange .gt. 0.0d0) then  
       
                atmMnXVxl = int(xVxl - xRange) + 1 
                atmMxXVxl = int(xVxl + xRange)  
                do l = max(atmMnXVxl,1), min(atmMxXVxl, mxXVxl)
                    stericGrid(l,m,n) = 1

                end do

            end if            

        end do

    end do
end subroutine placeInGrid

!*******************************************************************************
!
! Subroutine:  ligandCenterOfMass
!
! Description: This function computes mask CoM
!
!*******************************************************************************

subroutine ligandCenterOfMass(atm_cnt, mass, crd, mask, mX, mY, mZ, com_voxel)

    implicit none
    integer, intent(in)           :: atm_cnt
    double precision, intent(in)           :: mass(atm_cnt)
    double precision, intent(in)           :: crd(3, atm_cnt)
    integer, intent(in)           :: mask(atm_cnt)
    integer, intent(in)           :: mX, mY, mZ
    integer, intent(out)          :: com_voxel(3)
    double precision              :: com(3)
    double precision              :: fracx, fracy, fracz
    integer :: i
    double precision :: mass_total
    com(:)=0.0d0
    do i = 1, atm_cnt
        if(mask(i) .eq. 1) then
            com(1) = crd(1,i) * mass(i) + com(1)
            com(2) = crd(2,i) * mass(i) + com(2)
            com(3) = crd(3,i) * mass(i) + com(3)
            mass_total = mass(i) + mass_total
        end if
    end do
    com(1) = com(1)/mass_total
    com(2) = com(2)/mass_total
    com(3) = com(3)/mass_total
    call get_atom_fract_crd(com(1),com(2),com(3),fracx,fracy,fracz)
    com_voxel(1)=int(fracx*mX) + 1
    com_voxel(2)=int(fracy*mY) + 1
    com_voxel(3)=int(fracz*mZ) + 1
end subroutine

!*******************************************************************************
!
! Subroutine:  determineNumberOfEmptyVoxels
!
! Description: This function is actually not necessary. it it just if you want
! to know the number of empty voxels
!
!*******************************************************************************

subroutine determineNumberOfEmptyVoxels(stericGrid,mxXVxl,mxYVxl,mxZVxl,nmEmpVxls, voxelOffset)
    use mdin_ctrl_dat_mod
    implicit none
    integer, intent(in)           :: stericGrid(:,:,:)
    integer, intent(in)           :: mxXVxl, mxYVxl, mxZVxl  
    integer, intent(out)          :: nmEmpVxls
    integer                       :: voxelOffset(6)

 ! local variables
    integer                       :: l,m,n,counter 
    logical                       :: xbool, ybool, zbool

    counter = 0
    do n = 1, mxZVxl
        do m = 1, mxYVxl
            do l = 1, mxXVxl
              xbool=.false.
              ybool=.false.
              zbool=.false.
              if(voxelOffset(1) .gt. voxelOffset(2)) then
                if(l .le. voxelOffset(1) .or. l .ge. voxelOffset(2)) then
                  xbool=.true.
                end if
              else
                if(l .ge. voxelOffset(1) .and. l .le. voxelOffset(2)) then
                  xbool=.true.
                end if
              end if
              if(voxelOffset(3) .gt. voxelOffset(4)) then
                if(m .le. voxelOffset(3) .or. m .ge. voxelOffset(4)) then
                  ybool=.true.
                end if
              else
                if(m .ge. voxelOffset(3) .and. m .le. voxelOffset(4)) then
                  ybool=.true.
                end if
              end if
              if(voxelOffset(5) .gt. voxelOffset(6)) then
                if(n .le. voxelOffset(5) .or. n .ge. voxelOffset(6)) then
                  zbool=.true.
                end if
              else
                if(n .ge. voxelOffset(5) .and. n .le. voxelOffset(6)) then
                  zbool=.true.
                 end if
              end if
              if (xbool .and. ybool .and. zbool .and. stericGrid(l,m,n) .eq. 0) then
                    counter = counter + 1
              end if 
            end do
        end do
    end do
    nmEmpVxls = counter

end subroutine determineNumberOfEmptyVoxels


subroutine updateEmptyVoxelsList(stericGrid,mxXVxl,mxYVxl,mxZVxl,nmEmpVxls,listEmptyVoxels, voxelOffset)
    use mdin_ctrl_dat_mod

    implicit none
    integer, intent(in)           :: stericGrid(:,:,:)
    integer, intent(in)           :: mxXVxl, mxYVxl, mxZVxl
    integer, intent(out)          :: nmEmpVxls
    integer, intent(inout) :: listEmptyVoxels(:,:)
    integer                       :: voxelOffset(6)

 ! local variables
    integer                       :: l,m,n, counter, mnXVxl, mnYVxl, mnZVxl
    logical                       :: xbool, ybool, zbool
    counter = 0
    do n = 1, mxZVxl
        do m = 1, mxYVxl
            do l = 1, mxXVxl
              xbool=.false.
              ybool=.false.
              zbool=.false.
              if(voxelOffset(1) .gt. voxelOffset(2)) then
                if(l .le. voxelOffset(1) .or. l .ge. voxelOffset(2)) then
                  xbool=.true.
                end if
              else
                if(l .ge. voxelOffset(1) .and. l .le. voxelOffset(2)) then
                  xbool=.true.
                end if
              end if
              if(voxelOffset(3) .gt. voxelOffset(4)) then
                if(m .le. voxelOffset(3) .or. m .ge. voxelOffset(4)) then
                  ybool=.true.
                end if
              else
                if(m .ge. voxelOffset(3) .and. m .le. voxelOffset(4)) then
                  ybool=.true.
                end if
              end if
              if(voxelOffset(5) .gt. voxelOffset(6)) then
                if(n .le. voxelOffset(5) .or. n .ge. voxelOffset(6)) then
                  zbool=.true.
                end if
              else
                if(n .ge. voxelOffset(5) .and. n .le. voxelOffset(6)) then
                  zbool=.true.
                 end if
              end if
              if (xbool .and. ybool .and. zbool .and. stericGrid(l,m,n) .eq. 0) then
                 counter = counter + 1
                 listEmptyVoxels(counter,1) = l
                 listEmptyVoxels(counter,2) = m
                 listEmptyVoxels(counter,3) = n
              end if
            end do
        end do
    end do

    nmEmpVxls = counter

end subroutine updateEmptyVoxelsList


!*******************************************************************************
!
! Subroutine:  returnRandomVoxel
!
! Description: The function returns a random voxel that will be tested as a possible
! destination of the water molecule to be moved.
!
!*******************************************************************************

#ifdef CUDA
subroutine returnRandomVoxel(nmEmpVxls,listEmptyVoxels,xdim,ydim,zdim,rndmVxl)
#else
subroutine returnRandomVoxel(nmEmpVxls,listEmptyVoxels,rndmVxl)
#endif

    use pmemd_lib_mod
    use random_mod

    implicit none

    integer, intent(in)          :: nmEmpVxls
#ifdef CUDA
    integer, intent(in)          :: listEmptyVoxels(:)
    integer, intent(in)          :: xdim, ydim, zdim
#else
    integer, intent(in)          :: listEmptyVoxels(:,:)
#endif
    integer, intent(out)          :: rndmVxl(:)

    !local variable
    double precision              :: temp 
    integer                       :: tempVxl

    call amrand_gen(mcres_randgen, temp)

#ifdef CUDA
    tempVxl = listEmptyVoxels(ceiling(temp * nmEmpVxls))
    rndmVxl(1) = tempVxl/(xdim*ydim)+1
    tempVxl = tempVxl - ((rndmVxl(1)-1)*xdim*ydim)
    rndmVxl(2) = tempVxl/xdim+1
    rndmVxl(3) = mod(tempVxl,xdim)+1
#else
    tempVxl = ceiling(temp * nmEmpVxls)
    rndmVxl(1) = listEmptyVoxels(tempVxl,1)
    rndmVxl(2) = listEmptyVoxels(tempVxl,2)
    rndmVxl(3) = listEmptyVoxels(tempVxl,3)
#endif

end subroutine returnRandomVoxel

!*******************************************************************************
!
! Subroutine:  returnRandomWaterMolecule
!
! Description: the function returns a random water molecule for a move attempt.
!
!*******************************************************************************

subroutine returnRandomWaterMolecule(rndmWtr, fraction, mol_str,atm_cnt,mxXVxl,mxYVxl,mxZVxl, &
                                     voxelOffset)

    use pmemd_lib_mod 
    use mdin_ctrl_dat_mod
    use file_io_dat_mod
    use prmtop_dat_mod
    use mol_list_mod

    implicit none  

    integer, intent(in)           :: atm_cnt, mxXVxl,mxYVxl,mxZVxl
    double precision, intent(in)  :: fraction(3,atm_cnt)
    integer, intent(out)          :: rndmWtr
    character(4), intent(in)      :: mol_str
    integer                       :: voxelOffset(6)
! local variables 
    logical                       :: yesWat
    double precision              :: temp
    double precision              :: X,Y,Z
    integer                       :: vxlx,vxly,vxlz
    integer                       :: attempt
    logical                       :: xbool, ybool, zbool
    attempt=1
    yesWat = .false. 
    do while (yesWat .eqv. .false.)

        call amrand_gen(mcres_randgen, temp)
        rndmWtr = ceiling(temp*nres) 

        if(gbl_labres(rndmWtr) .eq. mol_str)  then
            X=fraction(1,gbl_res_atms(rndmWtr))-0.5
            Y=fraction(2,gbl_res_atms(rndmWtr))-0.5
            Z=fraction(3,gbl_res_atms(rndmWtr))-0.5
            if(X .lt. 0.0d0) then
                X=X+1.0d0
            end if
            if(Y .lt. 0.0d0) then
                Y=Y+1.0d0
            end if
            if(Z .lt. 0.0d0) then
                Z=Z+1.0d0
            end if
            vxlx=X*mxXVxl
            vxly=Y*mxYVxl
            vxlz=Z*mxZVxl
            xbool=.false.
            ybool=.false.
            zbool=.false.
            if(voxelOffset(1) .gt. voxelOffset(2)) then
               if(vxlx .le. voxelOffset(1) .or. vxlx .ge. voxelOffset(2)) then
                   xbool=.true.
               end if
            else
               if(vxlx .ge. voxelOffset(1) .and. vxlx .le. voxelOffset(2)) then
                   xbool=.true.
               end if
            end if
            if(voxelOffset(3) .gt. voxelOffset(4)) then
               if(vxly .le. voxelOffset(3) .or. vxly .ge. voxelOffset(4)) then
                   ybool=.true.
               end if
            else
               if(vxly .ge. voxelOffset(3) .and. vxly .le. voxelOffset(4)) then
                   ybool=.true.
               end if
            end if
            if(voxelOffset(5) .gt. voxelOffset(6)) then
               if(vxlz .le. voxelOffset(5) .or. vxlz .ge. voxelOffset(6)) then
                   zbool=.true.
               end if
            else
               if(vxlz .ge. voxelOffset(5) .and. vxlz .le. voxelOffset(6)) then
                   zbool=.true.
               end if
            end if
            if(xbool .and. ybool .and. zbool) then
               yesWat = .true.
            end if
        end if
        if(attempt .gt. mcwatretry) then
            rndmWtr=-1
            exit
        end if
        attempt=attempt+1
    end do 

end subroutine returnRandomWaterMolecule


subroutine buildCoarseGrid(fraction,coarseGrid,atm_cnt,coarseGSP, &
                           maxCoarseXVxl,maxCoarseYVxl,maxCoarseZVxl,mxX,mxY,mxZ,maxatoms)

   use mdin_ctrl_dat_mod
   use pbc_mod
   use ti_mod 

   implicit none

   integer, intent(in)                       :: atm_cnt
   double precision, intent(in)              :: fraction(3, atm_cnt)
   double precision, intent(in)              :: mxX,mxY,mxZ, coarseGSP
   integer, intent(inout)                    :: coarseGrid(:,:,:,:)
   integer, intent(in)                       :: maxCoarseXVxl,maxCoarseYVxl,maxCoarseZVxl
   integer, intent(out)                      :: maxatoms

!local variables
   double precision       :: X,Y,Z
   double precision  :: fx,fy,fz
   integer           :: tempX,tempY,tempZ, i
   integer   :: maxbktx,maxbkty,maxbktz
   maxatoms=0
   do i = 1, atm_cnt
       !if (ti_mode .ge. 1) then
       !     if (ti_lst(2,i) .eq. 1 .and. clambda .lt. 0.5) then
       !         cycle
       !     endif
       !     if (ti_lst(1,i) .eq. 1 .and. clambda .ge. 0.5) then
       !         cycle
       !     endif
       !endif
       fx = fraction(1,i) - 0.5d0
       fy = fraction(2,i) - 0.5d0
       fz = fraction(3,i) - 0.5d0
       if(fx .lt. 0.0d0) then
           fx = fx + 1.0d0
       end if
       if(fy .lt. 0.0d0) then
           fy = fy + 1.0d0
       end if
       if(fz .lt. 0.0d0) then
           fz = fz + 1.0d0
       end if
       tempX = fx*maxCoarseXVxl+1
       tempY = fy*maxCoarseYVxl+1
       tempZ = fz*maxCoarseZVxl+1
       coarseGrid(tempX,tempY,tempZ,1) = (coarseGrid(tempX,tempY,tempZ,1) + 1)
       if(coarseGrid(tempX,tempY,tempZ,1) .gt. maxatoms) then
           maxatoms=coarseGrid(tempX,tempY,tempZ,1)
           maxbktx=tempX-1
           maxbkty=tempY-1
           maxbktz=tempZ-1
       end if
       coarseGrid(tempX,tempY,tempZ,(coarseGrid(tempX,tempY,tempZ,1)+1)) = i
    end do
    ! Add one to max atoms since this is used for allocation and the first index is the size of the bucket
    maxatoms=maxatoms+1

end subroutine buildCoarseGrid


subroutine calculateCoarseEnergy(atm_cnt, mol_id, crd, fraction, shift, CoarseGrid, coarseGSP, mX, mY, mZ, &
                                 mxX, mxY, mxZ, coarseene)

    use mol_list_mod
    use pbc_mod
    use prmtop_dat_mod
    use mdin_ctrl_dat_mod
    use prmtop_dat_mod, only : ntypes, typ_ico, gbl_cn1, gbl_cn2, atm_qterm
    use ene_frc_splines_mod
    use ti_mod
 
    implicit none

    integer, intent(in)           :: atm_cnt, mol_id
    double precision, intent(in)  :: crd(3,atm_cnt)
    double precision, intent(in)  :: fraction(3,atm_cnt)
    double precision, intent(inout)  :: shift(3)
    integer, intent(inout)        :: coarseGrid(:,:,:,:)
    double precision, intent(in)  :: coarseGSP
    integer, intent(in)           :: mX, mY, mZ ! max of coarse grid
    double precision, intent(in)  :: mxX, mxY, mxZ ! max of pbc box
    double precision              :: coarseene

    integer                       :: grid(3)
    integer                       :: i,j,k, res_offset, cur_wat_atom
    integer                       :: bkt_x,bkt_y,bkt_z, bkt_offset,bkt_atms
    integer                       :: wrap_bkt_x, wrap_bkt_y, wrap_bkt_z
    integer                       :: iaci, ic
    double precision              :: crd_i(3)
    double precision              :: atm_i(3)
    double precision              :: atm_j(3)
    double precision              :: cn1, cn2
    double precision              :: dx,dy,dz, delx, dely, delz, delr2
    double precision              :: charge_i, charge_j
    double precision              :: delr2inv, r6, f12, f6
    double precision              :: elec_ene, lj_ene 
    double precision              :: dens_efs, del_efs, du, du2, du3, b0
    integer                       :: ind
    double precision              :: ele_weight, vdw_weight

    elec_ene = 0.0d0
    lj_ene = 0.0d0

    bkt_offset = 3

    res_offset = gbl_res_atms(mol_id)
    coarseene = 0.0
    do j = res_offset, res_offset + mcres_numatms-1
       crd_i(1) = crd(1,j) - shift(1)
       crd_i(2) = crd(2,j) - shift(2)
       crd_i(3) = crd(3,j) - shift(3)
       call get_atom_fract_crd(crd_i(1),crd_i(2),crd_i(3),atm_i(1),atm_i(2),atm_i(3))
       grid(1) = int(atm_i(1)*mX) + 1
       grid(2) = int(atm_i(2)*mY) + 1
       grid(3) = int(atm_i(3)*mZ) + 1

       iaci = ntypes * (atm_iac(j) - 1)

       do bkt_x=-bkt_offset,bkt_offset
           wrap_bkt_x = modulo(grid(1)+bkt_x-1,mX)+1
          
           do bkt_y=-bkt_offset,bkt_offset
               wrap_bkt_y = modulo(grid(2)+bkt_y-1,mY)+1
             
               do bkt_z=-bkt_offset,bkt_offset
                   wrap_bkt_z = modulo(grid(3)+bkt_z-1,mZ)+1
                   bkt_atms = coarseGrid(wrap_bkt_x,wrap_bkt_y,wrap_bkt_z,1)
                   do i=1,bkt_atms  ! We ignore a distance check here and do all calcs

                        k = coarseGrid(wrap_bkt_x, wrap_bkt_y, wrap_bkt_z, i+1)

                        cur_wat_atom=j-res_offset
                        if(k-j+cur_wat_atom .lt. mcres_numatms .and. k-j+cur_wat_atom .ge. 0) cycle

                        atm_j(1) = fraction(1,k) - 0.50d0
                        atm_j(2) = fraction(2,k) - 0.50d0 
                        atm_j(3) = fraction(3,k) - 0.50d0

                        delx = (atm_i(1) - atm_j(1)) - nint(atm_i(1)-atm_j(1))
                        dely = (atm_i(2) - atm_j(2)) - nint(atm_i(2)-atm_j(2))
                        delz = (atm_i(3) - atm_j(3)) - nint(atm_i(3)-atm_j(3))
                        
                        dx = delx*ucell(1,1) + dely*ucell(1,2) + delz*ucell(1,3)
                        dy = dely*ucell(2,2) + delz*ucell(2,3)
                        dz = delz*ucell(3,3)

                        delr2 = dx * dx + dy * dy + dz * dz
                        ic = typ_ico(iaci + atm_iac(k))

                        if(delr2 .le. vdw_cutoff*vdw_cutoff) then
                            ele_weight=1.0
                            vdw_weight=1.0
                            if(ti_mode .ge. 1) then
                              if(ti_lst(1,k) .gt. 0) then
                                ele_weight=ti_item_weights(4, 1)
                                vdw_weight=ti_item_weights(6, 1)
                              end if
                              if(ti_lst(2,k) .gt. 0) then
                                ele_weight=ti_item_weights(4, 2)
                                vdw_weight=ti_item_weights(6, 2)
                              end if
                            end if
                            cn1 = gbl_cn1(ic)
                            cn2 = gbl_cn2(ic)
                            delr2inv = 1.0d0/delr2
                            r6 = delr2inv * delr2inv *delr2inv
                            f6 = cn2 * r6
                            f12 = cn1 * r6 * r6
                            lj_ene = lj_ene + (f12 - f6) * vdw_weight

                            charge_i = atm_qterm(j)
                            charge_j = atm_qterm(k)
                            dens_efs = efs_tbl_dens
                            del_efs =  1.0d0/dens_efs
                            ind = int(dens_efs * delr2)
                            du = delr2 - dble(ind) * del_efs
                            du2 = du * du
                            du3 = du * du2
                            ind = ishft(ind, 3)             ! 8 * ind

                            elec_ene = elec_ene + charge_i*charge_j * ele_weight * (efs_tbl(1 + ind) + du * efs_tbl(2 + ind) + &
                                 du2 * efs_tbl(3 + ind) + du3 * efs_tbl(4 + ind))
                        end if

                    end do
                end do
            end do 
        end do
     end do

     coarseene = lj_ene + elec_ene
end subroutine calculateCoarseEnergy


!*******************************************************************************
!
! Subroutine:  mcres
!
! Description: This is the main subroutine, which orchestrates
! the other subroutines for the MC function.
! This routine uses the namelist variables mcwat and mcint. 
! If mcwat is > 0 this routine gets called.
! It is also only called if nstep routine it proceeds as follows:
!
!              1) Calculate initial energy (via call to pme_force) for current coordinates.
!              2) Check if residue to be move mcres (currently hardwired to 'WAT '
!                 exists.
!              3) Picks a random residue within the set of mcres.
!              4) Sets up a translate vector and moves the residue.
!              5) Calculates new energy.
!              6) Applies standard Metropolis criteria to accept or reject move.
!              7) If accepted it quits. If rejected it uses the negative of the
!                 translate vector to translate back to original state.
!
!*******************************************************************************

subroutine mcres(accept,mol_str, atm_cnt, crd, frc, &
                img_atm_map, atm_img_map, new_list, &
                need_pot_enes, need_virials, pot_ene, nstep, virial, &
                ekcmt, pme_err_est)

    use mol_list_mod
    use mdin_ctrl_dat_mod
    use pbc_mod
    use gbl_constants_mod, only : KB
    use energy_records_mod, only : pme_pot_ene_rec 
    use file_io_dat_mod
    use pme_force_mod
    use random_mod, only : random_state, amrand_gen, amrset_gen
    use prmtop_dat_mod
    use pmemd_lib_mod 
    use pbc_mod
    use pme_force_mod
    use file_io_dat_mod
    use pmemd_lib_mod
    use runfiles_mod
    use timers_mod
    use state_info_mod
    use bintraj_mod

    implicit none

    character(4)                  :: mol_str
    integer                       :: atm_cnt
    double precision              :: crd(3, atm_cnt)
    double precision              :: frc(3, atm_cnt)
    integer                       :: img_atm_map(atm_cnt)
    integer                       :: atm_img_map(atm_cnt)
    logical                       :: new_list
    logical                       :: need_pot_enes
    logical                       :: need_virials
    type(pme_pot_ene_rec)         :: pot_ene
    double precision              :: virial(3)            ! Only used for MD
    double precision              :: ekcmt(3)             ! Only used for MD
    double precision              :: pme_err_est          ! Only used for MD
    integer, intent(in)           :: nstep

    !local variables
    double precision translate(3)
    integer mol_id, res_offset, i,j,k, atm_id, res_id
    double precision random, rx, ry, rz
    double precision acceptance_criteria
    logical accept
    double precision old_ene, new_ene
    integer cyc_cnt, cnt
    integer mol_atm_cnt
    real current_time
    double precision :: coarseene    
    integer total_nstep
    integer l, rndmWtr
    integer full_energy_cnt
    integer maxCoarseAtoms
    double precision      :: fraction(3, atm_cnt)
 
!!! for the coarse energy
    integer              :: mxXVxl,mxYVxl,mxZVxl 
    integer, allocatable :: waterList(:)
    double precision     :: deltafrac(3)
    integer, allocatable :: stericGrid(:,:,:)  
    integer, allocatable :: coarseGrid(:,:,:,:)
    double precision                 :: coarse_before_ene, coarse_after_ene, delta_coarse
    double precision     :: boltz_sum, boltz_tracker
    integer              :: accept_idx
    integer              :: rndmWatAtomArray(mcrescyc)
    integer              :: rndmWatArray(mcrescyc)
    double precision     :: translate_array(3,mcrescyc)
    double precision     :: acceptance_array(mcrescyc)
    double precision     :: full_energies(mcrescyc+1)
    integer              :: full_to_coarse(mcrescyc+1)
#ifdef CUDA
    integer, allocatable :: listEmptyVoxels(:)  
#else
    integer, allocatable :: listEmptyVoxels(:,:)  
#endif
    double precision     :: save_crd(3,mcres_numatms)
    double precision     :: vel(3, atm_cnt)
    integer              :: ligoffset
    integer              :: com_voxel(3)
    integer              :: voxelOffset(6)

    !Outer loop over number of MC moves to try.
    !It will quit this loop if it gets an accepted move.

#ifdef CUDA
    ! Crds are only updated when nstep%ntwx=0.  This ensures the CPU has the latest
    ! coordinate copy for molecule movement.
    call gpu_download_crd(crd)
    !This call fixes translation vectors going to NaN
    call gpu_download_ucell(ucell)
#endif
#ifndef MPI
    call pme_force(atm_cnt, crd, frc, img_atm_map, atm_img_map, &
                         new_list, .true., .false., &
                         pot_ene, nstep, virial, ekcmt, pme_err_est)
#endif
    ! Store old energy
    old_ene = pot_ene%total

    call get_fract_crds(atm_cnt, crd, fraction)

    coarseGSP = (es_cutoff/2)
    mxXVxl = ceiling(pbc_box(1) / GSP)
    mxYVxl = ceiling(pbc_box(2) / GSP)
    mxZVxl = ceiling(pbc_box(3) / GSP)    ! exceeding the box into the solvent in one axis

    if (.not. allocated(steric_radius)) then
        allocate(steric_radius(atm_cnt))
        do i=1,atm_cnt
            call determineAtomRadius(i,steric_radius(i))
        end do 
    end if
    if(mcligshift .gt. 0.0) then
        call ligandCenterOfMass(atm_cnt, atm_mass, crd, mcwatregionmask, mxXVxl, mxYVxl, mxZVxl, com_voxel)
        ligoffset=ceiling(mcligshift/GSP)
        voxelOffset(1)=com_voxel(1)-ligoffset
        voxelOffset(2)=com_voxel(1)+ligoffset
        voxelOffset(3)=com_voxel(2)-ligoffset
        voxelOffset(4)=com_voxel(2)+ligoffset
        voxelOffset(5)=com_voxel(3)-ligoffset
        voxelOffset(6)=com_voxel(3)+ligoffset
        if(voxelOffset(1) .lt. 0) then
            voxelOffset(1) = mxXVxl - voxelOffset(1)
        end if
        if(voxelOffset(2) .gt. mxXVxl) then
            voxelOffset(2) = voxelOffset(2) - mxXVxl
        end if
        if(voxelOffset(3) .lt. 0) then
            voxelOffset(3) = mxYVxl - voxelOffset(3)
        end if
        if(voxelOffset(4) .gt. mxYVxl) then
            voxelOffset(4) = voxelOffset(4) - mxYVxl
        end if
        if(voxelOffset(5) .lt. 0) then
            voxelOffset(5) = mxZVxl - voxelOffset(5)
        end if
        if(voxelOffset(6) .gt. mxZVxl) then
            voxelOffset(6) = voxelOffset(6) - mxZVxl
        end if
    else
        voxelOffset(1)=1
        voxelOffset(2)=mxXVxl
        voxelOffset(3)=1
        voxelOffset(4)=mxYVxl
        voxelOffset(5)=1
        voxelOffset(6)=mxZVxl
    end if
#ifdef CUDA
    call gpu_steric_grid(atm_cnt, mxXVxl, mxYVxl, mxZVxl, GSP, steric_radius, nmEmpVxls, voxelOffset)
    allocate(listEmptyVoxels(int(nmEmpVxls*1.2)),stat=error)
    if (stat.ne.0) then
        write(0,*)"error: couldnt allocate memory for array, listEmptyVoxels=",stat
        stop 1
    endif
    call gpu_steric_grid_empty_voxels(nmEmpVxls,listEmptyVoxels)
#else
    allocate(stericGrid(mxXVxl, mxYVxl, mxZVxl),stat=error)

    if (stat.ne.0) then
        write(0,*)"error: couldnt allocate memory for array, g=",stat
        stop 1
    endif
    stericGrid(:,:,:)=0
    call buildGrid(fraction,stericGrid,atm_cnt,pbc_box(1),pbc_box(2),pbc_box(3),mxXVxl,mxYVxl,mxZVxl)

    call determineNumberOfEmptyVoxels(stericGrid,mxXVxl,mxYVxl,mxZVxl,nmEmpVxls,voxelOffset)
    allocate(listEmptyVoxels(int(nmEmpVxls*1.2),3),stat=error)

    if (stat.ne.0) then
        write(0,*)"error: couldnt allocate memory for array, listEmptyVoxels=",stat
        stop 1
    endif

    call updateEmptyVoxelsList(stericGrid,mxXVxl,mxYVxl,mxZVxl,nmEmpVxls,listEmptyVoxels,voxelOffset)
#endif
    maxCoarseXVxl = ceiling(pbc_box(1) /coarseGSP)
    maxCoarseYVxl = ceiling(pbc_box(2) /coarseGSP)
    maxCoarseZVxl = ceiling(pbc_box(3) /coarseGSP)

    allocate(coarseGrid(maxCoarseXVxl+3, maxCoarseYVxl+3, maxCoarseZVxl+3, nint(coarseGSP*15)),stat=error)

    if (stat.ne.0) then
        write(0,*)"error: couldnt allocate memory for array, coarseGrid=",stat
        stop 1
    endif
    coarseGrid(:,:,:,:)=0

    call buildCoarseGrid(fraction,coarseGrid,atm_cnt,coarseGSP, maxCoarseXVxl,maxCoarseYVxl,maxCoarseZVxl, &
                         pbc_box(1),pbc_box(2),pbc_box(3),maxCoarseAtoms)
    call returnRandomWaterMolecule(rndmWtr,fraction,mol_str,atm_cnt,mxXVxl,mxYVxl,mxZVxl, voxelOffset)
    !write(0,*)"Random water",rndmWtr,crd(1,gbl_res_atms(rndmWtr)),crd(2,gbl_res_atms(rndmWtr)),crd(3,gbl_res_atms(rndmWtr))
    !write(0,*)"Voxels",voxelOffset(1),voxelOffset(2),voxelOffset(3),voxelOffset(4),voxelOffset(5),voxelOffset(6)
    if(rndmWtr .eq. -1) then
        write(mdout,'(a/a/a, i8, a)') "Could not find a suitable water.", &
            "Skipping monte carlo water for this iteration.", &
            "Try raising mcwatretry from ", mcwatretry," if this continues."
        goto 9999
    end if
    do cyc_cnt= 1, mcrescyc
#ifdef CUDA
        call returnRandomVoxel(nmEmpVxls,listEmptyVoxels,mxXVxl,mxYVxl,mxZVxl,rndmVxl)
#else
        call returnRandomVoxel(nmEmpVxls,listEmptyVoxels,rndmVxl)
#endif
        call amrand_gen(mcres_randgen, rx)
        call amrand_gen(mcres_randgen, ry)
        call amrand_gen(mcres_randgen, rz)
        rndmWatArray(cyc_cnt)=rndmWtr
        rndmWatAtomArray(cyc_cnt)=gbl_res_atms(rndmWtr)-1 ! Starting atom is -1 on GPU side
        translate_array(1,cyc_cnt)=(rx+rndmVxl(1)-0.5)/float(mxXVxl)
        translate_array(2,cyc_cnt)=(ry+rndmVxl(2)-0.5)/float(mxYVxl)
        translate_array(3,cyc_cnt)=(rz+rndmVxl(3)-0.5)/float(mxZVxl)
        rx = translate_array(1,cyc_cnt)
        ry = translate_array(2,cyc_cnt)
        rz = translate_array(3,cyc_cnt)
        translate_array(1,cyc_cnt) = crd(1,gbl_res_atms(rndmWtr)) - rx*ucell(1,1) - ry*ucell(1,2) - rz*ucell(1,3)
        translate_array(2,cyc_cnt) = crd(2,gbl_res_atms(rndmWtr)) - rx*ucell(2,1) - ry*ucell(2,2) - rz*ucell(2,3)
        translate_array(3,cyc_cnt) = crd(3,gbl_res_atms(rndmWtr)) - rx*ucell(3,1) - ry*ucell(3,2) - rz*ucell(3,3)
    end do
    ! Uploads data, does coarse grid calculations, and outputs an acceptance array
    acceptance_array(:)=0.0
#ifdef CUDA
    call gpu_calc_coarsegrid(coarseGrid, coarseGSP, maxCoarseAtoms, maxCoarseXVxl, maxCoarseYVxl, maxCoarseZVxl, mcrescyc, &
                             rndmWatAtomArray, translate_array, acceptance_array, is_orthog, ti_mode)
#else
    deltafrac(:)=0.0
    do cyc_cnt=1, mcrescyc
        coarseene=0.0
        call calculateCoarseEnergy(atm_cnt, rndmWatArray(cyc_cnt), crd, fraction, deltafrac, coarseGrid, coarseGSP,&
                                   maxCoarseXVxl, maxCoarseYVxl, maxCoarseZVxl, &
                                   pbc_box(1), pbc_box(2), pbc_box(3), coarseene)
        coarse_before_ene = coarseene
        coarseene=0.0
        call calculateCoarseEnergy(atm_cnt, rndmWatArray(cyc_cnt), crd, fraction, translate_array(:,cyc_cnt), &
                                   coarseGrid, coarseGSP,&
                                   maxCoarseXVxl, maxCoarseYVxl, maxCoarseZVxl, &
                                   pbc_box(1),pbc_box(2),pbc_box(3), coarseene)
        coarse_after_ene = coarseene
        acceptance_array(cyc_cnt)=coarse_after_ene-coarse_before_ene
    end do
#endif

    ! First index is the base state of the system
    full_energy_cnt=1
    full_energies(1) = 1.0
    boltz_sum = 1.0

    do cyc_cnt=1, mcrescyc
        if(acceptance_array(cyc_cnt) .lt. mccoarsethresh) then
             full_energy_cnt=full_energy_cnt+1
             cnt = gbl_res_atms(rndmWatArray(cyc_cnt))
             do l=0, mcres_numatms-1
                 save_crd(1,l+1)=crd(1,cnt+l)
                 save_crd(2,l+1)=crd(2,cnt+l)
                 save_crd(3,l+1)=crd(3,cnt+l)
                 crd(1,cnt+l)=crd(1,cnt+l)-translate_array(1,cyc_cnt)
                 crd(2,cnt+l)=crd(2,cnt+l)-translate_array(2,cyc_cnt)
                 crd(3,cnt+l)=crd(3,cnt+l)-translate_array(3,cyc_cnt)
             end do
#ifdef CUDA
             call gpu_force_new_neighborlist()
             call gpu_upload_crd(crd)
#else
             new_list = .true.
#endif
#ifndef MPI
             ! Get new energy
             call pme_force(atm_cnt, crd, frc, img_atm_map, atm_img_map, &
                            .true., .true., .false., &
                            pot_ene, nstep, virial, ekcmt, pme_err_est)
#endif
             new_ene = pot_ene%total
             if(isnan(new_ene) .or. new_ene > 999999999.0d0 .or. new_ene < -999999999.0d0) then
                 new_ene=999999999.0d0
             end if
             !write(0,*)"Ene:",new_ene-old_ene,acceptance_array(cyc_cnt), pot_ene%vdw_tot,pot_ene%elec_tot
             full_energies(full_energy_cnt)=exp(-(new_ene-old_ene)/(temp0*KB))
             full_to_coarse(full_energy_cnt)=cyc_cnt
             boltz_sum=boltz_sum+full_energies(full_energy_cnt)
             ! Revert water position
             do l = 0, mcres_numatms-1 ! return the molecule
                 crd(1,cnt+l) = save_crd(1,1+l)
                 crd(2,cnt+l) = save_crd(2,1+l)
                 crd(3,cnt+l) = save_crd(3,1+l)
             end do
        end if
    end do
    ! Convert everything into probabilities.
    do cyc_cnt=1,full_energy_cnt
        full_energies(cyc_cnt)=full_energies(cyc_cnt)/boltz_sum
    end do
    ! Pick new state out of probabilities.
    call amrand_gen(mcres_randgen, random)
    boltz_tracker=0.0
    accept_idx=0
    do cyc_cnt=1,full_energy_cnt
        boltz_tracker=full_energies(cyc_cnt)+boltz_tracker
        if(random < boltz_tracker) then
            accept_idx=cyc_cnt
            exit
        end if
    end do
    if(accept_idx .ne. 1) then
        write(mdout,*)"Water move success.  Probability of movement: ", 1.0-full_energies(1)
        write(mdout,*)"Water moved: ",rndmWatArray(full_to_coarse(accept_idx))
        write(mdout,*)"Accepted probability: ",full_energies(accept_idx)
        ! Move coordinates
        cyc_cnt=full_to_coarse(accept_idx)
        cnt = gbl_res_atms(rndmWatArray(cyc_cnt))
        do l=0, mcres_numatms-1
            crd(1,cnt+l)=crd(1,cnt+l)-translate_array(1,cyc_cnt)
            crd(2,cnt+l)=crd(2,cnt+l)-translate_array(2,cyc_cnt)
            crd(3,cnt+l)=crd(3,cnt+l)-translate_array(3,cyc_cnt)
        end do
    end if
#ifdef CUDA
    call gpu_force_new_neighborlist()
    call gpu_upload_crd(crd)
#else
    new_list=.true.
#endif

9999 continue

    deallocate(steric_radius,stat=error)
    if (stat.ne.0) then
        write(0,*)"error in deallocating array steric_radius",stat
        stop 1
    endif

    deallocate(listEmptyVoxels,stat=error)
    if (stat.ne.0) then
        write(0,*)"error in deallocating array listEmptyVoxels=",stat
        stop 1
    endif

#ifndef CUDA
    deallocate(stericGrid,stat=error)
    if (error.ne.0) then
        write(0,*)"error in deallocating array stericGrid"
        stop 1
    endif
#endif

    deallocate(coarseGrid,stat=error)
    if (error.ne.0) then
        write(0,*)"error in deallocating array coarseGrid"
        stop 1
    endif

end subroutine mcres

end module mcres_mod
