!=============================================================================================
! This file contains all subroutines needed to construct CV, call CV value or CV force and all
! related utilities. Also a tutorial for how to build a new CV by youself is attached here.
! If you have any questions, please email to sagui@ncsu.edu or fpan3@ncsu.edu.
!=============================================================================================
!**********************************************************************************************
!**       Tutorial -- How to build a new customized Collective Variable in Amber16-PMEMD     **
!**                                                                                          **
!**              by Feng Pan, Christopher Roland and Celeste Sagui ( Jun. 2016 )             **
!**********************************************************************************************

!ABSTRACT
!--------
!This tutorial introduces how to build you own customized Collective Variable (CV) in NFE module 
!of PMEMD in Amber16. After the implementation of the new CV, you can use it for SMD and ABMD 
!simulations included in NFE module. In addition, they are also useable for the PMD and BBMD
!modules for umbrella sampling and H- or T-REMD simulations, which are two modules NOT currently
!released for PMEMD. The PMD and BBMD modules for PMEMD for `as is' usage may be obtained
!by writing to sagui@ncsu.edu and fpan3@ncsu.edu, and will also be made available via a future
!tutorial on the AMBER website

!Current and Future Status of CVs
!-------------------------------
!Current CVs implemented include:

!(i) DISTANCE : distance between two atoms
!(ii) COM_DISTANCE: distance between center of mass of two atom groups
!(iii) DF_COM_DISTANCE: difference of distances between center of mass of first two atom
!      groups and second two atom groups
!(iv) LCOD: linear combination of distances between pairs of atoms
!(v) ANGLE: angle between lines joining atoms
!(vi) COM_ANGLE: angle formed by center of mass of three atom groups
!(vii) TORSION: dihedral angles
!(viii) COM_TORSION: dihedral angle formed by the center of mass of four atom groups
!(ix) COS_OF_DIHEDRAL: sum of cosines of dihedral angles
!(x) SIN_OF_DIHEDRAL: sum of sines of dihedral angles
!(xi) PAIR_DIHEDRAL: sum of cosines of a list of angles
!(xii) PATTERN_DIHEDRAL: a pattern-recognizing function defined on a list of dihedral angles
!(xiii) R_OF_GYRATION: radius of gyration of a group of atoms
!(xiv) MULTI_RMSD: RMS or RMSDs of several groups of atoms with respect to a reference position
!(xv) N_OF_BONDS: cv for hydrogen-bonds
!(xvi) HANDEDNESS: for handedness of different helices
!(xvii) N_OF_STRUCTURES: this counts the number of structures that stay close in the
!       sense with respect to a group of reference structures

!For more details on these CVs, please see the AMBER16 manual.

!We are continuously adding new CVs to this list. To incorporate a newly defined and tested
!CV for use with ABMD and other modules, please contact sagui@ncsu.edu and/or fpan3@ncsu.edu.

!During the coming year, we aim to release a set of quaternion-based CVs for large-scale motions
!of a large number of atoms.

!CV Usage
!--------

!The following are a few examples of using CVs for an ABMD simulation:

!To use CVs in simulation, we need to add the namelist &abmd(for ABMD) or &smd(for SMD) in the mdin
!file, and also put the CV variables in another file to be read.

!For example, if we want to use DISTANCE between two atoms as the CV for SMD and ABMD simulation,
!we can create a new file cv.in with 

! &colvar
!      npath = 2, nharm = 1,
!      path = 8.0, 5.0,
!      harm = 10.0

!      cv_type = 'DISTANCE'
!      cv_min = 3.0, cv_max = 10.0, resolution = 1.0
!      
!      cv_ni = 2, 
!      cv_i = 5, 8,
! /

!To run SMD, we can add the following &smd section to mdin file

! &smd
!   output_file = 'smd.dat'
!   output_freq = 1
!   cv_file = 'cv.in'
! /
! 
!To run ABMD, we can add the following &abmd section to mdin file

! &abmd
!   mode = 'FLOODING'

!   monitor_file = 'monitor.txt'
!   monitor_freq = 1

!   umbrella_file = 'umbrella.nc'

!   timescale = 1.0

!   wt_temperature =  10000.0
!   wt_umbrella_file = 'wt_umbrella.nc'
!   cv_file = 'cv.in'
! /

!Also you can run a more complicated case with two CVs. In this case, you just need to add another
!&colvar section in the cv input file.
!For example, if you have another CV COS_OF_DIHEDRAL, you can write the cv.in file like:

! &colvar
!      npath = 2, nharm = 1,
!      path = 8.0, 5.0,
!      harm = 10.0

!      cv_type = 'DISTANCE'
!      cv_min = 3.0, cv_max = 10.0, resolution = 1.0
!      
!      cv_ni = 2, 
!      cv_i = 5, 8,
! /
! 
! &colvar
!      cv_type = 'COS_OF_DIHEDRAL' ! sum of cosines of dihedral angles

!      cv_ni = 8
!      cv_i =  2,  5,  7, 17,
!              17, 19, 21, 31,

!      path = -5.0, -4.0, -4.0, -5.0
!      npath = 4

!      harm = 100.0
!      nharm = 1
!    
!      cv_min = -2, cv_max = 2, resolution = 0.5,
! /

!If you want to use those CVs in your own codes, it is also not very hard. The basic aspects for
!usage includes the following steps:

!---read CV from a defined input file 
!---calculate the value of CV 
!---calculate the forces base on your own added potential 
!---print out CV information

!For example, if you want to use the CV COM_ANGLE in your code, a possible coding pattern could be:

!subroutine my_method(x,f,cv)
!  
!  use nfe_colvar_mod
!  
!  implicit none
!  
!  type(colvar_t) :: cv
!  double precision :: x(:), f(:)
!  
!  ...
!  call colvar_nlread(cv_unit, cv)
!  ...
!  value = v_COM_ANGLE(cv, x)
!  ...
!  U = ...   ! your defined restraint potential
!  frc = ... ! restraint force
!  call f_COM_ANGLE(cv, x, fcv, f)
!  ...
!  call p_COM_ANGLE(cv, mdout_unit)
!  ...
!end subroutine my_method


!Constructing a new CV
!---------------------
!Here we illustrate the process of implementing a new CV with a simple example. We define the 
!new CV as:

!DISTANCE_Z : projection of a distance vector on an axis.

!This CV has been implemented in NAMD but not in AMBER. Here we just pick the easiest case: the
!distance is between two atoms and the axis is given by a constant vector. It will is not very 
!hard to extend it to more complicated cases (like the distance is between two groups of atoms) 
!once we are familiar of the process of implementation. Perhaps we will finish a generalized
!version of this CV in the future.

!This tutorial consists of the following steps:

!  1.Define the input format of CV.
!  2.Decide the mathematical expression of the value and derivative of the new CV.
!  3.Build new related subroutines inside the source code.
!  4.Compile the new codes and get new executable.
!  5.Test the new CV.


!***************************************************************************************************
!1. Define the input format of CV.

!The first thing we should concern is that how can we read the CV from our input file. The default 
!CV type in the source code has two arrays to read: cv_i and cv_r, so it is straightforward that 
!we have the array cv_i to store the index of atoms and cv_r to store the unit vector of a given axis.

!In this way, an example of CV input file of the new CV could be:

! &colvar
!   cv_type = 'DISTANCE_Z'
!   cv_i = 1, 2 ! the distance between atom 1 and atom 2
!   cv_r = 0, 0, 1 ! calculate the projection to z-axis
! /

!***************************************************************************************************
!2. Decide the mathematical expression of the value and derivative of the new CV.

!Before we write the code, we need to know the mathematical expression to calculate the value of 
!DISTANCE_Z and the derivative with respect to coordinates. Then we can easily put it in Fortran code. 
!For this, the calculation is pretty easy:
!                             value = (r2 - r1) • v0
!                             derivative(r1) = - v0
!                             derivative(r2) = v0
!Here r1,r2 are the coordinates of the two atoms, and v0 is the unit vector of given axis. 

!***************************************************************************************************
!3. Build new related subroutines inside the source code.

!All CV related subroutines are included inside the module nfe_colvar_mod in the file nfe_colvar.F90. 
!So to build a new CV, we need to add certain subroutines in nfe_colvar.F90. First, you need to add 
!the new CV type to the parameters section:

!    integer, public, parameter :: COLVAR_DISTANCE_Z = 20 

!Then basically, four functions or subroutines are needed for a new CV:

!    v_DISTANCE_Z : to calculate the value of the CV
!    f_DISTANCE_Z : to calculate the force by the CV
!    b_DISTANCE_Z : a bootstrap subroutine to check CV before use
!    p_DISTANCE_Z : a print subroutine to print out CV information to mdout

!Here we have written the four subroutines, and attach them to the end of this file.

!(Comments: It is always not hard to write codes to get the value function. For the derivative, if 
!you find it too hard to code the force subroutine, you can use Maple or Matlab to get the fortran code.)

!Now you can freely call the four subprograms from outside to use the new CV. If you want it be well-
!integrated in the NFE module. You need to add the selection of those subroutines to the select case 
!sections in the following subroutines or functions:

!    colvar_value; colvar_force; colvar_bootstrap; colvar_print
!    
!For example, in colvar_value, you should do changes like this:

!  select case(cv%type)
!    case(COLVAR_ANGLE)
!      value = v_ANGLE(cv, x)
!      ......
!    case(COLVAR_DISTANCE_Z)
!      value = v_DISTANCE_Z(cv, x)
!    case default
!      nfe_assert_not_reached()
!      value = dble(0)
!  end select

!Similar changes should be made in the other three subprograms.

!***************************************************************************************************
!4. Compile the new codes and get new executable.

!After finishing the coding part, add new code to nfe_colvar.F90 file. Then recompile PMEMD to get 
!updated executable. For example, to get new pmemd.cuda, we do

!    cd $AMBERHOME && ./configure -cuda < compiler > 
!    cd src/pmemd && make install
!    
!***************************************************************************************************    
!5. Test the new CV.

!To make sure our new CV works well, one should of course do some testing.

!***************************************************************************************************


!***************************************************************************************************
!*********************************** S U B R O U T I N E S *****************************************
!***************************************************************************************************
!function v_DISTANCE_Z(cv, x) result(value)

!   use nfe_lib_mod
!   use parallel_dat_mod

!   implicit none

!   double precision :: value

!   type(colvar_t), intent(in) :: cv

!   double precision, intent(in) :: x(*)

!   integer :: a1, a2

!   nfe_assert(cv%type == COLVAR_DISTANCE_Z)

!   nfe_assert(associated(cv%i))
!   nfe_assert(size(cv%i) == 2)
!   nfe_assert(size(cv%r) == 3)

!   nfe_assert(cv%i(1) > 0 .and. cv%i(1) <= pmemd_natoms())
!   a1 = 3*cv%i(1) - 2

!   nfe_assert(cv%i(2) > 0 .and. cv%i(2) <= pmemd_natoms())
!   a2 = 3*cv%i(2) - 2

!   nfe_assert(a1 /= a2)

!#ifdef MPI
!   if (mytaskid.eq.0) then
!#endif /* MPI */
!      value = (x(a2)-x(a1))*cv%r(1)+(x(a2+1)-x(a1+1))*cv%r(2)+(x(a2+2)-x(a1+2))*cv%r(3)
!#ifdef MPI
!   else
!      value = ZERO
!   end if
!#endif /* MPI */

!end function v_DISTANCE_Z

!!=============================================================================

!subroutine f_DISTANCE_Z(cv, x, fcv, f)

!   use nfe_lib_mod
!   use parallel_dat_mod

!   implicit none

!   type(colvar_t), intent(in) :: cv

!   double precision, intent(in) :: x(*), fcv

!   double precision, intent(inout) :: f(*)

!   double precision :: d1(3), d2(3)

!   integer :: a1, a2
!#ifdef MPI
!   integer :: error
!#endif /* MPI */

!   nfe_assert(cv%type == COLVAR_DISTANCE_Z)

!   nfe_assert(associated(cv%i))
!   nfe_assert(size(cv%i) == 2)
!   nfe_assert(size(cv%r) == 3)

!   nfe_assert(cv%i(1) > 0 .and. cv%i(1) <= pmemd_natoms())
!   a1 = 3*cv%i(1) - 2

!   nfe_assert(cv%i(2) > 0 .and. cv%i(2) <= pmemd_natoms())
!   a2 = 3*cv%i(2) - 2

!   nfe_assert(a1 /= a2)

!   NFE_MASTER_ONLY_BEGIN
!   
!   d1(1) = -cv%r(1); d1(2) = -cv%r(2); d1(3) = -cv%r(3)
!   d1(2) =  cv%r(1); d2(2) =  cv%r(2); d2(3) =  cv%r(3)

!   f(a1:a1 + 2) = f(a1:a1 + 2) + fcv*d1
!   f(a2:a2 + 2) = f(a2:a2 + 2) + fcv*d2

!   NFE_MASTER_ONLY_END
!#ifdef MPI
!   f(1:3*pmemd_natoms()), 3*pmemd_natoms(), MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
!   nfe_assert(error.eq.0)
!#endif /* MPI */

!end subroutine f_DISTANCE_Z

!!=============================================================================

!subroutine b_DISTANCE_Z(cv, cvno)

!   use nfe_lib_mod
!   use parallel_dat_mod

!   implicit none

!   type(colvar_t), intent(inout) :: cv
!   integer,        intent(in)    :: cvno
!   integer :: error

!   nfe_assert(cv%type == COLVAR_DISTANCE_Z)
!   call check_i(cv%i, cvno, 'DISTANCE_Z', 2)

!   if ((size(cv%r).ne.3).or.(norm3(cv%(r)).ne.1)) then
!      NFE_MASTER_ONLY_BEGIN
!         write (unit = ERR_UNIT, fmt = '(a,a,'//pfmt(cvno)//',a)') &
!            NFE_ERROR, 'CV #', cvno, &
!            ' (DISTANCE_Z) : unexpected unit vector for axis'
!      NFE_MASTER_ONLY_END
!      call terminate()
!   end if

!end subroutine b_DISTANCE_Z

!!=============================================================================

!subroutine p_DISTANCE_Z(cv, lun)

!   use nfe_lib_mod

!   implicit none

!   type(colvar_t), intent(in) :: cv
!   integer, intent(in) :: lun

!   nfe_assert(is_master())
!   nfe_assert(cv%type == COLVAR_DISTANCE_Z)
!   nfe_assert(associated(cv%i))

!   write (unit = lun, fmt = '(a,a,'//pfmt(cv%r(1))//','//pfmt(cv%r(2))//','//pfmt(cv%r(3))//')') &
!      NFE_INFO, '      axis = ', cv%r(1), cv%r(2), cv%r(3)
!   call print_i(cv%i, lun)

!end subroutine p_DISTANCE_Z
!!============================================================================



!=============================================================================
! Here is the start of the code
!=============================================================================
! Pre-processing
#ifndef NFE_UTILS_H
#define NFE_UTILS_H

#ifndef NFE_DISABLE_ASSERT
#  define nfe_assert(stmt) if (.not.(stmt)) call afailed(__FILE__, __LINE__)
#  define nfe_assert_not_reached() call afailed(__FILE__, __LINE__)
#  define NFE_PURE_EXCEPT_ASSERT
#  define NFE_USE_AFAILED use nfe_lib_mod, only : afailed
#else
#  define nfe_assert(s) 
#  define nfe_assert_not_reached() 
#  define NFE_PURE_EXCEPT_ASSERT pure
#  define NFE_USE_AFAILED
#endif /* NFE_DISABLE_ASSERT */

#define NFE_OUT_OF_MEMORY call out_of_memory(__FILE__, __LINE__)

#ifdef MPI
#  define NFE_MASTER_ONLY_BEGIN if (mytaskid.eq.0) then
#  define NFE_MASTER_ONLY_END end if
#else
#  define NFE_MASTER_ONLY_BEGIN
#  define NFE_MASTER_ONLY_END
#endif /* MPI */

#define NFE_ERROR   ' ** NFE-Error ** : '
#define NFE_WARNING ' ** NFE-Warning ** : '
#define NFE_INFO    ' NFE : '

#endif /* NFE_UTILS_H */
!----------------------------------------------------------------

module nfe_colvar_mod

implicit none

private
!=================== BASIC DEFINITION ===========================
integer, public, parameter :: COLVAR_ANGLE           = 1
integer, public, parameter :: COLVAR_TORSION         = 2
integer, public, parameter :: COLVAR_DISTANCE        = 3
integer, public, parameter :: COLVAR_MULTI_RMSD      = 4
integer, public, parameter :: COLVAR_R_OF_GYRATION   = 5
integer, public, parameter :: COLVAR_HANDEDNESS      = 6
integer, public, parameter :: COLVAR_N_OF_BONDS      = 7
integer, public, parameter :: COLVAR_N_OF_STRUCTURES = 8
integer, public, parameter :: COLVAR_LCOD            = 9
integer, public, parameter :: COLVAR_COS_OF_DIHEDRAL = 10
integer, public, parameter :: COLVAR_COM_ANGLE       = 11
integer, public, parameter :: COLVAR_COM_TORSION     = 12
integer, public, parameter :: COLVAR_COM_DISTANCE    = 13
integer, public, parameter :: COLVAR_PCA             = 14
integer, public, parameter :: COLVAR_SIN_OF_DIHEDRAL = 15
integer, public, parameter :: COLVAR_PAIR_DIHEDRAL   = 16
integer, public, parameter :: COLVAR_PATTERN_DIHEDRAL= 17
integer, public, parameter :: COLVAR_DNA_RISE        = 18
integer, public, parameter :: COLVAR_DF_COM_DISTANCE = 19
integer, public, parameter :: COLVAR_ORIENTATION_ANGLE= 20
integer, public, parameter :: COLVAR_ORIENTATION_PROJ = 21
integer, public, parameter :: COLVAR_SPINANGLE        = 22
integer, public, parameter :: COLVAR_TILT             = 23
integer, public, parameter :: COLVAR_QUATERNION0      = 24
integer, public, parameter :: COLVAR_QUATERNION1      = 25
integer, public, parameter :: COLVAR_QUATERNION2      = 26
integer, public, parameter :: COLVAR_QUATERNION3      = 27

type, public :: colvar_t

   integer :: type = -1

   integer,   pointer :: i(:) => null()
   double precision, pointer :: r(:) => null()
   
   integer :: tag ! (see nfe_cv_priv.*)

   ! avgcrd : average crd of the trajectory 
   ! r      : reference crd 
   ! evec   : eigenvector from PCA 
   double precision,     pointer :: avgcrd(:) => null() 
   double precision,     pointer :: evec(:) => null() 
   double precision,     pointer :: axis(:) => null()
   integer,              pointer :: q_index => null()
     
   ! state(:) stores the sate of reference part of ref.crd
   
   integer,  pointer :: state_ref(:) => null()
   integer,  pointer :: state_pca(:) => null()
   integer,  pointer :: ipca_to_i(:) => null()
 
end type colvar_t

double precision, public                 :: cv_min,cv_max,resolution
character(len = 256), public             :: cv_type
integer,dimension(20000),public          :: cv_i 
double precision,dimension(20000),public :: cv_r
integer, public                          :: cv_ni, cv_nr, refcrd_len

double precision,dimension(20000),public :: path, harm
integer, public                          :: npath, nharm, q_index
character(len = 256), public             :: path_mode, harm_mode, refcrd_file

double precision,dimension(4), public    :: anchor_position
double precision,dimension(2), public    :: anchor_strength
double precision,dimension(3), public    :: axis

public  :: colvar
namelist / colvar /      cv_min, cv_max, resolution, cv_type, &
                         cv_i, cv_r, cv_ni, cv_nr, &
                         path, harm, npath, nharm, path_mode, harm_mode, &
                         anchor_position, anchor_strength, axis, q_index, &
                         refcrd_file 
!=========================== COLVAR UTILS ====================================

public :: check_i
public :: print_i
public :: print_pca

public :: com_check_i
public :: com_print_i

public :: com_init_weights

public :: group_com
public :: group_com_d

!=============================== COLVAR INTERFACE ============================

public :: colvar_value
public :: colvar_force

public :: colvar_difference
public :: colvar_interpolate

public :: colvar_is_periodic

public :: colvar_has_min
public :: colvar_min

public :: colvar_has_max
public :: colvar_max

public :: colvar_print
public :: colvar_cleanup

!public :: colvar_mdread  ! replaced by colvar_nlread
public :: colvar_bootstrap
public :: colvar_nlread

public :: colvar_is_quaternion
public :: colvar_has_axis
public :: colvar_has_refcrd
public :: colvar_read_refcrd

!=========================== PRELIM MATH =====================================
public :: distance
public :: distance_d

public :: angle
public :: angle_d

public :: torsion
public :: torsion_d

public :: dot3
public :: norm3
public :: cross3

public :: cosine
public :: cosine_d

private :: torsion3
private :: torsion3_d

! leaving the pucker stuff here -- may be useful in the future
#if 0
public :: pucker_z
public :: pucker_dz
#endif

private :: value4, derivative4  ! used in HANDEDNESS

private :: cos3, cos4, cos3_d, cos4_d ! used in COS_OF_DIHEDRAL,PATTERN_DIHEDRAL

private :: sin4, sin4_d ! used in SIN_OF_DIHEDRAL

private :: rise_v, rise_d ! used in DNA_RISE

#ifdef NFE_ENABLE_RMSD_CANNED
public :: rmsd_canned
#endif /* NFE_ENABLE_RMSD_CANNED */

public :: rmsd_q
public :: rmsd_q1
public :: rmsd_q2u
public :: rmsd_q3u
public :: orientation_q
!=========================== SPECIAL UTILS ===================================
! used by MULTI_RMSD
type, private :: group_t_MR  

   integer :: i0, i1, r0

   double precision, pointer :: mass(:) => null()

   double precision, pointer :: cm_crd(:) => null()
   double precision, pointer :: ref_crd(:) => null()

   double precision :: total_mass
   double precision :: ref_nrm
   double precision :: quaternion(4)
   double precision :: rmsd2

end type group_t_MR

type, private :: priv_t_MR
   double precision :: value
   double precision :: total_mass
   type(group_t_MR), pointer :: groups(:) => null()
   integer :: tag
   type(priv_t_MR), pointer :: next
end type priv_t_MR

private :: get_priv_MR
private :: new_priv_MR
private :: del_priv_MR

type(priv_t_MR), private, save, pointer :: priv_list_MR => null()
!-----------------------------------------------------------------------------
! used by R_OF_GYRATION
type, private :: priv_t_RG  
   double precision :: Rg, cm(3) ! value & center of mass
   double precision, pointer :: weights(:) ! m_i/total_mass
   integer :: tag
   type(priv_t_RG), pointer :: next
end type priv_t_RG

private :: get_priv_RG
private :: new_priv_RG
private :: del_priv_RG

type(priv_t_RG), private, save, pointer :: priv_list_RG => null()
!-----------------------------------------------------------------------------
! used by N_OF_STRUCTURES
type, private :: group_t_NS

   integer :: i0, i1, r0

   double precision, pointer :: mass(:) => null()

   double precision, pointer :: cm_crd(:) => null()
   double precision, pointer :: ref_crd(:) => null()

   double precision :: total_mass
   double precision :: ref_nrm
   double precision :: quaternion(4)
   double precision :: srmsd2
   double precision :: value

   double precision :: threshold

end type group_t_NS

type, private :: priv_t_NS
   type(group_t_NS), pointer :: groups(:) => null()
#ifdef MPI
   integer :: first_cpu
#endif /* MPI */
   integer :: tag
   type(priv_t_NS), pointer :: next
end type priv_t_NS

private :: get_priv_NS
private :: new_priv_NS
private :: del_priv_NS

type(priv_t_NS), private, save, pointer :: priv_list_NS => null()
!-----------------------------------------------------------------------------
! used by PCA
type, private :: priv_t_PCA
   
   double precision :: value, total_mass 
   double precision :: ref_cm(3), quaternion(4), U(3,3) ! COM of ref, Q, and rot matrix U  
   double precision, pointer :: mass(:) => null()
   double precision, pointer :: ref_crd(:) => null() ! translated ref  
   double precision, pointer :: cm_crd(:)  => null()  ! translated orig 
!  double precision, pointer :: fit_crd => null()  ! do we need this ?
    
   integer :: tag
   type(priv_t_PCA), pointer :: next
end type priv_t_PCA

private :: get_priv_PCA
private :: new_priv_PCA
private :: del_priv_PCA

type(priv_t_PCA), private, save, pointer :: priv_list_PCA => null()
!----------------------------------------------------------------------------
! used by QUATERNIONS
type, private :: priv_t
   double precision :: value
   double precision :: total_mass
   integer :: n_atoms
   double precision, pointer :: mass(:) => null()
   double precision, pointer :: cm_crd(:) => null()
   double precision, pointer :: ref_crd(:) => null()
   double precision :: ref_nrm
   double precision :: quaternion(4,4)
   double precision :: lambda(4)

   integer :: tag
   type(priv_t), pointer :: next
end type priv_t

private :: get_priv
private :: new_priv
private :: del_priv

type(priv_t), private, save, pointer :: priv_list => null()
!=============================================================================

contains

!=============================================================================

subroutine check_i(cvi, cvno, cvtype, expected_isize)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   integer, pointer :: cvi(:)
   integer, intent(in) :: cvno
   character(*), intent(in) :: cvtype
   integer, optional, intent(in) :: expected_isize

   integer :: a, b

   nfe_assert(cvno > 0)

   if (.not.associated(cvi)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a,a,a/)') &
            NFE_ERROR, 'CV #', cvno, ' (', cvtype, &
            ') : no integers found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cvi)

   if (present(expected_isize)) then
      if (size(cvi) /= expected_isize) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
            fmt = '(/a,a,'//pfmt(cvno)//',a,a,a,'//pfmt &
            (size(cvi))//',a,'//pfmt(expected_isize)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (', cvtype, ') : unexpected number of integers (', &
               size(cvi), ' instead of ', expected_isize, ')'
         NFE_MASTER_ONLY_END
         call terminate()
      end if ! size(cvi) /= isize
   end if ! present(expected_isize)

   do a = 1, size(cvi)
      if (cvi(a) < 1 .or. cvi(a) > pmemd_natoms()) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,a,a,'//pfmt(a)//',a,'//pfmt &
               (cvi(a))//',a,'//pfmt(pmemd_natoms())//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (', cvtype, ') : integer #', a, ' (', cvi(a), &
               ') is out of range [1, ', pmemd_natoms(), ']'
         NFE_MASTER_ONLY_END
         call terminate()
      end if
   end do

   ! check for duplicates

   do a = 1, size(cvi)
      do b = a + 1, size(cvi)
         if (cvi(a) == cvi(b)) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,a,a,'//pfmt &
                  (a)//',a,'//pfmt(b)//',a,'//pfmt(cvi(a))//',a/)') &
                  NFE_ERROR, 'CV #', cvno, ' (', cvtype, &
                  ') : integers #', a, ' and #', b, ' are equal (', cvi(a), ')'
            NFE_MASTER_ONLY_END
            call terminate()
         end if ! cvi(a) == cvi(b)
      end do
   end do

end subroutine check_i

!=============================================================================

subroutine print_i(cvi, lun)

   use nfe_lib_mod

   implicit none

   integer, pointer :: cvi(:)
   integer, intent(in) :: lun

   integer :: a
   character(4) :: aname

   nfe_assert(is_master())
   nfe_assert(associated(cvi))
   nfe_assert(size(cvi) > 0)

   write (unit = lun, fmt = '(a,a)', advance = 'NO') NFE_INFO, '  atoms = ('

   do a = 1, size(cvi)

      nfe_assert(cvi(a) > 0 .and. cvi(a) <= pmemd_natoms())
      aname = pmemd_atom_name(cvi(a))

      write (unit = lun, fmt = '('//pfmt(cvi(a))//',a,a,a)', advance = 'NO') &
         cvi(a), ' [', trim(aname), ']'

      if (a == size(cvi)) then
         write (unit = lun, fmt = '(a)') ')'
      else if (mod(a, 5) == 0) then
         write (unit = lun, fmt = '(a,/a,a)', advance = 'NO') &
            ',', NFE_INFO, '          '
      else
         write (unit = lun, fmt = '(a)', advance = 'NO') ', '
      end if

   end do

end subroutine print_i

!=============================================================================

subroutine com_check_i(cvi, cvno, cvtype, ngroups)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   integer, pointer :: cvi(:)
   integer, intent(in) :: cvno
   character(*), intent(in) :: cvtype
   integer, intent(out) :: ngroups

   integer :: a, a0

   nfe_assert(cvno > 0)

   if (.not.associated(cvi)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a,a,a/)') &
            NFE_ERROR, 'CV #', cvno, ' (', cvtype, &
            ') : no integers found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cvi)

   if (.not.size(cvi).ge.3) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a,a,a/)') &
            NFE_ERROR, 'CV #', cvno, ' (', cvtype, &
            ') : too few integers'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cvi)

   a0 = 1
   ngroups = 0

   do a = 1, size(cvi)
      if (cvi(a).eq.0) then
         ngroups = ngroups + 1
         if (a.eq.a0) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,a,a,'//pfmt(a)//',a/)') &
                  NFE_ERROR, 'CV #', cvno, &
                  ' (', cvtype, ') : unexpected zero (integer #', a, ')'
            NFE_MASTER_ONLY_END
            call terminate()
         end if ! a.eq.a0
         a0 = a + 1
      else if (cvi(a).lt.1.or.cvi(a).gt.pmemd_natoms()) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,a,a,'//pfmt(a)//',a,'//pfmt &
               (cvi(a))//',a,'//pfmt(pmemd_natoms())//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (', cvtype, ') : integer #', a, ' (', cvi(a), &
               ') is out of range [1, ', pmemd_natoms(), ']'
         NFE_MASTER_ONLY_END
         call terminate()
      end if
   end do

   if (cvi(size(cvi)).gt.0) &
      ngroups = ngroups + 1

end subroutine com_check_i

!=============================================================================

subroutine com_print_i(cvi, lun)

   use nfe_lib_mod

   implicit none

   integer, pointer :: cvi(:)
   integer, intent(in) :: lun

   integer :: a, g, c, ncvi
   character(4) :: aname

   nfe_assert(is_master())
   nfe_assert(associated(cvi))
   nfe_assert(size(cvi).gt.0)

   c = 1
   g = 1

   ncvi = size(cvi)

   if (cvi(ncvi).eq.0) &
      ncvi = ncvi - 1

   write (unit = lun, fmt = '(a,a,i1,a)', advance = 'NO') &
      NFE_INFO, '  group #', g, ' = ('

   do a = 1, ncvi
      if (cvi(a).eq.0) then
         g = g + 1
         c = 1

         write (unit = lun, fmt = '(a,a,i1,a)', advance = 'NO') &
            NFE_INFO, '  group #', g, ' = ('
      else ! cvi(a).ne.0
         nfe_assert(cvi(a).gt.0.and.cvi(a).le.pmemd_natoms())
         aname = pmemd_atom_name(cvi(a))

      write (unit = lun, fmt = '('//pfmt(cvi(a))//',a,a,a)', advance = 'NO') &
         cvi(a), ' [', trim(aname), ']'

         if (next_is_atom(a)) then
            if (mod(c, 5) == 0) then
               write (unit = lun, fmt = '(a,/a,a)', advance = 'NO') &
                  ',', NFE_INFO, '              '
            else
               write (unit = lun, fmt = '(a)', advance = 'NO') ', '
            end if
         else
            write (unit = lun, fmt = '(a)') ')'
         end if ! next_is_atom(a)
         c = c + 1
      end if ! cvi(a).eq.0
   end do

contains

logical function next_is_atom(ii)

   implicit none

   integer, intent(in) :: ii

   next_is_atom = .false.
   if (ii.lt.size(cvi)) &
      next_is_atom = cvi(ii + 1).gt.0

end function next_is_atom

end subroutine com_print_i

!=============================================================================

! cvr holds atomic masses of different groups padded by zeros
! com_init_weights() divides the individual masses by total group mass
subroutine com_init_weights(cvr)

   use nfe_lib_mod

   implicit none

   double precision, pointer :: cvr(:)

   integer :: a, a0, n
   double precision :: mass

   nfe_assert(associated(cvr))
   nfe_assert(size(cvr).gt.2)

   mass = ZERO
   a0 = 1

   do a = 1, size(cvr)
      mass = mass + cvr(a)
      if (abs(cvr(a)).lt.1.0d-10.or.a.eq.size(cvr)) then
         nfe_assert(a.ge.a0)
         nfe_assert(mass.gt.ZERO)
         do n = a0, a
            cvr(n) = cvr(n)/mass
         end do
         a0 = a + 1
         mass = ZERO
      end if
   end do

end subroutine com_init_weights

!=============================================================================

! see nfe_cv_COM_*.f
subroutine group_com(cv, x, pos, cm)

   use nfe_lib_mod
   
   implicit none

   type(colvar_t), intent(in) :: cv
   double precision, intent(in) :: x(3,pmemd_natoms())
   integer, intent(inout) :: pos
   double precision, intent(out) :: cm(3)

   integer :: a

   cm = ZERO

   do while(pos.le.size(cv%i))
      a = cv%i(pos)
      if (a.gt.0) then
         cm = cm + cv%r(pos)*x(1:3,a)
      else
         exit
      end if
      pos = pos + 1
   end do

   pos = pos + 1 ! skip zero

end subroutine group_com

!=============================================================================

! see nfe_cv_COM_*.f
subroutine group_com_d(cv, f, d, pos)

   implicit none

   type(colvar_t), intent(in) :: cv
   double precision, intent(inout) :: f(*)
   double precision, intent(in) :: d(3)
   integer, intent(inout) :: pos

   integer :: a

   do while(pos.le.size(cv%i))
      a = cv%i(pos)
      if (a.gt.0) then
         a = 3*a - 2
         f(a:a + 2) = f(a:a + 2) + cv%r(pos)*d(1:3)
      else
         exit
      end if
      pos = pos + 1
   end do

   pos = pos + 1 ! skip zero

end subroutine group_com_d

!=============================================================================

subroutine print_pca(cvi, lun)

   use nfe_lib_mod

   implicit none

   integer, pointer :: cvi(:)
   integer, intent(in) :: lun


   nfe_assert(is_master())
   nfe_assert(associated(cvi))
   nfe_assert(size(cvi) == 3)

   write (unit = lun, fmt = '(a,a,i5)') NFE_INFO, '     solute part total atom number = ', cvi(1)
   write (unit = lun, fmt = '(a,a,i5)') NFE_INFO, '        ref part total atom number = ', cvi(2)
   write (unit = lun, fmt = '(a,a,i5)') NFE_INFO, '        pca part total atom number = ', cvi(3)   

end subroutine print_pca

!=============================================================================
!=============================================================================

! the value is needed only on master
function colvar_value(cv, x) result(value)

   use nfe_lib_mod
   
   implicit none

   double precision :: value

   type(colvar_t) :: cv ! mutable
   double precision, intent(in) :: x(3,pmemd_natoms())

   select case(cv%type)
      case(COLVAR_ANGLE)
         value = v_ANGLE(cv, x)
      case(COLVAR_TORSION)
         value = v_TORSION(cv, x)
      case(COLVAR_DISTANCE)
         value = v_DISTANCE(cv, x)
      case(COLVAR_MULTI_RMSD)
         value = v_MULTI_RMSD(cv, x)
      case(COLVAR_R_OF_GYRATION)
         value = v_R_OF_GYRATION(cv, x)
      case(COLVAR_HANDEDNESS)
         value = v_HANDEDNESS(cv, x)
      case(COLVAR_N_OF_BONDS)
         value = v_N_OF_BONDS(cv, x)
      case(COLVAR_N_OF_STRUCTURES)
         value = v_N_OF_STRUCTURES(cv, x)
      case(COLVAR_LCOD)
         value = v_LCOD(cv, x)
      case(COLVAR_COS_OF_DIHEDRAL)
         value = v_COS_OF_DIHEDRAL(cv, x)
      case(COLVAR_COM_ANGLE)
         value = v_COM_ANGLE(cv, x)
      case(COLVAR_COM_TORSION)
         value = v_COM_TORSION(cv, x)
      case(COLVAR_COM_DISTANCE)
         value = v_COM_DISTANCE(cv, x)   
      case(COLVAR_PCA)
         value = v_PCA(cv,x)
      case(COLVAR_SIN_OF_DIHEDRAL)
         value = v_SIN_OF_DIHEDRAL(cv,x)
      case(COLVAR_PAIR_DIHEDRAL)
         value = v_PAIR_DIHEDRAL(cv,x)
      case(COLVAR_PATTERN_DIHEDRAL)
         value = v_PATTERN_DIHEDRAL(cv,x) 
      case(COLVAR_DNA_RISE)
         value = v_DNA_RISE(cv,x)
      case(COLVAR_DF_COM_DISTANCE)
         value = v_DF_COM_DISTANCE(cv, x)
      case(COLVAR_ORIENTATION_ANGLE)
         value = v_ORIENTATION_ANGLE(cv,x)
      case(COLVAR_ORIENTATION_PROJ)
         value = v_ORIENTATION_PROJ(cv,x)
      case(COLVAR_SPINANGLE)
         value = v_SPINANGLE(cv,x)
      case(COLVAR_TILT)
         value = v_TILT(cv,x)
      case(COLVAR_QUATERNION0)
         value = v_QUATERNION0(cv,x)
      case(COLVAR_QUATERNION1)
         value = v_QUATERNION1(cv,x)
      case(COLVAR_QUATERNION2)
         value = v_QUATERNION2(cv,x)
      case(COLVAR_QUATERNION3)
         value = v_QUATERNION3(cv,x)
      case default
         nfe_assert_not_reached()
         value = dble(0)
   end select

end function colvar_value

!=============================================================================

subroutine colvar_force(cv, x, fcv, f)

   use nfe_lib_mod
   
   implicit none

   type(colvar_t) :: cv ! mutable

   double precision, intent(in) :: x(3,pmemd_natoms()), fcv
   double precision, intent(inout) :: f(3,pmemd_natoms())

   select case(cv%type)
      case(COLVAR_ANGLE)
         call f_ANGLE(cv, x, fcv, f)
      case(COLVAR_TORSION)
         call f_TORSION(cv, x, fcv, f)
      case(COLVAR_DISTANCE)
         call f_DISTANCE(cv, x, fcv, f)
      case(COLVAR_MULTI_RMSD)
         call f_MULTI_RMSD(cv, x, fcv, f)
      case(COLVAR_R_OF_GYRATION)
         call f_R_OF_GYRATION(cv, x, fcv, f)
      case(COLVAR_HANDEDNESS)
         call f_HANDEDNESS(cv, x, fcv, f)
      case(COLVAR_N_OF_BONDS)
         call f_N_OF_BONDS(cv, x, fcv, f)
      case(COLVAR_N_OF_STRUCTURES)
         call f_N_OF_STRUCTURES(cv, x, fcv, f)
      case(COLVAR_LCOD)
         call f_LCOD(cv, x, fcv, f)
      case(COLVAR_COS_OF_DIHEDRAL)
         call f_COS_OF_DIHEDRAL(cv, x, fcv, f)
      case(COLVAR_COM_ANGLE)
         call f_COM_ANGLE(cv, x, fcv, f)
      case(COLVAR_COM_TORSION)
         call f_COM_TORSION(cv, x, fcv, f)
      case(COLVAR_COM_DISTANCE)
         call f_COM_DISTANCE(cv, x, fcv, f)     
      case(COLVAR_PCA)
         call f_PCA(cv, x, fcv, f)
      case(COLVAR_SIN_OF_DIHEDRAL)
         call f_SIN_OF_DIHEDRAL(cv, x, fcv, f)
      case(COLVAR_PAIR_DIHEDRAL)
         call f_PAIR_DIHEDRAL(cv, x, fcv, f)
      case(COLVAR_PATTERN_DIHEDRAL)
         call f_PATTERN_DIHEDRAL(cv, x, fcv, f)
      case(COLVAR_DNA_RISE)
         call f_DNA_RISE(cv,x,fcv,f)
      case(COLVAR_DF_COM_DISTANCE)
         call f_DF_COM_DISTANCE(cv, x, fcv, f)
      case(COLVAR_ORIENTATION_ANGLE)
         call f_ORIENTATION_ANGLE(cv, fcv, f)
      case(COLVAR_ORIENTATION_PROJ)
         call f_ORIENTATION_PROJ(cv, fcv, f)
      case(COLVAR_SPINANGLE)
         call f_SPINANGLE(cv, fcv, f)
      case(COLVAR_TILT)
         call f_TILT(cv, fcv, f)
      case(COLVAR_QUATERNION0)
         call f_QUATERNION0(cv, fcv, f)
      case(COLVAR_QUATERNION1)
         call f_QUATERNION1(cv, fcv, f)
      case(COLVAR_QUATERNION2)
         call f_QUATERNION2(cv, fcv, f)
      case(COLVAR_QUATERNION3)
         call f_QUATERNION3(cv, fcv, f)
      case default
         nfe_assert_not_reached()
   end select

end subroutine colvar_force

!=============================================================================

function colvar_difference(cv, v1, v2) result(diff)

   use gbl_constants_mod, only : PI
   use nfe_lib_mod

   implicit none

   double precision :: diff
   type(colvar_t), intent(in) :: cv
   double precision, intent(in) :: v1, v2

   double precision :: t1, t2

   t1 = fix_value(v1)
   t2 = fix_value(v2)

   diff = t1 - t2
  
   if (cv%type.eq.COLVAR_TORSION.or.cv%type.eq.COLVAR_COM_TORSION) then
      nfe_assert(- PI.le.t1.and.t1.le.PI)
      nfe_assert(- PI.le.t2.and.t2.le.PI)
      if (diff.gt.PI) then
         diff = diff - PI - PI
      else if (diff.lt.-PI) then
         diff = diff + PI + PI
      end if
   end if

contains

function fix_value(v) result(t)

   implicit none

   double precision :: t
   double precision, intent(in) :: v

   if (cv%type.eq.COLVAR_ANGLE.or.cv%type.eq.COLVAR_COM_ANGLE) then
      if (ZERO.le.v.and.v.le.PI) then
         t = v
      else
         t = acos(cos(v))
      end if
   else if (cv%type.eq.COLVAR_TORSION.or.cv%type.eq.COLVAR_COM_TORSION) then
      if (- PI.le.v.and.v.le.PI) then
        t = v
      else
         t = atan2(sin(v), cos(v))
      endif
   else
      t = v
   end if

end function fix_value

end function colvar_difference

!=============================================================================

function colvar_interpolate(cv, a1, v1, a2, v2) result(interp)

   implicit none

   double precision :: interp
   type(colvar_t), intent(in) :: cv
   double precision, intent(in) :: a1, v1, a2, v2

   double precision :: ts, tc

   if (cv%type.eq.COLVAR_TORSION.or.cv%type.eq.COLVAR_COM_TORSION) then
      ts = a1*sin(v1)+a2*sin(v2)
      tc = a1*cos(v1)+a2*cos(v2)
      interp = atan2(ts,tc)
   else
      interp = a1*v1+a2*v2
   end if

end function colvar_interpolate

!=============================================================================

logical function colvar_is_periodic(cv)

   implicit none

   type(colvar_t), intent(in) :: cv

   if (cv%type.eq.COLVAR_TORSION.or.cv%type.eq.COLVAR_COM_TORSION) then
      colvar_is_periodic = .true.
   else
      colvar_is_periodic = .false.
   end if
 
end function colvar_is_periodic

!=============================================================================

logical function colvar_has_axis(cv)

   implicit none

   type(colvar_t), intent(in) :: cv

   if (cv%type.eq.COLVAR_TILT.or.cv%type.eq.COLVAR_SPINANGLE) then
      colvar_has_axis = .true.
   else
      colvar_has_axis = .false.
   end if

end function colvar_has_axis

!=============================================================================

logical function colvar_is_quaternion(cv)

   implicit none

   type(colvar_t), intent(in) :: cv

   if (cv%type.eq.COLVAR_QUATERNION0.or.cv%type.eq.COLVAR_QUATERNION1 &
      .or.cv%type.eq.COLVAR_QUATERNION2.or.cv%type.eq.COLVAR_QUATERNION3) then
       colvar_is_quaternion = .true.
   else
      colvar_is_quaternion = .false.
   end if

end function colvar_is_quaternion

!=============================================================================

logical function colvar_has_refcrd(cv)

   implicit none

   type(colvar_t), intent(in) :: cv

   if (cv%type.eq.COLVAR_ORIENTATION_ANGLE.or.cv%type.eq.COLVAR_ORIENTATION_PROJ &
       .or.cv%type.eq.COLVAR_TILT.or.cv%type.eq.COLVAR_SPINANGLE &
       .or.cv%type.eq.COLVAR_QUATERNION0.or.cv%type.eq.COLVAR_QUATERNION1 &
       .or.cv%type.eq.COLVAR_QUATERNION2.or.cv%type.eq.COLVAR_QUATERNION3) then
      colvar_has_refcrd = .true.
   else
      colvar_has_refcrd = .false.
   end if

end function colvar_has_refcrd

!=============================================================================

logical function colvar_has_min(cv)

   implicit none

   type(colvar_t), intent(in) :: cv

   select case(cv%type)
      case(COLVAR_ANGLE:COLVAR_R_OF_GYRATION)
         colvar_has_min = .true.
      case(COLVAR_N_OF_BONDS:COLVAR_N_OF_STRUCTURES)
         colvar_has_min = .true.
      case(COLVAR_COM_ANGLE:COLVAR_COM_DISTANCE)
         colvar_has_min = .true.
      case(COLVAR_ORIENTATION_ANGLE)
         colvar_has_min = .true.
      case(COLVAR_ORIENTATION_PROJ)
         colvar_has_min = .true.
      case(COLVAR_SPINANGLE)
         colvar_has_min = .true.
      case(COLVAR_TILT)
         colvar_has_min = .true.
      case default
         colvar_has_min = .false.
   end select

end function colvar_has_min

!=============================================================================

double precision function colvar_min(cv)

   use gbl_constants_mod, only : PI
   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv

   nfe_assert(colvar_has_min(cv))

   select case(cv%type)
      case(COLVAR_ANGLE)
         colvar_min = ZERO
      case(COLVAR_TORSION)
         colvar_min = -PI
      case(COLVAR_DISTANCE:COLVAR_R_OF_GYRATION)
         colvar_min = ZERO
      case(COLVAR_N_OF_BONDS:COLVAR_N_OF_STRUCTURES)
         colvar_min = ZERO
      case(COLVAR_COM_ANGLE)
         colvar_min = ZERO
      case(COLVAR_COM_TORSION)
         colvar_min = -PI
      case(COLVAR_COM_DISTANCE)
         colvar_min = ZERO
      case(COLVAR_ORIENTATION_ANGLE)
         colvar_min = ZERO
      case(COLVAR_ORIENTATION_PROJ)
         colvar_min = -ONE
      case(COLVAR_SPINANGLE)
         colvar_min = -PI
      case(COLVAR_TILT)
         colvar_min = -ONE
      case default
         nfe_assert_not_reached()
         colvar_min = ZERO
   end select

end function colvar_min

!=============================================================================

logical function colvar_has_max(cv)

   implicit none

   type(colvar_t), intent(in) :: cv

   select case(cv%type)
      case(COLVAR_ANGLE:COLVAR_TORSION)
         colvar_has_max = .true.
      case(COLVAR_COM_ANGLE:COLVAR_COM_TORSION)
         colvar_has_max = .true.
      case(COLVAR_ORIENTATION_ANGLE)
         colvar_has_max = .true.
      case(COLVAR_ORIENTATION_PROJ)
         colvar_has_max = .true.
      case(COLVAR_SPINANGLE)
         colvar_has_max = .true.
      case(COLVAR_TILT)
         colvar_has_max = .true.
      case default
         colvar_has_max = .false.
   end select

end function colvar_has_max

!=============================================================================

double precision function colvar_max(cv)

   use gbl_constants_mod, only : PI
   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv

   nfe_assert(colvar_has_max(cv))

   select case(cv%type)
      case(COLVAR_ANGLE:COLVAR_TORSION)
         colvar_max = PI
      case(COLVAR_COM_ANGLE:COLVAR_COM_TORSION)
         colvar_max = PI
      case(COLVAR_ORIENTATION_ANGLE)
         colvar_max = PI
      case(COLVAR_ORIENTATION_PROJ)
         colvar_max = ONE
      case(COLVAR_SPINANGLE)
         colvar_max = PI
      case(COLVAR_TILT)
         colvar_max = ONE
      case default
         nfe_assert_not_reached()
         colvar_max = ZERO
   end select

end function colvar_max

!=============================================================================

subroutine colvar_print(cv, lun)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   integer,        intent(in) :: lun

   write (unit = lun, fmt = '(a,a)', advance = 'NO') NFE_INFO, '  type = '''

   select case(cv%type)
      case(COLVAR_ANGLE)
         write (unit = lun, fmt = '(a)') 'ANGLE'''
         call p_ANGLE(cv, lun)
      case(COLVAR_TORSION)
         write (unit = lun, fmt = '(a)') 'TORSION'''
         call p_TORSION(cv, lun)
      case(COLVAR_DISTANCE)
         write (unit = lun, fmt = '(a)') 'DISTANCE'''
         call p_DISTANCE(cv, lun)
      case(COLVAR_MULTI_RMSD)
         write (unit = lun, fmt = '(a)') 'MULTI_RMSD'''
         call p_MULTI_RMSD(cv, lun)
      case(COLVAR_R_OF_GYRATION)
         write (unit = lun, fmt = '(a)') 'R_OF_GYRATION'''
         call p_R_OF_GYRATION(cv, lun)
      case(COLVAR_HANDEDNESS)
         write (unit = lun, fmt = '(a)') 'HANDEDNESS'''
         call p_HANDEDNESS(cv, lun)
      case(COLVAR_N_OF_BONDS)
         write (unit = lun, fmt = '(a)') 'N_OF_BONDS'''
         call p_N_OF_BONDS(cv, lun)
      case(COLVAR_N_OF_STRUCTURES)
         write (unit = lun, fmt = '(a)') 'N_OF_STRUCTURES'''
         call p_N_OF_STRUCTURES(cv, lun)
      case(COLVAR_LCOD)
         write (unit = lun, fmt = '(a)') &
            'LCOD'' (Linear Combination Of Distances)'
         call p_LCOD(cv, lun)
      case(COLVAR_COS_OF_DIHEDRAL)
         write (unit = lun, fmt = '(a)') 'COS_OF_DIHEDRAL'''
         call p_COS_OF_DIHEDRAL(cv, lun)
      case(COLVAR_COM_ANGLE)
         write (unit = lun, fmt = '(a)') 'COM_ANGLE'''
         call p_COM_ANGLE(cv, lun)
      case(COLVAR_COM_TORSION)
         write (unit = lun, fmt = '(a)') 'COM_TORSION'''
         call p_COM_TORSION(cv, lun)
      case(COLVAR_COM_DISTANCE)
         write (unit = lun, fmt = '(a)') 'COM_DISTANCE'''
         call p_COM_DISTANCE(cv, lun)
      case(COLVAR_PCA) 
         write (unit = lun, fmt = '(a)') 'PRINCIPAL COMPONENT'''
         call p_PCA(cv, lun)
      case(COLVAR_SIN_OF_DIHEDRAL)
         write (unit = lun, fmt = '(a)') 'SIN_OF_DIHEDRAL'''
         call p_SIN_OF_DIHEDRAL(cv, lun)
      case(COLVAR_PAIR_DIHEDRAL)
         write (unit = lun, fmt = '(a)') 'PAIR_DIHEDRAL'''
         call p_PAIR_DIHEDRAL(cv, lun)
      case(COLVAR_PATTERN_DIHEDRAL)
         write (unit = lun, fmt = '(a)') 'PATTERN_DIHEDRAL'''
         call p_PATTERN_DIHEDRAL(cv, lun)
      case(COLVAR_DNA_RISE)
         write (unit = lun, fmt = '(a)') 'DNA_RISE'''
         call p_DNA_RISE(cv,lun)
      case(COLVAR_DF_COM_DISTANCE)
         write (unit = lun, fmt = '(a)') 'DF_COM_DISTANCE'''
         call p_DF_COM_DISTANCE(cv, lun)
      case(COLVAR_ORIENTATION_ANGLE)
         write (unit = lun, fmt = '(a)') 'ORIENTATION_ANGLE'''
         call p_ORIENTATION_ANGLE(cv, lun)
      case(COLVAR_ORIENTATION_PROJ)
         write (unit = lun, fmt = '(a)') 'ORIENTATION_PROJ'''
         call p_ORIENTATION_PROJ(cv, lun)
      case(COLVAR_SPINANGLE)
         write (unit = lun, fmt = '(a)') 'SPINANGLE'''
         call p_SPINANGLE(cv, lun)
      case(COLVAR_TILT)
         write (unit = lun, fmt = '(a)') 'TILT'''
         call p_TILT(cv, lun)
      case(COLVAR_QUATERNION0)
         write (unit = lun, fmt = '(a)') 'QUATERNION0'''
         call p_QUATERNION0(cv, lun)
      case(COLVAR_QUATERNION1)
         write (unit = lun, fmt = '(a)') 'QUATERNION1'''
         call p_QUATERNION1(cv, lun)
      case(COLVAR_QUATERNION2)
         write (unit = lun, fmt = '(a)') 'QUATERNION2'''
         call p_QUATERNION2(cv, lun)
      case(COLVAR_QUATERNION3)
         write (unit = lun, fmt = '(a)') 'QUATERNION3'''
         call p_QUATERNION3(cv, lun)
      case default
         nfe_assert_not_reached()
         continue
   end select

end subroutine colvar_print

!=============================================================================

subroutine colvar_cleanup(cv) 

   implicit none

   type(colvar_t), intent(inout) :: cv

   select case(cv%type)
      case(COLVAR_MULTI_RMSD)
         call c_MULTI_RMSD(cv)
      case(COLVAR_R_OF_GYRATION)
         call c_R_OF_GYRATION(cv)
      case(COLVAR_N_OF_STRUCTURES)
         call c_N_OF_STRUCTURES(cv)
      case (COLVAR_PCA) 
         call c_PCA(cv)
      case (COLVAR_ORIENTATION_ANGLE)
         call c_ORIENTATION_ANGLE(cv)
      case (COLVAR_ORIENTATION_PROJ)
         call c_ORIENTATION_PROJ(cv)
      case (COLVAR_SPINANGLE)
         call c_SPINANGLE(cv)
      case (COLVAR_TILT)
         call c_TILT(cv)
      case (COLVAR_QUATERNION0)
         call c_QUATERNION0(cv)
      case (COLVAR_QUATERNION1)
         call c_QUATERNION1(cv)
      case (COLVAR_QUATERNION2)
         call c_QUATERNION2(cv)
      case (COLVAR_QUATERNION3)
         call c_QUATERNION3(cv)
      case default
         continue
   end select

   if (associated(cv%i)) &
      deallocate(cv%i)

   if (associated(cv%r)) &
      deallocate(cv%r)

   if (associated(cv%avgcrd)) &
      deallocate(cv%avgcrd)

   if (associated(cv%evec)) &
      deallocate(cv%evec)
 
   if (associated(cv%state_ref)) &
      deallocate(cv%state_ref)
      
   if (associated(cv%state_pca)) &
      deallocate(cv%state_pca) 
   
   if (associated(cv%ipca_to_i)) &
      deallocate(cv%ipca_to_i) 

   if (associated(cv%q_index)) &
      deallocate(cv%q_index)

   if (associated(cv%axis)) &
      deallocate(cv%axis)
   
   cv%type = -1

end subroutine colvar_cleanup

!=============================================================================

subroutine colvar_read_refcrd(coor, refcrd)

  use nfe_lib_mod

  implicit none

  double precision, pointer :: coor(:) !=> null()

  character(len = *), intent(in)  :: refcrd

  integer ::i, j, error
  character :: dummy

  open(REF_UNIT1, FILE = refcrd, iostat = error, status = 'old')
     if (error.ne.0) then
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a,a/)') &
         NFE_ERROR, 'Failed to open reference coordinates file :''', trim(refcrd), ''''
      call terminate()
     end if
  ! skip the first line of the CRD files 
  read(REF_UNIT1, '(A1)') dummy
  read(REF_UNIT1, *) j                             ! Read total number of atoms
  read(REF_UNIT1, '(6F12.7)') (coor(i), i=1, 3*j)  ! Read coordinates 

  close(REF_UNIT1)

end subroutine colvar_read_refcrd

!=============================================================================

!subroutine colvar_mdread(cv, node, cvno)

!   use nfe_utils_mod
!   use nfe_value_mod
!   use nfe_cftree_mod
!   use nfe_constants_mod
!   use nfe_colvar_type_mod
!   use nfe_pmemd_proxy_mod
!   use nfe_read_pca_mod

!   implicit none

!   type(colvar_t), intent(inout) :: cv
!   type(node_t),   intent(in)    :: node
!   integer,        intent(in)    :: cvno

!   integer :: n, error
!   integer :: found, crdsize 
!   logical :: found2
!   integer :: nsolut
!   
!   type(value_node_t), pointer :: alist, aiter

!   character(len = STRING_LENGTH) :: type

!   ! declare three strings: ref_file avg_file evec_file 
!   character(len = STRING_LENGTH) :: ref_file
!   character(len = STRING_LENGTH) :: avg_file
!   character(len = STRING_LENGTH) :: evec_file
!   character(len = STRING_LENGTH) :: index_file 
!   integer :: first, last 

!   nfe_assert(is_master())

!   nfe_assert(.not. associated(cv%i))
!   nfe_assert(.not. associated(cv%r))
!   nfe_assert(.not. associated(cv%avgcrd)) 
!   nfe_assert(.not. associated(cv%evec))
!   nfe_assert(.not. associated(cv%state_ref))
!   nfe_assert(.not. associated(cv%state_pca)) 

!   nfe_assert(node_title(node) == 'variable')

!   !
!   ! type
!   !

!   if (.not.node_lookup_string(node, 'type', type)) then
!      write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//'/)') &
!            NFE_ERROR, 'type is not specified for CV #', cvno
!      call terminate()
!   end if

!   if (type == 'ANGLE') then
!      cv%type = COLVAR_ANGLE
!   else if (type == 'TORSION') then
!      cv%type = COLVAR_TORSION
!   else if (type == 'DISTANCE') then
!      cv%type = COLVAR_DISTANCE
!   else if (type == 'MULTI_RMSD') then
!      cv%type = COLVAR_MULTI_RMSD
!   else if (type == 'R_OF_GYRATION') then
!      cv%type = COLVAR_R_OF_GYRATION
!   else if (type == 'HANDEDNESS') then
!      cv%type = COLVAR_HANDEDNESS
!   else if (type == 'N_OF_BONDS') then
!      cv%type = COLVAR_N_OF_BONDS
!   else if (type == 'N_OF_STRUCTURES') then
!      cv%type = COLVAR_N_OF_STRUCTURES
!   else if (type == 'LCOD') then
!      cv%type = COLVAR_LCOD
!   else if (type == 'COS_OF_DIHEDRAL') then
!      cv%type = COLVAR_COS_OF_DIHEDRAL
!   else if (type == 'COM_ANGLE') then
!      cv%type = COLVAR_COM_ANGLE
!   else if (type == 'COM_TORSION') then
!      cv%type = COLVAR_COM_TORSION
!   else if (type == 'COM_DISTANCE') then
!      cv%type = COLVAR_COM_DISTANCE
!   else if (type == 'PCA') then
!      cv%type = COLVAR_PCA
!      nsolut = pmemd_nsolut()
!   else if (type == 'SIN_OF_DIHEDRAL') then
!      cv%type = COLVAR_SIN_OF_DIHEDRAL
!   else if (type == 'PAIR_DIHEDRAL') then
!      cv%type = COLVAR_PAIR_DIHEDRAL
!   else if (type == 'PATTERN_DIHEDRAL') then
!      cv%type = COLVAR_PATTERN_DIHEDRAL
!   else if (type == 'DNA_RISE') then
!      cv%type = COLVAR_DNA_RISE
!   else if (type == 'DF_COM_DISTANCE') then
!      cv%type = COLVAR_DF_COM_DISTANCE      
!   else
!      write (unit = ERR_UNIT, fmt = '(/a,a,a,a,'//pfmt(cvno)//',a/)') &
!            NFE_ERROR, 'CV type ''', trim(type), &
!            ''' is not supported so far (CV #', cvno, ')'
!      call terminate()
!   end if

!   !
!   ! cv%i
!   !

!   if (node_lookup_list(node, 'i', alist)) then

!      n = 0
!      aiter => alist

!      do while (associated(aiter))
!         n = n + 1
!         if (.not. value_is_integer(aiter%value)) then
!            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
!               NFE_ERROR, 'CV #', cvno, ' : unexpected &
!               &(not an integer) element of ''i'' list'
!            call terminate()
!         end if
!         aiter => aiter%next
!      end do

!      if (n > 0) then
!         allocate(cv%i(n), stat = error)
!         if (error /= 0) &
!            NFE_OUT_OF_MEMORY

!         n = 0
!         aiter => alist
!         do while (associated(aiter))
!            n = n + 1
!            cv%i(n) = value_get_integer(aiter%value)
!            aiter => aiter%next
!         end do
!      end if

!   end if ! node_lookup_list(vnode, 'i', alist))

!   !
!   ! cv%r
!   !

!  if (type /= 'PCA') then 
!   if (node_lookup_list(node, 'r', alist)) then

!      n = 0
!      aiter => alist

!      do while (associated(aiter))
!         n = n + 1
!         if (.not. value_is_real(aiter%value)) then
!            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
!               NFE_ERROR, 'CV #', cvno, ' : unexpected &
!               &(not a real number) element of ''r'' list'
!            call terminate()
!         end if
!         aiter => aiter%next
!      end do

!      if (n > 0) then
!         allocate(cv%r(n), stat = error)
!         if (error /= 0) &
!            NFE_OUT_OF_MEMORY

!         n = 0
!         aiter => alist
!         do while (associated(aiter))
!            n = n + 1
!            cv%r(n) = value_get_real(aiter%value)
!            aiter => aiter%next
!         end do
!      end if

!   end if ! node_lookup_list(vnode, 'r', alist))
!  endif  
!  
!   ! read refcrd from "refcrd" 
!   ! store in cv% r 
!   if (type == 'PCA') then
!!       first = 1
!!       last  = pmemd_nsolut()
!!       crdsize = (last-first+1)*3
!       
!!       write(*,*) "pmemd_nsolut=", pmemd_nsolut();
!       
!!       nfe_assert(cv%i(1)*3 == crdsize)
!         
!       found2 = node_lookup_string(node, 'refcrd', ref_file)
!       if (.not. found2) then
!            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(n)//'/)') &
!               NFE_ERROR, 'could not find ''reference coordinates'' for CV #', n
!            call terminate()
!       else
!       ! write (unit = OUT_UNIT, fmt = *) 'crdsize is', crdsize
!       ! allocate and read reference crd file into CV      
!       	   if(cv%i(1) > 0) then 
!              allocate(cv%r(cv%i(1)*3), stat = error)
!              if (error /= 0) &
!                 NFE_OUT_OF_MEMORY
!               call read_refcrd(cv, ref_file)
!               write (unit = OUT_UNIT, fmt = '(a,a,a,a)', advance = 'NO') NFE_INFO, &
!               'ref_file = ', trim(ref_file), ' ('
!               write (unit = OUT_UNIT, fmt = '(a)') 'loaded)'
!              end if
!       endif

!      ! cv % evec 
!      found2 = node_lookup_string(node, 'evec', evec_file)
!      if (.not. found2 ) then
!           write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(n)//'/)') &
!              NFE_ERROR, 'could not find ''reference coordinates'' for CV #', n
!            call terminate()
!      else
!      ! allocate and read evec crd into CV 
!!            if(crdsize > 0) then
!!              allocate(cv%evec(crdsize), stat = error)
!            if(cv%i(3)>0) then
!             	 allocate( cv%evec(cv%i(3)*3), stat = error)
!              ! write (unit = OUT_UNIT, fmt = *) 'crdsize is', crdsize
!              if (error /= 0) &
!                 NFE_OUT_OF_MEMORY
!!              call read_evec(cv, evec_file, first, last)
!              call read_evec(cv, evec_file)
!              write (unit = OUT_UNIT, fmt = '(a,a,a,a)', advance = 'NO') NFE_INFO, &
!              'evec_file = ', trim(evec_file), ' ('
!              write (unit = OUT_UNIT, fmt = '(a)') 'loaded)'
!            endif
!      endif

!     
!      ! cv % avgcrd
!      found2 = node_lookup_string(node, 'avgcrd', avg_file)
!      if (.not. found2) then
!            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(n)//'/)') &
!               NFE_ERROR, 'could not find ''average coordinates'' for CV #', n
!            call terminate()
!      else
!      	
!      	
!      ! allocate and read average crd into CV 
!!          if(crdsize > 0) then
!!              allocate(cv%avgcrd(crdsize), stat = error)
!           if(cv%i(3)>0) then
!              allocate( cv%avgcrd(cv%i(3)*3), stat = error)
!              if (error /= 0) &
!                 NFE_OUT_OF_MEMORY
!!             call read_avgcrd(cv, avg_file, first, last)
!              
!              call read_avgcrd(cv, avg_file)
!              write (unit = OUT_UNIT, fmt = '(a,a,a,a)', advance = 'NO') NFE_INFO, &
!              'avg_file = ', trim(avg_file), ' ('
!              write (unit = OUT_UNIT, fmt = '(a)') 'loaded)'
!           endif
!       endif
!   
!   
!      ! cv % index
!      found2 = node_lookup_string(node, 'index', index_file)
!      if (.not. found2) then
!            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(n)//'/)') &
!               NFE_ERROR, 'could not find ''index for ref and pca part'' for CV #', n
!            call terminate()
!      else
!      ! allocate and read average crd into CV 
!          if(cv%i(1) > 0 .and. cv%i(3) > 0) then
!              allocate(cv%state_ref(cv%i(1)), cv%state_pca(cv%i(1)), cv%ipca_to_i(cv%i(3)), stat = error)
!              if (error /= 0) &
!                 NFE_OUT_OF_MEMORY
!              call read_index(cv, index_file)
!              write (unit = OUT_UNIT, fmt = '(a,a,a,a)', advance = 'NO') NFE_INFO, &
!              'index_file = ', trim(index_file), ' ('
!              write (unit = OUT_UNIT, fmt = '(a)') 'loaded)'
!           endif
!       endif

!   
!  endif 
!     
! 
!end subroutine colvar_mdread

!=============================================================================

subroutine colvar_bootstrap(cv, cvno, amass)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)

#ifdef MPI  
   integer :: bcastdata(10), ierr

   !
   ! bcast type/i/r first
   !

#ifndef NFE_DISABLE_ASSERT
   if (mytaskid == 0) then
      nfe_assert(cv%type > 0)
   else
      nfe_assert(.not. associated(cv%i))
      nfe_assert(.not. associated(cv%r))
      if (cv%type.eq.COLVAR_QUATERNION0.or.cv%type.eq.COLVAR_QUATERNION1 &
          .or.cv%type.eq.COLVAR_QUATERNION2.or.cv%type.eq.COLVAR_QUATERNION3) then
          nfe_assert(.not. associated(cv%q_index))
      endif
      if (cv%type.eq.COLVAR_TILT.or.cv%type.eq.COLVAR_SPINANGLE) then
         nfe_assert(.not. associated(cv%axis))
      endif
      if (cv%type == COLVAR_PCA) then 
        nfe_assert(.not. associated(cv%avgcrd))
        nfe_assert(.not. associated(cv%evec))
        nfe_assert(.not. associated(cv%state_ref))
        nfe_assert(.not. associated(cv%state_pca))
        nfe_assert(.not. associated(cv%ipca_to_i)) 
      endif
   end if
#endif /* NFE_DISABLE_ASSERT */

   if (mytaskid == 0) then
      bcastdata(1) = cv%type

      bcastdata(2) = 0
      if (associated(cv%i)) &
         bcastdata(2) = size(cv%i)

      bcastdata(3) = 0
      if (associated(cv%r)) &
         bcastdata(3) = size(cv%r)
     
      bcastdata(4) = 0 
      bcastdata(5) = 0
      bcastdata(6) = 0
      bcastdata(7) = 0
      bcastdata(8) = 0
      bcastdata(9) = 0             !q_index
      bcastdata(10) = 0            !axis
     
     if (cv%type == COLVAR_PCA) then 
       
        if (associated(cv%avgcrd)) &
          bcastdata(4) = size(cv%avgcrd) 
        
        if (associated(cv%evec)) &
          bcastdata(5) = size(cv%evec)  
          
        if (associated(cv%state_ref)) &
          bcastdata(6) = size(cv%state_ref)
          
        if (associated(cv%state_pca)) &
          bcastdata(7) = size(cv%state_pca)
        
        if (associated(cv%ipca_to_i)) &
          bcastdata(8) = size(cv%ipca_to_i)  
          
     endif
     if (cv%type.eq.COLVAR_QUATERNION0.or.cv%type.eq.COLVAR_QUATERNION1 &
          .or.cv%type.eq.COLVAR_QUATERNION2.or.cv%type.eq.COLVAR_QUATERNION3) then
        if (associated(cv%q_index)) &
           bcastdata(9) = cv%q_index
     endif
     if (cv%type.eq.COLVAR_TILT.or.cv%type.eq.COLVAR_SPINANGLE) then
       if (associated(cv%axis)) &
         bcastdata(10) = size(cv%axis)
     endif
     
   end if ! mytaskid == 0

   call mpi_bcast(bcastdata, size(bcastdata), MPI_INTEGER, 0, pmemd_comm, ierr)
   nfe_assert(ierr == 0)

   if (mytaskid /= 0) &
      cv%type = bcastdata(1)

   !
   ! cv%i
   !

   if (bcastdata(2) > 0) then
      if (.not. associated(cv%i)) then
         nfe_assert(mytaskid > 0)
         allocate(cv%i(bcastdata(2)), stat = ierr)
         if (ierr /= 0) &
            NFE_OUT_OF_MEMORY
      end if ! .not. associated(cv%i)

      call mpi_bcast(cv%i, bcastdata(2), MPI_INTEGER, 0, pmemd_comm, ierr)
      nfe_assert(ierr == 0)
   else
      nullify(cv%i)
   end if ! bcastdata(2) > 0

   !
   ! cv%r
   !

   if (bcastdata(3) > 0) then
      if (.not. associated(cv%r)) then
         nfe_assert(mytaskid > 0)
         allocate(cv%r(bcastdata(3)), stat = ierr)
         if (ierr /= 0) &
            NFE_OUT_OF_MEMORY
      end if ! .not. associated(cv%r)

      call mpi_bcast(cv%r, bcastdata(3), MPI_DOUBLE_PRECISION, &
                     0, pmemd_comm, ierr)
      nfe_assert(ierr == 0)
   else
      nullify(cv%r)
   end if ! bcastdata(3) > 0

   ! 
   ! cv % avgcrd    
   ! 
   if (bcastdata(4) > 0) then
      if (.not. associated(cv%avgcrd)) then
         nfe_assert(mytaskid > 0)
         allocate(cv%avgcrd(bcastdata(4)), stat = ierr)
         if (ierr /= 0) &
            NFE_OUT_OF_MEMORY
      end if ! .not. associated(cv%avgcrd)

      call mpi_bcast(cv%avgcrd, bcastdata(4), MPI_DOUBLE_PRECISION, &
                     0, pmemd_comm, ierr)
      nfe_assert(ierr == 0)
   else
      nullify(cv%avgcrd)
   end if ! bcastdata(4) > 0
  
   ! 
   ! cv % evec 
   ! 
    if (bcastdata(5) > 0) then
      if (.not. associated(cv%evec)) then
         nfe_assert(mytaskid > 0)
         allocate(cv%evec(bcastdata(5)), stat = ierr)
         if (ierr /= 0) &
            NFE_OUT_OF_MEMORY
      end if ! .not. associated(cv%evec)

      call mpi_bcast(cv%evec, bcastdata(5), MPI_DOUBLE_PRECISION, &
                     0, pmemd_comm, ierr)
      nfe_assert(ierr == 0)
   else
      nullify(cv%evec)
   end if ! bcastdata(5) > 0


   ! 
   ! cv % state_ref 
   ! 
    if (bcastdata(6) > 0) then
      if (.not. associated(cv%state_ref)) then
         nfe_assert(mytaskid > 0)
         allocate(cv%state_ref(bcastdata(6)), stat = ierr)
         if (ierr /= 0) &
            NFE_OUT_OF_MEMORY
      end if ! .not. associated(cv%state_ref)

      call mpi_bcast(cv%state_ref, bcastdata(6), MPI_INTEGER, &
                     0, pmemd_comm, ierr)
      nfe_assert(ierr == 0)
   else
      nullify(cv%state_ref)
   end if ! bcastdata(6) > 0


   ! 
   ! cv % state_pca
   ! 
    if (bcastdata(7) > 0) then
      if (.not. associated(cv%state_pca)) then
         nfe_assert(mytaskid > 0)
         allocate(cv%state_pca(bcastdata(7)), stat = ierr)
         if (ierr /= 0) &
            NFE_OUT_OF_MEMORY
      end if ! .not. associated(cv%state_pca)

      call mpi_bcast(cv%state_pca, bcastdata(7), MPI_INTEGER, &
                     0, pmemd_comm, ierr)
      nfe_assert(ierr == 0)
   else
      nullify(cv%state_pca)
   end if ! bcastdata(7) > 0
   
   ! 
   ! cv % ipca_to_i()
   ! 
    if (bcastdata(8) > 0) then
      if (.not. associated(cv%ipca_to_i)) then
         nfe_assert(mytaskid > 0)
         allocate(cv%ipca_to_i(bcastdata(8)), stat = ierr)
         if (ierr /= 0) &
            NFE_OUT_OF_MEMORY
      end if ! .not. associated(cv%state_pca)

      call mpi_bcast(cv%ipca_to_i, bcastdata(8), MPI_INTEGER, &
                     0, pmemd_comm, ierr)
      nfe_assert(ierr == 0)
   else
      nullify(cv%ipca_to_i)
   end if ! bcastdata(8) > 0

   !
   ! cv%q_index
   !
   if (bcastdata(9) > 0) then
     if (.not. associated(cv%q_index)) then
         nfe_assert(mytaskid > 0)
         allocate(cv%q_index, stat = ierr)
         if (ierr /= 0) &
            NFE_OUT_OF_MEMORY
      end if ! .not. associated(cv%q_index)
      call mpi_bcast(cv%q_index, 1, MPI_INTEGER, &
                     0, pmemd_comm, ierr)
      nfe_assert(ierr == 0)
   else
      nullify(cv%q_index)
   end if ! bcastdata(9) > 0 

   !
   ! cv%axis
   !
   if (bcastdata(10) > 0) then
      if (.not. associated(cv%axis)) then
         nfe_assert(mytaskid > 0)
         allocate(cv%axis(bcastdata(10)), stat = ierr)
         if (ierr /= 0) &
            NFE_OUT_OF_MEMORY
      end if ! .not. associated(cv%axis)

      call mpi_bcast(cv%axis, bcastdata(10), MPI_DOUBLE_PRECISION, &
                     0, pmemd_comm, ierr)
      nfe_assert(ierr == 0)
   else
      nullify(cv%axis)
   end if ! bcastdata(10) > 0



   
#endif /* MPI */

   !
   ! dispatch according to the type
   !

   select case(cv%type)
      case(COLVAR_ANGLE)
         call b_ANGLE(cv, cvno, amass)
      case(COLVAR_TORSION)
         call b_TORSION(cv, cvno, amass)
      case(COLVAR_DISTANCE)
         call b_DISTANCE(cv, cvno, amass)
      case(COLVAR_MULTI_RMSD)
         call b_MULTI_RMSD(cv, cvno, amass)
      case(COLVAR_R_OF_GYRATION)
         call b_R_OF_GYRATION(cv, cvno, amass)
      case(COLVAR_HANDEDNESS)
         call b_HANDEDNESS(cv, cvno, amass)
      case(COLVAR_N_OF_BONDS)
         call b_N_OF_BONDS(cv, cvno, amass)
      case(COLVAR_N_OF_STRUCTURES)
         call b_N_OF_STRUCTURES(cv, cvno, amass)
      case(COLVAR_LCOD)
         call b_LCOD(cv, cvno, amass)
      case(COLVAR_COS_OF_DIHEDRAL)
         call b_COS_OF_DIHEDRAL(cv, cvno, amass)
      case(COLVAR_COM_ANGLE)
         call b_COM_ANGLE(cv, cvno, amass)
      case(COLVAR_COM_TORSION)
         call b_COM_TORSION(cv, cvno, amass)
      case(COLVAR_COM_DISTANCE)
         call b_COM_DISTANCE(cv, cvno, amass)
      case(COLVAR_PCA)
         call b_PCA(cv, cvno, amass)
      case(COLVAR_SIN_OF_DIHEDRAL)
         call b_SIN_OF_DIHEDRAL(cv, cvno)
      case(COLVAR_PAIR_DIHEDRAL)
         call b_PAIR_DIHEDRAL(cv, cvno)
      case(COLVAR_PATTERN_DIHEDRAL)
         call b_PATTERN_DIHEDRAL(cv, cvno)
      case(COLVAR_DNA_RISE)
         call b_DNA_RISE(cv, cvno, amass)
      case(COLVAR_DF_COM_DISTANCE)
         call b_DF_COM_DISTANCE(cv, cvno, amass)
      case(COLVAR_ORIENTATION_ANGLE)
         call b_ORIENTATION_ANGLE(cv, cvno, amass)
      case(COLVAR_ORIENTATION_PROJ)
         call b_ORIENTATION_PROJ(cv, cvno, amass)
      case(COLVAR_SPINANGLE)
         call b_SPINANGLE(cv, cvno, amass)
      case(COLVAR_TILT)
         call b_TILT(cv, cvno, amass)
      case(COLVAR_QUATERNION0)
         call b_QUATERNION0(cv, cvno, amass)
      case(COLVAR_QUATERNION1)
         call b_QUATERNION1(cv, cvno, amass)
      case(COLVAR_QUATERNION2)
         call b_QUATERNION2(cv, cvno, amass)
      case(COLVAR_QUATERNION3)
         call b_QUATERNION3(cv, cvno, amass)
      case default
         nfe_assert_not_reached()
         continue
   end select

end subroutine colvar_bootstrap

!=============================================================================

subroutine colvar_nlread(cv_unit,cv)

  use nfe_lib_mod
  use file_io_mod, only : nmlsrc, amopen
  use pmemd_lib_mod, only : upper

  implicit none

  integer, intent(in)           :: cv_unit
  type(colvar_t), intent(inout) :: cv

  integer                    :: ifind, nsolut, error, i, tmp
  
  call nmlsrc('colvar', cv_unit, ifind)

  if (ifind .eq. 0) then
     write(unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR, &
                 'cannot find collective variables info...'
     call terminate()
  end if
  
  ! initiate 
  cv_min = 0.0
  cv_max = 5.0
  resolution = 0.0
  cv_type = ''
  cv_i(:) = 0
  cv_r(:) = 0.0
  cv_ni = 0
  cv_nr = 0
  path(:) = 0.0
  harm(:) = 0.0
  path_mode = ' '
  harm_mode = ' '
  anchor_position(:) = 0.0
  anchor_strength(:) = 0.0
  q_index = 1
  axis = [0.0, 0.0, 1.0]
  refcrd_file = 'inpcrd'

  read(cv_unit,nml=colvar,err=666)
  
  if (cv_ni.le.0) then
     write(unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR, &
           'number of elements in cv_i array must be positive!'
     call terminate()
  else
     allocate(cv%i(cv_ni), stat = error)
     if (error /= 0) &
     NFE_OUT_OF_MEMORY
     
     i = 1
     do while (i.le.cv_ni)
       cv%i(i) = cv_i(i)
       i = i + 1
     end do
  end if
  
  if (cv_nr.lt.0) then
     write(unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR, &
           'number of elements in cv_r array cannot be negative!'
     call terminate()
  end if
  
  if (cv_nr.gt.0) then
     allocate(cv%r(cv_nr), stat = error)
     if (error /= 0) &
     NFE_OUT_OF_MEMORY
     
     i = 1
     do while (i.le.cv_nr)
       cv%r(i) = cv_r(i)
       i = i + 1
     end do
  end if
  
  call upper(cv_type)
  if (cv_type == 'ANGLE') then
      cv%type = COLVAR_ANGLE
  else if (cv_type == 'TORSION') then
      cv%type = COLVAR_TORSION
  else if (cv_type == 'DISTANCE') then
      cv%type = COLVAR_DISTANCE
  else if (cv_type == 'MULTI_RMSD') then
      cv%type = COLVAR_MULTI_RMSD
  else if (cv_type == 'R_OF_GYRATION') then
      cv%type = COLVAR_R_OF_GYRATION
  else if (cv_type == 'HANDEDNESS') then
      cv%type = COLVAR_HANDEDNESS
  else if (cv_type == 'N_OF_BONDS') then
      cv%type = COLVAR_N_OF_BONDS
  else if (cv_type == 'N_OF_STRUCTURES') then
      cv%type = COLVAR_N_OF_STRUCTURES
  else if (cv_type == 'LCOD') then
      cv%type = COLVAR_LCOD
  else if (cv_type == 'COS_OF_DIHEDRAL') then
      cv%type = COLVAR_COS_OF_DIHEDRAL
  else if (cv_type == 'COM_ANGLE') then
      cv%type = COLVAR_COM_ANGLE
  else if (cv_type == 'COM_TORSION') then
      cv%type = COLVAR_COM_TORSION
  else if (cv_type == 'COM_DISTANCE') then
      cv%type = COLVAR_COM_DISTANCE
  else if (cv_type == 'PCA') then
      cv%type = COLVAR_PCA
      nsolut = pmemd_nsolut()
  else if (cv_type == 'SIN_OF_DIHEDRAL') then
      cv%type = COLVAR_SIN_OF_DIHEDRAL
  else if (cv_type == 'PAIR_DIHEDRAL') then
      cv%type = COLVAR_PAIR_DIHEDRAL
  else if (cv_type == 'PATTERN_DIHEDRAL') then
      cv%type = COLVAR_PATTERN_DIHEDRAL
  else if (cv_type == 'DNA_RISE') then
      cv%type = COLVAR_DNA_RISE
  else if (cv_type == 'DF_COM_DISTANCE') then
      cv%type = COLVAR_DF_COM_DISTANCE
  else if (cv_type == 'ORIENTATION_ANGLE') then
      cv%type = COLVAR_ORIENTATION_ANGLE
  else if (cv_type == 'ORIENTATION_PROJ') then
      cv%type = COLVAR_ORIENTATION_PROJ
  else if (cv_type == 'SPINANGLE') then
      cv%type = COLVAR_SPINANGLE
  else if (cv_type == 'TILT') then
      cv%type = COLVAR_TILT
  else if (cv_type == 'QUATERNION0') then
      cv%type = COLVAR_QUATERNION0
  else if (cv_type == 'QUATERNION1') then
      cv%type = COLVAR_QUATERNION1
  else if (cv_type == 'QUATERNION2') then
      cv%type = COLVAR_QUATERNION2
  else if (cv_type == 'QUATERNION3') then
      cv%type = COLVAR_QUATERNION3
  else
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
            NFE_ERROR, 'CV type ''', trim(cv_type), &
            ''' is not supported so far '
      call terminate()
  end if
    
  return
666 write(unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR,'Cannot read &colvar namelist!'
    call terminate()
end subroutine colvar_nlread

!=============================================================================
!================================ ANGLE ======================================
!=============================================================================
function v_ANGLE(cv, x) result(value)

#ifndef NFE_DISABLE_ASSERT
   use nfe_lib_mod
#endif /* NFE_DISABLE_ASSERT */

   use parallel_dat_mod

   implicit none

   double precision :: value

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: x(3,pmemd_natoms())

   nfe_assert(cv%type.eq.COLVAR_ANGLE)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).eq.3)

   nfe_assert(cv%i(1) > 0 .and. cv%i(1) <= pmemd_natoms())

   nfe_assert(cv%i(2) > 0 .and. cv%i(2) <= pmemd_natoms())

   nfe_assert(cv%i(3) > 0 .and. cv%i(3) <= pmemd_natoms())

   nfe_assert(cv%i(1) /= cv%i(2) .and. cv%i(1) /= cv%i(3) .and. cv%i(2) /= cv%i(3))

#ifdef MPI
   if (mytaskid.eq.0) then
#endif /* MPI */
      value = angle(x(1:3,cv%i(1)), x(1:3,cv%i(2)), x(1:3,cv%i(3)))
#ifdef MPI
   else
      value = ZERO
   end if
#endif /* MPI */

end function v_ANGLE

!=============================================================================

subroutine f_ANGLE(cv, x, fcv, f)

#ifndef NFE_DISABLE_ASSERT
   use nfe_lib_mod
#endif /* NFE_DISABLE_ASSERT */

   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: x(3,pmemd_natoms()), fcv

   double precision, intent(inout) :: f(3,pmemd_natoms())

   double precision :: d1(3), d2(3), d3(3)

#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_ANGLE)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i) == 3)

   nfe_assert(cv%i(1) > 0 .and. cv%i(1) <= pmemd_natoms())

   nfe_assert(cv%i(2) > 0 .and. cv%i(2) <= pmemd_natoms())

   nfe_assert(cv%i(3) > 0 .and. cv%i(3) <= pmemd_natoms())

   nfe_assert(cv%i(1) /= cv%i(2) .and. cv%i(1) /= cv%i(3) .and. cv%i(2) /= cv%i(3))

   NFE_MASTER_ONLY_BEGIN

   call angle_d(x(1:3,cv%i(1)), x(1:3,cv%i(2)), x(1:3,cv%i(3)), d1, d2, d3)

   f(1:3,cv%i(1)) = f(1:3,cv%i(1)) + fcv*d1
   f(1:3,cv%i(2)) = f(1:3,cv%i(2)) + fcv*d2
   f(1:3,cv%i(3)) = f(1:3,cv%i(3)) + fcv*d3

   NFE_MASTER_ONLY_END

#ifdef MPI
   call mpi_bcast(f, 3*pmemd_natoms(), MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

end subroutine f_ANGLE

!=============================================================================

subroutine b_ANGLE(cv, cvno, amass)

   use nfe_lib_mod
   
   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)

   nfe_assert(cv%type == COLVAR_ANGLE)
   call check_i(cv%i, cvno, 'ANGLE', 3)

end subroutine b_ANGLE

!=============================================================================

subroutine p_ANGLE(cv, lun)

#ifndef NFE_DISABLE_ASSERT
   use nfe_lib_mod
#endif /* NFE_DISABLE_ASSERT */

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_ANGLE)
   nfe_assert(associated(cv%i))

   call print_i(cv%i, lun)

end subroutine p_ANGLE

!=============================================================================
!================================ TORSION ====================================
!=============================================================================
function v_TORSION(cv, x) result(value)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   double precision :: value

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: x(*)

   integer :: a1, a2, a3, a4

   nfe_assert(cv%type == COLVAR_TORSION)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i) == 4)

   nfe_assert(cv%i(1) > 0 .and. cv%i(1) <= pmemd_natoms())
   a1 = 3*cv%i(1) - 2

   nfe_assert(cv%i(2) > 0 .and. cv%i(2) <= pmemd_natoms())
   a2 = 3*cv%i(2) - 2

   nfe_assert(cv%i(3) > 0 .and. cv%i(3) <= pmemd_natoms())
   a3 = 3*cv%i(3) - 2

   nfe_assert(cv%i(4) > 0 .and. cv%i(4) <= pmemd_natoms())
   a4 = 3*cv%i(4) - 2

   nfe_assert(a1 /= a2 .and. a1 /= a3 .and. a1 /= a4)
   nfe_assert(a2 /= a3 .and. a2 /= a4)
   nfe_assert(a3 /= a4)

#ifdef MPI
   if (mytaskid.eq.0) then
#endif /* MPI */
      value = torsion(x(a1:a1 + 2), x(a2:a2 + 2), x(a3:a3 + 2), x(a4:a4 + 2))
#ifdef MPI
   else
      value = ZERO
   end if
#endif /* MPI */

end function v_TORSION

!=============================================================================

subroutine f_TORSION(cv, x, fcv, f)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t),   intent(in) :: cv
   double precision, intent(in) :: x(*), fcv
   double precision, intent(inout) :: f(*)

   double precision :: d1(3), d2(3), d3(3), d4(3)

   integer :: a1, a2, a3, a4
#ifdef MPI
   integer :: error
#endif

   nfe_assert(cv%type == COLVAR_TORSION)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i) == 4)

   nfe_assert(cv%i(1) > 0 .and. cv%i(1) <= pmemd_natoms())
   a1 = 3*cv%i(1) - 2

   nfe_assert(cv%i(2) > 0 .and. cv%i(2) <= pmemd_natoms())
   a2 = 3*cv%i(2) - 2

   nfe_assert(cv%i(3) > 0 .and. cv%i(3) <= pmemd_natoms())
   a3 = 3*cv%i(3) - 2

   nfe_assert(cv%i(4) > 0 .and. cv%i(4) <= pmemd_natoms())
   a4 = 3*cv%i(4) - 2

   nfe_assert(a1 /= a2 .and. a1 /= a3 .and. a1 /= a4)
   nfe_assert(a2 /= a3 .and. a2 /= a4)
   nfe_assert(a3 /= a4)

   NFE_MASTER_ONLY_BEGIN
   call torsion_d(x(a1:a1 + 2), x(a2:a2 + 2), &
      x(a3:a3 + 2), x(a4:a4 + 2), d1, d2, d3, d4)

   f(a1:a1 + 2) = f(a1:a1 + 2) + fcv*d1
   f(a2:a2 + 2) = f(a2:a2 + 2) + fcv*d2
   f(a3:a3 + 2) = f(a3:a3 + 2) + fcv*d3
   f(a4:a4 + 2) = f(a4:a4 + 2) + fcv*d4

   NFE_MASTER_ONLY_END
#ifdef MPI
   call mpi_bcast(f(1:3*pmemd_natoms()), 3*pmemd_natoms(), MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif
end subroutine f_TORSION

!=============================================================================

subroutine b_TORSION(cv, cvno, amass)
   
   use nfe_lib_mod
   
   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)

   nfe_assert(cv%type == COLVAR_TORSION)
   call check_i(cv%i, cvno, 'TORSION', 4)

end subroutine b_TORSION

!=============================================================================

subroutine p_TORSION(cv, lun)

#ifndef NFE_DISABLE_ASSERT
   use nfe_lib_mod
#endif /* NFE_DISABLE_ASSERT */

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_TORSION)
   nfe_assert(associated(cv%i))

   call print_i(cv%i, lun)

end subroutine p_TORSION
!=============================================================================
!=============================== DISTANCE ====================================
!=============================================================================
logical function do_PBC(cv)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv

   do_PBC = .false.
   if (0.ne.pmemd_ntb()) then
      nfe_assert(associated(cv%r))
      if (.not.cv%r(1).gt.ZERO) then
         do_PBC = .false.
      else
         do_PBC = .true.
      end if
   end if

end function do_PBC

!=============================================================================

subroutine PBC_r2f(r, f)

   use pbc_mod, only : recip 

   implicit none

   double precision, intent(in) :: r(3)
   double precision, intent(out) :: f(3)

   f(1) = r(1)*recip(1,1) + r(2)*recip(2,1) + r(3)*recip(3,1)
   f(2) = r(1)*recip(1,2) + r(2)*recip(2,2) + r(3)*recip(3,2)
   f(3) = r(1)*recip(1,3) + r(2)*recip(2,3) + r(3)*recip(3,3)

   f(1) = f(1) - floor(f(1))
   f(2) = f(2) - floor(f(2))
   f(3) = f(3) - floor(f(3))

end subroutine PBC_r2f

!=============================================================================

subroutine PBC_f2r(f, r, t1, t2, t3)

   use pbc_mod, only : ucell

   implicit none

   double precision, intent(in) :: f(3)
   double precision, intent(out) :: r(3)

   integer, intent(in) :: t1, t2, t3

   r(1) = &
      (f(1) + t1)*ucell(1,1) + (f(2) + t2)*ucell(1,2) + (f(3) + t3)*ucell(1,3)
   r(2) = &
      (f(1) + t1)*ucell(2,1) + (f(2) + t2)*ucell(2,2) + (f(3) + t3)*ucell(2,3)
   r(3) = &
      (f(1) + t1)*ucell(3,1) + (f(2) + t2)*ucell(3,2) + (f(3) + t3)*ucell(3,3)

end subroutine PBC_f2r

!=============================================================================

function PBC_distance(r1, r2, d0) result(value)

   use nfe_lib_mod

   implicit none

   double precision :: value
   double precision, intent(in) :: r1(*), r2(*), d0

   double precision :: f1(3), f2(3), d(-1:1,-1:1,-1:1)
   double precision :: x1(3), x2(3), d_min

   integer :: i, j, k

   nfe_assert(d0.gt.ZERO)

   ! get the fractionals

   call PBC_r2f(r1(1:3), f1)
   call PBC_r2f(r2(1:3), f2)

   ! wrap r1 to the primary cell

   call PBC_f2r(f1, x1, 0, 0, 0)
   d_min = ZERO

   ! wrap r2 to the primary/neighboring cells

   do i = -1, 1
      do j = -1, 1
         do k = -1, 1
            call PBC_f2r(f2, x2, i, j, k)
            d(i,j,k) = distance(x1, x2)
            if (d(i,j,k).lt.d_min) &
               d_min = d(i,j,k)
         end do
      end do
   end do

   value = ZERO

   do i = -1, 1
      do j = -1, 1
         do k = -1, 1
            value = value + exp((d_min - d(i,j,k))/d0)
         end do
      end do
   end do

   value = d_min - d0*log(value)

end function PBC_distance

!=============================================================================

subroutine PBC_distance_d(r1, r2, d0, d1, d2)

   use nfe_lib_mod

   implicit none

   double precision, intent(in) :: r1(*), r2(*), d0
   double precision, intent(out) :: d1(*), d2(*)

   double precision :: f1(3), f2(3), d(-1:1,-1:1,-1:1)
   double precision :: x1(3), x2(3), d_min, sum, ttt

   integer :: i, j, k, l

   nfe_assert(d0.gt.ZERO)

   ! get the fractionals

   call PBC_r2f(r1(1:3), f1)
   call PBC_r2f(r2(1:3), f2)

   ! wrap r1 to the primary cell

   call PBC_f2r(f1, x1, 0, 0, 0)

   ! find d_min
   d_min = ZERO

   do i = -1, 1
      do j = -1, 1
         do k = -1, 1
            call PBC_f2r(f2, x2, i, j, k)
            d(i,j,k) = distance(x1, x2)
            if (d(i,j,k).lt.d_min) &
               d_min = d(i,j,k)
         end do
      end do
   end do

   do l = 1, 3
      d1(l) = ZERO
      d2(l) = ZERO
   end do

   sum = ZERO

   do i = -1, 1
      do j = -1, 1
         do k = -1, 1
            call PBC_f2r(f2, x2, i, j, k)

            ttt = exp((d_min - d(i,j,k))/d0)
            sum = sum + ttt
            ttt = ttt/d(i,j,k)

            do l = 1, 3
               d1(l) = d1(l) + (x1(l) - x2(l))*ttt
               d2(l) = d2(l) - (x1(l) - x2(l))*ttt
            end do
         end do
      end do
   end do

   do l = 1, 3
      d1(l) = d1(l)/sum
      d2(l) = d2(l)/sum
   end do

end subroutine PBC_distance_d

!=============================================================================

function v_DISTANCE(cv, x) result(value)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   double precision :: value

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: x(*)

   integer :: a1, a2

   nfe_assert(cv%type == COLVAR_DISTANCE)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i) == 2)

   nfe_assert(cv%i(1) > 0 .and. cv%i(1) <= pmemd_natoms())
   a1 = 3*cv%i(1) - 2

   nfe_assert(cv%i(2) > 0 .and. cv%i(2) <= pmemd_natoms())
   a2 = 3*cv%i(2) - 2

   nfe_assert(a1 /= a2)

#ifdef MPI
   if (mytaskid.eq.0) then
#endif /* MPI */
   if (do_PBC(cv)) then
      value = PBC_distance(x(a1:a1 + 2), x(a2:a2 + 2), cv%r(1))
   else
      value = distance(x(a1:a1 + 2), x(a2:a2 + 2))
   end if ! do_PBC(cv)
#ifdef MPI
   else
      value = ZERO
   end if
#endif /* MPI */

end function v_DISTANCE

!=============================================================================

subroutine f_DISTANCE(cv, x, fcv, f)

#ifndef NFE_DISABLE_ASSERT
   use nfe_lib_mod
#endif /* NFE_DISABLE_ASSERT */

   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: x(*), fcv

   double precision, intent(inout) :: f(*)

   double precision :: d1(3), d2(3)

   integer :: a1, a2
#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_DISTANCE)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i) == 2)

   nfe_assert(cv%i(1) > 0 .and. cv%i(1) <= pmemd_natoms())
   a1 = 3*cv%i(1) - 2

   nfe_assert(cv%i(2) > 0 .and. cv%i(2) <= pmemd_natoms())
   a2 = 3*cv%i(2) - 2

   nfe_assert(a1 /= a2)

   NFE_MASTER_ONLY_BEGIN

   if (do_PBC(cv)) then
      call PBC_distance_d(x(a1:a1 + 2), x(a2:a2 + 2), cv%r(1), d1, d2)
   else
      call distance_d(x(a1:a1 + 2), x(a2:a2 + 2), d1, d2)
   end if ! do_PBC(cv)

   f(a1:a1 + 2) = f(a1:a1 + 2) + fcv*d1
   f(a2:a2 + 2) = f(a2:a2 + 2) + fcv*d2

   NFE_MASTER_ONLY_END
#ifdef MPI
   call mpi_bcast(f(1:3*pmemd_natoms()), 3*pmemd_natoms(), MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

end subroutine f_DISTANCE

!=============================================================================

subroutine b_DISTANCE(cv, cvno, amass)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)
   integer :: error

   nfe_assert(cv%type == COLVAR_DISTANCE)
   call check_i(cv%i, cvno, 'DISTANCE', 2)

   if (0.eq.pmemd_ntb()) &
      return

   if (.not.associated(cv%r)) then
      allocate(cv%r(1), stat = error)
      if (error.ne.0) &
         NFE_OUT_OF_MEMORY
      cv%r(1) = ONE
   end if

   if (size(cv%r).ne.1) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(a,a,'//pfmt(cvno)//',a)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (DISTANCE) : unexpected number of reals (at most 1 is possible)'
      NFE_MASTER_ONLY_END
      call terminate()
   end if

end subroutine b_DISTANCE

!=============================================================================

subroutine p_DISTANCE(cv, lun)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_DISTANCE)
   nfe_assert(associated(cv%i))

   if (0.ne.pmemd_ntb()) then
      nfe_assert(associated(cv%r))
      nfe_assert(size(cv%r).eq.1)
      if (cv%r(1).gt.ZERO) then
         write (unit = lun, fmt = '(a,a,'//pfmt(cv%r(1), 3)//',a)') &
         NFE_INFO, '    (between the closest PBC images smoothed with d0 = ', &
         cv%r(1), ' A)'
      end if ! cv%r(1).gt.ZERO
   end if

   call print_i(cv%i, lun)

end subroutine p_DISTANCE

!=============================================================================
!============================== MULTI_RMSD ===================================
!=============================================================================
! input:
!
! cv%i = (a1, a2, a3, 0, b1, b2, b3, b4, 0, ..., c1, c2, c3, 0)
!
!     (a[1-3] - 1st group, b[1-4] - 2nd group, etc; an atom may
!      enter a few groups simultaneously; last zero is optional;
!      empty groups [e.g., 2+ zeros in a row] are not allowed)
!
! cv%r = (a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z, b1x, ...)
!
!        (reference coordinates without '0' sentinel(s))
!
! value = sqrt[(M_1*rmsd1^2 + ... + M_N*rmsdN^2)/(M_1 + ... + M_N)]
!         M_1 - mass of group 1, ..., M_N - mass of group N
!


!=============================================================================

function v_MULTI_RMSD(cv, x) result(value)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   double precision :: value

   type(colvar_t), intent(inout) :: cv

   double precision, intent(in) :: x(*)

   type(priv_t_MR), pointer :: priv
   integer :: g

   nfe_assert(cv%type == COLVAR_MULTI_RMSD)

   priv => get_priv_MR(cv)
   nfe_assert(associated(priv%groups))
   nfe_assert(size(priv%groups).gt.0)

   priv%value = ZERO

   do g = 1, size(priv%groups)
      call group_evaluate(cv, priv%groups(g), x)
      priv%value = priv%value + priv%groups(g)%rmsd2
   end do

   priv%value = sqrt(priv%value/priv%total_mass)
   
   value = priv%value
   
contains
subroutine group_evaluate(cv, grp, x)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in)    :: cv
   type(group_t_MR),  intent(inout) :: grp

   double precision, intent(in) :: x(*)

   double precision :: cm(3), cm_nrm, lambda
   integer :: i, n, a, a3

   ! compute the center of mass (of the moving atoms)
   cm = ZERO
   n = 1
   do a = grp%i0, grp%i1
      a3 = 3*cv%i(a)
      cm = cm + grp%mass(n)*x(a3 - 2:a3)
      n = n + 1
   end do

   cm = cm/grp%total_mass

   ! populate cm_crd && compute the "norm"
   cm_nrm = ZERO
   n = 1
   do a = grp%i0, grp%i1
      a3 = 3*(n - 1)
      do i = 1, 3
         grp%cm_crd(a3 + i) = x(3*(cv%i(a) - 1) + i) - cm(i)
         cm_nrm = cm_nrm + grp%mass(n)*grp%cm_crd(a3 + i)**2
      end do
      n = n + 1
   end do

   call rmsd_q(grp%i1 - grp%i0 + 1, grp%mass, &
         grp%cm_crd, grp%ref_crd, lambda, grp%quaternion)

   grp%rmsd2 = max(ZERO, ((grp%ref_nrm + cm_nrm) - 2*lambda))

end subroutine group_evaluate

end function v_MULTI_RMSD

!=============================================================================

!
! assumes that atom positions have not been
!   changed since last call to 'value()'
!

subroutine f_MULTI_RMSD(cv, x, fcv, f)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: x(*), fcv
   double precision, intent(inout) :: f(*)
   
   integer :: n, g, a, a3, f3
   double precision :: U(3,3), tmp
   type(priv_t_MR), pointer :: priv
#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_MULTI_RMSD)
   nfe_assert(associated(cv%i))

   priv => get_priv_MR(cv)

   nfe_assert(associated(priv))
   nfe_assert(associated(priv%groups))
   nfe_assert(size(priv%groups).gt.0)

   NFE_MASTER_ONLY_BEGIN
   do g = 1, size(priv%groups)

      call rmsd_q2u(priv%groups(g)%quaternion, U)
      tmp = fcv/max(priv%value, dble(0.000001))/priv%total_mass

      n = 1
      do a = priv%groups(g)%i0, priv%groups(g)%i1
         a3 = 3*n
         f3 = 3*cv%i(a)

         f(f3 - 2:f3) = f(f3 - 2:f3) &
            + tmp*priv%groups(g)%mass(n)*(priv%groups(g)%cm_crd(a3 - 2:a3) &
               - matmul(U, priv%groups(g)%ref_crd(a3 - 2:a3)))
         n = n + 1
      end do
   end do
   NFE_MASTER_ONLY_END
#ifdef MPI
   call mpi_bcast(f(1:3*pmemd_natoms()), 3*pmemd_natoms(), MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

end subroutine f_MULTI_RMSD

!=============================================================================

subroutine b_MULTI_RMSD(cv, cvno, amass)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)

   integer :: i, i0, n_atoms, n_groups, error

   type(priv_t_MR), pointer :: priv

   nfe_assert(cv%type == COLVAR_MULTI_RMSD)

   ! very basic checks
   if (.not. associated(cv%i)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (MULTI_RMSD) : no integers found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cv%i)

   if (.not. associated(cv%r)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_INFO, 'CV #', cvno, &
            ' (MULTI_RMSD) : no reals (reference coordinates) found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if

   ! count the groups (number of zeros in the cv%i array)
   n_atoms = 0
   n_groups = 0
   i0 = 1

   do i = 1, size(cv%i)
      if (cv%i(i) == 0) then
         n_groups = n_groups + 1
         if (i .eq. i0) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(i)//',a/)') &
                  NFE_ERROR, 'CV #', cvno, &
                  ' (MULTI_RMSD) : unexpected zero (integer #', i, ')'
            NFE_MASTER_ONLY_END
            call terminate()
         end if ! i .eq. i0
         i0 = i + 1
      else if (cv%i(i) .lt. 0) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (MULTI_RMSD) : negative integer'
         NFE_MASTER_ONLY_END
         call terminate()
      else
         n_atoms = n_atoms + 1
      end if
   end do

   if (size(cv%i) .gt. 0 .and. cv%i(size(cv%i)) .gt. 0) &
      n_groups = n_groups + 1

   nfe_assert(n_groups .gt. 0)

   if (size(cv%r) /= 3*n_atoms) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (MULTI_RMSD) : wrong number of reals'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! size(cv%r) /= 3*n_atoms

   ! allocate priv_t instance for this variable
   priv => new_priv_MR(cv)

   ! allocate/setup groups
   allocate(priv%groups(n_groups), stat = error)
   if (error /= 0) &
      NFE_OUT_OF_MEMORY

   i0 = 1 ! first atom
   n_atoms = 1
   n_groups = 0

   priv%total_mass = ZERO

   do i = 1, size(cv%i)
      if (cv%i(i) .eq. 0) then
         n_groups = n_groups + 1
         nfe_assert(n_groups .le. size(priv%groups))
         call group_bootstrap(priv%groups(n_groups), &
                              cv, cvno, amass, i0, i - 1, n_atoms)
         n_atoms = n_atoms + 3*(i - i0)
         i0 = i + 1
         priv%total_mass = priv%total_mass &
            + priv%groups(n_groups)%total_mass
      end if
   end do

   if (size(cv%i) .gt. 0 .and. cv%i(size(cv%i)) .gt. 0) then
      n_groups = n_groups + 1
      nfe_assert(n_groups .le. size(priv%groups))
      call group_bootstrap(priv%groups(n_groups), &
                           cv, cvno, amass, i0, size(cv%i), n_atoms)
      priv%total_mass = priv%total_mass &
         + priv%groups(n_groups)%total_mass
   end if

contains
subroutine group_bootstrap(grp, cv, cvno, amass, i0, i1, r0)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(group_t_MR),  intent(inout) :: grp

   type(colvar_t), intent(in) :: cv
   integer,        intent(in) :: cvno
   double precision,      intent(in) :: amass(*)
   integer, intent(in) :: i0, i1, r0

   integer :: a, b, n_atoms, error
   double precision :: cm(3)

   nfe_assert(cv%type == COLVAR_MULTI_RMSD)
   nfe_assert(associated(cv%i).and.associated(cv%r))
   nfe_assert(i0.gt.0.and.i0.le.size(cv%i))
   nfe_assert(i1.gt.0.and.i1.le.size(cv%i))
   nfe_assert(r0.gt.0.and.r0.lt.size(cv%r))

   ! basic checks
   if (i0 + 2 > i1) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, &
            fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
            (i0)//',a,'//pfmt(i1)//',a/)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (MULTI_RMSD) : too few integers in group (', i0, ':', i1, ')'
      NFE_MASTER_ONLY_END
      call terminate()
   end if

   do a = i0, i1
      nfe_assert(a .le. size(cv%i))
      if (cv%i(a) < 1 .or. cv%i(a) > pmemd_natoms()) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
               (cv%i(a))//',a,'//pfmt(pmemd_natoms())//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (MULTI_RMSD) : integer #', a, ' (', cv%i(a), &
               ') is out of range [1, ', pmemd_natoms(), ']'
         NFE_MASTER_ONLY_END
         call terminate()
      end if
   end do

   ! check for duplicates
   do a = i0, i1
      do b = a + 1, i1
         if (cv%i(a) == cv%i(b)) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
                  (b)//',a,'//pfmt(cv%i(a))//',a/)') &
                  NFE_ERROR, 'CV #', cvno, &
                  ' (MULTI_RMSD) : integers #', a, ' and #', b, &
                  ' are equal (', cv%i(a), ')'
            NFE_MASTER_ONLY_END
            call terminate()
         end if ! cv%i(a) == cv%i(b)
      end do
   end do

   ! allocate/setup

   n_atoms = i1 - i0 + 1

   grp%i0 = i0
   grp%i1 = i1
   grp%r0 = r0

   allocate(grp%mass(n_atoms), grp%cm_crd(3*n_atoms), &
            grp%ref_crd(3*n_atoms), stat = error)

   if (error /= 0) &
      NFE_OUT_OF_MEMORY

   cm = ZERO
   grp%total_mass = ZERO

   do a = 1, n_atoms
      grp%mass(a) = amass(cv%i(a + grp%i0 - 1))
      grp%total_mass = grp%total_mass + grp%mass(a)
      cm = cm + grp%mass(a)*cv%r(r0 + 3*(a - 1):r0 + 3*a - 1)
   end do

   cm = cm/grp%total_mass

   ! translate reference coordinates to CM frame
   grp%ref_nrm = ZERO
   do a = 1, n_atoms
      do b = 1, 3
         grp%ref_crd(3*(a - 1) + b) = cv%r(r0 + 3*(a - 1) + b - 1) - cm(b)
         grp%ref_nrm = grp%ref_nrm + grp%mass(a)*grp%ref_crd(3*(a - 1) + b)**2
      end do
   end do

end subroutine group_bootstrap

end subroutine b_MULTI_RMSD

!=============================================================================

subroutine c_MULTI_RMSD(cv)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(inout) :: cv

   integer :: g
   type(priv_t_MR), pointer :: priv

   nfe_assert(cv%type == COLVAR_MULTI_RMSD)
   nfe_assert(associated(cv%i).and.size(cv%i).gt.0)

   priv => get_priv_MR(cv)
   nfe_assert(associated(priv%groups))

   do g = 1, size(priv%groups)
      call group_finalize(priv%groups(g))
   end do

   deallocate(priv%groups)
   call del_priv_MR(cv)

contains
subroutine group_finalize(grp)

   use nfe_lib_mod

   implicit none

   type(group_t_MR), intent(inout) :: grp

   nfe_assert(associated(grp%mass))
   nfe_assert(associated(grp%cm_crd))
   nfe_assert(associated(grp%ref_crd))

   deallocate(grp%mass, grp%cm_crd, grp%ref_crd)

end subroutine group_finalize
end subroutine c_MULTI_RMSD

!=============================================================================

subroutine p_MULTI_RMSD(cv, lun)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   integer :: g
   type(priv_t_MR), pointer :: priv

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_MULTI_RMSD)
   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   priv => get_priv_MR(cv)
   nfe_assert(associated(priv%groups))

   do g = 1, size(priv%groups)
      write (unit = lun, fmt = '(a,a,'//pfmt(g)//',a)') &
         NFE_INFO, '<> group <> #', g, ':'
      call group_print(cv, priv%groups(g), lun)
   end do

contains
subroutine group_print(cv, grp, lun)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   type(group_t_MR),  intent(in) :: grp
   integer,        intent(in) :: lun

   integer :: a, c
   character(4) :: aname

   nfe_assert(associated(cv%i).and.associated(cv%r))
   nfe_assert(grp%i0.gt.0.and.grp%i0.le.size(cv%i))
   nfe_assert(grp%i1.gt.0.and.grp%i1.le.size(cv%i))

   write (unit = lun, fmt = '(a,a)', advance = 'NO') NFE_INFO, 'atoms = ('

   c = 1
   do a = grp%i0, grp%i1

      nfe_assert(cv%i(a) > 0 .and. cv%i(a) <= pmemd_natoms())
      aname = pmemd_atom_name(cv%i(a))

      write (unit = lun, fmt = '('//pfmt(cv%i(a))//',a,a,a)', advance = 'NO') &
         cv%i(a), ' [', trim(aname), ']'

      if (a == grp%i1) then
         write (unit = lun, fmt = '(a)') ')'
      else if (mod(c, 5) == 0) then
         write (unit = lun, fmt = '(a,/a,3x)', advance = 'NO') ',', NFE_INFO
      else
         write (unit = lun, fmt = '(a)', advance = 'NO') ', '
      end if

      c = c + 1
   end do

   write (unit = lun, fmt = '(a,a)') &
         NFE_INFO, 'reference coordinates :'

   c = 0
   do a = grp%i0, grp%i1
      nfe_assert(grp%r0 + c + 2 <= size(cv%r))
      write (unit = lun, fmt = '(a,3x,i5,a,f8.3,a,f8.3,a,f8.3)') &
         NFE_INFO, cv%i(a), ' : ', &
         cv%r(grp%r0 + c), ', ', cv%r(grp%r0 + c + 1), &
         ', ', cv%r(grp%r0 + c + 2)
      c = c + 3
   end do

end subroutine group_print
end subroutine p_MULTI_RMSD
!=============================================================================
function get_priv_MR(cv) result(ptr)

   NFE_USE_AFAILED

   implicit none

   type(colvar_t), intent(in) :: cv
   type(priv_t_MR), pointer :: ptr

   ptr => priv_list_MR
   do while (associated(ptr))
      if (ptr%tag.eq.cv%tag) &
         exit
      ptr => ptr%next
   end do

   nfe_assert(associated(ptr))

end function get_priv_MR
!=============================================================================
! allocates and appends to the list
function new_priv_MR(cv) result(ptr)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(inout) :: cv

   type(priv_t_MR), pointer :: ptr
   type(priv_t_MR), pointer :: head

   integer :: error

   allocate(ptr, stat = error)
   if (error.ne.0) &
      NFE_OUT_OF_MEMORY

   ptr%next => null()

   if (.not.associated(priv_list_MR)) then
      ptr%tag = 0
      priv_list_MR => ptr
   else
      head => priv_list_MR
      do while (associated(head%next))
         head => head%next
      end do
      ptr%tag = head%tag + 1
      head%next => ptr
   end if

   cv%tag = ptr%tag

end function new_priv_MR
!=============================================================================
! removes from the list and deallocates
subroutine del_priv_MR(cv)

   NFE_USE_AFAILED

   implicit none

   type(colvar_t), intent(in) :: cv

   type(priv_t_MR), pointer :: curr, prev

   nfe_assert(associated(priv_list_MR))

   curr => priv_list_MR
   if (curr%tag.eq.cv%tag) then
      priv_list_MR => curr%next
   else
      prev => curr
      curr => curr%next
      do while (associated(curr))
         if (curr%tag.eq.cv%tag) then
            prev%next => curr%next
            exit
         end if
         prev => curr
         curr => curr%next
      end do
   end if

   nfe_assert(associated(curr))
   nfe_assert(curr%tag.eq.cv%tag)

   deallocate(curr)

end subroutine del_priv_MR
!=============================================================================
!============================ R_OF_GYRATION ==================================
!=============================================================================
!
! cv%i = (i1, ..., iN) -- list of participating atoms
!
function v_R_OF_GYRATION(cv, x) result(value)

   use nfe_lib_mod

   implicit none

   double precision :: value

   type(colvar_t), intent(inout) :: cv

   double precision, intent(in) :: x(*)

   integer :: natoms, a, a3, i
   double precision :: tmp

   type(priv_t_RG), pointer :: priv

   nfe_assert(cv%type == COLVAR_R_OF_GYRATION)
   nfe_assert(associated(cv%i))

   natoms = size(cv%i)
   nfe_assert(natoms > 2)

   priv => get_priv_RG(cv)

   priv%Rg = ZERO
   priv%cm(1:3) = ZERO

   do a = 1, natoms
      a3 = 3*(cv%i(a) - 1)
      tmp = ZERO
      do i = 1, 3
         tmp = tmp + x(a3 + i)**2
         priv%cm(i) = priv%cm(i) + priv%weights(a)*x(a3 + i)
      end do
      priv%Rg = priv%Rg + priv%weights(a)*tmp
   end do

   priv%Rg = sqrt(priv%Rg - priv%cm(1)**2 - priv%cm(2)**2 - priv%cm(3)**2)
   value = priv%Rg

end function v_R_OF_GYRATION

!=============================================================================

!
! assumes that atom positions have not been
!   changed since last call to 'value()'
!

subroutine f_R_OF_GYRATION(cv, x, fcv, f)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: x(*), fcv
   double precision, intent(inout) :: f(*)

   integer :: natoms, a, a3
   double precision :: tmp

   type(priv_t_RG), pointer :: priv
#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_R_OF_GYRATION)
   nfe_assert(associated(cv%i))

   natoms = size(cv%i)
   nfe_assert(natoms > 2)

   priv => get_priv_RG(cv)
   nfe_assert(priv%Rg > ZERO)

   tmp = fcv/priv%Rg
  
   NFE_MASTER_ONLY_BEGIN
   do a = 1, natoms

      a3 = 3*cv%i(a)

      f(a3 - 2:a3) = f(a3 - 2:a3) &
         + tmp*priv%weights(a)*(x(a3 - 2:a3) - priv%cm(1:3))
   end do
   NFE_MASTER_ONLY_END
   
#ifdef MPI
   call mpi_bcast(f(1:3*pmemd_natoms()), 3*pmemd_natoms(), MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

end subroutine f_R_OF_GYRATION

!=============================================================================

subroutine b_R_OF_GYRATION(cv, cvno, amass)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)

   integer   :: natoms, a, error
   double precision :: total_mass

   type(priv_t_RG), pointer :: priv

   nfe_assert(cv%type == COLVAR_R_OF_GYRATION)

   natoms = size(cv%i)

   call check_i(cv%i, cvno, 'R_OF_GYRATION')
   if (.not. natoms > 2) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(a,a,'//pfmt(cvno)//',a)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (R_OF_GYRATION) : too few integers'
      NFE_MASTER_ONLY_END
      call terminate()
   end if

   total_mass = ZERO

   do a = 1, natoms
      total_mass = total_mass + amass(cv%i(a))
   end do

   nfe_assert(total_mass > ZERO)

   priv => new_priv_RG(cv)

   allocate(priv%weights(natoms), stat = error)
   if (error.ne.0) &
      NFE_OUT_OF_MEMORY

   do a = 1, natoms
      priv%weights(a) = amass(cv%i(a))/total_mass
   end do

end subroutine b_R_OF_GYRATION

!=============================================================================

subroutine p_R_OF_GYRATION(cv, lun)

   NFE_USE_AFAILED

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   nfe_assert(cv%type == COLVAR_R_OF_GYRATION)
   nfe_assert(associated(cv%i))

   call print_i(cv%i, lun)

end subroutine p_R_OF_GYRATION

!=============================================================================

subroutine c_R_OF_GYRATION(cv)

   NFE_USE_AFAILED

   implicit none

   type(colvar_t), intent(inout) :: cv

   type(priv_t_RG), pointer :: priv

   nfe_assert(cv%type.eq.COLVAR_R_OF_GYRATION)

   priv => get_priv_RG(cv)
   nfe_assert(associated(priv))

   deallocate(priv%weights)
   call del_priv_RG(cv)

end subroutine c_R_OF_GYRATION
!=============================================================================
function get_priv_RG(cv) result(ptr)

   NFE_USE_AFAILED

   implicit none

   type(colvar_t), intent(in) :: cv
   type(priv_t_RG), pointer :: ptr

   ptr => priv_list_RG
   do while (associated(ptr))
      if (ptr%tag.eq.cv%tag) &
         exit
      ptr => ptr%next
   end do

   nfe_assert(associated(ptr))

end function get_priv_RG
!=============================================================================
! allocates and appends to the list
function new_priv_RG(cv) result(ptr)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(inout) :: cv

   type(priv_t_RG), pointer :: ptr
   type(priv_t_RG), pointer :: head

   integer :: error

   allocate(ptr, stat = error)
   if (error.ne.0) &
      NFE_OUT_OF_MEMORY

   ptr%next => null()

   if (.not.associated(priv_list_RG)) then
      ptr%tag = 0
      priv_list_RG => ptr
   else
      head => priv_list_RG
      do while (associated(head%next))
         head => head%next
      end do
      ptr%tag = head%tag + 1
      head%next => ptr
   end if

   cv%tag = ptr%tag

end function new_priv_RG
!=============================================================================
! removes from the list and deallocates
subroutine del_priv_RG(cv)

   NFE_USE_AFAILED

   implicit none

   type(colvar_t), intent(in) :: cv

   type(priv_t_RG), pointer :: curr, prev

   nfe_assert(associated(priv_list_RG))

   curr => priv_list_RG
   if (curr%tag.eq.cv%tag) then
      priv_list_RG => curr%next
   else
      prev => curr
      curr => curr%next
      do while (associated(curr))
         if (curr%tag.eq.cv%tag) then
            prev%next => curr%next
            exit
         end if
         prev => curr
         curr => curr%next
      end do
   end if

   nfe_assert(associated(curr))
   nfe_assert(curr%tag.eq.cv%tag)

   deallocate(curr)

end subroutine del_priv_RG
!=============================================================================
!============================== HANDEDNESS ===================================
!=============================================================================
!
! cv%i = (i1, ..., iN) -- list of participating atoms
!
function v_HANDEDNESS(cv, x) result(value)

   use nfe_lib_mod

   implicit none

   double precision :: value

   type(colvar_t), intent(inout) :: cv

   double precision, intent(in) :: x(*)

   integer :: natoms, a, j1, j2, j3, j4

   nfe_assert(cv%type.eq.COLVAR_HANDEDNESS)

   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   natoms = size(cv%i)
   nfe_assert(natoms.gt.3)

   value = ZERO

   do a = 1, natoms - 3

      j1 = 3*cv%i(a + 0)
      j2 = 3*cv%i(a + 1)
      j3 = 3*cv%i(a + 2)
      j4 = 3*cv%i(a + 3)

      value = value + value4(cv%r(1), &
         x(j1 - 2:j1), x(j2 - 2:j2), x(j3 - 2:j3), x(j4 - 2:j4))

   end do

end function v_HANDEDNESS

!=============================================================================

subroutine f_HANDEDNESS(cv, x, fcv, f)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: x(*), fcv
   double precision, intent(inout) :: f(*)

   integer :: natoms, a, j1, j2, j3, j4
   double precision :: d1(3), d2(3), d3(3), d4(3)
#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type.eq.COLVAR_HANDEDNESS)

   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   natoms = size(cv%i)
   nfe_assert(natoms.gt.3)

   NFE_MASTER_ONLY_BEGIN

   do a = 1, natoms - 3

      j1 = 3*cv%i(a + 0)
      j2 = 3*cv%i(a + 1)
      j3 = 3*cv%i(a + 2)
      j4 = 3*cv%i(a + 3)

      call derivative4(cv%r(1), &
         x(j1 - 2:j1), x(j2 - 2:j2), x(j3 - 2:j3), x(j4 - 2:j4), &
         d1, d2, d3, d4)

      f(j1 - 2:j1) = f(j1 - 2:j1) + fcv*d1
      f(j2 - 2:j2) = f(j2 - 2:j2) + fcv*d2
      f(j3 - 2:j3) = f(j3 - 2:j3) + fcv*d3
      f(j4 - 2:j4) = f(j4 - 2:j4) + fcv*d4

   end do

   NFE_MASTER_ONLY_END
#ifdef MPI
   call mpi_bcast(f(1:3*pmemd_natoms()), 3*pmemd_natoms(), MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

end subroutine f_HANDEDNESS

!=============================================================================

subroutine b_HANDEDNESS(cv, cvno, amass)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)

   integer :: natoms, error

   nfe_assert(cv%type == COLVAR_HANDEDNESS)

   natoms = size(cv%i)

   call check_i(cv%i, cvno, 'HANDEDNESS')
   if (.not.natoms.gt.3) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(a,a,'//pfmt(cvno)//',a)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (HANDEDNESS) : too few integers'
      NFE_MASTER_ONLY_END
      call terminate()
   end if

   if (.not.associated(cv%r)) then
      allocate(cv%r(1), stat = error)
      if (error.ne.0) &
         NFE_OUT_OF_MEMORY
      cv%r(1) = ZERO
   end if

   if (size(cv%r).ne.1) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(a,a,'//pfmt(cvno)//',a)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (HANDEDNESS) : unexpected number of reals (just 1 is needed)'
      NFE_MASTER_ONLY_END
      call terminate()
   end if

   cv%r(1) = max(ZERO, cv%r(1))
   cv%r(1) = min(ONE,  cv%r(1))

end subroutine b_HANDEDNESS

!=============================================================================

subroutine p_HANDEDNESS(cv, lun)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   nfe_assert(cv%type.eq.COLVAR_HANDEDNESS)

   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   write (unit = lun, fmt = '(a,a,'//pfmt(cv%r(1), 3)//')') &
      NFE_INFO, '      w = ', cv%r(1)
   call print_i(cv%i, lun)

end subroutine p_HANDEDNESS

!=============================================================================

!
! ** value4() & derivative4() are maple-generated & postprocessed by hands **
!
! dot3 := proc(u, v)
!     u[1]*v[1] + u[2]*v[2] + u[3]*v[3]
! end proc;
! 
! cross3 := proc(u, v)
!     [
!         u[2]*v[3] - u[3]*v[2],
!         u[3]*v[1] - u[1]*v[3],
!         u[1]*v[2] - v[1]*u[2]
!     ]
! end proc;
! 
! 
! r1 := array(1..3, [r1x, r1y, r1z]);
! r2 := array(1..3, [r2x, r2y, r2z]);
! r3 := array(1..3, [r3x, r3y, r3z]);
! r4 := array(1..3, [r4x, r4y, r4z]);
! 
! u1 := evalm(r2 - r1);
! u2 := evalm(r4 - r3);
! u3 := evalm((1 - w)*(r3 - r2) + w*(r4 - r1));
! 
! u1n := dot3(u1, u1) + TINY;
! u2n := dot3(u2, u2) + TINY;
! u3n := dot3(u3, u3) + TINY;
! 
! dofs := [r1x, r1y, r1z, r2x, r2y, r2z, r3x, r3y, r3z, r4x, r4y, r4z];
! 
! handedness := dot3(u3, cross3(u1, u2))/sqrt(u1n*u2n*u3n);
! handedness := simplify(handedness);
! 
! handedness_v := codegen[makeproc](handedness, dofs);
! handedness_v := codegen[optimize](handedness_v, tryhard);
! 
! handedness_d := codegen[GRADIENT]
!    (codegen[split](handedness_v, dofs), dofs,
!     function_value = false, mode = reverse);
! 
! handedness_v := codegen[packargs](handedness_v, [r1x, r1y, r1z], `x1`);
! handedness_v := codegen[packargs](handedness_v, [r2x, r2y, r2z], `x2`);
! handedness_v := codegen[packargs](handedness_v, [r3x, r3y, r3z], `x3`);
! handedness_v := codegen[packargs](handedness_v, [r4x, r4y, r4z], `x4`);
! 
! handedness_d := codegen[optimize](handedness_d, tryhard);
! handedness_d := codegen[packargs](handedness_d, [r1x, r1y, r1z], `x1`);
! handedness_d := codegen[packargs](handedness_d, [r2x, r2y, r2z], `x2`);
! handedness_d := codegen[packargs](handedness_d, [r3x, r3y, r3z], `x3`);
! handedness_d := codegen[packargs](handedness_d, [r4x, r4y, r4z], `x4`);
! 
! CodeGeneration[Fortran](handedness_v, deducetypes = false,
!                         optimize, defaulttype = numeric,
! 						output = `ncsu-handedness.maple.f`);
! 
! CodeGeneration[Fortran](handedness_d, deducetypes = false,
!                         optimize, defaulttype = numeric,
! 						output = `ncsu-handedness.maple.f`);
! 

double precision function value4 (w, x1, x2, x3, x4)

   use nfe_lib_mod, only : TINY
   
   implicit none

   double precision, intent(in) :: w
   double precision, intent(in) :: x1(3)
   double precision, intent(in) :: x2(3)
   double precision, intent(in) :: x3(3)
   double precision, intent(in) :: x4(3)

   double precision ::  t13
   double precision ::  t14
   double precision ::  t119
   double precision ::  t1
   double precision ::  t23
   double precision ::  t22
   double precision ::  t34
   double precision ::  t5
   double precision ::  t6
   double precision ::  t12
   double precision ::  t40
   double precision ::  t18
   double precision ::  t7
   double precision ::  t55
   double precision ::  t32
   double precision ::  t8
   double precision ::  t25
   double precision ::  t30
   double precision ::  t29
   double precision ::  t43
   double precision ::  t35
   double precision ::  t21
   double precision ::  t2
   double precision ::  t9
   double precision ::  t27
   double precision ::  t36
   double precision ::  t24
   double precision ::  t28
   double precision ::  t26
   double precision ::  t3
   double precision ::  t15
   double precision ::  t31
   double precision ::  t33
   double precision ::  t17
   double precision ::  t16
   double precision ::  t39
   double precision ::  t52
   double precision ::  t41
   double precision ::  t42
   double precision ::  t10
   double precision ::  t4
   double precision ::  t11
   double precision ::  t20
   double precision ::  t19
   double precision ::  t37
   double precision ::  t38

   t1 = x2(3)
   t7 = t1 ** 2
   t2 = x2(2)
   t9 = t2 ** 2
   t3 = x2(1)
   t11 = t3 ** 2
   t32 = t7 + t9 + t11
   t4 = x3(3)
   t8 = t4 ** 2
   t5 = x3(2)
   t10 = t5 ** 2
   t6 = x3(1)
   t12 = t6 ** 2
   t31 = t8 + t10 + t12
   t13 = 0.2D1 * t2
   t30 = -t13
   t33 = 0.2D1 * t3
   t29 = -t33
   t34 = 0.2D1 * t4
   t28 = -t34
   t35 = x4(1)
   t36 = x1(1)
   t27 = t35 - t36
   t37 = x4(2)
   t38 = x1(2)
   t26 = t37 - t38
   t39 = x4(3)
   t40 = x1(3)
   t25 = t39 - t40
   t24 = -t35 + t6
   t23 = -t37 + t5
   t22 = t4 - t39
   t21 = t33
   t20 = t13
   t19 = 0.2D1 * t1
   t41 = t40 ** 2
   t42 = t38 ** 2
   t43 = t36 ** 2
   t18 = t41 + t42 + t43 + t32
   t17 = -t6 * t37 + t5 * t35
   t16 = -t5 * t39 + t4 * t37
   t15 = t6 * t39 - t4 * t35
   t52 = 0.2D1 * t5
   t55 = 0.2D1 * t6
   t14 = t31 + (t28 + t39) * t39 + (-t52 + t37) * t37 + (-t55 + t35) * t35
   t119 = sqrt((t36 * t29 + t38 * t30 - 0.2D1 * t1 * t40 + TINY + t18) &
        * (TINY + t14) * (t1 * t28 + TINY + t6 * t29 + t5 * t30 + (-0.2D1 &
        * t12 - 0.2D1 * t11 - 0.2D1 * t7 - 0.2D1 * t9 - 0.2D1 * t10 - &
        0.2D1 * t8 + (t52 - t26) * t20 + (t55 - t27) * t21 + (t34 - t25) * &
       t19 + (t14 - t24 * t21 - t23 * t20 - t22 * t19 + t18) * w) * w + &
       0.2D1 * (t27 * t6 + t25 * t4 + t26 * t5 + ((-t3 + t24) * t36 + (-t1 &
       + t22) * t40 + (-t2 + t23) * t38) * w) * w + t31 + t32))

   value4 = (t17 * t1 + t15 * t2 + t16 * t3 + (-t24 * t2 + t23 * t3 - t17) &
      * t40 + (t24 * t1 - t22 * t3 - t15) * t38 + (-t23 * t1 + t22 * &
      t2 - t16) * t36) / t119

end function value4

subroutine derivative4(w, x1, x2, x3, x4, d1, d2, d3, d4)

   use nfe_lib_mod, only : TINY
   
   implicit none

   double precision, intent(in) :: w
   double precision, intent(in) :: x1(3)
   double precision, intent(in) :: x2(3)
   double precision, intent(in) :: x3(3)
   double precision, intent(in) :: x4(3)

   double precision, intent(out) :: d1(3)
   double precision, intent(out) :: d2(3)
   double precision, intent(out) :: d3(3)
   double precision, intent(out) :: d4(3)

   double precision ::  t216
   double precision ::  t41
   double precision ::  t42
   double precision ::  t43
   double precision ::  t40
   double precision ::  t128
   double precision ::  t184
   double precision ::  t191
   double precision ::  t20
   double precision ::  t93
   double precision ::  t307
   double precision ::  t19
   double precision ::  t18
   double precision ::  s1
   double precision ::  t17
   double precision ::  t106
   double precision ::  t83
   double precision ::  t85
   double precision ::  t124
   double precision ::  df(27)
   double precision ::  t55
   double precision ::  t16
   double precision ::  t15
   double precision ::  t38
   double precision ::  t39
   double precision ::  t14
   double precision ::  t110
   double precision ::  t100
   double precision ::  t59
   double precision ::  t60
   double precision ::  t87
   double precision ::  t52
   double precision ::  t137
   double precision ::  t1
   double precision ::  t2
   double precision ::  t114
   double precision ::  t117
   double precision ::  t125
   double precision ::  t122
   double precision ::  t146
   double precision ::  t295
   double precision ::  t13
   double precision ::  t33
   double precision ::  t22
   double precision ::  t21
   double precision ::  t34
   double precision ::  t35
   double precision ::  t321
   double precision ::  t3
   double precision ::  t174
   double precision ::  t36
   double precision ::  t37
   double precision ::  t95
   double precision ::  s0
   double precision ::  t6
   double precision ::  t7
   double precision ::  t123
   double precision ::  t9
   double precision ::  t11
   double precision ::  t4
   double precision ::  t5
   double precision ::  t32
   double precision ::  t8
   double precision ::  t10
   double precision ::  t12
   double precision ::  t31
   double precision ::  t30
   double precision ::  t29
   double precision ::  t28
   double precision ::  t179
   double precision ::  t27
   double precision ::  t237
   double precision ::  t201
   double precision ::  t204
   double precision ::  t206
   double precision ::  t79
   double precision ::  t291
   double precision ::  t293
   double precision ::  t26
   double precision ::  t25
   double precision ::  t58
   double precision ::  t24
   double precision ::  t23
   double precision ::  t226
   double precision ::  t310

   t1 = x2(3)
   t7 = t1 ** 2
   t2 = x2(2)
   t9 = t2 ** 2
   t3 = x2(1)
   t11 = t3 ** 2
   t32 = t7 + t9 + t11
   t4 = x3(3)
   t8 = t4 ** 2
   t5 = x3(2)
   t10 = t5 ** 2
   t6 = x3(1)
   t12 = t6 ** 2
   t31 = t8 + t10 + t12
   t13 = 0.2D1 * t2
   t30 = -t13
   t33 = 0.2D1 * t3
   t29 = -t33
   t34 = 0.2D1 * t4
   t28 = -t34
   t35 = x4(1)
   t36 = x1(1)
   t27 = t35 - t36
   t37 = x4(2)
   t38 = x1(2)
   t26 = t37 - t38
   t39 = x4(3)
   t40 = x1(3)
   t25 = t39 - t40
   t24 = -t35 + t6
   t23 = -t37 + t5
   t22 = t4 - t39
   t21 = t33
   t20 = t13
   t19 = 0.2D1 * t1
   t41 = t40 ** 2
   t42 = t38 ** 2
   t43 = t36 ** 2
   t18 = t41 + t42 + t43 + t32
   t17 = -t6 * t37 + t5 * t35
   t16 = -t5 * t39 + t4 * t37
   t15 = t6 * t39 - t4 * t35
   t52 = 0.2D1 * t5
   t55 = 0.2D1 * t6
   t14 = t31 + (t28 + t39) * t39 + (-t52 + t37) * t37 + (-t55 + t35) * t35
   t58 = t1 * t28
   t59 = t6 * t29
   t60 = t5 * t30
   t79 = (-0.2D1 * t12 - 0.2D1 * t11 - 0.2D1 * t7 - 0.2D1 * t9 - 0.2D1 * t10 &
   - 0.2D1 * t8 + (t52 - t26) * t20 + (t55 - t27) * t21 + &
   (t34 - t25) * t19 + (t14 - t24 * t21 - t23 * t20 - t22 * t19 + t18)* w) * w
   t83 = -t3 + t24
   t85 = -t1 + t22
   t87 = -t2 + t23
   t93 = 0.2D1 * (t27 * t6 + t25 * t4 + t26 * t5 + &
      (t83 * t36 + t85 * t40 + t87 * t38) * w) * w
   t95 = TINY + t14
   s0 = (t58 + TINY + t59 + t60 + t79 + t93 + t31 + t32) * t95
   t100 = t36 * t29 + t38 * t30 - 0.2D1 * t1 * t40 + TINY + t18
   s1 = t100 * s0
   t106 = -t24 * t2 + t23 * t3 - t17
   t110 = t24 * t1 - t22 * t3 - t15
   t114 = -t23 * t1 + t22 * t2 - t16
   t117 = sqrt(s1)
    df(27) = -(t17 * t1 + t15 * t2 + t16 * t3 + t106 * t40 + t110 * t38 &
       + t114 * t36) / t117 / s1 / 0.2D1
   t122 = df(27)
   df(26) = t122 * t100
   t123 = df(26)
   t124 = w ** 2
   t125 = t124 * t95
   df(25) = t123 * (t125 + t58 + TINY + t59 + t60 + t79 + t93 + t31 + t32)
   t128 = 0.1D1 / t117
   df(24) = (t2 - t38) * t128
   df(23) = (t3 - t36) * t128
   df(22) = (t1 - t40) * t128
   df(21) = t122 * s0 + t123 * t124 * t95
   t137 = w * t95
   df(20) = t123 * (t34 - t25 - t22 * w) * t137
   df(19) = t123 * (t52 - t26 - t23 * w) * t137
   df(18) = t123 * (t55 - t27 - t24 * w) * t137
   t146 = 0.2D1 * t40 * t124
   df(17) = t123 * (-t19 * t124 + t146) * t95 + (-t3 * t38 + t2 * t36) * t128
   df(16) = t123 * (-t20 * t124 + 0.2D1 * t38 * t124) * t95 + &
     (t3*t40 - t1 * t36) * t128
   df(15) = t123 * (-t21 * t124 + 0.2D1 * t36 * t124) * t95 + &
     (-t2 * t40 + t1 * t38) * t128
   t174 = t19 * w
   df(14) = t123 * (-t174 + 0.2D1 * t4 * w) * t95
   t179 = t20 * w
   df(13) = t123 * (-t179 + 0.2D1 * t5 * w) * t95
   t184 = t21 * w
   df(12) = t123 * (-t184 + 0.2D1 * t6 * w) * t95
   t191 = df(25)
   df(11) = t123 * t1 * t95 + t191 * t39
   df(10) = t122 * t36 * s0 + t123 * t6 * t95
   df(9) = t122 * t38 * s0 + t123 * t5 * t95
   t201 = t123 * t95
   df(8) = t201 + t191
   t204 = 0.2D1 * t123 * w * t95
   df(7) = -t204 + df(8)
   df(6) = df(7)
   df(5) = df(6)
   t206 = df(21)
   df(4) = t201 + t206
   df(3) = -t204 + df(4)
   df(2) = df(3)
   df(1) = df(2)
   t216 = df(12)
   t226 = df(13)
   t237 = df(14)
   t291 = df(24)
   t293 = df(22)
   t295 = df(15)
   t307 = df(23)
   t310 = df(16)
   t321 = df(17)
   d1(1) = t114 * t128 + t122 * t29 * s0 + 0.2D1 * t123 * t83 * t125 &
   + 0.2D1 * t206 * t36 - t216
   d1(2) = t110 * t128 + t122 * t30 * s0 + 0.2D1 * t123 * t87 * t125 &
   + 0.2D1 * t206 * t38 - t226
   d1(3) = t106 * t128 - 0.2D1 * t122 * t1 * s0 + 0.2D1 * t123 * t85 * t125 &
    + 0.2D1 * t206 * t40 - t237
   d2(1) = (t16 + t23 * t40 - t22 * t38) * t128 - 0.2D1 * t123 * t36 * t125 &
   + 0.2D1 * df(18) - 0.2D1 * df(10) + 0.2D1 * df(2) * t3
   d2(2) = (t15 - t24 * t40 + t22 * t36) * t128 - 0.2D1 * t123 * t38 * t125 &
   + 0.2D1 * df(19) - 0.2D1 * df(9) + 0.2D1 * df(1) * t2
   d2(3) = (t17 + t24 * t38 - t23 * t36) * t128 - 0.2D1 * t122 * &
     t40 * s0 + t123 * (t28 - t146) * t95 + 0.2D1 * df(20) + 0.2D1 * df(1) * t1
   d3(1) = t123 * (t29 + 0.2D1 * t184 + 0.2D1 * t27 * w) * t95 - &
     0.2D1 * t191 * t35 + t291 * t39 - t293 * t37 + t295 + 0.2D1 * df(6) * t6
   d3(2) = t123 * (t30 + 0.2D1 * t179 + 0.2D1 * t26 * w) * t95 - &
     0.2D1 * t191 * t37 - t307 * t39 + t293 * t35 + t310 + 0.2D1 * df(5) * t5
   d3(3) = t123 * (0.2D1 * t174 + 0.2D1 * t25 * w) * t95 - t291 &
      * t35 + t307 * t37 + t321 - 0.2D1 * df(11) + 0.2D1 * df(5) * t4
   d4(1) = -0.2D1 * t191 * t24 - t291 * t4 + t293 * t5 - t295 + t216
   d4(2) = -0.2D1 * t191 * t23 + t307 * t4 - t293 * t6 - t310 + t226
   d4(3) = t191 * (0.2D1 * t39 + t28) + t291 * t6 - t307 * t5 - t321 + t237

end subroutine derivative4
!=============================================================================
!============================= N_OF_BONDS ====================================
!=============================================================================
double precision pure function f6_12_value(x)

   use nfe_lib_mod, only : ONE

   implicit none

   double precision, intent(in) :: x

   f6_12_value = ONE/(ONE + x**3)

end function f6_12_value

!=============================================================================

double precision pure function f6_12_derivative(x)

   use nfe_lib_mod, only : ONE

   implicit none

   double precision, parameter :: MINUS_THREE = -3.000000000000000000000D0 ! dble(-3)

   double precision, intent(in) :: x

   double precision :: x2, tmp

   x2 = x*x
   tmp = ONE + x*x2

   f6_12_derivative = MINUS_THREE*x2/(tmp*tmp)

end function f6_12_derivative

!=============================================================================

subroutine partition(cv, first, last)

   NFE_USE_AFAILED

   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(out) :: first, last

#  ifdef MPI
   integer :: tmp
#  endif /* MPI */

   nfe_assert(cv%type == COLVAR_N_OF_BONDS)
   nfe_assert(associated(cv%i))
   nfe_assert(mod(size(cv%i), 2).eq.0)

#  ifdef MPI
   tmp = (size(cv%i)/2)/numtasks
   if (tmp.gt.0) then
      if (mytaskid.eq.(numtasks - 1)) then
         first = 2*tmp*mytaskid + 1
         last = size(cv%i) - 1
      else
         first = 2*tmp*mytaskid + 1
         last = 2*(mytaskid + 1)*tmp - 1
      end if
   else
      if (mytaskid.eq.(numtasks - 1)) then
         first = 1
         last = size(cv%i) - 1
      else
         first = 1
         last = 0
      end if
   end if
#  else
   first = 1
   last = size(cv%i) - 1
#  endif /* MPI */

end subroutine partition

!=============================================================================
!F Pan & H Macdermott-Opeskin
function v_N_OF_BONDS(cv, x) result(value)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   double precision :: value

   type(colvar_t), intent(inout) :: cv

   double precision, intent(in) :: x(*)

#  ifdef MPI
   integer :: error
   double precision :: accu
#  endif /* MPI */

   integer :: first, last
   integer :: i, i3, j3

   double precision :: r2
   double precision, parameter :: d0 = ONE

   nfe_assert(cv%type == COLVAR_N_OF_BONDS)

   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   nfe_assert(mod(size(cv%i), 2).eq.0)

   value = ZERO

   call partition(cv, first, last)
   
   do i = first, last, 2

      i3 = 3*cv%i(i) - 2
      j3 = 3*cv%i(i + 1) - 2

      !r2 = (x(i3 + 0) - x(j3 + 0))**2 &
      !   + (x(i3 + 1) - x(j3 + 1))**2 &
      !   + (x(i3 + 2) - x(j3 + 2))**2

      ! correct this by using PBC distance
      r2 = (PBC_distance(x(i3:i3+2), x(j3:j3+2), d0))**2

      r2 = r2/cv%r(1)**2

      value = value + f6_12_value(r2)

   end do

#  ifdef MPI
   call mpi_reduce(value, accu, 1, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
   value = accu
#  endif /* MPI */

end function v_N_OF_BONDS

!=============================================================================
!F Pan & H Macdermott-Opeskin
subroutine f_N_OF_BONDS(cv, x, fcv, f)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: x(*), fcv
   double precision, intent(inout) :: f(*)

   integer   :: i, i3, j3, k
   double precision :: r2, tmp, r
   double precision :: dxf1(3), dxf2(3), dx1(3), dx2(3)
   double precision, parameter :: d0 = ONE
#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_N_OF_BONDS)

   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   NFE_MASTER_ONLY_BEGIN
   do i = 1, size(cv%i), 2

      i3 = 3*cv%i(i) - 2
      j3 = 3*cv%i(i + 1) - 2
      !find the distance
      r = PBC_distance(x(i3:i3 + 2), x(j3:j3 + 2), d0)
      !find the unit displacement vectors
      call PBC_distance_d(x(i3:i3 + 2), x(j3:j3 + 2), d0, dx1, dx2)
      !calculate the displacement i to j 
      do k = 1,3
        dxf1(k) = r*dx1(k)
        dxf2(k) = r*dx2(k)
      end do

      !r2
      r2 = (r**2)/(cv%r(1)**2)

      tmp = (2*fcv/cv%r(1)**2)*f6_12_derivative(r2)

      dxf1(1) = tmp*dxf1(1)
      dxf2(1) = tmp*dxf2(1)
      f(i3 + 0) = f(i3 + 0) + dxf1(1)
      f(j3 + 0) = f(j3 + 0) + dxf2(1)

      dxf1(2) = tmp*dxf1(2)
      dxf2(2) = tmp*dxf2(2)
      f(i3 + 1) = f(i3 + 1) + dxf1(2)
      f(j3 + 1) = f(j3 + 1) + dxf2(2)

      dxf1(3) = tmp*dxf1(3)
      dxf2(3) = tmp*dxf2(3)
      f(i3 + 2) = f(i3 + 2) + dxf1(3)
      f(j3 + 2) = f(j3 + 2) + dxf2(3)

   end do
   NFE_MASTER_ONLY_END
   
#ifdef MPI
   call mpi_bcast(f(1:3*pmemd_natoms()), 3*pmemd_natoms(), MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */
end subroutine f_N_OF_BONDS

!=============================================================================

subroutine b_N_OF_BONDS(cv, cvno, amass)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)
   integer :: a

   nfe_assert(cv%type == COLVAR_N_OF_BONDS)

   if (.not.associated(cv%i)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (N_OF_BONDS) : no integers'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not.associated(cv%i)

   if (size(cv%i).lt.2) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (N_OF_BONDS) : too few integers'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! size(cv%i).lt.3

   if (mod(size(cv%i), 2).ne.0) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (N_OF_BONDS) : number of integers is odd'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! sep.eq.size(cv%i)

   do a = 1, size(cv%i) - 1, 2
      if (cv%i(a).lt.1.or.cv%i(a).gt.pmemd_natoms()) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
              (a)//',a,'//pfmt(cv%i(a))//',a,'//pfmt(pmemd_natoms())//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (N_OF_BONDS) : integer #', a, ' (', cv%i(a), &
               ') is out of range [1, ', pmemd_natoms(), ']'
         NFE_MASTER_ONLY_END
         call terminate()
      end if
      if (cv%i(a + 1).lt.1.or.cv%i(a + 1).gt.pmemd_natoms()) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
               (a + 1)//',a,'//pfmt(cv%i(a + 1))//',a,'//pfmt &
               (pmemd_natoms())//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (N_OF_BONDS) : integer #', a + 1, ' (', cv%i(a + 1), &
               ') is out of range [1, ', pmemd_natoms(), ']'
         NFE_MASTER_ONLY_END
         call terminate()
      end if
      if (cv%i(a).eq.cv%i(a + 1)) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
               (a)//',a,'//pfmt(a + 1)//',a,'//pfmt(cv%i(a))//',a/)') &
               NFE_WARNING, 'CV #', cvno, &
               ' (N_OF_BONDS) : integers #', a, ' and ', a + 1, &
               ' are equal (', cv%i(a), ')'
         NFE_MASTER_ONLY_END
      end if
   end do

   if (.not.associated(cv%r)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (N_OF_BONDS) : no reals found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cvi)

   if (size(cv%r).ne.1) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (N_OF_BONDS) : number of reals is not 1'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cvi)

   if (cv%r(1).le.ZERO) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (N_OF_BONDS) : r(1).le.0.0D0 is .true.'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cvi)

end subroutine b_N_OF_BONDS

!=============================================================================

subroutine p_N_OF_BONDS(cv, lun)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   integer :: a
   character(4) :: aname

   nfe_assert(cv%type == COLVAR_N_OF_BONDS)

   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   nfe_assert(mod(size(cv%i), 2).eq.0)

   write (unit = lun, fmt = '(a,a,'//pfmt(cv%r(1), 3)//')') &
      NFE_INFO, '    d0 = ', cv%r(1)
   write (unit = lun, fmt = '(a,a)', advance = 'NO') NFE_INFO, ' pairs = ('

   do a = 1, size(cv%i)

      nfe_assert(cv%i(a).gt.0.and.cv%i(a).le.pmemd_natoms())
      aname = pmemd_atom_name(cv%i(a))

      write (unit = lun, fmt = '('//pfmt(cv%i(a))//',a,a,a)', advance = 'NO') &
         cv%i(a), ' [', trim(aname), ']'

      if (a.eq.size(cv%i)) then
         write (unit = lun, fmt = '(a)') ')'
      else if (mod(a, 4).eq.0) then
         write (unit = lun, fmt = '(a,/a,10x)', advance = 'NO') ',', NFE_INFO
      else if (mod(a, 2).eq.1) then
         write (unit = lun, fmt = '(a)', advance = 'NO') ' <=> '
      else
         write (unit = lun, fmt = '(a)', advance = 'NO') ', '
      end if

   end do

end subroutine p_N_OF_BONDS

!=============================================================================
!=========================== N_OF_STRUCTURES =================================
!=============================================================================
!
! input:
!
! cv%i = (a1, a2, a3, 0, b1, b2, b3, b4, 0, ..., c1, c2, c3, 0)
!
!     (a[1-3] - 1st group, b[1-4] - 2nd group, etc; an atom may
!      enter a few groups simultaneously; last zero is optional;
!      empty groups [e.g., 2+ zeros in a row] are not allowed)
!
! cv%r = (a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z, Ra, b1x, ...)
!
!        (reference coordinates followed by a threshod distance
!                         without '0' sentinel(s))
!
! value = (1 - rmsd1^6)/(1 - rmsd1^12) + ... + (1 - rmsdN^6)/(1 - rmsdN^12)
!
function v_N_OF_STRUCTURES(cv, x) result(value)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   double precision :: value

   type(colvar_t), intent(inout) :: cv

   double precision, intent(in) :: x(*)

   type(priv_t_NS), pointer :: priv
   integer :: g

#ifdef MPI
   integer :: error
   double precision :: accu
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_N_OF_STRUCTURES)

   priv => get_priv_NS(cv)
   nfe_assert(associated(priv%groups))
   nfe_assert(size(priv%groups).gt.0)

#ifdef MPI
   accu = ZERO
#else
   value = ZERO
#endif /* MPI */

   do g = 1, size(priv%groups)
#ifdef MPI
      if (mod(g + priv%first_cpu - 1, numtasks).eq.mytaskid) then
#endif /* MPI */
      call group_evaluate(cv, priv%groups(g), x)
#ifdef MPI
      accu = accu + priv%groups(g)%value
      endif
#else
      value = value + priv%groups(g)%value
#endif /* MPI */
   end do

#ifdef MPI
   call mpi_reduce(accu, value, 1, MPI_DOUBLE_PRECISION, &
                   MPI_SUM, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

contains
subroutine group_evaluate(cv, grp, x)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in)    :: cv
   type(group_t_NS),  intent(inout) :: grp

   double precision, intent(in) :: x(*)

   double precision :: cm(3), cm_nrm, lambda, rmsd2
   integer :: i, n, a, a3

   ! compute the center of mass (of the moving atoms)
   cm = ZERO
   n = 1
   do a = grp%i0, grp%i1
      a3 = 3*cv%i(a)
      cm = cm + grp%mass(n)*x(a3 - 2:a3)
      n = n + 1
   end do

   cm = cm/grp%total_mass

   ! populate cm_crd && compute the "norm"
   cm_nrm = ZERO
   n = 1
   do a = grp%i0, grp%i1
      a3 = 3*(n - 1)
      do i = 1, 3
         grp%cm_crd(a3 + i) = x(3*(cv%i(a) - 1) + i) - cm(i)
         cm_nrm = cm_nrm + grp%mass(n)*grp%cm_crd(a3 + i)**2
      end do
      n = n + 1
   end do

   call rmsd_q(grp%i1 - grp%i0 + 1, grp%mass, &
         grp%cm_crd, grp%ref_crd, lambda, grp%quaternion)

   grp%srmsd2 = ((grp%ref_nrm + cm_nrm) - 2*lambda)/grp%total_mass
   grp%srmsd2 = grp%srmsd2/grp%threshold**2

   grp%value = ONE/(ONE + grp%srmsd2**3)

end subroutine group_evaluate

end function v_N_OF_STRUCTURES

!=============================================================================

!
! assumes that atom positions have not been
!   changed since last call to 'value()'
!

subroutine f_N_OF_STRUCTURES(cv, x, fcv, f)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: x(*), fcv
   double precision, intent(inout) :: f(*)

   integer :: n, g, a, a3, f3
   double precision :: U(3,3), tmp
   type(priv_t_NS), pointer :: priv

   nfe_assert(cv%type == COLVAR_N_OF_STRUCTURES)
   nfe_assert(associated(cv%i))

   priv => get_priv_NS(cv)

   nfe_assert(associated(priv))
   nfe_assert(associated(priv%groups))
   nfe_assert(size(priv%groups).gt.0)

   do g = 1, size(priv%groups)

      call rmsd_q2u(priv%groups(g)%quaternion, U)

      tmp = (-6)*(priv%groups(g)%srmsd2*priv%groups(g)%value)**2
      tmp = tmp/priv%groups(g)%threshold**2
      tmp = fcv*tmp/priv%groups(g)%total_mass

      n = 1
      do a = priv%groups(g)%i0, priv%groups(g)%i1
         a3 = 3*n
         f3 = 3*cv%i(a)

         f(f3 - 2:f3) = f(f3 - 2:f3) &
            + tmp*priv%groups(g)%mass(n)*(priv%groups(g)%cm_crd(a3 - 2:a3) &
               - matmul(U, priv%groups(g)%ref_crd(a3 - 2:a3)))
         n = n + 1
      end do
   end do

end subroutine f_N_OF_STRUCTURES

!=============================================================================
subroutine b_N_OF_STRUCTURES(cv, cvno, amass)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)

   integer :: i, i0, n_atoms, n_groups, error

   type(priv_t_NS), pointer :: priv

   nfe_assert(cv%type == COLVAR_N_OF_STRUCTURES)

   ! very basic checks
   if (.not.associated(cv%i)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (N_OF_STRUCTURES) : no integers found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cv%i)

   if (.not.associated(cv%r)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_INFO, 'CV #', cvno, &
            ' (N_OF_STRUCTURES) : no reals (reference coordinates) found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if

   ! count the groups (number of zeros in the cv%i array)
   n_atoms = 0
   n_groups = 0
   i0 = 1

   do i = 1, size(cv%i)
      if (cv%i(i).eq.0) then
         n_groups = n_groups + 1
         if (i.eq.i0) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(i)//',a/)') &
                  NFE_ERROR, 'CV #', cvno, &
                  ' (N_OF_STRUCTURES) : unexpected zero (integer #', i, ')'
            NFE_MASTER_ONLY_END
            call terminate()
         end if ! i .eq. i0
         i0 = i + 1
      else if (cv%i(i).lt.0) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (N_OF_STRUCTURES) : negative integer'
         NFE_MASTER_ONLY_END
         call terminate()
      else
         n_atoms = n_atoms + 1
      end if
   end do

   if (size(cv%i).gt.0.and.cv%i(size(cv%i)).gt.0) &
      n_groups = n_groups + 1

   nfe_assert(n_groups.gt.0)

   if (size(cv%r).ne.(3*n_atoms + n_groups)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (N_OF_STRUCTURES) : wrong number of reals'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! size(cv%r) /= 3*n_atoms + n_groups

   ! allocate priv_t instance for this variable
   priv => new_priv_NS(cv)

   ! allocate/setup groups
   allocate(priv%groups(n_groups), stat = error)
   if (error /= 0) &
      NFE_OUT_OF_MEMORY

   i0 = 1 ! first atom
   n_atoms = 1
   n_groups = 0

   do i = 1, size(cv%i)
      if (cv%i(i).eq.0) then
         n_groups = n_groups + 1
         nfe_assert(n_groups.le.size(priv%groups))
         call group_bootstrap(priv%groups(n_groups), &
                              cv, cvno, amass, i0, i - 1, n_atoms)
         n_atoms = n_atoms + 3*(i - i0) + 1
         i0 = i + 1
      end if
   end do

   if (size(cv%i).gt.0.and.cv%i(size(cv%i)).gt.0) then
      n_groups = n_groups + 1
      nfe_assert(n_groups.le.size(priv%groups))
      call group_bootstrap(priv%groups(n_groups), &
                           cv, cvno, amass, i0, size(cv%i), n_atoms)
   end if

#ifdef MPI
   priv%first_cpu = 0
#endif /* MPI */

contains
subroutine group_bootstrap(grp, cv, cvno, amass, i0, i1, r0)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(group_t_NS),  intent(inout) :: grp

   type(colvar_t), intent(in) :: cv
   integer,        intent(in) :: cvno
   double precision,      intent(in) :: amass(*)

   integer, intent(in) :: i0, i1, r0

   integer :: a, b, n_atoms, error
   double precision :: cm(3)

   nfe_assert(cv%type == COLVAR_N_OF_STRUCTURES)
   nfe_assert(associated(cv%i).and.associated(cv%r))
   nfe_assert(i0.gt.0.and.i0.le.size(cv%i))
   nfe_assert(i1.gt.0.and.i1.le.size(cv%i))
   nfe_assert(r0.gt.0.and.r0.lt.size(cv%r))

   ! basic checks
   if (i0 + 2 > i1) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, &
            fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
            (i0)//',a,'//pfmt(i1)//',a/)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (N_OF_STRUCTURES) : too few integers in group (', i0, ':', i1, ')'
      NFE_MASTER_ONLY_END
      call terminate()
   end if

   do a = i0, i1
      nfe_assert(a .le. size(cv%i))
      if (cv%i(a) < 1 .or. cv%i(a) > pmemd_natoms()) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
                  (a)//',a,'//pfmt(cv%i(a))//',a,'//pfmt(pmemd_natoms())//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (N_OF_STRUCTURES) : integer #', a, ' (', cv%i(a), &
               ') is out of range [1, ', pmemd_natoms(), ']'
         NFE_MASTER_ONLY_END
         call terminate()
      end if
   end do

   ! check for duplicates
   do a = i0, i1
      do b = a + 1, i1
         if (cv%i(a) == cv%i(b)) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
                  (a)//',a,'//pfmt(b)//',a,'//pfmt(cv%i(a))//',a/)') &
                  NFE_ERROR, 'CV #', cvno, &
                  ' (N_OF_STRUCTURES) : integers #', a, ' and #', b, &
                  ' are equal (', cv%i(a), ')'
            NFE_MASTER_ONLY_END
            call terminate()
         end if ! cv%i(a) == cv%i(b)
      end do
   end do

   ! allocate/setup

   n_atoms = i1 - i0 + 1

   grp%i0 = i0
   grp%i1 = i1
   grp%r0 = r0

   allocate(grp%mass(n_atoms), grp%cm_crd(3*n_atoms), &
            grp%ref_crd(3*n_atoms), stat = error)

   if (error /= 0) &
      NFE_OUT_OF_MEMORY

   cm = ZERO
   grp%total_mass = ZERO

   do a = 1, n_atoms
      grp%mass(a) = amass(cv%i(a + grp%i0 - 1))
      grp%total_mass = grp%total_mass + grp%mass(a)
      cm = cm + grp%mass(a)*cv%r(r0 + 3*(a - 1):r0 + 3*a - 1)
   end do

   nfe_assert(grp%total_mass.gt.ZERO)
   cm = cm/grp%total_mass

   ! translate reference coordinates to CM frame
   grp%ref_nrm = ZERO
   do a = 1, n_atoms
      do b = 1, 3
         grp%ref_crd(3*(a - 1) + b) = cv%r(r0 + 3*(a - 1) + b - 1) - cm(b)
         grp%ref_nrm = grp%ref_nrm + grp%mass(a)*grp%ref_crd(3*(a - 1) + b)**2
      end do
   end do

   ! save the threshold distance
   grp%threshold = cv%r(r0 + 3*n_atoms)
   if (grp%threshold.le.ZERO) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
            (r0 + 3*n_atoms)//',a,'//pfmt(grp%threshold, 2)//',a/)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (N_OF_STRUCTURES) : threshold distance (r[', &
            r0 + 3*n_atoms, '] = ', grp%threshold, ') must be positive'
      NFE_MASTER_ONLY_END
      call terminate()
   end if

end subroutine group_bootstrap

end subroutine b_N_OF_STRUCTURES

!=============================================================================
subroutine c_N_OF_STRUCTURES(cv)

   NFE_USE_AFAILED

   implicit none

   type(colvar_t), intent(inout) :: cv

   integer :: g
   type(priv_t_NS), pointer :: priv

   nfe_assert(cv%type == COLVAR_N_OF_STRUCTURES)
   nfe_assert(associated(cv%i).and.size(cv%i).gt.0)

   priv => get_priv_NS(cv)
   nfe_assert(associated(priv%groups))

   do g = 1, size(priv%groups)
      call group_finalize(priv%groups(g))
   end do

   deallocate(priv%groups)
   call del_priv_NS(cv)
   
contains
subroutine group_finalize(grp)

   NFE_USE_AFAILED

   implicit none

   type(group_t_NS), intent(inout) :: grp

   nfe_assert(associated(grp%mass))
   nfe_assert(associated(grp%cm_crd))
   nfe_assert(associated(grp%ref_crd))

   deallocate(grp%mass, grp%cm_crd, grp%ref_crd)

end subroutine group_finalize

end subroutine c_N_OF_STRUCTURES

!=============================================================================
subroutine p_N_OF_STRUCTURES(cv, lun)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   integer :: g
   type(priv_t_NS), pointer :: priv

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_N_OF_STRUCTURES)
   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   priv => get_priv_NS(cv)
   nfe_assert(associated(priv%groups))

   do g = 1, size(priv%groups)
      write (unit = lun, fmt = '(a,a,'//pfmt(g)//',a)') &
         NFE_INFO, '<> group <> #', g, ':'
      call group_print(cv, priv%groups(g), lun)
   end do

contains
subroutine group_print(cv, grp, lun)

   use nfe_lib_mod
   
   implicit none

   type(colvar_t), intent(in) :: cv
   type(group_t_NS),  intent(in) :: grp
   integer,        intent(in) :: lun

   integer :: a, c
   character(4) :: aname

   nfe_assert(associated(cv%i).and.associated(cv%r))
   nfe_assert(grp%i0.gt.0.and.grp%i0.le.size(cv%i))
   nfe_assert(grp%i1.gt.0.and.grp%i1.le.size(cv%i))

   write (unit = lun, fmt = '(a,a,'//pfmt(grp%threshold, 3)//',a)') &
      NFE_INFO, 'threshold = ', grp%threshold, ' A'
   write (unit = lun, fmt = '(a,a)', advance = 'NO') NFE_INFO, 'atoms = ('

   c = 1
   do a = grp%i0, grp%i1

      nfe_assert(cv%i(a) > 0 .and. cv%i(a) <= pmemd_natoms())
      aname = pmemd_atom_name(cv%i(a))

      write (unit = lun, fmt = '('//pfmt(cv%i(a))//',a,a,a)', advance = 'NO') &
         cv%i(a), ' [', trim(aname), ']'

      if (a == grp%i1) then
         write (unit = lun, fmt = '(a)') ')'
      else if (mod(c, 5) == 0) then
         write (unit = lun, fmt = '(a,/a,3x)', advance = 'NO') ',', NFE_INFO
      else
         write (unit = lun, fmt = '(a)', advance = 'NO') ', '
      end if

      c = c + 1
   end do

   write (unit = lun, fmt = '(a,a)') &
         NFE_INFO, 'reference coordinates :'

   c = 0
   do a = grp%i0, grp%i1
      nfe_assert(grp%r0 + c + 2 <= size(cv%r))
      write (unit = lun, fmt = '(a,3x,i5,a,f8.3,a,f8.3,a,f8.3)') &
         NFE_INFO, cv%i(a), ' : ', &
         cv%r(grp%r0 + c), ', ', cv%r(grp%r0 + c + 1), &
         ', ', cv%r(grp%r0 + c + 2)
      c = c + 3
   end do

end subroutine group_print

end subroutine p_N_OF_STRUCTURES
!=============================================================================
function get_priv_NS(cv) result(ptr)

   NFE_USE_AFAILED

   implicit none

   type(colvar_t), intent(in) :: cv
   type(priv_t_NS), pointer :: ptr

   ptr => priv_list_NS
   do while (associated(ptr))
      if (ptr%tag.eq.cv%tag) &
         exit
      ptr => ptr%next
   end do

   nfe_assert(associated(ptr))

end function get_priv_NS
!=============================================================================
! allocates and appends to the list
function new_priv_NS(cv) result(ptr)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(inout) :: cv

   type(priv_t_NS), pointer :: ptr
   type(priv_t_NS), pointer :: head

   integer :: error

   allocate(ptr, stat = error)
   if (error.ne.0) &
      NFE_OUT_OF_MEMORY

   ptr%next => null()

   if (.not.associated(priv_list_NS)) then
      ptr%tag = 0
      priv_list_NS => ptr
   else
      head => priv_list_NS
      do while (associated(head%next))
         head => head%next
      end do
      ptr%tag = head%tag + 1
      head%next => ptr
   end if

   cv%tag = ptr%tag

end function new_priv_NS
!=============================================================================
! removes from the list and deallocates
subroutine del_priv_NS(cv)

   NFE_USE_AFAILED

   implicit none

   type(colvar_t), intent(in) :: cv

   type(priv_t_NS), pointer :: curr, prev

   nfe_assert(associated(priv_list_NS))

   curr => priv_list_NS
   if (curr%tag.eq.cv%tag) then
      priv_list_NS => curr%next
   else
      prev => curr
      curr => curr%next
      do while (associated(curr))
         if (curr%tag.eq.cv%tag) then
            prev%next => curr%next
            exit
         end if
         prev => curr
         curr => curr%next
      end do
   end if

   nfe_assert(associated(curr))
   nfe_assert(curr%tag.eq.cv%tag)

   deallocate(curr)

end subroutine del_priv_NS
!=============================================================================
!================================ L_COD ======================================
!=============================================================================
!
! LCOD (Linear Combination Of Distances)
!
! cv%i = (a11, a12, a21, a22, ..., aN1, aN2)
!
!     indexes of the participating atoms
!
! cv%r = (r1, r2, ..., rN)
!
!     non-zero weights
!
! value = r1*d1 + r2*d2 + ... + rN*dN
!
!     d1 -- distance between atoms #a11 and #a12
!     d2 -- distance between atoms #a21 and #a22
!                       . . .
!     dN -- distance between atoms #aN1 and #aN2
!

function v_LCOD(cv, x) result(value)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   double precision :: value

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: x(*)

   integer :: a1, a2, n

   nfe_assert(cv%type.eq.COLVAR_LCOD)

   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   nfe_assert(size(cv%i).gt.0)
   nfe_assert(size(cv%r).gt.0)

   nfe_assert(mod(size(cv%i), 2).eq.0)
   nfe_assert((size(cv%i)/2).eq.size(cv%r))

   value = ZERO

   NFE_MASTER_ONLY_BEGIN

      do n = 1, size(cv%r)

         nfe_assert(cv%i(2*n - 1).gt.0.and.cv%i(2*n - 1).le.pmemd_natoms())
         a1 = 3*cv%i(2*n - 1) - 2

         nfe_assert(cv%i(2*n).gt.0.and.cv%i(2*n).le.pmemd_natoms())
         a2 = 3*cv%i(2*n) - 2

         nfe_assert(a1.ne.a2)

         value = value + cv%r(n)*distance(x(a1:a1 + 2), x(a2:a2 + 2))

      end do

   NFE_MASTER_ONLY_END

end function v_LCOD

!=============================================================================

subroutine f_LCOD(cv, x, fcv, f)

#  ifndef NFE_DISABLE_ASSERT
   use nfe_lib_mod
#  endif /* NFE_DISABLE_ASSERT */

   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: x(*), fcv

   double precision, intent(inout) :: f(*)

   double precision :: d1(3), d2(3)

   integer :: a1, a2, n
#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type.eq.COLVAR_LCOD)

   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   nfe_assert(size(cv%i).gt.0)
   nfe_assert(size(cv%r).gt.0)

   nfe_assert(mod(size(cv%i), 2).eq.0)
   nfe_assert((size(cv%i)/2).eq.size(cv%r))

   NFE_MASTER_ONLY_BEGIN

      do n = 1, size(cv%r)

         nfe_assert(cv%i(2*n - 1).gt.0.and.cv%i(2*n - 1).le.pmemd_natoms())
         a1 = 3*cv%i(2*n - 1) - 2

         nfe_assert(cv%i(2*n).gt.0.and.cv%i(2*n).le.pmemd_natoms())
         a2 = 3*cv%i(2*n) - 2

         nfe_assert(a1.ne.a2)

         call distance_d(x(a1:a1 + 2), x(a2:a2 + 2), d1, d2)

         f(a1:a1 + 2) = f(a1:a1 + 2) + cv%r(n)*fcv*d1
         f(a2:a2 + 2) = f(a2:a2 + 2) + cv%r(n)*fcv*d2

      end do

   NFE_MASTER_ONLY_END
#ifdef MPI
   call mpi_bcast(f(1:3*pmemd_natoms()), 3*pmemd_natoms(), MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

end subroutine f_LCOD

!=============================================================================

subroutine b_LCOD(cv, cvno, amass)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)
   
   integer :: n, j

   nfe_assert(cv%type.eq.COLVAR_LCOD)

   if (.not.associated(cv%i)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (LCOD) : no integers found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not.associated(cv%i)

   if (.not.associated(cv%r)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (LCOD) : no reals found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not.associated(cv%r)

   if (mod(size(cv%i), 2).ne.0) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (LCOD) : odd number of integers'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! mod(size(cv%i), 2).ne.0

   if (size(cv%i).lt.2) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (LCOD) : too few integers'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! size(cv%i).lt.2

   if (size(cv%r).ne.size(cv%i)/2) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (LCOD) : size(cv%r).ne.size(cv%i)/2'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! size(cv%r).ne.size(cv%i)/2

   do n = 1, size(cv%r)
      if (abs(cv%r(n)).lt.1.0D-8) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(n)//',a/)') &
               NFE_ERROR, 'CV #', cvno, ' (LCOD) : real number #', n, &
               ' is too small'
         NFE_MASTER_ONLY_END
         call terminate()
      end if ! abs(cv%r(n)).lt.1.0D-8

      do j = 0, 1
         if (cv%i(2*n - j).lt.1.or.cv%i(2*n - j).gt.pmemd_natoms()) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(2*n - j)//',a,'//pfmt &
                  (cv%i(2*n - j))//',a,'//pfmt(pmemd_natoms())//',a/)') &
                  NFE_ERROR, 'CV #', cvno, ' (LCOD) : integer #', &
                  (2*n - j), ' (', cv%i(2*n - j), ') is out of range [1, ', &
                  pmemd_natoms(), ']'
            NFE_MASTER_ONLY_END
            call terminate()
         end if
      end do

      if (cv%i(2*n - 1).eq.cv%i(2*n)) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(2*n - 1)//',a,'//pfmt &
               (2*n)//',a,'//pfmt(cv%i(2*n))//',a/)') NFE_ERROR, 'CV #', &
               cvno, ' (LCOD) : integers #', (2*n - 1), ' and #', (2*n), &
               ' are equal (', cv%i(2*n), ')'
         NFE_MASTER_ONLY_END
         call terminate()
      end if
   end do

end subroutine b_LCOD

!=============================================================================

subroutine p_LCOD(cv, lun)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   integer :: n, a1, a2
   character(4) :: aname1, aname2

   nfe_assert(is_master())
   nfe_assert(cv%type.eq.COLVAR_LCOD)

   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   nfe_assert(size(cv%i).gt.0)
   nfe_assert(size(cv%r).gt.0)

   nfe_assert(mod(size(cv%i), 2).eq.0)
   nfe_assert((size(cv%i)/2).eq.size(cv%r))

   do n = 1, size(cv%r)

      a1 = cv%i(2*n - 1)
      a2 = cv%i(2*n)

      aname1 = pmemd_atom_name(a1)
      aname2 = pmemd_atom_name(a2)

      write (unit = lun, fmt = '(a,4x,f8.3,a,'//pfmt &
         (a1)//',a,a,a,'//pfmt(a2)//',a,a,a)') NFE_INFO, cv%r(n), ' * (', &
         a1, ' [', trim(aname1), '] <=> ', a2, ' [', trim(aname2), '])'

   end do

end subroutine p_LCOD
!=============================================================================
!=========================== COS_OF_DIHEDRAL =================================
!=============================================================================
function v_COS_OF_DIHEDRAL(cv, x) result(value)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   double precision :: value

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: x(*)

   integer :: n, a1, a2, a3, a4

   nfe_assert(cv%type == COLVAR_COS_OF_DIHEDRAL)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).gt.0)
   nfe_assert(mod(size(cv%i), 4).eq.0)

   value = ZERO

   do n = 1, size(cv%i)/4

      nfe_assert(cv%i(4*n - 3).gt.0.and.cv%i(4*n - 3).le.pmemd_natoms())
      a1 = 3*cv%i(4*n - 3) - 2

      nfe_assert(cv%i(4*n - 2).gt.0.and.cv%i(4*n - 2).le.pmemd_natoms())
      a2 = 3*cv%i(4*n - 2) - 2

      nfe_assert(cv%i(4*n - 1).gt.0.and.cv%i(4*n - 1).le.pmemd_natoms())
      a3 = 3*cv%i(4*n - 1) - 2

      nfe_assert(cv%i(4*n - 0).gt.0.and.cv%i(4*n - 0).le.pmemd_natoms())
      a4 = 3*cv%i(4*n - 0) - 2

      nfe_assert(a1.ne.a2.and.a1.ne.a3.and.a1.ne.a4)
      nfe_assert(a2.ne.a3.and.a2.ne.a4)
      nfe_assert(a3.ne.a4)

      NFE_MASTER_ONLY_BEGIN
      value = value &
         + cos4(x(a1:a1 + 2), x(a2:a2 + 2), x(a3:a3 + 2), x(a4:a4 + 2))
      NFE_MASTER_ONLY_END

   end do

end function v_COS_OF_DIHEDRAL

!=============================================================================

subroutine f_COS_OF_DIHEDRAL(cv, x, fcv, f)

#ifndef NFE_DISABLE_ASSERT
   use nfe_lib_mod
#endif /* NFE_DISABLE_ASSERT */

   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: x(*), fcv

   double precision, intent(inout) :: f(*)

   double precision :: d1(3), d2(3), d3(3), d4(3)

   integer :: n, a1, a2, a3, a4
#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_COS_OF_DIHEDRAL)
   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).gt.0)
   nfe_assert(mod(size(cv%i), 4).eq.0)

   NFE_MASTER_ONLY_BEGIN 
   do n = 1, size(cv%i)/4

      nfe_assert(cv%i(4*n - 3).gt.0.and.cv%i(4*n - 3).le.pmemd_natoms())
      a1 = 3*cv%i(4*n - 3) - 2

      nfe_assert(cv%i(4*n - 2).gt.0.and.cv%i(4*n - 2).le.pmemd_natoms())
      a2 = 3*cv%i(4*n - 2) - 2

      nfe_assert(cv%i(4*n - 1).gt.0.and.cv%i(4*n - 1).le.pmemd_natoms())
      a3 = 3*cv%i(4*n - 1) - 2

      nfe_assert(cv%i(4*n - 0).gt.0.and.cv%i(4*n - 0).le.pmemd_natoms())
      a4 = 3*cv%i(4*n - 0) - 2

      nfe_assert(a1.ne.a2.and.a1.ne.a3.and.a1.ne.a4)
      nfe_assert(a2.ne.a3.and.a2.ne.a4)
      nfe_assert(a3.ne.a4)

      call cos4_d(x(a1:a1 + 2), x(a2:a2 + 2), &
         x(a3:a3 + 2), x(a4:a4 + 2), d1, d2, d3, d4)

      f(a1:a1 + 2) = f(a1:a1 + 2) + fcv*d1
      f(a2:a2 + 2) = f(a2:a2 + 2) + fcv*d2
      f(a3:a3 + 2) = f(a3:a3 + 2) + fcv*d3
      f(a4:a4 + 2) = f(a4:a4 + 2) + fcv*d4

   end do
   NFE_MASTER_ONLY_END
#ifdef MPI
   call mpi_bcast(f(1:3*pmemd_natoms()), 3*pmemd_natoms(), MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

end subroutine f_COS_OF_DIHEDRAL

!=============================================================================

subroutine b_COS_OF_DIHEDRAL(cv, cvno, amass)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)

   integer :: n, j, l

   nfe_assert(cv%type.eq.COLVAR_COS_OF_DIHEDRAL)

   if (.not.associated(cv%i)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (COS_OF_DIHEDRAL) : no integers found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not.associated(cv%i)

   if (mod(size(cv%i), 4).ne.0) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (COS_OF_DIHEDRAL) : number of integers is not a multiple of 4'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! mod(size(cv%i), 4).ne.0

   if (size(cv%i).lt.4) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (COS_OF_DIHEDRAL) : too few integers'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! size(cv%i).lt.4

   do n = 1, size(cv%i)/4
      do j = 0, 3
         if (cv%i(4*n - j).lt.1.or.cv%i(4*n - j).gt.pmemd_natoms()) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(4*n - j)//',a,'//pfmt &
                  (cv%i(4*n - j))//',a,'//pfmt(pmemd_natoms())//',a/)') &
                  NFE_ERROR, 'CV #', cvno, ' (COS_OF_DIHEDRAL) : integer #', &
                  (4*n - j), ' (', cv%i(4*n - j), ') is out of range [1, ', &
                  pmemd_natoms(), ']'
            NFE_MASTER_ONLY_END
            call terminate()
         end if

         do l = 0, 3
            if (l.ne.j.and.cv%i(4*n - l).eq.cv%i(4*n - j)) then
               NFE_MASTER_ONLY_BEGIN
                  write (unit = ERR_UNIT, &
                     fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
                     (4*n - l)//',a,'//pfmt(4*n - j)//',a,'//pfmt &
                     (cv%i(4*n - j))//',a/)') NFE_ERROR, 'CV #', &
                     cvno, ' (COS_OF_DIHEDRAL) : integers #', (4*n - l), &
                     ' and #', (4*n - j), ' are equal (', cv%i(4*n - j), ')'
               NFE_MASTER_ONLY_END
               call terminate()
            end if
         end do
      end do
   end do

end subroutine b_COS_OF_DIHEDRAL

!=============================================================================

subroutine p_COS_OF_DIHEDRAL(cv, lun)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   integer :: n, a1, a2, a3, a4
   character(4) :: aname1, aname2, aname3, aname4

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_COS_OF_DIHEDRAL)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).gt.0)
   nfe_assert(mod(size(cv%i), 4).eq.0)

   do n = 1, size(cv%i)/4

      a1 = cv%i(4*n - 3)
      a2 = cv%i(4*n - 2)
      a3 = cv%i(4*n - 1)
      a4 = cv%i(4*n - 0)

      aname1 = pmemd_atom_name(a1)
      aname2 = pmemd_atom_name(a2)
      aname3 = pmemd_atom_name(a3)
      aname4 = pmemd_atom_name(a4)

      write (unit = lun, fmt = '(a,8x,'//pfmt(a1)//',a,a,a,'//pfmt &
            (a2)//',a,a,a,'//pfmt(a3)//',a,a,a,'//pfmt(a4)//',a,a,a)') &
         NFE_INFO, a1, ' [', trim(aname1), '] ==> ', &
                    a2, ' [', trim(aname2), '] ==> ', &
                    a3, ' [', trim(aname3), '] ==> ', &
                    a4, ' [', trim(aname4), ']'
   end do

end subroutine p_COS_OF_DIHEDRAL
!=============================================================================

!
! for A == B == C == D, a torsion along B == C is arccos([ABxBC]*[CDx(-BC)])
!        and its sign is given by sign(BC*[[ABxBC]x[CDx(-BC)]])
!

pure double precision function cos4(r1, r2, r3, r4)

   implicit none

   double precision, intent(in) :: r1(3), r2(3), r3(3), r4(3)

   cos4 = cos3(r2 - r1, r3 - r2, r4 - r3)

end function cos4

!=============================================================================

subroutine cos4_d(r1, r2, r3, r4, d1, d2, d3, d4)

   implicit none

   double precision, intent(in) :: r1(3), r2(3), r3(3), r4(3)
   double precision, intent(out) :: d1(3), d2(3), d3(3), d4(3)

   double precision :: t(3)

   call cos3_d(r2 - r1, r3 - r2, r4 - r3, d1, t, d4)

   d2 = d1 - t
   d1 = - d1
   d3 = t - d4

end subroutine cos4_d

!=============================================================================

pure double precision function cos3(v1, v2, v3)

   implicit none

   double precision, intent(in) :: v1(3), v2(3), v3(3)

   double precision :: n1(3), n2(3)

   n1 = cross3(v1, v2)
   n2 = cross3(v2, v3)

   cos3 = cosine(n1, n2)

end function cos3

!=============================================================================

subroutine cos3_d(v1, v2, v3, d1, d2, d3)

   implicit none

   double precision, intent(in) :: v1(3), v2(3), v3(3)
   double precision, intent(out) :: d1(3), d2(3), d3(3)

   double precision :: n1(3), n2(3), dc1(3), dc2(3), c

   n1 = cross3(v1, v2)
   n2 = cross3(v2, v3)

   c = cosine_d(n1, n2, dc1, dc2)

   d1 = cross3(v2, dc1)
   d2 = cross3(dc1, v1) + cross3(v3, dc2)
   d3 = cross3(dc2, v2)

end subroutine cos3_d

!=============================================================================
!============================== COM_ANGLE ====================================
!=============================================================================
! vbabin-at-ncsu-dot-edu, 09/08/2010
!
! COM_ANGLE -- angle formed by the centers of mass of three atom groups
!
! input: i = (a1, ..., aN, 0, b1, ..., bM, 0, c1, ..., cK, 0)
! last zero is optional; a?/b?/c? -- indices of the participating atoms

function v_COM_ANGLE(cv, x) result(value)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   double precision :: value

   type(colvar_t), intent(in) :: cv
   double precision, intent(in) :: x(*)

   integer :: n

   double precision :: cm1(3), cm2(3), cm3(3)

   nfe_assert(cv%type == COLVAR_COM_ANGLE)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).ge.5)
   nfe_assert(associated(cv%r))
   nfe_assert(size(cv%r).eq.size(cv%i))

#ifdef MPI
   if (mytaskid.eq.0) then
#endif /* MPI */
   n = 1
   call group_com(cv, x, n, cm1)
   nfe_assert(n.le.size(cv%i))
   call group_com(cv, x, n, cm2)
   nfe_assert(n.le.size(cv%i))
   call group_com(cv, x, n, cm3)

   value = angle(cm1, cm2, cm3)
#ifdef MPI
   else
      value = ZERO
   end if
#endif /* MPI */

end function v_COM_ANGLE

!=============================================================================

subroutine f_COM_ANGLE(cv, x, fcv, f)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   double precision, intent(in) :: x(*), fcv
   double precision, intent(inout) :: f(*)

   double precision :: d1(3), d2(3), d3(3)
   double precision :: cm1(3), cm2(3), cm3(3)

   integer :: n
#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_COM_ANGLE)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).ge.5)
   nfe_assert(associated(cv%r))
   nfe_assert(size(cv%r).eq.size(cv%i))

   NFE_MASTER_ONLY_BEGIN

   n = 1
   call group_com(cv, x, n, cm1)
   call group_com(cv, x, n, cm2)
   call group_com(cv, x, n, cm3)

   if (distance(cm1, cm2).lt.1.0d-8.or.distance(cm2, cm3).lt.1.0d-8) then
      write (unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR, &
         'COM_ANGLE : centers of mass got too close to each other'
      call terminate()
   end if ! distance too small

   call angle_d(cm1, cm2, cm3, d1, d2, d3)

   d1 = fcv*d1
   d2 = fcv*d2
   d3 = fcv*d3

   n = 1
   call group_com_d(cv, f, d1, n)
   call group_com_d(cv, f, d2, n)
   call group_com_d(cv, f, d3, n)

   NFE_MASTER_ONLY_END
#ifdef MPI
   call mpi_bcast(f(1:3*pmemd_natoms()), 3*pmemd_natoms(), MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

end subroutine f_COM_ANGLE

!=============================================================================

subroutine b_COM_ANGLE(cv, cvno, amass)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)

   integer :: error, ngroups

   nfe_assert(cv%type == COLVAR_COM_ANGLE)

   call com_check_i(cv%i, cvno, 'COM_ANGLE', ngroups)
   if (ngroups.ne.3) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (COM_ANGLE) : number of atom groups is not three'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! ngroups.ne.3

   if (associated(cv%r)) &
      deallocate(cv%r)

   allocate(cv%r(size(cv%i)), stat = error)
   if (error.ne.0) &
      NFE_OUT_OF_MEMORY

   cv%r = ZERO
   do error = 1, size(cv%i)
      if (cv%i(error).gt.0) &
         cv%r(error) = amass(cv%i(error))
   end do

   call com_init_weights(cv%r)

end subroutine b_COM_ANGLE

!=============================================================================

subroutine p_COM_ANGLE(cv, lun)

#ifndef NFE_DISABLE_ASSERT
   use nfe_lib_mod
#endif /* NFE_DISABLE_ASSERT */

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_COM_ANGLE)
   nfe_assert(associated(cv%i))

   call com_print_i(cv%i, lun)

end subroutine p_COM_ANGLE

!=============================================================================
!============================= COM_TORSION ===================================
!=============================================================================
! vbabin-at-ncsu-dot-edu, 09/08/2010
!
! COM_TORSION -- dihedral angle formed by the centers
!                of mass of four atom groups
!
! input: i = (a1, ..., aN, 0, b1, ..., bM, 0, c1, ..., cK, 0, d1, ..., dL, 0)
! last zero is optional; a?/b?/c?/d? -- indices of the participating atoms

function v_COM_TORSION(cv, x) result(value)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   double precision :: value

   type(colvar_t), intent(in) :: cv
   double precision, intent(in) :: x(*)

   integer :: n

   double precision :: cm1(3), cm2(3), cm3(3), cm4(3)

   nfe_assert(cv%type == COLVAR_COM_TORSION)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).ge.7)
   nfe_assert(associated(cv%r))
   nfe_assert(size(cv%r).eq.size(cv%i))

#ifdef MPI
   if (mytaskid.eq.0) then
#endif /* MPI */
   n = 1
   call group_com(cv, x, n, cm1)
   nfe_assert(n.le.size(cv%i))
   call group_com(cv, x, n, cm2)
   nfe_assert(n.le.size(cv%i))
   call group_com(cv, x, n, cm3)
   nfe_assert(n.le.size(cv%i))
   call group_com(cv, x, n, cm4)

   value = torsion(cm1, cm2, cm3, cm4)
#ifdef MPI
   else
      value = ZERO
   end if
#endif /* MPI */

end function v_COM_TORSION

!=============================================================================

subroutine f_COM_TORSION(cv, x, fcv, f)

   NFE_USE_AFAILED

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   double precision, intent(in) :: x(*), fcv
   double precision, intent(inout) :: f(*)

   double precision :: d1(3), d2(3), d3(3), d4(3)
   double precision :: cm1(3), cm2(3), cm3(3), cm4(3)

   double precision :: d12, d23, d34

   integer :: n
#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_COM_TORSION)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).ge.7)
   nfe_assert(associated(cv%r))
   nfe_assert(size(cv%r).eq.size(cv%i))

   NFE_MASTER_ONLY_BEGIN

   n = 1
   call group_com(cv, x, n, cm1)
   call group_com(cv, x, n, cm2)
   call group_com(cv, x, n, cm3)
   call group_com(cv, x, n, cm4)

   d12 = distance(cm1, cm2)
   d23 = distance(cm2, cm3)
   d34 = distance(cm3, cm4)

   if (d12.lt.TINY.or.d23.lt.TINY.or.d34.lt.TINY) then
      write (unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR, &
         'COM_TORSION : centers of mass got too close to each other'
      call terminate()
   end if ! too close

   call torsion_d(cm1, cm2, cm3, cm4, d1, d2, d3,  d4)

   d1 = fcv*d1
   d2 = fcv*d2
   d3 = fcv*d3
   d4 = fcv*d4

   n = 1
   call group_com_d(cv, f, d1, n)
   call group_com_d(cv, f, d2, n)
   call group_com_d(cv, f, d3, n)
   call group_com_d(cv, f, d4, n)

   NFE_MASTER_ONLY_END
   
#ifdef MPI
   call mpi_bcast(f(1:3*pmemd_natoms()), 3*pmemd_natoms(), MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

end subroutine f_COM_TORSION

!=============================================================================

subroutine b_COM_TORSION(cv, cvno, amass)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)

   integer :: error, ngroups

   nfe_assert(cv%type == COLVAR_COM_TORSION)

   call com_check_i(cv%i, cvno, 'COM_TORSION', ngroups)
   if (ngroups.ne.4) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (COM_TORSION) : number of atom groups is not four'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! ngroups.ne.4

   if (associated(cv%r)) &
      deallocate(cv%r)

   allocate(cv%r(size(cv%i)), stat = error)
   if (error.ne.0) &
      NFE_OUT_OF_MEMORY

   cv%r = ZERO
   do error = 1, size(cv%i)
      if (cv%i(error).gt.0) &
         cv%r(error) = amass(cv%i(error))
   end do

   call com_init_weights(cv%r)

end subroutine b_COM_TORSION

!=============================================================================

subroutine p_COM_TORSION(cv, lun)

#ifndef NFE_DISABLE_ASSERT
   use nfe_lib_mod
#endif /* NFE_DISABLE_ASSERT */

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_COM_TORSION)
   nfe_assert(associated(cv%i))

   call com_print_i(cv%i, lun)

end subroutine p_COM_TORSION

!=============================================================================
!============================ COM_DISTANCE ===================================
!=============================================================================
! vbabin-at-ncsu-dot-edu, 09/08/2010
!
! COM_DISTANCE -- distance between centers of mass of two atom groups
!
! input: i = (a1, ..., aN, 0, b1, ..., bM, 0)
! last zero is optional; a?/b? -- indices of the participating atoms

function v_COM_DISTANCE(cv, x) result(value)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   double precision :: value

   type(colvar_t), intent(in) :: cv
   double precision, intent(in) :: x(*)

   integer :: n

   double precision :: cm1(3), cm2(3)

   nfe_assert(cv%type == COLVAR_COM_DISTANCE)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).gt.3)
   nfe_assert(associated(cv%r))
   nfe_assert(size(cv%r).eq.size(cv%i))

#ifdef MPI
   if (mytaskid.eq.0) then
#endif /* MPI */
   n = 1
   call group_com(cv, x, n, cm1)
   nfe_assert(n.le.size(cv%i))
   call group_com(cv, x, n, cm2)

   value = distance(cm1, cm2)
#ifdef MPI
   else
      value = ZERO
   end if
#endif /* MPI */

end function v_COM_DISTANCE

!=============================================================================

subroutine f_COM_DISTANCE(cv, x, fcv, f)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   double precision, intent(in) :: x(*), fcv
   double precision, intent(inout) :: f(*)

   double precision :: d1(3), d2(3)
   double precision :: cm1(3), cm2(3)

   integer :: n
#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_COM_DISTANCE)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).gt.3)
   nfe_assert(associated(cv%r))
   nfe_assert(size(cv%r).eq.size(cv%i))

   NFE_MASTER_ONLY_BEGIN

   n = 1
   call group_com(cv, x, n, cm1)
   call group_com(cv, x, n, cm2)

   if (distance(cm1, cm2).lt.1.0d-8) then
      write (unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR, &
         'COM_DISTANCE : the distance got too small'
      call terminate()
   end if ! distance too small

   call distance_d(cm1, cm2, d1, d2)

   d1 = fcv*d1
   d2 = fcv*d2

   n = 1
   call group_com_d(cv, f, d1, n)
   call group_com_d(cv, f, d2, n)

   NFE_MASTER_ONLY_END
   
#ifdef MPI
   call mpi_bcast(f(1:3*pmemd_natoms()), 3*pmemd_natoms(), MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

end subroutine f_COM_DISTANCE

!=============================================================================

subroutine b_COM_DISTANCE(cv, cvno, amass)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)

   integer :: error, ngroups

   nfe_assert(cv%type == COLVAR_COM_DISTANCE)

   call com_check_i(cv%i, cvno, 'COM_DISTANCE', ngroups)
   if (ngroups.ne.2) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (COM_DISTANCE) : number of atom groups is not two'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! ngroups.ne.2

   if (associated(cv%r)) &
      deallocate(cv%r)

   allocate(cv%r(size(cv%i)), stat = error)
   if (error.ne.0) &
      NFE_OUT_OF_MEMORY

   cv%r = ZERO
   do error = 1, size(cv%i)
      if (cv%i(error).gt.0) &
         cv%r(error) = amass(cv%i(error))
   end do

   call com_init_weights(cv%r)

end subroutine b_COM_DISTANCE

!=============================================================================

subroutine p_COM_DISTANCE(cv, lun)

#ifndef NFE_DISABLE_ASSERT
   use nfe_lib_mod
#endif /* NFE_DISABLE_ASSERT */

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_COM_DISTANCE)
   nfe_assert(associated(cv%i))

   call com_print_i(cv%i, lun)

end subroutine p_COM_DISTANCE
!=============================================================================
!================================= PCA =======================================
!=============================================================================
! By Sishi Tang and Lin Fu, 
! June 17, 2011 
!
! input:
!
! cv%i = (a1, a2)
!
!     (a1, a2 : first and last atom for PCA projection) 
!
! cv%r = (a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z, b1x, ...)
!
!        (reference coordinates)  
! cv%avgcrd = (b1x, b1y, b1z, b2x, b2y, b2z,... ) 
! cv%evec = (v1x, v1y, ... ) 
! 
! value = ((fitcrd-avgcrd)*transpose(evec)) 
!
! Note: for v2, all solute atoms are fitted, and the projection is only calculated 
!       for selected atoms from a1 - a2. 
! dimension of arrays (clarification): 
!                   evec:   3*n_solut 
!                   avgcrd: 3*n_solut 
!                   refcrd: 3*n_solut 
!                   cm_crd: 3*n_solut 
!                  fit_crd: 3*natoms   

! corrected one:    evec:   3*size(cv%i)
!                   avgcrd: 3*size(cv%i)
!                   refcrd: 3*size(cv%j) read from both j and ref.crd
!                   cm_crd: 3*size(cv%j)
!                  fit_crd: 3*size(cv%i)

function v_PCA(cv, x) result(value)

   use nfe_lib_mod
   use parallel_dat_mod
   implicit none

   double precision :: value

   type(colvar_t), intent(inout) :: cv

   double precision, intent(in) :: x(*)
   double precision :: x_cm(3), lambda 
   double precision, pointer :: fit_crd(:) => null() 
   integer :: natoms, a, a3, i1, i2, i3, n, i, error, nsolut 

   type(priv_t_PCA), pointer :: priv

   nfe_assert(cv%type == COLVAR_PCA)
   nfe_assert(associated(cv%i))

!  i1: nsolut
!  i2: nref
!  i3: npca

   i1 = cv%i(1) 
   i2 = cv%i(2)
   i3 = cv%i(3)
  
!   nsolut = pmemd_nsolut() 
 
!   natoms = i2 - i1 + 1
!   nfe_assert(natoms > 1)

    natoms = i1
    nfe_assert(natoms > 1)
    
!    nfe_assert(i1 == nsolut)

!   allocate(fit_crd(3*natoms), stat = error)   
   allocate(fit_crd(3*cv%i(3)), stat = error)
   
   if (error /= 0) &
      NFE_OUT_OF_MEMORY

   priv => get_priv_PCA(cv)

   priv%value = ZERO
   x_cm = ZERO
   
   ! compute COM of moving atoms for all solutes 
   n = 1 
   do a = 1, cv%i(1)
   	  if(cv%state_ref(a) == 0) cycle
      a3 = 3*a 
!      x_cm = x_cm + priv%mass(n)*x(a3 - 2:a3)  
      x_cm = x_cm + priv%mass(a)*x(a3 - 2:a3)  

      n = n + 1 
   end do

   nfe_assert(i2 == n-1)

   x_cm = x_cm/priv%total_mass 
   
!   write(*,*) "x_cm=, worldrank=", x_cm, worldrank
   
   ! computer translated moving atoms for all solutes     
   n = 1
   do a = 1, cv%i(1)
      a3 = 3*(n - 1)
      do i = 1, 3
         priv%cm_crd(a3 + i) = x(3*(a - 1) + i) - x_cm(i)
      end do
      n = n + 1
   end do

   ! Calculate quaternion 
   call rmsd_q1(natoms, cv%state_ref, priv%mass, priv%cm_crd,& 
               priv%ref_crd, lambda, priv%quaternion)
  
   ! Calculate transposed rotation matrix U 
   call rmsd_q3u(priv%quaternion, priv%U) 
   
!   write(*,*) "priv%U=, worldrank=", priv%U, worldrank
 
   ! Calculate fitted crd (from i1-i2 ONLY) wrt to original ref_crd     
   ! best fit X for ref_crd : U*cm_crd + ref_cm 
   ! then subject <X> to get zero mean 
   ! make it mass weighted 
   ! avgcrd includes i1:i2  
   ! cm_crd includes 1:nsolut
   ! fit_crd includes i1:i2  
   n = 1 
!   do a = i1, i2
    do a = 1, cv%i(1)
    	if(cv%state_pca(a) == 0) cycle
      a3 = 3*a
!      fit_crd(3*n-2:3*n) = (matmul(priv%U, priv%cm_crd(a3-2:a3)) &
!                        + priv%ref_cm - cv%avgcrd(a3-2:a3))  & 
!                        * sqrt(priv%mass(a)) 

       fit_crd(3*n-2:3*n) = (matmul(priv%U, priv%cm_crd(a3-2:a3)) &
                         + priv%ref_cm - cv%avgcrd(3*n-2:3*n))  & 
                         * sqrt(priv%mass(a)) 
      n = n + 1 

   end do

   nfe_assert(i3 == n-1)


   ! Calculate projection
   ! P = XT  
!   priv%value = dot_product(fit_crd,cv%evec(i1*3-2:i2*3)) 
   priv%value = dot_product(fit_crd,cv%evec(1:3*i3))
   value = priv%value  
   

   deallocate(fit_crd)

end function v_PCA

!=============================================================================

!
! assumes that atom positions have not been
!   changed since last call to 'value()'
!

subroutine f_PCA(cv, x, fcv, f)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: x(*), fcv
   double precision, intent(inout) :: f(*)

#ifdef MPI_OLD
   integer :: a_first, a_last
#endif

   integer :: natoms, npca, a, a3, i1, i2, i3  

   type(priv_t_PCA), pointer :: priv

   nfe_assert(cv%type == COLVAR_PCA)
   nfe_assert(associated(cv%i))

   i1 = cv%i(1) 
   i2 = cv%i(2)
   i3 = cv%i(3)
!   natoms = i2 - i1 + 1
!   natoms = i1
   npca = i3
   
!   nfe_assert(natoms > 1)
   nfe_assert(npca > 1)

   priv => get_priv_PCA(cv)
   !nfe_assert(priv%value > ZERO)
   nfe_assert(priv%total_mass > ZERO)

#ifdef MPI_OLD
!   a = natoms/worldsize
   a = npca/worldsize
   if (a.gt.0) then
      if (mytaskid.ne.(worldsize - 1)) then
         a_first = 1 + mytaskid*a
         a_last = (mytaskid + 1)*a
      else
         a_first = 1 + mytaskid*a
         ! a_last = natoms
         a_last = npca
         
      end if
   else
      if (mytaskid.eq.0) then
         a_first = 1
         ! a_last = natoms
         a_last = npca
      else
         a_first = 1
         a_last = 0
      end if
   end if
   do a = a_first, a_last
#else
!   do a = 1, natoms
   do a = 1, npca
#endif /* MPI_OLD */
!      a3 = 3*(i1 + a - 1)
       a3 = 3*cv%ipca_to_i(a)
!      if(state_pca(a) == 0) cycle
      ! dc/dx = transpose(T)*U 
!      f(a3 - 2:a3) = f(a3 - 2:a3) + fcv * matmul(cv%evec(a3 - 2:a3), priv%U)*sqrt(priv%mass(i1 + a - 1))! now evec has the same index as force 
       f(a3 - 2:a3) = f(a3 - 2:a3) + fcv * matmul(cv%evec(3*a - 2 : 3*a), priv%U)*sqrt(priv%mass(cv%ipca_to_i(a)))
       
   end do

end subroutine f_PCA

!=============================================================================

subroutine b_PCA(cv, cvno, amass)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)

   integer   :: a, b, error, crdsize, nsolut, natoms  
   double precision :: total_mass

   type(priv_t_PCA), pointer :: priv



   nfe_assert(cv%type == COLVAR_PCA)

!   natoms = cv%i(2) - cv%i(1) + 1 
!   nsolut = pmemd_nsolut()
!   nfe_assert(nsolut == cv%i(1))

!  make sure that there are only three integers 
   call check_i(cv%i, cvno, 'PCA', 3)
   if (.not. cv%i(3) >= 2) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(a,a,'//pfmt(cvno)//',a)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (PCA) : too few integers'
      NFE_MASTER_ONLY_END
      call terminate()
   end if

   priv => new_priv_PCA(cv)

   allocate(priv%mass(cv%i(1)), priv%ref_crd(3*cv%i(1)), &
            priv%cm_crd(3*cv%i(1)), stat = error)      
   

   if (error.ne.0) &
      NFE_OUT_OF_MEMORY

   priv%ref_cm = ZERO
   priv%total_mass = ZERO
  
   total_mass = ZERO

   ! calculate total mass and refcm 
   do a = 1, cv%i(1)
   	  priv%mass(a) = amass(a)
   	  if(cv%state_ref(a) == 0 ) cycle
      total_mass = total_mass + priv%mass(a)
      priv%ref_cm = priv%ref_cm + priv%mass(a)*cv%r(3*a-2:3*a)
   end do

   priv%total_mass = total_mass 
   priv%ref_cm = priv%ref_cm/total_mass

   ! translate reference coordinates to CM frame
   ! ref_crd and r should be the same size 
   do a = 1, cv%i(1)
     do b = 1, 3
        priv%ref_crd(3*(a - 1) + b) = cv%r(3*(a - 1) + b ) - priv%ref_cm(b)
     end do
   end do
   

end subroutine b_PCA

!=============================================================================

subroutine p_PCA(cv, lun)

   NFE_USE_AFAILED

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   nfe_assert(cv%type == COLVAR_PCA)
   nfe_assert(associated(cv%i))

!  call print_i(cv%i, lun)
   call print_pca(cv%i, lun)

end subroutine p_PCA

!=============================================================================

subroutine c_PCA(cv)

   NFE_USE_AFAILED

   implicit none

   type(colvar_t), intent(inout) :: cv

   type(priv_t_PCA), pointer :: priv

   nfe_assert(cv%type.eq.COLVAR_PCA)

   priv => get_priv_PCA(cv)
   nfe_assert(associated(priv))

   deallocate(priv%mass, priv%cm_crd, priv%ref_crd) 

   call del_priv_PCA(cv)

end subroutine c_PCA
!=============================================================================
function get_priv_PCA(cv) result(ptr)

   NFE_USE_AFAILED

   implicit none

   type(colvar_t), intent(in) :: cv
   type(priv_t_PCA), pointer :: ptr

   ptr => priv_list_PCA
   do while (associated(ptr))
      if (ptr%tag.eq.cv%tag) &
         exit
      ptr => ptr%next
   end do

   nfe_assert(associated(ptr))

end function get_priv_PCA
!=============================================================================
! allocates and appends to the list
function new_priv_PCA(cv) result(ptr)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(inout) :: cv

   type(priv_t_PCA), pointer :: ptr
   type(priv_t_PCA), pointer :: head

   integer :: error

   allocate(ptr, stat = error)
   if (error.ne.0) &
      NFE_OUT_OF_MEMORY

   ptr%next => null()

   if (.not.associated(priv_list_PCA)) then
      ptr%tag = 0
      priv_list_PCA => ptr
   else
      head => priv_list_PCA
      do while (associated(head%next))
         head => head%next
      end do
      ptr%tag = head%tag + 1
      head%next => ptr
   end if

   cv%tag = ptr%tag

end function new_priv_PCA
!=============================================================================
! removes from the list and deallocates
subroutine del_priv_PCA(cv)

   NFE_USE_AFAILED

   implicit none

   type(colvar_t), intent(in) :: cv

   type(priv_t_PCA), pointer :: curr, prev

   nfe_assert(associated(priv_list_PCA))

   curr => priv_list_PCA
   if (curr%tag.eq.cv%tag) then
      priv_list_PCA => curr%next
   else
      prev => curr
      curr => curr%next
      do while (associated(curr))
         if (curr%tag.eq.cv%tag) then
            prev%next => curr%next
            exit
         end if
         prev => curr
         curr => curr%next
      end do
   end if

   nfe_assert(associated(curr))
   nfe_assert(curr%tag.eq.cv%tag)

   deallocate(curr)

end subroutine del_priv_PCA
!=============================================================================
!=========================== SIN_OF_DIHEDRAL =================================
!=============================================================================
!
! written by mmoradi-at-ncsu-dot-edu (12/2008)
!
! SIN_OF_DIHEDRAL (sum of sines of dihedral angles)
!
! cv%i = (a11, a12, a13, a14,
!         a21, a22, a23, a24,
!         ...,
!         aN1, aN2, aN3, aN4)
!
!     indexes of the participating atoms
!
! value = sin(d1) + sin(d2) + ... + sin(dN)
!
!   where d1 is dihedral formed by the atoms #a11, #a13, #a13, #a14;
!         d2 is formed by #a21, #a22, #a23, #a24, and so on
!

function v_SIN_OF_DIHEDRAL(cv, x) result(value)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   double precision :: value

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: x(*)

   integer :: n, a1, a2, a3, a4

   nfe_assert(cv%type == COLVAR_SIN_OF_DIHEDRAL)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).gt.0)
   nfe_assert(mod(size(cv%i), 4).eq.0)

   value = ZERO

   do n = 1, size(cv%i)/4

      nfe_assert(cv%i(4*n - 3).gt.0.and.cv%i(4*n - 3).le.pmemd_natoms())
      a1 = 3*cv%i(4*n - 3) - 2

      nfe_assert(cv%i(4*n - 2).gt.0.and.cv%i(4*n - 2).le.pmemd_natoms())
      a2 = 3*cv%i(4*n - 2) - 2

      nfe_assert(cv%i(4*n - 1).gt.0.and.cv%i(4*n - 1).le.pmemd_natoms())
      a3 = 3*cv%i(4*n - 1) - 2

      nfe_assert(cv%i(4*n - 0).gt.0.and.cv%i(4*n - 0).le.pmemd_natoms())
      a4 = 3*cv%i(4*n - 0) - 2

      nfe_assert(a1.ne.a2.and.a1.ne.a3.and.a1.ne.a4)
      nfe_assert(a2.ne.a3.and.a2.ne.a4)
      nfe_assert(a3.ne.a4)

      NFE_MASTER_ONLY_BEGIN
      value = value &
         + sin4(x(a1:a1 + 2), x(a2:a2 + 2), x(a3:a3 + 2), x(a4:a4 + 2))
      NFE_MASTER_ONLY_END

   end do

end function v_SIN_OF_DIHEDRAL

!=============================================================================

subroutine f_SIN_OF_DIHEDRAL(cv, x, fcv, f)

#ifndef NFE_DISABLE_ASSERT
   use nfe_lib_mod
#endif /* NFE_DISABLE_ASSERT */

   use parallel_dat_mod

   implicit none

   type(colvar_t),   intent(in) :: cv
   double precision, intent(in) :: x(*), fcv
   double precision, intent(inout) :: f(*)

   double precision :: d1(3), d2(3), d3(3), d4(3)

   integer :: n, a1, a2, a3, a4
#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_SIN_OF_DIHEDRAL)
   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).gt.0)
   nfe_assert(mod(size(cv%i), 4).eq.0)

   NFE_MASTER_ONLY_BEGIN
   do n = 1, size(cv%i)/4

      nfe_assert(cv%i(4*n - 3).gt.0.and.cv%i(4*n - 3).le.pmemd_natoms())
      a1 = 3*cv%i(4*n - 3) - 2

      nfe_assert(cv%i(4*n - 2).gt.0.and.cv%i(4*n - 2).le.pmemd_natoms())
      a2 = 3*cv%i(4*n - 2) - 2

      nfe_assert(cv%i(4*n - 1).gt.0.and.cv%i(4*n - 1).le.pmemd_natoms())
      a3 = 3*cv%i(4*n - 1) - 2

      nfe_assert(cv%i(4*n - 0).gt.0.and.cv%i(4*n - 0).le.pmemd_natoms())
      a4 = 3*cv%i(4*n - 0) - 2

      nfe_assert(a1.ne.a2.and.a1.ne.a3.and.a1.ne.a4)
      nfe_assert(a2.ne.a3.and.a2.ne.a4)
      nfe_assert(a3.ne.a4)

      call sin4_d(x(a1:a1 + 2), x(a2:a2 + 2), &
         x(a3:a3 + 2), x(a4:a4 + 2), d1, d2, d3, d4)

      f(a1:a1 + 2) = f(a1:a1 + 2) + fcv*d1
      f(a2:a2 + 2) = f(a2:a2 + 2) + fcv*d2
      f(a3:a3 + 2) = f(a3:a3 + 2) + fcv*d3
      f(a4:a4 + 2) = f(a4:a4 + 2) + fcv*d4

   end do
   NFE_MASTER_ONLY_END
#ifdef MPI
   call mpi_bcast(f(1:3*pmemd_natoms()), 3*pmemd_natoms(), MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

end subroutine f_SIN_OF_DIHEDRAL

!=============================================================================

subroutine b_SIN_OF_DIHEDRAL(cv, cvno)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno

   integer :: n, j, l

   nfe_assert(cv%type.eq.COLVAR_SIN_OF_DIHEDRAL)

   if (.not.associated(cv%i)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (SIN_OF_DIHEDRAL) : no integers found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not.associated(cv%i)

   if (mod(size(cv%i), 4).ne.0) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (SIN_OF_DIHEDRAL) : number of integers is not a multiple of 4'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! mod(size(cv%i), 4).ne.0

   if (size(cv%i).lt.4) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (SIN_OF_DIHEDRAL) : too few integers'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! size(cv%i).lt.4

   do n = 1, size(cv%i)/4
      do j = 0, 3
         if (cv%i(4*n - j).lt.1.or.cv%i(4*n - j).gt.pmemd_natoms()) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(4*n - j)//',a,'//pfmt &
                  (cv%i(4*n - j))//',a,'//pfmt(pmemd_natoms())//',a/)') &
                  NFE_ERROR, 'CV #', cvno, ' (SIN_OF_DIHEDRAL) : integer #', &
                  (4*n - j), ' (', cv%i(4*n - j), ') is out of range [1, ', &
                  pmemd_natoms(), ']'
            NFE_MASTER_ONLY_END
            call terminate()
         end if

         do l = 0, 3
            if (l.ne.j.and.cv%i(4*n - l).eq.cv%i(4*n - j)) then
               NFE_MASTER_ONLY_BEGIN
                  write (unit = ERR_UNIT, &
                     fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
                     (4*n - l)//',a,'//pfmt(4*n - j)//',a,'//pfmt &
                     (cv%i(4*n - j))//',a/)') NFE_ERROR, 'CV #', &
                     cvno, ' (SIN_OF_DIHEDRAL) : integers #', (4*n - l), &
                     ' and #', (4*n - j), ' are equal (', cv%i(4*n - j), ')'
               NFE_MASTER_ONLY_END
               call terminate()
            end if
         end do
      end do
   end do

end subroutine b_SIN_OF_DIHEDRAL

!=============================================================================

subroutine p_SIN_OF_DIHEDRAL(cv, lun)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   integer :: n, a1, a2, a3, a4
   character(4) :: aname1, aname2, aname3, aname4

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_SIN_OF_DIHEDRAL)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).gt.0)
   nfe_assert(mod(size(cv%i), 4).eq.0)

   do n = 1, size(cv%i)/4

      a1 = cv%i(4*n - 3)
      a2 = cv%i(4*n - 2)
      a3 = cv%i(4*n - 1)
      a4 = cv%i(4*n - 0)

      aname1 = pmemd_atom_name(a1)
      aname2 = pmemd_atom_name(a2)
      aname3 = pmemd_atom_name(a3)
      aname4 = pmemd_atom_name(a4)

      write (unit = lun, fmt = '(a,8x,'//pfmt(a1)//',a,a,a,'//pfmt &
            (a2)//',a,a,a,'//pfmt(a3)//',a,a,a,'//pfmt(a4)//',a,a,a)') &
         NFE_INFO, a1, ' [', trim(aname1), '] ==> ', &
                    a2, ' [', trim(aname2), '] ==> ', &
                    a3, ' [', trim(aname3), '] ==> ', &
                    a4, ' [', trim(aname4), ']'
   end do

end subroutine p_SIN_OF_DIHEDRAL
!=============================================================================

!
! for A == B == C == D, a torsion along B == C is arccos([ABxBC]*[CDx(-BC)])
!        and its sign is given by sign(BC*[[ABxBC]x[CDx(-BC)]])
!

pure double precision function sin4(r1, r2, r3, r4)

   implicit none

   double precision, intent(in) :: r1(3), r2(3), r3(3), r4(3)

   sin4 = sin(torsion(r1, r2, r3, r4))

end function sin4

!=============================================================================

subroutine sin4_d(r1, r2, r3, r4, d1, d2, d3, d4)

   implicit none

   double precision, intent(in) :: r1(3), r2(3), r3(3), r4(3)
   double precision, intent(out) :: d1(3), d2(3), d3(3), d4(3)

   double precision :: cos_tor

   call torsion_d(r1, r2, r3, r4, d1, d2, d3, d4)
   
   cos_tor = cos4(r1, r2, r3, r4)

   d1 = d1 * cos_tor
   d2 = d2 * cos_tor
   d3 = d3 * cos_tor
   d4 = d4 * cos_tor

end subroutine sin4_d

!=============================================================================
!=========================== PAIR_DIHEDRAL ===================================
!=============================================================================
!
! written by mmoradi-at-ncsu-dot-edu (12/2008)
!
! PAIR_DIHEDRAL (Lambda Collective Variable: sum of cosines of sum of each pair of neighboring dihedral angles)
!
! cv%i = (a11, a12, a13, a14,
!         a21, a22, a23, a24,
!         ...,
!         aN1, aN2, aN3, aN4)
!
!     indexes of the participating atoms
!
! value = cos(d1+d2) + cos(d2+d3) + ... + cos(nN-1+dN)
!
!   where d1 is dihedral formed by the atoms #a11, #a13, #a13, #a14;
!         d2 is formed by #a21, #a22, #a23, #a24, and so on
!

function v_PAIR_DIHEDRAL(cv, x) result(value)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   double precision :: value

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: x(*)

   integer :: n, a1, a2, a3, a4
   double precision :: tor1, tor2

   nfe_assert(cv%type == COLVAR_PAIR_DIHEDRAL)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).gt.0)
   nfe_assert(mod(size(cv%i), 4).eq.0)

   value = ZERO

   do n = 1, size(cv%i)/4

      nfe_assert(cv%i(4*n - 3).gt.0.and.cv%i(4*n - 3).le.pmemd_natoms())
      a1 = 3*cv%i(4*n - 3) - 2

      nfe_assert(cv%i(4*n - 2).gt.0.and.cv%i(4*n - 2).le.pmemd_natoms())
      a2 = 3*cv%i(4*n - 2) - 2

      nfe_assert(cv%i(4*n - 1).gt.0.and.cv%i(4*n - 1).le.pmemd_natoms())
      a3 = 3*cv%i(4*n - 1) - 2

      nfe_assert(cv%i(4*n - 0).gt.0.and.cv%i(4*n - 0).le.pmemd_natoms())
      a4 = 3*cv%i(4*n - 0) - 2

      nfe_assert(a1.ne.a2.and.a1.ne.a3.and.a1.ne.a4)
      nfe_assert(a2.ne.a3.and.a2.ne.a4)
      nfe_assert(a3.ne.a4)

      NFE_MASTER_ONLY_BEGIN
      tor2 = torsion(x(a1:a1 + 2), x(a2:a2 + 2), x(a3:a3 + 2), x(a4:a4 + 2))

      if (n>1) then
         value = value &
            + cos(tor1 + tor2)
      end if 
      tor1 = tor2
      NFE_MASTER_ONLY_END

   end do

end function v_PAIR_DIHEDRAL

!=============================================================================

subroutine f_PAIR_DIHEDRAL(cv, x, fcv, f)

#ifndef NFE_DISABLE_ASSERT
   use nfe_lib_mod
#endif /* NFE_DISABLE_ASSERT */

   use parallel_dat_mod

   implicit none

   type(colvar_t),   intent(in) :: cv
   double precision, intent(in) :: x(*), fcv
   double precision, intent(inout) :: f(*)

   double precision :: d1(3), d2(3), d3(3), d4(3)
   double precision :: d1_(3), d2_(3), d3_(3), d4_(3)
   double precision :: dc, tor1, tor2
   integer :: n, a1, a2, a3, a4, b1, b2, b3, b4
#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_PAIR_DIHEDRAL)
   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).gt.0)
   nfe_assert(mod(size(cv%i), 4).eq.0)

   NFE_MASTER_ONLY_BEGIN
   do n = 1, size(cv%i)/4

      nfe_assert(cv%i(4*n - 3).gt.0.and.cv%i(4*n - 3).le.pmemd_natoms())
      a1 = 3*cv%i(4*n - 3) - 2

      nfe_assert(cv%i(4*n - 2).gt.0.and.cv%i(4*n - 2).le.pmemd_natoms())
      a2 = 3*cv%i(4*n - 2) - 2

      nfe_assert(cv%i(4*n - 1).gt.0.and.cv%i(4*n - 1).le.pmemd_natoms())
      a3 = 3*cv%i(4*n - 1) - 2

      nfe_assert(cv%i(4*n - 0).gt.0.and.cv%i(4*n - 0).le.pmemd_natoms())
      a4 = 3*cv%i(4*n - 0) - 2

      nfe_assert(a1.ne.a2.and.a1.ne.a3.and.a1.ne.a4)
      nfe_assert(a2.ne.a3.and.a2.ne.a4)
      nfe_assert(a3.ne.a4)

      call torsion_d(x(a1:a1 + 2), x(a2:a2 + 2), &
         x(a3:a3 + 2), x(a4:a4 + 2), d1, d2, d3, d4)

      tor2 = torsion(x(a1:a1 + 2), x(a2:a2 + 2), x(a3:a3 + 2), x(a4:a4 + 2))

      if (n>1) then

         dc = - sin(tor1 + tor2)

         f(a1:a1 + 2) = f(a1:a1 + 2) + fcv*d1 * dc
         f(a2:a2 + 2) = f(a2:a2 + 2) + fcv*d2 * dc
         f(a3:a3 + 2) = f(a3:a3 + 2) + fcv*d3 * dc
         f(a4:a4 + 2) = f(a4:a4 + 2) + fcv*d4 * dc

         f(b1:b1 + 2) = f(b1:b1 + 2) + fcv*d1_ * dc
         f(b2:b2 + 2) = f(b2:b2 + 2) + fcv*d2_ * dc
         f(b3:b3 + 2) = f(b3:b3 + 2) + fcv*d3_ * dc
         f(b4:b4 + 2) = f(b4:b4 + 2) + fcv*d4_ * dc

      end if

      tor1 = tor2

      b1 = a1
      b2 = a2
      b3 = a3
      b4 = a4

      d1_ = d1
      d2_ = d2
      d3_ = d3
      d4_ = d4

   end do
   NFE_MASTER_ONLY_END
#ifdef MPI
   call mpi_bcast(f(1:3*pmemd_natoms()), 3*pmemd_natoms(), MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

end subroutine f_PAIR_DIHEDRAL

!=============================================================================

subroutine b_PAIR_DIHEDRAL(cv, cvno)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   integer :: n, j, l

   nfe_assert(cv%type.eq.COLVAR_PAIR_DIHEDRAL)

   if (.not.associated(cv%i)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (PAIR_DIHEDRAL) : no integers found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not.associated(cv%i)

   if (mod(size(cv%i), 4).ne.0) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (PAIR_DIHEDRAL) : number of integers is not a multiple of 4'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! mod(size(cv%i), 4).ne.0

   if (size(cv%i).lt.4) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (PAIR_DIHEDRAL) : too few integers'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! size(cv%i).lt.4

   do n = 1, size(cv%i)/4
      do j = 0, 3
         if (cv%i(4*n - j).lt.1.or.cv%i(4*n - j).gt.pmemd_natoms()) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(4*n - j)//',a,'//pfmt &
                  (cv%i(4*n - j))//',a,'//pfmt(pmemd_natoms())//',a/)') &
                  NFE_ERROR, 'CV #', cvno, ' (PAIR_DIHEDRAL) : integer #', &
                  (4*n - j), ' (', cv%i(4*n - j), ') is out of range [1, ', &
                  pmemd_natoms(), ']'
            NFE_MASTER_ONLY_END
            call terminate()
         end if

         do l = 0, 3
            if (l.ne.j.and.cv%i(4*n - l).eq.cv%i(4*n - j)) then
               NFE_MASTER_ONLY_BEGIN
                  write (unit = ERR_UNIT, &
                     fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
                     (4*n - l)//',a,'//pfmt(4*n - j)//',a,'//pfmt &
                     (cv%i(4*n - j))//',a/)') NFE_ERROR, 'CV #', &
                     cvno, ' (PAIR_DIHEDRAL) : integers #', (4*n - l), &
                     ' and #', (4*n - j), ' are equal (', cv%i(4*n - j), ')'
               NFE_MASTER_ONLY_END
               call terminate()
            end if
         end do
      end do
   end do

end subroutine b_PAIR_DIHEDRAL

!=============================================================================

subroutine p_PAIR_DIHEDRAL(cv, lun)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   integer :: n, a1, a2, a3, a4
   character(4) :: aname1, aname2, aname3, aname4

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_PAIR_DIHEDRAL)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).gt.0)
   nfe_assert(mod(size(cv%i), 4).eq.0)

   do n = 1, size(cv%i)/4

      a1 = cv%i(4*n - 3)
      a2 = cv%i(4*n - 2)
      a3 = cv%i(4*n - 1)
      a4 = cv%i(4*n - 0)

      aname1 = pmemd_atom_name(a1)
      aname2 = pmemd_atom_name(a2)
      aname3 = pmemd_atom_name(a3)
      aname4 = pmemd_atom_name(a4)

      write (unit = lun, fmt = '(a,8x,'//pfmt(a1)//',a,a,a,'//pfmt &
            (a2)//',a,a,a,'//pfmt(a3)//',a,a,a,'//pfmt(a4)//',a,a,a)') &
         NFE_INFO, a1, ' [', trim(aname1), '] ==> ', &
                    a2, ' [', trim(aname2), '] ==> ', &
                    a3, ' [', trim(aname3), '] ==> ', &
                    a4, ' [', trim(aname4), ']'
   end do

end subroutine p_PAIR_DIHEDRAL

!=============================================================================
!============================ PATTERN_DIHEDRAL ===============================
!=============================================================================
!
! written by mmoradi-at-ncsu-dot-edu (12/2008)
!
! PATTERN_DIHEDRAL (Gamma Collective Variable: Labeling each cis-trans pattern by an integer)
!
! cv%i = (a11, a12, a13, a14,
!         a21, a22, a23, a24,
!         ...,
!         aN1, aN2, aN3, aN4)
!
!     indexes of the participating atoms
!
! value = (cos(d1)+1)/2 * 2^0 + (cos(d2)+1)/2 * 2^1 + ... + (cos(dN)+1)/2 * 2^(N-1)
!
!   where d1 is dihedral formed by the atoms #a11, #a13, #a13, #a14;
!         d2 is formed by #a21, #a22, #a23, #a24, and so on
!

function v_PATTERN_DIHEDRAL(cv, x) result(value)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   double precision :: value

   type(colvar_t),   intent(in) :: cv
   double precision, intent(in) :: x(*)

   integer :: n, a1, a2, a3, a4

   nfe_assert(cv%type == COLVAR_PATTERN_DIHEDRAL)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).gt.0)
   nfe_assert(mod(size(cv%i), 4).eq.0)

   value = ZERO

   do n = 1, size(cv%i)/4

      nfe_assert(cv%i(4*n - 3).gt.0.and.cv%i(4*n - 3).le.pmemd_natoms())
      a1 = 3*cv%i(4*n - 3) - 2

      nfe_assert(cv%i(4*n - 2).gt.0.and.cv%i(4*n - 2).le.pmemd_natoms())
      a2 = 3*cv%i(4*n - 2) - 2

      nfe_assert(cv%i(4*n - 1).gt.0.and.cv%i(4*n - 1).le.pmemd_natoms())
      a3 = 3*cv%i(4*n - 1) - 2

      nfe_assert(cv%i(4*n - 0).gt.0.and.cv%i(4*n - 0).le.pmemd_natoms())
      a4 = 3*cv%i(4*n - 0) - 2

      nfe_assert(a1.ne.a2.and.a1.ne.a3.and.a1.ne.a4)
      nfe_assert(a2.ne.a3.and.a2.ne.a4)
      nfe_assert(a3.ne.a4)

      NFE_MASTER_ONLY_BEGIN
      value = value &
         + ((cos4(x(a1:a1 + 2), x(a2:a2 + 2), x(a3:a3 + 2), x(a4:a4 + 2)) + 1)/2) * (2**(n-1))
      NFE_MASTER_ONLY_END

   end do

end function v_PATTERN_DIHEDRAL

!=============================================================================

subroutine f_PATTERN_DIHEDRAL(cv, x, fcv, f)

#ifndef NFE_DISABLE_ASSERT
   use nfe_lib_mod
#endif /* NFE_DISABLE_ASSERT */

   use parallel_dat_mod

   implicit none

   type(colvar_t),   intent(in) :: cv
   double precision, intent(in) :: x(*), fcv
   double precision, intent(inout) :: f(*)

   double precision :: d1(3), d2(3), d3(3), d4(3)

   integer :: n, a1, a2, a3, a4
#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_PATTERN_DIHEDRAL)
   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).gt.0)
   nfe_assert(mod(size(cv%i), 4).eq.0)

   NFE_MASTER_ONLY_BEGIN
   do n = 1, size(cv%i)/4

      nfe_assert(cv%i(4*n - 3).gt.0.and.cv%i(4*n - 3).le.pmemd_natoms())
      a1 = 3*cv%i(4*n - 3) - 2

      nfe_assert(cv%i(4*n - 2).gt.0.and.cv%i(4*n - 2).le.pmemd_natoms())
      a2 = 3*cv%i(4*n - 2) - 2

      nfe_assert(cv%i(4*n - 1).gt.0.and.cv%i(4*n - 1).le.pmemd_natoms())
      a3 = 3*cv%i(4*n - 1) - 2

      nfe_assert(cv%i(4*n - 0).gt.0.and.cv%i(4*n - 0).le.pmemd_natoms())
      a4 = 3*cv%i(4*n - 0) - 2

      nfe_assert(a1.ne.a2.and.a1.ne.a3.and.a1.ne.a4)
      nfe_assert(a2.ne.a3.and.a2.ne.a4)
      nfe_assert(a3.ne.a4)

      call cos4_d(x(a1:a1 + 2), x(a2:a2 + 2), &
         x(a3:a3 + 2), x(a4:a4 + 2), d1, d2, d3, d4)

      f(a1:a1 + 2) = f(a1:a1 + 2) + (fcv*d1/2) * 2**(n-1)
      f(a2:a2 + 2) = f(a2:a2 + 2) + (fcv*d2/2) * 2**(n-1)
      f(a3:a3 + 2) = f(a3:a3 + 2) + (fcv*d3/2) * 2**(n-1)
      f(a4:a4 + 2) = f(a4:a4 + 2) + (fcv*d4/2) * 2**(n-1)

   end do
   NFE_MASTER_ONLY_END
#ifdef MPI
   call mpi_bcast(f(1:3*pmemd_natoms()), 3*pmemd_natoms(), MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

end subroutine f_PATTERN_DIHEDRAL

!=============================================================================

subroutine b_PATTERN_DIHEDRAL(cv, cvno)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   
   integer :: n, j, l

   nfe_assert(cv%type.eq.COLVAR_PATTERN_DIHEDRAL)

   if (.not.associated(cv%i)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (PATTERN_DIHEDRAL) : no integers found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not.associated(cv%i)

   if (mod(size(cv%i), 4).ne.0) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (PATTERN_DIHEDRAL) : number of integers is not a multiple of 4'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! mod(size(cv%i), 4).ne.0

   if (size(cv%i).lt.4) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_ERROR, 'CV #', cvno, ' (PATTERN_DIHEDRAL) : too few integers'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! size(cv%i).lt.4

   do n = 1, size(cv%i)/4
      do j = 0, 3
         if (cv%i(4*n - j).lt.1.or.cv%i(4*n - j).gt.pmemd_natoms()) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(4*n - j)//',a,'//pfmt &
                  (cv%i(4*n - j))//',a,'//pfmt(pmemd_natoms())//',a/)') &
                  NFE_ERROR, 'CV #', cvno, ' (PATTERN_DIHEDRAL) : integer #', &
                  (4*n - j), ' (', cv%i(4*n - j), ') is out of range [1, ', &
                  pmemd_natoms(), ']'
            NFE_MASTER_ONLY_END
            call terminate()
         end if

         do l = 0, 3
            if (l.ne.j.and.cv%i(4*n - l).eq.cv%i(4*n - j)) then
               NFE_MASTER_ONLY_BEGIN
                  write (unit = ERR_UNIT, &
                     fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
                     (4*n - l)//',a,'//pfmt(4*n - j)//',a,'//pfmt &
                     (cv%i(4*n - j))//',a/)') NFE_ERROR, 'CV #', &
                     cvno, ' (PATTERN_DIHEDRAL) : integers #', (4*n - l), &
                     ' and #', (4*n - j), ' are equal (', cv%i(4*n - j), ')'
               NFE_MASTER_ONLY_END
               call terminate()
            end if
         end do
      end do
   end do

end subroutine b_PATTERN_DIHEDRAL

!=============================================================================

subroutine p_PATTERN_DIHEDRAL(cv, lun)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   integer :: n, a1, a2, a3, a4
   character(4) :: aname1, aname2, aname3, aname4

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_PATTERN_DIHEDRAL)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).gt.0)
   nfe_assert(mod(size(cv%i), 4).eq.0)

   do n = 1, size(cv%i)/4

      a1 = cv%i(4*n - 3)
      a2 = cv%i(4*n - 2)
      a3 = cv%i(4*n - 1)
      a4 = cv%i(4*n - 0)

      aname1 = pmemd_atom_name(a1)
      aname2 = pmemd_atom_name(a2)
      aname3 = pmemd_atom_name(a3)
      aname4 = pmemd_atom_name(a4)

      write (unit = lun, fmt = '(a,8x,'//pfmt(a1)//',a,a,a,'//pfmt &
            (a2)//',a,a,a,'//pfmt(a3)//',a,a,a,'//pfmt(a4)//',a,a,a)') &
         NFE_INFO, a1, ' [', trim(aname1), '] ==> ', &
                    a2, ' [', trim(aname2), '] ==> ', &
                    a3, ' [', trim(aname3), '] ==> ', &
                    a4, ' [', trim(aname4), ']'
   end do

end subroutine p_PATTERN_DIHEDRAL
!=============================================================================
!============================== DNA_RISE =====================================
!=============================================================================
!
! written by fpan3-at-ncsu-dot-edu in 01/2016
!
! DNA_RISE (the step rise of DNA/RNA duplex, the 
!           distance between two planes which are nearly parallel)
!
! cv%i = (a11, ..., a1m,
!         a21, ..., a2n,
!         ...,
!         a61, ..., a6k)
!
! indexes of the participating atoms : 
! a1?, a2? : the phosphate atom or group of atoms on the backbone of first plane
! a3?      : the atoms of base pairs of first plane
! a4?, a5? : the phosphate atom or group of atoms on the backbone of second plane
! a6?      : the atoms of base pairs of second plane
!
! value = norm(ga3-ga6)*cos((ga3-ga6),h12)
!
!   where ga3 is the center of mass of a1?, ga6 is the center of mass of a6?
!         h12 is the average unit vector of the normal vectors of two planes
!

function v_DNA_RISE(cv, x) result(value)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   double precision :: value

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: x(*)

   integer :: n
   double precision :: cm1(3),cm2(3),cm3(3),cm4(3),cm5(3),cm6(3)

   nfe_assert(cv%type == COLVAR_DNA_RISE)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).ge.11)

   value = ZERO

      NFE_MASTER_ONLY_BEGIN

      n = 1
      call group_com(cv, x, n, cm1)
      call group_com(cv, x, n, cm2)
      call group_com(cv, x, n, cm3)
      call group_com(cv, x, n, cm4)
      call group_com(cv, x, n, cm5)
      call group_com(cv, x, n, cm6)
      
      value = rise_v(cm1,cm2,cm3,cm4,cm5,cm6)
      
      NFE_MASTER_ONLY_END

end function v_DNA_RISE

!=============================================================================

subroutine f_DNA_RISE(cv, x, fcv, f)

#ifndef NFE_DISABLE_ASSERT
   use nfe_lib_mod
#endif /* NFE_DISABLE_ASSERT */

   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: x(*), fcv

   double precision, intent(inout) :: f(*)

   double precision :: cm1(3),cm2(3),cm3(3),cm4(3),cm5(3),cm6(3)
   double precision :: d1(3), d2(3), d3(3), d4(3), d5(3), d6(3)

   integer :: n
#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_DNA_RISE)
   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).ge.11)

   NFE_MASTER_ONLY_BEGIN 
      n = 1
      call group_com(cv, x, n, cm1)
      call group_com(cv, x, n, cm2)
      call group_com(cv, x, n, cm3)
      call group_com(cv, x, n, cm4)
      call group_com(cv, x, n, cm5)
      call group_com(cv, x, n, cm6)
     
      call rise_d(cm1,cm2,cm3,cm4,cm5,cm6,d1,d2,d3,d4,d5,d6)

      d1 = fcv*d1
      d2 = fcv*d2
      d3 = fcv*d3
      d4 = fcv*d4
      d5 = fcv*d5
      d6 = fcv*d6
   
      n = 1
      call group_com_d(cv, f, d1, n)
      call group_com_d(cv, f, d2, n)
      call group_com_d(cv, f, d3, n)
      call group_com_d(cv, f, d4, n)
      call group_com_d(cv, f, d5, n)  
      call group_com_d(cv, f, d6, n)
      
   NFE_MASTER_ONLY_END
#ifdef MPI
   call mpi_bcast(f(1:3*pmemd_natoms()), 3*pmemd_natoms(), MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

end subroutine f_DNA_RISE

!=============================================================================

subroutine b_DNA_RISE(cv, cvno, amass)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)

   integer :: error, ngroups

   nfe_assert(cv%type == COLVAR_DNA_RISE)

   call com_check_i(cv%i, cvno, 'DNA_RISE', ngroups)
   if (ngroups.ne.6) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (COM_TORSION) : number of atom groups is not six'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! ngroups.ne.6

   if (associated(cv%r)) &
      deallocate(cv%r)

   allocate(cv%r(size(cv%i)), stat = error)
   if (error.ne.0) &
      NFE_OUT_OF_MEMORY

   cv%r = ZERO
   do error = 1, size(cv%i)
      if (cv%i(error).gt.0) &
         cv%r(error) = amass(cv%i(error))
   end do

   call com_init_weights(cv%r)

end subroutine b_DNA_RISE

!=============================================================================

subroutine p_DNA_RISE(cv, lun)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_DNA_RISE)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).gt.11)

   call com_print_i(cv%i,lun)

end subroutine p_DNA_RISE

!=============================================================================
! the function rise_v and rise_d are based on MATLAB code and post-processing
!

pure double precision function rise_v(r1, r2, r3, r4, r5, r6)

   implicit none

   double precision, intent(in) :: r1(3), r2(3), r3(3), r4(3), r5(3), r6(3)
   double precision  :: t2
   double precision  :: t3
   double precision  :: t4
   double precision  :: t5
   double precision  :: t6
   double precision  :: t24
   double precision  :: t7
   double precision  :: t8
   double precision  :: t10
   double precision  :: t11
   double precision  :: t26
   double precision  :: t27
   double precision  :: t28
   double precision  :: t9
   double precision  :: t30
   double precision  :: t31
   double precision  :: t32
   double precision  :: t12
   double precision  :: t13
   double precision  :: t14
   double precision  :: t15
   double precision  :: t16
   double precision  :: t17
   double precision  :: t37
   double precision  :: t18
   double precision  :: t19
   double precision  :: t21
   double precision  :: t22
   double precision  :: t39
   double precision  :: t40
   double precision  :: t41
   double precision  :: t20
   double precision  :: t43
   double precision  :: t44
   double precision  :: t45
   double precision  :: t23
   double precision  :: t25
   double precision  :: t29
   double precision  :: t33
   double precision  :: t34
   double precision  :: t35
   double precision  :: t36
   double precision  :: t38
   double precision  :: t42
   double precision  :: t46
   double precision  :: t47
   double precision  :: t48
   double precision  :: t49
   double precision  :: t50
   double precision  :: t51
   double precision  :: t54
   double precision  :: t55
   double precision  :: t56
   double precision  :: t52
   double precision  :: t59
   double precision  :: t60
   double precision  :: t61
   double precision  :: t53
   double precision  :: t57
   double precision  :: t58
   double precision  :: t62
   double precision  :: t63
   double precision  :: t64
   double precision  :: t0
 
      t2 = r1(1)-r2(1)
      t3 = r1(2)-r3(2)
      t4 = t2*t3
      t5 = r1(2)-r2(2)
      t6 = r1(1)-r3(1)
      t24 = t5*t6
      t7 = t4-t24
      t8 = abs(t7)
      t10 = r1(3)-r3(3)
      t11 = r1(3)-r2(3)
      t26 = t2*t10
      t27 = t6*t11
      t28 = t26-t27
      t9 = abs(t28)
      t30 = t5*t10
      t31 = t3*t11
      t32 = t30-t31
      t12 = abs(t32)
      t13 = r4(1)-r5(1)
      t14 = r4(2)-r6(2)
      t15 = t13*t14
      t16 = r4(2)-r5(2)
      t17 = r4(1)-r6(1)
      t37 = t16*t17
      t18 = t15-t37
      t19 = abs(t18)
      t21 = r4(3)-r6(3)
      t22 = r4(3)-r5(3)
      t39 = t13*t21
      t40 = t17*t22
      t41 = t39-t40
      t20 = abs(t41)
      t43 = t16*t21
      t44 = t14*t22
      t45 = t43-t44
      t23 = abs(t45)
      t25 = t8**2
      t29 = t9**2
      t33 = t12**2
      t34 = t25+t29+t33
      t35 = 1.0D0/sqrt(t34)
      t36 = t7*t35*(1.0D0/2.0D0)
      t38 = t19**2
      t42 = t20**2
      t46 = t23**2
      t47 = t38+t42+t46
      t48 = 1.0D0/sqrt(t47)
      t49 = t18*t48*(1.0D0/2.0D0)
      t50 = t36+t49
      t51 = abs(t50)
      t54 = t28*t35*(1.0D0/2.0D0)
      t55 = t41*t48*(1.0D0/2.0D0)
      t56 = t54+t55
      t52 = abs(t56)
      t59 = t32*t35*(1.0D0/2.0D0)
      t60 = t45*t48*(1.0D0/2.0D0)
      t61 = t59+t60
      t53 = abs(t61)
      t57 = t51**2
      t58 = t52**2
      t62 = t53**2
      t63 = t57+t58+t62
      t64 = 1.0D0/sqrt(t63)
      t0 = abs(t50*t64*(r3(3)-r6(3))-t56*t64*(r3(2)-r6(2))+t61*t64*(r3(1)-r6(1)))
      
      rise_v = t0

end function rise_v

!=============================================================================

subroutine rise_d(r1, r2, r3, r4, r5, r6, d1, d2, d3, d4, d5, d6)

   implicit none

   double precision, intent(in)  :: r1(3), r2(3), r3(3), r4(3), r5(3), r6(3)
   double precision, intent(out) :: d1(3), d2(3), d3(3), d4(3), d5(3), d6(3)
   double precision  :: t2
   double precision  :: t3
   double precision  :: t4
   double precision  :: t5
   double precision  :: t6
   double precision  :: t24
   double precision  :: t7
   double precision  :: t8
   double precision  :: t10
   double precision  :: t11
   double precision  :: t26
   double precision  :: t27
   double precision  :: t28
   double precision  :: t9
   double precision  :: t30
   double precision  :: t31
   double precision  :: t32
   double precision  :: t12
   double precision  :: t13
   double precision  :: t14
   double precision  :: t15
   double precision  :: t16
   double precision  :: t17
   double precision  :: t37
   double precision  :: t18
   double precision  :: t19
   double precision  :: t21
   double precision  :: t22
   double precision  :: t39
   double precision  :: t40
   double precision  :: t41
   double precision  :: t20
   double precision  :: t43
   double precision  :: t44
   double precision  :: t45
   double precision  :: t23
   double precision  :: t25
   double precision  :: t29
   double precision  :: t33
   double precision  :: t34
   double precision  :: t35
   double precision  :: t36
   double precision  :: t38
   double precision  :: t42
   double precision  :: t46
   double precision  :: t47
   double precision  :: t48
   double precision  :: t49
   double precision  :: t50
   double precision  :: t51
   double precision  :: t54
   double precision  :: t55
   double precision  :: t56
   double precision  :: t52
   double precision  :: t59
   double precision  :: t60
   double precision  :: t61
   double precision  :: t53
   double precision  :: t57
   double precision  :: t58
   double precision  :: t62
   double precision  :: t63
   double precision  :: t64
   double precision  :: t65
   double precision  :: t66
   double precision  :: t67
   double precision  :: t68
   double precision  :: t69
   double precision  :: t70
   double precision  :: t71
   double precision  :: t72
   double precision  :: t73
   double precision  :: t74
   double precision  :: t75
   double precision  :: t80
   double precision  :: t76
   double precision  :: t77
   double precision  :: t83
   double precision  :: t78
   double precision  :: t79
   double precision  :: t81
   double precision  :: t82
   double precision  :: t84
   double precision  :: t85
   double precision  :: t89
   double precision  :: t86
   double precision  :: t87
   double precision  :: t88
   double precision  :: t90
   double precision  :: t91
   double precision  :: t108
   double precision  :: t92
   double precision  :: t93
   double precision  :: t94
   double precision  :: t95
   double precision  :: t98
   double precision  :: t99
   double precision  :: t96
   double precision  :: t97
   double precision  :: t103
   double precision  :: t100
   double precision  :: t101
   double precision  :: t102
   double precision  :: t104
   double precision  :: t105
   double precision  :: t107
   double precision  :: t106
   double precision  :: t109
   double precision  :: t110
   double precision  :: t111
   double precision  :: t114
   double precision  :: t112
   double precision  :: t116
   double precision  :: t113
   double precision  :: t115
   double precision  :: t117
   double precision  :: t119
   double precision  :: t118
   double precision  :: t120
   double precision  :: t121
   double precision  :: t122
   double precision  :: t123
   double precision  :: t127
   double precision  :: t124
   double precision  :: t125
   double precision  :: t129
   double precision  :: t126
   double precision  :: t128
   double precision  :: t130
   double precision  :: t132
   double precision  :: t131
   double precision  :: t133
   double precision  :: t136
   double precision  :: t134
   double precision  :: t135
   double precision  :: t140
   double precision  :: t137
   double precision  :: t138
   double precision  :: t139
   double precision  :: t141
   double precision  :: t142
   double precision  :: t144
   double precision  :: t143
   double precision  :: t145
   double precision  :: t146
   double precision  :: t147
   double precision  :: t150
   double precision  :: t148
   double precision  :: t152
   double precision  :: t149
   double precision  :: t151
   double precision  :: t153
   double precision  :: t155
   double precision  :: t154
   double precision  :: t156
   double precision  :: t157
   double precision  :: t158
   double precision  :: t159
   double precision  :: t163
   double precision  :: t160
   double precision  :: t161
   double precision  :: t165
   double precision  :: t162
   double precision  :: t164
   double precision  :: t166
   double precision  :: t168
   double precision  :: t167
   double precision  :: t169
   double precision  :: t172
   double precision  :: t170
   double precision  :: t171
   double precision  :: t176
   double precision  :: t173
   double precision  :: t174
   double precision  :: t175
   double precision  :: t177
   double precision  :: t178
   double precision  :: t180
   double precision  :: t179
   double precision  :: t181
   double precision  :: t182
   double precision  :: t183
   double precision  :: t186
   double precision  :: t184
   double precision  :: t188
   double precision  :: t185
   double precision  :: t187
   double precision  :: t189
   double precision  :: t191
   double precision  :: t190
   double precision  :: t192
   double precision  :: t193
   double precision  :: t194
   double precision  :: t195
   double precision  :: t196
   double precision  :: t197
   double precision  :: t198
   double precision  :: t199
   double precision  :: t200
   double precision  :: t204
   double precision  :: t201
   double precision  :: t202
   double precision  :: t206
   double precision  :: t203
   double precision  :: t205
   double precision  :: t207
   double precision  :: t209
   double precision  :: t208
   double precision  :: t210
   double precision  :: t211
   double precision  :: t214
   double precision  :: t215
   double precision  :: t212
   double precision  :: t213
   double precision  :: t219
   double precision  :: t216
   double precision  :: t217
   double precision  :: t218
   double precision  :: t220
   double precision  :: t221
   double precision  :: t223
   double precision  :: t222
   double precision  :: t224
   double precision  :: t225
   double precision  :: t226
   double precision  :: t229
   double precision  :: t227
   double precision  :: t231
   double precision  :: t228
   double precision  :: t230
   double precision  :: t232
   double precision  :: t234
   double precision  :: t233
   double precision  :: t235
   double precision  :: t236
   double precision  :: t237
   double precision  :: t238
   double precision  :: t242
   double precision  :: t239
   double precision  :: t240
   double precision  :: t244
   double precision  :: t241
   double precision  :: t243
   double precision  :: t245
   double precision  :: t247
   double precision  :: t246
   double precision  :: t248
   double precision  :: t251
   double precision  :: t249
   double precision  :: t250
   double precision  :: t255
   double precision  :: t252
   double precision  :: t253
   double precision  :: t254
   double precision  :: t256
   double precision  :: t257
   double precision  :: t259
   double precision  :: t258
   double precision  :: t260
   double precision  :: t261
   double precision  :: t262
   double precision  :: t265
   double precision  :: t263
   double precision  :: t267
   double precision  :: t264
   double precision  :: t266
   double precision  :: t268
   double precision  :: t270
   double precision  :: t269
   double precision  :: t271
   double precision  :: t272
   double precision  :: t273
   double precision  :: t274
   double precision  :: t278
   double precision  :: t275
   double precision  :: t276
   double precision  :: t280
   double precision  :: t277
   double precision  :: t279
   double precision  :: t281
   double precision  :: t283
   double precision  :: t282
   double precision  :: t284
   double precision  :: t285
   double precision  :: t288
   double precision  :: t286
   double precision  :: t287
   double precision  :: t292
   double precision  :: t289
   double precision  :: t290
   double precision  :: t291
   double precision  :: t293
   double precision  :: t294
   double precision  :: t296
   double precision  :: t295
   double precision  :: t297
   double precision  :: t298
   double precision  :: t299
   double precision  :: t300
   double precision  :: t303
   double precision  :: t301
   double precision  :: t305
   double precision  :: t302
   double precision  :: t304
   double precision  :: t306
   double precision  :: t308
   double precision  :: t307
 
      t2 = r1(1)-r2(1)
      t3 = r1(2)-r3(2)
      t4 = t2*t3
      t5 = r1(2)-r2(2)
      t6 = r1(1)-r3(1)
      t24 = t5*t6
      t7 = t4-t24
      t8 = abs(t7)
      t10 = r1(3)-r3(3)
      t11 = r1(3)-r2(3)
      t26 = t2*t10
      t27 = t6*t11
      t28 = t26-t27
      t9 = abs(t28)
      t30 = t5*t10
      t31 = t3*t11
      t32 = t30-t31
      t12 = abs(t32)
      t13 = r4(1)-r5(1)
      t14 = r4(2)-r6(2)
      t15 = t13*t14
      t16 = r4(2)-r5(2)
      t17 = r4(1)-r6(1)
      t37 = t16*t17
      t18 = t15-t37
      t19 = abs(t18)
      t21 = r4(3)-r6(3)
      t22 = r4(3)-r5(3)
      t39 = t13*t21
      t40 = t17*t22
      t41 = t39-t40
      t20 = abs(t41)
      t43 = t16*t21
      t44 = t14*t22
      t45 = t43-t44
      t23 = abs(t45)
      t25 = t8**2
      t29 = t9**2
      t33 = t12**2
      t34 = t25+t29+t33
      t35 = 1.0D0/sqrt(t34)
      t36 = t7*t35*(1.0D0/2.0D0)
      t38 = t19**2
      t42 = t20**2
      t46 = t23**2
      t47 = t38+t42+t46
      t48 = 1.0D0/sqrt(t47)
      t49 = t18*t48*(1.0D0/2.0D0)
      t50 = t36+t49
      t51 = abs(t50)
      t54 = t28*t35*(1.0D0/2.0D0)
      t55 = t41*t48*(1.0D0/2.0D0)
      t56 = t54+t55
      t52 = abs(t56)
      t59 = t32*t35*(1.0D0/2.0D0)
      t60 = t45*t48*(1.0D0/2.0D0)
      t61 = t59+t60
      t53 = abs(t61)
      t57 = t51**2
      t58 = t52**2
      t62 = t53**2
      t63 = t57+t58+t62
      t64 = 1.0D0/sqrt(t63)
      t65 = r2(2)-r3(2)
      t66 = r3(3)-r6(3)
      t67 = r2(3)-r3(3)
      t68 = (t7/abs(t7))
      t69 = t8*t65*t68*2.0D0
      t70 = (t28/abs(t28))
      t71 = t9*t67*t70*2.0D0
      t72 = t69+t71
      t73 = 1.0D0/t34**(3.0D0/2.0D0)
      t74 = r3(2)-r6(2)
      t75 = t35*t65*(1.0D0/2.0D0)
      t80 = t7*t72*t73*(1.0D0/4.0D0)
      t76 = t75-t80
      t77 = t35*t67*(1.0D0/2.0D0)
      t83 = t28*t72*t73*(1.0D0/4.0D0)
      t78 = t77-t83
      t79 = (t50/abs(t50))
      t81 = t51*t76*t79*2.0D0
      t82 = (t56/abs(t56))
      t84 = t52*t78*t82*2.0D0
      t85 = (t61/abs(t61))
      t89 = t32*t53*t72*t73*t85*(1.0D0/2.0D0)
      t86 = t81+t84-t89
      t87 = 1.0D0/t63**(3.0D0/2.0D0)
      t88 = r3(1)-r6(1)
      t90 = t50*t64*t66
      t91 = t61*t64*t88
      t108 = t56*t64*t74
      t92 = t90+t91-t108
      t93 = (t92/abs(t92))
      t94 = r2(1)-r3(1)
      t95 = (t32/abs(t32))
      t98 = t8*t68*t94*2.0D0
      t99 = t12*t67*t95*2.0D0
      t96 = t98-t99
      t97 = t35*t94*(1.0D0/2.0D0)
      t103 = t7*t73*t96*(1.0D0/4.0D0)
      t100 = t97-t103
      t101 = t32*t73*t96*(1.0D0/4.0D0)
      t102 = t77+t101
      t104 = t53*t85*t102*2.0D0
      t105 = t28*t52*t73*t82*t96*(1.0D0/2.0D0)
      t107 = t51*t79*t100*2.0D0
      t106 = t104+t105-t107
      t109 = t9*t70*t94*2.0D0
      t110 = t12*t65*t95*2.0D0
      t111 = t109+t110
      t114 = t28*t73*t111*(1.0D0/4.0D0)
      t112 = t97-t114
      t116 = t32*t73*t111*(1.0D0/4.0D0)
      t113 = t75-t116
      t115 = t52*t82*t112*2.0D0
      t117 = t53*t85*t113*2.0D0
      t119 = t7*t51*t73*t79*t111*(1.0D0/2.0D0)
      t118 = t115+t117-t119
      t120 = t3*t8*t68*2.0D0
      t121 = t9*t10*t70*2.0D0
      t122 = t120+t121
      t123 = t3*t35*(1.0D0/2.0D0)
      t127 = t7*t73*t122*(1.0D0/4.0D0)
      t124 = t123-t127
      t125 = t10*t35*(1.0D0/2.0D0)
      t129 = t28*t73*t122*(1.0D0/4.0D0)
      t126 = t125-t129
      t128 = t51*t79*t124*2.0D0
      t130 = t52*t82*t126*2.0D0
      t132 = t32*t53*t73*t85*t122*(1.0D0/2.0D0)
      t131 = t128+t130-t132
      t133 = t6*t8*t68*2.0D0
      t136 = t10*t12*t95*2.0D0
      t134 = t133-t136
      t135 = t6*t35*(1.0D0/2.0D0)
      t140 = t7*t73*t134*(1.0D0/4.0D0)
      t137 = t135-t140
      t138 = t32*t73*t134*(1.0D0/4.0D0)
      t139 = t125+t138
      t141 = t53*t85*t139*2.0D0
      t142 = t28*t52*t73*t82*t134*(1.0D0/2.0D0)
      t144 = t51*t79*t137*2.0D0
      t143 = t141+t142-t144
      t145 = t6*t9*t70*2.0D0
      t146 = t3*t12*t95*2.0D0
      t147 = t145+t146
      t150 = t28*t73*t147*(1.0D0/4.0D0)
      t148 = t135-t150
      t152 = t32*t73*t147*(1.0D0/4.0D0)
      t149 = t123-t152
      t151 = t52*t82*t148*2.0D0
      t153 = t53*t85*t149*2.0D0
      t155 = t7*t51*t73*t79*t147*(1.0D0/2.0D0)
      t154 = t151+t153-t155
      t156 = t5*t8*t68*2.0D0
      t157 = t9*t11*t70*2.0D0
      t158 = t156+t157
      t159 = t5*t35*(1.0D0/2.0D0)
      t163 = t7*t73*t158*(1.0D0/4.0D0)
      t160 = t159-t163
      t161 = t11*t35*(1.0D0/2.0D0)
      t165 = t28*t73*t158*(1.0D0/4.0D0)
      t162 = t161-t165
      t164 = t51*t79*t160*2.0D0
      t166 = t52*t82*t162*2.0D0
      t168 = t32*t53*t73*t85*t158*(1.0D0/2.0D0)
      t167 = t164+t166-t168
      t169 = t2*t8*t68*2.0D0
      t172 = t11*t12*t95*2.0D0
      t170 = t169-t172
      t171 = t2*t35*(1.0D0/2.0D0)
      t176 = t7*t73*t170*(1.0D0/4.0D0)
      t173 = t171-t176
      t174 = t32*t73*t170*(1.0D0/4.0D0)
      t175 = t161+t174
      t177 = t53*t85*t175*2.0D0
      t178 = t28*t52*t73*t82*t170*(1.0D0/2.0D0)
      t180 = t51*t79*t173*2.0D0
      t179 = t177+t178-t180
      t181 = t2*t9*t70*2.0D0
      t182 = t5*t12*t95*2.0D0
      t183 = t181+t182
      t186 = t28*t73*t183*(1.0D0/4.0D0)
      t184 = t171-t186
      t188 = t32*t73*t183*(1.0D0/4.0D0)
      t185 = t159-t188
      t187 = t52*t82*t184*2.0D0
      t189 = t53*t85*t185*2.0D0
      t191 = t7*t51*t73*t79*t183*(1.0D0/2.0D0)
      t190 = t187+t189-t191
      t192 = r5(2)-r6(2)
      t193 = r5(3)-r6(3)
      t194 = (t18/abs(t18))
      t195 = t19*t192*t194*2.0D0
      t196 = (t41/abs(t41))
      t197 = t20*t193*t196*2.0D0
      t198 = t195+t197
      t199 = 1.0D0/t47**(3.0D0/2.0D0)
      t200 = t48*t192*(1.0D0/2.0D0)
      t204 = t18*t198*t199*(1.0D0/4.0D0)
      t201 = t200-t204
      t202 = t48*t193*(1.0D0/2.0D0)
      t206 = t41*t198*t199*(1.0D0/4.0D0)
      t203 = t202-t206
      t205 = t51*t79*t201*2.0D0
      t207 = t52*t82*t203*2.0D0
      t209 = t45*t53*t85*t198*t199*(1.0D0/2.0D0)
      t208 = t205+t207-t209
      t210 = r5(1)-r6(1)
      t211 = (t45/abs(t45))
      t214 = t19*t194*t210*2.0D0
      t215 = t23*t193*t211*2.0D0
      t212 = t214-t215
      t213 = t48*t210*(1.0D0/2.0D0)
      t219 = t18*t199*t212*(1.0D0/4.0D0)
      t216 = t213-t219
      t217 = t45*t199*t212*(1.0D0/4.0D0)
      t218 = t202+t217
      t220 = t53*t85*t218*2.0D0
      t221 = t41*t52*t82*t199*t212*(1.0D0/2.0D0)
      t223 = t51*t79*t216*2.0D0
      t222 = t220+t221-t223
      t224 = t20*t196*t210*2.0D0
      t225 = t23*t192*t211*2.0D0
      t226 = t224+t225
      t229 = t41*t199*t226*(1.0D0/4.0D0)
      t227 = t213-t229
      t231 = t45*t199*t226*(1.0D0/4.0D0)
      t228 = t200-t231
      t230 = t52*t82*t227*2.0D0
      t232 = t53*t85*t228*2.0D0
      t234 = t18*t51*t79*t199*t226*(1.0D0/2.0D0)
      t233 = t230+t232-t234
      t235 = t14*t19*t194*2.0D0
      t236 = t20*t21*t196*2.0D0
      t237 = t235+t236
      t238 = t14*t48*(1.0D0/2.0D0)
      t242 = t18*t199*t237*(1.0D0/4.0D0)
      t239 = t238-t242
      t240 = t21*t48*(1.0D0/2.0D0)
      t244 = t41*t199*t237*(1.0D0/4.0D0)
      t241 = t240-t244
      t243 = t51*t79*t239*2.0D0
      t245 = t52*t82*t241*2.0D0
      t247 = t45*t53*t85*t199*t237*(1.0D0/2.0D0)
      t246 = t243+t245-t247
      t248 = t17*t19*t194*2.0D0
      t251 = t21*t23*t211*2.0D0
      t249 = t248-t251
      t250 = t17*t48*(1.0D0/2.0D0)
      t255 = t18*t199*t249*(1.0D0/4.0D0)
      t252 = t250-t255
      t253 = t45*t199*t249*(1.0D0/4.0D0)
      t254 = t240+t253
      t256 = t53*t85*t254*2.0D0
      t257 = t41*t52*t82*t199*t249*(1.0D0/2.0D0)
      t259 = t51*t79*t252*2.0D0
      t258 = t256+t257-t259
      t260 = t17*t20*t196*2.0D0
      t261 = t14*t23*t211*2.0D0
      t262 = t260+t261
      t265 = t41*t199*t262*(1.0D0/4.0D0)
      t263 = t250-t265
      t267 = t45*t199*t262*(1.0D0/4.0D0)
      t264 = t238-t267
      t266 = t52*t82*t263*2.0D0
      t268 = t53*t85*t264*2.0D0
      t270 = t18*t51*t79*t199*t262*(1.0D0/2.0D0)
      t269 = t266+t268-t270
      t271 = t16*t19*t194*2.0D0
      t272 = t20*t22*t196*2.0D0
      t273 = t271+t272
      t274 = t16*t48*(1.0D0/2.0D0)
      t278 = t18*t199*t273*(1.0D0/4.0D0)
      t275 = t274-t278
      t276 = t22*t48*(1.0D0/2.0D0)
      t280 = t41*t199*t273*(1.0D0/4.0D0)
      t277 = t276-t280
      t279 = t51*t79*t275*2.0D0
      t281 = t52*t82*t277*2.0D0
      t283 = t45*t53*t85*t199*t273*(1.0D0/2.0D0)
      t282 = t279+t281-t283
      t284 = t56*t64
      t285 = t13*t19*t194*2.0D0
      t288 = t22*t23*t211*2.0D0
      t286 = t285-t288
      t287 = t13*t48*(1.0D0/2.0D0)
      t292 = t18*t199*t286*(1.0D0/4.0D0)
      t289 = t287-t292
      t290 = t45*t199*t286*(1.0D0/4.0D0)
      t291 = t276+t290
      t293 = t53*t85*t291*2.0D0
      t294 = t41*t52*t82*t199*t286*(1.0D0/2.0D0)
      t296 = t51*t79*t289*2.0D0
      t295 = t293+t294-t296
      t297 = t50*t64
      t298 = t13*t20*t196*2.0D0
      t299 = t16*t23*t211*2.0D0
      t300 = t298+t299
      t303 = t41*t199*t300*(1.0D0/4.0D0)
      t301 = t287-t303
      t305 = t45*t199*t300*(1.0D0/4.0D0)
      t302 = t274-t305
      t304 = t52*t82*t301*2.0D0
      t306 = t53*t85*t302*2.0D0
      t308 = t18*t51*t79*t199*t300*(1.0D0/2.0D0)
      t307 = t304+t306-t308
      d1(1) = -t93*(-t64*t66*t76+t64*t74*t78+t50*t66*t86*t87*(1.0D0/2.&
     &0D0)-t56*t74*t86*t87*(1.0D0/2.0D0)+t61*t86*t87*t88*(1.0D0/2.0D0)+t&
     &32*t64*t72*t73*t88*(1.0D0/4.0D0))
      d1(2) = -t93*(t64*t66*t100-t64*t88*t102+t50*t66*t87*t106*(1.0D0/&
     &2.0D0)-t56*t74*t87*t106*(1.0D0/2.0D0)+t61*t87*t88*t106*(1.0D0/2.0D&
     &0)+t28*t64*t73*t74*t96*(1.0D0/4.0D0))
      d1(3) = t93*(t64*t74*t112-t64*t88*t113+t50*t66*t87*t118*(1.0D0/2&
     &.0D0)-t56*t74*t87*t118*(1.0D0/2.0D0)+t61*t87*t88*t118*(1.0D0/2.0D0&
     &)+t7*t64*t66*t73*t111*(1.0D0/4.0D0))
      d2(1) = t93*(-t64*t66*t124+t64*t74*t126+t50*t66*t87*t131*(1.0D0/&
     &2.0D0)-t56*t74*t87*t131*(1.0D0/2.0D0)+t61*t87*t88*t131*(1.0D0/2.0D&
     &0)+t32*t64*t73*t88*t122*(1.0D0/4.0D0))
      d2(2) = t93*(t64*t66*t137-t64*t88*t139+t50*t66*t87*t143*(1.0D0/2&
     &.0D0)-t56*t74*t87*t143*(1.0D0/2.0D0)+t61*t87*t88*t143*(1.0D0/2.0D0&
     &)+t28*t64*t73*t74*t134*(1.0D0/4.0D0))
      d2(3) = -t93*(t64*t74*t148-t64*t88*t149+t50*t66*t87*t154*(1.0D0/&
     &2.0D0)-t56*t74*t87*t154*(1.0D0/2.0D0)+t61*t87*t88*t154*(1.0D0/2.0D&
     &0)+t7*t64*t66*t73*t147*(1.0D0/4.0D0))
      d3(1) = -t93*(-t61*t64-t64*t66*t160+t64*t74*t162+t50*t66*t87*t16&
     &7*(1.0D0/2.0D0)-t56*t74*t87*t167*(1.0D0/2.0D0)+t61*t87*t88*t167*(1&
     &.0D0/2.0D0)+t32*t64*t73*t88*t158*(1.0D0/4.0D0))
      d3(2) = -t93*(t284+t64*t66*t173-t64*t88*t175+t50*t66*t87*t179*(1&
     &.0D0/2.0D0)-t56*t74*t87*t179*(1.0D0/2.0D0)+t61*t87*t88*t179*(1.0D0&
     &/2.0D0)+t28*t64*t73*t74*t170*(1.0D0/4.0D0))
      d3(3) = t93*(t297+t64*t74*t184-t64*t88*t185+t50*t66*t87*t190*(1.&
     &0D0/2.0D0)-t56*t74*t87*t190*(1.0D0/2.0D0)+t61*t87*t88*t190*(1.0D0/&
     &2.0D0)+t7*t64*t66*t73*t183*(1.0D0/4.0D0))
      d4(1) = -t93*(-t64*t66*t201+t64*t74*t203+t50*t66*t87*t208*(1.0D&
     &0/2.0D0)-t56*t74*t87*t208*(1.0D0/2.0D0)+t61*t87*t88*t208*(1.0D0/2.&
     &0D0)+t45*t64*t88*t198*t199*(1.0D0/4.0D0))
      d4(2) = -t93*(t64*t66*t216-t64*t88*t218+t50*t66*t87*t222*(1.0D0&
     &/2.0D0)-t56*t74*t87*t222*(1.0D0/2.0D0)+t61*t87*t88*t222*(1.0D0/2.0&
     &D0)+t41*t64*t74*t199*t212*(1.0D0/4.0D0))
      d4(3) = t93*(t64*t74*t227-t64*t88*t228+t50*t66*t87*t233*(1.0D0/&
     &2.0D0)-t56*t74*t87*t233*(1.0D0/2.0D0)+t61*t87*t88*t233*(1.0D0/2.0D&
     &0)+t18*t64*t66*t199*t226*(1.0D0/4.0D0))
      d5(1) = t93*(-t64*t66*t239+t64*t74*t241+t50*t66*t87*t246*(1.0D0&
     &/2.0D0)-t56*t74*t87*t246*(1.0D0/2.0D0)+t61*t87*t88*t246*(1.0D0/2.0&
     &D0)+t45*t64*t88*t199*t237*(1.0D0/4.0D0))
      d5(2) = t93*(t64*t66*t252-t64*t88*t254+t50*t66*t87*t258*(1.0D0/&
     &2.0D0)-t56*t74*t87*t258*(1.0D0/2.0D0)+t61*t87*t88*t258*(1.0D0/2.0D&
     &0)+t41*t64*t74*t199*t249*(1.0D0/4.0D0))
      d5(3) = -t93*(t64*t74*t263-t64*t88*t264+t50*t66*t87*t269*(1.0D0&
     &/2.0D0)-t56*t74*t87*t269*(1.0D0/2.0D0)+t61*t87*t88*t269*(1.0D0/2.0&
     &D0)+t18*t64*t66*t199*t262*(1.0D0/4.0D0))
      d6(1) = -t93*(t61*t64-t64*t66*t275+t64*t74*t277+t50*t66*t87*t28&
     &2*(1.0D0/2.0D0)-t56*t74*t87*t282*(1.0D0/2.0D0)+t61*t87*t88*t282*(1&
     &.0D0/2.0D0)+t45*t64*t88*t199*t273*(1.0D0/4.0D0))
      d6(2) = -t93*(-t284+t64*t66*t289-t64*t88*t291+t50*t66*t87*t295*&
     &(1.0D0/2.0D0)-t56*t74*t87*t295*(1.0D0/2.0D0)+t61*t87*t88*t295*(1.0&
     &D0/2.0D0)+t41*t64*t74*t199*t286*(1.0D0/4.0D0))
      d6(3) = t93*(-t297+t64*t74*t301-t64*t88*t302+t50*t66*t87*t307*(&
     &1.0D0/2.0D0)-t56*t74*t87*t307*(1.0D0/2.0D0)+t61*t87*t88*t307*(1.0D&
     &0/2.0D0)+t18*t64*t66*t199*t300*(1.0D0/4.0D0))

end subroutine rise_d
!=============================================================================
!============================ DF_COM_DISTANCE ================================
!=============================================================================
! fpan3-at-ncsu-dot-edu, 03/30/2016
!
! DF_COM_DISTANCE -- difference of distances between centers of mass of two atom groups
!
! input: i = (a1, ..., aN, 0, b1, ..., bM, 0, c1, ..., ck, d1, ..., dL)
! last zero is optional; a?/b?/c?/d? -- indices of the participating atoms

function v_DF_COM_DISTANCE(cv, x) result(value)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   double precision :: value

   type(colvar_t), intent(in) :: cv
   double precision, intent(in) :: x(*)

   integer :: n

   double precision :: cm1(3), cm2(3), cm3(3), cm4(4)

   nfe_assert(cv%type == COLVAR_DF_COM_DISTANCE)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).gt.6)
   nfe_assert(associated(cv%r))
   nfe_assert(size(cv%r).eq.size(cv%i))

#ifdef MPI
   if (mytaskid.eq.0) then
#endif /* MPI */
   n = 1
   call group_com(cv, x, n, cm1)
   nfe_assert(n.le.size(cv%i))
   call group_com(cv, x, n, cm2)
   nfe_assert(n.le.size(cv%i))   
   call group_com(cv, x, n, cm3)
   nfe_assert(n.le.size(cv%i))   
   call group_com(cv, x, n, cm4)

   value = distance(cm1, cm2) - distance(cm3, cm4)
#ifdef MPI
   else
      value = ZERO
   end if
#endif /* MPI */

end function v_DF_COM_DISTANCE

!=============================================================================

subroutine f_DF_COM_DISTANCE(cv, x, fcv, f)

   NFE_USE_AFAILED

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   double precision, intent(in) :: x(*), fcv
   double precision, intent(inout) :: f(*)

   double precision :: d1(3), d2(3), d3(3), d4(3)
   double precision :: cm1(3), cm2(3), cm3(3), cm4(3)

   integer :: n
#ifdef MPI
   integer :: error
#endif /* MPI */

  nfe_assert(cv%type == COLVAR_DF_COM_DISTANCE)

   nfe_assert(associated(cv%i))
   nfe_assert(size(cv%i).gt.6)
   nfe_assert(associated(cv%r))
   nfe_assert(size(cv%r).eq.size(cv%i))

   NFE_MASTER_ONLY_BEGIN

        n = 1
        call group_com(cv, x, n, cm1)
        call group_com(cv, x, n, cm2)
        call group_com(cv, x, n, cm3)
        call group_com(cv, x, n, cm4)

        if ((distance(cm1, cm2).lt.1.0d-8).or.(distance(cm3, cm4).lt.1.0d-8)) then
        write (unit = ERR_UNIT, fmt = '(/a,a/)') NFE_ERROR, &
        'DF_COM_DISTANCE : the distance got too small'
        call terminate()
        end if ! distance too small

        call distance_d(cm1, cm2, d1, d2)
        call distance_d(cm3, cm4, d3, d4)
        d3 = -d3
        d4 = -d4

        d1 = fcv*d1
        d2 = fcv*d2
        d3 = fcv*d3
        d4 = fcv*d4

        n = 1
        call group_com_d(cv, f, d1, n)
        call group_com_d(cv, f, d2, n)
        call group_com_d(cv, f, d3, n)
call group_com_d(cv, f, d4, n)

        NFE_MASTER_ONLY_END

#ifdef MPI
   call mpi_bcast(f(1:3*pmemd_natoms()), 3*pmemd_natoms(), MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

end subroutine f_DF_COM_DISTANCE

!=============================================================================

subroutine b_DF_COM_DISTANCE(cv, cvno, amass)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)

   integer :: error, ngroups

   nfe_assert(cv%type == COLVAR_DF_COM_DISTANCE)

   call com_check_i(cv%i, cvno, 'DF_COM_DISTANCE', ngroups)
   if (ngroups.ne.4) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (DF_COM_DISTANCE) : number of atom groups is not four'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! ngroups.ne.4

   if (associated(cv%r)) &
      deallocate(cv%r)

   allocate(cv%r(size(cv%i)), stat = error)
   if (error.ne.0) &
      NFE_OUT_OF_MEMORY

   cv%r = ZERO
   do error = 1, size(cv%i)
      if (cv%i(error).gt.0) &
         cv%r(error) = amass(cv%i(error))
   end do

   call com_init_weights(cv%r)

end subroutine b_DF_COM_DISTANCE

!=============================================================================

subroutine p_DF_COM_DISTANCE(cv, lun)

#ifndef NFE_DISABLE_ASSERT
   use nfe_lib_mod
#endif /* NFE_DISABLE_ASSERT */

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_DF_COM_DISTANCE)
   nfe_assert(associated(cv%i))

   call com_print_i(cv%i, lun)

end subroutine p_DF_COM_DISTANCE
!=============================================================================
!============================ ORIENTATION_ANGLE ==============================
!=============================================================================
!
! input:
!
! cv%i = (a1, a2, ..., aN)
!
!     (list of participating atoms)
!
! cv%r = (a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z, ...)
!
!        (reference coordinates)
!
! value = orientation angle: 2*acos(q0),
! in which (q0,q1,q2,q3) is orientation quaternion representing optimum rotation w.r.t. reference
!

function v_ORIENTATION_ANGLE(cv, x) result(value)

   NFE_USE_AFAILED

   use nfe_lib_mod
   use gbl_constants_mod
   use parallel_dat_mod

   implicit none

   double precision :: value

   type(colvar_t), intent(inout) :: cv

   double precision, intent(in) :: x(*)


   type(priv_t), pointer :: priv

   double precision :: cm(3), cm_nrm
   integer :: i, n

   integer :: a, a3

   nfe_assert(cv%type == COLVAR_ORIENTATION_ANGLE)

   priv => get_priv(cv)
   nfe_assert(associated(priv))
   nfe_assert(associated(priv%mass))
   nfe_assert(size(priv%mass).gt.0)

   priv%value = ZERO
   cm = ZERO
   n = 1
   do a = 1, priv%n_atoms
      a3 = 3*cv%i(a)
      cm = cm + priv%mass(n)*x(a3 - 2:a3)
      n = n + 1
   end do
   cm = cm/priv%total_mass
   cm_nrm = ZERO
   n = 1
   do a = 1, priv%n_atoms
      a3 = 3*(n - 1)
      do i = 1, 3
         priv%cm_crd(a3 + i) = x(3*(cv%i(a) - 1) + i) - cm(i)
         cm_nrm = cm_nrm + priv%mass(n)*priv%cm_crd(a3 + i)**2
      end do
      n = n + 1
   end do

   call orientation_q(priv%n_atoms, priv%mass, &
         priv%cm_crd, priv%ref_crd, priv%lambda, priv%quaternion)
   priv%value = (180.0/PI) * 2.0 * acos(abs(priv%quaternion(1,1)))

   value = priv%value
end function v_ORIENTATION_ANGLE
!=============================================================================
subroutine f_ORIENTATION_ANGLE(cv, fcv, f)

   use nfe_lib_mod
   use gbl_constants_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: fcv
   double precision, intent(inout) :: f(*)

   integer :: n, a, a3, f3
   double precision :: dSx(4,4), dSy(4,4), dSz(4,4), QL(3)
   double precision :: dSxq1(4), dSyq1(4), dSzq1(4)
   double precision :: a1x, a1y, a1z, tmp
   type(priv_t), pointer :: priv

#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_ORIENTATION_ANGLE)
   nfe_assert(associated(cv%i))
   priv => get_priv(cv)

   nfe_assert(associated(priv))
   nfe_assert(associated(priv%mass))
   nfe_assert(size(priv%mass).gt.0)

#ifdef MPI
   call mpi_bcast(priv%value, 1, MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

NFE_MASTER_ONLY_BEGIN
   if ( priv%quaternion(1,1)**2 < ONE ) then
      tmp = (180.0/PI)*(-2.0)* fcv / &
             max(sqrt(ONE - (priv%quaternion(1,1))**2),dble(0.00000001))
   else
      tmp = ZERO
   end if

      QL(1) = tmp * priv%quaternion(1,2) / &
                    max(priv%lambda(1)-priv%lambda(2),dble(0.00000001))
      QL(2) = tmp * priv%quaternion(1,3) / &
                    max(priv%lambda(1)-priv%lambda(3),dble(0.00000001))
      QL(3) = tmp * priv%quaternion(1,4) / &
                    max(priv%lambda(1)-priv%lambda(4),dble(0.00000001))
      n = 1
      do a = 1, priv%n_atoms
         a3 = 3*n
         f3 = 3*cv%i(a)

         a1x = priv%ref_crd(a3 - 2)
         a1y = priv%ref_crd(a3 - 1)
         a1z = priv%ref_crd(a3)

         dSx(1,1) = a1x
         dSy(1,1) = a1y
         dSz(1,1) = a1z

         dSx(2,1) = 0.0
         dSy(2,1) = -a1z
         dSz(2,1) = a1y

         dSx(1,2) = 0.0
         dSy(1,2) = -a1z
         dSz(1,2) = a1y

         dSx(3,1) = a1z
         dSy(3,1) = 0.0
         dSz(3,1) = -a1x

         dSx(1,3) = a1z
         dSy(1,3) = 0.0
         dSz(1,3) = -a1x

         dSx(4,1) = -a1y
         dSy(4,1) = a1x
         dSz(4,1) = 0.0

         dSx(1,4) = -a1y
         dSy(1,4) = a1x
         dSz(1,4) = 0.0

         dSx(2,2) = a1x
         dSy(2,2) = -a1y
         dSz(2,2) = -a1z

         dSx(3,2) = a1y
         dSy(3,2) = a1x
         dSz(3,2) = 0.0

         dSx(2,3) = a1y
         dSy(2,3) = a1x
         dSz(2,3) = 0.0

         dSx(4,2) = a1z
         dSy(4,2) = 0.0
         dSz(4,2) = a1x

         dSx(2,4) = a1z
         dSy(2,4) = 0.0
         dSz(2,4) = a1x

         dSx(3,3) = -a1x
         dSy(3,3) = a1y
         dSz(3,3) = -a1z

         dSx(4,3) = 0.0
         dSy(4,3) = a1z
         dSz(4,3) = a1y

         dSx(3,4) = 0.0
         dSy(3,4) = a1z
         dSz(3,4) = a1y

         dSx(4,4) = -a1x
         dSy(4,4) = -a1y
         dSz(4,4) = a1z

         dSxq1 = matmul(dSx,priv%quaternion(:,1))
         dSyq1 = matmul(dSy,priv%quaternion(:,1))
         dSzq1 = matmul(dSz,priv%quaternion(:,1))
         f(f3 - 2) = f(f3 - 2) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSxq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSxq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSxq1,priv%quaternion(:,4)) )
         f(f3 - 1) = f(f3 - 1) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSyq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSyq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSyq1,priv%quaternion(:,4)) )
         f(f3) = f(f3) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSzq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSzq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSzq1,priv%quaternion(:,4)) )
         n = n + 1
      end do
      NFE_MASTER_ONLY_END
end subroutine f_ORIENTATION_ANGLE
!============================================================================
subroutine b_ORIENTATION_ANGLE(cv, cvno, amass)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)
   
   double precision, pointer :: coor(:) => null()

   integer :: i, n_atoms, error
   integer :: a, b, a1
   double precision :: cm(3)
   
   type(priv_t), pointer :: priv


   nfe_assert(cv%type == COLVAR_ORIENTATION_ANGLE)

   ! very basic checks
   !call check_i(cv%i, cvno, 'COLVAR_ORIENTATION_ANGLE')
   if (.not. associated(cv%i)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION_ANGLE) : no integers found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cv%i)

   if (cv_nr.gt.0) then
    if (.not. associated(cv%r)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_INFO, 'CV #', cvno, &
            ' (ORIENTATION_ANGLE) : no reals (reference coordinates) found'
      NFE_MASTER_ONLY_END
      call terminate()
    end if
   end if 

   ! count the number of atoms
   n_atoms = 0

   do i = 1, size(cv%i)
      if (cv%i(i) == 0) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(i)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION_ANGLE) : unexpected zero (integer #', i, ')'
            NFE_MASTER_ONLY_END
            call terminate()
      else if (cv%i(i) .lt. 0) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION_ANGLE) : negative integer'
         NFE_MASTER_ONLY_END
         call terminate()
      else
         n_atoms = n_atoms + 1
      end if
   end do

   if (cv_nr.eq.0) then
      if (associated(coor)) &
       deallocate(coor)
       allocate(coor(3*pmemd_natoms()), stat = error)
        if (error.ne.0) &
         NFE_OUT_OF_MEMORY
         call colvar_read_refcrd(coor, refcrd_file)
         allocate(cv%r(3*size(cv%i)), stat = error)
          if (error.ne.0) &
           NFE_OUT_OF_MEMORY
            do a = 1, n_atoms
               a1=cv%i(a)
               cv%r(1 + 3*(a - 1): 3*a) = coor(1 + 3*(a1 - 1): 3*a1)
            end do
       deallocate(coor)
   end if

   if (size(cv%r) /= 3*n_atoms) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION_ANGLE) : wrong number of reals'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! size(cv%r) /= 3*n_atoms

   ! allocate priv_t instance for this variable
   priv => new_priv(cv)

   ! basic checks
   if (n_atoms < 3) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, &
            fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
            (1)//',a,'//pfmt(n_atoms)//',a/)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (ORIENTATION_ANGLE) : too few integers in reference (', 1, ':', n_atoms, ')'
      NFE_MASTER_ONLY_END
      call terminate()
   end if

   do a = 1, n_atoms
      if (cv%i(a) < 1 .or. cv%i(a) > pmemd_natoms()) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
               (cv%i(a))//',a,'//pfmt(pmemd_natoms())//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION_ANGLE) : integer #', a, ' (', cv%i(a), &
               ') is out of range [1, ', pmemd_natoms(), ']'
         NFE_MASTER_ONLY_END
         call terminate()
      end if
   end do

   ! check for duplicates
   do a = 1, n_atoms
      do b = a + 1, n_atoms
         if (cv%i(a) == cv%i(b)) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
                  (b)//',a,'//pfmt(cv%i(a))//',a/)') &
                  NFE_ERROR, 'CV #', cvno, &
                  ' (ORIENTATION_ANGLE) : integers #', a, ' and #', b, &
                  ' are equal (', cv%i(a), ')'
            NFE_MASTER_ONLY_END
            call terminate()
         end if ! cv%i(a) == cv%i(b)
      end do
   end do

   ! allocate/setup

   priv%n_atoms = n_atoms

   allocate(priv%mass(n_atoms), priv%cm_crd(3*n_atoms), &
            priv%ref_crd(3*n_atoms), stat = error)

   if (error /= 0) &
      NFE_OUT_OF_MEMORY

   cm = ZERO
   priv%total_mass = ZERO

   do a = 1, n_atoms
      priv%mass(a) = amass(cv%i(a))
      priv%total_mass = priv%total_mass + priv%mass(a)
      cm = cm + priv%mass(a)*cv%r(1 + 3*(a - 1): 3*a)
   end do
   cm = cm/priv%total_mass
   priv%ref_nrm = ZERO
   do a = 1, n_atoms
      do b = 1, 3
         priv%ref_crd(3*(a - 1) + b) = cv%r(3*(a - 1) + b) - cm(b)
         priv%ref_nrm = priv%ref_nrm + priv%mass(a)*priv%ref_crd(3*(a - 1) + b)**2
      end do
   end do

end subroutine b_ORIENTATION_ANGLE
!=============================================================================

subroutine c_ORIENTATION_ANGLE(cv)

   NFE_USE_AFAILED

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(inout) :: cv

   type(priv_t), pointer :: priv

   nfe_assert(cv%type == COLVAR_ORIENTATION_ANGLE)
   nfe_assert(associated(cv%i).and.size(cv%i).gt.0)

   priv => get_priv(cv)

   nfe_assert(associated(priv%mass))
   nfe_assert(associated(priv%cm_crd))
   nfe_assert(associated(priv%ref_crd))

   deallocate(priv%mass, priv%cm_crd, priv%ref_crd)

   call del_priv(cv)

end subroutine c_ORIENTATION_ANGLE

!=============================================================================
subroutine p_ORIENTATION_ANGLE(cv, lun)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   type(priv_t), pointer :: priv
   integer :: a, c
   character(4) :: aname

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_ORIENTATION_ANGLE)
   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   priv => get_priv(cv)
   nfe_assert(priv%n_atoms.gt.0.and.priv%n_atoms.eq.size(cv%i))

   write (unit = lun, fmt = '(a,a)', advance = 'NO') NFE_INFO, 'atoms = ('

   c = 1
   do a = 1, priv%n_atoms

      nfe_assert(cv%i(a) > 0 .and. cv%i(a) <= pmemd_natoms())
      aname = pmemd_atom_name(cv%i(a))

      write (unit = lun, fmt = '('//pfmt(cv%i(a))//',a,a,a)', advance = 'NO') &
         cv%i(a), ' [', trim(aname), ']'

      if (a == priv%n_atoms) then
         write (unit = lun, fmt = '(a)') ')'
      else if (mod(c, 5) == 0) then
         write (unit = lun, fmt = '(a,/a,3x)', advance = 'NO') ',', NFE_INFO
      else
         write (unit = lun, fmt = '(a)', advance = 'NO') ', '
      end if

      c = c + 1
   end do

   if (cv_nr.eq.0) &
      write (unit = lun, fmt = '(a,a,a)') &
         NFE_INFO, 'reference coordinates are loaded from : ', trim(refcrd_file)

   write (unit = lun, fmt = '(a,a)') &
         NFE_INFO, 'reference coordinates :'

   c = 1
   do a = 1, priv%n_atoms
      nfe_assert(c + 2 <= size(cv%r))
      write (unit = lun, fmt = '(a,3x,i5,a,f8.3,a,f8.3,a,f8.3)') &
         NFE_INFO, cv%i(a), ' : ', &
         cv%r(c), ', ', cv%r(c + 1), &
         ', ', cv%r(c + 2)
      c = c + 3
   end do

end subroutine p_ORIENTATION_ANGLE
!=========================================================================
!========================== ORIENTATION_PROJ =============================
!=========================================================================
!
! input:
!
! cv%i = (a1, a2, ..., aN)
!
!     (list of participating atoms)
!
! cv%r = (a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z, ...)
!
!        (reference coordinates)
!
! value = orientation proj: cos(theta)=2q0**2-1,
! in which (q0,q1,q2,q3) is orientation quaternion representing optimum rotation w.r.t. reference
!

!============================================================================
function v_ORIENTATION_PROJ(cv, x) result(value)
   NFE_USE_AFAILED

   use nfe_lib_mod
   use gbl_constants_mod
   use parallel_dat_mod

   implicit none

   double precision :: value

   type(colvar_t), intent(inout) :: cv

   double precision, intent(in) :: x(*)



   type(priv_t), pointer :: priv

#ifdef MPI
   integer :: error
#endif /* MPI */

   double precision :: cm(3), cm_nrm
   integer :: i, n
! gfortran 4.2 and lower does not support the volatile keyword, but this is
! necessary to work around an intel compiler optimization bug. So turn off
! the volatile keyword for qualifying compilers
#if !defined(__GNUC__) || (__GNUC__ >= 4 && __GNUC_MINOR__ > 2)
   integer, volatile :: a, a3 ! work around ifort optimization bug
#else
   integer :: a, a3
#endif /* __GFORTRAN__ */

   nfe_assert(cv%type == COLVAR_ORIENTATION_PROJ)

   priv => get_priv(cv)

   nfe_assert(associated(priv))
   nfe_assert(associated(priv%mass))
   nfe_assert(size(priv%mass).gt.0)

   priv%value = ZERO

   cm = ZERO
   n = 1
   do a = 1, priv%n_atoms
      a3 = 3*cv%i(a)
      cm = cm + priv%mass(n)*x(a3 - 2:a3)
      n = n + 1
   end do

   cm = cm/priv%total_mass

   ! populate cm_crd && compute the "norm"
   cm_nrm = ZERO
   n = 1
   do a = 1, priv%n_atoms
      a3 = 3*(n - 1)
      do i = 1, 3
         priv%cm_crd(a3 + i) = x(3*(cv%i(a) - 1) + i) - cm(i)
         cm_nrm = cm_nrm + priv%mass(n)*priv%cm_crd(a3 + i)**2
      end do
      n = n + 1
   end do

   call orientation_q(priv%n_atoms, priv%mass, &
         priv%cm_crd, priv%ref_crd, priv%lambda, priv%quaternion)

!#ifdef MPI
!   call mpi_reduce(priv%value, value, 1, MPI_DOUBLE_PRECISION, &
!                   MPI_SUM, 0, pmemd_comm, error)
!   nfe_assert(error.eq.0)
!#else
   priv%value = 2.0 * priv%quaternion(1,1) * priv%quaternion(1,1) - 1.0
!#endif /* MPI */

   value = priv%value

end function v_ORIENTATION_PROJ

!=============================================================================

!
! assumes that atom positions have not been
!   changed since last call to 'value()'
!

subroutine f_ORIENTATION_PROJ(cv, fcv, f)

   use nfe_lib_mod
   use gbl_constants_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: fcv
   double precision, intent(inout) :: f(*)



   integer :: n, a, a3, f3
   double precision :: dSx(4,4), dSy(4,4), dSz(4,4), QL(3)
   double precision :: dSxq1(4), dSyq1(4), dSzq1(4)
   double precision :: a1x, a1y, a1z, tmp
   type(priv_t), pointer :: priv

!#ifdef MPI
   integer :: error
!#endif /* MPI */

   nfe_assert(cv%type == COLVAR_ORIENTATION_PROJ)
   nfe_assert(associated(cv%i))

   priv => get_priv(cv)

   nfe_assert(associated(priv))
   nfe_assert(associated(priv%mass))
   nfe_assert(size(priv%mass).gt.0)

#ifdef MPI
   call mpi_bcast(priv%value, 1, MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

      tmp = 2.0 * 2.0 * (priv%quaternion(1,1)) * fcv                                         !w/o any condition!

      QL(1) = tmp * priv%quaternion(1,2) / &
                    max(priv%lambda(1)-priv%lambda(2),dble(0.00000001))
      QL(2) = tmp * priv%quaternion(1,3) / &
                    max(priv%lambda(1)-priv%lambda(3),dble(0.00000001))
      QL(3) = tmp * priv%quaternion(1,4) / &
                    max(priv%lambda(1)-priv%lambda(4),dble(0.00000001))
      n = 1
      do a = 1, priv%n_atoms
         a3 = 3*n
         f3 = 3*cv%i(a)

         a1x = priv%ref_crd(a3 - 2)
         a1y = priv%ref_crd(a3 - 1)
         a1z = priv%ref_crd(a3)

         dSx(1,1) = a1x
         dSy(1,1) = a1y
         dSz(1,1) = a1z

         dSx(2,1) = 0.0
         dSy(2,1) = -a1z
         dSz(2,1) = a1y

         dSx(1,2) = 0.0
         dSy(1,2) = -a1z
         dSz(1,2) = a1y

         dSx(3,1) = a1z
         dSy(3,1) = 0.0
         dSz(3,1) = -a1x

         dSx(1,3) = a1z
         dSy(1,3) = 0.0
         dSz(1,3) = -a1x

         dSx(4,1) = -a1y
         dSy(4,1) = a1x
         dSz(4,1) = 0.0

         dSx(1,4) = -a1y
         dSy(1,4) = a1x
         dSz(1,4) = 0.0

         dSx(2,2) = a1x
         dSy(2,2) = -a1y
         dSz(2,2) = -a1z

         dSx(3,2) = a1y
         dSy(3,2) = a1x
         dSz(3,2) = 0.0

         dSx(2,3) = a1y
         dSy(2,3) = a1x
         dSz(2,3) = 0.0

         dSx(4,2) = a1z
         dSy(4,2) = 0.0
         dSz(4,2) = a1x

         dSx(2,4) = a1z
         dSy(2,4) = 0.0
         dSz(2,4) = a1x

         dSx(3,3) = -a1x
         dSy(3,3) = a1y
         dSz(3,3) = -a1z

         dSx(4,3) = 0.0
         dSy(4,3) = a1z
         dSz(4,3) = a1y

         dSx(3,4) = 0.0
         dSy(3,4) = a1z
         dSz(3,4) = a1y

         dSx(4,4) = -a1x
         dSy(4,4) = -a1y
         dSz(4,4) = a1z

         dSxq1 = matmul(dSx,priv%quaternion(:,1))
         dSyq1 = matmul(dSy,priv%quaternion(:,1))
         dSzq1 = matmul(dSz,priv%quaternion(:,1))

         f(f3 - 2) = f(f3 - 2) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSxq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSxq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSxq1,priv%quaternion(:,4)) )

         f(f3 - 1) = f(f3 - 1) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSyq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSyq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSyq1,priv%quaternion(:,4)) )

         f(f3) = f(f3) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSzq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSzq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSzq1,priv%quaternion(:,4)) )

         n = n + 1
      end do

end subroutine f_ORIENTATION_PROJ

!=============================================================================

subroutine b_ORIENTATION_PROJ(cv, cvno, amass)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)
   
   double precision, pointer :: coor(:) => null()

   integer :: i, n_atoms, error
   integer :: a, b, a1
   double precision :: cm(3)

   type(priv_t), pointer :: priv



   nfe_assert(cv%type == COLVAR_ORIENTATION_PROJ)

   ! very basic checks
   if (.not. associated(cv%i)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : no integers found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cv%i)

   if (cv_nr.gt.0) then
    if (.not. associated(cv%r)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_INFO, 'CV #', cvno, &
            ' (ORIENTATION) : no reals (reference coordinates) found'
      NFE_MASTER_ONLY_END
      call terminate()
    end if
   end if 

   ! count the number of atoms
   n_atoms = 0

   do i = 1, size(cv%i)
      if (cv%i(i) == 0) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(i)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : unexpected zero (integer #', i, ')'
            NFE_MASTER_ONLY_END
            call terminate()
      else if (cv%i(i) .lt. 0) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : negative integer'
         NFE_MASTER_ONLY_END
         call terminate()
      else
         n_atoms = n_atoms + 1
      end if
   end do

   if (cv_nr.eq.0) then
      if (associated(coor)) &
       deallocate(coor)
       allocate(coor(3*pmemd_natoms()), stat = error)
        if (error.ne.0) &
         NFE_OUT_OF_MEMORY
         call colvar_read_refcrd(coor, refcrd_file)
         allocate(cv%r(3*size(cv%i)), stat = error)
          if (error.ne.0) &
           NFE_OUT_OF_MEMORY
            do a = 1, n_atoms
               a1=cv%i(a)
               cv%r(1 + 3*(a - 1): 3*a) = coor(1 + 3*(a1 - 1): 3*a1)
            end do
       deallocate(coor)
   end if

   if (size(cv%r) /= 3*n_atoms) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : wrong number of reals'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! size(cv%r) /= 3*n_atoms

   ! allocate priv_t instance for this variable
   priv => new_priv(cv)

   ! basic checks
   if (n_atoms < 3) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, &
            fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
            (1)//',a,'//pfmt(n_atoms)//',a/)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (ORIENTATION) : too few integers in reference (', 1, ':', n_atoms, ')'
      NFE_MASTER_ONLY_END
      call terminate()
   end if

   do a = 1, n_atoms
      if (cv%i(a) < 1 .or. cv%i(a) > pmemd_natoms()) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
               (cv%i(a))//',a,'//pfmt(pmemd_natoms())//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : integer #', a, ' (', cv%i(a), &
               ') is out of range [1, ', pmemd_natoms(), ']'
         NFE_MASTER_ONLY_END
         call terminate()
      end if
   end do

   ! check for duplicates
   do a = 1, n_atoms
      do b = a + 1, n_atoms
         if (cv%i(a) == cv%i(b)) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
                  (b)//',a,'//pfmt(cv%i(a))//',a/)') &
                  NFE_ERROR, 'CV #', cvno, &
                  ' (ORIENTATION) : integers #', a, ' and #', b, &
                  ' are equal (', cv%i(a), ')'
            NFE_MASTER_ONLY_END
            call terminate()
         end if ! cv%i(a) == cv%i(b)
      end do
   end do

   ! allocate/setup

   priv%n_atoms = n_atoms

   allocate(priv%mass(n_atoms), priv%cm_crd(3*n_atoms), &
            priv%ref_crd(3*n_atoms), stat = error)

   if (error /= 0) &
      NFE_OUT_OF_MEMORY

   cm = ZERO
   priv%total_mass = ZERO

   do a = 1, n_atoms
      priv%mass(a) = amass(cv%i(a))
      priv%total_mass = priv%total_mass + priv%mass(a)
      cm = cm + priv%mass(a)*cv%r(1 + 3*(a - 1): 3*a)
   end do

   cm = cm/priv%total_mass

   ! translate reference coordinates to CM frame
   priv%ref_nrm = ZERO
   do a = 1, n_atoms
      do b = 1, 3
         priv%ref_crd(3*(a - 1) + b) = cv%r(3*(a - 1) + b) - cm(b)
         priv%ref_nrm = priv%ref_nrm + priv%mass(a)*priv%ref_crd(3*(a - 1) + b)**2
      end do
   end do


end subroutine b_ORIENTATION_PROJ

!=============================================================================

subroutine c_ORIENTATION_PROJ(cv)

   NFE_USE_AFAILED

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(inout) :: cv

   type(priv_t), pointer :: priv

   nfe_assert(cv%type == COLVAR_ORIENTATION_PROJ)
   nfe_assert(associated(cv%i).and.size(cv%i).gt.0)

   priv => get_priv(cv)

   nfe_assert(associated(priv%mass))
   nfe_assert(associated(priv%cm_crd))
   nfe_assert(associated(priv%ref_crd))

   deallocate(priv%mass, priv%cm_crd, priv%ref_crd)

   call del_priv(cv)

end subroutine c_ORIENTATION_PROJ

!=============================================================================

subroutine p_ORIENTATION_PROJ(cv, lun)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   type(priv_t), pointer :: priv
   integer :: a, c
   character(4) :: aname

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_ORIENTATION_PROJ)
   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   priv => get_priv(cv)
   nfe_assert(priv%n_atoms.gt.0.and.priv%n_atoms.eq.size(cv%i))

   write (unit = lun, fmt = '(a,a)', advance = 'NO') NFE_INFO, 'atoms = ('

   c = 1
   do a = 1, priv%n_atoms

      nfe_assert(cv%i(a) > 0 .and. cv%i(a) <= pmemd_natoms())
      aname = pmemd_atom_name(cv%i(a))

      write (unit = lun, fmt = '('//pfmt(cv%i(a))//',a,a,a)', advance = 'NO') &
         cv%i(a), ' [', trim(aname), ']'

      if (a == priv%n_atoms) then
         write (unit = lun, fmt = '(a)') ')'
      else if (mod(c, 5) == 0) then
         write (unit = lun, fmt = '(a,/a,3x)', advance = 'NO') ',', NFE_INFO
      else
         write (unit = lun, fmt = '(a)', advance = 'NO') ', '
      end if

      c = c + 1
   end do

   if (cv_nr.eq.0) &
      write (unit = lun, fmt = '(a,a,a)') &
         NFE_INFO, 'reference coordinates are loaded from : ', trim(refcrd_file)

   write (unit = lun, fmt = '(a,a)') &
         NFE_INFO, 'reference coordinates :'

   c = 1
   do a = 1, priv%n_atoms
      nfe_assert(c + 2 <= size(cv%r))
      write (unit = lun, fmt = '(a,3x,i5,a,f8.3,a,f8.3,a,f8.3)') &
         NFE_INFO, cv%i(a), ' : ', &
         cv%r(c), ', ', cv%r(c + 1), &
         ', ', cv%r(c + 2)
      c = c + 3
   end do

end subroutine p_ORIENTATION_PROJ
!=========================================================================
!========================== TILT  ========================================
!=========================================================================

!
! input:
!
! cv%i = (a1, a2, ..., aN)
!
!     (list of participating atoms)
!
! cv%r = (a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z, ...)
!
!        (reference coordinates)
!
! axis = (e1, e2, e3)
!      (axis which defined by user)
! value = tilt[-1,1]: t = cos(w)= 2*(q0/cos(alpha/2))**2-1
! where alpha=2*atan(q.e)/q0
! in which q=(q1,q2,q3) is orientation vector and e=axis
!

!=============================================================================

function v_TILT(cv, x) result(value)

   NFE_USE_AFAILED

   use nfe_lib_mod
   use gbl_constants_mod
   use parallel_dat_mod
   !use nfe_colvar_math, only : norm3

   implicit none

   double precision :: value

   type(colvar_t), intent(inout) :: cv

   double precision, intent(in) :: x(*)
   double precision  :: e(3)
   double precision :: get_vec(3)
   double precision :: alpha, beta, alpha1, cos_spin, cos_theta2, cos_w

   type(priv_t), pointer :: priv

#ifdef MPI
   integer :: error
#endif /* MPI */

   double precision :: cm(3), cm_nrm
   integer :: i, n
! gfortran 4.2 and lower does not support the volatile keyword, but this is
! necessary to work around an intel compiler optimization bug. So turn off
! the volatile keyword for qualifying compilers
#if !defined(__GNUC__) || (__GNUC__ >= 4 && __GNUC_MINOR__ > 2)
   integer, volatile :: a, a3 ! work around ifort optimization bug
#else
   integer :: a, a3
#endif /* __GFORTRAN__ */

   nfe_assert(cv%type == COLVAR_TILT)

   priv => get_priv(cv)

   nfe_assert(associated(priv))
   nfe_assert(associated(priv%mass))
   nfe_assert(size(priv%mass).gt.0)

   priv%value = ZERO

   cm = ZERO
   n = 1
   do a = 1, priv%n_atoms
      a3 = 3*cv%i(a)
      cm = cm + priv%mass(n)*x(a3 - 2:a3)
      n = n + 1
   end do

   cm = cm/priv%total_mass

   ! populate cm_crd && compute the "norm"
   cm_nrm = ZERO
   n = 1
   do a = 1, priv%n_atoms
      a3 = 3*(n - 1)
      do i = 1, 3
         priv%cm_crd(a3 + i) = x(3*(cv%i(a) - 1) + i) - cm(i)
         cm_nrm = cm_nrm + priv%mass(n)*priv%cm_crd(a3 + i)**2
      end do
      n = n + 1
   end do

   call orientation_q(priv%n_atoms, priv%mass, &
         priv%cm_crd, priv%ref_crd, priv%lambda, priv%quaternion)

   get_vec = [priv%quaternion(2,1), priv%quaternion(3,1), priv%quaternion(4,1)]
   if (norm3(cv%axis).ne.ZERO) then
       e = cv%axis/norm3(cv%axis)
   else
       e = [0.0, 0.0, 1.0]                          !defualt value for the axis: [0,0, 0,0, 1.0]
   end if
   alpha1 = 2.0 * atan(dot_product(get_vec, e) / priv%quaternion(1,1))
   cos_spin = cos (0.5 * alpha1)
   if (cos_spin.ne.ZERO) then
      cos_theta2 = priv%quaternion(1,1) / cos_spin
   else
      cos_theta2 = ZERO
   end if

   cos_w = 2.0 * (cos_theta2 *cos_theta2) - 1.0
   priv%value = cos_w

   value = priv%value

end function v_TILT

!=============================================================================

!
! assumes that atom positions have not been
!   changed since last call to 'value()'
!

subroutine f_TILT(cv, fcv, f)

   use nfe_lib_mod
   use gbl_constants_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: fcv
   double precision, intent(inout) :: f(*)
   double precision :: e(3)



   integer :: n, a, a3, f3, j
   double precision :: dSx(4,4), dSy(4,4), dSz(4,4), QL(3)
   double precision :: dSxq1(4), dSyq1(4), dSzq1(4)
   double precision :: a1x, a1y, a1z, tmp, tmp1, gamma, phi    !phi: spinAngle
   double precision :: get_vec(3)
   double precision :: alpha, beta, alpha1, cos_spin, cos_theta2, cos_w

   type(priv_t), pointer :: priv

#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_TILT)
   nfe_assert(associated(cv%i))

   priv => get_priv(cv)

   nfe_assert(associated(priv))
   nfe_assert(associated(priv%mass))
   nfe_assert(size(priv%mass).gt.0)

#ifdef MPI
   call mpi_bcast(priv%value, 1, MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */
   get_vec = [priv%quaternion(2,1), priv%quaternion(3,1), priv%quaternion(4,1)]
   if (norm3(cv%axis).ne.ZERO) then
       e = cv%axis/norm3(cv%axis)
   else
       e = [0.0, 0.0, 1.0]
   end if
   phi = atan (dot_product(get_vec, e) / priv%quaternion(1,1))
   beta = dot_product(get_vec, e)
   tmp = 4.0 * priv%quaternion(1,1) * beta / (priv%quaternion(1,1)**2 + beta**2)
   tmp1 = tmp*priv%quaternion(1,1)
   gamma = 4*priv%quaternion(1,1) - tmp * beta
      QL(1) = fcv * (gamma * priv%quaternion(1,2) + &
                     tmp1 * (e(1) * priv%quaternion(2,2) + &
                             e(2) * priv%quaternion(3,2) + &
                             e(3) * priv%quaternion(4,2))) / &
                    max(priv%lambda(1)-priv%lambda(2),dble(0.00000001)) * &
                    cos(phi)**2

      QL(2) = fcv * (gamma * priv%quaternion(1,3) + &
                     tmp1 * (e(1) * priv%quaternion(2,3) + &
                             e(2) * priv%quaternion(3,3) + &
                             e(3) * priv%quaternion(4,3))) / &
                    max(priv%lambda(1)-priv%lambda(3),dble(0.00000001)) * &
                    cos(phi)**2

      QL(3) = fcv * (gamma * priv%quaternion(1,4) + &
                     tmp1 * (e(1) * priv%quaternion(2,4) + &
                             e(2) * priv%quaternion(3,4) + &
                             e(3) * priv%quaternion(4,4))) / &
                    max(priv%lambda(1)-priv%lambda(4),dble(0.00000001)) * &
                    cos(phi)**2
      n = 1
      do a = 1, priv%n_atoms
         a3 = 3*n
         f3 = 3*cv%i(a)

         a1x = priv%ref_crd(a3 - 2)
         a1y = priv%ref_crd(a3 - 1)
         a1z = priv%ref_crd(a3)

         dSx(1,1) = a1x
         dSy(1,1) = a1y
         dSz(1,1) = a1z

         dSx(2,1) = 0.0
         dSy(2,1) = -a1z
         dSz(2,1) = a1y

         dSx(1,2) = 0.0
         dSy(1,2) = -a1z
         dSz(1,2) = a1y

         dSx(3,1) = a1z
         dSy(3,1) = 0.0
         dSz(3,1) = -a1x

         dSx(1,3) = a1z
         dSy(1,3) = 0.0
         dSz(1,3) = -a1x

         dSx(4,1) = -a1y
         dSy(4,1) = a1x
         dSz(4,1) = 0.0

         dSx(1,4) = -a1y
         dSy(1,4) = a1x
         dSz(1,4) = 0.0

         dSx(2,2) = a1x
         dSy(2,2) = -a1y
         dSz(2,2) = -a1z

         dSx(3,2) = a1y
         dSy(3,2) = a1x
         dSz(3,2) = 0.0

         dSx(2,3) = a1y
         dSy(2,3) = a1x
         dSz(2,3) = 0.0

         dSx(4,2) = a1z
         dSy(4,2) = 0.0
         dSz(4,2) = a1x

         dSx(2,4) = a1z
         dSy(2,4) = 0.0
         dSz(2,4) = a1x

         dSx(3,3) = -a1x
         dSy(3,3) = a1y
         dSz(3,3) = -a1z

         dSx(4,3) = 0.0
         dSy(4,3) = a1z
         dSz(4,3) = a1y

         dSx(3,4) = 0.0
         dSy(3,4) = a1z
         dSz(3,4) = a1y

         dSx(4,4) = -a1x
         dSy(4,4) = -a1y
         dSz(4,4) = a1z

         dSxq1 = matmul(dSx,priv%quaternion(:,1))
         dSyq1 = matmul(dSy,priv%quaternion(:,1))
         dSzq1 = matmul(dSz,priv%quaternion(:,1))

          f(f3 - 2) = f(f3 - 2) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSxq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSxq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSxq1,priv%quaternion(:,4)) )

          f(f3 - 1) = f(f3 - 1) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSyq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSyq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSyq1,priv%quaternion(:,4)) )

          f(f3) = f(f3) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSzq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSzq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSzq1,priv%quaternion(:,4)) )
         n = n + 1
      end do

end subroutine f_TILT

!=============================================================================

subroutine b_TILT(cv, cvno, amass)

   use nfe_lib_mod
   use parallel_dat_mod
   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)

   double precision, pointer :: coor(:) => null()

   integer :: i, n_atoms, error
   integer :: a, b, a1 
   double precision :: cm(3)

   type(priv_t), pointer :: priv



   nfe_assert(cv%type == COLVAR_TILT)

   ! very basic checks
   if (.not. associated(cv%i)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : no integers found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cv%i)

   if (.not. associated(cv%r)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_INFO, 'CV #', cvno, &
            ' (ORIENTATION) : no reals (reference coordinates) found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if

   if (cv_nr.eq.0) then
      if (associated(coor)) &
          deallocate(coor)
          allocate(coor(3*pmemd_natoms()), stat = error)
          if (error.ne.0) &
           NFE_OUT_OF_MEMORY
           call colvar_read_refcrd(coor, refcrd_file)
           allocate(cv%r(3*size(cv%i)), stat = error)
             if (error.ne.0) &
              NFE_OUT_OF_MEMORY
              do a = 1, n_atoms
               a1=cv%i(a)
               cv%r(1 + 3*(a - 1): 3*a) = coor(1 + 3*(a1 - 1): 3*a1)
              end do
          deallocate(coor)
   end if

   ! count the number of atoms
   n_atoms = 0

   do i = 1, size(cv%i)
      if (cv%i(i) == 0) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(i)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : unexpected zero (integer #', i, ')'
            NFE_MASTER_ONLY_END
            call terminate()
      else if (cv%i(i) .lt. 0) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : negative integer'
         NFE_MASTER_ONLY_END
         call terminate()
      else
         n_atoms = n_atoms + 1
      end if
   end do

   if (size(cv%r) /= 3*n_atoms) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : wrong number of reals'
      NFE_MASTER_ONLY_END
      call terminate()
   end if 

   ! allocate priv_t instance for this variable
   priv => new_priv(cv)

   ! basic checks
   if (n_atoms < 3) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, &
            fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
            (1)//',a,'//pfmt(n_atoms)//',a/)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (ORIENTATION) : too few integers in reference (', 1, ':', n_atoms, ')'
      NFE_MASTER_ONLY_END
      call terminate()
   end if

   do a = 1, n_atoms
      if (cv%i(a) < 1 .or. cv%i(a) > pmemd_natoms()) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
               (cv%i(a))//',a,'//pfmt(pmemd_natoms())//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : integer #', a, ' (', cv%i(a), &
               ') is out of range [1, ', pmemd_natoms(), ']'
         NFE_MASTER_ONLY_END
         call terminate()
      end if
   end do

   ! check for duplicates
   do a = 1, n_atoms
      do b = a + 1, n_atoms
         if (cv%i(a) == cv%i(b)) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
                  (b)//',a,'//pfmt(cv%i(a))//',a/)') &
                  NFE_ERROR, 'CV #', cvno, &
                  ' (ORIENTATION) : integers #', a, ' and #', b, &
                  ' are equal (', cv%i(a), ')'
            NFE_MASTER_ONLY_END
            call terminate()
         end if 
      end do
   end do

   ! allocate/setup

   priv%n_atoms = n_atoms

   allocate(priv%mass(n_atoms), priv%cm_crd(3*n_atoms), &
            priv%ref_crd(3*n_atoms), stat = error)

   if (error /= 0) &
      NFE_OUT_OF_MEMORY

   cm = ZERO
   priv%total_mass = ZERO

   do a = 1, n_atoms
      priv%mass(a) = amass(cv%i(a))
      priv%total_mass = priv%total_mass + priv%mass(a)
      cm = cm + priv%mass(a)*cv%r(1 + 3*(a - 1): 3*a)
   end do

   cm = cm/priv%total_mass

   ! translate reference coordinates to CM frame
   priv%ref_nrm = ZERO
   do a = 1, n_atoms
      do b = 1, 3
         priv%ref_crd(3*(a - 1) + b) = cv%r(3*(a - 1) + b) - cm(b)
         priv%ref_nrm = priv%ref_nrm + priv%mass(a)*priv%ref_crd(3*(a - 1) + b)**2
      end do
   end do


end subroutine b_TILT

!=============================================================================

subroutine c_TILT(cv)

   NFE_USE_AFAILED

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(inout) :: cv

   type(priv_t), pointer :: priv

   nfe_assert(cv%type == COLVAR_TILT)
   nfe_assert(associated(cv%i).and.size(cv%i).gt.0)

   priv => get_priv(cv)

   nfe_assert(associated(priv%mass))
   nfe_assert(associated(priv%cm_crd))
   nfe_assert(associated(priv%ref_crd))

   deallocate(priv%mass, priv%cm_crd, priv%ref_crd)

   call del_priv(cv)

end subroutine c_TILT

!=============================================================================

subroutine p_TILT(cv, lun)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   type(priv_t), pointer :: priv
   integer :: a, c
   character(4) :: aname

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_TILT)
   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   priv => get_priv(cv)
   nfe_assert(priv%n_atoms.gt.0.and.priv%n_atoms.eq.size(cv%i))

   write (unit = lun, fmt = '(a,a)', advance = 'NO') NFE_INFO, 'atoms = ('

   c = 1
   do a = 1, priv%n_atoms

      nfe_assert(cv%i(a) > 0 .and. cv%i(a) <= pmemd_natoms())
      aname = pmemd_atom_name(cv%i(a))

      write (unit = lun, fmt = '('//pfmt(cv%i(a))//',a,a,a)', advance = 'NO') &
         cv%i(a), ' [', trim(aname), ']'

      if (a == priv%n_atoms) then
         write (unit = lun, fmt = '(a)') ')'
      else if (mod(c, 5) == 0) then
         write (unit = lun, fmt = '(a,/a,3x)', advance = 'NO') ',', NFE_INFO
      else
         write (unit = lun, fmt = '(a)', advance = 'NO') ', '
      end if

      c = c + 1
   end do

   if (cv_nr.eq.0) &
      write (unit = lun, fmt = '(a,a,a)') &
         NFE_INFO, 'reference coordinates are loaded from : ', trim(refcrd_file)

   write (unit = lun, fmt = '(a,a)') &
         NFE_INFO, 'reference coordinates :'

   c = 1
   do a = 1, priv%n_atoms
      nfe_assert(c + 2 <= size(cv%r))
      write (unit = lun, fmt = '(a,3x,i5,a,f8.3,a,f8.3,a,f8.3)') &
         NFE_INFO, cv%i(a), ' : ', &
         cv%r(c), ', ', cv%r(c + 1), &
         ', ', cv%r(c + 2)
      c = c + 3
   end do

end subroutine p_TILT
!=========================================================================
!========================== SPINANGLE ====================================
!=========================================================================

!
! input:
!
! cv%i = (a1, a2, ..., aN)
!
!     (list of participating atoms)
!
! cv%r = (a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z, ...)
!
!        (reference coordinates)
!
! axis = (e1, e2, e3)
!      (axis which defined by user)
! value = spinAngle[-180,180]: phi = 2.0 * atan (alpha)
! where alpha= q.e/q0
! in which q=(q1,q2,q3) is orientation vector and e=axis
!

!=============================================================================

function v_SPINANGLE(cv, x) result(value)

   NFE_USE_AFAILED

   use nfe_lib_mod
   use gbl_constants_mod
   use parallel_dat_mod
   !use nfe_colvar_math, only : norm3

   implicit none

   double precision :: value

   type(colvar_t), intent(inout) :: cv

   double precision, intent(in) :: x(*)
   double precision :: e(3)
   double precision :: get_vec(3)
   double precision :: alpha, beta, phi


   type(priv_t), pointer :: priv

#ifdef MPI
   integer :: error
#endif /* MPI */

   double precision :: cm(3), cm_nrm
   integer :: i, n
! gfortran 4.2 and lower does not support the volatile keyword, but this is
! necessary to work around an intel compiler optimization bug. So turn off
! the volatile keyword for qualifying compilers
#if !defined(__GNUC__) || (__GNUC__ >= 4 && __GNUC_MINOR__ > 2)
   integer, volatile :: a, a3 !work around ifort optimization bug      !VOLATILE attribute and statement specifies
                                                                        !that the value of an object is unpredictable.
#else
   integer :: a, a3
#endif /* __GFORTRAN__ */

   nfe_assert(cv%type == COLVAR_SPINANGLE)

   priv => get_priv(cv)

   nfe_assert(associated(priv))
   nfe_assert(associated(priv%mass))
   nfe_assert(size(priv%mass).gt.0)

   priv%value = ZERO

   cm = ZERO
   n = 1
   do a = 1, priv%n_atoms
      a3 = 3*cv%i(a)
      cm = cm + priv%mass(n)*x(a3 - 2:a3)
      n = n + 1
   end do

   cm = cm/priv%total_mass

   ! populate cm_crd && compute the "norm"
   cm_nrm = ZERO
   n = 1
   do a = 1, priv%n_atoms
      a3 = 3*(n - 1)
      do i = 1, 3
         priv%cm_crd(a3 + i) = x(3*(cv%i(a) - 1) + i) - cm(i)
         cm_nrm = cm_nrm + priv%mass(n)*priv%cm_crd(a3 + i)**2
      end do
      n = n + 1
   end do

   call orientation_q(priv%n_atoms, priv%mass, &
         priv%cm_crd, priv%ref_crd, priv%lambda, priv%quaternion)

   get_vec = [priv%quaternion(2,1), priv%quaternion(3,1), priv%quaternion(4,1)]
   if (norm3(cv%axis).ne.ZERO) then
       e = cv%axis/norm3(cv%axis)
   else
       e = [0.0, 0.0, 1.0] !defualt value for the axis
   end if
   alpha = (dot_product(get_vec, e)) / priv%quaternion(1,1)
   priv%value = 2.0 * atan (alpha) * 180.0/PI
   do while (priv%value > 180.0)    !check the value is between [-180,180] !XX!TO DO! FIX ME!
         priv%value = 360.0 - priv%value
   enddo
   do while (priv%value < - 180.0)
         priv%value = 360.0 + priv%value
   enddo

   value = priv%value

end function v_SPINANGLE

!=============================================================================

!
! assumes that atom positions have not been
!   changed since last call to 'value()'
!

subroutine f_SPINANGLE(cv, fcv, f)

   use nfe_lib_mod
   use gbl_constants_mod
   use parallel_dat_mod
   implicit none

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: fcv
   double precision, intent(inout) :: f(*)
   double precision :: e(3)



   integer :: n, a, a3, f3
   double precision :: dSx(4,4), dSy(4,4), dSz(4,4), QL(3)
   double precision :: dSxq1(4), dSyq1(4), dSzq1(4)
   double precision :: a1x, a1y, a1z, tmp
   double precision :: get_vec(3)
   double precision :: alpha, beta, phi
   type(priv_t), pointer :: priv

#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_SPINANGLE)
   nfe_assert(associated(cv%i))

   priv => get_priv(cv)

   nfe_assert(associated(priv))
   nfe_assert(associated(priv%mass))
   nfe_assert(size(priv%mass).gt.0)

#ifdef MPI
   call mpi_bcast(priv%value, 1, MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

   get_vec = [priv%quaternion(2,1), priv%quaternion(3,1), priv%quaternion(4,1)]
   if (norm3(cv%axis).ne.ZERO) then
       e = cv%axis/norm3(cv%axis)
   else
       e = [0.0, 0.0, 1.0]
   end if
   alpha = (dot_product(get_vec, e)) / priv%quaternion(1,1)                           
   beta = dot_product(get_vec, e)
   if ( priv%quaternion(1,1).ne.ZERO ) then                                           
      tmp = (2.0/priv%quaternion(1,1))*(1.0/(1.0 + alpha*alpha)) * 180.0/PI           
                                                                                     
      QL(1) = (fcv*tmp*(alpha*priv%quaternion(1,2) - &
                 (e(1) * priv%quaternion(2,2) - &
                  e(2) * priv%quaternion(3,2) - &
                  e(3) * priv%quaternion(4,2))))/ &
                 max(priv%lambda(1)-priv%lambda(2),dble(0.00000001))

      QL(2) = (fcv*tmp*(alpha*priv%quaternion(1,3) - &
                 (e(1) * priv%quaternion(2,3) - &
                  e(2) * priv%quaternion(3,3) - &
                  e(3) * priv%quaternion(4,3))))/ &
                 max(priv%lambda(1)-priv%lambda(3),dble(0.00000001))

      QL(3) = (fcv*tmp*(alpha*priv%quaternion(1,4) - &
                 (e(1) * priv%quaternion(2,4) - &
                  e(2) * priv%quaternion(3,4) - &
                  e(3) * priv%quaternion(4,4))))/ &
                 max(priv%lambda(1)-priv%lambda(4),dble(0.00000001))
   else
      tmp = 2.0/dot_product(get_vec, e)                                   
                                                                 
      QL(1) = fcv*tmp*priv%quaternion(1,2) / &
              max(priv%lambda(1)-priv%lambda(2),dble(0.00000001))
      QL(2) = fcv*tmp*priv%quaternion(1,3) / &
              max(priv%lambda(1)-priv%lambda(3),dble(0.00000001))
      QL(3) = fcv*tmp*priv%quaternion(1,4) / &
              max(priv%lambda(1)-priv%lambda(4),dble(0.00000001))  
   end if

      n = 1
      do a = 1, priv%n_atoms
         a3 = 3*n
         f3 = 3*cv%i(a)

         a1x = priv%ref_crd(a3 - 2)
         a1y = priv%ref_crd(a3 - 1)
         a1z = priv%ref_crd(a3)

         dSx(1,1) = a1x
         dSy(1,1) = a1y
         dSz(1,1) = a1z

         dSx(2,1) = 0.0
         dSy(2,1) = -a1z
         dSz(2,1) = a1y

         dSx(1,2) = 0.0
         dSy(1,2) = -a1z
         dSz(1,2) = a1y

         dSx(3,1) = a1z
         dSy(3,1) = 0.0
         dSz(3,1) = -a1x

         dSx(1,3) = a1z
         dSy(1,3) = 0.0
         dSz(1,3) = -a1x

         dSx(4,1) = -a1y
         dSy(4,1) = a1x
         dSz(4,1) = 0.0

         dSx(1,4) = -a1y
         dSy(1,4) = a1x
         dSz(1,4) = 0.0

         dSx(2,2) = a1x
         dSy(2,2) = -a1y
         dSz(2,2) = -a1z

         dSx(3,2) = a1y
         dSy(3,2) = a1x
         dSz(3,2) = 0.0

         dSx(2,3) = a1y
         dSy(2,3) = a1x
         dSz(2,3) = 0.0

         dSx(4,2) = a1z
         dSy(4,2) = 0.0
         dSz(4,2) = a1x

         dSx(2,4) = a1z
         dSy(2,4) = 0.0
         dSz(2,4) = a1x

         dSx(3,3) = -a1x
         dSy(3,3) = a1y
         dSz(3,3) = -a1z

         dSx(4,3) = 0.0
         dSy(4,3) = a1z
         dSz(4,3) = a1y

         dSx(3,4) = 0.0
         dSy(3,4) = a1z
         dSz(3,4) = a1y

         dSx(4,4) = -a1x
         dSy(4,4) = -a1y
         dSz(4,4) = a1z

         dSxq1 = matmul(dSx,priv%quaternion(:,1))
         dSyq1 = matmul(dSy,priv%quaternion(:,1))
         dSzq1 = matmul(dSz,priv%quaternion(:,1))

          f(f3 - 2) = f(f3 - 2) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSxq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSxq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSxq1,priv%quaternion(:,4)) )

          f(f3 - 1) = f(f3 - 1) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSyq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSyq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSyq1,priv%quaternion(:,4)) )

          f(f3) = f(f3) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSzq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSzq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSzq1,priv%quaternion(:,4)) )
         n = n + 1
      end do

end subroutine f_SPINANGLE

!=============================================================================

subroutine b_SPINANGLE(cv, cvno, amass)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)
   double precision, pointer :: coor(:) => null()

   integer :: i, n_atoms, error
   integer :: a, b, a1
   double precision :: cm(3)

   type(priv_t), pointer :: priv


   nfe_assert(cv%type == COLVAR_SPINANGLE)

   ! very basic checks
   if (.not. associated(cv%i)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : no integers found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cv%i)

   if (cv_nr.gt.0) then 
    if (.not. associated(cv%r)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_INFO, 'CV #', cvno, &
            ' (ORIENTATION) : no reals (reference coordinates) found'
      NFE_MASTER_ONLY_END
      call terminate()
    end if
   end if 

   ! count the number of atoms
   n_atoms = 0

   do i = 1, size(cv%i)
      if (cv%i(i) == 0) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(i)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : unexpected zero (integer #', i, ')'
            NFE_MASTER_ONLY_END
            call terminate()
      else if (cv%i(i) .lt. 0) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : negative integer'
         NFE_MASTER_ONLY_END
         call terminate()
      else
         n_atoms = n_atoms + 1
      end if
   end do

   if (cv_nr.eq.0) then
      if (associated(coor)) &
          deallocate(coor)
          allocate(coor(3*pmemd_natoms()), stat = error)
          if (error.ne.0) &
           NFE_OUT_OF_MEMORY
           call colvar_read_refcrd(coor, refcrd_file)
           allocate(cv%r(3*size(cv%i)), stat = error)
             if (error.ne.0) &
              NFE_OUT_OF_MEMORY
              do a = 1, n_atoms
               a1=cv%i(a)
               cv%r(1 + 3*(a - 1): 3*a) = coor(1 + 3*(a1 - 1): 3*a1)
              end do
          deallocate(coor)
   end if

   if (size(cv%r) /= 3*n_atoms) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : wrong number of reals'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! size(cv%r) /= 3*n_atoms

   ! allocate priv_t instance for this variable
   priv => new_priv(cv)

   ! basic checks
   if (n_atoms < 3) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, &
            fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
            (1)//',a,'//pfmt(n_atoms)//',a/)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (ORIENTATION) : too few integers in reference (', 1, ':', n_atoms, ')'
      NFE_MASTER_ONLY_END
      call terminate()
   end if

   do a = 1, n_atoms
      if (cv%i(a) < 1 .or. cv%i(a) > pmemd_natoms()) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
               (cv%i(a))//',a,'//pfmt(pmemd_natoms())//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : integer #', a, ' (', cv%i(a), &
               ') is out of range [1, ', pmemd_natoms(), ']'
         NFE_MASTER_ONLY_END
         call terminate()
      end if
   end do

   ! check for duplicates
   do a = 1, n_atoms
      do b = a + 1, n_atoms
         if (cv%i(a) == cv%i(b)) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
                  (b)//',a,'//pfmt(cv%i(a))//',a/)') &
                  NFE_ERROR, 'CV #', cvno, &
                  ' (ORIENTATION) : integers #', a, ' and #', b, &
                  ' are equal (', cv%i(a), ')'
            NFE_MASTER_ONLY_END
            call terminate()
         end if ! cv%i(a) == cv%i(b)
      end do
   end do

   ! allocate/setup

   priv%n_atoms = n_atoms

   allocate(priv%mass(n_atoms), priv%cm_crd(3*n_atoms), &
            priv%ref_crd(3*n_atoms), stat = error)

   if (error /= 0) &
      NFE_OUT_OF_MEMORY

   cm = ZERO
   priv%total_mass = ZERO

   do a = 1, n_atoms
      priv%mass(a) = amass(cv%i(a))
      priv%total_mass = priv%total_mass + priv%mass(a)
      cm = cm + priv%mass(a)*cv%r(1 + 3*(a - 1): 3*a)
   end do

   cm = cm/priv%total_mass

   ! translate reference coordinates to CM frame
   priv%ref_nrm = ZERO
   do a = 1, n_atoms
      do b = 1, 3
         priv%ref_crd(3*(a - 1) + b) = cv%r(3*(a - 1) + b) - cm(b)
         priv%ref_nrm = priv%ref_nrm + priv%mass(a)*priv%ref_crd(3*(a - 1) + b)**2
      end do
   end do


end subroutine b_SPINANGLE

!=============================================================================

subroutine c_SPINANGLE(cv)

   NFE_USE_AFAILED

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(inout) :: cv

   type(priv_t), pointer :: priv

   nfe_assert(cv%type == COLVAR_SPINANGLE)
   nfe_assert(associated(cv%i).and.size(cv%i).gt.0)

   priv => get_priv(cv)

   nfe_assert(associated(priv%mass))
   nfe_assert(associated(priv%cm_crd))
   nfe_assert(associated(priv%ref_crd))

   deallocate(priv%mass, priv%cm_crd, priv%ref_crd)

   call del_priv(cv)

end subroutine c_SPINANGLE

!=============================================================================

subroutine p_SPINANGLE(cv, lun)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   type(priv_t), pointer :: priv
   integer :: a, c
   character(4) :: aname

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_SPINANGLE)
   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   priv => get_priv(cv)
   nfe_assert(priv%n_atoms.gt.0.and.priv%n_atoms.eq.size(cv%i))

   write (unit = lun, fmt = '(a,a)', advance = 'NO') NFE_INFO, 'atoms = ('

   c = 1
   do a = 1, priv%n_atoms

      nfe_assert(cv%i(a) > 0 .and. cv%i(a) <= pmemd_natoms())
      aname = pmemd_atom_name(cv%i(a))

      write (unit = lun, fmt = '('//pfmt(cv%i(a))//',a,a,a)', advance = 'NO') &
         cv%i(a), ' [', trim(aname), ']'

      if (a == priv%n_atoms) then
         write (unit = lun, fmt = '(a)') ')'
      else if (mod(c, 5) == 0) then
         write (unit = lun, fmt = '(a,/a,3x)', advance = 'NO') ',', NFE_INFO
      else
         write (unit = lun, fmt = '(a)', advance = 'NO') ', '
      end if

      c = c + 1
   end do

   if (cv_nr.eq.0) &
      write (unit = lun, fmt = '(a,a,a)') &
         NFE_INFO, 'reference coordinates are loaded from : ', trim(refcrd_file)

   write (unit = lun, fmt = '(a,a)') &
         NFE_INFO, 'reference coordinates :'

   c = 1
   do a = 1, priv%n_atoms
      nfe_assert(c + 2 <= size(cv%r))
      write (unit = lun, fmt = '(a,3x,i5,a,f8.3,a,f8.3,a,f8.3)') &
         NFE_INFO, cv%i(a), ' : ', &
         cv%r(c), ', ', cv%r(c + 1), &
         ', ', cv%r(c + 2)
      c = c + 3
   end do

end subroutine p_SPINANGLE
!=========================================================================
!========================== QUATERNION0 ==================================
!=========================================================================

!
! input:
!
! cv%i = (a1, a2, ..., aN)
!
!     (list of participating atoms)
!
! cv%r = (a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z, ...)
!
!        (reference coordinates)
!
! value = q0
! in which (q0,q1,q2,q3) is orientation quaternion representing optimum rotation w.r.t. reference
!
! add normalization
!=============================================================================

function v_QUATERNION0(cv, x) result(value)

   NFE_USE_AFAILED

   use nfe_lib_mod
   use gbl_constants_mod
   use parallel_dat_mod

   implicit none

   double precision :: value

   type(colvar_t), intent(inout) :: cv

   double precision, intent(in) :: x(*)



   type(priv_t), pointer :: priv

#ifdef MPI
   integer :: error
#endif /* MPI */

   double precision :: cm(3), cm_nrm
   integer :: i, n
! gfortran 4.2 and lower does not support the volatile keyword, but this is
! necessary to work around an intel compiler optimization bug. So turn off
! the volatile keyword for qualifying compilers
#if !defined(__GNUC__) || (__GNUC__ >= 4 && __GNUC_MINOR__ > 2)
   integer, volatile :: a, a3 ! work around ifort optimization bug
#else
   integer :: a, a3
#endif /* __GFORTRAN__ */

   nfe_assert(cv%type == COLVAR_QUATERNION0)

   priv => get_priv(cv)

   nfe_assert(associated(priv))
   nfe_assert(associated(priv%mass))
   nfe_assert(size(priv%mass).gt.0)

   priv%value = ZERO

   ! compute the center of mass (of the moving atoms)
   cm = ZERO
   n = 1
   do a = 1, priv%n_atoms
      a3 = 3*cv%i(a)
      cm = cm + priv%mass(n)*x(a3 - 2:a3)
      n = n + 1
   end do

   cm = cm/priv%total_mass

   ! populate cm_crd && compute the "norm"
   cm_nrm = ZERO
   n = 1
   do a = 1, priv%n_atoms
      a3 = 3*(n - 1)
      do i = 1, 3
         priv%cm_crd(a3 + i) = x(3*(cv%i(a) - 1) + i) - cm(i)
         cm_nrm = cm_nrm + priv%mass(n)*priv%cm_crd(a3 + i)**2
      end do
      n = n + 1
   end do

   call orientation_q(priv%n_atoms, priv%mass, &
         priv%cm_crd, priv%ref_crd, priv%lambda, priv%quaternion)
!#ifdef MPI
   if (priv%quaternion(1,1) >= 0.0) then
     priv%value = priv%quaternion(1,1)
   else
     priv%value = -1.0 * priv%quaternion(1,1)
   endif
!#endif /* MPI */

   value = priv%value

end function v_QUATERNION0

!=============================================================================

!
! assumes that atom positions have not been
!   changed since last call to 'value()'
!

subroutine f_QUATERNION0(cv, fcv, f)

   use nfe_lib_mod
   use gbl_constants_mod
   use parallel_dat_mod
   implicit none

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: fcv
   double precision, intent(inout) :: f(*)



   integer :: n, a, a3, f3
   double precision :: dSx(4,4), dSy(4,4), dSz(4,4), QL(3)
   double precision :: dSxq1(4), dSyq1(4), dSzq1(4)
   double precision :: a1x, a1y, a1z, tmp
   type(priv_t), pointer :: priv

#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_QUATERNION0)
   nfe_assert(associated(cv%i))

   priv => get_priv(cv)

   nfe_assert(associated(priv))
   nfe_assert(associated(priv%mass))
   nfe_assert(size(priv%mass).gt.0)

#ifdef MPI
   call mpi_bcast(priv%value, 1, MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

      QL(1) = fcv * priv%quaternion(1,2) / &
                    max(priv%lambda(1)-priv%lambda(2),dble(0.00000001))
      QL(2) = fcv * priv%quaternion(1,3) / &
                    max(priv%lambda(1)-priv%lambda(3),dble(0.00000001))
      QL(3) = fcv * priv%quaternion(1,4) / &
                    max(priv%lambda(1)-priv%lambda(4),dble(0.00000001))

      n = 1
      do a = 1, priv%n_atoms
         a3 = 3*n
         f3 = 3*cv%i(a)

         a1x = priv%ref_crd(a3 - 2)
         a1y = priv%ref_crd(a3 - 1)
         a1z = priv%ref_crd(a3)

         dSx(1,1) = a1x
         dSy(1,1) = a1y
         dSz(1,1) = a1z

         dSx(2,1) = 0.0
         dSy(2,1) = -a1z
         dSz(2,1) = a1y

         dSx(1,2) = 0.0
         dSy(1,2) = -a1z
         dSz(1,2) = a1y

         dSx(3,1) = a1z
         dSy(3,1) = 0.0
         dSz(3,1) = -a1x

         dSx(1,3) = a1z
         dSy(1,3) = 0.0
         dSz(1,3) = -a1x

         dSx(4,1) = -a1y
         dSy(4,1) = a1x
         dSz(4,1) = 0.0

         dSx(1,4) = -a1y
         dSy(1,4) = a1x
         dSz(1,4) = 0.0

         dSx(2,2) = a1x
         dSy(2,2) = -a1y
         dSz(2,2) = -a1z

         dSx(3,2) = a1y
         dSy(3,2) = a1x
         dSz(3,2) = 0.0

         dSx(2,3) = a1y
         dSy(2,3) = a1x
         dSz(2,3) = 0.0

         dSx(4,2) = a1z
         dSy(4,2) = 0.0
         dSz(4,2) = a1x

         dSx(2,4) = a1z
         dSy(2,4) = 0.0
         dSz(2,4) = a1x

         dSx(3,3) = -a1x
         dSy(3,3) = a1y
         dSz(3,3) = -a1z

         dSx(4,3) = 0.0
         dSy(4,3) = a1z
         dSz(4,3) = a1y

         dSx(3,4) = 0.0
         dSy(3,4) = a1z
         dSz(3,4) = a1y

         dSx(4,4) = -a1x
         dSy(4,4) = -a1y
         dSz(4,4) = a1z

         dSxq1 = matmul(dSx,priv%quaternion(:,1))
         dSyq1 = matmul(dSy,priv%quaternion(:,1))
         dSzq1 = matmul(dSz,priv%quaternion(:,1))

         f(f3 - 2) = f(f3 - 2) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSxq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSxq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSxq1,priv%quaternion(:,4)) )

         f(f3 - 1) = f(f3 - 1) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSyq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSyq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSyq1,priv%quaternion(:,4)) )

         f(f3) = f(f3) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSzq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSzq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSzq1,priv%quaternion(:,4)) )

         n = n + 1
      end do

end subroutine f_QUATERNION0

!=============================================================================
subroutine b_QUATERNION0(cv, cvno, amass)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)
   double precision, pointer :: coor(:) => null()

   integer :: i, n_atoms, error
   integer :: a, b, a1 
   double precision :: cm(3)

   type(priv_t), pointer :: priv


   nfe_assert(cv%type == COLVAR_QUATERNION0)

   ! very basic checks
   if (.not. associated(cv%i)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : no integers found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cv%i)

   if (cv_nr.gt.0) then
    if (.not. associated(cv%r)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_INFO, 'CV #', cvno, &
            ' (ORIENTATION) : no reals (reference coordinates) found'
      NFE_MASTER_ONLY_END
      call terminate()
    end if
   end if 

   ! count the number of atoms
   n_atoms = 0

   do i = 1, size(cv%i)
      if (cv%i(i) == 0) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(i)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : unexpected zero (integer #', i, ')'
            NFE_MASTER_ONLY_END
            call terminate()
      else if (cv%i(i) .lt. 0) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : negative integer'
         NFE_MASTER_ONLY_END
         call terminate()
      else
         n_atoms = n_atoms + 1
      end if
   end do

   if (cv_nr.eq.0) then
      if (associated(coor)) &
          deallocate(coor)
          allocate(coor(3*pmemd_natoms()), stat = error)
          if (error.ne.0) &
           NFE_OUT_OF_MEMORY
           call colvar_read_refcrd(coor, refcrd_file)
           allocate(cv%r(3*size(cv%i)), stat = error)
             if (error.ne.0) &
              NFE_OUT_OF_MEMORY
              do a = 1, n_atoms
               a1=cv%i(a)
               cv%r(1 + 3*(a - 1): 3*a) = coor(1 + 3*(a1 - 1): 3*a1)
              end do
          deallocate(coor)
   end if

   if (size(cv%r) /= 3*n_atoms) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : wrong number of reals'
      NFE_MASTER_ONLY_END
      call terminate()
   end if 

   ! allocate priv_t instance for this variable
   priv => new_priv(cv)

   ! basic checks
   if (n_atoms < 3) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, &
            fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
            (1)//',a,'//pfmt(n_atoms)//',a/)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (ORIENTATION) : too few integers in reference (', 1, ':', n_atoms, ')'
      NFE_MASTER_ONLY_END
      call terminate()
   end if

   do a = 1, n_atoms
      if (cv%i(a) < 1 .or. cv%i(a) > pmemd_natoms()) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
               (cv%i(a))//',a,'//pfmt(pmemd_natoms())//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : integer #', a, ' (', cv%i(a), &
               ') is out of range [1, ', pmemd_natoms(), ']'
         NFE_MASTER_ONLY_END
         call terminate()
      end if
   end do

   ! check for duplicates
   do a = 1, n_atoms
      do b = a + 1, n_atoms
         if (cv%i(a) == cv%i(b)) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
                  (b)//',a,'//pfmt(cv%i(a))//',a/)') &
                  NFE_ERROR, 'CV #', cvno, &
                  ' (ORIENTATION) : integers #', a, ' and #', b, &
                  ' are equal (', cv%i(a), ')'
            NFE_MASTER_ONLY_END
            call terminate()
         end if
      end do
   end do

   ! allocate/setup

   priv%n_atoms = n_atoms

   allocate(priv%mass(n_atoms), priv%cm_crd(3*n_atoms), &
            priv%ref_crd(3*n_atoms), stat = error)

   if (error /= 0) &
      NFE_OUT_OF_MEMORY

   cm = ZERO
   priv%total_mass = ZERO

   do a = 1, n_atoms
      priv%mass(a) = amass(cv%i(a))
      priv%total_mass = priv%total_mass + priv%mass(a)
      cm = cm + priv%mass(a)*cv%r(1 + 3*(a - 1): 3*a)
   end do

   cm = cm/priv%total_mass

   ! translate reference coordinates to CM frame
   priv%ref_nrm = ZERO
   do a = 1, n_atoms
      do b = 1, 3
         priv%ref_crd(3*(a - 1) + b) = cv%r(3*(a - 1) + b) - cm(b)
         priv%ref_nrm = priv%ref_nrm + priv%mass(a)*priv%ref_crd(3*(a - 1) + b)**2
      end do
   end do


end subroutine b_QUATERNION0

!=============================================================================

subroutine c_QUATERNION0(cv)

   NFE_USE_AFAILED

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(inout) :: cv

   type(priv_t), pointer :: priv

   nfe_assert(cv%type == COLVAR_QUATERNION0)
   nfe_assert(associated(cv%i).and.size(cv%i).gt.0)

   priv => get_priv(cv)

   nfe_assert(associated(priv%mass))
   nfe_assert(associated(priv%cm_crd))
   nfe_assert(associated(priv%ref_crd))

   deallocate(priv%mass, priv%cm_crd, priv%ref_crd)

   call del_priv(cv)

end subroutine c_QUATERNION0

!=============================================================================

subroutine p_QUATERNION0(cv, lun)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   type(priv_t), pointer :: priv
   integer :: a, c
   character(4) :: aname

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_QUATERNION0)
   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   priv => get_priv(cv)
   nfe_assert(priv%n_atoms.gt.0.and.priv%n_atoms.eq.size(cv%i))

   write (unit = lun, fmt = '(a,a)', advance = 'NO') NFE_INFO, 'atoms = ('

   c = 1
   do a = 1, priv%n_atoms

      nfe_assert(cv%i(a) > 0 .and. cv%i(a) <= pmemd_natoms())
      aname = pmemd_atom_name(cv%i(a))

      write (unit = lun, fmt = '('//pfmt(cv%i(a))//',a,a,a)', advance = 'NO') &
         cv%i(a), ' [', trim(aname), ']'

      if (a == priv%n_atoms) then
         write (unit = lun, fmt = '(a)') ')'
      else if (mod(c, 5) == 0) then
         write (unit = lun, fmt = '(a,/a,3x)', advance = 'NO') ',', NFE_INFO
      else
         write (unit = lun, fmt = '(a)', advance = 'NO') ', '
      end if

      c = c + 1
   end do

   if (cv_nr.eq.0) &
      write (unit = lun, fmt = '(a,a,a)') &
         NFE_INFO, 'reference coordinates are loaded from : ', trim(refcrd_file)

   write (unit = lun, fmt = '(a,a)') &
         NFE_INFO, 'reference coordinates :'

   c = 1
   do a = 1, priv%n_atoms
      nfe_assert(c + 2 <= size(cv%r))
      write (unit = lun, fmt = '(a,3x,i5,a,f8.3,a,f8.3,a,f8.3)') &
         NFE_INFO, cv%i(a), ' : ', &
         cv%r(c), ', ', cv%r(c + 1), &
         ', ', cv%r(c + 2)
      c = c + 3
   end do

end subroutine p_QUATERNION0
!=========================================================================
!========================== QUATERNION1 ==================================
!=========================================================================

!
! input:
!
! cv%i = (a1, a2, ..., aN)
!
!     (list of participating atoms)
!
! cv%r = (a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z, ...)
!
!        (reference coordinates)
!
! value = q1
! in which (q0,q1,q2,q3) is orientation quaternion representing optimum rotation w.r.t. reference
!
!=============================================================================

function v_QUATERNION1(cv, x) result(value)

   NFE_USE_AFAILED

   use nfe_lib_mod
   use parallel_dat_mod
   use gbl_constants_mod

   implicit none

   double precision :: value

   type(colvar_t), intent(inout) :: cv

   double precision, intent(in) :: x(*)

!#  include "nfe-mpi.h"

   type(priv_t), pointer :: priv

#ifdef MPI
   integer :: error
#endif /* MPI */

   double precision :: cm(3), cm_nrm
   integer :: i, n
! gfortran 4.2 and lower does not support the volatile keyword, but this is
! necessary to work around an intel compiler optimization bug. So turn off
! the volatile keyword for qualifying compilers
#if !defined(__GNUC__) || (__GNUC__ >= 4 && __GNUC_MINOR__ > 2)
   integer, volatile :: a, a3 ! work around ifort optimization bug
#else
   integer :: a, a3
#endif /* __GFORTRAN__ */

   nfe_assert(cv%type == COLVAR_QUATERNION1)

   priv => get_priv(cv)

   nfe_assert(associated(priv))
   nfe_assert(associated(priv%mass))
   nfe_assert(size(priv%mass).gt.0)

   priv%value = ZERO

   cm = ZERO
   n = 1
   do a = 1, priv%n_atoms
      a3 = 3*cv%i(a)
      cm = cm + priv%mass(n)*x(a3 - 2:a3)
      n = n + 1
   end do

   cm = cm/priv%total_mass

   ! populate cm_crd && compute the "norm"
   cm_nrm = ZERO
   n = 1
   do a = 1, priv%n_atoms
      a3 = 3*(n - 1)
      do i = 1, 3
         priv%cm_crd(a3 + i) = x(3*(cv%i(a) - 1) + i) - cm(i)
         cm_nrm = cm_nrm + priv%mass(n)*priv%cm_crd(a3 + i)**2
      end do
      n = n + 1
   end do

   call orientation_q(priv%n_atoms, priv%mass, &
         priv%cm_crd, priv%ref_crd, priv%lambda, priv%quaternion)
!#ifdef MPI
   if (priv%quaternion(1,1) >= 0.0) then
     priv%value = priv%quaternion(2,1)
   else
     priv%value = -1.0 * priv%quaternion(2,1)
   endif
!#endif /* MPI */

   value = priv%value

end function v_QUATERNION1

!=============================================================================

!
! assumes that atom positions have not been
!   changed since last call to 'value()'
!

subroutine f_QUATERNION1(cv, fcv, f)

   use nfe_lib_mod
   use gbl_constants_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: fcv
   double precision, intent(inout) :: f(*)

!#  include "nfe-mpi.h"

   integer :: n, a, a3, f3
   double precision :: dSx(4,4), dSy(4,4), dSz(4,4), QL(3)
   double precision :: dSxq1(4), dSyq1(4), dSzq1(4)
   double precision :: a1x, a1y, a1z, tmp
   type(priv_t), pointer :: priv

#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_QUATERNION1)
   nfe_assert(associated(cv%i))

   priv => get_priv(cv)

   nfe_assert(associated(priv))
   nfe_assert(associated(priv%mass))
   nfe_assert(size(priv%mass).gt.0)

#ifdef MPI
   call mpi_bcast(priv%value, 1, MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

      QL(1) = fcv * priv%quaternion(2,2) / &
                    max(priv%lambda(1)-priv%lambda(2),dble(0.00000001))
      QL(2) = fcv * priv%quaternion(2,3) / &
                    max(priv%lambda(1)-priv%lambda(3),dble(0.00000001))
      QL(3) = fcv * priv%quaternion(2,4) / &
                    max(priv%lambda(1)-priv%lambda(4),dble(0.00000001))
      n = 1
      do a = 1, priv%n_atoms
         a3 = 3*n
         f3 = 3*cv%i(a)

         a1x = priv%ref_crd(a3 - 2)
         a1y = priv%ref_crd(a3 - 1)
         a1z = priv%ref_crd(a3)

         dSx(1,1) = a1x
         dSy(1,1) = a1y
         dSz(1,1) = a1z

         dSx(2,1) = 0.0
         dSy(2,1) = -a1z
         dSz(2,1) = a1y

         dSx(1,2) = 0.0
         dSy(1,2) = -a1z
         dSz(1,2) = a1y

         dSx(3,1) = a1z
         dSy(3,1) = 0.0
         dSz(3,1) = -a1x

         dSx(1,3) = a1z
         dSy(1,3) = 0.0
         dSz(1,3) = -a1x

         dSx(4,1) = -a1y
         dSy(4,1) = a1x
         dSz(4,1) = 0.0

         dSx(1,4) = -a1y
         dSy(1,4) = a1x
         dSz(1,4) = 0.0

         dSx(2,2) = a1x
         dSy(2,2) = -a1y
         dSz(2,2) = -a1z

         dSx(3,2) = a1y
         dSy(3,2) = a1x
         dSz(3,2) = 0.0

         dSx(2,3) = a1y
         dSy(2,3) = a1x
         dSz(2,3) = 0.0

         dSx(4,2) = a1z
         dSy(4,2) = 0.0
         dSz(4,2) = a1x

         dSx(2,4) = a1z
         dSy(2,4) = 0.0
         dSz(2,4) = a1x

         dSx(3,3) = -a1x
         dSy(3,3) = a1y
         dSz(3,3) = -a1z

         dSx(4,3) = 0.0
         dSy(4,3) = a1z
         dSz(4,3) = a1y

         dSx(3,4) = 0.0
         dSy(3,4) = a1z
         dSz(3,4) = a1y

         dSx(4,4) = -a1x
         dSy(4,4) = -a1y
         dSz(4,4) = a1z

         dSxq1 = matmul(dSx,priv%quaternion(:,1))
         dSyq1 = matmul(dSy,priv%quaternion(:,1))
         dSzq1 = matmul(dSz,priv%quaternion(:,1))

         f(f3 - 2) = f(f3 - 2) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSxq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSxq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSxq1,priv%quaternion(:,4)) )

         f(f3 - 1) = f(f3 - 1) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSyq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSyq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSyq1,priv%quaternion(:,4)) )

         f(f3) = f(f3) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSzq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSzq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSzq1,priv%quaternion(:,4)) )

         n = n + 1
      end do

end subroutine f_QUATERNION1

!=============================================================================
subroutine b_QUATERNION1(cv, cvno, amass)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)
   double precision, pointer :: coor(:) => null()

   integer :: i, n_atoms, error
   integer :: a, b, a1
   double precision :: cm(3)

   type(priv_t), pointer :: priv


   nfe_assert(cv%type == COLVAR_QUATERNION1)

   ! very basic checks
   if (.not. associated(cv%i)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : no integers found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cv%i)

   if (cv_nr.gt.0) then
    if (.not. associated(cv%r)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_INFO, 'CV #', cvno, &
            ' (ORIENTATION) : no reals (reference coordinates) found'
      NFE_MASTER_ONLY_END
      call terminate()
    end if
   end if 

   ! count the number of atoms
   n_atoms = 0

   do i = 1, size(cv%i)
      if (cv%i(i) == 0) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(i)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : unexpected zero (integer #', i, ')'
            NFE_MASTER_ONLY_END
            call terminate()
      else if (cv%i(i) .lt. 0) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : negative integer'
         NFE_MASTER_ONLY_END
         call terminate()
      else
         n_atoms = n_atoms + 1
      end if
   end do

   if (cv_nr.eq.0) then
      if (associated(coor)) &
          deallocate(coor)
          allocate(coor(3*pmemd_natoms()), stat = error)
          if (error.ne.0) &
           NFE_OUT_OF_MEMORY
           call colvar_read_refcrd(coor, refcrd_file)
           allocate(cv%r(3*size(cv%i)), stat = error)
             if (error.ne.0) &
              NFE_OUT_OF_MEMORY
              do a = 1, n_atoms
               a1=cv%i(a)
               cv%r(1 + 3*(a - 1): 3*a) = coor(1 + 3*(a1 - 1): 3*a1)
              end do
          deallocate(coor)
   end if

   if (size(cv%r) /= 3*n_atoms) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : wrong number of reals'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! size(cv%r) /= 3*n_atoms

   ! allocate priv_t instance for this variable
   priv => new_priv(cv)

   ! basic checks
   if (n_atoms < 3) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, &
            fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
            (1)//',a,'//pfmt(n_atoms)//',a/)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (ORIENTATION) : too few integers in reference (', 1, ':', n_atoms, ')'
      NFE_MASTER_ONLY_END
      call terminate()
   end if

   do a = 1, n_atoms
      if (cv%i(a) < 1 .or. cv%i(a) > pmemd_natoms()) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
               (cv%i(a))//',a,'//pfmt(pmemd_natoms())//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : integer #', a, ' (', cv%i(a), &
               ') is out of range [1, ', pmemd_natoms(), ']'
         NFE_MASTER_ONLY_END
         call terminate()
      end if
   end do

   ! check for duplicates
   do a = 1, n_atoms
      do b = a + 1, n_atoms
         if (cv%i(a) == cv%i(b)) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
                  (b)//',a,'//pfmt(cv%i(a))//',a/)') &
                  NFE_ERROR, 'CV #', cvno, &
                  ' (ORIENTATION) : integers #', a, ' and #', b, &
                  ' are equal (', cv%i(a), ')'
            NFE_MASTER_ONLY_END
            call terminate()
         end if 
      end do
   end do

   ! allocate/setup

   priv%n_atoms = n_atoms

   allocate(priv%mass(n_atoms), priv%cm_crd(3*n_atoms), &
            priv%ref_crd(3*n_atoms), stat = error)

   if (error /= 0) &
      NFE_OUT_OF_MEMORY

   cm = ZERO
   priv%total_mass = ZERO

   do a = 1, n_atoms
      priv%mass(a) = amass(cv%i(a))
      priv%total_mass = priv%total_mass + priv%mass(a)
      cm = cm + priv%mass(a)*cv%r(1 + 3*(a - 1): 3*a)
   end do

   cm = cm/priv%total_mass

   ! translate reference coordinates to CM frame
   priv%ref_nrm = ZERO
   do a = 1, n_atoms
      do b = 1, 3
         priv%ref_crd(3*(a - 1) + b) = cv%r(3*(a - 1) + b) - cm(b)
         priv%ref_nrm = priv%ref_nrm + priv%mass(a)*priv%ref_crd(3*(a - 1) + b)**2
      end do
   end do


end subroutine b_QUATERNION1

!=============================================================================

subroutine c_QUATERNION1(cv)

   NFE_USE_AFAILED

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(inout) :: cv

   type(priv_t), pointer :: priv

   nfe_assert(cv%type == COLVAR_QUATERNION1)
   nfe_assert(associated(cv%i).and.size(cv%i).gt.0)

   priv => get_priv(cv)

   nfe_assert(associated(priv%mass))
   nfe_assert(associated(priv%cm_crd))
   nfe_assert(associated(priv%ref_crd))

   deallocate(priv%mass, priv%cm_crd, priv%ref_crd)

   call del_priv(cv)

end subroutine c_QUATERNION1

!=============================================================================

subroutine p_QUATERNION1(cv, lun)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   type(priv_t), pointer :: priv
   integer :: a, c
   character(4) :: aname

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_QUATERNION1)
   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   priv => get_priv(cv)
   nfe_assert(priv%n_atoms.gt.0.and.priv%n_atoms.eq.size(cv%i))

   write (unit = lun, fmt = '(a,a)', advance = 'NO') NFE_INFO, 'atoms = ('

   c = 1
   do a = 1, priv%n_atoms

      nfe_assert(cv%i(a) > 0 .and. cv%i(a) <= pmemd_natoms())
      aname = pmemd_atom_name(cv%i(a))

      write (unit = lun, fmt = '('//pfmt(cv%i(a))//',a,a,a)', advance = 'NO') &
         cv%i(a), ' [', trim(aname), ']'

      if (a == priv%n_atoms) then
         write (unit = lun, fmt = '(a)') ')'
      else if (mod(c, 5) == 0) then
         write (unit = lun, fmt = '(a,/a,3x)', advance = 'NO') ',', NFE_INFO
      else
         write (unit = lun, fmt = '(a)', advance = 'NO') ', '
      end if

      c = c + 1
   end do

   if (cv_nr.eq.0) &
      write (unit = lun, fmt = '(a,a,a)') &
         NFE_INFO, 'reference coordinates are loaded from : ', trim(refcrd_file)

   write (unit = lun, fmt = '(a,a)') &
         NFE_INFO, 'reference coordinates :'

   c = 1
   do a = 1, priv%n_atoms
      nfe_assert(c + 2 <= size(cv%r))
      write (unit = lun, fmt = '(a,3x,i5,a,f8.3,a,f8.3,a,f8.3)') &
         NFE_INFO, cv%i(a), ' : ', &
         cv%r(c), ', ', cv%r(c + 1), &
         ', ', cv%r(c + 2)
      c = c + 3
   end do

end subroutine p_QUATERNION1
!=========================================================================
!========================== QUATERNION2 ==================================
!=========================================================================

!
! input:
!
! cv%i = (a1, a2, ..., aN)
!
!     (list of participating atoms)
!
! cv%r = (a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z, ...)
!
!        (reference coordinates)
!
! value = q2
! in which (q0,q1,q2,q3) is orientation quaternion representing optimum rotation w.r.t. reference
!

!=============================================================================

function v_QUATERNION2(cv, x) result(value)

   NFE_USE_AFAILED

   use nfe_lib_mod
   use parallel_dat_mod
   use gbl_constants_mod

   implicit none

   double precision :: value

   type(colvar_t), intent(inout) :: cv

   double precision, intent(in) :: x(*)

!#  include "nfe-mpi.h"

   type(priv_t), pointer :: priv

#ifdef MPI
   integer :: error
#endif /* MPI */

   double precision :: cm(3), cm_nrm
   integer :: i, n
! gfortran 4.2 and lower does not support the volatile keyword, but this is
! necessary to work around an intel compiler optimization bug. So turn off
! the volatile keyword for qualifying compilers
!#if !defined(__GNUC__) || (__GNUC__ >= 4 && __GNUC_MINOR__ > 2)
!   integer, volatile :: a, a3 ! work around ifort optimization bug
!#else
   integer :: a, a3
!#endif /* __GFORTRAN__ */

   nfe_assert(cv%type == COLVAR_QUATERNION2)

   priv => get_priv(cv)

   nfe_assert(associated(priv))
   nfe_assert(associated(priv%mass))
   nfe_assert(size(priv%mass).gt.0)

   priv%value = ZERO

   cm = ZERO
   n = 1
   do a = 1, priv%n_atoms
      a3 = 3*cv%i(a)
      cm = cm + priv%mass(n)*x(a3 - 2:a3)
      n = n + 1
   end do

   cm = cm/priv%total_mass

   ! populate cm_crd && compute the "norm"
   cm_nrm = ZERO
   n = 1
   do a = 1, priv%n_atoms
      a3 = 3*(n - 1)
      do i = 1, 3
         priv%cm_crd(a3 + i) = x(3*(cv%i(a) - 1) + i) - cm(i)
         cm_nrm = cm_nrm + priv%mass(n)*priv%cm_crd(a3 + i)**2
      end do
      n = n + 1
   end do

   call orientation_q(priv%n_atoms, priv%mass, &
         priv%cm_crd, priv%ref_crd, priv%lambda, priv%quaternion)

!#ifdef MPI
   if (priv%quaternion(1,1) >= 0.0) then
     priv%value = priv%quaternion(3,1)
   else
     priv%value = -1.0 * priv%quaternion(3,1)
   endif

!#endif /* MPI */

   value = priv%value

end function v_QUATERNION2

!=============================================================================

!
! assumes that atom positions have not been
!   changed since last call to 'value()'
!

subroutine f_QUATERNION2(cv, fcv, f)

   use nfe_lib_mod
   use gbl_constants_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: fcv
   double precision, intent(inout) :: f(*)

   integer :: n, a, a3, f3
   double precision :: dSx(4,4), dSy(4,4), dSz(4,4), QL(3)
   double precision :: dSxq1(4), dSyq1(4), dSzq1(4)
   double precision :: a1x, a1y, a1z, tmp
   type(priv_t), pointer :: priv

#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_QUATERNION2)
   nfe_assert(associated(cv%i))

   priv => get_priv(cv)

   nfe_assert(associated(priv))
   nfe_assert(associated(priv%mass))
   nfe_assert(size(priv%mass).gt.0)

#ifdef MPI
   call mpi_bcast(priv%value, 1, MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

      QL(1) = fcv * priv%quaternion(3,2) / &
                    max(priv%lambda(1)-priv%lambda(2),dble(0.00000001))
      QL(2) = fcv * priv%quaternion(3,3) / &
                    max(priv%lambda(1)-priv%lambda(3),dble(0.00000001))
      QL(3) = fcv * priv%quaternion(3,4) / &
                    max(priv%lambda(1)-priv%lambda(4),dble(0.00000001))

      n = 1
      do a = 1, priv%n_atoms
         a3 = 3*n
         f3 = 3*cv%i(a)

         a1x = priv%ref_crd(a3 - 2)
         a1y = priv%ref_crd(a3 - 1)
         a1z = priv%ref_crd(a3)

         dSx(1,1) = a1x
         dSy(1,1) = a1y
         dSz(1,1) = a1z

         dSx(2,1) = 0.0
         dSy(2,1) = -a1z
         dSz(2,1) = a1y

         dSx(1,2) = 0.0
         dSy(1,2) = -a1z
         dSz(1,2) = a1y

         dSx(3,1) = a1z
         dSy(3,1) = 0.0
         dSz(3,1) = -a1x

         dSx(1,3) = a1z
         dSy(1,3) = 0.0
         dSz(1,3) = -a1x

         dSx(4,1) = -a1y
         dSy(4,1) = a1x
         dSz(4,1) = 0.0

         dSx(1,4) = -a1y
         dSy(1,4) = a1x
         dSz(1,4) = 0.0

         dSx(2,2) = a1x
         dSy(2,2) = -a1y
         dSz(2,2) = -a1z

         dSx(3,2) = a1y
         dSy(3,2) = a1x
         dSz(3,2) = 0.0

         dSx(2,3) = a1y
         dSy(2,3) = a1x
         dSz(2,3) = 0.0

         dSx(4,2) = a1z
         dSy(4,2) = 0.0
         dSz(4,2) = a1x

         dSx(2,4) = a1z
         dSy(2,4) = 0.0
         dSz(2,4) = a1x

         dSx(3,3) = -a1x
         dSy(3,3) = a1y
         dSz(3,3) = -a1z

         dSx(4,3) = 0.0
         dSy(4,3) = a1z
         dSz(4,3) = a1y

         dSx(3,4) = 0.0
         dSy(3,4) = a1z
         dSz(3,4) = a1y

         dSx(4,4) = -a1x
         dSy(4,4) = -a1y
         dSz(4,4) = a1z

         dSxq1 = matmul(dSx,priv%quaternion(:,1))
         dSyq1 = matmul(dSy,priv%quaternion(:,1))
         dSzq1 = matmul(dSz,priv%quaternion(:,1))

         f(f3 - 2) = f(f3 - 2) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSxq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSxq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSxq1,priv%quaternion(:,4)) )

         f(f3 - 1) = f(f3 - 1) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSyq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSyq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSyq1,priv%quaternion(:,4)) )

         f(f3) = f(f3) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSzq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSzq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSzq1,priv%quaternion(:,4)) )

         n = n + 1
      end do

end subroutine f_QUATERNION2

!=============================================================================
subroutine b_QUATERNION2(cv, cvno, amass)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)
   double precision, pointer :: coor(:) => null()

   integer :: i, n_atoms, error
   integer :: a, b, a1
   double precision :: cm(3)

   type(priv_t), pointer :: priv

   nfe_assert(cv%type == COLVAR_QUATERNION2)

   ! very basic checks
   if (.not. associated(cv%i)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : no integers found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cv%i)

   if (cv_nr.gt.0) then
    if (.not. associated(cv%r)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_INFO, 'CV #', cvno, &
            ' (ORIENTATION) : no reals (reference coordinates) found'
      NFE_MASTER_ONLY_END
      call terminate()
    end if
   end if 

   ! count the number of atoms
   n_atoms = 0

   do i = 1, size(cv%i)
      if (cv%i(i) == 0) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(i)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : unexpected zero (integer #', i, ')'
            NFE_MASTER_ONLY_END
            call terminate()
      else if (cv%i(i) .lt. 0) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : negative integer'
         NFE_MASTER_ONLY_END
         call terminate()
      else
         n_atoms = n_atoms + 1
      end if
   end do

   if (cv_nr.eq.0) then
      if (associated(coor)) &
          deallocate(coor)
          allocate(coor(3*pmemd_natoms()), stat = error)
          if (error.ne.0) &
           NFE_OUT_OF_MEMORY
           call colvar_read_refcrd(coor, refcrd_file)
           allocate(cv%r(3*size(cv%i)), stat = error)
             if (error.ne.0) &
              NFE_OUT_OF_MEMORY
              do a = 1, n_atoms
               a1=cv%i(a)
               cv%r(1 + 3*(a - 1): 3*a) = coor(1 + 3*(a1 - 1): 3*a1)
              end do
          deallocate(coor)
   end if

   if (size(cv%r) /= 3*n_atoms) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : wrong number of reals'
      NFE_MASTER_ONLY_END
      call terminate()
   end if 

   ! allocate priv_t instance for this variable
   priv => new_priv(cv)

   ! basic checks
   if (n_atoms < 3) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, &
            fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
            (1)//',a,'//pfmt(n_atoms)//',a/)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (ORIENTATION) : too few integers in reference (', 1, ':', n_atoms, ')'
      NFE_MASTER_ONLY_END
      call terminate()
   end if

   do a = 1, n_atoms
      if (cv%i(a) < 1 .or. cv%i(a) > pmemd_natoms()) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
               (cv%i(a))//',a,'//pfmt(pmemd_natoms())//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : integer #', a, ' (', cv%i(a), &
               ') is out of range [1, ', pmemd_natoms(), ']'
         NFE_MASTER_ONLY_END
         call terminate()
      end if
   end do

   ! check for duplicates
   do a = 1, n_atoms
      do b = a + 1, n_atoms
         if (cv%i(a) == cv%i(b)) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
                  (b)//',a,'//pfmt(cv%i(a))//',a/)') &
                  NFE_ERROR, 'CV #', cvno, &
                  ' (ORIENTATION) : integers #', a, ' and #', b, &
                  ' are equal (', cv%i(a), ')'
            NFE_MASTER_ONLY_END
            call terminate()
         end if 
      end do
   end do

   ! allocate/setup

   priv%n_atoms = n_atoms

   allocate(priv%mass(n_atoms), priv%cm_crd(3*n_atoms), &
            priv%ref_crd(3*n_atoms), stat = error)

   if (error /= 0) &
      NFE_OUT_OF_MEMORY

   cm = ZERO
   priv%total_mass = ZERO

   do a = 1, n_atoms
      priv%mass(a) = amass(cv%i(a))
      priv%total_mass = priv%total_mass + priv%mass(a)
      cm = cm + priv%mass(a)*cv%r(1 + 3*(a - 1): 3*a)
   end do

   cm = cm/priv%total_mass

   ! translate reference coordinates to CM frame
   priv%ref_nrm = ZERO
   do a = 1, n_atoms
      do b = 1, 3
         priv%ref_crd(3*(a - 1) + b) = cv%r(3*(a - 1) + b) - cm(b)
         priv%ref_nrm = priv%ref_nrm + priv%mass(a)*priv%ref_crd(3*(a - 1) + b)**2
      end do
   end do


end subroutine b_QUATERNION2

!=============================================================================

subroutine c_QUATERNION2(cv)

   NFE_USE_AFAILED

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(inout) :: cv

   type(priv_t), pointer :: priv

   nfe_assert(cv%type == COLVAR_QUATERNION2)
   nfe_assert(associated(cv%i).and.size(cv%i).gt.0)

   priv => get_priv(cv)

   nfe_assert(associated(priv%mass))
   nfe_assert(associated(priv%cm_crd))
   nfe_assert(associated(priv%ref_crd))

   deallocate(priv%mass, priv%cm_crd, priv%ref_crd)

   call del_priv(cv)

end subroutine c_QUATERNION2

!=============================================================================

subroutine p_QUATERNION2(cv, lun)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   type(priv_t), pointer :: priv
   integer :: a, c
   character(4) :: aname

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_QUATERNION2)
   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   priv => get_priv(cv)
   nfe_assert(priv%n_atoms.gt.0.and.priv%n_atoms.eq.size(cv%i))

   write (unit = lun, fmt = '(a,a)', advance = 'NO') NFE_INFO, 'atoms = ('

   c = 1
   do a = 1, priv%n_atoms

      nfe_assert(cv%i(a) > 0 .and. cv%i(a) <= pmemd_natoms())
      aname = pmemd_atom_name(cv%i(a))

      write (unit = lun, fmt = '('//pfmt(cv%i(a))//',a,a,a)', advance = 'NO') &
         cv%i(a), ' [', trim(aname), ']'

      if (a == priv%n_atoms) then
         write (unit = lun, fmt = '(a)') ')'
      else if (mod(c, 5) == 0) then
         write (unit = lun, fmt = '(a,/a,3x)', advance = 'NO') ',', NFE_INFO
      else
         write (unit = lun, fmt = '(a)', advance = 'NO') ', '
      end if

      c = c + 1
   end do

   if (cv_nr.eq.0) &
      write (unit = lun, fmt = '(a,a,a)') &
         NFE_INFO, 'reference coordinates are loaded from : ', trim(refcrd_file)

   write (unit = lun, fmt = '(a,a)') &
         NFE_INFO, 'reference coordinates :'

   c = 1
   do a = 1, priv%n_atoms
      nfe_assert(c + 2 <= size(cv%r))
      write (unit = lun, fmt = '(a,3x,i5,a,f8.3,a,f8.3,a,f8.3)') &
         NFE_INFO, cv%i(a), ' : ', &
         cv%r(c), ', ', cv%r(c + 1), &
         ', ', cv%r(c + 2)
      c = c + 3
   end do

end subroutine p_QUATERNION2
!=========================================================================
!========================== QUATERNION3 ==================================
!=========================================================================
!
! input:
!
! cv%i = (a1, a2, ..., aN)
!
!     (list of participating atoms)
!
! cv%r = (a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z, ...)
!
!        (reference coordinates)
!
! value = q3
! in which (q0,q1,q2,q3) is orientation quaternion representing optimum rotation w.r.t. reference
!
!=============================================================================

function v_QUATERNION3(cv, x) result(value)

   NFE_USE_AFAILED

   use nfe_lib_mod
   use parallel_dat_mod
   use gbl_constants_mod

   implicit none

   double precision :: value

   type(colvar_t), intent(inout) :: cv

   double precision, intent(in) :: x(*)

!#  include "nfe-mpi.h"

   type(priv_t), pointer :: priv

!#ifdef MPI
   integer :: error
!#endif /* MPI */

   double precision :: cm(3), cm_nrm
   integer :: i, n
! gfortran 4.2 and lower does not support the volatile keyword, but this is
! necessary to work around an intel compiler optimization bug. So turn off
! the volatile keyword for qualifying compilers
!#if !defined(__GNUC__) || (__GNUC__ >= 4 && __GNUC_MINOR__ > 2)
!   integer, volatile :: a, a3 ! work around ifort optimization bug
!#else
   integer :: a, a3
!#endif /* __GFORTRAN__ */

   nfe_assert(cv%type == COLVAR_QUATERNION3)

   priv => get_priv(cv)

   nfe_assert(associated(priv))
   nfe_assert(associated(priv%mass))
   nfe_assert(size(priv%mass).gt.0)

   priv%value = ZERO

   cm = ZERO
   n = 1
   do a = 1, priv%n_atoms
      a3 = 3*cv%i(a)
      cm = cm + priv%mass(n)*x(a3 - 2:a3)
      n = n + 1
   end do

   cm = cm/priv%total_mass

   ! populate cm_crd && compute the "norm"
   cm_nrm = ZERO
   n = 1
   do a = 1, priv%n_atoms
      a3 = 3*(n - 1)
      do i = 1, 3
         priv%cm_crd(a3 + i) = x(3*(cv%i(a) - 1) + i) - cm(i)
         cm_nrm = cm_nrm + priv%mass(n)*priv%cm_crd(a3 + i)**2
      end do
      n = n + 1
   end do

   call orientation_q(priv%n_atoms, priv%mass, &
         priv%cm_crd, priv%ref_crd, priv%lambda, priv%quaternion)

!#ifdef MPI
   if (priv%quaternion(1,1) >= 0.0) then
     priv%value = priv%quaternion(4,1)
   else
     priv%value = -1.0 * priv%quaternion(4,1)
   endif
!#endif /* MPI */

   value = priv%value

end function v_QUATERNION3

!=============================================================================

!
! assumes that atom positions have not been
!   changed since last call to 'value()'
!

subroutine f_QUATERNION3(cv, fcv, f)

   use nfe_lib_mod
   use gbl_constants_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(in) :: cv

   double precision, intent(in) :: fcv
   double precision, intent(inout) :: f(*)


   integer :: n, a, a3, f3
   double precision :: dSx(4,4), dSy(4,4), dSz(4,4), QL(3)
   double precision :: dSxq1(4), dSyq1(4), dSzq1(4)
   double precision :: a1x, a1y, a1z, tmp
   type(priv_t), pointer :: priv

#ifdef MPI
   integer :: error
#endif /* MPI */

   nfe_assert(cv%type == COLVAR_QUATERNION3)
   nfe_assert(associated(cv%i))

   priv => get_priv(cv)

   nfe_assert(associated(priv))
   nfe_assert(associated(priv%mass))
   nfe_assert(size(priv%mass).gt.0)

#ifdef MPI
   call mpi_bcast(priv%value, 1, MPI_DOUBLE_PRECISION, 0, pmemd_comm, error)
   nfe_assert(error.eq.0)
#endif /* MPI */

      QL(1) = fcv * priv%quaternion(4,2) / &
                    max(priv%lambda(1)-priv%lambda(2),dble(0.00000001))
      QL(2) = fcv * priv%quaternion(4,3) / &
                    max(priv%lambda(1)-priv%lambda(3),dble(0.00000001))
      QL(3) = fcv * priv%quaternion(4,4) / &
                    max(priv%lambda(1)-priv%lambda(4),dble(0.00000001))
      n = 1
      do a = 1, priv%n_atoms
         a3 = 3*n
         f3 = 3*cv%i(a)

         a1x = priv%ref_crd(a3 - 2)
         a1y = priv%ref_crd(a3 - 1)
         a1z = priv%ref_crd(a3)

         dSx(1,1) = a1x
         dSy(1,1) = a1y
         dSz(1,1) = a1z

         dSx(2,1) = 0.0
         dSy(2,1) = -a1z
         dSz(2,1) = a1y

         dSx(1,2) = 0.0
         dSy(1,2) = -a1z
         dSz(1,2) = a1y

         dSx(3,1) = a1z
         dSy(3,1) = 0.0
         dSz(3,1) = -a1x

         dSx(1,3) = a1z
         dSy(1,3) = 0.0
         dSz(1,3) = -a1x

         dSx(4,1) = -a1y
         dSy(4,1) = a1x
         dSz(4,1) = 0.0

         dSx(1,4) = -a1y
         dSy(1,4) = a1x
         dSz(1,4) = 0.0

         dSx(2,2) = a1x
         dSy(2,2) = -a1y
         dSz(2,2) = -a1z

         dSx(3,2) = a1y
         dSy(3,2) = a1x
         dSz(3,2) = 0.0

         dSx(2,3) = a1y
         dSy(2,3) = a1x
         dSz(2,3) = 0.0

         dSx(4,2) = a1z
         dSy(4,2) = 0.0
         dSz(4,2) = a1x

         dSx(2,4) = a1z
         dSy(2,4) = 0.0
         dSz(2,4) = a1x

         dSx(3,3) = -a1x
         dSy(3,3) = a1y
         dSz(3,3) = -a1z

         dSx(4,3) = 0.0
         dSy(4,3) = a1z
         dSz(4,3) = a1y

         dSx(3,4) = 0.0
         dSy(3,4) = a1z
         dSz(3,4) = a1y

         dSx(4,4) = -a1x
         dSy(4,4) = -a1y
         dSz(4,4) = a1z

         dSxq1 = matmul(dSx,priv%quaternion(:,1))
         dSyq1 = matmul(dSy,priv%quaternion(:,1))
         dSzq1 = matmul(dSz,priv%quaternion(:,1))

         f(f3 - 2) = f(f3 - 2) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSxq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSxq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSxq1,priv%quaternion(:,4)) )

         f(f3 - 1) = f(f3 - 1) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSyq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSyq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSyq1,priv%quaternion(:,4)) )

         f(f3) = f(f3) &
            + priv%mass(n) &
            * ( QL(1) * dot_product(dSzq1,priv%quaternion(:,2)) &
              + QL(2) * dot_product(dSzq1,priv%quaternion(:,3)) &
              + QL(3) * dot_product(dSzq1,priv%quaternion(:,4)) )

         n = n + 1
      end do

end subroutine f_QUATERNION3

!=============================================================================
subroutine b_QUATERNION3(cv, cvno, amass)

   use nfe_lib_mod
   use parallel_dat_mod

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   double precision,      intent(in)    :: amass(*)
   double precision, pointer :: coor(:) => null()

   integer :: i, n_atoms, error
   integer :: a, b, a1 
   double precision :: cm(3)

   type(priv_t), pointer :: priv

   nfe_assert(cv%type == COLVAR_QUATERNION3)

   ! very basic checks
   if (.not. associated(cv%i)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : no integers found'
      NFE_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cv%i)

   if (cv_nr.gt.0) then
    if (.not. associated(cv%r)) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
            NFE_INFO, 'CV #', cvno, &
            ' (ORIENTATION) : no reals (reference coordinates) found'
      NFE_MASTER_ONLY_END
      call terminate()
    end if
   end if 

   ! count the number of atoms
   n_atoms = 0

   do i = 1, size(cv%i)
      if (cv%i(i) == 0) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(i)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : unexpected zero (integer #', i, ')'
            NFE_MASTER_ONLY_END
            call terminate()
      else if (cv%i(i) .lt. 0) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : negative integer'
         NFE_MASTER_ONLY_END
         call terminate()
      else
         n_atoms = n_atoms + 1
      end if
   end do

   if (cv_nr.eq.0) then
      if (associated(coor)) &
          deallocate(coor)
          allocate(coor(3*pmemd_natoms()), stat = error)
          if (error.ne.0) &
           NFE_OUT_OF_MEMORY
           call colvar_read_refcrd(coor, refcrd_file)
           allocate(cv%r(3*size(cv%i)), stat = error)
             if (error.ne.0) &
              NFE_OUT_OF_MEMORY
              do a = 1, n_atoms
               a1=cv%i(a)
               cv%r(1 + 3*(a - 1): 3*a) = coor(1 + 3*(a1 - 1): 3*a1)
              end do
          deallocate(coor)
   end if
                                     
   if (size(cv%r) /= 3*n_atoms) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : wrong number of reals'
      NFE_MASTER_ONLY_END
      call terminate()
   end if 

   ! allocate priv_t instance for this variable
   priv => new_priv(cv)

   ! basic checks
   if (n_atoms < 3) then
      NFE_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, &
            fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt &
            (1)//',a,'//pfmt(n_atoms)//',a/)') &
            NFE_ERROR, 'CV #', cvno, &
            ' (ORIENTATION) : too few integers in reference (', 1, ':', n_atoms, ')'
      NFE_MASTER_ONLY_END
      call terminate()
   end if

   do a = 1, n_atoms
      if (cv%i(a) < 1 .or. cv%i(a) > pmemd_natoms()) then
         NFE_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
               (cv%i(a))//',a,'//pfmt(pmemd_natoms())//',a/)') &
               NFE_ERROR, 'CV #', cvno, &
               ' (ORIENTATION) : integer #', a, ' (', cv%i(a), &
               ') is out of range [1, ', pmemd_natoms(), ']'
         NFE_MASTER_ONLY_END
         call terminate()
      end if
   end do

   ! check for duplicates
   do a = 1, n_atoms
      do b = a + 1, n_atoms
         if (cv%i(a) == cv%i(b)) then
            NFE_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
                  (b)//',a,'//pfmt(cv%i(a))//',a/)') &
                  NFE_ERROR, 'CV #', cvno, &
                  ' (ORIENTATION) : integers #', a, ' and #', b, &
                  ' are equal (', cv%i(a), ')'
            NFE_MASTER_ONLY_END
            call terminate()
         end if
      end do
   end do

   ! allocate/setup

   priv%n_atoms = n_atoms

   allocate(priv%mass(n_atoms), priv%cm_crd(3*n_atoms), &
            priv%ref_crd(3*n_atoms), stat = error)

   if (error /= 0) &
      NFE_OUT_OF_MEMORY

   cm = ZERO
   priv%total_mass = ZERO

   do a = 1, n_atoms
      priv%mass(a) = amass(cv%i(a))
      priv%total_mass = priv%total_mass + priv%mass(a)
      cm = cm + priv%mass(a)*cv%r(1 + 3*(a - 1): 3*a)
   end do

   cm = cm/priv%total_mass

   ! translate reference coordinates to CM frame
   priv%ref_nrm = ZERO
   do a = 1, n_atoms
      do b = 1, 3
         priv%ref_crd(3*(a - 1) + b) = cv%r(3*(a - 1) + b) - cm(b)
         priv%ref_nrm = priv%ref_nrm + priv%mass(a)*priv%ref_crd(3*(a - 1) + b)**2
      end do
   end do


end subroutine b_QUATERNION3

!=============================================================================

subroutine c_QUATERNION3(cv)

   NFE_USE_AFAILED

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(inout) :: cv

   type(priv_t), pointer :: priv

   nfe_assert(cv%type == COLVAR_QUATERNION3)
   nfe_assert(associated(cv%i).and.size(cv%i).gt.0)

   priv => get_priv(cv)

   nfe_assert(associated(priv%mass))
   nfe_assert(associated(priv%cm_crd))
   nfe_assert(associated(priv%ref_crd))

   deallocate(priv%mass, priv%cm_crd, priv%ref_crd)

   call del_priv(cv)

end subroutine c_QUATERNION3

!=============================================================================

subroutine p_QUATERNION3(cv, lun)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(in) :: cv
   integer, intent(in) :: lun

   type(priv_t), pointer :: priv
   integer :: a, c
   character(4) :: aname

   nfe_assert(is_master())
   nfe_assert(cv%type == COLVAR_QUATERNION3)
   nfe_assert(associated(cv%i))
   nfe_assert(associated(cv%r))

   priv => get_priv(cv)
   nfe_assert(priv%n_atoms.gt.0.and.priv%n_atoms.eq.size(cv%i))

   write (unit = lun, fmt = '(a,a)', advance = 'NO') NFE_INFO, 'atoms = ('

   c = 1
   do a = 1, priv%n_atoms

      nfe_assert(cv%i(a) > 0 .and. cv%i(a) <= pmemd_natoms())
      aname = pmemd_atom_name(cv%i(a))

      write (unit = lun, fmt = '('//pfmt(cv%i(a))//',a,a,a)', advance = 'NO') &
         cv%i(a), ' [', trim(aname), ']'

      if (a == priv%n_atoms) then
         write (unit = lun, fmt = '(a)') ')'
      else if (mod(c, 5) == 0) then
         write (unit = lun, fmt = '(a,/a,3x)', advance = 'NO') ',', NFE_INFO
      else
         write (unit = lun, fmt = '(a)', advance = 'NO') ', '
      end if

      c = c + 1
   end do

   if (cv_nr.eq.0) &
      write (unit = lun, fmt = '(a,a,a)') &
         NFE_INFO, 'reference coordinates are loaded from : ', trim(refcrd_file)

   write (unit = lun, fmt = '(a,a)') &
         NFE_INFO, 'reference coordinates :'

   c = 1
   do a = 1, priv%n_atoms
      nfe_assert(c + 2 <= size(cv%r))
      write (unit = lun, fmt = '(a,3x,i5,a,f8.3,a,f8.3,a,f8.3)') &
         NFE_INFO, cv%i(a), ' : ', &
         cv%r(c), ', ', cv%r(c + 1), &
         ', ', cv%r(c + 2)
      c = c + 3
   end do

end subroutine p_QUATERNION3
!============================================================================
! searches the list
function get_priv(cv) result(ptr)

   NFE_USE_AFAILED

   implicit none

   type(colvar_t), intent(in) :: cv
   type(priv_t), pointer :: ptr

   ptr => priv_list
   do while (associated(ptr))
      if (ptr%tag.eq.cv%tag) &
         exit
      ptr => ptr%next
   end do

   nfe_assert(associated(ptr))

end function get_priv
!============================================================================
! allocates and appends to the list
function new_priv(cv) result(ptr)

   use nfe_lib_mod

   implicit none

   type(colvar_t), intent(inout) :: cv

   type(priv_t), pointer :: ptr
   type(priv_t), pointer :: head

   integer :: error

   allocate(ptr, stat = error)
   if (error.ne.0) &
      NFE_OUT_OF_MEMORY

   ptr%next => null()

   if (.not.associated(priv_list)) then
      ptr%tag = 0
      priv_list => ptr
   else
      head => priv_list
      do while (associated(head%next))
         head => head%next
      end do
      ptr%tag = head%tag + 1
      head%next => ptr
   end if

   cv%tag = ptr%tag

end function new_priv
!=============================================================================
! removes from the list and deallocates
subroutine del_priv(cv)

   NFE_USE_AFAILED


   implicit none

   type(colvar_t), intent(in) :: cv

   type(priv_t), pointer :: curr, prev

   nfe_assert(associated(priv_list))

   curr => priv_list
   if (curr%tag.eq.cv%tag) then
      priv_list => curr%next
   else
      prev => curr
      curr => curr%next
      do while (associated(curr))
         if (curr%tag.eq.cv%tag) then
            prev%next => curr%next
            exit
         end if
         prev => curr
         curr => curr%next
      end do
   end if

   nfe_assert(associated(curr))
   nfe_assert(curr%tag.eq.cv%tag)

   deallocate(curr)

end subroutine del_priv
!=============================================================================
!============================ MATH FUNCTIONS =================================
!=============================================================================
pure double precision function distance(r1, r2)

   implicit none

   double precision, intent(in) :: r1(3), r2(3)

   distance = norm3(r1 - r2)

end function distance

!=============================================================================

subroutine distance_d(r1, r2, d1, d2)

   implicit none

   double precision, intent(in)  :: r1(3), r2(3)
   double precision, intent(out) :: d1(3), d2(3)

   double precision :: dr(3)

   dr = r1 - r2

   d1 = dr/norm3(dr)
   d2 = - d1

end subroutine distance_d

!=============================================================================

pure double precision function angle(r1, r2, r3)

   implicit none

   double precision, intent(in) :: r1(3), r2(3), r3(3)

   angle = acos(cosine(r1 - r2, r3 - r2))

end function angle

!=============================================================================

subroutine angle_d(r1, r2, r3, d1, d2, d3)

   use nfe_lib_mod, only : ONE

   implicit none

   double precision, intent(in)  :: r1(3), r2(3), r3(3)
   double precision, intent(out) :: d1(3), d2(3), d3(3)

   double precision :: c, d

   c = cosine_d(r1 - r2, r3 - r2, d1, d3)
   d = -ONE/sqrt(ONE - c*c)

   d1 = d*d1
   d3 = d*d3

   d2 = - (d1 + d3)

end subroutine angle_d

!=============================================================================

!
! for A == B == C == D, a torsion along B == C is arccos([ABxBC]*[CDx(-BC)])
!        and its sign is given by sign(BC*[[ABxBC]x[CDx(-BC)]])
!

pure double precision function torsion(r1, r2, r3, r4)

   implicit none

   double precision, intent(in) :: r1(3), r2(3), r3(3), r4(3)

   torsion = torsion3(r2 - r1, r3 - r2, r4 - r3)

end function torsion

!=============================================================================

subroutine torsion_d(r1, r2, r3, r4, d1, d2, d3, d4)

   implicit none

   double precision, intent(in)  :: r1(3), r2(3), r3(3), r4(3)
   double precision, intent(out) :: d1(3), d2(3), d3(3), d4(3)

   double precision :: t(3)

   call torsion3_d(r2 - r1, r3 - r2, r4 - r3, d1, t, d4)

   d2 = d1 - t
   d1 = - d1
   d3 = t - d4

end subroutine torsion_d

!=============================================================================

pure double precision function dot3(v1, v2)

   implicit none

   double precision, intent(in) :: v1(3), v2(3)

   dot3 = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)

end function dot3

!=============================================================================

pure double precision function norm3(v)

   implicit none

   double precision, intent(in) :: v(3)

   norm3 = sqrt(v(1)**2 + v(2)**2 + v(3)**2)

end function norm3

!=============================================================================

pure function cross3(v1, v2) result(cross)

   implicit none

   double precision             :: cross(3)
   double precision, intent(in) :: v1(3), v2(3)

   cross(1) = v1(2)*v2(3) - v1(3)*v2(2)
   cross(2) = v1(3)*v2(1) - v1(1)*v2(3)
   cross(3) = v1(1)*v2(2) - v1(2)*v2(1)

end function cross3

!=============================================================================

pure double precision function cosine(v1, v2)

   use nfe_lib_mod, only : ONE

   implicit none

   double precision, intent(in) :: v1(3), v2(3)

   cosine = dot3(v1, v2)/(norm3(v1)*norm3(v2))

   if (cosine > ONE) cosine = ONE
   if (cosine < -ONE) cosine = -ONE

end function cosine

!=============================================================================

double precision function cosine_d(v1, v2, d1, d2)

   use nfe_lib_mod, only : ONE

   implicit none

   double precision, intent(in)  :: v1(3), v2(3)
   double precision, intent(out) :: d1(3), d2(3)

   double precision :: n1, n2, p12

   n1  = norm3(v1)
   n2  = norm3(v2)
   p12 = dot3(v1, v2)

   cosine_d = p12/(n1*n2)

   d1 = (v2 - v1*p12/(n1**2))/(n1*n2)
   d2 = (v1 - v2*p12/(n2**2))/(n1*n2)

   if (cosine_d > ONE) cosine_d = ONE
   if (cosine_d < -ONE) cosine_d = -ONE

end function cosine_d

!=============================================================================

pure double precision function torsion3(v1, v2, v3)

   implicit none

   double precision, intent(in) :: v1(3), v2(3), v3(3)

   double precision :: n1(3), n2(3)

   n1 = cross3(v1, v2)
   n2 = cross3(v2, v3)

   torsion3 = sign(acos(cosine(n1, n2)), dot3(v2, cross3(n1, n2)))

end function torsion3

!=============================================================================

subroutine torsion3_d(v1, v2, v3, d1, d2, d3)

   use nfe_lib_mod, only : ONE

   implicit none

   double precision, intent(in)  :: v1(3), v2(3), v3(3)
   double precision, intent(out) :: d1(3), d2(3), d3(3)

   double precision, parameter :: epsilon = 1.0d-8
   double precision :: n1(3), n2(3), dc1(3), dc2(3), c, c2, s, d

   n1 = cross3(v1, v2)
   n2 = cross3(v2, v3)

   c = cosine_d(n1, n2, dc1, dc2)
   s = dot3(v2, cross3(n1, n2))

   c2 = c*c - epsilon
   d = sign(ONE/sqrt(ONE - c2), s)

   d1 = d*cross3(dc1, v2)
   d2 = d*(cross3(v1, dc1) + cross3(dc2, v3))
   d3 = d*cross3(v2, dc2)

end subroutine torsion3_d

!=============================================================================
#if 0
!
! computes out-of-plane (puckering) coordinates as described in
!        http://dx.doi.org/10.1021/ja00839a011
!

subroutine pucker_z(R, z)

   use nfe_lib_mod

   implicit none

   double precision, intent(in) :: R(:,:)
   double precision, intent(out) :: z(:)

   double precision :: R0(3), Rt(3), Rc(3), Rs(3), phase, N(3), pi
   integer :: nring, k

   nring = size(R, 2)
   nfe_assert(nring.gt.0)

   pi = 4*atan(dble(1))

   R0 = ZERO
   do k = 1, nring
      R0 = R0 + R(:, k)
   end do
   R0 = R0/nring

   Rc = ZERO
   Rs = ZERO

   do k = 1, nring
      Rt = R(:, k) - R0
      phase = 2*pi*(k - 1)/dble(nring)
      Rc = Rc + Rt*cos(phase)
      Rs = Rs + Rt*sin(phase)
   end do

   N = cross3(Rs, Rc)
   N = N/norm3(N)

   do k = 1, nring
      Rt = R(:, k) - R0
      z(k) = dot3(Rt, N)
   end do

end subroutine pucker_z

!-----------------------------------------------------------------------------

subroutine pucker_dz(R, z, dz)

   use nfe_lib_mod

   implicit none

   double precision, intent(in) :: R(:,:)
   double precision, intent(out) :: z(:), dz(:,:,:)

   double precision :: R0(3), Rt(3), Rc(3), Rs(3), pi
   double precision :: N(3), phi(3), z0(3), NN, phase

   integer :: nring, k, m

   nring = size(R, 2)
   nfe_assert(nring.gt.0)

   pi = 4*atan(dble(1))

   R0 = ZERO
   do k = 1, nring
      R0 = R0 + R(:, k)
   end do
   R0 = R0/nring

   Rc = ZERO
   Rs = ZERO

   do k = 1, nring
      Rt = R(:, k) - R0
      phase = 2*pi*(k - 1)/dble(nring)
      Rc = Rc + Rt*cos(phase)
      Rs = Rs + Rt*sin(phase)
   end do

   N = cross3(Rs, Rc)
   NN = norm3(N)
   N = N/NN

   do m = 1, nring
      Rt = R(:, m) - R0
      z(m) = dot3(Rt, N)
      Rt = Rt - z(m)*N
      z0 = ZERO
      do k = 1, nring
         phase = 2*pi*(k - 1)/dble(nring)
         phi = (Rs*cos(phase) - Rc*sin(phase))/NN
         dz(:, k, m) = cross3(Rt, phi)
         if (m.eq.k) &
            dz(:, k, m) = dz(:, k, m) + N
         z0 = z0 + dz(:, k, m)
      end do
      z0 = z0/nring
      do k = 1, nring
         dz(:, k, m) = dz(:, k, m) - z0
      end do
   end do

end subroutine pucker_dz
#endif

!=============================================================================
#ifdef NFE_RMSD_CANNED
double precision function rmsd_canned(m, x1, x2, g)

   use nfe_lib_mod

   implicit none

   double precision, intent(in)    :: m(:)
   double precision, intent(inout) :: x1(:), x2(:) ! destroyed upon return

   double precision, intent(out), optional :: g(:)

   integer   :: n, a, a3, i
   double precision :: mass, c1(3), c2(3), q(4), lambda, x1n, x2n, tmp, U(3,3)

   n = size(m)

   nfe_assert(n > 1)

   nfe_assert(3*n == size(x1))
   nfe_assert(3*n == size(x2))

   !
   ! find centers of mass
   !

   c1 = ZERO
   c2 = ZERO

   mass = ZERO

   do a = 1, n
      a3 = 3*(a - 1)
      mass = mass + m(a)
      do i = 1, 3
         c1(i) = c1(i) + m(a)*x1(a3 + i)
         c2(i) = c2(i) + m(a)*x2(a3 + i)
      end do
   end do

   nfe_assert(mass > ZERO)

   c1 = c1/mass
   c2 = c2/mass

   ! center x1, x2 && find "norms"

   x1n = ZERO
   x2n = ZERO

   do a = 1, n
      a3 = 3*(a - 1)
      do i = 1, 3
         x1(a3 + i) = x1(a3 + i) - c1(i)
         x1n = x1n + m(a)*x1(a3 + i)**2
         x2(a3 + i) = x2(a3 + i) - c2(i)
         x2n = x2n + m(a)*x2(a3 + i)**2
      end do
   end do

   call rmsd_q(n, m, x1, x2, lambda, q)

   rmsd_canned = sqrt(max(ZERO, ((x1n + x2n) - 2*lambda))/mass)

   ! g is w.r.t x1

   if (present(g)) then
      nfe_assert(3*n == size(g))
      call rmsd_q2u(q, U)

      tmp = ONE/mass/max(rmsd_canned, dble(0.000001))
      do a = 1, n
         a3 = 3*a
         g(a3 - 2:a3) = m(a)*tmp*(x1(a3 - 2:a3) - matmul(U, x2(a3 - 2:a3)))
      end do
   end if ! present g

end function rmsd_canned
#endif /* NFE_RMSD_CANNED */

!=============================================================================

!
! finds optimal rotation; size(w) == n; size(x?) == 3*n
!

subroutine rmsd_q(n, w, x1, x2, lambda, q)

   use nfe_lib_mod

   implicit none

   integer,   intent(in) :: n
   double precision, intent(in) :: w(*), x1(*), x2(*)

   double precision, intent(out) :: lambda, q(4)

   integer   :: a, a3, i, j
   double precision :: R(4,4), S(4,4)

   ! calculate the R matrix

   R = ZERO

   do a = 0, n - 1
      a3 = 3*a
      do i = 1, 3
         do j = 1, 3
            R(i,j) = R(i,j) + w(a + 1)*x1(a3 + i)*x2(a3 + j)
         end do
      end do
   end do

   ! S matrix

   S(1,1) = R(1,1) + R(2,2) + R(3,3)
   S(2,1) = R(2,3) - R(3,2)
   S(3,1) = R(3,1) - R(1,3)
   S(4,1) = R(1,2) - R(2,1)

   S(1,2) = S(2,1)
   S(2,2) = R(1,1) - R(2,2) - R(3,3)
   S(3,2) = R(1,2) + R(2,1)
   S(4,2) = R(1,3) + R(3,1)

   S(1,3) = S(3,1)
   S(2,3) = S(3,2)
   S(3,3) =-R(1,1) + R(2,2) - R(3,3)
   S(4,3) = R(2,3) + R(3,2)

   S(1,4) = S(4,1)
   S(2,4) = S(4,2)
   S(3,4) = S(4,3)
   S(4,4) =-R(1,1) - R(2,2) + R(3,3) 

   call dstmev(S, lambda, q)

end subroutine rmsd_q

!=============================================================================

! add variable control

subroutine rmsd_q1(n, state, w, x1, x2, lambda, q)

   use nfe_lib_mod

   implicit none

   integer,   intent(in) :: n, state(*)
   double precision, intent(in) :: w(*), x1(*), x2(*)

   double precision, intent(out) :: lambda, q(4)

   integer   :: a, a3, i, j
   double precision :: R(4,4), S(4,4)

   ! calculate the R matrix

   R = ZERO

   do a = 0, n - 1
   	  if(state(a+1)==0) cycle
      a3 = 3*a
      do i = 1, 3
         do j = 1, 3
            R(i,j) = R(i,j) + w(a + 1)*x1(a3 + i)*x2(a3 + j)
         end do
      end do
   end do

   ! S matrix

   S(1,1) = R(1,1) + R(2,2) + R(3,3)
   S(2,1) = R(2,3) - R(3,2)
   S(3,1) = R(3,1) - R(1,3)
   S(4,1) = R(1,2) - R(2,1)

   S(1,2) = S(2,1)
   S(2,2) = R(1,1) - R(2,2) - R(3,3)
   S(3,2) = R(1,2) + R(2,1)
   S(4,2) = R(1,3) + R(3,1)

   S(1,3) = S(3,1)
   S(2,3) = S(3,2)
   S(3,3) =-R(1,1) + R(2,2) - R(3,3)
   S(4,3) = R(2,3) + R(3,2)

   S(1,4) = S(4,1)
   S(2,4) = S(4,2)
   S(3,4) = S(4,3)
   S(4,4) =-R(1,1) - R(2,2) + R(3,3) 

   call dstmev(S, lambda, q)

end subroutine rmsd_q1

!=============================================================================


!
! This subroutine constructs (transposed) rotation matrix U from quaternion q.
! 

subroutine rmsd_q2u(q, U)

   use nfe_lib_mod

   implicit none

   double precision, intent(in)  :: q(*) ! 4
   double precision, intent(out) :: U(3,3)

   double precision :: b0, b1, b2, b3
   double precision :: q00, q01, q02, q03
   double precision :: q11, q12, q13, q22, q23, q33

   b0 = q(1) + q(1)
   b1 = q(2) + q(2)
   b2 = q(3) + q(3)
   b3 = q(4) + q(4)

   q00 = b0*q(1) - ONE
   q01 = b0*q(2)
   q02 = b0*q(3)
   q03 = b0*q(4)

   q11 = b1*q(2)
   q12 = b1*q(3)
   q13 = b1*q(4)

   q22 = b2*q(3)
   q23 = b2*q(4)

   q33 = b3*q(4)

   U(1,1) = q00 + q11
   U(2,1) = q12 - q03
   U(3,1) = q13 + q02

   U(1,2) = q12 + q03
   U(2,2) = q00 + q22
   U(3,2) = q23 - q01

   U(1,3) = q13 - q02
   U(2,3) = q23 + q01
   U(3,3) = q00 + q33

end subroutine rmsd_q2u

!=============================================================================
!
! This subroutine constructs (UNtransposed) rotation matrix U from quaternion q.
! 

subroutine rmsd_q3u(q, U)

   use nfe_lib_mod

   implicit none

   double precision, intent(in)  :: q(*) ! 4
   double precision, intent(out) :: U(3,3)

   double precision :: b0, b1, b2, b3
   double precision :: q00, q01, q02, q03
   double precision :: q11, q12, q13, q22, q23, q33

   b0 = q(1) + q(1)
   b1 = q(2) + q(2)
   b2 = q(3) + q(3)
   b3 = q(4) + q(4)

   q00 = b0*q(1) - ONE
   q01 = b0*q(2)
   q02 = b0*q(3)
   q03 = b0*q(4)

   q11 = b1*q(2)
   q12 = b1*q(3)
   q13 = b1*q(4)

   q22 = b2*q(3)
   q23 = b2*q(4)

   q33 = b3*q(4)

   U(1,1) = q00 + q11
   U(1,2) = q12 - q03
   U(1,3) = q13 + q02

   U(2,1) = q12 + q03
   U(2,2) = q00 + q22
   U(2,3) = q23 - q01

   U(3,1) = q13 - q02
   U(3,2) = q23 + q01
   U(3,3) = q00 + q33

end subroutine rmsd_q3u

!=============================================================================
! orientation_q using for quaternions

subroutine orientation_q(n, w, x1, x2, lambda, q)

   !use gbl_constants_mod
   use nfe_lib_mod

   implicit none

   integer,   intent(in) :: n
   double precision, intent(in) :: w(*), x1(*), x2(*)
   double precision, intent(out) :: lambda(4), q(4,4)
   double precision :: R(4,4), S(4,4)
   integer   :: a, a3, i, j

   ! calculate the R matrix
   R = ZERO
   do a = 0, n - 1
      a3 = 3*a
      do i = 1, 3
         do j = 1, 3
            !R(i,j) = R(i,j) + w(a+1)*x1(a3 + i)*x2(a3 + j)
             R(i,j) = R(i,j) + x1(a3 + i)*x2(a3 + j)
         end do
      end do
   end do
      ! S matrix

   S(1,1) = R(1,1) + R(2,2) + R(3,3)
   S(2,1) = R(2,3) - R(3,2)
   S(3,1) = R(3,1) - R(1,3)
   S(4,1) = R(1,2) - R(2,1)

   S(1,2) = S(2,1)
   S(2,2) = R(1,1) - R(2,2) - R(3,3)
   S(3,2) = R(1,2) + R(2,1)
   S(4,2) = R(1,3) + R(3,1)

   S(1,3) = S(3,1)
   S(2,3) = S(3,2)
   S(3,3) =-R(1,1) + R(2,2) - R(3,3)
   S(4,3) = R(2,3) + R(3,2)

   S(1,4) = S(4,1)
   S(2,4) = S(4,2)
   S(3,4) = S(4,3)
   S(4,4) =-R(1,1) - R(2,2) + R(3,3)
   call dstmev4(S, lambda, q)


end subroutine orientation_q


!=============================================================================

!
! a simple subroutine to compute the leading eigenvalue and eigenvector
! of a symmetric, traceless 4x4 matrix A by an inverse power iteration:
! (1) the matrix is converted to tridiagonal form by 3 Givens
! rotations;  V*A*V' = T
! (2) Gershgorin's theorem is used to estimate a lower
! bound for the leading negative eigenvalue:
! lambda_1 > g=min(T11-t12,-t21+T22-t23,-t32+T33-t34,-t43+T44)
!          =
! where tij=abs(Tij)
! (3) Form the positive definite matrix 
!     B = T-gI
! (4) Use svd (algorithm svdcmp from "Numerical Recipes")
!     to compute eigenvalues and eigenvectors for SPD matrix B
! (5) Shift spectrum back and keep leading singular vector
!     and largest eigenvalue.
! (6) Convert eigenvector to original matrix A, through 
!     multiplication by V'.  
!

subroutine dstmev(A, lambda, evec)

   implicit none

   double precision, intent(in)  :: A(4,4)
   double precision, intent(out) :: lambda, evec(4)

   double precision :: T(4,4), V(4,4), SV(4,4), SW(4)

   integer :: i, max_loc(1)

   ! (I) Convert to tridiagonal form, keeping similarity transform
   !            (a product of 3 Givens rotations)
   call givens4(A, T, V)

   ! (II) Estimate lower bound of smallest eigenvalue by Gershgorin's theorem
   lambda = min(T(1,1) - abs(T(1,2)), &
                T(2,2) - abs(T(2,1)) - abs(T(2,3)), &
                T(3,3) - abs(T(3,2)) - abs(T(3,4)), &
                T(4,4) - abs(T(4,3)))

   ! (III) Form positive definite matrix  T = lambda*I - T
   do i = 1, 4
      T(i,i) = T(i,i) - lambda
   end do

   ! (IV) Compute singular values/vectors of SPD matrix B
   call svdcmp(T, SW, SV)

   ! (V) Shift spectrum back
   max_loc = maxloc(SW) 
   lambda = lambda + SW(max_loc(1))

   ! (VI) Convert eigenvector to original coordinates: (V is transposed!)
   evec = matmul(V, SV(:, max_loc(1)))

end subroutine dstmev

!============================================================================

subroutine dstmev4(A, lambda, evec)

   implicit none

   double precision, intent(in)  :: A(4,4)
   double precision, intent(out) :: lambda(4), evec(4,4)

   double precision :: T(4,4), V(4,4), SV(4,4), SW(4)
   double precision :: U(4,4), SVT(4,4), work(20), norm1, norm2
   integer :: lwork, info, iwork, nrot

   integer :: i, ip, iq


   lwork=size(work)
   iwork=size(work)
   info=0

   call jacobi(A, 4, 4, SW, U, nrot)
   call eigsrt (SW, U, 4, 4)

   lambda(4) = SW(4)
   lambda(3) = SW(3)
   lambda(2) = SW(2)
   lambda(1) = SW(1)

   ! SV = trans(U)
   SV = transpose (U)

   ! Normalization
   do ip = 1, 4
      norm2 = 0.0
      do iq = 1, 4
         norm2 = norm2 + SV(ip,iq)*SV(ip,iq)
      enddo
      norm1 = sqrt(norm2)
      do iq = 1, 4
         SV(:,iq) = SV(:,iq)/norm1
      enddo
   enddo

   evec(:,1) = SV(:,1)
   evec(:,2) = SV(:,2)
   evec(:,3) = SV(:,3)
   evec(:,4) = SV(:,4)

end subroutine dstmev4

!=============================================================================

!
! performs givens rotations to reduce symmetric 4x4 matrix to tridiagonal
!

subroutine givens4(S, T, V)

   use nfe_lib_mod

   implicit none

   double precision, dimension(4,4), intent(in)  :: S
   double precision, dimension(4,4), intent(out) :: T,V

   double precision :: c1, c2, c3, s1, s2, s3, r1, r2, r3, c1c2, s1c2

   T = S
   V = ZERO

   ! zero out entries T(4,1) and T(1,4)
   ! compute cos and sin of rotation angle in the 3-4 plane

   r1 = pythag(T(3,1), T(4,1))

   if (r1 .ne. ZERO) then
      c1 = T(3,1)/r1
      s1 = T(4,1)/r1

      V(3,3) = c1
      V(3,4) = s1
      V(4,3) =-s1
      V(4,4) = c1

      T(3,1) = r1
      T(4,1) = ZERO

      T(3:4,2:4) = matmul(V(3:4,3:4), T(3:4,2:4))
      T(1:2,3:4) = transpose(T(3:4,1:2))
      T(3:4,3:4) = matmul(T(3:4,3:4), transpose(V(3:4,3:4)))
   else
      c1 = ONE
      s1 = ZERO
   end if

   ! zero out entries T(3,1) and T(1,3)
   ! compute cos and sin of rotation angle in the 2-3 plane

   r2 = pythag(T(3,1), T(2,1))

   if (r2 .ne. ZERO) then
      c2 = T(2,1)/r2
      s2 = T(3,1)/r2

      V(2,2) = c2
      V(2,3) = s2
      V(3,2) =-s2
      V(3,3) = c2

      T(2,1) = r2
      T(3,1) = ZERO

      T(2:3,2:4) = matmul(V(2:3,2:3), T(2:3,2:4))
      T(1,2:3)   = T(2:3,1)
      T(4,2:3)   = T(2:3,4)
      T(2:3,2:3) = matmul(T(2:3,2:3), transpose(V(2:3,2:3)))
   else
      c2 = ONE
      s2 = ZERO
   end if

   ! zero out entries T(4,2) and T(2,4)
   ! compute cos and sin of rotation angle in the 3-4 plane

   r3 = pythag(T(4,2), T(3,2))

   if (r3 .ne. ZERO) then
      c3 = T(3,2)/r3
      s3 = T(4,2)/r3

      V(3,3) = c3
      V(3,4) = s3
      V(4,3) =-s3
      V(4,4) = c3

      T(3,2) = r3
      T(4,2) = ZERO

      T(3:4,3:4) = matmul(V(3:4,3:4), T(3:4,3:4))
      T(1:2,3:4) = transpose(T(3:4,1:2))
      T(3:4,3:4) = matmul(T(3:4,3:4), transpose(V(3:4,3:4)))
   else
      c3 = ONE
      s3 = ZERO
   end if

   ! compute net rotation matrix (accumulate similarity for evec. computation)
   ! To save transposing later, This is the transpose!

   V(1,1)   = ONE
   V(1,2:4) = ZERO
   V(2:4,1) = ZERO

   V(2,2) = c2
   V(3,2) = c1*s2
   V(4,2) = s1*s2

   c1c2 = c1*c2
   s1c2 = s1*c2

   V(2,3) =-s2*c3
   V(3,3) = c1c2*c3 - s1*s3
   V(4,3) = s1c2*c3 + c1*s3
   V(2,4) = s2*s3
   V(3,4) =-c1c2*s3 - s1*c3
   V(4,4) =-s1c2*s3 + c1*c3

end subroutine givens4

!============================================================================

subroutine svdcmp(a, w, v)

   use nfe_lib_mod

   implicit none

   integer, parameter :: N = 4

   double precision, intent(inout) :: a(N,*)
   double precision, intent(out)   :: v(N,*), w(*)

   integer :: i, its, j, jj, k, l, nm

   double precision :: anorm, c, f, g, h, s, scale, x, y, z, rv1(2*N)

   g = ZERO
   scale = ZERO
   anorm = ZERO

   nm = 0 ! for g95

   do i = 1, N

      l = i + 1
      rv1(i) = scale*g

      g = ZERO
      s = ZERO
      scale = ZERO

      do k = i, N
         scale = scale + abs(a(k,i))
      end do

      if (scale .ne. ZERO) then
         do k = i, N
            a(k,i) = a(k,i)/scale
            s = s + a(k,i)*a(k,i)
         end do

         f = a(i,i)
         g =-sign(sqrt(s),f)
         h = f*g - s
         a(i,i) = f - g

         do j = l, N 
            s = ZERO
            do k = i, N
               s = s + a(k,i)*a(k,j)
            end do

            f = s/h
            do k = i, N
               a(k,j) = a(k,j) + f*a(k,i)
            end do
         end do

         do k = i, N
            a(k,i) = scale*a(k,i)
         end do
      endif ! scale .ne. ZERO

      w(i) = scale*g
      g = ZERO
      s = ZERO
      scale = ZERO

      if (i .ne. N) then
         do k = l, N
            scale = scale + abs(a(i,k))
         end do
         if (scale .ne. ZERO) then
            do k = l, N
               a(i,k) = a(i,k)/scale
               s = s + a(i,k)*a(i,k)
            end do
            f = a(i,l)
            g =-sign(sqrt(s),f)
            h = f*g - s
            a(i,l) = f - g
            do k = l, N
               rv1(k) = a(i,k)/h
            end do
            do j = l, N
               s = ZERO
               do k = l, N
                  s = s + a(j,k)*a(i,k)
               end do
               do k = l, N
                  a(j,k) = a(j,k) + s*rv1(k)
               end do
            end do
            do k = l, N
               a(i,k) = scale*a(i,k)
            end do
         endif
      endif
      anorm = max(anorm, (abs(w(i)) + abs(rv1(i))))
   end do

   do i = N, 1, -1
      if (i .lt. N) then
         if (g .ne. ZERO) then
            do j = l, N
               v(j,i) = (a(i,j)/a(i,l))/g
            end do
            do j = l, N
               s = ZERO
               do k = l, N
                  s = s + a(i,k)*v(k,j)
               end do
               do k = l, N
                  v(k,j) = v(k,j) + s*v(k,i)
               end do
            end do
         endif ! g .ne. ZERO
         do j = l, N
            v(i,j) = ZERO
            v(j,i) = ZERO
         end do
      endif
      v(i,i) = ONE
      g = rv1(i)
      l = i
   end do

   do i = N, 1, -1
      l = i + 1
      g = w(i)
      do j = l, N
         a(i,j) = ZERO
      end do
      if (g .ne. ZERO) then
         g = ONE/g
         do j = l, N
            s = ZERO
            do k = l, N
               s = s + a(k,i)*a(k,j)
            end do
            f = (s/a(i,i))*g
            do k = i, N
               a(k,j) = a(k,j) + f*a(k,i)
            end do
         end do
         do j = i, N
            a(j,i) = a(j,i)*g
         end do
      else
         do j = i, N
            a(j,i) = ZERO
         end do
      endif ! g .ne. ZERO
      a(i,i) = a(i,i) + ONE
   end do

   do k = N, 1, -1
      do its = 1, 30
         do l = k, 1, -1
            nm = l - 1
            if ((abs(rv1(l)) + anorm) .eq. anorm) &
               goto 2
            if ((abs(w(nm)) + anorm) .eq. anorm) &
               goto 1
         end do
1        c = ZERO
         s = ONE
         do i = l, k
            f = s*rv1(i)
            rv1(i) = c*rv1(i)
            if ((abs(f) + anorm) .eq. anorm) &
               goto 2
            g = w(i)
            h = pythag(f,g)
            w(i) = h
            h = ONE/h
            c = (g*h)
            s =-(f*h)
            do j = 1, N     
               y = a(j,nm)
               z = a(j,i)
               a(j,nm) = (y*c) + (z*s)
               a(j,i) = -(y*s) + (z*c)
            end do
         end do
2        z = w(k)
         if (l .eq. k) then
            if (z .lt. ZERO) then
               w(k) = -z
               do j = 1, N
                  v(j,k) = -v(j,k)
               end do
            end if
            goto 3
         end if
         if (its .eq. 30) &
            call fatal('nfe_rmsd: no convergence in svdcmp()')
         x = w(l)
         nm = k - 1
         y = w(nm)
         g = rv1(nm)
         h = rv1(k)
         f = ((y - z)*(y + z) + (g - h)*(g + h))/(2*h*y)
         g = pythag(f, ONE)
         f = ((x - z)*(x + z) + h*((y/(f + sign(g, f))) - h))/x
         c = ONE
         s = ONE
         do j = l, nm
            i = j + 1
            g = rv1(i)
            y = w(i)
            h = s*g
            g = c*g    
            z = pythag(f, h)
            rv1(j) = z
            c = f/z
            s = h/z
            f = (x*c) + (g*s)
            g =-(x*s) + (g*c)
            h = y*s
            y = y*c
            do jj = 1, N
               x = v(jj,j)
               z = v(jj,i)
               v(jj,j) = (x*c) + (z*s)
               v(jj,i) =-(x*s) + (z*c)
            end do
            z = pythag(f, h)
            w(j) = z
            if (z .ne. ZERO) then
               z = ONE/z
               c = f*z
               s = h*z
            end if
            f = (c*g) + (s*y)
            x =-(s*g) + (c*y)
            do jj = 1, N
               y = a(jj,j)
               z = a(jj,i)
               a(jj,j) = (y*c) + (z*s)
               a(jj,i) =-(y*s) + (z*c)
            end do
         end do
         rv1(l) = ZERO
         rv1(k) = f
         w(k) = x
      end do
3     continue
   end do

end subroutine svdcmp

!============================================================================

subroutine jacobi(a, n, np, d, v, nrot)

   use nfe_lib_mod
   implicit none

   integer :: n, np, nrot, NMAX
   double precision :: a(np,np), d(np), v(np,np)
   integer :: i, ip, iq, j
   double precision :: c, g, h, s, sm, t, tau, theta, tresh, b(n), z(n)

   do ip = 1, n                !initialze to the identity
      do iq = 1, n
         v(ip, iq) = 0.0
      enddo
      v(ip, ip) = 1.0
   enddo
   do ip = 1, n
      b(ip) = a(ip,ip)       ! initialze b and d to the diagnol of a
      d(ip) = b(ip)
      z(ip) = 0.0
   enddo
   nrot = 0
   do i = 1, 50
      sm = 0.0
      do ip = 1, n-1
         do iq = ip+1, n
            sm = sm + abs(a(ip,iq))
         enddo
      enddo
      if (sm .eq. ZERO) return
      if (i .lt. 4) then
         tresh = 0.2*sm/n**2
      else
         tresh = 0.0
      endif
      do ip =1, n-1
         do iq = ip+1, n
            g = 100.0*abs(a(ip,iq))
            if ((i .gt. 4) .and. (abs(d(ip))+ &
            g.eq.abs(d(ip))).and.(abs(d(iq))+ &
            g .eq. abs(d(iq)))) then
               a(ip,iq) = 0.0
            else if (abs(a(ip,iq)) .gt. tresh) then
                 h = d(iq) - d(ip)
                 if (abs(h) + g .eq. abs(h)) then
                    t = a(ip,iq)/h
                 else
                    theta = 0.5*h/a(ip,iq)
                    t = 1.0/(abs(theta) + sqrt(1.0+theta**2))
                    if (theta .lt. ZERO) t = -t
                 endif
                 c = 1.0/sqrt(1.0+t**2)
                 s = t*c
                 tau = s/(1.0+c)
                 h = t*a(ip,iq)
                 z(ip) = z(ip) - h
                 z(iq) = z(iq) + h
                 d(ip) = d(ip) - h
                 d(iq) = d(iq) + h
                 a(ip,iq) = 0.0
                 do j =1, ip-1
                    g = a(j,ip)
                    h = a(j,iq)
                    a(j,ip) = g - s*(h+g*tau)
                    a(j,iq) = h + s*(g-h*tau)
                 enddo
                 do j = ip+1, iq-1
                    g = a(ip,j)
                    h = a(j,iq)
                    a(ip,j) = g-s*(h+g*tau)
                    a(j, iq) = h+s*(g-h*tau)
                 enddo
                 do j = iq+1,n
                    g = a(ip,j)
                    h = a(iq,j)
                    a(ip,j) = g-s*(h+g*tau)
                    a(iq,j) = h+s*(g-h*tau)
                 enddo
                 do j =1, n
                 g = v(j,ip)
                 h = v(j,iq)
                 v(j,ip) = g-s*(h+g*tau)
                 v(j, iq) = h+s*(g-h*tau)
                 enddo
                 nrot = nrot + 1
            endif
         enddo
      enddo
      do ip = 1, n
         b(ip) = b(ip) + z(ip)
         d(ip) = b(ip)
         z(ip) = 0.0
      enddo
   enddo
   return
end subroutine jacobi

!=============================================================================

subroutine eigsrt(d, v, n, np)

   use nfe_lib_mod
   implicit none

   integer :: n, np, i, j, k
   double precision :: d(np), v(np,np), p

   do i = 1, n-1
      k = i
      p = d(i)
      do j = i+1, n
         if (d(j) .ge. p)then
            k = j
            p = d(j)
         endif
      enddo
      if (k .ne. i)then
          d(k) = d(i)
          d(i) = p
          do j = 1, n
             p = v(j,i)
             v(j,i) = v(j,k)
             v(j,k) = p
          enddo
      endif
   enddo
   return
end subroutine eigsrt


!============================================================================
!
! computes sqrt(a**2 + b**2) carefully
!

pure double precision function pythag(a, b)

   use nfe_lib_mod

   implicit none

   double precision, intent(in) :: a, b

   double precision :: absa, absb

   absa = abs(a)
   absb = abs(b)

   if (absa .gt. absb) then
      pythag = absa*sqrt(ONE + (absb/absa)**2)
   else
      if (absb .eq. ZERO) then
         pythag = ZERO
      else
         pythag = absb*sqrt(ONE + (absa/absb)**2)
     endif
   endif

end function pythag
!============================================================================
end module nfe_colvar_mod
