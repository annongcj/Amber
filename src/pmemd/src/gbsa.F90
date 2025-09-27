#include "copyright.i"
! GBSA code

!*******************************************************************************
!
! Module:  gbsa_mod
!
! Description: *Routines to calculate the GBSA contribution
!              *CPU call is located in gb_ene
!              *GPU call is located in gb_force and runs on the CPU for now,
!               a kernel needs to be writen to move the calculation to the GPU
!              *pmemd.cuda only takes nrespa=1, therefore the GBSA routine is 
!               called regardless of the value of need_pot_ene  
!              
!*******************************************************************************

module gbsa_mod

  use gbl_datatypes_mod

  implicit none

! GBSA data structures

  double precision, allocatable, save           :: gbsafrc(:)
  double precision, allocatable, save, private  :: vdwrad(:)
  double precision, allocatable, save, private  :: lcpo_p1(:)
  double precision, allocatable, save, private  :: lcpo_p2(:)
  double precision, allocatable, save, private  :: lcpo_p3(:)
  double precision, allocatable, save, private  :: lcpo_p4(:)

  integer, allocatable, save, private           :: ineighbor(:)
  integer, allocatable, save, private           :: j_vals(:)
  integer, allocatable, save, private           :: k_vals(:)

! End GBSA data structures

contains

!*******************************************************************************
!
! Subroutine:  gbsa_setup
!
! Description:
! GBSA Initialization
! Set up the GB/SA data structures if we're doing a GB/SA calculation
!              
!*******************************************************************************

subroutine gbsa_setup(atm_cnt, num_ints, num_reals)

  use gbl_constants_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer                       :: atm_cnt
  integer                       :: alloc_failed
  integer, allocatable          :: bond_partner(:)
  integer, allocatable          :: total_bond_partner(:)

! Local variables:
  integer                       :: i
  integer                       :: atm_1
  integer                       :: atm_2
  integer                       :: num_bonds
  integer                       :: total_bonds
  integer                       :: atomicnumber
  character(len=4)              :: isymbl ! holder for atom type

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals     

! ROMEE Add this before calling  if (gbsa .eq. 1) then

! First we allocate the main arrays

    allocate(gbsafrc(atm_cnt*3), &
             lcpo_p1(atm_cnt), &
             lcpo_p2(atm_cnt), &
             lcpo_p3(atm_cnt), &
             lcpo_p4(atm_cnt), &
             vdwrad(atm_cnt), &
             j_vals(atm_cnt), &
             k_vals(atm_cnt), &
             ineighbor(atm_cnt * 40), &
             bond_partner(atm_cnt), &
             total_bond_partner(atm_cnt), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    num_reals = num_reals + size(lcpo_p1) + &
                            size(lcpo_p2) + &
                            size(lcpo_p3) + &
                            size(lcpo_p4) + &
                            size(vdwrad)

    ! Initialize
    lcpo_p1(:) = 0.d0
    lcpo_p2(:) = 0.d0
    lcpo_p3(:) = 0.d0
    lcpo_p4(:) = 0.d0
    vdwrad(:)  = 0.d0
    
    ! I don't think we need to include bond_partner here, since we're
    ! deallocating it immediately after we set up other arrays

    num_ints = num_ints + size(j_vals) + &
                          size(k_vals) + &
                          size(ineighbor) + &
                          size(bond_partner) + &
                          size(total_bond_partner)

    ! Initialize
    j_vals(:) = 0
    k_vals(:) = 0
    ineighbor(:) = 0
    bond_partner(:) = 0

    ! For each atom index we find in the BONDS_WITHOUT_HYDROGEN array, add
    ! 1 to that index of bond_partner. That way, we can tell how many atoms
    ! are bonded to each atom

    do i = 0, nbona - 1
      atm_1 = gbl_bond(bonda_idx + i)%atm_i
      atm_2 = gbl_bond(bonda_idx + i)%atm_j
      bond_partner(atm_1) = bond_partner(atm_1) + 1
      bond_partner(atm_2) = bond_partner(atm_2) + 1
    end do
    
    total_bond_partner = bond_partner

    do i = 1, nbonh
      atm_1 = gbl_bond(i)%atm_i
      atm_2 = gbl_bond(i)%atm_j
      total_bond_partner(atm_1) = total_bond_partner(atm_1) + 1
      total_bond_partner(atm_2) = total_bond_partner(atm_2) + 1
    end do

    ! Construct the parameters for the SA calculation; note that the radii
    ! stored in vdwrad are augmented by 1.4 Angstroms

    do i = 1, atm_cnt

      if (loaded_atm_atomicnumber) then
        atomicnumber = atm_atomicnumber(i)
      else
        call get_atomic_number(atm_igraph(i), atm_mass(i), atomicnumber)
      end if
      isymbl = atm_isymbl(i)
      call upper(isymbl)
      num_bonds = bond_partner(i)
      total_bonds = total_bond_partner(i)

      if (atomicnumber .eq. 6) then
        if(total_bonds .eq. 4) then
          if (num_bonds .eq. 1) then
            vdwrad(i)  = 1.70d0 + 1.4d0
            lcpo_p1(i) = 0.77887d0
            lcpo_p2(i) = -0.28063d0
            lcpo_p3(i) = -0.0012968d0
            lcpo_p4(i) = 0.00039328d0
          else if (num_bonds .eq. 2) then
            vdwrad(i)  = 1.70d0 + 1.4d0
            lcpo_p1(i) = 0.56482d0
            lcpo_p2(i) = -0.19608d0
            lcpo_p3(i) = -0.0010219d0
            lcpo_p4(i) = 0.0002658d0
          else if (num_bonds .eq. 3) then
            vdwrad(i)  = 1.70d0 + 1.4d0
            lcpo_p1(i) = 0.23348d0
            lcpo_p2(i) = -0.072627d0
            lcpo_p3(i) = -0.00020079d0
            lcpo_p4(i) = 0.00007967d0
          else if (num_bonds .eq. 4) then
            vdwrad(i)  = 1.70d0 + 1.4d0
            lcpo_p1(i) = 0.00000d0
            lcpo_p2(i) = 0.00000d0
            lcpo_p3(i) = 0.00000d0
            lcpo_p4(i) = 0.00000d0
          else
            write(6,*) 'Unusual nbond for sp3 C:', i, num_bonds, &
              ' Using default carbon LCPO parameters'
            vdwrad(i)  = 1.70d0 + 1.4d0
            lcpo_p1(i) = 0.77887d0
            lcpo_p2(i) = -0.28063d0
            lcpo_p3(i) = -0.0012968d0
            lcpo_p4(i) = 0.00039328d0
          end if
        else
          if (num_bonds .eq. 2) then
            vdwrad(i)  = 1.70d0 + 1.4d0
            lcpo_p1(i) = 0.51245d0
            lcpo_p2(i) = -0.15966d0
            lcpo_p3(i) = -0.00019781d0
            lcpo_p4(i) = 0.00016392d0
          else if (num_bonds .eq. 3) then
            vdwrad(i)  = 1.70d0 + 1.4d0
            lcpo_p1(i) = 0.070344d0
            lcpo_p2(i) = -0.019015d0
            lcpo_p3(i) = -0.000022009d0
            lcpo_p4(i) = 0.000016875d0
          else
            write(6,*) 'Unusual nbond for C :', i, num_bonds, &
              ' Using default carbon LCPO parameters'
            vdwrad(i)  = 1.70d0 + 1.4d0
            lcpo_p1(i) = 0.77887d0
            lcpo_p2(i) = -0.28063d0
            lcpo_p3(i) = -0.0012968d0
            lcpo_p4(i) = 0.00039328d0
          end if
        end if
      else if (atomicnumber .eq. 8) then
        if (isymbl(1:2) .eq. 'O ') then
          vdwrad(i)  = 1.60d0 + 1.4d0
          lcpo_p1(i) = 0.68563d0
          lcpo_p2(i) = -0.1868d0
          lcpo_p3(i) = -0.00135573d0
          lcpo_p4(i) = 0.00023743d0
        else if (isymbl(1:2) .eq. 'O2') then
          vdwrad(i)  = 1.60d0 + 1.4d0
          lcpo_p1(i) = 0.88857d0
          lcpo_p2(i) = -0.33421d0
          lcpo_p3(i) = -0.0018683d0
          lcpo_p4(i) = 0.00049372d0
        else
          if (num_bonds .eq. 1) then
             vdwrad(i)  = 1.60d0 + 1.4d0
             lcpo_p1(i) = 0.77914d0
             lcpo_p2(i) = -0.25262d0
             lcpo_p3(i) = -0.0016056d0
             lcpo_p4(i) = 0.00035071d0
          else if (num_bonds .eq. 2) then
             vdwrad(i)  = 1.60d0 + 1.4d0
             lcpo_p1(i) = 0.49392d0
             lcpo_p2(i) = -0.16038d0
             lcpo_p3(i) = -0.00015512d0
             lcpo_p4(i) = 0.00016453d0
          else
             write(6,*) 'Unusual nbond for O:', i, num_bonds, &
                ' Using default oxygen LCPO parameters'
             vdwrad(i)  = 1.60d0 + 1.4d0
             lcpo_p1(i) = 0.77914d0
             lcpo_p2(i) = -0.25262d0
             lcpo_p3(i) = -0.0016056d0
             lcpo_p4(i) = 0.00035071d0
          end if
        end if
      else if (atomicnumber .eq. 7) then
        if (isymbl(1:2) .eq. 'N3') then
          if (num_bonds .eq. 1) then
             vdwrad(i)  = 1.65d0 + 1.4d0
             lcpo_p1(i) = 0.078602d0
             lcpo_p2(i) = -0.29198d0
             lcpo_p3(i) = -0.0006537d0
             lcpo_p4(i) = 0.00036247d0
          else if (num_bonds .eq. 2) then
             vdwrad(i)  = 1.65d0 + 1.4d0
             lcpo_p1(i) = 0.22599d0
             lcpo_p2(i) = -0.036648d0
             lcpo_p3(i) = -0.0012297d0
             lcpo_p4(i) = 0.000080038d0
          else if (num_bonds .eq. 3) then
             vdwrad(i)  = 1.65d0 + 1.4d0
             lcpo_p1(i) = 0.051481d0
             lcpo_p2(i) = -0.012603d0
             lcpo_p3(i) = -0.00032006d0
             lcpo_p4(i) = 0.000024774d0
          else
             write(6,*) 'Unusual nbond for N3:', i, num_bonds, &
                ' Using default nitrogen LCPO parameters'
             vdwrad(i)  = 1.65d0 + 1.4d0
             lcpo_p1(i) = 0.078602d0
             lcpo_p2(i) = -0.29198d0
             lcpo_p3(i) = -0.0006537d0
             lcpo_p4(i) = 0.00036247d0
          end if
        else
          if (num_bonds .eq. 1) then
             vdwrad(i)  = 1.65d0 + 1.4d0
             lcpo_p1(i) = 0.73511d0
             lcpo_p2(i) = -0.22116d0
             lcpo_p3(i) = -0.00089148d0
             lcpo_p4(i) = 0.0002523d0
          else if (num_bonds .eq. 2) then
             vdwrad(i)  = 1.65d0 + 1.4d0
             lcpo_p1(i) = 0.41102d0
             lcpo_p2(i) = -0.12254d0
             lcpo_p3(i) = -0.000075448d0
             lcpo_p4(i) = 0.00011804d0
          else if (num_bonds .eq. 3) then
             vdwrad(i)  = 1.65d0 + 1.4d0
             lcpo_p1(i) = 0.062577d0
             lcpo_p2(i) = -0.017874d0
             lcpo_p3(i) = -0.00008312d0
             lcpo_p4(i) = 0.000019849d0
          else
             write(6,*) 'Unusual nbond for N:', i, num_bonds, &
                ' Using default nitrogen LCPO parameters'
             vdwrad(i)  = 1.65d0 + 1.4d0
             lcpo_p1(i) = 0.078602d0
             lcpo_p2(i) = -0.29198d0
             lcpo_p3(i) = -0.0006537d0
             lcpo_p4(i) = 0.00036247d0
          end if
        end if
      else if (atomicnumber .eq. 16) then
        if (isymbl(1:2) .eq. 'SH') then
          vdwrad(i)  = 1.90d0 + 1.4d0
          lcpo_p1(i) = 0.7722d0
          lcpo_p2(i) = -0.26393d0
          lcpo_p3(i) = 0.0010629d0
          lcpo_p4(i) = 0.0002179d0
        else
          vdwrad(i)  = 1.90d0 + 1.4d0
          lcpo_p1(i) = 0.54581d0
          lcpo_p2(i) = -0.19477d0
          lcpo_p3(i) = -0.0012873d0
          lcpo_p4(i) = 0.00029247d0
        end if
      else if (atomicnumber .eq. 15) then
        if (num_bonds .eq. 3) then
          vdwrad(i)  = 1.90d0 + 1.4d0
          lcpo_p1(i) = 0.3865d0
          lcpo_p2(i) = -0.18249d0
          lcpo_p3(i) = -0.0036598d0
          lcpo_p4(i) = 0.0004264d0
        else if (num_bonds .eq. 4) then
          vdwrad(i)  = 1.90d0 + 1.4d0
          lcpo_p1(i) = 0.03873d0
          lcpo_p2(i) = -0.0089339d0
          lcpo_p3(i) = 0.0000083582d0
          lcpo_p4(i) = 0.0000030381d0
        else
          write(6,*) 'Unusual nbond for P:', i, num_bonds, &
            ' Using default phosphorus LCPO parameters'
          vdwrad(i)  = 1.90d0 + 1.4d0
          lcpo_p1(i) = 0.3865d0
          lcpo_p2(i) = -0.18249d0
          lcpo_p3(i) = -0.0036598d0
          lcpo_p4(i) = 0.0004264d0
        end if
      else if (isymbl(1:1) .eq. 'Z') then
        vdwrad(i)  = 0.00000d0 + 1.4d0
        lcpo_p1(i) = 0.00000d0
        lcpo_p2(i) = 0.00000d0
        lcpo_p3(i) = 0.00000d0
        lcpo_p4(i) = 0.00000d0
      else if (atomicnumber .eq. 1) then
        vdwrad(i)  = 0.00000d0 + 1.4d0
        lcpo_p1(i) = 0.00000d0
        lcpo_p2(i) = 0.00000d0
        lcpo_p3(i) = 0.00000d0
        lcpo_p4(i) = 0.00000d0
      else if (isymbl(1:2) .eq. 'MG') then
        !  Mg radius = 0.99A: ref. 21 in J. Chem. Phys. 1997, 107, 5422
        !  Mg radius = 1.18A: ref. 30 in J. Chem. Phys. 1997, 107, 5422
        !  Mg radius = 1.45A: Aqvist 1992
        vdwrad(i)  = 1.18d0 + 1.4d0
        !  The following values were taken from O.sp3 with two bonded 
        !  neighbors -> O has the smallest van der Waals radius 
        ! compared to all other elements which had been parametrized
        lcpo_p1(i) = 0.49392d0
        lcpo_p2(i) = -0.16038d0
        lcpo_p3(i) = -0.00015512d0
        lcpo_p4(i) = 0.00016453d0
      else
        ! write( 0,* ) 'bad atom type: ',atype
        ! call mexit( 6,1 )
        vdwrad(i)  = 1.70 + 1.4;
        lcpo_p1(i) = 0.51245;
        lcpo_p2(i) = -0.15966;
        lcpo_p3(i) = -0.00019781;
        lcpo_p4(i) = 0.00016392;
        write(6,'(a,a)') 'Using carbon SA parms for atom type', isymbl
      end if ! (isymbl)

    end do ! (i=1, atm_cnt)

    deallocate(bond_partner)
    deallocate(total_bond_partner)

    num_ints = num_ints - size(bond_partner) - size(total_bond_partner)

!ROMEE Add this to where this is called  end if ! (gbsa .eq. 1)
  return

end subroutine gbsa_setup


!*******************************************************************************
!
! Subroutine:  gbsa_ene
!
! Description: Calculate forces, energies based on Generalized Born SA.
!          Compute a surface-area dependent term if gbsa=1.
! 
!          Do these terms only at "nrespa" multiple-time step intervals;
!          (when igb=2 or 5, one may need to do this at every step)
! 
!
!
!*******************************************************************************

!subroutine gbsa_ene(crd, frc,rgbmaxpsmax2,max_i,frespa, esurf,atm_cnt,jj,r2x )
subroutine gbsa_ene(crd, frc, esurf,atm_cnt,jj,r2x,natbel )

  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use timers_mod

  implicit none
! Formal arguments:

  double precision      :: crd(*)
  double precision, intent(inout) :: frc(*)
  double precision, intent(out)   :: esurf
!  double precision, intent(in)    :: frespa
  !double precision, intent(in)    :: rgbmaxpsmax2                 !
  double precision                :: r2x(*)
  integer               :: natbel

  integer                         :: jj(*)
!  integer, intent(in)             :: max_i
  integer               :: atm_cnt

! Local variables:
  double precision      :: frespa                       !
  double precision      :: rgbmaxpsmax2                 !
  double precision      :: xij, yij, zij                !
  double precision      :: xi, yi, zi                   !
  double precision      :: r2                           !
  double precision      :: dij                          !
  integer               :: max_i
  integer               :: neibr_cnt
  integer               :: icount
  integer               :: i, j, k

  frespa = nrespa
  neibr_cnt = 0 ! we have no neighbors yet for LCPO gbsa
  rgbmaxpsmax2 = (rgbmax + gb_fs_max)**2
  max_i = atm_cnt
  if (natbel .gt. 0) max_i = natbel

!ROMEE check if for MPI this values max_i, numtasks are correct
#ifdef MPI
    do i = mytaskid + 1, max_i, numtasks
#else
    do i = 1, max_i
#endif
      xi = crd(3 * i - 2)
      yi = crd(3 * i - 1)
      zi = crd(3 * i)


      icount = 0
      do j = 1, atm_cnt
        if (i .eq. j) cycle
            
        xij = xi - crd(3 * j - 2)
        yij = yi - crd(3 * j - 1)
        zij = zi - crd(3 * j)
        r2 = xij * xij + yij * yij + zij * zij
        if (r2 .gt. rgbmaxpsmax2) cycle

        ! pairlist contains only atoms within rgbmax + safety margin
            
        icount = icount + 1
        jj(icount) = j
        r2x(icount) = r2
            
      end do
!from above I only need jj and r2x


!      if (gbsa .eq. 1) then
      
        do k = 1, icount

          j = jj(k)
          dij = sqrt(r2x(k))

          if (vdwrad(i) + vdwrad(j) .gt. dij) then
          
            ! Only count LCPO for non-Hydrogens

            if (vdwrad(i) .gt. 2.5 .and. vdwrad(j) .gt. 2.5) then
              neibr_cnt = neibr_cnt + 1
              ineighbor(neibr_cnt) = j
            end if

          end if

        end do

        neibr_cnt = neibr_cnt + 1
        ineighbor(neibr_cnt) = 0

!      end if
         
    end do   ! end loop over atom i = mytaskid + 1, max_i, numtasks
      
!    if (gbsa .eq. 1) 
    call gbsa_lcpo(crd, frc, neibr_cnt, max_i, frespa, esurf)


  return

end subroutine gbsa_ene

!*******************************************************************************
!
! Subroutine: gbsa_lcpo
!
! Description: Calculates the surface area contribution to GB/SA via LCPO approx
!              NOTE: The max_count passed in here was neibr_cnt from gb_ene!
!              neibr_cnt is re-used here for the same type of meaning as above,
!              but it is NOT the same variable as in gb_ene
!
!*******************************************************************************

subroutine gbsa_lcpo(crd, frc, max_count, max_i, frespa, esurf)

  use gbl_constants_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod

  implicit none

! Passed variables

  integer, intent(in)             :: max_count  ! how many neighbors there are
  integer, intent(in)             :: max_i

  double precision, intent(in)    :: frespa
  double precision, intent(in)    :: crd(*)
  double precision, intent(out)   :: esurf
  double precision, intent(inout) :: frc(*)

! Local variables

  integer          :: nei_count  ! neighbor count
  integer          :: i, j, k    ! atom indices
  integer          :: num_j_vals
  integer          :: num_k_vals
  integer          :: icount
  integer          :: count2
  integer          :: count2_fin
  integer          :: j_loop, k_loop

  double precision :: surf_i
  double precision :: xi, yi, zi
  double precision :: xj, yj, zj
  double precision :: xk, yk, zk
  double precision :: xji, yji, zji
  double precision :: xjk, yjk, zjk
  double precision :: aij, ajk
  double precision :: rij, one_rij
  double precision :: rjk, one_rjk
  double precision :: p3p4aij
  double precision :: vdw2dif
  double precision :: lastxj
  double precision :: lastyj
  double precision :: lastzj
  double precision :: ai
  double precision :: totsasa
  
! Running sums (i.e. sum a_ij)
  double precision :: sumaij
  double precision :: sumajk
  double precision :: sumajk2
  double precision :: sumaijajk
  double precision :: sumdaijddijdxi
  double precision :: sumdaijddijdyi
  double precision :: sumdaijddijdzi
  double precision :: sumdaijddijdxiajk
  double precision :: sumdaijddijdyiajk
  double precision :: sumdaijddijdziajk
  double precision :: sumdajkddjkdxj
  double precision :: sumdajkddjkdyj
  double precision :: sumdajkddjkdzj
! Derivatives of pairwise overlaps (i.e., d a_ij / d d_ij)
  double precision :: daijddij
  double precision :: daijddijdxj
  double precision :: daijddijdyj
  double precision :: daijddijdzj
  double precision :: daidxi
  double precision :: daidyi
  double precision :: daidzi
  double precision :: daidxj
  double precision :: daidyj
  double precision :: daidzj
  double precision :: dajkddjk
  double precision :: dajkddjkdxj
  double precision :: dajkddjkdyj
  double precision :: dajkddjkdzj

  nei_count = 1
  totsasa = 0

! Loop over i
#ifdef MPI
  do i = mytaskid + 1, max_i, numtasks
#else
  do i = 1, max_i
#endif
  
    if (ineighbor(nei_count) .eq. 0) then
      nei_count = nei_count + 1
    else
      surf_i = 4.d0 * PI * vdwrad(i) * vdwrad(i)
      sumaij            = 0.d0
      sumajk            = 0.d0
      sumaijajk         = 0.d0
      sumdaijddijdxi    = 0.d0
      sumdaijddijdyi    = 0.d0
      sumdaijddijdzi    = 0.d0
      sumdaijddijdxiajk = 0.d0
      sumdaijddijdyiajk = 0.d0
      sumdaijddijdziajk = 0.d0

      icount = nei_count
      num_j_vals = 0
      
      do
      
        j = ineighbor(nei_count)
        num_j_vals = num_j_vals + 1
        j_vals(num_j_vals) = j
        nei_count = nei_count + 1

        if (ineighbor(nei_count) .ne. 0) cycle
        
        nei_count = nei_count + 1
        exit

      end do

      xi = crd(3 * i - 2)
      yi = crd(3 * i - 1)
      zi = crd(3 * i    )

      ! Loop over j

      do j_loop = 1, num_j_vals

        j = j_vals(j_loop)

        xj = crd(3 * j - 2)
        yj = crd(3 * j - 1)
        zj = crd(3 * j    )
        
        xji = xj - xi
        yji = yj - yi
        zji = zj - zi
        rij = sqrt(xji * xji + yji * yji + zji * zji)
        one_rij = 1 / rij
         
        ! Equation 3 of Weiser, et. al.; J Comput Chem, v20 p217

        aij = 2.d0 * PI * vdwrad(i) * ( vdwrad(i) - rij * 0.5d0 - &
              (vdwrad(i) * vdwrad(i) - vdwrad(j) * vdwrad(j)) * &
              one_rij * 0.5d0)
        
!       write(0, '(a,F16.5)') 'aij is ', aij
        
        ! First derivatives: Appendix A of Weiser, et. al.

         daijddij = PI * vdwrad(i) * (one_rij *one_rij * &
                    (vdwrad(i) * vdwrad(i) - vdwrad(j) * vdwrad(j)) - 1)
         daijddijdxj = daijddij * xji * one_rij
         daijddijdyj = daijddij * yji * one_rij
         daijddijdzj = daijddij * zji * one_rij

         ! Add it up

         sumaij = sumaij + aij

         count2 = icount

         sumajk2 = 0.d0
         sumdajkddjkdxj = 0.d0
         sumdajkddjkdyj = 0.d0
         sumdajkddjkdzj = 0.d0

         p3p4aij = -surften * (lcpo_p3(i) + lcpo_p4(i) * aij) * frespa

         ! Loop over k

         do k_loop = count2, max_count
           if (ineighbor(k_loop) .eq. 0) then
             count2_fin = k_loop - 1
             exit
           end if
         end do

         num_k_vals = 0

         do k_loop = count2, count2_fin
           k = ineighbor(k_loop)

           if (k .ne. j) then
             num_k_vals = num_k_vals + 1
             k_vals(num_k_vals) = k
           end if
         end do

         count2 = icount

         do k_loop = 1, num_k_vals

           k = k_vals(k_loop)
           xk = crd(3 * k - 2)
           yk = crd(3 * k - 1)
           zk = crd(3 * k    )

           xjk = xj - xk
           yjk = yj - yk
           zjk = zj - zk
           rjk = sqrt(xjk * xjk + yjk * yjk + zjk * zjk)
           one_rjk = 1.d0 / rjk

           if (vdwrad(j) + vdwrad(k) .gt. rjk) then
             
             vdw2dif = vdwrad(j) * vdwrad(j) - vdwrad(k) * vdwrad(k)

             ajk = PI * vdwrad(j) * &
                  (2.d0 * vdwrad(j) - rjk - vdw2dif * one_rjk)

             sumajk = sumajk + ajk
             sumajk2 = sumajk2 + ajk

             ! First derivatives

               dajkddjk = PI * vdwrad(j) * one_rjk * &
                     (one_rjk * one_rjk * vdw2dif - 1.d0)
               
               dajkddjkdxj = dajkddjk * xjk
               dajkddjkdyj = dajkddjk * yjk
               dajkddjkdzj = dajkddjk * zjk
               
               frc(3*k-2) = frc(3*k-2) - dajkddjkdxj * p3p4aij
               frc(3*k-1) = frc(3*k-1) - dajkddjkdyj * p3p4aij
               frc(3*k  ) = frc(3*k  ) - dajkddjkdzj * p3p4aij
               
               sumdajkddjkdxj = sumdajkddjkdxj + dajkddjkdxj
               sumdajkddjkdyj = sumdajkddjkdyj + dajkddjkdyj
               sumdajkddjkdzj = sumdajkddjkdzj + dajkddjkdzj
             

           end if ! vdwrad(j) + vdwrad(k) .gt. rjk

         end do ! k_loop = 1, num_k_vals

         ! Finished looping over k

         sumaijajk = sumaijajk+aij*sumajk2

         ! First derivatives 

         sumdaijddijdxi = sumdaijddijdxi - daijddijdxj
         sumdaijddijdyi = sumdaijddijdyi - daijddijdyj
         sumdaijddijdzi = sumdaijddijdzi - daijddijdzj

         sumdaijddijdxiajk = sumdaijddijdxiajk - daijddijdxj * sumajk2
         sumdaijddijdyiajk = sumdaijddijdyiajk - daijddijdyj * sumajk2
         sumdaijddijdziajk = sumdaijddijdziajk - daijddijdzj * sumajk2
         
         lastxj = daijddijdxj * sumajk2 + aij * sumdajkddjkdxj
         lastyj = daijddijdyj * sumajk2 + aij * sumdajkddjkdyj
         lastzj = daijddijdzj * sumajk2 + aij * sumdajkddjkdzj
         
         daidxj = surften * (lcpo_p2(i) * daijddijdxj + lcpo_p3(i) * &
                             sumdajkddjkdxj + lcpo_p4(i) * lastxj)
         daidyj = surften * (lcpo_p2(i) * daijddijdyj + lcpo_p3(i) * &
                             sumdajkddjkdyj + lcpo_p4(i) * lastyj)
         daidzj = surften * (lcpo_p2(i) * daijddijdzj + lcpo_p3(i) * &
                             sumdajkddjkdzj + lcpo_p4(i) * lastzj)

         frc(3*j-2) = frc(3*j-2) - daidxj * frespa
         frc(3*j-1) = frc(3*j-1) - daidyj * frespa
         frc(3*j  ) = frc(3*j  ) - daidzj * frespa

      end do ! (loop_cnt = 1, num_j_vals)

      ai = lcpo_p1(i) * surf_i + lcpo_p2(i) * sumaij + lcpo_p3(i) * sumajk + &
           lcpo_p4(i) * sumaijajk

      daidxi = surften * (lcpo_p2(i) * sumdaijddijdxi + lcpo_p4(i) * &
                          sumdaijddijdxiajk)
      daidyi = surften * (lcpo_p2(i) * sumdaijddijdyi + lcpo_p4(i) * &
                          sumdaijddijdyiajk)
      daidzi = surften * (lcpo_p2(i) * sumdaijddijdzi + lcpo_p4(i) * &
                          sumdaijddijdziajk)
      
      frc(3*i-2) = frc(3*i-2) - daidxi*frespa
      frc(3*i-1) = frc(3*i-1) - daidyi*frespa
      frc(3*i  ) = frc(3*i  ) - daidzi*frespa
      
      totsasa = totsasa + ai
      
    end if ! ineighbor(nei_count) .eq. 0

  end do ! i = 1|mytaskid+1, max_i, [numtasks]

  esurf = surften * totsasa

end subroutine gbsa_lcpo

end module gbsa_mod
! End of Original GBSA code
