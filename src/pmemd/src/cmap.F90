#include "copyright.i"

!*******************************************************************************
!
! Module:  cmap_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module cmap_mod

  use gbl_datatypes_mod

  implicit none

! The following are derived from prmtop dihedral info:

  integer, save                             :: cit_cmap_term_count
  type(cmap_rec), allocatable, save         :: cit_cmap(:)

  !Local
  double precision, parameter :: zero      = 0.0d0
  double precision, parameter :: one       = 1.0d0
  double precision, parameter :: two       = 2.0d0
  double precision, parameter :: three     = 3.0d0
  double precision, parameter :: half      = one/two

  integer  :: gridOrigin=-180 !Where the 2D grid starts in degrees


contains

!*******************************************************************************
!
! Subroutine:  cmap_setup
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine cmap_setup(num_ints, num_reals, use_atm_map)

  use parallel_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals
  integer                       :: use_atm_map(natom)

! Local variables:

  integer               :: alloc_failed
  type(cmap_rec)        :: cmap_copy(cmap_term_count)
  integer               :: atm_i, atm_j, atm_k, atm_l, atm_m, cmap_idx
  integer               :: my_cmap_cnt

! This routine can handle reallocation, and thus can be called multiple
! times.

! Find all diheds for which this process owns either atom:

  my_cmap_cnt = 0

  do cmap_idx = 1, cmap_term_count

    atm_i = gbl_cmap(cmap_idx)%atm_i
    atm_j = gbl_cmap(cmap_idx)%atm_j
    atm_k = gbl_cmap(cmap_idx)%atm_k
    atm_l = gbl_cmap(cmap_idx)%atm_l
    atm_m = gbl_cmap(cmap_idx)%atm_m

#if defined(MPI) && !defined(CUDA)
    if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
      my_cmap_cnt = my_cmap_cnt + 1
      cmap_copy(my_cmap_cnt) = gbl_cmap(cmap_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
      use_atm_map(atm_l) = 1
      use_atm_map(atm_m) = 1
    else if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
      my_cmap_cnt = my_cmap_cnt + 1
      cmap_copy(my_cmap_cnt) = gbl_cmap(cmap_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
      use_atm_map(atm_l) = 1
      use_atm_map(atm_m) = 1
    else if (gbl_atm_owner_map(atm_k) .eq. mytaskid) then
      my_cmap_cnt = my_cmap_cnt + 1
      cmap_copy(my_cmap_cnt) = gbl_cmap(cmap_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
      use_atm_map(atm_l) = 1
      use_atm_map(atm_m) = 1
    else if (gbl_atm_owner_map(atm_l) .eq. mytaskid) then
      my_cmap_cnt = my_cmap_cnt + 1
      cmap_copy(my_cmap_cnt) = gbl_cmap(cmap_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
      use_atm_map(atm_l) = 1
      use_atm_map(atm_m) = 1
    else if (gbl_atm_owner_map(atm_m) .eq. mytaskid) then
      my_cmap_cnt = my_cmap_cnt + 1
      cmap_copy(my_cmap_cnt) = gbl_cmap(cmap_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
      use_atm_map(atm_l) = 1
      use_atm_map(atm_m) = 1
    end if

#else
    my_cmap_cnt = my_cmap_cnt + 1
    cmap_copy(my_cmap_cnt) = gbl_cmap(cmap_idx)
    use_atm_map(atm_i) = 1
    use_atm_map(atm_j) = 1
    use_atm_map(atm_k) = 1
    use_atm_map(atm_l) = 1
    use_atm_map(atm_m) = 1
#endif

  end do

  cit_cmap_term_count = my_cmap_cnt

  if (my_cmap_cnt .gt. 0) then
    if (allocated(cit_cmap)) then
      if (size(cit_cmap) .lt. my_cmap_cnt) then
        num_ints = num_ints - size(cit_cmap) * cmap_rec_ints
        deallocate(cit_cmap)
        allocate(cit_cmap(my_cmap_cnt), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
        num_ints = num_ints + size(cit_cmap) * cmap_rec_ints
      end if
    else
      allocate(cit_cmap(my_cmap_cnt), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(cit_cmap) * cmap_rec_ints
    end if
    cit_cmap(1:my_cmap_cnt) = cmap_copy(1:my_cmap_cnt)
  end if
  
#ifdef CUDA
  call gpu_cmap_setup(cmap_term_count, gbl_cmap, cmap_type_count, gbl_cmap_grid, gbl_cmap_dPhi, gbl_cmap_dPsi, gbl_cmap_dPhi_dPsi)
#endif

  return

end subroutine cmap_setup

!*******************************************************************************
!
! Subroutine:  get_dihed_energy
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine get_cmap_energy(cmap_cnt, cmap, x, frc, ep)

  use gbl_constants_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use charmm_mod, only : AM_VAL_GEOM_torsion
  use ti_mod

  implicit none

! Formal arguments:

  integer,             intent(in)     :: cmap_cnt
  type(cmap_rec),      intent(in)     :: cmap(*)
  double precision,    intent(in)     :: x(3, *)
  double precision,    intent(inout)  :: frc(3, *)
  double precision,    intent(inout)  :: ep

  !LOCAL
  double precision :: phi,psi !Value of repsective dihedral in degrees
  double precision :: dE_dPhi, dE_dPsi !Derivative of energy with
                             !respect to dihedral angles

  double precision :: crd_ijkl(12),dPhi_dijkl(12) !Coordinates and derivatives
  double precision :: cosPhi_ijkl,sinPhi_ijkl     !of first dihedral

  double precision :: crd_jklm(12),dPsi_djklm(12)
  double precision :: cosPsi_jklm,sinPsi_jklm


  double precision :: dE_dijkl(12), dE_djklm(12) !Final derivatives of energy wrt
                                       !to coordinate

  integer          :: jn, n
  integer          :: i, j, k, l, m
  integer          :: cmap_type_index
  double precision :: epw, epl
  double precision :: rad_to_deg_coeff

  epl = 0.d0

  do jn = 1, cmap_cnt

    i = cmap(jn)%atm_i
    j = cmap(jn)%atm_j
    k = cmap(jn)%atm_k
    l = cmap(jn)%atm_l
    m = cmap(jn)%atm_m

    !First dihedral (phi)
    do n = 1,3 !First loop assign x, second assign y, and finally z
      crd_ijkl(n)   = x(n, i)
      crd_ijkl(n+3) = x(n, j)
      crd_ijkl(n+6) = x(n, k)
      crd_ijkl(n+9) = x(n, l)
    enddo

    !Calculate the dihedral angle (phi) and the derivatives of the
    !four coordinates with respect to phi. Remember this subroutine is
    !operating in radians
    call AM_VAL_GEOM_torsion(crd_ijkl,dPhi_dijkl,cosphi_ijkl,sinphi_ijkl)
    phi = sign(acos(cosphi_ijkl),sinphi_ijkl) * RAD_TO_DEG
!DEBUGGING
!        write(6,'(10x,a)'),"Dihedral 1 "
!        write(6,'(10x,4I8)'), i,j,k,l
!
!        write(6,'( 10x,a)'),"coordinates"
!        write(6,'( 10x,a1,3(f10.6,x) )'),"i",crd_ijkl(1:3)
!        write(6,'( 10x,a1,3(f10.6,x) )'),"j",crd_ijkl(4:6)
!        write(6,'( 10x,a1,3(f10.6,x) )'),"k",crd_ijkl(7:9)
!        write(6,'( 10x,a1,3(f10.6,x) )'),"l",crd_ijkl(10:12)
!        write(6,'(a)'),""
!        write(6,'( 10x,a10,f8.3 )'),"Angle: ",phi
!        write(6,'(a)'),""

    !Second dihedral
    !First dihedral (phi)
    do n = 1,3 !First loop assign x, second assign y, and finally z
      crd_jklm(n)   = x(n, j)
      crd_jklm(n+3) = x(n, k)
      crd_jklm(n+6) = x(n, l)
      crd_jklm(n+9) = x(n, m)
    enddo

    !Calculate the dihedral angle (psi) and the derivatives of the
    !four coordinates with respect to psi. Remember this subroutine is
    !operating in radians
    call AM_VAL_GEOM_torsion(crd_jklm,dPsi_djklm,cospsi_jklm,sinpsi_jklm)
    psi = sign(acos(cospsi_jklm),sinpsi_jklm) * RAD_TO_DEG
!DEBUGGING
!        write(6,'(10x,a)'),"Dihedral 2 "
!        write(6,'(10x,4I8)'), j,k,l,m
!
!        write(6,'( 10x,a)'),"coordinates"
!        write(6,'( 10x,a1,3(f10.6,x) )'),"j",crd_jklm(1:3)
!        write(6,'( 10x,a1,3(f10.6,x) )'),"k",crd_jklm(4:6)
!        write(6,'( 10x,a1,3(f10.6,x) )'),"l",crd_jklm(7:9)
!        write(6,'( 10x,a1,3(f10.6,x) )'),"m",crd_jklm(10:12)
!        write(6,'(a)'),""
!        write(6,'( 10x,a10,f8.3 )'),"Angle: ",psi
!        write(6,'(a)'),""

     !Calculate the CMAP correction energy using the associated CMAP
     !parameter. Also calculate the partial derivatives:  dE/dPhi and dE/Psi

     !Remember the sixth element of cmap_index ( cmap_index(6,i) )
     !is the mapping to its corresponding parameter type:
     !
     !       11      13      15      21      23      1
     !        i       j       k       l       m      index to CMAP type

     cmap_type_index = cmap(jn)%parm_idx

     call charmm_calc_cmap_from_phi_psi(phi, psi,                      &
                                        cmap_type_index, &
                                        epw, dE_dPhi, dE_dPsi )

     !In AM_VAL_GEOM_torsion(), we've been working in radians,
     !and conversly charmm_calc_cmap_from_phi_psi(), in degrees.
     !We need to convert over to degrees per interval.

     rad_to_deg_coeff = RAD_TO_DEG * 1/gbl_cmap_grid_step_size(cmap_type_index)


     dPhi_dijkl(1:12) = dPhi_dijkl(1:12) * rad_to_deg_coeff
     dPsi_djklm(1:12) = dPsi_djklm(1:12) * rad_to_deg_coeff

     !Use chain rule to obtain the energy gradient wrt to coordinate

     !    dE_d(ijkl) = dE_dPhi*dPhi_d(ijkl)
     !    dE_d(jklm) = dE_dPsi*dPhi_d(jklm)

     !dE_d(ijkl)
     dE_dijkl(1:12) = dE_dPhi*dPhi_dijkl(1:12)

     !dE_d(jklm)
     dE_djklm(1:12) = dE_dPsi*dPsi_djklm(1:12)

!DEBUGGING
!        write(6,'( 10x,a1,i6)'),"jn",jn
!        write(6,'( 10x,a)'),"forces"
!        write(6,'( 10x,a9,3(f10.6,x) )'),"dE_dijkl",dE_dijkl(1:3)
!        write(6,'( 10x,a9,3(f10.6,x) )'),"dE_dijkl",dE_dijkl(4:6)
!        write(6,'( 10x,a9,3(f10.6,x) )'),"dE_dijkl",dE_dijkl(7:9)
!        write(6,'( 10x,a9,3(f10.6,x) )'),"dE_dijkl",dE_dijkl(10:12)
!
!        write(6,'( 10x,a9,3(f10.6,x) )'),"dE_djklm",dE_djklm(1:3)
!        write(6,'( 10x,a9,3(f10.6,x) )'),"dE_djklm",dE_djklm(4:6)
!        write(6,'( 10x,a9,3(f10.6,x) )'),"dE_djklm",dE_djklm(7:9)
!        write(6,'( 10x,a9,3(f10.6,x) )'),"dE_djklm",dE_djklm(10:12)
!
!        write(6,'(a)'),""

     !Now, update the forces on the five atoms in this CMAP term;
     !remember, it is the negative of the derivative

  ! Now add the force to the main force array.
  ! Remember, force is the negative of the gradient

    if (ti_mode .ne. 0) then
      ti_region = 0
      if (ti_lst(1,i) .eq. 1 .or. ti_lst(1,j) .eq. 1 &
          .or. ti_lst(1,k) .eq. 1 .or. ti_lst(1,l) .eq. 1 &
          .or. ti_lst(1,m) .eq. 1) then
        ti_region = 1
      else if (ti_lst(2,i) .eq. 1 .or. ti_lst(2,j) .eq. 1 &
               .or. ti_lst(2,k) .eq. 1 .or. ti_lst(2,l) .eq. 1 &
               .or. ti_lst(2,m) .eq. 1) then
        ti_region = 2
      end if

      if (ti_region .gt. 0) then
        if (ti_mode .ne. 1) then
          ti_lscale = (ti_sc_lst(i) .lt. 2 .and. &
                       ti_sc_lst(j) .lt. 2 .and. &
                       ti_sc_lst(k) .lt. 2 .and. &
                       ti_sc_lst(l) .lt. 2 .and. &
                       ti_sc_lst(m) .lt. 2)
        end if
#ifdef MPI
        if (gbl_atm_owner_map(i) .eq. mytaskid) then
#endif    
          if (ti_mode .eq. 1 .or. ti_lscale) then !linear scaling        
            call ti_update_ene(epw, si_cmap_ene, ti_region, 1)
            epw = epw * ti_weights(ti_region)
          else
            ti_ene(ti_region,si_cmap_ene) = &
                   ti_ene(ti_region,si_cmap_ene) + epw
            epw = 0.d0    
          end if
#ifdef MPI
        end if
#endif
        if (ti_mode .eq. 1 .or. ti_lscale) then
          dE_dijkl(:) = dE_dijkl(:) * ti_weights(ti_region)
          dE_djklm(:) = dE_djklm(:) * ti_weights(ti_region)
        end if
      end if
    end if

#ifdef MPI
! We use atm_i to determine who sums up the energy...
    if (gbl_atm_owner_map(i) .eq. mytaskid) then
      epl = epl  + epw
      do n = 1,3
        frc(n, i) = frc(n, i) - dE_dijkl(n) 
      enddo
    end if

    if (gbl_atm_owner_map(j) .eq. mytaskid) then
      do n = 1,3
        frc(n, j) = frc(n, j) - dE_dijkl(n+3) - dE_djklm(n)
      enddo
    end if

    if (gbl_atm_owner_map(k) .eq. mytaskid) then
      do n = 1,3
        frc(n, k) = frc(n, k) - dE_dijkl(n+6) - dE_djklm(n+3)
      enddo
    end if

    if (gbl_atm_owner_map(l) .eq. mytaskid) then
      do n = 1,3
        frc(n, l) = frc(n, l) - dE_dijkl(n+9) - dE_djklm(n+6)
      enddo
    end if

    if (gbl_atm_owner_map(m) .eq. mytaskid) then
      do n = 1,3
        frc(n, m) = frc(n, m)                 - dE_djklm(n+9) 
      enddo
    end if

#else
    do n = 1,3
      frc(n, i) = frc(n, i) - dE_dijkl(n)
      frc(n, j) = frc(n, j) - dE_dijkl(n+3) - dE_djklm(n) 
      frc(n, k) = frc(n, k) - dE_dijkl(n+6) - dE_djklm(n+3)
      frc(n, l) = frc(n, l) - dE_dijkl(n+9) - dE_djklm(n+6)
      frc(n, m) = frc(n, m) -                 dE_djklm(n+9)
    end do

    !Update the total cmap energy with calculated energy from the current
    !CMAP term
    epl = epl + epw

#endif

  end do

  ep = epl

  return
end subroutine get_cmap_energy


function cmapGridWrapper(value)
  !Use modular arithmetic to ensure that coordinates on the grid,
  !wrap around at the edges

  integer, intent(in)  :: value
  integer              :: cmapGridWrapper
  integer, parameter   :: resolution = 24 !I should take this from
                                          !cmapParameter%resolution

  cmapGridWrapper = modulo(value-1,resolution)+1

  !Debug
  !write(6,'(a,i8,a,i8)'),"cmapGridWrapper takes",value," returns ",cmapGridWrapper
end function


subroutine generate_cmap_derivatives
!Called *once* after CMAP parameters have been read. This populates
!various partial derivatives in the cmapParameter%{dPsi,dPhi,dPsi_dPhi}
!object. It does this using a cubic spline on the read in CMAP grid.
!
!Later, this information is
!used to evaluate an arbitary phi/phi angle for a given crossterm.

  use prmtop_dat_mod

  implicit none

  integer :: i,row,col,j,k,step_size
  integer :: res,halfRes,twoRes
  double precision, allocatable,dimension(:)   :: tmpy,tmpy2

  write(6,'(a)') ''
  do i=1,cmap_type_count

    !Assign these here to make the code more readable
    !in this section
    res       = gbl_cmap_res(i)
    halfRes   = gbl_cmap_res(i)/2
    twoRes    = gbl_cmap_res(i)*2
    step_size = gbl_cmap_grid_step_size(i)

    allocate( tmpy(twoRes) )
    allocate( tmpy2(twoRes) )

    !1) calculate dE/dPhi
    do row=1,res !Step up one row each cycle,
                 !splining across all columns


      !Fill an *extended* tmp array (tmpy) for with CMAP values
      !for the 1D splining
      !CHARMM does this to (I think) avoid edge issues

      k=1 !interal offset counter for tmp array
      do j=(1-halfRes),(res+halfRes) !-11 --> 36
        tmpy(k) = gbl_cmap_grid(i,row,cmapGridWrapper(j))
        k=k+1
      enddo


      !Calculate spline coeffients (tmpy2) for each of the 1D
      !horizontal rows in the CMAP table
      call generate_cubic_spline(twoRes, step_size, tmpy, tmpy2)
     
      !Calculate %dPhi for using each row
      !of energies and corresponding splines in tmpy and tmpy2
      do j=1,res
        !evaluate_cubic_spline(n,step_size, gridOrigin, y,y2,xin,yout)

        !offset array passed to evaluate_cubic_spline
        call evaluate_cubic_spline(step_size, tmpy, tmpy2, &
                                   j+12, gbl_cmap_dPhi(i,row,j))
      enddo
    enddo !row

    !2) calculate dE/dPsi
    do col=1,res !Step across one column each cycle,
                 !spinling up each column

      k=1
      do j=(1-halfRes),(res+halfRes) !-11 --> 36
        !DEBUG
        !write(6,'(i4,f15.6)'),j,cmap_types(i)%grid(cmapGridWrapper(j),col)
        tmpy(k) = gbl_cmap_grid(i,cmapGridWrapper(j), col)
        k=k+1
      enddo

      !Calculate spline coeffients (tmpy2) for each of the 1D
      !vertical columns in the CMAP table
      call generate_cubic_spline(twoRes, step_size, tmpy, tmpy2)
     
      !Calculate %dPsi for using each column of energies and
      !corresponding splines in tmpy and tmpy2
      do j=1,res
        !offset array passed to evaluate_cubic_spline
        call evaluate_cubic_spline(step_size, tmpy, tmpy2,&
                                   j+12,gbl_cmap_dPsi(i,j,col))
      enddo
    enddo !col


    !3) calculate d^2E/dPhidPsi
    !TODO Interpolate partitial derivative of psi; dE/dPhi ?

    do col=1,res !TODO

      k=1
      do j=(1-halfRes),(res+halfRes) !-11 --> 36
        tmpy(k) = gbl_cmap_dPhi(i,cmapGridWrapper(j),col)
        k=k+1
      enddo

      call generate_cubic_spline(twoRes, step_size, tmpy, tmpy2)
     
      !TODO
      do j=1,res
        call evaluate_cubic_spline(step_size, tmpy, tmpy2, &
                                   j+12,gbl_cmap_dPhi_dPsi(i,j,col))
      enddo
    enddo

    !release the arrays used by the splining code
    deallocate(tmpy)
    deallocate(tmpy2)

 enddo !cmap_type_count

end subroutine generate_cmap_derivatives

subroutine generate_cubic_spline(n, step_size, y, y2)
!Generates the set cubic splining coefficients from a 1D array with
!the assumption that the distance between the points is constant. This
!interpolation method ensures that the derivative of this spline is
!continuous across the boundary of two intervals.

!These coefficients are only calculated *once*, but are later used on
!multiple occassions by another subroutine to interpolate anywhere
!between two points in this 1D array

  implicit none

  integer, intent(in)                :: n !Number of values in below array
  integer, intent(in)                :: step_size !in degrees
!Since n and step_size are known, there is no need to pass an input
!array of x to this subroutine and many simplications can be made
!within this subroutine since we are using a fixed interval
  double precision,  intent(in), dimension(:)  :: y !1D array of y values
  double precision,  intent(out),dimension(:)  :: y2 !calculated coefficients

!Local
  double precision, allocatable,dimension(:)   :: tmp
  double precision                             :: p
  integer                                      :: i

  !tmp array used internally for the decomposition loop
  allocate( tmp(n) )
  
  !Lower and upper boundaries are natural
  y2(1)=zero
  tmp(1)=zero

  do i=2, n-1
    !y2 is used initially as a temp storage array
    p=half*y2(i-1) + two
    !Debug
    y2(i) = -HALF/p
    tmp(i) = ( y(i+1)-two*y(i) + y(i-1) ) / step_size
    tmp(i) = ( three* (tmp(i)/step_size) - HALF *tmp( i-1 ) )/p

  enddo

  !Set the upper boundary
  y2(n) = zero

  do i=n-1,1,-1
    y2(i)=y2(i)*y2(i+1)+tmp(i)
    !write(6, *) "n is :",n
    !write(6, *) "Step size is :",step_size
    !write(6, *) "y is :",y(i)
    !write(6, *) "y2 is :",y2(i)
  enddo


  deallocate(tmp)

end subroutine generate_cubic_spline

subroutine evaluate_cubic_spline(step_size, y, y2, grid_point, dyout)

!subroutine evaluate_cubic_spline(n,step_size, gridOrigin, y,y2,xin,dyout)
!This is reduced version of a typical cubic spline routine
!that one may find in a recipe book of a numerical flavour.

!It is "reduced" since it is used here to obtain smooth derivatives
!ONLY at discrete points on the CMAP grid, hence it never interpolates
!between the points, therefore the coefficient a is always 1.0 and b is
!always 0.0. In addition,  there is never a need to return the
!interpolated value since it will be aways the same as the value passed to it.


  implicit none

  integer, intent(in)                :: step_size !in degrees
                                                   !in degrees
!Since n and step_size are known, there is no need to pass an input
!array of x to this subroutine and many simplications can be made
!within this subroutine since we are using a fixed interval
  double precision, intent(in),dimension(:)    :: y  !1D array of y values
  double precision, intent(in),dimension(:)    :: y2 !calculated coefficients
                                           !from generate_cubic_spline()
!Hack to allow grid points to be passed
  integer, intent(in)                 :: grid_point !WIP
                                             !nearest GRID point
!  _REAL_, intent(out)                :: yout  !cubuic spline
                                              !interpolated value
  double precision, intent(out)                :: dyout !cubuic spline
                                              !interpolated gradient
  double precision                             :: a,b!,c,d
                                        !Cubic spline coefficents
  integer                            :: lo
!  _REAL_                             :: yout


  !Work out nearest complete grid point on the CMAP grid from xin
  !lo =  int( (xin - gridOrigin)/(step_size) ) + 1

  lo = grid_point
  !write(6,'(a,I4)'),"Lo is: ",lo

  !b = ( xin - ( (lo-1)*step_size + gridOrigin)  )/step_size
  !a = 1-b

  !write(6,'(a,f15.6)'),"a is : ",a
  !write(6,'(a,f15.6)'),"b is : ",b


  a = 1.0
  b = 0.0
  !DEBUG

!  yout =   a*y(lo)                                        &
!         + b*y(lo+1)                                      &
!         + (1/6)*(a*a*a-a)*(step_size*step_size)*y2(lo)   &
!         + (1/6)*(b*b*b-b)*(step_size*step_size)*y2(lo+1)


  dyout =  (y(lo+1)-y(lo))/step_size                    &
         - ((3*a*a-1)/6)*step_size*y2(lo)                 &
         + ((3*b*b-1)/6)*step_size*y2(lo+1)

  !DEBUG
  !write(6,'(a,f15.6)'),"y(lo) is :", y(lo)
  !write(6,'(a,f15.6)'),"y(lo+1) is :", y(lo+1)
  !write(6,'(a,f15.6)'),"y2(lo) is :", y2(lo)
  !write(6,'(a,f15.6)'),"y2(lo+1) is :", y2(lo+1)
  !!write(6,'(a,f15.6)'),"yout is :", yout
  !write(6,'(a)'),""
end subroutine evaluate_cubic_spline

subroutine weight_stencil(step_size,         &
                          E_stencil,         &
                          dPhi_stencil,      &
                          dPsi_stencil,      &
                          dPhi_dPsi_stencil, &
                          c )

  implicit none

  integer, intent(in)         :: step_size
  double precision,  dimension(4),&
           intent(in)         :: E_stencil,    dPhi_stencil,&
                                 dPsi_stencil, dPhi_dPsi_stencil
  double precision,  dimension(4,4),&
           intent(out)        :: c !The coefficient array to be be
                                   !used by bicubic_interpolation()


  double precision,  dimension(16)      :: tmp !temp vector built from
                                     !E,dPhi,dPsi,dPhi_dPsi
  integer, dimension(16,16)   :: wT !Weight matrix
!DEBUGGING
!  integer                     :: i
!  integer                     :: j

  double precision,  dimension(16)      :: c_tmp !Vector prior to matrix packing
  integer, dimension (1:2)    :: my_shape = (/ 4, 4 /) !Used for reshape
  integer, dimension (1:2)    :: my_order = (/ 2, 1 /) !Used for reshape
  wt(1 ,1:16) =   (/ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
  wt(2 ,1:16) =   (/ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0/)
  wt(3 ,1:16) =   (/-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0/)
  wt(4 ,1:16) =   (/ 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0/)
  wt(5 ,1:16) =   (/ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
  wt(6 ,1:16) =   (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0/)
  wt(7 ,1:16) =   (/ 0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1/)
  wt(8 ,1:16) =   (/ 0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1/)
  wt(9 ,1:16) =   (/-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
  wt(10,1:16) =   (/ 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0/)
  wt(11,1:16) =   (/ 9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2/)
  wt(12,1:16) =   (/-6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2/)
  wt(13,1:16) =   (/ 2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
  wt(14,1:16) =   (/ 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0/)
  wt(15,1:16) =   (/-6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1/)
  wt(16,1:16) =   (/ 4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1/)

  !Create the vector to be operated on
  tmp(1:4)   = E_stencil(1:4)
  tmp(5:8)   = dPhi_stencil(1:4)      *step_size
  tmp(9:12)  = dPsi_stencil(1:4)      *step_size
  tmp(13:16) = dPhi_dPsi_stencil(1:4) *step_size*step_size

  !Use the matmul intrinsic to operate on the vector
  !http://gcc.gnu.org/onlinedocs/gfortran/MATMUL.html
  c_tmp = matmul(wT,tmp)

  c = reshape(c_tmp, my_shape, order=my_order )

!DEBUGGING
!  write(6,'(a)'),"tmp vector is :"
!  write(6,'(f16.6)'),tmp
!  write(6,'(a)'),""
!  write(6,'(a)'),"c matrix is:"
!  do i=4,1,-1
!     write(6,'(4f16.6)'),c(i,1:4)
!  enddo

end subroutine weight_stencil

!------------------------------------------------------------

!------------------------------------------------------------
subroutine bicubic_interpolation(frac_x, frac_y, c, E, dx, dy)
  implicit none

  double precision, intent(in)         :: frac_x, frac_y
  double precision, dimension(4,4),&
          intent(in)         :: c
  double precision, intent(out)        :: E, dx, dy

  integer                    :: i

  E = 0.0d0
  dx = 0.0d0
  dy = 0.0d0

  do i=4,1,-1

     E =  E*frac_x +( (c(i,4)*frac_y + c(i,3))*frac_y + c(i,2) )*frac_y + c(i,1)
     dx = dx*frac_y +( three*c(4,i)*frac_x + two*c(3,i) )*frac_x + c(2,i)
     dy = dy*frac_x +( three*c(i,4)*frac_y + two*c(i,3) )*frac_y + c(i,2)

  enddo

!DEBUGGING
!  write(6,'(a)'),""
!
!  write(6,'(A,f10.5,f10.5)'),"frac_x,frac_y  ",frac_x,frac_y
!
!  write(6,'(a,f16.6)'),"E is ",E
!  write(6,'(a,f16.6)'),"dx is ",dx
!  write(6,'(a,f16.6)'),"dy is ",dy

end subroutine bicubic_interpolation

subroutine charmm_calc_cmap_from_phi_psi(phi,psi,cmap_idx,E,dPhi,dPsi)
  !Eventually, this will calculate the CMAP energy given
  !psi,phi and the cmap parameter

  use prmtop_dat_mod

  implicit none

  double precision, intent(in)             :: phi,psi
  integer ,intent(in)                      :: cmap_idx
  double precision, intent(out)            :: E !final CMAP energy
  double precision, intent(out)            :: dPhi,dPsi !partitial derivatives


  double precision                         :: phiFrac,psiFrac
  double precision, dimension(2,2)         :: E_stencil,    dPhi_stencil,&
                                    dPsi_stencil, dPhi_dPsi_stencil
  double precision, dimension(4)           :: E_stencil_1D,    dPhi_stencil_1D,&
                                    dPsi_stencil_1D, dPhi_dPsi_stencil_1D
  integer                        :: i,j
  integer                        :: x,y = 0
  double precision,  dimension(4,4)        :: c = 0



  !Work out nearest complete grid point on the CMAP grid from
  !phi and psi and use this to form a 2x2 stencil
  x =  int( (phi - gridOrigin)/(gbl_cmap_grid_step_size(cmap_idx)) ) + 1
  y =  int( (psi - gridOrigin)/(gbl_cmap_grid_step_size(cmap_idx)) ) + 1

  !We are currently at the bottom left of the inner 2x2 elements of the
  !stencil, now move to the bottom left hand corner of the 4x4 stencil
  !x = cmapGridWrapper(x-1)
  !y = cmapGridWrapper(y-1)


  !Work out the fraction of the CMAP grid step that the interpolated
  !point takes up.
  !This will give the remainder part and then it is divided by the step size
  phiFrac = modulo( (phi - gridOrigin), &
               dble(gbl_cmap_grid_step_size(cmap_idx)) ) / gbl_cmap_grid_step_size(cmap_idx)
  psiFrac = modulo( (psi - gridOrigin), &
               dble(gbl_cmap_grid_step_size(cmap_idx)) ) / gbl_cmap_grid_step_size(cmap_idx)

 !2x2 stencil
 !Set the stencils to the local CMAP grid region
 !determined by psi and phi
 do i=1,2
   do j=1,2
    E_stencil(i,j)    = gbl_cmap_grid(cmap_idx ,cmapGridWrapper( y+(i-1) ),&
                                                cmapGridWrapper( x+(j-1) ) )
    dPhi_stencil(i,j) = gbl_cmap_dPhi(cmap_idx,cmapGridWrapper( y+(i-1) ),&
                                                cmapGridWrapper( x+(j-1) ) )
    dPsi_stencil(i,j) = gbl_cmap_dPsi(cmap_idx,cmapGridWrapper( y+(i-1) ),&
                                                cmapGridWrapper( x+(j-1) ) )
    dPhi_dPsi_stencil(i,j) &
                      = gbl_cmap_dPhi_dPsi(cmap_idx, cmapGridWrapper( y+(i-1) ),&
                                                     cmapGridWrapper( x+(j-1) ) )
   enddo
  enddo
  
!DEBUGGING
!    write(6,'(A,f10.3,f10.3)'),"phi, psi ",phi,psi
!    write(6,'(A,f10.5,f10.5)'),"phiFrac,psiFrac  ",phiFrac,psiFrac
!    write(6,'(a,i2)'),"CMAP grid coordinate of psi : ",x
!    write(6,'(a,i2)'),"CMAP grid coordinate of phi : ",y
!
!
!    write(6,'(A)'),""
!    write(6,'(A)'),"2x2 grid stencils:"
!    write(6,'(A)'),"CMAP:"
!
!    write(6,'(4(f9.6,2x))'),E_stencil(2,1:2)
!    write(6,'(4(f9.6,2x))'),E_stencil(1,1:2)
!
!    write(6,'(A)'),""
!
!    write(6,'(A)'),"dPhi:"
!
!    write(6,'(4(f9.6,2x))'),dPhi_stencil(2,1:2)
!    write(6,'(4(f9.6,2x))'),dPhi_stencil(1,1:2)
!
!    write(6,'(A)'),""
!
!    write(6,'(A)'),"dPsi:"
!
!    write(6,'(4(f9.6,2x))'),dPsi_stencil(2,1:2)
!    write(6,'(4(f9.6,2x))'),dPsi_stencil(1,1:2)
!
!    write(6,'(A)'),""
!
!    write(6,'(A)'),"dPhi_dPsi:"
!
!    write(6,'(4(f9.6,2x))'),dPhi_dPsi_stencil(2,1:2)
!    write(6,'(4(f9.6,2x))'),dPhi_dPsi_stencil(1,1:2)
!
!    write(6,'(A)'),""

  !Convert the 2x2 stencils into a 1D array for processing
  !by weight_stencil.
  !The array starts at the bottom left of the 2x2, working
  !around counterclockwise:
  !
  !          4 3
  !          1 2
  !
  !There may be a cleaner/better way of doing this

  E_stencil_1D(1) = E_stencil(1,1)
  E_stencil_1D(2) = E_stencil(1,2)
  E_stencil_1D(3) = E_stencil(2,2)
  E_stencil_1D(4) = E_stencil(2,1)

  !DEBUG
  !write(6,'(f9.6)'),E_stencil_1D

  dPhi_stencil_1D(1) = dPhi_stencil(1,1)
  dPhi_stencil_1D(2) = dPhi_stencil(1,2)
  dPhi_stencil_1D(3) = dPhi_stencil(2,2)
  dPhi_stencil_1D(4) = dPhi_stencil(2,1)

  dPsi_stencil_1D(1) = dPsi_stencil(1,1)
  dPsi_stencil_1D(2) = dPsi_stencil(1,2)
  dPsi_stencil_1D(3) = dPsi_stencil(2,2)
  dPsi_stencil_1D(4) = dPsi_stencil(2,1)

  dPhi_dPsi_stencil_1D(1) = dPhi_dPsi_stencil(1,1)
  dPhi_dPsi_stencil_1D(2) = dPhi_dPsi_stencil(1,2)
  dPhi_dPsi_stencil_1D(3) = dPhi_dPsi_stencil(2,2)
  dPhi_dPsi_stencil_1D(4) = dPhi_dPsi_stencil(2,1)

  !Weight d{Psi,dPhi,dPhidPsi}_stencil
  call weight_stencil(gbl_cmap_grid_step_size(cmap_idx), &
                      E_stencil_1D,              &
                      dPhi_stencil_1D,           &
                      dPsi_stencil_1D,           &
                      dPhi_dPsi_stencil_1D,      &
                      c )

  call bicubic_interpolation(phiFrac, psiFrac, c, E, dphi, dpsi)

end subroutine charmm_calc_cmap_from_phi_psi




end module cmap_mod
