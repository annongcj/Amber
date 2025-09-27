#include "copyright.i"
!DG: subroutines are from sander neb code with minor modifications.
module neb_mod

  use file_io_dat_mod
  use parallel_dat_mod

  implicit none

Contains

!-------------------------------------------------------------------------------
#ifdef MPI
  subroutine neb_init()

    use prmtop_dat_mod, only : natom
    use nebread_mod

    implicit none

    integer :: alloc_failed

    neb_nbead = numgroups

    beadid = worldrank*neb_nbead/worldsize + 1

    allocate(neb_nrg_all(neb_nbead), stat=alloc_failed)

    neb_nrg_all(:)=0.d0

    if ( alloc_failed .ne. 0 ) &
      call alloc_error('neb_init()','Error allocating temporary arrays')

    allocate(tangents(3,natom), stat=alloc_failed)

    tangents(:,:)=0.d0

    if ( alloc_failed .ne. 0 ) &
      call alloc_error('neb_init()','Error allocating temporary arrays')

    allocate(springforce(3,natom), stat=alloc_failed)

    springforce(:,:)=0.d0

    if ( alloc_failed .ne. 0 ) &
      call alloc_error('neb_init()','Error allocating temporary arrays')

    allocate(neb_force(3,natom), stat=alloc_failed)

    neb_force(:,:)=0.d0

    if ( alloc_failed .ne. 0 ) &
      call alloc_error('neb_init()','Error allocating temporary arrays')

    allocate(xprev(3,natom), stat=alloc_failed)

    if ( alloc_failed .ne. 0 ) &
      call alloc_error('neb_init()','Error allocating temporary arrays')

    allocate(xnext(3,natom), stat=alloc_failed)

    if ( alloc_failed .ne. 0 ) &
      call alloc_error('neb_init()','Error allocating temporary arrays')

    allocate(fitmask(natom), stat=alloc_failed)

    fitmask(:)=0

    if ( alloc_failed .ne. 0 ) &
      call alloc_error('neb_init()','Error allocating temporary arrays')

    allocate(rmsmask(natom), stat=alloc_failed)

    rmsmask(:)=0

    if ( alloc_failed .ne. 0 ) &
      call alloc_error('neb_init()','Error allocating temporary arrays')


  end subroutine neb_init
#endif
!-------------------------------------------------------------------------------
#ifdef MPI
  subroutine neb_finalize()

    use nebread_mod

    deallocate(neb_nrg_all)
    deallocate(tangents)
    deallocate(springforce)
    deallocate(xprev)
    deallocate(xnext)
    deallocate(neb_force)
    deallocate(fitmask)
    deallocate(rmsmask)
    !deallocate(sortrms)
    !deallocate(sortfit)
    !deallocate(combined_mask)
    !deallocate(partial_mask)
    !deallocate (idx_mask)

  end subroutine neb_finalize
#endif
!-------------------------------------------------------------------------------
#ifdef MPI

#ifdef CUDA
  subroutine transfer_fit_neb_crd(irespa)

    use mdin_ctrl_dat_mod
    use nebread_mod

    implicit none

    integer, intent(in)                 :: irespa

     if ( irespa .eq. 1 .or. mod(irespa,nebfreq)==0 ) then
       call gpu_neb_exchange_crd()
       call gpu_neb_rmsfit()
     endif

  end subroutine transfer_fit_neb_crd
#endif

#ifdef CUDA
  subroutine full_neb_forces(irespa)
#else
  subroutine full_neb_forces(mass, x, f, epot, fitgp, rmsgp, irespa)
#endif

    use mdin_ctrl_dat_mod
    use prmtop_dat_mod
    use nebread_mod
    use pmemd_lib_mod,     only : mexit

    implicit none

    !passed variables
#ifndef CUDA
    double precision                    :: x(:,:)
    double precision                    :: mass(:), f(:,:)
    double precision                    :: epot
    integer                             :: rmsgp(:), fitgp(:)
#endif
    integer, intent(in)                 :: irespa

#ifndef CUDA
    !local variables
    logical                             :: rmsok
    integer                             :: i, j, rep, iatm, jatm
    integer, dimension(MPI_STATUS_SIZE) :: st
    double precision                    :: eplus, eminus, energy, ezero, emax
    double precision                    :: dotproduct, dotproduct2
    double precision                    :: rmsdvalue, spring, spring2
#endif

    if (mod(neb_nbead,2)/=0) then
       write (6,*) "must have even number of beads",neb_nbead
       call mexit(6,1)
    endif

!**********************************************************
! note that beadid runs 1 to neb_nbead and MPI nodes run 0 to #nodes-1
!
! note also that next_node and prev_node are circular,
!  i.e. first and last still have both neighbors defined.
!**********************************************************
#ifdef CUDA
 if(irespa .eq. 1 .or. mod(irespa,nebfreq)==0) then
  call gpu_set_neb_springs()
  call gpu_calculate_neb_frc(neb_force)
 else
  call gpu_calculate_neb_frc_nstep(neb_force)
 endif
#else
 if(irespa .eq. 1 .or. mod(irespa,nebfreq)==0) then
  !these ranks imply we use pmemd_master_comm for the communicator
  next_node = master_rank+1
  if(next_node>=master_size) next_node = 0
  prev_node = master_rank-1
  if(prev_node<0) prev_node = master_size-1

#  ifdef NEBDEBUG
   write(50+worldrank,'(a,i3,a,i3)') "NEBDBG: next_node=",next_node," prev_node=",prev_node
   write(50+worldrank,'(a,i8)') "NEBDBG: last_neb_atom=",last_neb_atom
   write(50+worldrank,'(a,i8)') "NEBDBG: nattgtrms=",nattgtrms
   write(0,'(a,i8)') "NEBDBG: nattgtrms=",nattgtrms
   write(50+worldrank,'(a,i3)') "NEBDBG: neb_nbead=",neb_nbead
#  endif

      !even send coords to rep+1
    if (mod(beadid,2)==0) then  !even, exchange with bead+1

       call mpi_sendrecv(x,3*last_neb_atom,MPI_DOUBLE_PRECISION,next_node,101010, &
                          xnext,3*last_neb_atom,MPI_DOUBLE_PRECISION,next_node,101010, &
                          pmemd_master_comm,st,err_code_mpi)

    else   !odd, exchange with bead-1

       call mpi_sendrecv(x,3*last_neb_atom,MPI_DOUBLE_PRECISION,prev_node,101010, &
                  xprev,3*last_neb_atom,MPI_DOUBLE_PRECISION,prev_node,101010, &
                  pmemd_master_comm,st,err_code_mpi)
    endif

    !even send coords to rep-1
    if (mod(beadid,2)==0) then  !even, exchange with bead-1

       call mpi_sendrecv(x,3*last_neb_atom,MPI_DOUBLE_PRECISION,prev_node,202020, &
                  xprev,3*last_neb_atom,MPI_DOUBLE_PRECISION,prev_node,202020, &
                  pmemd_master_comm,st,err_code_mpi)

    else   !odd, exchange with bead+1 (except bead neb_nbead)

       call mpi_sendrecv(x,3*last_neb_atom,MPI_DOUBLE_PRECISION,next_node,202020, &
                  xnext,3*last_neb_atom,MPI_DOUBLE_PRECISION,next_node,202020, &
                  pmemd_master_comm,st,err_code_mpi)

    endif

#  ifdef NEBDEBUG
   write(50+worldrank,'(a,f12.4)') "NEBDBG: epot=",epot
   write(0,'(a,f12.4)') "NEBDBG: epot=",epot
#  endif

    call mpi_allgather(epot,1,MPI_DOUBLE_PRECISION, &
                       neb_nrg_all,1,MPI_DOUBLE_PRECISION, &
                       pmemd_master_comm,err_code_mpi)

! fit neighbor coords to self

! first argument is refc (which gets fit to the second argument, the actual MD coords)
! refc in this case is the neighbor image
    call rmsfit( xprev, x, mass, fitgp,rmsgp, rmsdvalue, nattgtrms, nattgtfit, rmsok)
    call rmsfit( xnext, x, mass, fitgp,rmsgp, rmsdvalue, nattgtrms, nattgtfit, rmsok)

!ezero is the higher energy of the two end points
#  ifdef NEBDEBUG
   write(50+worldrank,'(a,i3)') "NEBDBG: Getting max (ezero) between 1 and ",neb_nbead
#  endif
   ezero = max( neb_nrg_all(1), neb_nrg_all(neb_nbead) )
#  ifdef NEBDEBUG
   write(50+worldrank,'(a,f12.4)') "NEBDBG: ezero=",ezero
   write(50+worldrank,'(a,f12.4)') "NEBDBG:   neb_nrg_all(1)=",neb_nrg_all(1)
#  endif

   !emax is the highest energy of all the points
   emax = neb_nrg_all(1)
   do rep = 2, neb_nbead
#  ifdef NEBDEBUG
   write(50+worldrank,'(a,i3,a,f12.4)') "NEBDBG:   neb_nrg_all(",rep,")=",neb_nrg_all(rep)
#  endif
   emax = max(emax, neb_nrg_all(rep))
   end do
#  ifdef NEBDEBUG
   write(50+worldrank,'(a,f12.4)') "NEBDBG: emax=",emax
#  endif

   !spring2 is the spring constant of the second point in path
   if (skmin.EQ.skmax) then
     spring2 = skmax
   else if (neb_nrg_all(2)>ezero.AND.emax/=ezero) then
     spring2 = skmax - skmin*((emax-max(neb_nrg_all(1),neb_nrg_all(2)))/  &
         (emax-ezero))
   else
     spring2 = skmax - skmin
   end if

   energy = 0.d0

   rep = beadid

   do i=1,nattgtrms
     do j=1,3
       neb_force(j,i)=0.d0
     end do
   enddo

! CHANGE LATER IF WE DON'T DO (FAKE) MD ON ENDPOINTS!
   if( beadid==1.or.beadid==neb_nbead) then
#  ifdef NEBDEBUG
   write(worldrank+50,'(a)') "-------------------------------------------------"
#  endif
        return
   endif
#  ifdef NEBDEBUG
   write(50+worldrank,'(a,i3,a,i3,a,i3)') "NEBDBG: Rep=",rep," rep+1=",rep+1," rep-1=",rep-1
#  endif
   !calculate spring constant for rep and rep+1

   spring = spring2
   if (skmin.EQ.skmax) then
      spring2 = skmax
   else if (neb_nrg_all(rep+1)>ezero.AND.emax/=ezero) then
      spring2 = skmax - skmin*((emax-max(neb_nrg_all(rep),neb_nrg_all(rep-1)))/ &
                (emax-ezero))
   else
      spring2 = skmax - skmin
   end if

! only fill the tangent elements for atoms in the rmsmask.
! this way the normalization will be ok, and we get to keep the
! tangent pointer and atom coordinate pointers in sync to make the
! code easier to follow as compared to a shorter vector for the
! tangent for just the rmsgroup atoms

! init tangent array since only some elements will be filled

   do i=1,last_neb_atom
     do j=1,3
       tangents(j,i)=0.d0
     end do
   enddo

   if (tmode.EQ.1) then
     !calculate the tangents (for all images except first and last)

     if (neb_nrg_all(rep+1)>neb_nrg_all(rep).AND.neb_nrg_all(rep)>neb_nrg_all(rep-1)) then
        do j=1,nattgtrms
           iatm=rmsgp(j)
           tangents(1,iatm) = xnext(1,iatm) - x(1,iatm)
           tangents(2,iatm) = xnext(2,iatm) - x(2,iatm)
           tangents(3,iatm) = xnext(3,iatm) - x(3,iatm)
        end do
     else if (neb_nrg_all(rep+1)<neb_nrg_all(rep).AND.neb_nrg_all(rep)<neb_nrg_all(rep-1)) then
        do j=1,nattgtrms
           iatm=rmsgp(j)
           tangents(1,iatm) =  x(1,iatm) - xprev(1,iatm)
           tangents(2,iatm) =  x(2,iatm) - xprev(2,iatm)
           tangents(3,iatm) =  x(3,iatm) - xprev(3,iatm)
        end do
     else if (neb_nrg_all(rep+1)>neb_nrg_all(rep-1)) then
        eplus = max( abs(neb_nrg_all(rep+1)-neb_nrg_all(rep)), &
           abs(neb_nrg_all(rep-1)-neb_nrg_all(rep)) )
        eminus = min(abs(neb_nrg_all(rep+1)-neb_nrg_all(rep)), &
           abs(neb_nrg_all(rep-1)-neb_nrg_all(rep)) )

        do j=1,nattgtrms
           iatm=rmsgp(j)
           tangents(1,iatm) = (xnext(1,iatm) - x(1,iatm))*eplus
           tangents(2,iatm) = (xnext(2,iatm) - x(2,iatm))*eplus
           tangents(3,iatm) = (xnext(3,iatm) - x(3,iatm))*eplus
           tangents(1,iatm) = tangents(1,iatm) + (x(1,iatm) - xprev(1,iatm))*eminus
           tangents(2,iatm) = tangents(2,iatm) + (x(2,iatm) - xprev(2,iatm))*eminus
           tangents(3,iatm) = tangents(3,iatm) + (x(3,iatm) - xprev(3,iatm))*eminus
        end do
     else !neb_nrg_all(rep+1)<=neb_nrg_all(rep-1)
        eplus = max(abs(neb_nrg_all(rep+1)-neb_nrg_all(rep)), &
           abs(neb_nrg_all(rep-1)-neb_nrg_all(rep)))
        eminus = min(abs(neb_nrg_all(rep+1)-neb_nrg_all(rep)), &
           abs(neb_nrg_all(rep-1)-neb_nrg_all(rep)))
        do j=1,nattgtrms
           iatm=rmsgp(j)
           tangents(1,iatm) = (xnext(1,iatm) - x(1,iatm))*eminus
           tangents(2,iatm) = (xnext(2,iatm) - x(2,iatm))*eminus
           tangents(3,iatm) = (xnext(3,iatm) - x(3,iatm))*eminus
           tangents(1,iatm) = tangents(1,iatm) + (x(1,iatm) - xprev(1,iatm))*eplus
           tangents(2,iatm) = tangents(2,iatm) + (x(2,iatm) - xprev(2,iatm))*eplus
           tangents(3,iatm) = tangents(3,iatm) + (x(3,iatm) - xprev(3,iatm))*eplus
        end do
     end if
   else !tmode.NE.1 so use strict tangents definition
     do j=1,nattgtrms
        iatm=rmsgp(j)
        tangents(1,iatm) = xnext(1,iatm) - x(1,iatm)
        tangents(2,iatm) = xnext(2,iatm) - x(2,iatm)
        tangents(3,iatm) = xnext(3,iatm) - x(3,iatm)
     end do
   end if

#  ifdef NEBDEBUG
   do i=1,nattgtrms
      iatm=rmsgp(i)
      write(worldrank+50,'(a,f12.4,a,f12.4,a,f12.4)') &
         "NEBDBG: t1=",tangents(1,iatm)," t2=",tangents(2,iatm)," t3=",tangents(3,iatm)
   enddo
#  endif
   call normalize(tangents, last_neb_atom)

   dotproduct = 0.d0
   dotproduct2 = 0.d0

   do i=1,nattgtrms
      iatm=rmsgp(i)

      ! this is "real" force along tangent
      dotproduct = dotproduct + f(1,iatm)*tangents(1,iatm)
      dotproduct = dotproduct + f(2,iatm)*tangents(2,iatm)
      dotproduct = dotproduct + f(3,iatm)*tangents(3,iatm)

      ! now we do the spring forces
      ! note use of two different spring constants

      springforce(1,iatm) = (xnext(1,iatm)-x(1,iatm))*spring2 - (x(1,iatm)-xprev(1,iatm))*spring
      springforce(2,iatm) = (xnext(2,iatm)-x(2,iatm))*spring2 - (x(2,iatm)-xprev(2,iatm))*spring
      springforce(3,iatm) = (xnext(3,iatm)-x(3,iatm))*spring2 - (x(3,iatm)-xprev(3,iatm))*spring

      dotproduct2 = dotproduct2 + springforce(1,iatm)*tangents(1,iatm)
      dotproduct2 = dotproduct2 + springforce(2,iatm)*tangents(2,iatm)
      dotproduct2 = dotproduct2 + springforce(3,iatm)*tangents(3,iatm)

      energy = energy + 0.5*spring2*(xnext(1,iatm)-x(1,iatm))*(xnext(1,iatm)-x(1,iatm))
      energy = energy + 0.5*spring2*(xnext(2,iatm)-x(2,iatm))*(xnext(2,iatm)-x(2,iatm))
      energy = energy + 0.5*spring2*(xnext(3,iatm)-x(3,iatm))*(xnext(3,iatm)-x(3,iatm))

      energy = energy + 0.5*spring*(x(1,iatm)-xprev(1,iatm))*(x(1,iatm)-xprev(1,iatm))
      energy = energy + 0.5*spring*(x(2,iatm)-xprev(2,iatm))*(x(2,iatm)-xprev(2,iatm))
      energy = energy + 0.5*spring*(x(3,iatm)-xprev(3,iatm))*(x(3,iatm)-xprev(3,iatm))
   enddo

   do jatm=1,nattgtrms
      iatm=rmsgp(jatm)

      ! leave modified forces packed in rmsgp, less to communicate
      ! unpack them in force() when they are added to normal forces
      ! note that tangents still use normal coordinate index iatm, while neb_force uses rmsgp pointer

      neb_force(1,jatm) = neb_force(1,jatm) - dotproduct*tangents(1,iatm)
      neb_force(2,jatm) = neb_force(2,jatm) - dotproduct*tangents(2,iatm)
      neb_force(3,jatm) = neb_force(3,jatm) - dotproduct*tangents(3,iatm)

      neb_force(1,jatm) = neb_force(1,jatm) + dotproduct2*tangents(1,iatm)
      neb_force(2,jatm) = neb_force(2,jatm) + dotproduct2*tangents(2,iatm)
      neb_force(3,jatm) = neb_force(3,jatm) + dotproduct2*tangents(3,iatm)

   end do

#  ifdef NEBDEBUG
   write(worldrank+50,'(a)') "-------------------------------------------------"
#  endif

 endif !nebfreq
#endif /* CUDA */
   return
  end subroutine full_neb_forces

#endif /* MPI */
!-------------------------------------------------------------------------------
  subroutine neb_energy_report(filenumber)

    use nebread_mod

    !local variables
    implicit none
    integer :: filenumber,ix
    double precision :: sum

#if defined(CUDA) && defined(MPI)
    call gpu_report_neb_ene(master_size, neb_nrg_all)
#endif

    sum = 0.d0
    do ix = 1, neb_nbead
       sum = sum + neb_nrg_all(ix)
    enddo

    write(filenumber,'(a)') "NEB replicate breakdown:"
    do ix = 1, neb_nbead
       write(filenumber,'(a,i3,a,f13.4)') "Energy for replicate ",ix," = ",neb_nrg_all(ix)
    enddo
    write(filenumber,'(a,f13.4)') "Total Energy of replicates = ",sum

  end subroutine neb_energy_report
!-------------------------------------------------------------------------------
  subroutine rmsfit(xc,x,mass,fitgroup,rmsgroup,rmsdvalue,ntgtatms,nfitatms,rmsok)

    use mdin_ctrl_dat_mod, only : ntr
    use nebread_mod, only : last_neb_atom

    implicit none

    !passed variables:
    double precision :: xc(:,:), x(:,:), mass(:), rmsdvalue
    integer :: fitgroup(:), rmsgroup(:)
    integer :: ntgtatms, nfitatms
    logical :: rmsok

    !local variables:
    double precision :: rotat(3,3), kabs(3,3), e(3), b(3,3)
    double precision :: xcm1, ycm1, zcm1, tmpx, tmpy, tmpz
    double precision :: xcm2, ycm2, zcm2, det, norm
    double precision :: small, totmass
    integer :: i, ii,j, k, ismall

    !following vars needed for lapack diagonalization
    integer :: ldum, info
    integer, parameter :: lwork = 48
    double precision :: kabs2(3,3), work(lwork), dum
    double precision :: wr(3), wi(3), a(3,3)

    rmsok = .true.
    small = 1.d20

    ! zero some variables
    do i=1,3
      e(i) = 0.d0
      do j=1,3
        a(i,j)     = 0.d0
        b(i,j)     = 0.d0
        kabs(i,j)  = 0.d0
        kabs2(i,j) = 0.d0
        rotat(i,j) = 0.d0
      end do
    end do

    ! do the fitting (=overlap) for ntr=0 (no restraints)
    if (ntr == 0) then
      ! STEPS:
      ! 1) center both coord sets on respective CM of fit regions
      ! 2) determine overlap (fit) matrix
      ! 3) rotate reference coords and get rmsdvalue
      ! 4) restore CM of both coords

      ! calculate center of mass of fit region, and center both
      ! molecules (whole molecules not just fit regions)
      xcm1 = 0.d0; ycm1 = 0.d0; zcm1 = 0.d0
      xcm2 = 0.d0; ycm2 = 0.d0; zcm2 = 0.d0

      totmass = 0.d0   ! totmass for fit region

      do ii=1, nfitatms
        i = fitgroup(ii)  ! 'i' is actual atom number
        totmass = totmass + mass(i)
        xcm1 = mass(i) * x(1,i)  + xcm1
        ycm1 = mass(i) * x(2,i)  + ycm1
        zcm1 = mass(i) * x(3,i)  + zcm1
        xcm2 = mass(i) * xc(1,i) + xcm2
        ycm2 = mass(i) * xc(2,i) + ycm2
        zcm2 = mass(i) * xc(3,i) + zcm2
      end do

      xcm1 = xcm1 / totmass
      ycm1 = ycm1 / totmass
      zcm1 = zcm1 / totmass
      xcm2 = xcm2 / totmass
      ycm2 = ycm2 / totmass
      zcm2 = zcm2 / totmass

      ! Move both molecules (all atoms) with respect to CM of fit regions.
      ! Don't modify xc (refc) coordinates, instead do all transformation
      ! that modify xc in a temporary array xctmp().
      do i=1, last_neb_atom !natom DG
        x(1,i) = x(1,i) - xcm1
        x(2,i) = x(2,i) - ycm1
        x(3,i) = x(3,i) - zcm1
        xc(1,i) = xc(1,i) - xcm2
        xc(2,i) = xc(2,i) - ycm2
        xc(3,i) = xc(3,i) - zcm2
      end do

      ! calculate the Kabsch matrix
      do ii=1, nfitatms
        i = fitgroup(ii)
        kabs(1,1) = kabs(1,1) + mass(i) * x(1,i)*xc(1,i)
        kabs(1,2) = kabs(1,2) + mass(i) * x(1,i)*xc(2,i)
        kabs(1,3) = kabs(1,3) + mass(i) * x(1,i)*xc(3,i)
        kabs(2,1) = kabs(2,1) + mass(i) * x(2,i)*xc(1,i)
        kabs(2,2) = kabs(2,2) + mass(i) * x(2,i)*xc(2,i)
        kabs(2,3) = kabs(2,3) + mass(i) * x(2,i)*xc(3,i)
        kabs(3,1) = kabs(3,1) + mass(i) * x(3,i)*xc(1,i)
        kabs(3,2) = kabs(3,2) + mass(i) * x(3,i)*xc(2,i)
        kabs(3,3) = kabs(3,3) + mass(i) * x(3,i)*xc(3,i)
      end do

      !calculate the determinant
      det = kabs(1,1)*kabs(2,2)*kabs(3,3) - kabs(1,1)*kabs(2,3)*kabs(3,2) -  &
            kabs(1,2)*kabs(2,1)*kabs(3,3) + kabs(1,2)*kabs(2,3)*kabs(3,1) +  &
            kabs(1,3)*kabs(2,1)*kabs(3,2) - kabs(1,3)*kabs(2,2)*kabs(3,1)

      ! check that the determinant is not zero
      if (abs(det) < 1.d-5) then
        rmsok = .false.
        write (6,*) "small determinant in rmsfit():", abs(det)
        goto 99
      end if

      ! construct a positive definite matrix by multiplying kabs
      ! matrix by its transpose
      do i=1, 3
        do j=1, 3
          do k=1, 3
            kabs2(i,j) = kabs2(i,j) + kabs(k,i)*kabs(k,j)
          end do
        end do
      end do

      ! Diagonalize kabs2 (matrix RtR in the above reference). On output
      ! kabs2 is destroyed and the eigenvectors are provided in 'a':
      ! a(i=1,3 ; j) is the j-th eigenvector, 'wr' has the eigenvalues.

      ldum = 1
      call dgeev('N','V',3,kabs2,3,wr,wi,dum,ldum,a,3,work,lwork,info)

      if (info /= 0) then
        write(6,'(" Error in diagonalization routine dgeev")')
        rmsok = .false.
        goto 99
      end if

      ! find the smallest eigenvalue
      do i=1, 3
        if (wr(i) < small) then
          ismall = i
          small = wr(i)
        end if
      end do

      ! generate the b vectors
      do j=1, 3
        ! 'wr' very small and negative sometimes -> abs()
        norm = 1.d0/sqrt(abs(wr(j)))
        do i=1, 3
          do k=1, 3
            b(i,j) = b(i,j) + kabs(i,k)*a(k,j)*norm
          end do
        end do
      end do


      ! calculate the rotation matrix rotat

  200 continue

      rotat(:, :) = 0.d0  ! array assignment

      do i=1, 3
        do j=1, 3
          do k=1,3
            rotat(i,j) = rotat(i,j) + b(i,k)*a(j,k)
          end do
        end do
      end do

      !calculate the determinant
      det = rotat(1,1)*rotat(2,2)*rotat(3,3) -  &
            rotat(1,1)*rotat(2,3)*rotat(3,2) -  &
            rotat(1,2)*rotat(2,1)*rotat(3,3) +  &
            rotat(1,2)*rotat(2,3)*rotat(3,1) +  &
            rotat(1,3)*rotat(2,1)*rotat(3,2) -  &
            rotat(1,3)*rotat(2,2)*rotat(3,1)

      !check if the determinant is not zero
      if (abs(det) < 1.d-10) then
        rmsok = .false.
        goto 99
      end if

      ! check if the determinant is negative
      ! If the determinant is negative, invert all elements to
      ! obtain rotation transformation (excluding inversion)
      if (det < 0) then
        do i=1,3
          b(i,ismall) = -b(i,ismall)
        end do
        goto 200
      end if

      ! Rotate reference coordinates: modify xc().
      do i=1, last_neb_atom !natom DG
        tmpx = rotat(1,1)*xc(1,i) + rotat(1,2)*xc(2,i) + rotat(1,3)*xc(3,i)
        tmpy = rotat(2,1)*xc(1,i) + rotat(2,2)*xc(2,i) + rotat(2,3)*xc(3,i)
        tmpz = rotat(3,1)*xc(1,i) + rotat(3,2)*xc(2,i) + rotat(3,3)*xc(3,i)
        xc(1,i) = tmpx
        xc(2,i) = tmpy
        xc(3,i) = tmpz
      end do

      ! calculate current rmsd for rms region (itgtrmsgp)
      totmass = 0.d0      ! totmass for rms region
      rmsdvalue = 0.d0

      do ii=1, ntgtatms
        i = rmsgroup(ii)

        totmass = totmass + mass(i)

        tmpx = xc(1,i) - x(1,i)
        tmpy = xc(2,i) - x(2,i)
        tmpz = xc(3,i) - x(3,i)

        rmsdvalue = rmsdvalue + mass(i)*(tmpx*tmpx + tmpy*tmpy + tmpz*tmpz)
      end do

      ! Now we need to translate both simulation coordinates X and reference
      ! coordinates XC back: remove the centering done previously. We need to
      ! do this to reference coordinates so they are still fit to simulation
      ! coordinates for energy/force calculation (in routines xconst() and
      ! xtgtmd()) - that means ADD the *simulation* coords CM.

      ! for neb we are translating back both the current and the neighbor beads
      ! so that we restore the current image CM (NOT the neighbor's CM).
      ! this should be done the same way for both rmsfit calls done in
      ! full_neb_forces (one for -1 neighbor and one for   1)
      ! as a result CM of both neighbors will be at the original CM of the MD run
      ! this makes no real difference to those images since the coordinates changed
      ! here are copies, and we do not actually modify the coordinates being
      ! used by the neighbor image MD (just the copies that we received from them)

      do i=1, last_neb_atom
        x(1,i) = x(1,i) + xcm1
        x(2,i) = x(2,i) + ycm1
        x(3,i) = x(3,i) + zcm1
        xc(1,i) = xc(1,i) + xcm1  ! yes, it's CM1 (not CM2)
        xc(2,i) = xc(2,i) + ycm1
        xc(3,i) = xc(3,i) + zcm1
      end do

    else
      ! ntr=1 -> don't fit just calculate rmsdvalue. This does not
      ! change X() or XC().

      ! calculate current rmsd for rms region (itgtrmsgp)
      rmsdvalue = 0.d0
      totmass = 0.d0

      do ii=1,ntgtatms
        i = rmsgroup(ii)
        totmass = totmass + mass(i)
        tmpx = xc(1,i) - x(1,i)
        tmpy = xc(2,i) - x(2,i)
        tmpz = xc(3,i) - x(3,i)
        rmsdvalue = rmsdvalue + mass(i)*(tmpx*tmpx + tmpy*tmpy + tmpz*tmpz)
      end do

    end if

    rmsdvalue = sqrt(rmsdvalue/totmass)

    return

  99 rmsdvalue=-99.
    rmsok = .false.
    return

  end subroutine rmsfit
!-------------------------------------------------------------------------------
  subroutine quench(f,v)

    use mdin_ctrl_dat_mod, only: vv, vfac
    use prmtop_dat_mod, only : natom

    implicit none

    !passed variables
    double precision :: f(:,:), v(:,:)
    !f is the force and v is the velocity

    !local variables
    double precision :: force, dotprod

    integer index, indx
    dotprod = 0.d0
    force = 0.d0
    do index=1,natom
      do indx=1,3
        force = force + f(indx,index)**2
        dotprod = dotprod + v(indx,index)*f(indx,index)
      end do
    enddo

    if (force/=0.0d0) then
      force = 1.0d0/sqrt(force)
      dotprod = dotprod*force
    end if

    if (dotprod>0.0d0) then
      v(1,1:natom) = dotprod*f(1,1:natom)*force
      v(2,1:natom) = dotprod*f(2,1:natom)*force
      v(3,1:natom) = dotprod*f(3,1:natom)*force
    else
      v(1,1:natom) = vfac*dotprod*f(1,1:natom)*force
      v(2,1:natom) = vfac*dotprod*f(2,1:natom)*force
      v(3,1:natom) = vfac*dotprod*f(3,1:natom)*force
    end if

  end subroutine quench
!-------------------------------------------------------------------------------
  subroutine normalize(array, ncol)

    implicit none

    !passed variables
    double precision :: array(:,:)
    integer :: ncol

    !local variables
    integer :: i, j
    double precision :: length

    length = 0.d0

    do i = 1, ncol
      do j = 1, 3
        length = length + array(j,i)*array(j,i)
      end do
    end do

    if (length>1.0d-6) then
      length = 1.0d0/sqrt(length)
      array(:,:) = array(:,:)*length
   else
      array(:,:) = 0.d0
   end if

 end subroutine normalize
!-------------------------------------------------------------------------------
  subroutine alloc_error(routine,message)

    use gbl_constants_mod, only : error_hdr, extra_line_hdr
    use pmemd_lib_mod,     only : mexit

    implicit none

    ! passed variables

    character(*), intent(in) :: routine, message

    write(mdout, '(a,a,a)') error_hdr, 'Error in ', routine
    write(mdout, '(a,a)') extra_line_hdr, message

    call mexit(mdout, 1)

  end subroutine alloc_error

end module neb_mod
