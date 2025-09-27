!*******************************************************************************
!
! Module: external_mod
!
! Description: This module houses the capabilities of calling external libraries
!              to compute energies and forces, instead of using the force field ones.
!
!              Adapted here by Vinicius Wilian D. Cruzeiro and Delaram Ghoreishi
!
!*******************************************************************************

module external_mod

  use external_dat_mod
  use pmemd_lib_mod, only : upcase

  IMPLICIT NONE

  ! The active namelist:

  private       :: extpot

  namelist /extpot/      extprog, json

  contains

  subroutine external_init()

    use prmtop_dat_mod, only  : natom, nres, gbl_labres, gbl_res_atms, atm_mass, atm_igraph, &
                                atm_atomicnumber, loaded_atm_atomicnumber
    use pmemd_lib_mod, only   : get_atomic_number
    use file_io_dat_mod, only : mdout, mdin
    use file_io_mod, only     : nmlsrc
    use pmemd_lib_mod, only   : mexit

    IMPLICIT NONE

    integer :: i, ifind, cnt
    double precision :: coord(3*natom)
    character(len=5), allocatable :: monomers(:)
    integer, allocatable :: nats(:)
    character(len=5), allocatable :: at_name(:)
    integer :: nmon, val
    integer :: atomicnumber(natom)
    logical :: isvalid

    extprog = ''
    json = ''

    ! Read input in namelist format:

    rewind(mdin)                  ! Insurance against maintenance mods.

    call nmlsrc('extpot', mdin, ifind)

    if (ifind .ne. 0) then        ! Namelist found. Read it:
      read(mdin, nml = extpot)
    else                          ! No namelist was present,
      write(mdout, '(a)') 'ERROR: Could not find the namelist extpot in the mdin file'
      call mexit(mdout, 1)
    end if

    if (upcase(extprog) .eq. 'KMMD') then
        call kmmd_init(trim(json)//CHAR(0))
        return
    end if


    do i = 1, natom
      if(loaded_atm_atomicnumber) then
        atomicnumber(i) = atm_atomicnumber(i)
      else
        call get_atomic_number(atm_igraph(i), atm_mass(i), atomicnumber(i))
      end if
    end do

    ! For MBX
    if (upcase(extprog) .eq. 'MBX') then
#ifdef MBX
      nmon = nres

      allocate(nats(nmon),monomers(nmon),at_name(natom))

      cnt = 1
      do i=1,nmon
         isvalid = .True.
         if (gbl_labres(i) == "WAT") then
            val = 3
            monomers(i) = "h2o"
            at_name(cnt+0) = "O"
            if (atomicnumber(cnt+0) .ne. 8) isvalid = .False.
            at_name(cnt+1) = "H"
            if (atomicnumber(cnt+1) .ne. 1) isvalid = .False.
            at_name(cnt+2) = "H"
            if (atomicnumber(cnt+2) .ne. 1) isvalid = .False.
            cnt = cnt + val
         else if (gbl_labres(i) == "N2O") then
            val = 7
            monomers(i) = "n2o5"
            at_name(cnt+0) = "O"
            if (atomicnumber(cnt+0) .ne. 8) isvalid = .False.
            at_name(cnt+1) = "N"
            if (atomicnumber(cnt+1) .ne. 7) isvalid = .False.
            at_name(cnt+2) = "N"
            if (atomicnumber(cnt+2) .ne. 7) isvalid = .False.
            at_name(cnt+3) = "O"
            if (atomicnumber(cnt+3) .ne. 8) isvalid = .False.
            at_name(cnt+4) = "O"
            if (atomicnumber(cnt+4) .ne. 8) isvalid = .False.
            at_name(cnt+5) = "O"
            if (atomicnumber(cnt+5) .ne. 8) isvalid = .False.
            at_name(cnt+6) = "O"
            if (atomicnumber(cnt+6) .ne. 8) isvalid = .False.
            cnt = cnt + val
         else
            write(mdout, '(a,a,a)') 'ERROR: The residue ',gbl_labres(i),' is not recognized by MBX!'
            call mexit(mdout, 1)
         end if
         if (val == gbl_res_atms(i+1)-gbl_res_atms(i)) then
            nats(i) = val
         else
            write(mdout, '(a,a,a)') 'ERROR: The number of atoms in residue ',gbl_labres(i),' does not match the expected by MBX!'
            call mexit(mdout, 1)
         end if
         if (.not. isvalid) then
            write(mdout, '(a,a,a)') 'ERROR: The order or type of the atoms in residue ',gbl_labres(i),&
                                    ' does not match the expected by MBX!'
            call mexit(mdout, 1)
         end if
      end do

      do i =1, natom
        at_name(i)=trim(at_name(i))//CHAR(0)
      end do
      do i=1,nmon
        monomers(i) = trim(monomers(i))//CHAR(0)
      end do

      if (json .ne. '') then
        call initialize_system(coord, nats, at_name, monomers, nmon, trim(json)//CHAR(0))
      else
        call initialize_system(coord, nats, at_name, monomers, nmon)
      end if
#endif
    else
      write(mdout, '(a,a,a)') 'ERROR: External program ',trim(extprog)//CHAR(0),&
             ' is not valid! Please set a valid value in the extprog flag'
      call mexit(mdout, 1)
    end if

  end subroutine

  subroutine gb_external(crd, frc, pot_ene)

    use prmtop_dat_mod,    only :  natom
    use gbl_constants_mod, only :  EV_TO_KCAL, AU_TO_EV, BOHRS_TO_A

    IMPLICIT NONE

    double precision           ::  crd(:,:)
    double precision           ::  frc(:,:)
    double precision           ::  pot_ene
    integer                    ::  i
    double precision           ::  coord(3*natom), grads(3*natom)

    ! For MBX
    if (upcase(extprog) .eq. 'MBX') then
#ifdef MBX
      do i = 1, natom
        coord(3*(i-1)+1) = crd(1,i)
        coord(3*(i-1)+2) = crd(2,i)
        coord(3*(i-1)+3) = crd(3,i)
      end do

      call get_energy_g(coord, natom, pot_ene, grads)

      do i = 1, natom
        frc(1,i) = -grads(3*(i-1)+1)
        frc(2,i) = -grads(3*(i-1)+2)
        frc(3,i) = -grads(3*(i-1)+3)
      end do
#endif
    end if
    
    if (upcase(extprog) .eq. 'KMMD') then
      do i = 1, natom
        coord(3*(i-1)+1) = crd(1,i)
        coord(3*(i-1)+2) = crd(2,i)
        coord(3*(i-1)+3) = crd(3,i)
        grads(3*(i-1)+1) = 0.0
        grads(3*(i-1)+2) = 0.0
        grads(3*(i-1)+3) = 0.0
      end do
      
      call kmmd_frccalc(coord, grads, pot_ene)
      
      do i = 1, natom
        frc(1,i) = frc(1,i) + grads(3*(i-1)+1)
        frc(2,i) = frc(2,i) + grads(3*(i-1)+2)
        frc(3,i) = frc(3,i) + grads(3*(i-1)+3)
      end do
    end if

  end subroutine

  subroutine pme_external(crd, frc, pot_ene)

    use prmtop_dat_mod,    only :  natom
    use pbc_mod,           only :  ucell, recip

    IMPLICIT NONE

    double precision           ::  crd(:,:)
    double precision           ::  frc(:,:)
    double precision           ::  pot_ene
    integer                    ::  i
    double precision           ::  coord(3*natom), grads(3*natom), box(9)

    ! For MBX
    if (upcase(extprog) .eq. 'MBX') then
#ifdef MBX
      do i = 1, natom
        coord(3*(i-1)+1) = crd(1,i)
        coord(3*(i-1)+2) = crd(2,i)
        coord(3*(i-1)+3) = crd(3,i)
      end do

      box(1) = ucell(1,1)
      box(2) = ucell(2,1)
      box(3) = ucell(3,1)
      box(4) = ucell(1,2)
      box(5) = ucell(2,2)
      box(6) = ucell(3,2)
      box(7) = ucell(1,3)
      box(8) = ucell(2,3)
      box(9) = ucell(3,3)

      call get_energy_pbc_g(coord, natom, box, pot_ene, grads)

      do i = 1, natom
        frc(1,i) = -grads(3*(i-1)+1)
        frc(2,i) = -grads(3*(i-1)+2)
        frc(3,i) = -grads(3*(i-1)+3)
      end do
#endif
    end if
    
    if (upcase(extprog) .eq. 'KMMD') then
      do i = 1, natom
        coord(3*(i-1)+1) = crd(1,i)
        coord(3*(i-1)+2) = crd(2,i)
        coord(3*(i-1)+3) = crd(3,i)
        grads(3*(i-1)+1) = 0.0
        grads(3*(i-1)+2) = 0.0
        grads(3*(i-1)+3) = 0.0
      end do

      call kmmd_frccalc(coord, grads, pot_ene)

      do i = 1, natom
        frc(1,i) = frc(1,i) + grads(3*(i-1)+1)
        frc(2,i) = frc(2,i) + grads(3*(i-1)+2)
        frc(3,i) = frc(3,i) + grads(3*(i-1)+3)
      end do
    end if

  end subroutine

  subroutine external_cleanup()

    IMPLICIT NONE

  end subroutine

end module external_mod
