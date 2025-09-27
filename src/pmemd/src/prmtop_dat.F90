#include "copyright.i"

!*******************************************************************************
!
! Module:  prmtop_dat_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module prmtop_dat_mod

  use gbl_datatypes_mod
  use file_io_dat_mod
  use AmberNetcdf_mod, only: NC_checkRestart

  implicit none

! Global parameter definitions:

! The following parameter MAY need to be increased for a very large system:

  integer, parameter    :: max_dihed_dups = 10000

! Global data.  This stuff should be broadcast:

  ! Starting with nttyp, the integer values are derived rather than read.
  ! nttyp is the number of 6-12 vdw parameters.

  integer, parameter    :: prmtop_int_cnt = 38

  integer               natom, ntypes, nbonh, ntheth, nphih, next, nres, &    !7
                        nbona, ntheta, nphia, numbnd, numang, nptra, nphb, &  !14
                        ifbox, ifcap, nspm, numextra, ncopy, nttyp, &         !20
                        bonda_idx, anglea_idx, diheda_idx, C4Pairwise, &    !24 New2022
                        gbl_bond_allocsize, gbl_angle_allocsize, &            !25
                        gbl_dihed_allocsize, next_mult_fac, &                 !27
                        nub, nubtypes, &                                      !29
                        nimphi, nimprtyp, &                                   !31
                        gbl_angle_ub_allocsize, gbl_dihed_imp_allocsize,&     !33
                        cmap_term_count, cmap_type_count, &                   !35
                        gbl_cmap_allocsize, nattgtrms, nattgtfit              !38

  common / prmtop_int / natom, ntypes, nbonh, ntheth, nphih, next, nres, &
                        nbona, ntheta, nphia, numbnd, numang, nptra, nphb, &
                        ifbox, ifcap, nspm, numextra, ncopy, nttyp, &
                        bonda_idx, anglea_idx, diheda_idx, C4Pairwise, & ! New2022
                        gbl_bond_allocsize, gbl_angle_allocsize, &
                        gbl_dihed_allocsize, next_mult_fac, &
                        nub, nubtypes, &
                        nimphi, nimprtyp, &
                        gbl_angle_ub_allocsize, gbl_dihed_imp_allocsize, &
                        cmap_term_count, cmap_type_count, &
                        gbl_cmap_allocsize, nattgtrms, nattgtfit

  save  :: / prmtop_int /

! next_mult_fac is a multiplier that is set to handle differences in gbl_natex,
! etc. storage requirements due to extra points nb14 redefinition.  It is set
! to 2 for non-ep code (handles list doubling), but bumped up to 9 for ep
! code; so far only glycam-associated ep really seem to need this...  This
! is a hack...

! Cap data not currently supported...

! integer, parameter    :: prmtop_dbl_cnt = 4

! double precision      cutcap, xcap, ycap, zcap

! common / prmtop_dbl / cutcap, xcap, ycap, zcap

  ! atm_qterm = the atom partial charge array for amber pme.
  ! atm_mass = the atom mass array.
  ! atm_iac = TBS
  ! typ_ico = TBS
  ! atm_nsp = the atom submolecule index array.
  ! atm_igraph = the atom name array.
  ! gbl_res_atms = The residue atoms index array.
  ! gbl_labres = the residue labels array.

  ! Atom and residue physical consts, names, indices into parameters,
  ! parameters, etc.

  double precision,     allocatable, save       :: atm_qterm(:)
  integer,              allocatable, save       :: atm_atomicnumber(:)
  double precision,     allocatable, save       :: atm_mass(:)
  integer,              allocatable, save       :: atm_iac(:)
  integer,              allocatable, save       :: typ_ico(:)
  double precision,     allocatable, save       :: gbl_cn1(:)
  double precision,     allocatable, save       :: gbl_cn2(:)
  double precision,     allocatable, save       :: gbl_cn6(:)
  integer,              allocatable, save       :: gbl_cn7(:) ! New2022
  double precision,     allocatable, save       :: gbl_cn8(:) ! New2022
  double precision,     allocatable, save       :: gbl_cn114(:)
  double precision,     allocatable, save       :: gbl_cn214(:)
  double precision,     allocatable, save       :: gbl_rk(:)
  double precision,     allocatable, save       :: gbl_req(:)
  double precision,     allocatable, save       :: gbl_rpullk(:)
  double precision,     allocatable, save       :: gbl_rpulleq(:)
  double precision,     allocatable, save       :: gbl_rpressk(:)
  double precision,     allocatable, save       :: gbl_rpresseq(:)
  double precision,     allocatable, save       :: gbl_ub_rk(:)
  double precision,     allocatable, save       :: gbl_ub_r0(:)
  double precision,     allocatable, save       :: gbl_tk(:)
  double precision,     allocatable, save       :: gbl_teq(:)
  double precision,     allocatable, save       :: gbl_asol(:)
  double precision,     allocatable, save       :: gbl_bsol(:)
  double precision,     allocatable, save       :: gbl_pk(:)
  double precision,     allocatable, save       :: gbl_pn(:)
  double precision,     allocatable, save       :: gbl_phase(:)
  double precision,     allocatable, save       :: gbl_imp_pk(:)
  double precision,     allocatable, save       :: gbl_imp_phase(:)
  double precision,     allocatable, save       :: gbl_one_scee(:)
  double precision,     allocatable, save       :: gbl_one_scnb(:)
  double precision,     allocatable, save       :: gbl_cmap_grid(:,:,:)
  double precision,     allocatable, save       :: gbl_cmap_dPhi(:,:,:)
  double precision,     allocatable, save       :: gbl_cmap_dPsi(:,:,:)
  double precision,     allocatable, save       :: gbl_cmap_dPhi_dPsi(:,:,:)

  ! Derived from atom paramaters:

  double precision,     allocatable, save       :: gbl_gamc(:)
  double precision,     allocatable, save       :: gbl_gams(:)
  double precision,     allocatable, save       :: gbl_fmn(:) 
  integer,              allocatable, save       :: gbl_ipn(:)
  integer,              allocatable, save       :: gbl_cmap_res(:)
  integer,              allocatable, save       :: gbl_cmap_grid_step_size(:)


  character(4),         allocatable, save       :: atm_igraph(:)
  integer,              allocatable, save       :: gbl_res_atms(:)
  character(4),         allocatable, save       :: gbl_labres(:)
  integer,              allocatable, save       :: atm_nsp(:)

  ! For Generalized Born:

  double precision,     allocatable, save       :: atm_gb_fs(:)
  double precision,     allocatable, save       :: atm_gb_radii(:)

  ! This is the original atom bond, angle, and dihedral information.
  ! It is used at initialization and when atoms are reassigned to make
  ! the arrays that are actually used in the bonds, angles, and dihedrals
  ! modules.

  type(bond_rec),       allocatable, save       :: gbl_bond(:)
  type(angle_rec),      allocatable, save       :: gbl_angle(:)
  type(angle_ub_rec),   allocatable, save       :: gbl_angle_ub(:)
  type(dihed_rec),      allocatable, save       :: gbl_dihed(:)
  type(dihed_imp_rec),  allocatable, save       :: gbl_dihed_imp(:)
  type(cmap_rec),       allocatable, save       :: gbl_cmap(:)

  ! Module-specific virtual sites frame datatype:
  type prm_ep_frame_rec
    sequence
    integer             :: extra_pnt
    integer             :: ep_cnt
    integer             :: type
    integer             :: parent_atm
    integer             :: frame_atm1
    integer             :: frame_atm2
    integer             :: frame_atm3
    double precision    :: d1
    double precision    :: d2
    double precision    :: d3
  end type prm_ep_frame_rec

  integer, parameter    :: ep_frame_rec_ints = 8 ! don't use for allocation!
  type(prm_ep_frame_rec), allocatable, save     :: prm_ep_frames(:)
  type(prm_ep_frame_rec), parameter :: null_prm_ep_frame_rec = &
                                       prm_ep_frame_rec(2*0,0,0,0,0,0,0,0.0d0,0.0d0,0.0d0)
  integer, save         :: num_custom_ep
  
  ! NOTE! These are deallocated after setup, unless we are using
  !       Generalized Born:

  ! atm_numex = TBS
  ! gbl_natex = TBS
  ! atm_isymbl = the atom symbol array.
  ! atm_itree = TBS

  integer,              allocatable, save       :: atm_numex(:)
  integer,              allocatable, save       :: gbl_natex(:)
  character(4),         allocatable, save       :: atm_isymbl(:)
  character(4),         allocatable, save       :: atm_itree(:)

  ! Main title from prmtop.  No need to broadcast, since only master does i/o.
  
  character(80), save   :: prmtop_ititl

! Old prmtop variables that must be printed out but that are otherwise unused:

  integer, private      :: mbona
  integer, private      :: mphia
  integer, private      :: mtheta
  integer, private      :: natyp
  integer, private      :: nhparm
  integer, private      :: nmxrs
  integer, private      :: nparm

! Have the atomic numbers been loaded in from the prmtop?

  logical,  save        :: loaded_atm_atomicnumber = .false.

! Copies of original values of prmtop variables, the original values of which
! need to be printed out, but which are modified prior to printout:

  integer, private      :: orig_nphih
  integer, private      :: orig_nphia

! Hide internal routines:

  private        alloc_amber_prmtop_mem, duplicate_dihedrals, calc_dihedral_parms

contains

!*******************************************************************************
!
! Subroutine:  init_prmtop_dat
!
! Description:  Read the molecular topology and parameters.
!
!*******************************************************************************

subroutine init_prmtop_dat(num_ints, num_reals, inpcrd_natom)

  use parallel_dat_mod
  use pmemd_lib_mod
  use file_io_mod
  use gbl_constants_mod
  use mdin_ctrl_dat_mod
  use nextprmtop_section_mod
  use charmm_mod, only : charmm_active, check_for_charmm

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals     
  integer, intent(in)           :: inpcrd_natom

! Local variables:

  integer               :: alloc_failed
  integer               :: mol_atm_cnt
  double precision      :: dielc_sqrt
  integer               :: i, j
  integer               :: ic
  character(len=80), allocatable, dimension(:)&
                        :: ffdesc
  character(len=80)     :: line
  integer               :: nlines, iscratch
  character(len=2)      :: word
  logical               :: safe_vs

  ! Stuff used in support of the amber 7 prmtop format:

  integer               :: errcode
  integer               :: atomicnumber
  character(80)         :: fmt
  character(80)         :: fmtin
  character(80)         :: afmt, ifmt, rfmt, mfmt, nfmt ! New2022
  character(80)         :: type

! Other unused prmtop variables:

  integer               :: ifpert, mbper, mdper, mgper, nbper, ndper, ngper
  
! BEGIN DBG
! integer               :: nxt_atm
! END DBG       

! Initialize stuff used in support of the amber 7 prmtop format:

  afmt = '(20A4)'
  ifmt = '(12I6)'
  rfmt = '(5E16.8)'
  mfmt = '(3I8)' ! New2022
  nfmt = '(1E16.8)' ! New2022

! Formatted input:

  call amopen(prmtop, prmtop_name, 'O', 'F', 'R')
  call nxtsec_reset()   ! insurance; the nextprmtop stuff caches the fileptr.

  ! Support both TITLE and CTITLE - CTITLE is used for a chamber
  ! prmtop file. Essentially to prevent earlier versions of the code
  ! from loading such files.

  fmtin = '(a80)'
  type = 'CTITLE'
  call nxtsec(prmtop, mdout, 1, fmtin, type, fmt, errcode)
  if ( errcode /= 0 ) then
    fmtin = '(a80)'
    type = 'TITLE'
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)
  end if

  ! NOTE the hack below (copied from sander) that prevents the hollerith format
  ! returned from the prmtop from screwing things up (the read format is
  ! fmtin, not the returned fmt).

  read(prmtop, fmtin) prmtop_ititl


  fmtin = ifmt
  type = 'POINTERS'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  numextra = 0
  ncopy = 0
  C4Pairwise = 0
  read(prmtop, fmt) natom, ntypes, nbonh, mbona, ntheth, mtheta, nphih, mphia, &
               nhparm, nparm, next, nres, nbona, ntheta, nphia, &
               numbnd, numang, nptra, natyp, nphb, ifpert, nbper, ngper, &
               ndper, mbper, mgper, mdper, ifbox, nmxrs, ifcap, numextra, &
               ncopy, C4Pairwise ! New2022
  
  if (natom .ne. inpcrd_natom) then
    write(mdout, '(a,a)') error_hdr, &
      'natom mismatch in inpcrd/restrt and prmtop files!'
    call mexit(6, 1)
  end if

  if (ifpert .ne. 0) then
    write(mdout, '(a,a)') error_hdr, &
      'PMEMD cannot process prmtop files with perturbation information!'
    call mexit(6, 1)
  end if

  if (nparm .eq. 1) then
    write(mdout, '(a,a)') error_hdr, &
      'PMEMD cannot process prmtop files created by ADDLES with nparm=1!'
    write(mdout, '(a,a)') extra_line_hdr, &
      'Please use sander compiled with -DLES.'
    call mexit(6, 1)
  end if

  if (mbona .ne. nbona .or. mtheta .ne. ntheta .or. mphia .ne. nphia) then
    write(mdout,'(a,a)') error_hdr, &
      'PMEMD no longer allows constraints in prmtop!'
    write(mdout,'(a,a)') extra_line_hdr, &
      'The prmtop must have mbona=nbona, mtheta=ntheta, mphia=nphia.'
    call mexit(6,1)
  end if

  if (numextra .gt. 0 .and. no_intermolecular_bonds .eq. 0) then
    write(mdout,'(a,a)') error_hdr, &
      'PMEMD does not allow intermolecular bonding with extra points!'
    write(mdout,'(a,a)') extra_line_hdr, &
      'You must use no_intermolecular_bonds = 1 (the default) for this prmtop.'
  end if

  ! BUGBUG:
  ! This is a heuristic for possible nb exclusions list growth, a total
  ! hack that may fail in some circumstances.  Need to reengineer a bunch
  ! of code to avoid this in the future...

  if (numextra .eq. 0) then
if(usemidpoint) then
    next_mult_fac = 2
else
    next_mult_fac = 3
endif
  else
    next_mult_fac = 16
  end if

  nttyp = ntypes * (ntypes + 1) / 2

  !    --- Read the force field information from the prmtop if available
  fmtin = '(i,a)'
  type = 'FORCE_FIELD_TYPE'
  call nxtsec(prmtop, mdout, 1, fmtin, type, fmt, errcode)
  if(errcode == 0) then
    ! We found a force field description. Should be 1 or more lines
    ! of text which will be echoed back to the user. In some case
    ! we also take some decisions based on this.
    write (mdout,'(a)') '| Force field information read from topology file: '
    ! format should be (nlines, text) where we read the first lines nlines
    ! to know how many more lines to read - we ignore nlines on subsequent
    ! lines but they still need to be there.
    ! e.g.
    ! 2 PARM99
    ! 2 GAFF
    read(prmtop,fmt) nlines,line
    allocate (ffdesc(nlines), stat = errcode)
    if (errcode /= 0) then
      write(mdout,*) 'Allocation of ffdesc in init_amber_prmtop_dat failed.'
      call mexit(6, 1)
    end if
    ffdesc(1) = line
    write(mdout,'(a,a)') '| ',ffdesc(1)
    do i = 2, nlines
      read(prmtop,fmt) iscratch,ffdesc(i)
      write(mdout,'(a,a)') '| ',ffdesc(i)
    end do

    !Test to see if any special force fields are in use.

    !1) Check for CHARMM
    !Sets charmm_active = .true. if the charmm force field is in use.
if(.not. usemidpoint) &
    call check_for_charmm(nlines,ffdesc)

    !End Test
    deallocate (ffdesc, stat = errcode)
    if (errcode /= 0) then
      write(mdout,*) 'Deallocation of ffdesc in init_amber_prmtop_dat failed.'
      call mexit(6, 1)
    end if

  end if

  if (charmm_active) then
    !Read UB counts so that an allocate can be done
    fmtin = '(2i8)'
    type = 'CHARMM_UREY_BRADLEY_COUNT'
    call nxtsec(prmtop, mdout, 1, fmtin, type, fmt, errcode)

    if (errcode == 0) then
      read(prmtop,fmt) nub, nubtypes
    end if


    !Read CHARMM impropers
    fmtin = '(10i8)'
    type = 'CHARMM_NUM_IMPROPERS'
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)
    read(prmtop,fmt) nimphi

    fmtin = '(i8)'
    type = 'CHARMM_NUM_IMPR_TYPES'
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)
    read(prmtop,fmt) nimprtyp

  end if

    ! if CMAP is active
    fmtin = '(2I8)'
    type = 'CHARMM_CMAP_COUNT'
    call nxtsec(prmtop, mdout, 1, fmtin, type, fmt, errcode)

    if (errcode .ne. 0) then
      type = 'CMAP_COUNT'
      call nxtsec(prmtop, mdout, 1, fmtin, type, fmt, errcode)
    end if

    if (errcode == 0) then
      read(prmtop,fmt) cmap_term_count, cmap_type_count
    end if





! Here we set nscm and ndfmin. We need natom to do this...
  if (ntr.eq.1) &
     nscm = 0


  if (nscm .gt. 0 .and. ntb .eq. 0) then
    ndfmin = 6   ! both translation and rotation com motion removed
    if (natom == 1) ndfmin = 3
    if (natom == 2) ndfmin = 5
  else if (nscm .gt. 0) then
    ndfmin = 3    ! just translation com will be removed
  else
    ndfmin = 0
  end if
  if (ibelly .gt. 0) then   ! No COM Motion Removal, ever.
    nscm = 0
    ndfmin = 0
  end if
  
  if (nscm .le. 0) nscm = 0

  if (ischeme == 1 .and. therm_par .gt. 0.0) ndfmin = 0

  ! For Langevin Dynamics...

  if (gamma_ln .gt. 0.0d0) ndfmin = 0 ! No COM motion removal for LD simulation
 
  call init_amber_prmtop_dat

  ! Reset cap data if necessary:

  if (ivcap .eq. 2) ifcap = 0

  if (ifcap .ne. 0) then
    write(mdout, '(a,a)') error_hdr, &
      'PMEMD cannot process prmtop files with a solvent cap (ifcap != 0)!'
    call mexit(6, 1)
  end if

  ! Check for single H residues, which are archaic, a pain, and no longer
  ! supported.

  do i = 1, nres
    if (gbl_res_atms(i + 1) - gbl_res_atms(i) .eq. 1) then
      if (atm_mass(gbl_res_atms(i)) .lt. 2.d0) then

        ! Check that this is not a virtual site
        safe_vs = .false. 
        do j = 1, num_custom_ep
          if (prm_ep_frames(j)%extra_pnt .eq. gbl_res_atms(i)) then
            safe_vs = .true.
            exit 
          end if
        end do
        if (safe_vs) then
          cycle
        end if
        
        write(mdout, '(a,a)') error_hdr, &
          'PMEMD does not support single H residues!'
        write(mdout, '(a,a,i6,a,i6)') extra_line_hdr, &
          'Atom number ', gbl_res_atms(i), ', Residue number ', i
        call mexit(6, 1)
      end if
    end if
  end do

  ! Check that all atoms belong to a molecule for constant P simulations.
  ! Failure leads to segfaults and chaos.

  if (ntp .gt. 0) then
    
    mol_atm_cnt = 0

    do i = 1, nspm
      mol_atm_cnt = mol_atm_cnt + atm_nsp(i)
    end do

    if (mol_atm_cnt .ne. natom) then
      write(mdout, '(a)') 'Error: Bad topology file. Sum of ATOMS_PER_MOLECULE &
                          &does not equal NATOM.'
      call mexit(6,1)
    end if

  end if

  ! Check Lennard-Jones 12-6-4 functionality
if(.not. usemidpoint) &
  call post_prm_validate_mdin_ctrl_dat

!BEGIN DBG
! nxt_atm = 1
! do i = 1, nspm
!   if (atm_nsp(i) .gt. 3) then
!     write(0,2000)'DBG: molecule, atm_cnt, atoms: ', &
!       i, atm_nsp(i), nxt_atm, '-', nxt_atm + atm_nsp(i) - 1
!   end if
!   nxt_atm = nxt_atm + atm_nsp(i)
! end do

! do i = 1, nres
!   if (gbl_res_atms(i+1) - gbl_res_atms(i) .gt. 3) then
!     write(0,2000)'DBG: res, atm_cnt, atoms: ', &
!       i, gbl_res_atms(i+1) - gbl_res_atms(i), &
!       gbl_res_atms(i),'-',gbl_res_atms(i+1) - 1
!   end if
! end do

!2000 format(a, 3i7, a1, i7)
!END DBG

  return

contains

!----------------------------------------------------------------------------------------------
! count_custom_ep_frames: read a number of custom-specified extra point (virtual site) frames.
!                         The number of frames specified in this section need not be the full
!                         number of extra points in the file.  If it is less, the remaining
!                         extra points will be inferred by the traditional rules.
!
! Arguments:
!     
!----------------------------------------------------------------------------------------------
subroutine count_custom_ep_frames(iunit, iout, max_frames, nfrm)

  ! Formal arguments
  integer, intent(in)  :: iunit, iout, max_frames
  integer, intent(out) :: nfrm

  ! Local variables
  integer              :: nc, tc, is_eof, ierr, nibuff
  character(128)       :: fline
  character            :: fli
  logical              :: line_is_integers, on_number 
  
  ! Read lines sequentially so long as they contain integers
  line_is_integers = .true.
  nibuff = 0
  is_eof = 0
  do while(line_is_integers .and. is_eof .ge. 0)
    fline(:) = '\n'
    on_number = .false.
    read(prmtop, '(A)', iostat=is_eof) fline
    if (is_eof .eq. 0) then
      nc = 0
      tc = 0
      do i = 1, 128
        fli = fline(i:i)
        if (fli .eq. '\n') then
          exit ! Exit the loop 
        else
          nc = nc + 1 
          if ((fli .ge. '0' .and. fli .le. '9') .or. fli .eq. '-') then
            tc = tc + 1
            if (.not. on_number) then
              on_number = .true.
            end if 
          else if (fli .eq. ' ') then
            tc = tc + 1
            if (on_number) then
              on_number = .false.
              nibuff = nibuff + 1
            end if
          end if
        end if
      end do
      if (tc .eq. nc) then
        if (on_number) then
          nibuff = nibuff + 1
        end if
      else
        line_is_integers = .false.
      end if
    end if
  end do
  nfrm = nibuff / 6
  
  ! Rewind the input prmtop
  rewind(iunit)
  
  return
end subroutine count_custom_ep_frames

!*******************************************************************************
!
! Internal Subroutine:  init_amber_prmtop_dat
!
! Description:  Amber forcefield-specific prmtop processing.
!
!*******************************************************************************

subroutine init_amber_prmtop_dat

  use charmm_mod, only          : charmm_active
  use gbl_constants_mod, only   : warn_hdr, info_hdr, extra_line_hdr
#ifndef NOXRAY
  use xray_globals_module, only : xray_active, bulk_solvent_model
#endif
  
  implicit none

  double precision      :: dumd                 ! scratch or trash variable
  integer               :: idum                 ! scratch or trash variable

  !Variables for testing if off diagonal VDW terms have been edited. - NBFIX
  double precision, allocatable :: sigma_tmp(:), epsilon_tmp(:)
  double precision :: c1_tmp, c2_tmp, sig_tmp, eps_tmp, aij_tmp, bij_tmp
  integer :: nbtype_tmp, loop_count, i, j, alloc_failed
  logical :: lj_offdiag
  
  ! BEWARE! The bond, angle, and dihedral arrays get repacked, and the amount
  ! of space used changes.  For this reason, we store the sizes allocated
  ! by the master, and non-master processes just use those sizes.  Also, the
  ! allocated sizes are used to broadcast these arrays.  It is a little
  ! wasteful, but the safest approach for now.  The values for nbonh, nbona,
  ! ntheth, ntheta, nphih, and nphia are all the values read from the prmtop
  ! file; they may later be adjusted in eliminating contstrained bonds or
  ! handling dihedral duplicates.

  bonda_idx = nbonh + 1                   ! Bond array non-hydrogens idx.
  anglea_idx = ntheth + 1                 ! Angle array non-hydrogens idx.
  diheda_idx = nphih + max_dihed_dups + 1 ! Dihedral array non-hydrogens idx.
if(usemidpoint) then
  gbl_bond_allocsize = nbonh + nbona ! nbonh + nbona !the input file for midpoint does
                              !not have nbonh an nbona, thus I artificially put it
  gbl_angle_allocsize = ntheth + ntheta !ntheth + ntheta
  gbl_dihed_allocsize = nphih + nphia + 2 * max_dihed_dups
  gbl_angle_ub_allocsize = nub
  gbl_dihed_imp_allocsize =  nimphi
  gbl_cmap_allocsize =  cmap_term_count
else
  gbl_bond_allocsize = nbonh + nbona
  gbl_angle_allocsize = ntheth + ntheta
  gbl_dihed_allocsize = nphih + nphia + 2 * max_dihed_dups
  gbl_angle_ub_allocsize = nub
  gbl_dihed_imp_allocsize = nimphi
  gbl_cmap_allocsize = cmap_term_count
endif ! usemidpoint

  call alloc_amber_prmtop_mem(num_ints, num_reals)

! Read the various prmtop arrays, making mods as required:

  fmtin = afmt
  type = 'ATOM_NAME'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (atm_igraph(i), i = 1, natom)

  fmtin = rfmt
  type = 'CHARGE'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (atm_qterm(i), i = 1, natom)

! Scale the charges if dielc .gt. 1.0d0:

  if (using_pme_potential .and. dielc .ne. 1.0d0) then
    dielc_sqrt = sqrt(dielc)
    do i = 1, natom
      atm_qterm(i) = atm_qterm(i) / dielc_sqrt
    end do
  end if

  fmtin = rfmt
  type = 'MASS'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (atm_mass(i), i = 1, natom)

  fmtin = ifmt
  type = 'ATOM_TYPE_INDEX'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (atm_iac(i), i = 1, natom)
  
  fmtin = ifmt
  type = 'NUMBER_EXCLUDED_ATOMS'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (atm_numex(i), i = 1, natom)

  fmtin = ifmt
  type = 'NONBONDED_PARM_INDEX'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (typ_ico(i), i = 1, ntypes * ntypes)

  fmtin = afmt
  type = 'RESIDUE_LABEL'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_labres(i), i = 1, nres)

  fmtin = ifmt
  type = 'RESIDUE_POINTER'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_res_atms(i), i = 1, nres)

  gbl_res_atms(nres + 1) = natom + 1

! Read the parameters:

  ! Throw away the solty array.  It is currently not used.

  fmtin = rfmt
  type = 'BOND_FORCE_CONSTANT'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_rk(i),    i = 1, numbnd)

  fmtin = rfmt
  type = 'BOND_EQUIL_VALUE'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_req(i),   i = 1, numbnd)

  ! Read bond adjustments if they exist.  There are two adjustments, each
  ! a piecewise parabolic curve that is flat on one side of its minimum.
  fmtin = rfmt
  type = 'BOND_STIFFNESS_PULL_ADJ'
  call nxtsec(prmtop, mdout, 1, fmtin, type, fmt, errcode)
  if (errcode == 0) then
    write (mdout, '(a)') ''
    write (mdout, '(a)') '| Note: bond stiffnesses will adjust when pulled &
                         &beyond normal lengths'
    read(prmtop, fmt) (gbl_rpullk(i),   i = 1, numbnd)
  else
    gbl_rpullk(:) = 0.0
  endif
  
  fmtin = rfmt
  type = 'BOND_EQUIL_PULL_ADJ'
  call nxtsec(prmtop, mdout, 1, fmtin, type, fmt, errcode)
  if (errcode == 0) then
    read(prmtop, fmt) (gbl_rpulleq(i),   i = 1, numbnd)
  else
    gbl_rpulleq(:) = 100.0
  endif
 
  fmtin = rfmt
  type = 'BOND_STIFFNESS_PRESS_ADJ'
  call nxtsec(prmtop, mdout, 1, fmtin, type, fmt, errcode)
  if (errcode == 0) then
    write (mdout, '(a)') ''
    write (mdout, '(a)') '| Note: bond stiffnesses will adjust when &
                         &compression is deep'
    read(prmtop, fmt) (gbl_rpressk(i),   i = 1, numbnd)
  else
    gbl_rpressk(:) = 0.0    
  endif

  fmtin = rfmt
  type = 'BOND_EQUIL_PRESS_ADJ'
  call nxtsec(prmtop, mdout, 1, fmtin, type, fmt, errcode)
  if (errcode == 0) then
    read(prmtop, fmt) (gbl_rpresseq(i),   i = 1, numbnd)
  else
    gbl_rpresseq(:) = -100.0
  endif
 
  fmtin = rfmt
  type = 'ANGLE_FORCE_CONSTANT'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_tk(i),    i = 1, numang)

  fmtin = rfmt
  type = 'ANGLE_EQUIL_VALUE'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_teq(i),   i = 1, numang)

  fmtin = rfmt
  type = 'DIHEDRAL_FORCE_CONSTANT'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_pk(i),    i = 1, nptra)

  fmtin = rfmt
  type = 'DIHEDRAL_PERIODICITY'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_pn(i),    i = 1, nptra)

  fmtin = rfmt
  type = 'DIHEDRAL_PHASE'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_phase(i), i = 1, nptra)

  ! Read variable SCEE and SCNB values if they exist.
  ! 1) SCEE
  fmtin = rfmt
  type = 'SCEE_SCALE_FACTOR'
  call nxtsec(prmtop,  mdout, 1, fmtin,  type,  fmt,  errcode)
  if (errcode == 0) then
    !We found the SCEE scale factor data so read it in, there should be
    !a value for each dihedral type. Note while there is one for each
    !dihedral type not all are actually used in the calculation since
    !1-4's are only done over unique dihedral terms. If multiple dihedrals
    !for a given set of atom types exist, with different pn's for example
    !then the last value of SCEE/SCNB should be used.

    write (mdout,'(a)') ''
    write (mdout,'(a)') '| Note: 1-4 EEL scale factors are being read &
                        &from the topology file.'

    read(prmtop,fmt) (gbl_one_scee(i), i = 1,nptra)
    ! We used to use vdinv here but this could in principal cause an issue since
    ! nptra includes improper dihedrals which have 0.0d0 in the prmtop file.
    ! call vdinv(nptra,gbl_one_scee,gbl_one_scee)
    do i = 1, nptra
      if (gbl_one_scee(i) /= 0.0d0) then
        gbl_one_scee(i) = 1.0d0/gbl_one_scee(i)
      !else it gets left at 0.0d0.
      end if
    end do
  else
    !We will use default scee of 1.2
    write (mdout,'(a)') ''
    write (mdout,'(a)') '| Note: 1-4 EEL scale factors were NOT found in the &
                         &topology file.'
    write (mdout,'(a)') '|       Using default value of 1.2.'
    do i=1,nptra
      gbl_one_scee(i)=1.0d0/1.2d0
    end do
  end if !if errcode==0

  ! 2) SCNB
  fmtin = rfmt
  type = 'SCNB_SCALE_FACTOR'
  call nxtsec(prmtop,  mdout,  1,fmtin,  type,  fmt,  errcode)
  if (errcode == 0) then
    write (mdout,'(a)') ''
    write (mdout,'(a)') '| Note: 1-4 VDW scale factors are being read from the &
                         &topology file.'

    read(prmtop,fmt) (gbl_one_scnb(i), i = 1,nptra)
    
    ! We used to use vdinv here but this could in principal cause an issue since
    ! nptra includes improper dihedrals which have 0.0d0 in the prmtop file.
    ! call vdinv(nptra,gbl_one_scnb,gbl_one_scnb)
    do i = 1, nptra
      if (gbl_one_scnb(i) /= 0.0d0) then
        gbl_one_scnb(i) = 1.0d0/gbl_one_scnb(i)
      !else it gets left at 0.0d0.
      end if
    end do
  else
    !We will use default scnb of 2.0d0
    write (mdout,'(a)') ''
    write (mdout,'(a)') '| Note: 1-4 VDW scale factors were NOT found in the &
                        &topology file.'
    write (mdout,'(a)') '|       Using default value of 2.0.'
    do i=1,nptra
      gbl_one_scnb(i)=1.0d0/2.0d0
    end do
  end if !if errcode==0
  ! End read variable SCEE and SCNB values.
  
  fmtin = rfmt
  type = 'SOLTY'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (dumd, i = 1, natyp)           ! solty

  fmtin = rfmt
  type = 'LENNARD_JONES_ACOEF'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_cn1(i),   i = 1, nttyp)
  
  fmtin = rfmt
  type = 'LENNARD_JONES_BCOEF'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_cn2(i),   i = 1, nttyp)

if(.not. usemidpoint) &
  gbl_cn6(:)=0.0
  gbl_cn7(:)=0 ! C4PairwiseCUDA
  gbl_cn8(:)=0.0
  
  if (lj1264 .ne. 0) then

    fmtin = rfmt
    type = 'LENNARD_JONES_CCOEF'
    call nxtsec(prmtop, mdout, 1, fmtin, type, fmt, errcode)

    if (errcode .ne. 0) then
      if (lj1264 .eq. -1) then
        ! Off by default if C-coefficients are not present
        lj1264 = 0
      else
        write(mdout, '(2a)') error_hdr, 'Could not find LENNARD_JONES_CCOEF &
                             &in topology file! Recreate topology.'
        call mexit(mdout, 1)
      end if
    else
      read(prmtop, fmt) (gbl_cn6(i),   i = 1, nttyp)
      ! It was either 1 to begin with, or is the default with C-coeffs present
      lj1264 = 1
    end if

  end if

  if (plj1264 .ne. 0) then ! New2022
    ! write(6, '(a)') 'hahaha' ! New2022Debug
    fmtin = mfmt
    type = 'LENNARD_JONES_DCOEF'
    call nxtsec(prmtop, mdout, 1, fmtin, type, fmt, errcode)
    if (errcode .ne. 0) then
      if (plj1264 .eq. -1) then
        ! Off by default if D-coefficients are not present
        plj1264 = 0
      else
        write(mdout, '(2a)') error_hdr, 'Could not find LENNARD_JONES_DCOEF &
                             &in topology file! Recreate topology.'
        call mexit(mdout, 1)
      end if
    else
      read(prmtop, fmt) (gbl_cn7(i),   i = 1, 3*C4Pairwise) ! New2022
      ! It was either 1 to begin with, or is the default with D-coeffs present
      plj1264 = 1
    end if
    if (plj1264 .ne. 0) then
    fmtin = nfmt
    type = 'LENNARD_JONES_DVALUE'
    call nxtsec(prmtop, mdout, 1, fmtin, type, fmt, errcode)
    if (errcode .ne. 0) then
      if (plj1264 .eq. -1) then
        ! Off by default if D-values are not present
        plj1264 = 0
      else
        write(mdout, '(2a)') error_hdr, 'Could not find LENNARD_JONES_DVALUE &
                             &in topology file! Recreate topology.'
        call mexit(mdout, 1)
      end if
    else
      read(prmtop, fmt) (gbl_cn8(i),   i = 1, C4Pairwise) ! New2022
      ! write(6, '(f8.1)') gbl_cn8(1) ! New2022Debug
      ! It was either 1 to begin with, or is the default with D-values present
      plj1264 = 1
      end if
    end if
  end if
! Read the bonding information:

  fmtin = ifmt
  type = 'BONDS_INC_HYDROGEN'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_bond(i)%atm_i, &
                gbl_bond(i)%atm_j, &
                gbl_bond(i)%parm_idx, &
                i = 1, nbonh)

  do i = 1, nbonh
    gbl_bond(i)%atm_i = gbl_bond(i)%atm_i / 3 + 1
    gbl_bond(i)%atm_j = gbl_bond(i)%atm_j / 3 + 1
  end do

  fmtin = ifmt
  type = 'BONDS_WITHOUT_HYDROGEN'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_bond(bonda_idx + i)%atm_i, &
                gbl_bond(bonda_idx + i)%atm_j, &
                gbl_bond(bonda_idx + i)%parm_idx, &
                i = 0, nbona - 1)

  do i = 0, nbona - 1
    gbl_bond(bonda_idx + i)%atm_i = gbl_bond(bonda_idx + i)%atm_i / 3 + 1
    gbl_bond(bonda_idx + i)%atm_j = gbl_bond(bonda_idx + i)%atm_j / 3 + 1
  end do

! Read the angle information:

  fmtin = ifmt
  type = 'ANGLES_INC_HYDROGEN'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_angle(i)%atm_i, &
                gbl_angle(i)%atm_j, &
                gbl_angle(i)%atm_k, &
                gbl_angle(i)%parm_idx, &
                i = 1, ntheth)

  do i = 1, ntheth
    gbl_angle(i)%atm_i = gbl_angle(i)%atm_i / 3 + 1
    gbl_angle(i)%atm_j = gbl_angle(i)%atm_j / 3 + 1
    gbl_angle(i)%atm_k = gbl_angle(i)%atm_k / 3 + 1
  end do

  fmtin = ifmt
  type = 'ANGLES_WITHOUT_HYDROGEN'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_angle(anglea_idx + i)%atm_i, &
                gbl_angle(anglea_idx + i)%atm_j, &
                gbl_angle(anglea_idx + i)%atm_k, &
                gbl_angle(anglea_idx + i)%parm_idx, &
                i = 0, ntheta - 1)

  do i = 0, ntheta - 1
    gbl_angle(anglea_idx + i)%atm_i = gbl_angle(anglea_idx + i)%atm_i / 3 + 1
    gbl_angle(anglea_idx + i)%atm_j = gbl_angle(anglea_idx + i)%atm_j / 3 + 1
    gbl_angle(anglea_idx + i)%atm_k = gbl_angle(anglea_idx + i)%atm_k / 3 + 1
  end do

  if (charmm_active) then
    !Specific 1-4 CHARMM scaling
    fmtin = rfmt 
    type = 'LENNARD_JONES_14_ACOEF'
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

    read(prmtop, fmt) (gbl_cn114(i),   i = 1, nttyp)
  
    fmtin = rfmt
    type = 'LENNARD_JONES_14_BCOEF' 
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

    read(prmtop, fmt) (gbl_cn214(i),   i = 1, nttyp)

    fmtin = ifmt
    type = 'CHARMM_UREY_BRADLEY'
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)
    read(prmtop, fmt) (gbl_angle_ub(i)%atm_i,&
                       gbl_angle_ub(i)%atm_k,&
                       gbl_angle_ub(i)%parm_idx,&
                  i = 1, nub)

    fmtin = rfmt
    type = 'CHARMM_UREY_BRADLEY_FORCE_CONSTANT'
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

    read(prmtop, fmt) (gbl_ub_rk(i),    i = 1, nubtypes)

    fmtin = rfmt
    type = 'CHARMM_UREY_BRADLEY_EQUIL_VALUE'
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

    read(prmtop, fmt) (gbl_ub_r0(i),    i = 1, nubtypes)


    ! Read CHARMM impropers
    fmtin = ifmt
    type = 'CHARMM_IMPROPERS'
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)
    read(prmtop, fmt) (gbl_dihed_imp(i)%atm_i,&
                       gbl_dihed_imp(i)%atm_j,&
                       gbl_dihed_imp(i)%atm_k,&
                       gbl_dihed_imp(i)%atm_l,&
                       gbl_dihed_imp(i)%parm_idx,&
                  i = 1, nimphi)

    fmtin = rfmt
    type = 'CHARMM_IMPROPER_FORCE_CONSTANT'
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

    read(prmtop, fmt) (gbl_imp_pk(i),    i = 1, nimprtyp)

    fmtin = rfmt
    type = 'CHARMM_IMPROPER_PHASE'
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

    read(prmtop, fmt) (gbl_imp_phase(i),    i = 1, nimprtyp)
  
  end if

    !CMAP

  if ( cmap_term_count > 0 ) then
    fmtin = '(9I8)'
    type = 'CHARMM_CMAP_INDEX'
    call nxtsec(prmtop, mdout, 1, fmtin, type, fmt, errcode)
    if (errcode .ne. 0) then
      type = 'CMAP_INDEX'
      call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)
    end if
    read(prmtop, fmt) (gbl_cmap(i)%atm_i,&
                       gbl_cmap(i)%atm_j,&
                       gbl_cmap(i)%atm_k,&
                       gbl_cmap(i)%atm_l,&
                       gbl_cmap(i)%atm_m,&
                       gbl_cmap(i)%parm_idx,&
                    i=1, cmap_term_count)

    fmtin = '(20I4)'
    type = 'CHARMM_CMAP_RESOLUTION'
    call nxtsec(prmtop, mdout, 1, fmtin, type, fmt, errcode)
    if (errcode .ne. 0) then
      type = 'CMAP_RESOLUTION'
      call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)
    end if
    read(prmtop, fmt) (gbl_cmap_res(i), i=1, cmap_type_count)

    do i=1,cmap_type_count
      write(word,'(i2.2)')i
      fmtin = '(8(F9.5))'
      type = "CHARMM_CMAP_PARAMETER_" // word
      call nxtsec(prmtop, mdout, 1, fmtin, type, fmt, errcode)
      if (errcode .ne. 0) then
        type = "CMAP_PARAMETER_" // word
        call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)
      end if

      read(prmtop, fmt) gbl_cmap_grid(i, 1:gbl_cmap_res(i), &
                                         1:gbl_cmap_res(i) )

      !Also take advantage of this loop to set the step size
      gbl_cmap_grid_step_size(i) = 360/gbl_cmap_res(i)
    enddo

  end if ! cmap_term_count > 0


! Read the dihedral information:

  fmtin = ifmt
  type = 'DIHEDRALS_INC_HYDROGEN'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_dihed(i)%atm_i, &
                gbl_dihed(i)%atm_j, &
                gbl_dihed(i)%atm_k, &
                gbl_dihed(i)%atm_l, &
                gbl_dihed(i)%parm_idx, &
                i = 1, nphih)

  do i = 1, nphih
    gbl_dihed(i)%atm_i = gbl_dihed(i)%atm_i / 3 + 1
    gbl_dihed(i)%atm_j = gbl_dihed(i)%atm_j / 3 + 1
    if (gbl_dihed(i)%atm_k .gt. 0) then
      gbl_dihed(i)%atm_k = gbl_dihed(i)%atm_k / 3 + 1
    else
      gbl_dihed(i)%atm_k = gbl_dihed(i)%atm_k / 3 - 1
    end if
    if (gbl_dihed(i)%atm_l .gt. 0) then
      gbl_dihed(i)%atm_l = gbl_dihed(i)%atm_l / 3 + 1
    else
      gbl_dihed(i)%atm_l = gbl_dihed(i)%atm_l / 3 - 1
    end if
  end do

  fmtin = ifmt
  type = 'DIHEDRALS_WITHOUT_HYDROGEN'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_dihed(diheda_idx + i)%atm_i, &
                gbl_dihed(diheda_idx + i)%atm_j, &
                gbl_dihed(diheda_idx + i)%atm_k, &
                gbl_dihed(diheda_idx + i)%atm_l, &
                gbl_dihed(diheda_idx + i)%parm_idx, &
                i = 0, nphia - 1)

  do i = 0, nphia - 1
    gbl_dihed(diheda_idx + i)%atm_i = gbl_dihed(diheda_idx + i)%atm_i / 3 + 1
    gbl_dihed(diheda_idx + i)%atm_j = gbl_dihed(diheda_idx + i)%atm_j / 3 + 1
    if (gbl_dihed(diheda_idx + i)%atm_k .gt. 0) then
      gbl_dihed(diheda_idx + i)%atm_k = gbl_dihed(diheda_idx + i)%atm_k / 3 + 1
    else
      gbl_dihed(diheda_idx + i)%atm_k = gbl_dihed(diheda_idx + i)%atm_k / 3 - 1
    end if
    if (gbl_dihed(diheda_idx + i)%atm_l .gt. 0) then
      gbl_dihed(diheda_idx + i)%atm_l = gbl_dihed(diheda_idx + i)%atm_l / 3 + 1
    else
      gbl_dihed(diheda_idx + i)%atm_l = gbl_dihed(diheda_idx + i)%atm_l / 3 - 1
    end if
  end do

  fmtin = ifmt
  type = 'EXCLUDED_ATOMS_LIST'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_natex(i), i = 1, next)

! Read the h-bond parameters:

  fmtin = rfmt
  type = 'HBOND_ACOEF'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_asol(i), i = 1, nphb)

  fmtin = rfmt
  type = 'HBOND_BCOEF'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (gbl_bsol(i), i = 1, nphb)

  fmtin = rfmt
  type = 'HBCUT'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (dumd, i = 1, nphb)    ! Throw away hbcut. It is not used.

! Read isymbl, itree, join and irotat arrays.

  fmtin = afmt
  type = 'AMBER_ATOM_TYPE'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (atm_isymbl(i), i = 1, natom)

  fmtin = afmt
  type = 'TREE_CHAIN_CLASSIFICATION'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (atm_itree(i), i = 1, natom)

  fmtin = ifmt
  type = 'JOIN_ARRAY'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (idum, i = 1, natom) ! Throw away the join array. Not used.

  fmtin = ifmt
  type = 'IROTAT'
  call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

  read(prmtop, fmt) (idum, i = 1, natom)! Throw away the irotat array. Not used.
  
! Read the boundary condition stuff if needed:

  if (ifbox .gt. 0) then

    fmtin = ifmt
    type = 'SOLVENT_POINTERS'
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

    read(prmtop, fmt) idum, nspm, idum

  else

    nspm = 1

  end if

  allocate(atm_nsp(nspm), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_ints = num_ints + size(atm_nsp)

  ! We read the prmtop box and angle data here but ignore it, per standard
  ! amber behaviour.  The ifbox value also influences atm_nsp() setup though.

  if (ifbox .gt. 0) then

    fmtin = ifmt
    type = 'ATOMS_PER_MOLECULE'
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

    read(prmtop, fmt) (atm_nsp(i), i = 1, nspm)

    fmtin = rfmt
    type = 'BOX_DIMENSIONS'
    call nxtsec(prmtop, mdout, 1, fmtin, type, fmt, errcode)

    if (errcode == 0) then
      read(prmtop, fmt) dumd, dumd, dumd, dumd
    end if


  else

    atm_nsp(1) = natom

  end if

! Load the cap information if needed:

  if (ifcap .gt. 0) then

    fmtin = '(I6)'
    type = 'CAP_INFO'
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

    read(prmtop, fmt) idum   ! natcap

    fmtin = '(4E16.8)'
    type = 'CAP_INFO2'
    call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

    ! At least for now, cap data is thrown away.  As far as I know, all
    ! cap simulations do not have PBC.

    read(prmtop, fmt) dumd, dumd, dumd, dumd ! cutcap, xcap, ycap, zcap

! else

!   natcap = 0
!   cutcap = 0.0d0
!   xcap = 0.0d0
!   ycap = 0.0d0
!   zcap = 0.0d0

  end if

! Generalized Born, if appropriate:

  if (using_gb_potential .or. icnstph .eq. 2 .or. icnste .eq. 2 .or.  hybridgb .gt. 0 ) then

    if (using_gb_potential .or. icnstph .eq. 2 .or. icnste .eq. 2 .or. hybridgb .gt. 0) then
      if (errcode .eq. -1) then
        write(mdout, '(a,a)') error_hdr, &
          'GB calculations now require a new-style prmtop file'
        call mexit(6, 1)
      end if

      fmtin = rfmt
      type = 'RADII'
      call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

      read(prmtop, fmt) (atm_gb_radii(i), i = 1, natom)

      fmtin = rfmt
      type = 'SCREEN'
      call nxtsec(prmtop, mdout, 0, fmtin, type, fmt, errcode)

      read(prmtop, fmt) (atm_gb_fs(i), i = 1, natom)
    end if

    fmtin = ifmt
    type = 'ATOMIC_NUMBER'
    call nxtsec(prmtop, mdout, 1, fmtin, type, fmt, errcode)
      
    if (errcode .eq. 0) then
      write(mdout,'(2a)') info_hdr, 'Reading atomic numbers from topology file.'
      read(prmtop, fmt) (atm_atomicnumber(i), i = 1, natom)
      loaded_atm_atomicnumber = .true.
    else
      atm_atomicnumber(:) = 0
      if (igb .eq. 7 .or. igb .eq. 8 .or. gbsa .eq. 1 .or. &
          hybridgb .eq. 7 .or. hybridgb .eq. 8) then
        write(mdout,'(2a)') warn_hdr, &
          'ATOMIC_NUMBER section not found. Guessing atomic numbers from'
        write(mdout,'(2a)') extra_line_hdr, &
          'masses for GB parameters. Remake topology file with AmberTools 12+'
        write(mdout, '(2a)') extra_line_hdr, &
          'or add atomic numbers with ParmEd to remove this warning.'
      end if
    end if
   
  end if

  ! Extra point frames specified by the topology, not inferred
  fmtin = '(10I8)'
  type = 'VIRTUAL_SITE_FRAMES'
  call nxtsec(prmtop, mdout, 1, fmtin, type, fmt, errcode)
  if (errcode .eq. 0) then
    allocate(prm_ep_frames(numextra), stat=alloc_failed)
    if (alloc_failed .ne. 0) then
      write(mdout, '(2A)') error_hdr, 'Allocation for virtual site frames failed.'
      call mexit(mdout, 1)
    end if
    prm_ep_frames(:) = null_prm_ep_frame_rec
    do i = 1, numextra
      prm_ep_frames(i)%ep_cnt = 1
    end do
    write(mdout, '(2A)') info_hdr, 'Reading virtual site frames from topology file.'
    call count_custom_ep_frames(prmtop, mdout, numextra, num_custom_ep)
    call nxtsec(prmtop, mdout, 1, fmtin, type, fmt, errcode)
    read(prmtop, fmt) (prm_ep_frames(i)%extra_pnt, prm_ep_frames(i)%parent_atm, &
                       prm_ep_frames(i)%frame_atm1, prm_ep_frames(i)%frame_atm2, &
                       prm_ep_frames(i)%frame_atm3, prm_ep_frames(i)%type, i = 1, &
                       num_custom_ep)
    write(mdout, '(2A,I8,A)') info_hdr, 'Read ', num_custom_ep, ' virtual site frames.'
    if (num_custom_ep .lt. numextra) then
      write(mdout, '(2A,I8,A)') info_hdr, 'Other virtual sites will be inferred from &
                                          &bonding patterns.'
    end if
    type = 'VIRTUAL_SITE_FRAME_DETAILS'
    fmtin = rfmt
    call nxtsec(prmtop, mdout, 1, fmtin, type, fmt, errcode)
    read(prmtop, fmt) (prm_ep_frames(i)%d1, prm_ep_frames(i)%d2, prm_ep_frames(i)%d3, &
                       i = 1, num_custom_ep)
  else if (numextra .gt. 0) then
    write(mdout, '(2A)') info_hdr, &
                         'Virtual site frames will be inferred from bonding patterns.'
  end if
 
! Polarizabilities would need to be read here if they are supported.

! Fix up the typ-ico array for more efficient processing of vdw and 10_12
! interactions.  Basically, if the ico array entry is 0, then vdw processing
! is unnecessary for amber pme.
! Note this prevents division by zero in other parts of the code for things
! like the H in water which has a zero VDW radii. Note it only deals with both
! cn1 and cn2 being zero - there is nothing for dealing with the case where 
! just one of these is zero.
  do i = 1, ntypes * ntypes
    ic = typ_ico(i)
    if (ic .gt. 0) then
      if (gbl_cn1(ic) .eq. 0.d0 .and. gbl_cn2(ic) .eq. 0.d0) typ_ico(i) = 0
    else if (ic .lt. 0) then
      ic = -ic
      if (gbl_asol(ic) .eq. 0.d0 .and. gbl_bsol(ic) .eq. 0.d0) typ_ico(i) = 0
    end if
  end do

  ! Test for modification of off diagonal elements of the LJ terms.a NBFIX

  ! Currently the CPU code and the GPU code support NBFIX but we will leave
  ! the detection in place here and print an informational message to the mdout
  ! file. It may be necessary to deactivate support for this at some point on
  ! older GPUs if there is a large performance regression or if the atom type
  ! counts exceed a certain value so it makes sense to leave the test here 
  ! for now.

  ! Note while the frcmod and parmXX dat files do not support of diagonal
  ! entires in the LJ terms a lot of people have been hacking the prmtop to
  ! attempt to work around this - rather than fixing the actual underlying
  ! parmXX and frcmod file formats as they should do.
  ! Leap also does not support pairwise LJ terms but the
  ! prmtop file for some reason has the LJ terms expanded over all atom type
  ! pairs rather than just storing the individual atom type values. The GPU code,
  ! to save shared memory, used to recalculate the LJ parameters for the
  ! individual atom types from the diagonal and then use them. Thus if someone
  ! modified the off-diagonal elements in the prmtop then the GPU code would not
  ! give the answer expected. This code was designed to trap for this by computing the LJ params
  ! from the diagonal and checking all the off diagonal elements. The GPU code
  ! now supports these off diagonal modifications so for now this test is just
  ! to provide an informational message that Off-diagonal LJ terms are detected.
  !
  ! Note the GPU code only supports the default AMBER combining rules. That is
  ! the Lorentz-Bethelot rules where sigma(ij)=(sigma(i)+sigma(j))/2 and
  ! epsilon(ij) = sqrt(epsilon(i) x epsilon(j))
  ! If the prmtop entries for A and B have been built using other combining rules
  ! then we trap it here and quit. - Currently quit for CPU as well since I
  ! don't know of any test examples that use other combining rules.

  lj_offdiag = .false.

  allocate (sigma_tmp(ntypes), stat = errcode)
  if (errcode /= 0) then
    write(mdout,*) 'Allocation of sigma_tmp in init_amber_prmtop_dat failed.'
    call mexit(mdout, 1)
  end if
  allocate (epsilon_tmp(ntypes), stat = errcode)
  if (errcode /= 0) then
    write(mdout,*) 'Allocation of sigma_tmp in init_amber_prmtop_dat failed.'
    call mexit(mdout, 1)
  end if

  !1) Find Sigma[i] and Epsilon[i] from diagonal for 1 -> ntypes
  do i = 1, ntypes
    nbtype_tmp = typ_ico(ntypes*(i-1) + i)

#ifdef CUDA
    !Check for -ve nbtype_tmp which would indicate a 10-12 interaction which is
    !not supported.
    if (nbtype_tmp < 0) then
       write(mdout, '(a)') 'CUDA (GPU): Prmtop appears to have 10-12 terms present.'
       write(mdout, '(a)') '            The CUDA code does not currently support this.'
       write(mdout, '(a,i4,a,i4)') '            Atom Type i = ', i, 'Atom Type Index ico = ', nbtype_tmp
       call mexit(mdout,1)
    end if
#endif

    if (nbtype_tmp == 0) then
      !The A and B coefficients are zero. To avoid division by zero we set sigma
      !to 0.5d0.
      sigma_tmp(i) = 0.5d0
      epsilon_tmp(i) = 0.0d0
    else
      c1_tmp = gbl_cn1(nbtype_tmp)
      c2_tmp = gbl_cn2(nbtype_tmp)
      sig_tmp = (c1_tmp / c2_tmp)**(1.0d0/6.0d0)
      eps_tmp = (c2_tmp * c2_tmp)/(4.0d0 * c1_tmp)
      sigma_tmp(i) = 0.5d0 * sig_tmp
      epsilon_tmp(i) = sqrt(eps_tmp)
    endif
  end do

  !2) Check that using these Sigma and Epsilon values that we get the A and B
  !   (cn1 and cn2) values in the prmtop.

  loop_count = 0
  do i = 1, ntypes
    do j = 1, i
      loop_count = loop_count + 1
   
      aij_tmp = 4.0d0 * (sigma_tmp(i) + sigma_tmp(j))**12 * &
                epsilon_tmp(i) * epsilon_tmp(j)
      bij_tmp = 4.0d0 * (sigma_tmp(i) + sigma_tmp(j))**6 * &
                epsilon_tmp(i) * epsilon_tmp(j)

      ! Compare to a small difference here as absolute comparisons are
      ! unreliable for floating point.

      if (abs(aij_tmp - gbl_cn1(loop_count)) > (1.0d-7*abs(aij_tmp)) ) then
        lj_offdiag=.true.
      end if

! OLD GPU FATAL CODE
!      if (abs(aij_tmp - gbl_cn1(loop_count)) > (1.0d-7*abs(aij_tmp)) ) then
!        write(mdout, '(a)') &
!         'CUDA (GPU): Prmtop appears to contain modified off-diagonal elements &
!         &for VDW A coef.'
!        write(mdout, '(a)') &
!         '            GPU code does NOT support such force field modifications.'
!        write(mdout, '(a,f22.10,a,f22.10)') &
!         '            Found Aij = ', aij_tmp, ' expected ', gbl_cn1(loop_count)
!        write(mdout, '(a,i4,a,i4,a,i6)') &
!         '            i = ', i, ' j = ', j, ' Matrix element = ', loop_count
!        call mexit(mdout,1)
!      end if
!      if (abs(bij_tmp - gbl_cn2(loop_count)) > (1.0d-7*abs(bij_tmp)) ) then
!        write(mdout, '(a)') &
!         'CUDA (GPU): Prmtop appears to contain modified off-diagonal elements &
!         &for VDW B coef.'
!        write(mdout, '(a)') &
!         '            GPU code does NOT support such force field modifications.'
!        write(mdout, '(a,f22.10,a,f18.10)') &
!         '            Found Bij = ', bij_tmp, ' expected ', gbl_cn2(loop_count)
!        write(mdout, '(a,i4,a,i4,a,i6)') &
!         '            i = ', i, ' j = ', j, ' Matrix element = ', loop_count
!        call mexit(mdout,1)
!      end if
    end do
  end do

  deallocate (sigma_tmp, stat = errcode)
  if (errcode /= 0) then
    write(mdout,*) 'Deallocation of sigma_tmp in init_amber_prmtop_dat failed.'
    call mexit(6, 1)
  end if
  deallocate (epsilon_tmp, stat = errcode)
  if (errcode /= 0) then
    write(mdout,*) 'Deallocation of sigma_tmp in init_amber_prmtop_dat failed.'
    call mexit(6, 1)
  end if

! Print infomational message if off diagonal LJ found
  if ( lj_offdiag ) then
    write(mdout, '(a)') '| INFO: Off Diagonal (NBFIX) LJ terms found in prmtop file.'
    write(mdout, '(a)') '|       The prmtop file has been modified to support atom'
    write(mdout, '(a)') '|       type based pairwise Lennard-Jones terms.'
  end if

! --- End testing for off diagonal VDW terms ---

! Duplicate dihedral pointers for vector ephi.  nphih and nphia may be increased
! in these calls.

  orig_nphih = nphih        ! Saved for reporting later
  orig_nphia = nphia        ! Saved for reporting later

  call duplicate_dihedrals(nphih, gbl_dihed, gbl_pn)
  call duplicate_dihedrals(nphia, gbl_dihed(diheda_idx), gbl_pn)

! Pre-calculate and save some parameters for vector ephi:

  call calc_dihedral_parms(nptra, gbl_pk, gbl_pn, gbl_phase, gbl_gamc, &
                           gbl_gams, gbl_ipn, gbl_fmn)

! Print prmtop parameters header ('resource use'):

  write(mdout, 1020)

  ! This is after-the-fact, but done here for output consistency with sander 8:
  if (ntb .gt. 0) then
    if( NC_checkRestart(inpcrd_name) ) then
      write(mdout, '(a)')' getting box info from netcdf restart file'
    else
      write(mdout, '(a)')' getting new box info from bottom of inpcrd'
    endif
  end if

! Print prmtop parameters:

!  print*,natom, ntypes, nbonh, mbona, ntheth, mtheta,&
!          orig_nphih, mphia, nhparm, nparm, next, &
!          nres, nbona, ntheta, orig_nphia, numbnd, numang, &
!          nptra, natyp, nphb, ifbox, nmxrs, ifcap,numextra, ncopy
  if (C4Pairwise .gt. 0) then
  write(mdout, 1031) natom, ntypes, nbonh, mbona, ntheth, mtheta, &
                     orig_nphih, mphia, nhparm, nparm, next, &
                     nres, nbona, ntheta, orig_nphia, numbnd, numang, &
                     nptra, natyp, nphb, ifbox, nmxrs, ifcap, numextra, C4Pairwise, ncopy ! New2022
  else
  write(mdout, 1030) natom, ntypes, nbonh, mbona, ntheth, mtheta, &
                     orig_nphih, mphia, nhparm, nparm, next, &
                     nres, nbona, ntheta, orig_nphia, numbnd, numang, &
                     nptra, natyp, nphb, ifbox, nmxrs, ifcap, numextra, ncopy ! New2022
  end if

  if (using_gb_potential .or. icnstph .eq. 2 .or. icnste .eq. 2 .or. hybridgb .gt. 0) then 

    fmtin = afmt
    type = 'RADIUS_SET'
    call nxtsec(prmtop, mdout, 1, fmtin, type, fmt, errcode)

    ! Write implicit solvent radius and screening info to mdout:

    if (errcode .eq. 0) then
      read(prmtop, fmt) type    ! note re-use of type to hold result, a hack...
      write(mdout, '(a,a)') &
        ' Implicit solvent radii are ', type
    end if

    ! If igb .eq. 7 use special S_x screening params; here we overwrite
    ! the tinker values read from the prmtop.

    if (igb .eq. 7) then

      write(mdout,'(a)') &
        ' Replacing prmtop screening parameters with GBn (igb=7) values'

      do i = 1, natom

        if (loaded_atm_atomicnumber) then
          atomicnumber = atm_atomicnumber(i)
        else
          call get_atomic_number(atm_igraph(i), atm_mass(i), atomicnumber)
        end if

        if (atomicnumber .eq. 6) then
          atm_gb_fs(i) = 4.84353823306d-1
        else if (atomicnumber .eq. 1) then
          atm_gb_fs(i) = 1.09085413633d0
        else if (atomicnumber .eq. 7) then
          atm_gb_fs(i) = 7.00147318409d-1
        else if (atomicnumber .eq. 8) then
          atm_gb_fs(i) = 1.06557401132d0
        else if (atomicnumber .eq. 16) then
          atm_gb_fs(i) = 6.02256336067d-1
        else
          atm_gb_fs(i) = 5.d-1 ! not optimized
        end if

      end do
    
    else if (igb .eq. 8) then
      
      write(mdout, '(a)') &
        ' Replacing prmtop screening parameters with GBn2 (igb=8) values'

      do i = 1, natom

        if(loaded_atm_atomicnumber) then
          atomicnumber = atm_atomicnumber(i)
        else
          call get_atomic_number(atm_igraph(i), atm_mass(i), atomicnumber)
        end if

        ! The screen_ variables are found in mdin_ctrl_dat_mod

        call isnucat(nucat,i,nres,60,gbl_res_atms(1:nres),gbl_labres(1:nres))
        ! check if atom belong to nucleic or protein residue
        ! 37 = the total number of different nucleic acid residue type in AMBER
        ! Need to update this number if adding new residue type
        ! update isnucat function in gb_ene.F90 too

        ! RCW - HACKY CODE ALERT!
        !       This (and the isnucat subroutine) is incredibly fragile and people will forget
        !       to update this. Plus one can in principal use any residue name they want. This
        !       is a total HACK and should NOT be here. Make it robust and move it to leap where
        !       it belogs please.

        if (nucat == 1) then
            ! If atom belongs to nuc, use GBNeck2nu pars
            if (atomicnumber .eq. 6) then
              atm_gb_fs(i) = screen_cnu
            else if (atomicnumber .eq. 1) then
              atm_gb_fs(i) = screen_hnu
            else if (atomicnumber .eq. 7) then
              atm_gb_fs(i) = screen_nnu
            else if (atomicnumber .eq. 8) then
              atm_gb_fs(i) = screen_onu
            else if (atomicnumber .eq. 15) then
              atm_gb_fs(i) = screen_pnu
            else
              atm_gb_fs(i) = 5.d-1 ! not optimized
            end if
        else
            !if atom does not belong to nuc, use protein pars
            if (atomicnumber .eq. 6) then
              atm_gb_fs(i) = screen_c
            else if (atomicnumber .eq. 1) then
              atm_gb_fs(i) = screen_h
            else if (atomicnumber .eq. 7) then
              atm_gb_fs(i) = screen_n
            else if (atomicnumber .eq. 8) then
              atm_gb_fs(i) = screen_o
            else if (atomicnumber .eq. 16) then
              atm_gb_fs(i) = screen_s
            else if (atomicnumber .eq. 15) then
              atm_gb_fs(i) = screen_p
            else
              atm_gb_fs(i) = 5.d-1 ! not optimized
            end if
         end if

      end do

    end if

    ! Put fs(i) * (rborn(i) - offset) into the "fs" array:

    gb_fs_max = 0.d0

    do i = 1, natom
      atm_gb_fs(i) = atm_gb_fs(i) * (atm_gb_radii(i) - offset)
      gb_fs_max = max(gb_fs_max, atm_gb_fs(i))
    end do

  end if

  return

 1020 format(80('-')/'   1.  RESOURCE   USE: ',/80('-')/)

 1030 format(t2, &
        'NATOM  = ', i7, ' NTYPES = ', i7, ' NBONH = ', i7, ' MBONA  = ', i7, &
      /' NTHETH = ', i7, ' MTHETA = ', i7, ' NPHIH = ', i7, ' MPHIA  = ', i7, &
      /' NHPARM = ', i7, ' NPARM  = ', i7, ' NNB   = ', i7, ' NRES   = ', i7, &
      /' NBONA  = ', i7, ' NTHETA = ', i7, ' NPHIA = ', i7, ' NUMBND = ', i7, &
      /' NUMANG = ', i7, ' NPTRA  = ', i7, ' NATYP = ', i7, ' NPHB   = ', i7, &
      /' IFBOX  = ', i7, ' NMXRS  = ', i7, ' IFCAP = ', i7, ' NEXTRA = ', i7, &
      /' NCOPY  = ', i7/) ! New2022
 
 1031 format(t2, &
        'NATOM  = ', i7, ' NTYPES = ', i7, ' NBONH = ', i7, ' MBONA  = ', i7, &
      /' NTHETH = ', i7, ' MTHETA = ', i7, ' NPHIH = ', i7, ' MPHIA  = ', i7, &
      /' NHPARM = ', i7, ' NPARM  = ', i7, ' NNB   = ', i7, ' NRES   = ', i7, &
      /' NBONA  = ', i7, ' NTHETA = ', i7, ' NPHIA = ', i7, ' NUMBND = ', i7, &
      /' NUMANG = ', i7, ' NPTRA  = ', i7, ' NATYP = ', i7, ' NPHB   = ', i7, &
      /' IFBOX  = ', i7, ' NMXRS  = ', i7, ' IFCAP = ', i7, ' NEXTRA = ', i7, &
      /' C4Pairwise  = ', i7, ' NCOPY  = ', i7/) ! New2022
end subroutine init_amber_prmtop_dat

end subroutine init_prmtop_dat

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  bcast_amber_prmtop_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine bcast_amber_prmtop_dat

  use pmemd_lib_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use charmm_mod, only : charmm_active

  implicit none

  integer               :: alloc_failed
  integer               :: bytes_per_unit
  integer               :: num_ints, num_reals  ! returned values discarded

!PHMD
  call mpi_bcast(iphmd, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)

  call mpi_bcast(natom, prmtop_int_cnt, mpi_integer, 0, pmemd_comm, &
                 err_code_mpi)

  call mpi_bcast(charmm_active, 1, mpi_logical, 0, pmemd_comm, &
                 err_code_mpi)

! There is only cap data, which we currently don't use, so we comment this
! out:
!
! call mpi_bcast(cutcap, prmtop_dbl_cnt, mpi_double_precision, 0, &
!                pmemd_comm, err_code_mpi)

  if (.not. master) then
    num_ints = 0
    num_reals = 0
    call alloc_amber_prmtop_mem(num_ints, num_reals)
  end if

  call mpi_bcast(atm_qterm, natom, mpi_double_precision, 0, &
                 pmemd_comm, err_code_mpi)

  if (using_gb_potential .or. icnstph .eq. 2 .or. icnste .eq. 2 .or. hybridgb .gt. 0) then
     call mpi_bcast(atm_atomicnumber, natom, mpi_integer, 0, &
                    pmemd_comm, err_code_mpi)
     call mpi_bcast(loaded_atm_atomicnumber, 1, mpi_logical, 0, &
                    pmemd_comm, err_code_mpi)
  end if

  call mpi_bcast(atm_mass, natom, mpi_double_precision, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(atm_iac, natom, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(typ_ico, ntypes * ntypes, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(gbl_res_atms, nres + 1, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(atm_igraph, natom, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(gbl_labres, nres, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  if (nphb > 0) then
    call mpi_bcast(gbl_asol, nphb, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    call mpi_bcast(gbl_bsol, nphb, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
  end if

  if (nttyp > 0) then
    call mpi_bcast(gbl_cn1, nttyp, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    call mpi_bcast(gbl_cn2, nttyp, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    call mpi_bcast(gbl_cn6, nttyp, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    if (C4Pairwise > 0) then
      call mpi_bcast(gbl_cn7, 3*C4Pairwise, mpi_integer, 0, & ! New2022
                     pmemd_comm, err_code_mpi)
      call mpi_bcast(gbl_cn8, C4Pairwise, mpi_double_precision, 0, &
                     pmemd_comm, err_code_mpi)
    end if
    if (charmm_active) then
      call mpi_bcast(gbl_cn114, nttyp, mpi_double_precision, 0, &
                     pmemd_comm, err_code_mpi)
      call mpi_bcast(gbl_cn214, nttyp, mpi_double_precision, 0, &
                     pmemd_comm, err_code_mpi)
    end if
  end if

  if (numbnd > 0) then
    call mpi_bcast(gbl_rk, numbnd, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    call mpi_bcast(gbl_req, numbnd, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    call mpi_bcast(gbl_rpullk, numbnd, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    call mpi_bcast(gbl_rpulleq, numbnd, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    call mpi_bcast(gbl_rpressk, numbnd, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    call mpi_bcast(gbl_rpresseq, numbnd, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
  end if

  if (numang > 0) then
    call mpi_bcast(gbl_tk, numang, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    call mpi_bcast(gbl_teq, numang, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
  end if

  if (nptra > 0) then
    call mpi_bcast(gbl_pk, nptra, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    call mpi_bcast(gbl_pn, nptra, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    call mpi_bcast(gbl_phase, nptra, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    call mpi_bcast(gbl_gamc, nptra, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    call mpi_bcast(gbl_gams, nptra, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    call mpi_bcast(gbl_fmn, nptra, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    call mpi_bcast(gbl_ipn, nptra, mpi_integer, 0, &
                   pmemd_comm, err_code_mpi)
    call mpi_bcast(gbl_one_scee, nptra, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    call mpi_bcast(gbl_one_scnb, nptra, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
  end if

  if (charmm_active) then
    !Urey Bradley
    if (nubtypes > 0) then
      call mpi_bcast(gbl_ub_rk, nubtypes, mpi_double_precision, 0, &
                     pmemd_comm, err_code_mpi)
      call mpi_bcast(gbl_ub_r0, nubtypes, mpi_double_precision, 0, &
                     pmemd_comm, err_code_mpi)
    end if

    !Charmm Impropers
    if (nimprtyp > 0) then
      call mpi_bcast(gbl_imp_pk, nimprtyp, mpi_double_precision, 0, &
                     pmemd_comm, err_code_mpi)
      call mpi_bcast(gbl_imp_phase, nimprtyp, mpi_double_precision, 0, &
                     pmemd_comm, err_code_mpi)
    end if
  end if
  
  ! if CMAP is active 
  ! CMAP TODO: 24 is a hard coded hack
    if ( cmap_type_count > 0 ) then
      call mpi_bcast(gbl_cmap_res, cmap_type_count, mpi_integer, 0, &
                     pmemd_comm, err_code_mpi)

      call mpi_bcast(gbl_cmap_grid_step_size, cmap_type_count, mpi_integer, 0, &
                     pmemd_comm, err_code_mpi)

      call mpi_bcast(gbl_cmap_grid, cmap_type_count*24*24, mpi_double_precision, 0, &
                     pmemd_comm, err_code_mpi)

      call mpi_bcast(gbl_cmap_dPhi, cmap_type_count*24*24, mpi_double_precision, 0, &
                     pmemd_comm, err_code_mpi)
      call mpi_bcast(gbl_cmap_dPsi, cmap_type_count*24*24, mpi_double_precision, 0, &
                     pmemd_comm, err_code_mpi)
      call mpi_bcast(gbl_cmap_dPhi_dPsi, cmap_type_count*24*24, mpi_double_precision, 0, &
                     pmemd_comm, err_code_mpi)

   end if 





  ! Generalized Born:

  if (using_gb_potential .or. icnstph .eq. 2 .or. icnste .eq. 2 .or. hybridgb .gt. 0) then
    ! For IGB = 8, we need to know about atm_isymbl
    call mpi_bcast(atm_gb_fs, natom, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    call mpi_bcast(atm_gb_radii, natom, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    call mpi_bcast(atm_isymbl, natom * 4, mpi_character, 0, &
                   pmemd_comm, err_code_mpi)
  end if

  ! Polarizabilities would need to be broadcast if they are supported.

  call get_bytesize(gbl_bond(1), gbl_bond(2), bytes_per_unit)

  call mpi_bcast(gbl_bond, size(gbl_bond) * bytes_per_unit, mpi_byte, 0, &
                 pmemd_comm, err_code_mpi)

  call get_bytesize(gbl_angle(1), gbl_angle(2), bytes_per_unit)

  call mpi_bcast(gbl_angle, size(gbl_angle) * bytes_per_unit, mpi_byte, 0, &
                 pmemd_comm, err_code_mpi)

  call get_bytesize(gbl_dihed(1), gbl_dihed(2), bytes_per_unit)

  call mpi_bcast(gbl_dihed, size(gbl_dihed) * bytes_per_unit, mpi_byte, 0, &
                 pmemd_comm, err_code_mpi)

  if (charmm_active) then
    call get_bytesize(gbl_angle_ub(1), gbl_angle_ub(2), bytes_per_unit)

    call mpi_bcast(gbl_angle_ub, size(gbl_angle_ub) * bytes_per_unit, mpi_byte, 0, &
                   pmemd_comm, err_code_mpi)

    call get_bytesize(gbl_dihed_imp(1), gbl_dihed_imp(2), bytes_per_unit)

    call mpi_bcast(gbl_dihed_imp, size(gbl_dihed_imp) * bytes_per_unit, mpi_byte, 0, &
                   pmemd_comm, err_code_mpi)
  end if
    ! if CMAP is active
  if (cmap_term_count > 0) then
    call get_bytesize(gbl_cmap(1), gbl_cmap(2), bytes_per_unit)

    call mpi_bcast(gbl_cmap, size(gbl_cmap) * bytes_per_unit, mpi_byte, 0, &
                   pmemd_comm, err_code_mpi)
  end if



  call mpi_bcast(atm_numex, natom, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(gbl_natex, next * next_mult_fac, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(atm_nsp, nspm, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  ! no need to bcast atm_tree. We need atm_isymbl only for some GB models

  return

end subroutine bcast_amber_prmtop_dat
#endif

!*******************************************************************************
!
! Subroutine:  alloc_amber_prmtop_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_amber_prmtop_mem(num_ints, num_reals)

  use pmemd_lib_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use charmm_mod, only : charmm_active
#ifndef NOXRAY
  use xray_globals_module, only : xray_active, bulk_solvent_model
#endif

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals     

! Local variables:

  integer                       :: alloc_failed

  ! Generic allocations required for just about everything:

  allocate(atm_qterm(natom), &
           atm_mass(natom), &
           atm_iac(natom), &
           typ_ico(ntypes * ntypes), &
           gbl_cn1(nttyp), &
           gbl_cn2(nttyp), &
           gbl_cn6(nttyp), &
           gbl_cn7(3*C4Pairwise), & ! New2022
           gbl_cn8(C4Pairwise), &   ! New2022
           gbl_rk(numbnd), &
           gbl_req(numbnd), &
           gbl_rpressk(numbnd), &
           gbl_rpresseq(numbnd), &
           gbl_rpullk(numbnd), &
           gbl_rpulleq(numbnd), &
           gbl_tk(numang), &
           gbl_teq(numang), &
           atm_igraph(natom), &
           gbl_res_atms(nres + 1), &
           gbl_labres(nres), &
           stat = alloc_failed)

  if (charmm_active) then
    allocate( &
             gbl_cn114(nttyp), &
             gbl_cn214(nttyp), &
             gbl_ub_rk(nubtypes),&
             gbl_ub_r0(nubtypes),&
             gbl_imp_pk(nimprtyp),& 
             gbl_imp_phase(nimprtyp),&
             stat = alloc_failed)
  end if 
  ! if CMAP is active
  if (cmap_term_count > 0) then
    allocate( &
             gbl_cmap_res(cmap_type_count), &
             gbl_cmap_grid_step_size(cmap_type_count), &
             gbl_cmap_grid(cmap_type_count,24,24),& !24 is a hardcoded hack
             gbl_cmap_dPhi(cmap_type_count,24,24),&
             gbl_cmap_dPsi(cmap_type_count,24,24),&
             gbl_cmap_dPhi_dPsi(cmap_type_count,24,24),&
             stat = alloc_failed)
  end if

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_reals = num_reals + size(atm_qterm) + size(atm_mass) + &
                          size(gbl_cn1) + size(gbl_cn2) + size(gbl_rk) + &
                          size(gbl_req) + size(gbl_rpullk) + size(gbl_rpulleq) + &
                          size(gbl_rpressk) + size(gbl_rpresseq) + size(gbl_tk) + &
                          size(gbl_teq) + size(gbl_cn6) + size(gbl_cn8) ! New2022

  num_ints = num_ints + size(atm_iac) + size(typ_ico) + &
                        size(atm_igraph) + size(gbl_res_atms) + &
                        size(gbl_labres) + size(gbl_cn7) ! New2022

  if (charmm_active) then
    num_reals = num_reals + size(gbl_cn114) + size(gbl_cn214) + &
                            size(gbl_ub_rk) + size(gbl_ub_r0) + size(gbl_imp_pk) + &
                            size(gbl_imp_phase)
  end if
  ! if CMAP is active
  if (cmap_term_count > 0) then
    num_reals = num_reals + &
                            size(gbl_cmap_res) + size(gbl_cmap_grid_step_size) + &
                            size(gbl_cmap_grid) + size(gbl_cmap_dPhi) + &
                            size(gbl_cmap_dPsi) + size(gbl_cmap_dPhi_dPsi)
  end if


  if (nphb .ne. 0) then
    allocate(gbl_asol(nphb), &
             gbl_bsol(nphb), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    num_reals = num_reals + size(gbl_asol) + size(gbl_bsol)

  end if

  ! Allocate space for atom bond, angle and dihedral information:

  allocate(gbl_bond(gbl_bond_allocsize), &
           gbl_angle(gbl_angle_allocsize), &
           gbl_dihed(gbl_dihed_allocsize), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  ! Zero the entire dihedral array since it is not entirly used as a result
  ! the uninitialized padding bytes cause false positivies with Valgrind
  gbl_dihed = null_dihed_rec
  
  if (charmm_active) then
    allocate( &
             gbl_angle_ub(gbl_angle_ub_allocsize), &
             gbl_dihed_imp(gbl_dihed_imp_allocsize), &
             stat = alloc_failed)
  
    if (alloc_failed .ne. 0) call setup_alloc_error
  end if
  ! if CHARMM or CMAP is active
  if (charmm_active .or. cmap_term_count > 0) then
    allocate( &
             gbl_cmap(gbl_cmap_allocsize), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error
  end if


  num_ints = num_ints + size(gbl_bond) * bond_rec_ints + &
             size(gbl_angle) * angle_rec_ints + &
             size(gbl_dihed) * dihed_rec_ints 

  if (charmm_active) then
    num_ints = num_ints + size(gbl_angle_ub) * angle_ub_rec_ints + &
               size(gbl_dihed_imp) * dihed_imp_rec_ints
  end if
  ! if CMAP is active
  if (cmap_term_count > 0) then
    num_ints = num_ints + &
               size(gbl_cmap) * cmap_rec_ints
  end if


  if (nptra > 0) then
    allocate(gbl_pk(nptra), &
             gbl_pn(nptra), &
             gbl_phase(nptra), &
             gbl_gamc(nptra), &
             gbl_gams(nptra), &
             gbl_fmn(nptra), &
             gbl_one_scee(nptra), &
             gbl_one_scnb(nptra), &
             gbl_ipn(nptra), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    num_reals = num_reals + size(gbl_pk) + &
                            size(gbl_pn) + &
                            size(gbl_phase) + &
                            size(gbl_gamc) + &
                            size(gbl_gams) + &
                            size(gbl_fmn) + &
                            size(gbl_one_scee) + &
                            size(gbl_one_scnb)

    num_ints = num_ints + size(gbl_ipn)

    ! We know these get initialized properly, so no need to zero.

  end if

  ! Allocations for Generalized Born: (isymbl, atomicnumber
  ! needed for igb = 8)

!RCW - why are these being done when igb /= 8? - The comment above
!      implies just for IGB=8.

  if (using_gb_potential .or. icnstph .eq. 2 .or. icnste .eq. 2 .or. hybridgb .gt. 0) then 

    allocate(atm_gb_fs(natom), atm_gb_radii(natom), atm_atomicnumber(natom), &
             atm_isymbl(natom), stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    num_reals = num_reals + size(atm_gb_fs) + size(atm_gb_radii)

    num_ints = num_ints + size(atm_isymbl) + size(atm_atomicnumber)

  end if
#ifndef NOXRAY
  if (xray_active .and. .not. allocated(atm_atomicnumber)) then
    allocate(atm_atomicnumber(natom), stat = alloc_failed)
    
    if (alloc_failed .ne. 0) call setup_alloc_error

    num_ints = num_ints + size(atm_atomicnumber)
 
  end if
#endif

  ! Polarizabilities would need to have space alloc'd if they are supported.

  ! Allocation of stuff that will later be deallocated:

  ! gbl_natex dimension is allocated as 2*next here since it can be expanded inside
  ! the extra point code to be more than next. From the prmtop though only next
  ! elements are read and subsequently broadcast.
  ! Bugfix note: gbl_natex dimension increased to provide workspace for
  ! extra points code...  This is a real hack that should be fixed...

  allocate(atm_numex(natom), &
           gbl_natex(next * next_mult_fac), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_ints = num_ints + size(atm_numex) + size(gbl_natex)

  gbl_natex(:) = 0

  ! Allocation of stuff that will later be deallocated that is only used in
  ! the master:

  if (master) then

    allocate(atm_itree(natom), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    num_ints = num_ints + size(atm_itree)

    ! In this case, we haven't allocated atm_isymbl yet, so do it here!

    if (.not. using_gb_potential .and. icnstph .ne. 2 .and. icnste .ne. 2 .and. hybridgb .le. 0) then 

      allocate(atm_isymbl(natom), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error

      num_ints = num_ints + size(atm_isymbl)

    end if

  else
  
    allocate(atm_nsp(nspm), &
           stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    num_ints = num_ints + size(atm_nsp)

  end if

  return

end subroutine alloc_amber_prmtop_mem

!*******************************************************************************
!
! Subroutine:  duplicate_dihedrals
!
! Description:  Duplicates pointers to multi-term dihedrals for vector ephi.
!               H-atom diheds are duplicated at the end of the h-atom lists,
!               but heavy atom diheds must be inserted between the end of the
!               heavy-atom diheds.  In order to use this code, extra space MUST
!               be allocated (2 * max_dihed_dups).
!
! Author:  George Seibel
!              
!*******************************************************************************

subroutine duplicate_dihedrals(nphi, dihed, pn)

  use gbl_constants_mod
  use parallel_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer               :: nphi         ! num real dihedrals
  type(dihed_rec)       :: dihed(*)     ! the dihedral record
  double precision      :: pn(*)        ! periodicity; negative if multi-term,
                                        ! read until + encountered:
! Local variables:

  integer               :: i
  integer               :: ic           ! working parameter pointer:
  integer               :: ndup         ! num of diheds duplicated:

  ndup = 0

  do i = 1, nphi

    ic = dihed(i)%parm_idx

    do while (pn(ic) .lt. 0)

      ndup = ndup + 1

      if (ndup .gt. max_dihed_dups) then

        write(mdout,'(a,a,i5,a)') error_hdr, &
          'max_dihed_dups = ',max_dihed_dups,' exceeded.'

        write(mdout,'(a,a,i5)') extra_line_hdr, &
          'Set max_dihed_dups =', ndup
        call mexit(6, 1)

      else 

     ! Duplicate pointers of multi-term dihedrals:

        dihed(nphi + ndup)%atm_i  = dihed(i)%atm_i
        dihed(nphi + ndup)%atm_j  = dihed(i)%atm_j
        dihed(nphi + ndup)%atm_k  = dihed(i)%atm_k
        dihed(nphi + ndup)%atm_l  = dihed(i)%atm_l

     ! But use the NEXT parameter pointer:

        dihed(nphi + ndup)%parm_idx = ic + 1

      end if

      ! Check for a third or higher term:

      ic = ic + 1

    end do

  end do

  nphi = nphi + ndup

  write(mdout,'(a,i5,a,/)') '| Duplicated',ndup,' dihedrals'

  return

end subroutine duplicate_dihedrals

!*******************************************************************************
!
! Subroutine:  calc_dihedral_parms
!
! Description: Routine to calc additional parameters for the vectorised ephi.
!
!*******************************************************************************

subroutine calc_dihedral_parms(numphi, pk, pn, phase, gamc, gams, ipn, fmn)

  use gbl_constants_mod

  implicit none

! Formal arguments:

  integer           numphi
  double precision  pk(*)
  double precision  pn(*)
  double precision  phase(*)
  double precision  gamc(*)
  double precision  gams(*)
  integer           ipn(*)
  double precision  fmn(*)

! Local variables:

  double precision  dum
  double precision  dumc
  double precision  dums
  double precision  four
  integer           i
  double precision  one
  double precision  pim
  double precision  tenm3
  double precision  tenm10
  double precision  zero

  data zero, one, tenm3, tenm10/0.0d+00, 1.0d+00, 1.0d-03, 1.0d-06/
  data four /4.0d+00/

  pim = four*atan(one)

  do i = 1, numphi
    dum = phase(i)
    if (abs(dum - PI) .le. tenm3) dum = sign(pim, dum)
    dumc = cos(dum)
    dums = sin(dum)
    if (abs(dumc) .le. tenm10) dumc = zero
    if (abs(dums) .le. tenm10) dums = zero
    gamc(i) = dumc*pk(i)
    gams(i) = dums*pk(i)
    fmn(i) = one
    if (pn(i) .le. zero) fmn(i) = zero
    pn(i) = abs(pn(i))
    ipn(i) = int(pn(i) + tenm3)
  end do

  return

end subroutine calc_dihedral_parms

! checking if an atom belongs to nucleic or not
! make sure to add new DNA/RNA residue name to this list
! anyway not to hardcoded nucnamenum? (=60)
! tag: gbneck2nu
subroutine isnucat(nucat,atomindex,nres,nucnamenum,ipres,lbres)
    implicit none
    integer, intent(in) :: atomindex,nres,nucnamenum,ipres(nres)
    integer :: j,k,nucat
    character(4) :: nucname(nucnamenum),lbres(nres)
    !nucnamenum = 60
    nucat = 0
    nucname = (/ &
     "A   ", & !1
     "A3  ", &
     "A5  ", &
     "AN  ", &
     "C   ", & !5
     "C3  ", &
     "C5  ", &
     "CN  ", &
     "DA  ", &
     "DA3 ", & !10
     "DA5 ", &
     "DAN ", &
     "DC  ", &
     "DC3 ", &
     "DC5 ", & !15
     "DCN ", &
     "DG  ", &
     "DG3 ", &
     "DG5 ", &
     "DGN ", & !20
     "DT  ", &
     "DT3 ", &
     "DT5 ", &
     "DTN ", &
     "G   ", & !25
     "G3  ", &
     "G5  ", &
     "8OG ", &
     "GN  ", &
     "U   ", & !30
     "U3  ", &
     "U5  ", &
     "UN  ", &
     "AP  ", &
     "DAP ", & !35
     "CP  ", &
     "AF2 ", &
     "AF5 ", &
     "AF3 ", &
     "GF2 ", & !40
     "GF3 ", &
     "GF5 ", &
     "CF2 ", &
     "CFZ ", &
     "CF3 ", & !45
     "CF5 ", &
     "UF2 ", &
     "UF5 ", &
     "UF3 ", &
     "UFT ", & !50
     "OMC ", &
     "OMG ", &
     "MC3 ", &
     "MC5 ", &
     "MG3 ", & !55
     "MG5 ", &
     "SIC ", &
     "NAM ", &
     "NM5 ", &
     "CAP "  /)!60

    do j=1,nres-1
        !find residue to which atom i belongs
        if (ipres(j)  <= atomindex .and. atomindex < ipres(j+1)) then
            do k=1,nucnamenum
                !if atom i belongs to nucleic, return 1
                if (lbres(j)(1:4) == nucname(k)(1:4)) then
                    nucat = 1
                end if
            end do
        end if
    end do
    if (atomindex >= ipres(nres)) then
       do k=1,nucnamenum
           !if atom i belongs to nucleic, return 1
           if (lbres(j)(1:4) == nucname(k)(1:4)) then
               nucat = 1
           end if
       end do
    end if
end subroutine isnucat

end module prmtop_dat_mod
