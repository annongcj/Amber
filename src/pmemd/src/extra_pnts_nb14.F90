#include "copyright.i"
#include "hybrid_datatypes.i"

!*******************************************************************************!
! Module:  extra_pnts_nb14_mod
!
! Description: <TBS>
!
!*******************************************************************************

module extra_pnts_nb14_mod

  use gbl_datatypes_mod
  use file_io_dat_mod, only : mdout
#ifdef _OPENMP_
  use omp_lib_kinds
#endif

  implicit none

! Frames are centered on an atom parent_atm that "owns" extra points
! ep_cnt are number of extra points attached to parent_atm
! extra_pnt are pointers to (at most 2) extra points attached to parent_atm
! type is 1 if atom has at least 2 other bonds, pref to heavy atoms
! type is 2 if atom has only one other bond e.g. O6 in guanine
! frame_atm1 frame_atm2 and frame_atm3 are atom nums defining frame
! loc_frame are local coords

! Original code by Tom Darden based on Jim Caldwell's ideas.
!
! All atom names beginning with EP are considered extra points.
!
! If you want bond angle and dihedral forces as usual for lone pairs
! (i.e. as if they are amber atoms) then set frameon = 0.
!
! If frameon .eq. 1, (DEFAULT) the bonds, angles and dihedral interactions
! involving the lone pairs / extra points are removed except for constraints
! added during parm. The lone pairs are kept in ideal geometry relative to
! local atoms, and resulting torques are transferred to these atoms.
!
! If chngmask .eq. 1 (DEFAULT), new 1-1, 1-2, 1-3 and 1-4 interactions are
! calculated. An extra point belonging to an atom has a 1-1 interaction with it,
! and participates in any 1-2, 1-3 or 1-4 interaction that atom has.
! For example, suppose (excusing the geometry)
! C1,C2,C3,C4 form a dihedral and each has 1 extra point attached as below
!
!           C1------C2------C3---------C4
!           |        |       |         |
!           |        |       |         |
!          Ep1      Ep2     Ep3       Ep4
!
! The 1-4 interactions include  C1&C4, Ep1&C4, C1&Ep4, and Ep1&Ep4
!
! To see a printout of all 1-1, 1-2, 1-3 and 1-4 interactions set verbose = 1.
! These interactions are masked out of nonbonds. Thus the amber mask list is
! rebuilt from these 1-1, 1-2, 1-3 and 1-4 pairs. Pairs that aren't in the
! union of these are not masked.
!
! A separate list of 1-4 nonbonds is then compiled. This list does not agree
! in general with the above 1-4, since a 1-4 could also be a 1-3 if its
! in a ring. The rules in ephi() are used to see who is included:
!
! Here is that code:
!
!             DO 700 JN = 1,MAXLEN
!               I3 = IP(JN + IST)
!               K3T = KP(JN + IST)
!               L3T = LP(JN + IST)
!               IC0 = ICP(JN + IST)
!               IDUMI = ISIGN(1,K3T)
!               IDUML = ISIGN(1,L3T)
!               KDIV = (2 + IDUMI + IDUML) / 4
!               L3 = IABS(L3T)
!               FMULN = FLOAT(KDIV) * FMN(IC0)
!   C
!               II = (I3 + 3) / 3
!               JJ = (L3 + 3) / 3
!               IA1 = IAC(II)
!               IA2 = IAC(JJ)
!               IBIG = MAX0(IA1,IA2)
!               ISML = MIN0(IA1,IA2)
!               IC = IBIG * (IBIG - 1) / 2 + ISML
!   C
!   C             ----- CALCULATE THE 14-EEL ENERGY -----
!   C
!               R2 = FMULN / CT(JN)
!               R1 = SQRT(R2)
!       ...........
!
! So a pair is included in the 1-4 list if kdiv .gt. 0 and FMN(ic0) .gt. 0.
! This is decided at startup. This decision logic is applied to the parent
! atoms, and if they are included, so are the extra points attached:
!
! That is, in the above situation, if C1 and C4 pass the test then
! C1&C4, Ep1&C4, C1&Ep4, and Ep1&Ep4 are included. The dihedrals involving the
! extra points are not tested since the decision is based solely on parent
! atoms.
!
! The list of 1-4 nonbonds is also spit out if verbose = 1.
!
! To scale 1-4 charge-dipole and dipole-dipole interactions the same as
! 1-4 charge-charge (i.e. divided by scee) set scaldip = 1 (DEFAULT).
! If scaldip .eq. 0 the 1-4 charge-dipole and dipole-dipole interactions
! are treated the same as other dipolar interactions (i.e. divided by 1)
!
!-----------------------------------------------------------------------------
! Using the bond list (the one not involving hydrogens), find the number
! of neighbors (heavy atom, hydrogens and extra points) attached to each
! atom. For example if atom i is heavy, numnghbr(1,i) is the number of heavy
! atoms attached to atom i, while numnghbr(2,i) is the number of hydrogens
! attached, and numnghbr(3,i) is the number of extra points attached to i.
! The identities of neighbors are packed back to back in nghbrlst. The attyp
! array holds the atom types, usded to distinguish extra points or lone pairs
! from regular atoms.
!-----------------------------------------------------------------------------

  integer, parameter    :: extra_pnts_nb14_int_cnt = 2

  integer                     gbl_nb14_cnt, gbl_frame_cnt

  common / extra_pnts_nb14_int / gbl_nb14_cnt, gbl_frame_cnt

  ! Stuff below here used for 1-4 nonbonded interactions:

  integer, save                                 :: cit_nb14_cnt = 0

  integer,              allocatable, save       :: gbl_nb14(:,:)

  ! This is allocated / reallocated as needed in nb14_setup():

  integer,              allocatable, save       :: cit_nb14(:,:)

  ! Nonbonded 1-4 interactions possibly increase by:
  integer, parameter                            :: nb14_mult_fac = 9

  ! Extra Points Frame datatype:

  type ep_frame_rec
    sequence
    integer             :: extra_pnt(2)
    integer             :: ep_cnt
    integer             :: type
    integer             :: parent_atm
    integer             :: frame_atm1
    integer             :: frame_atm2
    integer             :: frame_atm3
  end type ep_frame_rec

  integer, parameter    :: ep_frame_rec_ints  = 8 ! don't use for allocation!

  type(ep_frame_rec), parameter :: null_ep_frame_rec = &
                                   ep_frame_rec(2*0,0,0,0,0,0,0)

#ifdef _OPENMP_
  GBFloat, allocatable                      :: ee14_arr(:)
  GBFloat, allocatable                      :: enb14_arr(:)
  GBFloat, allocatable                      :: e14vir11(:)
  GBFloat, allocatable                      :: e14vir12(:)
  GBFloat, allocatable                      :: e14vir13(:)
  GBFloat, allocatable                      :: e14vir22(:)
  GBFloat, allocatable                      :: e14vir23(:)
  GBFloat, allocatable                      :: e14vir33(:)

  integer(kind=omp_lock_kind), allocatable, private  :: omplk(:)
#endif /*_OPENMP_*/

  ! Stuff below here used only for extra points:

  type(ep_frame_rec),   allocatable, save       :: ep_frames(:)
  integer,              allocatable, save       :: custom_ep(:)
  double precision,     allocatable, save       :: ep_lcl_crd(:,:,:)
#ifdef MPI
  integer, save                                 :: my_ep_frame_cnt
  integer, allocatable, save                    :: gbl_my_ep_frame_lst(:)
#endif /* MPI */

  character(len=4), parameter                   :: ep_symbl = 'EP  '
  integer, private                              :: maxa

  private       get_nghbrs, &
                define_frames, &
                fill_bonded, &
                redo_masked, &
                build_nb14, &
                copy_nb14, &
                trim_bonds, &
                trim_theta, &
                trim_phi, &
                do_bond_pairs, &
                do_angle_pairs, &
                do_dihed_pairs, &
                add_one_list_iblo, &
                add_one_list_inb, &
                do_14pairs

contains

#ifdef MPI

!*******************************************************************************
!
! Subroutine:  nb14_setup_midpoint
!
! Description:  Handle workload subdivision for nb14 list.  This is also called
!               once in uniprocessor code, but this is really not necessary -
!               just the convention borrowed from bond-angle-dihedral code. One
!               advantage, though, is that gbl_nb14 can then be deallocated in
!               uniprocessor code, and it is probably bigger than it needs to
!               be.
!*******************************************************************************
subroutine nb14_setup_midpoint(num_ints, num_reals, use_atm_map)

  use parallel_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
  use mdin_ctrl_dat_mod
  use processor_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals
  integer                       :: use_atm_map(natom)
!size if really big, can we reduce it?

! Local variables:

  integer               :: alloc_failed
  integer               :: nb14_copy(3, gbl_nb14_cnt)
  integer               :: atm_i, atm_j, nb14_idx

! This routine can handle reallocation, and thus can be called multiple
! times.

! Find all diheds for which this process owns either atom:

  cit_nb14_cnt = 0

  do nb14_idx = 1, gbl_nb14_cnt

    atm_i = gbl_nb14(1, nb14_idx)
    atm_j = gbl_nb14(2, nb14_idx)
!    nb14_parm_idx = gbl_nb14(3, nb14_idx)

    if(proc_atm_space(atm_i) /= 0 .and. proc_atm_space(atm_j) /= 0) then
      cit_nb14_cnt = cit_nb14_cnt + 1
      nb14_copy(1, cit_nb14_cnt) = proc_atm_space(gbl_nb14(1, nb14_idx))
      nb14_copy(2, cit_nb14_cnt) = proc_atm_space(gbl_nb14(2, nb14_idx))
      nb14_copy(3, cit_nb14_cnt) = gbl_nb14(3, nb14_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
    end if
  end do

  if (cit_nb14_cnt .gt. 0) then
    if (allocated(cit_nb14)) then
      if (size(cit_nb14) .lt. cit_nb14_cnt * 3) then
        num_ints = num_ints - size(cit_nb14) * 3
        deallocate(cit_nb14)
        allocate(cit_nb14(3, 2*cit_nb14_cnt), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
        num_ints = num_ints + size(cit_nb14) * 3 * 2
      end if
    else
      allocate(cit_nb14(3, cit_nb14_cnt * 2), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(cit_nb14) * 3 * 2
    end if
    cit_nb14(:, 1:cit_nb14_cnt) = nb14_copy(:, 1:cit_nb14_cnt)
  else
    if(.not. allocated(cit_nb14)) allocate(cit_nb14(3,2*natom/numtasks))
    !this is required since cit_nb14 is used in another rouitine
  end if

  return

end subroutine nb14_setup_midpoint

!*******************************************************************************
!
! Subroutine:   build_cit_nb14_midpoint
!
!*******************************************************************************

subroutine build_cit_nb14_midpoint(dihed_cnt_a, dihed_a, dihed_cnt_h, dihed_h)

  use parallel_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
  use mdin_ctrl_dat_mod
  use processor_mod

  implicit none

! Formal variables:
  integer                       :: dihed_cnt_a
  integer                       :: dihed_cnt_h
  type(dihed_rec)               :: dihed_a(*)
  type(dihed_rec)               :: dihed_h(*)


! Local variables:

  integer               :: nb14_copy(3, gbl_nb14_cnt)
  integer               :: nb14_idx
  double precision      :: xcrd, ycrd, zcrd

  cit_nb14_cnt = 0
  cit_nb14(:,:)=0

  do nb14_idx = 1, dihed_cnt_a
    xcrd = proc_atm_crd(1,dihed_a(nb14_idx)%atm_j)
    ycrd = proc_atm_crd(2,dihed_a(nb14_idx)%atm_j)
    zcrd = proc_atm_crd(3,dihed_a(nb14_idx)%atm_j)
    if (dihed_a(nb14_idx)%atm_k .gt. 0 .and. &
        dihed_a(nb14_idx)%atm_l .gt. 0 .and. &
        gbl_fmn(dihed_a(nb14_idx)%parm_idx) .ne. ZERO) then
      if(xcrd .lt. proc_max_x_crd .and. xcrd .ge. proc_min_x_crd .and. &
         ycrd .lt. proc_max_y_crd .and. ycrd .ge. proc_min_y_crd .and. &
         zcrd .lt. proc_max_z_crd .and. zcrd .ge. proc_min_z_crd) then
        cit_nb14_cnt = cit_nb14_cnt + 1
        nb14_copy(1, cit_nb14_cnt) = dihed_a(nb14_idx)%atm_i
        nb14_copy(2, cit_nb14_cnt) = dihed_a(nb14_idx)%atm_l
        nb14_copy(3, cit_nb14_cnt) = dihed_a(nb14_idx)%parm_idx
        index_14(1:3,cit_nb14_cnt) = mult_veca_dihed(1:3,(nb14_idx-1)*3+3)
      end if
    end if
  end do

  do nb14_idx = 1, dihed_cnt_h
    xcrd = proc_atm_crd(1,dihed_h(nb14_idx)%atm_j)
    ycrd = proc_atm_crd(2,dihed_h(nb14_idx)%atm_j)
    zcrd = proc_atm_crd(3,dihed_h(nb14_idx)%atm_j)
    if (dihed_h(nb14_idx)%atm_k .gt. 0 .and. &
        dihed_h(nb14_idx)%atm_l .gt. 0 .and. &
        gbl_fmn(dihed_h(nb14_idx)%parm_idx) .ne. ZERO) then
      if(xcrd .lt. proc_max_x_crd .and. xcrd .ge. proc_min_x_crd .and. &
         ycrd .lt. proc_max_y_crd .and. ycrd .ge. proc_min_y_crd .and. &
         zcrd .lt. proc_max_z_crd .and. zcrd .ge. proc_min_z_crd) then
        cit_nb14_cnt = cit_nb14_cnt+1
        nb14_copy(1, cit_nb14_cnt) = dihed_h(nb14_idx)%atm_i
        nb14_copy(2, cit_nb14_cnt) = dihed_h(nb14_idx)%atm_l
        nb14_copy(3, cit_nb14_cnt) = dihed_h(nb14_idx)%parm_idx
        index_14(1:3,cit_nb14_cnt) = mult_vech_dihed(1:3,(nb14_idx-1)*3+3)
     end if
    end if
  end do

  if (cit_nb14_cnt .gt. 0) then
    cit_nb14(:, 1:cit_nb14_cnt) = nb14_copy(:, 1:cit_nb14_cnt)
  end if

  return

end subroutine build_cit_nb14_midpoint
#endif

!*******************************************************************************
!
! Subroutine:  nb14_setup
!
! Description:  Handle workload subdivision for nb14 list.  This is also called
!               once in uniprocessor code, but this is really not necessary -
!               just the convention borrowed from bond-angle-dihedral code. One
!               advantage, though, is that gbl_nb14 can then be deallocated in
!               uniprocessor code, and it is probably bigger than it needs to
!               be.
!*******************************************************************************

subroutine nb14_setup(num_ints, num_reals, use_atm_map)

  use parallel_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
!  use gamd_mod,only :: igamd
#ifdef _OPENMP_
  use mdin_ctrl_dat_mod
#endif

#ifdef CUDA
  use mdin_ctrl_dat_mod
  use charmm_mod, only : charmm_active
#ifdef GTI
  use ti_mod
  use reaf_mod
#endif
#endif
  
  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals
  integer                       :: use_atm_map(natom)

! Local variables:

  integer               :: alloc_failed
  integer               :: nb14_copy(3, gbl_nb14_cnt)
  integer               :: atm_i, atm_j, nb14_idx
#ifdef _OPENMP_
  integer               :: i
#endif

#if defined(GTI) && defined (CUDA)
  integer, allocatable::iac_ti(:)
#endif

! This routine can handle reallocation, and thus can be called multiple
! times.

! Find all diheds for which this process owns either atom:
  
  cit_nb14_cnt = 0

  do nb14_idx = 1, gbl_nb14_cnt

    atm_i = gbl_nb14(1, nb14_idx)
    atm_j = gbl_nb14(2, nb14_idx)
!    nb14_parm_idx = gbl_nb14(3, nb14_idx)

#if defined(MPI) && !defined(CUDA)
#ifdef _OPENMP_
 if(using_gb_potential) then
   if(master) then
      cit_nb14_cnt = cit_nb14_cnt + 1
      nb14_copy(:, cit_nb14_cnt) = gbl_nb14(:, nb14_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
   end if
 else
#endif
    if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
      cit_nb14_cnt = cit_nb14_cnt + 1
      nb14_copy(:, cit_nb14_cnt) = gbl_nb14(:, nb14_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
    else if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
      cit_nb14_cnt = cit_nb14_cnt + 1
      nb14_copy(:, cit_nb14_cnt) = gbl_nb14(:, nb14_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
    end if
#ifdef _OPENMP_
 end if
#endif
#else
    cit_nb14_cnt = cit_nb14_cnt + 1
    nb14_copy(:, cit_nb14_cnt) = gbl_nb14(:, nb14_idx)
    use_atm_map(atm_i) = 1
    use_atm_map(atm_j) = 1
#endif

  end do

  if (cit_nb14_cnt .gt. 0) then
    if (allocated(cit_nb14)) then
      if (size(cit_nb14) .lt. cit_nb14_cnt * 3) then
        num_ints = num_ints - size(cit_nb14) * 3
        deallocate(cit_nb14)
        allocate(cit_nb14(3, cit_nb14_cnt), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
        num_ints = num_ints + size(cit_nb14) * 3
      end if
    else
      allocate(cit_nb14(3, cit_nb14_cnt), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(cit_nb14) * 3
    end if
    cit_nb14(:, 1:cit_nb14_cnt) = nb14_copy(:, 1:cit_nb14_cnt)
  end if

#ifdef CUDA

#ifdef GTI        
        allocate(iac_ti(natom))
        iac_ti(1:natom)=atm_iac(1:natom)

        if (charmm_active) then
            call gti_nb_setup( ntypes, iac_ti, typ_ico, &
                gbl_cn114, gbl_cn214, gbl_cn6, gbl_cn7, gbl_cn8, C4Pairwise) ! C4PairwiseCUDA
        else
            call gti_nb_setup( ntypes, iac_ti, typ_ico, &
                gbl_cn1, gbl_cn2, gbl_cn6, gbl_cn7, gbl_cn8, C4Pairwise) ! C4PairwiseCUDA
        end if 

        if (reaf_mode .ge. 0) then 
          call gti_reaf_nb_setup(ntypes, iac_ti, typ_ico, &
                cit_nb14_cnt, cit_nb14, gbl_one_scee, gbl_one_scnb, &
                gbl_cn1, gbl_cn2)
        endif

        if (ti_mode .ne. 0) then   
    
          !! Add TI-vdw mask here
          do atm_i=1, natom
              if (ti_vdw_mask_list(atm_i).ne.0) iac_ti(atm_i)=0
          end do

          if(.not.(igamd.ge.6.and.igamd.le.11))then
            if (charmm_active) then
                call gti_ti_nb_setup(ntypes, iac_ti, typ_ico, &
                    cit_nb14_cnt, cit_nb14, gbl_one_scee, gbl_one_scnb, &
                    gbl_cn114, gbl_cn214)
            else
                call gti_ti_nb_setup(ntypes, iac_ti, typ_ico, &
                    cit_nb14_cnt, cit_nb14, gbl_one_scee, gbl_one_scnb, &
                    gbl_cn1, gbl_cn2)
            end if
          else
            if (charmm_active) then
                call gti_ti_nb_gamd_setup(ntypes, iac_ti, typ_ico, &
                    cit_nb14_cnt, cit_nb14, gbl_one_scee, gbl_one_scnb, &
                    gbl_cn114, gbl_cn214)
            else
                call gti_ti_nb_gamd_setup(ntypes, iac_ti, typ_ico, &
                     cit_nb14_cnt, cit_nb14, gbl_one_scee, gbl_one_scnb, &
                     gbl_cn1, gbl_cn2)
            end if
           endif

           if (reaf_mode .lt. 0) then
            do atm_i=1, natom
              if (ti_lst(1,atm_i).ne.0  .or. ti_lst(2,atm_i).ne.0 ) iac_ti(atm_i)=0
            end do
           else
              do atm_i=1, natom
                if (ti_lst(2-reaf_mode,atm_i).ne.0 ) iac_ti(atm_i)=0
              end do
            endif

         end if

        if (lj1264 .ne. 0 .and. plj1264 .ne. 0) then 
            call gti_lj1264plj1264_nb_setup(atm_isymbl, ntypes, iac_ti, typ_ico, gbl_cn6, gbl_cn7, gbl_cn8, C4Pairwise);
            if (allocated(atm_isymbl)) deallocate(atm_isymbl)
        endif
        if (lj1264 .eq. 0 .and. plj1264 .ne. 0) then ! C4PairwiseCUDA
            call gti_plj1264_nb_setup(atm_isymbl, ntypes, iac_ti, typ_ico, gbl_cn7, gbl_cn8, C4Pairwise);
            if (allocated(atm_isymbl)) deallocate(atm_isymbl)
        endif
        if (lj1264 .ne. 0 .and. plj1264 .eq. 0) then ! C4PairwiseCUDA
            call gti_lj1264_nb_setup(atm_isymbl, ntypes, iac_ti, typ_ico, gbl_cn6);
            if (allocated(atm_isymbl)) deallocate(atm_isymbl)
        endif

       if(igamd.eq.12.or.igamd.eq.13.or.igamd.eq.14.or.igamd.eq.15.or.igamd.eq.16.or.igamd.eq.17 &
          .or.(igamd.ge.18.and.igamd.le.28).or.(igamd.ge.110.and.igamd.le.120))then
            if (charmm_active) then
                call gti_ti_nb_gamd_setup(ntypes, iac_ti, typ_ico, &
                    cit_nb14_cnt, cit_nb14, gbl_one_scee, gbl_one_scnb, &
                    gbl_cn114, gbl_cn214)
            else
                call gti_ti_nb_gamd_setup(ntypes, iac_ti, typ_ico, &
                    cit_nb14_cnt, cit_nb14, gbl_one_scee, gbl_one_scnb, &
                    gbl_cn1, gbl_cn2)
            end if
       endif
        

                    
        if (charmm_active) then
            call gpu_nb14_setup(cit_nb14_cnt, cit_nb14, gbl_one_scee, gbl_one_scnb, ntypes, lj1264, plj1264, iac_ti, typ_ico, & 
            gbl_cn1, gbl_cn2, gbl_cn6, gbl_cn7, gbl_cn8, C4Pairwise, gbl_cn114, gbl_cn214)  ! C4PairwiseCUDA
        else
            call gpu_nb14_setup(cit_nb14_cnt, cit_nb14, gbl_one_scee, gbl_one_scnb, ntypes, lj1264, plj1264, iac_ti, typ_ico, &
            gbl_cn1, gbl_cn2, gbl_cn6, gbl_cn7, gbl_cn8, C4Pairwise, gbl_cn1, gbl_cn2) ! C4PairwiseCUDA
        end if 
        
        deallocate (iac_ti)
#else
  if (charmm_active) then
    call gpu_nb14_setup(cit_nb14_cnt, cit_nb14, gbl_one_scee, gbl_one_scnb, ntypes, lj1264, plj1264, &
                        atm_iac, typ_ico, gbl_cn1, gbl_cn2, gbl_cn6, gbl_cn7, gbl_cn8, C4Pairwise, gbl_cn114, gbl_cn214) ! C4PairwiseCUDA
  else
    call gpu_nb14_setup(cit_nb14_cnt, cit_nb14, gbl_one_scee, gbl_one_scnb, ntypes, lj1264, plj1264, &
                        atm_iac, typ_ico, gbl_cn1, gbl_cn2, gbl_cn6, gbl_cn7, gbl_cn8, C4Pairwise, gbl_cn1, gbl_cn2) ! C4PairwiseCUDA
  end if
#endif 

#endif

#ifdef _OPENMP_
if(using_gb_potential) then
   if (flag_nb14 .eqv. .true.) then
     flag_nb14 = .false. ! 
     allocate(ee14_arr(natom), &
              enb14_arr(natom), &
              e14vir11(natom), &
              e14vir12(natom), &
              e14vir13(natom), &
              e14vir22(natom), &
              e14vir23(natom), &
              e14vir33(natom), &
              omplk(natom), &
           stat = alloc_failed)
     if (alloc_failed .ne. 0) call setup_alloc_error

     do i = 1, natom
       call omp_init_lock(omplk(i))
     end do
   end if
end if
#endif /*_OPENMP_*/
  return

end subroutine nb14_setup

!*******************************************************************************!
! Subroutine:  alloc_nb14_mem_only
!
! Description: <TBS>
!
!*******************************************************************************

subroutine alloc_nb14_mem_only(max_nb14, num_ints, num_reals)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: max_nb14

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed

  allocate(gbl_nb14(3, max_nb14), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_ints = num_ints + size(gbl_nb14)

  gbl_nb14(:,:) = 0
  
  return

end subroutine alloc_nb14_mem_only 

!*******************************************************************************!
! Subroutine:  alloc_extra_pnts_nb14_mem
!
! Description: <TBS>
!
!*******************************************************************************

subroutine alloc_extra_pnts_nb14_mem(max_frames, max_nb14, num_ints, num_reals)

  use pmemd_lib_mod
  use prmtop_dat_mod, only : natom
  implicit none

! Formal arguments:

  integer, intent(in)           :: max_frames
  integer, intent(in)           :: max_nb14

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed
  if (allocated(ep_frames)) deallocate(ep_frames)
  if (allocated(custom_ep)) deallocate(custom_ep)
  allocate(ep_frames(max_frames), custom_ep(natom), ep_lcl_crd(3, 2, max_frames), &
           gbl_nb14(3, max_nb14), stat = alloc_failed)
  
  if (alloc_failed .ne. 0) call setup_alloc_error

  num_ints = num_ints + natom + size(ep_frames) * ep_frame_rec_ints + size(gbl_nb14)

  num_reals = num_reals + size(ep_lcl_crd)

  ep_frames(:) = null_ep_frame_rec
  gbl_nb14(:,:) = 0
  ep_lcl_crd(:,:,:) = ZERO 
  
  return

end subroutine alloc_extra_pnts_nb14_mem 

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  bcast_extra_pnts_nb14_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine bcast_extra_pnts_nb14_dat

  use parallel_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod

  implicit none

! Local variables:

  integer       :: num_ints, num_reals  ! returned values discarded
  integer       :: bytes_per_unit
  integer       :: alloc_failed

  call mpi_bcast(gbl_nb14_cnt, extra_pnts_nb14_int_cnt, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  num_ints = 0
  num_reals = 0
  if (.not. master) then

    if (gbl_frame_cnt .eq. 0) then
      call alloc_nb14_mem_only(gbl_nb14_cnt, num_ints, num_reals)
    else
      call alloc_extra_pnts_nb14_mem(numextra, gbl_nb14_cnt, &
                                     num_ints, num_reals)
    end if
  end if

  ! All callers need the nb14 stuff...

  call mpi_bcast(gbl_nb14, gbl_nb14_cnt * 3, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)


  if (gbl_frame_cnt .ne. 0) then

    call get_bytesize(ep_frames(1), ep_frames(2), bytes_per_unit)

    call mpi_bcast(ep_frames, gbl_frame_cnt * bytes_per_unit, mpi_byte, 0, &
                   pmemd_comm, err_code_mpi)
    call mpi_bcast(ep_lcl_crd, 3 * 2 * gbl_frame_cnt, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    
    ! Here we allocate the process-local frame list in all tasks, including
    ! the master:

    allocate(gbl_my_ep_frame_lst(gbl_frame_cnt), stat = alloc_failed)
    if (alloc_failed .ne. 0) call setup_alloc_error
    num_ints = num_ints + size(gbl_my_ep_frame_lst)
    gbl_my_ep_frame_lst(:) = 0
    my_ep_frame_cnt = 0 ! until first atom division...

  end if

  return

end subroutine bcast_extra_pnts_nb14_dat
#endif

!*******************************************************************************!
! Subroutine:  init_nb14_only
!
! Description: <TBS>
!
!*******************************************************************************

subroutine init_nb14_only(num_ints, num_reals)

  use prmtop_dat_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables

  integer                       :: dihed_idx

  call alloc_nb14_mem_only(2 * (nphih + nphia), num_ints, num_reals)
   
  gbl_nb14_cnt = 0
  gbl_frame_cnt = 0

  do dihed_idx = 1, nphih + nphia

    if (gbl_dihed(dihed_idx)%atm_k .gt. 0 .and. &
        gbl_dihed(dihed_idx)%atm_l .gt. 0 .and. &
        gbl_fmn(gbl_dihed(dihed_idx)%parm_idx) .ne. ZERO) then
        gbl_nb14_cnt = gbl_nb14_cnt + 1
        gbl_nb14(1, gbl_nb14_cnt) = gbl_dihed(dihed_idx)%atm_i
        gbl_nb14(2, gbl_nb14_cnt) = gbl_dihed(dihed_idx)%atm_l
        gbl_nb14(3, gbl_nb14_cnt) = gbl_dihed(dihed_idx)%parm_idx
    end if
  end do  !  n = 1, nphih + nphia

  return

end subroutine init_nb14_only 

!*******************************************************************************!
! Subroutine:  init_extra_pnts_nb14
!
! Description: <TBS>
!
!*******************************************************************************

subroutine init_extra_pnts_nb14(num_ints, num_reals)

  use mdin_ewald_dat_mod
  use mdin_ctrl_dat_mod, only : using_gb_potential
  use charmm_mod, only : charmm_active
  use prmtop_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer, allocatable  :: epbtyp(:,:)
  integer, allocatable  :: nghbrs(:,:)
  integer, allocatable  :: hnghbrs(:,:)
  integer, allocatable  :: enghbrs(:,:)
  integer, allocatable  :: numnghbr(:,:)
  integer, allocatable  :: epowner(:)
  integer, allocatable  :: offset(:)
  integer, allocatable  :: test(:)
  integer, allocatable  :: i11(:,:)
  integer, allocatable  :: i12(:,:)
  integer, allocatable  :: i13(:,:)
  integer, allocatable  :: i14(:,:)
  integer, allocatable  :: nb_14_list(:,:)
  integer, allocatable  :: s3(:)
   
  integer               :: n
  integer               :: numextra_test
  integer               :: num11, num12, num13, num14
  integer               :: max11, max12, max13, max14
  integer               :: alloc_failed
  integer               :: i
   
  gbl_nb14_cnt = 0
  gbl_frame_cnt = 0

  call alloc_extra_pnts_nb14_mem(numextra, nb14_mult_fac * (nphih + nphia), &
                                 num_ints, num_reals)
   
  numextra_test = 0
  do n = 1, natom
    if (atm_isymbl(n) .eq. ep_symbl) numextra_test = numextra_test + 1
  end do

  if (numextra_test .ne. numextra) then
    write(mdout, *) 'Error in numextra_test'
    call mexit(mdout, 1)
  end if

  if (numextra .eq. 0) then
    frameon = 0
  end if
 
  max11 = natom + numextra
  max12 = 3 * (nbonh + nbona + 8*numextra)
  max13 = 3 * (ntheth + ntheta + 32*numextra)
  max14 = nb14_mult_fac * (nphih + nphia + 64*numextra)
  maxa = max(max11, max12, max13, max14)

  allocate(epbtyp(5, natom), &
           nghbrs(5, natom), &
           hnghbrs(5, natom), &
           enghbrs(5, natom), &
           numnghbr(3, natom), &
           epowner(natom), &
           offset(natom), &
           test(natom), &
           i11(2, max11), &
           i12(2, max12), &
           i13(2, max13), &
           i14(2, max14), &
           nb_14_list(3, max14), &
           s3(maxa), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  nghbrs(:,:) = 0
  hnghbrs(:,:) = 0
  enghbrs(:,:) = 0
  numnghbr(:,:) = 0
  epowner(:) = 0
  epbtyp(:,:) = 0

  i11(:,:) = 0
  i12(:,:) = 0
  i13(:,:) = 0
  i14(:,:) = 0
  nb_14_list(:,:) = 0

  call get_nghbrs(nghbrs, hnghbrs, enghbrs, numnghbr, epowner, epbtyp)

  call define_frames(natom, atm_isymbl, nghbrs, hnghbrs, enghbrs, &
                     numnghbr, ep_frames, ep_lcl_crd, &
                     gbl_frame_cnt, epbtyp, gbl_req, verbose)

  ! Activate the frame orientation of extra points if the count is nonzero and
  ! there is no box--assume that frames are needed to place these particles,
  ! even though the frameon flag may not have been set already.
  if (gbl_frame_cnt .gt. 0 .and. using_gb_potential .and. .not. charmm_active) then
    frameon = 1
    chngmask = 1
  end if
  
  call fill_bonded(max11, max12, max13, max14, num11, num12, num13, num14, &
                   i11, i12, i13, i14, enghbrs, numnghbr, epowner, verbose)

  if (chngmask .eq. 1) then
    call redo_masked(natom, atm_numex, gbl_natex, next, &
                     num11, num12, num13, num14, &
                     i11, i12, i13, i14, offset, test)
  end if

  ! Use nb_14_list and then copy to permanent:

  call build_nb14(nb_14_list, gbl_nb14_cnt, max14, epowner, numnghbr, enghbrs, &
                  natom, chngmask, verbose)
  
  call copy_nb14(nb_14_list, gbl_nb14, gbl_nb14_cnt)
  
  ! if frameon = 1 use frames and forces to define ep's;
  ! else use their masses and amber force params
   
  if (frameon .eq. 1) then
      
    ! Zero out mass and massinv for extra points:
      
    call fix_masses(natom, atm_mass, epowner)
      
    ! Now remove bonds etc involving extra points:
      
    call trim_bonds(nbonh, gbl_bond, epowner)
    call trim_bonds(nbona, gbl_bond(bonda_idx), epowner)

    ! Make the bond arrays sequential for shake and force routines:

    do i = 1, nbona
      gbl_bond(nbonh + i) = gbl_bond(bonda_idx + i - 1)
    end do

    bonda_idx = nbonh + 1

    call trim_theta(ntheth, gbl_angle, epowner)
    call trim_theta(ntheta, gbl_angle(anglea_idx), epowner)

    ! Make the angle arrays sequential:

    do i = 1, ntheta
      gbl_angle(ntheth + i) = gbl_angle(anglea_idx + i - 1)
    end do

    anglea_idx = ntheth + 1

    call trim_phi(nphih, gbl_dihed, epowner)
    call trim_phi(nphia, gbl_dihed(diheda_idx), epowner)

    do i = 1, nphia
      gbl_dihed(nphih + i) = gbl_dihed(diheda_idx + i - 1)
    end do

    diheda_idx = nphih + 1

  end if
   
  deallocate(epbtyp, nghbrs, hnghbrs, enghbrs, numnghbr, epowner, offset, &
             test, i11, i12, i13, i14, nb_14_list, s3)

  return

end subroutine init_extra_pnts_nb14 

!----------------------------------------------------------------------------------------------
! get_nghbrs: find atoms connected to one another, sorted into heavy, hydrogen, and EP classes.
!              
! Arguments:
!   ngbhrs:    two-dimensional integer array listing all of the neighbors of any particular
!              atom.  ngbhrs(i,k) gives the topological atom number of the ith neighbor of the
!              kth atom.
!   hnghbrs:   same as nghbrs, but all hydrogens in the topology are listed in this array.
!              Both arrays are obviously mostly blank.
!   enghbrs:   same as nghbrs, but all extra points in the topology are listed in this array.
!   numnghbr:  the number of heavy, hydrogen, and EP neighbors for each atom to guide indexing
!              into the arrays above.
!   epowner:   the owner of an extra point.  epowner(j) will state the topological number of
!              the atom that owns the extra point, or 0 otherwise.
!   epbtyp:    two dimensional array listing bond types connecting an atom to neighboring EPs
!----------------------------------------------------------------------------------------------
subroutine get_nghbrs(nghbrs, hnghbrs, enghbrs, numnghbr, epowner, epbtyp)

  use prmtop_dat_mod, only     : gbl_bond, bonda_idx, nbona, nbonh, atm_isymbl, &
                                 num_custom_ep, prm_ep_frames
  use file_io_dat_mod, only    : mdout
  use gbl_constants_mod, only  : error_hdr
  use pmemd_lib_mod, only      : mexit
  
  implicit none

  ! Formal arguments:
  integer               :: nghbrs(5, *)
  integer               :: hnghbrs(5, *)
  integer               :: enghbrs(5, *)
  integer               :: numnghbr(3, *)
  integer               :: epowner(*)
  integer               :: epbtyp(5, *)

  ! Local variables:
  integer               :: bond_idx, atm_i, atm_j, i
  logical               :: found
  
  do bond_idx = bonda_idx, bonda_idx + nbona - 1

    atm_i = gbl_bond(bond_idx)%atm_i
    atm_j = gbl_bond(bond_idx)%atm_j

    if (atm_isymbl(atm_i) .eq. ep_symbl) then
      numnghbr(3, atm_j) = numnghbr(3, atm_j) + 1
      enghbrs(numnghbr(3, atm_j), atm_j) = atm_i
      epowner(atm_i) = atm_j
      epbtyp(numnghbr(3, atm_j), atm_j) = gbl_bond(bond_idx)%parm_idx
    else if (atm_isymbl(atm_j) .eq. ep_symbl) then
      numnghbr(3, atm_i) = numnghbr(3, atm_i) + 1
      enghbrs(numnghbr(3, atm_i), atm_i) = atm_j
      epowner(atm_j) = atm_i
      epbtyp(numnghbr(3, atm_i), atm_i) = gbl_bond(bond_idx)%parm_idx
    else
      numnghbr(1, atm_i) = numnghbr(1, atm_i) + 1
      numnghbr(1, atm_j) = numnghbr(1, atm_j) + 1
      nghbrs(numnghbr(1, atm_i), atm_i) = atm_j
      nghbrs(numnghbr(1, atm_j), atm_j) = atm_i
    end if

  end do

  do bond_idx = 1, nbonh

    atm_i = gbl_bond(bond_idx)%atm_i
    atm_j = gbl_bond(bond_idx)%atm_j

    numnghbr(2, atm_i) = numnghbr(2, atm_i) + 1
    numnghbr(2, atm_j) = numnghbr(2, atm_j) + 1
    hnghbrs(numnghbr(2, atm_i), atm_i) = atm_j
    hnghbrs(numnghbr(2, atm_j), atm_j) = atm_i

  end do

  ! Custom extra points (virtual sites)
  do bond_idx = 1, num_custom_ep
    atm_i = prm_ep_frames(bond_idx)%extra_pnt
    atm_j = prm_ep_frames(bond_idx)%parent_atm

    ! Check that the neighbor relationship is not already logged
    found = .false.
    do i = 1, numnghbr(1, atm_i)
      if (nghbrs(i, atm_i) .eq. atm_j) then
        found = .true.
      end if
    end do
    do i = 1, numnghbr(2, atm_i)
      if (hnghbrs(i, atm_i) .eq. atm_j) then
        found = .true.
      end if
    end do
    do i = 1, numnghbr(3, atm_i)
      if (enghbrs(i, atm_i) .eq. atm_j) then
        found = .true.
      end if
    end do
    if (found) then
      write(mdout, '(2A)') error_hdr, 'Custom virtual site placement would disrupt &
                                      &pre-existing bonds.'
      call mexit(mdout, 1)
    end if
    if (.not. found) then
      numnghbr(3, atm_j) = numnghbr(3, atm_j) + 1
      enghbrs(numnghbr(3, atm_j), atm_j) = atm_i
      epowner(atm_i) = atm_j
      epbtyp(numnghbr(3, atm_j), atm_j) = 0
    end if
  end do
  
  return

end subroutine get_nghbrs 

!*******************************************************************************!
! Internal Subroutine:  define_frames
!
! Description:  Fix the positions of EP in the local frame / coord depending
!               on the kind of frame and atom types.
!
!*******************************************************************************

subroutine define_frames(natom, isymbl, nghbrs, hnghbrs, enghbrs, numnghbr, &
                         frames, ep_lcl_crd, frame_cnt, epbtyp, req, verbose)

  use gbl_constants_mod, only : DEG_TO_RAD, info_hdr, error_hdr
  use pmemd_lib_mod
  use prmtop_dat_mod, only    : num_custom_ep, prm_ep_frames
  
  implicit none
  
! Formal arguments:

  integer               :: natom
  character(len=4)      :: isymbl(*)
  integer               :: nghbrs(5, *)
  integer               :: hnghbrs(5, *)
  integer               :: enghbrs(5, *)
  integer               :: numnghbr(3, *)
  type(ep_frame_rec)    :: frames(*)
  double precision      :: ep_lcl_crd(3, 2, *)
  integer               :: frame_cnt
  integer               :: epbtyp(5, *)
  double precision      :: req(*)
  integer               :: verbose

! Local variables:

  integer               :: k, m, n, l, ierr
  double precision      :: tetcos, tetsin, angle, scos, ssin
  character(len=4)      :: sulf, sulfh
  integer, allocatable  :: custom_ep_tbl(:,:)
  integer, allocatable  :: custom_ep_count(:)
  
  sulf = 'S   '
  sulfh = 'SH  '
   
  ! Get half-angle for tetrahedral:
   
  angle = 54.735d0
  angle = angle * DEG_TO_RAD
  tetcos = cos(angle)
  tetsin = sin(angle)

  ! Allocate a table for custom EP placement
  allocate(custom_ep_tbl(5, natom), custom_ep_count(natom), stat=ierr)
  if (ierr .ne. 0) then
    write(mdout, '(2A)') error_hdr, 'define_frames: custom EP table allocation failed'
    call mexit(mdout, 1)
  end if
  custom_ep_tbl(:,:) = 0
  custom_ep_count(:) = 0
  do n = 1, num_custom_ep
    k = prm_ep_frames(n)%parent_atm
    custom_ep_count(k) = custom_ep_count(k) + 1
    custom_ep_tbl(custom_ep_count(k),k) = prm_ep_frames(n)%extra_pnt
  end do
     
  ! Get cos, sin for 60:
  scos = HALF
  ssin = sqrt(ONE - scos * scos)

  ! Deal with custom EP frames first
  do n = 1, num_custom_ep
    frames(n)%parent_atm   = prm_ep_frames(n)%parent_atm
    frames(n)%ep_cnt       = 1
    frames(n)%extra_pnt(1) = prm_ep_frames(n)%extra_pnt
    frames(n)%extra_pnt(2) = 0 
    frames(n)%frame_atm1   = prm_ep_frames(n)%frame_atm1
    frames(n)%frame_atm2   = prm_ep_frames(n)%frame_atm2
    frames(n)%frame_atm3   = prm_ep_frames(n)%frame_atm3
    frames(n)%type         = prm_ep_frames(n)%type
  end do
  frame_cnt = num_custom_ep
  
  ! Deal with hard-wired frame types
  do n = 1, natom
      
    if (numnghbr(3, n) .gt. 0) then

      ! Trap to deal with custom EP frames
      if (custom_ep_count(n) .gt. 0) then
        cycle
      end if

      ! Hardwired frame definitions are limited in scope
      if (numnghbr(1, n) + numnghbr(2, n) .gt. 2) then
        write(mdout, *) 'EXTRA_PTS: too many nghbrs!!'
        call mexit(mdout, 1)
      end if

      frame_cnt = frame_cnt + 1
      frames(frame_cnt)%parent_atm = n
      frames(frame_cnt)%ep_cnt = numnghbr(3, n)

      frames(frame_cnt)%extra_pnt(:) = 0

      do k = 1, frames(frame_cnt)%ep_cnt
        frames(frame_cnt)%extra_pnt(k) = enghbrs(k, n)
      end do
         
      if (numnghbr(1, n) .eq. 0 .and. numnghbr(2, n) .eq. 2 .and. &
          numnghbr(3, n) .eq. 1) then
            
        ! TIP4P water:
        ! (temporarily assign type to 3, will set to 1 in next section)
            
        frames(frame_cnt)%type = 3
        frames(frame_cnt)%frame_atm1 = hnghbrs(1, n)
        frames(frame_cnt)%frame_atm2 = n
        frames(frame_cnt)%frame_atm3 = hnghbrs(2, n)
            
      else if (numnghbr(1, n) .eq. 0 .and. numnghbr(2, n) .eq. 2 .and. &
               numnghbr(3, n) .eq. 2) then
            
        ! TIP5P water:
            
        frames(frame_cnt)%type = 1
        frames(frame_cnt)%frame_atm1 = hnghbrs(1, n)
        frames(frame_cnt)%frame_atm2 = n
        frames(frame_cnt)%frame_atm3 = hnghbrs(2, n)
            
      else if (numnghbr(1, n) .gt. 1) then
            
        ! "ordinary" type of frame defined by two other heavy atoms:
            
        frames(frame_cnt)%type = 1
        frames(frame_cnt)%frame_atm1 = nghbrs(1, n)
        frames(frame_cnt)%frame_atm2 = n
        frames(frame_cnt)%frame_atm3 = nghbrs(2, n)
            
      else if (numnghbr(1, n) .eq. 1 .and. numnghbr(2, n) .eq. 1) then
            
        ! frame defined by one heavy atom and one hydrogen:
            
        frames(frame_cnt)%type = 1
        frames(frame_cnt)%frame_atm1 = nghbrs(1, n)
        frames(frame_cnt)%frame_atm2 = n
        frames(frame_cnt)%frame_atm3 = hnghbrs(1, n)
            
      else if (numnghbr(1, n) .eq. 1 .and. numnghbr(2, n) .eq. 0) then
            
        ! Assume this is CARBONYL oxygen:
        ! (Need to use midpoints of other two bonds of the carbon
        ! for frame_atm1 and frame_atm3 in orient force; thus (mis)use
        ! frame_atm1 frame_atm2 frame_atm3 and parent_atm to store 4 atoms for
        ! this special case.)
            
        frames(frame_cnt)%type = 2
        m = nghbrs(1, n)
        frames(frame_cnt)%frame_atm2 = m
        if (numnghbr(1, m) .ne. 3 .or.numnghbr(2, m) .gt. 0) then
          write(mdout, *) 'EXTRA_PTS: frame type 2 Should not be here'
          write(mdout, *) n, m, numnghbr(1, m), numnghbr(2, m)
          call mexit(mdout, 1)
        end if
            
        ! numnghbr(1, m) = 3. One is n (the oxygen) and other 2 are
        ! other bonding partners of carbon:

        frames(frame_cnt)%frame_atm1 = 0
        frames(frame_cnt)%frame_atm3 = 0

        k = 1
        do while (k .lt. 4 .and. frames(frame_cnt)%frame_atm1 .eq. 0)

          if (nghbrs(k, m) .ne. n) then
            frames(frame_cnt)%frame_atm1 = nghbrs(k, m)
          end if
          k = k + 1

        end do

        k = 1
        do while (k .lt. 4 .and. frames(frame_cnt)%frame_atm3 .eq. 0)

          if (nghbrs(k, m) .ne. n .and. &
              nghbrs(k, m) .ne. frames(frame_cnt)%frame_atm1) then
            frames(frame_cnt)%frame_atm3 = nghbrs(k, m)
          end if
          k = k + 1

        end do

        if (frames(frame_cnt)%frame_atm1.eq. 0 .or. &
            frames(frame_cnt)%frame_atm3 .eq. 0) then
          write(mdout, *) 'EXTRA_PTS: cannot find first or third frame point '
          write(mdout, *) 'define: ', n, numnghbr(1, n), numnghbr(2, n), &
            numnghbr(3, n), frames(frame_cnt)%frame_atm1, &
            frames(frame_cnt)%frame_atm3
          call mexit(mdout, 1)
        endif

      else

        write(mdout, *) 'EXTRA_PTS: unexpected numnghbr array: '
        write(mdout, *) 'define: ', n, numnghbr(1, n), numnghbr(2, n), &
          numnghbr(3, n)
        call mexit(mdout, 1)

      end if  !  ( numnghbr(1, n) .eq. 0  .and. numnghbr(2, n) .eq. 2

    end if  ! ( numnghbr(3, n) .gt. 0 )

  end do ! 1, natom
  
  ! Get the local coords:
   
  ep_lcl_crd(1:3, 1:2, 1:frame_cnt) = ZERO

  do n = 1, frame_cnt

    l = frames(n)%parent_atm

    if (frames(n)%type .eq. 1 .or. frames(n)%type .eq. 3) then
         
      ! Z axis is along symmetry axis of second atom opposite
      ! bisector of frame_atm1, frame_atm3;
      ! X axis is along the diff vector frame_atm3 minus frame_atm1;
      ! Y axis is cross product.
         
      if (frames(n)%ep_cnt .eq. 1) then
            
        ! Extra point is along the z-direction: positive for ordinary
        ! lone pair, negative for TIP4P water extra point:
            
        ep_lcl_crd(3, 1, n) = req(epbtyp(1, l))
        if (frames(n)%type .eq. 3) then
          ep_lcl_crd(3, 1, n) = -req(epbtyp(1, l))
          frames(n)%type = 1
        end if
            
      else if (frames(n)%ep_cnt .eq. 2) then
            
        ! Extra points are in the z, y plane, tetrahedrally
        ! (unless frame_atm2 atom is sulfur, in which case they
        ! are opposite along y):
            
        m = frames(n)%frame_atm2
        if (isymbl(m) .eq. sulf .or. isymbl(m) .eq. sulfh) then
          ep_lcl_crd(2, 1, n) = req(epbtyp(1, l))
          ep_lcl_crd(2, 2, n) = -req(epbtyp(2, l))
        else
          ep_lcl_crd(3, 1, n) = tetcos * req(epbtyp(1, l))
          ep_lcl_crd(2, 1, n) = tetsin * req(epbtyp(1, l))
          ep_lcl_crd(3, 2, n) = tetcos * req(epbtyp(2, l))
          ep_lcl_crd(2, 2, n) = -tetsin * req(epbtyp(2, l))
        end if
            
      else

        write(mdout, *) 'EXTRA_PTS: unexpected ep_cnt value: ', frames(n)%ep_cnt
        call mexit(mdout, 1)

      end if  ! ( frames(n)%ep_cnt .eq. 1 )
         
    else if (frames(n)%type .eq. 2) then
         
      ! Z axis is along bond from frame_atm2 to parent_atm;
      ! X axis in plane of parent_atm and midpoints of frame_atm1, frame_atm2
      ! and frame_atm2, frame_atm3.
         
      if (frames(n)%ep_cnt .eq. 1) then
        ep_lcl_crd(3, 1, n) = req(epbtyp(1, l))
      else if (frames(n)%ep_cnt .eq. 2) then
        ep_lcl_crd(3, 1, n) = scos * req(epbtyp(1, l))
        ep_lcl_crd(1, 1, n) = ssin * req(epbtyp(1, l))
        ep_lcl_crd(3, 2, n) = scos * req(epbtyp(2, l))
        ep_lcl_crd(1, 2, n) = -ssin * req(epbtyp(2, l))
      else
        write(mdout, *) 'EXTRA_PTS: unexpected ep_cnt value: ', frames(n)%ep_cnt
        call mexit(mdout, 1)
      end if

    ! All frame types below must be custom definitions and will draw on the
    ! assumption that the first num_custom_ep defined frames line up with the
    ! prm_ep_frames array
    else if (frames(n)%type .eq. 4) then
      ep_lcl_crd(1, 1, n) = prm_ep_frames(n)%d1
    else if (frames(n)%type .ge. 5 .and. frames(n)%type .le. 7) then
      ep_lcl_crd(1, 1, n) = prm_ep_frames(n)%d1
      ep_lcl_crd(2, 1, n) = prm_ep_frames(n)%d2
    else if (frames(n)%type .eq. 8) then
      ep_lcl_crd(1, 1, n) = prm_ep_frames(n)%d1
      ep_lcl_crd(2, 1, n) = prm_ep_frames(n)%d2
      ep_lcl_crd(3, 1, n) = prm_ep_frames(n)%d3
    else if (frames(n)%type .eq. 9) then
      ep_lcl_crd(1, 1, n) = prm_ep_frames(n)%d1
    else if (frames(n)%type .eq. 10) then
      ep_lcl_crd(1, 1, n) = prm_ep_frames(n)%d1
      ep_lcl_crd(2, 1, n) = prm_ep_frames(n)%d2
      ep_lcl_crd(3, 1, n) = prm_ep_frames(n)%d3
    else
      write(mdout, *) 'EXTRA_PTS: unexpected frame type value: ', frames(n)%type
      call mexit(mdout, 1)
    end if  ! ( frames(n)%type .eq. 1 .or. frames(n)%type .eq. 3 )

  end do  !  n = 1, frame_cnt
   
   if (verbose .gt. 3) then
     write(mdout, *) 'frames:'
     do n = 1, frame_cnt
       write(mdout, 666) n, frames(n)%parent_atm, isymbl(frames(n)%parent_atm), &
         frames(n)%ep_cnt, &
         frames(n)%extra_pnt(1), frames(n)%extra_pnt(2), frames(n)%type, &
         frames(n)%frame_atm1, frames(n)%frame_atm2, frames(n)%frame_atm3
       write(mdout, 667) ep_lcl_crd(1,1,n), ep_lcl_crd(2,1,n), ep_lcl_crd(3,1,n), &
         ep_lcl_crd(1,2,n), ep_lcl_crd(2,2,n), ep_lcl_crd(3,2,n)
     end do
   end if
   
  return
  
  ! Free allocated memory
  deallocate(custom_ep_tbl)

666 format(1x, 2i7, 1x, a4, 7i7)
667 format(1x, 6(1x, f10.4))

end subroutine define_frames 

!*******************************************************************************!
! Internal Subroutine:  fill_bonded
!
! Description: <TBS>
!
!*******************************************************************************

subroutine fill_bonded(max11, max12, max13, max14, &
                       num11, num12, num13, num14, &
                       list11, list12, list13, list14, &
                       enghbrs, numnghbr, epowner, verbose)

  use pmemd_lib_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer       :: max11
  integer       :: max12
  integer       :: max13
  integer       :: max14
  integer       :: num11
  integer       :: num12
  integer       :: num13
  integer       :: num14
  integer       :: list11(2, max11)
  integer       :: list12(2, max12)
  integer       :: list13(2, max13)
  integer       :: list14(2, max14)
  integer       :: enghbrs(5, *)
  integer       :: numnghbr(3, *)
  integer       :: epowner(*)
  integer       :: verbose

! Local variables:

  integer       :: n, ifail

  num11 = 0
  num12 = 0
  num13 = 0
  num14 = 0

  do n = 1, natom

    if (numnghbr(3, n) .eq. 1) then
      num11 = num11 + 1
      if (num11 .gt. max11) goto 100
      list11(1, num11) = n
      list11(2, num11) = enghbrs(1, n)
    else if (numnghbr(3, n) .eq. 2) then
      if (num11  + 3 .gt. max11) goto 100
      num11 = num11 + 1
      if (num11 .gt. max11) goto 100
      list11(1, num11) = n
      list11(2, num11) = enghbrs(1, n)
      num11 = num11 + 1
      if (num11 .gt. max11) goto 100
      list11(1, num11) = n
      list11(2, num11) = enghbrs(2, n)
      num11 = num11 + 1
      if (num11 .gt. max11) goto 100
      list11(1, num11) = enghbrs(1, n)
      list11(2, num11) = enghbrs(2, n)
    else if (numnghbr(3, n) .eq. 3) then
      write(mdout, *) 'fill_bonded: should not be here!'
      call mexit(mdout, 1)
    end if

  end do

  ! Bonds:

  call do_bond_pairs(list12, num12, max12, epowner, numnghbr, enghbrs, ifail)

  if (ifail .eq. 1) then
    write(mdout, *) 'fill_bonded: max12 exceeded!!'
    call mexit(mdout, 1)
  end if

  call sort_pairs(list12, num12, natom)

  ! Angles:

  call do_angle_pairs(list13, num13, max13, epowner, numnghbr, enghbrs, ifail)

  if (ifail .eq. 1) then
    write(mdout, *) 'fill_bonded: max13 exceeded!!'
    call mexit(mdout, 1)
  end if

  call sort_pairs(list13, num13, natom)

  ! Dihedrals:

  call do_dihed_pairs(list14, num14, max14, epowner, numnghbr, enghbrs, ifail)

  if (ifail .eq. 1) then
    write(mdout, *) 'fill_bonded: max14 exceeded!!'
    call mexit(mdout, 1)
  end if

  call sort_pairs(list14, num14, natom)
   
  if (verbose .gt. 0) &
    write(mdout, '(a, 4i6)') '| EXTRA PNTS fill_bonded: num11-14 = ', &
      num11, num12, num13, num14

  if (verbose .gt. 3) then
    write(mdout, *) '$$$$$$$$$$$$$$$$$ 1-1 pairs $$$$$$$$$$$$$$$$$$$$$$$$'
    do n = 1, num11
      write(mdout, 666) n, list11(1, n), list11(2, n)
    end do
    write(mdout, *) '$$$$$$$$$$$$$$$$$ 1-2 pairs $$$$$$$$$$$$$$$$$$$$$$$$'
    do n = 1, num12
      write(mdout, 666) n, list12(1, n), list12(2, n)
    end do
    write(mdout, *) '$$$$$$$$$$$$$$$$$ 1-3 pairs $$$$$$$$$$$$$$$$$$$$$$$$'
    do n = 1, num13
      write(mdout, 666) n, list13(1, n), list13(2, n)
    end do
    write(mdout, *) '$$$$$$$$$$$$$$$$$ 1-4 pairs $$$$$$$$$$$$$$$$$$$$$$$$'
    do n = 1, num14
      write(mdout, 666) n, list14(1, n), list14(2, n)
    end do
    write(mdout, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
  end if

  return

100 write(mdout, *) 'fill_bonded: max11 exceeded!!'
  call mexit(mdout, 1)

666 format(1x, i5, ':', 3x, 'i, j = ', i5, 1x, i5)

end subroutine fill_bonded 

!*******************************************************************************!
! Internal Subroutine:  redo_masked
!
! Description: <TBS>
!
!*******************************************************************************

subroutine redo_masked(natom, iblo, inb, nnb, &
                       num11, num12, num13, num14, &
                       list11, list12, list13, list14, offset, test)

  use pmemd_lib_mod
  use prmtop_dat_mod, only : next_mult_fac

  implicit none

! Formal arguments:

  integer       :: natom
  integer       :: iblo(*)
  integer       :: inb(*)
  integer       :: nnb
  integer       :: num11
  integer       :: num12
  integer       :: num13
  integer       :: num14
  integer       :: list11(2, *)
  integer       :: list12(2, *)
  integer       :: list13(2, *)
  integer       :: list14(2, *)
  integer       :: offset(*)
  integer       :: test(*)
   
! Local variables:

  integer       :: j, n, m, ntot

  ! Build the mask list from list11-14. Make sure no duplication:
   
  do n = 1, natom
    iblo(n) = 0
    offset(n) = 0
    test(n) = 0
  end do
   
  ! Pass 1:  Fill iblo.
   
  call add_one_list_iblo(iblo, list11, num11)
  call add_one_list_iblo(iblo, list12, num12)
  call add_one_list_iblo(iblo, list13, num13)
  call add_one_list_iblo(iblo, list14, num14)
   
  ! Check totals while finding offsets, resetting iblo:
   
  ntot = 0
  do n = 1, natom
      offset(n) = ntot
      ntot = ntot + iblo(n)
      iblo(n) = 0
  end do

  if (ntot .gt. next_mult_fac * nnb) then
    write(mdout, *) 'EXTRA POINTS: nnb too small! '
    write(mdout, *) 'nnb, ntot = ', nnb, ntot
    call mexit(mdout, 1)
  end if
   
  ! Pass 2 fill inb, redo iblo:
   
  call add_one_list_inb(iblo, inb, offset, list11, num11)
  call add_one_list_inb(iblo, inb, offset, list12, num12)
  call add_one_list_inb(iblo, inb, offset, list13, num13)
  call add_one_list_inb(iblo, inb, offset, list14, num14)
   
  ! Pass 3 filter inb, remove duplicate entries:
   
  do n = 1, natom - 1
    do m = 1, iblo(n)
      j = inb(offset(n) + m)
      if (test(j) .ne. n ) then
        test(j) = n
      else
        inb(offset(n) + m) = 0
      end if
    end do
  end do

  return

end subroutine redo_masked 

!*******************************************************************************!
! Internal Subroutine:  build_nb14
!
! Description: <TBS>
!
!*******************************************************************************

subroutine build_nb14(nb14, nb14_cnt, maxnb14, epowner, numnghbr, enghbrs, &
                      atm_cnt, chngmask, verbose)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer               :: nb14(3, *)
  integer               :: nb14_cnt
  integer               :: maxnb14
  integer               :: epowner(*)
  integer               :: numnghbr(3, *)
  integer               :: enghbrs(5, *)
  integer               :: atm_cnt
  integer               :: chngmask
  integer               :: verbose

! Local variables:

  integer               :: ifail, n

  nb14_cnt = 0

  call do_14pairs(nb14, nb14_cnt, maxnb14, epowner, numnghbr, enghbrs, &
                  ifail, chngmask)

  if (ifail .eq. 1) then
    write(mdout, *) 'exceeded maxnb14 in build14: check extra_pnts_nb14.fpp'
    call mexit(mdout, 1)
  end if

  call sort_pairs_14_nb(nb14, nb14_cnt, atm_cnt)

  if (verbose .gt. 0) write(mdout, '(a, i6)') &
    '| EXTRA_PTS, build_nb14: num of 14 terms = ', nb14_cnt

  if (verbose .gt. 3) then
    write(mdout, *) '$$$$$$$$$$$$$$$$$$$$$$$  1-4 nb list $$$$$$$$$$'
    do n = 1, nb14_cnt
      write(mdout, 666) n, nb14(1, n), nb14(2, n)
    end do
    write(mdout, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
  end if

  return

666 format(1x, i5, ':', 3x, 'i, j = ', i5, 1x, i5)

end subroutine build_nb14 

!*******************************************************************************!
! Internal Subroutine:  copy_nb14
!
! Description: <TBS>
!
!*******************************************************************************

subroutine copy_nb14(from14, to14, nb14_cnt)

  implicit none

! Formal arguments:

  integer       :: from14(3, *)
  integer       :: to14(3, *)
  integer       :: nb14_cnt

! Local variables:

  integer       :: n

  do n = 1, nb14_cnt
    to14(1, n) = from14(1, n)
    to14(2, n) = from14(2, n)
    to14(3, n) = from14(3, n)
  end do

  return

end subroutine copy_nb14 

!*******************************************************************************!
! Internal Subroutine:  trim_bonds
!
! Description: <TBS>
!
!*******************************************************************************

subroutine trim_bonds(bond_cnt, bonds, epowner)
  
  implicit none

! Formal arguments:

  integer, intent(in out)               :: bond_cnt
  type(bond_rec), intent(in out)        :: bonds(*)
  integer, intent(in)                   :: epowner(*)

! Local variables:

  integer       :: atm_i, atm_j
  integer       :: n, m
   
  write(mdout, '(a, 2i6)') &
    '|      EXTRA_PTS, trim_bonds: num bonds BEFORE trim =', bond_cnt, 0
   
  m = 0
  
  do n = 1, bond_cnt

    atm_i = bonds(n)%atm_i
    atm_j = bonds(n)%atm_j
      
    ! Only keep if neither is extra:
      
    if (epowner(atm_i) .eq. 0 .and. epowner(atm_j) .eq. 0) then
      m = m + 1
      bonds(m) = bonds(n)
    end if

  end do

  bond_cnt = m
   
  write(mdout, '(a, 2i6)') &
    '|      EXTRA_PTS, trim_bonds: num bonds AFTER  trim =', bond_cnt, 0

  return

end subroutine trim_bonds 

!*******************************************************************************!
! Internal Subroutine:  trim_theta
!
! Description: <TBS>
!
!*******************************************************************************

subroutine trim_theta(angle_cnt, angles, epowner)

  implicit none

! Formal arguments:

  integer, intent(in out)               :: angle_cnt
  type(angle_rec), intent(in out)       :: angles(*)
  integer, intent(in)                   :: epowner(*)

! Local variables:

  integer       :: atm_i, atm_k
  integer       :: n, m

  m = 0

  write(mdout, '(a, 2i6)') &
    '|      EXTRA_PTS, trim_theta: num angle BEFORE trim =', angle_cnt, 0
   
  do n = 1, angle_cnt

    atm_i = angles(n)%atm_i
    atm_k = angles(n)%atm_k
      
    ! Only keep if neither is extra:
      
    if (epowner(atm_i) .eq. 0 .and. epowner(atm_k) .eq. 0) then
      m = m + 1
      angles(m) = angles(n)
    end if

  end do

  angle_cnt = m
   
  write(mdout, '(a, 2i6)') &
    '|      EXTRA_PTS, trim_theta: num angle AFTER  trim =', angle_cnt, 0

  return

end subroutine trim_theta 

!*******************************************************************************!
! Internal Subroutine:  trim_phi
!
! Description: <TBS>
!
!*******************************************************************************

subroutine trim_phi(dihed_cnt, dihedrals, epowner)

  implicit none

! Formal arguments:

  integer, intent(in out)               :: dihed_cnt
  type(dihed_rec), intent(in out)       :: dihedrals(*)
  integer, intent(in)                   :: epowner(*)

! Local variables:

  integer       :: atm_i, atm_l
  integer       :: n, m

  m = 0

  write(mdout, '(a, 2i6)') &
    '|      EXTRA_PTS, trim_phi:  num diheds BEFORE trim =', dihed_cnt, 0
   
  do n = 1, dihed_cnt

    atm_i = dihedrals(n)%atm_i
    atm_l = iabs(dihedrals(n)%atm_l)
      
    ! Only keep if neither is extra:
      
    if (epowner(atm_i) .eq. 0 .and. epowner(atm_l) .eq. 0) then
      m = m + 1
      dihedrals(m) = dihedrals(n)
    end if

  end do

  dihed_cnt = m
   
  write(mdout, '(a, 2i6)') &
    '|      EXTRA_PTS, trim_phi:  num diheds AFTER  trim =', dihed_cnt, 0

  return

end subroutine trim_phi 

!*******************************************************************************!
! Internal Subroutine:  do_bond_pairs
!
! Description: <TBS>
!
!*******************************************************************************

subroutine do_bond_pairs(list, num, maxp, epowner, numnghbr, enghbrs, ifail)

  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer       :: list(2, *)
  integer       :: num
  integer       :: maxp
  integer       :: epowner(*)
  integer       :: numnghbr(3, *)
  integer       :: enghbrs(5, *)
  integer       :: ifail

! Local variables:

  integer       :: atm_i, atm_j
  integer       :: i, j, k, l, n

  ifail = 0

  do n = 1, nbonh + nbona
      
    atm_i = gbl_bond(n)%atm_i
    atm_j = gbl_bond(n)%atm_j
      
    ! Check neither is extra; Count this bond and also extra attached:
      
    if (epowner(atm_i) .eq. 0 .and. epowner(atm_j) .eq. 0) then
      k = numnghbr(3, atm_i)
      l = numnghbr(3, atm_j)
      if (num + 1 + k + l + k * l .gt. maxp) then
        ifail = 1
        return
      end if
      num = num + 1
      list(1, num) = atm_i
      list(2, num) = atm_j
      do i = 1, k
        num = num + 1
        list(1, num) = enghbrs(i, atm_i)
        list(2, num) = atm_j
      end do
      do j = 1, l
        num = num + 1
        list(1, num) = atm_i
        list(2, num) = enghbrs(j, atm_j)
      end do
      do i = 1, k
        do j = 1, l
          num = num + 1
          list(1, num) = enghbrs(i, atm_i)
          list(2, num) = enghbrs(j, atm_j)
        end do
      end do
    end if
  end do  !  n = 1, nbonh + nbona

  return

end subroutine do_bond_pairs 

!*******************************************************************************!
! Internal Subroutine:  do_angle_pairs
!
! Description: <TBS>
!
!*******************************************************************************

subroutine do_angle_pairs(list, num, maxp, epowner, numnghbr, enghbrs, ifail)

  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer       :: list(2, *)
  integer       :: num
  integer       :: maxp
  integer       :: epowner(*)
  integer       :: numnghbr(3, *)
  integer       :: enghbrs(5, *)
  integer       :: ifail

! Local variables:

  integer       :: atm_i, atm_k
  integer       :: i, j, k, l, n

  ifail = 0

  do n = 1, ntheth + ntheta
      
    atm_i = gbl_angle(n)%atm_i
    atm_k = gbl_angle(n)%atm_k
      
    ! Check neither is extra; Count this bond and also extra attached:
      
    if (epowner(atm_i) .eq. 0 .and. epowner(atm_k) .eq. 0) then
      k = numnghbr(3, atm_i)
      l = numnghbr(3, atm_k)
      if (num + 1 + k + l + k * l .gt. maxp) then
        ifail = 1
        return
      end if
      num = num + 1
      list(1, num) = atm_i
      list(2, num) = atm_k
      do i = 1, k
        num = num + 1
        list(1, num) = enghbrs(i, atm_i)
        list(2, num) = atm_k
      end do
      do j = 1, l
        num = num + 1
        list(1, num) = atm_i
        list(2, num) = enghbrs(j, atm_k)
      end do
      do i = 1, k
        do j = 1, l
          num = num + 1
          list(1, num) = enghbrs(i, atm_i)
          list(2, num) = enghbrs(j, atm_k)
        end do
      end do
    end if
  end do  !  n = 1, ntheth + ntheta

  return

end subroutine do_angle_pairs 

!*******************************************************************************!
! Internal Subroutine:  do_dihed_pairs
!
! Description: <TBS>
!
!*******************************************************************************

subroutine do_dihed_pairs(list, num, maxp, epowner, numnghbr, enghbrs, ifail)

  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer       :: list(2, *)
  integer       :: num
  integer       :: maxp
  integer       :: epowner(*)
  integer       :: numnghbr(3, *)
  integer       :: enghbrs(5, *)
  integer       :: ifail

! Local variables:

  integer       :: atm_i, atm_l
  integer       :: i, j, k, l, n

  ifail = 0

  do n = 1, nphih + nphia
      
    ! Sometimes second index negative (improper dihedrals):
      
    atm_i = gbl_dihed(n)%atm_i
    atm_l = iabs(gbl_dihed(n)%atm_l)
      
    ! Check neither is extra; Count this bond and also extra attached:
      
    if (epowner(atm_i) .eq. 0 .and. epowner(atm_l) .eq. 0) then
      k = numnghbr(3, atm_i)
      l = numnghbr(3, atm_l)
      if (num + 1 + k + l + k * l .gt. maxp) then
        ifail = 1
        return
      end if
      num = num + 1
      list(1, num) = atm_i
      list(2, num) = atm_l
      do i = 1, k
        num = num + 1
        list(1, num) = enghbrs(i, atm_i)
        list(2, num) = atm_l
      end do
      do j = 1, l
        num = num + 1
        list(1, num) = atm_i
        list(2, num) = enghbrs(j, atm_l)
      end do
      do i = 1, k
        do j = 1, l
          num = num + 1
          list(1, num) = enghbrs(i, atm_i)
          list(2, num) = enghbrs(j, atm_l)
        end do
      end do
    end if
  end do  !  n = 1, nphih + nphia

  return

end subroutine do_dihed_pairs 

!*******************************************************************************!
! Internal Subroutine:  sort_pairs
!
! Description: <TBS>
!
!*******************************************************************************

subroutine sort_pairs(list, num, atm_cnt)

  implicit none

! Formal arguments:

  integer       :: list(2, *)
  integer       :: num
  integer       :: atm_cnt

! Local variables:

  integer       :: i, j, k, m, n, ntot
  integer       :: scr1(atm_cnt), scr2(atm_cnt), scr3(maxa)

  do n = 1, num
    i = list(1, n)
    j = list(2, n)
    if (i .lt. j) then
      list(1, n) = i
      list(2, n) = j
    else
      list(1, n) = j
      list(2, n) = i
    end if
  end do
   
  ! Now get rid of duplicates:
   
  ! First pass:
   
  do n = 1, atm_cnt
    scr1(n) = 0
  end do

  do n = 1, num
    i = list(1, n)
    scr1(i) = scr1(i) + 1
  end do

  scr2(1) = 0
  do n = 2, atm_cnt
    scr2(n) = scr2(n - 1) + scr1(n - 1)
    scr1(n - 1) = 0
  end do
  scr1(atm_cnt) = 0
   
  ! Second pass:
   
  do n = 1, num
    i = list(1, n)
    j = list(2, n)
    scr1(i) = scr1(i) + 1
    scr3(scr1(i) + scr2(i)) = j
  end do

  scr2(1:atm_cnt) = 0
   
  ! Now trim them:
   
  ntot = 0
  k = 0
  do n = 1, atm_cnt
    do m = 1, scr1(n)
      j = scr3(ntot + m)
      if (scr2(j) .ne. n) then
        k = k + 1
        list(1, k) = n
        list(2, k) = j
        scr2(j) = n
      end if
    end do
    ntot = ntot + scr1(n)
  end do
  num = k
   
  return
   
end subroutine sort_pairs 

!*******************************************************************************!
! Internal Subroutine:  sort_pairs_14_nb
!
! Description: <TBS>
!
!*******************************************************************************

subroutine sort_pairs_14_nb(nb_14_list, num, atm_cnt)

  implicit none

! Formal arguments:

  integer       :: nb_14_list(3, *)
  integer       :: num
  integer       :: atm_cnt

! Local variables:

  integer       :: i, j, k, m, n, ntot, parm_idx
  integer       :: scr1(atm_cnt), scr2(atm_cnt), scr3(maxa), scr4(maxa)

 !Indexes of nb_14_list = 1,j,parm_index
 !First pass, make sure the first index contains the lowest number. Ignore
 !index 3 here since this won't change anything.

  do n = 1, num
    i = nb_14_list(1, n)
    j = nb_14_list(2, n)
    if (i .lt. j) then
      nb_14_list(1, n) = i
      nb_14_list(2, n) = j
    else
      nb_14_list(1, n) = j
      nb_14_list(2, n) = i
    end if
  end do
   
  ! Now get rid of duplicates:
   
  ! First pass:
  scr1(1:atm_cnt) = 0 
  do n = 1, num
    i = nb_14_list(1, n)
    scr1(i) = scr1(i) + 1
  end do

  scr2(1) = 0
  do n = 2, atm_cnt
    scr2(n) = scr2(n - 1) + scr1(n - 1)
    scr1(n - 1) = 0
  end do
  scr1(atm_cnt) = 0
   
  ! Second pass:
   
  do n = 1, num
    i = nb_14_list(1, n)
    j = nb_14_list(2, n)
    parm_idx = nb_14_list(3, n)
    scr1(i) = scr1(i) + 1
    scr3(scr1(i) + scr2(i)) = j
    scr4(scr1(i) + scr2(i)) = parm_idx
  end do

  scr2(1:atm_cnt) = 0
   
  ! Now trim them:
   
  ntot = 0
  k = 0
  do n = 1, atm_cnt
    do m = 1, scr1(n)
      j = scr3(ntot + m)
      parm_idx = scr4(ntot + m)
      if (scr2(j) .ne. n) then
        k = k + 1
        nb_14_list(1, k) = n
        nb_14_list(2, k) = j
        nb_14_list(3, k) = parm_idx
        scr2(j) = n
      end if
    end do
    ntot = ntot + scr1(n)
  end do
  num = k
   
  return

end subroutine sort_pairs_14_nb

!*******************************************************************************!
! Internal Subroutine:  add_one_list_iblo
!
! Description: <TBS>
!
!*******************************************************************************

subroutine add_one_list_iblo(iblo, list, num)

  implicit none

! Formal arguments:

  integer       :: iblo(*)
  integer       :: list(2, *)
  integer       :: num

! Local variables:

  integer       :: n, i

  do n = 1, num
    i = list(1, n)
    iblo(i) = iblo(i) + 1
  end do

  return

end subroutine add_one_list_iblo 

!*******************************************************************************!
! Internal Subroutine:  add_one_list_inb
!
! Description: <TBS>
!
!*******************************************************************************

subroutine add_one_list_inb(iblo, inb, offset, list, num)

  implicit none

! Formal arguments:

  integer       :: iblo(*)
  integer       :: inb(*)
  integer       :: offset(*)
  integer       :: list(2, *)
  integer       :: num

! Local variables:

  integer       :: n, i, j, m

  do n = 1, num
    i = list(1, n)
    j = list(2, n)
    m = offset(i)
    iblo(i) = iblo(i) + 1
    inb(m + iblo(i)) = j
  end do

  return

end subroutine add_one_list_inb 

!*******************************************************************************!
! Internal Subroutine:  do_14pairs
!
! Description: <TBS>
!
!*******************************************************************************

subroutine do_14pairs(nb_14_list, num, maxnb14, epowner, numnghbr, enghbrs, &
                      ifail, chngmask)

  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer               :: nb_14_list(3, *)
  integer               :: num
  integer               :: maxnb14
  integer               :: epowner(*)
  integer               :: numnghbr(3, *)
  integer               :: enghbrs(5, *)
  integer               :: ifail
  integer               :: chngmask

! Local variables:

  integer               :: atm_i, atm_l, parm_idx
  integer               :: i, j, n, ni, nl

  ifail = 0

  do n = 1, nphih + nphia

    if (gbl_dihed(n)%atm_k .gt. 0 .and. &
        gbl_dihed(n)%atm_l .gt. 0 .and. &
        gbl_fmn(gbl_dihed(n)%parm_idx) .gt. ZERO) then

      atm_i = gbl_dihed(n)%atm_i
      atm_l = gbl_dihed(n)%atm_l
      parm_idx = gbl_dihed(n)%parm_idx
      if (chngmask .eq. 0) then
        if (num + 1 .gt. maxnb14) then
          ifail = 1
          return
        end if
        num = num + 1
        nb_14_list(1, num) = atm_i
        nb_14_list(2, num) = atm_l
        nb_14_list(3, num) = parm_idx  !used for lookup into scnb and scee arrays.
      else
        ! Check neither is extra. Count this bond and also extra attached:
        if (epowner(atm_i) .eq. 0 .and. epowner(atm_l) .eq. 0) then
          ni = numnghbr(3, atm_i)
          nl = numnghbr(3, atm_l)
          if (num + 1 + ni + nl + ni * nl .gt. maxnb14) then
            ifail = 1
            return
          end if
          num = num + 1
          nb_14_list(1, num) = atm_i
          nb_14_list(2, num) = atm_l
          nb_14_list(3, num) = parm_idx  !used for lookup into scnb and scee arrays.
          do i = 1, ni
            num = num + 1
            nb_14_list(1, num) = enghbrs(i, atm_i)
            nb_14_list(2, num) = atm_l
            nb_14_list(3, num) = parm_idx  !used for lookup into scnb and scee arrays.
          end do
          do j = 1, nl
            num = num + 1
            nb_14_list(1, num) = atm_i
            nb_14_list(2, num) = enghbrs(j, atm_l)
            nb_14_list(3, num) = parm_idx  !used for lookup into scnb and scee arrays.
          end do
          do i = 1, ni
            do j = 1, nl
              num = num + 1
              nb_14_list(1, num) = enghbrs(i, atm_i)
              nb_14_list(2, num) = enghbrs(j, atm_l)
              nb_14_list(3, num) = parm_idx  !used for lookup into scnb and scee arrays.
            end do
          end do
        end if
      end if  ! ( chngmask .eq. 0 )
    end if
  end do  !  n = 1, nphih + nphia

  return

end subroutine do_14pairs

!*******************************************************************************!
! Subroutine:  all_local_to_global
!
! Description:  Put EP in position in world coord system based on the
!               position of the frame and the local coordinates, for all
!               coordinates; done only in master.
!
! Frames from Stone and Alderton Mol Phys. 56, 5, 1047 (1985) plus weird
! modification for carbonyl to symmetrize.
!*******************************************************************************

subroutine all_local_to_global(crd, frames, lcl_crd, frame_cnt)
  
  implicit none

! Formal arguments:

  double precision      :: crd(3, *)
  type(ep_frame_rec)    :: frames(*)
  double precision      :: lcl_crd(3, 2, *)
  integer               :: frame_cnt
   
! Local variables:

  GBFloat      :: uvec(3)
  GBFloat      :: vvec(3)
  GBFloat      :: ave(3)
  GBFloat      :: diff(3)
  GBFloat      :: usiz
  GBFloat      :: vsiz
  GBFloat      :: asiz
  GBFloat      :: dsiz
  double precision      :: f(3, 3)
  double precision      :: a(3), b(3), c(3), d(3), dv(3), rab(3), rac(3), rbc(3), rPerp(3)
  double precision      :: abac(3), rad(3), rja(3), rjb(3), rm(3)
  double precision      :: magdvec, invab, invab2, abbc0abab, magrP2, rabfac, rPfac
  integer               :: j, k, m, n, j1, j2, j3, j4
   
  if (frame_cnt .eq. 0) return
  
  do n = 1, frame_cnt

    if (frames(n)%type .le. 3) then

      ! Hard-wired virtual site / extra point frames first
      if (frames(n)%type .eq. 1) then
        do m = 1, 3
          a(m) = crd(m, frames(n)%frame_atm1)
          b(m) = crd(m, frames(n)%frame_atm2)
          c(m) = crd(m, frames(n)%frame_atm3)
        end do
      else if (frames(n)%type .eq. 2) then
        do m = 1, 3
          a(m) = HALF * (crd(m, frames(n)%frame_atm1) + &
                 crd(m, frames(n)%frame_atm2))
          b(m) = crd(m, frames(n)%parent_atm)
          c(m) = HALF * (crd(m, frames(n)%frame_atm3) + &
                 crd(m, frames(n)%frame_atm2))
        end do
      end if
    
      ! Z-axis along symmmetry axis of b midway between
      ! unit vector to a and unit vector to c; points opposite:

      usiz = 0.d0
      vsiz = 0.d0
      do m = 1, 3
        uvec(m) = a(m) - b(m)
        usiz = usiz + uvec(m) * uvec(m)
        vvec(m) = c(m) - b(m)
        vsiz = vsiz + vvec(m) * vvec(m)
      end do
      usiz = sqrt(usiz)
      vsiz = sqrt(vsiz)
      asiz = 0.d0
      dsiz = 0.d0
      do m = 1, 3
        uvec(m) = uvec(m) / usiz
        vvec(m) = vvec(m) / vsiz
        ave(m) = (uvec(m) + vvec(m)) / TWO
        asiz = asiz + ave(m) * ave(m)
        diff(m) = (vvec(m) - uvec(m)) / TWO
        dsiz = dsiz + diff(m) * diff(m)
      end do
      asiz = sqrt(asiz)
      dsiz = sqrt(dsiz)
      do m = 1, 3
        f(m, 3) = -ave(m) / asiz
        f(m, 1) = diff(m) / dsiz
      end do
      f(1, 2) = f(2, 3) * f(3, 1) - f(3, 3) * f(2, 1)
      f(2, 2) = f(3, 3) * f(1, 1) - f(1, 3) * f(3, 1)
      f(3, 2) = f(1, 3) * f(2, 1) - f(2, 3) * f(1, 1)
      do k = 1, frames(n)%ep_cnt
        j = frames(n)%extra_pnt(k)
        do m = 1, 3
          crd(m, j) = crd(m, frames(n)%parent_atm) + &
                      lcl_crd(1, k, n) * f(m, 1) + &
                      lcl_crd(2, k, n) * f(m, 2) + &
                      lcl_crd(3, k, n) * f(m, 3)
        end do
      end do
    else

      ! Custom virtual site / extra point frames below
      if (frames(n)%type .eq. 4) then
         
        ! Type IV: mdgx style 1, EP on the line between two parent atoms
        j = frames(n)%extra_pnt(1)
        do m = 1, 3
          a(m) = crd(m, frames(n)%parent_atm)
          b(m) = crd(m, frames(n)%frame_atm1)
          crd(m, j) = a(m) + lcl_crd(1, 1, n)*(b(m) - a(m))
        end do
      else if (frames(n)%type .eq. 5) then

        ! Type V: mdgx style 2, EP on the line between two parent atoms
        j = frames(n)%extra_pnt(1)
        do m = 1, 3
          a(m) = crd(m, frames(n)%parent_atm)
          b(m) = crd(m, frames(n)%frame_atm1)
          c(m) = crd(m, frames(n)%frame_atm2)
          crd(m, j) = a(m) + lcl_crd(1, 1, n)*(b(m) - a(m)) + &
                      lcl_crd(2, 1, n)*(c(m) - a(m))
        end do
      else if (frames(n)%type .eq. 6) then

        ! Type VI: mdgx style 3, EP on the bisector of the parent atom and
        ! two frame atoms, the bisector being determined by an arbitrary
        ! fraction of the distance between the two frame atoms.
        j = frames(n)%extra_pnt(1)
        magdvec = 0.0d0 
        do m = 1, 3
          a(m) = crd(m, frames(n)%parent_atm)
          b(m) = crd(m, frames(n)%frame_atm1)
          c(m) = crd(m, frames(n)%frame_atm2)
          dv(m) = b(m) - a(m) + lcl_crd(2, 1, n)*(c(m) - b(m))
          magdvec = magdvec + dv(m)*dv(m)
        end do
        magdvec = 1.0d0 / sqrt(magdvec)
        do m = 1, 3
          crd(m, j) = a(m) + lcl_crd(1, 1, n) * dv(m) * magdvec; 
        end do
      else if (frames(n)%type .eq. 7) then

        ! Type VII: mdgx style 4, EP hanging off of one atom, in the plane of
        !           its parent and two other atoms, at a fixed angle and distance
        !           from the parent atom.
        j = frames(n)%extra_pnt(1)
        j1 = frames(n)%parent_atm
        j2 = frames(n)%frame_atm1
        j3 = frames(n)%frame_atm2
        do m = 1, 3
          a(m) = crd(m, j1)
          b(m) = crd(m, j2)
          c(m) = crd(m, j3)
          rab(m) = b(m) - a(m)
          rac(m) = c(m) - a(m)
          rbc(m) = c(m) - b(m)
        end do
        invab2 = 1.0 / (rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))
        abbc0abab = (rab(1)*rbc(1) + rab(2)*rbc(2) + rab(3)*rbc(3)) * invab2
        rPerp(1) = rbc(1) - abbc0abab*rab(1)
        rPerp(2) = rbc(2) - abbc0abab*rab(2)
        rPerp(3) = rbc(3) - abbc0abab*rab(3)
        magrP2 = rPerp(1)*rPerp(1) + rPerp(2)*rPerp(2) + rPerp(3)*rPerp(3)
        rabfac = lcl_crd(1, 1, n) * cos(lcl_crd(2, 1, n)) * sqrt(invab2)
        rPfac  = lcl_crd(1, 1, n) * sin(lcl_crd(2, 1, n)) / sqrt(magrP2)
        do m = 1, 3
          crd(m, j) = crd(m, j1) + rabfac*rab(m) + rPfac*rPerp(m)
        end do
      else if (frames(n)%type .eq. 8) then

        ! Type VIII: mdgx style 5, EP out of plane anchored by three atoms in plane
        j = frames(n)%extra_pnt(1)
        j1 = frames(n)%parent_atm
        j2 = frames(n)%frame_atm1
        j3 = frames(n)%frame_atm2
        do m = 1, 3
          rab(m) = crd(m, j2) - crd(m, j1)
          rac(m) = crd(m, j3) - crd(m, j1)
        end do
        abac(1) = rab(2)*rac(3) - rab(3)*rac(2)
        abac(2) = rab(3)*rac(1) - rab(1)*rac(3)
        abac(3) = rab(1)*rac(2) - rab(2)*rac(1)
        do m = 1, 3
          crd(m, j) = crd(m, j1) + rab(m)*lcl_crd(1, 1, n) + rac(m)*lcl_crd(2, 1, n) + &
                      abac(m)*lcl_crd(3, 1, n)
        end do
      else if (frames(n)%type .eq. 9) then

        ! Type IX: EP again on the line between two parent atoms,
        !          but with a normalization for absolute distance
        j = frames(n)%extra_pnt(1)
        invab = 0.0d0 
        do m = 1, 3
          a(m) = crd(m, frames(n)%parent_atm)
          rab(m) = crd(m, frames(n)%frame_atm1) - a(m)
          invab = invab + rab(m)*rab(m)
        end do
        invab = 1.0d0 / sqrt(invab)
        do m = 1, 3
          crd(m, j) = a(m) + lcl_crd(1, 1, n)* invab * rab(m)
        end do
      else if (frames(n)%type .eq. 10) then

        ! Type X: mdgx style 10, EP at a fixed distance on a pyramid based on four atoms
        j = frames(n)%extra_pnt(1)
        j1 = frames(n)%parent_atm
        j2 = frames(n)%frame_atm1
        j3 = frames(n)%frame_atm2
        j4 = frames(n)%frame_atm3
        do m = 1, 3
          rab(m) = crd(m, j2) - crd(m, j1)
          rac(m) = crd(m, j3) - crd(m, j1)
          rad(m) = crd(m, j4) - crd(m, j1)
          rja(m) = lcl_crd(1, 1, n)*rac(m) - rab(m)
          rjb(m) = lcl_crd(2, 1, n)*rad(m) - rab(m)
        end do
        rm(1) = rja(2)*rjb(3) - rja(3)*rjb(2)
        rm(2) = rja(3)*rjb(1) - rja(1)*rjb(3)
        rm(3) = rja(1)*rjb(2) - rja(2)*rjb(1)
        magdvec =  lcl_crd(3, 1, n) / sqrt(rm(1)*rm(1) + rm(2)*rm(2) + rm(3)*rm(3))
        do m = 1, 3
          crd(m, j) = crd(m, j1) + magdvec*rm(m)
        end do
      end if
    end if
  end do  !  n = 1, frame_cnt

  return

end subroutine all_local_to_global 

!*******************************************************************************!
! Subroutine:  local_to_global
!
! Description:  Put EP in position in world coord system based on the
!               position of the frame and the local coordinates.
!
! Frames from Stone and Alderton Mol Phys. 56, 5, 1047 (1985) plus weird
! modification for carbonyl to symmetrize.
!*******************************************************************************

subroutine local_to_global(crd, frames, lcl_crd, frame_cnt)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  double precision      :: crd(3, *)
  type(ep_frame_rec)    :: frames(*)
  double precision      :: lcl_crd(3, 2, *)
  integer               :: frame_cnt
   
! Local variables:

  GBFloat      :: uvec(3)
  GBFloat      :: vvec(3)
  GBFloat      :: ave(3)
  GBFloat      :: diff(3)
  GBFloat      :: usiz
  GBFloat      :: vsiz
  GBFloat      :: dsiz
  GBFloat      :: asiz
  double precision      :: f(3, 3)
  double precision      :: a(3), b(3), c(3), d(3), dv(3), rab(3), rac(3), rbc(3), rPerp(3)
  double precision      :: abac(3), rad(3), rja(3), rjb(3), rm(3)
  integer               :: i, j, k, m, j1, j2, j3, j4
  integer               :: frame_id
#ifdef MPI
  integer               :: my_lst_idx
#endif
  double precision      :: magdvec, invab, invab2, abbcOabab, magrP2, rabfac, rPfac
  
  if (frame_cnt .eq. 0) return

#ifdef MPI
  do my_lst_idx = 1, my_ep_frame_cnt
    frame_id = gbl_my_ep_frame_lst(my_lst_idx)
#else
  do frame_id = 1, frame_cnt
#endif

    if (frames(frame_id)%type .le. 3) then

      ! Hard-wired virtual site / extra point frames first
      if (frames(frame_id)%type .eq. 1) then
        do m = 1, 3
          a(m) = crd(m, frames(frame_id)%frame_atm1)
          b(m) = crd(m, frames(frame_id)%frame_atm2)
          c(m) = crd(m, frames(frame_id)%frame_atm3)
        end do
      else if (frames(frame_id)%type .eq. 2) then
        do m = 1, 3
          a(m) = HALF * (crd(m, frames(frame_id)%frame_atm1) + &
                 crd(m, frames(frame_id)%frame_atm2))
          b(m) = crd(m, frames(frame_id)%parent_atm)
          c(m) = HALF * (crd(m, frames(frame_id)%frame_atm3) + &
                 crd(m, frames(frame_id)%frame_atm2))
        end do
      end if

      ! Z-axis along symmmetry axis of b midway between
      ! unit vector to a and unit vector to c; points opposite:

      usiz = 0.d0
      vsiz = 0.d0
      do m = 1, 3
        uvec(m) = a(m) - b(m)
        usiz = usiz + uvec(m) * uvec(m)
        vvec(m) = c(m) - b(m)
        vsiz = vsiz + vvec(m) * vvec(m)
      end do
      usiz = sqrt(usiz)
      vsiz = sqrt(vsiz)
      asiz = 0.d0
      dsiz = 0.d0
      do m = 1, 3
        uvec(m) = uvec(m) / usiz
        vvec(m) = vvec(m) / vsiz
        ave(m) = (uvec(m) + vvec(m)) / TWO
        asiz = asiz + ave(m) * ave(m)
        diff(m) = (vvec(m) - uvec(m)) / TWO
        dsiz = dsiz + diff(m) * diff(m)
      end do
      asiz = sqrt(asiz)
      dsiz = sqrt(dsiz)
      do m = 1, 3
        f(m, 3) = -ave(m) / asiz
        f(m, 1) = diff(m) / dsiz
      end do
      f(1, 2) = f(2, 3) * f(3, 1) - f(3, 3) * f(2, 1)
      f(2, 2) = f(3, 3) * f(1, 1) - f(1, 3) * f(3, 1)
      f(3, 2) = f(1, 3) * f(2, 1) - f(2, 3) * f(1, 1)
      do k = 1, frames(frame_id)%ep_cnt
        j = frames(frame_id)%extra_pnt(k)
        do m = 1, 3
          crd(m, j) = crd(m, frames(frame_id)%parent_atm) + &
                      lcl_crd(1, k, frame_id) * f(m, 1) + &
                      lcl_crd(2, k, frame_id) * f(m, 2) + &
                      lcl_crd(3, k, frame_id) * f(m, 3)
        end do
      end do
    else

      ! Custom virtual site / extra point frames below
      if (frames(frame_id)%type .eq. 4) then

        ! Type IV: mdgx style 1, EP on the line between two parent atoms
        j = frames(frame_id)%extra_pnt(1)
        do m = 1, 3
          f(m, 1) = crd(m, frames(frame_id)%parent_atm)
          f(m, 2) = crd(m, frames(frame_id)%frame_atm1)
          crd(m, j) = f(m, 1) + lcl_crd(1, 1, frame_id)*(f(m, 2) - f(m, 1))
        end do
      else if (frames(frame_id)%type .eq. 5) then

        ! Type V: mdgx style 2, EP in the plane with three frame atoms
        j = frames(frame_id)%extra_pnt(1)
        do m = 1, 3
          a(m) = crd(m, frames(frame_id)%parent_atm)
          b(m) = crd(m, frames(frame_id)%frame_atm1)
          c(m) = crd(m, frames(frame_id)%frame_atm2)
          crd(m, j) = a(m) + lcl_crd(1, 1, frame_id)*(b(m) - a(m)) + &
                      lcl_crd(2, 1, frame_id)*(c(m) - a(m))
        end do
      else if (frames(frame_id)%type .eq. 6) then

        ! Type VI: mdgx style 3, EP on the bisector of the parent atom and
        ! two frame atoms, the bisector being determined by an arbitrary
        ! fraction of the distance between the two frame atoms.
        j = frames(frame_id)%extra_pnt(1)
        magdvec = 0.0d0
        do m = 1, 3
          a(m) = crd(m, frames(frame_id)%parent_atm)
          b(m) = crd(m, frames(frame_id)%frame_atm1)
          c(m) = crd(m, frames(frame_id)%frame_atm2)
          dv(m) = b(m) - a(m) + lcl_crd(2, 1, frame_id)*(c(m) - b(m))
          magdvec = magdvec + dv(m)*dv(m)
        end do
        magdvec = 1.0d0 / sqrt(magdvec)
        do m = 1, 3
          crd(m, j) = a(m) + lcl_crd(1, 1, frame_id) * dv(m) * magdvec;
        end do
      else if (frames(frame_id)%type .eq. 7) then

        ! Type VII: mdgx style 4, EP hanging off of one atom, in the plane of
        !           its parent and two other atoms, at a fixed angle and distance
        !           from the parent atom.
        j = frames(frame_id)%extra_pnt(1)
        j1 = frames(frame_id)%parent_atm
        j2 = frames(frame_id)%frame_atm1
        j3 = frames(frame_id)%frame_atm2
        do m = 1, 3
          a(m) = crd(m, j1)
          b(m) = crd(m, j2)
          c(m) = crd(m, j3)
          rab(m) = b(m) - a(m)
          rbc(m) = c(m) - b(m)
        end do
        invab2 = 1.0 / (rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))
        abbcOabab = (rab(1)*rbc(1) + rab(2)*rbc(2) + rab(3)*rbc(3)) * invab2
        rPerp(1) = rbc(1) - abbcOabab*rab(1)
        rPerp(2) = rbc(2) - abbcOabab*rab(2)
        rPerp(3) = rbc(3) - abbcOabab*rab(3)
        magrP2 = rPerp(1)*rPerp(1) + rPerp(2)*rPerp(2) + rPerp(3)*rPerp(3)
        rabfac = lcl_crd(1, 1, frame_id) * cos(lcl_crd(2, 1, frame_id)) * sqrt(invab2)
        rPfac  = lcl_crd(1, 1, frame_id) * sin(lcl_crd(2, 1, frame_id)) / sqrt(magrP2)
        do m = 1, 3
          crd(m, j) = crd(m, j1) + rabfac*rab(m) + rPfac*rPerp(m)
        end do
      else if (frames(frame_id)%type .eq. 8) then

        ! Type VIII: mdgx style 5, EP out of plane anchored by three atoms in plane
        j = frames(frame_id)%extra_pnt(1)
        j1 = frames(frame_id)%parent_atm
        j2 = frames(frame_id)%frame_atm1
        j3 = frames(frame_id)%frame_atm2
        do m = 1, 3
          a(m) = crd(m, j1)
          b(m) = crd(m, j2)
          c(m) = crd(m, j3)
          rab(m) = b(m) - a(m)
          rac(m) = c(m) - a(m)
        end do
        abac(1) = rab(2)*rac(3) - rab(3)*rac(2)
        abac(2) = rab(3)*rac(1) - rab(1)*rac(3)
        abac(3) = rab(1)*rac(2) - rab(2)*rac(1)
        do m = 1, 3
           crd(m, j) = crd(m, j1) + rab(m)*lcl_crd(1, 1, frame_id) + &
                       rac(m)*lcl_crd(2, 1, frame_id) + abac(m)*lcl_crd(3, 1, frame_id)
        end do
      else if (frames(frame_id)%type .eq. 9) then

        ! Type IX: EP again on the line between two parent atoms,
        !          but with a normalization for absolute distance
        j = frames(frame_id)%extra_pnt(1)
        invab = 0.0d0 
        do m = 1, 3
          a(m) = crd(m, frames(frame_id)%parent_atm)
          rab(m) = crd(m, frames(frame_id)%frame_atm1) - a(m)
          invab = invab + rab(m)*rab(m)
        end do
        invab = 1.0d0 / sqrt(invab)
        do m = 1, 3
          crd(m, j) = a(m) + (lcl_crd(1, 1, frame_id) * invab * rab(m))
        end do
      else if (frames(frame_id)%type .eq. 10) then

        ! Type X: mdgx style 10, EP at a fixed distance on a pyramid based on four atoms
        j = frames(frame_id)%extra_pnt(1)
        j1 = frames(frame_id)%parent_atm
        j2 = frames(frame_id)%frame_atm1
        j3 = frames(frame_id)%frame_atm2
        j4 = frames(frame_id)%frame_atm3
        do m = 1, 3
          rab(m) = crd(m, j2) - crd(m, j1)
          rac(m) = crd(m, j3) - crd(m, j1)
          rad(m) = crd(m, j4) - crd(m, j1)
          rja(m) = lcl_crd(1, 1, frame_id)*rac(m) - rab(m)
          rjb(m) = lcl_crd(2, 1, frame_id)*rad(m) - rab(m)
        end do
        rm(1) = rja(2)*rjb(3) - rja(3)*rjb(2)
        rm(2) = rja(3)*rjb(1) - rja(1)*rjb(3)
        rm(3) = rja(1)*rjb(2) - rja(2)*rjb(1)
        magdvec = lcl_crd(3, 1, frame_id) / sqrt(rm(1)*rm(1) + rm(2)*rm(2) + rm(3)*rm(3))
        do m = 1, 3
          crd(m, j) = crd(m, j1) + magdvec*rm(m)
        end do
      end if
    end if
  end do  !  frame_id = 1, frame_cnt

  return

end subroutine local_to_global 

!----------------------------------------------------------------------------------------------
! orient_frc: transfer forces from EP to main atoms in frame.  Frames from Stone and Alderton,
!             Mol Phys. 56, 5, 1047 (1985).
!
! Arguments:
!   crd:        atomic coordinates for the whole simulation
!   frc:        atomic forces for the whole simulation
!   frames:     atoms and frame style types that determine the locations of extra points
!   lcl_crd:    local coordinate information, originally used for placement of the 'old' EP
!               frame types, here used for managing the new EP frame types
!   frame_cnt:  the number of frames
!----------------------------------------------------------------------------------------------
subroutine orient_frc(crd, frc, framevir, frames, lcl_crd, frame_cnt)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  double precision      :: crd(3, *)
  double precision      :: frc(3, *)
  double precision      :: framevir(3, 3)
  type(ep_frame_rec)    :: frames(*)
  double precision      :: lcl_crd(3, 2, *)
  integer               :: frame_cnt
   
! Local variables:

  GBFloat      :: u(3), v(3), w(3)
  GBFloat      :: up(3), vp(3), diff(3)
  GBFloat      :: usiz, vsiz, wsiz, upsiz, vpsiz
  GBFloat      :: dotdu, dotdv, dphidu, dphidv, dphidw
  GBFloat      :: c, s, uvdis, vudis, du(3), dv(3)
  GBFloat      :: force(3), torque(3), rel(3)
  integer               :: i, j, k, j1, j2, j3, j4, m, n
  integer               :: frame_id
#ifdef MPI
  integer               :: my_lst_idx
#endif
  double precision      :: ap(3), bp(3), ab(3), cp(3), raEP(3), pbold(3), rPerp(3)
  double precision      :: rab(3), rac(3), rbc(3), rad(3), rja(3), rjb(3), rjab(3)
  double precision      :: fb(3), fc(3), fd(3), rm(3), rt(3) 
  double precision      :: magdvec, magdvec2, gamma, pbfac, f1fac, f2fac, g01, g02, g12
  double precision      :: nf1fac, nf2fac, abbcOabab, F1, F2, F3, invab, invab2, magab, magab2
  double precision      :: magPr, magPr2, invPr, invPr2, fproj, cfx, cfy, cfz

  ! Motions of the frame can be described in terms of rotations about the
  ! unit vectors u and v from frame_atm2 to frame_atm1, frame_atm3 respecively
  ! and rotation about the unit cross product w.
   
  if (frame_cnt .eq. 0) return
   
  framevir = 0.d0

#ifdef MPI
  do my_lst_idx = 1, my_ep_frame_cnt
    frame_id = gbl_my_ep_frame_lst(my_lst_idx)
#else
  do frame_id = 1, frame_cnt
#endif

    if (frames(frame_id)%type .le. 3) then

      ! As with placement, force transmission for hardwired frames goes first.
      ! Start by collecting the force and torque due to ALL extra points, then
      ! immediately zero the forces on extra point / virtaul site particles. 
      force  = 0.d0
      torque = 0.d0
      i = frames(frame_id)%parent_atm
      do k = 1, frames(frame_id)%ep_cnt
        j = frames(frame_id)%extra_pnt(k)
        force(1:3) = force(1:3) + frc(1:3, j)
        rel(1:3) = crd(1:3, j) - crd(1:3, i)
         
        ! Get transferred force component of virial:         
        do m = 1, 3
          framevir(:, m) = framevir(:, m) + frc(:, j) * rel(m)
        end do
         
        ! Torque is rel x frc:
        torque(1) = torque(1) + rel(2) * frc(3, j) - rel(3) * frc(2, j)
        torque(2) = torque(2) + rel(3) * frc(1, j) - rel(1) * frc(3, j)
        torque(3) = torque(3) + rel(1) * frc(2, j) - rel(2) * frc(1, j)
        frc(:, j) = 0.d0
      end do

      ! Collect information about the frame geometry
      if (frames(frame_id)%type .eq. 1) then
        do m = 1, 3
          ap(m) = crd(m, frames(frame_id)%frame_atm1)
          bp(m) = crd(m, frames(frame_id)%frame_atm2)
          cp(m) = crd(m, frames(frame_id)%frame_atm3)
        end do
      else if (frames(frame_id)%type .eq. 2) then
        do m = 1, 3
          ap(m) = HALF * (crd(m, frames(frame_id)%frame_atm1) + &
                  crd(m, frames(frame_id)%frame_atm2))
          bp(m) = crd(m, frames(frame_id)%parent_atm)
          cp(m) = HALF * (crd(m, frames(frame_id)%frame_atm3) + &
                  crd(m, frames(frame_id)%frame_atm2))
        end do
      end if
      usiz = 0.d0
      vsiz = 0.d0
      do m = 1, 3
        u(m) = ap(m) - bp(m)
        usiz = usiz + u(m) * u(m)
        v(m) = cp(m) - bp(m)
        vsiz = vsiz + v(m) * v(m)
      end do
      usiz = sqrt(usiz)
      vsiz = sqrt(vsiz)
      w(1) = u(2) * v(3) - u(3) * v(2)
      w(2) = u(3) * v(1) - u(1) * v(3)
      w(3) = u(1) * v(2) - u(2) * v(1)
      wsiz = sqrt(w(1) * w(1) + w(2) * w(2) + w(3) * w(3))
      dotdu = 0.d0
      dotdv = 0.d0
      do m = 1, 3
        u(m) = u(m) / usiz
        v(m) = v(m) / vsiz
        w(m) = w(m) / wsiz
        diff(m) = v(m) - u(m)
        dotdu = dotdu + u(m) * diff(m)
        dotdv = dotdv + v(m) * diff(m)
      end do
      
      ! Get perps to u, v to get direction of motion of u or v
      ! due to rotation about the cross product vector w:
      upsiz = 0.d0
      vpsiz = 0.d0
      do m = 1, 3
        up(m) = diff(m) - dotdu * u(m)
        vp(m) = diff(m) - dotdv * v(m)
        upsiz = upsiz + up(m) * up(m)
        vpsiz = vpsiz + vp(m) * vp(m)
      end do
      upsiz = sqrt(upsiz)
      vpsiz = sqrt(vpsiz)
      do m = 1, 3
        up(m) = up(m) / upsiz
        vp(m) = vp(m) / vpsiz
      end do
      
      ! Negative of dot product of torque with unit vectors
      ! along u, v and w.  Give result of infinitesmal rotation
      ! along these vectors, i.e. dphi / dtheta = dot product.
      dphidu = -(torque(1) * u(1) + torque(2) * u(2) + torque(3) * u(3))
      dphidv = -(torque(1) * v(1) + torque(2) * v(2) + torque(3) * v(3))
      dphidw = -(torque(1) * w(1) + torque(2) * w(2) + torque(3) * w(3))
      
      ! Get projected distances between vectors:
      c = u(1) * v(1) + u(2) * v(2) + u(3) * v(3)
      s = sqrt(ONE - c * c)
      uvdis = usiz * s
      vudis = vsiz * s
      
      !---------------------------------------------------------------------
      ! Frame formed by bisector of u, v, its perp, and w.
      ! movement of u by dz out of plane -> rotation about v of -dz / uvdis
      ! since positive rotation about v move u in negative dir. wrt w
      ! dphi / dz = dphi / dtheta dtheta / dz = -dotvt / uvdis
      ! movement of v by dz out of plane -> rotation about u of dz / vudis
      ! movement of u by dy along up -> rotation about w of 1/2 dy / usiz
      ! since bisector only rotates 1/2 as much as u or v in isolation
      ! movement of v by dy along vperp -> rotation about w of 1/2 dy / vsiz
      ! movement of u by dx along u doesn't change frame
      ! movement of v by dx along v doesn't change frame
      ! So... du_du = 0, du_dw = -dotvt / uvdis, du_dup = dotwt / (2.d0 * usiz)
      ! So... dv_dv = 0, dv_dw = dotut / vudis, du_dup = dotwt / (2.d0 * usiz)
      !---------------------------------------------------------------------
      if (frames(frame_id)%type .eq. 1) then
        j1 = frames(frame_id)%frame_atm1
        j2 = frames(frame_id)%frame_atm2
        j3 = frames(frame_id)%frame_atm3
        do m = 1, 3
          du(m) = -w(m) * dphidv / uvdis + up(m) * dphidw / (TWO * usiz)
          dv(m) = w(m) * dphidu / vudis + vp(m) * dphidw / (TWO * vsiz)
          frc(m, j1) = frc(m, j1) - du(m)
          frc(m, j3) = frc(m, j3) - dv(m)
          frc(m, j2) = frc(m, j2) + dv(m) + du(m) + force(m)
        end do

        ! Get torque contribution to virial:
        do m = 1, 3
          framevir(:, m) = framevir(:, m) + du(:) * (ap(m) - bp(m)) + &
                           dv(:) * (cp(m) - bp(m))
        end do
      else if (frames(frame_id)%type .eq. 2) then
      
        ! Need to transfer forces from midpoints to atoms:
        j1 = frames(frame_id)%frame_atm1
        j2 = frames(frame_id)%frame_atm2
        j3 = frames(frame_id)%frame_atm3
        j4 = frames(frame_id)%parent_atm
        do m = 1, 3
          du(m) = -w(m) * dphidv / uvdis + up(m) * dphidw / (TWO * usiz)
          dv(m) = w(m) * dphidu / vudis + vp(m) * dphidw / (TWO * vsiz)
          frc(m, j1) = frc(m, j1) - HALF * du(m)
          frc(m, j3) = frc(m, j3) - HALF * dv(m)
          frc(m, j2) = frc(m, j2) - HALF * (du(m) + dv(m))
          frc(m, j4) = frc(m, j4) + dv(m) + du(m) + force(m)
        end do

        ! Get torque contribution to virial:
        do m = 1, 3
          framevir(:, m) = framevir(:, m) + du(:) * (ap(m) - bp(m)) + &
                           dv(:) * (cp(m) - bp(m))
        end do
      end if  ! ( frames(frame_id)%type .eq. 1 )
      
      !---------------------------------------------------------------------
      ! OTHER TYPE FRAME; NOT SEEN YET
      ! Frame formed by  u, its perp, and w
      ! movement of v in plane doesn't change frame
      ! movement of u by dz out of plane -> rotation about v of -dz / uvdis
      ! since positive rotation about v move u in negative dir. wrt w
      ! dphi / dz = dphi / dtheta dtheta / dz = -dotvt / uvdis
      ! movement of v by dz out of plane -> rotation about u of dz / vudis
      ! movement of u by dy along up -> rotation about w of dy / usiz
      ! since frame rotates as much as u in isolation
      !---------------------------------------------------------------------
      !  do m = 1, 3
      !    du(m) = -w(m) * dphidv / uvdis + up(m) * dphidw / usiz
      !    dv(m) = w(m) * dphidu / vudis
      !    frc(m, j1) = frc(m, j1) - du(m) + force(m)
      !    frc(m, j3) = frc(m, j3) - dv(m)
      !    frc(m, j2) = frc(m, j2) + dv(m) + du(m)
      !  enddo
      !
      !  get torque contribution to virial:
      !
      !  do m = 1, 3
      !    do l = 1, 3
      !      framevir(l, m) = framevir(l, m) + &
      !                       du(l) * (crd(m, j1) - crd(m, j2)) + &
      !                       dv(l) * (crd(m, j3) - crd(m, j2))
      !    enddo
      !  enddo
    else

      ! Force transmission for custom virtual sites
      if (frames(frame_id)%type .eq. 4) then

        ! Type IV: mdgx style 1, EP on the line between two parent atoms
        j  = frames(frame_id)%extra_pnt(1)
        j1 = frames(frame_id)%parent_atm
        j2 = frames(frame_id)%frame_atm1
        do m = 1, 3
          frc(m, j1) = frc(m, j1) + ((1.0d0 - lcl_crd(1,1,frame_id)) * frc(m, j))
          frc(m, j2) = frc(m, j2) + (lcl_crd(1,1,frame_id) * frc(m, j))
          frc(m, j)  = 0.0d0
        end do
      else if (frames(frame_id)%type .eq. 5) then

        ! Type V: mdgx style 2, EP on the plane determined by three parent parent atoms
        j  = frames(frame_id)%extra_pnt(1)
        j1 = frames(frame_id)%parent_atm
        j2 = frames(frame_id)%frame_atm1
        j3 = frames(frame_id)%frame_atm2
        do m = 1, 3
          frc(m, j1) = frc(m, j1) + &
                       (1.0d0 - lcl_crd(1,1,frame_id) - lcl_crd(2,1,frame_id)) * frc(m, j)
          frc(m, j2) = frc(m, j2) + lcl_crd(1,1,frame_id) * frc(m, j)
          frc(m, j3) = frc(m, j3) + lcl_crd(2,1,frame_id) * frc(m, j)
          frc(m, j)  = 0.0d0
        end do
     else if (frames(frame_id)%type .eq. 6) then

        ! Type VI: mdgx style 3, EP on the bisector of the parent atom and
        ! two frame atoms, the bisector being determined by an arbitrary
        ! fraction of the distance between the two frame atoms.
        j  = frames(frame_id)%extra_pnt(1)
        j1 = frames(frame_id)%parent_atm
        j2 = frames(frame_id)%frame_atm1
        j3 = frames(frame_id)%frame_atm2
        gamma = 0.0d0
        pbfac = 0.0d0
        do m = 1, 3
          ap(m)   = crd(m, j2) - crd(m, j1)
          ab(m)   = crd(m, j3) - crd(m, j2)
          raEP(m) = crd(m, j) - crd(m, j1)
          dv(m) = ap(m) + lcl_crd(2,1,frame_id) * ab(m)
          gamma = gamma + dv(m)*dv(m)
          pbfac = pbfac + raEP(m)*frc(m, j)
        end do
        gamma = lcl_crd(1,1,frame_id) / sqrt(gamma)
        pbfac = pbfac / (raEP(1)*raEP(1) + raEP(2)*raEP(2) + raEP(3)*raEP(3))
        do m = 1, 3
          pbold(m) = frc(m, j) - pbfac*raEP(m)
          frc(m, j1) = frc(m, j1) + frc(m, j) - gamma*pbold(m)
          frc(m, j2) = frc(m, j2) + ((1.0d0 - lcl_crd(2,1,frame_id)) * gamma * pbold(m))
          frc(m, j3) = frc(m, j3) + (         lcl_crd(2,1,frame_id)  * gamma * pbold(m))
        end do

        ! There were no contributions to the virial form the first two custom frame types,
        ! but this fixed-distance frame, normalizing the parent :: virtual site distance,
        ! requires calculation of a virial contribution.
        do m = 1, 3
          do n = 1, 3
            framevir(n, m) = framevir(n, m) + raEP(m)*frc(n, j) + (ap(m) * gamma * pbold(n))
          end do
        end do

        ! Zero the force on the virtual particle (it was spared
        ! earlier, so that the virial could be calculated)
        do m = 1, 3
          frc(m, j) = 0.0d0
        end do
      else if (frames(frame_id)%type .eq. 7) then
         
        ! Type VII: mdgx style 4, EP on an arm out beyond the parent atom, in the plane of
        !           three frame atoms in all, at a fixed distance and EP-parent-next atom
        !           angle.
        j  = frames(frame_id)%extra_pnt(1)
        j1 = frames(frame_id)%parent_atm
        j2 = frames(frame_id)%frame_atm1
        j3 = frames(frame_id)%frame_atm2
        do m = 1, 3
          rab(m)  = crd(m, j2) - crd(m, j1)
          rbc(m)  = crd(m, j3) - crd(m, j2)
          raEP(m) = crd(m, j)  - crd(m, j1)
        end do
        magab2  = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
        magdvec = (rab(1)*rbc(1) + rab(2)*rbc(2) + rab(3)*rbc(3)) / magab2
        magab   = sqrt(magab2)
        invab   = 1.0 / magab
        invab2  = invab * invab
        do m = 1, 3
          rPerp(m) = rbc(m) - magdvec*rab(m)
        end do
        magPr2 = rPerp(1)*rPerp(1) + rPerp(2)*rPerp(2) + rPerp(3)*rPerp(3)
        magPr  = sqrt(magPr2)
        invPr  = 1.0 / magPr
        invPr2 = invPr * invPr
        f1fac  = (   rab(1)*frc(1, j) +    rab(2)*frc(2, j) +    rab(3)*frc(3, j)) * invab2
        f2fac  = (rPerp(1)*frc(1, j) + rPerp(2)*frc(2, j) + rPerp(3)*frc(3, j)) * invPr2
        nf1fac = lcl_crd(1,1,frame_id) * cos(lcl_crd(2,1,frame_id)) * invab
        nf2fac = lcl_crd(1,1,frame_id) * sin(lcl_crd(2,1,frame_id)) * invPr
        abbcOabab = (rab(1)*rbc(1) + rab(2)*rbc(2) + rab(3)*rbc(3)) * invab2
        do m = 1, 3
          F1 = frc(m, j) - f1fac*rab(m)
          F2 = F1 - f2fac*rPerp(m)
          F3 = f1fac * rPerp(m)
          frc(m, j1) = frc(m, j1) + frc(m, j) - nf1fac*F1 + nf2fac*(abbcOabab*F2 + F3)
          frc(m, j2) = frc(m, j2) + nf1fac*F1 - nf2fac*(F2 + abbcOabab*F2 + F3)
          frc(m, j3) = frc(m, j3) + nf2fac * F2

          ! Compute the virial contributions usign the quantities already obtained
          rac(m) = crd(m, j3) - crd(m, j1)
          do n = 1, 3
            F1 = frc(n, j) - f1fac*rab(n)
            F2 = F1 - f2fac*rPerp(n)
            F3 = f1fac * rPerp(n)
            framevir(n, m) = framevir(n, m) + (raEP(m) * frc(n, j)) - &
                             (rab(m) * (nf1fac*F1 + nf2fac*(abbcOabab*F2 + F3))) - &
                             (rac(m) * nf2fac * F2)
          end do
        end do

        ! Zero the force on the virtual particle (it was spared
        ! earlier, so that the virial could be calculated)
        do m = 1, 3
          frc(m, j)  = 0.0
        end do
      else if (frames(frame_id)%type .eq. 8) then

        ! Type VIII: mdgx style 5, EP out of plane anchored by three atoms in plane
        j = frames(frame_id)%extra_pnt(1)
        j1 = frames(frame_id)%parent_atm
        j2 = frames(frame_id)%frame_atm1
        j3 = frames(frame_id)%frame_atm2
        do m = 1, 3
          rab(m) = crd(m, j2) - crd(m, j1)
          rac(m) = crd(m, j3) - crd(m, j1)
        end do
        g01 = lcl_crd(3,1,frame_id) * rac(3)
        g02 = lcl_crd(3,1,frame_id) * rac(2)
        g12 = lcl_crd(3,1,frame_id) * rac(1)
        fb(1) =  lcl_crd(1,1,frame_id)*frc(1, j) - g01*frc(2, j) + g02*frc(3, j)
        fb(2) =  g01*frc(1, j) + lcl_crd(1,1,frame_id)*frc(2, j) - g12*frc(3, j)
        fb(3) = -g02*frc(1, j) + g12*frc(2, j) + lcl_crd(1,1,frame_id)*frc(3, j)
        g01 = lcl_crd(3,1,frame_id) * rab(3)
        g02 = lcl_crd(3,1,frame_id) * rab(2)
        g12 = lcl_crd(3,1,frame_id) * rab(1)
        fc(1) =  lcl_crd(2,1,frame_id)*frc(1, j) + g01*frc(2, j) - g02*frc(3, j)
        fc(2) = -g01*frc(1, j) + lcl_crd(2,1,frame_id)*frc(2, j) + g12*frc(3, j)
        fc(3) =  g02*frc(1, j) - g12*frc(2, j) + lcl_crd(2,1,frame_id)*frc(3, j)
        do m = 1, 3
          frc(m, j1) = frc(m, j1) + frc(m, j) - fb(m) - fc(m)
          frc(m, j2) = frc(m, j2) + fb(m)
          frc(m, j3) = frc(m, j3) + fc(m)
        end do

        ! Compute the virial contributions
        do m = 1, 3
          raEP(m) = crd(m, j) - crd(m, j1)
          do n = 1, 3
            framevir(n, m) = framevir(n, m) + (raEP(m) * frc(n, j)) - &
                             (rab(m) * fb(n)) - (rac(m) * fc(n))
          end do
        end do

        ! Remove force on the virtual site
        do m = 1, 3
          frc(m,  j) = 0.0d0
        end do
      else if (frames(frame_id)%type .eq. 9) then

        ! Type IX: virtual site along the line between two points
        j = frames(frame_id)%extra_pnt(1)
        j1 = frames(frame_id)%parent_atm
        j2 = frames(frame_id)%frame_atm1
        invab = 0.0d0
        fproj = 0.0d0
        do m = 1, 3
          rab(m) = crd(m, j2) - crd(m, j1)
          invab = invab + rab(m)*rab(m)
          raEP(m) = crd(m, j) - crd(m, j1)
          fproj = fproj + raEP(m)*frc(m, j)
        end do
        invab = 1.0d0 / sqrt(invab)
        fproj = fproj / (raEP(1)*raEP(1) + raEP(2)*raEP(2) + raEP(3)*raEP(3))
        do m = 1, 3
          fb(m) = invab * lcl_crd(1,1,frame_id) * (frc(m, j) - fproj*raEP(m))
          frc(m, j1) = frc(m, j1) + frc(m, j) - fb(m)
          frc(m, j2) = frc(m, j2) + fb(m)
        end do

        ! Compute the virial contributions
        do m = 1, 3
          do n = 1, 3
            framevir(n, m) = framevir(n, m) + (raEP(m) * frc(n, j)) - (rab(m) * fb(n))
          end do
        end do

        ! Remove force on the virtual site
        do m = 1, 3
          frc(m, j) = 0.0d0
        end do
      else if (frames(frame_id)%type .eq. 10) then

        ! Type X: virtual site atop a pyramid determined by four frame atoms.
        ! This is the lengthiest, most complicated chain rule to execute.
        j = frames(frame_id)%extra_pnt(1)
        j1 = frames(frame_id)%parent_atm
        j2 = frames(frame_id)%frame_atm1
        j3 = frames(frame_id)%frame_atm2
        j4 = frames(frame_id)%frame_atm3
        do m = 1, 3
          rab(m)  = crd(m, j2) - crd(m, j1)
          rac(m)  = crd(m, j3) - crd(m, j1)
          rad(m)  = crd(m, j4) - crd(m, j1)
          rja(m)  = lcl_crd(1,1,frame_id) * rac(m) - rab(m)
          rjb(m)  = lcl_crd(2,1,frame_id) * rad(m) - rab(m)
          rjab(m) = rjb(m) - rja(m)
        end do
        rm(1) = rja(2)*rjb(3) - rja(3)*rjb(2)
        rm(2) = rja(3)*rjb(1) - rja(1)*rjb(3)
        rm(3) = rja(1)*rjb(2) - rja(2)*rjb(1)
        magdvec = 1.0d0 / sqrt(rm(1)*rm(1) + rm(2)*rm(2) + rm(3)*rm(3))
        magdvec2 = magdvec * magdvec
        cfx = lcl_crd(3,1,frame_id) * magdvec * frc(1, j)
        cfy = lcl_crd(3,1,frame_id) * magdvec * frc(2, j)
        cfz = lcl_crd(3,1,frame_id) * magdvec * frc(3, j)
        rt(1) = (rm(2)*rjab(3) - rm(3)*rjab(2)) * magdvec2
        rt(2) = (rm(3)*rjab(1) - rm(1)*rjab(3)) * magdvec2
        rt(3) = (rm(1)*rjab(2) - rm(2)*rjab(1)) * magdvec2
        fb(1) = -rm(1)*rt(1)*cfx + (rjab(3) - rm(2)*rt(1))*cfy - (rjab(2) + rm(3)*rt(1))*cfz
        fb(2) = -(rjab(3) + rm(1)*rt(2))*cfx - rm(2)*rt(2)*cfy + (rjab(1) - rm(3)*rt(2))*cfz
        fb(3) = (rjab(2) - rm(1)*rt(3))*cfx - (rjab(1) + rm(2)*rt(3))*cfy - rm(3)*rt(3)*cfz
        rt(1) = (rjb(2)*rm(3) - rjb(3)*rm(2)) * magdvec2 * lcl_crd(1,1,frame_id)
        rt(2) = (rjb(3)*rm(1) - rjb(1)*rm(3)) * magdvec2 * lcl_crd(1,1,frame_id)
        rt(3) = (rjb(1)*rm(2) - rjb(2)*rm(1)) * magdvec2 * lcl_crd(1,1,frame_id)
        fc(1) = -rm(1)*rt(1)*cfx - (lcl_crd(1,1,frame_id)*rjb(3) + rm(2)*rt(1))*cfy + &
                (lcl_crd(1,1,frame_id)*rjb(2) - rm(3)*rt(1))*cfz
        fc(2) = (lcl_crd(1,1,frame_id)*rjb(3) - rm(1)*rt(2))*cfx - rm(2)*rt(2)*cfy - &
                (lcl_crd(1,1,frame_id)*rjb(1) + rm(3)*rt(2))*cfz
        fc(3) = -(lcl_crd(1,1,frame_id)*rjb(2) + rm(1)*rt(3))*cfx + &
                (lcl_crd(1,1,frame_id)*rjb(1) - rm(2)*rt(3))*cfy - rm(3)*rt(3)*cfz
        rt(1) = (rm(2)*rja(3) - rm(3)*rja(2)) * magdvec2 * lcl_crd(2,1,frame_id)
        rt(2) = (rm(3)*rja(1) - rm(1)*rja(3)) * magdvec2 * lcl_crd(2,1,frame_id)
        rt(3) = (rm(1)*rja(2) - rm(2)*rja(1)) * magdvec2 * lcl_crd(2,1,frame_id)
        fd(1) = -rm(1)*rt(1)*cfx + (lcl_crd(2,1,frame_id)*rja(3) - rm(2)*rt(1))*cfy - &
                (lcl_crd(2,1,frame_id)*rja(2) + rm(3)*rt(1))*cfz
        fd(2) = -(lcl_crd(2,1,frame_id)*rja(3) + rm(1)*rt(2))*cfx - rm(2)*rt(2)*cfy + &
                (lcl_crd(2,1,frame_id)*rja(1) - rm(3)*rt(2))*cfz
        fd(3) = (lcl_crd(2,1,frame_id)*rja(2) - rm(1)*rt(3))*cfx + &
                -(lcl_crd(2,1,frame_id)*rja(1) + rm(2)*rt(3))*cfy - rm(3)*rt(3)*cfz
        do m = 1, 3
          frc(m, j1) = frc(m, j1) + frc(m, j) - fb(m) - fc(m) - fd(m)
          frc(m, j2) = frc(m, j2) + fb(m)
          frc(m, j3) = frc(m, j3) + fc(m)
          frc(m, j4) = frc(m, j4) + fd(m)
        end do

        ! Compute virial contributions
        do m = 1, 3
          raEP(m) = crd(m, j) - crd(m, j1)
          do n = 1, 3
            framevir(n, m) = framevir(n, m) + (raEP(m) * frc(n, j)) - (rab(m) * fb(n)) - &
                             (rac(m) * fc(n)) - (rad(m) * fd(n))
          end do
        end do

        ! Remove force on the virtual site        
        do m = 1, 3
          frc(m, j) = 0.0d0
        end do
      end if
    end if
  end do  !  frame_id = 1, frame_cnt

  return

end subroutine orient_frc 

!*******************************************************************************!
! Subroutine:   zero_extra_pnts_vec
!
! Description:  Set extra points vector entries to 0.0.
!
!*******************************************************************************

subroutine zero_extra_pnts_vec(vec, frames, frame_cnt)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  double precision      :: vec(3, *)
  type(ep_frame_rec)    :: frames(*)
  integer               :: frame_cnt
   
! Local variables:

  integer               :: i
  integer               :: frame_id
#ifdef MPI
  integer               :: my_lst_idx
#endif

  if (frame_cnt .eq. 0) return
   
#ifdef MPI
  do my_lst_idx = 1, my_ep_frame_cnt
    frame_id = gbl_my_ep_frame_lst(my_lst_idx)
#else
  do frame_id = 1, frame_cnt
#endif
    do i = 1, frames(frame_id)%ep_cnt
      vec(:, frames(frame_id)%extra_pnt(i)) = 0.d0
    end do
  end do
  return

end subroutine zero_extra_pnts_vec 

!*******************************************************************************!
! Internal Subroutine:  get_nb14_energy
!
! Description: <TBS>
!
!*******************************************************************************
subroutine get_nb14_energy(charge, crd, frc, iac, ico, cn1, cn2, &
                           nb14, nb14_cnt, ee14, enb14, e14vir)
#include "extra_pnts_nb14.i"
end subroutine get_nb14_energy 

!*******************************************************************************!
! Internal Subroutine:  get_nb14_energy_midpoint
!
! Description: <TBS>
!
!*******************************************************************************
#ifdef MPI
subroutine get_nb14_energy_midpoint(charge, crd, frc, iac, ico, cn1, cn2, &
                           nb14, nb14_cnt, ee14, enb14, e14vir)
#include "extra_pnts_nb14_midpoint.i"
end subroutine get_nb14_energy_midpoint 
#endif


#ifdef _OPENMP_
subroutine get_nb14_energy_gb(atm_cnt,charge, crd, frc, iac, ico, cn1, cn2, &
                           nb14, nb14_cnt, ee14, enb14, e14vir)
#define GBorn
#include "extra_pnts_nb14.i"
#undef GBorn
end subroutine get_nb14_energy_gb
#endif /*_OPENMP_*/

!*******************************************************************************!
! Subroutine:  fix_masses
!
! Description: We just fix up the atom masses; everything else is derived later
!              from them...
!
!*******************************************************************************

subroutine fix_masses(atm_cnt, amass, epowner)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: amass(*)
  integer               :: epowner(*)

! Local variables:

  integer               :: n
   
  ! Zero out mass for extra points;
   
  do n = 1, atm_cnt
    if (epowner(n) .ne. 0) then
      amass(n) = 0.d0
    end if
  end do
   
  return

end subroutine fix_masses 

end module extra_pnts_nb14_mod
