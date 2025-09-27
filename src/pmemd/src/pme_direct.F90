#include "copyright.i"
#include "mpi_err_check.i"

#include "include_precision.i"
!*******************************************************************************
!
! Module:  pme_direct_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module pme_direct_mod

  use gbl_datatypes_mod
  implicit none

  ! The following storage is per-process common; ie., it SHOULD be
  ! broadcast from the master to the other processes!

  integer, parameter    :: pme_direct_int_cnt = 1

  integer                  mxeedtab

  common / pme_direct_int / mxeedtab

  save  :: / pme_direct_int /

  double precision, allocatable, save   :: gbl_eed_cub(:)

  ! originally -1, I did +1 to accomodate 1 ghost bkt each side
  integer, allocatable  :: proc_crd_idx_tbl(:,:,:)
  integer, allocatable  :: tail_idx_tbl(:,:,:)
  type(atm_lst_rec), allocatable     :: proc_atm_lst(:) ! *2 to accomodate local and ghost atoms
  !these are related to pme exchange code in pme_direct
  logical, save                         :: alloc_flag = .false.
  !integer, dimension(2, 0:numtasks-1) :: b_send_map, a_send_map, d_send_map
!  integer, allocatable,save            :: b_send_map(:,:), a_send_map(:,:), d_send_map(:,:)
contains

!*******************************************************************************
!
! Subroutine:  init_pme_direct_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine init_pme_direct_dat(num_ints, num_reals)

  use gbl_constants_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use prmtop_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  double precision      :: dxdr
  integer               :: i

! Setup dxdr map from r to x in table lookup of eed

  dxdr = ew_coeff
  
! For eed table, assume all nonbond distances are less than 1.5 * es_cutoff
! between nonbond updates; i.e. no excess motion; this is enforced by
! es_cutoff.

  mxeedtab = int(dxdr * eedtbdns * es_cutoff * 1.5d0)

  call alloc_pme_direct_mem(num_ints, num_reals)

  return

end subroutine init_pme_direct_dat

!*******************************************************************************
!
! Subroutine:  alloc_pme_direct_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_pme_direct_mem(num_ints, num_reals)

  use pmemd_lib_mod
  use processor_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed

  allocate(gbl_eed_cub(4 * mxeedtab), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_reals = num_reals + size(gbl_eed_cub)

  gbl_eed_cub(:) = 0.d0

  alloc_flag = .true.

  return

end subroutine alloc_pme_direct_mem

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  bcast_pme_direct_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine bcast_pme_direct_dat

  use parallel_dat_mod

  implicit none

! Local variables:

  integer       :: num_ints, num_reals  ! returned values discarded

  call mpi_bcast(mxeedtab, pme_direct_int_cnt, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  if (.not. master) then
    num_ints = 0
    num_reals = 0
    call alloc_pme_direct_mem(num_ints, num_reals)
  end if

  ! The allocated data is not initialized from the master node.

  return

end subroutine bcast_pme_direct_dat

!*******************************************************************************
!
! Subroutine:  pme_list / pme_list_midpoint
!
! Description:  Handles set-up and error checking for calling of get_nb_list
!               which creates the nonbond list.
!
!*******************************************************************************
subroutine pme_list_midpoint(atm_cnt, crd, atm_maskdata, atm_mask &
                , proc_atm_maskdata, proc_atm_mask,nstep)

  use mdin_ctrl_dat_mod
  use cit_mod
  use constraints_mod
  use pmemd_lib_mod
  use pbc_mod
  use nb_pairlist_mod
  use img_mod
  use pme_recip_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use parallel_mod
  use pbc_mod
  use prmtop_dat_mod
  use timers_mod
  use nbips_mod
  use ti_mod
  use nb_exclusions_mod
  use processor_mod
  use  parallel_processor_mod
  use extra_pnts_nb14_mod
  use bonds_midpoint_mod
  use angles_midpoint_mod
  use dihedrals_midpoint_mod
  use shake_mod
  use mol_list_mod
  use dynamics_mod
  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: crd(3, atm_cnt)
  type(listdata_rec)    :: atm_maskdata(*)
  integer               :: atm_mask(*)
  type(listdata_rec)    :: proc_atm_maskdata(*)
  integer               :: proc_atm_mask(*)

! Local variables:
  integer               :: ifail
  integer               :: alloc_failed
  logical               :: dont_skip_belly_pairs
  integer               :: i,j,k
  ! We 0-base the following array for efficiency, but don't use atm_lst(0) so
  ! we can use 0 as the distinguished non-entry (NULL).

  double precision      :: fraction(3, atm_cnt)
!  type(cit_tbl_rec)     :: crd_idx_tbl(0 : cit_tbl_x_dim - 1, &
!                                       0 : cit_tbl_y_dim - 1, &
!                                       0 : cit_tbl_z_dim - 1)
  integer                   :: lc, gc, ac, lmax, lmin, lsum, gmax, gmin, gsum,amax, amin, asum
  integer                   :: bkt, atm_i_idx, atm_i, x,y,z
!debug
  integer                   :: atmid
  integer, save             :: counts=0
  integer                    :: nstep
  ! Comm debugging arrays
  double precision, allocatable :: dbg_atm_vel(:,:)
  Integer                       :: ghost_cnt
  double precision, allocatable :: dbg_atm_mass(:)
  double precision, allocatable :: dbg_atm_qterm(:)
  Integer, allocatable :: dbg_atm_iac(:)


  ! Under amber 7, when both members of an atom pair are in the belly (ie.,
  ! fixed), we don't add them to the nonbonded pairlist.  This saves time but
  ! makes the energies look different (in a manner that actually does not
  ! affect the simulation).
!  print*, "pme_list is called ", nstep, mytaskid

#ifdef USE_VTUNE
  if (mytaskid == 0) then
    call ITT_RESUME()
  end if
#endif
  counts = counts + 1
  if(allocated(proc_crd_idx_tbl))deallocate(proc_crd_idx_tbl)
  allocate(proc_crd_idx_tbl(0 : int_bkt_dimx + 1, &
                            0 : int_bkt_dimy + 1, &
                            0 : int_bkt_dimz + 1))
  if(allocated(tail_idx_tbl))deallocate(tail_idx_tbl)
  allocate(tail_idx_tbl(0 : int_bkt_dimx + 1, &
                        0 : int_bkt_dimy + 1, &
                        0 : int_bkt_dimz + 1))
  if(.not. allocated(bkt_atm_lst)) allocate(bkt_atm_lst(256, proc_numbkts))
  if(.not. allocated(bkt_atm_cnt)) allocate(bkt_atm_cnt(proc_numbkts))
 
!  write(15+mytaskid,*)"New List"

  dont_skip_belly_pairs = ibelly .eq. 0 .and. ti_mode .eq. 0
   
  gbl_saved_box(:) = pbc_box(:)     ! Needed when pressure scaling.

!we use counts to make sure we call it from 2nd nb list build
  if(counts .gt. 1) then
    call comm_3dimensions_3dbls(proc_atm_vel)

  end if
  call update_time(fcve_dist_time)
!!The following two calls probably should be called from 2nd call of the nb list
!if(counts .gt. 1) then
  !alternative
  do i=1, atm_cnt
    j=proc_atm_space(i)
  end do

!   old_proc_atm_to_full_list = proc_atm_to_full_list
  if(reorient_flag) then
   call resetup_processor(pbc_box, bkt_size)
   call proc_reorient_volume(counts, nstep)
!  do i = proc_num_atms+1, proc_num_atms + proc_ghost_num_atms
!         atmid = proc_atm_to_full_list(i)
!         proc_atm_space(atmid) = i
!  end do
!  do i = 1, proc_num_atms
!         atmid = proc_atm_to_full_list(i)
!         proc_atm_space(atmid) = i
!  end do
  end if
  if(.not. reorient_flag) then
   call proc_add_delete_atms(counts,nstep)
  end if
  reorient_flag = .false.
  call update_time(nonbond_time)

  call comm_setup
  if(allocated(proc_atm_lst)) deallocate(proc_atm_lst)
  allocate(proc_atm_lst(2*proc_atm_alloc_size)) 

  ! Communicate the ghost region's parameters
  call comm_3dimensions_1int(proc_atm_to_full_list)
  call comm_3dimensions_1int(proc_iac)
  call comm_3dimensions_1dbl(proc_atm_mass)
  call comm_3dimensions_1dbl(proc_atm_qterm)
  call comm_3dimensions_3ints(proc_atm_wrap, .true.)

  call update_time(fcve_dist_time)

  call proc_get_fract_crds(proc_num_atms, proc_fraction,int_bkt_dimx, int_bkt_dimy, int_bkt_dimz)
  call proc_setup_crd_idx_tbl(proc_num_atms, proc_fraction,proc_crd_idx_tbl,tail_idx_tbl,& 
                               proc_atm_lst, int_bkt_dimx, int_bkt_dimy, int_bkt_dimz)
  call proc_get_fract_crds_ghost(proc_num_atms, proc_ghost_num_atms, &
                                 proc_fraction, int_bkt_dimx, int_bkt_dimy, int_bkt_dimz)
  call proc_setup_crd_idx_tbl_ghost(proc_num_atms, proc_ghost_num_atms, proc_fraction, &
                proc_crd_idx_tbl, tail_idx_tbl, proc_atm_lst, int_bkt_dimx, int_bkt_dimy, int_bkt_dimz)

  call update_pme_time(cit_setup_timer)

  ! NOTE - It is important to note that no image mapping should occur prior
  !        to running get_nb_list().  We now exit setup_cit() without the
  !        used image map set also.  This is all being done to facilitate
  !        speed in image mapping in get_nb_list(), which should always be
  !        called immediately after setup_cit().


  ! convert the linked list data structure to SoA
  do bkt = 1,proc_numbkts
    atm_i_idx = 0
    call int_unflatten(bkt,x,y,z)
    atm_i = proc_crd_idx_tbl(x,y,z)
    do while(atm_i .ne. 0)
      atm_i_idx = atm_i_idx + 1
      bkt_atm_lst(atm_i_idx,bkt) = atm_i
      atm_i = proc_atm_lst(atm_i)%nxt ! next atom_i
    enddo
    bkt_atm_cnt(bkt) = atm_i_idx
!    if(atm_i_idx .eq. 0) then
!      write(*,"(A,I4,I3,I3,I3)") "#### empty bucket: ",mytaskid, proc_bkt_minx+x,proc_bkt_miny+y,proc_bkt_minz+z
!    end if
  enddo


   call atom_sorting()
!   do i=1,proc_num_atms
!     proc_sort_map(i) = i
!   end do

      ! but later we need to optimize so that we can call just one time
!  call build_atm_space() !now calling everytime nb list builds,&
if(counts .gt. 1) then !need to call for the 1st time
 call bonded_sorting()   
end if
if(counts .eq. 1) then !need to call for the 1st time
   call bonds_midpoint_setup()
   call angles_midpoint_setup()
   call dihedrals_midpoint_setup()
!  call mol_mass_calc(gbl_mol_cnt, proc_atm_mass)
endif


  do
    ifail = 0
    if (ifail .eq. 0) exit

    ! Deallocate the old array, grow it by 10%, reallocate,
    ! and go back up to the top to try again.

    deallocate(gbl_ipairs)

    ipairs_maxsize = 11 * ipairs_maxsize / 10       

    allocate(gbl_ipairs(ipairs_maxsize), stat = alloc_failed)
    if (alloc_failed .ne. 0) then
      call alloc_error('get_nb_list', 'ipairs array reallocation failed!');
    end if

! code that will handle ifail_dev .eq. 0

#ifdef PAIRLST_DBG
    write(mdout, '(/,a,5x,a,i12)') &
          '|', 'Nonbonded Pairs Reallocation:', ipairs_maxsize
#endif

  end do

  do
    ifail = 0

  call update_pme_time(build_list_timer)

! if(mytaskid .eq. 0) write(0,*) "PME list build list Time step", nstep
  call get_nb_list_midpoint( proc_num_atms+proc_ghost_num_atms, proc_atm_crd, &
                               typ_ico, ntypes, ibelly, &
                     es_cutoff, vdw_cutoff, skinnb, verbose, ifail, proc_corner_neighborbuckets, atm_maskdata, atm_mask)

!  call update_pme_time(build_nb_list_timer)
!need excluding list (excl_img_flags), currently they use atm_maskdata, and
!atm_mask
!need ico or typ_ico
!how do we handle ntypes and ibelly -- shall we reuse


    if (ifail .eq. 0) exit

    ! Deallocate the old array, grow it by 10%, reallocate,
    ! and go back up to the top to try again.
    
    deallocate(proc_ipairs)

    proc_ipairs_maxsize = 12 * proc_ipairs_maxsize / 10

    allocate(proc_ipairs(proc_ipairs_maxsize), stat = alloc_failed)

    if (alloc_failed .ne. 0) then
      call alloc_error('get_nb_list_midpoint', 'ipairs array reallocation failed!')
    end if

  end do

  call update_time(nonbond_time)

#if 1

#ifdef BEXDBG
  if (new_list .and. nstep .gt. 0) then
    call pme_exchange_midpoint(proc_atm_crd, (cit_nbonh+cit_nbona), (cit_ntheth+cit_ntheta), (cit_nphih+cit_nphia), neighbor_mpi_cnt, nstep)
  end if
#else
  if(counts .gt. 1) then ! we only call the routine after the 1st pme_list call
    call pme_exchange_midpoint(proc_atm_crd, (cit_nbonh+cit_nbona), (cit_ntheth+cit_ntheta), &
                               (cit_nphih+cit_nphia), neighbor_mpi_cnt,counts,nstep)
  end if  

  call update_time(fcve_dist_time)
!building the mult_vector for bonds
  if (ntf .le. 1) then
    if (cit_nbonh .gt. 0) then
        call  build_mult_vector_bond(cit_nbonh, cit_h_bond, proc_atm_crd,mult_vech_bond)
    end if
  end if

  if (ntf .le. 2) then
    if (cit_nbona .gt. 0) then
        call  build_mult_vector_bond(cit_nbona, cit_a_bond, proc_atm_crd,mult_veca_bond)
    end if
  end if
 
  call update_time(bond_time)

!building the mult_vector for angles
  if (ntf .le. 3) then
    if (cit_ntheth .gt. 0) then
       call  build_mult_vector_angle(cit_ntheth, cit_h_angle, proc_atm_crd,mult_vech_angle)
    end if
  end if

  if (ntf .le. 4) then
    if (cit_ntheta .gt. 0) then
        call  build_mult_vector_angle(cit_ntheta, cit_a_angle, proc_atm_crd,mult_veca_angle)
    end if
  end if

  call update_time(angle_time)

!building the mult_vector for dihedrals
  if (ntf .le. 5) then
    if (cit_nphih .gt. 0) then
        call  build_mult_vector_dihed(cit_nphih, cit_h_dihed, proc_atm_crd,mult_vech_dihed)
    end if
  end if

  if (ntf .le. 6) then
    if (cit_nphia .gt. 0) then
       call  build_mult_vector_dihed(cit_nphia, cit_a_dihed, proc_atm_crd,mult_veca_dihed)
    end if
  end if

  call update_time(dihedral_time)

  call make_nb_adjust_pairlst_midpoint(gbl_nb_adjust_pairlst, counts)


  call build_cit_nb14_midpoint(cit_nphia, cit_a_dihed,cit_nphih, cit_h_dihed)

if (ntc .ne. 1) then
  proc_shake_space(:) = 0
  do i = proc_num_atms + proc_ghost_num_atms,1,-1
         j = proc_atm_to_full_list(i)
         proc_shake_space(j) = i
  end do
 call shake_claim_midpoint()
end if
#endif

#endif

  log_used_img_cnt = dble(gbl_used_img_cnt)

  call update_pme_time(build_list_timer)
  call update_time(nonbond_time)

  call update_time(fcve_dist_time)

  call proc_save_all_atom_crds(proc_atm_crd, proc_saved_atm_crd)

  call zero_pme_time()

#ifdef USE_VTUNE
  if (mytaskid == 0) then
    call ITT_PAUSE()
  end if
#endif

! print sub-domain size stats
#if 0
!if(nstep > 5000) then
  if(mytaskid == 0) write(0,*) "Time step:", nstep
  lc=proc_num_atms
  gc=proc_ghost_num_atms
  ac=proc_num_atms+proc_ghost_num_atms
  write(0,"('['I3'] local: 'I5 ' ghost: 'I5' all: 'I5' --- Buckets:('I2','I2','I2')')") mytaskid, lc, gc, ac, (proc_bkt_maxx-proc_bkt_minx+1), (proc_bkt_maxy-proc_bkt_miny+1), (proc_bkt_maxz-proc_bkt_minz+1)
!  write(0,*) mytaskid, lc, gc, ac
  call mpi_allreduce(ac, amax, 1,mpi_int, mpi_max, pmemd_comm, err_code_mpi)
  call mpi_allreduce(ac, amin, 1,mpi_int, mpi_min, pmemd_comm, err_code_mpi)
  call mpi_allreduce(ac, asum, 1,mpi_int, mpi_sum, pmemd_comm, err_code_mpi)
  if(ac == amax) write(0,*) mytaskid, "max", amax
  if(ac == amin) write(0,*) mytaskid, "min", amin
  if(mytaskid == 0) write(0,*) mytaskid, "*** avg ***", asum/numtasks
!endif
#endif

  return

end subroutine pme_list_midpoint

#endif /*MPI*/


subroutine pme_list(atm_cnt, crd, atm_maskdata, atm_mask)

  use mdin_ctrl_dat_mod
  use cit_mod
  use constraints_mod
  use pmemd_lib_mod
  use pbc_mod
  use nb_pairlist_mod
  use img_mod
  use pme_blk_recip_mod
  use pme_slab_recip_mod
  use pme_recip_dat_mod
  use loadbal_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use parallel_mod
  use pbc_mod
  use prmtop_dat_mod
  use timers_mod
  use nbips_mod
  use ti_mod
  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: crd(3, atm_cnt)
  type(listdata_rec)    :: atm_maskdata(*)
  integer               :: atm_mask(*)

! Local variables:

  integer               :: ifail
  integer               :: alloc_failed
  logical               :: dont_skip_belly_pairs
  integer               :: i
  ! We 0-base the following array for efficiency, but don't use atm_lst(0) so
  ! we can use 0 as the distinguished non-entry (NULL).

  double precision      :: fraction(3, atm_cnt)
  type(cit_tbl_rec)     :: crd_idx_tbl(0 : cit_tbl_x_dim - 1, &
                                       0 : cit_tbl_y_dim - 1, &
                                       0 : cit_tbl_z_dim - 1)

  ! Under amber 7, when both members of an atom pair are in the belly (ie.,
  ! fixed), we don't add them to the nonbonded pairlist.  This saves time but
  ! makes the energies look different (in a manner that actually does not
  ! affect the simulation).

  dont_skip_belly_pairs = ibelly .eq. 0 .and. ti_mode .eq. 0
   
  gbl_saved_box(:) = pbc_box(:)     ! Needed when pressure scaling.
#ifndef MPI
  call save_all_atom_crds(atm_cnt, crd, gbl_atm_saved_crd)
#endif

  if(igb .ne. 0 .or. use_pme .ne. 0 .or. ips .ne. 0) then
    call get_fract_crds(atm_cnt, crd, fraction)
 
    call setup_cit(atm_cnt, fraction, crd_idx_tbl, &
                   gbl_atm_img_map, gbl_img_atm_map)
  end if

#ifdef MPI
  call init_used_img_map(gbl_used_img_map, gbl_used_img_cnt, gbl_used_img_lst)
  if( ips > 0 ) then
    call init_used_img_ips_map(gbl_used_img_ips_map, gbl_used_img_ips_cnt, gbl_used_img_ips_lst)
  endif
#endif /* MPI */

  call update_pme_time(cit_setup_timer)

  ! NOTE - It is important to note that no image mapping should occur prior
  !        to running get_nb_list().  We now exit setup_cit() without the
  !        used image map set also.  This is all being done to facilitate
  !        speed in image mapping in get_nb_list(), which should always be
  !        called immediately after setup_cit().

#ifdef MPI
  call map_pairlist_imgs(atm_cnt, crd_idx_tbl, &
                         fraction, atm_qterm, atm_iac, gbl_img_crd, &
                         gbl_img_qterm, gbl_img_iac, gbl_mapped_img_cnt, &
                         gbl_mapped_img_lst, gbl_img_atm_map)
#else
  if(igb .ne. 0 .or. use_pme .ne. 0 .or. ips .ne. 0) then
    call map_pairlist_imgs(atm_cnt, fraction, atm_qterm, atm_iac, gbl_img_crd, &
                          gbl_img_qterm, gbl_img_iac, gbl_img_atm_map)
  else
    gbl_img_crd(:,:)=crd(:,:)
    do i=1,atm_cnt
        gbl_img_qterm(i)=atm_qterm(i)
        gbl_img_atm_map(i) = i
        gbl_atm_img_map(i) = i
        gbl_img_iac(i) = atm_iac(i)
    end do
  end if
#endif

  do

    ifail = 0
  if( ips > 0 ) then
    call get_nb_ips_list(atm_cnt, crd_idx_tbl, gbl_img_crd, gbl_img_atm_map, &
#ifdef MPI
                     gbl_used_img_map, gbl_used_img_cnt, gbl_used_img_lst, &
                     gbl_used_img_ips_map, gbl_used_img_ips_cnt, gbl_used_img_ips_lst, &
#endif
                     typ_ico, fraction, gbl_tranvec, atm_maskdata, &
                     atm_mask, gbl_atm_img_map, gbl_excl_img_flags, &
                     gbl_img_iac, atm_igroup, ntypes, ibelly, &
                     es_cutoff, vdw_cutoff, skinnb, verbose, ifail)
  else

    if(use_pme .ne. 0) then
      call get_nb_list(atm_cnt, crd_idx_tbl, gbl_img_crd, gbl_img_atm_map, &

#ifdef MPI
                     gbl_used_img_map, gbl_used_img_cnt, gbl_used_img_lst, &
#endif
                     typ_ico, fraction, gbl_tranvec, atm_maskdata, &
                     atm_mask, gbl_atm_img_map, gbl_excl_img_flags, &
                     gbl_img_iac, atm_igroup, ntypes, ibelly, &
                     es_cutoff, vdw_cutoff, skinnb, verbose, ifail)
    else
      call get_nb_list_gp(atm_cnt, crd_idx_tbl, gbl_img_crd, gbl_img_atm_map, &
#ifdef MPI
                     gbl_used_img_map, gbl_used_img_cnt, gbl_used_img_lst, &
#endif
                     typ_ico, fraction, gbl_tranvec, atm_maskdata, &
                     atm_mask, gbl_atm_img_map, gbl_excl_img_flags, &
                     gbl_img_iac, atm_igroup, ntypes, ibelly, &
                     es_cutoff, vdw_cutoff, skinnb, verbose, ifail)
    end if
  endif
    if (ifail .eq. 0) exit

    ! Deallocate the old array, grow it by 10%, reallocate,
    ! and go back up to the top to try again.

    deallocate(gbl_ipairs)

    ipairs_maxsize = 11 * ipairs_maxsize / 10       

    allocate(gbl_ipairs(ipairs_maxsize), stat = alloc_failed)
    if (alloc_failed .ne. 0) then
      call alloc_error('get_nb_list', 'ipairs array reallocation failed!');
    end if

! code that will handle ifail_dev .eq. 0

#ifdef PAIRLST_DBG
    write(mdout, '(/,a,5x,a,i12)') &
          '|', 'Nonbonded Pairs Reallocation:', ipairs_maxsize
#endif

  end do

#ifdef MPI
  if (block_fft .eq. 0) then
    if (is_orthog .ne. 0) then
      call claim_slab_recip_imgs(atm_cnt, fraction, pbc_box, crd_idx_tbl, &
                                 gbl_img_crd, gbl_img_qterm, gbl_img_iac, &
                                 gbl_mapped_img_cnt, gbl_mapped_img_lst, &
                                 gbl_img_atm_map, gbl_used_img_map, &
                                 gbl_used_img_cnt, gbl_used_img_lst)
    else
      call claim_slab_recip_imgs_nonorthog(atm_cnt, fraction, crd_idx_tbl, &
                                           gbl_img_crd, gbl_img_qterm, &
                                           gbl_img_iac, &
                                           gbl_mapped_img_cnt, &
                                           gbl_mapped_img_lst, &
                                           gbl_img_atm_map, &
                                           gbl_used_img_map, &
                                           gbl_used_img_cnt, gbl_used_img_lst)
    end if
  else
    if (is_orthog .ne. 0) then
      call claim_blk_recip_imgs(atm_cnt, fraction, pbc_box, crd_idx_tbl, &
                                gbl_img_crd, gbl_img_qterm, gbl_img_iac, &
                                gbl_mapped_img_cnt, gbl_mapped_img_lst, &
                                gbl_img_atm_map, gbl_used_img_map, &
                                gbl_used_img_cnt, gbl_used_img_lst)
    else
      call claim_blk_recip_imgs_nonorthog(atm_cnt, fraction, crd_idx_tbl, &
                                          gbl_img_crd, gbl_img_qterm, &
                                          gbl_img_iac, gbl_mapped_img_cnt, &    
                                          gbl_mapped_img_lst, &
                                          gbl_img_atm_map, &
                                          gbl_used_img_map, &
                                          gbl_used_img_cnt, gbl_used_img_lst)
    end if
  end if

  log_used_img_cnt = dble(gbl_used_img_cnt)

  call update_pme_time(build_list_timer)
  call update_time(nonbond_time)

  if( ips > 0 ) then
    call get_send_ips_atm_lst(atm_cnt)
  endif
  call get_send_atm_lst(atm_cnt)

  call update_time(fcve_dist_time)
  call zero_pme_time()

  if (imin .eq. 0) then
    call save_used_atom_crds(atm_cnt, crd, gbl_atm_saved_crd)
  else
    ! Minimizations require all coords for skin check...
    call save_all_atom_crds(atm_cnt, crd, gbl_atm_saved_crd)
  end if
#else
  call update_pme_time(build_list_timer)
  call update_time(nonbond_time)
#endif /* MPI */


  return

end subroutine pme_list

#ifdef MPI

!*******************************************************************************
!
! Subroutine:  bonded_sorting
!
! Description: <TBS>
!*******************************************************************************
subroutine bonded_sorting()

  use pbc_mod
  use processor_mod 
  use bonds_midpoint_mod
  use angles_midpoint_mod
  use dihedrals_midpoint_mod
  implicit none
  integer        :: i,j,k , atmid
 
#if 1 
! proc_atm_space(:) = 0 !very expensive, need to only reset the previous nb steps
! proc_atm_space_ghosts(:) = 0!very expensive, need to only reset the previous nb steps
  do i = proc_num_atms+1, proc_num_atms + proc_ghost_num_atms
         atmid = proc_atm_to_full_list(i)
         proc_atm_space(atmid) = i
  end do
  do i = 1, proc_num_atms
         atmid = proc_atm_to_full_list(i)
         proc_atm_space(atmid) = i
  end do
  
  do j = 1, cit_nbonh
    cit_h_bond(j)%atm_i = proc_atm_space(old_cit_h_bond(j)%atm_i)
    cit_h_bond(j)%atm_j = proc_atm_space(old_cit_h_bond(j)%atm_j)
!  if(cit_h_bond(j)%atm_i==0 .or. cit_h_bond(j)%atm_j==0) print*, "bonded_sorting H", mytaskid, j, cit_h_bond(j)%atm_i, cit_h_bond(j)%atm_j, &
!   old_cit_h_bond(j)%atm_i, old_cit_h_bond(j)%atm_j, &
!   proc_atm_space(old_cit_h_bond(j)%atm_i), proc_atm_space(old_cit_h_bond(j)%atm_j)
  end do
  do j = 1, cit_nbona
    cit_a_bond(j)%atm_i = proc_atm_space(old_cit_a_bond(j)%atm_i)
    cit_a_bond(j)%atm_j = proc_atm_space(old_cit_a_bond(j)%atm_j)
  end do
  !save the angle index as global id
  do j = 1, cit_ntheth
    cit_h_angle(j)%atm_i = proc_atm_space(old_cit_h_angle(j)%atm_i)
    cit_h_angle(j)%atm_j = proc_atm_space(old_cit_h_angle(j)%atm_j)
    cit_h_angle(j)%atm_k = proc_atm_space(old_cit_h_angle(j)%atm_k)
  end do
  do j = 1, cit_ntheta
    cit_a_angle(j)%atm_i = proc_atm_space(old_cit_a_angle(j)%atm_i)
    cit_a_angle(j)%atm_j = proc_atm_space(old_cit_a_angle(j)%atm_j)
    cit_a_angle(j)%atm_k = proc_atm_space(old_cit_a_angle(j)%atm_k)
  end do
  !save the dihedrals index as global id
  do j = 1, cit_nphih
    cit_h_dihed(j)%atm_i = proc_atm_space(old_cit_h_dihed(j)%atm_i)
    cit_h_dihed(j)%atm_j = proc_atm_space(old_cit_h_dihed(j)%atm_j)
    cit_h_dihed(j)%atm_k =proc_atm_space(iabs(old_cit_h_dihed(j)%atm_k))*sign(1,old_cit_h_dihed(j)%atm_k)
    cit_h_dihed(j)%atm_l =proc_atm_space(iabs(old_cit_h_dihed(j)%atm_l))*sign(1,old_cit_h_dihed(j)%atm_l)
  end do
  do j = 1, cit_nphia
    cit_a_dihed(j)%atm_i = proc_atm_space(old_cit_a_dihed(j)%atm_i)
    cit_a_dihed(j)%atm_j = proc_atm_space(old_cit_a_dihed(j)%atm_j)
    cit_a_dihed(j)%atm_k = proc_atm_space(iabs(old_cit_a_dihed(j)%atm_k))*sign(1,old_cit_a_dihed(j)%atm_k)
    cit_a_dihed(j)%atm_l = proc_atm_space(iabs(old_cit_a_dihed(j)%atm_l))*sign(1,old_cit_a_dihed(j)%atm_l)
  end do

#endif

end subroutine bonded_sorting

!*******************************************************************************
!
! Subroutine:  proc_add_delete_atms
! This routine is called in the time step when "new_list" is true. after the "crd_distribution" 
! and before "atm_distribution"
!
! Description: when a new list is to be created, we need to check if any atoms
! moves in or out my local space. in this step, crd already transferred
! (crd_distribution) . Based on this crd (local and ghost crd) we find out which 
! atoms moved, and then update the proc_atm_crd() array. after that atm_distribution 
! will happen. Currently we are doing full ghost atm list transfer. Later we will 
! optimize so that only moved atoms are sent 
!*******************************************************************************

subroutine proc_add_delete_atms(counts,nstep)

  use pbc_mod
  use processor_mod !, only: proc_num_atms, proc_atm_crd, &
                     !   proc_min_x_crd, proc_max_x_crd, &
                     !   proc_min_y_crd, proc_max_y_crd, &
                     !   proc_min_z_crd, proc_max_z_crd, &
                     !   proc_ghost_num_atms, proc_atm_vel, &
                     !   proc_atm_last_vel, &
                     !   proc_atm_mass, proc_atm_qterm, &
                     !   proc_atm_to_full_list,&
                     !   proc_iac,myblockid, &  
                     !   sys_bkt_dims, sys_bkt_dimy, sys_bkt_dimz
  use bonds_midpoint_mod
  use angles_midpoint_mod
  use dihedrals_midpoint_mod
  use parallel_processor_mod
  implicit none

! Formal arguments:
  integer             :: counts

! Local variables:
!optimization note: please move the arrays to one time dynamic allocation in the
!setup routines
  integer               :: i,n
  integer               :: add_index(proc_num_atms_min_bound)
  integer               :: delete_index(proc_num_atms_min_bound)
  integer               :: deleted, added 
  integer               :: index
  double precision      :: tmp1, tmp2, tmp3 
  double precision      :: tmp4, tmp5, tmp6 
  double precision      :: tmp7, tmp8
  character(4)          :: tmp9 
  integer               :: tmp10, tmp11
  integer               :: tmp12, tmp13, tmp14
  double precision      :: add_crd(3, proc_num_atms_min_bound)
  double precision      :: add_vel(3, proc_num_atms_min_bound)
  integer               :: add_wrap(3,proc_num_atms_min_bound)
  double precision      :: add_mass(proc_num_atms_min_bound)
  double precision      :: add_qterm(proc_num_atms_min_bound)
  integer               :: add_atm_to_full_list(proc_num_atms_min_bound)
  integer               :: add_iac(proc_num_atms_min_bound)
  !For pme exchange code, bonded data transfer
  double Precision :: crd_x, crd_y, crd_z
!  double precision :: bkt_size_x, bkt_size_y, bkt_size_z
  double precision :: recip_x, recip_y, recip_z
  double precision :: fract_x, fract_y, fract_z
  integer          :: bkt_x, bkt_y, bkt_z
  integer          :: nstep,j 
  integer          :: proc(3), neighbor_proc(3), dest
  
  ! We will need these variables outside the scope of this subroutine while
  ! exchanging bonds between domains, hence copying in global variables
  !the following codes necessary for calculating dest_sdmns() 
!  bkt_size_x = pbc_box(1) / dble(sys_bkt_dimx)
!  bkt_size_y = pbc_box(2) / dble(sys_bkt_dimy)
!  bkt_size_z = pbc_box(3) / dble(sys_bkt_dimz)
  recip_x = 1.d0 / dble(pbc_box(1))
  recip_y = 1.d0 / dble(pbc_box(2))
  recip_z = 1.d0 / dble(pbc_box(3))

  n = 0
  !find which atoms move away from my local space
  do i = 1, proc_num_atms
    if(proc_atm_crd(1, i) .lt. proc_min_x_crd .or. proc_atm_crd(1,i) .ge. proc_max_x_crd &
       .or. proc_atm_crd(2, i) .lt. proc_min_y_crd .or. proc_atm_crd(2,i) .ge. proc_max_y_crd &
       .or. proc_atm_crd(3, i) .lt. proc_min_z_crd .or. proc_atm_crd(3,i) .ge. proc_max_z_crd) then
        n = n + 1
        delete_index(n) = i 
        !calculating dest_sdmns() array

      ! New approach through coordinates information
        call destinationmpi(dest_sdmns(n), proc_atm_crd(:,i))
    end if
  end do
  
  deleted = n
              !since after that proc_atom_to_full_list array has been changed 
  del_cnt = deleted
  atms_to_del = delete_index

 call bonded_packing(proc_atm_crd, neighbor_mpi_cnt,counts,nstep) ! Calling this here is needed
 call proc_save_bonded(counts) ! saving the bonded index's global id, to use later

  n = 0
  !find which atoms came in local space from ghost region
  do i = proc_num_atms + 1, proc_num_atms + proc_ghost_num_atms
    if(proc_atm_crd(1, i) .ge. proc_min_x_crd .and. proc_atm_crd(1,i) .lt. proc_max_x_crd &
       .and. proc_atm_crd(2, i) .ge. proc_min_y_crd .and. proc_atm_crd(2,i) .lt. proc_max_y_crd &
       .and. proc_atm_crd(3, i) .ge. proc_min_z_crd .and. proc_atm_crd(3,i) .lt. proc_max_z_crd) then
        n = n + 1
        add_crd(1,n) = proc_atm_crd(1,i)
        add_crd(2,n) = proc_atm_crd(2,i)
        add_crd(3,n) = proc_atm_crd(3,i)
        add_vel(1,n) = proc_atm_vel(1,i)
        add_vel(2,n) = proc_atm_vel(2,i)
        add_vel(3,n) = proc_atm_vel(3,i)
        add_wrap(1,n) = proc_atm_wrap(1,i)
        add_wrap(2,n) = proc_atm_wrap(2,i)
        add_wrap(3,n) = proc_atm_wrap(3,i)
        ! Because comm send for proc_atm_wrap also updates the wrap counter we
        ! need to double compensate here.
        if(proc_saved_atm_crd(1,i) .ge. pbc_box(1)) add_wrap(1,n) = add_wrap(1,n)-2
        if(proc_saved_atm_crd(2,i) .ge. pbc_box(2)) add_wrap(2,n) = add_wrap(2,n)-2
        if(proc_saved_atm_crd(3,i) .ge. pbc_box(3)) add_wrap(3,n) = add_wrap(3,n)-2
        if(proc_saved_atm_crd(1,i) .lt. 0.d0) add_wrap(1,n) = add_wrap(1,n)+2
        if(proc_saved_atm_crd(2,i) .lt. 0.d0) add_wrap(2,n) = add_wrap(2,n)+2
        if(proc_saved_atm_crd(3,i) .lt. 0.d0) add_wrap(3,n) = add_wrap(3,n)+2
        add_mass(n) = proc_atm_mass(i)
        add_qterm(n) = proc_atm_qterm(i)
        add_atm_to_full_list(n) = proc_atm_to_full_list(i)
        add_iac(n) = proc_iac(i)
    end if
  end do
  added = n
  !delete from my local list and fill the whole
  !by the last atom of local list, and adjust the proc_num_atms
  !do i = 1, deleted 
  !the reason the do loop below is backward, for example: for 10 atoms. 
  !need to delete 3rd and 10th. 
  !now if scan from beginning, 3rd will be replaced by 10th, and 10th will be
  !lost. If you scan from end, you delete 10th, and then 3rd will be replaced by
  !9th. So you have to scan from backward
  do i =  deleted, 1 , -1 
    index = delete_index(i)
    !if the deleted atom is at the end, just remove it
    !by reducing the proc_num_atms
    if(index .eq. proc_num_atms) then
      proc_num_atms = proc_num_atms - 1
      cycle 
    end if
    !if the deleted atom is not the last atom, copy the last atom
    !and paste it at the place of deleted atom(index)
    tmp1 = proc_atm_crd(1, proc_num_atms)
    tmp2 = proc_atm_crd(2, proc_num_atms)
    tmp3 = proc_atm_crd(3, proc_num_atms)
    tmp4 = proc_atm_vel(1, proc_num_atms)
    tmp5 = proc_atm_vel(2, proc_num_atms)
    tmp6 = proc_atm_vel(3, proc_num_atms)
    tmp7 = proc_atm_mass(proc_num_atms)
    tmp8 = proc_atm_qterm(proc_num_atms)
    tmp10 = proc_atm_to_full_list(proc_num_atms)
    tmp11 = proc_iac(proc_num_atms)
    tmp12 = proc_atm_wrap(1, proc_num_atms)
    tmp13 = proc_atm_wrap(2, proc_num_atms)
    tmp14 = proc_atm_wrap(3, proc_num_atms)

    proc_atm_crd(1, index) = tmp1
    proc_atm_crd(2, index) = tmp2
    proc_atm_crd(3, index) = tmp3

    proc_atm_vel(1, index) = tmp4
    proc_atm_vel(2, index) = tmp5
    proc_atm_vel(3, index) = tmp6
    !in this routine proc_atm_vel and proc_atm_last_vel contain same value
    proc_atm_last_vel(1, index) = tmp4
    proc_atm_last_vel(2, index) = tmp5
    proc_atm_last_vel(3, index) = tmp6

    proc_atm_mass(index) = tmp7
    proc_atm_qterm(index) = tmp8
    proc_atm_to_full_list(index) = tmp10
    proc_iac(index) = tmp11
    proc_atm_wrap(1,index) = tmp12
    proc_atm_wrap(2,index) = tmp13
    proc_atm_wrap(3,index) = tmp14
    !since we moved the last atom, recude it 
    proc_num_atms = proc_num_atms - 1 
  end do
  !add the new atoms to the end of my list
  !this may override the ghost list, dont care, 
  !since a fresh ghost list
  ! will be created after atm_distribution
  do i = 1, added
    proc_atm_crd(1,proc_num_atms+i) = add_crd(1, i)
    proc_atm_crd(2,proc_num_atms+i) = add_crd(2, i)
    proc_atm_crd(3,proc_num_atms+i) = add_crd(3, i)
    proc_atm_vel(1,proc_num_atms+i) = add_vel(1, i)
    proc_atm_vel(2,proc_num_atms+i) = add_vel(2, i)
    proc_atm_vel(3,proc_num_atms+i) = add_vel(3, i)
    !in this routine proc_atm_vel and proc_atm_last_vel contain same value
    proc_atm_last_vel(1,proc_num_atms+i) = add_vel(1, i)
    proc_atm_last_vel(2,proc_num_atms+i) = add_vel(2, i)
    proc_atm_last_vel(3,proc_num_atms+i) = add_vel(3, i)
    proc_atm_wrap(1,proc_num_atms+i) = add_wrap(1,i)
    proc_atm_wrap(2,proc_num_atms+i) = add_wrap(2,i)
    proc_atm_wrap(3,proc_num_atms+i) = add_wrap(3,i)
    proc_atm_mass(proc_num_atms+i) = add_mass(i)
    proc_atm_qterm(proc_num_atms+i) = add_qterm(i)
    proc_atm_to_full_list(proc_num_atms+i) = add_atm_to_full_list(i)
    proc_iac(proc_num_atms+i) = add_iac(i)
  end do
  ! updating my local list cnt after atoms are added
  proc_num_atms = proc_num_atms + added
  if(proc_num_atms .lt. proc_num_atms_min_size) then
    proc_num_atms_min_bound = proc_num_atms_min_size
  else
    proc_num_atms_min_bound = proc_num_atms 
  endif

  return

end subroutine proc_add_delete_atms 

!*******************************************************************************
!
! Subroutine:  pme_exchange_midpoint
!
! Description: Move bonds, angles and dihedrals from this subdomain to neighbor
!              subdomain in case thelead atom have moved to the neighbor
!              subdomain.
!              
!*******************************************************************************
! Note:
! 1) Now that dest_sdmn calculation is being done in nb_pairlist.F90, this
! subroutine does not need crd argument.
!
! 2) Send and recv buffers are 2D arrays where each column corresponds to an
! individual send/recv buffer which is mapped to one of the neighbor sub-domain.
! According to Fortran array notation ARRAY(ROW, COL), the base address of send
! and recv buffers should be (1, i).
!*******************************************************************************
#ifdef BEXDBG
subroutine pme_exchange_midpoint(crd, proc_bond_cnt, proc_angle_cnt, proc_dihed_cnt, nbr_sdmn_cnt, nstep)
#else
subroutine pme_exchange_midpoint(crd, proc_bond_cnt, proc_angle_cnt, proc_dihed_cnt, nbr_sdmn_cnt, counts,nstep)
#endif
  use pbc_mod
  use processor_mod
  use bonds_midpoint_mod
  use angles_midpoint_mod
  use dihedrals_midpoint_mod
!  use nb_pairlist_mod, only: del_cnt, atms_to_del, dest_sdmns
  use prmtop_dat_mod, only : natom
  
  implicit none
  
  ! Formal arguments
  ! Array of coordinates
  double precision :: crd(3, *)
  ! Number of neighboring sub-domains
  integer          :: nbr_sdmn_cnt
  ! Number of bonds, angles and dihedrals in this sub-domain
  integer          :: proc_bond_cnt, proc_angle_cnt, proc_dihed_cnt
  integer          :: counts
#ifdef BEXDBG
  ! Only for debug purpose, iteration count
  integer          :: nstep
#endif
  
  ! Local variables
  integer, parameter :: b_mpi_tag = 101, a_mpi_tag = 102, d_mpi_tag = 103
  logical :: del_flag
  integer :: i, j, k, num_ints, atm_to_del
  integer :: b_lead_atm_h, b_lead_atm_a
  integer :: a_lead_atm_h, a_lead_atm_a
  integer :: d_lead_atm_h, d_lead_atm_a
  integer :: src_sdmn, dest_sdmn, src_recv_cnt, dest_send_cnt, src_idx, dest_idx
  integer, dimension(nbr_sdmn_cnt) :: b_recv_req, b_send_req
  integer, dimension(nbr_sdmn_cnt) :: a_recv_req, a_send_req
  integer, dimension(nbr_sdmn_cnt) :: d_recv_req, d_send_req
 ! integer, dimension(2, 0:numtasks-1) :: b_send_map, a_send_map, d_send_map
  integer, dimension(mpi_status_size) :: b_recv_stat, a_recv_stat, d_recv_stat
  integer, dimension(mpi_status_size, nbr_sdmn_cnt) :: b_send_stat, a_send_stat, d_send_stat
! type(bond_trans), dimension(proc_bond_cnt, nbr_sdmn_cnt) :: b_recv_buf, b_send_buf
! type(angle_trans), dimension(proc_angle_cnt, nbr_sdmn_cnt) :: a_recv_buf, a_send_buf
! type(dihed_trans), dimension(proc_dihed_cnt, nbr_sdmn_cnt) :: d_recv_buf, d_send_buf
  double Precision :: crd_x, crd_y, crd_z
  double precision :: bkt_size_x, bkt_size_y, bkt_size_z
  double precision :: recip_x, recip_y, recip_z
  double precision :: fract_x, fract_y, fract_z
  integer          :: bkt_x, bkt_y, bkt_z
  integer          :: atmid,src_mpi
  integer          :: nstep
#ifdef BEXDBG
  integer          :: non_lead_atm_del_cnt=0, self_dest_cnt=0
  integer          :: b_send_cnt=0, a_send_cnt=0, d_send_cnt=0
  integer          :: fileunit
  character(len=32):: dbgfile
#endif
  !NOTE: No need of receieve map because recv_req arrays themselves act as a map
  ! Initialize the send map to zero
!  b_send_map(:,:) = 0
!  a_send_map(:,:) = 0
!  d_send_map(:,:) = 0
!Experimental
  proc_atm_space(:) = 0 !very expensive, need to only reset the previous nb steps
! proc_atm_space_ghosts(:) = 0!very expensive, need to only reset the previous nb steps
  do i = proc_num_atms+1, proc_num_atms + proc_ghost_num_atms
         atmid = proc_atm_to_full_list(i)
         proc_atm_space(atmid) = i
  end do
  do i = 1, proc_num_atms
         atmid = proc_atm_to_full_list(i)
         proc_atm_space(atmid) = i
  end do

  ! Post asynchronous receives for incoming bonds, angles and dihedrals.
  ! NOTE: Each rank has a receive buffer just enough to store as many
  ! bonds, angles and dihedrals as it already has.
  ! Consider, this rank has X bonds and the neighbor rank has to send Y
  ! bonds to this rank. It's highly unlikey that Y > X, but this design
  ! won't be able to handle such scenario.
  do i = 1, nbr_sdmn_cnt
    src_sdmn = neighbor_mpi_map(i)
    num_ints = max(proc_bond_cnt,proc_num_atms_min_bound) * bond_trans_ints
    call mpi_irecv(b_recv_buf(1, i), num_ints, mpi_integer, src_sdmn,     &
                   b_mpi_tag, pmemd_comm, b_recv_req(i), err_code_mpi); MPI_ERR_CHK("pme_exchange irecv bond")
    
    num_ints = max (proc_angle_cnt,proc_num_atms_min_bound) * angle_trans_ints
    call mpi_irecv(a_recv_buf(1, i), num_ints, mpi_integer, src_sdmn,      &
                   a_mpi_tag, pmemd_comm, a_recv_req(i), err_code_mpi); MPI_ERR_CHK("pme_exchange irecv angle")
    
    num_ints = max(proc_dihed_cnt,proc_num_atms_min_bound) * dihed_trans_ints
    call mpi_irecv(d_recv_buf(1, i), num_ints, mpi_integer, src_sdmn,      &
                   d_mpi_tag, pmemd_comm, d_recv_req(i), err_code_mpi); MPI_ERR_CHK("pme_exchange irecv dihed")
    
  end do
  
  ! Go through the lead atoms moving out of this sub-domain
  ! Post asynchronous sends for outgoing bonds, angles and dihedrals
  do i = 1, nbr_sdmn_cnt
    dest_sdmn = neighbor_mpi_map(i)
    dest_send_cnt = b_send_map(2, dest_sdmn)
    num_ints = dest_send_cnt * bond_trans_ints
    call mpi_isend(b_send_buf(1, i), num_ints, mpi_integer, dest_sdmn,         &
                   b_mpi_tag, pmemd_comm, b_send_req(i), err_code_mpi); MPI_ERR_CHK("pme_exchange isend bond")
    dest_send_cnt = a_send_map(2, dest_sdmn)
    num_ints = dest_send_cnt * angle_trans_ints
    call mpi_isend(a_send_buf(1, i), num_ints, mpi_integer, dest_sdmn,         &
                   a_mpi_tag, pmemd_comm, a_send_req(i), err_code_mpi); MPI_ERR_CHK("pme_exchange isend angle")
    
    dest_send_cnt = d_send_map(2, dest_sdmn)
    num_ints = dest_send_cnt * dihed_trans_ints
    call mpi_isend(d_send_buf(1, i), num_ints, mpi_integer, dest_sdmn,         &
                   d_mpi_tag, pmemd_comm, d_send_req(i), err_code_mpi); MPI_ERR_CHK("pme_exchange isend dihed")
  end do
  
  ! Wait for all receives to finish
  do i = 1, nbr_sdmn_cnt
    ! Wait to receive bonds, unpack and add
    call mpi_waitany(nbr_sdmn_cnt, b_recv_req, src_idx, b_recv_stat, err_code_mpi); MPI_ERR_CHK("pme_exchange wait bond")
    call mpi_get_count(b_recv_stat, mpi_integer, num_ints, err_code_mpi); MPI_ERR_CHK("pme_exchange get count bond")
    src_recv_cnt = int(num_ints / bond_trans_ints)
    src_mpi = b_recv_stat(mpi_source)

    ! Unpack bonds and add to list after converting global indices to local
    do j = 1, src_recv_cnt
      if (b_recv_buf(j, src_idx)%hatype .eq. 0) then
        cit_nbonh = cit_nbonh + 1
        cit_h_bond(cit_nbonh)%atm_i = proc_atm_space(b_recv_buf(j, src_idx)%atm_i)
        cit_h_bond(cit_nbonh)%atm_j = proc_atm_space(b_recv_buf(j, src_idx)%atm_j)
        cit_h_bond(cit_nbonh)%parm_idx = b_recv_buf(j, src_idx)%parm_idx
        my_hbonds_leads(cit_nbonh) = cit_h_bond(cit_nbonh)%atm_i
      else
        cit_nbona = cit_nbona + 1
        cit_a_bond(cit_nbona)%atm_i = proc_atm_space(b_recv_buf(j, src_idx)%atm_i)
        cit_a_bond(cit_nbona)%atm_j = proc_atm_space(b_recv_buf(j, src_idx)%atm_j)
        cit_a_bond(cit_nbona)%parm_idx = b_recv_buf(j, src_idx)%parm_idx
        my_abonds_leads(cit_nbona) = cit_a_bond(cit_nbona)%atm_i
      end if
    end do
    
    ! Wait to receive angles, unpack and add
    call mpi_waitany(nbr_sdmn_cnt, a_recv_req, src_idx, a_recv_stat, err_code_mpi); MPI_ERR_CHK("pme_exchange wait angle")
    call mpi_get_count(a_recv_stat, mpi_integer, num_ints, err_code_mpi); MPI_ERR_CHK("pme_exchange cout angle")
    src_recv_cnt = int(num_ints / angle_trans_ints)
    src_mpi = a_recv_stat(mpi_source)
    ! Unpack and add to list after converting global indices to local
    do j = 1, src_recv_cnt
      if (a_recv_buf(j, src_idx)%hatype .eq. 0) then
        cit_ntheth = cit_ntheth + 1
        cit_h_angle(cit_ntheth)%atm_i = proc_atm_space(a_recv_buf(j, src_idx)%atm_i)
        cit_h_angle(cit_ntheth)%atm_j = proc_atm_space(a_recv_buf(j, src_idx)%atm_j)
        cit_h_angle(cit_ntheth)%atm_k = proc_atm_space(a_recv_buf(j, src_idx)%atm_k)
        cit_h_angle(cit_ntheth)%parm_idx = a_recv_buf(j, src_idx)%parm_idx
        my_hangles_leads(cit_ntheth) = cit_h_angle(cit_ntheth)%atm_j
      else
        cit_ntheta = cit_ntheta + 1
        cit_a_angle(cit_ntheta)%atm_i = proc_atm_space(a_recv_buf(j, src_idx)%atm_i)
        cit_a_angle(cit_ntheta)%atm_j = proc_atm_space(a_recv_buf(j, src_idx)%atm_j)
        cit_a_angle(cit_ntheta)%atm_k = proc_atm_space(a_recv_buf(j, src_idx)%atm_k)
        cit_a_angle(cit_ntheta)%parm_idx = a_recv_buf(j, src_idx)%parm_idx
        my_aangles_leads(cit_ntheta) = cit_a_angle(cit_ntheta)%atm_j
      end if
    end do
    
    ! Wait to receive dihedrals, unpack and add
    call mpi_waitany(nbr_sdmn_cnt, d_recv_req, src_idx, d_recv_stat, err_code_mpi); MPI_ERR_CHK("pme_exchange wait dihed")
    call mpi_get_count(d_recv_stat, mpi_integer, num_ints, err_code_mpi); MPI_ERR_CHK("pme_exchange count dihed")
    src_recv_cnt = int(num_ints / dihed_trans_ints)
    src_mpi = d_recv_stat(mpi_source)
    
    ! Unpack and add to list after converting global indices to local
    do j = 1, src_recv_cnt
      if (d_recv_buf(j, src_idx)%hatype .eq. 0) then
        cit_nphih = cit_nphih + 1
        cit_h_dihed(cit_nphih)%atm_i = proc_atm_space(d_recv_buf(j, src_idx)%atm_i)
        cit_h_dihed(cit_nphih)%atm_j = proc_atm_space(d_recv_buf(j, src_idx)%atm_j)
        cit_h_dihed(cit_nphih)%atm_k = proc_atm_space(iabs(d_recv_buf(j, src_idx)%atm_k))*sign(1,d_recv_buf(j,src_idx)%atm_k)
        cit_h_dihed(cit_nphih)%atm_l = proc_atm_space(iabs(d_recv_buf(j, src_idx)%atm_l))*sign(1,d_recv_buf(j,src_idx)%atm_l)
        cit_h_dihed(cit_nphih)%parm_idx = d_recv_buf(j, src_idx)%parm_idx
        my_hdiheds_leads(cit_nphih) = cit_h_dihed(cit_nphih)%atm_j
      else
        cit_nphia = cit_nphia + 1
        cit_a_dihed(cit_nphia)%atm_i = proc_atm_space(d_recv_buf(j, src_idx)%atm_i)
        cit_a_dihed(cit_nphia)%atm_j = proc_atm_space(d_recv_buf(j, src_idx)%atm_j)
        cit_a_dihed(cit_nphia)%atm_k = proc_atm_space(iabs(d_recv_buf(j, src_idx)%atm_k))*sign(1,d_recv_buf(j,src_idx)%atm_k)
        cit_a_dihed(cit_nphia)%atm_l = proc_atm_space(iabs(d_recv_buf(j, src_idx)%atm_l))*sign(1,d_recv_buf(j,src_idx)%atm_l)
        cit_a_dihed(cit_nphia)%parm_idx = d_recv_buf(j, src_idx)%parm_idx
        my_adiheds_leads(cit_nphia) = cit_a_dihed(cit_nphia)%atm_j
      end if
    end do
  end do

  ! Wait for all sends to finish
  call mpi_waitall(nbr_sdmn_cnt, b_send_req, b_send_stat, err_code_mpi); MPI_ERR_CHK("pme_exchange wait send bond")
  call mpi_waitall(nbr_sdmn_cnt, a_send_req, a_send_stat, err_code_mpi); MPI_ERR_CHK("pme_exchange wait send angle")
  call mpi_waitall(nbr_sdmn_cnt, d_send_req, d_send_stat, err_code_mpi); MPI_ERR_CHK("pme_exchange wait send dihed")
#ifdef BEXDBG
  write(fileunit, "('Exchange['I2', 'I5']: All send requests complete.')"), mytaskid, nstep
#endif
end subroutine pme_exchange_midpoint


!*******************************************************************************
!
! Subroutine:   atom_sorting
!
! Description:  Routine sorting local atoms, crd, charge, iac, full list
!
!*******************************************************************************
subroutine atom_sorting()

  use cit_mod
  use gbl_datatypes_mod
  use img_mod
  use parallel_dat_mod
  use pbc_mod
  use ti_mod
  use processor_mod

  implicit none

!local variables
  integer               :: x, y, z !3d co-ordinate of bucket 
  double precision      :: tmp_atm_crd(3,proc_num_atms_min_bound), tmp_qterm(proc_num_atms_min_bound)
  double precision      :: tmp_atm_frc(3,proc_num_atms_min_bound)
  integer               :: tmp_send_indices(proc_num_atms_min_bound)
  integer               :: tmp_old_local_id(proc_num_atms_min_bound)
  double precision      :: tmp_atm_vel(3,proc_num_atms_min_bound),  tmp_atm_last_vel(3,proc_num_atms_min_bound)
  integer               :: tmp_iac(proc_num_atms_min_bound),tmp_atm_to_full_list(proc_num_atms_min_bound)
  double precision      :: tmp_mass(proc_num_atms_min_bound)
  double precision      :: tmp_fraction(3,proc_num_atms_min_bound) 
  integer               :: tmp_wrap(3,proc_num_atms_min_bound)
  integer               :: i,j, bkt, atm_i_idx, atm_i
  integer               :: tmp_bkt_atm_lst(size(bkt_atm_lst,1))

  do i = 1, proc_num_atms
    tmp_atm_crd(1:3,i) = proc_atm_crd(1:3,i)
    tmp_wrap(1:3,i) = proc_atm_wrap(1:3,i)
    tmp_atm_frc(1:3,i) = proc_atm_frc(1:3,i)
    tmp_atm_vel(1:3,i) = proc_atm_vel(1:3,i)
    tmp_atm_last_vel(1:3,i) = proc_atm_last_vel(1:3,i)
    tmp_mass(i)        = proc_atm_mass(i)
    tmp_qterm(i)       = proc_atm_qterm(i)
    tmp_old_local_id(i)       = proc_old_local_id(i)
    tmp_iac(i)         = proc_iac(i)
    tmp_atm_to_full_list(i)= proc_atm_to_full_list(i)
    tmp_fraction(1:3,i) = proc_fraction(1:3,i)
  end do

  j = 1
  do z = 1, int_bkt_dimz
    do y= 1, int_bkt_dimy
      do x= 1, int_bkt_dimx

        bkt=int_flatten(x,y,z)
        do atm_i_idx =1,bkt_atm_cnt(bkt)
          tmp_bkt_atm_lst(atm_i_idx) = bkt_atm_lst(atm_i_idx,bkt)
        enddo
 
        do atm_i_idx =1,bkt_atm_cnt(bkt)
          i = tmp_bkt_atm_lst(atm_i_idx)
          if(i .le. proc_num_atms) then
            proc_atm_crd(1:3,j) = tmp_atm_crd(1:3,i) 
            proc_atm_wrap(1:3,j) =  tmp_wrap(1:3,i) 
            proc_atm_frc(1:3,j) = tmp_atm_frc(1:3,i) 
            proc_atm_vel(1:3,j) = tmp_atm_vel(1:3,i) 
            proc_atm_last_vel(1:3,j) = tmp_atm_last_vel(1:3,i) 
            proc_atm_qterm(j) = tmp_qterm(i)
            proc_atm_mass(j) = tmp_mass(i) 
            proc_old_local_id(j) = tmp_old_local_id(i) 
            proc_iac(j) = tmp_iac(i)
            proc_atm_to_full_list(j) = tmp_atm_to_full_list(i)
            proc_sort_map(i) = j
            proc_fraction(1:3,j) = tmp_fraction(1:3,i)
            bkt_atm_lst(atm_i_idx, bkt) = j
            j=j+1
          else
            proc_sort_map(i) = i
          endif
        end do

    end do
  end do
end do


  ! Sort the new communication method ghost atoms indices
  do i=1,6
    tmp_send_indices(1:send_buf_cnt(i)) = send_buf_indices(1:send_buf_cnt(i), i)
    do j=1,send_buf_cnt(i)
      if(tmp_send_indices(j) .gt. proc_num_atms) exit ! we do not sort ghost atoms
!if(proc_sort_map(tmp_send_indices(j)) == 0 ) print*, "index map 0:", i, j, tmp_send_indices(j)
      send_buf_indices(j, i) = proc_sort_map(tmp_send_indices(j))
    end do
  end do

  return
end subroutine atom_sorting

#include "get_nb_energy_midpoint.i"
#endif /*MPI*/

#include "get_nb_energy.i"
#include "get_nb_energy_offload.i"
#include "get_nb_energy_ti_AUTO.i"
#include "get_nb_energy_ips.i"
#include "get_nb_energy_lj1264.i"
#include "get_nb_energy_plj1264.i" ! New2022
#include "get_nb_energy_lj1264plj1264.i" ! New2022

end module pme_direct_mod
