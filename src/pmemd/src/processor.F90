#include "copyright.i"
#include "include_precision.i"

!*******************************************************************************
!
! Module: processor_mod
!
! Description: Generates processor specific information and processor grid
!              
!*******************************************************************************

module processor_mod

#ifdef MPI

use parallel_dat_mod
use energy_records_mod, only : pme_pot_ene_rec

#ifdef _OPENMP_
use omp_lib
#endif

  implicit none

  !Global Variables
  !Block variables
  !# of MPI ranks = number of blocks=sudomains=proc_dimx*proc_dimy*proc_dimz
  !proc_dimx, proc_dimy, proc_dimz divides the domain in 3 dimensions
  !these are the division, not length
  double Precision, save              :: new_bkt_size(3)
  integer, save              :: proc_dimx, proc_dimy, proc_dimz, proc_max_num_ghost_atm
  integer, save              :: proc_dimx_recip, proc_dimy_recip, proc_dimz_recip
  integer, save              :: myblockid

  !Bucket (Bkt) Variables
  !CIT
  !Bucket dimensions for entire system
  !Individual
  !beginning and ending bucket numbers (not length) in each dim owned by one MPI
  integer, save              :: proc_bkt_minx, proc_bkt_maxx, proc_bkt_miny, &
                                proc_bkt_maxy, proc_bkt_minz, proc_bkt_maxz
  !internal bucket numbering without including ghost portion
  integer, save              :: int_bkt_minx, int_bkt_maxx, int_bkt_miny, & 
                                int_bkt_maxy, int_bkt_minz, int_bkt_maxz
  integer, save              :: bkt_dim(3)
  integer, save              :: int_bkt_dimx, int_bkt_dimy, int_bkt_dimz
  integer                    :: proc_numbkts, proc_numghostbkts! per MPI 
  integer                    :: number_buckets_in_block !The number of buckets in the min max range
  integer, allocatable, save :: proc_range_ghostbkts(:,:)
  integer, allocatable, save :: proc_temp_ghostbkts(:)
  integer, allocatable, save :: proc_corner_neighborbuckets(:,:) !Added, contains the range of corner bucketids
  integer, allocatable, save :: proc_corner_neighborbuckets_cnt(:) !Added, contains the number of corner bucketids
  integer, allocatable, save :: send_bucket_mpi_list(:,:)

  integer, save              :: neighbor_mpi_map(27),neighbor_mpi_cnt
#if 0
  integer, allocatable, save :: block_minmax_array(:,:) !Contains the minmax ranges for each block(task)
#endif
  integer, allocatable, save :: int_ext_bkt(:) !Contains int to ext bkt mapping
  integer, allocatable, save :: int_mpi_map(:) !Contains if internal bucket is self or not
! each bkt is perfectly cube, with a cut_lst lenght in each 3 dimension
  
  !for bucketing fractiona_proc is needed
  !populating proc_fraction is written in pbc.F90 file
  double precision, allocatable           :: proc_fraction(:,:)
#if defined(pmemd_SPDP)
  real, allocatable           :: proc_fraction_sp(:,:)
#endif

  ! Atom min max
  ! acutal co-ordinate of each block
  double precision           :: proc_crd_bound(2,3) !(min/max, x/y/z)
  double precision           :: proc_crd_len(3), ghost_crd_len(3)
  double precision, allocatable :: crd_boundx_recip(:,:),crd_boundy_recip(:,:),crd_boundz_recip(:,:) !(min/max, MPI rank)
  double precision           :: proc_min_x_crd, proc_max_x_crd
  double precision           :: proc_min_y_crd, proc_max_y_crd
  double precision           :: proc_min_z_crd, proc_max_z_crd
  integer                    :: proc_num_atms, proc_num_atms_min_bound
  integer                    :: proc_atm_alloc_size
  integer,parameter          :: proc_num_atms_min_size=10

  ! Standard Atom Arrays
  double precision, allocatable :: proc_atm_crd(:,:)

#if defined(pmemd_SPDP)
  real, allocatable :: proc_crd_q_sp(:,:)
  real, allocatable, save, public :: proc_gbl_cn1_sp(:)
  real, allocatable, save, public :: proc_gbl_cn2_sp(:)
#endif /*pmemd_SPDP*/
  double precision, allocatable :: proc_saved_atm_crd(:,:)
  double precision, allocatable :: proc_atm_vel(:,:)
  double precision, allocatable :: proc_atm_last_vel(:,:)
  double precision, allocatable :: proc_atm_frc(:,:)
  double precision, allocatable :: proc_atm_mass(:)
  !DIR$ ATTRIBUTES ALIGN : 64 :: proc_atm_qterm
  double precision, allocatable :: proc_atm_qterm(:)
  character(4), allocatable     :: proc_atm_igraph(:)
  !DIR$ ATTRIBUTES ALIGN : 64 :: proc_atm_to_full_list
  integer, allocatable          :: proc_atm_to_full_list(:)
  integer, allocatable          :: old_proc_atm_to_full_list(:)
! proc_atm_to_full_list takes in the proc index for the atom and gives the index
! that the inpcrd has for the atom
  integer, allocatable          :: proc_atm_wrap(:,:)
! The idea of atm_wrap_counter is if you proc_atm_crd +
! proc_atm_wrap you'll get the coordinate in real space.
  integer, allocatable          :: proc_old_local_id(:)
  integer, allocatable          :: proc_old_mpi(:) 
  integer, allocatable          :: proc_sort_map(:)
  
  !Energy and related arrays ! later we can reuse from the actual code
  
  ! Parameter arrays
  ! proc_iac allows us to find the ids needed to find the parameters for an atom
  ! pair.  So assume 2 atms that we need to calculate force for proc_atm_i and proc_atm_j.
  ! iaci = ntypes_stk * (proc_iac(proc_atm_i)-1) (flattens one to be length
  ! ntypes)
  ! ic = ico(iaci  + proc_iac(proc_atm_j))
  ! Then we can use cn1(ic) and cn2(ic) to get the LJ parameters for this atom
  ! pair. cn1, cn2, and ico do not need separate per processor arrays because
  ! they're sized ntypes x ntypes.
  integer, allocatable           :: proc_iac(:)

  !GB Atom Arrays
  integer, allocatable          :: proc_atm_atomicnumber(:)
  character(4), allocatable     :: proc_atm_isymbl(:)
  double precision, allocatable :: proc_atm_gb_radii(:)
  double precision, allocatable :: proc_atm_gb_fs(:)

  !NEW Ghost Atom Arrays
  ! TERMENOLOGY: The arrays related to the 6 sides of the box are
  !   stored in the order: x-, x+, y-, y+, z-, z+
  integer, allocatable :: dbg_atm_to_full_list(:), dbg_ghost_permute(:)
  double precision, allocatable :: dbg_atm_crd(:,:)
  integer, allocatable :: send_buf_indices(:,:)
  integer              :: recv_buf_idx_range(2,6), proc_neighbor_rank(6)
  integer              :: send_buf_cnt(6), recv_buf_cnt(6)
  double precision     :: proc_send_bounds(6), proc_side_wrap_offset(3,6)
  double precision, allocatable :: mpi_send_buf_dble(:,:),mpi_recv_buf_dble(:,:)
  integer, allocatable :: mpi_send_buf_int(:,:),mpi_recv_buf_int(:,:)

  integer              :: ghost_bkt_max_size, ghost_bkt_max_size_max_proc

  integer              :: proc_ghost_num_atms

  ! Output MPI Transfer
  double precision, allocatable   :: send_output_dbl(:,:) ! send int array
  integer, allocatable            :: send_output_int(:,:) ! send double array
  double precision, allocatable   :: recv_output_dbl(:,:) ! recv int array
  integer, allocatable            :: recv_output_int(:,:) ! recv  double array

  integer                         :: frc_mult   ! multiplier needed for array
                                                ! size due to ghosts for force

  ! PME Variables & Arrays
  ! Gives min max crds for pme
  double precision                :: proc_pme_minx, proc_pme_maxx, &
                                     proc_pme_miny, proc_pme_maxy, &
                                     proc_pme_minz, proc_pme_maxz
  integer                         :: proc_min_nfft1, proc_max_nfft1, &
                                     proc_min_nfft2, proc_max_nfft2, &
                                     proc_min_nfft3, proc_max_nfft3
  integer                         :: fft_halo_size
  integer, allocatable            :: send_xfft_cnts(:)
  integer, allocatable            :: send_xfft_offset(:)
  integer, allocatable            :: recv_xfft_cnts(:)
  integer, allocatable            :: recv_xfft_offset(:)
  integer, allocatable            :: send_yfft_cnts(:)
  integer, allocatable            :: send_yfft_offset(:)
  integer, allocatable            :: recv_yfft_cnts(:)
  integer, allocatable            :: recv_yfft_offset(:)
  integer, allocatable            :: send_zfft_cnts(:)
  integer, allocatable            :: send_zfft_offset(:)
  integer, allocatable            :: recv_zfft_cnts(:)
  integer, allocatable            :: recv_zfft_offset(:)
  integer, allocatable            :: x_offset(:)
  integer, allocatable            :: y_offset(:)
  integer, allocatable            :: z_offset(:)
!max cnts so that mpi_buf can be calculated
  integer                         :: max_xfft_cnts, max_yfft_cnts, max_zfft_cnts
  integer                         :: my_xfft_cnt, my_yfft_cnt, my_zfft_cnt
  integer                         :: my_mpix, my_mpiy, my_mpiz !my 3d co-ordinates
  
  ! Gives blocks (mpi tasks) in x, xy, or z directions not including self block
  integer, allocatable            :: proc_x_comm_lst(:), proc_xy_comm_lst(:), &
                                     proc_y_comm_lst(:), proc_z_comm_lst(:)
  integer, allocatable            :: proc_charge_list(:)
  integer                         :: proc_num_charge_atms
  !bond related arrays
  integer, allocatable, save            :: proc_atm_space(:)
  integer, allocatable, save            :: proc_shake_space(:)
  integer, allocatable, save            :: proc_atm_space_ghosts(:)
  integer, allocatable, save            :: mult_vech_bond(:,:)
  integer, allocatable, save            :: mult_veca_bond(:,:)
  integer, allocatable, save            :: mult_vech_angle(:,:)
  integer, allocatable, save            :: mult_veca_angle(:,:)
  integer, allocatable, save            :: mult_vech_dihed(:,:)
  integer, allocatable, save            :: mult_veca_dihed(:,:)
  integer, allocatable, save            :: index_14(:,:)
  integer, allocatable, save            :: mult_vec_adjust(:,:)
  integer, save                         :: del_cnt
  integer, allocatable, save            :: atms_to_del(:)
  integer, allocatable, save            :: dest_sdmns(:)
  integer, allocatable,save             :: b_send_map(:,:), a_send_map(:,:), d_send_map(:,:)
  double precision, allocatable, save   :: box_len(:)!other names for pbc_box()
  integer,allocatable,save              :: excl_img_flags(:) !it can be byte size
  integer,allocatable,save              :: start_index(:) !it can be byte size
#ifdef _OPENMP_
  ! Used for OMP in direct force calculations.
  double precision,allocatable,save     :: frc_thread(:,:) 
  integer, save                         ::  nthreads
#endif /*_OPENMP_ */
  pme_float, save                       :: neg2InvSqrtPi
  pme_float, parameter                  :: PI = 3.1415926535897930

  !DIR$ ATTRIBUTES ALIGN : 64 :: atm_j_excl_list, bkt_atm_lst
  integer,allocatable,save              :: bkt_atm_lst(:,:) ! list of atoms in each bucket
  integer,allocatable,save              :: bkt_atm_cnt(:) ! number if atoms in each bucket
  integer,allocatable,save              :: atm_j_excl_list(:)
  logical                               :: reorient_flag

contains

!*******************************************************************************
!
! Subroutine: adjust_pbc
!
! Description: adjusts atoms and bucket widths
!
!*******************************************************************************

!*******************************************************************************
!
! Subroutine: proc_data_dump
!
! Description: debug routine
!
!*******************************************************************************

subroutine proc_data_dump(pbc_box, natom)

  implicit none

  double precision, intent(in) ::  pbc_box(3)
  integer debugout
  integer i
  double precision r11,r22,r33
  double precision f1, f2, f3
  integer, intent(in) :: natom
  debugout = 15+myblockid  

  write(debugout, *) "Proc Dim: ",proc_dimx,' ',proc_dimy,' ', proc_dimz," Myblockid: ",myblockid

  write(debugout, '(A,I5)')'numatoms ',proc_num_atms
  write(debugout, '(A,F6.3,1x,F6.3,1x,F6.3)')"PBC (A): ",pbc_box(1),pbc_box(2),pbc_box(3)
  write(debugout, '(A,F6.3,1x,F6.3,A,F6.3,1x,F6.3,A,F6.3,1x,F6.3)') &
       'X: ', proc_min_x_crd, proc_max_x_crd, ' Y: ', &
       proc_min_y_crd, proc_max_y_crd, ' Z: ', &
       proc_min_z_crd, proc_max_z_crd
!  write(debugout, '(A,I5,1x,I5,1x,I5)') 'Sys bkt dim: ',sys_bkt_dimx,sys_bkt_dimy,sys_bkt_dimz
  write(debugout, '(A,I5,1x,I5,1x,I5)') 'Proc bkt dim: ',int_bkt_dimx,int_bkt_dimy,int_bkt_dimz
  write(debugout, '(A,I5,1x,I5,1x,I5,1x,I5,1x,I5,1x,I5)')'bucket ranges: ', &
   proc_bkt_minx, proc_bkt_maxx, proc_bkt_miny, proc_bkt_maxy, proc_bkt_minz, proc_bkt_maxz

!  r11=1.d0/(dble(int_bkt_dimx + 2)*pbc_box(1)/dble(sys_bkt_dimx))
!  r22=1.d0/(dble(int_bkt_dimy + 2)*pbc_box(2)/dble(sys_bkt_dimy))
!  r33=1.d0/(dble(int_bkt_dimz + 2)*pbc_box(3)/dble(sys_bkt_dimz))

!  do i=1, natom
!        f1=r11 * (proc_atm_crd(1, i) - proc_min_x_crd + pbc_box(1)/dble(sys_bkt_dimx))
!        f2=r22 * (proc_atm_crd(2, i) - proc_min_y_crd + pbc_box(2)/dble(sys_bkt_dimy))
!        f3=r33 * (proc_atm_crd(3, i) - proc_min_z_crd + pbc_box(3)/dble(sys_bkt_dimz))
!        write(debugout, '(I5,1x,F6.3,1x,F6.3,1x,F6.3,1x,A,1x,I5)') i, &
!          proc_atm_frc(1,i), proc_atm_frc(2,i),proc_atm_frc(3,i),'fullid',proc_atm_to_full_list(i)
!  end do

end subroutine proc_data_dump


!*******************************************************************************
!
! Subroutine:  processor_cleanup
!
! Description: Deallocates stuff 
!              
!*******************************************************************************

subroutine processor_cleanup

  implicit none
  integer ier_alloc
  
  if(allocated(proc_atm_crd)) deallocate (proc_atm_crd,stat=ier_alloc)
#if defined(pmemd_SPDP)
  if(allocated(proc_crd_q_sp)) deallocate (proc_crd_q_sp,stat=ier_alloc)
#endif
  if(allocated(proc_atm_vel)) deallocate (proc_atm_vel,stat=ier_alloc)
  if(allocated(proc_atm_frc)) deallocate (proc_atm_frc,stat=ier_alloc)
  if(allocated(proc_atm_last_vel)) deallocate (proc_atm_last_vel,stat=ier_alloc)
  if(allocated(proc_atm_to_full_list)) &
             deallocate(proc_atm_to_full_list,stat=ier_alloc)
  if(allocated(proc_iac)) deallocate (proc_iac,stat=ier_alloc)
  if(allocated(proc_atm_mass)) deallocate (proc_atm_mass,stat=ier_alloc)
  if(allocated(proc_atm_qterm)) deallocate (proc_atm_qterm,stat=ier_alloc)
  if(allocated(proc_atm_igraph)) deallocate (proc_atm_igraph,stat=ier_alloc)
  if(allocated(proc_atm_atomicnumber)) deallocate (proc_atm_atomicnumber,stat=ier_alloc)
  if(allocated(proc_atm_isymbl)) deallocate (proc_atm_isymbl,stat=ier_alloc)
  if(allocated(proc_atm_gb_radii)) deallocate (proc_atm_gb_radii,stat=ier_alloc)
  if(allocated(proc_atm_gb_fs)) deallocate (proc_atm_gb_fs,stat=ier_alloc)
 
 
  return

end subroutine

!*******************************************************************************
!
! Subroutine: ensure_arrays_with_ghost
!
! Description: moves atm_crd, atm_vel, atm_last_vel, atm_frc, atm_iac, atm_mass,
!              atm_atomicnumber, atm_qterm, atm_igraph
!
!*******************************************************************************
subroutine ensure_arrays_with_ghost(alloc_size)
  use ensure_alloc_mod
  implicit none

  integer alloc_size

  integer,save :: numDouEle=9  ! 3 coord, one mass one qterm. Change these as needed for allocation
  integer atm_cnt

  proc_atm_alloc_size = alloc_size

 
#if defined(pmemd_SPDP)
  call ensure_alloc_float4(proc_crd_q_sp, alloc_size, .true.)
  call ensure_alloc_float3(proc_fraction_sp, alloc_size, .true.)
#endif
  call ensure_alloc_dble3(proc_saved_atm_crd, alloc_size, .true.)
  call ensure_alloc_dble3(proc_atm_frc, alloc_size, .true.)
  call ensure_alloc_dble3(proc_atm_last_vel, alloc_size, .true.)
  call ensure_alloc_dble3(proc_atm_vel, alloc_size, .true.)
  call ensure_alloc_int(proc_atm_to_full_list, alloc_size, .true.)
  call ensure_alloc_int(old_proc_atm_to_full_list, alloc_size, .true.)
  call ensure_alloc_char4(proc_atm_igraph, alloc_size, .true.)
  call ensure_alloc_char4(proc_atm_isymbl, alloc_size, .true.)
  call ensure_alloc_dble(proc_atm_mass, alloc_size, .true.)
  call ensure_alloc_dble(proc_atm_qterm, alloc_size, .true.)
  call ensure_alloc_int3(proc_atm_wrap, alloc_size, .true.)
  call ensure_alloc_dble3(proc_fraction, alloc_size, .true.)
  call ensure_alloc_int(proc_iac, alloc_size, .true.)
  call ensure_alloc_int(proc_old_local_id, alloc_size, .true.)
  call ensure_alloc_int(proc_old_mpi, alloc_size, .true.)
  call ensure_alloc_int(proc_sort_map, alloc_size, .true.)

  call ensure_alloc_int3(mult_vech_bond, alloc_size, .true.)
  call ensure_alloc_int3(mult_veca_bond, alloc_size, .true.)
  call ensure_alloc_int3(mult_vech_angle, alloc_size, .true.)
  call ensure_alloc_int3(mult_veca_angle, alloc_size, .true.)
  call ensure_alloc_int3(mult_vech_dihed, 3*alloc_size, .true.)
  call ensure_alloc_int3(mult_veca_dihed, 3*alloc_size, .true.)
  call ensure_alloc_int3(index_14, alloc_size, .true.)
  call ensure_alloc_int3(mult_vec_adjust, alloc_size*10)
  call ensure_alloc_int(excl_img_flags, numtasks*alloc_size)
  call ensure_alloc_int(start_index, numtasks*alloc_size)
  !output mpi transfer related
  if(master) then
    call ensure_alloc_int2d(recv_output_int, numDouEle*alloc_size+2, numtasks)
    call ensure_alloc_dble2d(recv_output_dbl, numDouEle*alloc_size+2, numtasks)
  end if

end subroutine ensure_arrays_with_ghost


!*******************************************************************************
!
! Subroutine: alloc_full_list_to_proc
!
! Description: moves atm_crd, atm_vel, atm_last_vel, atm_frc, atm_iac, atm_mass,
!              atm_atomicnumber, atm_qterm, atm_igraph
!
!*******************************************************************************

subroutine alloc_full_list_to_proc(atm_crd,atm_frc,atm_vel,atm_last_vel,natoms,pbc_box)

  use prmtop_dat_mod
  use mdin_ctrl_dat_mod
  use ensure_alloc_mod
  implicit none
  
  integer i,j,natoms
  double precision, allocatable :: atm_crd(:,:)
  double precision, allocatable :: atm_frc(:,:)
  double precision, allocatable :: atm_vel(:,:)
  double precision, allocatable :: atm_last_vel(:,:)
  double precision temp(3)  
  integer temp_wrap(3)
  double precision, intent(in) ::  pbc_box(3)
  integer,save :: numDouEle=9  ! 3 coord, one mass one qterm. Change these as needed for allocation
  proc_num_atms = 0

! count atoms

!  print *, "processor: numtasks",numtasks
  do i=1, natoms
    temp(:) = atm_crd(:,i)
    do while(temp(1) .gt. pbc_box(1)) ! pbc_box is the entire box
      temp(1)=temp(1)-pbc_box(1)      
    end do
    do while(temp(1) .lt. 0.d0)
      temp(1)=temp(1)+pbc_box(1)      
    end do
    do while(temp(2) .gt. pbc_box(2))
      temp(2)=temp(2)-pbc_box(2)      
    end do
    do while(temp(2) .lt. 0.d0)
      temp(2)=temp(2)+pbc_box(2)      
    end do
    do while(temp(3) .gt. pbc_box(3))
      temp(3)=temp(3)-pbc_box(3)      
    end do
    do while(temp(3) .lt. 0.d0)
      temp(3)=temp(3)+pbc_box(3)      
    end do
    if(temp(1) .lt. proc_max_x_crd .and. temp(1) .ge. proc_min_x_crd .and. &
       temp(2) .lt. proc_max_y_crd .and. temp(2) .ge. proc_min_y_crd .and. &
       temp(3) .lt. proc_max_z_crd .and. temp(3) .ge. proc_min_z_crd) then
      proc_num_atms = proc_num_atms + 1 !this atom is in the block, this MPI
                                        !rank own
    end if
  end do
  if(proc_num_atms .lt. proc_num_atms_min_size) then  ! temporary solution for empty sub domaind
    proc_num_atms_min_bound = proc_num_atms_min_size
!    write(0, "(I7' Has 'I7' local atoms. Setting allocation assuming 'I7' atoms')") mytaskid, proc_num_atms, proc_num_atms_min_bound
  else
    proc_num_atms_min_bound = proc_num_atms 
  endif

  call ensure_alloc_dble3(proc_atm_crd, proc_num_atms_min_bound)
  call ensure_alloc_int2d(send_buf_indices, proc_num_atms_min_bound, 6)
  call ensure_arrays_with_ghost(proc_num_atms_min_bound)

 if(.not. allocated (proc_atm_space)) allocate(proc_atm_space(natom))
 if(.not. allocated (proc_shake_space)) allocate(proc_shake_space(natom))
 if(.not. allocated (proc_atm_space_ghosts)) allocate(proc_atm_space_ghosts(natom))
 if(.not. allocated(atms_to_del)) allocate(atms_to_del(proc_num_atms_min_bound))
 if(.not. allocated(dest_sdmns)) allocate(dest_sdmns(proc_num_atms_min_bound))
 if(.not. allocated(b_send_map)) allocate(b_send_map(2, 0:numtasks-1))
 if(.not. allocated(a_send_map)) allocate(a_send_map(2, 0:numtasks-1))
 if(.not. allocated(d_send_map)) allocate(d_send_map(2, 0:numtasks-1))
 if(.not. allocated(box_len)) allocate(box_len(3))
 if(.not. allocated(atm_j_excl_list))allocate(atm_j_excl_list(1024*1024))

!output MPI transfer related
  if(.not. master) then
    if(.not. allocated(send_output_int)) allocate(send_output_int(numDouEle*proc_num_atms_min_bound*2+2,1))
    if(.not. allocated(send_output_dbl)) allocate(send_output_dbl(numDouEle*proc_num_atms_min_bound*2+2,1))
  end if

  start_index = 0 ! used for OpenMP in nb list and direct Force
  excl_img_flags = 0
  atm_j_excl_list = 0

#ifdef pmemd_SPDP
 neg2InvSqrtPi = -2.0 / sqrt(PI) 
#else
 neg2InvSqrtPi = -2.d0 / dsqrt(PI) 
#endif /* pmemd_SPDP */
!call mpi_barrier(pmemd_comm, err_code_mpi); if(mytaskid .eq. 0)  print*, mytaskid, "### Here 1120"
#ifdef _OPENMP_
!$OMP PARALLEL
        nthreads = omp_get_num_threads()
!$OMP END PARALLEL
 if(.not. allocated(frc_thread))allocate(frc_thread(3,proc_num_atms_min_bound*frc_mult*nthreads))

!call mpi_barrier(pmemd_comm, err_code_mpi); if(mytaskid .eq. 0)  print*, mytaskid, "### Here 1130"
#endif /* _OPENMP_ */
#if defined(pmemd_SPDP)
  if(.not. allocated(proc_gbl_cn1_sp)) allocate(proc_gbl_cn1_sp(1:nttyp)) 
    allocate(proc_gbl_cn2_sp(1:nttyp)) 
    proc_gbl_cn1_sp(1:nttyp) = gbl_cn1(1:nttyp)
    proc_gbl_cn2_sp(1:nttyp) = gbl_cn2(1:nttyp)
#endif /* pmemd_SPDP */
!call mpi_barrier(pmemd_comm, err_code_mpi); if(mytaskid .eq. 0)  print*, mytaskid, "### Here 1140"
  proc_atm_wrap(3,:)=0
  mult_vech_bond(:,:)=0 !used in the bond code to unwrap atom
  mult_veca_bond(:,:)=0 !used in the bond code to unwrap atom
  mult_vech_angle(:,:)=0 !used in the bond code to unwrap atom
  mult_veca_angle(:,:)=0 !used in the bond code to unwrap atom
  mult_vech_dihed(:,:)=0 !used in the bond code to unwrap atom
  mult_veca_dihed(:,:)=0 !used in the bond code to unwrap atom
  index_14(:,:) = 0 ! used in extra_pnts....F90 and i
  mult_vec_adjust(:,:) = 0 ! used in extra_pnts....F90 and i
  box_len(1:3) = pbc_box(1:3)

  j=0
  do i=1, natoms
    temp_wrap(:)=0
    temp(:) = atm_crd(:,i)
    do while(temp(1) .ge. pbc_box(1))
      temp(1)=temp(1)-pbc_box(1)
      temp_wrap(1) = temp_wrap(1) + 1      
    end do
    do while(temp(1) .lt. 0.d0)
      temp(1)=temp(1)+pbc_box(1)      
      temp_wrap(1) = temp_wrap(1) - 1
    end do
    do while(temp(2) .ge. pbc_box(2))
      temp(2)=temp(2)-pbc_box(2)      
      temp_wrap(2) = temp_wrap(2) + 1      
    end do
    do while(temp(2) .lt. 0.d0)
      temp(2)=temp(2)+pbc_box(2)      
      temp_wrap(2) = temp_wrap(2) - 1      
    end do
    do while(temp(3) .ge. pbc_box(3))
      temp(3)=temp(3)-pbc_box(3)      
      temp_wrap(3) = temp_wrap(3) + 1      
    end do
    do while(temp(3) .lt. 0.d0)
      temp(3)=temp(3)+pbc_box(3)      
      temp_wrap(3) = temp_wrap(3) - 1      
    end do
    if(temp(1) .lt. proc_max_x_crd .and. temp(1) .ge. proc_min_x_crd .and. &
       temp(2) .lt. proc_max_y_crd .and. temp(2) .ge. proc_min_y_crd .and. &
       temp(3) .lt. proc_max_z_crd .and. temp(3) .ge. proc_min_z_crd) then
      j=j+1
      proc_atm_to_full_list(j)=i
      proc_atm_space(j) = i
      proc_atm_wrap(:,j)=temp_wrap(:)
      proc_atm_crd(:,j)=temp(:) 
#if defined(pmemd_SPDP)
      proc_crd_q_sp(1:3,j)=real(temp(:)) 
      proc_crd_q_sp(4,j)=real(atm_qterm(i)) 
#endif
      proc_atm_frc(:,j)=atm_frc(:,i)
      proc_atm_vel(:,j)=atm_vel(:,i)
      proc_atm_last_vel(:,j)=atm_last_vel(:,i)
      proc_iac(j) = atm_iac(i)
      if(allocated(atm_mass)) proc_atm_mass(j)=atm_mass(i)
      if(allocated(atm_qterm)) proc_atm_qterm(j)=atm_qterm(i)
      if(allocated(atm_igraph)) proc_atm_igraph(j)=atm_igraph(i)
      if(allocated(atm_isymbl)) proc_atm_isymbl(j)=atm_isymbl(i)
    end if
  end do

  deallocate(atm_frc, atm_last_vel)
  allocate(atm_frc(3,0), atm_last_vel(3,0))
  if(allocated(atm_mass)) deallocate(atm_mass)
  if(allocated(atm_igraph)) deallocate(atm_igraph)
  if(allocated(atm_isymbl)) deallocate(atm_isymbl)
  allocate(atm_mass(0))
  allocate(atm_igraph(0))
  allocate(atm_isymbl(0))


!call mpi_barrier(pmemd_comm, err_code_mpi); if(mytaskid .eq. 0)  print*, mytaskid, "### Here 1150"
!Everything: atm_qterm, atm_mass, atm_iac, atm_igraph

!GB: atm_gb_fs, atm_gb_radii, atm_atomicnumber, atm_isymbl
  if(igb .gt. 0 .or. icnstph .eq. 2) then
    allocate(proc_atm_atomicnumber(proc_num_atms_min_bound))
    allocate(proc_atm_gb_radii(proc_num_atms_min_bound))
    allocate(proc_atm_gb_fs(proc_num_atms_min_bound))
    j=1
    do i=0,natoms
      if(proc_atm_to_full_list(j) .eq. i) then
        proc_atm_atomicnumber(j)=atm_atomicnumber(i)
        proc_atm_gb_radii(j)=atm_gb_radii(i)
        proc_atm_gb_fs(j)=atm_gb_fs(i)
        j=j+1
      end if
    end do
  end if
!  print *, proc_num_atms_min_bound, myblockid
!call mpi_barrier(pmemd_comm, err_code_mpi); if(mytaskid .eq. 0)  print*, mytaskid, "### Here 1160"

end subroutine alloc_full_list_to_proc

!*******************************************************************************
!
! Subroutine:  resetup_processor
!
! Description: Sets up processor information
!
!*******************************************************************************

subroutine resetup_processor(pbc_box, bkt_size)

  implicit none
  double precision, intent(in)  :: pbc_box(3), bkt_size(3)

  call assign_cit_buckets(pbc_box)

! call processor_ewald_setup(pbc_box)

  ! We find the neighborbuckets for entire space (including ghost space) that we
  ! on
  call cornerbuckets_for_range(2,int_bkt_minx-1,int_bkt_miny-1,int_bkt_minz-1,int_bkt_maxx+1,int_bkt_maxy+1,int_bkt_maxz+1)

! write(15+mytaskid,*) "dimensions",proc_bkt_minx, proc_bkt_miny, &
!proc_bkt_minz, proc_bkt_maxx, proc_bkt_maxy, proc_bkt_maxz

end subroutine resetup_processor

!*******************************************************************************
!
! Subroutine:  setup_processor
!
! Description: Sets up processor information
!              
!*******************************************************************************

subroutine setup_processor(pbc_box, bkt_size)

  implicit none
  double precision, intent(in)  :: pbc_box(3), bkt_size(3)
  integer p(3), np(3,6), i

  new_bkt_size = bkt_size

  ! Get the MPI topology
  call setup_processor_decomposition(pbc_box, proc_dimx, proc_dimy, proc_dimz, numtasks)
  ! Get the MPI ranks of the neighbors 6 faces of the domain
  np=0
  proc_neighbor_rank=0
  call blockunflatten(mytaskid,p(1),p(2),p(3))
  np(1,:) = p(1)
  np(2,:) = p(2)
  np(3,:) = p(3)
  np(1,1) = mod(proc_dimx+p(1)-1,proc_dimx)
  np(1,2) = mod(proc_dimx+p(1)+1,proc_dimx)
  np(2,3) = mod(proc_dimy+p(2)-1,proc_dimy)
  np(2,4) = mod(proc_dimy+p(2)+1,proc_dimy)
  np(3,5) = mod(proc_dimz+p(3)-1,proc_dimz)
  np(3,6) = mod(proc_dimz+p(3)+1,proc_dimz)
  do i =1,6
    proc_neighbor_rank(i) = np(1,i) + ( np(2,i) + np(3,i)*proc_dimy )*proc_dimx
  end do
  ! Setup the offset required for the sides of the MPI tasks at the edge of the domain
  proc_side_wrap_offset = 0.0
  if(p(1) .eq. 0) proc_side_wrap_offset(1,1) = -1.0 !box_len(1)
  if(p(2) .eq. 0) proc_side_wrap_offset(2,3) = -1.0 !box_len(2)
  if(p(3) .eq. 0) proc_side_wrap_offset(3,5) = -1.0 !box_len(3)
  if(p(1) .eq. proc_dimx-1) proc_side_wrap_offset(1,2) = 1.0 !-1*box_len(1)
  if(p(2) .eq. proc_dimy-1) proc_side_wrap_offset(2,4) = 1.0 !-1*box_len(2)
  if(p(3) .eq. proc_dimz-1) proc_side_wrap_offset(3,6) = 1.0 !-1*box_len(3)


  call setup_processor_decomposition(pbc_box, proc_dimx_recip, proc_dimy_recip, proc_dimz_recip, numtasks_recip)

  proc_max_num_ghost_atm = max(proc_dimx, proc_dimy, proc_dimz)
  proc_max_num_ghost_atm = 100*proc_max_num_ghost_atm*proc_max_num_ghost_atm !assuming max atoms per bucket = 100

  ! These should be called if CIT changes
  call assign_cit_buckets(pbc_box)
 
  call processor_ewald_setup(pbc_box)
 
  ! We find the neighborbuckets for entire space (including ghost space) that we
  ! on
  call cornerbuckets_for_range(2,int_bkt_minx-1,int_bkt_miny-1,int_bkt_minz-1,int_bkt_maxx+1,int_bkt_maxy+1,int_bkt_maxz+1)

  call setup_neighbor_mpi_map

  frc_mult=50
 
end subroutine setup_processor

!*******************************************************************************
!
! Subroutine:  setup_processor_decomposition
!
! Description: Compute the processers coordinates bounds
!              
!*******************************************************************************


subroutine setup_processor_decomposition(pbc_box, dimx, dimy, dimz, tasks)

!  use file_io_dat_mod
!  use parallel_dat_mod
!  use pmemd_lib_mod,     only : mexit

  implicit none

  double precision, intent(in)  :: pbc_box(3)
  integer                       :: dimx, dimy, dimz, tasks

  ! local variables
  integer facX,facY,facZ
  double precision xySurface,xzSurface, yzSurface,netSurface,newSurface
  integer ModulosX,ModulosY;
! we are using # of buckets in 3D, not the box length
! later if the following breaks, then we need to use box length 
!currently it is assumed that always box length is divisible by cut-off,
!otherwise the code will break????


  xySurface= pbc_box(1)*pbc_box(2)
  xzSurface= pbc_box(1)*pbc_box(3)
  yzSurface= pbc_box(2)*pbc_box(3)
  !Surface Area = SA
  netSurface= 2 * (xySurface+ xzSurface + yzSurface)  ! total surface

  !Find the first factor of numTasks
  !facx is the number of division of doman in x direction
  do facX=1, tasks,1 
    ModulosX= mod(tasks,facX)
    if(ModulosX .eq. 0) then
      !Find the second factor of numTasks(which is the first factor of
      !numTasks/facX
      !facy is the number of division of domain in y direction
      do facY=1, tasks/facX,1
        ModulosY = mod(tasks/facX,facY)
        !Find the last factor of numTasks, safe to divide because modulo always is 0
        !so it divides evenly
        if(ModulosY .eq. 0) then
          !facz is the number of division of domanin in z direction
          facZ=((tasks/facX)/facY)
          !If the SA we get is the smallest, then grid is being
          !divided most evenly(larger block sizes)
          newSurface = (xySurface/facX)/facY+(xzSurface/facX)/facZ+(yzSurface/facZ)/facY
          if(newSurface .lt. netSurface) then
            netSurface= newSurface
            dimx=facX
            dimy=facY
            dimz=facZ
          end if
        end if
      end do
    end if
  end do


!  if(pbc_box(1)/new_bkt_size(1) .lt. 1 .or. pbc_box(2)/new_bkt_size(2) .lt. 1 .or. pbc_box(3)/new_bkt_size(3) .lt. 1) then
!     if(master) then
!       write (mdout,'(a,3I8,a)') 'THERE ARE MORE PROCESSORS IN A DIMENSION THAN THERE &
!                               ARE BUCKETS.  TRY A DIFFERENT MPI COUNT.  THE &
!                               MAX MPI COUNT SUPPORTED BY YOUR SYSTEM IN EACH DIMENSION: ', &
!                               int(pbc_box(1)/new_bkt_size(1)), int(pbc_box(2)/new_bkt_size(2)), &
!                               int(pbc_box(3)/new_bkt_size(3)),'. Try something with more factors.'
!       call mexit(mdout, 1)
!     end if
!  end if
end subroutine setup_processor_decomposition


!*******************************************************************************
!
! Subroutine:  assign_cit_buckets
!
! Description: Assign cit ID of processors. Sets proc_minx, proc_maxx,
! proc_miny, proc_maxy, proc_minz, proc_maxz, proc_ownedbkts
!
!*******************************************************************************

!Need to test the following routine
subroutine assign_cit_buckets(pbc_box)

  use parallel_dat_mod
!  use cit_mod, only: bkt_size

  implicit none

  double precision, intent(in)  :: pbc_box(3)
  real :: real_proc_bkt_minx,real_proc_bkt_maxx,real_proc_bkt_miny,real_proc_bkt_maxy,real_proc_bkt_minz,real_proc_bkt_maxz
  integer flattened
  integer x,y,z,i,j,k,mpiid
  integer tmpx,tmpy,tmpz

  flattened=mytaskid
  myblockid=mytaskid

  
  !Decompose owned block into x,y,z indicies
  call blockunflatten(flattened,x,y,z)

  int_bkt_minx=1
  int_bkt_miny=1
  int_bkt_minz=1
 

  if(allocated(crd_boundx_recip)) deallocate(crd_boundx_recip, crd_boundy_recip, crd_boundz_recip)
  if(.not. allocated(crd_boundx_recip)) allocate(crd_boundx_recip(2,proc_dimx))
  if(.not. allocated(crd_boundy_recip)) allocate(crd_boundy_recip(2,proc_dimy))
  if(.not. allocated(crd_boundz_recip)) allocate(crd_boundz_recip(2,proc_dimz))

  do i = 1, proc_dimx
    crd_boundx_recip(1,i) = pbc_box(1)/proc_dimx * (i-1)
    crd_boundx_recip(2,i) = pbc_box(1)/proc_dimx * i
  end do
  do i = 1, proc_dimy
    crd_boundy_recip(1,i) = pbc_box(2)/proc_dimy * (i-1)
    crd_boundy_recip(2,i) = pbc_box(2)/proc_dimy * i
  end do
  do i = 1, proc_dimz
    crd_boundz_recip(1,i) = pbc_box(3)/proc_dimz * (i-1)
    crd_boundz_recip(2,i) = pbc_box(3)/proc_dimz * i
  end do

  proc_crd_bound(1,1) = pbc_box(1)/proc_dimx * x
  proc_crd_bound(2,1) = pbc_box(1)/proc_dimx * (x+1)
  proc_crd_bound(1,2) = pbc_box(2)/proc_dimy * y
  proc_crd_bound(2,2) = pbc_box(2)/proc_dimy * (y+1)
  proc_crd_bound(1,3) = pbc_box(3)/proc_dimz * z
  proc_crd_bound(2,3) = pbc_box(3)/proc_dimz * (z+1)

  bkt_dim(1:3) = ceiling((proc_crd_bound(2,1:3)-proc_crd_bound(1,1:3))/new_bkt_size(1:3))


  proc_min_x_crd = proc_crd_bound(1,1)
  proc_max_x_crd = proc_crd_bound(2,1)
  proc_min_y_crd = proc_crd_bound(1,2)
  proc_max_y_crd = proc_crd_bound(2,2)
  proc_min_z_crd = proc_crd_bound(1,3)
  proc_max_z_crd = proc_crd_bound(2,3)

  int_bkt_dimx=bkt_dim(1)
  int_bkt_dimy=bkt_dim(2)
  int_bkt_dimz=bkt_dim(3)

  int_bkt_maxx=bkt_dim(1)
  int_bkt_maxy=bkt_dim(2)
  int_bkt_maxz=bkt_dim(3)

  proc_crd_len(1) = proc_max_x_crd-proc_min_x_crd
  proc_crd_len(2) = proc_max_y_crd-proc_min_y_crd
  proc_crd_len(3) = proc_max_z_crd-proc_min_z_crd

  ghost_crd_len = new_bkt_size*2.d0+proc_crd_len

  proc_numbkts=(int_bkt_dimx+2)*(int_bkt_dimy+2)*(int_bkt_dimz+2)
end subroutine assign_cit_buckets

!*******************************************************************************
!
! Subroutine:  processor_ewald_setup
!
! Description: Processor setup for ewald.  Sets up communicators and the block
! dimensions for ewald
!
!*******************************************************************************

subroutine processor_ewald_setup(pbc_box)

  use mdin_ewald_dat_mod
  use parallel_dat_mod

  implicit none
  
  double precision, intent(in)  :: pbc_box(3)
  integer flattened
  integer x,y,z,i
  double precision blockmincrd, blockmaxcrd
  double precision nfft1len, nfft2len, nfft3len, minlenx, minleny, minlenz
  double precision blockcrdstart, blockcrdend, pmecrdstart, pmecrdend
  integer pmelenx, pmeleny, pmelenz, nfftstart, nfftend
  integer proc_other_minx, proc_other_maxx !for other MPI ranks
  integer proc_other_miny, proc_other_maxy  !for other MPI ranks
  integer proc_other_minz, proc_other_maxz  !for other MPI ranks
  integer proc_rangex,proc_rangey,proc_rangez !for other MPI ranks


  !getting 3d MPI co-ordinates, my_mpix, my_mpiy, my_mpiz from myblockid
  call blockunflatten_recip(myblockid,my_mpix,my_mpiy,my_mpiz) 

! Assigns min, max for block
  !allocate the buffers needed for MPI_Alltoallv
  if(allocated(send_xfft_cnts)) deallocate(send_xfft_cnts  ,send_xfft_offset &
                   ,recv_xfft_cnts  ,recv_xfft_offset  ,send_yfft_cnts &
                   ,send_yfft_offset  ,recv_yfft_cnts  ,recv_yfft_offset &
                   ,send_zfft_cnts  ,send_zfft_offset  ,recv_zfft_cnts &
                   ,recv_zfft_offset  ,x_offset  ,y_offset  ,z_offset)


  if(.not. allocated(send_xfft_cnts)) allocate(send_xfft_cnts(0:proc_dimx_recip-1))
  if(.not. allocated(send_xfft_offset)) allocate(send_xfft_offset(0:proc_dimx_recip-1))
  if(.not. allocated(recv_xfft_cnts)) allocate(recv_xfft_cnts(0:proc_dimx_recip-1))
  if(.not. allocated(recv_xfft_offset)) allocate(recv_xfft_offset(0:proc_dimx_recip-1))
  if(.not. allocated(send_yfft_cnts)) allocate(send_yfft_cnts(0:proc_dimy_recip-1))
  if(.not. allocated(send_yfft_offset)) allocate(send_yfft_offset(0:proc_dimy_recip-1))
  if(.not. allocated(recv_yfft_cnts)) allocate(recv_yfft_cnts(0:proc_dimy_recip-1))
  if(.not. allocated(recv_yfft_offset)) allocate(recv_yfft_offset(0:proc_dimy_recip-1))
  if(.not. allocated(send_zfft_cnts)) allocate(send_zfft_cnts(0:proc_dimz_recip-1))
  if(.not. allocated(send_zfft_offset)) allocate(send_zfft_offset(0:proc_dimz_recip-1))
  if(.not. allocated(recv_zfft_cnts)) allocate(recv_zfft_cnts(0:proc_dimz_recip-1))
  if(.not. allocated(recv_zfft_offset)) allocate(recv_zfft_offset(0:proc_dimz_recip-1))
  if(.not. allocated(x_offset)) allocate(x_offset(0:proc_dimx_recip-1))
  if(.not. allocated(y_offset)) allocate(y_offset(0:proc_dimy_recip-1))
  if(.not. allocated(z_offset)) allocate(z_offset(0:proc_dimz_recip-1))

  ! Grabs length of one nfft grid
  nfft1len=pbc_box(1)/dble(nfft1)
  nfft2len=pbc_box(2)/dble(nfft2)
  nfft3len=pbc_box(3)/dble(nfft3)
  ! Grabs maximum number of nfft grids per block
  minlenx = dble(nfft1)/dble(proc_dimx_recip)
  minleny = dble(nfft2)/dble(proc_dimy_recip)
  minlenz = dble(nfft3)/dble(proc_dimz_recip)
  pmelenx = ceiling(minlenx)*2 ! max number of grid in x for 1 block
  pmeleny = ceiling(minleny)*2 ! max number of grid in y for 1 block
  pmelenz = ceiling(minlenz)*2 ! max number of grid in z for 1 block
  ! Length of an individual bucket
  ! Grabs my blockids
  call blockunflatten_recip(myblockid,x,y,z)! find x,y,z of my block

  ! Does a check in x if the pmelenx works
  pmecrdstart=0.d0 ! crd  of block 
  nfftstart=1 ! integer crd of fft grid
  do i=0,proc_dimx_recip-1 ! run through all blocks in x direction
!    pmelenx = ceiling (size of the block)
    pmecrdend=pmecrdstart+nfft1len*pmelenx
    nfftend=nfftstart+pmelenx-1 ! To make sure the start and end wroks

    blockmincrd = crd_boundx_recip(1,i+1)
    blockmaxcrd = crd_boundx_recip(2,i+1)

    blockcrdstart=blockmincrd ! Block start without ghost
    blockcrdend=blockmaxcrd ! Block end without ghost
    do while(blockcrdend .lt. pmecrdend) !1st check. blockcrdend includes the ghost space, if crossed the block, reduce the nfftend
      pmecrdend=pmecrdend-nfft1len
      nfftend=nfftend-1 ! make sure we are inside
    end do
    ! Shifts it to be one out of block (overestimate)
    pmecrdend=pmecrdend+nfft1len
    nfftend=nfftend+1
    !2nd check is needed since the above blockcrdend includes ghost space 
    if(nfftend .gt. nfft1) then ! 2nd check if we are insidei nfft1, if not we reduce the nfftend
      do while(nfftend .ne. nfft1)
        pmecrdend=pmecrdend-nfft1len
        nfftend=nfftend-1 !make sure, we are inside
      end do
    end if
    if(i .eq. x) then !if this is myblock
      proc_pme_minx = pmecrdstart
      proc_pme_maxx = pmecrdend
      proc_min_nfft1 = nfftstart ! this is the one we use for fft range
      proc_max_nfft1 = nfftend! this is the one we use for fft range
      !send_xfft_cnts(i) = 0 ! send 0 to myself 
      !keep my own fft count also, will be easy later in fft code
      send_xfft_cnts(i) = nfftend - nfftstart + 1 ! remvoe the send_ word 
      my_xfft_cnt = nfftend - nfftstart + 1
      x_offset(i) = nfftstart-1
    else
     send_xfft_cnts(i) = nfftend - nfftstart + 1 
     x_offset(i) = nfftstart-1
    end if
    pmecrdstart=pmecrdend !ready for the next block
    nfftstart=nfftend+1 !ready for the next block
  end do
  ! Does a check in y if pmeleny works

!the explanation of the following y and z are same as the above x
!thus we dont have much explanataiton in the below
  pmecrdstart=0.d0
  nfftstart=1
  do i=0,proc_dimy_recip-1
    pmecrdend=pmecrdstart+pmeleny*nfft2len
    nfftend=nfftstart+pmeleny-1

    blockmincrd = crd_boundy_recip(1,i+1)
    blockmaxcrd = crd_boundy_recip(2,i+1)

    blockcrdstart=blockmincrd
    blockcrdend=blockmaxcrd
    do while(blockcrdend .lt. pmecrdend)
      pmecrdend=pmecrdend-nfft2len
      nfftend=nfftend-1
    end do
    pmecrdend=pmecrdend+nfft2len
    nfftend=nfftend+1
    if(nfftend .gt. nfft2) then
      do while(nfftend .ne. nfft2)
        pmecrdend=pmecrdend-nfft2len
        nfftend=nfftend-1
      end do
    end if
    if(i .eq. y) then
      proc_pme_miny=pmecrdstart
      proc_pme_maxy=pmecrdend
      proc_min_nfft2=nfftstart
      proc_max_nfft2=nfftend
      !send_yfft_cnts(i) = 0 ! send 0 to myself 
      my_yfft_cnt = nfftend - nfftstart + 1
      !keep my own fft count also, will be easy later in fft code
      send_yfft_cnts(i) = nfftend - nfftstart + 1 !change the name to only
                                                  !yfft_cnts make more sense 
      y_offset(i) = nfftstart-1
     !write(15+mytaskid_recip,*)"Selected"
    else
     send_yfft_cnts(i) = nfftend - nfftstart + 1 
     y_offset(i) = nfftstart-1
    end if
    pmecrdstart=pmecrdend !ready for the next block
    nfftstart=nfftend+1 !ready for the next block
  end do
  pmecrdstart=0.d0
  nfftstart=1
  do i=0,proc_dimz_recip-1
    pmecrdend=pmecrdstart+pmelenz*nfft3len
    nfftend=nfftstart+pmelenz-1

    blockmincrd = crd_boundz_recip(1,i+1)
    blockmaxcrd = crd_boundz_recip(2,i+1)

    blockcrdstart=blockmincrd
    blockcrdend=blockmaxcrd
    do while(blockcrdend .lt. pmecrdend)
      pmecrdend=pmecrdend-nfft3len
      nfftend=nfftend-1
    end do
    pmecrdend=pmecrdend+nfft3len
    nfftend=nfftend+1
    if(nfftend .gt. nfft3) then
      do while(nfftend .ne. nfft3)
        pmecrdend=pmecrdend-nfft3len
        nfftend=nfftend-1
      end do
    end if
    if(i .eq. z) then
      proc_pme_minz=pmecrdstart
      proc_pme_maxz=pmecrdend
      proc_min_nfft3=nfftstart
      proc_max_nfft3=nfftend
      !send_zfft_cnts(i) = 0 ! send 0 to myself 
      !keep my own fft count also, will be easy later in fft code
      send_zfft_cnts(i) = nfftend - nfftstart + 1  ! remvoe the send_ word 
      my_zfft_cnt = nfftend - nfftstart + 1
      z_offset(i) = nfftstart-1
     !write(15+mytaskid_recip,*)"Selected"
    else
     send_zfft_cnts(i) = nfftend - nfftstart + 1 
     z_offset(i) = nfftstart-1
    end if

    pmecrdstart=pmecrdend
    nfftstart=nfftend+1
  end do


!find the max fft_x, fft_y, fft_z
!Use this value later to allocate the MPI buffer in pme_recip_midpoint file
  max_xfft_cnts = max(maxval(send_xfft_cnts(0:proc_dimx_recip-1)),(proc_max_nfft1-proc_min_nfft1+1))
  max_yfft_cnts = max(maxval(send_yfft_cnts(0:proc_dimy_recip-1)),(proc_max_nfft2-proc_min_nfft2+1))
  max_zfft_cnts = max(maxval(send_zfft_cnts(0:proc_dimz_recip-1)),(proc_max_nfft3-proc_min_nfft3+1))
     
! Creates a list of all blocks in x, xy, and z directions not including self.

  if(allocated(proc_x_comm_lst)) deallocate(proc_x_comm_lst)
  if(allocated(proc_y_comm_lst)) deallocate(proc_y_comm_lst)
  if(allocated(proc_xy_comm_lst)) deallocate(proc_xy_comm_lst)
  if(allocated(proc_z_comm_lst)) deallocate(proc_z_comm_lst)

  allocate(proc_x_comm_lst(0:proc_dimx_recip-1))
  allocate(proc_y_comm_lst(0:proc_dimy_recip-1))
  allocate(proc_xy_comm_lst(0:proc_dimx_recip*proc_dimy_recip-1))
  allocate(proc_z_comm_lst(0:proc_dimz_recip-1))

  call blockunflatten_recip(myblockid,x,y,z) ! same y and z will be used
  i=1
  do x=0,proc_dimx_recip-1
    flattened=x+proc_dimx_recip*y+z*proc_dimx_recip*proc_dimy_recip 
  !  if(flattened .ne. myblockid) then
      !proc_x_comm_lst(i)=flattened ! y and z are same for my peers in x directions
      proc_x_comm_lst(x)=flattened ! y and z are same for my peers in x directions
            !we include even myblockid. Later when use it, we can excluded
            !sending myself
  !    i=i+1
  !  end if
  end do

  call blockunflatten_recip(myblockid,x,y,z) ! same x and z will be used
  i=1
  do y=0,proc_dimy_recip-1
    flattened=x+proc_dimx_recip*y+z*proc_dimx_recip*proc_dimy_recip
    !if(flattened .ne. myblockid) then
      !proc_y_comm_lst(i)=flattened ! x and z are same for my peers in y directions
      proc_y_comm_lst(y)=flattened ! x and z are same for my peers in y directions
            !we include even myblockid. Later when use it, we can excluded
            !sending myself
    !  i=i+1
    !end if
  end do
  call blockunflatten_recip(myblockid,x,y,z) !same z will be used, used in 
                          !fft comm in x and y direction
  i=1
  !do x=0,proc_dimx_recip-1
  !  do y=0,proc_dimy_recip-1
  do y=0,proc_dimy_recip-1
    do x=0,proc_dimx_recip-1
      flattened=x+proc_dimx_recip*y+z*proc_dimx_recip*proc_dimy_recip ! same z
      !if(flattened .ne. myblockid) then
        !proc_xy_comm_lst(i)=flattened
        proc_xy_comm_lst(x+y*proc_dimx_recip)=flattened
      !  i=i+1
      !end if
    end do
  end do

  call blockunflatten_recip(myblockid,x,y,z)
  i=1
  do z=0,proc_dimz_recip-1
    flattened=x+proc_dimx_recip*y+z*proc_dimx_recip*proc_dimy_recip
    !if(flattened .ne. myblockid) then
      !proc_z_comm_lst(i)=flattened ! x and y are same for my peers in z directions
      proc_z_comm_lst(z)=flattened ! x and y are same for my peers in z directions
            !we include even myblockid. Later when use it, we can excluded
            !sending myself
   !   i=i+1
   ! end if
  end do 

end subroutine processor_ewald_setup

!*******************************************************************************
! function: int_unflatten
! Description: takes a 1d bucketid and converts it to a 3d coordinate
! 
! Parameters: mybucketid is the id to be coordinates to be flattened.
! 
! Can call for example by 
! !print *,"x,y,z is ",x,y,z," Flattened is ",flatten(x,y,z)
!*******************************************************************************
subroutine int_unflatten(mybucketid,x,y,z)
  integer x,y,z
  integer flattened,mybucketid,minx,maxx,miny,maxy,minz,maxz
  integer proc_bkt_minx,proc_bkt_miny,proc_bkt_minz,proc_bkt_maxx,proc_bkt_maxy,proc_bkt_maxz


  flattened=mybucketid-1
  z=flattened/((int_bkt_dimx+2)*(int_bkt_dimy+2))
  flattened=flattened - (z*(int_bkt_dimx+2)*(int_bkt_dimy+2))
  y=flattened/(int_bkt_dimx+2)
  x=mod(flattened  , (int_bkt_dimx+2))

end subroutine int_unflatten


!*******************************************************************************
! function: int_flatten
! Description: takes a 3d coordinate and flattens it into a single value
! 
! Parameters: x,y,z are the coordinates to be flattened.
! 
! Can call for example by 
! !print *,"x,y,z is ",x,y,z," Flattened is ",flatten(x,y,z)
!*******************************************************************************

integer function int_flatten(x,y,z)

  integer x,y,z

  int_flatten=x+y*(int_bkt_dimx+2)+z*(int_bkt_dimx+2)*(int_bkt_dimy+2)+1

  !flatten=x+(y*sys_bkt_dimx*proc_dimx)+(z*sys_bkt_dimx*proc_dimx*sys_bkt_dimy*proc_dimy)

end function int_flatten

!*******************************************************************************
!
! Subroutine:  simple_i_neighbors
!
! Description: Given a distance i and a bucketid
!
! Variables:
!   mybucketid: the bucketid that from which the bottom left neighborbuckets are taken from
!   index: is single_bucket_neighborlist's array index
!   i: the distance of the neighborbuckets being taken from
!
!*******************************************************************************

subroutine simple_i_neighbors (mybucketid,i,index,single_bucket_neighborlist)
  implicit none
  integer :: i,index
  integer :: x,y,z,xid,yid,zid,flattened,mybucketid
  integer, intent(INOUT):: single_bucket_neighborlist(:)
  integer  :: bucket

  call int_unflatten(mybucketid,xid,yid,zid)

  !Transfer value of z,y,x to zid,yid,xid, because z,y,x will be used as iteration variables

  do z = zid -i , zid !zid, zid+i 
    do y = yid -i , yid + i
      do x = xid -i , xid + i
        ! If we're already a ghost bucket we don't want to grab imaginary
        ! buckets.
        if((z .lt. zid) .or. (y .lt. yid) .or. ((y .eq. yid) .and. (x .lt. xid))) then !lower left

          if(x .lt. 0 .or. x .gt. int_bkt_dimx+1 .or. &
             y .lt. 0 .or. y .gt. int_bkt_dimy+1 .or. &
             z .lt. 0 .or. z .gt. int_bkt_dimz+1) cycle
          bucket = int_flatten(x, y, z)
          single_bucket_neighborlist(index) = bucket
          index = index + 1
        end if
      end do
    end do
  end do

end subroutine simple_i_neighbors

!*******************************************************************************
!
! Subroutine:  cornerbuckets_for_range
!
! Description: Given a distance i, and a block defined by min x,y,z and max x,y,z,
! return all neighborbuckets and itself within i in the range
! Give all the corner neighborbuckets into the 2d array proc_corner_neighborbuckets
!
! Note: Take care to insure that the minimums are indeed smaller than the maximums
!
!*******************************************************************************

subroutine cornerbuckets_for_range(i,minx,miny,minz,maxx,maxy,maxz)

  integer x,y,z,flattened,i
  integer minx,miny,minz,maxx,maxy,maxz,tmpx,tmpy,tmpz
  integer xindex,yindex,numElements,numBuckets,counter,mybucketid
  integer, allocatable :: single_bucket_neighborlist(:)
  integer :: bucket_cnt(minx:maxx,miny:maxy,minz:maxz)

  integer bkt, swap, j, cnt

  !xindex=1;
  yindex=1;
  numElements=0
  numBuckets=0

  !numBuckets is the range of indicies needed for the system since we don't
  !search upper we don't need the buffer zone for the z axis
  numBuckets=(maxx-minx+3)*(maxy-miny+3)*(maxz-minz+2)

  !numElements in the number of neighborbuckets
  !(i*((2*i+1)**2)) is the layers of the block that are taken entirely
  !(i*(2*i+1)) is the half layer taken at the center bucket layer
  ! +i is for the half column of buckets
  ! +1 is for the center bucket id itself
  numElements=(i*((2*i+1)**2))+(i*(2*i+1))+i+1
  numElements=250
! write(0,*)"Elements", numElements
  !Allocate the 2d array
  if(allocated(proc_corner_neighborbuckets))deallocate(proc_corner_neighborbuckets)
  allocate(proc_corner_neighborbuckets(numElements,numBuckets))

  proc_corner_neighborbuckets(:,:)=-1

  !Allocate the temporary 1d array
  allocate(single_bucket_neighborlist(numElements))
  !For each bucket in the range as defined by the min and max
  do z=minz,maxz
    do y=miny,maxy
      do x=minx,maxx

        !Flatten the current iteration(which is the coordinates of the bucket in the block)
        mybucketid=int_flatten(x,y,z)
        !Add the center bucket(the origin) to the 1d array
        single_bucket_neighborlist(:)=-1
        single_bucket_neighborlist(yindex)=int_flatten(x,y,z)
        yindex=yindex+1
        !Then add all the neighborbuckets in the range
        !ashraf wrote the simpler one
        call simple_i_neighbors(mybucketid,i,yindex,single_bucket_neighborlist)
        !Then set single_bucket_neighborlist at the xindex to the values in the proc_corner_neighborbuckets
        !set the xindex as the mybucketid, so that all the neighborbucket can be
        !found at the mybicketid index of the array. Thus it will be straight
        !forward to find the neighbor buckets in the neighbor list code
        xindex = mybucketid 
        do counter=1,size(single_bucket_neighborlist)
          proc_corner_neighborbuckets(counter,xindex)=single_bucket_neighborlist(counter)
        end do
        !Reset the 1d array index
        yindex=1;
        !Increment the 2d array index
        !xindex=xindex+1
      end do
    end do
  end do

  !Deallocate the 1d array
  deallocate(single_bucket_neighborlist)

  ! Compute the neighbors count for each bucket
  if(allocated(proc_corner_neighborbuckets_cnt)) &
            deallocate(proc_corner_neighborbuckets_cnt)
  allocate(proc_corner_neighborbuckets_cnt(numBuckets))
  do bkt = 1, numBuckets
    j=1
    do while ((proc_corner_neighborbuckets(j, bkt) .ne. -1) .and. (j<= size(proc_corner_neighborbuckets, 1)) )
      j=j+1
    end do
    proc_corner_neighborbuckets_cnt(bkt)=j-1
  end do
  ! put the current bucket at the end of the list for more efficient neighbor build list
  do bkt = 1, numBuckets
    cnt = proc_corner_neighborbuckets_cnt(bkt)
    if(cnt .eq. 0) cycle
    swap = proc_corner_neighborbuckets(cnt,bkt)
    proc_corner_neighborbuckets(cnt,bkt) = proc_corner_neighborbuckets(1,bkt)
    proc_corner_neighborbuckets(1,bkt) = swap
  end do


  !Following statement prints out the values in the 2d array, uncomment to debug
! bucket_cnt(:,:,:)=0
!  do y=1,numBuckets
!    do x=1,numElements
!       if(proc_corner_neighborbuckets(x,y) .ne. -1) then
!         call int_unflatten(proc_corner_neighborbuckets(x,y),tmpx,tmpy,tmpz)
!         write(15+mytaskid,*)"tmp",tmpx,tmpy,tmpz
!         write(15+mytaskid,*)x,' ',int_ext_bkt(y)+1,' ',int_ext_bkt(proc_corner_neighborbuckets(x,y)+1)
!         bucket_cnt(tmpx,tmpy,tmpz)=bucket_cnt(tmpx,tmpy,tmpz)+1
!       end if
!    end do
!  end do
!  do x=0,maxx
!    do y=0,maxy
!      do z=0,maxz
!        write(15+mytaskid,*)"Bkt cnt: ",x,y,z,bucket_cnt(x,y,z)
!      end do
!    end do
!  end do
end subroutine cornerbuckets_for_range

!*******************************************************************************
! function: blockunflatten_recip
! Description: takes a flattened 3d recip block coordinate and unflattens to 3d coordinates
! 
! Parameters: myblockid is the value to be unflattened.
! 
!*******************************************************************************
subroutine blockunflatten_recip(myblockid,x,y,z)
  integer x,y,z
  integer flattened,myblockid
  flattened=myblockid
  z=flattened/(proc_dimx_recip*proc_dimy_recip)
  flattened=flattened - (z*proc_dimx_recip*proc_dimy_recip)
  y=flattened/proc_dimx_recip
  x=mod(flattened  , proc_dimx_recip)
end subroutine blockunflatten_recip

!*******************************************************************************
! function: destinationmpi
! Description: Given atom crds get calculates the destination processor mpi task
!
! Parameters: Atm crd x, atm crd y, atm crd z returns mpi task id
!
!********************************************************************************
subroutine destinationmpi(neighborid, crd)
   double precision, intent(in)   :: crd(3)
   integer, intent(out)           :: neighborid
   integer                        :: neighbor_proc(3), proc(3)
   ! New approach through coordinates information
   call blockunflatten(mytaskid,proc(1),proc(2),proc(3))
   neighbor_proc = proc
   if(crd(1) .lt. proc_min_x_crd) neighbor_proc(1) = mod(proc_dimx+proc(1)-1,proc_dimx)
   if(crd(1) .ge. proc_max_x_crd) neighbor_proc(1) = mod(proc_dimx+proc(1)+1,proc_dimx)
   if(crd(2) .lt. proc_min_y_crd) neighbor_proc(2) = mod(proc_dimy+proc(2)-1,proc_dimy)
   if(crd(2) .ge. proc_max_y_crd) neighbor_proc(2) = mod(proc_dimy+proc(2)+1,proc_dimy)
   if(crd(3) .lt. proc_min_z_crd) neighbor_proc(3) = mod(proc_dimz+proc(3)-1,proc_dimz)
   if(crd(3) .ge. proc_max_z_crd) neighbor_proc(3) = mod(proc_dimz+proc(3)+1,proc_dimz)
   neighborid = neighbor_proc(1) + ( neighbor_proc(2)+neighbor_proc(3)*proc_dimy )*proc_dimx
end subroutine

!*******************************************************************************
! function: blockunflatten
! Description: takes a flattened 3d block coordinate and unflattens to 3d coordinates
! 
! Parameters: myblockid is the value to be unflattened.
! 
!*******************************************************************************
subroutine blockunflatten(myblockid,x,y,z)
  integer x,y,z
  integer flattened,myblockid
  flattened=myblockid
  z=flattened/(proc_dimx*proc_dimy)
  flattened=flattened - (z*proc_dimx*proc_dimy)
  y=flattened/proc_dimx
  x=mod(flattened  , proc_dimx)
end subroutine blockunflatten

!*******************************************************************************
!
! Subroutine:  setup_neighbor_mpi_map
!
! Description: sets up the list of neighbors MPI ranks
!              
!*******************************************************************************
subroutine setup_neighbor_mpi_map()
  implicit none

  integer p(3),np(3),neighbor_rank, i, j,k

  neighbor_mpi_map = -1
  neighbor_mpi_cnt = 0
  call blockunflatten(mytaskid,p(1),p(2),p(3))
  do i=-1,1
    do j=-1,1
      do k=-1,1
        if(i==0 .and. j==0. .and. k==0) cycle
        np(1) = mod(proc_dimx+p(1)+i,proc_dimx)
        np(2) = mod(proc_dimy+p(2)+j,proc_dimy)
        np(3) = mod(proc_dimz+p(3)+k,proc_dimz)
        neighbor_rank = np(1) + ( np(2) + np(3)*proc_dimy )*proc_dimx
        if(all(neighbor_mpi_map(:neighbor_mpi_cnt) .ne. neighbor_rank) ) then
          neighbor_mpi_cnt = neighbor_mpi_cnt+1
          neighbor_mpi_map(neighbor_mpi_cnt)= neighbor_rank
        end if
      end do
    end do
  end do

end subroutine setup_neighbor_mpi_map


!*******************************************************************************
!
! Subroutine:  fill_ghost_tranvec
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine fill_ghost_tranvec(tranvec)

  implicit none

  integer      :: tranvec(3, 27) ! this is (1:3,0:17) externally.
  integer               :: iv, i0, i1, i2, i3

! This works for both orthogonal and nonorthogonal unit cells.

  iv = 0

  do i3 = -1, 1
    do i2 = -1, 1
      do i1 = -1, 1
        iv = iv + 1
          tranvec(1, iv) = i1
          tranvec(2, iv) = i2
          tranvec(3, iv) = i3
      end do
    end do
  end do

  return

end subroutine fill_ghost_tranvec

#endif /*MPI*/
end module processor_mod
