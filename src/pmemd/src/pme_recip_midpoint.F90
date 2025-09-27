#include "copyright.i"

!*******************************************************************************
!
! Module: pme_recip_midpoint_mod
!
! Description: This file contains the setup for reciprocal routines and
!              fft code, all the sub calls for all the subroutines
!
! NOTE NOTE NOTE: This code assumes frc_int .eq. 0 and should only be used
!                 under these conditions!!!
!*******************************************************************************

module pme_recip_midpoint_mod

#ifdef MPI

  use processor_mod
  use pme_recip_dat_mod
#include "include_precision.i"

  implicit none

  integer, save     :: mpi_comm_x, mpi_comm_xy, mpi_comm_z
  integer           :: fft_x_dim, fft_y_dim, fft_z_dim
  integer*8         :: forward_plan(1)
  integer*8         :: backward_plan(1)    
! double precision, allocatable  ::q_x(:),q_y(:) 
! double precision, allocatable  ::q_z(:) 
! double precision, allocatable  :: mpi_send_buf(:), mpi_recv_buf(:)
!we dont need the following
! double precision, allocatable  :: dbl_mpi_send_bufx(:), dbl_mpi_recv_bufx(:)
! double precision, allocatable  :: dbl_mpi_send_bufy(:), dbl_mpi_recv_bufy(:)
! double precision, allocatable  :: dbl_mpi_send_bufz(:), dbl_mpi_recv_bufz(:)
  integer, allocatable           :: new_fft_x_cnts(:) 
  integer, allocatable           :: new_fft_x_offs(:) 
  integer, allocatable           :: new_fft_y_cnts(:)
  integer, allocatable           :: new_fft_y_offs(:)
  integer, allocatable           :: new_fft_y_offs2(:)
  integer, allocatable           :: i_tbl(:), j_tbl(:), k_tbl(:)
#ifdef pmemd_SPDP
  real,allocatable,save             :: xyz_q(:,:,:)
#else
  double precision,allocatable,save             :: xyz_q(:,:,:)
#endif
  integer           :: x_cnt, x_off, y_cnt, y_off, z_cnt, z_off

contains

!*******************************************************************************
!
! Subroutine:  pme_recip_setup
!
! Description:
!              
!*******************************************************************************

subroutine pme_recip_setup()


  use processor_mod, only: proc_dimx_recip, proc_dimy_recip, proc_dimz_recip, &
                           max_xfft_cnts, max_yfft_cnts, max_zfft_cnts, &
                           blockunflatten_recip, y_offset, &
                           proc_min_nfft1, proc_max_nfft1, &
                           proc_min_nfft2, proc_max_nfft2, &
                           proc_min_nfft3, proc_max_nfft3
  use pme_fft_midpoint_mod
  use gbl_constants_mod
  use pmemd_lib_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use pbc_mod
  implicit none
!for FFTW use
!#include <fftw3.f> 

  
  integer          :: color, key
  integer          :: nproc_x, nproc_xy, nproc_z
  integer          :: my_x_rank, my_xy_rank, my_z_rank
  integer          :: ierror
  integer          :: nlocalx, nlocaly, nlocalz

!the following should come from Charlies code
! the following variables are for fft grid range for my MPI rank
  integer           :: ngrid(3)
  integer           :: x_start, x_end, x_local
  integer           :: y_start, y_end, y_local
  integer           :: z_start, z_end, z_local
  integer           :: localmax
  integer           :: index_x, index_y, index_z
  integer           :: nfftmax,nx,ny,nz 
  integer*8         :: plan(3), plan_x, plan_y, plan_z
  integer           :: siz_buf1, siz_buf2, siz_buf3
  integer           :: max_proc_cnt
!For new fft x division
  integer           :: quot,rem, cnts_idx
  integer           :: new_fft_x_cntx(proc_dimx_recip) 
  integer           :: my_x, my_y, my_z
  integer           :: base_y_fft
  integer           :: i,j,k,y_idx,cnt_idx
!Index table for fill charge grid & grad sum
!  integer           :: x_cnt, x_off, y_cnt, y_off, z_cnt, z_off
  integer           :: jbot, jtop, kbot, ktop, order
#ifndef try
  integer           :: num_int=100, num_real=100  ! this is for the call pme_dat_setup_midpoint
#endif

! if(mytaskid_recip .eq. 0) print *, "proc dimx",proc_dimx_recip, "proc_dimy_recip",proc_dimy_recip , "proc_dimz_recip",proc_dimz_recip
#if 0
!MPI sub communicator in x direction, total division = x * y
  color = mytaskid_recip / (proc_dimz_recip) ! # of colors=# of rows ! (proc_dimy_recip * proc_dimz_recip) `
  key = mytaskid_recip
  call mpi_comm_split(mpi_comm_world, color, key, mpi_comm_x, ierror) 
  call mpi_comm_size (mpi_comm_x, nproc_x, ierror)
  call mpi_comm_rank (mpi_comm_x, my_x_rank, ierror)
!MPI sub communicator in xy direction
  color = mytaskid_recip / (proc_dimx_recip * proc_dimz_recip) ! may need to verify)
  key = mytaskid_recip
  call mpi_comm_split(mpi_comm_world, color, key, mpi_comm_xy, ierror) 
  call mpi_comm_size (mpi_comm_xy, nproc_xy, ierror)
  call mpi_comm_rank (mpi_comm_xy, my_xy_rank, ierror)
!MPI sub communicator in z direction
  color = mytaskid_recip / (proc_dimx_recip *  proc_dimy_recip)
  key = mytaskid_recip
  call mpi_comm_split(mpi_comm_world, color, key, mpi_comm_z, ierror)
  call mpi_comm_size (mpi_comm_z, nproc_z, ierror)
  call mpi_comm_rank (mpi_comm_z, my_z_rank, ierror)
#endif /*all the split code is commented*/
  
  allocate(new_fft_x_cnts(0:proc_dimx_recip*proc_dimy_recip - 1))
  allocate(new_fft_x_offs(0:proc_dimx_recip*proc_dimy_recip - 1))

  allocate(new_fft_y_cnts(0:proc_dimx_recip*proc_dimy_recip-1))!for initial fft transfers in x dim
  allocate(new_fft_y_offs(0:proc_dimx_recip*proc_dimy_recip-1))!for fft transfers in x dim
  allocate(new_fft_y_offs2(0:proc_dimx_recip*proc_dimy_recip-1))!for fft transfers in x dim
  fft_x_dim = (nfft1+2)/2
  fft_y_dim = nfft2
  fft_z_dim = nfft3

  fft_halo_size= max(es_cutoff,vdw_cutoff)+skinnb+bspl_order+1 ! use equation here
  if(mod(fft_halo_size,2) .eq. 0) fft_halo_size = fft_halo_size + 1 ! make it odd

  if(.not. allocated(xyz_q)) allocate(xyz_q(-1*fft_halo_size:2 *fft_x_dim+fft_halo_size+1, &
                                 -1*fft_halo_size:proc_max_nfft2-proc_min_nfft2+fft_halo_size+1, &
                                 -1*fft_halo_size:proc_max_nfft3-proc_min_nfft3+fft_halo_size+1))
  xyz_q(-1*fft_halo_size:2 * fft_x_dim+fft_halo_size+1, &
                                 -1*fft_halo_size:proc_max_nfft2-proc_min_nfft2+fft_halo_size+1, &
                                 -1*fft_halo_size:proc_max_nfft3-proc_min_nfft3+fft_halo_size+1) = 0.d0

!the following division in fft_x_dim is necessary in kspace
!dividing the x dim using fft_x_dim 
!here we are dividing fft_x_dim into proc_dimx_recip*proc_dimy_recip division 
!So before the FFT for x_dim, we use the division in processor.F90, after the fft start (kspace) we
!use the followng new didivison
!z we use the previous division 
!new y division is given below

! danger if fft wrong size

if(fft_x_dim .le. proc_dimx_recip*proc_dimy_recip) print *, "fft_x_dim is less than processor in xy dimension"
if(fft_y_dim .le. proc_dimx_recip*proc_dimy_recip) print *, "fft_y_dim is less than processor in xy dimension"

  quot = fft_x_dim / (proc_dimx_recip*proc_dimy_recip)
  rem = fft_x_dim - quot * proc_dimx_recip*proc_dimy_recip

  new_fft_x_cnts(0:proc_dimx_recip*proc_dimy_recip-1) = quot

  if (rem .ne. 0) then
    do i = 1, rem
      cnts_idx = proc_dimx_recip*proc_dimy_recip * i / rem -1 ! -1 to make it 0:proc_dimx_recip-1 
      new_fft_x_cnts(cnts_idx) = quot + 1
    end do
  end if
  offset = 0
  do i = 0, proc_dimx_recip*proc_dimy_recip - 1 
    new_fft_x_offs(i) = offset
    offset = offset + new_fft_x_cnts(i)
  end do

!the following y_fft_dim division will be used in k space
!in real space y_fft_dim is already divided into proc_dimy_recip in processor.F90
!but in kspace y_fft_dim will be divided by (proc_dimx_recip*proc_dimy_recip) division

!we will use already proc_dimy_recip divided fft_y_dim parts to further divided each part into proc_dimx_recip
!division
!the following division is little bit complicated than the above x division,
!since, each MPI in y directon only has a particular amount of y fft, so we take
!this amount and offset into account.  But in the
!above x fft division in k-spce we have full dim of x
  call blockunflatten_recip(mytaskid_recip, my_x,my_y,my_z)
  do y_idx = 0, proc_dimy_recip - 1
    base_y_fft = y_offset(y_idx)
    quot =  send_yfft_cnts(y_idx)/ proc_dimx_recip
    rem = send_yfft_cnts(y_idx) - quot * proc_dimx_recip
    !reasonable to rename new_fft_y_cnts to kspce_y_fft_cnts
    new_fft_y_cnts(y_idx*proc_dimx_recip+0:y_idx*proc_dimx_recip+proc_dimx_recip-1) = quot 
    !each division y_idx is further divided into proc_dimx_recip division in y direction
    if (rem .ne. 0 ) then
      do i =1,rem
        cnts_idx = y_idx* proc_dimx_recip + proc_dimx_recip * i / rem -1 ! -1 to make MPI index in y between 0:proc_dimx_recip-1 
        new_fft_y_cnts(cnts_idx) = quot + 1   
      end do 
    end if 
    offset = 0 ! this is a local offset for this MPI 
    !this local offset is requred, since berore I had a bigger chunk of y,
    ! it was not the full y dim. Now I need a smaller chunk of y 
    !that is my_yfft_cnt/proc_dimx_recip chunk. The array we build 
    !new_fft_y_cnts() is 
    !for full fft_y_dim. But my xyz_data only contain
    !my_yfft_cnt( = fft_y_dim/proc_dimy_recip) from the real space. but in 
    !k space I will even divide the my_yfft_cnt by proc_dimx_recip so that 
    ! no processor is idle, I will have one chunk=my_yfft_cnt/proc_dimx_recip
    !so offset just to make sure in my local space (1:my_yfft_cnt) from 
    !which value my smaller chunk of actual yfft starts 
    do i = 0, proc_dimx_recip-1 
      new_fft_y_offs(i+y_idx*proc_dimx_recip) = offset
      offset = offset + new_fft_y_cnts(i+y_idx*proc_dimx_recip)
    end do
  end do

  x_off = proc_min_nfft1 - 1
  x_cnt = proc_max_nfft1 - proc_min_nfft1 + 1

  y_off = proc_min_nfft2 - 1
  y_cnt = proc_max_nfft2 - proc_min_nfft2 + 1

  z_off = proc_min_nfft3 - 1
  z_cnt = proc_max_nfft3 - proc_min_nfft3 + 1

  ! Initialize the indexing tables.  It actually produces faster code to
  ! do this here rather than caching the info remotely.

  ! The following i_tbl, j_tbl, k_tbl take an input from -order to nfft+order.
  ! The output is a value from 1 to nfft (it corresponds to the grid space).

  ! We set a table that tells us how our i value corresponds to the grid space.
  order = bspl_order
  allocate(i_tbl(-order:nfft1+order),j_tbl(-order:nfft2+order),k_tbl(-order:nfft3+order))

  do i= -order, nfft1
    i_tbl(i) = i + 1
    if(i_tbl(i) .ge. proc_min_nfft1 .and. i_tbl(i) .le. proc_max_nfft1) then
      i_tbl(i) = i + 1
    else
      i_tbl(i) = 0
    end if
  end do
  ! If we go over nfft1 we wrap it back to nfft1-1, nfft1-2, etc (because
  ! ifracts starts from -order).  So for example if nfft1 is 45, when i =
  ! nfft-1, i is 44, when i=nfft1, i is 45, and when i is nfft+1, i is 44. This
  ! is important because the upper right portion of the code needs to contribute
  ! back into the grid charge space.
  do i= nfft1+1, nfft1+order
!   i_tbl(i) = nfft1 - (i - nfft1) + 1
    i_tbl(i) = 0
  end do

  jbot = 1
  jtop = y_cnt
  ! We set a table that moves our real space into our y space (because we only
  ! take a portion of y space we need an offset)
  do j = -order, nfft2
    j_tbl(j) = j - y_off + 1
    if(j_tbl(j) .ge. 1 .and. j_tbl(j) .le. y_cnt) then
    else
      j_tbl(j) = 0
    end if
  end do
  ! We set our y values to go in reverse pass nfft3
  do j=nfft2+1, nfft2+order
    j_tbl(j) = 0
  end do

  kbot = 1
  ktop = z_cnt
  ! We set up a table that moves our real space into our z space (because we
  ! only takea portion of the z space there's an offset)
  do k = -order, nfft3
    k_tbl(k) = k - z_off + 1
    if(k_tbl(k) .ge. 1 .and. k_tbl(k) .le. z_cnt) then
    else
      k_tbl(k) = 0
    end if
  end do
  ! We set our z values to go in reverse pass nfft3
  do k = nfft3+1, nfft3+order
    k_tbl(k) = 0
  end do

  !interface of midpoint real space to old code k space:
  !setup the data structure and variables to use old k space code
  !call blk_fft_setup_midpoint(x_dim, y_dim, z_dim,new_fft_x_cnts, &
  call blk_fft_setup_midpoint(fft_x_dim, fft_y_dim, fft_z_dim,new_fft_x_cnts, &
                          new_fft_x_offs, new_fft_y_cnts, new_fft_y_offs)
#ifndef try
  call pme_recip_dat_setup(num_int, num_real)
  call fft_dat_setup(num_int, num_real)
#endif
#if 0
  if(mytaskid_recip .eq. 0) print *, "fft total dims", nfft1, nfft2, nfft2
  if(mytaskid_recip .eq. 0) print *, "send_x cnts",send_xfft_cnts,"y cnts", send_yfft_cnts, "z cnts",send_zfft_cnts
  if(mytaskid_recip .eq. 0) print *, "send_offs",x_offset,"y offs",y_offset,"z offs", z_offset

  if(mytaskid_recip .eq. 0) print *, "new_x cnts",new_fft_x_cnts,"y cnts",new_fft_y_cnts,"z cnts",send_zfft_cnts
  if(mytaskid_recip .eq. 0) print *, "new_offs",new_fft_x_offs,"y offs",new_fft_y_offs,"z offs",z_offset
#endif
  return

end subroutine pme_recip_setup

!*******************************************************************************
!
! Subroutine:  do_kspace_midpoint
!
! Description:  <TBS>
!
! INPUT:
!
! OUTPUT:
!
! frc:      Forces incremented by k-space sum.
! virial:       Virial due to k-space sum (valid for atomic scaling;
!               rigid molecule virial needs a correction term not computed here.
!*******************************************************************************

subroutine do_kspace_midpoint(nb_frc, eer, virial, need_pot_enes, need_virials)

  use pme_fft_midpoint_mod
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use pbc_mod
  use prmtop_dat_mod
  use timers_mod
  use processor_mod
  use pme_direct_mod

  implicit none

  ! Formal arguments:

  double precision      :: nb_frc(3,proc_num_atms+proc_ghost_num_atms)
  double precision      :: eer
  double precision      :: virial(3, 3)
  logical               :: need_pot_enes
  logical               :: need_virials

  ! Local variables:

  integer           :: my_chgs(proc_atm_alloc_size)
  integer           :: my_chg_cnt
#ifdef pmemd_SPDP
  real              :: theta(bspl_order * 3 * proc_atm_alloc_size)
  real              :: dtheta(bspl_order * 3 * proc_atm_alloc_size)
  !real      :: xyz_q(2 * fft_x_dim, &
  !                   proc_max_nfft2-proc_min_nfft2+1, &
  !                   proc_max_nfft3-proc_min_nfft3+1)
  real      :: yxz_q(2 * fft_y_dim, &
                     fft_x_cnts(my_grid_idx1), &
                     fft_z_cnts(my_grid_idx2))
  real      :: zxy_q(2*fft_z_dim, &
                     fft_x_cnts(my_grid_idx1), &
                     fft_y_cnts2(my_grid_idx2))
#else
  double precision      :: theta(bspl_order * 3 * proc_atm_alloc_size)
  double precision      :: dtheta(bspl_order * 3 * proc_atm_alloc_size)
 ! double precision      :: xyz_q(-11:2 * fft_x_dim+12, &
 !                                -11:proc_max_nfft2-proc_min_nfft2+1+11, &
 !                                -11:proc_max_nfft3-proc_min_nfft3+1+11)
  double precision      :: yxz_q(2 * fft_y_dim, &
                                fft_x_cnts(my_grid_idx1), &
                                fft_z_cnts(my_grid_idx2))
                                ! proc_max_nfft1-proc_min_nfft1+1, &
                                ! proc_max_nfft3-proc_min_nfft3+1)
  double precision      :: zxy_q(2*fft_z_dim, &
                             fft_x_cnts(my_grid_idx1), &
                             fft_y_cnts2(my_grid_idx2))
#endif
  integer               :: ifracts(3, proc_atm_alloc_size)
  integer               :: i,j,k, my_x, my_y, my_z

   if (ntp .gt. 0 .and. is_orthog .ne. 0) call load_m_exp_tbls 

  call get_grid_weights(ifracts, theta, dtheta, bspl_order, my_chgs, my_chg_cnt)

  call update_pme_time(bspline_timer)

  call fill_charge_grid(xyz_q, theta, bspl_order, ifracts, &
                        my_chgs, my_chg_cnt)

  call update_pme_time(grid_charges_timer)

#if 0
call blockunflatten_recip(mytaskid_recip,my_x,my_y,my_z)
do k=1, my_zfft_cnt
  do j=1, my_yfft_cnt
    do i=proc_min_nfft1, proc_max_nfft1
      write(200+mytaskid_recip,'(3I0.4,ES20.10)') i,y_offset(my_y)+j,z_offset(my_z)+k,xyz_q(i,j,k)
    end do
  end do
end do
call mpi_barrier(pmemd_comm, err_code_mpi)
stop 1
#endif

 call fft3d_forward_midpoint(xyz_q,yxz_q,zxy_q,new_fft_y_cnts, new_fft_y_offs, &
                                    fft_x_dim,fft_y_dim, fft_z_dim, &
                                    new_fft_x_cnts, new_fft_x_offs)
  call update_pme_time(fft_timer)

!After MPI comm , verifying output
!call blockunflatten_recip(mytaskid_recip,my_x,my_y,my_z)

#if 0
!print *, "new after for z",fft_z_dim,fft_x_offs(my_grid_idx1),fft_x_cnts(my_grid_idx1), fft_y_offs2(my_grid_idx2),fft_y_cnts2(my_grid_idx2), myblockid
  do k=1, fft_y_cnts2(my_grid_idx2)
    do j=1, fft_x_cnts(my_grid_idx1)
      do i=1, 2*fft_z_dim
          write(300+myblockid,'(3I0.4,ES18.8)') i,fft_x_offs(my_grid_idx1)+j,fft_y_offs2(my_grid_idx2)+k,zxy_q(i,j,k) 
      end do
     end do
   end do
call mpi_barrier(pmemd_comm, err_code_mpi)
stop 1
#endif
  if (is_orthog .ne. 0) then

    if (need_virials) then
      call scalar_sumrc_uv(zxy_q, ew_coeff, uc_volume, &
                           gbl_prefac1, gbl_prefac2, gbl_prefac3, &
                           fft_x_dim, fft_y_dim, fft_z_dim, eer, virial)
    else if (need_pot_enes) then
!      call scalar_sumrc_u(zxy_q, ew_coeff, uc_volume, &
!                          gbl_prefac1, gbl_prefac2, gbl_prefac3, &
!                          fft_x_dim, fft_y_dim, fft_z_dim, eer)
! Ross Walker - 7/8/2014 - bug bug - The PME Error estimate requires the vir
! since it compares virial vs energy.
! pme_err_est =   vir_vs_ene = vir%elec_recip(1, 1) + &
!               vir%elec_recip(2, 2) + &
!               vir%elec_recip(3, 3) + &
!               vir%eedvir + &
!               vir%elec_nb_adjust(1, 1) + &
!               vir%elec_nb_adjust(2, 2) + &
!               vir%elec_nb_adjust(3, 3)
! So we can't just naively skip it here unless we choose not to print the
! Ewald error estimate anymore (like in the GPU code). For now we will calculate
! the reciprocal virial whenever need_pot_enes is true.
      call scalar_sumrc_uv(zxy_q, ew_coeff, uc_volume, &
                           gbl_prefac1, gbl_prefac2, gbl_prefac3, &
                           fft_x_dim, fft_y_dim, fft_z_dim, eer, virial)
    else
      call scalar_sumrc(zxy_q, ew_coeff, uc_volume, &
                        gbl_prefac1, gbl_prefac2, gbl_prefac3, &
                        fft_x_dim, fft_y_dim, fft_z_dim)
    end if
  else
    if (need_virials) then

      call scalar_sumrc_nonorthog_uv(zxy_q, ew_coeff, uc_volume, &
                                     gbl_prefac1, gbl_prefac2, gbl_prefac3, &
                                     fft_x_dim, fft_y_dim, fft_z_dim, &
                                     eer, virial)
    else if (need_pot_enes) then
!      call scalar_sumrc_nonorthog_u(zxy_q, ew_coeff, uc_volume, &
!                                    fft_x_dim, fft_y_dim, fft_z_dim, &
!                                    eer)
! Ross Walker - 7/8/2014 - bug bug - The PME Error estimate requires the vir
! since it compares virial vs energy.
! pme_err_est =   vir_vs_ene = vir%elec_recip(1, 1) + &
!               vir%elec_recip(2, 2) + &
!               vir%elec_recip(3, 3) + &
!               vir%eedvir + &
!               vir%elec_nb_adjust(1, 1) + &
!               vir%elec_nb_adjust(2, 2) + &
!               vir%elec_nb_adjust(3, 3)
! So we can't just naively skip it here unless we choose not to print the
! Ewald error estimate anymore (like in the GPU code). For now we will calculate
! the reciprocal virial whenever need_pot_enes is true.
      call scalar_sumrc_nonorthog_uv(zxy_q, ew_coeff, uc_volume, &
                                     gbl_prefac1, gbl_prefac2, gbl_prefac3, &
                                     fft_x_dim, fft_y_dim, fft_z_dim, &
                                     eer, virial)
    else
      call scalar_sumrc_nonorthog(zxy_q, ew_coeff, uc_volume, &
                                  gbl_prefac1, gbl_prefac2, gbl_prefac3, &
                                  fft_x_dim, fft_y_dim, fft_z_dim)
    end if
  end if
  call update_pme_time(scalar_sum_timer)
#if 0
!scalar sum debug
!call blockunflatten_recip(mytaskid_recip,my_x,my_y,my_z)
  do k=1, nfft3*2
    do j=1, fft_y_cnts2(my_grid_idx2)
      do i=1, fft_x_cnts(my_grid_idx1)
         write(500+mytaskid_recip,'(3I0.4,ES13.3)')fft_x_offs(my_grid_idx1)+i,fft_y_offs2(my_grid_idx2)+j,k,zxy_q(k,i,j)
      end do
    end do
 end do
#endif

  call fft3d_back_midpoint(zxy_q, xyz_q, fft_x_dim, fft_y_dim, fft_z_dim, &
                                          new_fft_y_cnts, new_fft_y_offs)
  call update_pme_time(fft_timer)
#if 0
  do k=1, fft_z_cnts(my_grid_idx2)
    do j=1, fft_y_cnts1(my_grid_idx1)
      do i=1, 2*fft_x_dim
!          write(DBG_ARRAYS+mytaskid_recip,'(I0.2,3(" ",I0.4),":",F)') 1, i,fft_y_offs1(my_grid_idx1)+j,fft_z_offs(my_grid_idx2)+k, xyz_q(i,j,k) 
      end do
     end do
   end do
#endif

  call grad_sum(nb_frc, xyz_q, theta, dtheta, bspl_order, ifracts, my_chgs, my_chg_cnt)
  call update_pme_time(grad_sum_timer)

!#if 0
!!if (counts .eq. 1) then
!  do i=1, proc_num_atms+proc_ghost_num_atms
!    write(15+mytaskid_recip,'(A,I5,3x,F12.8,3x,F12.8,3x,F12.8,3x,A,1x,I5)')"new",i,&
!    !    write(15+mytaskid_recip,'(A,I5,4x,F10.6,4x,F10.6,4x,F10.6,4x,A,1x,I5)')"new",i,&
!           proc_atm_frc(1,i), proc_atm_frc(2,i),proc_atm_frc(3,i),'fullid',proc_atm_to_full_list(i)
!  end do
!!end if
!#endif

!   call proc_data_dump(pbc_box, proc_num_atms+proc_ghost_num_atms)

  return

end subroutine do_kspace_midpoint

!*******************************************************************************
!
! Subroutine:  get_grid_weights
!
! Description: <TBS>
!              
!*******************************************************************************
subroutine get_grid_weights(ifracts, theta, dtheta, order,      my_chgs, my_chg_cnt)

  use bspline_mod
  use pme_fft_midpoint_mod
  use gbl_datatypes_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use pbc_mod
  use processor_mod

  implicit none

! Formal arguments:

  integer               :: ifracts(3, proc_atm_alloc_size)
  integer               :: order
#ifdef pmemd_SPDP
!  real      :: atm_crd(4, *) !it is a packed crd+q array, so 4
  real      :: theta(order, 3, proc_atm_alloc_size)
  real      :: dtheta(order, 3, proc_atm_alloc_size)
#else
!  double precision      :: atm_crd(3, *) !it is a packed crd+q array, so 4
  double precision      :: theta(order, 3, proc_atm_alloc_size)
  double precision      :: dtheta(order, 3, proc_atm_alloc_size)
#endif
  integer               :: my_chgs(proc_atm_alloc_size), my_chg_cnt

! Local variables:

#ifdef pmemd_SPDP
  real      :: box_x, box_y, box_z
  real      :: crd_x, crd_y, crd_z
  real      :: factor1, factor2, factor3
  real      :: fract(3)
  real      :: weight(3)
  real      :: ucellorder1, ucellorder2, ucellorder3
#else
  double precision      :: box_x, box_y, box_z
  double precision      :: crd_x, crd_y, crd_z
  double precision      :: factor1, factor2, factor3
  double precision      :: fract(3)
  double precision      :: weight(3)
  double precision      :: ucellorder1, ucellorder2, ucellorder3
#endif
  integer               :: i, j, k, x, y, z
#ifdef pmemd_SPDP
  real, parameter   :: one_third = 1.0 / 3.0
  real, parameter   :: one_half=0.50, one=1.0, two=2.0, three=3.0
#else
  double precision, parameter   :: one_third = 1.d0 / 3.d0
  double precision, parameter   :: one_half=0.5d0, one=1.d0, two=2.d0, three=3.d0
#endif

  integer               :: lst_idx
  integer               :: i_base, j_base, k_base

  my_chg_cnt = 0

  box_x = ucell(1, 1)
  box_y = ucell(2, 2)
  box_z = ucell(3, 3)

  ! Scaling factors to get from cit table indexes to nfft indexes:

  factor1 = dble(nfft1) * recip(1, 1)
  factor2 = dble(nfft2) * recip(2, 2)
  factor3 = dble(nfft3) * recip(3, 3)

  ! This sets an upper limit for our nfft grid which is box size + order*nfft
  ! grid space
  ucellorder1 = ucell(1,1) + order*ucell(1,1)/dble(nfft1)
  ucellorder2 = ucell(2,2) + order*ucell(2,2)/dble(nfft2)
  ucellorder3 = ucell(3,3) + order*ucell(3,3)/dble(nfft3)
  do  i = 1, proc_num_atms + proc_ghost_num_atms

#ifdef pmemd_SPDP
    crd_x = proc_crd_q_sp(1, i)
    crd_y = proc_crd_q_sp(2, i)
    crd_z = proc_crd_q_sp(3, i)
#else
    crd_x = proc_atm_crd(1, i)
    crd_y = proc_atm_crd(2, i)
    crd_z = proc_atm_crd(3, i)
#endif

    ! This if statement makes it so we only look at atoms that start from 0,0,0
    ! to ucell + nfft length * order.  This makes it so we deal with a buffer
    ! that is the order.
    if(crd_x .ge. 0 .and. crd_x .lt. ucellorder1 .and. &
       crd_y .ge. 0 .and. crd_y .lt. ucellorder2 .and. &
       crd_z .ge. 0 .and. crd_z .lt. ucellorder3) then
      !   i=atm_lst(i)%nxt
      !   cycle
      !end if
      fract(1) = factor1 * crd_x
      fract(2) = factor2 * crd_y
      fract(3) = factor3 * crd_z
    
      !    k = int(fract(3))
    
      my_chg_cnt = my_chg_cnt + 1
      my_chgs(my_chg_cnt) = i
      weight(:) = fract(:) - int(fract(:))
    
      ifracts(:, my_chg_cnt) = int(fract(:)) - order
            
      theta(1, 1:3, my_chg_cnt ) = one - weight(:)
      theta(2, 1:3, my_chg_cnt ) = weight(:)
      do j = 1, 3
        ! One pass to order 3:
        theta(3, j, my_chg_cnt ) = one_half * weight(j) * theta(2, j,my_chg_cnt )
        theta(2, j,my_chg_cnt ) = one_half * ((weight(j) + one) * theta(1, j,my_chg_cnt ) + &
                      (two - weight(j)) * theta(2, j,my_chg_cnt ))
        theta(1, j,my_chg_cnt ) = one_half * (one - weight(j)) * theta(1, j,my_chg_cnt )

        ! Diff to get dtheta:
        dtheta(1, j ,my_chg_cnt ) = -theta(1, j ,my_chg_cnt )
        dtheta(2, j,my_chg_cnt ) = theta(1, j,my_chg_cnt ) - theta(2, j,my_chg_cnt )
        dtheta(3, j,my_chg_cnt ) = theta(2, j,my_chg_cnt ) - theta(3, j,my_chg_cnt )
        dtheta(4, j,my_chg_cnt ) = theta(3, j,my_chg_cnt )

        ! One final pass to order 4:
        theta(4, j,my_chg_cnt ) = one_third * weight(j) * theta(3, j,my_chg_cnt )
        theta(3, j,my_chg_cnt ) = one_third * ((weight(j) + one) * theta(2, j,my_chg_cnt ) + &
                      (three - weight(j)) * theta(3, j,my_chg_cnt ))
        theta(2, j,my_chg_cnt ) = one_third * ((weight(j) + two) * theta(1, j,my_chg_cnt ) + &
                      (two - weight(j)) * theta(2, j,my_chg_cnt ))
        theta(1, j,my_chg_cnt ) = one_third * (one - weight(j)) * theta(1, j,my_chg_cnt )
      end do

    end if
  end do


!#endif
  return

end subroutine get_grid_weights

!*******************************************************************************
!
! Subroutine:  fill_charge_grid
!
! Description: <TBS>
!
! INPUT:
!
! theta1, theta2, theta3:       Spline coeff arrays.
! ifracts:                      int(scaled fractional coords).
! nfft1, nfft2, nfft3:          Logical charge grid dimensions.
!
! fft_x_dim, fft_y_dim, fft_z_dim: Physical charge grid dims.
!
! order:                    Order of spline interpolation.
!
! OUTPUT:
!
! q:                            Charge grid.
!              
!*******************************************************************************

subroutine fill_charge_grid(q, theta, order, ifracts, &
                            my_chgs, my_chg_cnt)

  use mdin_ewald_dat_mod
  use pme_fft_midpoint_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use processor_mod, only : proc_atm_qterm, proc_min_nfft1, proc_max_nfft1, &
                            proc_min_nfft2, proc_max_nfft2, proc_min_nfft3, &
                            proc_max_nfft3
! BEGIN DBG
  use pmemd_lib_mod
! END DBG

  implicit none

! Formal arguments:

  integer               :: order
  integer               :: ifracts(3, proc_atm_alloc_size)
#ifdef pmemd_SPDP
  real      :: theta(order, 3, *)
  real      :: q(-1*fft_halo_size : 2 * fft_x_dim+fft_halo_size+1, &
                           -1*fft_halo_size:proc_max_nfft2-proc_min_nfft2+fft_halo_size+1,&
                            -1*fft_halo_size:proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)
#else
  double precision      :: theta(order, 3, *)
  double precision      :: q(-1*fft_halo_size : 2 * fft_x_dim + fft_halo_size+1, &
                           -1*fft_halo_size:proc_max_nfft2-proc_min_nfft2+fft_halo_size+1,&
                            -1*fft_halo_size:proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)
#endif
  integer               :: my_chgs(proc_atm_alloc_size), my_chg_cnt
! Local variables:

  integer               :: i, j, k
  integer               :: i_base, j_base, k_base
  integer               :: ith1, ith2, ith3
  integer               :: my_chg_idx
#ifdef pmemd_SPDP
  real      :: charge
  real      :: i_term, k_term, j_term
  real      :: tmp_theta(1:4)
#else
  double precision      :: charge
  double precision      :: i_term, k_term, j_term
  double precision      :: tmp_theta(1:4)
#endif
!  integer               :: x_off, y_off, z_off
!  integer               :: x_cnt, y_cnt, z_cnt
  integer               :: ibot, itop
  integer               :: jbot, jtop
  integer               :: kbot, ktop

  q(1:2*fft_x_dim , 1:y_cnt , 1:z_cnt) = 0.0
   !since the q starts from -ve co-ordinate

  ! We special-case order 4, the default.
  if (order .eq. 4) then
    do my_chg_idx = 1, my_chg_cnt

#ifdef pmemd_SPDP
      charge = proc_crd_q_sp(4,my_chgs((my_chg_idx)))
#else
      charge = proc_atm_qterm(my_chgs((my_chg_idx)))
#endif

      i_base = ifracts(1, my_chg_idx) + 1
      j_base = ifracts(2, my_chg_idx) - y_off + 1
      k_base = ifracts(3, my_chg_idx) !- z_off + 1
      tmp_theta(1:4) =   theta(1:4, 1, my_chg_idx) 

!omp simd collapse(2) private(k,k_term,j_term) reduction(+:q)
      do ith3 = 1, 4
        k = k_tbl(k_base + ith3)
       if(k .gt. 0) then
        k_term =  theta(ith3, 3, my_chg_idx) * charge
   !     !$dir  simd
   !     !$dir  vector always
        do ith2 = 1, 4
          j_term = k_term * theta(ith2, 2, my_chg_idx)
        !do ith2 = j_base+1, j_base+4
        !  j_term = k_term * theta(ith2-j_base, 2, my_chg_idx)
           q(i_base+1:i_base+4, j_base+ith2, k) = &
                            q(i_base+1:i_base+4,j_base+ith2, k) +  tmp_theta(1:4)*j_term
        end do
        end if ! for k
      end do
    end do
  else
    do my_chg_idx = 1, my_chg_cnt

      charge = proc_atm_qterm(my_chgs(my_chg_idx))

      i_base = ifracts(1, my_chg_idx)
      j_base = ifracts(2, my_chg_idx)
      k_base = ifracts(3, my_chg_idx)

      do ith3 = 1, order

        k = k_tbl(k_base + ith3)

        if (k .ne. 0) then

          k_term =  theta(ith3, 3, my_chg_idx) * charge
          do ith2 = 1, order
            j = j_tbl(j_base + ith2)
            if (j .ne. 0) then
              j_term = k_term * theta(ith2, 2, my_chg_idx)
              do ith1 = 1, order
                i = i_tbl(i_base + ith1)
                q(i, j, k) = q(i, j, k) + theta(ith1, 1, my_chg_idx) * j_term
              end do
            end if
          end do
        end if
      end do
    end do
  end if
  return

end subroutine fill_charge_grid

!*******************************************************************************
!
! Subroutine:  grad_sum
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine grad_sum(frc, q, theta, dtheta, order, ifracts, &
                    my_chgs, my_chg_cnt)

  use img_mod
  use mdin_ewald_dat_mod
  use prmtop_dat_mod
  use pbc_mod
  use pme_fft_midpoint_mod
  use ti_mod
  use processor_mod, only : proc_atm_qterm, proc_min_nfft1, proc_max_nfft1, &
                            proc_min_nfft2, proc_max_nfft2, proc_min_nfft3, &
                            proc_max_nfft3, proc_num_atms, proc_ghost_num_atms

  implicit none

! Formal arguments:

  double precision      :: frc(3, proc_num_atms+proc_ghost_num_atms)
  integer               :: order
#ifdef pmemd_SPDP
  real      :: q(-1*fft_halo_size : 2 * fft_x_dim+fft_halo_size+1, &
                           -1*fft_halo_size:proc_max_nfft2-proc_min_nfft2+fft_halo_size+1,&
                            -1*fft_halo_size:proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)
  real      :: theta(order, 3, *)
  real      :: dtheta(order, 3, *)
#else
  double precision      :: q(-1*fft_halo_size : 2 * fft_x_dim+fft_halo_size+1, &
                           -1*fft_halo_size:proc_max_nfft2-proc_min_nfft2+fft_halo_size+1,&
                            -1*fft_halo_size:proc_max_nfft3-proc_min_nfft3+fft_halo_size+1)
  double precision      :: theta(order, 3, *)
  double precision      :: dtheta(order, 3, *)
#endif
  integer               :: ifracts(3, *)

  integer               :: my_chgs(*), my_chg_cnt

! Local variables:

  integer               :: i, j, k
  integer               :: i_base, j_base, k_base
  integer               :: ith1, ith2, ith3, jth
  integer               :: my_chg_idx
#ifdef pmemd_SPDP
  real      :: f1, f2, f3
  real      :: charge, qterm(4), temp1, temp2
  real      :: f1_term(4), f2_term(4), f3_term(4)
  double precision      :: dfx, dfy, dfz
  real      :: recip_11, recip_22, recip_33
  real      :: dnfft1, dnfft2, dnfft3
#else
  double precision      :: f1, f2, f3
  double precision      :: charge, qterm(4),temp1, temp2
  double precision      :: f1_term(4), f2_term(4), f3_term(4)
  double precision      :: dfx, dfy, dfz
  double precision      :: recip_11, recip_22, recip_33
  double precision      :: dnfft1, dnfft2, dnfft3
#endif
  integer               :: atm_i

  recip_11 = recip(1, 1)
  recip_22 = recip(2, 2)
  recip_33 = recip(3, 3)

  dnfft1 = dble(nfft1)
  dnfft2 = dble(nfft2)
  dnfft3 = dble(nfft3)

!zero them for not using i .gt. 0 and not using the i_tbl
   q(proc_min_nfft1 -1*fft_halo_size : proc_min_nfft1-1, &
     1:proc_max_nfft2-proc_min_nfft2+1,&
     1:proc_max_nfft3-proc_min_nfft3+1) = 0.d0
  q(proc_max_nfft1+1: proc_max_nfft1+fft_halo_size+1, &
    1:proc_max_nfft2-proc_min_nfft2+1,&
   1:proc_max_nfft3-proc_min_nfft3+1) = 0.d0
!zero them for not using j .gt. 0 and not using the i_tbl
   q(proc_min_nfft1 -1*fft_halo_size : proc_max_nfft1+fft_halo_size+1, &
    -1*fft_halo_size:0,&
     1:proc_max_nfft3-proc_min_nfft3+1) = 0.d0
   q(proc_min_nfft1 -1*fft_halo_size : proc_max_nfft1+fft_halo_size+1, &
    proc_max_nfft2-proc_min_nfft2+2:proc_max_nfft2-proc_min_nfft2+fft_halo_size+1,&
     1:proc_max_nfft3-proc_min_nfft3+1) = 0.d0

  do my_chg_idx = 1, my_chg_cnt
    atm_i = my_chgs(my_chg_idx)

#ifdef pmemd_SPDP
      charge = proc_crd_q_sp(4,atm_i)
#else
      charge = proc_atm_qterm(atm_i)
#endif

    i_base = ifracts(1, my_chg_idx)+1
    j_base = ifracts(2, my_chg_idx)-y_off+1
    k_base = ifracts(3, my_chg_idx)

    f1 = 0.d0
    f2 = 0.d0
    f3 = 0.d0

    ! We special-case order 4, the default.

    if (order .eq. 4) then
      do ith3 = 1, 4

        k = k_tbl(k_base + ith3)
        if( k .gt. 0) then

            f1_term(1:4) = theta(1:4, 2, my_chg_idx) * theta(ith3, 3, my_chg_idx)
            f2_term(1:4) = dtheta(1:4, 2, my_chg_idx) * theta(ith3, 3, my_chg_idx)
            f3_term(1:4) = theta(1:4, 2, my_chg_idx) * dtheta(ith3, 3, my_chg_idx)

          do ith2 = 1, 4
            qterm(1:4) = q(i_base+1:i_base+4, j_base+ith2, k)
            f1 = f1 - f1_term(ith2) * SUM(qterm(1:4) * dtheta(1:4, 1, my_chg_idx))
            temp2 = SUM(qterm(1:4) * theta(1:4, 1, my_chg_idx))
            f2 = f2 - f2_term(ith2) * temp2
            f3 = f3 - f3_term(ith2) * temp2
          end do
        end if
      end do
    else
      do ith3 = 1, order

        k = k_tbl(k_base + ith3)

        if (k .ne. 0) then

          do ith2 = 1, order

            j = j_tbl(j_base + ith2)

            if (j .ne. 0) then

              f1_term(1) = theta(ith2, 2, my_chg_idx) * theta(ith3, 3, my_chg_idx)
              f2_term(1) = dtheta(ith2, 2, my_chg_idx) * theta(ith3, 3, my_chg_idx)
              f3_term(1) = theta(ith2, 2, my_chg_idx) * dtheta(ith3, 3, my_chg_idx)

              do ith1 = 1, order
                qterm(1) = q(i_tbl(i_base+ith1), j, k)

                ! Force is negative of grad:
                f1 = f1 - qterm(1) * dtheta(ith1, 1, my_chg_idx) * f1_term(1)
                f2 = f2 - qterm(1) * theta(ith1, 1, my_chg_idx) * f2_term(1)
                f3 = f3 - qterm(1) * theta(ith1, 1, my_chg_idx) * f3_term(1)
  
              end do
            end if

          end do

        end if
      end do
    end if

    f1 = f1 * dnfft1 * charge
    f2 = f2 * dnfft2 * charge
    f3 = f3 * dnfft3 * charge

    if (is_orthog .ne. 0) then
      dfx = recip_11 * f1
      dfy = recip_22 * f2
      dfz = recip_33 * f3
    else
      dfx = recip(1, 1) * f1 + recip(1, 2) * f2 + recip(1, 3) * f3
      dfy = recip(2, 1) * f1 + recip(2, 2) * f2 + recip(2, 3) * f3
      dfz = recip(3, 1) * f1 + recip(3, 2) * f2 + recip(3, 3) * f3
    end if

    frc(1, atm_i) = frc(1, atm_i) + dfx
    frc(2, atm_i) = frc(2, atm_i) + dfy
    frc(3, atm_i) = frc(3, atm_i) + dfz

  end do
  return

end subroutine grad_sum

#define POT_ENES 1
#define VIRIAL 1
#include "scalar_sum_midpoint.i"
#undef VIRIAL
#define VIRIAL 0
#include "scalar_sum_midpoint.i"
#undef POT_ENES
#define POT_ENES 0
#include "scalar_sum_midpoint.i"
#undef POT_ENES
#undef VIRIAL

#endif /*MPI*/

end module pme_recip_midpoint_mod
