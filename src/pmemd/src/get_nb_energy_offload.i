
#ifdef MIC_offload
!*******************************************************************************
!
! Subroutine:  get_nb_energy_offload
!
! Description:
!              
! The main routine for non bond energy (vdw and hbond) as well as direct part
! of ewald sum.  It is structured for parallelism.
!
!*******************************************************************************
subroutine get_nb_energy_offload(img_frc, img_crd, img_qterm, eed_cub, &
                         ipairs, tranvec, need_pot_enes, need_virials, &
                         eed, evdw, ehb, eedvir, virial, atm_cnt, new_list, min_img, max_img)
  use img_mod
  use timers_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use omp_lib 
  use offload_allocation_mod 
  use ene_frc_splines_mod

  implicit none

! Formal arguments:

  double precision, intent(in out) :: img_frc(3, *)
  double precision, intent(in)     :: img_crd(3, *)
  double precision , intent(in)    :: img_qterm(*)
  double precision, intent(in)     :: eed_cub(*)
  integer                          :: ipairs(*)
  double precision, intent(in)     :: tranvec(1:3, 0:17)
  logical, intent(in)              :: need_pot_enes
  logical, intent(in)              :: need_virials
  double precision, intent(out)    :: eed
  double precision, intent(out)    :: evdw
  double precision, intent(out)    :: ehb
  double precision, intent(out)    :: eedvir
  double precision, intent(out)    :: virial(3, 3)
! Local variables and parameters:

/* Newly added Formal argument for OpenMP and offload efficiency*/
integer		:: min_img, max_img /*range of neighbors*/
integer		:: atm_cnt
logical		:: new_list 

!DIR$ OPTIONS /offload_attribute_target=mic
integer		:: my_thread_id, j

/*OpenMP*/

  double precision,save ::      del
  double precision,save ::      dxdr
  double precision      eedtbdns_stk
  double precision      eedvir_stk, eed_stk, evdw_stk, ehb_stk
  double precision     max_nb_cut2, es_cut2, es_cut
  double precision      x_i, y_i, z_i
  double precision      x_tran(1:3, 0:17)
  double precision      vxx, vxy, vxz, vyy, vyz, vzz
  integer               i
  integer               ipairs_idx
  integer,save   ::               ntypes_stk
  integer               img_i
  integer               ee_eval_cnt
  integer               full_eval_cnt

!offload related variables
  integer		start_ipairs_host, end_ipairs_host
  integer		start_ipairs_dev, end_ipairs_dev
  integer	        start_img_crd_host, end_img_crd_host
  integer	        start_img_crd_dev, end_img_crd_dev	
  !double precision	::	mic_fraction
  integer		:: sig1, sig2 ! asynchronous variables
  integer			:: offset_atm
  

  logical			:: need_pot_enes_dev !in
  logical			:: need_virials_dev !in 
  double precision		:: eed_dev !inout
  double precision		:: evdw_dev !inout
  double precision		:: ehb_dev !inout
  double precision		:: eedvir_dev
  double precision		:: virial_dev(3, 3)
  double precision 		:: vxx_dev, vxy_dev, vxz_dev, vyy_dev, vyz_dev, vzz_dev
  double precision 		:: eedvir_stk_dev, eed_stk_dev, evdw_stk_dev, ehb_stk_dev

  integer               common_tran    ! flag - 1 if translation not needed
  logical               cutoffs_equal
  integer			:: max_local, min_local, atm_local ! This is done for offload coi buffer efficiency
!DIR$ END OPTIONS

  if (my_img_lo .gt. my_img_hi) return
!print *, "get nb energy starts",mytaskid

!offload asynch related signal tag
sig1=1
sig2=2


  ntypes_stk = ntypes

  eedvir_stk = 0.d0
  eed_stk = 0.d0
  evdw_stk = 0.d0
  ehb_stk = 0.d0

  dxdr = ew_coeff
  eedtbdns_stk = eedtbdns
  del = 1.d0 / eedtbdns_stk
  max_nb_cut2 = vdw_cutoff * vdw_cutoff
  es_cut = es_cutoff
  es_cut2 = es_cut * es_cut
  cutoffs_equal = (vdw_cutoff .eq. es_cutoff)

  vxx = 0.d0
 ! vxy = 0.d0
 ! vxz = 0.d0
  vyy = 0.d0
 ! vyz = 0.d0
  vzz = 0.d0

  ipairs_idx = 1

/*The following is required for OpenMP to work, this is an overhead, but happens one for each MPI process*/

/*divide the atom range into two parts: 1st part to MIC and 2nd part to host*/

  ipairs_idx = 1
/*the following are done due to compiler inefficiency for offloading*/
  atm_local = atm_cnt

if (offload_tasks(mytaskid)) then 
  max_local = max_img
  min_local = min_img
  start_img_crd_host = my_img_lo+ int(dble(my_img_hi - my_img_lo)*mic_fraction) 
  end_img_crd_host = my_img_hi ! end crd for host
  start_img_crd_dev = my_img_lo ! start crd for mic
  end_img_crd_dev = start_img_crd_host -1

  start_ipairs_dev = 1
  end_ipairs_dev = num_packed_dev /*not necessary*/
else
  start_img_crd_host = my_img_lo
  end_img_crd_host = my_img_hi ! end crd for host
end if
  start_ipairs_host = 1
  end_ipairs_host = num_packed 

if(new_list .eq. .true.) then
if (offload_tasks(mytaskid)) then 
!dir$ offload begin mandatory target (mic:0) &
 in(img_qterm(min_local:max_local) : alloc_if(.false.) free_if(.false.) into(img_qterm_dev(min_local:max_local))) & 
 in(gbl_img_iac(min_local:max_local) : alloc_if(.false.) free_if(.false.) into(gbl_img_iac_dev(min_local:max_local))) & 
  nocopy(thread_num_packed_dev(:) : alloc_if(.false.) free_if(.false.)) &
  nocopy(thread_ipairs_dev(0:ipairs_maxsize_dev) : alloc_if(.false.) free_if(.false.)) &
!  nocopy(ipairs_dev(:): alloc_if(.false.) free_if(.false.)) &
  nocopy(start_indexes_dev(:): alloc_if(.false.) free_if(.false.)) &
  nocopy(ipairs_idx,img_i,ee_eval_cnt,full_eval_cnt) &
  in( start_img_crd_dev,end_img_crd_dev, i, num_threads_offload, num_packed_dev,ipairs_maxsize_dev )
  ipairs_idx = 0
  num_packed_dev = 0
  i=0
	
  do img_i = start_img_crd_dev, end_img_crd_dev
      start_indexes_dev(img_i - start_img_crd_dev + 1) = ipairs_idx  ! starting from correct start
      ipairs_idx = ipairs_idx + 1
      ee_eval_cnt = thread_ipairs_dev(ipairs_idx)
      full_eval_cnt = thread_ipairs_dev(ipairs_idx + 1)
      ipairs_idx = ipairs_idx + 2
      ipairs_idx = ipairs_idx + ee_eval_cnt + full_eval_cnt
      /*ipairs_idx traverses with "nblist of each thread", and then jump to start from the nblist of next thread*/
      if (ipairs_idx .eq. num_packed_dev + thread_num_packed_dev(i)  .and. (i .lt. num_threads_offload)) then
        ipairs_idx =  ipairs_maxsize_dev*(i+1)/num_threads_offload 
        /*num_packed_dev helps to track when traversing of one thread is done*/
	num_packed_dev = ipairs_maxsize_dev*(i+1)/num_threads_offload 
	i = i + 1
     end if
      /*Need to take care if any OMP threads are not used in nb list,i.e.,thread_num_packed_dev(i) = 0*/
      do while(thread_num_packed_dev(i) .eq. 0 .and. (i .lt. num_threads_offload))
        ipairs_idx =  ipairs_maxsize_dev*(i+1)/num_threads_offload 
        /*num_packed_dev helps to track when traversing of one thread is done*/
	num_packed_dev = ipairs_maxsize_dev*(i+1)/num_threads_offload 
        i = i + 1
      end do
  end do
!dir$ end offload
 end if
end if

if (need_pot_enes) then
!1st part to the device
if (offload_tasks(mytaskid)) then 
!dir$ offload begin  mandatory signal(sig1) target (mic:0)  &

  out(img_frc_dev(:,1:atm_cnt): alloc_if(.false.) free_if(.false.))  &
  nocopy(img_frc_thread_dev(:,:): alloc_if(.false.) free_if(.false.))  &
  in(img_crd(:, 1:atm_cnt):  alloc_if(.false.) free_if(.false.) into(img_crd_dev(:,1:atm_cnt))) &
  nocopy(thread_ipairs_dev(0:ipairs_maxsize_dev): alloc_if(.false.) free_if(.false.)) &
  nocopy(start_indexes_dev(:): alloc_if(.false.) free_if(.false.)) &
  nocopy(x_tran_dev(:,: ): alloc_if(.false.) free_if(.false.)) &
  in(tranvec(:,0:17): alloc_if(.false.) free_if(.false.) into(tranvec_dev)) & ! may preloaded
  nocopy(img_qterm_dev(:) : alloc_if(.false.) free_if(.false.)) & ! may preloaded 
  nocopy(eed_cub_dev(:) : alloc_if(.false.) free_if(.false.)) &
  nocopy(gbl_cn1_dev(:) : alloc_if(.false.) free_if(.false.)) &
  nocopy(gbl_cn2_dev(:) : alloc_if(.false.) free_if(.false.)) & 
  nocopy(efs_tbl_dev(:) : alloc_if(.false.) free_if(.false.)) &
  nocopy(typ_ico_dev(:) : alloc_if(.false.) free_if(.false.)) &
  nocopy(gbl_img_iac_dev(:) : alloc_if(.false.) free_if(.false.)) & 
  in(img_i, start_img_crd_dev, end_img_crd_dev, ipairs_idx, &
  my_thread_id,num_threads_offload, my_img_hi, my_img_lo,ee_eval_cnt,full_eval_cnt,x_i, y_i, z_i)  &
  in(lowest_efs_delr2, es_cut2, max_nb_cut2, efs_tbl_dens, eedtbdns_stk, dxdr,del,ntypes_stk, & 
  common_tran, cutoffs_equal, atm_local,i,min_local,max_local) &
  out(eedvir_stk_dev,eed_stk_dev,evdw_stk_dev,ehb_stk_dev, vxx_dev,  &
  vyy_dev, vzz_dev, j, offset_atm)

  eedvir_stk_dev = 0.d0
  eed_stk_dev = 0.d0
  evdw_stk_dev = 0.d0
  ehb_stk_dev = 0.d0

  vxx_dev = 0.d0
  !vxy_dev = 0.d0
  !vxz_dev = 0.d0
  vyy_dev = 0.d0
  !vyz_dev = 0.d0
  vzz_dev = 0.d0
	!$OMP  parallel  
        !$omp DO
  	do i=min_local*3-2,max_local*3
    		img_frc_dev(i,1)=0.d0
  	enddo
	!$omp end do nowait
        !DIR$ VECTOR NONTEMPORAL
	!$OMP DO
	do i=0,num_threads_offload-1
	!dir$ simd
  	do j=min_local*3-2,max_local*3
    		img_frc_thread_dev(i*3*atm_local+j,1)=0.d0
  	enddo
	enddo
        !$omp end parallel

!$omp parallel default(shared) private(img_i, common_tran, ipairs_idx, ee_eval_cnt, &
!$omp&  full_eval_cnt, x_tran_dev, x_i, y_i, z_i, i, my_thread_id) &
!$omp&  reduction(+: eedvir_stk_dev) &
!$omp&  reduction(+: eed_stk_dev, evdw_stk_dev, ehb_stk_dev) &
!$omp&  reduction(+: vxx_dev, vyy_dev, vzz_dev) 

   my_thread_id = omp_get_thread_num()

!$omp do schedule(static)
    !do img_i = my_img_lo, my_img_hi
    do img_i = start_img_crd_dev, end_img_crd_dev

	ipairs_idx = start_indexes_dev(img_i - start_img_crd_dev  + 1)
	!ipairs_idx = start_indexes_dev(img_i - my_img_lo  + 1)

      ! Common translation (ie. no translation) flag is packed at
      ! the front of each sublist followed by the count(s) of sublist
      ! image pair entries.
      common_tran = thread_ipairs_dev(ipairs_idx)
      ipairs_idx = ipairs_idx + 1

      ! Electrostatic evaluation-only count followed by
      ! full evaluation count packed at the front of each pair sublist.

      ee_eval_cnt = thread_ipairs_dev(ipairs_idx)
      full_eval_cnt = thread_ipairs_dev(ipairs_idx + 1)
      ipairs_idx = ipairs_idx + 2
      if (ee_eval_cnt + full_eval_cnt .gt. 0) then

        x_i = img_crd_dev(1, img_i)
        y_i = img_crd_dev(2, img_i)
        z_i = img_crd_dev(3, img_i)

        if (common_tran .eq. 0) then
          ! We need all the translation vectors:
          do i = 0, 17
            x_tran_dev(1, i) = tranvec_dev(1, i) - x_i
            x_tran_dev(2, i) = tranvec_dev(2, i) - y_i
            x_tran_dev(3, i) = tranvec_dev(3, i) - z_i
          end do
        else
          ! Just put the x,y,z values in the middle cell
          x_tran_dev(1, 13) = - x_i
          x_tran_dev(2, 13) = - y_i
          x_tran_dev(3, 13) = - z_i
        end if

        ! We always need virials from this routine if we need energies,
        ! because the virials are used in estimating the pme error.
        if (cutoffs_equal) then
        call pairs_calc_efv_offload(img_frc_dev, img_crd_dev, img_qterm_dev, efs_tbl_dev, &
                              eed_cub_dev, typ_ico_dev, thread_ipairs_dev(ipairs_idx), &
                              gbl_img_iac_dev, gbl_cn1_dev, gbl_cn2_dev, x_tran_dev, img_i,ee_eval_cnt, full_eval_cnt,common_tran,evdw_stk_dev,eed_stk_dev, eedvir_stk_dev, ehb_stk_dev, vxx_dev, vyy_dev, vzz_dev,img_frc_thread_dev,my_thread_id, atm_local, es_cut2, max_nb_cut2,eedtbdns_stk)

        else
          call pairs_calc_efv_2cut_offload(img_frc_dev, img_crd_dev, img_qterm_dev, efs_tbl_dev, &
                                   eed_cub_dev, typ_ico_dev, thread_ipairs_dev(ipairs_idx), &
                                   gbl_img_iac_dev, gbl_cn1_dev, gbl_cn2_dev, x_tran_dev, img_i,ee_eval_cnt, full_eval_cnt,common_tran,evdw_stk_dev,eed_stk_dev, eedvir_stk_dev, ehb_stk_dev, vxx_dev, vyy_dev, vzz_dev,img_frc_thread_dev,my_thread_id, atm_local, es_cut2, max_nb_cut2,eedtbdns_stk)
        end if
       ! ipairs_idx = ipairs_idx + ee_eval_cnt + full_eval_cnt !for OpenMP this is taken care of befor the loop

      end if
    end do
!$omp end parallel
#if 1
	do j=0, num_threads_offload-1
		offset_atm = j*3*atm_local
!		!dir$ vector nontemporal
		!$omp parallel do
		do i=min_local*3 - 2, max_local*3 
			img_frc_dev(i,1) = img_frc_dev(i,1) + img_frc_thread_dev(i+offset_atm ,1)
		end do
	end do
#endif

!dir$ end offload 
end if !mytaskid

!2nd half on host if need energies
  ipairs_idx = start_ipairs_host
    do img_i = start_img_crd_host, end_img_crd_host

      ! Common translation (ie. no translation) flag is packed at
      ! the front of each sublist followed by the count(s) of sublist
      ! image pair entries.

      common_tran = ipairs(ipairs_idx)
      ipairs_idx = ipairs_idx + 1
    
      ! Electrostatic evaluation-only count followed by
      ! full evaluation count packed at the front of each pair sublist.

      ee_eval_cnt = ipairs(ipairs_idx)
      full_eval_cnt = ipairs(ipairs_idx + 1)
      ipairs_idx = ipairs_idx + 2

      if (ee_eval_cnt + full_eval_cnt .gt. 0) then

        x_i = img_crd(1, img_i)
        y_i = img_crd(2, img_i)
        z_i = img_crd(3, img_i)

        if (common_tran .eq. 0) then
          ! We need all the translation vectors:
          do i = 0, 17
            x_tran(1, i) = tranvec(1, i) - x_i
            x_tran(2, i) = tranvec(2, i) - y_i
            x_tran(3, i) = tranvec(3, i) - z_i
          end do
        else
          ! Just put the x,y,z values in the middle cell
          x_tran(1, 13) = - x_i
          x_tran(2, 13) = - y_i
          x_tran(3, 13) = - z_i
        end if

        ! We always need virials from this routine if we need energies,
        ! because the virials are used in estimating the pme error.

        if (cutoffs_equal) then
          call pairs_calc_efv(img_frc, img_crd, img_qterm, efs_tbl, &
                              eed_cub, typ_ico, ipairs(ipairs_idx), &
                              gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
        else
          call pairs_calc_efv_2cut(img_frc, img_crd, img_qterm, efs_tbl, &
                                   eed_cub, typ_ico, ipairs(ipairs_idx), &
                                   gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
        end if

        ipairs_idx = ipairs_idx + ee_eval_cnt + full_eval_cnt !for OpenMP this is taken care of befor the loop

      end if

    end do

if (offload_tasks(mytaskid)) then 
!now make sure the offload compute is complete before the final reduction
!dir$ offload begin mandatory target(mic:0) wait(sig1)
!dir$ end offload 

!Finally Reducing the frc from host and device
	!dir$ simd
	do j=1, atm_cnt
		img_frc(1,j) = img_frc(1,j) +  img_frc_dev(1,j)
		img_frc(2,j) = img_frc(2,j) +  img_frc_dev(2,j)
		img_frc(3,j) = img_frc(3,j) +  img_frc_dev(3,j)
	end do

  eedvir_stk = eedvir_stk + eedvir_stk_dev
  eed_stk = eed_stk + eed_stk_dev
  evdw_stk = evdw_stk + evdw_stk_dev 
  ehb_stk = ehb_stk + ehb_stk_dev 

  vxx = vxx + vxx_dev
!  vxy = vxy + vxy_dev
!  vxz = vxz + vxz_dev 
  vyy = vyy + vyy_dev
!  vyz = vyz + vyz_dev
  vzz = vzz + vzz_dev
end if 

  else !! do not need energies


!1st part on MIC if do not need energies
if (offload_tasks(mytaskid)) then 
!dir$ offload begin mandatory signal(sig2) target (mic:0)  & !changed

  out(img_frc_dev(:,min_local:max_local): alloc_if(.false.) free_if(.false.))  &
  nocopy(img_frc_thread_dev(:,:): alloc_if(.false.) free_if(.false.))  &
  in(img_crd(:, min_local:max_local):  alloc_if(.false.) free_if(.false.) into(img_crd_dev(:,min_local:max_local))) &
  nocopy(thread_ipairs_dev(0:ipairs_maxsize_dev): alloc_if(.false.) free_if(.false.)) &
  nocopy(start_indexes_dev(:): alloc_if(.false.) free_if(.false.)) &
  nocopy(x_tran_dev(:,: ): alloc_if(.false.) free_if(.false.)) &
  in(tranvec(:,0:17): alloc_if(.false.) free_if(.false.) into(tranvec_dev)) & ! may preloaded
  nocopy(img_qterm_dev(:) : alloc_if(.false.) free_if(.false.)) & ! may preloaded 
!  in(eed_cub(1:mxeedtab) : alloc_if(.false.) free_if(.false.) into(eed_cub_dev)) &
  nocopy(eed_cub_dev(:) : alloc_if(.false.) free_if(.false.)) &
  nocopy(gbl_cn1_dev(:) : alloc_if(.false.) free_if(.false.)) &
  nocopy(gbl_cn2_dev(:) : alloc_if(.false.) free_if(.false.)) & 
  nocopy(fs_tbl_dev(:) : alloc_if(.false.) free_if(.false.)) &
  nocopy(typ_ico_dev(:) : alloc_if(.false.) free_if(.false.)) &
  nocopy(gbl_img_iac_dev(:) : alloc_if(.false.) free_if(.false.)) & 
  in(img_i, start_img_crd_dev, end_img_crd_dev, ipairs_idx, &
  my_thread_id,num_threads_offload, my_img_hi, my_img_lo,ee_eval_cnt,full_eval_cnt,x_i, y_i, z_i,min_local, max_local)  &
  in(lowest_efs_delr2, es_cut2, max_nb_cut2, efs_tbl_dens, eedtbdns_stk, dxdr,del,ntypes_stk, & 
  common_tran, cutoffs_equal, atm_local,i) &
  out(eedvir_stk_dev, vxx_dev, &
  vyy_dev, vzz_dev, j, offset_atm) 
  eedvir_stk_dev = 0.d0
  vxx_dev = 0.d0
!  vxy_dev = 0.d0
!  vxz_dev = 0.d0
  vyy_dev = 0.d0
!  vyz_dev = 0.d0
  vzz_dev = 0.d0
	!$OMP  parallel  
        !$omp DO
  	do i=min_local*3-2,max_local*3
    		img_frc_dev(i,1)=0.d0
  	enddo
	!$omp end do nowait
        !DIR$ VECTOR NONTEMPORAL
	!$OMP DO
	do i=0,num_threads_offload-1
	!dir$ simd
  	do j=min_local*3-2,max_local*3
    		img_frc_thread_dev(i*3*atm_local+j,1)=0.d0
  	enddo
	enddo
        !$omp end parallel

!$omp parallel default(shared) private(img_i, common_tran, ipairs_idx, ee_eval_cnt, &
!$omp&  full_eval_cnt, x_tran_dev, x_i, y_i, z_i, i, my_thread_id) &
!$omp&  reduction(+: eedvir_stk_dev) &
!$omp&  reduction(+: vxx_dev, vyy_dev, vzz_dev) 

    my_thread_id = omp_get_thread_num()
    !$omp do schedule(static)
    do img_i = start_img_crd_dev, end_img_crd_dev
	/* The following is the OpenMP to work*/
	ipairs_idx = start_indexes_dev(img_i - start_img_crd_dev + 1)

      ! Common translation (ie. no translation) flag is packed at
      ! the front of each sublist followed by the count(s) of sublist
      ! image pair entries.

      common_tran = thread_ipairs_dev(ipairs_idx)
      ipairs_idx = ipairs_idx + 1
    
      ! Electrostatic evaluation-only count followed by
      ! full evaluation count packed at the front of each pair sublist.

      ee_eval_cnt = thread_ipairs_dev(ipairs_idx)
      full_eval_cnt = thread_ipairs_dev(ipairs_idx + 1)
      ipairs_idx = ipairs_idx + 2
      if (ee_eval_cnt + full_eval_cnt .gt. 0) then

        x_i = img_crd_dev(1, img_i)
        y_i = img_crd_dev(2, img_i)
        z_i = img_crd_dev(3, img_i)

        if (common_tran .eq. 0) then
          ! We need all the translation vectors:
          do i = 0, 17
            x_tran_dev(1, i) = tranvec_dev(1, i) - x_i
            x_tran_dev(2, i) = tranvec_dev(2, i) - y_i
            x_tran_dev(3, i) = tranvec_dev(3, i) - z_i
          end do
        else
          ! Just put the x,y,z values in the middle cell
          x_tran_dev(1, 13) = - x_i
          x_tran_dev(2, 13) = - y_i
          x_tran_dev(3, 13) = - z_i
        end if

        if (need_virials) then

          if (cutoffs_equal) then

            call pairs_calc_fv_offload(img_frc_dev, img_crd_dev, img_qterm_dev, fs_tbl_dev, &
                               eed_cub_dev, typ_ico_dev, thread_ipairs_dev(ipairs_idx), &
                               gbl_img_iac_dev, gbl_cn1_dev, gbl_cn2_dev, x_tran_dev, img_i,ee_eval_cnt, full_eval_cnt,common_tran, eedvir_stk_dev, vxx_dev, vyy_dev, vzz_dev,img_frc_thread_dev,my_thread_id, atm_local, es_cut2, max_nb_cut2,eedtbdns_stk)
 
         else
            call pairs_calc_fv_2cut_offload(img_frc_dev, img_crd_dev, img_qterm_dev, fs_tbl_dev, &
                                    eed_cub_dev, typ_ico_dev, thread_ipairs_dev(ipairs_idx), &
                                    gbl_img_iac_dev, gbl_cn1_dev, gbl_cn2_dev, x_tran_dev, img_i,ee_eval_cnt, full_eval_cnt,common_tran,eedvir_stk_dev, vxx_dev, vyy_dev, vzz_dev,img_frc_thread_dev,my_thread_id, atm_local, es_cut2, max_nb_cut2,eedtbdns_stk)
          end if

        else

          if (cutoffs_equal) then

            call pairs_calc_f_offload(img_frc_dev, img_crd_dev, img_qterm_dev, fs_tbl_dev, &
                              eed_cub_dev, typ_ico_dev, thread_ipairs_dev(ipairs_idx), &
                              gbl_img_iac_dev, gbl_cn1_dev, gbl_cn2_dev, x_tran_dev, img_i,ee_eval_cnt, full_eval_cnt,common_tran, eedvir_stk_dev, vxx_dev, vyy_dev, vzz_dev,img_frc_thread_dev,my_thread_id, atm_local, es_cut2, max_nb_cut2,eedtbdns_stk)
          else
            call pairs_calc_f_2cut_offload(img_frc_dev, img_crd_dev, img_qterm_dev, fs_tbl_dev, &
                                   eed_cub_dev, typ_ico_dev, thread_ipairs_dev(ipairs_idx), &
                                   gbl_img_iac_dev, gbl_cn1_dev, gbl_cn2_dev, x_tran_dev, img_i,ee_eval_cnt, full_eval_cnt,common_tran, eedvir_stk_dev, vxx_dev, vyy_dev, vzz_dev,img_frc_thread_dev,my_thread_id, atm_local, es_cut2, max_nb_cut2,eedtbdns_stk)
          end if

        end if

        !ipairs_idx = ipairs_idx + ee_eval_cnt + full_eval_cnt

      end if
    end do
!$omp end parallel

	do j=0, num_threads_offload-1
		offset_atm = j*3*atm_local
!		!dir$ vector nontemporal
		!$omp parallel do
		do i=min_local*3 - 2, max_local*3 
			img_frc_dev(i,1) = img_frc_dev(i,1) + img_frc_thread_dev(i+offset_atm ,1)
		end do
	end do

!dir$ end offload 
end if !mytaskid

!2nd part on host if do not need energies
  ipairs_idx = start_ipairs_host
    !do img_i = my_img_lo, my_img_hi
    do img_i = start_img_crd_host, end_img_crd_host
      ! Common translation (ie. no translation) flag is packed at
      ! the front of each sublist followed by the count(s) of sublist
      ! image pair entries.

      common_tran = ipairs(ipairs_idx)
      ipairs_idx = ipairs_idx + 1
    
      ! Electrostatic evaluation-only count followed by
      ! full evaluation count packed at the front of each pair sublist.

      ee_eval_cnt = ipairs(ipairs_idx)
      full_eval_cnt = ipairs(ipairs_idx + 1)
      ipairs_idx = ipairs_idx + 2

      if (ee_eval_cnt + full_eval_cnt .gt. 0) then

        x_i = img_crd(1, img_i)
        y_i = img_crd(2, img_i)
        z_i = img_crd(3, img_i)

        if (common_tran .eq. 0) then
          ! We need all the translation vectors:
          do i = 0, 17
            x_tran(1, i) = tranvec(1, i) - x_i
            x_tran(2, i) = tranvec(2, i) - y_i
            x_tran(3, i) = tranvec(3, i) - z_i
          end do
        else
          ! Just put the x,y,z values in the middle cell
          x_tran(1, 13) = - x_i
          x_tran(2, 13) = - y_i
          x_tran(3, 13) = - z_i
        end if

        if (need_virials) then

          if (cutoffs_equal) then
            call pairs_calc_fv(img_frc, img_crd, img_qterm, fs_tbl, &
                               eed_cub, typ_ico, ipairs(ipairs_idx), &
                               gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
          else
            call pairs_calc_fv_2cut(img_frc, img_crd, img_qterm, fs_tbl, &
                                    eed_cub, typ_ico, ipairs(ipairs_idx), &
                                    gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
          end if

        else

          if (cutoffs_equal) then
            call pairs_calc_f(img_frc, img_crd, img_qterm, fs_tbl, &
                              eed_cub, typ_ico, ipairs(ipairs_idx), &
                              gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
          else
            call pairs_calc_f_2cut(img_frc, img_crd, img_qterm, fs_tbl, &
                                   eed_cub, typ_ico, ipairs(ipairs_idx), &
                                   gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
          end if

        end if

        ipairs_idx = ipairs_idx + ee_eval_cnt + full_eval_cnt

      end if

    end do

if (offload_tasks(mytaskid)) then 
!now make sure the offload compute is complete before the final reduction
	!dir$ offload_wait mandatory target(mic:0) wait(sig2) !changed
end if
if (offload_tasks(mytaskid)) then 
!Finally Reducing the frc from host and device
	do j=min_img, max_img
		img_frc(1,j) = img_frc(1,j) +  img_frc_dev(1,j)
		img_frc(2,j) = img_frc(2,j) +  img_frc_dev(2,j)
		img_frc(3,j) = img_frc(3,j) +  img_frc_dev(3,j)
	end do
  	eedvir_stk = eedvir_stk + eedvir_stk_dev

  vxx = vxx + vxx_dev
 !vxy = vxy + vxy_dev
 !vxz = vxz + vxz_dev 
  vyy = vyy + vyy_dev
 !vyz = vyz + vyz_dev
  vzz = vzz + vzz_dev

end if

  end if
!end of the force calculation loop for need_energy and dont_need_energy
  ! Save the energies:
  eedvir = eedvir_stk
  eed = eed_stk
  evdw = evdw_stk
  ehb = ehb_stk

  ! Save the virials.

  virial(1, 1) = vxx
!  virial(1, 2) = vxy
!  virial(2, 1) = vxy
!  virial(1, 3) = vxz
!  virial(3, 1) = vxz
  virial(2, 2) = vyy
!  virial(2, 3) = vyz
!  virial(3, 2) = vyz
  virial(3, 3) = vzz

!print *, "get nb energy ends",mytaskid
  return

contains

#define BUILD_PAIRS_CALC_EFV
#define NEED_ENE
#define NEED_VIR
#include "pairs_calc_offload.i"
#undef NEED_VIR
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_EFV

#define BUILD_PAIRS_CALC_FV
#define NEED_VIR
#include "pairs_calc_offload.i"
#undef NEED_VIR
#undef BUILD_PAIRS_CALC_FV

#define BUILD_PAIRS_CALC_F
#include "pairs_calc_offload.i"
#undef BUILD_PAIRS_CALC_F

#define BUILD_PAIRS_CALC_EFV
#define NEED_ENE
#define NEED_VIR
#include "pairs_calc.i"
#undef NEED_VIR
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_EFV

#define BUILD_PAIRS_CALC_FV
#define NEED_VIR
#include "pairs_calc.i"
#undef NEED_VIR
#undef BUILD_PAIRS_CALC_FV

#define BUILD_PAIRS_CALC_F
#include "pairs_calc.i"
#undef BUILD_PAIRS_CALC_F

#define BUILD_PAIRS_CALC_EFV
#define NEED_ENE
#define NEED_VIR
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#include "pairs_calc_offload.i"
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef NEED_VIR
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_EFV

#define BUILD_PAIRS_CALC_EFV
#define NEED_ENE
#define NEED_VIR
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#include "pairs_calc.i"
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef NEED_VIR
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_EFV

#define BUILD_PAIRS_CALC_FV
#define NEED_VIR
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#include "pairs_calc_offload.i"
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_FV

#define BUILD_PAIRS_CALC_FV
#define NEED_VIR
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#include "pairs_calc.i"
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_FV

#define BUILD_PAIRS_CALC_F
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#include "pairs_calc_offload.i"
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef BUILD_PAIRS_CALC_F

#define BUILD_PAIRS_CALC_F
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#include "pairs_calc.i"
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef BUILD_PAIRS_CALC_F

end subroutine get_nb_energy_offload
#endif /* MIC_offload */
