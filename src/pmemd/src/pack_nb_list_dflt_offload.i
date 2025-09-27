
#ifdef MIC_offload
!*******************************************************************************
!
! Internal Subroutine:  pack_nb_list_offload
!
! Description: See get_nb_list comments above.
!              
!*******************************************************************************
!DIR$ ATTRIBUTES OFFLOAD : MIC :: pack_nb_list_offload
subroutine pack_nb_list_offload(ee_eval_cnt, img_j_ee_eval, &
#ifdef MPI
	full_eval_cnt, img_j_full_eval, ipairs, num_packed, common_tran, ipairs_maxsize,ifail_thread, my_thread_id,total_threads)! ipairs_maxsize if per thread if openmp is enabled
#else
	full_eval_cnt, img_j_full_eval, ipairs, num_packed)
#endif
                        

  implicit none

! Formal arguments:

  integer               :: ee_eval_cnt
  integer               :: img_j_ee_eval(*)
  integer               :: full_eval_cnt
  integer               :: img_j_full_eval(*)
  integer               :: ipairs(*)
  integer               :: num_packed
#ifdef MPI
  logical		:: common_tran
  integer		:: ipairs_maxsize
  integer		:: total_threads
  integer		:: ifail_thread(0:total_threads)
  !integer		:: ifail_thread(*)
 integer		:: my_thread_id
#endif
! Local variables:

  integer               :: i

! Check for enough room on the pairs list.  We also allow space for the
! two pairlist counters at the front of the list, and for one blank integer at
! the end of the list because some pairs calc implementations may prefetch one
! integer and could thus run off the end of the list.  If DIRFRC_COMTRANS is 
! defined we also allow space for one translation flag at the front of the list.


  if (num_packed + ee_eval_cnt + full_eval_cnt + 4 .le. ipairs_maxsize-1) then
    if (common_tran) then
      ipairs(num_packed + 1) = 1
    else
      ipairs(num_packed + 1) = 0
    end if

    ipairs(num_packed + 2) = ee_eval_cnt
    ipairs(num_packed + 3) = full_eval_cnt

    if (ee_eval_cnt .gt. 0) then
      i = num_packed + 4
      ipairs(i:i + ee_eval_cnt - 1) = img_j_ee_eval(1:ee_eval_cnt)
    end if

    if (full_eval_cnt .gt. 0) then
      i = num_packed + 4 + ee_eval_cnt
      ipairs(i:i + full_eval_cnt - 1) = img_j_full_eval(1:full_eval_cnt)
    end if

    num_packed = num_packed + ee_eval_cnt + full_eval_cnt + 3

  else
    ifail_thread(my_thread_id) =  1 !all the openmp thread  adds to it
    !ifail =  1

  end if


  return

end subroutine pack_nb_list_offload
#endif
