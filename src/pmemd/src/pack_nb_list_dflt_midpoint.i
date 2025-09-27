
!*******************************************************************************
!
! Internal Subroutine:  pack_nb_list_midpoint
!
! Description: See get_nb_list comments above.
!              
!*******************************************************************************

subroutine pack_nb_list_midpoint(ee_eval_cnt, atm_j_ee_eval, &
                        full_eval_cnt, atm_j_full_eval, ipairs, & 
                        num_packed)
  implicit none

! Formal arguments:

  integer               :: ee_eval_cnt
  integer               :: atm_j_ee_eval(*)
  integer               :: full_eval_cnt
  integer               :: atm_j_full_eval(*)
  integer               :: ipairs(*)
  integer               :: num_packed

! Local variables:

  integer               :: i
  integer               :: dummy_common_tran

! Check for enough room on the pairs list.  We also allow space for the
! two pairlist counters at the front of the list, and for one blank integer at
! the end of the list because some pairs calc implementations may prefetch one
! integer and could thus run off the end of the list. 
! We also allow space for one translation flag at the front of the list.
  
  dummy_common_tran = 1 ! in case we use common_tran, we modify 
  if (num_packed + ee_eval_cnt + full_eval_cnt + 4 .le. proc_ipairs_maxsize) then

    ipairs(num_packed + 1) = dummy_common_tran!may need to modify

    ipairs(num_packed + 2) = ee_eval_cnt
    ipairs(num_packed + 3) = full_eval_cnt

    if (ee_eval_cnt .gt. 0) then
      i = num_packed + 4
      ipairs(i:i + ee_eval_cnt - 1) = atm_j_ee_eval(1:ee_eval_cnt)
    end if

    if (full_eval_cnt .gt. 0) then
      i = num_packed + 4 + ee_eval_cnt
      ipairs(i:i + full_eval_cnt - 1) = atm_j_full_eval(1:full_eval_cnt)
    end if

    num_packed = num_packed + ee_eval_cnt + full_eval_cnt + 3

  else

    ifail = 1

  end if

  return

end subroutine pack_nb_list_midpoint
