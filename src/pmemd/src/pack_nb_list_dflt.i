!*******************************************************************************
!
! Internal Subroutine:  pack_nb_list
!
! Description: See get_nb_list comments above.
!              
!*******************************************************************************

subroutine pack_nb_list(ee_eval_cnt, img_j_ee_eval, &
                        full_eval_cnt, img_j_full_eval, ipairs, & 
                        num_packed)
  implicit none

! Formal arguments:

  integer               :: ee_eval_cnt
  integer               :: img_j_ee_eval(*)
  integer               :: full_eval_cnt
  integer               :: img_j_full_eval(*)
  integer               :: ipairs(*)
  integer               :: num_packed

! Local variables:

  integer               :: i

! Check for enough room on the pairs list.  We also allow space for the
! two pairlist counters at the front of the list, and for one blank integer at
! the end of the list because some pairs calc implementations may prefetch one
! integer and could thus run off the end of the list. 
! We also allow space for one translation flag at the front of the list.

  if (num_packed + ee_eval_cnt + full_eval_cnt + 4 .le. ipairs_maxsize) then

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

    ifail = 1

  end if

  return

end subroutine pack_nb_list

!*******************************************************************************
!
! Internal Subroutine:  pack_nb_list_skip_belly_pairs
!
! Description: See get_nb_list comments above.
!              
!*******************************************************************************

subroutine pack_nb_list_skip_belly_pairs(ee_eval_cnt, img_j_ee_eval, &
                                         full_eval_cnt, img_j_full_eval, &
                                         ipairs, igroup, img_atm_map, &
                                       num_packed)

  implicit none

! Formal arguments:

  integer               :: ee_eval_cnt
  integer               :: img_j_ee_eval(*)
  integer               :: full_eval_cnt
  integer               :: img_j_full_eval(*)
  integer               :: ipairs(*)
  integer               :: igroup(*)
  integer               :: img_atm_map(*)
  integer               :: num_packed

! Local variables:

  integer               :: i
  integer               :: img_j
  integer               :: new_ee_eval_cnt
  integer               :: new_full_eval_cnt
  integer, parameter    :: mask27 = int(Z"07FFFFFF")


  if (igroup(atm_i) .eq. 0) then

    new_ee_eval_cnt = 0
    new_full_eval_cnt = 0

    if (common_tran) then

      do i = 1, ee_eval_cnt
        if (igroup(img_atm_map(img_j_ee_eval(i))) .ne. 0) then
          new_ee_eval_cnt = new_ee_eval_cnt + 1
          img_j_ee_eval(new_ee_eval_cnt) = img_j_ee_eval(i)
        end if
      end do

      do i = 1, full_eval_cnt
        if (igroup(img_atm_map(img_j_full_eval(i))) .ne. 0) then
          new_full_eval_cnt = new_full_eval_cnt + 1
          img_j_full_eval(new_full_eval_cnt) = img_j_full_eval(i)
        end if
      end do

    else

      do i = 1, ee_eval_cnt
        img_j = iand(img_j_ee_eval(i), mask27)
        if (igroup(img_atm_map(img_j)) .ne. 0) then
          new_ee_eval_cnt = new_ee_eval_cnt + 1
          img_j_ee_eval(new_ee_eval_cnt) = img_j_ee_eval(i)
        end if
      end do

      do i = 1, full_eval_cnt
        img_j = iand(img_j_full_eval(i), mask27)
        if (igroup(img_atm_map(img_j)) .ne. 0) then
          new_full_eval_cnt = new_full_eval_cnt + 1
          img_j_full_eval(new_full_eval_cnt) = img_j_full_eval(i)
        end if
      end do

    end if

    ee_eval_cnt = new_ee_eval_cnt
    full_eval_cnt = new_full_eval_cnt

  end if

  call pack_nb_list(ee_eval_cnt, img_j_ee_eval, &
                    full_eval_cnt, img_j_full_eval, &
                    gbl_ipairs, num_packed)

  return

end subroutine pack_nb_list_skip_belly_pairs

!*******************************************************************************
!
! Internal Subroutine:  ti_pack_nb_list
!
! Description: See get_nb_list comments above. Modified for TI to seperate
! calculation into common atms, interactions with TI atoms, and TI-TI 
! interactions. Any interaction between V0/V1 TI atoms should NOT be included.
!              
!*******************************************************************************

subroutine ti_pack_nb_list(ee_eval_cnt, img_j_ee_eval, &
                        full_eval_cnt, img_j_full_eval, &
                        ipairs, num_packed)

  implicit none

! Formal arguments:
  integer               :: ee_eval_cnt(:)
  integer               :: img_j_ee_eval(:,:)
  integer               :: full_eval_cnt(:)
  integer               :: img_j_full_eval(:,:)
  
  integer               :: ipairs(*)
  integer               :: num_packed

! Local variables:
  integer               :: i, j
  integer               :: header_size
! Check for enough room on the pairs list.  We also allow space for the
! pairlist counters at the front of the list, and for one blank integer at
! the end of the list because some pairs calc implementations may prefetch one
! integer and could thus run off the end of the list. 
! We also allow space for one translation flag at the front of the list.

  header_size = 1 + nb_list_cnt * 2 ! + 1 for common_tran

  if (num_packed + header_size + &
      sum(ee_eval_cnt) + sum(full_eval_cnt) .le. ipairs_maxsize) then

    if (common_tran) then
      ipairs(num_packed + 1) = 1
    else
      ipairs(num_packed + 1) = 0
    end if

    j = num_packed + header_size + 1
    do i = 1, nb_list_cnt   
      ipairs(num_packed + 2*(i-1) + 2) = ee_eval_cnt(i)    
      ipairs(num_packed + 2*(i-1) + 3) = full_eval_cnt(i)  
      
      if (ee_eval_cnt(i) .gt. 0) then
        ipairs(j:j + ee_eval_cnt(i) - 1) = img_j_ee_eval(i, 1:ee_eval_cnt(i))        
        j = j + ee_eval_cnt(i)
      end if
      if (full_eval_cnt(i) .gt. 0) then
        ipairs(j:j + full_eval_cnt(i) - 1) = img_j_full_eval(i, 1:full_eval_cnt(i))        
        j = j + full_eval_cnt(i)
      end if      
    end do     

    num_packed = num_packed + header_size + &
                 sum(ee_eval_cnt) + sum(full_eval_cnt)
  else

    ifail = 1

  end if

  return

end subroutine ti_pack_nb_list

!*******************************************************************************
!
! Internal Subroutine:  pack_nb_list_skip_ti_pairs
!
! Description: See get_nb_list comments above. Modified to prevent cross terms
! between the non-interacting TI regions.
!              
!*******************************************************************************

subroutine pack_nb_list_skip_ti_pairs(ee_eval_cnt, img_j_ee_eval, &
                                         full_eval_cnt, img_j_full_eval, &
                                         ipairs, img_atm_map, &
                                         num_packed)
  use ti_mod
  implicit none

! Formal arguments:

  integer               :: ee_eval_cnt
  integer               :: img_j_ee_eval(*)
  integer               :: full_eval_cnt
  integer               :: img_j_full_eval(*)
  integer               :: ipairs(*)
  integer               :: img_atm_map(*)
  integer               :: num_packed
  
! Local variables: 
  integer               :: i
  integer               :: img_j
  integer               :: atm_j
  integer               :: ti_ee_eval_cnt(nb_list_cnt)
  integer               :: ti_full_eval_cnt(nb_list_cnt)
  integer               :: ti_img_j_ee_eval(nb_list_cnt, ee_eval_cnt) !max size is ee_eval_cnt
  integer               :: ti_img_j_full_eval(nb_list_cnt, full_eval_cnt)
  
  ! For TI:
  integer               :: cur_nb_list
  integer               :: sc_lcl_sum
  integer, parameter    :: mask27 = int(Z"07FFFFFF")
  
  ti_ee_eval_cnt(:) = 0
  ti_full_eval_cnt(:) = 0
  ti_img_j_ee_eval(:,:) = 0
  ti_img_j_full_eval(:,:) = 0
  
  if (common_tran) then

    do i = 1, ee_eval_cnt
      atm_j = img_atm_map(img_j_ee_eval(i))
      sc_lcl_sum = ti_lst(1,atm_i) + ti_lst(2,atm_i) + &
                   ti_lst(1,atm_j) + ti_lst(2,atm_j)

      cur_nb_list = nb_list_null      
      if (sc_lcl_sum .eq. 0) then ! 0 TI atoms
        cur_nb_list = nb_list_common
      else if (sc_lcl_sum .eq. 1) then ! 1 TI atom
        if (ti_mode .eq. 1 .or. (ti_sc_lst(atm_i) .eq. 0 .and. ti_sc_lst(atm_j) .eq. 0)) then    
          if (ti_lst(1,atm_i)+ti_lst(1,atm_j) .ge. 1) then
            cur_nb_list = nb_list_linear_V0
          else
            cur_nb_list = nb_list_linear_V1
          end if    
        else if (ti_sc_lst(atm_i)+ti_sc_lst(atm_j) .eq. 2) then
          if (ti_lst(1,atm_i)+ti_lst(1,atm_j) .ge. 1) then
            cur_nb_list = nb_list_sc_common_V0
          else
            cur_nb_list = nb_list_sc_common_V1
          end if           
        end if
      else if (ti_lst(1,atm_i)+ti_lst(2,atm_j) .eq. 1) then ! 2 TI atoms, no cross term V0 to V1
        if (ti_sc_lst(atm_i)+ti_sc_lst(atm_j) .eq. 4) then ! both atoms softcore
          if (ti_lst(1,atm_i)+ti_lst(1,atm_j) .ge. 1) then
            cur_nb_list = nb_list_sc_sc_V0
          else
            cur_nb_list = nb_list_sc_sc_V1
          end if      
        else if (ti_sc_lst(atm_i)+ti_sc_lst(atm_j) .eq. 2) then ! one atom softcore
          if (ti_lst(1,atm_i)+ti_lst(1,atm_j) .ge. 1) then
            cur_nb_list = nb_list_sc_common_V0
          else
            cur_nb_list = nb_list_sc_common_V1
          end if 
        else ! only linearly scaled atoms
          if (ti_lst(1,atm_i)+ti_lst(1,atm_j) .ge. 1) then
            cur_nb_list = nb_list_linear_V0
          else
            cur_nb_list = nb_list_linear_V1
          end if 
        end if      
      end if

      if (cur_nb_list .ne. nb_list_null) then
        ti_ee_eval_cnt(cur_nb_list) = ti_ee_eval_cnt(cur_nb_list) + 1
        ti_img_j_ee_eval(cur_nb_list, ti_ee_eval_cnt(cur_nb_list)) = img_j_ee_eval(i)
      end if
    end do

    do i = 1, full_eval_cnt
      atm_j = img_atm_map(img_j_full_eval(i))
      sc_lcl_sum = ti_lst(1,atm_i) + ti_lst(2,atm_i) + &
                   ti_lst(1,atm_j) + ti_lst(2,atm_j)

      cur_nb_list = nb_list_null      
      if (sc_lcl_sum .eq. 0) then ! 0 TI atoms
        cur_nb_list = nb_list_common
      else if (sc_lcl_sum .eq. 1) then ! 1 TI atom
        if (ti_mode .eq. 1 .or. (ti_sc_lst(atm_i) .eq. 0 .and. ti_sc_lst(atm_j) .eq. 0)) then    
          if (ti_lst(1,atm_i)+ti_lst(1,atm_j) .ge. 1) then
            cur_nb_list = nb_list_linear_V0
          else
            cur_nb_list = nb_list_linear_V1
          end if    
        else if (ti_sc_lst(atm_i)+ti_sc_lst(atm_j) .eq. 2) then
          if (ti_lst(1,atm_i)+ti_lst(1,atm_j) .ge. 1) then
            cur_nb_list = nb_list_sc_common_V0
          else
            cur_nb_list = nb_list_sc_common_V1
          end if           
        end if
      else if (ti_lst(1,atm_i)+ti_lst(2,atm_j) .eq. 1) then ! 2 TI atoms, no cross term V0 to V1
        if (ti_sc_lst(atm_i)+ti_sc_lst(atm_j) .eq. 4) then ! both atoms softcore
          if (ti_lst(1,atm_i)+ti_lst(1,atm_j) .ge. 1) then
            cur_nb_list = nb_list_sc_sc_V0
          else
            cur_nb_list = nb_list_sc_sc_V1
          end if      
        else if (ti_sc_lst(atm_i)+ti_sc_lst(atm_j) .eq. 2) then ! one atom softcore
          if (ti_lst(1,atm_i)+ti_lst(1,atm_j) .ge. 1) then
            cur_nb_list = nb_list_sc_common_V0
          else
            cur_nb_list = nb_list_sc_common_V1
          end if 
        else ! only linearly scaled atoms
          if (ti_lst(1,atm_i)+ti_lst(1,atm_j) .ge. 1) then
            cur_nb_list = nb_list_linear_V0
          else
            cur_nb_list = nb_list_linear_V1
          end if 
        end if      
      end if
   
      if (cur_nb_list .ne. nb_list_null) then
        ti_full_eval_cnt(cur_nb_list) = ti_full_eval_cnt(cur_nb_list) + 1
        ti_img_j_full_eval(cur_nb_list, ti_full_eval_cnt(cur_nb_list)) = img_j_full_eval(i)
      end if
    end do

  else

    do i = 1, ee_eval_cnt
      img_j = iand(img_j_ee_eval(i), mask27)
      atm_j = img_atm_map(img_j)
      sc_lcl_sum = ti_lst(1,atm_i) + ti_lst(2,atm_i) + &
                   ti_lst(1,atm_j) + ti_lst(2,atm_j)

      cur_nb_list = nb_list_null      
      if (sc_lcl_sum .eq. 0) then ! 0 TI atoms
        cur_nb_list = nb_list_common
      else if (sc_lcl_sum .eq. 1) then ! 1 TI atom
        if (ti_mode .eq. 1 .or. (ti_sc_lst(atm_i) .eq. 0 .and. ti_sc_lst(atm_j) .eq. 0)) then    
          if (ti_lst(1,atm_i)+ti_lst(1,atm_j) .ge. 1) then
            cur_nb_list = nb_list_linear_V0
          else
            cur_nb_list = nb_list_linear_V1
          end if    
        else if (ti_sc_lst(atm_i)+ti_sc_lst(atm_j) .eq. 2) then
          if (ti_lst(1,atm_i)+ti_lst(1,atm_j) .ge. 1) then
            cur_nb_list = nb_list_sc_common_V0
          else
            cur_nb_list = nb_list_sc_common_V1
          end if           
        end if
      else if (ti_lst(1,atm_i)+ti_lst(2,atm_j) .eq. 1) then ! 2 TI atoms, no cross term V0 to V1
        if (ti_sc_lst(atm_i)+ti_sc_lst(atm_j) .eq. 4) then ! both atoms softcore
          if (ti_lst(1,atm_i)+ti_lst(1,atm_j) .ge. 1) then
            cur_nb_list = nb_list_sc_sc_V0
          else
            cur_nb_list = nb_list_sc_sc_V1
          end if      
        else if (ti_sc_lst(atm_i)+ti_sc_lst(atm_j) .eq. 2) then ! one atom softcore
          if (ti_lst(1,atm_i)+ti_lst(1,atm_j) .ge. 1) then
            cur_nb_list = nb_list_sc_common_V0
          else
            cur_nb_list = nb_list_sc_common_V1
          end if 
        else ! only linearly scaled atoms
          if (ti_lst(1,atm_i)+ti_lst(1,atm_j) .ge. 1) then
            cur_nb_list = nb_list_linear_V0
          else
            cur_nb_list = nb_list_linear_V1
          end if 
        end if      
      end if

      if (cur_nb_list .ne. nb_list_null) then      
        ti_ee_eval_cnt(cur_nb_list) = ti_ee_eval_cnt(cur_nb_list) + 1
        ti_img_j_ee_eval(cur_nb_list, ti_ee_eval_cnt(cur_nb_list)) = img_j_ee_eval(i)
      end if
    end do

    do i = 1, full_eval_cnt
      img_j = iand(img_j_full_eval(i), mask27)
      atm_j = img_atm_map(img_j)
      sc_lcl_sum = ti_lst(1,atm_i) + ti_lst(2,atm_i) + &
                   ti_lst(1,atm_j) + ti_lst(2,atm_j)
                   
      cur_nb_list = nb_list_null      
      if (sc_lcl_sum .eq. 0) then ! 0 TI atoms
        cur_nb_list = nb_list_common
      else if (sc_lcl_sum .eq. 1) then ! 1 TI atom
        if (ti_mode .eq. 1 .or. (ti_sc_lst(atm_i) .eq. 0 .and. ti_sc_lst(atm_j) .eq. 0)) then    
          if (ti_lst(1,atm_i)+ti_lst(1,atm_j) .ge. 1) then
            cur_nb_list = nb_list_linear_V0
          else
            cur_nb_list = nb_list_linear_V1
          end if    
        else if (ti_sc_lst(atm_i)+ti_sc_lst(atm_j) .eq. 2) then
          if (ti_lst(1,atm_i)+ti_lst(1,atm_j) .ge. 1) then
            cur_nb_list = nb_list_sc_common_V0
          else
            cur_nb_list = nb_list_sc_common_V1
          end if           
        end if
      else if (ti_lst(1,atm_i)+ti_lst(2,atm_j) .eq. 1) then ! 2 TI atoms, no cross term V0 to V1
        if (ti_sc_lst(atm_i)+ti_sc_lst(atm_j) .eq. 4) then ! both atoms softcore
          if (ti_lst(1,atm_i)+ti_lst(1,atm_j) .ge. 1) then
            cur_nb_list = nb_list_sc_sc_V0
          else
            cur_nb_list = nb_list_sc_sc_V1
          end if      
        else if (ti_sc_lst(atm_i)+ti_sc_lst(atm_j) .eq. 2) then ! one atom softcore
          if (ti_lst(1,atm_i)+ti_lst(1,atm_j) .ge. 1) then
            cur_nb_list = nb_list_sc_common_V0
          else
            cur_nb_list = nb_list_sc_common_V1
          end if 
        else ! only linearly scaled atoms
          if (ti_lst(1,atm_i)+ti_lst(1,atm_j) .ge. 1) then
            cur_nb_list = nb_list_linear_V0
          else
            cur_nb_list = nb_list_linear_V1
          end if 
        end if      
      end if

      if (cur_nb_list .ne. nb_list_null) then
        ti_full_eval_cnt(cur_nb_list) = ti_full_eval_cnt(cur_nb_list) + 1
        ti_img_j_full_eval(cur_nb_list, ti_full_eval_cnt(cur_nb_list)) = img_j_full_eval(i)
      end if
    end do

  end if

  call ti_pack_nb_list(ti_ee_eval_cnt, ti_img_j_ee_eval, &
                       ti_full_eval_cnt, ti_img_j_full_eval, &
                       gbl_ipairs, num_packed)
  
  return

end subroutine pack_nb_list_skip_ti_pairs
