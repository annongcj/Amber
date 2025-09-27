#include "mpi_err_check.i"

!*******************************************************************************
!
! Module: parallel_processor_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module parallel_processor_mod

#ifdef MPI

  !use gbl_datatypes_mod
  use parallel_dat_mod
  implicit none

contains


!*******************************************************************************
!
! Subroutine:  comm_ensure_ghost_alloc
!
! Description: get the ghost atoms count  and expands arrays as needed 
!                before the actual data receive
!
!*******************************************************************************
subroutine comm_ensure_ghost_alloc(side)
  use processor_mod
  use ensure_alloc_mod
  ! Parameters
  integer side
  
  ! local variables
  integer        :: recv_req(2), send_req(2),recvd_idx
  integer       :: irecv_stat(mpi_status_size)
  integer       :: isend_stat(mpi_status_size, 2)
  integer i, j, k, s, atm_idx, direction, sides(2)
  integer, save :: call_cnt=0
  integer alloc_size

  ! with values 1 to 6 in the order: x-, x+, y-, y+, z-, z+
  sides(1)=(side-1)*2+1
  sides(2)=sides(1)+1
  !post MPI receive
  do s =1,2
    direction=sides(s)
    call mpi_irecv(recv_buf_cnt(direction), 1, MPI_INTEGER, proc_neighbor_rank(direction),&
                    3-s, pmemd_comm, recv_req(s), err_code_mpi)
  enddo
  do s =1,2
    direction=sides(s)
    !post MPI send
    call mpi_isend(send_buf_cnt(direction), 1, MPI_INTEGER, proc_neighbor_rank(direction),&
                   s, pmemd_comm, send_req(s), err_code_mpi)
  enddo
  !wait MPI receive and count
  do s =1,2
    call mpi_waitany(2, recv_req, recvd_idx, irecv_stat, err_code_mpi)
  end do

  !wait MPI send
  call mpi_waitall(2, send_req, isend_stat, err_code_mpi)


  ! update the allocation size if needed
  alloc_size = proc_num_atms_min_bound+sum(recv_buf_cnt(1:2*side))
  if(alloc_size .gt. proc_num_atms_min_bound+proc_ghost_num_atms) then  
    call ensure_alloc_dble3(proc_atm_crd, alloc_size, .true.)
    call ensure_alloc_int2d(send_buf_indices, alloc_size, 6, .true.)
    call ensure_alloc_int2d(mpi_send_buf_int, alloc_size, 2)
    call ensure_alloc_int2d(mpi_recv_buf_int, alloc_size, 2)
    call ensure_alloc_dble2d(mpi_send_buf_dble, 3*alloc_size, 2)
    call ensure_alloc_dble2d(mpi_recv_buf_dble, 3*alloc_size, 2)
    proc_atm_alloc_size = size(send_buf_indices)
  endif

#if 0
  call_cnt=call_cnt+1
  do s =1,2
    direction=sides(s)
    write(0, "(I3,':', I1,'-', I0.3,' With neighbor ',I1,' rank:',I3,' Sent ',I,' Recvd ',I)") &
                  call_cnt/3+1, side, mytaskid,direction, proc_neighbor_rank(direction), send_buf_cnt(direction), recv_buf_cnt(direction)
  end do
#endif
end subroutine comm_ensure_ghost_alloc



!*******************************************************************************
!
! Subroutine:  comm_3dimensions_1int
!
! Description: Performs the communication in the three dimensions of the domain
!              for a 1 components Integer array
!
!*******************************************************************************
subroutine comm_3dimensions_1int(int1_array)
  use processor_mod
  use pmemd_lib_mod,     only : mexit

  ! Parameters
  integer int1_array(proc_atm_alloc_size)
  
  ! local variables
  integer mpi_recv_cnt(2), mpi_send_cnt(2)! 2 elements for two directions
  integer        :: recv_req(2), send_req(2),recvd_idx
  integer       :: irecv_stat(mpi_status_size)
  integer       :: isend_stat(mpi_status_size, 2)
  integer i, j, k, s, atm_idx, direction, sides(2), ghost_cnt, dimensions

  ghost_cnt=0

  do dimensions =1,3
    sides(1)=(dimensions-1)*2+1
    sides(2)=sides(1)+1
    !post MPI receive
    do s =1,2
      direction=sides(s)
      mpi_recv_cnt(s) = recv_buf_cnt(direction)
      call mpi_irecv(mpi_recv_buf_int(:,s), mpi_recv_cnt(s), mpi_integer, &
             proc_neighbor_rank(direction), 3-s, pmemd_comm, recv_req(s), err_code_mpi)
    enddo
    do s =1,2
      direction=sides(s)
      !pack the send data
      do i =1, send_buf_cnt(direction)
        atm_idx=send_buf_indices(i,direction)
        mpi_send_buf_int(i,s) = int1_array(atm_idx)
      enddo
      !post MPI send
      mpi_send_cnt(s) = send_buf_cnt(direction)
      call mpi_isend(mpi_send_buf_int(:,s), mpi_send_cnt(s), mpi_integer, &
             proc_neighbor_rank(direction), s, pmemd_comm, send_req(s), err_code_mpi)
    enddo
    !wait MPI receive, count, and unpack
    do s =1,2
      call mpi_waitany(2, recv_req, recvd_idx, irecv_stat, err_code_mpi)
      call mpi_get_count(irecv_stat, mpi_integer, mpi_recv_cnt(recvd_idx), err_code_mpi)
      direction=sides(recvd_idx)
      if(recv_buf_cnt(direction) .ne. mpi_recv_cnt(recvd_idx)) then
           write(0,"('ERROR: 1int recv count does not match the expected. Expected: ',I9,' Recvd: ', I9)")&
                 recv_buf_cnt(direction),  mpi_recv_cnt(recvd_idx) 
           call mexit(6,1)
      end if
    end do
    !unpack
    do s =1,2
      do j= 1,mpi_recv_cnt(s)
        ghost_cnt = ghost_cnt+1
        int1_array(proc_num_atms+ghost_cnt) = mpi_recv_buf_int(j,s)
      enddo
    end do
    !wait MPI send
    call mpi_waitall(2, send_req, isend_stat, err_code_mpi)
#if 0
    do s =1,2
      direction=sides(s)
      write(0, "(I1,'-', I0.3,' With neighbor ',I1,' rank:',I3,' Sent ',I9,' Recvd ',I9)") &
                    dimensions, mytaskid,direction, proc_neighbor_rank(direction), mpi_send_cnt(s),  mpi_recv_cnt(s)
    end do
#endif
  end do
end subroutine comm_3dimensions_1int

!*******************************************************************************
!
! Subroutine:  comm_3dimensions_1dbl
!
! Description: Performs the communication in the three dimensions of the domain
!              for a 1D Double Precision array
!
!*******************************************************************************
! Aggregate the 3D comm calls in one call
subroutine comm_3dimensions_1dbl(dbl1_array)
  use processor_mod
  use pmemd_lib_mod,     only : mexit ! Parameters
  double precision dbl1_array(proc_atm_alloc_size)
  
  ! local variables
  integer          mpi_recv_cnt(2), mpi_send_cnt(2)! 2 elements for two directions
  integer       :: recv_req(2), send_req(2),recvd_idx
  integer       :: irecv_stat(mpi_status_size)
  integer       :: isend_stat(mpi_status_size, 2)
  integer          i, j, k, s, atm_idx, direction, sides(2), ghost_cnt, dimensions

  ghost_cnt=0

  do dimensions=1,3
    sides(1)=(dimensions-1)*2+1
    sides(2)=sides(1)+1
    !post MPI receive
    do s =1,2
      direction=sides(s)
      mpi_recv_cnt(s) = recv_buf_cnt(direction)
      call mpi_irecv(mpi_recv_buf_dble(:,s), mpi_recv_cnt(s), mpi_double, proc_neighbor_rank(direction),&
                                3-s, pmemd_comm, recv_req(s), err_code_mpi)
    enddo
    do s =1,2
      direction=sides(s)
      !pack the send data
      do i =1, send_buf_cnt(direction)
        atm_idx=send_buf_indices(i,direction)
        mpi_send_buf_dble(i,s) = dbl1_array(atm_idx)
      enddo
      !post MPI send
      mpi_send_cnt(s) = send_buf_cnt(direction)
      call mpi_isend(mpi_send_buf_dble(:,s), mpi_send_cnt(s), mpi_double, proc_neighbor_rank(direction),&
                                 s, pmemd_comm, send_req(s), err_code_mpi)
    enddo
    !wait MPI receive, count, and unpack
    do s =1,2
      call mpi_waitany(2, recv_req, recvd_idx, irecv_stat, err_code_mpi)
      call mpi_get_count(irecv_stat, mpi_double, mpi_recv_cnt(recvd_idx), err_code_mpi)
      direction=sides(recvd_idx)
      if(recv_buf_cnt(direction) .ne. mpi_recv_cnt(recvd_idx)) then
           write(0,"('ERROR: 1dble recv count does not match the expected. Expected: ',I9,' Recvd: ', I9)")&
                      recv_buf_cnt(direction),  mpi_recv_cnt(recvd_idx) 
           call mexit(6,1)
      end if
    end do
    !unpack
    do s =1,2
      do j= 1,mpi_recv_cnt(s)
        ghost_cnt = ghost_cnt+1
        dbl1_array(proc_num_atms+ghost_cnt) = mpi_recv_buf_dble(j,s)
      enddo
    end do
    !wait MPI send
    call mpi_waitall(2, send_req, isend_stat, err_code_mpi)
#if 0
    do s =1,2
      direction=sides(s)
      write(0, "(I1,'-', I0.3,' With neighbor ',I1,' rank:',I3,' Sent ',I9,' Recvd ',I9)") &
                    dimensions, mytaskid,direction, proc_neighbor_rank(direction), mpi_send_cnt(s),  mpi_recv_cnt(s)
    end do
#endif
  end do
end subroutine comm_3dimensions_1dbl





!*******************************************************************************
!
! Subroutine:  comm_frc
!
! Description: Performs the force arrays updates across MPI tasks
!
!*******************************************************************************
subroutine comm_frc(atm_frc)
  use processor_mod
  use pmemd_lib_mod,     only : mexit
  ! Parameters
  double precision atm_frc(3,proc_atm_alloc_size)
  
  ! local variables
  integer mpi_recv_cnt(2), mpi_send_cnt(2)! 2 elements for two directions
  integer        :: recv_req(2), send_req(2),recvd_idx
  integer       :: irecv_stat(mpi_status_size)
  integer       :: isend_stat(mpi_status_size, 2)
  integer          ind, i, j, k, s, atm_idx, direction, sides(2)
  integer side

  ! process the dimensions in the order Z,Y, then X
  do side=3,1,-1

    ! with values 1 to 6 in the order: x-, x+, y-, y+, z-, z+
    sides(1)=(side-1)*2+1
    sides(2)=sides(1)+1
    !post MPI receive
    do s =1,2
      direction=sides(s)
      mpi_recv_cnt(s) = 3*proc_num_atms! max estimate
      call mpi_irecv(mpi_recv_buf_dble(:,s), mpi_recv_cnt(s), mpi_double, proc_neighbor_rank(direction),&
                      3-s, pmemd_comm, recv_req(s), err_code_mpi)
    enddo
    do s =1,2
      direction=sides(s)
      ind=0
      do i = recv_buf_idx_range(1, direction), recv_buf_idx_range(2, direction)
        ind=ind+1
        mpi_send_buf_dble((ind-1)*3+1,s) = atm_frc(1,i)
        mpi_send_buf_dble((ind-1)*3+2,s) = atm_frc(2,i)
        mpi_send_buf_dble((ind-1)*3+3,s) = atm_frc(3,i)
      enddo
      !post MPI send
      mpi_send_cnt(s) =3*(recv_buf_idx_range(2, direction)-recv_buf_idx_range(1, direction)+1) !3*send_buf_cnt(direction)
      call mpi_isend(mpi_send_buf_dble(:,s), mpi_send_cnt(s), mpi_double, proc_neighbor_rank(direction),&
                     s, pmemd_comm, send_req(s), err_code_mpi)
    enddo
    !wait MPI receive and count
    do s =1,2
      call mpi_waitany(2, recv_req, recvd_idx, irecv_stat, err_code_mpi)
      call mpi_get_count(irecv_stat, mpi_double, mpi_recv_cnt(recvd_idx), err_code_mpi)
      if(3*proc_num_atms .lt. mpi_recv_cnt(recvd_idx)) then
           write(0,"('ERROR: Force recv count does not match the expected. Expected: bound ',I9,' Recvd: ', I9)") &
                      3*proc_num_atms,  mpi_recv_cnt(recvd_idx) 
           call mexit(6,1)
      end if
    end do

    !unpack
    do s =1,2
      direction=sides(s)
      do j= 1,mpi_recv_cnt(s)/3
        atm_idx=send_buf_indices(j,direction)
        atm_frc(1,atm_idx) = atm_frc(1,atm_idx) + mpi_recv_buf_dble((j-1)*3+1,s)
        atm_frc(2,atm_idx) = atm_frc(2,atm_idx) + mpi_recv_buf_dble((j-1)*3+2,s)
        atm_frc(3,atm_idx) = atm_frc(3,atm_idx) + mpi_recv_buf_dble((j-1)*3+3,s)
      enddo
    end do
    !wait MPI send
    call mpi_waitall(2, send_req, isend_stat, err_code_mpi)
#if 0
    do s =1,2
      direction=sides(s)
      write(0, "(I1,'-', I0.3,' With neighbor ',I1,' rank:',I3,' Sent ',I9,' Recvd ',I9)") &
                    side, mytaskid,direction, proc_neighbor_rank(direction), mpi_send_cnt(s),  mpi_recv_cnt(s)
    end do
#endif

  end do
end subroutine comm_frc


!*******************************************************************************
!
! Subroutine:  comm_3dimensions_3ints
!
! Description: Performs the communication in the three directions of the domain
! of 3
!              components double array wrap_array_param has to be set to true
!              in the case of crd data
!
!*******************************************************************************
subroutine comm_3dimensions_3ints(int3_array, wrap_array_param)
  use processor_mod
  ! Parameters
  integer :: int3_array(3, proc_atm_alloc_size)
  Logical, optional :: wrap_array_param

  integer ghost_cnt
  Logical          is_wrapped

  ! array wrapping switch setup
  is_wrapped=.false.
  if(present(wrap_array_param)) then
    is_wrapped = wrap_array_param
  end if

  call comm_1dimension_3ints(int3_array, ghost_cnt, 1, is_wrapped)
  call comm_1dimension_3ints(int3_array, ghost_cnt, 2, is_wrapped)
  call comm_1dimension_3ints(int3_array, ghost_cnt, 3, is_wrapped)
end subroutine comm_3dimensions_3ints

!*******************************************************************************
!
! Subroutine:  comm_1dimension_3ints
!
! Description: Performs the communication in one direction of the domain of 3
!              components double array
!
!*******************************************************************************
subroutine comm_1dimension_3ints(int3_array, ghost_cnt, side, is_wrapped)
  use processor_mod
  use pbc_mod
  use pmemd_lib_mod,     only : mexit
  ! Parameters
  Logical   is_wrapped
  integer ghost_cnt, side
  integer int3_array(3, proc_atm_alloc_size)

  ! local variables
  integer mpi_recv_cnt(2), mpi_send_cnt(2)! 2 elements for two directions
  integer        :: recv_req(2), send_req(2),recvd_idx
  integer       :: irecv_stat(mpi_status_size)
  integer       :: isend_stat(mpi_status_size, 2)
  integer i, j, k, s, atm_idx, direction, sides(2)

  if(side .eq. 1) ghost_cnt=0
  ! with values 1 to 6 in the order: x-, x+, y-, y+, z-, z+
  sides(1)=(side-1)*2+1
  sides(2)=sides(1)+1
  !post MPI receive
  do s =1,2
    direction=sides(s)
    mpi_recv_cnt(s) = 3*recv_buf_cnt(direction)
    call mpi_irecv(mpi_recv_buf_int(:,s), mpi_recv_cnt(s), mpi_integer, &
                   proc_neighbor_rank(direction), 3-s, pmemd_comm,&
                   recv_req(s), err_code_mpi)
  enddo
  do s =1,2
    direction=sides(s)
    !pack the send data
    do i =1, send_buf_cnt(direction)
      atm_idx=send_buf_indices(i,direction)
      mpi_send_buf_int((i-1)*3+1,s) = int3_array(1,atm_idx)
      mpi_send_buf_int((i-1)*3+2,s) = int3_array(2,atm_idx)
      mpi_send_buf_int((i-1)*3+3,s) = int3_array(3,atm_idx)
    enddo
    !post MPI send
    mpi_send_cnt(s) = 3*send_buf_cnt(direction)
    call mpi_isend(mpi_send_buf_int(:,s), mpi_send_cnt(s), mpi_integer,&
                   proc_neighbor_rank(direction), s, pmemd_comm,&
                   send_req(s), err_code_mpi)
  enddo
  !wait MPI receive and count
  do s =1,2
    call mpi_waitany(2, recv_req, recvd_idx, irecv_stat, err_code_mpi)
    call mpi_get_count(irecv_stat, mpi_integer, mpi_recv_cnt(recvd_idx), err_code_mpi)
    direction=sides(recvd_idx)
    if(3*recv_buf_cnt(direction) .ne. mpi_recv_cnt(recvd_idx)) then
         write(0,"('ERROR: 3int recv count does not match the expected. Expected: ',I9,' Recvd: ', I9)")&
                   recv_buf_cnt(direction),  mpi_recv_cnt(recvd_idx)
         call mexit(6,1)
    end if
  end do

  !unpack
  do s =1,2
    direction=sides(s)
    if(is_wrapped) then
      do j= 1,mpi_recv_cnt(s)/3
        ghost_cnt = ghost_cnt+1
        int3_array(1,proc_num_atms+ghost_cnt) = mpi_recv_buf_int((j-1)*3+1,s) + proc_side_wrap_offset(1,direction)
        int3_array(2,proc_num_atms+ghost_cnt) = mpi_recv_buf_int((j-1)*3+2,s) + proc_side_wrap_offset(2,direction)
        int3_array(3,proc_num_atms+ghost_cnt) = mpi_recv_buf_int((j-1)*3+3,s) + proc_side_wrap_offset(3,direction)
      enddo
    else
      do j= 1,mpi_recv_cnt(s)/3
        ghost_cnt = ghost_cnt+1
        int3_array(1,proc_num_atms+ghost_cnt) = mpi_recv_buf_int((j-1)*3+1,s)
        int3_array(2,proc_num_atms+ghost_cnt) = mpi_recv_buf_int((j-1)*3+2,s)
        int3_array(3,proc_num_atms+ghost_cnt) = mpi_recv_buf_int((j-1)*3+3,s)
      enddo
    endif
  end do
  !wait MPI send
  call mpi_waitall(2, send_req, isend_stat, err_code_mpi)
#if 0
  do s =1,2
    direction=sides(s)
    write(0, "(I1,'-', I0.3,' With neighbor ',I1,' rank:',I3,' Sent ',I9,' Recvd ',I9)") &
                  side, mytaskid,direction, proc_neighbor_rank(direction), mpi_send_cnt(s),  mpi_recv_cnt(s)
  end do
#endif
end subroutine comm_1dimension_3ints

!*******************************************************************************
!
! Subroutine:  comm_3dimensions_3dbls
!
! Description: Performs the communication in the three directions of the domain of 3
!              components double array wrap_array_param has to be set to true
!              in the case of crd data
!
!*******************************************************************************
subroutine comm_3dimensions_3dbls(dbl3_array, wrap_array_param)
  use processor_mod
  ! Parameters
  double precision dbl3_array(3,proc_atm_alloc_size)
  Logical, optional :: wrap_array_param

  integer ghost_cnt
  Logical          is_wrapped

  ! array wrapping switch setup
  is_wrapped=.false.
  if(present(wrap_array_param)) then
    is_wrapped = wrap_array_param
  end if

  call comm_1dimension_3dbls(dbl3_array, ghost_cnt, 1, is_wrapped)
  call comm_1dimension_3dbls(dbl3_array, ghost_cnt, 2, is_wrapped)
  call comm_1dimension_3dbls(dbl3_array, ghost_cnt, 3, is_wrapped)
end subroutine comm_3dimensions_3dbls

!*******************************************************************************
!
! Subroutine:  comm_1dimension_3dbls
!
! Description: Performs the communication in one direction of the domain of 3
!              components double array 
!
!*******************************************************************************
subroutine comm_1dimension_3dbls(dbl3_array, ghost_cnt, side, is_wrapped)
  use processor_mod
  use pmemd_lib_mod,     only : mexit

  ! Parameters
  Logical   is_wrapped
  integer ghost_cnt, side
  double precision dbl3_array(3,proc_atm_alloc_size)
  
  ! local variables
  integer mpi_recv_cnt(2), mpi_send_cnt(2)! 2 elements for two directions
  integer        :: recv_req(2), send_req(2),recvd_idx
  integer       :: irecv_stat(mpi_status_size)
  integer       :: isend_stat(mpi_status_size, 2)
  integer i, j, k, s, atm_idx, direction, sides(2)

  if(side .eq. 1) ghost_cnt=0
  ! with values 1 to 6 in the order: x-, x+, y-, y+, z-, z+
  sides(1)=(side-1)*2+1
  sides(2)=sides(1)+1
  !post MPI receive
  do s =1,2
    direction=sides(s)
    mpi_recv_cnt(s) = 3*recv_buf_cnt(direction)
    call mpi_irecv(mpi_recv_buf_dble(:,s), mpi_recv_cnt(s), mpi_double, proc_neighbor_rank(direction),&
                    3-s, pmemd_comm, recv_req(s), err_code_mpi)
  enddo
  do s =1,2
    direction=sides(s)
    !pack the send data
    do i =1, send_buf_cnt(direction)
      atm_idx=send_buf_indices(i,direction)
      mpi_send_buf_dble((i-1)*3+1,s) = dbl3_array(1,atm_idx)
      mpi_send_buf_dble((i-1)*3+2,s) = dbl3_array(2,atm_idx)
      mpi_send_buf_dble((i-1)*3+3,s) = dbl3_array(3,atm_idx)
    enddo
    !post MPI send
    mpi_send_cnt(s) = 3*send_buf_cnt(direction)
    call mpi_isend(mpi_send_buf_dble(:,s), mpi_send_cnt(s), mpi_double, proc_neighbor_rank(direction),&
                   s, pmemd_comm, send_req(s), err_code_mpi)
  enddo
  !wait MPI receive and count
  do s =1,2
    call mpi_waitany(2, recv_req, recvd_idx, irecv_stat, err_code_mpi)
    call mpi_get_count(irecv_stat, mpi_double, mpi_recv_cnt(recvd_idx), err_code_mpi)
    direction=sides(recvd_idx)
    if(3*recv_buf_cnt(direction) .ne. mpi_recv_cnt(recvd_idx)) then
         write(0,"('ERROR: 3dbl recv count does not match the expected. Expected: ',I9,' Recvd: ', I9)")&
                  recv_buf_cnt(direction),  mpi_recv_cnt(recvd_idx) 
         call mexit(6,1)
    end if
  end do

  !unpack
  do s =1,2
    direction=sides(s)
    if(is_wrapped) then ! Preserve the recived atoms count from each side for frc communication
      recv_buf_idx_range(1,direction) = proc_num_atms+ghost_cnt+1
      recv_buf_idx_range(2,direction) = proc_num_atms+ghost_cnt+mpi_recv_cnt(s)/3
    endif
    if(is_wrapped) then
      do j= 1,mpi_recv_cnt(s)/3
        ghost_cnt = ghost_cnt+1
        dbl3_array(1,proc_num_atms+ghost_cnt) = mpi_recv_buf_dble((j-1)*3+1,s) + proc_side_wrap_offset(1,direction)*box_len(1)
        dbl3_array(2,proc_num_atms+ghost_cnt) = mpi_recv_buf_dble((j-1)*3+2,s) + proc_side_wrap_offset(2,direction)*box_len(2)
        dbl3_array(3,proc_num_atms+ghost_cnt) = mpi_recv_buf_dble((j-1)*3+3,s) + proc_side_wrap_offset(3,direction)*box_len(3)
      enddo
    else
      do j= 1,mpi_recv_cnt(s)/3
        ghost_cnt = ghost_cnt+1
        dbl3_array(1,proc_num_atms+ghost_cnt) = mpi_recv_buf_dble((j-1)*3+1,s)
        dbl3_array(2,proc_num_atms+ghost_cnt) = mpi_recv_buf_dble((j-1)*3+2,s)
        dbl3_array(3,proc_num_atms+ghost_cnt) = mpi_recv_buf_dble((j-1)*3+3,s)
      enddo
    endif
  end do
  !wait MPI send
  call mpi_waitall(2, send_req, isend_stat, err_code_mpi)
#if 0
  do s =1,2
    direction=sides(s)
    write(0, "(I1,'-', I0.3,' With neighbor ',I1,' rank:',I3,' Sent ',I9,' Recvd ',I9)") &
                  side, mytaskid,direction, proc_neighbor_rank(direction), mpi_send_cnt(s),  mpi_recv_cnt(s)
  end do
#endif
end subroutine comm_1dimension_3dbls

!*******************************************************************************
!
! Subroutine:  comm_setup
!
! Description: Setup the communicates information of the neighbor MPI tasks (proc_neighbor_rank)
!                and the communication region (proc_send_bounds) at first call
!              Also, initialize the list of atoms contributing in the communication at every call
!                send_buf_indices contains the indices to be communicated to
!                each of the 6 sides and the size is stored in send_buf_cnt
!              The ghost coordinates are exchanged here as part of the initialization
!
!*******************************************************************************
subroutine comm_setup()
! TERMENOLOGY: The arrays related to the 6 sides of the box are
!   stored in the order: x-, x+, y-, y+, z-, z+
  use processor_mod
  use ensure_alloc_mod

  ! Local variables
  integer i, j, k, direction
  double precision ghost_len(3)
  integer local_send_buf_cnt(6) ! for debugging
  integer old_alloc_size

  ! Get communication coordinates bounds
  ghost_len(1) = new_bkt_size(1)
  ghost_len(2) = new_bkt_size(2)
  ghost_len(3) = new_bkt_size(3)

  proc_send_bounds(1) = proc_min_x_crd + ghost_len(1)
  proc_send_bounds(2) = proc_max_x_crd - ghost_len(1)
  proc_send_bounds(3) = proc_min_y_crd + ghost_len(2)
  proc_send_bounds(4) = proc_max_y_crd - ghost_len(2)
  proc_send_bounds(5) = proc_min_z_crd + ghost_len(3)
  proc_send_bounds(6) = proc_max_z_crd - ghost_len(3)


  ! find the atoms index for sending data
  send_buf_cnt = 0
  do i = 1,proc_num_atms
    k=1
    do j=1,5,2 ! lower index bounds
      if(proc_atm_crd(k,i) .le. proc_send_bounds(j)) then
        send_buf_cnt(j) = send_buf_cnt(j) + 1
        send_buf_indices(send_buf_cnt(j), j) = i
      end if
      k=k+1
    end do
    k=1
    do j=2,6,2 ! higher index bounds
      if(proc_atm_crd(k,i) .ge. proc_send_bounds(j)) then
        send_buf_cnt(j) = send_buf_cnt(j) + 1
        send_buf_indices(send_buf_cnt(j), j) = i
      end if
      k=k+1
    end do
  end do
  local_send_buf_cnt = send_buf_cnt 

! Setup the ghost indices by communicating the ghost coordinates

! exchange the local crd in x-/x+ sides
  call comm_ensure_ghost_alloc(1) 
  call comm_1dimension_3dbls(proc_atm_crd, proc_ghost_num_atms, 1, .true.)

! updated the y-/y+ buffers with the new x-/x+ ghost indices
  do i = proc_num_atms+1,proc_num_atms+proc_ghost_num_atms
    if(proc_atm_crd(2,i) .le. proc_send_bounds(3)) then
      send_buf_cnt(3) = send_buf_cnt(3) + 1
      send_buf_indices(send_buf_cnt(3), 3) = i
    end if
    if(proc_atm_crd(2,i) .ge. proc_send_bounds(4)) then
      send_buf_cnt(4) = send_buf_cnt(4) + 1
      send_buf_indices(send_buf_cnt(4), 4) = i
    end if
  end do
! exchage the local/ghost crd in y-/y+ sides
  call comm_ensure_ghost_alloc(2)
  call comm_1dimension_3dbls(proc_atm_crd, proc_ghost_num_atms, 2, .true.)

! updated the z-/z+ buffers with the new x-/x+ and y-/y+ ghost indices
  do i = proc_num_atms+1,proc_num_atms+proc_ghost_num_atms
    if(proc_atm_crd(3,i) .le. proc_send_bounds(5)) then
      send_buf_cnt(5) = send_buf_cnt(5) + 1
      send_buf_indices(send_buf_cnt(5), 5) = i
    end if
    if(proc_atm_crd(3,i) .ge. proc_send_bounds(6)) then
      send_buf_cnt(6) = send_buf_cnt(6) + 1
      send_buf_indices(send_buf_cnt(6), 6) = i
    end if
  end do

! exchage the local/ghost crd in z-/z+ sides (this is to complete the data)
  call comm_ensure_ghost_alloc(3)
  call comm_1dimension_3dbls(proc_atm_crd, proc_ghost_num_atms, 3, .true.)

  ! Update the arrays size if the new local+ghost atoms > current allocation
  old_alloc_size = size(proc_atm_qterm)
  if(old_alloc_size .lt. proc_num_atms_min_bound+proc_ghost_num_atms) &
    call ensure_arrays_with_ghost(proc_num_atms_min_bound+proc_ghost_num_atms)

#if 0 /* Debugging send buffers data*/
! Debug the index ranges of the ghost data from the 6 sides of the box
#if 0 
    write(0, "(I0.3,'-0',' Ghost total: ',I,' Sides ranges: ',6(2(I8,','),'|'))") mytaskid, proc_ghost_num_atms, recv_buf_idx_range
    write(0, "(I0.3,'-1',' Sides size send: ',6(I8,' |'), ' recv: ',6(I8,' |'))") mytaskid, send_buf_cnt, recv_buf_idx_range(2,:)-recv_buf_idx_range(1,:)+1
!    call MPI_Barrier(pmemd_comm)
!    stop
#endif

! Debug current and neighbots corrdinates and their MPI ranks
#if 0 
    write(0, "(I0.3,' pcrd:(',3I3,')',' nb:(',6(3I3,','),')', 6I3)") mytaskid,p(:), np(:,:), proc_neighbor_rank
    call MPI_Barrier(pmemd_comm)
    stop
#endif
! Debug boundary processors for crd offset
#if 0
    write(0, "('000-',I0.3,' pcrd:(',3I3,')',' wrapping:',6(3F11.5,','))") mytaskid,p(:), proc_side_wrap_offset
!    call MPI_Barrier(pmemd_comm)
!    stop
#endif
#if 0
  write(0, "('1-', I0.3,' Ghost size: ',3F7.3,' Proc bounds: ',3(2F11.3, ','))") mytaskid, ghost_len, proc_min_x_crd, proc_max_x_crd, proc_min_y_crd, proc_max_y_crd, proc_min_z_crd, proc_max_z_crd
  write(0, "('2-', I0.3, ' Local region comm bounds: ', 6(F8.3, ','))") mytaskid, proc_send_bounds
!  call MPI_Barrier(pmemd_comm)
!  stop
#endif
! Check if the send buffer sizes
#if 0
  write(0,"('0-',I0.3,' Local boundary buf sizes:', 6I)") mytaskid, send_buf_cnt
#endif

#endif

end subroutine comm_setup

!*******************************************************************************
!
! Subroutine:  proc_pack_output
!
! Description: Uses proc_atm_cnt to build output array
!
!*******************************************************************************

subroutine proc_pack_output
  use processor_mod
  use parallel_dat_mod, only : master, master_rank  
  !integer            ::  mpi_atm_cnt(:)
  integer j
  integer,save :: numIntEle=1  ! global id, 333 proc wrap counter. Change these as needed for allocation
  integer,save :: numDouEle=9  ! 3 coord, 3 vel. Change these as needed for allocation
  integer,save :: numMPI=26    ! max number of MPI to send
  !the col and row are switched, since in MPI fortran sends and recvs in col
  !Major, so fast changing part is row. Now, row = atom lists, col= MPI ranks
  !send_atm_dbl = previous pack_atm_array_double
  !send_atm_int = previous pack_atm_array_int

  send_output_dbl(1,1) = dble(master_rank)  
  send_output_int(1,1) = master_rank 
  ! Second element is size of the mpi array 
  send_output_dbl(2,1) = dble(proc_num_atms*numDouEle)
  send_output_int(2,1) = proc_num_atms*numIntEle

! write(0,*)"Atms:",proc_num_atms

  do j=1, proc_num_atms
    ! Double arrays
    send_output_dbl(3*(j-1)+3,1)=proc_atm_crd(1,j)
    send_output_dbl(3*(j-1)+4,1)=proc_atm_crd(2,j)
    send_output_dbl(3*(j-1)+5,1)=proc_atm_crd(3,j)
    send_output_dbl(3*(j-1)+3*proc_num_atms+3,1)=proc_atm_vel(1,j)
    send_output_dbl(3*(j-1)+3*proc_num_atms+4,1)=proc_atm_vel(2,j)
    send_output_dbl(3*(j-1)+3*proc_num_atms+5,1)=proc_atm_vel(3,j)
    send_output_dbl(3*(j-1)+3*proc_num_atms+3,1)=proc_atm_wrap(1,j)
    send_output_dbl(3*(j-1)+3*proc_num_atms+4,1)=proc_atm_wrap(2,j)
    send_output_dbl(3*(j-1)+3*proc_num_atms+5,1)=proc_atm_wrap(3,j)
        
    ! Integer arrays
    send_output_int(2+j,1)=proc_atm_to_full_list(j)

  end do 

!do i=1,size(recv_task_map)
!      if(recv_task_map(i) .ne. -1) then
!        write(15+mytaskid,*)"recv_task_map(",i,"):", recv_task_map(i)
!      end if
!      if(send_task_map(i) .ne. -1) then
!        write(15+mytaskid,*)"send_task_map(",i,"):", send_task_map(i)
!      end if
! end do
end subroutine proc_pack_output

!*******************************************************************************
!
! Subroutine:  proc_unpack_output
!
! Description: Uses proc_atm_cnt
!
!*******************************************************************************

subroutine proc_unpack_output(crd, vel)
  use processor_mod

  integer i,j,k
  integer numAtms
  integer,save :: numIntEle=1  ! global id, wrap counter, Change these as needed for allocation
  integer,save :: numDouEle=9  ! 3 coord, 3 velocity. Change these as needed for allocation
  integer,save :: numMPI=26    ! max number of MPI to receive, Change to 13 later
  double precision :: crd(3,*)
  double precision :: vel(3,*)

  do i=1,numtasks-1
!    if(i-1 .eq. master_rank) cycle
    numAtms = int(recv_output_dbl(2,i)/numDouEle) !# of atoms received
    do j=1,numAtms !2nd indice is how many atoms
      k = recv_output_int(2+j,i)
!     write(0,*)k,numAtms,j
      crd(1,k) =recv_output_dbl(3*(j-1)+3,i)
      crd(2,k) =recv_output_dbl(3*(j-1)+4,i)
      crd(3,k) =recv_output_dbl(3*(j-1)+5,i)
      vel(1,k) =recv_output_dbl(3*(j-1)+3*numAtms+3,i)
      vel(2,k) =recv_output_dbl(3*(j-1)+3*numAtms+4,i)
      vel(3,k) =recv_output_dbl(3*(j-1)+3*numAtms+5,i)
    end do
  end do
    
  ! Unpack self
  do i=1, proc_num_atms
      crd(1,proc_atm_to_full_list(i)) = proc_atm_crd(1,i)
      crd(2,proc_atm_to_full_list(i)) = proc_atm_crd(2,i)
      crd(3,proc_atm_to_full_list(i)) = proc_atm_crd(3,i)
      vel(1,proc_atm_to_full_list(i)) = proc_atm_vel(1,i)
      vel(2,proc_atm_to_full_list(i)) = proc_atm_vel(2,i)
      vel(3,proc_atm_to_full_list(i)) = proc_atm_vel(3,i)
  end do

! write(15+mytaskid,*)"ghost atm cnt ",proc_ghost_num_atms
!  do i=1,recv_atm_dbl(2,1)+2
!    write(15+mytaskid,*)"recv: ", i, &
!       recv_atm_dbl(i,1)
!  end do
! do i=proc_num_atms+1,proc_num_atms+proc_ghost_num_atms
!   write(15+mytaskid,*)"iac: ",proc_iac(i), "mass: ", &
!        proc_atm_mass(i),proc_old_mpi(i)
! end do
end subroutine proc_unpack_output

!*******************************************************************************
!
! Subroutine:  output_distribution
!
! Description: Send and Recieve atm crd and relevant parameters called during
! the output process.
!
!*******************************************************************************

subroutine output_distribution()

  use processor_mod, only: send_output_dbl, send_output_int, &
                           recv_output_dbl, recv_output_int, &
                           proc_num_atms, &
                           proc_atm_alloc_size, proc_num_atms_min_bound
  use parallel_dat_mod, only: mytaskid, pmemd_comm,err_code_mpi,&
                              numtasks,master,master_rank
  implicit none

! Local common variables:
  integer, parameter    :: gifd_tag_int = 12 
  integer, parameter    :: gifd_tag_dbl = 13 
                                                        !recv_tasks= max tasks recv from
  integer       :: i
  integer        :: wait_call
  integer        :: taskmap_idx, taskmap_idx_map
  integer        :: node
!Recv variables
  integer       :: recv_atm_cnt
  integer       :: recv_atm_cnts_sum
  integer       :: receivers_cnt
  integer       :: recv_task_dbl, recv_task_int
  !currently 26 is the max MPI ranks that can recv or send to one MPI
  !thus hard coded, but we can generalize later
  integer       :: recv_atm_cnts_int(26) ! may not required since
                                                 ! send_atm_int and
                                                 ! _double contains the count
  integer        :: recv_req_int(numtasks-1)
  integer        :: recv_req_dbl(numtasks-1)
  logical, save  :: recv_atm_cnts_zeroed = .false.
  integer       :: irecv_stat(mpi_status_size)
!Send variables 
  integer       :: send_atm_cnt
  integer       :: send_task_dbl, send_task_int
  integer       :: send_req_dbl
  integer       :: send_req_int
  integer       :: isend_stat(mpi_status_size, 26)
  integer       ::  recv_taskmap(numtasks-1)

  !integer       :: problem ; it will catch if recv_buf size is smaller than
  !recv bytes in the MPI_Irecv call. Actual recv bytes is much smaller though. 
  !probably need to address this, but currently working without fixing this.
  
  ! Post the asynchronous receives first:
  if(master) then
    do taskmap_idx = 1, numtasks
      if(taskmap_idx-1 .eq. master_rank) cycle
      recv_task_dbl = taskmap_idx - 1
      recv_task_int = taskmap_idx - 1
      recv_output_int(1,taskmap_idx) = recv_task_int
      recv_output_dbl(1,taskmap_idx) = recv_task_dbl
      taskmap_idx_map = taskmap_idx
      if(taskmap_idx_map .gt. master_rank) taskmap_idx_map = taskmap_idx_map-1
        call mpi_irecv(recv_output_int(3, taskmap_idx_map), proc_atm_alloc_size, mpi_integer, &
                     recv_task_int, gifd_tag_int, pmemd_comm, &
                     recv_req_int(taskmap_idx_map), err_code_mpi); MPI_ERR_CHK("in output_distribution irecv int")
        call mpi_irecv(recv_output_dbl(3,taskmap_idx_map), 9*proc_atm_alloc_size, mpi_double, &
                     recv_task_dbl, gifd_tag_dbl, pmemd_comm, &
                     recv_req_dbl(taskmap_idx_map), err_code_mpi); MPI_ERR_CHK("in output_distribution irecv dbl")
    end do
  end if

  ! Now set up and post the asynchronous sends:
  if(.not. master) then
      send_task_dbl = int(send_output_dbl(1, 1)) ! the mpi rank to send
      send_task_int = send_output_int(1,1) ! the mpi rank to send

      !sending the integers
      send_atm_cnt = send_output_int(2,1) ! # of integers to send 

      call mpi_isend(send_output_int(3,1), send_atm_cnt, &
                     mpi_integer, send_task_int, gifd_tag_int, pmemd_comm, &
                     send_req_int, err_code_mpi); MPI_ERR_CHK("in output_distribution isend int")

      !sending the doubles
      send_atm_cnt = int(send_output_dbl(2,1)) ! # of doubles to send

      call mpi_isend(send_output_dbl(3,1), send_atm_cnt, &
                     mpi_double, send_task_dbl, gifd_tag_dbl, pmemd_comm, &
                     send_req_dbl, err_code_mpi); MPI_ERR_CHK("in output_distribution isend dbl")
  end if 

  ! Wait on and process the pending receive requests:
  if(master) then
    do wait_call = 1, numtasks -1
      !Now for integer
      call mpi_waitany(numtasks-1, recv_req_int, taskmap_idx, irecv_stat, &
                       err_code_mpi); MPI_ERR_CHK("in output_distribution wait recv int")
      !taskmap idx get the value index from which recv_task_dbl derived
      !recv_task_int = recv_task_map(taskmap_idx)
      call mpi_get_count(irecv_stat, mpi_integer, recv_atm_cnt, err_code_mpi); MPI_ERR_CHK("in output_distribution get recv count int")

      !recv_atm_cnts_int(recv_task_int) = recv_atm_cnt ! may not required
      recv_output_int(2,taskmap_idx) = recv_atm_cnt !storing the count
     ! print *, "from recv waitcall int ", recv_task_int, taskmap_idx, mytaskid, recv_atm_cnt
      !if (problem .eq. 1) then
  
      !end if

      !Now for double
      call mpi_waitany(numtasks-1, recv_req_dbl, taskmap_idx, irecv_stat, &
                       err_code_mpi); MPI_ERR_CHK("in output_distribution wait recv dbl")
      !taskmap idx find the index from which recv_task_dbl derived
      !recv_task_dbl = recv_task_map(taskmap_idx)

      call mpi_get_count(irecv_stat, mpi_double, recv_atm_cnt, err_code_mpi); MPI_ERR_CHK("in output_distribution get recv count dbl")

      recv_output_dbl(2,taskmap_idx) = dble(recv_atm_cnt) !storing the count
!   write(15+mytaskid,*) "received atoms", recv_atm_cnt,"recv_task_dbl ", recv_task_dbl,"mytaskid", mytaskid
    end do
  end if

  if(.not. master) then
    ! Wait for all sends to complete:
    call mpi_waitall(1, send_req_int, isend_stat, err_code_mpi); MPI_ERR_CHK("in output_distribution wait send int")
    call mpi_waitall(1, send_req_dbl, isend_stat, err_code_mpi); MPI_ERR_CHK("in output_distribution wait send dbl")
  end if
 return

end subroutine output_distribution

subroutine mpi_output_transfer(crd, vel)
 
  use parallel_dat_mod, only : master
 
  implicit none

  double precision :: crd(3,*)
  double precision :: vel(3,*)

  if(.not. master) call proc_pack_output
  call output_distribution
  if(master)  call proc_unpack_output(crd,vel)

end subroutine mpi_output_transfer

!
! This subroutine moves atoms and calls 3 routines
!
!
subroutine proc_reorient_volume(counts, nstep)

  use processor_mod
  use file_io_dat_mod
  implicit none

  integer :: counts
  integer :: nstep
  integer delete_cnt
  integer :: deleted_array(proc_num_atms_min_bound)
  integer :: dest_mpi(proc_num_atms_min_bound)
  integer :: pack_array_int(proc_num_atms_min_bound*6,neighbor_mpi_cnt)
  double precision :: pack_array_dbl(proc_num_atms_min_bound*14, neighbor_mpi_cnt)
  integer :: unpack_array_int(proc_num_atms_min_bound*6, neighbor_mpi_cnt)
  double precision :: unpack_array_dbl(proc_num_atms_min_bound*14, neighbor_mpi_cnt)
  integer i
  double precision :: crd(3)
  reorient_flag = .true.
  
!  write(15+mytaskid,*)"proc_num_pre_orient",proc_num_atms
!  write(0,*)"proc_num_pre_orient",proc_num_atms
  delete_cnt = 0
  call proc_pack_reorient(deleted_array, delete_cnt, pack_array_int, pack_array_dbl, dest_mpi)
  del_cnt = delete_cnt
  atms_to_del(1:del_cnt) = deleted_array(1:del_cnt)
  dest_sdmns(1:del_cnt)=dest_mpi(1:del_cnt)
  call bonded_packing(proc_atm_crd, neighbor_mpi_cnt,counts,nstep)! Calling this here is needed
  call proc_save_bonded(counts) ! saving the bonded index's global id, to use later
  call reorient_distribution(pack_array_int, pack_array_dbl, unpack_array_int, unpack_array_dbl)
  call proc_unpack_reorient(deleted_array, delete_cnt, unpack_array_int, unpack_array_dbl)
!  write(15+mytaskid,*)"proc_num_post_orient",proc_num_atms
!  write(0,*)"proc_num_post_orient",proc_num_atms

end subroutine proc_reorient_volume

subroutine proc_pack_reorient(deleted_array, delete_cnt, pack_array_int, pack_array_dbl,dest_mpi)

  use pbc_mod
  use processor_mod
  implicit none

  integer :: deleted_array(proc_num_atms_min_bound)
  integer :: dest_mpi(proc_num_atms_min_bound)
  integer :: delete_cnt
  integer :: pack_array_int(proc_num_atms_min_bound*6,neighbor_mpi_cnt)
  double precision :: pack_array_dbl(proc_num_atms_min_bound*14,neighbor_mpi_cnt)
  double precision   :: crd(3)
  double precision   :: recipbox(3)
  integer i, j, bkt_x, bkt_y, bkt_z, n, blockid
  integer              :: numIntEle = 6 
  integer              :: numDouEle = 14

  delete_cnt=0
  recipbox(:) = 1.d0 / dble(pbc_box(:))
  pack_array_dbl=0.0d0
  pack_array_int=0

  do j=1,neighbor_mpi_cnt
    ! First element is the mpi sending to
    pack_array_dbl(1,j) = dble(neighbor_mpi_map(j)) !i - 1 !since actual range is 0 to numtasks -1
    pack_array_int(1,j) = neighbor_mpi_map(j) !i - 1 !since actual range is 0 to numtasks-1
  end do
  do i=1,proc_num_atms
    crd(:)=mod(proc_atm_crd(:,i)*recipbox(:),1.0d0)*pbc_box(:)
    if(crd(1) .lt. 0.0d0) crd(1)=crd(1)+pbc_box(1)
    if(crd(2) .lt. 0.0d0) crd(2)=crd(2)+pbc_box(2)
    if(crd(3) .lt. 0.0d0) crd(3)=crd(3)+pbc_box(3)

    call destinationmpi(blockid, proc_atm_crd(:,i))
    if(blockid .ne. mytaskid) then
      delete_cnt=delete_cnt+1
      deleted_array(delete_cnt) = i
      dest_mpi(delete_cnt)=blockid
      ! Packing atom
      do j=1,neighbor_mpi_cnt
        if(neighbor_mpi_map(j) .eq. blockid) then
          ! Double arrays
          pack_array_dbl(int(pack_array_dbl(2,j))+3,j)= crd(1)
          pack_array_dbl(int(pack_array_dbl(2,j))+4,j)= crd(2)
          pack_array_dbl(int(pack_array_dbl(2,j))+5,j)= crd(3)
          pack_array_dbl(int(pack_array_dbl(2,j))+6,j)= proc_atm_vel(1,i)
          pack_array_dbl(int(pack_array_dbl(2,j))+7,j)= proc_atm_vel(2,i)
          pack_array_dbl(int(pack_array_dbl(2,j))+8,j)= proc_atm_vel(3,i)
          pack_array_dbl(int(pack_array_dbl(2,j))+9,j)=proc_atm_mass(i)
          pack_array_dbl(int(pack_array_dbl(2,j))+10,j)=proc_atm_qterm(i)
          pack_array_dbl(int(pack_array_dbl(2,j))+11,j)=proc_atm_last_vel(1,i)
          pack_array_dbl(int(pack_array_dbl(2,j))+12,j)=proc_atm_last_vel(2,i)
          pack_array_dbl(int(pack_array_dbl(2,j))+13,j)=proc_atm_last_vel(3,i)
          pack_array_dbl(int(pack_array_dbl(2,j))+14,j)=proc_atm_frc(1,i)
          pack_array_dbl(int(pack_array_dbl(2,j))+15,j)=proc_atm_frc(2,i)
          pack_array_dbl(int(pack_array_dbl(2,j))+16,j)=proc_atm_frc(3,i)
        
          ! Integer arrays
          pack_array_int(pack_array_int(2,j)+3,j)=proc_atm_to_full_list(i)
          pack_array_int(pack_array_int(2,j)+4,j)=i  ! Its good to send the local index of the real atom too
          pack_array_int(pack_array_int(2,j)+5,j)=proc_iac(i)
          pack_array_int(pack_array_int(2,j)+6,j)=proc_atm_wrap(1,i)
          pack_array_int(pack_array_int(2,j)+7,j)=proc_atm_wrap(2,i)
          pack_array_int(pack_array_int(2,j)+8,j)=proc_atm_wrap(3,i)

          if(proc_atm_crd(1,i) .ge. pbc_box(1)) then
            pack_array_int(pack_array_int(2,j)+6,j)=proc_atm_wrap(1,i)+1
          end if
          if(proc_atm_crd(2,i) .ge. pbc_box(2)) then
            pack_array_int(pack_array_int(2,j)+7,j)=proc_atm_wrap(2,i)+1
          end if
          if(proc_atm_crd(3,i) .ge. pbc_box(3)) then
            pack_array_int(pack_array_int(2,j)+8,j)=proc_atm_wrap(3,i)+1
          end if
          if(proc_atm_crd(1,i) .lt. 0.d0) then
            pack_array_int(pack_array_int(2,j)+6,j)=proc_atm_wrap(1,i)-1
          end if
          if(proc_atm_crd(2,i) .lt. 0.d0) then
            pack_array_int(pack_array_int(2,j)+7,j)=proc_atm_wrap(2,i)-1
          end if
          if(proc_atm_crd(3,i) .lt. 0.d0) then
            pack_array_int(pack_array_int(2,j)+8,j)=proc_atm_wrap(3,i)-1
          end if

          ! Second element is size of the mpi array
          pack_array_dbl(2,j) = pack_array_dbl(2,j) + dble(numDouEle)
          pack_array_int(2,j) = pack_array_int(2,j) + numIntEle
        end if
      end do
    else
      if(proc_atm_crd(1,i) .ge. pbc_box(1)) then
        proc_atm_wrap(1,i)=proc_atm_wrap(1,i)+1
      end if
      if(proc_atm_crd(2,i) .ge. pbc_box(2)) then
        proc_atm_wrap(2,i)=proc_atm_wrap(2,i)+1
      end if
      if(proc_atm_crd(3,i) .ge. pbc_box(3)) then
        proc_atm_wrap(3,i)=proc_atm_wrap(3,i)+1
      end if
      if(proc_atm_crd(1,i) .lt. 0.d0) then
        proc_atm_wrap(1,i)=proc_atm_wrap(1,i)-1
      end if
      if(proc_atm_crd(2,i) .lt. 0.d0) then
        proc_atm_wrap(2,i)=proc_atm_wrap(2,i)-1
      end if
      if(proc_atm_crd(3,i) .lt. 0.d0) then
        proc_atm_wrap(3,i)=proc_atm_wrap(3,i)-1
      end if
      proc_atm_crd(1,i)=crd(1)
      proc_atm_crd(2,i)=crd(2)
      proc_atm_crd(3,i)=crd(3)
    end if
  end do

end subroutine proc_pack_reorient

subroutine proc_unpack_reorient(deleted_array, delete_cnt, unpack_array_int, unpack_array_dbl)

  use processor_mod, only : proc_atm_crd, proc_atm_vel, proc_atm_mass, &
                            proc_atm_qterm, proc_atm_to_full_list, &
                            proc_old_local_id, proc_iac, proc_old_mpi, &
                            neighbor_mpi_cnt, proc_num_atms, &
                            proc_num_atms_min_bound, proc_num_atms_min_size, &
                            proc_atm_wrap, proc_atm_last_vel, proc_atm_frc

  use pbc_mod
  implicit none

  integer :: deleted_array(proc_num_atms_min_bound)
  integer :: delete_cnt
  integer :: total_added, added_cnt
  integer :: i,j,n
  integer :: atm_idx, proc_num_atms_idx
  integer, intent(in) :: unpack_array_int(proc_num_atms_min_bound*6,neighbor_mpi_cnt)
  double precision, intent(in) :: unpack_array_dbl(proc_num_atms_min_bound*14,neighbor_mpi_cnt)
  integer, save :: numIntEle = 6
  integer, save :: numDouEle = 14
  added_cnt = 0
  total_added = 0

!  write(25+mytaskid,*)"Pre unpack", proc_num_atms

  ! First move end of list to beginning of list to prevent overwriting important
  ! data.
!   write(25+mytaskid,*)"total added delete cnt",total_added,delete_cnt
  proc_num_atms_idx=proc_num_atms
  do i=1,delete_cnt
    atm_idx = deleted_array(i)

    do while(any(proc_num_atms_idx .eq. deleted_array(1:delete_cnt)))
      proc_num_atms_idx=proc_num_atms_idx-1
    end do
    proc_atm_crd(1,atm_idx) = proc_atm_crd(1,proc_num_atms_idx)
    proc_atm_crd(2,atm_idx) = proc_atm_crd(2,proc_num_atms_idx)
    proc_atm_crd(3,atm_idx) = proc_atm_crd(3,proc_num_atms_idx)
    proc_atm_vel(1,atm_idx) = proc_atm_vel(1,proc_num_atms_idx)
    proc_atm_vel(2,atm_idx) = proc_atm_vel(2,proc_num_atms_idx)
    proc_atm_vel(3,atm_idx) = proc_atm_vel(3,proc_num_atms_idx)
    proc_atm_last_vel(1,atm_idx) = proc_atm_last_vel(1,proc_num_atms_idx)
    proc_atm_last_vel(2,atm_idx) = proc_atm_last_vel(2,proc_num_atms_idx)
    proc_atm_last_vel(3,atm_idx) = proc_atm_last_vel(3,proc_num_atms_idx)
    proc_atm_frc(1,atm_idx) = proc_atm_frc(1,proc_num_atms_idx)
    proc_atm_frc(2,atm_idx) = proc_atm_frc(2,proc_num_atms_idx)
    proc_atm_frc(3,atm_idx) = proc_atm_frc(3,proc_num_atms_idx)
    proc_atm_mass(atm_idx) = proc_atm_mass(proc_num_atms_idx)
    proc_atm_qterm(atm_idx) = proc_atm_qterm(proc_num_atms_idx)
    proc_atm_to_full_list(atm_idx) = proc_atm_to_full_list(proc_num_atms_idx)
    proc_old_local_id(atm_idx) = proc_old_local_id(proc_num_atms_idx)
    proc_iac(atm_idx) = proc_iac(proc_num_atms_idx)
    proc_atm_wrap(1,atm_idx) = proc_atm_wrap(1,proc_num_atms_idx)
    proc_atm_wrap(2,atm_idx) = proc_atm_wrap(2,proc_num_atms_idx)
    proc_atm_wrap(3,atm_idx) = proc_atm_wrap(3,proc_num_atms_idx)
    proc_old_mpi(atm_idx) = proc_old_mpi(proc_num_atms_idx)
    proc_num_atms_idx=proc_num_atms_idx-1
    proc_num_atms = proc_num_atms - 1
  end do

  ! Now fill the new atoms to the end of the list

  do i=1,neighbor_mpi_cnt
    do j=1,unpack_array_int(2,i)/numIntEle !2nd indice is how many atoms

      ! If we've filled in all empty array indicies make proc_num_atms bigger
      atm_idx = proc_num_atms + 1
      proc_num_atms = proc_num_atms + 1

      added_cnt = added_cnt + 1
      proc_atm_crd(1,atm_idx) = unpack_array_dbl((j-1)*numDouEle+3,i)
      proc_atm_crd(2,atm_idx) = unpack_array_dbl((j-1)*numDouEle+4,i)
      proc_atm_crd(3,atm_idx) = unpack_array_dbl((j-1)*numDouEle+5,i)
      proc_atm_vel(1,atm_idx) = unpack_array_dbl((j-1)*numDouEle+6,i)
      proc_atm_vel(2,atm_idx) = unpack_array_dbl((j-1)*numDouEle+7,i)
      proc_atm_vel(3,atm_idx) = unpack_array_dbl((j-1)*numDouEle+8,i)
      proc_atm_mass(atm_idx) = unpack_array_dbl((j-1)*numDouEle+9,i)
      proc_atm_qterm(atm_idx) = unpack_array_dbl((j-1)*numDouEle+10,i)
      proc_atm_last_vel(1,atm_idx) = unpack_array_dbl((j-1)*numDouEle+11,i)
      proc_atm_last_vel(2,atm_idx) = unpack_array_dbl((j-1)*numDouEle+12,i)
      proc_atm_last_vel(3,atm_idx) = unpack_array_dbl((j-1)*numDouEle+13,i)
      proc_atm_frc(1,atm_idx) = unpack_array_dbl((j-1)*numDouEle+14,i)
      proc_atm_frc(2,atm_idx) = unpack_array_dbl((j-1)*numDouEle+15,i)
      proc_atm_frc(3,atm_idx) = unpack_array_dbl((j-1)*numDouEle+16,i)
      proc_atm_to_full_list(atm_idx) = unpack_array_int((j-1)*numIntEle+3,i)
      proc_old_local_id(atm_idx)= unpack_array_int((j-1)*numIntEle+4,i)
      proc_iac(atm_idx) = unpack_array_int((j-1)*numIntEle+5,i)
      proc_old_mpi(atm_idx) = unpack_array_int(1,i) 
      proc_atm_wrap(1,atm_idx) = unpack_array_int((j-1)*numIntEle+6,i)
      proc_atm_wrap(2,atm_idx) = unpack_array_int((j-1)*numIntEle+7,i)
      proc_atm_wrap(3,atm_idx) = unpack_array_int((j-1)*numIntEle+8,i)
!      write(15+mytaskid,*)"unpacked",proc_atm_to_full_list(atm_idx)
    end do
  end do

!  write(25+mytaskid,*)"total added",total_added

  if(proc_num_atms .lt. proc_num_atms_min_size) then  ! temporary solution for empty sub domain
    proc_num_atms_min_bound = proc_num_atms_min_size
  else
    proc_num_atms_min_bound = proc_num_atms
  endif

end subroutine proc_unpack_reorient

subroutine reorient_distribution(pack_array_int, pack_array_dbl, unpack_array_int, unpack_array_dbl)

  use processor_mod, only: proc_num_atms, neighbor_mpi_map, &
                           neighbor_mpi_cnt, proc_num_atms_min_bound
  use parallel_dat_mod, only: mytaskid, pmemd_comm,err_code_mpi,numtasks

  implicit none

  integer :: pack_array_int(proc_num_atms_min_bound*6,neighbor_mpi_cnt)
  double precision :: pack_array_dbl(proc_num_atms_min_bound*14,neighbor_mpi_cnt)
  integer :: unpack_array_int(proc_num_atms_min_bound*6, neighbor_mpi_cnt)
  double precision :: unpack_array_dbl(proc_num_atms_min_bound*14, neighbor_mpi_cnt)

! Local common variables:
  integer, parameter    :: gifd_tag_int = 12 
  integer, parameter    :: gifd_tag_dbl = 13 
                                                        !neighbor_mpi_cnt= max tasks recv from
  integer       :: i
  integer        :: wait_call
  integer        :: taskmap_idx
  integer        :: node
!Recv variables
  integer       :: recv_atm_cnt
  integer       :: recv_atm_cnts_sum
  integer       :: receivers_cnt
  integer       :: recv_task_dbl, recv_task_int
  !currently 26 is the max MPI ranks that can recv or send to one MPI
  !thus hard coded, but we can generalize later
  integer       :: recv_atm_cnts_int(26) ! may not required since
                                                 ! send_atm_int and
                                                 ! _double contains the count
  integer       :: recv_atm_cnts_dbl(26)!may not be required 
  integer        :: recv_req_int(26)
  integer        :: recv_req_dbl(26)
  logical, save  :: recv_atm_cnts_zeroed = .false.
  integer       :: irecv_stat(mpi_status_size)
!Send variables 
  integer       :: send_atm_cnt
  integer       :: send_task_dbl, send_task_int
  integer       :: send_req_dbl(26)
  integer       :: send_req_int(26)
  integer       :: isend_stat(mpi_status_size, 26)
  !integer      :: owned_atm_cnts(0:tasks)
  !integer      :: off_tbl(0:tasks)
  integer       ::  recv_taskmap(numtasks -1)
  !integer       :: problem ; it will catch if recv_buf size is smaller than
  !recv bytes in the MPI_Irecv call. Actual recv bytes is much smaller though. 
  !probably need to address this, but currently working without fixing this.
  
 ! problem = 0

  ! Post the asynchronous receives first:

    do taskmap_idx = 1, neighbor_mpi_cnt 

      recv_task_dbl = neighbor_mpi_map(taskmap_idx) ! created by Allan's code 
      recv_task_int = neighbor_mpi_map(taskmap_idx) ! created by Allan's code
      unpack_array_int(1,taskmap_idx) = recv_task_int
      unpack_array_dbl(1,taskmap_idx) = recv_task_dbl
       
        call mpi_irecv(unpack_array_int(3, taskmap_idx),proc_num_atms_min_bound*6, mpi_integer, &
                     recv_task_int, gifd_tag_int, pmemd_comm, &
                     recv_req_int(taskmap_idx), err_code_mpi); MPI_ERR_CHK("in reorient_distribution irecv int")
  
        call mpi_irecv(unpack_array_dbl(3,taskmap_idx), proc_num_atms_min_bound*14, mpi_double, &
                     recv_task_dbl, gifd_tag_dbl, pmemd_comm, &
                     recv_req_dbl(taskmap_idx), err_code_mpi);  MPI_ERR_CHK("in reorient_distribution irecv dbl")
    end do

  ! Now set up and post the asynchronous sends:

  do taskmap_idx = 1, neighbor_mpi_cnt 

    send_task_dbl = int(pack_array_dbl(1, taskmap_idx)) ! the mpi rank to send
    send_task_int = pack_array_int(1,taskmap_idx) ! the mpi rank to send

      send_atm_cnt = pack_array_int(2,taskmap_idx) ! # of integers to send 

      call mpi_isend(pack_array_int(3,taskmap_idx), send_atm_cnt, &
                     mpi_integer, send_task_int, gifd_tag_int, pmemd_comm, &
                     send_req_int(taskmap_idx), err_code_mpi); MPI_ERR_CHK("in reorient_distribution isend int")

      !sending the doubles
      send_atm_cnt = int(pack_array_dbl(2,taskmap_idx)) ! # of doubles to send

      call mpi_isend(pack_array_dbl(3,taskmap_idx), send_atm_cnt, &
                     mpi_double, send_task_dbl, gifd_tag_dbl, pmemd_comm, &
                     send_req_dbl(taskmap_idx), err_code_mpi); MPI_ERR_CHK("in reorient_distribution isend dbl")

!      write(25+mytaskid,*) "sending atoms", send_atm_cnt,"send_task_dbl ", send_task_dbl,"mytaskid", mytaskid

  end do

  ! Wait on and process the pending receive requests:

    do wait_call = 1, neighbor_mpi_cnt 
      !Now for integer
      call mpi_waitany(neighbor_mpi_cnt , recv_req_int, taskmap_idx, irecv_stat, &
                       err_code_mpi); MPI_ERR_CHK("in reorient_distribution wait recv int")
      !taskmap idx get the value index from which recv_task_dbl derived
      recv_task_int = neighbor_mpi_map(taskmap_idx)

      call mpi_get_count(irecv_stat, mpi_integer, recv_atm_cnt, err_code_mpi); MPI_ERR_CHK("in reorient_distribution wait recv dbl")

      unpack_array_int(2,taskmap_idx) = recv_atm_cnt !storing the count
     ! print *, "from recv waitcall int ", recv_task_int, taskmap_idx, mytaskid, recv_atm_cnt

      !Now for double
      call mpi_waitany(neighbor_mpi_cnt, recv_req_dbl, taskmap_idx, irecv_stat, &
                       err_code_mpi); MPI_ERR_CHK("in reorient_distributioni wait recv dbl")
      !taskmap idx find the index from which recv_task_dbl derived
      recv_task_dbl = neighbor_mpi_map(taskmap_idx)

      call mpi_get_count(irecv_stat, mpi_double, recv_atm_cnt, err_code_mpi); MPI_ERR_CHK("in reorient_distribution get recv count")

      !recv_atm_cnts_dbl(recv_task_dbl) = recv_atm_cnt ! may not required
      unpack_array_dbl(2,taskmap_idx) = dble(recv_atm_cnt) !storing the count
!       write(25+mytaskid,*) "received atoms", recv_atm_cnt,"recv_task_dbl ", recv_task_dbl,"mytaskid", mytaskid
    end do

  ! Wait for all sends to complete:
  call mpi_waitall(neighbor_mpi_cnt , send_req_int, isend_stat, err_code_mpi); MPI_ERR_CHK("in reorient_distribution wait send int")
  call mpi_waitall(neighbor_mpi_cnt , send_req_dbl, isend_stat, err_code_mpi); MPI_ERR_CHK("in reorient_distribution wait send dbl")

end subroutine reorient_distribution

!*******************************************************************************
!
! Subroutine:   proc_save_bonded
!
! Description: Saving bonds, angle and dihed global ids in arrays, later after atom_sorting
! and atm_distributiion it will be restored, so that the bonded arrays also sorted
!              
!*******************************************************************************
subroutine  proc_save_bonded(counts)
  use processor_mod
  use bonds_midpoint_mod
  use angles_midpoint_mod
  use dihedrals_midpoint_mod
  use parallel_dat_mod
  implicit none
  
  integer       :: counts
  integer       :: j
  !save the bonds index as global id

  do j = 1, cit_nbonh
    old_cit_h_bond(j)%atm_i = proc_atm_to_full_list(cit_h_bond(j)%atm_i)
    old_cit_h_bond(j)%atm_j = proc_atm_to_full_list(cit_h_bond(j)%atm_j)
  end do
  do j = 1, cit_nbona
    old_cit_a_bond(j)%atm_i = proc_atm_to_full_list(cit_a_bond(j)%atm_i)
    old_cit_a_bond(j)%atm_j = proc_atm_to_full_list(cit_a_bond(j)%atm_j)
  end do
  !save the angle index as global id
  do j = 1, cit_ntheth
    old_cit_h_angle(j)%atm_i = proc_atm_to_full_list(cit_h_angle(j)%atm_i)
    old_cit_h_angle(j)%atm_j = proc_atm_to_full_list(cit_h_angle(j)%atm_j)
    old_cit_h_angle(j)%atm_k = proc_atm_to_full_list(cit_h_angle(j)%atm_k)
  end do
  do j = 1, cit_ntheta
    old_cit_a_angle(j)%atm_i = proc_atm_to_full_list(cit_a_angle(j)%atm_i)
    old_cit_a_angle(j)%atm_j = proc_atm_to_full_list(cit_a_angle(j)%atm_j)
    old_cit_a_angle(j)%atm_k = proc_atm_to_full_list(cit_a_angle(j)%atm_k)
  end do
  !save the dihedrals index as global id
  do j = 1, cit_nphih
    old_cit_h_dihed(j)%atm_i = proc_atm_to_full_list(cit_h_dihed(j)%atm_i)
    old_cit_h_dihed(j)%atm_j = proc_atm_to_full_list(cit_h_dihed(j)%atm_j)
    old_cit_h_dihed(j)%atm_k =proc_atm_to_full_list(iabs(cit_h_dihed(j)%atm_k))*sign(1,cit_h_dihed(j)%atm_k)
    old_cit_h_dihed(j)%atm_l =proc_atm_to_full_list(iabs(cit_h_dihed(j)%atm_l))*sign(1,cit_h_dihed(j)%atm_l)
  end do
  do j = 1, cit_nphia
    old_cit_a_dihed(j)%atm_i = proc_atm_to_full_list(cit_a_dihed(j)%atm_i)
    old_cit_a_dihed(j)%atm_j = proc_atm_to_full_list(cit_a_dihed(j)%atm_j)
    old_cit_a_dihed(j)%atm_k = proc_atm_to_full_list(iabs(cit_a_dihed(j)%atm_k))*sign(1,cit_a_dihed(j)%atm_k)
    old_cit_a_dihed(j)%atm_l = proc_atm_to_full_list(iabs(cit_a_dihed(j)%atm_l))*sign(1,cit_a_dihed(j)%atm_l)
  end do
end subroutine proc_save_bonded

!*******************************************************************************
!
! Subroutine:  bonded_packing
!
! Description: pack bonds, angles and dihedrals that are escape from my domein
! and goes to another domain
!              
!*******************************************************************************
subroutine bonded_packing(crd, nbr_sdmn_cnt,counts,nstep)
  use pbc_mod
  use processor_mod
  use bonds_midpoint_mod
  use angles_midpoint_mod
  use dihedrals_midpoint_mod
  !use nb_pairlist_mod, only: del_cnt, atms_to_del, dest_sdmns
  use prmtop_dat_mod, only : natom
  
  implicit none
  
  ! Formal arguments
  ! Array of coordinates
  double precision :: crd(3, *)
  ! Number of neighboring sub-domains
  integer          :: nbr_sdmn_cnt
  integer          :: counts
  
  ! Local variables
  logical :: del_flag
  integer :: i, j, k, num_ints, atm_to_del
  integer :: b_lead_atm_h, b_lead_atm_a
  integer :: a_lead_atm_h, a_lead_atm_a
  integer :: d_lead_atm_h, d_lead_atm_a
  integer :: src_sdmn, dest_sdmn, src_recv_cnt, dest_send_cnt, src_idx, dest_idx
  !integer, dimension(2, 0:numtasks-1) :: b_send_map, a_send_map, d_send_map
  double Precision :: crd_x, crd_y, crd_z
  integer          :: atmid,src_mpi
  integer          :: new_cit_nbonh, new_cit_nbona
  integer          :: new_cit_ntheth, new_cit_ntheta
  integer          :: new_cit_nphih, new_cit_nphia
#ifdef BEXDBG
  integer          :: non_lead_atm_del_cnt=0, self_dest_cnt=0
  integer          :: b_send_cnt=0, a_send_cnt=0, d_send_cnt=0
  integer          :: fileunit
  character(len=32):: dbgfile
#endif
  integer         :: nstep
  
  !NOTE: No need of receieve map because recv_req arrays themselves act as a map
  ! Initialize the send map to zero
  b_send_map(:,:) = 0
  a_send_map(:,:) = 0
  d_send_map(:,:) = 0
    ! H type bonds
    new_cit_nbonh = cit_nbonh
    new_cit_nbona = cit_nbona
    new_cit_ntheth = cit_ntheth 
    new_cit_ntheta = cit_ntheta 
    new_cit_nphih = cit_nphih
    new_cit_nphia = cit_nphia
! if(mytaskid==0) write(0,*) "____________________"
    do j = cit_nbonh,1, -1
      b_lead_atm_h = cit_h_bond(j)%atm_i 
      do i = 1, del_cnt
        atm_to_del = atms_to_del(i)
        dest_sdmn = dest_sdmns(i) !the index i is consecutive based on how it was build 
        del_flag = .false.
        if (b_lead_atm_h .eq. atm_to_del) then
          ! Update the send_map for bonds
          ! send map's first row stores the index of neighbor sub-domain's
          ! rank in neighbor_mpi_map array; and second row stores the count of
          ! bonds to be sent to that sub-domain.
          do k = 1, nbr_sdmn_cnt
            if (neighbor_mpi_map(k) .eq. dest_sdmn) then
              b_send_map(1, dest_sdmn) = k
              b_send_map(2, dest_sdmn) = b_send_map(2, dest_sdmn) + 1
              dest_idx = b_send_map(1, dest_sdmn)
              dest_send_cnt = b_send_map(2, dest_sdmn)
              exit
            end if
          end do
!if(mytaskid==0) write(0,*) j, i,dest_send_cnt,k,dest_sdmn
          ! Pack the bond info in send buffer; convert local index to global
          b_send_buf(dest_send_cnt, dest_idx)%atm_i = proc_atm_to_full_list(cit_h_bond(j)%atm_i)
          b_send_buf(dest_send_cnt, dest_idx)%atm_j = proc_atm_to_full_list(cit_h_bond(j)%atm_j)
          b_send_buf(dest_send_cnt, dest_idx)%parm_idx = cit_h_bond(j)%parm_idx
          b_send_buf(dest_send_cnt, dest_idx)%hatype = 0

        ! Overwrite last bond record at location j (implies deletion)
          my_hbonds_leads(j) = my_hbonds_leads(new_cit_nbonh)
          cit_h_bond(j)%atm_i = cit_h_bond(new_cit_nbonh)%atm_i
          cit_h_bond(j)%atm_j = cit_h_bond(new_cit_nbonh)%atm_j
          cit_h_bond(j)%parm_idx = cit_h_bond(new_cit_nbonh)%parm_idx
          new_cit_nbonh = new_cit_nbonh - 1
          del_flag = .true.
          exit
        end if
      end do  ! i loop
    end do  ! j loop
    ! Non-H type bonds
!    if (.not. del_flag) then
    do j = cit_nbona, 1, -1
      b_lead_atm_a = cit_a_bond(j)%atm_i 
      do i = 1, del_cnt
        atm_to_del = atms_to_del(i)
        dest_sdmn = dest_sdmns(i) !the index i is consecutive based on how it was build 
        del_flag = .false.
        if (b_lead_atm_a .eq. atm_to_del) then
          ! Update the send_map for bonds
          do k = 1, nbr_sdmn_cnt
            if (neighbor_mpi_map(k) .eq. dest_sdmn) then
              b_send_map(1, dest_sdmn) = k
              b_send_map(2, dest_sdmn) = b_send_map(2, dest_sdmn) + 1
              dest_idx = b_send_map(1, dest_sdmn)
              dest_send_cnt = b_send_map(2, dest_sdmn)
              exit
            end if
          end do
          ! Pack the bond info in send buffer; convert local index to global
          b_send_buf(dest_send_cnt, dest_idx)%atm_i = proc_atm_to_full_list(cit_a_bond(j)%atm_i)
          b_send_buf(dest_send_cnt, dest_idx)%atm_j = proc_atm_to_full_list(cit_a_bond(j)%atm_j)
          b_send_buf(dest_send_cnt, dest_idx)%parm_idx = cit_a_bond(j)%parm_idx
          b_send_buf(dest_send_cnt, dest_idx)%hatype = 1
          ! Overwrite last bond record at location j (implies deletion)
          my_abonds_leads(j) = my_abonds_leads(new_cit_nbona)
          cit_a_bond(j)%atm_i = cit_a_bond(new_cit_nbona)%atm_i
          cit_a_bond(j)%atm_j = cit_a_bond(new_cit_nbona)%atm_j
          cit_a_bond(j)%parm_idx = cit_a_bond(new_cit_nbona)%parm_idx
          new_cit_nbona = new_cit_nbona - 1
          exit
        end if
      end do  ! i loop
  end do  ! j loop
    
    ! Packing of angles info in send buffer and deletion from local lists
    ! H type angles
  del_flag = .false.
  do j = cit_ntheth, 1, -1
    a_lead_atm_h = cit_h_angle(j)%atm_j
    do i = 1, del_cnt
      atm_to_del = atms_to_del(i)
      dest_sdmn = dest_sdmns(i) !the index i is consecutive based on how it was build 
      del_flag = .false.
      if (a_lead_atm_h .eq. atm_to_del) then
        ! Update the send_map for angles
        do k = 1, nbr_sdmn_cnt
          if (neighbor_mpi_map(k) .eq. dest_sdmn) then
            a_send_map(1, dest_sdmn) = k
            a_send_map(2, dest_sdmn) = a_send_map(2, dest_sdmn) + 1
            dest_idx = a_send_map(1, dest_sdmn) !the subscripts were switched
            dest_send_cnt = a_send_map(2, dest_sdmn)
            exit
          end if
        end do
        ! Pack the angle info in send buffer; convert local index to global
        a_send_buf(dest_send_cnt, dest_idx)%atm_i = proc_atm_to_full_list(cit_h_angle(j)%atm_i)
        a_send_buf(dest_send_cnt, dest_idx)%atm_j = proc_atm_to_full_list(cit_h_angle(j)%atm_j)
        a_send_buf(dest_send_cnt, dest_idx)%atm_k = proc_atm_to_full_list(cit_h_angle(j)%atm_k)
        a_send_buf(dest_send_cnt, dest_idx)%parm_idx = cit_h_angle(j)%parm_idx
        a_send_buf(dest_send_cnt, dest_idx)%hatype = 0
        ! Overwrite last angle record at location j (implies deletion)
        my_hangles_leads(j) = my_hangles_leads(new_cit_ntheth)
        cit_h_angle(j)%atm_i = cit_h_angle(new_cit_ntheth)%atm_i
        cit_h_angle(j)%atm_j = cit_h_angle(new_cit_ntheth)%atm_j
        cit_h_angle(j)%atm_k = cit_h_angle(new_cit_ntheth)%atm_k
        cit_h_angle(j)%parm_idx = cit_h_angle(new_cit_ntheth)%parm_idx
        new_cit_ntheth = new_cit_ntheth - 1
        del_flag = .true.
        exit
      end if
    end do  ! i loop
  end do  ! j loop
    ! Non-H type angles
!    if (.not. del_flag) then
  do j =  cit_ntheta, 1, -1
    a_lead_atm_a = cit_a_angle(j)%atm_j
    do i = 1, del_cnt
      atm_to_del = atms_to_del(i)
      dest_sdmn = dest_sdmns(i) !the index i is consecutive based on how it was build 
      del_flag = .false.
      if (a_lead_atm_a .eq. atm_to_del) then
        ! Update the send_map for angles
        do k = 1, nbr_sdmn_cnt
          if (neighbor_mpi_map(k) .eq. dest_sdmn) then
            a_send_map(1, dest_sdmn) = k
            a_send_map(2, dest_sdmn) = a_send_map(2, dest_sdmn) + 1
            dest_idx = a_send_map(1, dest_sdmn) !the subscripts were switched
            dest_send_cnt = a_send_map(2, dest_sdmn) !the subscripts were switched
            exit
          end if
        end do
        ! Pack the angle info in send buffer; convert local index to global
        a_send_buf(dest_send_cnt, dest_idx)%atm_i = proc_atm_to_full_list(cit_a_angle(j)%atm_i)
        a_send_buf(dest_send_cnt, dest_idx)%atm_j = proc_atm_to_full_list(cit_a_angle(j)%atm_j)
        a_send_buf(dest_send_cnt, dest_idx)%atm_k = proc_atm_to_full_list(cit_a_angle(j)%atm_k)
        a_send_buf(dest_send_cnt, dest_idx)%parm_idx = cit_a_angle(j)%parm_idx
        a_send_buf(dest_send_cnt, dest_idx)%hatype = 1
        ! Overwrite last angle record at location j (implies deletion)
        my_aangles_leads(j) = my_aangles_leads(new_cit_ntheta)
        cit_a_angle(j)%atm_i = cit_a_angle(new_cit_ntheta)%atm_i
        cit_a_angle(j)%atm_j = cit_a_angle(new_cit_ntheta)%atm_j
        cit_a_angle(j)%atm_k = cit_a_angle(new_cit_ntheta)%atm_k
        cit_a_angle(j)%parm_idx = cit_a_angle(new_cit_ntheta)%parm_idx
        new_cit_ntheta = new_cit_ntheta - 1
        exit
      end if
    end do  ! i loop
  end do  ! j loop
!    end if ! del_flag
    
    ! Packing of dihedrals info in send buffer and deletion from local lists
    ! H type dihedrals
    del_flag = .false.
  do j = cit_nphih, 1, -1
    d_lead_atm_h = cit_h_dihed(j)%atm_j 

    do i = 1, del_cnt
      atm_to_del = atms_to_del(i)
      dest_sdmn = dest_sdmns(i) !the index i is consecutive based on how it was build 
      del_flag = .false.
      if (d_lead_atm_h .eq. atm_to_del) then
        ! Update the send_map for dihedrals
        do k = 1, nbr_sdmn_cnt
          if (neighbor_mpi_map(k) .eq. dest_sdmn) then
            d_send_map(1, dest_sdmn) = k
            d_send_map(2, dest_sdmn) = d_send_map(2, dest_sdmn) + 1
            dest_idx = d_send_map(1, dest_sdmn)
            dest_send_cnt = d_send_map(2, dest_sdmn)
            exit
          end if
        end do
        ! Pack the dihedral info in send buffer; convert local index to global
        d_send_buf(dest_send_cnt, dest_idx)%atm_i = proc_atm_to_full_list(cit_h_dihed(j)%atm_i)
        d_send_buf(dest_send_cnt, dest_idx)%atm_j = proc_atm_to_full_list(cit_h_dihed(j)%atm_j)
        d_send_buf(dest_send_cnt, dest_idx)%atm_k = proc_atm_to_full_list(iabs(cit_h_dihed(j)%atm_k)) * sign(1,cit_h_dihed(j)%atm_k)
        d_send_buf(dest_send_cnt, dest_idx)%atm_l = proc_atm_to_full_list(iabs(cit_h_dihed(j)%atm_l)) * sign(1,cit_h_dihed(j)%atm_l)
        d_send_buf(dest_send_cnt, dest_idx)%parm_idx = cit_h_dihed(j)%parm_idx
        d_send_buf(dest_send_cnt, dest_idx)%hatype = 0
        ! Overwrite last dihedral record at location j (implies deletion)
        my_hdiheds_leads(j) = my_hdiheds_leads(new_cit_nphih)
        cit_h_dihed(j)%atm_i = cit_h_dihed(new_cit_nphih)%atm_i
        cit_h_dihed(j)%atm_j = cit_h_dihed(new_cit_nphih)%atm_j
        cit_h_dihed(j)%atm_k = cit_h_dihed(new_cit_nphih)%atm_k
        cit_h_dihed(j)%atm_l = cit_h_dihed(new_cit_nphih)%atm_l
        cit_h_dihed(j)%parm_idx = cit_h_dihed(new_cit_nphih)%parm_idx
        new_cit_nphih = new_cit_nphih - 1
        del_flag = .true.
        exit
      end if
    end do  ! i loop
  end do  ! j loop
    
    ! Non-H type dihedrals
!    if (.not. del_flag) then
  do j = cit_nphia, 1, -1
    d_lead_atm_a = cit_a_dihed(j)%atm_j 
    do i = 1, del_cnt
      atm_to_del = atms_to_del(i)
      dest_sdmn = dest_sdmns(i) !the index i is consecutive based on how it was build 
      del_flag = .false.
        if (d_lead_atm_a .eq. atm_to_del) then
          ! Update the send_map for dihedrals
          do k = 1, nbr_sdmn_cnt
            if (neighbor_mpi_map(k) .eq. dest_sdmn) then
              d_send_map(1, dest_sdmn) = k
              d_send_map(2, dest_sdmn) = d_send_map(2, dest_sdmn) + 1
              dest_idx = d_send_map(1, dest_sdmn)
              dest_send_cnt = d_send_map(2, dest_sdmn)
              exit
            end if
          end do
          ! Pack the dihedral info in send buffer; convert local index to global
          d_send_buf(dest_send_cnt, dest_idx)%atm_i = proc_atm_to_full_list(cit_a_dihed(j)%atm_i)
          d_send_buf(dest_send_cnt, dest_idx)%atm_j = proc_atm_to_full_list(cit_a_dihed(j)%atm_j)
          d_send_buf(dest_send_cnt, dest_idx)%atm_k = proc_atm_to_full_list(iabs(cit_a_dihed(j)%atm_k)) &
                                                      * sign(1,cit_a_dihed(j)%atm_k)
          d_send_buf(dest_send_cnt, dest_idx)%atm_l = proc_atm_to_full_list(iabs(cit_a_dihed(j)%atm_l)) &
                                                      * sign(1,cit_a_dihed(j)%atm_l)
          d_send_buf(dest_send_cnt, dest_idx)%parm_idx = cit_a_dihed(j)%parm_idx
          d_send_buf(dest_send_cnt, dest_idx)%hatype = 1
          ! Overwrite last dihedral record at location j (implies deletion)
          my_adiheds_leads(j) = my_adiheds_leads(new_cit_nphia)
          cit_a_dihed(j)%atm_i = cit_a_dihed(new_cit_nphia)%atm_i
          cit_a_dihed(j)%atm_j = cit_a_dihed(new_cit_nphia)%atm_j
          cit_a_dihed(j)%atm_k = cit_a_dihed(new_cit_nphia)%atm_k
          cit_a_dihed(j)%atm_l = cit_a_dihed(new_cit_nphia)%atm_l
          cit_a_dihed(j)%parm_idx = cit_a_dihed(new_cit_nphia)%parm_idx
          new_cit_nphia = new_cit_nphia - 1
          exit
        end if
  !  end if ! del_flag
    end do  ! i loop
  end do  ! j loop
    cit_nbonh = new_cit_nbonh
    cit_nbona = new_cit_nbona
    cit_ntheth = new_cit_ntheth 
    cit_ntheta = new_cit_ntheta 
    cit_nphih = new_cit_nphih
    cit_nphia = new_cit_nphia
  return

end subroutine bonded_packing



#endif /*MPI*/

end module parallel_processor_mod
