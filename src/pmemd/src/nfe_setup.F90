! Pre-processing
#ifndef NFE_UTILS_H
#define NFE_UTILS_H

#ifndef NFE_DISABLE_ASSERT
#  define nfe_assert(stmt) if (.not.(stmt)) call afailed(__FILE__, __LINE__)
#  define nfe_assert_not_reached() call afailed(__FILE__, __LINE__)
#  define NFE_PURE_EXCEPT_ASSERT
#  define NFE_USE_AFAILED use nfe_lib_mod, only : afailed
#else
#  define nfe_assert(s) 
#  define nfe_assert_not_reached() 
#  define NFE_PURE_EXCEPT_ASSERT pure
#  define NFE_USE_AFAILED
#endif /* NFE_DISABLE_ASSERT */

#define NFE_OUT_OF_MEMORY call out_of_memory(__FILE__, __LINE__)

#ifdef MPI
#  define NFE_MASTER_ONLY_BEGIN if (mytaskid.eq.0) then
#  define NFE_MASTER_ONLY_END end if
#else
#  define NFE_MASTER_ONLY_BEGIN
#  define NFE_MASTER_ONLY_END
#endif /* MPI */

#define NFE_ERROR   ' ** NFE-Error ** : '
#define NFE_WARNING ' ** NFE-Warning ** : '
#define NFE_INFO    ' NFE : '

#endif /* NFE_UTILS_H */
!----------------------------------------------------------------

module nfe_setup_mod

! old REMD subroutines
! currently not used, only for records
#ifdef MPI
public :: on_delta
public :: on_exchange
#endif /* MPI */

! old subroutine, not used
public :: on_multipmemd_exit

! pmemd.F90
public :: on_pmemd_init
public :: on_pmemd_exit

! gb/pme_force.F90
public :: on_force

! runmd.F90, write biasing potential
public :: on_mdwrit

! runmd.F90
#ifdef MPI
public :: on_mdstep
#endif /* MPI */

! ---------------------------------------------------------------

contains

! ---------------------------------------------------------------

#ifdef MPI
subroutine on_delta(o_master_rank, need_U_xx, U_mm, U_mo, U_om, U_oo)

   use nfe_lib_mod, only : ZERO
   use nfe_pmd_mod, only : pmd_on_delta => on_delta
   use nfe_abmd_mod, only : abmd_on_delta => on_delta

   implicit none

   integer, intent(in) :: o_master_rank
   logical, intent(in) :: need_U_xx

   double precision, intent(out) :: U_mm, U_mo, U_om, U_oo

   U_mm = ZERO
   U_mo = ZERO
   U_om = ZERO
   U_oo = ZERO

   call pmd_on_delta(o_master_rank, need_U_xx, U_mm, U_mo, U_om, U_oo)
   call abmd_on_delta(o_master_rank, need_U_xx, U_mm, U_mo, U_om, U_oo)

end subroutine on_delta

!-----------------------------------------------------------------------------

subroutine on_exchange(o_master_rank)

   use nfe_pmd_mod, only : pmd_on_exchange => on_exchange
   use nfe_abmd_mod, only : abmd_on_exchange => on_exchange

   implicit none

   integer, intent(in) :: o_master_rank

   call pmd_on_exchange(o_master_rank)
   call abmd_on_exchange(o_master_rank)

end subroutine on_exchange
#endif /* MPI */

!-----------------------------------------------------------------------------

subroutine on_multipmemd_exit()

   use parallel_dat_mod
   use nfe_pmd_mod, only : pmd_on_multipmemd_exit => on_multipmemd_exit
   use nfe_abmd_mod, only : abmd_on_multipmemd_exit => on_multipmemd_exit
   use nfe_lib_mod

   implicit none

   call pmd_on_multipmemd_exit()
   call abmd_on_multipmemd_exit()

   NFE_MASTER_ONLY_BEGIN
   call proxy_finalize()
   NFE_MASTER_ONLY_END

end subroutine on_multipmemd_exit

!-----------------------------------------------------------------------------

!
! 'on_pmemd_init()' is called at the point when
! MPI is initialized and MDIN/PRMTOP/INPCRD are loaded
!

subroutine on_pmemd_init(ih, amass, acrds, rem)

   use file_io_dat_mod, only :mdin
   use parallel_dat_mod
   use nfe_lib_mod
   use nfe_smd_mod, only : smd_on_pmemd_init => on_pmemd_init
   use nfe_pmd_mod, only : pmd_on_pmemd_init => on_pmemd_init
   use nfe_abmd_mod, only : abmd_on_pmemd_init => on_pmemd_init
   use nfe_stsm_mod, only : stsm_on_pmemd_init => on_pmemd_init
#ifdef MPI
   use nfe_bbmd_mod, only : bbmd_on_pmemd_init => on_pmemd_init, bbmd_active => active
#endif /* MPI */

   implicit none

   character(len = 4), intent(in) :: ih(*)

   double precision, intent(in) :: amass(*)
   double precision, intent(in) :: acrds(*)

   integer, intent(in) :: rem

   integer, save :: my_unit = mdin

#ifdef MPI
   integer :: total_nfe_atm_cnt, err, i
   integer, allocatable :: all_nfe_atm_cnt(:), all_nfe_atm_lst(:), displs(:)
#endif

   call remember_rem(rem)

   NFE_MASTER_ONLY_BEGIN
   call remember_atom_names(ih)
   NFE_MASTER_ONLY_END

#ifdef MPI
   if (multipmemd_rem().eq.0) &
      call smd_on_pmemd_init(my_unit, amass, acrds)
#else
   call smd_on_pmemd_init(my_unit, amass, acrds)
#endif

   call pmd_on_pmemd_init(my_unit, amass)
   call abmd_on_pmemd_init(my_unit, amass)
   call stsm_on_pmemd_init(my_unit, amass)

#ifdef MPI      
   if (multipmemd_rem().eq.0.and.multipmemd_numgroup().gt.1) &
      call bbmd_on_pmemd_init(my_unit, amass)
#endif 

#ifdef CUDA
#ifdef MPI
   ! allgather the atm_lst if there is replica exchange
   NFE_MASTER_ONLY_BEGIN
   if (bbmd_active.eq.1 .or. multipmemd_rem().ne.0) then
      allocate(all_nfe_atm_cnt(master_size))
      call mpi_allgather(nfe_atm_cnt, 1, MPI_INTEGER, all_nfe_atm_cnt, 1, MPI_INTEGER, pmemd_master_comm, err)
      nfe_assert(err.eq.0)
      call mpi_allreduce(nfe_atm_cnt, total_nfe_atm_cnt, 1, MPI_INTEGER, MPI_SUM, pmemd_master_comm, err)
      nfe_assert(err.eq.0)

      allocate(all_nfe_atm_lst(total_nfe_atm_cnt), stat=err)

      allocate(displs(master_size))
      displs(1) = 0
      do i=2,master_size
         displs(i) = displs(i-1) + all_nfe_atm_cnt(i-1)
      end do

      call mpi_allgatherv(nfe_atm_lst, nfe_atm_cnt, MPI_INTEGER, all_nfe_atm_lst, all_nfe_atm_cnt, displs, MPI_INTEGER, &
                         pmemd_master_comm, err)
      nfe_assert(err.eq.0)

      nfe_atm_cnt = total_nfe_atm_cnt
      deallocate(nfe_atm_lst)
      call move_alloc(all_nfe_atm_lst, nfe_atm_lst)
   end if
   NFE_MASTER_ONLY_END
#endif
   if (allocated(nfe_atm_lst)) then
      call bubble_sort(nfe_atm_lst,nfe_atm_cnt)
      call remove_zero_duplicate(nfe_atm_lst,nfe_atm_cnt)
      nfe_atm_cnt = size(nfe_atm_lst)
   end if

contains
subroutine bubble_sort(a,n)
  implicit none
  integer :: n,a(:)
  integer i,j,temp
  do i=n-1,1,-1
    do j=1,i
      if ( a(j) > a(j+1) ) then
        temp=a(j)
        a(j)=a(j+1)
        a(j+1)=temp
      end if
    end do
  end do
  return
end subroutine

subroutine remove_zero_duplicate(a,n)
  implicit none
  integer :: n
  integer, allocatable :: tmp(:)
  integer, allocatable, intent(inout) :: a(:)
  integer :: prev = 0, i

  do i=1,n
    if (a(i).ne.0 .and. a(i).ne.prev) then
       call AddToList(tmp, a(i))
       prev = a(i)
    end if
  end do
  deallocate(a)
  call move_alloc(tmp, a)
end subroutine
#endif

end subroutine on_pmemd_init

!-----------------------------------------------------------------------------

subroutine on_pmemd_exit()

   use nfe_smd_mod, only : smd_on_pmemd_exit => on_pmemd_exit
   use nfe_pmd_mod, only : pmd_on_pmemd_exit => on_pmemd_exit
   use nfe_abmd_mod, only : abmd_on_pmemd_exit => on_pmemd_exit
   use nfe_stsm_mod, only : stsm_on_pmemd_exit => on_pmemd_exit
#ifdef MPI
   use nfe_bbmd_mod, only : bbmd_on_pmemd_exit => on_pmemd_exit
#endif /* MPI */

   implicit none

   call smd_on_pmemd_exit()
   call pmd_on_pmemd_exit()
   call abmd_on_pmemd_exit()
   call stsm_on_pmemd_exit()

#ifdef MPI
   call bbmd_on_pmemd_exit()
#endif /* MPI */

end subroutine on_pmemd_exit

!-----------------------------------------------------------------------------

subroutine on_force(x, f, pot)

   use nfe_lib_mod
   use nfe_smd_mod, only : smd_on_force => on_force
   use nfe_pmd_mod, only : pmd_on_force => on_force
   use nfe_abmd_mod, only : abmd_on_force => on_force
   use nfe_stsm_mod, only : stsm_on_force => on_force
#ifdef MPI
   use nfe_bbmd_mod, only : bbmd_on_force => on_force
#endif /* MPI */

   implicit none

   double precision, intent(in) :: x(3,pmemd_natoms())
   double precision, intent(inout) :: f(3,pmemd_natoms())
   double precision, intent(inout) :: pot
! Modified by M Moradi
! for driven ABMD
   double precision :: wdriven = ZERO
   double precision :: udriven = ZERO
! Moradi end
   
   nfe_pot_ene = null_nfe_pot_ene_rec

   call smd_on_force(x, f, wdriven,udriven, nfe_pot_ene%smd)
   call pmd_on_force(x, f, nfe_pot_ene%pmd)
   call abmd_on_force(x, f, wdriven,udriven, nfe_pot_ene%abmd)
   call stsm_on_force(x, f, nfe_pot_ene%stsm)
#ifdef MPI
   call bbmd_on_force(x, f, wdriven, udriven, nfe_pot_ene%bbmd)
#endif /* MPI */
   
   nfe_pot_ene%total = nfe_pot_ene%smd + nfe_pot_ene%pmd + nfe_pot_ene%abmd + &
                       nfe_pot_ene%bbmd + nfe_pot_ene%stsm
   pot = nfe_pot_ene%total

   nfe_real_mdstep = .true.  ! reset real_mdstep to true

end subroutine on_force

!-----------------------------------------------------------------------------

subroutine on_mdwrit()

   use nfe_abmd_mod, only : abmd_on_mdwrit => on_mdwrit
#ifdef MPI
   use nfe_bbmd_mod, only : bbmd_on_mdwrit => on_mdwrit
#endif /* MPI */

   implicit none

   call abmd_on_mdwrit()
#ifdef MPI
   call bbmd_on_mdwrit()
#endif /* MPI */

end subroutine on_mdwrit

!-----------------------------------------------------------------------------
! Modified by M Moradi
! for selection algorithm (runmd.F90 was modified accordingly as well)
#ifdef MPI
subroutine on_mdstep(eptot, x, v)

   use nfe_lib_mod
   use nfe_bbmd_mod, only : bbmd_on_mdstep => on_mdstep
   use nfe_abmd_mod, only : abmd_on_mdstep => on_mdstep

   implicit none

   double precision, intent(in) :: eptot
   double precision, intent(inout) :: x(3,pmemd_natoms())
   double precision, intent(inout) :: v(3,pmemd_natoms())
   double precision :: trs_x(3*pmemd_natoms()), trs_v(3*pmemd_natoms())

   trs_x = reshape(x,shape(trs_x))
   trs_v = reshape(v,shape(trs_v))
   call bbmd_on_mdstep(eptot, trs_v)
   call abmd_on_mdstep(trs_x, trs_v)
   v = reshape(trs_v,shape(v))
   x = reshape(trs_x,shape(x))

end subroutine on_mdstep
#endif /* MPI */
! Moradi end

end module nfe_setup_mod
