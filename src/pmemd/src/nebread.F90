#include "copyright.i"

module nebread_mod

  use prmtop_dat_mod
  use inpcrd_dat_mod, only : atm_crd
  use mdin_ctrl_dat_mod

  implicit none

  type :: atm_info
    integer :: atm_num     !atom number
    integer :: atm_idx     !atom position index
  end type atm_info

  double precision, allocatable :: tangents(:,:), springforce(:,:), neb_force(:,:)
  double precision, allocatable :: xprev(:,:), xnext(:,:), neb_nrg_all(:)
  integer :: neb_nbead, last_neb_atom, next_node, prev_node
  integer :: beadid  !mybeadid in sander neb
  integer, allocatable :: fitmask(:), rmsmask(:)
  type(atm_info), allocatable :: sortfit(:), sortrms(:)
  type(atm_info), allocatable :: combined_mask(:)   !DG: combined and sorted mask from fitmask and rmsmask
  integer, allocatable :: partial_mask(:)    !DG: same as combined_mask with no end zeroes, needed for partial coordinates download from gpu.
  integer, allocatable :: idx_mask(:)
  integer :: pcnt                            !DG: partial counts, partial_mask number of elements.
  integer :: alloc_stat

contains

subroutine nebread()

   use findmask_mod,      only : atommask
   use pmemd_lib_mod,     only : mexit
   use parallel_dat_mod
   implicit none

   integer :: inerr, n

#if !defined(MPI)
   if(ineb>0) then
     write(6,*) 'NEB requires MPI'
     inerr=1
   end if
#endif

   if (ineb>0) then
      last_neb_atom = 0

      if (ntr.ne.0) then
        write(6,'(/2x,a)') 'cannot use NEB with ntr restraints'
        call mexit(6,1)
      else
         ! read in atom group for fitting (=overlap region)
         call atommask( natom, nres, 0, atm_igraph, atm_isymbl, &
            gbl_res_atms, gbl_labres, atm_crd, tgtfitmask, fitmask )
         ! see comments above (for ntr) for the following reduction cycle
         nattgtfit = 0 ! number of atoms for tgtmd fitting (=overlap)
         do n=1,natom
           if( fitmask(n) <= 0 ) cycle
           nattgtfit = nattgtfit + 1
           fitmask(nattgtfit) = n
           if (n.gt.last_neb_atom) last_neb_atom = n
         end do
         write(6,'(a)') "The following selection will be used for NEB structure fitting"
         write(6,'(a,a,a,i5,a)')  &
         '     Mask "', tgtfitmask(1:len_trim(tgtfitmask)-1),  &
         '" matches ',nattgtfit,' atoms'
      end if
      ! read in atom group for tgtmd rmsd calculation
      call atommask( natom, nres, 0, atm_igraph, atm_isymbl, &
         gbl_res_atms, gbl_labres, atm_crd, tgtrmsmask, rmsmask )
      nattgtrms = 0 ! number of atoms for tgtmd rmsd calculation
      do n=1,natom
        if( rmsmask(n) <= 0 ) cycle
        nattgtrms = nattgtrms + 1
        rmsmask(nattgtrms) = n

        if (n.gt.last_neb_atom) last_neb_atom = n
      end do
      write(6,'(a)') "The following selection will be used for NEB force application"
      write(6,'(a,a,a,i5,a)')  &
      '     Mask "', tgtrmsmask(1:len_trim(tgtrmsmask)-1),  &
      '" matches ',nattgtrms,' atoms'
      write(6,'(/2x,a,i6)') "Last atom in NEB fitmask or rmsmask is ",last_neb_atom

      if (nattgtrms<=0 .or. nattgtfit <= 0) then
        write(6,'(/2x,a)') 'NEB requires use of tgtfitmask and tgtrmsmask'
        call mexit(6,1)
      endif

      allocate(sortfit(nattgtfit), stat=alloc_stat)
      if ( alloc_stat .ne. 0 ) &
        call allocation_error('nebread','Error allocating temporary arrays')
      sortfit(:)%atm_idx = 0

      allocate(sortrms(nattgtrms), stat=alloc_stat)
      if ( alloc_stat .ne. 0 ) &
        call allocation_error('nebread','Error allocating temporary arrays')
      sortrms(:)%atm_idx = 0

      allocate(combined_mask(nattgtrms+nattgtfit), stat=alloc_stat)
      if ( alloc_stat .ne. 0 ) &
        call allocation_error('nebread','Error allocating temporary arrays')
      combined_mask(:)%atm_num = 0
      combined_mask(:)%atm_idx = 0

      do n=1,nattgtfit
        sortfit(n)%atm_num = fitmask(n)
      end do

      do n=1,nattgtrms
        sortrms(n)%atm_num = rmsmask(n)
      end do

       call combine_sort(combined_mask, sortfit, sortrms, nattgtfit, nattgtrms)

       pcnt=0
       do n=1,size(combined_mask)
         if(combined_mask(n)%atm_num .ne. 0) then
           pcnt = pcnt + 1
         endif
       enddo

       allocate(partial_mask(pcnt), stat=alloc_stat)
       if ( alloc_stat .ne. 0 ) &
        call allocation_error('nebread','Error allocating temporary arrays')
       partial_mask(:) = 0

       allocate(idx_mask(pcnt), stat=alloc_stat)
       if ( alloc_stat .ne. 0 ) &
        call allocation_error('nebread','Error allocating temporary arrays')
       idx_mask(:) = 0

       partial_mask(1:pcnt) = combined_mask(1:pcnt)%atm_num
       idx_mask(1:pcnt) = combined_mask(1:pcnt)%atm_idx
   endif
end subroutine nebread
!-------------------------------------------------------------------------
subroutine combine_sort(combined, sortfit, sortrms, nattgtfit, nattgtrms)

  implicit none

  !Passed variables
  type(atm_info) :: combined(:)
  type(atm_info) :: sortfit(:), sortrms(:)
  integer :: nattgtfit, nattgtrms

  !Local variables
  integer :: fitcnt, rmscnt, cntr
  integer :: i, m

  fitcnt = 1
  rmscnt = 1
  cntr = 0

  do i=1,nattgtfit+nattgtrms
    if(sortfit(fitcnt)%atm_num .lt. sortrms(rmscnt)%atm_num) then
      combined(i)%atm_num = sortfit(fitcnt)%atm_num
      sortfit(fitcnt)%atm_idx = i
      combined(i)%atm_idx = 1  !atom belongs to sortfit
      fitcnt = fitcnt + 1
    else
      combined(i)%atm_num = sortrms(rmscnt)%atm_num
      sortrms(rmscnt)%atm_idx = i !leave combined(i)%atm_idx = 0, atom belongs to sortrms
      if(sortfit(fitcnt)%atm_num .eq. sortrms(rmscnt)%atm_num) then
        sortfit(fitcnt)%atm_idx = i
        combined(i)%atm_idx = 2 !atom belongs to both sortfit and sortrms
        fitcnt = fitcnt + 1
      endif
      rmscnt = rmscnt + 1
    endif
    cntr = cntr + 1
    if(fitcnt .gt. nattgtfit .or. rmscnt .gt. nattgtrms) exit
  enddo

  if(rmscnt .gt. nattgtrms) then
    m = nattgtfit - fitcnt
    do i=0,m
      combined(cntr + i + 1)%atm_num = sortfit(fitcnt + i)%atm_num
      sortfit(fitcnt +i)%atm_idx = cntr + i + 1
      combined(cntr + i + 1)%atm_idx = 1  !atom belongs to sortfit
    enddo
  else if(fitcnt .gt. nattgtfit) then
    m = nattgtrms - rmscnt
    do i=0,m
      combined(cntr + i + 1)%atm_num = sortrms(rmscnt + i)%atm_num
      sortrms(rmscnt + i)%atm_idx = cntr + i + 1 !leave combined(i)%atm_idx = 0, atom belongs to sortrms
    enddo
  endif

end subroutine combine_sort
!-------------------------------------------------------------------------
subroutine allocation_error(routine,message)

  use gbl_constants_mod, only : error_hdr, extra_line_hdr
  use pmemd_lib_mod,     only : mexit

  implicit none

  ! passed variables

  character(*), intent(in) :: routine, message

  write(mdout, '(a,a,a)') error_hdr, 'Error in ', routine
  write(mdout, '(a,a)') extra_line_hdr, message

  call mexit(mdout, 1)

end subroutine allocation_error

end module nebread_mod
