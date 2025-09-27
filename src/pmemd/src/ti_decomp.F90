#ifdef TIDECOMP
#include "copyright.i"

!*******************************************************************************
!
! Module: ti_decomp_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module ti_decomp_mod

  implicit none

  ! Atomic Decomposition
  double precision, allocatable :: ti_decomp(:)
  integer, allocatable :: pme_cell_cnt(:,:,:,:)
  integer, allocatable :: pme_cell_atom(:,:,:,:)
  integer, allocatable :: pme_cell_weight(:,:,:,:)
  integer, private     :: total_atom
  integer              :: decomp_region
  integer, allocatable :: outputmask(:)
  integer :: protstart, protend,cofacstart,cofacend

  ! Ligand Decomposition
  double precision :: lig_dvdl_decomp(30)
  integer, allocatable :: id_atom_region(:)
  double precision, allocatable :: lig_ti_decomp(:,:)

  contains

!*******************************************************************************
!
! Subroutine: init_region_decomp
!
! Description: <TBS>
!
!*******************************************************************************

  subroutine init_region_decomp(atm_cnt, nres, igraph, isymbl, res_atms, labres, crd)

      use mdin_ctrl_dat_mod
      use file_io_dat_mod
      use findmask_mod
#if defined(CUDA)
      use pmemd_lib_mod !I think we just need setup_alloc_error
#endif

      implicit none

      ! Formal arguments:
      integer, intent(in)             :: atm_cnt, nres
      integer, intent(in)             :: res_atms(nres)
      character(len=4), intent(in)    :: igraph(atm_cnt)
      character(len=4), intent(in)    :: isymbl(atm_cnt)
      character(len=4), intent(in)    :: labres(nres)
      double precision, intent(in)    :: crd(3 * atm_cnt)
       
      integer :: i
#if defined(CUDA)
      integer                         :: alloc_failed
#endif
      integer :: cofac_region(atm_cnt)
      integer :: protein_region(atm_cnt)
      integer :: lig_region(atm_cnt)
      
      cofac_region(:)=0
      protein_region(:)=0
      lig_region(:)=0 

      allocate(id_atom_region(atm_cnt))

      if(ligmask .ne. '') then
        call atommask(atm_cnt, nres, 0, igraph, isymbl, res_atms, labres, &
               crd, ligmask, lig_region)
      end if
      if(cofacmask .ne. '') then
        call atommask(atm_cnt, nres, 0, igraph, isymbl, res_atms, labres, &
               crd, cofacmask, cofac_region)
      end if
      if(proteinmask .ne. '') then
        call atommask(atm_cnt, nres, 0, igraph, isymbl, res_atms, labres, &
               crd, proteinmask, protein_region)
      end if

      ! For ligand decomposition we use 30 because we identify the regions using
      ! a set of n numbers where the subset of any sum of 1 or 2 numbers in set n
      ! is unique.
      ! Ligand - 1
      ! Protein - 3
      ! Cofactor - 7
      ! Water - 15
      ! Then ligand-ligand interactions are 2, lig-protein 4, lig-cofactor-8, lig-water-16
      protstart=0
      protend=0

      do i=1,atm_cnt
          id_atom_region(i) = 7*cofac_region(i) + lig_region(i) + 3*protein_region(i)
          if(id_atom_region(i) .eq. 0) id_atom_region(i)=15
          if(protstart .eq. 0 .and. protein_region(i) .eq. 1) protstart=i
          if(protein_region(i) .eq. 1) protend=i
          if(cofacstart .eq. 0 .and. cofac_region(i) .eq. 1) cofacstart=i
          if(cofac_region(i) .eq. 1) cofacend=i
      end do 

  end subroutine init_region_decomp

!*******************************************************************************
!
! Subroutine: init_ti_decomp
!
! Description: <TBS>
!              
!*******************************************************************************

  subroutine init_ti_decomp(natom, output)
 
      use file_io_mod
      use ti_mod

      implicit none
   
      integer,intent(in)        :: natom
      integer,intent(in)        :: output
      integer :: i

      allocate(ti_decomp(natom))
      allocate(outputmask(natom))
      allocate(lig_ti_decomp(30,natom))
      total_atom = natom
      decomp_region = 1

      ti_decomp(:)=0.0d0
      lig_ti_decomp(:,:)=0.0d0

      if(output .gt. 0) call amopen(decomp_unit,decomp_name,'R', 'F', 'W')
 
  end subroutine init_ti_decomp

!*******************************************************************************
!
! Subroutine:  init_ti_decomp_mask
!
! Description: Fill in the ti_lst array, based on the amber mask.
!
!*******************************************************************************

subroutine init_ti_decomp_mask(atm_cnt, nres, igraph, isymbl, res_atms, labres, &
                            crd)
  use file_io_dat_mod
  use findmask_mod
#if defined(CUDA)
  use pmemd_lib_mod !I think we just need setup_alloc_error
#endif
  use mdin_ctrl_dat_mod, only : decompmask
  implicit none

! Formal arguments:
  integer, intent(in)             :: atm_cnt, nres
  integer, intent(in)             :: res_atms(nres)
  character(len=4), intent(in)    :: igraph(atm_cnt)
  character(len=4), intent(in)    :: isymbl(atm_cnt)
  character(len=4), intent(in)    :: labres(nres)
  double precision, intent(in)    :: crd(3 * atm_cnt)

#if defined(CUDA)
  integer                         :: alloc_failed
#endif

       outputmask(:)=1
       if(decompmask .ne. '') then
         call atommask(atm_cnt, nres, 0, igraph, isymbl, res_atms, labres, &
                crd, decompmask, outputmask)
       end if

end subroutine init_ti_decomp_mask

!*******************************************************************************
!
! Subroutine: init_pme_decomp_bookkeeping
!
! Description: <TBS>
!              
!*******************************************************************************

  subroutine init_pme_decomp_bookkeeping

      use mdin_ewald_dat_mod
   
      implicit none

      ! xdim, ydim, zdim 1 stores cnt
      ! xdim, ydim, zdim 2 stores total parts of the order
      allocate(pme_cell_cnt(nfft1,nfft2,nfft3,2))
      pme_cell_cnt(:,:,:,:) = 0
      ! xdim, ydim, zdim, 1..n stores atoms in xdim, ydim, zdim cell
      allocate(pme_cell_atom(nfft1,nfft2,nfft3,(bspl_order*2+1)**3))
      ! xdim, ydim, zdim, 1..n stores weights.
      ! Divide weight/cnt for contribution.
      allocate(pme_cell_weight(nfft1,nfft2,nfft3,(bspl_order*2+1)**3))

  end subroutine init_pme_decomp_bookkeeping

!*******************************************************************************
!
! Subroutine: zero_decomp_energies
!
! Description: zeroes decomp arrays if they're allocated
!              
!*******************************************************************************

  subroutine zero_decomp_energies

      implicit none

      if(allocated(pme_cell_cnt))pme_cell_cnt(:,:,:,:)=0
      if(allocated(pme_cell_atom))pme_cell_atom(:,:,:,:)=0
      if(allocated(pme_cell_weight))pme_cell_weight(:,:,:,:)=0

!      lig_dvdl_decomp(:)=0

  end subroutine zero_decomp_energies

!*******************************************************************************
!
! Subroutine: mddecomp
!
! Description: Writes data to output file
!              
!*******************************************************************************

  subroutine mddecomp(total_nstep)

      use file_io_dat_mod
      use ti_mod, only : ti_mode, ti_ti_atm_cnt, ti_latm_cnt

      implicit none

      integer, intent(in)         :: total_nstep
      integer                     :: i

      write(decomp_unit,*) "Steps ran:", total_nstep

      ! Atomic Decomp
      if(ti_mode .gt. 0) then
          !Region Decomp
          write(decomp_unit,*)"Ligand-Ligand DV/DL", lig_dvdl_decomp(2)+lig_dvdl_decomp(1)
          write(decomp_unit,*)"Ligand-Protein DV/DL", lig_dvdl_decomp(4)+lig_dvdl_decomp(3)
          write(decomp_unit,*)"Ligand-Cofactor DV/DL", lig_dvdl_decomp(8)+lig_dvdl_decomp(7)
          write(decomp_unit,*)"Ligand-Water/Bulk DV/DL", lig_dvdl_decomp(16)+lig_dvdl_decomp(15)+&
                  lig_dvdl_decomp(30)+lig_dvdl_decomp(18)  !Water-water / Water-protein
          write(decomp_unit,*)"Atom   ","LigDecomp   ","ProtDecomp   ","WatDecomp   ","DV/DL"
      end if

      do i=1, total_atom
        if(outputmask(i) .eq. 1) then
          if(ti_mode .gt. 0) then
              write(decomp_unit,*)i,' ', lig_ti_decomp(2,i)+lig_ti_decomp(1,i), &
                                         lig_ti_decomp(4,i)+lig_ti_decomp(3,i), &
                                         lig_ti_decomp(16,i)+lig_ti_decomp(15,i), &
                                         ti_decomp(i)
          end if
        end if
      end do

  end subroutine mddecomp

end module
#endif
