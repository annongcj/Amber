! Another module needed to avoid cyclic dependencies

module constantph_dat_mod

  implicit none

  public

  logical, save :: on_cpstep = .false.
  double precision, allocatable, save :: proposed_qterm(:)
  double precision, allocatable, save :: cnstph_frc(:,:)

contains

subroutine allocate_cnstph_dat(natom, num_reals, ierror)

  implicit none

  integer, intent(in)     :: natom
  integer, intent(in out) :: num_reals
  integer, intent(out)    :: ierror

  ierror = 0
  allocate(proposed_qterm(natom), cnstph_frc(3, natom), &
           stat=ierror)

  num_reals = num_reals + size(proposed_qterm) + size(cnstph_frc)

  return

end subroutine allocate_cnstph_dat

subroutine cleanup_cnstph_dat

  implicit none

  if (allocated(proposed_qterm)) deallocate(proposed_qterm)
  if (allocated(cnstph_frc)    ) deallocate(cnstph_frc)

  return

end subroutine cleanup_cnstph_dat

end module constantph_dat_mod
