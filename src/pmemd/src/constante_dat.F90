! Another module needed to avoid cyclic dependencies

module constante_dat_mod

  use constantph_dat_mod, only : proposed_qterm
  ! For consistency, we need to update the proposed changes on a single variable.
  ! For this reason we are only going to using proposed_qterm from constantph_dat.F90

  implicit none

  public

  logical, save :: on_cestep = .false.
  double precision, allocatable, save :: cnste_frc(:,:)

contains

subroutine allocate_cnste_dat(natom, num_reals, ierror)

  implicit none

  integer, intent(in)     :: natom
  integer, intent(in out) :: num_reals
  integer, intent(out)    :: ierror

  ierror = 0
  if (allocated(proposed_qterm)) then
    allocate(cnste_frc(3, natom), stat=ierror)
  else
    allocate(proposed_qterm(natom), cnste_frc(3, natom), &
             stat=ierror)
  end if

  num_reals = num_reals + size(proposed_qterm) + size(cnste_frc)

  return

end subroutine allocate_cnste_dat

subroutine cleanup_cnste_dat

  implicit none

  if (allocated(proposed_qterm)) deallocate(proposed_qterm)
  if (allocated(cnste_frc)    ) deallocate(cnste_frc)

  return

end subroutine cleanup_cnste_dat

end module constante_dat_mod
