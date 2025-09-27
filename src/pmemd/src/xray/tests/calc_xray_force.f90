! Primarily used to profile/debug xray modules
program calc_xray_force
  
  use xray_pure_utils, only: real_kind
  use xray_unit_cell_module, only: unit_cell_t
  use xray_debug_dump_module, only: load
  use xray_interface2_module
  
  implicit none
  
  character(len=*), parameter :: xray_params_filename = "dump.txt"
  character(len=*), parameter :: coordinates_filename = "coords.txt"

  integer, allocatable :: hkl(:, :)
  complex(real_kind), allocatable :: Fobs(:)
  real(real_kind), allocatable :: sigma_Fobs(:)
  logical, allocatable :: work_flag(:)
  type(unit_cell_t) :: unit_cell
  real(real_kind), allocatable :: scatter_coefficients(:, :, :) ! Fourier coefficients (2,n_scatter_coeffs,n_scatter_types)
  real(real_kind), allocatable :: atom_b_factor(:)
  real(real_kind), allocatable :: atom_occupancy(:)
  integer, allocatable :: atom_scatter_type(:)
  logical, allocatable :: atom_is_not_bulk(:)
  integer, allocatable :: atom_atomic_number(:)

  real(real_kind), allocatable :: xyz(:, :)
  real(real_kind), allocatable :: force(:, :)
  real(real_kind) :: energy
  real(real_kind) :: weight = 1.0
  real(real_kind) :: r_work, r_free
  real(real_kind) :: k_sol, b_sol

  integer, parameter :: n_steps = 5
  integer :: step
  real(real_kind) :: gamma = 1e-6
  real :: t1, t2 ! time measurement
  
  call load(xray_params_filename, &
      & hkl, Fobs, sigma_Fobs, work_flag, unit_cell, scatter_coefficients, &
      & atom_b_factor, atom_occupancy, atom_scatter_type, atom_is_not_bulk, atom_atomic_number &
      &)
  
  call init(&
      & "ml", "afonine-2013", &
      & hkl, Fobs, sigma_Fobs, work_flag, unit_cell, scatter_coefficients, &
      & atom_b_factor, atom_occupancy, atom_scatter_type, atom_is_not_bulk, atom_atomic_number, &
      & 1, 1, 1, &
      & k_sol, b_sol, &
      & 1.11_real_kind, 0.9_real_kind &
      &)
  
  allocate(xyz(3, size(atom_b_factor)))
  allocate(force(3, size(atom_b_factor)))
  
  call load_coords()
  
  do step = 0, n_steps - 1

    force = 0
    call cpu_time(t1)
    call calc_force(xyz, step, weight, force, energy)
    call cpu_time(t2)
    call get_r_factors(r_work, r_free)
    xyz = xyz + force * gamma

    ! print *, "Step / Energy / R-work / R-free / Max force / Force calc time = ", step, energy, r_work, r_free, maxval(norm2(force, 1)), t2 - t1
  end do
  
  call finalize()
  
contains
  
  subroutine load_coords()
    integer, parameter :: unit_id = 19
    integer :: i
    open(unit=unit_id, action="READ", file=coordinates_filename)
    do i = 1, size(xyz, 2)
      read (unit_id, *) xyz(:, i)
    end do
    close(unit=unit_id)
  end subroutine load_coords

end program calc_xray_force