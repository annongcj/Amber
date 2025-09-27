#include "copyright.i"

!*******************************************************************************
!
! Module: remd_exchg_mod
!
! Description:
!
! This module holds all of the subroutines that conduct exchanges.
!
! This module needed to be introduced external to the REMD module because many
! of the more advanced exchange subroutines introduced dependencies for remd_mod
! that resulted in cyclic dependencies that prevented compiling. It's probably
! good to have them separate, anyway.
!
! IMPORTANT: If you add new exchange types (in the form of subroutines) here,
! you MUST make sure that, if you must call gb/pme_force, you do so with the
! temporary force array (frc_temp), and in the event of a successful exchange,
! you *must* replace frc with frc_temp so our forces are right for our first
! step after this exchange. Note if the exchange attempt fails, DO NOT perform
! this update. Also, any stand-alone exchange type should write data to the
! remlog if and *only* if it is not called from multid_exchange. Return before
! any print-outs are done if remd_method .eq. -1
!
!*******************************************************************************

module remd_exchg_mod

  use parallel_dat_mod
  use random_mod, only : random_state, amrset_gen, amrand_gen
  use remd_mod

  implicit none

#ifdef MPI
! Constants

  double precision, parameter :: ONEKB = 503.01d0 ! 1/Kb in internal units
  double precision, parameter :: ONE = 1.d0 ! For Non-Boltzmann R-REMD
  double precision, parameter :: TINY = 1.d-6

! To get pH-REMD and E-REMD to agree with sander

  logical, save :: jumpright_ph = .true.
  logical, save :: jumpright_e = .true.

! Random number generator

  type(random_state), save, private :: remd_rand_gen

contains

!*******************************************************************************
!
! Subroutine: multid_exchange
!
! Description:
!
! This subroutine is responsible for coordinating the setting up of the REMD
! comms to facilitate multi-dimensional replica exchange attempts. You'll need
! to edit the select case() sections to add new exchange types
!
!*******************************************************************************

subroutine multid_exchange(atm_cnt, crd, vel, frc, remd_ptot, &
                           actual_temperature, my_atm_lst, gbl_img_atm_map, &
                           gbl_atm_img_map, mdloop)

  use file_io_dat_mod, only : remd_dimension_name, remd_file, max_fn_len
  use pmemd_lib_mod, only   : strip
  use mdin_ctrl_dat_mod, only : iphmd

  implicit none

! Passed variables

  integer, intent(in) :: atm_cnt
  integer, intent(in) :: mdloop
  integer             :: my_atm_lst(*)
  integer             :: gbl_img_atm_map(*)
  integer             :: gbl_atm_img_map(*)

  double precision, intent(in) :: remd_ptot
  double precision, intent(in) :: actual_temperature
  double precision             :: crd(3,atm_cnt)
  double precision             :: frc(3,atm_cnt)
  double precision             :: vel(3,atm_cnt)

! Local variables

  double precision :: l_fe
  double precision :: r_fe

  integer   :: my_dim
  integer   :: i, j, k
  integer   :: rep
  integer   :: gid ! group id

  logical   :: group_master

  character :: torf_char

  character(len=max_fn_len) :: filename
  character(len=max_fn_len) :: buf

! Variable descriptions
!
!  my_dim       : which dimension we choose to exchange between. Alternate
!                 b/w different dimensions
!  i            : counter/iterator
!  buf          : holder string for the remlog filename extension
!  group_master : master of the group_master_comm
!  filename     : the full file name of the remlog we need to write to
!  torf_char    : character for if we succeeded (T) or not (F) in H-REMD output

  ! Determine my_dim, starting with the lowest dimension. We want my_dim to
  ! cycle between 1 - remd_dimension, starting at 1

  my_dim = mod(mdloop - 1, remd_dimension) + 1

  ! Nobody is group master yet, this will be set later

  group_master = .false.

  ! Set up REMD communicators, and free them if they're still formed

  if (master) then

    if (remd_comm .ne. mpi_comm_null) &
      call mpi_comm_free(remd_comm, err_code_mpi)
    remd_comm = mpi_comm_null

    call mpi_comm_split(pmemd_master_comm, group_num(my_dim), &
                        master_rank, remd_comm, err_code_mpi)
    call mpi_comm_size(remd_comm, remd_size, err_code_mpi)
    call mpi_comm_rank(remd_comm, remd_rank, err_code_mpi)

    remd_master = remd_rank .eq. 0

    if (group_master_comm .ne. mpi_comm_null) &
      call mpi_comm_free(group_master_comm, err_code_mpi)
    group_master_comm = mpi_comm_null

    if (remd_master) then
      call mpi_comm_split(pmemd_master_comm, 0, master_rank, &
                          group_master_comm, err_code_mpi)
      call mpi_comm_rank(group_master_comm, group_master_rank, err_code_mpi)
      call mpi_comm_size(group_master_comm, num_rem_grp, err_code_mpi)
      group_master = group_master_rank .eq. 0
    else
      ! If we are not remd_master then call split with MPI_UNDEFINEDs.
      call mpi_comm_split(pmemd_master_comm, MPI_UNDEFINED, MPI_UNDEFINED, &
                          group_master_comm, err_code_mpi)
    end if

  end if ! master

  select case (remd_types(my_dim))

    case(1)
      if(rremd_type .gt. 0 .and. hybridgb .le. 0) then
         call reservoir_exchange(atm_cnt, crd, vel, remd_ptot, my_dim, remd_size, &
                                actual_temperature, .true., mdloop)
      else if(hybridgb .gt. 0) then
         call hybridsolvent_exchange(atm_cnt, crd, vel, frc, numwatkeep, &
                                remd_ptot, my_dim, remd_size, &
                                actual_temperature, .true., mdloop)
      else
         call temperature_exchange(atm_cnt, vel, remd_ptot, my_dim, remd_size, &
                                actual_temperature, .true., mdloop)
      end if
    case(3)
      call hamiltonian_exchange(atm_cnt, my_atm_lst, crd, vel, frc, remd_ptot, &
                                gbl_img_atm_map, gbl_atm_img_map, my_dim, &
                                remd_size, .true., mdloop)
    case(4)
      call ph_remd_exchange(my_dim, remd_size, mdloop)
    case(5)
      call e_remd_exchange(my_dim, remd_size, mdloop)

  end select

  ! Now collect all of the multid_print_data's for all of the various
  ! remd_masters to the group_master

  if (remd_master) then
    call mpi_allgather(multid_print_data, SIZE_REMLOG_DATA*numgroups, &
                       mpi_double_precision, multid_print_data_buf, &
                       SIZE_REMLOG_DATA*numgroups, mpi_double_precision, &
                       group_master_comm, err_code_mpi)

   ! Update the free energy data structures if this dimension is H-REMD

    if (remd_types(my_dim) .eq. 3) then

      do i = 1, num_rem_grp

        gid = int(multid_print_data_buf(1,i)%group_num)

        do j = 1, int(multid_print_data_buf(1,i)%num_rep)

          rep = int(multid_print_data_buf(j,i)%repnum)
          num_left_exchg(my_dim,gid,rep) = num_left_exchg(my_dim,gid,rep) &
                 + int(multid_print_data_buf(j,i)%left_exchg)
          num_right_exchg(my_dim,gid,rep) = num_right_exchg(my_dim,gid,rep) &
                 + int(multid_print_data_buf(j,i)%right_exchg)
          total_left_fe(my_dim,gid,rep) = total_left_fe(my_dim,gid,rep) &
                 + multid_print_data_buf(j,i)%left_fe
          total_right_fe(my_dim,gid,rep) = total_right_fe(my_dim,gid,rep) &
                 + multid_print_data_buf(j,i)%right_fe

        end do ! j = 1, num_rep

      end do ! i = 1, num_rem_grp

    end if ! remd_types(my_dim) .eq. 3

  end if

  ! Now that we have all of the data on group_master, have *that* process
  ! open up the appropriate remlog and write out the data to that file. We
  ! call the intrinsic "open" so we can append to this old file (it should
  ! have already been opened via amopen in the setup routine). This is the
  ! approach taken for the hard-flushing of mdout/mdinfo as well (so this has
  ! the side-effect of flushing this rem.log every time we reach this point

  if (group_master) then

    write(buf, '(i5)') my_dim
    call strip(buf)
    filename = trim(remlog_name) // '.' // trim(buf)

    open(unit=remlog, file=filename, status='OLD', position='APPEND')

    ! Modify your particular REMD's exchange method's printout here:

    select case(remd_types(my_dim))

      case(1) ! TEMPERATURE
        do i = 1, num_rem_grp

          gid = int(multid_print_data_buf(1,i)%group_num)

          write(remlog, '(2(a,i8))') '# exchange ', mdloop, ' REMD group ', gid

          do j = 1, int(multid_print_data_buf(1,i)%num_rep)
            write(remlog, '(i2,6f10.2,i8)') &
               j, &
               multid_print_data_buf(j,i)%scaling,       &
               multid_print_data_buf(j,i)%real_temp,     &
               multid_print_data_buf(j,i)%pot_ene_tot,   &
               multid_print_data_buf(j,i)%temp0,         &
               multid_print_data_buf(j,i)%new_temp0,     &
               multid_print_data_buf(j,i)%success_ratio, &
               int(multid_print_data_buf(j,i)%struct_num)
          end do ! j = 1, %num_rep

        end do ! i = 1, num_rem_grp

      case(3) ! HAMILTONIAN
        ! We need to update our running tally of number of left/right exchanges
        ! as well as our total left and right free energies
        do i = 1, num_rem_grp

          gid = int(multid_print_data_buf(1,i)%group_num)

          write(remlog, '(2(a,i8))') '# exchange ', mdloop, ' REMD group ', gid

          do j = 1, int(multid_print_data_buf(1,i)%num_rep)

            rep = int(multid_print_data_buf(j,i)%repnum)

            if (num_right_exchg(my_dim, i, rep) .eq. 0) then
              r_fe = 0.d0
            else
              r_fe = multid_print_data_buf(j,i)%temp0 / ONEKB * &
                     log(total_right_fe(my_dim, gid, rep) / &
                         num_right_exchg(my_dim, gid, rep))
            end if

            if (num_left_exchg(my_dim, i, rep) .eq. 0) then
              l_fe = 0.d0
            else
              l_fe = multid_print_data_buf(j,i)%temp0 / ONEKB * &
                     log(total_left_fe(my_dim, gid, rep) / &
                         num_left_exchg(my_dim, gid, rep))
            end if

            if (multid_print_data_buf(j,i)%success .eq. 1.d0) then
              torf_char = 'T'
            else
              torf_char = 'F'
            end if

            write(remlog, '(2i6,5f10.2,4x,a,2x,f10.2)')      &
               int(multid_print_data_buf(j,i)%repnum),       &
               int(multid_print_data_buf(j,i)%neighbor_rep), &
               multid_print_data_buf(j,i)%temp0,             &
               multid_print_data_buf(j,i)%pot_ene_tot,       &
               multid_print_data_buf(j,i)%nei_pot_ene,       &
               l_fe, &
               r_fe, &
               torf_char, &
               multid_print_data_buf(j,i)%success_ratio

          end do ! j = 1, %num_rep

        end do ! i = 1, num_rem_grp

      case(4) ! PH
        do i = 1, num_rem_grp

          gid = int(multid_print_data_buf(1,i)%group_num)

          write(remlog, '(2(a,i8))') '# exchange ', mdloop, ' REMD group ', gid

          do j = 1, int(multid_print_data_buf(1,i)%num_rep)
            !PHMD
            if(iphmd .ne. 0) then
               write(remlog, '(i6,x,3f7.3,x,f8.4)') &
                  j, &
                  multid_print_data_buf(j,i)%nprot,   &
                  multid_print_data_buf(j,i)%my_ph,         &
                  multid_print_data_buf(j,i)%nei_ph,        &
                  multid_print_data_buf(j,i)%success_ratio
            else
               write(remlog, '(i6,x,i7,x,2f7.3,x,f8.4)') &
                  j, &
                  nint(multid_print_data_buf(j,i)%nprot),   &
                  multid_print_data_buf(j,i)%my_ph,         &
                  multid_print_data_buf(j,i)%nei_ph,        &
                  multid_print_data_buf(j,i)%success_ratio
            end if
          end do ! j = 1, %num_rep

        end do ! i = 1, num_rem_grp

      case(5) ! REDOX
        do i = 1, num_rem_grp

          gid = int(multid_print_data_buf(1,i)%group_num)

          write(remlog, '(2(a,i8))') '# exchange ', mdloop, ' REMD group ', gid

          do j = 1, int(multid_print_data_buf(1,i)%num_rep)
            write(remlog, '(i6,x,i7,x,2f7.3,x,f8.4)') &
               j, &
               nint(multid_print_data_buf(j,i)%nelec),   &
               multid_print_data_buf(j,i)%my_e,          &
               multid_print_data_buf(j,i)%nei_e,         &
               multid_print_data_buf(j,i)%success_ratio
          end do ! j = 1, %num_rep

        end do ! i = 1, num_rem_grp

    end select

    close(remlog)

  end if ! group_master

  return

end subroutine multid_exchange

!*******************************************************************************
!
! Subroutine: pv_correction 
!
! Description: Calculate the pressure/volume correction to exchange delta. 
!
!*******************************************************************************
function pv_correction(ourtemp, nbrtemp, neighbor, comm_rep_master)
  ! USE STATEMENTS
  use gbl_constants_mod, only : KB
  ! ARGUMENTS
  implicit none
  double precision pv_correction
  double precision, intent(in)  :: ourtemp, nbrtemp
  integer, intent(in) :: neighbor
  integer, intent(in) :: comm_rep_master
  ! INCLUDES
  include 'mpif.h'
  ! LOCAL VARIABLES
  double precision, dimension(2) :: our, nbr
  integer ierr, istat(mpi_status_size)

  ! Pressure in position 1, volume in position 2.
  our(1) = remd_pressure
  our(2) = remd_volume
  call mpi_sendrecv(our, 2, mpi_double_precision, neighbor, 52, &
                    nbr, 2, mpi_double_precision, neighbor, 52, &
                    comm_rep_master, istat, ierr)
# ifdef VERBOSE_REMD
  write(mdout,'(3(a,f16.8))') '| REMD: OurPressure= ', our(1), ' atm, OurVolume= ', our(2), &
          ' Ang^3, OurTemp= ', ourtemp
  write(mdout,'(3(a,f16.8))') '| REMD: NbrPressure= ', nbr(1), ' atm, NbrVolume= ', nbr(2), &
          ' Ang^3, NbrTemp= ', nbrtemp
# endif
  pv_correction = (((1.d0/(KB*nbrtemp))*nbr(1)) - ((1.d0/(KB*ourtemp))*our(1))) &
                  * (nbr(2) - our(2))
  ! Need to convert from units of atm*Ang^3 to kcal/mol
  ! 1 atm*Ang^3 * 101325 Pa/atm * E-30 m^3/Ang^3 * J/Pa*m^3 * kcal/4184 J * 6.02E23/mol
  pv_correction = pv_correction * 1.4584E-5
# ifdef VERBOSE_REMD
  write(mdout,'(a,E16.8)') '| REMD: PvCorrection= ', pv_correction
# endif
end function pv_correction

!*******************************************************************************
!
! Subroutine: temperature_exchange
!
! Description: Performs the temperature replica exchange attempt. We don't need
!              any extra energies here.
!
!*******************************************************************************

subroutine temperature_exchange(atm_cnt, vel, my_pot_ene_tot, t_dim, &
                                num_replicas, actual_temperature,    &
                                print_exch_data, mdloop)

  use mdin_ctrl_dat_mod, only : temp0
  use pmemd_lib_mod,     only : mexit

  implicit none

! Passed variables

  integer, intent(in)             :: atm_cnt
  integer, intent(in)             :: num_replicas
  integer, intent(in)             :: mdloop

  double precision, intent(in)    :: my_pot_ene_tot
  double precision, intent(inout) :: vel(3, atm_cnt)
  double precision, intent(in)    :: actual_temperature

  integer, intent(in)             :: t_dim ! which dimension T-exchanges are
  logical, intent(in)             :: print_exch_data

! Local variables

  type :: exchange_data ! facilitates gather for remlog
    sequence
    double precision :: scaling
    double precision :: real_temp
    double precision :: pot_ene_tot
    double precision :: temp0
    double precision :: new_temp0
    double precision :: struct_num
  end type exchange_data

  integer, parameter :: SIZE_EXCHANGE_DATA = 6 ! for mpi_gather

  type (exchange_data) :: exch_buffer(num_replicas) ! gather buffer

  type (exchange_data) :: my_exch_data ! stores data for this replica
  type (exchange_data) :: exch_data_tbl(num_replicas) ! exch_data for all reps

  double precision  :: delta
  double precision  :: metrop
  double precision  :: pvterm
  double precision  :: random_value
  double precision  :: success_ratio

  integer           :: neighbor_rank
  integer           :: i
  integer           :: success_array(num_replicas)
  integer           :: success_buf(num_replicas)

  ! For exchanging replica positions in ladders upon successful exchanges
  integer           :: group_num_buffer(remd_dimension)
  integer           :: replica_indexes_buffer(remd_dimension)

  logical           :: success
  logical           :: i_do_exchg

  integer, dimension(mpi_status_size) :: stat_array ! needed for mpi_recv

! Explanation of local variables:
!
! delta: The calculated DELTA value used as part of the detailed balance calc.
!
! metrop: The actual value compared to a random_value to decide success
!
! success: Did this exchange attempt succeed?
!
! success_array: buffer for keeping track of which temperatures exchanged
!                successfully. Used in a mpi_reduce at the end
!
! success_buf: buffer for mpi_reduce call on success_array
!
! neighbor_rank: the remd_rank of the replica with whom we need to be doing ALL
!                of our communication. Trying to do everything with the
!                'partners' array is a bit too confusing.

! Set the variables that we know at the beginning

  my_exch_data%real_temp   = actual_temperature
  my_exch_data%temp0       = temp0
  my_exch_data%new_temp0   = temp0 ! changed later if exch. succeeds
  my_exch_data%struct_num  = -1    ! not implemented yet
  my_exch_data%pot_ene_tot = my_pot_ene_tot
  my_exch_data%scaling     = -1.d0
  success_array(:)         = 0
  pvterm = 0.d0

  if (master) then

    call set_partners(t_dim, num_replicas, remd_random_partner)

    ! Call amrand_gen here to generate a random number. Sander does this
    ! to pick out structures from a reservoir which pmemd doesn't implement
    ! (yet). To keep the random #s synched, though, and to ensure proper
    ! regression testing, this must be done to keep the random #s synched
    ! b/w sander and pmemd

    call amrand_gen(remd_rand_gen, random_value)

#ifdef VERBOSE_REMD
    write(mdout, '(a,26("="),a,26("="))') '| ', 'REMD EXCHANGE CALCULATION'
    write(mdout, '(a,i10,a,i1)') '| Exch= ', mdloop, ' RREMD= ', 0
#endif

    ! We alternate exchanges: temp up followed by temp down
    if (even_exchanges(t_dim)) then
      i_do_exchg = control_exchange(t_dim)
    else
      i_do_exchg = .not. control_exchange(t_dim)
    end if

    even_exchanges(t_dim) = .not. even_exchanges(t_dim)

    ! Find out what rank our neighbor is (-1 for indexing from 0), and
    ! If we are controlling the exchange, then our partner will be our
    ! 2nd partner (above). Otherwise it will be our 1st partner (below)

    if (i_do_exchg) then
      neighbor_rank = partners(2) - 1
    else
      neighbor_rank = partners(1) - 1
    end if

    ! Collect potential energies/temperatures

    call mpi_allgather(my_exch_data, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                exch_data_tbl, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                remd_comm, err_code_mpi)

#ifdef VERBOSE_REMD
    write(mdout, '(a,f6.2,2(a7,i2),a7,f10.2)') '| Replica          Temp= ', &
      my_exch_data%temp0, ' Indx= ', replica_indexes(t_dim), ' Rep#= ', &
      remd_rank+1, ' EPot= ', my_pot_ene_tot
    write(mdout,'(a7,i6,a8,i6)') '| RepIdx=', remd_repidx, ' CrdIdx=', remd_crdidx
#endif

    ! Calculate pressure/volume correction if necessary
    if (use_pv) pvterm = pv_correction(my_exch_data%temp0, &
                                       exch_data_tbl(neighbor_rank+1)%temp0, &
                                       neighbor_rank, remd_comm)

    ! Calculate the exchange probability. Note that if delta is negative, then
    ! metrop will be > 1, and it will always succeed (so we don't need an extra
    ! conditional)

    if (i_do_exchg) then
#ifdef VERBOSE_REMD
      write(mdout, '(a,f6.2,2(a7,i2),a7,f10.2)') '| Partner          Temp= ', &
        exch_data_tbl(neighbor_rank+1)%temp0, ' Indx= ', neighbor_rank+1, &
        ' Rep#= ', neighbor_rank + 1, ' EPot= ', &
        exch_data_tbl(neighbor_rank+1)%pot_ene_tot
#endif

      delta = (my_exch_data%pot_ene_tot - &
               exch_data_tbl(neighbor_rank+1)%pot_ene_tot) * &
              (my_exch_data%temp0 - &
               exch_data_tbl(neighbor_rank+1)%temp0) * ONEKB /  &
              (my_exch_data%temp0 * exch_data_tbl(neighbor_rank+1)%temp0)

      metrop = exp(-delta)

      call amrand_gen(remd_rand_gen, random_value)

      success = random_value .lt. metrop

      ! Let my neighbor know about our success

      call mpi_send(success, 1, mpi_logical, neighbor_rank, &
                    remd_tag, remd_comm, err_code_mpi)

      if (success) then
        success_array(replica_indexes(t_dim)) = 1
        my_exch_data%new_temp0 = exch_data_tbl(neighbor_rank+1)%temp0
        my_exch_data%scaling = sqrt(my_exch_data%new_temp0 / my_exch_data%temp0)
      end if

#ifdef VERBOSE_REMD
      write(mdout,'(a8,E16.6,a8,E16.6,a12,f10.2)') &
               "| Metrop= ",metrop," delta= ",delta," o_scaling= ", &
               1 / my_exch_data%scaling
#endif

    else

#ifdef VERBOSE_REMD
      ! Even if we don't control exchange, still print REMD info to mdout

      write(mdout, '(a,f6.2,2(a7,i2),a7,f10.2)') '| Partner          Temp= ', &
        exch_data_tbl(neighbor_rank+1)%temp0, ' Indx= ', neighbor_rank+1, &
        ' Rep#= ', neighbor_rank + 1, ' EPot= ', &
        exch_data_tbl(neighbor_rank+1)%pot_ene_tot
      write(mdout, '(a)') '| Not controlling exchange.'
#endif

      call amrand_gen(remd_rand_gen, random_value) ! to stay synchronized

      ! Get the message from the exchanging replica about our success
      call mpi_recv(success, 1, mpi_logical, neighbor_rank, &
                    remd_tag, remd_comm, stat_array, err_code_mpi)

      ! We scale velocities only if we succeed
      if (success) then
        my_exch_data%new_temp0 = exch_data_tbl(neighbor_rank+1)%temp0
        my_exch_data%scaling = sqrt(my_exch_data%new_temp0 / my_exch_data%temp0)
      end if

    end if

#ifdef VERBOSE_REMD
    write(mdout,'(a8,E16.6,a12,f10.2,a10,L1)') &
      '| Rand=   ', random_value, ' MyScaling= ', my_exch_data%scaling, &
      ' Success= ', success

    write(mdout,'(a,24("="),a,24("="))') '| ', "END REMD EXCHANGE CALCULATION"
#endif

  end if ! master

  ! We'll broadcast all of the my_exch_data here once. From that, we can tell
  ! if the exchange succeeded based on whether or not the temperatures have
  ! changed. Then make sure we update temp0 so the simulation reflects that.

  call mpi_bcast(my_exch_data, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                 0, pmemd_comm, err_code_mpi)

  temp0 = my_exch_data%new_temp0

  if (.not. master) &
    success = abs(my_exch_data%temp0 - my_exch_data%new_temp0) .gt. TINY

#ifdef VERBOSE_REMD
  if (master) &
    write(mdout, '(2(a,f6.2))') &
      '| REMD: checking to see if bath T has changed: ', &
      my_exch_data%temp0, '->', my_exch_data%new_temp0
#endif

  if (success) call rescale_velocities(atm_cnt, vel, my_exch_data%scaling)

  ! We want to keep track of our exchange successes, so gather them
  ! from all of the masters here, and increment them in exchange_successes
  ! We do this whether or not we print so we can keep track of exchange
  ! success ratios. For the index and exchange data, only gather those to
  ! the remd_master if we actually plan to print. This should be much
  ! more efficient for high exchange attempts.

  ! Then recollect the temperatures as some may have changed, and reset our
  ! partners.

  if (master) then

    if (print_exch_data) then
      call mpi_gather(my_exch_data, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                      exch_buffer, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                      0, remd_comm, err_code_mpi)
    end if

    call mpi_reduce(success_array, success_buf, num_replicas, mpi_integer, &
                    mpi_sum, 0, remd_comm, err_code_mpi)

    if (remd_master) &
      exchange_successes(t_dim,:) = exchange_successes(t_dim,:) + success_buf(:)

    ! Increment exchange counter, and swap our replica ranks and numbers with
    ! our neighbor if our attempt succeeded, since we effectively swapped places
    ! in this array. That's because we swap temperatures and not structures due
    ! to performance considerations.

    if (success) then
      call mpi_sendrecv(group_num, remd_dimension, mpi_integer, &
                        neighbor_rank, remd_tag, &
                        group_num_buffer, remd_dimension, mpi_integer, &
                        neighbor_rank, remd_tag, remd_comm, &
                        stat_array, err_code_mpi)
      call mpi_sendrecv(replica_indexes, remd_dimension, mpi_integer, &
                        neighbor_rank, remd_tag, &
                        replica_indexes_buffer, remd_dimension, mpi_integer, &
                        neighbor_rank, remd_tag, remd_comm, &
                        stat_array, err_code_mpi)
      call mpi_sendrecv_replace(remd_repidx, 1, mpi_integer, &
                                neighbor_rank, remd_tag, &
                                neighbor_rank, remd_tag, &
                                remd_comm, stat_array, err_code_mpi)
      group_num = group_num_buffer
      replica_indexes = replica_indexes_buffer
    end if
  end if

  ! If we're doing NMROPT, then we need to make sure we re-read in the TEMP0,
  ! since this subroutine overwrites TEMP0 to the original since you can modify
  ! temperature with NMROPT settings.

  remd_modwt = .true.

  ! Write the data to the remlog

  ! remd_method == -1 means we're doing multi-D REMD.
  if (remd_method .ne. -1) then
    if (print_exch_data .and. remd_master) then
      write(remlog, '(a,i8)') '# exchange ', mdloop
      do i = 1, num_replicas
        success_ratio = dble(exchange_successes(t_dim, index_list(i))) / &
                        dble(mdloop) * dble(remd_dimension) * 2
        write(remlog, '(i2,f12.4,f10.2,f15.2,3f10.2,i8)') &
              i, &
              exch_buffer(i)%scaling, &
              exch_buffer(i)%real_temp, &
              exch_buffer(i)%pot_ene_tot, &
              exch_buffer(i)%temp0, &
              exch_buffer(i)%new_temp0, &
              success_ratio, &
              int(exch_buffer(i)%struct_num)
      end do
    end if
  else if (remd_master) then
    do i = 1, num_replicas
      success_ratio = dble(exchange_successes(t_dim, index_list(i))) / &
                      dble(mdloop) * dble(remd_dimension) * 2
      multid_print_data(i)%scaling       = exch_buffer(i)%scaling
      multid_print_data(i)%real_temp     = exch_buffer(i)%real_temp
      multid_print_data(i)%pot_ene_tot   = exch_buffer(i)%pot_ene_tot
      multid_print_data(i)%temp0         = exch_buffer(i)%temp0
      multid_print_data(i)%new_temp0     = exch_buffer(i)%new_temp0
      multid_print_data(i)%group_num     = group_num(t_dim)
      multid_print_data(i)%success_ratio = success_ratio
      multid_print_data(i)%num_rep       = num_replicas
    end do
  end if ! remd_method .ne. -1

  return

end subroutine temperature_exchange

!*******************************************************************************
!
! Subroutine: remd_bcast_cell
!
! Description: Broadcast new unit cell information in box1 and update unit cell.
!
!*******************************************************************************
subroutine remd_bcast_cell(box0, box1, commpmemd)
  ! USE STATEMENTS
  use pbc_mod, only : pbc_box, pressure_scale_pbc_data, ucell, recip, &
                      uc_volume, cut_factor
  use mdin_ctrl_dat_mod, only : vdw_cutoff
  use mdin_ewald_dat_mod, only : skinnb, verbose
  use cit_mod, only : set_cit_tbl_dims
  use pmemd_lib_mod,     only : mexit

  ! ARGUMENTS
  implicit none
  double precision, intent(in), dimension(:) :: box0       ! Original box
  double precision, intent(in), dimension(:) :: box1       ! New box
  integer, intent(in)                        :: commpmemd ! COMM for broadcast
  ! INCLUDES
  include 'mpif.h'
  ! LOCAL VARS
  integer ierr
  double precision rmu(3)

  call mpi_bcast(box1, 3, mpi_double_precision, 0, commpmemd, ierr)
  ! Calculate scale from old box to new
  rmu(1:3) = box1(1:3) / box0(1:3)
  if (abs(rmu(1)/rmu(2)-1.0) .gt. 1e-5 .or. abs(rmu(1)/rmu(3)-1.0) .gt. 1e-5) then
    write(mdout, '(a)') '| ERROR: The replica boxes do not have the same shape. '
    write(mdout, '(a)') '| ERROR: They can only differ by a common scaling factor in all dimensions. '
    call mexit(mdout,1)
  endif 
  !write(6, '(a,6f8.3,3(x,E10.4))') '| DBG: Box {old, new, scale} = ', box0(1:3), box1(1:3), rmu(1:3)
  call pressure_scale_pbc_data(rmu, vdw_cutoff + skinnb, verbose)
#ifdef CUDA
  !call gpu_pressure_scale(ucell, recip, uc_volume)
  call gpu_ucell_set(ucell, recip, uc_volume)
#else
  call set_cit_tbl_dims(pbc_box, vdw_cutoff + skinnb, cut_factor)
#endif
end subroutine remd_bcast_cell

!*******************************************************************************
!
! Subroutine: remd_update_mol_com 
!
! Description: Ensure molecule CoM are consistent with given coords. 
!
!*******************************************************************************
subroutine remd_update_mol_com(crd)
  ! USE STATEMENTS
  use mol_list_mod,      only : gbl_mol_cnt
  use dynamics_dat_mod,  only : my_mol_cnt, gbl_mol_mass_inv, gbl_my_mol_lst, gbl_mol_com
  use prmtop_dat_mod,    only : atm_mass
  use dynamics_mod,      only : get_mol_com, get_mol_com_midpoint
  use mdin_ctrl_dat_mod, only : usemidpoint
  use pbc_mod,           only : old_pbc
  ! ARGUMENTS
  implicit none
  double precision :: crd(3,*)
  ! INCLUDES
  ! LOCAL VARS

  if (usemidpoint) then
    call get_mol_com_midpoint(gbl_mol_cnt, crd, atm_mass, gbl_mol_mass_inv, gbl_mol_com, old_pbc)
  else
    call get_mol_com(my_mol_cnt, crd, atm_mass, gbl_mol_mass_inv, gbl_my_mol_lst, gbl_mol_com)
  endif
end subroutine remd_update_mol_com

!*******************************************************************************
!
! Subroutine: hamiltonian_exchange
!
! Description: Performs the Hamiltonian replica exchange attempt. Additional
!              energies need to be computed here.
!
!*******************************************************************************

subroutine hamiltonian_exchange(atm_cnt, my_atm_lst, crd, vel, frc, &
                                my_pot_ene_tot, gbl_img_atm_map, &
                                gbl_atm_img_map, h_dim, &
                                num_replicas, print_exch_data, mdloop)

  use gb_force_mod
  use mdin_ctrl_dat_mod
  use nmr_calls_mod, only : nmrdcp, skip_print
  use pme_force_mod
  ! Below needed when exchanging unit cell information
  use axis_optimize_mod
  use pbc_mod,          only : pbc_box
  use reservoir_mod
! NFE
  use nfe_lib_mod,   only : nfe_real_mdstep
  
  use ti_mod

  implicit none

! Passed variables

  integer          :: atm_cnt
  integer          :: num_replicas
  integer          :: mdloop
  integer          :: my_atm_lst(*)
  integer          :: gbl_img_atm_map(*)
  integer          :: gbl_atm_img_map(*)
  integer          :: my_idx
  logical          :: print_exch_data
  double precision :: crd(3,atm_cnt)
  double precision :: frc(3,atm_cnt)
  double precision :: vel(3,atm_cnt)
  integer          :: h_dim
  double precision :: my_pot_ene_tot

! Local variables

  type :: ene_temp ! Useful type to limit communication
    sequence
    double precision :: neighbor_rep ! neighbor replica # (for remlog printing)
    double precision :: energy_1     ! energy with MY coordinates
    double precision :: energy_2     ! energy with THEIR coordinates
    double precision :: temperature  ! my temperature
    double precision :: repnum       ! replica number
    double precision :: left_fe      ! My left free energy for this step
    double precision :: right_fe     ! My right free energy for this step
    double precision :: right_exchg  ! 1 if we exchanged to the right, 0 left
    double precision :: left_exchg   ! 1 if we exchanged to the left, 0 right
  end type ene_temp
  integer, parameter :: SIZE_ENE_TEMP = 9

  double precision :: my_pot_ene_tot_2   ! MY pot ene with THEIR coordinates
  double precision :: delta
  double precision :: metrop
  double precision :: pvterm
  double precision :: original_box(3)    ! Box XYZ before exchange
  double precision :: new_box(3)         ! Neighbor box XYZ
  double precision :: boxtmp(3)          ! Temp storage for box when exchanging
  double precision :: random_value
  double precision :: success_ratio
  double precision :: virial(3)          ! For pme_force_call
  double precision :: ekcmt(3)           ! For pme_force_call
  double precision :: pme_err_est        ! For pme_force_call
  double precision :: vel_scale_coff     ! Velocity scaling factor
  double precision :: l_fe, r_fe         ! Holders for FE calcs
  double precision :: frc_bak(3,atm_cnt)

  logical,save          :: firstTime=.true.
  logical          :: i_do_exchg
  logical          :: silent
  logical          :: do_reservoir
  logical          :: rank_success(num_replicas) ! log of all replica's success
  logical          :: success
  integer          :: success_array(num_replicas) ! to keep track of successes
  integer          :: success_buf(num_replicas) ! reduce buffer for successes
  character        :: success_char ! to print T or F if we succeeded in remlog
  integer          :: neighbor_rank
  integer          :: istat(mpi_status_size)
  integer          :: ier, i,j
  integer          :: rep ! replica counter
  integer          :: ord1, ord2, ord3 ! For obtaining cell lengths
  integer,save          :: tiRegion
  double precision, parameter :: energy_cut=500.0

  integer, dimension(mpi_status_size) :: stat_array ! needed for mpi_recv
    
  type(gb_pot_ene_rec)  :: my_new_gb_ene  ! gb energy record for THEIR coords
  type(pme_pot_ene_rec) :: my_new_pme_ene ! pme energy record for THEIR coords
  type(ene_temp)   :: my_ene_temp
  type(ene_temp)   :: neighbor_ene_temp
  type(ene_temp)   :: ene_temp_buffer(num_replicas) ! to collect/print REMD data

  if(firstTime .and. (rremd_type .ge. 4 .and. rremd_type .le. 6) ) then
    
    tiRegion=-1
    if ( (rremd_type .eq. 4 .or. rremd_type .eq. 6) .and. replica_indexes(1) .eq. 1) then
      tiRegion=1
    else if ( (rremd_type .eq. 5 .or. rremd_type .eq. 6) .and. replica_indexes(1) .eq. num_replicas) then
      tiRegion=2
    endif 

    allocate(rremd_crd(3,atm_cnt))

  end if
  firstTime=.false.
  
  ! Get coordinates and forces from GPU
#ifdef CUDA
  call gpu_download_crd(crd)
  call gpu_download_frc(frc)
#endif

  pvterm = 0.d0

  ! Set up some stuff for pme_force call:

  virial(:) = 0.d0
  ekcmt(:) = 0.d0
  
  vel_scale_coff=1.0
  ! First determine if we're controlling the exchange or not

  silent=.false.
  if (master) then

    ! Figure out our partners

    call set_partners(h_dim, num_replicas, remd_random_partner)

    ! Initialize some data

    success_array(:) = 0
    my_ene_temp%energy_1 = my_pot_ene_tot
    my_ene_temp%energy_2 = 0;
    my_ene_temp%temperature = temp0
    my_ene_temp%repnum = replica_indexes(h_dim)
    my_ene_temp%left_fe = 0.d0
    my_ene_temp%right_fe = 0.d0
    my_ene_temp%right_exchg = 0.d0
    my_ene_temp%left_exchg = 0.d0

#ifdef VERBOSE_REMD
    write(mdout, '(a,24("="),a,24("="))') '| ', ' H-REMD EXCHANGE CALCULATION '
    write(mdout, '(a,i10,a,i1)') '| Exch= ', mdloop, ' RREMD= ', 0
#endif

    ! Alternate exchanges

    if (even_exchanges(h_dim)) then
      i_do_exchg = control_exchange(h_dim)
    else
      i_do_exchg = .not. control_exchange(h_dim)
    end if

    even_exchanges(h_dim) = .not. even_exchanges(h_dim)

    ! Set partner rank
    if (i_do_exchg) then
      my_ene_temp%neighbor_rep = partners(2)
      neighbor_rank = partners(2) - 1
    else
      my_ene_temp%neighbor_rep = partners(1)
      neighbor_rank = partners(1) - 1
    end if

    silent=(neighbor_rank .lt. 0)
  endif ! master

  ! Make sure non-master processes know about silent
  call mpi_bcast(silent, 1, mpi_logical, 0, pmemd_comm, err_code_mpi)
  if (silent) i_do_exchg = .false.
  call mpi_bcast(i_do_exchg, 1, mpi_logical, 0, pmemd_comm, err_code_mpi)

  ! Back to master-only stuff
  if (master) then
    
    ! Exchange coordinates with your neighbor. crd_temp must be a backup since
    ! the reciprocal space sum needs the global array to be updated with the new
    ! coordinates (and crd is a reference to the global array)
    crd_temp(1,:) = crd(1,:)
    crd_temp(2,:) = crd(2,:)
    crd_temp(3,:) = crd(3,:)
    if (.not. silent) call mpi_sendrecv_replace(crd, atm_cnt * 3, mpi_double_precision, &
               neighbor_rank, remd_tag, &
               neighbor_rank, remd_tag, &
               remd_comm, istat, err_code_mpi)
    
    ! Exchange unit cells if necessary TODO are there other things that need exchanging
    if (exchange_ucell) then
      ! To ensure axes order remains consistent, exchange box parameters in order XYZ
      ord1 = axis_flipback_ords(1)
      ord2 = axis_flipback_ords(2)
      ord3 = axis_flipback_ords(3)
      original_box(1) = pbc_box(ord1)
      original_box(2) = pbc_box(ord2)
      original_box(3) = pbc_box(ord3)
      new_box(ord1) = original_box(1)
      new_box(ord2) = original_box(2)
      new_box(ord3) = original_box(3)
      if (.not. silent) then
        call mpi_sendrecv(original_box, 3, mpi_double_precision, neighbor_rank, 11, &
                          boxtmp,       3, mpi_double_precision, neighbor_rank, 11, &
                          remd_comm, istat, err_code_mpi)
        new_box(ord1) = boxtmp(1)
        new_box(ord2) = boxtmp(2)
        new_box(ord3) = boxtmp(3)
      endif
    end if
           
    !! reservoir H-REMD; only valid at the end states
    do_reservoir=.false. 
    enforce_reset_velocities=.false.
  
    if(silent .and. rremd_type.ge.4 .and. rremd_type.le.6) then
      call amrand_gen(remd_rand_gen, random_value) ! keep in sync
      if (tiRegion .gt. 0 .and. neighbor_rank .lt. 0) then
        j=1
        if (mod(num_replicas,2) .eq.1 .and. tiRegion .eq. 2) j=0
        if(mod(mdloop,reservoir_exchange_step).eq.j) then
            do_reservoir=.true. 
            call ti_mergeReservoir(atm_cnt, crd, vel, boxtmp, random_value, tiRegion)
            if (exchange_ucell) then
              new_box(ord1) = boxtmp(1)
              new_box(ord2) = boxtmp(2)
              new_box(ord3) = boxtmp(3)
            endif
        endif
      endif
    endif
    
   end if ! master

  ! backup forces since forces may be reset during the load_balancing of
  ! pme_force() below
  frc_bak(1,:) = frc(1,:)
  frc_bak(2,:) = frc(2,:)
  frc_bak(3,:) = frc(3,:)

  ! Now broadcast the temporary coordinates to the rest of pmemd_comm. I don't
  ! know if this is necessary or not (i.e. is it done in the load balancing?)

  if (.not. silent) call mpi_bcast(crd, atm_cnt * 3, mpi_double_precision, &
                 0, pmemd_comm, err_code_mpi)
  if (.not. silent) call mpi_bcast(crd_temp, atm_cnt * 3, mpi_double_precision, &
                 0, pmemd_comm, err_code_mpi)

#ifdef CUDA
  ! Upload the temporary coordinates
  call gpu_upload_crd(crd)
#endif

  ! Update the unit cell if necessary
  if ( (.not. silent .or. do_reservoir).and. exchange_ucell) &
    call remd_bcast_cell(pbc_box, new_box, pmemd_comm)
  !write(6,'(a,3f8.3)') '| DBG: New box: ', pbc_box(1:3)

  ! We do not want to print out any DUMPAVE values for this force call

  if (nmropt .ne. 0) skip_print = .true.

  ! Don't print for NFE
  if (infe .ne. 0) nfe_real_mdstep = .false.

  ! Now it's time to get the energy of the other structure:
  if (using_gb_potential) then
    call gb_force(atm_cnt, crd, frc_temp, my_new_gb_ene, irespa, .true.)

    my_ene_temp%energy_2 = my_new_gb_ene%total

  else if (using_pme_potential) then
#ifdef CUDA
    ! CUDA pairlist building is not controlled by new_list, it is handled
    ! internally. Since we are changing coordinates, we must build a new
    ! pairlist
    call gpu_force_new_neighborlist()
#endif
    ! First .true. is new_list, next is need_pot_enes, .false. is need_virials
    ! Last value is nstep - set to -1 here since we don't have nstep in this
    ! routine and it may??? not be relevant here?

    if (.not. silent .or. do_reservoir) then
      call pme_force(atm_cnt, crd, frc_temp, gbl_img_atm_map, &
                     gbl_atm_img_map, my_atm_lst, .true., .true., .false., &
                     my_new_pme_ene, -1, virial, ekcmt, pme_err_est)
      my_ene_temp%energy_2 = my_new_pme_ene%total
      ! Since we made an 'extra' force call just now, decrement the nmropt counter
      if (nmropt .ne. 0) call nmrdcp
    endif
  end if

  neighbor_ene_temp = my_ene_temp
  ! Now trade that potential energy with your partner
  if (master) then

    if (.not. silent) call mpi_sendrecv(my_ene_temp, SIZE_ENE_TEMP, mpi_double_precision, &
                      neighbor_rank, remd_tag, &
                      neighbor_ene_temp, SIZE_ENE_TEMP, mpi_double_precision, &
                      neighbor_rank, remd_tag, &
                      remd_comm, istat, err_code_mpi)
    ! Calculate the PV correction if necessary
    if (.not.silent .and. exchange_ucell .and. use_pv) &
      pvterm = pv_correction(temp0, neighbor_ene_temp%temperature, &
                             neighbor_rank, remd_comm)

  ! Now it's time to calculate the exchange criteria. This detailed balance was
  ! derived under the assumption that *if* we have different T's, then they are
  ! NOT exchanged. Swapping temperatures requires a slightly different equation
  ! and then requires you to update temp0 at the end with an mpi_sendrecv call.

    success = .false.
    if (.not. silent) then
    if (i_do_exchg) then
      if (abs(my_ene_temp%energy_1-my_ene_temp%energy_2)<energy_cut .and. &
          abs(neighbor_ene_temp%energy_1-neighbor_ene_temp%energy_2)<energy_cut)  then ! exclude extreme cases
 
          delta = -ONEKB / temp0 * (my_ene_temp%energy_2 - my_ene_temp%energy_1) - &
                  ONEKB / neighbor_ene_temp%temperature * &
                  (neighbor_ene_temp%energy_2 - neighbor_ene_temp%energy_1)

          metrop = exp(delta + pvterm)
          call amrand_gen(remd_rand_gen, random_value)
          success = random_value .lt. metrop
      end if

      ! Let my neighbor know about our success
      call mpi_send(success, 1, mpi_logical, neighbor_rank, &
           remd_tag, remd_comm, err_code_mpi)
      if (success) then
#ifdef CUDA
        call gpu_download_vel(vel)
#endif
        success_array(replica_indexes(h_dim)) = 1
        ! swapping velocities after exchange succeed
        ! and scaling velocities   (DSD 09/12)
        call mpi_sendrecv_replace(vel, atm_cnt * 3, mpi_double_precision, &
                   neighbor_rank, remd_tag, &
                   neighbor_rank, remd_tag, &
                   remd_comm, istat, err_code_mpi)
        vel_scale_coff = sqrt(temp0 / neighbor_ene_temp%temperature)

        ! Exchange coordinate indices (control replica call)
        call mpi_sendrecv_replace(remd_crdidx, 1, mpi_integer, &
                                  neighbor_rank, remd_tag, &
                                  neighbor_rank, remd_tag, &
                                  remd_comm, istat, err_code_mpi)
      end if  ! if (success) 

    else
      
      call amrand_gen(remd_rand_gen, random_value) ! keep in sync
      call mpi_recv(success, 1, mpi_logical, neighbor_rank, &
                    remd_tag, remd_comm, istat, err_code_mpi)
      if (success) then
#ifdef CUDA
        call gpu_download_vel(vel)
#endif
        call mpi_sendrecv_replace(vel, atm_cnt * 3, mpi_double_precision, &
                      neighbor_rank, remd_tag, &
                      neighbor_rank, remd_tag, &
                      remd_comm, istat, err_code_mpi)
        vel_scale_coff = sqrt(temp0 / neighbor_ene_temp%temperature)
        ! Exchange coordinate indices (non-control replica call)
        call mpi_sendrecv_replace(remd_crdidx, 1, mpi_integer, &
                                  neighbor_rank, remd_tag, &
                                  neighbor_rank, remd_tag, &
                                  remd_comm, istat, err_code_mpi)
      end if

    end if ! i_do_exchg
    end if ! silent
    
    if (silent) then
        if (do_reservoir) then
            call amrand_gen(remd_rand_gen, random_value)
            if (gti_resv_accept .gt. 0) then 
               success=.true.
            else
                delta = -ONEKB / temp0 * (my_ene_temp%energy_2 - my_ene_temp%energy_1) 
                metrop = exp(delta + pvterm)
                success = random_value .lt. metrop
            endif
            if(success) enforce_reset_velocities=.true.

            enforce_reset_velocities=.false.
       endif
    endif
    
#ifdef VERBOSE_REMD
      write(mdout, '(a,f16.6)') '| My Eptot_1:       ', my_ene_temp%energy_1
      write(mdout, '(a,f16.6)') '| My Eptot_2:       ', my_ene_temp%energy_2
      write(mdout, '(a,f16.6)') '| Neighbor Eptot_1: ', neighbor_ene_temp%energy_1
      write(mdout, '(a,f16.6)') '| Neighbor Eptot_2: ', neighbor_ene_temp%energy_2
      if (i_do_exchg) then
        write(mdout, '(3(a,f16.6))') '| Delta= ', delta, ' Metrop= ', metrop, &
                                     ' Random #= ', random_value
        write(mdout,'(a7,i6,a8,i6)') '| RepIdx=', remd_repidx, ' CrdIdx=', remd_crdidx
      else
        write(mdout, '(3(a))') '| Not controlling exchange.'
      endif
      if (success) then
        write(mdout, '(a)') '| Exchange Succeeded!'
      else
        write(mdout, '(a)') '| Exchange Failed!'
      end if
#endif    

  end if ! master

  ! Now broadcast our success for all to hear! This is the ONLY part that is not
  ! for masters only

  call mpi_bcast(success, 1, mpi_logical, 0, pmemd_comm, err_code_mpi)
    
#ifdef CUDA
  ! If the exchange succeeded, then the coordinates and forces on the gpu
  ! are already correct and we don't need to do anything. Otherwise, we
  ! need to send back the originals.
  if (.not. success) then
    ! Upload the orginal coordinates again and force a new neighborlist if we're
    ! running PME
    if (using_pme_potential) call gpu_force_new_neighborlist()
    call gpu_upload_crd(crd_temp)
    call gpu_upload_frc(frc_bak)
  else
    if (.not.do_reservoir) then
      call mpi_bcast(vel, atm_cnt * 3, mpi_double_precision, 0, &
                     pmemd_comm, err_code_mpi)
      call mpi_bcast(vel_scale_coff, 1, mpi_double_precision, 0, &
                     pmemd_comm, err_code_mpi)
    else 
      if(.not. reservoir_velocity(tiRegion)) enforce_reset_velocities=.true.
    endif

    if(.not.enforce_reset_velocities) then
#ifdef CUDA
      call gpu_upload_vel(vel)
#endif
      call rescale_velocities(atm_cnt, vel, vel_scale_coff)
    endif 
  end if
#else
  ! Update our forces if we succeeded. Note that the exchange probability
  ! equation used currently satisfies detailed balance only for NOT exchanging
  ! temperatures. We also update our force array, since we won't have to run
  ! our force call on the first pass through our next MD loop in runmd
  if (success) then
    frc(1, :) = frc_temp(1, :)
    frc(2, :) = frc_temp(2, :)
    frc(3, :) = frc_temp(3, :)
    if (.not.do_reservoir) then
        call mpi_bcast(vel, atm_cnt * 3, mpi_double_precision, 0, &
                       pmemd_comm, err_code_mpi)
        call mpi_bcast(vel_scale_coff, 1, mpi_double_precision, 0, &
                       pmemd_comm, err_code_mpi)
        call rescale_velocities(atm_cnt, vel, vel_scale_coff)
    endif
  else
    ! Restore the coordinates and forces
    crd(1, :) = crd_temp(1, :)
    crd(2, :) = crd_temp(2, :)
    crd(3, :) = crd_temp(3, :)

    frc(1, :) = frc_bak(1, :)
    frc(2, :) = frc_bak(2, :)
    frc(3, :) = frc_bak(3, :)
  end if
#endif /* CUDA */
  if (exchange_ucell) then
    if (.not.success) then
      ! Restore original unit cell
      call remd_bcast_cell(pbc_box, original_box, pmemd_comm)
      !write(6,'(a,3f8.3)') '| DBG: Original box: ', pbc_box(1:3)
    else
      ! Update molecule centers of mass for new coords
      call remd_update_mol_com( crd )
    endif
  endif

  ! Calculate the free energies for jumping right and jumping left. This follows
  ! free energy perturbation (FEP). The equation for FEP is
  !
  !                          /        -(Eb - Ea)     \
  ! D G_(a->b) = -Kb * T ln | < exp (-------------) > |
  !                          \          Kb * T       /  A
  !
  ! Where the exponential is averaged. Thus, we have to keep a running total of
  ! the exponential Eb - Ea term (total_right_fe/total_left_fe), which can be
  ! combined into a total D G for this step via the above equation. These are
  ! stored in the ene_temp type to minimize communications. They can be
  ! calculated regardless of whether we succeeded or not!
  ! For the layout of the rem.log file: left = up, right = down.
  ! Then calculate the 'averaged' free energy right after, keeping in mind that
  ! right/left alternate. Therefore, we introduce a counter for the number of
  ! left-moves and number of right-moves to make sure we divide by the right
  ! number of samples.
  
  if(rremd_type.ge.4) then
    if(repnum.eq.num_replicas) &
    call mpi_send(rremd_idx(2), 1, mpi_integer, 0, &
                    remd_tag, remd_comm, stat_array, err_code_mpi)
    if(repnum.eq.1) &
    call mpi_recv(rremd_idx(2), 1, mpi_integer, num_replicas-1, &
                    remd_tag, remd_comm, stat_array, err_code_mpi)
  endif
  
  if (master) then
    if (i_do_exchg) then

      my_ene_temp%right_fe = &
                   exp((my_ene_temp%energy_1 - neighbor_ene_temp%energy_2) * &
                   ONEKB / my_ene_temp%temperature)
      my_ene_temp%right_exchg = 1.d0
      my_ene_temp%left_fe = 0.d0

    else

      my_ene_temp%left_fe = &
                   exp((my_ene_temp%energy_1 - neighbor_ene_temp%energy_2) * &
                   ONEKB / my_ene_temp%temperature)
      my_ene_temp%left_exchg = 1.d0
      my_ene_temp%right_fe = 0.d0

    end if

    ! Collect data and write it out to the rem.log

    call mpi_reduce(success_array, success_buf, num_replicas, mpi_integer, &
                    mpi_sum, 0, remd_comm, err_code_mpi)

    call mpi_gather(success, 1, mpi_logical, rank_success, 1, mpi_logical, &
                    0, remd_comm, err_code_mpi)

    call mpi_gather(my_ene_temp, SIZE_ENE_TEMP, mpi_double_precision, &
                    ene_temp_buffer, SIZE_ENE_TEMP, mpi_double_precision, &
                    0, remd_comm, err_code_mpi)

    if (remd_master) then

      exchange_successes(h_dim,:) = exchange_successes(h_dim,:) + success_buf(:)

      if (remd_method .ne. -1 .and. print_exch_data) then

        do i = 1, num_replicas
          num_left_exchg(h_dim, group_num(h_dim), i)  = &
                num_left_exchg(h_dim, group_num(h_dim), i) &
              + int(ene_temp_buffer(i)%left_exchg)

          num_right_exchg(h_dim, group_num(h_dim), i) = &
                num_right_exchg(h_dim, group_num(h_dim), i) &
              + int(ene_temp_buffer(i)%right_exchg)
        end do

        write(remlog, '(a,i8)') '# exchange ', mdloop

        do rep = 1, num_replicas

          ! Did we succeed or not? We need to print a character
          if (rank_success(rep)) then
            success_char = 'T'
          else
            success_char = 'F'
          end if

          i = ene_temp_buffer(rep)%repnum

          total_left_fe(h_dim, group_num(h_dim), i) = &
             total_left_fe(h_dim, group_num(h_dim), i) + &
             ene_temp_buffer(rep)%left_fe

          total_right_fe(h_dim, group_num(h_dim), i) = &
             total_right_fe(h_dim, group_num(h_dim), i) + &
             ene_temp_buffer(rep)%right_fe

          ! Calculate the success ratio:

          success_ratio = dble(exchange_successes(h_dim, &
                          int(ene_temp_buffer(rep)%repnum))) / &
                          dble(mdloop) * dble(remd_dimension) * 2

          ! Print info to the remlog

          if (num_left_exchg(h_dim, group_num(h_dim), rep) > 0) then
            l_fe = ene_temp_buffer(rep)%temperature / ONEKB * &
                   log(total_left_fe(h_dim, group_num(h_dim), rep) / &
                   num_left_exchg(h_dim, group_num(h_dim), rep))
          else
            l_fe = 0.d0
          end if

          if (num_right_exchg(h_dim, group_num(h_dim), rep) > 0) then
            r_fe = ene_temp_buffer(rep)%temperature / ONEKB * &
                   log(total_right_fe(h_dim, group_num(h_dim), rep) / &
                   num_right_exchg(h_dim, group_num(h_dim), rep))
          else
            r_fe = 0.d0
          end if

          if (rremd_type.ge.4) then
            i=0
            j=mod(mdloop,reservoir_exchange_step)
            if ((rremd_type.eq.4 .or. rremd_type.eq.6) .and. rep.eq.1 .and. j.eq.1) i=1
            if ((rremd_type.eq.5 .or. rremd_type.eq.6) .and. rep.eq.num_replicas) then
              if (mod(num_replicas,2) .eq. (1-j) ) i=2
            end if

            if (i.gt.0) then
              if (rremd_idx(i).gt.0) then
                write(remlog, '(2i6,5f15.7,4x,a,2x,f13.5, i8)') &
                  index_list(rep), &
                  int(ene_temp_buffer(rep)%neighbor_rep), &
                  ene_temp_buffer(rep)%temperature, &
                  ene_temp_buffer(rep)%energy_1,    &
                  ene_temp_buffer(rep)%energy_2,    &
                  l_fe,         &
                  r_fe,         &
                  success_char, &
                  success_ratio, rremd_idx(i)
              end if
            else
              write(remlog, '(2i6,5f15.7,4x,a,2x,f13.5)') &
                  index_list(rep), &
                  int(ene_temp_buffer(rep)%neighbor_rep), &
                  ene_temp_buffer(rep)%temperature, &
                  ene_temp_buffer(rep)%energy_1,    &
                  ene_temp_buffer(rep)%energy_2,    &
                  l_fe,         &
                  r_fe,         &
                  success_char, &
                  success_ratio
            endif
          else
            write(remlog, '(2i6,5f10.2,4x,a,2x,f10.2)') &
                  index_list(rep), &
                  int(ene_temp_buffer(rep)%neighbor_rep), &
                  ene_temp_buffer(rep)%temperature, &
                  ene_temp_buffer(rep)%energy_1,    &
                  ene_temp_buffer(rep)%energy_2,    &
                  l_fe,         &
                  r_fe,         &
                  success_char, &
                  success_ratio

          endif

        end do
      else
        do rep = 1, num_replicas

          success_ratio = dble(exchange_successes(h_dim, &
                          int(ene_temp_buffer(rep)%repnum))) /  &
                          dble(mdloop) * dble(remd_dimension) * 2
          if (rank_success(rep)) then
            multid_print_data(rep)%success = 1.d0
          else
            multid_print_data(rep)%success = 0.d0
          end if

          multid_print_data(rep)%neighbor_rep= ene_temp_buffer(rep)%neighbor_rep
          multid_print_data(rep)%temp0       = ene_temp_buffer(rep)%temperature
          multid_print_data(rep)%pot_ene_tot = ene_temp_buffer(rep)%energy_1
          multid_print_data(rep)%nei_pot_ene = ene_temp_buffer(rep)%energy_2
          multid_print_data(rep)%left_fe     = ene_temp_buffer(rep)%left_fe
          multid_print_data(rep)%right_fe    = ene_temp_buffer(rep)%right_fe
          multid_print_data(rep)%success_ratio= success_ratio
          multid_print_data(rep)%num_rep     = num_replicas
          multid_print_data(rep)%group_num   = group_num(h_dim)
          multid_print_data(rep)%repnum      = ene_temp_buffer(rep)%repnum
          multid_print_data(rep)%left_exchg  = ene_temp_buffer(rep)%left_exchg
          multid_print_data(rep)%right_exchg = ene_temp_buffer(rep)%right_exchg
        end do
      end if

    end if ! remd_master

#ifdef VERBOSE_REMD
    write(mdout, '(a,26("="),a,26("="))') '| ', 'END H-REMD CALCULATION'
#endif

  end if ! master

  remd_modwt = .true.

  return

end subroutine hamiltonian_exchange

!*******************************************************************************
!
! Subroutine: ph_remd_exchange
!
! Description: Performs the pH replica exchange attempt.
!
!*******************************************************************************

subroutine ph_remd_exchange(ph_dim, num_replicas, mdloop)

  use constantph_mod, only    : total_protonation
  use gbl_constants_mod, only : LN_TO_LOG
  use mdin_ctrl_dat_mod, only : solvph
!PHMD
  use mdin_ctrl_dat_mod, only : iphmd,solvph
  use phmd_mod, only          : double_protonation,reset_park

  use parallel_dat_mod

  implicit none

! Passed variables

  integer, intent(in) :: ph_dim
  integer, intent(in) :: num_replicas
  integer, intent(in) :: mdloop

! Local variables

  double precision :: delta                  ! DELTA value for MC transition
  double precision :: randval                ! Random #
  double precision :: all_ph(num_replicas)   ! pH table
  double precision :: new_ph(num_replicas)   ! table of pHs after exchange
  double precision :: exchfrac(num_replicas) ! fraction of successes
  double precision :: o_ph                   ! pH of neighbor replica
  double precision :: success_ratio          ! success rate
!PHMD
  double precision :: double_prot            ! sum of (1-lambda)
  double precision :: double_prot_table(num_replicas) ! table of double_prot
  double precision :: double_o_prot          ! sum of (1-lambda) in other replica

  integer :: i    ! counter
  integer :: prot ! # of protons in this replica
  integer :: suc_arry(num_replicas)   ! array with success values
  integer :: suc_buf(num_replicas)    ! reduce buffer for the above
  integer :: prot_table(num_replicas) ! table of `prot' for each replica
  integer :: my_index                 ! my position in pH table
  integer :: o_index                  ! index of replica we are exchanging with
  integer :: o_repnum                 ! replica number we're exchanging with
  integer :: o_prot                   ! # of protons in other replica

  logical :: success                  ! did we succeed?
  logical :: i_do_exchg               ! Do I do the exchange?

  integer, dimension(mpi_status_size) :: stat_array ! needed for mpi_recv

  ! Initialize
  suc_arry(:) = 0
  suc_buf(:) = 0
  success = .false.
  delta = 0.d0
  randval = 0.d0

  if (master) then

    call set_partners(ph_dim, num_replicas, remd_random_partner)

#ifdef VERBOSE_REMD
    write(6,'(a)') '| =============== REMD ==============='
#endif

    call mpi_allgather(solvph, 1, mpi_double_precision, &
                       all_ph, 1, mpi_double_precision, &
                       remd_comm, err_code_mpi)
    !PHMD
    if(iphmd .ne. 0) then
        double_prot = double_protonation()
        call mpi_allgather(double_prot, 1, mpi_double_precision, &
                           double_prot_table, 1, mpi_double_precision, &
                           remd_comm, err_code_mpi)
    else
       prot = total_protonation()

       call mpi_allgather(prot, 1, mpi_integer, prot_table, &
                          1, mpi_integer, remd_comm, err_code_mpi)
    endif

    my_index = replica_indexes(ph_dim)

    if (jumpright_ph) then

      if (mod(my_index, 2) .eq. 0) then
        o_index = my_index + 1
        if (o_index > num_replicas) o_index = 1
        o_repnum = partners(2)
      else
        o_index = my_index - 1
        if (o_index < 1) o_index = num_replicas
        o_repnum = partners(1)
      end if

    else ! not jumpright_ph

      if (mod(my_index, 2) .eq. 0) then
        o_index = my_index - 1
        if (o_index < 1) o_index = num_replicas
        o_repnum = partners(1)
      else
        o_index = my_index + 1
        if (o_index > num_replicas) o_index = 1
        o_repnum = partners(2)
      end if

    end if ! (jumpright_ph)

    o_ph = all_ph(o_repnum)
    !PHMD
    if(iphmd .ne. 0) then
       double_o_prot = double_prot_table(o_repnum)
    else
       o_prot = prot_table(o_repnum)
    endif

#ifdef VERBOSE_REMD
    !PHMD
    if(iphmd .NE. 0) then
      write(6,'(a,f5.2,a,f5.2,a,i3)') '| Found partner information: &
             &protonation = ', double_o_prot, ' pH = ', o_ph, ' repnum = ', o_repnum
    else
       write(6,'(a,i3,a,f5.2,a,i3)') '| Found partner information: &
                &protonation = ', o_prot, ' pH = ', o_ph, ' repnum = ', o_repnum
    endif

#endif

    if (mod(my_index, 2) .eq. 0) then

      call amrand_gen(remd_rand_gen, randval)
      !PHMD
      if(iphmd .NE. 0) then
         delta = LN_TO_LOG * (double_prot - double_o_prot) * (o_ph - solvph)
      else
         delta = LN_TO_LOG * (prot - o_prot) * (o_ph - solvph)
      endif

      success = randval < exp(-delta)

      ! Tell our neighbor if we succeeded

      call mpi_send(success, 1, mpi_logical, o_repnum-1, &
                    22, remd_comm, err_code_mpi)

      ! Increment our success

      if (success) then
        if (jumpright_ph) then
          suc_arry(my_index) = 1
        else
          suc_arry(o_index) = 1
        end if
      end if

#ifdef VERBOSE_REMD
       !PHMD
       if(iphmd .NE. 0) then
         write(6, '(2(a,f12.3))') '| Proton count: ', double_prot, ' --> ', double_o_prot
       else
          write(6, '(2(a,i3))') '| Proton count: ', prot, ' --> ', o_prot
       endif

       write(6, '(2(a,f8.3))') '| pH transition:', solvph, ' --> ', o_ph
       if (success) then
         write(6, '(a,E16.8)') '| Success! delta = ', delta
       else
         write(6, '(a,E16.8)') '| Failure. delta = ', delta
       end if
#endif

    else ! I do not exchange

      call mpi_recv(success, 1, mpi_logical, o_repnum-1, &
                    22, remd_comm, stat_array, err_code_mpi)
#ifdef VERBOSE_REMD
       !PHMD
       if(iphmd .NE. 0) then
         write(6, '(2(a,f12.3))') '| Proton count: ', double_prot, ' --> ', double_o_prot
       else
         write(6, '(2(a,i3))') '| Proton count: ', prot, ' --> ', o_prot
       endif
         write(6, '(2(a,f8.3))') '| pH transition:', solvph, ' --> ', o_ph
         if (success) then
            write(6, '(a)') '| Success!'
         else
            write(6, '(a)') '| Failure.'
         end if
#endif
    end if ! (mod(my_index, 2) .eq. 0)

    ! Update our solvph if necessary

    if (success) solvph = o_ph

      ! Swap the replica_indexes and group_num if our exchange was successful
    if (success .and. master) then
      call mpi_sendrecv_replace(replica_indexes, remd_dimension, &
                mpi_integer, o_repnum-1, 22, o_repnum-1, 22, &
                remd_comm, stat_array, err_code_mpi)
      call mpi_sendrecv_replace(group_num, remd_dimension, mpi_integer, &
                o_repnum-1, 23, o_repnum-1, 23, remd_comm, stat_array, &
                err_code_mpi)
      call mpi_sendrecv_replace(remd_repidx, 1, mpi_integer, &
                                o_repnum-1, 24, o_repnum-1, 24, &
                                remd_comm, stat_array, err_code_mpi)
    end if

    ! Reduce our success array to master for printing

    call mpi_reduce(suc_arry, suc_buf, num_replicas, mpi_integer, &
                    mpi_sum, 0, remd_comm, err_code_mpi)

    ! Generate our new pH table

    call mpi_gather(solvph, 1, mpi_double_precision, new_ph, 1, &
                    mpi_double_precision, 0, remd_comm, err_code_mpi)

  end if ! (master)

  ! Tell the rest of our replica about our (maybe) new pH

  call mpi_bcast(solvph, 1, mpi_double_precision, 0, pmemd_comm, err_code_mpi)

  !PHMD
  if(iphmd .NE. 0) then
     call reset_park()
  endif

  if (remd_master) &
    exchange_successes(ph_dim,:) = exchange_successes(ph_dim,:) + suc_buf(:)

  ! Write to the remlog

  ! remd_method == -1 means we're doing multi-D REMD.
  if (remd_method .ne. -1) then
    if (remd_master) then
      write(remlog, '(a,i8)') '# exchange ', mdloop
      do i = 1, num_replicas
        success_ratio = dble(exchange_successes(ph_dim, index_list(i))) / &
                        dble(mdloop) * dble(remd_dimension) * 2
        !PHMD
        if(iphmd .ne.0) then
           write(remlog, '(i6,x,3f7.3,x,f8.4)') &
                 i, &
                 double_prot_table(i), &
                 all_ph(i), &
                 new_ph(i), &
                 success_ratio
        else
           write(remlog, '(i6,x,i7,x,2f7.3,x,f8.4)') &
                 i, &
                 prot_table(i), &
                 all_ph(i), &
                 new_ph(i), &
                 success_ratio
        end if
      end do
    end if
  else if (remd_master) then
    do i = 1, num_replicas
      success_ratio = dble(exchange_successes(ph_dim, index_list(i))) / &
                      dble(mdloop) * dble(remd_dimension) * 2
      !PHMD
      if(iphmd .ne. 0) then
         multid_print_data(i)%nprot      = double_prot_table(i)
      else
         multid_print_data(i)%nprot      = prot_table(i)
      endif

      multid_print_data(i)%my_ph         = all_ph(i)
      multid_print_data(i)%nei_ph        = new_ph(i)
      multid_print_data(i)%group_num     = group_num(ph_dim)
      multid_print_data(i)%success_ratio = success_ratio
      multid_print_data(i)%num_rep       = num_replicas
    end do
  end if ! remd_method .ne. -1

  jumpright_ph = .not. jumpright_ph

end subroutine ph_remd_exchange
!*******************************************************************************
!
! Subroutine: e_remd_exchange
!
! Description: Performs the Redox potential replica exchange attempt.
!
!*******************************************************************************

subroutine e_remd_exchange(e_dim, num_replicas, mdloop)

  use constante_mod, only     : total_reduction
  use constantph_mod, only    : total_reduction2
  use gbl_constants_mod, only : KB, FARADAY
  use mdin_ctrl_dat_mod, only : solve, temp0
  use parallel_dat_mod
  use get_cmdline_mod, only   : cpein_specified

  implicit none

! Passed variables

  integer, intent(in) :: e_dim
  integer, intent(in) :: num_replicas
  integer, intent(in) :: mdloop

! Local variables

  double precision :: delta                  ! DELTA value for MC transition
  double precision :: randval                ! Random #
  double precision :: all_e(num_replicas)   ! Redox potential table
  double precision :: new_e(num_replicas)   ! table of Redox potentials after exchange
  double precision :: exchfrac(num_replicas) ! fraction of successes
  double precision :: o_e                   ! Redox potential of neighbor replica
  double precision :: success_ratio          ! success rate

  integer :: i    ! counter
  integer :: elec ! # of electrons in this replica
  integer :: suc_arry(num_replicas)   ! array with success values
  integer :: suc_buf(num_replicas)    ! reduce buffer for the above
  integer :: elec_table(num_replicas) ! table of `elec' for each replica
  integer :: my_index                 ! my position in Redox potential table
  integer :: o_index                  ! index of replica we are exchanging with
  integer :: o_repnum                 ! replica number we're exchanging with
  integer :: o_elec                   ! # of electrons in other replica

  logical :: success                  ! did we succeed?
  logical :: i_do_exchg               ! Do I do the exchange?

  integer, dimension(mpi_status_size) :: stat_array ! needed for mpi_recv

  ! Initialize
  suc_arry(:) = 0
  suc_buf(:) = 0
  success = .false.
  delta = 0.d0
  randval = 0.d0

  if (master) then

    call set_partners(e_dim, num_replicas, remd_random_partner)

#ifdef VERBOSE_REMD
    write(6,'(a)') '| =============== REMD ==============='
#endif

    call mpi_allgather(solve, 1, mpi_double_precision, &
                       all_e, 1, mpi_double_precision, &
                       remd_comm, err_code_mpi)

    if (.not. cpein_specified) then
      elec = total_reduction()
    else
      elec = total_reduction2()
    end if

    call mpi_allgather(elec, 1, mpi_integer, elec_table, &
                       1, mpi_integer, remd_comm, err_code_mpi)

    my_index = replica_indexes(e_dim)

    if (jumpright_e) then

      if (mod(my_index, 2) .eq. 0) then
        o_index = my_index + 1
        if (o_index > num_replicas) o_index = 1
        o_repnum = partners(2)
      else
        o_index = my_index - 1
        if (o_index < 1) o_index = num_replicas
        o_repnum = partners(1)
      end if

    else ! not jumpright_e

      if (mod(my_index, 2) .eq. 0) then
        o_index = my_index - 1
        if (o_index < 1) o_index = num_replicas
        o_repnum = partners(1)
      else
        o_index = my_index + 1
        if (o_index > num_replicas) o_index = 1
        o_repnum = partners(2)
      end if

    end if ! (jumpright_e)

    o_e = all_e(o_repnum)
    o_elec = elec_table(o_repnum)

#ifdef VERBOSE_REMD
    write(6,'(a,i3,a,f5.2,a,i3)') '| Found partner information: &
             &reduction = ', o_elec, ' E = ', o_e, ' repnum = ', o_repnum
#endif

    if (mod(my_index, 2) .eq. 0) then

      call amrand_gen(remd_rand_gen, randval)
      delta = ( FARADAY / ( KB * temp0 ) ) * (elec - o_elec) * (o_e - solve)
      success = randval < exp(-delta)

      ! Tell our neighbor if we succeeded

      call mpi_send(success, 1, mpi_logical, o_repnum-1, &
                    22, remd_comm, err_code_mpi)

      ! Increment our success

      if (success) then
        if (jumpright_e) then
          suc_arry(my_index) = 1
        else
          suc_arry(o_index) = 1
        end if
      end if

#ifdef VERBOSE_REMD
       write(6, '(2(a,i3))') '| Electron count: ', elec, ' --> ', o_elec
       write(6, '(2(a,f8.3))') '| E transition:', solve, ' --> ', o_e
       if (success) then
         write(6, '(a,E16.8)') '| Success! delta = ', delta
       else
         write(6, '(a,E16.8)') '| Failure. delta = ', delta
       end if
#endif

    else ! I do not exchange

      call mpi_recv(success, 1, mpi_logical, o_repnum-1, &
                    22, remd_comm, stat_array, err_code_mpi)
#ifdef VERBOSE_REMD
         write(6, '(2(a,i3))') '| Electron count: ', elec, ' --> ', o_elec
         write(6, '(2(a,f8.3))') '| E transition:', solve, ' --> ', o_e
         if (success) then
            write(6, '(a)') '| Success!'
         else
            write(6, '(a)') '| Failure.'
         end if
#endif
    end if ! (mod(my_index, 2) .eq. 0)

    ! Update our solve if necessary

    if (success) solve = o_e

      ! Swap the replica_indexes and group_num if our exchange was successful
    if (success .and. master) then
      call mpi_sendrecv_replace(replica_indexes, remd_dimension, &
                mpi_integer, o_repnum-1, 22, o_repnum-1, 22, &
                remd_comm, stat_array, err_code_mpi)
      call mpi_sendrecv_replace(group_num, remd_dimension, mpi_integer, &
                o_repnum-1, 23, o_repnum-1, 23, remd_comm, stat_array, &
                err_code_mpi)
      call mpi_sendrecv_replace(remd_repidx, 1, mpi_integer, &
                                o_repnum-1, 24, o_repnum-1, 24, &
                                remd_comm, stat_array, err_code_mpi)
    end if

    ! Reduce our success array to master for printing

    call mpi_reduce(suc_arry, suc_buf, num_replicas, mpi_integer, &
                    mpi_sum, 0, remd_comm, err_code_mpi)

    ! Generate our new Redox potential table

    call mpi_gather(solve, 1, mpi_double_precision, new_e, 1, &
                    mpi_double_precision, 0, remd_comm, err_code_mpi)

  end if ! (master)

  ! Tell the rest of our replica about our (maybe) new Redox potential

  call mpi_bcast(solve, 1, mpi_double_precision, 0, pmemd_comm, err_code_mpi)

  if (remd_master) &
    exchange_successes(e_dim,:) = exchange_successes(e_dim,:) + suc_buf(:)

  ! Write to the remlog

  ! remd_method == -1 means we're doing multi-D REMD.
  if (remd_method .ne. -1) then
    if (remd_master) then
      write(remlog, '(a,i8)') '# exchange ', mdloop
      do i = 1, num_replicas
        success_ratio = dble(exchange_successes(e_dim, index_list(i))) / &
                        dble(mdloop) * dble(remd_dimension) * 2
        write(remlog, '(i6,x,i7,x,2f7.3,x,f8.4)') &
              i, &
              elec_table(i), &
              all_e(i), &
              new_e(i), &
              success_ratio
      end do
    end if
  else if (remd_master) then
    do i = 1, num_replicas
      success_ratio = dble(exchange_successes(e_dim, index_list(i))) / &
                      dble(mdloop) * dble(remd_dimension) * 2
      multid_print_data(i)%nelec         = elec_table(i)
      multid_print_data(i)%my_e          = all_e(i)
      multid_print_data(i)%nei_e         = new_e(i)
      multid_print_data(i)%group_num     = group_num(e_dim)
      multid_print_data(i)%success_ratio = success_ratio
      multid_print_data(i)%num_rep       = num_replicas
    end do
  end if ! remd_method .ne. -1

  jumpright_e = .not. jumpright_e

end subroutine e_remd_exchange

!*******************************************************************************
!
! Subroutine: reservoir_exchange
!
! Description: Performs the reservoir temperature replica exchange attempt. We don't need
!              any extra energies here.
!
!*******************************************************************************

subroutine reservoir_exchange(atm_cnt, crd, vel, my_pot_ene_tot, t_dim, &
                                num_replicas, actual_temperature,    &
                                print_exch_data, mdloop)

  use mdin_ctrl_dat_mod, only : temp0, reservoir_exchange_step, using_pme_potential
  use pmemd_lib_mod,     only : mexit
  use reservoir_mod

  implicit none

! Passed variables

  integer, intent(in)             :: atm_cnt
  integer, intent(in)             :: num_replicas
  integer, intent(in)             :: mdloop

  double precision, intent(in)    :: my_pot_ene_tot
  double precision, intent(inout) :: crd(3, atm_cnt)
  double precision, intent(inout) :: vel(3, atm_cnt)
  double precision, intent(in)    :: actual_temperature
  
  integer, intent(in)             :: t_dim ! which dimension T-exchanges are
  logical, intent(in)             :: print_exch_data

! Local variables

  type :: exchange_data ! facilitates gather for remlog
    sequence
    double precision :: scaling
    double precision :: real_temp
    double precision :: pot_ene_tot
    double precision :: temp0
    double precision :: new_temp0
    double precision :: struct_num
  end type exchange_data

  integer, parameter :: SIZE_EXCHANGE_DATA = 6 ! for mpi_gather

  type (exchange_data) :: exch_buffer(num_replicas) ! gather buffer

  type (exchange_data) :: my_exch_data ! stores data for this replica
  type (exchange_data) :: exch_data_tbl(num_replicas) ! exch_data for all reps

  double precision  :: delta
  double precision  :: metrop
  double precision  :: random_value
  double precision  :: success_ratio

  integer           :: neighbor_rank
  integer           :: i
  integer           :: success_array(num_replicas)
  integer           :: success_buf(num_replicas)

  ! For exchanging replica positions in ladders upon successful exchanges
  integer           :: group_num_buffer(remd_dimension)
  integer           :: replica_indexes_buffer(remd_dimension)

  logical           :: success
  logical           :: i_do_exchg

  integer, dimension(mpi_status_size) :: stat_array ! needed for mpi_recv

  integer           :: my_idx ! Replica index 

  ! Used for broadcasting coordinates if exchange with reservoir is successful
  logical           :: highest_replica_reservoir_exchange_success

! Explanation of local variables:
!
! delta: The calculated DELTA value used as part of the detailed balance calc.
!
! metrop: The actual value compared to a random_value to decide success
!
! success: Did this exchange attempt succeed?
!
! success_array: buffer for keeping track of which temperatures exchanged
!                successfully. Used in a mpi_reduce at the end
!
! success_buf: buffer for mpi_reduce call on success_array
!
! neighbor_rank: the remd_rank of the replica with whom we need to be doing ALL
!                of our communication. Trying to do everything with the
!                'partners' array is a bit too confusing.

! Set the variables that we know at the beginning

  my_exch_data%real_temp   = actual_temperature
  my_exch_data%temp0       = temp0
  my_exch_data%new_temp0   = temp0 ! changed later if exch. succeeds
  my_exch_data%struct_num  = -1  
  my_exch_data%pot_ene_tot = my_pot_ene_tot
  my_exch_data%scaling     = -1.d0
  success_array(:)         = 0

  ! Need coordinates for swapping structures from the reservoir
  if (master) then
#ifdef CUDA
        call gpu_download_crd(crd)
        call gpu_download_vel(vel)
#endif
  end if

  ! Set below variable to .false. initially. It will be set to .true. if
  ! exchange with reservoir is successful. If .true., new coordinates and
  ! velocities will be broadcasted. Force will be calculated for the new 
  ! structure.
  highest_replica_reservoir_exchange_success = .FALSE.
  call mpi_bcast(highest_replica_reservoir_exchange_success, 1, mpi_logical, &
                 0, pmemd_comm, err_code_mpi)

  if (master) then

    call set_partners(t_dim, num_replicas, remd_random_partner)

    ! Call amrand_gen here to generate a random number. Sander does this
    ! to pick out structures from a reservoir which pmemd doesn't implement
    ! (yet). To keep the random #s synched, though, and to ensure proper
    ! regression testing, this must be done to keep the random #s synched
    ! b/w sander and pmemd

    call amrand_gen(remd_rand_gen, random_value)

#ifdef VERBOSE_REMD
    write(mdout, '(26("="),a,26("="))') 'REMD EXCHANGE CALCULATION'
    write(mdout, '(a,i10,a,i1)') 'Exch= ', mdloop, ' RREMD= ', rremd_type
#endif

    ! We alternate exchanges: temp up followed by temp down
    if (even_exchanges(t_dim)) then
      i_do_exchg = control_exchange(t_dim)
    else
      i_do_exchg = .not. control_exchange(t_dim)
    end if

    even_exchanges(t_dim) = .not. even_exchanges(t_dim)

    ! Find out what rank our neighbor is (-1 for indexing from 0), and
    ! If we are controlling the exchange, then our partner will be our 
    ! 2nd partner (above). Otherwise it will be our 1st partner (below)

    if (i_do_exchg) then
      neighbor_rank = partners(2) - 1
    else
      neighbor_rank = partners(1) - 1
    end if

    ! Collect potential energies/temperatures

    call mpi_allgather(my_exch_data, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                exch_data_tbl, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                remd_comm, err_code_mpi)

#ifdef VERBOSE_REMD
    write(mdout, '(a,f6.2,2(a7,i2),a7,f10.2)') 'Replica          Temp= ', &
      my_exch_data%temp0, ' Indx= ', replica_indexes(t_dim), ' Rep#= ', &
      remd_rank+1, ' EPot= ', my_pot_ene_tot
    write(mdout,'(a7,i6,a8,i6)') 'RepIdx=', remd_repidx, ' CrdIdx=', remd_crdidx
#endif

    ! Get the replica indexes. my_idx will be equal to 1 for replica with 
    ! lowest T and equal to num_replicas for replica with highest T
    my_idx = replica_indexes(remd_dimension)
    
    ! Initialize rremd_idx, so that it is -1 for any replica other than the 
    ! one that will use the reservoir coordinates.
    rremd_idx = -1

    ! Calculate the exchange probability. Note that if delta is negative, then
    ! metrop will be > 1, and it will always succeed (so we don't need an extra
    ! conditional)

    ! Attempt exchange with reservoir if mod(mdloop, reservoir_exchange_step) is 1.
    ! Because the number of replicas should be even, the highest replica will always 
    ! have an even index. The even replicas control the exchange (i_do_exchg) when 
    ! the exchange number is odd. Therefore for the highest replica to exchange with 
    ! the reservoir, the exchange number has to be odd and hence the modulus should
    ! be equal to 1 and not zero. Since we attempt exchange with the reservoir every
    ! alternate exchange, the reservoir_exchange_step has to be a multiple of 2.
    if(mod(mdloop, reservoir_exchange_step) .eq. 1) then
       if (i_do_exchg) then
       ! Reservoir REMD: If replica is doing exchange and has highest T, then attempt an
       ! exchange with a random structure in the reservoir.
          if (rremd_type .gt. 0) then
              if (my_idx .eq. num_replicas) then
                ! Store reservoir temperature
                exch_data_tbl(neighbor_rank+1)%temp0 = reservoir_temperature(1)
                ! Get random structure index
                ! amrand has already been called at the start of the routine so
                ! that all threads stay in sync.
                rremd_idx = int(random_value * reservoir_size(1)) + 1
                ! Random number should not be outside reservoir size!
                ! By definition it won't be unless something breaks in amrand()
                if (rremd_idx(1) .gt. reservoir_size(1)) then
                    rremd_idx(1) = reservoir_size(1)
                end if
                ! Get Potential Energy for the corresponding reservoir structure index (rremd_idx)
                ! The routine load_reservoir_files() in remd.F90 file has already loaded the 
                ! energies in reservoir_structure_energies_array. So, we just read it from there.
                exch_data_tbl(neighbor_rank+1)%pot_ene_tot = &
                   reservoir_structure_energies_array(rremd_idx(1))
              end if
          end if
#ifdef VERBOSE_REMD
          if (rremd_type .gt. 0) then
            if (my_idx .eq. num_replicas) then
            ! Partner is Reservoir
              write(mdout, '(a16, a7, f6.2, a10, i8, a7, f10.2)') &
                 "Reservoir      ", " Temp= ", exch_data_tbl(neighbor_rank+1)%temp0, &
                 " Struct#= ", rremd_idx, &
                 " EPot= ", exch_data_tbl(neighbor_rank+1)%pot_ene_tot
            else
              write(mdout, '(a,f6.2,2(a7,i2),a7,f10.2)') 'Partner          Temp= ', &
                exch_data_tbl(neighbor_rank+1)%temp0, ' Indx= ', neighbor_rank+1, &
                ' Rep#= ', neighbor_rank + 1, ' EPot= ', &
                exch_data_tbl(neighbor_rank+1)%pot_ene_tot
            end if
          else ! Not Reservoir REMD
            write(mdout, '(a,f6.2,2(a7,i2),a7,f10.2)') 'Partner          Temp= ', &
              exch_data_tbl(neighbor_rank+1)%temp0, ' Indx= ', neighbor_rank+1, &
              ' Rep#= ', neighbor_rank + 1, ' EPot= ', &
              exch_data_tbl(neighbor_rank+1)%pot_ene_tot
          end if
#endif

          if(rremd_type .eq. 2) then
             if(my_idx .eq. num_replicas) then
                ! No delta beta for exchange with reservoir for highest replica
                ! assign weight of 1/N to each structure in reservoir.
                ! The exchange criterion becomes:
                ! metrop = exp[-beta_replica * (E_reservoir-E_replica)]
                ! Need -1 to switch sign for Non-Boltzmann
                delta = -ONE * (my_exch_data%pot_ene_tot - &
                         exch_data_tbl(neighbor_rank+1)%pot_ene_tot) * ONEKB /  &
                         (my_exch_data%temp0)
             else
                delta = (my_exch_data%pot_ene_tot - &
                         exch_data_tbl(neighbor_rank+1)%pot_ene_tot) * &
                        (my_exch_data%temp0 - &
                         exch_data_tbl(neighbor_rank+1)%temp0) * ONEKB /  &
                        (my_exch_data%temp0 * exch_data_tbl(neighbor_rank+1)%temp0)
             end if

          else if(rremd_type .eq. 1) then
            ! All replicas use delta beta including the highest temperature replica.
            delta = (my_exch_data%pot_ene_tot - &
                     exch_data_tbl(neighbor_rank+1)%pot_ene_tot) * &
                    (my_exch_data%temp0 - &
                     exch_data_tbl(neighbor_rank+1)%temp0) * ONEKB /  &
                    (my_exch_data%temp0 * exch_data_tbl(neighbor_rank+1)%temp0)

          else
            ! Not Reservoir REMD, same as rremd_type .eq. 1. Separated here for clarity.
            ! This should not be used anyway, because Reservoir REMD requires either
            ! rremd_type .eq. 1 or rremd_type .eq. 2.
            delta = (my_exch_data%pot_ene_tot - &
                     exch_data_tbl(neighbor_rank+1)%pot_ene_tot) * &
                    (my_exch_data%temp0 - &
                     exch_data_tbl(neighbor_rank+1)%temp0) * ONEKB /  &
                    (my_exch_data%temp0 * exch_data_tbl(neighbor_rank+1)%temp0)
          end if

          metrop = exp(-delta)

          call amrand_gen(remd_rand_gen, random_value)

          success = random_value .lt. metrop

          ! Let my neighbor know about our success

          ! In sander, for rremd_type .gt. 0, the "success" is only sent for
          ! the replicas in the middle. The replicas with highest and lowest 
          ! T do not communicate. Here, we use a different scheme. We send 
          ! "success" for all replicas and then change the temperature and scaling 
          ! for replicas with highest and lowest T to keep their current values
          ! (see below). 

          call mpi_send(success, 1, mpi_logical, neighbor_rank, &
                        remd_tag, remd_comm, err_code_mpi)

          if (success) then
            success_array(replica_indexes(t_dim)) = 1
            my_exch_data%new_temp0 = exch_data_tbl(neighbor_rank+1)%temp0
            my_exch_data%scaling = sqrt(my_exch_data%new_temp0 / my_exch_data%temp0)

            ! Do not change temperature of the replica with highest T if it 
            ! exchanged with the reservoir. So new_temp0 should be equal to temp0.
            if (rremd_type .gt. 0 .and. my_idx .eq. num_replicas) then

               highest_replica_reservoir_exchange_success = .TRUE.

               my_exch_data%new_temp0 = my_exch_data%temp0
            ! If Reservoir REMD and exchanging with the reservoir we are changing 
            ! structure and not temperature, so scaling is inverted.
               if (reservoir_velocity(1)) then
                   my_exch_data%scaling = 1.d0/my_exch_data%scaling
            ! Don't scale if we are not reading velocities in the reservoir.
               else
                   my_exch_data%scaling = 1.d0
               end if
            ! Load the reservoir structure
               call load_reservoir_structure(crd, vel, atm_cnt)
            ! Make success .FALSE.. Otherwise, when mpi_gather() is called below 
            ! for my_exch_data, having success = .TRUE. for the highest T
            ! effects the temperature and replica indices for the 
            ! subsequent exchanges. This is a simple fix. The same is done for 
            ! the replica with the lowest T (see below).
               success = .FALSE.
            end if

          end if

#ifdef VERBOSE_REMD
          write(mdout,'(a8,E16.6,a8,E16.6,a12,f10.2)') &
                   "Metrop= ",metrop," delta= ",delta," o_scaling= ", &
                   1 / my_exch_data%scaling
#endif

       else

          call amrand_gen(remd_rand_gen, random_value) ! to stay synchronized

#ifdef VERBOSE_REMD
          ! Even if we don't control exchange, still print REMD info to mdout

          if (rremd_type .gt. 0) then
            if (my_idx .eq. 1) then
            ! No partner for lowest temperature replica
              write(mdout, '(a)') &
                 ' No partner. Highest T replica exchanging with reservoir'
              write(mdout, '(a)') 'Not controlling exchange.'
            else
              write(mdout, '(a,f6.2,2(a7,i2),a7,f10.2)') 'Partner          Temp= ', &
                exch_data_tbl(neighbor_rank+1)%temp0, ' Indx= ', neighbor_rank+1, &
                ' Rep#= ', neighbor_rank + 1, ' EPot= ', &
                exch_data_tbl(neighbor_rank+1)%pot_ene_tot
              write(mdout, '(a)') 'Not controlling exchange.'
            end if
          else ! Not Reservoir REMD
            write(mdout, '(a,f6.2,2(a7,i2),a7,f10.2)') 'Partner          Temp= ', &
              exch_data_tbl(neighbor_rank+1)%temp0, ' Indx= ', neighbor_rank+1, &
              ' Rep#= ', neighbor_rank + 1, ' EPot= ', &
              exch_data_tbl(neighbor_rank+1)%pot_ene_tot
            write(mdout, '(a)') 'Not controlling exchange.'

          end if
#endif

          ! Get the message from the exchanging replica about our success

          ! For reservoir REMD, even the replica with lowest T receives 
          ! success when replica with highest T exchanges with reservoir.
          ! We modify this success for the replica with lowest T below.
          ! This is different from sander. In sander, the replica with
          ! lowest T does not receive "success" at all.
          call mpi_recv(success, 1, mpi_logical, neighbor_rank, &
                        remd_tag, remd_comm, stat_array, err_code_mpi)

          ! We scale velocities only if we succeed
          if (success) then
            my_exch_data%new_temp0 = exch_data_tbl(neighbor_rank+1)%temp0
            my_exch_data%scaling = sqrt(my_exch_data%new_temp0 / my_exch_data%temp0)
            ! When replica with highest T exchanges with reservoir, the replica with 
            ! lowest T should retain the same temperature.
            if(rremd_type .gt. 0 .and. my_idx .eq. 1) then
               my_exch_data%new_temp0 = my_exch_data%temp0
               my_exch_data%scaling = 1.d0
              ! Make success .FALSE.. Otherwise, when mpi_gather() is called below 
              ! for my_exch_data, having success = .TRUE. for the highest T
              ! effects the temperature and replica indices for the 
              ! subsequent exchanges. This is a simple fix. The same is done for 
              ! the replica with the lowest T (see below).
              success = .FALSE.
            end if
          end if

       end if

#ifdef VERBOSE_REMD
       if (rremd_type .gt. 0) then
         if (my_idx .eq. num_replicas) then
           if(i_do_exchg) then
             if (highest_replica_reservoir_exchange_success .eqv. .TRUE.) then
               write(mdout, '(a)') 'Reservoir Exchange = T'
             else
               write(mdout, '(a)') 'Reservoir Exchange = F'
             end if

             write(mdout,'(a8,E16.6,a12,f10.2,a10,L1)') &
               'Rand=   ', random_value, ' MyScaling= ', my_exch_data%scaling, &
               ' Success= ', highest_replica_reservoir_exchange_success
           else
             write(mdout,'(a8,E16.6,a12,f10.2,a10,L1)') &
               'Rand=   ', random_value, ' MyScaling= ', my_exch_data%scaling, &
               ' Success= ', success
           end if

         else
           write(mdout,'(a8,E16.6,a12,f10.2,a10,L1)') &
             'Rand=   ', random_value, ' MyScaling= ', my_exch_data%scaling, &
             ' Success= ', success

         end if
       end if

       write(mdout,'(24("="),a,24("="))') "END RESERVOIR REMD EXCHANGE CALCULATION"
#endif

       my_exch_data%struct_num = rremd_idx(1) ! Store rremd_idx in struct_num

    else ! Don't attempt exchange with reservoir
       if (i_do_exchg) then
#ifdef VERBOSE_REMD
          write(mdout, '(a,f6.2,2(a7,i2),a7,f10.2)') 'Partner          Temp= ', &
            exch_data_tbl(neighbor_rank+1)%temp0, ' Indx= ', neighbor_rank+1, &
            ' Rep#= ', neighbor_rank + 1, ' EPot= ', &
            exch_data_tbl(neighbor_rank+1)%pot_ene_tot
#endif

          delta = (my_exch_data%pot_ene_tot - &
                   exch_data_tbl(neighbor_rank+1)%pot_ene_tot) * &
                  (my_exch_data%temp0 - &
                   exch_data_tbl(neighbor_rank+1)%temp0) * ONEKB /  &
                  (my_exch_data%temp0 * exch_data_tbl(neighbor_rank+1)%temp0)

          metrop = exp(-delta)

          call amrand_gen(remd_rand_gen, random_value)

          success = random_value .lt. metrop

          ! Let my neighbor know about our success

          call mpi_send(success, 1, mpi_logical, neighbor_rank, &
                        remd_tag, remd_comm, err_code_mpi)

          if (success) then
            success_array(replica_indexes(t_dim)) = 1
            my_exch_data%new_temp0 = exch_data_tbl(neighbor_rank+1)%temp0
            my_exch_data%scaling = sqrt(my_exch_data%new_temp0 / my_exch_data%temp0)
          end if

#ifdef VERBOSE_REMD
          write(mdout,'(a8,E16.6,a8,E16.6,a12,f10.2)') &
                   "Metrop= ",metrop," delta= ",delta," o_scaling= ", &
                   1 / my_exch_data%scaling
#endif

       else

#ifdef VERBOSE_REMD
      ! Even if we don't control exchange, still print REMD info to mdout

          write(mdout, '(a,f6.2,2(a7,i2),a7,f10.2)') 'Partner          Temp= ', &
            exch_data_tbl(neighbor_rank+1)%temp0, ' Indx= ', neighbor_rank+1, &
            ' Rep#= ', neighbor_rank + 1, ' EPot= ', &
            exch_data_tbl(neighbor_rank+1)%pot_ene_tot
          write(mdout, '(a)') 'Not controlling exchange.'
#endif

          call amrand_gen(remd_rand_gen, random_value) ! to stay synchronized

          ! Get the message from the exchanging replica about our success
          call mpi_recv(success, 1, mpi_logical, neighbor_rank, &
                        remd_tag, remd_comm, stat_array, err_code_mpi)

          ! We scale velocities only if we succeed
          if (success) then
            my_exch_data%new_temp0 = exch_data_tbl(neighbor_rank+1)%temp0
            my_exch_data%scaling = sqrt(my_exch_data%new_temp0 / my_exch_data%temp0)
          end if

       end if

#ifdef VERBOSE_REMD
       write(mdout,'(a8,E16.6,a12,f10.2,a10,L1)') &
         'Rand=   ', random_value, ' MyScaling= ', my_exch_data%scaling, &
         ' Success= ', success

       write(mdout,'(24("="),a,24("="))') "END RESERVOIR REMD EXCHANGE CALCULATION"
#endif

    end if ! End if(mod(mdloop,reservoir_exchange_step) .eq. 1)

  end if ! master


  ! We'll broadcast all of the my_exch_data here once. From that, we can tell
  ! if the exchange succeeded based on whether or not the temperatures have
  ! changed. Then make sure we update temp0 so the simulation reflects that.

  call mpi_bcast(my_exch_data, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                 0, pmemd_comm, err_code_mpi)

  temp0 = my_exch_data%new_temp0

  if (.not. master) &
    success = abs(my_exch_data%temp0 - my_exch_data%new_temp0) .gt. TINY

#ifdef VERBOSE_REMD
  if (master) then
     if (highest_replica_reservoir_exchange_success .eqv. .TRUE.) then
         write(mdout, '(a)') &
         'Exchange with reservoir successful. Do not change T'
     else
        write(mdout, '(2(a,f6.2))') &
        'REMD: checking to see if bath T has changed: ', &
        my_exch_data%temp0, '->', my_exch_data%new_temp0
     end if
  end if
#endif

  ! Upload/bcast new coordinates and velocities if exchanged with the reservoir.
  call mpi_bcast(highest_replica_reservoir_exchange_success, 1, mpi_logical, &
                  0, pmemd_comm, err_code_mpi)
  if(highest_replica_reservoir_exchange_success .eqv. .TRUE.) then        
#ifdef CUDA
        call gpu_upload_crd(crd)
        call gpu_upload_vel(vel)
        if (using_pme_potential) call gpu_force_new_neighborlist()
#endif
      ! The coordinates of the structure obtained
      ! from the reservoir should be broadcasted to the slave process. Otherwise,
      ! the slave process for each replica will be using the second-half of 
      ! the coordinates that were present before the exchange with reservoir
      ! and the master process for each replica will be using the first-half
      ! of coordinates of the structure obtained from the reservoir. If this
      ! is not corrected, it will lead to vlimit errors during MD part of REMD.

      ! Now master needs to broadcast coordinates and velocity for reservoir
      ! to all processes in its group.

        call mpi_bcast(crd, atm_cnt*3, mpi_double_precision, 0, &
                       pmemd_comm, err_code_mpi)
        call mpi_bcast(vel, atm_cnt*3, mpi_double_precision, 0, &
                       pmemd_comm, err_code_mpi)
  end if

  if (success) call rescale_velocities(atm_cnt, vel, my_exch_data%scaling)

  ! We want to keep track of our exchange successes, so gather them
  ! from all of the masters here, and increment them in exchange_successes
  ! We do this whether or not we print so we can keep track of exchange
  ! success ratios. For the index and exchange data, only gather those to
  ! the remd_master if we actually plan to print. This should be much
  ! more efficient for high exchange attempts.

  ! Then recollect the temperatures as some may have changed, and reset our
  ! partners.

  if (master) then

    if (print_exch_data) then
      call mpi_gather(my_exch_data, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                      exch_buffer, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                      0, remd_comm, err_code_mpi)
    end if

    call mpi_reduce(success_array, success_buf, num_replicas, mpi_integer, &
                    mpi_sum, 0, remd_comm, err_code_mpi)

    if (remd_master) &
      exchange_successes(t_dim,:) = exchange_successes(t_dim,:) + success_buf(:)

    ! Increment exchange counter, and swap our replica ranks and numbers with 
    ! our neighbor if our attempt succeeded, since we effectively swapped places
    ! in this array. That's because we swap temperatures and not structures due
    ! to performance considerations.

    if (success) then
      call mpi_sendrecv(group_num, remd_dimension, mpi_integer, &
                        neighbor_rank, remd_tag, &
                        group_num_buffer, remd_dimension, mpi_integer, &
                        neighbor_rank, remd_tag, remd_comm, &
                        stat_array, err_code_mpi)
      call mpi_sendrecv(replica_indexes, remd_dimension, mpi_integer, &
                        neighbor_rank, remd_tag, &
                        replica_indexes_buffer, remd_dimension, mpi_integer, &
                        neighbor_rank, remd_tag, remd_comm, &
                        stat_array, err_code_mpi)
      call mpi_sendrecv_replace(remd_repidx, 1, mpi_integer, &
                                neighbor_rank, remd_tag, &
                                neighbor_rank, remd_tag, &
                                remd_comm, stat_array, err_code_mpi)
      group_num = group_num_buffer
      replica_indexes = replica_indexes_buffer
    end if
  end if

  ! If we're doing NMROPT, then we need to make sure we re-read in the TEMP0,
  ! since this subroutine overwrites TEMP0 to the original since you can modify
  ! temperature with NMROPT settings.

  remd_modwt = .true.

  ! Write the data to the remlog

  ! remd_method == -1 means we're doing multi-D REMD.
  if (remd_method .ne. -1) then
    if (print_exch_data .and. remd_master) then
      write(remlog, '(a,i8)') '# exchange ', mdloop
      do i = 1, num_replicas
        success_ratio = dble(exchange_successes(t_dim, index_list(i))) / &
                        dble(mdloop) * dble(remd_dimension) * 2
        write(remlog, '(i2,f12.4,f10.2,f15.2,3f10.2,i8)') &
              i, &
              exch_buffer(i)%scaling, &
              exch_buffer(i)%real_temp, &
              exch_buffer(i)%pot_ene_tot, &
              exch_buffer(i)%temp0, &
              exch_buffer(i)%new_temp0, &
              success_ratio, &
              int(exch_buffer(i)%struct_num)
      end do
    end if
  else if (remd_master) then
    do i = 1, num_replicas
      success_ratio = dble(exchange_successes(t_dim, index_list(i))) / &
                      dble(mdloop) * dble(remd_dimension) * 2
      multid_print_data(i)%scaling       = exch_buffer(i)%scaling
      multid_print_data(i)%real_temp     = exch_buffer(i)%real_temp
      multid_print_data(i)%pot_ene_tot   = exch_buffer(i)%pot_ene_tot
      multid_print_data(i)%temp0         = exch_buffer(i)%temp0
      multid_print_data(i)%new_temp0     = exch_buffer(i)%new_temp0
      multid_print_data(i)%group_num     = group_num(t_dim)
      multid_print_data(i)%success_ratio = success_ratio
      multid_print_data(i)%num_rep       = num_replicas
    end do
  end if ! remd_method .ne. -1

  return

end subroutine reservoir_exchange

!*******************************************************************************
!
! Subroutine: hybridsolvent_exchange
!
! Description: Calculate energy of "solute + numwatkeep waters" system in GB
!              and save it in GB.
!
!              First, it strips the molecule to keep the desired number of atoms.
!              Then, calls gb_hybridsolvent_remd_ene() in gb_force.F90 to 
!              calculate the energy.
!              Also, if VERBOSE_REMD is set, calls hybridsolvent_remd_writeenergies() 
!              and hybridsolvent_remd_writetraj() to write energies and coordinates,
!              respectively.
!
!*******************************************************************************

subroutine hybridsolvent_exchange(atm_cnt, crd, vel, frc, numwatkeep, &
                                my_pot_ene_tot, t_dim, num_replicas, &
                                actual_temperature, print_exch_data, mdloop)

  use pmemd_lib_mod,     only : mexit
  use gb_force_mod
  use runfiles_mod,      only : prntmd, corpac
  use hybridsolvent_remd_mod
  use file_io_mod

  implicit none

! Passed variables

  integer, intent(in)             :: atm_cnt
  integer, intent(in)             :: num_replicas
  integer, intent(in)             :: mdloop
  integer, intent(in)             :: numwatkeep

  double precision, intent(in)    :: my_pot_ene_tot
  double precision, intent(inout) :: crd(3, atm_cnt)
  double precision, intent(inout) :: vel(3, atm_cnt)
  double precision, intent(inout) :: frc(3, atm_cnt)
  double precision, intent(in)    :: actual_temperature
  
  integer, intent(in)             :: t_dim ! which dimension T-exchanges are
  logical, intent(in)             :: print_exch_data

! Local variables
  type(gb_pot_ene_rec) :: my_gb_pot_ene_tot
  double precision  :: hybridsolvent_remd_crd(3, atm_cnt) ! Coordinates of stripped system
  double precision  :: hybridsolvent_remd_frc(3, atm_cnt) ! Coordinates of stripped system

  integer :: i, j

! Variables for writing coordinates of stripped system  
  integer                   :: hybridsolvent_remd_crd_cnt ! 3 * hybridsolvent_remd_atm_cnt
  integer                   :: local_total_nstep
  integer                   :: my_mdloop


  ! Need coordinates for building stripped system
  if (master) then
#ifdef CUDA
        call gpu_download_crd(crd)
#endif
  end if

  ! The code below strips the water molecules, gets energy of stripped system
  ! 1. Copying coordinates to hybridsolvent_crd and forces to hybridsolvent_frc. We use these
  ! two arrays for the rest of the calculation so that the original coordinates
  ! are kept intact.
  do i = 1, 3
    do j =1, atm_cnt
      hybridsolvent_remd_crd(i,j) = crd(i,j)
      hybridsolvent_remd_frc(i,j) = 0.d0
    end do
  end do

#ifdef VERBOSE_REMD
  if (master) then
    write(mdout,'(a)')
    write(mdout,'(20("="),a,20("="))') "BEGIN HYBRID SOLVENT REMD EXCHANGE"
    write(mdout,'(a)') 'Hybrid Solvent REMD: Stripping waters'
  end if
#endif

  ! 2. Strip the water molecules.
  hybridsolvent_remd_atm_cnt = atm_cnt
  call hybridsolvent_remd_strip_water(hybridsolvent_remd_atm_cnt, hybridsolvent_remd_crd, numwatkeep)

#ifdef VERBOSE_REMD    
  if (master) then
    write(6,'(a,i8)') "Done stripping atoms. Hybrid Solvent REMD: &
                       Stripped system atm_cnt= ", hybridsolvent_remd_atm_cnt
  end if
#endif

  ! 3. Broadcast or update the coordinates of the new system
!#ifdef CUDA 
!  ! Upload the coordinates and forces of the stripped system onto the GPU
!  call gpu_update_natoms(hybridsolvent_remd_atm_cnt, .false.)
!  call gpu_upload_crd_gb_cph(hybridsolvent_remd_crd)
!  call gpu_upload_frc(hybridsolvent_remd_frc)
!#else
  ! Broadcast hybridsolvent_remd_atm_cnt first so that all threads know how many coordinates to receive
  call mpi_bcast(hybridsolvent_remd_atm_cnt, 1, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)
  ! Broadcast the coordinates for the stripped system
  call mpi_bcast(hybridsolvent_remd_crd, hybridsolvent_remd_atm_cnt*3, mpi_double_precision, 0, &
                 pmemd_comm, err_code_mpi)
  ! Broadcast the forces for the stripped system
  call mpi_bcast(hybridsolvent_remd_frc, hybridsolvent_remd_atm_cnt*3, mpi_double_precision, 0, &
                 pmemd_comm, err_code_mpi)
!#endif

  ! 4. Calculate the energy of the stripped system. The new energy is saved in
  ! my_gb_pot_ene_tot. It's a energy data structure. Using my_gb_pot_ene_tot%total
  ! will give the total potential energy for the stripped system.

#ifdef VERBOSE_REMD
  if (master) then
    write(mdout,'(a)') 'Hybrid Solvent REMD: Calling gb_hybridsolvent_remd_ene() in gb_force.F90'
  end if
#endif

  call gb_hybridsolvent_remd_ene(hybridsolvent_remd_atm_cnt, hybridsolvent_remd_crd, &
                            hybridsolvent_remd_frc, my_gb_pot_ene_tot, 1, .true.)

#ifdef VERBOSE_REMD
  if (master) then
          write(mdout,'(a)') 'Hybrid Solvent REMD: Done calling gb_hybridsolvent_remd_ene()'
  end if
#endif

!#ifdef CUDA
!  call gpu_update_natoms(atm_cnt, .true.)
!  call gpu_upload_crd(crd)
!#endif

#ifdef VERBOSE_REMD
  if (master) then
    write(mdout,'(a)')
    write(mdout, '(30("="),a,30("="))') 'ENERGY FROM HYBRID SOLVENT'
    write (mdout, '(a,f13.2)')' EPtot      = ', my_gb_pot_ene_tot%total
    write (mdout, '(a,f14.4,a,f14.4,a,f14.4)') ' BOND   = ', my_gb_pot_ene_tot%bond, &
       '  ANGLE   = ', my_gb_pot_ene_tot%angle, '  DIHED      = ', my_gb_pot_ene_tot%dihedral
    write (mdout, '(a,f14.4,a,f14.4,a,f14.4)') ' 1-4 NB = ', my_gb_pot_ene_tot%vdw_14,&
       '  1-4 EEL = ', my_gb_pot_ene_tot%elec_14, '  VDWAALS    = ', my_gb_pot_ene_tot%vdw_tot
    write (mdout, '(a,f14.4,a,f14.4,a,f14.4)') ' EELEC  = ', my_gb_pot_ene_tot%elec_tot, &
       '  EGB     = ', my_gb_pot_ene_tot%gb, '  ESURF  = ', my_gb_pot_ene_tot%surf
    write (mdout, '(a,f14.4,a,f14.4)') '  RESTRAINT  = ', my_gb_pot_ene_tot%restraint, &
       ' EAMBER (non-restraint)  = ', my_gb_pot_ene_tot%total - my_gb_pot_ene_tot%restraint
    write(mdout, '(a)')
  end if
#endif

!! 6. Print stripped coordinates if asked for.
  if(hybridsolvent_remd_traj_name .ne. ' ') then
    hybridsolvent_remd_crd_cnt = 3 * hybridsolvent_remd_atm_cnt
    if (master) then
      ! Need to decrement mdloop because REMD exchanges happen on first step until
      ! total_steps - nstlim steps.
      my_mdloop = mdloop - 1
      local_total_nstep = (my_mdloop * nstlim)
      if(ioutfm .eq. 0) then
        if(my_mdloop .eq. 0) then
             call amopen(hybridsolvent_remd_traj, hybridsolvent_remd_traj_name, owrite, 'F', 'W')
             write(hybridsolvent_remd_traj, '(a)') 'Hybrid Solvent REMD stripped trajectory'
             write(hybridsolvent_remd_traj, '(a,3(1x,i8),1x,f8.3)') 'REMD ', repnum, &
                 my_mdloop, local_total_nstep, temp0
             call corpac(hybridsolvent_remd_crd_cnt, hybridsolvent_remd_crd, 1, hybridsolvent_remd_traj)
        else
             write(hybridsolvent_remd_traj, '(a,3(1x,i8),1x,f8.3)') 'REMD ', repnum, &
                 my_mdloop, local_total_nstep, temp0
             call corpac(hybridsolvent_remd_crd_cnt, hybridsolvent_remd_crd, 1, hybridsolvent_remd_traj)
        end if
      else ! ioutfm .eq. 1
#ifdef BINTRAJ
        !! write hybrid traj in netcdf format
        if(my_mdloop .eq. 0) then
          call open_hybridsolvent_remd_traj_binary_file(hybridsolvent_remd_atm_cnt)
          ! Need to pass remd arguments to avoid cyclic dependencies
          call hybridsolvent_remd_traj_setup_remd_indices(remd_method, remd_dimension, &
                                                    remd_types, group_num)
          ! Need to pass remd arguments to avoid cyclic dependencies
          call write_hybridsolvent_remd_traj_binary_crds(hybridsolvent_remd_crd, hybridsolvent_remd_atm_cnt, &
                                 remd_method, replica_indexes, remd_dimension, &
                                 remd_repidx, remd_crdidx)
        else
          ! Need to pass remd arguments to avoid cyclic dependencies
          call write_hybridsolvent_remd_traj_binary_crds(hybridsolvent_remd_crd, hybridsolvent_remd_atm_cnt, &
                                 remd_method, replica_indexes, remd_dimension, &
                                 remd_repidx, remd_crdidx)
        end if
#endif
      end if
    end if
  end if

  ! Do the exchange.
  if(rremd_type .le. 0) then
         call temperature_exchange(atm_cnt, vel, my_gb_pot_ene_tot%total, t_dim, &
                                num_replicas, actual_temperature, print_exch_data, mdloop)
  else
  ! If exchanging with reservoir, make sure the reservoir
  ! has coordinates and velocities for the complete (not stripped) system and
  ! has energies for the stripped (not complete) system.
         call reservoir_exchange(atm_cnt, crd, vel, my_gb_pot_ene_tot%total, t_dim, &
                                num_replicas, actual_temperature, print_exch_data, mdloop)
  end if

#ifdef VERBOSE_REMD
  if(master) then
    write(mdout,'(20("="),a,20("="))') "END HYBRID SOLVENT REMD EXCHANGE"
  end if
#endif

  return

end subroutine hybridsolvent_exchange

!*******************************************************************************
!
! Subroutine: setup_remd_randgen
!
! Description:
!
! Initializes the REMD random number stream. This *should* go in remd.fpp
! (subroutine remd_setup), but the Intel compilers fail when it's taken from a
! different module. I think this is caused because the size of random_state is
! not a multiple of its largest member (causes a warning), so ifort confuses the
! types. This seems to me to be a compiler bug, but oh well...
!
!*******************************************************************************

subroutine setup_remd_randgen

  use mdin_ctrl_dat_mod, only : ig
  use reservoir_mod,     only : reservoir_iseed
  use remd_mod,     only : remd_types

  implicit none

  integer rand_seed
  integer my_ig
  integer orepnum ! other replica number

  integer, dimension(mpi_status_size) :: stat_array ! needed for mpi_recv

  ! Only master nodes do exchanges, so seed the generator to protect against
  ! calling amrand on non-masters, but don't go through the logic below since
  ! the necessary communicator is not set up

  ! Set random # generator based on the reservoir_iseed from reservoir_name file.

  ! DAN ROE's comments from Sander: This reservoir_iseed can be used to do independent runs with different sequence of structures
  ! from the reservoir. Here we have no concept of the ig seed in md.in and we don't want to
  ! hardcode. So, the reservoir_iseed in reservoir related files is a convenient location.
  if (.not. master) then
    if(rremd_type .gt. 0) then
      call amrset_gen(remd_rand_gen, ig + reservoir_iseed + rand_seed )
      return
    else
      call amrset_gen(remd_rand_gen, ig + rand_seed )
      return
    end if
  end if

  ! sander and pmemd use different strategies when it comes to exchanging. While
  ! sander has only even replicas exchange 'up' and 'down', pmemd has every
  ! replica exchange 'up', and alternates which replica exchanges. Therefore,
  ! the random number streams between adjacent replicas here that are used for
  ! exchange attempt evaluatinos should be identical so that the random numbers
  ! used in REMD evaluations are synchronized between pmemd and sander. This way
  ! all results will be identical b/w the two programs and they will be easier
  ! to validate

  ! Set random # generator based on the reservoir_iseed from reservoir_name file.

  ! DAN ROE's comments from Sander: This reservoir_iseed can be used to do independent runs with different sequence of structures
  ! from the reservoir. Here we have no concept of the ig seed in md.in and we don't want to
  ! hardcode. So, the reservoir_iseed in reservoir related files is a convenient location.
  if(rremd_type .gt. 0 .and. remd_types(1).eq.1) then
    my_ig = ig + reservoir_iseed
  else
    my_ig = ig
  end if

  if (mod(repnum, 2) .eq. 1) then
    ! This is actually an 'even' replica when indexing from 0
    if (repnum .eq. 1) then
      orepnum = numgroups
    else
      orepnum = repnum - 1
    end if
    call mpi_send(my_ig, 1, mpi_integer, orepnum-1, 10, &
                  pmemd_master_comm, err_code_mpi)
  else
    if (repnum .eq. numgroups) then
      orepnum = 1
    else
      orepnum = repnum + 1
    end if
    call mpi_recv(my_ig, 1, mpi_integer, orepnum-1, 10, pmemd_master_comm, &
                  stat_array, err_code_mpi)
  end if
  rand_seed = mod(repnum-1, numgroups) / 2 * 2

  call amrset_gen(remd_rand_gen, my_ig + rand_seed)

  return

end subroutine setup_remd_randgen

!*******************************************************************************
!
! Subroutine: rxsgld_exchange
!
! Description: Performs the SGLD replica exchange attempt. It allow replicas with
!     different guiding effects and/or different temperatures to exchange.  It is
!     an extended version of temperature_exchange.
!
!*******************************************************************************
subroutine rxsgld_exchange(atm_cnt, my_atm_lst, amass, crd, vel, my_pot_ene_tot, t_dim, &
                                num_replicas, actual_temperature,    &
                                print_exch_data, mdloop)

  use mdin_ctrl_dat_mod, only : temp0
  use pmemd_lib_mod,     only : mexit
  use sgld_mod, only :  tsgset, sgfti,sgffi,psgldg,epotlf,epotllf, &
                     templf,rxsgld_scale

  implicit none

! Passed variables

  integer, intent(in)             :: atm_cnt
  integer                         :: my_atm_lst(*)
  integer, intent(in)             :: num_replicas
  integer, intent(in)             :: mdloop

  double precision, intent(in)    :: my_pot_ene_tot
  double precision, intent(inout) :: amass(atm_cnt)
  double precision, intent(inout) :: crd(3, atm_cnt)
  double precision, intent(inout) :: vel(3, atm_cnt)
  double precision, intent(in)    :: actual_temperature

  integer, intent(in)             :: t_dim ! which dimension T-exchanges are
  logical, intent(in)             :: print_exch_data

! Local variables

  type :: exchange_data ! facilitates gather for remlog
    sequence
    double precision :: temp0
    double precision :: sgft
    double precision :: sgff
    double precision :: psgldg
    double precision :: new_temp0
    double precision :: new_sgft
    double precision :: new_sgff
    double precision :: new_psgldg
    double precision :: real_temp
    double precision :: pot_ene_tot
    double precision :: epotlf
    double precision :: epotllf
    double precision :: templf
    double precision :: scaling
    double precision :: scalsg
    double precision :: exchange
  end type exchange_data

  integer, parameter :: SIZE_EXCHANGE_DATA = 16 ! for mpi_gather

  type (exchange_data) :: exch_buffer(num_replicas) ! gather buffer

  type (exchange_data) :: my_exch_data ! stores data for this replica
  type (exchange_data) :: exch_data_tbl(num_replicas) ! exch_data for all reps

   double precision  :: myscaling,o_scaling,myscalsg, o_scalsg
   double precision  :: o_sglf, sglf, o_sghf, sghf

   double precision  :: delta,dbeta,dbetamiu
   double precision  :: depot,depotlf,depotllf
  double precision  :: metrop
  double precision  :: random_value
  double precision  :: success_ratio

  integer           :: neighbor_rank
  integer           :: i,index_base
  integer           :: success_array(num_replicas)
  integer           :: success_buf(num_replicas)

  ! For exchanging replica positions in ladders upon successful exchanges
  integer           :: group_num_buffer(remd_dimension)
  integer           :: replica_indexes_buffer(remd_dimension)

  logical           :: success
  logical           :: i_do_exchg

  integer, dimension(mpi_status_size) :: stat_array ! needed for mpi_recv

! Explanation of local variables:
!
! delta: The calculated DELTA value used as part of the detailed balance calc.
!
! metrop: The actual value compared to a random_value to decide success
!
! success: Did this exchange attempt succeed?
!
! success_array: buffer for keeping track of which temperatures exchanged
!                successfully. Used in a mpi_reduce at the end
!
! success_buf: buffer for mpi_reduce call on success_array
!
! neighbor_rank: the remd_rank of the replica with whom we need to be doing ALL
!                of our communication. Trying to do everything with the
!                'partners' array is a bit too confusing.

! Set the variables that we know at the beginning

  my_exch_data%temp0       = temp0
  my_exch_data%sgft   = sgfti
  my_exch_data%sgff   = sgffi
  my_exch_data%psgldg   = psgldg
  my_exch_data%new_temp0       = temp0
  my_exch_data%new_sgft   = sgfti
  my_exch_data%new_sgff   = sgffi
  my_exch_data%new_psgldg   = psgldg
  my_exch_data%real_temp   = actual_temperature
  my_exch_data%pot_ene_tot = my_pot_ene_tot
  my_exch_data%epotlf     = epotlf
  my_exch_data%epotllf     = epotllf
  my_exch_data%templf     = templf
  my_exch_data%scaling     = -1.0d0
  my_exch_data%scalsg     = -1.0d0
  my_exch_data%exchange   = -1.0d0
  success_array(:)         = 0

  if (master) then
    my_exch_data%exchange   = replica_indexes(t_dim)

    call set_partners(t_dim, num_replicas, remd_random_partner)

    ! Call amrand_gen here to generate a random number. Sander does this
    ! to pick out structures from a reservoir which pmemd doesn't implement
    ! (yet). To keep the random #s synched, though, and to ensure proper
    ! regression testing, this must be done to keep the random #s synched
    ! b/w sander and pmemd

    call amrand_gen(remd_rand_gen, random_value)

#ifdef VERBOSE_REMD
    write(mdout, '(a,26("="),a,26("="))') '| ', 'RXSGLD EXCHANGE CALCULATION'
    write(mdout, '(a,i10,a,i1)') '| Exch= ', mdloop, ' RREMD= ', 0
    write(mdout, *) '| exdata= ', my_exch_data
#endif

    ! We alternate exchanges: temp up followed by temp down
    if (even_exchanges(t_dim)) then
      i_do_exchg = control_exchange(t_dim)
    else
      i_do_exchg = .not. control_exchange(t_dim)
    end if

    even_exchanges(t_dim) = .not. even_exchanges(t_dim)

    ! Find out what rank our neighbor is (-1 for indexing from 0), and
    ! If we are controlling the exchange, then our partner will be our
    ! 2nd partner (above). Otherwise it will be our 1st partner (below)

    if (i_do_exchg) then
      neighbor_rank = partners(2) - 1
    else
      neighbor_rank = partners(1) - 1
    end if

    ! Collect potential energies/temperatures
    call mpi_allgather(my_exch_data, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                exch_data_tbl, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                remd_comm, err_code_mpi)


#ifdef VERBOSE_REMD
    write(mdout, '(a,f8.2,f8.4,f8.4,2(a7,i2),a7,f10.2)') &
    '| Temp,sgft,sgff= ', &
      temp0,sgfti,sgffi, ' Indx= ', replica_indexes(t_dim), ' Rep#= ', &
      remd_rank+1, ' EPot= ', my_pot_ene_tot
    write(mdout, *) '| exdatn= ', exch_data_tbl(neighbor_rank+1)
#endif

    ! Calculate the exchange probability. Note that if delta is negative, then
    ! metrop will be > 1, and it will always succeed (so we don't need an extra
    ! conditional)

    if (i_do_exchg) then
#ifdef VERBOSE_REMD
      write(mdout, '(a,f6.2,2f6.4,2(a7,i2),a7,f10.2)') '| Partner Temp,sgft,sgff= ', &
        exch_data_tbl(neighbor_rank+1)%temp0,exch_data_tbl(neighbor_rank+1)%sgft,&
        exch_data_tbl(neighbor_rank+1)%sgff, ' Indx= ', neighbor_rank+1, &
        ' Rep#= ', neighbor_rank + 1, ' EPot= ', &
        exch_data_tbl(neighbor_rank+1)%pot_ene_tot
#endif


            ! Replica exchange self-guided Langevin dynamics
            !  Works also for temperature-based replica exchange
   ! delta = D(beta)D(Ep)+D(beta*miu)*(D(Eplf)-D(Epllf))
   !  D(beta)=1/kT(neighbor)-1/kT
   !  D(betamiu)=(sgff(neighbor)-psgldg(neighbor))/kT(neighbor)-(sgff-psgldg)/kT
   !  D(Ep)=Eppt(neighbor)-Eppt
   !  D(Eplf)=Epotlf(neighbor)-Epotlf
   !  D(Epllf)=Epotllf(neighbor)-Epotllf
   
    dbeta = 503.01d0/ exch_data_tbl(neighbor_rank+1)%temp0-503.01d0 / my_exch_data%temp0
    dbetamiu=503.01d0*(exch_data_tbl(neighbor_rank+1)%sgff-exch_data_tbl(neighbor_rank+1)%psgldg)/ exch_data_tbl(neighbor_rank+1)%temp0-  &
             503.01d0*(my_exch_data%sgff-my_exch_data%psgldg) / my_exch_data%temp0
    depot=exch_data_tbl(neighbor_rank+1)%pot_ene_tot- my_exch_data%pot_ene_tot 
    depotlf=exch_data_tbl(neighbor_rank+1)%epotlf - my_exch_data%epotlf 
    depotllf=exch_data_tbl(neighbor_rank+1)%epotllf - my_exch_data%epotllf 

    delta = dbeta * depot + dbetamiu*(depotlf-depotllf) 
                 

      metrop = exp(delta)

      call amrand_gen(remd_rand_gen, random_value)

      success = random_value .lt. metrop

      ! Let my neighbor know about our success

      call mpi_send(success, 1, mpi_logical, neighbor_rank, &
                    remd_tag, remd_comm, err_code_mpi)

      if (success) then
        success_array(replica_indexes(t_dim)) = 1
        myscaling = sqrt(exch_data_tbl(neighbor_rank+1)%temp0 / my_exch_data%temp0)
        myscalsg = sqrt(exch_data_tbl(neighbor_rank+1)%templf / my_exch_data%templf)
        my_exch_data%new_temp0 = exch_data_tbl(neighbor_rank+1)%temp0
        my_exch_data%new_sgft = exch_data_tbl(neighbor_rank+1)%sgft
        my_exch_data%new_sgff = exch_data_tbl(neighbor_rank+1)%sgff
        my_exch_data%new_psgldg = exch_data_tbl(neighbor_rank+1)%psgldg
        my_exch_data%scaling = myscaling
        my_exch_data%scalsg = myscalsg
        o_scaling = 1.0d0 / myscaling
        o_scalsg = 1.0d0 / myscalsg
      end if

#ifdef VERBOSE_REMD
      write(mdout,'(a8,E16.6,a8,E16.6,a12,f10.2,f10.2)') &
               "| Metrop= ",metrop," delta= ",delta," o_scaling,o_scalsg= ", &
               o_scaling,o_scalsg
#endif

    else

#ifdef VERBOSE_REMD
      ! Even if we don't control exchange, still print REMD info to mdout

      write(mdout, '(a,f6.2,2f6.4,2(a7,i2),a7,f10.2)') '| Partner Temp,sgft,sgff= ', &
        exch_data_tbl(neighbor_rank+1)%temp0,exch_data_tbl(neighbor_rank+1)%sgft, &
        exch_data_tbl(neighbor_rank+1)%sgff, ' Indx= ', neighbor_rank+1, &
        ' Rep#= ', neighbor_rank + 1, ' EPot= ', &
        exch_data_tbl(neighbor_rank+1)%pot_ene_tot
      write(mdout, '(a)') '| Not controlling exchange.'
#endif

      call amrand_gen(remd_rand_gen, random_value) ! to stay synchronized

      ! Get the message from the exchanging replica about our success
      call mpi_recv(success, 1, mpi_logical, neighbor_rank, &
                    remd_tag, remd_comm, stat_array, err_code_mpi)

      ! We scale velocities only if we succeed
      if (success) then
        myscaling = sqrt(exch_data_tbl(neighbor_rank+1)%temp0 / my_exch_data%temp0)
        myscalsg = sqrt(exch_data_tbl(neighbor_rank+1)%templf / my_exch_data%templf)
        my_exch_data%new_temp0 = exch_data_tbl(neighbor_rank+1)%temp0
        my_exch_data%new_sgft = exch_data_tbl(neighbor_rank+1)%sgft
        my_exch_data%new_sgff = exch_data_tbl(neighbor_rank+1)%sgff
        my_exch_data%new_psgldg = exch_data_tbl(neighbor_rank+1)%psgldg
        my_exch_data%scaling = myscaling
        my_exch_data%scalsg = myscalsg
        o_scaling = 1.0d0 / myscaling
        o_scalsg = 1.0d0 / myscalsg
      end if

    end if

#ifdef VERBOSE_REMD
    write(mdout,'(a8,E16.6,a12,f10.2,f10.2,a10,L1)') &
      '| Rand=   ', random_value, ' MyScaling,myscalsg= ', myscaling, &
      myscalsg,' Success= ', success

    write(mdout,'(a,24("="),a,24("="))') '| ', "END RXSGLD EXCHANGE CALCULATION"
#endif

  end if ! master

  ! We'll broadcast all of the my_exch_data here once. From that, we can tell
  ! if the exchange succeeded based on whether or not the temperatures have
  ! changed. Then make sure we update temp0 so the simulation reflects that.

  call mpi_bcast(my_exch_data, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                 0, pmemd_comm, err_code_mpi)
  if (.not. master ) &
                 success = my_exch_data%scaling.gt.TINY
             

#ifdef VERBOSE_REMD
  if (master) &
    write(mdout,*)'| rank, sucess: ', &
       my_exch_data%exchange,replica_indexes(t_dim),neighbor_rank+1
   if (master) &
     write(mdout, '(2(a,i4,1x,f6.2,f6.4,f6.4))') &
      '| RXSGLD: checking to see if TEMP0,SGFT,sgff has changed: ', &
      replica_indexes(t_dim),temp0,sgfti,sgffi, '->', &
     int(my_exch_data%exchange),my_exch_data%temp0,  &
      my_exch_data%sgft,my_exch_data%sgff
#endif

  if (success)then
   ! ---=== RXSGLD property SCALING ===---
     temp0 = my_exch_data%new_temp0
     tsgset=temp0
     sgfti = my_exch_data%new_sgft
     sgffi = my_exch_data%new_sgff
     psgldg = my_exch_data%new_psgldg
     call rxsgld_scale(atm_cnt,my_atm_lst, my_exch_data%scaling,my_exch_data%scalsg, amass, crd, vel)
  endif
  ! We want to keep track of our exchange successes, so gather them
  ! from all of the masters here, and increment them in exchange_successes
  ! We do this whether or not we print so we can keep track of exchange
  ! success ratios. For the index and exchange data, only gather those to
  ! the remd_master if we actually plan to print. This should be much
  ! more efficient for high exchange attempts.

  ! Then recollect the temperatures as some may have changed, and reset our
  ! partners.
  if (master) then

    if (print_exch_data) then
      call mpi_gather(my_exch_data, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                      exch_buffer, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                      0, remd_comm, err_code_mpi)
    end if

    call mpi_reduce(success_array, success_buf, num_replicas, mpi_integer, &
                    mpi_sum, 0, remd_comm, err_code_mpi)

    if (remd_master) &
      exchange_successes(t_dim,:) = exchange_successes(t_dim,:) + success_buf(:)

    ! Increment exchange counter, and swap our replica ranks and numbers with
    ! our neighbor if our attempt succeeded, since we effectively swapped places
    ! in this array. That's because we swap temperatures and not structures due
    ! to performance considerations.

    if (success) then
      call mpi_sendrecv(group_num, remd_dimension, mpi_integer, &
                        neighbor_rank, remd_tag, &
                        group_num_buffer, remd_dimension, mpi_integer, &
                        neighbor_rank, remd_tag, remd_comm, &
                        stat_array, err_code_mpi)
      call mpi_sendrecv(replica_indexes, remd_dimension, mpi_integer, &
                        neighbor_rank, remd_tag, &
                        replica_indexes_buffer, remd_dimension, mpi_integer, &
                        neighbor_rank, remd_tag, remd_comm, &
                        stat_array, err_code_mpi)
      call mpi_sendrecv_replace(remd_repidx, 1, mpi_integer, &
                                neighbor_rank, remd_tag, &
                                neighbor_rank, remd_tag, &
                                remd_comm, stat_array, err_code_mpi)
      group_num = group_num_buffer
      replica_indexes = replica_indexes_buffer
    end if
  end if

  ! If we're doing NMROPT, then we need to make sure we re-read in the TEMP0,
  ! since this subroutine overwrites TEMP0 to the original since you can modify
  ! temperature with NMROPT settings.

  remd_modwt = .true.

  ! Write the data to the remlog

  ! remd_method == -1 means we're doing multi-D REMD.
  if (remd_method .ne. -1) then
    if (print_exch_data .and. remd_master) then
      write(remlog, '(a,i8)') '# exchange ', mdloop
      do i = 1, num_replicas
        success_ratio = dble(exchange_successes(t_dim, index_list(i))) / &
                        dble(mdloop) * dble(remd_dimension) * 2
        write(remlog, '(i4,i4,2f8.4,2f8.2,e14.6,f8.4)') &
              i, int(exch_buffer(i)%exchange),&
              exch_buffer(i)%scaling, &
              exch_buffer(i)%scalsg, &
              exch_buffer(i)%real_temp, &
              exch_buffer(i)%templf, &
              exch_buffer(i)%pot_ene_tot, &
              success_ratio
      end do
      flush(mdout)

    end if
  else if (remd_master) then
    do i = 1, num_replicas
      success_ratio = dble(exchange_successes(t_dim, index_list(i))) / &
                      dble(mdloop) * dble(remd_dimension) * 2
      multid_print_data(i)%scaling       = exch_buffer(i)%scaling
      multid_print_data(i)%real_temp     = exch_buffer(i)%real_temp
      multid_print_data(i)%pot_ene_tot   = exch_buffer(i)%pot_ene_tot
      multid_print_data(i)%temp0         = exch_buffer(i)%temp0
      multid_print_data(i)%new_temp0     = exch_buffer(i)%new_temp0
      multid_print_data(i)%group_num     = group_num(t_dim)
      multid_print_data(i)%success_ratio = success_ratio
      multid_print_data(i)%num_rep       = num_replicas
    end do
  end if ! remd_method .ne. -1
#ifdef VERBOSE_REMD
  if (master) &
    write(mdout, '(2(a,i8))') &
      '| RXSGLD: checking  ', &
      replica_indexes(t_dim), '->', &
     int(my_exch_data%exchange)
#endif

  return

end subroutine rxsgld_exchange

#endif /* MPI */

end module remd_exchg_mod
