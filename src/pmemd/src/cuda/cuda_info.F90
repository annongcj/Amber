#include "copyright_f90.i"

!----------------------------------------------------------------------------------------------
! gpu_write_cuda_info: Routines to print info regarding CUDA based GPU Calculations.
!                      (By Ross C. Walker)
!
! Arguments:
!   iout_unit:           unit to write output to (master_setup supplies the value mdout)
!   mytaskid:            MPI task ID
!   numtasks:            the number of MPI tasks
!   using_pme_potential: integer that functions as a logical (it IS a logical up in the
!                        Fortran code, unlike the other flags fed to this subroutine)
!   iamd:                flag to denote accelerated molecular dynamics (NOT Advanced Micro
!                        Devices).
!   igamd:               flag to activate Gaussian accelerated molecular dynamics
!                        (akin to meta-dynamics).
!   icnstph:             flag to activate constant pH molecular dynamics (with either implicit
!                        or explicit solvent).
!   icnste:              flag to activate constant Redox Potential molecular dynamics (with
!                        either implicit or explicit solvent).
!   icfe:                flag to activate alchemical free energy (either TI, MBAR, BAR or FEP).
!----------------------------------------------------------------------------------------------
#ifdef CUDA
subroutine gpu_write_cuda_info(iout_unit, mytaskid, numtasks, using_pme_potential, &
                               iamd, igamd, icnstph, icnste, icfe)

  implicit none

  ! Formal arguments
  integer, intent(in) :: iout_unit, mytaskid, numtasks, using_pme_potential, &
                         iamd, igamd, icnstph, icnste, icfe

  ! Local variables
  integer :: gpu_dev_count, gpu_dev_id, gpu_dev_mem
  integer :: gpu_num_multi, gpu_num_proc, name_len, i
  double precision :: gpu_core_freq
  character(len=80) :: gpu_name
  logical :: b_single_node, b_p2p_enabled, b_nccl_enabled
  character(len=512) :: char_tmp_512
  integer :: env_var_status

  write(iout_unit,'(a)') '|--------------------- INFORMATION ----------------------'
  write(iout_unit,'(a)') '| GPU (CUDA) Version of PMEMD in use: NVIDIA GPU IN USE.'
  write(iout_unit,'(a)') '|                    Version 18.0.0'
  write(iout_unit,'(a)') '| '
  write(iout_unit,'(a)') '|                      03/25/2018'
  write(iout_unit,'(a)') '| '
  write(iout_unit,'(a)') '| Implementation by:'
  write(iout_unit,'(a)') '|                    Ross C. Walker     (SDSC)'
  write(iout_unit,'(a)') '|                    Scott Le Grand     (nVIDIA)'
  write(iout_unit,'(a)') '| '
  write(iout_unit,'(a)') '| Version 18 performance extensions by:'
  write(iout_unit,'(a)') '|                    David Cerutti     (Rutgers)'
  write(iout_unit,'(a)') '| '
  write(iout_unit,'(a)') '| Precision model in use:'
#ifdef use_SPFP
  write(iout_unit,'(a)') '|      [SPFP] - Single Precision Forces, 64-bit Fixed Point'
  write(iout_unit,'(a)') '|               Accumulation. (Default)'
#elif use_DPFP
  write(iout_unit,'(a)') '|      [DPFP] - Double Precision Forces, 64-bit Fixed point'
  write(iout_unit,'(a)') '|               Accumulation.'
#else
  write(iout_unit,'(a)') '|      ERROR ERROR ERROR - Unknown GPU Precision Model'
  call mexit(iout_unit, 1)
#endif
  write(iout_unit,'(a)') '| '
  write(iout_unit,'(a)') '|--------------------------------------------------------'
  write(iout_unit,'(a)') ' '

  !Write citation information for the GPU code to mdout.
  call write_cuda_citation(iout_unit, using_pme_potential, iamd, igamd, icnstph, icnste, icfe)

  ! Get info about GPU(s) in use.
  write(iout_unit,'(a)') '|------------------- GPU DEVICE INFO --------------------'
#ifdef MPI
  do i = 0, numtasks - 1
    if (i == 0) then
#endif
      call gpu_get_device_info(gpu_dev_count, gpu_dev_id, gpu_dev_mem, gpu_num_proc, &
                               gpu_core_freq, gpu_name, name_len)
#ifdef MPI
    else
      call gpu_get_slave_device_info(i, gpu_dev_count, gpu_dev_id, gpu_dev_mem, &
                                     gpu_num_proc, gpu_core_freq, gpu_name, name_len)
    end if
#endif
    write(iout_unit,'(a)')      '|'
#ifdef MPI
    write(iout_unit,'(a,i6)')   '|                         Task ID: ',i
#endif
    call get_environment_variable("CUDA_VISIBLE_DEVICES", char_tmp_512, env_var_status)
    if (env_var_status == 0) then
      write(iout_unit,'(a)')    '|            CUDA_VISIBLE_DEVICES: not set'
    else
      write(iout_unit,'(a,a)') '|            CUDA_VISIBLE_DEVICES: ', trim(char_tmp_512)
    end if
    write(iout_unit,'(a,i6)')   '|   CUDA Capable Devices Detected: ', gpu_dev_count
    write(iout_unit,'(a,i6)')   '|           CUDA Device ID in use: ', gpu_dev_id
    write(iout_unit,'(a,a)')    '|                CUDA Device Name: ', gpu_name(1:name_len)
    write(iout_unit,'(a,i6,a)') '|     CUDA Device Global Mem Size: ', gpu_dev_mem, ' MB'
    write(iout_unit,'(a,i6)')   '| CUDA Device Num Multiprocessors: ', gpu_num_proc
    write(iout_unit,'(a,f6.2,a)') '|           CUDA Device Core Freq: ', gpu_core_freq, ' GHz'
    write(iout_unit,'(a)')      '|'
#ifdef MPI
  end do
#endif
  write(iout_unit,'(a)') '|--------------------------------------------------------'
  write(iout_unit,'(a)')   ' '

#ifdef MPI
  write(iout_unit,'(a)') '|---------------- GPU PEER TO PEER INFO -----------------'
  write(iout_unit,'(a)') '|'
  call gpu_get_p2p_info(b_single_node, b_p2p_enabled, b_nccl_enabled)
  if (.not. b_single_node) then
    write(iout_unit,'(a)') '|   Peer to Peer support: DISABLED'
    write(iout_unit,'(a)') '|'
    write(iout_unit,'(a)') '|              (MPI tasks span multiple nodes)'
  else if (.not. b_p2p_enabled) then
    write(iout_unit,'(a)') '|   Peer to Peer support: DISABLED'
    write(iout_unit,'(a)') '|'
    write(iout_unit,'(a)') '|         (Selected GPUs cannot communicate over P2P)'
  else  !both b_single_node and b_p2p_enabled are true.
    write(iout_unit,'(a)') '|   Peer to Peer support: ENABLED'
  end if
  write(iout_unit,'(a)') '|'
  if (b_nccl_enabled) then
     write(iout_unit,'(a)') '|   NCCL support: ENABLED'
  else
     write(iout_unit,'(a)') '|   NCCL support: DISABLED'
  end if
  write(iout_unit,'(a)') '|'
  write(iout_unit,'(a)') '|--------------------------------------------------------'
  write(iout_unit,'(a)') ' '
#endif

  return

end subroutine gpu_write_cuda_info

!----------------------------------------------------------------------------------------------
! write_cuda_citation: write the pmemd.CUDA citations so that users can give us press.
!
! Arguments:
!   See gpu_write_cuda_info() above for descriptions.
!----------------------------------------------------------------------------------------------
subroutine write_cuda_citation(iout_unit, using_pme_potential, iamd, igamd, icnstph, icnste, icfe)

  implicit none

  ! Formal arguments
  integer, intent(in) :: iout_unit, using_pme_potential, iamd, igamd, icnstph, icnste, icfe

  write(iout_unit,'(a)') '|----------------- CITATION INFORMATION -----------------'
  write(iout_unit,'(a)') '|'
  write(iout_unit,'(a)') '|    When publishing work that utilized the CUDA version'
  write(iout_unit,'(a)') '|    of AMBER, please cite the following in addition to'
  write(iout_unit,'(a)') '|    the regular AMBER citations:'
  write(iout_unit,'(a)') '|'
  if (using_pme_potential .ne. 0) then
    write(iout_unit,'(a)') '|  - Romelia Salomon-Ferrer; Andreas W. Goetz; Duncan'
    write(iout_unit,'(a)') '|    Poole; Scott Le Grand; Ross C. Walker "Routine'
    write(iout_unit,'(a)') '|    microsecond molecular dynamics simulations with'
    write(iout_unit,'(a)') '|    AMBER - Part II: Particle Mesh Ewald", J. Chem.'
    write(iout_unit,'(a)') '|    Theory Comput., 2013, 9 (9), pp3878-3888,'
    write(iout_unit,'(a)') '|    DOI: 10.1021/ct400314y.'
    write(iout_unit,'(a)') '|'
  endif
  write(iout_unit,'(a)') '|  - Andreas W. Goetz; Mark J. Williamson; Dong Xu;'
  write(iout_unit,'(a)') '|    Duncan Poole; Scott Le Grand; Ross C. Walker'
  write(iout_unit,'(a)') '|    "Routine microsecond molecular dynamics simulations'
  write(iout_unit,'(a)') '|    with AMBER - Part I: Generalized Born", J. Chem.'
  write(iout_unit,'(a)') '|    Theory Comput., 2012, 8 (5), pp1542-1555.'
#ifdef use_SPFP
  write(iout_unit,'(a)') '|'
  write(iout_unit,'(a)') '|  - Scott Le Grand; Andreas W. Goetz; Ross C. Walker'
  write(iout_unit,'(a)') '|    "SPFP: Speed without compromise - a mixed precision'
  write(iout_unit,'(a)') '|    model for GPU accelerated molecular dynamics'
  write(iout_unit,'(a)') '|    simulations.", Comp. Phys. Comm., 2013, 184'
  write(iout_unit,'(a)') '|    pp374-380, DOI: 10.1016/j.cpc.2012.09.022'
#endif
  if (iamd .ne. 0) then
    write(iout_unit,'(a)') '|'
    write(iout_unit,'(a)') '|    When publishing work that utilized the CUDA version'
    write(iout_unit,'(a)') '|    of AMD, please cite the following in addition to'
    write(iout_unit,'(a)') '|    the regular AMBER citations:'
    write(iout_unit,'(a)') '|'
    write(iout_unit,'(a)') '|  - Levi C. T. Pierce; Romelia Salomon-Ferrer; '
    write(iout_unit,'(a)') '|    Cesar Augusto F de Oliveira; J. Andrew McCammon'
    write(iout_unit,'(a)') '|    and Ross C. Walker "Routine access to milli-second '
    write(iout_unit,'(a)') '|    time scales with accelerated molecular dynamics".'
    write(iout_unit,'(a)') '|    J. Chem. Theory Comput., 2012, 8(9), pp2997-3002.'
    write(iout_unit,'(a)') '|    DOI: 10.1021/ct300284c.'
    write(iout_unit,'(a)') '|'
  endif
  if (igamd .ne. 0) then
    write(iout_unit,'(a)') '|'
    write(iout_unit,'(a)') '|    When publishing work that utilized the CUDA version'
    write(iout_unit,'(a)') '|    of Gaussian Accelerated Molecular Dynamics(GaMD), '
    write(iout_unit,'(a)') '|    please cite the following in addition to'
    write(iout_unit,'(a)') '|    the regular AMBER GPU citations:'
    write(iout_unit,'(a)') '|'
    write(iout_unit,'(a)') '|  - Yinglong Miao; Victoria A. Feher; J. Andrew McCammon'
    write(iout_unit,'(a)') '|    "Gaussian Accelerated Molecular Dynamics: Unconstrained '
    write(iout_unit,'(a)') '|    Enhanced Sampling and Free Energy Calculation".'
    write(iout_unit,'(a)') '|    J. Chem. Theory Comput., 2015, 11(8):3584-95.'
    write(iout_unit,'(a)') '|    DOI: 10.1021/acs.jctc.5b00436.'
    write(iout_unit,'(a)') '|'
  endif
  if (icnstph .ne. 0) then
    write(iout_unit,'(a)') '|'
    write(iout_unit,'(a)') '|    When publishing work that utilized the CUDA version'
    write(iout_unit,'(a)') '|    of Constant pH MD please cite the following in'
    write(iout_unit,'(a)') '|    addition to the regular AMBER GPU citations:'
    write(iout_unit,'(a)') '|'
    write(iout_unit,'(a)') '|  - Daniel J. Mermelstein; J. Andrew McCammon; Ross C. Walker'
    write(iout_unit,'(a)') '|    "pH dependent conformational dynamics of Beta-secretase 1:'
    write(iout_unit,'(a)') '|    a molecular dynamics study".'
    write(iout_unit,'(a)') '|    J. Chem. Theory Comput., 2018, in review.'
    write(iout_unit,'(a)') '|'
  endif
  if (icnste .ne. 0) then
    write(iout_unit,'(a)') '|'
    write(iout_unit,'(a)') '|    When publishing work that utilized the CUDA version'
    write(iout_unit,'(a)') '|    of Constant Redox Potential MD please cite the'
    write(iout_unit,'(a)') '|    following in addition to the regular AMBER GPU citations:'
    write(iout_unit,'(a)') '|'
    write(iout_unit,'(a)') '|  - Vinicius Wilian D. Cruzeiro; Marcos S. Amaral; Adrian E. Roitberg'
    write(iout_unit,'(a)') '|    "Redox potential replica exchange molecular dynamics at constant'
    write(iout_unit,'(a)') '|    pH in AMBER: Implementation and validation".'
    write(iout_unit,'(a)') '|    J. Chem. Phys., 2018, 149, 072338.'
    write(iout_unit,'(a)') '|    DOI: 10.1063/1.5027379'
    write(iout_unit,'(a)') '|'
  endif
  if (icfe .ne. 0) then
    write(iout_unit,'(a)') '|'
    write(iout_unit,'(a)') '|    When publishing work that utilized the CUDA version'
    write(iout_unit,'(a)') '|    of TI, BAR, MBAR or FEP please cite the following '
    write(iout_unit,'(a)') '|    publications in addition to the regular AMBER '
    write(iout_unit,'(a)') '|    GPU citations:'
    write(iout_unit,'(a)') '|'
    write(iout_unit,'(a)') '|  - Daniel J. Mermelstein; Charles Lin; Gard Nelson; '
    write(iout_unit,'(a)') '|    Rachael Kretsch; J. Andrew McCammon; Ross C. Walker'
    write(iout_unit,'(a)') '|    "Fast and Flexible GPU Accelerated Binding '
    write(iout_unit,'(a)') '|    Free Energy Calculations within the AMBER Molecular'
    write(iout_unit,'(a)') '|    Dynamics Package" J. Comp. Chem., 2018,'
    write(iout_unit,'(a)') '|    DOI: 10.1002/jcc.25187'
    write(iout_unit,'(a)') '|'
    write(iout_unit,'(a)') '|  - Tai-Sung Lee; Yuan Hu; Brad Sherborne; Zhuyan Guo;'
    write(iout_unit,'(a)') '|    Darrin M. York'
    write(iout_unit,'(a)') '|    "Toward Fast and Accurate Binding Affinity Prediction with'
    write(iout_unit,'(a)') '|    pmemdGTI: An Efficient Implementation of GPU-Accelerated'
    write(iout_unit,'(a)') '|    Thermodynamic Integration"'
    write(iout_unit,'(a)') '|    J. Chem. Theory Comput., 2017, 13 (7), 3077'
    write(iout_unit,'(a)') '|'
  endif
  write(iout_unit,'(a)') '|'
  write(iout_unit,'(a)') '|--------------------------------------------------------'
  write(iout_unit,'(a)') ' '

  return

end subroutine write_cuda_citation

!----------------------------------------------------------------------------------------------
! gpu_write_memory_info: write information about GPU memory usage to the mdout file
!
! Arguments:
!   iout_unit:    theunit to write to (given as mdout in the pmemd main program)
!----------------------------------------------------------------------------------------------
subroutine gpu_write_memory_info(iout_unit)

  implicit none

  ! Formal arguments:
  integer, intent(in)   :: iout_unit !Unit to write info to.

  ! Local variables:
  integer :: gpumemory, cpumemory

  call gpu_get_memory_info(gpumemory, cpumemory)

  write(iout_unit,'(a)') '| GPU memory information (estimate):'
  write(iout_unit,'(a,i9)') '| KB of GPU memory in use: ', gpumemory
  write(iout_unit,'(a,i9/)') '| KB of CPU memory in use: ', cpumemory

  return

end subroutine gpu_write_memory_info

#endif
