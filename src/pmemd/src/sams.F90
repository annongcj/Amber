!*******************************************************************************
!
! Module: sams_mod
!
! SAMS implementaion by 
! Taisung Lee and Darrin York
!
! Ref: Tan, Z. (2017) Optimally adjusted mixture sampling and locally weighted histogram analysis, 
! Journal of Computational and Graphical Statistics, 26, 54-65.
!
!*******************************************************************************

module sams_mod

  use gbl_constants_mod
  use energy_records_mod
  use state_info_mod
  use mdin_ctrl_dat_mod
  use ti_mod
    
  use file_io_mod
  use file_io_dat_mod
  use inpcrd_dat_mod
  use pbc_mod
  use mdin_ewald_dat_mod, only : skinnb
  use reservoir_mod

  use random_mod,      only : random_state, amrand_gen, amrset_gen, amrand
  implicit none  

  !! To be put into namelist
  integer::sams_direction=0
  double precision::sams_delta=0.1 
  
  ! Random number generator
  !type(random_state), save, private :: sams_rand_gen
  
  double precision, save :: onset_time
  double precision, dimension(:), allocatable, save:: logP, logPi, logPZ, &
      logZ_drive, logZ_report, logWeight
  double precision, dimension(:), allocatable :: jumpProb
  double precision, dimension(:), allocatable, save:: targetProb;
  double precision, dimension(:), allocatable, save:: dvdl_drive, dvdl_report
  double precision, dimension(:), allocatable, save:: dvdl_square

  integer, dimension(:, :), allocatable, save:: transition_counter
  double precision, dimension(:, :), allocatable, save:: transition_matrix
  double precision, dimension(:), allocatable, save:: eigval
  ! These are scratch/dummy arrays for LAPACK.
  double precision, dimension(:), allocatable, save:: eigval_imag
  double precision, dimension(:,:), allocatable, save:: eigvec
  double precision, dimension(:), allocatable, save:: eigscratch
  
  integer, dimension(:), allocatable, save::lambda_counter, sample_counter
  integer, dimension(:), allocatable, save::resident_time
  integer, save:: currentIteration=0
  integer, save:: currentLambdaIndex=0
  integer, save:: burnin_stage=0 ! 0: initial; 1:optimal
  integer, save:: scan_stage=0
  integer, save:: n_considered_states
  integer, save::local_atm_cnt
  
  integer, save:: sams_restart_step=0
  double precision, dimension(:,:), allocatable, save:: q_ave;
  logical, save:: sams_onstep, sams_do_average=.false.
  
  double precision, save:: current_dvdl=0.0, dvdl_ave=0.0
  
  type(random_state), save :: sams_rand_gen
  double precision  :: my_random_value
  
  ! formatters
  character(len=160), parameter :: form_iter = &
    "('#Step: ',I8, '  Iter#:',I8,'  Sample#:',I8 &
    ,'  Lambda Index:', I4, ' Lambda:', F8.5 &
            ,'  Stage:', I3)"
  character(len=160), parameter :: form_estimate = &
    "('#Estimate deltaF:',F12.6, ' Smoothness:', &
    F12.6,'  dvdl:', F12.6, '  dvdl_ave:', F12.6 )"  
  character(len=160), parameter :: form_counter = "(10I10)"
  character(len=160), parameter :: form_logZ = "(10F10.4)"
  character(len=160), parameter :: form_dvdl = "(10F10.4)"
  character(len=160), parameter :: form_mbar = "(10F12.4)"
  character(len=160), parameter :: form_prob = "(20F7.4)"
  
contains

function smoothness(array, n) result(sm)

  implicit none

  integer, intent(in)::n
  integer, intent(in), dimension(n)::array
  double precision::sm

  double precision::t0, tt

  t0=sum(array)  
  tt=0.0
  sm=0.0
  if(t0.gt.0) then
    t0=t0/(n*1.0)
    tt=sum((array(1:n)-t0)**2) 
    sm=sqrt(tt/(n*1.0))/t0
  endif

end 

subroutine sams_setup(irest, error, atm_cnt)
 
  implicit none

  integer, intent(in)::irest, atm_cnt
  integer, intent(out)::error

  character(10) :: char_date
  character(10) :: char_time
  character(5)  :: char_zone

  integer       :: dt(8)
  integer       :: year, month, day, hour, minute, sec, msec
  integer       :: feb_days
  integer       :: i,j
 
  logical       :: read_init
  
  double precision :: box(3), factor(3)
  
  !make sure parameters are in the right range
  if (bar_states < 0) then
     error=-1
     return
  endif 
  
  sams_alpha=max(0.1, sams_alpha)  
  sams_beta=max(0.00001, sams_beta)  
  sams_gamma=max(0.00001, sams_gamma)
  sams_scan=max(sams_scan, 1)
  
  i=min(sams_lambda_start, sams_lambda_end)
  j=max(sams_lambda_start, sams_lambda_end)
  sams_lambda_start=max(1, i)
  sams_lambda_end=min(bar_states, j)
  n_considered_states=(sams_lambda_end-sams_lambda_start+1)
  
  sams_jump_range=max(1, sams_jump_range)
  sams_update_freq=max(1, sams_update_freq)  
  

  
  call date_and_time(char_date, char_time, char_zone, dt)
  call amopen(sams_log, sams_log_name, 'U', 'F', 'W')
  write(sams_log,'("SAMS output: ", A10, A10, A5)')char_date, &
    char_time, char_zone
  write(sams_log,'("SAMS output verbose:", I10)') sams_verbose 
  write(sams_log,*)"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
  write(sams_log,'("SAMS Type = ", I4)') sams_type
  write(sams_log,'("SAMS Scan = ", I4)') sams_scan
  write(sams_log,'("alpha =", F10.6)') sams_alpha  
  write(sams_log,'("beta =", F10.6)') sams_beta  
  write(sams_log,'("gamma =", F10.6)') sams_gamma
  write(sams_log,'("Lambda range =", 2I8)')  &
    sams_lambda_start, sams_lambda_end    
  write(sams_log,'("jump range =", I10)') sams_jump_range
  write(sams_log,'("update Freq =", I10)') sams_update_freq  
  write(sams_log,'("SAMS starts at:", I10, " time step")') sams_init_start 
  write(sams_log,'("Statisics starts at:", I10, " time step")') sams_stat_start  
  write(sams_log,'("SAMS optimal estimator ends at:", I10, " SAMS iteration")') &
    sams_opt_end 
  write(sams_log,'("SAMS Variance control:",I4)') sams_variance 
  
  if (sams_opt_end .lt.0) sams_opt_end=huge(sams_opt_end)
  onset_time=0.0  
   
  allocate(logP(bar_states))
  allocate(logPZ(bar_states))
  allocate(logPi(bar_states))
  allocate(logZ_drive(bar_states))
  allocate(logZ_report(bar_states))  
  allocate(logWeight(bar_states))
  allocate(targetProb(bar_states))  
  allocate(resident_time(bar_states))  
  allocate(dvdl_drive(bar_states))
  allocate(dvdl_report(bar_states))  
  allocate(dvdl_square(bar_states))  

  allocate(jumpProb(bar_states))  
  allocate(lambda_counter(bar_states))  
  allocate(sample_counter(bar_states))    

  allocate(transition_counter(bar_states, bar_states))
  allocate(transition_matrix(bar_states, bar_states))
  allocate(eigval(bar_states))
  allocate(eigval_imag(bar_states))
  allocate(eigvec(bar_states, bar_states))
  ! Note that this same size MUST be passed to LAPACK calls.
  allocate(eigscratch(4*bar_states))

  !!TBT
  !!allocate(q_ave(bar_states, sams_update_freq))
  allocate(q_ave(bar_states, max(1000, sams_update_freq)))
  q_ave=0.0
 
  logP(:)=0.0
  logPZ(:)=0.0
  logPi(:)=0.0  
  targetProb(:)=0.0
  resident_time(:)=0
  logWeight(:)=0.0
  targetProb(sams_lambda_start:sams_lambda_end)=1.0/(n_considered_states)
  logWeight(sams_lambda_start:sams_lambda_end) = &
    log(targetProb(sams_lambda_start:sams_lambda_end))

  read_init=.false.;
  if (sams_restart .eq. 1 .and. irest .eq. 0) then
      write(*,*) "SAMS can only restart when irest=1. Start normal SAMS."
      sams_restart=0
  endif
  
  if (sams_type.le.2) then  !! SAMS run
    
    if (sams_restart .ne. 0) then
      call sams_read_init(read_init)
      if (.not. read_init) then
        write(*,*) "SAMS init file read failed.  Start normal SAMS." 
      else 
        call sams_output(sams_restart_step)
      endif
    end if
  
    if (.not. read_init) then
      currentLambdaIndex=sams_lambda_start
      currentIteration=0
      burnin_stage=0 ! initial  
      if (sams_type .eq. 2) then 
        scan_stage=0
      else 
        scan_stage=1
      endif
      logZ_drive(:)=0.0
      logZ_report(:)=0.0    
      dvdl_drive(:)=0.0
      dvdl_report(:)=0.0   
      dvdl_square(:)=0.0
      lambda_counter(:)=0
      lambda_counter(currentLambdaIndex)=1
      sample_counter(:)=0
      transition_counter(:,:)=0
    end if
    
  else !! NFE run
    
      local_atm_cnt=atm_cnt
      call load_reservoir_files(1)
      
      CALL SYSTEM_CLOCK(i)
      call amrset_gen(sams_rand_gen, ig+reservoir_iseed+i)
      call amrand_gen(sams_rand_gen, my_random_value) 
      rremd_idx(1) = min(reservoir_size(1),  &
        int(my_random_value * reservoir_size(1)) + 1)
        !!TBT
      !!rremd_idx(1) =4147
      call load_reservoir_structure(atm_crd, atm_vel, atm_cnt,  &
        box(1), box(2), box(3), 1)
      atm_last_vel(:,:)=atm_vel(:,:)
      
      factor(:)=box(:)/pbc_box(:)
      
      call pressure_scale_pbc_data(factor, vdw_cutoff + skinnb, 0)
      
      write(sams_log,'("NFE: Init Load Reservoir # ",I6)') rremd_idx(1) 
      
      if (sams_direction .eq. 0 ) then
        currentLambdaIndex=sams_lambda_start
      else
        currentLambdaIndex=sams_lambda_end
      endif   

      lambda_counter(:)=0
      lambda_counter(currentLambdaIndex)=1
      currentIteration=0
      scan_stage=0
      dvdl_drive(:)=0.0
  endif

  error=0
  close (sams_log)
  
end subroutine

subroutine sams_update_weight()

  implicit none

  if (sams_variance.eq.2) then
    logWeight(:)=-logZ_drive(:)
  else
    logWeight(:)=log(targetProb(:))-logZ_drive(:)
  endif   
 
end subroutine

subroutine sams_update_Z_estimate()
 
  implicit none
  double precision::pi_star, gamma, gamma_report
  double precision:: t0, tt
  integer::i,j
  logical, save::updateTarget=.false.

  j=currentIteration/sams_update_freq
  if(sams_variance.ne.0 .and. j .le. sams_opt_end &
    .and. burnin_stage.eq.1  .and. .not. updateTarget) then
    updateTarget=.true.
    do i=sams_lambda_start, sams_lambda_end
      updateTarget= updateTarget .and. (sample_counter(i)>10) 
    end do
  endif
  
  if(updateTarget .and. mod(j,10).eq.0 ) then
    tt=0
    do i=sams_lambda_start, sams_lambda_end
        tt=tt+(dvdl_square(i)-dvdl_report(i)*dvdl_report(i))/sample_counter(i)
    enddo
    t0=-100
    do i=sams_lambda_start, sams_lambda_end
        targetProb(i)=(dvdl_square(i)-dvdl_report(i)*dvdl_report(i))/(sample_counter(i)*tt)
        if (targetProb(i) > t0) t0=targetProb(i)
    enddo  

    tt=0
    if (sams_delta .gt. 1e-4) then 
      do i=sams_lambda_start, sams_lambda_end
          if (targetProb(i) < t0*sams_delta) targetProb(i)=sams_delta*t0
          tt=tt+targetProb(i)
      enddo  
    end if

    do i=sams_lambda_start, sams_lambda_end
        targetProb(i) = targetProb(i)/tt
    enddo  

    do i=sams_lambda_start, sams_lambda_end
      resident_time(i)=targetProb(i)/t0*20+1
    enddo
  endif
  
  logPi(sams_lambda_start:sams_lambda_end)= &
    log(targetProb(sams_lambda_start:sams_lambda_end))
  pi_star=MINVAL(exp(logPi(sams_lambda_start:sams_lambda_end)))

  tt=min(j,sams_opt_end)   
  if(burnin_stage .eq. 0) then
    gamma = sams_gamma * min(pi_star, tt**(-sams_beta))
    gamma_report = sams_gamma * min(pi_star, j**(-sams_beta))
  else
    if (sams_type .eq. 0 ) then
      gamma = sams_gamma &
          * min(pi_star, (tt - onset_time + onset_time**sams_beta)**(-1))
      gamma_report = sams_gamma &
          * min(pi_star, (j - onset_time + onset_time**sams_beta)**(-1))    
    endif
  endif

  if( ((burnin_stage .eq. 0) .or. (sams_type .eq. 0)) &
            .and. (sams_type .ne. 2)  ) then

    logZ_drive(sams_lambda_start:sams_lambda_end) =  &
          logZ_drive(sams_lambda_start:sams_lambda_end) + & 
          gamma * exp(logPZ(sams_lambda_start:sams_lambda_end)- &
                      logPi(sams_lambda_start:sams_lambda_end))

    logZ_report(sams_lambda_start:sams_lambda_end) =  &
          logZ_report(sams_lambda_start:sams_lambda_end) + & 
          gamma_report * exp(logPZ(sams_lambda_start:sams_lambda_end)- &
                      logPi(sams_lambda_start:sams_lambda_end))

    tt = logZ_drive(sams_lambda_start)
    logZ_drive(sams_lambda_start:sams_lambda_end) = &
        logZ_drive(sams_lambda_start:sams_lambda_end) - tt

    tt = logZ_report(sams_lambda_start)
    logZ_report(sams_lambda_start:sams_lambda_end) = &
        logZ_report(sams_lambda_start:sams_lambda_end) - tt

  else ! TI-updating scheme
    t0=-1.0/((bar_states-1)*1.0)/temp0/KB
    logZ_drive(1)=0.0
    logZ_drive(2)=(dvdl_drive(1) + dvdl_drive(2))*t0*0.5
    do i=3, bar_states
        tt= ( dvdl_drive(i-2) + 4*dvdl_drive(i-1) + dvdl_drive(i) )/3.0 * t0
        logZ_drive(i)=logZ_drive(i-2)+tt
    enddo

    logZ_report(1)=0.0
    logZ_report(2)=(dvdl_report(1) + dvdl_report(2))*t0*0.5
    do i=3, bar_states
        tt= ( dvdl_report(i-2) + 4*dvdl_report(i-1) + dvdl_report(i) )/3.0 * t0
        logZ_report(i)=logZ_report(i-2)+tt
    enddo

  endif
 
end subroutine

subroutine sams_update_stage()

  implicit none

  logical::advance
  integer::i

  double precision::t0, tt
  integer, save::delay_counter=0
 
  if (burnin_stage.eq.1) return
  
  advance=.true.
 
  tt= smoothness(lambda_counter(sams_lambda_start:sams_lambda_end) &
      , n_considered_states) 

  if (scan_stage .eq. 1) delay_counter = delay_counter+1

  if (sams_type .eq. 2) then
    advance = (scan_stage .eq. 1)
  else 
    advance = ( (tt .le. sams_alpha) &
          .and. (scan_stage .eq. 1) &
          .and. (delay_counter.gt. n_considered_states*10) )
  endif

  if (advance) then
    burnin_stage=1
    onset_time=((1.0*currentIteration)/(1.0*sams_update_freq)-1)*1.0
  end if
  
  end subroutine

subroutine sams_update(input_dvdl, new_list)

  implicit none

  double precision, parameter:: MBAR_limit=1.0d100

  double precision,intent(in)::input_dvdl
  logical,intent(inout)::new_list
  integer::i, j, k, counter
  integer::j_start, j_end
  
  integer,save::test_counter=0, window_counter=0, rest_counter=0
  
  integer,save::update_direction=1
  
  double precision:: tt, t0
  double precision:: lambda

  !! Store the current MBAR energies
  counter=mod((currentIteration),sams_update_freq)+1
  currentIteration=currentIteration+1
  
  current_dvdl=input_dvdl

  t0=temp0*KB
  t0=log(MBAR_limit)*temp0*KB

  if (sams_type .le. 2) then
    do i=1, bar_states
      if ( abs(bar_cont(i)) .gt. t0) then
        tt=t0
      else
        tt=bar_cont(i)
      endif
      q_ave(i, counter)=exp(-tt/temp0/KB)
    end do
  endif

  sams_onstep=( counter .eq. sams_update_freq )
  if (.not. sams_onstep) return  

  !! Update du/dl average for driving du/dl for each lambda
  if (sams_type .ge. 2) then
      tt=lambda_counter(currentLambdaIndex)
      dvdl_drive(currentLambdaIndex) = dvdl_drive(currentLambdaIndex)/tt * (tt-1.0) &
        + current_dvdl/tt
  endif
  
  !! Update du/dl average for reporting du/dl for each lambda
  if (sams_do_average .and. (scan_stage.eq.1)  ) then
    sample_counter(currentLambdaIndex)=sample_counter(currentLambdaIndex)+1
    tt=sample_counter(currentLambdaIndex)
    dvdl_report(currentLambdaIndex) = dvdl_report(currentLambdaIndex)/tt * (tt-1.0)
    dvdl_report(currentLambdaIndex) = dvdl_report(currentLambdaIndex) + current_dvdl/tt

    dvdl_square(currentLambdaIndex) = &
      dvdl_square(currentLambdaIndex)/tt * (tt-1.0) &
      + current_dvdl*current_dvdl/tt   

  endif
  

  if (sams_type .le. 2) then  !! Update overall du/dl estimation  (only for SAMS run)
     
    if (mod(n_considered_states,2).eq.1) then
      tt=0.0
      tt=dvdl_report(sams_lambda_start)
      do i=sams_lambda_start+1, sams_lambda_end-1, 2
        tt=tt+dvdl_report(i)*4.0
      end do 
      do i=sams_lambda_start+2, sams_lambda_end-2, 2
        tt=tt+dvdl_report(i)*2.0
      end do 
      tt=tt+dvdl_report(sams_lambda_end)
      dvdl_ave=tt/3.0/(n_considered_states*1.0-1.0)
    
    else
      dvdl_ave=sum(dvdl_report(sams_lambda_start:sams_lambda_end)) &
        /(n_considered_states*1.0)
    endif

    dvdl_ave = dvdl_ave*(n_considered_states*1.0-1.0)/(bar_states*1.0-1.0)

    tt=sams_update_freq*1.0
    t0=0
 
    do i=1, bar_states
       logPZ(i)=logWeight(i)+log(sum(q_ave(i,:)/tt))
       logP(i)=logWeight(i)+log(q_ave(i,sams_update_freq))
    end do
    logP(:)=logP(:)-t0 
 
    !! normalize logP 
    tt=0.0
    t0=0.0
    do i=1, bar_states
       tt=tt+exp(logP(i))
       t0=t0+exp(logPZ(i))
    end do
    tt=log(tt)
    t0=log(t0)
    logP(:)=logP(:)-tt
    logPZ(:)=logPZ(:)-t0
  
  endif

  if (scan_stage .eq. 0) then
  
    if (sams_type .le. 2) then  !! normal SAMS run
        if (currentLambdaIndex .eq. sams_lambda_end) then 
          test_counter=test_counter+1
          if (test_counter>sams_scan) scan_stage=1
        else
            window_counter=window_counter+1
            if (logP(currentLambdaIndex+1)-logP(currentLambdaIndex) .gt. -2.0 )&
                test_counter=test_counter+1
            if (window_counter.gt.(sams_scan+100)) &
                test_counter=sams_scan+1
        endif
        if (test_counter>=sams_scan) then 
            window_counter=0
            test_counter=0
            currentLambdaIndex=currentLambdaIndex+1
        endif
    else 
      
        !!NFE case
        test_counter=test_counter+1
        if (sams_direction .eq.0) then
          k=currentLambdaIndex+1
        else
          k=currentLambdaIndex-1
        endif
        tt=bar_cont(k)-bar_cont(currentLambdaIndex)
        q_ave(k, test_counter)=exp(-tt/temp0/KB)
       
        if ((sams_direction.eq.0 .and. currentLambdaIndex .eq. sams_lambda_end-1) &
          .or. (sams_direction.eq.1 .and. currentLambdaIndex .eq. sams_lambda_start+1) ) then 
          if (test_counter>=sams_scan)  then
            call sams_update_NFE()
            new_list=.true.
            test_counter=0
          endif
        else
          if (test_counter>=sams_scan) then 
              window_counter=0
              test_counter=0
              if (sams_direction.eq.0) then 
                currentLambdaIndex=currentLambdaIndex+1
              else
                currentLambdaIndex=currentLambdaIndex-1
              endif
          endif
        endif
    endif

  else 
    
    if(sams_variance .eq.2 .and. rest_counter.le. resident_time(currentLambdaIndex) ) then
      rest_counter = rest_counter+1
    else
      rest_counter=0 
   
      j_start=max(currentLambdaIndex-sams_jump_range, sams_lambda_start)
      j_end=min(currentLambdaIndex+sams_jump_range, sams_lambda_end)    
      
      if (sams_direction .ne. 0) then 
        if (update_direction .eq. 1) j_start=max(j_start, currentLambdaIndex)
        if (update_direction .eq. -1) j_end=min(j_end, currentLambdaIndex)        
      endif
      !! normalize jumpProb
      tt=0.0
      do i=j_start, j_end
        jumpProb(i)=exp(logP(i)) 
        if (i .eq. currentLambdaIndex) jumpProb(i)=0
        tt=tt+jumpProb(i)
      end do
      jumpProb(j_start:j_end)=jumpProb(j_start:j_end)/tt

      call amrand(tt)
      j=-1
      do i=j_start, j_end
        if (tt<jumpProb(i)) then
          j=i
          exit
        end if
        tt=tt-jumpProb(i)
      end do  
      if (j<sams_lambda_start .or. j>sams_lambda_end) then
        j=currentLambdaIndex
      endif
      
      if (sams_direction .ne. 0) then 
        if (update_direction .eq. 1) j=max(j, currentLambdaIndex)
        if (update_direction .eq. -1) j=min(j, currentLambdaIndex)        
        if (sams_direction .eq. 1) then 
          if ((currentLambdaIndex .ne. sams_lambda_start) .and. (j.eq. sams_lambda_start)) &
            update_direction=1
          if ((currentLambdaIndex .ne. sams_lambda_end) .and. (j.eq. sams_lambda_end)) &
            update_direction=-1     
        else if (sams_direction .eq. 2) then 
          if ( (j-sams_lambda_start) .le. sams_jump_range/5 ) &
            update_direction=1
          if ( (sams_lambda_end-j) .le. sams_jump_range/5 ) &
            update_direction=-1     
        endif
      end if
      
      ! Note that we only collect transitions due to MC moves.
      transition_counter(currentLambdaIndex, j)=transition_counter(currentLambdaIndex, j)+1
 
      currentLambdaIndex=j
    endif
  endif
  
  lambda=mbar_lambda(currentLambdaIndex)
  call ti_update_lambda(lambda, currentLambdaIndex)
  
  lambda_counter(currentLambdaIndex)=lambda_counter(currentLambdaIndex)+1

  if (sams_type .le. 2) then  !! normal SAMS run
    call sams_update_stage()
    call sams_update_Z_estimate()
    call sams_update_weight()
  endif

end subroutine

subroutine sams_output(currentStep)

  implicit none

  integer, intent(in)::currentStep
  integer::outUnit, endUnit 
  integer::i, j, tot, info
  double precision::t0,tt

  call amopen(sams_log, sams_log_name, 'U', 'F', 'A')

  if (burnin_stage .eq. 0) then 
    endUnit=sams_log
  else 
    endUnit=sams_rest
    call amopen(sams_rest, sams_restart_name, 'U', 'F', 'W') 
  endif
  
  do outUnit = sams_log, endUnit, max(1,(endUnit-sams_log))
 
    if (burnin_stage .eq. 0) then 
      t0=sum(lambda_counter)  
    else
      t0=sum(sample_counter)  
    end if

    tt=smoothness(sample_counter, n_considered_states)

    write(outUnit, *)
    write(outUnit, *)  
    write(outUnit, form_iter) currentStep, currentIteration, &
        sum(sample_counter(1:bar_states)) , currentLambdaIndex, &
        mbar_lambda(currentLambdaIndex), burnin_stage     
    if (sams_type .le. 2) then
        if (sams_type .eq. 2) then
          write(outUnit, form_estimate) &
             (logZ_drive(sams_lambda_start)-logZ_drive(sams_lambda_end))*temp0*KB, tt, &
             current_dvdl, dvdl_ave
        else
          write(outUnit, form_estimate) &
           (logZ_report(sams_lambda_start)-logZ_report(sams_lambda_end))*temp0*KB, tt, &
             current_dvdl, dvdl_ave
        endif
           
        if ((sams_verbose.ge.2) .and. (outUnit.ne.sams_rest)) then  
          write(outUnit, *) "#MBAR Energy: (Relative to current state)"
          write(outUnit, form_mbar) bar_cont(1:bar_states)    
        endif
     
        if ((sams_verbose.ge.2) .or. (outUnit.eq.sams_rest)) then  

          write(outUnit, *) "#Lambda counter:"            
          write(outUnit, form_counter) lambda_counter(1:bar_states)
    
          write(outUnit, *) "#Sample counter:", tt             
          write(outUnit, form_counter) sample_counter(1:bar_states)    

          if (burnin_stage.eq.1 .and. sams_verbose.ge.3) then
            ! Compute empirical transition matrix under long-term reversibility
            ! assumption. This is well-defined for a given trajectory, but may
            ! be a (rapidly) moving target for shorter simulations.
            !
            write(outUnit, *) "#Transition matrix:"
            do i=1, bar_states
              do j=1, bar_states
                tot=sum(transition_counter(i, :)) + sum(transition_counter(:, i))
                if ((tot.le.0)) then
                  transition_matrix(i, j)=0.0
                else
                  tt=dble(transition_counter(i, j)) + dble(transition_counter(j, i))
                  transition_matrix(i, j)=tt / dble(tot)
                endif
              end do
              write(outUnit, form_prob) transition_matrix(i, 1:bar_states)
            end do
            call dgeev('N','N',bar_states,transition_matrix,bar_states, &
                       eigval,eigval_imag,eigvec,bar_states,eigvec, &
                       bar_states,eigscratch,4*bar_states,info)
    !        write(outUnit, *) "#Eigenvalues:"
    !        write(outUnit, form_prob) eigval(1:bar_states)
            ! Find the largest non-unit eigenvalue.
            t0=0.0
            do i=1, bar_states
              ! exponential time (parameter in correlation function)
              !tt = -1 / log(abs(eigval(i)))
              ! long-time approximation (best near ~1, else upper bound)
              !tt = 1 / (1 - eigval(i))
              ! integrated time (damped integral of correlation function)
              tt = 0.5*(1 + eigval(i)) / (1 - eigval(i))
              ! 1e5 should filter out values > 0.9999
              if ((t0.lt.tt) .and. (tt.lt.1e5)) then
                t0=tt
              endif
            end do
            t0=t0*bar_intervall*dt
            write(outUnit, "('#Correlation time: ',F10.3,' ps')") t0
          endif

          write(outUnit, *) "#logZ:"
          write(outUnit, form_logZ) logZ_drive(1:bar_states)

          write(outUnit, *) "#Current dV/dL average:", dvdl_ave
          write(outUnit, form_dvdl) dvdl_report(1:bar_states)   
    
        endif
    endif   
  
    close(outUnit)
  
  end do

end subroutine

subroutine sams_read_init(read_init)

  implicit none

  logical, intent(inout)::read_init
  double precision:: delta=1e-4, t1, t2, t3
  integer::inUnit 
  integer::stat_counter
  
  call amopen(sams_init, sams_init_name, 'O', 'F', 'R') 
  
  inUnit = sams_init
  
  read_init=.true.
  
  read(inUnit, *)
  read(inUnit, *)  
  read(inUnit, form_iter, err=100) sams_restart_step, currentIteration,  &
                     stat_counter, currentLambdaIndex, t1, burnin_stage     
  read_init = read_init .and. (abs(t1-mbar_lambda(currentLambdaIndex)) < delta)   
  
  if (burnin_stage .eq. 1) scan_stage=1
 
  read(inUnit, form_estimate, err=100) &
              t1, t2, current_dvdl, dvdl_ave    
  
  read(inUnit, *, err=100) !"#Lambda counter:"             
  read(inUnit, form_counter, err=100) lambda_counter(1:bar_states)
    
  read(inUnit, *, err=100) !"#Sample counter:"             
  read(inUnit, form_counter, err=100) sample_counter(1:bar_states)
  read_init = read_init .and. &
      (stat_counter .eq. sum(sample_counter(1:bar_states)) )  
       
  read(inUnit, *, err=100) !"#logZ:"
  read(inUnit, form_logZ) logZ_drive(1:bar_states)
    
  read(inUnit, *, err=100) !"#Current dV/dL average:"
  read(inUnit, form_dvdl, err=100) dvdl_report(1:bar_states)     
 
  close (sams_init)
  return

100 read_init=.false.
 
end subroutine

subroutine sams_update_NFE()
  
  implicit none
  integer       :: i,j,k,l,m
  double precision :: tt, t0, t1, t2, box(3), ratio(3)

  logical       :: read_init

  if (sams_direction .eq. 0) then
    i=2
    j=bar_states
    k=1
  else
    i=bar_states-1
    j=1
    k=-1
  endif
    
  t0=0.0
  t1=0.0
  call amopen(sams_log, sams_log_name, 'U', 'F', 'A')
  
  do l=i, j, k
    tt=-KB*temp0*log( sum(q_ave(l, 1:sams_scan))/sams_scan )
    t0=t0+tt
  
    write(sams_log,'("NFE: Window #",I5, "  :  work = ",F10.6, "     Ave. du/dl = ", F10.6 )') l, tt, dvdl_drive(l)
  enddo
  
  t1=1.0/(bar_states-1.0)
  do i=3, bar_states,2
        t2= t2+( dvdl_drive(i-2) + 4*dvdl_drive(i-1) + dvdl_drive(i) )/3.0 * t1
  enddo

  write(sams_log,'("NFE: Total:  work = ",F10.6, "  Work2  ", F10.6, "    Ave. du/dl = ", F10.6 )')  &
    t0, -KB*temp0*sum(log(q_ave(2:bar_states,sams_scan))), t2

  if (sams_direction.eq.0) then 
    currentLambdaIndex=sams_lambda_start
  else
    currentLambdaIndex=sams_lambda_end
  endif
  dvdl_drive(:)=0.0
  lambda_counter(:)=0
  
  call amrand_gen(sams_rand_gen, my_random_value) 
  rremd_idx(1) = min(reservoir_size(1),  &
      int(my_random_value * reservoir_size(1)) + 1)
  !!write(*,*) "Loading reservoir", rremd_idx, local_atm_cnt
  !!TBT
  !!rremd_idx(1) =4147
  call load_reservoir_structure(atm_crd, atm_vel, local_atm_cnt, &
      box(1), box(2), box(3), 1)
  atm_last_vel=atm_vel
  
  ratio(:)=box(:)/pbc_box(:)
  
  if (abs(ratio(1)-ratio(2))>1e-6 .or. abs(ratio(2)-ratio(3))>1e-6 .or. &
    abs(ratio(1)-ratio(3))>1e-6 ) then
    write(*,*)' Wrong box info in the reservoir'
    call mexit()
  else
    call pressure_scale_pbc_data(ratio, vdw_cutoff + skinnb , 0)
  endif  
  
#ifdef CUDA  
  call gpu_upload_crd(atm_crd)
  call gpu_upload_vel(atm_vel)
  call gpu_upload_last_vel(atm_last_vel)
  if (ntb .ne. 1) then
    call gpu_ucell_set(ucell, recip, uc_volume)
  endif
  call gpu_force_new_neighborlist()
#endif  

  write(sams_log,'("NFE: Load Reservoir # ",I6)') rremd_idx(1) 
  
  close(sams_log)
  
end subroutine

end module


