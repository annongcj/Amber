
        use mdin_ctrl_dat_mod
        use parallel_dat_mod
        use prmtop_dat_mod
        use gbl_datatypes_mod
        use timers_mod
        use ti_mod

        implicit none

        ! Formal arguments:

#ifdef __INTEL_COMPILER
#define LENGTH atm_cnt
#else
#define LENGTH *
#endif /* INTEL_COMPILER */

        integer               :: atm_cnt
        double precision      :: crd(*)
        double precision      :: rborn(LENGTH)
        double precision      :: fs(LENGTH)
        double precision      :: charge(LENGTH)
        integer               :: iac(LENGTH)
        integer               :: ico(LENGTH)
#ifdef __INTEL_COMPILER
        integer               :: numex(*) 
        integer               :: natex(*)
#else
        integer               :: numex(LENGTH) 
        integer               :: natex(LENGTH)
#endif
        integer               :: natbel
        integer, intent(in)   :: irespa
        logical               :: skip_radii ! for converting the optional argument

        double precision      :: egb, eelt, evdw, esurf
        double precision      :: frc(*)
        logical, optional, intent(in)         :: skip_radii_

        ! Local variables:

        GBFloat               :: cut2, cut_inner2
        GBFloat               :: extdiel_inv
        GBFloat               :: intdiel_inv
        GBFloat               :: ri, rj
        GBFloat               :: ri1i
        GBFloat               :: xij, yij, zij
        GBFloat               :: dij1i, dij2i, dij3i
        GBFloat               :: r2
        GBFloat               :: dij
        GBFloat               :: sj, sj2
        GBFloat               :: frespa
        GBFloat               :: qi, qiqj
        GBFloat               :: dumx, dumy, dumz
        GBFloat               :: sumi            
        integer               :: x, counter, start_iter, end_iter, Jmin
        integer               :: icount
        GBFloat               :: fgbi
        GBFloat               :: rinv, r2inv, r6inv, r10inv
        GBFloat               :: fgbk
        GBFloat               :: expmkf
        GBFloat               :: dl
        GBFloat               :: de
        GBFloat               :: e
        GBFloat               :: temp1
        GBFloat               :: temp4, temp5, temp6, temp7
        GBFloat               :: eel
        GBFloat               :: f6, f12, f10
        GBFloat               :: dedx, dedy, dedz
        GBFloat               :: qi2h, qid2h
        GBFloat               :: datmp
        GBFloat               :: thi, thi2
        GBFloat               :: f_x, f_y, f_z
        GBFloat               :: f_xi, f_yi, f_zi
        GBFloat               :: xi, yi, zi
        GBFloat               :: dumbo
        GBFloat               :: tmpsd
        integer               :: thread_id

        ! Variables needed for smooth integration cutoff in Reff:

        GBFloat               :: rgbmax1i
        GBFloat               :: rgbmax2i
        GBFloat               :: rgbmaxpsmax2

#ifdef SOFTCORE_TI
        ! Variables needed for softcore LJ terms in TI:

        double precision      :: sc_vdw                       !
        double precision      :: r6                           !

        ! Variables needed for softcore eel terms in TI:
        double precision      :: sc_eel                       !
        double precision      :: sc_eel_denom                 ! 
        double precision      :: sc_eel_denom_sqrt            !
        integer               :: sceeorderinv

        integer               :: num_ti_atms_cntd_i 
        integer               :: num_ti_atms_cntd_j 

        logical       :: i_is_sc 
        logical       :: j_is_sc 
#endif ! SOFTCORE_TI

        ! Scratch variables used for calculating neck correction:

        GBFloat               ::  mdist
        GBFloat               ::  mdist2
        GBFloat               ::  mdist3
        GBFloat               ::  mdist5
        GBFloat               ::  mdist6

        ! Stuff for alpb:

        double precision      :: alpb_beta
        double precision      :: one_arad_beta
        double precision      :: gb_kappa_inv

        ! Alpha prefactor for alpb_alpha:
        double precision, parameter   :: alpb_alpha = 0.571412d0

        integer               :: neibr_cnt
        integer               :: i, j, k
        integer               :: kk1,kod
        integer               :: max_i
        integer               :: iaci
        integer               :: iexcl, jexcl
        integer               :: jexcl_last
        integer               :: jjv
        integer               :: ic
        integer               :: j3
        logical               :: onstep
        GBFloat               :: si, si2
        GBFloat               :: theta
        GBFloat               :: uij
        GBFloat               :: reff_i
        integer               :: nstart 

        ! Mask variables to remove gb_kappa and alpb conditional equations

        integer               :: outer_k
        integer               :: gb_kappa_zero
        integer               :: alpb_zero
        GBFloat               :: alpb_factor
        GBFloat               :: gb_kappa_factor
        GBFloat               :: use_nrespa

        integer               :: index_icj

        ! Scratch variables used for calculating neck correction:

        double precision      ::  neck


        ! Variables for blocking
        
        integer               :: outer_i, jend

        !Variables for atom distribution


        integer               :: temp_icount
        GBFloat               :: r2x_offdiag(BATCH_SIZE)

        ! FGB taylor coefficients follow
        ! from A to H :
        ! 1/3 , 2/5 , 3/7 , 4/9 , 5/11
        ! 4/3 , 12/5 , 24/7 , 40/9 , 60/11

        GBFloat               :: vt1, vt2, vt3, vt4, vt5 
        
        skip_radii = (irespa .gt. 1 .and. mod(irespa, nrespai) .ne. 0)
  
        ! For gas phase calculations gbradii are not needed so set skip to true.
        if (igb == 6) skip_radii = .true.
 
       if (present(skip_radii_)) &
        skip_radii = skip_radii_ .or. skip_radii

        if (mod(irespa, nrespai) .ne. 0) return

        egb = 0.d0
        eelt = 0.d0
        evdw = 0.d0
        esurf = 0.d0

        cut2 = gb_cutoff * gb_cutoff
        cut_inner2 = cut_inner * cut_inner
        onstep = mod(irespa, nrespa) .eq. 0

        if (gb_kappa .eq. 0.d0) then
           gb_kappa_zero = 1 
        else
           gb_kappa_zero = 0 
        end if

        if (alpb .eq. 0) then
        ! Standard Still's GB
          alpb_zero = 1 
          extdiel_inv = ONE / extdiel
          intdiel_inv = ONE / intdiel
          one_arad_beta = 0.d0
          gb_kappa_inv = 0.d0
        else
           alpb_zero = 0 
        ! Sigalov Onufriev ALPB (epsilon-dependent GB):
           alpb_beta = alpb_alpha * (intdiel / extdiel)
           extdiel_inv = ONE / (extdiel * (ONE + alpb_beta))
           intdiel_inv = ONE / (intdiel * (ONE + alpb_beta))
           one_arad_beta = alpb_beta / arad
           if (gb_kappa .ne. 1.d0) gb_kappa_inv = 1.d0 / gb_kappa
        end if

        max_i = atm_cnt

        if (natbel .gt. 0) max_i = natbel

#ifdef TIMODE
        do j = 1, atm_cnt
         if(ti_lst(ti_mask_piece,j) .ne. 0) then
            ti_mode_skip(j) = ZERO 
         else 
           ti_mode_skip(j) = ONE
         end if
       end do  
#endif
        ! Smooth "cut-off" in calculating GB effective radii.
        ! Implemented by Andreas Svrcek-Seiler and Alexey Onufriev.
        ! The integration over solute is performed up to rgbmax and includes
        ! parts of spheres; that is an atom is not just "in" or "out", as
        ! with standard non-bonded cut.  As a result, calculated effective
        ! radii are less than rgbmax. This saves time, and there is no
        ! discontinuity in dReff / drij.

        ! Only the case rgbmax > 5*max(sij) = 5*gb_fs_max ~ 9A is handled; this is
        ! enforced in mdread().  Smaller values would not make much physical
        ! sense anyway.

        rgbmax1i = ONE / rgbmax
        rgbmax2i = rgbmax1i * rgbmax1i
        rgbmaxpsmax2 = (rgbmax + gb_fs_max)**2

        !------------------------------------------------------------------
        ! Check if some MPI ranks are completing before the others
        ! Frequency of this check is currently set to 10 seconds 
        ! This should ideally be 10 seconds or time to complete around 20
        ! 20 iterations - whichever is more
        ! If there is imbalance, the distribute_atoms routine is called 
        !-----------------------------------------------------------------

        call get_wall_time(wall_s, wall_u)
        strt_time_ms = dble(wall_s) * 1000.d0 + dble(wall_u) / 1000.d0

        if(numtasks .gt. 1) then
              if((strt_time_ms - last_balance_time) .ge. 10000.d0) then
                  time_to_balance = .true.
              else
                  time_to_balance = .false.
              end if

              call mpi_bcast(time_to_balance, 1, mpi_logical, 0, pmemd_comm,&
                           err_code_mpi)

             if(time_to_balance) then
                 call check_imbalance
                 last_balance_time = strt_time_ms
              end if

              if (atm_cnt .ne. prev_atm_cnt) then
                 atmcnt_change = .true.
                 last_balance_time = strt_time_ms
              end if

              if (atmcnt_change .or. imbalance_radii .or. imbalance_offd) then
                call distribute_atoms(atm_cnt)
                prev_atm_cnt = atm_cnt
                imbalance_radii = .false.
                imbalance_offd = .false.
                atmcnt_change = .false.
                time_to_balance = .false.
              end if
        end if

        !---------------------------------------------------------------------------
        ! Step 1: loop over pairs of atoms to compute the effective Born radii.
        !---------------------------------------------------------------------------

        ! OpenMP region starts for Calc-Born Radii

        !$omp parallel default(shared) &
        !$omp& private(xi,yi,zi,reff_i,nstart,temp_jj,temp_r2,&
        !$omp& ri,ri1i,si,si2,temp_icount,xij,yij,zij,r2,vectmp1,sj,sj2,dij1i,dij,uij,dij2i,tmpsd,dumbo,vt2,&             
        !$omp& vt4,theta,mdist,mdist2,mdist3,mdist6,neck)
       
        if(.not. skip_radii) then
          !$omp do
          do i= 1, atm_cnt
             reff(i) = 0.d0
          end do
          !$omp end do
        end if
          
       !$omp do
       do i= 1, atm_cnt
             crdx(i) = ToGBFloat(crd( 3 * i - 2 ))
             crdy(i) = ToGBFloat(crd( 3 * i - 1 ))
             crdz(i) = ToGBFloat(crd( 3 * i ))
              rborn_lcl(i) = ToGBFloat(rborn(i))
              fs_lcl(i) = ToGBFloat(fs(i))
              charge_lcl(i) = ToGBFloat(charge(i))
        end do
        !$omp end do

        !$omp do schedule(dynamic)
        do outer_i = start_atm_radii, end_atm_radii, BLKSIZE
           do i = outer_i, MIN0(outer_i+BLKSIZE, end_atm_radii+1)-1
#ifdef TIMODE
                if (ti_lst(ti_mask_piece,i) .ne. 0) cycle
#endif
              nstart = (i -start_atm_radii )* max_nlist_count 
              xi = crdx(i)
              yi = crdy(i)
              zi = crdz(i)

              reff_i = reff(i)
              ri = rborn_lcl(i) - offset ! convert to float
              ri1i = ONE / ri
              si = fs_lcl(i)
              si2 = si * si

              temp_icount = 0

        ! Here, reff_i will sum the contributions to the inverse effective
        ! radius from all of the atoms surrounding atom "i"; later the
        ! inverse of its own intrinsic radius will be added in

#ifndef TIMODE
                do j=1, atm_cnt
                   xij = xi - crdx(j)
                   yij = yi - crdy(j)
                   zij = zi - crdz(j)
                   r2 = xij * xij + yij * yij + zij * zij

                   if( j .eq. i ) cycle
                   if (r2 .gt. rgbmaxpsmax2) cycle
                   temp_icount = temp_icount + 1
                   temp_jj(temp_icount) = j
                   temp_r2(temp_icount) = r2
                end do
#else
                do j = 1, atm_cnt
                   if (ti_lst(ti_mask_piece,j) .ne. 0) cycle
                   xij = xi - crdx(j)
                   yij = yi - crdy(j)
                   zij = zi - crdz(j)
                   r2 = xij * xij + yij * yij + zij * zij

                   if( j .eq. i ) cycle
                   if (r2 .gt. rgbmaxpsmax2) cycle
                   temp_icount = temp_icount + 1
                   temp_jj(temp_icount) = j
                   temp_r2(temp_icount) = r2
                end do
#endif

             if(temp_icount .gt. max_nlist_count) then
               print *, " neighbour list size error"
             end if
             if (calc_nlist) then
                nlist_count(i) = temp_icount
                nlist_jj(nstart+1:nstart+temp_icount) = temp_jj(1:temp_icount)
                nlist_r2(nstart+1:nstart+temp_icount) = temp_r2(1:temp_icount)
             end if
      
             if (.not. skip_radii) then
               vectmp1(1:temp_icount) = ONE/sqrt(temp_r2(1:temp_icount))
               !dir$ simd reduction(+:reff_i)
                do k=1, temp_icount
                   j = temp_jj(k) 
                   r2 = temp_r2(k) 
                   sj = fs_lcl(j)
                   sj2 = sj * sj

        ! don't fill the remaining vectmp arrays if atoms don't see each other:

                  dij1i = vectmp1(k)
                  dij = r2 * dij1i

                  if (dij .le. rgbmax + sj) then
                     if ((dij .gt. rgbmax - sj)) then
                         uij = ONE / (dij -sj)
                         reff_i = reff_i - EIGHTH * dij1i * (ONE + &
                                  TWO * dij *uij +  rgbmax2i * &
                                  (r2 - FOUR * rgbmax * dij - sj2) + &
                                  TWO * log((dij - sj) * rgbmax1i))
                     else if (dij .gt. FOUR * sj) then
                        dij2i = dij1i * dij1i
                        tmpsd = sj2 * dij2i
                        dumbo = ta + tmpsd *  (tb + tmpsd * (tc + tmpsd * &
                                 (td + tmpsd * tdd)))
                        reff_i = reff_i - tmpsd * sj * dij2i * dumbo
                        else
                          vt2 = ONE/(dij + sj)
                          if (dij .gt. ri + sj) then
                            vt4 = log(vt2 * (dij -sj) )
                            reff_i =  reff_i - HALF * (sj / (r2 - sj2) + &
                                      HALF * dij1i * vt4)
                          else if (dij .gt. abs(ri - sj)) then
                            vt4 = log(vt2 * ri)
                            theta = HALF * ri1i * dij1i * (r2 + ri * ri - sj2)
                            reff_i = reff_i - QUARTER * (ri1i * &
                                    (TWO - theta) -  vt2 + dij1i * vt4)
                          else if (ri .lt. sj) then
                            vt4 = log(vt2 * (sj - dij) )
                            reff_i = reff_i - HALF * (sj / (r2 - sj2) + &
                                   TWO * ri1i + HALF * dij1i * vt4)
                           end if !  if (dij .gt. ri + sj) then
                        end if  ! (dij .gt. 4.d0 * sj)

                      if (igb78) then
                         if (dij .lt. rborn_lcl(i) + rborn_lcl(j) + gb_neckcut) then
                            mdist = dij - neckMaxPos(neck_idx(i), neck_idx(j))
                            mdist2 = mdist * mdist
                            mdist3 = mdist2 * mdist
                            mdist6 = mdist3 * mdist3
                            neck = neckMaxVal(neck_idx(i), neck_idx(j)) / &
                            (ONE + mdist2 + THREE * TENTH * mdist6)
                            reff_i = reff_i - gb_neckscale * neck
                          end if
                      end if ! igb .eq. 7 .or. 8
                  end if
               end do ! k = 1, temp_icount
               reff(i) = reff_i
            end if ! skip_radii
         end do ! i = outer_i , outer_i + BLKSIZE
      end do  !  outer_i = 1, atm_cnt
!$omp end do
!$omp end parallel

! OpenMP region ends for Calc-Born Radii

    call get_wall_time(wall_s, wall_u)
    loop_time = dble(wall_s) * 1000.d0 + dble(wall_u) / 1000.d0 - &
                strt_time_ms
#ifdef GBTimer
    print *, "Radii =", loop_time 
#endif
    radii_time(mytaskid + 1) = radii_time(mytaskid +1 ) + loop_time


if (.not. skip_radii) then
     
#ifdef MPI
       call update_gb_time(calc_gb_rad_timer)

       ! Collect the (inverse) effective radii from other nodes:
if (numtasks .gt. 1) then
        call mpi_allreduce(reff, red_buf, atm_cnt, mpi_double_precision, &
                        mpi_sum, pmemd_comm, err_code_mpi)

        reff(1:atm_cnt) = red_buf(1:atm_cnt)
end if        

        call update_gb_time(dist_gb_rad_timer)
#endif

        if (igb2578) then

        ! apply the new Onufriev "gbalpha, gbbeta, gbgamma" correction:
        
        do i = 1, atm_cnt
#ifdef TIMODE
                if (ti_lst(ti_mask_piece,i) .ne. 0) cycle
#endif
                ri = rborn_lcl(i) - offset
                ri1i = 1.d0 / ri
                psi(i) = -ri * reff(i)
                reff(i) = ri1i - tanh((gb_alpha_arry(i) + gb_gamma_arry(i) * psi(i) * &
                                        psi(i) - gb_beta_arry(i) * psi(i)) * psi(i)) / rborn_lcl(i)

                if (reff(i) .lt. 0.d0) reff(i) = THIRTEETH 

                reff(i) = 1.d0 / reff(i)
        end do

        else

        ! "standard" GB, including the "diagonal" term here:
        do i = 1, atm_cnt
#ifdef TIMODE
                if (ti_lst(ti_mask_piece,i) .ne. 0) cycle
#endif
                ri = rborn_lcl(i) - offset
                ri1i = 1.d0 / ri
                reff(i) = 1.d0 / (reff(i) + ri1i)
        end do
        end if

        if (rbornstat .eq. 1) then
        do i = 1, atm_cnt
#ifdef TIMODE
                if (ti_lst(ti_mask_piece,i) .ne. 0) cycle
#endif
                gbl_rbave(i) = gbl_rbave(i) + reff(i)
                gbl_rbfluct(i) = gbl_rbfluct(i) + reff(i) * reff(i)
                if (gbl_rbmax(i) .le. reff(i)) gbl_rbmax(i) = reff(i)
                if (gbl_rbmin(i) .ge. reff(i)) gbl_rbmin(i) = reff(i)
        end do
        end if

        call update_gb_time(calc_gb_rad_timer)

end if /*skip_radii*/
        !--------------------------------------------------------------------------

        !--------------------------------------------------------------------------
        !
        ! Step 2: Loop over all pairs of atoms, computing the gas-phase
        !         electrostatic energies, the LJ terms, and the off-diagonal
        !         GB terms.  Also accumulate the derivatives of these off-
        !         diagonal terms with respect to the inverse effective radii,
        !         sumdeijda(k) will hold  sum over i, j>i (deij / dak),  where
        !         "ak" is the inverse of the effective radius for atom "k".
        !
        !         Update the forces with the negative derivatives of the
        !         gas-phase terms, plus the derivatives of the explicit
        !         distance dependence in Fgb, i.e. the derivatives of the
        !         GB energy terms assuming that the effective radii are constant.
        !
        !--------------------------------------------------------------------------

        ! Note: this code assumes that the belly atoms are the first natbel
        !       atoms...this is checked in mdread.

        iexcl = 1
        iexcl_arr(1) = 1
        do k = 2, atm_cnt
             iexcl_arr(k) = iexcl_arr(k-1) + numex(k-1)
        end do
        
        sumdeijda(1:atm_cnt) = 0.d0
        call get_wall_time(wall_s, wall_u)
        strt_time_ms = dble(wall_s) * 1000.d0 + dble(wall_u) / 1000.d0

! OpenMP region starts for Off-Diagonal

        !$omp parallel default(shared) &
        !$omp& private(iexcl,xi,yi,zi,qi,ri,iaci,jexcl,jexcl_last,dumx,dumy,dumz,sumi,skipv,r2,thread_id,&
        !$omp& xij,yij,zij,x,jmin,jend,counter,r2x_offdiag,vectmp1,vectmp3,vectmp2,jj_od,r2x_od, rjx_od,icount,&
        !$omp& vectmp4,vectmp5,kk1,index_icj,rj,qiqj,expmkf,fgbk,dl,fgbi,e,temp4,temp6,temp1,de,temp5,rinv,r2inv, &
#ifdef SOFTCORE_TI
        !$omp& r6,num_ti_atms_cntd_i,num_ti_atms_cntd_j,i_is_sc,j_is_sc, sc_eel_denom,sc_eel_denom_sqrt, &
#endif ! SOFTCORE_TI
        !$omp& r6inv,eel,ic,f6,f12,dedx,dedy,dedz,start_iter,vt1,vt2,vt3,vt4,vt5,end_iter) &
#ifdef __INTEL_COMPILER
        !$omp& firstprivate(gb_kappa_factor,gb_kappa_zero, alpb_zero, alpb_factor,nrespa, nrespai, use_nrespa)
#else
        !$omp& firstprivate(gb_kappa_factor,gb_kappa_zero, alpb_zero, alpb_factor) 
#endif

        thread_id = omp_get_thread_num() + 1

        frcx_red(1:atm_cnt,thread_id) = ZERO
        frcy_red(1:atm_cnt,thread_id) = ZERO 
        frcz_red(1:atm_cnt,thread_id) = ZERO 

#ifdef SOFTCORE_TI
  !if (ifsc .gt. 0) then
  !print ti_region
    if (ti_mask_piece .eq. 2) then
      sc_vdw = scalpha * clambda !vdw
      sc_eel = scbeta * clambda !+ r2 !eel
    else if (ti_mask_piece .eq. 1) then
      sc_vdw = scalpha * (1.d0 - clambda)
      sc_eel = scbeta * (1.d0 - clambda) !+ r2 
    end if

    num_ti_atms_cntd_i = 0 
#endif ! SOFTCORE_TI

#ifdef SOFTCORE_TI
        !$omp do reduction(+:eelt,evdw) schedule(dynamic)
#else
#ifdef GB_ENERGY
        !$omp do reduction(+:eelt,egb,sumdeijda,evdw) schedule(dynamic)
#else
        !$omp do reduction(+:sumdeijda) schedule(dynamic)
#endif
#endif ! SOFTCORE_TI
       do outer_i=start_atm_offd, end_atm_offd,S2_BATCH_SIZE
          do i = outer_i, MIN0(outer_i+S2_BATCH_SIZE, end_atm_offd+1)-1
            iexcl = iexcl_arr(i)

#ifdef TIMODE
               if (ti_lst(ti_mask_piece,i) .ne. 0) then
                  iexcl = iexcl + numex(i)
                  cycle
               end if
#endif

#ifdef SOFTCORE_TI
        if (ti_lst(ti_region,i) .ne. 0) then 
          num_ti_atms_cntd_i = num_ti_atms_cntd_i + 1
          if (ti_sc_lst(num_ti_atms_cntd_i) .eq. 1) i_is_sc = .true. 
        end if         

        num_ti_atms_cntd_j = 0 

#endif ! SOFTCORE_TI

             xi = crdx(i)
             yi = crdy(i)
             zi = crdz(i)
          icount = 0      
          do j = i + 1, atm_cnt
#ifdef TIMODE
               if (ti_lst(ti_mask_piece,j) .ne. 0) cycle
#endif

#ifdef SOFTCORE_TI
        if (ti_lst(ti_region,j) .ne. 0) then
          num_ti_atms_cntd_j = num_ti_atms_cntd_j + 1
          if (ti_sc_lst(num_ti_atms_cntd_j) .ne. 0) j_is_sc = .true.
        end if
#endif ! SOFTCORE_TI

               xij = xi - crdx(j)
               yij = yi - crdy(j)
               zij = zi - crdz(j)
               r2 = xij * xij + yij * yij + zij * zij
               if (r2 .gt. cut2) cycle
               if (.not. onstep .and. r2 .gt. cut_inner2) cycle
               icount = icount + 1
               jj_od(icount) = j
               r2x_od(icount) = r2
               !if ( igb /= 6) rjx_od(icount) = reff(j)
               rjx_od(icount) = reff(j)
#ifdef SOFTCORE_TI
        if (ifsc .gt. 0) then
          sc_eel_denom = 1/(sc_eel + r2)
          sc_eel_denom_sqrt = sc_eel_denom ** sceeorderinv 
        end if
#endif ! SOFTCORE_TI
          end do
             
             qi = charge_lcl(i)
             !if ( igb /= 6 ) ri = reff(i)
             ri = reff(i)
             iaci = ntypes * iac(i) 
             jexcl = iexcl
             jexcl_last = iexcl + numex(i) - 1

             dumx = ZERO
             dumy = ZERO
             dumz = ZERO
             sumi = ZERO

             ! check the exclusion list for eel and vdw:
             !$dir ivdep
             do k = i + 1, atm_cnt
                skipv(k) = 1.0
             end do

             !$dir ivdep
             do jjv = jexcl, jexcl_last
                skipv(natex(jjv)) = 0.0
             end do
!dir$ ivdep
        do k = 1, icount
           j = jj_od(k)
           r2 = r2x_od(k)
           rj = rjx_od(k)
        !TODO_IGB6
           kk1 = iaci + iac(j)

           xij = xi - crdx(j)
           yij = yi - crdy(j)
           zij = zi - crdz(j)
!           r2 = xij * xij + yij * yij + zij * zij
 
#ifndef SOFTCORE_TI
    if (igb /= 6) then
           !vt1 = FOUR * !ri * rj
           ! !   vt1 = 1.d0 / vt1
           !vt1 = -r2 / vt1
           !vt1 = exp(vt1)
           ! Rewritting the above calculations to remove frequent writebacks
           ! to the main memory
           vt1 = exp(-r2 / (FOUR * ri * rj))
            ! vectmp1 now contains exp(-rij^2/[4*ai*aj])
           vt3 = r2 + rj * ri * vt1
            ! vectmp3 now contains fij
           vt2 = ONE/sqrt(vt3)
            ! vectmp2 now contains 1/fij
    end if 

           vt4 = ZERO

           if (gb_kappa .ne. 0.d0) then
             vt3 = -ToGBFloat(gb_kappa / vt2)
             vt4 = exp(vt3)
           end if
#endif ! SOFTCORE_TI
           !if (.not. skipv(j)) then
           !    rinv = ONE/sqrt(r2) ! 1/rij
           !else
           !     rinv = ZERO 
           !end if
           ! Rewriting the above calculation to avid if-else statement
           !rinv = skipv(j) / sqrt(r2) ! 1/rij or 0
           ! Moved this computation ahead

!TODO 
!           if (onstep .and. r2 .gt. cut_inner2) then
!                  use_nrespa = 1.d0
!           else
                  use_nrespa = ZERO 
!           end if
           qiqj = qi * charge_lcl(j)
#ifdef SOFTCORE_TI
           de = 0.0d0 !because igb == 6 
           num_ti_atms_cntd_j = 0 
#endif ! SOFTCORE_TI


           alpb_factor = ToGBFloat(1 - alpb_zero)
           gb_kappa_factor = ToGBFloat(1 - gb_kappa_zero)

#ifndef SOFTCORE_TI
      if (igb /= 6) then
           expmkf = extdiel_inv * (ToGBFloat(gb_kappa_zero)  + (gb_kappa_factor * vt4 ))

           fgbk = gb_kappa_factor  * vt3 * expmkf
           fgbk = fgbk + gb_kappa_factor *  alpb_factor * &
                          (fgbk * one_arad_beta * &
                          ToGBFloat(-vt3 * gb_kappa_inv)) 

           dl = intdiel_inv - expmkf

           fgbi = vt2 ! 1.d0/fij
            ! Energy : egb partial sum of ij pair

           e = -qiqj * dl * (fgbi + alpb_factor * one_arad_beta)

#ifdef GB_ENERGY
#ifdef TIMODE
           egb = egb + e * ti_mode_skip(j)
#else
           egb = egb + e 
#endif
#endif

           temp4 = fgbi * fgbi * fgbi ! 1.d0/fij^3

             ! [here, and in the gas-phase part, "de" contains -(1/r)(dE/dr)]
           temp6 = -qiqj * temp4 * (dl + fgbk)
              ! -qiqj/fij^3*[1/Ein - e(-Kfij)/Eout) -kappa*fij*
              ! exp(-kappa*fij)(1 + fij*a*b/A ) /Eout]

           temp1 = vt1 ! exp(-rij^2/[4*ai*aj])
           de = temp6 * (ONE - QUARTER * temp1)
           temp5 = HALF * temp1 * temp6 * (ri * rj + QUARTER * r2)

#ifdef TIMODE
           sumi = sumi + ri * temp5 * ti_mode_skip(j)
           sumdeijda(j) = sumdeijda(j) + rj * temp5 * ti_mode_skip(j)
#else
           sumi = sumi + ri * temp5 
           sumdeijda(j) = sumdeijda(j) + rj * temp5 
#endif
      else
        de = 0.0d0
      end if !igb/=6
#endif ! SOFTCORE_TI

        ! Computing skipv here
        rinv = skipv(j) / sqrt(r2) ! 1/rij or 0
        ! skip exclusions for remaining terms:
             ! gas-phase Coulomb energy:
           r2inv = rinv * rinv
#ifdef SOFTCORE_TI 
          if (i_is_sc .or. j_is_sc) then
            eel = qiqj *intdiel_inv * sc_eel_denom_sqrt
            eelt = eelt + eel
            de = de - qiqj * sc_eel_denom + sc_eel_denom_sqrt &
                  ** (1 + sceeorder)
          else
            eel = intdiel_inv * qiqj * rinv
            eelt = eelt + eel 
            de = de + eel * r2inv 
          end if
#else
           eel = intdiel_inv * qiqj * rinv
#ifdef GB_ENERGY
#ifdef TIMODE
           eelt = eelt + eel * ti_mode_skip(j)
#else
           eelt = eelt + eel 
#endif
#endif ! GB_ENERGY
           de = de + eel * r2inv
#endif ! SOFTCORE_TI

          ! 6-12 potential:
            r6inv = r2inv * r2inv * r2inv
#ifdef SOFTCORE_TI
          ic = ico(iaci + iac(j))
          if (ic .gt. 0) then
            r6 = r2 * r2 * r2 
            if (i_is_sc .or. j_is_sc) then 
              f6 = 1.0d0/(sc_vdw + r6 * ti_sigma6(ic))
              f12 = f6 * f6
              evdw = evdw + ti_foureps(ic) * ( f12 - f6)
              de = de + ti_foureps(ic) * r2 * r2 * f12 * &
                  ti_sigma6(ic) * (12.0d0 * f6 - 6.0d0)
            else 
              f6 = gbl_cn2(ic) * r6     
              f12 = gbl_cn1(ic) * (r6 * r6)
              evdw = evdw + (f12 - f6)
              de = de + (12.d0 * f12 - 6.d0 * f6) * r2inv
            end if !either atom is sc
          end if  ! (ic .gt. 0)
#else
            f6 = data_cn2(kk1) * r6inv
            f12 =data_cn1(kk1)  * (r6inv * r6inv)
#ifdef GB_ENERGY
#ifdef TIMODE
            evdw = evdw + (f12*TWELVTH - f6*SIXTH) * ti_mode_skip(j)
#else
            evdw = evdw + (f12*TWELVTH - f6*SIXTH) 
#endif
#endif
            de = de + (f12 - f6) * r2inv
#endif ! SOFTCORE_TI


#if defined (HAS_10_12) && !defined(SOFTCORE_TI)
          ! The following could be commented out if the Cornell et al.
          ! force field was always used, since then all hbond terms are zero.
          ic   = ico(ntypes * (iac(i) - 1) + iac(j))
        if (ic .lt. 0 ) then
          ! 10-12 potential:
          r10inv = r2inv * r2inv * r2inv * r2inv * r2inv
          f10 = gbl_bsol(-ic) * r10inv f12 = gbl_asol(-ic) * r10inv * r2inv
#ifdef GB_ENERGY
#ifdef TIMODE
          evdw = evdw + (f12 - f10) * ti_mode_skip(j)
#else
          evdw = evdw + (f12 - f10) 
#endif
#endif
          de = de + (TWELVE * f12 - TEN * f10) * r2inv
        end if
#endif

             ! derivatives

          de = de * ( use_nrespa * nrespa + ( ONE - use_nrespa) * nrespai)

#ifdef TIMODE
           dedx = de * xij * ti_mode_skip(j)
           dedy = de * yij * ti_mode_skip(j)
           dedz = de * zij * ti_mode_skip(j)
#else
           dedx = de * xij 
           dedy = de * yij
           dedz = de * zij
#endif

           frcx_red(j,thread_id) = frcx_red(j,thread_id) - dedx
           frcy_red(j,thread_id) = frcy_red(j,thread_id) - dedy
           frcz_red(j,thread_id) = frcz_red(j,thread_id) - dedz
        
           dumx = dumx + dedx
           dumy = dumy + dedy
           dumz = dumz + dedz
       end do
       frcx_red(i,thread_id) = frcx_red(i,thread_id) + dumx
       frcy_red(i,thread_id) = frcy_red(i,thread_id) + dumy
       frcz_red(i,thread_id) = frcz_red(i,thread_id) + dumz
     
#ifndef SOFTCORE_TI
       sumdeijda(i) = sumdeijda(i) + sumi
#endif ! SOFTCORE_TI

       iexcl = iexcl + numex(i)
    end do
  end do  !  i = 1, max_i
!$omp end do
!$omp end parallel

!OpenMP region ends for Off-Diagonal

    call get_wall_time(wall_s, wall_u)
    loop_time = dble(wall_s) * 1000.d0 + dble(wall_u) / 1000.d0 - &
           strt_time_ms
#ifdef GBTimer
    print *, "Off Diagonal =", loop_time 
#endif
    offd_time(mytaskid + 1) = offd_time(mytaskid +1 ) + loop_time
    call get_wall_time(wall_s, wall_u)
    strt_time_ms = dble(wall_s) * 1000.d0 + dble(wall_u) / 1000.d0

    call update_gb_time(calc_gb_offdiag_timer)
        !--------------------------------------------------------------------------
        !
        ! Step 3:  Finally, do the reduction over the sumdeijda terms:, adding
        !          into the forces those terms that involve derivatives of
        !          the GB terms (including the diagonal or "self" terms) with
        !          respect to the effective radii.  This is done by computing
        !          the vector dai / dxj, and using the chain rule with the
        !          previously-computed sumdeijda vector.
        !
        !          Also compute a surface-area dependent term if gbsa=1. This
        !          perhaps should be moved to a new subroutine, but it relies
        !          on some other stuff going on in this loop, so we may(?) save
        !          time putting everything here and occasionally evaluating a
        !          conditional if (gbsa .eq. 1)
        !
        !          Do these terms only at "nrespa" multiple-time step intervals;
        !          (when igb=2 or 5, one may need to do this at every step)
        !
        !--------------------------------------------------------------------------

  if (igb /=6) then !GAS Phase Check
       if (onstep) then
           neibr_cnt = 0 ! we have no neighbors yet for LCPO gbsa
#ifdef MPI
        !first, collect all the sumdeijda terms:
           if (numtasks .gt. 1) then
              call mpi_allreduce(sumdeijda, red_buf, atm_cnt, mpi_double_precision, &
                        mpi_sum, pmemd_comm, err_code_mpi)
              sumdeijda(1:atm_cnt) = red_buf(1:atm_cnt)

           end if

        call update_gb_time(dist_gb_rad_timer)
#endif
        
        frespa = nrespa

! OpenMP region starts for Diagonal

        !$omp parallel default(private) &
        !$omp& shared(outer_i,start_atm_radii,end_atm_radii,atm_cnt, iac, crdx, crdy, crdz,fs_lcl,rborn_lcl,charge_lcl,reff, &
        !$omp& gb_kappa,sumdeijda,egb,frcx_red,frcy_red,frcz_red,frespa,igb,igb78,igb2578,&
        !$omp& rgbmax,rgbmaxpsmax2,gb_neckscale,neck_idx,offset,frc,&
        !$omp& nlist_r2, nlist_jj, nlist_count, max_nlist_count,&
        !$omp& cut2,cut_inner2,onstep,extdiel_inv,intdiel_inv,alpb_beta,alpb, &
#ifdef __INTEL_COMPILER
        !$omp& gb_alpha_arry,gb_beta_arry,psi, neckMaxPos,neckMaxVal, num_threads, &
#else
        !$omp& gb_alpha_arry,gb_beta_arry,psi,num_threads, &
#endif
        !$omp& gb_gamma_arry, natbel,one_arad_beta,max_i,ti_mode,ti_lst,ti_mask_piece, &
        !$omp& ntypes,esurf,gbsa,rgbmax2i,calc_nlist) 

        thread_id = omp_get_thread_num() + 1

#ifdef GB_ENERGY
        !$omp do reduction (+:egb) schedule(dynamic)
#else
        !$omp do schedule(dynamic)
#endif
        do i = start_atm_radii, end_atm_radii
#ifdef TIMODE
              if (ti_lst(ti_mask_piece,i) .ne. 0) cycle
#endif

           f_xi = ZERO
           f_yi = ZERO
           f_zi = ZERO 
           qi = charge_lcl(i)
           expmkf = exp(-gb_kappa * reff(i)) * extdiel_inv
           dl = intdiel_inv - expmkf
           qi2h = HALF * qi * qi
           qid2h = qi2h * dl

           if (alpb .eq. 0) then
#ifdef GB_ENERGY
             egb = egb - qid2h / reff(i)
#endif
             temp7 = -sumdeijda(i) + qid2h - gb_kappa * qi2h * expmkf * reff(i)
           else
#ifdef GB_ENERGY
             egb = egb - qid2h * (ONE/reff(i) + one_arad_beta)
#endif
             temp7 = -sumdeijda(i) + qid2h - gb_kappa * qi2h * expmkf * reff(i) * &
               (ONE + one_arad_beta * reff(i))
           end if

           xi = crdx(i)
           yi = crdy(i)
           zi = crdz(i)

           ri = rborn_lcl(i) - offset
           ri1i = ONE / ri

           if (igb2578) then

        ! new onufriev: we have to later scale values by a
        !               alpha,beta,gamma -dependent factor:

              ri = rborn_lcl(i) - offset
              thi = tanh((gb_alpha_arry(i) + gb_gamma_arry(i) * psi(i) * psi(i) - &
                                gb_beta_arry(i) * psi(i)) * psi(i))
              thi2 = (gb_alpha_arry(i) + THREE * gb_gamma_arry(i) * psi(i) * psi(i) - &
                        TWO * gb_beta_arry(i) * psi(i)) * (ONE - thi * thi) * ri / &
                     rborn_lcl(i)
           end if

             if (calc_nlist) then
                temp_icount = nlist_count(i) 
                nstart = (i - start_atm_radii) * max_nlist_count
                temp_jj(1:temp_icount) = nlist_jj(nstart + 1:nstart + temp_icount)
                temp_r2(1:temp_icount) = nlist_r2(nstart + 1:nstart + temp_icount)
             else
                   temp_icount = 0
                do j=1, atm_cnt
                   xij = xi - crdx(j)
                   yij = yi - crdy(j)
                   zij = zi - crdz(j)
                   r2 = xij * xij + yij * yij + zij * zij

                   if( j .eq. i ) cycle
                   if (r2 .gt. rgbmaxpsmax2) cycle
                   temp_icount = temp_icount + 1
                   temp_jj(temp_icount) = j
                   temp_r2(temp_icount) = r2
                end do
             end if

           vectmp1(1:temp_icount) = ONE/sqrt(temp_r2(1:temp_icount))

           do kod = 1, temp_icount, DIAG_BATCH_SIZE
             kk1 = 0

             do k = kod, MIN0(kod+DIAG_BATCH_SIZE,temp_icount+1) -1 
                j = temp_jj(k) 
                r2 = temp_r2(k)
                sj =  fs_lcl(j)

                dij1i = vectmp1(k)
                dij = r2 * dij1i
                sj2 = sj * sj

                if (dij .gt. FOUR * sj) cycle
                kk1 = kk1 + 1
                vectmp3(kk1) = dij + sj
                if (dij .gt. ri + sj) then
                  vectmp2(kk1) = r2 - sj2
                  vectmp4(kk1) = dij - sj
                else if (dij .gt. abs(ri - sj)) then
                  vectmp2(kk1) = dij + sj
                  vectmp4(kk1) = ri
                else if (ri .lt. sj) then
                  vectmp2(kk1) = r2 - sj2
                  vectmp4(kk1) = sj - dij
                else
                  vectmp2(kk1) = ONE 
                  vectmp4(kk1) = ONE 
                end if
             end do ! K Loop END

             vectmp2(1:kk1) = ONE/vectmp2(1:kk1)
             vectmp3(1:kk1) = ONE/vectmp3(1:kk1)
             vectmp4(1:kk1) = vectmp4(1:kk1) * vectmp3(1:kk1)
             vectmp4(1:kk1) = log(vectmp4(1:kk1))

             kk1 = 0
        
             do outer_k = kod, MIN0(kod+DIAG_BATCH_SIZE,temp_icount+1) -1 ,kBATCH_SIZE
                do k = outer_k, MIN0(outer_k+kBATCH_SIZE, temp_icount+1)-1
                   j  = temp_jj(k)
                   r2 = temp_r2(k)
                   j3 = 3 * j

                   xij = xi - crdx(j)
                   yij = yi - crdy(j)
                   zij = zi - crdz(j)

                   dij1i = vectmp1(k)
                   dij = r2 * dij1i
                   sj = fs_lcl(j)
                   if (dij .gt. rgbmax + sj) cycle
                   sj2 = sj * sj

        ! datmp will hold (1/r)(dai/dr):

                  dij2i = dij1i * dij1i
                  dij3i = dij2i * dij1i

                  if (dij .gt. rgbmax - sj) then
                     temp1 = ONE / (dij - sj)
                     datmp = EIGHTH * dij3i * ((r2 + sj2) * &
                        (temp1 * temp1 - rgbmax2i) - TWO * log(rgbmax * temp1))

                  else if (dij .gt. FOUR * sj) then
                     tmpsd = sj2 * dij2i
                     dumbo = te + tmpsd * (tf + tmpsd * (tg + tmpsd * (th + tmpsd * thh)))
                     datmp = tmpsd * sj * dij2i * dij2i * dumbo

                  else if (dij .gt. ri + sj) then
                     kk1 = kk1 + 1
                     datmp = vectmp2(kk1) * sj * (-HALF * dij2i + vectmp2(kk1)) + &
                     QUARTER * dij3i * vectmp4(kk1)

                  else if (dij .gt. abs(ri - sj)) then
                     kk1 = kk1 + 1
                     datmp = -QUARTER * (-HALF * (r2 - ri * ri + sj2) * &
                        dij3i * ri1i * ri1i + dij1i * vectmp2(kk1) * &
                        (vectmp2(kk1) - dij1i) - dij3i * vectmp4(kk1))

                  else if (ri .lt. sj) then
                     kk1 = kk1 + 1
                     datmp = -HALF * (sj * dij2i * vectmp2(kk1) - &
                        TWO * sj * vectmp2(kk1) * vectmp2(kk1) - &
                        HALF * dij3i * vectmp4(kk1))

                  else
                     kk1 = kk1 + 1
                    datmp = ZERO 
                  end if  ! (dij .gt. 4.d0 * sj)

                  if (igb78) then
                     if (dij .lt. rborn_lcl(i) + rborn_lcl(j) + gb_neckcut) then

        ! Derivative of neck with respect to dij is:
        !                     5
        !              9 mdist
        !   (2 mdist + --------) neckMaxVal gb_neckscale
        !                 5
        ! -(------------------------)
        !                        6
        !             2   3 mdist  2
        !   (1 + mdist  + --------)
        !                    10

                        mdist = dij - neckMaxPos(neck_idx(i), neck_idx(j))
                        mdist2 = mdist * mdist
                        mdist3 = mdist2 * mdist
                        mdist5 = mdist2 * mdist3
                        mdist6 = mdist3 * mdist3

                        temp1 = ONE + mdist2 + (THREE * TENTH) * mdist6
                        temp1 = temp1 * temp1 * dij

        ! (Note "+" means subtracting derivative, since above
                        !     expression has leading "-")

                        datmp = datmp + ((TWO * mdist + (NINE/FIVE) * mdist5) * &
                        neckMaxVal(neck_idx(i), neck_idx(j)) * &
                        gb_neckscale) / temp1

                    end if ! if (dij < rborn(i) +rborn(j) + gb_neckcut)
                  end if ! (igb .eq. 7 .or. igb .eq. 8)

               datmp = -datmp * frespa * temp7

               if (igb2578) &
                 datmp = datmp * thi2

                 f_x = xij * datmp
                 f_y = yij * datmp
                 f_z = zij * datmp
  
                 frcx_red(j,thread_id) = frcx_red(j,thread_id) + f_x
                 frcy_red(j,thread_id) = frcy_red(j,thread_id) + f_y
                 frcz_red(j,thread_id) = frcz_red(j,thread_id) + f_z
          
                 f_xi = f_xi - f_x
                 f_yi = f_yi - f_y
                 f_zi = f_zi - f_z
  
               end do
           end do  !  k = 1, icount
        end do        

        frcx_red(i,thread_id) = frcx_red(i,thread_id) + f_xi
        frcy_red(i,thread_id) = frcy_red(i,thread_id) + f_yi
        frcz_red(i,thread_id) = frcz_red(i,thread_id) + f_zi
        
      end do   ! end loop over atom i = mytaskid + 1, max_i, numtasks
      !$omp end do
      
      ! Putting back forces from reduction arrays to frc array
      !$omp do
      do i = 1, atm_cnt
        frc(3 * i - 2) =  sum(ToDBLE(frcx_red(i,1:num_threads)))
        frc(3 * i - 1) =  sum(ToDBLE(frcy_red(i,1:num_threads)))
        frc(3 * i)     =  sum(ToDBLE(frcz_red(i,1:num_threads)))
      end do
      !$omp end do
      !$omp end parallel

! OpenMP region ends for Diagonal

    call get_wall_time(wall_s, wall_u)
    loop_time = dble(wall_s) * 1000.d0 + dble(wall_u) / 1000.d0 - &
               strt_time_ms
#ifdef GBTimer
    print *, "Diagonal =", loop_time 
#endif
    diag_time(mytaskid + 1) = diag_time(mytaskid +1 ) + loop_time

      call update_gb_time(calc_gb_diag_timer)
      if (gbsa .eq. 1) then
         call gbsa_ene(crd, frc, esurf ,atm_cnt,jj,r2x,natbel)
      end if

      ! Define neighbor list ineighbor for calculating LCPO areas
      call update_gb_time(calc_gb_lcpo_timer)
   end if  !  onstep
  end if  ! if igb/=6 - GB gas phase check 

if (igb == 6) then
      ! Putting back forces from reduction arrays to frc array
      do i = 1, atm_cnt
        frc(3 * i - 2) =  sum(ToDBLE(frcx_red(i,1:num_threads)))
        frc(3 * i - 1) =  sum(ToDBLE(frcy_red(i,1:num_threads)))
        frc(3 * i)     =  sum(ToDBLE(frcz_red(i,1:num_threads)))
      end do
end if

return

