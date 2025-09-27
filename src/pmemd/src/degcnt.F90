#include "copyright.i"

!*******************************************************************************
!
! Module: degcnt_mod
!
! Description: Count solute and solvent degrees of freedom.
!              
!*******************************************************************************

module degcnt_mod

  implicit none

contains

!*******************************************************************************
!
! Subroutine:   degcnt (DEGree CouNT)
!
! Description:
!
! This routine determines how many degrees of freedom can be ascribed
! to the "solute", and how many to the "solvent". Unlike any previous
! version of AMBER, it works correctly with the belly option, with
! separate solute/solvent coupling, etc.
!
! Essentially, this routine loops through the atoms in a belly run and 
! determines how many of the moving atoms are in the "solute" (IBELSL), 
! and how many are in the "solvent" (IBELSV). 
! 
! It then determines the number of degrees of freedom in the solute and in
! the solvent which are lost to SHAKE and TORCON constraints.
!
! INPUT:
! -----
! IBELLY > 0 if belly run is being performed.
! NAT is the number of atoms in the system.
! IGRP(I) is > 0 if atom I is part of moving belly.
! NSOLUT is the number of the last solute atom.
!
!   MAINTENANCE NOTE: NSOLUT was a separate variable, but was always set to
!                     NATOM.  Hence we now just pass NATOM in to this procedure.
!                     I have not confirmed that this is correct, but it is
!                     effectively what was done with obfuscation in sander!
!
! NBONH, bondh() are the NBONH pointers to bonds to hydrogens.
! NBONA, bonda() are the NBONA pointers to bonds to non-hydrogens
! NTC is shake flag (1=no shake; 2,4=shake on bonds to H; 3=shake on all bonds)
! 
! Output:
! RNDFP = net number of degrees of freedom for the solute
! RNDFS = net number of degrees of freedom for the solvent
! 
! Removed from call interface:
!
! IBELSV = number of moving atoms in solvent, if IBELLY > 0.
!
! Author: David Pearlman
! Date: 5/91
!
!*******************************************************************************

subroutine degcnt(ibelly, nat, igrp, nsolut, bondh, bonda, ntc, rndfp, rndfs)

  use prmtop_dat_mod
  use ti_mod

  implicit none

! Formal arguments:

  integer               :: ibelly
  integer               :: nat
  integer               :: igrp(*)
  integer               :: nsolut
  type(bond_rec)        :: bondh(*)
  type(bond_rec)        :: bonda(*)
  integer               :: ntc
  double precision      :: rndfp
  double precision      :: rndfs

! Local variables:

  double precision      :: fract
  integer               :: i, j, m
  integer               :: ibelsl
  integer               :: ibelsv
  integer               :: ib, jb
  integer               :: im
  integer               :: inum
  double precision      :: rstssl
  double precision      :: rstssv
  ! For TI
  integer               :: iscsl(2) !number of dof in the TI region
  integer               :: iscsv(2)
  double precision      :: rstssl_sc(2)
  double precision      :: rstssv_sc(2)
  integer               :: iscnum
  integer               :: nsolut_tot
  integer               :: nsolvent_tot
  logical               :: is_sc_atm
! ibelsl = number of moving atoms in solute, if ibelly > 0.
! rstssl = number of degrees of freedom in the solute lost to shake.
! rstssv = number of degrees of freedom in the solvent lost to shake.
!         (rstssl and rstssv omit degrees of freedom already lost in the
!          belly case because the bond corresponds to two non-moving atoms)

! count up degrees of freedom in the belly, if any.

  ibelsl = 0
  ibelsv = 0
  if (ibelly .gt. 0) then
    do i = 1, nat
      if (igrp(i) .gt. 0) then
        if (i .le. nsolut) then
          ibelsl = ibelsl + 1
        else
          ibelsv = ibelsv + 1
        end if
      end if
    end do
  end if

! count up degrees of freedom that are associated with TI atoms

  iscsl(:) = 0
  iscsv(:) = 0
  if (ti_mode .ne. 0) then
    do i = 1, nat
      if (ti_lst(1,i) .ne. 0) then
        iscnum = 1
      else if (ti_lst(2,i) .ne. 0) then
        iscnum = 2
      else
        cycle
      end if
      if (i .le. nsolut) then
        iscsl(iscnum) = iscsl(iscnum) + 1
      else
        iscsv(iscnum) = iscsv(iscnum) + 1
      end if
    end do
    rstssl_sc(:) = 0.d0
    rstssv_sc(:) = 0.d0
    ti_region = 0
  end if

! Now, if shake is one, loop over the appropriate bond. Add up all bonds
! in the solute/solvent which cannot move because of shake. 
! For a belly run, do not count bonds which already cannot move because both
! atoms are part of the non-moving section.
! If a bond is 1/2 in solute, 1/2 in solvent, assign 1/2 to both counters.

  rstssl = 0.0d0
  rstssv = 0.0d0

! bonds to H

  if (ntc .ge. 2) then
    do i = 1, nbonh
      ib = bondh(i)%atm_i
      jb = bondh(i)%atm_j

      if (ti_mode .ne. 0) then
        if (ti_lst(1,ib)+ti_lst(1,jb) .ge. 1) then
          ti_region = 1
        else if (ti_lst(2,ib)+ti_lst(2,jb) .ge. 1) then    
          ti_region = 2
        else
          ti_region = 0
        end if
        is_sc_atm = ti_sc_lst(ib) .ne. 0 .or. ti_sc_lst(jb) .ne. 0
      end if

      if (ibelly .le. 0) then
         if (ib .le. nsolut .and. jb .le. nsolut) then
            rstssl = rstssl + 1.0d0
            if (ti_mode .ne. 0) then
               if (ti_region .eq. 0) then
                  rstssl_sc(1) = rstssl_sc(1) + 1.d0
                  rstssl_sc(2) = rstssl_sc(2) + 1.d0
               else
                  rstssl_sc(ti_region) = rstssl_sc(ti_region) + 1.d0      
                  if (is_sc_atm) then
                    ti_dof_shaked(ti_region) = ti_dof_shaked(ti_region) + 1                  
                  end if
               end if
            end if
         else if (ib .gt. nsolut .and. jb .gt. nsolut) then
            rstssv = rstssv + 1.0d0
            if (ti_mode .ne. 0) then
               if (ti_region .eq. 0) then
                  rstssv_sc(1) = rstssv_sc(1) + 1.d0
                  rstssv_sc(2) = rstssv_sc(2) + 1.d0
               else
                  rstssv_sc(ti_region) = rstssv_sc(ti_region) + 1.d0
                  if (is_sc_atm) then                  
                    ti_dof_shaked(ti_region) = ti_dof_shaked(ti_region) + 1
                  end if
               end if
            end if
         else
            rstssl = rstssl + 0.5d0
            rstssv = rstssv + 0.5d0
            if (ti_mode .ne. 0) then
               if (ti_region .eq. 0) then
                  rstssv_sc(1) = rstssv_sc(1) + 0.5d0
                  rstssv_sc(2) = rstssv_sc(2) + 0.5d0
                  rstssl_sc(1) = rstssl_sc(1) + 0.5d0
                  rstssl_sc(2) = rstssl_sc(2) + 0.5d0
               else
                  rstssv_sc(ti_region) = rstssv_sc(ti_region) + 0.5d0
                  rstssl_sc(ti_region) = rstssl_sc(ti_region) + 0.5d0         
                  if (is_sc_atm) then                  
                    ti_dof_shaked(ti_region) = ti_dof_shaked(ti_region) + 1
                  end if                 
               end if
            end if
         end if
      else if (igrp(ib) .gt. 0 .or. igrp(jb) .gt. 0) then
         if (ib .le. nsolut .and. jb .le. nsolut) then
            rstssl = rstssl + 1.0d0
         else if (ib .gt. nsolut .and. jb .gt. nsolut) then
            rstssv = rstssv + 1.0d0
         else
            rstssl = rstssl + 0.5d0
            rstssv = rstssv + 0.5d0
         end if
      end if
    end do
  end if

! bonds to heavy atoms

  if (ntc .eq. 3) then
    do i = 1, nbona
      ib = bonda(i)%atm_i
      jb = bonda(i)%atm_j
      if (ibelly .le. 0) then
        if (ib .le. nsolut .and. jb .le. nsolut) then
          rstssl = rstssl + 1.0d0
        else if (ib .gt. nsolut .and. jb .gt. nsolut) then
          rstssv = rstssv + 1.0d0
        else
          rstssl = rstssl + 0.5d0
          rstssv = rstssv + 0.5d0
        end if
      else if (igrp(ib) .gt. 0 .or. igrp(jb) .gt. 0) then
        if (ib .le. nsolut .and. jb .le. nsolut) then
          rstssl = rstssl + 1.0d0
        else if (ib .gt. nsolut .and. jb .gt. nsolut) then
          rstssv = rstssv + 1.0d0
        else
          rstssl = rstssl + 0.5d0
          rstssv = rstssv + 0.5d0
        end if
      end if
    end do
  end if

! Now determine RNDFP (net number of degrees of freedom for solute) and
!               RNDFS (net number of degrees of freedom for solvent).

  if (ibelly .le. 0) then
    rndfp = 3*nsolut - rstssl
    rndfs = 3*(nat - nsolut) - rstssv
  else
    rndfp = 3*ibelsl - rstssl
    rndfs = 3*ibelsv - rstssv
  end if

  if (ti_mode .ne. 0) then 
    nsolut_tot = nsolut - iscsl(1) - iscsl(2)
    nsolvent_tot = nat - nsolut - iscsv(1) - iscsv(2)
    rndfp = 3*nsolut_tot - rstssl
    rndfs = 3*nsolvent_tot - rstssv
    do i = 1, 2      
      nsolut_tot = nsolut - iscsl(2-i+1)
      nsolvent_tot = nat - nsolut - iscsv(2-i+1)
      ti_rndfp(i) = 3*nsolut_tot - rstssl_sc(i)
      ti_rndfs(i) = 3*nsolvent_tot - rstssv_sc(i)
    end do
  end if

  return

end subroutine degcnt

end module degcnt_mod
