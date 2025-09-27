#include "copyright.i"

!*******************************************************************************
!
! Module:  nmr_calls_mod
!
! Description:
!
!*******************************************************************************

module nmr_calls_mod

  use gbl_constants_mod

  implicit none

! Global data definitions.

  ! The following storage is per-process common; ie., it SHOULD be
  ! broadcast from the master to the other processes!

  integer, parameter    :: iredir_cnt = 8

  ! added jrsttyp so changed to 60
  integer, parameter    :: nmr_dat_int_cnt = iredir_cnt + 60

  ! nimprp is the number of improper torsional parameters.
  ! (nptra - nimprp is the number of regular torsional parameters).

  integer, private      :: intreq, irlreq, iredir(iredir_cnt), iuse, nstepl, &
                           nstep, j1nmr, j2nmr, j3nmr, j4nmr, jk2nmr, jk3nmr, &
                           jnmrat, jnmrst, jcomdf, jshrtd, jfntyp, jgravt, &
                           jwtnrg, jbave, jaave, jtave1, jtave2, jbav0, jaav0, &
                           jtav01, jtav02, jjcoef, jwtstp, jwttyp, jxpk, &
                           ichgwt, nmrnum, ishrtb, nave(3), nexact(3), &
                           ipower(3), navint(3), iavtyp(3), idmpav, maxrst, &
                           itimav, maxwt, maxgrp, jaltd, nimprp, jfxyz, &
                           joutxyz, jdxyz, jrsttyp

  common / nmr_dat_int /   intreq, irlreq, iredir, iuse, nstepl, &
                           nstep, j1nmr, j2nmr, j3nmr, j4nmr, jk2nmr, jk3nmr, &
                           jnmrat, jnmrst, jcomdf, jshrtd, jfntyp, jgravt, &
                           jwtnrg, jbave, jaave, jtave1, jtave2, jbav0, jaav0, &
                           jtav01, jtav02, jjcoef, jwtstp, jwttyp, jxpk, &
                           ichgwt, nmrnum, ishrtb, nave, nexact, &
                           ipower, navint, iavtyp, idmpav, maxrst, &
                           itimav, maxwt, maxgrp, jaltd, nimprp, jfxyz, &
                           joutxyz, jdxyz, jrsttyp

  ! JMS: add a flag to skip printing to the DUMPAVE file. Since this file is
  ! written during the nmrnrg call, simply calling nmrdcp will not prevent this
  ! data from being written even for frames where you don't want it to

  logical, public :: skip_print = .false.

  integer, parameter    :: nmr_dat_dbl_cnt = 14

  save  :: / nmr_dat_int /

  double precision, private     :: rwell, stpmlt, drjar, wtls(2), tauave(3), &
                                   taufin(3), ravein(3)

  common / nmr_dat_dbl /           rwell, stpmlt, drjar, wtls, tauave, &
                                   taufin, ravein

  save  :: / nmr_dat_dbl /

  ! nmr_work = the nmr reals work array.
  ! nmr_iwork = the nmr integers work array.

  double precision, private, allocatable, save  :: nmr_work(:)
  integer,          private, allocatable, save  :: nmr_iwork(:)

! redir and iredir contain information regarding listin, listout, readnmr,
! noesy, shifts, dumpave, pcshift and dipole respectively. If iredir(i) > 0,
! then that input/output has been redirected. The iredir array is currently
! broadcast (it is in the integer common block above).  I am not sure this is
! necessary.

  character(80), save   :: redir(iredir_cnt)

#ifdef MPI
  ! If eadev/ebdev are ever handled under mpi, they will need to be summed.
#else
  double precision, save        :: eadev = 0.d0
  double precision, save        :: ebdev = 0.d0
#endif

! RMS angle and bond energy deviations from ideal:


! Module private data that is NOT broadcast.  This is possible because it
! is either parameters or 0.d0-initialized stuff.

  integer, parameter, private           :: iscop1 = 33
  integer, parameter, private           :: iscop2 = 34
  integer, parameter, private           :: iscop3 = 35

  double precision, save, private       :: dvdis(2, 4) = 0.d0
  double precision, save, private       :: dvang(2, 4) = 0.d0
  double precision, save, private       :: dvtor(2, 4) = 0.d0
  double precision, save, private       :: eenmr(2, 3) = 0.d0

! Hide internal routines:

  private       alloc_nmr_mem, restal, nmrred, nmrnrg, modwt, nmrprt, ndvprt

contains

!*******************************************************************************
!
! Subroutine:  init_nmr_dat
!
! Description: <TBS>
!
!*******************************************************************************

subroutine init_nmr_dat(num_ints, num_reals)

  use mdin_ctrl_dat_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: ier
  integer               :: numgrp

  rewind(mdin)

  call restal(iscop1, mdin, ichgwt, nmrnum, itimav, mdout, ier, numgrp)


  maxrst = nmrnum + 1 + 5 * ier
  maxwt = ichgwt + 1 + 5 * ier
  maxgrp = numgrp + 1
  ! If imin .eq. 1 (minimization run), override any "time-averaged" requests.

  if (imin .eq. 1) itimav = 0

  intreq = 36 * maxrst + 7 * maxwt + 8 * maxgrp
  irlreq = 18 * maxrst + 12 * maxrst * itimav + 2 * maxwt

  call alloc_nmr_mem(num_ints, num_reals)

! Initialize nstepl and set the work/iwork storage partition.

  nstepl = 0
  nstep = 0

! work array:
  j1nmr = 1
  j2nmr = j1nmr + 2 * maxrst
  j3nmr = j2nmr + 2 * maxrst
  j4nmr = j3nmr + 2 * maxrst
  jk2nmr = j4nmr + 2 * maxrst
  jk3nmr = jk2nmr + 2 * maxrst
  jwtnrg = jk3nmr + 2 * maxrst
  jbave = jwtnrg + 2 * maxwt
  jbav0 = jbave + 2 * maxrst * itimav
  jaave = jbav0 + maxrst * itimav
  jaav0 = jaave + 2 * maxrst * itimav
  jtave1 = jaav0 + maxrst * itimav
  jtav01 = jtave1 + 2 * maxrst * itimav
  jtave2 = jtav01 + maxrst * itimav
  jtav02 = jtave2 + 2 * maxrst * itimav
  jjcoef = jtav02 + maxrst * itimav
  jdxyz  = jjcoef + 3 * maxrst
! jdxyz requires 3 * maxrst storage.

! iwork array:

  jnmrat = 1
  jnmrst = jnmrat + 16 * maxrst !16 nmrat
  jcomdf = jnmrst + 3 * maxrst
  jshrtd = jcomdf + 2 * maxgrp
  jfntyp = jshrtd + maxrst
  jgravt = jfntyp + maxrst
  jaltd  = jgravt + maxrst
  jwtstp = jaltd  + maxrst
  jwttyp = jwtstp + 3 * maxwt
  joutxyz = jwttyp + maxwt
  jfxyz  = joutxyz +  maxrst
  jxpk   = jfxyz  + 5 * maxrst
  jrsttyp = jxpk + 6 * maxrst !6 nmropt cases, currently only using 3
! jxpk requires 2 * maxrst of storage

  rwell = 1.0d0
  wtls(:) = 1.0d0
  tauave(:) = 1.0d+7

  return

end subroutine init_nmr_dat

!*******************************************************************************
!
! Subroutine:  alloc_nmr_mem
!
! Description: <TBS>
!
!*******************************************************************************

subroutine alloc_nmr_mem(num_ints, num_reals)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer                       :: alloc_failed

  allocate(nmr_work(irlreq), &
           nmr_iwork(intreq), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_reals = num_reals + size(nmr_work)

  num_ints = num_ints + size(nmr_iwork)

  nmr_work(:) = 0.d0
  nmr_iwork(:) = 0

  return

end subroutine alloc_nmr_mem

#ifdef MPI
!*******************************************************************************
!
! Subroutine:   bcast_nmr_dat
!
! Description:  When using the current AMBER/MPI implementation, only the master
!               reads and processes the input file.  This implies that only the
!               master sets anything in the nmrcal "READ" step or in the restal
!               subroutine call.  Therefore, we need to give the other
!               processors the data.
!*******************************************************************************

subroutine bcast_nmr_dat

  use parallel_dat_mod

  implicit none

! Local variables:

  integer               :: num_ints, num_reals  ! returned values discarded

  call mpi_bcast(intreq, nmr_dat_int_cnt, mpi_integer, 0, pmemd_comm, &
                 err_code_mpi)

  call mpi_bcast(rwell, nmr_dat_dbl_cnt, mpi_double_precision, 0, &
                 pmemd_comm, err_code_mpi)


  if (.not. master) then
    num_ints = 0
    num_reals = 0
    call alloc_nmr_mem(num_ints, num_reals)
  end if

  call mpi_bcast(nmr_work, irlreq, mpi_double_precision, &
                   0, pmemd_comm, err_code_mpi)

  call mpi_bcast(nmr_iwork, intreq, mpi_integer, &
                   0, pmemd_comm, err_code_mpi)

  return

end subroutine bcast_nmr_dat
#endif

!*******************************************************************************
!
! NOTE NOTE NOTE: The following subroutine header is being kept for historical
!                 purposes.  The nmrcal() subroutine and its various entry
!                 points have all been replaced by individual subroutine
!                 calls, with module private data taking the place of any
!                 shared components, and initialization activities moved into
!                 init_nmr_dat() and bcast_nmr_dat().
!
! Subroutine:  nmrcal
!
! Description:
!
! This is the top-level calling routine for the NMR/MD refinement suite of
! programs.
!
! This suite was written by
!
!            David A. Pearlman,
!            Department of Pharmaceutical Chemistry
!            University of California, San Francisco, 94143-0446
!
! Subroutine NMR calls. This routine makes the calls to the following
! NMR-MD refinement related routines:
!       restal : Does a cursory read of the NMR restraint files to
!                determine how much memory will be required. This
!                can be used in dynamic memory allocation schemes.
!       nmrred : Reads information defining NMR restraints & changes
!                in the relative weights of various energy terms.
!       nmrnrg : Calculates the energy/derivatives resulting from the
!                NMR restraints
!       modwt  : Modifies the relative weights of the various energy terms.
!
!       nmrprt : Prints the energy and average deviations for the restraints
!       ndvprt : Prints the deviation and energy contribution of each restraint.
!
! Calls (except to restal) are made thorugh a call to nmrcal
! with the appropriate value of caltyp (see below). Additionally,
! calls to nmrprt can be made through the entry-point nmrptx,
! and calls to ndvprt can be made through the entry point ndvptx.
!
! Calls to restal are always made through the entry point restlx.
!
! The entry point nmrdcp simply decrements the step-number counter and returns.
!
! Author: David A. Pearlman
! Date: 7/89
!
! INPUT:
!      caltyp: (Character*4 variable)
!              caltyp = "READ": Call to nmrred to read weight-change
!                               and restraint info is made.
!              caltyp = "CALC": Call to nmrnrg is made.
!              caltyp = "WEIT": Call to modwt is made.
!        f(i): Force array.
!atm_igraph(i): Name array. name(i) is the name of atom i. Packed as integers.
!gbl_labres(i): Residue name array. rsnam(i) is name of residue i. Packed as i.
!   rimass(i): The inverse mass of atom i.
!       temp0: (MD only) On the first call to modwt, should be the target
!              temperature as specified in the normal input. On this and
!              subsequent calls to modwt, the appropriate new target temperature
!              temp0 is returned.
!       tautp: Same behavior as temp0, but for the Berendsen scaling
!              parmameter.
! nmr_work(i): Real work array.
!              If averaging of restraints has not been requested,
!              work must be at least
!              12*maxrst + 12*maxrst*itimav + 2*maxwt elements long.
!              The work storage must not be modified between calls to nmrcal.
!
!              Note that  maxrst, itimav, maxwt, and maxgrp are all set
!              during the initialization call to restal. itimav is a flag
!              which is one if time-averaged restraints were requested.
!nmr_iwork(i): Integer work array.
!              iwork must be at least 12*maxrst+5*maxwt+2*maxgrp elements long.
!              The iwork storage must not be modified between calls to nmrcal.
!      maxrst: The maximum number of NMR restraints allowed.
!       maxwt: The maximum number of relative weight change definitions allowed.
!          in: The unit from which the NMR restraints/weight change definitions
!              will be read.
!        iout: The unit for informational prints.
!
!     The following are used only in modwt calls. They are modified if the
!     relative weights of the appropriate terms change:
!
!       rk(i): Force constant array for bonds.
!       tk(i): Force constant array for angles.
!       pk(i): Force constant array for torsions. Regular torsions, followed
!              by impropers.
!      cn1(i): "A" (r**12) coefficient array for vdw term.
!      cn2(i): "B" (r**6) coefficient array for vdw term.
!       ag(i): "A" (r**12) coefficient array for h-bond term.
!       bg(i): "B" (r**10) coefficient array for h-bond term.
!      numbnd: Number of bond force constants.
!      numang: Number of angle force constants.
!      numphi: Number of torsional force constants.
!      nimprp: Number of "improper" torsional force constants (not included
!              in numphi).
!         nhb: Number of h-bond parameters.
!      ncharg: Number of charge parameters. Equal to the number of atoms for
!              minimization and normal MD. Equal to 2*number_of_atoms for
!              perturbed systems (GIBBS).
!       natom: The number of atoms.
!      ntypes: The number of types of non-bonded (6-12) parameters.
!        nres: The number of residues in the system
!       rwell: The force constant used for non-bond interactions.
!
! The following are obtained from gbl_dyn_mem storage:
!
! atm_qterm(i): Charge on atom I.
!      ipres(i): Residue pointer array. ipres(i) points to the first atom
!                of residue i.
!
! The following are only passed in a call to restlx:
!
!      itotst: If > 0, then time-averaged requests will be honored.
!              Set to 0 to ignore time-averaged requests (e.g. for MIN).
!
! OUTPUT:
! ------
! (The following are set when a call to nmrnrg [caltyp="CALC"] is made).
!
!     enmr(3): (1) is the total energy contribution from distance restraints.
!              (2) is the total energy contribution from angle restraints.
!              (3) is the total energy contribution from torsion restraints.
!
! (devdis, devang, and devtor are not being used externally, so they are now
!  local to this routine)
!
!   devdis(4): (1) is the average deviation of the distance restraints from
!                  the "target" value. The average of r2 and r3 (see routine
!                  nmrnrg) is used as the "target distance"
!              (2) is the rms deviation of the distance restraints from the
!                  "target" value. The "target" value is defined as for (1).
!              (3) is the average deviation of the distance restraints from
!                  the "target". But in this case, a deviation of 0 is used
!                  when restraints fall between r2 and r3. For restraints
!                  outside this range, either r2 or r3 (whichever is closer)
!                  is used as the target value.
!              (4) is the rms deviation, calculated using the same deviations
!                  defined for (3).
!   devang(4): Same as devdis(4), but for angles. In radians.
!   devtor(4): Same as devdis(4), but for torsions. In radians.
!
! The following are set when a call to restal [entry point restlx] is made:
!
!     nmrnum:
! itmal(2,4):  Elements (i,1)->(i,3) contain the amount of storage allocated
!              for time-averaging of bonds, angles, and torsions (0 if
!              time-averaging not requested). itmal(1,4) contains the sum total
!              amount of space required. Generally calculated by routine
!              restal from the number of MD steps to be run, the numbers
!              of restraints, and the values of nave and nexact (see the
!              headers to routines nmrred and restal for more details).
!    maxrst :  The maximum allowable number of restraints.
!    itimav :  = 0 if no time-averaged restraints requested.
!              = 1 if time-averaged restraints requested.
!    maxwt  :  The maximum allowable number of weight-change cards.
!    maxgrp  : The maximum number of ranges of atoms which can be stored
!              as group information (only needed when using group center-of-mass
!              coordinates in distance restraints). This value is
!              computed in the initial call to restal.
!    intreq :  The amount of integer work array (iwork) storage required.
!    irlreq :  The amount of real work arary (work) storage required.
!
! Partitioning the work arrays:
!   The work storage is partitioned here. The following variables
!   are set to indicate the storage locations of the various arrays
!   in the work vectors:
!
!         j1nmr:  r1nmr(2,i)
!         j2nmr:  r2nmr(2,i)
!         j3nmr:  r3nmr(2,i)
!         j4nmr:  r4nmr(2,i)
!        jk2nmr:  rk2nmr(2,i)
!        jk3nmr:  rk3nmr(2,i)
!        jnmrat:  nmrat(4,i)
!        jnmrst:  nmrst(3,i)
!        jcomdf:  nmrcom(2,i)
!        jshrtd:  nmrshb(i)
!        jfntyp:  nmrfty(i)
!        jgravt:  igravt(i)
!         jaltd:  ialtdis(i)
!        jwtnrg:  wtnrg(2,i)
!         jbave:  bave(i)
!         jaave:  aave(i)
!        jtave1:  tave1(i)
!        jtave2:  tave2(i)
!         jbav0:  bave0(i)
!         jaav0:  aave0(i)
!        jtav01:  tave01(i)
!        jtav02:  tave02(i)
!        jjcoef:  ajcoef(3,i)
!        jwtstp:  iwtstp(3,i)
!        jwttyp:  iwttyp(i)
!          jxpk:  kxpk(i)
!
! Other Local Common:
!         nstepl: Local accumulator which keeps track of the number of
!         calls to calnrg. On each call, nstepl is incremented by 1.
!         The effective value of nstep used is nint(nstepl*stpmlt),
!         where stpmlt=1.0 if not reset by the user (resulting in
!         nstepl = nstep).
!         Since nstepl is incrmented after calnrg is called, modwt
!         should be called *before* calnrg on any cycle.
!
!         stpmlt modifies the effective increment in nstep per call.
!         By default, stpmlt is 1.0.
!
!         dvdis(2,4),dvang(2,4),dvtor(2,4) and eenmr(2,3) keep running totals
!         of devdis,devang,devtor and enmr, respectively, in the second
!         elements (e.g. dvdis(2,1), dvdis(2,2),...). The first elements
!         save the values of these quantities on the last call to nmrnrg.
!
! iscop1 and iscop2 are scratch unit numbers, used when input/output
! redirection is requested. iscop3 is a dedicated unit number, used for
! opening the file which contains periodic "dumps" of the restraint values.
!
!*******************************************************************************

!*******************************************************************************
!
! Subroutine:   nmr_read
!
! Description:  <TBS>
!*******************************************************************************

subroutine nmr_read(crd, in, iout)

  use prmtop_dat_mod
  use ti_mod

  implicit none

! Formal arguments:

  double precision      :: crd(*)
  integer               :: in
  integer               :: iout

! Local variables:

  integer               :: nstep0

  call nmrred(crd, atm_igraph, nmr_work(j1nmr), nmr_work(j2nmr), &
              nmr_work(j3nmr), nmr_work(j4nmr), nmr_work(jk2nmr), &
              nmr_work(jk3nmr), nmr_work(jjcoef), nmr_iwork(jnmrat), &
              nmr_iwork(jnmrst), nmr_iwork(jcomdf), nmr_iwork(jshrtd), &
              nmr_iwork(jfntyp), nmr_iwork(jgravt), nmr_iwork(jaltd), &
              nmrnum, nmr_work(jwtnrg), nmr_iwork(jwtstp), &
              nmr_iwork(jwttyp), nmr_iwork(jxpk), ichgwt, ishrtb, &
              nstep0, stpmlt, in, iout, iscop1, iscop2, &
              maxrst, maxwt, maxgrp, nave, nexact, navint, ipower, &
              ravein, taufin, iavtyp, idmpav, itimav, nmr_iwork(jfxyz), &
              nmr_iwork(joutxyz), nmr_work(jdxyz), nmr_iwork(jrsttyp))

  if (ti_mode .ne. 0) then
    call ti_check_res_int(nmrnum, nmr_iwork(jnmrat), nmr_iwork(jcomdf), iout)
  end if

  if (nstep0 .ne. 0) nstepl = nstep0
  if (nstep0 .ne. 0) nstep = nint(nstep0 * stpmlt)

#ifdef CUDA
  call gpu_nmr_set_nstep(nstep)
#endif

  return

end subroutine nmr_read

!*******************************************************************************
!
! Subroutine:   nmr_calc
!
! Description:  <TBS>
!*******************************************************************************

subroutine nmr_calc(crd, frc, enmr, iout)

  implicit none

! Formal arguments:

  double precision      :: crd(*)
  double precision      :: frc(*)
  double precision      :: enmr(3)
  integer               :: iout

! Local variables:

  double precision      :: devdis(4)
  double precision      :: devang(4)
  double precision      :: devtor(4)

  call nmrnrg(crd, frc, nmrnum, nstep, nmr_iwork(jnmrat), &
              nmr_iwork(jnmrst), nmr_work(j1nmr), nmr_work(j2nmr), &
              nmr_work(j3nmr), nmr_work(j4nmr), nmr_work(jk2nmr), &
              nmr_work(jk3nmr), nmr_iwork(jcomdf), nmr_iwork(jfntyp), &
              nmr_iwork(jgravt), nmr_iwork(jaltd), nmr_work(jbave), &
              nmr_work(jbav0), nmr_work(jaave), nmr_work(jaav0), &
              nmr_work(jtave1), nmr_work(jtav01), nmr_work(jtave2), &
              nmr_work(jtav02), nmr_work(jjcoef), nave, ipower, &
              tauave, ravein, navint, iavtyp, idmpav, iscop3, &
              enmr, devdis, devang, devtor, iout, nmr_iwork(jfxyz), &
              nmr_iwork(joutxyz), nmr_work(jdxyz), nmr_iwork(jrsttyp))

  dvdis(1, :) = devdis(:)
  dvang(1, :) = devang(:)
  dvtor(1, :) = devtor(:)
  eenmr(1, :) = enmr(:)
  dvdis(2, :) = dvdis(2, :) + devdis(:)
  dvang(2, :) = dvang(2, :) + devang(:)
  dvtor(2, :) = dvtor(2, :) + devtor(:)
  eenmr(2, :) = eenmr(2, :) + enmr(:)
  nstepl = nstepl + 1
  nstep = nint(nstepl*stpmlt)

#ifdef CUDA
  call gpu_nmr_set_nstep(nstep)
#endif

  return

end subroutine nmr_calc

!*******************************************************************************
!
! Subroutine:   nmr_weight
!
! Description:  See modwt for parameter descriptions.
!*******************************************************************************

subroutine nmr_weight(ncharg, crd, iout)

  use mdin_ctrl_dat_mod
  use nmr_lib_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer               :: ncharg
  double precision      :: crd(*)
  integer               :: iout

! Local variables:

  integer               :: numphi

  numphi = nptra - nimprp       ! "proper" dihedrals

! First check/reset short-range distance interactions flags:

  call nmrsht(crd, nmrnum, nstep, nmr_iwork(jnmrat), &
              nmr_iwork(jcomdf), ishrtb, nmr_iwork(jwtstp), &
              nmr_work(jwtnrg), nmr_work(jk2nmr), nmr_work(jk3nmr), &
              nmr_iwork(jnmrst), nmr_iwork(jshrtd), wtls, nmr_iwork(jfntyp), &
              nmr_iwork(jgravt), nmr_iwork(jaltd), nmr_work(jbave), &
              nmr_work(jbav0), nave, &
              ipower, tauave, navint, ravein, dt, iout, 1, nmr_iwork(jfxyz), &
              nmr_work(jdxyz), nmr_iwork(jrsttyp))

  call modwt(nmr_work(jwtnrg), nmr_iwork(jwtstp), nmr_iwork(jwttyp), ichgwt, &
             ishrtb, nstep, gbl_rk, gbl_tk, gbl_pk, &
             gbl_cn1, gbl_cn2, gbl_asol, gbl_bsol, &
             nmr_work(jk2nmr), nmr_work(jk3nmr), &
             nmr_iwork(jshrtd), nmr_iwork(jnmrst), numphi, &
             nimprp, nphb, ncharg, rwell, tauave, wtls, nmrnum, tautp, temp0)

  return

end subroutine nmr_weight

!*******************************************************************************
!
! Subroutine:   nmrptx
!
! Description:  The following subroutine allows a call to nmrprt without having
!               to specify lots of unnecessary stuff.
!*******************************************************************************

subroutine nmrptx(iout)

  implicit none

! Formal arguments:

  integer               :: iout

  call nmrprt(dvdis, dvang, dvtor, eenmr, nstep, iout)

  return

end subroutine nmrptx

!*******************************************************************************
!
! Subroutine:   ndvptx
!
! Description:  The following subroutine  allows a call to ndvprt without
!               having to specify lots of unnecessary stuff.
!*******************************************************************************

subroutine ndvptx(crd, frc, iout)

  use prmtop_dat_mod

  implicit none

! Formal arguments:

  double precision      :: crd(*)
  double precision      :: frc(*)
  integer               :: iout

  call ndvprt(crd, frc, atm_igraph, gbl_labres, nmrnum, nstep, &
              nmr_iwork(jnmrat), nmr_iwork(jxpk), nmr_iwork(jnmrst), &
              nmr_work(j1nmr), nmr_work(j2nmr), nmr_work(j3nmr), &
              nmr_work(j4nmr), nmr_work(jk2nmr), nmr_work(jk3nmr), &
              nmr_iwork(jcomdf), nmr_iwork(jfntyp), nmr_iwork(jgravt), &
              nmr_iwork(jaltd), nmr_work(jbave), nmr_work(jbav0), &
              nmr_work(jaave), nmr_work(jaav0), nmr_work(jtave1), &
              nmr_work(jtav01), nmr_work(jtave2), nmr_work(jtav02), &
              nmr_work(jjcoef), nave, ipower, taufin, ravein, navint, iscop1, &
              iout, nmr_iwork(jfxyz), nmr_iwork(joutxyz), nmr_work(jdxyz), &
              nmr_iwork(jrsttyp))

  return

end subroutine ndvptx

!*******************************************************************************
!
! Subroutine:   nmrdcp
!
! Description:  The following subroutine allows the local step number counter,
!               nstepl, to be decremented by one. This is necessary when a call
!               to force is made which does! not correspond to a dynamics step.
!               The accumuated deviation and energy counters are also
!               decremented by the distances on the last call.
!*******************************************************************************

subroutine nmrdcp

  implicit none

  nstepl = nstepl - 1
  nstep = nint(nstepl * stpmlt)

#ifdef CUDA
  call gpu_nmr_set_nstep(nstep)
#endif

  dvdis(2, :) = dvdis(2, :) - dvdis(1, :)
  dvang(2, :) = dvang(2, :) - dvang(1, :)
  dvtor(2, :) = dvtor(2, :) - dvtor(1, :)
  eenmr(2, :) = eenmr(2, :) - eenmr(1, :)

  return

end subroutine nmrdcp

!*******************************************************************************
!
! Subroutine:  restal (RESTraint Allocation)
!
! Description:
!
! This subroutine does a quick read of the weight-change/restraint
! cards used for NMR restraints, and returns the numbers of each (nmrnum
! and ichgwt, respectively).
!
! These quantities are required in setting up the dynamic memory
! allocation in routine locmem.
!
! It is expected that this routine will be called before any other of the
! NMR suite.
!
! Author: David A. Pearlman
! Date: 5/90
!
! INPUT:
!
!      iscop1 : Scratch units for informational reads/writes. Used when
!               redirection of input is requested.
!          in : Unit from which the restraints/weight changes will be read.
!        iout : Unit for informational prints.
!
! OUTPUT:
!
!      ichgwt : Number of change parameters read.
!      nmrnum : The total number of NMR restraints defined.
!      numgrp : The total number of NMR groups defined.
!      itimav : = 1 If any time-averaged restraints were requested.
!               = 0 otherwise.
!        ierr : Returned as 0 if no problems reading input. Returned as
!               1 if there were problems. Actual informational messages to
!               user about problems will be generated during the call to
!               nmrred.
!
!*******************************************************************************

subroutine restal(iscop1, in, ichgwt, nmrnum, itimav, iout, ierr, numgrp)

  use file_io_mod
  use file_io_dat_mod
  use nmr_lib_mod
  use prmtop_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer           iscop1
  integer           in
  integer           ichgwt
  integer           nmrnum
  integer           itimav
  integer           iout
  integer           ierr
  integer           numgrp

! Local variables:

! maxigr is maximum number of entries in a COM group
  integer           maxigr

  character(4)      atnam(8)
  character(80)     aline
  character(8)      flag(3)
  character(4), allocatable :: grnam1(:), grnam2(:), grnam3(:), grnam4(:), &
                               grnam5(:), grnam6(:), grnam7(:), grnam8(:), &
                               grnamarr(:,:), grnampass(:)

  character(10)     redirc(8)
  double precision  rjcoef(3)
  character(8)      type

  integer, allocatable :: igr1(:), igr2(:), igr3(:), igr4(:), igr5(:),igr6(:), &
                          igr7(:), igr8(:) !to hold 8 groups

  integer           outxyz, alloc_failed

  integer           istep1, istep2, iinc, imult

  double precision  value1, value2

  namelist /wt/     type, istep1, istep2, value1, value2, iinc, imult

  integer           iat(8)

  integer           nstep1, nstep2, irstyp, ninc, iresid, ir6, ifntyp, &
                    ifvari, ixpk, nxpk, ialtd, &
                    iin, ibeg, i, ifind, iformw, j, m, n, &
                    imat, istrt, ilengt, iothrd, ig

  double precision  r1, r2, r3, r4, rk2, rk3, r1a, r2a, r3a, r4a, rk2a, rk3a

  integer           fxyz(5)

  namelist /rst/        iat, nstep1, nstep2, irstyp, ninc, iresid, imult, &
                        atnam, igr1, igr2, igr3, igr4, igr5, igr6, igr7, igr8, &
                        grnam1, grnam2, grnam3, grnam4, grnam5, grnam6, grnam7, &
                        grnam8, r1, r2, r3, r4, rk2, rk3, r1a, r2a, r3a, r4a, &
                        rk2a, rk3a, ir6, ifntyp, ifvari, rjcoef, ixpk, nxpk, &
                        ialtd, fxyz, outxyz

  data redirc /'LISTIN    ' , 'LISTOUT   ' ,'DISANG    ', 'NOESY     ' , &
               'SHIFTS    ' , 'DUMPAVE   ', 'PCSHIFT   ' ,'DIPOLE    '/

  ! MAINTENANCE NOTE: NOESY, SHIFTS, PCSHIFT, and DIPOLE file input,
  !                   corresponding to iredir(4,5,7, and 8) are invalid in
  !                   pmemd.  In theory this input is only valid with
  !                   nmropt = 2, which is screened out, but in reality we
  !                   can get into restal if nmropt .gt. 0.  So we need to
  !                   filter these out here...

  data flag /'DISAVE  ' , 'ANGAVE  ' , 'TORAVE  '/

  ! Set maxigr value to the number of atoms set in the prmtop file
  maxigr = natom
  ! Allocate variables that depend on maxigr
  allocate(grnam1(maxigr), grnam2(maxigr), grnam3(maxigr), grnam4(maxigr), &
           grnam5(maxigr), grnam6(maxigr), grnam7(maxigr), grnam8(maxigr), &
           grnamarr(maxigr,8), grnampass(maxigr), &
           igr1(maxigr), igr2(maxigr), igr3(maxigr), igr4(maxigr), &
           igr5(maxigr), igr6(maxigr), igr7(maxigr), igr8(maxigr), &
           stat = alloc_failed)
  if (alloc_failed .ne. 0) call setup_alloc_error

  ichgwt = 0
  nmrnum = 0
  numgrp = 0
  itimav = 0
  iat(1:8)= 0
  ierr = 0
  iin = in
  ibeg = 1
  fxyz (1:5) = 1
  outxyz = 0

  do i = 1, 8
      iredir(i) = 0
  end do

! Look for a card starting with the string "&formwt". If one is found, then
! the subsequent cards, up until one with TYPE="END" is encountered, are
! the weight cards, as formatted input. If this string is not found,
! assume input will be in namelist format.

  iformw = 0
  call nmlsrc('formwt', in, ifind)
  if (ifind .eq. 1) then
    read(in, *)
    iformw = 1
  end if

! Read the weight changes:

9 do 10 i = ibeg, 99999

! Branch depending on whether formatted weight changes (iformw = 1) or
! namelist weight changes are being read:

    if (iformw .ne. 1) then

! Namelist-type read:

! Set the default values of the parameters which can be set by the namelist:

      istep1 = 0
      istep2 = 0
      value1 = 0.0
      value2 = 0.0
      iinc = 0
      imult = 0

! Look for "wt" namelist. If not found, stop with error:

      call nmlsrc('wt', iin, ifind)
      if (ifind .eq. 0) goto 1001

      read(iin, wt)

! Formatted-type weight card read:

    else
      read(iin, 9005, end = 1001, err = 1005) type, istep1, istep2, &
           value1, value2, iinc, imult
    end if

! End of weight card read. Now process given options:

! Do not count comment lines:

    if (type(1:1) .eq. '#') then
          goto 10

! If itype="END", then goto redirection reading section:

    else if (type(1:3) .eq. 'END') then
          goto 20
    else

! Check for cards requesting time-averaging:

      do j = 1, 3
        if (type .eq. flag(j)) then
          itimav = 1
        end if
      end do

    end if

    ichgwt = ichgwt + 1

10 continue

! Now read the redirection cards, if any. The first non-blank line encountered
! is assumed to be the start of the redirection section. If it is not a known
! option, assume no redirection cards where given.

20 read(iin, 9004, end = 25) aline

! Skip blank lines:

    do j = 1, len(aline)
      if (aline(j:j) .ne. ' ') goto 29
    end do

    goto 20

29  do m = 1, 8  ! 8 redirection flag, only 3 used
      do n = 10, 1, - 1
          if (redirc(m)(n:n) .ne. ' ') goto 24
      end do

24 continue

      call chklin(aline, redirc(m), n, imat, istrt, ilengt, iout)

      if (imat .ne. 0) then
        redir(m) = aline(istrt:istrt + ilengt - 1)
        iredir(m) = ilengt
        goto 20
      end if
    end do

    backspace(iin)

25 continue

! Check to be sure we did not read a card for an nmr option that has been
! dropped in pmemd.  These include NOESY, SHIFTS, PCSHIFT and DIPOLE.
! Note that the indices checked are hardwired - what a kludge!

    if (iredir(4) .ne. 0 .or. iredir(5) .ne. 0 .or. iredir(7) .ne. 0 .or. &
        iredir(8) .ne. 0) then
      write(mdout,*) &
        'No PMEMD support for NOESY, SHIFTS, PCSHIFT or DIPOLE options!'
      write(mdout,*)'Please use SANDER instead.'

      call mexit(6, 1)
    end if

! Now read the restraints:

! If restraint input has been redirected, open the appropriate file:

    iothrd = 0
    if (iredir(3) .ne. 0) then
      call opnmrg(redir(3)(1:iredir(3)), iscop1, 2, iout, ierr)
      if (ierr .eq. 1) then
        call mexit(mdout, 1)
      end if
      iin = iscop1
      iothrd = iin
    end if

! Count the restraints:
    ifind = 1
    do while (ifind .ne. 0)

! Look for "rst" namelist. If not found, assume we are done reading them:

      call nmlsrc('rst', iin, ifind)
      if (ifind .eq. 0) goto 100

! Read the rst namelist.
      iat(1:8)= 0   ! need to reset these before every namelist read
      igr1(:) = 0
      igr2(:) = 0
      igr3(:) = 0 !add 3 to 8; 3 and 4 for COM angle, torsion grp
      igr4(:) = 0 ! 5 to 8 for future COM uses
      igr5(:) = 0
      igr6(:) = 0
      igr7(:) = 0
      igr8(:) = 0
      fxyz(:) = 1
      outxyz = 0
      read(iin, rst, end = 100)

! If iat(5)!=0, this is a restraint unknown to pmemd. Warn and bail.
! for iat(5), iat(6), iat(7), iat(8) make changes here

      if (iat(5) .ne. 0) then
        write(mdout, '(a)') 'PMEMD only supports distance, angle, and dihedral &
                            &restraints!'
        write(mdout, '(a)') 'Please use SANDER instead!'
        call mexit(mdout, 1)
      end if

#ifdef CUDA
      do n = 1, 8
        if ( iat(n) < 0 ) then
          if ( n > 4 ) then
            write(mdout, '(a)')   'CUDA (GPU): Implementation does not currently support iat(n) > 4'
                write(mdout, '(a,i0,a)') '            Require iat(',n,') >= 0.'
                call mexit(mdout, 1)
          end if
        endif
     enddo

#endif

! Set fxyz(4) to be single direction or not.  0 is multidirection.  1 is single
! in x, 2 is single in y, 3 is single in z.  Set fxyz(5) to 0.


      fxyz(5) = 0

! Restricting restraints to the x-, y-, or z-direction is only allowed for
! distance and COM distance restraints and must be between 0.0 and 1.0.
! Catch errors relating to this.
     if (fxyz(1) .ne. 1 .or. fxyz(2) .ne. 1 .or. fxyz(3) .ne. 1) then
       if (iat(3) .ne. 0) then
            write(mdout, '(a)') "Restricting restraints in the x-, y-, and/or z-direction(s) using fxyz &
                                & is only allowed for distance and COM distance restraints"
            call mexit(mdout,1)
       end if
       if (fxyz(1) .ne. 1 .and. fxyz(1) .ne. 0) then
            write(mdout, '(a)') "Distance Restraint Error: fxyz x-dimension must be 0 or 1"
            call mexit(mdout,1)
       end if
       if (fxyz(2) .ne. 1 .and. fxyz(2) .ne. 0) then
            write(mdout, '(a)') "Distance Restraint Error: fxyz y-dimension must be 0 or 1"
            call mexit(mdout,1)
       end if
       if (fxyz(3) .ne. 1 .and. fxyz(3) .ne. 0) then
            write(mdout, '(a)') "Distance Restraint Error: fxyz z-dimension must be 0 or 1"
            call mexit(mdout,1)
       end if
     end if

     if (outxyz .ne. 0 .and. iat(3) .ne. 0) then
         write(mdout, '(a)') "XYZ Component Output: outxyz only available for distance and COM distance restraints"
         call mexit(mdout,1)
     end if

! Now we have support for COM restraints for angle and torsions
! erase comments below or ignore
! The only COM restraints allowed are for distances. Catch any attempted COM
! restraints that are not distances and bail out with an error

      if (iat(1) .lt. 0 .or. iat(2) .lt. 0 .or. iat(3) .lt. 0 .or. iat(4) .lt. 0) then
        if (iat(5) .ne. 0) then
          write(mdout, '(a)') 'PMEMD only supports atom groups in distance &
                              &angle, and torsion restraints only!'
          write(mdout, '(a,4i5)') 'iat1->4 are: ', iat(1),iat(2),iat(3),iat(4)
          write(mdout, '(a)') 'Please use SANDER instead!'
          call mexit(mdout, 1)
        end if
      end if

! If iat(1)=0, assume this a blank line & the last restraint has been read.

      if (iat(1) .lt. 0) then
        numgrp = numgrp + 1
        do ig = 2, maxigr
          if (igr1(ig) .eq. 0) goto 32
#ifdef CUDA
          numgrp = numgrp + 1
#else
          if (igr1(ig) .ne. igr1(ig - 1) + 1) numgrp = numgrp + 1
#endif
        end do
      end if

32 continue

      if (iat(2) .lt. 0) then
        numgrp = numgrp + 1
        do ig = 2, maxigr
          if (igr2(ig) .eq. 0) goto 34
#ifdef CUDA
          numgrp=numgrp + 1
#else
          if (igr2(ig) .ne. igr2(ig - 1) + 1) numgrp = numgrp + 1
#endif
        end do
      end if

34 continue

      if (iat(3) .lt. 0) then          ! COM angle
        numgrp = numgrp + 1
        do ig = 2, maxigr
          if (igr3(ig) .eq. 0) goto 36
          if (igr3(ig) .ne. igr3(ig - 1) + 1) numgrp = numgrp + 1
        end do
      end if

36 continue

      if (iat(4) .lt. 0) then          !  COM torsion
        numgrp = numgrp + 1
        do ig = 2, maxigr
          if (igr4(ig) .eq.0) goto 38
          if (igr4(ig) .ne. igr4(ig - 1) + 1) numgrp = numgrp + 1
        end do
      end if

38 continue

      nmrnum = nmrnum + 1

      if (iat(1) .eq. 0) goto 100

    end do ! do while (ifind .ne. 0)

100 continue

  if (iothrd .ne. 0) call opnmrg(redir(3), iin, 3, iout, ierr)

  return

! Errors:
! Do not report these here, just return with ierr=1 flag. The actual reporting
! and program termination will occur when routine nmrred is called.

! End-of-file read errors:

1001 continue

  ierr = 1
  return

! Other read errors (probably format mismatch):
! (If "#" is in first column, this is a comment line):

1005 continue

  backspace(iin)
  read(iin, 9004) aline
  if (aline(1:1) .eq. '#') then
    ibeg = i + 1
    goto 9
  else
    ierr = 1
    return
  end if

9004 format(a80)
9005 format(a8, 2i7, 2f12.6, 2i7)

end subroutine restal

!*******************************************************************************
!
! Subroutine:   nmrred (NMR REaD)
!
! Description:
!
! This subroutine reads: 1) The NMR restraints;
!                        2) Flags describing changes to the relative
!                           weights of the energy terms.
!
! Author: David A. Pearlman
! Date: 7/89
!
! Information is read in three sections.
!
! SECTION 1 -- "WEIGHT CHANGE" cards:
!
! Information can either be read formatted, in the format described below,
! or as part of the "wt" namelist.  To signal that formatted input is
! desired, you must preceed the weight-change card section by the flag
! "&formwt".  This flag can be in any column, but must be the first not
! blank set of characters on the line.
!
! Namelist input follows standard namelist rules...
!
! Control lines of the format:
!
!        type , istep1 , istep2 , value1 , value2 , iinc , imult
!
!            format(a8,2i5,2f12.6,2i5)
!
! After "END" is read in the previous section, input/output redirection
! information can be read as described here.  The inclusion of cards for
! redirection is optional.  The format of the redirection cards is:
!
!         TYPE = filename
!
! where TYPE is any valid redirection keyword (see below), and filename
! is any character string. The equals sign ("=") is required.
!
! The input/output redirection cards (if any) are followed by the restraints.
!
! INPUT VARIABLES:
! ---------------
!
!     name(i) : The name corresponding to atom i.
!   rimass(i) : The inverse mass of atom i.
!          in : Unit from which the restraints/weight changes will be read.
!        iout : Unit for informational prints.
!      iscop1
!      iscop2 : Scratch units for informational reads/writes. Used when
!               redirection of input/output is requested.
!       natom : The number of atoms.
!        nres : The number of residues.
!      maxrst : The maximum number of restraints for which array space is avail.
!       maxwt : The maximum number of weight changes
!      maxgrp : The maximum number of atom range definitions relating to
!               group definitions (for center-of-mass distance restraints only)
!               for which array space is available.
!
! Input variables from global memory:
!
!gbl_res_atms(i) : Residue pointer array. gbl_res_atms(i) points to the first
!                  atom of residue i.
!
! OUTPUT:
!
! (re: restraints)
!
!  r1nmr(2,i) : The slope and intercept, respectively, defining the dependency
!               of r1 on step number for restraint i.
!  r2nmr(2,i) : The slope and intercept for r2.
!  r3nmr(2,i) : The slope and intercept for r3.
!  r4nmr(2,i) : The slope and intercept for r4.
! rk2nmr(2,i) : The slope and intercept for k2.
! rk3nmr(2,i) : The slope and intercept for k3.
! ajcoef(3,i) : For torsional restraints, if rjcoef(3,i) .ne. 0, then
!               a J-coupling restraint based on the Karplus relationship
!               will be used.
!  nmrat(16,i) : The atoms defining restraint i. Stored as (3*(iat-1)), where
!               iat is the absolute atom number.
!  RESTTYPE(I) :  An integer identifying the type of restraint
!                           = 1, this is a distance restraint
!                           = 2, this is an angle restraint
!                           = 3, this is a torsion or j-coupling restraint
!                           = 4, this is a plane-point restraint
!                           = 5, this is a plane-plane restraint
!                           = 6, this is a generalized distance coordinate restraint
!                               There are plans to implement generalized angle and torsion
!                               restraints, but not for this version of AMBER

!  nmrst(3,i) : The beginning (1) and ending (2) steps over which this
!               restraint is to be applied. nmrst(3,i) gives the increment
!               between changes in the target values. If nmrst(3,i)=0,
!               the target values are varied continuously.
! ajcoef(3,i) : The coefficients to be used in the Karplus equation for
!               any torsional restraints to be imposed using J-coupling.
!   nmrshb(i) : Flag for each restraint, which keeps track of whether that
!               restraint is part of a defined "short-range" interactions
!               group. nmrshb(i) = 0 (not short-range or short range not
!               defined); = 1 (short-range based on residue sequence);
!               =2 (short-range based on distance criteria); =3
!               (short-range based on both residue sequence and distance
!               criteria).
!   nmrfty(i) : Flag for each restraint which indicates what functional
!               form is used for the restraint. See description of ifntyp
!               above.
!      nmrnum : The total number of NMR restraints defined.
!
! (re: weight/temperature changes)
!
!  wtnrg(2,i) : The slope and intercept, respectively, defining the dependency
!               of weight/temperature/etc. Change i on step number.
! iwtstp(3,i) : The beginning and ending steps over which change i is to be
!               applied (elements 1 and 2). If iwtstp(3,i)>0, the change
!               in target values is done as a step function every iwtstp(3,i)
!               steps. Otherwise, the change is done continuously.
!   iwttyp(i) : The type(s) of variable(s) being changed.
!        kxpk : Peak identifier.
!      ichgwt : Number of change parameters read.
!      ishrtb : If >0, then short-range interactions have been defined
!               (using the SHORT command card). The definitions of the
!               interactions are stored in iwtstp(,ishrtb) and
!               wtnrg(,ishrtb)>
!
! (re: other)
!
!      nstep0 : If an nstep0 card is read, the initial value of the nstep0
!               counter will be set to the given value. Otherwise, nstep0
!               is returned with a value of 0.
!      stpmlt : If a stpmlt card is read, the step counter multiplier will
!               be set to the given value. Otherwise, stpmlt is returned
!               with a value of 1.0.
!     nave(3) : For time-averaged restraints. If time averaged restraints
!               have been requested, nave(i) > 0.
!               elements 1->3 correspond to bonds, angles, torsions.
!   nexact(3) : Not currently used.
!   navint(3) : If > 0, then averaged restraints are updated every
!               navint(i) steps.
!   ipower(3) : For time-averaged restraints. Gives the exponent used in
!               developing the time average.
!   ravein(3) : For time-averaged values, the value for the integral
!               on the first step will be the instantaneous value+ravein(i)
!               (i=1,3 for bonds, angles, torsions).
!   taufin(3) : The value of the exponential decay constant used in calculating
!               the final time-averaged values reported by ndvprt.
!   iavtyp(3) : Determines how forces will be calculated when time-averaged
!               retraints are used. iavtyp(i) = 1 --> calculate forces
!               according to the standard energy expression. iavtyp(i) = 2 -->
!               calculate forces as dE/d(r_ave) dr(t)/d(x).
!   idmpav    : If > 0, then "dump" the values of the restraints to a file
!               every idmpav steps. The dump file is specified by a
!               dumpave redirection card.
!   itimav    : If = 0, then requests for time-averaged restraints will be
!               ignored.
!
!*******************************************************************************

subroutine nmrred(crd, name, r1nmr, r2nmr, r3nmr, r4nmr, rk2nmr, rk3nmr, &
                  ajcoef, nmrat, nmrst, nmrcom, nmrshb, nmrfty, igravt, &
                  ialtdis, nmrnum, wtnrg, iwtstp, iwttyp, kxpk, ichgwt, &
                  ishrtb, nstep0, stpmlt, in, iout, iscop1, iscop2, maxrst, &
                  maxwt, maxgrp, nave, nexact, navint, ipower, ravein, &
                  taufin, iavtyp, idmpav, itimav, ifxyz, ioutxyz, idxyz, &
                  resttype)

  use file_io_mod
  use file_io_dat_mod
  use nmr_lib_mod
  use prmtop_dat_mod
  use pmemd_lib_mod
  use mdin_ctrl_dat_mod, only : jar, nstlim

  implicit none

! Formal arguments:

  double precision  crd(*)
  character(4)      name(*)
  double precision  r1nmr(2, *)
  double precision  r2nmr(2, *)
  double precision  r3nmr(2, *)
  double precision  r4nmr(2, *)
  double precision  rk2nmr(2, *)
  double precision  rk3nmr(2, *)
  double precision  ajcoef(3, *)
  integer           nmrat(16, *)
  integer           nmrst(3, *)
  integer           nmrcom(2, *)
  integer           nmrshb(*)
  integer           nmrfty(*)
  integer           igravt(*)
  integer           ialtdis(*)
  integer           nmrnum
  double precision  wtnrg(2, *)
  integer           iwtstp(3, *)
  integer           iwttyp(*)
  integer           kxpk(2, *)
  integer           ichgwt
  integer           ishrtb
  integer           nstep0
  double precision  stpmlt
  integer           in
  integer           iout
  integer           iscop1
  integer           iscop2
  integer           maxrst
  integer           maxwt
  integer           maxgrp
  integer           nave(3)
  integer           nexact(3)
  integer           navint(3)
  integer           ipower(3)
  double precision  ravein(3)
  double precision  taufin(3)
  integer           iavtyp(3)
  integer           idmpav
  integer           itimav
  integer           ii
  integer           resttype(*) !added for case implementation
  integer           resttypetmp !added for case implementation

  double precision  idxyz(*)
  integer           ioutxyz(*)
  integer           ifxyz(*)

! Local variables:

  character(80)     aline

  ! maxigr is maximum number of entries in a COM group
  integer ixpk, nxpk, maxigr, alloc_failed

! Added 9/2007 by Matthew Seetin from sander to enable natural language input of restraints
  character(4)      atnam(8)
  character(4), allocatable :: grnam1(:), grnam2(:), grnam3(:), grnam4(:), & ! grnam3 and grnam4 for COM angle and torsions
                               grnam5(:), grnam6(:), grnam7(:), grnam8(:), & ! groups 5-8 for future-proofing if anyone wants to add more restraints
                               grnamarr(:,:), grnampass(:)
  character(4)      nametmp(8) ! holds the COM or iat info
  logical jcoupl
  logical usecom    ! add for com
  ! iflag is the number of "TYPE" flags recognized (make sure "END" is last):              |  ! iflag is the number o
! iflag is the number of "TYPE" flags recognized (make sure "END" is last):

  integer           iflag

  parameter (iflag = 27)

! zerner is a value near zero used for weights set to 0:

  double precision  zerner

  parameter (zerner = 1.0d-7)

  character(8)      type, flag(iflag)

  character(10)     redirc(8)

  double precision  rjcoef(3)

  integer, allocatable :: igr1(:), igr2(:), igr3(:), igr4(:), igr5(:), igr6(:), &
                          igr7(:), igr8(:), igrarr(:,:), igrpass(:)

  integer           idumm(3)
  double precision  xcom(12), dfr(12), rmstot(4), wtls(2)
  character(4)      iatnm(8)

! The following is all related to the use of namelist reads:

  integer           istep1, istep2, iinc, imult

  double precision  value1, value2

  namelist /wt/ type, istep1, istep2, value1, value2, iinc, imult

  integer           iat(8), nstep1, nstep2, irstyp, ninc, iresid, &
                    ir6, ifntyp, ifvari, ialtd, outxyz

  double precision  r1, r2, r3, r4, rk2, rk3, r1a, r2a, r3a, &
                    r4a, rk2a, rk3a

  integer           fxyz(5)


  namelist /rst/ iat, nstep1, nstep2, irstyp, ninc, iresid, imult, &
                 atnam, igr1, igr2, igr3, igr4, igr5, igr6, igr7, igr8, &
                 grnam1, grnam2, grnam3, grnam4, grnam5, grnam6, grnam7,  &
                 grnam8, r1, r2, r3, r4, rk2, rk3, r1a, r2a, r3a, r4a, rk2a, &
                 rk3a, ir6, ifntyp, ifvari, rjcoef, ixpk, nxpk, ialtd, fxyz, &
                 outxyz

  ! NOTE - The following flag data array is referenced by number for
  !        filtering purposes, which means that one had better not move
  !        anything without fixing things up!!!
  ! Options not supported in PMEMD: CUT (17)

  data flag /'BOND    ' ,  'ANGLE   ' ,  'TORSION ' ,  'VDW     ',  &
             'HB      ' ,  'ELEC    ' ,  'NB      ' ,  'ATTRACT ',  &
             'REPULSE ' ,  'INTERN  ' ,  'ALL     ' ,  'REST    ',  &
             'TEMP0   ' ,  'RSTAR   ' ,  'IMPROP  ' ,  'SOFTR   ',  &
             'CUT     ' ,  'SHORT   ' ,  'RESTS   ' ,  'RESTL   ',  &
             'NOESY   ' ,  'SHIFTS  ' ,  'TAUTP   ' ,  'DISAVE  ',  &
             'ANGAVE  ' ,  'TORAVE  ' ,  'END     '/

  data redirc /'LISTIN    ' ,  'LISTOUT   ' ,  'DISANG    ',  &
               'NOESY     ' ,  'SHIFTS    ' ,  'DUMPAVE   ',  &
               'PCSHIFT   ' ,  'DIPOLE    '/

  ! MAINTENANCE NOTE: NOESY, SHIFTS, PCSHIFT, and DIPOLE file input,
  !                   corresponding to iredir(4,5,7, and 8) are invalid in
  !                   pmemd.  In theory this input is only valid with
  !                   nmropt = 2, which is screened out, but in reality we
  !                   can get into restal if nmropt .gt. 0.  So we need to
  !                   filter these out here...

  integer           lastpk, ibeg, ihits, i, ifind, iin, iformw, j, jj, &
                    jdx, ntu, m, n, imat, istrt, ilengt, iothrd, ierr, kk

  double precision  slope, rint, rntu, step1, step2, frc_dummy(1), &
                    e, dumm, convrt, rprt, apmean, dev2, rcurr, dumm_array(3)

  ! Set maxigr value to the number of atoms set in the prmtop file
  maxigr = natom
  ! Allocate variables that depend on maxigr
  allocate(grnam1(maxigr), grnam2(maxigr), grnam3(maxigr), grnam4(maxigr), &
           grnam5(maxigr), grnam6(maxigr), grnam7(maxigr), grnam8(maxigr), &
           grnamarr(maxigr,8), grnampass(maxigr), &
           igr1(maxigr), igr2(maxigr), igr3(maxigr), igr4(maxigr), &
           igr5(maxigr), igr6(maxigr), igr7(maxigr), igr8(maxigr), &
           igrarr(maxigr,8), igrpass(maxigr), &
           stat = alloc_failed)
  if (alloc_failed .ne. 0) call setup_alloc_error

  nstep0 = 0
  stpmlt = 1.0d0
  lastpk = 1
  ishrtb = 0
  idmpav = 0
  ibeg = 1
  iat(1:8) = 0
  wtls(1) = 1.0d0
  wtls(2) = 1.0d0
  ravein(1) = 0.0d0
  ravein(2) = 0.0d0
  ravein(3) = 0.0d0
  ihits = 0
  fxyz(:) = 1
  outxyz = 0

  do i = 1, 3
    iredir(i) = 0
    nave(i) = 0
    nexact(i) = 0
    taufin(i) = 1.0d+7
    navint(i) = 1
    iavtyp(i) = 1
  end do

  iredir(4) = 0
  iredir(5) = 0
  iredir(6) = 0
  iredir(7) = 0
  iredir(8) = 0

  write(iout, 9000)

! Look for a card starting with the string "&formwt". If one is found, then
! the subsequent cards, up until one with TYPE="END" is encountered, are
! the weight cards, as formatted input. If this string is not found,
! assume input will be in namelist format.

  iformw = 0
  call nmlsrc('formwt', in, ifind)
  if (ifind .eq. 1) then
    read(in, *)
    iformw = 1
  end if

! Read the weight changes:

  iin = in

9 continue

  do 10 i = ibeg, 99999

! Branch depending on whether formatted weight changes (iformw = 1) or
! namelist weight changes are being read:

11 continue

    if (iformw .ne. 1) then

! NAMELIST-type read:

! Set the default values of the parameters which can be set by the namelist:

      istep1 = 0
      istep2 = 0
      value1 = 0.0
      value2 = 0.0
      iinc = 0
      imult = 0

! Look for "wt" namelist. If not found, stop with error:

      call nmlsrc('wt', iin, ifind)

      if (ifind .eq. 0) goto 1002

      read(iin, wt)

! Formatted-type weight card read:

    else
      read(iin, 9005, end = 1001, err = 1005) type, istep1, istep2, &
      value1, value2, iinc, imult
    end if

! End of weight card read. Now process given options:

! ITYPE = "#" is a spacer card:

    if (type(1:1) .eq. '#') then
      write(iout, 9006) type, istep1, istep2, value1, value2, iinc, imult
      goto 11
    end if

! If ITYPE="END", then goto redirection reading section:

    if (type .eq. flag(iflag) .or. type(1:3) == 'end') then
      if (i .eq. 1) write(iout, 9001)
      write(iout, 9002)
      ichgwt = i - 1
      goto 20
    end if

! If ITYPE="NSTEP0", then reset the initial value of the step counter
! and goto the next restraint:

    if (type .eq. 'NSTEP0') then
      write(iout, 9006) type, istep1, istep2, value1, value2, iinc, imult
      nstep0 = istep1
      goto 11
    end if

! If ITYPE="STPMLT", then reset the value of the step counter multiplier.
! By default, this multiplier has a value of 1.0.

    if (type .eq. 'STPMLT') then
      write(iout, 9006) type, istep1, istep2, value1, value2, iinc, imult
      if (value1 .gt. zerner) stpmlt = value1
      goto 11
    end if

! If ITYPE="DISAVI", "ANGAVI", or "TORAVI", then set various values related
! to the time-integral. Also set the frequency of "dumping" to the restraint
! value dump file.

    if (type .eq. 'DISAVI') then
      navint(1) = 1
      if (idmpav .eq. 0) idmpav = istep2
      ravein(1) = value1
      if (value2 .gt. zerner) taufin(1) = value2
      if (iinc .eq. 1) iavtyp(1) = 2
      write(iout, 9006) type, istep1, istep2, value1, value2, iinc, imult
      goto 11
    else if (type .eq. 'ANGAVI') then
      navint(2) = 1
      if (idmpav .eq. 0) idmpav = istep2
      ravein(2) = value1
      if (value2 .gt. zerner) taufin(2) = value2
      if (iinc .eq. 1) iavtyp(2) = 2
      write(iout, 9006) type, istep1, istep2, value1, value2, iinc, imult
      goto 11
    else if (type .eq. 'TORAVI') then
      navint(3) = 1
      if (idmpav .eq. 0) idmpav = istep2
      ravein(3) = value1
      if (value2 .gt. zerner) taufin(3) = value2
      if (iinc .eq. 1) iavtyp(3) = 2
      write(iout, 9006) type, istep1, istep2, value1, value2, iinc, imult
      goto 11
    else if (type == 'DUMPFREQ') then
      idmpav = istep1
      write(iout, 9006) type, istep1, istep2, value1, value2, iinc, imult
      goto 11
    end if

! Check to make sure maximum allowable number of weight cards not exceeded:

    if (i .gt. maxwt) then
      write(iout, 9008) maxwt
      call mexit(iout, 1)
    end if

! Loop through valid types:

    do j = 1, iflag
      if (type .eq. flag(j)) then
        ! Check for types not now supported under pmemd:
        if (j .eq. 17) then
          write(iout, *)' Error: CUT TYPE not supported by pmemd; line:'
          write(iout, 9006) type, istep1, istep2, value1, value2, iinc, imult
          call mexit(iout, 1)
        end if
#ifdef CUDA
        ! JMS: We do not support scaling potential terms on the GPU, so there
        ! are many options we don't support: BOND ANGLE TORSION VDW HB ELEC NB
        ! ATTRACT REPULSE INTERN ALL RSTAR IMPROP SOFTR(??what is this??)
        if ((j .ge. 1 .and. j .le. 11) .or. (j .ge. 14 .and. j .le. 17)) then
          write(iout, *) 'Error: Cannot vary potential terms on the GPU'
          call mexit(iout, 1)
        end if
#endif /* CUDA */
        goto 14
      endif
    end do

    write(iout, 9010)
    write(iout, 9006) type, istep1, istep2, value1, value2, iinc, imult
    call mexit(iout, 1)

! The type is valid. Store an integer pointer to the type of flag
! in iwttyp(i). Determine the slope and intercept of changes in value
! with step number, and store these in wtnrg. Store the step range
! over which this change is valid in iwtstp.

14 continue

    write(iout, 9006) type, istep1, istep2, value1, value2, iinc, imult
    iwttyp(i) = j
    if (istep1 .gt. istep2) write(iout, 9015)

! Special handling if averaging requested:

    if (j .eq. 24 .or. j .eq. 25 .or. j .eq. 26) then

      jdx = j - 23

! If itimav=0, do not turn on averaging, even if requesting (e.g.
! if this is a min run).

      if (itimav .eq. 0) then
        write(iout, 9013)
        goto 11
      end if

      if (nave(jdx) .eq. 0) nave(jdx) = max(iinc, 1)
      if (ipower(jdx) .eq. 0) ipower(jdx) = nint(value2)
      imult = 0
      if (value1 .lt. zerner) value1 = 1.0d+7
      value2 = value1
    end if

    if (value1 .eq. 0.0d0) value1 = zerner
    if (value2 .eq. 0.0d0) value2 = zerner
    if (value1 .lt. 0.0d0) value1 = 0.0d0
    if (value2 .lt. 0.0d0) value2 = 0.0d0

    if (istep1 .le. istep2 .or. istep2 .le. 0) then
      iwtstp(1, i) = istep1
      iwtstp(2, i) = istep2
    else
      iwtstp(1, i) = istep2
      iwtstp(2, i) = istep1
    end if
    iwtstp(3, i) = max(1, iinc)

! Store special information if TYPE = SHORT:

    if (type(1:5) .eq. 'SHORT') then
      slope = value1
      rint = value2
      iwtstp(1, i) = istep1
      iwtstp(2, i) = istep2
      if (iinc .le. 0) iwtstp(3, i) = 0
      ishrtb = i

    else if (istep2 .le. 0 .or. istep1 .eq. istep2) then
      slope = 0.0d0
      rint = value1

! If imult > 0, then use a constant multipier every iinc steps, instead
! of a linear interpolation. Set iwtstp(3,i) = -iinc, to signal this type
! of change, store the multiplicative factor in wtnrg(1,i), and store
! the initial value in wtnrg(2,i).

    else if (imult .gt. 0) then
      if (abs(value1) .lt. zerner) then
        write(iout, 9011)
        call mexit(iout, 1)
      else if (value2/value1 .lt. 0.0d0) then
        write(iout, 9012)
        call mexit(iout, 1)
      end if
      ntu = (iwtstp(2, i) - iwtstp(1, i))/iwtstp(3, i)
      rntu = dble(max(1, ntu))
      slope = (value2/value1)**(1.0d0/rntu)
      rint = value1
      iwtstp(3, i) = - iwtstp(3, i)
    else
      step1 = dble(istep1)
      step2 = dble(istep2)
      slope = (value2 - value1) / (step2 - step1)
      rint = value2 - slope*step2
    end if
    wtnrg(1, i) = slope
    wtnrg(2, i) = rint

10 continue

! Now read the redirection cards, if any. The first non-blank line encountered
! is assumed to be the start of the redirection section. If it is not a known
! option, assume no redirection cards where given.

20 continue

    read(iin, 9004, end = 25) aline

! Skip blank lines:

    do j = 1, len(aline)
      if (aline(j:j) .ne. ' ') goto 29
    end do

    goto 20

! Now see if this line starts with a valid redirection name:

29 continue

    do m = 1, 8
      do n = 10, 1, - 1
        if (redirc(m)(n:n) .ne. ' ') goto 24
      end do

24 continue

      call chklin(aline, redirc(m), n, imat, istrt, ilengt, iout)
      if (imat .ne. 0) then
        redir(m) = aline(istrt:istrt + ilengt - 1)
        iredir(m) = ilengt
        ihits = ihits + 1
        if (ihits .eq. 1) write(iout, 9069)
        write(iout, 9070) redirc(m), redir(m)(1:ilengt)
        goto 20
      end if
    end do

    if (ihits .eq. 0) write(iout, 9071)

    backspace(iin)

! Now read the restraints:

25 continue

! If restraint input has been redirected, open the appropriate file:

    iothrd = 0
    if (iredir(3) .ne. 0) then
      call opnmrg(redir(3)(1:iredir(3)), iscop1, 2, iout, ierr)
      if (ierr .eq. 1) then
        call mexit(iout, 1)
      end if
      iin = iscop1
      iothrd = iin
      write(iout, 9031) redir(3)(1:iredir(3))

! Read and echo comment cards from this file:

      write(mdout, '(a)') 'Here are comments from the DISANG input file:'

42 continue

      read(iin, '(a)', end = 43) aline

      if (aline(1:1) .eq. '#') then
        write(mdout, '(a)') aline
        goto 42
      end if
      backspace (iin)
      write(mdout, *)

43 continue

    end if

! If echoing of restraints read has been redirected, open the appropriate file:

    iuse = 0
    if (iredir(1) .ne. 0) then
      if (redir(1)(1:iredir(1)) .eq. 'POUT') then
        iuse = iout
      else
        call opnmrg(redir(1)(1:iredir(1)), iscop2, 0, iout, ierr)
        if (ierr .eq. 1) then
          call mexit(iout, 1)
        end if
        iuse = iscop2
        write(iuse, 9032)
        write(iout, 9033) redir(1)(1:iredir(1))
      end if
    end if

! Set the 0 defaults for many variables which can be set by the
! "rst" namelist read. The variables which appear before the "do 30"
! loop take the value of the last &rst namelist that set them. The remainder
! revert to the defaults (0) in the "do 30 loop":

    ir6 = 0
    ialtd = 0
    ninc = 0
    imult = 0

    r1 = 0.0d0
    r2 = 0.0d0
    r3 = 0.0d0
    r4 = 0.0d0
    rk2 = 0.0d0
    rk3 = 0.0d0
    r1a = 0.0d0
    r2a = 0.0d0
    r3a = 0.0d0
    r4a = 0.0d0
    rk2a = 0.0d0
    rk3a = 0.0d0

    do 30 i = 1, 999999
      ! provide ability to initiate all 8 groups
      do ii = 1, 8
        iat(ii) = 0
        atnam(ii) = '    '
      end do

      do ii = 1, 3
        rjcoef(ii) = 0.0d0
      end do


      do ii = 1, maxigr
        igr1(ii) = 0
        igr2(ii) = 0
        igr3(ii) = 0
        igr4(ii) = 0
        igr5(ii) = 0
        igr6(ii) = 0
        igr7(ii) = 0
        igr8(ii) = 0
        grnam1(ii) = '    '
        grnam2(ii) = '    '
        grnam3(ii) = '    '
        grnam4(ii) = '    '
        grnam5(ii) = '    '
        grnam6(ii) = '    '
        grnam7(ii) = '    '
        grnam8(ii) = '    '
      end do

      ixpk = 0
      nxpk = 0
      nstep1 = 0
      nstep2 = 0
      irstyp = 0
      ifvari = 0
      iresid = 0
      ifntyp = 0

      fxyz(:) = 1
      outxyz = 0
! Look for "rst" namelist. If not found, assume we are done reading them:

      call nmlsrc('rst', iin, ifind)

      if (ifind .eq. 0) goto 227

! Read the rst namelist.

      read(iin, rst, end = 227)

! taken from sander
      do ii = 1, maxigr
        igrarr(ii,1) = igr1(ii)
        igrarr(ii,2) = igr2(ii)
        igrarr(ii,3) = igr3(ii)
        igrarr(ii,4) = igr4(ii)
        igrarr(ii,5) = igr5(ii)
        igrarr(ii,6) = igr6(ii)
        igrarr(ii,7) = igr7(ii)
        igrarr(ii,8) = igr8(ii)
        grnamarr(ii,1) = grnam1(ii)
        grnamarr(ii,2) = grnam2(ii)
        grnamarr(ii,3) = grnam3(ii)
        grnamarr(ii,4) = grnam4(ii)
        grnamarr(ii,5) = grnam5(ii)
        grnamarr(ii,6) = grnam6(ii)
        grnamarr(ii,7) = grnam7(ii)
        grnamarr(ii,8) = grnam8(ii)
      end do


! If iat(1)=0, assume the last restraint has been read.

227 continue

      if (iat(1) .eq. 0) then
        nmrnum = i - 1
        if (nmrnum .le. 0) write(iout, 9025)
        if (nmrnum .gt. 0) write(iout, 9026) nmrnum
        write(iout, 9027)
        goto 100
      end if

! figure out the type of restraint
      resttype(i) = 0
      if (iat(3) .eq. 0) then        ! Distance restraint
        resttype(i) = 1
      else if (iat(4) .eq. 0) then   ! Angle restraint
        resttype(i) = 2
      else if (iat(5) .eq. 0) then   ! torsion and j coupling restraint
        resttype(i) = 3
      end if
      ! add the other 3 restraints for future pmemd kelB

      ! Assignment of some variables if jar run (roit. 02/27/05)
      ! We are using the restraint potential as defined for an umbrella sampling
      ! calculation, where nstep2=nstlim, and the minimum of the well is
      ! linearly changed from r2=r3 to r2a=r3a. Constant is rk2=rk3 and
      ! r1 and r4 are taken large enough [r1(r4)=r2(r3)-(+)100]

      if (jar .eq. 1) then

        write(mdout, '(a)') 'jar option running '

        if (nmrnum .ne. 1) goto 1004

        nstep2 = nstlim
        ifvari = 1

        if (r2 .ne. r3) then
          if (r2 .eq. 0.) then
            r2 = r3
          else if (r3 .eq. 0.) then
            r3 = r2
          else
            write(6,'(/2x,a)') &
              'r2 and r3 must be equal, or only one of them should be defined'
          endif
        endif

        if (r2a .ne. r3a) then
          if (r2a .eq. 0.) then
            r2a = r3a
          else if (r3a .eq. 0.) then
            r3a = r2a
          else
            write(6,'(/2x,a)') &
              'r2a and r3a must be equal, or only one of them should be defined'
          endif
        endif

        if (rk2 .ne. rk3) then
          if (rk2 .eq. 0.) then
            rk2 = rk3
          else if (rk3 .eq. 0.) then
            rk3 = rk2
          else
            write(6,'(/2x,a)') &
              'rk2 and rk3 must be equal, or only one of them should be defined'
          endif
        endif

        if (rk2a .ne. rk2) then
          if (rk2a .eq. 0.) then
            rk2a = rk2
          else if (rk2 .eq. 0.) then
            rk2 = rk2a
          else
            write(6,'(/2x,a)') &
              'rk2 and rk2a must be equal, or only one of them should be defined'
          endif
        endif

        if (rk2a .ne. rk3a) then
          if (rk2a .eq. 0.) then
            rk2a = rk3a
          else if (rk3a .eq. 0.) then
            rk3a = rk2a
          else
            write(6,'(/2x,a)') &
              'rk2a and rk3a must be equal, or only one of them should be defined'
          endif
        endif

        if (rk2 .eq. 0.) then
          write(6,'(/2x,a)') 'rk2 should not be 0'
        endif

        r1 = r2 - 100
        r4 = r3 + 100
        r1a = r2a - 100
        r4a = r3a + 100
        drjar = (r2a - r2) / (nstep2 - nstep1)

      endif !jar

! Check to make sure maximum allowable number of weight cards not exceeded

      if (i .gt. maxrst) then
        write(iout, 9028) maxrst
        call mexit(iout, 1)
      end if

! If iresid = 1, read another line and determine the absolute atom numbers
! of the atoms making up the restraint. The values in iat(1)->iat(4) will
! be replaced by these absolute atom numbers.

      if (iresid .eq. 1) then

! Must transfer the character name to an integer variable, if necessary:
        do ii = 1, 8
          read(atnam(ii), '(a4)') iatnm(ii)
          ierr = 0
          call getnat(iat(ii), iatnm(ii), name, gbl_res_atms, nres, iout, ierr)
        end do

        if (ierr .eq. 1) then
          call mexit(iout, 1)
        end if
      end if

! Set flag igravt (GRoup AVerage Type) to ir6:

      igravt(i) = ir6

! Set flag ialtdis (Use alternative restraint function) to ialtd:

      ialtdis(i) = ialtd

! Set the function-type flag to ifntyp:

      nmrfty(i) = ifntyp

! Store the atom pointers (3*(iat-1)) and effective ranges. The atom
! pointers will be changed with the group info., if iat(1) or iat(2) < 0
! (see below).
! similar to sander
      do ii = 1, 8 ! Added atoms pointer for all 8 type restraints
        nmrat(ii, i) = 3*(iat(ii) - 1)

! Store the fxyz coefficient in a global array. Default value is 1.d0
        ifxyz(5*(i-1) + 1) = fxyz(1)
        ifxyz(5*(i-1) + 2) = fxyz(2)
        ifxyz(5*(i-1) + 3) = fxyz(3)
        ! Check if we only have single direction in either x, y or z 
        ! if we do then set the 4th element to single direction index
 
        if(fxyz(1) + fxyz(2) + fxyz(3) .eq. 1) then
           if(fxyz(1) .eq. 1) then
             ifxyz(5*(i-1) + 4) = 1
           end if
           if(fxyz(2) .eq. 1) then
             ifxyz(5*(i-1) + 4) = 2
           end if
           if(fxyz(3) .eq. 1) then
             ifxyz(5*(i-1) + 4) = 3
           end if
        ! else we have multidirection so set 4th element to zero
        else
          ifxyz(5*(i-1)+4) = 0
        end if
        ifxyz(5*(i-1) + 5) = 0
     end do
!RCW idxyz used to be defined above as idxyz(1) implying it was
!    1 element long despite being an array. The following zero
!    would therefore, I believe, just zero the first element and
!    not the entire array. Now it is correctly set to be assumed size
!    so we manually just zero the first element and assume this is what
!    the original author intended. Who knows??? How did we get away
!    with having all these arrays defined as (1) for so long?
!
!    As far as I can tell there is no requirement to zero the first
!    element of this array.
!
!      idxyz(:) = 0.d0

      ioutxyz(i) = outxyz

! Initialize nmrat(9, i) to nmrat(16,i) to 0, these eventually hold the last
! atom index for groups
      do ii = 9, 16
        nmrat(ii, i) = 0
      end do

      if (nstep2 .gt. nstep1 .or. nstep2 .le. 0) then
        nmrst(1, i) = nstep1
        nmrst(2, i) = nstep2
      else
        nmrst(1, i) = nstep2
        nmrst(2, i) = nstep1
      end if
      nmrst(3, i) = max(1, ninc)
      kxpk(1, i) = ixpk
      kxpk(2, i) = nxpk

! Store the Karplus equation constants:

      ajcoef(1, i) = rjcoef(1)
      ajcoef(2, i) = rjcoef(2)
      ajcoef(3, i) = rjcoef(3)

      usecom = .false. !added this logical to track whether a group is involved in the restraint

! Since option was added for iat(1) to iat(4) to  be less than 0, then we do
! not need to be checking if iat(3)=iat(4) = 0, DO a loop in the future
! should we will still check for to ensure iat(5) = 0, just incase?

      do ii = 1, 8
        if (iat(ii) .lt. 0) then
          do jj =1, maxigr
            igrpass(jj) = igrarr(jj, ii)
            grnampass(jj) = grnamarr(jj, ii)
          end do
          call nmrgrp(i, ii-1, lastpk, natom, nmrat, nmrcom, iout, maxgrp, igrpass, &
                    grnampass, iresid, name, gbl_res_atms, nres, maxigr)
        end if
        ! usecom tracks whether a group is involved in the restraint
        usecom = usecom .or. nmrat(ii+8, i) .lt. 0
      end do

! Calculate the current value of the internal. If irstyp=1, this
! will be used to calculate the target distance. If irstyp=0, it
! will be used simply to report the deviation between the current & target
! values.

      jcoupl = .false.
      resttypetmp = resttype(i)
      select case (resttypetmp)
        !for distance stuff
        case (1)  ! distance
          if (.not. usecom) then
            call disnrg(crd, frc_dummy, dfr, rcurr, nmrat(1, i), nmrat(2, i), e, &
                        r1, r2, r3, r4, rk2, rk3, dumm_array, dumm_array, &
                        ipower(1), dumm, ravein(1), dumm, 0, 0, 0, 0, 1, &
                        ialtdis(i), ifxyz(5*(i - 1) + 1), idxyz(3*(i - 1) +1))

          else if (igravt(i) .eq. 1) then !r6 distance
            !call r6ave(crd, nmrat, nmrcom, rcurr, i)
            call r6ave(crd, nmrat, nmrcom, rcurr, i, ifxyz(5*(i - 1) + 1), &
                      idxyz(3*(i - 1) + 1))
            call disnrg(crd, frc_dummy, dfr, rcurr, nmrat(1, i), nmrat(2, i), e, &
                        r1, r2, r3, r4, rk2, rk3, dumm_array, dumm_array, &
                        ipower(1), dumm, ravein(1), dumm, 0, 0, 0, 0, 3, &
                        ialtdis(i), ifxyz(5*(i - 1) + 1), idxyz(3*(i - 1) + 1))
          else ! COM distance
            call nmrcms(crd, xcom, nmrat, nmrcom, rmstot, i)
            call disnrg(xcom, frc_dummy, dfr, rcurr, 0, 3, e, r1, r2, r3, r4, &
                        rk2, rk3, dumm_array, dumm_array, ipower(1), dumm, &
                        ravein(1), dumm, 0, 0, 0, 0, 1, ialtdis(i), ifxyz(5*(i - 1) + 1), &
                        idxyz(3*(i - 1) + 1))
          end if ! end distance type
        case (2) ! angle stuff
          if (usecom) then !COM angles
            call nmrcms(crd, xcom, nmrat, nmrcom, rmstot, i)
            call angnrg(xcom, frc_dummy, dfr, rcurr, 0, 3,            &
                        6, e, r1, r2, r3, r4, rk2, rk3,                 &
                        dumm_array, dumm_array, ipower(2), dumm,        &
                        ravein(2), dumm, 0, 0, 0, 0, 1)
          else ! regular angles
            call angnrg(crd, frc_dummy, dfr, rcurr, nmrat(1, i), nmrat(2, i), &
                        nmrat(3, i), e, r1, r2, r3, r4, rk2, rk3, &
                        dumm_array, dumm_array, ipower(2), dumm, ravein(2), &
                        dumm, 0, 0, 0, 0, 1)

          end if  ! end of angle type
        case (3) ! Torsion and Jcoupling stuff
! If J coefficients are 0 for this torsion, it is a normal torsional
! restraint. If they are .ne.0, then it is a J-coupling factor restraint.
          jcoupl = abs(ajcoef(1,i)) .gt. zerner .or. abs(ajcoef(2,i)) .gt. zerner
          if (.not. jcoupl) then
            if (usecom) then !COM torsion
              call nmrcms(crd, xcom, nmrat, nmrcom, rmstot, i)
              call tornrg(xcom , frc_dummy, dfr, rcurr, 0,        3,     &
                          6, 9, e, r1, r2, r3, r4, rk2, rk3,             &
                          dumm_array, dumm_array, dumm_array, dumm_array, &
                          ipower(3), dumm, ravein(3), dumm,               &
                          0, 0, 0, 0, 1)
            else ! regular torsion
              call tornrg(crd, frc_dummy, dfr, rcurr, nmrat(1, i), nmrat(2,i), &
                          nmrat(3, i), nmrat(4, i), e, r1, r2, r3, r4, rk2, &
                          rk3, dumm_array, dumm_array, dumm_array, dumm_array, &
                          ipower(3), dumm, ravein(3), dumm, 0, 0, 0, 0, 1)
            end if  ! end com and regular torsion
          else
            if (usecom) then !j coupling
             ! Since atom groups are not supported for j-coupling restraints,
             ! print out an error message and terminate execution if the user
             ! specifies atom groups AND j-coupling
              write(iout,9074)
              call mexit(iout, 1)
            else
              call jnrg(crd, frc_dummy, rcurr, nmrat(1, i), nmrat(2, i), &
                        nmrat(3, i), nmrat(4, i), e, r1, r2, r3, r4, &
                        rk2, rk3, dumm_array, dumm_array, ipower(1), dumm, &
                        ravein(1), dumm, 0, 0, 0, 0, ajcoef(1, i), 1)
            end if ! end j coupling or not torsion
          end if ! end torsion types
        !case (4)
           ! add functionality for plane-point angle restraint, see sander
        !case (5)
           ! add functionality for generalized distance restraint, see sander
        !case (6)
           ! add functionality for plane plane angle restraint, see sander
        case default
        ! Improper restraint specified
        write(iout,9075)
        call mexit(iout, 1)
      end select

      rint = 0.0d0
      if (irstyp .eq. 1) rint = rcurr

! Echo the restraint to the user, if the user has requested this:
      do ii = 1, 8
        if (iat(ii) .lt. 0) then
          nametmp(ii) = 'COM '
        else if (iat(ii) .gt. 0) then
          nametmp(ii) = name(iat(ii))
        else
          nametmp(ii) = '????'
        end if
      end do

      if (iuse .ne. 0) then
        resttypetmp = resttype(i)
        select case (resttypetmp)
          case (1) ! distances
            write(iuse,9040) nametmp(1),iat(1),nametmp(2),iat(2), &
                  nstep1,nstep2
          case (2)  !angles
            write(iuse,9041) nametmp(1),iat(1),nametmp(2),iat(2), &
                 nametmp(3),iat(3),nstep1,nstep2
          case (3) ! torsions
            write(iuse,9042) nametmp(1),iat(1),nametmp(2),iat(2), &
                 nametmp(3),iat(3),nametmp(4),iat(4),nstep1,nstep2
          ! case (4) case (5) case (6)
        end select

        if (ninc .gt. 1) write(iuse, 9043) ninc

! If a group-type distance restraint was defined, echo the sorted group:
        do ii = 1, 8 ! only four restraint type but still use 8 for future
          if (nmrat(ii, i) .lt. 0 .and. nmrat(ii+8, i) .lt. 0 .and. ii .eq. 1) then
            write(iuse, 9065)
            write(iuse, 9066) (nmrcom(1, kk), nmrcom(2, kk), &
                  kk=-nmrat(ii, i),-nmrat(ii+8, i))

          else if (nmrat(ii, i) .lt. 0 .and. nmrat(ii+8, i) .lt. 0 .and. ii .eq. 2) then
            write(iuse, 9067)
            write(iuse, 9066) (nmrcom(1, kk), nmrcom(2, kk), &
                  kk=-nmrat(ii, i),-nmrat(ii+8, i))

          else if (nmrat(ii, i) .lt. 0 .and. nmrat(ii+8, i) .lt. 0 .and. ii .eq. 3) then
            write(iuse, 9076)
            write(iuse, 9066) (nmrcom(1, kk), nmrcom(2, kk), &
                  kk=-nmrat(ii, i),-nmrat(ii+8, i))

          else if (nmrat(ii, i) .lt. 0 .and. nmrat(ii+8, i) .lt. 0 .and. ii .eq. 4) then
            write(iuse, 9077) ii
            write(iuse, 9066) (nmrcom(1, kk), nmrcom(2, kk), &
                  kk=-nmrat(ii, i),-nmrat(ii+8, i))
          end if
        end do
      end if !(iuse /= 0)

      convrt = 1.0d0


      if ((iat(3) .ne. 0 .or. iat(4) .ne. 0) .and. .not. jcoupl)  convrt = 180./PI
      rprt = rint*convrt
      rcurr = rcurr*convrt

! If this is a torsional restraint, translate the value of the torsion
! (by +- N*360) to bring it as close as possible to one of the two central
! "cutoff" points (r2,r3). Use this as the value of the torsion in the report.

      if (iat(4) .ne. 0 .and. .not.jcoupl) then
        apmean = (r2 + rprt + r3 + rprt)*0.5
   18   if (rcurr - apmean .gt. 180.0d0) then
          rcurr = rcurr - 360.0d0
          goto 18
        else if (apmean - rcurr .gt. 180.0d0) then
          rcurr = rcurr + 360.0d0
          goto 18
        end if
      end if

      if (iuse .ne. 0) then
        write(iuse, 9045) r1 + rprt, r2 + rprt, r3 + rprt, r4 + rprt, rk2, rk3
        if (ifvari .gt. 0) &
          write(iuse, 9046) r1a + rprt, r2a + rprt, r3a + rprt, r4a + rprt, &
                            rk2a, rk3a

        dev2 = min(abs(rcurr - r2 - rprt), abs(rcurr - r3 - rprt))
        if (rcurr .ge. r2 + rprt .and. rcurr .le. r3 + rprt) dev2 = 0
        write(iuse, 9047) rcurr, abs(rcurr - (r2 + rprt + r3 + rprt)/2.0d0), &
                          dev2
      end if

! The energy routines only work correctly if r1->r4 increase monotonically.
! If not, flag an error -- dap (8/91)

      if (r1 .gt. r2 .or. r2 .gt. r3 .or. r3 .gt. r4) then
        goto 1003
      end if
      if (nstep2 .ne. 0 .and. ifvari .gt. 0 .and. nstep1 .ne. nstep2) then
         if (r1a .gt. r2a .or. r2a .gt. r3a .or. r3a .gt. r4a) then
           goto 1003
         end if
      end if

! Convert r1->r4 to radians, if it is an angle:



      if ((iat(3) .ne. 0 .or. iat(4) .ne. 0) .and. .not. jcoupl) then
        r1 = r1/convrt
        r2 = r2/convrt
        r3 = r3/convrt
        r4 = r4/convrt
        r1a = r1a/convrt
        r2a = r2a/convrt
        r3a = r3a/convrt
        r4a = r4a/convrt
      end if

! Now calculate the slope/intercepts for the r and k values with step
! number and store these in the appropriate arrays:

      if (nstep2 .eq. 0 .or. ifvari .le. 0 .or. nstep1 .eq. nstep2) then
        r1nmr(1, i) = 0.0d0
        r1nmr(2, i) = r1 + rint
        r2nmr(1, i) = 0.0d0
        r2nmr(2, i) = r2 + rint
        r3nmr(1, i) = 0.0d0
        r3nmr(2, i) = r3 + rint
        r4nmr(1, i) = 0.0d0
        r4nmr(2, i) = r4 + rint
        rk2nmr(1, i) = 0.0d0
        rk2nmr(2, i) = rk2
        rk3nmr(1, i) = 0.0d0
        rk3nmr(2, i) = rk3
      else
        step1 = dble(nstep1)
        step2 = dble(nstep2)
        r1nmr(1, i) = (r1a - r1) / (step2 - step1)
        r1nmr(2, i) = (r1a + rint) - r1nmr(1, i) * step2
        r2nmr(1, i) = (r2a - r2) / (step2 - step1)
        r2nmr(2, i) = (r2a + rint) - r2nmr(1, i) * step2
        r3nmr(1, i) = (r3a - r3) / (step2 - step1)
        r3nmr(2, i) = (r3a + rint) - r3nmr(1, i) * step2
        r4nmr(1, i) = (r4a - r4) / (step2 - step1)
        r4nmr(2, i) = (r4a + rint) - r4nmr(1, i) * step2

! If imult > 0, then use a constant multipier every iinc steps, instead
! of a linear interpolation for the force constants:

        if (imult .gt. 0) then
          ntu = (nmrst(2, i) - nmrst(1, i))/nmrst(3, i)
          rntu = dble(max(1, ntu))
          nmrst(3, i) = - nmrst(3, i)
          rk2nmr(1, i) = (rk2a/rk2)**(1.0d0/rntu)
          rk2nmr(2, i) = rk2
          rk3nmr(1, i) = (rk3a/rk3)**(1.0d0/rntu)
          rk3nmr(2, i) = rk3
        else
          rk2nmr(1, i) = (rk2a - rk2) / (step2 - step1)
          rk2nmr(2, i) = rk2a - rk2nmr(1, i) * step2
          rk3nmr(1, i) = (rk3a - rk3) / (step2 - step1)
          rk3nmr(2, i) = rk3a - rk3nmr(1, i) * step2
        end if
      end if !(nstep2 == 0 .or. ifvari <= 0 .or. nstep1 == nstep2)

! Unconvert r1->r4 to degrees, if it is an angle (this is done in case the
! next namelist does not reset these variables, with the user expecting the
! same values to be repeated for several constraints).

! changing .gt. to .ne. to account for COM angle + torsion restraints

    if ((iat(3) .ne. 0 .or. iat(4) .ne. 0) .and. .not.jcoupl) then
      r1 = r1*convrt
      r2 = r2*convrt
      r3 = r3*convrt
      r4 = r4*convrt
      r1a = r1a*convrt
      r2a = r2a*convrt
      r3a = r3a*convrt
      r4a = r4a*convrt
    end if

30 continue
! Call nmrsht to set up short-range definitions. set nave(1) (idumm(1) in
! call list) to zero, so that time-averaged distances are not used on this
! (initial) call.

    idumm(1) = 0
    idumm(2) = 0
    idumm(3) = 0

100 continue

  call nmrsht(crd, nmrnum, 0, nmrat, nmrcom, ishrtb, iwtstp, wtnrg, rk2nmr, &
              rk3nmr, nmrst, nmrshb, wtls, nmrfty, igravt, idumm, &
              dumm_array, dumm_array, idumm, idumm, dumm_array, idumm, &
              dumm_array, dumm, iout, 0, ifxyz, idxyz, resttype)

  if (iothrd .ne. 0) call opnmrg(redir(3), iin, 3, iout, ierr)
  if (iuse .ne. 0 .and. iuse .ne. iout)  &
      call opnmrg(redir(1), iuse, 3, iout, ierr)
  return

! End-of-file read errors:

1001 continue

  write(iout, 2000)
  call mexit(iout, 1)

! No "wt" namelist with type=END read:

1002 continue

  write(iout, 2002)
  call mexit(iout, 1)

! r1->r4 were not monotonic error:

1003 continue

  write(iout, 2003)
    do ii = 1, 8
      if (iat(ii) .lt. 0) then
        nametmp(ii) = 'COM '
      else
        nametmp(ii) = name(iat(ii))
      end if
    end do
    resttypetmp = resttype(i)
    select case (resttypetmp)
      case (1) ! distances
        write(iout,9040) nametmp(1),iat(1),nametmp(2),iat(2), &
                 nstep1,nstep2
      case (2) ! angles
        write(iout,9041) nametmp(1),iat(1),nametmp(2),iat(2), &
                 nametmp(3),iat(3),nstep1,nstep2
      case (3) ! torsions
        write(iout,9042) nametmp(1),iat(1),nametmp(2),iat(2), &
                 nametmp(3),iat(3),nametmp(4),iat(4),nstep1,nstep2
      !other cases here
    end select

  write(iout, 9045) r1 + rprt, r2 + rprt, r3 + rprt, r4 + rprt, rk2, rk3
  if (ifvari .gt. 0) &
    write(iout, 9046) r1a + rprt, r2a + rprt, r3a + rprt, r4a + rprt, rk2a, rk3a
  call mexit(iout, 1)

! Other read errors (probably format mismatch):
! (If "#" is in first column, this is a comment line):

1004 continue

  write(iout, 9072) !jarzynski error exit (roit. 02/27/05)
  call mexit(iout, 1)

1005 continue

  backspace(iin)
  read(iin, 9004) aline
  if (aline(1:1) .eq. '#') then
    write(iout, 9068) aline
    ibeg = i + 1
    goto 9
  else
    write(iout, 2001) aline
    call mexit(iout, 1)
  end if

2000 format(' ERROR: End of file read during formatted weight card ', 'read.')
2001 format(' ERROR reading the line:', /, a80, /, &
            ' during formatted weight-card read.')
2002 format(' ERROR: No "wt" namelist with TYPE=END found')
2003 format(' ERROR: r1 -> r4 (and r1a -> r4a) must be monotonically ', &
            'increasing; Offending restraint:')
9000 format(/ /, t12, 'Begin reading energy term weight changes/NMR ', &
            'restraints', /, ' WEIGHT CHANGES:')
9001 format(t26, '** No weight changes given **')
9002 format(/, ' RESTRAINTS:')
9004 format(a80)
9005 format(a8, 2i7, 2f12.6, 2i7)
9006 format(1x, a8, 2i7, 2f12.6, 2i7)
9008 format(' Error: Maximum allowable number of weight change cards', &
            ' (', i5, ') exceeded.')
9010 format(' Error: Invalid TYPE flag in line:')
9011 format(' Error: Multiplicative-type change (IMULT > 0) ', &
            'invalid when VALUE1 = 0.0.')
9012 format(' Error: Multiplicative-type change (IMULT > 0) invalid when ', &
            'VALUE1 and', /, t8, 'VALUE2 do not have the same sign.')
9013 format(' Warning: Averaging request ignored.')
9015 format(' Warning: ISTEP1 >= ISTEP2')
9025 format(t27, '** No restraint defined **')
9026 format(t23, ' Number of restraints read = ', i5)
9027 format(/, t19, 'Done reading weight changes/NMR restraints', / /)
9028 format(' Error: Maximum allowable number of restraints (',i5,') exceeded.')
9031 format(' Restraints will be read from file: ', a)
9032 format(' INITIAL restraints/deviations/energy contributions:')
9033 format(' Initial restraints echoed to file: ', a)
9040 format('******', /, &
            1x, a4, '(', i5, ')-', a4, '(', i5, ')', t53, 'NSTEP1=', i6, &
            ' NSTEP2=', i6)
9041 format('******', /, &
            1x, a4, '(', i5, ')-', a4, '(', i5, ')-', a4, '(', i5, ')',  &
            t53, 'NSTEP1=', i6, ' NSTEP2=', i6)
9042 format('******', /, &
            1x, a4, '(', i5, ')-', a4, '(', i5, ')-', a4, '(', i5, ')-', &
            a4, '(', i5, ')', t53, 'NSTEP1=', i6, ' NSTEP2=', i6)
9043 format(1x, 'NINC=', i6)
9045 format('R1 =', f8.3, ' R2 =', f8.3, ' R3 =', f8.3, ' R4 =', f8.3, &
            ' RK2 =', f8.3, ' RK3 = ', f8.3)
9046 format('R1A=', f8.3, ' R2A=', f8.3, ' R3A=', f8.3, ' R4A=', f8.3, &
            ' RK2A=', f8.3, ' RK3A= ', f8.3)

9047 format(1x, 'Rcurr: ', f8.3, '  Rcurr-(R2+R3)/2: ', f8.3, &
            '  MIN(Rcurr-R2,Rcurr-R3): ', f8.3)

9065 format(' Atom ranges defining Center of Mass Group in &
             &the first position: ')
9066 format(5(1x, i5, ' -> ', i5, '/'))
9067 format(' Atom ranges defining Center of Mass Group in &
             &the second position: ')
9068 format(1x, a79)
9069 format(' Requested file redirections:')
9070 format(t3, a, t13, '= ', a)
9071 format(t3, 'No valid redirection requests found')
9072 format('In Jarzynsky runs, there must only one restrain, stopping program')
9073 format('Error: Invalid residue number or atom name in restraint.')
9074 format('Error: Atom groupings not supported for j-coupling restraints.')
9075 format('Error: Improper restraint specified.')
9076 format(' Atom ranges defining the Center of Mass Group in the third position: ')
9077 format(' Atom ranges defining the Center of Mass Group in the ',i1,'th position: ')
9078 format('In Jarzynsky runs using COM restraints, GPU support is not available, stopping program')

end subroutine nmrred

!*******************************************************************************
!
! Subroutine:   modwt (MODify WeighTs)
!
! Description:
!
! This routine modifies the various energy coefficients as requested.
!
! Author: David A. Pearlman
! Date: 7/89
!
! INPUT
! -----
!  wtnrg(2,i): Slope/intercept describing how the relative weights change
!              with step number.
! iwtstp(3,i): The beginning and ending steps over which a modification
!              is to be active. If iwtstp(3,i)>0, it gives the step increment
!              between changes in the target value.
!   iwttyp(i): The type of modification I.
!      ichgwt: Number of relative weight change instructions.
!      ishrtb: If ishrtb >0, short-range interactions have been defined;
!              Relevant information is then stored in iwtstp(,ishrtB) and
!              wtnrg(,ishrtb).
!       nstep: Current iteration/step number
!       temp0: Target temperature (molecular dynamics only)
!       tautp: Berendsen temperature scaling parameter (md only)
!       rk(i): Force constants for bonds
!       tk(i): Force constants for angles
!       pk(i): Force constants for torsions. Regular torsions, followed by
!              impropers.
!      cn1(i): "a" (r**12) coefficient array for vdw term
!      cn2(i): "b" (r**6) coefficient array for vdw term
!       ag(i): "a" (r**12) coefficient for h-bond term
!       bg(i): "b" (r**10) coefficient for h-bond term
! rk2nmr(2,i): Slope/intercept of line describing how k2 force constants for NMR
!              restraints vary with step number.
! rk3nmr(2,i): Slope/intercept of line describing how k3 force constants for NMR
!              restraints vary with step number.
!   nmrshb(i): If nmrshb(i) >0, then restraint i has been defined as a
!              "short-range" restraint.
!  nmrst(3,i): nmrst(3,i)>0 if force constants for NMR restraint i are linearly
!              interpolated. nmrst(3,i) < 0 if force constants for NMR restraint
!              i are generated by constant factor multiplication. If the
!              latter is true, do not modify rk?nmr(1,i) when restraint weights
!              are modified.
!
!      numbnd: Number of bond parmameters
!      numang: Number of angle parameters
!      numphi: Number of torsional parameters
!      nimprp: Number of "improper" torsional parameters (not included in numphi
!        nphb: Number of hydrogen bond-parameters
!      ncharg: Number of charge parameters. Equal to the number of atoms
!              for minimization and regular molecular dynamics. Equal to
!              2*Number of atoms for perturbed systems (GIBBS).
!      ntypes: Number of non-bonded 6-12 parameter types.
!       rwell: The force constant used in the "soft-repulsion" function
!   tauave(3): The exponential decay constants used if time-averaged restraints
!              are implemented
!      nmrnum: Number of added NMR restraints.
!
! Global variables:
!
!       atm_qterm(i): Charge on atom i.
!
! Output:
!     wtls(2): wtls(1) is set = weight(10)  (used in resetting short-range)
!              wtls(2) is set = weight(11)  (used in resetting short-range)
!
! Local Common:
! Array weight--
!  weight(1) : The current weight of the bonds
!  weight(2) : The current weight of the angles
!  weight(3) : The current weight of the torsions
!  weight(4) : The current weight of the attractive 6 (vdw) term
!  weight(5) : The current weight of the replusive 12 (vdw) term
!  weight(6) : The current weight of the attractive 10 h-bond term
!  weight(7) : The current weight of the replusive 12 h-bond term
!  weight(8) : The current weight of the electrostatic terms
!  weight(9) : The current weight of the "improper" torsions
!  weight(10) : The current weight of the "short-range" NMR restraints
!  weight(11) : The current weight of the non-"short-range" NMR restraints
!  weight(12) : The current weight of the NOESY restraints
!  weight(13) : The current weight of the chemical shift restraints
!  weight(14) : Temperature scaling parameter on the first call
!  weight(15): The "soft-repulsion" force constant as specified in the input
!  weight(16) : The value of the target temperature on the first call
!  weight(17) : NO LONGER USED IN PMEMD.  Was square of cutoff.
!  weight(18) : Exponential decay time const for time-averaged bond restraints
!  weight(19) : Exponential decay time const for time-averaged angle restraints
!  weight(20) : Exponential decay time const for time-averaged torsion restraints
!
!*******************************************************************************

subroutine modwt(wtnrg, iwtstp, iwttyp, ichgwt, ishrtb, nstep, &
                 rk, tk, pk, cn1, cn2, ag, bg, rk2nmr, rk3nmr, &
                 nmrshb, nmrst, numphi, nimprp_arg, nhb, ncharg, &
                 rwell, tauave, wtls, nmrnum, tautp, temp0)

  use file_io_mod
  use file_io_dat_mod
  use nmr_lib_mod
  use prmtop_dat_mod
#ifdef MPI
  use remd_mod, only : remd_modwt
#endif

  implicit none

! Formal arguments:

  double precision  wtnrg(2, *)
  integer           iwtstp(3, *)
  integer           iwttyp(*)
  integer           ichgwt
  integer           ishrtb
  integer           nstep
  double precision  rk(*)
  double precision  tk(*)
  double precision  pk(*)
  double precision  cn1(*)
  double precision  cn2(*)
  double precision  ag(*)
  double precision  bg(*)
  double precision  rk2nmr(2, *)
  double precision  rk3nmr(2, *)
  integer           nmrshb(*)
  integer           nmrst(3, *)
  integer           numphi
  integer           nimprp_arg
  integer           nhb
  integer           ncharg
  double precision  rwell
  double precision  tauave(3)
  double precision  wtls(2)
  integer           nmrnum
  double precision  tautp
  double precision  temp0

! Local variables:

  double precision small, wt, wt2, wt3

  logical fixed

  integer       jweit, ionce, i, itorfc, irstyp, itype, nstepu, irstr, j, indx

  double precision  rmult

  parameter (jweit = 20)
  parameter (small = 1.0d-12)

  integer           ichang(jweit)
  integer           ichold(jweit)
  integer           ireset(jweit)
  double precision  weight(jweit)


  data ionce/0/
  save weight, ichold

! On the first call, initialize the ichold indicator. On subsequent calls,
! ichold(i)=0, if weight of term was set/remained at 1.0 on previous call.
! ichold(i)=1, if weight of term was set to non-unity as of the previous call.
! ichold(i)=2, if weight of term has been set to 0, or if the weight
!              of the term is constant throughout the run (nstep1=nstep2=0).

#ifdef MPI
  if (ionce .ne. 99 .or. remd_modwt) then
#else
  if (ionce .ne. 99) then
#endif
    do i = 1, jweit
      ichold(i) = 0
      weight(i) = 1.0d0
    end do
    weight(14) = tautp
    weight(15) = rwell
    weight(16) = temp0
    weight(17) = 0.d0           ! Cutoff no longer used!
    weight(18) = tauave(1)
    weight(19) = tauave(2)
    weight(20) = tauave(3)
    ionce = 99
  end if
  itorfc = 0

! Initialize the change indicator. Any energy term not modified by one
! of the change instructions should be set to the default weight of 1.0
! at the end of this routine.

  do i = 1, jweit
    ichang(i) = 0
    ireset(i) = 0
  end do

! Loop over the requested changes. Use an implied loop, rather than
! an explicit one, so that this code can be conveniently used to
! reset unchanged parameter weights to 1.0 later on.

  if (ichgwt .le. 0) goto 75

  irstyp = 0
  i = 1

20 continue

  fixed = (iwtstp(2, i) .eq. 0)

! Skip, if the current step is outside of the desired range:

  if (nstep .lt. iwtstp(1, i) .or. (nstep .gt. iwtstp(2, i).and. &
      iwtstp(2, i) .ne. 0)) goto 71

  itype = iwttyp(i)

! Modify the weight of the appropriate term. Changes are step-wise,
! every iwtstp(3,i) steps (iwtstp(3,i)=1, if user did not specify it).

! If iwtstp(3,i) > 0, then the weight is linearly interpolated.
! If iwtstp(3,i) < 0, then the weight is modified by a multplicative factor.

  if (iwtstp(3, i) .gt. 0) then
    nstepu = nstep - mod(nstep - iwtstp(1, i), iwtstp(3, i))
    wt = dble(nstepu)*wtnrg(1, i) + wtnrg(2, i)
  else if (iwtstp(3, i) .lt. 0) then
    nstepu = (nstep - iwtstp(1, i))/abs(iwtstp(3, i))
    wt = wtnrg(2, i) * wtnrg(1, i)**nstepu
  end if
  wt2 = wt

! TYPE = RSTAR:
! In this case, set wt to the appropriate R**12 weight; set wt2 to the
! appropriate attractive R**6 weight, and set ITYPE=4 (VDW). Also set
! irstr = 1 and wt3 = wt**10. This will ensure h-bond parameters will
! also be appropriately modified. The actual weight changes will be handled
! by the subsequent loops.

  irstr = 0
  if (itype .eq. 14) then
    irstr = 1
    itype = 4
    wt2 = wt**6.0d0
    wt3 = wt**10.0d0
    wt = wt2*wt2
  end if


! TYPE = BOND/INTERN/ALL:

25 continue

  if ((itype .eq. 1 .or. itype .eq. 10 .or.  &
       itype .eq. 11 .or. ireset(1) .eq. 1) .and.ichold(1) .ne. 2) then
    rmult = wt/weight(1)
    do j = 1, numbnd
      rk(j) = rk(j)*rmult
    end do
    weight(1) = wt
    ichang(1) = 1
    if (abs(wt) .le. small) ichold(1) = 2
    if (fixed) ichold(1) = 2
  end if

! TYPE = ANGLE/INTERN/ALL:

  if ((itype .eq. 2 .or. itype .eq. 10 .or.  &
       itype .eq. 11.or. ireset(2) .eq. 1) .and. ichold(2) .ne. 2) then
    rmult = wt/weight(2)
    do j = 1, numang
      tk(j) = tk(j)*rmult
    end do
    weight(2) = wt
    ichang(2) = 1
    if (abs(wt) .le. small) ichold(2) = 2
    if (fixed) ichold(2) = 2
  end if

! TYPE = TORSION/INTERN/ALL:

  if ((itype .eq. 3 .or. itype .eq. 10 .or.  &
       itype .eq. 11.or. ireset(3) .eq. 1) .and. ichold(3) .ne. 2) then
    rmult = wt/weight(3)
    do j = 1, numphi
      pk(j) = pk(j)*rmult
    end do
    weight(3) = wt
    ichang(3) = 1
    if (abs(wt) .le. small) ichold(3) = 2
    if (fixed) ichold(3) = 2
    itorfc = 1
  end if

! TYPE = VDW/ATTRACT/NB/ALL:

  if ((itype .eq. 4 .or. itype .eq. 7 .or. itype .eq. 8 .or. &
       itype .eq. 11.or. ireset(4) .eq. 1) .and. ichold(4) .ne. 2) then
    rmult = wt2/weight(4)
    do j = 1, nttyp
      cn2(j) = cn2(j)*rmult
    end do
    weight(4) = wt2
    ichang(4) = 1
    if (abs(wt2) .le. small) ichold(4) = 2
    if (fixed) ichold(4) = 2
  end if

! TYPE = VDW/REPULSE/NB/ALL:

  if ((itype .eq. 4 .or. itype .eq. 7 .or. itype .eq. 9 .or. &
       itype .eq. 11.or. ireset(5) .eq. 1) .and. ichold(5) .ne. 2) then
    rmult = wt/weight(5)
    do j = 1, nttyp
      cn1(j) = cn1(j)*rmult
    end do
    weight(5) = wt
    ichang(5) = 1
    if (abs(wt) .le. small) ichold(5) = 2
    if (fixed) ichold(5) = 2
  end if

! TYPE = HB/ATTRACT/NB/ALL:

  if ((irstr .eq. 1 .or. itype .eq. 5 .or. itype .eq. 7 .or.  &
       itype .eq. 8 .or. itype .eq. 11 .or. ireset(6) .eq. 1) .and.  &
       ichold(6) .ne. 2) then

    rmult = wt/weight(6)
    if (irstr .eq. 1) rmult = wt3/weight(6)
    do j = 1, nhb
      bg(j) = bg(j)*rmult
    end do
    weight(6) = wt
    ichang(6) = 1
    if (abs(wt) .le. small) ichold(6) = 2
    if (fixed) ichold(6) = 2
  end if

! TYPE = HB/REPULSE/NB/ALL:

  if ((irstr .eq. 1 .or. itype .eq. 5 .or. itype .eq. 7 .or.  &
       itype .eq. 9 .or. itype .eq. 11 .or. ireset(7) .eq. 1) .and.  &
       ichold(7) .ne. 2) then

    rmult = wt/weight(7)
    do j = 1, nhb
      ag(j) = ag(j)*rmult
    end do
    weight(7) = wt
    ichang(7) = 1
    if (abs(wt) .le. small) ichold(7) = 2
    if (fixed) ichold(7) = 2
  end if

! TYPE = ELEC/NB/ALL:

  if ((itype .eq. 6 .or. itype .eq. 7 .or.  &
       itype .eq. 11.or. ireset(8) .eq. 1) .and. ichold(8) .ne. 2) then
    rmult = sqrt(wt/weight(8))
    do j = 1, ncharg
      atm_qterm(j) = atm_qterm(j) * rmult
    end do
    weight(8) = wt
    ichang(8) = 1
    if (abs(wt) .le. small) ichold(8) = 2
    if (fixed) ichold(8) = 2
  end if

! TYPE = IMPROP:

  if ((itype .eq. 15 .or. ireset(9) .eq. 1) .and. ichold(9) .ne. 2) then
    rmult = wt/weight(9)
    do j = 1, nimprp_arg
      pk(numphi + j) = pk(numphi + j)*rmult
    end do
    weight(9) = wt
    ichang(9) = 1
    if (abs(wt) .le. small) ichold(9) = 2
    if (fixed) ichold(9) = 2
    itorfc = 1
  end if

! TYPE = REST/RESTS:

  if ((itype .eq. 12 .or. itype .eq. 19 .or. ireset(10) .eq. 1) .and. &
       ishrtb .gt. 0 .and. ichold(10) .ne. 2) then
    rmult = wt/weight(10)
    do j = 1, nmrnum
      if (nmrshb(j) .gt. 0) then
        rk2nmr(2, j) = rk2nmr(2, j)*rmult
        rk3nmr(2, j) = rk3nmr(2, j)*rmult
        if (nmrst(3, j) .gt. 0) then
          rk2nmr(1, j) = rk2nmr(1, j)*rmult
          rk3nmr(1, j) = rk3nmr(1, j)*rmult
        end if
      end if
    end do
    weight(10) = wt
    ichang(10) = 1
    if (abs(wt) .le. small) ichold(10) = 2
    if (fixed) ichold(10) = 2
  end if

! TYPE = REST/RESTL:

  if ((itype .eq. 12 .or. itype .eq. 20 .or. ireset(11) .eq. 1) .and. &
       ichold(11) .ne. 2) then
    rmult = wt/weight(11)
    do j = 1, nmrnum
      if (nmrshb(j) .le. 0) then
        rk2nmr(2, j) = rk2nmr(2, j)*rmult
        rk3nmr(2, j) = rk3nmr(2, j)*rmult
        if (nmrst(3, j) .gt. 0) then
          rk2nmr(1, j) = rk2nmr(1, j)*rmult
          rk3nmr(1, j) = rk3nmr(1, j)*rmult
        end if
      end if
    end do

    weight(11) = wt
    ichang(11) = 1
    if (abs(wt) .le. small) ichold(11) = 2
    if (fixed) ichold(11) = 2
  end if

! TYPE = NOESY:

  if ((itype .eq. 21 .or. ireset(12) .eq. 1) .and. ichold(12) .ne. 2) then
    weight(12) = wt
    ichang(12) = 1
    if (abs(wt) .lt. small .or. fixed) ichold(12) = 2
  end if

! TYPE = SHIFTS:

  if ((itype .eq. 22 .or. ireset(13) .eq. 1) .and. ichold(13) .ne. 2) then
    weight(13) = wt
    ichang(13) = 1
    if (abs(wt) .lt. small .or. fixed) ichold(13) = 2
  end if

! TYPE = SOFTR:

  if (itype .eq. 16) then
    rwell = wt
    ichang(15) = 1
  end if


! TYPE = TEMP0:

  if (itype .eq. 13) then
    temp0 = wt
    ichang(16) = 1
  end if

! TYPE = TAUTP:

  if (itype .eq. 23) then
    tautp = wt
    ichang(14) = 1
  end if

! TYPE = CUTOFF:

  ! No longer supported and filtered out on input.  Since this routine does not
  ! do file i/o, we assume we got it on input, and just comment out this block.

! if (itype .eq. 17) then
!   cut = wt**2
!   ichang(17) = 1
! end if

! TYPE = DISAVE, ANGAVE, or TORAVE:

  if (itype .ge. 24 .and. itype .le. 26) then
    indx = itype-23
    tauave(indx) = wt
    ichang(17 + indx) = 1
  end if

! End of implied loop over ichgwt instructions:

71 continue

  i = i + 1

  if (i .le. ichgwt .and. irstyp .eq. 0) goto 20

! Now see what weights were modified. Any that were not modified
! on this call, and that were not set to 1.0 previously, should be reset
! to 1.0.

75 continue

  irstyp = 1
  do 80 i = 1, 14
    if (ireset(i) .eq. 1) then
      goto 80
    else if (ichang(i) .eq. 0 .and. ichold(i) .eq. 1) then
      ireset(i) = 1
      wt = 1.0d0
      wt2 = 1.0d0
      ichold(i) = 0
      goto 25
    else if (ichang(i) .eq. 1 .and. ichold(i) .eq. 0) then
      ichold(i) = 1
    end if
   80 continue

! wnoesy = weight(12)
! wshift = weight(13)
  if (ichang(14) .ne. 1) tautp = weight(14)
  if (ichang(15) .ne. 1) rwell = weight(15)
  if (ichang(16) .ne. 1) temp0 = weight(16)
! cut = weight(17)
  if (ichang(18) .ne. 1) tauave(1) = weight(18)
  if (ichang(19) .ne. 1) tauave(2) = weight(19)
  if (ichang(20) .ne. 1) tauave(3) = weight(20)

  wtls(1) = weight(10)
  wtls(2) = weight(11)

! setgms sets a couple of arrays which depend on the torsional term
! force constants, if those force constants were changed here.

  if (itorfc .eq. 1) call setgms(numphi + nimprp_arg)

  return

end subroutine modwt

!*******************************************************************************
!
! Subroutine:   nmrprt  (NMR PRinT)
!
! Description:
!
! This routine prints the energies/deviations due to the NMR restraints
! to the user-formatted file.
!
! Author: David A. Pearlman
! Date: 7/89
!
! dvdis,dvang,dvtor,eenmr(1,i) are the values on the last call to nmrnrg.
! dvdis,dvang,dvtor,eenmr(2,i) are the accumulated totals over the entire run.
!
!*******************************************************************************

subroutine nmrprt(dvdis, dvang, dvtor, eenmr, nstep, iout)

  use file_io_mod
  use file_io_dat_mod
  use mdin_ctrl_dat_mod
  use nmr_lib_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  double precision  dvdis(2, 4)
  double precision  dvang(2, 4)
  double precision  dvtor(2, 4)
  double precision  eenmr(2, 3)
  integer           nstep
  integer           iout

! Local variables:

  integer               :: j

  if (nmropt .lt. 1) return
  if (nstep .lt. 0) return

  write(iout, 20) (eenmr(1, j), j = 1, 3)
  write(iout, 40)
  return

20 format(' NMR restraints: Bond =',f9.3,3x,'Angle = ',f9.3,3x, &
          'Torsion = ',f9.3)
40 format(79('='))

end subroutine nmrprt

!*******************************************************************************
!
! Subroutine:   ndvprt (Nmr DeViations PRinT)
!
! This routine prints the deviations of each restraint from its
! target value, and the energy contribution for that restraint.
!
! Author: David A. Pearlman
! Date: 7/89
!
! INPUT:
!
!  crd(i)     : Coordinates array.
!  frc(i)     : Force array; Updated during this call.
!  name(i)    : Atom name array
!  irsnam()   : Residue name array
!  nres       : Number of residues in the system
!  natom      : Number of atoms in the system.
!  rimass(i)  : Array of inverse masses
!  nmrnum     : The number of NMR restraints
!  nstep      : The step/iteration number in the main calling program
!  nmrat(16,i) : The 2-4 atoms defining restraint i.
!  kxpk(2,i)  : Integer peak identifier, just for printing 1:peak 2:dataset
!  nmrst(3,i) : Restraint i is imposed over the range of steps
!               nmrst(1,i) -> nmrst(2,i). if nmrst(3,i) .ne. 0, change
!               is implemented as a step function every abs(nmrst(3,i)) steps.
!               Otherwise, changes in target values occur continuously.
!  r1nmr(2,i) : The slope/intercept of the dependence of r1 on step
!               number (see routine disnrg,angnrg or tornrg for description
!               of r1, r2, r3, r4, k2, and k3).
!  r2nmr(2,i) : The slope/int. of the dependence of r2 on step number.
!  r3nmr(2,i) : The slope/int. of the dependence of r3 on step number.
!  r4nmr(2,i) : The slope/int. of the dependence of r4 on step number.
!  rk2nmr(2,i) : The slope/int. of the dependence of k2 on step number.
!  rk3nmr(2,i) : The slope/int. of the dependence of k3 on step number.
!  nmrcom(2,i) : The ranges of atoms defining the center-of-mass group,
!                for applicable distance restraints
!  nmrfty(i)  : Flag for each restraint indicating what functional form should
!               be used for that restraint. If time-averaged restraints
!               are being applied to a particular class of internal
!               (bond, angle, torsion) and nmrfty(i)=1, then an instantaneous
!               restraint will be applied to internal i.
!  igravt(i)  : For averaged group positions (if requested)
!               = 0 for center of mass
!               = 1 for r**-6 average
!  ialtdis(i) : flag to use alternate restraint functional form
!    bave(i)
!   bave0(i)
!    aave(i)
!   aave0(i)    For time-averaged restraints. _ave(i) contains the time-averaged
!   tave1(i)    value of the restraint, using the exponential decay value in
!   tave01(i)   tauave(i). _ave0(i) contains the time-averaged value of the
!   tave2(i)    restraint, using no exponential decay (real average). For
!   tave02(i)   torsions, two values are required for average, corresponding
!               to the cos(tau) and sin(tau)
!
! ajcoef(3,i) : The coefficients to be used in the Karplus equation for
!               any torsional restraints to be imposed using J-coupling.
!
!     nave(3) : For time-averaged restraints. >0 for time-averaged restraints.
!               Elements 1->3 correspond to values for bonds, angles, torsions.
!   ipower(3) : For time-averaged restraints. Gives the exponent used in
!               developing the time average.
!   tauave(3) : For time-averaged restraints. Gives the time constant used
!               in the exponential decay multiplier.
!   ravein(3) : For time-averaged restraints. On the first step, the value
!               returned is the current value+ravein(i) (i=1,3 for bonds,
!               angles, torsions).
!   navint(3) : The averaged value is only recalculated every navint() steps.
!          dt : Integration timestep (ps). Used for time-averaged restraints.
!  iscopn     : Scratch unit for redirected prints (if requested)
!  iout       : Unit for prints (if requested to pout)
!
! Input variables in global storage:
!
!  gbl_res_atms()    : Residue i contains atoms
!                      gbl_res_atms(i) -> gbl_res_atms(i+1)-1
!
!*******************************************************************************

subroutine ndvprt(crd, frc, name, irsnam, nmrnum, nstep, nmrat, kxpk, nmrst, &
                  r1nmr, r2nmr, r3nmr, r4nmr, rk2nmr, rk3nmr, nmrcom, nmrfty, &
                  igravt,  ialtdis, bave, bave0, aave, aave0, tave1, tave01, &
                  tave2, tave02, ajcoef, nave, ipower, tauave, ravein, navint, &
                  iscopn, iout, ifxyz, ioutxyz, idxyz, resttype)

  use file_io_mod
  use file_io_dat_mod
  use mdin_ctrl_dat_mod
  use nmr_lib_mod
  use prmtop_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:

  double precision  crd(*)
  double precision  frc(*)
  character(4)      name(*)
  character(4)      irsnam(*)
  integer           nmrnum
  integer           nstep
  integer           nmrat(16, *) ! change from 4 to 16 to hold CoM  ang/tor
  integer           kxpk(2, *)
  integer           nmrst(3, *)
  double precision  r1nmr(2, *)
  double precision  r2nmr(2, *)
  double precision  r3nmr(2, *)
  double precision  r4nmr(2, *)
  double precision  rk2nmr(2, *)
  double precision  rk3nmr(2, *)
  integer           nmrcom(2, *)
  integer           nmrfty(*)
  integer           igravt(*)
  integer           incflg(3)
  integer           ialtdis(*)
  double precision  bave(*)
  double precision  bave0(*)
  double precision  aave(*)
  double precision  aave0(*)
  double precision  tave1(*)
  double precision  tave01(*)
  double precision  tave2(*)
  double precision  tave02(2)
  double precision  ajcoef(3, *)
  integer           nave(3)
  integer           ipower(3)
  double precision  tauave(3)
  double precision  ravein(3)
  integer           navint(3)
  integer           iscopn
  integer           iout
  integer           ifirst
  integer           ioutxyz(*)
  integer           ifxyz(*)
  integer           resttype(*)
  integer           resttypetmp

  double precision  idxyz(*)
  integer           ii
! Local variables:

  ! Scratch array of length at least natom:

  integer, dimension(natom) :: iscrth

  double precision  bound
  double precision  convrt
  double precision  dev2
  double precision  dfr(12)
  double precision  e
  double precision  eang
  double precision  edis
  double precision  etor
  integer           i
  integer           iat(8)
  integer           iave
  integer           ierr
  integer           itimes
  integer           j
  logical           jcoupl
  logical           usecom
  integer           nstepu
  double precision  r1, r2, r3, r4
  double precision  rint
  double precision  rk2, rk3
  double precision  rmstot(4)
  double precision  small
  double precision  step
  character(14)     tmpnam(4)
  double precision  xcom(12)

  data itimes/0/
  data small/1.0d-7/

  save itimes, small

! First, figure out where this information should be output (based on
! user-defined redirection).

  iuse = iout
    if (iredir(2) .eq. 0) then
      return
    else
      if (redir(2)(1:iredir(2)) .eq. 'POUT') then
        go to 1

! If LISTIN file = LISTOUT file, or if more than one LISTOUT type dump
! is being done, open the LISTOUT file as "old" and position the pointer
! at the bottom of the file:

      else if ((iredir(1) .eq. iredir(2) .and. &
                redir(1)(1:iredir(1)) .eq. redir(2)(1:iredir(2))) .or. &
                itimes .gt. 0) then

        call opnmrg(redir(2)(1:iredir(2)), iscopn, 2, iout, ierr)
        if (ierr .eq. 1) return
      else
        call opnmrg(redir(2)(1:iredir(2)), iscopn, 0, iout, ierr)
        if (ierr .eq. 1) return
      end if
    end if

    write(iout, 9020) redir(2)(1:iredir(2))

    iuse = iscopn

1 continue

  write(iuse, 9032)
  write(iuse, 9077) restrt_name(1:40)
  write(iuse, 9030) pencut
  write(iuse, 9032)
  write(iuse, 9031)
  write(iuse, 9032)

! Set up pointers to the residue number corresponding to each atom in
! the ISCRTH array.

  do i = 1, nres
     do j = gbl_res_atms(i), gbl_res_atms(i + 1) - 1
        iscrth(j) = i
     end do
  end do

! Main loop over all the restraints:

  edis = 0.0
  eang = 0.0
  etor = 0.0
  do 10 i = 1, nmrnum

! Skip restraints if the current step # is outside the applicable range:

    if (nstep .lt. nmrst(1, i) .or. (nstep .gt. nmrst(2, i) .and. &
        nmrst(2, i) .gt. 0)) go to 10

! Calculate the values of r1,r2,r3,r4,k2, and k3 for this step & restraint:

! Vary the step used in calculating the values only every NMRST(3,I) steps.
! If the user did not specify NINC in the input, NMRST(3,I)=1.

    nstepu = nstep - mod(nstep - nmrst(1, i), abs(nmrst(3, i)))
    step = dble(nstepu)
    r1 = r1nmr(1, i)*step + r1nmr(2, i)
    r2 = r2nmr(1, i)*step + r2nmr(2, i)
    r3 = r3nmr(1, i)*step + r3nmr(2, i)
    r4 = r4nmr(1, i)*step + r4nmr(2, i)

! If NMRST(3,I) > 0, then the weights are linearly interpolated.
! If NMRST(3,I) < 0, then the weights are modified by a multplicative factor.

    if (nmrst(3, i) .gt. 0) then
      rk2 = rk2nmr(1, i)*step + rk2nmr(2, i)
      rk3 = rk3nmr(1, i)*step + rk3nmr(2, i)
    else
      nstepu = (nstep - nmrst(1, i)) / abs(nmrst(3, i))
      rk2 = rk2nmr(2, i) * rk2nmr(1, i)**nstepu
      rk3 = rk3nmr(2, i) * rk3nmr(1, i)**nstepu
    end if

! Determine what type of restraint (distance,angle,torsion) this
! is, and call the appropriate routine to calculate its contribution.

    usecom = .false.
    do ii = 9, 16
      usecom = usecom .or. nmrat(ii, i) .lt. 0
    end do

    iave = 0
    jcoupl = .false.
    resttypetmp = resttype(i)
    select case (resttypetmp)
      case (1)  !distances
        if (nave(1) .gt. 0 .and. nmrfty(i) .eq. 0) iave = 1
        if (.not. usecom) then  ! regular distances
          call disnrg(crd, frc, dfr, rint, nmrat(1, i), nmrat(2, i), e, &
                      r1, r2, r3, r4, rk2, rk3, bave(2*i - 1), bave0(i), &
                      ipower(1), tauave(1), ravein(1), dt, navint(1), iave, &
                      0, 0, 2, ialtdis(i), ifxyz(5*(i - 1) + 1), &
                      idxyz(3*(i - 1) + 1))
        else if (igravt(i) .eq. 1) then  ! r6 distances
          !call r6ave(crd, nmrat, nmrcom, rint, i)
          call r6ave(crd, nmrat, nmrcom, rint, i, ifxyz(5*(i - 1) + 1), &
                      idxyz(3*(i - 1) + 1))
          call disnrg(crd, frc, dfr, rint, nmrat(1, i), nmrat(2, i), e, &
                      r1, r2, r3, r4, rk2, rk3, bave(2*i - 1), bave0(i), &
                      ipower(1), tauave(1), ravein(1), dt, navint(1), iave, &
                      0, 0, 3, ialtdis(i), ifxyz(5*(i - 1) + 1), &
                      idxyz(3*(i - 1) + 1))
        else   !com distances
          call nmrcms(crd, xcom, nmrat, nmrcom, rmstot, i)
          call disnrg(xcom, frc, dfr, rint, 0, 3, e, r1, r2, r3, r4, rk2, rk3, &
                      bave(2*i - 1), bave0(i), ipower(1), tauave(1), &
                      ravein(1), dt, navint(1), iave, 0, 0, 2, ialtdis(i), &
                      ifxyz(5*(i - 1) + 1), idxyz(3*(i - 1) + 1))
        end if
        edis = edis + e
      case (2) ! angles
        if (nave(2) .gt. 0 .and. nmrfty(i) .eq. 0) iave = 1
        if (usecom) then  ! com angles
          call nmrcms(crd, xcom, nmrat, nmrcom, rmstot, i)
          call angnrg(xcom, frc, dfr, rint, 0, 3, 6, e,                       &
                    r1, r2, r3, r4, rk2, rk3, aave(2*i - 1), aave0(i), ipower(2), &
                    tauave(2), ravein(2), dt, navint(2), iave, 0, 0, 2)
        else  ! regular angles
          call angnrg(crd, frc, dfr, rint, nmrat(1, i), nmrat(2, i), nmrat(3, i), e, &
                    r1, r2, r3, r4, rk2, rk3, aave(2*i - 1), aave0(i), ipower(2), &
                    tauave(2), ravein(2), dt, navint(2), iave, 0, 0, 2)
        end if
        eang = eang + e
      case (3) ! torsions and j coupling
        if (nave(3) .gt. 0 .and. nmrfty(i) .eq. 0) iave = 1
        ! If J coefficients are 0 for this torsion, it is a normal torsional
        ! restraint. If they are .NE.0, then it is a J-coupling factor restraint.
        jcoupl = abs(ajcoef(1, i)) .gt. small .or. abs(ajcoef(2, i)) .gt. small
        if (.not.jcoupl) then
          if (usecom) then ! COM torsions
            call nmrcms(crd, xcom, nmrat, nmrcom, rmstot, i)
            call tornrg(xcom, frc, dfr, rint, 0, 3, 6, 9,   &
                       e, r1, r2, r3, r4, rk2, rk3,   &
                       tave1(2*i - 1), tave01(i), tave2(2*i - 1), &
                       tave02(i), ipower(3),  &
                       tauave(3), ravein(3), dt, navint(3), iave, 0, 0, 2)
          else ! regular torsions
            call tornrg(crd, frc, dfr, rint, nmrat(1, i), nmrat(2, i), nmrat(3, i), &
                       nmrat(4, i), e, r1, r2, r3, r4, rk2, rk3, &
                       tave1(2*i - 1), tave01(i), tave2(2*i - 1), &
                       tave02(i), ipower(3), &
                       tauave(3), ravein(3), dt, navint(3), iave, 0, 0, 2)
          end if
        else
          if (usecom) then
            ! Since atom groups are not supported for j-coupling restraints,
            ! print out an error message and terminate execution if the user
            ! specifies atom groups AND j-coupling
            write(iout, 9086)
            call mexit(iout, 1)
          else ! jcoupling
            call jnrg(crd, frc, rint, nmrat(1, i), nmrat(2, i), nmrat(3, i), &
                     nmrat(4, i), e, r1, r2, r3, r4, rk2, rk3, &
                     tave1(2*i - 1), tave01(i), ipower(3), &
                     tauave(3), ravein(3), dt, navint(3), iave, 0, 0, &
                     ajcoef(1, i), 2)
          end if
        end if
        etor = etor + e
        ! case 4 5 6
      case default
        write(iout, 9087)
        call mexit(iout, 1)
    end select


! If energy of this restraint is less than the cutoff, do not print anything

    if (e .lt. pencut) go to 10

! Calculate the deviation of the current value of the internal from its target

      if (rint .ge. r2 .and. rint .le. r3) then
        dev2 = 0.0d0
! if no penalty, print upper bound:
        bound = r3
      else if (rint .lt. r2) then
        dev2 = r2 - rint
        bound = r2
      else
        dev2 = rint - r3
        bound = r3
      end if

! Print the internal, its value, its deviations from the target,
! and its energy contribution:

      do ii = 1, 8
        iat(ii) = (nmrat(ii, i) / 3) + 1
        if (nmrat(ii, i) .lt. 0 .and. nmrat(ii+8, i) .lt. 0) &
          iat(ii) = nmrcom(2, -nmrat(ii+8, i))
      end do

      usecom = .false.

! convention for printing: if a group is defined, identify it
!        by its final atom, and put a "*" in the first column.
!        This is usually sufficient, e.g. if
!        the group is the three atoms of a methyl group, or two
!        equivalent delta or epsilon protons on an aromatic ring,
!        etc.  For more exotic groupings of atoms, this will not
!        be as informative as it could be, but we will cross that
!        bridge if it turns out to be a real problem.

!        Also add a one-letter code at the end of each line:
!              d       for distance constraint
!              a       for angle constraint
!              t       for torsion constraint
!              j       for j-coupling constraint


      convrt = 1.0d0
      ! 6 for when we port generalized distance coordinate
      if (resttype(i) .ne. 1 .and. resttype(i) .ne. 6 .and. &
            .not.jcoupl) convrt  = 180.0d0 / PI


! for torsions, print the two central atoms and the deviations:
!       (for J-couplings, print the two end atoms, since this is generally
!        more informative.)



      resttypetmp = resttype(i)
      select case (resttypetmp)
        case (1)
          write(tmpnam(1)(2:14), '(A4,1x,A4,I4)') &
              name(iat(1)), irsnam(iscrth(iat(1))), iscrth(iat(1))
          write(tmpnam(2)(2:14), '(A4,1x,A4,I4)') &
              name(iat(2)), irsnam(iscrth(iat(2))), iscrth(iat(2))
          tmpnam(1)(1:1) = ' '
          tmpnam(2)(1:1) = ' '
          if (nmrat(1, i) .lt. 0) tmpnam(1)(1:1) = '*'
          if (nmrat(2, i) .lt. 0) tmpnam(2)(1:1) = '*'

          write(iuse, 9073) tmpnam(1), tmpnam(2), &
                rint, bound, abs(dev2), e, kxpk(1, i), kxpk(2, i)

        case (2)
          write(tmpnam(1)(2:14), '(A4,1x,A4,I4)') &
               name(iat(1)), irsnam(iscrth(iat(1))), iscrth(iat(1))
          write(tmpnam(3)(2:14), '(A4,1x,A4,I4)') &
               name(iat(3)), irsnam(iscrth(iat(3))), iscrth(iat(3))
          tmpnam(1)(1:1) = ' '
          tmpnam(3)(1:1) = ' '
          if (nmrat(1,i) .lt. 0) tmpnam(1)(1:1) = '*'
          if (nmrat(3,i) .lt. 0) tmpnam(3)(1:1) = '*'
          write(iuse, 9041) tmpnam(1), tmpnam(3), &
                convrt*rint, convrt*bound, convrt*abs(dev2), e

        case (3)
          ! for torsions, print the two central atoms and the deviations:
          !       (for J-couplings, print the two end atoms, since this is generally
          !        more informative.)
          write(tmpnam(1)(2:14), '(A4,1x,A4,I4)') &
                name(iat(1)), irsnam(iscrth(iat(1))), iscrth(iat(1))
          write(tmpnam(2)(2:14), '(A4,1x,A4,I4)') &
                name(iat(2)), irsnam(iscrth(iat(2))), iscrth(iat(2))
          write(tmpnam(3)(2:14), '(A4,1x,A4,I4)') &
                name(iat(3)), irsnam(iscrth(iat(3))), iscrth(iat(3))
          write(tmpnam(4)(2:14), '(A4,1x,A4,I4)') &
                name(iat(4)), irsnam(iscrth(iat(4))), iscrth(iat(4))
          tmpnam(1)(1:1) = ' '
          tmpnam(2)(1:1) = ' '
          tmpnam(3)(1:1) = ' '
          tmpnam(4)(1:1) = ' '
          if (nmrat(1,i) .lt. 0) tmpnam(1)(1:1) = '*'
          if (nmrat(2,i) .lt. 0) tmpnam(2)(1:1) = '*'
          if (nmrat(3,i) .lt. 0) tmpnam(3)(1:1) = '*'
          if (nmrat(4,i) .lt. 0) tmpnam(4)(1:1) = '*'
          if (jcoupl) then
            write(iuse, 9082) tmpnam(1), tmpnam(4), &
                  convrt*rint, convrt*bound, convrt*abs(dev2), e
          else
            write(iuse, 9083) tmpnam(2), tmpnam(3), &
                  convrt*rint, convrt*bound, convrt*abs(dev2), e
          end if
      end select

! If a group-type distance restraint was defined, echo the sorted group:
10 continue

  if (edis .gt. 0.0) write(iuse, 9074) edis
  if (eang .gt. 0.0) write(iuse, 9075) eang
  if (etor .gt. 0.0) write(iuse, 9076) etor
#ifndef MPI
# ifndef CUDA
  !Averaging not working for nmropt on CUDA at present.
  write(iuse, 9084) ebdev
  write(iuse, 9085) eadev
# endif
#endif
  write(iuse, 9032)

  itimes = itimes + 1

  if (iuse .eq. iscopn) call opnmrg(redir(1), iuse, 3, iout, ierr)

  return

9020 format(' Restraints/deviations being written to file: ', a)
9030 format(' Restraints, deviations, and energy contributions:', &
            '    pencut = ', f7.2 /)
9031 format('     First atom        Last atom    curr. value target', &
            ' deviation  penalty')
9032 format(' ', 78('-'))
9041 format(' ', a14, ' -- ', a14, ':', 4f9.3, ' a')
9073 format(' ', a14, ' -- ', a14, ':', 4f9.3, ' d', i5, ':', i2)
9074 format(39x, 'Total distance penalty: ', f10.3)
9075 format(39x, 'Total angle    penalty: ', f10.3)
9076 format(39x, 'Total torsion  penalty: ', f10.3)
9077 format(/ / ' Final Restraint Analysis for coords: ', a40 / /)
9082 format(' ', a14, ' -- ', a14, ':', 4f9.3, ' j')
9083 format(' ', a14, ' -- ', a14, ':', 4f9.3, ' t')
#ifndef MPI
9084 format('| ', 30x, 'RMS deviation from ideal bonds : ', f11.4)
9085 format('| ', 30x, 'RMS deviation from ideal angles: ', f10.3)
#endif
9086 format('Error: Atom groupings not supported for j-coupling restraints.')
9087 format('Error: Improper restraint specified.')

end subroutine ndvprt

!*******************************************************************************
! Subroutine:   cuda_nmr_setup()
!
! Description:
!
! This routine passesa arrays to the gpu code - gpu_nmr_setup_ - so forces
!  and energies can all be evaluated and coordinates integrated on the device.
!
! Author: Ross Walker .and. Scott Le Grand
! Date: 10/21/2015
#ifdef CUDA
subroutine cuda_nmr_setup()
use mdin_ctrl_dat_mod, only : jar
    !added the restraint type array to transfer the restraint type to the GPU
    call gpu_nmr_setup(nmrnum, nmr_iwork(jnmrat), &
              nmr_iwork(jnmrst), nmr_work(j1nmr), nmr_work(j2nmr), &
              nmr_work(j3nmr), nmr_work(j4nmr), nmr_work(jk2nmr), &
              nmr_work(jk3nmr), jar, drjar, nmr_iwork(jcomdf), maxgrp,&
              nmr_iwork(jgravt), nmr_iwork(jfxyz), nmr_work(jdxyz), &
              nmr_iwork(jrsttyp))
end subroutine cuda_nmr_setup
#endif

!*******************************************************************************
!
! Subroutine:   nmrnrg (NMR eNeRGy)
!
! Description:
!
! This routine calculates the contributions
! to the energy/forces due to distance/valence angle/torsion restraints.
!
! Author: David A. Pearlman
! Date: 7/89
!
! INPUT:
!
!  crd(i)     : Coordinate array
!  frc(i)     : Force array; Updated during this call.
!  rimass(i)  : Array of inverse masses
!  nmrnum     : The number of NMR restraints
!  nstep      : The step/iteration number in the main calling program
!  nmrat(16,i) : The 2-6 atoms defining restraint i.
!                           NMRAT(9,I)-NMRAT(16,I) are used for atom groups
!  nmrst(3,i) : Restraint i is imposed over the range of steps
!               nmrst(1,i) -> nmrst(2,i). if nmrst(3,i) .ne. 0, change
!               is implemented as a step function every abs(nmrst(3,i)) steps.
!               Otherwise, changes in target values occur continuously.
!  r1nmr(2,i) : The slope/intercept of the dependence of r1 on step
!               number (see routine disnrg,angnrg or tornrg for description
!               of r1, r2, r3, r4, k2, and k3).
!  r2nmr(2,i) : The slope/int. of the dependence of r2 on step number.
!  r3nmr(2,i) : The slope/int. of the dependence of r3 on step number.
!  r4nmr(2,i) : The slope/int. of the dependence of r4 on step number.
!  rk2nmr(2,i) : The slope/int. of the dependence of k2 on step number.
!  rk3nmr(2,i) : The slope/int. of the dependence of k3 on step number.
!  nmrcom(2,i) : The ranges of atoms defining the center-of-mass group,
!                for applicable distance restraints
!  nmrfty(i)  : Flag for each restraint indicating what functional form should
!               be used for that restraint. If time-averaged restraints
!               are being applied to a particular class of internal
!               (bond, angle, torsion) and nmrfty(i)=1, then an instantaneous
!               restraint will be applied to internal i.
!  igravt(i)  : For averaged group restraint positions.
!               = 0, use center of mass
!               = 1, use r**-6 averaged position.
!  ialtdis(i) : Flag to use alternative restraint functional form
!    bave(i)
!   bave0(i)
!    aave(i)
!   aave0(i)    For time-averaged restraints. _ave(i) contains the time-averaged
!   tave1(i)    value of the restraint, using the exponential decay value in
!   tave01(i)   tauave(i). _ave0(i) contains the time-averaged value of the
!   tave2(i)    restraint, using no exponential decay (real average). For
!   tave02(i)   torsions, two values are required for average, corresponding
!               to the cos(tau) and sin(tau)
! ajcoef(3,i) : The coefficients to be used in the Karplus equation for
!               any torsional restraints to be imposed using J-coupling.
!     nave(3) : If time-averaged restraint have been requested, nave(i)>0.
!               Elements 1->3 correspond to values for bonds, angles, torsions.
!   ipower(3) : For time-averaged restraints. Gives the exponent used in
!               developing the time average.
!   tauave(3) : For time-averaged restraints. Gives the time constant used
!               in the exponential decay multiplier for the energies/forces
!               reported back to the calling program.
!   ravein(3) : For time-averaged restraints. On the first step, the value
!               returned is the current value+ravein(i) (i=1,3 for bonds,
!               angles, torsions).
!          dt : Integration timestep (ps). Used for time-averaged restraints.
!      navint : Time-averaged restraints are only updated every navint
!               steps. Typically, navint=1.
!   iavtyp(3) : Determines how forces will be calculated when time-averaged
!               retraints are used. iavtyp(i) = 1 --> Calculate forces
!               according to the standard energy expression. iavtyp(i) = 2 -->
!               calculate forces as dE/d(r_ave) dr(t)/d(x). Latter form
!               integrates into a non-intuitive pseudo-energy, but avoids
!               large forces from (1+ipower) exponentiation.
!      idmpav : Values of retraints will be dumped to the file redir(6) every
!               idmpav steps, if idmpav > 0. A time-averaged value for a
!               restraint will be reported, if time-averaging is being performed
!               for that retraint.
!      idumpu : Unit to be used for restraint dumps if idmpav > 0. This
!               unit must be dedicated for the entire run.
!
! OUTPUT:
!
! enmr(1)     : Total energy contribution from distance restraints.
! enmr(2)     : Total energy contribution from angle restraints.
! enmr(3)     : Total energy contribution from torsion restraints.
! devdis(1): The average deviation of the distance restraints. The
!            average of r2 and r3 is used as the "target" distance.
!       (2): The rms deviation of the distance restraints. The average
!            of r2 and r3 is used as the "target" distance.
!       (3): Average deviation, but a deviation of 0 is used for restraints
!            falling between r2 and r3, and r2 or r3 (whichever is
!            closer) is used outside this range.
!       (4): The rms deviation, but same definitions of deviations as (3).
! devang(3): same as devdis, but for angles
! devtor(3): same as devdis, but for torsions
!
!*******************************************************************************

subroutine nmrnrg(crd, frc, nmrnum, nstep, nmrat, nmrst, r1nmr, r2nmr, r3nmr, &
                  r4nmr, rk2nmr, rk3nmr, nmrcom, nmrfty, igravt, ialtdis, &
                  bave, bave0, aave, aave0, tave1, tave01, tave2, tave02, &
                  ajcoef, nave, ipower, tauave, ravein,  navint, &
                  iavtyp, idmpav, idumpu, enmr, devdis, devang, devtor, iout, ifxyz, &
                  ioutxyz, idxyz, resttype)

  use file_io_mod
  use file_io_dat_mod
  use mdin_ctrl_dat_mod
  use nmr_lib_mod
  use prmtop_dat_mod
  use parallel_dat_mod
  use ti_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:

  double precision  crd(*)
  double precision  frc(*)
#ifdef CUDA
#ifdef MPI
  double precision  cmsfrc(3*natom)
#endif
#endif
  integer           nmrnum
  integer           nstep
  integer           ioutxyz(*)
  integer           nmrat(16, *)
  integer           nmrst(3, *)
  double precision  r1nmr(2, *)
  double precision  r2nmr(2, *)
  double precision  r3nmr(2, *)
  double precision  r4nmr(2, *)
  double precision  rk2nmr(2, *)
  double precision  rk3nmr(2, *)
  integer           nmrcom(2, *)
  integer           nmrfty(*)
  integer           igravt(*)
  integer           ialtdis(*)
  double precision  bave(*)
  double precision  bave0(*)
  double precision  aave(*)
  double precision  aave0(*)
  double precision  tave1(*)
  double precision  tave01(*)
  double precision  tave2(*)
  double precision  tave02(*)
  double precision  ajcoef(3, *)
  integer           nave(3)
  integer           ipower(3)
  double precision  tauave(3)
  double precision  ravein(3)
  integer           navint(3)
  integer           iavtyp(3)
  integer           idmpav
  integer           resttype(*) ! added for case
  integer           resttypetmp !added for case implementation
  integer           idumpu
  double precision  enmr(3)
  double precision  devdis(4)
  double precision  devang(4)
  double precision  devtor(4)
  integer           iout
  integer           ifxyz(*)
  double precision  idxyz(*)

! Local variables:

  double precision  convrt
  double precision  dev
  double precision  dev2
  double precision  deviat(3, 4)
  double precision  dfr(12)
  double precision  e
  integer           i
  integer           iave
  integer           ifirst
  integer           incflg(3)
  integer           inum(3)
  integer           j
  logical           jcoupl
  logical           usecom
  integer           jj
  integer           ii
  integer           jtyp
  integer           nstepu
  double precision  r1
  double precision  r2
  double precision  r3
  double precision  r4
  double precision  rint
  double precision  rk2
  double precision  rk3
  double precision  rmstot(4)
  double precision  rnum
  double precision  step
  double precision  target
  double precision  xcom(12)

  integer, save                 :: ionce = 0
  double precision, save        :: small = 1.d-7

  ! For jar:

  double precision              :: fcurr
  double precision, save        :: fold = 0.d0
  double precision, save        :: work = 0.d0
  integer loop
  ! For CUDA MPI
#ifdef CUDA
#ifdef MPI
  integer         buf_size
#endif
#endif
#ifdef CUDA
    if (.not. master) return
    if (idmpav .le. 0) return
    if (mod(nstep,idmpav) .ne. 0) return
    call gpu_download_crd(crd)
#endif

  ifirst = 0
  if (ionce .ne. 9899) then
    ifirst = 1
    ionce = 9899
  end if

! If incremental "dumps" of restraint values have been requested,
! and this is the first call to nmrnrg, open the appropriate file:

  if (master .and. ifirst .eq. 1 .and. iredir(6) .gt. 0) &
    call amopen(idumpu, redir(6)(1:iredir(6)), owrite, 'F', 'W')

! Zero the accumulators:

  do i = 1, 3
    do j = 1, 4
      deviat(i, j) = 0.0d0
    end do
    inum(i) = 0
    incflg(i) = 1
    enmr(i) = 0.0d0
  end do
  jcoupl = .false.

! Main loop over all the restraints:
  do i = 1, nmrnum

! Skip restraints if the current step # is outside the applicable range:

    rint = 0.0d0
    if (nstep .lt. nmrst(1, i) .or. (nstep .gt. nmrst(2, i) .and. &
        nmrst(2, i) .gt. 0)) goto 12

! Determine if the restraint includes softcore atoms and add to energy to the
! correct array if it does. Note: If the new restraints in sander are added
! to pmemd, this routine will need to be updated too.
    if (ti_mode .ne. 0) then
       call ti_check_res(nmrat, nmrcom, i)
       ti_lscale = .false.
       if (ti_mode .ne. 1 .and. res_ti_region .ne. 0) then
         call ti_check_res_lscale(nmrat, nmrcom, i)
       end if
       if (emil_sc_lcl .eq. 1) ti_lscale = .true.
    end if

! Calculate the values of r1,r2,r3,r4,k2, and k3 for this step & restraint:

! Vary the step used in calculating the values only every abs(nmrst(3,i)) steps.
! If the user did not specify ninc in the input, abs(nmrst(3,i))=1.

    nstepu = nstep - mod(nstep - nmrst(1, i), abs(nmrst(3, i)))
    step = dble(nstepu)
    r1 = r1nmr(1, i)*step + r1nmr(2, i)
    r2 = r2nmr(1, i)*step + r2nmr(2, i)
    r3 = r3nmr(1, i)*step + r3nmr(2, i)
    r4 = r4nmr(1, i)*step + r4nmr(2, i)

! If nmrst(3,i) > 0, then the weights are linearly interpolated.
! If nmrst(3,i) < 0, then the weights are modified by a multplicative factor.

    if (nmrst(3, i) .gt. 0) then
      rk2 = rk2nmr(1, i)*step + rk2nmr(2, i)
      rk3 = rk3nmr(1, i)*step + rk3nmr(2, i)
    else
      nstepu = (nstep - nmrst(1, i))/abs(nmrst(3, i))
      rk2 = rk2nmr(2, i) * rk2nmr(1, i)**nstepu
      rk3 = rk3nmr(2, i) * rk3nmr(1, i)**nstepu
    end if

! Take the average of r2 & r3 as the average "target" value:

    target = (r2 + r3)/2.0d0

! Determine what type of restraint (distance,angle,torsion) this
! is, and call the appropriate routine to calculate its contribution.
    usecom = .false.
    do ii = 9, 16
      usecom = usecom .or. nmrat(ii, i) .lt. 0
    end do

    iave = 0
    jcoupl = .false.
    resttypetmp = resttype(i)
    select case (resttypetmp)
      case (1) ! Distance restrain
        if (nave(1) .gt. 0 .and. nmrfty(i) .eq. 0) iave = iavtyp(1)
        if (.not. usecom) then  ! regular distance restraints
          call disnrg(crd, frc, dfr, rint, nmrat(1, i), nmrat(2, i), e, r1, &
                      r2, r3, r4, rk2, rk3, bave(2*i - 1), bave0(i), ipower(1), &
                      tauave(1), ravein(1), dt, navint(1), iave, &
                      incflg(1), ifirst, 0, ialtdis(i), ifxyz(5*(i - 1) + 1), &
                      idxyz(3*(i - 1) + 1))
        else if (igravt(i) .eq. 1) then ! r6 averaging distances
          !call r6ave(crd, nmrat, nmrcom, rint, i)
          call r6ave(crd, nmrat, nmrcom, rint, i, ifxyz(5*(i - 1) + 1), &
                      idxyz(3*(i - 1) + 1))
          call disnrg(crd, frc, dfr, rint, nmrat(1, i), nmrat(2, i), e, r1, &
                      r2, r3, r4, rk2, rk3, bave(2*i - 1), &
                      bave0(i), ipower(1), tauave(1), ravein(1), dt, &
                      navint(1), iave, incflg(1), ifirst, 3, ialtdis(i), ifxyz(5*(i - 1) + 1), &
                      idxyz(3*(i - 1) + 1))
          !call r6drv(frc, dfr(1))
          call r6drv(frc, dfr(1), ifxyz(5*(i - 1) + 1), idxyz(3*(i - 1) + 1))
        else  !COM distances
          call nmrcms(crd, xcom, nmrat, nmrcom, rmstot, i)
          call disnrg(xcom, frc, dfr, rint, 0, 3, e, r1, r2, r3, &
                    r4, rk2, rk3, bave(2*i - 1), bave0(i), &
                    ipower(1), tauave(1), ravein(1), dt, navint(1), &
                    iave, incflg(1), ifirst, 2, ialtdis(i), ifxyz(5*(i - 1) + 1), &
                    idxyz(3*(i - 1) + 1))
          call nmrcmf(frc, dfr, nmrat, nmrcom, rmstot, i)
        end if
        jtyp = 1
      case (2)
        if (nave(2) .gt. 0 .and. nmrfty(i) .eq. 0) iave = iavtyp(2)
        if (usecom) then ! com angles
           call nmrcms(crd, xcom, nmrat, nmrcom, rmstot, i)
           call angnrg(xcom, frc, dfr, rint, 0, 3, 6,         &
                      e, r1, r2, r3, r4, rk2, rk3, aave(2*i - 1),   &
                      aave0(i), ipower(2), tauave(2), ravein(2), dt, navint(2),  &
                      iave, incflg(2), ifirst, 3)
           call nmrcmf(frc, dfr, nmrat, nmrcom, rmstot, i)
        else ! regular angles
           call angnrg(crd, frc, dfr, rint, nmrat(1, i), nmrat(2, i), nmrat(3, i), &
                       e, r1, r2, r3, r4, rk2, rk3, aave(2*i - 1), &
                       aave0(i), ipower(2), tauave(2), ravein(2), dt, navint(2), &
                       iave, incflg(2), ifirst, 0)
        end if
        jtyp = 2
      case (3)
        if (nave(3) .gt. 0 .and. nmrfty(i) .eq. 0) iave = iavtyp(3)

! If J coefficients are 0 for this torsion, it is a normal torsional
! restraint. If they are .ne.0, then it is a J-coupling factor restraint.
        jcoupl = abs(ajcoef(1, i)) .gt. small .or. abs(ajcoef(2, i)) .gt. small
        if (.not.jcoupl) then
          if (usecom) then ! com torsion
            call nmrcms(crd, xcom, nmrat, nmrcom, rmstot, i)
            call tornrg(xcom, frc, dfr, rint, 0, 3, 6, 9,            &
                        e, r1, r2, r3, r4, rk2, rk3,                 &
                        tave1(2*i - 1), tave01(i), tave2(2*i - 1),   &
                        tave02(i), ipower(3), tauave(3), ravein(3), dt, &
                        navint(3), iave, incflg(3), ifirst, 3)
            call nmrcmf(frc, dfr, nmrat, nmrcom, rmstot, i)
          else ! regular torsion
            call tornrg(crd, frc, dfr, rint, nmrat(1, i), nmrat(2, i), nmrat(3, i), &
                    nmrat(4, i), e, r1, r2, r3, r4, rk2, rk3, &
                    tave1(2*i - 1), tave01(i), tave2(2*i - 1), &
                    tave02(i), ipower(3), tauave(3), ravein(3), dt, &
                    navint(3), iave, incflg(3), ifirst, 0)
          end if
        else
          if (usecom) then ! jcoupling
             ! Since atom groups are not supported for j-coupling restraints,
             ! print out an error message and terminate execution if the user
             ! specifies atom groups AND j-coupling
             ! This exit condition shouldn't be used, since it should have been
             ! checked in nmrred.  However, I have it here just for safety.

             write(iout,9074)
             call mexit(iout, 1)
          else
              call   jnrg(crd, frc, rint, nmrat(1, i), nmrat(2, i), nmrat(3, i), &
                         nmrat(4, i), e, r1, r2, r3, r4, rk2, rk3, &
                         tave1(2*i - 1), tave01(i), ipower(3), &
                         tauave(3), ravein(3), dt, navint(3), iave, &
                         incflg(3), ifirst, ajcoef(1, i), 0)
          end if
        end if
        jtyp = 3
      case default
        ! Improper restraint specified.  This should have also been checked by nmrred,
        ! but again I have the exit condition here just in case.
        write(iout,9075)
        call mexit(iout, 1)
    end select
    !kelb
    !write(*,*) 'line4194--nmr_calls.F90--nmrat=', nmrat(1:16,i)
    if (ti_mode .eq. 0) then
      enmr(jtyp) = enmr(jtyp) + e
    else if (ti_mode .eq. 1 .or. ti_lscale) then
      if (res_ti_region .eq. 0) then
        enmr(jtyp) = enmr(jtyp) + e
      else
        enmr(jtyp) = enmr(jtyp) + e * ti_weights(res_ti_region)
#if !defined(CUDA) || !defined(GTI)
!because dvdl is dealt with in a weird way from GPU to CPU
!it's easier to just not calculate the dvdl on GPU

!TL: GTI implementation keeps the CPU-way
        call ti_update_ene(e, si_restraint_ene, res_ti_region, 1)
#endif
      end if
    else
      if(res_ti_region .eq. 0) then
        enmr(jtyp) = enmr(jtyp) + e
      else if(res_ti_region .eq. 1) then !jtyp = 1,2,3 for dist,ang,tors
        ti_ene_aug(1,ti_rest_dist_ene - 1 + jtyp) = &
          ti_ene_aug(1,ti_rest_dist_ene - 1 + jtyp) + e
      else if(res_ti_region .eq. 2) then
        ti_ene_aug(2,ti_rest_dist_ene - 1 + jtyp) = &
          ti_ene_aug(2,ti_rest_dist_ene - 1 + jtyp) + e
      end if
   end if

! Sum the appropriate quantities into the deviation accumulators.
! deviat(1,1->4) stores bonds; (2,1->4) stores angles; (3,1->4) stores torsions
! If this was a coupling constant restraint (jcoupl= .true. ), scale the
! deviation by 3.14/180. here, as it will be scaled by 180./3.14 before
! printing...

    dev = abs(target - rint)

    if (rint .ge. r2 .and. rint .le. r3) then
      dev2 = 0.0d0
    else if (rint .lt. r2) then
      dev2 = r2 - rint
    else
      dev2 = rint - r3
    end if

    if (jcoupl) dev  = dev  * (PI/180.0d0)
    if (jcoupl) dev2 = dev2 * (PI/180.0d0)

    deviat(jtyp, 1) = deviat(jtyp, 1) + dev
    deviat(jtyp, 2) = deviat(jtyp, 2) + dev**2
    deviat(jtyp, 3) = deviat(jtyp, 3) + dev2
    deviat(jtyp, 4) = deviat(jtyp, 4) + dev2**2
    inum(jtyp) = inum(jtyp) + 1

! Set incflg to zero after the first time-averaged restraint of a
! particular type has been stored:
    if (iave .gt. 0) incflg(resttype(i)) = 0

! If requested, write out restraint values, all on a single line:

12 continue

    ! Adding jar option to write correctly (roit. 02/27/05)
    ! if jar=0, the behaviour is as usual, if jar=1,
    ! writting is done in the same file, but in the way:
    ! r2(t), rint(t), force(t), accumulated work(t)

    if (jar .eq. 0 .and. .not. skip_print) then

      if (master .and. idmpav .gt. 0) then
         if (mod(nstep, idmpav) .eq. 0) then
            if (i .eq. 1) write(idumpu, '(i8,1x)', advance='no') nstep
            convrt = 1.d0
            if (ioutxyz(i) .eq. 1) then
                write(idumpu, '(4x,A,f9.3,4x)', advance='no')'x:', idxyz(3*(i-1)+1) * convrt
                write(idumpu, '(A,1x,f9.3,4x)', advance='no')'y:', idxyz(3*(i-1)+2) * convrt
                write(idumpu, '(A,1x,f9.3,4x)', advance='no')'z:', idxyz(3*(i-1)+3) * convrt
            end if
            if (resttype(i) .ne. 1 .and. .not. jcoupl) convrt = RAD_TO_DEG
            if(ifxyz(5*(i-1)+1) .ne. 1 .or. ifxyz(5*(i-1)+2) .ne. 1 .or. &
               ifxyz(5*(i-1)+3) .ne. 1) then
              write(idumpu, '(f9.3,1x)', advance='no')sqrt(idxyz(3*(i-1)+1)**2 &
                + idxyz(3*(i-1)+2)**2 + idxyz(3*(i-1)+3)**2)
            else
              write(idumpu, '(f9.3,1x)', advance='no') rint * convrt
            end if
         end if
      end if

    else if (jar .eq. 1 .and. .not. skip_print) then

      if (master .and. idmpav .gt. 0) then

#ifdef CUDA
         if (mod(nstep, idmpav) .eq. 0) then
#endif
            convrt = 1.d0
            if (resttype(i) .ne. 1 .and. .not. jcoupl) convrt = RAD_TO_DEG

#ifdef CUDA
            call gpu_nmr_get_jar_data(r2, rint, fcurr, work)
#else
               ! Set jar variables to the initial values (roit. 02/27/05)

            if (ifirst .eq. 1) then
                work = 0.d0
                fold = 2.d0 * rk2 * (rint - r2)
            endif

            ! Calculate the jar force and work (roit. 02/27/05)

            fcurr = -2.d0 * rk2 * (rint - r2)
            work = work + (fcurr + fold) * drjar / 2.d0
            fold = fcurr
         if (mod(nstep, idmpav) .eq. 0) then
#endif


            write(idumpu, '(f12.5,2x)', advance='no') r2 * convrt
            write(idumpu, '(f12.5,2x)', advance='no') rint * convrt
            write(idumpu, '(f12.5,2x)', advance='no') fcurr / convrt
            write(idumpu, '(f12.5,2x)', advance='no') work / convrt
         end if !mod(nstep,...)
      endif  !master
    end if !jar

  end do !end of loop over all restraints

  ! End the dump record:

  if (master .and. idmpav .gt. 0 .and. .not. skip_print) then
    if (mod(nstep, idmpav) .eq. 0) write(idumpu, '(a)')
  end if

  ! skip_print will only skip one record. After this, it's turned back off

  skip_print = .false.

  ! Determine the RMS averages and return

  do jtyp = 1, 3
    if (inum(jtyp) .gt. 0) then
      rnum = dble(inum(jtyp))
      deviat(jtyp, 1) = deviat(jtyp, 1)/rnum
      deviat(jtyp, 2) = sqrt(deviat(jtyp, 2)/rnum)
      deviat(jtyp, 3) = deviat(jtyp, 3)/rnum
      deviat(jtyp, 4) = sqrt(deviat(jtyp, 4)/rnum)
    end if
  end do

  devdis(1) = deviat(1, 1)
  devdis(2) = deviat(1, 2)
  devdis(3) = deviat(1, 3)
  devdis(4) = deviat(1, 4)
  devang(1) = deviat(2, 1)
  devang(2) = deviat(2, 2)
  devang(3) = deviat(2, 3)
  devang(4) = deviat(2, 4)
  devtor(1) = deviat(3, 1)
  devtor(2) = deviat(3, 2)
  devtor(3) = deviat(3, 3)
  devtor(4) = deviat(3, 4)

  return

9074 format('Error: Atom groupings not supported for j-coupling restraints.')
9075 format('Error: Improper restraint specified.')

end subroutine nmrnrg

!*******************************************************************************
!
! Subroutine:   set_num_improp_dihed  (set improper number)
!
! Description:
!
! This simple subroutine looks through the lp() array for improper torsions.
! The number of improper torsional parameters (not torsions themselves)
! is returned. This value is used in modifying the relative
! weight of regular torsion terms vs. improper torsion terms.
!
! Author: David A. Pearlman
! Date: 8/89
!
! INPUT
!
! dihedh(i) and diheda(i), containing:
!
! dihedh(i)%atm_l : Array containing the fourth atom of torsion i. For torsions
!                   with hydrogens at one/both ends. If < 0, then this is an
!                   improper torsion.
! diheda(i)%atm_l : Array containing fourth atom of torsion i. For torsions with
!                   no hydrogens at ends. For perturbation calculations
!                   (nfpert>0), these are followed by 2*nfpert perturbed
!                   torsions.
! dihedh(i)%parm_idx: Pointer into the force constant array for torsion i, for
!                     torsions with hydrogens at one/both ends.
! diheda(i)%parm_idx: Pointer corresponding to the torsions of diheda(i)%atm_l.
! nphih:   Number of torsion angles with hydrogens at one/both ends.
! nphia:   Number of torsion angles without hydrogens at the ends.
! numphi: Total number of torsional parameters (regular+improper)
!
! OUTPUT
!
! nimprp: The number of improper torsional parameters. the lowest
!                   value icp(i) corresponding to an improper is used to
!                   calculate this value, assuming improper torsion parameters
!                   are packed after the regular torsions.
!*******************************************************************************

subroutine set_num_improp_dihed(nphih, dihedh, nphia, diheda, numphi)

  use gbl_datatypes_mod

  implicit none

! Formal arguments:

  integer               :: nphih
  type(dihed_rec)       :: dihedh(*)
  integer               :: nphia
  type(dihed_rec)       :: diheda(*)
  integer               :: numphi

! Local variables:

  integer           i

  nimprp = numphi + 1

  do i = 1, nphih
    if (dihedh(i)%atm_l .lt. 0) then
      if (dihedh(i)%parm_idx .lt. nimprp) &
        nimprp = dihedh(i)%parm_idx
    end if
  end do

  do i = 1, nphia
    if (diheda(i)%atm_l .lt. 0) then
      if (diheda(i)%parm_idx .lt. nimprp) &
        nimprp = diheda(i)%parm_idx
    end if
  end do

  nimprp = numphi - nimprp + 1

  return

end subroutine set_num_improp_dihed

end module nmr_calls_mod
