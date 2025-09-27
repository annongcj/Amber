! <compile=optimized>
! This file contains non standard compliant code that must be compiled with
! optimizations in order to avoid SIGSEGV, invalid memory references.
! Search on SIGSEGV for details.

#include "copyright.i"
! Original GB code

!*******************************************************************************
!
! Module:  gb_ene_mod
!
! Description: <TBS>
!
!*******************************************************************************

module gb_ene_mod

  use gbl_datatypes_mod
  use gbsa_mod

  implicit none

! Global data definitions:

  double precision, allocatable, save           :: gbl_rbmax(:)
  double precision, allocatable, save           :: gbl_rbmin(:)
  double precision, allocatable, save           :: gbl_rbave(:)
  double precision, allocatable, save           :: gbl_rbfluct(:)
  double precision, allocatable, save           :: r2x(:)
  double precision, allocatable, save           :: reff(:)
  integer,          allocatable, save           :: jj(:)

! Atom-exclusion array for continuous pHMD;
! this is only allocated if iphmd .ne. 0.
  integer,          allocatable, save           :: excludeatoms(:)

! Private data definitions:

  double precision, allocatable, save, private  :: psi(:)
  double precision, allocatable, save, private  :: rjx(:)
  double precision, allocatable, save, private  :: sumdeijda(:)
  double precision, allocatable, save, private  :: vectmp1(:)
  double precision, allocatable, save, private  :: vectmp2(:)
  double precision, allocatable, save, private  :: vectmp3(:)
  double precision, allocatable, save, private  :: vectmp4(:)
  double precision, allocatable, save, private  :: vectmp5(:)

  double precision, allocatable, save, private  :: gb_alpha_arry(:)
  double precision, allocatable, save, private  :: gb_beta_arry(:)
  double precision, allocatable, save, private  :: gb_gamma_arry(:)

! Integer mask specifying if atom is in ion and which ion. 0:not ion;
! >0: res name is ion_target_res_arr(atom_ion_mask(i))
  integer, allocatable, save, private, dimension(:) :: atom_ion_mask
  double precision, save, private, dimension(0:5)  :: atom_ion_coef

  ! target ion residue names array
  character(len=4), save, private, dimension(5) :: ion_target_res_arr = (/ 'CO6 ', 'Na+ ', 'Cl- ', 'K+  ', 'Mg2+' /)

  ! For GPU-SASA private data definitions:
  double precision, allocatable, save, private  :: AgsuminvRn(:)
  double precision, allocatable, save, private  :: sigmaNP(:)
  double precision, allocatable, save, private  :: epsilonNP(:)
  double precision, allocatable, save, private  :: radiusNP(:)
  double precision, allocatable, save, private  :: maxSASANP(:)
  integer,          allocatable, save, private  :: atm_residue(:)
  integer,          allocatable, save, private  :: bond_partner(:)
  integer,          allocatable, save, private  :: total_bond_partner(:)
  ! this radius is the internal radius of each element augmented by 1.4Angstroms

  logical,          allocatable, save, private  :: skipv(:)

  integer,          allocatable, save, private  :: neck_idx(:)

  ! gb_neckcut: 2.8d0 (diameter of water) is "correct" value but
  ! larger values give smaller discontinuities at the cut:

  double precision, parameter    :: gb_neckcut = 6.8d0

  ! Lookup tables for position (atom separation, r) and value of the maximum
  ! of the neck function for given atomic radii ri and rj. Values of neck
  ! maximum are already divided by 4*Pi to save time. Values are given
  ! for each 0.05 angstrom between 1.0 and 2.0 (inclusive), so map to index
  ! with dnint((r-1.0)*20)).  Values were numerically determined in
  ! Mathematica; note FORTRAN column-major array storage, so the data below
  ! may be transposed from how you might expect it.

   double precision, parameter  :: neckMaxPos(0:20,0:20) = reshape((/ &
     2.26685,2.32548,2.38397,2.44235,2.50057,2.55867,2.61663,2.67444, &
     2.73212,2.78965,2.84705,2.9043,2.96141,3.0184,3.07524,3.13196, &
     3.18854,3.24498,3.30132,3.35752,3.4136, &
     2.31191,2.37017,2.4283,2.48632,2.5442,2.60197,2.65961,2.71711, &
     2.77449,2.83175,2.88887,2.94586,3.00273,3.05948,3.1161,3.1726, &
     3.22897,3.28522,3.34136,3.39738,3.45072, &
     2.35759,2.41549,2.47329,2.53097,2.58854,2.646,2.70333,2.76056, &
     2.81766,2.87465,2.93152,2.98827,3.0449,3.10142,3.15782,3.21411, &
     3.27028,3.32634,3.3823,3.43813,3.49387, &
     2.4038,2.46138,2.51885,2.57623,2.63351,2.69067,2.74773,2.80469, &
     2.86152,2.91826,2.97489,3.0314,3.08781,3.1441,3.20031,3.25638, &
     3.31237,3.36825,3.42402,3.4797,3.53527, &
     2.45045,2.50773,2.56492,2.62201,2.679,2.7359,2.7927,2.8494,2.90599, &
     2.9625,3.0189,3.07518,3.13138,3.18748,3.24347,3.29937,3.35515, &
     3.41085,3.46646,3.52196,3.57738, &
     2.4975,2.5545,2.61143,2.66825,2.72499,2.78163,2.83818,2.89464, &
     2.95101,3.00729,3.06346,3.11954,3.17554,3.23143,3.28723,3.34294, &
     3.39856,3.45409,3.50952,3.56488,3.62014, &
     2.54489,2.60164,2.6583,2.71488,2.77134,2.8278,2.88412,2.94034, &
     2.9965,3.05256,3.10853,3.16442,3.22021,3.27592,3.33154,3.38707, &
     3.44253,3.49789,3.55316,3.60836,3.66348, &
     2.59259,2.6491,2.70553,2.76188,2.81815,2.87434,2.93044,2.98646, &
     3.04241,3.09827,3.15404,3.20974,3.26536,3.32089,3.37633,3.4317, &
     3.48699,3.54219,3.59731,3.65237,3.70734, &
     2.64054,2.69684,2.75305,2.80918,2.86523,2.92122,2.97712,3.03295, &
     3.0887,3.14437,3.19996,3.25548,3.31091,3.36627,3.42156,3.47677, &
     3.5319,3.58695,3.64193,3.69684,3.75167, &
     2.68873,2.74482,2.80083,2.85676,2.91262,2.96841,3.02412,3.07976, &
     3.13533,3.19082,3.24623,3.30157,3.35685,3.41205,3.46718,3.52223, &
     3.57721,3.63213,3.68696,3.74174,3.79644, &
     2.73713,2.79302,2.84884,2.90459,2.96027,3.01587,3.0714,3.12686, &
     3.18225,3.23757,3.29282,3.34801,3.40313,3.45815,3.51315,3.56805, &
     3.6229,3.67767,3.73237,3.78701,3.84159, &
     2.78572,2.84143,2.89707,2.95264,3.00813,3.06356,3.11892,3.17422, &
     3.22946,3.28462,3.33971,3.39474,3.44971,3.5046,3.55944,3.61421, &
     3.66891,3.72356,3.77814,3.83264,3.8871, &
     2.83446,2.89,2.94547,3.00088,3.05621,3.11147,3.16669,3.22183, &
     3.27689,3.33191,3.38685,3.44174,3.49656,3.55132,3.60602,3.66066, &
     3.71523,3.76975,3.82421,3.8786,3.93293, &
     2.88335,2.93873,2.99404,3.04929,3.10447,3.15959,3.21464,3.26963, &
     3.32456,3.37943,3.43424,3.48898,3.54366,3.5983,3.65287,3.70737, &
     3.76183,3.81622,3.87056,3.92484,3.97905, &
     2.93234,2.9876,3.04277,3.09786,3.15291,3.20787,3.26278,3.31764, &
     3.37242,3.42716,3.48184,3.53662,3.591,3.64551,3.69995,3.75435, &
     3.80867,3.86295,3.91718,3.97134,4.02545, &
     2.98151,3.0366,3.09163,3.14659,3.20149,3.25632,3.3111,3.36581, &
     3.42047,3.47507,3.52963,3.58411,3.63855,3.69293,3.74725,3.80153, &
     3.85575,3.90991,3.96403,4.01809,4.07211, &
     3.03074,3.08571,3.14061,3.19543,3.25021,3.30491,3.35956,3.41415, &
     3.46869,3.52317,3.57759,3.63196,3.68628,3.74054,3.79476,3.84893, &
     3.90303,3.95709,4.01111,4.06506,4.11897, &
     3.08008,3.13492,3.1897,3.2444,3.29905,3.35363,3.40815,3.46263, &
     3.51704,3.57141,3.62572,3.67998,3.73418,3.78834,3.84244,3.8965, &
     3.95051,4.00447,4.05837,4.11224,4.16605, &
     3.12949,3.18422,3.23888,3.29347,3.348,3.40247,3.45688,3.51124, &
     3.56554,3.6198,3.674,3.72815,3.78225,3.83629,3.8903,3.94425, &
     3.99816,4.05203,4.10583,4.15961,4.21333, &
     3.17899,3.23361,3.28815,3.34264,3.39706,3.45142,3.50571,3.55997, &
     3.61416,3.66831,3.72241,3.77645,3.83046,3.8844,3.93831,3.99216, &
     4.04598,4.09974,4.15347,4.20715,4.26078, &
     3.22855,3.28307,3.33751,3.39188,3.4462,3.50046,3.55466,3.6088, &
     3.6629,3.71694,3.77095,3.82489,3.8788,3.93265,3.98646,4.04022, &
     4.09395,4.14762,4.20126,4.25485,4.3084 &
     /), (/21,21/))

   double precision, parameter  :: neckMaxVal(0:20,0:20) = reshape((/ &
     0.0381511,0.0338587,0.0301776,0.027003,0.0242506,0.0218529, &
     0.0197547,0.0179109,0.0162844,0.0148442,0.0135647,0.0124243, &
     0.0114047,0.0104906,0.00966876,0.008928,0.0082587,0.00765255, &
     0.00710237,0.00660196,0.00614589, &
     0.0396198,0.0351837,0.0313767,0.0280911,0.0252409,0.0227563, &
     0.0205808,0.0186681,0.0169799,0.0154843,0.014155,0.0129696, &
     0.0119094,0.0109584,0.0101031,0.00933189,0.0086348,0.00800326, &
     0.00742986,0.00690814,0.00643255, &
     0.041048,0.0364738,0.0325456,0.0291532,0.0262084,0.0236399, &
     0.0213897,0.0194102,0.0176622,0.0161129,0.0147351,0.0135059, &
     0.0124061,0.0114192,0.0105312,0.00973027,0.00900602,0.00834965, &
     0.0077535,0.00721091,0.00671609, &
     0.0424365,0.0377295,0.0336846,0.0301893,0.0271533,0.0245038, &
     0.0221813,0.0201371,0.018331,0.0167295,0.0153047,0.014033, &
     0.0128946,0.0118727,0.0109529,0.0101229,0.00937212,0.00869147, &
     0.00807306,0.00751003,0.00699641, &
     0.0437861,0.0389516,0.0347944,0.0311998,0.0280758,0.0253479, &
     0.0229555,0.0208487,0.0189864,0.0173343,0.0158637,0.0145507, &
     0.0133748,0.0123188,0.0113679,0.0105096,0.0097329,0.00902853, &
     0.00838835,0.00780533,0.0072733, &
     0.0450979,0.0401406,0.0358753,0.0321851,0.0289761,0.0261726, &
     0.0237125,0.0215451,0.0196282,0.017927,0.0164121,0.0150588, &
     0.0138465,0.0127573,0.0117761,0.0108902,0.0100882,0.00936068, &
     0.00869923,0.00809665,0.00754661, &
     0.0463729,0.0412976,0.0369281,0.0331456,0.0298547,0.026978, &
     0.0244525,0.0222264,0.0202567,0.0185078,0.0169498,0.0155575, &
     0.0143096,0.0131881,0.0121775,0.0112646,0.010438,0.00968781, &
     0.00900559,0.00838388,0.00781622, &
     0.0476123,0.0424233,0.0379534,0.034082,0.0307118,0.0277645, &
     0.0251757,0.0228927,0.0208718,0.0190767,0.0174768,0.0160466, &
     0.0147642,0.0136112,0.0125719,0.0116328,0.0107821,0.0100099, &
     0.00930735,0.00866695,0.00808206, &
     0.0488171,0.0435186,0.038952,0.0349947,0.0315481,0.0285324, &
     0.0258824,0.0235443,0.0214738,0.0196339,0.0179934,0.0165262, &
     0.0152103,0.0140267,0.0129595,0.0119947,0.0111206,0.0103268, &
     0.00960445,0.00894579,0.00834405, &
     0.0499883,0.0445845,0.0399246,0.0358844,0.032364,0.0292822, &
     0.0265729,0.0241815,0.0220629,0.0201794,0.0184994,0.0169964, &
     0.0156479,0.0144345,0.0133401,0.0123504,0.0114534,0.0106386, &
     0.00989687,0.00922037,0.00860216, &
     0.0511272,0.0456219,0.040872,0.0367518,0.0331599,0.0300142, &
     0.0272475,0.0248045,0.0226392,0.0207135,0.0189952,0.0174574, &
     0.0160771,0.0148348,0.0137138,0.0126998,0.0117805,0.0109452, &
     0.0101846,0.00949067,0.00885636, &
     0.0522348,0.0466315,0.0417948,0.0375973,0.0339365,0.030729, &
     0.0279067,0.0254136,0.023203,0.0212363,0.0194809,0.0179092, &
     0.016498,0.0152275,0.0140807,0.013043,0.012102,0.0112466, &
     0.0104676,0.00975668,0.00910664, &
     0.0533123,0.0476145,0.042694,0.0384218,0.0346942,0.0314268, &
     0.0285507,0.026009,0.0237547,0.0217482,0.0199566,0.018352, &
     0.0169108,0.0156128,0.0144408,0.0133801,0.0124179,0.011543, &
     0.010746,0.0100184,0.00935302, &
     0.0543606,0.0485716,0.04357,0.0392257,0.0354335,0.0321082, &
     0.02918,0.0265913,0.0242943,0.0222492,0.0204225,0.0187859, &
     0.0173155,0.0159908,0.0147943,0.0137111,0.0127282,0.0118343, &
     0.0110197,0.0102759,0.00959549, &
     0.0553807,0.0495037,0.0444239,0.0400097,0.0361551,0.0327736, &
     0.0297949,0.0271605,0.0248222,0.0227396,0.0208788,0.0192111, &
     0.0177122,0.0163615,0.0151413,0.0140361,0.013033,0.0121206, &
     0.0112888,0.0105292,0.00983409, &
     0.0563738,0.0504116,0.0452562,0.0407745,0.0368593,0.0334235, &
     0.0303958,0.0277171,0.0253387,0.0232197,0.0213257,0.0196277, &
     0.0181013,0.0167252,0.0154817,0.0143552,0.0133325,0.0124019, &
     0.0115534,0.0107783,0.0100688, &
     0.0573406,0.0512963,0.0460676,0.0415206,0.0375468,0.0340583, &
     0.030983,0.0282614,0.0258441,0.0236896,0.0217634,0.020036, &
     0.0184826,0.017082,0.0158158,0.0146685,0.0136266,0.0126783, &
     0.0118135,0.0110232,0.0102998, &
     0.0582822,0.0521584,0.0468589,0.0422486,0.038218,0.0346784, &
     0.0315571,0.0287938,0.0263386,0.0241497,0.0221922,0.0204362, &
     0.0188566,0.0174319,0.0161437,0.0149761,0.0139154,0.0129499, &
     0.0120691,0.0112641,0.0105269, &
     0.0591994,0.0529987,0.0476307,0.042959,0.0388734,0.0352843, &
     0.0321182,0.0293144,0.0268225,0.0246002,0.0226121,0.0208283, &
     0.0192232,0.0177751,0.0164654,0.015278,0.0141991,0.0132167, &
     0.0123204,0.0115009,0.0107504, &
     0.0600932,0.053818,0.0483836,0.0436525,0.0395136,0.0358764, &
     0.0326669,0.0298237,0.0272961,0.0250413,0.0230236,0.0212126, &
     0.0195826,0.0181118,0.0167811,0.0155744,0.0144778,0.0134789, &
     0.0125673,0.0117338,0.0109702, &
     0.0609642,0.0546169,0.0491183,0.0443295,0.0401388,0.036455, &
     0.0332033,0.030322,0.0277596,0.0254732,0.0234266,0.0215892, &
     0.0199351,0.018442,0.0170909,0.0158654,0.0147514,0.0137365, &
     0.0128101,0.0119627,0.0111863 &
     /), (/21,21/))

contains

!*******************************************************************************
!
! Subroutine:  final_gb_setup
!
! Description: <TBS>
!
!*******************************************************************************

subroutine final_gb_setup(atm_cnt, num_ints, num_reals, my_igb)

  use gbl_constants_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: atm_cnt
  integer, intent(in)           :: my_igb

!RCW Why is there my_igb here instead of igb - this gets passed in as cph_igb
!but I can't figure why one did not just modify the igb in mdin_ctrl_dat.
!This seems to be ripe for mistakes. And it is only used in the allocation why
!not in the gb_ene routine itself which just references igb?

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer                       :: alloc_failed
  integer                       :: i
  integer                       :: j
  integer                       :: atm_1
  integer                       :: atm_2
  integer                       :: num_bonds
  integer                       :: total_bonds

  integer                       :: atomicnumber
  character(len=4)              :: isymbl ! holder for atom type
  
  if (my_igb /= 6) then
    allocate(reff(atm_cnt), &
             psi(atm_cnt), &
             rjx(atm_cnt), &
             r2x(atm_cnt), &
             sumdeijda(atm_cnt), &
             vectmp1(atm_cnt), &
             vectmp2(atm_cnt), &
             vectmp3(atm_cnt), &
             vectmp4(atm_cnt), &
             vectmp5(atm_cnt), &
             jj(atm_cnt), &
             skipv(0:atm_cnt), &
             gb_alpha_arry(atm_cnt), &
             gb_beta_arry(atm_cnt), &
             gb_gamma_arry(atm_cnt), &
  ! we must allocate memories for our added dimensions:
             AgsuminvRn(atm_cnt), &
             sigmaNP(atm_cnt), &
             epsilonNP(atm_cnt), &
             radiusNP(atm_cnt), &
             maxSASANP(atm_cnt),&
             atm_residue(atm_cnt),&
             bond_partner(atm_cnt),&
             total_bond_partner(atm_cnt),&
             stat = alloc_failed)
  else
    allocate( &
             r2x(atm_cnt), &
             vectmp5(atm_cnt), &
             jj(atm_cnt), &
             skipv(0:atm_cnt), &
             stat = alloc_failed)
  end if

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_reals = num_reals + size(r2x) + size(vectmp5)

  if (my_igb /= 6) then
    num_reals = num_reals + size(psi) + &
                            size(rjx) + &
                            size(sumdeijda) + &
                            size(vectmp1) + &
                            size(vectmp2) + &
                            size(vectmp3) + &
                            size(vectmp4) + &
                            size(gb_alpha_arry) + &
                            size(gb_beta_arry) + &
                            size(gb_gamma_arry)
  end if

  num_ints = num_ints + size(jj) + &
                        size(skipv)     ! logical assumed same size as integer.

  if (rbornstat .ne. 0) then

    allocate(gbl_rbmax(atm_cnt), &
             gbl_rbmin(atm_cnt), &
             gbl_rbave(atm_cnt), &
             gbl_rbfluct(atm_cnt), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    num_reals = num_reals + size(gbl_rbmax) + &
                            size(gbl_rbmin) + &
                            size(gbl_rbave) + &
                            size(gbl_rbfluct)


    gbl_rbmax(:) = 0.d0
    gbl_rbmin(:) = 999.d0
    gbl_rbave(:) = 0.d0
    gbl_rbfluct(:) = 0.d0

  end if

  if (my_igb .eq. 7 .or. my_igb .eq. 8) then

    allocate(neck_idx(atm_cnt), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    num_ints = num_ints + size(neck_idx)

    ! Some final error checking before run start for my_igb 7

    do i = 1, atm_cnt

      neck_idx(i) = dnint((atm_gb_radii(i) - 1.d0) * 20.d0)

      if (neck_idx(i) .lt. 0 .or. neck_idx(i) .gt. 20) then

        if (master) then
          write(mdout, '(a,a,i6,a,f7.3,a)') error_hdr, 'Atom ', i, &
            ' has a radius (', atm_gb_radii(i), ') outside the allowed range of'
          write(mdout, '(a,a,a)') extra_line_hdr, &
            '1.0 - 2.0 angstrom for igb=7 and 8. ', &
            'Regenerate prmtop with bondi radii.'
        end if

        call mexit(mdout, 1)

      end if

    end do

  end if

  ! gbion options(gbion=1); initialization
  if (gbion .eq. 1) then
    allocate( atom_ion_mask(atm_cnt), &
              stat = alloc_failed )
    if (alloc_failed .ne. 0) call setup_alloc_error
    num_ints = num_ints + size(atom_ion_mask)
    atom_ion_mask = 0
    atom_ion_coef(0) = 1.0d0
    atom_ion_coef(1:5) = gi_coef

    ! update values in atom_ion_mask() and atom_ion_coef(); structured to handle
    ! ion-specific gi_coef later
    do i=1,nres
      if ( ion_target_res_arr(1) == gbl_labres(i)(1:4) ) then
        do j=gbl_res_atms(i),gbl_res_atms(i+1)-1
          atom_ion_mask(j) = 1
        end do
      else if ( ion_target_res_arr(2) == gbl_labres(i)(1:4) ) then 
        do j=gbl_res_atms(i),gbl_res_atms(i+1)-1
          atom_ion_mask(j) = 2
        end do
      else if ( ion_target_res_arr(3) == gbl_labres(i)(1:4) ) then
        do j=gbl_res_atms(i),gbl_res_atms(i+1)-1
          atom_ion_mask(j) = 3
        end do
      else if ( ion_target_res_arr(4) == gbl_labres(i)(1:4) ) then
        do j=gbl_res_atms(i),gbl_res_atms(i+1)-1
          atom_ion_mask(j) = 4
        end do
      else if ( ion_target_res_arr(5) == gbl_labres(i)(1:4) ) then
        do j=gbl_res_atms(i),gbl_res_atms(i+1)-1
          atom_ion_mask(j) = 5
        end do
      end if
    end do

  end if

  ! gbion options(gbion=2); initialization
  if (gbion .eq. 2) then
    allocate( atom_ion_mask(atm_cnt), &
              stat = alloc_failed )
    if (alloc_failed .ne. 0) call setup_alloc_error
    num_ints = num_ints + size(atom_ion_mask)
    atom_ion_mask = 0
    atom_ion_coef(0) = gi_coef_1
    atom_ion_coef(1:5) = gi_coef_2

    ! update values in atom_ion_mask() and atom_ion_coef(); structured to handle
    ! ion-specific gi_coef later
    do i=1,nres
      if ( ion_target_res_arr(1) == gbl_labres(i)(1:4) ) then
        do j=gbl_res_atms(i),gbl_res_atms(i+1)-1
          atom_ion_mask(j) = 1
        end do
      else if ( ion_target_res_arr(2) == gbl_labres(i)(1:4) ) then 
        do j=gbl_res_atms(i),gbl_res_atms(i+1)-1
          atom_ion_mask(j) = 2
        end do
      else if ( ion_target_res_arr(3) == gbl_labres(i)(1:4) ) then
        do j=gbl_res_atms(i),gbl_res_atms(i+1)-1
          atom_ion_mask(j) = 3
        end do
      else if ( ion_target_res_arr(4) == gbl_labres(i)(1:4) ) then
        do j=gbl_res_atms(i),gbl_res_atms(i+1)-1
          atom_ion_mask(j) = 4
        end do
      else if ( ion_target_res_arr(5) == gbl_labres(i)(1:4) ) then
        do j=gbl_res_atms(i),gbl_res_atms(i+1)-1
          atom_ion_mask(j) = 5
        end do
      end if
    end do

  end if

  ! gbion options(gbion=3); initialization
  if (gbion .eq. 3) then
    allocate( atom_ion_mask(atm_cnt), &
              stat = alloc_failed )
    if (alloc_failed .ne. 0) call setup_alloc_error
    num_ints = num_ints + size(atom_ion_mask)
    atom_ion_mask = 0
    ! atom_ion_coef(0) = gi_coef_1
    ! atom_ion_coef(1:5) = gi_coef_2
    ! no longer need atom_ion_coef() for gbion=3

    ! update values in atom_ion_mask() and atom_ion_coef(); structured to handle
    ! ion-specific gi_coef later
    do i=1,nres
      if ( ion_target_res_arr(1) == gbl_labres(i)(1:4) ) then
        do j=gbl_res_atms(i),gbl_res_atms(i+1)-1
          atom_ion_mask(j) = 1
        end do
      else if ( ion_target_res_arr(2) == gbl_labres(i)(1:4) ) then 
        do j=gbl_res_atms(i),gbl_res_atms(i+1)-1
          atom_ion_mask(j) = 2
        end do
      else if ( ion_target_res_arr(3) == gbl_labres(i)(1:4) ) then
        do j=gbl_res_atms(i),gbl_res_atms(i+1)-1
          atom_ion_mask(j) = 3
        end do
      else if ( ion_target_res_arr(4) == gbl_labres(i)(1:4) ) then
        do j=gbl_res_atms(i),gbl_res_atms(i+1)-1
          atom_ion_mask(j) = 4
        end do
      else if ( ion_target_res_arr(5) == gbl_labres(i)(1:4) ) then
        do j=gbl_res_atms(i),gbl_res_atms(i+1)-1
          atom_ion_mask(j) = 5
        end do
      end if
    end do

  end if



  ! Since the development of igb = 8, the gb_alpha, gb_beta, and gb_gamma
  ! parameters have to be stored in arrays (one for each atom), since that
  ! model introduced atom-dependent values for those parameters. We fill
  ! those arrays here. If we're using igb .eq. 1, 2, 5, or 7, then we'll fill
  ! every value with gb_alpha, gb_beta, and gb_gamma set in init_mdin_ctrl_dat

  if (my_igb .eq. 1 .or. my_igb .eq. 2 .or. my_igb .eq. 5 .or. my_igb .eq. 7) then

    gb_alpha_arry(:) = gb_alpha
    gb_beta_arry(:) = gb_beta
    gb_gamma_arry(:) = gb_gamma

  else if (my_igb .eq. 8) then

    do i = 1, atm_cnt

      if (loaded_atm_atomicnumber) then
        atomicnumber = atm_atomicnumber(i)
      else
        call get_atomic_number(atm_igraph(i), atm_mass(i), atomicnumber)
      end if

      call isnucat(nucat,i,nres,60,gbl_res_atms(1:nres),gbl_labres(1:nres))
      ! update GBNeck2nu
      if (nucat == 1) then
          !if atom belong to nucleic part, use nuc pars
          if (atomicnumber .eq. 1) then
            gb_alpha_arry(i) = gb_alpha_hnu
            gb_beta_arry(i) = gb_beta_hnu
            gb_gamma_arry(i) = gb_gamma_hnu
          else if (atomicnumber .eq. 6) then
            gb_alpha_arry(i) = gb_alpha_cnu
            gb_beta_arry(i) = gb_beta_cnu
            gb_gamma_arry(i) = gb_gamma_cnu
          else if (atomicnumber .eq. 7) then
            gb_alpha_arry(i) = gb_alpha_nnu
            gb_beta_arry(i) = gb_beta_nnu
            gb_gamma_arry(i) = gb_gamma_nnu
          else if (atomicnumber .eq. 8) then
            gb_alpha_arry(i) = gb_alpha_osnu
            gb_beta_arry(i) = gb_beta_osnu
            gb_gamma_arry(i) = gb_gamma_osnu
          else if (atomicnumber .eq. 16) then
            gb_alpha_arry(i) = gb_alpha_osnu
            gb_beta_arry(i) = gb_beta_osnu
            gb_gamma_arry(i) = gb_gamma_osnu
          else if (atomicnumber .eq. 15) then
            gb_alpha_arry(i) = gb_alpha_pnu
            gb_beta_arry(i) = gb_beta_pnu
            gb_gamma_arry(i) = gb_gamma_pnu
          else
            ! Use GB^OBC (II) (igb = 5) set for other atoms
            gb_alpha_arry(i) = 1.0d0
            gb_beta_arry(i) = 0.8d0
            gb_gamma_arry(i) = 4.851d0
          end if
      else
          !if not nucleic part, use protein pars
          if (atomicnumber .eq. 1) then
            gb_alpha_arry(i) = gb_alpha_h
            gb_beta_arry(i) = gb_beta_h
            gb_gamma_arry(i) = gb_gamma_h
          else if (atomicnumber .eq. 6) then
            gb_alpha_arry(i) = gb_alpha_c
            gb_beta_arry(i) = gb_beta_c
            gb_gamma_arry(i) = gb_gamma_c
          else if (atomicnumber .eq. 7) then
            gb_alpha_arry(i) = gb_alpha_n
            gb_beta_arry(i) = gb_beta_n
            gb_gamma_arry(i) = gb_gamma_n
          else if (atomicnumber .eq. 8) then
            gb_alpha_arry(i) = gb_alpha_os
            gb_beta_arry(i) = gb_beta_os
            gb_gamma_arry(i) = gb_gamma_os
          else if (atomicnumber .eq. 16) then
            gb_alpha_arry(i) = gb_alpha_os
            gb_beta_arry(i) = gb_beta_os
            gb_gamma_arry(i) = gb_gamma_os
          else if (atomicnumber .eq. 15) then
            gb_alpha_arry(i) = gb_alpha_p
            gb_beta_arry(i) = gb_beta_p
            gb_gamma_arry(i) = gb_gamma_p
          else
            ! Use GB^OBC (II) (igb = 5) set for other atoms
            gb_alpha_arry(i) = 1.0d0
            gb_beta_arry(i) = 0.8d0
            gb_gamma_arry(i) = 4.851d0
          end if
        end if !testing if atom belongs to nuc

    end do

  end if

! Set up the GPU-SASA data structures if we're doing a GB/SA calculation
if (gbsa .eq. 3 ) then
! Built atm_residue(atm_cnt) array
  do i=1,nres
     do j=gbl_res_atms(i),gbl_res_atms(i+1)-1
        atm_residue(j) = i
     end do
  end do

    bond_partner(:) = 0

    ! For each atom index we find in the BONDS_WITHOUT_HYDROGEN array, add
    ! 1 to that index of bond_partner. That way, we can tell how many atoms
    ! are bonded to each atom

    do i = 0, nbona - 1
      atm_1 = gbl_bond(bonda_idx + i)%atm_i
      atm_2 = gbl_bond(bonda_idx + i)%atm_j
      bond_partner(atm_1) = bond_partner(atm_1) + 1
      bond_partner(atm_2) = bond_partner(atm_2) + 1
    end do

    total_bond_partner = bond_partner

    do i = 1, nbonh
      atm_1 = gbl_bond(i)%atm_i
      atm_2 = gbl_bond(i)%atm_j
      total_bond_partner(atm_1) = total_bond_partner(atm_1) + 1
      total_bond_partner(atm_2) = total_bond_partner(atm_2) + 1
    end do

! Look up SASA type by "bonded partner" bond_partner(i),
! "total bonded partner" total_bond_partner(i),  "residue name"
! gbl_labres(atm_residue(i)) and  "atom name" atm_igraph(i)
  do i = 1, atm_cnt
     if ( atm_atomicnumber(i) == 6 ) then
        radiusNP(i) = 3.10d0 ! 1.7d0 + 1.4d0
        if ( total_bond_partner(i) == 4 ) then
           if ( bond_partner(i) == 1) then
              if ( gbl_labres(atm_residue(i))(1:3) == 'MET') then
                 ! SASA type '1CS'
                 sigmaNP(i) = 1.20831904d0
                 epsilonNP(i) = 14.4050341d0
                 maxSASANP(i) = 101.172789338d0
              else if ( gbl_labres(atm_residue(i))(1:3) == 'NME') then
                 ! SASA type '1CN'
                 sigmaNP(i) = 0.317053639d0
                 epsilonNP(i) = 7.14828051d0
                 maxSASANP(i) = 93.3032786658d0
              else
                 ! SASA type '1CC'
                 sigmaNP(i) = 4.37011638d0
                 epsilonNP(i) = 19.5924803d0
                 maxSASANP(i) = 89.5418522949d0
              end if
           else if ( bond_partner(i) == 2) then
              if ( gbl_labres(atm_residue(i))(1:3) == 'SER') then
                 ! SASA type '2CCO'
                 sigmaNP(i) = 1.56800573d0
                 epsilonNP(i) = 4.73805582d0
                 maxSASANP(i) = 62.3149039685d0
              else if ( (gbl_labres(atm_residue(i))(1:3) == 'MET' .and. &
                         atm_igraph(i)(1:2) == 'CG') .or. &
                         gbl_labres(atm_residue(i))(1:2) == 'CY' ) then
                 ! SASA type '2CCS'
                 sigmaNP(i) = 3.85601041d0
                 epsilonNP(i) = 13.7926170d0
                 maxSASANP(i) = 72.0652500364d0
              else if ( (gbl_labres(atm_residue(i))(1:3) == 'ARG' .and. &
                         atm_igraph(i)(1:2) == 'CD') .or. &
                         (gbl_labres(atm_residue(i))(1:3) == 'LYS' .and. &
                         atm_igraph(i)(1:2) == 'CE') .or. &
                         (gbl_labres(atm_residue(i))(1:3) == 'PRO' .and. &
                         atm_igraph(i)(1:2) == 'CD') .or. &
                         gbl_labres(atm_residue(i))(1:3) == 'GLY' ) then
                 ! SASA type '2CCN'
                 sigmaNP(i) = 5.49283775d0
                 epsilonNP(i) = 19.2142035d0
                 maxSASANP(i) = 78.9232674787d0
              else
                 ! SASA type '2CCC'
                 sigmaNP(i) = 7.24965924d0
                 epsilonNP(i) = 18.7931978d0
                 maxSASANP(i) = 67.848137034d0
              end if
           else !   bond_partner(i) == 3
              if ( gbl_labres(atm_residue(i))(1:3) == 'THR' .and. &
                   atm_igraph(i)(1:2) == 'CB') then
                 ! SASA type '4CCO'
                 sigmaNP(i) = 1.84288147d0
                 epsilonNP(i) = 1.60001261d0
                 maxSASANP(i) = 33.1334096093d0
              else if ( (gbl_labres(atm_residue(i))(1:3) == 'ILE' .and. &
                         atm_igraph(i)(1:2) == 'CB') .or. &
                        (gbl_labres(atm_residue(i))(1:3) == 'LEU' .and. &
                         atm_igraph(i)(1:2) == 'CG') .or. &
                        (gbl_labres(atm_residue(i))(1:3) == 'VAL' .and. &
                         atm_igraph(i)(1:2) == 'CB') ) then
                 ! SASA type '4CCC'
                 sigmaNP(i) = 2.46334544d0
                 epsilonNP(i) = 1.58675633d0
                 maxSASANP(i) = 29.2869985239d0
              else
                 ! SASA type '4CCN'
                 sigmaNP(i) = 0.100000d0
                 epsilonNP(i) = 0.328277427d0
                 maxSASANP(i) = 22.5608806495d0
              end if
           end if
        else ! Carbon  total_bond_partner(i) == 3
           if (  bond_partner(i) == 2 ) then
              if ( (gbl_labres(atm_residue(i))(1:3) == 'TRP' .and. &
                    atm_igraph(i)(1:3) == 'CD1') .or. &
                   (gbl_labres(atm_residue(i))(1:2) == 'HI' .and. &
                    atm_igraph(i)(1:3) == 'CD2') ) then
                 ! SASA type '3CCN'
                 sigmaNP(i) = 1.13740114d0
                 epsilonNP(i) = 4.91448211d0
                 maxSASANP(i) = 70.6469194565d0
              else if ( (gbl_labres(atm_residue(i))(1:2) == 'HI' .and. &
                         atm_igraph(i)(1:3) == 'CE1') ) then
                 ! SASA type '3CNN'
                 sigmaNP(i) = 4.79302099d0
                 epsilonNP(i) = 17.7860582d0
                 maxSASANP(i) = 85.8858602953d0
              else
                 ! SASA type '3CC'
                 sigmaNP(i) = 5.52380670d0
                 epsilonNP(i) = 13.1799067d0
                 maxSASANP(i) = 72.2402965204d0
              end if
            else ! bond_partner(i) == 3
               if ( (gbl_labres(atm_residue(i))(1:3) == 'TYR' .and. &
                     atm_igraph(i)(1:2) == 'CG') .or. &
                    (gbl_labres(atm_residue(i))(1:3) == 'TRP' .and. &
                     atm_igraph(i)(1:2) == 'CG') .or. &
                    (gbl_labres(atm_residue(i))(1:3) == 'TRP' .and. &
                     atm_igraph(i)(1:3) == 'CD2') .or. &
                    (gbl_labres(atm_residue(i))(1:3) == 'PHE' .and. &
                     atm_igraph(i)(1:2) == 'CG') ) then
                 ! SASA type '3CCC'
                 sigmaNP(i) = 7.92563682d0
                 epsilonNP(i) = 0.700404609d0
                 maxSASANP(i) = 14.8528474421d0
               else if ( (gbl_labres(atm_residue(i))(1:3) == 'TYR' .and. &
                          atm_igraph(i)(1:2) == 'CZ') ) then
                 ! SASA type '3CCO'
                 sigmaNP(i) = 2.85247115d0
                 epsilonNP(i) = 0.690432536d0
                 maxSASANP(i) = 20.2921312847d0
               else if ( (gbl_labres(atm_residue(i))(1:3) == 'TRP' .and. &
                     atm_igraph(i)(1:3) == 'CE2') ) then
                 ! SASA type '5CCN1'
                 sigmaNP(i) = 3.53275857d0
                 epsilonNP(i) = 0.371096858d0
                 maxSASANP(i) = 16.5093987597d0
               else if ( (gbl_labres(atm_residue(i))(1:2) == 'HI' .and. &
                          atm_igraph(i)(1:2) == 'CG') ) then
                 ! SASA type '5CCN2'
                 sigmaNP(i) = 0.902827736d0
                 epsilonNP(i) = 0.0212409206d0
                 maxSASANP(i) = 6.90712731848d0
               else if ( (gbl_labres(atm_residue(i))(1:3) == 'ARG' .and. &
                          atm_igraph(i)(1:2) == 'CZ') ) then
                 ! SASA type '5CNN'
                 sigmaNP(i) = 6.51644207d0
                 epsilonNP(i) = 3.68184016d0
                 maxSASANP(i) = 32.0850765017d0
               else if ( (gbl_labres(atm_residue(i))(1:3) == 'ASP' .and. &
                          atm_igraph(i)(1:2) == 'CG') .or. &
                         (gbl_labres(atm_residue(i))(1:3) == 'GLU' .and. &
                          atm_igraph(i)(1:2) == 'CD') .or. &
                         (atm_residue(i) == nres .and. &
                          gbl_labres(atm_residue(i))(1:3) /= 'NME' .and. &
                          gbl_labres(atm_residue(i))(1:3) /= 'NHE') ) then
                 ! SASA type '5COO'
                 sigmaNP(i) = 9.77659495d0
                 epsilonNP(i) = 1.43809938d0
                 maxSASANP(i) = 16.6436516422d0
               else
                 ! SASA type '5CNO'
                 sigmaNP(i) = 5.99708166d0
                 epsilonNP(i) = 0.739936182d0
                 maxSASANP(i) = 16.1244401843d0
               end if
            end if
        end if
     else if ( atm_atomicnumber(i) == 7 ) then
        radiusNP(i) = 2.950d0 ! 1.55d0 + 1.4d0
        if ( total_bond_partner(i) == 4 ) then
                 ! SASA type '1NC2'
                 sigmaNP(i) = 4.29095538d0
                 epsilonNP(i) = 34.2745747d0
                 maxSASANP(i) = 95.1108847805d0
        else !  total_bond_partner(i) == 3
           if ( bond_partner(i) == 1 ) then
                 ! SASA type '1NC1'
                 sigmaNP(i) = 3.00848540d0
                 epsilonNP(i) = 23.5119766d0
                 maxSASANP(i) = 94.0970695867d0
           else if ( bond_partner(i) == 2 ) then
              if ( (gbl_labres(atm_residue(i))(1:3) == 'TRP' .and. &
                    atm_igraph(i)(1:3) == 'NE1') .or. &
                   (gbl_labres(atm_residue(i))(1:2) == 'HI' .and. &
                    atm_igraph(i)(1:3) == 'ND1') .or. &
                   (gbl_labres(atm_residue(i))(1:2) == 'HI' .and. &
                    atm_igraph(i)(1:3) == 'NE2') ) then
                 ! SASA type '3NCC'
                 sigmaNP(i) = 5.58999754d0
                 epsilonNP(i) = 16.2639366d0
                 maxSASANP(i) = 64.7488801386d0
              else
                 ! SASA type '2NCC'
                 sigmaNP(i) = 3.29603051d0
                 epsilonNP(i) = 0.919202164d0
                 maxSASANP(i) = 22.8610751485d0
              end if
           else ! bond_partner(i) == 3 for PRO amide
                 ! SASA type '4NCC'
                 sigmaNP(i) = 1.0d0
                 epsilonNP(i) = 0.0d0
                 maxSASANP(i) = 0.180783832809d0
           end if
        end if
     else if ( atm_atomicnumber(i) == 8 ) then
        radiusNP(i) = 2.90d0 ! 1.5d0 + 1.4d0
        if ( total_bond_partner(i) == 2 ) then
           if ( gbl_labres(atm_residue(i))(1:3) == 'TYR' ) then
                 ! SASA type '1OC3'
                 sigmaNP(i) = 2.82723035d0
                 epsilonNP(i) = 11.2361169d0
                 maxSASANP(i) = 80.1047057149d0
           else ! SER, THR hydroxyl
                 ! SASA type '1OC2', same parameter as '1OC3'
                 sigmaNP(i) = 2.82723035d0
                 epsilonNP(i) = 11.2361169d0
                 maxSASANP(i) = 69.8105657286d0
           end if
        else
                 ! SASA type '1OC1'
                 sigmaNP(i) = 6.76485820d0
                 epsilonNP(i) = 12.6706340d0
                 maxSASANP(i) = 58.3692979586d0
        end if
     else if ( atm_atomicnumber(i) == 16 ) then
        radiusNP(i) = 3.20d0 ! 1.8d0 + 1.4d0
        if ( bond_partner(i) == 1 ) then
                 ! SASA type '1SC'
                 sigmaNP(i) = 2.52036225d0
                 epsilonNP(i) = 16.7889854d0
                 maxSASANP(i) = 105.113824567d0
        else ! bond_partner(i) == 2
                 ! SASA type '2SCC'
                 sigmaNP(i) = 1.13372522d0
                 epsilonNP(i) = 5.82867046d0
                 maxSASANP(i) = 75.8598197695d0
        end if
     else ! All the hydrogen atoms
                 radiusNP(i) = 1.4d0
                 sigmaNP(i) = 1.0d0
                 epsilonNP(i) = 0.0d0
                 maxSASANP(i) = 0.0d0
     end if
  end do

endif !(gbsa == 3)

  if (gbsa .eq. 1) then
    call gbsa_setup(atm_cnt, num_ints, num_reals)
  end if ! (gbsa .eq. 1)

  ! We needed atm_isymbl, etc. for igb=8 and gbsa initialization, but now
  ! we have no more use for it. So deallocate it.

  num_ints = num_ints - size(atm_isymbl) - &
               size(atm_atomicnumber)
  if (allocated(atm_isymbl)) deallocate(atm_isymbl)
  if (allocated(atm_isymbl)) deallocate(atm_atomicnumber)



#ifdef CUDA
  call gpu_create_outputbuffers()
  if (my_igb .ne. 6) then !DJM
    call gpu_upload_rborn(atm_gb_radii)
    call gpu_upload_fs(atm_gb_fs)
  end if
  call gpu_final_gb_setup(my_igb, alpb, saltcon, rgbmax, gb_neckcut, gb_fs_max,&
                    atm_gb_radii, atm_gb_fs, neckMaxVal, neckMaxPos)

  if (gbion .eq. 2 .or. gbion .eq. 3) then ! gpu setup for gbion=2 or 3
    call gpu_gb_gbion_setup(gbion, gi_coef_1, gi_coef_2, intdiel_ion_1, intdiel_ion_2,&
                      gb_neckscale_ion_1, gb_neckscale_ion_2, atom_ion_mask, &
                      gi_coef_1_p, gi_coef_1_n, gi_coef_2_pp, gi_coef_2_pn, gi_coef_2_nn, &
                      gb_neckscale_ion_1_p, gb_neckscale_ion_1_n, &
                      gb_neckscale_ion_2_pp, gb_neckscale_ion_2_pn, gb_neckscale_ion_2_nn, &
                      intdiel_ion_1_p, intdiel_ion_1_n, & 
                      intdiel_ion_2_pp, intdiel_ion_2_pn, intdiel_ion_2_nn)
  end if

  if (my_igb .eq. 8) then
    call gpu_gb_igb8_setup(gb_alpha_arry, gb_beta_arry, gb_gamma_arry)
  end if
  call gpu_build_threadblock_work_list(atm_numex, gbl_natex)
  if (gbsa .eq. 1) then
    call gpu_gbsa_setup()
  end if ! (gbsa .eq. 1)
  !gpu_gbsa3 _setup()
  if (gbsa .eq. 3) then
    call gpu_gbsa3_setup(sigmaNP, epsilonNP, radiusNP, maxSASANP)
  end if
#endif

  return

end subroutine final_gb_setup

!*******************************************************************************
!
! Subroutine:  gb_cleanup
!
! Description: Deallocate GB arrays etc.
!
!*******************************************************************************
subroutine gb_cleanup(num_ints,num_reals,my_igb)

  use mdin_ctrl_dat_mod, only : rbornstat, gbion
  use pmemd_lib_mod, only : setup_dealloc_error

  implicit none

  ! num_ints and num_reals are used to return allocation counts. Don't zero.
  integer, intent(in out)       :: num_ints, num_reals
  integer, intent(in)           :: my_igb

  integer cleanup_alloc_failed

  num_reals = num_reals - size(r2x) - size(vectmp5)

  if (my_igb /= 6) then
    num_reals = num_reals - size(psi) - &
                            size(rjx) - &
                            size(sumdeijda) - &
                            size(vectmp1) - &
                            size(vectmp2) - &
                            size(vectmp3) - &
                            size(vectmp4) - &
                            size(gb_alpha_arry) - &
                            size(gb_beta_arry) - &
                            size(gb_gamma_arry)
  end if

  num_ints = num_ints - size(jj) - &
                        size(skipv)     ! logical assumed same size as integer.

  if (my_igb /= 6) then
    deallocate(reff, &
             psi, &
             rjx, &
             r2x, &
             sumdeijda, &
             vectmp1, &
             vectmp2, &
             vectmp3, &
             vectmp4, &
             vectmp5, &
             jj, &
             skipv, &
             gb_alpha_arry, &
             gb_beta_arry, &
             gb_gamma_arry, &
             stat = cleanup_alloc_failed)
  else
    deallocate( &
             r2x, &
             vectmp5, &
             jj, &
             skipv, &
             stat = cleanup_alloc_failed)
  end if

  if (cleanup_alloc_failed .ne. 0) call setup_dealloc_error

  if (rbornstat .ne. 0) then
    num_reals = num_reals - size(gbl_rbmax) - &
                            size(gbl_rbmin) - &
                            size(gbl_rbave) - &
                            size(gbl_rbfluct)

    deallocate(gbl_rbmax, &
             gbl_rbmin, &
             gbl_rbave, &
             gbl_rbfluct, &
             stat = cleanup_alloc_failed)

    if (cleanup_alloc_failed .ne. 0) call setup_dealloc_error

  end if

  if (my_igb .eq. 7 .or. my_igb .eq. 8) then

    num_ints = num_ints - size(neck_idx)

    deallocate(neck_idx, &
             stat = cleanup_alloc_failed)

    if (cleanup_alloc_failed .ne. 0) call setup_dealloc_error

  end if

  ! gbion option
  if (gbion .eq. 1 .or. gbion .eq. 2 .or. gbion .eq. 3) then
    num_ints = num_ints - size(atom_ion_mask)
    deallocate( atom_ion_mask, & 
              stat = cleanup_alloc_failed )
    if (cleanup_alloc_failed .ne. 0) call setup_dealloc_error
  end if


  return

end subroutine gb_cleanup

!*******************************************************************************
!
! Subroutine:  gb_ene
!
! Description: Calculate forces, energies based on Generalized Born.
!
!   Compute nonbonded interactions with a generalized Born model,
!   getting the "effective" Born radii via the approximate pairwise method
!   Use Eqs 9-11 of Hawkins, Cramer, Truhlar, J. Phys. Chem. 100:19824
!   (1996).  Aside from the scaling of the radii, this is the same
!   approach developed in Schaefer and Froemmel, JMB 216:1045 (1990).
!
!   The input coordinates are in the "x" array, and the forces in "f"
!   get updated; energy components are returned in "egb", "eelt" and
!   "evdw".
!
!   Input parameters for the generalized Born model are "rborn(i)", the
!   intrinsic dielectric radius of atom "i", and "fs(i)", which is
!   set (in init_prmtop_dat()) to (rborn(i) - offset)*si.
!
!   Input parameters for the "gas-phase" electrostatic energies are
!   the charges, in the "charge()" array.
!
!   Input parameters for the van der Waals terms are "cn1()" and "cn2()",
!   containing LJ 12-6 parameters, and "asol" and "bsol" containing
!   LJ 12-10 parameters.  (The latter are not used in 1994 and later
!   forcefields.)  The "iac" and "ico" arrays are used to point into
!   these matrices of coefficients.
!
!   The "numex" and "natex" arrays are used to find "excluded" pairs of
!   atoms, for which gas-phase electrostatics and LJ terms are skipped;
!   note that GB terms are computed for all pairs of atoms.
!
!   The code also supports a multiple-time-step facility in which:
!
!   Pairs closer than sqrt(cut_inner) are evaluated every nrespai steps, pairs
!   between sqrt(cut_inner) and sqrt(cut) are evaluated every nrespa steps,
!   and pairs beyond sqrt(cut) are ignored
!
!   The forces arising from the derivatives of the GB terms with respect
!   to the effective Born radii are evaluated every nrespa steps.
!
!   The surface-area dependent term is evaluated every nrespa steps.
!
!   The effective radii are only updated every nrespai steps
!
!   (Be careful with the above: what seems to work is dt=0.001,
!    nrespai=2, nrespa=4; anything beyond this seems dangerous.)
!
!   Written 1999-2000, primarily by D.A. Case, with help from C. Brooks,
!   T. Simonson, R. Sinkovits  and V. Tsui.  The LCPO implementation
!   was written by V. Tsui.
!
!   Vectorization and optimization 1999-2000, primarily by C. P. Sosa,
!   T. Hewitt, and D. A. Case.  Work presented at CUG Fall of 2000.
!
!   NOTE - in the old sander code, the Generalized Born energy was calc'd
!          and returned as epol; here we rename this egb and will pass it
!          all the way out to the run* routines with this name.
!
!*******************************************************************************
subroutine gb_ene(crd, frc, rborn, fs, charge, iac, ico, numex, &
                  natex, atm_cnt, natbel, egb, eelt, evdw, esurf, &
                  irespa, skip_radii_)
!#endif

  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use timers_mod
  use ti_mod

  implicit none

! Formal arguments:

  double precision      :: crd(*)
  double precision      :: frc(*)
  double precision      :: rborn(*)
  double precision      :: fs(*)
  double precision      :: charge(*)
  integer               :: iac(*)
  integer               :: ico(*)
  integer               :: numex(*)
  integer               :: natex(*)
  integer               :: atm_cnt
  integer               :: natbel
  double precision      :: egb, eelt, evdw, esurf, totsasa
  integer, intent(in)   :: irespa
  logical, optional, intent(in) :: skip_radii_

  logical :: skip_radii ! for converting the optional argument

! Local variables:

  double precision      :: cut2, cut_inner2             !
  double precision      :: extdiel_inv                  !
  double precision      :: intdiel_inv                  !
  double precision      :: intdiel_inv_ion_1            ! for ion-solute
  double precision      :: intdiel_inv_ion_2            ! for ion-ion
  double precision      :: intdiel_inv_ion_1_p          ! cation-solute
  double precision      :: intdiel_inv_ion_1_n          ! anion-solute
  double precision      :: intdiel_inv_ion_2_pp         ! cation-cation
  double precision      :: intdiel_inv_ion_2_pn         ! cation-anion
  double precision      :: intdiel_inv_ion_2_nn         ! anion-anion
  double precision      :: ri, rj                       !
  double precision      :: ri1i
  double precision      :: xij, yij, zij                !
  double precision      :: dij1i, dij2i, dij3i          !
  double precision      :: r2                           !
  double precision      :: dij                          !
  double precision      :: sj, sj2                      !
  double precision      :: frespa                       !
  double precision      :: qi, qiqj                     !
  double precision      :: dumx, dumy, dumz             !
  double precision      :: fgbi                         !
  double precision      :: rinv, r2inv, r6inv, r10inv   !
  double precision      :: fgbk                         !
  double precision      :: expmkf                       !
  double precision      :: dl                           !
  double precision      :: de                           !
  double precision      :: e                            !
  double precision      :: temp1                        !
  double precision      :: temp4, temp5, temp6, temp7   !
  double precision      :: eel                          !
  double precision      :: f6, f12, f10                 !
  double precision      :: dedx, dedy, dedz             !
  double precision      :: qi2h, qid2h                  !
  double precision      :: datmp                        !
  double precision      :: thi, thi2                    !
  double precision      :: f_x, f_y, f_z                !
  double precision      :: f_xi, f_yi, f_zi             !
  double precision      :: xi, yi, zi                   !
  double precision      :: dumbo                        !
  double precision      :: tmpsd                        !

  !Variables added for gbsa == 3
  double precision      :: paraA, paraB
  double precision      :: Cij,dist,CijinvS,invS
  double precision      :: tempsum,tempsum1
  double precision      :: tempsum2,tempsum3
  double precision      :: tempsum4,tempsum5
  double precision      :: desurfi,desurfj
  integer               :: orderm, ordern

  ! Variables needed for smooth integration cutoff in Reff:

  double precision      :: rgbmax1i                     !
  double precision      :: rgbmax2i                     !
  double precision      :: rgbmaxpsmax2                 !

  ! Scratch variables used for calculating neck correction:

  double precision      ::  mdist
  double precision      ::  mdist2
  double precision      ::  mdist3
  double precision      ::  mdist5
  double precision      ::  mdist6

  ! Stuff for alpb:

  double precision              :: alpb_beta
  double precision              :: one_arad_beta
  double precision              :: gb_kappa_inv
  ! Alpha prefactor for alpb_alpha:
  double precision, parameter   :: alpb_alpha = 0.571412d0

  integer               :: icount
  integer               :: neibr_cnt
  integer               :: i, j, k
  integer               :: kk1
  integer               :: max_i
  integer               :: iaci
  integer               :: iexcl, jexcl
  integer               :: jexcl_last
  integer               :: jjv
  integer               :: ic
  integer               :: j3
  logical               :: onstep

  ! atom i, j is cation or anion; parameter index 0-8 in 3*3 matrix
  ! parameter matrix
  ! | 1.0(0)   1_n(1)   1_p(2)  | 
  ! | 1_n(3)   2_nn(4)  2_pn(5) |
  ! | 1_p(6)   2_pn(7)  2_pp(8) |
  integer               :: ion_mask_pn_i, ion_mask_pn_j, ion_mask_pn_index

  ! FGB taylor coefficients follow
  ! from A to H :
  ! 1/3 , 2/5 , 3/7 , 4/9 , 5/11
  ! 4/3 , 12/5 , 24/7 , 40/9 , 60/11

  double precision, parameter  :: te = 4.d0 / 3.d0
  double precision, parameter  :: tf = 12.d0 / 5.d0
  double precision, parameter  :: tg = 24.d0 / 7.d0
  double precision, parameter  :: th = 40.d0 / 9.d0
  double precision, parameter  :: thh = 60.d0 / 11.d0

  ! Determine if we want to skip the calculation of the effective radii. There
  ! is an optional argument to force skipping it

  skip_radii = (irespa .gt. 1 .and. mod(irespa, nrespai) .ne. 0)

  ! For gas phase calculations gbradii are not needed so set skip to true.
  if (igb == 6) skip_radii = .true.

  !RCW - What is this code doing? What is skip_radii_?
  if (present(skip_radii_)) &
    skip_radii = skip_radii_ .or. skip_radii

  egb = 0.d0
  eelt = 0.d0
  evdw = 0.d0
  esurf = 0.d0
  totsasa = 0.d0
  max_i = atm_cnt

  ! GPU-SASA: scaling factor multiplied to atomic SASA
  if ( gbsa .eq. 3 ) then
     AgsuminvRn(1:atm_cnt) = 0.d0
#ifdef MPI
     if (mytaskid == 0) then
        totsasa = 361.108307897d0 ! initial maxSASA for molecular SASA transformation
     end if
     do i = mytaskid + 1, max_i, numtasks
#else
     totsasa = 361.108307897d0 ! initial maxSASA for molecular SASA transformation
     do i = 1, max_i
#endif
        AgsuminvRn(i) = 0.681431329392d0 * maxSASANP(i) !MaxSasa * slope
     end do
  end if

  if (mod(irespa, nrespai) .ne. 0) return

  cut2 = gb_cutoff * gb_cutoff
  ! GPU-SASA: scaling factor multiplied to atomic SASA
  cut_inner2 = cut_inner * cut_inner
  onstep = mod(irespa, nrespa) .eq. 0

  if (alpb .eq. 0) then
    ! Standard Still's GB
    extdiel_inv = 1.d0 / extdiel
    intdiel_inv = 1.d0 / intdiel
    if (gbion .eq. 2) then  ! intdiel_ion_1 and _2 only for alpb=1
      intdiel_inv_ion_1 = 1.d0 / intdiel_ion_1
      intdiel_inv_ion_2 = 1.d0 / intdiel_ion_2
    else if (gbion .eq. 3) then  ! separate cation and anion
      intdiel_inv_ion_1_p = 1.d0 / intdiel_ion_1_p
      intdiel_inv_ion_1_n = 1.d0 / intdiel_ion_1_n
      intdiel_inv_ion_2_pp = 1.d0 / intdiel_ion_2_pp
      intdiel_inv_ion_2_pn = 1.d0 / intdiel_ion_2_pn
      intdiel_inv_ion_2_nn = 1.d0 / intdiel_ion_2_nn
    end if
  else
    ! Sigalov Onufriev ALPB (epsilon-dependent GB):
    alpb_beta = alpb_alpha * (intdiel / extdiel)
    extdiel_inv = 1.d0 / (extdiel * (1.d0 + alpb_beta))
    intdiel_inv = 1.d0 / (intdiel * (1.d0 + alpb_beta))
    one_arad_beta = alpb_beta / arad
    if (gb_kappa .ne. 0.d0) gb_kappa_inv = 1.d0 / gb_kappa
  end if

  if (natbel .gt. 0) max_i = natbel

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

  rgbmax1i = 1.d0 / rgbmax
  rgbmax2i = rgbmax1i * rgbmax1i
  rgbmaxpsmax2 = (rgbmax + gb_fs_max)**2

  !---------------------------------------------------------------------------
  ! Step 1: loop over pairs of atoms to compute the effective Born radii.
  !---------------------------------------------------------------------------

  if (.not. skip_radii) call calc_born_radii(atm_cnt, crd, fs, rborn)

  !--------------------------------------------------------------------------
  !
  ! Step 2: Loop over all pairs of atoms, computing the gas-phase
  !         electrostatic energies, the LJ terms, and the off-diagonal
  !         GB terms.  Also accumulate the derivatives of these off-diagonal
  !         terms with respect to the effective radii,
  !         sumdeijda(k) will hold  sum over i, j>i (deij / dak),  where
  !         "ak" is the inverse of the effective radius for atom "k".
  !
  !         Update the forces with the negative derivatives of the
  !         gas-phase terms, plus the derivatives of the explicit
  !         distance dependence in Fgb, i.e. the derivatives of the
  !         GB energy terms assuming that the effective radii are constant.
  !
  !--------------------------------------------------------------------------

  if ( igb /= 6) sumdeijda(1:atm_cnt) = 0.d0

  ! Note: this code assumes that the belly atoms are the first natbel
  !       atoms...this is checked in mdread.

  iexcl = 1

#ifdef MPI
    do i = 1, mytaskid
      iexcl = iexcl + numex(i)
    end do
#endif

#ifdef MPI
  do i = mytaskid + 1, max_i, numtasks
#else
  do i = 1, max_i
#endif
    if (ti_mode .ne. 0) then !non softcore ti
      if (ti_lst(ti_mask_piece,i) .ne. 0) then
      ! We still need to take into account the exclusions for this atom
#ifdef MPI
        do k = i, min(i + numtasks - 1, atm_cnt)
          iexcl = iexcl + numex(k)
        end do
#else
        iexcl = iexcl + numex(i)
#endif
        cycle
      end if
    end if

    xi = crd(3 * i - 2)
    yi = crd(3 * i - 1)
    zi = crd(3 * i)
    qi = charge(i)
    if ( igb /= 6) ri = reff(i)
    iaci = ntypes * (iac(i) - 1)
    jexcl = iexcl
    jexcl_last = iexcl + numex(i) - 1

    dumx = 0.d0
    dumy = 0.d0
    dumz = 0.d0

    ! check the exclusion list for eel and vdw:

    do k = i + 1, atm_cnt
      skipv(k) = .false.
    end do
    do jjv = jexcl, jexcl_last
      skipv(natex(jjv)) = .true.
    end do

    icount = 0
    if (ti_mode .eq. 0) then
      do j = i + 1, atm_cnt
        xij = xi - crd(3 * j - 2)
        yij = yi - crd(3 * j - 1)
        zij = zi - crd(3 * j)
        r2 = xij * xij + yij * yij + zij * zij !calculates distance^2 between
                                               !atoms i and j
#ifdef CUDA
  !CUDA GB Code does not currently support cutoffs in GB.
  !      if (r2 .gt. cut2) cycle
  !      if (.not. onstep .and. r2 .gt. cut_inner2) cycle
#else
        if (r2 .gt. cut2) cycle      !if outside the cutoff, no calculation
        if (.not. onstep .and. r2 .gt. cut_inner2) cycle
#endif

        icount = icount + 1
        jj(icount) = j
        r2x(icount) = r2
        if ( igb /= 6) rjx(icount) = reff(j) !reff undefined at all times for
      end do
    else
      do j = i + 1, atm_cnt
        if (ti_lst(ti_mask_piece,j) .ne. 0) cycle
        xij = xi - crd(3 * j - 2)
        yij = yi - crd(3 * j - 1)
        zij = zi - crd(3 * j)
        r2 = xij * xij + yij * yij + zij * zij
#ifdef CUDA
  !CUDA GB Code does not currently support cutoffs in GB.
  !      if (r2 .gt. cut2) cycle
  !      if (.not. onstep .and. r2 .gt. cut_inner2) cycle
#else
        if (r2 .gt. cut2) cycle
        if (.not. onstep .and. r2 .gt. cut_inner2) cycle
#endif
        icount = icount + 1
        jj(icount) = j
        r2x(icount) = r2
        if ( igb /= 6) rjx(icount) = reff(j)
      end do
    end if

    if (igb /= 6) then
        ! gbion option:
        if (gbion .eq. 1) then
          do k = 1, icount
            vectmp1(k) = 4.d0 * ri * rjx(k) * & 
                   atom_ion_coef(atom_ion_mask(i))*atom_ion_coef(atom_ion_mask(jj(k)))
          end do
        else if (gbion .eq. 2) then
          if (atom_ion_mask(i) .eq. 0) then ! atom i is solute
            do k = 1, icount
              if (atom_ion_mask(jj(k)) .gt. 0) then ! atom j is ion
                vectmp1(k) = 4.d0 * ri * rjx(k) * gi_coef_1
              else  ! atom j is solute
                vectmp1(k) = 4.d0 * ri * rjx(k)
              end if
            end do
          else ! atom i is ion
            do k = 1, icount
              if (atom_ion_mask(jj(k)) .gt. 0) then ! atom j is ion
                vectmp1(k) = 4.d0 * ri * rjx(k) * gi_coef_2
              else ! atom j is solute
                vectmp1(k) = 4.d0 * ri * rjx(k) * gi_coef_1
              end if
            end do 
          end if
        else if (gbion .eq. 3) then
          if (atom_ion_mask(i) .eq. 0) then ! atom i is solute
            ion_mask_pn_i=0
          else if (atom_ion_mask(i) .eq. 3) then ! atom i is anion
            ion_mask_pn_i=1
          else ! atom i is cation
            ion_mask_pn_i=2
          end if
          do k = 1, icount
            if (atom_ion_mask(jj(k)) .eq. 0) then ! atom j is solute
              ion_mask_pn_j=0
            else if (atom_ion_mask(jj(k)) .eq. 3) then ! atom j is anion
              ion_mask_pn_j=1
            else ! atom j is cation
              ion_mask_pn_j=2
            end if
            ion_mask_pn_index = ion_mask_pn_i*3+ion_mask_pn_j
              ! parameters matrix
              ! | 1.0(0)   1_n(1)   1_p(2)  | 
              ! | 1_n(3)   2_nn(4)  2_pn(5) |
              ! | 1_p(6)   2_pn(7)  2_pp(8) |
            if (ion_mask_pn_index .eq. 0) then ! sulote-solute pair
              vectmp1(k) = 4.d0 * ri * rjx(k)
            else if (ion_mask_pn_index .eq. 1 .or. ion_mask_pn_index .eq. 3) then ! solute-anion pair
              vectmp1(k) = 4.d0 * ri * rjx(k) * gi_coef_1_n
            else if (ion_mask_pn_index .eq. 2 .or. ion_mask_pn_index .eq. 6) then ! solute-cation pair
              vectmp1(k) = 4.d0 * ri * rjx(k) * gi_coef_1_p
            else if (ion_mask_pn_index .eq. 5 .or. ion_mask_pn_index .eq. 7) then ! cation-anion pair
              vectmp1(k) = 4.d0 * ri * rjx(k) * gi_coef_2_pn
            else if (ion_mask_pn_index .eq. 4) then ! anion-anion pair
              vectmp1(k) = 4.d0 * ri * rjx(k) * gi_coef_2_nn
            else ! ion_mask_pn_index .eq. 8 -- cation-cation pair
              vectmp1(k) = 4.d0 * ri * rjx(k) * gi_coef_2_pp
            end if
          end do
        else ! (gbion .eq. 0)
          vectmp1(1:icount) = 4.d0 * ri * rjx(1:icount)  !4 * reff(i) * reff(j)
        end if
        call vdinv(icount, vectmp1, vectmp1)
        vectmp1(1:icount) = -r2x(1:icount) * vectmp1(1:icount) !invert born radii
        call vdexp(icount, vectmp1, vectmp1)
    ! vectmp1 now contains exp(-rij^2/[4*ai*aj])
        vectmp3(1:icount) = r2x(1:icount) + rjx(1:icount) * ri * vectmp1(1:icount)
    ! vectmp3 now contains fij
        call vdinvsqrt(icount, vectmp3, vectmp2)
    ! vectmp2 now contains 1/fij
    end if

    if (gb_kappa .ne. 0.d0) then !
      call vdinv(icount, vectmp2, vectmp3)
      vectmp3(1:icount) = -gb_kappa * vectmp3(1:icount)
      call vdexp(icount, vectmp3, vectmp4)
      ! vectmp4 now contains exp(-kappa*fij)
    end if

    call vdinvsqrt(icount, r2x, vectmp5) ! 1/rij

    ! Current vector array contents
    ! vectmp1 = exp(-rij^2/[4*ai*aj])
    ! vectmp2 = 1/fij
    ! vectmp3 = -kappa*fij - if kappa .ne. 0.d0, otherwise .eq. fij
    ! vectmp4 = exp(-kappa*fij)
    ! vectmp5 = 1/rij

    ! Start first outer loop
    ! dir$ ivdep
    do k = 1, icount

      j = jj(k)
      xij = xi - crd(3 * j - 2)
      yij = yi - crd(3 * j - 1)
      zij = zi - crd(3 * j)
      r2 = r2x(k)
      qiqj = qi * charge(j)

      if (igb /= 6) then
        if (gb_kappa .eq. 0.d0) then
          fgbk = 0.d0
          expmkf = extdiel_inv
        else
          expmkf = vectmp4(k) * extdiel_inv
          fgbk = vectmp3(k)*expmkf !-kappa*fij*exp(-kappa*fij)/Eout
          if (alpb .eq. 1) &
            fgbk = fgbk + (fgbk * one_arad_beta * (-vectmp3(k) * gb_kappa_inv))
            ! (-kappa*fij*exp(-kappa*fij)(1 + fij*ab/A)/Eout)*(1/fij+ab/A)
            ! Note: -vectmp2(k)*kappa_inv = fij
        end if

        if (gbion .eq. 2) then ! gb w/ ion, gbion=2
          if ( (atom_ion_mask(i) .eq. 0 .and. atom_ion_mask(j) .gt. 0) &
          .or. (atom_ion_mask(i) .gt. 0 .and. atom_ion_mask(j) .eq. 0) ) then ! ion-solute
            dl = intdiel_inv_ion_1 - expmkf
          else if ( atom_ion_mask(i) .gt. 0 .and. atom_ion_mask(j) .gt. 0 ) then ! ion-ion
            dl = intdiel_inv_ion_2 - expmkf
          else  ! solute-solute
            dl = intdiel_inv - expmkf
          end if
        else if (gbion .eq. 3) then ! gb w/ ion, gbion=3, separate cation and anion
          if (atom_ion_mask(i) .eq. 0) then ! atom i is solute
            ion_mask_pn_i=0
          else if (atom_ion_mask(i) .eq. 3) then ! atom i is anion
            ion_mask_pn_i=1
          else ! atom i is cation
            ion_mask_pn_i=2
          end if
          if (atom_ion_mask(j) .eq. 0) then ! atom j is solute
            ion_mask_pn_j=0
          else if (atom_ion_mask(j) .eq. 3) then ! atom j is anion
            ion_mask_pn_j=1
          else ! atom j is cation
            ion_mask_pn_j=2
          end if
          ion_mask_pn_index = ion_mask_pn_i*3+ion_mask_pn_j
          if (ion_mask_pn_index .eq. 0) then ! sulote-solute pair
            dl = intdiel_inv - expmkf
          else if (ion_mask_pn_index .eq. 1 .or. ion_mask_pn_index .eq. 3) then ! solute-anion pair
            dl = intdiel_inv_ion_1_n - expmkf
          else if (ion_mask_pn_index .eq. 2 .or. ion_mask_pn_index .eq. 6) then ! solute-cation pair
            dl = intdiel_inv_ion_1_p - expmkf
          else if (ion_mask_pn_index .eq. 5 .or. ion_mask_pn_index .eq. 7) then ! cation-anion pair
            dl = intdiel_inv_ion_2_pn - expmkf
          else if (ion_mask_pn_index .eq. 4) then ! anion-anion pair
            dl = intdiel_inv_ion_2_nn - expmkf
          else ! ion_mask_pn_index .eq. 8 -- cation-cation pair
            dl = intdiel_inv_ion_2_pp - expmkf
          end if
        else  ! default gb, gbion!=2 or 3
          dl = intdiel_inv - expmkf
        end if
        fgbi = vectmp2(k) ! 1.d0/fij

        if (alpb .eq. 0) then
          e = -qiqj * dl * fgbi
        else
          e = -qiqj * dl * (fgbi + one_arad_beta)
        end if

        egb = egb + e

        temp4 = fgbi * fgbi * fgbi ! 1.d0/fij^3

        ! [here, and in the gas-phase part, "de" contains -(1/r)(dE/dr)]

        temp6 = -qiqj * temp4 * (dl + fgbk)

        ! -qiqj/fij^3*[1/Ein - e(-Kfij)/Eout) -kappa*fij*
        ! exp(-kappa*fij)(1 + fij*a*b/A ) /Eout]

        temp1 = vectmp1(k) ! exp(-rij^2/[4*ai*aj])

        de = temp6 * (1.d0 - 0.25d0 * temp1)

        rj = rjx(k)

        temp5 = 0.5d0 * temp1 * temp6 * (ri * rj + 0.25d0 * r2)


        sumdeijda(i) = sumdeijda(i) + ri * temp5
        sumdeijda(j) = sumdeijda(j) + rj * temp5
      else
        de = 0.0d0
      end if !igb/=6

      ! GPU-SASA calculations
      if (gbsa .eq. 3) then
        if (epsilonNP(i) /= 0 .and. epsilonNP(j) /= 0) then
          rinv = vectmp5(k) ! 1.d0/rij
          dist = 1.0d0/rinv
          if (dist .lt. radiusNP(i)+radiusNP(j)) then
            orderm = 10
            ordern = 4
            tempsum = (sqrt(epsilonNP(i)*epsilonNP(j)))/(orderm-ordern)
            !paraA=ordern*tempsum*(sigmaNP(i) + sigmaNP(j))**orderm
            !paraB=orderm*tempsum*(sigmaNP(i) + sigmaNP(j))**ordern
            paraA=ordern*tempsum
            paraB=orderm*tempsum
            !Sij = radiusNP(i)+radiusNP(j)+sigmaNP(i)+sigmaNP(j)
            !tempsum1 = paraA*(Sij-dist)**(-1*orderm)
            !tempsum2 = paraB*(Sij-dist)**(-1*ordern)
            !tempsum3 = (Sij-dist)**(-1)
            !tempsum4 = tempsum2 - tempsum1 - sqrt(epsilonNP(i)*epsilonNP(j))
            Cij = radiusNP(i)+radiusNP(j)-dist
            invS = 1.0d0/(sigmaNP(i)+sigmaNP(j))
            CijinvS = 1 + Cij * invS
            tempsum1 = paraA*(CijinvS)**(-1*orderm)
            tempsum2 = paraB*(CijinvS)**(-1*ordern)
            tempsum3 = (CijinvS)**(-1)
            tempsum4 = tempsum2 - tempsum1 - tempsum*(orderm-ordern)
            AgsuminvRn(i) = AgsuminvRn(i) + tempsum4 * 0.6d0
            AgsuminvRn(j) = AgsuminvRn(j) + tempsum4 * 0.6d0

            tempsum5 = (orderm*tempsum1*tempsum3-ordern*tempsum2*tempsum3)*rinv*surften*0.6d0*invS *2.0d0
            desurfi = tempsum5
            desurfj = tempsum5


            dumx = dumx + desurfi * xij
            dumy = dumy + desurfi * yij
            dumz = dumz + desurfi * zij

            !accffromf6x(i) = accffromf6x(i) + desurfi * xij
            !accffromf6y(i) = accffromf6y(i) + desurfi * yij
            !accffromf6z(i) = accffromf6z(i) + desurfi * zij

            frc(3*j-2) = frc(3*j-2) - desurfj * xij
            frc(3*j-1) = frc(3*j-1) - desurfj * yij
            frc(3*j  ) = frc(3*j  ) - desurfj * zij
            !accffromf6x(j) = accffromf6x(j) - desurfj * xij
            !accffromf6y(j) = accffromf6y(j) - desurfj * yij
            !accffromf6z(j) = accffromf6z(j) - desurfj * zij
            !write (*,*) "i,j,dist,desurfi,paraA,paraB",i,j,dist,desurfi,paraA,paraB
          end if !( dist < 5.8-6.4 Angstroms)
        end if ! (tempsum not equal to 0)
      end if ! ( gbsa == 3 )


      ! skip exclusions for remaining terms:

      if (.not. skipv(j)) then

        ! gas-phase Coulomb energy:

        rinv = vectmp5(k) ! 1.d0/rij
        r2inv = rinv * rinv
        if (gbion .eq. 2) then  ! gb w/ ion gbion = 2
          if ( (atom_ion_mask(i) .eq. 0 .and. atom_ion_mask(j) .gt. 0) &
          .or. (atom_ion_mask(i) .gt. 0 .and. atom_ion_mask(j) .eq. 0) ) then ! ion-solute
            eel = intdiel_inv_ion_1 * qiqj * rinv
          else if ( atom_ion_mask(i) .gt. 0 .and. atom_ion_mask(j) .gt. 0 ) then ! ion-ion
            eel = intdiel_inv_ion_2 * qiqj * rinv
          else ! solute-solute
            eel = intdiel_inv * qiqj * rinv
          end if
        else if (gbion .eq. 3) then ! gb w/ ion, gbion=3, separate cation and anion
          if (atom_ion_mask(i) .eq. 0) then ! atom i is solute
            ion_mask_pn_i=0
          else if (atom_ion_mask(i) .eq. 3) then ! atom i is anion
            ion_mask_pn_i=1
          else ! atom i is cation
            ion_mask_pn_i=2
          end if
          if (atom_ion_mask(j) .eq. 0) then ! atom j is solute
            ion_mask_pn_j=0
          else if (atom_ion_mask(j) .eq. 3) then ! atom j is anion
            ion_mask_pn_j=1
          else ! atom j is cation
            ion_mask_pn_j=2
          end if
          ion_mask_pn_index = ion_mask_pn_i*3+ion_mask_pn_j
          if (ion_mask_pn_index .eq. 0) then ! sulote-solute pair
            eel = intdiel_inv * qiqj * rinv
          else if (ion_mask_pn_index .eq. 1 .or. ion_mask_pn_index .eq. 3) then ! solute-anion pair
            eel = intdiel_inv_ion_1_n * qiqj * rinv
          else if (ion_mask_pn_index .eq. 2 .or. ion_mask_pn_index .eq. 6) then ! solute-cation pair
            eel = intdiel_inv_ion_1_p * qiqj * rinv
          else if (ion_mask_pn_index .eq. 5 .or. ion_mask_pn_index .eq. 7) then ! cation-anion pair
            eel = intdiel_inv_ion_2_pn * qiqj * rinv
          else if (ion_mask_pn_index .eq. 4) then ! anion-anion pair
            eel = intdiel_inv_ion_2_nn * qiqj * rinv
          else ! ion_mask_pn_index .eq. 8 -- cation-cation pair
            eel = intdiel_inv_ion_2_pp * qiqj * rinv
          end if
        else ! gb w/o ion, gbion!=2 or 3
          eel = intdiel_inv * qiqj * rinv
        end if
        eelt = eelt + eel
        de = de + eel * r2inv

        ! van der Waals energy:

        ic = ico(iaci + iac(j))
        if (ic .gt. 0) then
          ! 6-12 potential:
          r6inv = r2inv * r2inv * r2inv
          f6 = gbl_cn2(ic) * r6inv
          f12 = gbl_cn1(ic) * (r6inv * r6inv)
          evdw = evdw + (f12 - f6)
          de = de + (12.d0 * f12 - 6.d0 * f6) * r2inv

#ifdef HAS_10_12
          ! The following could be commented out if the Cornell et al.
          ! force field was always used, since then all hbond terms are zero.

        else if (ic .lt. 0 ) then
          ! 10-12 potential:
          r10inv = r2inv * r2inv * r2inv * r2inv * r2inv
          f10 = gbl_bsol(-ic) * r10inv
          f12 = gbl_asol(-ic) * r10inv * r2inv
          evdw = evdw + f12 - f10
          de = de + (12.d0 * f12 - 10.d0 * f10) * r2inv
#endif
        end if  ! (ic .gt. 0)
      end if  ! (.not. skipv(j))

      ! derivatives:

      if (onstep .and. r2 .gt. cut_inner2) then
        de = de * nrespa
      else
        de = de * nrespai
      end if

      dedx = de * xij
      dedy = de * yij
      dedz = de * zij
      dumx = dumx + dedx
      dumy = dumy + dedy
      dumz = dumz + dedz
      frc(3 * j - 2) = frc(3 * j - 2) - dedx
      frc(3 * j - 1) = frc(3 * j - 1) - dedy
      frc(3 * j) = frc(3 * j) - dedz
    end do

    frc(3 * i - 2) = frc(3 * i - 2) + dumx
    frc(3 * i - 1) = frc(3 * i - 1) + dumy
    frc(3 * i) = frc(3 * i) + dumz
#ifdef MPI
    do k = i, min(i + numtasks - 1, atm_cnt)
      iexcl = iexcl + numex(k)
    end do
#else
    iexcl = iexcl + numex(i)
#endif
  end do  !  i = 1, max_i
!#endif
  if (gbsa .eq. 3) then
#ifdef MPI
    call mpi_allreduce(AgsuminvRn, vectmp1, atm_cnt, mpi_double_precision, &
                       mpi_sum, pmemd_comm, err_code_mpi)
    AgsuminvRn(1:atm_cnt) = vectmp1(1:atm_cnt)
#endif

#ifdef MPI
    do i = mytaskid + 1, max_i, numtasks
#else
    do i=1, max_i  !"natom" in sander was replaced by "atm_cnt = max_i"
#endif

      totsasa = totsasa + AgsuminvRn(i)
    end do 
  end if !( gbsa == 3  )

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
  if (igb .ne. 6) then
    if (onstep) then

      neibr_cnt = 0 ! we have no neighbors yet for LCPO gbsa

#ifdef MPI

     ! first, collect all the sumdeijda terms:

      call mpi_allreduce(sumdeijda, vectmp1, atm_cnt, mpi_double_precision, &
                         mpi_sum, pmemd_comm, err_code_mpi)

      sumdeijda(1:atm_cnt) = vectmp1(1:atm_cnt)

      call update_gb_time(dist_gb_rad_timer)
#endif

      frespa = nrespa

      ! diagonal egb term, plus off-diag derivs wrt alpha .eq. reff^-1:

#ifdef MPI
      do i = mytaskid + 1, max_i, numtasks
#else
      do i = 1, max_i
#endif
        if (ti_mode .ne. 0) then
          if (ti_lst(ti_mask_piece,i) .ne. 0) cycle
        end if

        f_xi = 0.d0
        f_yi = 0.d0
        f_zi = 0.d0
        qi = charge(i)
        expmkf = exp(-gb_kappa * reff(i)) * extdiel_inv
        if (gbion .eq. 2) then  ! gb w/ ion, gbion = 2
          if ( atom_ion_mask(i) .gt. 0 ) then ! ion-ion
            dl = intdiel_inv_ion_2 - expmkf
          else ! solute-solute
            dl = intdiel_inv - expmkf
          end if
        else if (gbion .eq. 3) then ! gb w/ ion, gbion=3, separate cation and anion
          if (atom_ion_mask(i) .eq. 0) then ! atom i is solute
            ! ion_mask_pn_i=0
            dl = intdiel_inv - expmkf
          else if (atom_ion_mask(i) .eq. 3) then ! atom i is anion
            ! ion_mask_pn_i=1
            dl = intdiel_inv_ion_2_nn - expmkf
          else ! atom i is cation
            ! ion_mask_pn_i=2
            dl = intdiel_inv_ion_2_pp - expmkf
          end if
        else ! gb w/o ion, gbion != 2 or 3
          dl = intdiel_inv - expmkf
        end if
        qi2h = 0.5d0 * qi * qi
        qid2h = qi2h * dl

        if (alpb .eq. 0) then
          egb = egb - qid2h / reff(i)
          temp7 = -sumdeijda(i) + qid2h - gb_kappa * qi2h * expmkf * reff(i)
        else
          egb = egb - qid2h * (1.d0/reff(i) + one_arad_beta)
          temp7 = -sumdeijda(i) + qid2h - gb_kappa * qi2h * expmkf * reff(i) * &
                  (1.d0 + one_arad_beta * reff(i))
        end if

        xi = crd(3 * i - 2)
        yi = crd(3 * i - 1)
        zi = crd(3 * i)
        ri = rborn(i) - offset
        ri1i = 1.d0 / ri
        iaci = ntypes * (iac(i) - 1)

        if (igb .eq. 2 .or. igb .eq. 5 .or. igb .eq. 7 .or. igb .eq. 8) then
        !entire region is skipped for igb=6

          ! new onufriev: we have to later scale values by a
          !               alpha,beta,gamma -dependent factor:

          ri = rborn(i) - offset
          thi = tanh((gb_alpha_arry(i) + gb_gamma_arry(i) * psi(i) * psi(i) - &
                gb_beta_arry(i) * psi(i)) * psi(i))

          thi2 = (gb_alpha_arry(i) + 3.d0 * gb_gamma_arry(i) * psi(i) * psi(i) - &
                 2.d0 * gb_beta_arry(i) * psi(i)) * (1.d0 - thi * thi) * ri / &
                 rborn(i)
        end if

        icount = 0
        if (ti_mode .eq. 0) then
          do j = 1, atm_cnt
            if (i .eq. j) cycle
            xij = xi - crd(3 * j - 2)
            yij = yi - crd(3 * j - 1)
            zij = zi - crd(3 * j)
            r2 = xij * xij + yij * yij + zij * zij
            if (r2 .gt. rgbmaxpsmax2) cycle

            ! pairlist contains only atoms within rgbmax + safety margin

            icount = icount + 1
            jj(icount) = j
            r2x(icount) = r2
          end do
        else
          do j = 1, atm_cnt
            if (i .eq. j) cycle
            if (ti_lst(ti_mask_piece,j) .ne. 0) cycle
            xij = xi - crd(3 * j - 2)
            yij = yi - crd(3 * j - 1)
            zij = zi - crd(3 * j)
            r2 = xij * xij + yij * yij + zij * zij
            if (r2 .gt. rgbmaxpsmax2) cycle

            ! pairlist contains only atoms within rgbmax + safety margin

            icount = icount + 1
            jj(icount) = j
            r2x(icount) = r2
          end do
        end if

        call vdinvsqrt(icount, r2x, vectmp1)

        kk1 = 0
        do k = 1, icount
          j = jj(k)
          r2 = r2x(k)
          sj =  fs(j) !is sj scaling factor? DJM

          dij1i = vectmp1(k)
          dij = r2 * dij1i
          sj2 = sj * sj

          if (dij .gt. 4.d0 * sj) cycle
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
            vectmp2(kk1) = 1.d0
            vectmp4(kk1) = 1.d0
          end if
        end do

        call vdinv(kk1, vectmp2, vectmp2)
        call vdinv(kk1, vectmp3, vectmp3)
        vectmp4(1:kk1) = vectmp4(1:kk1) * vectmp3(1:kk1)
        call vdln(kk1, vectmp4, vectmp4)

        kk1 = 0
        do k = 1, icount
          j = jj(k)
          j3 = 3 * j
          r2 = r2x(k)
          xij = xi - crd(j3 - 2)
          yij = yi - crd(j3 - 1)
          zij = zi - crd(j3)

          dij1i = vectmp1(k)
          dij = r2 * dij1i
          sj = fs(j)
          if (dij .gt. rgbmax + sj) cycle
          sj2 = sj * sj

          ! datmp will hold (1/r)(dai/dr):

          dij2i = dij1i * dij1i
          dij3i = dij2i * dij1i

          if (dij .gt. rgbmax - sj) then

            temp1 = 1.d0 / (dij - sj)
            datmp = 0.125d0 * dij3i * ((r2 + sj2) * &
                    (temp1 * temp1 - rgbmax2i) - 2.d0 * log(rgbmax * temp1))

          else if (dij .gt. 4.d0 * sj) then

            tmpsd = sj2 * dij2i
            dumbo = te + tmpsd * (tf + tmpsd * (tg + tmpsd * (th + tmpsd * thh)))
            datmp = tmpsd * sj * dij2i * dij2i * dumbo

          else if (dij .gt. ri + sj) then

            kk1 = kk1 + 1
            datmp = vectmp2(kk1) * sj * (-0.5d0 * dij2i + vectmp2(kk1)) + &
                    0.25d0 * dij3i * vectmp4(kk1)

          else if (dij .gt. abs(ri - sj)) then

            kk1 = kk1 + 1
            datmp = -0.25d0 * (-0.5d0 * (r2 - ri * ri + sj2) * &
                    dij3i * ri1i * ri1i + dij1i * vectmp2(kk1) * &
                    (vectmp2(kk1) - dij1i) - dij3i * vectmp4(kk1))

          else if (ri .lt. sj) then

            kk1 = kk1 + 1
            datmp = -0.5d0 * (sj * dij2i * vectmp2(kk1) - &
                    2.d0 * sj * vectmp2(kk1) * vectmp2(kk1) - &
                    0.5d0 * dij3i * vectmp4(kk1))

          else

            kk1 = kk1 + 1
            datmp = 0.d0

          end if  ! (dij .gt. 4.d0 * sj)

          if (igb .eq. 7 .or. igb .eq. 8) then

            if (dij .lt. rborn(i) + rborn(j) + gb_neckcut) then

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

              ! temp1 will be divisor of above fraction * dij
              ! (datmp is deriv * 1/r)

              temp1 = 1.d0 + mdist2 + (0.3d0) * mdist6
              temp1 = temp1 * temp1 * dij

              ! (Note "+" means subtracting derivative, since above
              !     expression has leading "-")

              if (gbion .eq. 2) then  ! gb w/ ion gbion = 2
                if ( (atom_ion_mask(i) .eq. 0 .and. atom_ion_mask(j) .gt. 0) &
                .or. (atom_ion_mask(i) .gt. 0 .and. atom_ion_mask(j) .eq. 0) ) then ! ion-solute
                  datmp = datmp + ((2.d0 * mdist + (9.d0/5.d0) * mdist5) * &
                    neckMaxVal(neck_idx(i), neck_idx(j)) * &
                    gb_neckscale * gb_neckscale_ion_1) / temp1
                else if ( atom_ion_mask(i) .gt. 0 .and. atom_ion_mask(j) .gt. 0 ) then ! ion-ion
                  datmp = datmp + ((2.d0 * mdist + (9.d0/5.d0) * mdist5) * &
                    neckMaxVal(neck_idx(i), neck_idx(j)) * &
                    gb_neckscale * gb_neckscale_ion_2) / temp1
                else ! solute-solute
                  datmp = datmp + ((2.d0 * mdist + (9.d0/5.d0) * mdist5) * &
                    neckMaxVal(neck_idx(i), neck_idx(j)) * &
                    gb_neckscale) / temp1
                end if
              else if (gbion .eq. 3) then ! gb w/ ion, gbion=3, separate cation and anion
                if (atom_ion_mask(i) .eq. 0) then ! atom i is solute
                  ion_mask_pn_i=0
                else if (atom_ion_mask(i) .eq. 3) then ! atom i is anion
                  ion_mask_pn_i=1
                else ! atom i is cation
                  ion_mask_pn_i=2
                end if
                if (atom_ion_mask(j) .eq. 0) then ! atom j is solute
                  ion_mask_pn_j=0
                else if (atom_ion_mask(j) .eq. 3) then ! atom j is anion
                  ion_mask_pn_j=1
                else ! atom j is cation
                  ion_mask_pn_j=2
                end if
                ion_mask_pn_index = ion_mask_pn_i*3+ion_mask_pn_j
                if (ion_mask_pn_index .eq. 0) then ! sulote-solute pair
                  datmp = datmp + ((2.d0 * mdist + (9.d0/5.d0) * mdist5) * &
                    neckMaxVal(neck_idx(i), neck_idx(j)) * &
                    gb_neckscale) / temp1
                else if (ion_mask_pn_index .eq. 1 .or. ion_mask_pn_index .eq. 3) then ! solute-anion pair
                  datmp = datmp + ((2.d0 * mdist + (9.d0/5.d0) * mdist5) * &
                    neckMaxVal(neck_idx(i), neck_idx(j)) * &
                    gb_neckscale * gb_neckscale_ion_1_n) / temp1
                else if (ion_mask_pn_index .eq. 2 .or. ion_mask_pn_index .eq. 6) then ! solute-cation pair
                  datmp = datmp + ((2.d0 * mdist + (9.d0/5.d0) * mdist5) * &
                    neckMaxVal(neck_idx(i), neck_idx(j)) * &
                    gb_neckscale * gb_neckscale_ion_1_p) / temp1
                else if (ion_mask_pn_index .eq. 5 .or. ion_mask_pn_index .eq. 7) then ! cation-anion pair
                  datmp = datmp + ((2.d0 * mdist + (9.d0/5.d0) * mdist5) * &
                    neckMaxVal(neck_idx(i), neck_idx(j)) * &
                    gb_neckscale * gb_neckscale_ion_2_pn) / temp1
                else if (ion_mask_pn_index .eq. 4) then ! anion-anion pair
                  datmp = datmp + ((2.d0 * mdist + (9.d0/5.d0) * mdist5) * &
                    neckMaxVal(neck_idx(i), neck_idx(j)) * &
                    gb_neckscale * gb_neckscale_ion_2_nn) / temp1
                else ! ion_mask_pn_index .eq. 8 -- cation-cation pair
                  datmp = datmp + ((2.d0 * mdist + (9.d0/5.d0) * mdist5) * &
                    neckMaxVal(neck_idx(i), neck_idx(j)) * &
                    gb_neckscale * gb_neckscale_ion_2_pp) / temp1
                end if
              else ! gb w/o ion gbion=0 or w/ ion gbion=1 (old imple, remove later)
                datmp = datmp + ((2.d0 * mdist + (9.d0/5.d0) * mdist5) * &
                      neckMaxVal(neck_idx(i), neck_idx(j)) * &
                      gb_neckscale) / temp1
              end if

            end if ! if (dij < rborn(i) +rborn(j) + gb_neckcut)

          end if ! (igb .eq. 7 .or. igb .eq. 8)

          datmp = -datmp * frespa * temp7

          if (igb .eq. 2 .or. igb .eq. 5 .or. igb .eq. 7 .or. igb .eq. 8) &
            datmp = datmp * thi2

          f_x = xij * datmp
          f_y = yij * datmp
          f_z = zij * datmp
          frc(j3 - 2) = frc(j3 - 2) + f_x
          frc(j3 - 1) = frc(j3 - 1) + f_y
          frc(j3) = frc(j3) + f_z
          f_xi = f_xi - f_x
          f_yi = f_yi - f_y
          f_zi = f_zi - f_z

        end do  !  k = 1, icount

        frc(3 * i - 2) = frc(3 * i - 2) + f_xi
        frc(3 * i - 1) = frc(3 * i - 1) + f_yi
        frc(3 * i) = frc(3 * i) + f_zi

      end do   ! end loop over atom i = mytaskid + 1, max_i, numtasks

      call update_gb_time(calc_gb_diag_timer)

      ! Define neighbor list ineighbor for calculating LCPO areas

! BEGIN GBSA implementation to work with the CPU code.
! The gpu call is located in gb_force
      if (gbsa .eq. 1) then
        call gbsa_ene(crd, frc, esurf ,atm_cnt,jj,r2x,natbel)
      end if
! END GBSA implementation to work with the CPU code.
      if (gbsa .eq. 3) then
         esurf = surften*totsasa
      end if

      call update_gb_time(calc_gb_lcpo_timer)

    end if  !  onstep
  end if  ! if igb/=6 - GB gas phase check
  return

end subroutine gb_ene

!*******************************************************************************
!
! Subroutine:  gb_ene_sc
!
! Description: Calculate forces, energies in "gas phase" for alchemical free
!              energy with softcore potentials.
!
!   igb = 6 is equivalent to no continuum solvent being used. This
!   us to run alchemical free energy calculations in gas phase with
!   softcore potentials (ifsc .gt. 0).
!
!   NOTE: this subroutine is based on the above gb_ene subroutine, with
!   modifications to accommodate softcore potentials and topologies. since
!   softcore only works for igb = 6, steps 1 and 3 are completely removed
!   from this subroutine
!
!   The input coordinates are in the "x" array, and the forces in "f"
!   get updated; energy components are returned in "egb", "eelt" and
!   "evdw".
!
!   Input parameters for the "gas-phase" electrostatic energies are
!   the charges, in the "charge()" array.
!
!   Input parameters for the van der Waals terms are "cn1()" and "cn2()",
!   containing LJ 12-6 parameters, and "asol" and "bsol" containing
!   LJ 12-10 parameters.  (The latter are not used in 1994 and later
!   forcefields.)  The "iac" and "ico" arrays are used to point into
!   these matrices of coefficients.
!
!   The "numex" and "natex" arrays are used to find "excluded" pairs of
!   atoms, for which gas-phase electrostatics and LJ terms are skipped;
!
!   The code also supports a multiple-time-step facility in which:
!
!   Pairs closer than sqrt(cut_inner) are evaluated every nrespai steps, pairs
!   between sqrt(cut_inner) and sqrt(cut) are evaluated every nrespa steps,
!   and pairs beyond sqrt(cut) are ignored
!
!   The forces arising from the derivatives of the GB terms with respect
!   to the effective Born radii are evaluated every nrespa steps.
!
!   The surface-area dependent term is evaluated every nrespa steps.
!
!   The effective radii are only updated every nrespai steps
!
!   (Be careful with the above: what seems to work is dt=0.001,
!    nrespai=2, nrespa=4; anything beyond this seems dangerous.)
!
!   Written 2015-2017, primarily by D.J. Mermelstein, with help
!   from R.C Walker, H. Loeffler, and J.W. Kaus.
!
!*******************************************************************************

subroutine gb_ene_sc(crd, frc, rborn, fs, charge, iac, ico, numex, &
                  natex, atm_cnt, natbel, egb, eelt, evdw, esurf, &
                  scdvdl, sceel, scevdw, irespa, skip_radii_)

  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use timers_mod
  use ti_mod

  implicit none

! Formal arguments:

  double precision      :: crd(*)
  double precision      :: frc(*)
  double precision      :: rborn(*)
  double precision      :: fs(*)
  double precision      :: charge(*)
  integer               :: iac(*)
  integer               :: ico(*)
  integer               :: numex(*)
  integer               :: natex(*)
  integer               :: atm_cnt
  integer               :: natbel
  double precision      :: egb, eelt, evdw, esurf
!#ifdef SCTI
  double precision      :: scdvdl, sceel, scevdw
!#endif
  integer, intent(in)   :: irespa
  logical, optional, intent(in) :: skip_radii_

  logical :: skip_radii ! for converting the optional argument

! Local variables:

  double precision      :: cut2, cut_inner2             !
  double precision      :: extdiel_inv                  !
  double precision      :: intdiel_inv                  !
  double precision      :: ri, rj                       !
  double precision      :: ri1i
  double precision      :: xij, yij, zij                !
  double precision      :: dij1i, dij2i, dij3i          !
  double precision      :: r2                           !
  double precision      :: dij                          !
  double precision      :: sj, sj2                      !
  double precision      :: frespa                       !
  double precision      :: qi, qiqj                     !
  double precision      :: dumx, dumy, dumz             !
  double precision      :: fgbi                         !
  double precision      :: rinv, r2inv, r6inv, r10inv   !
  double precision      :: fgbk                         !
  double precision      :: expmkf                       !
  double precision      :: dl                           !
  double precision      :: de                           !
  double precision      :: e                            !
  double precision      :: temp1                        !
  double precision      :: temp4, temp5, temp6, temp7   !
  double precision      :: eel                          !
  double precision      :: f6, f12, f10                 !
  double precision      :: dedx, dedy, dedz             !
  double precision      :: qi2h, qid2h                  !
  double precision      :: datmp                        !
  double precision      :: thi, thi2                    !
  double precision      :: f_x, f_y, f_z                !
  double precision      :: f_xi, f_yi, f_zi             !
  double precision      :: xi, yi, zi                   !
  double precision      :: dumbo                        !
  double precision      :: tmpsd                        !

  ! Variables needed for smooth integration cutoff in Reff:

  double precision      :: rgbmax1i                     !
  double precision      :: rgbmax2i                     !
  double precision      :: rgbmaxpsmax2                 !

  ! Variables needed for softcore LJ terms in TI:

  double precision      :: sc_vdw                       !
  double precision      :: r6                           !

  ! Variables needed for softcore eel terms in TI/MBAR:
  double precision      :: sc_eel                       !
  double precision      :: sc_eel_denom                 !
  double precision      :: sc_eel_denom_sqrt            !
  double precision      :: mbar_eel_denom               !
  double precision      :: mbar_f6                      !
  double precision      :: sceeorderinv
  double precision      :: ti_sign
  integer               :: num_ti_atms_cntd_i
  integer               :: num_ti_atms_cntd_j
  integer               :: ti_region_i
  integer               :: ti_region_j


  logical       :: i_is_sc
  logical       :: j_is_sc

  ! Scratch variables used for calculating neck correction:

  double precision      ::  mdist
  double precision      ::  mdist2
  double precision      ::  mdist3
  double precision      ::  mdist5
  double precision      ::  mdist6

  ! Stuff for alpb:

  double precision              :: alpb_beta
  double precision              :: one_arad_beta
  double precision              :: gb_kappa_inv
  ! Alpha prefactor for alpb_alpha:
  double precision, parameter   :: alpb_alpha = 0.571412d0

  integer               :: icount
  integer               :: neibr_cnt
  integer               :: i, j, k
  integer               :: kk1
  integer               :: max_i
  integer               :: iaci
  integer               :: iexcl, jexcl
  integer               :: jexcl_last
  integer               :: jjv
  integer               :: ic
  integer               :: j3
  logical               :: onstep

  ! FGB taylor coefficients follow
  ! from A to H :
  ! 1/3 , 2/5 , 3/7 , 4/9 , 5/11
  ! 4/3 , 12/5 , 24/7 , 40/9 , 60/11

  double precision, parameter  :: te = 4.d0 / 3.d0
  double precision, parameter  :: tf = 12.d0 / 5.d0
  double precision, parameter  :: tg = 24.d0 / 7.d0
  double precision, parameter  :: th = 40.d0 / 9.d0
  double precision, parameter  :: thh = 60.d0 / 11.d0

  ! Determine if we want to skip the calculation of the effective radii. There
  ! is an optional argument to force skipping it

  skip_radii = (irespa .gt. 1 .and. mod(irespa, nrespai) .ne. 0)

  ! For gas phase calculations gbradii are not needed so set skip to true.
  if (igb == 6) skip_radii = .true.

  !RCW - What is this code doing? What is skip_radii_?
  if (present(skip_radii_)) &
    skip_radii = skip_radii_ .or. skip_radii

  egb = 0.d0
  eelt = 0.d0
  evdw = 0.d0
  esurf = 0.d0
  sceel = 0.d0
  scevdw = 0.d0

  if (mod(irespa, nrespai) .ne. 0) return

  cut2 = gb_cutoff * gb_cutoff
  cut_inner2 = cut_inner * cut_inner
  onstep = mod(irespa, nrespa) .eq. 0

  if (alpb .eq. 0) then
    ! Standard Still's GB
    extdiel_inv = 1.d0 / extdiel
    intdiel_inv = 1.d0 / intdiel
  else
    ! Sigalov Onufriev ALPB (epsilon-dependent GB):
    alpb_beta = alpb_alpha * (intdiel / extdiel)
    extdiel_inv = 1.d0 / (extdiel * (1.d0 + alpb_beta))
    intdiel_inv = 1.d0 / (intdiel * (1.d0 + alpb_beta))
    one_arad_beta = alpb_beta / arad
    if (gb_kappa .ne. 0.d0) gb_kappa_inv = 1.d0 / gb_kappa
  end if

  max_i = atm_cnt
  if (natbel .gt. 0) max_i = natbel

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

  rgbmax1i = 1.d0 / rgbmax
  rgbmax2i = rgbmax1i * rgbmax1i
  rgbmaxpsmax2 = (rgbmax + gb_fs_max)**2

  !---------------------------------------------------------------------------
  ! Step 1: loop over pairs of atoms to compute the effective Born radii.
  !---------------------------------------------------------------------------

  if (.not. skip_radii) call calc_born_radii(atm_cnt, crd, fs, rborn)

  !--------------------------------------------------------------------------
  !
  ! Step 2: Loop over all pairs of atoms, computing the gas-phase
  !         electrostatic energies, the LJ terms, and the off-diagonal
  !         GB terms.  Also accumulate the derivatives of these off-diagonal
  !         terms with respect to the effective radii,
  !         sumdeijda(k) will hold  sum over i, j>i (deij / dak),  where
  !         "ak" is the inverse of the effective radius for atom "k".
  !
  !         Update the forces with the negative derivatives of the
  !         gas-phase terms, plus the derivatives of the explicit
  !         distance dependence in Fgb, i.e. the derivatives of the
  !         GB energy terms assuming that the effective radii are constant.
  !
  !--------------------------------------------------------------------------

  if ( igb /= 6 ) sumdeijda(1:atm_cnt) = 0.d0

  ! Note: this code assumes that the belly atoms are the first natbel
  !       atoms...this is checked in mdread.

  iexcl = 1

#ifdef MPI
    do i = 1, mytaskid
      iexcl = iexcl + numex(i)
    end do
#endif

    if (ti_mask_piece .eq. 2) then
  !region 1
      sc_vdw = scalpha * clambda !ti_wt_stk(2) !vdw
      sc_eel = scbeta * clambda !ti_wt_stk(2) !+ r2 !eel
      ti_sign = -1.0d0
    else if (ti_mask_piece .eq. 1) then
  !region 2
      sc_vdw = scalpha * (1.d0 - clambda)
      sc_eel = scbeta * (1.d0 - clambda) !+ r2
      ti_sign = 1.0d0
    end if


    num_ti_atms_cntd_i = 0

    sceeorderinv = 0.5d0
#ifdef MPI
    do i = mytaskid + 1, max_i, numtasks
#else
    do i = 1, max_i
#endif
      if (ti_lst(ti_mask_piece,i) .ne. 0) then
#ifdef MPI
        do k = i, min(i + numtasks - 1, atm_cnt)
          iexcl = iexcl + numex(k)
        end do
#else
!if atom is not in this region, count the exclusions, and then skip the rest
        iexcl = iexcl + numex(i)
#endif
        cycle
      end if
      i_is_sc = .false.
      if (ti_mask_piece .eq. 2) then
        ti_region_i = 1
      else
        ti_region_i = 2
      endif
      if (ti_lst(ti_region_i,i) .ne. 0) then
        num_ti_atms_cntd_i = num_ti_atms_cntd_i + 1
        if (ti_sc_lst(ti_atm_lst(ti_region_i, num_ti_atms_cntd_i)) .gt. 0) i_is_sc = .true.
      end if

      xi = crd(3 * i - 2)
      yi = crd(3 * i - 1)
      zi = crd(3 * i)
      qi = charge(i)
      iaci = ntypes * (iac(i) - 1)
      jexcl = iexcl
      jexcl_last = iexcl + numex(i) - 1

      dumx = 0.d0
      dumy = 0.d0
      dumz = 0.d0

    ! check the exclusion list for eel and vdw:

      do k = i + 1, atm_cnt
        skipv(k) = .false.
      end do
      do jjv = jexcl, jexcl_last
        skipv(natex(jjv)) = .true.
      end do

!set num_ti_atms_cntd_j = num_ti_atms_cntd_i to account for any atoms we've
!already calculated for in our upper diagonal loop
      num_ti_atms_cntd_j = num_ti_atms_cntd_i

      icount = 0
      do j = i + 1, atm_cnt
        if (ti_lst(ti_mask_piece,j) .ne. 0) cycle
        j_is_sc = .false.
        if (ti_mask_piece .eq. 2) then
          ti_region_j = 1
        else
          ti_region_j = 2
        endif
        if (ti_lst(ti_region_j,j) .ne. 0) then
          num_ti_atms_cntd_j = num_ti_atms_cntd_j + 1
          if (ti_sc_lst(ti_atm_lst(ti_region_j, num_ti_atms_cntd_j)) .gt. 0) j_is_sc = .true.
        end if
!working through here
        xij = xi - crd(3 * j - 2)
        yij = yi - crd(3 * j - 1)
        zij = zi - crd(3 * j)
        r2 = xij * xij + yij * yij + zij * zij
#ifdef CUDA
  !CUDA GB Code does not currently support cutoffs in GB.
  !      if (r2 .gt. cut2) cycle
  !      if (.not. onstep .and. r2 .gt. cut_inner2) cycle
#else
        if (r2 .gt. cut2) cycle
        if (.not. onstep .and. r2 .gt. cut_inner2) cycle
#endif
        icount = icount + 1
        jj(icount) = j
        r2x(icount) = r2
        !if (ifsc .gt. 0) then
        sc_eel_denom = 1/(sc_eel + r2)
!sceeorderinv = 0.5
        sc_eel_denom_sqrt = sc_eel_denom ** sceeorderinv
        !end if
      end do

      call vdinvsqrt(icount, r2x, vectmp5) ! 1/rij

    ! Current vector array contents
    ! vectmp1 = exp(-rij^2/[4*ai*aj])
    ! vectmp2 = 1/fij
    ! vectmp3 = -kappa*fij - if kappa .ne. 0.d0, otherwise .eq. fij
    ! vectmp4 = exp(-kappa*fij)
    ! vectmp5 = 1/rij

    ! Start first outer loop
    ! dir$ ivdep
      do k = 1, icount

        j = jj(k)
        xij = xi - crd(3 * j - 2)
        yij = yi - crd(3 * j - 1)
        zij = zi - crd(3 * j)
        r2 = r2x(k)
        qiqj = qi * charge(j)
        de = 0.0d0 !because igb == 6

      ! skip exclusions for remaining terms:

        if (.not. skipv(j)) then

        ! gas-phase Coulomb energy:
          rinv = vectmp5(k) ! 1.d0/rij
          r2inv = rinv * rinv
          if (.NOT. (i_is_sc .or. j_is_sc)) then !c-c or c-cv
            eel = intdiel_inv * qiqj * rinv
            eelt = eelt + eel
            de = de + eel * r2inv
            if (ifmbar_lcl .ne. 0 .and. do_mbar) then
              if (ti_mask_piece .eq. 2) then
                do bar_i = 1, bar_states
                  bar_cont(bar_i) = bar_cont(bar_i) + eel * &
                                    (bar_lambda(1,bar_i,1) - 1.0d0 + clambda)
                end do
              else if (ti_mask_piece .eq. 1) then
                do bar_i = 1, bar_states
                  bar_cont(bar_i) = bar_cont(bar_i) + eel * &
                                    (bar_lambda(2,bar_i,1) - clambda)
                end do
              end if
            end if
          else if (XOR(i_is_sc, j_is_sc)) then !sc-c
            eel = qiqj *intdiel_inv * sc_eel_denom_sqrt
            eelt = eelt + eel
            scdvdl = scdvdl + (eel * scbeta * ti_sign * sceeorderinv * &
                     sc_eel_denom)
            de = de - qiqj * sc_eel_denom + sc_eel_denom_sqrt &
                  ** (1 + sceeorder)
            if (ifmbar_lcl .ne. 0 .and. do_mbar) then
              if (ti_mask_piece .eq. 2) then
!for mbar, we need energy at each lambda point - energy at lambda for region 1
!so we have to calculate energy at all lambda points which means recalculating
!any lambda dependent portions with all different lambda values
                do bar_i = 1, bar_states
                  mbar_eel_denom = 1.d0 / (r2 + scbeta * (1.d0 - &
                                   bar_lambda(1,bar_i,1))) ** sceeorderinv
                  bar_cont(bar_i) = bar_cont(bar_i) - eel * (1.d0 - clambda) + &
                                    mbar_eel_denom * bar_lambda(1,bar_i,1)
                end do
              else !ti_mask_piece .eq. 1, because only ti_mode > 1 calls this
                do bar_i = 1, bar_states
                  mbar_eel_denom = 1.d0 / (r2 + scbeta * (1.d0 - &
                                   bar_lambda(2,bar_i,1))) ** sceeorderinv
                  bar_cont(bar_i) = bar_cont(bar_i) - eel * clambda + &
                                    mbar_eel_denom * bar_lambda(2,bar_i,1)
                end do
              end if
            end if
          else !sc-sc
            eel = qiqj *intdiel_inv * sc_eel_denom_sqrt
            sceel = sceel + eel
            de = de + eel * r2inv !this should be force, which is the same
          end if

        ! van der Waals energy:

          ic = ico(iaci + iac(j))
          if (ic .gt. 0) then
            r6inv = r2inv * r2inv * r2inv
            r6 = r2 * r2 * r2
            if (.NOT. (i_is_sc .or. j_is_sc)) then !c-c or c-cv
              f6 = gbl_cn2(ic) * r6inv
              f12 = gbl_cn1(ic) * (r6inv * r6inv)
              evdw = evdw + (f12 - f6)
              de = de + (12.d0 * f12 - 6.d0 * f6) * r2inv
              if (ifmbar_lcl .ne. 0 .and. do_mbar) then
                if (ti_mask_piece .eq. 2) then
                  do bar_i = 1, bar_states
                    bar_cont(bar_i) = bar_cont(bar_i) + (f12 - f6) * &
                                      (bar_lambda(1,bar_i,1) - 1.0d0 + clambda)
                  end do
                else if (ti_mask_piece .eq. 1) then
                  do bar_i = 1, bar_states
                    bar_cont(bar_i) = bar_cont(bar_i) + (f12 - f6) * &
                                      (bar_lambda(2,bar_i,1) - clambda)
                  end do
                end if
              end if
            else if (XOR(i_is_sc, j_is_sc)) then !sc-c
              f6 = 1.0d0/(sc_vdw + r6 * ti_sigma6(ic))
              f12 = f6 * f6
              evdw = evdw + ti_foureps(ic) * ( f12 - f6)
              scdvdl = scdvdl + ti_foureps(ic) * scalpha * &
                       ti_signs(ti_region_i) * f12 * (2.0d0 * f6 - 1.0d0)
              de = de + ti_foureps(ic) * r2 * r2 * f12 * &
                  ti_sigma6(ic) * (12.0d0 * f6 - 6.0d0)
              if (ifmbar_lcl .ne. 0 .and. do_mbar) then
                if (ti_mask_piece .eq. 2) then
!for mbar, we need energy at each lambda point - energy at lambda for region 1
!so we have to calculate energy at all lambda points which means recalculating
!any lambda dependent portions with all different lambda values
                  do bar_i = 1, bar_states
                    bar_cont(bar_i) = bar_cont(bar_i) - &
                      ti_foureps(ic) * ( f12 - f6 ) * (1.0d0 - clambda)
                  end do
                  do bar_i = 1, bar_states
                    mbar_f6 = 1.d0 / (scalpha * (1.0d0 - &
                              bar_lambda(1,bar_i,1)) + r6 + ti_sigma6(ic))
!evdw would be scaled later in pme, but we scale here instead for simplicity
                    bar_cont(bar_i) = bar_cont(bar_i) - evdw * (1.d0 - clambda) + &
                                      ti_foureps(ic) * (mbar_f6 * mbar_f6 - &
                                      mbar_f6) * bar_lambda(1,bar_i,1)
                  end do
                else !ti_mask_piece .eq. 1, because only ti_mode > 1 calls this
                  do bar_i = 1, bar_states
                    bar_cont(bar_i) = bar_cont(bar_i) - &
                      ti_foureps(ic) * ( f12 - f6 ) * clambda
                  end do
                  do bar_i = 1, bar_states
                    mbar_f6 = 1.d0 / (scalpha * bar_lambda(2,bar_i,1) + r6 + &
                              ti_sigma6(ic))
                    bar_cont(bar_i) = bar_cont(bar_i) - evdw * clambda + &
                                      ti_foureps(ic) * (mbar_f6 * mbar_f6 - &
                                      mbar_f6) * bar_lambda(2,bar_i,1)
                  end do
                end if
              end if
            else !sc-sc
              f6 = gbl_cn2(ic) * r6inv
              f12 = gbl_cn1(ic) * (r6inv * r6inv)
              scevdw = scevdw + (f12 - f6)
              de = de + (12.d0 * f12 - 6.d0 * f6) * r2inv
            end if !sc determination
          end if  ! (ic .gt. 0)

        end if  ! (.not. skipv(j))

      ! derivatives:

        if (onstep .and. r2 .gt. cut_inner2) then
          de = de * nrespa
        else
          de = de * nrespai
        end if

        dedx = de * xij
        dedy = de * yij
        dedz = de * zij
        dumx = dumx + dedx
        dumy = dumy + dedy
        dumz = dumz + dedz
        frc(3 * j - 2) = frc(3 * j - 2) - dedx
        frc(3 * j - 1) = frc(3 * j - 1) - dedy
        frc(3 * j) = frc(3 * j) - dedz
      end do

      frc(3 * i - 2) = frc(3 * i - 2) + dumx
      frc(3 * i - 1) = frc(3 * i - 1) + dumy
      frc(3 * i) = frc(3 * i) + dumz
#ifdef MPI
      do k = i, min(i + numtasks - 1, atm_cnt)
        iexcl = iexcl + numex(k)
      end do
#else
      iexcl = iexcl + numex(i)
#endif
    end do  !  i = 1, max_i
  call update_gb_time(calc_gb_offdiag_timer)
  return

end subroutine gb_ene_sc

!*******************************************************************************
!
! Subroutine:  calc_born_radii
!
! Description: <TBS>
!
!*******************************************************************************

subroutine calc_born_radii(atm_cnt, crd, fs, rborn)

  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use timers_mod
  use ti_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: crd(*)
  double precision      :: fs(atm_cnt)
  double precision      :: rborn(atm_cnt)

! Local variables:

  double precision      :: ri, rj
  double precision      :: ri1i, rj1i
  double precision      :: xij, yij, zij
  double precision      :: dij1i, dij2i
  double precision      :: r2
  double precision      :: dij
  double precision      :: si, si2
  double precision      :: sj, sj2
  double precision      :: theta
  double precision      :: uij
  double precision      :: xi, yi, zi
  double precision      :: reff_i
  double precision      :: dumbo
  double precision      :: tmpsd

  ! Variables needed for smooth integration cutoff in Reff:

  double precision      :: rgbmax1i
  double precision      :: rgbmax2i
  double precision      :: rgbmaxpsmax2

  ! Scratch variables used for calculating neck correction:

  double precision      ::  mdist
  double precision      ::  mdist2
  double precision      ::  mdist3
  double precision      ::  mdist6
  double precision      ::  neck

  integer               :: icount

  integer               :: i, j, k
  integer               :: kk1, kk2

  ! atom i, j is cation or anion; parameter index 0-8 in 3*3 matrix
  ! parameter matrix
  ! | 1.0(0)   1_n(1)   1_p(2)  | 
  ! | 1_n(3)   2_nn(4)  2_pn(5) |
  ! | 1_p(6)   2_pn(7)  2_pp(8) |
  integer               :: ion_mask_pn_i, ion_mask_pn_j, ion_mask_pn_index 

  ! FGB taylor coefficients follow
  ! from A to H :
  ! 1/3 , 2/5 , 3/7 , 4/9 , 5/11
  ! 4/3 , 12/5 , 24/7 , 40/9 , 60/11

  double precision, parameter  :: ta = 1.d0 / 3.d0
  double precision, parameter  :: tb = 2.d0 / 5.d0
  double precision, parameter  :: tc = 3.d0 / 7.d0
  double precision, parameter  :: td = 4.d0 / 9.d0
  double precision, parameter  :: tdd = 5.d0 / 11.d0

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

  rgbmax1i = 1.d0 / rgbmax
  rgbmax2i = rgbmax1i * rgbmax1i
  rgbmaxpsmax2 = (rgbmax + gb_fs_max)**2

  reff(1:atm_cnt) = 0.d0

#ifdef MPI
  do i = mytaskid + 1, atm_cnt, numtasks
#else
  do i = 1, atm_cnt
#endif
    if (ti_mode .ne. 0) then
      if (ti_lst(ti_mask_piece,i) .ne. 0) cycle
    end if
    xi = crd(3 * i - 2)
    yi = crd(3 * i - 1)
    zi = crd(3 * i)

    reff_i = reff(i)
    ri = rborn(i) - offset
    ri1i = 1.d0 / ri
    si = fs(i)
    si2 = si * si

    ! Here, reff_i will sum the contributions to the inverse effective
    ! radius from all of the atoms surrounding atom "i"; later the
    ! inverse of its own intrinsic radius will be added in

    icount = 0
    if (ti_mode .eq. 0) then
      do j = i + 1, atm_cnt
        xij = xi - crd(3 * j - 2)
        yij = yi - crd(3 * j - 1)
        zij = zi - crd(3 * j)
        r2 = xij * xij + yij * yij + zij * zij
        if (r2 .gt. rgbmaxpsmax2) cycle
        icount = icount + 1
        jj(icount) = j
        r2x(icount) = r2
      end do
    else
      do j = i + 1, atm_cnt
        if (ti_lst(ti_mask_piece,j) .ne. 0) cycle
        xij = xi - crd(3 * j - 2)
        yij = yi - crd(3 * j - 1)
        zij = zi - crd(3 * j)
        r2 = xij * xij + yij * yij + zij * zij
        if (r2 .gt. rgbmaxpsmax2) cycle
        icount = icount + 1
        jj(icount) = j
        r2x(icount) = r2
      end do
    end if
    call vdinvsqrt(icount, r2x, vectmp1)

    kk1 = 0
    kk2 = 0
    !dir$ ivdep
    do k = 1, icount

      j = jj(k)
      r2 = r2x(k)
      sj = fs(j)

      ! don't fill the remaining vectmp arrays if atoms don't see each other:

      dij1i = vectmp1(k)
      dij = r2 * dij1i
      if (dij .gt. rgbmax + si .and. dij .gt. rgbmax + sj) cycle
      rj = rborn(j) - offset

      if (dij .le. 4.d0 * sj) then
        kk1 = kk1 + 1
        vectmp2(kk1) = dij + sj
        if (dij .gt. ri + sj) then
          vectmp4(kk1) = dij - sj
        else if (dij .gt. abs(ri - sj)) then
          vectmp4(kk1) = ri
        else if (ri .lt. sj) then
          vectmp4(kk1) = sj - dij
        else
          vectmp4(kk1) = 1.d0
        end if
      end if

      if (dij .le. 4.d0 * si) then
        kk2 = kk2 + 1
        vectmp3(kk2) = dij + si
        if (dij .gt. rj + si) then
          vectmp5(kk2) = dij - si
        else if (dij .gt. abs(rj - si)) then
          vectmp5(kk2) = rj
        else if (rj .lt. si) then
          vectmp5(kk2) = si - dij
        else
          vectmp5(kk2) = 1.d0
        end if
      end if

    end do  !  k = 1, icount

    call vdinv(kk1, vectmp2, vectmp2)
    call vdinv(kk2, vectmp3, vectmp3)
    vectmp4(1:kk1) = vectmp2(1:kk1) * vectmp4(1:kk1)
    vectmp5(1:kk2) = vectmp3(1:kk2) * vectmp5(1:kk2)
    call vdln(kk1, vectmp4, vectmp4)
    call vdln(kk2, vectmp5, vectmp5)

    kk1 = 0
    kk2 = 0
    do k = 1, icount

      j = jj(k)
      r2 = r2x(k)

      rj = rborn(j) - offset
      rj1i = 1.d0 / rj
      sj = fs(j)

      sj2 = sj * sj

      xij = xi - crd(3 * j - 2)
      yij = yi - crd(3 * j - 1)
      zij = zi - crd(3 * j)

      dij1i = vectmp1(k)
      dij = r2 * dij1i


      if (dij .le. rgbmax + sj) then

        if ((dij .gt. rgbmax - sj)) then

          uij = 1.d0 / (dij - sj)
! Fortran does not guarantee short circuit logical expressions.
! This if-statement is not standard compliant; it must be compiled with
! optimizations in order to avoid SIGSEGV, invalid memory references.
          if ((iphmd .eq. 0) .or. excludeatoms(j) .eq. 0) then
            reff_i = reff_i - 0.125d0 * dij1i * (1.d0 + 2.d0 * dij *uij + &
                     rgbmax2i * (r2 - 4.d0 * rgbmax * dij - sj2) + &
                     2.d0 * log((dij - sj) * rgbmax1i))
          end if

        else if (dij .gt. 4.d0 * sj) then

          dij2i = dij1i * dij1i
          tmpsd = sj2 * dij2i
          dumbo = ta + tmpsd *  (tb + tmpsd * (tc + tmpsd * (td + tmpsd * tdd)))

! Fortran does not guarantee short circuit logical expressions.
! This if-statement is not standard compliant; it must be compiled with
! optimizations in order to avoid SIGSEGV, invalid memory references.
          if ((iphmd .eq. 0) .or. excludeatoms(j) .eq. 0) then
            reff_i = reff_i - tmpsd * sj * dij2i * dumbo
          end if

          !     ---following are from the Appendix of Schaefer and Froemmel,
          !        J. Mol. Biol. 216:1045-1066, 1990, divided by (4*Pi):

        else if (dij .gt. ri + sj) then

          kk1 = kk1 + 1
! Fortran does not guarantee short circuit logical expressions.
! This if-statement is not standard compliant; it must be compiled with
! optimizations in order to avoid SIGSEGV, invalid memory references.
          if ((iphmd .eq. 0) .or. excludeatoms(j) .eq. 0) then
            reff_i = reff_i - 0.5d0 * (sj / (r2 - sj2) + 0.5d0 * dij1i * &
                   vectmp4(kk1))
          end if

          !-----------------------------------------------------------------

        else if (dij .gt. abs(ri - sj)) then

          kk1 = kk1 + 1
          theta = 0.5d0 * ri1i * dij1i * (r2 + ri * ri - sj2)
! Fortran does not guarantee short circuit logical expressions.
! This if-statement is not standard compliant; it must be compiled with
! optimizations in order to avoid SIGSEGV, invalid memory references.
          if ((iphmd .eq. 0) .or. excludeatoms(j) .eq. 0) then
            reff_i = reff_i - 0.25d0 * (ri1i * (2.d0 - theta) - &
                     vectmp2(kk1) + dij1i * vectmp4(kk1))
          end if

          !-----------------------------------------------------------------

        else if (ri .lt. sj) then

          kk1 = kk1 + 1
! Fortran does not guarantee short circuit logical expressions.
! This if-statement is not standard compliant; it must be compiled with
! optimizations in order to avoid SIGSEGV, invalid memory references.
          if ((iphmd .eq. 0) .or. excludeatoms(j) .eq. 0) then
            reff_i = reff_i - 0.5d0 * (sj / (r2 - sj2) + 2.d0 * ri1i + &
                     0.5d0 * dij1i * vectmp4(kk1))
          end if

          !-----------------------------------------------------------------

        else

          kk1 = kk1 + 1

        end if  ! (dij .gt. 4.d0 * sj)

        if (igb .eq. 7 .or. igb .eq. 8) then

          if (dij .lt. rborn(i) + rborn(j) + gb_neckcut) then
            mdist = dij - neckMaxPos(neck_idx(i), neck_idx(j))
            mdist2 = mdist * mdist
            mdist3 = mdist2 * mdist
            mdist6 = mdist3 * mdist3
            neck = neckMaxVal(neck_idx(i), neck_idx(j)) / &
                   (1.d0 + mdist2 + 0.3d0 * mdist6)
! Fortran does not guarantee short circuit logical expressions.
! This if-statement is not standard compliant; it must be compiled with
! optimizations in order to avoid SIGSEGV, invalid memory references.
            if ((iphmd .eq. 0) .or. excludeatoms(j) .eq. 0) then
              if (gbion .eq. 2) then  ! gb w/ ion, gbion = 2
                if ( (atom_ion_mask(i) .eq. 0 .and. atom_ion_mask(j) .gt. 0) &
                .or. (atom_ion_mask(i) .gt. 0 .and. atom_ion_mask(j) .eq. 0) ) then ! ion-solute
                  reff_i = reff_i - gb_neckscale_ion_1 * gb_neckscale * neck
                else if ( atom_ion_mask(i) .gt. 0 .and. atom_ion_mask(j) .gt. 0 ) then ! ion-ion
                  reff_i = reff_i - gb_neckscale_ion_2 * gb_neckscale * neck
                else ! solute-solute
                  reff_i = reff_i - gb_neckscale * neck
                end if
              else if (gbion .eq. 3) then ! gb w/ ion, gbion=3, separate cation and anion
                if (atom_ion_mask(i) .eq. 0) then ! atom i is solute
                  ion_mask_pn_i=0
                else if (atom_ion_mask(i) .eq. 3) then ! atom i is anion
                  ion_mask_pn_i=1
                else ! atom i is cation
                  ion_mask_pn_i=2
                end if
                if (atom_ion_mask(j) .eq. 0) then ! atom j is solute
                  ion_mask_pn_j=0
                else if (atom_ion_mask(j) .eq. 3) then ! atom j is anion
                  ion_mask_pn_j=1
                else ! atom j is cation
                  ion_mask_pn_j=2
                end if
                ion_mask_pn_index = ion_mask_pn_i*3+ion_mask_pn_j
                if (ion_mask_pn_index .eq. 0) then ! sulote-solute pair
                  reff_i = reff_i - gb_neckscale * neck
                else if (ion_mask_pn_index .eq. 1 .or. ion_mask_pn_index .eq. 3) then ! solute-anion pair
                  reff_i = reff_i - gb_neckscale_ion_1_n * gb_neckscale * neck
                else if (ion_mask_pn_index .eq. 2 .or. ion_mask_pn_index .eq. 6) then ! solute-cation pair
                  reff_i = reff_i - gb_neckscale_ion_1_p * gb_neckscale * neck
                else if (ion_mask_pn_index .eq. 5 .or. ion_mask_pn_index .eq. 7) then ! cation-anion pair
                  reff_i = reff_i - gb_neckscale_ion_2_pn * gb_neckscale * neck
                else if (ion_mask_pn_index .eq. 4) then ! anion-anion pair
                  reff_i = reff_i - gb_neckscale_ion_2_nn * gb_neckscale * neck
                else ! ion_mask_pn_index .eq. 8 -- cation-cation pair
                  reff_i = reff_i - gb_neckscale_ion_2_pp * gb_neckscale * neck
                end if
              else ! gb w/o ion
                reff_i = reff_i - gb_neckscale * neck
              end if
            end if
          end if

        end if

      end if

      ! --- Now the same thing, but swap i and j:

      if (dij .gt. rgbmax + si) cycle

      if (dij .gt. rgbmax - si) then

        uij = 1.d0 / (dij - si)
! Fortran does not guarantee short circuit logical expressions.
! This if-statement is not standard compliant; it must be compiled with
! optimizations in order to avoid SIGSEGV, invalid memory references.
        if ((iphmd .eq. 0) .or. excludeatoms(i) .eq. 0) then
          reff(j) = reff(j) - 0.125d0 * dij1i * (1.d0 + 2.d0 * dij * uij + &
                    rgbmax2i * (r2 - 4.d0 * rgbmax * dij - si2) + &
                    2.d0 * log((dij - si) * rgbmax1i))
        end if

      else if (dij .gt. 4.d0 * si) then

        dij2i = dij1i * dij1i
        tmpsd = si2 * dij2i
        dumbo = ta + tmpsd * (tb + tmpsd * (tc + tmpsd * (td + tmpsd * tdd)))
! Fortran does not guarantee short circuit logical expressions.
! This if-statement is not standard compliant; it must be compiled with
! optimizations in order to avoid SIGSEGV, invalid memory references.
        if ((iphmd .eq. 0) .or. excludeatoms(i) .eq. 0) then
          reff(j) = reff(j) - tmpsd * si * dij2i * dumbo
        end if

      else if (dij .gt. rj + si) then

        kk2 = kk2 + 1
! Fortran does not guarantee short circuit logical expressions.
! This if-statement is not standard compliant; it must be compiled with
! optimizations in order to avoid SIGSEGV, invalid memory references.
        if ((iphmd .eq. 0) .or. excludeatoms(i) .eq. 0) then
          reff(j) = reff(j) - 0.5d0 * (si / (r2 - si2) + &
                    0.5d0 * dij1i * vectmp5(kk2))
        end if

        !-----------------------------------------------------------------

      else if (dij .gt. abs(rj - si)) then

        kk2 = kk2 + 1
        theta = 0.5d0 * rj1i * dij1i * (r2 + rj * rj - si2)
! Fortran does not guarantee short circuit logical expressions.
! This if-statement is not standard compliant; it must be compiled with
! optimizations in order to avoid SIGSEGV, invalid memory references.
        if ((iphmd .eq. 0) .or. excludeatoms(i) .eq. 0) then
          reff(j) = reff(j) - 0.25d0 * (rj1i * (2.d0 - theta) - &
                    vectmp3(kk2) + dij1i * vectmp5(kk2))
        end if

        !-----------------------------------------------------------------

      else if (rj .lt. si) then

        kk2 = kk2 + 1
! Fortran does not guarantee short circuit logical expressions.
! This if-statement is not standard compliant; it must be compiled with
! optimizations in order to avoid SIGSEGV, invalid memory references.
        if ((iphmd .eq. 0) .or. excludeatoms(i) .eq. 0) then
          reff(j) = reff(j) - 0.5d0 * (si / (r2 - si2) + 2.d0 * rj1i + &
                    0.5d0 * dij1i * vectmp5(kk2))
        end if

        !-----------------------------------------------------------------

      else

        kk2 = kk2 + 1

      end if  ! (dij .gt. 4.d0 * si)

      if (igb .eq. 7 .or. igb .eq. 8) then
        if (dij .lt. rborn(j) + rborn(i) + gb_neckcut) then
          mdist = dij - neckMaxPos(neck_idx(j), neck_idx(i))
          mdist2 = mdist * mdist
          mdist3 = mdist2 * mdist
          mdist6 = mdist3 * mdist3
          neck = neckMaxVal(neck_idx(j), neck_idx(i)) / &
                 (1.d0 + mdist2 + 0.3d0 * mdist6)
! Fortran does not guarantee short circuit logical expressions.
! This if-statement is not standard compliant; it must be compiled with
! optimizations in order to avoid SIGSEGV, invalid memory references.
          if ((iphmd .eq. 0) .or. excludeatoms(i) .eq. 0) then
            if (gbion .eq. 2) then  ! gb w/ ion, gbion = 2
              if ( (atom_ion_mask(i) .eq. 0 .and. atom_ion_mask(j) .gt. 0) &
              .or. (atom_ion_mask(i) .gt. 0 .and. atom_ion_mask(j) .eq. 0) ) then ! ion-solute
                reff(j) = reff(j) - gb_neckscale_ion_1 * gb_neckscale * neck
              else if ( atom_ion_mask(i) .gt. 0 .and. atom_ion_mask(j) .gt. 0 ) then ! ion-ion
                reff(j) = reff(j) - gb_neckscale_ion_2 * gb_neckscale * neck
              else ! solute-solute
                reff(j) = reff(j) - gb_neckscale * neck
              end if
            else if (gbion .eq. 3) then ! gb w/ ion, gbion=3, separate cation and anion
              if (atom_ion_mask(i) .eq. 0) then ! atom i is solute
                ion_mask_pn_i=0
              else if (atom_ion_mask(i) .eq. 3) then ! atom i is anion
                ion_mask_pn_i=1
              else ! atom i is cation
                ion_mask_pn_i=2
              end if
              if (atom_ion_mask(j) .eq. 0) then ! atom j is solute
                ion_mask_pn_j=0
              else if (atom_ion_mask(j) .eq. 3) then ! atom j is anion
                ion_mask_pn_j=1
              else ! atom j is cation
                ion_mask_pn_j=2
              end if
              ion_mask_pn_index = ion_mask_pn_i*3+ion_mask_pn_j
              if (ion_mask_pn_index .eq. 0) then ! sulote-solute pair
                reff(j) = reff(j) - gb_neckscale * neck
              else if (ion_mask_pn_index .eq. 1 .or. ion_mask_pn_index .eq. 3) then ! solute-anion pair
                reff(j) = reff(j) - gb_neckscale_ion_1_n * gb_neckscale * neck
              else if (ion_mask_pn_index .eq. 2 .or. ion_mask_pn_index .eq. 6) then ! solute-cation pair
                reff(j) = reff(j) - gb_neckscale_ion_1_p * gb_neckscale * neck
              else if (ion_mask_pn_index .eq. 5 .or. ion_mask_pn_index .eq. 7) then ! cation-anion pair
                reff(j) = reff(j) - gb_neckscale_ion_2_pn * gb_neckscale * neck
              else if (ion_mask_pn_index .eq. 4) then ! anion-anion pair
                reff(j) = reff(j) - gb_neckscale_ion_2_nn * gb_neckscale * neck
              else ! ion_mask_pn_index .eq. 8 -- cation-cation pair
                reff(j) = reff(j) - gb_neckscale_ion_2_pp * gb_neckscale * neck
              end if
            else ! gb w/o ion
              reff(j) = reff(j) - gb_neckscale * neck
            end if
          end if
        end if
      end if

    end do  !  k = 1, icount

    ! we are ending the do-i-loop, reassign the scalar to the original array:

    reff(i) = reff_i
  end do  !  i = 1, atm_cnt


#ifdef MPI
  call update_gb_time(calc_gb_rad_timer)

  ! Collect the (inverse) effective radii from other nodes:

  call mpi_allreduce(reff, vectmp1, atm_cnt, mpi_double_precision, &
                     mpi_sum, pmemd_comm, err_code_mpi)

  reff(1:atm_cnt) = vectmp1(1:atm_cnt)

  call update_gb_time(dist_gb_rad_timer)
#endif

  if (igb .eq. 2 .or. igb .eq. 5 .or. igb .eq. 7 .or. igb .eq. 8) then

    ! apply the new Onufriev "gbalpha, gbbeta, gbgamma" correction:

    do i = 1, atm_cnt
      if (ti_mode .ne. 0) then
        if (ti_lst(ti_mask_piece,i) .ne. 0) cycle
      end if
      ri = rborn(i) - offset
      ri1i = 1.d0 / ri
      psi(i) = -ri * reff(i)
      reff(i) = ri1i - tanh((gb_alpha_arry(i) + gb_gamma_arry(i) * psi(i) * &
                psi(i) - gb_beta_arry(i) * psi(i)) * psi(i)) / rborn(i)

      if (reff(i) .lt. 0.d0) reff(i) = 1.d0/30.d0

      reff(i) = 1.d0 / reff(i)
    end do

  else

    ! "standard" GB, including the "diagonal" term here:

    do i = 1, atm_cnt
      if (ti_mode .ne. 0) then
        if (ti_lst(ti_mask_piece,i) .ne. 0) cycle
      end if
      ri = rborn(i) - offset
      ri1i = 1.d0 / ri
      reff(i) = 1.d0 / (reff(i) + ri1i)
    end do
  end if

  if (rbornstat .eq. 1) then
    do i = 1, atm_cnt
      if (ti_mode .ne. 0) then
        if (ti_lst(ti_mask_piece,i) .ne. 0) cycle
      end if
      gbl_rbave(i) = gbl_rbave(i) + reff(i)
      gbl_rbfluct(i) = gbl_rbfluct(i) + reff(i) * reff(i)
      if (gbl_rbmax(i) .le. reff(i)) gbl_rbmax(i) = reff(i)
      if (gbl_rbmin(i) .ge. reff(i)) gbl_rbmin(i) = reff(i)
    end do
  end if

  call update_gb_time(calc_gb_rad_timer)

  return

end subroutine calc_born_radii

end module gb_ene_mod
! End of Original GB code
