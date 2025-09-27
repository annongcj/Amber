#ifdef _OPENMP_
#include "copyright.i"

#define MAX_NLIST_LIMIT 60000000
!********************************************************************
!
! Module:  gb_ene_hybrid_mod
!
! Description: <TBS>
!
!*******************************************************************************
#include "hybrid_datatypes.i"

module gb_ene_hybrid_mod

use gbl_datatypes_mod
use gbsa_mod
use omp_lib

implicit none

 
        ! Global data definitions:
        double precision, allocatable, save           :: gbl_rbmax(:)
        double precision, allocatable, save           :: gbl_rbmin(:)
        double precision, allocatable, save           :: gbl_rbave(:) 
        double precision, allocatable, save           :: gbl_rbfluct(:)
        double precision, allocatable, save           :: r2x(:),r2x_od(:),rjx_od(:)
        integer,          allocatable, save           :: jj(:),jj_od(:)
        GBFloat, allocatable, save                    :: gbl_asol_sp(:)
        GBFloat, allocatable, save                    :: gbl_bsol_sp(:)

        ! Private data definitions:
        GBFloat, allocatable, save,private            :: data_cn1(:)
        GBFloat, allocatable, save,private            :: data_cn2(:)

        GBFloat, allocatable, save,private            :: fs_lcl(:)
        GBFloat, allocatable, save,private            :: rborn_lcl(:)
        GBFloat, allocatable, save,private            :: charge_lcl(:)

        double precision, allocatable, save, private  :: reff(:)
        GBFloat, allocatable, save, private           :: psi(:)

        double precision, allocatable, save, private  :: sumdeijda(:)
        double precision, allocatable, save, private  :: red_buf(:)
        GBFloat, allocatable, save, private           :: vectmp1(:) 
        GBFloat, allocatable, save, private           :: vectmp2(:) 
        GBFloat, allocatable, save, private           :: vectmp3(:) 
        GBFloat, allocatable, save, private           :: vectmp4(:) 
        GBFloat, allocatable, save, private           :: vectmp5(:) 

        double precision, allocatable, save, private  :: gb_alpha_arry(:) 
        double precision, allocatable, save, private  :: gb_beta_arry(:)
        double precision, allocatable, save, private  :: gb_gamma_arry(:) 

        real, allocatable, save, private              :: skipv(:)
        integer, allocatable, save, private           :: neck_idx(:) 

       ! MPI atom distribution and autobalance variables: 
        integer,save,private                          :: start_atm_offd
        integer,save,private                          :: end_atm_offd
        integer,save,private                          :: start_atm_radii
        integer,save,private                          :: end_atm_radii
        integer,save,private                          :: prev_atm_cnt
        real*8, allocatable,save, private             :: radii_time(:)
        real*8, allocatable,save, private             :: offd_time(:)
        real*8, allocatable,save, private             :: diag_time(:)
        integer, allocatable,save,private             :: radii_num_atoms(:)
        integer, allocatable,save, private            :: offd_num_atoms(:)
        real*8,save,private                           :: last_balance_time
        real*8,save, private                          :: loop_time 
        logical,save, private                         :: imbalance_radii
        logical,save, private                         :: imbalance_offd
        logical,save, private                         :: atmcnt_change
        logical,save, private                         :: time_to_balance
        integer,save, private                         :: wall_s, wall_u
        real*8, save, private                         :: strt_time_ms

       ! timode and exclusions: 
        GBFloat, allocatable,private,save             :: ti_mode_skip(:)
        integer,allocatable,save,private              :: iexcl_arr(:)

       ! Reduction Arrays
        GBFloat, allocatable,save, private            :: frcx_red(:,:)
        GBFloat, allocatable,save, private            :: frcy_red(:,:)
        GBFloat, allocatable,save, private            :: frcz_red(:,:)

       ! Neighbor List
        integer, allocatable, private                 :: count_nlist(:) 
        integer, allocatable, private,save            :: nlist_count(:) 
        integer, private, save                        :: max_nlist_count
        GBFloat, allocatable, private                 :: temp_r2(:)
        integer, allocatable,  private                :: temp_jj(:)
        GBFloat, allocatable, save,private            :: nlist_r2(:)
        integer, allocatable, save, private           :: nlist_jj(:)

        GBFloat, allocatable,save, private            :: crdx(:)
        GBFloat, allocatable,save, private            :: crdy(:)
        GBFloat, allocatable,save, private            :: crdz(:)
        integer, save, private                        :: num_threads

        !Precomputed condition variables
        logical,save,private                          :: igb78, igb2578

        ! gb_neckcut: 2.8d0 (diameter of water) is "correct" value but
        ! larger values give smaller discontinuities at the cut:

#ifdef SPDP
        GBFloat, parameter  :: gb_neckcut = 6.8
#else
        GBFloat, parameter  :: gb_neckcut = 6.8d0
#endif
        GBFloat, parameter  :: ta = ONE / THREE
        GBFloat, parameter  :: tb = TWO / FIVE
        GBFloat, parameter  :: tc = THREE / SEVEN 
        GBFloat, parameter  :: td = FOUR / NINE
        GBFloat, parameter  :: tdd = FIVE / ELEVEN

        GBFloat, parameter  :: te = FOUR / THREE
        GBFloat, parameter  :: tf = (THREE * FOUR) / FIVE 
        GBFloat, parameter  :: tg = (SIX * FOUR) / SEVEN
        GBFloat, parameter  :: th = (TEN * FOUR)  / NINE 
        GBFloat, parameter  :: thh = (TEN * SIX) / ELEVEN


        ! Lookup tables for position (atom separation, r) and value of the maximum
        ! of the neck function for given atomic radii ri and rj. Values of neck
        ! maximum are already divided by 4*Pi to save time. Values are give
        ! for each 0.05 angstrom between 1.0 and 2.0 (inclusive), so map to index
        ! with dnint((r-1.0)*20)).  Values were numerically determined in
        ! Mathematica; note FORTRAN column-major array storage, so the data below
        ! may be transposed from how you might expect it.

        GBFloat, parameter  :: neckMaxPos(0:20,0:20) = reshape((/ &
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

        GBFloat, parameter  :: neckMaxVal(0:20,0:20) = reshape((/ &
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
        ! Subroutine:  final_gb_setup_hybrid
        !
        ! Description: <TBS>
        !
        !*******************************************************************************

subroutine final_gb_setup_hybrid(atm_cnt, crd,num_ints, num_reals, my_igb)

        use gbl_constants_mod
        use mdin_ctrl_dat_mod
        use parallel_dat_mod
        use pmemd_lib_mod
        use prmtop_dat_mod

        implicit none

        ! Formal arguments:

        integer, intent(in)           :: atm_cnt
        double precision              :: crd(*)
        integer, intent(in)           :: my_igb

        ! num_ints and num_reals are used to return allocation counts. Don't zero.

        integer, intent(in out)       :: num_ints, num_reals

        ! Local variables:

        integer                       :: alloc_failed
        integer                       :: i
        integer                       :: atomicnumber
        integer                       :: outer_i,j , temp_icount
        double precision              :: rgbmaxpsmax2
        double precision              :: r2, xi,yi,zi, xij,yij,zij 
        integer                       :: pos, o1, o2
        integer                       :: total_nlist_count

        ! Getting number of openmp threads 
        !$omp parallel 
             num_threads = omp_get_num_threads()
        !$omp end parallel 

        allocate(reff(atm_cnt), &
               sumdeijda(atm_cnt), &
               psi(atm_cnt), &
               jj(atm_cnt), &
               r2x(atm_cnt), &
               jj_od(atm_cnt),&
               r2x_od(atm_cnt),&
               rjx_od(atm_cnt),&
               vectmp1(atm_cnt), &
               red_buf(atm_cnt), &
               vectmp2(DIAG_BATCH_SIZE), &
               vectmp3(DIAG_BATCH_SIZE), &
               vectmp4(DIAG_BATCH_SIZE), &
               vectmp5(BATCH_SIZE), &
               iexcl_arr(atm_cnt), &
               frcx_red(atm_cnt,num_threads), &
               frcy_red(atm_cnt,num_threads), &
               frcz_red(atm_cnt,num_threads), &
               count_nlist(atm_cnt), &
               nlist_count(atm_cnt), &
               crdx(atm_cnt), &
               crdy(atm_cnt), &
               crdz(atm_cnt), &
               ti_mode_skip(atm_cnt), &
               skipv(0:atm_cnt), &
               gb_alpha_arry(atm_cnt), &
               gb_beta_arry(atm_cnt), &
               gb_gamma_arry(atm_cnt), &
               fs_lcl(atm_cnt), &
               rborn_lcl(atm_cnt), &
               charge_lcl(atm_cnt), &
               data_cn1(ntypes:(ntypes+1)*(ntypes+1)), &
               data_cn2(ntypes:(ntypes+1)*(ntypes+1)), &
               radii_time(numtasks), &
               diag_time(numtasks), &
               radii_num_atoms(numtasks), &
               offd_num_atoms(numtasks), &
               offd_time(numtasks), &
              stat = alloc_failed)

        if (alloc_failed .ne. 0) call setup_alloc_error

        num_reals = num_reals + size(reff) + &
        size(psi) + &
        size(r2x) + &
        size(crdx) + &
        size(crdy) + &
        size(crdz) + &
        size(skipv) + &  ! skipv now allocated as real
        size(r2x_od) + &
        size(rjx_od) + &
        size(fs_lcl) + &
        size(vectmp1) + &
        size(vectmp2) + &
        size(vectmp3) + &
        size(vectmp4) + &
        size(vectmp5) + &
        size(red_buf) + &
        size(data_cn1) + &
        size(data_cn2) + &
        size(frcx_red) + &
        size(frcy_red) + &
        size(frcz_red) + &
        size(sumdeijda) + &
        size(rborn_lcl) + &
        size(offd_time) + &
        size(diag_time) + &
        size(radii_time) + &
        size(charge_lcl) + &
        size(gb_beta_arry) + &
        size(gb_alpha_arry) + &
        size(gb_gamma_arry) + &
        size(ti_mode_skip)

        num_ints = num_ints + size(iexcl_arr) + &
        size(jj) + &
        size(jj_od) + &
        size(nlist_count) + &
        size(offd_num_atoms) + &
        size(radii_num_atoms) + &
        size(count_nlist)

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


        gbl_rbmax(:) = 0.0
        gbl_rbmin(:) = 999.0
        gbl_rbave(:) = 0.0
        gbl_rbfluct(:) = 0.0

        end if

        if (my_igb .eq. 2 .or. my_igb .eq. 5 .or. my_igb .eq. 7 .or. my_igb .eq.  8) then
              igb2578 = .true.
        end if

        if (my_igb .eq. 7 .or. my_igb .eq. 8) then
             igb78 = .true.

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
               '1.0 - 2.0 angstrom for igb=7. ', &
               'Regenerate prmtop with bondi radii.'
             end if

            call mexit(mdout, 1)

          end if

        end do

        end if

        ! Since the development of igb = 8, the gb_alpha, gb_beta, and gb_gamma
        ! parameters have to be stored in arrays (one for each atom), since that
        ! model introduced atom-dependent values for those parameters. We fill
        ! those arrays here. If we're using igb .eq. 1, 2, 5, or 7, then we'll fill
        ! every value with gb_alpha, gb_beta, and gb_gamma set in init_mdin_ctrl_dat

        if (my_igb .eq. 1 .or. my_igb .eq. 2 .or. my_igb .eq. 5 .or. &
                        my_igb .eq. 7) then

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
      !Hai Nguyen: update GBNeck2nu
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
            gb_alpha_arry(i) = 1.d0
            gb_beta_arry(i) = 0.8d0
            gb_gamma_arry(i) = 4.851d0
          end if
        end if !testing if atom belongs to nuc

        end do

        end if

        ! Set up the GB/SA data structures if we're doing a GB/SA calculation

        if (gbsa .eq. 1) then
                call gbsa_setup(atm_cnt, num_ints, num_reals)
        end if ! (gbsa .eq. 1)

        ! We needed atm_isymbl, etc. for igb=8 and gbsa initialization, but now
        ! we have no more use for it. So deallocate it.

        num_ints = num_ints - size(atm_isymbl) - &
        size(atm_atomicnumber)
        if (allocated(atm_isymbl)) deallocate(atm_isymbl)
        if (allocated(atm_isymbl)) deallocate(atm_atomicnumber)

        rgbmaxpsmax2 = (rgbmax + gb_fs_max)**2
        do outer_i = 1, atm_cnt, BLKSIZE
           do i = outer_i, MIN0(outer_i+BLKSIZE, atm_cnt+1)-1
            xi = crd(3*i - 2)
            yi = crd(3*i -1 )
            zi = crd(3*i)
      
            temp_icount = 0 
            do j=1, atm_cnt
              xij = xi - crd(3*j - 2)
              yij = yi - crd(3*j - 1)
              zij = zi - crd(3*j)
              r2 = xij * xij + yij * yij + zij * zij
              if (i .eq. j) cycle
              if (r2 .gt. rgbmaxpsmax2) cycle
              temp_icount = temp_icount + 1
            end do
            count_nlist(i) = temp_icount
          end do
        end do

       !************************************************* 
       ! Checking whether neighbour list can be reused.
       ! Allocating suitable arrays if it can be reused.
       !************************************************* 

       max_nlist_count = CEILING(1.3 * MAXVAL(count_nlist(1:atm_cnt))/64.0)*64
       if (max_nlist_count .gt. atm_cnt)  max_nlist_count = atm_cnt - 1
       total_nlist_count = max_nlist_count * atm_cnt

       if (total_nlist_count < MAX_NLIST_LIMIT) then
          allocate(nlist_r2(max_nlist_count*atm_cnt), &
                   nlist_jj(max_nlist_count*atm_cnt), &
                   temp_jj(max_nlist_count), &
                   temp_r2(max_nlist_count), &
                   stat = alloc_failed)
          if (alloc_failed .ne. 0) call setup_alloc_error

          num_reals = num_reals + &
          size(nlist_r2) + & 
          size(temp_r2)

          num_ints = num_ints + &
          size(nlist_jj) + &
          size(temp_jj)

          calc_nlist = .true. ! turning neighbour list reuse off.

      else 
        allocate(temp_jj(atm_cnt), &
                 temp_r2(atm_cnt), &
                 stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error

        num_reals = num_reals + &
        size(temp_r2) 

        num_ints = num_ints + &
        size(temp_jj)

        calc_nlist = .false. ! turning neighbour list reuse on.

      end if !  MAX_NLIST_LIMIT

      ! Initializing data_cn1,2 arrays to avoid indirect memory access gbl_cn1,2
      ! arrays.

      data_cn1(:) = ZERO 
      data_cn2(:) = ZERO 
      do i = 1, ntypes
         do j = 1, i
            pos = typ_ico(ntypes * (i -1) + j)
            o1  = i * (ntypes) + j
            o2  = j * (ntypes) + i
            data_cn1(o1) = gbl_cn1(pos) * THREE * FOUR
            data_cn1(o2) = gbl_cn1(pos) * THREE* FOUR
            data_cn2(o1) = gbl_cn2(pos) * SIX
            data_cn2(o2) = gbl_cn2(pos) * SIX 
         end do ! j = 1, i
      end do !do i = 1, ntypes

      ! Initializing atom distribution for load balancing

      atmcnt_change = .true.
      imbalance_radii = .false.
      imbalance_offd = .false.
      time_to_balance = .false.
      radii_time(1:numtasks) = 0.d0
      offd_time(1:numtasks) = 0.d0
      diag_time(1:numtasks) = 0.d0      
      radii_num_atoms(1:numtasks) = 0
      offd_num_atoms(1:numtasks) = 0
      prev_atm_cnt = atm_cnt
      call get_wall_time(wall_s, wall_u)
      last_balance_time = dble(wall_s) * 1000.d0 + dble(wall_u) / 1000.d0
      call distribute_atoms(atm_cnt)

   return

   end subroutine final_gb_setup_hybrid

   !*******************************************************************************
   !
   ! Subroutine:  gb_cleanup_hybrid
   !
   ! Description: Deallocate GB arrays etc.
   !
   !*******************************************************************************

   subroutine gb_cleanup_hybrid(num_ints,num_reals,my_igb)
   
     use mdin_ctrl_dat_mod, only : rbornstat
     use pmemd_lib_mod, only : setup_dealloc_error
   
     implicit none
   
     ! num_ints and num_reals are used to return allocation counts. Don't zero.
     integer, intent(in out)       :: num_ints, num_reals
     integer, intent(in)           :: my_igb
   
     integer cleanup_alloc_failed
   
     num_reals = num_reals - size(reff) - &
        size(psi) - &
        size(r2x) - &
        size(crdx) - &
        size(crdy) - &
        size(crdz) - &
        size(skipv) - &  ! skipv now allocated as real
        size(r2x_od) - &
        size(rjx_od) - &
        size(fs_lcl) - &
        size(vectmp1) - &
        size(vectmp2) - &
        size(vectmp3) - &
        size(vectmp4) - &
        size(vectmp5) - &
        size(red_buf) - &
        size(data_cn1) - &
        size(data_cn2) - &
        size(frcx_red) - &
        size(frcy_red) - &
        size(frcz_red) - &
        size(sumdeijda) - &
        size(rborn_lcl) - &
        size(offd_time) - &
        size(diag_time) - &
        size(radii_time) - &
        size(charge_lcl) - &
        size(gb_beta_arry) - &
        size(gb_alpha_arry) - &
        size(gb_gamma_arry) - &
        size(ti_mode_skip)

        num_ints = num_ints - size(iexcl_arr) - &
        size(jj) - &
        size(jj_od) - &
        size(nlist_count) - &
        size(offd_num_atoms) - &
        size(radii_num_atoms) - &
        size(count_nlist)

        deallocate(reff, &
                sumdeijda, &
                psi, &
                jj, &
                r2x, &
                jj_od, &
                r2x_od, &
                rjx_od, &
                vectmp1, &
                red_buf, &
                vectmp2, & 
                vectmp3, & 
                vectmp4, & 
                vectmp5, & 
                iexcl_arr, &
                frcx_red, & 
                frcy_red, & 
                frcz_red, & 
                count_nlist, &
                nlist_count, &
                crdx, &
                crdy, &
                crdz, &
                ti_mode_skip, &
                skipv, &
                gb_alpha_arry, &
                gb_beta_arry, &
                gb_gamma_arry, &
                fs_lcl, &
                rborn_lcl, &
                charge_lcl, &
                data_cn1, &
                data_cn2, &
                radii_time, &
                diag_time, &
                radii_num_atoms, &
                offd_num_atoms, &
                offd_time, &
                stat = cleanup_alloc_failed)

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
   
     return
   
   end subroutine gb_cleanup_hybrid

   !*******************************************************************************
   !
   ! Subroutine:  gb_ene_hybrid
   !
   ! Description: Calculate forces, energies based on Generalized Born.
   !
   !   This is OpenMP+MPI implementation of Generalized Born algorithm.
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

!************************************************************************
! Using multiple #defines to generate suitable gb_ene_hybrid subroutines
! such as support for ti_mode, softcore etc.
!************************************************************************

#define SOFTCORE_TI
subroutine gb_ene_hyb_energy_ifsc(crd, frc, rborn, fs, charge, iac, ico, numex, &
                 natex, atm_cnt, natbel, egb, eelt, evdw, esurf, irespa, &
                 skip_radii_)
#include "gb_ene_hybrid.i"
end subroutine gb_ene_hyb_energy_ifsc
#undef SOFTCORE_TI

#define GB_ENERGY
#define TIMODE 
subroutine gb_ene_hyb_energy_timode(crd, frc, rborn, fs, charge, iac, ico, numex, &
                 natex, atm_cnt, natbel, egb, eelt, evdw, esurf, irespa, &
                 skip_radii_)
#include "gb_ene_hybrid.i"
end subroutine gb_ene_hyb_energy_timode

#undef TIMODE 
subroutine gb_ene_hyb_energy(crd, frc, rborn, fs, charge, iac, ico, numex, &
                 natex, atm_cnt, natbel, egb, eelt, evdw, esurf, irespa, &
                 skip_radii_)
#include "gb_ene_hybrid.i"
end subroutine gb_ene_hyb_energy

#undef GB_ENERGY
subroutine gb_ene_hyb_force(crd, frc, rborn, fs, charge, iac, ico, numex, &
                 natex, atm_cnt, natbel, egb, eelt, evdw, esurf, irespa, &
                 skip_radii_)
#include "gb_ene_hybrid.i"
end subroutine gb_ene_hyb_force

#define TIMODE 
subroutine gb_ene_hyb_force_timode(crd, frc, rborn, fs, charge, iac, ico, numex, &
                 natex, atm_cnt, natbel, egb, eelt, evdw, esurf, irespa, &
                 skip_radii_)
#include "gb_ene_hybrid.i"
end subroutine gb_ene_hyb_force_timode
#undef TIMODE 


!************************************************************************
!
! Subroutine:  check_imbalance 
!
! Description: 
!
!   Checks whether there is imbalance in wall time taken by the different
!   MPI ranks for the 3 computations (radii, off-diagonal and diagonal) 
!   The % difference that decides if there is imbalance between ranks is
!   now hardcoded to 3% for off-diagonal computation and 5% for sum of 
!   radii + diagonal
!************************************************************************

subroutine check_imbalance()
  use parallel_dat_mod
!   double precision        :: average_radii_time
!   double precision        :: average_diag_time
   double precision        :: average_RD_time
   double precision        :: average_offd_time
   double precision        :: time_diff
   double precision        :: tmp1(numtasks)
   integer                 :: i 

  imbalance_radii = .false.
  imbalance_offd = .false.
  call mpi_allreduce(radii_time, tmp1, numtasks, mpi_double_precision, &
                        mpi_sum, pmemd_comm, err_code_mpi)
  radii_time(1:numtasks) = tmp1(1:numtasks)

  call mpi_allreduce(offd_time, tmp1, numtasks, mpi_double_precision, &
                        mpi_sum, pmemd_comm, err_code_mpi)
  offd_time(1:numtasks) = tmp1(1:numtasks)

  call mpi_allreduce(diag_time, tmp1, numtasks, mpi_double_precision, &
                        mpi_sum, pmemd_comm, err_code_mpi)
  diag_time(1:numtasks) = tmp1(1:numtasks)

!  average_radii_time = sum(radii_time(1:numtasks))/DBLE(numtasks)
!  average_diag_time = sum(diag_time(1:numtasks))/DBLE(numtasks)
 average_RD_time = sum(radii_time(1:numtasks))/DBLE(numtasks) + &
                    sum(diag_time(1:numtasks))/DBLE(numtasks)
  average_offd_time = sum(offd_time(1:numtasks))/DBLE(numtasks)

   if(average_RD_time .gt. 0.d0 .and.  average_offd_time .gt. 0.d0) then
      do i=1,numtasks
        time_diff =  ABS((radii_time(i) + diag_time(i)- average_RD_time)/average_RD_time)
        if (time_diff .ge. 0.05) imbalance_radii = .true.
        time_diff =  ABS((offd_time(i) - average_offd_time)/average_offd_time)
        if (time_diff .ge. 0.03) imbalance_offd = .true.
      end do 
   end if

  if(.not. imbalance_radii) then 
       radii_time(1:numtasks) = 0.d0
       diag_time(1:numtasks) = 0.d0
  end if
  if(.not. imbalance_offd) offd_time(1:numtasks) = 0.d0
   
end subroutine check_imbalance
    

!************************************************************************
!
! Subroutine:  distribute_atoms 
!
! Description: 
!
!   The three computations loops in in gb_ene are the radii, off-diagonal 
!   and diagonal.
!   This subroutine distributes range of atoms taken by the outer loop of
!   the the 3 computations loops to each MPI rank.
!   When distributing atoms for the first time, it is assumed that each
!   rank has equal compute power.
!   Atoms are equally divided to for radii and diagonal computations. 
!   For off-diagonal the calculation is triangular (inner loop starts 
!   from j+1) and hence distribution is based on equal i-j interactions.
!   When this routine is invoked again, the distrubution is based on 
!   time taken by each rank to process 1 outer atom. So this means that
!   faster ranks will be given more atoms. This is useful when ranks
!   reside in non-homogeneous clusters. Note that since same range of
!   atoms are used by radii and diagonal computation, the time taken
!   per atom is sum of radii and diagonal wall times
!
!************************************************************************

subroutine distribute_atoms(atm_cnt)
   use parallel_dat_mod
      implicit none
   ! Formal arguments
     integer, intent(in)           :: atm_cnt
   ! local variables
   double precision    :: num_interactions_atm_cnt
   double precision    :: a,b,c  ! Quadratic constants 
   double precision    :: d, root 
   integer             :: i,j 
   integer             :: task_per_node
   double precision    :: time_per_atm(numtasks)
   integer             :: tmp1(numtasks)
   double precision    :: atm_cnt_dble


   if(atmcnt_change) then
      atm_cnt_dble = dble(atm_cnt)
      num_interactions_atm_cnt = (atm_cnt_dble * (atm_cnt_dble - 1.d0))/2.d0

      a = 1
      b = -1
      c =  -2 * ((numtasks - (mytaskid)) * num_interactions_atm_cnt / numtasks)
      d = b*b - 4 * a * c
      d = sqrt(d)
      root = (-b + d)/(2.0*a)     ! first root __ Only consider positive root

      start_atm_offd = (atm_cnt - root) + 1 

      a = 1
      b = -1
      c =  -2 * ((numtasks - (mytaskid + 1)) * num_interactions_atm_cnt / numtasks) 
      d = b*b - 4 * a * c
      d = sqrt(d)
      root = (-b + d)/(2.0*a)     ! first root __ Only consider positive root

      end_atm_offd = (atm_cnt - root ) 
      if(end_atm_offd .eq. (atm_cnt - 1)) end_atm_offd = atm_cnt

      task_per_node = atm_cnt / numtasks
      start_atm_radii = task_per_node * mytaskid + 1 
      end_atm_radii = task_per_node * (mytaskid+1)
        
      if(mytaskid .eq. (numtasks-1)) end_atm_radii = atm_cnt
      radii_num_atoms(mytaskid + 1) = (end_atm_radii - start_atm_radii ) + 1 
      offd_num_atoms(mytaskid + 1) = (end_atm_offd - start_atm_offd ) + 1 
#ifdef GBtimer
      print *, "atm count change"
      print *, "taskid=",mytaskid,"radii range",start_atm_radii,end_atm_radii
      print *, "taskid=",mytaskid,"offd range",start_atm_offd,end_atm_offd
#endif
    end if

   if(imbalance_radii) then
#ifdef GBtimer
      print *, "radii before balance", mytaskid, start_atm_radii, end_atm_radii
#endif
      call mpi_allreduce(radii_num_atoms, tmp1, numtasks, mpi_integer, &
                        mpi_sum, pmemd_comm, err_code_mpi)
      radii_num_atoms(1:numtasks) = tmp1(1:numtasks)

      do j=1,numtasks
#ifdef GBtimer
         if (master) then
             print *,"rank", j,"radii_time=",radii_time(j),"atoms=",&
                     radii_num_atoms(j)
         end if
#endif
         if (radii_num_atoms(j) .le. 0) radii_num_atoms(j) = 1
         time_per_atm(j) = (radii_time(j)+diag_time(j))/DBLE(radii_num_atoms(j))
         if (time_per_atm(j) .le. 0) time_per_atm(j) = 99999999.0d0
      end do
     
      radii_num_atoms(1:numtasks) = 0 
      end_atm_radii = 0 

      do i=1,mytaskid+1
        root = 0.d0
        do j=1,numtasks
           root = root +  time_per_atm(i)/time_per_atm(j)
        end do
        start_atm_radii = end_atm_radii + 1 
        end_atm_radii =  start_atm_radii + floor(DBLE(atm_cnt)/root) -1
      end do

      if(mytaskid .eq. (numtasks-1)) end_atm_radii = atm_cnt
      radii_num_atoms(mytaskid + 1) = (end_atm_radii - start_atm_radii ) + 1 
#ifdef GBTimer
      print *, "radii balanced", mytaskid, start_atm_radii, end_atm_radii
#endif
      radii_time(1:numtasks) = 0.d0
      diag_time(1:numtasks) = 0.d0
   end if

   if(imbalance_offd) then
#ifdef GBTimer
      print *, "offd before balance", mytaskid, start_atm_offd, end_atm_offd
#endif
      call mpi_allreduce(offd_num_atoms, tmp1, numtasks, mpi_integer, &
                        mpi_sum, pmemd_comm, err_code_mpi)
      offd_num_atoms(1:numtasks) = tmp1(1:numtasks)

      do j=1,numtasks
#ifdef GBTimer
         if (master) then
             print *,"rank", j,"offd_time=",offd_time(j),"atoms=",&
              offd_num_atoms(j)
         end if
#endif
         if (offd_num_atoms(j) .le. 0) offd_num_atoms(j) = 1
         time_per_atm(j) = offd_time(j)/DBLE(offd_num_atoms(j))
         if (time_per_atm(j) .le. 0) time_per_atm(j) = 99999999.0d0
      end do
     
      offd_num_atoms(1:numtasks) = 0 
      end_atm_offd = 0 
      do i=1,mytaskid+1
        root = 0.d0
        do j=1,numtasks
           root = root +  time_per_atm(i)/time_per_atm(j)
        end do
        start_atm_offd = end_atm_offd + 1 
        end_atm_offd =  start_atm_offd + floor(DBLE(atm_cnt)/root) -1
      end do
      if(mytaskid .eq. (numtasks-1)) end_atm_offd = atm_cnt
      offd_num_atoms(mytaskid + 1) = (end_atm_offd - start_atm_offd ) + 1 
#ifdef GBTimer
      print *, "offd balanced", mytaskid, start_atm_offd, end_atm_offd
#endif
      offd_time(1:numtasks) = 0.d0
   end if

end subroutine

end module gb_ene_hybrid_mod
! End of gb_ene hybrid code.
#endif /* gb_ene_hybrid works only in _OPENMP_ currently */
