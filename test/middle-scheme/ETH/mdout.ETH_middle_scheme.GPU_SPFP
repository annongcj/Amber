
          -------------------------------------------------------
          Amber 20 PMEMD                              2020
          -------------------------------------------------------

| PMEMD implementation of SANDER, Release 18

|  Compiled date/time: Fri Mar 11 15:03:13 2022
| Run on 03/12/2022 at 13:14:15

|   Executable path: /home/zhf/work/amberfold/amber20/bin/pmemd.cuda_SPFP
| Working directory: /home/zhf/work/amberfold/amber20/test/middle-scheme/ETH
|          Hostname: Unknown
  [-O]verwriting output

File Assignments:
|   MDIN: mdin                                                                  
|  MDOUT: mdout.ETH_middle_scheme                                               
| INPCRD: inpcrd                                                                
|   PARM: prmtop                                                                
| RESTRT: restrt                                                                
|   REFC: refc                                                                  
|  MDVEL: mdvel                                                                 
|   MDEN: mden                                                                  
|  MDCRD: mdcrd                                                                 
| MDINFO: mdinfo                                                                
|  MDFRC: mdfrc                                                                 


 Here is the input file:

 MD for Ethaline (1:2 ratio Chcl and EG)                                       
 &cntrl                                                                        
   ntx=5, irest=1,                   ! read coordinates                        
   ntb = 1,                          ! constant volume                         
   nstlim=10,                        ! MD for 10 steps                         
   ntpr=1, ntwr=10,                  ! output settings                         
   dt=0.0005,                        ! timestep in ps                          
   ig=71277,                         ! random seed                             
   cut = 10.0,                       ! non-bond cut off                        
   temp0 = 298.15, tempi = 298.15,   ! temerature settings                     
   ischeme = 1,                      !! Leapfrog middle scheme                 
   ithermostat = 1,                  !! Langevin thermostat                    
   therm_par = 100,                  !! thermostat parameter in ps^-1          
 /                                                                             
 &ewald                                                                        
   nfft1=60, nfft2=60, nfft3=60, netfrc=0, skinnb = 2.0                        
 /                                                                             


| irandom = 1, using AMBER's internal random number generator (default).
 
|--------------------- INFORMATION ----------------------
| GPU (CUDA) Version of PMEMD in use: NVIDIA GPU IN USE.
|                    Version 18.0.0
| 
|                      03/25/2018
| 
| Implementation by:
|                    Ross C. Walker     (SDSC)
|                    Scott Le Grand     (nVIDIA)
| 
| Version 18 performance extensions by:
|                    David Cerutti     (Rutgers)
| 
| Precision model in use:
|      [SPFP] - Single Precision Forces, 64-bit Fixed Point
|               Accumulation. (Default)
| 
|--------------------------------------------------------
 
|----------------- CITATION INFORMATION -----------------
|
|    When publishing work that utilized the CUDA version
|    of AMBER, please cite the following in addition to
|    the regular AMBER citations:
|
|  - Romelia Salomon-Ferrer; Andreas W. Goetz; Duncan
|    Poole; Scott Le Grand; Ross C. Walker "Routine
|    microsecond molecular dynamics simulations with
|    AMBER - Part II: Particle Mesh Ewald", J. Chem.
|    Theory Comput., 2013, 9 (9), pp3878-3888,
|    DOI: 10.1021/ct400314y.
|
|  - Andreas W. Goetz; Mark J. Williamson; Dong Xu;
|    Duncan Poole; Scott Le Grand; Ross C. Walker
|    "Routine microsecond molecular dynamics simulations
|    with AMBER - Part I: Generalized Born", J. Chem.
|    Theory Comput., 2012, 8 (5), pp1542-1555.
|
|  - Scott Le Grand; Andreas W. Goetz; Ross C. Walker
|    "SPFP: Speed without compromise - a mixed precision
|    model for GPU accelerated molecular dynamics
|    simulations.", Comp. Phys. Comm., 2013, 184
|    pp374-380, DOI: 10.1016/j.cpc.2012.09.022
|
|--------------------------------------------------------
 
|------------------- GPU DEVICE INFO --------------------
|
|            CUDA_VISIBLE_DEVICES: not set
|   CUDA Capable Devices Detected:      1
|           CUDA Device ID in use:      0
|                CUDA Device Name: NVIDIA GeForce GTX 1050 Ti
|     CUDA Device Global Mem Size:   4038 MB
| CUDA Device Num Multiprocessors:      6
|           CUDA Device Core Freq:   1.42 GHz
|
|--------------------------------------------------------
 
 
| Conditional Compilation Defines Used:
| PUBFFT
| BINTRAJ
| CUDA
| EMIL

| Largest sphere to fit in unit cell has radius =    29.170

| New format PARM file being parsed.
| Version =    1.000 Date = 05/16/21 Time = 18:16:03

| Note: 1-4 EEL scale factors are being read from the topology file.

| Note: 1-4 VDW scale factors are being read from the topology file.
| Duplicated    0 dihedrals

| Duplicated    0 dihedrals

--------------------------------------------------------------------------------
   1.  RESOURCE   USE: 
--------------------------------------------------------------------------------

 getting box info from netcdf restart file
 NATOM  =   21504 NTYPES =       7 NBONH =   13312 MBONA  =    6144
 NTHETH =   27136 MTHETA =    6144 NPHIH =   43520 MPHIA  =    4608
 NHPARM =       0 NPARM  =       0 NNB   =   96256 NRES   =    2048
 NBONA  =    6144 NTHETA =    6144 NPHIA =    4608 NUMBND =       6
 NUMANG =      10 NPTRA  =       8 NATYP =       7 NPHB   =       0
 IFBOX  =       1 NMXRS  =      21 IFCAP =       0 NEXTRA =       0
 NCOPY  =       0

| Coordinate Index Table dimensions:     9    9    9
| Direct force subcell size =     6.4823    6.4824    6.4850

     BOX TYPE: RECTILINEAR

--------------------------------------------------------------------------------
   2.  CONTROL  DATA  FOR  THE  RUN
--------------------------------------------------------------------------------

default_name                                                                    

General flags:
     imin    =       0, nmropt  =       0

Nature and format of input:
     ntx     =       5, irest   =       1, ntrx    =       1

Nature and format of output:
     ntxo    =       2, ntpr    =       1, ntrx    =       1, ntwr    =      10
     iwrap   =       0, ntwx    =       0, ntwv    =       0, ntwe    =       0
     ioutfm  =       1, ntwprt  =       0, idecomp =       0, rbornstat=      0

Potential function:
     ntf     =       1, ntb     =       1, igb     =       0, nsnb    =      25
     ipol    =       0, gbsa    =       0, iesp    =       0
     dielc   =   1.00000, cut     =  10.00000, intdiel =   1.00000

Frozen or restrained atoms:
     ibelly  =       0, ntr     =       0

Molecular dynamics:
     nstlim  =        10, nscm    =      1000, nrespa  =         1
     t       =   0.00000, dt      =   0.00050, vlimit  =  -1.00000

| Intermolecular bonds treatment:
|     no_intermolecular_bonds =       1

| Energy averages sample interval:
|     ene_avg_sampling =       1

Ewald parameters:
     verbose =       0, ew_type =       0, nbflag  =       1, use_pme =       1
     vdwmeth =       1, eedmeth =       1, netfrc  =       0
     Box X =   58.340   Box Y =   58.341   Box Z =   58.365
     Alpha =   90.000   Beta  =   90.000   Gamma =   90.000
     NFFT1 =   60       NFFT2 =   60       NFFT3 =   60
     Cutoff=   10.000   Tol   =0.100E-04
     Ewald Coefficient =  0.27511
     Interpolation order =    4

--------------------------------------------------------------------------------
   3.  ATOMIC COORDINATES AND VELOCITIES
--------------------------------------------------------------------------------

default_name                                                                    
 begin time read from input coords = 38000.000 ps

 
 Number of triangulated 3-point waters found:        0

     Sum of charges from parm topology file =   0.00153593
     Forcing neutrality...

| Dynamic Memory, Types Used:
| Reals              681172
| Integers          2055745

| Nonbonded Pairs Initial Allocation:     8425267

| GPU memory information (estimate):
| KB of GPU memory in use:    191575
| KB of CPU memory in use:     45751

--------------------------------------------------------------------------------
   4.  RESULTS
--------------------------------------------------------------------------------


 NSTEP =        1   TIME(PS) =   38000.000  TEMP(K) =   288.48  PRESS =     0.0
 Etot   =     -7051.4731  EKtot   =     18491.1543  EPtot      =    -25542.6274
 BOND   =      5989.3618  ANGLE   =      8603.1911  DIHED      =      3954.6724
 1-4 NB =       565.3615  1-4 EEL =    -31517.5611  VDWAALS    =    -11575.8080
 EELEC  =     -1561.8451  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        2   TIME(PS) =   38000.001  TEMP(K) =   286.89  PRESS =     0.0
 Etot   =     -7107.1862  EKtot   =     18389.2891  EPtot      =    -25496.4753
 BOND   =      6011.6966  ANGLE   =      8625.6181  DIHED      =      3954.1501
 1-4 NB =       565.6489  1-4 EEL =    -31517.7163  VDWAALS    =    -11573.7784
 EELEC  =     -1562.0942  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        3   TIME(PS) =   38000.001  TEMP(K) =   286.48  PRESS =     0.0
 Etot   =     -7105.2057  EKtot   =     18363.3066  EPtot      =    -25468.5123
 BOND   =      6019.1449  ANGLE   =      8643.4680  DIHED      =      3953.7223
 1-4 NB =       566.0671  1-4 EEL =    -31518.7673  VDWAALS    =    -11571.1784
 EELEC  =     -1560.9690  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        4   TIME(PS) =   38000.002  TEMP(K) =   286.59  PRESS =     0.0
 Etot   =     -7093.7605  EKtot   =     18370.1816  EPtot      =    -25463.9421
 BOND   =      6009.1779  ANGLE   =      8654.5000  DIHED      =      3953.6560
 1-4 NB =       566.6098  1-4 EEL =    -31520.1291  VDWAALS    =    -11568.2351
 EELEC  =     -1559.5215  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        5   TIME(PS) =   38000.002  TEMP(K) =   286.44  PRESS =     0.0
 Etot   =     -7107.5750  EKtot   =     18360.2930  EPtot      =    -25467.8680
 BOND   =      5998.3595  ANGLE   =      8659.5569  DIHED      =      3953.8408
 1-4 NB =       567.2082  1-4 EEL =    -31521.9193  VDWAALS    =    -11565.7479
 EELEC  =     -1559.1663  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        6   TIME(PS) =   38000.003  TEMP(K) =   286.66  PRESS =     0.0
 Etot   =     -7100.7909  EKtot   =     18374.4160  EPtot      =    -25475.2069
 BOND   =      5991.3395  ANGLE   =      8660.0689  DIHED      =      3954.0974
 1-4 NB =       567.9024  1-4 EEL =    -31524.1315  VDWAALS    =    -11564.7103
 EELEC  =     -1559.7733  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        7   TIME(PS) =   38000.003  TEMP(K) =   286.51  PRESS =     0.0
 Etot   =     -7120.3306  EKtot   =     18364.8398  EPtot      =    -25485.1705
 BOND   =      5987.4114  ANGLE   =      8655.5912  DIHED      =      3954.4874
 1-4 NB =       568.6711  1-4 EEL =    -31526.3100  VDWAALS    =    -11564.9006
 EELEC  =     -1560.1210  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        8   TIME(PS) =   38000.004  TEMP(K) =   285.97  PRESS =     0.0
 Etot   =     -7166.9443  EKtot   =     18330.6387  EPtot      =    -25497.5830
 BOND   =      5987.0327  ANGLE   =      8644.6836  DIHED      =      3954.7501
 1-4 NB =       569.3997  1-4 EEL =    -31527.9773  VDWAALS    =    -11565.2807
 EELEC  =     -1560.1911  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        9   TIME(PS) =   38000.004  TEMP(K) =   286.86  PRESS =     0.0
 Etot   =     -7124.4868  EKtot   =     18387.1562  EPtot      =    -25511.6430
 BOND   =      5987.9151  ANGLE   =      8630.2284  DIHED      =      3954.8050
 1-4 NB =       570.0090  1-4 EEL =    -31528.6795  VDWAALS    =    -11565.9690
 EELEC  =     -1559.9520  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =       10   TIME(PS) =   38000.005  TEMP(K) =   287.88  PRESS =     0.0
 Etot   =     -7075.4840  EKtot   =     18452.5566  EPtot      =    -25528.0407
 BOND   =      5986.2526  ANGLE   =      8615.6231  DIHED      =      3954.9482
 1-4 NB =       570.5803  1-4 EEL =    -31528.8794  VDWAALS    =    -11566.9368
 EELEC  =     -1559.6287  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


      A V E R A G E S   O V E R      10 S T E P S


 NSTEP =       10   TIME(PS) =   38000.005  TEMP(K) =   286.88  PRESS =     0.0
 Etot   =     -7105.3237  EKtot   =     18388.3832  EPtot      =    -25493.7069
 BOND   =      5996.7692  ANGLE   =      8639.2529  DIHED      =      3954.3130
 1-4 NB =       567.7458  1-4 EEL =    -31523.2071  VDWAALS    =    -11568.2545
 EELEC  =     -1560.3262  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


      R M S  F L U C T U A T I O N S


 NSTEP =       10   TIME(PS) =   38000.005  TEMP(K) =     0.71  PRESS =     0.0
 Etot   =        28.8836  EKtot   =        45.2940  EPtot      =        25.5341
 BOND   =        11.5573  ANGLE   =        18.7574  DIHED      =         0.4553
 1-4 NB =         1.7678  1-4 EEL =         4.3385  VDWAALS    =         3.7685
 EELEC  =         0.9398  EHBOND  =         0.0000  RESTRAINT  =         0.0000
|E(PBS) =        27.8862
 ------------------------------------------------------------------------------

--------------------------------------------------------------------------------
   5.  TIMINGS
--------------------------------------------------------------------------------

|  NonSetup CPU Time in Major Routines:
|
|     Routine           Sec        %
|     ------------------------------
|     Nonbond           0.03   68.03
|     Bond              0.00    0.00
|     Angle             0.00    0.00
|     Dihedral          0.00    0.00
|     Shake             0.00    0.00
|     RunMD             0.01   31.87
|     Other             0.00    0.11
|     ------------------------------
|     Total             0.04

|  PME Nonbond Pairlist CPU Time:
|
|     Routine              Sec        %
|     ---------------------------------
|     Set Up Cit           0.00    0.00
|     Build List           0.00    0.00
|     ---------------------------------
|     Total                0.00    0.00

|  PME Direct Force CPU Time:
|
|     Routine              Sec        %
|     ---------------------------------
|     NonBonded Calc       0.00    0.00
|     Exclude Masked       0.00    0.00
|     Other                0.00    0.02
|     ---------------------------------
|     Total                0.00    0.02

|  PME Reciprocal Force CPU Time:
|
|     Routine              Sec        %
|     ---------------------------------
|     1D bspline           0.00    0.00
|     Grid Charges         0.00    0.00
|     Scalar Sum           0.00    0.00
|     Gradient Sum         0.00    0.00
|     FFT                  0.00    0.00
|     ---------------------------------
|     Total                0.00    0.00

|  Final Performance Info:
|     -----------------------------------------------------
|     Average timings for last       9 steps:
|     Elapsed(s) =       0.03 Per Step(ms) =       3.14
|         ns/day =      13.76   seconds/ns =    6281.30
|
|     Average timings for all steps:
|     Elapsed(s) =       0.04 Per Step(ms) =       4.35
|         ns/day =       9.93   seconds/ns =    8705.18
|     -----------------------------------------------------

|  Setup CPU time:            0.44 seconds
|  NonSetup CPU time:         0.04 seconds
|  Total CPU time:            0.48 seconds     0.00 hours

|  Setup wall time:           0    seconds
|  NonSetup wall time:        0    seconds
|  Total wall time:           0    seconds     0.00 hours
