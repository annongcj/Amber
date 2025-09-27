
          -------------------------------------------------------
          Amber 22 PMEMD                              2022
          -------------------------------------------------------

| PMEMD implementation of SANDER, Release 22

|  Compiled date/time: Tue Apr 26 08:53:16 2022
| Run on 04/27/2022 at 10:47:06

|   Executable path: /home/mahdieh/apps/amber_gitlab2/amber/dist/bin/pmemd.cuda_DPFP
| Working directory: /home/mahdieh/apps/amber_gitlab2/amber/test/middle-scheme/ETH
|          Hostname: TheraC02

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
|      [DPFP] - Double Precision Forces, 64-bit Fixed point
|               Accumulation.
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
|--------------------------------------------------------
 
|------------------- GPU DEVICE INFO --------------------
|
|            CUDA_VISIBLE_DEVICES: not set
|   CUDA Capable Devices Detected:      8
|           CUDA Device ID in use:      0
|                CUDA Device Name: 
|     CUDA Device Global Mem Size:  32752 MB
| CUDA Device Num Multiprocessors:    120
|           CUDA Device Core Freq:   1.50 GHz
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
| KB of GPU memory in use:    218417
| KB of CPU memory in use:     57137

--------------------------------------------------------------------------------
   4.  RESULTS
--------------------------------------------------------------------------------


 NSTEP =        1   TIME(PS) =   38000.000  TEMP(K) =   289.56  PRESS =     0.0
 Etot   =     -6982.2748  EKtot   =     18560.3555  EPtot      =    -25542.6303
 BOND   =      5989.3618  ANGLE   =      8603.1911  DIHED      =      3954.6724
 1-4 NB =       565.3616  1-4 EEL =    -31517.5613  VDWAALS    =    -11575.8095
 EELEC  =     -1561.8464  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        2   TIME(PS) =   38000.001  TEMP(K) =   288.78  PRESS =     0.0
 Etot   =     -6980.6083  EKtot   =     18510.3187  EPtot      =    -25490.9270
 BOND   =      6020.1872  ANGLE   =      8623.1748  DIHED      =      3954.1303
 1-4 NB =       565.6934  1-4 EEL =    -31517.3866  VDWAALS    =    -11573.7741
 EELEC  =     -1562.9519  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        3   TIME(PS) =   38000.001  TEMP(K) =   288.03  PRESS =     0.0
 Etot   =     -6984.5421  EKtot   =     18462.3613  EPtot      =    -25446.9034
 BOND   =      6050.2454  ANGLE   =      8636.4442  DIHED      =      3953.6816
 1-4 NB =       566.0289  1-4 EEL =    -31517.6836  VDWAALS    =    -11571.4566
 EELEC  =     -1564.1633  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        4   TIME(PS) =   38000.002  TEMP(K) =   287.82  PRESS =     0.0
 Etot   =     -6968.3564  EKtot   =     18448.8875  EPtot      =    -25417.2439
 BOND   =      6071.0722  ANGLE   =      8643.9696  DIHED      =      3953.5387
 1-4 NB =       566.3041  1-4 EEL =    -31518.2342  VDWAALS    =    -11568.8029
 EELEC  =     -1565.0913  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        5   TIME(PS) =   38000.002  TEMP(K) =   287.26  PRESS =     0.0
 Etot   =     -6981.4955  EKtot   =     18412.7511  EPtot      =    -25394.2465
 BOND   =      6084.6352  ANGLE   =      8650.4978  DIHED      =      3953.7618
 1-4 NB =       566.7108  1-4 EEL =    -31519.1029  VDWAALS    =    -11565.7003
 EELEC  =     -1565.0488  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        6   TIME(PS) =   38000.003  TEMP(K) =   286.97  PRESS =     0.0
 Etot   =     -6992.3755  EKtot   =     18394.6373  EPtot      =    -25387.0128
 BOND   =      6084.6564  ANGLE   =      8654.9347  DIHED      =      3954.0925
 1-4 NB =       567.2124  1-4 EEL =    -31520.6335  VDWAALS    =    -11563.0715
 EELEC  =     -1564.2038  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        7   TIME(PS) =   38000.003  TEMP(K) =   286.99  PRESS =     0.0
 Etot   =     -7003.5186  EKtot   =     18395.5291  EPtot      =    -25399.0477
 BOND   =      6070.3501  ANGLE   =      8654.7556  DIHED      =      3954.2599
 1-4 NB =       567.7274  1-4 EEL =    -31522.6126  VDWAALS    =    -11560.7961
 EELEC  =     -1562.7319  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        8   TIME(PS) =   38000.004  TEMP(K) =   287.04  PRESS =     0.0
 Etot   =     -7032.4589  EKtot   =     18399.1947  EPtot      =    -25431.6536
 BOND   =      6040.9694  ANGLE   =      8649.5812  DIHED      =      3954.1626
 1-4 NB =       568.1704  1-4 EEL =    -31524.9070  VDWAALS    =    -11558.6438
 EELEC  =     -1560.9864  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        9   TIME(PS) =   38000.004  TEMP(K) =   288.31  PRESS =     0.0
 Etot   =     -6997.9300  EKtot   =     18480.0574  EPtot      =    -25477.9874
 BOND   =      6001.0213  ANGLE   =      8642.1665  DIHED      =      3953.9263
 1-4 NB =       568.5502  1-4 EEL =    -31526.5368  VDWAALS    =    -11557.2309
 EELEC  =     -1559.8840  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =       10   TIME(PS) =   38000.005  TEMP(K) =   289.08  PRESS =     0.0
 Etot   =     -6989.4297  EKtot   =     18529.4547  EPtot      =    -25518.8843
 BOND   =      5965.4146  ANGLE   =      8635.4920  DIHED      =      3953.7582
 1-4 NB =       569.0547  1-4 EEL =    -31527.7175  VDWAALS    =    -11556.4629
 EELEC  =     -1558.4235  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


      A V E R A G E S   O V E R      10 S T E P S


 NSTEP =       10   TIME(PS) =   38000.005  TEMP(K) =   287.98  PRESS =     0.0
 Etot   =     -6991.2990  EKtot   =     18459.3547  EPtot      =    -25450.6537
 BOND   =      6037.7913  ANGLE   =      8639.4207  DIHED      =      3953.9984
 1-4 NB =       567.0814  1-4 EEL =    -31521.2376  VDWAALS    =    -11565.1749
 EELEC  =     -1562.5331  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


      R M S  F L U C T U A T I O N S


 NSTEP =       10   TIME(PS) =   38000.005  TEMP(K) =     0.89  PRESS =     0.0
 Etot   =        16.6010  EKtot   =        56.8444  EPtot      =        51.8368
 BOND   =        39.9616  ANGLE   =        15.2670  DIHED      =         0.3172
 1-4 NB =         1.2002  1-4 EEL =         3.7462  VDWAALS    =         6.6723
 EELEC  =         2.1286  EHBOND  =         0.0000  RESTRAINT  =         0.0000
|E(PBS) =        16.8315
 ------------------------------------------------------------------------------

--------------------------------------------------------------------------------
   5.  TIMINGS
--------------------------------------------------------------------------------

|  NonSetup CPU Time in Major Routines:
|
|     Routine           Sec        %
|     ------------------------------
|     Nonbond           0.01   67.46
|     Bond              0.00    0.00
|     Angle             0.00    0.00
|     Dihedral          0.00    0.00
|     Shake             0.00    0.00
|     RunMD             0.01   32.31
|     Other             0.00    0.23
|     ------------------------------
|     Total             0.02

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
|     Other                0.00    0.04
|     ---------------------------------
|     Total                0.00    0.04

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
|     Elapsed(s) =       0.02 Per Step(ms) =       2.70
|         ns/day =      16.00   seconds/ns =    5400.66
|
|     Average timings for all steps:
|     Elapsed(s) =       0.03 Per Step(ms) =       2.97
|         ns/day =      14.54   seconds/ns =    5943.80
|     -----------------------------------------------------

|  Setup CPU time:            3.05 seconds
|  NonSetup CPU time:         0.02 seconds
|  Total CPU time:            3.07 seconds     0.00 hours

|  Setup wall time:           3    seconds
|  NonSetup wall time:        0    seconds
|  Total wall time:           3    seconds     0.00 hours
