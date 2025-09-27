
          -------------------------------------------------------
          Amber 22 PMEMD                              2022
          -------------------------------------------------------

| PMEMD implementation of SANDER, Release 22

|  Compiled date/time: Tue Apr 26 08:53:16 2022
| Run on 04/27/2022 at 11:17:31

|   Executable path: /home/mahdieh/apps/amber_gitlab2/amber/dist/bin/pmemd.cuda_DPFP
| Working directory: /home/mahdieh/apps/amber_gitlab2/amber/test/middle-scheme/DNA7
|          Hostname: TheraC02

  [-O]verwriting output

File Assignments:
|   MDIN: mdin                                                                  
|  MDOUT: mdout.DNA_middle_scheme                                               
| INPCRD: com-sol_equil_0.rst                                                   
|   PARM: prmtop                                                                
| RESTRT: restrt                                                                
|   REFC: refc                                                                  
|  MDVEL: mdvel                                                                 
|   MDEN: mden                                                                  
|  MDCRD: mdcrd                                                                 
| MDINFO: mdinfo                                                                
|  MDFRC: mdfrc                                                                 


 Here is the input file:

 MD for DNA                                                                    
 &cntrl                                                                        
   ntx=5, irest=1,                   ! read coordinates                        
   ntc=2, ntf=2, tol=0.0000001,      ! constrain lengths of bonds              
                                                                               
   ntb = 1,                          ! constant volume                         
   nstlim=10,                        ! MD for 10 steps                         
   ntpr=1, ntwr=10000,               ! output settings                         
   dt=0.001,                         ! timestep in ps                          
   ig=71277,                         ! random seed                             
   cut = 9.0,                        ! non-bond cut off                        
   temp0 = 300, tempi = 300,         ! temerature settings                     
   ischeme = 1,                      !! Leapfrog middle scheme                 
   ithermostat = 1,                  !! Langevin thermostat                    
   therm_par = 5.0,                  !! thermostat parameter in ps^-1          
 /                                                                             
 &ewald                                                                        
   nfft1=60, nfft2=60, nfft3=60, netfrc=0,                                     
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

| Largest sphere to fit in unit cell has radius =    21.045

| New format PARM file being parsed.
| Version =    1.000 Date = 02/20/19 Time = 20:21:43

| Note: 1-4 EEL scale factors are being read from the topology file.

| Note: 1-4 VDW scale factors are being read from the topology file.
| Duplicated    0 dihedrals

| Duplicated    0 dihedrals

--------------------------------------------------------------------------------
   1.  RESOURCE   USE: 
--------------------------------------------------------------------------------

 getting box info from netcdf restart file
 NATOM  =   10511 NTYPES =      17 NBONH =   10218 MBONA  =     314
 NTHETH =     389 MTHETA =     481 NPHIH =     763 MPHIA  =     918
 NHPARM =       0 NPARM  =       0 NNB   =   15889 NRES   =    3377
 NBONA  =     314 NTHETA =     481 NPHIA =     918 NUMBND =      41
 NUMANG =      75 NPTRA  =      56 NATYP =      29 NPHB   =       1
 IFBOX  =       2 NMXRS  =      33 IFCAP =       0 NEXTRA =       0
 NCOPY  =       0

| Coordinate Index Table dimensions:     9    9    9
| Direct force subcell size =     5.7277    5.7277    5.7277

     BOX TYPE: TRUNCATED OCTAHEDRON

--------------------------------------------------------------------------------
   2.  CONTROL  DATA  FOR  THE  RUN
--------------------------------------------------------------------------------

default_name                                                                    

General flags:
     imin    =       0, nmropt  =       0

Nature and format of input:
     ntx     =       5, irest   =       1, ntrx    =       1

Nature and format of output:
     ntxo    =       2, ntpr    =       1, ntrx    =       1, ntwr    =   10000
     iwrap   =       0, ntwx    =       0, ntwv    =       0, ntwe    =       0
     ioutfm  =       1, ntwprt  =       0, idecomp =       0, rbornstat=      0

Potential function:
     ntf     =       2, ntb     =       1, igb     =       0, nsnb    =      25
     ipol    =       0, gbsa    =       0, iesp    =       0
     dielc   =   1.00000, cut     =   9.00000, intdiel =   1.00000

Frozen or restrained atoms:
     ibelly  =       0, ntr     =       0

Molecular dynamics:
     nstlim  =        10, nscm    =      1000, nrespa  =         1
     t       =   0.00000, dt      =   0.00100, vlimit  =  -1.00000

SHAKE:
     ntc     =       2, jfastw  =       0
     tol     =   0.00000

| Intermolecular bonds treatment:
|     no_intermolecular_bonds =       1

| Energy averages sample interval:
|     ene_avg_sampling =       1

Ewald parameters:
     verbose =       0, ew_type =       0, nbflag  =       1, use_pme =       1
     vdwmeth =       1, eedmeth =       1, netfrc  =       0
     Box X =   51.549   Box Y =   51.549   Box Z =   51.549
     Alpha =  109.471   Beta  =  109.471   Gamma =  109.471
     NFFT1 =   60       NFFT2 =   60       NFFT3 =   60
     Cutoff=    9.000   Tol   =0.100E-04
     Ewald Coefficient =  0.30768
     Interpolation order =    4

--------------------------------------------------------------------------------
   3.  ATOMIC COORDINATES AND VELOCITIES
--------------------------------------------------------------------------------

default_name                                                                    
 begin time read from input coords =    50.000 ps

 
 Number of triangulated 3-point waters found:     3351

     Sum of charges from parm topology file =   0.00000006
     Forcing neutrality...

| Dynamic Memory, Types Used:
| Reals              396757
| Integers           330055

| Nonbonded Pairs Initial Allocation:     2396508

| GPU memory information (estimate):
| KB of GPU memory in use:     89290
| KB of CPU memory in use:     33418

--------------------------------------------------------------------------------
   4.  RESULTS
--------------------------------------------------------------------------------


 NSTEP =        1   TIME(PS) =      50.001  TEMP(K) =   295.48  PRESS =     0.0
 Etot   =    -29236.0128  EKtot   =      6257.9186  EPtot      =    -35493.9314
 BOND   =       109.0883  ANGLE   =       238.3825  DIHED      =       351.3869
 1-4 NB =       117.1807  1-4 EEL =      -239.3030  VDWAALS    =      4481.9885
 EELEC  =    -40552.6553  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        2   TIME(PS) =      50.002  TEMP(K) =   296.14  PRESS =     0.0
 Etot   =    -29232.8601  EKtot   =      6271.8622  EPtot      =    -35504.7223
 BOND   =       109.7376  ANGLE   =       235.6698  DIHED      =       351.6545
 1-4 NB =       116.4693  1-4 EEL =      -238.5945  VDWAALS    =      4476.6359
 EELEC  =    -40556.2949  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        3   TIME(PS) =      50.003  TEMP(K) =   296.82  PRESS =     0.0
 Etot   =    -29229.0650  EKtot   =      6286.1922  EPtot      =    -35515.2572
 BOND   =       109.7991  ANGLE   =       232.0853  DIHED      =       352.0084
 1-4 NB =       115.6901  1-4 EEL =      -238.0099  VDWAALS    =      4471.5057
 EELEC  =    -40558.3358  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        4   TIME(PS) =      50.004  TEMP(K) =   297.05  PRESS =     0.0
 Etot   =    -29233.8324  EKtot   =      6290.9713  EPtot      =    -35524.8037
 BOND   =       109.3022  ANGLE   =       228.4242  DIHED      =       352.3890
 1-4 NB =       114.9536  1-4 EEL =      -237.6826  VDWAALS    =      4467.0322
 EELEC  =    -40559.2222  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        5   TIME(PS) =      50.005  TEMP(K) =   297.26  PRESS =     0.0
 Etot   =    -29237.4324  EKtot   =      6295.4341  EPtot      =    -35532.8665
 BOND   =       108.4580  ANGLE   =       225.0978  DIHED      =       352.6826
 1-4 NB =       114.3356  1-4 EEL =      -237.6423  VDWAALS    =      4463.1638
 EELEC  =    -40558.9621  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        6   TIME(PS) =      50.006  TEMP(K) =   297.54  PRESS =     0.0
 Etot   =    -29237.4012  EKtot   =      6301.4884  EPtot      =    -35538.8897
 BOND   =       108.0691  ANGLE   =       222.2229  DIHED      =       352.7803
 1-4 NB =       113.8519  1-4 EEL =      -237.8633  VDWAALS    =      4459.1169
 EELEC  =    -40557.0675  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        7   TIME(PS) =      50.007  TEMP(K) =   297.54  PRESS =     0.0
 Etot   =    -29240.1171  EKtot   =      6301.4892  EPtot      =    -35541.6063
 BOND   =       108.6527  ANGLE   =       220.0932  DIHED      =       352.8122
 1-4 NB =       113.4707  1-4 EEL =      -238.2578  VDWAALS    =      4454.9477
 EELEC  =    -40553.3250  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        8   TIME(PS) =      50.008  TEMP(K) =   297.65  PRESS =     0.0
 Etot   =    -29236.4785  EKtot   =      6303.8615  EPtot      =    -35540.3400
 BOND   =       110.0470  ANGLE   =       219.2253  DIHED      =       352.8078
 1-4 NB =       113.2570  1-4 EEL =      -238.9358  VDWAALS    =      4450.7395
 EELEC  =    -40547.4807  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        9   TIME(PS) =      50.009  TEMP(K) =   297.45  PRESS =     0.0
 Etot   =    -29237.1992  EKtot   =      6299.6441  EPtot      =    -35536.8433
 BOND   =       111.8173  ANGLE   =       219.5155  DIHED      =       352.6480
 1-4 NB =       113.2263  1-4 EEL =      -239.9006  VDWAALS    =      4446.5605
 EELEC  =    -40540.7103  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =       10   TIME(PS) =      50.010  TEMP(K) =   297.54  PRESS =     0.0
 Etot   =    -29230.8497  EKtot   =      6301.4478  EPtot      =    -35532.2975
 BOND   =       114.0114  ANGLE   =       220.5322  DIHED      =       352.3399
 1-4 NB =       113.3393  1-4 EEL =      -240.9916  VDWAALS    =      4442.4863
 EELEC  =    -40534.0150  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


      A V E R A G E S   O V E R      10 S T E P S


 NSTEP =       10   TIME(PS) =      50.010  TEMP(K) =   297.05  PRESS =     0.0
 Etot   =    -29235.1249  EKtot   =      6291.0309  EPtot      =    -35526.1558
 BOND   =       109.8983  ANGLE   =       226.1249  DIHED      =       352.3509
 1-4 NB =       114.5774  1-4 EEL =      -238.7181  VDWAALS    =      4461.4177
 EELEC  =    -40551.8069  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


      R M S  F L U C T U A T I O N S


 NSTEP =       10   TIME(PS) =      50.010  TEMP(K) =     0.68  PRESS =     0.0
 Etot   =         3.2272  EKtot   =        14.3970  EPtot      =        15.5482
 BOND   =         1.6924  ANGLE   =         6.7631  DIHED      =         0.4835
 1-4 NB =         1.3656  1-4 EEL =         1.0318  VDWAALS    =        12.3750
 EELEC  =         8.1120  EHBOND  =         0.0000  RESTRAINT  =         0.0000
|E(PBS) =         3.2832
 ------------------------------------------------------------------------------

--------------------------------------------------------------------------------
   5.  TIMINGS
--------------------------------------------------------------------------------

|  NonSetup CPU Time in Major Routines:
|
|     Routine           Sec        %
|     ------------------------------
|     Nonbond           0.01   56.28
|     Bond              0.00    0.00
|     Angle             0.00    0.00
|     Dihedral          0.00    0.00
|     Shake             0.00    0.47
|     RunMD             0.00   43.00
|     Other             0.00    0.25
|     ------------------------------
|     Total             0.01

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
|     Other                0.00    0.06
|     ---------------------------------
|     Total                0.00    0.06

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
|     Elapsed(s) =       0.01 Per Step(ms) =       1.55
|         ns/day =      55.91   seconds/ns =    1545.22
|
|     Average timings for all steps:
|     Elapsed(s) =       0.02 Per Step(ms) =       1.80
|         ns/day =      47.91   seconds/ns =    1803.30
|     -----------------------------------------------------

|  Setup CPU time:            2.93 seconds
|  NonSetup CPU time:         0.01 seconds
|  Total CPU time:            2.94 seconds     0.00 hours

|  Setup wall time:           3    seconds
|  NonSetup wall time:        0    seconds
|  Total wall time:           3    seconds     0.00 hours
