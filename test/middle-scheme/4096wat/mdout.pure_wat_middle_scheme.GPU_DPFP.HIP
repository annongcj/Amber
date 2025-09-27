
          -------------------------------------------------------
          Amber 22 PMEMD                              2022
          -------------------------------------------------------

| PMEMD implementation of SANDER, Release 22

|  Compiled date/time: Tue Apr 26 08:53:16 2022
| Run on 04/27/2022 at 11:17:27

|   Executable path: /home/mahdieh/apps/amber_gitlab2/amber/dist/bin/pmemd.cuda_DPFP
| Working directory: /home/mahdieh/apps/amber_gitlab2/amber/test/middle-scheme/4096wat
|          Hostname: TheraC02

  [-O]verwriting output

File Assignments:
|   MDIN: mdin                                                                  
|  MDOUT: mdout.pure_wat_middle_scheme                                          
| INPCRD: eq1.x                                                                 
|   PARM: prmtop                                                                
| RESTRT: restrt                                                                
|   REFC: refc                                                                  
|  MDVEL: mdvel                                                                 
|   MDEN: mden                                                                  
|  MDCRD: mdcrd                                                                 
| MDINFO: mdinfo                                                                
|  MDFRC: mdfrc                                                                 


 Here is the input file:

 short md, nve ensemble                                                        
 &cntrl                                                                        
   ntx=5, irest=1,                   ! read coordinates                        
   ntc=2, ntf=2, tol=0.0000001,      ! constrain lengths of bonds              
                                     ! having hydrogen atoms                   
   nstlim=10,                        ! MD for 10 steps                         
   ntpr=1, ntwr=10000,               ! output settings                         
   dt=0.001,                         ! timestep in ps                          
   ig=71277,                         ! random seed                             
   cut = 7.0,                        ! non-bond cut off                        
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

| Largest sphere to fit in unit cell has radius =    24.800

|  INFO: Old style PARM file read


| Note: 1-4 EEL scale factors were NOT found in the topology file.
|       Using default value of 1.2.

| Note: 1-4 VDW scale factors were NOT found in the topology file.
|       Using default value of 2.0.
| Duplicated    0 dihedrals

| Duplicated    0 dihedrals

--------------------------------------------------------------------------------
   1.  RESOURCE   USE: 
--------------------------------------------------------------------------------

 getting new box info from bottom of inpcrd
 NATOM  =   12288 NTYPES =       2 NBONH =   12288 MBONA  =       0
 NTHETH =       0 MTHETA =       0 NPHIH =       0 MPHIA  =       0
 NHPARM =       0 NPARM  =       0 NNB   =   16384 NRES   =    4096
 NBONA  =       0 NTHETA =       0 NPHIA =       0 NUMBND =       2
 NUMANG =       0 NPTRA  =       0 NATYP =       2 NPHB   =       1
 IFBOX  =       1 NMXRS  =       3 IFCAP =       0 NEXTRA =       0
 NCOPY  =       0

| Coordinate Index Table dimensions:    12   12   12
| Direct force subcell size =     4.1333    4.1333    4.1333

     BOX TYPE: RECTILINEAR

--------------------------------------------------------------------------------
   2.  CONTROL  DATA  FOR  THE  RUN
--------------------------------------------------------------------------------

                                                                                

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
     dielc   =   1.00000, cut     =   7.00000, intdiel =   1.00000

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
     Box X =   49.600   Box Y =   49.600   Box Z =   49.600
     Alpha =   90.000   Beta  =   90.000   Gamma =   90.000
     NFFT1 =   60       NFFT2 =   60       NFFT3 =   60
     Cutoff=    7.000   Tol   =0.100E-04
     Ewald Coefficient =  0.40167
     Interpolation order =    4

--------------------------------------------------------------------------------
   3.  ATOMIC COORDINATES AND VELOCITIES
--------------------------------------------------------------------------------

                                                                                
 begin time read from input coords =     1.000 ps

 
 Number of triangulated 3-point waters found:     4096

     Sum of charges from parm topology file =   0.00000000
     Forcing neutrality...

| Dynamic Memory, Types Used:
| Reals              421811
| Integers           331783

| Nonbonded Pairs Initial Allocation:     1452441

| GPU memory information (estimate):
| KB of GPU memory in use:     73610
| KB of CPU memory in use:     27674

--------------------------------------------------------------------------------
   4.  RESULTS
--------------------------------------------------------------------------------


 NSTEP =        1   TIME(PS) =       1.001  TEMP(K) =   297.39  PRESS =     0.0
 Etot   =    -32081.2193  EKtot   =      7261.7995  EPtot      =    -39343.0189
 BOND   =         0.0000  ANGLE   =         0.0000  DIHED      =         0.0000
 1-4 NB =         0.0000  1-4 EEL =         0.0000  VDWAALS    =      6029.0685
 EELEC  =    -45372.0874  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        2   TIME(PS) =       1.002  TEMP(K) =   296.97  PRESS =     0.0
 Etot   =    -32083.6524  EKtot   =      7251.4959  EPtot      =    -39335.1484
 BOND   =         0.0000  ANGLE   =         0.0000  DIHED      =         0.0000
 1-4 NB =         0.0000  1-4 EEL =         0.0000  VDWAALS    =      6030.1349
 EELEC  =    -45365.2833  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        3   TIME(PS) =       1.003  TEMP(K) =   296.92  PRESS =     0.0
 Etot   =    -32079.5168  EKtot   =      7250.4137  EPtot      =    -39329.9305
 BOND   =         0.0000  ANGLE   =         0.0000  DIHED      =         0.0000
 1-4 NB =         0.0000  1-4 EEL =         0.0000  VDWAALS    =      6031.1486
 EELEC  =    -45361.0792  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        4   TIME(PS) =       1.004  TEMP(K) =   296.77  PRESS =     0.0
 Etot   =    -32080.4156  EKtot   =      7246.7530  EPtot      =    -39327.1686
 BOND   =         0.0000  ANGLE   =         0.0000  DIHED      =         0.0000
 1-4 NB =         0.0000  1-4 EEL =         0.0000  VDWAALS    =      6032.2604
 EELEC  =    -45359.4290  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        5   TIME(PS) =       1.005  TEMP(K) =   296.23  PRESS =     0.0
 Etot   =    -32092.6656  EKtot   =      7233.4374  EPtot      =    -39326.1030
 BOND   =         0.0000  ANGLE   =         0.0000  DIHED      =         0.0000
 1-4 NB =         0.0000  1-4 EEL =         0.0000  VDWAALS    =      6033.1115
 EELEC  =    -45359.2145  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        6   TIME(PS) =       1.006  TEMP(K) =   295.98  PRESS =     0.0
 Etot   =    -32098.8509  EKtot   =      7227.3261  EPtot      =    -39326.1770
 BOND   =         0.0000  ANGLE   =         0.0000  DIHED      =         0.0000
 1-4 NB =         0.0000  1-4 EEL =         0.0000  VDWAALS    =      6033.5964
 EELEC  =    -45359.7735  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        7   TIME(PS) =       1.007  TEMP(K) =   295.97  PRESS =     0.0
 Etot   =    -32099.8123  EKtot   =      7227.2939  EPtot      =    -39327.1062
 BOND   =         0.0000  ANGLE   =         0.0000  DIHED      =         0.0000
 1-4 NB =         0.0000  1-4 EEL =         0.0000  VDWAALS    =      6033.7271
 EELEC  =    -45360.8333  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        8   TIME(PS) =       1.008  TEMP(K) =   296.38  PRESS =     0.0
 Etot   =    -32091.2513  EKtot   =      7237.1354  EPtot      =    -39328.3867
 BOND   =         0.0000  ANGLE   =         0.0000  DIHED      =         0.0000
 1-4 NB =         0.0000  1-4 EEL =         0.0000  VDWAALS    =      6033.7138
 EELEC  =    -45362.1005  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        9   TIME(PS) =       1.009  TEMP(K) =   296.81  PRESS =     0.0
 Etot   =    -32081.3289  EKtot   =      7247.7828  EPtot      =    -39329.1117
 BOND   =         0.0000  ANGLE   =         0.0000  DIHED      =         0.0000
 1-4 NB =         0.0000  1-4 EEL =         0.0000  VDWAALS    =      6033.8319
 EELEC  =    -45362.9436  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =       10   TIME(PS) =       1.010  TEMP(K) =   296.60  PRESS =     0.0
 Etot   =    -32087.1620  EKtot   =      7242.5489  EPtot      =    -39329.7109
 BOND   =         0.0000  ANGLE   =         0.0000  DIHED      =         0.0000
 1-4 NB =         0.0000  1-4 EEL =         0.0000  VDWAALS    =      6033.7782
 EELEC  =    -45363.4891  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


      A V E R A G E S   O V E R      10 S T E P S


 NSTEP =       10   TIME(PS) =       1.010  TEMP(K) =   296.60  PRESS =     0.0
 Etot   =    -32087.5875  EKtot   =      7242.5987  EPtot      =    -39330.1862
 BOND   =         0.0000  ANGLE   =         0.0000  DIHED      =         0.0000
 1-4 NB =         0.0000  1-4 EEL =         0.0000  VDWAALS    =      6032.4371
 EELEC  =    -45362.6233  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


      R M S  F L U C T U A T I O N S


 NSTEP =       10   TIME(PS) =       1.010  TEMP(K) =     0.44  PRESS =     0.0
 Etot   =         7.2628  EKtot   =        10.6504  EPtot      =         4.9530
 BOND   =         0.0000  ANGLE   =         0.0000  DIHED      =         0.0000
 1-4 NB =         0.0000  1-4 EEL =         0.0000  VDWAALS    =         1.6485
 EELEC  =         3.6511  EHBOND  =         0.0000  RESTRAINT  =         0.0000
|E(PBS) =         7.0376
 ------------------------------------------------------------------------------

--------------------------------------------------------------------------------
   5.  TIMINGS
--------------------------------------------------------------------------------

|  NonSetup CPU Time in Major Routines:
|
|     Routine           Sec        %
|     ------------------------------
|     Nonbond           0.00   49.03
|     Bond              0.00    0.00
|     Angle             0.00    0.00
|     Dihedral          0.00    0.00
|     Shake             0.00    0.53
|     RunMD             0.00   50.12
|     Other             0.00    0.31
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
|     Other                0.00    0.08
|     ---------------------------------
|     Total                0.00    0.08

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
|     Elapsed(s) =       0.01 Per Step(ms) =       1.52
|         ns/day =      56.85   seconds/ns =    1519.88
|
|     Average timings for all steps:
|     Elapsed(s) =       0.02 Per Step(ms) =       1.78
|         ns/day =      48.51   seconds/ns =    1781.10
|     -----------------------------------------------------

|  Setup CPU time:            2.92 seconds
|  NonSetup CPU time:         0.01 seconds
|  Total CPU time:            2.93 seconds     0.00 hours

|  Setup wall time:           3    seconds
|  NonSetup wall time:        0    seconds
|  Total wall time:           3    seconds     0.00 hours
