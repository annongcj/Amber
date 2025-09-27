
          -------------------------------------------------------
          Amber 22 PMEMD                              2022
          -------------------------------------------------------

| PMEMD implementation of SANDER, Release 22

|  Compiled date/time: Wed May 25 08:34:30 2022
| Run on 06/07/2022 at 10:54:44

|   Executable path: /home/mahdieh/apps/amber_patch/amber/dist/bin/pmemd.cuda_SPFP
| Working directory: /home/mahdieh/apps/amber_patch/amber/test/middle-scheme/DNA7
|          Hostname: TheraC05

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
| KB of GPU memory in use:     84150
| KB of CPU memory in use:     28278

--------------------------------------------------------------------------------
   4.  RESULTS
--------------------------------------------------------------------------------


 NSTEP =        1   TIME(PS) =      50.001  TEMP(K) =   295.48  PRESS =     0.0
 Etot   =    -29236.0038  EKtot   =      6257.9185  EPtot      =    -35493.9223
 BOND   =       109.0883  ANGLE   =       238.3825  DIHED      =       351.3869
 1-4 NB =       117.1807  1-4 EEL =      -239.3030  VDWAALS    =      4481.9899
 EELEC  =    -40552.6477  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        2   TIME(PS) =      50.002  TEMP(K) =   296.14  PRESS =     0.0
 Etot   =    -29232.8471  EKtot   =      6271.8613  EPtot      =    -35504.7085
 BOND   =       109.7376  ANGLE   =       235.6698  DIHED      =       351.6546
 1-4 NB =       116.4693  1-4 EEL =      -238.5945  VDWAALS    =      4476.6364
 EELEC  =    -40556.2817  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        3   TIME(PS) =      50.003  TEMP(K) =   296.82  PRESS =     0.0
 Etot   =    -29229.0480  EKtot   =      6286.1919  EPtot      =    -35515.2399
 BOND   =       109.7991  ANGLE   =       232.0853  DIHED      =       352.0085
 1-4 NB =       115.6901  1-4 EEL =      -238.0098  VDWAALS    =      4471.5079
 EELEC  =    -40558.3208  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        4   TIME(PS) =      50.004  TEMP(K) =   297.05  PRESS =     0.0
 Etot   =    -29233.8202  EKtot   =      6290.9702  EPtot      =    -35524.7904
 BOND   =       109.3022  ANGLE   =       228.4242  DIHED      =       352.3890
 1-4 NB =       114.9535  1-4 EEL =      -237.6826  VDWAALS    =      4467.0345
 EELEC  =    -40559.2112  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        5   TIME(PS) =      50.005  TEMP(K) =   297.26  PRESS =     0.0
 Etot   =    -29237.4239  EKtot   =      6295.4331  EPtot      =    -35532.8570
 BOND   =       108.4580  ANGLE   =       225.0978  DIHED      =       352.6827
 1-4 NB =       114.3356  1-4 EEL =      -237.6423  VDWAALS    =      4463.1659
 EELEC  =    -40558.9548  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        6   TIME(PS) =      50.006  TEMP(K) =   297.54  PRESS =     0.0
 Etot   =    -29237.3910  EKtot   =      6301.4878  EPtot      =    -35538.8788
 BOND   =       108.0691  ANGLE   =       222.2229  DIHED      =       352.7804
 1-4 NB =       113.8519  1-4 EEL =      -237.8633  VDWAALS    =      4459.1193
 EELEC  =    -40557.0591  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        7   TIME(PS) =      50.007  TEMP(K) =   297.54  PRESS =     0.0
 Etot   =    -29240.1076  EKtot   =      6301.4883  EPtot      =    -35541.5958
 BOND   =       108.6527  ANGLE   =       220.0932  DIHED      =       352.8123
 1-4 NB =       113.4707  1-4 EEL =      -238.2578  VDWAALS    =      4454.9503
 EELEC  =    -40553.3173  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        8   TIME(PS) =      50.008  TEMP(K) =   297.65  PRESS =     0.0
 Etot   =    -29236.4669  EKtot   =      6303.8594  EPtot      =    -35540.3262
 BOND   =       110.0470  ANGLE   =       219.2254  DIHED      =       352.8079
 1-4 NB =       113.2569  1-4 EEL =      -238.9358  VDWAALS    =      4450.7422
 EELEC  =    -40547.4698  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =        9   TIME(PS) =      50.009  TEMP(K) =   297.45  PRESS =     0.0
 Etot   =    -29237.1875  EKtot   =      6299.6436  EPtot      =    -35536.8310
 BOND   =       111.8173  ANGLE   =       219.5156  DIHED      =       352.6481
 1-4 NB =       113.2262  1-4 EEL =      -239.9006  VDWAALS    =      4446.5640
 EELEC  =    -40540.7015  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


 NSTEP =       10   TIME(PS) =      50.010  TEMP(K) =   297.54  PRESS =     0.0
 Etot   =    -29230.8384  EKtot   =      6301.4463  EPtot      =    -35532.2846
 BOND   =       114.0114  ANGLE   =       220.5322  DIHED      =       352.3399
 1-4 NB =       113.3393  1-4 EEL =      -240.9915  VDWAALS    =      4442.4894
 EELEC  =    -40534.0053  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


      A V E R A G E S   O V E R      10 S T E P S


 NSTEP =       10   TIME(PS) =      50.010  TEMP(K) =   297.05  PRESS =     0.0
 Etot   =    -29235.1134  EKtot   =      6291.0300  EPtot      =    -35526.1435
 BOND   =       109.8983  ANGLE   =       226.1249  DIHED      =       352.3510
 1-4 NB =       114.5774  1-4 EEL =      -238.7181  VDWAALS    =      4461.4200
 EELEC  =    -40551.7969  EHBOND  =         0.0000  RESTRAINT  =         0.0000
 ------------------------------------------------------------------------------


      R M S  F L U C T U A T I O N S


 NSTEP =       10   TIME(PS) =      50.010  TEMP(K) =     0.68  PRESS =     0.0
 Etot   =         3.2291  EKtot   =        14.3967  EPtot      =        15.5483
 BOND   =         1.6924  ANGLE   =         6.7631  DIHED      =         0.4835
 1-4 NB =         1.3656  1-4 EEL =         1.0318  VDWAALS    =        12.3743
 EELEC  =         8.1116  EHBOND  =         0.0000  RESTRAINT  =         0.0000
|E(PBS) =         3.2852
 ------------------------------------------------------------------------------

--------------------------------------------------------------------------------
   5.  TIMINGS
--------------------------------------------------------------------------------

|  NonSetup CPU Time in Major Routines:
|
|     Routine           Sec        %
|     ------------------------------
|     Nonbond           0.00   44.48
|     Bond              0.00    0.00
|     Angle             0.00    0.00
|     Dihedral          0.00    0.00
|     Shake             0.00    0.62
|     RunMD             0.00   54.42
|     Other             0.00    0.48
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
|     Elapsed(s) =       0.01 Per Step(ms) =       1.26
|         ns/day =      68.66   seconds/ns =    1258.44
|
|     Average timings for all steps:
|     Elapsed(s) =       0.01 Per Step(ms) =       1.48
|         ns/day =      58.57   seconds/ns =    1475.20
|     -----------------------------------------------------

|  Setup CPU time:            2.84 seconds
|  NonSetup CPU time:         0.01 seconds
|  Total CPU time:            2.85 seconds     0.00 hours

|  Setup wall time:           3    seconds
|  NonSetup wall time:        0    seconds
|  Total wall time:           3    seconds     0.00 hours
