#ifndef _BASE_CUDA_SIMULATION
#define _BASE_CUDA_SIMULATION

#include "gputypes.h"

class base_simulationConst {

public:

  void InitSimulationConst();

public:
  int grid;                 // Grid size
  int gridBits;             // Grid bit count
  int atoms;                // Total number of atoms
  int nQQatoms;             // Number of atoms bearing partial charges
  int nLJatoms;             // Number of atoms capable of any van-der Waals interaction
  int exclAtomIDLimit;      // The first atom in the master topology (original order) known
                            //   to have no excluded Lennard-Jones interactions
  int paddedNumberOfAtoms;  // Atom count padded out to fit grid size
  int SMPCount;             // The number of SMPs on the GPU (needed for initializing certain
                            //   counters)

  // AMBER parameters
  int ntp;                          // AMBER constant pressure setting
  int barostat;                     // AMBER constant pressure barostat
  int alpb;                         // Analytical Linearized Poisson Boltzmann setting
  int igb;                          // Generalized Born overall setting
  int icnstph;                      // constant pH flag
  int icnste;                       // constant Redox potential flag
  int hybridgb;                     // Hybrid Solvent REMD flag
  int ti_mode;                      // Softcore TI mode
  double scnb;                      // 1-4 nonbond scale factor
  double scee;                      // 1-4 electrostatic scale factor
  PMEFloat cut;                     // Nonbond interaction cutoff
  PMEFloat cut2;                    // Nonbond interaction cutoff squared
  PMEFloat cutinv;                  // Nonbonded interaction cutoff inverse
  PMEFloat cut2inv;                 // Nonbonded cutoff inverse squared
  PMEFloat cut3inv;                 // Nonbonded cutoff inverse cubed
  PMEFloat cut6inv;                 // Nonbonded cutoff inverse to the sixth power
  PMEFloat cut3;                    // Nonbonded cutoff, cubed
  PMEFloat cut6;                    // Nonbonded cutoff, sixth power
  PMEFloat es_cutoff;               //
  PMEFloat es_cutoff2;              // Electrostatic and Lennard-Jones non-bonded cutoffs,
  PMEFloat lj_cutoff;               //   possibly squared
  PMEFloat lj_cutoff2;              //
  PMEFloat esCutPlusSkin;           //
  PMEFloat ljCutPlusSkin;           // Electrostatic or Lennard-Jones non-bonded cutoff plus
  PMEFloat esCutPlusSkin2;          //   pair list margins, possibly squared
  PMEFloat ljCutPlusSkin2;          //
  PMEFloat fswitch;                 // VDW force switch constant
  PMEFloat fswitch2;                // VDW force switch constant squared
  PMEFloat fswitch3;                // VDW force switch constant cubed
  PMEFloat fswitch6;                // VDW force switch constant to the sixth power
  PMEFloat invfswitch6cut6;         // VDW force switch constant, 1 / (fswitch^6*cut^6)
  PMEFloat invfswitch3cut3;         // VDW force switch constant, 1 / (fswitch^3*cut^3)
  PMEFloat cut6invcut6minfswitch6;  // VDW force switch constant, cut6 / (cut6-fswitch6)
  PMEFloat cut3invcut3minfswitch3;  // VDW force switch constant, cut3 / (cut3-fswitch3)
  PMEFloat cutPlusSkin;             // Nonbond cutoff plus skin
  PMEFloat cutPlusSkin2;            // Nonbond cutoff plus skin, squared
  int efn;                          // Normalize electric field vectors
  PMEFloat efx;                     // Electric field x vector
  PMEFloat efy;                     // Electric field y vector
  PMEFloat efz;                     // Electric field z vector
  PMEFloat efphase;                 // Electric field spatial phase
  PMEFloat effreq;                  // Electric field time frequency
  double dielc;                     // Dielectric constant
  double gamma_ln;                  // Langevin integration parameter
  double c_ave;                     // Langevin integration parameter
  double c_implic;                  // Langevin integration parameter
  double c_explic;                  // Langevin integration parameter
  bool bUseVlimit;                  // Use vlimit flag
  double vlimit;                    // vlimit
  double tol;                       // SHAKE tolerance
  double massH;                     // Hydrogen mass
  aligned_double invMassH;          // Inverse hydrogen mass
  PMEFloat gb_alpha;                // Generalized Born (GB) parameter
  PMEFloat gb_beta;                 // GB parameter
  PMEFloat gb_gamma;                // GB parameter
  PMEFloat gb_kappa;                // GB parameter
  PMEFloat gb_kappa_inv;            // GB derived parameter
  PMEFloat gb_cutoff;               // GB cutoff
  PMEFloat gb_fs_max;               // GB something or other
  PMEFloat gb_neckscale;            // GB neck scaling factor
  PMEFloat gb_neckcut;              // GB neck cutoff
  PMEFloat gb_neckoffset;           // GB neck offset for LUT index
  PMEFloat rgbmax;                  // GB Born radius cutoff
  PMEFloat rgbmax1i;                // Inverse GB Radius cutoff
  PMEFloat rgbmax2i;                // Inverse GB Radius cutoff
  PMEFloat rgbmaxpsmax2;            // GB derived quantity
  PMEFloat intdiel;                 // GB interior dielectric
  PMEFloat extdiel;                 // GB exterior dielectric
  PMEFloat alpb_alpha;              // GB parameter
  PMEFloat alpb_beta;               // GB derived parameter
  PMEFloat extdiel_inv;             // GB derived parameter
  PMEFloat intdiel_inv;             // GB derived parameter
  PMEFloat saltcon;                 // GB salt conductivity
  PMEFloat surften;                 // GB surface tension
  PMEFloat offset;                  // GB Born radius offset
  PMEFloat arad;                    // GB atomic radius
  PMEFloat one_arad_beta;           // GB derived parameters
  int gbion;                        // GB with explict ion 0 or 2
  PMEFloat gi_coef_1;               // gbion=2, exp(-dij^2/(4*Ri*Rj*gi_coef_1)), ion-solute
  PMEFloat gi_coef_2;               // ion-ion
  PMEFloat intdiel_ion_1;           // intdiel coef for ion-solute
  PMEFloat intdiel_ion_2;           // ion-ion
  PMEFloat gb_neckscale_ion_1;      // gb_neckscale for ion-solute
  PMEFloat gb_neckscale_ion_2;      // ion-ion
  PMEFloat gi_coef_1_p;             // gbion=3, exp(-dij^2/(4*Ri*Rj*gi_coef_1)), cation-solute
  PMEFloat gi_coef_1_n;             // gbion=3, exp(-dij^2/(4*Ri*Rj*gi_coef_1)), anion-solute
  PMEFloat gi_coef_2_pp;            // cation-cation
  PMEFloat gi_coef_2_pn;            // cation-anion
  PMEFloat gi_coef_2_nn;            // anion-anion
  PMEFloat intdiel_ion_1_p;         // intdiel coef for cation-solute
  PMEFloat intdiel_ion_1_n;         // intdiel coef for anion-solute
  PMEFloat intdiel_ion_2_pp;        // cation-cation
  PMEFloat intdiel_ion_2_pn;        // cation-anion
  PMEFloat intdiel_ion_2_nn;        // anion-anion
  PMEFloat gb_neckscale_ion_1_p;    // gb_neckscale for cation-solute
  PMEFloat gb_neckscale_ion_1_n;    // gb_neckscale for anion-solute
  PMEFloat gb_neckscale_ion_2_pp;   // cation-cation
  PMEFloat gb_neckscale_ion_2_pn;   // cation-anion
  PMEFloat gb_neckscale_ion_2_nn;   // anion-anion

  // PME parameters
  PMEDouble a;                      // PBC x cell length
  PMEDouble b;                      // PBC y cell length
  PMEDouble c;                      // PBC z cell length
  PMEFloat af;                      // PBC single-precision x cell length
  PMEFloat bf;                      // PBC single-precision y cell length
  PMEFloat cf;                      // PBC single-precision z cell length
  PMEFloat alpha;                   // PBC Alpha
  PMEFloat beta;                    // PBC Beta
  PMEFloat gamma;                   // PBC Gamma
  PMEFloat pbc_box[3];              // PBC Box
  PMEFloat reclng[3];               // PBC reciprocal cell lengths
  PMEFloat cut_factor[3];           // PBC cut factors
  aligned_double ucell[3][3];       // PBC cell coordinate system
  PMEFloat ucellf[3][3];            // Single precision PBC cell coordinate system
  aligned_double recip[3][3];       // PBC reciprocal cell coordinate system
  PMEFloat recipf[3][3];            // Single precision PBC reciprocal cell coordinate system
  PMEFloat uc_volume;               // PBC cell volume
  PMEFloat uc_sphere;               // PBC bounding sphere
  PMEFloat pi_vol_inv;              // PBC/PME constant
  PMEFloat fac;                     // PME Ewald factor
  PMEFloat fac2;                    // PME Ewald factor x 2
  bool is_orthog;                   // PBC cell orthogonality flag
  int nfft1;                        // x real-space mesh grid size
  int nfft2;                        // y real-space mesh grid size
  int nfft3;                        // z real-space mesh grid size
  double dnfft1;                    //
  double dnfft2;                    // Double-precision variants of the mesh grid sizes
  double dnfft3;                    //
  int nMeshSubRegion;               // Number of PME grid subsections
  int xSmemWidth;                   // Widths of the grid actually stored by the interpolation
  int ySmemWidth;                   //   kernel in __shared__ memory, measured along the x, y,
  int zSmemWidth;                   //   and z dimensions
  int xySmemWidth;                  // Useable size plus padding for each xy slab of the
                                    //   __shared__ memory grid used by the interpolation
                                    //   kernel.  Each slab is padded by an additional amount
                                    //   to line up memory banks for interpolating each atom.
  int nfft1xnfft2;                  // Product of nfft1 and nfft2
  int nfft1xnfft2xnfft3;            // Product of nfft1, nfft2 and nfft3
  int fft_x_dim;                    // x fft dimension
  int fft_y_dim;                    // y fft dimension
  int fft_z_dim;                    // z fft dimension
  int fft_quarter_z_dim_m1;         // One quarter of the z fft dimension, minus one
  int fft_y_dim_times_x_dim;        // Product of x and y dimension in the complex grid
  int fft_x_y_z_dim;                // Product of x, y, and z dimensions in the complex grid
  int fft_x_y_z_quarter_dim;        // One quarter of fft_x_y_z_dim (this works because the
                                    //   z dimension, as well as x and y dimensions, are all
                                    //   guaranteed to be multiples of 4)
  int nf1;                          // x scalar sum coefficient
  int nf2;                          // y scalar sum coefficient
  int nf3;                          // z scalar sum coefficient
  int n2Offset;                     // Offset to y prefac data
  int n3Offset;                     // Offset to z prefac data
  int nSum;                         // Total prefac data
  int pmeOrder;                     // PME interpolation order
  int orderMinusOne;                // Interpolation order - 1
  int XYZStride;                    // Stride of PME buffers
  PMEFloat ew_coeffSP;              // Ewald coefficient (fp32 for SPFP mode)
  double ew_coeff;                  // Ewald coefficient in double precision
  PMEFloat ew_coeff2;               // Ewald coefficient squared
  PMEFloat negTwoEw_coeffRsqrtPI;   // Switching function constant
  PMEFloat4* pErfcCoeffsTable;      // d/dr [erfc(a * r)/r] spline coefficients

  // IPS paramers
  bool bIPSActive;                  // Flag to indicate IPS is active
  PMEFloat rips;                    // Radius of IPS local region
  PMEFloat rips2;                   // rips^2
  PMEFloat ripsr;                   // 1/rips
  PMEFloat rips2r;                  // 1/rips^2
  PMEFloat rips6r;                  // 1/rips^6
  PMEFloat rips12r;                 // 1/rips^12
  PMEDouble eipssnb;                // IPS self VDW energy
  PMEDouble eipssel;                // IPS self-electrostatic energy
  PMEDouble virips;                 // IPS self-virial energy
  PMEDouble EIPSEL;                 // IPS exclusion self electrostatic energy
  PMEDouble EIPSNB;                 // IPS exclusion self Nonbond energy

  // Electrostatic IPS parameters:
  PMEFloat aipse0;
  PMEFloat aipse1;
  PMEFloat aipse2;
  PMEFloat aipse3;
  PMEFloat pipsec;
  PMEFloat pipse0;
  PMEFloat bipse1;
  PMEFloat bipse2;
  PMEFloat bipse3;

  // Dispersion IPS parameters:
  PMEFloat aipsvc0;
  PMEFloat aipsvc1;
  PMEFloat aipsvc2;
  PMEFloat aipsvc3;
  PMEFloat pipsvcc;
  PMEFloat pipsvc0;
  PMEFloat bipsvc1;
  PMEFloat bipsvc2;
  PMEFloat bipsvc3;

  // Repulsion IPS parameters:
  PMEFloat aipsva0;
  PMEFloat aipsva1;
  PMEFloat aipsva2;
  PMEFloat aipsva3;
  PMEFloat pipsvac;
  PMEFloat pipsva0;
  PMEFloat bipsva1;
  PMEFloat bipsva2;
  PMEFloat bipsva3;

  // AMD parameters:
  int  iamd;                        // Flag to activate accelerated Molecular Dynamics (aMD)
                                    //   of various flavors
  int  w_amd;                       // Weight to use in the aMD
  int  iamdlag;                     // Commence aMD after a brief equilibration period (?)
  int  amd_print_interval;          // Interval to print aMD energies (?)
  int  AMDNumRecs;                  //
  int  AMDNumLag;                   //
  PMEDouble AMDtboost;              // Total aMD boost being applied (varies step to step,
                                    //   used to calculate the force to apply)
  PMEDouble AMDfwgt;                // GPU analog of fwgt (calculated from alphaD and EV)
  PMEUllInt* pAMDEDihedral;         //
  PMEDouble* pAMDfwgtd;             //
  PMEDouble amd_EthreshP;           // Average total potential energy threshold
  PMEDouble amd_alphaP;             // Inverse strength boost factor for total potential energy
  PMEDouble amd_EthreshD;           // Average dihedral energy threshold
  PMEDouble amd_alphaD;             // Inverse strength boost factor for the dihedral energy
  PMEDouble amd_temp0;              // Starting temperature for aMD
  PMEDouble amd_EthreshP_w;         //
  PMEDouble amd_alphaP_w;           // Weighted versions of the above
  PMEDouble amd_EthreshD_w;         //
  PMEDouble amd_alphaD_w;           // Sign of the aMD weight (?)
  PMEDouble amd_w_sign;

  // GaMD parameters:
  int igamd;                        // Flag to activate Gaussian-accelerated MD
  int igamd0;                       // Flag to active potential collected in GaMD equilibration
  int tspanP;                       //
  int tspanD;                       //
  int tspan;                        //
  int igamdlag;                     //
  int gamd_print_interval;          //
  int GaMDNumRecs;                  //
  int GaMDNumLag;                   // All of these variables have analogs in the aMD control
  PMEDouble GaMDtboost;             //   variables listed above.  See above for descriptions.
  PMEDouble GaMDfwgt;               //
  PMEDouble GaMDfwgtb;
  PMEUllInt* pGaMDEDihedral;        //
  PMEDouble* pGaMDfwgtd;            //
  PMEDouble gamd_EthreshP;          //
  PMEDouble gamd_kP;                //
  PMEDouble gamd_EthreshD;          //
  PMEDouble gamd_kD;                //
  PMEDouble gamd_EthreshB;          //
  PMEDouble gamd_kB;                //
  PMEDouble gamd_temp0;             //

  // scaledMD parameters:
  int scaledMD;                        // Flag to activate scaled molecular dynamics
  PMEDouble scaledMD_lambda;           //
  PMEDouble scaledMD_energy;           //
  PMEDouble scaledMD_weight;           //
  PMEDouble scaledMD_unscaled_energy;  //

  // For relaxation dynamics
  int first_update_atom;           // First atom whose crds should be updated


    // SGLD parameters:
    int                         isgld;
    int                         nsgatom;
    PMEDouble                   sgavg;
    PMEDouble                   sgavp;
    PMEDouble                   sgft;
    PMEDouble                   sgff;
    PMEDouble                   sgfg;
    PMEDouble                   fsgldg;
    PMEDouble                   tempsg;
    PMEDouble                   tsgset;
    PMEDouble                   sgmsum;
    PMEDouble                   sgrndf;
    // SGLD pointers
    int*                        pWsg;                              // SGLD weighter
    double*                     pX0sg;                              // local average X 
    double*                     pY0sg;                              // local average Y 
    double*                     pZ0sg;                              // local average Z 
    double*                     pX1sg;                              // local average X 
    double*                     pY1sg;                              // local average Y 
    double*                     pZ1sg;                              // local average Z 
    double*                     pX2sg;                              // local average X1 
    double*                     pY2sg;                              // local average Y1 
    double*                     pZ2sg;                              // local average Z1 
    double*                     pRsgX;                              // local average X random forces
    double*                     pRsgY;                              // local average Y random forces
    double*                     pRsgZ;                              // local average Z random forces
    double*                     pFPsg;                              // local average of momentums time forces
    double*                     pPPsg;                              // local average of momentums time momentums
    int*                        pImageWsg;                              // SGLD weighter
    double*                     pImageX0sg;                         // Image x local average 
    double*                     pImageY0sg;                         // Image y local average 
    double*                     pImageZ0sg;                         // Image z local average 
    double*                     pImageX1sg;                         // Image x local average 
    double*                     pImageY1sg;                         // Image y local average 
    double*                     pImageZ1sg;                         // Image z local average 
    double*                     pImageX2sg;                         // Image x1 local average 
    double*                     pImageY2sg;                         // Image y1 local average 
    double*                     pImageZ2sg;                         // Image z1 local average 
    double*                     pImageRsgX;                         // Image local average X random forces
    double*                     pImageRsgY;                         // Image local average Y random forces
    double*                     pImageRsgZ;                         // Image local average Z random forces
    double*                     pImageFPsg;                         // Image local average of forces time momentum
    double*                     pImagePPsg;                         // Image local average of momentum time momentum
    int*                        pImageWsg2;                              // SGLD weighter
    double*                     pImageX0sg2;                         // Image x 
    double*                     pImageY0sg2;                         // Image y 
    double*                     pImageZ0sg2;                         // Image z 
    double*                     pImageX1sg2;                         // Image x local average 
    double*                     pImageY1sg2;                         // Image y local average 
    double*                     pImageZ1sg2;                         // Image z local average 
    double*                     pImageX2sg2;                         // Image x1 local average 
    double*                     pImageY2sg2;                         // Image y1 local average 
    double*                     pImageZ2sg2;                         // Image z1 local average 
    double*                     pImageRsgX2;                         // Image local average X random forces
    double*                     pImageRsgY2;                         // Image local average Y random forces
    double*                     pImageRsgZ2;                         // Image local average Z random forces
    double*                     pImageFPsg2;                         // Image local average of forces time momentum
    double*                     pImagePPsg2;                         // Image local average of momentum time momentum
   // SGLD calculation pointer
    SGLDAverage*                pSGLDAverage;   

  // NTP stuff
  NTPData* pNTPData;               // PME NTP mutable values

  // Constant pH stuff
  double* pChargeRefreshBuffer;    // Pointer to new charges

  // Softcore TI stuff
  PMEDouble AFElambda[2];          // The value of lambda (mixing term) to use in (soft-core)
                                   //   Thermodynamic Integration
  PMEFloat AFElambdaSP[2];         //
  PMEDouble scalpha;               // Taking L = lambda, E = Lennard Jones epsilon, and
                                   //   S = Lennard Jones sigma, the equation for the softcore
                                   //   potential is:
                                   //
                                   // V(sc) = 4E(1-L)[ (1/(alpha*L + (r/S)^6)^2) -
                                   //                  (1/(alpha*L + (r/S)^6)) ]
                                   //
  PMEDouble scbeta;                // Equivalent scaling factor for electrostatic interactions
                                   //   in softcore TI
  PMEFloat scalphaSP;              // Floating-point equivalents of scalpha and scbeta, used
  PMEFloat scbetaSP;               //   in single-precision mode
  PMEFloat TIsigns[2];             //
  int bar_intervall;               // Frequency of calculating mbar values when ifmbar = 1
  int bar_states;                  // Number of mbars to calculate
  int bar_stride;
  int TIPaddedLinearAtmCnt;
  int TIlinearAtmCnt;
  int AFE_recip_region;
  int ifmbar;

#ifdef MPI
  // MPI stuff
  int  minLocalAtom;               // First local atom
  int  maxLocalAtom;               // Last local atom
#endif

  // Atom stuff
  double* pAtomX;                  // Atom X coordinates
  double* pAtomY;                  // Atom Y coordinates
  double* pAtomZ;                  // Atom Z coordinates
  double* pOldAtomX;               // Old Atom X position
  double* pOldAtomY;               // Old Atom Y position
  double* pOldAtomZ;               // Old Atom Z position
  double* pShakeOldAtomX;          // Old Atom X position added by zhf, used by middle-scheme
  double* pShakeOldAtomY;          // Old Atom Y position velocity recalculation
  double* pShakeOldAtomZ;          // Old Atom Z position
  PMEFloat2* pAtomXYSP;            // Single Precision Atom X and Y coordinates
  PMEFloat* pAtomZSP;              // Single Precision Atom Z coordinates
  PMEFloat2* pAtomSigEps;          // Atom nonbond parameters
  unsigned int* pAtomLJID;         // Atom Lennard-Jones index
  PMEFloat* pAtomS;                // Atom S Parameter
  PMEFloat* pAtomRBorn;            // Atom Born Radius
  double* pAtomCharge;             // Atom charges
  PMEFloat* pAtomChargeSP;         // Single precision atom charges
  PMEFloat2* pAtomChargeSPLJID;    // Single precision atom charges, Lennard Jones atom types
  double* pAtomMass;               // Atom masses
  double* pAtomInvMass;            // Atom inverse masses
  int* pTIRegion;                  // Softcore TI Region
  int* pTILinearAtmID;             // Softcore TI linear atom IDs for use in vector exchange
  unsigned long long int* pBarTot; // Holds final AFE MBar values
  double* pBarLambda;              // AFE Holds MBar Lambda constants
  PMEDouble* pReff;                // Effective Born Radius
  PMEFloat* pPsi;                  // GB intermediate psi value
  PMEFloat* pReffSP;               // Single Precision Effective Born Radius
  PMEFloat* pTemp7;                // Single Precision Born Force
  unsigned int* pAtomIonMask;               // Atom is ion (>=1) or not (=0)

#ifdef MPI
  // Remote data accumulator
  PMEAccumulator* pPeerAccumulator;
#endif
  double* pVelX;                 // Atom X velocities
  double* pVelY;                 // Atom Y velocities
  double* pVelZ;                 // Atom Z velocities
  double* pLVelX;                // Atom X last velocities
  double* pLVelY;                // Atom Y last velocities
  double* pLVelZ;                // Atom Z last velocities
  int  maxNonbonds;              // Maximum nonbond interaction buffers
  PMEFloat* pXMax;               // Recentering maximum x
  PMEFloat* pYMax;               // Recentering maximum y
  PMEFloat* pZMax;               // Recentering maximum z
  PMEFloat* pXMin;               // Recentering minimum x
  PMEFloat* pYMin;               // Recentering minimum y
  PMEFloat* pZMin;               // Recentering minimum z

  // Neighbor List stuff
  PMEFloat skinnb;                // Input Nonbond skin
  PMEFloat skinPermit;            // Input skin permittivity (the fraction of the pair list
                                  //   margin over which particles may travel before triggering
                                  //   a pair list rebuild)
  unsigned int ljexist[4];        // Bit mask to indicate whether Lennard-Jones potentials
                                  //   exist between an atom type and any others
  unsigned int* pImageIndex;      // Image # for spatial sort
  unsigned int* pImageIndex2;     // Image # for spatial sort
  unsigned int* pImageHash;       // Image hash for spatial sort
  unsigned int* pImageHash2;      // Image hash for spatial sort
  unsigned int* pImageAtom;       // Image atom #
  unsigned int* pImageAtom2;      // Image atom #
  unsigned int* pImageAtomLookup; // Original atom lookup table
  PMEFloat2* pAtomXYSaveSP;       // Saved atom coordinates from neighbor list generation
  PMEFloat* pAtomZSaveSP;         // Saved atom coordinates from neighbor list generation
  double* pImageX;                // Image x coordinates
  double* pImageY;                // Image y coordinates
  double* pImageZ;                // Image z coordinates
  double* pImageVelX;             // Image x velocities
  double* pImageVelY;             // Image y velocities
  double* pImageVelZ;             // Image z velocities
  double* pImageLVelX;            // Image last x velocities
  double* pImageLVelY;            // Image last y velocities
  double* pImageLVelZ;            // Image last z velocities
  double* pImageMass;             // Image masses
  double* pImageInvMass;          // Image inverse masses
  double* pImageCharge;           // Image charges
  PMEFloat2* pImageSigEps;        // Image sigma/epsilon data for nonbond interactions
  unsigned int* pImageLJID;       // Image Lennard-Jones index
  unsigned int* pImageCellID;     // Image cell ID for calculating local coordinate system
  int* pImageTIRegion;            // Image Softcore TI data
  int* pImageTILinearAtmID;       // Atom indexes for the linear atoms used in syncing vectors
                                  //   each step
  unsigned int* pUpdateIndex;     // The updated position for each index for use with syncing
                                  //   vectors at each step
  double* pImageX2;               // Image x coordinates
  double* pImageY2;               // Image y coordinates
  double* pImageZ2;               // Image z coordinates
  double* pImageVelX2;            // Image x velocities
  double* pImageVelY2;            // Image y velocities
  double* pImageVelZ2;            // Image z velocities
  double* pImageLVelX2;           // Image last x velocities
  double* pImageLVelY2;           // Image last y velocities
  double* pImageLVelZ2;           // Image last z velocities
  double* pImageMass2;            // Image masses
  double* pImageInvMass2;         // Image inverse masses
  double* pImageCharge2;          // Image charges
  PMEFloat2* pImageSigEps2;       // Image sigma/epsilon data for nonbond interactions
  unsigned int* pImageLJID2;      // Image Lennard-Jones index
  unsigned int* pImageCellID2;    // Image cell ID for calculating local coordinate system
  int* pImageTIRegion2;           // Image Softcore TI data
  int* pImageTILinearAtmID2;      // Atom indexes for the linear atoms used in syncing coords
  int2* pImageNMRCOMDistanceID;        // Remapped NMR COM Distance i, j
  int2* pImageNMRCOMDistanceCOM;       // Remapped NMR Distance COM range .x to .y
  int2* pImageNMRCOMDistanceCOMGrp;    // Remapped NMR Distance COM group range .x to .y
  int2* pImageNMRr6avDistanceID;       // Remapped NMR r6av Distance i, j
  int2* pImageNMRr6avDistancer6av;     // Remapped NMR Distance r6av range .x to .y
  int2* pImageNMRr6avDistancer6avGrp;  // Remapped NMR Distance r6av group range .x to .y

  int2* pImageNMRCOMAngleID1;          // Remapped NMR COM Angle i, j
  int* pImageNMRCOMAngleID2;           // Remapped NMR COM Angle k
  int2* pImageNMRCOMAngleCOM;         // Remapped NMR Angle COM range .x to .y
  int2* pImageNMRCOMAngleCOMGrp;       // Remapped NMR Angle COM group range .x to .y

  int4* pImageNMRCOMTorsionID1;          // Remapped NMR COM Torsion i, j, k, l
  int2* pImageNMRCOMTorsionCOM;         // Remapped NMR COM Torsion COM range .x to .y
  int2* pImageNMRCOMTorsionCOMGrp;       // Remapped NMR COM Torsion  COM group range .x to .y
  int4* pImageShakeID;                 // Remapped Shake Atom ID
  int4* pImageFastShakeID;             // Remapped Fast Shake Atom ID
  int* pImageSlowShakeID1;             // Remapped Slow Shake central Atom ID
  int4* pImageSlowShakeID2;            // Remapped Slow Shake hydrogen Atom IDs
  int* pImageSlowTIShakeID1;           // Remapped Slow Shake central Atom ID
  int4* pImageSlowTIShakeID2;          // Remapped Slow Shake hydrogen Atom IDs
  int4* pImageSlowTIShakeID3;          // Remapped Slow Shake hydrogen Atom IDs
  int4* pImageSolventAtomID;           // Remapped solvent molecules/ions
  int* pImageSoluteAtomID;             // Remapped solute atoms

  // Bonded term work units
  int bondWorkUnits;               // The number of bond work units
  int clearQBWorkUnits;            // The number of charge buffer clearing work units (the
                                   //   charge buffer is cleared while bonded interactions are
                                   //   being computed to backfill the GPU)
  int bondWorkUnitClusterSize;     // The number of consecutive work units that each thread
                                   //   block will take.  This number is determined by the
                                   //   CPU so that a maximum of 65536 blocks can be launched
                                   //   and still accommodate all of the work units.  Different
                                   //   work units will take different amounts of computation,
                                   //   but with up to 65536 clusters the granularity should
                                   //   offer very good GPU efficiency.
  int BwuCnstCount;                // The length of the atom positional constraints parameter
                                   //   array, as kept by the bond work units.  This array is
                                   //   padded, rounding up to the nearest 32 for every bond
                                   //   work unit that handles such constraints.  Other counts
                                   //   might need to be stored in future development.
  unsigned int *pBwuInstructions;  //
  unsigned int *pBwuBondID;        //
  unsigned int *pBwuAnglID;        // Pointers into the bond work unit unsigned integer array:
  unsigned int *pBwuDiheID;        //   the instructions are found at the head of the array,
  unsigned int *pBwuCmapID;        //   other lists of bit packed ID numbers later.  See the
  unsigned int *pBwuQQxcID;        //   descriptions in bondRemap.cpp for more explanations.
  unsigned int *pBwuNB14ID;        //   The pointer pBwuCmapID leads into twice as much data
  unsigned int *pBwuNMR2ID;        //   as the others--64 unsigned ints as opposed to just 32,
  unsigned int *pBwuNMR3ID;        //   read as 32 threads' worth of IJKL indices followed by
  unsigned int *pBwuNMR4ID;        //   32 threads' worth of energy map / M atom indices.
  unsigned int *pBwuUreyID;        //
  unsigned int *pBwuCImpID;        //
  unsigned int *pBwuCnstID;        //
  unsigned int *pBwuCnstUpdateIdx; // Pointer, also into the bond work unit integer data array,
                                   //   detailing the correspondence of parameters for each
                                   //   bond work unit's constraint objects.  Other origin
                                   //   arrays may be needed in future development.
  unsigned int *pBwuBondStatus;    //
  unsigned int *pBwuAnglStatus;    //
  unsigned int *pBwuDiheStatus;    // Pointers into the bond work unit unsigned integer data
  unsigned int *pBwuCmapStatus;    //   array, detailing the status of each object with respect
  unsigned int *pBwuQQxcStatus;    //   to Thermodynamic Integration.  Each is a packed string
  unsigned int *pBwuNB14Status;    //   of bits, the high eight being discarded, next eight
  unsigned int *pBwuNMR2Status;    //   giving the TI region, next SCBond, and the low eight
  unsigned int *pBwuNMR3Status;    //   giving the CVBond status.  These status integers start
  unsigned int *pBwuNMR4Status;    //   life as signed ints, but since they're all positive
  unsigned int *pBwuUreyStatus;    //   the conversion is safe.
  unsigned int *pBwuCimpStatus;    //
  unsigned int *pBwuCnstStatus;    //
  PMEDouble *pBwuBondLambda;       //
  PMEDouble *pBwuAnglLambda;       //
  PMEDouble *pBwuDiheLambda;       //
  PMEDouble *pBwuCmapLambda;       //
  PMEDouble *pBwuQQxcLambda;       // Pointers into the bond work unit double data array.
  PMEDouble *pBwuNB14Lambda;       //   These encode lambda values for Thermodynamic
  PMEDouble *pBwuNMR2Lambda;       //   integration.  Dynamic lambdas requrie that this array
  PMEDouble *pBwuNMR3Lambda;       //   be updated.
  PMEDouble *pBwuNMR4Lambda;       //
  PMEDouble *pBwuUreyLambda;       //
  PMEDouble *pBwuCimpLambda;       //
  PMEDouble *pBwuCnstLambda;       //
  PMEDouble2 *pBwuBond;            // Pointers into the bond work unit double2 data array: this
  PMEDouble2 *pBwuAngl;            //   contains details of all the bonded terms to be computed
  PMEDouble2 *pBwuNMR2;            //   and most attributes have analogs (change pBwu to "p" to
  PMEDouble2 *pBwuNMR3;            //   find them, with the exceptions that NMR{2,3,4} replace
  PMEDouble2 *pBwuNMR4;            //   NMRDistance, NMRAngle, and NMRTorsion, respectively.
  PMEDouble2 *pBwuUrey;            //   Also, impDihedral has been changed to CImp (CHARMM
  PMEDouble2 *pBwuCImp;            //   improper).  The NMR{2,3,4} pointers lead into three
  PMEDouble2 *pBwuCnst;            //   times as much space as any of the others, which by and
                                   //   large hold harmonic stiffness and equilibrium
                                   //   constants.  The NMR details contain R1/R2, then R3/R4,
                                   //   then K2/K3, but the indexing from pBwuNMR??ID is done
                                   //   so that the first read gets R1/R2, a second read offset
                                   //   by GRID elements gets R3/R4, and a third offset by
                                   //   2*GRID gets K2/K3.
  PMEFloat *pBwuCharges;           // Atomic partial charges for all atoms in the work units.
                                   //   This is somewhat redundant, as the electrostatic
                                   //   subimage array will contain the charges as well, but
                                   //   loads will not always proceed from there: work units
                                   //   requiring double precision coordinates will need to
                                   //   import charges separately--might as well be from a
                                   //   pre-arranged list.  Pointer into a larger array.
  PMEFloat *pBwuDihe3;             // In SP*P modes, these pointers hold fp32 representations
  PMEFloat2 *pBwuDihe12;           //   of the dihedral parameters--SP*P will do dihedrals in
                                   //   FP32 after reading double precision coordinates.  These
                                   //   point into the same larger array as pBwuLJnb14.
                                   //   In DPFP mode, all the PMEFloats are doubles anyway.
  PMEFloat *pBwuEEnb14;            // Charge-charge non-bonded 1:4 scaling factors--both
                                   //   charges and the attenuations are folded in here.
  PMEFloat2 *pBwuLJnb14;           // Parameters for computing the Lennard-Jones interactions
                                   //   of non-bonded 1:4 interactions: ACOEF and BCOEF.

  // Atom shuttling
  int nShuttle;            // The maximum number of atoms that can be shuttled.  The
                           //   shuttle arrays will be sized based on this number, but
                           //   their sizes may have other factors.
  int ShuttleType;         // The type of shuttling to perform
  double* pDataShuttle;    // Array to hold information on atoms to be shuttled between
                           //   the CPU and GPU.  This array will have eight-byte
                           //   elements in double, as the most "official" coordinates are
                           //   given in double precision and forces (stored in eight-byte
                           //   long long ints) willbe converted to double precision on the
                           //   GPU to fit into the Fortran code.
  int* pShuttleTickets;    // List of atoms that will get tickets on the shuttle to and
                           //   from the CPU.

  //MCWAT info
  int stericBitCellWidth;   // bit grid width (x dim)
  int stericBitCellHeight;  // bit grid width (y dim)
  int stericBitCellDepth;   // bit grid width (z dim) x y z should multiply to 32 for int32
  int stericMaxXVxl;        // steric grid in x direction
  int stericMaxYVxl;        // steric grid in y direction
  int stericMaxZVxl;        // steric grid in z direction
  PMEFloat stericGridSpacing; // steric grid spacing size
  int waterAtomCount;      // number of atoms in water.
  int coarseMaxXvxl;        // Coarse grid size in x direction
  int coarseMaxYvxl;        // Coarse grid size in y direction
  int coarseMaxZvxl;        // Coarse grid size in z direction
  int coarseMaxAtoms;        // Max atoms per grid cell
  PMEFloat coarseGridSpacing; // Coarse grid spacing size
  int mcressize;            // Number of samples to do for coarse grid
  int* pCoarseGrid;          // x,y,z,atom coarse grid
  int* pRandomWater;        // List of water molecules
  int* pStericEmptyVoxels;  // List of empty voxels in steric grid
  unsigned int* pStericGrid;        // List of sterically mapped atoms
  int* pStericThreadCount;  // List of how many open cells in steric grid for a thread
  PMEFloat* pStericRadius;   // List of radius of atoms
  int* pVoxelOffset;        // List of radius of atoms
  PMEFloat* pCoarseTrans;     // Translation for crd in coarse grid
  PMEFloat* pCoarseGridEne;     // Stores the coarse grid energies


  //NEB info
  int neb_nbead;                       // Total number of NEB replicas.
  int beadid;                          // NEB replica ID.
  int nattgtrms;                       // Number of atoms used for ENB force calculation.
  int nattgtfit;                       // Number of atoms used for NEB structure fitting.
  int last_neb_atm;                    // Last atom involved in NEB force calculation or
                                       //   structure fitting.
  int skmin;                           // Spring constants for NEB.
  int skmax;                           // Spring constants for NEB.
  int tmode;                           // Tangent mode for NEB.
  int* pRMSMask;                       // Mask of atoms used for ENB force calculation.
  int* pFitMask;                       // Mask of atoms used for NEB structure fitting.
  int* pAtmIdx;
  PMEDouble* pNEBEnergyAll;            // Energy of the replicas in NEB.
  double* pTangents;                   // Tangents to the path.
  double* pSpringForce;                // Artificial spring forces.
  double* pNEBForce;                   // NEB forces.
  double* pNextDataShuttle;            // Next neighboring replica's coordinate.
  double* pPrevDataShuttle;            // Previous neighboring replica's coordinate.
  double* pKabschCOM;
  double* pPrevKabsch;
  double* pNextKabsch;
  double* pSelfCOM;
  double* pNextCOM;
  double* pPrevCOM;
  double* pKabsch2;
  double* pPrevKabsch2;
  double* pNextKabsch2;
  double* pDataSPR;
  double* ptotFitMass;                 // Pointer to total mass of atoms in RMSFit array.
  double* pNorm;
  double* pDotProduct1;
  double* pDotProduct2;
  double* pSpring1;
  double* pSpring2;
  double* pRotAtm;
  double* pPrevRotAtm;
  double* pNextRotAtm;

  // Direct space decomposition
  uint2* pNLNonbondCellStartEnd;           // Nonbond cell boundaries pointer
  int  srTableLength;                      // Length of the list describing the reciprocal
                                           //   space mesh buffers
  int  cells;                              // Total number of nonbond cells
  int  xcells;                             // Number of x cells
  int  ycells;                             // Number of y cells
  int  zcells;                             // Number of z cells
  PMEDouble xcell;                         // x cell dimension
  PMEDouble ycell;                         // y cell dimension
  PMEDouble zcell;                         // z cell dimension
  PMEDouble minCellX;                      // Minimum x cell coordinate
  PMEDouble minCellY;                      // Minimum y cell coordinate
  PMEDouble minCellZ;                      // Minimum z cell coordinate
  PMEDouble maxCellX;                      // Maximum x cell coordinate
  PMEDouble maxCellY;                      // Maximum y cell coordinate
  PMEDouble maxCellZ;                      // Maximum z cell coordinate
  aligned_double oneOverXcells;            // Fractional x cell dimension
  aligned_double oneOverYcells;            // Fractional y cell dimension
  aligned_double oneOverZcells;            // Fractional z cell dimension
  PMEFloat oneOverXcellsf;                 // Single precision fractional x cell dimension
  PMEFloat oneOverYcellsf;                 // Single precision fractional y cell dimension
  PMEFloat oneOverZcellsf;                 // Single precision fractional z cell dimension
  PMEFloat cell;                           // Minimum cell dimension
  PMEFloat nonbond_skin;                   // Effective nonbond skin
  unsigned int* pFinishedBlocksCounters;   // Counters of finished blocks in
                                           // cudaDeviceSynchronize-less kernels
  unsigned int* pKernelCompletionFlags;    // Completion flags for
                                           // cudaDeviceSynchronize-less kernels
  PMEFloat one_half_nonbond_skin_squared;  // Skin test atom movement threshold
  bool* pNLbSkinTestFail;                  // Skin test result buffer
  unsigned int* pNLCellHash;               // Spatial ordering hash for within cells
  unsigned int* pNLTotalOffset;            // Pointer to total offset
  int NLMaxTotalOffset;                    // Maximum available exclusion masks
  unsigned int* pNLAtomList;               // Pointer to neighbor list atoms
  unsigned int* pNLExclusionList;          // Pointer to list of nonbond exclusions
  uint2* pNLExclusionStartCount;           // Pointer to per-atom exclusions
  unsigned int NLExclusions;               // Total number of exclusions
  unsigned int NLAtoms;                    // Total number of neighbor list atoms
  NLRecord* pNLRecord;                     // Pointer to neighbor list records
  NLEntry* pNLEntry;                       // Active neighbor list
  unsigned int* pNLEntries;                // Number of entries in current neighbor list
  unsigned int* pNLOffset;                 // Pointer to atom/exclusion offsets
  unsigned int* pFrcBlkCounters;           // Positions in building/traversing neighbor lists:
                                           //
                                           //   Idx    Description
                                           //    0   Hilbert Space neighbor list building,
                                           //        Hilbert Space NB tiles
                                           //    1   Bond work unit progress
                                           //    2   Sub-image expansion block progress
                                           //    3   Neutral Territory mapping progress, tracks
                                           //        the list of home pencils and each stencil
                                           //        being split into NT regions
                                           //    4   Counter for Neutral Territory regions as
                                           //        they are written out.  Stores the total
                                           //        number of NT regions by the time the
                                           //        pair list is being written and evaluated.
                                           //    5   Counter for traversing the NT region list
                                           //        during pair list construction
                                           //    6   Tracks pair list storage for singletons
                                           //        during pair list construction
  unsigned int NLBuildWarps;               // Number of warps in neighbor list build
  unsigned int NLNonbondEnergyWarps;       // Number of warps in nonbond energy calculation
  unsigned int NLNonbondForcesWarps;       // Number of warps in nonbond forces calculation
  unsigned int NLCellBuffers;              // Number of nonbond cell buffers
  unsigned int NLRecords;                  // Number of neighbor list records
  unsigned int NLMaxEntries;               // Maximum allow entries in neighbor list
                                           //   (~5x expected)
  unsigned int NLXEntryWidth;              // Number of x atoms per NL entry
  unsigned int NLYDivisor;                 // Number of y subdivisions per neighor list record
  unsigned int NLYStride;                  // Neighbor List build vertical atom stride
                                           //   (NLYDivisor * NLAtomsPerWarp)
  unsigned int NLEntryTypes;               // Neighbor List entry types (used to determine
                                           //   output buffers offset)
  unsigned int NLHomeCellBuffer;           // Register atom buffer
  unsigned int NLAtomsPerWarp;             // Number of atoms to process in each warp's
                                           //   registers
  unsigned int NLAtomsPerWarpBits;         // Number of bits in atoms to toprocess in each
                                           //   warp's registers
  unsigned int NLAtomsPerWarpBitsMask;     // NLAtomsPerWarp - 1
  unsigned int NLAtomsPerWarpMask;         // First AtomsPerWarp worth of bits set to 1
  unsigned int NLOffsetPerWarp;            // Number of entries needed to process 1 warp
                                           //   iteration's worth of nonbond forces
  unsigned int NLExclusionBufferSize;      // Exclusion buffer size for L1 cache-capable GPUs
  unsigned int forceBufferStride;          // Total force buffer stride for clearing forces 
  unsigned int* pBNLExclusionBuffer;       // Per-warp GMEM exclusion buffer

  // LJ data
  PMEFloat2* pLJTerm;             // Lennard-Jones terms
  int LJTerms;                    // Number of Lennard-Jones terms
  int LJTypes;                    // Number of Lennard-Jones types
  int LJOffset;                   // Offset to SCTI LJ Terms
  int nR4sources;                 // The number of particles projecting 1/r4 potentials.  Also
                                  //   serves as a flag to activate Lennard-Jones 12-6-4.
  int2* pR4sources;               // Atom IDs and Lennard-Jones types of the 1/r4 potential
                                  //   sources, as indexed in the master topology

  // GB Data
  PMEFloat2* pNeckMaxValPos;      // GB Born Radii and energy correction data
  PMEFloat* pgb_alpha;            // Pointer to per-atom GB alpha
  PMEFloat* pgb_beta;             // Pointer to per-atom GB beta
  PMEFloat* pgb_gamma;            // Pointer to per-atom GB gamma

  // GBSA3 data pwsasa
  PMEFloat* pgbsa_sigma;          // Pointer to per-atom GBSA sigma
  PMEFloat* pgbsa_epsilon;        // Pointer to per-atom GBSA epsilon
  PMEFloat* pgbsa_radius;         // Pointer to per-atom GBSA radius
  PMEFloat* pgbsa_maxsasa;        // Pointer to per-atom GBSA maxSASA

  // PME data
#ifdef use_DPFP
  long long int* plliXYZ_q;       // Input PME charge grid
#else
  int* plliXYZ_q;                 // Input PME charge grid (for single-precision calculations,
                                  //   the charges will still be stored to greater accuracy
                                  //   than the B-spline coefficients and the positions would
                                  //   allow us to calculate)
#endif
  PMEFloat* pXYZ_q;               // PME charge grid/buffer
  PMEComplex* pXYZ_qt;            // FFTed PME charge grid
  PMEFloat* pPrefac1;             // PME nfft1 pre-factors
  PMEFloat* pPrefac2;             // PME nfft2 pre-factors
  PMEFloat* pPrefac3;             // PME nfft3 pre-factors
  PMEFloat* pFractX;              // PME unit cell fractional x coordinates
  PMEFloat* pFractY;              // PME unit cell fractional y coordinates
  PMEFloat* pFractZ;              // PME unit cell fractional z coordinates

  // NTP molecule data
  int  soluteMolecules;           // Total solute molecules
  int  soluteMoleculeStride;      // Total solute molecule stride
  int  soluteAtoms;               // Total solute atoms
  int  soluteAtomsOffset;         // Used for remapping neighbor list data
  int  solventMolecules;          // Total solvent molecules
  int  solventMoleculeStride;     // Total solvent molecules padded to warp width
  int* pSoluteAtomMoleculeID;     // List of solute molecule IDs
  int* pSoluteAtomID;             // List of solute atom IDs
  PMEDouble* pSoluteAtomMass;     // Solute atom masses
  PMEDouble* pSoluteCOMX;         // X Last centers of mass for each solute molecule
  PMEDouble* pSoluteCOMY;         // Y Last centers of mass for each solute molecule
  PMEDouble* pSoluteCOMZ;         // Z Last centers of mass for each solute molecule
  PMEDouble* pSoluteDeltaCOMX;    // X change in center of mass for each solute molecule
  PMEDouble* pSoluteDeltaCOMY;    // Y change in center of mass for each solute molecule
  PMEDouble* pSoluteDeltaCOMZ;    // Z change in center of mass for each solute molecule
  PMEUllInt* pSoluteUllCOMX;      // X Current center of mass for each solute molecule
  PMEUllInt* pSoluteUllCOMY;      // Y Current center of mass for each solute molecule
  PMEUllInt* pSoluteUllCOMZ;      // Z Current center of mass for each solute molecule
  PMEUllInt* pSoluteUllEKCOMX;    // Pointer to x component of COM Kinetic energy buffer
  PMEUllInt* pSoluteUllEKCOMY;    // Pointer to x component of COM Kinetic energy buffer
  PMEUllInt* pSoluteUllEKCOMZ;    // Pointer to x component of COM Kinetic energy buffer
  PMEDouble* pSoluteInvMass;      // Total Inverse mass for each solute molecule
  int4* pSolventAtomID;           // List of solvent molecules/ions of 4 or fewer atoms
  PMEDouble* pSolventAtomMass1;   // First solvent atom mass
  PMEDouble* pSolventAtomMass2;   // Second solvent atom mass
  PMEDouble* pSolventAtomMass3;   // Third solvent atom mass
  PMEDouble* pSolventAtomMass4;   // Fourth solvent atom mass
  PMEDouble* pSolventCOMX;        // X Last centers of mass for each solvent molecule
  PMEDouble* pSolventCOMY;        // Y Last centers of mass for each solvent molecule
  PMEDouble* pSolventCOMZ;        // Z Last centers of mass for each solvent molecule
  PMEDouble* pSolventInvMass;     // Total inverse mass for eache solvent molecule
  int* pAFEMolType;               // Molecule type for COM calculation for alchemical free
                                  //   energy
  int* pAFEMolPartner;            // Partner molecule ID for COM calculations for alchemical
                                  //   free energy

  // NTP molecule constraint data
  PMEDouble* pConstraintAtomX;    // Pre-centered constraint atom x
  PMEDouble* pConstraintAtomY;    // Pre-centered constraint atom y
  PMEDouble* pConstraintAtomZ;    // Pre-centered constraint atom z
  PMEDouble* pConstraintCOMX;     // Original center of mass X, Y, and Z fractional
  PMEDouble* pConstraintCOMY;     //   coordinates for constraint atoms
  PMEDouble* pConstraintCOMZ;     //

  // Energy and Virial Buffers
  unsigned long long int* pEnergyBuffer;     // Generic energy buffer pointer
  unsigned long long int* pEELT;             // Pointer to electrostatic energy
  unsigned long long int* pEVDW;             // Pointer to vdw energy
  unsigned long long int* pEGB;              // Pointer to Generalized Born energy
  unsigned long long int* pEBond;            // Pointer to bond energy
  unsigned long long int* pEAngle;           // Pointer to bond angle energy
  unsigned long long int* pEDihedral;        // Pointer to dihedral energy
  unsigned long long int* pEEL14;            // Pointer to 1-4 electrostatic energy
  unsigned long long int* pENB14;            // Pointer to 1-4 vdw energy
  unsigned long long int* pEConstraint;      // Pointer to restraint energy
  unsigned long long int* pEER;              // Pointer to PME reciprocal space
                                             //   electrostatic energy
  unsigned long long int* pEED;              // Pointer to PME direct space
                                             //   electrostatic energy
  unsigned long long int* pEAngle_UB;        // Pointer to CHARMM Urey Bradley energy
  unsigned long long int* pEImp;             // Pointer to CHARMM improper dihedral energy
  unsigned long long int* pECmap;            // Pointer to CHARMM cmap energy
  unsigned long long int* pENMRDistance;     // Pointer to NMR distance energy
  unsigned long long int* pENMRCOMDistance;  // Pointer to NMR COM distance energy
  unsigned long long int* pENMRr6avDistance; // Pointer to NMR r6av distance energy
  unsigned long long int* pENMRAngle;        // Pointer to NMR angle energy
  unsigned long long int* pENMRCOMAngle;     // Pointer to NMR COM angle energy
  unsigned long long int* pENMRTorsion;      // Pointer to NMR torsion energy
  unsigned long long int* pENMRCOMTorsion;  // Pointer to NMR COM torsion energy
  unsigned long long int* pESurf;            // Pointer to GBSA surface energy pwsasa

  unsigned long long int* pEEField;          // Pointer to Electric Field energy
  unsigned long long int* pEEDE;             // Pointer to PME direct space electrostatic energy extended accumulator
  unsigned long long int* pEVDWE;            // Pointer to VDW energy extended accumulator 
  unsigned long long int* pEELTE;            // Pointer to GB electrostatic energy extended accumulator   
  unsigned long long int* pAFEBuffer;        // Pointer to alchemical free energy specific
                                             //   terms
  unsigned long long int* pDVDL;             // Pointer to TI deriv. of pot ene w.r.t. lambda
  unsigned long long int* pSCBondR1;         // Pointer to bond energies unique to region 1
  unsigned long long int* pSCBondR2;         // Pointer to bond energies unique to region 2
  unsigned long long int* pSCBondAngleR1;    // Pointer to angle energies unique to region 1
  unsigned long long int* pSCBondAngleR2;    // Pointer to angle energies unique to region 2
  unsigned long long int* pSCDihedralR1;     // Pointer to dihedral energies unique to region 1
  unsigned long long int* pSCDihedralR2;     // Pointer to dihedral energies unique to region 2
  unsigned long long int* pESCNMRDistanceR1;     // Pointer to softcore NMR distance restraint
  unsigned long long int* pESCNMRDistanceR2;     //   energies unique to region 1 or 2
  unsigned long long int* pESCNMRCOMDistanceR1;  // Pointer to softcore NMR COM restraint
  unsigned long long int* pESCNMRCOMDistanceR2;  //   energies unique to region 1 or 2
  unsigned long long int* pESCNMRr6avDistanceR1; // Pointer to softcore NMR r6av restraint
  unsigned long long int* pESCNMRr6avDistanceR2; //    energies unique to region 1 or 2
  unsigned long long int* pESCNMRAngleR1;        // Pointer to softcore NMR angle restraint
  unsigned long long int* pESCNMRAngleR2;        //   energies unique to region 1 or 2
  unsigned long long int* pESCNMRTorsionR1;      // Pointer to softcore NMR torsion restraint
  unsigned long long int* pESCNMRTorsionR2;      //  energies unique to region 1 or 2
  unsigned long long int* pSCVDWDirR1;       // Pointers to vdw direct energies unique to
  unsigned long long int* pSCVDWDirR2;       //      region 1 or 2
  unsigned long long int* pSCEELDirR1;       // Pointer to eel direct energies unique to
  unsigned long long int* pSCEELDirR2;       //      region 1 or 2
  unsigned long long int* pSCVDW14R1;        // Pointer to vdw 14 energies unique to region 1
  unsigned long long int* pSCVDW14R2;        // Pointer to vdw 14 energies unique to region 2
  unsigned long long int* pSCEEL14R1;        // Pointer to eel 14 energies unique to region 1
  unsigned long long int* pSCEEL14R2;        // Pointer to eel 14 energies unique to region 2
  unsigned long long int* pSCVDWDerR1;       // Pointer to vdw derivatives energies unique to
  unsigned long long int* pSCVDWDerR2;       //   region 1 or 2
  unsigned long long int* pSCEELDerR1;       // Pointer to eel derivatives energies unique to
  unsigned long long int* pSCEELDerR2;       //   region 1 or 2
  unsigned long long int* pVirial;           // Pointer to PME virial
  unsigned long long int* pVirial_11;        // Pointer to PME virial component
  unsigned long long int* pVirial_22;        // Pointer to PME virial component
  unsigned long long int* pVirial_33;        // Pointer to PME virial component
  unsigned long long int* pVirial_11E;       // Pointer to PME virial component extended precision
  unsigned long long int* pVirial_22E;       // Pointer to PME virial component extended precision
  unsigned long long int* pVirial_33E;       // Pointer to PME virial component extended precision
  unsigned long long int* pEKCOMX;           // Pointers to the X, Y, and Z components of
  unsigned long long int* pEKCOMY;           //   PME center of mass kinetic energy
  unsigned long long int* pEKCOMZ;           //
  unsigned int EnergyTerms;                  // Total energy terms
  unsigned int AFETerms;                     // Total alchemical free energy terms

  // Kinetic Energy Buffers
  KineticEnergy* pKineticEnergy;        // Pointer to per-block Kinetic Energy entries
  AFEKineticEnergy* pAFEKineticEnergy;  // Pointer to per-block extra terms for Alchemical
                                        //   Free Energy kinetic energy entries

  // Random Number stuff
  unsigned int randomSteps;             // Number of steps between RNG calls
  unsigned int randomNumbers;           // Number of randoms to generate per atom per RNG call
  double* pRandom;                      // Pointer to overall RNG buffer
  double* pRandomX;                     // Pointer to x random numbers
  double* pRandomY;                     // Pointer to y random numbers
  double* pRandomZ;                     // Pointer to z random numbers

  // Extra points stuff
  int EPs;                          // Counts of extra points.  These refer to the total number
  int EP11s;                        //   of extra points and then the numbers of EP type I with
  int EP12s;                        //   one extra point on the parent atom, two extra points
  int EP21s;                        //   on the parent atom, then type II with one extra point
  int EP22s;                        //   on the parent atom, or two.  Finally, the number of
  int EPCustomCount;                //   custom EP frames (all have one point) is given.
  int EP11Offset;                   //
  int EP12Offset;                   // Offsets in the GPU data storage arrays for EP frame
  int EP21Offset;                   //   data.  These numbers are cumulative sums of the above,
  int EP22Offset;                   //   padded forward by GRID.
  int EPCustomOffset;               //
  int4*   pExtraPoint11Frame;       //
  int*    pExtraPoint11Index;       //
  double* pExtraPoint11X;           //
  double* pExtraPoint11Y;           //
  double* pExtraPoint11Z;           //
  int4*   pExtraPoint12Frame;       //
  int*    pExtraPoint12Index;       //
  double* pExtraPoint12X;           //
  double* pExtraPoint12Y;           //
  double* pExtraPoint12Z;           //
  int4*   pExtraPoint21Frame;       //
  int2*   pExtraPoint21Index;       //
  double* pExtraPoint21X1;          //
  double* pExtraPoint21Y1;          //
  double* pExtraPoint21Z1;          //
  double* pExtraPoint21X2;          //
  double* pExtraPoint21Y2;          //
  double* pExtraPoint21Z2;          //
  int4*   pExtraPoint22Frame;       //
  int2*   pExtraPoint22Index;       //
  double* pExtraPoint22X1;          //
  double* pExtraPoint22Y1;          //
  double* pExtraPoint22Z1;          //
  double* pExtraPoint22X2;          //
  double* pExtraPoint22Y2;          //
  double* pExtraPoint22Z2;          //
  int4*   pExtraPointCustomFrame;   //
  int*    pExtraPointCustomIndex;   //
  double* pExtraPointCustomD1;      //
  double* pExtraPointCustomD2;      //
  double* pExtraPointCustomD3;      //
  int4*   pImageExtraPoint11Frame;  //
  int*    pImageExtraPoint11Index;  //
  int4*   pImageExtraPoint12Frame;  //
  int*    pImageExtraPoint12Index;  //
  int4*   pImageExtraPoint21Frame;  //
  int2*   pImageExtraPoint21Index;  //
  int4*   pImageExtraPoint22Frame;  //
  int2*   pImageExtraPoint22Index;  //
  int4*   pImageEPCustomFrame;      //
  int*    pImageEPCustomIndex;      //

  // Shake constraint stuff
  unsigned int shakeConstraints;     // Traditional SHAKE constraints
  unsigned int shakeOffset;          // Offset to end of traditional SHAKE constraints
  unsigned int fastShakeConstraints; // Fast SHAKE constraints (H2O molecules)
  unsigned int fastShakeOffset;      // Offset to end of fast SHAKE constraints
  unsigned int slowShakeConstraints; // XH4 (really slow) SHAKE constraints
  unsigned int slowShakeOffset;      // Offset to end of slow SHAKE constraints
  unsigned int slowTIShakeConstraints; // XH4 (really slow) SHAKE constraints
  unsigned int slowTIShakeOffset;      // Offset to end of slow SHAKE constraints
  int tishake;                       // TI Shake toggle to determine which shake variant to do for TI
  int4* pShakeID;                    // SHAKE central atom plus up to 3 hydrogens
  double2* pShakeParm;               // SHAKE central atom mass and equilibrium bond length
  double* pShakeInvMassH;            // SHAKE HMR inverse hydrogen mass
  int4* pFastShakeID;                // H2O oxygen plus two hydrogens atom ID
  int* pSlowShakeID1;                // Central atom of XH4 Shake constraint
  int4* pSlowShakeID2;               // XH4 SHAKE constraint hydrogens
  double2* pSlowShakeParm;           // XH4 SHAKE central atom mass and equilibrium bond length
  double* pSlowShakeInvMassH;        // XH4 SHAKE HMR inverse hydrogen mass
  int* pSlowTIShakeID1;              // Central atom of XH8 Shake constraint
  int4* pSlowTIShakeID2;             // XH8 SHAKE constraint hydrogens
  int4* pSlowTIShakeID3;             // XH8 SHAKE constraint hydrogens
  double2* pSlowTIShakeParm;         // XH8 SHAKE central atom mass and equilibrium bond length
  double* pSlowTIShakeInvMassH;      // XH8 SHAKE HMR inverse hydrogen mass
  double2* pShakeParm2;              // central atom mass and equilibrium bond length
  double2* pSlowShakeParm2;          // XH4 SHAKE central atom mass and equilibrium bond length
  double2* pSlowTIShakeParm2;        // XH8 SHAKE central atom mass and equilibrium bond length
  aligned_double ra;                 // Fast SHAKE parameter
  aligned_double ra_inv;             // Fast SHAKE parameter
  aligned_double rb;                 // Fast SHAKE parameter
  aligned_double rc;                 // Fast SHAKE parameter
  aligned_double rc2;                // Fast SHAKE parameter
  aligned_double hhhh;               // Fast SHAKE parameter
  aligned_double wo_div_wohh;        // Fast SHAKE parameter
  aligned_double wh_div_wohh;        // Fast SHAKE parameter
  aligned_double wo;                 // Fast RATTLE parameter added by zhf used for
  aligned_double wh;                 // Fast RATTLE parameter middle-scheme RATTLE

  // Bonded interaction stuff
  int bonds;                // Total number of bonds
  int bondOffset;           // Offset to end of bondss
  int bondAngles;           // Total number of bond angles
  int bondAngleOffset;      // Offset to end of bond angles
  int dihedrals;            // Total number of dihedrals
  int dihedralOffset;       // Offset to end of dihedrals
  int nb14s;                // Total number of 1-4 nonbond interactions
  int nb14Offset;           // Offset to end of 1-4 nonbond interactions
  int constraints;          // Number of positional constraints
  int constraintOffset;     // Offset to end of positional constraints
  int UBAngles;             // Total number of Urey Bradley angles
  int UBAngleOffset;        // Offset to end of Urey Bradley Angles
  int impDihedrals;         // Total number of (CHARMM) improper dihedrals
  int impDihedralOffset;    // Offset to end of (CHARMM) improper dihedrals
  int cmaps;                // Total number of cmap terms
  int cmapOffset;           // Offset to end of cmap terms
  int qqxcs;                // Number of charge-charge interaction exclusions (these are merely
                            //   subtracted from the current forces and energies, unlike
                            //   Lennard-Jones interactions which we prevent from being
                            //   calculated)
  PMEDouble2* pBond;        // Bond rk, req
  int2* pBondID;            // Bond i, j
  int2* pQQxcID;            // Electrostatic exclusion i, j
  PMEDouble2* pBondAngle;   // Bond Angle Kt, teq
  int2* pBondAngleID1;      // Bond Angle i, j
  int* pBondAngleID2;       // Bond Angle k;
  PMEFloat2* pDihedral1;    // Dihedral Ipn, pn
  PMEFloat2* pDihedral2;    // Dihedral pk, gamc
  PMEFloat* pDihedral3;     // Dihedral gams
  int4* pDihedralID1;       // Dihedral i, j, k, l
  PMEDouble2* pNb141;       // 1-4 nonbond scee, scnb0
  PMEDouble2* pNb142;       // 1-4 nonbond cn1, cn2
  int2* pNb14ID;            // 1-4 nonbond i, j
  PMEDouble2* pConstraint1; // Constraint weight and xc
  PMEDouble2* pConstraint2; // Constraint yc and zc
  int* pConstraintID;       // Atom constraint ID

  // CHARMM interaction stuff
  PMEDouble2* pUBAngle;     // Urey Bradley Angle rk, r0
  int2* pUBAngleID;         // Urey Bradley Angle i, j
  PMEDouble2* pImpDihedral; // Improper Dihedral pk, phase
  int4* pImpDihedralID1;    // Improper Dihedral i, j, k, l
  int4* pCmapID1;           // Cmap i, j, k, l
  int* pCmapID2;            // Cmap m
  int* pCmapType;           // Cmap type

  // NMR refinement stuff
  bool NMRR6av;                        // Whether or not R6av restraints are on or not
  int NMRDistances;                    // Number of NMR distance constraints
  int NMRCOMDistances;                 // Number of NMR COM distance constraints
  int NMRr6avDistances;                // Number of NMR r6av distance constraints
  int NMRMaxgrp;                       // Number of NMR distance COM groupings
  int NMRDistanceOffset;               // Offset to end of distance constraints
  int NMRCOMDistanceOffset;            // Offset to end of COM distance constraints
  int NMRr6avDistanceOffset;           // Offset to end of r6av distance constraints
  int NMRAngles;                       // Number of NMR angle constraints
  int NMRCOMAngles;                 // number of nmr com angle constraints
  int NMRCOMAngleOffset;            // offset to end of com angle constraints
  int NMRAngleOffset;                  // Offset to end of angle constraints
  int NMRTorsions;                     // Number of NMR torsion constraints
  int NMRTorsionOffset;                // Offset to end of torsion constraints
  int NMRCOMTorsions;                 // Number of NMR COM torsion constraints
  int NMRCOMTorsionOffset;            // Offset to end of COM torsion constraints

  bool bJar;                           // Determines whether Jarzynski MD is active
  double drjar;                        // Jar increment
  double* pNMRJarData;                 // Jarzynski accumulated work data
  int2* pNMRDistanceID;                // NMR distance i, j
  PMEDouble2* pNMRDistanceR1R2;        // NMR distance computed r1, r2
  PMEDouble2* pNMRDistanceR3R4;        // NMR distance computed r3, r4
  PMEDouble2* pNMRDistanceK2K3;        // NMR distance computed k2, k3
  PMEDouble* pNMRDistanceAve;          // NMR distance restraint linear and
                                       //   exponential averaged value
  PMEDouble2* pNMRDistanceTgtVal;      // NMR distance target and value for current step
  int2* pNMRDistanceStep;              // NMR distance first and last step for
                                       //   application of restraint
  int* pNMRDistanceInc;                // NMR distance increment for step weighting
  PMEDouble2* pNMRDistanceR1R2Slp;     // NMR distance r1, r2 slope
  PMEDouble2* pNMRDistanceR3R4Slp;     // NMR distance r3, r4 slope
  PMEDouble2* pNMRDistanceK2K3Slp;     // NMR distance k2, k3 slope
  PMEDouble2* pNMRDistanceR1R2Int;     // NMR distance r1, r2 intercept
  PMEDouble2* pNMRDistanceR3R4Int;     // NMR distance r3, r4 intercept
  PMEDouble2* pNMRDistanceK2K3Int;     // NMR distance k2, k3 intercept

  int2* pNMRCOMDistanceID;             // NMR COM distance i, j
  int2* pNMRCOMDistanceCOM;            // NMR COM distance COM ranges from .x to .y
  int2* pNMRCOMDistanceCOMGrp;         // NMR COM distance COM indexing for no. of
                                       //   atom ranges in a COM group from .x to .y
  PMEDouble2* pNMRCOMDistanceR1R2;     // NMR COM distance computed r1, r2
  PMEDouble2* pNMRCOMDistanceR3R4;     // NMR COM distance computed r3, r4
  PMEDouble2* pNMRCOMDistanceK2K3;     // NMR COM distance computed k2, k3
  PMEDouble* pNMRCOMDistanceAve;       // NMR COM distance restraint linear and
                                       //   exponential averaged value
  PMEDouble2* pNMRCOMDistanceTgtVal;   // NMR COM distance target and value for current step
  int2* pNMRCOMDistanceStep;           // NMR COM distance first and last step for
                                       //   application of restraint
  int* pNMRCOMDistanceInc;             // NMR COM distance increment for step weighting
  PMEDouble2* pNMRCOMDistanceR1R2Slp;  // NMR COM distance r1, r2 slope
  PMEDouble2* pNMRCOMDistanceR3R4Slp;  // NMR COM distance r3, r4 slope
  PMEDouble2* pNMRCOMDistanceK2K3Slp;  // NMR COM distance k2, k3 slope
  PMEDouble2* pNMRCOMDistanceR1R2Int;  // NMR COM distance r1, r2 intercept
  PMEDouble2* pNMRCOMDistanceR3R4Int;  // NMR COM distance r3, r4 intercept
  PMEDouble2* pNMRCOMDistanceK2K3Int;  // NMR COM distance k2, k3 intercept
  int* pNMRCOMDistanceWeights;         // NMR COM distance weights
  PMEDouble* pNMRCOMDistanceXYZ;       // NMR X,Y,Z COM Distance Components
  int2* pNMRr6avDistanceID;            // NMR r6av distance i, j
  int2* pNMRr6avDistancer6av;          // NMR r6av distance r6av ranges from .x to .y
  int2* pNMRr6avDistancer6avGrp;       // NMR r6av distance r6av indexing for no. of
                                       //   atom ranges in a r6av group from .x to .y
  PMEDouble2* pNMRr6avDistanceR1R2;    // NMR r6av distance computed r1, r2
  PMEDouble2* pNMRr6avDistanceR3R4;    // NMR r6av distance computed r3, r4
  PMEDouble2* pNMRr6avDistanceK2K3;    // NMR r6av distance computed k2, k3
  PMEDouble* pNMRr6avDistanceAve;      // NMR r6av distance restraint linear and
                                       //   exponential averaged value
  PMEDouble2* pNMRr6avDistanceTgtVal;  // NMR r6av distance target and value for current step
  int2* pNMRr6avDistanceStep;          // NMR r6av distance first and last step for
                                       //   application of restraint
  int* pNMRr6avDistanceInc;            // NMR r6av distance increment for step weighting
  PMEDouble2* pNMRr6avDistanceR1R2Slp; // NMR r6av distance r1, r2 slope
  PMEDouble2* pNMRr6avDistanceR3R4Slp; // NMR r6av distance r3, r4 slope
  PMEDouble2* pNMRr6avDistanceK2K3Slp; // NMR r6av distance k2, k3 slope
  PMEDouble2* pNMRr6avDistanceR1R2Int; // NMR r6av distance r1, r2 intercept
  PMEDouble2* pNMRr6avDistanceR3R4Int; // NMR r6av distance r3, r4 intercept
  PMEDouble2* pNMRr6avDistanceK2K3Int; // NMR r6av distance k2, k3 intercept
  int* pNMRr6avDistanceWeights;         // NMR r6av  distance weights
  PMEDouble* pNMRr6avDistanceXYZ;       // NMR X,Y,Z r6av Distance Components
  int2* pNMRAngleID1;                  // NMR angle i, j
  int* pNMRAngleID2;                   // NMR angle k
  PMEDouble2* pNMRAngleR1R2;           // NMR angle computed r1, r2
  PMEDouble2* pNMRAngleR3R4;           // NMR angle computed r3, r4
  PMEDouble2* pNMRAngleK2K3;           // NMR angle computed k2, k3
  PMEDouble* pNMRAngleAve;             // NMR angle restraint linear and exponential
                                       //   averaged value
  PMEDouble2* pNMRAngleTgtVal;         // NMR angle target and value for current step
  int2* pNMRAngleStep;                 // NMR angle first and last step for application
                                       //   of restraint
  int* pNMRAngleInc;                   // NMR angle increment for step weighting
  PMEDouble2* pNMRAngleR1R2Slp;        // NMR angle r1, r2 slope
  PMEDouble2* pNMRAngleR3R4Slp;        // NMR angle r3, r4 slope
  PMEDouble2* pNMRAngleK2K3Slp;        // NMR angle k2, k3 slope
  PMEDouble2* pNMRAngleR1R2Int;        // NMR angle r1, r2 intercept
  PMEDouble2* pNMRAngleR3R4Int;        // NMR angle r3, r4 intercept
  PMEDouble2* pNMRAngleK2K3Int;        // NMR angle k2, k3 intercept

  int2* pNMRCOMAngleID1;               // NMR COM angle i, j
  int*  pNMRCOMAngleID2;               // NMR COM angle k
  int2* pNMRCOMAngleCOM;               // NMR COM angle COM ranges from .x to .y
  int2* pNMRCOMAngleCOMGrp;            // NMR COM angle COM indexing for no. of
                                       //   atom ranges in a COM group from .x to .y
  PMEDouble2* pNMRCOMAngleR1R2;        // NMR COM angle computed r1, r2
  PMEDouble2* pNMRCOMAngleR3R4;        // NMR COM angle computed r3, r4
  PMEDouble2* pNMRCOMAngleK2K3;        // NMR COM angle computed k2, k3
  PMEDouble* pNMRCOMAngleAve;          // NMR COM angle restraint linear and
                                       //   exponential averaged value
  PMEDouble2* pNMRCOMAngleTgtVal;   // NMR COM angle target and value for current step
  int2* pNMRCOMAngleStep;           // NMR COM angle first and last step for
                                       //   application of restraint
  int* pNMRCOMAngleInc;             // NMR COM angle increment for step weighting
  PMEDouble2* pNMRCOMAngleR1R2Slp;  // NMR COM angle r1, r2 slope
  PMEDouble2* pNMRCOMAngleR3R4Slp;  // NMR COM angle r3, r4 slope
  PMEDouble2* pNMRCOMAngleK2K3Slp;  // NMR COM angle k2, k3 slope
  PMEDouble2* pNMRCOMAngleR1R2Int;  // NMR COM angle r1, r2 intercept
  PMEDouble2* pNMRCOMAngleR3R4Int;  // NMR COM angle r3, r4 intercept
  PMEDouble2* pNMRCOMAngleK2K3Int;  // NMR COM angle k2, k3 intercept


  int4* pNMRTorsionID1;                // NMR torsion i, j, k, l
  PMEDouble2* pNMRTorsionR1R2;         // NMR torsion computed r1, r2
  PMEDouble2* pNMRTorsionR3R4;         // NMR torsion computed r3, r4
  PMEDouble2* pNMRTorsionK2K3;         // NMR torsion computed k2, k3
  PMEDouble* pNMRTorsionAve1;          // NMR torsion restraint linear and exponential
                                       //   averaged value
  PMEDouble* pNMRTorsionAve2;          // NMR torsion restraint linear and exponential
                                       //   averaged value
  PMEDouble2* pNMRTorsionTgtVal;       // NMR torsion target and value for current step
  int2* pNMRTorsionStep;               // NMR torsion first and last step for application
                                       //   of restraint
  int* pNMRTorsionInc;                 // NMR torsion increment for step weighting
  PMEDouble2* pNMRTorsionR1R2Slp;      // NMR torsion r1, r2 slope
  PMEDouble2* pNMRTorsionR3R4Slp;      // NMR torsion r3, r4 slope
  PMEDouble2* pNMRTorsionK2K3Slp;      // NMR torsion k2, k3 slope
  PMEDouble* pNMRTorsionK4Slp;         // NMR torsion k4 slope
  PMEDouble2* pNMRTorsionR1R2Int;      // NMR torsion r1, r2 intercept
  PMEDouble2* pNMRTorsionR3R4Int;      // NMR torsion r3, r4 intercept
  PMEDouble2* pNMRTorsionK2K3Int;      // NMR torsion k2, k3 intercept

  int4* pNMRCOMTorsionID1;             // NMR COM torsion i, j,k,l
  int2* pNMRCOMTorsionCOM;            // NMR COM torsion COM ranges from .x to .y
  int2* pNMRCOMTorsionCOMGrp;         // NMR COM torsion COM indexing for no. of
                                       //   atom ranges in a COM group from .x to .y
  PMEDouble2* pNMRCOMTorsionR1R2;     // NMR COM torsion computed r1, r2
  PMEDouble2* pNMRCOMTorsionR3R4;     // NMR COM torsion computed r3, r4
  PMEDouble2* pNMRCOMTorsionK2K3;     // NMR COM torsion computed k2, k3
  PMEDouble* pNMRCOMTorsionAve;       // NMR COM torsion restraint linear and
                                       //   exponential averaged value
  PMEDouble2* pNMRCOMTorsionTgtVal;   // NMR COM torsion target and value for current step
  int2* pNMRCOMTorsionStep;           // NMR COM torsion first and last step for
                                       //   application of restraint
  int* pNMRCOMTorsionInc;             // NMR COM torsion increment for step weighting
  PMEDouble2* pNMRCOMTorsionR1R2Slp;  // NMR COM torsion r1, r2 slope
  PMEDouble2* pNMRCOMTorsionR3R4Slp;  // NMR COM torsion r3, r4 slope
  PMEDouble2* pNMRCOMTorsionK2K3Slp;  // NMR COM torsion k2, k3 slope
  PMEDouble2* pNMRCOMTorsionR1R2Int;  // NMR COM torsion r1, r2 intercept
  PMEDouble2* pNMRCOMTorsionR3R4Int;  // NMR COM torsion r3, r4 intercept
  PMEDouble2* pNMRCOMTorsionK2K3Int;  // NMR COM torsion k2, k3 intercept



  // Cmap lookup table
  int cmapTermStride;                  // Per term stride
  int cmapRowStride;                   // Per row of terms stride
  PMEFloat4* pCmapEnergy;              // Pointer to Cmap LUT data (E, dPhi, dPsi, dPhi_dPsi)

  // Accumulation buffers
  PMEAccumulator* pForceAccumulator;         // Bare fixed point force accumulator
  PMEAccumulator* pForceXAccumulator;        // Bare fixed point x force accumulator
  PMEAccumulator* pForceYAccumulator;        // Bare fixed point y force accumulator
  PMEAccumulator* pForceZAccumulator;        // Bare fixed point z force accumulator
  PMEAccumulator* pNBForceAccumulator;       // Bare fixed point nonbond force accumulator
  PMEAccumulator* pNBForceXAccumulator;      // Fixed point nonbond x force accumulator
  PMEAccumulator* pNBForceYAccumulator;      // Fixed point nonbond y force accumulator
  PMEAccumulator* pNBForceZAccumulator;      // Fixed point nonbond z force accumulator
  PMEAccumulator* pIntForceAccumulator;      // Integer force for minimizations
  PMEAccumulator* pIntForceXAccumulator;     // Integer force x component for minimizations
  PMEAccumulator* pIntForceYAccumulator;     // Integer force y component for minimizations
  PMEAccumulator* pIntForceZAccumulator;     // Integer force z component for minimizations    
  PMEAccumulator* pReffAccumulator;          // Effective Born Radius buffer
  PMEAccumulator* pSumdeijdaAccumulator;     // Atom Sumdeijda buffer
  PMEAccumulator* pBondedForceAccumulator;   // Bare fixed point bonded force accumulator
  PMEAccumulator* pBondedForceXAccumulator;  // Bare fixed point bonded x force accumulator
  PMEAccumulator* pBondedForceYAccumulator;  // Bare fixed point bonded y force accumulator
  PMEAccumulator* pBondedForceZAccumulator;  // Bare fixed point bonded z force accumulator
  int stride;                                // Atom quantity stride
  int stride2;                               // Atom quantity 2x stride
  int stride3;                               // Atom quantity 3x stride
  int stride4;                               // Atom quantity 4x stride
  int stride5;                               // Atom quantity 5x stride
  int imageStride;                           // Neighor list index stride (1K-aligned) for
                                             // Duane Merrill's radix sort

  // Kernel call stuff
  unsigned int workUnits;                    // Total work units
  unsigned int excludedWorkUnits;            // Total work units with exclusions
  PMEMask* pExclusion;                       // Exclusion masks
  unsigned int* pWorkUnit;                   // Work unit list
  unsigned int* pGBBRPosition;               // Generalized Born Born Radius workunit position
  unsigned int* pGBNB1Position;              // Generalized Born Nonbond 1 workunit position
  unsigned int* pGBNB2Position;              // Generalized Born Nonbond 2 workunit position
  unsigned int GBTotalWarps[3];              // Total warps in use for nonbond kernels
  unsigned int maxForceBuffers;              // Total Force buffer count
  unsigned int nonbondForceBuffers;          // Nonbond Force buffer count
  PMEFloat cellOffset[NEIGHBOR_CELLS][3];    // Local coordinate system offsets

  // texture object
  cudaTextureObject_t texImageX;
  cudaTextureObject_t texImageX2;
  cudaTextureObject_t texAtomX;
  cudaTextureObject_t texXYZ_q;
  cudaTextureObject_t texOldAtomX;
#if !defined(use_DPFP)
  cudaTextureObject_t texAtomXYSP;
  cudaTextureObject_t texAtomChargeSPLJID;
  cudaTextureObject_t texLJTerm;
  cudaTextureObject_t texAtomZSP;
  cudaTextureObject_t texErfcCoeffsTable;
#endif
  // Data for continuous pHMD
  bool qsetlam;                              // Are the titration variables fixed?
  bool use_forces;                           // Are the forces used for the titration variables?
  int ntitr;                                 // Number of titrating residues
  int iphmd;                                 // Version of phmd we are running
  int ntitratoms;                            // Number of titrating atoms
  int phmd_stride;                           // Number of titrating atoms per thread in phmd
  int phmd_atom_stride;                      // Number of titrating atoms per thread in phmd
  int phmd_padded_ntitr;                     // Padded number of titrating residues
  PMEFloat qmass_phmd;                         // Mass of titration particles
  PMEFloat temp_phmd;                          // Temperature of titration particles
  PMEFloat phbeta;                             // Friction constant in titration integrator
  PMEFloat solvph;                            // Solution pH
  int randomNumbersPHMD;                     // Number of random numbers to generate for 
                                             // pHMD Langevin integrators
  double* pRandomPHMD;                      // Pointer to overall RNG buffer for pHMD
  PMEAccumulator* pdph_accumulator;   // Pointer to titration force accumulators
  PMEAccumulator* pdph_plus_accumulator;   // Pointer to titration force accumulators
  unsigned long long int* pphmd_ene;   // Pointer to energy accumulator
  PMEFloat* pdph_theta;                       // Pointer to theta forces
  PMEFloat* pvph_theta;                       // Pointer to the velocities of the titration coordinates
  PMEFloat* pdudls;                           // Pointer to the dudls in cphmd
  PMEFloat* pvphold;                          // Pointer to the hold velocities of the titration coordinates
  PMEFloat* pph_theta;                        // Pointer to the titration variables
  PMEFloat* pdphold;                          // Pointer to the hold derivatives of the titration coordinates
  int2*   presbounds;                       // Pointer to the resbounds
  bool* pbaseres;                           // Pointer to whether it is a base residue
  int* psp_grp;                         // Pointer to the sp_grp ids
  int* ptitratoms;                            //Pointer to the array of titrating atom indices
  PMEFloat2* pqstate1;                       // Pointer to the qstate1 array
  PMEFloat2* pqstate2;                       // Pointer to the qstate2 array
  PMEFloat2* pImageQstate1;                  // Pointer to the qstate1 array, image
  PMEFloat2* pImageQstate2;                  // Pointer to the qstate2 array, image
  PMEFloat2* pImageQstate12;                  // Pointer to the qstate1 array, image
  PMEFloat2* pImageQstate22;                  // Pointer to the qstate2 array, image
  PMEFloat2* pqstate1_md;                   // Pointer to the qstate1_md array
  PMEFloat2* pqstate2_md;                    // Pointer to the qstate2_md array
  PMEFloat2* pvstate1;                       // Pointer to the vstate1 array
  PMEFloat2* pvstate2;                       // Pointer to the vstate2 array
  PMEFloat2* pImageVstate1;                  // Pointer to the vstate1 array, image
  PMEFloat2* pImageVstate2;                 // Pointer to the vstate2 array, image
  PMEFloat2* pImageVstate12;                  // Pointer to the vstate1 array, image
  PMEFloat2* pImageVstate22;                 // Pointer to the vstate2 array, image
  PMEFloat* pcharge_phmd;                     // Pointer to the charges for phmd
  PMEFloat* pImageCharge_phmd;                // Pointer to the charges for phmd, image
  PMEFloat* pImageCharge_phmd2;               // Pointer to the charges for phmd, image
  PMEFloat* pparmod;                          // Pointer to parmod parameters
  PMEFloat* pbarr;                            // Pointer to barrier heights
  PMEFloat* ppark;                            // Pointer to park parameters
  PMEFloat* ppara;                            // Pointer to para parameters
  PMEFloat* pparb;                            // Pointer to parb parameters
  PMEFloat2* psp_par;                         // Pointer to sp_par parameters
  int padded_numch;                         // Number of threads to launch for each residue
  int* pexcludeatoms;                        // Pointer to excluded atoms
  int* pgrplist;                             // Pointer to grplist
  int* plinkgrplist;                         // Pointer to list of linked atoms
  int* pImageGrplist;                        // Pointer to image grplist
  int* pImageGrplist2;                        // Pointer to image grplist
  int* pImageIndexPHMD;                       //Pointer to image indices
  int* pImageIndexPHMD2;                      //Pointer to image indices
};

#endif // _BASE_CUDA_SIMULATION
