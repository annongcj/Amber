#ifndef BOND_REMAP_STRUCTS
#define BOND_REMAP_STRUCTS

//TL: enforce correct usage of const in different env.(CPU or GPU)
#ifdef GL_CONST
#undef GL_CONST
#endif
#ifdef  __CUDA_ARCH__
#define GL_CONST  static const __constant__
#else
#define GL_CONST  static const
#endif

// Object codes--positive integers from 0 to 255 are acceptable
GL_CONST int BOND_CODE                        =    2;
GL_CONST int ANGL_CODE                        =    3;
GL_CONST int DIHE_CODE                        =    4;
GL_CONST int CMAP_CODE                        =    5;
GL_CONST int QQXC_CODE                        =    6;
GL_CONST int NB14_CODE                        =    7;
GL_CONST int NMR2_CODE                        =    8;
GL_CONST int NMR3_CODE                        =    9;
GL_CONST int NMR4_CODE                        =   10;
GL_CONST int UREY_CODE                        =   11;
GL_CONST int CIMP_CODE                        =   12;
GL_CONST int CNST_CODE                        =   13;

// Force accumulator codes (destination for each unit's results)
GL_CONST int BOND_FORCE_ACCUMULATOR           =    0;
GL_CONST int NB_FORCE_ACCUMULATOR             =    1;

// Indexing into the GPU instruction sets
GL_CONST int ATOM_IMPORT_COUNT_IDX            =    0;
GL_CONST int ATOM_IMPORT_SOURCE_IDX           =    1;
GL_CONST int WARP_INSTRUCTION_COUNT_IDX       =    2;
GL_CONST int FORCE_ACCUMULATOR_IDX            =    3;
GL_CONST int WARP_INSTRUCTION_OFFSET          =    8;
GL_CONST int BOND_WORK_UNIT_BLOCKS_MULTIPLIER =    6;
GL_CONST int BOND_WORK_UNIT_WARPS_PER_BLOCK   =    4;
#ifdef AMBER_PLATFORM_AMD_WARP64
GL_CONST int BOND_WORK_UNIT_THREADS_PER_BLOCK =  256;
GL_CONST int CHARGE_BUFFER_CYCLES             =   12;
#else
GL_CONST int BOND_WORK_UNIT_THREADS_PER_BLOCK =  128;
GL_CONST int CHARGE_BUFFER_CYCLES             =   24;
#endif
GL_CONST int CHARGE_BUFFER_STRIDE             = BOND_WORK_UNIT_THREADS_PER_BLOCK *
                                                CHARGE_BUFFER_CYCLES;


// Indexing into kernel energy accumulators
GL_CONST int BOND_EACC_OFFSET = 0;
GL_CONST int ANGL_EACC_OFFSET = BOND_EACC_OFFSET + BOND_WORK_UNIT_WARPS_PER_BLOCK;
GL_CONST int DIHE_EACC_OFFSET = ANGL_EACC_OFFSET + BOND_WORK_UNIT_WARPS_PER_BLOCK;
GL_CONST int CMAP_EACC_OFFSET = DIHE_EACC_OFFSET + BOND_WORK_UNIT_WARPS_PER_BLOCK;
GL_CONST int QQXC_EACC_OFFSET = CMAP_EACC_OFFSET + BOND_WORK_UNIT_WARPS_PER_BLOCK;
GL_CONST int SCNB_EACC_OFFSET = QQXC_EACC_OFFSET + BOND_WORK_UNIT_WARPS_PER_BLOCK;
GL_CONST int SCEE_EACC_OFFSET = SCNB_EACC_OFFSET + BOND_WORK_UNIT_WARPS_PER_BLOCK;
GL_CONST int NMR2_EACC_OFFSET = SCEE_EACC_OFFSET + BOND_WORK_UNIT_WARPS_PER_BLOCK;
GL_CONST int NMR3_EACC_OFFSET = NMR2_EACC_OFFSET + BOND_WORK_UNIT_WARPS_PER_BLOCK;
GL_CONST int NMR4_EACC_OFFSET = NMR3_EACC_OFFSET + BOND_WORK_UNIT_WARPS_PER_BLOCK;
GL_CONST int UREY_EACC_OFFSET = NMR4_EACC_OFFSET + BOND_WORK_UNIT_WARPS_PER_BLOCK;
GL_CONST int CIMP_EACC_OFFSET = UREY_EACC_OFFSET + BOND_WORK_UNIT_WARPS_PER_BLOCK;
GL_CONST int CNST_EACC_OFFSET = CIMP_EACC_OFFSET + BOND_WORK_UNIT_WARPS_PER_BLOCK;
GL_CONST int EACC_TOTAL_SIZE  = CNST_EACC_OFFSET + BOND_WORK_UNIT_WARPS_PER_BLOCK;

// Scaling factors needed by bond work units to cast results to
// 32/64-bit ints in SP*P or DPFP precision modes, respectively.
#define BWU_SCALING_FACTOR    (1ll << 18)
#define BWU_PROMOTION_FACTOR  (1ll << 22)
GL_CONST double BSCALE                     = (double)BWU_SCALING_FACTOR;
GL_CONST double BPROMOTE                   = (double)BWU_PROMOTION_FACTOR;
GL_CONST float BSCALEF                     = (float)BWU_SCALING_FACTOR;
GL_CONST float BPROMOTEF                   = (float)BWU_PROMOTION_FACTOR;

// Indexing for GPU kernels - TODO appears not to be used, check where it might need to be applied
GL_CONST int BWU_YCRD_OFFSET                  =    BOND_WORK_UNIT_THREADS_PER_BLOCK + 8;
GL_CONST int BWU_ZCRD_OFFSET                  =  2*BOND_WORK_UNIT_THREADS_PER_BLOCK + 16;

// Codes for atom import sources
GL_CONST int ATOM_DP_SOURCE_CODE              =    0;
GL_CONST int ATOM_FP_SOURCE_CODE              =    1;

struct BondedWorkUnit {
  int natom;           // The number of unique atoms that must be imported for this unit
  int nbond;           // The number of bonds in this work unit (this, and any of the other
                       //   counts, could be zero)
  int nqqxc;           // The number of charge-charge interaction exclusions to be processed
                       //   by this work unit
  int nnb14;           // The number of 1-4 exclusions to be processed
  int nangl;           // The number of bond angle terms in this work unit
  int ndihe;           // The number of dihedral terms in this work unit
  int nnmr2;           // The number of two-body (bond) NMR restraints in this work unit
  int nnmr3;           // The number of three-body (angle) NMR restraints in this work unit
  int nnmr4;           // The number of four-body (dihedral) NMR restraints in this work unit
  int ncmap;           // The number of CMAP terms in this work unit
  int nurey;           // The number of (CHARMM) Urey-Bradley angle terms in this work unit
  int ncimp;           // The number of CHARMM improper dihedrals in this work unit
  int ncnst;           //
  int frcAcc;          // The force accumulator (0 = bonded, 1 = non-bonded) to dump into
  int bondDbl2Idx;     //
  int anglDbl2Idx;     //
  int nmr2Dbl2Idx;     //
  int nmr3Dbl2Idx;     //
  int nmr4Dbl2Idx;     //
  int ureyDbl2Idx;     // Indices at which this work unit starts reading various parameters
  int cimpDbl2Idx;     //   from the common arrays kept on the GPU to serve all work units.
  int cnstDbl2Idx;     //
  int dihePFloat2Idx;  //
  int dihePFloatIdx;   //
  int nb14PFloat2Idx;  //
  int nb14PFloatIdx;   //
  int qPFloatIdx;      //
  int *atomList;       // List of atoms to be imported (pointer to a region of data)
  int *bondList;       // List of bonds that this work unit will cover
  int *qqxcList;       // List of electrostatic exclusions that this work unit will cover
  int *nb14List;       // List of 1:4 exclusions that this work unit will cover
  int *bondAngleList;  // List of bond angles that this work unit will cover
  int *dihedralList;   // List of dihedrals (and impropers) that this work unit will cover
  int *cmapList;       // List of CMAP terms that this work unit will cover
  int *nmr2List;       // List of two-body NMR terms
  int *nmr3List;       // List of three-body NMR terms
  int *nmr4List;       // List of four-body NMR terms
  int *ureyList;       // List of Urey-Bradley terms
  int *cimpList;       // List of CHARMM improper dihedral terms
  int *cnstList;       // List of atom positional constraints
  int *bondMap;        //
  int *qqxcMap;        //
  int *nb14Map;        //
  int *bondAngleMap;   //
  int *dihedralMap;    // Maps of atoms in the imported array detailing which atom positions
  int *cmapMap;        //   the GPU shall reference in order to fulfill each interaction, and
  int *nmr2Map;        //   where to assign forces once the computations are complete.
  int *nmr3Map;        //
  int *nmr4Map;        //
  int *ureyMap;        //
  int *cimpMap;        //
  int *cnstMap;        //
  int *bondStatus;     //
  int *anglStatus;     //
  int *diheStatus;     //
  int *cmapStatus;     //
  int *qqxcStatus;     // Thermodynamic integration status of each object.  These are packed
  int *nb14Status;     //   bit strings, the high eight bits being blank, next highest eight
  int *nmr2Status;     //   bits being the TI region, after that soft-core status and then
  int *nmr3Status;     //   CV status.
  int *nmr4Status;     //
  int *ureyStatus;     //
  int *cimpStatus;     //
  int *cnstStatus;     //
  int* data;           // Data held by this work unit.  This is the actual array.  All the
                       //   pointers above just reference it.
};
typedef struct BondedWorkUnit bondwork;

struct BondedWorkAllocation {
  int nbondwarps;      // Number of warps devoted to bond stretching computations
  int nanglwarps;      // Number of warps devoted to bond angle computations
  int ndihewarps;      // Number of warps devoted to dihedral term computations
  int ncmapwarps;      // Number of warps devoted to CMAP computations
  int nqqxcwarps;      // Number of warps devoted to electrostatic exclusions
  int nnb14warps;      // Number of warps devoted to 1:4 non-bonded terms
  int nnmr2warps;      // Number of warps devoted to NMR distance restraints
  int nnmr3warps;      // Number of warps devoted to NMR angle restraints
  int nnmr4warps;      // Number of warps devoted to NMR torsion restraints
  int nureywarps;      // Number of warps devoted to Urey-Bradley angle terms
  int ncimpwarps;      // Number of warps devoted to Charmm improper dihedrals
  int ncnstwarps;      // Number of warps devoted to atom positional constraints
  int nUint;           // Aggregate unsigned int memory to allocate in a GpuBuffer
  int nDbl;            // Aggregate PMEDouble memory to allocate in a GpuBuffer
  int nDbl2;           // Aggregate PMEDouble2 memory to allocate in a GpuBuffer
  int nPFloat;         // Aggregate PMEFloat memory to allocate in a GpuBuffer
  int nPFloat2;        // Aggregate PMEFloat2 memory to allocate in a GpuBuffer
};
typedef BondedWorkAllocation bwalloc;

struct ParameterBuffer {
  int    i0;           // Integer value, the only use currently being to record the energy map
                       //   needed by a particular CMAP term
  double p1;           //
  double p2;           //
  double p3;           //
  double p4;           //
  double p5;           // Real-valued parameters, taken from PMEDouble, PMEDouble2, PMEFloat,
  double p6;           //   or PMEFloat2 types for bonds, angles, dihedrals, etc.  All values
  double p7;           //   get promoted to double precision--comparisons must be made within
  double p8;           //   tolerances.
  double p9;           //
  double p10;          //
  double p11;          //
  double p12;          //
};
typedef ParameterBuffer prmbuff;

struct AtomIDBuffer {
  int refI;            //
  int refJ;            // Atoms I, J, K, L, and M (as per Amber name conventions) referencing
  int refK;            //   positions in the master topology
  int refL;            //
  int refM;            //
};
typedef struct AtomIDBuffer aidbuff;

#endif
