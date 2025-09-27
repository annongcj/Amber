#ifndef GPU_HONEYCOMB_STRUCTS
#define GPU_HONEYCOMB_STRUCTS

// Enforce correct usage of const in different env.(CPU or GPU)
#ifdef  __CUDA_ARCH__
#define GL_CONST  static const __constant__
#else
#define GL_CONST  static const
#endif

GL_CONST int NTMAP_SPACE          = 160;
GL_CONST int NTMAP_TYPE_IDX       = 132;
GL_CONST int NTMAP_ATOM_COUNT_IDX = 133;
GL_CONST int NTMAP_PAIR_START_IDX = 134;
GL_CONST int NTMAP_PAIR_COUNT_IDX = 135;
GL_CONST int NTMAP_TILE_START_IDX = 136;
GL_CONST int NTMAP_TILE_COUNT_IDX = 137;
GL_CONST int NTMAP_SIZE           = 138;
GL_CONST int NT_MAX_TILES         = 768;
GL_CONST int HCMB_PL_STORAGE_SIZE = 262144;
GL_CONST int HCMB_PL_PADDING_SIZE = 40960;
#ifdef use_DPFP
GL_CONST int NT_MAX_ATOMS         = 784;
GL_CONST int NT_PAD_ATOMS         = 800;
GL_CONST int HCMB_PL_BLANK        = ((0x318 << 16) | 0x310);
GL_CONST int HCMB_TILE_BLANK      = ((0x318 << 16) | 0x310);
#else
GL_CONST int NT_MAX_ATOMS         = 1008;
GL_CONST int NT_PAD_ATOMS         = 1024;
GL_CONST int HCMB_PL_BLANK        = ((0x3f8 << 16) | 0x3f0);
GL_CONST int HCMB_TILE_BLANK      = ((0x3f8 << 16) | 0x3f0);
#endif

#include "matrixDS.h"

struct neighborListKit {
  bool yshift;            // Flag to indicate hash cells are staggered along the Y-axis
  int ny;                 // Number of pencil cells along the Y-axis
  int nz;                 // Number of pencil cells along the Z-axis
  int TotalQQSpace;       // Total allocation needed for the expanded electrostatic or
  int TotalLJSpace;       //   Lennard-Jones images (stored here as they are not needed on
                          //   the device)
  int TotalPairSpace;     // Total allocation needed for the combined pair list storage.  The
                          //   pair list is stored as bit packed pairs of integers which refer
                          //   to the imported atom numbers of a specific Neutral Territory
                          //   import region.  Because the NT region could contain
                          //   electrostatic or Lennard-Jones particles and will point to a
                          //   specific index of the pair list array in global memory, one
                          //   array will be able to store all of the relevant pairs data.
  imat qqPencils;         // Pencil cells that interact with the cell at the origin to cover
  imat ljPencils;         //   all electrostatic and Lennard-Jones interactions, respectively
  dmat qqRanges;          // Ranges between pencil cells in the electrostatic and
  dmat ljRanges;          //   Lennard-Jones lists above
};
typedef struct neighborListKit nlkit;

struct CellCoordinates {
  int ycell;              // Y cell index to which the atom will be assigned
  int zcell;              // Z cell index to which the atom will be assigned
  double newx;            // New coordinates for the atom placed as close as possible to the
  double newy;            //   central axis of its home cell such that its X coordinate lies
  double newz;            //   in the primary unit cell
};
typedef struct CellCoordinates cellcrd;

struct ImageCellContents {
  int ypos;         // Y position of the hash cell in the expanded grid
  int zpos;         // Z position of the hash cell in the expanded grid
  int ysrc;         // Y and Z positions of the primary image hash cell that will fill this 
  int zsrc;         //   hash cell (ysrc and zsrc are indexed in the expanded image but must
                    //   point to a cell not in the extended, padding regions of the grid)
  int coreStart;    // Start and end of the core contents--all atoms in the cell excluding
  int coreEnd;      //   those in capping groups
  int lcapStart;    // Start of the lower capping group and end of the upper capping group
  int hcapEnd;      //   for ghost atoms in the cell replicated along the X axis
  int nntmap;       // The number of Nuetral Territory maps owned by this hash cell
  int *atoms;       // List of atoms in the hash cell, indexed according to the resorted list
                    //   in HcmbIndex
  int *absid;       // List of atom ID numbers according to the master topology but tracking
                    //   the order in the 'atoms' array
  int* ntmaps;      // The Neutral Territory maps owned by this hash cell
  int* pairlimits;  // The limits of each Neutral Territory region's pair list within the
                    //   pairs attribute (there are nntmap + 1 elements in this array, the kth
                    //   Neutral Territory region extends from limits defined in indices k to
                    //   k + 1)
  int* pairs;       // The pair list handled by this Neutral Territory region.  As on the GPU,
                    //   the entries of this list are stored as packed bit strings.  There is
                    //   only the array of pairs, no tiles, on the CPU.
  PMEFloat4 *crdq;  // Coordinates and properties of atoms in the hash cell, listed according
                    //   to the resorted list in HcmbIndex and tracking the 'atoms' array
};
typedef struct ImageCellContents icontent;

struct ExpandedImage {
  int xpadding;        // Padding in the X direction--this is a number of ATOMS, not a number 
                       //   of hash cells
  int ypadding;        // Hash cell padding in the Y direction
  int zpadding;        // Hash cell padding in the Z direction
  int ydimFull;        // The complete width of the expanded hash cell grid in the Y direction
  int zdimFull;        // The complete width of the expanded hash cell grid in the Z direction
  int ydimPrimary;     // The width of the primary image's hash cell grid in the Y direction
  int zdimPrimary;     // The width of the primary image's hash cell grid in the Z direction
  int cellspace;       // The maximum number of atoms that can be held in each cell
  int* atoms;          // Aggregated list of atoms for every hash cell
  int* absid;          // List of atom ID numbers according to the master topology but tracking
                       //   the order in the 'atoms' array
  PMEFloat4* crdq;     // Coordinates and properties of atoms in the expanded representation
  icontent* cells;     // The list of pencil hash cells into which atoms are packed
};
typedef struct ExpandedImage eimage;

struct PairlistChecker {
  int natom;             // The number of atoms in the system, copied over for convenience 
  int* ngbrCounts;       // Neighbor counts for each atom, a prefix sum denoting starting
                         //   points to read in ngbrList and ngbrRanges
  int* ngbrList;         // The list of neighbors for each atom
  PMEFloat* ngbrRanges;  // The ranges of all neighbors to each atom
};
typedef struct PairlistChecker pairlist;

struct PairTileScoreCard {
  int n1x1;
  int n4x4;
  int n8x8;
  int nTrace;
};
typedef struct PairTileScoreCard tilescore;

#endif
