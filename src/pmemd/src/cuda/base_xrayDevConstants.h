#ifndef _BASE_XRAY_DEV_CONSTANTS
#define _BASE_XRAY_DEV_CONSTANTS

#include "gputypes.h"

class base_xrayDevConstants {

public:
  void InitXrayDevConstants();
  
public:
  int nHKL;                      // Number of H-K-L reflections
  int natom;                     // The number of atoms in the simulation (same as in the
                                 //   simulationConst struct)
  int nres;                      // The number of residues in the simulation (same as in the
                                 //   simulationConst struct)
  int nSelectedAtom;             // The number of atoms selected for X-ray structure factor
                                 //   calculations
  int tarfunc;                   // Target function for linking deviations in calculated
                                 //   X-ray structure factors to forces on atoms.  0 = overall
                                 //   harmonic restraint, 1 = vector-based RMSD restraint
  
  // General parameters for structure factor calculations
  int spaceGroupID;              // Identification number for the space group
  int nScatterCoeff;             // The number of scattering coefficients for each atom type
  int nScatterTypes;             // Number of atom types with different scattering properties
  int SolventModel;              // Identifier code for the solvent model in use (default 0,
                                 //   no solvent, 1 = 'simple' bulk solvent mask)
  PMEFloat lowResolution;        // Lower resolution limit
  PMEFloat highResolution;       // Upper resolution limit
  PMEFloat xrayGeneralWt;        // Weight term for forces derived from X-ray restraints
  PMEFloat BFacRefineInterval;   // B-factor refinement interval
  PMEFloat minBFactor;           // Minimum value for B-factors (default 1.0)
  PMEFloat maxBFactor;           // Maximum value for B-factors (default 999.0)
  PMEFloat NormScale;            // Sum of squares of observed structure factors
  PMEFloat SumFObsWork;          // Sum of absolute values of observed structure factors for
                                 //   unmasked reflections
  PMEFloat SumFObsFree;          // Sum of absolute values of observed structure factors for
                                 //   masked reflections
  
  // Solvent mask parameters
  int bsMaskModel;               // Index code for the type of bulk solvent employed.
                                 //   0 : no solvent
                                 //   1 : "simple" solvent with exlcusion-based grid coloring
  int bsMaskUpdateInterval;      // Update interval for computing the bulk solvent mask
  int bsFFTMethod;               // Method for computing the mask
  int bsMaskXdim;                //
  int bsMaskYdim;                // Bulk solvent mask X, Y, and Z grid dimensions
  int bsMaskZdim;                //
  PMEFloat solventProbeRadius;   // Probe radius
  PMEFloat solventMaskExpand;    // Expansion thickness for the second phase of grid coloring
  PMEFloat solventScale;         // The scaling coefficient for solvent contributions to the
                                 //   structure factors
  PMEFloat solventBfactor;       // Assumed B-factor for solvent atoms
  
  // Pointers to data in the overarching xrayHostContext class instance
  int* pSelectionMask;           // Mask of selected atoms out of all particles in the system
  int* pSelectedAtoms;           // List of selected atoms
  int* pHKLaidx;                 //
  int* pHKLbidx;                 // First, second, and third H-K-L reflection indices
  int* pHKLcidx;                 //
  int* pHKLMask;                 // Mask for H-K-L reflections to apply
  int* pScatterType;             // Atom scattering types
  PMEFloat* pAtomBFactor;        // Isotropic atomic thermal B-factors
  PMEFloat* pOccupancy;          // Atom occupancies
  PMEFloat* pAbsFObs;            // Absolute values of the observed structure factors
  PMEFloat2* pFObs;              // Observed structure factors
  PMEFloat* pSigFObs;            // Sigmas for observed structure factor calculation
  PMEFloat* pSFScaleFactors;     // Scaling factors derived from structure factor sums across
                                 //   all reflections
  PMEFloat* pMSS4;               // Pointer to reflection data
  PMEFloat*  pSolventKMask;      // Pointer to solvent contributions, static reflection mask
  PMEFloat2* pSolventFMask;      // Pointer to solvent contributions, periodically updated for
                                 //   each reflection
  PMEFloat* pScatterCoeffs;      // Scattering coefficients (this array is linearized, in
                                 //   contrast to the Fortran code's implementation)
  PMEFloat2* pFcalc;             // Complex representations of the computed structure factors
  PMEFloat* pAbsFcalc;           // Magnitudes of the computed structure factors
  PMEFloat2* pReflDeriv;         // Derivative of the energy function with respect to each
                                 //   structure factor
  PMEAccumulator* pForceX;       // Forces acting on atoms in the X direction.  This points to
                                 //   the beginning of a larger array holding all X, Y, and Z
                                 //   forces in contiguous memory.
  PMEAccumulator* pForceY;       // Forces acting on atoms in the Y direction
  PMEAccumulator* pForceZ;       // Forces acting on atoms in the Z direction

  // Results (individual floating point numbers) derived from structure factor calculations
  PMEAccumulator* pSumFoFc;      //
  PMEAccumulator* pSumFoFo;      // Sums of the Fobs and Fcalc products
  PMEAccumulator* pSumFcFc;      //
  PMEAccumulator* pSumFobs;      // Sum of Fobs
  PMEAccumulator* pResidual;     // Residual from R-factor calculations (this is r_work)
  PMEAccumulator* pFreeResidual; // Residual from R-factor calculations (this is r_free)
  PMEAccumulator* pXrayEnergy;   // Energy due to X-ray restraints
};

#endif
