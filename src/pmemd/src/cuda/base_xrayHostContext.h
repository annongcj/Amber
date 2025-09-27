#ifndef _BASE_XRAY_HOST_CONTEXT
#define _BASE_XRAY_HOST_CONTEXT

#include "gputypes.h"
#include "gpuBuffer.h"
#include "xrayDevConstants.h"

class base_xrayHostContext {

public:
  virtual void Init();
  
  base_xrayHostContext();
  virtual ~base_xrayHostContext();
  
public:
  xrayDevConstants devInfo;                 // Bundle of parameters and pointers to ship to the
                                            //   GPU, to be referenced from constant memory
  GpuBuffer<int>* pbSelectionMask;          // Atom selection mask (array of zeros and ones,
                                            //   one for atoms that have been selected)
  GpuBuffer<int>* pbSelectedAtoms;          // Compact list of atoms selected
  GpuBuffer<int>* pbHKLaidx;                //
  GpuBuffer<int>* pbHKLbidx;                // H-K-L indices for reflections
  GpuBuffer<int>* pbHKLcidx;                //
  GpuBuffer<int>* pbHKLMask;                // Mask for H-K-L indices to apply
  GpuBuffer<int>* pbTestSelection;          //
  GpuBuffer<int>* pbWorkSelection;          //
  GpuBuffer<int>* pbScatterType;            // Scattering types for each atom (there are not
                                            //   very many different values--one for each
                                            //   element type at most)
  GpuBuffer<int>* pbResNum;                 // Number of the residue to which each atom belongs
  GpuBuffer<unsigned int>* pbResChain;      // Integer-encoded chain IDs (four characters) for
                                            //   all atoms in the simulation
  GpuBuffer<unsigned int>* pbResICode;      // Residue integer codes, similar to the
                                            //   4-character encoded format above
  GpuBuffer<unsigned int>* pbIElement;      // Element integer codes for all atoms
  GpuBuffer<unsigned int>* pbAltLoc;        // 
  GpuBuffer<PMEFloat>* pbAtomBFactor;       // B-factors for all atoms
  GpuBuffer<PMEFloat>* pbOccupancy;         // Occupancies for all atoms
  GpuBuffer<PMEFloat>* pbAbsFObs;           // Absolute values of observed structure factors
                                            //   for each reflection
  GpuBuffer<PMEFloat2>* pbFObs;             // Observed structure factors for each reflection
  GpuBuffer<PMEFloat>* pbSigFObs;           // Sigmas for FObs
  GpuBuffer<PMEFloat>* pbMSS4;              // 
  GpuBuffer<PMEFloat>* pbSFScaleFactors;    // Storage for structure factor normalizations,
                                            //   needed to match CPU code behavior
  GpuBuffer<PMEFloat>* pbSolventKMask;      // Bulk solvent mask contribution (calculated once
                                            //   at the outset of the run)
  GpuBuffer<PMEFloat2>* pbSolventFMask;     // Bulk solvent structure factor contributions
                                            //   (updated periodically during the simulation)
  GpuBuffer<PMEFloat>* pbScatterCoeffs;     // Scattering coefficients (array is linearized,
                                            //   rather than the three-dimensional array used
                                            //   in the Fortran code)
  GpuBuffer<PMEFloat>* pbScatterFactorLkp;  // Lookup table for pre-computed scattering factors
                                            //   of each reflection, for a given set of
                                            //   box dimensions
  GpuBuffer<PMEFloat2>* pbFcalc;            // Complex representations of the computed
                                            //   structure factors for each reflection
  GpuBuffer<PMEFloat>* pbAbsFcalc;          // Magnitudes of the computed structure factors
  GpuBuffer<PMEFloat2>* pbReflDeriv;        // Complex derivatives of the structure factors
  GpuBuffer<PMEAccumulator>* pbXrayResults; // Array to hold accumulated, system-wide results
  GpuBuffer<PMEAccumulator>* pbXrayForces;  // Array holding X, Y, and Z forces in contiguous
                                            //   blocks
};

#endif
