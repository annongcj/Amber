#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// kPMECondenseChargeInfo: condense the information about each particle's position and charge
//                         into a single PMEfloat4 vector type.  This will convert each
//                         position into an unsigned integer, the first 10 bits of which
//                         indicate the grid bin it goes in and the last 22 of which indicate
//                         the position in the bin, using fractional coordinates.  That
//                         representation is four to eight times as precise as the floating
//                         point representation of the same coordinates in the neighbor list
//                         from which this routine draws.  The charge is represented as an
//                         unsigned integer to fit into the save uint4 vector type, but can be
//                         simply read as a float.  With the charge and position information
//                         condensed in this manner, the bulk of the work of converting each
//                         atom to a representation suitable for the interpolation is only done
//                         once per atom.  It is then simple to read the information from a
//                         texture to form lists of atoms specific to each subsection of the
//                         grid.
//---------------------------------------------------------------------------------------------
{
#ifdef PME_VIRIAL
  __shared__ PMEFloat sRecip[9];
  if (threadIdx.x < 9) {
    sRecip[threadIdx.x] = cSim.pNTPData->recip[threadIdx.x];
  }
  __syncthreads();
#endif

  // Standard blockwise striding through the list of atoms in the neighbor list hash cell
  uint2 clim = tex1Dfetch(nlboundstexref, blockIdx.x);
  int pos  = clim.x + threadIdx.x;

  // Compute the origin of the cell, and also the grid spacings,
  // in terms of fractional coordinates.  Every thread computes
  // these to have them stored in registers.
  int xyslab = cSim.xcells * cSim.ycells;
#if defined PME_IS_ORTHOGONAL && !defined(PME_VIRIAL)
  int zcellid = blockIdx.x / xyslab;
  int ycellid = (blockIdx.x - zcellid*xyslab) / cSim.xcells;
  int xcellid = blockIdx.x - zcellid*xyslab - ycellid*cSim.xcells;
  double xcfract = (double)xcellid / (double)cSim.xcells;
  double ycfract = (double)ycellid / (double)cSim.ycells;
  double zcfract = (double)zcellid / (double)cSim.zcells;
#endif
  while (pos < clim.y) {

    // Load x, y, and z coordinates from global memory
    PMEFloat2 xy = cSim.pAtomXYSP[pos];
    double xcrd = xy.x;
    double ycrd = xy.y;
    double zcrd = cSim.pAtomZSP[pos];

#ifdef PME_IS_ORTHOGONAL
    // Convert to fractional coordinates
#  ifdef PME_VIRIAL
    xcrd *= sRecip[0];
    ycrd *= sRecip[3];
    zcrd *= sRecip[6];
#  else
    xcrd = xcrd*cSim.recip[0][0] + xcfract;
    ycrd = ycrd*cSim.recip[1][1] + ycfract;
    zcrd = zcrd*cSim.recip[2][2] + zcfract;
#  endif
#endif
    // Compute the first grid cell to which this atom shall contribute.
    xcrd += (PMEDouble)((xcrd < 0.0) - (xcrd >= 1.0));
    ycrd += (PMEDouble)((ycrd < 0.0) - (ycrd >= 1.0));
    zcrd += (PMEDouble)((zcrd < 0.0) - (zcrd >= 1.0));
    xcrd *= cSim.nfft1d;
    ycrd *= cSim.nfft2d;
    zcrd *= cSim.nfft3d;
    PMEUint ixcrd = xcrd;
    PMEUint iycrd = ycrd;
    PMEUint izcrd = zcrd;

    // Finally, load charge information from global memory (hopefully,
    // the latency can be masked by the arithmetic work above)
    PMEFloat2 qljid = cSim.pAtomChargeSPLJID[pos];
#if defined(use_SPFP)
    unsigned int xrem = (xcrd - (PMEDouble)ixcrd)*PMEGRID_SCALEF;
    unsigned int yrem = (ycrd - (PMEDouble)iycrd)*PMEGRID_SCALEF;
    unsigned int zrem = (zcrd - (PMEDouble)izcrd)*PMEGRID_SCALEF;
    uint4 crdq;
    crdq.x = (ixcrd << 22) | xrem;
    crdq.y = (iycrd << 22) | yrem;
    crdq.z = (izcrd << 22) | zrem;
    crdq.w = __float_as_uint(qljid.x);
    cSim.pAtomCRDQ[pos] = crdq;
#else  // use_DPFP
    unsigned long long int xrem = (xcrd - (PMEDouble)ixcrd)*PMEGRID_SCALE;
    unsigned long long int yrem = (ycrd - (PMEDouble)iycrd)*PMEGRID_SCALE;
    unsigned long long int zrem = (zcrd - (PMEDouble)izcrd)*PMEGRID_SCALE;
    llconstruct Lx, Ly, Lz, Lq;
    Lx.ulli = (ixcrd << 54) | xrem;
    Ly.ulli = (iycrd << 54) | yrem;
    Lz.ulli = (izcrd << 54) | zrem;
    Lq.d = qljid.x;
    uint4 crdq0, crdq1;
    crdq0.x = Lx.ui[0];
    crdq0.y = Ly.ui[0];
    crdq0.z = Lz.ui[0];
    crdq0.w = Lq.ui[0];
    crdq1.x = Lx.ui[1];
    crdq1.y = Ly.ui[1];
    crdq1.z = Lz.ui[1];
    crdq1.w = Lq.ui[1];
    cSim.pAtomCRDQ[2*pos    ] = crdq0;
    cSim.pAtomCRDQ[2*pos + 1] = crdq1;
#endif
    pos += blockDim.x;
  }
}
