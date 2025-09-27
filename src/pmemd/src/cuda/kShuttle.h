#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// This file is included by kForceUpdate.cu to construct a number of kernels for downloading
// and uploading data on a subset of the atoms.  Each thread block is to contain 96 threads to
// manage transfers of 32 atoms per block.
//
// Defines: TAKE_IMAGE (used when namelists are in effect and atoms are re-arranged into an
//          image of their initial order), PULL and POST (direction of the transfer), COORD
//          and FORCE (type of transfer)
//---------------------------------------------------------------------------------------------
{
  __shared__ int atidx[GRID];
  int warpIdx = threadIdx.x >> GRID_BITS;
 unsigned int tgx = threadIdx.x & GRID_BITS_MASK;
  int pos = tgx + (GRID * blockIdx.x);
#ifdef PULL_FORCE
  bool frcsplit = (cSim.ntp > 0) && (cSim.barostat == 1);
#endif
  while (pos < cSim.nShuttle) {

    // Read the ID numbers in the shuttle tickets for a warp-sized batch of atoms
    if (warpIdx == 0) {
#ifdef TAKE_IMAGE
      int shidx = cSim.pShuttleTickets[pos];
      atidx[tgx] = cSim.pImageAtomLookup[shidx];
#else
      atidx[tgx] = cSim.pShuttleTickets[pos];
#endif
    }
    __syncthreads();
#ifdef PULL_COORD

    // Pull coordinates
#  ifdef TAKE_IMAGE
    double atmcrd = cSim.pImageX[atidx[tgx] + warpIdx*cSim.stride];
#  else
    double atmcrd = cSim.pAtomX[atidx[tgx] + warpIdx*cSim.stride];
#  endif
    cSim.pDataShuttle[warpIdx*cSim.nShuttle + pos] = atmcrd;
#elif defined (PULL_FORCE)

    // Pull forces.  If the forces are split (bonded and non-bonded values),
    // fold them together on the GPU to reduce download volume and CPU work.
    PMEAccumulator fr = cSim.pForceAccumulator[atidx[tgx] + warpIdx*cSim.stride];
    if (frcsplit) {
      fr += cSim.pForceAccumulator[atidx[tgx] + warpIdx*cSim.stride + cSim.stride3];
    }
    cSim.pDataShuttle[warpIdx*cSim.nShuttle + pos] = (double)fr * ONEOVERFORCESCALE;
#elif defined (PULL_BOND_FORCE)

    // Pull forces from the first of what could be two accumulators.
    // If the forces are NOT split between two accumulators this will pull
    // down the entirety of the force (bonded and non-bonded values).
    // This is really only intended to be used in cases where the forces
    // are split between two accumulators, to pull down the bonded term
    // contributions to the forces.
    PMEAccumulator fr = cSim.pForceAccumulator[atidx[tgx] + warpIdx*cSim.stride];
    cSim.pDataShuttle[warpIdx*cSim.nShuttle + pos] = (double)fr * ONEOVERFORCESCALE;
#elif defined (PULL_NONBOND_FORCE)

    // Pull forces from the second of what must be two accumulators.
    // If the forces are NOT split between two accumulators this will fail.
    // This is only intended to be used in cases where the forces are split
    // between two accumulators, to pull down the non-bonded term
    // contributions to the forces.
    PMEAccumulator fr = cSim.pForceAccumulator[atidx[tgx] + warpIdx*cSim.stride +
                                               cSim.stride3];
    cSim.pDataShuttle[warpIdx*cSim.nShuttle + pos] = (double)fr * ONEOVERFORCESCALE;
#elif defined(POST_FORCE)

    // Post forces directly into GPU arrays
    PMEAccumulator ifr = cSim.pDataShuttle[warpIdx*cSim.nShuttle + pos] * FORCESCALE;
    cSim.pForceAccumulator[atidx[tgx] + warpIdx*cSim.stride] = ifr;
#elif defined(POST_NONBOND_FORCE)

    // Post forces directly into GPU arrays
    PMEAccumulator ifr = cSim.pDataShuttle[warpIdx*cSim.nShuttle + pos] * FORCESCALE;
    cSim.pForceAccumulator[atidx[tgx] + warpIdx*cSim.stride + cSim.stride3] = ifr;
#elif defined(POST_FORCE_ADD)

    // Post forces and add them to the existing values
    PMEAccumulator ifr = cSim.pDataShuttle[warpIdx*cSim.nShuttle + pos] * FORCESCALE;
    cSim.pForceAccumulator[atidx[tgx] + warpIdx*cSim.stride] += ifr;
#elif defined(POST_NONBOND_FORCE_ADD)

    // Post forces and add them to the existing values
    PMEAccumulator ifr = cSim.pDataShuttle[warpIdx*cSim.nShuttle + pos] * FORCESCALE;
    cSim.pForceAccumulator[atidx[tgx] + warpIdx*cSim.stride + cSim.stride3] += ifr;
#elif defined(POST_COORD)

    // Post coordinates directly into the GPU arrays
    double atmcrd = cSim.pDataShuttle[warpIdx*cSim.nShuttle + pos];
#  ifdef TAKE_IMAGE
    cSim.pImageX[atidx[tgx] + warpIdx*cSim.stride] = atmcrd;
#  else
    cSim.pAtomX[atidx[tgx] + warpIdx*cSim.stride] = atmcrd;
#  endif
#endif

    // Increment the position--three warps in each block
    // collectively operate on GRID (GRID = 32) atoms
    pos += gridDim.x * GRID;
  }
}
