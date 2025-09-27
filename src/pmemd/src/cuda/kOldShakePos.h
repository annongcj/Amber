#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// Dec 2021, by zhf
//---------------------------------------------------------------------------------------------
//
{
#ifdef SHAKE_NEIGHBORLIST
#define ATOMX(i) cSim.pImageX[i]
#define ATOMY(i) cSim.pImageY[i]
#define ATOMZ(i) cSim.pImageZ[i]
#else
#define ATOMX(i) cSim.pAtomX[i]
#define ATOMY(i) cSim.pAtomY[i]
#define ATOMZ(i) cSim.pAtomZ[i]
#endif

    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;

    if (pos < cSim.atoms) {
        cSim.pShakeOldAtomX[pos] = ATOMX(pos);
        cSim.pShakeOldAtomY[pos] = ATOMY(pos);
        cSim.pShakeOldAtomZ[pos] = ATOMZ(pos);
    }

#undef ATOMX
#undef ATOMY
#undef ATOMZ
}

