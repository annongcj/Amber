#ifndef BOND_REMAP_FUNCS
#define BOND_REMAP_FUNCS

#include "bondRemapDS.h"
#include "gputypes.h"
#include "gpuContext.h"

void DestroyBondWorkUnit(bondwork *bw);

bondwork* AssembleBondWorkUnits(gpuContext gpu, int *nbwunits);

bwalloc CalcBondBlockSize(gpuContext gpu, bondwork* bwunits, int nbwunits);

void MakeBondedUnitDirections(gpuContext gpu, bondwork *bwunits, int nbwunits,
                              bwalloc *bwdims);

void CheckWorkUnitTranslations(gpuContext gpu, bondwork* bwunits, const char* objtype);

#endif
