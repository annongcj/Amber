#ifndef HONEYCOMB_FUNCTIONS
#define HONEYCOMB_FUNCTIONS

#include "matrixDS.h"
#include "gpuContext.h"
#include "gpuHoneycombDS.h"

int CountQQAtoms(gpuContext gpu);

int CountLJAtoms(gpuContext gpu);

int GetNTInstructionEst(gpuContext gpu);

nlkit GetHoneycombDesign(gpuContext gpu, double atmcrd[][3], double es_cutoff,
                         double vdw_cutoff, double nbskin);

void DestroyNLKit(nlkit *tNL);

void PackPencilRelationship(gpuContext gpu, imat *P, dmat *R, int span, imat *ntPencils,
                            imat *ntInteractions);

void MakeSegmentRelays(gpuContext gpu, imat *ntInteractions);

int SortInt2ByX(const void *dualA, const void *dualB);

int FindSubImgOrderingPartition(int natom, int nsmp, int maxcell);

int FindExpansionBatchCount(gpuContext gpu);

void CheckHoneycombExpansion(gpuContext gpu, bool printAll=false, bool checkPairs=true);

#endif
