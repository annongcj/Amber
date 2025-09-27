#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------
#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include <string.h>
#include "matrix.h"
#include "bondRemap.h"
#include "gputypes.h"

using namespace std;

//---------------------------------------------------------------------------------------------
// codify_object_type: codify the type of object and set pointers for atom ID lookup.  Macro
//                     FTW because pointer arithmetic would be horrible otherwise.
//---------------------------------------------------------------------------------------------
#ifndef codify_object_type
#define codify_object_type(objtype, gpu, p1atm, p2atm, p4atm, objcode) \
{ \
  if (strcmp(objtype, "BOND") == 0) { \
    p2atm = gpu->pbBondID->_pSysData; \
    objcode = BOND_CODE; \
  } \
  else if (strcmp(objtype, "ANGLE") == 0) { \
    p2atm = gpu->pbBondAngleID1->_pSysData; \
    p1atm = gpu->pbBondAngleID2->_pSysData; \
    objcode = ANGL_CODE; \
  } \
  else if (strcmp(objtype, "DIHEDRAL") == 0) { \
    p4atm = gpu->pbDihedralID1->_pSysData; \
    objcode = DIHE_CODE; \
  } \
  else if (strcmp(objtype, "CMAP") == 0) { \
    p4atm = gpu->pbCmapID1->_pSysData; \
    p1atm = gpu->pbCmapID2->_pSysData; \
    objcode = CMAP_CODE; \
  } \
  else if (strcmp(objtype, "QQXC") == 0) { \
    p2atm = gpu->pQQxcID; \
    objcode = QQXC_CODE; \
  } \
  else if (strcmp(objtype, "NB14") == 0) { \
    p2atm = gpu->pbNb14ID->_pSysData; \
    objcode = NB14_CODE; \
  } \
  else if (strcmp(objtype, "NMR2") == 0) { \
    p2atm = gpu->pbNMRDistanceID->_pSysData; \
    objcode = NMR2_CODE; \
  } \
  else if (strcmp(objtype, "NMR3") == 0) { \
    p2atm = gpu->pbNMRAngleID1->_pSysData; \
    p1atm = gpu->pbNMRAngleID2->_pSysData; \
    objcode = NMR3_CODE; \
  } \
  else if (strcmp(objtype, "NMR4") == 0) { \
    p4atm = gpu->pbNMRTorsionID1->_pSysData; \
    objcode = NMR4_CODE; \
  } \
  else if (strcmp(objtype, "UREY") == 0) { \
    p2atm = gpu->pbUBAngleID->_pSysData; \
    objcode = UREY_CODE; \
  } \
  else if (strcmp(objtype, "CIMP") == 0) { \
    p4atm = gpu->pbImpDihedralID1->_pSysData; \
    objcode = CIMP_CODE; \
  } \
  else if (strcmp(objtype, "CNST") == 0) { \
    p1atm = gpu->pbConstraintID->_pSysData; \
    objcode = CNST_CODE; \
  } \
}
#endif

//---------------------------------------------------------------------------------------------
// FindAtomInUnit: determine whether a particular atom is within a given work unit by
//                 referencing tables of atoms as they are included in all work units.
//---------------------------------------------------------------------------------------------
static int FindAtomInUnit(int unitID, imat *unitMap, int* unitMapCounts, int atomID)
{
  int i;

  int imax = unitMapCounts[atomID];
  int *itmp;
  itmp = unitMap->map[atomID];
  int found = 0;
  for (i = 0; i < imax; i++) {
    found += (itmp[i] == unitID);
  }

  return found;
}

//---------------------------------------------------------------------------------------------
// GetWorkUnitObjectSize: get the size of a work unit object based on its code.
//
// Arguments:
//   objcode:  identifies the type of object (2 = bonds, 3 = bond angles, etc.)
//---------------------------------------------------------------------------------------------
static int GetWorkUnitObjectSize(int objcode)
{
  int objsize = 0;

  if (objcode == CNST_CODE) {
    objsize = 1;
  }
  else if (objcode == BOND_CODE || objcode == QQXC_CODE || objcode == NB14_CODE ||
    objcode == NMR2_CODE || objcode == UREY_CODE) {
    objsize = 2;
  }
  else if (objcode == ANGL_CODE || objcode == NMR3_CODE) {
    objsize = 3;
  }
  else if (objcode == DIHE_CODE || objcode == NMR4_CODE || objcode == CIMP_CODE) {
    objsize = 4;
  }
  else if (objcode == CMAP_CODE) {
    objsize = 5;
  }

  return objsize;
}

//---------------------------------------------------------------------------------------------
// FindObjectAtomsInUnit: for any type of object, find whether the atoms that make up the
//                        object are already present in a bond work unit.  Return the number of
//                        atoms already in the work unit and the number of atoms that would
//                        have to be added in order to include this object in the work unit.
//
//---------------------------------------------------------------------------------------------
static int2 FindObjectAtomsInUnit(int unitID, imat *unitMap, imat *unitMapCounts, int cnobj,
                                  int objcode, int* p1atm, int2* p2atm, int4* p4atm,
                                  int* ispresent)
{
  int i, objsize;

  // Find the object size
  objsize = GetWorkUnitObjectSize(objcode);

  // Initalize the atom presence
  for (i = 0; i < 5; i++) {
    ispresent[i] = 0;
  }

  // Determine the footprint in the unit as it stands.  Special-case unitMapCounts == NULL,
  // meaning that unitMap is merely a 1 x (number of simulation atoms) vector that stores
  // the ID of the last unit that any atom was placed in.  The special case is instigated by
  // BaseNGroups >> TraceBranch, when new groups are being built one at a time.  The more
  // general case is instigated by ExtendWorkUnits, when many units may have each atom.
  int *itmp;
  itmp = unitMap->data;
  if (unitMapCounts == NULL) {
    if (objcode == CNST_CODE) {
      ispresent[0] += (itmp[p1atm[cnobj]] == unitID);
    }
    if (objcode == BOND_CODE || objcode == ANGL_CODE || objcode == QQXC_CODE ||
        objcode == NB14_CODE || objcode == NMR2_CODE || objcode == UREY_CODE) {
      ispresent[0] += (itmp[p2atm[cnobj].x] == unitID);
      ispresent[1] += (itmp[p2atm[cnobj].y] == unitID);
    }
    if (objcode == ANGL_CODE || objcode == NMR3_CODE) {
      ispresent[2] += (itmp[p1atm[cnobj]] == unitID);
    }
    if (objcode == DIHE_CODE || objcode == CMAP_CODE || objcode == NMR4_CODE ||
        objcode == CIMP_CODE) {
      ispresent[0] += (itmp[p4atm[cnobj].x] == unitID);
      ispresent[1] += (itmp[p4atm[cnobj].y] == unitID);
      ispresent[2] += (itmp[p4atm[cnobj].z] == unitID);
      ispresent[3] += (itmp[p4atm[cnobj].w] == unitID);
    }
    if (objcode == CMAP_CODE) {
      ispresent[4] += (itmp[p1atm[cnobj]] == unitID);
    }
  }
  else {
    if (objcode == CNST_CODE) {
      ispresent[0] = FindAtomInUnit(unitID, unitMap, unitMapCounts->data, p1atm[cnobj]);
    }
    if (objcode == BOND_CODE || objcode == ANGL_CODE || objcode == QQXC_CODE ||
        objcode == NB14_CODE || objcode == NMR2_CODE || objcode == UREY_CODE) {
      ispresent[0] = FindAtomInUnit(unitID, unitMap, unitMapCounts->data, p2atm[cnobj].x);
      ispresent[1] = FindAtomInUnit(unitID, unitMap, unitMapCounts->data, p2atm[cnobj].y);
    }
    if (objcode == ANGL_CODE || objcode == NMR3_CODE) {
      ispresent[2] = FindAtomInUnit(unitID, unitMap, unitMapCounts->data, p1atm[cnobj]);
    }
    if (objcode == DIHE_CODE || objcode == CMAP_CODE || objcode == NMR4_CODE ||
        objcode == CIMP_CODE) {
      ispresent[0] = FindAtomInUnit(unitID, unitMap, unitMapCounts->data, p4atm[cnobj].x);
      ispresent[1] = FindAtomInUnit(unitID, unitMap, unitMapCounts->data, p4atm[cnobj].y);
      ispresent[2] = FindAtomInUnit(unitID, unitMap, unitMapCounts->data, p4atm[cnobj].z);
      ispresent[3] = FindAtomInUnit(unitID, unitMap, unitMapCounts->data, p4atm[cnobj].w);
    }
    if (objcode == CMAP_CODE) {
      ispresent[4] = FindAtomInUnit(unitID, unitMap, unitMapCounts->data, p1atm[cnobj]);
    }
  }

  // Score the result
  int2 presence;
  presence.x = 0;
  for (i = 0; i < 5; i++) {
    ispresent[i] = (ispresent[i] > 0);
    presence.x += ispresent[i];
  }
  presence.y = objsize - presence.x;

  return presence;
}

//---------------------------------------------------------------------------------------------
// ScoreObjectPlacementInUnit: compute a score for placing an object in a unit.  The score is
//                             increased for every atom that would need to be added to the unit
//                             in order for it to cover the object, but decreased for every
//                             terminal atom that the object brings with it.
//
// Arguments:
//   bw:             the bond work unit to which the object will be added
//   unitID:         the ID number of the bond work unit in the alrger bwunits array
//   unitMap:        map of all the units that each atom is a part of
//   unitMapCounts:  the total number of work units in which each atom is included
//   cnobj:          ID number of the object in question
//   objcode:        code number for the object type
//   p1atm:          array of single atoms included in certain object types (i.e. all atom Ks
//                   in bond angles)
//   p2atm:          array of two-atom tuples included in certain object types (i.e. both atoms
//                   in any given bond)
//   p4atm:          array of four-atom tuples included in certain object types (i.e. all four
//                   atoms in a dihedral, or atoms IJKL of a CMAP term)
//   atomEndPts:     table indicating which atoms connect to their molecules by only one bond
//   objTerminates:  table indicating which objects contain terminal atoms
//   creditTermini:  score to credit for objects containing terminal atoms
//---------------------------------------------------------------------------------------------
static int ScoreObjectPlacementInUnit(bondwork *bw, int unitID, imat *unitMap,
                                      imat *unitMapCounts, int cnobj, int objcode, int* p1atm,
                                      int2* p2atm, int4* p4atm, int* atomEndPts,
                                      int* objTerminates, int creditTermini)
{
  int i;
  int ispresent[5], isterminal[5];

  // Get the number of atoms that would be added to the group if this object were included
  int2 score = FindObjectAtomsInUnit(unitID, unitMap, unitMapCounts, cnobj, objcode, p1atm,
                                     p2atm, p4atm, ispresent);

  // This is not a candidate if it would add more atoms than the group can currently support.
  if (score.y + bw->natom > BOND_WORK_UNIT_THREADS_PER_BLOCK) {
    return score.y + 1000;
  }

  // Legitimate candidate if we're still here.  If the object will add a
  // terminal atom that is NOT already in the group, credit the score.
  if (creditTermini != 0 && objTerminates[cnobj] == 1) {
    for (i = 0; i < 5; i++) {
      isterminal[i] = 0;
    }
    if (objcode == CNST_CODE) {
      isterminal[0] = atomEndPts[p1atm[cnobj]];
    }
    if (objcode == BOND_CODE || objcode == ANGL_CODE || objcode == QQXC_CODE ||
        objcode == NB14_CODE || objcode == NMR2_CODE || objcode == UREY_CODE) {
      isterminal[0] = atomEndPts[p2atm[cnobj].x];
      isterminal[1] = atomEndPts[p2atm[cnobj].y];
    }
    if (objcode == ANGL_CODE || objcode == NMR3_CODE) {
      isterminal[2] = atomEndPts[p1atm[cnobj]];
    }
    if (objcode == DIHE_CODE || objcode == CMAP_CODE || objcode == NMR4_CODE ||
        objcode == CIMP_CODE) {
      isterminal[0] = atomEndPts[p4atm[cnobj].x];
      isterminal[1] = atomEndPts[p4atm[cnobj].y];
      isterminal[2] = atomEndPts[p4atm[cnobj].z];
      isterminal[3] = atomEndPts[p4atm[cnobj].w];
    }
    if (objcode == CMAP_CODE) {
      isterminal[4] = atomEndPts[p1atm[cnobj]];
    }
    for (i = 0; i < 5; i++) {
      if (isterminal[i] == 1 && ispresent[i] == 0) {
        score.y -= creditTermini;
      }
    }
  }

  return score.y;
}

//---------------------------------------------------------------------------------------------
// AddObjectToUnit: add an object to a unit.
//
// Arguments:
//   bw:             the bond work unit to which the object will be added
//   unitID:         the ID number of the bond work unit in the alrger bwunits array
//   unitMap:        map of all the units that each atom is a part of
//   unitMapCounts:  the total number of work units in which each atom is included
//   objID:          ID number of the object to add
//   objcode:        code number for the object type
//   p1atm:          array of single atoms included in certain object types (i.e. all atom Ks
//                   in bond angles)
//   p2atm:          array of two-atom tuples included in certain object types (i.e. both atoms
//                   in any given bond)
//   p4atm:          array of four-atom tuples included in certain object types (i.e. all four
//                   atoms in a dihedral, or atoms IJKL of a CMAP term)
//   objAssigned:    table indicating which objects have been assigned to work units
//   nassigned:      counter of all the objects that have been assigned
//---------------------------------------------------------------------------------------------
static void AddObjectToUnit(bondwork *bw, int unitID, imat *unitMap, imat *unitMapCounts,
                            int objID, int objcode, int* p1atm, int2* p2atm,
                            int4* p4atm, int* objAssigned, int *nassigned)
{
  // Update the atom list and unit map
  int nunitatom = bw->natom;
  int ispresent[5];
  int *itmp;
  FindObjectAtomsInUnit(unitID, unitMap, unitMapCounts, objID, objcode, p1atm, p2atm, p4atm,
                        ispresent);
  if (objcode == CNST_CODE) {
    if (ispresent[0] == 0) {
      bw->atomList[nunitatom] = p1atm[objID];
      nunitatom++;
    }
    if (unitMapCounts == NULL) {
      itmp = unitMap->data;
      itmp[p1atm[objID]] = unitID;
    }
    else {
      unitMap->map[p1atm[objID]][unitMapCounts->data[p1atm[objID]]] = unitID;
      unitMapCounts->data[p1atm[objID]] += 1;
    }
  }
  if (objcode == BOND_CODE || objcode == QQXC_CODE || objcode == NB14_CODE ||
      objcode == NMR2_CODE || objcode == ANGL_CODE || objcode == NMR3_CODE ||
      objcode == UREY_CODE) {
    if (ispresent[0] == 0) {
      bw->atomList[nunitatom] = p2atm[objID].x;
      nunitatom++;
    }
    if (ispresent[1] == 0) {
      bw->atomList[nunitatom] = p2atm[objID].y;
      nunitatom++;
    }
    if (unitMapCounts == NULL) {
      itmp = unitMap->data;
      itmp[p2atm[objID].x] = unitID;
      itmp[p2atm[objID].y] = unitID;
    }
    else {
      unitMap->map[p2atm[objID].x][unitMapCounts->data[p2atm[objID].x]] = unitID;
      unitMapCounts->data[p2atm[objID].x] += 1;
      unitMap->map[p2atm[objID].y][unitMapCounts->data[p2atm[objID].y]] = unitID;
      unitMapCounts->data[p2atm[objID].y] += 1;
    }
  }
  if (objcode == DIHE_CODE || objcode == CMAP_CODE || objcode == NMR4_CODE ||
      objcode == CIMP_CODE) {
    if (ispresent[0] == 0) {
      bw->atomList[nunitatom] = p4atm[objID].x;
      nunitatom++;
    }
    if (ispresent[1] == 0) {
      bw->atomList[nunitatom] = p4atm[objID].y;
      nunitatom++;
    }
    if (ispresent[2] == 0) {
      bw->atomList[nunitatom] = p4atm[objID].z;
      nunitatom++;
    }
    if (ispresent[3] == 0) {
      bw->atomList[nunitatom] = p4atm[objID].w;
      nunitatom++;
    }
    if (unitMapCounts == NULL) {
      itmp = unitMap->data;
      itmp[p4atm[objID].x] = unitID;
      itmp[p4atm[objID].y] = unitID;
      itmp[p4atm[objID].z] = unitID;
      itmp[p4atm[objID].w] = unitID;
    }
    else {
      unitMap->map[p4atm[objID].x][unitMapCounts->data[p4atm[objID].x]] = unitID;
      unitMapCounts->data[p4atm[objID].x] += 1;
      unitMap->map[p4atm[objID].y][unitMapCounts->data[p4atm[objID].y]] = unitID;
      unitMapCounts->data[p4atm[objID].y] += 1;
      unitMap->map[p4atm[objID].z][unitMapCounts->data[p4atm[objID].z]] = unitID;
      unitMapCounts->data[p4atm[objID].z] += 1;
      unitMap->map[p4atm[objID].w][unitMapCounts->data[p4atm[objID].w]] = unitID;
      unitMapCounts->data[p4atm[objID].w] += 1;
    }
  }
  if (objcode == ANGL_CODE || objcode == NMR3_CODE || objcode == CMAP_CODE) {
    if (((objcode == ANGL_CODE || objcode == NMR3_CODE) && ispresent[2] == 0) ||
        (objcode == CMAP_CODE && ispresent[4] == 0)) {
      bw->atomList[nunitatom] = p1atm[objID];
      nunitatom++;
    }
    if (unitMapCounts == NULL) {
      itmp = unitMap->data;
      itmp[p1atm[objID]] = unitID;
    }
    else {
      unitMap->map[p1atm[objID]][unitMapCounts->data[p1atm[objID]]] = unitID;
      unitMapCounts->data[p1atm[objID]] += 1;
    }
  }
  bw->natom = nunitatom;

  // Mark the object as having been assigned
  objAssigned[objID] = 1;
  *nassigned += 1;

  // Add the object ID to the appropriate list
  if (objcode == BOND_CODE) {
    bw->bondList[bw->nbond] = objID;
    bw->nbond += 1;
  }
  else if (objcode == ANGL_CODE) {
    bw->bondAngleList[bw->nangl] = objID;
    bw->nangl += 1;
  }
  else if (objcode == DIHE_CODE) {
    bw->dihedralList[bw->ndihe] = objID;
    bw->ndihe += 1;
  }
  else if (objcode == CMAP_CODE) {
    bw->cmapList[bw->ncmap] = objID;
    bw->ncmap += 1;
  }
  else if (objcode == QQXC_CODE) {
    bw->qqxcList[bw->nqqxc] = objID;
    bw->nqqxc += 1;
  }
  else if (objcode == NB14_CODE) {
    bw->nb14List[bw->nnb14] = objID;
    bw->nnb14 += 1;
  }
  else if (objcode == NMR2_CODE) {
    bw->nmr2List[bw->nnmr2] = objID;
    bw->nnmr2 += 1;
  }
  else if (objcode == NMR3_CODE) {
    bw->nmr3List[bw->nnmr3] = objID;
    bw->nnmr3 += 1;
  }
  else if (objcode == NMR4_CODE) {
    bw->nmr4List[bw->nnmr4] = objID;
    bw->nnmr4 += 1;
  }
  else if (objcode == UREY_CODE) {
    bw->ureyList[bw->nurey] = objID;
    bw->nurey += 1;
  }
  else if (objcode == CIMP_CODE) {
    bw->cimpList[bw->ncimp] = objID;
    bw->ncimp += 1;
  }
  else if (objcode == CNST_CODE) {
    bw->cnstList[bw->ncnst] = objID;
    bw->ncnst += 1;
  }
}

//---------------------------------------------------------------------------------------------
// TraceBranch: function for tracing interconnected objects and making groups of them that
//              share as many atoms as possible
//
// Arguments:
//   objMap:         integer table mapping out all the objects in which any atom plays a part
//   objMapCounts:   integer vector giving the length of relevant data in each row of objMap
//   objAssigned:    table indicating whether each object has been assigned to a group
//   objTerminates:  table indicating whether each object contains a terminal atom
//   atomEndPts:     table to indciate whether each atom is an end point, bonded to exactly
//                   one other atom in a larger structure
//   tobj:           index of the object that seeds the branching
//   objcode:        numerical translation of the object type (2 = bonds, 3 = angles, etc.)
//   nassigned:      the number of objects currently assigned to groups
//   minFreeID:      the index of the first object that is yet unassigned
//   p[124]atm:      arrays giving the atom indices of atoms referenced by each object of
//                   type objcode in the simulation (i.e. p2atm holds i and j atoms for all
//                   bonds)
//   nobj:           total number of such objects in the simulation, assigned or not
//   scores:         scratch space for the keeping scores of each candidate object to add to
//                   the growing work unit
//   bw:             the bond work unit being built
//---------------------------------------------------------------------------------------------
static void TraceBranch(imat *objMap, int* objMapCounts, int* objAssigned, int* objTerminates,
                        int* atomEndPts, imat *unitMap, int tobj, int objcode, int *nassigned,
                        int *minFreeID, int *minFreeTermID, int* p1atm, int2* p2atm,
                        int4* p4atm, int nobj, imat *scores, bondwork *bw, int unitID)
{
  int i, j, k, ncandidate;
  int earmarks[BOND_WORK_UNIT_THREADS_PER_BLOCK];

  // The group begins with the atoms of the seed object.  Counts of atoms, bonds,
  // and other objects were initialized when the work unit was allocated.
  AddObjectToUnit(bw, unitID, unitMap, NULL, tobj, objcode, p1atm, p2atm, p4atm, objAssigned,
                  nassigned);

  // Permit up to 128 non-bonded electrostatic exclusions in a single
  // bond work unit.  Any other interactions have a limit of 128.
  int maxobj = BOND_WORK_UNIT_THREADS_PER_BLOCK;
  if (objcode == QQXC_CODE || objcode == NB14_CODE) {
    maxobj *= 8;
  }
  else if (objcode == DIHE_CODE) {
    maxobj *= 2;
  }

  // Repeat until all objects are exhausted, or the group fills up in terms
  // of either atoms or objects.  The number of objects in the unit is tracked
  // separately so that this function has a convenient metric without having
  // to figure out which of the work unit's counters is relevant.
  int nobjInUnit = 1;
  int objsize = GetWorkUnitObjectSize(objcode);
  int mfID = *minFreeID;
  int mftID = *minFreeTermID;
  bool filled = false;
  while (filled == false && nobjInUnit < maxobj && *nassigned < nobj) {

    // Find all objects yet unmapped that touch atoms of the current group
    ncandidate = 0;
    int nunitatom = bw->natom;
    for (i = 0; i < nunitatom; i++) {
      int iatm = bw->atomList[i];
      for (j = 0; j < objMapCounts[iatm]; j++) {
        int cnobj = objMap->map[iatm][j];
        if (objAssigned[cnobj] == 0 && ncandidate < BOND_WORK_UNIT_THREADS_PER_BLOCK) {

          // Compute a score based on the number of atoms that would
          // be added to the work unit if this object were included.
          scores->map[0][ncandidate] = cnobj;
          scores->map[1][ncandidate] = ScoreObjectPlacementInUnit(bw, unitID, unitMap, NULL,
                                                                  cnobj, objcode, p1atm, p2atm,
                                                                  p4atm, atomEndPts,
                                                                  objTerminates, 0);

          // Increment the number of candidates only if the score tells us that the
          // new object shares at least one atom with the rest of the growing unit.
          // Temporarily give this object a different assignment to avoid including
          // it twice in the list of candidates.
          if (scores->map[1][ncandidate] < objsize) {
            objAssigned[cnobj] = 2;
            earmarks[ncandidate] = cnobj;
            ncandidate++;
          }
        }
      }
    }

    // Revert earmarked objects
    for (i = 0; i < ncandidate; i++) {
      objAssigned[earmarks[i]] = 0;
    }

    // If there have not been any objects found thus far, there must be at least
    // enough space in the unit to add a completely new one to the group.  If not,
    // we're done buidling this group.
    if (ncandidate == 0 && nunitatom + objsize > BOND_WORK_UNIT_THREADS_PER_BLOCK) {
      filled = true;
      continue;
    }

    // If there are no candidates, we need to start searching other objects
    // that may not be connected to the current group of atoms.  Same deal:
    // give things with terminal atoms priority.
    if (ncandidate == 0) {
      for (i = mftID; i < nobj; i++) {
        mftID += (objAssigned[mftID] == 1 || objTerminates[mftID] == 0);
        if (objAssigned[i] == 0 && objTerminates[i] == 1) {
          scores->map[0][0] = i;
          scores->map[1][0] = ScoreObjectPlacementInUnit(bw, unitID, unitMap, NULL, i, objcode,
                                                         p1atm, p2atm, p4atm, atomEndPts,
                                                         objTerminates, 1);
          if (scores->map[1][0] <= objsize) {
            ncandidate = 1;
            break;
          }
        }
      }
    }
    if (ncandidate == 0) {
      for (i = mfID; i < nobj; i++) {
        mfID += (objAssigned[mfID] == 1);
        if (objAssigned[i] == 0) {
          scores->map[0][0] = i;
          scores->map[1][0] = ScoreObjectPlacementInUnit(bw, unitID, unitMap, NULL, i, objcode,
                                                         p1atm, p2atm, p4atm, atomEndPts,
                                                         objTerminates, 1);
          if (scores->map[1][0] <= objsize) {
            ncandidate = 1;
            break;
          }
        }
      }
    }

    // No suitable candidates?  We're done building this group.
    filled = (ncandidate == 0);

    // Select candidates with the lowest score (the score is given by the number
    // of atoms in the dihedral candidate not already found in other dihedrals of
    // the group, less the number of terminal atoms in the candidate).
    if (ncandidate > 0) {
      int bscore = scores->map[1][0];
      for (i = 1; i < ncandidate; i++) {
        if (scores->map[1][i] < bscore) {
          bscore = scores->map[1][i];
        }
      }

      // Add the best candidate, and its atoms, to the group
      i = 0;
      int nadded = 0;
      while (i < ncandidate && nobjInUnit < maxobj && *nassigned < nobj &&
             (nadded == 0 || bw->natom < BOND_WORK_UNIT_THREADS_PER_BLOCK - objsize)) {
        if (scores->map[1][i] == bscore) {
          AddObjectToUnit(bw, unitID, unitMap, NULL, scores->map[0][i], objcode, p1atm, p2atm,
                          p4atm, objAssigned, nassigned);

          // Increment the number of objects in the group,
          // and the number of objects thus far assigned
          nobjInUnit++;
          nadded++;
        }
        i++;
      }
    }
  }

  // These local variables were assigned and incremented to avoid constant de-referencing.
  *minFreeID = mfID;
  *minFreeTermID = mftID;
}

//---------------------------------------------------------------------------------------------
// PartitionIntegerArray: function for setting a pointer to a segment of an integer data array
//                        based on an input position, then incrementing the input position to
//                        reflect the amount of data that has bene spoken for.
//
// Arguments:
//   V:         the data array (integer type)
//   pos:       input position to start the new segment (also returned)
//   len:       length of the segment (strictly input)
//---------------------------------------------------------------------------------------------
static int* PartitionIntegerArray(int* V, int *pos, int len)
{
  int *iptr;

  iptr = &V[*pos];
  *pos += len;

  return iptr;
}

//---------------------------------------------------------------------------------------------
// CreateBondWorkUnit: create a bonded interaction work unit--that is, allocated memory for
//                     all aspects of its data, and initialize the counts of each type of
//                     object it might contain to zero.
//
// Arguments:
//   ti_mode:   passed in from the eponymous attribute of a cudaSimulation (as found in host
//              memory).  This is the Thermodynamic Integration mode
//---------------------------------------------------------------------------------------------
static bondwork CreateBondWorkUnit(int ti_mode)
{
  bondwork bw;

  // Initialize the numbers of each object
  bw.natom = 0;
  bw.nbond = 0;
  bw.nqqxc = 0;
  bw.nnb14 = 0;
  bw.nangl = 0;
  bw.ndihe = 0;
  bw.ncmap = 0;
  bw.nnmr2 = 0;
  bw.nnmr3 = 0;
  bw.nnmr4 = 0;
  bw.nurey = 0;
  bw.ncimp = 0;
  bw.ncnst = 0;

  // Set the accumulator to the bonded force array by default
  bw.frcAcc = BOND_FORCE_ACCUMULATOR;

  // Allocate data array and set pointers into it.  Start with lists of atoms and objects.
  int nsets = 94 + 28*(ti_mode > 0);
  bw.data = (int*)malloc(BOND_WORK_UNIT_THREADS_PER_BLOCK * nsets * sizeof(int));
  int pos = 0;
  bw.atomList      = PartitionIntegerArray(bw.data, &pos,    BOND_WORK_UNIT_THREADS_PER_BLOCK);
  bw.bondList      = PartitionIntegerArray(bw.data, &pos,    BOND_WORK_UNIT_THREADS_PER_BLOCK);
  bw.qqxcList      = PartitionIntegerArray(bw.data, &pos,  8*BOND_WORK_UNIT_THREADS_PER_BLOCK);
  bw.nb14List      = PartitionIntegerArray(bw.data, &pos,  8*BOND_WORK_UNIT_THREADS_PER_BLOCK);
  bw.bondAngleList = PartitionIntegerArray(bw.data, &pos,    BOND_WORK_UNIT_THREADS_PER_BLOCK);
  bw.dihedralList  = PartitionIntegerArray(bw.data, &pos,  2*BOND_WORK_UNIT_THREADS_PER_BLOCK);
  bw.cmapList      = PartitionIntegerArray(bw.data, &pos,    BOND_WORK_UNIT_THREADS_PER_BLOCK);
  bw.nmr2List      = PartitionIntegerArray(bw.data, &pos,    BOND_WORK_UNIT_THREADS_PER_BLOCK);
  bw.nmr3List      = PartitionIntegerArray(bw.data, &pos,    BOND_WORK_UNIT_THREADS_PER_BLOCK);
  bw.nmr4List      = PartitionIntegerArray(bw.data, &pos,    BOND_WORK_UNIT_THREADS_PER_BLOCK);
  bw.ureyList      = PartitionIntegerArray(bw.data, &pos,    BOND_WORK_UNIT_THREADS_PER_BLOCK);
  bw.cimpList      = PartitionIntegerArray(bw.data, &pos,    BOND_WORK_UNIT_THREADS_PER_BLOCK);
  bw.cnstList      = PartitionIntegerArray(bw.data, &pos,    BOND_WORK_UNIT_THREADS_PER_BLOCK);

  // Space to store how each object maps into the imported atoms table for this work unit
  bw.bondMap       = PartitionIntegerArray(bw.data, &pos,  2*BOND_WORK_UNIT_THREADS_PER_BLOCK);
  bw.qqxcMap       = PartitionIntegerArray(bw.data, &pos, 16*BOND_WORK_UNIT_THREADS_PER_BLOCK);
  bw.nb14Map       = PartitionIntegerArray(bw.data, &pos, 16*BOND_WORK_UNIT_THREADS_PER_BLOCK);
  bw.bondAngleMap  = PartitionIntegerArray(bw.data, &pos,  3*BOND_WORK_UNIT_THREADS_PER_BLOCK);
  bw.dihedralMap   = PartitionIntegerArray(bw.data, &pos,  8*BOND_WORK_UNIT_THREADS_PER_BLOCK);
  bw.cmapMap       = PartitionIntegerArray(bw.data, &pos,  5*BOND_WORK_UNIT_THREADS_PER_BLOCK);
  bw.nmr2Map       = PartitionIntegerArray(bw.data, &pos,  2*BOND_WORK_UNIT_THREADS_PER_BLOCK);
  bw.nmr3Map       = PartitionIntegerArray(bw.data, &pos,  3*BOND_WORK_UNIT_THREADS_PER_BLOCK);
  bw.nmr4Map       = PartitionIntegerArray(bw.data, &pos,  4*BOND_WORK_UNIT_THREADS_PER_BLOCK);
  bw.ureyMap       = PartitionIntegerArray(bw.data, &pos,  2*BOND_WORK_UNIT_THREADS_PER_BLOCK);
  bw.cimpMap       = PartitionIntegerArray(bw.data, &pos,  4*BOND_WORK_UNIT_THREADS_PER_BLOCK);
  bw.cnstMap       = PartitionIntegerArray(bw.data, &pos,    BOND_WORK_UNIT_THREADS_PER_BLOCK);

  // Space to store the identities of each object vis-a-vis TI
  if (ti_mode > 0) {
    bw.bondStatus  = PartitionIntegerArray(bw.data, &pos,    BOND_WORK_UNIT_THREADS_PER_BLOCK);
    bw.anglStatus  = PartitionIntegerArray(bw.data, &pos,    BOND_WORK_UNIT_THREADS_PER_BLOCK);
    bw.diheStatus  = PartitionIntegerArray(bw.data, &pos,  2*BOND_WORK_UNIT_THREADS_PER_BLOCK);
    bw.cmapStatus  = PartitionIntegerArray(bw.data, &pos,    BOND_WORK_UNIT_THREADS_PER_BLOCK);
    bw.qqxcStatus  = PartitionIntegerArray(bw.data, &pos,  8*BOND_WORK_UNIT_THREADS_PER_BLOCK);
    bw.nb14Status  = PartitionIntegerArray(bw.data, &pos,  8*BOND_WORK_UNIT_THREADS_PER_BLOCK);
    bw.nmr2Status  = PartitionIntegerArray(bw.data, &pos,    BOND_WORK_UNIT_THREADS_PER_BLOCK);
    bw.nmr3Status  = PartitionIntegerArray(bw.data, &pos,    BOND_WORK_UNIT_THREADS_PER_BLOCK);
    bw.nmr4Status  = PartitionIntegerArray(bw.data, &pos,    BOND_WORK_UNIT_THREADS_PER_BLOCK);
    bw.ureyStatus  = PartitionIntegerArray(bw.data, &pos,    BOND_WORK_UNIT_THREADS_PER_BLOCK);
    bw.cimpStatus  = PartitionIntegerArray(bw.data, &pos,    BOND_WORK_UNIT_THREADS_PER_BLOCK);
    bw.cnstStatus  = PartitionIntegerArray(bw.data, &pos,    BOND_WORK_UNIT_THREADS_PER_BLOCK);
  }

  return bw;
}

//---------------------------------------------------------------------------------------------
// DestroyBondWorkUnit: destructor for the bond work unit data struct.  This is not C++
//                      because of how dynamic I want the bond work units array to be.
//
// Arguments:
//   bw:      the bond work unit to destroy
//---------------------------------------------------------------------------------------------
void DestroyBondWorkUnit(bondwork *bw)
{
  if(bw) if(bw->data) free(bw->data);
}

//---------------------------------------------------------------------------------------------
// GetTerminalObjects: get a table of objects containing terminal atoms.
//
// Arguments:
//   nobj:        the number of objects listed in p1atm, p2atm, and p4atm
//   p1atm:
//   p2atm:       arrays detailing the atom IDs of particles in a list of interactions
//   p4atm:
//   atomEndPts:  array of atoms already known to bond only to one other atom in their
//                respective structures
//---------------------------------------------------------------------------------------------
static imat GetTerminalObjects(int nobj, int* p1atm, int2* p2atm, int4* p4atm, int* atomEndPts)
{
  int i;
  imat objTerminates;

  objTerminates = CreateImat(1, nobj);
  if (p1atm != NULL) {
    for (i = 0; i < nobj; i++) {
      objTerminates.data[i] = atomEndPts[p1atm[i]];
    }
  }
  if (p2atm != NULL) {
    for (i = 0; i < nobj; i++) {
      objTerminates.data[i] = (objTerminates.data[i] ||
                               atomEndPts[p2atm[i].x] == 1 || atomEndPts[p2atm[i].y] == 1);
    }
  }
  if (p4atm != NULL) {
    for (i = 0; i < nobj; i++) {
      objTerminates.data[i] = (objTerminates.data[i] ||
                               atomEndPts[p4atm[i].x] == 1 || atomEndPts[p4atm[i].y] == 1 ||
                               atomEndPts[p4atm[i].z] == 1 || atomEndPts[p4atm[i].w] == 1);
    }
  }

  return objTerminates;
}


//---------------------------------------------------------------------------------------------
// GetSimulationObjectCount: get a count of the number of objects based on the object code or
//                           object type
//
// Arguments:
//   gpu:       overarching data structure containing simulation information, here used for
//              object counts
//   objcode:   code for the object to seek out--this, or objtype, must be supplied
//   objtype:   type name of the object to search for
//---------------------------------------------------------------------------------------------
static int GetSimulationObjectCount(gpuContext gpu, int objcode, const char* objtype)
{
  // Work off of the type name, if it is supplied
  if (strcmp(objtype, "BOND") == 0 || objcode == BOND_CODE) {
    return gpu->sim.bonds;
  }
  else if (strcmp(objtype, "ANGLE") == 0 || objcode == ANGL_CODE) {
    return gpu->sim.bondAngles;
  }
  else if (strcmp(objtype, "DIHEDRAL") == 0 || objcode == DIHE_CODE) {
    return gpu->sim.dihedrals;
  }
  else if (strcmp(objtype, "CMAP") == 0 || objcode == CMAP_CODE) {
    return gpu->sim.cmaps;
  }
  else if (strcmp(objtype, "QQXC") == 0 || objcode == QQXC_CODE) {
    return gpu->sim.qqxcs;
  }
  else if (strcmp(objtype, "NB14") == 0 || objcode == NB14_CODE) {
    return gpu->sim.nb14s;
  }
  else if (strcmp(objtype, "NMR2") == 0 || objcode == NMR2_CODE) {
    return gpu->sim.NMRDistances;
  }
  else if (strcmp(objtype, "NMR3") == 0 || objcode == NMR3_CODE) {
    return gpu->sim.NMRAngles;
  }
  else if (strcmp(objtype, "NMR4") == 0 || objcode == NMR4_CODE) {
    return gpu->sim.NMRTorsions;
  }
  else if (strcmp(objtype, "UREY") == 0 || objcode == UREY_CODE) {
    return gpu->sim.UBAngles;
  }
  else if (strcmp(objtype, "CIMP") == 0 || objcode == CIMP_CODE) {
    return gpu->sim.impDihedrals;
  }
  else if (strcmp(objtype, "CNST") == 0 || objcode == CNST_CODE) {
    return gpu->sim.constraints;
  }

  // Return 0 if nothing else
  return 0;
}

//---------------------------------------------------------------------------------------------
// GetWorkUnitObjectCount: get the number of objects in a work unit based on the object code.
//
// Arguments:
//   bw:       the work unit
//   objcode:  the object code (2 = bonds, 4 = dihedrals, etc.)
//---------------------------------------------------------------------------------------------
static int GetWorkUnitObjectCount(bondwork *bw, int objcode)
{
  if (objcode == BOND_CODE) {
    return bw->nbond;
  }
  else if (objcode == ANGL_CODE) {
    return bw->nangl;
  }
  else if (objcode == DIHE_CODE) {
    return bw->ndihe;
  }
  else if (objcode == CMAP_CODE) {
    return bw->ncmap;
  }
  else if (objcode == QQXC_CODE) {
    return bw->nqqxc;
  }
  else if (objcode == NB14_CODE) {
    return bw->nnb14;
  }
  else if (objcode == NMR2_CODE) {
    return bw->nnmr2;
  }
  else if (objcode == NMR3_CODE) {
    return bw->nnmr3;
  }
  else if (objcode == NMR4_CODE) {
    return bw->nnmr4;
  }
  else if (objcode == UREY_CODE) {
    return bw->nurey;
  }
  else if (objcode == CIMP_CODE) {
    return bw->ncimp;
  }
  else if (objcode == CNST_CODE) {
    return bw->ncnst;
  }

  // Return 0 if nothing else
  return 0;
}

//---------------------------------------------------------------------------------------------
// GetObjectMaps: create tables showing how objects of this type map to system atoms and
//                vice-versa.
//
// Arguments:
//   gpu:           data structure containing bonded interactions and their atom indices
//   objtype:       name specifying the object type (i.e. "DIHEDRALS")
//   atomEndPts:    table indicating which atoms are bonded to one and only one other atom
//   objMapCounts:  (returned) number of objects of type objtype that each atom is found in
//   objMap:        (returned) objects of type objtype that each atom is found in
//   objTerminates: (returned) table indicating which objects contain terminal atoms
//---------------------------------------------------------------------------------------------
static void GetObjectMaps(gpuContext gpu, const char* objtype, imat *atomEndPts,
                          imat *objMapCounts, imat *objMap, imat *objTerminates)
{
  int i, j, objcode=0;

  // Return immediately if there are no objects to map--don't
  // get caught mucking around with unassigned pointers.
  int nobj = GetSimulationObjectCount(gpu, 0, objtype);
  if (nobj == 0) {
    return;
  }

  // Set pointers and codify the object type
  int *p1atm = NULL;
  int2 *p2atm = NULL;
  int4 *p4atm = NULL;
  codify_object_type(objtype, gpu, p1atm, p2atm, p4atm, objcode);

  // Get the number of objects each atom is a part of
  *objMapCounts = CreateImat(1, gpu->sim.atoms);
  for (i = 0; i < nobj; i++) {
    if (p1atm != NULL) {
      objMapCounts->data[p1atm[i]] += 1;
    }
    if (p2atm != NULL) {
      objMapCounts->data[p2atm[i].x] += 1;
      objMapCounts->data[p2atm[i].y] += 1;
    }
    if (p4atm != NULL) {
      objMapCounts->data[p4atm[i].x] += 1;
      objMapCounts->data[p4atm[i].y] += 1;
      objMapCounts->data[p4atm[i].z] += 1;
      objMapCounts->data[p4atm[i].w] += 1;
    }
  }
  int maxincl = IVecExtreme(objMapCounts->data, objMapCounts->col);
  *objMap = CreateImat(objMapCounts->col, maxincl);
  SetIVec(objMapCounts->data, objMapCounts->col, 0);
  for (i = 0; i < nobj; i++) {
    if (p1atm != NULL) {
      j = p1atm[i];
      objMap->map[j][objMapCounts->data[j]] = i;
      objMapCounts->data[j] += 1;
    }
    if (p2atm != NULL) {
      j = p2atm[i].x;
      objMap->map[j][objMapCounts->data[j]] = i;
      objMapCounts->data[j] += 1;
      j = p2atm[i].y;
      objMap->map[j][objMapCounts->data[j]] = i;
      objMapCounts->data[j] += 1;
    }
    if (p4atm != NULL) {
      j = p4atm[i].x;
      objMap->map[j][objMapCounts->data[j]] = i;
      objMapCounts->data[j] += 1;
      j = p4atm[i].y;
      objMap->map[j][objMapCounts->data[j]] = i;
      objMapCounts->data[j] += 1;
      j = p4atm[i].z;
      objMap->map[j][objMapCounts->data[j]] = i;
      objMapCounts->data[j] += 1;
      j = p4atm[i].w;
      objMap->map[j][objMapCounts->data[j]] = i;
      objMapCounts->data[j] += 1;
    }
  }

  // Get a table of whether each object contains terminal atoms
  *objTerminates = GetTerminalObjects(nobj, p1atm, p2atm, p4atm, atomEndPts->data);
}

//---------------------------------------------------------------------------------------------
// MaskBlankObjects: create a mask (by marking the objects as already assigned) to prevent
//                   objects with zero force constants from being included in the work units.
//                   Only certain objects will be masked: for example, the NMR terms are not
//                   masked because they have to be set explicitly in the input file and we do
//                   not assume that anyone would inject large numbers of blank NMR restraints.
//                   Electrostatic exclusions are already masked by the routine that created
//                   them.
//
// Arguments:
//   gpu:          overarching data structure containing all simulation parameters and terms
//   objcode:      the object code, obtained from the codify_object_type macro above
//   objAssigned:  array that will take the mask
//---------------------------------------------------------------------------------------------
static void MaskBlankObjects(gpuContext gpu, int objcode, imat *objAssigned)
{
  int i;

  int nobj = GetSimulationObjectCount(gpu, objcode, " ");
  if (objcode == BOND_CODE) {
    for (i = 0; i < nobj; i++) {
      if (fabs(gpu->pbBond->_pSysData[i].x) < 1.0e-8) {
        objAssigned->data[i] = 1;
      }
    }
  }
  else if (objcode == ANGL_CODE) {
    for (i = 0; i < nobj; i++) {
      if (fabs(gpu->pbBondAngle->_pSysData[i].x) < 1.0e-8) {
        objAssigned->data[i] = 1;
      }
    }
  }
  else if (objcode == DIHE_CODE) {
    for (i = 0; i < nobj; i++) {
      if (fabs(gpu->pbDihedral2->_pSysData[i].x) < 1.0e-8) {
        objAssigned->data[i] = 1;
      }
    }
  }
  else if (objcode == UREY_CODE) {
    for (i = 0; i < nobj; i++) {
      if (fabs(gpu->pbUBAngle->_pSysData[i].x) < 1.0e-8) {
        objAssigned->data[i] = 1;
      }
    }
  }
}

//---------------------------------------------------------------------------------------------
// LocateObjectAmongUnits: locate an object among all work units.
//
// Arguments:
//   bwunits:   the list of all bond work units
//   nbwunits:  the length of bwunits
//   objid:     ID number of the object to seek out
//   objcode:   code number of the object to seek out
//   msg:       message to print to help with traceback
//
// This is a debugging function.
//---------------------------------------------------------------------------------------------
static void LocateObjectAmongUnits(bondwork *bwunits, int nbwunits, int objid, int objcode,
                                   const char* msg)
{
  int i, j;

  int *listPtr;
  bool found = false;
  for (i = 0; i < nbwunits; i++) {
    if (objcode == BOND_CODE) {
      listPtr = bwunits[i].bondList;
    }
    else if (objcode == ANGL_CODE) {
      listPtr = bwunits[i].bondAngleList;
    }
    else if (objcode == DIHE_CODE) {
      listPtr = bwunits[i].dihedralList;
    }
    else if (objcode == CMAP_CODE) {
      listPtr = bwunits[i].cmapList;
    }
    else if (objcode == QQXC_CODE) {
      listPtr = bwunits[i].qqxcList;
    }
    else if (objcode == NB14_CODE) {
      listPtr = bwunits[i].nb14List;
    }
    else if (objcode == NMR2_CODE) {
      listPtr = bwunits[i].nmr2List;
    }
    else if (objcode == NMR3_CODE) {
      listPtr = bwunits[i].nmr3List;
    }
    else if (objcode == NMR4_CODE) {
      listPtr = bwunits[i].nmr4List;
    }
    else if (objcode == UREY_CODE) {
      listPtr = bwunits[i].ureyList;
    }
    else if (objcode == CIMP_CODE) {
      listPtr = bwunits[i].cimpList;
    }
    else if (objcode == CNST_CODE) {
      listPtr = bwunits[i].cnstList;
    }
    int nobj = GetWorkUnitObjectCount(&bwunits[i], objcode);
    for (j = 0; j < nobj; j++) {
      if (listPtr[j] == objid) {
        printf("%s :: Object %6d (code %3d) is found in unit %6d [%4d]\n", msg, objid,
               objcode, i, j);
        found = true;
      }
    }
  }
  if (found == false) {
    printf("%s :: Object %6d (code %3d) was not found.\n", msg, objid, objcode);
  }
}

//---------------------------------------------------------------------------------------------
// BaseNGroups: function to encapsulate the creation of bonded groups starting with objects of
//              N atoms.  For most contexts, N will imply the sort of object: bonds = 2,
//              angles = 3, dihedrals/impropers = 4, CMAP = 5.  Additional modifiers may be
//              added later.  The number of unique atoms in each group may not exceed 128--for
//              most systems this will still permit 128 dihedrals--not sure about CMAP terms.
//
// Arguments:
//   gpu:            overarching structure holding all simulation information, critically in
//                   this case bonded term atom IDs
//   objtype:        four-letter code for the type of object under consideration
//   atomEndPts:     array indicating which atoms are terminal atoms
//   objAssigned:    array indicating which objects have already found a work unit home
//   objmapCounts:   the number of objects currently assigned to each work unit
//   objTerminates:  array indicating whether each object contains one or more terminal atoms
//   unitMap:        map of the atoms needed by each work unit
//   bwunits:        array of bond work units stored in terms of CPU-accessible information
//   nbwunits:       the length of bwunits
//---------------------------------------------------------------------------------------------
static bondwork* BaseNGroups(gpuContext gpu, const char* objtype, imat *atomEndPts,
                             imat *objAssigned, imat* objMapCounts, imat *objMap,
                             imat *objTerminates, imat *unitMap, bondwork* bwunits,
                             int *nbwunits)
{
  int i, j, objcode=0;
  int lnbwunits = *nbwunits;

  // Return immediately if there are no objects to map--don't
  // get caught mucking around with unassigned pointers.
  int nobj = GetSimulationObjectCount(gpu, 0, objtype);
  if (nobj == 0) {
    return bwunits;
  }

  // Set pointers and codify the object type
  int *p1atm = NULL;
  int2 *p2atm = NULL;
  int4 *p4atm = NULL;
  codify_object_type(objtype, gpu, p1atm, p2atm, p4atm, objcode);

  // Permit up to 1024 non-bonded electrostatic exclusions in a single
  // bond work unit.  Any other interactions have a limit of 128.
  int maxobj = BOND_WORK_UNIT_THREADS_PER_BLOCK;
  if (objcode == QQXC_CODE || objcode == NB14_CODE) {
    maxobj *= 8;
  }
  else if (objcode == DIHE_CODE) {
    maxobj *= 2;
  }

  // Allocate new work units as necessary
  int nnew = ((nobj - ISum(objAssigned->data, nobj) + maxobj - 1) / maxobj) + 1;
  bwunits = (bondwork*)realloc(bwunits, (lnbwunits + nnew)*sizeof(bondwork));
  for (i = lnbwunits; i < lnbwunits + nnew; i++) {
    bwunits[i] = CreateBondWorkUnit(gpu->sim.ti_mode);
  }

  // Make assignments to groups until there are no more objects of this type to assign.
  int nassigned = ISum(objAssigned->data, objAssigned->col);
  int minFreeID = 0;
  int minFreeTermID = 0;
  int currunit = lnbwunits;
  lnbwunits += nnew;
  imat scores = CreateImat(2, BOND_WORK_UNIT_THREADS_PER_BLOCK);
  while (nassigned < nobj) {

    // Preferably, find an unassigned object that involves a terminal atom
    int tobj = -1;
    i = minFreeTermID;
    while (i < nobj && tobj == -1) {
      if (objAssigned->data[i] == 0 && objTerminates->data[i] == 1) {
        tobj = i;
      }
      i++;
    }

    // Choose the first available object
    if (tobj == -1) {
      i = minFreeID;
      while (i < nobj && tobj == -1) {
        if (objAssigned->data[i] == 0) {
          tobj = i;
        }
        i++;
      }
    }

    // Loop through other objects and expand the group as permissible
    TraceBranch(objMap, objMapCounts->data, objAssigned->data, objTerminates->data,
                atomEndPts->data, unitMap, tobj, objcode, &nassigned, &minFreeID,
                &minFreeTermID, p1atm, p2atm, p4atm, nobj, &scores, &bwunits[currunit],
                currunit);
    currunit++;

    // If the number of work units has been hit (some of them didn't fill up completely
    // with bonded interaction objects because the atom limits were reached), extend the
    // array.
    if (currunit == lnbwunits) {
      nnew = (lnbwunits > GRID * 2) ? (lnbwunits / 2) : GRID;
      bwunits = (bondwork*)realloc(bwunits, (lnbwunits + nnew)*sizeof(bondwork));
      for (i = lnbwunits; i < lnbwunits + nnew; i++) {
        bwunits[i] = CreateBondWorkUnit(gpu->sim.ti_mode);
      }
      lnbwunits += nnew;
    }
  }

  // Free allocated memory
  DestroyImat(&scores);

  // Update the number of bond work units (a local, non-pointer placeholder was in use)
  *nbwunits = lnbwunits;

  return bwunits;
}

//---------------------------------------------------------------------------------------------
// MergeIntegerLists: function for adding all unique integers in one list to another.  This is
//                    here for the time being but is really a much more general function.
//
// Arguments:
//   listA:       the list which will receive the union of both lists
//   nA:          on input, the length of listA--on output, the combined list length
//   listB:       the list to incorporate
//   nB:          the length of listB
//---------------------------------------------------------------------------------------------
static void MergeIntegerLists(int* listA, int *nA, int* listB, int nB)
{
  int i, j;

  int lnA = *nA;
  for (i = 0; i < nB; i++) {
    int redundant = 0;
    for (j = 0; j < lnA; j++) {
      redundant += (listA[j] == listB[i]);
    }
    if (redundant == 0) {
      listA[lnA] = listB[i];
      lnA++;
    }
  }
  *nA = lnA;
}

//---------------------------------------------------------------------------------------------
// ExtendWorkUnits: this function will add objects of a certain type to EXISTING bond work
//                  units under either of two circumstances: if the bonded interaction object
//                  involves one or more atoms that are already part of the work unit (also
//                  requiring that there be space in the work unit to hold any additional
//                  atoms that would be needed by the interaction), or if the work unit has
//                  enough space to hold all of the atoms and empty slots in a warp that
//                  could take the interaction.
//
// Arguments:
//   gpu:            overarching structure holding all simulation information, critically in
//                   this case bonded term atom IDs
//   objtype:        four-letter code for the type of object under consideration
//   atomEndPts:     array indicating which atoms are terminal atoms
//   objAssigned:    array indicating which objects have already found a work unit home
//   objmapCounts:   the number of objects currently assigned to each work unit
//   objMap:         map of the work units that any one object could be placed in, along with
//                   scores indicating which choice would be preferable
//   objTerminates:  array indicating whether each object contains one or more terminal atoms
//   unitMap:        map of the atoms needed by each work unit
//   unitMapCounts:  number of the atoms needed by each work unit
//   bwunits:        array of bond work units stored in terms of CPU-accessible information
//   nbwunits:       the length of bwunits
//---------------------------------------------------------------------------------------------
static void ExtendWorkUnits(gpuContext gpu, const char* objtype, imat *atomEndPts,
                            imat *objAssigned, imat *objMapCounts, imat *objMap,
                            imat *objTerminates, imat *unitMap, imat *unitMapCounts,
                            bondwork* bwunits, int nbwunits)
{
  int i, j, objcode=0;
  int ispresent[5];

  // Return immediately if there are no objects to map--don't
  // get caught mucking around with unassigned pointers.
  int nobj = GetSimulationObjectCount(gpu, -1, objtype);
  if (nobj == 0 || nbwunits == 0) {
    return;
  }

  // Set pointers and codify the object type
  int *p1atm = NULL;
  int2 *p2atm = NULL;
  int4 *p4atm = NULL;
  codify_object_type(objtype, gpu, p1atm, p2atm, p4atm, objcode);

  // Mask objects that will have no consequence on the system
  MaskBlankObjects(gpu, objcode, objAssigned);

  // Proceed one object at a time.  Make a map of all units that have
  // one or more atoms needed by this object.  Find the best fit.
  int objsize = GetWorkUnitObjectSize(objcode);
  int nassigned = ISum(objAssigned->data, nobj);
  imat unitFoothold;
  unitFoothold = CreateImat(2, nbwunits);

  // Permit up to 1024 non-bonded electrostatic exclusions in a single
  // bond work unit.  Any other interactions have a limit of 128.
  int maxobj = BOND_WORK_UNIT_THREADS_PER_BLOCK;
  if (objcode == QQXC_CODE || objcode == NB14_CODE) {
    maxobj *= 8;
  }
  else if (objcode == DIHE_CODE) {
    maxobj *= 2;
  }

  // Loop over all objects and make assignments
  for (i = 0; i < nobj; i++) {

    // Make a list of all work units this object could go in
    int ufhc = 0;
    if (objcode == CNST_CODE) {
      MergeIntegerLists(unitFoothold.map[0], &ufhc, unitMap->map[p1atm[i]],
                        unitMapCounts->data[p1atm[i]]);
    }
    if (objcode == BOND_CODE || objcode == ANGL_CODE || objcode == QQXC_CODE ||
        objcode == NB14_CODE || objcode == NMR2_CODE || objcode == NMR3_CODE ||
        objcode == UREY_CODE) {
      MergeIntegerLists(unitFoothold.map[0], &ufhc, unitMap->map[p2atm[i].x],
                        unitMapCounts->data[p2atm[i].x]);
      MergeIntegerLists(unitFoothold.map[0], &ufhc, unitMap->map[p2atm[i].y],
                        unitMapCounts->data[p2atm[i].y]);
    }
    if (objcode == DIHE_CODE || objcode == CMAP_CODE || objcode == NMR4_CODE ||
        objcode == CIMP_CODE) {
      MergeIntegerLists(unitFoothold.map[0], &ufhc, unitMap->map[p4atm[i].x],
                        unitMapCounts->data[p4atm[i].x]);
      MergeIntegerLists(unitFoothold.map[0], &ufhc, unitMap->map[p4atm[i].y],
                        unitMapCounts->data[p4atm[i].y]);
      MergeIntegerLists(unitFoothold.map[0], &ufhc, unitMap->map[p4atm[i].z],
                        unitMapCounts->data[p4atm[i].z]);
      MergeIntegerLists(unitFoothold.map[0], &ufhc, unitMap->map[p4atm[i].w],
                        unitMapCounts->data[p4atm[i].w]);
    }
    if (objcode == ANGL_CODE || objcode == NMR3_CODE || objcode == CMAP_CODE) {
      MergeIntegerLists(unitFoothold.map[0], &ufhc, unitMap->map[p1atm[i]],
                        unitMapCounts->data[p1atm[i]]);
    }
    for (j = 0; j < ufhc; j++) {
      unitFoothold.map[1][j] = ScoreObjectPlacementInUnit(&bwunits[unitFoothold.map[0][j]],
                                                          unitFoothold.map[0][j], unitMap,
                                                          unitMapCounts, i, objcode, p1atm,
                                                          p2atm, p4atm, atomEndPts->data,
                                                          objTerminates->data, 1);
      int ncurrobj = GetWorkUnitObjectCount(&bwunits[unitFoothold.map[0][j]], objcode);
      int warpsnap = ((ncurrobj & GRID_BITS_MASK) == 0) * 10;
      int unitfull = (ncurrobj >= maxobj) * 1000;
      unitFoothold.map[1][j] += warpsnap + unitfull;
    }

    // Check that there is a valid placement in at least one of the units
    int valid = 0;
    for (j = 0; j < ufhc; j++) {
      valid += (unitFoothold.map[1][j] <= objsize + 10);
    }

    // No candidates?  Just move on.
    if (ufhc == 0 || valid == 0) {
      continue;
    }

    // Select the best home for this object
    int bestunit = unitFoothold.map[0][0];
    int bestscore = unitFoothold.map[1][0];
    for (j = 0; j < ufhc; j++) {
      if (unitFoothold.map[1][j] < bestscore) {
        bestunit = unitFoothold.map[0][j];
        bestscore = unitFoothold.map[1][j];
      }
    }

    // Add the object to the work unit
    AddObjectToUnit(&bwunits[bestunit], bestunit, unitMap, unitMapCounts, i, objcode, p1atm,
                    p2atm, p4atm, objAssigned->data, &nassigned);
  }

  // Free allocated memory
  DestroyImat(&unitFoothold);
}

//---------------------------------------------------------------------------------------------
// GetUnitAtomMap: make a map of all units to which any atom contributes.  This will anticipate
//                 further additions from the object type at hand, and also use a customized,
//                 packed imat that allocates different amounts of space to each row of the
//                 array.
//
// Arguments:
//   bwunits:        the list of bond work units currently allocated
//   nbwunits:       the number of bond work units currently allocated
//   objMapCounts:   counts of each atom's inclusion into objects of the type to be added
//                   to the work units
//   unitMapCounts:  counts of each atom's inclusion into bond work units (input initialized
//                   to zeros, returned with values)
//---------------------------------------------------------------------------------------------
static imat GetUnitAtomMap(bondwork *bwunits, int nbwunits, imat *objMapCounts,
                           imat *unitMapCounts)
{
  int i, j;

  // Get the number of atoms in the system from one of the imported matrices, to avoid
  // having to import the number (or the entire gpuContext) as an argument.
  int natom = objMapCounts->col;

  int *uMCounts, *itmp;
  uMCounts = unitMapCounts->data;
  for (i = 0; i < nbwunits; i++) {
    itmp = bwunits[i].atomList;
    for (j = 0; j < bwunits[i].natom; j++) {
      uMCounts[itmp[j]] += 1;
    }
  }
  itmp = objMapCounts->data;
  for (i = 0; i < natom; i++) {
    uMCounts[i] += itmp[i];
  }

  // Make the unit atom map
  imat unitMap;
  int ndata = ISum(uMCounts, natom);
  unitMap.data = (int*)malloc(ndata*sizeof(int));
  unitMap.map = (int**)malloc(natom*sizeof(int*));
  j = 0;
  for (i = 0; i < natom; i++) {
    unitMap.map[i] = &unitMap.data[j];
    j += uMCounts[i];
  }
  SetIVec(uMCounts, natom, 0);
  for (i = 0; i < nbwunits; i++) {
    itmp = bwunits[i].atomList;
    for (j = 0; j < bwunits[i].natom; j++) {
      unitMap.map[itmp[j]][uMCounts[itmp[j]]] = i;
      uMCounts[itmp[j]] += 1;
    }
  }

  return unitMap;
}

//---------------------------------------------------------------------------------------------
// CheckWorkUnitAssignments: this will check the accumulated work units for completeness.  All
//                           terms of a particular object type must have been assigned to a
//                           work unit or have been masked by MaskBlankObjects() in order for
//                           this to pass.
//
// Arguments:
//   gpu:       overarching structure holding all simulation information, critically the atom
//              IDs, parameters, and quantities of bonded interactions
//   objtype:   the type of object (i.e. bonds, dihedrals, NB 1-4s)
//   bwunits:   growing list of bond work units (also returned)
//   nbwunits:  number of bond work units (returned, passed by pointer)
//
// This is a debugging function.
//---------------------------------------------------------------------------------------------
static void CheckWorkUnitAssignments(gpuContext gpu, const char* objtype, imat *objAssigned,
                                     bondwork* bwunits, int nbwunits)
{
  int i, j, objcode=0;

  int nobj = GetSimulationObjectCount(gpu, -1, objtype);
  if (strcmp(objtype, "BOND") == 0) {
    objcode = BOND_CODE;
  }
  else if (strcmp(objtype, "ANGLE") == 0) {
    objcode = ANGL_CODE;
  }
  else if (strcmp(objtype, "DIHEDRAL") == 0) {
    objcode = DIHE_CODE;
  }
  else if (strcmp(objtype, "CMAP") == 0) {
    objcode = CMAP_CODE;
  }
  else if (strcmp(objtype, "QQXC") == 0) {
    objcode = QQXC_CODE;
  }
  else if (strcmp(objtype, "NB14") == 0) {
    objcode = NB14_CODE;
  }
  else if (strcmp(objtype, "NMR2") == 0) {
    objcode = NMR2_CODE;
  }
  else if (strcmp(objtype, "NMR3") == 0) {
    objcode = NMR3_CODE;
  }
  else if (strcmp(objtype, "NMR4") == 0) {
    objcode = NMR4_CODE;
  }
  else if (strcmp(objtype, "UREY") == 0) {
    objcode = UREY_CODE;
  }
  else if (strcmp(objtype, "CIMP") == 0) {
    objcode = CIMP_CODE;
  }
  else if (strcmp(objtype, "CNST") == 0) {
    objcode = CNST_CODE;
  }
  int* objcovered;
  objcovered = (int*)calloc(nobj, sizeof(int));
  imat objMasked;
  objMasked = CreateImat(1, nobj);
  MaskBlankObjects(gpu, objcode, &objMasked);
  int *listPtr;
  for (i = 0; i < nbwunits; i++) {
    if (objcode == BOND_CODE) {
      listPtr = bwunits[i].bondList;
    }
    else if (objcode == ANGL_CODE) {
      listPtr = bwunits[i].bondAngleList;
    }
    else if (objcode == DIHE_CODE) {
      listPtr = bwunits[i].dihedralList;
    }
    else if (objcode == CMAP_CODE) {
      listPtr = bwunits[i].cmapList;
    }
    else if (objcode == QQXC_CODE) {
      listPtr = bwunits[i].qqxcList;
    }
    else if (objcode == NB14_CODE) {
      listPtr = bwunits[i].nb14List;
    }
    else if (objcode == NMR2_CODE) {
      listPtr = bwunits[i].nmr2List;
    }
    else if (objcode == NMR3_CODE) {
      listPtr = bwunits[i].nmr3List;
    }
    else if (objcode == NMR4_CODE) {
      listPtr = bwunits[i].nmr4List;
    }
    else if (objcode == UREY_CODE) {
      listPtr = bwunits[i].ureyList;
    }
    else if (objcode == CIMP_CODE) {
      listPtr = bwunits[i].cimpList;
    }
    else if (objcode == CNST_CODE) {
      listPtr = bwunits[i].cnstList;
    }
    int jlim = GetWorkUnitObjectCount(&bwunits[i], objcode);
    for (j = 0; j < jlim; j++) {
      objcovered[listPtr[j]] += 1;
    }
  }
  for (i = 0; i < nobj; i++) {
    if (objAssigned->data[i] != 1) {
      printf("CheckWorkUnitAssignments :: Error.  No assignment for %s %4d.\n", objtype, i);
    }
    if (objcovered[i] != 1 && objMasked.data[i] == 0) {
      printf("CheckWorkUnitAssignments :: Error.  Coverage for %s %4d is %4d:", objtype, i,
             objcovered[i]);
      for (j = 0; j < nbwunits; j++) {
        if (objcode == BOND_CODE) {
          listPtr = bwunits[j].bondList;
        }
        else if (objcode == ANGL_CODE) {
          listPtr = bwunits[j].bondAngleList;
        }
        else if (objcode == DIHE_CODE) {
          listPtr = bwunits[j].dihedralList;
        }
        else if (objcode == CMAP_CODE) {
          listPtr = bwunits[j].cmapList;
        }
        else if (objcode == QQXC_CODE) {
          listPtr = bwunits[j].qqxcList;
        }
        else if (objcode == NB14_CODE) {
          listPtr = bwunits[j].nb14List;
        }
        else if (objcode == NMR2_CODE) {
          listPtr = bwunits[j].nmr2List;
        }
        else if (objcode == NMR3_CODE) {
          listPtr = bwunits[j].nmr3List;
        }
        else if (objcode == NMR4_CODE) {
          listPtr = bwunits[j].nmr4List;
        }
        else if (objcode == UREY_CODE) {
          listPtr = bwunits[j].ureyList;
        }
        else if (objcode == CIMP_CODE) {
          listPtr = bwunits[j].cimpList;
        }
        else if (objcode == CNST_CODE) {
          listPtr = bwunits[j].cnstList;
        }
        int k;
        int klim = GetWorkUnitObjectCount(&bwunits[j], objcode);
        for (k = 0; k < klim; k++) {
          if (listPtr[k] == i) {
            printf("%4d ", j);
          }
        }
      }
      printf("\n");
    }
  }
  DestroyImat(&objMasked);
}

//---------------------------------------------------------------------------------------------
// AccumulateWorkUnits: function for encapsulating what would otherwise be repetitious code.
//                      In order to commit all objects such as bonds, dihedrals, or nonbonded
//                      1-4 exclusions to the work units, we allocate some tables to track
//                      assignments and atom footprints, then attempt to commit the objects to
//                      any existing work units.  In a second pass, we make new work units as
//                      needed to hold any remaining objects.  Finally, we clean up memory that
//                      is no longer needed.
//
// Arguments:
//   gpu:         overarching structure holding all simulation information, critically the atom
//                IDs, parameters, and quantities of bonded interactions
//   objtype:     the type of object (i.e. bonds, dihedrals, NB 1-4s)
//   atomEndPts:  table indicating whether each atom is a terminal atom, connected to the rest
//                of some molecule by only one bond
//   bwunits:     growing list of bond work units (also returned)
//   nbwunits:    number of bond work units (returned, passed by pointer)
//   nbpivot:     start testing existing work units to accommodate each object at this point
//                (no earlier in the list)
//---------------------------------------------------------------------------------------------
static bondwork* AccumulateWorkUnits(gpuContext gpu, const char* objtype, imat *atomEndPts,
                                     bondwork* bwunits, int *nbwunits, int nbpivot)
{
  int nobj = GetSimulationObjectCount(gpu, -1, objtype);
  if (nobj > 0) {
    imat objAssigned, objMapCounts, objMap, objTerminates, unitMap, unitMapCounts;

    // Create tables
    objAssigned = CreateImat(1, nobj);
    GetObjectMaps(gpu, objtype, atomEndPts, &objMapCounts, &objMap, &objTerminates);

    // First pass to assign objects to existing work units
    unitMapCounts = CreateImat(1, gpu->sim.atoms);
    unitMap = GetUnitAtomMap(&bwunits[nbpivot], *nbwunits - nbpivot, &objMapCounts,
                             &unitMapCounts);
    ExtendWorkUnits(gpu, objtype, atomEndPts, &objAssigned, &objMapCounts, &objMap,
                    &objTerminates, &unitMap, &unitMapCounts, &bwunits[nbpivot],
                    *nbwunits - nbpivot);
    DestroyImat(&unitMap);
    DestroyImat(&unitMapCounts);

    // A second pass to create new work units as necessary
    unitMap = CreateImat(1, gpu->sim.atoms);
    SetIVec(unitMap.data, gpu->sim.atoms, -1);
    bwunits = BaseNGroups(gpu, objtype, atomEndPts, &objAssigned, &objMapCounts, &objMap,
                          &objTerminates, &unitMap, bwunits, nbwunits);
    DestroyImat(&unitMap);

    // Clear out units that have nothing in them
    int i, j;
    int* unitSurvives;
    unitSurvives = (int*)calloc(*nbwunits, sizeof(int));
    for (i = 0; i < *nbwunits; i++) {
      if (bwunits[i].nbond + bwunits[i].nangl + bwunits[i].ndihe + bwunits[i].ncmap +
          bwunits[i].nqqxc + bwunits[i].nnb14 + bwunits[i].nnmr2 + bwunits[i].nnmr3 +
          bwunits[i].nnmr4 + bwunits[i].nurey + bwunits[i].ncimp + bwunits[i].ncnst == 0) {
        DestroyBondWorkUnit(&bwunits[i]);
      }
      else {
        unitSurvives[i] = 1;
      }
    }
    j = 0;
    int nblank = 0;
    for (i = 0; i < *nbwunits; i++) {
      if (unitSurvives[i] == 1) {
        if (j < i) {
          bwunits[j] = bwunits[i];
        }
        j++;
      }
      else {
        nblank++;
      }
    }
    *nbwunits -= nblank;

    // Free allocated memory
    DestroyImat(&objAssigned);
    DestroyImat(&objMap);
    DestroyImat(&objMapCounts);
    DestroyImat(&objTerminates);
    free(unitSurvives);
  }

  return bwunits;
}

//---------------------------------------------------------------------------------------------
// SetWorkUnitObjectPointers: macro for setting pointers to a particular work unit's object
//                            list and imported atom map arrays.
//
// Arguments:
//   bwunits:  the array of all bond work units
///  i:        index into bwunits to focus on
//   objcode:  codified object identifier
//   listPtr:  pointer into the object ID list held by the work unit
//   mapPtr:   pointer into the list of imported atoms drawn upon by each object in the work
//             unit (as enumerated by contents of listPtr)
//   tiPtr:    pointer into the TI status array (not always needed, but present, so set anyway)
//---------------------------------------------------------------------------------------------
#ifndef SetWorkUnitObjectPointers
#define SetWorkUnitObjectPointers(bwunits, i, objcode, listPtr, mapPtr, tiPtr) \
{ \
  if (objcode == BOND_CODE) { \
    listPtr = bwunits[i].bondList; \
    mapPtr = bwunits[i].bondMap; \
    tiPtr = bwunits[i].bondStatus; \
  } \
  else if (objcode == ANGL_CODE) { \
    listPtr = bwunits[i].bondAngleList; \
    mapPtr = bwunits[i].bondAngleMap; \
    tiPtr = bwunits[i].anglStatus; \
  } \
  else if (objcode == DIHE_CODE) { \
    listPtr = bwunits[i].dihedralList; \
    mapPtr = bwunits[i].dihedralMap; \
    tiPtr = bwunits[i].diheStatus; \
  } \
  else if (objcode == CMAP_CODE) { \
    listPtr = bwunits[i].cmapList; \
    mapPtr = bwunits[i].cmapMap; \
    tiPtr = bwunits[i].cmapStatus; \
  } \
  else if (objcode == QQXC_CODE) { \
    listPtr = bwunits[i].qqxcList; \
    mapPtr = bwunits[i].qqxcMap; \
    tiPtr = bwunits[i].qqxcStatus; \
  } \
  else if (objcode == NB14_CODE) { \
    listPtr = bwunits[i].nb14List; \
    mapPtr = bwunits[i].nb14Map; \
    tiPtr = bwunits[i].nb14Status; \
  } \
  else if (objcode == NMR2_CODE) { \
    listPtr = bwunits[i].nmr2List; \
    mapPtr = bwunits[i].nmr2Map; \
    tiPtr = bwunits[i].nmr2Status; \
  } \
  else if (objcode == NMR3_CODE) { \
    listPtr = bwunits[i].nmr3List; \
    mapPtr = bwunits[i].nmr3Map; \
    tiPtr = bwunits[i].nmr3Status; \
  } \
  else if (objcode == NMR4_CODE) { \
    listPtr = bwunits[i].nmr4List; \
    mapPtr = bwunits[i].nmr4Map; \
    tiPtr = bwunits[i].nmr4Status; \
  } \
  else if (objcode == UREY_CODE) { \
    listPtr = bwunits[i].ureyList; \
    mapPtr = bwunits[i].ureyMap; \
    tiPtr = bwunits[i].ureyStatus; \
  } \
  else if (objcode == CIMP_CODE) { \
    listPtr = bwunits[i].cimpList; \
    mapPtr = bwunits[i].cimpMap; \
    tiPtr = bwunits[i].cimpStatus; \
  } \
  else if (objcode == CNST_CODE) { \
    listPtr = bwunits[i].cnstList; \
    mapPtr = bwunits[i].cnstMap; \
    tiPtr = bwunits[i].cnstStatus; \
  } \
}
#endif

//---------------------------------------------------------------------------------------------
// RemapObjectsToImportedAtoms: function for creating
//
// Arguments:
//   gpu:       overarching structure holding all simulation information, critically the atom
//              IDs, parameters, and quantities of bonded interactions
//   bwunits:   array of all bond work units
//   nbwunits:  length of bwunits
//   objtype:   type of object to remap
//---------------------------------------------------------------------------------------------
static void RemapObjectsToImportedAtoms(gpuContext gpu, bondwork* bwunits, int nbwunits,
                                        const char* objtype)
{
  int i, j, objcode=0;

  // Return immediately if there are no objects to map--don't
  // get caught mucking around with unassigned pointers.
  int nobj = GetSimulationObjectCount(gpu, 0, objtype);
  if (nobj == 0 || nbwunits == 0) {
    return;
  }

  // Set up pointers, as in other functions
  int *p1atm = NULL;
  int2 *p2atm = NULL;
  int4 *p4atm = NULL;
  codify_object_type(objtype, gpu, p1atm, p2atm, p4atm, objcode);

  // Allocate a scratch array for recording the positions of each atom in the work unit list
  // accoring to their ID numbers in the master topology.  The vast majority of atoms will
  // not be in the work unit, and the list will not be cleaned after each work unit, but
  // because every object in the unit is guaranteed to make use of atoms imported by the unit,
  // the array can be safely referenced to find the locations of each atom in the unit's list.
  int* unitAtomPos;
  unitAtomPos = (int*)malloc(gpu->sim.atoms*sizeof(int));

  // Remap the objects, unit by unit
  for (i = 0; i < nbwunits; i++) {

    // Remap all atoms of the work unit
    int *listPtr, *mapPtr, *tiPtr;
    listPtr = bwunits[i].atomList;
    for (j = 0; j < bwunits[i].natom; j++) {
      unitAtomPos[listPtr[j]] = j;
    }

    // Set up additional pointers into the work unit's data arrays
    SetWorkUnitObjectPointers(bwunits, i, objcode, listPtr, mapPtr, tiPtr);
    int nobj = GetWorkUnitObjectCount(&bwunits[i], objcode);
    for (j = 0; j < nobj; j++) {
      int objid = listPtr[j];
      if (objcode == CNST_CODE) {
        mapPtr[j] = unitAtomPos[p1atm[objid]];
      }
      else if (objcode == BOND_CODE || objcode == QQXC_CODE || objcode == NB14_CODE ||
               objcode == NMR2_CODE || objcode == UREY_CODE) {
        mapPtr[2*j    ] = unitAtomPos[p2atm[objid].x];
        mapPtr[2*j + 1] = unitAtomPos[p2atm[objid].y];
      }
      else if (objcode == ANGL_CODE || objcode == NMR3_CODE) {
        mapPtr[3*j    ] = unitAtomPos[p2atm[objid].x];
        mapPtr[3*j + 1] = unitAtomPos[p2atm[objid].y];
        mapPtr[3*j + 2] = unitAtomPos[p1atm[objid]];
      }
      else if (objcode == DIHE_CODE || objcode == NMR4_CODE || objcode == CIMP_CODE) {
        mapPtr[4*j    ] = unitAtomPos[p4atm[objid].x];
        mapPtr[4*j + 1] = unitAtomPos[p4atm[objid].y];
        mapPtr[4*j + 2] = unitAtomPos[p4atm[objid].z];
        mapPtr[4*j + 3] = unitAtomPos[p4atm[objid].w];
      }
      else if (objcode == CMAP_CODE) {
        mapPtr[5*j    ] = unitAtomPos[p4atm[objid].x];
        mapPtr[5*j + 1] = unitAtomPos[p4atm[objid].y];
        mapPtr[5*j + 2] = unitAtomPos[p4atm[objid].z];
        mapPtr[5*j + 3] = unitAtomPos[p4atm[objid].w];
        mapPtr[5*j + 4] = unitAtomPos[p1atm[objid]];
      }
    }
  }

  // Free allocated memory
  free(unitAtomPos);
}

//---------------------------------------------------------------------------------------------
// FindBestWarpPlacement: find the best warp in which to place an object.  This function can
//                        consider "flipping" the order of atoms that define an object so as
//                        to make it fit better.
//
// Arguments:
//   objidx:       index of the object in the work unit
//   nwarps:       the number of warps that this work unit will devote to objects of this type
//   flipOK:       flag indicating whether it is all right to flip the order of atoms in this
//                 object (not all objects are flippable)
//   occupancy:    the occupancy of the work unit's imported atoms, accumulated from placing
//                 other objects of this type
//---------------------------------------------------------------------------------------------
static int2 FindBestWarpPlacement(bondwork *bw, int objidx, int objcode, int nwarp,
                                  bool flipOK, imat *occupancy, imat *warpcounts, int *mapPtr)
{
  int i, ncf;
  int2 result;

  int bestcf = -1;
  for (i = 0; i < nwarp; i++) {
    if (warpcounts->data[i] >= GRID) {
      continue;
    }
    if (objcode == CNST_CODE) {
      ncf = occupancy->map[i][mapPtr[objidx]];
    }
    else if (objcode == BOND_CODE || objcode == QQXC_CODE || objcode == NB14_CODE ||
             objcode == NMR2_CODE || objcode == UREY_CODE) {
      int i1 = occupancy->map[i          ][mapPtr[2*objidx    ]];
      int i2 = occupancy->map[i +   nwarp][mapPtr[2*objidx + 1]];
      ncf = i1*i1 + i2*i2;
    }
    else if (objcode == ANGL_CODE || objcode == NMR3_CODE) {
      int i1 = occupancy->map[i          ][mapPtr[3*objidx    ]];
      int i2 = occupancy->map[i +   nwarp][mapPtr[3*objidx + 1]];
      int i3 = occupancy->map[i + 2*nwarp][mapPtr[3*objidx + 2]];
      ncf = i1*i1 + i2*i2 + i3*i3;
    }
    else if (objcode == DIHE_CODE || objcode == NMR4_CODE) {
      int i1 = occupancy->map[i          ][mapPtr[4*objidx    ]];
      int i2 = occupancy->map[i +   nwarp][mapPtr[4*objidx + 1]];
      int i3 = occupancy->map[i + 2*nwarp][mapPtr[4*objidx + 2]];
      int i4 = occupancy->map[i + 3*nwarp][mapPtr[4*objidx + 3]];
      ncf = i1*i1 + i2*i2 + i3*i3 + i4*i4;
    }
    else if (objcode == CMAP_CODE) {
      int i1 = occupancy->map[i          ][mapPtr[5*objidx    ]];
      int i2 = occupancy->map[i +   nwarp][mapPtr[5*objidx + 1]];
      int i3 = occupancy->map[i + 2*nwarp][mapPtr[5*objidx + 2]];
      int i4 = occupancy->map[i + 3*nwarp][mapPtr[5*objidx + 3]];
      int i5 = occupancy->map[i + 4*nwarp][mapPtr[5*objidx + 4]];
      ncf = i1*i1 + i2*i2 + i3*i3 + i4*i4 + i5*i5;
    }
    if (bestcf < 0 || ncf < bestcf) {
      bestcf = ncf;
      result.x = i;
      result.y = 0;
    }
    if (flipOK == false) {
      continue;
    }

    // Check the reverse order of object atoms
    if (objcode == BOND_CODE || objcode == QQXC_CODE || objcode == NB14_CODE ||
        objcode == NMR2_CODE || objcode == UREY_CODE) {
      int i1 = occupancy->map[i          ][mapPtr[2*objidx + 1]];
      int i2 = occupancy->map[i +   nwarp][mapPtr[2*objidx    ]];
      ncf = i1*i1 + i2*i2;
    }
    else if (objcode == ANGL_CODE || objcode == NMR3_CODE) {
      int i1 = occupancy->map[i          ][mapPtr[3*objidx + 2]];
      int i2 = occupancy->map[i +   nwarp][mapPtr[3*objidx + 1]];
      int i3 = occupancy->map[i + 2*nwarp][mapPtr[3*objidx    ]];
      ncf = i1*i1 + i2*i2 + i3*i3;
    }
    else if (objcode == DIHE_CODE || objcode == NMR4_CODE) {
      int i1 = occupancy->map[i          ][mapPtr[4*objidx + 3]];
      int i2 = occupancy->map[i +   nwarp][mapPtr[4*objidx + 2]];
      int i3 = occupancy->map[i + 2*nwarp][mapPtr[4*objidx + 1]];
      int i4 = occupancy->map[i + 3*nwarp][mapPtr[4*objidx    ]];
      ncf = i1*i1 + i2*i2 + i3*i3 + i4*i4;
    }
    if (bestcf < 0 || ncf < bestcf) {
      bestcf = ncf;
      result.x = i;
      result.y = 1;
    }
  }

  return result;
}

//---------------------------------------------------------------------------------------------
// PlotAtomOccupancy: plot the occupancy of an object's atoms against the bond work unit's
//                    list of imported atoms.
//
// Arguments:
//
//---------------------------------------------------------------------------------------------
static void PlotAtomOccupancy(int *mapPtr, int objcode, int objidx, bool flip, int warpidx,
                              int nwarp, imat *occupancy)
{
  if (objcode == CNST_CODE) {
    occupancy->map[warpidx][mapPtr[objidx]] += 1;
  }
  else if (objcode == BOND_CODE || objcode == QQXC_CODE || objcode == NB14_CODE ||
           objcode == NMR2_CODE || objcode == UREY_CODE) {
    if (flip) {
      occupancy->map[warpidx          ][mapPtr[2*objidx + 1]] += 1;
      occupancy->map[warpidx +   nwarp][mapPtr[2*objidx    ]] += 1;
    }
    else {
      occupancy->map[warpidx          ][mapPtr[2*objidx    ]] += 1;
      occupancy->map[warpidx +   nwarp][mapPtr[2*objidx + 1]] += 1;
    }
  }
  else if (objcode == ANGL_CODE || objcode == NMR3_CODE) {
    if (flip) {
      occupancy->map[warpidx          ][mapPtr[3*objidx + 2]] += 1;
      occupancy->map[warpidx +   nwarp][mapPtr[3*objidx + 1]] += 1;
      occupancy->map[warpidx + 2*nwarp][mapPtr[3*objidx    ]] += 1;
    }
    else {
      occupancy->map[warpidx          ][mapPtr[3*objidx    ]] += 1;
      occupancy->map[warpidx +   nwarp][mapPtr[3*objidx + 1]] += 1;
      occupancy->map[warpidx + 2*nwarp][mapPtr[3*objidx + 2]] += 1;
    }
  }
  else if (objcode == DIHE_CODE || objcode == NMR4_CODE) {
    if (flip) {
      occupancy->map[warpidx          ][mapPtr[4*objidx + 3]] += 1;
      occupancy->map[warpidx +   nwarp][mapPtr[4*objidx + 2]] += 1;
      occupancy->map[warpidx + 2*nwarp][mapPtr[4*objidx + 1]] += 1;
      occupancy->map[warpidx + 3*nwarp][mapPtr[4*objidx    ]] += 1;
    }
    else {
      occupancy->map[warpidx          ][mapPtr[4*objidx    ]] += 1;
      occupancy->map[warpidx +   nwarp][mapPtr[4*objidx + 1]] += 1;
      occupancy->map[warpidx + 2*nwarp][mapPtr[4*objidx + 2]] += 1;
      occupancy->map[warpidx + 3*nwarp][mapPtr[4*objidx + 3]] += 1;
    }
  }
  else if (objcode == CMAP_CODE) {
    if (flip) {
      occupancy->map[warpidx          ][mapPtr[5*objidx + 4]] += 1;
      occupancy->map[warpidx +   nwarp][mapPtr[5*objidx + 3]] += 1;
      occupancy->map[warpidx + 2*nwarp][mapPtr[5*objidx + 2]] += 1;
      occupancy->map[warpidx + 3*nwarp][mapPtr[5*objidx + 1]] += 1;
      occupancy->map[warpidx + 4*nwarp][mapPtr[5*objidx    ]] += 1;
    }
    else {
      occupancy->map[warpidx          ][mapPtr[5*objidx    ]] += 1;
      occupancy->map[warpidx +   nwarp][mapPtr[5*objidx + 1]] += 1;
      occupancy->map[warpidx + 2*nwarp][mapPtr[5*objidx + 2]] += 1;
      occupancy->map[warpidx + 3*nwarp][mapPtr[5*objidx + 3]] += 1;
      occupancy->map[warpidx + 4*nwarp][mapPtr[5*objidx + 4]] += 1;
    }
  }
}

//---------------------------------------------------------------------------------------------
// PlaceObjectInWarp: place an object in the occupancy matrix to reflect its inclusion in the
//                    set of objects handled by one warp.
//
// Arguments:
//
//---------------------------------------------------------------------------------------------
static void PlaceObjectInWarp(int *mapPtr, int objcode, int objidx, int warpidx, bool flip,
                              int nwarp, imat *occupancy, imat *neworder, imat *warpcounts)
{
  // Contribute to occupancy
  PlotAtomOccupancy(mapPtr, objcode, objidx, flip, warpidx, nwarp, occupancy);

  // Assign in the ordering
  neworder->map[0][warpidx*GRID + warpcounts->data[warpidx]] = objidx;
  neworder->map[1][warpidx*GRID + warpcounts->data[warpidx]] = (flip);
  warpcounts->data[warpidx] += 1;
}

//---------------------------------------------------------------------------------------------
// CountWarpAtomConflicts:
//
// This is a debugging and performance checking function.
//---------------------------------------------------------------------------------------------
static int CountWarpAtomConflicts(int *mapPtr, imat *order, const char* objtype, int objcode,
                                  int nwarp)
{
  int i, j, k;

  // Allocate an occupancy map (do this locally so as not to blow away anything else)
  imat occupancy;
  int maxwarps = 16 * BOND_WORK_UNIT_WARPS_PER_BLOCK;
  occupancy = CreateImat(maxwarps, BOND_WORK_UNIT_THREADS_PER_BLOCK);

  // Loop over all warps and plot the occupancy
  for (i = 0; i < nwarp; i++) {
    for (j = 0; j < GRID; j++) {
      if (order->map[0][i*GRID + j] >= 0) {
        PlotAtomOccupancy(mapPtr, objcode, order->map[0][i*GRID + j],
                          (order->map[1][i*GRID + j] == 1), i, nwarp, &occupancy);
      }
    }
  }

  // Get the number of reads for this type of object and total the conflicts
  int nsets;
  if (objcode == CNST_CODE) {
    nsets = 1;
  }
  else if (objcode == BOND_CODE || objcode == QQXC_CODE || objcode == NB14_CODE ||
           objcode == NMR2_CODE || objcode == UREY_CODE) {
    nsets = 2;
  }
  else if (objcode == ANGL_CODE || objcode == NMR3_CODE) {
    nsets = 3;
  }
  else if (objcode == DIHE_CODE || objcode == NMR4_CODE) {
    nsets = 4;
  }
  else if (objcode == CMAP_CODE) {
    nsets = 5;
  }
  printf("Conflicts for %s :: ", objtype);
  for (i = 0; i < nwarp; i++) {
    int ncf = 0;
    int ncf2 = 0;
    for (j = 0; j < nsets; j++) {
      int *itmp = occupancy.map[i + j*nwarp];
      for (k = 0; k < BOND_WORK_UNIT_THREADS_PER_BLOCK; k++) {
        int itkm1 = itmp[k] - (itmp[k] > 0);
        ncf += itkm1;
        ncf2 += itkm1*itkm1;
      }
    }
    printf("%4d [ %4d ]", ncf, ncf2);
  }
  printf("\n");

  // Free allocated memory
  DestroyImat(&occupancy);
}

//---------------------------------------------------------------------------------------------
// ReorderWorkUnitObjects: re-order objects of the stated type in a work unit to avoid
//                         processing objects requiring the same atoms in the same warps.
//
// Arguments:
//   gpu:         overarching structure holding all simulation information, critically the atom
//                IDs, parameters, and quantities of bonded interactions
//   bwunits:     array of all bond work units
//   nbwunits:    length of bwunits
//   objtype:     the type of object (i.e. bonds, dihedrals, NB 1-4s)
//   flipOK:      flag to indicate that it is all right flip the order of atoms in an object
//---------------------------------------------------------------------------------------------
static void ReorderWorkUnitObjects(gpuContext gpu, bondwork* bwunits, int nbwunits,
                                   const char* objtype, bool flipOK)
{
  int i, j, k, objcode=0;

  // Return immediately if there are no objects to map--don't
  // get caught mucking around with unassigned pointers.
  int nobj = GetSimulationObjectCount(gpu, -1, objtype);
  if (nobj == 0 || nbwunits == 0) {
    return;
  }

  // Set pointers and codify the object type
  int *p1atm = NULL;
  int2 *p2atm = NULL;
  int4 *p4atm = NULL;
  codify_object_type(objtype, gpu, p1atm, p2atm, p4atm, objcode);
  int objsize = GetWorkUnitObjectSize(objcode);

  // Allocate a table for re-ordering the objects, then a
  // table for the occupancy across all imported atoms
  imat neworder, occupancy, warpcounts, newmap;
  int maxwarps = 16 * BOND_WORK_UNIT_WARPS_PER_BLOCK;
  neworder = CreateImat(2, maxwarps * GRID);
  occupancy = CreateImat(maxwarps, BOND_WORK_UNIT_THREADS_PER_BLOCK);
  warpcounts = CreateImat(1, maxwarps);
  newmap = CreateImat(1, 2 * maxwarps * GRID);

  // Each unit gets handled separately
  for (i = 0; i < nbwunits; i++) {

    // Initialize the re-ordering and occupancy tables
    SetIVec(neworder.data, 8*BOND_WORK_UNIT_THREADS_PER_BLOCK, -1);
    SetIVec(occupancy.data, maxwarps * BOND_WORK_UNIT_THREADS_PER_BLOCK, 0);
    SetIVec(warpcounts.data, maxwarps, 0);

    // Set the object list and atom map pointers
    int *listPtr, *mapPtr, *tiPtr;
    SetWorkUnitObjectPointers(bwunits, i, objcode, listPtr, mapPtr, tiPtr);
    int nobj = GetWorkUnitObjectCount(&bwunits[i], objcode);
    int nwarp = (nobj + GRID_BITS_MASK) / GRID;

    // Set all unfilled slots in the work unit's object list to -1 to indicate they are blanks
    int maxobj = BOND_WORK_UNIT_THREADS_PER_BLOCK;
    if (objcode == QQXC_CODE || objcode == NB14_CODE) {
      maxobj *= 8;
    }
    else if (objcode == DIHE_CODE) {
      maxobj *= 2;
    }
    for (j = nobj; j < maxobj; j++) {
      listPtr[j] = -1;
    }

    // Loop over all objects, assign them to warps.
    for (j = 0; j < nobj; j++) {
      int2 plc = FindBestWarpPlacement(&bwunits[i], j, objcode, nwarp, flipOK, &occupancy,
                                       &warpcounts, mapPtr);
      PlaceObjectInWarp(mapPtr, objcode, j, plc.x, (plc.y == 1), nwarp, &occupancy, &neworder,
                        &warpcounts);
    }

    // Make a new map of the imported atoms for the new order of objects
    for (j = 0; j < nwarp*GRID; j++) {

      // Skip blanks
      if (neworder.map[0][j] < 0) {
        continue;
      }
      if (neworder.map[1][j] == 0) {
        for (k = 0; k < objsize; k++) {
          newmap.data[j*objsize + k] = mapPtr[(objsize * neworder.map[0][j]) + k];
        }
      }
      else {
        for (k = 0; k < objsize; k++) {
          newmap.data[j*objsize + k] = mapPtr[(objsize * neworder.map[0][j]) +
                                              objsize - k - 1];
        }
      }

      // Replace the new ordering with the object ID from the master topology,
      // preparing to replace the list in the bond work unit from neworder.
      if (neworder.map[0][j] >= 0) {
        neworder.map[0][j] = listPtr[neworder.map[0][j]];
      }
    }
    for (j = 0; j < nwarp*GRID*objsize; j++) {
      mapPtr[j] = newmap.data[j];
    }
    for (j = 0; j < nwarp*GRID; j++) {
      if (neworder.map[0][j] < 0) {
        listPtr[j] = -1;
      }
      else {
        listPtr[j] = neworder.map[0][j];
      }
    }

    // Set the work unit object count to reflect the warp padding
    if (objcode == BOND_CODE) {
      bwunits[i].nbond = nwarp * GRID;
    }
    else if (objcode == ANGL_CODE) {
      bwunits[i].nangl = nwarp * GRID;
    }
    else if (objcode == DIHE_CODE) {
      bwunits[i].ndihe = nwarp * GRID;
    }
    else if (objcode == CMAP_CODE) {
      bwunits[i].ncmap = nwarp * GRID;
    }
    else if (objcode == QQXC_CODE) {
      bwunits[i].nqqxc = nwarp * GRID;
    }
    else if (objcode == NB14_CODE) {
      bwunits[i].nnb14 = nwarp * GRID;
    }
    else if (objcode == NMR2_CODE) {
      bwunits[i].nnmr2 = nwarp * GRID;
    }
    else if (objcode == NMR3_CODE) {
      bwunits[i].nnmr3 = nwarp * GRID;
    }
    else if (objcode == NMR4_CODE) {
      bwunits[i].nnmr4 = nwarp * GRID;
    }
    else if (objcode == UREY_CODE) {
      bwunits[i].nurey = nwarp * GRID;
    }
    else if (objcode == CIMP_CODE) {
      bwunits[i].ncimp = nwarp * GRID;
    }
    else if (objcode == CNST_CODE) {
      bwunits[i].ncnst = nwarp * GRID;
    }
  }

  // Free allocated memory
  DestroyImat(&neworder);
  DestroyImat(&occupancy);
  DestroyImat(&warpcounts);
  DestroyImat(&newmap);
}

//---------------------------------------------------------------------------------------------
// GetAtomEndPts: make a table of all atoms that are connect to their parent molecules by only
//                one bond.  This needs to be done in a couple of places, hence encapsulation.
//
// Arguments:
//   gpu:   overarching data type of simulation parameters, here used for atom and bond counts
//---------------------------------------------------------------------------------------------
static imat GetAtomEndPts(gpuContext gpu)
{
  int i;
  imat atomEndPts;

  // Unpack the gpu data structure a bit
  int natom = gpu->sim.atoms;
  int nbond = gpu->sim.bonds;

  atomEndPts = CreateImat(1, natom);
  int2 *pBondID;
  pBondID = gpu->pbBondID->_pSysData;
  for (i = 0; i < nbond; i++) {
    atomEndPts.data[pBondID[i].x] += 1;
    atomEndPts.data[pBondID[i].y] += 1;
  }
  for (i = 0; i < natom; i++) {
    atomEndPts.data[i] = (atomEndPts.data[i] == 1);
  }

  return atomEndPts;
}

//---------------------------------------------------------------------------------------------
// GetBwuObjectStatus: get the Thermodynamic Integration status of an object in a bond work
//                     unit.  This encapsulates the tests that used to be done at every step
//                     on the GPU.
//
// Arguments:
//   gpu:       overarching structure holding all simulation information, critically the atom
//              IDs, parameters, and quantities of bonded interactions
//   bwunits:   complete array of existing bond work units (this is a "capstone" function,
//              intended to be called at the end of the process)
//   nbwunits:  the number of work units allocated (returned)
//   objtype:   type of object
//---------------------------------------------------------------------------------------------
static void GetBwuObjectStatus(gpuContext gpu, bondwork* bwunits, int nbwunits,
                               const char* objtype)
{
  int i, j, objcode=0;

  // Return immediately if there are no objects to map--don't
  // get caught mucking around with unassigned pointers.
  int nobj = GetSimulationObjectCount(gpu, -1, objtype);
  if (nobj == 0 || nbwunits == 0) {
    return;
  }

  // Set up pointers, as in other functions
  int *p1atm = NULL;
  int2 *p2atm = NULL;
  int4 *p4atm = NULL;
  codify_object_type(objtype, gpu, p1atm, p2atm, p4atm, objcode);
  int objsize = GetWorkUnitObjectSize(objcode);

  for (i = 0; i < nbwunits; i++) {
    int nobj = GetWorkUnitObjectCount(&bwunits[i], objcode);
    int *listPtr, *mapPtr, *tiPtr;
    SetWorkUnitObjectPointers(bwunits, i, objcode, listPtr, mapPtr, tiPtr);
    for (j = 0; j < nobj; j++) {

      // Search all atoms in the object, seeking three results: whether the
      // object is classified as "CV" or "soft core," and the TI region to
      // which it belongs.  If one or more atoms are in the CV region, it
      // is "CV" so long as no atoms are in the SC (soft core) region.  If
      // any atoms are in the soft core region, it is classified as SC.
      int CVterm, SCterm, TIObjRegion;
      int objid = listPtr[j];
      if (objid == 0xffffffff) {
        continue;
      }
      if (objsize == 1) {
        int stti = gpu->pbTIRegion->_pSysData[p1atm[objid]];
        CVterm = ((stti >> 1) > 0);
        SCterm = stti & 0x1;
        TIObjRegion = (stti > 0) ? (stti >> 1) : 0;
      }
      else if (objsize == 2) {
        int stti = gpu->pbTIRegion->_pSysData[p2atm[objid].x];
        int sttj = gpu->pbTIRegion->_pSysData[p2atm[objid].y];
        CVterm = ((stti >> 1) > 0) || ((sttj >> 1) > 0);
        SCterm = (stti | sttj) & 0x1;
        TIObjRegion = ((stti | sttj) > 0) ? ((stti | sttj) >> 1) : 0;
      }
      else if (objsize == 3) {
        int stti = gpu->pbTIRegion->_pSysData[p2atm[objid].x];
        int sttj = gpu->pbTIRegion->_pSysData[p2atm[objid].y];
        int sttk = gpu->pbTIRegion->_pSysData[p1atm[objid]];
        CVterm = ((stti >> 1) > 0) || ((sttj >> 1) > 0) || ((sttk >> 1) > 0);
        SCterm = (stti | sttj | sttk) & 0x1;
        TIObjRegion = ((stti | sttj | sttk) > 0) ? ((stti | sttj | sttk) >> 1) : 0;
      }
      else if (objsize == 4) {
        int stti = gpu->pbTIRegion->_pSysData[p4atm[objid].x];
        int sttj = gpu->pbTIRegion->_pSysData[p4atm[objid].y];
        int sttk = gpu->pbTIRegion->_pSysData[p4atm[objid].z];
        int sttl = gpu->pbTIRegion->_pSysData[p4atm[objid].w];
        CVterm = ((stti >> 1) > 0) || ((sttj >> 1) > 0) || ((sttk >> 1) > 0) ||
                 ((sttl >> 1) > 0);
        SCterm = (stti | sttj | sttk | sttl) & 0x1;
        TIObjRegion = ((stti | sttj | sttk | sttl) > 0) ?
                      ((stti | sttj | sttk | sttl) >> 1) : 0;
      }
      else if (objsize == 5) {
        int stti = gpu->pbTIRegion->_pSysData[p4atm[objid].x];
        int sttj = gpu->pbTIRegion->_pSysData[p4atm[objid].y];
        int sttk = gpu->pbTIRegion->_pSysData[p4atm[objid].z];
        int sttl = gpu->pbTIRegion->_pSysData[p4atm[objid].w];
        int sttm = gpu->pbTIRegion->_pSysData[p1atm[objid]];
        CVterm = ((stti >> 1) > 0) || ((sttj >> 1) > 0) || ((sttk >> 1) > 0) ||
                 ((sttl >> 1) > 0) || ((sttm >> 1) > 0);
        SCterm = (stti | sttj | sttk | sttl | sttm) & 0x1;
        TIObjRegion = ((stti | sttj | sttk | sttl | sttm) > 0) ?
                      ((stti | sttj | sttk | sttl | sttm) >> 1) : 0;
      }
      if (SCterm == 1) {
        CVterm = 0;
      }

      // Pack all results into one integer
      tiPtr[j] = (TIObjRegion << 16) | (SCterm << 8) | CVterm;
    }
  }
}

//---------------------------------------------------------------------------------------------
// AssembleBondWorkUnits: create collections of atoms which can be imported as a group to cover
//                        bonds, electrostatic exclusions, bond angles, dihedrals, and CMAP
//                        terms.  Results are passed back up to the calling function in gpu.cpp
//                        to allow new GpuBuffers to get allocated before more functions in
//                        this library work on them.
//
// Arguments:
//   gpu:       overarching structure holding all simulation information, critically the atom
//              IDs, parameters, and quantities of bonded interactions
//   nbwunits:  the number of work units allocated (returned)
//---------------------------------------------------------------------------------------------
bondwork* AssembleBondWorkUnits(gpuContext gpu, int *nbwunits)
{
  // Start the bond work untis array as a tiny bit of memory--it will be expanded as needed
  bondwork* bwunits;
  *nbwunits = 0;
  bwunits = (bondwork*)malloc(sizeof(bondwork));

  imat atomEndPts;
  atomEndPts = GetAtomEndPts(gpu);

  // Map all types of objects and allocate assignment tables.  Once BaseNGroups has been
  // called, all objects of that type will have been placed in work units, so a lot of
  // memory can be free'd immediately.
  bwunits = AccumulateWorkUnits(gpu, "CMAP", &atomEndPts, bwunits, nbwunits, 0);
  bwunits = AccumulateWorkUnits(gpu, "DIHEDRAL", &atomEndPts, bwunits, nbwunits, 0);
  bwunits = AccumulateWorkUnits(gpu, "CIMP", &atomEndPts, bwunits, nbwunits, 0);
  bwunits = AccumulateWorkUnits(gpu, "ANGLE", &atomEndPts, bwunits, nbwunits, 0);
  bwunits = AccumulateWorkUnits(gpu, "BOND", &atomEndPts, bwunits, nbwunits, 0);
  bwunits = AccumulateWorkUnits(gpu, "UREY", &atomEndPts, bwunits, nbwunits, 0);
  bwunits = AccumulateWorkUnits(gpu, "NMR2", &atomEndPts, bwunits, nbwunits, 0);
  bwunits = AccumulateWorkUnits(gpu, "NMR3", &atomEndPts, bwunits, nbwunits, 0);
  bwunits = AccumulateWorkUnits(gpu, "NMR4", &atomEndPts, bwunits, nbwunits, 0);
  bwunits = AccumulateWorkUnits(gpu, "CNST", &atomEndPts, bwunits, nbwunits, 0);

  // In the case of NTP with a Berendsen barostat, bond work units for 1-4 non-bonded and
  // electrostatic non-bonded exclusions must have their forces mapped back to the non-bonded
  // force array, and so must be handled separately from objects that accumulate into the
  // bonded force array.
  if (gpu->sim.ntp > 0 && gpu->sim.barostat == 1) {
    int nbpivot = *nbwunits;
    bwunits = AccumulateWorkUnits(gpu, "NB14", &atomEndPts, bwunits, nbwunits, nbpivot);
    for (int i = nbpivot; i < *nbwunits; i++) {
      bwunits[i].frcAcc = NB_FORCE_ACCUMULATOR;
    }
  }
  else {
    bwunits = AccumulateWorkUnits(gpu, "NB14", &atomEndPts, bwunits, nbwunits, 0);
  }

  // Print an almanac of the bond work units
#ifdef GVERBOSE
  printf("Work units:\n");
  for (int i = 0; i < *nbwunits; i++) {
    printf("%4d :: %4d   %4d %4d %4d %4d   %4d %4d %4d   %4d %4d %4d   %4d %4d\n", i,
           bwunits[i].natom, bwunits[i].nbond, bwunits[i].nangl, bwunits[i].ndihe,
           bwunits[i].ncnst, bwunits[i].nnmr2, bwunits[i].nnmr3, bwunits[i].nnmr4,
           bwunits[i].nurey, bwunits[i].ncimp, bwunits[i].ncmap, bwunits[i].nqqxc,
           bwunits[i].nnb14);
  }
  printf("\n");
#endif

  // Map each work unit's objects to the array of atoms it will import
  RemapObjectsToImportedAtoms(gpu, bwunits, *nbwunits, "BOND");
  RemapObjectsToImportedAtoms(gpu, bwunits, *nbwunits, "ANGLE");
  RemapObjectsToImportedAtoms(gpu, bwunits, *nbwunits, "DIHEDRAL");
  RemapObjectsToImportedAtoms(gpu, bwunits, *nbwunits, "CMAP");
  RemapObjectsToImportedAtoms(gpu, bwunits, *nbwunits, "QQXC");
  RemapObjectsToImportedAtoms(gpu, bwunits, *nbwunits, "NB14");
  RemapObjectsToImportedAtoms(gpu, bwunits, *nbwunits, "NMR2");
  RemapObjectsToImportedAtoms(gpu, bwunits, *nbwunits, "NMR3");
  RemapObjectsToImportedAtoms(gpu, bwunits, *nbwunits, "NMR4");
  RemapObjectsToImportedAtoms(gpu, bwunits, *nbwunits, "UREY");
  RemapObjectsToImportedAtoms(gpu, bwunits, *nbwunits, "CIMP");
  RemapObjectsToImportedAtoms(gpu, bwunits, *nbwunits, "CNST");

  // Re-organize the objects in each work unit so as not to create bank conflicts
  ReorderWorkUnitObjects(gpu, bwunits, *nbwunits, "BOND", true);
  ReorderWorkUnitObjects(gpu, bwunits, *nbwunits, "ANGLE", true);
  ReorderWorkUnitObjects(gpu, bwunits, *nbwunits, "DIHEDRAL", true);
  ReorderWorkUnitObjects(gpu, bwunits, *nbwunits, "CMAP", false);
  ReorderWorkUnitObjects(gpu, bwunits, *nbwunits, "QQXC", true);
  ReorderWorkUnitObjects(gpu, bwunits, *nbwunits, "NB14", true);
  ReorderWorkUnitObjects(gpu, bwunits, *nbwunits, "NMR2", true);
  ReorderWorkUnitObjects(gpu, bwunits, *nbwunits, "NMR3", true);
  ReorderWorkUnitObjects(gpu, bwunits, *nbwunits, "NMR4", true);
  ReorderWorkUnitObjects(gpu, bwunits, *nbwunits, "UREY", true);
  ReorderWorkUnitObjects(gpu, bwunits, *nbwunits, "CIMP", false);
  ReorderWorkUnitObjects(gpu, bwunits, *nbwunits, "CNST", false);

  // For TI modes, encode the status of each object--does it count as a CV object,
  // SC object, or a standard object?
  if (gpu->sim.ti_mode > 0) {
    GetBwuObjectStatus(gpu, bwunits, *nbwunits, "BOND");
    GetBwuObjectStatus(gpu, bwunits, *nbwunits, "ANGLE");
    GetBwuObjectStatus(gpu, bwunits, *nbwunits, "DIHEDRAL");
    GetBwuObjectStatus(gpu, bwunits, *nbwunits, "CMAP");
    GetBwuObjectStatus(gpu, bwunits, *nbwunits, "QQXC");
    GetBwuObjectStatus(gpu, bwunits, *nbwunits, "NB14");
    GetBwuObjectStatus(gpu, bwunits, *nbwunits, "NMR2");
    GetBwuObjectStatus(gpu, bwunits, *nbwunits, "NMR3");
    GetBwuObjectStatus(gpu, bwunits, *nbwunits, "NMR4");
    GetBwuObjectStatus(gpu, bwunits, *nbwunits, "UREY");
    GetBwuObjectStatus(gpu, bwunits, *nbwunits, "CIMP");
    GetBwuObjectStatus(gpu, bwunits, *nbwunits, "CNST");
  }

  // Free allocated memory
  DestroyImat(&atomEndPts);

  return bwunits;
}

//---------------------------------------------------------------------------------------------
// CalcBondBlockSize: calculate the total amount of memory needed to store the bond work units'
//                    representations on the GPU (and allocate that amount on both sides of
//                    the fence so that the CPU can fill it up to port up to the device).  See
//                    the description for MakeBondedUnitDirections() below to get a better
//                    sense of what all the memory will hold.  This will return a tuple of
//                    ints encoding how much of various data types to allocate and how many
//                    warps will be needed to handle each type of bonded object.
//
// Arguments:
//   gpu:       overarching structure holding all simulation information, critically the atom
//              IDs, parameters, and quantities of bonded interactions
//   bwunits:   the bond work units
//   nbwunits:  the length of bwunits
//---------------------------------------------------------------------------------------------
bwalloc CalcBondBlockSize(gpuContext gpu, bondwork* bwunits, int nbwunits)
{
  int i;
  bwalloc bwdims;

  // The integer allocations serve instruction sets (12 sets of 32 ints per unit) and
  // Lennard-Jones parameter indices (4 sets of 32 ints per unit).
  bwdims.nUint = nbwunits * 3 * BOND_WORK_UNIT_THREADS_PER_BLOCK;
  bwdims.nbondwarps = 0;
  bwdims.nanglwarps = 0;
  bwdims.ndihewarps = 0;
  bwdims.ncmapwarps = 0;
  bwdims.nqqxcwarps = 0;
  bwdims.nnb14warps = 0;
  bwdims.nnmr2warps = 0;
  bwdims.nnmr3warps = 0;
  bwdims.nnmr4warps = 0;
  bwdims.nureywarps = 0;
  bwdims.ncimpwarps = 0;
  bwdims.ncnstwarps = 0;
  for (i = 0; i < nbwunits; i++) {
    bwdims.nbondwarps += (bwunits[i].nbond + GRID_BITS_MASK) >> GRID_BITS;
    bwdims.nanglwarps += (bwunits[i].nangl + GRID_BITS_MASK) >> GRID_BITS;
    bwdims.ndihewarps += (bwunits[i].ndihe + GRID_BITS_MASK) >> GRID_BITS;
    bwdims.ncmapwarps += (bwunits[i].ncmap + GRID_BITS_MASK) >> GRID_BITS;
    bwdims.nqqxcwarps += (bwunits[i].nqqxc + GRID_BITS_MASK) >> GRID_BITS;
    bwdims.nnb14warps += (bwunits[i].nnb14 + GRID_BITS_MASK) >> GRID_BITS;
    bwdims.nnmr2warps += (bwunits[i].nnmr2 + GRID_BITS_MASK) >> GRID_BITS;
    bwdims.nnmr3warps += (bwunits[i].nnmr3 + GRID_BITS_MASK) >> GRID_BITS;
    bwdims.nnmr4warps += (bwunits[i].nnmr4 + GRID_BITS_MASK) >> GRID_BITS;
    bwdims.nureywarps += (bwunits[i].nurey + GRID_BITS_MASK) >> GRID_BITS;
    bwdims.ncimpwarps += (bwunits[i].ncimp + GRID_BITS_MASK) >> GRID_BITS;
    bwdims.ncnstwarps += (bwunits[i].ncnst + GRID_BITS_MASK) >> GRID_BITS;
  }
  bwdims.nUint += (  bwdims.nbondwarps +   bwdims.nanglwarps +   bwdims.ndihewarps +
                   2*bwdims.ncmapwarps +   bwdims.nqqxcwarps +   bwdims.nnb14warps +
                   4*bwdims.nnmr2warps + 4*bwdims.nnmr3warps + 4*bwdims.nnmr4warps +
                     bwdims.nureywarps +   bwdims.ncimpwarps + 2*bwdims.ncnstwarps) * GRID;

  // Double-precision parameter data
  bwdims.nDbl2 = (  bwdims.nbondwarps +   bwdims.nanglwarps + 6*bwdims.nnmr2warps +
                  6*bwdims.nnmr3warps + 6*bwdims.nnmr4warps +   bwdims.nureywarps +
                    bwdims.ncimpwarps + 2*bwdims.ncnstwarps) * GRID;

  // PMEFloats for partial charges and dihedral parameters
  bwdims.nPFloat = (nbwunits * BOND_WORK_UNIT_THREADS_PER_BLOCK) +
                   (bwdims.ndihewarps + bwdims.nnb14warps)*GRID;
  bwdims.nPFloat2 = (2*bwdims.ndihewarps + bwdims.nnb14warps) * GRID;

  // Adjustments to accommodate Thermodynamic Integration
  if (gpu->sim.ti_mode > 0) {
    int nTIval = (bwdims.nbondwarps + bwdims.nanglwarps + bwdims.ndihewarps +
                  bwdims.ncmapwarps + bwdims.nqqxcwarps + bwdims.nnb14warps +
                  bwdims.nnmr2warps + bwdims.nnmr3warps + bwdims.nnmr4warps +
                  bwdims.nureywarps + bwdims.ncimpwarps + bwdims.ncnstwarps) * GRID;
    bwdims.nUint += nTIval;
  }

  return bwdims;
}

//---------------------------------------------------------------------------------------------
// MakeBondedUnitDirections: lay out set of directions needed for processing a unit of bonded
//                           interactions into the system-side array of the gpuContext.  Once
//                           all units have been laid out, the unit is ready for upload to the
//                           device.
//
// The directions take the form of a massive stream of integers, with some degree of order
// if broken down by 32-int sets.  The goal is to read a small number of atoms (20-40) that
// comprise all the coordinates needed by the 32 dihedral computations (it COULD be 768
// unique atoms, but this would only happen in particular circumstances, like a simulation of
// a fluid wherein each molecule has precisely one dihedral).  Then, the small number of atoms
// gets translated and copied ("expanded") into an array of 128 where each dihedral will work
// on its own four atoms.  Electrostatic exclusions, bonds, angles, and even CMAP terms then
// work off the expanded array of 128, computing forces and energies using indices pre-comuted
// on the CPU and loaded with the directions to avoid memory bank conflicts whenever possible.
//
// Each block will run 128 threads: four warps.  All threads will be used in reading
// directions, and then reading atoms from the sub-images into __shared__.  The warps will then
// specialize according to the directions they read.  Finally, all warps will cooperate in
// merging forces and writing results back to global memory.
//
//    Set of 32 ints   Description
//     SP*P    DPFP
//    --------------   ------------------------------------------------------------------------
//     1- 4    1- 4    Critical counts for the entire thread block (these will be read into
//                     __shared__):
//                     - Number of atom imports
//                     - Source of atom imports (pImage for double precision, pQQSubImageSP
//                         for single precision)
//                     - Number of warp instructions (each instruction is an unsigned int,
//                         packed with the type of interaction to perform and the starting
//                         point in the list of such interactions)
//                     - List of warp instructions (up to 120 available)
//     5- 8    5- 8    Atom ID number imports (these are 'absolute' IDs, referencing the order
//                     in the original topology.  They will be cross-referenced against the
//                     pSubImageAtomLkp array.  The IDs are not read into __shared__, but
//                     coordinates are.
//     9-12    9-12    Atom ID number imports indexing into the pImage (for double precision)
//                     or pSubImage (for single precision) arrays.  The single precision arrays
//                     are accessed if the work unit contains only electrostatic exclusions.
//                     The double precision arrays are not only more data to access, they
//                     require separate reads for X, Y, and Z coordinates and possibly charges
//                     as well.
//
// That map is repeated for each of N work units--the first 12*N sets of 32 numbers in the very
// large array filled by this function are object counts, warp instructions, and imported atom
// IDs for each of the work units.  The warp instructions are themselves bit-packed integers:
//
// [ 8 bits = object type ] [ 24 bits = warp number within the list of such objects ]
//
// The number encoded in the last 27 bits must be multiplied by 32 to get their actual start
// positions--like the atom IDs, the data substrates for each warp instruction are all
// multiples of 32 ints.
//
// Object type                              Warp Descriptions
// -----------     -------------------------------------------------------------------------
//    BOND         Atom IDs (uint), stretching and equilibrium length constants (double2)
//    ANGLE        Atom IDs (uint), stretching and equilibrium length constants (double2)
//   DIHEDRAL      Atom IDs (uint), multiplicities, amplitudes, and phase angles
//                   (three PMEFloats)
//    CMAP         Atom IDs I-J-K-L (uint) and energy map index plus atom ID M (uint)
//    QQXC         Atom IDs (uint), charges for all imported atoms (PMEFloat)
//    NB14         Atom IDs (uint), charges referenced from fields used by QQXC,
//                   Lennard-Jones ID numbers for all imported atoms (uint)
//  NMR(2,3,4)     Atom IDs (uint), r1 and r2 (double2), r3 and r4 (double2),
//                   k2 and k3 (double2)
//    UREY         Atom IDs (uint), stretching and equilibrium length constants (double2)
//    CIMP         Atom IDs (uint), stretching and equilibrium length constants (double2)
//
// It is worth pointing out that sets 1-4 are the only one read directly into __shared__ and
// kept there.  The others are all read as needed, 32 unsigned ints at a time, by the warp
// that needs them.
//
// Arguments:
//   gpu:        overarching data type containing all simulation parameters, including bond
//               work unit instruction arrays
//   bwunits:    the array of bond work units (this will be translated into arrays in gpu)
//   nbwunits:   length of bwunits
//   bwdims:     dimensions of the arrays allocated as well as other counts for GPU activities
//---------------------------------------------------------------------------------------------
void MakeBondedUnitDirections(gpuContext gpu, bondwork *bwunits, int nbwunits, bwalloc *bwdims)
{
  int i, j;
  unsigned int *uiptr, *ljid;
  PMEDouble2 *d2ptr;
  PMEFloat *fptr;
  PMEFloat2 *f2ptr;

  // Offsets for each type of object
  int bondUintPos    = nbwunits * 3 * BOND_WORK_UNIT_THREADS_PER_BLOCK;
  int anglUintPos    = bondUintPos +    (GRID * bwdims->nbondwarps);
  int diheUintPos    = anglUintPos +    (GRID * bwdims->nanglwarps);
  int cmapUintPos    = diheUintPos +    (GRID * bwdims->ndihewarps);
  int qqxcUintPos    = cmapUintPos +  2*(GRID * bwdims->ncmapwarps);
  int nb14UintPos    = qqxcUintPos +    (GRID * bwdims->nqqxcwarps);
  int nmr2UintPos    = nb14UintPos +    (GRID * bwdims->nnb14warps);
  int nmr3UintPos    = nmr2UintPos +  4*(GRID * bwdims->nnmr2warps);
  int nmr4UintPos    = nmr3UintPos +  4*(GRID * bwdims->nnmr3warps);
  int ureyUintPos    = nmr4UintPos +  4*(GRID * bwdims->nnmr4warps);
  int cimpUintPos    = ureyUintPos +    (GRID * bwdims->nureywarps);
  int cnstUintPos    = cimpUintPos +    (GRID * bwdims->ncimpwarps);
  int cnstUpdatePos  = cnstUintPos +    (GRID * bwdims->ncnstwarps);
  int bondDbl2Pos    = 0;
  int anglDbl2Pos    = bondDbl2Pos +    (GRID * bwdims->nbondwarps);
  int nmr2Dbl2Pos    = anglDbl2Pos +    (GRID * bwdims->nanglwarps);
  int nmr3Dbl2Pos    = nmr2Dbl2Pos +  6*(GRID * bwdims->nnmr2warps);
  int nmr4Dbl2Pos    = nmr3Dbl2Pos +  6*(GRID * bwdims->nnmr3warps);
  int ureyDbl2Pos    = nmr4Dbl2Pos +  6*(GRID * bwdims->nnmr4warps);
  int cimpDbl2Pos    = ureyDbl2Pos +    (GRID * bwdims->nureywarps);
  int cnstDbl2Pos    = cimpDbl2Pos +    (GRID * bwdims->ncimpwarps);
  int qPmefPos       = 0;
  int dihePmefPos    = nbwunits * BOND_WORK_UNIT_THREADS_PER_BLOCK;
  int nb14PmefPos    = dihePmefPos +    (GRID * bwdims->ndihewarps);
  int dihePmef2Pos   = 0;
  int nb14Pmef2Pos   = dihePmef2Pos + 2*(GRID * bwdims->ndihewarps);
  int cnstDbl2Base   = cnstDbl2Pos;
  int bondStatusPos = 0, anglStatusPos = 0, diheStatusPos = 0, cmapStatusPos = 0,
      qqxcStatusPos = 0, nb14StatusPos = 0, nmr2StatusPos = 0, nmr3StatusPos = 0,
      nmr4StatusPos = 0, ureyStatusPos = 0, cimpStatusPos = 0, cnstStatusPos = 0;
  if (gpu->sim.ti_mode > 0) {
    bondStatusPos = cnstUpdatePos + (GRID * bwdims->ncnstwarps);
    anglStatusPos = bondStatusPos + (GRID * bwdims->nbondwarps);
    diheStatusPos = anglStatusPos + (GRID * bwdims->nanglwarps);
    cmapStatusPos = diheStatusPos + (GRID * bwdims->ndihewarps);
    qqxcStatusPos = cmapStatusPos + (GRID * bwdims->ncmapwarps);
    nb14StatusPos = qqxcStatusPos + (GRID * bwdims->nqqxcwarps);
    nmr2StatusPos = nb14StatusPos + (GRID * bwdims->nnb14warps);
    nmr3StatusPos = nmr2StatusPos + (GRID * bwdims->nnmr2warps);
    nmr4StatusPos = nmr3StatusPos + (GRID * bwdims->nnmr3warps);
    ureyStatusPos = nmr4StatusPos + (GRID * bwdims->nnmr4warps);
    cimpStatusPos = ureyStatusPos + (GRID * bwdims->nureywarps);
    cnstStatusPos = cimpStatusPos + (GRID * bwdims->ncimpwarps);
  }

  // Counters for the total number of warps operating on each object across all work units
  unsigned int Tbondwarp = 0;
  unsigned int Tanglwarp = 0;
  unsigned int Tdihewarp = 0;
  unsigned int Tcmapwarp = 0;
  unsigned int Tqqxcwarp = 0;
  unsigned int Tnb14warp = 0;
  unsigned int Tnmr2warp = 0;
  unsigned int Tnmr3warp = 0;
  unsigned int Tnmr4warp = 0;
  unsigned int Tureywarp = 0;
  unsigned int Tcimpwarp = 0;
  unsigned int Tcnstwarp = 0;

  // Set pointers
  uiptr = gpu->pbBondWorkUnitUINT->_pSysData;
  d2ptr = gpu->pbBondWorkUnitDBL2->_pSysData;
  fptr  = gpu->pbBondWorkUnitPFLOAT->_pSysData;
  f2ptr = gpu->pbBondWorkUnitPFLOAT2->_pSysData;

  // Loop over all bond work units, incrementing offsets in the process.
  for (i = 0; i < nbwunits; i++) {

    // Load the atom list array.  Start with the basic instructions,
    // then load atom IDs.  Load 0xffffffff for every blank in the
    // atom list so that the device will know not to go looking for
    // anything in that slot during the neighbor list remapping.
    int ioffset = (3*i + 1) * BOND_WORK_UNIT_THREADS_PER_BLOCK;
    for (j = 0; j < bwunits[i].natom; j++) {
      uiptr[ioffset + j] = bwunits[i].atomList[j];
    }
    for (j = bwunits[i].natom; j < BOND_WORK_UNIT_THREADS_PER_BLOCK; j++) {
      uiptr[ioffset + j] = 0xffffffff;
    }

    // Bit-pack the atom IDs and statuses for each object
    for (j = 0; j < bwunits[i].nbond; j++) {
      if (bwunits[i].bondList[j] < 0) {
        continue;
      }
      unsigned int packID = (bwunits[i].bondMap[2*j] << 8) | (bwunits[i].bondMap[2*j+1]);
      uiptr[bondUintPos + j] = packID;
      if (gpu->sim.ti_mode > 0) {
        uiptr[bondStatusPos + j] = bwunits[i].bondStatus[j];
      }
    }
    for (j = 0; j < bwunits[i].nangl; j++) {
      if (bwunits[i].bondAngleList[j] < 0) {
        continue;
      }
      unsigned int packID = (bwunits[i].bondAngleMap[3*j    ] << 16) |
                            (bwunits[i].bondAngleMap[3*j + 1] << 8) |
                            (bwunits[i].bondAngleMap[3*j + 2]);
      uiptr[anglUintPos + j] = packID;
      if (gpu->sim.ti_mode > 0) {
        uiptr[anglStatusPos + j] = bwunits[i].anglStatus[j];
      }
    }
    for (j = 0; j < bwunits[i].ndihe; j++) {
      if (bwunits[i].dihedralList[j] < 0) {
        continue;
      }
      unsigned int packID = (bwunits[i].dihedralMap[4*j    ] << 24) |
                            (bwunits[i].dihedralMap[4*j + 1] << 16) |
                            (bwunits[i].dihedralMap[4*j + 2] << 8) |
                            (bwunits[i].dihedralMap[4*j + 3]);
      uiptr[diheUintPos + j] = packID;
      if (gpu->sim.ti_mode > 0) {
        uiptr[diheStatusPos + j] = bwunits[i].diheStatus[j];
      }
    }
    for (j = 0; j < bwunits[i].ncmap; j++) {
      if (bwunits[i].cmapList[j] < 0) {
        continue;
      }
      unsigned int packID = (bwunits[i].cmapMap[5*j    ] << 24) |
                            (bwunits[i].cmapMap[5*j + 1] << 16) |
                            (bwunits[i].cmapMap[5*j + 2] << 8) |
                            (bwunits[i].cmapMap[5*j + 3]);
     unsigned int tgx = j & GRID_BITS_MASK;
      int warpStep = 2 * (j >> GRID_BITS) * GRID;
      uiptr[cmapUintPos + tgx + warpStep] = packID;

      // Note that this second bit-packed unsigned int takes the energy map ID first, then
      // the atom M ID--this allows virtually unlimited numbers of unique CMAP surfaces.
      packID = (gpu->pbCmapType->_pSysData[bwunits[i].cmapList[j]] << 8) |
               (bwunits[i].cmapMap[5*j + 4]);
      uiptr[cmapUintPos + tgx + warpStep + GRID] = packID;
      if (gpu->sim.ti_mode > 0) {
        uiptr[cmapStatusPos + j] = bwunits[i].cmapStatus[j];
      }
    }
    for (j = 0; j < bwunits[i].nqqxc; j++) {
      if (bwunits[i].qqxcList[j] < 0) {
        continue;
      }
      unsigned int packID = (bwunits[i].qqxcMap[2*j] << 8) | (bwunits[i].qqxcMap[2*j + 1]);
      uiptr[qqxcUintPos + j] = packID;
      if (gpu->sim.ti_mode > 0) {
        uiptr[qqxcStatusPos + j] = bwunits[i].qqxcStatus[j];
      }
    }
    for (j = 0; j < bwunits[i].nnb14; j++) {
      if (bwunits[i].nb14List[j] < 0) {
        continue;
      }
      unsigned int packID = (bwunits[i].nb14Map[2*j] << 8) | (bwunits[i].nb14Map[2*j + 1]);
      uiptr[nb14UintPos + j] = packID;
      if (gpu->sim.ti_mode > 0) {
        uiptr[nb14StatusPos + j] = bwunits[i].nb14Status[j];
      }
    }
    for (j = 0; j < bwunits[i].nnmr2; j++) {
      if (bwunits[i].nmr2List[j] < 0) {
        continue;
      }
      unsigned int packID = (bwunits[i].nmr2Map[2*j] << 8) | (bwunits[i].nmr2Map[2*j + 1]);
     unsigned int tgx = j & GRID_BITS_MASK;
      int warpStep = nmr2UintPos + (4 * (j >> GRID_BITS) * GRID) + tgx;
      uiptr[warpStep         ] = packID;
      uiptr[warpStep +   GRID] = gpu->pbNMRDistanceInc->_pSysData[bwunits[i].nmr2List[j]];
      uiptr[warpStep + 2*GRID] = gpu->pbNMRDistanceStep->_pSysData[bwunits[i].nmr2List[j]].x;
      uiptr[warpStep + 3*GRID] = gpu->pbNMRDistanceStep->_pSysData[bwunits[i].nmr2List[j]].y;
      if (gpu->sim.ti_mode > 0) {
        uiptr[nmr2StatusPos + j] = bwunits[i].nmr2Status[j];
      }
    }
    for (j = 0; j < bwunits[i].nnmr3; j++) {
      if (bwunits[i].nmr3List[j] < 0) {
        continue;
      }
      unsigned int packID = (bwunits[i].nmr3Map[3*j    ] << 16) |
                            (bwunits[i].nmr3Map[3*j + 1] << 8) |
                            (bwunits[i].nmr3Map[3*j + 2]);
     unsigned int tgx = j & GRID_BITS_MASK;
      int warpStep = nmr3UintPos + (4 * (j >> GRID_BITS) * GRID) + tgx;
      uiptr[warpStep         ] = packID;
      uiptr[warpStep +   GRID] = gpu->pbNMRAngleInc->_pSysData[bwunits[i].nmr3List[j]];
      uiptr[warpStep + 2*GRID] = gpu->pbNMRAngleStep->_pSysData[bwunits[i].nmr3List[j]].x;
      uiptr[warpStep + 3*GRID] = gpu->pbNMRAngleStep->_pSysData[bwunits[i].nmr3List[j]].y;
      if (gpu->sim.ti_mode > 0) {
        uiptr[nmr3StatusPos + j] = bwunits[i].nmr3Status[j];
      }
    }
    for (j = 0; j < bwunits[i].nnmr4; j++) {
      if (bwunits[i].nmr4List[j] < 0) {
        continue;
      }
      unsigned int packID = (bwunits[i].nmr4Map[4*j    ] << 24) |
                            (bwunits[i].nmr4Map[4*j + 1] << 16) |
                            (bwunits[i].nmr4Map[4*j + 2] << 8) |
                            (bwunits[i].nmr4Map[4*j + 3]);
     unsigned int tgx = j & GRID_BITS_MASK;
      int warpStep = nmr4UintPos + (4 * (j >> GRID_BITS) * GRID) + tgx;
      uiptr[warpStep         ] = packID;
      uiptr[warpStep +   GRID] = gpu->pbNMRTorsionInc->_pSysData[bwunits[i].nmr4List[j]];
      uiptr[warpStep + 2*GRID] = gpu->pbNMRTorsionStep->_pSysData[bwunits[i].nmr4List[j]].x;
      uiptr[warpStep + 3*GRID] = gpu->pbNMRTorsionStep->_pSysData[bwunits[i].nmr4List[j]].y;
      if (gpu->sim.ti_mode > 0) {
        uiptr[nmr4StatusPos + j] = bwunits[i].nmr4Status[j];
      }
    }
    for (j = 0; j < bwunits[i].nurey; j++) {
      if (bwunits[i].ureyList[j] < 0) {
        continue;
      }
      unsigned int packID = (bwunits[i].ureyMap[2*j] << 8) | (bwunits[i].ureyMap[2*j + 1]);
      uiptr[ureyUintPos + j] = packID;
      if (gpu->sim.ti_mode > 0) {
        uiptr[ureyStatusPos + j] = bwunits[i].ureyStatus[j];
      }
    }
    for (j = 0; j < bwunits[i].ncimp; j++) {
      if (bwunits[i].cimpList[j] < 0) {
        continue;
      }
      unsigned int packID = (bwunits[i].cimpMap[4*j    ] << 24) |
                            (bwunits[i].cimpMap[4*j + 1] << 16) |
                            (bwunits[i].cimpMap[4*j + 2] << 8) |
                            (bwunits[i].cimpMap[4*j + 3]);
      uiptr[cimpUintPos + j] = packID;
      if (gpu->sim.ti_mode > 0) {
        uiptr[cimpStatusPos + j] = bwunits[i].cimpStatus[j];
      }
    }
    for (j = 0; j < bwunits[i].ncnst; j++) {
      if (bwunits[i].cnstList[j] < 0) {
        continue;
      }
      unsigned int packID = bwunits[i].cnstMap[j];
      uiptr[cnstUintPos + j] = packID;

      // Special case!  Record the index of the bond work units' compiled PMEDouble2 array
      // where the parameters for this atom position constraint are stored (the constraint's
      // ID in the master topology is stored in bwunits[i].cnstList[j]).  This permits us to
      // access pBwuCnstUpdateIdx[constraint number] during rescaling in order to update the
      // constraints' reference positions in pBwuCnst.  Note also that the local variable
      // cnstUpdatePos is not incremented along with the other counters--it is a fixed
      // reference point for all work units.
     unsigned int tgx = j & GRID_BITS_MASK;
      int warpStep = 2 * (j >> GRID_BITS) * GRID;
      uiptr[cnstUpdatePos + bwunits[i].cnstList[j]] = cnstDbl2Pos - cnstDbl2Base +
                                                      warpStep + tgx;
      if (gpu->sim.ti_mode > 0) {
        uiptr[cnstStatusPos + j] = bwunits[i].cnstStatus[j];
      }
    }

    // Now pack the PMEDouble2 parameter details of each object
    for (j = 0; j < bwunits[i].nbond; j++) {
      if (bwunits[i].bondList[j] < 0) {
        continue;
      }
      d2ptr[bondDbl2Pos + j] = gpu->pbBond->_pSysData[bwunits[i].bondList[j]];
    }
    for (j = 0; j < bwunits[i].nangl; j++) {
      if (bwunits[i].bondAngleList[j] < 0) {
        continue;
      }
      d2ptr[anglDbl2Pos + j] = gpu->pbBondAngle->_pSysData[bwunits[i].bondAngleList[j]];
    }
    for (j = 0; j < bwunits[i].nnmr2; j++) {
      if (bwunits[i].nmr2List[j] < 0) {
        continue;
      }
     unsigned int tgx = j & GRID_BITS_MASK;
      int warpStep = nmr2Dbl2Pos + (6 * (j >> GRID_BITS) * GRID) + tgx;
      int objid = bwunits[i].nmr2List[j];
      d2ptr[warpStep         ] = gpu->pbNMRDistanceR1R2Slp->_pSysData[objid];
      d2ptr[warpStep +   GRID] = gpu->pbNMRDistanceR3R4Slp->_pSysData[objid];
      d2ptr[warpStep + 2*GRID] = gpu->pbNMRDistanceK2K3Slp->_pSysData[objid];
      d2ptr[warpStep + 3*GRID] = gpu->pbNMRDistanceR1R2Int->_pSysData[objid];
      d2ptr[warpStep + 4*GRID] = gpu->pbNMRDistanceR3R4Int->_pSysData[objid];
      d2ptr[warpStep + 5*GRID] = gpu->pbNMRDistanceK2K3Int->_pSysData[objid];
    }
    for (j = 0; j < bwunits[i].nnmr3; j++) {
      if (bwunits[i].nmr3List[j] < 0) {
        continue;
      }
     unsigned int tgx = j & GRID_BITS_MASK;
      int warpStep = nmr3Dbl2Pos + (6 * (j >> GRID_BITS) * GRID) + tgx;
      int objid = bwunits[i].nmr3List[j];
      d2ptr[warpStep         ] = gpu->pbNMRAngleR1R2Slp->_pSysData[objid];
      d2ptr[warpStep +   GRID] = gpu->pbNMRAngleR3R4Slp->_pSysData[objid];
      d2ptr[warpStep + 2*GRID] = gpu->pbNMRAngleK2K3Slp->_pSysData[objid];
      d2ptr[warpStep + 3*GRID] = gpu->pbNMRAngleR1R2Int->_pSysData[objid];
      d2ptr[warpStep + 4*GRID] = gpu->pbNMRAngleR3R4Int->_pSysData[objid];
      d2ptr[warpStep + 5*GRID] = gpu->pbNMRAngleK2K3Int->_pSysData[objid];
    }
    for (j = 0; j < bwunits[i].nnmr4; j++) {
      if (bwunits[i].nmr4List[j] < 0) {
        continue;
      }
     unsigned int tgx = j & GRID_BITS_MASK;
      int warpStep = nmr4Dbl2Pos + (6 * (j >> GRID_BITS) * GRID) + tgx;
      int objid = bwunits[i].nmr4List[j];
      d2ptr[warpStep         ] = gpu->pbNMRTorsionR1R2Slp->_pSysData[objid];
      d2ptr[warpStep +   GRID] = gpu->pbNMRTorsionR3R4Slp->_pSysData[objid];
      d2ptr[warpStep + 2*GRID] = gpu->pbNMRTorsionK2K3Slp->_pSysData[objid];
      d2ptr[warpStep + 3*GRID] = gpu->pbNMRTorsionR1R2Int->_pSysData[objid];
      d2ptr[warpStep + 4*GRID] = gpu->pbNMRTorsionR3R4Int->_pSysData[objid];
      d2ptr[warpStep + 5*GRID] = gpu->pbNMRTorsionK2K3Int->_pSysData[objid];
    }
    for (j = 0; j < bwunits[i].nurey; j++) {
      if (bwunits[i].ureyList[j] < 0) {
        continue;
      }
      d2ptr[ureyDbl2Pos + j] = gpu->pbUBAngle->_pSysData[bwunits[i].ureyList[j]];
    }
    for (j = 0; j < bwunits[i].ncimp; j++) {
      if (bwunits[i].cimpList[j] < 0) {
        continue;
      }
      d2ptr[cimpDbl2Pos + j] = gpu->pbImpDihedral->_pSysData[bwunits[i].cimpList[j]];
    }
    for (j = 0; j < bwunits[i].ncnst; j++) {
      if (bwunits[i].cnstList[j] < 0) {
        continue;
      }
     unsigned int tgx = j & GRID_BITS_MASK;
      int warpStep = 2 * (j >> GRID_BITS) * GRID;
      int objid = bwunits[i].cnstList[j];
      d2ptr[cnstDbl2Pos + warpStep + tgx       ] = gpu->pbConstraint1->_pSysData[objid];
      d2ptr[cnstDbl2Pos + warpStep + tgx + GRID] = gpu->pbConstraint2->_pSysData[objid];
    }

    // Pack the PMEFloat details of each work unit
    for (j = 0; j < bwunits[i].natom; j++) {
      fptr[qPmefPos + j] = gpu->pbAtomCharge->_pSysData[bwunits[i].atomList[j]];
    }
    for (j = 0; j < bwunits[i].ndihe; j++) {
      if (bwunits[i].dihedralList[j] < 0) {
        continue;
      }
      fptr[dihePmefPos + j] = gpu->pbDihedral3->_pSysData[bwunits[i].dihedralList[j]];
    }
    for (j = 0; j < bwunits[i].nnb14; j++) {
      if (bwunits[i].nb14List[j] < 0) {
        continue;
      }
      fptr[nb14PmefPos + j] = gpu->pbNb141->_pSysData[bwunits[i].nb14List[j]].x;
    }

    // Pack the PMEFloat2 details of each work unit
    for (j = 0; j < bwunits[i].nnb14; j++) {
      if (bwunits[i].nb14List[j] < 0) {
        continue;
      }
      f2ptr[nb14Pmef2Pos + j].x = gpu->pbNb142->_pSysData[bwunits[i].nb14List[j]].x;
      f2ptr[nb14Pmef2Pos + j].y = gpu->pbNb142->_pSysData[bwunits[i].nb14List[j]].y;
      double dampener = gpu->pbNb141->_pSysData[bwunits[i].nb14List[j]].y;
      f2ptr[nb14Pmef2Pos + j].x *= dampener;
      f2ptr[nb14Pmef2Pos + j].y *= dampener;
    }
    for (j = 0; j < bwunits[i].ndihe; j++) {
      if (bwunits[i].dihedralList[j] < 0) {
        continue;
      }
     unsigned int tgx = j & GRID_BITS_MASK;
      int warpStep = 2 * (j >> GRID_BITS) * GRID;
      int objid = bwunits[i].dihedralList[j];
      f2ptr[dihePmef2Pos + warpStep + tgx       ] = gpu->pbDihedral1->_pSysData[objid];
      f2ptr[dihePmef2Pos + warpStep + tgx + GRID] = gpu->pbDihedral2->_pSysData[objid];
    }

    // Finally, write the warp instructions
    int nWUbondwarp = (bwunits[i].nbond + GRID_BITS_MASK) >> GRID_BITS;
    int nWUanglwarp = (bwunits[i].nangl + GRID_BITS_MASK) >> GRID_BITS;
    int nWUdihewarp = (bwunits[i].ndihe + GRID_BITS_MASK) >> GRID_BITS;
    int nWUcmapwarp = (bwunits[i].ncmap + GRID_BITS_MASK) >> GRID_BITS;
    int nWUqqxcwarp = (bwunits[i].nqqxc + GRID_BITS_MASK) >> GRID_BITS;
    int nWUnb14warp = (bwunits[i].nnb14 + GRID_BITS_MASK) >> GRID_BITS;
    int nWUnmr2warp = (bwunits[i].nnmr2 + GRID_BITS_MASK) >> GRID_BITS;
    int nWUnmr3warp = (bwunits[i].nnmr3 + GRID_BITS_MASK) >> GRID_BITS;
    int nWUnmr4warp = (bwunits[i].nnmr4 + GRID_BITS_MASK) >> GRID_BITS;
    int nWUureywarp = (bwunits[i].nurey + GRID_BITS_MASK) >> GRID_BITS;
    int nWUcimpwarp = (bwunits[i].ncimp + GRID_BITS_MASK) >> GRID_BITS;
    int nWUcnstwarp = (bwunits[i].ncnst + GRID_BITS_MASK) >> GRID_BITS;

    // Just put big things at the front.  There are a lot of counters in use
    // here: WUtask is a counter of the number of tasks in this particular work
    // unit.  T<object>warp is a counter over all <object> dealt out across all
    // work units. (T<object>warp is needed to index into the contiguous arrays
    // of imported atom IDs and parameters that all work units draw upon, so
    // that each work unit can access IDs of imported atoms that make sense
    // and get parameters applicable to them.) Meanwhile, nWU<object>warp is
    // the total number of warps that the work unit must devote to <object>,
    // limits that various counters will eventually reach.
    int WUtask = 0;
    ioffset = (3 * i * BOND_WORK_UNIT_THREADS_PER_BLOCK) + WARP_INSTRUCTION_OFFSET;
    for (j = 0; j < nWUcmapwarp; j++) {
      uiptr[ioffset + WUtask] = (CMAP_CODE << 24) | Tcmapwarp;
      Tcmapwarp++;
      WUtask++;
    }
    for (j = 0; j < nWUdihewarp; j++) {
      uiptr[ioffset + WUtask] = (DIHE_CODE << 24) | Tdihewarp;
      Tdihewarp++;
      WUtask++;
    }
    for (j = 0; j < nWUcimpwarp; j++) {
      uiptr[ioffset + WUtask] = (CIMP_CODE << 24) | Tcimpwarp;
      Tcimpwarp++;
      WUtask++;
    }
    for (j = 0; j < nWUnmr4warp; j++) {
      uiptr[ioffset + WUtask] = (NMR4_CODE << 24) | Tnmr4warp;
      Tnmr4warp++;
      WUtask++;
    }
    for (j = 0; j < nWUanglwarp; j++) {
      uiptr[ioffset + WUtask] = (ANGL_CODE << 24) | Tanglwarp;
      Tanglwarp++;
      WUtask++;
    }
    for (j = 0; j < nWUnmr3warp; j++) {
      uiptr[ioffset + WUtask] = (NMR3_CODE << 24) | Tnmr3warp;
      Tnmr3warp++;
      WUtask++;
    }
    for (j = 0; j < nWUnmr2warp; j++) {
      uiptr[ioffset + WUtask] = (NMR2_CODE << 24) | Tnmr2warp;
      Tnmr2warp++;
      WUtask++;
    }
    for (j = 0; j < nWUnb14warp; j++) {
      uiptr[ioffset + WUtask] = (NB14_CODE << 24) | Tnb14warp;
      Tnb14warp++;
      WUtask++;
    }
    for (j = 0; j < nWUqqxcwarp; j++) {
      uiptr[ioffset + WUtask] = (QQXC_CODE << 24) | Tqqxcwarp;
      Tqqxcwarp++;
      WUtask++;
    }
    for (j = 0; j < nWUbondwarp; j++) {
      uiptr[ioffset + WUtask] = (BOND_CODE << 24) | Tbondwarp;
      Tbondwarp++;
      WUtask++;
    }
    for (j = 0; j < nWUureywarp; j++) {
      uiptr[ioffset + WUtask] = (UREY_CODE << 24) | Tureywarp;
      Tureywarp++;
      WUtask++;
    }
    for (j = 0; j < nWUcnstwarp; j++) {
      uiptr[ioffset + WUtask] = (CNST_CODE << 24) | Tcnstwarp;
      Tcnstwarp++;
      WUtask++;
    }

    // Write out the preliminary instructions for the work unit
    ioffset = 3 * i * BOND_WORK_UNIT_THREADS_PER_BLOCK;
    uiptr[ioffset + ATOM_IMPORT_COUNT_IDX]      = bwunits[i].natom;
    uiptr[ioffset + ATOM_IMPORT_SOURCE_IDX]     = ATOM_DP_SOURCE_CODE;
    uiptr[ioffset + WARP_INSTRUCTION_COUNT_IDX] = WUtask;
    uiptr[ioffset + FORCE_ACCUMULATOR_IDX]      = bwunits[i].frcAcc;

    // Record where in the common arrays this work unit will read
    bwunits[i].bondDbl2Idx    = bondDbl2Pos;
    bwunits[i].anglDbl2Idx    = anglDbl2Pos;
    bwunits[i].nmr2Dbl2Idx    = nmr2Dbl2Pos;
    bwunits[i].nmr3Dbl2Idx    = nmr3Dbl2Pos;
    bwunits[i].nmr4Dbl2Idx    = nmr4Dbl2Pos;
    bwunits[i].ureyDbl2Idx    = ureyDbl2Pos;
    bwunits[i].cimpDbl2Idx    = cimpDbl2Pos;
    bwunits[i].cnstDbl2Idx    = cnstDbl2Pos;
    bwunits[i].dihePFloatIdx  = dihePmefPos;
    bwunits[i].nb14PFloatIdx  = nb14PmefPos;
    bwunits[i].dihePFloat2Idx = dihePmef2Pos;
    bwunits[i].nb14PFloat2Idx = nb14Pmef2Pos;
    bwunits[i].qPFloatIdx     = qPmefPos;

    // Increment the offsets
    bondUintPos   +=     nWUbondwarp * GRID;
    anglUintPos   +=     nWUanglwarp * GRID;
    diheUintPos   +=     nWUdihewarp * GRID;
    cmapUintPos   += 2 * nWUcmapwarp * GRID;
    qqxcUintPos   +=     nWUqqxcwarp * GRID;
    nb14UintPos   +=     nWUnb14warp * GRID;
    nmr2UintPos   += 4 * nWUnmr2warp * GRID;
    nmr3UintPos   += 4 * nWUnmr3warp * GRID;
    nmr4UintPos   += 4 * nWUnmr4warp * GRID;
    ureyUintPos   +=     nWUureywarp * GRID;
    cimpUintPos   +=     nWUcimpwarp * GRID;
    cnstUintPos   +=     nWUcnstwarp * GRID;
    bondDbl2Pos   +=     nWUbondwarp * GRID;
    anglDbl2Pos   +=     nWUanglwarp * GRID;
    nmr2Dbl2Pos   += 6 * nWUnmr2warp * GRID;
    nmr3Dbl2Pos   += 6 * nWUnmr3warp * GRID;
    nmr4Dbl2Pos   += 6 * nWUnmr4warp * GRID;
    ureyDbl2Pos   +=     nWUureywarp * GRID;
    cimpDbl2Pos   +=     nWUcimpwarp * GRID;
    cnstDbl2Pos   += 2 * nWUcnstwarp * GRID;
    qPmefPos      +=     BOND_WORK_UNIT_THREADS_PER_BLOCK;
    dihePmefPos   +=     nWUdihewarp * GRID;
    nb14PmefPos   +=     nWUnb14warp * GRID;
    dihePmef2Pos  += 2 * nWUdihewarp * GRID;
    nb14Pmef2Pos  +=     nWUnb14warp * GRID;
    bondStatusPos +=     nWUbondwarp * GRID;
    anglStatusPos +=     nWUanglwarp * GRID;
    diheStatusPos +=     nWUdihewarp * GRID;
    cmapStatusPos +=     nWUcmapwarp * GRID;
    qqxcStatusPos +=     nWUqqxcwarp * GRID;
    nb14StatusPos +=     nWUnb14warp * GRID;
    nmr2StatusPos +=     nWUnmr2warp * GRID;
    nmr3StatusPos +=     nWUnmr3warp * GRID;
    nmr4StatusPos +=     nWUnmr4warp * GRID;
    ureyStatusPos +=     nWUureywarp * GRID;
    cimpStatusPos +=     nWUcimpwarp * GRID;
    cnstStatusPos +=     nWUcnstwarp * GRID;
  }
}

//---------------------------------------------------------------------------------------------
// MatchObjectByWorkUnitData: accept the details of a bonded term object--its type, ID numbers
//                            of its atoms, and force constants--and match those to the first
//                            available instance of such as object, then check that off of a
//                            list of object assignments.  There are edge cases here, namely
//                            that there could be some force field that happens to assign two
//                            identical yet separate terms to the same tuple of atoms.  In this
//                            case, the goal is merely to make sure that all terms have been
//                            assigned, so always checking off the first unassigned term while
//                            looping over all work units will find them all if the work units
//                            have been built properly.
//
// Arguments:
//   gpu:           overarching data structure containing the tranlsated bond work units and
//                  all other simulation parameters
//   objcode:       numerical code for each bonded term type
//   atomIDBuffer:  the ID numbers of all atoms in the object--holds up to five
//   prmBuffer:     parameters for the object--holds all sorts of details.  The object type
//                  will dictate which attributes must be matched.
//---------------------------------------------------------------------------------------------
static int MatchObjectByWorkUnitData(gpuContext gpu, int objcode, aidbuff atomIDBuffer,
                                     prmbuff prmBuffer, imat *objAssigned, int* p1atm,
                                     int2* p2atm, int4* p4atm)
{
  int i;
  int *paramI;
  PMEDouble *paramD;
  PMEDouble2 *paramD2a;
  PMEDouble2 *paramD2b;
  PMEDouble2 *paramD2c;
  PMEDouble2 *paramD2d;
  PMEDouble2 *paramD2e;
  PMEDouble2 *paramD2f;
  PMEFloat *paramF;
  PMEFloat2 *paramF2a;
  PMEFloat2 *paramF2b;

  // Set pointers to the appropriate details.  No parameters will be checked for electrostatic
  // exclusions--if the atom IDs are correct, that's it for this test (charges are checked
  // elsewhere).
  if (objcode == BOND_CODE) {
    paramD2a = gpu->pbBond->_pSysData;
  }
  else if (objcode == ANGL_CODE) {
    paramD2a = gpu->pbBondAngle->_pSysData;
  }
  else if (objcode == DIHE_CODE) {
    paramF2a = gpu->pbDihedral1->_pSysData;
    paramF2b = gpu->pbDihedral2->_pSysData;
    paramF = gpu->pbDihedral3->_pSysData;
  }
  else if (objcode == CMAP_CODE) {
    paramI = gpu->pbCmapType->_pSysData;
  }
  else if (objcode == NB14_CODE) {
    paramD2a = gpu->pbNb141->_pSysData;
    paramD2b = gpu->pbNb142->_pSysData;
  }
  else if (objcode == NMR2_CODE) {
    paramD2a = gpu->pbNMRDistanceR1R2Slp->_pSysData;
    paramD2b = gpu->pbNMRDistanceR3R4Slp->_pSysData;
    paramD2c = gpu->pbNMRDistanceK2K3Slp->_pSysData;
    paramD2d = gpu->pbNMRDistanceR1R2Int->_pSysData;
    paramD2e = gpu->pbNMRDistanceR3R4Int->_pSysData;
    paramD2f = gpu->pbNMRDistanceK2K3Int->_pSysData;
  }
  else if (objcode == NMR3_CODE) {
    paramD2a = gpu->pbNMRAngleR1R2Slp->_pSysData;
    paramD2b = gpu->pbNMRAngleR3R4Slp->_pSysData;
    paramD2c = gpu->pbNMRAngleK2K3Slp->_pSysData;
    paramD2d = gpu->pbNMRAngleR1R2Int->_pSysData;
    paramD2e = gpu->pbNMRAngleR3R4Int->_pSysData;
    paramD2f = gpu->pbNMRAngleK2K3Int->_pSysData;
  }
  else if (objcode == NMR4_CODE) {
    paramD2a = gpu->pbNMRTorsionR1R2Slp->_pSysData;
    paramD2b = gpu->pbNMRTorsionR3R4Slp->_pSysData;
    paramD2c = gpu->pbNMRTorsionK2K3Slp->_pSysData;
    paramD2d = gpu->pbNMRTorsionR1R2Int->_pSysData;
    paramD2e = gpu->pbNMRTorsionR3R4Int->_pSysData;
    paramD2f = gpu->pbNMRTorsionK2K3Int->_pSysData;
  }
  else if (objcode == UREY_CODE) {
    paramD2a = gpu->pbUBAngle->_pSysData;
  }
  else if (objcode == CIMP_CODE) {
    paramD2a = gpu->pbImpDihedral->_pSysData;
  }
  else if (objcode == CNST_CODE) {
    paramD2a = gpu->pbConstraint1->_pSysData;
    paramD2b = gpu->pbConstraint2->_pSysData;
  }
  int nobj = GetSimulationObjectCount(gpu, objcode, " ");
  bool matched;
  for (i = 0; i < nobj; i++) {

    // Check for previous assignment
    if (objAssigned->data[i] == 1) {
      continue;
    }

    // Check that the atom IDs clear
    matched = false;
    if (objcode == BOND_CODE || objcode == NMR2_CODE || objcode == NMR3_CODE ||
        objcode == QQXC_CODE || objcode == NB14_CODE || objcode == UREY_CODE) {
      matched = ((atomIDBuffer.refI == p2atm[i].x && atomIDBuffer.refJ == p2atm[i].y) ||
                 (atomIDBuffer.refI == p2atm[i].y && atomIDBuffer.refJ == p2atm[i].x));
    }
    else if (objcode == ANGL_CODE || objcode == NMR3_CODE) {
      matched = (((atomIDBuffer.refI == p2atm[i].x && atomIDBuffer.refK == p1atm[i]) ||
                  (atomIDBuffer.refI == p1atm[i] && atomIDBuffer.refK == p2atm[i].x)) &&
                 atomIDBuffer.refJ == p2atm[i].y);
    }
    else if (objcode == DIHE_CODE || objcode == NMR4_CODE) {
      matched = ((atomIDBuffer.refI == p4atm[i].x && atomIDBuffer.refJ == p4atm[i].y &&
                  atomIDBuffer.refK == p4atm[i].z && atomIDBuffer.refL == p4atm[i].w) ||
                 (atomIDBuffer.refI == p4atm[i].w && atomIDBuffer.refJ == p4atm[i].z &&
                  atomIDBuffer.refK == p4atm[i].y && atomIDBuffer.refL == p4atm[i].x));
    }
    else if (objcode == CIMP_CODE || objcode == CMAP_CODE) {
      matched = (atomIDBuffer.refI == p4atm[i].x && atomIDBuffer.refJ == p4atm[i].y &&
                 atomIDBuffer.refK == p4atm[i].z && atomIDBuffer.refL == p4atm[i].w);
    }
    else if (objcode == CNST_CODE) {
      matched = (atomIDBuffer.refI == p1atm[i]);
    }
    if (matched && objcode == CMAP_CODE) {
      matched = (atomIDBuffer.refM == p1atm[i] && prmBuffer.i0 == paramI[i]);
    }
    if (matched == false) {
      continue;
    }

    // If the atom IDs have cleared, check the parameter details
    if (objcode == BOND_CODE || objcode == ANGL_CODE || objcode == NMR2_CODE ||
        objcode == NMR3_CODE || objcode == NMR4_CODE || objcode == UREY_CODE ||
        objcode == CIMP_CODE || objcode == CNST_CODE) {
      matched = ((fabs(paramD2a[i].x - prmBuffer.p1) < 1.0e-8 ||
                  fabs((paramD2a[i].x - prmBuffer.p1) / prmBuffer.p1) < 1.0e-6) &&
                 (fabs(paramD2a[i].y - prmBuffer.p2) < 1.0e-8 ||
                  fabs((paramD2a[i].y - prmBuffer.p2) / prmBuffer.p2) < 1.0e-6));
    }
    else if (objcode == DIHE_CODE) {
      matched = ((fabs(paramF2a[i].x - prmBuffer.p1) < 1.0e-8 ||
                  fabs((paramF2a[i].x - prmBuffer.p1) / prmBuffer.p1) < 1.0e-6) &&
                 (fabs(paramF2a[i].y - prmBuffer.p2) < 1.0e-8 ||
                  fabs((paramF2a[i].y - prmBuffer.p2) / prmBuffer.p2) < 1.0e-6) &&
                 (fabs(paramF2b[i].x - prmBuffer.p3) < 1.0e-8 ||
                  fabs((paramF2b[i].x - prmBuffer.p3) / prmBuffer.p3) < 1.0e-6) &&
                 (fabs(paramF2b[i].y - prmBuffer.p4) < 1.0e-8 ||
                  fabs((paramF2b[i].y - prmBuffer.p4) / prmBuffer.p4) < 1.0e-6) &&
                 (fabs(paramF[i] - prmBuffer.p5) < 1.0e-8 ||
                  fabs((paramF[i] - prmBuffer.p5) / prmBuffer.p5) < 1.0e-6));
    }
    if (matched == false) {
      continue;
    }
    if (objcode == CNST_CODE) {
      matched = ((fabs(paramD2b[i].x - prmBuffer.p3) < 1.0e-8 ||
                  fabs((paramD2b[i].x - prmBuffer.p3) / prmBuffer.p3) < 1.0e-6) &&
                 (fabs(paramD2b[i].y - prmBuffer.p4) < 1.0e-8 ||
                  fabs((paramD2b[i].y - prmBuffer.p4) / prmBuffer.p4) < 1.0e-6));
    }

    // Note the correspondence--allow for the possibility that one term
    // has been assigned to multiple work units simultaneously.
    if (matched) {
      objAssigned->data[i] += 1;
      break;
    }
  }

  // Return the result
  if (matched) {
    return 1;
  }
  else {
    return 0;
  }
}

//---------------------------------------------------------------------------------------------
// CheckWorkUnitTranslations: check the translations of each work unit, before they are
//                            uploaded to the GPU.  This will check that every bonded
//                            interaction object is assigned to a work units and that the
//                            atoms it access are the correct ones, after the rearrangements
//                            that have happened.  If this check passes, the GPU should be
//                            ready to process the work units and get the right answers.
//
// Arguments:
//   gpu:      overarching data structure containing the tranlsated bond work units and all
//             other simulation parameters (this function is checking whether all of that
//             data is consistent)
//   bwunits:  the bond work units as originally formulated (imported for CalcBondBlockSize--
//             otherwise this data is not referenced)
//   objtype:  the object type
//
// This is a debugging function.
//---------------------------------------------------------------------------------------------
void CheckWorkUnitTranslations(gpuContext gpu, bondwork* bwunits, const char* objtype)
{
  int i, j, k, objcode=0;

  // Return immediately if there are no objects to map--don't
  // get caught mucking around with unassigned pointers.
  int nobj = GetSimulationObjectCount(gpu, -1, objtype);
  if (nobj == 0 || gpu->sim.bondWorkUnits == 0) {
    return;
  }

  // Set up pointers, as in other functions.
  int *p1atm = NULL;
  int2 *p2atm = NULL;
  int4 *p4atm = NULL;
  codify_object_type(objtype, gpu, p1atm, p2atm, p4atm, objcode);

  // Set up more pointers, here a mock-up of the pointers that exist to GPU device
  // data regarding atom IDs and details of the bond work units.  These pointers
  // will point to host data, the image of the bond work units on the CPU RAM.
  bwalloc bwdims = CalcBondBlockSize(gpu, bwunits, gpu->sim.bondWorkUnits);
  int insrUintPos  = 0;
  int bondUintPos  = insrUintPos + (3 * gpu->sim.bondWorkUnits *
                                    BOND_WORK_UNIT_THREADS_PER_BLOCK);
  int anglUintPos  = bondUintPos +   (GRID * bwdims.nbondwarps);
  int diheUintPos  = anglUintPos +   (GRID * bwdims.nanglwarps);
  int cmapUintPos  = diheUintPos +   (GRID * bwdims.ndihewarps);
  int qqxcUintPos  = cmapUintPos + 2*(GRID * bwdims.ncmapwarps);
  int nb14UintPos  = qqxcUintPos +   (GRID * bwdims.nqqxcwarps);
  int nmr2UintPos  = nb14UintPos +   (GRID * bwdims.nnb14warps);
  int nmr3UintPos  = nmr2UintPos + 4*(GRID * bwdims.nnmr2warps);
  int nmr4UintPos  = nmr3UintPos + 4*(GRID * bwdims.nnmr3warps);
  int ureyUintPos  = nmr4UintPos + 4*(GRID * bwdims.nnmr4warps);
  int cimpUintPos  = ureyUintPos +   (GRID * bwdims.nureywarps);
  int cnstUintPos  = cimpUintPos +   (GRID * bwdims.ncimpwarps);
  int bondDbl2Pos  = 0;
  int anglDbl2Pos  = bondDbl2Pos +   (GRID * bwdims.nbondwarps);
  int nmr2Dbl2Pos  = anglDbl2Pos +   (GRID * bwdims.nanglwarps);
  int nmr3Dbl2Pos  = nmr2Dbl2Pos + 6*(GRID * bwdims.nnmr2warps);
  int nmr4Dbl2Pos  = nmr3Dbl2Pos + 6*(GRID * bwdims.nnmr3warps);
  int ureyDbl2Pos  = nmr4Dbl2Pos + 6*(GRID * bwdims.nnmr4warps);
  int cimpDbl2Pos  = ureyDbl2Pos +   (GRID * bwdims.nureywarps);
  int cnstDbl2Pos  = cimpDbl2Pos +   (GRID * bwdims.ncimpwarps);
  int qPmefPos     = 0;
  int dihePmefPos  = gpu->sim.bondWorkUnits * BOND_WORK_UNIT_THREADS_PER_BLOCK;
  int nb14PmefPos  = dihePmefPos +   (GRID * bwdims.ndihewarps);
  int dihePmef2Pos = 0;
  int nb14Pmef2Pos = dihePmef2Pos + 2*(GRID * bwdims.ndihewarps);
  unsigned int *pBwuInstructions = &gpu->pbBondWorkUnitUINT->_pSysData[insrUintPos];
  unsigned int *pBwuBondID  = &gpu->pbBondWorkUnitUINT->_pSysData[bondUintPos];
  unsigned int *pBwuAnglID  = &gpu->pbBondWorkUnitUINT->_pSysData[anglUintPos];
  unsigned int *pBwuDiheID  = &gpu->pbBondWorkUnitUINT->_pSysData[diheUintPos];
  unsigned int *pBwuCmapID  = &gpu->pbBondWorkUnitUINT->_pSysData[cmapUintPos];
  unsigned int *pBwuQQxcID  = &gpu->pbBondWorkUnitUINT->_pSysData[qqxcUintPos];
  unsigned int *pBwuNB14ID  = &gpu->pbBondWorkUnitUINT->_pSysData[nb14UintPos];
  unsigned int *pBwuNMR2ID  = &gpu->pbBondWorkUnitUINT->_pSysData[nmr2UintPos];
  unsigned int *pBwuNMR3ID  = &gpu->pbBondWorkUnitUINT->_pSysData[nmr3UintPos];
  unsigned int *pBwuNMR4ID  = &gpu->pbBondWorkUnitUINT->_pSysData[nmr4UintPos];
  unsigned int *pBwuUreyID  = &gpu->pbBondWorkUnitUINT->_pSysData[ureyUintPos];
  unsigned int *pBwuCImpID  = &gpu->pbBondWorkUnitUINT->_pSysData[cimpUintPos];
  unsigned int *pBwuCnstID  = &gpu->pbBondWorkUnitUINT->_pSysData[cnstUintPos];
  PMEDouble2 *pBwuBond      = &gpu->pbBondWorkUnitDBL2->_pSysData[bondDbl2Pos];
  PMEDouble2 *pBwuAngl      = &gpu->pbBondWorkUnitDBL2->_pSysData[anglDbl2Pos];
  PMEDouble2 *pBwuNMR2      = &gpu->pbBondWorkUnitDBL2->_pSysData[nmr2Dbl2Pos];
  PMEDouble2 *pBwuNMR3      = &gpu->pbBondWorkUnitDBL2->_pSysData[nmr3Dbl2Pos];
  PMEDouble2 *pBwuNMR4      = &gpu->pbBondWorkUnitDBL2->_pSysData[nmr4Dbl2Pos];
  PMEDouble2 *pBwuUrey      = &gpu->pbBondWorkUnitDBL2->_pSysData[ureyDbl2Pos];
  PMEDouble2 *pBwuCImp      = &gpu->pbBondWorkUnitDBL2->_pSysData[cimpDbl2Pos];
  PMEDouble2 *pBwuCnst      = &gpu->pbBondWorkUnitDBL2->_pSysData[cnstDbl2Pos];
  PMEFloat  *pBwuCharges    = &gpu->pbBondWorkUnitPFLOAT->_pSysData[qPmefPos];
  PMEFloat2 *pBwuLJnb14     = &gpu->pbBondWorkUnitPFLOAT2->_pSysData[nb14Pmef2Pos];
  PMEFloat  *pBwuSCEEnb14   = &gpu->pbBondWorkUnitPFLOAT->_pSysData[nb14PmefPos];
  PMEFloat  *pBwuDihe3      = &gpu->pbBondWorkUnitPFLOAT->_pSysData[dihePmefPos];
  PMEFloat2 *pBwuDihe12     = &gpu->pbBondWorkUnitPFLOAT2->_pSysData[dihePmef2Pos];

  // Make a new assignment table for the object, as well as a mask
  imat atomEndPts, objAssigned, objMasked, objMapCounts, objMap, objTerminates;
  objAssigned = CreateImat(1, nobj);
  objMasked = CreateImat(1, nobj);
  MaskBlankObjects(gpu, objcode, &objMasked);
  atomEndPts = GetAtomEndPts(gpu);
  GetObjectMaps(gpu, objtype, &atomEndPts, &objMapCounts, &objMap, &objTerminates);

  // Loop over all bond work units, striding through the unsigned integer array first.
  // The CPU will know how many bond work units there are--that's a kernel launch
  // parameter.  The pointer into the unsigned integer array must be set for each bond
  // work unit at a regular interval.  Otherwise, the unsigned integers combined with
  // object demarcations (pointers) indicate where in each of the other arrays to
  // access their details.
  unsigned int atomImports[BOND_WORK_UNIT_THREADS_PER_BLOCK];
  unsigned int instructions[BOND_WORK_UNIT_THREADS_PER_BLOCK];
  PMEFloat atomCharges[BOND_WORK_UNIT_THREADS_PER_BLOCK];
  PMEDouble2 *d2ptr;
  PMEFloat2 *f2ptr;
  PMEFloat *fptr;
  prmbuff wuParams;
  aidbuff wuAtomID;
  for (i = 0; i < gpu->sim.bondWorkUnits; i++) {
    unsigned int *uiptr;
    uiptr = &pBwuInstructions[3 * i * BOND_WORK_UNIT_THREADS_PER_BLOCK];
    for (j = 0; j < BOND_WORK_UNIT_THREADS_PER_BLOCK; j++) {

      // The first BOND_WORK_UNIT_THREADS_PER_BLOCK (currently four warps) are devoted
      // to critical instructions.  The next set is the absolute IDs of atoms to import.
      // The final set are the neighbor-list imaged IDs of atoms to import (that's what
      // is read by the GPU, but here we'll make use of absolute IDs).
      instructions[j] = uiptr[j];
      atomImports[j]  = uiptr[BOND_WORK_UNIT_THREADS_PER_BLOCK + j];
      atomCharges[j]  = pBwuCharges[j];
    }

    // Loop over instruction sets and do what they say
    for (j = 0; j < instructions[WARP_INSTRUCTION_COUNT_IDX]; j++) {
      unsigned int warpinsr = instructions[WARP_INSTRUCTION_OFFSET + j];

      // The high eight bits contain the object code
      unsigned int warpcode = warpinsr >> 24;

      // The low 24 bits contain the (initial) warp read position:
      // multiply by GRID in order to get the array index.
      unsigned int startidx = (warpinsr << 8) >> (8 - GRID_BITS);

      // Continue if the warp code does not match the object under scrutiny
      // (inefficient, but this is for debugging purposes).
      if (warpcode != objcode) {
        continue;
      }

      // Investigate each object (of the indicated type) that the work unit deals with.
      // General procedure: get the raw atom ID information (usually one unsigned int,
      // two in the case of CMAPs), unpack it to find the imported atom IDs to reference,
      // look up what atoms those came from in the master topology, then cross-reference
      // the atoms against lists of all as-yet unassigned objects those atoms are found
      // in to make a list of candidates for what object this could be.  Compare the
      // details of the candidate objects against the ones in the bond work unit and
      // check off the first match.
      for (k = 0; k < GRID; k++) {
        unsigned int rawID;
        if (warpcode == BOND_CODE) {
          rawID = pBwuBondID[startidx + k];
          d2ptr = pBwuBond;
        }
        else if (warpcode == ANGL_CODE) {
          rawID = pBwuAnglID[startidx + k];
          d2ptr = pBwuAngl;
        }
        else if (warpcode == DIHE_CODE) {
          rawID = pBwuDiheID[startidx + k];
          f2ptr = pBwuDihe12;
          fptr = pBwuDihe3;
        }
        else if (warpcode == CMAP_CODE) {
          rawID = pBwuCmapID[2*startidx + k];
        }
        else if (warpcode == QQXC_CODE) {
          rawID = pBwuQQxcID[startidx + k];
        }
        else if (warpcode == NB14_CODE) {
          rawID = pBwuNB14ID[startidx + k];
          f2ptr = pBwuLJnb14;
          fptr  = pBwuSCEEnb14;
        }
        else if (warpcode == NMR2_CODE) {
          rawID = pBwuNMR2ID[4*startidx + k];
          d2ptr = pBwuNMR2;
        }
        else if (warpcode == NMR3_CODE) {
          rawID = pBwuNMR3ID[4*startidx + k];
          d2ptr = pBwuNMR3;
        }
        else if (warpcode == NMR4_CODE) {
          rawID = pBwuNMR4ID[4*startidx + k];
          d2ptr = pBwuNMR4;
        }
        else if (warpcode == UREY_CODE) {
          rawID = pBwuUreyID[startidx + k];
          d2ptr = pBwuUrey;
        }
        else if (warpcode == CIMP_CODE) {
          rawID = pBwuCImpID[startidx + k];
          d2ptr = pBwuCImp;
        }
        else if (warpcode == CNST_CODE) {
          rawID = pBwuCnstID[startidx + k];
          d2ptr = pBwuCnst;
        }

        // Values of 0xffffffff in the imported atom IDs package
        // will indicate that there is nothing to compute.
        if (rawID == 0xffffffff) {
          continue;
        }
        if (warpcode == BOND_CODE || warpcode == QQXC_CODE || warpcode == NB14_CODE ||
            warpcode == NMR2_CODE || warpcode == UREY_CODE) {
          wuAtomID.refI = atomImports[(rawID >> 8) & 0xff];
          wuAtomID.refJ = atomImports[rawID & 0xff];
        }
        if (warpcode == ANGL_CODE || warpcode == NMR3_CODE) {
          wuAtomID.refI = atomImports[(rawID >> 16) & 0xff];
          wuAtomID.refJ = atomImports[(rawID >>  8) & 0xff];
          wuAtomID.refK = atomImports[rawID & 0xff];
        }
        if (warpcode == DIHE_CODE || warpcode == NMR4_CODE || warpcode == CIMP_CODE ||
            warpcode == CMAP_CODE) {
          wuAtomID.refI = atomImports[(rawID >> 24) & 0xff];
          wuAtomID.refJ = atomImports[(rawID >> 16) & 0xff];
          wuAtomID.refK = atomImports[(rawID >>  8) & 0xff];
          wuAtomID.refL = atomImports[rawID & 0xff];
        }
        if (warpcode == CNST_CODE) {
          wuAtomID.refI = atomImports[rawID & 0xff];
        }
        if (warpcode == CMAP_CODE) {
          wuParams.i0 = pBwuCmapID[2*startidx + GRID + k] >> 8;
          wuAtomID.refM = atomImports[pBwuCmapID[2*startidx + GRID + k] & 0xff];
        }
        if (warpcode == BOND_CODE || warpcode == ANGL_CODE || warpcode == NMR2_CODE ||
            warpcode == NMR3_CODE || warpcode == NMR4_CODE || warpcode == UREY_CODE ||
            warpcode == CIMP_CODE) {
          wuParams.p1 = d2ptr[startidx + k].x;
          wuParams.p2 = d2ptr[startidx + k].y;
        }
        if (warpcode == NMR2_CODE || warpcode == NMR3_CODE || warpcode == NMR4_CODE) {
          wuParams.p3  = d2ptr[6*startidx +   GRID + k].x;
          wuParams.p4  = d2ptr[6*startidx +   GRID + k].y;
          wuParams.p5  = d2ptr[6*startidx + 2*GRID + k].x;
          wuParams.p6  = d2ptr[6*startidx + 2*GRID + k].y;
          wuParams.p7  = d2ptr[6*startidx + 3*GRID + k].x;
          wuParams.p8  = d2ptr[6*startidx + 3*GRID + k].y;
          wuParams.p9  = d2ptr[6*startidx + 4*GRID + k].x;
          wuParams.p10 = d2ptr[6*startidx + 4*GRID + k].y;
          wuParams.p11 = d2ptr[6*startidx + 5*GRID + k].x;
          wuParams.p12 = d2ptr[6*startidx + 5*GRID + k].y;
        }
        if (warpcode == DIHE_CODE) {
          wuParams.p1 = f2ptr[2*startidx +        k].x;
          wuParams.p2 = f2ptr[2*startidx +        k].y;
          wuParams.p3 = f2ptr[2*startidx + GRID + k].x;
          wuParams.p4 = f2ptr[2*startidx + GRID + k].y;
          wuParams.p5 = fptr[startidx + k];
        }
        if (warpcode == NB14_CODE) {
          wuParams.p1 = f2ptr[startidx + k].x;
          wuParams.p2 = f2ptr[startidx + k].y;
          wuParams.p3 = fptr[startidx + k];
        }
        if (warpcode == CNST_CODE) {
          wuParams.p1 = d2ptr[2*startidx + k].x;
          wuParams.p2 = d2ptr[2*startidx + k].y;
          wuParams.p3 = d2ptr[2*startidx + GRID + k].x;
          wuParams.p4 = d2ptr[2*startidx + GRID + k].y;
        }
        if (MatchObjectByWorkUnitData(gpu, objcode, wuAtomID, wuParams, &objAssigned,
                                      p1atm, p2atm, p4atm) == 0) {
          printf("CheckWorkUnitTranslations :: Work Unit %4d %s warp %2d-%2d matches no known "
                 "terms.\n", i, objtype, j, k);
        }
      }
    }
  }

  // Check all assignments
  bool problem = false;
  for (i = 0; i < nobj; i++) {
    if (objAssigned.data[i] == 0 && objMasked.data[i] == 0) {
      printf("CheckWorkUnitTranslations :: %s %6d is not covered.\n", objtype, i);
      problem = true;
    }
    if (objAssigned.data[i] > 1) {
      printf("CheckWorkUnitTranslations :: %s %6d is covered by %d work units.\n", objtype,
             i, objAssigned.data[i]);
      problem = true;
    }
  }
  if (problem == false) {
    printf("CheckWorkUnitTranslations :: All %s have been correctly assigned.\n", objtype);
  }

  // Free allocated memory
  DestroyImat(&objMapCounts);
  DestroyImat(&objMap);
  DestroyImat(&objAssigned);
  DestroyImat(&objMasked);
  DestroyImat(&objTerminates);
  DestroyImat(&atomEndPts);
}
