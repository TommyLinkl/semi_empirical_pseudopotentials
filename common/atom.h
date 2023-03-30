/****************************************************************************/
//
//
//
/****************************************************************************/

#ifndef ATOM_H
#define ATOM_H

/****************************************************************************/
/* These are the library functions that are used within atoms.c */

#include <string.h>
#include "vector.h"

/****************************************************************************/
/* Structure declarations */

typedef struct atom_ {
  int type, index;
  char symbol[100];
  double mass;
  vector pos;
  vector *neighborPos;
  int *neighborList;
} atom;

/****************************************************************************/
/* Macro definitions: common unit conversion and multiplication schemes 
  that help make the program more readable */




/****************************************************************************/
/* Function declarations - public interface */

/* Functions for atom structures - in atom.c */
int removeDuplicatedAtoms(atom *atoms, int nAtoms, int nMaxBonds);
void translateAllAtoms(atom *atoms, int nAtoms, vector shiftVector);
int shallowCopyAllAtoms(atom *newAtoms, atom *origAtoms, int numAtomsToCopy);
int deepCopyAllAtoms(atom *newAtoms, atom *origAtoms, int numAtomsToCopy, int nMaxBonds);
void deepCopySingleAtom(atom *newAtoms, int newIndex, atom *origAtoms, int origIndex, int nMaxBonds);
void deepZeroAllAtoms(atom *atoms, int numAtomsToZero, int nMaxBonds);
void deepZeroSingleAtom(atom *atoms, int index, int nMaxBonds);
//int callocEntireAtomMemory(atom *atoms, int nAtoms, int nMaxBonds);
//void freeEntireAtomMemory(atom *atoms, int nAtoms);
//int callocAtomMemoryNoNeighbors(atom *atoms, int nAtoms);
//int callocAtomNeighborMemory(aomt *aotms, int nAtoms, int nMaxBonds);

/****************************************************************************/

#endif

/****************************************************************************/