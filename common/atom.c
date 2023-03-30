/****************************************************************************/
/* This file contains the transformations for the atom structure */

#include "nc.h"
#include "atom.h"
#include "vector.h"

/****************************************************************************/
// Replaces the first nUniqueAtoms of atoms with just the unique atoms

int removeDuplicatedAtoms(atom *atoms, int nAtoms, int nMaxBonds) {
	int i, j, uniqueAtomFlag, nUniqueAtoms = 0;
	vector diff;
	atom *uniqueAtoms;

	// Dynamically allocate memory
	uniqueAtoms = (atom *) calloc(nAtoms, sizeof(atom));
	for (i = 0; i < nAtoms; i++) {
		uniqueAtoms[i].neighborList = (int *) calloc(nMaxBonds, sizeof(int));
		uniqueAtoms[i].neighborPos = (vector *) calloc(nMaxBonds, sizeof(vector));
	}

	// Find unique atoms (if two atoms within EPS then they are the same atom)
	for (i = 0; i < nAtoms; i++) {
		uniqueAtomFlag = 1;
		for (j = i+1; j < nAtoms; j++) {
			diff = retSubtractedVectors(atoms[i].pos, atoms[j].pos);
			if (diff.mag < EPS) {
				uniqueAtomFlag = 0;
				break;
			}
		}
		if (uniqueAtomFlag) {
			deepCopySingleAtom(uniqueAtoms, nUniqueAtoms, atoms, i, nMaxBonds);
			nUniqueAtoms++;
		}
	}

	// Deep copy finalAtoms into the first nUniqueAtoms spots of atoms
	calcNearestNeighborLists(uniqueAtoms, nMaxBonds, nUniqueAtoms);
	deepCopyAllAtoms(atoms, uniqueAtoms, nUniqueAtoms, nMaxBonds);

	// Free dynamically allocated memory
	for (i = 0; i < nAtoms; i++) {
		free(uniqueAtoms[i].neighborPos);
		free(uniqueAtoms[i].neighborList);
	}
	free(uniqueAtoms);

	return nUniqueAtoms;
}

/****************************************************************************/
// Translates the atom position of all nAtoms in atoms by shiftVector
// shiftVector will be the new origin of the atoms

void translateAllAtoms(atom *atoms, int nAtoms, vector shiftVector) {
	int i;

	for (i = 0; i < nAtoms; i++) {
		atoms[i].pos = retSubtractedVectors(atoms[i].pos, shiftVector);
	}

	return;
}

/****************************************************************************/
// Shallow copy all the origAtoms pointers to newAtoms pointers
// newAtoms will now point to the same place as origAtoms

int shallowCopyAllAtoms(atom *newAtoms, atom *origAtoms, int numAtomsToCopy) {
	int i, numAtomsCopied = 0;

	for (i = 0; i < numAtomsToCopy; i++) {
		newAtoms[i] = origAtoms[i];
		numAtomsCopied++;
	}

	return numAtomsCopied;
}

/****************************************************************************/
// Deep copies all the origAtoms information to newAtoms
// newAtoms must have separately allocated memory for neighborPos/List
// newAtoms and origAtoms point to different locations now

int deepCopyAllAtoms(atom *newAtoms, atom *origAtoms, int numAtomsToCopy, int nMaxBonds) {
	int i, numAtomsCopied = 0;

	for (i = 0; i < numAtomsToCopy; i++) {
		deepCopySingleAtom(newAtoms, i, origAtoms, i, nMaxBonds);
		numAtomsCopied++;
	}

	return numAtomsCopied;
}

/****************************************************************************/
// Deep copy origAtoms[origIndex] to newAtoms[newIndex] 

void deepCopySingleAtom(atom *newAtoms, int newIndex, atom *origAtoms, int origIndex, int nMaxBonds) {
	int i;

	newAtoms[newIndex].type = origAtoms[origIndex].type;
	strcpy(newAtoms[newIndex].symbol, origAtoms[origIndex].symbol);
	newAtoms[newIndex].mass = origAtoms[origIndex].mass;
	newAtoms[newIndex].pos = origAtoms[origIndex].pos;
	for (i = 0; i < nMaxBonds; i++) {
		newAtoms[newIndex].neighborPos[i].x = origAtoms[origIndex].neighborPos[i].x;
		newAtoms[newIndex].neighborPos[i].y = origAtoms[origIndex].neighborPos[i].y;
		newAtoms[newIndex].neighborPos[i].z = origAtoms[origIndex].neighborPos[i].z;
		newAtoms[newIndex].neighborPos[i].mag = origAtoms[origIndex].neighborPos[i].mag;
		newAtoms[newIndex].neighborList[i] = origAtoms[origIndex].neighborList[i];
	}

	return;
}

/****************************************************************************/
// Deep zeroing/ nullifying of all atoms 

void deepZeroAllAtoms(atom *atoms, int numAtomsToZero, int nMaxBonds) {
	int i;

	for (i = 0; i < numAtomsToZero; i++) deepZeroSingleAtom(atoms, i, nMaxBonds);

	return;
}

/****************************************************************************/
// Deep zeroing/ nullifying of a single atom 

void deepZeroSingleAtom(atom *atoms, int index, int nMaxBonds) {
	int i;

	atoms[index].type = 0;
	strcpy(atoms[index].symbol, '\0');
	atoms[index].mass = 0.0;
	atoms[index].pos = retZeroVector();
	for (i = 0; i < nMaxBonds; i++) {
		atoms[index].neighborPos[i] = retZeroVector();
		atoms[index].neighborList[i] = -1; // never an atom index
	}

	return;
}

/****************************************************************************/

// int callocEntireAtomMemory(atom *atoms, int nAtoms, int nMaxBonds) {
// 	int i;

// 	// Dynamically allocate memory for the overall atom structure
// 	atoms = (atom *) calloc(nAtoms, sizeof(atom));
// 	if (! atoms && nAtoms) {
// 		memoryError("Error in allocating memory for an atom structure\n");
// 	}

// 	// Dynamically allocate memory for the neighbor lists for each atom structure
// 	for (i = 0; i < nAtoms; i++) {
// 		atoms[i].neighborList = (int *) calloc(nMaxBonds, sizeof(int));
// 		atoms[i].neighborPos = (vector *) calloc(nMaxBonds, sizeof(vector));
// 		if (! atoms[i].neighborPos && nMaxBonds) {
// 			memoryError("Error in allocating memory for an atom neighbor list\n");
// 		}
// 		if (! atoms[i].neighborList && nMaxBonds) {
// 			memoryError("Error in allocating memory for an atom neighbor list\n");
// 		}
// 	}


// 	return 0; // successful allocation
// }

// /****************************************************************************/

// void freeEntireAtomMemory(atom *atoms, int nAtoms) {
// 	int i;

// 	// Free dynamically allocated memory
// 	for (i = 0; i < nAtoms; i++) {
// 		free(atoms[i].neighborPos);
// 		free(atoms[i].neighborList);
// 	}
// 	free(atoms);

// 	return; 
// }

/****************************************************************************/
