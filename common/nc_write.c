/****************************************************************************/
/* This file does the printing for the program */

#include "nc.h"

/****************************************************************************/
//  Writes a configuration file that can be used in an electronic structure
//  calculation -> should utilize unitChange to print in desired units 
//  atomic units are required for the electronic structure calculations

void writeConfFile(atom *atoms, int nAtoms, double unitChange, FILE *pf) {
  int i;

  fprintf(pf, "%d\n", nAtoms);
  for (i = 0; i < nAtoms; i++) fprintf(pf, "%s % .8f % .8f % .8f\n", 
    atoms[i].symbol, atoms[i].pos.x*unitChange, 
    atoms[i].pos.y*unitChange, atoms[i].pos.z*unitChange);

  return;
}

/****************************************************************************/
// Writes a vmd, vesta, ... readable file -> prints H for passivation ligands 
// TODO: pass in character array instead of FILE * and add in .xyz 
// TODO: pass in character array for an optional comment

void writeXYZFile(atom *atoms, int nAtoms, double unitChange, FILE *pf) {
	int i;

	fprintf(pf, "%d\n", nAtoms);
	fprintf(pf, "# Comment line\n");
	for (i = 0; i < nAtoms; i++) { 
    if (isAPassivationSymbol(atoms[i].symbol)) fprintf(pf, "%s ", "H");
    else fprintf(pf, "%s ", atoms[i].symbol);
    fprintf(pf, "% .8f % .8f % .8f\n",
			atoms[i].pos.x*unitChange, atoms[i].pos.y*unitChange, atoms[i].pos.z*unitChange);
  }

	return;
}

/****************************************************************************/
// Writes a lammps configuration file. The atom positions must be in Angstroms

void writeLammpsConfFile(atom *atoms, vector ncDimensions, int nAtoms, int nAtomTypes, FILE *pf) {
  int i, j;

  fprintf(pf, "LAMMPS configuration data file\n\n");
  fprintf(pf, "  %d atoms\n\n", nAtoms);
  fprintf(pf, "  %d atom types\n\n", nAtomTypes);
  fprintf(pf, "  % .6f % .6f  xlo xhi\n", -0.5*ncDimensions.x-20.0, 0.5*ncDimensions.x+20.0);
  fprintf(pf, "  % .6f % .6f  ylo yhi\n", -0.5*ncDimensions.y-20.0, 0.5*ncDimensions.y+20.0);
  fprintf(pf, "  % .6f % .6f  zlo zhi\n\n", -0.5*ncDimensions.z-20.0, 0.5*ncDimensions.z+20.0);
  fprintf(pf, " Masses\n\n");
  for (i = 0; i < nAtomTypes; i++) {
    for (j = 0; j < nAtoms; j++) {
      if (atoms[j].type == i) {
        fprintf(pf, "  %d  % .3f\n", i+1, retAtomMass(atoms[j].symbol));
        break;
      }
    }
  }
  fprintf(pf, "\n Atoms\n\n");
  for (i = 0; i < nAtoms; i++) { 
    fprintf(pf, "  %7d  %2d  ", i+1, atoms[i].type+1);
    fprintf(pf, " % 5.8f % 5.8f % 5.8f \n", atoms[i].pos.x, atoms[i].pos.y, atoms[i].pos.z);
  }

  return;
}

/****************************************************************************/
// 

void writeNearestNeighbors(atom *atoms, param params, FILE *pf) {
	int i, k, numBonds;

	writeSeparation(pf);
	// Loop over each atom and its neighbors
	for (i = 0; i < params.nAtoms; i++) {
		numBonds = 0;
		for (k = 0; k < params.nMaxBonds;  k++) if (atoms[i].neighborPos[k].mag > EPS) numBonds++;
		fprintf(pf, "%s atom %d with the number of bonds = %d\n", atoms[i].symbol, i, numBonds);
		for (k = 0; k < numBonds;  k++) fprintf(pf, "Bond %d with neighbor %s has length = %.4f\n",
			k, atoms[atoms[i].neighborList[k]].symbol, atoms[i].neighborPos[k].mag);
		fprintf(pf, "\n");
	}
	writeSeparation(pf);

	return;
}

/****************************************************************************/
//

void writeSizeResults(vector *maxDimensions, atom *atoms, param params, FILE *pf) {
  int i, j, selectAtomIndex, numSelectAtoms = 0;

  // Approximate number of grid points to be used in our real space 
  // electronic structure calculations if crystal is passivated
  vector gridPointVector = retScaledVector(maxDimensions[0], ANGTOAU);
  gridPointVector.x += 10.0; gridPointVector.y += 10.0; gridPointVector.z += 10.0;
  gridPointVector = retScaledVector(gridPointVector, 1.25); // 1.25 = 1.0/0.8 = desired grid density
  gridPointVector.mag = gridPointVector.x*gridPointVector.y*gridPointVector.z;
  fprintf(pf, "The approximate number of grid points (grid density 0.8):\n");
  writeVectorShort(gridPointVector, pf);

  for (j = 0; j < params.nAtomTypes+2; j++) {
    if (! j) fprintf(pf, "\nAll atom types - including passivation\n");
    else if (j == 1) fprintf(pf, "All atom types - excluding passivation\n");
    else {
      numSelectAtoms = 0;
      for (i = 0; i < params.nAtoms; i++) if (atoms[i].type == j-2) {
        numSelectAtoms += 1;
        selectAtomIndex = i;
      }
      fprintf(pf, "Atom type = %d, atom symbol = %s, number of atoms = %d\n", 
        j-2, atoms[selectAtomIndex].symbol, numSelectAtoms);
    }
    writeVectorShort(maxDimensions[j], pf);
    fprintf(pf, "\n");
  }

  return;
}

/****************************************************************************/

void writeSystemInfo(atom *atoms, param params, FILE *pf) {
  int i, j, numSelectAtoms = 0;

  fprintf(pf, "System Information\n\n");
  fprintf(pf, "Total Number of Atoms = %d\n", params.nAtoms);
  fprintf(pf, "Number of atom types = %d\n", params.nAtomTypes);
  for (j = 0; j < params.nAtomTypes; j++) {
    for (i = 0; i < params.nAtoms; i++) if (atoms[i].type == j) {
      numSelectAtoms += 1;
      fprintf(pf, "Atom type %d has symbol %s\n", j, atoms[i].symbol);
      break;
    }
  }
    
  return;
}

/****************************************************************************/

void writeIntegerArray(int *intArray, int arrayLength, char *fileName) {
  FILE *pf;
  int i;

  pf = fopen(fileName, "w");
  for (i = 0; i < arrayLength; i++) fprintf(pf, "% d % d\n", i, intArray[i]);
  fclose(pf);

  return;
}

/****************************************************************************/

void writeVector(vector vect, FILE *pf) {
  fprintf(pf, "%14.10f %14.10f %14.10f %14.10f\n", vect.x, vect.y, vect.z, vect.mag);

  return;
}

/****************************************************************************/

void writeVectorShort(vector vect, FILE *pf) {
  fprintf(pf, "%10.3f %10.3f %10.3f %10.3f\n", vect.x, vect.y, vect.z, vect.mag);

  return;
}

/****************************************************************************/
// prints out current time to stdout 
 
void writeCurrentTime(FILE *pf) {
  time_t startTime;

  startTime = time(NULL);
  fprintf(pf, ctime(&startTime));

  return;
}

/****************************************************************************/

void writeSeparation(FILE *pf) {
  fprintf(pf, "\n******************************************************************************\n\n");
  
  return;
}

/****************************************************************************/
