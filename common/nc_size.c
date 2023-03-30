/****************************************************************************/
/* This file calculates statistics realted to the size of a nanocrystal */

#include "nc.h"

/***************************************************************************/

void calcSizeStatistics(atom *atoms, param params) {
  FILE *pf;
  int i, j, numSelectAtoms = 0;
  char currentSymbol[100];
  double volume, area, diameter;
  vector *atomPositions, *maxDimensions;

  // Assign atom types
  params.nAtomTypes = assignAtomTypes(atoms, params.nAtoms);

  // Open the file that will contain the output/ resulting statistics
  pf = fopen("nc_size.dat", "w");
  writeCurrentTime(pf);
  writeSeparation(pf);
  writeSystemInfo(atoms, params, pf);
  writeSeparation(pf);
  fprintf(pf, "Nanocrystal Size Statistics [Angstroms]\n\n");

  // Allocate memory for the positions of all atoms 
  atomPositions = (vector *) calloc(params.nAtoms, sizeof(vector));
  maxDimensions = (vector *) calloc(params.nAtomTypes+2, sizeof(vector));

  // Copy all atom positions into vector only atomPositions and 
  // calculate the size statistics for the entire nanocrystal
  for (i = 0; i < params.nAtoms; i++) atomPositions[i] = atoms[i].pos;
  calcNanocrystalDimensions(&maxDimensions[0], atomPositions, params.nAtoms);
  
  // Do not include passivation ligands in dimensions counting
  for (i = 0; i < params.nAtoms; i++) {
    if (! (atoms[i].symbol[0] == 'P')) {
      atomPositions[i] = atoms[i].pos;
      numSelectAtoms += 1;
    }
  }
  calcNanocrystalDimensions(&maxDimensions[1], atomPositions, numSelectAtoms);

  // Print (in nm) volume, area, radius, length, width, thickness to nc_size.dat
  if (! strcmp(params.ncType, "attachedNCs")) strcpy(params.ncType, "NPL"); 
  volume = retNanocrystalVolume(atomPositions, numSelectAtoms, params.ncType);
  area = retNanocrystalArea(atomPositions, numSelectAtoms, params.ncType);
  fprintf(pf, "The volume of the %s = %.3f nm^3\n", params.ncType, 0.001*volume);
  fprintf(pf, "The area of the %s = %.3f nm^2\n", params.ncType, 0.01*area);
  if (! strcmp(params.ncType, "QD")) {
    fprintf(pf, "The radius of the QD = %.3f nm\n", 0.05*maxDimensions[1].mag);
    fprintf(pf, "The diameter of the QD = %.3f nm\n\n", 0.1*maxDimensions[1].mag);
  } else if (! strcmp(params.ncType, "NR")) {
    diameter = retMaxPlaneProjectedLength(atomPositions, numSelectAtoms, "xy");
    fprintf(pf, "The radius of the NR = %.3f nm\n", 0.05*diameter);
    fprintf(pf, "The diameter of the NR = %.3f nm\n", 0.1*diameter);
    fprintf(pf, "The length of the NR = %.3f nm\n\n", 0.1*maxDimensions[1].z);
  } else if (! strcmp(params.ncType, "NPL")) {
    fprintf(pf, "The length of the NPL = %.3f nm\n", 0.1*maxDimensions[1].x);
    fprintf(pf, "The width of the NPL = %.3f nm\n", 0.1*maxDimensions[1].y);
    fprintf(pf, "The thickness of the NPL = %.3f nm\n\n", 0.1*maxDimensions[1].z);
  }

  // Copy only atoms of one type and calculate the size statistics 
  for (j = 0; j < params.nAtomTypes; j++) {
    numSelectAtoms = 0;
    for (i = 0; i < params.nAtoms; i++) {
      if (atoms[i].type == j) {
        strcpy(currentSymbol, atoms[i].symbol); 
        atomPositions[numSelectAtoms] = atoms[i].pos;
        numSelectAtoms += 1;
      }
    }
    calcNanocrystalDimensions(&maxDimensions[j+2], atomPositions, numSelectAtoms);
  }

  // Print end of file separator and close nc_size.dat
  writeSizeResults(maxDimensions, atoms, params, pf);
  writeSeparation(pf);
  fclose(pf);

  // Print out a lammps configuration file
  pf = fopen("lammpsconf.par", "w");
  if (! params.passivate) {
    writeLammpsConfFile(atoms, maxDimensions[1], params.nAtoms, params.nAtomTypes, pf);
  }
  fclose(pf);

  // Free dynamically allocated memory
  free(atomPositions); free(maxDimensions);

  return;
}

/***************************************************************************/

void calcNanocrystalDimensions(vector *maxDimensions, vector *atomPositions, int numVectors) {
  
  maxDimensions[0].x = retMaxLength(atomPositions, numVectors, 'x');
  maxDimensions[0].y = retMaxLength(atomPositions, numVectors, 'y');
  maxDimensions[0].z = retMaxLength(atomPositions, numVectors, 'z');
  maxDimensions[0].mag = retMaxDistance(atomPositions, numVectors);

  return;
}

/***************************************************************************/
// Returns 4/3*pi*r^3 for QD, pi*r^2*length(z) for NR 
// and length(x)*width(y)*thickness(z) for NPLs

double retNanocrystalVolume(vector *atomPositions, int numVectors, char *ncType) {

  if (! strcmp(ncType, "QD")) 
    return 4.0/3.0*PIE*cube(0.5*retMaxDistance(atomPositions, numVectors));
  else if (! strcmp(ncType, "NR"))
    return PIE*sqr(0.5*retMaxPlaneProjectedLength(atomPositions, numVectors, "xy"))
              *retMaxLength(atomPositions, numVectors, 'z');
  else if (! strcmp(ncType, "NPL"))
    return retMaxLength(atomPositions, numVectors, 'x')*retMaxLength(atomPositions, numVectors, 'y')
              *retMaxLength(atomPositions, numVectors, 'z');

  return 0.0;
}

/***************************************************************************/
// Returns pi*r^2 for QD and NR and length(x)*width(y) for NPL

double retNanocrystalArea(vector *atomPositions, int numVectors, char *ncType) {

  if (! strcmp(ncType, "QD")) 
    return PIE*sqr(0.5*retMaxDistance(atomPositions, numVectors));
  else if (! strcmp(ncType, "NR"))
    return PIE*sqr(0.5*retMaxPlaneProjectedLength(atomPositions, numVectors, "xy"));
  else if (! strcmp(ncType, "NPL"))
    return retMaxLength(atomPositions, numVectors, 'x')*retMaxLength(atomPositions, numVectors, 'y');

  return 0.0;
}

/***************************************************************************/
// Returns max length between two vectors that have been projected onto
// a plane specified by plane -> allowed are 'xy', 'xz' and 'yz'

double retMaxPlaneProjectedLength(vector *vectors, int numVectors, char *plane) {
  int i;
  double maxProjectedLength = 0.0;
  vector *projectedVectors;

  // Allocate memory for the positions of all atoms 
  projectedVectors = (vector *) calloc(numVectors, sizeof(vector));

  // Copy vectors into projectedVectors
  deepCopyVectors(projectedVectors, vectors, numVectors);

  // zero the direction that is orthogonal to the desired plane
  if (! strcmp(plane, "xy")) for (i = 0; i < numVectors; i++) projectedVectors[i].z = 0.0;
  else if (! strcmp(plane, "xz")) for (i = 0; i < numVectors; i++) projectedVectors[i].y = 0.0;
  else if (! strcmp(plane, "yz")) for (i = 0; i < numVectors; i++) projectedVectors[i].x = 0.0;

  // Can use retMaxDistance now to the get max plane projected length
  maxProjectedLength = retMaxDistance(projectedVectors, numVectors);

  // Free memory allocated within this function
  free(projectedVectors);

  return maxProjectedLength;
}

/***************************************************************************/

double retMaxLength(vector *vectors, int numVectors, char axis) {
  long i, j;
  double testLen2, maxLen2 = 0.0;

  for (i = 0; i < numVectors-1; i++) {
    for (j = i+1; j < numVectors; j++) {
      if (axis == 'x') testLen2 = sqr(vectors[i].x - vectors[j].x);
      else if (axis == 'y') testLen2 = sqr(vectors[i].y - vectors[j].y);
      else if (axis == 'z') testLen2 = sqr(vectors[i].z - vectors[j].z);
      else {
        printf("Must provide a char of 'x', 'y', or 'z' when calling retMaxLength\n");
        exit(EXIT_FAILURE);
      }
      if (testLen2 > maxLen2) maxLen2 = testLen2;
    }
  }
  
  return (sqrt(maxLen2));
}

/***************************************************************************/

double retMaxDistance(vector *vectors, int numVectors) {
  int i, j;
  double testDistance, maxDistance = 0.0;
  
  for (i = 0; i < numVectors-1; i++) {
    for (j = i+1; j < numVectors; j++) {
      testDistance = retDistanceBetweenPoints(vectors[i], vectors[j]);
      if (testDistance > maxDistance) maxDistance = testDistance;
    }
  }
  
  return maxDistance;
}

/***************************************************************************/
