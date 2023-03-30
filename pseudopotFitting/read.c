/*****************************************************************************************/

#include <string.h>
#include "fit.h"

/*****************************************************************************************/

void readInput(param *params) {
  FILE *pf;
  int i;
  char field[1000], tmp[1000];

  // Set default / initial values
  params->nAtomTypes = 0;
  params->kineticEnergyScaling = 1.00;

  pf = fopen("input.par", "r");
  while (1) {
    fscanf(pf, "%s", field);
    if (! strcmp(field, "maxKE")) fscanf(pf, "%s %lg", tmp, &params->maxKE);
    else if (! strcmp(field, "nBandStructures")) fscanf(pf, "%s %d", tmp, &params->nBandStructures);
    else if (! strcmp(field, "nAtoms")) fscanf(pf, "%s %d", tmp, &params->nAtoms);
    else if (! strcmp(field, "nElectrons")) fscanf(pf, "%s %d", tmp, &params->nElectrons);
    else if (! strcmp(field, "nIterations")) fscanf(pf, "%s %d", tmp, &params->nIterations);
    else if (! strcmp(field, "nKPoints")) fscanf(pf, "%s %d", tmp, &params->nKPoints);
    else if (! strcmp(field, "beta")) fscanf(pf, "%s %lg", tmp, &params->beta);
    else if (! strcmp(field, "stepSize")) fscanf(pf, "%s %lg", tmp, &params->stepSize);
    else if (! strcmp(field, "kineticEnergyScaling")) fscanf(pf, "%s %lg", tmp, &params->kineticEnergyScaling);
    else {
       printf("Invalid input field and/ or format of input\n");
       exit(EXIT_FAILURE);
    }
    i++;
    if (i > 9) break;
  }
  fclose(pf);

  //
  if (params->nBandStructures > 1) readMultiInput(params);
  else copyBaseParamsToArrays(params);
 
  return;
}

/*****************************************************************************************/

void copyBaseParamsToArrays(param *params) {
  params->na[0] = params->nAtoms;
  params->ne[0] = params->nElectrons;
  params->nkp[0] = params->nKPoints;
  params->mke[0] = params->maxKE;
  params->nKPE = params->nkp[0]*params->ne[0];

  return;
}

/*****************************************************************************************/

void readMultiInput(param *params) {
  FILE *pf;
  int i;
  char field[1000], tmp[1000];

  params->nKPE = 0;

  pf = fopen("multiInput.par", "r");
  for (i = 0; i < params->nBandStructures; i++) {
    fscanf(pf, "%s %s %d", field, tmp, &params->na[i]);
    fscanf(pf, "%s %s %d", field, tmp, &params->ne[i]);
    fscanf(pf, "%s %s %d", field, tmp, &params->nkp[i]);
    fscanf(pf, "%s %s %lg", field, tmp, &params->mke[i]);
    params->nKPE += params->nkp[i]*params->ne[i];
  }
  fclose(pf);

  return;
}

/*****************************************************************************************/

void readExpBandStructure(double *expBandStructure, int nKPoints, int nElectrons, int fileIndex) {
  FILE *pf;
  int i, j;
  double tmp;
  char fileName[100], tmpChar[2];  
  
  strcpy(fileName, "expBandStruct_");
  sprintf(tmpChar, "%d.par", fileIndex);
  strcat(fileName, tmpChar);

  pf = fopen(fileName, "r");
  for (i = 0; i < nKPoints; i++) {
    fscanf(pf, "%g", &tmp);
    for (j = 0; j < nElectrons; j++) fscanf(pf, "%lg", &expBandStructure[i*nElectrons+j]); 
  }
  fclose(pf);  
  
  return;
}

/*****************************************************************************************/

void readKPoints(vector *kPoint, int nKPoints, double scale, int fileIndex) {
  FILE *pf;
  int i;
  char fileName[100], tmpChar[2];  
  
  strcpy(fileName, "kpoints_");
  sprintf(tmpChar, "%d.par", fileIndex);
  strcat(fileName, tmpChar);

  pf = fopen(fileName, "r"); 
  for (i = 0; i < nKPoints; i++) {
    fscanf(pf, "%lg %lg %lg", &kPoint[i].x, &kPoint[i].y, &kPoint[i].z);
    kPoint[i] = retScaledVector(kPoint[i], (TWOPI / scale));
  }
  fclose(pf);


  
  return;
}

/*****************************************************************************************/

void readSystem(param *params, vector *unitCellvectors, atom *atoms, int fileIndex) {
  FILE *pf;
  int i;
  char buf[1000]; 
  char tmp[1000];
  char fileName[100], tmpChar[2];  
  
  strcpy(fileName, "system_");
  sprintf(tmpChar, "%d.par", fileIndex);
  strcat(fileName, tmpChar);
 
  pf = fopen(fileName, "r");
  while (fscanf(pf, "%s", buf) == 1) {
    if (! strcmp(buf, "atoms")) {
      for (i = 0; i < params->na[fileIndex]; i++) {
        fscanf (pf, "%s %lg %lg %lg", atoms[i].symbol, &atoms[i].pos.x, &atoms[i].pos.y, &atoms[i].pos.z);
      if (newAtomType(atoms, i)) params->nAtomTypes++;
      }
    }
    else if (! strcmp(buf, "cell")) {
      for (i = 0; i < NDIM; i++) {
        fscanf(pf, "%lg %lg %lg", &unitCellvectors[i].x, &unitCellvectors[i].y, &unitCellvectors[i].z);
      }
    }
    else if (! strcmp(buf, "scale")) fscanf(pf, "%s %lg", tmp, &(params->s[fileIndex]));
  }
  fclose(pf);
 
  // Scale the unit cell vectors by the input parameter scale
  for (i = 0; i < NDIM; i++) {
    unitCellvectors[i] = retScaledVector(unitCellvectors[i], params->s[fileIndex]);
  }
 
  // Scale the atom positions by the unit cell vectors 
  for (i = 0; i < params->na[fileIndex]; i++) {
    atoms[i].pos.x *= (unitCellvectors[0].x + unitCellvectors[1].x + unitCellvectors[2].x);
    atoms[i].pos.y *= (unitCellvectors[0].y + unitCellvectors[1].y + unitCellvectors[2].y);
    atoms[i].pos.z *= (unitCellvectors[0].z + unitCellvectors[1].z + unitCellvectors[2].z);
    atoms[i].pos = retScaledVector(atoms[i].pos, 1.0);
  }

  // Store the volume of the unit cell to be used later in scaling the first ppParam
  params->vol[fileIndex] = getCellVolume(unitCellvectors);
  params->volPerAtom[fileIndex] = (params->vol[fileIndex] / (double)(params->na[fileIndex]));
  params->normalizedVolPerAtom[fileIndex] = params->volPerAtom[0] / params->volPerAtom[fileIndex];
  printf("\nCell %d volume = %.8f\n", fileIndex, params->vol[fileIndex]);
  printf("Cell %d volPerAtom = %.8f\n", fileIndex, params->volPerAtom[fileIndex]);
  printf("Cell %d normalizedVolPerAtom = %.8f\n", fileIndex, params->normalizedVolPerAtom[fileIndex]);

  return;
}

/*****************************************************************************************/
