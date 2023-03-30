/* This file does the printing for the program */

#include "fit.h"

/*****************************************************************************************/

void writeResults(double *bandStructure, vector *kPoint, double bestMSE, atom *atoms, param params) {
  FILE *pf;
  int i, j, k;
  double bandGap[params.nBandStructures+1];
  double magBandGapKPoint[params.nBandStructures+1];
  char fileName[100], tmpChar[100];  

  // Write kpoints.dat files
  int ikp = 0;
  for (k = 0; k < params.nBandStructures; k++) {
    strcpy(fileName, "kpoints_");
    sprintf(tmpChar, "%d.dat", k);
    strcat(fileName, tmpChar);
    pf = fopen(fileName, "w");
    for (i = 0; i < params.nkp[k]; i++) {
      writeVector(kPoint[ikp + i], pf);
    }
    fclose(pf);
    ikp += params.nkp[k];
  }

  // Final results to screen
  writeSeparation(stdout);
  printf("The best mean squared error found was: % .5f\n", bestMSE);
  writePseudoParams(atoms, params, stdout);
  writeSeparation(stdout);

  // Write best pseudoparams to new files
  writeBestPseudoParams(atoms, params);

  // Final results to output.dat
  pf = fopen("output.dat", "w");
  writeCurrentTime(pf);
  writeSeparation(pf);
  fprintf(pf, "The best mean squared error found was: % .5f\n", bestMSE);
  writePseudoParams(atoms, params, pf);
  writeSeparation(pf);
  writeInput(atoms, params, pf); 
  fclose(pf);

  // Write bandStructure.dat files
  ikp = 0;
  int ikpe = 0;
  for (k = 0; k < params.nBandStructures; k++) {
    strcpy(fileName, "bandStruct_");
    sprintf(tmpChar, "%d.dat", k);
    strcat(fileName, tmpChar);
    pf = fopen(fileName, "w");
    for (i = 0; i < params.nkp[k]; i++) {
      fprintf(pf, "% .5f ", kPoint[ikp + i].mag);
      for (j = 0; j < params.ne[k]; j++) {
        fprintf(pf, "% .5f ", bandStructure[ikpe + i*params.ne[k] + j]);
      }
      fprintf(pf, "\n");
    }
    fclose(pf);
    ikp += params.nkp[k];
    ikpe += params.nkp[k]*params.ne[k];
  }

  return;
}

/*****************************************************************************************/

void writeIteration(double MSE, int iteration, atom *atoms, param params, FILE *pf) {
  if (iteration != 0) writeSeparation(pf);
  fprintf(pf, "Iteration %d\n", iteration);
  fprintf(pf, "MSE % .5f\n", MSE);
  writePseudoParams(atoms, params, pf);
  if (iteration == params.nIterations) writeSeparation(pf);  

  return;
}

/*****************************************************************************************/

void writeBestPseudoParams(atom *atoms, param params) {
	FILE *pf;
	int i, j, k;
  char fileName[100], tmpChar[10];

  int ina = 0;
  for (k = 0; k < params.nBandStructures; k++) {
    for (i = 0; i < params.na[k]; i++) {
      strcpy(fileName, "best_");
      sprintf(tmpChar, "%s_%d_Params.dat", atoms[ina].symbol, k);
      strcat(fileName, tmpChar);
      pf = fopen(fileName, "w");
      for (j = 0; j < NPParams; j++) {
        fprintf(pf, "% .8f\n", atoms[ina].ppParams[j]);
      }
      fclose(pf);
      ina++;
    }
  }

	return;
}

/*****************************************************************************************/

void writePseudoParams(atom *atoms, param params, FILE *pf) {
  int i, j, k;
  
  int ina = 0;
  for (k = 0; k < params.nBandStructures; k++) {
    if (k) fprintf(pf, "\n");
    fprintf(pf, "\nBand structure %d\n", k);
    for (i = 0; i < params.na[k]; i++) {
      fprintf(pf, "\n%s", atoms[ina].symbol);
      for (j = 0; j < NPParams; j++) fprintf(pf, " % .8f", atoms[ina].ppParams[j]);
      ina++;
    }
  }

  return;
}

/*****************************************************************************************/

void writeInput(atom *atoms, param params, FILE *pf) {
  int i;

  fprintf(pf, "INPUT\n\n");
  fprintf(pf, "The number of band structures fit to was: %d\n", params.nBandStructures);
  fprintf(pf, "The kinetic energy scaling was: % .4f\n", params.kineticEnergyScaling);
  fprintf(pf, "The total number of iterations was: %d\n", params.nIterations);
  fprintf(pf, "The beta used was: % .1f\n\n", params.beta);
  fprintf(pf, "The atoms were place at the following positions:\n\n");
  for (i = 0; i < params.nAtoms; i++) {
    fprintf(pf, "%s %.5f %.5f %.5f\n", atoms[i].symbol, atoms[i].pos.x, atoms[i].pos.y, atoms[i].pos.z);
  } 
  for (i = 0; i < params.nBandStructures; i++) {
    fprintf(pf, "\nBand structure %d\n", i);
    fprintf(pf, "The number of atoms in the cell was: %d\n", params.na[i]);
    fprintf(pf, "The number of bands calculated was: %d\n", params.ne[i]);
    fprintf(pf, "The number of k-points that were used to fit the pseudopotentials to was: %d\n", params.nkp[i]);
    fprintf(pf, "The maximum kinetic energy was: % .3f\n", params.mke[i]);
    fprintf(pf, "The number of basis functions used was: %d\n", params.nbv[i]);
  }
  writeSeparation(pf);  

  return;
}

/*****************************************************************************************/

void writeVector(vector vect, FILE *pf) {
  fprintf(pf, "% .10f % .10f % .10f % .10f\n", vect.x, vect.y, vect.z, vect.mag);

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

/*****************************************************************************************/

void writeSeparation(FILE *pf) {
  fprintf(pf, "\n\n***************************************************************\n\n");
  
  return;
}

/*****************************************************************************************/
