/****************************************************************************/

#include "fit.h"

/****************************************************************************/

void readInput(param *params) {
  FILE *pf;
  int i, j = 0;
  char field[1000], tmp[1000];

  // Set default
  params->kineticEnergyScaling = 1.000;

  pf = fopen("input.par", "r");
  while (1) {
    fscanf(pf, "%s", field);
    if (! strcmp(field, "maxke")) fscanf(pf, "%s %lg", tmp, &params->maxKE);
    else if (! strcmp(field, "natoms")) fscanf(pf, "%s %d", tmp, &params->nAtoms);
    else if (! strcmp(field, "nbands")) fscanf(pf, "%s %d", tmp, &params->nBands);
    else if (! strcmp(field, "niterations")) fscanf(pf, "%s %d", tmp, &params->nIterations);
    else if (! strcmp(field, "nkpoints")) fscanf(pf, "%s %d", tmp, &params->nKPoints);
    else if (! strcmp(field, "beta")) fscanf(pf, "%s %lg", tmp, &params->beta);
    else if (! strcmp(field, "stepsize")) fscanf(pf, "%s %lg", tmp, &params->stepSize);
    else if (! strcmp(field, "kineticEnergyScaling")) fscanf(pf, "%s %lg", tmp, &params->kineticEnergyScaling);
    else {
       printf("Invalid input field and/ or format of input\n");
       exit(EXIT_FAILURE);
    }
    i++;
    if (i > 8) break;
  }
  fflush(pf);
 
  return;
}

/****************************************************************************/

void readSystem(param *params, vector *unitCellvectors, atom *atoms) {
  FILE *pfile;
  int i;
  char buf[1000]; char tmp[1000];
 
  params->nAtomTypes = 0;
  pfile = fopen("system.par", "r");
  while (fscanf(pfile, "%s", buf) == 1) {
    if (! strcmp(buf, "atoms")) {
      for (i = 0; i < params->nAtoms; i++) {
        fscanf (pfile, "%s %lg %lg %lg", atoms[i].symbol, &atoms[i].pos.x, &atoms[i].pos.y, &atoms[i].pos.z);
      if (newAtomType(atoms, i)) params->nAtomTypes++;
      }
    }
    else if (! strcmp(buf, "cell")) {
      for (i = 0; i < NDIM; i++) {
        fscanf(pfile, "%lg %lg %lg", &unitCellvectors[i].x, &unitCellvectors[i].y, &unitCellvectors[i].z);
      }
    }
    else if (! strcmp(buf, "scale")) fscanf(pfile, "%s %lg", tmp, &params->scale);
  }
  fflush(pfile);
 
  /* Scale the unit cell vectors and atom position vectors in the input file by the input parameter scale */
  for (i = 0; i < NDIM; i++) {
    unitCellvectors[i] = retScaledVector(unitCellvectors[i], params->scale);
  }
 
  /* Scale the atom positions by the unit cell vectors */
  for (i = 0; i < params->nAtoms; i++) {
    atoms[i].pos.x *= (unitCellvectors[0].x + unitCellvectors[1].x + unitCellvectors[2].x);
    atoms[i].pos.y *= (unitCellvectors[0].y + unitCellvectors[1].y + unitCellvectors[2].y);
    atoms[i].pos.z *= (unitCellvectors[0].z + unitCellvectors[1].z + unitCellvectors[2].z);
    atoms[i].pos = retScaledVector(atoms[i].pos, 1.0);
  }

  return;
}

/****************************************************************************/

void readExpBandStructure(double *expBandStructure, param *params) {
  FILE *pf;
  int i, j;
  double tmp, eg, kp, cbmEnergy, vbmEnergy, splitOffBandEnergy;  

  eg = 10000.0;
  kp = 10000.0;
  pf = fopen("so_expBandStruct.par", "r");
  for (i = 0; i < params->nKPoints; i++) {
    fscanf(pf, "%lg", &tmp);
    for (j = 0; j < params->nBands; j++) fscanf(pf, "%lg", &expBandStructure[i*params->nBands+j]); 
    if ((expBandStructure[i*params->nBands+(params->nBands/2)] - expBandStructure[i*params->nBands+(params->nBands/2-1)]) < eg) {
      cbmEnergy = expBandStructure[i*params->nBands+(params->nBands/2)];
      vbmEnergy = expBandStructure[i*params->nBands+(params->nBands/2-1)];
      eg = cbmEnergy - vbmEnergy;
      kp = tmp;
    }
    if (tmp < EPS) {
      splitOffBandEnergy = expBandStructure[i*params->nBands+(params->nBands/2) - 5];
    }
  }
  params->expEg = eg;
  params->expKp = kp;
  params->expCBM = cbmEnergy;
  params->expVBM = vbmEnergy;
  params->expSOB = splitOffBandEnergy;
  
  fclose(pf);  
  
  return;
}

/****************************************************************************/

void readKPoints(vector *kPoint, param params) {
  FILE *pf;
  int i;
 
  pf = fopen("kpoints.par", "r"); 
  for (i = 0; i < params.nKPoints; i++) {
    fscanf(pf, "%lg %lg %lg", &kPoint[i].x, &kPoint[i].y, &kPoint[i].z);
    kPoint[i] = retScaledVector(kPoint[i], (TWOPI / params.scale));
  }
  
  return;
}

/****************************************************************************/


