/****************************************************************************/
/* This files does the required initialization before calculating the hamiltonian */

#include "fit.h"

/****************************************************************************/

int newAtomType(atom *atoms, int i) {
  int j;
 
  for (j = 0; j < i; j++) {
    if (! strcmp(atoms[j].symbol, atoms[i].symbol)) return 0;
  }
 
  return 1;
}

/****************************************************************************/

double getCellVolume(vector *cellVectors) {
  double cellVolume;
  vector temp;
  
  temp = retCrossProduct(cellVectors[1], cellVectors[2]);
  cellVolume = retDotProduct(cellVectors[0], temp);
  
  return (cellVolume);
}

/****************************************************************************/

void getRecipSpaceBasisVectors(vector *unitCellvectors, vector *gVectors, param *params) {
  params->cellVolume = getCellVolume(&unitCellvectors[0]);
  double preFactor = TWOPI / params->cellVolume;

  gVectors[0] = retCrossProduct(unitCellvectors[1], unitCellvectors[2]);
  gVectors[0] = retScaledVector(gVectors[0], preFactor);
  gVectors[1] = retCrossProduct(unitCellvectors[2], unitCellvectors[0]);
  gVectors[1] = retScaledVector(gVectors[1], preFactor);
  gVectors[2] = retCrossProduct(unitCellvectors[0], unitCellvectors[1]);
  gVectors[2] = retScaledVector(gVectors[2], preFactor); 
  
  /* overestimate of the number of basis vectors */
  double minGsqrMag = 100000.0; int n, i;
  for (i = 0; i < NDIM; i++) {
    if (sqr(gVectors[i].mag) < minGsqrMag) {
      minGsqrMag = sqr(gVectors[i].mag);
      n = (int) (sqrt(params->maxKE / minGsqrMag)) + 1;
      params->nMaxBV = (2*n+1) * (2*n+1) * (2*n+1);
    }
  }
  return;
}

/****************************************************************************/

void initAtomicFormFactors(param params, atom *atoms) {
  FILE *pf;
  int i, j;
  char fileName[1000];
 
  for (i = 0; i < params.nAtoms; i++) {
    memset(&fileName[0], 0, sizeof(fileName));
    strcpy(fileName, atoms[i].symbol);
    strcat(fileName, "Params.par");
    pf = fopen(fileName, "r");
    for (j = 0; j < 5; j++) {
      if (j < 4) { 
        fscanf(pf, "%lg", &atoms[i].ppParams[j]); 
      }
      else { 
        fscanf(pf, "%lg", &atoms[i].spinOrbit); 
      }
    }
    fclose(pf);
  }

  return;
} 

/****************************************************************************/
