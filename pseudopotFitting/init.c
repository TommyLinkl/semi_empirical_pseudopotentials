/* This files does the required initialization before calculating the hamiltonian */

/*****************************************************************************************/

#include "fit.h"

/*****************************************************************************************/

int newAtomType(atom *atoms, int i) {
  int j;
 
  for (j = 0; j < i; j++) {
    if (! strcmp(atoms[j].symbol, atoms[i].symbol)) return 0;
  }
 
  return 1;
}

/*****************************************************************************************/

double getCellVolume(vector *cellVectors) {
  double cellVolume;
  vector temp;
  
  temp = retCrossProduct(cellVectors[1], cellVectors[2]);
  cellVolume = retDotProduct(cellVectors[0], temp);

  return (cellVolume);
}

/*****************************************************************************************/

int getRecipSpaceBasisVectors(vector *unitCellvectors, vector *gVectors, double maxKE) {
  double realSpacevolume = getCellVolume(&unitCellvectors[0]);
  double preFactor = TWOPI / realSpacevolume;

  gVectors[0] = retCrossProduct(unitCellvectors[1], unitCellvectors[2]);
  gVectors[0] = retScaledVector(gVectors[0], preFactor);
  gVectors[1] = retCrossProduct(unitCellvectors[2], unitCellvectors[0]);
  gVectors[1] = retScaledVector(gVectors[1], preFactor);
  gVectors[2] = retCrossProduct(unitCellvectors[0], unitCellvectors[1]);
  gVectors[2] = retScaledVector(gVectors[2], preFactor); 
  
  /* overestimate of the number of basis vectors */
  int i, n;
  double minGmag = 100000.0; 
  int nMaxBV = 0;
  for (i = 0; i < NDIM; i++) {
    if (gVectors[i].mag < minGmag) {
      minGmag = gVectors[i].mag;
      n = (int) (sqrt(maxKE / sqr(minGmag))) + 1;
      nMaxBV = (2*n+1) * (2*n+1) * (2*n+1);
    }
  }
  return nMaxBV;
}

/*****************************************************************************************/

void initAtomicFormFactors(atom *atoms, param params) {
  FILE *pf;
  int i, j, k, bandCopiedFrom;
  char fileName[1000];

  int ina = 0;
  for (k = 0; k < params.nBandStructures; k++) {
    for (i = 0; i < params.na[k]; i++) {
      if (newAtomType(atoms, ina)) {
        memset(&fileName[0], 0, sizeof(fileName));
        strcpy(fileName, atoms[ina].symbol);
        strcat(fileName, "Params.par");
        pf = fopen(fileName, "r");
        for (j = 0; j < NPParams; j++) fscanf(pf, "%lg", &atoms[ina].ppParams[j]);
        fclose(pf);
      }
      else {
        bandCopiedFrom = copyPreviousAtomParams(atoms, ina, params);
        if (bandCopiedFrom != k) {
          atoms[ina].ppParams[0] *= (params.normalizedVolPerAtom[k] / params.normalizedVolPerAtom[bandCopiedFrom]); 
        } 
      }
      ina++;
    }
  }
  writeSeparation(stdout);
  fprintf(stdout, "Initial pseudopotentials:\n");
  writePseudoParams(atoms, params, stdout);
  fflush(stdout);

  return;
} 

/*****************************************************************************************/
