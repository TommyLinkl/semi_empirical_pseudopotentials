/* This is the beginning of this program that calculates the band structure using purely local pseudopotentials */

/*****************************************************************************************/

#include "fit.h"

/*****************************************************************************************/

/* Initiates the main functions */

int main() {

  /* Defines the variables that are used throughout this function */
  int i;
  double *energies, *expBandStructure, *bandStructure;
  dcomplex *hamiltonian, *potential;
  vector *kPoint, *basisStates, *gVectors, *unitCellvectors, *tmpBV;
  atom *atoms;
  param params;

  /* Reads input parameters that determine the size of the system */
  writeCurrentTime(stdout);
  readInput(&params);
  int nbs = params.nBandStructures;

  /* Allocates the input dependent memory */
  unitCellvectors = (vector *) calloc(NDIM*nbs, sizeof(vector));
  gVectors = (vector *) calloc(NDIM*nbs, sizeof(vector));
  atoms = (atom *) calloc(params.nAtoms, sizeof(atom));
  expBandStructure = (double *) calloc(params.nKPE, sizeof(double));
  bandStructure = (double *) calloc(params.nKPE, sizeof(double));
  kPoint = (vector *) calloc(params.nKPoints, sizeof(vector));
 
  /* Initializes the atom positions, initial pseudopotential parameters, and calculates the basis states */
  int ina = 0;
  int ine = 0;
  int ikp = 0;
  int iebs = 0;
  for (i = 0; i < nbs; i++) {
    readExpBandStructure(&expBandStructure[iebs], params.nkp[i], params.ne[i], i);
    readSystem(&params, &unitCellvectors[i*NDIM], &atoms[ina], i);
    readKPoints(&kPoint[ikp], params.nkp[i], params.s[i], i);
    params.nmbv[i] = getRecipSpaceBasisVectors(&unitCellvectors[i*NDIM], &gVectors[i*NDIM], params.mke[i]);
    ina += params.na[i];
    ine += params.ne[i];
    ikp += params.nkp[i];
    iebs += params.nkp[i]*params.ne[i];
  }
  initAtomicFormFactors(atoms, params);
  
  int inmbv = 0;
  for (i = 0; i < nbs; i++) inmbv += params.nmbv[i];
  tmpBV = (vector *) calloc(inmbv, sizeof(vector));

  int inbv = 0;
  int inbv_squared = 0;
  inmbv = 0;
  for (i = 0; i < nbs; i++) {
    params.nbv[i] = calcBasisVectors(&tmpBV[inmbv], &gVectors[i*NDIM], params.mke[i]);
    inmbv += params.nmbv[i];
    inbv += params.nbv[i];
    inbv_squared += params.nbv[i]*params.nbv[i];
  }
  
  basisStates = (vector *) calloc(inbv, sizeof(vector));
  potential = (dcomplex *) calloc(inbv_squared, sizeof(dcomplex));
  hamiltonian = (dcomplex *) calloc(inbv_squared, sizeof(dcomplex));
  energies = (double *) calloc(inbv, sizeof(double));
  
  inbv = 0;
  inmbv = 0;
  for (i = 0; i < nbs; i++) {
    finalizeBasisStates(&basisStates[inbv], &tmpBV[inmbv], params.nmbv[i], params.mke[i]);
    inmbv += params.nmbv[i];
    inbv += params.nbv[i];
  }

  /* Begins the computational work */
  runMonteCarlo(hamiltonian, potential, basisStates, expBandStructure, kPoint, bandStructure, energies, atoms, params);
  writeCurrentTime(stdout);

  /* Free dynamically allocated memory */ 
  free(atoms);  
  free(energies); 
  free(unitCellvectors); free(gVectors);
  free(basisStates); free(potential); free(hamiltonian);
  free(tmpBV);
  free(bandStructure); free(expBandStructure); free(kPoint);

  return 0;
}

/*****************************************************************************************/
