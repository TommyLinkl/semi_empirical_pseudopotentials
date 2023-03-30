/******************************************************************************/
//
// This is the beginning of this program that calculates the band structure 
// using local and nonlocal (spin-orbit) pseudopotentials 
//
/******************************************************************************/

#include "fit.h"

/******************************************************************************/
/* Initiates the main functions */

int main() {

  /****************************************************************************/
  /* Defines the variables that are used throughout this function */
  double *energies, homo, lumo, *expBandStructure, *bandStructure, *vSpinOrbit;
  complexnumber *hamiltonian, *kineticEnergyMat, *localPotential, *spinOrbitPot;
  vector *kPoint, *basisStates, *gVectors, *unitCellvectors, *tmpBV;
  atom *atoms;
  param params;

  /****************************************************************************/
  /* Allocates initial  memory */
  unitCellvectors = (vector *) calloc(NDIM, sizeof(vector));
  gVectors = (vector *) calloc(NDIM, sizeof(vector));

  /****************************************************************************/ 
  /* Reads input parameters that determine the size of the system */
  readInput(&params);
  
  /****************************************************************************/
  /* Allocates the input dependent memory */
  atoms = (atom *) calloc(params.nAtoms, sizeof(atom));
  expBandStructure = (double *) calloc(params.nKPoints * params.nBands, sizeof(double));
  bandStructure = (double *) calloc(params.nKPoints * params.nBands, sizeof(double));
  kPoint = (vector *) calloc(params.nKPoints, sizeof(vector));
 
  /****************************************************************************/
  /* Initializes the atom positions, initial pseudopotential parameters, and calculates the basis states */
  readExpBandStructure(expBandStructure, &params);
  readSystem(&params, unitCellvectors, atoms);
  readKPoints(kPoint, params);
  getRecipSpaceBasisVectors(unitCellvectors, gVectors, &params);
  tmpBV = (vector *) calloc(params.nMaxBV, sizeof(vector));
  calcBasisVectors(tmpBV, gVectors, &params);
  basisStates = (vector *) calloc(params.nBasisVectors, sizeof(vector));
  
  /****************************************************************************/
  // The matrices are 4 times as large when spin-orbit coupling is included
  vSpinOrbit = (double *) calloc(params.nKPoints*sqr(params.nBasisVectors), sizeof(double));
  kineticEnergyMat = (complexnumber *) calloc(4*sqr(params.nBasisVectors), sizeof(complexnumber));
  localPotential = (complexnumber *) calloc(4*sqr(params.nBasisVectors), sizeof(complexnumber));
  spinOrbitPot = (complexnumber *) calloc(4*sqr(params.nBasisVectors), sizeof(complexnumber));
  hamiltonian = (complexnumber *) calloc(4*sqr(params.nBasisVectors), sizeof(complexnumber));
  energies = (double *) calloc(2*params.nBasisVectors, sizeof(double));
  finalizeBasisStates(basisStates, tmpBV, &params);
  initAtomicFormFactors(params, atoms);
  
  /****************************************************************************/
  /* The computational work */
  calcVSpinOrbit(vSpinOrbit, basisStates, kPoint, params);
  runMonteCarlo(hamiltonian, kineticEnergyMat, localPotential, spinOrbitPot, vSpinOrbit, basisStates, expBandStructure, kPoint, bandStructure, energies, atoms, params);

  /****************************************************************************/
  /* Free dynamically allocated memory */ 
  free(atoms);  
  free(energies);
  free(unitCellvectors); free(gVectors);
  free(basisStates); free(tmpBV); free(vSpinOrbit);
  free(localPotential); free(kineticEnergyMat); free(hamiltonian); free(spinOrbitPot); 
  free(bandStructure); free(expBandStructure); free(kPoint);

  return 0;
}

/****************************************************************************/
