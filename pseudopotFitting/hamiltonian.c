/*****************************************************************************************/

/* This file calculates the all the matrix elements for the Hamiltonian matrix */

/*****************************************************************************************/

#include "fit.h"

/*****************************************************************************************/

void calcHamiltonianMatrix(dcomplex *hamiltonian, dcomplex *potential, vector *basisStates, 
                            vector *kVector, int nbv, double kineticEnergyScaling) {
  int i, j;

  // Make sure hamiltonian is zero
  zeroComplexMatrix(hamiltonian, nbv, nbv);

  // Fill hamiltonian with the diagonal kinetic energy elements
  calcKineticEnergyMatrix(hamiltonian, basisStates, &kVector[0], nbv, kineticEnergyScaling);

  // Add the potential energy elements to the hamiltonian elements
  for (i = 0; i < nbv; i++) { 
    for (j = 0; j < nbv; j++) {
      hamiltonian[i*nbv+j].real += potential[i*nbv+j].real;
      hamiltonian[i*nbv+j].imag += potential[i*nbv+j].imag;
    }
  }

  return;
}

/*****************************************************************************************/

void calcKineticEnergyMatrix(dcomplex *hamiltonian, vector *basisStates, vector *kVector, 
                              int nbv, double kineticEnergyScaling) {
  int i;
  double preFactor;
  vector kPlusG;

  preFactor = ( kineticEnergyScaling * sqr(HBAR) / (2.0 * MASS) );

  // Fill hamiltonian with the kinetic energy matrix elements 
  // Diagonal elements = (hbar^2/2m)*(k+G)^2 and off-diagonal elements = 0
  for (i = 0; i < nbv; i++) {
    kPlusG = retAddedVectors(basisStates[i], kVector[0]);
    hamiltonian[i*nbv+i].real = kPlusG.mag*kPlusG.mag*preFactor;
  }    

  return;
}

/*****************************************************************************************/

void calcPotentialEnergyMatrix(dcomplex *potential, vector *basisStates, atom *atoms, int nbv, int nAtoms) { 
  int i, j, k, l;
  double *gDiffDotTau;
  vector gDiff; 
  dcomplex *structFact; 
  
  // Make sure potential is zero
  zeroComplexMatrix(potential, nbv, nbv);

  // Dynamically allocate memory
  gDiffDotTau = (double *) calloc(nAtoms, sizeof(double));
  structFact = (dcomplex *) calloc(nAtoms, sizeof(dcomplex));
 
  // Loop over all basisState pairs and each atom in the unit cell
  for (i = 0; i < nbv; i++) {
    for (j = 0; j < nbv; j++) {
      gDiff = retSubtractedVectors(basisStates[i], basisStates[j]);
      for (k = 0; k < nAtoms; k++) { 
        gDiffDotTau[k] = retDotProduct(gDiff, atoms[k].pos);
        structFact[k].real = (1.0 / (double) (nAtoms)) * cos(gDiffDotTau[k]);
        structFact[k].imag = (-1.0 / (double) (nAtoms)) * sin(gDiffDotTau[k]);
        atoms[k].fF = calcPot(gDiff.mag, atoms[k].ppParams);
        potential[i*nbv+j].real += atoms[k].fF * structFact[k].real;
        potential[i*nbv+j].imag += atoms[k].fF * structFact[k].imag; 
      }
    }  
  }

  // Free dynamically allocated memory
  free(gDiffDotTau);
  free(structFact);
  
  return;
}

/*****************************************************************************************/

double calcPot(double q, double *param) {
  double potential = (param[0]*(q*q - param[1]) / (param[2] * exp(param[3]*q*q) + 1.0));
  return potential;
}

/*****************************************************************************************/

void zeroComplexMatrix(dcomplex *mat, int n1, int n2) {
  int i, j;

  for (i = 0; i < n1; i++) {
    for (j = 0; j < n2; j++) {
      mat[i*n2+j].real = mat[i*n2+j].imag = 0.0;
    }
  }

  return;
}

/*****************************************************************************************/