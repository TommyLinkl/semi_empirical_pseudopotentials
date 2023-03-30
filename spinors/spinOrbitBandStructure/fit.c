/****************************************************************************/
/* This files calculates the band structure along specified directions in k-space */

#include "fit.h"

/****************************************************************************/

void runMonteCarlo(complexnumber *hamiltonian, complexnumber *kineticEnergyMat, complexnumber *localPotential, complexnumber *spinOrbitMatrix, double *vSpinOrbit, vector *basisStates, double *expBandStructure, vector *kPoint, double *bandStructure, double *energies, atom *atoms, param params) {
  FILE *pf;
  int i, j;
  double bestMSE, currentMSE, newMSE;
  double *tmpPseudos; /* Used to hold current pseudos while testing new ones */
  double *bestPseudos; 
  tmpPseudos = (double *) calloc(5*params.nAtoms, sizeof(double));
  bestPseudos = (double *) calloc(5*params.nAtoms, sizeof(double));
 
  pf = fopen("iterations.dat", "w");
  writeCurrentTime(pf);
  writeSeparation(pf); 

  calcBandStructure(bandStructure, hamiltonian, kineticEnergyMat, localPotential, spinOrbitMatrix, vSpinOrbit, basisStates, kPoint, energies, atoms, params);    
  bestMSE = currentMSE = calcMSE(bandStructure, expBandStructure, params.nKPoints * params.nBands);
  savePseudoParams(tmpPseudos, atoms, params); savePseudoParams(bestPseudos, atoms, params);
  writeIteration(currentMSE, 0, atoms, params, pf); 
  
  for (i = 0; i < params.nIterations; i++) {
    updatePseudoParams(atoms, params);
    calcBandStructure(bandStructure, hamiltonian, kineticEnergyMat, localPotential, spinOrbitMatrix, vSpinOrbit, basisStates, kPoint, energies, atoms, params);
    newMSE = calcMSE(bandStructure, expBandStructure, params.nKPoints * params.nBands);
    if ((i+1) % 100 == 0) writeIteration(newMSE, i+1, atoms, params, pf);
    if (newMSE < bestMSE) { 
      bestMSE = currentMSE = newMSE;
      savePseudoParams(bestPseudos, atoms, params); savePseudoParams(tmpPseudos, atoms, params);
    }
    else if (exp(-params.beta * (sqrt(newMSE) - sqrt(currentMSE))) > (0.5*(randNumber()+1.0))) {
       currentMSE = newMSE; 
       savePseudoParams(tmpPseudos, atoms, params);
    } 
    else { revertPseudoParams(atoms, tmpPseudos, params); } 
  }
  revertPseudoParams(atoms, bestPseudos, params);
  calcBandStructure(bandStructure, hamiltonian, kineticEnergyMat, localPotential, spinOrbitMatrix, vSpinOrbit, basisStates, kPoint, energies, atoms, params); 
  bestMSE = calcMSE(bandStructure, expBandStructure, params.nKPoints * params.nBands);
  writeResults(bandStructure, kPoint, bestMSE, atoms, params);
  writeCurrentTime(pf);
  fclose(pf);

  // Calculate and print the final r-space local pseudopotentials
  calcRealSpacePotential(atoms, params);

  // Free dynamically allocated memory
  free(tmpPseudos);
  free(bestPseudos);

  return;
}

/****************************************************************************/

void calcBandStructure(double *bandStructure, complexnumber *hamiltonian, complexnumber *kineticEnergyMat, complexnumber *localPotential, complexnumber *spinOrbitMatrix, double *vSpinOrbit, vector *basisStates, vector *kPoint, double *energies, atom *atoms, param params) {
  int i, j;
  vector kVect;

  calcLocalPotentialEnergyMatrix(localPotential, basisStates, atoms, params);
  for (i = 0; i < params.nKPoints; i++) {
    kVect = kPoint[i];
    calcSpinOrbitPotentialMatrix(spinOrbitMatrix, vSpinOrbit, basisStates, kVect, i, atoms, params); 
    calcKineticEnergyMatrix(kineticEnergyMat, basisStates, kVect, params);
    calcHamiltonianMatrix(hamiltonian, kineticEnergyMat, localPotential, spinOrbitMatrix, params);
    diagnolizeHamiltonian(hamiltonian, energies, params);
    storeEnergies(bandStructure, energies, i, params.nBands);
  }

  return;
}

/****************************************************************************/

void updatePseudoParams(atom *atoms, param params) {
  int i, j;

  for (i = 0; i < params.nAtoms; i++) {
    for (j = 0; j < 4; j++) {
      if (newAtomType(atoms, i)) {
        atoms[i].ppParams[j] += params.stepSize*randNumber(); /* randNumber returns a number between -1.0 and 1.0 */
        if (j == 0) atoms[i].spinOrbit += params.stepSize*randNumber();
      } 
      else { 
        copyPreviousAtomParams(atoms, i); 
      }
    }
  }

  return;
}

/****************************************************************************/

void copyPreviousAtomParams(atom *atoms, int i) {
  int j, k;

  for (j = 0; j < i; j++) {
    if (! strcmp(atoms[j].symbol, atoms[i].symbol)) {
      for (k = 0; k < 4; k++) {
        atoms[i].ppParams[k] = atoms[j].ppParams[k];
      }
      atoms[i].spinOrbit = atoms[j].spinOrbit;
    }
  }

  return;
}

/****************************************************************************/

void revertPseudoParams(atom *atoms, double *tmpPseudos, param params) {
  int i, j;

  for (i = 0; i < params.nAtoms; i++) {
    for (j = 0; j < 4; j++) atoms[i].ppParams[j] = tmpPseudos[i*5+j];
    atoms[i].spinOrbit = tmpPseudos[i*5+j];
  }

  return;
}

/****************************************************************************/

void savePseudoParams(double *tmpPseudos, atom *atoms, param params) {
  int i, j;

  for (i = 0; i < params.nAtoms; i++) {
    for (j = 0; j < 4; j++)  tmpPseudos[i*5+j] = atoms[i].ppParams[j];
    tmpPseudos[i*5+j] = atoms[i].spinOrbit;
  }

  return;
}

/****************************************************************************/

double calcMSE(double *array1, double *array2, int arrayLength) {
  int i; double mse = 0.0;
  
  for (i = 0; i < arrayLength; i++) mse += sqr((array1[i] - array2[i]));
  
  return mse / (double) (arrayLength);
}

/****************************************************************************/

void storeEnergies(double *bandStructure, double *energies, int index, int length) {
  int j;

  for (j = 0; j < length; j++) {
    bandStructure[index*length+j] = energies[j]*AUTOEV;
  }

  return;
}

/****************************************************************************/
