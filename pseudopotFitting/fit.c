/* This files calculates the band structure along specified directions in k-space */

/*****************************************************************************************/

#include "fit.h"

/*****************************************************************************************/

void runMonteCarlo(dcomplex *hamiltonian, dcomplex *potential, vector *basisStates, double *expBandStructure, 
  vector *kPoint, double *bandStructure, double *energies, atom *atoms, param params) {
  FILE *pf;
  int i, j;
  double bestMSE, currentMSE, newMSE;
  double *tmpPseudos; /* Used to hold current pseudos while testing new ones */
  double *bestPseudos; 
  tmpPseudos = (double *) calloc(NPParams*params.nAtoms, sizeof(double));
  bestPseudos = (double *) calloc(NPParams*params.nAtoms, sizeof(double));

  pf = fopen("iterations.dat", "w"); 
  writeCurrentTime(pf);
  writeSeparation(pf);
  calcBandStructure(bandStructure, hamiltonian, potential, basisStates, kPoint, energies, atoms, params);  
  bestMSE = currentMSE = calcAllIndividualWeightedMSE(bandStructure, expBandStructure, params);
  savePseudoParams(tmpPseudos, atoms, params); 
  savePseudoParams(bestPseudos, atoms, params);
  writeIteration(currentMSE, 0, atoms, params, pf);  
  fclose(pf);

  for (i = 0; i < params.nIterations; i++) {
    updatePseudoParams(atoms, params);
    calcBandStructure(bandStructure, hamiltonian, potential, basisStates, kPoint, energies, atoms, params);
    newMSE = calcAllIndividualWeightedMSE(bandStructure, expBandStructure, params);
    if ((i+1) % 10 == 0) {
	  pf = fopen("iterations.dat", "a");
	  writeIteration(newMSE, i+1, atoms, params, pf);
      fclose(pf);
	}
	if (newMSE < bestMSE) { 
      bestMSE = currentMSE = newMSE;
      savePseudoParams(bestPseudos, atoms, params); 
      savePseudoParams(tmpPseudos, atoms, params);
    }
    else if (exp(-params.beta * (sqrt(newMSE) - sqrt(currentMSE))) > (0.5*(randNumber()+1.0))) {
      currentMSE = newMSE; 
      savePseudoParams(tmpPseudos, atoms, params);
    } 
    else { // reject the move
      revertPseudoParams(atoms, tmpPseudos, params); 
    } 
  }
  revertPseudoParams(atoms, bestPseudos, params);
  calcBandStructure(bandStructure, hamiltonian, potential, basisStates, kPoint, energies, atoms, params); 
  bestMSE = calcAllIndividualWeightedMSE(bandStructure, expBandStructure, params);
  pf = fopen("iterations.dat", "a");
  writeCurrentTime(pf);
  fclose(pf);

  // Write the results and print the real space pseudopotentials
  writeResults(bandStructure, kPoint, bestMSE, atoms, params);
  calcRealSpacePotential(atoms, params);

  // Free dynamically allocated memory
  free(tmpPseudos);
  free(bestPseudos);

  return;
}

/*****************************************************************************************/

void calcBandStructure(double *bandStructure, dcomplex *hamiltonian, dcomplex *potential, 
                      vector *basisStates, vector *kPoint, double *energies, atom *atoms, param params) {
  int i, j, k;

  int inbv = 0;
  int inbv_squared = 0;
  int ikp = 0;
  int ina = 0;
  int ikpe = 0;
  for (i = 0; i < params.nBandStructures; i++) {
    calcPotentialEnergyMatrix(&potential[inbv_squared], &basisStates[inbv], &atoms[ina], params.nbv[i], params.na[i]);
    for (j = 0; j < params.nkp[i]; j++) {
      calcHamiltonianMatrix(&hamiltonian[inbv_squared], &potential[inbv_squared], &basisStates[inbv], 
                              &kPoint[j+ikp], params.nbv[i], params.kineticEnergyScaling);
      diagonalizeHamiltonian(&hamiltonian[inbv_squared], &energies[0], params.nbv[i]);
      // store the lowest ne (number of electrons) energies
      for (k = 0; k < params.ne[i]; k++) bandStructure[ikpe + j*params.ne[i] + k] = energies[k]*AUTOEV;
    }
    inbv += params.nbv[i];
    inbv_squared += params.nbv[i]*params.nbv[i];
    ikp += params.nkp[i];
    ina += params.na[i];
    ikpe += params.nkp[i]*params.ne[i];
  }

  return;
}

/*****************************************************************************************/

double calcAllIndividualWeightedMSE(double *bandStructure, double *expBandStructure, param params) {
  int i, j; 
  double mse[params.nBandStructures], totalMSE = 0.0;
  double bandWeight[params.nElectrons];
  
  int ikpe = 0;
  for (i = 0; i < params.nBandStructures; i++) {
    // Weight bands closest to band-edge more than those far apart
    for (j = 0; j < params.ne[i]/2; j++) {
      if (j < params.ne[i]/8) bandWeight[j] = bandWeight[params.ne[i]-j-1] = 0.06125; 
      else if (j < params.ne[i]/4) bandWeight[j] = bandWeight[params.ne[i]-j-1] = 2.0;
      else bandWeight[j] = bandWeight[params.ne[i]-j-1] = 4.0;
    }
    // Calculated weighted mean squared error between bandStructure and expBandStructure
    mse[i] = calcWeightedMSE(&bandStructure[ikpe], &expBandStructure[ikpe], &bandWeight[0], params.ne[i], params.nkp[i]);
    totalMSE += mse[i];
    ikpe += params.ne[i]*params.nkp[i];
  }
  
  return (totalMSE / (double)(params.nBandStructures));
}

/*****************************************************************************************/

double calcWeightedMSE(double *array1, double *array2, double *weights, int nElectrons, int nKPoints) {
  int i, j, index, n = nElectrons*nKPoints; 
  double mse = 0.0;
  
  for (i = 0; i < nKPoints; i++) {
    for (j = 0; j < nElectrons; j++) {
      index = i*nElectrons+j; 
      if ((array2[index] > EPS) || (array2[index] < -EPS)) {
        mse += weights[j] * sqr((array1[index] - array2[index]));
      }
      else {
        n--;
      }
    }
  }
  
  return (mse / (double)(n));
}


/*****************************************************************************************/

double calcAllIndividualMSE(double *bandStructure, double *expBandStructure, param params) {
  int i; 
  double mse[params.nBandStructures], totalMSE = 0.0;
  
  int ikpe = 0;
  for (i = 0; i < params.nBandStructures; i++) {
    mse[i] = calcMSE(&bandStructure[ikpe], &expBandStructure[ikpe], params.ne[i]*params.nkp[i]);
    totalMSE += mse[i];
    ikpe += params.ne[i]*params.nkp[i];
  }
  
  return (totalMSE / (double)(params.nBandStructures));
}

/*****************************************************************************************/

double calcMSE(double *array1, double *array2, int arrayLength) {
  int i; 
  double mse = 0.0;
  
  for (i = 0; i < arrayLength; i++) {
  	if ((array2[i] > EPS) || (array2[i] < -EPS)) {
  	  mse += sqr((array1[i] - array2[i]));
    }
  	else {
  	  arrayLength--;
  	}
  }
  
  return (mse / (double)(arrayLength));
}

/*****************************************************************************************/

void storeEnergies(double *bandStructure, double *energies, int index, int length) {
  int i;

  for (i = 0; i < length; i++) {
    bandStructure[index*length+i] = energies[i]*AUTOEV;
  }

  return;
}

/*****************************************************************************************/

void updatePseudoParams(atom *atoms, param params) {
  int i, j, k, bandCopiedFrom;

  int ina = 0;
  for (k = 0; k < params.nBandStructures; k++) {
    for (i = 0; i < params.na[k]; i++) {
      for (j = 0; j < NPParams; j++) {
        if (newAtomType(atoms, ina)) {
          // randomly change one of the parameters, randNumber returns a number between -1.0 and 1.0 
          atoms[ina].ppParams[j] += params.stepSize*randNumber();
        } 
        else { 
          bandCopiedFrom = copyPreviousAtomParams(atoms, ina, params);
          if (bandCopiedFrom != k) {
            atoms[ina].ppParams[0] *= (params.normalizedVolPerAtom[k] / params.normalizedVolPerAtom[bandCopiedFrom]); 
          }
        }
      }
      ina++;
    }
  }

  return;
}

/*****************************************************************************************/

int copyPreviousAtomParams(atom *atoms, int i, param params) {
  int j, k, l;

  for (j = 0; j < i; j++) {
    if (! strcmp(atoms[j].symbol, atoms[i].symbol)) {
      for (k = 0; k < NPParams; k++) {
        atoms[i].ppParams[k] = atoms[j].ppParams[k];
      }
      break;
    }
  }

  // determine which band structure the previous atom was a part of 
  // and return the corresponding band index
  int ina = 0;
  for (l = 0; l < params.nBandStructures; l++) {
    ina += params.na[l];
    if (j < ina) {
      return l;
    }
  }

  return (-1); // should never reach here, -1 signifies unsuccesful
}

/*****************************************************************************************/

void revertPseudoParams(atom *atoms, double *tmpPseudos, param params) {
  int i, j;

  for (i = 0; i < params.nAtoms; i++) {
    for (j = 0; j < NPParams; j++) {
      atoms[i].ppParams[j] = tmpPseudos[i*NPParams+j];
    }
  }

  return;
}

/*****************************************************************************************/

void savePseudoParams(double *tmpPseudos, atom *atoms, param params) {
  int i, j;

  for (i = 0; i < params.nAtoms; i++) {
    for (j = 0; j < NPParams; j++)  {
      tmpPseudos[i*NPParams+j] = atoms[i].ppParams[j];
    }
  }

  return;
}

/*****************************************************************************************/
