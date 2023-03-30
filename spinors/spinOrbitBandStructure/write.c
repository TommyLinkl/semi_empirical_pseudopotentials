/****************************************************************************/
/* This file does the printing for the program */

#include "fit.h"

/****************************************************************************/

void writeResults(double *bandStructure, vector *kPoint, double bestMSE, atom *atoms, param params) {
  FILE *pf;
  int i, j, iNBands;
  double vbmEnergy, cbmEnergy, splitOffBandEnergy, gammaPointSpinOrbitSplitting, bandGap = 1000.0;
  double magBandGapKPoint = 1000.0;

  pf = fopen("bestBandStructure.dat", "w");
  for (i = 0; i < params.nKPoints; i++) {
    iNBands = i*params.nBands;
    fprintf(pf, "%.7f ", kPoint[i].mag);
    for (j = 0; j < params.nBands; j++) {
      fprintf(pf, "% .4f ", bandStructure[iNBands+j]);
    }
    fprintf(pf, "\n");
    // Determine the band gap and the k-point at which the band gap occurs    
    if ((bandStructure[iNBands + (params.nBands/2)] - bandStructure[iNBands + (params.nBands/2) - 1]) < bandGap) {
      vbmEnergy = bandStructure[iNBands + (params.nBands/2) - 1];
      cbmEnergy = bandStructure[iNBands + (params.nBands/2)];
      bandGap = (cbmEnergy - vbmEnergy);
      magBandGapKPoint = kPoint[i].mag;
    }
    // Determine the spin-orbit splitting at the gamma (k = (0,0,0)) point
    if (kPoint[i].mag < EPS) {
      splitOffBandEnergy = bandStructure[iNBands + (params.nBands/2) - 5];
      gammaPointSpinOrbitSplitting = (bandStructure[iNBands + (params.nBands/2) - 1] - splitOffBandEnergy);
    }
  }
  fclose(pf);

  pf = fopen("output.dat", "w");
  fprintf(pf, "RESULTS\n\n");
  fprintf(pf, "The best mean squared error found was: % .5f\n\n", bestMSE);
  fprintf(pf, "The magnitude of the k-point where the experimental band gap is: % .5f\n", params.expKp);
  fprintf(pf, "The experimental band gap is: % .3f eV\n", params.expEg);
  fprintf(pf, "The experimental energy of the conduction band minimum is: % .5f\n", params.expCBM);
  fprintf(pf, "The experimental energy of the valence band maximum is: % .5f\n", params.expVBM);
  fprintf(pf, "The experimental energy of the split off band maximum is: % .5f\n", params.expSOB);
  fprintf(pf, "The experimental spin-orbit splitting at the gamma point is: % .5f\n\n", params.expVBM-params.expSOB);  
  fprintf(pf, "The magnitude of the k-point where the calculated band gap is: % .5f\n", magBandGapKPoint);
  fprintf(pf, "The calculated band gap is: % .3f eV\n", bandGap);
  fprintf(pf, "The calculated energy of the conduction band minimum is: % .5f\n", cbmEnergy);
  fprintf(pf, "The calculated energy of the valence band maximum is: % .5f\n", vbmEnergy);
  fprintf(pf, "The calculated energy of the split off band maximum is: % .5f\n", splitOffBandEnergy);
  fprintf(pf, "The calculated spin-orbit splitting at the gamma point is: % .5f\n\n", gammaPointSpinOrbitSplitting);
  fprintf(pf, "The corresponding, optimized, parameters are:\n"); 
  writePseudoParams(atoms, params, pf);
  writeSeparation(pf);
  writeInput(atoms, params, pf);
  fclose(pf);

  return;
}

/****************************************************************************/

void writeIteration(double MSE, int iteration, atom *atoms, param params, FILE *pf) {
  
  if (iteration != 0) writeSeparation(pf);
  fprintf(pf, "Iteration %d\n", iteration);
  fprintf(pf, "MSE % .5f\n", MSE);
  writePseudoParams(atoms, params, pf);
  if (iteration == params.nIterations) writeSeparation(pf);
  fprintf(pf, "\n");  
  fflush(pf);

  return;
}

/****************************************************************************/

void writePseudoParams(atom *atoms, param params, FILE *pf) {
  int i, j;
  
  for (i = 0; i < params.nAtoms; i++) {
    if (newAtomType(atoms, i)) {
      fprintf(pf, "\n%s", atoms[i].symbol);
      for (j = 0; j < 4; j++) fprintf(pf, " % .8f", atoms[i].ppParams[j]);
      fprintf(pf, " % .8f", atoms[i].spinOrbit);
    }
  }

  return;
}

/****************************************************************************/

void writeVector(vector vect, FILE *pf) {
  fprintf(pf, "% .10f % .10f % .10f % .10f\n", vect.x, vect.y, vect.z, vect.mag);

  return;
}

/****************************************************************************/

void writeInput(atom *atoms, param params, FILE *pf) {
  int i;

  fprintf(pf, "INPUT\n\n");
  fprintf(pf, "The total number of iterations was: %d\n", params.nIterations);
  fprintf(pf, "The beta used was: % .1f\n\n", params.beta);
  fprintf(pf, "The cell volume was: % .3f\n", params.cellVolume);
  fprintf(pf, "The number of atoms in the cell was: %d\n", params.nAtoms);
  fprintf(pf, "The atoms were place at the following positions:\n\n");
  for (i = 0; i < params.nAtoms; i++) {
    fprintf(pf, "%s %.5f %.5f %.5f\n", atoms[i].symbol, atoms[i].pos.x, atoms[i].pos.y, atoms[i].pos.z);
  } fprintf(pf, "\n");
  fprintf(pf, "The maximum kinetic energy was: % .5f\n", params.maxKE);
  fprintf(pf, "The kinetic energy scaling was: % .5f\n", params.kineticEnergyScaling);
  fprintf(pf, "The number of basis functions used was: %d\n\n", params.nBasisVectors);
  fprintf(pf, "The number of bands calculated was: %d\n", params.nBands);
  fprintf(pf, "The number of k-points that were used to fit the pseudopotentials to was: %d\n", params.nKPoints);
  writeSeparation(pf);  

  return;
}

/****************************************************************************/
// prints out current time to pf 
 
void writeCurrentTime(FILE *pf) {
  time_t startTime;

  startTime = time(NULL);
  fprintf(pf, ctime(&startTime));

  return;
}

/****************************************************************************/

void writeSeparation(FILE *pf) {
  fprintf(pf, "\n\n******************************************************************************\n\n");
  
  return;
}

/****************************************************************************/
