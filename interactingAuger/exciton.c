#include "fd.h"

/*****************************************************************************/

double retExcitonEnergy(double *Cbs, double *Hbs, long numBasisStates) {
  double energy = 0.0; 
  long i, j;

  for (i = 0; i < numBasisStates; i++) {
    for (j = 0; j < numBasisStates; j++) {
      energy += Cbs[j] * Cbs[i] * Hbs[i*numBasisStates + j];
	}
  }
  
  return (energy);
}

/*****************************************************************************/
