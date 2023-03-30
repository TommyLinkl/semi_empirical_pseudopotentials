#include "fd.h"

/*****************************************************************************/
// Normalizes to 1 all numStates states in psi each of which is of length ngrid 
// with a grid volume of dr 
// works with a complex wavefunctions (complex spinors as well)

void normalize_all(zomplex *psi, double dr, long numStates, long ngrid)
{
  long i;

  for (i = 0; i < numStates; i++) {
    normalize(&(psi[i]), dr, ngrid);
  }

  return;
}

/*****************************************************************************/
// Normalizes a single state of length ngrid and grid volume of dr
// works with complex a wavefunction (complex spinors as well)

double normalize(zomplex *psi, double dr, long ngrid)
{
  long i;
  double N = norm(psi,dr,ngrid);
  
  for (i = 0; i < ngrid; i++){
    psi[i].re /= N;
    psi[i].im /= N;
  }

  return (N);
}

/*****************************************************************************/
// Returns the norm of a single state of length ngrid and grid volume of dr
// works with complex a wavefunction (complex spinors as well)

double norm(zomplex *psi, double dr, long ngrid)
{
  long i;
  double norm = 0.0;

  for (i = 0; i < ngrid; i++) norm += (sqr(psi[i].re) + sqr(psi[i].im));
  
  return (sqrt(norm * dr));
}

/*****************************************************************************/