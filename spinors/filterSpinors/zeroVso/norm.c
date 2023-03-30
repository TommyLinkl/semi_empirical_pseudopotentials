/*****************************************************************************/
// File contains functions that are involved with normalizing wavefunctions

#include "fd.h"

/*****************************************************************************/
// Returns the norm of a single state of length ngrid and grid volume of dr
// works with complex a wavefunction (complex spinors as well)

double norm(zomplex *psi, double dr,long ngrid,long nthreads)
{
  long i;
  double norm = 0.0;

  omp_set_dynamic(0);
  omp_set_num_threads(nthreads);
#pragma omp parallel for reduction(+:norm)
  for (i = 0; i < ngrid; i++) norm += (sqr(psi[i].re) + sqr(psi[i].im));
  
  return (sqrt(norm * dr));
}

/*****************************************************************************/
// Normalizes a single state of length ngrid and grid volume of dr
// works with complex a wavefunction (complex spinors as well)

double normalize(zomplex *psi,double dr,long ngrid,long nthreads)
{
  long k;
  double N = norm(psi,dr,ngrid,nthreads);
  
  omp_set_dynamic(0);
  omp_set_num_threads(nthreads);
#pragma omp parallel for 
  for (k = 0; k < ngrid; k++){
    psi[k].re /= N;
    psi[k].im /= N;
  }

  return (N);
}

/*****************************************************************************/
// Normalizes to 1 all numStates states in psi each of which is of length ngrid 
// with a grid volume of dr works with a complex wavefunctions (complex spinors as well)

void normalize_all(zomplex *psi,double dr,long ms,long ngrid,long nthreads)
{
  long i, ie, ieg;
  double nrm;

  omp_set_dynamic(0);
  omp_set_num_threads(nthreads);
#pragma omp parallel for private(ie,nrm,ieg,i)
  for (ie = 0; ie < ms; ie++){
    for (ieg = ie * ngrid, nrm = 0.0, i = 0; i < ngrid; i++) nrm += sqr(psi[ieg+i].re) + sqr(psi[ieg+i].im);
    if (nrm != 0.0) nrm = 1.0 / sqrt(nrm * dr);
    for (ieg = ie * ngrid, i = 0; i < ngrid; i++) {
      psi[ieg+i].re *= nrm;
      psi[ieg+i].im *= nrm;
    }
  }

  return;
}

/*****************************************************************************/