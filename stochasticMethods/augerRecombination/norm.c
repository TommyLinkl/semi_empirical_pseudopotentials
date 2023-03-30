#include "ar.h"

/*****************************************************************************/

double norm(zomplex *psi, double dr,long nGridPoints)
{
  long i;
  double norm = 0.0;

  for (i = 0; i < nGridPoints; i++) norm += (sqr(psi[i].re) + sqr(psi[i].im));
  return (sqrt(norm * dr));
}

/*****************************************************************************/

double normalize(zomplex *psi,double dr,long nGridPoints)
{
  long k;
  double N_1 = 1.0 / norm(psi,dr,nGridPoints);
  
  for (k = 0; k < nGridPoints; k++){
    psi[k].re *= N_1;
    psi[k].im *= N_1;
  }
  return (1.0/N_1);
}

/*****************************************************************************/

void normalize_all(double *psi, double dr, long ms, long nGridPoints)
{
  FILE *pf;
  long i, ie, ieg;
  double nrm;

  pf = fopen("norm.dat" , "w");
  for (ie = 0; ie < ms; ie++){
    for (ieg = ie * nGridPoints, nrm = 0.0, i = 0; i < nGridPoints; i++)
      nrm += sqr(psi[ieg+i]);
    if (nrm != 0.0){
      nrm = 1.0 / sqrt(nrm * dr);
      fprintf (pf,"%ld norm = %g\n",ie,nrm); fflush(0);
    }
    for (ieg = ie * nGridPoints, i = 0; i < nGridPoints; i++)
      psi[ieg+i] *= nrm;
  }
  fclose(pf);
  return;
}

/*****************************************************************************/

void normalize_all_omp(double *psi, double dr, long ms, long nGridPoints, long nThreads)
{
  long i, ie, ieg;
  double nrm;

  omp_set_dynamic(0);
  omp_set_num_threads(nThreads);
#pragma omp parallel for private(ie,nrm,ieg,i)
  for (ie = 0; ie < ms; ie++) {
    for (ieg = ie * nGridPoints, nrm = 0.0, i = 0; i < nGridPoints; i++) { 
		nrm += sqr(psi[ieg+i]);
	}
	if (nrm != 0.0) { 
		nrm = 1.0 / sqrt(nrm * dr);
	}
	for (ieg = ie * nGridPoints, i = 0; i < nGridPoints; i++) { 
		psi[ieg+i] *= nrm;
	}
  }

  return;
}

/*****************************************************************************/
