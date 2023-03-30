/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "fd.h"

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

double normd(double *psi, double dr,long ngrid,long nthreads)
{
	long i;
	double norm = 0.0;

	omp_set_dynamic(0);
	omp_set_num_threads(nthreads);
#pragma omp parallel for reduction(+:norm)
	for (i = 0; i < ngrid; i++) {
		norm += (sqr(psi[i]));
	}
	return (sqrt(norm * dr));
}

/*****************************************************************************/

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

void normalize_all(double *psi,double dr,long ms,long ngrid,long nthreads)
{
	long i, ie, ieg;
	double nrm;

	omp_set_dynamic(0);
	omp_set_num_threads(nthreads);
#pragma omp parallel for private(ie,nrm,ieg,i)
	for (ie = 0; ie < ms; ie++){
		for (ieg = ie * ngrid, nrm = 0.0, i = 0; i < ngrid; i++) nrm += sqr(psi[ieg+i]);
		if (nrm != 0.0) {
			nrm = 1.0 / sqrt(nrm * dr);
		}
		for (ieg = ie * ngrid, i = 0; i < ngrid; i++) psi[ieg+i] *= nrm;
	}
	return;
}
