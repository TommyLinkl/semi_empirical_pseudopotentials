#include "fd.h"

/*****************************************************************************/

void hnorm(zomplex *psim1,zomplex *psin,double *potl,double *ksqr,par_st par,nlc_st *nlc,long *nl,double *Elkb,double zm1,long_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  long i;
  hamiltonian(psin,psim1,potl,ksqr,ist,par,nlc,nl,Elkb,planfw,planbw,fftwpsi);

  for (i = 0; i < ist.nspinngrid; i++){
    /*** par.dE_1 = 4.0 / par.dE and therefore I don't multiply by 4 ***/
    psin[i].re = par.dE_1 * psin[i].re - (2.0 + zm1 + par.Vmin * par.dE_1) * psim1[i].re;
    psin[i].im = par.dE_1 * psin[i].im - (2.0 + zm1 + par.Vmin * par.dE_1) * psim1[i].im;
  }
  
  return;
}

/*****************************************************************************/

void calcAllSpinUpDownRatios(zomplex *psitot, long nGrid, long nStates) {
	FILE *pf;
	long iGrid, iState, iStateGrid;
	double stateNorm, spinUpPercentage, spinDownPercentage;

	pf = fopen("spinPercentages.dat", "w");
	for (iState = 0; iState < nStates; iState++) {
		iStateGrid = iState*nGrid;

		spinUpPercentage = 0.0;
		for (iGrid = 0; iGrid < (nGrid/2); iGrid++) {
			spinUpPercentage += (sqr(psitot[iStateGrid + iGrid].re) + sqr(psitot[iStateGrid + iGrid].im));
		}

		spinDownPercentage = 0.0;
		for (iGrid = (nGrid/2); iGrid < nGrid; iGrid++) {
			spinDownPercentage += (sqr(psitot[iStateGrid + iGrid].re) + sqr(psitot[iStateGrid + iGrid].im));
		}

		stateNorm = spinUpPercentage+spinDownPercentage;
		spinUpPercentage /= stateNorm;
		spinDownPercentage /= stateNorm;
		
		fprintf(pf, "%ld %lg %lg\n", iState, spinUpPercentage, spinDownPercentage);
		fflush(pf);
	}
	fclose(pf);

	return;
}

/*****************************************************************************/
