#include "fd.h"

/****************************************************************************/
// Calculates H|phi> = (T+V_local)|phi>=T|phi>+V_local|phi>
// and stores the result in phi

void hamiltonian(zomplex *phi,zomplex *psi,double *potl,double *ksqr,lng_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  long i;

  memcpy(&psi[0],&phi[0],ist.ngrid*sizeof(psi[0]));
  kinetic(phi,ksqr,planfw,planbw,fftwpsi,ist);
  for (i = 0; i < ist.ngrid; i++){
    phi[i].re += (potl[i] * psi[i].re);
    phi[i].im += (potl[i] * psi[i].im);
  }

  return;
}

/****************************************************************************/
// Calculates T|psi> via FFT and stores result in psi 


void kinetic(zomplex *psi,double *ksqr,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,lng_st ist)
{
  long i;

  memcpy(&fftwpsi[0],&psi[0],ist.ngrid*sizeof(fftwpsi[0]));
  fftw_execute(planfw);

  for (i = 0; i < ist.ngrid; i++){
    fftwpsi[i][0] *= ksqr[i];
    fftwpsi[i][1] *= ksqr[i];
  }
  fftw_execute(planbw);

  memcpy(&psi[0],&fftwpsi[0],ist.ngrid*sizeof(psi[0]));
  
  return;
}

/****************************************************************************/