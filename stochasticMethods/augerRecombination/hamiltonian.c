#include "ar.h"

/****************************************************************************/

void hamiltonian(zomplex *phi,zomplex *psi,double *potl,double *ksqr,lng_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  long i;

  //memcpy(&fftwpsi[0],&phi[0],ist.nGridPoints*sizeof(fftwpsi[0]));
  memcpy(&psi[0],&phi[0],ist.nGridPoints*sizeof(psi[0]));
  kinetic(phi,ksqr,planfw,planbw,fftwpsi,ist);
  for (i = 0; i < ist.nGridPoints; i++){
    phi[i].re += (potl[i] * psi[i].re);
    phi[i].im += (potl[i] * psi[i].im);
  }
  return;
}

/****************************************************************************/

void kinetic(zomplex *psi,double *ksqr,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,lng_st ist)
{
  long j, flags=0;

  memcpy(&fftwpsi[0],&psi[0],ist.nGridPoints*sizeof(fftwpsi[0]));
  fftw_execute(planfw);
  //fftwnd_one(planfw,psi,fftwpsi);

  for (j = 0; j < ist.nGridPoints; j++){
    fftwpsi[j][0] *= ksqr[j];
    fftwpsi[j][1] *= ksqr[j];
  }
  fftw_execute(planbw);
  //fftwnd_one(planbw,fftwpsi,psi);

  memcpy(&psi[0],&fftwpsi[0],ist.nGridPoints*sizeof(psi[0]));
  return;
}

/****************************************************************************/