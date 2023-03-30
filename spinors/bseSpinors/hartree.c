#include "fd.h"

/*****************************************************************************/
// Calculates the hatree potential via FFT for the integral over r 
// of rho(r)*potq(r-r') and stores result in poth

// TODO: Check to see if poth needs to be a zomplex - John
void hartree(zomplex *rho, zomplex *potq, double *poth, lng_st ist, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi)
{  
  long j;
  zomplex tmp;

  memcpy(&fftwpsi[0], &rho[0], ist.ngrid*sizeof(fftwpsi[0]));
  fftw_execute(planfw);

  for (j = 0; j < ist.ngrid; j++) {
    tmp.re = fftwpsi[j][0];
    tmp.im = fftwpsi[j][1];
    
    fftwpsi[j][0] = tmp.re * potq[j].re - tmp.im * potq[j].im;
    fftwpsi[j][1] = tmp.re * potq[j].im + tmp.im * potq[j].re;
  }

  fftw_execute(planbw);
  for (j = 0; j < ist.ngrid; j++) {
    poth[j] = fftwpsi[j][0];
  }

  return;
}

/*****************************************************************************/