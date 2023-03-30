/*****************************************************************************/
//
// This file deals with calculating Hartree potentials 
//
/*****************************************************************************/

#include "qp.h"

/*****************************************************************************/

void hartree(zomplex *rho, zomplex *potq, double *poth, long nGridPoints, 
              fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi) {
  long iGrid;
  zomplex tmp;
  
  // FFT to q-space
  memcpy(&fftwpsi[0], &rho[0], nGridPoints*sizeof(fftwpsi[0]));
  fftw_execute(planfw);

  // Do integration in q-space 
  for (iGrid = 0; iGrid < nGridPoints; iGrid++) {
    tmp.re = fftwpsi[iGrid][0];
    tmp.im = fftwpsi[iGrid][1];
    
    fftwpsi[iGrid][0] = tmp.re * potq[iGrid].re - tmp.im * potq[iGrid].im;
    fftwpsi[iGrid][1] = tmp.re * potq[iGrid].im + tmp.im * potq[iGrid].re;
  }

  // FFT to r-space
  fftw_execute(planbw);

  for (iGrid = 0; iGrid < nGridPoints; iGrid++) { 
    poth[iGrid] = fftwpsi[iGrid][0];
  }
  
  return;
}

/*****************************************************************************/
