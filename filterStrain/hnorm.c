/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "fd.h"

/*****************************************************************************/

void hnorm(zomplex *psim1,zomplex *psin,double *potl,double *ksqr,par_st par,double zm1,long_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)

{
  long i;
  hamiltonian_norm(psim1,psin,potl,ksqr,ist,par,planfw,planbw,fftwpsi);

  for (i = 0; i < ist.ngrid; i++){
    /*** par.dE_1 = 4.0 / par.dE and therefore I don't multiply by 4 ***/
    psin[i].re = par.dE_1 * psin[i].re -
      (2.0 + zm1 + par.Vmin * par.dE_1) * psim1[i].re;
    psin[i].im = par.dE_1 * psin[i].im -
      (2.0 + zm1 + par.Vmin * par.dE_1) * psim1[i].im;
  }
  return;
}

/*****************************************************************************/

void hamiltonian_norm(zomplex *psi,zomplex *phi,double *potl,double *ksqr,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  long i;
  
  kinetic(phi,ksqr,planfw,planbw,fftwpsi,ist);
  for (i = 0; i < ist.ngrid; i++){
    phi[i].re += (potl[i] * psi[i].re);
    phi[i].im += (potl[i] * psi[i].im);
  }
  return;
}

/*****************************************************************************/
