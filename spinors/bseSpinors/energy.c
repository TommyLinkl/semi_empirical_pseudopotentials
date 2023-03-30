#include "fd.h"

/***************************************************************************/
// Calculates <psi|H|psi> (the energy) by first calculating H|psi> and then
// doing the inner product with <psi| and returning the resulting energy

double energy(zomplex *psi,zomplex *phi,double *potl,double *ksqr,lng_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  long     i;
  zomplex ene;

  memcpy(&phi[0],&psi[0],ist.ngrid*sizeof(phi[0]));
  hamiltonian(phi,psi,potl,ksqr,ist,par,planfw,planbw,fftwpsi);

  for (ene.re = ene.im = 0.0, i = 0; i < ist.ngrid; i++){
    ene.re += (psi[i].re * phi[i].re + psi[i].im * phi[i].im);
    /*ene.im += (psi[i].re * phi[i].im - psi[i].im * phi[i].re);*/
  }
  
  return (ene.re*par.dv);
}

/***************************************************************************/