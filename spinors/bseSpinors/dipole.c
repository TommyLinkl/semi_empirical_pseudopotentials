#include "fd.h"

/***************************************************************************************/
// Calculates the dipole matrix elements between the quasielectron and hole spinor 
// states and stores them in mux - the square of which is used to calculate the oscillator
// strength here in the non-interacting (in OS0.dat) and later the matrix elements are used 
// in a coherent sum to calculate the optical absorption spectra (OS.dat)

void dipole(double *vx,double *vy,double *vz,zomplex *psi,zomplex *mux,zomplex *muy,zomplex *muz,double *eval,lng_st ist,par_st par)
{
  FILE *pf1, *pf2, *pf3, *pf; 
  long a, i, jx, jy, jz, jgrid, jyz, pairIndex; 
  double dz, dx, dy, r, ev, os;
  zomplex sumX, sumY, sumZ, tmp;
  
  pf = fopen("OS0.dat" , "w"); pf1 = fopen("ux.dat", "w"); pf2 = fopen("uy.dat", "w"); pf3 = fopen("uz.dat", "w");
  for (i = 0; i < ist.ms2; i++) mux[i].re = mux[i].im = muy[i].re = muy[i].im = muz[i].re = muz[i].im = 0.0;

  for (i = 0; i < ist.totalhomo; i++, fprintf(pf1,"\n"),fprintf(pf2,"\n"),fprintf(pf3,"\n")) {
    for (a = ist.nlumo; a < ist.nlumo+ist.totallumo; a++) {
      for (jz = 0; jz < ist.nz; jz++) {
        sumX.re = sumX.im = sumY.re = sumY.im = sumZ.re = sumZ.im = 0.0;
        dz = vz[jz];
        for (jy = 0; jy < ist.ny; jy++) {
          dy = vy[jy];
          jyz = ist.nx * (ist.ny * jz + jy);
          for (jx = 0; jx < ist.nx; jx++) {
            dx = vx[jx];
            jgrid = jyz + jx;
      	    tmp.re = ((psi[i*ist.nspinngrid+jgrid].re * psi[a*ist.nspinngrid+jgrid].re 
              + psi[i*ist.nspinngrid+jgrid].im * psi[a*ist.nspinngrid+jgrid].im
              + psi[i*ist.nspinngrid+ist.ngrid+jgrid].re * psi[a*ist.nspinngrid+ist.ngrid+jgrid].re
              + psi[i*ist.nspinngrid+ist.ngrid+jgrid].im * psi[a*ist.nspinngrid+ist.ngrid+jgrid].im) * par.dv);
            tmp.im = ((psi[i*ist.nspinngrid+jgrid].re * psi[a*ist.nspinngrid+jgrid].im 
              - psi[i*ist.nspinngrid+jgrid].im * psi[a*ist.nspinngrid+jgrid].re
              + psi[i*ist.nspinngrid+ist.ngrid+jgrid].re * psi[a*ist.nspinngrid+ist.ngrid+jgrid].im 
              - psi[i*ist.nspinngrid+ist.ngrid+jgrid].im * psi[a*ist.nspinngrid+ist.ngrid+jgrid].re) * par.dv);
            sumX.re += dx*tmp.re;
      	    sumX.im += dx*tmp.im;
            sumY.re += dy*tmp.re;
            sumY.im += dy*tmp.im;
            sumZ.re += dz*tmp.re;
      	    sumZ.im += dz*tmp.im;
          }
        }
      }
      pairIndex = i*ist.totallumo+(a-ist.nlumo);
      mux[pairIndex].re =  sumX.re; mux[pairIndex].im =  sumX.im;
      muy[pairIndex].re =  sumY.re; muy[pairIndex].im =  sumY.im;
      muz[pairIndex].re =  sumZ.re; muz[pairIndex].im =  sumZ.im;
      
      fprintf(pf1, "%ld %ld %.8f %.8f\n", i, a, sumX.re, sumX.im);
      fprintf(pf2, "%ld %ld %.8f %.8f\n", i, a, sumY.re, sumY.im);
      fprintf(pf3, "%ld %ld %.8f %.8f\n", i, a, sumZ.re, sumZ.im);

      os=(sqr(mux[pairIndex].re) + sqr(mux[pairIndex].im) + sqr(muy[pairIndex].re) + sqr(muy[pairIndex].im)  
        + sqr(muz[pairIndex].re) + sqr(muz[pairIndex].im));
      ev = eval[a] - eval[i];
      fprintf(pf,"%.8f %.12f %.8f %.8f %.8f %.8f\n",sqrt(os),ev,(4.0/3.0)*ev*os,
	       sqr(mux[pairIndex].re) + sqr(mux[pairIndex].im),
	       sqr(muy[pairIndex].re) + sqr(muy[pairIndex].im),
	       sqr(muz[pairIndex].re) + sqr(muz[pairIndex].im));
    }
  }

  fclose(pf); fclose(pf1); fclose(pf2); fclose(pf3);

  return;
}

/***************************************************************************************/