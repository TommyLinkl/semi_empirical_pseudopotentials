/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "fd.h"

/*****************************************************************************/

double energy(zomplex *psi,zomplex *phi,double *potl,double *ksqr,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  long     i;
  zomplex ene;

  memcpy(&phi[0], &psi[0], ist.ngrid*sizeof(phi[0]));
  hamiltonian(phi, psi, potl, ksqr, ist, par, planfw, planbw, fftwpsi);

  for (ene.re = ene.im = 0.0, i = 0; i < ist.ngrid; i++){
    ene.re += (psi[i].re * phi[i].re + psi[i].im * phi[i].im);
    /*ene.im += (psi[i].re * phi[i].im - psi[i].im * phi[i].re);*/
  }
  return (ene.re*par.dv);
}

/*****************************************************************************/

void energy_all(zomplex *psi,zomplex *phi,double *psi0,double *potl,double *ksqr,double *ene,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  long i, ie, ieg;

  for (ie = 0; ie < ist.ms; ie++){
    ieg = ie*ist.ngrid;
    for (i = 0; i < ist.ngrid; i++) {
      psi[i].re = psi0[ieg+i]; 
      psi[i].im = 0.0;
    }

    memcpy(&phi[0], &psi[0], ist.ngrid*sizeof(phi[0]));
    hamiltonian(phi, psi, potl, ksqr, ist, par, planfw, planbw, fftwpsi);
    
    for (ene[ie] = 0.0, i = 0; i < ist.ngrid; i++) {
      ene[ie] += (psi[i].re * phi[i].re + psi[i].im * phi[i].im);
    }
    ene[ie] *= par.dv;
  }

  return;
}

/*****************************************************************************/

void get_energy_range(zomplex *psi,zomplex *phi,double *potl,double *vx,double *vy,double *vz,double *ksqr,par_st *par,long_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  FILE *pf;
  long i; 
  long idum = -874917403;
  double ene, ene_old, norma;
  
  /*** initial wave function ***/
  init_psi(psi,ist,*par,&idum);
  
  /*** calculate the energy range ***/
  ene = (ene_old = 0.0) + 0.1;
  pf = fopen("ene-init.dat" , "w");
  for (i = 0; (fabs((ene-ene_old)/ene)>1.0e-4) && (i < 500); i++) {
    hamiltonian(psi, phi, potl, ksqr, ist, *par, planfw, planbw, fftwpsi);
    norma = normalize(psi, par->dv, ist.ngrid, ist.nthreads);
    ene_old = ene;
    ene = energy(psi, phi, potl, ksqr, ist, *par, planfw, planbw, fftwpsi);
    fprintf(pf, "%ld %.16g %.16g %.16g\n", i, ene_old, ene, norma);
    fflush(pf);
  }
  fclose(pf);

  par->dE = 1.1*(ene - par->Vmin);
  par->dE_1 = (4.0 / par->dE);

  printf("dE = %g\n", par->dE);
  fflush(stdout);

  return;
}


/*****************************************************************************/

void calc_sigma_E(zomplex *psi,zomplex *phi,double *psitot,double *potl,double *ksqr,double *eval2,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  long jgrid, ims;
  double eval;
  
  for (ims = 0; ims < ist.mstot; ims++){
    for (jgrid = 0; jgrid < ist.ngrid; jgrid++) {
      psi[jgrid].re = psitot[ims*ist.ngrid+jgrid];
      psi[jgrid].im = 0.0;
    }
    memcpy(&phi[0],&psi[0],ist.ngrid*sizeof(phi[0]));
    hamiltonian(phi,psi,potl,ksqr,ist,par,planfw,planbw,fftwpsi);

    for (eval = 0.0, jgrid = 0; jgrid < ist.ngrid; jgrid++) {
      eval += (psitot[ims*ist.ngrid+jgrid] * phi[jgrid].re);
    }
    eval *= par.dv;
    /*eval2[ims] = eval;*/
    
    memcpy(&psi[0],&phi[0],ist.ngrid*sizeof(psi[0]));
    hamiltonian(phi,psi,potl,ksqr,ist,par,planfw,planbw,fftwpsi);

    for (eval2[ims] = 0.0, jgrid = 0; jgrid < ist.ngrid; jgrid++) {
      eval2[ims] += (psitot[ims*ist.ngrid+jgrid] * phi[jgrid].re);
    }
    eval2[ims] *= par.dv;
    
    eval2[ims] -= sqr(eval);
    eval2[ims] = sqrt(fabs(eval2[ims]));
  }
  return;
}

/*****************************************************************************/