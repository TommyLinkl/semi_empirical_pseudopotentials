/****************************************************************************/

#include "fd.h"

/****************************************************************************/

double energy(zomplex *psi,zomplex *phi,double *potl,double *ksqr,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
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

double energy_norm(zomplex *psi,zomplex *phi,double *potl,double *ksqr,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  long     i;
  double  sum;
  zomplex ene;

  memcpy(&phi[0],&psi[0],ist.ngrid*sizeof(phi[0]));
  hamiltonian(phi,psi,potl,ksqr,ist,par,planfw,planbw,fftwpsi);

  for (ene.re = ene.im = 0.0, i = 0; i < ist.ngrid; i++){
    ene.re += (psi[i].re * phi[i].re + psi[i].im * phi[i].im);
    sum += (psi[i].re * psi[i].re + psi[i].im * psi[i].im);
    /*ene.im += (psi[i].re * phi[i].im - psi[i].im * phi[i].re);*/
  }
  return (ene.re/sum);
}

/***************************************************************************/

void energy_all(double *psi0,double *potl,double *ksqr,double *ene,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,long ms)
{
  long jgrid, jms, jmsg;
  zomplex *psi, *phi;

  if ((psi = (zomplex*)calloc(ist.ngrid,sizeof(zomplex)))==NULL)nerror("psi");
  if ((phi = (zomplex*)calloc(ist.ngrid,sizeof(zomplex)))==NULL)nerror("phi");
    
  for (jms = 0; jms < ms; jms++){
    jmsg = jms * ist.ngrid;
    for (jgrid = 0; jgrid < ist.ngrid; jgrid++) {psi[jgrid].re = psi0[jmsg+jgrid]; psi[jgrid].im = 0.0;}
    memcpy(&phi[0],&psi[0],ist.ngrid*sizeof(phi[0]));
    hamiltonian(phi,psi,potl,ksqr,ist,par,planfw,planbw,fftwpsi);
    
    for (ene[jms] = 0.0, jgrid = 0; jgrid < ist.ngrid; jgrid++)
      ene[jms] += (psi[jgrid].re * phi[jgrid].re + psi[jgrid].im * phi[jgrid].im);
    ene[jms] *= par.dv;
  }

  free(psi); free(phi);
  return;
}

/***************************************************************************/

void get_energy_range(double *vx,double *vy,double *vz,double *ksqr,double *potl,par_st *par,long_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  FILE *pf; long i;
  double ene, ene_old, norma;
  zomplex *psi, *phi;
  long idum = -21059783947;
  
  if ((psi = (zomplex*)calloc(ist.ngrid,sizeof(zomplex)))==NULL)nerror("psi");
  if ((phi = (zomplex*)calloc(ist.ngrid,sizeof(zomplex)))==NULL)nerror("phi");
  
  /*** initial wave function ***/
  init_psi(psi,vx,vy,vz,ist,*par,&idum);
  
  /*** calculate the energy range ***/
  ene = (ene_old = 0.0) + 0.1;
  pf = fopen("ene-init.dat" , "w");
  for (i = 0; (fabs((ene-ene_old)/ene)>1.0e-4) & (i < 500); i++){
    hamiltonian(psi,phi,potl,ksqr,ist,*par,planfw,planbw,fftwpsi);
    //norma = normalize(psi,par->dv,ist.ngrid);
    ene_old = ene;
    ene = energy(psi,phi,potl,ksqr,ist,*par,planfw,planbw,fftwpsi);
    fprintf (pf,"%ld %.16g %.16g %.16g\n",i,ene_old,ene,norma);
    fflush(pf);
  }
  fclose(pf);

  /*par->dE = 1.1 * (ene - par->Vmin);
  par->dE_1 = 4.0 / par->dE;

  printf ("dE = %g\n",par->dE);*/

  free(psi); free(phi);
  return;
}


/****************************************************************************************/

void calc_sigma_E(zomplex *psi,zomplex *phi,double *psitot,double *potl,double *ksqr,double *eval2,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  long jgrid, ims;
  double eval;
  
  for (ims = 0; ims < ist.ms2; ims++){
    for (jgrid = 0; jgrid < ist.ngrid; jgrid++) {
      psi[jgrid].re = psitot[ims*ist.ngrid+jgrid];
      psi[jgrid].im = 0.0;
    }
    memcpy(&phi[0],&psi[0],ist.ngrid*sizeof(phi[0]));
    hamiltonian(phi,psi,potl,ksqr,ist,par,planfw,planbw,fftwpsi);

    for (eval = 0.0, jgrid = 0; jgrid < ist.ngrid; jgrid++)
      eval += (psitot[ims*ist.ngrid+jgrid] * phi[jgrid].re);
    eval *= par.dv;
    /*eval2[ims] = eval;*/
    
    memcpy(&psi[0],&phi[0],ist.ngrid*sizeof(psi[0]));
    hamiltonian(phi,psi,potl,ksqr,ist,par,planfw,planbw,fftwpsi);

    for (eval2[ims] = 0.0, jgrid = 0; jgrid < ist.ngrid; jgrid++)
      eval2[ims] += (psitot[ims*ist.ngrid+jgrid] * phi[jgrid].re);
    eval2[ims] *= par.dv;
    
    eval2[ims] -= sqr(eval);
    eval2[ims] = sqrt(fabs(eval2[ims]));
  }
  return;
}

/****************************************************************************/