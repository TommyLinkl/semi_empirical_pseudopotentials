#include "ar.h"

/*****************************************************************************/

double energy(zomplex *psi,zomplex *phi,double *potl,double *ksqr,lng_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  long     i;
  zomplex ene;

  memcpy(&phi[0],&psi[0],ist.nGridPoints*sizeof(phi[0]));
  hamiltonian(phi,psi,potl,ksqr,ist,par,planfw,planbw,fftwpsi);

  for (ene.re = ene.im = 0.0, i = 0; i < ist.nGridPoints; i++){
    ene.re += (psi[i].re * phi[i].re + psi[i].im * phi[i].im);
    /*ene.im += (psi[i].re * phi[i].im - psi[i].im * phi[i].re);*/
  }
  return (ene.re*par.dv);
}

/*****************************************************************************/

void energy_all(zomplex *psi,zomplex *phi,double *psi0,double *potl,double *ksqr,double *ene,lng_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,long ms)
{
  long i, ie, ieg;

  for (ie = 0; ie < ms; ie++){
    ieg = ie * ist.nGridPoints;
    for (i = 0; i < ist.nGridPoints; i++) {psi[i].re = psi0[ieg+i]; psi[i].im = 0.0;}
    memcpy(&phi[0],&psi[0],ist.nGridPoints*sizeof(phi[0]));
    hamiltonian(phi,psi,potl,ksqr,ist,par,planfw,planbw,fftwpsi);
    
    for (ene[ie] = 0.0, i = 0; i < ist.nGridPoints; i++)
      ene[ie] += (psi[i].re * phi[i].re + psi[i].im * phi[i].im);
    ene[ie] *= par.dv;
  }
  return;
}

/*****************************************************************************/

void get_energy_range(double *vx,double *vy,double *vz,double *ksqr,double *potl,par_st *par,lng_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  FILE *pf; long i;
  double ene, ene_old, norma;
  zomplex *psi, *phi;
  long idum = -21059783947;
  
  if ((psi   = (zomplex*)calloc(ist.nGridPoints,sizeof(zomplex)))==NULL)nerror("psi");
  if ((phi   = (zomplex*)calloc(ist.nGridPoints,sizeof(zomplex)))==NULL)nerror("phi");
  
  /*** initial wave function ***/
  init_psi(psi,ist,*par,&idum);
  
  /*** calculate the energy range ***/
  ene = (ene_old = 0.0) + 0.1;
  pf = fopen("ene-init.dat" , "w");
  for (i = 0; (fabs((ene-ene_old)/ene)>1.0e-4) & (i < 500); i++){
    hamiltonian(psi,phi,potl,ksqr,ist,*par,planfw,planbw,fftwpsi);
    norma = normalize(psi,par->dv,ist.nGridPoints);
    ene_old = ene;
    ene = energy(psi,phi,potl,ksqr,ist,*par,planfw,planbw,fftwpsi);
    fprintf (pf,"%ld %.16g %.16g %.16g\n",i,ene_old,ene,norma);
    fflush(pf);
  }
  fclose(pf);

  par->dE = 1.1 * (ene - par->Vmin);
  par->dE_1 = 4.0 / par->dE;

  printf ("dE = %g\n",par->dE);

  free(psi); free(phi);
  return;
}


/*****************************************************************************/

void calc_sigma_E(double *psitot,double *potl,double *ksqr,double *eval2,lng_st ist,par_st par)
{
  long jgrid, flags=0; double eval; zomplex *psi, *phi;
  fftw_plan_loc planfw, planbw; fftw_complex *fftwpsi;

  fftwpsi = fftw_malloc(sizeof (fftw_complex )*ist.nGridPoints);
  if ((psi = (zomplex*)calloc(ist.nGridPoints,sizeof(zomplex)))==NULL)nerror("psi");
  if ((phi = (zomplex*)calloc(ist.nGridPoints,sizeof(zomplex)))==NULL)nerror("psi");

  planfw = fftw_plan_dft_3d(ist.nz,ist.ny,ist.nx,fftwpsi,fftwpsi,FFTW_FORWARD,flags);
  planbw = fftw_plan_dft_3d(ist.nz,ist.ny,ist.nx,fftwpsi,fftwpsi,FFTW_BACKWARD,flags);
  
  for (jgrid = 0; jgrid < ist.nGridPoints; jgrid++) {
    psi[jgrid].re = psitot[jgrid];
    psi[jgrid].im = 0.0;
  }
  memcpy(&phi[0],&psi[0],ist.nGridPoints*sizeof(phi[0]));
  hamiltonian(phi,psi,potl,ksqr,ist,par,planfw,planbw,fftwpsi);

  for (eval = 0.0, jgrid = 0; jgrid < ist.nGridPoints; jgrid++)
    eval += (psitot[jgrid] * phi[jgrid].re);
  eval *= par.dv;
    
  memcpy(&psi[0],&phi[0],ist.nGridPoints*sizeof(psi[0]));
  hamiltonian(phi,psi,potl,ksqr,ist,par,planfw,planbw,fftwpsi);

  for (*eval2 = 0.0, jgrid = 0; jgrid < ist.nGridPoints; jgrid++)
    (*eval2) += (psitot[jgrid] * phi[jgrid].re);
  (*eval2) *= par.dv;
    
  (*eval2) -= sqr(eval);
  (*eval2) = sqrt(fabs((*eval2)));
  
  free(psi); free(phi);
  fftw_destroy_plan(planfw);
  fftw_destroy_plan(planbw);
  fftw_free(fftwpsi);
  return;
}

/*****************************************************************************/
