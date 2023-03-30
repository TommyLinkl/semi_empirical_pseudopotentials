#include "fd.h"

/*****************************************************************************/

void generate_filter_states(double *psitot,double *eval,double *sige,double *evalbe,double *ksqr,double *potl,double *vx,double *vy,double *vz,lng_st *ist,par_st *par,fftw_plan_loc *planfw,fftw_plan_loc *planbw,fftw_complex *fftwpsi)
{
  FILE *pf; long idum;   zomplex *rho, *an; 
  long igrid, jmc, jms, tid, count1, count2, flags=0;  
  double wind = 1.0, rate1, rate2;
  double *el, *elfil, *zn, dE, dt, sum, tmp;
  
  idum = ist->seed;
  if (idum == 0) {
    Randomize(); idum = -random();
    printf("seed = %ld\n",idum); 
  }
  fflush(0);

  // Changed to include a larger energy window since the energy of the initial biexciton 
  // ranges from a minimum of 2.0*par->Egap (all carriers in homo/lumo states) to 2.0*par->Egap+3kbT
  par->Eamin = evalbe[0] + par->minInitE - par->DeltaE; // min electron energy when both excitons are in the ground state 
  par->Eamax = evalbe[ist->homoIndex] + par->maxInitE + par->DeltaE; // increased to include temperature effects 
  printf("Eamin = %g Eamax %g\n", par->Eamin, par->Eamax);
  
  // Changed to include a larger energy window since the energy of the initial biexciton
  // ranges from a minimum of 2.0*par->Egap (all carriers in homo/lumo states) to 2.0*par->Egap+3kbT 
  par->Eimin = evalbe[ist->lumoIndex] - (par->maxInitE + par->DeltaE); // decreased to include temperature effects 
  par->Eimax = evalbe[ist->numBandEdgeStates-1] - (par->minInitE - par->DeltaE);	// max hole energy when both excitons are in the ground state
  printf("Eimin = %g Eimax %g\n", par->Eimin, par->Eimax);
  
  if ((an = (zomplex*)calloc(ist->nc*ist->two,sizeof(zomplex)))==NULL)nerror("an");
  if ((zn = (double*)calloc(ist->nc,sizeof(double)))==NULL)nerror("zn");
  if ((el = (double*)calloc(ist->two,sizeof(double)))==NULL)nerror("el");
  if ((elfil = (double*)calloc(ist->two,sizeof(double)))==NULL)nerror("elfil");
  
  dt = sqr((double)(ist->nc) / (2.5 * par->dE));
  printf("nc = %ld dt = %g dE = %g\n",ist->nc,dt,par->dE); fflush(0);

  /*** generate the filter coeff ***/
  double dela = (par->Eamax - par->Eamin) / ((double)(ist->two) / 2.0 - 1.0);
  for (jms = 0; jms < ist->two / 2; jms++)
    el[jms] = par->Eamin + (double)(jms) * dela;

  double deli = (par->Eimax - par->Eimin) / ((double)(ist->two) / 2.0 - 1.0);
  for (jms = 0; jms < ist->two/2; jms++)
    el[jms+ist->two/2] = par->Eimin + (double)(jms) * deli;

  for (jms = 0; jms < ist->two; jms++)
    printf ("el = %g\n",el[jms]);
 
  coefficient(an,zn,ist->nc,ist->two,dt,par->dE,par->Vmin,el);

  /*** initial random states ***/
  if ((rho = (zomplex*)calloc(ist->ngrid,sizeof(zomplex)))==NULL)nerror("rho");
  for (jmc = 0; jmc < ist->nmc; jmc++){
    init_psi(rho,*ist,*par,&idum);
    for (jms = 0; jms < ist->two; jms++)
      for (igrid = 0; igrid < ist->ngrid; igrid++)
	psitot[jmc*ist->two*ist->ngrid+ist->ngrid*jms+igrid] = rho[igrid].re;
  }
  free(rho);
  
  /*** filter the states ***/
  omp_set_dynamic(0);
  omp_set_num_threads(ist->nthreads);
#pragma omp parallel for private(jmc,tid)
  for (jmc = 0; jmc < ist->nmc; jmc++){
    tid = omp_get_thread_num();	
    filtering(&psitot[jmc*ist->two*ist->ngrid],potl,ksqr,an,zn,el,*ist,*par,tid,jmc,&idum,planfw[tid],planbw[tid],&fftwpsi[tid*ist->ngrid]);
  } 

  /*** orthogonalize and normalize the filtered states using an svd routine ***/
  ist->ms = portho(psitot,par->dv,*ist);
  printf ("mstot ortho = %ld\n",ist->ms); fflush(0);
  normalize_all_omp(psitot,par->dv,ist->ms,ist->ngrid,ist->nthreads);

  /*** diagonalize the hamiltonian in the subspace spanned by the ***/
  Hmatreal(psitot,potl,ksqr,eval,*ist,*par,planfw[0],planbw[0],&fftwpsi[0]);
  normalize_all_omp(psitot,par->dv,ist->ms,ist->ngrid,ist->nthreads);

  /*** calculate the standard deviation of these states ***/
  pf = fopen("evalai.dat" , "w");
#pragma omp parallel for private(jms)
  for (jms = 0; jms < ist->ms; jms++)
    calc_sigma_E(&psitot[jms*ist->ngrid],potl,ksqr,&sige[jms],*ist,*par);

  for (jms = 0; jms < ist->ms; jms++) fprintf (pf,"%ld %g %g\n",jms,eval[jms],sige[jms]);
  fclose(pf);
  
  for (ist->itot = ist->atot = 0, jms = 0; jms < ist->ms; jms++) {
    if ((eval[jms] >= par->Eamin) && (eval[jms] <= par->Eamax) && (sige[jms] < par->deps)) ist->atot++;
    if ((eval[jms] >= par->Eimin) && (eval[jms] <= par->Eimax) && (sige[jms] < par->deps)) ist->itot++;
  }
  printf ("Number of filtered excited hole eigenstates in energy range = %ld\n", ist->itot);
  printf ("Number of filtered excited electron eigenstates in energy range = %ld\n", ist->atot);
  
  free(elfil); free(el); free(zn); free(an); 
  return;
}

/*****************************************************************************/
