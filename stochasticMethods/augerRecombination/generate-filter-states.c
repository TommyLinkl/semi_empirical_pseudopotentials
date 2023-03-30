/*****************************************************************************/

//

/*****************************************************************************/

#include "ar.h"

/*****************************************************************************/

void generate_filter_states(double *psitot, double *eval, double *sige, double *evalbe, double *ksqr, double *potl,
                            double *vx, double *vy, double *vz, lng_st *ist, par_st *par, 
                            fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwpsi)
{
  FILE *pf; 
  long igrid, jmc, jms, tid, idum = ist->seed;;  
  double dela, deli;
  double *el, *zn, dt;
  zomplex *rho, *an;

  // Write beginning of function
  writeSeparation(stdout);
  writeCurrentTime(stdout);
  fprintf(stdout, "Beginnning the calculation of the hot hole and electron final states\n\n");
  fflush(stdout);

  // Dynamically allocate memory
  if ((an = (zomplex *) calloc(ist->nNewtonIntSteps*ist->nStatesPerFilter, sizeof(zomplex))) == NULL) nerror("an");
  if ((zn = (double *) calloc(ist->nNewtonIntSteps, sizeof(double))) == NULL) nerror("zn");
  if ((el = (double *) calloc(ist->nStatesPerFilter, sizeof(double))) == NULL) nerror("el");

  // Determine the maximum and minimum allowed energy for the hot electron (a) and hole (a) states
  par->Eamin = evalbe[0] + par->minInitE - par->maxDeltaE;
  par->Eamax = evalbe[ist->homoIndex] + par->maxInitE + par->maxDeltaE;
  par->Eimin = evalbe[ist->lumoIndex] - (par->maxInitE + par->maxDeltaE);  
  par->Eimax = evalbe[ist->nHolesPlusElecs-1] - (par->minInitE - par->maxDeltaE);

  // Determine the energy to center the filtered states
  dela = (par->Eamax - par->Eamin) / ((double)(ist->nStatesPerFilter) / 2.0 - 1.0);
  for (jms = 0; jms < ist->nStatesPerFilter / 2; jms++) {
    el[jms] = par->Eamin + (double)(jms) * dela;
  }
  deli = (par->Eimax - par->Eimin) / ((double)(ist->nStatesPerFilter) / 2.0 - 1.0);
  for (jms = 0; jms < ist->nStatesPerFilter/2; jms++) {
    el[jms+ist->nStatesPerFilter/2] = par->Eimin + (double)(jms) * deli;
  }

  // Print some filter related information to stdout
  writeFilterInfo(el, *ist, *par, stdout);
  fflush(stdout);
  
  // Generate the filter coefficients
  dt = sqr((double)(ist->nNewtonIntSteps) / (2.5 * par->dE));
  coefficient(an, zn, ist->nNewtonIntSteps, ist->nStatesPerFilter, dt, par->dE, par->Vmin, el);

  // Generate initial random states 
  if ((rho = (zomplex *) calloc(ist->nGridPoints, sizeof(zomplex))) == NULL) nerror("rho");
  for (jmc = 0; jmc < ist->nFilterCycles; jmc++) {
    init_psi(rho, *ist, *par, &idum);
    for (jms = 0; jms < ist->nStatesPerFilter; jms++) {
      for (igrid = 0; igrid < ist->nGridPoints; igrid++) {
        psitot[jmc*ist->nStatesPerFilter*ist->nGridPoints + ist->nGridPoints*jms+igrid] = rho[igrid].re;
      }
    }
  }
  free(rho);
  
  // Filter the states
  printf("Maximum number of filtered hot carrier states that can be calculated = %ld\n", ist->nFilteredStates);
  fflush(stdout);
  omp_set_dynamic(0);
  omp_set_num_threads(ist->nThreads);
#pragma omp parallel for private(jmc,tid)
  for (jmc = 0; jmc < ist->nFilterCycles; jmc++) {
    tid = omp_get_thread_num();	
    filtering(&psitot[jmc*ist->nStatesPerFilter*ist->nGridPoints], potl, ksqr, an, zn, el, *ist, *par,
              tid, jmc, &idum, planfw[tid], planbw[tid], &fftwpsi[tid*ist->nGridPoints]);
  } 

  // Orthogonalize and normalize the filtered states using an svd routine 
  ist->nFilteredStates = portho(psitot, par->dv, *ist);
  printf("Total number of filtered hot carrier states, nFilteredStates = %ld\n", ist->nFilteredStates); 
  fflush(stdout);
  normalize_all_omp(psitot, par->dv, ist->nFilteredStates, ist->nGridPoints, ist->nThreads);

  // Diagonalize the hamiltonian in the subspace spanned by the psitot vectors
  Hmatreal(psitot, potl, ksqr, eval, *ist, *par, planfw[0], planbw[0], &fftwpsi[0]);
  normalize_all_omp(psitot, par->dv, ist->nFilteredStates, ist->nGridPoints, ist->nThreads);

  // Calculate the standard deviation of these states 
#pragma omp parallel for private(jms)
  for (jms = 0; jms < ist->nFilteredStates; jms++) {
    calc_sigma_E(&psitot[jms*ist->nGridPoints], potl, ksqr, &sige[jms], *ist, *par);
  }

  // Write eval.dat file containing the energies and sigma values of the hot hole (i) and electron (a) states
  pf = fopen("evalai.dat", "w");
  for (jms = 0; jms < ist->nFilteredStates; jms++) {
    fprintf(pf, "%ld % .10f %g\n", jms, eval[jms], sige[jms]);
  }
  fclose(pf);
  
  // Determine the number of electron and hole eigenstates that can be part of a final excitonic state
  ist->itot = ist->atot = 0;
  for (jms = 0; jms < ist->nFilteredStates; jms++) {
    if ((eval[jms] >= par->Eamin) && (eval[jms] <= par->Eamax) && (sige[jms] < par->sigmaCutoff)) {
      ist->atot++;
    }
    else if ((eval[jms] >= par->Eimin) && (eval[jms] <= par->Eimax) && (sige[jms] < par->sigmaCutoff)) {
      ist->itot++;
    }
  }
  printf("The number of filtered excited hole eigenstates in energy range, itot = %ld\n", ist->itot);
  printf("The number of filtered excited electron eigenstates in energy range, atot = %ld\n", ist->atot);
  // TODO: replace itot and atot with nHotHoles and nHotElecs
  ist->nHotHoles = ist->itot;
  ist->nHotElecs = ist->atot;
  ist->nHotNonIntStates = ist->nFilteredStates;

  // Write ending of function
  fprintf(stdout, "\nFinished the calculation of the hot hole and electron final states\n");
  writeCurrentTime(stdout);
  fflush(stdout);

  // Free dynamically allocated memory
  free(el); free(zn); free(an); 
  
  return;
}

/*****************************************************************************/
