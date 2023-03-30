/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "ar.h"

/*****************************************************************************/

int main(long argc, char *argv[]) {
  FILE *pf; 
  par_st par;  
  lng_st ist;
  long i, flags = 0; 
  double *ksqr, *vx, *vy, *vz, *poth, *potl, *rx, *ry, *rz;
  double *psibe, *evalbe, *sigebe;
  double *psiai, *evalai, *sigeai;
  double *vijck, *vabck;
  double *Cbs, *Ebs; 
  zomplex *potq, *psi, *phi;
  fftw_plan_loc *planfw, *planbw; fftw_complex *fftwpsi;

  /*************************************************************************/

  // Write time the program began and read input
  writeCurrentTime(stdout);
  writeSeparation(stdout);
  init_size(&par, &ist);
    
  /*************************************************************************/

  // Dynamically allocate memory
  if ((fftwpsi = fftw_malloc(sizeof (fftw_complex )*ist.nGridPoints*ist.nThreads)) == NULL) nerror("fftwpsi");
  if ((potq  = (zomplex *) calloc(ist.nGridPoints, sizeof(zomplex))) == NULL) nerror("potq");
  if ((poth = (double *) calloc(ist.nGridPoints*ist.nThreads, sizeof(double))) == NULL) nerror("poth");
  if ((potl = (double *) calloc(ist.nGridPoints, sizeof(double))) == NULL) nerror("potl");
  if ((ksqr = (double *) calloc(ist.nGridPoints, sizeof(double))) == NULL) nerror("ksqr");
  if ((vx = (double *) calloc(ist.nx, sizeof(double))) == NULL) nerror("vx");
  if ((vy = (double *) calloc(ist.ny, sizeof(double))) == NULL) nerror("vy");
  if ((vz = (double *) calloc(ist.nz, sizeof(double))) == NULL) nerror("vz");
  if ((rx = (double *) calloc(ist.nAtoms, sizeof(double))) == NULL) nerror("rx");
  if ((ry = (double *) calloc(ist.nAtoms, sizeof(double))) == NULL) nerror("ry");
  if ((rz = (double *) calloc(ist.nAtoms, sizeof(double))) == NULL) nerror("rz");
  if ((psi  = (zomplex *) calloc(ist.nGridPoints, sizeof(zomplex))) == NULL) nerror("psi");
  if ((phi  = (zomplex *) calloc(ist.nGridPoints, sizeof(zomplex))) == NULL) nerror("phi");
  if ((psibe = (double *) calloc(ist.nHolesPlusElecs*ist.nGridPoints, sizeof(double))) == NULL) nerror("psibe");
  if ((evalbe = (double *) calloc(ist.nHolesPlusElecs, sizeof(double))) == NULL) nerror("evalbe");
  if ((sigebe = (double *) calloc(ist.nHolesPlusElecs, sizeof(double))) == NULL) nerror("sigebe");
  if (! ist.readInHotQPs) { // new hot hole and elec states will be calculated
    if ((psiai = (double *) calloc(ist.nFilteredStates*ist.nGridPoints, sizeof(double))) == NULL) nerror("psiai");
    if ((evalai = (double *) calloc(ist.nFilteredStates, sizeof(double))) == NULL) nerror("evalai");
    if ((sigeai = (double *) calloc(ist.nFilteredStates, sizeof(double))) == NULL) nerror("sigeai");    
  }
  else { // hot hole and elec states will be read in
    if ((psiai = (double *) calloc(ist.nHotNonIntStates*ist.nGridPoints, sizeof(double))) == NULL) nerror("psiai");
    if ((evalai = (double *) calloc(ist.nHotNonIntStates, sizeof(double))) == NULL) nerror("evalai");
    if ((sigeai = (double *) calloc(ist.nHotNonIntStates, sizeof(double))) == NULL) nerror("sigeai");    
  }

  if (ist.intAR) {
    if ((Cbs = (double *) calloc(ist.nNonIntExcitons*ist.nIntExcitons, sizeof(double))) == NULL) nerror("Cbs");
    if ((Ebs = (double *) calloc(ist.nIntExcitons, sizeof(double))) == NULL) nerror("Ebs");
  }
  
  /**************************************************************************/
  if (ist.intAR) {
    readBetheSalpeterResults(Cbs, Ebs, ist);
  }
  init(vx, vy, vz, ksqr, potl, rx, ry, rz, &par, &ist);

  /**************************************************************************/

  // Initialize the FFT and allocate memory
  fftw_plan_with_nthreads(ist.nThreads);
  planfw = (fftw_plan_loc *) calloc(ist.nThreads, sizeof(fftw_plan_loc));
  planbw = (fftw_plan_loc *) calloc(ist.nThreads, sizeof(fftw_plan_loc));
  for (i = 0; i < ist.nThreads; i++){ 
    planfw[i] = fftw_plan_dft_3d(ist.nz, ist.ny, ist.nx, &fftwpsi[i*ist.nGridPoints], &fftwpsi[i*ist.nGridPoints], FFTW_FORWARD,  flags);
    planbw[i] = fftw_plan_dft_3d(ist.nz, ist.ny, ist.nx, &fftwpsi[i*ist.nGridPoints], &fftwpsi[i*ist.nGridPoints], FFTW_BACKWARD, flags);
  }
  
  /**************************************************************************/

  // Calculate the 
  init_pot(vx, vy, vz, potq, par, ist, planfw[0], planbw[0], &fftwpsi[0]);

  /**************************************************************************/
  // 
  gauss_test(vx, vy, vz, potq, &poth[0], par, ist, planfw[0], planbw[0], &fftwpsi[0]);

  /**************************************************************************/
  //
  get_energy_range(vx, vy, vz, ksqr, potl, &par, ist, planfw[0], planbw[0], &fftwpsi[0]);
  
  /**************************************************************************/

  // Store the wavefunctions and energies of the noninteracting hole and elec eigenstates
  storeNonIntEigenstates(psibe, evalbe, sigebe, par, ist);
  ist.homoIndex = ist.nHoles-1;  ist.lumoIndex = ist.nHoles;
  if (ist.intAR) {
    ist.nBiexcitons = calcNumIntBiexcStates(Ebs, par.maxInitE, ist.nIntExcitons, ist.nTrappedStates);
  }
  else if (ist.nonIntAR) {
    ist.nBiexcitons = calcNumNonIntBiexcStates(evalbe, par.maxInitE, ist.nHoles, ist.nElecs, ist.nTrappedStates);  
  }

  /**************************************************************************/

  // Obtain the hot hole and electron states
  if (ist.readInHotQPs) {
    // Read in the hot hole (i) and hot electron (a) states 
    ist.nHotNonIntStates = readNonIntEigenstates(psiai, evalai, sigeai, "psiai.par", "evalai.par", par.sigmaCutoff, ist.nGridPoints);
  }
  else {
    // Calculate the hot hole (i) and hot electron (a) states
    generate_filter_states(psiai, evalai, sigeai, evalbe, ksqr, potl, vx, vy, vz, &ist, &par, planfw, planbw, fftwpsi);

  }

  // Remove states from psiai, evalai and sigeai that are not needed
  //ist.nHotNonIntStates = remUnnecessaryHotNonIntStates(psiai, evalai, sigeai, par, ist);
  
  // Realloc memory for the smaller psiai, evalai and sigeai array
  if (! ist.readInHotQPs) {
    // Remove states from psiai, evalai and sigeai that are not needed
    ist.nHotNonIntStates = remUnnecessaryHotNonIntStates(psiai, evalai, sigeai, par, ist);
    psiai = realloc(psiai, (ist.nHotNonIntStates+1)*ist.nGridPoints*sizeof(psiai[0]));
    evalai = realloc(evalai, (ist.nHotNonIntStates+1)*sizeof(evalai[0]));
    sigeai = realloc(sigeai, (ist.nHotNonIntStates+1)*sizeof(sigeai[0]));
  }

  /**************************************************************************/

  // Allocate memory for the Coulomb matrix elements
  if ((vijck = (double *) calloc(ist.nHoles*ist.nHoles*ist.nElecs*ist.nHotHoles, sizeof(double))) == NULL) nerror("vijck");
  if ((vabck = (double *) calloc(ist.nElecs*ist.nElecs*ist.nHoles*ist.nHotElecs, sizeof(double))) == NULL) nerror("vabck");

  /**************************************************************************/

  // Calculate the AR lifetimes 
  if (ist.deterministic) {
    calcAllCoulombMatrixElements(vijck, vabck, psiai, evalai, psibe, potq, poth, ist, par, planfw, planbw, fftwpsi);
    if (ist.intAR) {
      calcDeterministicIntAR(Cbs, Ebs, vijck, vabck, evalbe, evalai, ist, par);
    }
    else if (ist.nonIntAR) {
      calcDeterministicNonIntAR(vijck, vabck, evalbe, evalai, ist, par);
    }
  }
  if (ist.stochastic) {
    if (ist.intAR) {
      if (ist.stochastic == 1) {
        calcStochasticFinalStatesIntAR(Cbs, Ebs, vijck, vabck, psibe, evalbe, psiai, evalai, potq, poth, ist, par, planfw, planbw, fftwpsi);
      }
      else if (ist.stochastic == 2) {
        calcStochasticCoulombIntAR(Cbs, Ebs, psibe, evalbe, psiai, evalai, potq, vx, vy, vz, ist, par, planfw, planbw, fftwpsi);
      }
      else if (ist.stochastic == 3) {
        calcDoublyStochasticIntAR(Cbs, Ebs, psibe, evalbe, psiai, evalai, potq, ist, par, planfw, planbw, fftwpsi);
      }    
    }
    else if (ist.nonIntAR) {
      calcStochasticNonIntAR(vijck, vabck, psibe, evalbe, psiai, evalai, potq, ist, par, planfw, planbw, fftwpsi);
    } 
  }
      
  /**************************************************************************/

  // Free dynamically allocated memory
  if (ist.intAR) {
    free(Ebs); free(Cbs);
  }
  free(potq); free(potl);  free(psibe);  free(psiai); free(evalbe);
  free(ksqr); free(vx); free(vy);  free(vz);  free(vabck); free(vijck);
  free(rx); free(ry); free(rz); free(poth); free(evalai); free(sigeai);
  for (i = 0; i < ist.nThreads; i++){ 
    fftw_destroy_plan(planfw[i]);
    fftw_destroy_plan(planbw[i]);
  }
  fftw_free(fftwpsi);

  /**************************************************************************/
  writeSeparation(stdout);
  writeCurrentTime(stdout);

  exit(EXIT_SUCCESS);
}

/**************************************************************************/
