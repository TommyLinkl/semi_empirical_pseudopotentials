/*****************************************************************************/
//
// This program calculates correlated two quasiparticle states. The options are:
//   1. elec-hole (i.e. solves the bethe-salpeter equation for excitons)
//   2. elec-elec (correlated two quasielectron states)
//   3. hole-hole (correlated two hole states)
//
// Input files:
//   input.par    -> required
//   eval.par     -> required
//   psi.par      -> required
//   conf.par     -> optional
//   grid.par     -> optional
//   spinEval.par -> optional
//
// Output files: 
//   stdout ->
//    ->
//    ->    
//
// Author: John P. Philbin
// Began: January 10th, 2019
// Last modified: August 28th, 2019
//
/*****************************************************************************/

#include "qp.h"

/*****************************************************************************/

int main(int argc, char *argv[]) {
  long i, flags = 0; 
  double *psiHoles, *psiElecs;
  double *Wrsut, *Vrsut, *h0Matrix, *hMatrix;
  zomplex *rSpaceCoulombPot, *rSpaceScreenedCoulombPot;
  zomplex *qSpaceCoulombPot, *qSpaceScreenedCoulombPot;
  nonintQP *holeQP, *elecQP;
  nonintTwoQPState *nonintTwoQP;
  intTwoQPState *intTwoQP;
  nonintExc *nonintExciton;
  intExc *intExciton;
  gridPoint *rSpaceGP, *qSpaceGP;
  grid rSpaceGrid, qSpaceGrid;
  dParams dPar;  
  lParams lPar;
  fftw_plan_loc *planfw, *planbw; 
  fftw_complex *fftwPsi;

  /*************************************************************************/  
  // Write time the program began and read input.par 
  writeCurrentTime(stdout); writeSeparation(stdout); fflush(stdout);
  readInputParFile(&rSpaceGrid, &dPar, &lPar);
    
  /*************************************************************************/
  // Dynamically allocate memory

  // Allocate memory for the noninteracting quasiparticle states
  if ((holeQP = (nonintQP *) calloc(lPar.nHoles, sizeof(nonintQP))) == NULL) memoryError("holeQP");
  if ((elecQP = (nonintQP *) calloc(lPar.nElecs, sizeof(nonintQP))) == NULL) memoryError("elecQP");
  if ((nonintTwoQP = (nonintTwoQPState *) calloc(lPar.nNonintTwoQPStates, sizeof(nonintTwoQPState))) == NULL) memoryError("nonintTwoQP");

  // Allocate memory for the wavefunctions
  if ((psiHoles = (double *) calloc(lPar.nHoles*rSpaceGrid.nGridPoints, sizeof(double))) == NULL) memoryError("psiHoles");
  if ((psiElecs = (double *) calloc(lPar.nElecs*rSpaceGrid.nGridPoints, sizeof(double))) == NULL) memoryError("psiElecs");

  // Allocate memory for the 3D Fast Fourier transforms using openMP
  if ((fftwPsi = fftw_malloc(sizeof (fftw_complex )*rSpaceGrid.nGridPoints*lPar.nThreads)) == NULL) memoryError("fftwPsi");
  fftw_plan_with_nthreads(lPar.nThreads);
  planfw = (fftw_plan_loc *) calloc(lPar.nThreads, sizeof(fftw_plan_loc));
  planbw = (fftw_plan_loc *) calloc(lPar.nThreads, sizeof(fftw_plan_loc));
  for (i = 0; i < lPar.nThreads; i++) { 
    planfw[i] = fftw_plan_dft_3d(rSpaceGrid.nGridPointsZ, rSpaceGrid.nGridPointsY, rSpaceGrid.nGridPointsX, 
                                  &fftwPsi[i*rSpaceGrid.nGridPoints], &fftwPsi[i*rSpaceGrid.nGridPoints], FFTW_FORWARD,  flags);
    planbw[i] = fftw_plan_dft_3d(rSpaceGrid.nGridPointsZ, rSpaceGrid.nGridPointsY, rSpaceGrid.nGridPointsX, 
                                  &fftwPsi[i*rSpaceGrid.nGridPoints], &fftwPsi[i*rSpaceGrid.nGridPoints], FFTW_BACKWARD, flags);
  }

  // Allocate memory for the grid points
  if ((rSpaceGP = (gridPoint *) calloc(rSpaceGrid.nGridPoints, sizeof(gridPoint))) == NULL) memoryError("rSpaceGP");
  if ((qSpaceGP = (gridPoint *) calloc(rSpaceGrid.nGridPoints, sizeof(gridPoint))) == NULL) memoryError("qSpaceGP");

  // Allocate memory for the Coulomb potentials
  if ((rSpaceCoulombPot = (zomplex *) calloc(rSpaceGrid.nGridPoints, sizeof(zomplex))) == NULL) memoryError("rSpaceCoulombPot");
  if ((qSpaceCoulombPot = (zomplex *) calloc(rSpaceGrid.nGridPoints, sizeof(zomplex))) == NULL) memoryError("qSpaceCoulombPot");
  if ((rSpaceScreenedCoulombPot = (zomplex *) calloc(rSpaceGrid.nGridPoints, sizeof(zomplex))) == NULL) memoryError("rSpaceScreenedCoulombPot");
  if ((qSpaceScreenedCoulombPot = (zomplex *) calloc(rSpaceGrid.nGridPoints, sizeof(zomplex))) == NULL) memoryError("qSpaceScreenedCoulombPot");

  // Allocate memory for the Coulomb matrix elements and Hamiltonians
  if ((Wrsut = (double *) calloc(lPar.nNonintTwoQPStates*lPar.nNonintTwoQPStates, sizeof(double))) == NULL) memoryError("Wrsut");
  if ((Vrsut = (double *) calloc(lPar.nNonintTwoQPStates*lPar.nNonintTwoQPStates, sizeof(double))) == NULL) memoryError("Vrsut");
  if ((h0Matrix = (double *) calloc(lPar.nNonintTwoQPStates*lPar.nNonintTwoQPStates, sizeof(double))) == NULL) memoryError("h0Matrix");
  if ((hMatrix = (double *) calloc(lPar.nNonintTwoQPStates*lPar.nNonintTwoQPStates, sizeof(double))) == NULL) memoryError("hMatrix");

  // Allocate memory for the interacting (i.e. correlated) two quasiparticle states
  if ((intTwoQP = (intTwoQPState *) calloc(lPar.nIntTwoQPStates, sizeof(intTwoQPState))) == NULL) memoryError("intTwoQP");

  /**************************************************************************/
  // Read in the noninteracting quasiparticle states  
  readEvalPsiParFiles(psiHoles, psiElecs, holeQP, elecQP, dPar.sigmaCutoff, lPar);
  
  /**************************************************************************/
  // Initialize the noninteracting two quasiparticle states
  fillNonintTwoQPStructure(nonintTwoQP, holeQP, elecQP, lPar);

  /**************************************************************************/
  // Initialize the r-space grid 
  initRSpaceGrid(&rSpaceGrid, rSpaceGP, lPar);
  initQSpaceGrid(&qSpaceGrid, qSpaceGP, rSpaceGrid);

  /**************************************************************************/
  // Normalize the wavefunctions
  normalizeAllDoubleWavefunctions(psiHoles, rSpaceGrid.nGridPoints, rSpaceGrid.dV, lPar.nHoles);
  normalizeAllDoubleWavefunctions(psiElecs, rSpaceGrid.nGridPoints, rSpaceGrid.dV, lPar.nElecs);

  /**************************************************************************/
  // Calculate the noninteracting (single-particle) quasiparticle probability densities
  // calcNonintAbsorptionProperties(nonintTwoQP, rSpaceGrid, qSpaceGrid, dPar, lPar, planfw, planbw, fftwPsi);
  calcNonintQPDensities(holeQP, elecQP, rSpaceGrid, dPar, lPar);
  if (lPar.calcIntDensitiesOnly && (! strcmp(lPar.calcType, "elec-hole"))) {
    readIntTwoQuasiparticleParFiles(intTwoQP, hMatrix, nonintTwoQP, dPar, lPar);
    calcAbsorptionProperties(intTwoQP, nonintTwoQP, rSpaceGrid, qSpaceGrid, dPar, lPar, planfw, planbw, fftwPsi);
    calcExcitonSizeStats(intTwoQP, nonintTwoQP, psiHoles, psiElecs, rSpaceGrid, dPar, lPar);
    calcIntTwoQPProjDensities(intTwoQP, rSpaceGrid, dPar, lPar);
    exit(EXIT_SUCCESS);
  }

  /**************************************************************************/
  // Calculate the r- and q-space Coulomb potentials on their respective grids
  calcRSpaceCoulombPotential(rSpaceCoulombPot, rSpaceGrid, dPar.epsilon, 0);        
  calcRSpaceCoulombPotential(rSpaceScreenedCoulombPot, rSpaceGrid, dPar.epsilon, 1); 
  calcQSpaceCoulombPotential(qSpaceCoulombPot, qSpaceGrid, rSpaceCoulombPot, rSpaceGrid, 
                              dPar.epsilon, 0, planfw, planbw, fftwPsi);             
  calcQSpaceCoulombPotential(qSpaceScreenedCoulombPot, qSpaceGrid, rSpaceScreenedCoulombPot, rSpaceGrid, 
                              dPar.epsilon, 1, planfw, planbw, fftwPsi);  
  
  /**************************************************************************/
  // Calculate the Coulomb matrix elements between the two quasiparticles
  calcAllCoulombMatrixElements(Wrsut, Vrsut, nonintTwoQP, psiHoles, psiElecs, qSpaceCoulombPot, qSpaceScreenedCoulombPot,
                                qSpaceGrid, rSpaceGrid, lPar, planfw, planbw, fftwPsi);

  /**************************************************************************/
  if (! strcmp(lPar.calcType, "auger-decay-ni") || ! strcmp(lPar.calcType, "auger-decay-i")) {
    calcAugerDecayNoninteracting(nonintTwoQP, Wrsut, Vrsut, dPar, lPar);
  }
  else {
    // Fill the Hamiltonian matrix with h0 and the direct and indirect Coulomb matrix elements 
    calcTwoQPH0Matrix(h0Matrix, nonintTwoQP, lPar.nNonintTwoQPStates);
    calcTwoQPHamiltonianWrapper(hMatrix, h0Matrix, Wrsut, Vrsut, nonintTwoQP, lPar);
    
    /**************************************************************************/
    // Diagonalize the Hamiltonian matrix to obtain the correlated two quasiparticle states
    diagonalizeTwoQPHamiltonain(intTwoQP, hMatrix, nonintTwoQP, lPar);

    /**************************************************************************/
    // Calculate exciton related properties results
    if (! strcmp(lPar.calcType, "elec-hole") || ! strcmp(lPar.calcType, "sp-elec-hole")) {
      calcAbsorptionProperties(intTwoQP, nonintTwoQP, rSpaceGrid, qSpaceGrid, dPar, lPar, planfw, planbw, fftwPsi);
      calcExcitonSizeStats(intTwoQP, nonintTwoQP, psiHoles, psiElecs, rSpaceGrid, dPar, lPar);
      calcIntTwoQPProjDensities(intTwoQP, rSpaceGrid, dPar, lPar);
    }
    if (lPar.nFixedQPPoints) {
      calcQPDensitiesFixedOtherQP(intTwoQP, psiHoles, psiElecs, rSpaceGrid, dPar, lPar);        
    }
  }

  /**************************************************************************/
  // Free the dynamically allocated memory

  // Free memory related to the interacting two quasiparticle states
  free(intTwoQP); 

  // Free memory related to the Coulomb matrix elements and Hamiltonians
  free(hMatrix); free(h0Matrix); 
  free(Vrsut); free(Wrsut);

  // Free memory related to the Coulomb potentials
  free(qSpaceScreenedCoulombPot); free(qSpaceCoulombPot);
  free(rSpaceScreenedCoulombPot); free(rSpaceCoulombPot);

  // Free memory related to the grid points
  free(rSpaceGP); free(qSpaceGP);

  // Free memory related to the Fast Fourier Transforms
  for (i = 0; i < lPar.nThreads; i++) { 
    fftw_destroy_plan(planfw[i]);
    fftw_destroy_plan(planbw[i]);
  }
  fftw_free(fftwPsi);

  // Free memory related to the quasiparticle states  
  free(psiElecs); free(psiHoles);
  free(elecQP); free(holeQP);
  free(nonintTwoQP);

  /**************************************************************************/
  // Write program timings 
  writeSeparation(stdout);
  writeCurrentTime(stdout);

  /**************************************************************************/
  // Exit the program
  exit(EXIT_SUCCESS);

}

/**************************************************************************/
