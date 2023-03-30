/*****************************************************************************/
//
// This program calculates properties related to the size of excitons 
//
// Input files:
//   input.par   ->  
//   eval.par    -> a list 
//   psi.par     -> 
//   exciton.par -> 
//
// Output files: 
//   stdout ->
//    ->
//    ->    
//
// Author: John P. Philbin
// Last modified: Jan. 2nd, 2019
//
/*****************************************************************************/

#include "es.h"

/*****************************************************************************/

int main(long argc, char *argv[]) {
  long i, flags = 0; 
  double *psiHoles, *psiElecs, *Cai, *Ebs;
  nonintQP *holeQP, *elecQP;
  nonintExc *nonintExciton;
  intExc *intExciton;
  gridPoint *rSpaceGP, *qSpaceGP;
  grid rSpaceGrid, qSpaceGrid;
  dParams dPar;  
  lParams lPar;

  /*************************************************************************/  
  // Write time the program began and read input.par 

  writeCurrentTime(stdout);
  writeSeparation(stdout);
  readInputParFile(&rSpaceGrid, &dPar, &lPar);
    
  /*************************************************************************/
  // Dynamically allocate memory

  // Allocate memory for the noninteracting quasiparticles and excitons
  if ((holeQP = (nonintQP *) calloc(lPar.nHoles, sizeof(nonintQP))) == NULL) memoryError("holeQP");
  if ((elecQP = (nonintQP *) calloc(lPar.nElecs, sizeof(nonintQP))) == NULL) memoryError("elecQP");
  if ((nonintExciton = (nonintExc *) calloc(lPar.nNonintExcitons, sizeof(nonintExc))) == NULL) memoryError("nonintExciton");

  // Allocate memory for the wavefunctions
  if ((psiHoles = (double *) calloc(lPar.nHoles*rSpaceGrid.nGridPoints, sizeof(double))) == NULL) memoryError("psiHoles");
  if ((psiElecs = (double *) calloc(lPar.nElecs*rSpaceGrid.nGridPoints, sizeof(double))) == NULL) memoryError("psiElecs");

  // Allocate memory for the grid points
  if ((rSpaceGP = (gridPoint *) calloc(rSpaceGrid.nGridPoints, sizeof(gridPoint))) == NULL) memoryError("rSpaceGP");
  if ((qSpaceGP = (gridPoint *) calloc(rSpaceGrid.nGridPoints, sizeof(gridPoint))) == NULL) memoryError("qSpaceGP");

  // Allocate memory for the coefficients and interacting excitons from the BSE results
  if ((Ebs = (double *) calloc(lPar.nIntExcitons, sizeof(double))) == NULL) memoryError("Ebs");
  if ((Cai = (double *) calloc(lPar.nNonintExcitons*lPar.nIntExcitons, sizeof(double))) == NULL) memoryError("Cai");
  if ((intExciton = (intExc *) calloc(lPar.nIntExcitons, sizeof(intExc))) == NULL) memoryError("intExciton");
    
  /**************************************************************************/
  // Read in and initialize the noninteracting QP states and excitonic states 

  // Read in the noninteracting QP states
  readEvalPsiParFiles(psiHoles, psiElecs, holeQP, elecQP, dPar.sigmaCutoff, lPar);
  lPar.iHomo = lPar.nHoles-1; lPar.iLumo = lPar.nHoles;
  fillNonintExcitonStructure(nonintExciton, holeQP, elecQP, lPar);

  // Read in interacting exciton information 
  if (lPar.intExcitons) {
    readBetheSalpeterResults(Cai, Ebs, lPar);
    fillIntExcitonStructure(intExciton, Cai, Ebs, lPar);
  }

  /**************************************************************************/
  // Initialize the r-space grid 
  initRSpaceGrid(&rSpaceGrid, rSpaceGP, lPar);
 
  /**************************************************************************/
  // Calculate properties related to the size of the excitonic states 
  calcExcitonSizeStats(intExciton, nonintExciton, psiHoles, psiElecs, rSpaceGrid, dPar, lPar);

  /**************************************************************************/
  // Free the dynamically allocated memory

  // Interacting exciton related memory
  free(intExciton); free(Cai); free(Ebs);

  // Grid point related memory
  free(rSpaceGP); free(qSpaceGP);

  // Quasiparticle related memory  
  free(psiElecs); free(psiHoles);
  free(nonintExciton); 
  free(elecQP); free(holeQP); 

  /**************************************************************************/
  // Write program timings 

  writeSeparation(stdout);
  writeCurrentTime(stdout);

  /**************************************************************************/
  // Exit the program

  exit(EXIT_SUCCESS);
}

/**************************************************************************/