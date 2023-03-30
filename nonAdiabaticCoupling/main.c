/*****************************************************************************/
// TODO: update this description
//
// This program calculates correlated two quasiparticle states. The options are:
//   1. elec-hole (i.e. solves the bethe-salpeter equation for excitons)
//   2. elec-elec (correlated two quasielectron states)
//   3. hole-hole (correlated two hole states)
//
// Input files:
//   input.par   ->  
//   eval.par    -> a list 
//   psi.par     -> 
//
// Output files: 
//   stdout ->
//    ->
//    ->    
//
// Began: May 13th, 2019
// Last modified: May 13th, 2019
//
/*****************************************************************************/

#include "qp.h"

/*****************************************************************************/
int main(long argc, char *argv[]) {
  double *psiHoles, *psiElecs;
  nonintQP *holeQP, *elecQP;
  nonintTwoQPState *nonintTwoQP;
  intTwoQPState *intTwoQP;
  vector *dPotdAtomPos;
  vector *Vab, *Vij;
  gridPoint3d *rSpaceGP; 
  grid3d rSpaceGrid;
  grid1d atomicPPGrid, dAtomicPPGrid;
  gridPoint1d *atomicPPGP, *dAtomicPPGP; 
  atom *atoms;
  dParams dPar;  
  lParams lPar;

  /*************************************************************************/  
  // Write time the program began and read input.par 

  writeCurrentTime(stdout);
  writeSeparation(stdout);
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

  // Allocate memory for the grid points
  if ((rSpaceGP = (gridPoint3d *) calloc(rSpaceGrid.nGridPoints, sizeof(gridPoint3d))) == NULL) memoryError("rSpaceGP");

  // Allocate memory for the interacting (i.e. correlated) two quasiparticle states
  if ((intTwoQP = (intTwoQPState *) calloc(lPar.nIntTwoQPStates, sizeof(intTwoQPState))) == NULL) memoryError("intTwoQP");
   
  // Allocate memory for the derivative of the potential energy 
  if ((dPotdAtomPos = (vector *) calloc(rSpaceGrid.nGridPoints, sizeof(vector))) == NULL) memoryError("dPotdAtomPos"); 

  // Allocate memory for the derivative coupling matrix elements
  if ((Vab = (vector *) calloc(lPar.nElecs*lPar.nElecs*lPar.nAtoms, sizeof(vector))) == NULL) memoryError("Vab"); 
  if ((Vij = (vector *) calloc(lPar.nHoles*lPar.nHoles*lPar.nAtoms, sizeof(vector))) == NULL) memoryError("Vij");   

  // Allocate memory for the atomic configuration
  if ((atoms = (atom *) calloc(lPar.nAtoms, sizeof(atom))) == NULL) memoryError("atoms");

  /**************************************************************************/
  // Read in the noninteracting quasiparticle states  
  readEvalPsiParFiles(psiHoles, psiElecs, holeQP, elecQP, dPar.sigmaCutoff, lPar);

  /**************************************************************************/
  // Read in the atom positions
  fillAtomStructureFromConfFile(atoms, &lPar);

  // Allocate memory for the atomic pseudopotentials and their derivatives
  atomicPPGrid.nGridPoints = dAtomicPPGrid.nGridPoints = 1024; // TODO: have more flexibility here
  if ((atomicPPGP = (gridPoint1d *) calloc(atomicPPGrid.nGridPoints*lPar.nSCAtomTypes, sizeof(gridPoint1d))) == NULL) memoryError("atomicPPGP");
  if ((dAtomicPPGP = (gridPoint1d *) calloc(dAtomicPPGrid.nGridPoints*lPar.nSCAtomTypes, sizeof(gridPoint1d))) == NULL) memoryError("dAtomicPPGP");

  // Read in the pseudopotentials for all atom types
  readAtomicPseudopotentials(&atomicPPGrid, atomicPPGP, atoms, lPar);

  // Calculate the derivative of the pseudopotentials for all atom types
  calcPseudopotentialDerivative(&dAtomicPPGrid, dAtomicPPGP, atomicPPGrid, lPar.nSCAtomTypes);

  /**************************************************************************/
  // Initialize the noninteracting two quasiparticle states
  //fillNonintTwoQPStructure(nonintTwoQP, holeQP, elecQP, lPar);

  /**************************************************************************/
  // Initialize the r-space grid 
  initRSpaceGrid(&rSpaceGrid, rSpaceGP, lPar);

  /**************************************************************************/ 
  if (! strcmp(lPar.stateType, "adiabatic")) {
    calcDerivativeCouplings(Vab, Vij, dPotdAtomPos, atoms, dAtomicPPGrid, holeQP, elecQP, rSpaceGrid, lPar);
  }
  else if (! strcmp(lPar.stateType, "diabatic")) {
    calcElPh(Vab, Vij, dPotdAtomPos, atoms, dAtomicPPGrid, holeQP, elecQP, rSpaceGrid, lPar);
  }

  /**************************************************************************/
  // Free the dynamically allocated memory

  // Free memory related to the interacting two quasiparticle states
  free(intTwoQP); 

  // Free memory related to the grid points
  free(rSpaceGP); 

  // Free memory related to the nonadiabatic coupling matrix elements
  free(Vij); free(Vab);
  free(dPotdAtomPos);

  // Free memory related to the atomic configuration
  free(atoms);

  // Free memory related to the atomic pseudopotentials
  free(atomicPPGP); free(dAtomicPPGP);

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
