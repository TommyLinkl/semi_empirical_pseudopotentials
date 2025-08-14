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
  double *psiHoles, *psiElecs, *a4Params, *a5Params;
  nonintQP *holeQP, *elecQP;
  nonintTwoQPState *nonintTwoQP;
  intTwoQPState *intTwoQP;
  vector *dPotdAtomPos;
  vector *Vab, *Vij;
  gridPoint3d *rSpaceGP; 
  grid3d rSpaceGrid;
  grid1d atomicPPGrid, dAtomicPPGrid;
  gridPoint1d *atomicPPGP, *dAtomicPPGP; 
  atom *atoms, *atomNearestNeighbors;
  dParams dPar;  
  lParams lPar;
  double *refTetrahedronVol, *tetrahedronVol, *strainScale;
  vector *tetrahedronVolDerivatives, *strainScaleDerivatives;
  int crystalStructureInt, outmostMaterialInt;

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
  if ((a4Params = (double *) calloc(lPar.nSCAtomTypes, sizeof(double))) == NULL) memoryError("a4Params");
  if ((a5Params = (double *) calloc(lPar.nSCAtomTypes, sizeof(double))) == NULL) memoryError("a5Params");

  // Read in the pseudopotentials for all atom types
  readAtomicPseudopotentials(&atomicPPGrid, atomicPPGP, atoms, lPar, a4Params, a5Params);

  // Calculate the derivative of the pseudopotentials for all atom types
  calcPseudopotentialDerivative(&dAtomicPPGrid, dAtomicPPGP, atomicPPGrid, lPar.nSCAtomTypes);

  /**************************************************************************/
  /* Strain-dependent related computation */
  if ((atomNearestNeighbors = (atom *) calloc(4*lPar.nAtoms, sizeof(atom))) == NULL) memoryError("atomNearestNeighbors");
  if ((refTetrahedronVol = (double *) calloc(lPar.nAtoms, sizeof(double))) == NULL) memoryError("refTetrahedronVol");
  if ((tetrahedronVol = (double *) calloc(lPar.nAtoms, sizeof(double))) == NULL) memoryError("tetrahedronVol");
  if ((tetrahedronVolDerivatives = (vector *) calloc(4*lPar.nAtoms, sizeof(vector))) == NULL) memoryError("tetrahedronVolDerivatives");
  if ((strainScale = (double *) calloc(lPar.nAtoms, sizeof(double))) == NULL) memoryError("strainScale");
  if ((strainScaleDerivatives = (vector *) calloc(4*lPar.nAtoms, sizeof(vector))) == NULL) memoryError("strainScaleDerivatives");

  if (! strcmp(lPar.crystalStructure, "wurtzite")) {
      crystalStructureInt = 0;
  }
  else if (! strcmp(lPar.crystalStructure, "zincblende")) {
      crystalStructureInt = 1;
  }
  else {
      printf("\n\nCrystal structure type %s not recognized -- the program is exiting!!!\n\n", lPar.crystalStructure);
      fflush(stdout);
      exit(EXIT_FAILURE);
  }

  /*** Assign outmostMaterialInt ***/
  if (! strcmp(lPar.outmostMaterial, "CdS")) {
    outmostMaterialInt = 0;
  }
  else if (! strcmp(lPar.outmostMaterial, "CdSe")) {
    outmostMaterialInt = 1;
  }
  else if (! strcmp(lPar.outmostMaterial, "InP")) {
    outmostMaterialInt = 2;
  }
  else if (! strcmp(lPar.outmostMaterial, "InAs")) {
    outmostMaterialInt = 3;
  }
  else if (! strcmp(lPar.outmostMaterial, "alloyInGaP")) {  //cation terminated surface only
    outmostMaterialInt = 4;
  }
  else if (! strcmp(lPar.outmostMaterial, "alloyInGaAs")) {  // cation terminated surface only
    outmostMaterialInt = 5;
  }
  else if (! strcmp(lPar.outmostMaterial, "GaAs")) {  // cation terminated surface only
    outmostMaterialInt = 6;
  }
  else {
    printf("\n\nOutmostMaterial type %s not recognized -- the program is exiting!!!\n\n", lPar.outmostMaterial);
    fflush(stdout);
    exit(EXIT_FAILURE);
  }

  readNearestNeighbors(lPar.nAtoms, atomNearestNeighbors);
  calculateRefTetrahedronVol(lPar.nAtoms, crystalStructureInt, outmostMaterialInt, atoms, atomNearestNeighbors, refTetrahedronVol);
  calculateTetrahedronVol(lPar.nAtoms, atoms, atomNearestNeighbors, tetrahedronVol);
  calculateTetrahedronVolDeriv(lPar.nAtoms, atoms, atomNearestNeighbors, tetrahedronVol, tetrahedronVolDerivatives);
  calculateStrainScale(lPar.nAtoms, atoms, refTetrahedronVol, tetrahedronVol, a4Params, a5Params, strainScale);
  calculateStrainScaleDeriv(lPar.nAtoms, atoms, atomNearestNeighbors, refTetrahedronVol, tetrahedronVol, tetrahedronVolDerivatives,
          a4Params, a5Params, strainScaleDerivatives);

  printf("Finished with strain-related calculations\n");
  fflush(stdout);
  
  free(refTetrahedronVol); free(tetrahedronVol); free(tetrahedronVolDerivatives); free(a4Params); free(a5Params);

  /**************************************************************************/
  // Initialize the noninteracting two quasiparticle states
  //fillNonintTwoQPStructure(nonintTwoQP, holeQP, elecQP, lPar);

  /**************************************************************************/
  // Initialize the r-space grid 
  initRSpaceGrid(&rSpaceGrid, rSpaceGP, lPar);
  printf("Finished with initializing r-space grid\n"); 

  /**************************************************************************/ 
  calcElPh(Vab, Vij, dPotdAtomPos, atoms, atomNearestNeighbors, atomicPPGrid, dAtomicPPGrid, 
          holeQP, elecQP, rSpaceGrid, strainScale, strainScaleDerivatives, lPar);

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
 
  // Free memory related to strain 
  free(atomNearestNeighbors); 
  free(strainScale); free(strainScaleDerivatives);

  /**************************************************************************/
  // Write program timings 
  writeSeparation(stdout);
  writeCurrentTime(stdout);

  /**************************************************************************/
  // Exit the program
  exit(EXIT_SUCCESS);

}

/**************************************************************************/
