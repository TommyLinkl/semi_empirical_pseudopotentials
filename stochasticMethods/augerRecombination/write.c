/****************************************************************************/

// This file does the printing for the program 

/*****************************************************************************/

#include "ar.h"

/****************************************************************************/
// Writes the input parameters that were used in AR lifetime calculation

void writeInputParameters(lng_st ist, par_st par, FILE *pf) {

  fprintf(pf, "Calculation type = %s\n", ist.calcType);
  if (ist.deterministic) {
    fprintf(pf, "Calculate AR lifetimes deterministically\n");
  }
  if (ist.stochastic) {
    fprintf(pf, "Calculate AR lifetimes stochastically\n");
    if (ist.stochastic == 2 || ist.stochastic == 3) {
      fprintf(pf, "The number of stochastic orbitals used to appox. the Coulomb operator = %ld\n", ist.nStochOrbitals);
    }
    if (ist.stochastic == 1 || ist.stochastic == 3) {
      fprintf(pf, "The number of stochastic hot hole orbitals = %ld\n", ist.nStochHotHoles);
      fprintf(pf, "The number of stochastic hot elec orbitals = %ld\n", ist.nStochHotElecs);
    }
  }
  fprintf(pf, "\nThe number of openMP threads used = %ld\n", ist.nThreads);
  if (ist.readInHotQPs) {
    fprintf(pf, "Using previously calculated hot hole and electron states, readInHotQPs = %ld\n", ist.readInHotQPs);
  }
  if (ist.readAaiMatrices) {
    fprintf(pf, "Using previously calculated AaiHole and AaiElec matrices, readAaiMatrices = %ld\n", ist.readAaiMatrices);
  }
  fprintf(pf, "The seed used = %ld\n", ist.seed);
  fprintf(pf, "Maximum length of pseudopotential files, nPseudoPot = %ld\n", ist.nPseudoPot);
  fprintf(pf, "The number of atoms in the system, nAtoms = %ld\n", ist.nAtoms);
  fprintf(pf, "The number of grid points used: nx = %ld  ny = %ld  nz = %ld\n", ist.nx, ist.ny, ist.nz);
  fprintf(pf, "The number of filter cycles, nFilterCycles = %ld\n", ist.nFilterCycles);
  fprintf(pf, "The number of states per filter, nStatesPerFilter = %ld\n", ist.nStatesPerFilter);
  fprintf(pf, "Total number of filtered states, nFilteredStates = %ld\n", ist.nFilteredStates);
  fprintf(pf, "Length of Newton interpolation used, nNewtonIntSteps = %ld\n", ist.nNewtonIntSteps);
  fprintf(pf, "Kinetic energy maximum = %.2f\n", par.Ekinmax);
  fprintf(pf, "Electronic temperature = %.1f\n", par.temp);
  fprintf(pf, "kbT = %.8f gives an energy range of, 3kbT = %.8f\n", par.kbT, par.boltzEnergyRange);
  fprintf(pf, "Sigma cutoff for eigenstate determination, sigmaCutoff = %.3f\n\n", par.sigmaCutoff);

  fprintf(pf, "Total number of hole eigenstates = %ld\n", ist.totalHomo);
  fprintf(pf, "Total number of electron eigenstates = %ld\n", ist.totalLumo);
  fprintf(pf, "The index of the HOMO state in eval.par = %ld\n", ist.homoIndex);
  fprintf(pf, "The index of the LUMO state in eval.par = %ld\n", ist.lumoIndex);
  fprintf(pf, "HOMO energy = %.10f\n", par.homoEnergy);
  fprintf(pf, "LUMO energy = %.10f\n", par.lumoEnergy);
  fprintf(pf, "Fundamental gap = %.10f %.8f\n", par.fundamentalGap, par.fundamentalGap*AUTOEV);
  if (ist.nonIntAR) {
    fprintf(pf, "The number of holes within 3kbT of the HOMO = %ld\n", ist.nHoles);
    fprintf(pf, "The number of elecs within 3kbT of the LUMO = %ld\n", ist.nElecs);
    fprintf(pf, "Total number of noninteracting states within boltzEnergyRange = %ld\n", ist.nHolesPlusElecs);
  }
  else if (ist.intAR) {
    fprintf(pf, "Optical gap = %.10f %.8f\n", par.opticalGap, par.opticalGap*AUTOEV);
    fprintf(pf, "Exciton binding energy = %.4f\n", (par.fundamentalGap-par.opticalGap)*1000.0*AUTOEV);    
    fprintf(pf, "The number of holes used in the BSE = %ld\n", ist.nHoles);
	  fprintf(pf, "The number of elecs used in the BSE = %ld\n", ist.nElecs);
    fprintf(pf, "The number of noninteracting excitons = %ld\n", ist.nNonIntExcitons);
    fprintf(pf, "The number of interacting excitons within 3kbT of the optical gap = %ld\n", ist.nIntExcitons);
  }
  fprintf(pf, "Min energy of the initial biexciton = %.10f %.8f\n", par.minInitE, par.minInitE*AUTOEV);
  fprintf(pf, "Max energy of the initial biexciton = %.10f %.8f\n\n", par.maxInitE, par.maxInitE*AUTOEV);

  fprintf(pf, "The maximum energy conservation window used = %.6f\n", par.maxDeltaE);
  fprintf(pf, "The total number of energy conservation windows = %ld\n\n", ist.nConsWindows);

  return;
}

/****************************************************************************/
// prints out information regarding a generate-filter-states function call

void writeFilterInfo(double *targetEnergies, lng_st ist, par_st par, FILE *pf) {
  long jms;
  double dt = sqr((double)(ist.nNewtonIntSteps) / (2.5*par.dE));

  fprintf(pf, "Newton interpolation length = %ld\n", ist.nNewtonIntSteps);
  fprintf(pf, "Energy range of the Hamiltonian, dE = % .4f\n", par.dE);
  fprintf(pf, "dt = % .4f\n", dt);
  fprintf(pf, "The minimum hot elec energy = % .6f and maximum = % .6f\n", par.Eamin, par.Eamax);
  fprintf(pf, "The %ld elec filter target energies are = ", ist.nStatesPerFilter/2);
  for (jms = 0; jms < ist.nStatesPerFilter/2; jms++) {
    fprintf(pf, "% .6f ", targetEnergies[jms]); 
  }
  fprintf(pf, "\nThe minimum hot hole energy = % .6f and maximum = % .6f\n", par.Eimin, par.Eimax);
  fprintf(pf, "The %ld hole filter target energies are = ", ist.nStatesPerFilter/2);
  for (jms = ist.nStatesPerFilter/2; jms < ist.nStatesPerFilter; jms++) {
    fprintf(pf, "% .6f ", targetEnergies[jms]);
  }
  fprintf(pf, "\n");

  return;
}

/****************************************************************************/
// prints out the noninteracting biexciton to pf
 
void writeNonIntBiexciton(nonIntBiexc biExciton, FILE *pf) {
  
  fprintf(pf, "\nBiexciton index = %d\n", biExciton.index);
  fprintf(pf, "Biexciton energy = % .6f % .4f\n", biExciton.energy, biExciton.energy*AUTOEV);
  fprintf(pf, "Elec 1 index, b = %d and energy = % .6f % .4f\n", biExciton.b, biExciton.bEnergy, biExciton.bEnergy*AUTOEV);
  fprintf(pf, "Elec 2 index, c = %d and energy = % .6f % .4f\n", biExciton.c, biExciton.cEnergy, biExciton.cEnergy*AUTOEV);
  fprintf(pf, "Hole 1 index, j = %d and energy = % .6f % .4f\n", biExciton.j, biExciton.jEnergy, biExciton.jEnergy*AUTOEV);
  fprintf(pf, "Hole 2 index, k = %d and energy = % .6f % .4f\n", biExciton.k, biExciton.kEnergy, biExciton.kEnergy*AUTOEV);

  return;
}

/****************************************************************************/
// writes the array to a file corresponding to the string in fileName
 
void writeLongArray(long *longArray, long arrayLength, char *fileNameBase) {
  FILE *pf;
  long i;
  char fileName[100];

  //if ()
  sprintf(fileName, "%s.dat", fileNameBase);

  pf = fopen(fileName, "w");
  for (i = 0; i < arrayLength; i++) {
    fprintf(pf, "%ld %ld\n", i, longArray[i]);
  }
  fclose(pf);

  return;
}

/****************************************************************************/
// writes the array to a file corresponding to the string in fileName
 
void writeDoubleArray(double *doubleArray, long arrayLength, char *fileNameBase) {
  FILE *pf;
  long i;
  char fileName[100];

  //if ()
  sprintf(fileName, "%s.dat", fileNameBase);

  pf = fopen(fileName, "w");
  for (i = 0; i < arrayLength; i++) {
    fprintf(pf, "%ld %g\n", i, doubleArray[i]);
  }
  fclose(pf);

  return;
}

/****************************************************************************/
// writes the array to a file corresponding to the string in fileName
 
void writeComplexArray(zomplex *zomplexArray, long arrayLength, char *fileNameBase) {
  FILE *pf;
  long i;
  char fileName[100];

  //if ()
  sprintf(fileName, "%s.dat", fileNameBase);

  pf = fopen(fileName, "w");
  for (i = 0; i < arrayLength; i++) {
    fprintf(pf, "%ld %g %g\n", i, zomplexArray[i].re, zomplexArray[i].im);
  }
  fclose(pf);

  return;
}

/****************************************************************************/
// prints out current time to pf 
 
void writeCurrentTime(FILE *pf) {
  time_t startTime;

  startTime = time(NULL);
  fprintf(pf, ctime(&startTime));

  return;
}

/****************************************************************************/

void writeSeparation(FILE *pf) {
  fprintf(pf, "\n******************************************************************************\n\n");
  
  return;
}

/****************************************************************************/



