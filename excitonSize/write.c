/****************************************************************************/

// This file does the printing for the program 

/*****************************************************************************/

#include "es.h"

/****************************************************************************/
// Writes the input parameters that were used in AR lifetime calculation

long writeInputParameters(dParams dPar, lParams lPar, FILE *pf) {

  fprintf(pf, "The number of openMP threads used = %ld\n", lPar.nThreads);
  //fprintf(pf, "The number of atoms in the system, nAtoms = %ld\n", lPar.nAtoms);
  //fprintf(pf, "The number of grid points used: nx = %ld  ny = %ld  nz = %ld\n", ist.nx, ist.ny, ist.nz);
  fprintf(pf, "Electronic temperature = %.1f\n", dPar.temp);
  fprintf(pf, "kbT = %.8f gives an energy range of, 3kbT = %.8f\n", dPar.kbT, dPar.boltzEnergyRange);
  fprintf(pf, "Sigma cutoff for eigenstate determination, sigmaCutoff = %.3f\n", dPar.sigmaCutoff);
  fprintf(pf, "Total number of hole eigenstates = %ld\n", lPar.nHoles);
  fprintf(pf, "Total number of electron eigenstates = %ld\n", lPar.nElecs);
  fprintf(pf, "Total number of noninteracting eigenstates = %ld\n", lPar.nHolesPlusElecs);
  fprintf(pf, "The total number of noninteracting excitons = %ld\n", lPar.nNonintExcitons);  
  fprintf(pf, "The index of the HOMO state in eval.dat = %ld\n", lPar.iHomo);
  fprintf(pf, "The index of the LUMO state in eval.dat = %ld\n", lPar.iLumo);
  fprintf(pf, "HOMO energy = %.10f\n", dPar.homoEnergy);
  fprintf(pf, "LUMO energy = %.10f\n", dPar.lumoEnergy);
  fprintf(pf, "Fundamental gap = %.10f %.8f\n", dPar.fundamentalGap, dPar.fundamentalGap*AUTOEV);
  if (lPar.intExcitons) {
    fprintf(pf, "Optical gap = %.10f %.8f\n", dPar.opticalGap, dPar.opticalGap*AUTOEV);
    fprintf(pf, "Exciton binding energy = %.4f\n", (dPar.fundamentalGap-dPar.opticalGap)*1000.0*AUTOEV);    
    fprintf(pf, "The number of holes used in the BSE = %ld\n", lPar.nHoles);
	  fprintf(pf, "The number of elecs used in the BSE = %ld\n", lPar.nElecs);
    fprintf(pf, "The number of interacting excitons  = %ld\n", lPar.nIntExcitons);
  }
  fflush(pf);

  return 0;
}

/****************************************************************************/
// writes

long writeNonintExcitonState(nonintExc niExc, FILE *pf) {
  fprintf(pf, "%ld %ld %ld %.10f %.10f\n", niExc.index, niExc.h->index, niExc.e->index, 
                                            niExc.energy, niExc.energy*AUTOEV);

  return 0;
}

/****************************************************************************/
// writes

long writeQuasiparticleState(nonintQP qp, FILE *pf) {
  fprintf(pf, "%ld %s %.10f %.10f\n", qp.index, qp.type, qp.energy, qp.sigma);

  return 0;
}

/****************************************************************************/
// writes the array to a file corresponding to the string in fileName
 
long writeLongArray(long *longArray, long arrayLength, char *fileNameBase) {
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

  return 0;
}

/****************************************************************************/
// writes the array to a file corresponding to the string in fileName
 
long writeDoubleArray(double *doubleArray, long arrayLength, char *fileNameBase) {
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

  return 0;
}

/****************************************************************************/
// writes the array to a file corresponding to the string in fileName
 
long writeComplexArray(zomplex *zomplexArray, long arrayLength, char *fileNameBase) {
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

  return 0;
}

/****************************************************************************/
// prints out current time to pf 
 
long writeCurrentTime(FILE *pf) {
  time_t startTime;

  startTime = time(NULL);
  fprintf(pf, ctime(&startTime));
  fflush(pf);

  return 0;
}

/****************************************************************************/

long writeSeparation(FILE *pf) {
  fprintf(pf, "\n******************************************************************************\n\n");
  fflush(pf);

  return 0;
}

/****************************************************************************/



