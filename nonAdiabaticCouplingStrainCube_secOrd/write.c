/****************************************************************************/

#include "qp.h"

/****************************************************************************/
/* Writes the input parameters that were used */
/****************************************************************************/
long writeInputParameters(dParams dPar, lParams lPar, FILE *pf) {

  fprintf(pf, "The number of openMP threads used = %ld\n", lPar.nThreads);
  fprintf(pf, "Sigma cutoff for eigenstate determination, sigmaCutoff = %.3f\n", dPar.sigmaCutoff);
  fprintf(pf, "Total number of hole eigenstates = %ld\n", lPar.nHoles);
  fprintf(pf, "Total number of electron eigenstates = %ld\n", lPar.nElecs);
  fprintf(pf, "Total number of noninteracting eigenstates = %ld\n", lPar.nQPStates);
  fprintf(pf, "The total number of noninteracting excitons = %ld\n", lPar.nNonintTwoQPStates);  
  fprintf(pf, "The index of the HOMO state in eval.dat = %ld\n", lPar.iHomo);
  fprintf(pf, "The index of the LUMO state in eval.dat = %ld\n", lPar.iLumo);
  fprintf(pf, "HOMO energy = %.8f %.8f\n", dPar.homoEnergy, dPar.homoEnergy*AUTOEV);
  fprintf(pf, "LUMO energy = %.8f %.8f\n", dPar.lumoEnergy, dPar.lumoEnergy*AUTOEV);
  fprintf(pf, "Fundamental gap = %.10f %.8f\n", dPar.fundamentalGap, dPar.fundamentalGap*AUTOEV);
  fflush(pf);

  return 0;
}

/****************************************************************************/
// writes quasiparticle state
/****************************************************************************/
long writeQuasiparticleState(nonintQP qp, FILE *pf) {
  fprintf(pf, "%ld %s %.10f %.10f\n", qp.index, qp.type, qp.energy, qp.sigma);

  return 0;
}

/****************************************************************************/
/* writes the array to a file corresponding to the string in fileName */
/****************************************************************************/
long writeLongArray(long *longArray, long arrayLength, char *fileNameBase) {
  FILE *pf;
  long i;
  char fileName[100];

  sprintf(fileName, "%s.dat", fileNameBase);

  pf = fopen(fileName, "w");
  for (i = 0; i < arrayLength; i++) {
    fprintf(pf, "%ld %ld\n", i, longArray[i]);
  }
  fclose(pf);

  return 0;
}

/****************************************************************************/
/* writes the array to a file corresponding to the string in fileName */
/****************************************************************************/
long writeDoubleArray(double *doubleArray, long arrayLength, char *fileNameBase) {
  FILE *pf;
  long i;
  char fileName[100];

  sprintf(fileName, "%s.dat", fileNameBase);

  pf = fopen(fileName, "w");
  for (i = 0; i < arrayLength; i++) {
    fprintf(pf, "%ld %g\n", i, doubleArray[i]);
  }
  fclose(pf);

  return 0;
}

/****************************************************************************/
/* writes the array to a file corresponding to the string in fileName */
/****************************************************************************/
long writeComplexArray(zomplex *zomplexArray, long arrayLength, char *fileNameBase) {
  FILE *pf;
  long i;
  char fileName[100];

  sprintf(fileName, "%s.dat", fileNameBase);

  pf = fopen(fileName, "w");
  for (i = 0; i < arrayLength; i++) {
    fprintf(pf, "%ld %g %g\n", i, zomplexArray[i].re, zomplexArray[i].im);
  }
  fclose(pf);

  return 0;
}

/****************************************************************************/
/* prints current time to pf */
/****************************************************************************/
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
