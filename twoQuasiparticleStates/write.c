/****************************************************************************/

// This file does the printing for the program 

/*****************************************************************************/

#include "qp.h"

/****************************************************************************/
// Writes the input parameters that were used in AR lifetime calculation

long writeInputParameters(dParams dPar, lParams lPar, FILE *pf) {

  fprintf(pf, "The number of openMP threads used = %ld\n", lPar.nThreads);
  fprintf(pf, "Fermi energy, fermiEnergy = %.3f\n", dPar.fermiEnergy);
  fprintf(pf, "Sigma cutoff for eigenstate determination, sigmaCutoff = %.3f\n", dPar.sigmaCutoff);
  fprintf(pf, "Total number of hole eigenstates = %ld\n", lPar.nHoleEigenstates);
  fprintf(pf, "Total number of electron eigenstates = %ld\n", lPar.nElecEigenstates); 
  fprintf(pf, "Total number of hole eigenstates used = %ld\n", lPar.nHoles);
  fprintf(pf, "Total number of electron eigenstates used = %ld\n", lPar.nElecs);
  fprintf(pf, "Total number of noninteracting eigenstates used = %ld\n", lPar.nQPStates);
  fprintf(pf, "The total number of noninteracting two quasiparticle states used = %ld\n", lPar.nNonintTwoQPStates);  
  fprintf(pf, "The index of the HOMO state in eval.dat = %ld\n", lPar.iHomo);
  fprintf(pf, "The index of the LUMO state in eval.dat = %ld\n", lPar.iLumo);
  fprintf(pf, "HOMO energy = %.8f %.8f\n", dPar.homoEnergy, dPar.homoEnergy*AUTOEV);
  fprintf(pf, "LUMO energy = %.8f %.8f\n", dPar.lumoEnergy, dPar.lumoEnergy*AUTOEV);
  fprintf(pf, "Fundamental gap = %.10f %.8f\n", dPar.fundamentalGap, dPar.fundamentalGap*AUTOEV);
  fprintf(pf, "Electronic temperature = %.1f\n", dPar.electronicTemperature);
  if (lPar.calcNonintDensitiesOnly) {
    fprintf(pf, "\nOnly calculating the noninteracting quasiparticle densities\n\n");
  }
  else if (lPar.calcIntDensitiesOnly) {
    fprintf(pf, "The number of interacting two quasiparticle states to be read in = %ld\n", lPar.nIntTwoQPStates);
  }
  else {
     fprintf(pf, "The number of interacting two quasiparticle states to be calculated = %ld\n", lPar.nIntTwoQPStates);
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
// Writes a vector property for the specified number of noninteracting two
// quasiparticle states

long writeNonintTwoQPVectorProperty(vector *vObservable, nonintTwoQPState *nonintTwoQP, 
                                    long nNonintTwoQPStates, char *fileName) {
  FILE *pf;
  long iNonintTwoQP;

  pf = fopen(fileName, "w");
  fprintf(pf, "#iTwoQP iH  iE   eTwoQP      |O|^2        <Ox>        <Oy>        <Oz>\n");
  for (iNonintTwoQP = 0; iNonintTwoQP < nNonintTwoQPStates; iNonintTwoQP++) {
    fprintf(pf, "%5ld %3ld %3ld % .8f % .8f % .8f % .8f % .8f\n", iNonintTwoQP, nonintTwoQP[iNonintTwoQP].qp1->index,
        nonintTwoQP[iNonintTwoQP].qp2->index, nonintTwoQP[iNonintTwoQP].energy, sqr(vObservable[iNonintTwoQP].mag), 
        vObservable[iNonintTwoQP].x, vObservable[iNonintTwoQP].y, vObservable[iNonintTwoQP].z);
  }
  fclose(pf);

  return 0;
}

/****************************************************************************/
// Writes a vector property for the specified number of noninteracting two
// quasiparticle states including the spin quantum numbers of the two quasiparticle state

long writeSpinPolarizedNonintTwoQPVectorProperty(vector *vObservable, nonintTwoQPState *nonintTwoQP, 
                                                 long nNonintTwoQPStates, char *fileName) {
  FILE *pf;
  long iNonintTwoQP;

  pf = fopen(fileName, "w");
  fprintf(pf, "#iTwoQP iH  iE   eTwoQP     S   Sz   Sz1  Sz2    |O|^2        <Ox>        <Oy>        <Oz>\n");
  for (iNonintTwoQP = 0; iNonintTwoQP < nNonintTwoQPStates; iNonintTwoQP++) {
    fprintf(pf, "%5ld %3ld %3ld % .8f % .1f % .1f % .1f % .1f % .8f % .8f % .8f % .8f\n", iNonintTwoQP, nonintTwoQP[iNonintTwoQP].qp1->index,
        nonintTwoQP[iNonintTwoQP].qp2->index, nonintTwoQP[iNonintTwoQP].energy, nonintTwoQP[iNonintTwoQP].spin,   
        nonintTwoQP[iNonintTwoQP].spinZ, nonintTwoQP[iNonintTwoQP].qp1->spinZ, nonintTwoQP[iNonintTwoQP].qp2->spinZ,
        sqr(vObservable[iNonintTwoQP].mag), vObservable[iNonintTwoQP].x, vObservable[iNonintTwoQP].y, vObservable[iNonintTwoQP].z);
  }
  fclose(pf);

  return 0;
}

/****************************************************************************/
// Writes a vector property for the specified number of 
// interacting two quasiparticle states

long writeIntTwoQPVectorProperty(vector *vObservable, intTwoQPState *intTwoQP, 
                                    long nIntTwoQPStates, char *fileName) {
  FILE *pf;
  long iIntTwoQP;

  pf = fopen(fileName, "w");
  fprintf(pf, "#iTwoQP   eTwoQP      |O|^2        <Ox>        <Oy>        <Oz>\n");
  for (iIntTwoQP = 0; iIntTwoQP < nIntTwoQPStates; iIntTwoQP++) {
    fprintf(pf, "%5ld  % .8f % .8f % .8f % .8f % .8f\n", iIntTwoQP, 
                intTwoQP[iIntTwoQP].energy, sqr(vObservable[iIntTwoQP].mag), 
                vObservable[iIntTwoQP].x, vObservable[iIntTwoQP].y, vObservable[iIntTwoQP].z);
  }
  fclose(pf);

  return 0;
}

/****************************************************************************/
// Writes a vector property for the specified number of interacting 
// two quasiparticle states including the spin quantum numbers of the two quasiparticle state

long writeSpinPolarizedIntTwoQPVectorProperty(vector *vObservable, intTwoQPState *intTwoQP, 
                                              long nIntTwoQPStates, char *fileName) {
  FILE *pf;
  long iIntTwoQP;

  pf = fopen(fileName, "w");
  fprintf(pf, "#iTwoQP   eTwoQP     S   Sz     |O|^2        <Ox>        <Oy>        <Oz>\n");
  for (iIntTwoQP = 0; iIntTwoQP < nIntTwoQPStates; iIntTwoQP++) {
    // TODO: remove setting these to 0 and calculate them after diagonalization
    intTwoQP[iIntTwoQP].spin  = 0.0;
    intTwoQP[iIntTwoQP].spinZ = 0.0;
    fprintf(pf, "%5ld  % .8f % .1f % .1f % .8f % .8f % .8f % .8f\n", iIntTwoQP, 
                intTwoQP[iIntTwoQP].energy, intTwoQP[iIntTwoQP].spin, intTwoQP[iIntTwoQP].spinZ,
                sqr(vObservable[iIntTwoQP].mag), vObservable[iIntTwoQP].x, 
                vObservable[iIntTwoQP].y, vObservable[iIntTwoQP].z);
  }
  fclose(pf);

  return 0;
}

/****************************************************************************/
//

long writeIntTwoQPStructuresToFiles(intTwoQPState *intTwoQP, long nNonintTwoQPStates, long nIntTwoQPStates) {
  FILE *pf;
  long iNonintTwoQP, iIntTwoQP;

  // Print energy related fields of the interacting two quasiparticle states
  pf = fopen("intTwoParticleEnergies.dat", "w");
  for (iIntTwoQP = 0; iIntTwoQP < nIntTwoQPStates; iIntTwoQP++) {
    fprintf(pf, "%4ld %.16f % .16f % .16f\n", iIntTwoQP, intTwoQP[iIntTwoQP].energy, intTwoQP[iIntTwoQP].bindingEnergy,
                                             intTwoQP[iIntTwoQP].correlationEnergy);
  }
  fclose(pf);

  // Print out the coefficients of the eigenstates
  pf = fopen("intTwoParticleCoefficients.dat", "w");
  for (iIntTwoQP = 0; iIntTwoQP < nIntTwoQPStates; iIntTwoQP++) {
    for (iNonintTwoQP = 0; iNonintTwoQP < nNonintTwoQPStates; iNonintTwoQP++) {
     fprintf(pf, "%ld %.16f\n", iNonintTwoQP, intTwoQP[iIntTwoQP].Crs[iNonintTwoQP]);
    }
  }
  fclose(pf);

  return 0;
}

/****************************************************************************/
// 

void writeCubeFile(double *rho, grid rSpaceGrid, char *fileName) {
  FILE *pf, *pConfFile;
  long iGrid, iX, iY, iZ, iYZ, iXY, nAtoms, atomType;
  double x, y, z;
  char line[80], atomSymbol[10];


  pConfFile = fopen("conf.par", "r");
  fscanf(pConfFile, "%ld", &nAtoms);
  pf = fopen(fileName, "w");
  fprintf(pf, "CUBE FILE\n");
  //fprintf(pf, "OUTER LOOP: Z, MIDDLE LOOP: Y, INNER LOOP: X\n");
  fprintf(pf, "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n");
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", nAtoms, rSpaceGrid.minPos.x, rSpaceGrid.minPos.y, rSpaceGrid.minPos.z);
  //fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", nAtoms, 0.0, 0.0, 0.0);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", rSpaceGrid.nGridPointsX, rSpaceGrid.stepSize.x, 0.0, 0.0);
  //fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", rSpaceGrid.nGridPointsZ, 0.0, 0.0, rSpaceGrid.stepSize.z);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", rSpaceGrid.nGridPointsY, 0.0, rSpaceGrid.stepSize.y, 0.0);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", rSpaceGrid.nGridPointsZ, 0.0, 0.0, rSpaceGrid.stepSize.z);
  //fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", rSpaceGrid.nGridPointsX, rSpaceGrid.stepSize.x, 0.0, 0.0);
  fgets(line, 80, pConfFile); 
  while(fgets(line, 80, pConfFile) != NULL) {
    sscanf(line, "%2s %lf %lf %lf", &atomSymbol, &x, &y, &z);
    if (! strcmp(atomSymbol, "Cd")) { 
      atomType = 48;
    }
    else if (! strcmp(atomSymbol, "S")) {  
      atomType = 16;
    }
    else if (! strcmp(atomSymbol, "Se")) { 
      atomType = 34;
    }
    else if (! strcmp(atomSymbol, "Zn")) {
      atomType = 30;
    }
    else if (! strcmp(atomSymbol, "Te")) {
      atomType = 52;
    }
	else if (! strcmp(atomSymbol, "C")) {
	  atomType = 6;
	}
    else if (! strcmp(atomSymbol, "Si")) {
	  atomType = 14;
    }
    else if (! strcmp(atomSymbol, "P")) {
	  atomType = 15;
    }
    else if (! strcmp(atomSymbol, "As")) {
	  atomType = 33;
    }
    else if (! strcmp(atomSymbol, "In")) {
	  atomType = 49;
    }
    else if (! strcmp(atomSymbol, "Ga")) {
    atomType = 31;
    }
    else { 
      atomType = 1; 
    }
	fprintf(pf, "%5ld%12.6f%12.6f%12.6f%12.6f\n", atomType, 0.0, x, y, z);
	//fprintf(pf, "%5i%12.6f%12.6f%12.6f%12.6f\n", atomType, 0.0, x-rSpaceGrid.minPos.x, 
	//		y-rSpaceGrid.minPos.y, z-rSpaceGrid.minPos.z);
  }

  for (iX = 0; iX < rSpaceGrid.nGridPointsX; iX++) {
    for (iY = 0; iY < rSpaceGrid.nGridPointsY; iY++) {
      iXY = rSpaceGrid.nGridPointsZ * (rSpaceGrid.nGridPointsY * iX + iY);
	    for (iZ = 0; iZ < rSpaceGrid.nGridPointsZ; iZ++) {
        iGrid = iXY + iZ;
        fprintf(pf, "%g ", rho[iGrid]);
        if (iX % 6 == 5) {
          fprintf(pf, "\n");
        }
      }
      fprintf(pf, "\n");
    }
  }
  fclose(pConfFile);
  fclose(pf);

  return;
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



