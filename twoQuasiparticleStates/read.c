/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "qp.h"

/*****************************************************************************/
// 

long readInputParFile(grid *rSpaceGrid, dParams *dPar, lParams *lPar) {
  FILE *pf; 
  long i, j, ieof, lenEvalParFile;
  double intExcEnergy, qpEnergy, sigma, *eval, *sige, a;
  char field[100], tmp[100], tmp2[100];
  
  // Set defaults
  dPar->fermiEnergy = -0.175; // eigenstates lower in energy are holes and higher in energy are elecs
  dPar->sigmaCutoff = 0.01;   // sigma=sqrt(<E^2>-<E>^2) -> values lower than sigmaCutoff are eigenstates of H
  lPar->nThreads = 1;
  lPar->calcSinglets = 1; // default to calculating the singlet two quasiparticle states
  lPar->calcTriplets = 0; // default to not calculating the triplet two quasiparticle states
  lPar->calcNonintDensitiesOnly = 0; // set to 1 (true) makes program exit after calculating nonint densities
  lPar->calcIntDensitiesOnly = 0; // set to 1 (true) makes program exit after calculating int densities
  lPar->maxElecStates = 150;
  lPar->maxHoleStates = 150;
  dPar->maxElecDeltaE = 0.2;
  dPar->maxHoleDeltaE = 0.2;
  dPar->maxHoleEnergyToPrint = dPar->minElecEnergyToPrint = dPar->fermiEnergy;
  lPar->maxElecStatesToPrint = lPar->maxElecStates;
  lPar->maxHoleStatesToPrint = lPar->maxHoleStates;
  dPar->minIntTwoQPEnergyToPrint = 0.0;
  lPar->maxIntTwoQPStatesToPrint = 10;
  lPar->maxIntTwoQPStatesExcSize = lPar->maxIntTwoQPStatesToPrint;
  lPar->nFixedQPPoints = 0;
  lPar->nIntTwoQPStates = 0;
  lPar->readInCoulombMatrixElements = 0;
  lPar->readInGrid = 0;
  lPar->readPsiValueByValue = 0;
  lPar->initialStateAverage = 1;
  dPar->electronicTemperature = 298.0;
  dPar->energyLevelSigma = 0.025/AUTOEV;
  dPar->maxEnergyConservation = 0.025/AUTOEV;

  // Read input.par if it exists - exit program otherwise 
  if ( access("input.par", F_OK) != -1 ) {
    pf = fopen("input.par", "r");
    i = 0;
    while (fscanf(pf, "%s", field) != EOF && i < 33) {
      if (! strcmp(field, "calcType")) fscanf(pf, "%s %s", tmp, &(lPar->calcType));
      else if (! strcmp(field, "gridSize")) {
        fscanf(pf, "%s %ld %ld %ld", tmp, &(rSpaceGrid->nGridPointsX), 
                                          &(rSpaceGrid->nGridPointsY), 
                                          &(rSpaceGrid->nGridPointsZ));
      }
      else if (! strcmp(field, "fermiEnergy")) fscanf(pf, "%s %lg", tmp, &(dPar->fermiEnergy));
      else if (! strcmp(field, "sigmaCutoff")) fscanf(pf, "%s %lg", tmp, &(dPar->sigmaCutoff));
      else if (! strcmp(field, "maxElecStates")) fscanf(pf, "%s %ld", tmp, &(lPar->maxElecStates));
      else if (! strcmp(field, "maxHoleStates")) fscanf(pf, "%s %ld", tmp, &(lPar->maxHoleStates));
      else if (! strcmp(field, "maxElecDeltaE")) fscanf(pf, "%s %lg", tmp, &(dPar->maxElecDeltaE));
      else if (! strcmp(field, "maxHoleDeltaE")) fscanf(pf, "%s %lg", tmp, &(dPar->maxHoleDeltaE));
      else if (! strcmp(field, "calcSinglets")) fscanf(pf, "%s %ld", tmp, &(lPar->calcSinglets));
      else if (! strcmp(field, "calcTriplets")) fscanf(pf, "%s %ld", tmp, &(lPar->calcTriplets));
      else if (! strcmp(field, "calcNonintDensitiesOnly")) fscanf(pf, "%s %ld", tmp, &(lPar->calcNonintDensitiesOnly));
      else if (! strcmp(field, "maxHoleEnergyToPrint")) fscanf(pf, "%s %lg", tmp, &(dPar->maxHoleEnergyToPrint));
      else if (! strcmp(field, "minElecEnergyToPrint")) fscanf(pf, "%s %lg", tmp, &(dPar->minElecEnergyToPrint));
      else if (! strcmp(field, "maxHoleStatesToPrint")) fscanf(pf, "%s %ld", tmp, &(lPar->maxHoleStatesToPrint));
      else if (! strcmp(field, "maxElecStatesToPrint")) fscanf(pf, "%s %ld", tmp, &(lPar->maxElecStatesToPrint));
      else if (! strcmp(field, "calcIntDensitiesOnly")) fscanf(pf, "%s %ld", tmp, &(lPar->calcIntDensitiesOnly));
      else if (! strcmp(field, "nIntTwoQPStates")) fscanf(pf, "%s %ld", tmp, &(lPar->nIntTwoQPStates));
      else if (! strcmp(field, "minIntTwoQPEnergyToPrint")) fscanf(pf, "%s %lg", tmp, &(dPar->minIntTwoQPEnergyToPrint));
      else if (! strcmp(field, "maxIntTwoQPStatesToPrint")) fscanf(pf, "%s %ld", tmp, &(lPar->maxIntTwoQPStatesToPrint));
      else if (! strcmp(field, "maxIntTwoQPStatesExcSize")) fscanf(pf, "%s %ld", tmp, &(lPar->maxIntTwoQPStatesExcSize));
      else if (! strcmp(field, "readInCoulombMatrixElements")) fscanf(pf, "%s %ld", tmp, &(lPar->readInCoulombMatrixElements));
      else if (! strcmp(field, "readInGrid")) fscanf(pf, "%s %ld", tmp, &(lPar->readInGrid));
      else if (! strcmp(field, "readPsiValueByValue")) fscanf(pf, "%s %ld", tmp, &(lPar->readPsiValueByValue));
      else if (! strcmp(field, "nFixedQPPoints")) {
        fscanf(pf, "%s %ld", tmp, &(lPar->nFixedQPPoints));
        for (j = 0; j < lPar->nFixedQPPoints; j++) {
          fscanf(pf, "%s %s %lg %lg %lg", tmp, tmp2, &(dPar->fixedQPPoints[j].x), 
                                                     &(dPar->fixedQPPoints[j].y), 
                                                     &(dPar->fixedQPPoints[j].z));
        }
      }
      else if (! strcmp(field, "epsilon")) {
        fscanf(pf, "%s %lg %lg %lg", tmp, &(dPar->epsilon.x), 
                                          &(dPar->epsilon.y), 
                                          &(dPar->epsilon.z));
      }
      else if (! strcmp(field, "electronicTemperature")) fscanf(pf, "%s %lg", tmp, &(dPar->electronicTemperature));
      else if (! strcmp(field, "maxEnergyConservation")) fscanf(pf, "%s %lg", tmp, &(dPar->maxEnergyConservation));
      else if (! strcmp(field, "energyLevelSigma")) fscanf(pf, "%s %lg", tmp, &(dPar->energyLevelSigma));
      else if (! strcmp(field, "nThreads")) fscanf(pf, "%s %ld", tmp, &(lPar->nThreads));
      else if (! strcmp(field, "initialStateAverage")) fscanf(pf, "%s %ld", tmp, &(lPar->initialStateAverage));
      else if (! strcmp(field, "iElecIndex")) fscanf(pf, "%s %ld", tmp, &(lPar->iElecIndex));
      else if (! strcmp(field, "iHoleIndex")) fscanf(pf, "%s %ld", tmp, &(lPar->iHoleIndex));
      else if (! strcmp(field, "fElecIndex")) fscanf(pf, "%s %ld", tmp, &(lPar->fElecIndex));
      else if (! strcmp(field, "fHoleIndex")) fscanf(pf, "%s %ld", tmp, &(lPar->fHoleIndex));
      else {
        printf("Invalid input field and/or format - equal sign required after each field\n");
        printf("Only allowed fields are (case-sensitive):\n\n");
        printf("calcType = elec-hole (required, elec-elec, hole-hole and sp-elec-hole also allowed)\n");
        printf("gridSize = nx ny nz (required, all integers)\n");
        printf("epsilon = eps.x eps.y eps.z (required, dielectric constant in the x y and z directions)\n");
        printf("sigmaCutoff = 0.01 (optional, 0.01 default, criteria to determine if eigenstate or not)\n");
        printf("nThreads = 1 (optional, serial default, number of openmp threads)\n");
        printf("\nThis list is not up to date - see readInputParFile function in read.c to see allowed fields\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
      }
      i++;
    }
    fclose(pf);
  }
  else {
    printf("\n\nNo input.par file detected in current working directory - the program is exiting!!!\n\n");
    fflush(stdout);
    exit(EXIT_FAILURE);
  }

  // Set useful parameter based on input.par file
  rSpaceGrid->nGridPoints = rSpaceGrid->nGridPointsX * rSpaceGrid->nGridPointsY * rSpaceGrid->nGridPointsZ;
  lPar->nGridPoints = rSpaceGrid->nGridPoints; 

  // Determine number of eigenstates if eval.par exits or exit the program if eval.par does not exist in cwd
  lPar->nHoles = lPar->nElecs = 0;
  if ( access("eval.par", F_OK) != -1 ) {
    pf = fopen("eval.par" , "r");
    for (i = ieof = 0; ieof != EOF; i++) {
      ieof = fscanf(pf, "%ld %lg %lg", &j, &qpEnergy, &sigma);
      // TODO: add in maxEnergy in addition to fermiEnergy 
      if (sigma < dPar->sigmaCutoff && qpEnergy < dPar->fermiEnergy) {
        dPar->homoEnergy = qpEnergy;
        lPar->iHomo = i;
        lPar->nHoles++;
      }
      // TODO: add in minEnergy in addition to fermiEnergy 
      else if (sigma < dPar->sigmaCutoff && qpEnergy > dPar->fermiEnergy && ieof != EOF) {
        if (! lPar->nElecs) {
          dPar->lumoEnergy = qpEnergy;
          lPar->iLumo = i;  
        }
        lPar->nElecs++;
      } 
    }
    fclose(pf);
    lenEvalParFile = i-1; // length of eval.par
  }
  else {
    printf("\n\nNo eval.par file detected in current working directory - the program is exiting!!!\n\n");
    fflush(stdout);
    exit(EXIT_FAILURE);
  }
   
  lPar->nHoleEigenstates = lPar->nHoles;
  lPar->nElecEigenstates = lPar->nElecs;
  if (lPar->nHoles > lPar->maxHoleStates) lPar->nHoles = lPar->maxHoleStates;
  if (lPar->nElecs > lPar->maxElecStates) lPar->nElecs = lPar->maxElecStates;

  // TODO: eventually remove this
  lPar->nHolesPlusElecs = lPar->nHoles + lPar->nElecs;
  lPar->nNonintExcitons = lPar->nHoles*lPar->nElecs;
  dPar->fundamentalGap = dPar->lumoEnergy - dPar->homoEnergy;

  // Set useful parameters based on eval.par file and the calcType
  if (! strcmp(lPar->calcType, "elec-hole") || ! strcmp(lPar->calcType, "sp-elec-hole") ||
      ! strcmp(lPar->calcType, "auger-decay-ni") || ! strcmp(lPar->calcType, "auger-decay-i")) {
    dPar->fundamentalGap = dPar->lumoEnergy - dPar->homoEnergy;
    lPar->nHolesPlusElecs = lPar->nHoles + lPar->nElecs;
    lPar->nNonintExcitons = lPar->nHoles*lPar->nElecs;
    lPar->nQPStates = lPar->nHoles + lPar->nElecs;
    lPar->nNonintTwoQPStates = lPar->nHoles*lPar->nElecs;
  }
  else if (! strcmp(lPar->calcType, "elec-elec")) {
    lPar->nQPStates = lPar->nElecs;
    if (lPar->calcSinglets) {
      lPar->nNonintTwoQPStates = (lPar->nElecs*lPar->nElecs + lPar->nElecs)/2;
    }
    else {
      lPar->nNonintTwoQPStates = (lPar->nElecs*lPar->nElecs - lPar->nElecs)/2;
    }
  }
  else if (! strcmp(lPar->calcType, "hole-hole")) {
    lPar->nQPStates = lPar->nHoles;
    if (lPar->calcSinglets) {
      lPar->nNonintTwoQPStates = (lPar->nHoles*lPar->nHoles + lPar->nHoles)/2;
    }
    else {
      lPar->nNonintTwoQPStates = (lPar->nHoles*lPar->nHoles - lPar->nHoles)/2;
    }
  }
  else {
    printf("\n\nInvalid calcType! Allowed calcTypes are: elec-hole, elec-elec, hole-hole and sp-elec-hole\n");
    fflush(stdout);
    exit(EXIT_FAILURE);
  }
  if (! lPar->calcIntDensitiesOnly || ! lPar->nIntTwoQPStates) {
	  lPar->nIntTwoQPStates = lPar->nNonintTwoQPStates;
  }

  // Write the input parameters to stdout
  writeInputParameters(*dPar, *lPar, stdout);
  fflush(stdout);

  return 0;
}

/****************************************************************************/
// 

long readEvalPsiParFiles(double *psiHoles, double *psiElecs, nonintQP *holeQP, nonintQP *elecQP, 
                          double sigmaCutoff, lParams lPar) {
  FILE *pf;
  long i, iGrid, ieof, nStatesInEval, qpIndex, eigenstateList[lPar.nHolesPlusElecs]; 
  long nStoredHoles, nStoredElecs, nStoredStates;
  double qpEnergy, qpSigma, qpSpinZ, *tmpQPEnergy, *tmpQPSigma, *tmpQPSpinZ;

  // Get the number of states in eval.par
  nStatesInEval = 0;
  if ( access("eval.par", F_OK) != -1 ) {
    pf = fopen("eval.par" , "r");
    for (i = ieof = 0; ieof != EOF; i++) {
      ieof = fscanf(pf, "%ld %lg %lg", &qpIndex, &qpEnergy, &qpSigma);
    }
    nStatesInEval = i-1;
    fclose(pf);
  }
  else {
    printf("\n\nNo eval.par file detected in current working directory - the program is exiting!!!\n\n");
    exit(EXIT_FAILURE);
  }
  
  // Dynamically allocate memory
  if ((tmpQPEnergy = (double *) calloc(nStatesInEval, sizeof(double))) == NULL) memoryError("tmpQPEnergy");
  if ((tmpQPSigma = (double *) calloc(nStatesInEval, sizeof(double))) == NULL) memoryError("tmpQPSigma");
  if ((tmpQPSpinZ = (double *) calloc(nStatesInEval, sizeof(double))) == NULL) memoryError("tmpQPSpinZ");

  // Store all the eigenvalues
  pf = fopen("eval.par", "r");
  for (i = 0; i < nStatesInEval; i++) {
    fscanf(pf, "%ld %lg %lg", &qpIndex, &tmpQPEnergy[i], &tmpQPSigma[i]);
  }
  fclose(pf);

  // Store the hole and elec energies in the fields and create a list of their indexes
  nStoredHoles = 0;
  for (i = lPar.iHomo; i >= 0; i--) {
    if (tmpQPSigma[i] < sigmaCutoff && nStoredHoles < lPar.maxHoleStates) {
      holeQP[nStoredHoles].index = nStoredHoles;
      strcpy(holeQP[nStoredHoles].type, "h");
      holeQP[nStoredHoles].energy = tmpQPEnergy[i];
      holeQP[nStoredHoles].sigma = tmpQPSigma[i];
      nStoredHoles++;
      eigenstateList[lPar.nHoles-nStoredHoles] = i;
    }
    if (nStoredHoles == lPar.nHoles) {
      break;
    }
  }
  nStoredElecs = 0;
  for (i = lPar.iLumo; i < nStatesInEval; i++) {
    if (tmpQPSigma[i] < sigmaCutoff && nStoredElecs < lPar.maxElecStates) {
      elecQP[nStoredElecs].index = nStoredElecs;
      strcpy(elecQP[nStoredElecs].type, "e");
      elecQP[nStoredElecs].energy = tmpQPEnergy[i];
      elecQP[nStoredElecs].sigma = tmpQPSigma[i];
      eigenstateList[lPar.nHoles+nStoredElecs] = i;
      nStoredElecs++;
    }
    if (nStoredElecs == lPar.nElecs) {
      break;
    }
  }

  // Store the eigenstates from psi.par in psiHoles and psiElecs
  nStoredStates = 0;
  if ( access("psi.par", F_OK) != -1 ) {
    pf = fopen("psi.par", "r");
    for (i = 0; i < nStatesInEval; i++) {
      if (nStoredStates < lPar.nHoles) {
        if (lPar.readPsiValueByValue) {
          for (iGrid = 0; iGrid < lPar.nGridPoints; iGrid++) {
            fscanf(pf, "%lg", &(psiHoles[(lPar.nHoles-1-nStoredStates)*lPar.nGridPoints + iGrid]));
          }
        }
        else {
          fread(&psiHoles[(lPar.nHoles-1-nStoredStates)*lPar.nGridPoints], sizeof(double), lPar.nGridPoints, pf);
        }
        holeQP[lPar.nHoles-1-nStoredStates].psi = &(psiHoles[(lPar.nHoles-1-nStoredStates)*lPar.nGridPoints]);
      }
      else {
        if (lPar.readPsiValueByValue) {
          for (iGrid = 0; iGrid < lPar.nGridPoints; iGrid++) {
            fscanf(pf, "%lg", &(psiElecs[(nStoredStates-lPar.nHoles)*lPar.nGridPoints + iGrid]));
          }
        }
        else {
          fread(&psiElecs[(nStoredStates-lPar.nHoles)*lPar.nGridPoints], sizeof(double), lPar.nGridPoints, pf);
        }
        elecQP[nStoredStates-lPar.nHoles].psi = &(psiElecs[(nStoredStates-lPar.nHoles)*lPar.nGridPoints]);
      }
      if (i == eigenstateList[nStoredStates]) {
        nStoredStates++;
      }
      if (nStoredStates == lPar.nHolesPlusElecs) {
        break; // as all desired eigenstates have been stored in psiHoles and psiElecs
      }
    }
    fclose(pf);
  }
  else {
    printf("\n\nNo psi.par file detected in current working directory - the program is exiting!!!\n\n");
    exit(EXIT_FAILURE);
  }

  // Read in spin states for the quasiparticle states
  if (! strcmp(lPar.calcType, "sp-elec-hole")) {
    nStoredStates = 0;
    if ( access("spinEval.par", F_OK) != -1 ) {
      pf = fopen("spinEval.par", "r");
      for (i = 0; i < nStatesInEval; i++) {
        fscanf(pf, "%ld %lg %lg %lg", &qpIndex, &tmpQPSpinZ[i], &qpEnergy, &qpSigma);
      }
      fclose(pf);
      for (i = 0; i < nStatesInEval; i++) {
        if (i == eigenstateList[nStoredStates]) {
          if (nStoredStates < lPar.nHoles) {
            holeQP[lPar.nHoles-1-nStoredStates].spinZ = -1.0*tmpQPSpinZ[i];
            holeQP[lPar.nHoles-1-nStoredStates].spin = fabs(tmpQPSpinZ[i]);
          }
          else {
            elecQP[nStoredStates-lPar.nHoles].spinZ = tmpQPSpinZ[i];
            elecQP[nStoredStates-lPar.nHoles].spin = fabs(tmpQPSpinZ[i]);
          }
          nStoredStates++;
          if (nStoredStates == lPar.nHolesPlusElecs) {
            break; // as all desired eigenstates have been stored in psiHoles and psiElecs
          }
        }
      }
    }
    else {
      printf("\n\nNo spinEval.par file detected in current working directory - the program is exiting!!!\n\n");
      exit(EXIT_FAILURE);
    }
  }

  // Write eval.dat
  pf = fopen("holeElecEigenstates.dat", "w");
  for (i = 0; i < lPar.nHoles; i++) writeQuasiparticleState(holeQP[i], pf);
  for (i = 0; i < lPar.nElecs; i++) writeQuasiparticleState(elecQP[i], pf);
  fclose(pf);

  // Sanity check
  if (nStoredStates != (nStoredHoles+nStoredElecs)) {
    printf("\n\nThere is a discripency between the eval.par and psi.par - the program is exiting!!!\n\n");
    exit(EXIT_FAILURE);   
  }

  // Free dynamically allocated memory
  free(tmpQPSpinZ); free(tmpQPSigma); free(tmpQPEnergy); 

  return (nStoredHoles+nStoredElecs);
}

/*****************************************************************************/
//

long readNonIntEigenstates(double *psi, double *eval, double *sige, char *psiFileName, char *evalFileName, 
                          double sigmaCutoff, long nGridPoints) {
  FILE *pfPsi, *pfEval;
  long i, ieofPsi, ieofEval, qpIndex, nEigenstates = 0;
  double qpEnergy, qpSigma, qpPsi[nGridPoints];

  // Make sure that both the wavefunction (psiFileName) and eigenvalue (evalFileName) both exist
  if ( ! (access(evalFileName, F_OK) != -1) ) {
    fprintf(stdout, "\n\nThe progam is existing because %s is not in current working directory!!!\n\n", evalFileName); fflush(stdout);
    exit(EXIT_FAILURE);
  }
  else if ( ! (access(psiFileName, F_OK) != -1) ) {
    fprintf(stdout, "\n\nThe progam is existing because %s is not in current working directory!!!\n\n", psiFileName); fflush(stdout);
    exit(EXIT_FAILURE);
  }
  // Since both files exist we can read in the eigenstates and eigenvalues now
  else {
    pfEval = fopen(evalFileName, "r");
    pfPsi = fopen(psiFileName, "r");
    for (i = ieofEval = 0; ieofEval != EOF; i++) {
      ieofEval = fscanf(pfEval, "%ld %lg %lg", &qpIndex, &qpEnergy, &qpSigma);
      ieofPsi = fread(&qpPsi[0], sizeof(double), nGridPoints, pfPsi);
      if (qpSigma < sigmaCutoff && ieofEval != EOF) {
        eval[nEigenstates] = qpEnergy;
        sige[nEigenstates] = qpSigma;
        memcpy(&(psi[nEigenstates*nGridPoints]), &(qpPsi[0]), nGridPoints*sizeof(double));
        nEigenstates++;
      }
    }
    fclose(pfPsi);
    fclose(pfEval);
  }

  // Return the number of eigenstates read in
  return nEigenstates;
}

/*****************************************************************************/

long readIntTwoQuasiparticleParFiles(intTwoQPState *intTwoQP, double *hMatrix, nonintTwoQPState *nonintTwoQP, 
                                  dParams dPar, lParams lPar) {
  long nReadTwoQPStates = 0;
  long iIntTwoQP, iNonintTwoQP, iR, iS, index, tmpIndex;
  double *tmpCrs, *E;

  // Dynamically allocate memory
  if ((E = (double *) calloc(lPar.nIntTwoQPStates, sizeof(double))) == NULL) memoryError("E");
  if ((tmpCrs = (double *) calloc(lPar.nIntTwoQPStates*lPar.nNonintTwoQPStates, sizeof(double))) == NULL) memoryError("tmpCrs");

  // Read input files
  if (! strcmp(lPar.calcType, "elec-hole") || ! strcmp(lPar.calcType, "auger-decay-i")) {
    readBetheSalpeterResults(tmpCrs, E, lPar.nNonintTwoQPStates, lPar.nIntTwoQPStates);
    // Store reordered tmpCrs in hMatrix
    // Reordering is necessary because BSEcoeff.par is in different order than in the nonintTwoQP structure 
    for (iIntTwoQP = 0; iIntTwoQP < lPar.nIntTwoQPStates; iIntTwoQP++) {
      tmpIndex = iIntTwoQP*lPar.nNonintTwoQPStates;
      for (iNonintTwoQP = 0; iNonintTwoQP < lPar.nNonintTwoQPStates; iNonintTwoQP++) {
        iR = ( lPar.nHoles - 1 - nonintTwoQP[iNonintTwoQP].qp1->index );
        iS = ( nonintTwoQP[iNonintTwoQP].qp2->index );
        hMatrix[tmpIndex + iNonintTwoQP] = tmpCrs[tmpIndex + iS*lPar.nHoles + iR];
      }
    }
  }
  else {
    fprintf(stdout, "readIntTwoQuasiparticleParFiles has only been implemented for calcType = elec-hole\n");
    fflush(stdout);
    exit(EXIT_FAILURE);
  }

  // Store results of the input file reading in intTwoQP
  for (iIntTwoQP = 0; iIntTwoQP < lPar.nIntTwoQPStates; iIntTwoQP++) {
    intTwoQP[iIntTwoQP].index = iIntTwoQP;
    intTwoQP[iIntTwoQP].nNonintTwoQPStates = lPar.nNonintTwoQPStates;
    intTwoQP[iIntTwoQP].energy = E[iIntTwoQP];
    intTwoQP[iIntTwoQP].Crs = &(hMatrix[iIntTwoQP*lPar.nNonintTwoQPStates]);
    intTwoQP[iIntTwoQP].niTwoQP = nonintTwoQP;
    intTwoQP[iIntTwoQP].bindingEnergy = nonintTwoQP[0].energy-intTwoQP[iIntTwoQP].energy;
    intTwoQP[iIntTwoQP].correlationEnergy = 0.0;
    //intTwoQP[iIntTwoQP].correlationEnergy = (nonintTwoQP[0].energy + nonintTwoQP[0].hartreeEnergy +
    //                                        nonintTwoQP[0].exchangeEnergy - intTwoQP[iIntTwoQP].energy);
    nReadTwoQPStates++;
  } 

  // Write interacting two quasiparticle structure to files
  writeIntTwoQPStructuresToFiles(intTwoQP, lPar.nNonintTwoQPStates, lPar.nIntTwoQPStates);

  // Free dynamicallly allocated memory
  free(tmpCrs); free(E);

  return nReadTwoQPStates;
}


/*****************************************************************************/

void readBetheSalpeterResults(double *Cai, double *Ebs, long nNonintExcitons, long nIntExcitons) {
  FILE *pf;
  long iIntExc, iNonintExc, tmpL;
  double tmpD;

  // Read the BSE coefficients from BScoeff.par - exit program if file doesn't exist
  if ( access("BSEcoeff.par", F_OK) != -1 ) { 
    pf = fopen("BSEcoeff.par", "r");
    for (iIntExc = 0; iIntExc < nIntExcitons; iIntExc++) {
      for (iNonintExc = 0; iNonintExc < nNonintExcitons; iNonintExc++) {
        fscanf(pf, "%lg", &(Cai[iNonintExc + iIntExc*nNonintExcitons]));
      }
    } 
    fclose(pf);
  }
  else {
    printf("\n\nNo BSEcoeff.par file detected in current working directory - the program is exiting!!!\n\n");
    exit(EXIT_FAILURE);
  }  

  // Read the interacting excitons energies from exciton.par - exit program if file doesn't exist
  if ( access("exciton.par", F_OK) != -1 ) { 
    pf = fopen("exciton.par", "r");
    for (iIntExc = 0; iIntExc < nIntExcitons; iIntExc++) {
      fscanf(pf, "%ld %lg %lg %lg %lg", &tmpL, &(Ebs[iIntExc]), &tmpD, &tmpD, &tmpD);
    } 
    fclose(pf);
  }
  else {
    printf("\n\nNo exciton.par file detected in current working directory - the program is exiting!!!\n\n");
    exit(EXIT_FAILURE);
  }  

  return;
}

/*****************************************************************************/
//

long readCoulombMatrixElements(coulombMatrixElement *cme, long *rsutIndexToCMEIndexList, 
                                  cmeController VrsutController, char *fileName) {
  FILE *pf;
  long ieof, rsutIndex, iR, iS, iU, iT, lTmp1, lTmp2, nReadInCME = 0;
  double dTmp1, dTmp2, value;

  // Useful constants
  const long nMaxCME = VrsutController.nCME;
  const long nS = VrsutController.nS;
  const long nT = VrsutController.nT;
  const long nU = VrsutController.nU;
  const long nUT = nU*nT; 
  const long nSUT = nS*nUT;

  // Read the interacting excitons energies from exciton.par - exit program if file doesn't exist
  if ( access(fileName, F_OK) != -1 ) { 
    pf = fopen(fileName, "r");
    ieof = 0;
    while (ieof != EOF && nReadInCME < nMaxCME) {
      // Based off of the ordering from Wrsut.dat being printed in coulomb.c
      ieof = fscanf(pf, "%ld %ld %ld %ld %ld %ld %lg %lg %lg", 
                          &iU, &iR, &iT, &iS, &lTmp1, &lTmp2, 
                          &dTmp1, &dTmp2, &value);
      rsutIndex = ( iR*nSUT + iS*nUT + iU*nT + iT );
      if (ieof == EOF) {
        break;
      }
      else if (rsutIndex >= nMaxCME) {
        continue;
      }
      else { // store the read in Coulomb matrix element value
        cme[rsutIndexToCMEIndexList[rsutIndex]].value = value;
        cme[rsutIndexToCMEIndexList[rsutIndex]].alreadyCalculatedFlag = 1;
        nReadInCME++;
      }
    } 
    fclose(pf);
  }
  else {
    printf("\n\nNo %s file detected in current working directory!!!\n\n", fileName);
  }

  return nReadInCME;
}

/*****************************************************************************/
//

long retMatchingCoulombMatrixElementIndex(coulombMatrixElement *coulombME, long nMaxCME, 
                                          long iR, long iS, long iU, long iT) {
  long i;

  for (i = 0; i < nMaxCME; i++) {
    if ( (coulombME[i].iR == iR) && (coulombME[i].iS == iS) && 
                                    (coulombME[i].iU == iU) && 
                                    (coulombME[i].iT == iT)) {
      return (i);
    }
  }

  return (-1);
}

/*****************************************************************************/
//

long readGridParFile(grid *rSpaceGrid, char *gridFileName) {
  FILE *pf;
  long i;
  char field[1000], tmp[1000];

  // Write beginning of function
  writeSeparation(stdout); fprintf(stdout, "Reading in %s file\n\n", gridFileName); fflush(stdout);

  // Read in the grid input parameters from the with name equal to gridFileName 
  if ( access(gridFileName, F_OK) != -1) {
    pf = fopen(gridFileName, "r");
    i = 0;
    while (fscanf(pf, "%s", field) != EOF && i < 1) {
      if (! strcmp(field, "minPos")) {
        fscanf(pf, "%s %lg %lg %lg", tmp, &(rSpaceGrid->minPos.x),
                                          &(rSpaceGrid->minPos.y),
                                          &(rSpaceGrid->minPos.z));
      }
      i++;
    }
    fclose(pf);
  }
  else {
    fprintf(stdout, "\n\nThe progam is existing because %s is not in current working directory!!!\n\n", gridFileName); 
    fflush(stdout); exit(EXIT_FAILURE);
  }

  // Set important rSpaceGrid fields based on what was read in as input
  // Box dimensions
  rSpaceGrid->maxPos.x = -rSpaceGrid->minPos.x;
  rSpaceGrid->maxPos.y = -rSpaceGrid->minPos.y;
  rSpaceGrid->maxPos.z = -rSpaceGrid->minPos.z;
  rSpaceGrid->volume = 8.0*rSpaceGrid->maxPos.x*rSpaceGrid->maxPos.y*rSpaceGrid->maxPos.z;
  // dx, dy, dz, dr and dV=dx*dy*dz -> spacing between the grid points and the volume element
  rSpaceGrid->stepSize.x  = (rSpaceGrid->maxPos.x - rSpaceGrid->minPos.x) / (double)(rSpaceGrid->nGridPointsX);
  rSpaceGrid->stepSize.y  = (rSpaceGrid->maxPos.y - rSpaceGrid->minPos.y) / (double)(rSpaceGrid->nGridPointsY);
  rSpaceGrid->stepSize.z  = (rSpaceGrid->maxPos.z - rSpaceGrid->minPos.z) / (double)(rSpaceGrid->nGridPointsZ);
  rSpaceGrid->stepSize.mag = retVectorMagnitude(rSpaceGrid->stepSize);
  rSpaceGrid->dV = rSpaceGrid->stepSize.x*rSpaceGrid->stepSize.y*rSpaceGrid->stepSize.z;

  // Print grid related statistics
  fprintf(stdout, "Box (quadrant) dimensions: xMax = %.2f ymax = %.2f zMax = %.2f\n", 
            rSpaceGrid->maxPos.x, rSpaceGrid->maxPos.y, rSpaceGrid->maxPos.z);
  fprintf(stdout, "Box volume = %.2f\n", rSpaceGrid->volume); 
  fprintf(stdout, "Grid point spacing: dx = %.6f dy = %.6f dz = %.6f dr = %.6f\n", 
            rSpaceGrid->stepSize.x, rSpaceGrid->stepSize.y, rSpaceGrid->stepSize.z, rSpaceGrid->stepSize.mag);
  fprintf(stdout, "Grid point volume, dV = %.6f\n", rSpaceGrid->dV); fflush(stdout);

  // Write ending of function
  fprintf(stdout, "\nFinished reading in %s\n", gridFileName);   writeSeparation(stdout); fflush(stdout); 

  return 0;
}

/*****************************************************************************/

long readConfFile(double *rx, double *ry, double *rz, atm_st *atm, long nAtoms, FILE *pf) {
  FILE *pw;
  long i, mfermi; 
  double xd, yd, zd;
  
  for (xd = yd = zd = 0.0, mfermi = i = 0; i < nAtoms; i++) {
    fscanf(pf, "%s %lf %lf %lf", atm[i].atyp, &rx[i], &ry[i], &rz[i]);
	atm[i].natyp = assign_atom_number(atm[i].atyp);

    if ((atm[i].natyp == 1) || (atm[i].natyp == 3) || (atm[i].natyp == 4) || (atm[i].natyp == 7)) mfermi++;
    
    xd += rx[i];
    yd += ry[i];
    zd += rz[i];
  }
  xd /= (double)(nAtoms);
  yd /= (double)(nAtoms);
  zd /= (double)(nAtoms);

  for (i = 0; i < nAtoms; i++){
    rx[i] -= xd;
    ry[i] -= yd;
    rz[i] -= zd;
  }  

  pw = fopen("conf.dat" , "w");
  fprintf(pw, "%ld\n", nAtoms);
  for (i = 0; i < nAtoms; i++) {
    fprintf(pw, "%s %lf %lf %lf\n", atm[i].atyp, rx[i], ry[i], rz[i]);
  }
  fclose(pw);
  
  return (4*mfermi);
}

/*****************************************************************************/

void read_pot(double *vr, double *pot, long *nPot, double *dr, atm_st *atm, long n, long ntype) {
  FILE *pf;  
  long i, j, iscan, nloc; 
  double *a, *b;
  char str[100], atype[3];

  if ((a = (double *) calloc(ntype, sizeof(double))) == NULL) memoryError("a");
  if ((b = (double *) calloc(ntype, sizeof(double))) == NULL) memoryError("b");

  a[8] = 0.64;
  a[9] = -0.384;
  a[10] = 0.64;
  a[11] = -0.2;
  b[8] = b[9] = b[10] = b[11] = 2.2287033;

  /*
      Cd = 0
      Se = 1
      In = 2
      As = 3
      Si = 4
      H  = 5
      Zn = 6
      S  = 7
      P1 = 8
      P2 = 9
      P3 = 10
      P4 = 11
      Te = 12
      Cdz = 13
      Sez = 14
      P = 15
      Ga = 16
  */

  for (j = 0; j < ntype*n; j++) pot[j] = vr[j] = 0;
  for (j = 0; j < ntype; j++) nPot[j] = 0;
  
  for (j = 0; j < ntype; j++){
    assign_atom_type(atype,j);
    if ((j==5) || (j==7) || (j==15)) sprintf (str,"pot%c.par",atype[0]);
    else if ((j <= 12) || (j==16)) sprintf (str,"pot%c%c.par",atype[0],atype[1]);
    else sprintf (str,"pot%c%c%c.par",atype[0],atype[1],atype[2]);
    pf = fopen(str , "r");
    if (pf != NULL){
      for (nPot[j] = iscan = 0; iscan != EOF; nPot[j]++)
	iscan = fscanf (pf,"%lg %lg",&vr[j*n+nPot[j]],&pot[j*n+nPot[j]]);
      fclose(pf);
      nPot[j]--;
    }      
    else {
      nPot[j] = nPot[0];
      for (i = 0; i < nPot[j]; i++) {
	vr[j*n+i] = vr[i];
	pot[j*n+i] = (a[j] * exp(-sqr(vr[j*n+i]) / b[j]));
      }
    }

    /*** shift the potentials and get the r-spacing ***/
    for (i = 0; i < nPot[j]; i++) {
      pot[j*n+i] -= pot[j*n+nPot[j]-1];
      dr[j] = vr[j*n+1] - vr[j*n+0];
    }
      
    /*** print shifted pot ***/
    if ((j==5) || (j==7) || (j==15)) sprintf (str,"pot%c.dat",atype[0]);
    else if ((j <= 12) || (j==16)) sprintf (str,"pot%c%c.dat",atype[0],atype[1]);
    else sprintf (str,"pot%c%c%c.dat",atype[0],atype[1],atype[2]);
    //sprintf (str,"pot%s.dat",atype);
    pf = fopen(str , "w");
    for (i = 0; i < nPot[j]; i++) fprintf (pf,"%g %g\n",vr[j*n+i],pot[j*n+i]);
    fclose(pf);
  }
  free(a); free(b);
  return;
}

/*****************************************************************************/

long assign_atom_number(char atyp[3])
{
  char strerror[100];
  
  if ((atyp[0] == 'C') && (atyp[1] == 'd')  && (atyp[2] == '\0')) return(0);
  else if ((atyp[0] == 'S') && (atyp[1] == 'e') && (atyp[2] == '\0')) return(1);
  else if ((atyp[0] == 'I') && (atyp[1] == 'n') && (atyp[2] == '\0')) return(2);
  else if ((atyp[0] == 'A') && (atyp[1] == 's') && (atyp[2] == '\0')) return(3);
  else if ((atyp[0] == 'S') && (atyp[1] == 'i') && (atyp[2] == '\0')) return(4);
  else if ((atyp[0] == 'H') && (atyp[1] == '\0') && (atyp[2] == '\0'))  return(5);
  else if ((atyp[0] == 'Z') && (atyp[1] == 'n') && (atyp[2] == '\0'))  return(6);
  else if ((atyp[0] == 'S') && (atyp[1] == '\0') && (atyp[2] == '\0'))  return(7);
  else if ((atyp[0] == 'S') && (atyp[1] == '\0')) return(7);
  else if ((atyp[0] == 'P') && (atyp[1] == '1') && (atyp[2] == '\0'))  return(8);
  else if ((atyp[0] == 'P') && (atyp[1] == '2') && (atyp[2] == '\0'))  return(9);
  else if ((atyp[0] == 'P') && (atyp[1] == '3') && (atyp[2] == '\0'))  return(10);
  else if ((atyp[0] == 'P') && (atyp[1] == '4') && (atyp[2] == '\0'))  return(11);
  else if ((atyp[0] == 'T') && (atyp[1] == 'e') && (atyp[2] == '\0')) return(12);
  else if ((atyp[0] == 'C') && (atyp[1] == 'd') && (atyp[2] == 'z')) return(13);
  else if ((atyp[0] == 'S') && (atyp[1] == 'e') && (atyp[2] == 'z')) return(14);
  else if ((atyp[0] == 'P') && (atyp[1] == '\0') && (atyp[2] == '\0'))  return(15);
  else if ((atyp[0] == 'P') && (atyp[1] == '\0'))  return(15);
  else if ((atyp[0] == 'G') && (atyp[1] == 'a') && (atyp[2] == '\0'))  return(16);
  else {
    sprintf (strerror,"atom type %s not in current list",atyp);
    memoryError (strerror);
  }
  return(0);
}

/*****************************************************************************/

void assign_atom_type(char *atyp, long j)
{
  if (j == 0) {atyp[0] = 'C'; atyp[1] = 'd'; atyp[2] = '\0';}
  else if (j == 1) {atyp[0] = 'S'; atyp[1] = 'e'; atyp[2] = '\0';}
  else if (j == 2) {atyp[0] = 'I'; atyp[1] = 'n'; atyp[2] = '\0';}
  else if (j == 3) {atyp[0] = 'A'; atyp[1] = 's'; atyp[2] = '\0';}
  else if (j == 4) {atyp[0] = 'S'; atyp[1] = 'i'; atyp[2] = '\0';}
  else if (j == 5) {atyp[0] = 'H'; atyp[1] = '\0'; atyp[2] = '\0';}
  else if (j == 6) {atyp[0] = 'Z'; atyp[1] = 'n'; atyp[2] = '\0';}
  else if (j == 7) {atyp[0] = 'S'; atyp[1] = '\0'; atyp[2] = '\0';}
  else if (j == 8) {atyp[0] = 'P'; atyp[1] = '1'; atyp[2] = '\0';}
  else if (j == 9) {atyp[0] = 'P'; atyp[1] = '2'; atyp[2] = '\0';}
  else if (j == 10) {atyp[0] = 'P'; atyp[1] = '3'; atyp[2] = '\0';}
  else if (j == 11) {atyp[0] = 'P'; atyp[1] = '4'; atyp[2] = '\0';}
  else if (j == 12) {atyp[0] = 'T'; atyp[1] = 'e'; atyp[2] = '\0';}
  else if (j == 13) {atyp[0] = 'C'; atyp[1] = 'd'; atyp[2] = 'z';}
  else if (j == 14) {atyp[0] = 'S'; atyp[1] = 'e'; atyp[2] = 'z';}
  else if (j == 15) {atyp[0] = 'P'; atyp[1] = '\0'; atyp[2] = '\0';}
  else if (j == 16) {atyp[0] = 'G'; atyp[1] = 'a'; atyp[2] = '\0';}
  return;
}

/*****************************************************************************/
