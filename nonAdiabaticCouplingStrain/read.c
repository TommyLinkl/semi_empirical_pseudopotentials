/*****************************************************************************/

#include "qp.h"

/*****************************************************************************/
long readInputParFile(grid3d *rSpaceGrid, dParams *dPar, lParams *lPar) {
  FILE *pf; 
  long i, j, ieof, lenEvalParFile;
  double intExcEnergy, qpEnergy, sigma, *eval, *sige, a;
  char field[100], tmp[100];
  
  // Set defaults
  lPar->nThreads = 1;
  lPar->nPhonons = 1;
  lPar->maxElecStates = 1000;
  lPar->maxHoleStates = 1000;
  dPar->maxElecDeltaE = 0.2;
  dPar->maxHoleDeltaE = 0.2;
  dPar->sigmaCutoff = 0.01;   // sigma=sqrt(<E^2>-<E>^2) -> values lower than sigmaCutoff are eigenstates of H
  dPar->fermiEnergy = -0.175; // eigenstates lower in energy aree holes and higheer in energy are elecs

  // Most likely not needed: May, 13th 2019
  lPar->calcSinglets = 1; // default to calculating the singlet two quasiparticle states
  lPar->calcTriplets = 0; // default to not calculating the triplet two quasiparticle states
  dPar->epsilon.x = dPar->epsilon.y = dPar->epsilon.z = 1.0;

  // Read input.par if it exists - exit program otherwise 
  if ( access("input.par", F_OK) != -1 ) {
    pf = fopen("input.par", "r");
    i = 0;
    while (fscanf(pf, "%s", field) != EOF && i < 9) {
      if (! strcmp(field, "calcType")) fscanf(pf, "%s %s", tmp, &(lPar->calcType));
      else if (! strcmp(field, "gridSize")) {
        fscanf(pf, "%s %ld %ld %ld", tmp, &(rSpaceGrid->nGridPointsX), 
                                          &(rSpaceGrid->nGridPointsY), 
                                          &(rSpaceGrid->nGridPointsZ));
      }
      else if (! strcmp(field, "stateType")) fscanf(pf, "%s %s", tmp, &(lPar->stateType));
      else if (! strcmp(field, "maxElecStates")) fscanf(pf, "%s %ld", tmp, &(lPar->maxElecStates));
      else if (! strcmp(field, "maxHoleStates")) fscanf(pf, "%s %ld", tmp, &(lPar->maxHoleStates));
      else if (! strcmp(field, "maxElecDeltaE")) fscanf(pf, "%s %lg", tmp, &(dPar->maxElecDeltaE));
      else if (! strcmp(field, "maxHoleDeltaE")) fscanf(pf, "%s %lg", tmp, &(dPar->maxHoleDeltaE));
      else if (! strcmp(field, "sigmaCutoff")) fscanf(pf, "%s %lg", tmp, &(dPar->sigmaCutoff));
      else if (! strcmp(field, "fermiEnergy")) fscanf(pf, "%s %lg", tmp, &(dPar->fermiEnergy));
      else if (! strcmp(field, "nPhonons")) fscanf(pf, "%s %ld", tmp, &(lPar->nPhonons));
      else if (! strcmp(field, "nThreads")) fscanf(pf, "%s %ld", tmp, &(lPar->nThreads));
      else if (! strcmp(field, "crystalStructure")) fscanf(pf, "%s %s", tmp, &(lPar->crystalStructure));
      else {
        printf("Invalid input field and/or format - equal sign required after each field\n");
        printf("Only allowed fields are (case-sensitive):\n\n");
        printf("calcType = elec-hole (required, elec-elec and hole-hole also allowed)\n");
        printf("gridSize = nx ny nz (required, all integers)\n");
        printf("epsilon = eps.x eps.y eps.z (required, dielectric constant in the x y and z directions)\n");
        printf("sigmaCutoff = 0.01 (optional, 0.01 default, criteria to determine if eigenstate or not)\n");
        printf("nThreads = 1 (optional, serial default, number of openmp threads)\n");
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

  // 
  // Determine the number of atoms
  if ( access("conf.par", F_OK) != -1 ) {  
    pf = fopen("conf.par" , "r");
    fscanf(pf, "%ld", &(lPar->nAtoms));
    fclose(pf);
  }
  else {
    printf("\n\nNo conf.par file detected in current working directory - the program is exiting!!!\n\n");
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
      if (sigma < dPar->sigmaCutoff && qpEnergy < dPar->fermiEnergy) {
        dPar->homoEnergy = qpEnergy;
        lPar->iHomo = i;
        lPar->nHoles++;
      }
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
   
  // TODO: do this after writeInputParameters so the total number of eigenstates is correct 
  if (lPar->nHoles > lPar->maxHoleStates) lPar->nHoles = lPar->maxHoleStates;
  if (lPar->nElecs > lPar->maxElecStates) lPar->nElecs = lPar->maxElecStates;

  // TODO: eventually remove this
  lPar->nHolesPlusElecs = lPar->nHoles + lPar->nElecs;
  lPar->nNonintExcitons = lPar->nHoles*lPar->nElecs;
  dPar->fundamentalGap = dPar->lumoEnergy - dPar->homoEnergy;

  // Set useful parameters based on eval.par file and the calcType
  if (! strcmp(lPar->calcType, "nonInt")) {
    dPar->fundamentalGap = dPar->lumoEnergy - dPar->homoEnergy;
    lPar->nHolesPlusElecs = lPar->nHoles + lPar->nElecs;
    lPar->nQPStates = lPar->nHoles + lPar->nElecs;
    lPar->nNonintExcitons = lPar->nHoles*lPar->nElecs;
    lPar->nNonintTwoQPStates = lPar->nHoles*lPar->nElecs;
  }
  // TODO: eventually remove this, May 13th 2019
  else if (! strcmp(lPar->calcType, "elec-hole")) {
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
    printf("\n\nInvalid calcType! Allowed calcTypes are: nonInt, elec-hole, elec-elec and hole-hole\n");
    fflush(stdout);
    exit(EXIT_FAILURE);
  }
  lPar->nIntTwoQPStates = lPar->nNonintTwoQPStates;

  // Write the input parameters to stdout
  writeInputParameters(*dPar, *lPar, stdout);
  fflush(stdout);

  return 0;
}

/****************************************************************************/
long fillAtomStructureFromConfFile(atom *atoms, lParams *lPar) {
  FILE *pf;
  long iAtom, tmp;
  double xd, yd, zd;

  lPar->nAtomTypes = 0;
  lPar->nSCAtomTypes = 0;
  lPar->nSCAtoms = 0;
  xd = yd = zd = 0.;

  fprintf(stdout, "Reading in the atomic positions\n"); fflush(stdout);

  if ( access("conf.par", F_OK) != -1 ) {  
    pf = fopen("conf.par" , "r");
    fscanf(pf, "%ld", &(tmp));

    for (iAtom = 0; iAtom < lPar->nAtoms; iAtom++) {
      atoms[iAtom].index = iAtom;

      fscanf(pf, "%s", atoms[iAtom].symbol);
      fscanf(pf, "%lg %lg %lg", &atoms[iAtom].pos.x, &atoms[iAtom].pos.y, &atoms[iAtom].pos.z);
      atoms[iAtom].pos.mag = retVectorMagnitude(atoms[iAtom].pos);
      xd += atoms[iAtom].pos.x;
      yd += atoms[iAtom].pos.y;
      zd += atoms[iAtom].pos.z;

      if (isNewAtomType(atoms, iAtom)) {
        atoms[iAtom].type = lPar->nAtomTypes; // gives new integer for each new atom type
        lPar->nAtomTypes++;
        if (! isAPassivationSymbol(atoms[iAtom].symbol)) {
          atoms[iAtom].scType = lPar->nSCAtomTypes;
          lPar->nSCAtomTypes++;
        }
      }
      // fprintf(stdout, "hi %ld %ld\n", iAtom, atoms[iAtom].type);
      if (! isAPassivationSymbol(atoms[iAtom].symbol)) {
        lPar->nSCAtoms++;
      }
    }
    fclose(pf);
  }
  else {
    exit(EXIT_FAILURE);
  }

  xd /= (double)(lPar->nAtoms);
  yd /= (double)(lPar->nAtoms);
  zd /= (double)(lPar->nAtoms);
  for (iAtom = 0; iAtom < lPar->nAtoms; iAtom++) {
      atoms[iAtom].pos.x -= xd;
      atoms[iAtom].pos.y -= yd;
      atoms[iAtom].pos.z -= zd;
  }

  fprintf(stdout, "Number of atoms = %ld\n", lPar->nAtoms);
  fprintf(stdout, "Number of atom types = %ld\n", lPar->nAtomTypes);
  fprintf(stdout, "Number of SC atoms = %ld\n", lPar->nSCAtoms);
  fprintf(stdout, "Number of SC atom types = %ld\n", lPar->nSCAtomTypes);
  fflush(stdout);

  pf = fopen("Conf.dat", "w");
  for (iAtom = 0; iAtom < lPar->nAtoms; iAtom++) {
    fprintf(pf, "%s %ld %ld %lg %lg %lg\n", atoms[iAtom].symbol, atoms[iAtom].index, atoms[iAtom].type,
                    atoms[iAtom].pos.x, atoms[iAtom].pos.y, atoms[iAtom].pos.z);
  }
  fclose(pf);
    
  return 0;
}

/****************************************************************************/
// Return 1 (True) if the atomSymbol signifies a passivation ligand (P1, P2, P3 or P4)
// and return 0 (False) otherwise
/****************************************************************************/
long isAPassivationSymbol(char *atomSymbol) {
  // check if a the symbol is P1, P2, P3 or P4 -> return 1/True if so
  if (! strcmp(atomSymbol, "P1") || ! strcmp(atomSymbol, "P2") ||
    ! strcmp(atomSymbol, "P3") || ! strcmp(atomSymbol, "P4")) return 1;
  else return 0; // not a passivation ligand symbol -> return 0/False
}

/****************************************************************************/
// returns 0/ false if atom symbol is already in the list or 1/true if 
// it is a new symbol in atoms[currIndex].symbol
/****************************************************************************/
long isNewAtomType(atom *atoms, long currIndex) {
  long i;

  for (i = 0; i < currIndex; i++) {
    if (! strcmp(atoms[i].symbol, atoms[currIndex].symbol)) {
      atoms[currIndex].type = atoms[i].type; // make types equal since same symbols
      return 0; // not a new atom type
    }
  }

  return 1; // new atom type
}

/****************************************************************************/
void readAtomicPseudopotentials(grid1d *atomicPPGrid, gridPoint1d *atomicPPGP, atom *atoms, lParams lPar, double *a4Params) {
  FILE *pf;
  long iAtom, iAtomType, iGrid;
  long nGridPoints = atomicPPGrid->nGridPoints;
  char fileName[100];
  char fileName2[100];

  fprintf(stdout, "Reading in the atomic pseudopotentials\n"); fflush(stdout);

  // r-space grid point will be indexed by 0
  atomicPPGrid->index = 1;
  atomicPPGrid->gP = atomicPPGP;

  for (iAtomType = 0; iAtomType < lPar.nSCAtomTypes; iAtomType++) {
    for (iAtom = 0; iAtom < lPar.nAtoms; iAtom++) {
      if (atoms[iAtom].type == iAtomType && ! isAPassivationSymbol(atoms[iAtom].symbol)) {
        memset(&fileName[0], 0, sizeof(fileName));
        strcpy(fileName, "pot");
        strcat(fileName, atoms[iAtom].symbol);
        strcat(fileName, ".par");
        fprintf(stdout, "Reading in %s file\n", fileName);
        pf = fopen(fileName, "r");
        for (iGrid = 0; iGrid < nGridPoints; iGrid++) {
          atomicPPGP[iAtomType*nGridPoints + iGrid].index = iGrid;
          fscanf(pf, "%lg %lg", &(atomicPPGP[iAtomType*nGridPoints + iGrid].pos), 
                                &(atomicPPGP[iAtomType*nGridPoints + iGrid].localVal)); 
        }
        fclose(pf);
        // read in a4 parameter
        memset(&fileName2[0], 0, sizeof(fileName2));
        strcpy(fileName2, "pot");
        strcat(fileName2, atoms[iAtom].symbol);
        strcat(fileName2, "_a4.par");
        fprintf(stdout, "Reading in %s file\n", fileName2);
        pf = fopen(fileName2, "r");
        if (pf != NULL) {
            fscanf(pf, "%lg", &a4Params[iAtomType]);
            fclose(pf);
        }
        else {
            a4Params[iAtomType] = 0.;
            printf("Warning: no %s file... setting a4 to 0.\n", fileName2);
        }
        break;
      }
      else if (isAPassivationSymbol(atoms[iAtom].symbol)) {
        iAtomType--;
        break;
      }
    }
  }

  atomicPPGrid->stepSize = (atomicPPGP[1].pos - atomicPPGP[0].pos);

  fprintf(stdout, "Finished reading in the atomic pseudopotentials\n"); fflush(stdout);
  
  return;
}

/****************************************************************************/
long readEvalPsiParFiles(double *psiHoles, double *psiElecs, nonintQP *holeQP, nonintQP *elecQP, 
                          double sigmaCutoff, lParams lPar) {
  FILE *pf;
  long i, ieof, nStatesInEval, qpIndex, eigenstateList[lPar.nHolesPlusElecs]; 
  long nStoredHoles, nStoredElecs, nStoredStates;
  double qpEnergy, qpSigma, *tmpQPEnergy, *tmpQPSigma;

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

  // Store all the eigenvalues
  pf = fopen("eval.par", "r");
  for (i = 0; i < nStatesInEval; i++) {
    fscanf(pf, "%ld %lg %lg", &qpIndex, &tmpQPEnergy[i], &tmpQPSigma[i]);
  }
  fclose(pf);

  // Store the hole and elec energies in evalbe and create a list of their indexes
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

  // Store the eigenstates from psi.par in psibe
  nStoredStates = 0;
  if ( access("psi.par", F_OK) != -1 ) {
    pf = fopen("psi.par" , "r");
    for (i = 0; i < nStatesInEval; i++) {
      if (nStoredStates < lPar.nHoles) {
        fread(&psiHoles[(lPar.nHoles-1-nStoredStates)*lPar.nGridPoints], sizeof(double), lPar.nGridPoints, pf);
        holeQP[lPar.nHoles-1-nStoredStates].psi = &psiHoles[(lPar.nHoles-1-nStoredStates)*lPar.nGridPoints];
      }
      else {
        fread(&psiElecs[(nStoredStates-lPar.nHoles)*lPar.nGridPoints], sizeof(double), lPar.nGridPoints, pf);
        elecQP[nStoredStates-lPar.nHoles].psi = &psiElecs[(nStoredStates-lPar.nHoles)*lPar.nGridPoints];
      }
      if (i == eigenstateList[nStoredStates]) {
        nStoredStates++;
      }
      if (nStoredStates == lPar.nHolesPlusElecs) {
        break; // as all desired eigenstates have been stored in psibe
      }
    }
    fclose(pf);
  }
  else {
    printf("\n\nNo psi.par file detected in current working directory - the program is exiting!!!\n\n");
    exit(EXIT_FAILURE);
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
  free(tmpQPSigma); free(tmpQPEnergy);

  return (nStoredHoles+nStoredElecs);
}

/*****************************************************************************/
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
void readBetheSalpeterResults(double *Cai, double *Ebs, lParams ist) {
  FILE *pf;
  long iIntExc, iNonintExc, tmpL;
  double tmpD;

  // Read the BSE coefficients from BScoeff.par - exit program if file doesn't exist
  if ( access("BSEcoeff.par", F_OK) != -1 ) { 
    pf = fopen("BSEcoeff.par", "r");
    for (iIntExc = 0; iIntExc < ist.nIntExcitons; iIntExc++) {
      for (iNonintExc = 0; iNonintExc < ist.nNonintExcitons; iNonintExc++) {
        fscanf(pf, "%lg", &(Cai[iNonintExc + iIntExc*ist.nNonintExcitons]));
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
    for (iIntExc = 0; iIntExc < ist.nIntExcitons; iIntExc++) {
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
  for (i = 0; i < nAtoms; i++) {
    fprintf(pw, "%s %g %g %g %ld\n", atm[i].atyp, rx[i], ry[i], rz[i], atm[i].natyp);
  }
  fclose(pw);
  
  return (4*mfermi);
}

/*****************************************************************************/
/* read all neighbors of each atom */
void readNearestNeighbors(long nAtoms, atom *atomNeighbors) {

    FILE *pf;
    long iAtom;
    char at[4];
    int tmpi;
    double tmpx, tmpy, tmpz;

    if ( access("allNeighborBonds.par", F_OK) != -1) {

        pf = fopen("allNeighborBonds.par", "r");
        for (iAtom = 0; iAtom < nAtoms; iAtom++) {
            fscanf(pf, "%s %i %lg %lg %lg %s %i %lg %lg %lg %s %i %lg %lg %lg %s %i %lg %lg %lg %s %i %lg %lg %lg",
                    at, &tmpi, &tmpx, &tmpy, &tmpz,
                    atomNeighbors[4*iAtom].symbol, &atomNeighbors[4*iAtom].index,
                    &atomNeighbors[4*iAtom].pos.x, &atomNeighbors[4*iAtom].pos.y, &atomNeighbors[4*iAtom].pos.z,
                    atomNeighbors[4*iAtom+1].symbol, &atomNeighbors[4*iAtom+1].index,
                    &atomNeighbors[4*iAtom+1].pos.x, &atomNeighbors[4*iAtom+1].pos.y, &atomNeighbors[4*iAtom+1].pos.z,
                    atomNeighbors[4*iAtom+2].symbol, &atomNeighbors[4*iAtom+2].index,
                    &atomNeighbors[4*iAtom+2].pos.x, &atomNeighbors[4*iAtom+2].pos.y, &atomNeighbors[4*iAtom+2].pos.z,
                    atomNeighbors[4*iAtom+3].symbol, &atomNeighbors[4*iAtom+3].index,
                    &atomNeighbors[4*iAtom+3].pos.x, &atomNeighbors[4*iAtom+3].pos.y, &atomNeighbors[4*iAtom+3].pos.z);
        }
    }
    else {
        printf("\n\nNo allNeighborBonds.par file in current working directory -- the program is exiting!!!\n\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }

    return;

}

/*****************************************************************************/
/* calculate reference tetrahedron volume for each atom */
void calculateRefTetrahedronVol(long nAtoms, int crystalStructure, atom *atoms, atom *atomNeighbors, double *refTetrahedronVol) {

    long iAtom;
    int numSe;

    for (iAtom = 0; iAtom < nAtoms; iAtom++) {

        if (isAPassivationSymbol(atoms[iAtom].symbol)) {
            refTetrahedronVol[iAtom] = 0.;
        }
        else {
            numSe = 0;
            // count number of Se neighbors
            if (! strcmp(atomNeighbors[4*iAtom].symbol, "Se")) numSe++;
            if (! strcmp(atomNeighbors[4*iAtom+1].symbol, "Se")) {
                numSe++;
                if (! strcmp(atomNeighbors[4*iAtom+2].symbol, "P1")) numSe++;
                if (! strcmp(atomNeighbors[4*iAtom+3].symbol, "P1")) numSe++;
            }
            if (! strcmp(atomNeighbors[4*iAtom+2].symbol, "Se")) numSe++;
            if (! strcmp(atomNeighbors[4*iAtom+3].symbol, "Se")) numSe++;

            // printf("%ld %s %i\n", iAtom, atoms[iAtom].symbol, numSe);

            // Se center
            if (! strcmp(atoms[iAtom].symbol, "Se")) {
                // printf("here Se!\n");
                if (!crystalStructure) refTetrahedronVol[iAtom] = 9.36357042145021;
                if (crystalStructure) refTetrahedronVol[iAtom] = 9.264686495055233;
            }
            // S center
            else if (! strcmp(atoms[iAtom].symbol, "S")) {
                // printf("here S!\n");
                if (!crystalStructure) refTetrahedronVol[iAtom] = 8.303025459832083;
                if (crystalStructure) refTetrahedronVol[iAtom] = 8.205905479438373;
            }
            // Cd center
            else if (! strcmp(atoms[iAtom].symbol, "Cd")) {
                // printf("here Cd!\n");
                // 4 Se neighbors
                if (numSe == 4) {
                    if (!crystalStructure) refTetrahedronVol[iAtom] = 9.36357042145021;
                    if (crystalStructure) refTetrahedronVol[iAtom] = 9.264686495055233;
                }
                // 3 Se neighbors, 1 S neighbor
                else if (numSe == 3) {
                    if (!crystalStructure) refTetrahedronVol[iAtom] = 9.087742385978258;
                    if (crystalStructure) refTetrahedronVol[iAtom] = 8.989214552610866;
                }
                // 2 Se neighbors, 2 S neighbors
                else if (numSe == 2) {
                    if (!crystalStructure) refTetrahedronVol[iAtom] = 8.819136771653406;
                    if (crystalStructure) refTetrahedronVol[iAtom] = 8.721023282605593;
                }
                // 1 Se neighbor,  3 S neighbors
                else if (numSe == 1) {
                    if (!crystalStructure) refTetrahedronVol[iAtom] = 8.557611741822427;
                    if (crystalStructure) refTetrahedronVol[iAtom] = 8.459968364920925;
                }
                // 4 S neighbors
                else {
                    if (!crystalStructure) refTetrahedronVol[iAtom] = 8.303025459832083;
                    if (crystalStructure) refTetrahedronVol[iAtom] = 8.205905479438373;
                }
            }
            // otherwise
            else {
                printf("\n\nTetrahedron reference volumes not implemented for atom type %s -- the program is exiting!!!\n\n",
                        atoms[iAtom].symbol);
                fflush(stdout);
                exit(EXIT_FAILURE);
            }
        }
    }

    return;
}

/*****************************************************************************/
/* calculate actual tetrahedron volume for each atom */
void calculateTetrahedronVol(long nAtoms, atom *atoms, atom *atomNeighbors, double *tetrahedronVol) {

    FILE *pf;
    long iAtom;
    vector v1, v2, v3;

    for (iAtom = 0; iAtom < nAtoms; iAtom++) {

        if (! isAPassivationSymbol(atoms[iAtom].symbol)) {

            // rescale positions of passivation ligands
            if (! strcmp(atomNeighbors[4*iAtom+2].symbol, "P1")) {
                atomNeighbors[4*iAtom+2].pos = retScaledVector(atomNeighbors[4*iAtom+2].pos, 1.0/0.55);
            }
            if (! strcmp(atomNeighbors[4*iAtom+3].symbol, "P1")) {
                atomNeighbors[4*iAtom+3].pos = retScaledVector(atomNeighbors[4*iAtom+3].pos, 1.0/0.55);
            }
            if (! strcmp(atomNeighbors[4*iAtom+2].symbol, "P2") && ! strcmp(atomNeighbors[4*iAtom+3].symbol, "P2")) {
                atomNeighbors[4*iAtom+2].pos = retScaledVector(atomNeighbors[4*iAtom+2].pos, 1./0.25);
                atomNeighbors[4*iAtom+3].pos = retScaledVector(atomNeighbors[4*iAtom+3].pos, 1./0.25);
            }
            else if (! strcmp(atomNeighbors[4*iAtom+3].symbol, "P2")) {
                atomNeighbors[4*iAtom+3].pos = retScaledVector(atomNeighbors[4*iAtom+3].pos, 1.0/0.30);
            }

            // calculate tetrahedron volume
            v1 = retSubtractedVectors(atomNeighbors[4*iAtom].pos, atomNeighbors[4*iAtom+3].pos);
            v2 = retSubtractedVectors(atomNeighbors[4*iAtom+1].pos, atomNeighbors[4*iAtom+3].pos);
            v3 = retSubtractedVectors(atomNeighbors[4*iAtom+2].pos, atomNeighbors[4*iAtom+3].pos);

            tetrahedronVol[iAtom] = retDotProduct(v1, retCrossProduct(v2, v3))/6.;
        }
    }

    return;
}

/*****************************************************************************/
/* calculate derivative of neighbor's tetrahedron volume with respect to atom position */
void calculateTetrahedronVolDeriv(long nAtoms, atom *atoms, atom *atomNeighbors, double *tetrahedronVol, vector *tetrahedronVolDerivatives) {

    long iAtom, iNeighbor, neighborIndex;
    int nsign;
    vector v1, v2, v3;

    for (iAtom = 0; iAtom < nAtoms; iAtom++) {

        if (isAPassivationSymbol(atoms[iAtom].symbol)) {
            for (iNeighbor = 0; iNeighbor < 4; iNeighbor++) {
                tetrahedronVolDerivatives[4*iAtom+iNeighbor] = retZeroVector();
            }
        }
        else {
            // iterate through iAtom's neighbors
            for (iNeighbor = 0; iNeighbor < 4; iNeighbor++) {

                neighborIndex = atomNeighbors[4*iAtom+iNeighbor].index;
                // printf("iAtom: %ld, iNeighbor: %ld, neighborIndex: %ld\n", iAtom, iNeighbor, neighborIndex);

                if (isAPassivationSymbol(atoms[neighborIndex].symbol)) {
                    tetrahedronVolDerivatives[4*iAtom+iNeighbor] = retZeroVector();
                }
                else {

                    if (tetrahedronVol[neighborIndex] > 0.) {
                        nsign = 1.;
                    }
                    else {
                        nsign = -1.;
                    }
                    // printf("%ld %ld %i\n", iAtom, neighborIndex, nsign);


                    // printf("iAtom: %ld, iNeighbor: %ld, neighborIndex: %ld\n", iAtom, iNeighbor, neighborIndex);
                    // printf("%ld %ld %ld %ld\n", atomNeighbors[4*neighborIndex].index, atomNeighbors[4*neighborIndex+1].index, atomNeighbors[4*neighborIndex+2].index, atomNeighbors[4*neighborIndex+3].index);

                    // if iAtom is neighborIndex's a atom, take derivatives accordingly
                    if (iAtom == atomNeighbors[4*neighborIndex].index) {
                        v1 = retSubtractedVectors(atomNeighbors[4*neighborIndex+1].pos, atomNeighbors[4*neighborIndex+3].pos);
                        v2 = retSubtractedVectors(atomNeighbors[4*neighborIndex+2].pos, atomNeighbors[4*neighborIndex+3].pos);

                        tetrahedronVolDerivatives[4*iAtom+iNeighbor] = retScaledVector(retCrossProduct(v1, v2), nsign/6.);
                        // v3 = retCrossProduct(v1, v2);
                        // printf("%i\n", nsign);
                        // printf("in here a! %.8f %.8f %.8f\n", v1.x, v1.y, v1.z);
                        // printf("in here a! %.8f %.8f %.8f\n", v2.x, v2.y, v2.z);
                        // printf("in here a! %.8f %.8f %.8f\n", v3.x, v3.y, v3.z);
                        // printf("in here a! %.8f %.8f %.8f\n", tetrahedronVolDerivatives[4*iAtom+iNeighbor].x, tetrahedronVolDerivatives[4*iAtom+iNeighbor].y, tetrahedronVolDerivatives[4*iAtom+iNeighbor].z);
                    }
                    // if iAtom is neighborIndex's b atom
                    else if (iAtom == atomNeighbors[4*neighborIndex+1].index) {
                        v1 = retSubtractedVectors(atomNeighbors[4*neighborIndex].pos, atomNeighbors[4*neighborIndex+3].pos);
                        v2 = retSubtractedVectors(atomNeighbors[4*neighborIndex+2].pos, atomNeighbors[4*neighborIndex+3].pos);

                        tetrahedronVolDerivatives[4*iAtom+iNeighbor].x = nsign*retDotProduct(v1, retCrossProduct(retIVector(), v2))/6.;
                        tetrahedronVolDerivatives[4*iAtom+iNeighbor].y = nsign*retDotProduct(v1, retCrossProduct(retJVector(), v2))/6.;
                        tetrahedronVolDerivatives[4*iAtom+iNeighbor].z = nsign*retDotProduct(v1, retCrossProduct(retKVector(), v2))/6.;
                        // printf("in here b! %.8f %.8f %.8f\n", tetrahedronVolDerivatives[4*iAtom+iNeighbor].x, tetrahedronVolDerivatives[4*iAtom+iNeighbor].x, tetrahedronVolDerivatives[4*iAtom+iNeighbor].x);
                    }
                    // if iAtom is neighborIndex's c atom
                    else if (iAtom == atomNeighbors[4*neighborIndex+2].index) {
                        v1 = retSubtractedVectors(atomNeighbors[4*neighborIndex].pos, atomNeighbors[4*neighborIndex+3].pos);
                        v2 = retSubtractedVectors(atomNeighbors[4*neighborIndex+1].pos, atomNeighbors[4*neighborIndex+3].pos);

                        tetrahedronVolDerivatives[4*iAtom+iNeighbor].x = nsign*retDotProduct(v1, retCrossProduct(v2, retIVector()))/6.;
                        tetrahedronVolDerivatives[4*iAtom+iNeighbor].y = nsign*retDotProduct(v1, retCrossProduct(v2, retJVector()))/6.;
                        tetrahedronVolDerivatives[4*iAtom+iNeighbor].z = nsign*retDotProduct(v1, retCrossProduct(v2, retKVector()))/6.;
                        // printf("in here c! %.8f %.8f %.8f\n", tetrahedronVolDerivatives[4*iAtom+iNeighbor].x, tetrahedronVolDerivatives[4*iAtom+iNeighbor].x, tetrahedronVolDerivatives[4*iAtom+iNeighbor].x);
                    }
                    // if iAtom is neighborIndex's d atom
                    else if (iAtom == atomNeighbors[4*neighborIndex+3].index) {
                        v1 = retSubtractedVectors(atomNeighbors[4*neighborIndex].pos, atomNeighbors[4*neighborIndex+3].pos);
                        v2 = retSubtractedVectors(atomNeighbors[4*neighborIndex+1].pos, atomNeighbors[4*neighborIndex+3].pos);
                        v3 = retSubtractedVectors(atomNeighbors[4*neighborIndex+2].pos, atomNeighbors[4*neighborIndex+3].pos);

                        tetrahedronVolDerivatives[4*iAtom+iNeighbor] = retScaledVector(retCrossProduct(v2, v3), -nsign/6.);

                        tetrahedronVolDerivatives[4*iAtom+iNeighbor].x -= nsign*retDotProduct(v1, retCrossProduct(retIVector(), v3))/6.;
                        tetrahedronVolDerivatives[4*iAtom+iNeighbor].y -= nsign*retDotProduct(v1, retCrossProduct(retJVector(), v3))/6.;
                        tetrahedronVolDerivatives[4*iAtom+iNeighbor].z -= nsign*retDotProduct(v1, retCrossProduct(retKVector(), v3))/6.;

                        tetrahedronVolDerivatives[4*iAtom+iNeighbor].x -= nsign*retDotProduct(v1, retCrossProduct(v2, retIVector()))/6.;
                        tetrahedronVolDerivatives[4*iAtom+iNeighbor].y -= nsign*retDotProduct(v1, retCrossProduct(v2, retJVector()))/6.;
                        tetrahedronVolDerivatives[4*iAtom+iNeighbor].z -= nsign*retDotProduct(v1, retCrossProduct(v2, retKVector()))/6.;
                    }
                    else {
                        printf("\n\nError calculating tetrahedron volume derivatives. Atom %ld is not a neighbor of atom %ld -- exiting!!!\n\n", iAtom, iNeighbor);
                        fflush(stdout);
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
    }

    return;
}

/*****************************************************************************/
/* compute strain scaling factor for each atom
 * strainScale = 1 + a4 * (Omega/Omega0 - 1),
 * where Omega is volume of tetrahedron formed by atom's nearest neighbors
 * and Omega0 is volume of reference tetrahedron */
void calculateStrainScale(long nAtoms, atom *atoms, double *tetrahedronVolRef, double *tetrahedronVol, double *a4Params, double *strainScale) {

    FILE *pf;
    long iAtom;

    for (iAtom = 0; iAtom < nAtoms; iAtom++) {

        if (isAPassivationSymbol(atoms[iAtom].symbol)) {
            strainScale[iAtom] = 1.;
        }
        else {
            strainScale[iAtom] = (1. + a4Params[atoms[iAtom].type]*(fabs(tetrahedronVol[iAtom])/tetrahedronVolRef[iAtom] - 1.));
        }
    }

    pf = fopen("strain.dat", "w");
    for (iAtom = 0; iAtom < nAtoms; iAtom++) {
        fprintf(pf, "%ld %s %.8f %.8f %.8f %.8f\n", iAtom, atoms[iAtom].symbol, a4Params[atoms[iAtom].type],
                fabs(tetrahedronVol[iAtom]), tetrahedronVolRef[iAtom], strainScale[iAtom]);
    }
    fclose(pf);

    return;
}

/*****************************************************************************/
/* compute strain scaling factor for each atom
 * strainScale = 1 + a4 * (Omega/Omega0 - 1),
 * where Omega is volume of tetrahedron formed by atom's nearest neighbors
 * and Omega0 is volume of reference tetrahedron */
void calculateStrainScaleDeriv(long nAtoms, atom *atoms, atom *atomNeighbors, 
        double *tetrahedronVolRef, vector *tetrahedronVolDeriv, double *a4Params, vector *strainScaleDeriv) {

    FILE *pf;
    long iAtom, iNeighbor, neighborIndex;

    for (iAtom = 0; iAtom < nAtoms; iAtom++) {

        if (isAPassivationSymbol(atoms[iAtom].symbol)) {
            for (iNeighbor = 0; iNeighbor < 4; iNeighbor++) {
                strainScaleDeriv[4*iAtom+iNeighbor] = retZeroVector();
            }
        }
        else {

            for (iNeighbor = 0; iNeighbor < 4; iNeighbor++) {
                
                neighborIndex = atomNeighbors[4*iAtom+iNeighbor].index;

                if (isAPassivationSymbol(atoms[neighborIndex].symbol)) {
                    strainScaleDeriv[4*iAtom+iNeighbor] = retZeroVector();
                }
                else {
                    strainScaleDeriv[4*iAtom+iNeighbor] = retScaledVector(tetrahedronVolDeriv[4*iAtom+iNeighbor], 
                            a4Params[atoms[neighborIndex].type]/tetrahedronVolRef[neighborIndex]);
                }
            }
        }
    }

    pf = fopen("tetrahedronVolDeriv.dat", "w");
    for (iAtom = 0; iAtom < nAtoms; iAtom++) {
        for (iNeighbor = 0; iNeighbor < 4; iNeighbor++) {
            fprintf(pf, "%ld %s %ld %ld %s %.8f %.8f %.8f\n", iAtom, atoms[iAtom].symbol, iNeighbor, 
                    atomNeighbors[4*iAtom+iNeighbor].index, atomNeighbors[4*iAtom+iNeighbor].symbol,
                    tetrahedronVolDeriv[4*iAtom+iNeighbor].x, tetrahedronVolDeriv[4*iAtom+iNeighbor].y, tetrahedronVolDeriv[4*iAtom+iNeighbor].z);
        }
    }
    fclose(pf);
    
    pf = fopen("strainDeriv.dat", "w");
    for (iAtom = 0; iAtom < nAtoms; iAtom++) {
        for (iNeighbor = 0; iNeighbor < 4; iNeighbor++) {
            fprintf(pf, "%ld %s %ld %ld %s %.8f %.8f %.8f %.8f\n", iAtom, atoms[iAtom].symbol, iNeighbor, 
                    atomNeighbors[4*iAtom+iNeighbor].index, atomNeighbors[4*iAtom+iNeighbor].symbol, a4Params[atoms[atomNeighbors[4*iAtom+iNeighbor].index].type],
                    strainScaleDeriv[4*iAtom+iNeighbor].x, strainScaleDeriv[4*iAtom+iNeighbor].y, strainScaleDeriv[4*iAtom+iNeighbor].z);
        }
    }

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
  else if ((atyp[0] == 'P') && (atyp[1] == '1') && (atyp[2] == '\0'))  return(8);
  else if ((atyp[0] == 'P') && (atyp[1] == '2') && (atyp[2] == '\0'))  return(9);
  else if ((atyp[0] == 'P') && (atyp[1] == '3') && (atyp[2] == '\0'))  return(10);
  else if ((atyp[0] == 'P') && (atyp[1] == '4') && (atyp[2] == '\0'))  return(11);
  else if ((atyp[0] == 'T') && (atyp[1] == 'e') && (atyp[2] == '\0')) return(12);
  else if ((atyp[0] == 'C') && (atyp[1] == 'd') && (atyp[2] == 'z')) return(13);
  else if ((atyp[0] == 'S') && (atyp[1] == 'e') && (atyp[2] == 'z')) return(14);
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
  return;
}

/*****************************************************************************/
