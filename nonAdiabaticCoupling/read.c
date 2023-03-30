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
void readAtomicPseudopotentials(grid1d *atomicPPGrid, gridPoint1d *atomicPPGP, atom *atoms, lParams lPar) {
  FILE *pf;
  long iAtom, iAtomType, iGrid;
  long nGridPoints = atomicPPGrid->nGridPoints;
  char fileName[100];

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
