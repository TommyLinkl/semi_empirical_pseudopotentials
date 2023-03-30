/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "es.h"

/*****************************************************************************/
// 

long readInputParFile(grid *rSpaceGrid, dParams *dPar, lParams *lPar) {
  FILE *pf; 
  long i, j, ieof, lenEvalParFile;
  double intExcEnergy, qpEnergy, sigma, *eval, *sige, a;
  char field[100], tmp[100];
  
  // Set defaults
  dPar->sigmaCutoff = 0.01;   // sigma=sqrt(<E^2>-<E>^2) -> values lower than sigmaCutoff are eigenstates of H
  dPar->temp = 298.0;
  dPar->kbT = dPar->temp*KB;
  dPar->boltzEnergyRange = 2.0*dPar->kbT/AUTOEV;
  lPar->intExcitons = 1;
  lPar->nThreads = 1;

  // Read input.par if it exists - exit program otherwise 
  if ( access("input.par", F_OK) != -1 ) {
    pf = fopen("input.par", "r");
    i = 0;
    while (fscanf(pf, "%s", field) != EOF && i < 4) {
      if (! strcmp(field, "gridSize")) {
        fscanf(pf, "%s %ld %ld %ld", tmp, &(rSpaceGrid->nGridPointsX), 
                                          &(rSpaceGrid->nGridPointsY), 
                                          &(rSpaceGrid->nGridPointsZ));
      }
      else if (! strcmp(field, "nThreads"))  fscanf(pf, "%s %ld", tmp, &(lPar->nThreads));
      else if (! strcmp(field, "intExcitons")) fscanf(pf, "%s %s", tmp, &(lPar->intExcitons));
      else if (! strcmp(field, "sigmaCutoff"))  fscanf(pf, "%s %lg", tmp, &(dPar->sigmaCutoff));
      else {
        printf("Invalid input field and/or format - equal sign required after each field\n");
        printf("Only allowed fields are (case-sensitive):\n\n");
        printf("gridSize = nx ny nz (required, all integers)\n");
        printf("intExcitons = 1 (default = 0 which means no interacting exciton properties will be calculated)\n");
        printf("nThreads = 1 (optional, serial default, number of openmp threads)\n");
        printf("sigmaCutoff = 0.01 (optional, 0.01 default, criteria to determine if eigenstate or not)\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
      }
      i++;
    }
    fclose(pf);
  }
  else {
    printf("\n\nNo input.par file detected in current working directory - the program is exiting!!!\n\n");
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
      if (sigma < dPar->sigmaCutoff && qpEnergy < -0.175) {
        dPar->homoEnergy = qpEnergy;
        lPar->iHomo = i;
        lPar->nHoles++;
      }
      else if (sigma < dPar->sigmaCutoff && qpEnergy > -0.175 && ieof != EOF) {
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
  
  // Set useful parameters based on eval.par file
  lPar->nHolesPlusElecs = lPar->nHoles+lPar->nElecs;
  lPar->nNonintExcitons = lPar->nHoles*lPar->nElecs; 
  dPar->fundamentalGap = dPar->lumoEnergy - dPar->homoEnergy;
  
  // Determine number of interacting excitonic states from the exciton.par file
  lPar->nIntExcitons = 0;
  if (lPar->intExcitons) {
    if ( access("exciton.par", F_OK) != -1 ) {
      pf = fopen("exciton.par", "r");
      for (i = ieof = 0; ieof != EOF; i++) {
        ieof = fscanf(pf, "%ld %lg %lg %lg %lg", &j, &intExcEnergy, &a, &a, &a);
        if (! i) {
          dPar->opticalGap = intExcEnergy;
          lPar->nIntExcitons++;
        }
        else if (intExcEnergy <= (dPar->opticalGap + dPar->boltzEnergyRange) && lPar->nIntExcitons < 100) {
          lPar->nIntExcitons++;
        }
        else {
          break; // all other excitonic states are of too high energy
        }
      } 
      fclose(pf);
    }
    else {
      printf("\n\nNo exciton.par file detected in current working directory - the program is exiting!!!\n\n");
      exit(EXIT_FAILURE);
    }
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
    if (tmpQPSigma[i] < sigmaCutoff) {
      holeQP[nStoredHoles].index = nStoredHoles;
      strcpy(holeQP[nStoredHoles].type, "h");
      holeQP[nStoredHoles].energy = tmpQPEnergy[i];
      holeQP[nStoredHoles].sigma = tmpQPSigma[i];
      nStoredHoles++;
      // evalbe[lPar.nHoles-nStoredHoles] = tmpQPEnergy[i];
      // sigebe[lPar.nHoles-nStoredHoles] = tmpQPSigma[i];
      eigenstateList[lPar.nHoles-nStoredHoles] = i;
    }
    if (nStoredHoles == lPar.nHoles) {
      break;
    }
  }
  nStoredElecs = 0;
  for (i = lPar.iLumo; i < nStatesInEval; i++) {
    if (tmpQPSigma[i] < sigmaCutoff) {
      elecQP[nStoredElecs].index = nStoredElecs;
      strcpy(elecQP[nStoredElecs].type, "e");
      elecQP[nStoredElecs].energy = tmpQPEnergy[i];
      elecQP[nStoredElecs].sigma = tmpQPSigma[i];
      // evalbe[lPar.nHoles+nStoredElecs] = tmpQPEnergy[i];
      // sigebe[lPar.nHoles+nStoredElecs] = tmpQPSigma[i];
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
  */

  for (j = 0; j < ntype*n; j++) pot[j] = vr[j] = 0;
  for (j = 0; j < ntype; j++) nPot[j] = 0;
  
  for (j = 0; j < ntype; j++){
    assign_atom_type(atype,j);
    if ((j==5) || (j==7)) sprintf (str,"pot%c.par",atype[0]);
    else if (j <= 12) sprintf (str,"pot%c%c.par",atype[0],atype[1]);
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
    if ((j==5) || (j==7)) sprintf (str,"pot%c.dat",atype[0]);
    else if (j <= 12) sprintf (str,"pot%c%c.dat",atype[0],atype[1]);
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
