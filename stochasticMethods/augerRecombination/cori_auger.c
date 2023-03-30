/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "ar.h"

/*****************************************************************************/

void calcDoublyStochasticIntAR(double *Cbs, double *Ebs, double *psibe, double *evalbe, double *psiai, double *evalai, zomplex *potq,
                                lng_st ist, par_st par, fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwpsi)
{
  FILE *pf;
  long a, i, b, j, nSO, nSI, nSA;
  long iDeltaE, iBiexc, iExc1, iExc2, iBiexcDeltaE, iZeta, idum, counter, iZetanThreads, nThetaZeta;
  long randomHoleIndex, randomElecIndex, randomHotHoleIndex, randomHotElecIndex;
  long *globalThetaIList, *globalThetaAList, *globalThetaIElecList, *globalThetaAHoleList;
  long *elecHotHoleList, *holeHotElecList, *nHotElecForHole, *nHotHoleForElec;
  long *thetaIList, *thetaIElecList, *thetaAList, *thetaAHoleList;
  double sum, tmp, *nHotHnSIRatio, *nHotEnSARatio;  
  double *eRate, *hRate, *arRate, *biexcitonEnergies;
  double deltaE[ist.nConsWindows];
  char fileName[50];
  zomplex *RckZeta, *RabZeta, *RijZeta, *TZeta;
  zomplex *ChiaiZetaHole, *ChiaiZetaElec;
  zomplex *AaiHole, *AaiElec, *AaiHole1, *AaiElec1, *AaiHole2, *AaiElec2;
  zomplex *thetaZeta;

  // Useful integers
  idum = ist.seed; 
  nSI = ist.nStochHotHoles;
  nSA = ist.nStochHotElecs;
  nSO = ist.nStochOrbitals;  
  const long THL = ist.nNonIntExcitons;
  const long TH2L = THL * ist.nHoles;
  const long THL2 = THL * ist.nElecs;

  // Write beginning of function
  writeSeparation(stdout);
  writeCurrentTime(stdout);
  fprintf(stdout, "Beginnning the calculation of the interacting AR lifetimes using the doubly stochastic method\n\n");
  fprintf(stdout, "The number of energy conservation windows = %d\n", ist.nConsWindows);
  fprintf(stdout, "The largest energy conservation window = %.6f %.6f\n", par.maxDeltaE, par.maxDeltaE*AUTOEV);
  fprintf(stdout, "The number of deterministic hot hole states = %d\n", ist.itot);
  fprintf(stdout, "The number of deterministic hot elec states = %d\n", ist.atot);
  fprintf(stdout, "The number of stochastic hot hole states sampled = %ld\n", nSI);
  fprintf(stdout, "The number of stochastic hot elec states sampled = %ld\n", nSA);
  fprintf(stdout, "The number of stochastic orbitals used to approximate the Coulomb operator = %ld\n", nSO);
  fflush(stdout);

  // Dynamically allocate memory
  if ((eRate = (double *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(double))) == NULL) nerror("eRate");
  if ((hRate = (double *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(double))) == NULL) nerror("hRate");
  if ((arRate = (double *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(double))) == NULL) nerror("arRate");
  if ((biexcitonEnergies = (double *) calloc(ist.nBiexcitons, sizeof(double))) == NULL) nerror("biexcitonEnergies");
  if ((globalThetaIList = (long *) calloc(nSI, sizeof(long))) == NULL) nerror("globalThetaIList");
  if ((globalThetaIElecList = (long *) calloc(nSI, sizeof(long))) == NULL) nerror("globalThetaIElecList");
  if ((globalThetaAList = (long *) calloc(nSA, sizeof(long))) == NULL) nerror("globalThetaAList");
  if ((globalThetaAHoleList = (long *) calloc(nSA, sizeof(long))) == NULL) nerror("globalThetaAHoleList");
  if ((elecHotHoleList = (long *) calloc(ist.nBiexcitons*ist.nConsWindows*ist.nElecs*ist.itot, sizeof(long))) == NULL) nerror("elecHotHoleList");
  if ((nHotHoleForElec = (long *) calloc(ist.nBiexcitons*ist.nConsWindows*ist.nElecs, sizeof(long))) == NULL) nerror("nHotHoleForElec");
  if ((holeHotElecList = (long *) calloc(ist.nBiexcitons*ist.nConsWindows*ist.nHoles*ist.atot, sizeof(long))) == NULL) nerror("holeHotElecList");
  if ((nHotElecForHole = (long *) calloc(ist.nBiexcitons*ist.nConsWindows*ist.nHoles, sizeof(long))) == NULL) nerror("nHotElecForHole");
  if ((thetaIList = (long *) calloc(nSI*ist.nBiexcitons*ist.nConsWindows, sizeof(long))) == NULL) nerror("thetaIList");
  if ((thetaIElecList = (long *) calloc(nSI*ist.nBiexcitons*ist.nConsWindows, sizeof(long))) == NULL) nerror("thetaIElecList");
  if ((thetaAList = (long *) calloc(nSA*ist.nBiexcitons*ist.nConsWindows, sizeof(long))) == NULL) nerror("thetaAList");
  if ((thetaAHoleList = (long *) calloc(nSA*ist.nBiexcitons*ist.nConsWindows, sizeof(long))) == NULL) nerror("thetaAHoleList");
  if ((nHotHnSIRatio = (double *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(double))) == NULL) nerror("nHotHnSIRatio");
  if ((nHotEnSARatio = (double *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(double))) == NULL) nerror("nHotEnSARatio");

  // Create multiple energy conservation windows so one does not have to repeat the entire calculation
  fillEnergyConservationWindows(deltaE, par.maxDeltaE, ist.nConsWindows);

  // Fill in biexcitonEnergies
  fillIntBiexcEnergies(biexcitonEnergies, Ebs, ist, par);

  // Create lists and counts of the allowed final states for each of the possible initial states and energyConsWindow
  storeAllowedFinalStatesIntAR(elecHotHoleList, nHotHoleForElec, holeHotElecList, nHotElecForHole, biexcitonEnergies, 
                                Ebs, evalbe, evalai, deltaE, ist, par);

  // Sampling the final states part of the function
  // Generate the stochastic orbitals by:
  // 1) randomly select a band-edge state 
  // 2) randomly select a hot carrier from the list of hot carriers such that the elec-hole pair
  //    satisfies energy conservation with at least one initial biexcitonic state
  if (ist.readAaiMatrices) {
    counter = 0;
    //readGlobalThetaLists(globalThetaIElecList, globalThetaIList, nSI, globalThetaAHoleList, globalThetaAList, nSA);
  }
  else {
    counter = 0;  
    for (i = 0; i < nSI; ) {
      randomElecIndex = (long)((double)(ist.nElecs)*ran_nrc(&idum)) + ist.nHoles; 
      randomHotHoleIndex = (long)((double)(ist.itot)*ran_nrc(&idum));
      counter++;
      if (ist.itot == 0 || counter > nSI*100) break; // test to leave loop if no final states 
      // Make sure this final state satisfies energy conservation for at least 1 initial state
      for (iBiexc = 0; iBiexc < ist.nBiexcitons; iBiexc++) {
        if (fabs(evalbe[randomElecIndex] - evalai[randomHotHoleIndex] - biexcitonEnergies[iBiexc]) <= par.maxDeltaE) {
          globalThetaIElecList[i] = randomElecIndex; // nHoles to nHoles+nElecs
          globalThetaIList[i] = randomHotHoleIndex;  // 0 to itot
          i++;
          break;
        }
      }
    }
    counter = 0;
    for (a = 0; a < nSA; ) {
      randomHoleIndex = (long)((double)(ist.nHoles)*ran_nrc(&idum));
      randomHotElecIndex = (long)((double)(ist.atot)*ran_nrc(&idum)) + ist.itot;
      counter++;
      if (ist.atot == 0 || counter > nSA*100) break; // test to leave loop if no final states
      // Make sure this final state satisfies energy conservation for at least 1 initial state
      for (iBiexc = 0; iBiexc < ist.nBiexcitons; iBiexc++) {
        if (fabs(evalai[randomHotElecIndex] - evalbe[randomHoleIndex] - biexcitonEnergies[iBiexc]) <= par.maxDeltaE) {
          globalThetaAHoleList[a] = randomHoleIndex; // 0 to nHoles
          globalThetaAList[a] = randomHotElecIndex;  // itot to itot+atot
          a++;
          break;
        }
      }
    }
  }

  // Write globalTheta lists so s3I calcutions can be restarted
  // Hole channel
  pf = fopen("holeChannelThetaLists.dat", "w");
  fprintf(pf, "%ld", nSI);
  for (i = 0; i < nSI; i++) {
    fprintf(pf, "%ld %ld %ld", i, globalThetaIElecList[i], globalThetaIList[i]);
  }
  fclose(pf);
  // Elec channel
  pf = fopen("elecChannelThetaLists.dat", "w");
  fprintf(pf, "%ld", nSA);
  for (a = 0; a < nSA; a++) {
    fprintf(pf, "%ld %ld %ld", a, globalThetaAHoleList[a], globalThetaAList[a]);
  }
  fclose(pf);

  // Fill thetaIList and thetaAList based on the globalTheta Lists (energy conservation check here)
  for (iDeltaE = 0; iDeltaE < ist.nConsWindows; iDeltaE++) {
    for (iBiexc = 0; iBiexc < ist.nBiexcitons; iBiexc++) {
      iBiexcDeltaE = iBiexc + iDeltaE*ist.nBiexcitons;
      // Hole channel
      for (i = 0; i < nSI; i++) {
        if (fabs(evalbe[globalThetaIElecList[i]] - evalai[globalThetaIList[i]] - biexcitonEnergies[iBiexc]) <= deltaE[iDeltaE]) {
          thetaIElecList[i + nSI*iBiexcDeltaE] = (globalThetaIElecList[i] - ist.nHoles); // only to be used in w2h indexing, 0 to nElecs
          thetaIList[i + nSI*iBiexcDeltaE] = globalThetaIList[i]; // 0 to itot
        }
        else { // used to signify it does not satisfy energy conservation
          thetaIElecList[i + nSI*iBiexcDeltaE] = -1; 
          thetaIList[i + nSI*iBiexcDeltaE] = -1;
        }
      }
      // Elec channel
      for (a = 0; a < nSA; a++) {
        if (fabs(evalai[globalThetaAList[a]] - evalbe[globalThetaAHoleList[a]] - biexcitonEnergies[iBiexc]) <= deltaE[iDeltaE]) {
          thetaAHoleList[a + nSA*iBiexcDeltaE] = globalThetaAHoleList[a]; // only to be used in w2e indexing, 0 to nHoles
          thetaAList[a + nSA*iBiexcDeltaE] = globalThetaAList[a];  // itot to itot+atot
        }
        else { // used to signify it does not satisfy energy conservation
          thetaAHoleList[a + nSA*iBiexcDeltaE] = -1; 
          thetaAList[a + nSA*iBiexcDeltaE] = -1;
        }
      }
    }
  }

  // Determine the required scalings for each pair of initial state and energy conservation window
  double sumSto, sumDet;
  for (iDeltaE = 0; iDeltaE < ist.nConsWindows; iDeltaE++) {
    iBiexc = 0;
    for (iExc1 = 0; iExc1 < ist.nIntExcitons; iExc1++) {
      for (iExc2 = 0; iExc2 < ist.nIntExcitons; iExc2++) { 
        if ((Ebs[iExc1] + Ebs[iExc2]) < par.maxInitE+EPS) {
          iBiexcDeltaE = iBiexc + iDeltaE*ist.nBiexcitons;
          // Hole channel
          for (sumSto = 0.0, i = 0; i < nSI; i++) {
            a = thetaIElecList[i + nSI*iBiexcDeltaE]; // a goes from 0 to ist.nElecs
            if (a != -1) { // final states that do not satisfy energy conservation have a = i = -1
              for (tmp = 0.0, j = 0; j < ist.nHoles; j++) { 
                tmp += Cbs[iExc1*ist.nNonIntExcitons + a*ist.nHoles + j];
              }
              sumSto += tmp*tmp;
            }
          }
          for (sumDet = 0.0, a = 0; a < ist.nElecs; a++) { // sum over a,j for c_a,j, delta_ab
            for (tmp = 0.0, j = 0; j < ist.nHoles; j++) {
              tmp += Cbs[iExc1*ist.nNonIntExcitons + a*ist.nHoles + j];
            }
            sumDet += tmp*tmp*(double)(nHotHoleForElec[a + iBiexcDeltaE*ist.nElecs]); 
          }
          nHotHnSIRatio[iBiexcDeltaE] = sumDet/sumSto;
          // Elec channel
          for (sumSto = 0.0, a = 0; a < nSA; a++) {
            i = thetaAHoleList[a + nSA*iBiexcDeltaE]; // a goes from 0 to ist.nHoles
            if (i != -1) { // final states that do not satisfy energy conservation have i = a = -1
              for (tmp = 0.0, b = 0; b < ist.nElecs; b++) {
                tmp += Cbs[iExc1*ist.nNonIntExcitons + b*ist.nHoles + i];
              }
              sumSto += tmp*tmp;
            }
          }
          for (sumDet = 0.0, i = 0; i < ist.nHoles; i++) { // sum over b,i for c_b,i, delta_ij
            for (tmp = 0.0, b = 0; b < ist.nElecs; b++) {
              tmp += Cbs[iExc1*ist.nNonIntExcitons + b*ist.nHoles + i];
            }
            sumDet += tmp*tmp*(double)(nHotElecForHole[i + iBiexcDeltaE*ist.nHoles]);
          }
          nHotEnSARatio[iBiexcDeltaE] = sumDet/sumSto;          
          iBiexc++;
        }
      }
    }
  }

  // Remove states from psiai/evalai that are not being sampled 
  remStatesNotInThetaLists(psiai, evalai, thetaIList, thetaAList, &ist);

  // Realloc memory for the smaller psiai, evalai and allocate memory for the stochastic Coulomb part of this function
  psiai = realloc(psiai, ist.nGridPoints*ist.nHotNonIntStates*sizeof(psiai[0]));
  evalai = realloc(evalai, ist.nHotNonIntStates*sizeof(evalai[0]));

  // Stochastic Coulomb part of the function
  if ((thetaZeta = (zomplex *) calloc(ist.nGridPoints*ist.nThreads, sizeof(zomplex))) == NULL) nerror("thetaZeta");
  if ((RckZeta = (zomplex *) calloc(ist.nElecs*ist.nHoles*ist.nThreads, sizeof(zomplex))) == NULL) nerror("RckZeta");  
  if ((TZeta = (zomplex *) calloc(ist.nIntExcitons*ist.nThreads, sizeof(zomplex))) == NULL) nerror("TZeta");  
  if ((RijZeta = (zomplex *) calloc(ist.itot*ist.nHoles*ist.nThreads, sizeof(zomplex))) == NULL) nerror("RijZeta");  
  if ((ChiaiZetaHole = (zomplex *) calloc(ist.itot*ist.nElecs*ist.nIntExcitons*ist.nThreads, sizeof(zomplex))) == NULL) nerror("ChiaiZetaHole");  
  if ((RabZeta = (zomplex *) calloc(ist.atot*ist.nElecs*ist.nThreads, sizeof(zomplex))) == NULL) nerror("RabZeta");  
  if ((ChiaiZetaElec = (zomplex *) calloc(ist.atot*ist.nHoles*ist.nIntExcitons*ist.nThreads, sizeof(zomplex))) == NULL) nerror("ChiaiZetaElec");
  if ((AaiHole = (zomplex *) calloc(ist.itot*ist.nElecs*ist.nIntExcitons*ist.nIntExcitons, sizeof(zomplex))) == NULL) nerror("AaiHole");  
  if ((AaiElec = (zomplex *) calloc(ist.atot*ist.nHoles*ist.nIntExcitons*ist.nIntExcitons, sizeof(zomplex))) == NULL) nerror("AaiElec"); 
  if ((AaiHole1 = (zomplex *) calloc(ist.itot*ist.nElecs*ist.nIntExcitons*ist.nIntExcitons, sizeof(zomplex))) == NULL) nerror("AaiHole1");  
  if ((AaiElec1 = (zomplex *) calloc(ist.atot*ist.nHoles*ist.nIntExcitons*ist.nIntExcitons, sizeof(zomplex))) == NULL) nerror("AaiElec1");  
  if ((AaiHole2 = (zomplex *) calloc(ist.itot*ist.nElecs*ist.nIntExcitons*ist.nIntExcitons, sizeof(zomplex))) == NULL) nerror("AaiHole2");  
  if ((AaiElec2 = (zomplex *) calloc(ist.atot*ist.nHoles*ist.nIntExcitons*ist.nIntExcitons, sizeof(zomplex))) == NULL) nerror("AaiElec2"); 

  // Check if previously calculated AaiMatrices should be read in
  long iInitZetanThreads, nAaiHoleMatElements, nAaiElecMatElements;
  if (ist.readAaiMatrices) {
    for (a = 0; a < 2; a++) {
      sprintf(fileName, "AaiMatrix_%ld.par", a);
      if ( access(fileName, F_OK) != -1 ) {
        pf = fopen(fileName, "r");
        fscanf(pf, "%ld %ld %ld %ld", &iInitZetanThreads, &(ist.seed), &nAaiHoleMatElements, &nAaiElecMatElements);
        printf("Reading in AaiMatrix%ld\n", a);
        printf("Number of previously calculated stochastic orbitals = %ld\n", iInitZetanThreads);
        ist.seed *= -1;
        printf("Seed set to = %ld\n", ist.seed); 
        fflush(stdout);
        for (i = 0; i < nAaiHoleMatElements; i++) {
          fscanf(pf, "%lg %lg", &(AaiHole[i].re), &(AaiHole[i].im));
        }
        for (i = 0; i < nAaiElecMatElements; i++) {
          fscanf(pf, "%lg %lg", &(AaiElec[i].re), &(AaiElec[i].im));
        }
        fclose(pf);
        if (! a && (iInitZetanThreads == (nSO/2))) {
          for (i = 0; i < nAaiHoleMatElements; i++) {
            AaiHole1[i].re = AaiHole[i].re;
            AaiHole1[i].im = AaiHole[i].im;
          }
          for (i = 0; i < nAaiElecMatElements; i++) {
            AaiElec1[i].re = AaiElec[i].re;
            AaiElec1[i].im = AaiElec[i].im;
          }
        }
      }
      else {
        printf("\nNo %s file detected in current working directory.\n", fileName); 
        if (! a) {
          iInitZetanThreads = 0;
          break; // leave loop, do not try to read AaiMatrix_1.par
        }
        printf("The number of stochastic orbitals previously calculated = %ld\n", iInitZetanThreads); 
        fflush(stdout);
      }
    }
  }
  else {
    iInitZetanThreads = 0;
  }
  if (iInitZetanThreads >= (nSO/2)) {
    for (i = 0; i < ist.itot*ist.nElecs*ist.nIntExcitons*ist.nIntExcitons; i++) {
      AaiHole1[i].re /= ((double)(nSO/2));
      AaiHole1[i].im /= ((double)(nSO/2)); 
    }
    for (i = 0; i < ist.atot*ist.nHoles*ist.nIntExcitons*ist.nIntExcitons; i++) {
      AaiElec1[i].re /= ((double)(nSO/2));
      AaiElec1[i].im /= ((double)(nSO/2)); 
    }            
  }
  if (iInitZetanThreads == (nSO/2)) {
    for (i = 0; i < ist.itot*ist.nElecs*ist.nIntExcitons*ist.nIntExcitons; i++) {
      AaiHole[i].re = 0.0;
      AaiHole[i].im = 0.0;
    }
    for (i = 0; i < ist.atot*ist.nHoles*ist.nIntExcitons*ist.nIntExcitons; i++) {
      AaiElec[i].re = 0.0;
      AaiElec[i].im = 0.0;
    }          
  }

  // Calculate stochastic orbitals that represent the Coulomb operator and their overlaps with the quasiparticle states
  fprintf(stdout, "\nBeginning the calculations of the RrsZeta, TZeta, ChiaiZeta and Aai Matrices\n");
  writeCurrentTime(stdout);
  fprintf(stdout, "The number of 3D integrals being performed = %ld\n", 
                    (nSO-iInitZetanThreads)*(ist.nElecs*ist.nHoles+ist.itot*ist.nElecs+ist.atot*ist.nHoles));
  fflush(stdout);

  for (iZetanThreads = iInitZetanThreads; iZetanThreads < nSO; iZetanThreads += ist.nThreads) {
    if ( ((iZetanThreads + ist.nThreads) > nSO) && (nSO % ist.nThreads)) {
      nThetaZeta = (nSO % ist.nThreads);
    }
    else {
      nThetaZeta = ist.nThreads;
    }
    // Generate stochastic orbitals used to approximate the Coulomb operator (the thetaZeta orbitals)
    ist.seed = calcRandomCoulombStates(thetaZeta, nThetaZeta, potq, ist, par, planfw, planbw, fftwpsi);

    // Calculate RckZeta (elec and hole) 
    calcRrsZetaMatrix(RckZeta, &(psibe[ist.nHoles*ist.nGridPoints]), ist.nElecs, psibe, ist.nHoles,
                      thetaZeta, nThetaZeta, ist.nGridPoints, par.dv);
    for (i = 0; i < ist.nElecs*ist.nHoles*nThetaZeta; i++) { RckZeta[i].im *= -1.0; }
    
    // Calculate RijZeta (hot hole and hole - hole channel) 
    calcRrsZetaMatrix(RijZeta, psiai, ist.itot, psibe, ist.nHoles, 
                      thetaZeta, nThetaZeta, ist.nGridPoints, par.dv);
    
    // Calculate RabZeta (hot elec and elec - elec channel) 
    calcRrsZetaMatrix(RabZeta, &(psiai[ist.itot*ist.nGridPoints]), ist.atot, &(psibe[ist.nHoles*ist.nGridPoints]), 
                      ist.nElecs, thetaZeta, nThetaZeta, ist.nGridPoints, par.dv);

    // Calculate TZeta for all excitonic states
    calcTZetaMatrix(TZeta, RckZeta, nThetaZeta, Cbs, ist);

    // Calculate ChiaiZetaHole (hole channel) for all excitonic states
    calcChiaiZetaMatrix(ChiaiZetaHole, RijZeta, nThetaZeta, Cbs, ist, 0);
    
    // Calculate ChiaiZetaElec (elec channel) for all excitonic states
    calcChiaiZetaMatrix(ChiaiZetaElec, RabZeta, nThetaZeta, Cbs, ist, 1);

    // Calculate Aai (hole channel) Aai = Ave over zeta (XaiZeta*TZeta) for all pairs of excitonic states
    calcAaiMatrix(AaiHole, ChiaiZetaHole, TZeta, nThetaZeta, ist.nIntExcitons, ist.nElecs, ist.itot);
    
    // Calculate Aai (elec channel) Aai = Ave over zeta (XaiZeta*TZeta) for all pairs of excitonic states
    calcAaiMatrix(AaiElec, ChiaiZetaElec, TZeta, nThetaZeta, ist.nIntExcitons, ist.atot, ist.nHoles);

    // AaiMatrices can be used to restart the calculation - rewrites the file each time
    // First line is the original seed entering calcRandomCoulombStates, nStochasticOrbitals performed so far and
    // the length of AaiHole and AaiElec which is needed to read them in properly when restarting here
    if (iZetanThreads < (nSO/2)) {
      sprintf(fileName, "AaiMatrix_%ld.dat", 0);
    }
    else {
      sprintf(fileName, "AaiMatrix_%ld.dat", 1);
    }
    pf = fopen(fileName, "w");
    fprintf(pf, "%ld %ld %ld %ld\n", iZetanThreads+nThetaZeta, ist.seed, 
                      ist.itot*ist.nElecs*ist.nIntExcitons*ist.nIntExcitons, 
                      ist.atot*ist.nHoles*ist.nIntExcitons*ist.nIntExcitons);
    // Print current values for restarting purposes
    for (i = 0; i < ist.itot*ist.nElecs*ist.nIntExcitons*ist.nIntExcitons; i++) {
      fprintf(pf, "%.16g %.16g\n", AaiHole[i].re, AaiHole[i].im);
    }
    for (i = 0; i < ist.atot*ist.nHoles*ist.nIntExcitons*ist.nIntExcitons; i++) {
      fprintf(pf, "%.16g %.16g\n", AaiElec[i].re, AaiElec[i].im);
    }
    fclose(pf);
    // Quick and dirty way to test if Calculating two sets of Aai Matrices is better
    if ((iZetanThreads == (nSO/2)) && (fabs(AaiHole1[i].re) < EPS)) { // Copy Aai to Aai1
      for (i = 0; i < ist.itot*ist.nElecs*ist.nIntExcitons*ist.nIntExcitons; i++) {
        AaiHole1[i].re = ( AaiHole[i].re / ((double)(nSO/2)) );
        AaiHole1[i].im = ( AaiHole[i].im / ((double)(nSO/2)) ); 
        AaiHole[i].re = 0.0;
        AaiHole[i].im = 0.0;
      }
      for (i = 0; i < ist.atot*ist.nHoles*ist.nIntExcitons*ist.nIntExcitons; i++) {
        AaiElec1[i].re = ( AaiElec[i].re / ((double)(nSO/2)) );
        AaiElec1[i].im = ( AaiElec[i].im / ((double)(nSO/2)) ); 
        AaiElec[i].re = 0.0;
        AaiElec[i].im = 0.0;
      }
      printf("Set AaiHole1 and AaiElec1 to Aai with iZetanThreads = %ld\n", iZetanThreads); fflush(stdout);
    }
  }
  for (i = 0; i < ist.itot*ist.nElecs*ist.nIntExcitons*ist.nIntExcitons; i++) {
    AaiHole2[i].re = ( AaiHole[i].re / ((double)(nSO/2)) );
    AaiHole2[i].im = ( AaiHole[i].im / ((double)(nSO/2)) ); 
  }
  for (i = 0; i < ist.atot*ist.nHoles*ist.nIntExcitons*ist.nIntExcitons; i++) {
    AaiElec2[i].re = ( AaiElec[i].re / ((double)(nSO/2)) );
    AaiElec2[i].im = ( AaiElec[i].im / ((double)(nSO/2)) ); 
  }      
  printf("Set AaiHole2 and AaiElec2 to Aai with iZetanThreads = %ld\n", iZetanThreads); fflush(stdout);
  writeCurrentTime(stdout);
  fprintf(stdout, "Finished the calculations of the RrsZeta, TZeta, ChiaiZeta and Aai Matrices\n");
  fflush(stdout);

  // Calculate the AR Lifetimes (energy conservation enforced here)
  long nAIHole = ist.nElecs*ist.itot;
  long nAIElec = ist.atot*ist.nHoles;
  long anIHole, anIElec;
  pf = fopen("doublyStochasticIntARLifetimes.dat", "w");
  for (iDeltaE = 0; iDeltaE < ist.nConsWindows; iDeltaE++) {
    iBiexc = 0;
    for (iExc1 = 0; iExc1 < ist.nIntExcitons; iExc1++) {
      for (iExc2 = 0; iExc2 < ist.nIntExcitons; iExc2++) {
        if ((Ebs[iExc1] + Ebs[iExc2]) < (par.maxInitE+EPS)) {
          iBiexcDeltaE = iBiexc + iDeltaE*ist.nBiexcitons; 
          // Hole channel
          for (i = 0; i < nSI; i++) { 
            if (thetaIList[i + nSI*iBiexcDeltaE] != -1) {
              anIHole = thetaIElecList[i + nSI*iBiexcDeltaE]*ist.itot;
              // Use two separate Aai's in order to decrease the bias
              hRate[iBiexcDeltaE] += (AaiHole1[(iExc1*ist.nIntExcitons+iExc2)*nAIHole + anIHole + thetaIList[i + nSI*iBiexcDeltaE]].re
                                        * AaiHole2[(iExc1*ist.nIntExcitons+iExc2)*nAIHole + anIHole + thetaIList[i + nSI*iBiexcDeltaE]].re)
                                    + (AaiHole1[(iExc1*ist.nIntExcitons+iExc2)*nAIHole + anIHole + thetaIList[i + nSI*iBiexcDeltaE]].im
                                        * AaiHole2[(iExc1*ist.nIntExcitons+iExc2)*nAIHole + anIHole + thetaIList[i + nSI*iBiexcDeltaE]].im);             
            }
          }
          // Elec channel
          for (a = 0; a < nSA; a++) { 
            if (thetaAList[a + nSA*iBiexcDeltaE] != -1) {
              anIElec = (thetaAList[a + nSA*iBiexcDeltaE]-ist.itot)*ist.nHoles;
              // Use two separate Aai's in order to decrease the bias
              eRate[iBiexcDeltaE] += (AaiElec1[(iExc1*ist.nIntExcitons+iExc2)*nAIElec + anIElec + thetaAHoleList[a + nSA*iBiexcDeltaE]].re
                                        * AaiElec2[(iExc1*ist.nIntExcitons+iExc2)*nAIElec + anIElec + thetaAHoleList[a + nSA*iBiexcDeltaE]].re)
                                    + (AaiElec1[(iExc1*ist.nIntExcitons+iExc2)*nAIElec + anIElec + thetaAHoleList[a + nSA*iBiexcDeltaE]].im
                                        * AaiElec2[(iExc1*ist.nIntExcitons+iExc2)*nAIElec + anIElec + thetaAHoleList[a + nSA*iBiexcDeltaE]].im);              
            }
          }
          // Scale the hRate and eRate and convert them to units of (ps)^-1
          hRate[iBiexcDeltaE] *= (nHotHnSIRatio[iBiexcDeltaE] * AUTOPS * TWOPI / (2.0 * deltaE[iDeltaE]));
          eRate[iBiexcDeltaE] *= (nHotEnSARatio[iBiexcDeltaE] * AUTOPS * TWOPI / (2.0 * deltaE[iDeltaE]));
          arRate[iBiexcDeltaE] = (hRate[iBiexcDeltaE] + eRate[iBiexcDeltaE]);
          // Print out AR lifetimes for each initial biexcitonic state for the current energy window
          fprintf(pf, "%ld %.6f %ld %.6f ", iDeltaE, deltaE[iDeltaE], iBiexc, biexcitonEnergies[iBiexc]);
          fprintf(pf, "%.8f %.8f ", hRate[iBiexcDeltaE], eRate[iBiexcDeltaE]);
          fprintf(pf, "%.8f %.8f\n", arRate[iBiexcDeltaE], 1.0/(arRate[iBiexcDeltaE]+EPS));
          iBiexc++;
        }
      }
    }
    writeSeparation(pf);  
  }
  fclose(pf);

  // Calculate the Boltzmann weighted lifetimes for each energy conservation window
  pf = fopen("IntARLifetimesDoublyStochastic.dat", "w");
  for (iDeltaE = 0; iDeltaE < ist.nConsWindows; iDeltaE++) {
    fprintf(pf, "Energy conservation window = %.6f\n\n", deltaE[iDeltaE]);
    calcBoltzmannWeightedRates(biexcitonEnergies, &(eRate[iDeltaE*ist.nBiexcitons]), 
            &(hRate[iDeltaE*ist.nBiexcitons]), par.temp, ist.nBiexcitons, pf); 
    writeSeparation(pf);
  }
  fclose(pf);

  // Write ending of function
  fprintf(stdout, "\nFinished the calculation of the interacting AR lifetimes using the doubly stochastic method\n");
  writeCurrentTime(stdout);
  fflush(stdout);

  // Free dynamically allocated memory
  free(biexcitonEnergies);
  free(eRate); free(hRate); free(arRate);
  free(AaiHole); free(AaiElec);
  free(TZeta); free(ChiaiZetaHole); free(ChiaiZetaElec);
  free(RckZeta); free(RijZeta); free(RabZeta);
  free(thetaZeta);
  free(nHotHoleForElec); free(nHotElecForHole);
  free(elecHotHoleList); free(holeHotElecList);
  free(thetaIList); free(thetaAList);
  free(thetaIElecList); free(thetaAHoleList);
  free(nHotHnSIRatio); free(nHotEnSARatio);
  free(globalThetaIElecList); free(globalThetaAHoleList);
  free(globalThetaAList); free(globalThetaIList);
  free(AaiHole2); free(AaiElec2);
  free(AaiHole1); free(AaiElec1);

  return;
}

/*****************************************************************************/

void calcStochasticCoulombIntAR(double *Cbs, double *Ebs, double *psibe, double *evalbe, double *psiai, double *evalai,
                                     zomplex *potq, double *vx, double *vy, double *vz, lng_st ist, par_st par, 
                                    fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwpsi)
{
  FILE *pf;
  long a, i, nSO;
  long iDeltaE, iBiexc, iExc1, iExc2, iBiexcDeltaE, iZeta, iZetanThreads, nThetaZeta;
  double *eRate, *hRate, *arRate, *biexcitonEnergies;
  double deltaE[ist.nConsWindows];  
  char fileName[50];
  zomplex *RckZeta, *RabZeta, *RijZeta, *TZeta;
  zomplex *ChiaiZetaHole, *ChiaiZetaElec;
  zomplex *AaiHole, *AaiElec, *AaiHole1, *AaiElec1, *AaiHole2, *AaiElec2;
  zomplex *thetaZeta;

  // Useful integers 
  nSO = ist.nStochOrbitals;  
  const long THL = ist.nNonIntExcitons;
  const long TH2L = THL * ist.nHoles;
  const long THL2 = THL * ist.nElecs;

  // Write beginning of function
  writeSeparation(stdout);
  writeCurrentTime(stdout);
  fprintf(stdout, "Beginnning the calculation of the interacting AR lifetimes by stochastically representing the Coulomb operator\n\n");
  fprintf(stdout, "The number of energy conservation windows = %d\n", ist.nConsWindows);
  fprintf(stdout, "The largest energy conservation window = %.6f %.6f\n", par.maxDeltaE, par.maxDeltaE*AUTOEV);
  fprintf(stdout, "The number of deterministic hot hole states = %d\n", ist.itot);
  fprintf(stdout, "The number of deterministic hot elec states = %d\n", ist.atot);
  fprintf(stdout, "The number of stochastic orbitals used to approximate the Coulomb operator = %ld\n", nSO);
  fflush(stdout);

  // Dynamically allocate memory
  if ((thetaZeta = (zomplex *) calloc(ist.nGridPoints*ist.nThreads, sizeof(zomplex))) == NULL) nerror("thetaZeta");
  if ((RckZeta = (zomplex *) calloc(ist.nElecs*ist.nHoles*ist.nThreads, sizeof(zomplex))) == NULL) nerror("RckZeta");  
  if ((TZeta = (zomplex *) calloc(ist.nIntExcitons*ist.nThreads, sizeof(zomplex))) == NULL) nerror("TZeta");  
  if ((RijZeta = (zomplex *) calloc(ist.itot*ist.nHoles*ist.nThreads, sizeof(zomplex))) == NULL) nerror("RijZeta");  
  if ((ChiaiZetaHole = (zomplex *) calloc(ist.itot*ist.nElecs*ist.nIntExcitons*ist.nThreads, sizeof(zomplex))) == NULL) nerror("ChiaiZetaHole");  
  if ((RabZeta = (zomplex *) calloc(ist.atot*ist.nElecs*ist.nThreads, sizeof(zomplex))) == NULL) nerror("RabZeta");  
  if ((ChiaiZetaElec = (zomplex *) calloc(ist.atot*ist.nHoles*ist.nIntExcitons*ist.nThreads, sizeof(zomplex))) == NULL) nerror("ChiaiZetaElec");
  if ((AaiHole = (zomplex *) calloc(ist.itot*ist.nElecs*ist.nIntExcitons*ist.nIntExcitons, sizeof(zomplex))) == NULL) nerror("AaiHole");  
  if ((AaiElec = (zomplex *) calloc(ist.atot*ist.nHoles*ist.nIntExcitons*ist.nIntExcitons, sizeof(zomplex))) == NULL) nerror("AaiElec");  
  if ((AaiHole1 = (zomplex *) calloc(ist.itot*ist.nElecs*ist.nIntExcitons*ist.nIntExcitons, sizeof(zomplex))) == NULL) nerror("AaiHole1");  
  if ((AaiElec1 = (zomplex *) calloc(ist.atot*ist.nHoles*ist.nIntExcitons*ist.nIntExcitons, sizeof(zomplex))) == NULL) nerror("AaiElec1");  
  if ((AaiHole2 = (zomplex *) calloc(ist.itot*ist.nElecs*ist.nIntExcitons*ist.nIntExcitons, sizeof(zomplex))) == NULL) nerror("AaiHole2");  
  if ((AaiElec2 = (zomplex *) calloc(ist.atot*ist.nHoles*ist.nIntExcitons*ist.nIntExcitons, sizeof(zomplex))) == NULL) nerror("AaiElec2"); 
  if ((eRate = (double *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(double))) == NULL) nerror("eRate");
  if ((hRate = (double *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(double))) == NULL) nerror("hRate");
  if ((arRate = (double *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(double))) == NULL) nerror("arRate");
  if ((biexcitonEnergies = (double *) calloc(ist.nBiexcitons, sizeof(double))) == NULL) nerror("biexcitonEnergies");

  // Create multiple energy conservation windows so one does not have to repeat the entire calculation
  fillEnergyConservationWindows(deltaE, par.maxDeltaE, ist.nConsWindows);

  // Fill in biexcitonEnergies
  fillIntBiexcEnergies(biexcitonEnergies, Ebs, ist, par);

  // Check if previously calculated AaiMatrices should be read in
  long iInitZetanThreads, nAaiHoleMatElements, nAaiElecMatElements;
  if (ist.readAaiMatrices) {
    for (a = 0; a < 2; a++) {
      sprintf(fileName, "AaiMatrix_%ld.par", a);
      if ( access(fileName, F_OK) != -1 ) {
        pf = fopen(fileName, "r");
        fscanf(pf, "%ld %ld %ld %ld", &iInitZetanThreads, &(ist.seed), &nAaiHoleMatElements, &nAaiElecMatElements);
        printf("Reading in AaiMatrix%ld\n", a);
        printf("Number of previously calculated stochastic orbitals = %ld\n", iInitZetanThreads);
        ist.seed *= -1;
        printf("Seed set to = %ld\n", ist.seed); 
        fflush(stdout);
        for (i = 0; i < nAaiHoleMatElements; i++) {
          fscanf(pf, "%lg %lg", &(AaiHole[i].re), &(AaiHole[i].im));
        }
        for (i = 0; i < nAaiElecMatElements; i++) {
          fscanf(pf, "%lg %lg", &(AaiElec[i].re), &(AaiElec[i].im));
        }
        fclose(pf);
        if (! a && (iInitZetanThreads == (nSO/2))) {
          for (i = 0; i < nAaiHoleMatElements; i++) {
            AaiHole1[i].re = AaiHole[i].re;
            AaiHole1[i].im = AaiHole[i].im;
          }
          for (i = 0; i < nAaiElecMatElements; i++) {
            AaiElec1[i].re = AaiElec[i].re;
            AaiElec1[i].im = AaiElec[i].im;
          }
        }
      }
      else {
        printf("\nNo %s file detected in current working directory.\n", fileName); 
        if (! a) {
          iInitZetanThreads = 0;
          break; // leave loop, do not try to read AaiMatrix_1.par
        }
        printf("The number of stochastic orbitals previously calculated = %ld\n", iInitZetanThreads); 
        fflush(stdout);
      }
    }
  }
  else {
    iInitZetanThreads = 0;
  }
  if (iInitZetanThreads >= (nSO/2)) {
    for (i = 0; i < ist.itot*ist.nElecs*ist.nIntExcitons*ist.nIntExcitons; i++) {
      AaiHole1[i].re /= ((double)(nSO/2));
      AaiHole1[i].im /= ((double)(nSO/2)); 
    }
    for (i = 0; i < ist.atot*ist.nHoles*ist.nIntExcitons*ist.nIntExcitons; i++) {
      AaiElec1[i].re /= ((double)(nSO/2));
      AaiElec1[i].im /= ((double)(nSO/2)); 
    }            
  }
  if (iInitZetanThreads == (nSO/2)) {
    for (i = 0; i < ist.itot*ist.nElecs*ist.nIntExcitons*ist.nIntExcitons; i++) {
      AaiHole[i].re = 0.0;
      AaiHole[i].im = 0.0;
    }
    for (i = 0; i < ist.atot*ist.nHoles*ist.nIntExcitons*ist.nIntExcitons; i++) {
      AaiElec[i].re = 0.0;
      AaiElec[i].im = 0.0;
    }          
  }

  // Calculate stochastic orbitals that represent the Coulomb operator and their overlaps with the quasiparticle states
  fprintf(stdout, "\nBeginning the calculations of the RrsZeta, TZeta, ChiaiZeta and Aai Matrices\n");
  writeCurrentTime(stdout);
  fprintf(stdout, "The number of 3D integrals being performed = %ld\n", 
                    (nSO-iInitZetanThreads)*(ist.nElecs*ist.nHoles+ist.itot*ist.nElecs+ist.atot*ist.nHoles));
  fflush(stdout);
  
  // Begin stochastic Coulomb orbital calculation 
  for (iZetanThreads = iInitZetanThreads; iZetanThreads < nSO; iZetanThreads += ist.nThreads) {
    if ( ((iZetanThreads + ist.nThreads) > nSO) && (nSO % ist.nThreads)) {
      nThetaZeta = (nSO % ist.nThreads);
    }
    else {
      nThetaZeta = ist.nThreads;
    }
    // Generate stochastic orbitals used to approximate the Coulomb operator (the thetaZeta orbitals)
    ist.seed = calcRandomCoulombStates(thetaZeta, nThetaZeta, potq, ist, par, planfw, planbw, fftwpsi);

    // Calculate RckZeta (elec and hole) 
    calcRrsZetaMatrix(RckZeta, &(psibe[ist.nHoles*ist.nGridPoints]), ist.nElecs, psibe, ist.nHoles,
                      thetaZeta, nThetaZeta, ist.nGridPoints, par.dv);
    for (i = 0; i < ist.nElecs*ist.nHoles*nThetaZeta; i++) { RckZeta[i].im *= -1.0; }
    
    // Calculate RijZeta (hot hole and hole - hole channel) 
    calcRrsZetaMatrix(RijZeta, psiai, ist.itot, psibe, ist.nHoles, 
                      thetaZeta, nThetaZeta, ist.nGridPoints, par.dv);
    
    // Calculate RabZeta (hot elec and elec - elec channel) 
    calcRrsZetaMatrix(RabZeta, &(psiai[ist.itot*ist.nGridPoints]), ist.atot, &(psibe[ist.nHoles*ist.nGridPoints]), 
                      ist.nElecs, thetaZeta, nThetaZeta, ist.nGridPoints, par.dv);

    // Calculate TZeta for all excitonic states
    calcTZetaMatrix(TZeta, RckZeta, nThetaZeta, Cbs, ist);

    // Calculate ChiaiZetaHole (hole channel) for all excitonic states
    calcChiaiZetaMatrix(ChiaiZetaHole, RijZeta, nThetaZeta, Cbs, ist, 0);
    
    // Calculate ChiaiZetaElec (elec channel) for all excitonic states
    calcChiaiZetaMatrix(ChiaiZetaElec, RabZeta, nThetaZeta, Cbs, ist, 1);

    // Calculate Aai (hole channel) Aai = Ave over zeta (XaiZeta*TZeta) for all pairs of excitonic states
    calcAaiMatrix(AaiHole, ChiaiZetaHole, TZeta, nThetaZeta, ist.nIntExcitons, ist.nElecs, ist.itot);
    
    // Calculate Aai (elec channel) Aai = Ave over zeta (XaiZeta*TZeta) for all pairs of excitonic states
    calcAaiMatrix(AaiElec, ChiaiZetaElec, TZeta, nThetaZeta, ist.nIntExcitons, ist.atot, ist.nHoles);

    // AaiMatrices can be used to restart the calculation - rewrites the file each time
    // First line is the original seed entering calcRandomCoulombStates, nStochasticOrbitals performed so far and
    // the length of AaiHole and AaiElec which is needed to read them in properly when restarting here
    if (iZetanThreads < (nSO/2)) {
      sprintf(fileName, "AaiMatrix_%ld.dat", 0);
    }
    else {
      sprintf(fileName, "AaiMatrix_%ld.dat", 1);
    }
    pf = fopen(fileName, "w");
    fprintf(pf, "%ld %ld %ld %ld\n", iZetanThreads+nThetaZeta, ist.seed, 
                      ist.itot*ist.nElecs*ist.nIntExcitons*ist.nIntExcitons, 
                      ist.atot*ist.nHoles*ist.nIntExcitons*ist.nIntExcitons);
    // Print current values for restarting purposes
    for (i = 0; i < ist.itot*ist.nElecs*ist.nIntExcitons*ist.nIntExcitons; i++) {
      fprintf(pf, "%.16g %.16g\n", AaiHole[i].re, AaiHole[i].im);
    }
    for (i = 0; i < ist.atot*ist.nHoles*ist.nIntExcitons*ist.nIntExcitons; i++) {
      fprintf(pf, "%.16g %.16g\n", AaiElec[i].re, AaiElec[i].im);
    }
    fclose(pf);
    // Quick and dirty way to test if Calculating two sets of Aai Matrices is better
    if ((iZetanThreads == (nSO/2)) && (fabs(AaiHole1[i].re) < EPS)) { // Copy Aai to Aai1
      for (i = 0; i < ist.itot*ist.nElecs*ist.nIntExcitons*ist.nIntExcitons; i++) {
        AaiHole1[i].re = ( AaiHole[i].re / ((double)(nSO/2)) );
        AaiHole1[i].im = ( AaiHole[i].im / ((double)(nSO/2)) ); 
        AaiHole[i].re = 0.0;
        AaiHole[i].im = 0.0;
      }
      for (i = 0; i < ist.atot*ist.nHoles*ist.nIntExcitons*ist.nIntExcitons; i++) {
        AaiElec1[i].re = ( AaiElec[i].re / ((double)(nSO/2)) );
        AaiElec1[i].im = ( AaiElec[i].im / ((double)(nSO/2)) ); 
        AaiElec[i].re = 0.0;
        AaiElec[i].im = 0.0;
      }
      printf("Set AaiHole1 and AaiElec1 to Aai with iZetanThreads = %ld\n", iZetanThreads); fflush(stdout);
    }
  }
  for (i = 0; i < ist.itot*ist.nElecs*ist.nIntExcitons*ist.nIntExcitons; i++) {
    AaiHole2[i].re = ( AaiHole[i].re / ((double)(nSO/2)) );
    AaiHole2[i].im = ( AaiHole[i].im / ((double)(nSO/2)) ); 
  }
  for (i = 0; i < ist.atot*ist.nHoles*ist.nIntExcitons*ist.nIntExcitons; i++) {
    AaiElec2[i].re = ( AaiElec[i].re / ((double)(nSO/2)) );
    AaiElec2[i].im = ( AaiElec[i].im / ((double)(nSO/2)) ); 
  }      
  printf("Set AaiHole2 and AaiElec2 to Aai with iZetanThreads = %ld\n", iZetanThreads); fflush(stdout);
  writeCurrentTime(stdout);
  fprintf(stdout, "Finished the calculations of the RrsZeta, TZeta, ChiaiZeta and Aai Matrices\n");
  fflush(stdout);

  // Calculate the AR Lifetimes (energy conservation enforced here)
  long nAIHole = ist.nElecs*ist.itot;
  long nAIElec = ist.atot*ist.nHoles;
  long anIHole, anIElec;
  pf = fopen("stochasticCoulombIntARLifetimes.dat", "w");
  for (iDeltaE = 0; iDeltaE < ist.nConsWindows; iDeltaE++) {
    iBiexc = 0;
    for (iExc1 = 0; iExc1 < ist.nIntExcitons; iExc1++) {
      for (iExc2 = 0; iExc2 < ist.nIntExcitons; iExc2++) {
        if ((Ebs[iExc1] + Ebs[iExc2]) < (par.maxInitE+EPS)) {
          iBiexcDeltaE = iBiexc + iDeltaE*ist.nBiexcitons; 
          // Hole channel
          for (a = 0; a < ist.nElecs; a++) {
            anIHole = a*ist.itot;
            for (i = 0; i < ist.itot; i++) {
              if (fabs(evalbe[a+ist.nHoles] - evalai[i] - biexcitonEnergies[iBiexc]) <= deltaE[iDeltaE]) {
                // Use two separate Aai's in order to decrease the bias
                hRate[iBiexcDeltaE] += (AaiHole1[(iExc1*ist.nIntExcitons+iExc2)*nAIHole + anIHole + i].re*AaiHole2[(iExc1*ist.nIntExcitons+iExc2)*nAIHole + anIHole + i].re)
                                      + (AaiHole1[(iExc1*ist.nIntExcitons+iExc2)*nAIHole + anIHole + i].im*AaiHole2[(iExc1*ist.nIntExcitons+iExc2)*nAIHole + anIHole + i].im);
              }
            }
          }
          // Elec channel
          for (a = 0; a < ist.atot; a++) {
            anIElec = a*ist.nHoles;
            for (i = 0; i < ist.nHoles; i++) {
              if (fabs(evalai[a + ist.itot] - evalbe[i] - biexcitonEnergies[iBiexc]) <= deltaE[iDeltaE]) {
                // Use two separate Aai's in order to decrease the bias
                eRate[iBiexcDeltaE] += (AaiElec1[(iExc1*ist.nIntExcitons+iExc2)*nAIElec + anIElec + i].re*AaiElec2[(iExc1*ist.nIntExcitons+iExc2)*nAIElec + anIElec + i].re)
                                      + (AaiElec1[(iExc1*ist.nIntExcitons+iExc2)*nAIElec + anIElec + i].im*AaiElec2[(iExc1*ist.nIntExcitons+iExc2)*nAIElec + anIElec + i].im);
              }
            }
          }
          // Scale the hRate and eRate and convert them to units of (ps)^-1
          hRate[iBiexcDeltaE] *= (AUTOPS * TWOPI / (2.0 * deltaE[iDeltaE]));
          eRate[iBiexcDeltaE] *= (AUTOPS * TWOPI / (2.0 * deltaE[iDeltaE]));
          arRate[iBiexcDeltaE] = (hRate[iBiexcDeltaE] + eRate[iBiexcDeltaE]);
          // Print out AR lifetimes for each initial biexcitonic state for the current energy window
          fprintf(pf, "%ld %.6f %ld %.6f ", iDeltaE, deltaE[iDeltaE], iBiexc, biexcitonEnergies[iBiexc]);
          fprintf(pf, "%.8f %.8f ", hRate[iBiexcDeltaE], eRate[iBiexcDeltaE]);
          fprintf(pf, "%.8f %.8f\n", arRate[iBiexcDeltaE], 1.0/(arRate[iBiexcDeltaE]+EPS));
          iBiexc++;
        }
      }
    }
    writeSeparation(pf);  
  }
  fclose(pf);

  // Calculate the Boltzmann weighted lifetimes for each energy conservation window
  pf = fopen("IntARLifetimesStochasticCoulomb.dat", "w");
  for (iDeltaE = 0; iDeltaE < ist.nConsWindows; iDeltaE++) {
    fprintf(pf, "Energy conservation window = %.6f\n\n", deltaE[iDeltaE]);
    calcBoltzmannWeightedRates(biexcitonEnergies, &(eRate[iDeltaE*ist.nBiexcitons]), 
            &(hRate[iDeltaE*ist.nBiexcitons]), par.temp, ist.nBiexcitons, pf); // TODO: remove -1 from ist.nBiexcitons
    writeSeparation(pf);
  }
  fclose(pf);

  // Write ending of function
  fprintf(stdout, "\nFinished the calculation of the interacting AR lifetimes by stochastically representing the Coulomb operator\n");
  writeCurrentTime(stdout);
  fflush(stdout);

  // Free dynamically allocated memory
  free(biexcitonEnergies);
  free(eRate); free(hRate); free(arRate);
  free(AaiHole); free(AaiElec);
  free(TZeta); free(ChiaiZetaHole); free(ChiaiZetaElec);
  free(RckZeta); free(RijZeta); free(RabZeta);
  free(thetaZeta);
  free(AaiHole2); free(AaiElec2);
  free(AaiHole1); free(AaiElec1);

  return;
}


/*****************************************************************************/

void calcStochasticFinalStatesIntAR(double *Cbs, double *Ebs, double *vijck, double *vabck, double *psibe, double *evalbe, 
                                    double *psiai, double *evalai, zomplex *potq, double *poth, lng_st ist, par_st par, 
                                    fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwpsi)
{
  FILE *pf, *fp;
  long a, i, b, c, j, k, nSA, nSI;
  long iDeltaE, iBiexc, iExc1, iExc2, iBiexcDeltaE, idum, counter;
  long *globalThetaIList, *globalThetaAList, *globalThetaIElecList, *globalThetaAHoleList;
  long *elecHotHoleList, *holeHotElecList, *nHotElecForHole, *nHotHoleForElec;
  long randomHoleIndex, randomElecIndex, randomHotHoleIndex, randomHotElecIndex;
  long *thetaIList, *thetaIElecList, *thetaAList, *thetaAHoleList;
  double sum, tmp, *nHotHnSIRatio, *nHotEnSARatio, *w2h, *w2e;
  double *eRate, *hRate, *arRate, *biexcitonEnergies;
  double deltaE[ist.nConsWindows]; 

  // Useful integers
  idum = ist.seed; 
  nSA = ist.nStochHotElecs;  
  nSI = ist.nStochHotHoles; 
  const long THL = ist.nNonIntExcitons;
  const long TH2L = THL * ist.nHoles;
  const long THL2 = THL * ist.nElecs;

  // Write beginning of function
  writeSeparation(stdout);
  writeCurrentTime(stdout);
  fprintf(stdout, "Beginnning the calculation of the interacting AR lifetimes by sampling the final states\n\n");
  fprintf(stdout, "The number of energy conservation windows = %d\n", ist.nConsWindows);
  fprintf(stdout, "The largest energy conservation windows = %.6f %.6f\n", par.maxDeltaE, par.maxDeltaE*AUTOEV);
  fprintf(stdout, "The maximum number of deterministic hot hole states = %d\n", ist.itot);
  fprintf(stdout, "The maximum number of deterministic hot elec states = %d\n", ist.atot);
  fprintf(stdout, "The number of stochastic hot hole states sampled = %ld\n", nSI);
  fprintf(stdout, "The number of stochastic hot elec states sampled = %ld\n", nSA);
  fflush(stdout);

  // Dynamically allocate memory
  if ((globalThetaIList = (long *) calloc(nSI, sizeof(long))) == NULL) nerror("globalThetaIList");
  if ((globalThetaIElecList = (long *) calloc(nSI, sizeof(long))) == NULL) nerror("globalThetaIElecList");
  if ((globalThetaAList = (long *) calloc(nSA, sizeof(long))) == NULL) nerror("globalThetaAList");
  if ((globalThetaAHoleList = (long *) calloc(nSA, sizeof(long))) == NULL) nerror("globalThetaAHoleList");
  if ((elecHotHoleList = (long *) calloc(ist.nBiexcitons*ist.nConsWindows*ist.nElecs*ist.itot, sizeof(long))) == NULL) nerror("elecHotHoleList");
  if ((nHotHoleForElec = (long *) calloc(ist.nBiexcitons*ist.nConsWindows*ist.nElecs, sizeof(long))) == NULL) nerror("nHotHoleForElec");
  if ((holeHotElecList = (long *) calloc(ist.nBiexcitons*ist.nConsWindows*ist.nHoles*ist.atot, sizeof(long))) == NULL) nerror("holeHotElecList");
  if ((nHotElecForHole = (long *) calloc(ist.nBiexcitons*ist.nConsWindows*ist.nHoles, sizeof(long))) == NULL) nerror("nHotElecForHole");
  if ((thetaIList = (long *) calloc(nSI*ist.nBiexcitons*ist.nConsWindows, sizeof(long))) == NULL) nerror("thetaIList");
  if ((thetaIElecList = (long *) calloc(nSI*ist.nBiexcitons*ist.nConsWindows, sizeof(long))) == NULL) nerror("thetaIElecList");
  if ((thetaAList = (long *) calloc(nSA*ist.nBiexcitons*ist.nConsWindows, sizeof(long))) == NULL) nerror("thetaAList");
  if ((thetaAHoleList = (long *) calloc(nSA*ist.nBiexcitons*ist.nConsWindows, sizeof(long))) == NULL) nerror("thetaAHoleList");
  if ((eRate = (double *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(double))) == NULL) nerror("eRate");
  if ((hRate = (double *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(double))) == NULL) nerror("hRate");
  if ((arRate = (double *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(double))) == NULL) nerror("arRate");
  if ((nHotHnSIRatio = (double *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(double))) == NULL) nerror("nHotHnSIRatio");
  if ((nHotEnSARatio = (double *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(double))) == NULL) nerror("nHotEnSARatio");
  if ((biexcitonEnergies = (double *) calloc(ist.nBiexcitons, sizeof(double))) == NULL) nerror("biexcitonEnergies");

  // Create multiple energy conservation windows so one does not have to repeat the entire calculation
  fillEnergyConservationWindows(deltaE, par.maxDeltaE, ist.nConsWindows);

  // Fill in biexcitonEnergies
  fillIntBiexcEnergies(biexcitonEnergies, Ebs, ist, par);
  
  // Create lists and counts of the allowed final states for each of the possible initial states and energyConsWindow
  storeAllowedFinalStatesIntAR(elecHotHoleList, nHotHoleForElec, holeHotElecList, nHotElecForHole, biexcitonEnergies, 
                                Ebs, evalbe, evalai, deltaE, ist, par);

  // Generate the stochastic orbitals by:
  // 1) randomly select a band-edge state 
  // 2) randomly select a hot carrier from the list of hot carriers such that the elec-hole pair
  //    satisfies energy conservation with at least one initial biexcitonic state
  counter = 0;
  for (i = 0; i < nSI; ) {
    randomElecIndex = (long)((double)(ist.nElecs)*ran_nrc(&idum)) + ist.nHoles; 
    randomHotHoleIndex = (long)((double)(ist.itot)*ran_nrc(&idum));
    counter++;
    if (ist.itot == 0 || counter > nSI*100) break; // test to leave loop if no final states 
    // Make sure this final state satisfies energy conservation for at least 1 initial state
    for (iBiexc = 0; iBiexc < ist.nBiexcitons; iBiexc++) {
      if (fabs(evalbe[randomElecIndex] - evalai[randomHotHoleIndex] - biexcitonEnergies[iBiexc]) <= par.maxDeltaE) {
        globalThetaIElecList[i] = randomElecIndex; // nHoles to nHoles+nElecs
        globalThetaIList[i] = randomHotHoleIndex;  // 0 to itot
        i++;
        break;
      }
    }
  }
  counter = 0;
  for (a = 0; a < nSA; ) {
    randomHoleIndex = (long)((double)(ist.nHoles)*ran_nrc(&idum));
    randomHotElecIndex = (long)((double)(ist.atot)*ran_nrc(&idum)) + ist.itot;
    counter++;
    if (ist.atot == 0 || counter > nSA*100) break; // test to leave loop if no final states
    // Make sure this final state satisfies energy conservation for at least 1 initial state
    for (iBiexc = 0; iBiexc < ist.nBiexcitons; iBiexc++) {
      if (fabs(evalai[randomHotElecIndex] - evalbe[randomHoleIndex] - biexcitonEnergies[iBiexc]) <= par.maxDeltaE) {
        globalThetaAHoleList[a] = randomHoleIndex; // 0 to nHoles
        globalThetaAList[a] = randomHotElecIndex;  // itot to itot+atot
        a++;
        break;
      }
    }
  }

  // Fill thetaIList and thetaAList based on the globalTheta Lists (energy conservation check here)
  pf = fopen("sHotElecStatesIntAR.dat", "w");
  fp = fopen("sHotHoleStatesIntAR.dat", "w");
  for (iDeltaE = 0; iDeltaE < ist.nConsWindows; iDeltaE++) {
    for (iBiexc = 0; iBiexc < ist.nBiexcitons; iBiexc++) {
      iBiexcDeltaE = iBiexc + iDeltaE*ist.nBiexcitons;
      // Hole channel
      for (i = 0; i < nSI; i++) {
        if (fabs(evalbe[globalThetaIElecList[i]] - evalai[globalThetaIList[i]] - biexcitonEnergies[iBiexc]) <= deltaE[iDeltaE]) {
          thetaIElecList[i + nSI*iBiexcDeltaE] = (globalThetaIElecList[i] - ist.nHoles); // only to be used in w2h indexing, 0 to nElecs
          thetaIList[i + nSI*iBiexcDeltaE] = globalThetaIList[i]; // 0 to itot
        }
        else { // used to signify it does not satisfy energy conservation
          thetaIElecList[i + nSI*iBiexcDeltaE] = -1; 
          thetaIList[i + nSI*iBiexcDeltaE] = -1;
        }
        // Print out information to see if maps are working as desired
        fprintf(fp, "%5d %.6f %4d %.6f ", i, deltaE[iDeltaE], iBiexc, biexcitonEnergies[iBiexc]);
        fprintf(fp, "%5d % .6f ", thetaIList[i + nSI*iBiexcDeltaE], evalai[globalThetaIList[i]]);
        fprintf(fp, "%5d % .6f ", thetaIElecList[i + nSI*iBiexcDeltaE] , evalbe[globalThetaIElecList[i]]); 
        fprintf(fp, "% .6f ", evalbe[globalThetaIElecList[i]] - evalai[globalThetaIList[i]]); 
        fprintf(fp, "% .6f\n", biexcitonEnergies[iBiexc] - (evalbe[globalThetaIElecList[i]] - evalai[globalThetaIList[i]]));
      }
      // Elec channel
      for (a = 0; a < nSA; a++) {
        if (fabs(evalai[globalThetaAList[a]] - evalbe[globalThetaAHoleList[a]] - biexcitonEnergies[iBiexc]) <= deltaE[iDeltaE]) {
          thetaAHoleList[a + nSA*iBiexcDeltaE] = globalThetaAHoleList[a]; // only to be used in w2e indexing, 0 to nHoles
          thetaAList[a + nSA*iBiexcDeltaE] = globalThetaAList[a];  // itot to itot+atot
        }
        else { // used to signify it does not satisfy energy conservation
          thetaAHoleList[a + nSA*iBiexcDeltaE] = -1; 
          thetaAList[a + nSA*iBiexcDeltaE] = -1;
        }
        // Print out information to see if maps are working as desired        
        fprintf(pf, "%5d %.6f %4d %.6f ", a, deltaE[iDeltaE], iBiexc, biexcitonEnergies[iBiexc]);
        fprintf(pf, "%5d % .6f ", thetaAHoleList[a + nSA*iBiexcDeltaE], evalbe[globalThetaAHoleList[a]]);
        fprintf(pf, "%5d % .6f ", thetaAList[a + nSA*iBiexcDeltaE], evalai[globalThetaAList[a]]); 
        fprintf(pf, "% .6f ", evalai[globalThetaAList[a]] - evalbe[globalThetaAHoleList[a]]); 
        fprintf(pf, "% .6f\n", biexcitonEnergies[iBiexc] - (evalai[globalThetaAList[a]] - evalbe[globalThetaAHoleList[a]]));  
      }
    }
  }
  fclose(fp);
  fclose(pf);

  // Determine the scaling for sampling final states
  double sumSto, sumDet;
  for (iDeltaE = 0; iDeltaE < ist.nConsWindows; iDeltaE++) {
    iBiexc = 0;
    for (iExc1 = 0; iExc1 < ist.nIntExcitons; iExc1++) {
      for (iExc2 = 0; iExc2 < ist.nIntExcitons; iExc2++) { 
        if ((Ebs[iExc1] + Ebs[iExc2]) < par.maxInitE+EPS) {
          iBiexcDeltaE = iBiexc + iDeltaE*ist.nBiexcitons;
          // Hole channel
          for (sumSto = 0.0, i = 0; i < nSI; i++) {
            a = thetaIElecList[i + nSI*iBiexcDeltaE]; // a goes from 0 to ist.nElecs
            if (a != -1) { // final states that do not satisfy energy conservation have a = i = -1
              for (tmp = 0.0, j = 0; j < ist.nHoles; j++) { 
                tmp += Cbs[iExc1*ist.nNonIntExcitons + a*ist.nHoles + j];
              }
              sumSto += tmp*tmp;
            }
          }
          for (sumDet = 0.0, a = 0; a < ist.nElecs; a++) { // sum over a,j for c_a,j, delta_ab
            for (tmp = 0.0, j = 0; j < ist.nHoles; j++) {
              tmp += Cbs[iExc1*ist.nNonIntExcitons + a*ist.nHoles + j];
            }
            sumDet += tmp*tmp*(double)(nHotHoleForElec[a + iBiexcDeltaE*ist.nElecs]); 
          }
          nHotHnSIRatio[iBiexcDeltaE] = sumDet/sumSto;
          // Elec channel
          for (sumSto = 0.0, a = 0; a < nSA; a++) {
            i = thetaAHoleList[a + nSA*iBiexcDeltaE]; // a goes from 0 to ist.nHoles
            if (i != -1) { // final states that do not satisfy energy conservation have i = a = -1
              for (tmp = 0.0, b = 0; b < ist.nElecs; b++) {
                tmp += Cbs[iExc1*ist.nNonIntExcitons + b*ist.nHoles + i];
              }
              sumSto += tmp*tmp;
            }
          }
          for (sumDet = 0.0, i = 0; i < ist.nHoles; i++) { // sum over b,i for c_b,i, delta_ij
            for (tmp = 0.0, b = 0; b < ist.nElecs; b++) {
              tmp += Cbs[iExc1*ist.nNonIntExcitons + b*ist.nHoles + i];
            }
            sumDet += tmp*tmp*(double)(nHotElecForHole[i + iBiexcDeltaE*ist.nHoles]);
          }
          nHotEnSARatio[iBiexcDeltaE] = sumDet/sumSto;          
          iBiexc++;
        }
      }
    }
  }

  // Remove states from psiai/evalai that are not being sampled 
  remStatesNotInThetaLists(psiai, evalai, thetaIList, thetaAList, &ist);

  // Realloc memory for the smaller psiai, evalai, vijck, vabck and allocate memory for coherent sums
  psiai = realloc(psiai, ist.nGridPoints*ist.nHotNonIntStates*sizeof(psiai[0]));
  evalai = realloc(evalai, ist.nHotNonIntStates*sizeof(evalai[0]));
  vijck = realloc(vijck, ist.nHoles*ist.nHoles*ist.nElecs*ist.itot*sizeof(vijck[0]));
  vabck = realloc(vabck, ist.nHoles*ist.nElecs*ist.nElecs*ist.atot*sizeof(vabck[0]));
  if ((w2h = (double *) calloc(ist.nElecs*ist.itot*ist.nBiexcitons, sizeof(double))) == NULL) nerror("w2h");
  if ((w2e = (double *) calloc(ist.atot*ist.nHoles*ist.nBiexcitons, sizeof(double))) == NULL) nerror("w2e");

  // Calculate required noninteracting matrix elements
  calcAllCoulombMatrixElements(vijck, vabck, psiai, evalai, psibe, potq,
                                poth, ist, par, planfw, planbw, fftwpsi);

  // Calculate w2h and w2e for all initial and final states by performing coherent summations
  calcCoherentSumsIntAR(w2h, w2e, Cbs, Ebs, vijck, vabck, par, ist);

  // Calculate the AR lifetimes and correctly scale based on number of stochastic orbitals used
  pf = fopen("samplingFinalStatesIntARLifetimes.dat", "w");
  for (iDeltaE = 0; iDeltaE < ist.nConsWindows; iDeltaE++) {
    for (iBiexc = 0; iBiexc < ist.nBiexcitons; iBiexc++) {
      iBiexcDeltaE = iBiexc + iDeltaE*ist.nBiexcitons;
      hRate[iBiexcDeltaE] = 0.0;
      eRate[iBiexcDeltaE] = 0.0;
      // Hole channel
      for (i = 0; i < nSI; i++) { 
        if (thetaIList[i + nSI*iBiexcDeltaE] != -1) {
          hRate[iBiexcDeltaE] += w2h[iBiexc*ist.itot*ist.nElecs + thetaIElecList[i + nSI*iBiexcDeltaE]*ist.itot + thetaIList[i + nSI*iBiexcDeltaE]];
        }
      }
      // Elec channel
      for (a = 0; a < nSA; a++) { 
        if (thetaAList[a + nSA*iBiexcDeltaE] != -1) {
          eRate[iBiexcDeltaE] += w2e[iBiexc*ist.atot*ist.nHoles + (thetaAList[a + nSA*iBiexcDeltaE]-ist.itot)*ist.nHoles + thetaAHoleList[a + nSA*iBiexcDeltaE]];
        }
      }
      // Scale the hRate and eRate and convert them to units of (ps)^-1
      hRate[iBiexcDeltaE] *= (nHotHnSIRatio[iBiexcDeltaE] * AUTOPS * TWOPI / (2.0 * deltaE[iDeltaE]));
      eRate[iBiexcDeltaE] *= (nHotEnSARatio[iBiexcDeltaE] * AUTOPS * TWOPI / (2.0 * deltaE[iDeltaE]));
      arRate[iBiexcDeltaE] = (hRate[iBiexcDeltaE] + eRate[iBiexcDeltaE]);
      // Print out AR lifetimes for each initial biexcitonic state for the current energy window
      fprintf(pf, "%ld %.6f %ld %.6f ", iDeltaE, deltaE[iDeltaE], iBiexc, biexcitonEnergies[iBiexc]);
      fprintf(pf, "%.8f %.8f ", hRate[iBiexcDeltaE], eRate[iBiexcDeltaE]);
      fprintf(pf, "%.8f %.8f\n", arRate[iBiexcDeltaE], 1.0/(arRate[iBiexcDeltaE]+EPS));
    }
    writeSeparation(pf);  
  }
  fclose(pf);

  // Calculate the Boltzmann weighted lifetimes for each energy conservation window
  pf = fopen("IntARLifetimesSamplingFinalStates.dat", "w");
  for (iDeltaE = 0; iDeltaE < ist.nConsWindows; iDeltaE++) {
    fprintf(pf, "Energy conservation window = %.6f\n\n", deltaE[iDeltaE]);
    calcBoltzmannWeightedRates(biexcitonEnergies, &(eRate[iDeltaE*ist.nBiexcitons]), 
            &(hRate[iDeltaE*ist.nBiexcitons]), par.temp, ist.nBiexcitons, pf); 
    writeSeparation(pf);
  }
  fclose(pf);

  // The next to lines are for convergence testing purposes (assumes nSI >= nSA)
  double sumDetHole, sumStoHole, holeChannelRate;
  double sumDetElec, sumStoElec, elecChannelRate;
  pf = fopen("ConvHoleElecSampFinalStatesIntARLifetimes.dat", "w");
  for (iDeltaE = 0; iDeltaE < ist.nConsWindows; iDeltaE++) {
    iBiexc = 0;
    for (iExc1 = 0; iExc1 < ist.nIntExcitons; iExc1++) {
      for (iExc2 = 0; iExc2 < ist.nIntExcitons; iExc2++) { 
        if ((Ebs[iExc1] + Ebs[iExc2]) < par.maxInitE+EPS) {
          holeChannelRate = 0.0;
          elecChannelRate = 0.0;
          iBiexcDeltaE = iBiexc + iDeltaE*ist.nBiexcitons;
          // Deterministic hole channel
          for (sumDetHole = 0.0, a = 0; a < ist.nElecs; a++) {
            for (tmp = 0.0, j = 0; j < ist.nHoles; j++) {
              tmp += Cbs[iExc1*ist.nNonIntExcitons + a*ist.nHoles + j];
            }
            sumDetHole += tmp*tmp*(double)(nHotHoleForElec[a + iBiexcDeltaE*ist.nElecs]);
          }
          // Deterministic elec channel
          for (sumDetElec = 0.0, i = 0; i < ist.nHoles; i++) {
            for (tmp = 0.0, b = 0; b < ist.nElecs; b++) {
              tmp += Cbs[iExc1*ist.nNonIntExcitons + b*ist.nHoles + i];
            }
            sumDetElec += tmp*tmp*(double)(nHotElecForHole[i + iBiexcDeltaE*ist.nHoles]);
          }
          sumStoHole = 0.0; sumStoElec = 0.0;
          for (i = 0; i < nSI; i++) {
            // Stochastic hole channel
            if (thetaIElecList[i + nSI*iBiexcDeltaE] != -1) {
              for (tmp = 0.0, j = 0; j < ist.nHoles; j++) {
                tmp += Cbs[iExc1*ist.nNonIntExcitons + thetaIElecList[i + nSI*iBiexcDeltaE]*ist.nHoles + j];
              }
              sumStoHole += tmp*tmp;
              holeChannelRate += w2h[iBiexc*ist.itot*ist.nElecs + thetaIElecList[i + nSI*iBiexcDeltaE]*ist.itot + thetaIList[i + nSI*iBiexcDeltaE]];
            }
            // Stochastic elec channel
            if (thetaAHoleList[i + nSA*iBiexcDeltaE] != -1) {
              for (tmp = 0.0, b = 0; b < ist.nElecs; b++) {
                tmp += Cbs[iExc1*ist.nNonIntExcitons + b*ist.nHoles + thetaAHoleList[i + nSA*iBiexcDeltaE]];
              }
              sumStoElec += tmp*tmp;
              elecChannelRate += w2e[iBiexc*ist.atot*ist.nHoles + (thetaAList[i + nSA*iBiexcDeltaE]-ist.itot)*ist.nHoles + thetaAHoleList[i + nSA*iBiexcDeltaE]];
            }
            // Scale and print the results every 5 sampling steps
            if (! (i % 5) && sumStoElec > EPS && sumStoHole > EPS) {
              fprintf(pf, "%ld %.12f  ", i, (holeChannelRate*(sumDetHole/sumStoHole)*AUTOPS*TWOPI/(2.0*deltaE[iDeltaE])));
              fprintf(pf, "%.12f\n", (elecChannelRate*(sumDetElec/sumStoElec)*AUTOPS*TWOPI/(2.0*deltaE[iDeltaE])));
            }
          }
          iBiexc++;
        }
      }
    }
  }
  fclose(pf);  

  // Write ending of function
  fprintf(stdout, "\nFinished the calculation of the interacting AR lifetimes by sampling the final states\n");
  writeCurrentTime(stdout);
  fflush(stdout);

  // Free dynamically allocated memory 
  free(w2h); free(w2e); free(biexcitonEnergies);
  free(eRate); free(hRate); free(arRate); 
  free(nHotHoleForElec); free(nHotElecForHole);
  free(elecHotHoleList); free(holeHotElecList);
  free(thetaIList); free(thetaAList);
  free(thetaIElecList); free(thetaAHoleList);
  free(nHotHnSIRatio); free(nHotEnSARatio);
  free(globalThetaIElecList); free(globalThetaAHoleList);
  free(globalThetaAList); free(globalThetaIList);

  return;
}

/*****************************************************************************/

void calcDeterministicIntAR(double *Cbs, double *Ebs, double *vijck, double *vabck, 
                     double *evalbe, double *evalai, lng_st ist, par_st par)
{
  FILE *pf;
  long a, i;
  long iDeltaE, iBiexc, iBiexcDeltaE;
  double *w2h, *w2e;
  double *eRate, *hRate, *arRate, *biexcitonEnergies;
  double deltaE[ist.nConsWindows]; // energy conservation windowns -> equally spaced from par.maxDeltaE to 0.1*par.maxDeltaE

  // Write beginning of function
  writeSeparation(stdout);
  writeCurrentTime(stdout);
  fprintf(stdout, "Beginnning the calculation of the deterministic interacting AR lifetimes\n\n");
  fprintf(stdout, "The number of energy conservation windows = %ld\n", ist.nConsWindows);
  fprintf(stdout, "The largest energy conservation windows = %.6f %.6f\n", par.maxDeltaE, par.maxDeltaE*AUTOEV);
  fprintf(stdout, "The maximum number of deterministic hot hole states = %ld\n", ist.itot);
  fprintf(stdout, "The maximum number of deterministic hot elec states = %ld\n", ist.atot);
  fflush(stdout);

  // Useful integers 
  const long THL = ist.nNonIntExcitons;
  const long TH2L = THL * ist.nHoles;
  const long THL2 = THL * ist.nElecs;

  // Dynamically allocate memory
  if ((biexcitonEnergies = (double *) calloc(ist.nBiexcitons, sizeof(double))) == NULL) nerror("biexcitonEnergies");
  if ((eRate = (double *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(double))) == NULL) nerror("eRate");
  if ((hRate = (double *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(double))) == NULL) nerror("hRate");
  if ((arRate = (double *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(double))) == NULL) nerror("arRate");
  if ((w2h = (double *) calloc(ist.nElecs*ist.itot*ist.nBiexcitons, sizeof(double))) == NULL) nerror("w2h");
  if ((w2e = (double *) calloc(ist.atot*ist.nHoles*ist.nBiexcitons, sizeof(double))) == NULL) nerror("w2e");

  // Create multiple energy conservation windows so one does not have to repeat the entire calculation
  fillEnergyConservationWindows(deltaE, par.maxDeltaE, ist.nConsWindows);

  // Fill in biexcitonEnergies
  fillIntBiexcEnergies(biexcitonEnergies, Ebs, ist, par);
  
  // Calculate w2h and w2e for all initial and final states
  calcCoherentSumsIntAR(w2h, w2e, Cbs, Ebs, vijck, vabck, par, ist);

  // Calculate the elec and hole channel rates
  pf = fopen("deterministicIntARLifetimes.dat", "w");
  for (iDeltaE = 0; iDeltaE < ist.nConsWindows; iDeltaE++) {
    for (iBiexc = 0; iBiexc < ist.nBiexcitons; iBiexc++) {
      iBiexcDeltaE = iBiexc + iDeltaE*ist.nBiexcitons;
      // hole channel
      for (a = 0; a < ist.nElecs; a++) {
        for (i = 0; i < ist.itot; i++) {
          // check if a,i (elec, hot hole) pair conserve energy
          if (fabs(evalbe[a+ist.nHoles] - evalai[i] - biexcitonEnergies[iBiexc]) <= deltaE[iDeltaE]) { 
            hRate[iBiexcDeltaE] += w2h[iBiexc*ist.itot*ist.nElecs + a*ist.itot + i];
          }
        }
      }
      // elec channel
      for (a = 0; a < ist.atot; a++) {
        for (i = 0; i < ist.nHoles; i++) {
          // check if a,i (hot elec, hole) pair conserve energy
          if (fabs(evalai[a+ist.itot] - evalbe[i] - biexcitonEnergies[iBiexc]) <= deltaE[iDeltaE]) { 
            eRate[iBiexcDeltaE] += w2e[iBiexc*ist.atot*ist.nHoles + a*ist.nHoles + i];
          }
        }
      }
      // scale the rates
      hRate[iBiexcDeltaE] *= (AUTOPS*TWOPI/(2.0*deltaE[iDeltaE]));
      eRate[iBiexcDeltaE] *= (AUTOPS*TWOPI/(2.0*deltaE[iDeltaE]));
      arRate[iBiexcDeltaE] = (hRate[iBiexcDeltaE] + eRate[iBiexcDeltaE]);
      // print out AR lifetimes for each initial biexcitonic state fir the current energy window
      fprintf(pf, "%ld %.6f %ld %.6f ", iDeltaE, deltaE[iDeltaE], iBiexc, biexcitonEnergies[iBiexc]);
      fprintf(pf, "%.8f %.8f ", hRate[iBiexcDeltaE], eRate[iBiexcDeltaE]);
      fprintf(pf, "%.8f %.8f\n", arRate[iBiexcDeltaE], 1.0/(arRate[iBiexcDeltaE]+EPS));
      fflush(pf);
    }
    writeSeparation(pf);  
  }
  fclose(pf);

  // Calculate the Boltzmann weighted lifetimes for each energy conservation window
  pf = fopen("IntARLifetimesDeterministic.dat", "w");
  for (iDeltaE = 0; iDeltaE < ist.nConsWindows; iDeltaE++) {
    fprintf(pf, "Energy conservation window = %.6f\n\n", deltaE[iDeltaE]);
    calcBoltzmannWeightedRates(biexcitonEnergies, &(eRate[iDeltaE*ist.nBiexcitons]), 
            &(hRate[iDeltaE*ist.nBiexcitons]), par.temp, ist.nBiexcitons, pf); 
    writeSeparation(pf);
  }
  fclose(pf);

  // Write ending of function
  fprintf(stdout, "\nFinished the calculation of the deterministic interacting AR lifetimes\n");
  writeCurrentTime(stdout);
  fflush(stdout);

  // Free dynamically allocated memory 
  free(w2h); free(w2e); 
  free(eRate); free(hRate); free(arRate); 
  free(biexcitonEnergies); 

  return;
}

/*****************************************************************************/

void calcCoherentSumsIntAR(double *w2h, double *w2e, double *Cbs, double *Ebs, 
                          double *vijck, double *vabck, par_st par, lng_st ist) 
{
  long a, i, b, j, c, k, iExc1, iExc2, iBiexc;  
  double sum;   

  // Useful integers 
  const long THL = ist.nNonIntExcitons;
  const long TH2L = THL * ist.nHoles;
  const long THL2 = THL * ist.nElecs;

  // Calculate w2h and w2e for all initial and final states
  iBiexc = 0;
  for (iExc1 = 0; iExc1 < ist.nIntExcitons; iExc1++) {
    for (iExc2 = 0; iExc2 < ist.nIntExcitons; iExc2++) { 
      if ((Ebs[iExc1] + Ebs[iExc2]) < par.maxInitE+EPS) {
        // Loop over final states for the hole channel
        for (a = 0; a < ist.nElecs; a++) { // spectator, delta_ab
          for (i = 0; i < ist.itot; i++) {
            // Coherent sum
            sum = 0.0;
            for (j = 0; j < ist.nHoles; j++) { 
              for (c = 0; c < ist.nElecs; c++) { 
                for (k = 0; k < ist.nHoles; k++) {
                	sum += vijck[i*TH2L + j*THL + c*ist.nHoles + k]              // Vijck
                          * Cbs[iExc1*ist.nNonIntExcitons + a*ist.nHoles + j]  // Caj
                          * Cbs[iExc2*ist.nNonIntExcitons + c*ist.nHoles + k]; // Cck
                }
              }
            }
            w2h[iBiexc*ist.nElecs*ist.itot + a*ist.itot + i] = sum*sum;
          }
        }
        // Loop over final states for the elec channel
        for (a = 0; a < ist.atot; a++) {
          for (i = 0; i < ist.nHoles; i++) { // spectator, delta_ij
            // Coherent sum
            sum = 0.0;
            for (b = 0; b < ist.nElecs; b++) {
              for (c = 0; c < ist.nElecs; c++) {
                for (k = 0; k < ist.nHoles; k++) {
                	sum += vabck[a*THL2 + b*THL + c*ist.nHoles + k]              // Vabck
                          * Cbs[iExc1*ist.nNonIntExcitons + b*ist.nHoles + i]  // Cbi
                          * Cbs[iExc2*ist.nNonIntExcitons + c*ist.nHoles + k]; // Cck
                }
              }
            }     
            w2e[iBiexc*ist.atot*ist.nHoles + a*ist.nHoles + i] = sum*sum;
          }
        }
        iBiexc++;
      }
    }
  }

  return;
}

/*****************************************************************************/

void calcStochasticNonIntAR(double *vijck, double *vabck, double *psibe, double *evalbe, double *psiai, double *evalai, 
                    zomplex *potq, lng_st ist, par_st par, fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwpsi) 
{
  FILE *pf, *fp; 
  long a, i, b, c, j, k, idum, tid;  
  long bIndex, jIndex, cIndex, kIndex;
  long aIndex, iIndex, jckIndex, bckIndex;
  long iBiexc, iDeltaE, iBiexcDeltaE, iBiexcDeltaEHotHole, iBiexcDeltaEHotElec;
  long numHotHoleEigenstates, numHotElecEigenstates;
  long *finalHoleList, *finalElecList, *nFinalHoles, *nFinalElecs;
  double tmp, *eRate, *hRate, *arRate;
  long nSA, nSI, *thetaAList, *thetaIList;
  long randomIndex, nMatrixElementsCalculated, coulombMatrixElementIndex;
  long *prevCalcHoleList, *prevCalcElecList;
  double deltaE[ist.nConsWindows]; // energy conservation windowns -> equally spaced from par.maxDeltaE to 0.1*par.maxDeltaE
  double nHotEnSARatio, nHotHnSIRatio, *biexcitonEnergies;
  double *psiB, *psiJ, *psiC, *psiK, *psiA, *psiI; // pointers only, will not allocate memory to them
  nonIntBiexc *initBiexcStates;

  // Useful integers 
  nSA = ist.nStochHotElecs; // number of stochastic hot elec orbitals 
  nSI = ist.nStochHotHoles; // number of stochastic hot hole orbitals
  nMatrixElementsCalculated = (nSA+nSI)*ist.nConsWindows*ist.nBiexcitons;
  const long THL = ist.nNonIntExcitons;
  const long TH2L = THL * ist.nHoles;
  const long THL2 = THL * ist.nElecs;

  // Write beginning of function
  writeSeparation(stdout);
  writeCurrentTime(stdout);
  fprintf(stdout, "Beginnning the calculation of the stochastic noninteracting AR lifetimes\n\n");
  fprintf(stdout, "The number of energy conservation windows = %ld\n", ist.nConsWindows);
  fprintf(stdout, "The largest energy conservation windows = %.6f %.6f\n", par.maxDeltaE, par.maxDeltaE*AUTOEV);
  fprintf(stdout, "The maximum number of deterministic hot hole states = %ld\n", ist.itot);
  fprintf(stdout, "The maximum number of deterministic hot elec states = %ld\n", ist.atot);
  fprintf(stdout, "The number of stochastic hot hole states sampled = %ld\n", nSI);
  fprintf(stdout, "The number of stochastic hot elec states sampled = %ld\n", nSA);
  fflush(stdout);

  // Dynamically allocate memory
  if ((nFinalHoles = (long *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(long))) == NULL) nerror("nFinalHoles");
  if ((nFinalElecs = (long *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(long))) == NULL) nerror("nFinalElecs");
  if ((thetaAList = (long *) calloc(nSA*ist.nBiexcitons*ist.nConsWindows, sizeof(long))) == NULL) nerror("thetaAList");
  if ((thetaIList = (long *) calloc(nSI*ist.nBiexcitons*ist.nConsWindows, sizeof(long))) == NULL) nerror("thetaIList");
  if ((prevCalcHoleList = (long *) calloc(ist.itot*TH2L+1, sizeof(long))) == NULL) nerror("prevCalcHoleList");
  if ((prevCalcElecList = (long *) calloc(ist.atot*THL2+1, sizeof(long))) == NULL) nerror("prevCalcElecList");
  if ((eRate = (double *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(double))) == NULL) nerror("eRate");
  if ((hRate = (double *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(double))) == NULL) nerror("hRate");
  if ((arRate = (double *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(double))) == NULL) nerror("arRate");
  if ((biexcitonEnergies = (double *) calloc(ist.nBiexcitons, sizeof(double))) == NULL) nerror("biexcitonEnergies");
  if ((initBiexcStates = (nonIntBiexc *) calloc(ist.nBiexcitons, sizeof(nonIntBiexc))) == NULL) nerror("initBiexcStates");
  if ((finalHoleList = (long *) calloc(ist.nBiexcitons*ist.nConsWindows*ist.itot+1, sizeof(long))) == NULL) nerror("finalHoleList");
  if ((finalElecList = (long *) calloc(ist.nBiexcitons*ist.nConsWindows*ist.atot+1, sizeof(long))) == NULL) nerror("finalElecList");
  
  // Create multiple energy conservation windows so one does not have to repeat the entire calculation
  fillEnergyConservationWindows(deltaE, par.maxDeltaE, ist.nConsWindows);

  // Fill in the initial states into the biexc structure - loop over all possible (b,j,k,c) pairs
  iBiexc = fillNonIntBiexcStruct(initBiexcStates, evalbe, par.maxInitE, ist);
  for (iBiexc = 0; iBiexc < ist.nBiexcitons; iBiexc++) {
    biexcitonEnergies[iBiexc] = initBiexcStates[iBiexc].energy;
  }
  fprintf(stdout, "The number of initial biexcitonic states = %d\n", iBiexc); // ist.nBiexcitons should also have same information
  fflush(stdout);

  // Create lists and counts of the allowed final states for each of the possible initial state and energyConsWindow
  storeAllowedFinalStatesNonIntAR(finalHoleList, nFinalHoles, finalElecList, nFinalElecs, evalai, ist.itot, 
              ist.atot, initBiexcStates, ist.nBiexcitons, deltaE, ist.nConsWindows);

  // Generate the stochastic orbitals (sampling the hot carriers) -> just need a list of the indices
  idum = ist.seed;  
  pf = fopen("sHotElecStates.dat", "w");
  fp = fopen("sHotHoleStates.dat", "w");
  for (iDeltaE = 0; iDeltaE < ist.nConsWindows; iDeltaE++) {
    for (iBiexc = 0; iBiexc < ist.nBiexcitons; iBiexc++) {
      iBiexcDeltaE = iBiexc + iDeltaE*ist.nBiexcitons;
      for (i = 0; i < nSI; i++) {
        randomIndex = (long) ((double)(nFinalHoles[iBiexcDeltaE])*ran_nrc(&idum)); // ranges from 0 to nFinalHoles[iBiexcDeltaE]
        thetaIList[i + nSI*iBiexcDeltaE] = finalHoleList[randomIndex + iBiexcDeltaE*ist.itot];
        fprintf(fp, "%ld %ld %ld\n", i, randomIndex, thetaIList[i + nSI*iBiexcDeltaE]);
      }
      for (a = 0; a < nSA; a++) {
        randomIndex = (long) ((double)(nFinalElecs[iBiexcDeltaE])*ran_nrc(&idum)); // ranges from 0 to nFinalElecs[iBiexcDeltaE] 
        thetaAList[a + nSA*iBiexcDeltaE] = finalElecList[randomIndex + iBiexcDeltaE*ist.atot];
        fprintf(pf, "%ld %ld %ld\n", a, randomIndex, thetaAList[a + nSA*iBiexcDeltaE]);
      }
    }
  }
  fclose(fp);
  fclose(pf);

  // Use biexciton structure to create lists of the Coulomb matrix elements that must be calculated
  //nMatrixElementsCalculated = determineHoleMatrixElementsToCalc(initBiexcStates, )  

  // Calculating all of the required stochastic matrix elements
  //calcAllStochasticCoulombMatrixElements(vijck, vabck, );

  // Calculate AR lifetimes for the hole (1/hRate) and electron (1/eRate) channels
  pf = fopen("stochasticNonIntARLifetimes.dat", "w");
  for (iDeltaE = 0; iDeltaE < ist.nConsWindows; iDeltaE++) {
    for (iBiexc = 0; iBiexc < ist.nBiexcitons; iBiexc++) {
      iBiexcDeltaE = iBiexc + iDeltaE*ist.nBiexcitons;
      bIndex = initBiexcStates[iBiexc].b - ist.lumoIndex;
      psiB = &(psibe[initBiexcStates[iBiexc].b*ist.nGridPoints]);
      jIndex = initBiexcStates[iBiexc].j;
      psiJ = &(psibe[initBiexcStates[iBiexc].j*ist.nGridPoints]);
      cIndex = initBiexcStates[iBiexc].c - ist.lumoIndex;
      psiC = &(psibe[initBiexcStates[iBiexc].c*ist.nGridPoints]);
      kIndex = initBiexcStates[iBiexc].k;
      psiK = &(psibe[initBiexcStates[iBiexc].k*ist.nGridPoints]);
      jckIndex = jIndex*THL + cIndex*ist.nHoles + kIndex;
      bckIndex = bIndex*THL + cIndex*ist.nHoles + kIndex;
      // Calculate the required hot hole Coulomb matrix elements and resulting hRate
    //#pragma omp parallel for private(i, iIndex, coulombMatrixElementIndex, psiI, tid)
      for (i = 0; i < nSI; i++) {
        iIndex = thetaIList[i + nSI*iBiexcDeltaE];
        coulombMatrixElementIndex = iIndex*TH2L + jckIndex;
        if (prevCalcHoleList[coulombMatrixElementIndex] == 0) {
          psiI = &(psiai[iIndex*ist.nGridPoints]);
          //tid = omp_get_thread_num();
          tid = 0;
          vijck[coulombMatrixElementIndex] = calcOneCoulombMatrixElement(psiC, psiK, psiJ, psiI, potq, 
                                              par.dv, ist.nGridPoints, planfw[tid], planbw[tid], &fftwpsi[tid*ist.nGridPoints]);
          prevCalcHoleList[coulombMatrixElementIndex] = 1; // to prevent same matrix element being calculated again
        }
        else {
          nMatrixElementsCalculated--;
        }
      //#pragma omp critical
        hRate[iBiexcDeltaE] += vijck[coulombMatrixElementIndex]*vijck[coulombMatrixElementIndex];
      }
      // Calculate the required hot elec Coulomb matrix elements and resulting eRate
    //#pragma omp parallel for private(a, aIndex, coulombMatrixElementIndex, psiA, tid)
      for (a = 0; a < nSA; a++) {
        aIndex = thetaAList[a + nSA*iBiexcDeltaE];
        coulombMatrixElementIndex = (aIndex-ist.itot)*THL2 + bckIndex;
        if (prevCalcElecList[coulombMatrixElementIndex] == 0) {
          psiA = &(psiai[aIndex*ist.nGridPoints]);
          //tid = omp_get_thread_num();
          tid = 0;
          vabck[coulombMatrixElementIndex] = calcOneCoulombMatrixElement(psiC, psiK, psiB, psiA, potq,
                                              par.dv, ist.nGridPoints, planfw[tid], planbw[tid], &fftwpsi[tid*ist.nGridPoints]);
          prevCalcElecList[coulombMatrixElementIndex] = 1; // to prevent same matrix element being calculated again
        }
        else {
          nMatrixElementsCalculated--;
        }
      //#pragma omp critical
        eRate[iBiexcDeltaE] += vabck[coulombMatrixElementIndex]*vabck[coulombMatrixElementIndex];
      }
      // Scale the hRate and eRate and convert them to units of (ps)^-1
      nHotHnSIRatio = ((double)(nFinalHoles[iBiexcDeltaE]) / (double)(nSI));
      nHotEnSARatio = ((double)(nFinalElecs[iBiexcDeltaE]) / (double)(nSA));
      hRate[iBiexcDeltaE] *= (nHotHnSIRatio * AUTOPS * TWOPI / (2.0 * deltaE[iDeltaE]));
      eRate[iBiexcDeltaE] *= (nHotEnSARatio * AUTOPS * TWOPI / (2.0 * deltaE[iDeltaE]));
      arRate[iBiexcDeltaE] = (hRate[iBiexcDeltaE] + eRate[iBiexcDeltaE]);
      fprintf(pf, "%ld %.6f %ld %.6f ", iDeltaE, deltaE[iDeltaE], initBiexcStates[iBiexc].index, initBiexcStates[iBiexc].energy);
      fprintf(pf, "%ld %ld %ld %ld %ld %ld ", bIndex+ist.lumoIndex, jIndex, cIndex+ist.lumoIndex, kIndex, nFinalHoles[iBiexcDeltaE], nFinalElecs[iBiexcDeltaE]);
      fprintf(pf, "%.8f %.8f %.8f ", hRate[iBiexcDeltaE], eRate[iBiexcDeltaE], arRate[iBiexcDeltaE]);
      fprintf(pf, "%.8f\n", 1.0/(arRate[iBiexcDeltaE]+EPS));  
    }
    writeSeparation(pf);
  }
  fclose(pf);
  fprintf(stdout, "The number of Coulomb matrix elements that were calculated = %ld\n", nMatrixElementsCalculated);
  fflush(stdout);

  // Calculate the Boltzmann weighted lifetimes for each energy conservation window
  pf = fopen("NonIntARLifetimesStochastic.dat", "w");
  for (iDeltaE = 0; iDeltaE < ist.nConsWindows; iDeltaE++) {
    fprintf(pf, "Energy conservation window = %.6f\n\n", deltaE[iDeltaE]);
    calcBoltzmannWeightedRates(biexcitonEnergies, &(eRate[iDeltaE*ist.nBiexcitons]), &(hRate[iDeltaE*ist.nBiexcitons]), par.temp, ist.nBiexcitons, pf);
    writeSeparation(pf);
  }
  fclose(pf);

  // The next to lines are for convergence testing purposes (assumes nSI >= nSA)
  pf = fopen("ConvHoleElecNonIntARLifetimes.dat", "w");
  for (iDeltaE = 0; iDeltaE < ist.nConsWindows; iDeltaE++) {
    for (iBiexc = 0; iBiexc < ist.nBiexcitons; iBiexc++) {
      iBiexcDeltaE = iBiexc + iDeltaE*ist.nBiexcitons;
      bIndex = initBiexcStates[iBiexc].b - ist.lumoIndex;
      jIndex = initBiexcStates[iBiexc].j;
      cIndex = initBiexcStates[iBiexc].c - ist.lumoIndex;
      kIndex = initBiexcStates[iBiexc].k;
      jckIndex = jIndex*THL + cIndex*ist.nHoles + kIndex;
      bckIndex = bIndex*THL + cIndex*ist.nHoles + kIndex;
      hRate[iBiexcDeltaE] = 0.0;
      eRate[iBiexcDeltaE] = 0.0;
      for (i = 0; i < nSI; i++) {
        a = i;
        iIndex = thetaIList[i + nSI*iBiexcDeltaE];
        coulombMatrixElementIndex = iIndex*TH2L + jckIndex;
        hRate[iBiexcDeltaE] += vijck[coulombMatrixElementIndex]*vijck[coulombMatrixElementIndex];
        nHotHnSIRatio = ((double)(nFinalHoles[iBiexcDeltaE]) / (double)(i+1));
        fprintf(pf, "%ld %ld %ld %.12f  ", i, iIndex, i+nSI*iBiexcDeltaE, (hRate[iBiexcDeltaE]*nHotHnSIRatio*AUTOPS*TWOPI/(2.0*deltaE[iDeltaE])));
        aIndex = thetaAList[a + nSA*iBiexcDeltaE];
        coulombMatrixElementIndex = (aIndex-ist.itot)*THL2 + bckIndex;
        eRate[iBiexcDeltaE] += vabck[coulombMatrixElementIndex]*vabck[coulombMatrixElementIndex];
        nHotEnSARatio = ((double)(nFinalElecs[iBiexcDeltaE]) / (double)(i+1));
        fprintf(pf, "%ld %ld %.12f\n", i, aIndex, (eRate[iBiexcDeltaE]*nHotEnSARatio*AUTOPS*TWOPI/(2.0*deltaE[iDeltaE])));
      }
      // Scale the hRate and eRate and convert them to units of (ps)^-1
      nHotHnSIRatio = ((double)(nFinalHoles[iBiexcDeltaE]) / (double)(nSI));
      nHotEnSARatio = ((double)(nFinalElecs[iBiexcDeltaE]) / (double)(nSA));
      hRate[iBiexcDeltaE] *= (nHotHnSIRatio * AUTOPS * TWOPI / (2.0 * deltaE[iDeltaE]));
      eRate[iBiexcDeltaE] *= (nHotEnSARatio * AUTOPS * TWOPI / (2.0 * deltaE[iDeltaE]));
      arRate[iBiexcDeltaE] = (hRate[iBiexcDeltaE] + eRate[iBiexcDeltaE]);
    }
  }
  fclose(pf);  

  // Write ending of function
  fprintf(stdout, "\nFinished the calculation of the stochastic noninteracting AR lifetimes\n");
  writeCurrentTime(stdout);
  fflush(stdout);

  // Free dynamically allocated memory
  free(biexcitonEnergies);
  free(thetaAList); free(thetaIList);
  free(initBiexcStates);
  free(finalHoleList); free(finalElecList); free(nFinalHoles); free(nFinalElecs);
  free(eRate); free(hRate); free(arRate);

  return;
}

/*****************************************************************************/
//

void calcDeterministicNonIntAR(double *vijck, double *vabck, double *evalbe, double *evalai, lng_st ist, par_st par) {
  FILE *pf; 
  long a, i, b, c, j, k, idum;
  long bIndex, jIndex, cIndex, kIndex, aIndex, iIndex;
  long iBiexc, iDeltaE, iBiexcDeltaE, iBiexcDeltaEHotHole, iBiexcDeltaEHotElec;
  long numHotHoleEigenstates, numHotElecEigenstates;
  long *finalHoleList, *finalElecList, *nFinalHoles, *nFinalElecs;
  double tmp, *eRate, *hRate, *arRate, *biexcitonEnergies;
  double deltaE[ist.nConsWindows]; // energy conservation windowns -> equally spaced from par.maxDeltaE to 0.1*par.maxDeltaE
  nonIntBiexc *initBiexcStates;

  // Useful integers 
  const long THL = ist.nNonIntExcitons;
  const long TH2L = THL * ist.nHoles;
  const long THL2 = THL * ist.nElecs;

  // Write beginning of function
  writeSeparation(stdout);
  writeCurrentTime(stdout);
  fprintf(stdout, "Beginnning the calculation of the deterministic noninteracting AR lifetimes\n\n");
  fprintf(stdout, "The number of energy conservation windows = %d\n", ist.nConsWindows);
  fprintf(stdout, "The largest energy conservation windows = %.6f %.6f\n", par.maxDeltaE, par.maxDeltaE*AUTOEV);
  fprintf(stdout, "The maximum number of deterministic hot hole states = %d\n", ist.itot);
  fprintf(stdout, "The maximum number of deterministic hot elec states = %d\n", ist.atot);
  fflush(stdout);

  // Dynamically allocate memory
  if ((nFinalHoles = (long *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(long))) == NULL) nerror("nFinalHoles");
  if ((nFinalElecs = (long *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(long))) == NULL) nerror("nFinalElecs");
  if ((eRate = (double *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(double))) == NULL) nerror("eRate");
  if ((hRate = (double *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(double))) == NULL) nerror("hRate");
  if ((arRate = (double *) calloc(ist.nBiexcitons*ist.nConsWindows, sizeof(double))) == NULL) nerror("arRate");
  if ((biexcitonEnergies = (double *) calloc(ist.nBiexcitons, sizeof(double))) == NULL) nerror("biexcitonEnergies");
  if ((initBiexcStates = (nonIntBiexc *) calloc(ist.nBiexcitons, sizeof(nonIntBiexc))) == NULL) nerror("initBiexcStates");
  if ((finalHoleList = (long *) calloc(ist.nBiexcitons*ist.nConsWindows*ist.itot+1, sizeof(long))) == NULL) nerror("finalHoleList");
  if ((finalElecList = (long *) calloc(ist.nBiexcitons*ist.nConsWindows*ist.atot+1, sizeof(long))) == NULL) nerror("finalElecList");
  
  // Create multiple energy conservation windows so one does not have to repeat the entire calculation
  fillEnergyConservationWindows(deltaE, par.maxDeltaE, ist.nConsWindows);

  // Fill in the initial states into the biexc structure - loop over all possible (b,j,k,c) pairs
  iBiexc = fillNonIntBiexcStruct(initBiexcStates, evalbe, par.maxInitE, ist);
  fprintf(stdout, "The number of initial biexcitonic states = %d\n", iBiexc); // ist.nBiexcitons should also have same information
  fflush(stdout);

  // Create lists and counts of the allowed final states for each of the possible initial state and energyConsWindow
  storeAllowedFinalStatesNonIntAR(finalHoleList, nFinalHoles, finalElecList, nFinalElecs, evalai, ist.itot, 
              ist.atot, initBiexcStates, ist.nBiexcitons, deltaE, ist.nConsWindows);

  // Calculate AR lifetimes for the hole (1/hRate) and electron (1/eRate) channels
  pf = fopen("deterministicNonIntARLifetimes.dat", "w");
  for (iDeltaE = 0; iDeltaE < ist.nConsWindows; iDeltaE++) {
    for (iBiexc = 0; iBiexc < ist.nBiexcitons; iBiexc++) {
      iBiexcDeltaE = iBiexc + iDeltaE*ist.nBiexcitons;
      iBiexcDeltaEHotHole = iBiexcDeltaE*ist.itot;
      iBiexcDeltaEHotElec = iBiexcDeltaE*ist.atot;
      bIndex = initBiexcStates[iBiexc].b - ist.lumoIndex;
      jIndex = initBiexcStates[iBiexc].j;
      cIndex = initBiexcStates[iBiexc].c - ist.lumoIndex;
      kIndex = initBiexcStates[iBiexc].k;
      for (i = 0; i < nFinalHoles[iBiexcDeltaE]; i++) {
        iIndex = finalHoleList[i + iBiexcDeltaEHotHole];
        tmp = vijck[iIndex*TH2L + jIndex*THL + cIndex*ist.nHoles + kIndex];
        hRate[iBiexcDeltaE] += tmp*tmp;
      }
      for (a = 0; a < nFinalElecs[iBiexcDeltaE]; a++) {
        aIndex = finalElecList[a + iBiexcDeltaEHotElec];
        tmp = vabck[(aIndex-ist.itot)*THL2 + bIndex*THL + cIndex*ist.nHoles + kIndex];
        eRate[iBiexcDeltaE] += tmp*tmp;
      }
      hRate[iBiexcDeltaE] *= (AUTOPS*TWOPI/(2.0*deltaE[iDeltaE]));
      eRate[iBiexcDeltaE] *= (AUTOPS*TWOPI/(2.0*deltaE[iDeltaE]));
      arRate[iBiexcDeltaE] = (hRate[iBiexcDeltaE] + eRate[iBiexcDeltaE]);
      fprintf(pf, "%ld %.6f %ld %.6f ", iDeltaE, deltaE[iDeltaE], initBiexcStates[iBiexc].index, initBiexcStates[iBiexc].energy);
      fprintf(pf, "%ld %ld %ld %ld %ld %ld ", bIndex+ist.lumoIndex, jIndex, cIndex+ist.lumoIndex, kIndex, nFinalHoles[iBiexcDeltaE], nFinalElecs[iBiexcDeltaE]);
      fprintf(pf, "%.8f %.8f %.8f ", hRate[iBiexcDeltaE], eRate[iBiexcDeltaE], arRate[iBiexcDeltaE]);
      fprintf(pf, "%.8f\n", 1.0/(arRate[iBiexcDeltaE]+EPS));  
    }
    writeSeparation(pf);
  }
  fclose(pf);

  // Calculate the Boltzmann weighted lifetimes for each energy conservation window
  pf = fopen("NonIntARLifetimesDeterministic.dat", "w");
  for (iDeltaE = 0; iDeltaE < ist.nConsWindows; iDeltaE++) {
    fprintf(pf, "Energy conservation window = %.6f\n\n", deltaE[iDeltaE]);
    calcBoltzmannWeightedRates(biexcitonEnergies, &(eRate[iDeltaE*ist.nBiexcitons]), &(hRate[iDeltaE*ist.nBiexcitons]), par.temp, ist.nBiexcitons, pf);
    writeSeparation(pf);
  }
  fclose(pf);

  // Write ending of function
  fprintf(stdout, "\nFinished the calculation of the deterministic noninteracting AR lifetimes\n");
  writeCurrentTime(stdout);
  fflush(stdout);

  // Free dynamically allocated memory
  free(biexcitonEnergies); free(initBiexcStates);
  free(finalHoleList); free(finalElecList); free(nFinalHoles); free(nFinalElecs);
  free(eRate); free(hRate); free(arRate);

  return;
}

/*****************************************************************************/
// calculates and prints the boltzmann weighted AR, k_x-, k_x+ rates and 
// the corresponding lifetimes and boltzmann populations

void calcBoltzmannWeightedRates(double *energies, double *eRate, double *hRate, double temp, long nStates, FILE *pf) {
  FILE *fp;
  long i;
  double pF, bwAugerRate, bwElecRate, bwHoleRate;
  double stateProb, e0 = energies[0];  // will be energy of the lowest state
  double beta = AUTOEV/(KB*(temp+EPS));  

  bwAugerRate = bwElecRate = bwHoleRate = 0.0;
  fp = fopen("boltzmannStats.dat", "w");

  pF = calcPartitionFunction(energies, temp, nStates);
  for (i = 0; i < nStates; i++) if (energies[i] < e0) e0 = energies[i]; // in case energies isn't ordered
  for (i = 0; i < nStates; i++) {
    stateProb = exp(-beta*(energies[i]-e0));
    bwElecRate += eRate[i]*stateProb;
    bwHoleRate += hRate[i]*stateProb;
    fprintf(fp, "%.10f %.8f %.12f %.12f %.12f\n", energies[i], stateProb/pF, 
      1.0/(eRate[i]+EPS), 1.0/(hRate[i]+EPS), 1.0/(eRate[i]+hRate[i]+EPS));
  }
  fclose(fp);

  bwElecRate /= pF;
  bwHoleRate /= pF;
  bwAugerRate = bwElecRate + bwHoleRate;

  fprintf(pf, "Lowest energy biexciton = %.10f\n", e0);
  fprintf(pf, "Partition function = %.4f\n", pF);
  fprintf(pf, "Boltzmann weighted k_x- = %.6f 1/ps\n", bwElecRate);
  fprintf(pf, "Boltzmann weighted k_x+ = %.6f 1/ps\n", bwHoleRate);
  fprintf(pf, "Boltzmann weighted kAR  = %.6f 1/ps\n", bwAugerRate);
  fprintf(pf, "Boltzmann weighted x- lifetime = %.8f ps\n", 1.0/(bwElecRate+EPS));
  fprintf(pf, "Boltzmann weighted x+ lifetime = %.8f ps\n", 1.0/(bwHoleRate+EPS));
  fprintf(pf, "Boltzmann weighted Auger lifetime = %.8f ps\n", 1.0/(bwAugerRate+EPS));

  return;
}

/*****************************************************************************/
// returns the partition function given an array of energies, a temperature,
// and the number of states (i.e., length of the energies array)

double calcPartitionFunction(double *energies, double temp, long nStates) {
  long i;
  double pF = 0.0;
  double e0 = energies[0];  // will be energy of the lowest state
  double beta = AUTOEV/(KB*(temp+EPS)); 

  for (i = 0; i < nStates; i++) if (energies[i] < e0) e0 = energies[i]; // in case energies isn't ordered
  for (i = 0; i < nStates; i++) pF += exp(-beta*(energies[i]-e0));

  return pF;
}

/*****************************************************************************/
//

void fillEnergyConservationWindows(double *deltaE, double maxDeltaE, long nConsWindows) {
  long iDeltaE;

  for (iDeltaE = 0; iDeltaE < nConsWindows; iDeltaE++) {
    deltaE[iDeltaE] = (maxDeltaE*(1.0 - 0.1*(double)(iDeltaE)));
  }  

  return;
}

/*****************************************************************************/
