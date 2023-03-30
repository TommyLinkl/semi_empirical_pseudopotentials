/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "ar.h"

/****************************************************************************/
// Returns the number of initial interacting biexcitonic states 
// that have an energy lower than maxBiexcEnergy

long calcNumIntBiexcStates(double *Ebs, double maxBiexcEnergy, long nIntExcitons) {
	long nBiexcitons = 0; // keeps track of the number of allowed initial states
	long iIntExc1, iIntExc2; // indices of the 1st and 2nd initial exciton
	double biexcEnergy; // test biexciton energy to see if it is within allowed energy range

	for (iIntExc1 = 0; iIntExc1 < nIntExcitons; iIntExc1++) {
		for (iIntExc2 = 0; iIntExc2 < nIntExcitons; iIntExc2++) {
			biexcEnergy = (Ebs[iIntExc1] + Ebs[iIntExc2]);
			if (biexcEnergy < (maxBiexcEnergy+EPS)) {
				nBiexcitons++;
			}
		}
	}

	// Print number of initial biexcitonic states
	fprintf(stdout, "The number of initial interacting biexcitonic states = %ld\n", nBiexcitons);
	fflush(stdout);

	return nBiexcitons;
}

/****************************************************************************/
// Returns the number of initial noninteracting biexcitonic states 
// that have an energy lower than maxBiexcEnergy

long calcNumNonIntBiexcStates(double *evalbe, double maxBiexcEnergy, long nHoles, long nElecs) {
	long nBiexcitons = 0; // keeps track of the number of allowed initial states
	long b, j; // indices of the 1st initial exciton
	long c, k; // indices of the 2nd initial exciton
	double biexcEnergy; // test biexciton energy to see if it is within allowed energy range

	for (b = 0; b < nElecs; b++) {
		for (j = 0; j < nHoles; j++) {
			for (c = 0; c < nElecs; c++) {
				for (k = 0; k < nHoles; k++) {
					biexcEnergy = (evalbe[nHoles+b] + evalbe[nHoles+c] - evalbe[j] - evalbe[k]);
					if (biexcEnergy < (maxBiexcEnergy+EPS)) {
						nBiexcitons++;
					}
				}
			}
		}
	}

	// Print number of initial biexcitonic states
	fprintf(stdout, "The number of initial noninteracting biexcitonic states = %ld\n", nBiexcitons);
	fflush(stdout);

	return nBiexcitons;
}

/****************************************************************************/
// Fills the noninteracting biexictonic state structure and returns the number 
// that were stored. Assumes all hole eigenstates from eval[0]-eval[nHoles]
// are eigenstates and the electrons are listed from eval[lumoIndex]-eval[nElecs]

long fillNonIntBiexcStruct(nonIntBiexc *biExciton, double *eval, double maxEnergy, lng_st ist) {
 	FILE *pf;
 	long b, j, k, c, nBiexcitons = 0;

	pf = fopen("biexcitonStates.dat", "w");
	for (b = 0; b < ist.nElecs; b++) { 
		for (j = 0; j < ist.nHoles; j++) {
			for (c = 0; c < ist.nElecs; c++) {
				for (k = 0; k < ist.nHoles; k++) {
					biExciton[nBiexcitons].energy = eval[ist.lumoIndex + b] + eval[ist.lumoIndex + c] - eval[j] - eval[k];
			  		if (biExciton[nBiexcitons].energy < (maxEnergy+EPS)) { // checks to see if initial biexciton is in allowed energy range
						biExciton[nBiexcitons].index = nBiexcitons;
						biExciton[nBiexcitons].b = b + ist.lumoIndex; // elec of exciton 1
						biExciton[nBiexcitons].bEnergy = eval[ist.lumoIndex + b]; // energy of e1
						biExciton[nBiexcitons].j = j; // hole of exciton 1
						biExciton[nBiexcitons].jEnergy = -eval[j]; // energy of h1 
						biExciton[nBiexcitons].c = c + ist.lumoIndex; // elec of exciton 2
						biExciton[nBiexcitons].cEnergy = eval[ist.lumoIndex + c]; // energy of e2
						biExciton[nBiexcitons].k = k; // hole of exciton 2
						biExciton[nBiexcitons].kEnergy = -eval[k]; // energy of h2
						writeNonIntBiexciton(biExciton[nBiexcitons], pf);
						nBiexcitons++; // increments biexciton index energy range of biexciton was satisfied
			  		}
				}
			}
		}
	}
	fclose(pf);

	return nBiexcitons;
}

/****************************************************************************/
// Fills the noninteracting biexictonic state structure and returns the number 
// that were stored. Assumes all hole eigenstates from eval[0]-eval[nHoles]
// are eigenstates and the electrons are listed from eval[lumoIndex]-eval[nElecs]

long fillIntBiexcEnergies(double *biexcitonEnergies, double *Ebs, lng_st ist, par_st par) {
 	long iExc1, iExc2, nBiexcitons = 0;
 	double biexcitonEnergy;

    for (iExc1 = 0; iExc1 < ist.nIntExcitons; iExc1++) {
    	for (iExc2 = 0; iExc2 < ist.nIntExcitons; iExc2++) { 
    		biexcitonEnergy = Ebs[iExc1] + Ebs[iExc2];
			if (biexcitonEnergy < par.maxInitE+EPS) {
				biexcitonEnergies[nBiexcitons] = biexcitonEnergy;
				nBiexcitons++;
			}
		}
	}

	return nBiexcitons;
}

/****************************************************************************/
// Create lists and counts of the allowed final states for each 
// of the possible initial state and energyConsWindow

void storeAllowedFinalStatesIntAR(long *elecHotHoleList, long *nHotHoleForElec, long *holeHotElecList, long *nHotElecForHole, 
				double *biexcitonEnergies, double *Ebs, double *evalbe, double *evalai, double *deltaE, lng_st ist, par_st par) {
	FILE *pf;
	long a, i, iExc1, iExc2, iDeltaE, iBiexc, iBiexcDeltaE;
	double biexcitonEnergy;

	pf = fopen("finalExcStatesIntAR.dat", "w");
	for (iDeltaE = 0; iDeltaE < ist.nConsWindows; iDeltaE++) {
		if (iDeltaE) writeSeparation(pf);
		fprintf(pf, "DeltaE index = %d has window = % .6f\n", iDeltaE, deltaE[iDeltaE]);
	    iBiexc = 0;
	    for (iExc1 = 0; iExc1 < ist.nIntExcitons; iExc1++) {
	    	for (iExc2 = 0; iExc2 < ist.nIntExcitons; iExc2++) { 
	    		biexcitonEnergy = Ebs[iExc1] + Ebs[iExc2];
				if (biexcitonEnergy < par.maxInitE+EPS) {
					biexcitonEnergies[iBiexc] = biexcitonEnergy;
	        		iBiexcDeltaE = iBiexc + iDeltaE*ist.nBiexcitons;
	        		// Loop over final states for the hole channel
	        		for (a = 0; a < ist.nElecs; a++) { // spectator, delta_ab
	        			for (i = 0; i < ist.itot; i++) {
	            			if (fabs(evalbe[a+ist.nHoles] - evalai[i] - biexcitonEnergies[iBiexc]) <= deltaE[iDeltaE]) { 
	            				elecHotHoleList[nHotHoleForElec[a + iBiexcDeltaE*ist.nElecs] + (a + iBiexcDeltaE*ist.nElecs)*ist.itot] = i;
	            				nHotHoleForElec[a + iBiexcDeltaE*ist.nElecs]++;
	            			}
	        			}
	        		}
	        		// Loop over final states for the elec channel
	          		for (a = 0; a < ist.atot; a++) {
	            		for (i = 0; i < ist.nHoles; i++) { // spectator, delta_ij
	              			if (fabs(evalai[a+ist.itot] - evalbe[i] - biexcitonEnergies[iBiexc]) <= deltaE[iDeltaE]) {
	                			holeHotElecList[nHotElecForHole[i + iBiexcDeltaE*ist.nHoles] + (i + iBiexcDeltaE*ist.nHoles)*ist.atot] = a+ist.itot;
	                			nHotElecForHole[i + iBiexcDeltaE*ist.nHoles]++;
	            			}
	            		} 
	          		}
					fprintf(pf, "\nBiexciton index = %d has energy = %.6f\n", iBiexc, biexcitonEnergies[iBiexc]);
					fprintf(pf, "The number of hot holes satisfying energy conservation = %d\n", nHotHoleForElec[iBiexcDeltaE*ist.nElecs]);
					fprintf(pf, "The number of hot elecs satisfying energy conservation = %d\n", nHotElecForHole[iBiexcDeltaE*ist.nHoles]);
	          		iBiexc++;
	        	}
	      	}
	    }  
	}

	return;
}

/****************************************************************************/
// Create lists and counts of the allowed final states for each 
// of the possible initial state and energyConsWindow

void storeAllowedFinalStatesNonIntAR(long *finalHoleList, long *nFinalHoles, long *finalElecList, long *nFinalElecs, 
				double *evalai, long nHotHoles, long nHotElecs, 
				nonIntBiexc *initBiexcStates, long nBiexcitons, double *deltaE, long nConsWindows) {
	FILE *pf;
	long i, a, iDeltaE, iBiexc, iBiexcDeltaE;

	pf = fopen("finalExcStates.dat", "w");
	for (iDeltaE = 0; iDeltaE < nConsWindows; iDeltaE++) {
		fprintf(pf, "DeltaE index = %d has window = % .6f\n", iDeltaE, deltaE[iDeltaE]);
		for (iBiexc = 0; iBiexc < nBiexcitons; iBiexc++) {
			iBiexcDeltaE = iBiexc + iDeltaE*nBiexcitons;
			for (i = 0; i < nHotHoles; i++) { // loop over all final hot hole states (delta_ab so just over i)
				// check to see if the final hole state is in energy window for the given initial biexciton energy
				if (fabs(-evalai[i] + initBiexcStates[iBiexc].bEnergy - initBiexcStates[iBiexc].energy) <= deltaE[iDeltaE]) {
					finalHoleList[nFinalHoles[iBiexcDeltaE] + iBiexcDeltaE*nHotHoles] = i;
					nFinalHoles[iBiexcDeltaE]++;
				}
			}
			for (a = 0; a < nHotElecs; a++) { // loop over all final hot elec states (delta_ij so just over a)
				// check to see if the final elec state is in energy window for the given initial biexciton energy
				if (fabs(evalai[a+nHotHoles] + initBiexcStates[iBiexc].jEnergy - initBiexcStates[iBiexc].energy) <= deltaE[iDeltaE]) {
					finalElecList[nFinalElecs[iBiexcDeltaE] + iBiexcDeltaE*nHotElecs] = a+nHotHoles; 
					nFinalElecs[iBiexcDeltaE]++;
				}
			}
			fprintf(pf, "\nBiexciton index = %d\n", initBiexcStates[iBiexc].index);
			fprintf(pf, "The number of hot holes satisfying energy conservation = %d\n", nFinalHoles[iBiexcDeltaE]);
			fprintf(pf, "The number of hot elecs satisfying energy conservation = %d\n", nFinalElecs[iBiexcDeltaE]);  
		}
		writeSeparation(pf);
	}
	fclose(pf);

  return;
}

/****************************************************************************/
// 

long storeNonIntEigenstates(double *psibe, double *evalbe, double *sigebe, par_st par, lng_st ist) {
	FILE *pf;
	long i, ieof, nStatesInEval, qpIndex, eigenstateList[ist.nHolesPlusElecs]; 
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
	if ((tmpQPEnergy = (double *) calloc(nStatesInEval, sizeof(double))) == NULL) nerror("tmpQPEnergy");
	if ((tmpQPSigma = (double *) calloc(nStatesInEval, sizeof(double))) == NULL) nerror("tmpQPSigma");

	// Store all the eigenvalues
	pf = fopen("eval.par", "r");
	for (i = 0; i < nStatesInEval; i++) {
		fscanf(pf, "%ld %lg %lg", &qpIndex, &tmpQPEnergy[i], &tmpQPSigma[i]);
	}
	fclose(pf);

	// Store the hole and elec energies in evalbe and create a list of their indexes
	nStoredHoles = 0;
	for (i = ist.homoIndex; i >= 0; i--) {
		if (tmpQPSigma[i] < par.sigmaCutoff) {
			nStoredHoles++;
			evalbe[ist.nHoles-nStoredHoles] = tmpQPEnergy[i];
			sigebe[ist.nHoles-nStoredHoles] = tmpQPSigma[i];
			eigenstateList[ist.nHoles-nStoredHoles] = i;
		}
		if (nStoredHoles == ist.nHoles) {
			break;
		}
	}
	nStoredElecs = 0;
	for (i = ist.lumoIndex; i < nStatesInEval; i++) {
		if (tmpQPSigma[i] < par.sigmaCutoff) {
			evalbe[ist.nHoles+nStoredElecs] = tmpQPEnergy[i];
			sigebe[ist.nHoles+nStoredElecs] = tmpQPSigma[i];
			eigenstateList[ist.nHoles+nStoredElecs] = i;
			nStoredElecs++;
		}
		if (nStoredElecs == ist.nElecs) {
			break;
		}
	}

	// Store the eigenstates from psi.par in psibe
	nStoredStates = 0;
	if ( access("psi.par", F_OK) != -1 ) {
		pf = fopen("psi.par" , "r");
		for (i = 0; i < nStatesInEval; i++) {
			fread(&psibe[nStoredStates*ist.nGridPoints], sizeof(double), ist.nGridPoints, pf);
			if (i == eigenstateList[nStoredStates]) {
				nStoredStates++;
			}
			if (nStoredStates == ist.nHolesPlusElecs) {
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
 	pf = fopen("evalEigOnly.dat", "w");
 	for (i = 0; i < ist.nHolesPlusElecs; i++) {
		fprintf(pf, "%ld % .12f %g\n", i, evalbe[i], sigebe[i]);
	}
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

/****************************************************************************/
// Removes the noninteracting states from psiai that are not eigenstates
// and orders them from lowest energy to highest energy starting at 0
// and returns the new total number of states

long remNonEigenstates(double *psi, double *eval, double *sige, double sigmaCutoff, long nGridPoints, long nOrigStates) {
	long i, nNewStates = 0;

	for (i = 0; i < nOrigStates; i++) {
		if (sige[i] < sigmaCutoff) {
			if (i == nNewStates) {
				nNewStates++;
			}
			else {
				sige[nNewStates] = sige[i];
				eval[nNewStates] = eval[i];
				memcpy(&(psi[nNewStates*nGridPoints]), &(psi[i*nGridPoints]), nGridPoints*sizeof(psi[0]));
				nNewStates++;
			}
		}
	}

	return nNewStates;
}

/****************************************************************************/
// Removes the noninteracting states from psiai that are not eigenstates
// and orders them from lowest energy to highest energy starting at 0
// and returns the new total number of states

long remNonIntStatesOutsideERange(double *psi, double *eval, double *sige, double minE, double maxE, long nGridPoints, long nOrigStates) {
	long i, nNewStates = 0;

	for (i = 0; i < nOrigStates; i++) {
		if (eval[i] > minE && eval[i] < maxE) {
			if (i == nNewStates) {
				nNewStates++;
			}
			else {
				sige[nNewStates] = sige[i];
				eval[nNewStates] = eval[i];
				memcpy(&(psi[nNewStates*nGridPoints]), &(psi[i*nGridPoints]), nGridPoints*sizeof(psi[0]));
				nNewStates++;
			}
		}
	}

	return nNewStates;
}

/****************************************************************************/
// Removes the noninteracting states from psiai that are not eigenstates
// and orders them from lowest energy to highest energy starting at 0
// and returns the new total number of states

long remNonIntStatesInsideERange(double *psi, double *eval, double *sige, double minE, double maxE, long nGridPoints, long nOrigStates) {
	long i, nNewStates = 0;

	for (i = 0; i < nOrigStates; i++) {
		if (eval[i] < minE || eval[i] > maxE) {
			if (i == nNewStates) {
				nNewStates++;
			}
			else {
				sige[nNewStates] = sige[i];
				eval[nNewStates] = eval[i];
				memcpy(&(psi[nNewStates*nGridPoints]), &(psi[i*nGridPoints]), nGridPoints*sizeof(psi[0]));
				nNewStates++;
			}
		}
	}

	return nNewStates;
}

/****************************************************************************/
// Removes the noninteracting states from psiai that are not eigenstates
// and that are not in the needed energy range and orders them 
// from lowest energy to highest energy and returns the new total number of states

long remUnnecessaryHotNonIntStates(double *psiai, double *evalai, double *sigeai, par_st par, lng_st ist) {
	FILE *pf;
	long i, nNewStates;

	// Write beginning of function
	writeSeparation(stdout);
	fprintf(stdout, "Removing unnessary hot quasiparticle states\n\n");
	fprintf(stdout, "The original number of hot hole and elec states = %ld\n", ist.nHotNonIntStates);
 
 	// Remove the states that are not eigenstates (i.e. states with sigmas above par.sigmaCutoff)
	nNewStates = remNonEigenstates(psiai, evalai, sigeai, par.sigmaCutoff, ist.nGridPoints, ist.nHotNonIntStates);
	printf("The number of hot hole and elec states after removing non eigenstates = %ld\n", nNewStates);

	// Remove states that are outside of Eimin and Eamax
	nNewStates = remNonIntStatesOutsideERange(psiai, evalai, sigeai, par.Eimin, par.Eamax, ist.nGridPoints, nNewStates);
	printf("The number of hot hole and elec states after removing states outside energy region = %ld\n", nNewStates);

	// Remove states that are inside Eimax and Eamin
	nNewStates = remNonIntStatesInsideERange(psiai, evalai, sigeai, par.Eimax, par.Eamin, ist.nGridPoints, nNewStates);
	printf("The number of hot hole and elec states after removing states inside energy region = %ld\n", nNewStates);
	fflush(stdout);

	// Write the resulting smaller evalai file
	pf = fopen("evalaiEigOnly.dat", "w");
	for (i = 0; i < nNewStates; i++) {
		fprintf(pf, "%ld % .12f %g\n", i, evalai[i], sigeai[i]);
	}
	fclose(pf);

	// Write a psi file containing just these required eigenstates
	pf = fopen("psiaiEigOnly.dat", "w");
	fwrite(psiai, sizeof(double), nNewStates*ist.nGridPoints, pf);
	fclose(pf);

	// Return the new number of hot hole and elec eigenstates
	return nNewStates;
}

/****************************************************************************/
// Removes the noninteracting states from psiai and evalai for which Coulomb
// matrix elements do not need to be calculated to get an AR lifetime. 
// Used when stochastically sampling the final hole and elec states

long remStatesNotInThetaLists(double *psiai, double *evalai, long *thetaIList, long *thetaAList, lng_st *ist) {
	FILE *pf;
	long i, a, iThetaI, iThetaA, nNewHotHoles, nNewHotElecs;
	long *uniqueThetaIList, *uniqueThetaAList;

	// Dynamically allocate memory
	if ((uniqueThetaIList = (long *) calloc(ist->itot, sizeof(long))) == NULL) nerror("uniqueThetaIList");
	if ((uniqueThetaAList = (long *) calloc(ist->atot, sizeof(long))) == NULL) nerror("uniqueThetaAList");

	// Create unique list of the indexes of the hot holes in the thetaIList
	nNewHotHoles = 0;
	for (i = 0; i < ist->itot; i++) {
		for (iThetaI = 0; iThetaI < ist->nStochHotHoles*ist->nBiexcitons*ist->nConsWindows; iThetaI++) {
			if (i == thetaIList[iThetaI]) {
				uniqueThetaIList[nNewHotHoles] = i;
				nNewHotHoles++;
				break;
			}
		}
	}

	// Create unique list of the indexes of the hot elecs in the thetaAList
	nNewHotElecs = 0;
	for (a = ist->itot; a < ist->itot+ist->atot; a++) {
		for (iThetaA = 0; iThetaA < ist->nStochHotElecs*ist->nBiexcitons*ist->nConsWindows; iThetaA++) {
			if (a == thetaAList[iThetaA]) {
				uniqueThetaAList[nNewHotElecs] = a;	
				nNewHotElecs++;
				break;
			}
		}
	}		

	// Go through uniquethetaIList and remove hot holes from psiai/evalai that are not in the list
	for (i = 0; i < nNewHotHoles; i++) {
		if (i == uniqueThetaIList[i]) {
			continue;
		}
		else {
			evalai[i] = evalai[uniqueThetaIList[i]];
			memcpy(&(psiai[i*ist->nGridPoints]), &(psiai[uniqueThetaIList[i]*ist->nGridPoints]), ist->nGridPoints*sizeof(psiai[0]));
		}
		// Update thetaIList to reflect new position of the hot hole in psiai/ evalai
		for (iThetaI = 0; iThetaI < ist->nStochHotHoles*ist->nBiexcitons*ist->nConsWindows; iThetaI++) {
			if (thetaIList[iThetaI] == uniqueThetaIList[i]) {
				thetaIList[iThetaI] = i;
			}
		}
	}
	// Go through uniquethetaAList and remove hot elecs from psiai/evalai that are not in the list
	for (a = 0; a < nNewHotElecs; a++) {
		if ((nNewHotHoles + a) == uniqueThetaAList[a]) {
			continue;
		}
		else {
			evalai[nNewHotHoles+a] = evalai[uniqueThetaAList[a]];
			memcpy(&(psiai[(nNewHotHoles+a)*ist->nGridPoints]), &(psiai[uniqueThetaAList[a]*ist->nGridPoints]), ist->nGridPoints*sizeof(psiai[0]));
		}
		// Update thetaAList to reflect new position of the hot elec in psiai/ evalai
		for (iThetaA = 0; iThetaA < ist->nStochHotElecs*ist->nBiexcitons*ist->nConsWindows; iThetaA++) {
			if (thetaAList[iThetaA] == uniqueThetaAList[a]) {
				thetaAList[iThetaA] = (a + nNewHotHoles);
			}
		}
	}

	// Print indexes of hot carrier states that Coulomb matrix elements will be calculated 
	pf = fopen("uniqueHotNonIntStates.dat", "w");
	for (i = 0; i < nNewHotHoles+nNewHotElecs; i++) {
		if (i < nNewHotHoles) a = uniqueThetaIList[i];
		else a = uniqueThetaAList[i-nNewHotHoles];
		fprintf(pf, "%ld %ld % .6f\n", i, a, evalai[i]);
	}
	fclose(pf);
	printf("The number of unique hot holes stochastically sampled = %ld\n", nNewHotHoles);
	printf("The number of unique hot elecs stochastically sampled = %ld\n", nNewHotElecs);
	fflush(stdout);

	// Update ist structure
	ist->itot = nNewHotHoles;
	ist->atot = nNewHotElecs;
	ist->nHotNonIntStates = ist->itot+ist->atot;

	// Free dynamically allocated memory
	free(uniqueThetaAList); free(uniqueThetaIList);

	return (nNewHotHoles+nNewHotElecs);
}

/****************************************************************************/