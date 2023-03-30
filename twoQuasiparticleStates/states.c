/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "qp.h"

/****************************************************************************/
// 

long fillNonintTwoQPStructure(nonintTwoQPState *nonintTwoQP, nonintQP *holeQP, nonintQP *elecQP, lParams lPar) {
	FILE *pf;
	long iHole1, iHole2, iElec1, iElec2, iTwoQP = 0;
	long singletFlag;

	// Determine if iElec2 (iHole2) can equal iElec1 (iHole1)
	if (lPar.calcSinglets) {
		singletFlag = 0;
	}
	else {
		singletFlag = 1; // only used in elec-elec and hole-hole
	}

	// Fill the noninteracting two quasiparticle structures
	if (! strcmp(lPar.calcType, "elec-hole") || ! strcmp(lPar.calcType, "sp-elec-hole") ||
		! strcmp(lPar.calcType, "auger-decay-ni") || ! strcmp(lPar.calcType, "auger-decay-i")) {
		for (iHole1 = 0; iHole1 < lPar.nHoles; iHole1++) {
			for (iElec1 = 0; iElec1 < lPar.nElecs; iElec1++) {
				nonintTwoQP[iTwoQP].index = iTwoQP;
				nonintTwoQP[iTwoQP].energy = (elecQP[iElec1].energy - holeQP[iHole1].energy);
				nonintTwoQP[iTwoQP].qp1 = &holeQP[iHole1];
				nonintTwoQP[iTwoQP].qp2 = &elecQP[iElec1];
				if (! strcmp(lPar.calcType, "sp-elec-hole")) {
					nonintTwoQP[iTwoQP].spinZ = (elecQP[iElec1].spinZ + holeQP[iHole1].spinZ);
					if (fabs(nonintTwoQP[iTwoQP].spinZ) > EPS) {
						nonintTwoQP[iTwoQP].spin = 1.0; // spin triplet
					}
					else {
						nonintTwoQP[iTwoQP].spin = 0.0; // spin singlet or triplet
					}
				}
				iTwoQP++;
			}
		}
	}
	else if (! strcmp(lPar.calcType, "elec-elec")) {
		for (iElec1 = 0; iElec1 < lPar.nElecs; iElec1++) {
			for (iElec2 = iElec1+singletFlag; iElec2 < lPar.nElecs; iElec2++) {
				nonintTwoQP[iTwoQP].index = iTwoQP;
				nonintTwoQP[iTwoQP].energy = (elecQP[iElec2].energy + elecQP[iElec1].energy);
				nonintTwoQP[iTwoQP].qp1 = &elecQP[iElec1];
				nonintTwoQP[iTwoQP].qp2 = &elecQP[iElec2];
				iTwoQP++;
			}
		}
	}
	else if (! strcmp(lPar.calcType, "hole-hole")) {
		for (iHole1 = 0; iHole1 < lPar.nHoles; iHole1++) {
			for (iHole2 = iHole1+singletFlag; iHole2 < lPar.nHoles; iHole2++) {
				nonintTwoQP[iTwoQP].index = iTwoQP;
				nonintTwoQP[iTwoQP].energy = -(holeQP[iHole1].energy + holeQP[iHole2].energy);
				nonintTwoQP[iTwoQP].qp1 = &holeQP[iHole1];
				nonintTwoQP[iTwoQP].qp2 = &holeQP[iHole2];
				iTwoQP++;
			}
		}
	} 

	// Print the noninteracting two quasiparticle pair state information
	pf = fopen("nonintTwoQuasiparticleStates.dat", "w");
	for (iTwoQP = 0; iTwoQP < lPar.nNonintTwoQPStates; iTwoQP++) {
		fprintf(pf, "%ld %.6f %ld %.6f %ld %.6f\n", nonintTwoQP[iTwoQP].index, nonintTwoQP[iTwoQP].energy, 
				nonintTwoQP[iTwoQP].qp1->index, nonintTwoQP[iTwoQP].qp1->energy, 
				nonintTwoQP[iTwoQP].qp2->index, nonintTwoQP[iTwoQP].qp2->energy);
	}
	fclose(pf);

	// Print the noninteracting two quasiparticle pair state information
	if (! strcmp(lPar.calcType, "sp-elec-hole")) {
		pf = fopen("spNonintTwoQuasiparticleStates.dat", "w");
		for (iTwoQP = 0; iTwoQP < lPar.nNonintTwoQPStates; iTwoQP++) {
			fprintf(pf, "%ld %.2f % .2f %.6f  %ld %.2f % .2f %.6f  %ld %.2f % .2f %.6f\n", 
					nonintTwoQP[iTwoQP].index, nonintTwoQP[iTwoQP].spin, nonintTwoQP[iTwoQP].spinZ, nonintTwoQP[iTwoQP].energy, 
					nonintTwoQP[iTwoQP].qp1->index, nonintTwoQP[iTwoQP].qp1->spin, nonintTwoQP[iTwoQP].qp1->spinZ, nonintTwoQP[iTwoQP].qp1->energy,
					nonintTwoQP[iTwoQP].qp2->index, nonintTwoQP[iTwoQP].qp2->spin, nonintTwoQP[iTwoQP].qp2->spinZ, nonintTwoQP[iTwoQP].qp2->energy);
		}
		fclose(pf);
	}

	return 0;
}

/****************************************************************************/
// 

long fillIntExcitonStructure(intExc *intExciton, double *Cai, double *Ebs, lParams lPar) {
	long i;

	for (i = 0; i < lPar.nIntExcitons; i++) {
		intExciton[i].index = i;
		intExciton[i].energy = Ebs[i];
		intExciton[i].Cai = &Cai[i*lPar.nNonintExcitons];
	}

	return 0;
}

/****************************************************************************/
// 

long fillNonintExcitonStructure(nonintExc *nonintExciton, nonintQP *holeQP, nonintQP *elecQP, lParams lPar) {
	FILE *pf;
	long iHole, iElec, iExciton = 0;

	// Fill the noninteracting exciton structure
	for (iHole = 0; iHole < lPar.nHoles; iHole++) {
		for (iElec = 0; iElec < lPar.nElecs; iElec++) {
			nonintExciton[iExciton].index = iExciton;
			nonintExciton[iExciton].energy = (elecQP[iElec].energy - holeQP[iHole].energy);
			nonintExciton[iExciton].h = &holeQP[iHole];
			nonintExciton[iExciton].e = &elecQP[iElec];
			iExciton++;
		}
	}

	// Print out information summarizing the noninteracting excitonic states
	pf = fopen("nonintExcitons.dat", "w");
	for (iExciton = 0; iExciton < lPar.nNonintExcitons; iExciton++) {
		writeNonintExcitonState(nonintExciton[iExciton], pf);
	}

	fclose(pf);


	return 0;
}

/****************************************************************************/
// 

long storeNonintEigenstates(double *psibe, double *evalbe, double *sigebe, dParams par, lParams ist) {
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
	for (i = ist.iHomo; i >= 0; i--) {
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
	for (i = ist.iLumo; i < nStatesInEval; i++) {
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

// long remUnnecessaryHotNonIntStates(double *psiai, double *evalai, double *sigeai, par_st par, lng_st ist) {
// 	FILE *pf;
// 	long i, nNewStates;

// 	// Write beginning of function
// 	writeSeparation(stdout);
// 	fprintf(stdout, "Removing unnessary hot quasiparticle states\n\n");
// 	fprintf(stdout, "The original number of hot hole and elec states = %ld\n", ist.nHotNonIntStates);
 
//  	// Remove the states that are not eigenstates (i.e. states with sigmas above par.sigmaCutoff)
// 	nNewStates = remNonEigenstates(psiai, evalai, sigeai, par.sigmaCutoff, ist.nGridPoints, ist.nHotNonIntStates);
// 	printf("The number of hot hole and elec states after removing non eigenstates = %ld\n", nNewStates);

// 	// Remove states that are outside of Eimin and Eamax
// 	nNewStates = remNonIntStatesOutsideERange(psiai, evalai, sigeai, par.Eimin, par.Eamax, ist.nGridPoints, nNewStates);
// 	printf("The number of hot hole and elec states after removing states outside energy region = %ld\n", nNewStates);

// 	// Remove states that are inside Eimax and Eamin
// 	nNewStates = remNonIntStatesInsideERange(psiai, evalai, sigeai, par.Eimax, par.Eamin, ist.nGridPoints, nNewStates);
// 	printf("The number of hot hole and elec states after removing states inside energy region = %ld\n", nNewStates);
// 	fflush(stdout);

// 	// Write the resulting smaller evalai file
// 	pf = fopen("evalaiEigOnly.dat", "w");
// 	for (i = 0; i < nNewStates; i++) {
// 		fprintf(pf, "%ld % .12f %g\n", i, evalai[i], sigeai[i]);
// 	}
// 	fclose(pf);

// 	// Write a psi file containing just these required eigenstates
// 	pf = fopen("psiaiEigOnly.dat", "w");
// 	fwrite(psiai, sizeof(double), nNewStates*ist.nGridPoints, pf);
// 	fclose(pf);

// 	// Return the new number of hot hole and elec eigenstates
// 	return nNewStates;
// }

/****************************************************************************/
