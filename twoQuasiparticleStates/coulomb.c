/*****************************************************************************/
//
// This file deals with Coulomb potential 
//
/*****************************************************************************/

#include "qp.h"

/****************************************************************************/
// 

long calcAllCoulombMatrixElements(double *Wrsut, double *Vrsut, nonintTwoQPState *nonintTwoQP, double *psiHoles, double *psiElecs, 
					zomplex *qSpaceCoulombPot, zomplex *qSpaceScreenedCoulombPot, grid qSpaceGrid, grid rSpaceGrid, 
					lParams lPar, fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwPsi) { 
	FILE *pf;
	long iNonintTwoQP1, iNonintTwoQP2, iGrid, iThread, index;
	long nGridPoints, nReadInCME, nCME, nNonintTwoQPStatesSquared;
	long iR, iS, iU, iT, nR, nS, nU, nT;
	long *rsutIndexToCMEIndexList;
	double *psiR, *psiS, *psiU, *psiT;
	double *hartreePotential, sum;
	zomplex *rho, *tmpCoulombPot;
	coulombMatrixElement *coulombMatElement;
	cmeController WrsutController, VrsutController;

	// Useful integers
	if (! strcmp(lPar.calcType, "elec-elec")) {
		nR = nS = nU = nT = lPar.nElecs;
		nCME = lPar.nElecs*lPar.nElecs*lPar.nElecs*lPar.nElecs;
	}
	else if (! strcmp(lPar.calcType, "hole-hole")) {
		nR = nS = nU = nT = lPar.nHoles;
		nCME = lPar.nHoles*lPar.nHoles*lPar.nHoles*lPar.nHoles;
	}
	else {
		nCME = lPar.nNonintTwoQPStates*lPar.nNonintTwoQPStates; // = nHoles*nElecs*nHoles*nElecs
	} 
	nNonintTwoQPStatesSquared = lPar.nNonintTwoQPStates*lPar.nNonintTwoQPStates;
	nGridPoints = rSpaceGrid.nGridPoints;

	// Write beginning of function
	writeSeparation(stdout); writeCurrentTime(stdout);
	fprintf(stdout, "Beginnning the calculation of all the Coulomb matrix elements\n\n");
	fprintf(stdout, "The maximum number of Coulomb matrix elements that must be calculated = %ld\n", 2*nCME); fflush(stdout);

	// Dynamically allocate memory
	if ((hartreePotential = (double *) calloc(lPar.nThreads*nGridPoints, sizeof(double))) == NULL) memoryError("hartreePotential");
	if ((rho = (zomplex *) calloc(lPar.nThreads*nGridPoints, sizeof(zomplex))) == NULL) memoryError("rho");
	if ((tmpCoulombPot = (zomplex *) calloc(rSpaceGrid.nGridPoints, sizeof(zomplex))) == NULL) memoryError("tmpCoulombPot");
	if ((coulombMatElement = (coulombMatrixElement *) calloc(nCME, sizeof(coulombMatrixElement))) == NULL) memoryError("coulombMatElement");
	if ((rsutIndexToCMEIndexList = (long *) calloc(nCME, sizeof(long))) == NULL) memoryError("rsutIndexToCMEIndexList");

	// Initialize which Coulomb matrix elements structure, sets alreadyCalculatedFlag equal to 0 (False)
	initCoulombMatrixElementStruct(coulombMatElement, nCME);

	// Fill the Coulomb matrix element structures
	for (iNonintTwoQP1 = 0; iNonintTwoQP1 < lPar.nNonintTwoQPStates; iNonintTwoQP1++) {
		for (iNonintTwoQP2 = 0; iNonintTwoQP2 < lPar.nNonintTwoQPStates; iNonintTwoQP2++) {
			index = (iNonintTwoQP1*lPar.nNonintTwoQPStates + iNonintTwoQP2);
			coulombMatElement[index].index = index;
			coulombMatElement[index].iU   = nonintTwoQP[iNonintTwoQP1].qp1->index; // hole for elec-hole and hole-hole
			coulombMatElement[index].psiU = nonintTwoQP[iNonintTwoQP1].qp1->psi;
			coulombMatElement[index].iR   = nonintTwoQP[iNonintTwoQP1].qp2->index; // elec for elec-hole and elec-elec
			coulombMatElement[index].psiR = nonintTwoQP[iNonintTwoQP1].qp2->psi; 
			coulombMatElement[index].iT   = nonintTwoQP[iNonintTwoQP2].qp1->index; // hole for elec-hole and hole-hole
			coulombMatElement[index].psiT = nonintTwoQP[iNonintTwoQP2].qp1->psi;
			coulombMatElement[index].iS   = nonintTwoQP[iNonintTwoQP2].qp2->index; // elec for elec-hole and elec-elec
			coulombMatElement[index].psiS = nonintTwoQP[iNonintTwoQP2].qp2->psi; 
		}
	}

	if (! strcmp(lPar.calcType, "elec-elec") || ! strcmp(lPar.calcType, "hole-hole")) {
		coulombMatElement[index+1].iU = coulombMatElement[index+1].iR = coulombMatElement[index+1].iT = coulombMatElement[index+1].iS = 0;
		coulombMatElement[index+1].psiU = coulombMatElement[index+1].psiR = coulombMatElement[index+1].psiT = coulombMatElement[index+1].psiS = nonintTwoQP[0].qp1->psi;
		coulombMatElement[index+2].iU = coulombMatElement[index+2].iR = coulombMatElement[index+2].iT = coulombMatElement[index+2].iS = (nR-1);
		coulombMatElement[index+2].psiU = coulombMatElement[index+2].psiR = coulombMatElement[index+2].psiT = coulombMatElement[index+2].psiS = nonintTwoQP[nR-1].qp2->psi;
	}

	// Initialize the Coulomb matrix element controller structure
	initCoulombMatrixElementController(&WrsutController, coulombMatElement, nCME);

	// Calculate and store the rsutIndex for all the Coulomb matrix elements
	calcAllCoulombMatrixElementsRSUTIndices(coulombMatElement, WrsutController.nCME, WrsutController.nR, 
											WrsutController.nS, WrsutController.nU, WrsutController.nT);

	// Calculate list that maps rsutIndex to the order of Coulomb matrix elements 
	calcRSUTIndexToCMEIndexList(rsutIndexToCMEIndexList, coulombMatElement, WrsutController.nCME);

	// Determine if matrix elements have been previously calculated and read them in if so
	if (lPar.readInCoulombMatrixElements) {
		nReadInCME = readCoulombMatrixElements(coulombMatElement, rsutIndexToCMEIndexList, WrsutController, "Wrsut.par");
		fprintf(stdout, "The number of hartree-like Coulomb matrix elements that were read in = %ld\n", nReadInCME);
	}
	else {
		nReadInCME = 0;
	}

	// Calculated all the missing screened hartree-like Coulomb matrix elements using openMP
	calcCoulombMatrixElementsOpenMP(coulombMatElement, nCME, rsutIndexToCMEIndexList, qSpaceScreenedCoulombPot, nGridPoints, 
									rSpaceGrid.dV, lPar.nThreads, planfw, planbw, fftwPsi, "RS", 1, "Wrsut.dat");

	// Set openMP settings
	omp_set_dynamic(0);
	omp_set_num_threads(lPar.nThreads);

	// Calculate Wrsut, direct Coulomb term
	pf = fopen("subset_Wrsut.dat", "w");
#pragma omp parallel for private(iNonintTwoQP1, iNonintTwoQP2, iThread, index, psiR, psiS, psiU, psiT, iGrid, sum)
//#pragma omp parallel for private(iNonintTwoQP1, iNonintTwoQP2, index)
	for (iNonintTwoQP1 = 0; iNonintTwoQP1 < lPar.nNonintTwoQPStates; iNonintTwoQP1++) {
		psiU = nonintTwoQP[iNonintTwoQP1].qp1->psi; // hole for elec-hole and hole-hole
		psiR = nonintTwoQP[iNonintTwoQP1].qp2->psi; // elec for elec-hole and elec-elec
		for (iNonintTwoQP2 = 0; iNonintTwoQP2 < lPar.nNonintTwoQPStates; iNonintTwoQP2++) {
			iThread = omp_get_thread_num();
			index = iNonintTwoQP1*lPar.nNonintTwoQPStates + iNonintTwoQP2;
			psiT = nonintTwoQP[iNonintTwoQP2].qp1->psi; // hole for elec-hole and hole-hole
			psiS = nonintTwoQP[iNonintTwoQP2].qp2->psi; // elec for elec-hole and elec-elec
			if (coulombMatElement[index].alreadyCalculatedFlag) { // this should always be true for (sp-)elec-hole calcs
				Wrsut[index] = coulombMatElement[index].value;
			}
			else {
				for (iGrid = 0; iGrid < nGridPoints; iGrid++) {
					rho[iThread*nGridPoints + iGrid].re = (psiR[iGrid] * psiS[iGrid]);
					rho[iThread*nGridPoints + iGrid].im = 0.0;
				}
				// Screened Coulomb matrix elements, Wrsut
				hartree(&rho[iThread*nGridPoints], qSpaceScreenedCoulombPot, &hartreePotential[iThread*nGridPoints], 
						nGridPoints, planfw[iThread], planbw[iThread], &fftwPsi[iThread*nGridPoints]);
				sum = 0.0;
				for (iGrid = 0; iGrid < nGridPoints; iGrid++) {
					sum += (psiU[iGrid] * psiT[iGrid] * hartreePotential[iThread*nGridPoints + iGrid]);
				}
				Wrsut[index] = sum*rSpaceGrid.dV;
			}
			fprintf(pf, "%ld %ld %ld %ld %ld %ld %.16f %.16f %.16f\n", 
							nonintTwoQP[iNonintTwoQP1].qp1->index, nonintTwoQP[iNonintTwoQP1].qp2->index, 
							nonintTwoQP[iNonintTwoQP2].qp1->index, nonintTwoQP[iNonintTwoQP2].qp2->index, iNonintTwoQP1, 
							iNonintTwoQP2, nonintTwoQP[iNonintTwoQP1].energy, nonintTwoQP[iNonintTwoQP2].energy, Wrsut[index]);
			fflush(pf);
			if (iNonintTwoQP1 == iNonintTwoQP2) {
				nonintTwoQP[iNonintTwoQP1].hartreeEnergy = Wrsut[index];
				nonintTwoQP[iNonintTwoQP1].exchangeEnergy = 0.0; // filled in below
				if (! strcmp(lPar.calcType, "elec-hole") || ! strcmp(lPar.calcType, "sp-elec-hole")) {
					nonintTwoQP[iNonintTwoQP1].hartreeEnergy *= -1.0;
				}
			}
		}
	}
	fclose(pf);

	// Calculate Vrsut, indirect Coulomb term (exchange-like)
	if ((lPar.calcSinglets && ! strcmp(lPar.calcType, "elec-hole")) || ! strcmp(lPar.calcType, "sp-elec-hole") ||
		(lPar.calcTriplets && (! strcmp(lPar.calcType, "elec-elec") || ! strcmp(lPar.calcType, "hole-hole")))) {
		if (! strcmp(lPar.calcType, "elec-hole") || ! strcmp(lPar.calcType, "sp-elec-hole")) { 
			pf = fopen("subset_Vrsut.dat", "w"); // bare Coulomb
			memcpy(&tmpCoulombPot[0], &qSpaceCoulombPot[0], rSpaceGrid.nGridPoints*sizeof(zomplex));
		}
		else { // elec-elec or hole-hole calculations
			// TODO: add in calcType = "auger-decay-ni" and "auger-decay-i" 
			//		 and decide if the indirect term should be screened or not
			pf = fopen("subset_Wrust.dat", "w"); // screened Coulomb
			memcpy(&tmpCoulombPot[0], &qSpaceScreenedCoulombPot[0], rSpaceGrid.nGridPoints*sizeof(zomplex));
		}
		resetAlreadyCalculatedFlagCME(coulombMatElement, nCME);
		// Determine if matrix elements have been previously calculated and read them in if so
		if (lPar.readInCoulombMatrixElements) {
			nReadInCME = readCoulombMatrixElements(coulombMatElement, rsutIndexToCMEIndexList, WrsutController, "Vrsut.par");
			fprintf(stdout, "The number of exchange-like Coulomb matrix elements that were read in = %ld\n", nReadInCME);
		}
		calcCoulombMatrixElementsOpenMP(coulombMatElement, nCME, rsutIndexToCMEIndexList, tmpCoulombPot, nGridPoints, 
											rSpaceGrid.dV, lPar.nThreads, planfw, planbw, fftwPsi, "RU", 1, "Vrsut.dat");
	#pragma omp parallel for private(iNonintTwoQP1, iNonintTwoQP2, iThread, index, psiR, psiS, psiU, psiT, iGrid, sum)
	//#pragma omp parallel for private(iNonintTwoQP1, iNonintTwoQP2, index)
		for (iNonintTwoQP1 = 0; iNonintTwoQP1 < lPar.nNonintTwoQPStates; iNonintTwoQP1++) {
			psiU = nonintTwoQP[iNonintTwoQP1].qp1->psi; // hole for elec-hole and hole-hole
			psiR = nonintTwoQP[iNonintTwoQP1].qp2->psi; // elec for elec-hole and elec-elec
			for (iNonintTwoQP2 = 0; iNonintTwoQP2 < lPar.nNonintTwoQPStates; iNonintTwoQP2++) {
				iThread = omp_get_thread_num();
				index = iNonintTwoQP1*lPar.nNonintTwoQPStates + iNonintTwoQP2;
				psiT = nonintTwoQP[iNonintTwoQP2].qp1->psi; // hole for elec-hole and hole-hole
				psiS = nonintTwoQP[iNonintTwoQP2].qp2->psi; // elec for elec-hole and elec-elec
				if (coulombMatElement[index].alreadyCalculatedFlag) {
					Vrsut[index] = coulombMatElement[index].value;
				}
				else {
					for (iGrid = 0; iGrid < nGridPoints; iGrid++) {
						rho[iThread*nGridPoints + iGrid].re = (psiR[iGrid] * psiU[iGrid]);
						rho[iThread*nGridPoints + iGrid].im = 0.0;
					}
					// Coulomb matrix elements, Vrsut
					hartree(&rho[iThread*nGridPoints], tmpCoulombPot, &hartreePotential[iThread*nGridPoints], 
							nGridPoints, planfw[iThread], planbw[iThread], &fftwPsi[iThread*nGridPoints]);
					sum = 0.0;
					for (iGrid = 0; iGrid < nGridPoints; iGrid++) {
						sum += (psiT[iGrid] * psiS[iGrid] * hartreePotential[iThread*nGridPoints + iGrid]);
					}
					Vrsut[index] = sum*rSpaceGrid.dV;
				}
				fprintf(pf, "%ld %ld %ld %ld %ld %ld %.16f %.16f %.16f\n", nonintTwoQP[iNonintTwoQP1].qp1->index,
					nonintTwoQP[iNonintTwoQP1].qp2->index, nonintTwoQP[iNonintTwoQP2].qp1->index, 
					nonintTwoQP[iNonintTwoQP2].qp2->index, iNonintTwoQP1, iNonintTwoQP2, 
					nonintTwoQP[iNonintTwoQP1].energy, nonintTwoQP[iNonintTwoQP2].energy, Vrsut[index]);
				fflush(pf);
				if (iNonintTwoQP1 == iNonintTwoQP2) {
					nonintTwoQP[iNonintTwoQP1].exchangeEnergy = -1.0*Vrsut[index];
					if (! strcmp(lPar.calcType, "elec-hole") || ! strcmp(lPar.calcType, "sp-elec-hole")) {
						nonintTwoQP[iNonintTwoQP1].exchangeEnergy *= -2.0;
					}
				}
			}
		}
		fclose(pf);
	}

	// Write the hartree and exchange energies for all the noninteracting two quasiparticle states
	pf = fopen("nonintTwoQPHartreeExchangeEnergies.dat", "w");
	fprintf(pf, "# iRS   eR+eS    iR     eR     iS     eS       eHartree       eExchange\n");
	for (iNonintTwoQP1 = 0; iNonintTwoQP1 < lPar.nNonintTwoQPStates; iNonintTwoQP1++) {
		fprintf(pf, "%5ld %.6f %3ld %.6f %3ld %.6f %.12f %.12f\n", iNonintTwoQP1, nonintTwoQP[iNonintTwoQP1].energy, 
					nonintTwoQP[iNonintTwoQP1].qp1->index, nonintTwoQP[iNonintTwoQP1].qp1->energy, 
					nonintTwoQP[iNonintTwoQP1].qp2->index, nonintTwoQP[iNonintTwoQP1].qp2->energy,
					nonintTwoQP[iNonintTwoQP1].hartreeEnergy, nonintTwoQP[iNonintTwoQP1].exchangeEnergy);
	}
	fclose(pf);

	// // Free dynamically allocated memory
	free(hartreePotential); free(rho);
	free(coulombMatElement); free(tmpCoulombPot);
	free(rsutIndexToCMEIndexList);

	// Write ending of function
	fprintf(stdout, "\nFinished the calculation of the Coulomb matrix elements\n");
	writeCurrentTime(stdout); fflush(stdout);

	return 0;
}

/****************************************************************************/
//

long calcCoulombMatrixElementsOpenMP(coulombMatrixElement *coulombME, long nCME, long *rsutIndexToCMEIndexList, 
									zomplex *qSpaceCoulombPot, long nGridPoints, double rSpaceVolumeElement, long nThreads, 
									fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwPsi, 
									char *iRhoPair, long printCMEValuesFlag, char *cmeValuesFileName) {
	FILE *pf;
	long i, j, index, rsutIndex, iThread, iGrid, iMax, lTmp;
	long *tmpCoulombList; 
	long iR, iS, iU, iT, nR, nS, nU, nT, *nCalculatedCMEs, *nStoredCMEs;
	double *psiR, *psiS, *psiU, *psiT;
	double sum, *hartreePotential, *cmeValues;
	zomplex *rho;

	// Dynamically allocate memory
	if ((tmpCoulombList = (long *) calloc(nCME, sizeof(long))) == NULL) memoryError("tmpCoulombList");
	if ((nCalculatedCMEs = (long *) calloc(nThreads, sizeof(long))) == NULL) memoryError("nCalculatedCMEs");
	if ((nStoredCMEs = (long *) calloc(nThreads, sizeof(long))) == NULL) memoryError("nStoredCMEs");
	if ((cmeValues = (double *) calloc(nCME, sizeof(double))) == NULL) memoryError("cmeValues");
	if ((hartreePotential = (double *) calloc(nThreads*nGridPoints, sizeof(double))) == NULL) memoryError("hartreePotential");
	if ((rho = (zomplex *) calloc(nThreads*nGridPoints, sizeof(zomplex))) == NULL) memoryError("rho");

	// Determine which matrix elements need to be calculated
	nR = nS = nU = nT = 0;
	for (i = 0; i < nCME; i++) {
			if (coulombME[i].iR > nR) nR = coulombME[i].iR;
			if (coulombME[i].iS > nS) nS = coulombME[i].iS;
			if (coulombME[i].iU > nU) nU = coulombME[i].iU;
			if (coulombME[i].iT > nT) nT = coulombME[i].iT;
	}
	nR++; nS++; nU++; nT++;
	
	// Set all array of pointers to NULL 
	double *psiRList[nR+nS+nU+nT], *psiSList[nR+nS+nU+nT], *psiUList[nR+nS+nU+nT], *psiTList[nR+nS+nU+nT];
	double *tmpPsiList[nR+nS+nU+nT];
	for (i = 0; i < nR+nS+nU+nT; i++) {
		psiRList[i] = 0; psiSList[i] = 0; psiUList[i] = 0; psiTList[i] = 0; tmpPsiList[i] = 0;
	}

	// Determine order in which the matrix elements will be calculated
	for (i = 0; i < nCME; i++) { // Set up RS default
		if (coulombME[i].iR >= 0 &&  coulombME[i].iS >= 0 && coulombME[i].iU >= 0 && coulombME[i].iT >= 0) {
			if (! psiRList[coulombME[i].iR]) psiRList[coulombME[i].iR] = coulombME[i].psiR;
			if (! psiSList[coulombME[i].iS]) psiSList[coulombME[i].iS] = coulombME[i].psiS;
			if (! psiUList[coulombME[i].iU]) psiUList[coulombME[i].iU] = coulombME[i].psiU;
			if (! psiTList[coulombME[i].iT]) psiTList[coulombME[i].iT] = coulombME[i].psiT;
		}
	}

	// Create mapping from order of coulombME structure to how the matrix elements will be calculated
	// TODO: replace with memcpy
	for (iR = 0; iR < nR; iR++) {
 		for (iS = 0; iS < nS; iS++) {
 			for (iU = 0; iU < nU; iU++) {
 				for (iT = 0; iT < nT; iT++) {
 					index = iR*(nS*nU*nT) + iS*(nU*nT) + iU*(nT) + iT;
 					tmpCoulombList[index] = rsutIndexToCMEIndexList[index];
	 			}
	 		}
 		}
	}

	if (! strcmp(iRhoPair, "RU")) { // Switch U and S for exchange-like matrix elements -> RU pairing
		for (iU = 0; iU < nU; iU++) tmpPsiList[iU] = psiUList[iU];
		for (iS = 0; iS < nS; iS++) psiUList[iS] = psiSList[iS];
		for (iU = 0; iU < nU; iU++) psiSList[iU] = tmpPsiList[iU];
		lTmp = nU;  nU = nS;   nS = lTmp;
		for (iR = 0; iR < nR; iR++) {
 			for (iS = 0; iS < nS; iS++) {
 				for (iU = 0; iU < nU; iU++) {
 					for (iT = 0; iT < nT; iT++) {
 						index = iR*(nS*nU*nT) + iS*(nU*nT) + iU*(nT) + iT;
 						tmpCoulombList[index] = rsutIndexToCMEIndexList[iR*(nS*nU*nT) + iU*(nS*nT) + iS*(nT) + iT];
	 				}
	 			}
 			}
		}	
	}

	// Set openMP threads for the main computational part of the program	
    omp_set_dynamic(0);
    omp_set_num_threads(nThreads);
	
    // Main computational part of the program
	for (iR = 0; iR < nR; iR++) {
		psiR = psiRList[iR];
	#pragma omp parallel for private(sum, iThread, iGrid, index, rsutIndex, iS, psiS, iU, psiU, iT, psiT)
		for (iS = 0; iS < nS; iS++) {
			psiS = psiSList[iS];
			iThread = omp_get_thread_num();
			for (iGrid = 0; iGrid < nGridPoints; iGrid++) {
				rho[iThread*nGridPoints + iGrid].re = (psiR[iGrid] * psiS[iGrid]);
				rho[iThread*nGridPoints + iGrid].im = 0.0;
			}
			hartree(&rho[iThread*nGridPoints], qSpaceCoulombPot, &hartreePotential[iThread*nGridPoints], 
					nGridPoints, planfw[iThread], planbw[iThread], &fftwPsi[iThread*nGridPoints]);
			for (iU = 0; iU < nU; iU++) {
				psiU = psiUList[iU];
				for (iT = 0; iT < nT; iT++) {
					rsutIndex = iR*(nS*nU*nT) + iS*(nU*nT) + iU*(nT) + iT;
					index = tmpCoulombList[rsutIndex];
					if (index >= 0) {
						if (! coulombME[index].alreadyCalculatedFlag) {
							psiT = psiTList[iT];
							sum = 0.0;
							for (iGrid = 0; iGrid < nGridPoints; iGrid++) {
								sum += (psiU[iGrid] * psiT[iGrid] * hartreePotential[iThread*nGridPoints + iGrid]);
							}
							coulombME[index].value = sum*rSpaceVolumeElement;
							coulombME[index].alreadyCalculatedFlag = 1;
							nCalculatedCMEs[iThread] += 1;
							nStoredCMEs[iThread] += 1;
						}
						cmeValues[rsutIndex] = coulombME[index].value;
					}
					else {
						psiT = psiTList[iT];
						sum = 0.0;
						for (iGrid = 0; iGrid < nGridPoints; iGrid++) {
							sum += (psiU[iGrid] * psiT[iGrid] * hartreePotential[iThread*nGridPoints + iGrid]);
						}
						cmeValues[rsutIndex] = sum*rSpaceVolumeElement;
						nCalculatedCMEs[iThread] += 1;
					}
				}
			}
		}
	}

	// Print all cmeValues matrix
	if (printCMEValuesFlag) {
		pf = fopen(cmeValuesFileName, "w");
//	#pragma omp parallel for private(rsutIndex, iR, iS, iU, iT)	
		for (iR = 0; iR < nR; iR++) {
	 		for (iS = 0; iS < nS; iS++) {
	 			for (iU = 0; iU < nU; iU++) {
	 				for (iT = 0; iT < nT; iT++) {
	 					rsutIndex = iR*(nS*nU*nT) + iS*(nU*nT) + iU*(nT) + iT;
						fprintf(pf, "%ld %ld %ld %ld %.16f\n", iU, iR, iT, iS, cmeValues[rsutIndex]); fflush(pf);
		 			}
		 		}
	 		}
		}
		fclose(pf);
	}

	// Print the number of matrix elements that were calculated and stored 
	for (i = 0, j = 0, iThread = 0; iThread < nThreads; iThread++) {
		i += nCalculatedCMEs[iThread];
		j += nStoredCMEs[iThread];
	}
	fprintf(stdout, "%ld Coulomb matrix elements were calculated\n", i);
	fprintf(stdout, "%ld Coulomb matrix elements were stored\n", j); fflush(stdout);

	// Free dynamically allocated memory
	free(tmpCoulombList); free(cmeValues);
	free(nCalculatedCMEs); free(nStoredCMEs);
	free(rho); free(hartreePotential); 

	return 0;
}

/****************************************************************************/
//

long initCoulombMatrixElementController(cmeController *VrsutController, coulombMatrixElement *cme, long nCME) {
	long i;
	long iMinR, iMinS, iMinU, iMinT;
	long iMaxR, iMaxS, iMaxU, iMaxT; 

	// Set basic properties of the controller
	VrsutController->nCME = nCME;
	VrsutController->cme = cme;

	iMinR = iMinS = iMinU = iMinT = 1000000;
	iMaxR = iMaxS = iMaxU = iMaxT = -1;
	for (i = 0; i < nCME; i++) {
		if (cme[i].iR >= 0 && cme[i].iS >= 0 && cme[i].iU >= 0 && cme[i].iT >= 0) {
			if (cme[i].iR > iMaxR) { iMaxR = cme[i].iR; }
			if (cme[i].iR < iMinR) { iMinR = cme[i].iR; }
			if (cme[i].iS > iMaxS) { iMaxS = cme[i].iS; }
			if (cme[i].iS < iMinS) { iMinS = cme[i].iS; }
			if (cme[i].iU > iMaxU) { iMaxU = cme[i].iU; }
			if (cme[i].iU < iMinU) { iMinU = cme[i].iU; }
			if (cme[i].iT > iMaxT) { iMaxT = cme[i].iT; }
			if (cme[i].iT < iMinT) { iMinT = cme[i].iT; }
		}
	}
	VrsutController->nR = (iMaxR - iMinR + 1);
	VrsutController->nS = (iMaxS - iMinS + 1);
	VrsutController->nU = (iMaxU - iMinU + 1);
	VrsutController->nT = (iMaxT - iMinT + 1);

	return 0;
}

/****************************************************************************/
// cme[i].rsutIndex = ( iR*(nS*nU*nT) + iS*(nU*nT) + iU*(nT) + iT )

long calcAllCoulombMatrixElementsRSUTIndices(coulombMatrixElement *cme, long nCME, 
											 long nR, long nS, long nU, long nT) {
	long i;

	// Useful constants
	const long nUT = nU*nT; 
	const long nSUT = nS*nUT;

	// Store the rsutIndex field of all the coulombMatrixElement structures
	for (i = 0; i < nCME; i++) {
		cme[i].rsutIndex = ( cme[i].iR*nSUT + cme[i].iS*nUT + cme[i].iU*nT + cme[i].iT );
	}

	return 0;
}

/****************************************************************************/
//

long calcRSUTIndexToCMEIndexList(long *rsutIndexToCMEIndexList, coulombMatrixElement *cme, long nCME) {
	long i, nMissingIndices = 0;

	// Initialize all rsutIndexToCMEIndexList to -1
	for (i = 0; i < nCME; i++) {
		rsutIndexToCMEIndexList[i] = -1;
	}

	// Store the position index, i, of the cme structure in the list[cme[i].rsutIndex]  
	for (i = 0; i < nCME; i++) {
		if (cme[i].rsutIndex < 0) {
			continue;
		} 
		else {
			rsutIndexToCMEIndexList[cme[i].rsutIndex] = i;
		}
	}

	return 0;
}

/****************************************************************************/
//

long deepCopySingleCoulombMatrixElement(coulombMatrixElement *destCME, long iDestCME, 
										coulombMatrixElement *origCME, long iOrigCME) {
	
	destCME[iDestCME].index = origCME[iOrigCME].index;
	destCME[iDestCME].alreadyCalculatedFlag = origCME[iOrigCME].alreadyCalculatedFlag;
	destCME[iDestCME].value = origCME[iOrigCME].value;
	destCME[iDestCME].psiS  = origCME[iOrigCME].psiS;
	destCME[iDestCME].iU    = origCME[iOrigCME].iU; 
	destCME[iDestCME].psiU  = origCME[iOrigCME].psiU;
	destCME[iDestCME].iR    = origCME[iOrigCME].iR; 
	destCME[iDestCME].psiR  = origCME[iOrigCME].psiR; 
	destCME[iDestCME].iT    = origCME[iOrigCME].iT; 
	destCME[iDestCME].psiT  = origCME[iOrigCME].psiT;
	destCME[iDestCME].iS    = origCME[iOrigCME].iS; 
	destCME[iDestCME].psiS  = origCME[iOrigCME].psiS; 

	return 0;
}

/****************************************************************************/
//

long resetAlreadyCalculatedFlagCME(coulombMatrixElement *coulombMatElement, long nCoulombMatrixElements) {
	long i;

	for (i = 0; i < nCoulombMatrixElements; i++) {
		coulombMatElement[i].alreadyCalculatedFlag = 0;
	}

	return i;
}

/****************************************************************************/
//

long initCoulombMatrixElementStruct(coulombMatrixElement *coulombMatElement, long nCoulombMatrixElements) {
	long i;

	for (i = 0; i < nCoulombMatrixElements; i++) {
		coulombMatElement[i].index = i;
		coulombMatElement[i].alreadyCalculatedFlag = 0;
		coulombMatElement[i].value = 0.0;
		coulombMatElement[i].type = -1;
		coulombMatElement[i].screenedFlag = -1;
		coulombMatElement[i].iR = -1;
		coulombMatElement[i].iS = -1;
		coulombMatElement[i].iU = -1;
		coulombMatElement[i].iT = -1;
	}

	return i;
}


/****************************************************************************/
// 

long calcQSpaceCoulombPotential(zomplex *qSpaceCoulombPot, grid qSpaceGrid, zomplex *rSpaceCoulombPot, grid rSpaceGrid, 
								vector epsilon, long screenedFlag, fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwPsi) { 
	FILE *pf;
	long iGrid, iKx, iKy, iKz, iKyz;
	double dr2, alpha, cosAlpha, sinAlpha, gamma2, sqrk0, minBoxLength, nGridPoints_1;
	zomplex tmp;
	vector pos2;

	// Determine minimum box length and gamma2
	if (rSpaceGrid.maxPos.x < rSpaceGrid.maxPos.y && rSpaceGrid.maxPos.x < rSpaceGrid.maxPos.z) {
		minBoxLength = (rSpaceGrid.maxPos.x - rSpaceGrid.minPos.x);
	}
	else if (rSpaceGrid.maxPos.y < rSpaceGrid.maxPos.x && rSpaceGrid.maxPos.y < rSpaceGrid.maxPos.z) {
		minBoxLength = (rSpaceGrid.maxPos.y - rSpaceGrid.minPos.y);
	}
	else {
		minBoxLength = (rSpaceGrid.maxPos.z - rSpaceGrid.minPos.z);
	}
	gamma2 = sqr(7.0 / minBoxLength);
	sqrk0 = gamma2 * ((epsilon.x + epsilon.y + epsilon.z) / 3.0);

	// FFT the r-space Coulomb potential to get the first part of the qSpaceCoulomb potential
	memcpy(&fftwPsi[0], &rSpaceCoulombPot[0], rSpaceGrid.nGridPoints*sizeof(fftwPsi[0]));
	fftw_execute(planfw[0]);
	memcpy(&qSpaceCoulombPot[0], &fftwPsi[0], rSpaceGrid.nGridPoints*sizeof(qSpaceCoulombPot[0]));	

	// Add the second part of the Coulomb potential to get the final qSpaceCoulomb potential
	nGridPoints_1 = (1.0 / (double)(rSpaceGrid.nGridPoints));
	for (iKz = 0; iKz < qSpaceGrid.nGridPointsZ; iKz++) {
		for (iKy = 0; iKy < qSpaceGrid.nGridPointsY; iKy++) {
			iKyz = qSpaceGrid.nGridPointsX * (qSpaceGrid.nGridPointsY * iKz + iKy);
			for (iKx = 0; iKx < qSpaceGrid.nGridPointsX; iKx++) {
				iGrid = iKyz + iKx;

				alpha = PIE * (double)(iKx + iKy + iKz + qSpaceGrid.nGridPoints / 2);
				cosAlpha = cos(alpha);
				sinAlpha = sin(alpha);

				tmp.re = qSpaceCoulombPot[iGrid].re;
				tmp.im = qSpaceCoulombPot[iGrid].im;

				qSpaceCoulombPot[iGrid].re = (tmp.re * cosAlpha - tmp.im * sinAlpha) * rSpaceGrid.dV;
				qSpaceCoulombPot[iGrid].im = (tmp.re * sinAlpha + tmp.im * cosAlpha) * rSpaceGrid.dV;

				pos2 = retElementWiseVectorMultiplication(qSpaceGrid.gP[iGrid].pos, qSpaceGrid.gP[iGrid].pos);
				if (! screenedFlag) {
					dr2 = (pos2.x + pos2.y + pos2.z + gamma2);
					qSpaceCoulombPot[iGrid].re += (FOURPI / dr2);
				}
				else {
					dr2 = (epsilon.x * pos2.x + epsilon.y * pos2.y + epsilon.z * pos2.z);
					qSpaceCoulombPot[iGrid].re += (FOURPI * (1.0 - exp(-0.25 * dr2 / sqrk0)) / (dr2 + EPSR));
				}
				qSpaceCoulombPot[iGrid].re *= nGridPoints_1;
				qSpaceCoulombPot[iGrid].im *= nGridPoints_1;
			}
		}
	}

	return 0;
}

/****************************************************************************/
// 

long calcRSpaceCoulombPotential(zomplex *rSpaceCoulombPot, grid rSpaceGrid, vector epsilon, long screenedFlag) {
	FILE *pf;
	long iGrid;
	double dr, gamma, gammaEpsilon, sqrtEpsXYZ_1, minBoxLength;
	vector pos2, epsilon_1;

	// Determine minimum box length and gamma
	if (rSpaceGrid.maxPos.x < rSpaceGrid.maxPos.y && rSpaceGrid.maxPos.x < rSpaceGrid.maxPos.z) {
		minBoxLength = (rSpaceGrid.maxPos.x - rSpaceGrid.minPos.x);
	}
	else if (rSpaceGrid.maxPos.y < rSpaceGrid.maxPos.x && rSpaceGrid.maxPos.y < rSpaceGrid.maxPos.z) {
		minBoxLength = (rSpaceGrid.maxPos.y - rSpaceGrid.minPos.y);
	}
	else {
		minBoxLength = (rSpaceGrid.maxPos.z - rSpaceGrid.minPos.z);
	}
	gamma = 7.0 / minBoxLength;
	if (! screenedFlag) {
		fprintf(stdout, "Gamma = %lg\n", gamma); fflush(stdout);
	}

	// Set gammaEpsilon and epsilon_1 depending on if screened or not
	if (! screenedFlag) { // unscreened Coulomb
		gammaEpsilon = gamma;
		epsilon_1.x = epsilon_1.y = epsilon_1.z = 1.0;
		sqrtEpsXYZ_1 = 1.0;
	}
	else { // screened Coulomb
		gammaEpsilon = gamma * sqrt((epsilon.x + epsilon.y + epsilon.z) / 3.0);
		epsilon_1.x = 1.0 / epsilon.x;
		epsilon_1.y = 1.0 / epsilon.y;
		epsilon_1.z = 1.0 / epsilon.z;
		sqrtEpsXYZ_1 = 1.0 / sqrt(epsilon.x * epsilon.y *epsilon.z);
	}

	// Fill in the grid point structures for all grid points
	// The grid points go in the order of z->y->x (i.e. x changes the quickest)
	for (iGrid = 0; iGrid < rSpaceGrid.nGridPoints; iGrid++) {
		rSpaceCoulombPot[iGrid].re = rSpaceCoulombPot[iGrid].im = 0.0;
		dr = rSpaceGrid.gP[iGrid].pos.mag;
		if (dr < minBoxLength) { // used in bse/init.c
			if (screenedFlag) {
				pos2 = retElementWiseVectorMultiplication(rSpaceGrid.gP[iGrid].pos, rSpaceGrid.gP[iGrid].pos);
				dr = sqrt(epsilon_1.x * pos2.x + epsilon_1.y * pos2.y + epsilon_1.z * pos2.z);
			}
			rSpaceCoulombPot[iGrid].re = sqrtEpsXYZ_1 * screenedCoulomb(dr, gammaEpsilon);
		}
	}

	return 0;
}

/****************************************************************************/
//

double screenedCoulomb(double dr, double gamma) {
	if (dr < EPSR) {
		return ( 2.0*gamma / (sqrt(PIE) ));
	}
	else {
		return ( erf(gamma * dr) / dr );
	}
}

/****************************************************************************/
