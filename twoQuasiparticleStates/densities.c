/*****************************************************************************/
//
// This file calculates carrier density related properties  
//
/*****************************************************************************/

#include "qp.h"

/*****************************************************************************/
//

void calcQPDensitiesFixedOtherQP(intTwoQPState *intTwoQP, double *psiHoles, double *psiElecs, 
									grid rSpaceGrid, dParams dPar, lParams lPar) {
	FILE *pf;
	long iIntState, iNonIntState1, iNonIntState2, iR, iS, iU, iT, iFixedQP, iRiUiFQP, iSiTiFQP;
	long iGrid, nIntTwoQPStatesToPrint, nPrintedStates, nR, nS, nU, nT, nF, nGrid, nFPnIntTwoQPs; 
	long index, iFixedGridPoint[lPar.nFixedQPPoints];
	double *psiR, *psiS, *psiU, *psiT, Crs, Cut, CrsCut;
	double *psiSpsiT, *psiRpsiU, *Cruf, *Cstf, *qp1RhoFixedQP2, *qp2RhoFixedQP1;  
	// double qp1RhoFixedQP2[rSpaceGrid.nGridPoints*lPar.nFixedQPPoints];
	// double qp2RhoFixedQP1[rSpaceGrid.nGridPoints*lPar.nFixedQPPoints]; 
	char fileName[100];
	char qp1Type[100], qp2Type[100];

	// Determine number of excitonic states for which projected densities will be printed
	if (lPar.maxIntTwoQPStatesToPrint > lPar.nIntTwoQPStates) {
		nIntTwoQPStatesToPrint = lPar.nIntTwoQPStates;
	}
	else {
		nIntTwoQPStatesToPrint = lPar.maxIntTwoQPStatesToPrint;
	}
	nF = lPar.nFixedQPPoints;
	nFPnIntTwoQPs = nIntTwoQPStatesToPrint*nF;
	nGrid = rSpaceGrid.nGridPoints;

	// Write beginning of function
	writeSeparation(stdout); writeCurrentTime(stdout);
	fprintf(stdout, "Beginning the calculation of the QP densities with other QP at fixed location\n");
	fprintf(stdout, "Number of densities to be calculated = %ld\n", nFPnIntTwoQPs); 
	fflush(stdout);

	// Determine the number of QP1 and QP2 states
	if (! strcmp(lPar.calcType, "elec-hole")) {
		nR = nU = lPar.nHoles;
		psiR = psiU = psiHoles;
		nS = nT = lPar.nElecs;
		psiS = psiT = psiElecs;
	}
	else if (! strcmp(lPar.calcType, "hole-hole")) {
		nR = nS = nU = nT = lPar.nHoles;
		psiR = psiU = psiS = psiT = psiHoles;
	}
	else { // elec-elec calculation
		nR = nS = nU = nT = lPar.nElecs;
		psiR = psiU = psiS = psiT = psiElecs;
	}

	// Dynamically allocate memory
	if ((psiRpsiU = (double *) calloc(nR*nU*nF, sizeof(double))) == NULL) memoryError("psiRpsiU");
	if ((psiSpsiT = (double *) calloc(nS*nT*nF, sizeof(double))) == NULL) memoryError("psiSpsiT");
	if ((Cruf = (double *) calloc(nR*nU*nFPnIntTwoQPs, sizeof(double))) == NULL) memoryError("Cruf");
	if ((Cstf = (double *) calloc(nS*nT*nFPnIntTwoQPs, sizeof(double))) == NULL) memoryError("Cstf");
	if ((qp1RhoFixedQP2 = (double *) calloc(nF*nGrid, sizeof(double))) == NULL) memoryError("qp1RhoFixedQP2");
	if ((qp2RhoFixedQP1 = (double *) calloc(nF*nGrid, sizeof(double))) == NULL) memoryError("qp2RhoFixedQP1");

	// Find and store the index of the nearest grid point for each fixedQPPoint
	for (iFixedQP = 0; iFixedQP < nF; iFixedQP++) {
		fprintf(stdout, "Fixed point %ld = %.4f %.4f %.4f\n", iFixedQP, dPar.fixedQPPoints[iFixedQP].x, 
										dPar.fixedQPPoints[iFixedQP].y, dPar.fixedQPPoints[iFixedQP].z);
		iFixedGridPoint[iFixedQP] = retIndexOfClosestGridPoint(rSpaceGrid, dPar.fixedQPPoints[iFixedQP]);
		fprintf(stdout, "Closest grid point = %.4f %.4f %.4f\n", rSpaceGrid.gP[iFixedGridPoint[iFixedQP]].pos.x, 
							rSpaceGrid.gP[iFixedGridPoint[iFixedQP]].pos.y, 
							rSpaceGrid.gP[iFixedGridPoint[iFixedQP]].pos.z); fflush(stdout);
	}

	// Calculate psiS[iFP]*psiT[iFP] and psiR[iFP]*psiU[iFP] for all fixed points
	for (iS = 0; iS < nS; iS++) {
		for (iT = 0; iT < nT; iT++) {
			for (iFixedQP = 0; iFixedQP < nF; iFixedQP++) {
				index = iS*(nT+nF) + iT*nF + iFixedQP;
				iGrid = iFixedGridPoint[iFixedQP];
				psiSpsiT[index] = psiS[iS*nGrid + iGrid]*psiT[iT*nGrid + iGrid];
			}
		}
	}
	// only needed for elec-hole since psiSpsiT=psiRpsiU for hole-hole and elec-elec calculations
	//if (! strcmp(lPar.calcType, "elec-hole")) { 
	for (iR = 0; iR < nR; iR++) {
		for (iU = 0; iU < nU; iU++) {
			for (iFixedQP = 0; iFixedQP < nF; iFixedQP++) {
				index = iR*(nU+nF) + iU*nF + iFixedQP;
				iGrid = iFixedGridPoint[iFixedQP];
				psiRpsiU[index] = psiR[iR*nGrid + iGrid]*psiU[iU*nGrid + iGrid];
			}
		}
	}
	//}

	// Main computational part of the program
	nPrintedStates = 0;
	for (iIntState = 0; iIntState < lPar.nIntTwoQPStates && nPrintedStates < nIntTwoQPStatesToPrint; iIntState++) {
		if (intTwoQP[iIntState].energy < dPar.minIntTwoQPEnergyToPrint) {
			continue;
		}
		nPrintedStates++;
		// Calculate Cruf and Cstf matrices
		zeroDoubleArray(Cruf, nR*nU*nF);
		zeroDoubleArray(Cstf, nS*nT*nF);
		// TODO: have separate loops
		for (iNonIntState1 = 0; iNonIntState1 < lPar.nNonintTwoQPStates; iNonIntState1++) {
			iR = intTwoQP[iIntState].niTwoQP[iNonIntState1].qp1->index;
			psiR = intTwoQP[iIntState].niTwoQP[iNonIntState1].qp1->psi; // hole for elec-hole and hole-hole
			iS = intTwoQP[iIntState].niTwoQP[iNonIntState1].qp2->index;
			psiS = intTwoQP[iIntState].niTwoQP[iNonIntState1].qp2->psi; // elec for elec-hole and elec-elec
			Crs = intTwoQP[iIntState].Crs[iNonIntState1];
			for (iNonIntState2 = 0; iNonIntState2 < lPar.nNonintTwoQPStates; iNonIntState2++) {
				iU = intTwoQP[iIntState].niTwoQP[iNonIntState2].qp1->index;
				psiU = intTwoQP[iIntState].niTwoQP[iNonIntState2].qp1->psi; // hole for elec-hole and hole-hole
				iT = intTwoQP[iIntState].niTwoQP[iNonIntState2].qp2->index;
				psiT = intTwoQP[iIntState].niTwoQP[iNonIntState2].qp2->psi; // elec for elec-hole and elec-elec
				Cut = intTwoQP[iIntState].Crs[iNonIntState2];
				CrsCut = Crs*Cut;
				for (iFixedQP = 0; iFixedQP < nF; iFixedQP++) {
					index = iS*(nT+nF) + iT*nF + iFixedQP;
					Cruf[iR*(nU+nF) + iU*nF + iFixedQP] += CrsCut*psiSpsiT[index];
				}
				for (iFixedQP = 0; iFixedQP < nF; iFixedQP++) {
					index = iR*(nU+nF) + iU*nF + iFixedQP;
					Cstf[iS*(nT+nF) + iT*nF + iFixedQP] += CrsCut*psiRpsiU[index];
				}
			}
		}
		// Calculate qpRho
		zeroDoubleArray(qp1RhoFixedQP2, nGrid*nF);
		zeroDoubleArray(qp2RhoFixedQP1, nGrid*nF);
		// TEMP 
		psiR = psiU = psiHoles;
		psiS = psiT = psiElecs;
		// END TEMP
		for (iR = 0; iR < nR; iR++) {
			for (iU = 0; iU < nU; iU++) {
				for (iFixedQP = 0; iFixedQP < nF; iFixedQP++) {
					index = iR*(nU+nF) + iU*nF + iFixedQP;
					for (iGrid = 0; iGrid < nGrid; iGrid++) {
						qp1RhoFixedQP2[iFixedQP*nGrid + iGrid] += (Cruf[index]*
																		psiR[iR*nGrid + iGrid]*
																		psiU[iU*nGrid + iGrid]);
					}
				}
			}
		}
		for (iS = 0; iS < nS; iS++) {
			for (iT = 0; iT < nT; iT++) {
				for (iFixedQP = 0; iFixedQP < nF; iFixedQP++) {
					index = iS*(nT+nF) + iT*nF + iFixedQP;
					for (iGrid = 0; iGrid < nGrid; iGrid++) {
						qp2RhoFixedQP1[iFixedQP*nGrid + iGrid] += (Cstf[index]*
																		psiS[iR*nGrid + iGrid]*
																		psiT[iU*nGrid + iGrid]);
					}
				}
			}
		}
		// Normalize distributions
		for (iFixedQP = 0; iFixedQP < nF; iFixedQP++) {
			normalizeDoubleArray(&(qp1RhoFixedQP2[iFixedQP*nGrid]), nGrid);
			normalizeDoubleArray(&(qp2RhoFixedQP1[iFixedQP*nGrid]), nGrid);
		}
		// Print cube files and projected densities
		if (! strcmp(lPar.calcType, "elec-hole") || ! strcmp(lPar.calcType, "hole-hole")) {
			sprintf(qp1Type, "Hole-i-%ld", iIntState);
			sprintf(qp2Type, "Elec-i-%ld", iIntState);
		}
		else { // calcType == "elec-elec"
			sprintf(qp1Type, "Elec-i-%ld", iIntState);
			sprintf(qp2Type, "Elec2-i-%ld", iIntState); // TODO: remove eventually
		}
		for (iFixedQP = 0; iFixedQP < nF; iFixedQP++) {
			sprintf(fileName, "fp%ld-rho%s%s", iFixedQP, qp1Type, ".cube");
			writeCubeFile(&(qp1RhoFixedQP2[iFixedQP*nGrid]), rSpaceGrid, fileName);
			// Calculate and write x, y and z projected probability densities
			sprintf(fileName, "fp%ld-projRho%s%s", iFixedQP, qp1Type, ".dat");
			calcXYZProjectedProbDensities(&(qp1RhoFixedQP2[iFixedQP*nGrid]), rSpaceGrid, fileName);
			//// For elec-hole need to print projected density for second (elec) quasiparticle too
			//if (! strcmp(lPar.calcType, "elec-hole")) {
			sprintf(fileName, "fp%ld-rho%s%s", iFixedQP, qp2Type, ".cube");
			writeCubeFile(&(qp2RhoFixedQP1[iFixedQP*nGrid]), rSpaceGrid, fileName);
			sprintf(fileName, "fp%ld-projRho%s%s", iFixedQP, qp2Type, ".dat");
			calcXYZProjectedProbDensities(&(qp2RhoFixedQP1[iFixedQP*nGrid]), rSpaceGrid, fileName);
			//}
		}
	}
	
	// Free dynamically allocated memory
	free(qp1RhoFixedQP2); free(qp2RhoFixedQP1);
	free(Cruf); free(Cstf);
	free(psiSpsiT); free(psiRpsiU);

	// Write ending of function
	fprintf(stdout, "\nFinished  the calculation of the QP densities with other QP at fixed location\n");
	writeCurrentTime(stdout); writeSeparation(stdout); fflush(stdout);

	return;
}

/*****************************************************************************/
// amplitude for e-h = psiR(xH)*psiS(xE) = psiH(xH)*psiE(xE)
// amplitude for h-h = psiR(xH1)*psiS(xH2)-psiR(xH2)*psiS(xH1) = psiH(xH1)*psiH(xH2)-psiH(xH2)*psiH(xH1)
// amplitude for e-e = psiR(xE1)*psiS(xE2)-psiR(xE2)*psiS(xE1) = psiE(xE1)*psiE(xE2)-psiE(xE2)*psiE(xE1)

void calcIntTwoQPProjDensities(intTwoQPState *intTwoQP, grid rSpaceGrid, dParams dPar, lParams lPar) {
        
	FILE *pf;
	long iIntState, iNonIntState1, iNonIntState2, iR, iS, iU, iT;
	long iGrid, nIntTwoQPStatesToPrint, nPrintedStates; 
	double *psiR, *psiS, *psiU, *psiT, Crs, Cut, CrsCut;
	// double rhoQP1[rSpaceGrid.nGridPoints], rhoQP2[rSpaceGrid.nGridPoints];
        double *rhoQP1, *rhoQP2;
	char fileName[100];
	char qp1Type[100], qp2Type[100];

	// Determine number of excitonic states for which projected densities will be printed
	if (lPar.maxIntTwoQPStatesToPrint > lPar.nIntTwoQPStates) {
		nIntTwoQPStatesToPrint = lPar.nIntTwoQPStates;
	}
	else {
		nIntTwoQPStatesToPrint = lPar.maxIntTwoQPStatesToPrint;
	}

	// Write beginning of function
	writeSeparation(stdout); writeCurrentTime(stdout);
	fprintf(stdout, "Beginning the calculation of the interacting two QP state projected densities\n");
	fprintf(stdout, "Maximum number of states to be calculated = %ld\n", nIntTwoQPStatesToPrint);
	fprintf(stdout, "Mininum energy of state  to be calculated = %.8f\n", dPar.minIntTwoQPEnergyToPrint); fflush(stdout);

        // added by Dipti
        // dynamically allocate memory
        if ((rhoQP1 = (double *) calloc(rSpaceGrid.nGridPoints, sizeof(double))) == NULL) memoryError("rhoQP1");
        if ((rhoQP2 = (double *) calloc(rSpaceGrid.nGridPoints, sizeof(double))) == NULL) memoryError("rhoQP2");

	// Main computational part of program
	nPrintedStates = 0;
	for (iIntState = 0; iIntState < lPar.nIntTwoQPStates && nPrintedStates < nIntTwoQPStatesToPrint; iIntState++) {
		if (intTwoQP[iIntState].energy < dPar.minIntTwoQPEnergyToPrint) {
			continue;
		}
		nPrintedStates++;
		for (iGrid = 0; iGrid < rSpaceGrid.nGridPoints; iGrid++) {
			rhoQP1[iGrid] = 0.0;
			rhoQP2[iGrid] = 0.0;
		}
		for (iNonIntState1 = 0; iNonIntState1 < lPar.nNonintTwoQPStates; iNonIntState1++) {
			iR = intTwoQP[iIntState].niTwoQP[iNonIntState1].qp1->index;
			psiR = intTwoQP[iIntState].niTwoQP[iNonIntState1].qp1->psi; // hole for elec-hole and hole-hole
			iS = intTwoQP[iIntState].niTwoQP[iNonIntState1].qp2->index;
			psiS = intTwoQP[iIntState].niTwoQP[iNonIntState1].qp2->psi; // elec for elec-hole and elec-elec
			Crs = intTwoQP[iIntState].Crs[iNonIntState1];
			for (iNonIntState2 = 0; iNonIntState2 < lPar.nNonintTwoQPStates; iNonIntState2++) {
				iU = intTwoQP[iIntState].niTwoQP[iNonIntState2].qp1->index;
				psiU = intTwoQP[iIntState].niTwoQP[iNonIntState2].qp1->psi; // hole for elec-hole and hole-hole
				iT = intTwoQP[iIntState].niTwoQP[iNonIntState2].qp2->index;
				psiT = intTwoQP[iIntState].niTwoQP[iNonIntState2].qp2->psi; // elec for elec-hole and elec-elec
				Cut = intTwoQP[iIntState].Crs[iNonIntState2];
				// TODO: edit this to work with elec-elec and hole-hole
				if (iS == iT) { // delta_ST = integral(psiS*psiT)
					// Calculate rhoQP1(r) += Crs*Cus*psiR(r)*psiU(r)
					CrsCut = Crs*Cut;
				#pragma omp parallel for private(iGrid)
					for (iGrid = 0; iGrid < rSpaceGrid.nGridPoints; iGrid++) {
						rhoQP1[iGrid] += CrsCut*psiR[iGrid]*psiU[iGrid];
					}
				}
				if (iR == iU && ! strcmp(lPar.calcType, "elec-hole")) { // delta_RU = integral(psiR*psiI)
					// Calculate rhoQP2(r) += Crs*Crt*psiS(r)*psiT(r)
					CrsCut = Crs*Cut;
				#pragma omp parallel for private(iGrid)
					for (iGrid = 0; iGrid < rSpaceGrid.nGridPoints; iGrid++) {
						rhoQP2[iGrid] += CrsCut*psiS[iGrid]*psiT[iGrid];
					}
				}
			}
		}
		// Normalize rhoQP1 and rhoQP2
		for (iGrid = 0; iGrid < rSpaceGrid.nGridPoints; iGrid++) {
			rhoQP1[iGrid] *= rSpaceGrid.dV;
			rhoQP2[iGrid] *= rSpaceGrid.dV;
			// TODO: have functions calls to normalizeDoubleArray
		}
		// Print rhoQP1 and rhoQP2
		if (! strcmp(lPar.calcType, "elec-hole") || ! strcmp(lPar.calcType, "hole-hole")) {
			sprintf(qp1Type, "Hole-i-%ld", iIntState);
			sprintf(qp2Type, "Elec-i-%ld", iIntState);
		}
		else { // calcType == "elec-elec"
			sprintf(qp1Type, "Elec-i-%ld", iIntState);
		}
		sprintf(fileName, "rho%s%s", qp1Type, ".cube");
		writeCubeFile(rhoQP1, rSpaceGrid, fileName);
		// Calculate and write x, y and z projected probability densities
		sprintf(fileName, "projRho%s%s", qp1Type, ".dat");
		calcXYZProjectedProbDensities(rhoQP1, rSpaceGrid, fileName);
		// For elec-hole need to print projected density for second (elec) quasiparticle too
		if (! strcmp(lPar.calcType, "elec-hole")) {
			sprintf(fileName, "rho%s%s", qp2Type, ".cube");
			writeCubeFile(rhoQP2, rSpaceGrid, fileName);
			sprintf(fileName, "projRho%s%s", qp2Type, ".dat");
			calcXYZProjectedProbDensities(rhoQP2, rSpaceGrid, fileName);
		}
	}

        // Free dynamically allocated memory
        free(rhoQP1); free(rhoQP2);

	// Write ending of function
	fprintf(stdout, "\nFinished  the calculation of the interacting two QP state projected densities\n");
	writeCurrentTime(stdout);
	writeSeparation(stdout); fflush(stdout);

	return;
}

/*****************************************************************************/
//

void calcNonintQPDensities(nonintQP *holeQP, nonintQP *elecQP, grid rSpaceGrid, dParams dPar, lParams lPar) {
        
	FILE *pf;
	long iQP, iElec, iHole, nPrintedElecs, nPrintedHoles, printQPStateFlag;
	double *psiQP;
	// double rho[rSpaceGrid.nGridPoints];
        double *rho;
	char fileName[100];
	char qpType[100];

	// Write beginning of function
	writeSeparation(stdout);
	writeCurrentTime(stdout);
	fprintf(stdout, "Beginning the calculation of the noninteracting QP densities\n"); 
	fprintf(stdout, "Maximum number of hole densities to be calculated = %ld\n", lPar.nHoles);
	fprintf(stdout, "Maximum energy of hole state to be calculated = % .8f\n", dPar.maxHoleEnergyToPrint);
	fprintf(stdout, "Maximum number of elec densities to be calculated = %ld\n", lPar.nElecs); 
	fprintf(stdout, "Minimum energy of elec state to be calculated = % .8f\n", dPar.minElecEnergyToPrint); fflush(stdout);

        // Dynamically allocate memory
        if ((rho = (double *) calloc(rSpaceGrid.nGridPoints, sizeof(double))) == NULL) memoryError("rho");

	// Main computational part of function
	iElec = iHole = nPrintedElecs = nPrintedHoles = 0;
	for (iQP = 0; iQP < lPar.nHolesPlusElecs; iQP++) {
		printQPStateFlag = 0;
		if (iQP < lPar.nHoles && nPrintedHoles < lPar.maxHoleStatesToPrint) {
			sprintf(qpType, "Hole-ni-%ld", iHole);
			psiQP = holeQP[iHole].psi;
			if (holeQP[iHole].energy < dPar.maxHoleEnergyToPrint) {
				printQPStateFlag = 1;
				nPrintedHoles++;
			}
			else {
				printQPStateFlag = 0;
			}
			iHole++;
		}
		else if (iQP >= lPar.nHoles && nPrintedElecs < lPar.maxElecStatesToPrint) {
			sprintf(qpType, "Elec-ni-%ld", iElec);
			psiQP = elecQP[iElec].psi;
			if (elecQP[iElec].energy > dPar.minElecEnergyToPrint) {
				printQPStateFlag = 1;
				nPrintedElecs++;
			}
			else {
				printQPStateFlag = 0;
			}
			iElec++;
		}
		if (! printQPStateFlag) {
			continue;
		}
		else {
		 	// Write the wavefunction cube file
		 	sprintf(fileName, "psi%s%s", qpType, ".cube");
		 	writeCubeFile(psiQP, rSpaceGrid, fileName);
		 	
		 	// Calculate and write probability density file
		 	calcProbDensity(rho, psiQP, rSpaceGrid.nGridPoints, rSpaceGrid.dV);
		 	sprintf(fileName, "rho%s%s", qpType, ".cube");
			writeCubeFile(rho, rSpaceGrid, fileName);
			
			// Calculate and write x, y and z projected probability densities
			sprintf(fileName, "projRho%s%s", qpType, ".dat");
			calcXYZProjectedProbDensities(rho, rSpaceGrid, fileName);
		}
	}

        // Free dynamically allocated memory
        free(rho);

	// Write ending of function
	fprintf(stdout, "\nFinished  the calculation of the noninteracting QP densities\n");
	writeCurrentTime(stdout); fflush(stdout);
	
	if (! lPar.calcIntDensitiesOnly) {
		writeSeparation(stdout); fflush(stdout);
	}

	if (lPar.calcNonintDensitiesOnly) {
		writeSeparation(stdout); fflush(stdout);
		exit(EXIT_SUCCESS);
	}

	return;
}

/*****************************************************************************/
// Calculates the projected probability densities of rho onto the x, y and z
// assumes rho is ordered: z,y,x (i.e. x is faster changing index)

void calcXYZProjectedProbDensities(double *rho, grid rSpaceGrid, char *fileName) {
	FILE *pf;
	long i, iGrid, iX, iY, iZ, iYZ, nMaxGridPoints;
	double *x, *y, *z;
	double *rhoX, *rhoY, *rhoZ;

	// Determine largest grid point dimension
	if (rSpaceGrid.nGridPointsZ >= rSpaceGrid.nGridPointsX 
								&& rSpaceGrid.nGridPointsZ >= rSpaceGrid.nGridPointsY) {
		nMaxGridPoints = rSpaceGrid.nGridPointsZ;
	}
	else if (rSpaceGrid.nGridPointsY >= rSpaceGrid.nGridPointsX) {
		nMaxGridPoints = rSpaceGrid.nGridPointsY;
	}
	else {
		nMaxGridPoints = rSpaceGrid.nGridPointsX;
	}

	// Dynamically allocate memory for projected probability densities
	if ((rhoX = (double *) calloc(nMaxGridPoints, sizeof(double))) == NULL) memoryError("rhoX");
	if ((x = (double *) calloc(nMaxGridPoints, sizeof(double))) == NULL) memoryError("x");
	if ((rhoY = (double *) calloc(nMaxGridPoints, sizeof(double))) == NULL) memoryError("rhoY");
	if ((y = (double *) calloc(nMaxGridPoints, sizeof(double))) == NULL) memoryError("y");
	if ((rhoZ = (double *) calloc(nMaxGridPoints, sizeof(double))) == NULL) memoryError("rhoZ");
	if ((z = (double *) calloc(nMaxGridPoints, sizeof(double))) == NULL) memoryError("z");

	// Calculate rhoX, rhoY and rhoZ
	for (iZ = 0; iZ < rSpaceGrid.nGridPointsZ; iZ++) {
		for (iY = 0; iY < rSpaceGrid.nGridPointsY; iY++) {
			iYZ = rSpaceGrid.nGridPointsX * (rSpaceGrid.nGridPointsY * iZ + iY);
			for (iX = 0; iX < rSpaceGrid.nGridPointsX; iX++) {
				iGrid = iYZ + iX;
        		rhoX[iX] += rho[iGrid];
        		rhoY[iY] += rho[iGrid];
        		rhoZ[iZ] += rho[iGrid];
  			}
  		}
  		rhoZ[iZ] /= rSpaceGrid.stepSize.z;
  	}
  	for (iY = 0; iY < rSpaceGrid.nGridPointsY; iY++) { rhoY[iY] /= rSpaceGrid.stepSize.y; }
  	for (iX = 0; iX < rSpaceGrid.nGridPointsX; iX++) { rhoX[iX] /= rSpaceGrid.stepSize.x; }

  	// Store x, y and z positions for printing 
  	for (iZ = 0; iZ < rSpaceGrid.nGridPointsZ; iZ++) { 
  		z[iZ] = rSpaceGrid.gP[iZ*rSpaceGrid.nGridPointsY*rSpaceGrid.nGridPointsX].pos.z;
  	}
  	for (iY = 0; iY < rSpaceGrid.nGridPointsY; iY++) { 
  		y[iY] = rSpaceGrid.gP[iY*rSpaceGrid.nGridPointsX].pos.y; 
  	}
  	for (iX = 0; iX < rSpaceGrid.nGridPointsX; iX++) { 
  		x[iX] = rSpaceGrid.gP[iX].pos.x; 
  	}

	// Write rhoX, rhoY and rhoZ
	pf = fopen(fileName, "w");
	fprintf(pf, "#   x        rhoX        y         rhoY         z        rhoZ\n");
	for (i = 0; i < nMaxGridPoints; i++) {
		fprintf(pf, "%8.3f  %.8f  %8.3f  %.8f  %8.3f  %.8f\n", 
					x[i], rhoX[i], y[i], rhoY[i], z[i], rhoZ[i]);
	}
	fclose(pf);

	// Free dynamically allocated memory
	free(z); free(y); free(x);
	free(rhoZ); free(rhoY); free(rhoX);

	return;
}


/*****************************************************************************/
// Calculates rho[iGrid] = psi[iGrid]*psi[iGrid]*dV 

void calcProbDensity(double *rho, double *psi, long nGridPoints, double dV) {
	long iGrid;

	// Calculate |psi[iGrid]|^2
	calcPsiSquared(rho, psi, nGridPoints);

	// Multiple by volume element to get probability density
	for (iGrid = 0; iGrid < nGridPoints; iGrid++) {
		rho[iGrid] *= dV;
	}

	return;
}

/*****************************************************************************/
// Calculates psiSquared[iGrid] = psi[iGrid]*psi[iGrid] 

void calcPsiSquared(double *psiSquared, double *psi, long nGridPoints) {
	long iGrid;

	for (iGrid = 0; iGrid < nGridPoints; iGrid++) {
		psiSquared[iGrid] = psi[iGrid]*psi[iGrid];
	}

	return;
}

/*****************************************************************************/
// Zeros a double array

void zeroDoubleArray(double *dArray, long nElements) {
	long i;

	for (i = 0; i < nElements; i++) {
		dArray[i] = 0.0;
	}

	return;
}

/*****************************************************************************/
