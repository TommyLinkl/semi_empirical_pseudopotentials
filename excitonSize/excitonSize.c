/*****************************************************************************/
//
// File does the main computations for calculating statistics related to the
// size of excitonic (noninteracting and interacting) states in nanostructures
//
/*****************************************************************************/

#include "es.h"

/****************************************************************************/
// 

long calcExcitonSizeStats(intExc *intExciton, nonintExc *nonintExciton, double *psiHoles, 
							double *psiElecs, grid rSpaceGrid, dParams dPar, lParams lPar) {
	FILE *pf;
	long iGrid;
	vector *Oij, *Oab; // <psi_i|x^2|psi_j> and <psi_a|x^2|psi_b>
	vector *Uij, *Uab; //  <psi_i|x|psi_j>  and  <psi_a|x|psi_b>

	// Allocate memory
	if ((Oij = (vector *) calloc(lPar.nHoles*lPar.nHoles, sizeof(vector))) == NULL) memoryError("Oij");
	if ((Oab = (vector *) calloc(lPar.nElecs*lPar.nElecs, sizeof(vector))) == NULL) memoryError("Oab");
	if ((Uij = (vector *) calloc(lPar.nHoles*lPar.nHoles, sizeof(vector))) == NULL) memoryError("Uij");
	if ((Uab = (vector *) calloc(lPar.nElecs*lPar.nElecs, sizeof(vector))) == NULL) memoryError("Uab");

	// Calculate the matrices required for <psi(r_e, r_h)|r^2|psi(r_e, r_h)> matrix elements
	calcOrsMatrix(Oij, psiHoles, lPar.nHoles, psiHoles, lPar.nHoles, rSpaceGrid);
	calcOrsMatrix(Oab, psiElecs, lPar.nElecs, psiElecs, lPar.nElecs, rSpaceGrid);
	calcUrsMatrix(Uij, psiHoles, lPar.nHoles, psiHoles, lPar.nHoles, rSpaceGrid);
	calcUrsMatrix(Uab, psiElecs, lPar.nElecs, psiElecs, lPar.nElecs, rSpaceGrid);

	// Calculate final exciton size properties for the noninteracting excitons
	calcExcSizeNoninteractingExcitons(nonintExciton, Oij, Oab, Uij, Uab, lPar);

	// Calculate final exciton size properties for the interacting excitons
	if (lPar.intExcitons) {
		calcExcSizeInteractingExcitons(intExciton, Oij, Oab, Uij, Uab, lPar);
	}

	// Free dynamically allocated memory
	free(Uij); free(Uab);
	free(Oij); free(Oab);

	return 0;
}

/****************************************************************************/
//

long calcExcSizeInteractingExcitons(intExc *intExciton, vector *Oij, vector *Oab, vector *Uij, 
										vector *Uab, lParams lPar) {
	FILE *pf;
	long iIntExc, i, j, a, b, iHoleI, iHoleJ, iElecA, iElecB;
	//long   ; // TODO: optimize 
	double Cai, Cbj, CaiCbj;
	vector tmpUijabVector, tmpOijabVector, tmpVector, elecHoleSeparation[lPar.nIntExcitons];

	pf = fopen("sizeIntExcitons.dat", "w");
	fprintf(pf, "#iExc  energy   <rx>    <ry>    <rxy>   <rz>    <r>\n");
	for (iIntExc = 0; iIntExc < lPar.nIntExcitons; iIntExc++) {
		elecHoleSeparation[iIntExc] = retZeroVector();
		for (i = 0; i < lPar.nHoles; i++) {
			for (a = 0; a < lPar.nElecs; a++) {
				Cai = intExciton[iIntExc].Cai[a*lPar.nHoles + (lPar.nHoles-1-i)];
				for (j = 0; j < lPar.nHoles; j++) {
					for (b = 0; b < lPar.nElecs; b++) {
						Cbj = intExciton[iIntExc].Cai[b*lPar.nHoles + (lPar.nHoles-1-j)];
						CaiCbj = Cai*Cbj;
						if (i == j && a == b) {
							tmpOijabVector = retAddedVectors(Oij[i*lPar.nHoles + j], Oab[a*lPar.nElecs + b]);
						}
						else if (i == j) {
							tmpOijabVector = Oab[a*lPar.nElecs + b];
						}
						else if (a == b) {
							tmpOijabVector = Oij[i*lPar.nHoles + j];
						}
						else {
							tmpOijabVector = retZeroVector();
						}
						tmpOijabVector = retScaledVector(tmpOijabVector, CaiCbj);
						tmpUijabVector = retElementWiseVectorMultiplication(Uij[i*lPar.nHoles + j], Uab[a*lPar.nElecs + b]);
						tmpUijabVector = retScaledVector(tmpUijabVector, -2.0*CaiCbj);
						tmpVector = retAddedVectors(tmpOijabVector, tmpUijabVector);
						elecHoleSeparation[iIntExc] = retAddedVectors(elecHoleSeparation[iIntExc], tmpVector);
					}
				}
			}
		}
		// convert from atomic units to nanometers
		elecHoleSeparation[iIntExc] = retScaledVector(elecHoleSeparation[iIntExc], sqr(AUTONM));
		fprintf(pf, "%5d  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n", iIntExc, intExciton[iIntExc].energy*AUTOEV, sqrt(elecHoleSeparation[iIntExc].x), 
					sqrt(elecHoleSeparation[iIntExc].y), sqrt((elecHoleSeparation[iIntExc].x  + elecHoleSeparation[iIntExc].y)),
					sqrt(elecHoleSeparation[iIntExc].z), sqrt(elecHoleSeparation[iIntExc].mag));
		fflush(pf);
	}
	fclose(pf);

	return 0;
}

/****************************************************************************/
//

long calcExcSizeNoninteractingExcitons(nonintExc *nonintExciton, vector *Oij, vector *Oab, vector *Uij, 
										vector *Uab, lParams lPar) {
	FILE *pf;
	long iNonintExc, i, a, iDiagonal, aDiagonal;
	vector tmpVector, elecHoleSeparation[lPar.nNonintExcitons];

	// i == j and a == b for the elec-hole density in noninteracting exciton
	pf = fopen("sizeNonintExcitons.dat", "w");
	fprintf(pf, "#iExc iHole iElec  energy   <rx>    <ry>    <rxy>   <rz>    <r>\n");
	for (iNonintExc = 0; iNonintExc < lPar.nNonintExcitons; iNonintExc++) {
		i = nonintExciton[iNonintExc].h->index;
		iDiagonal = i*lPar.nHoles + i; // used to index Uii and Oii
		a = nonintExciton[iNonintExc].e->index;
		aDiagonal = a*lPar.nElecs + a; // used to index Uaa and Oaa
		tmpVector = retElementWiseVectorMultiplication(Uij[iDiagonal], Uab[aDiagonal]);
		tmpVector = retScaledVector(tmpVector, -2.0);
		// now calculate Oii + Oaa - 2.0*Uii*Uaa
		elecHoleSeparation[iNonintExc] = retAddedVectors(Oij[iDiagonal], Oab[aDiagonal]);
		elecHoleSeparation[iNonintExc] = retAddedVectors(elecHoleSeparation[iNonintExc], tmpVector);
		// convert from atomic units to nanometers
		elecHoleSeparation[iNonintExc] = retScaledVector(elecHoleSeparation[iNonintExc], sqr(AUTONM));
		fprintf(pf, "%5d %5d %5d  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n", iNonintExc, i, a, nonintExciton[iNonintExc].energy*AUTOEV,
					sqrt(elecHoleSeparation[iNonintExc].x), sqrt(elecHoleSeparation[iNonintExc].y), 
					sqrt((elecHoleSeparation[iNonintExc].x + elecHoleSeparation[iNonintExc].y)),
					sqrt(elecHoleSeparation[iNonintExc].z), sqrt(elecHoleSeparation[iNonintExc].mag));
		fflush(pf);
	}
	fclose(pf);


	return 0;
}

/****************************************************************************/
//

long calcUrsMatrix(vector *Urs, double *psiR, long nRStates, double *psiS, long nSStates, grid rSpaceGrid) {
	long r, s, iGrid;
	long rnGridPoints, snGridPoints, rnSStates, snRStates;
	double tmp;
	vector sum;

	// Calculate Urs matrix elements: <psi_r|x|psi_s>, <psi_r|y|psi_s>, <psi_r|z|psi_s>
#pragma omp parallel for private(r, s, sum, tmp, iGrid,  rnGridPoints, snGridPoints, rnSStates, snRStates)
	for (r = 0; r < nRStates; r++) {
		rnGridPoints = r*rSpaceGrid.nGridPoints;
		rnSStates = r*nSStates;
		for (s = r; s < nSStates; s++) {
			snGridPoints = s*rSpaceGrid.nGridPoints;
			snRStates = s*nRStates;
			sum = retZeroVector();
			for (iGrid = 0; iGrid < rSpaceGrid.nGridPoints; iGrid++) {
				tmp = psiR[rnGridPoints + iGrid] * psiS[snGridPoints + iGrid];
				sum = retAddedVectors(sum, retScaledVector(rSpaceGrid.gP[iGrid].pos, tmp));
			}
			Urs[rnSStates + s] = retScaledVector(sum, rSpaceGrid.dV);
			Urs[r + snRStates] = retScaledVector(sum, rSpaceGrid.dV); // since x,y,z can act either way
		}
	}

	return 0;
}

/****************************************************************************/
//

long calcOrsMatrix(vector *Ors, double *psiR, long nRStates, double *psiS, long nSStates, grid rSpaceGrid) {
	long r, s, iGrid;
	long rnGridPoints, snGridPoints, rnSStates, snRStates;
	double tmp;
	vector sum;

	// Calculate Ors matrix elements: <psi_r|x^2|psi_s>, <psi_r|y^2|psi_s>, <psi_r|z^2|psi_s>
#pragma omp parallel for private(r, s, sum, tmp, iGrid,  rnGridPoints, snGridPoints, rnSStates, snRStates)
	for (r = 0; r < nRStates; r++) {
		rnGridPoints = r*rSpaceGrid.nGridPoints;
		rnSStates = r*nSStates;
		for (s = r; s < nSStates; s++) {
			snGridPoints = s*rSpaceGrid.nGridPoints;
			snRStates = s*nRStates;
			sum = retZeroVector();
			for (iGrid = 0; iGrid < rSpaceGrid.nGridPoints; iGrid++) {
				tmp = psiR[rnGridPoints + iGrid] * psiS[snGridPoints + iGrid];
				sum.x += tmp * sqr(rSpaceGrid.gP[iGrid].pos.x);
				sum.y += tmp * sqr(rSpaceGrid.gP[iGrid].pos.y);
				sum.z += tmp * sqr(rSpaceGrid.gP[iGrid].pos.z);
			}
			Ors[rnSStates + s] = retScaledVector(sum, rSpaceGrid.dV);
			Ors[r + snRStates] = retScaledVector(sum, rSpaceGrid.dV); // since x^2,y^2,z^2 can act either way
		}
	}

	return 0;
}

/****************************************************************************/
