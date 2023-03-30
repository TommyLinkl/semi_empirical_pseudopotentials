/*****************************************************************************/
//
// This file deals with the Hamiltonian for correlated two quasiparticle states
//
/*****************************************************************************/

#include "qp.h"

/****************************************************************************/
//

long calcTwoQPHamiltonianWrapper(double *hMatrix, double *h0Matrix, double *Wrsut, double *Vrsut, 
						  nonintTwoQPState *nonintTwoQP, lParams lPar) {

	// Call correct two QP Hamiltonian function depending on the calcType
	if (! strcmp(lPar.calcType, "sp-elec-hole")) {
		calcSpinPolarizedTwoQPHamiltonian(hMatrix, h0Matrix, Wrsut, Vrsut, nonintTwoQP, lPar);
	}
	else {
		calcTwoQPHamiltonian(hMatrix, h0Matrix, Wrsut, Vrsut, nonintTwoQP, lPar);
	}

	return 0;
}

/****************************************************************************/
//

long calcSpinPolarizedTwoQPHamiltonian(double *hMatrix, double *h0Matrix, double *Wrsut, double *Vrsut, 
										nonintTwoQPState *nonintTwoQP, lParams lPar) {
	FILE *pf;
	long iNonintTwoQP1, iNonintTwoQP2, index;
	double deltaSpin, coulombSign;

	// Coulomb term is attractive and exchange is repulsive for sp-elec-hole calculations 
	coulombSign = -1.0; 
	deltaSpin = 1.0; // 1.0 since working in non-diagonal spin basis (2000 PRB, Louie, BSE Theory paper)

	// Calculate the Hamiltonian matrix (hMatrix)
	// h = single-particle term + direct Coulomb term - indirect(exchange-like) Coulomb term
	// For elec-hole: hMatrix = h0Matrix - Wrsut + 2.0*Vrsut*delta_sigma_sigmaPrime
	// qp1 = holeQP -> has spinZ that is the negative of that in the valence band electrons
	// qp2 = elecQP -> has spinZ that is equal to that of the conduction band electrons
	for (iNonintTwoQP1 = 0; iNonintTwoQP1 < lPar.nNonintTwoQPStates; iNonintTwoQP1++) {
		for (iNonintTwoQP2 = 0; iNonintTwoQP2 < lPar.nNonintTwoQPStates; iNonintTwoQP2++) {
			index = iNonintTwoQP1*lPar.nNonintTwoQPStates + iNonintTwoQP2;
			// 2 of the 3 spin triplet eigenstates: |1,1> and |1,-1>:
			// |0.5,0.5>|0.5,0.5> or |-0.5,-0.5>|-0.5,-0.5>
			if ( (fabs(nonintTwoQP[iNonintTwoQP1].spin - 1.0) < EPS) && 
				 (fabs(nonintTwoQP[iNonintTwoQP2].spin - 1.0) < EPS) &&
			     (fabs(nonintTwoQP[iNonintTwoQP1].spinZ - nonintTwoQP[iNonintTwoQP2].spinZ) < EPS) ) {
				hMatrix[index] = ( h0Matrix[index] + coulombSign*Wrsut[index] );
			}
			// diagonal part of |1,0> and |0,0> spin eigenstates Hilbert space:
			// |-0.5,0.5>|-0.5,0.5> or |0.5,-0.5>|0.5,-0.5>
			else if ( ( (fabs(nonintTwoQP[iNonintTwoQP1].spinZ) < EPS) && (fabs(nonintTwoQP[iNonintTwoQP2].spinZ) < EPS) ) && 
					    (fabs(nonintTwoQP[iNonintTwoQP1].qp1->spinZ - nonintTwoQP[iNonintTwoQP2].qp1->spinZ) < EPS) &&
					    (fabs(nonintTwoQP[iNonintTwoQP1].qp2->spinZ - nonintTwoQP[iNonintTwoQP2].qp2->spinZ) < EPS) ) { 
				hMatrix[index] = (h0Matrix[index] + coulombSign*(Wrsut[index] - deltaSpin*Vrsut[index]));
			}
			// off diagonal part of |1,0> and |0,0> spin eigenstates Hilbert space:
			// |-0.5,0.5>|0.5,-0.5> or |0.5,-0.5>|-0.5,0.5>
			else if ( ( (fabs(nonintTwoQP[iNonintTwoQP1].spinZ) < EPS) && (fabs(nonintTwoQP[iNonintTwoQP2].spinZ) < EPS) ) && 
					    (fabs(nonintTwoQP[iNonintTwoQP1].qp1->spinZ - nonintTwoQP[iNonintTwoQP2].qp2->spinZ) < EPS) &&
					    (fabs(nonintTwoQP[iNonintTwoQP1].qp2->spinZ - nonintTwoQP[iNonintTwoQP2].qp1->spinZ) < EPS) ) { 	
				hMatrix[index] = ( deltaSpin*Vrsut[index] );
			}
			else { // J != J' for |J,M>|J',M'>
				hMatrix[index] = 0.0;
			}
		}
	}

	return 0;
}

/****************************************************************************/
//

long calcTwoQPHamiltonian(double *hMatrix, double *h0Matrix, double *Wrsut, double *Vrsut, 
						  nonintTwoQPState *nonintTwoQP, lParams lPar) {
	FILE *pf;
	long iNonintTwoQP1, iNonintTwoQP2, index;
	double deltaSpin, coulombSign;

	// Determine is spin singlets or triplets are being calculated
	if (lPar.calcSinglets && ! strcmp(lPar.calcType, "elec-hole")) {
		deltaSpin = 2.0; 
	}
	else if (lPar.calcTriplets && (! strcmp(lPar.calcType, "elec-elec") || ! strcmp(lPar.calcType, "hole-hole"))) {
		deltaSpin = 1.0;
	}
	else {
		deltaSpin = 0.0;
	}

	// Determine is the Coulomb interaction is attractive or repulsive
	if (! strcmp(lPar.calcType, "elec-hole")) {
		coulombSign = -1.0; // elec-hole Coulomb interaction is attractive
	}
	else {
		coulombSign = 1.0;  // elec-elec & hole-hole Coulomb interaction is repulsive
	}

	// Calculate the Hamiltonian matrix (hMatrix)
	// h = single-particle term + direct Coulomb term - indirect(exchange-like) Coulomb term
	// For elec-hole: hMatrix = h0Matrix - Wrsut + 2.0*Vrsut*delta_sigma_sigmaPrime
	// For elec-elec: hMatrix = h0Matrix + Wrsut - 1.0*Wusrt*delta_sigma_sigmaPrime
	// For hole-hole: hMatrix = h0Matrix + Wrsut - 1.0*Wusrt*delta_sigma_sigmaPrime  
	for (iNonintTwoQP1 = 0; iNonintTwoQP1 < lPar.nNonintTwoQPStates; iNonintTwoQP1++) {
		for (iNonintTwoQP2 = 0; iNonintTwoQP2 < lPar.nNonintTwoQPStates; iNonintTwoQP2++) {
			index = iNonintTwoQP1*lPar.nNonintTwoQPStates + iNonintTwoQP2;
			hMatrix[index] = (h0Matrix[index] + coulombSign*(Wrsut[index] - deltaSpin*Vrsut[index]));
		}
	}

	// Print Hamiltonian matrix
	pf = fopen("Hmatrix.dat", "w");
	for (iNonintTwoQP1 = 0; iNonintTwoQP1 < lPar.nNonintTwoQPStates; iNonintTwoQP1++) {
		for (iNonintTwoQP2 = 0; iNonintTwoQP2 < lPar.nNonintTwoQPStates; iNonintTwoQP2++) {
			index = iNonintTwoQP1*lPar.nNonintTwoQPStates + iNonintTwoQP2;
			fprintf(pf, "% .16f ", hMatrix[index]);
		}
		fprintf(pf, "\n");
	}	
	fclose(pf);

	return 0;
}

/****************************************************************************/
//

long calcTwoQPH0Matrix(double *h0Matrix, nonintTwoQPState *nonintTwoQP, long nNonintTwoQPStates) {
	FILE *pf;
	long iNonintTwoQP1, iNonintTwoQP2, index;

	// Calculate h0Matrix
	for (iNonintTwoQP1 = 0; iNonintTwoQP1 < nNonintTwoQPStates; iNonintTwoQP1++) {
		for (iNonintTwoQP2 = 0; iNonintTwoQP2 < nNonintTwoQPStates; iNonintTwoQP2++) {
			index = iNonintTwoQP1*nNonintTwoQPStates + iNonintTwoQP2;
			if (iNonintTwoQP1 == iNonintTwoQP2) { // diagonal only
				h0Matrix[index] = nonintTwoQP[iNonintTwoQP1].energy;
			}
			else { // off-diagonals are zero for h0
				h0Matrix[index] = 0.0;
			}
		}
	}	

	// Print h0 matrix
	pf = fopen("h0matrix.dat", "w");
	for (iNonintTwoQP1 = 0; iNonintTwoQP1 < nNonintTwoQPStates; iNonintTwoQP1++) {
		for (iNonintTwoQP2 = 0; iNonintTwoQP2 < nNonintTwoQPStates; iNonintTwoQP2++) {
			index = iNonintTwoQP1*nNonintTwoQPStates + iNonintTwoQP2;
			fprintf(pf, "%.8f ", h0Matrix[index]);
		}
		fprintf(pf, "\n");
	}	
	fclose(pf);

	return 0;
}

/****************************************************************************/
