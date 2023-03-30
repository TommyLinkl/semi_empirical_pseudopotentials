/*****************************************************************************/
//
// This file deals with normalizing arrays of various types  
//
/*****************************************************************************/

#include "qp.h"

/*****************************************************************************/
// 

void normalizeAllDoubleWavefunctions(double *dPsi, long nBasisStates, double dV, long nPsis) {
	long iPsi;

	for (iPsi = 0; iPsi < nPsis; iPsi++){
		normalizeDoubleWavefunction(&(dPsi[iPsi*nBasisStates]), nBasisStates, dV);
	}

	return;
}

/*****************************************************************************/
// 

void normalizeDoubleWavefunction(double *dPsi, long nBasisStates, double dV) {
	long i;
	double initNorm, oneOverSqrtInitNorm;

	initNorm = retNormOfDoubleWavefunction(dPsi, nBasisStates, dV);
	oneOverSqrtInitNorm = 1.0/sqrt(initNorm);
	for (i = 0; i < nBasisStates; i++) {
		dPsi[i] *= oneOverSqrtInitNorm;
	}

	return;
}

/*****************************************************************************/
//

double retNormOfDoubleWavefunction(double *dPsi, long nBasisStates, double dV) {
	long i;
	double norm = 0.0;

	for (i = 0; i < nBasisStates; i++) {
		norm += sqr(dPsi[i]);
	}
	norm *= dV;

	return norm;
}

/*****************************************************************************/
// 

void normalizeDoubleArray(double *dArray, long nElements) {
	long i;
	double initNorm, oneOverInitNorm;

	initNorm = retNormOfDoubleArray(dArray, nElements);
	oneOverInitNorm = 1.0/initNorm;
	for (i = 0; i < nElements; i++) {
		dArray[i] *= oneOverInitNorm;
	}

	return;
}

/*****************************************************************************/
//

double retNormOfDoubleArray(double *dArray, long nElements) {
	long i;
	double norm = 0.0;

	for (i = 0; i < nElements; i++) {
		norm += dArray[i];
	}

	return norm;
}

/*****************************************************************************/