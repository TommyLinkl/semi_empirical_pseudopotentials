/*****************************************************************************/
//
// This file calculates auger decay rates
//
/*****************************************************************************/

#include "qp.h"

/****************************************************************************/
// 

void calcAugerDecayNoninteracting(nonintTwoQPState *nonintTwoQP, double *Wrust, double *Vrsut, 
									dParams dPar, lParams lPar) {
	FILE *pFile;
	long index, iIS, iFS, nFinalStates = 0;
	long nInitialStates = 0;
	long initialStateList[lPar.nNonintTwoQPStates];
	double initialStateEnergy, boltzmannWeight; 
	double zeroOfEnergy = nonintTwoQP[lPar.iHoleIndex*lPar.nElecs + lPar.iElecIndex].energy;
	double pF = 0.0;
	double augerDecayRate = 0.0;

	// Write beginning of function
	writeSeparation(stdout); writeCurrentTime(stdout);
	fprintf(stdout, "Beginning calculation of the noninteracting auger decay rate\n"); fflush(stdout);

	// Set constants
	double temperature = dPar.electronicTemperature;
	double beta = AUTOEV/(KB*temperature);

	// DEBUG printing
	fprintf(stdout, "iHoleIndex = %ld\n", lPar.iHoleIndex);
	fprintf(stdout, "iElecIndex = %ld\n", lPar.iElecIndex);
	fprintf(stdout, "zeroOfEnergy = %.16f\n", zeroOfEnergy);
	fprintf(stdout, "temperature  = %.16f\n", temperature);
	fprintf(stdout, "beta = %.16f\n", beta); fflush(stdout);

	// 
	for (iIS = 0; iIS < lPar.nNonintTwoQPStates; iIS++) {
		if (fabs(nonintTwoQP[iIS].energy - zeroOfEnergy) > 3.0/beta) {
			initialStateList[iIS] = 0;
		}
		else if (nonintTwoQP[iIS].qp1->index == lPar.iHoleIndex &&
			nonintTwoQP[iIS].qp2->index == lPar.iElecIndex) {
			pF += exp(-beta*(nonintTwoQP[iIS].energy - zeroOfEnergy));
			initialStateList[iIS] = 1;
			nInitialStates++;
		}
		else if (nonintTwoQP[iIS].qp1->index >= lPar.iHoleIndex && 
				 nonintTwoQP[iIS].qp2->index >= lPar.iElecIndex &&
				 lPar.initialStateAverage) {
			pF += exp(-beta*(nonintTwoQP[iIS].energy - zeroOfEnergy));
			initialStateList[iIS] = 1;
			nInitialStates++;
		}
		else {
			initialStateList[iIS] = 0;
		}
	}

	// 
	for (iIS = 0; iIS < lPar.nNonintTwoQPStates; iIS++) {
		if (initialStateList[iIS]) {
			initialStateEnergy = nonintTwoQP[iIS].energy;
			boltzmannWeight = exp(-beta*(initialStateEnergy - zeroOfEnergy));
			for (iFS = 0; iFS < lPar.nNonintTwoQPStates; iFS++) {
				if (iIS != iFS && fabs(nonintTwoQP[iFS].energy - initialStateEnergy) < dPar.maxEnergyConservation) {
					if ( (lPar.fHoleIndex < 0 && nonintTwoQP[iFS].qp2->index == lPar.fElecIndex) ||
						 (lPar.fElecIndex < 0 && nonintTwoQP[iFS].qp1->index == lPar.fHoleIndex) ) {
						index = iIS*lPar.nNonintTwoQPStates + iFS;
						// TODO: check screening and final states - add filter capabilities
						augerDecayRate += boltzmannWeight*sqr(Wrust[index]-Vrsut[index]);
						if (iIS == (lPar.iHoleIndex*lPar.nElecs + lPar.iElecIndex)) {
							nFinalStates++;
						}
						//nFinalStates++;
					}
				}
			}
		}
	}

	// Scale auger decay rate
	augerDecayRate *= TWOPI/(2.0*dPar.maxEnergyConservation)/pF;

	// Print out the partition function
	fprintf(stdout, "Partition function = %.16f\n", pF);
	fprintf(stdout, "Number of initial states = %ld\n", nInitialStates); 
	fprintf(stdout, "Number of final states = %ld\n", nFinalStates);
	fprintf(stdout, "Auger decay rate = %.16f %.16f\n", augerDecayRate, augerDecayRate*AUTOPS); 
	fflush(stdout);

	return;
}

/*****************************************************************************/