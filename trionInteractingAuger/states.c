#include "fd.h"

/****************************************************************************/
// returns the number of initial biexcitonic states that are in the desired
// energy range and stores their respective energies

long calc_bse_biexcitonic_states(double *Ebs, double *Cbs, double *Hbs, lng_st ist, par_st par) {
	FILE *pf;
	long numBiexcitons = 0; // keeps track of the number of allowed initial states
	long i, iExc1, iExc2; // indices of the two excitons of the initial state
	double biexcEnergy; // test biexciton energy to see if it is within allowed energy range

	// calculate and store the individual exciton energies - should be the same as in exciton.par
	for (i = 0; i < ist.numExcitons; i++) {
		Ebs[i] = retExcitonEnergy(&Cbs[i*ist.msbs2], Hbs, ist.msbs2);
		printf("Excitonic state %ld has energy = %.10f\n", i, Ebs[i]);
	}

	// check all possible combinations of two excitonic states are in allowed initial state
	// energy range (less than maximum allowed initial state energy) 
	pf = fopen("biexcStates.dat","w");
	for (iExc1 = 0; iExc1 < ist.numExcitons; iExc1++) {
		for (iExc2 = 0; iExc2 < ist.numExcitons; iExc2++) {
			biexcEnergy = Ebs[iExc1] + Ebs[iExc2];
			if (biexcEnergy < par.maxInitE+EPS) {
				fprintf(pf, "%ld %ld %ld %.10f %.10f %.10f\n", numBiexcitons, iExc1, iExc2, 
					Ebs[iExc1], Ebs[iExc2], biexcEnergy);
				numBiexcitons++;
			}
		}
	}
	fclose(pf);

	return numBiexcitons;
}

/****************************************************************************/
