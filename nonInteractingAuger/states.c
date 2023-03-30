#include "fd.h"

/****************************************************************************/
// stores specific non-interacting eigenstates that are a desired energy range

// void store_sp_eigenstates(double *psi, double *eval, double *de, ... ) {
// 	long i;




// 	return;
// }


/****************************************************************************/
// returns the number of initial biexcitonic states that are in the desired
// energy range and stores their respective energies

long calc_sp_biexcitonic_states(double *eval, double *de, lng_st ist, par_st par) {
	FILE *pf;
	long numBiexcitons = 0; // keeps track of the number of allowed initial states
	long b, j; // indices of the 1st initial exciton
	long c, k; // indices of the 2nd initial exciton
	double biexcEnergy; // test biexciton energy to see if it is within allowed energy range

	pf = fopen("biexcStates.dat","w");
	for (b = 0; b < ist.numBandEdgeElectrons; b++) {
		for (j = 0; j < ist.numBandEdgeHoles; j++) {
			for (c = 0; c < ist.numBandEdgeElectrons; c++) {
				for (k = 0; k < ist.numBandEdgeHoles; k++) {
					biexcEnergy = eval[ist.lumoIndex+b]+eval[ist.lumoIndex+c]-eval[j]-eval[k];
					fprintf(pf, "%ld %ld %ld %ld %.10f\n", ist.lumoIndex+b, j,
						ist.lumoIndex+c, k, biexcEnergy);
					if (biexcEnergy < (par.maxInitE+EPS)) numBiexcitons++;
				}
			}
		}
	}
	fclose(pf);

	return numBiexcitons;
}

/****************************************************************************/
