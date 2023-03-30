/************************************************************************/
//
// Sets up the initial populations for each system
//
/************************************************************************/

#include "monteCarlo.h"

/************************************************************************/
// Initializes the populations of the systems to satisfy 
// the Poisson distribution in the limit of an infinite # of systems

void initPopulations(syst *sys, event *events, param params) {
	int i, j, initRodPop, noPopCount = 0;
	randomize();

	for (i = 0; i < params.nSystems; i++) {
		if (! strcmp(params.initPopDist, "constant")) {
			initRodPop = params.aveInitialPop;
		}
		else if (! strcmp(params.initPopDist, "poisson")) {
			initRodPop = randPoissonDist(params.aveInitialPop, params.maxPop, params.seed);
		}
		else {
			fprintf(stdout, "Invalid initPopDist. constant or poisson allowed\n"); fflush(stdout);
			exit(EXIT_FAILURE);
		}
		if (! initRodPop) noPopCount++;
		for (j = 0; j < params.nPopulations; j++) {
			if (sys[i].pop[j].initPopBoolean) {
				sys[i].pop[j].nParticles = initRodPop;
			}
			else sys[i].pop[j].nParticles = 0;
		}
	}

	//printf("The number of systems with no initial populations = %d\n", noPopCount);	

	return;
}

/************************************************************************/

