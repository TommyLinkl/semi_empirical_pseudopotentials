/************************************************************************/
//
// This file contains the code that calculates the trajectories for all 
// the systems using the Gillespie algorithm, solving the classical 
// master equation
//
/************************************************************************/

#include "monteCarlo.h"
 
/************************************************************************/

void calcTrajectories(syst *sys, event *events, param params) {
	int i, j, tStep, nBins = 1 + (int) (params.totalTime / params.printTimeStepSize);
	double t, *avePops, *aveEvents, totalRate;
	int eventID;

	// Dynamically allocate memory
	avePops = (double *) calloc(nBins*params.nPopulations, sizeof(double));
	aveEvents = (double *) calloc(nBins*params.nEvents, sizeof(double));

	for (i = 0; i < params.nSystems; i++) {
		tStep = 0;
		t = 0.0;
		while (t < params.totalTime) {
			tStep++;
			totalRate = calcTotalRate(&sys[i], events, params);
			// advance the time using an exponential time dependence
			t += (1.0/(totalRate+EPS) * log(1.0/randZeroToOne(&(params.seed))));
			calcEnsembleAverage(avePops, &sys[i], t, params);
			if (t > params.totalTime) {
				break;
			}
			// select which event occurs, update the population, and keep track of statistics
			eventID = selectEvent(&sys[i], events, params);
			if (eventID) {
				for (j = 0; j < params.nPopulations; j++) {
					sys[i].pop[j].nParticles += events[eventID-1].popChanges[j];
				}
				//sys[i].pop[0].nParticles = (sys[i].pop[1].nParticles + sys[i].pop[2].nParticles);
				//sys[i].pop[0].nParticles = MIN(sys[i].pop[1].nParticles, sys[i].pop[2].nParticles); // since exciton pop depends on carriers
				//sys[i].pop[3].nParticles = MIN(sys[i].pop[4].nParticles, sys[i].pop[5].nParticles); // since exciton pop depends on carriers
				events[eventID-1].nOccurences++;
				// TODO: fill correct aveEvents bin 				
			}
		}
	}

	writeAveragedTrajectory(avePops, sys[0], params);

	// Free dynamically allocated memory
	free(aveEvents); free(avePops);

	return;
}

/************************************************************************/
// returns the event id of the process that will be occur 

int selectEvent(syst *s, event *events, param params) {
	int i, j; 
	double randNumber, pTest, totalRate = 0.0;

	totalRate = calcTotalRate(s, events, params);
	if (totalRate < EPS) return 0; // no events are possible

	randNumber= randZeroToOne(&(params.seed));
	pTest = 0.0;
	for (i = 0; i < params.nEvents; i++) {
		pTest += (events[i].baseRate*calcRateScaling(s, &events[i], params)) / totalRate;
		if (randNumber <= pTest) {
			return (i+1); 
		}
	}

	return 0;	// should never get to here
}

/************************************************************************/
// calculates the rate that any event occurs

double calcTotalRate(syst *s, event *events, param params) {
	int i, j, k, count = 0;
	double totalRate = 0.0;

	for (i = 0; i < params.nEvents; i++) {
			totalRate += events[i].baseRate*calcRateScaling(s, &events[i], params);
	}

	return totalRate;
}

/************************************************************************/

double calcRateScaling(syst *s, event *e, param params) {
	int i, scaling = 1;
	int count = 0;

	// check to see if event is allowed to occur (min particles present and no max pop restrictions)
	for (i = 0; i < params.nPopulations; i++) {
		if (s->pop[i].nParticles >= e->minParticles[i]) {
			if (! s->pop[i].maxPop) count++; // maxPop = 0 means there is not a max restrictions on pop[j]
			else if ((s->pop[i].nParticles+e->popChanges[i]) <= s->pop[i].maxPop) count++;
		}
	}
	// if event is allowed to occur then calculates its scaling
	if (count == params.nPopulations) {
		for (i = 0; i < params.nPopulations; i++) {
			if (e->rateScaling[i]) scaling *= (calcNChooseK(s->pop[i].nParticles, e->rateScaling[i]));
		}	
	}
	else scaling = 0; // event cannot occur

	return (double)(scaling);
}

/************************************************************************/

int calcNChooseK(int n, int k) {
	if (n >= k) {
		if (k == 1) {
			return n;
		}
		else if (k == 2) {
			return ((n*(n-1))/2);
		}
		else if (k == 3) {
			return ((n*(n-1)/6)*(n-2));
		}
		else {
			return (calcFactorial(n) / (calcFactorial(k)*calcFactorial(n-k)));
		}
	}
	else return 0;
}

/************************************************************************/

int calcFactorial(int n) {
	int i, a = 1;

	if (! n) return 1;
	else {
		for (i = n; i > 0; i--) a=a*i;
		return a;
	}
}

/************************************************************************/

void calcEnsembleAverage(double *ave, syst *s, double t, param params) {
	int i;
	static double lastTime = 0.0;
	static long startIndex = 0;
	
	while (lastTime < t && lastTime < params.totalTime && startIndex < 1 + (int)(params.totalTime / params.printTimeStepSize)) {
		for (i = 0; i < params.nPopulations; i++) {
			ave[startIndex*params.nPopulations+i] += s->pop[i].nParticles;	
		}
		startIndex++;
		lastTime += params.printTimeStepSize;
	}
	lastTime = t;

	// signifies that a new systems trajectory will be sent here when 
	// this function is called next 
	if (t > params.totalTime) {
		startIndex = 0;
		lastTime = 0.0;
		return;
	}

	return;
}

/************************************************************************/
