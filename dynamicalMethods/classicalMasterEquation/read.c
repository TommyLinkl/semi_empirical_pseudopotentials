/************************************************************************/
//
// This file reads the input parameters from params.par and stores them 
// in a param structure and from system.par and stores them population,
// events, and system structures
//
/************************************************************************/

#include "monteCarlo.h"

/************************************************************************/

void readParams(param *params) {
	FILE *pf;
	int i = 0;
	char field[1000], tmp[1000];

	pf = fopen("params.par", "r");
	while (1) {
		fscanf(pf, "%s", field);
		if (! strcmp(field, "totalTime")) fscanf(pf, "%s %lg", tmp, &(params->totalTime));
		else if (! strcmp(field, "printTimeStepSize")) fscanf(pf, "%s %lg", tmp, &params->printTimeStepSize);
		else if (! strcmp(field, "aveInitialPop")) fscanf(pf, "%s %lg", tmp, &params->aveInitialPop);
		else if (! strcmp(field, "nEvents")) fscanf(pf, "%s %d", tmp, &params->nEvents);
		else if (! strcmp(field, "nPopulations")) fscanf(pf, "%s %d", tmp, &params->nPopulations);		
		else if (! strcmp(field, "nSystems")) fscanf(pf, "%s %d", tmp, &params->nSystems);
		else if (! strcmp(field, "seed")) fscanf(pf, "%s %ld", tmp, &params->seed);
		else if (! strcmp(field, "initPopDist")) fscanf(pf, "%s %s", tmp, &params->initPopDist);
		else {
   			printf("Invalid input field and/ or format of input\n");
   			exit(EXIT_FAILURE);
		}	
		i++;
		if (i > 8) break;
	}
	fclose(pf);

	if (params->aveInitialPop > 1.0) params->maxPop = 10*(int)params->aveInitialPop;
	else params->maxPop = 10; 

	return;
}

/************************************************************************/

void readSystem(syst *sys, event *events, param params) {
	FILE *pf;
	int i, j, k;
	char field[1000], tmp[1000];

	pf = fopen("system.par", "r");
	for (i = 0; i < params.nPopulations; i++) {
		fscanf(pf, "%s %d", tmp, &sys[0].pop[i].id);
		for (j = 0; j < 4; j++) {
			fscanf(pf, "%s", field);
			if (! strcmp(field, "name")) fscanf(pf, "%s %s", tmp, &sys[0].pop[i].name);
			else if (! strcmp(field, "initPopBoolean")) fscanf(pf, "%s %d", tmp, &sys[0].pop[i].initPopBoolean);
			else if (! strcmp(field, "maxPop")) fscanf(pf, "%s %d", tmp, &sys[0].pop[i].maxPop);
			else if (! strcmp(field, "energy")) fscanf(pf, "%s %lg", tmp, &sys[0].pop[i].energy);
			else break;
		}
	}
	for (i = 0; i < params.nEvents; i++) {
		fscanf(pf, "%s %d", tmp, &events[i].id);
		for (j = 0; j < 5; j++) {
			fscanf(pf, "%s", field);
			if (! strcmp(field, "name")) fscanf(pf, "%s %s", tmp, &events[i].name);
			else if (! strcmp(field, "baseRate")) fscanf(pf, "%s %lg", tmp, &events[i].baseRate);
			else if (! strcmp(field, "minParticles")) for (k = 0; k < params.nPopulations; k++) fscanf(pf, "%d", &events[i].minParticles[k]);
			else if (! strcmp(field, "rateScaling")) for (k = 0; k < params.nPopulations; k++) fscanf(pf, "%d", &events[i].rateScaling[k]);
			else if (! strcmp(field, "popChanges")) for (k = 0; k < params.nPopulations; k++) fscanf(pf, "%d", &events[i].popChanges[k]);
		}
	}
	fclose(pf);

	for (i = 1; i < params.nSystems; i++) {
		for (j = 0; j < params.nPopulations; j++) {
			strcpy(sys[i].pop[j].name, sys[0].pop[j].name); 
			sys[i].pop[j].id = sys[0].pop[j].id;
			sys[i].pop[j].initPopBoolean = sys[0].pop[j].initPopBoolean;
			sys[i].pop[j].maxPop = sys[0].pop[j].maxPop;
			sys[i].pop[j].energy = sys[0].pop[j].energy;
		}
	}
	
	return;
}

/************************************************************************/
