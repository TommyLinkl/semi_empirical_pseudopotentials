/************************************************************************/
//
// This program solves a coupled set of differential equations 
// (i.e. classical master equations) using the Gillespie algorithm.
//
//
// Input files:
// 		params.par -> a list of 'key = value' pairs on each line
//					  see readParams() in read.c to see required 
//					  'key = value' pairs
// 		system.par -> defines all the populations and the events
//					  see readSystem() in read.c to see required 
//					  'key = value' pairs
//
// 	Output files:
// 		output.dat -> summarizes the calculations (e.g. average number of 
//					  each type of event occured for the given initial 
//					  conditions)
//		trajectories.dat -> solution of the classical master equation
//						    first column contains time, 2nd-2nd to last 
//  						contain the average population for each   
// 							population at the corresponding time and the  
//						    last column contains the energy of the system 
//
// 	Example run commands:
//		$ ./kinMCvDt.x &
//
/************************************************************************/

#include "monteCarlo.h"

/************************************************************************/
/* Initiates the main functions */

int main() {

	/* Defines the variables that are used throughout this function */
	int i; param params; syst *sys; event *events;

	/* Reads input parameters that determine the size of the system */
	readParams(&params);
 	
	/* Allocates the input dependent memory */
	sys = (syst *) calloc(params.nSystems, sizeof(syst));
	for (i = 0; i < params.nSystems; i++) sys[i].pop = (population *) calloc(params.nPopulations, sizeof(population));
	events = (event *) calloc(params.nEvents, sizeof(event));
	for (i = 0; i < params.nEvents; i++) { 
		events[i].minParticles = (int *) calloc(params.nPopulations, sizeof(int));
		events[i].rateScaling = (int *) calloc(params.nPopulations, sizeof(int));
		events[i].popChanges = (int *) calloc(params.nPopulations, sizeof(int));
	}
 
	/* Initializes the system */
	readSystem(sys, events, params);
	initPopulations(sys, events, params);

	/* Performs the computational work */
	calcTrajectories(sys, events, params);
  
	/* Write the output files */
	writeResults(events, params);

	/* Free dynamically allocated memory */ 
	for (i = 0; i < params.nSystems; i++) free(sys[i].pop);
	free(sys);
	for (i = 0; i < params.nEvents; i++) { 
		free(events[i].minParticles);
		free(events[i].rateScaling);
		free(events[i].popChanges);
	}
	free(events);

	return 0;
}

/************************************************************************/
