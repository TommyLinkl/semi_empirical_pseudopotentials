/************************************************************************/
//
// This file does the printing for the program 
//
//
/************************************************************************/

#include "monteCarlo.h"

/************************************************************************/

void writeAveragedTrajectory(double *ave, syst s, param params) {
  FILE *pf;
  int i, j, tStep = 0;
  int np = params.nPopulations;
  double aveEnergy, t = 0.0;

  pf = fopen("trajectories.dat", "w");
  while (t < params.totalTime + EPS) {
	aveEnergy = 0.0;
    fprintf(pf, "%.5f ", t);
    for (j = 0; j < np; j++) {
	  fprintf(pf, "%.8f ", (ave[np*tStep+j]/(double)params.nSystems));
      aveEnergy += ave[np*tStep+j]*s.pop[j].energy;
	}
	aveEnergy /= ((double)params.nSystems);
	fprintf(pf, "%.5f\n", aveEnergy);
	//fprintf(pf, "\n");
    t += params.printTimeStepSize;
    tStep++;
  }
  fclose(pf);

  return;
}

/************************************************************************/

void writeResults(event *events, param params) {
  FILE *pf;
  int i;

  pf = fopen("output.dat", "w");
  fprintf(pf, "RESULTS\n\n");
  writeTime(pf); 

  for (i = 0; i < params.nEvents; i++) {
    fprintf(pf, "The average number of %s events was = %.5f\n", events[i].name, (double)events[i].nOccurences/(double)params.nSystems);
  }

  writeSeparation(pf);
  writeInput(events, params, pf);
  fclose(pf);

  return;
}

/************************************************************************/

void writeInput(event *events, param params, FILE *pf) {
  int i;

  fprintf(pf, "INPUT PARAMETERS\n\n");
  fprintf(pf, "The number of independent system trajectories was = %d\n", params.nSystems);
  fprintf(pf, "The average initial population was = %.3f\n", params.aveInitialPop);
  fprintf(pf, "The maximum population allowed in a system was = %d\n\n", params.maxPop);
  fprintf(pf, "The total time was: %.5f\n", params.totalTime);
  fprintf(pf, "The base rates used were (in inverse units of time):\n");
  for (i = 0; i < params.nEvents; i++) fprintf(pf, "Base rate for %s was = %.6f\n", events[i].name, events[i].baseRate);
  writeSeparation(pf);  

  return;
}

/************************************************************************/
// used to organize output files 

void writeSeparation(FILE *pf) {
  fprintf(pf, "\n***************************************************************\n\n");
  
  return;
}

/************************************************************************/
// Print the formatted string of time to given file

void writeTime(FILE *pf) {
  time_t rawtime;
  struct tm * timeinfo;
  
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );

  fprintf(pf, "The calculation finished on %d %d %d at %d:%d:%d\n\n",timeinfo->tm_mday, timeinfo->tm_mon + 1, timeinfo->tm_year + 1900, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);

  return;
}

/************************************************************************/
