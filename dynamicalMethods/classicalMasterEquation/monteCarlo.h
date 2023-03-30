/************************************************************************/
//
// Header file 
//
/************************************************************************/
/* C library packages that are used */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <sys/time.h>
#include <assert.h>
#include <time.h>

/************************************************************************/
/* These are some common unit conversion and multiplication schemes */

#define sqr(x)       ((x) * (x))
#define cube(x)      ((x) * (x) * (x))
#define forth(x)     ((x) * (x) * (x) * (x))
#define ISWAP(a,b)   {double temp = (a); (a) = (b); (b) = -temp;}
#define MIN(x,y)     (((x) < (y)) ? (x) : (y))
#define AUTOEV    27.2114
#define AUTONM	  0.05291772108
#define PIE       3.14159265358979323846
#define TWOPI     6.28318530717958647692
#define EPS		  0.0000000001

/************************************************************************/
/* Structure declarations */

// parameter structure - contains input parameters read from input.par
typedef struct param {
	double totalTime, printTimeStepSize;
	double aveInitialPop;
	int maxPop, nEvents, nSystems, nPopulations;
	long seed;
	char initPopDist[100];
} param;

// population structure - contains the name of the particle type along with number of particles of that type in the system
typedef struct population {
	char name[50];
	int id;				// unique id to identify population/ loop over 
	int nParticles;
	int initPopBoolean; // 0 (init pop is zero) or 1 (initial poisson distribution)
	int maxPop;	// maximum population that a single system can have at any time, 0 indicates no limit
	double energy;
} population;

// system structure - contains the state (i.e., particle populations) and the system's name
typedef struct system {
	char name[50];
	double energy;
	population *pop; 	// dynamically allocated array of pointers to the population structure
} syst;

// event structure - contains the event properties (e.g., how the event changes the population when it occurs)
typedef struct event {
	char name[50];
	int id; 
	long nOccurences;	// keeps track of how many times an event occurs
	double baseRate;
	int *minParticles;	// defines minumum populations of the particle types for process to occur
	int *rateScaling;	// how rate depends on number of each type of particle present in system
	int *popChanges;	// how the population of each particle type changes when event occurs
} event;

/************************************************************************/
/* Function declarations */

/* Functions that read input */
void readParams(param *params);
void readSystem(syst *sys, event *events, param params);

/* Functions that perform the computational work */
void initPopulations(syst *sys, event *events, param params);
void calcTrajectories(syst *sys, event *events, param params);
int selectEvent(syst *sys, event *events, param params);
double calcTotalRate(syst *s, event *events, param params);
double calcRateScaling(syst *s, event *e, param params);
void calcEnsembleAverage(double *ave, syst *s, double time, param params);

/* Functions that write the output files */
void writeAveragedTrajectory(double *ave, syst s, param params);
void writeResults(event *events, param params);
void writeInput(event *events, param params, FILE *pf);
void writeSeparation(FILE *pf);
void writeTime(FILE *pf);

/* General functions */
double randZeroToOne(long *idum);
void randomize(void);
int randPoissonDist(double ave, int max, long seed);
int calcNChooseK(int n, int k);
int calcFactorial(int n);

/************************************************************************/
