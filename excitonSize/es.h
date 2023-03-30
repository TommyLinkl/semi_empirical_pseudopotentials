/*****************************************************************************/
//
//
//
/*****************************************************************************/

#ifndef ES_H
#define ES_H

/*****************************************************************************/
// Library functions 

// Standard library 
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include <time.h>
#include <assert.h>
#include <fftw3.h>
#include <stdlib.h>
#include <omp.h>
#include "unistd.h"
#include "mkl.h"

// Personal headers
#include "vector.h"

/****************************************************************************/
// Macro definitions

// Common multiplication schemes
#define sqr(x)       ((x) * (x))
#define cube(x)      ((x) * (x) * (x))
#define forth(x)     ((x) * (x) * (x) * (x))

// Common unit conversions
#define AUTOEV    27.2114
#define AUTONM    0.0529177249
#define AUTOANG   0.529177249
#define KB        8.6173303e-5  // units of eV/K
#define PIE       3.14159265358979323846
#define TWOPI     6.28318530717958647692
#define EPS       1.0e-10

/****************************************************************************/
// Structure declarations

// Complex number structure
typedef struct zomplex {
  double re, im;
} zomplex;

// Atom structure
typedef struct st9 {
  long natyp;
  char atyp[3];
} atm_st;

// FFTW3 structure
typedef fftw_plan fftw_plan_loc;

// Double parameters structure
typedef struct dParams_ {
  double lumoEnergy, homoEnergy, fundamentalGap, opticalGap;
  double sigmaCutoff;
  double temp, kbT, boltzEnergyRange;
} dParams;

// Long parameters structure
typedef struct lParams_ {
  long debug, intExcitons;
  long nThreads, nGridPoints;
  long iHomo, iLumo, nTotalNonintQPs;
  long nHoles, nElecs, nHolesPlusElecs; 
  long nNonintExcitons, nIntExcitons;
} lParams;

// Grid point structure
typedef struct gridPoint_ {
  long index;       // unique identifier for each grid point
  vector pos;       // position and magnitude from origin
  double localPot;  // the local potential at the grid point
} gridPoint;

// Overall grid structure
typedef struct grid_ {
  long index;             // unique identifier for each grid
  long nGridPointsX, nGridPointsY, nGridPointsZ;
  long nGridPoints;       // total number of grid point structures gridPoint points to
  long iOrigin;           // index of the origin grid point structure
  vector origin;          // position of the origin for this grid
  double volume, dV;      // volume of the grid (assumes cuboid) and stepSize.x*stepSize.y*stepSize.z -> the volume element
  vector minPos, maxPos;  // minimum and maximum grid point positions, respectively
  vector stepSize;        // dx, dy, dz, dr
  gridPoint *gP;          // pointers to the grid points
} grid;

// A noninteracting quasiparticle structure
typedef struct nonintQP_ {
  long index;       // unique index of the quasiparticle
  char type[10];    // electron (e) or hole (h) 
  double energy;    // particle energy
  double sigma;     // sqrt(<E^2>-<E>^2)
  double *psi;      // pointer to beginning of the state
} nonintQP;

// A noninteracting exciton structure
typedef struct nonintExc_ {
  long index;       // unique index of the exciton
  double energy;    // exciton energy
  nonintQP *e, *h;  // pointer to electron (e) and hole (h)
} nonintExc;

// An interacting exciton structure
typedef struct intExc_ {
  long index;      // unique index of the exciton
  double energy;   // exciton energy
  double *Cai;     // coefficients (i*nElecs+a to access Cai, there are nHoles i states and nElecs a states)
//  nonintQP *e, *h; // pointers to the first electron (e) and hole (h) states
//  nonintExc *iNonintExc; //
} intExc;

/****************************************************************************/

// Functions in excitonSize.c
long calcExcitonSizeStats(intExc *intExciton, nonintExc *nonintExciton, double *psiHoles, double *psiElecs, 
                          grid rSpaceGrid, dParams dPar, lParams lPar);
long calcExcSizeNoninteractingExcitons(nonintExc *nonintExciton, vector *Oij, vector *Oab, vector *Uij, 
                                      vector *Uab, lParams lPar);
long calcExcSizeInteractingExcitons(intExc *intExciton, vector *Oij, vector *Oab, vector *Uij, 
                                    vector *Uab, lParams lPar);
long calcUrsMatrix(vector *Urs, double *psiR, long nRStates, double *psiS, long nSStates, grid rSpaceGrid);
long calcOrsMatrix(vector *Ors, double *psiR, long nRStates, double *psiS, long nSStates, grid rSpaceGrid);

// Functions in read.c
long readInputParFile(grid *rSpaceGrid, dParams *dPar, lParams *lPar);
long readEvalPsiParFiles(double *psiHoles, double *psiElecs, nonintQP *holeQP, nonintQP *elecQP, 
                          double sigmaCutoff, lParams lPar);
long readNonIntEigenstates(double *psi, double *eval, double *sige, char *psiFileName, char *evalFileName, 
                          double sigmaCutoff, long nGridPoints); 
void readBetheSalpeterResults(double *Cai, double *Ebs, lParams ist);
long readConfFile(double *rx, double *ry, double *rz, atm_st *atm, long nAtoms, FILE *pf);
long assign_atom_number(char atyp[2]);
void assign_atom_type(char *atype,long j);

// Functions in states.c
long fillIntExcitonStructure(intExc *intExciton, double *Cai, double *Ebs, lParams lPar);
long fillNonintExcitonStructure(nonintExc *nonintExciton, nonintQP *holeQP, nonintQP *elecQP, lParams lPar);
long storeNonintEigenstates(double *psibe, double *evalbe, double *sigebe, dParams par, lParams ist);
long removeNonEigenstates(double *psi, double *eval, double *sige, double sigmaCutoff, long nGridPoints, long nOrigStates);
long remNonIntStatesOutsideERange(double *psi, double *eval, double *sige, double minE, double maxE, long nGridPoints, long nOrigStates);
long remNonIntStatesInsideERange(double *psi, double *eval, double *sige, double minE, double maxE, long nGridPoints, long nOrigStates);

// Functions in grid.c
long initRSpaceGrid(grid *rSpaceGrid, gridPoint *rSpaceGP, lParams lPar);
long determineRSpaceGridSize(grid *rSpaceGrid, long nAtoms, FILE *confFilePointer);

// Functions in errorHandling.c
void memoryError(char *str);

// Functions in size.c
double get_dot_ligand_size_z(double *rz, long n);

// Functions in write.c
long writeInputParameters(dParams dPar, lParams lPar, FILE *pf);
long writeNonintExcitonState(nonintExc niExc, FILE *pf);
long writeQuasiparticleState(nonintQP qp, FILE *pf);
long writeLongArray(long *longArray, long arrayLength, char *fileNameBase);
long writeDoubleArray(double *doubleArray, long arrayLength, char *fileNameBase);
long writeComplexArray(zomplex *zomplexArray, long arrayLength, char *fileNameBase);
long writeCurrentTime(FILE *pf);
long writeSeparation(FILE *pf);

/****************************************************************************/

#endif

/****************************************************************************/
