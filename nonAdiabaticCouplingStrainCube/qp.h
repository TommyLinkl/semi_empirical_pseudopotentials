/*****************************************************************************/
//
//
//
/*****************************************************************************/

#ifndef QP_H
#define QP_H

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
#define AUTOLAMBDANM  45.5640   // nm = AUTOLAMBDANM / energyInAU
#define AUTONM    0.0529177249
#define AUTOANG   0.529177249
#define AUTONS    2.41884e-8    
#define KB        8.6173303e-5  // units of eV/K
#define PIE       3.14159265358979323846
#define TWOPI     6.28318530717958647692
#define FOURPI    (4.0*3.14159265358979323846)
#define EPS       1.0e-10
#define EPSR      1.0e-10

/****************************************************************************/
// Structure declarations

// Complex number structure
typedef struct zomplex {
  double re, im;
} zomplex;

// Atom structure - old -> Remove at some point
typedef struct st9 {
  long natyp;
  char atyp[3];
} atm_st;

// FFTW3 structure
typedef fftw_plan fftw_plan_loc;

// Atom structure
typedef struct atom_ {
  long type, index;
  long scType;
  char symbol[100];
  double mass;
  vector pos;
} atom;

// Double parameters structure
typedef struct dParams_ {
  double lumoEnergy, homoEnergy, fundamentalGap, opticalGap;
  double sigmaCutoff;
  double fermiEnergy;
  double maxElecDeltaE, maxHoleDeltaE;
  vector epsilon;
} dParams;

// Long parameters structure
typedef struct lParams_ {
  char calcType[20];
  char stateType[20];
  char crystalStructure[20];
  char outmostMaterial[20]; 
  long calcSinglets, calcTriplets;
  long debug;
  long nAtoms, nAtomTypes;
  long nSCAtoms, nSCAtomTypes;
  long nThreads, nGridPoints;
  long iHomo, iLumo, nTotalNonintQPs;
  long nHoles, nElecs, nHolesPlusElecs;
  long maxElecStates, maxHoleStates; 
  long nNonintExcitons, nIntExcitons;
  long nQPStates, nNonintTwoQPStates, nIntTwoQPStates;
  long nPhonons;
} lParams;

// grid point structure for 1D grid
typedef struct gridPoint1d_ {
  long index;         // unique identifier for each grid point
  double pos;         // position
  double localVal;    // value of the function (i.e. potential) at that grid point
} gridPoint1d;

// grid structure for 1D
typedef struct grid1d_ {
  long index;                 // unique identifier for each 1D grid
  long nGridPoints;           // number of grid point structures that gP points to
  long iOrigin;               // index of the origin grid point structure
  double origin;              // position of the origin for this grid
  double minPos, maxPos;      // minimum and maximum grid point positions
  double stepSize;            // dr
  gridPoint1d *gP;              // pointers to the grid points
} grid1d;

// grid point structure for 3D grid
typedef struct gridPoint3d_ {
  long index;         // unique identifier for each grid point
  vector pos;         // position
  double localVal;    // value of the function (i.e. potential) at that grid point
} gridPoint3d;

// grid structure for 3D
typedef struct grid3d_ {
  long index;                 // unique identifier for each 1D grid
  long nGridPoints;           // total number of grid point structures that gP points to
  long nGridPointsX, nGridPointsY, nGridPointsZ;
  long iOrigin;               // index of the origin grid point structure
  vector origin;              // position of the origin for this grid
  vector minPos, maxPos;      // minimum and maximum grid point positions
  vector stepSize;            // dx, dy, dz
  double volume, dV;          // volume of the grid and volume element (stepSize.x*stepSize.y*stepSize.z)
  gridPoint3d *gP;            // pointers to the grid points
} grid3d;

// A noninteracting quasiparticle structure
typedef struct nonintQP_ {
  long index;       // unique index of the quasiparticle
  char type[10];    // electron (e) or hole (h) 
  double spin;      // both electrons and holes are spin 0.5 particles
  double spinZ;     // spin up = 0.5, spin down = -0.5
  double energy;    // particle energy
  double sigma;     // sqrt(<E^2>-<E>^2)
  double *psi;      // pointer to beginning of the state
} nonintQP;

// A noninteracting two quasiparticle state structure
typedef struct nonintTwoQPState_ {
  long index;       // unique index of the exciton
  char type[20];    // elec-hole or elec-elec or hole-hole
  double spin;      // singlet = 0.0 and triplet = 1.0
  double spinZ;     // -spin, -spin+1, ... , spin-1, spin
  double energy;    // exciton energy
  double hartreeEnergy; // Coulomb (i.e. Hartree) energy for the two QP state
  double exchangeEnergy;
  nonintQP *qp1, *qp2;  // pointer to electron (e) and hole (h)
} nonintTwoQPState;

// An interacting (i.e. correlated) two quasiparticle state structure
typedef struct intTwoQPState_ {
  long index;       // unique index of the exciton
  char type[20];    // elec-hole or elec-elec or hole-hole
  double spin;      // singlet = 0.0 and triplet = 1.0
  double spinZ;     // -spin, -spin+1, ... , spin-1, spin
  double energy;    // energy or the correlated two quasiparticle state
  double correlationEnergy; // fundamental gap + coulombEnergy - energy
  double bindingEnergy;     // exciton binding energy for elec-hole pairs
  long nNonintTwoQPStates;
  double *Crs;      // coefficients 
  nonintTwoQPState *niTwoQP; //
} intTwoQPState;

/****************************************************************************/
// functions in nonadiabatic.c
// void calcDerivativeCouplings(vector *Vab, vector *Vij, vector *dPotdR, atom *atoms,
//         grid1d dAtomicPPGrid, nonintQP *holeQP, nonintQP *elecQP, grid3d rSpaceGrid, 
//         lParams lPar);
void calcElPh(vector *Vab, vector *Vij, vector *dPotdR, atom *atoms, atom *atomNeighbors,
        grid1d atomicPPGrid, grid1d dAtomicPPGrid, nonintQP *holeQP, nonintQP *elecQP, 
        grid3d rSpaceGrid, double *strainScale, vector *strainScaleDeriv, lParams lPar);

/****************************************************************************/
// functions in pseudopotential.c
void calcPseudopotentialDerivative(grid1d *dAtomicPPGrid, gridPoint1d *dAtomicPPGP,
        grid1d atomicPPGrid, long nAtomTypes);
void calcPotential(vector *pot, grid3d rSpaceGrid, vector atomPosition,
        double *potentialValues, double *potentialGridPointValues,
        double potGridStepSize, long nPotValues, vector strainScaleDeriv);
void calcdPotentialdR(vector *dPotdR, grid3d rSpaceGrid, vector atomPosition, 
        double *potentialValues, double *potentialGridPointValues, 
        double potGridStepSize, long nPotValues, double strainScale);
double interpolate(double evaluationPoint, double stepSize, double *potentialValues, 
        double *potentialGrid, long nPotValues);

/****************************************************************************/
// functions in grid.c
long initRSpaceGrid(grid3d *rSpaceGrid, gridPoint3d *rSpaceGP, lParams lPar);
long determineRSpaceGridSize(grid3d *rSpaceGrid, long nAtoms, FILE *confFilePointer);

/****************************************************************************/
// Functions in read.c
long readInputParFile(grid3d *rSpaceGrid, dParams *dPar, lParams *lPar);
long readEvalPsiParFiles(double *psiHoles, double *psiElecs, nonintQP *holeQP, 
        nonintQP *elecQP, double sigmaCutoff, lParams lPar);
long readNonIntEigenstates(double *psi, double *eval, double *sige, char *psiFileName, 
        char *evalFileName, double sigmaCutoff, long nGridPoints); 
void readBetheSalpeterResults(double *Cai, double *Ebs, lParams ist);
long readConfFile(double *rx, double *ry, double *rz, atm_st *atm, long nAtoms, FILE *pf);
long assign_atom_number(char atyp[2]);
void assign_atom_type(char *atype,long j);
double retIdealBondLength(long natyp_1, long natyp_2, int crystalStructureInt); 
long fillAtomStructureFromConfFile(atom *atoms, lParams *lPar);
long isAPassivationSymbol(char *atomSymbol);
long isNewAtomType(atom *atoms, long currIndex);
void readAtomicPseudopotentials(grid1d *atomicPPGrid, gridPoint1d *atomicPPGP, atom *atoms, lParams lPar, double *a4Params, double *a5Params);
void readNearestNeighbors(long nAtoms, atom *atomNeighbors);
void calculateRefTetrahedronVol(long nAtoms, int crystalStructure, int outmostMaterial, atom *atoms, atom *atomNeighbors, double *refTetrahedronVol); 
void calculateTetrahedronVol(long nAtoms, atom *atoms, atom *atomNeighbors, double *tetrahedronVol);
void calculateTetrahedronVolDeriv(long nAtoms, atom *atoms, atom *atomNeighbors, double *tetrahedronVol, vector *tetrahedronVolDerivatives);
void calculateStrainScale(long nAtoms, atom *atoms, double *tetrahedronVolRef, double *tetrahedronVol, double *a4Params, double *a5Params, double *strainScale);
void calculateStrainScaleDeriv(long nAtoms, atom *atoms, atom *atomNeighbors, double *tetrahedronVolRef, double *tetrahedronVol, vector *tetrahedronVolDeriv,
        double *a4Params, double *a5Params, vector *strainScaleDeriv);

/****************************************************************************/
// Functions in errorHandling.c
void memoryError(char *str);

/****************************************************************************/
// Functions in size.c
double get_dot_ligand_size_z(double *rz, long n);

/****************************************************************************/
// Functions in write.c
long writeInputParameters(dParams dPar, lParams lPar, FILE *pf);
long writeQuasiparticleState(nonintQP qp, FILE *pf);
long writeLongArray(long *longArray, long arrayLength, char *fileNameBase);
long writeDoubleArray(double *doubleArray, long arrayLength, char *fileNameBase);
long writeComplexArray(zomplex *zomplexArray, long arrayLength, char *fileNameBase);
long writeCurrentTime(FILE *pf);
long writeSeparation(FILE *pf);

/****************************************************************************/

#endif

/****************************************************************************/
