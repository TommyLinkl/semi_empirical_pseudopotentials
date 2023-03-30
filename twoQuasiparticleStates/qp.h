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
#define AUTOMEV   27211.4
#define AUTOLAMBDANM  45.5640   // nm = AUTOLAMBDANM / energyInAU
#define AUTONM    0.0529177249
#define AUTOANG   0.529177249
#define AUTONS    2.41884e-8
#define AUTOPS    2.41884e-5    
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
  double fermiEnergy, sigmaCutoff;
  double maxHoleEnergyToPrint, minElecEnergyToPrint;
  double minIntTwoQPEnergyToPrint;
  double maxElecDeltaE, maxHoleDeltaE;
  double electronicTemperature, energyLevelSigma;
  double maxEnergyConservation; 
  vector epsilon;
  vector fixedQPPoints[10];
} dParams;

// Long parameters structure
typedef struct lParams_ {
  char calcType[20];
  long calcSinglets, calcTriplets;
  long calcNonintDensitiesOnly, calcIntDensitiesOnly, nFixedQPPoints;
  long readInCoulombMatrixElements, readInGrid;
  long readPsiValueByValue, debug;
  long nThreads, nGridPoints;
  long iHomo, iLumo, nTotalNonintQPs;
  long nHoles, nElecs, nHolesPlusElecs;
  long maxElecStates, maxHoleStates;
  long nHoleEigenstates, nElecEigenstates;
  long maxElecStatesToPrint, maxHoleStatesToPrint; 
  long maxIntTwoQPStatesToPrint, maxIntTwoQPStatesExcSize;
  long nNonintExcitons, nIntExcitons;
  long nQPStates, nNonintTwoQPStates, nIntTwoQPStates;
  long iElecIndex, iHoleIndex, fElecIndex, fHoleIndex;
  long initialStateAverage;
} lParams;

// Coulomb matrix index structure
typedef struct coulombMatrixElement_ {
  long index;
  long type;
  long screenedFlag;
  long alreadyCalculatedFlag;
  long iR, iS, iU, iT;
  long rsutIndex;
  double *psiR, *psiS, *psiU, *psiT;
  double value;
} coulombMatrixElement;

// Overall Coulomb matrix element structure
typedef struct cmeController_ {
	long index;
	long nCME;
	long nR, nS, nU, nT;
	long nGridPointsForPsis;
	long realPsiFlag;
	coulombMatrixElement *cme;
} cmeController;

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
  double energy;    // sum of qp1 and qp1 energies
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

// A noninteracting exciton structure
// TODO: remove or make it contain nonintTwoQPState
typedef struct nonintExc_ {
  long index;       // unique index of the exciton
  double energy;    // exciton energy
  nonintQP *e, *h;  // pointer to electron (e) and hole (h)
} nonintExc;

// An interacting exciton structure
// TODO: remove or make it contain intTwoQPState
typedef struct intExc_ {
  long index;      // unique index of the exciton
  double energy;   // exciton energy
  double *Cai;     // coefficients (i*nElecs+a to access Cai, there are nHoles i states and nElecs a states)
//  nonintQP *e, *h; // pointers to the first electron (e) and hole (h) states
//  nonintExc *iNonintExc; //
} intExc;

/****************************************************************************/

// Functions in twoQPHamiltonian.c
long calcTwoQPHamiltonianWrapper(double *hMatrix, double *h0Matrix, double *Wrsut, double *Vrsut, 
                                 nonintTwoQPState *nonintTwoQP, lParams lPar);
long calcSpinPolarizedTwoQPHamiltonian(double *hMatrix, double *h0Matrix, double *Wrsut, double *Vrsut, 
                                       nonintTwoQPState *nonintTwoQP, lParams lPar);
long calcTwoQPHamiltonian(double *hMatrix, double *h0Matrix, double *Wrsut, double *Vrsut, 
                          nonintTwoQPState *nonintTwoQP, lParams lPar);
long calcTwoQPH0Matrix(double *h0Matrix, nonintTwoQPState *nonintTwoQP, long nNonintTwoQPStates);

// Functions in coulomb.c
long calcAllCoulombMatrixElements(double *Wrsut, double *Vrsut, nonintTwoQPState *nonintTwoQP, double *psiHoles, double *psiElecs, 
          zomplex *qSpaceCoulombPot, zomplex *qSpaceScreenedCoulombPot, grid qSpaceGrid, grid rSpaceGrid, 
          lParams lPar, fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwPsi);
long calcCoulombMatrixElementsOpenMP(coulombMatrixElement *coulombME, long nCME, long *rsutIndexToCMEIndexList, 
                                    zomplex *qSpaceCoulombPot, long nGridPoints, double rSpaceVolumeElement, long nThreads, 
                                    fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwPsi, 
                                    char *iRhoPair, long printCMEValuesFlag, char *cmeValuesFileName);
long deepCopySingleCoulombMatrixElement(coulombMatrixElement *destCME, long iDestCME, 
                    coulombMatrixElement *origCME, long iOrigCME);
long resetAlreadyCalculatedFlagCME(coulombMatrixElement *coulombMatElement, long nCoulombMatrixElements);
long initCoulombMatrixElementStruct(coulombMatrixElement *coulombMatElement, long nCoulombMatrixElements);
long calcQSpaceCoulombPotential(zomplex *qSpaceCoulombPot, grid qSpaceGrid, zomplex *rSpaceCoulombPot, grid rSpaceGrid, 
                vector epsilon, long screenedFlag, fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwPsi);
long calcRSpaceCoulombPotential(zomplex *rSpaceCoulombPot, grid rSpaceGrid, vector epsilon, long screenedFlag);
double screenedCoulomb(double r, double gamma);
long initCoulombMatrixElementController(cmeController *VrsutController, coulombMatrixElement *cme, long nCME);
long calcAllCoulombMatrixElementsRSUTIndices(coulombMatrixElement *cme, long nCME, 
                                            long nR, long nS, long nU, long nT);
long calcRSUTIndexToCMEIndexList(long *rsutIndexToCMEIndexList, coulombMatrixElement *cme, long nCME);

// Functions in augerDecay.c
void calcAugerDecayNoninteracting(nonintTwoQPState *nonintTwoQP, double *Wrust, double *Vrsut, 
                                  dParams dPar, lParams lPar);

// Functions in spin.c
double retSpinZProjNonintTwoQPState(nonintTwoQPState nonintTwoQP);
double retTotalSpinNonintTwoQPState(nonintTwoQPState nonintTwoQP);
double retSpinZProjIntTwoQPState(intTwoQPState intTwoQP, long nNonintTwoQPStates);
double retTotalSpinIntTwoQPState(intTwoQPState intTwoQP, long nNonintTwoQPStates);

// Functions in diag.c
long diagonalizeTwoQPHamiltonain(intTwoQPState *intTwoQP, double *hMatrix, 
                nonintTwoQPState *nonintTwoQP, lParams lPar);
void diag(const long long n, int nthreads, double *mat, double *eval);

// Functions in hartree.c
void hartree(zomplex *rho, zomplex *potq, double *poth, long nGridPoints, 
              fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi);

// Functions in absorption.c
long calcAbsorptionProperties(intTwoQPState *intTwoQP, nonintTwoQPState *nonintTwoQP, grid rSpaceGrid, grid qSpaceGrid,
                dParams dPar, lParams lPar, fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwPsi);
long calcIntMagneticMoments(vector *intMagDipMoment, vector *magDipoleMoment, intTwoQPState *intTwoQP, 
                long nNonintStates, long nIntStates);
long calcNonintAbsorptionProperties(nonintTwoQPState *nonintTwoQP, grid rSpaceGrid, grid qSpaceGrid, dParams dPar, lParams lPar,
                  fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwPsi);
long calcNonintMagneticMoments(vector *magDipoleMoment, nonintTwoQPState *nonintTwoQP, grid rSpaceGrid, grid qSpaceGrid, 
                long nStates, fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwPsi);
long calcRadiativeLifetimesElecDipoleApprox(vector *intElecDipMoment, vector *elecDipoleMoment, intTwoQPState *intTwoQP,
                nonintTwoQPState *nonintTwoQP, dParams dPar, lParams lPar);
long calcIntElecDipoleMoments(vector *intElecDipMoment, vector *elecDipoleMoment, intTwoQPState *intTwoQP, 
                long nNonintStates, long nIntStates);
long calcNonintElecDipoleMoments(vector *elecDipoleMoment, nonintTwoQPState *nonintTwoQP, grid rSpaceGrid, long nStates);

// Functions in densities.c
void calcQPDensitiesFixedOtherQP(intTwoQPState *intTwoQP, double *psiHoles, double *psiElecs, 
                                  grid rSpaceGrid, dParams dPar, lParams lPar);
void calcIntTwoQPProjDensities(intTwoQPState *intTwoQP, grid rSpaceGrid, dParams dPar, lParams lPar);
void calcNonintQPDensities(nonintQP *holeQP, nonintQP *elecQP, grid rSpaceGrid, dParams dPar, lParams lPar);
void calcXYZProjectedProbDensities(double *rho, grid rSpaceGrid, char *fileName);
void calcProbDensity(double *rho, double *psi, long nGridPoints, double dV);
void calcPsiSquared(double *rho, double *psi, long nGridPoints);
void zeroDoubleArray(double *dArray, long nElements);

// Functions in read.c
long readInputParFile(grid *rSpaceGrid, dParams *dPar, lParams *lPar);
long readEvalPsiParFiles(double *psiHoles, double *psiElecs, nonintQP *holeQP, nonintQP *elecQP, 
                          double sigmaCutoff, lParams lPar);
long readNonIntEigenstates(double *psi, double *eval, double *sige, char *psiFileName, char *evalFileName, 
                          double sigmaCutoff, long nGridPoints); 
long readIntTwoQuasiparticleParFiles(intTwoQPState *intTwoQP, double *hMatrix, nonintTwoQPState *nonintTwoQP, 
                                  dParams dPar, lParams lPar);
void readBetheSalpeterResults(double *Cai, double *Ebs, long nNonintExcitons, long nIntExcitons);
long readCoulombMatrixElements(coulombMatrixElement *cme, long *rsutIndexToCMEIndexList, 
                                  cmeController VrsutController, char *fileName);
long retMatchingCoulombMatrixElementIndex(coulombMatrixElement *coulombME, long nMaxCME, 
                                          long iR, long iS, long iU, long iT);
long readConfFile(double *rx, double *ry, double *rz, atm_st *atm, long nAtoms, FILE *pf);
long readGridParFile(grid *rSpaceGrid, char *gridFileName);
long assign_atom_number(char atyp[2]);
void assign_atom_type(char *atype,long j);

// Functions in excitonSize.c
long calcExcitonSizeStats(intTwoQPState *intTwoQP, nonintTwoQPState *nonintTwoQP, double *psiHoles, 
              double *psiElecs, grid rSpaceGrid, dParams dPar, lParams lPar);
long calcExcSizeInteractingExcitons(intTwoQPState *intTwoQP, vector *Oij, vector *Oab, vector *Uij, 
                    vector *Uab, lParams lPar);
long calcExcSizeNoninteractingExcitons(nonintTwoQPState *nonintTwoQP, vector *Oij, vector *Oab, vector *Uij, 
                    vector *Uab, lParams lPar);
long calcUrsMatrix(vector *Urs, double *psiR, long nRStates, double *psiS, long nSStates, grid rSpaceGrid);
long calcOrsMatrix(vector *Ors, double *psiR, long nRStates, double *psiS, long nSStates, grid rSpaceGrid);

// Functions in states.c
long fillNonintTwoQPStructure(nonintTwoQPState *nonintTwoQP, nonintQP *holeQP, nonintQP *elecQP, lParams lPar);
long fillIntExcitonStructure(intExc *intExciton, double *Cai, double *Ebs, lParams lPar);
long fillNonintExcitonStructure(nonintExc *nonintExciton, nonintQP *holeQP, nonintQP *elecQP, lParams lPar);
long storeNonintEigenstates(double *psibe, double *evalbe, double *sigebe, dParams par, lParams ist);
long removeNonEigenstates(double *psi, double *eval, double *sige, double sigmaCutoff, long nGridPoints, long nOrigStates);
long remNonIntStatesOutsideERange(double *psi, double *eval, double *sige, double minE, double maxE, long nGridPoints, long nOrigStates);
long remNonIntStatesInsideERange(double *psi, double *eval, double *sige, double minE, double maxE, long nGridPoints, long nOrigStates);

// Functions in grid.c
long initRSpaceGrid(grid *rSpaceGrid, gridPoint *rSpaceGP, lParams lPar);
long determineRSpaceGridSize(grid *rSpaceGrid, long nAtoms, FILE *confFilePointer);
long initQSpaceGrid(grid *qSpaceGrid, gridPoint *qSpaceGP, grid rSpaceGrid);
long retIndexOfClosestGridPoint(grid rSpaceGrid, vector targetPoint);

// Functions in errorHandling.c
void memoryError(char *str);

// Functions in norm.c
void normalizeAllDoubleWavefunctions(double *dPsi, long nBasisStates, double dV, long nPsis);
void normalizeDoubleWavefunction(double *dPsi, long nBasisStates, double dV);
double retNormOfDoubleWavefunction(double *dPsi, long nBasisStates, double dV);
void normalizeDoubleArray(double *dArray, long nElements);
double retNormOfDoubleArray(double *dArray, long nElements);

// Functions in size.c
double retRangeOfDoubleArray(double *dArray, long nElements);

// Functions in write.c
long writeInputParameters(dParams dPar, lParams lPar, FILE *pf);
long writeNonintExcitonState(nonintExc niExc, FILE *pf);
long writeNonintTwoQPVectorProperty(vector *vObservable, nonintTwoQPState *nonintTwoQP, 
                                    long nNonintTwoQPStates, char *fileName);
long writeSpinPolarizedNonintTwoQPVectorProperty(vector *vObservable, nonintTwoQPState *nonintTwoQP, 
                                                 long nNonintTwoQPStates, char *fileName);
long writeIntTwoQPVectorProperty(vector *vObservable, intTwoQPState *intTwoQP, 
                                  long nIntTwoQPStates, char *fileName);
long writeSpinPolarizedIntTwoQPVectorProperty(vector *vObservable, intTwoQPState *intTwoQP, 
                                              long nIntTwoQPStates, char *fileName);
long writeIntTwoQPStructuresToFiles(intTwoQPState *intTwoQP, long nNonintTwoQPStates, long nIntTwoQPStates);
long writeQuasiparticleState(nonintQP qp, FILE *pf);
void writeCubeFile(double *rho, grid rSpaceGrid, char *fileName);
long writeLongArray(long *longArray, long arrayLength, char *fileNameBase);
long writeDoubleArray(double *doubleArray, long arrayLength, char *fileNameBase);
long writeComplexArray(zomplex *zomplexArray, long arrayLength, char *fileNameBase);
long writeCurrentTime(FILE *pf);
long writeSeparation(FILE *pf);

/****************************************************************************/

#endif

/****************************************************************************/
