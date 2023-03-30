/****************************************************************************/
/* These are the library functions that are used within the calculation of the band structure */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <assert.h>

/****************************************************************************/
/* These are some common unit conversion and multiplication schemes that help make the program more readable */

#define sqr(x)       ((x) * (x))
#define cube(x)      ((x) * (x) * (x))
#define forth(x)     ((x) * (x) * (x) * (x))
#define ISWAP(a,b)   {double temp = (a); (a) = (b); (b) = -temp;}
#define AUTOEV    27.2114
#define AUTONM	  0.05291772108
#define PIE       3.14159265358979323846
#define TWOPI     6.28318530717958647692
#define HBAR	  1.0
#define NDIM 	  3
#define MASS 	  1.0
#define EPS     1e-10

/****************************************************************************/
/* Structure declarations */

typedef struct complexnumber {
  double re, im;
} complexnumber;

typedef struct vector {
  double x, y, z;
  double mag;
} vector;

typedef struct atom {
  char symbol[1000];
  int type;
  vector pos;
  double ppParams[4];
  double spinOrbit;
  double fF;
} atom;

typedef struct param {
  int nAtoms, nAtomTypes, nBands;
  int nIterations, nKPoints;
  int nQpoints, nBasisVectors, nMaxBV;
  double scale, cellVolume, maxKE, beta, stepSize;
  double expEg, expKp, expCBM, expVBM, expSOB;
  double kineticEnergyScaling;
} param;

/****************************************************************************/
/* Function declarations */

/* Functions that initialize the system */
void readInput(param *params);
void readSystem(param *params, vector *unitCellvectors, atom *atoms);
void readExpBandStructure(double *expBandStructure, param *params);
void readKPoints(vector *kPoint, param params);
int newAtomType(atom *atoms, int i);
double getCellVolume(vector *cellVectors);
void getRecipSpaceBasisVectors(vector *unitCellvectors, vector *gVectors, param *params);
void initAtomicFormFactors(param params, atom *atoms);
double calcPot(double q, double *param);
void calcBasisVectors(vector *tmpBV, vector *gVectors, param *params);
void finalizeBasisStates(vector *basisStates, vector *tmpBV, param *params);
void quickSortVectors(vector *vectors, int l, int r);
int partitionVectors(vector *vectors, int l, int r);
void calcVSpinOrbit(double *vSpinOrbit, vector *basisStates, vector *kPoint, param params);

/* Functions that involve the Hamiltonian matrix */
void calcLocalPotentialEnergyMatrix(complexnumber *localPotential, vector *basisStates, atom *atoms, param params);
void calcKineticEnergyMatrix(complexnumber *hamiltonian, vector *basisStates, vector kVector, param params);
void calcHamiltonianMatrix(complexnumber *hamiltonian, complexnumber *kineticEnergyMat, complexnumber *localPotential, complexnumber *spinOrbitMatrix, param params);
void calcSpinOrbitPotentialMatrix(complexnumber *spinOrbitMatrix, double *vSpinOrbit, vector *basisStates, vector kVector, int kindex, atom *atoms, param params);
void diagnolizeHamiltonian(complexnumber *hamiltonian, double *energies, param params);
void calcBandStructure(double *bandStructure, complexnumber *hamiltonian, complexnumber *kineticEnergyMat, complexnumber *localPotential, complexnumber *spinOrbitMatrix, double *vSpinOrbit, vector *basisStates, vector *kPoint, double *energies, atom *atoms, param params);
void zeroComplexMatrix(complexnumber *mat, int n1, int n2);

/* Functions needed for the Monte Carlo fitting process */
void runMonteCarlo(complexnumber *hamiltonian, complexnumber *kineticEnergyMat, complexnumber *localPotential, complexnumber *spinOrbitMatrix, double *vSpinOrbit, vector *basisStates, double *expBandStructure, vector *kPoint, double *bandStructure, double *energies, atom *atoms, param params);
double calcMSE(double *array1, double *array2, int arrayLength);
void savePseudoParams(double *tmpPseudos, atom *atoms, param params);
void updatePseudoParams(atom *atoms, param params);
void copyPreviousAtomParams(atom *atoms, int i);
void revertPseudoParams(atom *atoms, double *tmpPseudos, param params);
void storeEnergies(double *bandStructure, double *energies, int index, int length);
double randNumber(void);

/* Functions that write the output */
void writeResults(double *bandStructure, vector *kPoint, double bestMSE, atom *atoms, param params);
void writeIteration(double MSE, int iteration, atom *atoms, param params, FILE *pf);
void writeInput(atom *atoms, param params, FILE *pf);
void writePseudoParams(atom *atoms, param params, FILE *pf);
void writeVector(vector vect, FILE *pf);
void writeCurrentTime(FILE *pf);
void writeSeparation(FILE *pf);

/* Functions for vector structures */
vector retAddedVectors(vector vect1, vector vect2);
vector retSubtractedVectors(vector vect1, vector vect2);
vector retScaledVector(vector vect, double scale);
vector retCrossProduct(vector vect1, vector vect2);
vector retZeroVector(void);
double retDotProduct(vector vect1, vector vect2);
int compareVectorMagnitudes(vector vect1, vector vect2);
int countDistinctMagnitudes(vector *vect, int nVectors);

/****************************************************************************/