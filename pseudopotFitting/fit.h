/*****************************************************************************************/
/* These are the library functions that are used within the calculation of the band structure */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <sys/time.h>
#include <assert.h>
#include "mkl.h"

/*****************************************************************************************/
/* These are some common unit conversion and multiplication schemes that help make the program more readable */

#define sqr(x)       ((x) * (x))
#define cube(x)      ((x) * (x) * (x))
#define forth(x)     ((x) * (x) * (x) * (x))
#define ISWAP(a,b)   {double temp = (a); (a) = (b); (b) = -temp;}
#define AUTOEV    27.2114
#define AUTONM	  0.05291772108
#define PIE       3.14159265358979323846
#define TWOPI     6.28318530717958647692
#define HBAR	    1.0
#define NDIM 	    3
#define MASS 	    1.0
#define EPS       1e-10
#define NPParams  4

/*****************************************************************************************/
/* Structure declarations */

typedef struct dcomplex_ {
  double real, imag;
} dcomplex;

typedef struct vector_ {
  double x, y, z;
  double mag;
} vector;

typedef struct atom_ {
  char symbol[10];
  int type;
  vector pos;
  double ppParams[4];
  double fF;
} atom;

typedef struct param_ {
  int nBandStructures, nIterations;
  int nAtoms, nAtomTypes, nElectrons;
  int nKPoints, nKPE; 
  int nkp[20], ne[20], na[20];
  int nqp[20], nbv[20], nmbv[20];
  double s[20], mke[20], vol[20], volPerAtom[20];
  double normalizedVolPerAtom[20];
  int nQpoints, nBasisVectors, nMaxBV;
  double scale, maxKE, beta, stepSize;
  double kineticEnergyScaling;
} param;

/*****************************************************************************************/
/* Function declarations */

/* Functions that read input - read.c */
void readInput(param *params);
void readMultiInput(param *params);
void copyBaseParamsToArrays(param *params);
void readExpBandStructure(double *expBandStructure, int nKPoints, int nElectrons, int index);
void readKPoints(vector *kPoint, int nKPoints, double scale, int index);
void readSystem(param *params, vector *unitCellvectors, atom *atoms, int index);

/* Functions that initialize the system - init.c */
int newAtomType(atom *atoms, int i);
double getCellVolume(vector *cellVectors);
int getRecipSpaceBasisVectors(vector *unitCellvectors, vector *gVectors, double maxKE);
void initAtomicFormFactors(atom *atoms, param params);

/* Functions that deal with the basis function generation - states.c */
int calcBasisVectors(vector *tmpBV, vector *gVectors, double maxKE);
void finalizeBasisStates(vector *basisStates, vector *tmpBV, int nMaxBV, double maxKE);
void quickSortVectors(vector *vectors, int l, int r);
int partitionVectors(vector *vectors, int l, int r);

/* Functions needed for the Monte Carlo fitting process */
void runMonteCarlo(dcomplex *hamiltonian, dcomplex *potential, vector *basisStates, 
  double *expBandStructure, vector *kPoint, double *bandStructure, double *energies, atom *atoms, param params);
void calcBandStructure(double *bandStructure, dcomplex *hamiltonian, dcomplex *potential, 
  vector *basisStates, vector *kPoint, double *energies, atom *atoms, param params);
double calcAllIndividualWeightedMSE(double *bandStructure, double *expBandStructure, param params);
double calcWeightedMSE(double *array1, double *array2, double *weights, int nElectrons, int nKPoints);
double calcAllIndividualMSE(double *bandStructure, double *expBandStructure, param params);
double calcMSE(double *array1, double *array2, int arrayLength);
void savePseudoParams(double *tmpPseudos, atom *atoms, param params);
void updatePseudoParams(atom *atoms, param params);
int copyPreviousAtomParams(atom *atoms, int i, param params);
void revertPseudoParams(atom *atoms, double *tmpPseudos, param params);
void storeEnergies(double *bandStructure, double *energies, int index, int length);

/* Functions related to random number generation */
double randNumber(void);

/* Functions that involve the Hamiltonian matrix */
void calcPotentialEnergyMatrix(dcomplex *potential, vector *basisStates, atom *atoms, int nbv, int nAtoms);
void calcKineticEnergyMatrix(dcomplex *hamiltonian, vector *basisStates, vector *kVector, 
                              int nbv, double kineticEnergyScaling);
void calcHamiltonianMatrix(dcomplex *hamiltonian, dcomplex *potential, vector *basisStates, 
                            vector *kVector, int nbv, double kineticEnergyScaling);
double calcPot(double q, double *param);
void zeroComplexMatrix(dcomplex *mat, int n1, int n2);

/* Functions that diagonalize the Hamiltonian matrix - diag.c */
void diagonalizeHamiltonian(dcomplex *hamiltonian, double *energies, int nBasisVectors);

/* Functions that calculate the real space pseudopotentials - pot.c */
void calcRealSpacePotential(atom *atoms, param params);

/* Functions that write the output */
void writeResults(double *bandStructure, vector *kPoint, double bestMSE, atom *atoms, param params);
void writeIteration(double MSE, int iteration, atom *atoms, param params, FILE *pf);
void writeInput(atom *atoms, param params, FILE *pf);
void writePseudoParams(atom *atoms, param params, FILE *pf);
void writeBestPseudoParams(atom *atoms, param params);
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

/*****************************************************************************************/
