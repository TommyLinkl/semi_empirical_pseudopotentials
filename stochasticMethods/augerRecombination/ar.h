/*****************************************************************************/
//
//
//
/*****************************************************************************/

#ifndef AR_H
#define AR_H

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

// Common multiplication and swap schemes
#define sqr(x)       ((x) * (x))
#define cube(x)      ((x) * (x) * (x))
#define forth(x)     ((x) * (x) * (x) * (x))
#define ISWAP(a,b)   {double tmp = (a); (a) = (b); (b) = -tmp;}
#define ISWAP2(a,b)  {double tmp = (a); (a) = 2.0 * (b); (b) = -2.0 * tmp;}

// Common unit conversions
#define AUTOEV    27.2114
#define AUTONS    41341374.575751
#define AUTOPS    41341.374575751
#define KB        8.6173303e-5  // units of eV/K
#define PIE       3.14159265358979323846
#define TWOPI     6.28318530717958647692
#define FOURPI    (4.0*3.14159265358979323846)
#define SVDEPS    1.0e-10
#define EPSR      1.0e-10
#define EPSCHI    1.0e-8
#define DENE      1.0e-10
#define NDOS      500
#define EPS       1.0e-10

/****************************************************************************/
// Structure declarations

// Complex number structure
typedef struct st0 {
  double re, im;
} zomplex;

// FFTW3 structure
typedef fftw_plan fftw_plan_loc;

// Doulbe parameter/ common constant structure
typedef struct st1 {
  double dx, dy, dz, dr, dkx, dky, dkz, dv;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double Ekinmax, gamma, gamma2, Elmin, Elmax;
  double lumoEnergy, homoEnergy, fundamentalGap, opticalGap;
  double Vmin, Vmax, dE, dE_1;
  double sigmaCutoff, maxDeltaE;
  double Eamin, Eamax, Eimin, Eimax;
  double fermiEnergy;
  double temp, kbT, boltzEnergyRange, minInitE, maxInitE;
} par_st;

// Long parameter structure
typedef struct st4 {
  long mfermi, debug;
  long readInPsiEvalFiltFiles, readInHotQPs, readAaiMatrices;
  long nx, ny, nz, nGridPoints;
  double nx_1, ny_1, nz_1, nGridPoints_1;
  long nAtoms, nAtomTypes, seed, nThreads;
  long nNewtonIntSteps, nPseudoPot;
  long nFilterCycles, nStatesPerFilter, nFilteredStates;
  long homoIndex, lumoIndex, totalHomo, totalLumo, msbs, msbs2;
  long nHoles, nElecs, nHolesPlusElecs, nNonIntExcitons, nIntExcitons;
  long nStochOrbitals, nStochHotHoles, nStochHotElecs;
  long nHotHoles, nHotElecs, nHotNonIntStates;
  long imin, imax, amin, amax, atot, itot;
  long nMaxIntExcitons;
  long nConsWindows, nBiexcitons;
  long nTrappedStates;
  long deterministic, stochastic, nonIntAR, intAR;
  char calcType[10]; // dNI, sNI, testNI, dI, s1I, s2I, test1I, test2I
} lng_st;

// Atom structure
typedef struct st9 {
  long natyp;
  char atyp[3];
} atm_st;

// External potential structure
typedef struct st10 {
  double x0, t, e;
} pot_st;

// Grid point structure
typedef struct gridPoint_ {
  long index;       // unique identifier for each grid point
  vector pos;       // position and magnitude from origin
  double localPot;  // the local potential at the grid point
} gridPoint;

// Overall grid structure
typedef struct grid_ {
  long index;       // unique identifier for each grid
  long nGridPoints; // total number of grid point structures gridPoint points to
  long iOrigin;     // index of the origin grid point structure
  vector origin;    // position of the origin for this grid
  double dV;        // stepSize.x*stepSize.y*stepSize.z -> the volume element
  vector stepSize;  // dx, dy, dz, dr
  gridPoint *gP;    // pointers to the grid points
} grid;

// A noninteracting quasiparticle structure
typedef struct nonIntQP_ {
  long index;       // unique index of the quasiparticle
  char type[10];    // electron (e) or hole (h) 
  double energy;    // particle energy
  double sigma;     // sqrt(<E^2>-<E>^2)
  double *psi;      // pointer to beginning of state
} nonIntQP;

// A noninteracting exciton structure
typedef struct nonIntExc_ {
  long index;       // unique index of the exciton
  double energy;    // exciton energy
  nonIntQP *e, *h;  // pointer to electron (e) and hole (h)
} nonIntExc;

// A noninteracting biexciton structure
typedef struct nonIntBiexc_ {
  long index;
  long b, j, c, k;
  double *psiB, *psiJ, *psiC, *psiK;
  double bEnergy, jEnergy, cEnergy, kEnergy;
  double energy;      // biexciton energy, sum of 4 nonIntQP energies
  nonIntExc *bj, *ck; // pointers to exciton 1 (b,j) and 2 (c,k) structures
} nonIntBiexc;

// An interacting exciton structure
typedef struct intExc_ {
  long index;      // unique index of the exciton
  double energy;   // exciton energy
  double *Cai;     // coefficients (i*nElecs+a to access Cai)
  nonIntQP *e, *h; // pointers to the first electron (e) and hole (h) states
} intExc;

// An interacting biexciton structure
typedef struct intBiexc_ {
  long index;      // unique index of the biexciton
  double energy;   // sum of the two exicton energies
  intExc *bj, *ck; // pointers to exciton 1 (b,j) and 2 (c,k) structures
} intBiexc;

/****************************************************************************/

// Functions in init.c 
void init_size(par_st *par, lng_st *ist);
void setCalcTypeParameters(lng_st *ist);
void init(double *vx, double *vy, double *vz, double *ksqr, double *potl, 
          double *rx, double *ry, double *rz, par_st *par, lng_st *ist);
void init_pot(double *vx ,double *vy, double *vz, zomplex *potq, par_st par, lng_st ist, 
              fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi);
void init_psi(zomplex *psi, lng_st ist, par_st par, long *idum);
double expot(double r, pot_st ppar);
double screenedcoulomb(double dr, double gamma);

// Functions in interpolate.c
double interpolate(double r, double dr, double *vr, double *pot, long npot, long n, long j);

// Functions in read.c
long readNonIntEigenstates(double *psi, double *eval, double *sige, char *psiFileName, char *evalFileName, 
                          double sigmaCutoff, long nGridPoints); 
void readBetheSalpeterResults(double *Cbs, double *Ebs, lng_st ist);
long read_conf(double *rx, double *ry, double *rz, atm_st *atm, long n, FILE *);
void read_pot(double *vr, double *pot, long *npot, double *dr, atm_st *atm, long n, long ntype);
long assign_atom_number(char atyp[2]);
void assign_atom_type(char *atype,long j);

// Functions in norm.c
double norm(zomplex *, double,long);
double normalize(zomplex *, double,long);
void normalize_all(double *, double,long,long);
void normalize_all_omp(double *psi, double dr, long ms, long nGridPoints, long nThreads);

// Functions in rand.c
double ran_nrc(long *idum);
double rint(double); 
long random(void);
void srandom(unsigned int);
void Randomize();
double ran();

// Functions in size.c
double get_dot_ligand_size(double *, double *, double *, long);
double get_dot_ligand_size_z(double *rz, long n);

// Functions in gauss.c 
void gauss_test(double *vx, double *vy, double *vz, zomplex *potq, double *poth, par_st par, lng_st ist, 
                fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi);

// Functions in hartree.c
void hartree(zomplex *rho, zomplex *potq, double *poth, long nGrid, 
              fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi);

// Functions in hamiltonian.c
void hamiltonian(zomplex *phi, zomplex *psi, double *potl, double *ksqr, lng_st ist, par_st par, 
                  fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi);
void kinetic(zomplex *psi, double *ksqr, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi, lng_st ist);

// Functions in filter.c
void filtering(double *psitot, double *potl, double *ksqr, zomplex *an, double *zn, double *el, lng_st ist, par_st par, 
               long tid, long jns, long *idum, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi);
void filter(zomplex *psin, zomplex *psim1, double *psi0, double *potl, double *ksqr, par_st par, zomplex *an, double *zn, 
            lng_st ist, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi, long tid, long jns);

// Functions in hnorm.c
void hnorm(zomplex *psim1, zomplex *psin, double *potl, double *ksqr, par_st par, double zm1, lng_st ist, 
            fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi);
void hamiltonian_norm(zomplex *psi, zomplex *phi, double *potl, double *ksqr, lng_st ist, par_st par, 
                      fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi);

// Functions in energy.c
double energy(zomplex *psi, zomplex *phi, double *potl, double *ksqr, lng_st ist, par_st par, 
              fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi);
void energy_all(zomplex *psi, zomplex *phi, double *psi0, double *potl, double *ksqr, double *ene, lng_st ist, par_st par, 
                fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi, long ms);
void get_energy_range(double *vx, double *vy, double *vz, double *ksqr, double *potl, par_st *par, lng_st ist, 
                      fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi);
void calc_sigma_E(double *psitot, double *potl, double *ksqr, double *eval2, lng_st ist, par_st par);

// Functions in states.c
long calcNumIntBiexcStates(double *Ebs, double maxBiexcEnergy, long nIntExcitons, long nTrappedStates);
long calcNumNonIntBiexcStates(double *evalbe, double maxBiexcEnergy, long nHoles, long nElecs, long nTrappedStates);
long storeNonIntEigenstates(double *psibe, double *evalbe, double *sigebe, par_st par, lng_st ist);
long fillNonIntBiexcStruct(nonIntBiexc *biExciton, double *eval, double maxInitE, lng_st ist);
long fillIntBiexcEnergies(double *biexcitonEnergies, double *Ebs, lng_st ist, par_st par);
void storeAllowedFinalStatesNonIntAR(long *finalHoleList, long *nFinalHoles, long *finalElecList, 
        long *nFinalElecs, double *evalai, long nHotHoles, long nHotElecs, 
        nonIntBiexc *initBiexcStates, long nBiexcitons, double *deltaE, long nConsWindows);
void storeAllowedFinalStatesIntAR(long *elecHotHoleList, long *nHotHoleForElec, long *holeHotElecList, long *nHotElecForHole, 
        double *biexcitonEnergies, double *Ebs, double *evalbe, double *evalai, double *deltaE, lng_st ist, par_st par);
long removeNonEigenstates(double *psi, double *eval, double *sige, double sigmaCutoff, long nGridPoints, long nOrigStates);
long remNonIntStatesOutsideERange(double *psi, double *eval, double *sige, double minE, double maxE, long nGridPoints, long nOrigStates);
long remNonIntStatesInsideERange(double *psi, double *eval, double *sige, double minE, double maxE, long nGridPoints, long nOrigStates);
long remUnnecessaryHotNonIntStates(double *psiai, double *evalai, double *sigeai, par_st par, lng_st ist);
long remStatesNotInThetaLists(double *psiai, double *evalai, long *thetaIList, long *thetaAList, lng_st *ist);

// Functions in coeff.c
void coefficient(zomplex *an, double *samp, long nc, long M, double dt, double dE, double Emin, double *El);
void chebyshev_reordered(double *polong, double min, double max, long nc);
double samp_points_ashkenazy(zomplex *polong, double min, double max, long nc);

// Functions in auger.c
void calcDoublyStochasticIntAR(double *Cbs, double *Ebs, double *psibe, double *evalbe, double *psiai, double *evalai, zomplex *potq,
                                lng_st ist, par_st par, fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwpsi);
void calcStochasticCoulombIntAR(double *Cbs, double *Ebs, double *psibe, double *evalbe, double *psiai, double *evalai,
                                zomplex *potq, double *vx, double *vy, double *vz, lng_st ist, par_st par, 
                                fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwpsi);
void calcStochasticFinalStatesIntAR(double *Cbs, double *Ebs, double *vijck, double *vabck, double *psibe, double *evalbe, 
                      double *psiai, double *evalai, zomplex *potq, double *poth, lng_st ist, par_st par, 
                      fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwpsi);
void calcDeterministicIntAR(double *Cbs, double *Ebs, double *vijck, double *vabck, 
                     double *evalbe, double *evalai, lng_st ist, par_st par);
void calcCoherentSumsIntAR(double *w2h, double *w2e, double *Cbs, double *Ebs, 
                          double *vijck, double *vabck, par_st par, lng_st ist);
void calcStochasticNonIntAR(double *vijck, double *vabck, double *psibe, double *evalbe, double *psiai, double *evalai, 
                    zomplex *potq, lng_st ist, par_st par, fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwpsi);
void calcDeterministicNonIntAR(double *vijck, double *vabck, double *evalbe, double *evalai, lng_st ist, par_st par);
void calcBoltzmannWeightedRates(double *energies, double *eRate, double *hRate, double temp, long nStates, FILE *pf, long iD);
double calcPartitionFunction(double *energies, double temp, long nStates);
void fillEnergyConservationWindows(double *deltaE, double maxDeltaE, long nConsWindows);

// Functions in stochasticCoulomb.c
void calcRrsZetaMatrix(zomplex *RrsZeta, double *psiR, long nR, double *psiS, long nS, 
                      zomplex *psiZeta, long nZeta, long nGridPoints, double dV);
void calcChiaiZetaMatrix(zomplex *ChiaiZeta, zomplex *RrsZeta, long nZeta, double *Cbs, lng_st ist, long elecHoleFlag);
void calcTZetaMatrix(zomplex *TZeta, zomplex *RckZeta, long nZeta, double *Cbs, lng_st ist);
void calcAaiMatrix(zomplex *Aai, zomplex *ChiaiZeta, zomplex *TZeta, long nZeta, long nIntExc, long nA, long nI);
long calcRandomCoulombStates(zomplex *thetaZeta, long nZeta, zomplex *potq, lng_st ist, par_st par, 
            fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwpsi);
void calcStochRealSpaceCoulomb(zomplex *thetaZeta, zomplex *potq, double *vx, double *vy, double *vz, lng_st ist, par_st par, 
                    fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwpsi);
void printCoulombAndRMatricesConvergenceInfo(zomplex *RckZeta, zomplex *RabZeta, lng_st ist);

// Functions in Hmat.c
double dotpreal(zomplex *psi, double *phi, long n, long m, long nGridPoints, double dv);
void Hmatreal(double *psitot, double *potl, double *ksqr, double *eval, lng_st ist, par_st par, 
              fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi);

// Functins in ortho.c
long portho(double *psi, double dv, lng_st ist);

// Functions in generate-filter-states.c
void generate_filter_states(double *psitot, double *eval, double *sige, double *evalbe, double *ksqr, double *potl,
                            double *vx, double *vy, double *vz, lng_st *ist, par_st *par, 
                            fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwpsi);

// Functions in coulomb.c
void calcAllCoulombMatrixElements(double *vijck, double *vabck, double *psiai, double *evalai,
            double *psibs, zomplex *potq, double *pothhl, lng_st ist, par_st par, 
            fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwpsi);
double calcOneCoulombMatrixElement(double *psiR, double *psiS, double *psiU, double *psiT, zomplex *potq,
            double dV, long nGrid, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi);

// Functions in nerror.c
void nerror(char *);

// Functions in write.c
void writeInputParameters(lng_st ist, par_st par, FILE *pf);
void writeFilterInfo(double *targetEnergies, lng_st ist, par_st par, FILE *pf);
void writeNonIntBiexciton(nonIntBiexc biExciton, FILE *pf);
void writeLongArray(long *longArray, long arrayLength, char *fileNameBase);
void writeDoubleArray(double *doubleArray, long arrayLength, char *fileNameBase);
void writeComplexArray(zomplex *zomplexArray, long arrayLength, char *fileNameBase);
void writeCurrentTime(FILE *pf);
void writeSeparation(FILE *pf);

/****************************************************************************/

#endif

/****************************************************************************/
