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

/****************************************************************************/

typedef struct st0 {
  double re, im;
} zomplex;

typedef fftw_plan fftw_plan_loc;

typedef struct st1 {
  double dx, dy, dz, dr, dkx, dky, dkz, dv;
  double xmin, xmax, ymin, ymax, zmin, zmax, kxmin, kymin, kzmin;
  double Ekinmax, gamma, gamma2, Elmin, Elmax, lumoEnergy, homoEnergy, Egap;
  double Vmin, Vmax, dE, dE_1, DeltaE;
  double deps, Ex;
  double Eamin, Eamax, Eimin, Eimax;  
  double temp, kbT, boltzEnergyRange, minInitE, maxInitE;
} par_st;

typedef struct st4 {
  long n1, n2, n12, natom;
  long ms, ms2, niter, nc, npot, npsi;
  long nx, ny, nz, ngrid, natomtype;
  long mfermi, ndos, nmc;
  long nthreads, two;
  long homoIndex, lumoIndex, totalHomo, totalLumo, msbs, msbs2;
  long numBandEdgeHoles, numBandEdgeElectrons, numBandEdgeStates;
  double nx_1, ny_1, nz_1, ngrid_1;
  long seed;
  long imin, imax, amin, amax, atot, itot;
  long nexcitons, numBiexcitons;
} lng_st;

typedef struct st9 {
  long natyp;
  char atyp[3];
} atm_st;

typedef struct st10 {
  double x0, t, e;
} pot_st;

/****************************************************************************/

#define sqr(x)       ((x) * (x))
#define cube(x)      ((x) * (x) * (x))
#define forth(x)     ((x) * (x) * (x) * (x))

#define ISWAP(a,b)   {double temp = (a); (a) = (b); (b) = -temp;}
#define ISWAP2(a,b)  {double temp = (a); (a) = 2.0 * (b); (b) = -2.0 * temp;}

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

long assign_atom_number(char atyp[2]);
void assign_atom_type(char *atype,long j);

void init(double *vx,double *vy,double *vz,double *ksqr,double *potl,double *rx,double *ry,double *rz,par_st *par,lng_st *ist);
void init_size(long, char *argv[],par_st *,lng_st *);
void init_pot(double *vx,double *vy,double *vz,zomplex *potq,par_st par,lng_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void init_psi(zomplex *psi,lng_st ist,par_st par,long *idum);
double expot(double r, pot_st ppar);

double screenedcoulomb(double dr,double gamma);
double interpolate(double r,double dr,double *vr,double *pot,long npot,long n,long j);

long read_conf(double *rx,double *ry,double *rz,atm_st *atm,long n,FILE *);

double norm(zomplex *, double,long);
double normalize(zomplex *,double,long);
void normalize_all(double *,double,long,long);
void normalize_all_omp(double *psi,double dr,long ms,long ngrid,long nthreads);

void write_psi(double *,double *,double *,double *,double *,lng_st,par_st);
void write_pot(double *,double *,double *,double *,lng_st);

void nerror(char *);

double rint(double); 
long random(void);
void srandom(unsigned int);
double ran_nrc(long *idum);
double ran();

double get_dot_ligand_size(double *,double *,double *,long);
double get_dot_ligand_size_z(double *rz,long n);

void gauss_test(double *vx,double *vy,double *vz,zomplex *potq,double *poth,par_st par,lng_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);

void hartree(zomplex *rho,zomplex *potq,double *poth,lng_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);

void hamiltonian(zomplex *phi,zomplex *psi,double *potl,double *ksqr,lng_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void kinetic(zomplex *psi,double *ksqr,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,lng_st ist);

void filtering(double *psitot,double *potl,double *ksqr,zomplex *an,double *zn,double *el,lng_st ist,par_st par,long tid,long jns,long *idum,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void filter(zomplex *psin,zomplex *psim1,double *psi0,double *potl,double *ksqr,par_st par,zomplex *an,double *zn,lng_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,long tid,long jns);

void hnorm(zomplex *psim1,zomplex *psin,double *potl,double *ksqr,par_st par,double zm1,lng_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void hamiltonian_norm(zomplex *psi,zomplex *phi,double *potl,double *ksqr,lng_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);

double energy(zomplex *psi,zomplex *phi,double *potl,double *ksqr,lng_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void energy_all(zomplex *psi,zomplex *phi,double *psi0,double *potl,double *ksqr,double *ene,lng_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,long ms);
void get_energy_range(double *vx,double *vy,double *vz,double *ksqr,double *potl,par_st *par,lng_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void calc_sigma_E(double *psitot,double *potl,double *ksqr,double *eval2,lng_st ist,par_st par);
long calc_bse_biexcitonic_states(double *eval, double *de, lng_st ist, par_st par);

void chebyshev_reordered(double *polong,double min,double max,long nc);
double samp_points_ashkenazy(zomplex *polong,double min,double max,long nc);

void calculate_auger(double *vkijb,double *vkcab,double *evalbe,double *evalai,double *sigeai,lng_st ist,par_st par);
void calc_boltzmann_weighted_rates(double *energies, double *eRate, double *hRate, double temp, int numStates);
double calc_partition_function(double *energies, double temp, int numStates);

double dotpreal(zomplex *psi,double *phi,long n,long m,long ngrid,double dv);
void Hmatreal(double *psitot,double *potl,double *ksqr,double *eval,lng_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
long portho(double *psi,double dv,lng_st ist);

void generate_filter_states(double *psitot,double *eval,double *sige,double *evalbe,double *ksqr,double *potl,double *vx,double *vy,double *vz,lng_st *ist,par_st *par,fftw_plan_loc *planfw,fftw_plan_loc *planbw,fftw_complex *fftwpsi);
void generate_coulomb_matrix_elements(double *vkijb,double *vkcab,double *psiai,double *evalai,double *sigeai,double *psibs,zomplex *potq,double *pothhl,lng_st ist,par_st par,fftw_plan_loc *planfw,fftw_plan_loc *planbw,fftw_complex *fftwps);

/****************************************************************************/
