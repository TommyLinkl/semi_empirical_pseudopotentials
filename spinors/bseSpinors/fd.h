/****************************************************************************/
// Required libraries

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include <assert.h>
#include <fftw3.h>
#include <omp.h>
#include "mkl.h"

/****************************************************************************/
// Application specific structures

typedef struct st0 {
  double re, im;
} zomplex;

typedef struct st1 {
  double dx, dy, dz, dr, dkx, dky, dkz, dv, epsx, epsy, epsz, boxl, minr;
  double xmin, xmax, ymin, ymax, zmin, zmax, kxmin, kymin, kzmin;
  double Ekinmax, gamma, gamma2, Elmin, Elmax, Elumo, Ehomo, Vmin, Vmax;
  double deltae, deltah, deps;
} par_st;

typedef struct st4 {
  long n1, n2, n12, natom, nthreads;
  long ms, ms2, niter, nc, npot, npsi;
  long nx, ny, nz, ngrid, natomtype;
  long nspin, nspinngrid;
  long nhomo, nlumo, totalhomo, totallumo;
  double nx_1, ny_1, nz_1, ngrid_1;
} lng_st;

typedef struct st9 {
  long natyp;
  char atyp[3];
} atm_st;

typedef fftw_plan fftw_plan_loc;

/****************************************************************************/
// Macro definitions

#define rlong(x)      (rint(x)) 
#define sqr(x)       ((x) * (x))
#define cube(x)      ((x) * (x) * (x))
#define forth(x)     ((x) * (x) * (x) * (x))
#define ISWAP(a,b)   {double temp = (a); (a) = (b); (b) = -temp;}
#define ISWAP2(a,b)  {double temp = (a); (a) = 2.0 * (b); (b) = -2.0 * temp;}
#define AUTOEV       27.2114
#define PIE       3.14159265358979323846
#define TWOPI     6.28318530717958647692
#define FOURPI    (4.0*3.14159265358979323846)
#define SVDEPS    1.0e-10
#define EPSR      1.0e-10
#define EPSCHI    1.0e-8
#define DENE      1.0e-10

/****************************************************************************/
// Function declarations

// initialization functions
void init(double *vx,double *vy,double *vz,double *ksqr,par_st *par,lng_st *ist);
void init_size(long, char *argv[],par_st *,lng_st *);
void init_conf(double *rx,double *ry,double *rz,par_st *par,lng_st *ist);
void init_pot(double *vx,double *vy,double *vz,zomplex *potq,zomplex *potqx,par_st par,lng_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
double screenedcoulomb(double dr,double gamma);
void read_conf(double *rx,double *ry,double *rz,atm_st *atm,long n,FILE *);
long assign_atom_number(char atyp[2]);
double get_dot_ligand_size(double *rx, double *ry, double *rz, long n);
double get_dot_ligand_size_z(double *rz,long n);
void normalize_all(zomplex *psi, double dr, long numStates, long ngrid);
double normalize(zomplex *psi, double dr, long ngrid);
double norm(zomplex *psi, double dr, long ngrid);

// main computational work functions
void hartree(zomplex *rho,zomplex *potq,double *poth,lng_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void single_coulomb(zomplex *psi,zomplex *potq,zomplex *potqx,double *poth,double *eval,lng_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,double *bsmat,double *h0mat);
void hamiltonian(zomplex *phi,zomplex *psi,double *potl,double *ksqr,lng_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void kinetic(zomplex *psi,double *ksqr,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,lng_st ist);
void single_coulomb_openmp(zomplex *psi,zomplex *potq,zomplex *potqx,double *poth,double *eval,lng_st ist,par_st par,fftw_plan_loc *planfw,fftw_plan_loc *planbw,fftw_complex *fftwpsi,double *bsmat,double *h0mat);
void bethe_salpeter(double *bsmat,double *h0mat,zomplex *psi,double *vz,zomplex *mux,zomplex *muy,zomplex *muz,lng_st ist,par_st par);
void diag(int n,int nthreads,double *mat,double *eval);
void dipole(double *vx,double *vy,double *vz,zomplex *psi,zomplex *mux,zomplex *muy,zomplex *muz,double *eval,lng_st ist,par_st par);

// printing functions
void print_sp_pz(zomplex *psi, double *sige, double *vz, par_st par, lng_st ist);
void print_pz_one(double *density, double *densityUp, double *densityDown, double *vz,par_st par,lng_st ist,char *str);
void print_cube(double *pgrid,lng_st ist,par_st par,char *fName);

// general, simple complex number functions
zomplex sumComplex(zomplex c1, zomplex c2);
zomplex multiplyComplex(zomplex c1, zomplex c2);
double normComplex(zomplex c);
zomplex compConjComplex(zomplex c);

// general functions
double ran_nrc(long *idum);
void nerror(char *str);

/****************************************************************************/
