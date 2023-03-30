/*****************************************************************************/
// Required libraries

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include <time.h>
#include <assert.h>
#include <fftw3.h>
#include <omp.h>
#include "mkl.h"

/*****************************************************************************/
// Application specific structures

typedef struct st0 {
  double re, im;
} zomplex;

typedef fftw_plan fftw_plan_loc;

typedef struct st1 {
  double dx, dy, dz, dr, dkx, dky, dkz, dv, dt;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double Vmin, Vmax, dE, dE_1, Elmin, Elmax, Ekinmax;
  double Rnlcut2, sigma, sigma_1;
} par_st;

typedef struct st4 {
  long natom, lumo, homo, flaghomo, nsemicond;
  long ms, ns, nc, npot, mstot, mstotngrid, nspinngrid;
  long nx, ny, nz, ngrid, nspin, nnlc;
  long natomtype;
  long nthreads;
  long flagkb;
  double nx_1, ny_1, nz_1;
} long_st;

typedef struct st5 {
  double x, y, z;
} xyz_st;

typedef struct st9 {
  long natyp;
  char atyp[3];
  double Vso, initE;
} atm_st;

typedef struct st11 {
  long jxyz;
  zomplex y10, y11, y1_1;
  double Vs0G1, Vs0G_2;
  double r, r2_1, r2, Vr;
} nlc_st;

typedef struct st12{
  zomplex G1y10up, G1y11up, G1y1_1up;
  zomplex G1y10dn, G1y11dn, G1y1_1dn;

  zomplex G_2y10up, G_2y11up, G_2y1_1up;
  zomplex G_2y10dn, G_2y11dn, G_2y1_1dn;
} ing_st;

/*****************************************************************************/
// Macro definitions

// double macros
#define sqr(x)       ((x) * (x))
#define cube(x)      ((x) * (x) * (x))
#define forth(x)     ((x) * (x) * (x) * (x))
#define ISWAP(a,b)   {double temp = (a); (a) = (b); (b) = -temp;}
#define ISWAP2(a,b)  {double temp = (a); (a) = 2.0 * (b); (b) = -2.0 * temp;}

// complex number macros
#define cplus(a,b,c)  {zomplex tmp; (tmp).re=(a).re+(b).re; (tmp).im=(a).im+(b).im; c=tmp;}
#define cminus(a,b,c) {zomplex tmp; (tmp).re=(a).re-(b).re; (tmp).im=(a).im-(b).im; c=tmp;}
#define cmul(a,b,c)   {zomplex tmp; (tmp).re=(a).re*(b).re-(a).im*(b).im; (tmp).im=(a).re*(b).im+(a).im*(b).re; c=tmp;}
#define cmuls(a,b,c)  {zomplex tmp; (tmp).re=(a).re*(b).re+(a).im*(b).im; (tmp).im=(a).re*(b).im-(a).im*(b).re; c=tmp;}
#define cdev(a,b,c)   {zomplex tmp; double mechane; mechane = (1.0 / ((b).re*(b).re+(b).im*(b).im)); (tmp).re=((a).re*(b).re+(a).im*(b).im) * mechane; (tmp).im = ((a).im * (b).re - (a).re*(b).im)*mechane; c=tmp;}
#define cexp(a,c)     {double texp = exp((a).re); zomplex tmp; (tmp).re=texp * cos((a).im); (tmp).im = texp * sin((a).im); c= tmp;}
#define cexpminx(a,c) {double texp = exp(-(a).re); (c).re=texp * cos((a).im); (c).im = -texp * sin((a).im);}

// constant macros
#define AUTOEV    27.211385
#define PIE       3.14159265358979323846
#define TWOPI     6.28318530717958647692
#define SVDEPS    1.0e-10
#define EPSR      1.0e-10
#define EPSE      1.0e-10
#define EPSR02	  1.0e-2
#define EPSR0     1.0e-4
#define EPSCHI    1.0e-8
#define DENE      1.0e-10
#define EPSDX     1.0e-20

/*****************************************************************************/
// Function declarations

void init(double *potl,double *vx,double *vy,double *vz,double *ksqr,double *rx,double *ry,double *rz,atm_st *atm,par_st *par,double *eval,long_st *ist,double *dr,double *vr,double *potatom,long *npot);
void init_psi(zomplex *psi,long_st ist,par_st par,long *idum);
void init_size(long, char *argv[],par_st *,long_st *);
void init_conf(double *rx,double *ry,double *rz,atm_st *atm,par_st *par,long_st *ist);
void init_list(nlc_st *nlc,long *nl,double *vx,double *vy,double *vz,double *rx,double *ry,double *rz,atm_st *atm,par_st par,long_st ist);
void init_kb(nlc_st *nlc,long *nl,double *Elkb,double *vx,double *vy,double *vz,double *rx,double *ry,double *rz,atm_st *atm,par_st par,long_st ist,double *dr,double *vr,double *potatom,long *npot);
void solve_for_ur(double *r,double *ur,double *pot,double dr,long nr,double kappa,double initE, char *str);

// General functions needed for numerical integration
// used in this program to solve the radial Dirac equation
void gen_numerov_solver(double *u, double *p, double *q, double *s, double *r, double u0, double u1, int istart, int iend);
void numerov_solver(double *y, double *k, double *b, double y0, double y1, double h, int istart, int iend);
void runge_kutta_solver(double *a, double *p, double *f, double *r, double a0, double h, int istart, int iend);

// Functions used to wrap dirac equation into general ODE formats
// to be used as inputs into the functions in num_int.c
// these functions are in init.c
void calc_P_rk_dirac(double *P, double *p, int istart, int iend);
double energy_corrector_pt(double *urIn, double *urOut, double dr, int jrMid);
double min_element(double *a, long len);
long count_nodes(double *a, long len);
double calc_average(double *a, int istart, int iend);


long assign_atom_number(char atyp[2]);
void assign_atom_type(char *atype,long j);
long get_number_of_atom_types(atm_st *atm,long_st ist,long *list);

long read_conf(double *rx,double *ry,double *rz,atm_st *atm,long n,FILE *);
void read_pot(double *vr,double *pot,long *npot,double *dr,atm_st *atm,long n,long ntype);
double get_dot_ligand_size_z(double *rz, long n);
double interpolate(double r,double dr,double *vr,double *pot,long npot,long n,long j);

double ran_nrc(long *idum);
void Randomize();
double ran();

void kinetic(zomplex *psi,double *ksqr,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,long_st ist);
void nonlocal_potential(zomplex *phi,zomplex *psi,long_st ist,par_st par,nlc_st *nlc,long *nl);
void nonlocal_potential_kb(zomplex *psinl,zomplex *psi,long_st ist,par_st par,nlc_st *nlc,long *nl,double *Elkb);
void hamiltonian(zomplex *phi,zomplex *psi,double *potl,double *ksqr,long_st ist,par_st par,nlc_st *nlc,long *nl,double *Elkb,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);

void filtering(zomplex *psi,double *potl,double *ksqr,zomplex *an,double *zn,double *el,long_st ist,par_st par,nlc_st *nlc,long *nl,double *Elkb,long tid,long jns);
void filter(zomplex *psin,zomplex *psim1,zomplex *psims,double *potl,double *ksqr,par_st par,nlc_st *nlc,long *nl,double *Elkb,zomplex *an,double *zn,long_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,long tid,long jns);
void hnorm(zomplex *psim1,zomplex *psin,double *potl,double *ksqr,par_st par,nlc_st *nlc,long *nl,double *Elkb,double zm1,long_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void calcAllSpinUpDownRatios(zomplex *psitot, long nGrid, long nStates);

double norm(zomplex *, double,long,long);
double normalize(zomplex *,double,long,long);
void normalize_all(zomplex *psi,double dr,long ms,long ngrid,long nthreads);

double energy(zomplex *psi,zomplex *phi,double *potl,double *ksqr,long_st ist,par_st par,nlc_st *nlc,long *nl,double *Elkb,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void energy_all(zomplex *psi,zomplex *phi,zomplex *psims,double *potl,double *ksqr,double *ene,long_st ist,par_st par,nlc_st *nlc,long *nl,double *Elkb,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void get_energy_range(zomplex *psi,zomplex *phi,double *potl,double *vx,double *vy,double *vz,double *ksqr,long_st ist,par_st *par,nlc_st *nlc,long *nl,double *Elkb,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void calc_sigma_E(zomplex *psi,zomplex *phi,zomplex *psitot,double *potl,double *ksqr,double *eval,long_st ist,par_st par,nlc_st *nlc,long *nl,double *Elkb,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);

void coefficient(zomplex *an,double *samp,double *Elkb,par_st par,long_st ist);
void chebyshev_reordered(double *,double,double,long);
double samp_points_ashkenazy(zomplex *point,double min,double max,long nc);
void check_function(zomplex *an,zomplex *samp,long_st ist,par_st par,double Elkb);

void reorder(long *,long);

void Hmat(zomplex *psi,zomplex *phi,MKL_Complex16 *psi0,double *potl,double *ksqr,double *eval,long_st ist,par_st par,nlc_st *nlc,long *nl,double *Elkb,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
MKL_Complex16 dotp(zomplex *psi,MKL_Complex16 *phi,long m,long ngrid,double dv);

void nerror(char *);

long ortho(double *,double,long_st);
long portho(MKL_Complex16 *,double,long_st);
//void dgesvd_(char *,char *,long *,long *, double *,long *,double *,void *,long *,void *,long *,double *,long *,long *);

double rint(double); 
/*long int random(void);
  void srandom(unsigned long);*/
//void dsyev_(char *,char *,long *,double *,long *,double *,double *,long *,long *);

/*****************************************************************************/
