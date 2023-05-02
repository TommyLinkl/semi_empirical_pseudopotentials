/****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include <time.h>
#include <time.h>
#include <assert.h>
#include <fftw3.h>
#include <omp.h>
#include "unistd.h"
#include "mkl.h"
#include "mkl_lapacke.h"

#include "vector.h"

/****************************************************************************/

typedef struct st0 {
  double re, im;
} zomplex;

typedef fftw_plan fftw_plan_loc;

typedef struct st1 {
  double dx, dy, dz, dr, dkx, dky, dkz, dv, dt;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double Vmin, Vmax, dE, dE_1, Elmin, Elmax, Ekinmax;
} par_st;

typedef struct st4 {
  long natom, lumo, homo;
  long ms, ns, nc, npot, mstot;
  long nx, ny, nz, ngrid;
  long natomtype;
  long nthreads;
  double nx_1, ny_1, nz_1;
  char crystalStructure[30];
  char outmostMaterial[30];
} long_st;

typedef struct st7 {
  double nr00, nr20, nr21, nr2_1, nr22, nr2_2;
} nrm_st;

typedef struct st5 {
  double x, y, z;
} xyz_st;

typedef struct st9 {
  long natyp;
  char atyp[3];
} atm_st;

typedef struct st10 {
  double x0,t,e;
} pot_st;

/****************************************************************************/

#define sqr(x)       ((x) * (x))
#define cube(x)      ((x) * (x) * (x))
#define forth(x)     ((x) * (x) * (x) * (x))

#define ISWAP(a,b)   {double temp = (a); (a) = (b); (b) = -temp;}
#define ISWAP2(a,b)  {double temp = (a); (a) = 2.0 * (b); (b) = -2.0 * temp;}

#define cplus(a,b,c)  {zomplex tmp; (tmp).re=(a).re+(b).re; (tmp).im=(a).im+(b).im; c=tmp;}
#define cminus(a,b,c) {zomplex tmp; (tmp).re=(a).re-(b).re; (tmp).im=(a).im-(b).im; c=tmp;}
#define cmul(a,b,c)   {zomplex tmp; (tmp).re=(a).re*(b).re-(a).im*(b).im; (tmp).im=(a).re*(b).im+(a).im*(b).re; c=tmp;}
#define cmuls(a,b,c)  {zomplex tmp; (tmp).re=(a).re*(b).re+(a).im*(b).im; (tmp).im=(a).re*(b).im-(a).im*(b).re; c=tmp;}
#define cdev(a,b,c)   {zomplex tmp; double mechane; mechane = (1.0 / ((b).re*(b).re+(b).im*(b).im)); (tmp).re=((a).re*(b).re+(a).im*(b).im) * mechane; (tmp).im = ((a).im * (b).re - (a).re*(b).im)*mechane; c=tmp;}
#define cexp(a,c)     {double texp = exp((a).re); zomplex tmp; (tmp).re=texp * cos((a).im); (tmp).im = texp * sin((a).im); c= tmp;}
#define cexpminx(a,c) {double texp = exp(-(a).re); (c).re=texp * cos((a).im); (c).im = -texp * sin((a).im);}

#define PIE       3.14159265358979323846
#define TWOPI     6.28318530717958647692
#define SVDEPS    1.0e-10
#define EPSR      1.0e-10
#define EPSCHI    1.0e-8
#define DENE      1.0e-10

/****************************************************************************/

void init(double *potl,double *vx,double *vy,double *vz,double *ksqr,double *rx,double *ry,double *rz,atm_st *atm,par_st *par,double *eval,long_st *ist,fftw_plan_loc *planfw,fftw_plan_loc *planbw,fftw_complex *fftwpsi);
void readNearestNeighbors(long nAtoms, int crystalStructure, vector *atomNeighbors, double *tetrahedronVolRef, int outmostMaterial);
void calculateStrainScale(long nAtoms, double *tetrahedronVolRef, atm_st *atm,
        vector *atomNeighbors, double *a4Params, double *a5Params, double *strainScale);
void init_psi(zomplex *psi,long_st ist,par_st par,long *idum);
void init_size(long, char *argv[],par_st *,long_st *);

long assign_atom_number(char atyp[2]);
void assign_atom_type(char *atype,long j);
double retIdealBondLength(long natyp_1, long natyp_2, int crystalStructureInt);

void read_conf(double *rx,double *ry,double *rz,atm_st *atm,long n,FILE *);
void read_pot(double *vr,double *pot,long *npot,double *dr,atm_st *atm,long n,long ntype, double *a4Params, double *a5Params);
double get_dot_ligand_size_z(double *rz,long n);
double interpolate(double r,double dr,double *vr,double *pot,long npot,long n,long j, double strainScaleFactor);

double ran_nrc(long *idum);
void Randomize();
double ran();

void kinetic(zomplex *psi,double *ksqr,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,long_st ist);
void hamiltonian(zomplex *phi,zomplex *psi,double *potl,double *ksqr,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);

void filtering(double *psitot,double *potl,double *ksqr,zomplex *an,double *zn,double *el,long_st ist,par_st par,long tid,long jns,long *idum);
void filter(zomplex *psin,zomplex *psim1,double *psi0,double *potl,double *ksqr,par_st par,zomplex *an,double *zn,long_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,long tid,long jns);
void hnorm(zomplex *psim1,zomplex *psin,double *potl,double *ksqr,par_st par,double zm1,long_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void hamiltonian_norm(zomplex *psi,zomplex *phi,double *potl,double *ksqr,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);


double norm(zomplex *, double,long,long);
double normd(double *, double,long,long);
double normalize(zomplex *,double,long,long);
void normalize_all(double *psi,double dr,long ms,long ngrid,long nthreads);

double energy(zomplex *psi,zomplex *phi,double *potl,double *ksqr,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void energy_all(zomplex *psi,zomplex *phi,double *psi0,double *potl,double *ksqr,double *ene,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void get_energy_range(zomplex *psi,zomplex *phi,double *potl,double *vx,double *vy,double *vz,double *ksqr,par_st *par,long_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);


void coefficient(zomplex *an,double *samp,double *El,par_st par,long_st ist);
void chebyshev_reordered(double *,double,double,long);
double samp_points_ashkenazy(zomplex *point,double min,double max,long nc);
void check_function(zomplex *an,zomplex *samp,long_st ist,par_st par,double El);

void reorder(long *,long);
long write_psi(double *,double *,double *,double *,double *,long,long,par_st);
void write_pot(double *,double *,double *,double *);
void writeCubeFile(double *rho, double xmin, double ymin, double zmin, double dx, double dy, double dz, long nx, long ny, long nz); 
void writeCurrentTime(FILE *pf);
void writeSeparation(FILE *pf);

void scalar_product(zomplex *,zomplex *,zomplex *,double,long,long);

void Hmatreal(zomplex *psi,zomplex *phi,double *psi0,double *potl,double *ksqr,double *eval,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
double dotpreal(zomplex *psi,double *phi,long n,long m,long ngrid,double dv);

void nerror(char *);

long ortho(double *,double,long_st);
long portho(double *,double,long_st);
//void dgesvd_(char *,char *,long *,long *, double *,long *,double *,void *,long *,void *,long *,double *,long *,long *);

double rint(double); 
/*long int random(void);
  void srandom(unsigned long);*/
//void dsyev_(char *,char *,long *,double *,long *,double *,double *,long *,long *);

void calc_sigma_E(zomplex *psi,zomplex *phi,double *psitot,double *potl,double *ksqr,double *eval,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);

xyz_st calc_dipole(long i,long a,double *vx,double *vy,double *vz,double *psitot,long_st ist,par_st par);

double expot(double r,pot_st ppar);
double retRegularTetrahedronVolume(double bondLength1, double bondLength2, double bondLength3, double bondLength4);
/****************************************************************************/
