#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include <time.h>
#include <assert.h>
#include <fftw3.h>
#include <omp.h>
#include "/opt/intel/compilers_and_libraries/linux/mkl/include/mkl.h"

/****************************************************************************/

typedef struct st0 {
  double re, im;
} zomplex;

//typedef fftw_complex zomplex;
//typedef fftwnd_plan fftw_plan_loc;
typedef fftw_plan fftw_plan_loc;

typedef struct st1 {
  double dx, dy, dz, dr, dkx, dky, dkz, dv, dt;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double Vmin, Vmax, dE, dE_1, Elmin, Elmax, Ekinmax;
  double electronMin, holeMax, holeFilts;
  long hprint, eprint;
  double sigma_e;
} par_st;

typedef struct st4 {
  long natom, lumo, homo, flaghomo;
  long ms, ns, nc, npot, mstot;
  long nx, ny, nz, ngrid;
  long natomtype;
  long nthreads;
  double nx_1, ny_1, nz_1;
} long_st;

typedef struct st7 {
  double nr00, nr20, nr21, nr2_1, nr22, nr2_2;
} nrm_st;

// typedef struct st5 {
//   double x, y, z;
// } xyz_st;

typedef struct st9 {
  long natyp;
  char atyp[3];
} atm_st;

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
#define AUTOEV    27.2114
#define AUTONS    2.418884254e-8
#define SVDEPS    1.0e-10
#define EPSR      1.0e-10
#define EPSCHI    1.0e-8
#define DENE      1.0e-10
#define DEPS      1.0e-3
#define DIM       3

/****************************************************************************/

// init.c functions
void init(double *potl,double *vx,double *vy,double *vz,double *ksqr,double *rx,double *ry,double *rz,atm_st *atm,par_st *par,double *eval,long_st *ist,fftw_plan_loc *planfw,fftw_plan_loc *planbw,fftw_complex *fftwpsi);
void init_psi(zomplex *psi,long_st ist,par_st par,long *idum);
void init_size(long, char *argv[],atm_st **atm,par_st *,long_st *,double **rx,double **ry,double **rz);
int rndToThreads(int gridpts, int threads);


// read.c functions
long assign_atom_number(char atyp[2]);
void assign_atom_type(char *atype,long j);
void read_conf(double *rx,double *ry,double *rz,atm_st *atm,long n,FILE *);
void read_pot(double *vr,double *pot,long *npot,double *dr,atm_st *atm,long n,long ntype);

// size.c functions
double get_dot_ligand_size_z(double *rz,long n);
double get_dot_ligand_size(double *rx,double *ry,double *rz,long n);

// interpolate.c functions
double interpolate(double r,double dr,double *vr,double *pot,long npot,long n,long j);

// hamiltonian.c functions
void kinetic(zomplex *psi,double *ksqr,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,long_st ist);
void hamiltonian(zomplex *phi,zomplex *psi,double *potl,double *ksqr,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);

// filter.c functions
void filtering(double *psitot,double *potl,double *ksqr,zomplex *an,zomplex *zn,double *el,long_st ist,par_st par,long tid,long jns,long *idum);
void filter(zomplex *psin,zomplex *psim1,double *psi0,double *potl,double *ksqr,par_st par,zomplex *an,zomplex *zn,long_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,long tid,long jns);

// hnorm.c functions
void hnorm(zomplex *psim1,zomplex *psin,double *potl,double *ksqr,par_st par,double zm1,long_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void hamiltonian_norm(zomplex *psi,zomplex *phi,double *potl,double *ksqr,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);

// norm.c functions
double norm(zomplex *, double,long,long);
double normalize(zomplex *,double,long,long);
void normalize_all(double *psi,double dr,long ms,long ngrid,long nthreads);

// energy.c functions
double energy(zomplex *psi,zomplex *phi,double *potl,double *ksqr,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void energy_all(zomplex *psi,zomplex *phi,double *psi0,double *potl,double *ksqr,double *ene,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void get_energy_range(zomplex *psi,zomplex *phi,double *potl,double *vx,double *vy,double *vz,double *ksqr,par_st *par,long_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void calc_sigma_E(zomplex *psi,zomplex *phi,double *psitot,double *potl,double *ksqr,double *eval,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);

// coeff.c functions
void coefficient(zomplex *an,zomplex *sampoints,double *target,long ncheby,long ntarget,long nthreads,double dE,double Emin,double dt);
double samp_points(zomplex *point,double min,double max,long nc);

// Hmat.c functions
void Hmatreal(zomplex *psi,zomplex *phi,double *psi0,double *potl,double *ksqr,double *eval,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
double dotpreal(zomplex *psi,double *phi,long n,long m,long ngrid,double dv);

// ortho.c functions
long ortho(double *,double,long_st);

// absorption.c functions
void absorption_spectrum(double *vx, double *vy, double *vz, double *psi, double *eval, double *sige, long_st ist, par_st par);
void calc_cd_spectrum(double *mux, double *muy, double *muz, double *eval, long *eigenstate_index_list, long num_holes, long num_excitons, long_st ist);
void fill_direction_averaged_array(zomplex *ave_array);

// states.c functions
long num_eigenstates_energy_range(double *eval, double *sige, double eMin, double eMax, long maxNumStates);
long get_eigenstate_list(long *eigenstate_index_list, double *eval, double *sige, double eMin, double eMax, long maxNumStates);

// rand.c functions
double ran_nrc(long *idum);
void Randomize();
double ran();

// nerror.c functions
void nerror(char *);

// ortho.c functions
long portho(double *psi,double dv,long_st ist);

// unused/ old functions
// double rint(double); 
// void reorder(long *,long);
// long write_psi(double *,double *,double *,double *,double *,long,long,par_st);
// void write_pot(double *,double *,double *,double *);
// long int random(void);
// void srandom(unsigned long);
// void dsyev_(char *,char *,long *,double *,long *,double *,double *,long *,long *);
// void dgesvd_(char *,char *,long *,long *, double *,long *,double *,void *,long *,void *,long *,double *,long *,long *);

/****************************************************************************/
