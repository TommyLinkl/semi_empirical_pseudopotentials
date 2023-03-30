/*****************************************************************************/

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

/*****************************************************************************/

struct pzdata {
    double z;
    double pz;
};

typedef struct st0 {
  double re, im;
} zomplex;

typedef struct st1 {
  double dx, dy, dz, dr, dkx, dky, dkz, dv, epsx, epsy, epsz, boxl, minr;
  double xmin, xmax, ymin, ymax, zmin, zmax, kxmin, kymin, kzmin;
  double Ekinmax, gamma, gamma2, Elmin, Elmax, Elumo, Ehomo, Vmin, Vmax;
  double deltae, deltah, deps, fermiEnergy;
} par_st;

typedef struct st4 {
  long n1, n2, n12, natom, nthreads;
  long ms, ms2, niter, nc, npot, npsi;
  long nx, ny, nz, ngrid, natomtype;
  long nhomo, nlumo, totalhomo, totallumo;
  double nx_1, ny_1, nz_1, ngrid_1;
  long maxElecStates, maxHoleStates;
  long printFPDensity; // 0 = False (default) or 1 = True
  long calcDarkStates; // 0 = False (default) or 1 = True
} long_st;

typedef struct st9 {
  long natyp;
  char atyp[3];
} atm_st;

typedef fftw_plan fftw_plan_loc;

/*****************************************************************************/

#define rlong(x)      (rint(x)) 

#define sqr(x)       ((x) * (x))
#define cube(x)      ((x) * (x) * (x))
#define forth(x)     ((x) * (x) * (x) * (x))

#define ISWAP(a,b)   {double temp = (a); (a) = (b); (b) = -temp;}
#define ISWAP2(a,b)  {double temp = (a); (a) = 2.0 * (b); (b) = -2.0 * temp;}

#define AUTOEV    27.2114
#define PIE       3.14159265358979323846
#define TWOPI     6.28318530717958647692
#define FOURPI    (4.0*3.14159265358979323846)
#define SVDEPS    1.0e-10
#define EPSR      1.0e-10
#define EPSCHI    1.0e-8
#define DENE      1.0e-10

/*****************************************************************************/

/*#define DEPS  0.02
  #define DENERGY 0.01*/
long assign_atom_number(char atyp[2]);
void assign_atom_type(char *atype,long j);

void init(double *potl,double *vx,double *vy,double *vz,double *ksqr,double *rx,double *ry,double *rz,par_st *par,long_st *ist);
void init_size(long, char *argv[],par_st *,long_st *);
void init_pot(double *vx,double *vy,double *vz,zomplex *potq,zomplex *potqx,par_st par,long_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void init_psi(zomplex *psi,double *vx,double *vy,double *vz,long_st ist,par_st par,long *idum);
double screenedcoulomb(double dr, double gamma);
double longerpolate(double r,double dr,double *vr,double *pot,long npot,long n,long j);

void read_conf(double *rx,double *ry,double *rz,atm_st *atm,long n,FILE *);

double norm(zomplex *, double,long);
double normalize_zomplex(zomplex *psi, double dr, long ngrid);
void normalize(double *vector, double dV, long ngrid);
void normalize_all(double *,double,long,long);
void norm_vector(double *vector, double dV, long length);
double norm_rho(zomplex *rho,double dr,long ngrid);


void write_psi(double *,double *,double *,double *,double *,long_st,par_st);
void write_pot(double *,double *,double *,double *,long_st);

void scalar_product(zomplex *,zomplex *,zomplex *,double,long,long);

void nerror(char *);

double rlong(double); 
/*
long int random(void);
void srandom(unsigned long);
*/
double ran_nrc(long *idum);
void Randomize();
double ran();

double get_dot_ligand_size(double *,double *,double *,long);
double get_dot_ligand_size_z(double *rz,long n);

void hartree(zomplex *rho,zomplex *potq,double *poth,long_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void single_coulomb(double *psi,zomplex *potq,zomplex *potqx,double *poth,double *eval,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,double *bsmat,double *h0mat);



double energy(zomplex *psi,zomplex *phi,double *potl,double *ksqr,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
double energy_norm(zomplex *psi,zomplex *phi,double *potl,double *ksqr,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void energy_all(double *psi0,double *potl,double *ksqr,double *ene,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,long ms);
void get_energy_range(double *vx,double *vy,double *vz,double *ksqr,double *potl,par_st *par,long_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void calc_sigma_E(zomplex *psi,zomplex *phi,double *psitot,double *potl,double *ksqr,double *eval2,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);


void prlong_pz(double *psi,double *sige,double *vz,par_st par,long_st ist);

void diag(const long long n, int nthreads, double *mat, double *eval);
void bethe_salpeter(double *bsmat, double *h0mat, double *psi, double *vz, double *mux, double *muy, double * muz,
					double *mx, double *my, double *mz, long_st ist,par_st par);

void psi_rnd(zomplex *psi,long ngrid,double dv,long *idum);

double findmaxabsre(zomplex *dwmat,long n);
double findmaxabsim(zomplex *dwmat,long n);

void dipole(double *vx,double *vy,double *vz,double *psi,double *mux,double *muy,double *muz,double *eval,long_st ist,par_st par);
void mag_dipole(double *vx, double *vy, double *vz, double *psi, double *mx, double *my, double *mz, 
  double *eval, fftw_plan_loc *planfw,fftw_plan_loc *planbw,fftw_complex *fftwpsi, long_st ist, par_st par);
void rotational_strength(double *rs, double *mux, double *muy, double *muz, double *mx, 
  double *my, double *mz, double *eval, long_st ist);

void hamiltonian(zomplex *phi,zomplex *psi,double *potl,double *ksqr,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void kinetic(zomplex *psi,double *ksqr,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,long_st ist);

void prlong_cube(double *pgrid,long_st ist,par_st par);

void single_coulomb_openmp(double *psi,zomplex *potq,zomplex *potqx,double *poth,double *eval,long_st ist,par_st par,fftw_plan_loc *planfw,fftw_plan_loc *planbw,fftw_complex *fftwpsi,double *bsmat,double *h0mat);

void print_pz_one(double *psi,double *vz,par_st par,long_st ist,char *str);
void print_pz(double *psi,double *sige,double *vz,par_st par,long_st ist);
int z_project(double *vector, double *vz, par_st par, long_st ist, char *fname);
void print_cube(double *pgrid,long_st ist,par_st par,char *fName);
void print_fixed_qp_density(double *psi, double *Cbs, double *vz, long_st ist, par_st par);

// Functions that write input or output - write.c
void writeCurrentTime(FILE *pf);
void writeSeparation(FILE *pf);

/*****************************************************************************/
