#include "fd.h"

/*****************************************************************************/

void init_size(long argc, char *argv[], par_st *par, lng_st *ist)
{
  FILE *pf;
 
  // read input and store parameters in ist and par structures
  pf = fopen("input.par" , "r");
  // TODO: allow for the xmin, ymin and zmin to be set in input.par
  fscanf(pf,"%ld",&ist->nx);  /*** number of grid points in x ***/
  fscanf(pf,"%ld",&ist->ny);  /*** number of grid points in y ***/
  fscanf(pf,"%ld",&ist->nz);  /*** number of grid points in z ***/
  fscanf(pf,"%lg",&par->dx);  /*** minimum grid spacing allowed ***/
  fscanf(pf,"%lg %lg %lg", &par->xmin, &par->ymin, &par->zmin);  /*** minimum x, y and z values of box ***/
  fscanf(pf,"%lg %lg %lg",&par->epsx,&par->epsy,&par->epsz); 
  fscanf(pf,"%lg %lg %lg",&par->deltae,&par->deltah,&par->deps); 
  fscanf(pf,"%ld",&ist->nthreads);  /*** number of threads as program uses openmp ***/
  fclose(pf);

  // read in the number of atoms in the system
  pf = fopen("conf.par" , "r");
  fscanf(pf,"%ld",&ist->natom);
  fclose(pf);

  par->Ekinmax = 20.0; /*** highest allowed value of ksqr ***/
  ist->npot = 8192; /*** length of pseudopotential files ***/

  par->dy = par->dz = par->dx;

  return;
}

/*****************************************************************************/

void init_conf(double *rx,double *ry,double *rz,par_st *par,lng_st *ist)
{
  FILE *pf;
  long ieof, i, a, nval, ntot, ntmp;
  double evalloc, deloc, *eval, *de;
  atm_st *atm;
  double xd, yd, zd;
  
  /*** read the pasivated nanocrystal configuration ***/
  pf = fopen("conf.par" , "r");
  fscanf(pf,"%ld",&ntot);
  assert(fabs((double)(ntot - ist->natom)) < 1.0e-15);
  if ((atm = (atm_st*)calloc(ist->natom,sizeof(atm_st)))==NULL)nerror("atm");
  read_conf(rx,ry,rz,atm,ntot,pf);
  fclose (pf);
  
  // Added July 2nd 2019 -> can now set xmin, ymin and zmin in input.par, set to negative to not use it
  if (par->xmin < EPS) {
    xd = rint(0.5 * get_dot_ligand_size_z(rx,ntot) + 5.0);
    yd = rint(0.5 * get_dot_ligand_size_z(ry,ntot) + 5.0);
    zd = rint(0.5 * get_dot_ligand_size_z(rz,ntot) + 5.0);
  }
  else {
    xd = par->xmin;
    yd = par->ymin;
    zd = par->zmin;
  }
  fprintf(stdout, "xd = %g yd = %g zd = %g\n", xd, yd, zd); fflush(stdout);

  /***initial parameters for the pot reduce mass, etc. in the x direction ***/
  par->xmin = -xd;
  par->xmax = xd;
  ntmp  = (long)((par->xmax - par->xmin) / par->dx);
  if (ntmp > ist->nx) ist->nx = ntmp;
  par->xmin = -((double)(ist->nx) * par->dx) / 2.0;
  par->xmax = ((double)(ist->nx) * par->dx) / 2.0;
  /*par->dx  = (par->xmax - par->xmin) / (double)(ist->nx);*/
  par->dkx = TWOPI / ((double)ist->nx * par->dx);
  
  /***initial parameters for the pot reduce mass, etc. in the y direction ***/
  par->ymin = -yd;
  par->ymax = yd;
  ntmp  = (long)((par->ymax - par->ymin) / par->dy);
  if (ntmp > ist->ny) ist->ny = ntmp;
  /*par->dy  = (par->ymax - par->ymin) / (double)(ist->ny);*/
  par->ymin = -((double)(ist->ny) * par->dy) / 2.0;
  par->ymax = ((double)(ist->ny) * par->dy) / 2.0;
  par->dky = TWOPI / ((double)ist->ny * par->dy);

  /***initial parameters for the pot reduce mass, etc. in the z direction ***/
  par->zmin = -zd;
  par->zmax = zd;
  ntmp  = (long)((par->zmax - par->zmin) / par->dz) - 1;
  if (ntmp > ist->nz) ist->nz = ntmp;
  /*par->dz  = (par->zmax - par->zmin) / (double)(ist->nz);*/
  par->zmin = -((double)(ist->nz) * par->dz) / 2.0;
  par->zmax = ((double)(ist->nz) * par->dz) / 2.0;
  par->dkz = TWOPI / ((double)ist->nz * par->dz);
  printf("new xd = %g yd = %g zd = %g\n",par->xmax,par->ymax,par->zmax);
  printf("nx = %ld  ny = %ld  nz = %ld npot = %ld\n", ist->nx,ist->ny,ist->nz,ist->npot);

  ist->nx_1 = 1.0 / (double)(ist->nx);
  ist->ny_1 = 1.0 / (double)(ist->ny);
  ist->nz_1 = 1.0 / (double)(ist->nz);
  ist->ngrid = ist->nx * ist->ny * ist->nz;
  ist->ngrid_1 = 1.0 / (double)(ist->ngrid);

  ist->nspin = 2;
  ist->natomtype = 15;
  ist->nspinngrid = ist->nspin*ist->ngrid;

  // make sure dielectric constant values are nonnegative
  if (par->epsx < 0.0) nerror("error in epsilonx");
  if (par->epsy < 0.0) nerror("error in epsilony");
  if (par->epsz < 0.0) nerror("error in epsilonz");

  // calculate total # of hole (totalhomo) and electron (totallumo) states 
  // and the indices of the homo (nhomo) and lumo (nlumo) states
  ist->nhomo = ist->nlumo = ist->totalhomo = ist->totallumo = 0;
  pf = fopen("eval.par" , "r");
  for (i = ieof = 0; ieof != EOF; i++){
    ieof = fscanf(pf,"%ld %lg %lg",&a,&evalloc,&deloc);
    if (deloc < par->deps && evalloc < -0.180) ist->nhomo = i;
    if (deloc < par->deps && evalloc < -0.180) ist->totalhomo++; 
    if (deloc < par->deps && evalloc > -0.180) ist->totallumo++; 
  }
  fclose(pf);
  nval = i - 1;
  pf = fopen("eval.par" , "r");
  for (i = 0; i <= ist->nhomo; i++) fscanf(pf,"%ld %lg %lg",&a,&evalloc,&deloc);
  for (i = ist->nhomo+1; i < nval; i++) {
    fscanf(pf,"%ld %lg %lg",&a,&evalloc,&deloc);
    if (deloc < par->deps) {
      ist->nlumo = i;
      break;
    }
  }
  fclose(pf);
  printf("The total number of hole eigenstates = %ld\n",ist->totalhomo);
  printf("The total number of electron eigenstates = %ld\n",ist->totallumo);
  printf("The index of the homo state = %ld\n",ist->nhomo);
  printf("The index of the lumo state = %ld\n",ist->nlumo);
  printf("nspinngrid = %ld\n", ist->nspinngrid);

  // store all energies and sigmas from the filter calculation
  // in eval and de, respectively
  eval = (double*)calloc(nval,sizeof(double)); 
  de = (double*)calloc(nval,sizeof(double)); 
  pf = fopen("eval.par" , "r");
  for (i = 0; i < nval; i++){
    fscanf(pf,"%ld %lg %lg",&a,&eval[i],&de[i]);
  }
  fclose(pf);
  printf("The homo energy = %g\n",eval[ist->nhomo]);
  printf("The lumo energy = %g\n",eval[ist->nlumo]);

  // check to see if the hole and electron states are within
  // deltah and deltae of the homo and lumo states, respectively
  for (ist->totalhomo = 0, i = ist->nhomo; i >= 0; i--)
    if ((eval[ist->nhomo] - eval[i] <= par->deltah) && de[i] < par->deps) ist->totalhomo++;
  for (ist->totallumo = 0, i = ist->nlumo; i < nval; i++)
    if ((eval[i] - eval[ist->nlumo] <= par->deltae) && de[i] < par->deps) ist->totallumo++;
  printf("The # of hole states within deltah of the homo energy = %ld\n", ist->totalhomo);
  printf("The # of elec states within deltae of the lumo energy = %ld\n", ist->totallumo);

  free(eval);  free(de);

  ist->ms = ist->totalhomo + ist->totallumo;
  printf("deltae = %g deltah = %g\ndeps = %g\n", par->deltae, par->deltah, par->deps);
  printf("nx = %ld  ny = %ld  nz = %ld npot = %ld\n", ist->nx,ist->ny,ist->nz,ist->npot);
  printf("ms = %ld  natom = %ld\n",ist->ms,ist->natom);
  printf("Ekinmax = %g\n",par->Ekinmax);
  printf("Epsilonx = %g epsilony = %g epsilonz = %g\n", par->epsx, par->epsy, par->epsz);

  free(atm);

  return;
}

/*****************************************************************************/

void init(double *vx,double *vy,double *vz,double *ksqr,par_st *par,lng_st *ist)
{
  long jx, jy, jz, jyz, jxyz, ntot, ntmp;
  double mx, my, mz, dx, dy, dz, *ksqrx, *ksqry, *ksqrz;

  if ((ksqrx  = (double*)calloc(ist->nx,sizeof(double)))==NULL)nerror("ksqrx");
  if ((ksqry  = (double*)calloc(ist->ny,sizeof(double)))==NULL)nerror("ksqry");
  if ((ksqrz  = (double*)calloc(ist->nz,sizeof(double)))==NULL)nerror("ksqrz");

  if ((par->xmax < par->ymax) && (par->xmax < par->zmax)) {
    par->minr = par->xmax;
    par->boxl = (double)(ist->nx) * par->dx;
  }
  else if ((par->ymax < par->xmax) && (par->ymax < par->zmax)) {
    par->minr = par->ymax;
    par->boxl = (double)(ist->ny) * par->dy;
  }
  else {
    par->minr = par->zmax;
    par->boxl = (double)(ist->nz) * par->dz;
  }    
  
  par->gamma = 7.0 / (2.0 * par->minr);
  par->gamma2 = sqr(par->gamma);
  
  par->Vmin = 1.0e10;
  par->Vmax = -1.0e10;

  par->dv = par->dx * par->dy * par->dz;
  par->dr = sqrt(sqr(par->dx) + sqr(par->dy) + sqr(par->dz));
  printf("dx = %g dy = %g dz = %g dv = %g dr = %g\n",
	  par->dx,par->dy,par->dz,par->dv,par->dr);

  /***initializing the ksqr vectors ***/
  for (ksqrx[0] = 0.0, jx = 1; jx <= ist->nx / 2; jx++)
    ksqrx[jx] = (ksqrx[ist->nx-jx] = 0.5 * sqr((double)(jx) * par->dkx) *
		ist->nx_1 * ist->ny_1 * ist->nz_1 / mx);

  for (ksqry[0] = 0.0, jy = 1; jy <= ist->ny / 2; jy++)
    ksqry[jy] = (ksqry[ist->ny-jy] = 0.5 * sqr((double)(jy) * par->dky) *
		ist->ny_1 * ist->nx_1 * ist->nz_1 / my);

  for (ksqrz[0] = 0.0, jz = 1; jz <= ist->nz / 2; jz++)
    ksqrz[jz] = (ksqrz[ist->nz-jz] = 0.5 * sqr((double)(jz) * par->dkz) *
		ist->nz_1 * ist->nx_1 * ist->ny_1 / mz);

  par->Ekinmax *= (ist->ny_1 * ist->nx_1 * ist->nz_1);
  for (jz = 0; jz < ist->nz; jz++) for (jy = 0; jy < ist->ny; jy++) {
    for (jyz = ist->nx * (ist->ny * jz + jy), jx = 0; jx < ist->nx; jx++) {
      jxyz = jyz + jx;
      ksqr[jxyz] = ksqrx[jx] + ksqry[jy] + ksqrz[jz];
      if (ksqr[jxyz] > par->Ekinmax) ksqr[jxyz] = par->Ekinmax;
    }
  }
  free(ksqrx); free(ksqry);  free(ksqrz);
  
  /***initializing the potential vector  ***/
  for (jx = 0, dx = par->xmin; jx < ist->nx; jx++, dx += par->dx) vx[jx] = dx;
  for (jy = 0, dy = par->ymin; jy < ist->ny; jy++, dy += par->dy) vy[jy] = dy;
  for (jz = 0, dz = par->zmin; jz < ist->nz; jz++, dz += par->dz) vz[jz] = dz;

  return;
}


/****************************************************************************/
// Calculates W(r-r') and stores it in potqx and V(r-r') and stores it in potq

void init_pot(double *vx,double *vy,double *vz,zomplex *potq,zomplex *potqx,par_st par,lng_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  long jx, jy, jz, jyz, jxyz, sx, sy, sz;
  double dr, dr2, x2, y2, z2, *kx2, *ky2, *kz2, alpha, cosa, sina;
  double ex_1, ey_1, ez_1, sqrtexeyez_1, sqrtaveps, boxl2 = sqr(par.boxl), sqrk0;
  double gammaeps;
  zomplex *potr, *potrx, tmp;

  ex_1 = 1.0 / par.epsx;
  ey_1 = 1.0 / par.epsy;
  ez_1 = 1.0 / par.epsz;
  sqrtexeyez_1 = 1.0 / sqrt(par.epsx * par.epsy * par.epsz);
  sqrtaveps = sqrt((par.epsx + par.epsy + par.epsz) / 3.0);
  gammaeps = par.gamma * sqrtaveps;

  /*** no yukawa screening for the exchange ***/
  sqrk0 = par.gamma2 * sqr(sqrtaveps);
  
  if ((kx2  = (double*)calloc(ist.nx,sizeof(double)))==NULL)nerror("kx2");
  if ((ky2  = (double*)calloc(ist.ny,sizeof(double)))==NULL)nerror("ky2");
  if ((kz2  = (double*)calloc(ist.nz,sizeof(double)))==NULL)nerror("kz2");
  if ((potr = (zomplex*)calloc(ist.ngrid,sizeof(zomplex)))==NULL)nerror("potr");
  if ((potrx = (zomplex*)calloc(ist.ngrid,sizeof(zomplex)))==NULL)nerror("potrx");
  
  for (jxyz = 0; jxyz < ist.ngrid; jxyz++) potr[jxyz].re = potr[jxyz].im = 0.0;
  for (jxyz = 0; jxyz < ist.ngrid; jxyz++) potrx[jxyz].re = potrx[jxyz].im = 0.0;

  for (jz = 0; jz < ist.nz; jz++) {
    z2 =sqr(vz[jz]);
    for (jy = 0; jy < ist.ny; jy++) {
      y2 = sqr(vy[jy]);
      jyz = ist.nx * (ist.ny * jz + jy);
      for (jx = 0; jx < ist.nx; jx++) {
      	x2 = sqr(vx[jx]);
      	jxyz = jyz + jx;
      	dr2 = (x2 + y2 + z2);
      	if (dr2 < boxl2) {
      	  dr = sqrt(x2 + y2 + z2);
      	  potr[jxyz].re = screenedcoulomb(dr,par.gamma);

      	  dr = sqrt(ex_1 * x2 + ey_1 * y2 + ez_1 * z2);
      	  potrx[jxyz].re = sqrtexeyez_1 * screenedcoulomb(dr,gammaeps);
      	}
      }
    }
  }

  for (jxyz = 0; jxyz < ist.ngrid; jxyz++) 
    potq[jxyz].re = potq[jxyz].im = potqx[jxyz].re = potqx[jxyz].im = 0.0;
  memcpy(&fftwpsi[0],&potr[0],ist.ngrid*sizeof(fftwpsi[0]));
  fftw_execute(planfw);
  memcpy(&potq[0],&fftwpsi[0],ist.ngrid*sizeof(potq[0]));

  memcpy(&fftwpsi[0],&potrx[0],ist.ngrid*sizeof(fftwpsi[0]));
  fftw_execute(planfw);
  memcpy(&potqx[0],&fftwpsi[0],ist.ngrid*sizeof(potqx[0]));
  
  for (kx2[0] = 0.0, jx = 1; jx <= ist.nx / 2; jx++)
    kx2[jx] = (kx2[ist.nx-jx] = sqr((double)(jx) * par.dkx));
  for (ky2[0] = 0.0, jy = 1; jy <= ist.ny / 2; jy++)
    ky2[jy] = (ky2[ist.ny-jy] = sqr((double)(jy) * par.dky));
  for (kz2[0] = 0.0, jz = 1; jz <= ist.nz / 2; jz++)
    kz2[jz] = (kz2[ist.nz-jz] = sqr((double)(jz) * par.dkz));

  for (sz = 1.0, jz = 0; jz < ist.nz; jz++, sz = -sz) {
    z2 = kz2[jz];
    for (sy = 1.0, jy = 0; jy < ist.ny; jy++, sy = -sy) {
      y2 = ky2[jy];
      jyz = ist.nx * (ist.ny * jz + jy);
      for (sx = 1.0, jx = 0; jx < ist.nx; jx++, sx = -sx) {
      	x2 = kx2[jx];
      	jxyz = jyz + jx;
      	alpha = PIE * (double)(jx + jy + jz + ist.ngrid / 2);
      	cosa = cos(alpha);
      	sina = sin(alpha);

      	/*** hatree term ***/
      	tmp.re = potq[jxyz].re;
      	tmp.im = potq[jxyz].im;
      	
      	potq[jxyz].re = (tmp.re * cosa - tmp.im * sina) * par.dv;
      	potq[jxyz].im = (tmp.re * sina + tmp.im * cosa) * par.dv;
      	potq[jxyz].re += (FOURPI / (x2 + y2 + z2 + par.gamma2));
      	potq[jxyz].re *= ist.ngrid_1;
      	potq[jxyz].im *= ist.ngrid_1;

      	/*** screened term ***/
      	tmp.re = potqx[jxyz].re;
      	tmp.im = potqx[jxyz].im;

      	potqx[jxyz].re = (tmp.re * cosa - tmp.im * sina) * par.dv;
      	potqx[jxyz].im = (tmp.re * sina + tmp.im * cosa) * par.dv;
      	potqx[jxyz].re += FOURPI / (par.epsx * x2 + par.epsy * y2 + par.epsz * z2 + sqrk0);
      	potqx[jxyz].re *= ist.ngrid_1;
      	potqx[jxyz].im *= ist.ngrid_1;
      }
    }
  }

  free(potr);  free(potrx); free(kx2); free(ky2); free(kz2);

  return;
}
  
/************************************************************************/

double screenedcoulomb(double dr, double gamma)
{
  if (dr < EPSR) return (gamma);
  return ((1.0 - exp(-gamma * dr)) / dr);
}

/*****************************************************************************/
