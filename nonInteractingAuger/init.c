#include "fd.h"

/*****************************************************************************/

void init_size(long argc, char *argv[], par_st *par, lng_st *ist)
{
  FILE* pf = fopen("input.par" , "r");
  double evalloc, deloc, *eval, *de, a;
  int i, ieof, nval;
  
  fscanf (pf,"%ld",&ist->nx);  /*** number of grid point in x ***/
  fscanf (pf,"%ld",&ist->ny);  /*** number of grid point in y ***/
  fscanf (pf,"%ld",&ist->nz);  /*** number of grid point in z ***/
  fscanf (pf,"%ld",&ist->nmc);  /*** number of iteration ***/
  fscanf (pf,"%ld",&ist->nc);    /*** length of newton interpolation ***/
  fscanf (pf,"%lg",&par->DeltaE); /*** DeltaE window used in enforcing energy conservation ***/
  fscanf (pf,"%ld",&ist->seed);    
  fscanf (pf,"%ld",&ist->nthreads);
  fscanf (pf,"%lg",&par->deps); /*** are you an eigenstate? ***/
  fscanf (pf,"%lg",&par->temp); /*** temperature for boltzmann averaging over initial states ***/
  fclose(pf);

  pf = fopen("conf.par" , "r");
  fscanf(pf,"%ld",&ist->natom);
  fclose(pf);
  
  /*** set default parameters - no longer specified in input.par - John ***/
  ist->npot = 8192;	/*** size of pseudopotential files ***/
  par->Elmin = -1;	/*** minimum energy ***/
  par->Elmax = 1;	/*** maximum energy ***/
  par->Ekinmax = 10;	/*** kinetic energy cutoff ***/
  ist->two = 8;
  ist->ms = ist->two * ist->nmc;	/*** number of states per filter ***/
  ist->nexcitons = 10;
  par->kbT = KB*par->temp/AUTOEV;
  par->boltzEnergyRange = 3.0*par->kbT;

  ist->nx_1 = 1.0 / (double)(ist->nx);
  ist->ny_1 = 1.0 / (double)(ist->ny);
  ist->nz_1 = 1.0 / (double)(ist->nz);
  ist->ngrid = ist->nx * ist->ny * ist->nz;
  ist->ngrid_1 = 1.0 / (double)(ist->ngrid);
  ist->natomtype = 15;
  
  printf("nx = %ld  ny = %ld  nz = %ld npot = %ld\n",
	  ist->nx, ist->ny, ist->nz, ist->npot);
  printf("ms = %ld  nmc = %ld natom = %ld nc = %ld nthreads = %ld\n",
	  ist->ms, ist->nmc, ist->natom, ist->nc, ist->nthreads);
  printf("Elmin = %g  Elmax = %g  Ekinmax = %g\n",
	  par->Elmin,par->Elmax,par->Ekinmax);
  printf("Temperature = %.3f kbT = %.8f boltzEnergyRange = %.8f\n", 
	  par->temp, par->kbT, par->boltzEnergyRange);

  ist->homoIndex = ist->lumoIndex = ist->totalHomo = ist->totalLumo = 0;
  pf = fopen("eval.par" , "r");
  for (i = ieof = 0; ieof != EOF; i++) {
    ieof = fscanf(pf,"%ld %lg %lg",&a,&evalloc,&deloc);
    if (deloc < par->deps && evalloc < -0.2) ist->homoIndex = i;
    if (deloc < par->deps && evalloc < -0.2) ist->totalHomo++; 
    if (deloc < par->deps && evalloc > -0.2 && ieof != EOF) ist->totalLumo++; 
  }
  fclose(pf);
  nval = i - 1; // the total number of states in original eval.par/psi.par
  pf = fopen("eval.par" , "r");
  for (i = 0; i <= ist->homoIndex; i++) fscanf (pf,"%ld %lg %lg",&a,&evalloc,&deloc);
  for (i = ist->homoIndex+1; i < nval; i++) {
    fscanf(pf,"%ld %lg %lg",&a,&evalloc,&deloc);
    if (deloc < par->deps) {
      ist->lumoIndex = i;
      break;
    }
  }
  fclose(pf);

  printf ("Total number of hole eigenstates = %ld\n",ist->totalHomo);
  printf ("Total number of electron eigenstates = %ld\n",ist->totalLumo);
  printf ("The index of the homo state = %ld\n",ist->homoIndex);
  printf ("The index of the lumo state = %ld\n",ist->lumoIndex);

  eval = (double*)calloc(nval,sizeof(double)); 
  de = (double*)calloc(nval,sizeof(double)); 
  pf = fopen("eval.par" , "r");
  for (i = 0; i < nval; i++){
    fscanf(pf,"%ld %lg %lg", &a, &eval[i], &de[i]);
  }
  fclose(pf);

  par->homoEnergy = eval[ist->homoIndex];
  par->lumoEnergy = eval[ist->lumoIndex];
  par->Egap = par->lumoEnergy - par->homoEnergy;
  par->minInitE = 2.0*par->Egap;
  par->maxInitE = 2.0*par->Egap + par->boltzEnergyRange;
 
  printf("Homo energy = %.10f\n", par->homoEnergy);
  printf("Lumo energy = %.10f\n", par->lumoEnergy);
  printf("Energy gap = %.10f %.8f\n", par->Egap, par->Egap*AUTOEV);
  printf("Min energy of the initial ground state biexciton =  %.10f\n", par->minInitE);
  printf("Max allowed energy of the initial biexciton =  %.10f\n", par->maxInitE);

  // Count the number of hole and electron band-edge states within 
  // boltzEnergyRange (= 3kbT) of the homo and lumo states, respectively
  // count hole states
  ist->numBandEdgeHoles = 0;
  for (i = ist->homoIndex; i >= 0; i--)
    if ((eval[i] >= par->homoEnergy-par->boltzEnergyRange) && de[i] < par->deps) ist->numBandEdgeHoles++;
  // count electron states
  ist->numBandEdgeElectrons = 0;
  for (i = ist->lumoIndex; i < nval; i++)
    if ((eval[i] <= par->lumoEnergy+par->boltzEnergyRange) && de[i] < par->deps) ist->numBandEdgeElectrons++;
  ist->numBandEdgeStates = ist->numBandEdgeHoles+ist->numBandEdgeElectrons;  

  printf("Number of holes within boltzEnergyRange of the homo energy = %ld\n", ist->numBandEdgeHoles);
  printf("Number of electrons within boltzEnergyRange of the lumo energy = %ld\n", ist->numBandEdgeElectrons);
  printf("Total number of band-edge states within boltzEnergyRange = %ld\n", ist->numBandEdgeStates);

  free(eval);  free(de);

  // used in BSE based AR calculation
  ist->msbs = ist->totalHomo + ist->totalLumo;
  ist->msbs2 = ist->totalHomo * ist->totalLumo;

  printf("Total number of non-interacting states used in the BSE %d %d\n", ist->msbs, ist->msbs2);

  return;
}

/*****************************************************************************/

void init(double *vx,double *vy,double *vz,double *ksqr,double *potl,double *rx,double *ry,double *rz,par_st *par,lng_st *ist)
{
  FILE *pf; atm_st *atm;
  long jexc, jx, jy, jz, jyz, jxyz, iatom, ntot, jp, *npot, flags=0;
  double del, mx, my, mz, xd, yd, zd, dx, dy, dz, *ksqrx, *ksqry, *ksqrz;
  double *vr, *potatom, *dr, sum, rex, potEx; 
  pot_st ppar;

  if ((ksqrx  = (double*)calloc(ist->nx,sizeof(double)))==NULL)nerror("ksqrx");
  if ((ksqry  = (double*)calloc(ist->ny,sizeof(double)))==NULL)nerror("ksqry");
  if ((ksqrz  = (double*)calloc(ist->nz,sizeof(double)))==NULL)nerror("ksqrz");
  if ((atm = (atm_st*)calloc(ist->natom,sizeof(atm_st)))==NULL)nerror("atm");
  
  // Read in the external potential and store it in ppar if it exists
  if ( access( "expot.par", F_OK) != -1 ) {
    pf = fopen("expot.par", "r");
    fscanf (pf,"%lg",&ppar.e);  /*** external potential parameter: Energy Scale ***/
    fscanf (pf,"%lg",&ppar.x0); /*** external potential parameter: Ef ***/
    fscanf (pf,"%lg",&ppar.t);  /*** external potential parameter: kT ***/
    fclose(pf);
    printf("expot.par was detected\n");
    printf("External potential energy scale, e = % .6f\n", ppar.e);
    printf("External potential Fermi radius, x0 = % .6f\n", ppar.x0);
    printf("External potential steepness (Temp), kT = % .6f\n", ppar.t);
  } else {
    ppar.e = ppar.x0 = 0.0; // setting ppar.e = 0 results in ext pot being 0
    ppar.t = 1.0;
    printf("No expot.par file detected in cwd - setting external potential to 0!\n");
  }

  /*** read the pasivated nanocrystal configuration ***/
  pf = fopen("conf.par" , "r");
  fscanf(pf,"%ld",&ntot);
  assert(fabs((double)(ntot - ist->natom)) < 1.0e-15);
  ist->mfermi = read_conf(rx, ry, rz, atm, ntot, pf);
  printf ("mfermi = %ld\n",ist->mfermi);
  fclose(pf);

  xd = rint(0.5 * get_dot_ligand_size_z(rx,ntot) + 5.0);
  yd = rint(0.5 * get_dot_ligand_size_z(ry,ntot) + 5.0);
  zd = rint(0.5 * get_dot_ligand_size_z(rz,ntot) + 5.0);
  printf("xd = %g yd = %g zd = %g\n", xd, yd, zd);

  /*xd = rint(0.5 * get_dot_ligand_size(rx,ry,rz,ntot) + 1.0);
  printf ("radius = %g\n",xd);
  yd = zd = xd;*/
  
  /***initial parameters for the pot reduce mass, etc. in the x direction ***/
  par->xmin = -xd;
  par->xmax = xd;
  mx = 1.0;
  par->dx  = (par->xmax - par->xmin) / (double)(ist->nx);
  par->dkx = TWOPI / ((double)ist->nx * par->dx);
  
  /***initial parameters for the pot reduce mass, etc. in the y direction ***/
  par->ymin = -yd;
  par->ymax = yd;
  my = 1.0;
  par->dy  = (par->ymax - par->ymin) / (double)(ist->ny);
  par->dky = TWOPI / ((double)ist->ny * par->dy);

  /***initial parameters for the pot reduce mass, etc. in the z direction ***/
  par->zmin = -zd;
  par->zmax = zd;
  mz = 1.0;
  par->dz  = (par->zmax - par->zmin) / (double)(ist->nz);
  par->dkz = TWOPI / ((double)ist->nz * par->dz);

  par->gamma = 7.0 / (2.0 * xd);
  par->gamma2 = sqr(par->gamma);
  
  par->dv = par->dx * par->dy * par->dz;
  par->dr = sqrt(sqr(par->dx) + sqr(par->dy) + sqr(par->dz));
  printf ("dx = %g dy = %g dz = %g dV = %g dr = %g\n",
	  par->dx, par->dy, par->dz, par->dv, par->dr);

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
  for (jz = 0; jz < ist->nz; jz++) for (jy = 0; jy < ist->ny; jy++){
    for (jyz = ist->nx * (ist->ny * jz + jy), jx = 0; jx < ist->nx; jx++){
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
  
  /*** read pseudopotentials ***/
  dr  = (double*)calloc(ist->natomtype,sizeof(double));
  vr  = (double*)calloc(ist->npot*ist->natomtype,sizeof(double));
  potatom = (double*)calloc(ist->npot*ist->natomtype,sizeof(double));
  npot = (long*)calloc(ist->natomtype,sizeof(long));
  read_pot(vr,potatom,npot,dr,atm,ist->npot,ist->natomtype);

  par->Vmin = 1.0e10;
  par->Vmax = -1.0e10;

  omp_set_dynamic(0);
  omp_set_num_threads(ist->nthreads);
#pragma omp parallel for private(dx,dy,dz,del,jy,jx,jyz,jxyz,sum,iatom)
  for (jz = 0; jz < ist->nz; jz++) {
    for (jy = 0; jy < ist->ny; jy++) {
      jyz = ist->nx * (ist->ny * jz + jy);
      for (jx = 0; jx < ist->nx; jx++) {
        rex = sqrt(vx[jx]*vx[jx]+vy[jy]*vy[jy]+vz[jz]*vz[jz]);
        potEx = expot(rex, ppar);
      	jxyz = jyz + jx;
      	for (sum = 0.0, iatom = 0; iatom < ntot; iatom++) {
      	  dx = vx[jx] - rx[iatom];
      	  dy = vy[jy] - ry[iatom];
      	  dz = vz[jz] - rz[iatom];
      	  del = sqrt(dx * dx + dy * dy + dz * dz);
      	  sum += interpolate(del,dr[atm[iatom].natyp],vr,potatom,ist->npot,npot[atm[iatom].natyp],atm[iatom].natyp);
      	}
      	potl[jxyz] = sum+potEx;
      }
    }
  }

  // Print out the external potential as a function of the absolute value of the z coordinate
  pf = fopen("pS.dat", "w");
  for (jz = 0; jz < ist->nz; jz++) {
    potEx = expot(fabs(vz[jz]), ppar);
    fprintf(pf, "% .6f % .8f\n", vz[jz], potEx);
  }
  fclose(pf);  

  for (jxyz= 0; jxyz < ist->ngrid; jxyz++){
    if (par->Vmax < potl[jxyz]) par->Vmax = potl[jxyz];
    if (par->Vmin > potl[jxyz]) par->Vmin = potl[jxyz];
  }

  printf ("dV = %g vmin = %g Vmax = %g\n",
	  par->Vmax-par->Vmin,par->Vmin,par->Vmax);

  par->dE = 0.5 * sqr(PIE) / (mx*par->dx*par->dx) +
    0.5 * sqr(PIE) / (my*par->dy*par->dy) +    
    0.5 * sqr(PIE) / (mz*par->dz*par->dz);
  printf ("dT = %g\n",par->dE);

  free(dr); free(vr); free (atm);  free(npot);

  return;
}

/****************************************************************************/

void init_pot(double *vx,double *vy,double *vz,zomplex *potq,par_st par,lng_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  long jx, jy, jz, jyz, jxyz, sx, sy, sz;
  double dr, x2, y2, z2, *kx2, *ky2, *kz2, alpha, cosa, sina;
  double boxl = (double)(ist.nx) * par.dx;
  zomplex *potr, tmp;

  if ((kx2  = (double*)calloc(ist.nx,sizeof(double)))==NULL)nerror("kx2");
  if ((ky2  = (double*)calloc(ist.ny,sizeof(double)))==NULL)nerror("ky2");
  if ((kz2  = (double*)calloc(ist.nz,sizeof(double)))==NULL)nerror("kz2");
  if ((potr = (zomplex*)calloc(ist.ngrid,sizeof(zomplex)))==NULL)nerror("potr");
  
  for (jx = 0; jx < ist.ngrid; jx++) potr[jx].re = potr[jx].im = 0.0;
  for (jz = 0; jz < ist.nz; jz++){
    z2 = sqr(vz[jz]);
    for (jy = 0; jy < ist.ny; jy++){
      y2 = sqr(vy[jy]);
      jyz = ist.nx * (ist.ny * jz + jy);
      for (jx = 0; jx < ist.nx; jx++){
	x2 = sqr(vx[jx]);
	jxyz = jyz + jx;
	dr = sqrt(x2 + y2 + z2);
	if (dr < boxl) potr[jxyz].re = screenedcoulomb(dr,par.gamma);
      }
    }
  }

  
  memcpy(&fftwpsi[0],&potr[0],ist.ngrid*sizeof(fftwpsi[0]));
  fftw_execute(planfw);
  memcpy(&potq[0],&fftwpsi[0],ist.ngrid*sizeof(potq[0]));
  /*fftwnd_one(planfw,potr,potq);*/

  /*for (jx = 0; jx < ist.nx; jx++)
    kx2[jx] = sqr(par.kxmin + (double)(jx) * par.dkx);
  for (jy = 0; jy < ist.ny; jy++)
    ky2[jy] = sqr(par.kymin + (double)(jy) * par.dky);
  for (jz = 0; jz < ist.nz; jz++)
  kz2[jz] = sqr(par.kzmin + (double)(jz) * par.dkz);*/
  
  for (kx2[0] = 0.0, jx = 1; jx <= ist.nx / 2; jx++)
    kx2[jx] = (kx2[ist.nx-jx] = sqr((double)(jx) * par.dkx));
  for (ky2[0] = 0.0, jy = 1; jy <= ist.ny / 2; jy++)
    ky2[jy] = (ky2[ist.ny-jy] = sqr((double)(jy) * par.dky));
  for (kz2[0] = 0.0, jz = 1; jz <= ist.nz / 2; jz++)
    kz2[jz] = (kz2[ist.nz-jz] = sqr((double)(jz) * par.dkz));

  for (sz = 1.0, jz = 0; jz < ist.nz; jz++, sz = -sz){
    z2 = kz2[jz];
    for (sy = 1.0, jy = 0; jy < ist.ny; jy++, sy = -sy){
      y2 = ky2[jy];
      jyz = ist.nx * (ist.ny * jz + jy);
      for (sx = 1.0, jx = 0; jx < ist.nx; jx++, sx = -sx){
	x2 = kx2[jx];
	jxyz = jyz + jx;
	alpha = PIE * (double)(jx + jy + jz + ist.ngrid / 2);
	cosa = cos(alpha);
	sina = sin(alpha);
	
	/*potq[jxyz].re *= (double)(sx * sy * sz) * par.dv;
	  potq[jxyz].im *= (double)(sx * sy * sz) * par.dv;*/

	tmp.re = potq[jxyz].re;
	tmp.im = potq[jxyz].im;

	potq[jxyz].re = (tmp.re * cosa - tmp.im * sina) * par.dv;
	potq[jxyz].im = (tmp.re * sina + tmp.im * cosa) * par.dv;
	
	double  q2 = x2 + y2 + z2 + 1.0e-20;
	/*potq[jxyz].re += (FOURPI / (q2 + par.gamma2));*/
	potq[jxyz].re += (FOURPI * (1.0 - exp(-0.25* q2 / par.gamma2)) / q2);


	potq[jxyz].re *= ist.ngrid_1;
	potq[jxyz].im *= ist.ngrid_1;

	/*printf ("%g %g %g %g\n",sqrt(x2+y2+z2),potq[jxyz].re,potq[jxyz].im,
	  (FOURPI / (x2 + y2 + z2 + par.gamma2)));*/
      }
    }
  }
  free(potr);  free(kx2); free(ky2); free(kz2);
  return;
}

/************************************************************************/

#define SQRTPI       (sqrt(3.14159265358979323846))

double screenedcoulomb(double dr,double gamma)
{
  /*if (dr < EPSR) return (gamma);*/
  if (dr < EPSR) return (2.0*gamma/SQRTPI);
  return (erf(gamma * dr) / dr);
  /*return ((1.0 - exp(-gamma * dr)) / dr);*/
}

/************************************************************************/

void init_psi(zomplex *psi,lng_st ist,par_st par,long *idum)
{
  long jx, jy, jz, jzy, jxyz;
  long tidum = (*idum);

  for (jz = 0; jz < ist.nz; jz++) for (jy = 0; jy < ist.ny; jy++){
    for (jzy = ist.nx * (ist.ny * jz + jy), jx = 0; jx < ist.nx; jx++){
      jxyz = jzy + jx;
      psi[jxyz].re = (-1.0 + 2.0 * ran_nrc(&tidum));
      psi[jxyz].im = 0.0;
    }
  }
  normalize(psi,par.dv,ist.ngrid);
  (*idum) = tidum;
  return;
}

/************************************************************************/

double expot(double r, pot_st ppar) {
  return (ppar.e*(1.0/(exp((r-ppar.x0)/ppar.t)+1.0)-1.0));
}

/************************************************************************/
