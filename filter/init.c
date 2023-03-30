/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "fd.h"
#include <time.h>

/****************************************************************************/

void init_size(long argc, char *argv[],par_st *par,long_st *ist)
{
  FILE* pf = fopen("input.par" , "r");
  fscanf (pf,"%ld",&ist->nx);  /*** number of grid point in x ***/
  fscanf (pf,"%ld",&ist->ny);  /*** number of grid point in y ***/
  fscanf (pf,"%ld",&ist->nz);  /*** number of grid point in z ***/
  fscanf (pf,"%ld",&ist->ms);  /*** number of state per filter ***/
  fscanf (pf,"%ld",&ist->ns);  /*** number of filter cycles ***/
  fscanf (pf,"%ld",&ist->nc);    /*** length of newton interpolation ***/
  fscanf (pf,"%lg",&par->Elmin); /***minimum energy ***/
  fscanf (pf,"%lg",&par->Elmax); /*** maximum energy ***/
  fscanf (pf,"%ld",&ist->nthreads);
  fclose(pf);

  pf = fopen("conf.par" , "r");
  fscanf(pf,"%ld",&ist->natom);
  fclose(pf);

  ist->npot = 8192;
  par->Ekinmax = 10.0;  
  ist->nx_1 = 1.0 / (double)(ist->nx);
  ist->ny_1 = 1.0 / (double)(ist->ny);
  ist->nz_1 = 1.0 / (double)(ist->nz);
  ist->ngrid = ist->nx * ist->ny * ist->nz;
  ist->natomtype = 17;

  ist->mstot = ist->ms * ist->ns; // number of states per filter * number of filters

  printf("Number of atoms in the system, natom = %ld\n", ist->natom);
  printf("Number of grid points used: nx = %ld  ny = %ld  nz = %ld\n", ist->nx, ist->ny, ist->nz);
  printf("Number of states per filter, ms = %ld\nNumber of filter cycles, ns = %ld\n", ist->ms, ist->ns);
  printf("Length of Newton interpolation used, nc = %ld\n", ist->nc);
  printf("Target energy for HOMO state = % .4f\n", par->Elmin);
  printf("Target energy for LUMO state = % .4f\n", par->Elmax);
  printf("Kinetic energy maximum = % .2f\n", par->Ekinmax); 
  printf("Number of openMP threads used = %ld\n", ist->nthreads);

  return;
}

/*************************************************************************/

void init(double *potl,double *vx,double *vy,double *vz,double *ksqr,double *rx,double *ry,double *rz,atm_st *atm,par_st *par,double *eval,long_st *ist,fftw_plan_loc *planfw,fftw_plan_loc *planbw,fftw_complex *fftwpsi)
{
  FILE *pf, *pfs; pot_st ppar;
  long jx, jy, jz, jyz, jxyz, ie, ntot, jp, *npot, nn, flags=0, ii;
  double del, mx, my, mz, xd, yd, zd, dx, dy, dz, *ksqrx, *ksqry, *ksqrz;
  double *vr, *potatom, *dr, sum, rex, potEx;

  if ((ksqrx  = (double*)calloc(ist->nx,sizeof(double)))==NULL)nerror("ksqrx");
  if ((ksqry  = (double*)calloc(ist->ny,sizeof(double)))==NULL)nerror("ksqry");
  if ((ksqrz  = (double*)calloc(ist->nz,sizeof(double)))==NULL)nerror("ksqrz");

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
  } 
  else {
    ppar.e = ppar.x0 = 0.0; // setting ppar.e = 0 results in ext pot being 0
    ppar.t = 1.0;
    printf("No expot.par file detected in cwd - setting external potential to 0!\n");
  }
  
  /*** read the pasivated nanocrystal configuration ***/
  pf = fopen("conf.par" , "r");
  fscanf(pf,"%ld",&ntot);
  assert(fabs((double)(ntot - ist->natom)) < 1.0e-15);
  read_conf(rx,ry,rz,atm,ntot,pf);
  fclose(pf);

  xd = rint(0.5 * get_dot_ligand_size_z(rx,ntot) + 5.0);
  yd = rint(0.5 * get_dot_ligand_size_z(ry,ntot) + 5.0);
  zd = rint(0.5 * get_dot_ligand_size_z(rz,ntot) + 5.0);
  printf("Box (quadrant) dimensions: xd = %.2f yd = %.2f zd = %.2f\n", xd, yd, zd);
  
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

  par->Vmin = 1.0e10;
  par->Vmax = -1.0e10;

  par->dv = par->dx * par->dy * par->dz;
  par->dr = sqrt(sqr(par->dx) + sqr(par->dy) + sqr(par->dz));
  printf("Grid point spacing: dx = %.4f dy = %.4f dz = %.4f dv = %.6f dr = %.6f\n",
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
  
  /*** read pseudopotentials ***/
  fprintf(stdout, "Number of atom types = %ld\n", ist->natomtype); fflush(stdout);

  dr  = (double *) calloc(ist->natomtype, sizeof(double));
  vr  = (double *) calloc(ist->npot*ist->natomtype, sizeof(double));
  potatom = (double *) calloc(ist->npot*ist->natomtype, sizeof(double));
  npot = (long *) calloc(ist->natomtype, sizeof(long));
  read_pot(vr, potatom, npot, dr, atm, ist->npot, ist->natomtype);
  
  /***initializing the potential vector  ***/
  for (jx = 0, dx = par->xmin; jx < ist->nx; jx++, dx += par->dx) vx[jx] = dx;
  for (jy = 0, dy = par->ymin; jy < ist->ny; jy++, dy += par->dy) vy[jy] = dy;
  for (jz = 0, dz = par->zmin; jz < ist->nz; jz++, dz += par->dz) vz[jz] = dz;

  omp_set_dynamic(0);
  omp_set_num_threads(ist->nthreads);
#pragma omp parallel for private(dx,dy,dz,del,jy,jx,jyz,jxyz,sum,ie)
  for (jz = 0; jz < ist->nz; jz++) {
    for (jy = 0; jy < ist->ny; jy++) {
      jyz = ist->nx * (ist->ny * jz + jy);
      for (jx = 0; jx < ist->nx; jx++) {
        rex = sqrt(vx[jx]*vx[jx]+vy[jy]*vy[jy]+vz[jz]*vz[jz]);
		// TODO: have a more flexible potEx / maybe have precalculated potEx at each grid point 
        potEx = expot(rex,ppar);
        jxyz = jyz + jx;
  	    for (sum = 0.0, ie = 0; ie < ntot; ie++) {
  	      dx = vx[jx] - rx[ie];
  	      dy = vy[jy] - ry[ie];
  	      dz = vz[jz] - rz[ie];
  	      del = sqrt(dx * dx + dy * dy + dz * dz);
  	      sum += interpolate(del,dr[atm[ie].natyp],vr,potatom,ist->npot,npot[atm[ie].natyp],atm[ie].natyp);
  	    }
        potl[jxyz] = sum+potEx;
      }
    }
  }

  // Print out the external potential as a function of the distance from the origin 
  pfs = fopen("pS.dat", "w");
  for (jz = 0; jz < ist->nz; jz++) {
    potEx=expot(fabs(vz[jz]),ppar);
    fprintf(pfs,"% .6f % .8f\n", vz[jz], potEx);
  }
  fclose(pfs);

  writeCubeFile(potl, par->xmin, par->ymin, par->zmin, par->dx, par->dy, par->dz, ist->nx, ist->ny, ist->nz);

  free(vr); free(potatom); free(dr); free(npot);  

  for (jxyz= 0; jxyz < ist->ngrid; jxyz++){
    if (par->Vmax < potl[jxyz]) par->Vmax = potl[jxyz];
    if (par->Vmin > potl[jxyz]) par->Vmin = potl[jxyz];
  }
  /*write_pot(vx,vy,vz,potl);*/
  printf("Potential energy statistics: dV = %.5f Vmin = % .5f Vmax = %.5f\n",
	  par->Vmax-par->Vmin,par->Vmin,par->Vmax);

  par->dE = 0.5 * sqr(PIE) / (mx*par->dx*par->dx) +
    0.5 * sqr(PIE) / (my*par->dy*par->dy) +    
    0.5 * sqr(PIE) / (mz*par->dz*par->dz);
  printf("dT = %g\n", par->dE);
  
  /*** setting the energy grid El ***/
  del = (par->Elmax - par->Elmin) / (double)(ist->ms-1);
  for (jx = 0; jx < ist->ms; jx++) eval[jx] = par->Elmin + (double)(jx) * del;
  for (pf = fopen("Egrid.dat","w"),jx = 0; jx < ist->ms; jx++)
    fprintf(pf,"%g\n", eval[jx]);
  fclose(pf);
  
  /*** initialization for the fast Fourier transform ***/
  //(*planfw) = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_FORWARD,flags);
  //(*planbw) = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_BACKWARD,flags);
  //(*planfw) = fftw3d_create_plan(ist->nz,ist->ny,ist->nx,FFTW_FORWARD,flags);
  //(*planbw) = fftw3d_create_plan(ist->nz,ist->ny,ist->nx,FFTW_BACKWARD,flags);
  
  for (nn = ie = 0; ie < ntot; ie++) {
    if (((atm[ie].natyp < 8) || (atm[ie].natyp > 11)) && (atm[ie].natyp % 2)) {
      nn++;
    }
  }
  printf("nn = %ld\n",nn);
  ist->homo = 4*nn-1;
  ist->lumo = ist->homo+1;
  printf("homo = %ld lumo = %ld\n",ist->homo,ist->lumo);

  return;
}
/****************************************************************************/

void init_psi(zomplex *psi,long_st ist,par_st par,long *idum)
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
  normalize(psi,par.dv,ist.ngrid,ist.nthreads);
  (*idum) = tidum;
  return;
}

/****************************************************************************/

double expot(double r, pot_st ppar) {
  return (ppar.e*(1.0/(exp((r-ppar.x0)/ppar.t)+1.0)-1.0));
}

/****************************************************************************/
