#include "fd.h"

/*************************************************************************/

void init_size(long argc, char *argv[],atm_st **atm,par_st *par,long_st *ist, double **rx, double **ry, double **rz)
{
  long ntot, ie;
  double mx, my, mz, xd, yd, zd;
  FILE* pf = fopen("input.par" , "r");
  fscanf (pf,"%lg",&par->dx);  /*** dx ***/
  fscanf (pf,"%lg",&par->dy);  /*** dy ***/
  fscanf (pf,"%lg",&par->dz);  /*** dz ***/
  fscanf (pf,"%ld",&ist->ms);  /*** number of state per filter ***/
  fscanf (pf,"%ld",&ist->ns);  /*** number of filter cycles ***/
  fscanf (pf,"%ld",&ist->nc);    /*** length of newton interpolation ***/
  fscanf (pf,"%ld",&ist->npot);    /*** size of pseudopotential files ***/
  fscanf (pf,"%lg %lg",&par->Elmin, &par->holeMax); /*** hole energy window ***/
  fscanf (pf,"%lg %lg",&par->electronMin,&par->Elmax); /*** electron energy window ***/
  fscanf (pf,"%lg",&par->holeFilts); /*** Hole filter fraction ***/
  fscanf (pf,"%lg",&par->Ekinmax); /*** kinetic energy cutoff ***/
  fscanf (pf,"%ld",&ist->flaghomo); /*** flaghomo ***/
  fscanf (pf,"%ld %ld %lg",&par->hprint,&par->eprint,&par->sigma_e); /*** holes to print, electrons to print, maximum eps ***/
  fscanf (pf,"%ld",&ist->nthreads);
  fclose(pf);

  printf ("Done reading input\n");

  pf = fopen("conf.par" , "r");
  fscanf(pf,"%ld",&ist->natom);
  //printf ("ist->natom=%ld, sizeof(double*)=%ld\n",ist->natom,sizeof(double));
  fclose(pf);
  if ((*rx = (double*)calloc(ist->natom,sizeof(double)))==NULL)nerror("rx");
  if ((*ry = (double*)calloc(ist->natom,sizeof(double)))==NULL)nerror("ry");
  if ((*rz = (double*)calloc(ist->natom,sizeof(double)))==NULL)nerror("rz");
  if ((*atm = (atm_st*)calloc(ist->natom,sizeof(atm_st)))==NULL)nerror("atm");
  //for (ie = 0; ie < ist->natom; ie++){printf("rx[%ld]=%g\n",ie,*rx[ie]);}
  pf = fopen("conf.par" , "r");
  ist->natomtype = 16;
  fscanf(pf,"%ld",&ntot);
  assert(fabs((double)(ntot - ist->natom)) < 1.0e-15);
  printf ("ntot=%ld\n",ntot);
  read_conf(*rx,*ry,*rz,*atm,ntot,pf);
  fclose (pf);
  //for (ie = 0; ie < ist->natom; ie++){printf("rx[%ld]=%g\n",ie,*rx[ie]);}
  printf ("Done reading conf\n");
  //printf ("dot radius = %g\n",rd);
  //rd = rint(0.5 * get_dot_ligand_size(rx,ry,rz,ntot) + 3.0);
   
  xd = rint(0.5 * get_dot_ligand_size_z(*rx,ntot) + 5.0);
  yd = rint(0.5 * get_dot_ligand_size_z(*ry,ntot) + 5.0);
  zd = rint(0.5 * get_dot_ligand_size_z(*rz,ntot) + 5.0);
  printf ("xd = %g yd = %g zd = %g\n",xd,yd,zd);
  
  /***initial parameters for the pot reduce mass, etc. in the x direction ***/
  par->xmin = -xd;
  par->xmax = xd;
  mx = 1.0;
  ist->nx = rndToThreads((int)((par->xmax - par->xmin) / (double)(par->dx)), ist->nthreads); //Calculate and round up
  //ist->nx = (int)((par->xmax - par->xmin) / (par->dx));
  //printf ("nx = %ld\n",ist->nx);
  par->dkx = TWOPI / ((double)ist->nx * par->dx);
   
  /***initial parameters for the pot reduce mass, etc. in the y direction ***/
  par->ymin = -yd;
  par->ymax = yd;
  my = 1.0;
  //ist->ny = (int)((par->ymax - par->ymin) / (double)(par->dy));
  ist->ny = rndToThreads((int)((par->ymax - par->ymin) / (double)(par->dy)), ist->nthreads);
  par->dky = TWOPI / ((double)ist->ny * par->dy);

  /***initial parameters for the pot reduce mass, etc. in the z direction ***/
  par->zmin = -zd;
  par->zmax = zd;
  mz = 1.0;
  ist->nz = rndToThreads((int)((par->zmax - par->zmin) / (double)(par->dz)), ist->nthreads);
  //ist->nz = (int)((par->zmax - par->zmin) / (par->dz));
  par->dkz = TWOPI / ((double)ist->nz * par->dz);
 
  ist->nx_1 = 1.0 / (double)(ist->nx);
  ist->ny_1 = 1.0 / (double)(ist->ny);
  ist->nz_1 = 1.0 / (double)(ist->nz);
  ist->ngrid = ist->nx * ist->ny * ist->nz;
     
  ist->mstot = ist->ms * ist->ns;
  
  printf ("nx = %ld  ny = %ld  nz = %ld npot = %ld\n",
	  ist->nx,ist->ny,ist->nz,ist->npot);
  printf ("ms = %ld  ns = %ld natom = %ld nc = %ld\n",
	  ist->ms,ist->ns,ist->natom,ist->nc);
  printf ("Elmin = %g  Elmax = %g  Ekinmax = %g\n",
	  par->Elmin,par->Elmax,par->Ekinmax);
  printf ("threads = %ld\n",ist->nthreads);
  printf ("dx = %g dy = %g dz = %g\n",
           par->dx,par->dy,par->dz);
  
return;
}

/*************************************************************************/

void init(double *potl,double *vx,double *vy,double *vz,double *ksqr,double *rx,double *ry,double *rz,atm_st *atm,par_st *par,double *eval,long_st *ist,fftw_plan_loc *planfw,fftw_plan_loc *planbw,fftw_complex *fftwpsi)
{
  FILE *pf;
  long jx, jy, jz, jyz, jxyz, ie, jp, *npot, nn, flags=0;
  double del, del_h, del_e, mx, my, mz, xd, yd, zd, dx, dy, dz, *ksqrx, *ksqry, *ksqrz;
  double *vr, *potatom, *dr, sum;

  if ((ksqrx  = (double*)calloc(ist->nx,sizeof(double)))==NULL)nerror("ksqrx");
  if ((ksqry  = (double*)calloc(ist->ny,sizeof(double)))==NULL)nerror("ksqry");
  if ((ksqrz  = (double*)calloc(ist->nz,sizeof(double)))==NULL)nerror("ksqrz");
  //for (ie = 0; ie < ist->natom; ie++){printf("rx[%ld]=%lg\n",ie,rx[ie]);} 
  mx = 1.0;
  my = 1.0;
  mz = 1.0;
  par->Vmin = 1.0e10;
  par->Vmax = -1.0e10;

  par->dv = par->dx * par->dy * par->dz;
  par->dr = sqrt(sqr(par->dx) + sqr(par->dy) + sqr(par->dz));
  printf ("In init():\n\ndx = %g dy = %g dz = %g dv = %g dr = %g\n",
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
  for (jz = 0; jz < ist->nz; jz++) for (jy = 0; jy < ist->ny; jy++){
    for (jyz = ist->nx * (ist->ny * jz + jy), jx = 0; jx < ist->nx; jx++){
      jxyz = jyz + jx;
      ksqr[jxyz] = ksqrx[jx] + ksqry[jy] + ksqrz[jz];
      if (ksqr[jxyz] > par->Ekinmax) ksqr[jxyz] = par->Ekinmax;
    }
  }
  free(ksqrx); free(ksqry);  free(ksqrz);
  
  /*** read pseudopotentials ***/
  dr  = (double*)calloc(ist->natomtype,sizeof(double));
  vr  = (double*)calloc(ist->npot*ist->natomtype,sizeof(double));
  potatom = (double*)calloc(ist->npot*ist->natomtype,sizeof(double));
  npot = (long*)calloc(ist->natomtype,sizeof(long));
  read_pot(vr,potatom,npot,dr,atm,ist->npot,ist->natomtype);
 printf("ist->natom=%ld\n",ist->natom); 
  /***initializing the potential vector  ***/
  for (jx = 0, dx = par->xmin; jx < ist->nx; jx++, dx += par->dx) vx[jx] = dx;
  for (jy = 0, dy = par->ymin; jy < ist->ny; jy++, dy += par->dy) vy[jy] = dy;
  for (jz = 0, dz = par->zmin; jz < ist->nz; jz++, dz += par->dz) vz[jz] = dz;
  omp_set_dynamic(0);
  omp_set_num_threads(ist->nthreads);
#pragma omp parallel for private(dx,dy,dz,del,jy,jx,jyz,jxyz,sum,ie)
  for (jz = 0; jz < ist->nz; jz++){
  //printf("jz=%d\n",jz);
    for (jy = 0; jy < ist->ny; jy++){
      jyz = ist->nx * (ist->ny * jz + jy);
      for (jx = 0; jx < ist->nx; jx++){
	jxyz = jyz + jx;
	for (sum = 0.0, ie = 0; ie < ist->natom; ie++){
	  //printf("dz=vz[%ld]-rz[%ld]=%lg\n",jz,ie,dz); 
          dx = vx[jx] - rx[ie];
	  dy = vy[jy] - ry[ie];
	  dz = vz[jz] - rz[ie];
	  del = sqrt(dx * dx + dy * dy + dz * dz);
	  //printf("del=%lg\n",del);
	  sum += interpolate(del,dr[atm[ie].natyp],vr,potatom,ist->npot,npot[atm[ie].natyp],atm[ie].natyp);
	  //printf("Sum=%lg\n",sum);
	}
	potl[jxyz] = sum;
      }
    }
  }
  free(vr);  free(potatom);  free(dr); free(npot);
  printf("Finished potl\n");
  for (jxyz= 0; jxyz < ist->ngrid; jxyz++){
    if (par->Vmax < potl[jxyz]) par->Vmax = potl[jxyz];
    if (par->Vmin > potl[jxyz]) par->Vmin = potl[jxyz];
  }
  /*write_pot(vx,vy,vz,potl);*/
  /*par->Vmin = -5.0;*/
  printf ("dV = %g vmin = %g Vmax = %g\n",
	  par->Vmax-par->Vmin,par->Vmin,par->Vmax);

  par->dE = 0.5 * sqr(PIE) / (mx*par->dx*par->dx) +
    0.5 * sqr(PIE) / (my*par->dy*par->dy) +    
    0.5 * sqr(PIE) / (mz*par->dz*par->dz);
  printf ("dT = %g\n",par->dE);
  
  /*** setting the energy grid El ***/
/*
  del = (par->Elmax - par->Elmin) / (double)(ist->ms-1);
  for (jx = 0; jx < ist->ms; jx++) eval[jx] = par->Elmin + (double)(jx) * del;
  for (pf = fopen("Egrid.dat","w"),jx = 0; jx < ist->ms; jx++)
    fprintf (pf,"%g\n",eval[jx]);
*/
  del_h = (par->holeMax - par->Elmin) / (double)(par->holeFilts*(ist->ms)-1);
  del_e = (par->Elmax - par->electronMin) / (double)((1-par->holeFilts)*(ist->ms)-1);
  for (jx = 0; jx < par->holeFilts*ist->ms; jx++) eval[jx] = par->Elmin + (double)(jx) * del_h;
  for (jx = 0; jx < (1-par->holeFilts)*ist->ms; jx++) eval[(int)(par->holeFilts*(ist->ms)+jx)] = par->electronMin + (double)(jx) * del_e; 
  for (pf = fopen("Egrid.dat","w"),jx = 0; jx < ist->ms; jx++)
    fprintf (pf,"%g\n",eval[jx]);
  fclose(pf);
  
  /*** initialization for the fast Fourier transform ***/
  //(*planfw) = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_FORWARD,flags);
  //(*planbw) = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_BACKWARD,flags);
  //(*planfw) = fftw3d_create_plan(ist->nz,ist->ny,ist->nx,FFTW_FORWARD,flags);
  //(*planbw) = fftw3d_create_plan(ist->nz,ist->ny,ist->nx,FFTW_BACKWARD,flags);
  
  /*ist->lumo = 4*n2;
    ist->homo = 4*n2-1;*/

  for (nn = ie = 0; ie < ist->natom; ie++)
    if (((atm[ie].natyp < 8) || (atm[ie].natyp > 11)) && (atm[ie].natyp % 2)) nn++;
  printf ("nn = %ld\n",nn);
  ist->homo = 4*nn-1;
  ist->lumo = ist->homo+1;
  printf ("homo = %ld lumo = %ld\n",ist->homo,ist->lumo);
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

int rndToThreads(int gridpts, int threads)
{
    int remainder = gridpts % threads;
    if (remainder == 0)
        return gridpts;

    return (gridpts + threads - remainder);
}
