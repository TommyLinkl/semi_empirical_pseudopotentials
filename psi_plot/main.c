/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "fd.h"

/*****************************************************************************/

int main(int argc, char *argv[]) {
  FILE *ppsi;  zomplex *psi, *phi, *an; long idum;
  fftw_plan_loc planfw, planbw; fftw_complex *fftwpsi;
  par_st par;   long_st  ist; xyz_st dipole; atm_st *atm;
  double *eval, *ksqr, *vx, *vy, *vz, *zn, *potl;
  double *psitot, *el, *sige, *rx, *ry, *rz, *rho, tci, twi;
  long jgrid, jms, jns, mssav, mstot, flags=0, tid;
  int n,m; // n: starting states; m: include m states
  n = 487;  // 487-InP, 471-GaP
  m = 6;  //6


  /*** read initial setup from input.par ***/
  writeCurrentTime(stdout);
  writeSeparation(stdout);
  init_size(argc,argv,&par,&ist);
  
  /*************************************************************************/
  /*** allocating memory ***/

  fftwpsi = fftw_malloc(sizeof (fftw_complex )*ist.ngrid);
  //if ((fftwpsi  = (fftw_complex*)calloc(ist.ngrid,sizeof(fftw_complex)))==NULL)nerror("fftwpsi");
  if ((psi   = (zomplex*)calloc(ist.ngrid,sizeof(zomplex)))==NULL)nerror("psi");
  if ((phi   = (zomplex*)calloc(ist.ngrid,sizeof(zomplex)))==NULL)nerror("phi");
  /*** the pseudopotential stored on the grid ***/
  if ((potl  = (double*)calloc(ist.ngrid,sizeof(double)))==NULL)nerror("potl");
  /*** the kinetic energy stored on the grid ***/
  if ((ksqr = (double*)calloc(ist.ngrid,sizeof(double)))==NULL)nerror("ksqr");
  /*** the grid in the x, y, and z directions ***/
  if ((vx = (double*)calloc(ist.nx,sizeof(double)))==NULL)nerror("vx");
  if ((vy = (double*)calloc(ist.ny,sizeof(double)))==NULL)nerror("vy");
  if ((vz = (double*)calloc(ist.nz,sizeof(double)))==NULL)nerror("vz");
  /*** the positions of the atoms in the x, y, and z directions ***/
  if ((rx = (double*)calloc(ist.natom,sizeof(double)))==NULL)nerror("rx");
  if ((ry = (double*)calloc(ist.natom,sizeof(double)))==NULL)nerror("ry");
  if ((rz = (double*)calloc(ist.natom,sizeof(double)))==NULL)nerror("rz");
  if ((atm = (atm_st*)calloc(ist.natom,sizeof(atm_st)))==NULL)nerror("atm");
  if ((rho  = (double*)calloc(ist.ngrid,sizeof(double)))==NULL)nerror("rho");
  /*** all wavefunctions in the energy range  ***/
  if ((psitot = (double*)calloc(ist.mstot*ist.ngrid,sizeof(double)))==NULL)nerror("psitot");
  /*** the filtered energies ***/
  if ((eval  = (double*)calloc(ist.mstot,sizeof(double)))==NULL)nerror("eval");
  if ((sige = (double*)calloc(ist.mstot,sizeof(double)))==NULL)nerror("sige");
  /*** the energy of the filters ***/
  if ((el   = (double*)calloc(ist.ms,sizeof(double)))==NULL)nerror("el");
   
  /**************************************************************************/
  /*** initialize the pseudopotential, the kinetic energy, ***/
  /*** the fourier transorm, the grid, and the energy window ***/
  init(potl,vx,vy,vz,ksqr,rx,ry,rz,atm,&par,el,&ist,&planfw,&planbw,fftwpsi);
  /*** initialization for the fast Fourier transform ***/
  planfw = fftw_plan_dft_3d(ist.nz,ist.ny,ist.nx,fftwpsi,fftwpsi,FFTW_FORWARD,flags);
  planbw = fftw_plan_dft_3d(ist.nz,ist.ny,ist.nx,fftwpsi,fftwpsi,FFTW_BACKWARD,flags);
  
  /**************************************************************************/
  /*** calcualte the energy range of the hamitonian ***/
  tci = (double)clock(); twi = (double)time(NULL);
  //get_energy_range(psi,phi,potl,vx,vy,vz,ksqr,&par,ist,planfw,planbw,fftwpsi);
  printf("done calculate energy range, CPU time (sec) %g, wall run time (sec) %g\n",
    ((double)clock()-tci)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-twi); fflush(0);

  /**************************************************************************/
  /*** set parameters for the newton interpolation ***/
  par.dt = sqr((double)(ist.nc) / (2.5 * par.dE));
  printf ("nc = %ld dt = %g dE = %g\n",ist.nc,par.dt,par.dE); fflush(0);
  an = (zomplex*)calloc(ist.nc*ist.ms,sizeof(zomplex));
  zn = (double*)calloc(ist.nc,sizeof(double));
  // coefficient(an,zn,el,par,ist);
  
  slice_psi(ist,par,n,m); // put psitot as the argument and save directly
  ppsi = fopen("psitot.dat","r");
  printf ("reading the first state %d from psitot.dat\n",n);
  for (int i=0; i < m; i++) fread (&psitot[i*ist.ngrid],sizeof(double),ist.ngrid,ppsi); // original correct input;

    /***Convert to .cube files for LMO***/
  char filename2[20];
  for (int j = 0; j<m; j++)
  {
    for (int i = 0; i<ist.ngrid; i++)
    {
      rho[i] = psitot[i+j*ist.ngrid];
    } 
  sprintf(filename2, "psi%i.cube", j+n);

  writeCubeFile(rho, par, ist, filename2);

  //Check for normalization
  double normal=0.0;
  int ix, iy, iz, iyz, ixyz;
  for (iz = 0; iz < ist.nz; iz++){
    for (iy = 0; iy < ist.ny; iy++){
      iyz = ist.nx * (ist.ny * iz + iy);
      for (ix = 0; ix < ist.nx; ix++){
        ixyz = iyz + ix;

        normal += rho[ixyz] * rho[ixyz] * par.dv;
      }
    }
  }
  printf("Normalization of the wavefunction No %d is: %.4f .....\n", j+n, normal);




  char filename3[20];
  double grid_x, grid_y, grid_z, dist_r; 
  long iGrid, iX, iY, iZ, iYZ;
  double box_length = 8;
  /***
  for (int j = 0; j<m; j++) {
    for (iX = 0; iX < ist.nx; iX++) {
      for (iY = 0; iY < ist.ny; iY++) {
        for (iZ = 0; iZ < ist.nz; iZ++) {
          grid_x = par.dx * (iX-0.5*ist.nx);
          grid_y = par.dy * (iY-0.5*ist.ny);
          grid_z = par.dz * (iZ-0.5*ist.nz);
          dist_r = sqrt(sqr(grid_x) + sqr(grid_y) + sqr(grid_z));
          iYZ = ist.nx * (ist.ny * iZ + iY);
          iGrid = iYZ + iX;
          // rho[iGrid] = psitot[iGrid+j*ist.ngrid] / exp(-dist_r);
          if (grid_x <= box_length && grid_y <= box_length && grid_z <= box_length && dist_r<=box_length) {
            rho[iGrid] = psitot[iGrid+j*ist.ngrid];
          }
          else {rho[iGrid] = 0;}
        }
      }
    }
  }
  ***/
  sprintf(filename3, "psi%i_unitCell.cube", j+n);

  //writeCubeFile(rho, par, ist, filename3);
  writeCubeFile_cubicUnitCell(rho, par, ist, filename3, -10.0, -10.0, -10.0, 20.0);   //-7.0, -7.0, -7.0, 14.0

  /***
  char filename4[20];
  for (int j = 0; j<m; j++) {
    for (iX = 0; iX < ist.nx; iX++) {
      for (iY = 0; iY < ist.ny; iY++) {
        for (iZ = 0; iZ < ist.nz; iZ++) {
          grid_x = par.dx * (iX-0.5*ist.nx);
          grid_y = par.dy * (iY-0.5*ist.ny);
          grid_z = par.dz * (iZ-0.5*ist.nz);
          dist_r = sqrt(sqr(grid_x) + sqr(grid_y) + sqr(grid_z));
          iYZ = ist.nx * (ist.ny * iZ + iY);
          iGrid = iYZ + iX;
          rho[iGrid] = psitot[iGrid+j*ist.ngrid] / exp(-dist_r);
          if (rho[iGrid]>100 || rho[iGrid]<-100) {rho[iGrid]=0;}
        }
      }
    }
  }

  sprintf(filename4, "psi%i_periodic.cube", j+n);

  writeCubeFile(rho, par, ist, filename4);
  ***/
  

  }

  
  
  /*************************************************************************/
  /*** free memeory ***/
  free(psitot); free(phi); free(psi); free(potl); free(eval); 
  free(el); free(ksqr); free(vx); free(vy); free(an); free(zn);
  free(vz); free(rx); free(ry); free(rz); free(sige); free(atm);

  fftw_destroy_plan(planfw);
  fftw_destroy_plan(planbw);
  fftw_free(fftwpsi);

  /*************************************************************************/
  /*** Print out the time that the computation finished ***/
  writeSeparation(stdout);
  writeCurrentTime(stdout);

  /*************************************************************************/
  /*** Exit the program ***/
  exit(0);
}
