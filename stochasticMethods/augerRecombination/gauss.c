/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "ar.h"

/*****************************************************************************/

void gauss_test(double *vx, double *vy, double *vz, zomplex *potq, double *poth, par_st par, lng_st ist,
                fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi)
{
  FILE *pf;
  zomplex *rho;
  long jx, jy, jz, jyz, jxyz;
  double x2, y2, z2, ax = 0.0, ay = 0.0, az = 0.0;
  double sigma = 3.0 * par.dx;

  //printf("ax = %g  ay = %g  az = %g   sigma = %g\n", ax, ay, az, sigma);
  
  // Dynamically allocated memory
  rho = (zomplex *) calloc(ist.nGridPoints, sizeof(zomplex));

  for (jz = 0; jz < ist.nz; jz++) {
    z2 = sqr((vz[jz]-az)/sigma);
    for (jy = 0; jy < ist.ny; jy++) {
      y2 = sqr((vy[jy]-ay)/sigma);
      jyz = ist.nx * (ist.ny * jz + jy);
      for (jx = 0; jx < ist.nx; jx++) {
      	x2 = sqr((vx[jx]-ax)/sigma);
      	jxyz = jyz + jx;
      	rho[jxyz].re = exp(-0.5*(x2 + y2 + z2));
      	rho[jxyz].im = 0.0;
      }
    }
  }
  
  hartree(rho,potq,poth,ist.nGridPoints,planfw,planbw,fftwpsi);

  pf = fopen("gauss.dat" , "w");
  for (jz = 0; jz < ist.nz; jz++) {
    z2 = sqr(vz[jz]);
    for (jy = 0; jy < ist.ny; jy++) {
      y2 = sqr(vy[jy]);
      jyz = ist.nx * (ist.ny * jz + jy);
      for (jx = 0; jx < ist.nx; jx++) {
      	x2 = sqr(vx[jx]);
      	jxyz = jyz + jx;
      	fprintf(pf, "%g %g\n", sqrt(x2+y2+z2), poth[jxyz]);
      }
    }
  }
  fclose(pf);
  
  // Free dynamically allocated memory
  free(rho);

  return;
}

/*****************************************************************************/
