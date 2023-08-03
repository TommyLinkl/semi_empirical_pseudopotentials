/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "fd.h"

/****************************************************************************/

xyz_st calc_dipole(long i,long a,double *vx,double *vy,double *vz,double *psi,long_st ist,par_st par)
{
  long jx, jy, jz, jyz, jxyz;
  xyz_st dip;
  double rhoia;

  dip.x = dip.y = dip.z = 0.0;
  for (jz = 0; jz < ist.nz; jz++) {
    for (jy = 0; jy < ist.ny; jy++) {
      for (jyz = ist.nx * (ist.ny * jz + jy), jx = 0; jx < ist.nx; jx++) {
      	jxyz = jyz + jx;
      	rhoia = psi[i*ist.ngrid+jxyz] * psi[a*ist.ngrid+jxyz];
      	dip.x += vx[jx] * rhoia;
      	dip.y += vy[jy] * rhoia;
      	dip.z += vz[jz] * rhoia;
      }
    }
  }
  dip.x *= par.dv;
  dip.y *= par.dv;
  dip.z *= par.dv;

  return (dip);
}

/****************************************************************************/