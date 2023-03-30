#include "fd.h"

/****************************************************************************/
// prints the z-projected electron density (total, spin up, and spin down) for
// a complex wavefunction psi 

void print_sp_pz(zomplex *psi, double *sige, double *vz, par_st par, lng_st ist)
{ 
  FILE *pf;
  char str[100];
  long jx, jy, jz, jyz, jxyz, jns;
  double pz, pzUp, pzDown; // total density, spin up density, and spin down density

  for (jns = 0; jns <= ist.nhomo; jns++) {
    if (sige[jns] < par.deps) {
      sprintf (str,"pzv%ld.dat",jns);
      pf = fopen (str , "w"); 
      for (jz = 0; jz < ist.nz; jz++) {
        pz = pzUp = pzDown = 0.0;
        for (jy = 0; jy < ist.ny; jy++) {
          for (jyz = ist.nx * (ist.ny * jz + jy), jx = 0; jx < ist.nx; jx++) {
            jxyz = jyz + jx;
            pzUp += (sqr(psi[jns*ist.nspinngrid+jxyz].re) + sqr(psi[jns*ist.nspinngrid+jxyz].im));
            pzDown += (sqr(psi[jns*ist.nspinngrid+ist.ngrid+jxyz].re) + sqr(psi[jns*ist.nspinngrid+ist.ngrid+jxyz].im));      
          }
        }
        pz = pzUp+pzDown;
        fprintf (pf,"%.4f %.10f %.10f %.10f\n",vz[jz], pz*par.dx*par.dy, pzUp*par.dx*par.dy, pzDown*par.dx*par.dy);
      }
      fclose(pf);
    }
  }

  for (jns = ist.nlumo; jns < ist.ms; jns++) {
    if (sige[jns] < par.deps) {
      sprintf (str,"pzc%ld.dat",jns);
      pf = fopen (str , "w");
      for (jz = 0; jz < ist.nz; jz++) {
        pz = pzUp = pzDown = 0.0;
        for (jy = 0; jy < ist.ny; jy++) {
          for (jyz = ist.nx * (ist.ny * jz + jy), jx = 0; jx < ist.nx; jx++) {
            jxyz = jyz + jx;
            pzUp += (sqr(psi[jns*ist.nspinngrid+jxyz].re) + sqr(psi[jns*ist.nspinngrid+jxyz].im));
            pzDown += (sqr(psi[jns*ist.nspinngrid+ist.ngrid+jxyz].re) + sqr(psi[jns*ist.nspinngrid+ist.ngrid+jxyz].im));
          }
        }
        pz = pzUp+pzDown;
        fprintf (pf,"%.4f %.10f %.10f %.10f\n",vz[jz], pz*par.dx*par.dy, pzUp*par.dx*par.dy, pzDown*par.dx*par.dy);
      }
      fclose(pf);
    }
  }

  return;
}

/****************************************************************************/
// calculates the density of a complex wavefunction psi over the first 
// ngrid points each of volume dv and stores result in the density pointer


/****************************************************************************/

void print_pz_one(double *density, double *densityUp, double *densityDown, double *vz, par_st par, lng_st ist, char *str)
{
  FILE *pf; 
  long jx, jy, jz, jyz, jxyz;
  double pz, pzUp, pzDown;

  pf = fopen (str , "w");    
  for (jz = 0; jz < ist.nz; jz++) {
    for (jy = 0; jy < ist.ny; jy++) {
      pz = 0.0; pzUp = 0.0; pzDown = 0.0; 
      for (jyz = ist.nx * (ist.ny * jz + jy), jx = 0; jx < ist.nx; jx++) {
      	jxyz = jyz + jx;
        pzUp += densityUp[jxyz];
        pzDown += densityDown[jxyz];
      	pz += density[jxyz];
      }
    }
    fprintf (pf,"%.4f %.10f %.10f %.10f\n", vz[jz], pz*par.dx*par.dy, pzUp*par.dx*par.dy, pzDown*par.dx*par.dy);
  }
  fclose(pf);

  return;
}

/****************************************************************************/

void print_cube(double *pgrid, lng_st ist, par_st par, char *fName)
{
  FILE *pf, *fr; 
  int jgrid, jx, jy, jz, jyz, ityp;
  int ncub = 10;
  double x, y, z;
  char line[80], atyp[10];
 
  if(ist.ms2<ncub) ncub = ist.ms2; 
  
  fr = fopen("conf.par","r");
  pf = fopen(fName , "w");
  fprintf(pf,"CUBE FILE\n");
  fprintf(pf,"OUTER LOOP: Z, MIDDLE LOOP: Y, INNER LOOP: X\n");
  fprintf(pf,"%5li%12.6f%12.6f%12.6f\n",ist.natom,par.xmin,par.ymin,par.zmin);
  fprintf(pf,"%5li%12.6f%12.6f%12.6f\n",ist.nz,0.0,0.0,par.dz);
  fprintf(pf,"%5li%12.6f%12.6f%12.6f\n",ist.ny,0.0,par.dy,0.0);
  fprintf(pf,"%5li%12.6f%12.6f%12.6f\n",ist.nx,par.dx,0.0,0.0);
  fgets(line, 80, fr);
  while(fgets(line, 80, fr) != NULL) {
    sscanf (line,"%2s %lf %lf %lf",&atyp,&x,&y,&z);
    if(!strcmp(atyp,"Cd")) { 
      ityp=48;
    }
    else { 
      if(!strcmp(atyp,"S")) {  
        ityp=16;
      }
      else { 
        if(!strcmp(atyp,"Se")){ ityp=34;}
        else{ityp=1;}
      }
    }
    fprintf(pf,"%5i%12.6f%12.6f%12.6f%12.6f\n", ityp, 0.0, x, y, z);
  }
  
  for (jz = 0; jz < ist.nz; jz++) {
    for (jy = 0; jy < ist.ny; jy++) {
      jyz = ist.nx * (ist.ny * jz + jy);
      for (jx = 0; jx < ist.nx; jx++) {
      	jgrid = jyz + jx;
      	fprintf(pf,"%g ", pgrid[jgrid]);
      	if (jx % 6 == 5) fprintf(pf,"\n");
      }
      fprintf(pf,"\n");
    }
  }
  
  fclose(pf); fclose(fr);
  
  return;
}

/****************************************************************************/