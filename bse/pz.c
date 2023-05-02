/****************************************************************************/

#include "fd.h"
#include <float.h>

/****************************************************************************/

int z_project(double *vector, double *vz, par_st par, long_st ist, char *fname)
/* This function prints the z-projection file and returns whether or not 
 * the state was surface trapped or not
 *
 * psi vector has form: 
 * v[(x0,y0,z0)..(xN-1,y0,z0)..(x0..y1..z0)..(xN-1,yN-1,z0)..(xN-1,yN-1,zN-1)]
 */
{
    FILE *fp;
    long i, j, k, xyz = 0;
    long nx = ist.nx, 
         ny = ist.ny, 
         nz = ist.nz;
    double dx = par.dx,
           dy = par.dy,
           pz;
    struct pzdata *zproj; 
    int flag = 0;

    if ((zproj = calloc(nz, sizeof(struct pzdata))) == NULL) {
        fprintf(stderr, "Memory failure in z_project\n");
        exit(EXIT_FAILURE);
    }
    fp = fopen(fname, "w");
    if (fp) {
        for (k = 0; k < nz; k++) {
            pz = 0.0;
            for (j = 0; j < ny; j++) {
                for (i = 0; i < nx; i++) {
                    xyz = nx * (ny * k + j) + i;
                    pz += vector[xyz] * vector[xyz];
                }
            }
            zproj[k].z = vz[k];
            zproj[k].pz = pz * dx * dy;
            fprintf(fp, "%.*g %.*g\n", DBL_DIG, zproj[k].z, DBL_DIG, zproj[k].pz);
        }
        flag = istrapped(zproj, nz);
    } else {
        free(zproj);
        fprintf(stderr, "Failed in z_project\n");
        exit(EXIT_FAILURE);
    }
    fclose(fp);
    free(zproj);
    if (flag) {
        remove(fname);
    }

    return flag;
}

/****************************************************************************/

int istrapped(struct pzdata *data, long len)
/* This is the same logic I used in my python scripts to select the state for 
 * psi-hl.x.
 *
 * NOTE: There is a slight difference that I've noticed in that trapped is
 * defaulted to false here where in psi-hl.py it is defaulted to true.
 *
 * Also changed abs to fabs. I wish the compiler would tell me when I'm doing
 * dumb things.
 */
{
    int i;
    double mean, var, stddev, norm, 
           xi_pi, xisq_pi, abs_xi, 
           maxz, maxpz;

    mean = var = stddev = norm = 0; 
    xi_pi = xisq_pi = abs_xi = 0;

    maxz = data[0].z;
    maxpz = data[0].pz;
    for (i = 0; i < len; i++) {
        if (data[i].pz > maxpz) {
            maxpz = data[i].pz;
            maxz = data[i].z;
        }
        norm += data[i].pz;
        xi_pi += data[i].z * data[i].pz;
        xisq_pi += data[i].z * data[i].z * data[i].pz;
        abs_xi += abs(data[i].z);
    }
    mean = xi_pi / norm;
    var = (xisq_pi / norm) - (mean * mean);
    stddev = sqrt(var);

    if (fabs(mean) >= (fabs(data[0].z) / 2.0)) {
        return 0; // Should return 1 since trapped but not using this feature: John - Aug 28 2018
        //return 1;
    }

    return 0;
}


/****************************************************************************/

void print_cube(double *pgrid, long_st ist, par_st par, char *fName) {
  FILE *pf, *fr; 
  int i, j, k, l, jgamma,jgrid,jx, jy, jz, jyz,jxyz,ityp;
  int ncub = 10;
  double dz, dx, dy,x,y,z;
  char line[80],atyp[10];
  if(ist.ms2<ncub) ncub = ist.ms2; 
  
  fr = fopen("conf.par", "r");
  pf = fopen(fName , "w");
  fprintf(pf, "CUBE FILE\n");
  fprintf(pf, "OUTER LOOP: Z, MIDDLE LOOP: Y, INNER LOOP: X\n");
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", ist.natom, par.xmin, par.ymin, par.zmin);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", ist.nz, 0.0, 0.0, par.dz);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", ist.ny, 0.0, par.dy, 0.0);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", ist.nx, par.dx, 0.0, 0.0);
  fgets(line, 80, fr);
  while(fgets(line, 80, fr) != NULL) {
    sscanf(line, "%2s %lf %lf %lf", &atyp, &x, &y, &z);
    if (! strcmp(atyp,"Cd")) { 
      ityp=48;
    }
    else if (! strcmp(atyp,"S")) {  
	  ityp=16;
    }
    else if (! strcmp(atyp,"Se")) { 
		ityp=34;
	}
	else if (! strcmp(atyp,"Zn")) {
		ityp=30;
	}
	else if (! strcmp(atyp,"Te")) {
		ityp=52;
	}
	else { 
		ityp=1; 
	}
    fprintf(pf, "%5i%12.6f%12.6f%12.6f%12.6f\n", ityp, 0.0, x, y, z);
  }
  for (jz = 0; jz < ist.nz; jz++) {
    for (jy = 0; jy < ist.ny; jy++) {
      jyz = ist.nx * (ist.ny * jz + jy);
      for (jx = 0; jx < ist.nx; jx++) {
	    jgrid = jyz + jx;
	    fprintf(pf, "%g ", pgrid[jgrid]);
	    if (jx % 6 == 5) {
		  fprintf(pf, "\n");
		}
      }
      fprintf(pf, "\n");
    }
  }
  fclose(pf); fclose(fr);
  
  return;
}

/****************************************************************************/
// TODO: write below function
// void print_all_axes_projections()

/****************************************************************************/

void print_pz_one(double *psi2, double *vz, par_st par, long_st ist, char *str) {
  FILE *pf; 
  long jx, jy, jz, jyz, jxyz;
  double pz;

  pf = fopen(str, "w");    
  for (jz = 0; jz < ist.nz; jz++) {
    for (pz = 0.0, jy = 0; jy < ist.ny; jy++) {
      for (jyz = ist.nx * (ist.ny * jz + jy), jx = 0; jx < ist.nx; jx++) {
        jxyz = jyz + jx;
        pz += psi2[jxyz];
      }
    }
    fprintf(pf, "% .4f % .8f\n", vz[jz], pz*par.dx*par.dy);
  }
  fclose(pf);

  return;
}

/****************************************************************************/
// Prints cube and z-projection files for the num_excitons=1 lowest energy excitons 
// for which the other quasiparticle is fixed to a small volume around the origin

void print_fixed_qp_density(double *psi, double *Cbs, double *vz, long_st ist, par_st par) {
  int jx, jy, jz, jyz, jgrid, fixedGrid;
  int a, i, b, j, ibs, jbs, exc_index;
  int ilist, *list;
  double x, y, z;
  double coeff, *psie, *psih, *pegrid, *phgrid;
  char fileName[100];

  // Print function initiation
  writeSeparation(stdout);
  writeCurrentTime(stdout);
  printf("\nPrinting quasiparticle densities with other particle at a fixed location:\n");

  // Parameters
  int num_excitons = 1;
  double origin_cutoff = 2.0;
  int numTestPoints= 1;
  int numListPoints[numTestPoints];
  double xTest[numTestPoints], yTest[numTestPoints], zTest[numTestPoints];

  // Dynamically allocate memory
  psih = (double *) calloc(ist.ngrid, sizeof(double));
  psie = (double *) calloc(ist.ngrid, sizeof(double));
  phgrid = (double *) calloc(ist.ngrid*numTestPoints, sizeof(double));
  pegrid = (double *) calloc(ist.ngrid*numTestPoints, sizeof(double));
  list = (int *) calloc(ist.ngrid*numTestPoints, sizeof(int));

  // The test points
  xTest[0] = 0.0; yTest[0] = 0.0; zTest[0] = 0.0; // the origin: (0, 0, 0)
  //xTest[1] = 0.0; yTest[1] = 0.0; zTest[1] = -0.75*par.zmin; // (0, 0, -0.75*zmin)
  //xTest[2] = 0.0; yTest[2] = -0.75*par.ymin; zTest[2] = 0.0; // (0, -0.75*ymin, 0)
  //xTest[3] = -0.75*par.xmin; yTest[3] = 0.0; zTest[3] = 0.0; // (-0.75*xmin, 0, 0)

  // create list of points within origin_cutoff of 
  // TODO: make the points based on the shape of the NC
  int totalNumListPoints = 0;
  for (i = 0; i < numTestPoints; i++) {
    numListPoints[i] = 0;
	  for (z = par.zmin, jz = 0; jz < ist.nz; jz++) {
      z += par.dz;
      for (y = par.ymin, jy = 0; jy < ist.ny; jy++) {
        y += par.dy;
        jyz = ist.nx * (ist.ny * jz + jy);
        for (x = par.xmin, jx = 0; jx < ist.nx; jx++) {
          x += par.dx;
          // Calculate distance from each of the test points and check if
          // the grid point is within a sphere surrounding the test point 
          if ((sqr((z-zTest[i]))+sqr((y-yTest[i]))+sqr((x-xTest[i]))) > origin_cutoff) { 
            continue;
          }
          else {
            list[totalNumListPoints] = jyz + jx;
            numListPoints[i] += 1;
			      totalNumListPoints += 1;
          }
        }
      }
    }
  }
  printf("The number of fixed points = %d\n", numTestPoints);
  for (i = 0; i < numTestPoints; i++) {
	printf("Fixed point %d:\n", i);
    printf("Fixed point location = % .5f % .5f % .5f\n", xTest[i], yTest[i], zTest[i]);
    printf("Number of points within sphere near the fixed point = %d\n", numListPoints[i]);
  }
  fflush(stdout);

  // Calculate pegrid and phgrid for the lowest num_excitons 
  int testPointIndex = 0;
  for (exc_index = 0; exc_index < num_excitons; exc_index++) {
    for (jgrid = 0; jgrid < ist.ngrid*numTestPoints; jgrid++) {
      pegrid[jgrid] = 0.0;
      phgrid[jgrid] = 0.0;
    }  
    // Now loop over states involved in the BSE calculation 
    for (ibs = 0, a = ist.nlumo; a < ist.nlumo+ist.totallumo; a++) {
      for (i = 0; i < ist.totalhomo; i++, ibs++) {
        for (jbs = 0, b = ist.nlumo; b < ist.nlumo+ist.totallumo; b++) {
          for (j = 0; j < ist.totalhomo; j++, jbs++) {
            coeff = Cbs[exc_index*ist.ms2 + ibs] * Cbs[exc_index*ist.ms2 + jbs];
            // Store cia*cjb*psii*psij and cia*cjb*psia*psib loop over all grid points
            // Can do these with just loops over just i and j (a and b) separately
            for (jgrid = 0; jgrid < ist.ngrid; jgrid++) {
              psih[jgrid] = coeff * psi[i*ist.ngrid + jgrid] * psi[j*ist.ngrid + jgrid];
              psie[jgrid] = coeff * psi[a*ist.ngrid + jgrid] * psi[b*ist.ngrid + jgrid];
            }
            // TODO: add loop over numTestPoints
			      totalNumListPoints = 0;
			      // for (testPointIndex = 0; testPointIndex < numTestPoints; testPointIndex++) {
            for (fixedGrid = 0; fixedGrid < numListPoints[testPointIndex]; fixedGrid++) {
              ilist = list[totalNumListPoints];
			  totalNumListPoints++;
            #pragma omp parallel for private(jgrid)
              for (jgrid = 0; jgrid < ist.ngrid; jgrid++) {
                phgrid[jgrid + testPointIndex*ist.ngrid] += psie[ilist] * psi[i*ist.ngrid + jgrid] * psi[j*ist.ngrid + jgrid];
                pegrid[jgrid + testPointIndex*ist.ngrid] += psih[ilist] * psi[a*ist.ngrid + jgrid] * psi[b*ist.ngrid + jgrid];
              }
            }
			      //}
          }
        }
      }
    }
    norm_vector(&pegrid[testPointIndex*ist.ngrid], par.dv, ist.ngrid);
    norm_vector(&phgrid[testPointIndex*ist.ngrid], par.dv, ist.ngrid);
    sprintf(fileName, "pe-fp%d-%d.cub", testPointIndex, exc_index);
    print_cube(&pegrid[testPointIndex*ist.ngrid], ist, par, fileName);
    sprintf(fileName, "pcz%d-fp%d.dat", exc_index, testPointIndex);
    print_pz_one(&pegrid[testPointIndex*ist.ngrid], vz, par, ist, fileName);
    sprintf(fileName, "ph-fp%d-%d.cub", testPointIndex, exc_index);
    print_cube(&phgrid[testPointIndex*ist.ngrid], ist, par, fileName);  
    sprintf(fileName, "pvz%d-fp%d.dat", exc_index, testPointIndex);
    print_pz_one(&phgrid[testPointIndex*ist.ngrid], vz, par, ist, fileName);
  }

  // Free dynamically allocated memory
  free(psih); free(psie);
  free(pegrid); free(phgrid);
  free(list); 

  return;
}

/****************************************************************************/
