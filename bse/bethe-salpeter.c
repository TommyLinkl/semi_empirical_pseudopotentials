/*****************************************************************************/

#include <float.h>
#include "fd.h"

/*****************************************************************************/

void bethe_salpeter(double *bsmat, double *h0mat, double *psi, double *vz, double *mux, double *muy, double * muz,
					double *mx, double *my, double *mz, long_st ist, par_st par)
{
  FILE *pf, *pf1, *pf2; 
  long a, b, i, j, k, l, ibs, jbs, jgamma, jgrid;
  double *h, *u, *eval, *ev, *mat, sum, sumx, sumy, sumz, os, msumx, msumy, msumz, mos, *pgrid;
  char str[100]; 

  mat = (double *) calloc(ist.ms2*ist.ms2, sizeof(double));
  h = (double *) calloc(ist.ms2*ist.ms2, sizeof(double));
  u = (double *) calloc(ist.ms2*ist.ms2, sizeof(double));
  eval = (double *) calloc(ist.ms2, sizeof(double));

  printf("The number of electron-hole pairs in the exciton basis = %d\n", ist.ms2);

#pragma omp parallel for private(i)
  for (i = 0; i < ist.ms2*ist.ms2; i++) {
    h[i] = u[i] = h0mat[i] - bsmat[i];
  }  
  diag(ist.ms2, ist.nthreads, u, eval);

  pf = fopen("HBSmat.dat", "w");
  for (i = 0; i < ist.ms2; i++, fprintf(pf,"\n")) {
    for (j = 0; j < ist.ms2; j++) {
      fprintf (pf,"%.*g ", DBL_DIG, h[i*ist.ms2+j]);
	}
  }
  fclose(pf);
  
  // Prints the coefficients for the 100 (or ist.ms2) lowest energy excitonic states
  long numExcStatesToPrint = ist.ms2; //100;
  if (ist.ms2 < 100) numExcStatesToPrint = ist.ms2;
  pf = fopen("BSEcoeff.dat", "w");
  for (j = 0; j < numExcStatesToPrint; j++) { 
    for (i = 0; i < ist.ms2; i++) {
	  fprintf(pf, "%.*g\n", DBL_DIG, u[i+j*ist.ms2]);
    } 	
  }
  fclose(pf);

  /*** lowest exciton energy ***/
  for (sum = 0.0, i = 0; i < ist.ms2; i++) {
    for (j = 0; j < ist.ms2; j++) {
      sum += u[j] * u[i] * h[i*ist.ms2+j];
    }
  }
  printf("Ground state exciton has energy = %.10f %.10f\n", sum, sum*AUTOEV);
  
#pragma omp parallel for private(l,j,k,sum)
  for (l = 0; l < ist.ms2; l++) {
    for (j = 0; j < ist.ms2; j++) {
      for (sum = 0, k = 0; k < ist.ms2; k++) {
      	sum +=  h[l*ist.ms2+k] * u[j*ist.ms2+k];
      }
      mat[l*ist.ms2+j] = sum;
    }
  }

#pragma omp parallel for private(i,j,k,sum)
  for (i = 0; i < ist.ms2; i++) {
    for (j = 0; j < ist.ms2; j++) {
      for (sum = 0, l = 0; l < ist.ms2; l++) {
      	sum +=  u[i*ist.ms2+l] * mat[l*ist.ms2+j];
      }
      h[i*ist.ms2+j] = sum;
    }
  }

#pragma omp parallel for private(l,j,k,sum)
  for (l = 0; l < ist.ms2; l++) {
    for (j = 0; j < ist.ms2; j++) {
      for (sum = 0, k = 0; k < ist.ms2; k++) {
        sum +=  h0mat[l*ist.ms2+k] * u[j*ist.ms2+k];
      }
      mat[l*ist.ms2+j] = sum;
    }
  }

#pragma omp parallel for private(i,l,sum)
  for (i = 0; i < ist.ms2; i++) {
    for (sum = 0, l = 0; l < ist.ms2; l++) {
      sum +=  u[i*ist.ms2+l] * mat[l*ist.ms2+i];
    }
    h[i*ist.ms2+i] = sum;
  }

  // Print out the energies of the excitonic states
  pf = fopen("exciton.dat" , "w");
  for (i = 0; i < ist.ms2; i++) {  
    fprintf(pf,"%ld % .12f % .12f  % .12f  % .12f\n", i, eval[i], h[i*ist.ms2+i], h0mat[i*ist.ms2+i],
	     (eval[i]-h[i*ist.ms2+i])*AUTOEV);
  }
  fclose(pf);

  // Calculate and print the electric and magnetic dipole strengths and the rotational strength
  pf = fopen("OS.dat", "w");
  pf1 = fopen("M.dat", "w");
  pf2 = fopen("rs.dat", "w");
  for (jgamma = 0; jgamma < ist.ms2; jgamma++) {
    sumx = 0.0; sumy = 0.0; sumz = 0.0;
    msumx = 0.0; msumy = 0.0; msumz = 0.0;
    for (ibs = 0, a = ist.nlumo; a < ist.nlumo+ist.totallumo; a++) {
      for (i = 0; i < ist.totalhomo; i++, ibs++) {
        sumx += u[jgamma*ist.ms2 + ibs] * mux[i*ist.totallumo + (a - ist.nlumo)];
        sumy += u[jgamma*ist.ms2 + ibs] * muy[i*ist.totallumo + (a - ist.nlumo)];
        sumz += u[jgamma*ist.ms2 + ibs] * muz[i*ist.totallumo + (a - ist.nlumo)];
        msumx += u[jgamma*ist.ms2 + ibs] * mx[i*ist.totallumo + (a - ist.nlumo)];
		    msumy += u[jgamma*ist.ms2 + ibs] * my[i*ist.totallumo + (a - ist.nlumo)];
		    msumz += u[jgamma*ist.ms2 + ibs] * mz[i*ist.totallumo + (a - ist.nlumo)];
      }
    } 
    os  = (sumx*sumx + sumy*sumy + sumz*sumz);
    mos = (msumx*msumx + msumy*msumy + msumz*msumz);
    fprintf(pf,  "%ld %.8f %.8f % .8f % .12f % .12f % .12f\n", jgamma, sqrt(os), eval[jgamma], (4.0/3.0)*eval[jgamma]*os, 
      				sumx*sumx, sumy*sumy, sumz*sumz);
    fprintf(pf1, "%ld %.8f %.8f % .8f % .12f % .12f % .12f\n", jgamma, sqrt(mos), eval[jgamma], (4.0/3.0)*eval[jgamma]*mos,
					    msumx*msumx, msumy*msumy, msumz*msumz);
    fprintf(pf2, "%ld %.8f % .16f\n", jgamma, eval[jgamma], (sumx*msumx + sumy*msumy + sumz*msumz));
  }
  fclose(pf); fclose(pf1); fclose(pf2);

  // Electron and hole carrier densities 
  pgrid = (double*)calloc(ist.ngrid,sizeof(double));
  for (jgamma = 0; jgamma < 1; jgamma++) {
	// Electron densities
    for (jgrid = 0; jgrid < ist.ngrid; jgrid++) { 
	  pgrid[jgrid] = 0.0;
	}
    for (ibs = 0, a = ist.nlumo; a < ist.nlumo+ist.totallumo; a++) {
      for (i = 0; i < ist.totalhomo; i++, ibs++) {
      	for (jbs = 0, b = ist.nlumo; b < ist.nlumo+ist.totallumo; b++) {
      	  for (j = 0; j < ist.totalhomo; j++, jbs++) {
      	    if (i == j) { 
      	      sum = u[jgamma*ist.ms2+ibs] * u[jgamma*ist.ms2+jbs];
            #pragma omp parallel for private(jgrid)
      	      for (jgrid = 0; jgrid < ist.ngrid; jgrid++)
            		pgrid[jgrid] += sum * psi[a*ist.ngrid+jgrid] * psi[b*ist.ngrid+jgrid];
      	    }
      	  }
      	}
      }
    }
    sprintf(str, "pcz%d-bs.dat", jgamma);
    print_pz_one(pgrid, vz, par, ist, str);
    sprintf(str, "pe-bs-%d.cub", jgamma);
    print_cube(pgrid, ist, par, str);
    
	// Hole densities
    for (jgrid = 0; jgrid < ist.ngrid; jgrid++) {
	  pgrid[jgrid] = 0.0;
    }
	for (ibs = 0, a = ist.nlumo; a < ist.nlumo+ist.totallumo; a++) {
      for (i = 0; i < ist.totalhomo; i++, ibs++) {
      	for (jbs = 0, b = ist.nlumo; b < ist.nlumo+ist.totallumo; b++) {
      	  for (j = 0; j < ist.totalhomo; j++, jbs++) {
      	    if (a == b) { 
      	      sum = u[jgamma*ist.ms2+ibs] * u[jgamma*ist.ms2+jbs];
            #pragma omp parallel for private(jgrid)
      	      for (jgrid = 0; jgrid < ist.ngrid; jgrid++)
            		pgrid[jgrid] += sum * psi[i*ist.ngrid+jgrid] * psi[j*ist.ngrid+jgrid];
      	    }
      	  }
      	}
      }
    }
    sprintf(str, "pvz%d-bs.dat", jgamma);
    print_pz_one(pgrid, vz, par, ist, str);
    sprintf(str, "ph-bs-%d.cub", jgamma);
    print_cube(pgrid, ist, par, str);
  }

  // Calculate pegrid and phgrid for the lowest 15 electron and hole state, respectively 
  // this is a noninteracting result -> test locations does not influence the result
  long nNonIntStates = 20;
  if (ist.totalhomo < nNonIntStates) nNonIntStates = ist.totalhomo;
  if (ist.totallumo < nNonIntStates) nNonIntStates = ist.totallumo;
  for (i = 0; i < nNonIntStates; i++) {
  #pragma omp parallel for private(jgrid)
    for (jgrid = 0; jgrid < ist.ngrid; jgrid++) { // electron first
      pgrid[jgrid] = par.dv *  psi[(ist.nlumo+i)*ist.ngrid + jgrid] * psi[(ist.nlumo+i)*ist.ngrid + jgrid];
    }
    norm_vector(pgrid, par.dv, ist.ngrid);
    sprintf(str, "pe-ni-%d.cub", i);
    print_cube(pgrid, ist, par, str);
    sprintf(str, "pcz%d-ni.dat", i); 
    print_pz_one(pgrid, vz, par, ist, str);
  #pragma omp parallel for private(jgrid)
    for (jgrid = 0; jgrid < ist.ngrid; jgrid++) { // hole second
      pgrid[jgrid] = par.dv *  psi[(ist.nlumo-1-i)*ist.ngrid + jgrid] * psi[(ist.nlumo-1-i)*ist.ngrid + jgrid];
    }
    norm_vector(pgrid, par.dv, ist.ngrid);
    sprintf(str, "ph-ni-%d.cub", i);
    print_cube(pgrid, ist, par, str);
    sprintf(str, "pvz%d-ni.dat", i); 
    print_pz_one(pgrid, vz, par, ist, str);
  }  

  // Calculate the pegrid and phgrid for the lowest excitonic state 
  // for fixed locations of the other carrier - an interacting result
  if (ist.printFPDensity) {
    print_fixed_qp_density(psi, u, vz, ist, par);
  } 
  
  free(pgrid);  free(u); free(h); free(eval); free(mat);

  return;
}

/*****************************************************************************/
