#include "fd.h"

/****************************************************************************/
// Diagnolizes the BSE matrix, u, <Sai|Hbse|Sbj> to obtain 
// the excitonic eigenstates and eigenenergies 

void bethe_salpeter(double *bsmat, double *h0mat, zomplex *psi, double *vz,
                    zomplex *mux,zomplex *muy,zomplex *muz,lng_st ist,par_st par) {

  FILE *pf; char str[100]; long a, b, i, j, k, l, ibs, jbs, jgamma,jgrid;
  double *h, *u, *eval, *ev, *mat, sum, os, *pgrid, *pgridUp, *pgridDown;
  zomplex sumX, sumY, sumZ;
  
  mat = (double *) calloc(ist.ms2*ist.ms2, sizeof(double));
  h = (double *) calloc(ist.ms2*ist.ms2, sizeof(double));
  u = (double *) calloc(ist.ms2*ist.ms2, sizeof(double));
  eval = (double *) calloc(ist.ms2, sizeof(double));

  // calculate the full BSE matrix, u 
  // u = h0 - bsmat
  // h0 = Ec-Ev and bsmat = W - v 
#pragma omp parallel for private(i)
  for (i = 0; i < ist.ms2 * ist.ms2; i++) {
    h[i] = u[i] = h0mat[i] - bsmat[i];
  }

  // diagnolize the full BSE matrix, u and store
  // the eigenvectors in u and the eigenvalues in eval
  diag(ist.ms2, ist.nthreads, u, eval);

  // write out the prediagnolized full BSE matrix
  // HBSmat.dat that contains h0.dat added with bs.dat
  pf = fopen("HBSmat.dat", "w");
  for (i = 0; i < ist.ms2; i++, fprintf(pf, "\n"))
    for (j = 0; j < ist.ms2; j++)
      fprintf(pf,"%g ", h[i*ist.ms2+j]);
  fclose(pf);
  
  // write out the coefficients defining the lowest energy exciton
  pf = fopen("Uvec.dat", "w");
  for (i = 0; i < ist.ms2; i++)
    fprintf (pf,"%g\n",u[i]);
  fclose(pf);

  // calculate the energy of the lowest energy exciton
  for (sum = 0.0, i = 0; i < ist.ms2; i++)
    for (j = 0; j < ist.ms2; j++)
      sum += u[j] * u[i] * h[i*ist.ms2+j];
  printf ("The ground state exciton has energy = %.10f %.8f\n", sum, sum*AUTOEV);

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
  
  pf = fopen("test.dat", "w");
  for (i = 0; i < ist.ms2; i++, fprintf(pf,"\n"))
    for (j = 0; j < ist.ms2; j++)
      fprintf (pf,"%g ",h[i*ist.ms2+j]);
  fclose(pf);

  // print out all coefficients/ eigenvectors obtained from
  // diagnolizing the full BSE matrix
  pf = fopen("Umat.dat", "w");
  for (i = 0; i < ist.ms2; i++, fprintf(pf,"\n"))
    for (j = 0; j < ist.ms2; j++)
      fprintf (pf,"%g ",u[i*ist.ms2+j]);
  fclose(pf);
  
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

  // print the excitonic state energies
  pf = fopen("exciton.dat" , "w");
  for (i = 0; i < ist.ms2; i++)
    fprintf (pf,"%ld %.16g %.16g  %.16g  %.16g\n",i,eval[i],h[i*ist.ms2+i],h0mat[i*ist.ms2+i],
	     (eval[i]-h[i*ist.ms2+i])*27.2114);
  fclose(pf);

  // print the oscillator strengths
  pf = fopen("OS.dat" , "w");
  for (jgamma = 0; jgamma < ist.ms2; jgamma++) {
    sumX.re = sumX.im = sumY.re = sumY.im = sumZ.re = sumZ.im = 0.0;
    for (ibs = 0, a = ist.nlumo; a < ist.nlumo+ist.totallumo; a++) {
      for (i = 0; i < ist.totalhomo; i++, ibs++) {
        sumX.re += u[jgamma*ist.ms2+ibs] * mux[i*ist.totallumo+(a-ist.nlumo)].re;
        sumX.im += u[jgamma*ist.ms2+ibs] * mux[i*ist.totallumo+(a-ist.nlumo)].im;
        sumY.re += u[jgamma*ist.ms2+ibs] * muy[i*ist.totallumo+(a-ist.nlumo)].re;
        sumY.im += u[jgamma*ist.ms2+ibs] * muy[i*ist.totallumo+(a-ist.nlumo)].im;
        sumZ.re += u[jgamma*ist.ms2+ibs] * muz[i*ist.totallumo+(a-ist.nlumo)].re;
        sumZ.im += u[jgamma*ist.ms2+ibs] * muz[i*ist.totallumo+(a-ist.nlumo)].im;
      }
    } 
    os=(sqr(sumX.re)+sqr(sumX.im)+sqr(sumY.re)+sqr(sumY.im)+sqr(sumZ.re)+sqr(sumZ.im));
    fprintf (pf,"%g %g %g %g %g %g\n",sqrt(os),eval[jgamma],(4.0/3.0)*eval[jgamma]*os,
      sqr(sumX.re)+sqr(sumX.im),sqr(sumY.re)+sqr(sumY.im),sqr(sumZ.re)+sqr(sumZ.im));
  }
  fclose(pf);

  // calculate and print the electron and z-projected
  // electron densities of the excitonic states
  pgrid = (double*)calloc(ist.ngrid,sizeof(double));
  pgridUp = (double*)calloc(ist.ngrid,sizeof(double));
  pgridDown = (double*)calloc(ist.ngrid,sizeof(double));
  double densSpinUp, densSpinDown;
  for (jgamma = 0; jgamma < 10; jgamma++) {
    for (jgrid = 0; jgrid < ist.ngrid; jgrid++) {
      pgrid[jgrid] = 0.0;
      pgridUp[jgrid] = 0.0;
      pgridDown[jgrid] = 0.0;
    }
    for (ibs = 0, a = ist.nlumo; a < ist.nlumo+ist.totallumo; a++) {
      for (i = 0; i < ist.totalhomo; i++, ibs++) {
        for (jbs = 0, b = ist.nlumo; b < ist.nlumo+ist.totallumo; b++) {
      	  for (j = 0; j < ist.totalhomo; j++, jbs++) {
      	    if (i == j) { 
      	      sum = u[jgamma*ist.ms2+ibs] * u[jgamma*ist.ms2+jbs];
            #pragma omp parallel for private(jgrid)
              for (jgrid = 0; jgrid < ist.ngrid; jgrid++) {
                // TODO: make sure complex conjugates are correct as well
                densSpinUp = sum * (psi[a*ist.nspinngrid+jgrid].re * psi[b*ist.nspinngrid+jgrid].re
                  + psi[a*ist.nspinngrid+jgrid].im * psi[b*ist.nspinngrid+jgrid].im);
                densSpinDown = sum * (psi[a*ist.nspinngrid+ist.ngrid+jgrid].re * psi[b*ist.nspinngrid+ist.ngrid+jgrid].re
                  + psi[a*ist.nspinngrid+ist.ngrid+jgrid].im * psi[b*ist.nspinngrid+ist.ngrid+jgrid].im);
		            // original: pgrid[jgrid] += sum * psi[a*ist.ngrid+jgrid] * psi[b*ist.ngrid+jgrid];
                pgridUp[jgrid] += densSpinUp;
                pgridDown[jgrid] += densSpinDown;
                pgrid[jgrid] += (densSpinUp+densSpinDown); 
              } 
            }
      	  }
      	}
      }
    }
    sprintf (str,"pcz%d-bs.dat",jgamma);
    print_pz_one(pgrid,pgridUp,pgridDown,vz,par,ist,str);

    sprintf(str,"pe-%d.cub",jgamma);
    print_cube(pgrid,ist,par,str);
    
    for (jgrid = 0; jgrid < ist.ngrid; jgrid++) {
      pgrid[jgrid] = 0.0;
      pgridUp[jgrid] = 0.0;
      pgridDown[jgrid] = 0.0;
    }
    for (ibs = 0, a = ist.nlumo; a < ist.nlumo+ist.totallumo; a++) {
      for (i = 0; i < ist.totalhomo; i++, ibs++) {
      	for (jbs = 0, b = ist.nlumo; b < ist.nlumo+ist.totallumo; b++) {
      	  for (j = 0; j < ist.totalhomo; j++, jbs++) {
      	    if (a == b) { 
      	      sum = u[jgamma*ist.ms2+ibs] * u[jgamma*ist.ms2+jbs];
            #pragma omp parallel for private(jgrid)
      	      for (jgrid = 0; jgrid < ist.ngrid; jgrid++) {
                // TODO: make sure complex conjugates are correct as well
                densSpinUp = sum * (psi[i*ist.nspinngrid+jgrid].re * psi[j*ist.nspinngrid+jgrid].re
                  + psi[i*ist.nspinngrid+jgrid].im * psi[j*ist.nspinngrid+jgrid].im);
                densSpinDown = sum * (psi[i*ist.nspinngrid+ist.ngrid+jgrid].re * psi[j*ist.nspinngrid+ist.ngrid+jgrid].re
                  + psi[i*ist.nspinngrid+ist.ngrid+jgrid].im * psi[j*ist.nspinngrid+ist.ngrid+jgrid].im);
      		      // original: pgrid[jgrid] += sum * psi[i*ist.ngrid+jgrid].re * psi[j*ist.ngrid+jgrid].re;
                pgridUp[jgrid] += densSpinUp;
                pgridDown[jgrid] += densSpinDown;
                pgrid[jgrid] += (densSpinUp+densSpinDown);              
              }
            }
      	  }
      	}
      }
    }
    sprintf (str,"pvz%d-bs.dat",jgamma);
    print_pz_one(pgrid,pgridUp,pgridDown,vz,par,ist,str);

    sprintf(str,"ph-%d.cub",jgamma);
    print_cube(pgrid,ist,par,str);
  }
  
  free(pgrid); free(pgridUp); free(pgridDown); 
  free(u); free(h); free(eval); free(mat);

  return;
}

/****************************************************************************/