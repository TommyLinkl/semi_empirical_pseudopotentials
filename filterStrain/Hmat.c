/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "fd.h"

/*****************************************************************************/

void Hmatreal(zomplex *psi,zomplex *phi,double *psitot,double *potl,double *ksqr,double *eval,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  FILE *pg;
  long ims, jms, jgrid;
  long long mstot = (long long)ist.mstot, info, lwk = 3*ist.mstot;
  double *H, *work, *tpsi, sum;

  // Dynamically allocate memory
  H = (double *) calloc(ist.mstot*ist.mstot, sizeof(double));
  work = (double *) calloc(lwk, sizeof(double));
  tpsi = (double *) calloc(ist.mstot, sizeof(double));

  // TODO: use openmp
  //omp_set_dynamic(0);
  //omp_set_num_threads(ist.nthreads);

  // Calculate H|psi_i> 
  pg = fopen("hmat.dat" , "w");
  fprintf(pg, "calculate the H matrix\n"); fflush(pg);
  for (ims = 0; ims < ist.mstot; ims++) {
    //#pragma omp parallel for private(jgrid)
    for (jgrid = 0; jgrid < ist.ngrid; jgrid++) {
      psi[jgrid].re = psitot[ims*ist.ngrid + jgrid];
      psi[jgrid].im = 0.0;
    }
    memcpy(&phi[0], &psi[0], ist.ngrid*sizeof(phi[0]));
    hamiltonian(phi, psi, potl, ksqr, ist, par, planfw, planbw, fftwpsi);

    // Calculate <psi_j|H|psi_i> 
    //#pragma omp parallel for private(jms)
    for (jms = 0; jms < ist.mstot; jms++) {
      H[ims*ist.mstot + jms] = dotpreal(phi, psitot, 0, jms, ist.ngrid, par.dv);
      //fprintf(pg, "%ld %ld %g\n", ims, jms, H[ims*ist.mstot+jms]);
    }
    fprintf(pg, "finshed row %ld\n", ims); fflush(pg);
  }

  // Diagonalize the Hamiltonian, H
  // TODO: make sure/ make dysev_ utilizes all physical cores
  // mkl_set_dynamic(0); 
  // mkl_set_num_threads(ist.nthreads);
  // omp_set_nested(int val) -> enables nested parallelism if val != 0 and disables it if val == 0
  // omp_set_nested(1); // enable nested parallelism 
  fprintf(pg, "diagonalize the new Hamiltonian H\n"); fflush(pg);
  dsyev_("V", "U", &mstot, &H[0], &mstot, &eval[0], &work[0], &lwk, &info);
  if (info) { 
    nerror("error in dsyev_ H");
  }    
  // Copy the new function into psitot 
  fprintf(pg, "copy the new function into psitot\n"); fflush(pg);
  for (jgrid = 0; jgrid < ist.ngrid; jgrid++) { 
    for (jms = 0; jms < ist.mstot; jms++) {
      tpsi[jms] = psitot[jms*ist.ngrid + jgrid];
      psitot[jms*ist.ngrid + jgrid] = 0.0;
    }
    for (jms = 0; jms < ist.mstot; jms++) {
      sum = 0.0; 
      //#pragma omp parallel for reduction(+:sum)
      for (ims = 0; ims < ist.mstot; ims++) {
	    sum += H[jms*ist.mstot + ims]*tpsi[ims];
      }
      //#pragma omp critical
      psitot[jms*ist.ngrid+jgrid] = sum;
    }
    if (!(jgrid % 1000)) { 
	  fprintf(pg, "finished grid point %ld\n", jgrid); fflush(pg);
    }
  }
  fprintf(pg, "free memory\n");
  fclose(pg);
  
  // Free dynamically allocated memory
  free(work); free(H); free(tpsi);

  return;
}

/*****************************************************************************/

double dotpreal(zomplex *psi,double *phi,long n,long m,long ngrid,double dv)
{
  long i;
  double sum;

  for (sum = 0.0, i = 0; i < ngrid; i++)
    sum += (psi[n*ngrid+i].re * phi[m*ngrid+i]);
  sum *= dv;
  return(sum);
}

/*****************************************************************************/
