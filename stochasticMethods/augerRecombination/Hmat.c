#include "ar.h"

void Hmatreal(double *psitot,double *potl,double *ksqr,double *eval,lng_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  FILE *pg;
  long ims, jms, jgrid;
  long long ms = (long long)ist.nFilteredStates, info, lwk = 3*ist.nFilteredStates;
  double *H, *work, *tpsi, sum;
  zomplex *psi,*phi;
  
  if ((psi = (zomplex*)calloc(ist.nGridPoints,sizeof(zomplex)))==NULL)nerror("psi");
  if ((phi = (zomplex*)calloc(ist.nGridPoints,sizeof(zomplex)))==NULL)nerror("phi");

  H = (double*)calloc(ist.nFilteredStates*ist.nFilteredStates,sizeof(double));
  work = (double*)calloc(lwk,sizeof(double));
  tpsi = (double*)calloc(ist.nFilteredStates,sizeof(double));

  /*omp_set_dynamic(0);
    omp_set_num_threads(ist.nThreads);*/

  /*** calculate H|psi_i> ***/
  pg = fopen("hmat.dat" , "w");
  fprintf (pg,"calculate the H matrix\n"); fflush(pg);
  for (ims = 0; ims < ist.nFilteredStates; ims++){
    //#pragma omp parallel for private(jgrid)
    for (jgrid = 0; jgrid < ist.nGridPoints; jgrid++) {
      psi[jgrid].re = psitot[ims*ist.nGridPoints+jgrid];
      psi[jgrid].im = 0.0;
    }
    memcpy(&phi[0],&psi[0],ist.nGridPoints*sizeof(phi[0]));
    hamiltonian(phi,psi,potl,ksqr,ist,par,planfw,planbw,fftwpsi);

    /*** calculate <psi_j|H|psi_i> ***/
    //#pragma omp parallel for private(jms)
    for (jms = 0; jms < ist.nFilteredStates; jms++){
      H[ims*ist.nFilteredStates+jms] = dotpreal(phi,psitot,0,jms,ist.nGridPoints,par.dv);
      //fprintf (pg,"%ld %ld %g\n",ims,jms,H[ims*ist.nFilteredStates+jms]);
    }
    fprintf (pg,"finshed row %ld\n",ims); fflush(pg);
  }

  /*** diagonalize the Hamiltonian H ***/
  fprintf (pg,"diagonalize the new Hamiltonian H\n"); fflush(pg);
  dsyev_("V","U",&ms,&H[0],&ms,&eval[0],&work[0],&lwk,&info);
  if (info) nerror("error in dsyev_ H");
    
  /*** copy the new function into psitot ***/
  fprintf (pg,"copy the new function into psitot\n"); fflush(pg);

  for (jgrid = 0; jgrid < ist.nGridPoints; jgrid++){
#pragma omp parallel for private (jms)
    for (jms = 0; jms < ist.nFilteredStates; jms++){
      tpsi[jms] = psitot[jms*ist.nGridPoints+jgrid];
      psitot[jms*ist.nGridPoints+jgrid] = 0.0;
    }

#pragma omp parallel for private (jms,sum)
    for (jms = 0; jms < ist.nFilteredStates; jms++) {
      sum = 0.0; 
      for (ims = 0; ims < ist.nFilteredStates; ims++)
	sum += H[jms*ist.nFilteredStates+ims]*tpsi[ims];

      //#pragma omp critical
      psitot[jms*ist.nGridPoints+jgrid] = sum;
    }
    if (!(jgrid % 1000)) fprintf (pg,"finished grid point %ld\n",jgrid); fflush(pg);
  }
  fprintf (pg,"free memory\n");
  fclose(pg);
  
  free(work); free(H); free(tpsi);
  free(psi); free(phi);
  return;
}

/****************************************************************************/

double dotpreal(zomplex *psi,double *phi,long n,long m,long nGridPoints,double dv)
{
  long i;
  double sum;

  for (sum = 0.0, i = 0; i < nGridPoints; i++)
    sum += (psi[n*nGridPoints+i].re * phi[m*nGridPoints+i]);
  sum *= dv;
  return(sum);
}
