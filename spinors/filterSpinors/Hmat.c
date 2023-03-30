#include "fd.h"

/*****************************************************************************/

void Hmat(zomplex *psi,zomplex *phi,MKL_Complex16 *psitot,double *potl,double *ksqr,double *eval,long_st ist,par_st par,nlc_st *nlc,long *nl,double *Elkb,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  FILE *pg;
  long ims, jms, jgrid;
  long long mstot = (long long)ist.mstot, info, lwk = 3*ist.mstot;
  double *rwork, sumre, sumim; MKL_Complex16 *H, *work; zomplex *tpsi;
  
  H = (MKL_Complex16*)calloc(ist.mstot*ist.mstot,sizeof(MKL_Complex16));
  work = (MKL_Complex16*)calloc(lwk,sizeof(MKL_Complex16));
  rwork = (double*)calloc(3*ist.mstot,sizeof(double));
  tpsi = (zomplex*)calloc(ist.mstot,sizeof(zomplex));

  /*omp_set_dynamic(0);
    omp_set_num_threads(ist.nthreads);*/

  /*** calculate H|psi_i> ***/
  pg = fopen("hmat.dat" , "w");
  fprintf (pg,"calculate the H matrix\n"); fflush(pg);
  for (ims = 0; ims < ist.mstot; ims++){
    //#pragma omp parallel for private(jgrid)
    for (jgrid = 0; jgrid < ist.nspinngrid; jgrid++) {
      psi[jgrid].re = psitot[ims*ist.nspinngrid+jgrid].real;
      psi[jgrid].im = psitot[ims*ist.nspinngrid+jgrid].imag;
    }
    memcpy(&phi[0],&psi[0],ist.nspinngrid*sizeof(phi[0]));
    hamiltonian(phi,psi,potl,ksqr,ist,par,nlc,nl,Elkb,planfw,planbw,fftwpsi);

    /*** calculate <psi_j|H|psi_i> ***/
    //#pragma omp parallel for private(jms)
    for (jms = 0; jms < ist.mstot; jms++){
      H[ims*ist.mstot+jms] = dotp(phi,psitot,jms,ist.nspinngrid,par.dv);
      //fprintf (pg,"%ld %ld %g %g\n",ims,jms,H[ims*ist.mstot+jms].real,H[ims*ist.mstot+jms].imag);
    }
    fprintf (pg,"finshed row %ld\n",ims); fflush(pg);
  }

  /*** diagonalize the Hamiltonian H ***/
  fprintf (pg,"diagonalize the new Hamiltonian H\n"); fflush(pg);
  zheev_("V","U",&mstot,&H[0],&mstot,&eval[0],&work[0],&lwk,&(rwork[0]),&info);
  if (info) nerror("error in zheev H");
    
  /*** copy the new function into psitot ***/
  fprintf (pg,"copy the new function into psitot\n"); fflush(pg);

  for (jgrid = 0; jgrid < ist.nspinngrid; jgrid++){
#pragma omp parallel for private (jms)
    for (jms = 0; jms < ist.mstot; jms++){
      tpsi[jms].re = psitot[jms*ist.nspinngrid+jgrid].real;
      tpsi[jms].im = psitot[jms*ist.nspinngrid+jgrid].imag;
      psitot[jms*ist.nspinngrid+jgrid].real = psitot[jms*ist.nspinngrid+jgrid].imag = 0.0;
    }

#pragma omp parallel for private (jms,sumre,sumim)
    for (jms = 0; jms < ist.mstot; jms++) {      
      for (sumre = sumim = 0.0, ims = 0; ims < ist.mstot; ims++){
	sumre += (H[jms*ist.mstot+ims].real * tpsi[ims].re + H[jms*ist.mstot+ims].imag * tpsi[ims].im);
	sumim += (H[jms*ist.mstot+ims].real * tpsi[ims].im - H[jms*ist.mstot+ims].imag * tpsi[ims].re);
      }

      //#pragma omp critical
      psitot[jms*ist.nspinngrid+jgrid].real = sumre;
      psitot[jms*ist.nspinngrid+jgrid].imag = sumim;
    }
    if (!(jgrid % 1000)) fprintf (pg,"finished grid point %ld\n",jgrid); fflush(pg);
  }
  fprintf (pg,"free memory\n");
  fclose(pg);

  free(work); free(H); free(tpsi); free(rwork);
  return;
}

/****************************************************************************/

MKL_Complex16 dotp(zomplex *psi,MKL_Complex16 *phi,long m,long ngrid,double dv)
{
  long i;
  MKL_Complex16 sum;

  for (sum.real = sum.imag = 0.0, i = 0; i < ngrid; i++){
    sum.real += (psi[i].re * phi[m*ngrid+i].real + psi[i].im * phi[m*ngrid+i].imag);
    sum.imag += (psi[i].re * phi[m*ngrid+i].imag - psi[i].im * phi[m*ngrid+i].real);
  }
  sum.real *= dv;
  sum.imag *= dv;
  return(sum);
}

/*****************************************************************************/