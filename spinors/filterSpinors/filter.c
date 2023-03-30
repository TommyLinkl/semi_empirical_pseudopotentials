#include "fd.h"

/*****************************************************************************/

void filtering(zomplex *psims,double *potl,double *ksqr,zomplex *an,double *zn,double *el,long_st ist,par_st par,nlc_st *nlc,long *nl,double *Elkb,long tid,long jns)
{
  FILE *pf; char str[100]; long flags = 0, ispn, jms, jgrid;
  fftw_plan_loc planfw, planbw; fftw_complex *fftwpsi;
  zomplex *psi, *phi;  double *eval;
  
  fftwpsi = fftw_malloc(sizeof (fftw_complex )*ist.ngrid);
  if ((psi = (zomplex*)calloc(ist.nspinngrid,sizeof(zomplex)))==NULL)nerror("psi");
  if ((phi = (zomplex*)calloc(ist.nspinngrid,sizeof(zomplex)))==NULL)nerror("phi");
  if ((eval = (double*)calloc(ist.ms,sizeof(double)))==NULL)nerror("eval");

  planfw = fftw_plan_dft_3d(ist.nz,ist.ny,ist.nx,fftwpsi,fftwpsi,FFTW_FORWARD,flags);
  planbw = fftw_plan_dft_3d(ist.nz,ist.ny,ist.nx,fftwpsi,fftwpsi,FFTW_BACKWARD,flags);

  for (jgrid = 0; jgrid < ist.nspinngrid; jgrid++) {
    psi[jgrid].re = psims[jgrid].re;
    psi[jgrid].im = psims[jgrid].im;
  }
  
  /***********************************************************************/
  /*** filter the states and normalize them ***/
  filter(psi,phi,psims,potl,ksqr,par,nlc,nl,Elkb,an,zn,ist,planfw,planbw,fftwpsi,tid,jns); 
  normalize_all(psims,par.dv,ist.ms,ist.nspinngrid,ist.nthreads);
  
  /*** calculate and print the energy of the filtered states ***/
  energy_all(psi,phi,psims,potl,ksqr,eval,ist,par,nlc,nl,Elkb,planfw,planbw,fftwpsi);

  sprintf (str,"eval-filt-%ld-%ld.dat",tid,jns);
  pf = fopen(str , "w");
  for (jms = 0; jms < ist.ms; jms++)
    fprintf (pf,"%ld %.16g %.16g\n",jms,eval[jms],el[jms]);
  fclose(pf);

  /*** write the filtered states to a file ***/
  sprintf (str,"psi-filt-%ld-%ld.dat",tid,jns);
  pf = fopen(str , "w");
  fwrite (psims,sizeof(double),ist.ms*ist.nspinngrid,pf);
  fclose(pf);
  
  free(psi); free(eval); free(phi);
  fftw_destroy_plan(planfw);
  fftw_destroy_plan(planbw);
  fftw_free(fftwpsi);
  return;
}

/*****************************************************************************/

void filter(zomplex *psin,zomplex *psim1,zomplex *psims,double *potl,double *ksqr,par_st par,nlc_st *nlc,long *nl,double *Elkb,zomplex *an,double *zn,long_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,long tid,long jns)
{
  FILE *pf; char str[100];
  long ispn, jgrid, jc, jms, ncjms, jmsg;

  for (jms = 0; jms < ist.ms; jms++){
    ncjms = ist.nc*jms; jmsg = jms * ist.nspinngrid;
    for (jgrid = 0; jgrid < ist.nspinngrid; jgrid++){
      psims[jmsg+jgrid].re = an[ncjms+0].re * psin[jgrid].re - an[ncjms+0].im * psin[jgrid].im;
      psims[jmsg+jgrid].im = an[ncjms+0].re * psin[jgrid].im + an[ncjms+0].im * psin[jgrid].re;
    }
  }

  sprintf (str,"prop%ld.dat",tid);
  for (pf = fopen(str , "w"), jc = 1; jc < ist.nc; jc++){
    memcpy(&psim1[0],&psin[0],ist.nspinngrid*sizeof(psim1[0]));
    hnorm(psim1,psin,potl,ksqr,par,nlc,nl,Elkb,zn[jc-1],ist,planfw,planbw,fftwpsi);
    
    for (jms = 0; jms < ist.ms; jms++){
      ncjms = ist.nc*jms; jmsg = jms * ist.nspinngrid;
      for (jgrid = 0; jgrid < ist.nspinngrid; jgrid++){
	psims[jmsg+jgrid].re += (an[ncjms+jc].re * psin[jgrid].re - an[ncjms+jc].im * psin[jgrid].im);
	psims[jmsg+jgrid].im += (an[ncjms+jc].re * psin[jgrid].im + an[ncjms+jc].im * psin[jgrid].re);
      }
    }
    if (!(jc % 100)) {fprintf (pf,"%ld %ld\n",jc,jns); fflush(pf);}
  }
  fclose(pf);
  return;
}

/*****************************************************************************/