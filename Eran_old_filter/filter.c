#include "fd.h"

void filtering(double *psitot,double *potl,double *ksqr,zomplex *an,zomplex *zn,double *el,long_st ist,par_st par,long tid,long jns,long *idum)
{
  FILE *pf; char str[100]; long flags = 0, jms, jgrid;
  fftw_plan_loc planfw, planbw; fftw_complex *fftwpsi;
  zomplex *psi, *phi;  double *eval;
  /*
  int err = fftw_init_threads();
  if(err==0){
    printf("Error when initializing multi-thread FFTW3D!\n");
    fflush(stdout);
    exit(0);
  }
  */
  
  fftwpsi = fftw_malloc(sizeof (fftw_complex )*ist.ngrid);
  if ((psi = (zomplex*)calloc(ist.ngrid,sizeof(zomplex)))==NULL)nerror("psi");
  if ((phi = (zomplex*)calloc(ist.ngrid,sizeof(zomplex)))==NULL)nerror("psi");
  if ((eval = (double*)calloc(ist.ms,sizeof(double)))==NULL)nerror("eval");

  planfw = fftw_plan_dft_3d(ist.nz,ist.ny,ist.nx,fftwpsi,fftwpsi,FFTW_FORWARD,flags);
  planbw = fftw_plan_dft_3d(ist.nz,ist.ny,ist.nx,fftwpsi,fftwpsi,FFTW_BACKWARD,flags);
  
  /*** start from an initial random state ***/  
  //init_psi(psi,ist,par,&idum[tid]);
  
  for (jgrid = 0; jgrid < ist.ngrid; jgrid++) psi[jgrid].re = psitot[jgrid];
  
  /***********************************************************************/
  /*** filter the states and normalize them ***/
  filter(psi,phi,psitot,potl,ksqr,par,an,zn,ist,planfw,planbw,fftwpsi,tid,jns); 
  normalize_all(psitot,par.dv,ist.ms,ist.ngrid,ist.nthreads);
  
  /*** calculate and print the energy of the filtered states ***/
  energy_all(psi,phi,psitot,potl,ksqr,eval,ist,par,planfw,planbw,fftwpsi);

  sprintf (str,"eval-filt-%ld-%ld.dat",tid,jns);
  pf = fopen(str , "w");
  for (jms = 0; jms < ist.ms; jms++)
    fprintf (pf,"%ld %.16g %.16g\n",jms,eval[jms],el[jms]);
  fclose(pf);

  /*** write the filtered states to a file ***/
  sprintf (str,"psi-filt-%ld-%ld.dat",tid,jns);
  pf = fopen(str , "w");
  fwrite (psitot,sizeof(double),ist.ms*ist.ngrid,pf);
  fclose(pf);
  
  free(psi); free(eval); free(phi);
  fftw_destroy_plan(planfw);
  fftw_destroy_plan(planbw);
  fftw_free(fftwpsi);
  return;
}

/****************************************************************************************/

void filter(zomplex *psin,zomplex *psim1,double *psitot,double *potl,double *ksqr,par_st par,zomplex *an,zomplex *zn,long_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,long tid,long jns)
{
  FILE *pf; char str[100];
  long jcheby, jgrid, jtarget;
  long ntarget = ist.ms, ngrid = ist.ngrid, ncheby = ist.nc;

  for (jtarget = 0; jtarget < ntarget; jtarget++)
    for (jgrid = 0; jgrid < ngrid; jgrid++)
      psitot[jtarget*ngrid+jgrid] = (an[ncheby*jtarget+0].re*psin[jgrid].re);

  sprintf (str,"prop%ld.dat",tid);
  for (pf = fopen(str , "w"), jcheby = 1; jcheby < ncheby; jcheby++){
    memcpy(&psim1[0],&psin[0],ngrid*sizeof(psim1[0]));
    hnorm(psim1,psin,potl,ksqr,par,zn[jcheby-1].re,ist,planfw,planbw,fftwpsi);
    
    for (jtarget = 0; jtarget < ntarget; jtarget++)
      for (jgrid = 0; jgrid < ist.ngrid; jgrid++)
	psitot[jtarget*ngrid+jgrid] += (an[ncheby*jtarget+jcheby].re*psin[jgrid].re);
    
    if (!(jcheby % 100)) {fprintf (pf,"%ld %ld\n",jcheby,jns); fflush(pf);}
  }
  fclose(pf);
  return;
}
