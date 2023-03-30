#include "ar.h"

/*****************************************************************************/

void filtering(double *psitot,double *potl,double *ksqr,zomplex *an,double *zn,double *el,lng_st ist,par_st par,long tid,long jns,long *idum,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  FILE *pf; char str[100]; long flags = 0, jms, jgrid;
  zomplex *psi, *phi;  double *eval;
  
  if ((psi = (zomplex*)calloc(ist.nGridPoints,sizeof(zomplex)))==NULL)nerror("psi");
  if ((phi = (zomplex*)calloc(ist.nGridPoints,sizeof(zomplex)))==NULL)nerror("psi");
  if ((eval = (double*)calloc(ist.nStatesPerFilter,sizeof(double)))==NULL)nerror("eval");

  /*** start from an initial random state ***/  
  //init_psi(psi,ist,par,&idum[tid]);
  
  for (jgrid = 0; jgrid < ist.nGridPoints; jgrid++) psi[jgrid].re = psitot[jgrid];
  
  /***********************************************************************/
  /*** filter the states and normalize them ***/
  filter(psi,phi,psitot,potl,ksqr,par,an,zn,ist,planfw,planbw,fftwpsi,tid,jns); 
  normalize_all(psitot,par.dv,ist.nStatesPerFilter,ist.nGridPoints);
  
  /*** calculate and print the energy of the filtered states ***/
  energy_all(psi,phi,psitot,potl,ksqr,eval,ist,par,planfw,planbw,fftwpsi,ist.nStatesPerFilter);

  sprintf (str,"eval-filt-%ld-%ld.dat",tid,jns);
  pf = fopen(str , "w");
  for (jms = 0; jms < ist.nStatesPerFilter; jms++)
    fprintf (pf,"%ld %.16g %.16g\n",jms,eval[jms],el[jms]);
  fclose(pf);

  /*** store the filtered states ***/
  sprintf (str,"psi-filt-%ld-%ld.dat",tid,jns);
  pf = fopen(str , "w");
  fwrite (psitot,sizeof(double),ist.nStatesPerFilter*ist.nGridPoints,pf);
  fclose(pf);
  
  free(psi); free(eval); free(phi);
  return;
}

/*****************************************************************************/

void filter(zomplex *psin,zomplex *psim1,double *psi0,double *potl,double *ksqr,par_st par,zomplex *an,double *zn,lng_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,long tid,long jns)
{
  FILE *pf; char str[100];
  long i, j, ie, ncie, ieg;

  for (ie = 0; ie < ist.nStatesPerFilter; ie++){
    ncie = ist.nNewtonIntSteps*ie; ieg = ie * ist.nGridPoints;
    for (i = 0; i < ist.nGridPoints; i++) psi0[ieg+i] = (an[ncie+0].re*psin[i].re);
  }

  sprintf (str,"prop%ld.dat",tid);
  for (pf = fopen(str , "w"), j = 1; j < ist.nNewtonIntSteps; j++){
    memcpy(&psim1[0],&psin[0],ist.nGridPoints*sizeof(psim1[0]));
    hnorm(psim1,psin,potl,ksqr,par,zn[j-1],ist,planfw,planbw,fftwpsi);
    
    for (ie = 0; ie < ist.nStatesPerFilter; ie++){
      ncie = ist.nNewtonIntSteps*ie; ieg = ie * ist.nGridPoints;
      for (i = 0; i < ist.nGridPoints; i++)
	psi0[ieg+i] += (an[ncie+j].re*psin[i].re);
    }
    if (!(j % 100)) {fprintf (pf,"%ld %ld\n",j,jns); fflush(pf);}
  }
  fclose(pf);
  return;
}

/*****************************************************************************/
