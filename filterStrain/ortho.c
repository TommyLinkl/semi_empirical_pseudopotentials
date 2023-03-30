/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "fd.h"

/*****************************************************************************/

long portho(double *psi,double dv,long_st ist)
{
  long long lwork;  long long info, one=1, i, cutoff;
  long long ngrid = (long long)ist.ngrid, mstot = (long long)ist.mstot;
  double *S, *work;
  
  lwork = 5*(long long)(ist.mstot*ist.mstot+ist.ngrid);
  S = (double*) malloc(ist.mstot * sizeof(double));
  work = (double*) malloc(lwork * sizeof(double));

  if (ngrid*mstot<2147483647){
    mkl_set_dynamic(0);
    mkl_set_num_threads(ist.nthreads);
    omp_set_nested(1);
  }
  else {
    mkl_set_dynamic(0);
    mkl_set_num_threads(1);
    omp_set_nested(1);
  }
  {
    dgesvd_("O","N",&(ngrid),&(mstot),&(psi[0]),&(ngrid),&(S[0]),
	    NULL,&one,NULL,&one,&(work[0]),&(lwork),&info);
  }
  
  if (info != 0) {
    printf("error in dgesvd(2) %ld, exiting",info);
    exit(0);
  }
  
  for (cutoff = ist.mstot, i=0; i<ist.mstot; i++) {
    if ((S[i] / S[0]) < 1.0e-10) {
      cutoff = i;
      break;
    }
  }
  printf("cutoff is %ld\n",cutoff);
  
  free(work); 
  free(S);
  return (cutoff);
}

/*****************************************************************************/

