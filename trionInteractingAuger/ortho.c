#include "fd.h"

/*****************************************************************************/

long portho(double *psi,double dv,lng_st ist)
{
  long long lwork;  long long info=0, one=1, i, cutoff;
  long long ngrid = (long long)ist.ngrid, ms = (long long)ist.ms;
  double *S, *work;
  
  lwork = 5*(long long)(ist.ms*ist.ms+ist.ngrid);
  S = (double*) malloc(ist.ms * sizeof(double));
  work = (double*) malloc(lwork * sizeof(double));

  if (ngrid*ms<2147483647){
    mkl_set_dynamic(0);
    mkl_set_num_threads(ist.nthreads);
    omp_set_max_active_levels(1);
    //omp_set_nested(1);
  }
  else {
    mkl_set_dynamic(0);
    mkl_set_num_threads(1);
    omp_set_max_active_levels(1);
    //omp_set_nested(1);
  }
  {
    dgesvd_("O","N",&(ngrid),&(ms),&(psi[0]),&(ngrid),&(S[0]),
	    NULL,&one,NULL,&one,&(work[0]),&(lwork),&info);
  }
  
  if (info != 0) {
    printf("error in dgesvd(2) %lld, exiting",info);
    exit(0);
  }
  
  for (cutoff = ist.ms, i=0; i<ist.ms; i++) {
    if ((S[i] / S[0]) < 1.0e-10) {
      cutoff = i;
      break;
    }
  }
  printf("cutoff is %lld\n",cutoff);
  
  free(work); 
  free(S);
  return (cutoff);
}

/*****************************************************************************/
