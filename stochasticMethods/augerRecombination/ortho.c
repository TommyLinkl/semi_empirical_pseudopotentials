#include "ar.h"

/*****************************************************************************/

long portho(double *psi,double dv,lng_st ist)
{
  long long lwork;  long long info, one=1, i, cutoff;
  long long nGridPoints = (long long)ist.nGridPoints, ms = (long long)ist.nFilteredStates;
  double *S, *work;
  
  // Dynamically allocate memory
  lwork = 5*(long long)(ist.nFilteredStates*ist.nFilteredStates+ist.nGridPoints);
  S = (double*) malloc(ist.nFilteredStates * sizeof(double));
  work = (double*) malloc(lwork * sizeof(double));

  if (nGridPoints*ms<2147483647){
    mkl_set_dynamic(0);
    mkl_set_num_threads(ist.nThreads);
    omp_set_nested(1);
  }
  else {
    mkl_set_dynamic(0);
    mkl_set_num_threads(1);
    omp_set_nested(1);
  }
  {
    dgesvd_("O","N",&(nGridPoints),&(ms),&(psi[0]),&(nGridPoints),&(S[0]),
	    NULL,&one,NULL,&one,&(work[0]),&(lwork),&info);
  }
  
  if (info != 0) {
    printf("error in dgesvd(2) %ld, exiting",info);
    exit(0);
  }
  
  for (cutoff = ist.nFilteredStates, i=0; i<ist.nFilteredStates; i++) {
    if ((S[i] / S[0]) < 1.0e-10) {
      cutoff = i;
      break;
    }
  }

  // Free dynamically allocated memory  
  free(work); 
  free(S);
  
  return (cutoff);
}

/*****************************************************************************/
