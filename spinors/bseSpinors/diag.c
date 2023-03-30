#include "fd.h"

/****************************************************************************/
// Diagnolizes the double matrix mat by calling dysev_ (LAPACK) and stores  
// the eigenvectors (coefficients in BSE) in mat and the eigenvalues in eval

void diag(int n, int nthreads, double *mat, double *eval)
{
  int lwork = 3*n, info;
  char jobz, uplo;
  double *work;

  mkl_set_dynamic(0);
  mkl_set_num_threads(nthreads);
  omp_set_nested(1);

  work=(double*)calloc(lwork,sizeof(double));
  jobz = 'V'; uplo = 'U';
  dsyev_(&jobz,&uplo,&n,&mat[0],&n,&eval[0],&work[0],&lwork,&info);
  if (info) exit(0);

  free(work);

  return;
}

/****************************************************************************/