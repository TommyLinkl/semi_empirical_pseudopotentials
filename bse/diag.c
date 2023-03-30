/*****************************************************************************/

#include "fd.h"

/*****************************************************************************/
// Uses MKL/LAPACK dysev_ subroutine to diagonalize the real, symmetric matrix mat
// the eigenvalues are stored in eval and eigenvectors go into mat if jobz='V'

void diag(const long long n, int nthreads, double *mat, double *eval) {
  // n is order of the matrix mat (also leading dimension of mat)
  const long long lwork = 3*n; // length of the array work. lwork >= max(1,3*n-1)
  long long info; // signifies successful exit or not -> used to handle errors
  char jobz; 
  char uplo;
  double *work; // dimension (MAX(1,lwork))

  // Use multiple threads
  mkl_set_dynamic(0);
  mkl_set_num_threads(nthreads);
  omp_set_nested(1);

  // Allocate memmory for work array
  work = (double *) calloc(lwork, sizeof(double));

  // Set what will be computed along with if matrix is upper or lower traingular
  jobz = 'V'; // Compute eigenvalues and eigenvectors ('N' for eigenvalues only) 
  uplo = 'U'; // mat matrix is upper triangular

  // Diagonalize the real, symmetric matrix
  // TODO: uncomment and test the line below
  // dsyev_(&jobz, &uplo, &n, mat, &n, eval, work, &lwork, &info);
  dsyev_(&jobz, &uplo, &n, &mat[0], &n, &eval[0], &work[0], &lwork, &info);

  // Error handling
  if (info < 0) {
    printf("Argument %d of dysev_ in diag.c had an illegal argument\n", abs(info));
    printf("The program is exiting due to this error!!!\n");
    exit(EXIT_FAILURE);
  } else if (info > 0) {
    printf("dysev_ failed to converge -> the program is exiting!\n");
    exit(EXIT_FAILURE);
  }

  // Free dynamically allocated memory that was only internally used by MKL 
  free(work);

  return;
}

/*****************************************************************************/
