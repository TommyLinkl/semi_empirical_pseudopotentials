/*****************************************************************************************/
/* This file diagnolizes the hamiltonian matrix and stores the eigenvalues in energies */

/*****************************************************************************************/

#include "fit.h"

/*****************************************************************************************/
// Uses MKL/LAPACK zheev_ subroutine to diagonalize the complex, Hermitian matrix Hamiltonian
// the eigenvalues are stored in energies and eigenvectors go into Hamiltonian if jobz='V'

void diagonalizeHamiltonian(dcomplex *hamiltonian, double *energies, int nBasisVectors) {
  int i;
  const long long n = nBasisVectors; // order of the matrix A (n >= 0)
  const long long lda = n; // leading dimension of array the array A (lda >= max(1, N))
  const long long lwork = 3*n;
  long long info = 0;
  double *rwork;

  MKL_Complex16 *work;
  MKL_Complex16 *A;

  // Allocate memory 
  rwork  = (double *) calloc(3*n-2, sizeof(double));
  work  = (MKL_Complex16 *) calloc(lwork, sizeof(MKL_Complex16));
  A  = (MKL_Complex16 *) calloc(n*n, sizeof(MKL_Complex16));

  for (i = 0; i < n*n; i++) {
    A[i].real = hamiltonian[i].real;
    A[i].imag = hamiltonian[i].imag;
  }

  // zheev_ computes all eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix A.
  // Compute eigenvalues and eigenvectors 'V' or 'N' for eigenvalues only
  // hamiltonian matrix is upper triangular 'U' or lower triangular 'L'
  zheev_("N", "L", &n, A, &lda, energies, work, &lwork, rwork, &info);
  if (info) { // 0 gets returned on successful exit
    printf("Error in zheev_ in diag.c. Paramater %d is incorrect!\n", -info);
    printf("Program is exiting!!!\n");
    exit(EXIT_FAILURE);
  }

  // Free dynamically allocated memory
  free(rwork);
  free(work);
  free(A);

  return;
}

/*****************************************************************************************/
