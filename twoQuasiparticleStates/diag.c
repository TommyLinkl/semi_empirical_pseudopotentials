/*****************************************************************************/
//
// This file deals with diagonalizing matrices
//
/*****************************************************************************/

#include "qp.h"

/*****************************************************************************/
// Wraps diag and stores / prints out some interacting two quasiparticle state information

long diagonalizeTwoQPHamiltonain(intTwoQPState *intTwoQP, double *hMatrix, 
                                 nonintTwoQPState *nonintTwoQP, lParams lPar) {
  FILE *pf;
  long iIntTwoQP, iNonintTwoQP;
  double *eigenvalues;

  // Write beginning of function
  writeSeparation(stdout);
  writeCurrentTime(stdout);
  fprintf(stdout, "Beginnning the diagonalization of the two quasiparticle Hamiltonian\n\n");
  fflush(stdout);

  // Allocate memory
  if ((eigenvalues = (double *) calloc(lPar.nNonintTwoQPStates, sizeof(double))) == NULL) memoryError("eigenvalues");

  // Set constants needed for diag function call
  const long long orderOfHamiltonian = lPar.nNonintTwoQPStates;
  int nThreads = lPar.nThreads;

  // Diagonalize the two quasiparticle hamiltonian matrix
  diag(orderOfHamiltonian, nThreads, hMatrix, eigenvalues);

  // Store results of the diagonalization
  for (iIntTwoQP = 0; iIntTwoQP < lPar.nIntTwoQPStates; iIntTwoQP++) {
    intTwoQP[iIntTwoQP].index = iIntTwoQP;
    intTwoQP[iIntTwoQP].nNonintTwoQPStates = lPar.nNonintTwoQPStates;
    intTwoQP[iIntTwoQP].energy = eigenvalues[iIntTwoQP];
    intTwoQP[iIntTwoQP].Crs = &(hMatrix[iIntTwoQP*lPar.nNonintTwoQPStates]);
    intTwoQP[iIntTwoQP].niTwoQP = nonintTwoQP;
    intTwoQP[iIntTwoQP].bindingEnergy = nonintTwoQP[0].energy-intTwoQP[iIntTwoQP].energy;
    intTwoQP[iIntTwoQP].correlationEnergy = (nonintTwoQP[0].energy + nonintTwoQP[0].hartreeEnergy +
                                            nonintTwoQP[0].exchangeEnergy - intTwoQP[iIntTwoQP].energy);
  }    

  // Print out the energies and coefficients of the interacting two quasiparticle eigenstates
  writeIntTwoQPStructuresToFiles(intTwoQP, lPar.nNonintTwoQPStates, lPar.nIntTwoQPStates);

  // Write ending of function
  fprintf(stdout, "\nFinished the diagonalization of the two quasiparticle Hamiltonian\n");
  writeCurrentTime(stdout);
  fflush(stdout);

  // Free dynamically allocated memory
  free(eigenvalues);

  return 0;
}


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
  dsyev_(&jobz, &uplo, &n, &mat[0], &n, &eval[0], &work[0], &lwork, &info);

  // Error handling
  if (info < 0) {
    printf("Argument %d of dysev_ in diag.c had an illegal argument\n", abs((int)info));
    printf("The program is exiting due to this error!!!\n");
    exit(EXIT_FAILURE);
  } 
  else if (info > 0) {
    printf("dysev_ failed to converge -> the program is exiting!\n");
    exit(EXIT_FAILURE);
  }

  // Free dynamically allocated memory that was only internally used by MKL 
  free(work);

  return;
}

/*****************************************************************************/
