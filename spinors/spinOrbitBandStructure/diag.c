/****************************************************************************/
/* This file diagnolizes the hamiltonian matrix and stores the eigenvalues in energies */

#include "fit.h"

/****************************************************************************/

void diagnolizeHamiltonian(complexnumber *hamiltonian, double *energies, param params) {
  int n = 2*params.nBasisVectors;
  int lwork = 2 * n -1;
  int info;
  char  jobz, uplo;
  double *rwork;
  complexnumber *work;
 
  rwork  = (double *) calloc(3*n-2, sizeof(double));
  work  = (complexnumber *) calloc(lwork, sizeof(complexnumber));
 
  /* Compute only eigenvalues */
  jobz = 'N';
  uplo = 'U';
  info = 0;
 
  zheev_(&jobz, &uplo, &n, &hamiltonian[0], &n,
      &energies[0], &work[0], &lwork, &rwork[0], &info);
 
  free(rwork);
  free(work);

  return;
}

/****************************************************************************/
