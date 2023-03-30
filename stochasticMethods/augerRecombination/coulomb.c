/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "ar.h"

/*****************************************************************************/

void calcAllCoulombMatrixElements(double *vijck, double *vabck, double *psiai, double *evalai,
            double *psibs, zomplex *potq, double *pothhl, lng_st ist, par_st par, 
            fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwpsi)
{
  FILE *pf; zomplex *rho; 
  long jms, k, i, j, b, a, c, iGrid, tid;
  double sum;

  // Useful integers 
  const long THL = ist.nNonIntExcitons;
  const long TH2L = THL * ist.nHoles;
  const long THL2 = THL * ist.nElecs;

  // Write beginning of function
  writeSeparation(stdout);
  writeCurrentTime(stdout);
  fprintf(stdout, "Beginnning the calculation of all the Coulomb matrix elements\n\n");
  fprintf(stdout, "The number of Coulomb matrix elements that must be calculated = %ld\n", ist.itot*TH2L + ist.atot*THL2);
  fflush(stdout);

  // Dynamically allocate memory
  if ((rho = (zomplex *) calloc(ist.nGridPoints*ist.nThreads, sizeof(zomplex))) == NULL) nerror("rho");
  
  // Calculate the Coulomb matrix elements for the hot hole
#pragma omp parallel for private(i, j, iGrid, k, c, sum, tid)
  for (i = 0; i < ist.itot; i++) {
    tid = omp_get_thread_num();
    for (j = 0; j < ist.nHoles; j++) {
      for (iGrid = 0; iGrid < ist.nGridPoints; iGrid++) {
        rho[tid*ist.nGridPoints+iGrid].re = psiai[i*ist.nGridPoints+iGrid] * psibs[j*ist.nGridPoints+iGrid];      	
        rho[tid*ist.nGridPoints+iGrid].im = 0.0;
      }
      hartree(&rho[tid*ist.nGridPoints],potq,&pothhl[tid*ist.nGridPoints],ist.nGridPoints,planfw[tid],planbw[tid],&fftwpsi[tid*ist.nGridPoints]);
      
      for (c = 0; c < ist.nElecs; c++) {
        for (k = 0; k < ist.nHoles; k++) {
      	  sum = 0.0;   
      	  for (iGrid = 0; iGrid < ist.nGridPoints; iGrid++) {
      	    sum += psibs[k*ist.nGridPoints+iGrid] * psibs[(c+ist.nHoles)*ist.nGridPoints+iGrid] * pothhl[tid*ist.nGridPoints+iGrid];
          }
      	  vijck[i*TH2L + j*THL + c*ist.nHoles + k] = sum*par.dv;
      	}
      }
    }
  }
  // Print the Coulomb matrix elements for the hot hole
  pf = fopen("vijck.dat" , "w");
  for (i = 0; i < ist.itot; i++) {
    for (j = 0; j < ist.nHoles; j++) {
      for (c = 0; c < ist.nElecs; c++) {
        for (k = 0; k < ist.nHoles; k++) {
      	  fprintf(pf, "%ld %ld %ld %ld %g %g\n", i, j, c, k, vijck[i*TH2L + j*THL + c*ist.nHoles + k],
            vijck[i*TH2L + j*THL + c*ist.nHoles + k]*vijck[i*TH2L + j*THL + c*ist.nHoles + k]);
        }
      }
    }
  }
  fclose(pf);

  // Calculate the Coulomb matrix elements for the hot elec
#pragma omp parallel for private(a, b, iGrid, c, k, sum, tid)
  for (a = 0; a < ist.atot; a++) {
    tid = omp_get_thread_num();
    for (b = 0; b < ist.nElecs; b++) {
      for (iGrid = 0; iGrid < ist.nGridPoints; iGrid++) {
      	rho[tid*ist.nGridPoints+iGrid].re = psiai[(a+ist.itot)*ist.nGridPoints+iGrid] * psibs[(b+ist.nHoles)*ist.nGridPoints+iGrid];
        rho[tid*ist.nGridPoints+iGrid].im = 0.0;
      }
      hartree(&rho[tid*ist.nGridPoints],potq,&pothhl[tid*ist.nGridPoints],ist.nGridPoints,planfw[tid],planbw[tid],&fftwpsi[tid*ist.nGridPoints]);
      
      for (c = 0; c < ist.nElecs; c++) {
        for (k = 0; k < ist.nHoles; k++) {
      	  sum = 0.0;
      	  for (iGrid = 0; iGrid < ist.nGridPoints; iGrid++) {
      	    sum += psibs[k*ist.nGridPoints+iGrid] * psibs[(c+ist.nHoles)*ist.nGridPoints+iGrid] * pothhl[tid*ist.nGridPoints+iGrid];
          }
      	  vabck[a*THL2 + b*THL + c*ist.nHoles + k] = sum*par.dv;
      	}
      }
    }
  }
  // Print the Coulomb matrix elements for the hot elec
  pf = fopen("vabck.dat" , "w");
  for (a = 0; a < ist.atot; a++) {
    for (b = 0; b < ist.nElecs; b++) {
      for (c = 0; c < ist.nElecs; c++) {
        for (k = 0; k < ist.nHoles; k++) {
      	  fprintf(pf, "%ld %ld %ld %ld %g %g\n", a, b, c, k, vabck[a*THL2 + b*THL + c*ist.nHoles + k],
            vabck[a*THL2 + b*THL + c*ist.nHoles + k]*vabck[a*THL2 + b*THL + c*ist.nHoles + k]);
        }
      }
    }
  }
  fclose(pf);
  
  // Write ending of function
  fprintf(stdout, "\nFinished the calculation of the Coulomb matrix elements\n");
  writeCurrentTime(stdout);
  fflush(stdout);

  // Free dynamically allocated memory
  free(rho);
  
  return;
}

/*****************************************************************************/
// Calculate and return a single Coulomb matrix element: Vrsut = <rs|V|ut> 
// Args: psiX -> pointers to the first element of each state
//       potq -> q-space potential of the Coulomb potential
//       dV -> volume element, 
//       nGrid -> number of grid points

double calcOneCoulombMatrixElement(double *psiR, double *psiS, double *psiU, double *psiT, zomplex *potq,
            double dV, long nGrid, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi) {
  long iGrid;
  double *hatreePot, Vrsut, sum = 0.0;
  zomplex *rho;

  // Dynamically allocate memory
  if ((hatreePot = (double *) calloc(nGrid, sizeof(double))) == NULL) nerror("hatreePot");
  if ((rho = (zomplex *) calloc(nGrid, sizeof(zomplex))) == NULL) nerror("rho");

  // Calculate the psiT*psiU then call hartree to perform forward and backward FTs
#pragma omp parallel for private(iGrid)
  for (iGrid = 0; iGrid < nGrid; iGrid++) {
    rho[iGrid].re = psiU[iGrid] * psiT[iGrid]; // product nStatesPerFilter of the states
    rho[iGrid].im = 0.0;
  }
  hartree(rho, potq, hatreePot, nGrid, planfw, planbw, &fftwpsi[0]);

  // Use the result of hartree, hatreePot, to calculate the matrix element 
#pragma omp parallel for private(iGrid)
  for (iGrid = 0; iGrid < nGrid; iGrid++) {
  #pragma omp critical
    sum += psiR[iGrid] * psiS[iGrid] * hatreePot[iGrid];
  }
  Vrsut = sum * dV;

  // Testing
  //printf("Vrsut = %g\n", Vrsut);

  // Free dynamically allocate memory
  free(rho); free(hatreePot);

  return (Vrsut);
}

/*****************************************************************************/
