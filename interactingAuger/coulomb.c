#include "fd.h"

/*****************************************************************************/

void generate_coulomb_matrix_elements(double *vkijb,double *vkcab,double *psiai,double *evalai,double *sigeai,double *psibs,zomplex *potq,double *pothhl,lng_st ist,par_st par,fftw_plan_loc *planfw,fftw_plan_loc *planbw,fftw_complex *fftwpsi)
{
  FILE *pf; zomplex *rho; double sum;
  long jms, k, i, j, b, a, c, igrid, tid, *list, thl, th2l;

  // Write time this function began
  writeFunctionStartTime("Program is now calculating the Coulomb matrix elements:", stdout);
  
  // Allocate memory used within this function
  if ((rho = (zomplex *) calloc(ist.ngrid*ist.nthreads, sizeof(zomplex))) == NULL) nerror("rho");

  thl = ist.numBandEdgeHoles * ist.numBandEdgeElectrons;
  th2l = thl * ist.numBandEdgeHoles;
  if ((list = (long *) calloc(ist.itot, sizeof(long))) == NULL) nerror("list i");

  for (i = 0, jms = 0; jms < ist.ms; jms++) {
    if ((evalai[jms] >= par.Eimin) && (evalai[jms] <= par.Eimax) && (sigeai[jms] < par.deps)) {
      list[i] = jms; 
      i++;
    }
  }
  
#pragma omp parallel for private(i,k,igrid,j,b,sum,tid)
  for (i = 0; i < ist.itot; i++) {
    tid = omp_get_thread_num();
    for (k = 0; k < ist.numBandEdgeHoles; k++) {
      for (igrid = 0; igrid < ist.ngrid; igrid++) {
      	rho[tid*ist.ngrid+igrid].re = psiai[list[i]*ist.ngrid+igrid] * psibs[k*ist.ngrid+igrid];
      	rho[tid*ist.ngrid+igrid].im = 0.0;
      }
      hartree(&rho[tid*ist.ngrid],potq,&pothhl[tid*ist.ngrid],ist,planfw[tid],planbw[tid],&fftwpsi[tid*ist.ngrid]);
      
      for (j = 0; j < ist.numBandEdgeHoles; j++) {
      	for (b = 0; b < ist.numBandEdgeElectrons; b++) {
      	  for (sum = 0.0, igrid = 0; igrid < ist.ngrid; igrid++) {
      	    sum += psibs[j*ist.ngrid+igrid] * psibs[(b+ist.numBandEdgeHoles)*ist.ngrid+igrid] * pothhl[tid*ist.ngrid+igrid];
          }
      	  vkijb[i*th2l + j*thl + k*ist.numBandEdgeElectrons + b] = sum * par.dv;
      	}
      }
    }
  }
  free(list);
  
  pf = fopen("vkijb.dat" , "w");
  for (i = 0; i < ist.itot; i++)
    for (k = 0; k < ist.numBandEdgeHoles; k++)
      for (j = 0; j < ist.numBandEdgeHoles; j++)
      	for (b = 0; b < ist.numBandEdgeElectrons; b++)
      	  fprintf (pf,"%ld %ld %ld %ld %g\n",i,k,j,b,vkijb[i*th2l + j*thl + k*ist.numBandEdgeElectrons + b]);
  fclose(pf);

  thl = ist.numBandEdgeHoles * ist.numBandEdgeElectrons;
  th2l = thl * ist.numBandEdgeElectrons;
  if ((list = (long *) calloc(ist.atot, sizeof(long))) == NULL) nerror("list a");

  for (a = 0, jms = 0; jms < ist.ms; jms++) {
    if ((evalai[jms] >= par.Eamin) && (evalai[jms] <= par.Eamax) && (sigeai[jms] < par.deps)) {
      list[a] = jms;
      a++;
    }
  }

#pragma omp parallel for private(a,c,igrid,j,b,sum,tid)
  for (a = 0; a < ist.atot; a++) {
    tid = omp_get_thread_num();
    for (b = 0; b < ist.numBandEdgeElectrons; b++) {
      for (igrid = 0; igrid < ist.ngrid; igrid++) {
      	rho[tid*ist.ngrid+igrid].re = psiai[list[a]*ist.ngrid+igrid] * psibs[(b+ist.numBandEdgeHoles)*ist.ngrid+igrid];
      	rho[tid*ist.ngrid+igrid].im = 0.0;
      }
      hartree(&rho[tid*ist.ngrid],potq,&pothhl[tid*ist.ngrid],ist,planfw[tid],planbw[tid],&fftwpsi[tid*ist.ngrid]);
      
      for (k = 0; k < ist.numBandEdgeHoles; k++) {
      	for (c = 0; c < ist.numBandEdgeElectrons; c++) {
      	  for (sum = 0.0, igrid = 0; igrid < ist.ngrid; igrid++) {
      	    sum += psibs[k*ist.ngrid+igrid] * psibs[(c+ist.numBandEdgeHoles)*ist.ngrid+igrid] * pothhl[tid*ist.ngrid+igrid];
          }
      	  vkcab[a*th2l + b*thl + k*ist.numBandEdgeElectrons + c] = sum * par.dv;
      	}
      }
    }
  }
  free(list);
  
  pf = fopen("vkcab.dat" , "w");
  for (a = 0; a < ist.atot; a++) 
    for (b = 0; b < ist.numBandEdgeElectrons; b++)
      for (k = 0; k < ist.numBandEdgeHoles; k++)
      	for (c = 0; c < ist.numBandEdgeElectrons; c++)
      	  fprintf (pf,"%ld %ld %ld %ld %g\n",a,b,k,c,vkcab[a*th2l + b*thl + k*ist.numBandEdgeElectrons + c]);
  fclose(pf);
  
  free(rho);
  
  return;
}

/*****************************************************************************/
