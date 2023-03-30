#include "fd.h"

/***************************************************************************************/
// Calculates bsmat = <Sai|Keh|Sbj> for excitonic states Sai and Sbj where Keh is the 
// electon-hole kernel and contains screened direct coulomb terms and unscreened 
// direct exchange terms and calculates h0mat as the diagnol matrix elements (Ea-Ei) 


// Function needs to be majorally edited to work with complex spinors
void single_coulomb_openmp(zomplex *psi,zomplex *potq,zomplex *potqx,double *poth,double *eval,lng_st ist,par_st par,fftw_plan_loc *planfw,fftw_plan_loc *planbw,fftw_complex *fftwpsi,double *bsmat,double *h0mat)
{
  FILE *pf;  
  long flag, i, j, a, b, ibs, jbs, igrid, *listibs; 
  int tid; zomplex *rho;
  double ene, ene1, ene2, sum1, sum2;
	  
  rho = (zomplex*)calloc(ist.ngrid*ist.nthreads,sizeof(zomplex));
  listibs = (long*)calloc(ist.ms2,sizeof(long));

  for (ibs = 0, a = ist.nlumo; a < ist.nlumo+ist.totallumo; a++) {
    for (i = 0; i < ist.totalhomo; i++, ibs++) {
      listibs[(a-ist.nlumo)*ist.totalhomo+i] = ibs; 
    }
  }

  omp_set_dynamic(0);
  omp_set_num_threads(ist.nthreads);

  pf = fopen("single-coulomb.dat" , "w");
  /*** vabji direct ***/
  for (a = ist.nlumo; a < ist.nlumo+ist.totallumo; a++) {
#pragma omp parallel for private(sum1,ibs,jbs,ene1,ene2,ene,tid,igrid,b,i,j)
    for (b = ist.nlumo; b < ist.nlumo+ist.totallumo; b++) {
      tid = omp_get_thread_num();
      for (igrid = 0; igrid < ist.ngrid; igrid++) {
        // TODO: just put .re here for now - needs to be corrected
        // TODO: make sure complex conjugates are correct as wellll
        rho[tid*ist.ngrid+igrid].re = (psi[a*ist.ngrid+igrid].re*psi[b*ist.ngrid+igrid].re);
        rho[tid*ist.ngrid+igrid].im = 0.0;
      }
      hartree(&rho[tid*ist.ngrid], potqx, &poth[tid*ist.ngrid], ist, planfw[tid], planbw[tid], &fftwpsi[tid*ist.ngrid]);
      for (i = 0; i < ist.totalhomo; i++) {
        for (j = 0; j < ist.totalhomo; j++) {
      	  ene1 = eval[a] - eval[i];
      	  ene2 = eval[b] - eval[j];
      	  ene = ene1 - ene2;

          for (sum1 = 0.0, igrid = 0; igrid < ist.ngrid; igrid++) {
            // TODO: just put .re here for now - needs to be corrected
            // TODO: make sure complex conjugates are correct as well
            sum1 += (poth[tid*ist.ngrid+igrid] * psi[j*ist.ngrid+igrid].re * psi[i*ist.ngrid+igrid].re);
          }
      	  sum1 *= par.dv;
      	  ibs = listibs[(a-ist.nlumo)*ist.totalhomo+i];
      	  jbs = listibs[(b-ist.nlumo)*ist.totalhomo+j];
      	  
      	  //bsmat[ibs*ist.ms2+jbs] = sum1;
      	  if (ibs == jbs) h0mat[ibs*ist.ms2+jbs] = eval[a] - eval[i];
      	  else h0mat[ibs*ist.ms2+jbs] = 0.0;
      	  fprintf (pf,"%ld %ld %ld %ld %ld %ld %g %g %g\n",a,i,b,j,
      		   listibs[(a-ist.nlumo)*ist.totalhomo+i],
      		   listibs[(b-ist.nlumo)*ist.totalhomo+j],
      		   ene1,ene2,sum1);
      	  fflush(0);
        }
      }
    }
  }

  
  /*** vjbai exchange ***/
  for (a = ist.nlumo; a < ist.nlumo+ist.totallumo; a++) {
#pragma omp parallel for private(sum2,ibs,jbs,ene1,ene2,ene,tid,igrid,b,i,j)
    for (i = 0; i < ist.totalhomo; i++) {
      tid = omp_get_thread_num();	
      ene1 = eval[a] - eval[i];

      for (igrid = 0; igrid < ist.ngrid; igrid++) {
        // TODO: just put .re here for now - needs to be corrected
        // TODO: make sure complex conjugates are correct as well
        rho[tid*ist.ngrid+igrid].re = psi[i*ist.ngrid+igrid].re * psi[a*ist.ngrid+igrid].re;
        rho[tid*ist.ngrid+igrid].im = 0.0;
      }
      hartree(&rho[tid*ist.ngrid], potq, &poth[tid*ist.ngrid], ist, planfw[tid], planbw[tid], &fftwpsi[tid*ist.ngrid]);

      for (b = ist.nlumo; b < ist.nlumo+ist.totallumo; b++) {
        for (j = 0; j < ist.totalhomo; j++) {
          ene2 = eval[b] - eval[j];
          ene = ene1 - ene2;
	  
          for (sum2 = 0.0, igrid = 0; igrid < ist.ngrid; igrid++) {
            // TODO: just put .re here for now - needs to be corrected
            // TODO: make sure complex conjugates are correct as well
            sum2 += (poth[tid*ist.ngrid+igrid] * psi[b*ist.ngrid+igrid].re * psi[j*ist.ngrid+igrid].re);
          }
          sum2 *= par.dv;
          ibs = listibs[(a-ist.nlumo)*ist.totalhomo+i];
          jbs = listibs[(b-ist.nlumo)*ist.totalhomo+j];
          //bsmat[ibs*ist.ms2+jbs] -= 2.0 * sum2;
          fprintf (pf,"%ld %ld %ld %ld %ld %ld %g %g %g\n",a,i,b,j,ibs,jbs,ene1,ene2,sum2);
          fflush(0);
        }
      }
    }
  }

  fclose(pf); free(rho); free(listibs);
  
  return;
}

/***************************************************************************************/

