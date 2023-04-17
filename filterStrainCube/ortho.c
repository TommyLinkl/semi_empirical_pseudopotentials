/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "fd.h"

/*****************************************************************************/

long portho(double *psi,double dv,long_st ist)
{
	long long lwork;  long long one=1, i, cutoff;
	long long ngrid = (long long)ist.ngrid, mstot = (long long)ist.mstot;
	double *work, *S, *vt;
	long *iwork;
	// MKL_INT info=0;
	long long info=0;

	// ORIGINAL memory allocations
	lwork = 5*(long long)(ist.mstot*ist.mstot+ist.ngrid);
	lwork = (long long)(3*mstot + ngrid + 5*mstot);
	S = (double*) calloc(mstot, sizeof(double));
	work = (double*) calloc(lwork, sizeof(double));

	// DGESDD_ memory allocations
	// lwork = (long long)(3*mstot + ngrid + 5*mstot*mstot + 4*mstot);
	// S = (double*) malloc((mstot+ngrid)* sizeof(double));
	// iwork = (long*) malloc(8*mstot* sizeof(long));
	// work = (double*) malloc(lwork* sizeof(double));
	// vt = (double*) malloc(ngrid*mstot* sizeof(double));

	// DUMMY MATRIX memory allocations
	// lwork = 5*(long long)(3*3 + 3 + 5*3);
	// S = (double*) calloc(3, sizeof(double));
	// work = (double*) calloc(lwork, sizeof(double));
	// double a[3*3] = {1, 0, 4, 2, 1, 0, 3, 1, 5};
	// int M_test=3, N_test=3;

	// LAPACKE memory allocations
	// MKL_INT lwork, info, one=1, i, cutoff;
	// MKL_INT ngrid = (lapack_int)ist.ngrid, mstot = (lapack_int)ist.mstot;
	// double superb[ngrid-1]; 
	// S = (double*) malloc(ist.mstot* sizeof(double));

	// printf("lwork = %ld \n", lwork);

	if (ngrid*mstot<2147483647){
		mkl_set_dynamic(0);
		mkl_set_num_threads(ist.nthreads);
		// omp_set_nested(1);
		omp_set_max_active_levels(1); 
	}
	else {
		mkl_set_dynamic(0);
		mkl_set_num_threads(1);
		// omp_set_nested(1);
		omp_set_max_active_levels(1); 
	}
	
	{
		// ORIGINAL 
		dgesvd_("O","N",&(ngrid),&(mstot),&(psi[0]),&(ngrid),&(S[0]),
		NULL,&one,NULL,&one,&(work[0]),&(lwork),&info);

		// DUMMY MATRIX 
		// dgesvd("O","N",&M_test,&N_test,&(a[0]),&M_test,&(S[0]),NULL,&one,NULL,&one,&(work[0]),&lwork,&info);
		// dgesvd("O","N",&M_test,&N_test, a,&M_test,S,NULL,&one,NULL,&one,work,&lwork,&info);
		
		// LAPACKE
		// info = LAPCKE_dgesvd(LAPACK_ROW_MAJOR, "O","N",ngrid,mstot,&(psi[0]),ngrid,&(S[0]),
		// NULL,one,NULL,one,&(superb)); 

		// DGESDD_ 
		// dgesdd_("O",&(ngrid),&(mstot),&(psi[0]),&(ngrid),&(S[0]),
		// NULL,&(ngrid),&(vt[0]),&(ngrid),&(work[0]),&(lwork),&(iwork[0]),&info);
	}
	
	printf("info = %ld\n",info);

	if (info != 0) {
		printf("error in dgesvd(2) %ld, exiting\n",info);
		exit(0);
	}
	
	for (cutoff = ist.mstot, i=0; i<ist.mstot; i++) {
		// printf("i, S[i] = %d, %g\n", i, S[i]);
		if ((S[i] / S[0]) < 1.0e-10) {
			cutoff = i;
			break;
		}
	}
	printf("cutoff is %ld\n",cutoff);
	
	free(work); 
	free(S);
	return (cutoff);
}

/*****************************************************************************/

