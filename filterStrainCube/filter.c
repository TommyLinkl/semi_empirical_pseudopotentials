/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "fd.h"

/****************************************************************************/

void filtering(double *psitot,double *potl,double *ksqr,zomplex *an,double *zn,double *el,long_st ist,par_st par,long tid,long jns,long *idum)
{
	FILE *pf; char str[100]; long flags = 0, jms, jgrid;
	fftw_plan_loc planfw, planbw; fftw_complex *fftwpsi;
	zomplex *psi, *phi;  double *eval;
	
	fftwpsi = fftw_malloc(sizeof (fftw_complex )*ist.ngrid);
	if ((psi = (zomplex*)calloc(ist.ngrid,sizeof(zomplex)))==NULL)nerror("psi");
	if ((phi = (zomplex*)calloc(ist.ngrid,sizeof(zomplex)))==NULL)nerror("phi");
	if ((eval = (double*)calloc(ist.ms,sizeof(double)))==NULL)nerror("eval");

	planfw = fftw_plan_dft_3d(ist.nz,ist.ny,ist.nx,fftwpsi,fftwpsi,FFTW_FORWARD,flags);
	planbw = fftw_plan_dft_3d(ist.nz,ist.ny,ist.nx,fftwpsi,fftwpsi,FFTW_BACKWARD,flags);
	
	/*** start from an initial random state ***/  
	//init_psi(psi,ist,par,&idum[tid]);
	
	for (jgrid = 0; jgrid < ist.ngrid; jgrid++) psi[jgrid].re = psitot[jgrid];
	
	/***********************************************************************/
	/*** filter the states and normalize them ***/
	filter(psi,phi,psitot,potl,ksqr,par,an,zn,ist,planfw,planbw,fftwpsi,tid,jns); 
	normalize_all(psitot,par.dv,ist.ms,ist.ngrid,ist.nthreads);
	
	/*** calculate and print the energy of the filtered states ***/
	energy_all(psi,phi,psitot,potl,ksqr,eval,ist,par,planfw,planbw,fftwpsi);

	sprintf (str,"eval-filt-%ld-%ld.dat",tid,jns);
	pf = fopen(str , "w");
	for (jms = 0; jms < ist.ms; jms++)
		fprintf (pf,"%ld %.16g %.16g\n",jms,eval[jms],el[jms]);
	fclose(pf);

	/*** write the filtered states to a file ***/
	sprintf (str,"psi-filt-%ld-%ld.dat",tid,jns);
	pf = fopen(str , "w");
	fwrite (psitot,sizeof(double),ist.ms*ist.ngrid,pf);
	fclose(pf);

	/***
	sprintf (str,"NEW_psi-filt-%ld-%ld.dat",tid,jns);
	pf = fopen(str , "w");
	printf("number of grid points: %ld \n", ist.ngrid);
	for (jgrid = 0; jgrid < ist.ngrid; jgrid++) {
		fprintf(pf,"%g \n",psitot[jgrid]);
	}
	fclose(pf);
	***/
	
	free(psi); free(eval); free(phi);
	fftw_destroy_plan(planfw);
	fftw_destroy_plan(planbw);
	fftw_free(fftwpsi);
	return;
}

/****************************************************************************/

void filter(zomplex *psin,zomplex *psim1,double *psi0,double *potl,double *ksqr,par_st par,zomplex *an,double *zn,long_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,long tid,long jns)
{
	FILE *pf; char str[100];
	long i, j, ie, ncie, ieg;

	for (ie = 0; ie < ist.ms; ie++){
		ncie = ist.nc*ie; ieg = ie * ist.ngrid;
		for (i = 0; i < ist.ngrid; i++) psi0[ieg+i] = (an[ncie+0].re*psin[i].re);
	}

	sprintf (str,"prop%ld.dat",tid);
	for (pf = fopen(str , "w"), j = 1; j < ist.nc; j++){
		memcpy(&psim1[0],&psin[0],ist.ngrid*sizeof(psim1[0]));
		hnorm(psim1,psin,potl,ksqr,par,zn[j-1],ist,planfw,planbw,fftwpsi);
		
		for (ie = 0; ie < ist.ms; ie++){
			ncie = ist.nc*ie; ieg = ie * ist.ngrid;
			for (i = 0; i < ist.ngrid; i++)
				psi0[ieg+i] += (an[ncie+j].re*psin[i].re);
		}
		if (!(j % 100)) {fprintf (pf,"%ld %ld\n",j,jns); fflush(pf);}
	}
	fclose(pf);
	return;
}

/****************************************************************************/