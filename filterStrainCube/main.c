/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "fd.h"

/*****************************************************************************/

int main(int argc, char *argv[]) {
	FILE *ppsi;  zomplex *psi, *phi, *an; long idum;
	fftw_plan_loc planfw, planbw; fftw_complex *fftwpsi;
	par_st par;   long_st  ist; xyz_st dipole; atm_st *atm;
	double *eval, *ksqr, *vx, *vy, *vz, *zn, *potl;
	double *psitot, *el, *sige, *rx, *ry, *rz, tci, twi;
	long jgrid, jms, jns, mssav, mstot, flags=0, tid;
	FILE *test_psi, *test_psitot;
	char str[100]; 

	/*** read initial setup from input.par ***/
	writeCurrentTime(stdout);
	writeSeparation(stdout);
	init_size(argc,argv,&par,&ist);
	
	/*************************************************************************/
	/*** allocating memory ***/

	fftwpsi = fftw_malloc(sizeof (fftw_complex )*ist.ngrid);
	//if ((fftwpsi  = (fftw_complex*)calloc(ist.ngrid,sizeof(fftw_complex)))==NULL)nerror("fftwpsi");
	if ((psi   = (zomplex*)calloc(ist.ngrid,sizeof(zomplex)))==NULL)nerror("psi");
	if ((phi   = (zomplex*)calloc(ist.ngrid,sizeof(zomplex)))==NULL)nerror("phi");
	/*** the pseudopotential stored on the grid ***/
	if ((potl  = (double*)calloc(ist.ngrid,sizeof(double)))==NULL)nerror("potl");
	/*** the kinetic energy stored on the grid ***/
	if ((ksqr = (double*)calloc(ist.ngrid,sizeof(double)))==NULL)nerror("ksqr");
	/*** the grid in the x, y, and z directions ***/
	if ((vx = (double*)calloc(ist.nx,sizeof(double)))==NULL)nerror("vx");
	if ((vy = (double*)calloc(ist.ny,sizeof(double)))==NULL)nerror("vy");
	if ((vz = (double*)calloc(ist.nz,sizeof(double)))==NULL)nerror("vz");
	/*** the positions of the atoms in the x, y, and z directions ***/
	if ((rx = (double*)calloc(ist.natom,sizeof(double)))==NULL)nerror("rx");
	if ((ry = (double*)calloc(ist.natom,sizeof(double)))==NULL)nerror("ry");
	if ((rz = (double*)calloc(ist.natom,sizeof(double)))==NULL)nerror("rz");
	if ((atm = (atm_st*)calloc(ist.natom,sizeof(atm_st)))==NULL)nerror("atm");

	/*** all wavefunctions in the energy range  ***/
	if ((psitot = (double*)calloc(ist.mstot*ist.ngrid,sizeof(double)))==NULL)nerror("psitot");
	/*** the filtered energies ***/
	if ((eval  = (double*)calloc(ist.mstot,sizeof(double)))==NULL)nerror("eval");
	if ((sige = (double*)calloc(ist.mstot,sizeof(double)))==NULL)nerror("sige");
	/*** the energy of the filters ***/
	if ((el   = (double*)calloc(ist.ms,sizeof(double)))==NULL)nerror("el");
	 
	/**************************************************************************/
	/*** initialize the pseudopotential, the kinetic energy, ***/
	/*** the fourier transorm, the grid, and the energy window ***/
	init(potl,vx,vy,vz,ksqr,rx,ry,rz,atm,&par,el,&ist,&planfw,&planbw,fftwpsi);
	/*** initialization for the fast Fourier transform ***/
	planfw = fftw_plan_dft_3d(ist.nz,ist.ny,ist.nx,fftwpsi,fftwpsi,FFTW_FORWARD,flags);
	planbw = fftw_plan_dft_3d(ist.nz,ist.ny,ist.nx,fftwpsi,fftwpsi,FFTW_BACKWARD,flags);
	
	/**************************************************************************/
	/*** calcualte the energy range of the hamitonian ***/
	tci = (double)clock(); twi = (double)time(NULL);
	get_energy_range(psi,phi,potl,vx,vy,vz,ksqr,&par,ist,planfw,planbw,fftwpsi);
	printf("done calculate energy range, CPU time (sec) %g, wall run time (sec) %g\n",
		((double)clock()-tci)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-twi); fflush(0);
	
	/**************************************************************************/
	/*** set parameters for the newton interpolation ***/
	par.dt = sqr((double)(ist.nc) / (2.5 * par.dE));
	printf ("nc = %ld dt = %g dE = %g\n",ist.nc,par.dt,par.dE); fflush(0);
	an = (zomplex*)calloc(ist.nc*ist.ms,sizeof(zomplex));
	zn = (double*)calloc(ist.nc,sizeof(double));
	coefficient(an,zn,el,par,ist);

	/**************************************************************************/
	/*** start filtering loop.  we run over ns cycles and calculate ***/
	/*** ms filtered states at each cycle ***/
	Randomize();  idum = -random();
	printf("seed = %ld\n",idum);  fflush(0);

	for (jns = 0; jns < ist.ns; jns++) {
		init_psi(psi,ist,par,&idum);
		for (jms = 0; jms < ist.ms; jms++) {
			for (jgrid = 0; jgrid < ist.ngrid; jgrid++) {
				psitot[jns*ist.ms*ist.ngrid+ist.ngrid*jms+jgrid] = psi[jgrid].re;
			}
		}
	}
	
	tci = (double)clock(); twi = (double)time(NULL);
	omp_set_dynamic(0);
	omp_set_num_threads(ist.nthreads);
#pragma omp parallel for private(jns)
	for (jns = 0; jns < ist.ns; jns++){
		tid = omp_get_thread_num();	
		filtering(&psitot[jns*ist.ms*ist.ngrid],potl,ksqr,an,zn,el,ist,par,tid,jns,&idum);
	} 
	printf("done calculating filter, CPU time (sec) %g, wall run time (sec) %g\n",
		((double)clock()-tci)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-twi); fflush(0);

	/*************************************************************************/
	/*** read all filtered states ***/
	/*ppsi = fopen("psi-filt.dat" , "r");
	fread (psitot,sizeof(double),ist.mstot*ist.ngrid,ppsi);
	fclose(ppsi);*/

	normalize_all(psitot,par.dv,ist.mstot,ist.ngrid,ist.nthreads);
	// printf("par.dv = %g \n", par.dv);

	/*** Tommy debug ***/
	/***
	test_psitot = fopen("DEBUG_psitot.dat", "w");
	for (jns = 0; jns < ist.ns; jns++) {
		for (jms = 0; jms < ist.ms; jms++) {
			printf("Norm of filtered state ns=%d, ms=%d is: %g\n", jns, jms, normd(&psitot[jns*ist.ms*ist.ngrid+ist.ngrid*jms], par.dv, ist.ngrid, ist.nthreads));
			
			sprintf (str,"DEBUG_psitot-%ld-%ld.dat",jns,jms);
			test_psi = fopen(str, "w");
			for (jgrid = 0; jgrid < ist.ngrid; jgrid++) {
				fprintf(test_psitot,"%g  ",psitot[jns*ist.ms*ist.ngrid + jms*ist.ngrid + jgrid]);
				fprintf(test_psi,"%g \n", psitot[jns*ist.ms*ist.ngrid + jms*ist.ngrid + jgrid]);
			}
			fclose(test_psi);
			fprintf(test_psitot,"\n");
		}
	}
	fclose(test_psitot);
	***/

	/*** orthogonalize and normalize the filtered states using an svd routine ***/
	tci = (double)clock(); twi = (double)time(NULL);
	ist.mstot = portho(psitot,par.dv,ist);
	printf("mstot ortho = %ld\n", ist.mstot); fflush(0);
	normalize_all(psitot,par.dv,ist.mstot,ist.ngrid,ist.nthreads);
	printf("done calculating ortho, CPU time (sec) %g, wall run time (sec) %g\n",
		((double)clock()-tci)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-twi); fflush(0);
		
	/***********************************************************************/
	/*** diagonalize the hamiltonian in the subspace spanned by the ***/
	/*** orthogonal filtered states. generate the eigenstates of the ***/
	/*** hamitonian within the desired energy reange ***/
	tci = (double)clock(); twi = (double)time(NULL);
	Hmatreal(psi, phi, psitot, potl, ksqr, eval,ist, par, planfw, planbw, fftwpsi);
	normalize_all(psitot, par.dv, ist.mstot, ist.ngrid, ist.nthreads);
	jms = ist.mstot;
	printf("done calculating Hmat, CPU time (sec) %g, wall run time (sec) %g\n",
		((double)clock()-tci)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-twi);
	fflush(0);

	/*** write the eigenstates to a file ***/
	ppsi = fopen("psi.dat", "w");
	fwrite(psitot, sizeof(double), jms*ist.ngrid, ppsi);
	fclose(ppsi);

	/*** calculate the standard deviation of these states ***/
	/*** this is used to check if there are ghost states ***/
	calc_sigma_E(psi, phi, psitot, potl, ksqr, sige, ist, par, planfw, planbw, fftwpsi);

	/*** write the eigenvalues and the standard deviation of the eigenstates ***/
	ppsi = fopen("eval.dat", "w");
	for (jms = 0; jms < ist.mstot; jms++) { 
		fprintf(ppsi, "%ld %.16g %g\n", jms, eval[jms], sige[jms]); 
	}
	fclose(ppsi);  
	
	/*************************************************************************/
	/*** free memeory ***/
	free(psitot); free(phi); free(psi); free(potl); free(eval); 
	free(el); free(ksqr); free(vx); free(vy); free(an); free(zn);
	free(vz); free(rx); free(ry); free(rz); free(sige); free(atm);

	fftw_destroy_plan(planfw);
	fftw_destroy_plan(planbw);
	fftw_free(fftwpsi);

	/*************************************************************************/
	/*** Print out the time that the computation finished ***/
	writeSeparation(stdout);
	writeCurrentTime(stdout);

	/*************************************************************************/
	/*** Exit the program ***/
	exit(0);
}
