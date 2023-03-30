#include "fd.h"

/*****************************************************************************/

int main(int argc, char *argv[])
{
  FILE *ppsi;  zomplex *psi, *phi, *an; zomplex *psitot;
  fftw_plan_loc planfw, planbw; fftw_complex *fftwpsi;
  par_st par;   long_st  ist; xyz_st dipole; atm_st *atm;
  double *eval, *ksqr, *vx, *vy, *vz, *zn, *potl;
  double *el, *sige, *rx, *ry, *rz, tci, twi, *Elkb;
  double *dr, *vr, *potatom; 
  long *npot, idum;
  long jgrid, ispn, jms, jns, mssav, mstot, flags=0, tid;
  long *nl;
  time_t currentTime = time(NULL);
  nlc_st *nlc; 

  printf("This calculation began at: %s", ctime(&currentTime)); 
  fflush(stdout);

  /*** read initial setup from input.par ***/
  init_size(argc, argv, &par, &ist);

  /*** allocating memory ***/
  /*** the positions of the atoms in the x, y, and z directions ***/
  if ((rx = (double *) calloc(ist.natom, sizeof(double))) == NULL) nerror("rx");
  if ((ry = (double *) calloc(ist.natom, sizeof(double))) == NULL) nerror("ry");
  if ((rz = (double *) calloc(ist.natom, sizeof(double))) == NULL) nerror("rz");
  if ((atm = (atm_st *) calloc(ist.natom, sizeof(atm_st))) == NULL) nerror("atm");
  init_conf(rx, ry, rz, atm, &par, &ist);
  
  /*************************************************************************/
  /*** allocating memory ***/
  fftwpsi = fftw_malloc(sizeof (fftw_complex )*ist.ngrid);
  if ((psi = (zomplex *)calloc(ist.nspinngrid, sizeof(zomplex))) == NULL) nerror("psi");
  if ((phi = (zomplex *)calloc(ist.nspinngrid, sizeof(zomplex))) == NULL) nerror("phi");
  /*** the pseudopotential stored on the grid ***/
  if ((potl = (double *) calloc(ist.ngrid, sizeof(double))) == NULL) nerror("potl");
  /*** the kinetic energy stored on the grid ***/
  if ((ksqr = (double *) calloc(ist.ngrid, sizeof(double))) == NULL) nerror("ksqr");
  /*** the grid in the x, y, and z directions ***/
  if ((vx = (double *) calloc(ist.nx, sizeof(double))) == NULL) nerror("vx");
  if ((vy = (double *) calloc(ist.ny, sizeof(double))) == NULL) nerror("vy");
  if ((vz = (double *) calloc(ist.nz, sizeof(double))) == NULL) nerror("vz");
  /*** nonlocal pseudopotential parameters ***/

  /*** all wavefunctions in the energy range  ***/
  if ((psitot = (zomplex *) calloc(ist.nspin*ist.mstotngrid, sizeof(zomplex))) == NULL) nerror("psitot");
  /*** the filtered energies ***/
  if ((eval = (double *) calloc(ist.mstot, sizeof(double))) == NULL) nerror("eval");
  if ((sige = (double *) calloc(ist.mstot, sizeof(double))) == NULL) nerror("sige");
  /*** the energy of the filters ***/
  if ((el = (double *) calloc(ist.ms, sizeof(double))) == NULL) nerror("el");
   
  /**************************************************************************/
  /*** initialize the pseudopotential, the kinetic energy, ***/
  dr = (double *) calloc(ist.natomtype, sizeof(double));
  vr = (double *) calloc(ist.npot*ist.natomtype, sizeof(double));
  potatom = (double *) calloc(ist.npot*ist.natomtype, sizeof(double));
  npot = (long *) calloc(ist.natomtype, sizeof(long));
  init(potl, vx, vy, vz, ksqr, rx, ry, rz, atm, &par, el, &ist, dr, vr, potatom, npot);

  /*** initialization for the non-local potential ***/
  if ((nlc = (nlc_st *) calloc(ist.nsemicond*ist.nnlc, sizeof(nlc_st))) == NULL) nerror("nlc");
  if ((nl = (long *) calloc(ist.nsemicond, sizeof(long))) == NULL) nerror("nl");
  init_list(nlc, nl, vx, vy, vz, rx, ry, rz, atm, par, ist);

  /*** generate the KB projections ***/
  if (ist.flagkb == 1) {
    printf("Kleinman-Bylander method used for nonlocal spin-orbit pseudopotential\n");
    if ((Elkb = (double *) calloc(2*ist.nsemicond, sizeof(double))) == NULL) nerror("Elkb");
    init_kb(nlc, nl, Elkb, vx, vy, vz, rx, ry, rz, atm, par, ist, dr, vr, potatom, npot);
  }
  else { 
    printf("Matrix method for the nonlocal spin-orbit pseudopotential\n");
  }
  free(vr); free(dr); free(potatom); free(npot);
  
  /*** initialization for the fast Fourier transform ***/
  planfw = fftw_plan_dft_3d(ist.nz, ist.ny, ist.nx, fftwpsi, fftwpsi, FFTW_FORWARD, flags);
  planbw = fftw_plan_dft_3d(ist.nz, ist.ny, ist.nx, fftwpsi, fftwpsi, FFTW_BACKWARD, flags);
  
  /**************************************************************************/
  /*** calculate the energy range of the hamitonian ***/
  tci = (double)clock(); 
  twi = (double)time(NULL);
  get_energy_range(psi, phi, potl, vx, vy, vz, ksqr, ist, &par, nlc, nl, Elkb, planfw, planbw, fftwpsi);
  printf("done calculate energy range, CPU time (sec) %g, wall run time (sec) %g\n",
            ((double)clock()-tci)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-twi); 
  fflush(stdout);
 
  /**************************************************************************/
  /*** set parameters for the newton interpolation ***/
  par.dt = sqr((double)(ist.nc) / (2.5*par.dE));
  printf("nc = %ld dt = %g dE = %g\n",ist.nc,par.dt,par.dE); fflush(0);
  an = (zomplex *) calloc(ist.nc*ist.ms, sizeof(zomplex));
  zn = (double *) calloc(ist.nc, sizeof(double));
  coefficient(an, zn, el, par, ist);

  /**************************************************************************/
  /*** start filtering loop.  we run over ns cycles and calculate ***/
  /*** ms filtered states at each cycle ***/
  Randomize();  idum = -random();
  printf("seed = %ld\n", idum);  
  fflush(stdout);

  for (jns = 0; jns < ist.ns; jns++) {
    for (ispn = 0; ispn < ist.nspin; ispn++) {
      init_psi(&psi[ispn*ist.ngrid], ist, par, &idum);
    }
    for (jms = 0; jms < ist.ms; jms++) {
    	for (jgrid = 0; jgrid < ist.nspinngrid; jgrid++) {
    	  psitot[jns*ist.ms*ist.nspinngrid+ist.nspinngrid*jms+jgrid].re = psi[jgrid].re;
    	  psitot[jns*ist.ms*ist.nspinngrid+ist.nspinngrid*jms+jgrid].im = psi[jgrid].im;
    	}
    }
  }
  
  tci = (double)clock(); 
  twi = (double)time(NULL);
  omp_set_dynamic(0);
  omp_set_num_threads(ist.nthreads);
#pragma omp parallel for private(jns)
  for (jns = 0; jns < ist.ns; jns++) {
    tid = omp_get_thread_num();	
    filtering(&psitot[jns*ist.ms*ist.nspinngrid], potl, ksqr, an, zn, el, ist, par, nlc, nl, Elkb, tid, jns);
  } 
  printf("done calculating filter, CPU time (sec) %g, wall run time (sec) %g\n",
            ((double)clock()-tci)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-twi); 
  fflush(stdout);

  /*************************************************************************/
  /*** read all filtered states ***/
  /*ppsi = fopen("psi-filt.dat" , "r");
  fread (psitot,sizeof(zomplex),ist.mstot*ist.ngrid,ppsi);
  fclose(ppsi);*/

  /*** orthogonalize and normalize the filtered states using an svd routine ***/
  tci = (double)clock(); twi = (double)time(NULL);
  ist.mstot = portho((MKL_Complex16*)psitot,par.dv,ist);
  printf("mstot ortho = %ld\n",ist.mstot); fflush(0);
  normalize_all(&psitot[0],par.dv,ist.mstot,ist.nspinngrid,ist.nthreads);
  printf("done calculating ortho, CPU time (sec) %g, wall run time (sec) %g\n",
            ((double)clock()-tci)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-twi); 
  fflush(stdout);
    
  /***********************************************************************/
  /*** diagonalize the hamiltonian in the subspace spanned by the ***/
  /*** orthogonal filtered states. generate the eigenstates of the ***/
  /*** hamitonian within the desired energy reange ***/
  tci = (double)clock(); twi = (double)time(NULL);
  Hmat(psi,phi,(MKL_Complex16*)psitot,potl,ksqr,eval,ist,par,nlc,nl,Elkb,planfw,planbw,fftwpsi);
  normalize_all(&psitot[0],par.dv,ist.mstot,ist.nspinngrid,ist.nthreads);
  jms = ist.mstot;
  printf("done calculating Hmat, CPU time (sec) %g, wall run time (sec) %g\n",
            ((double)clock()-tci)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-twi);
  fflush(stdout);

  /*** write the eigenstates to a file ***/
  ppsi = fopen("psi.dat" , "w");
  fwrite (&psitot[0],sizeof(zomplex),jms*ist.nspinngrid,ppsi);
  fclose(ppsi);

  /*** calculate the standard deviation of these states ***/
  /*** this is used to check if there are ghost states ***/
  calc_sigma_E(psi, phi, psitot, potl, ksqr, sige, ist, par, nlc, nl, Elkb, planfw, planbw, fftwpsi);

  /*** write the eigenvalues and the standard deviation of the eigenstates ***/
  ppsi = fopen("eval.dat" , "w");
  for (jms = 0; jms < ist.mstot; jms++) fprintf (ppsi,"%ld %.16g %g\n", jms, eval[jms], sige[jms]); 
  fclose(ppsi);  

  /*** calculate the spin purity ***/
  calcAllSpinUpDownRatios(psitot, ist.nspinngrid, ist.mstot);

  /*************************************************************************/
  /*** free memeory ***/
  free(psitot); free(phi); free(psi); free(potl); free(eval); 
  free(el); free(ksqr); free(vx); free(vy); free(an); free(zn);
  free(vz); free(rx); free(ry); free(rz); free(sige); free(atm);
  free(nlc); free(nl);
  if (ist.flagkb == 1) free(Elkb);

  fftw_destroy_plan(planfw);
  fftw_destroy_plan(planbw);
  fftw_free(fftwpsi);
 
  currentTime = time(NULL);
  printf("This calculation ended at: %s", ctime(&currentTime)); fflush(0);
 
  exit(0);
}

/*****************************************************************************/
