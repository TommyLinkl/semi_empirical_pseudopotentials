#include "fd.h"

int main(int argc, char *argv[])
{
  FILE *ppsi, *peval, *pf;  par_st par;
  long i, a, j, thomo, tlumo, flags=0;
  long indexfirsthomo;
  zomplex *potq, *potqx, *mux, *muy, *muz;
  lng_st ist;
  zomplex *psi, *psidummy;
  fftw_plan_loc *planfw, *planbw; fftw_complex *fftwpsi;
  double *eval, *de, *ksqr, *vx, *vy, *vz, *rx, *ry, *rz, *poth;
  double egap, tmp, *potl, *bsmat, *h0mat;

  /*************************************************************************/

  init_size(argc, argv, &par, &ist);
  
  /*************************************************************************/

  rx = (double*)calloc(ist.natom,sizeof(double));
  ry = (double*)calloc(ist.natom,sizeof(double));
  rz = (double*)calloc(ist.natom,sizeof(double));

  /*************************************************************************/  

  init_conf(rx, ry, rz, &par, &ist);

  /*************************************************************************/

  fftwpsi = fftw_malloc(sizeof (fftw_complex )*ist.ngrid*ist.nthreads);
  potl  = (double*)calloc(ist.ngrid,sizeof(double));
  potq  = (zomplex*)calloc(ist.ngrid,sizeof(zomplex));
  potqx  = (zomplex*)calloc(ist.ngrid,sizeof(zomplex));
  poth = (double*)calloc(ist.ngrid*ist.nthreads,sizeof(double));
  ksqr = (double*)calloc(ist.ngrid,sizeof(double));
  vx = (double*)calloc(ist.nx,sizeof(double));
  vy = (double*)calloc(ist.ny,sizeof(double));
  vz = (double*)calloc(ist.nz,sizeof(double));
  rx = (double*)calloc(ist.natom,sizeof(double));
  ry = (double*)calloc(ist.natom,sizeof(double));
  rz = (double*)calloc(ist.natom,sizeof(double));
  
  /**************************************************************************/

  init(vx,vy,vz,ksqr,&par,&ist);

  /*** initialization for the fast Fourier transform ***/
  fftw_plan_with_nthreads(ist.nthreads);
  
  planfw = (fftw_plan_loc*)calloc(ist.nthreads,sizeof(fftw_plan_loc));
  planbw = (fftw_plan_loc*)calloc(ist.nthreads,sizeof(fftw_plan_loc));
  for (i = 0; i < ist.nthreads; i++) { 
    planfw[i] = fftw_plan_dft_3d(ist.nz,ist.ny,ist.nx,&fftwpsi[i*ist.ngrid],&fftwpsi[i*ist.ngrid],FFTW_FORWARD,flags);
    planbw[i] = fftw_plan_dft_3d(ist.nz,ist.ny,ist.nx,&fftwpsi[i*ist.ngrid],&fftwpsi[i*ist.ngrid],FFTW_BACKWARD,flags);
  }
  init_pot(vx,vy,vz,potq,potqx,par,ist,planfw[0],planbw[0],&fftwpsi[0]);

  /**************************************************************************/

  eval = (double*)calloc(ist.nhomo+1,sizeof(double)); 
  de = (double*)calloc(ist.nhomo+1,sizeof(double)); 
  peval = fopen("eval.par", "r");
  for (i = 0; i < ist.nhomo+1; i++)
    fscanf(peval,"%ld %lg %lg",&a,&eval[i],&de[i]);
  fclose(peval);

  peval = fopen("eval.dat", "w");
  for (i = 0; i < ist.nhomo+1; i++) {
    fprintf(peval,"%ld %g %g\n",i,eval[i],de[i]);
  }
  fclose(peval);
  
  for (thomo = 0, i = ist.nhomo; i >= 0; i--) {
    if (de[i] < par.deps) thomo++;
    if (thomo == ist.totalhomo) {
      indexfirsthomo = i;
      break;
    }
  }
  free(eval);
  free(de);

  // TODO: Move this logic outside of main
  psi = (zomplex*)calloc(ist.ms*ist.nspinngrid,sizeof(zomplex));
  psidummy = (zomplex*)calloc(ist.nspinngrid,sizeof(zomplex));
  eval = (double*)calloc(ist.ms,sizeof(double)); 
  de = (double*)calloc(ist.ms,sizeof(double)); 
  peval = fopen("eval.par" , "r");
  ppsi = fopen("psi.par" , "r");

  // skip the beginning states in psi.par that are either not eigenstates
  // or too far away from the band-edge
  for (i = 0; i < indexfirsthomo; i++) {
    fscanf(peval,"%ld %lg %lg",&a,&tmp,&tmp);
    for (j = 0; j < ist.nspinngrid; j++) {
      fread(&(psidummy[j]), sizeof(zomplex), 1, ppsi);
    }
  }

  // read in the hole eigenstates and store them in psi
  for (a = thomo = 0; (thomo < ist.totalhomo) && (a <= ist.nhomo); a++) {
    fscanf(peval,"%ld %lg %lg",&a, &eval[thomo], &de[thomo]);
    for (j = 0; j < ist.nspinngrid; j++) fread(&(psidummy[j]), sizeof(zomplex), 1, ppsi);
    if (de[thomo] < par.deps) {
      for (j = 0; j < ist.nspinngrid; j++) {
        psi[thomo*ist.nspinngrid+j].re = psidummy[j].re;
        psi[thomo*ist.nspinngrid+j].im = psidummy[j].im;
      }
      thomo++;
    }
  }
  // skip the parts of psidummy that are not eigenstates
  for (i = ist.nhomo; i < ist.nlumo-1; i++) {
    fscanf(peval,"%ld %lg %lg",&a,&tmp,&tmp);
    for (j = 0; j < ist.nspinngrid; j++) fread(&(psidummy[j]), sizeof(zomplex), 1, ppsi);
  }
  // read in the quasielectron eigenstates and store them in psi
  for (tlumo = 0; tlumo < ist.totallumo; ) {
    fscanf(peval,"%ld %lg %lg",&a,&eval[thomo+tlumo],&de[thomo+tlumo]);
    for (j = 0; j < ist.nspinngrid; j++) fread(&(psidummy[j]), sizeof(zomplex), 1, ppsi);
    if (de[thomo+tlumo] < par.deps) {
      for (j = 0; j < ist.nspinngrid; j++) {
        psi[(thomo+tlumo)*ist.nspinngrid+j].re = psidummy[j].re;
        psi[(thomo+tlumo)*ist.nspinngrid+j].im = psidummy[j].im;
      }
      tlumo++;
    }
  }
  fclose(peval);
  fclose(ppsi);
  free(psidummy);
  
  printf("The number of hole states used in the BSE calculation = %ld\n", thomo);
  printf("The number of elec states used in the BSE calculation = %ld\n", tlumo);
  printf("The total number of band-edge states used in the BSE calculation = %ld\n", thomo + tlumo);
  printf("HOMO energy = %.10f %.10f\n", eval[thomo-1], eval[thomo-1]*AUTOEV);
  printf("LUMO energy = %.10f %.10f\n", eval[thomo], eval[thomo]*AUTOEV);
  printf("Fundamental gap  =  %.10f  %.10f\n", eval[thomo]-eval[thomo-1], (eval[thomo]-eval[thomo-1])*AUTOEV);
  fflush(stdout);

  /*************************************************************************/

  ist.ms = thomo+tlumo;
  ist.nlumo = thomo;
  ist.nhomo = thomo-1;
  ist.totallumo = tlumo;
  ist.totalhomo = thomo;
  
  /*************************************************************************/
  // normalize then print out (future densities and) z-projected 
  // electron densities for the quasiparticle states
  normalize_all(psi,par.dv,ist.ms,ist.nspinngrid);
  print_sp_pz(psi, de, vz, par, ist);

  /**************************************************************************/
  /*** this routine computes the coulomb coupling between
       single excitons.  On input - it requires the eigenstates stored in psi,
       the eigenvalues stored in eval and poth computed in init_pot.
       On output it stores the coulomb matrix elements on the disk
       in the following format: a, i, b, j, ene_ai, ene_bj, vjbai, vabji.
       a - the index of the electron in exciton Sai.
       i - the index of the hole in exciton Sai.
       b - the index of the electron in exciton Sbj.
       j - the index of the hole in exciton Sbj.
       ene_ai - the energy of exciton Sai.
       ene_bj - the energy of exciton Sbj.
       vjbai and vabji are the coulomb matrix elements needed to be used to
       generate the spin-depedent matrix elements as described by
       the last equation in our document.  ***/

  ist.ms2 = ist.totallumo * ist.totalhomo;
  bsmat = (double*)calloc(ist.ms2*ist.ms2,sizeof(double)); 
  h0mat = (double*)calloc(ist.ms2*ist.ms2,sizeof(double)); 

  single_coulomb_openmp(psi,potq,potqx,poth,eval,ist,par,planfw,planbw,fftwpsi,bsmat,h0mat);
  printf("made it through single_coulomb_openmp\n");

  ppsi = fopen("bs.dat", "w");
  for (i = 0; i < ist.ms2; i++, fprintf(ppsi,"\n"))
    for (j = 0; j < ist.ms2; j++)
      fprintf(ppsi,"%g ",bsmat[i*ist.ms2+j]);
  fclose(ppsi);

  ppsi = fopen("h0.dat", "w");
  for (i = 0; i < ist.ms2; i++, fprintf(ppsi,"\n"))
    for (j = 0; j < ist.ms2; j++)
      fprintf (ppsi,"%g ",h0mat[i*ist.ms2+j]);
  fclose(ppsi);

  mux = (zomplex*)calloc(ist.ms2,sizeof(zomplex));
  muy = (zomplex*)calloc(ist.ms2,sizeof(zomplex));
  muz = (zomplex*)calloc(ist.ms2,sizeof(zomplex));

  dipole(vx,vy,vz,psi,mux,muy,muz,eval,ist,par);
  printf("made it through dipole\n");
  bethe_salpeter(bsmat,h0mat,psi,vz,mux,muy,muz,ist,par);
  printf("made it through bethe_salpeter\n");
  
  /***********************************************************************/

  free(psi); free(potq); free(potqx); free(eval);  free(de); 
  free(ksqr); free(vx); free(vy);  free(vz);
  free(rx); free(ry); free(rz); free(poth);  free(potl);
  free(bsmat); free(h0mat); free(mux); free(muy); free(muz);
  free(planfw); free(planbw);

  printf("made it to the end of the program!\n");
  return 0;
}

/***********************************************************************/
